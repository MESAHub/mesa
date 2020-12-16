! ***********************************************************************
!
!   Copyright (C) 2009  Bill Paxton, Frank Timmes
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      program create_eosDT_files
		use eos_def
		use helm
		use helm_alloc
		use helm_opal_scvh_driver
		use const_def
      use const_lib
		
      implicit none



      real(dp) ::
     >      del_logT = 0.02d0, 
     >      del_logQ = 0.03d0, 
     >      logT_min = 2.1d0, 
     >      logT_rad_max = 8.2d0,  
     >      logT_no_rad_max = 7.7d0,  
     >      logQ_min = -10.09d0, 
     >      logQ_max = 5.69d0, 
     >      logRho_min = -14d0

      real(dp) :: logT_max


		integer, parameter :: version_number = 51 ! update this to force rebuiding of caches
		! update min_version in eosDT_load_tables to force rebuild of data files

      integer :: ix, io_unit, ios, info, irad
		real(dp) :: whichz

      ! control what is done by saving the Z in whichz.txt and, if necessary, setting the Xs below
      
      integer, parameter :: num_Xs = 5, num_Zs = 1 ! currently limited to a single Z data file per run
      real(dp) :: Xs(num_Xs)
      
      Xs(1:num_Xs) = (/ 0.80d0, 0.00d0, 0.20d0, 0.40d0, 0.60d0 /)
		
      io_unit = 40
      open(UNIT=io_unit, FILE=trim("whichz.txt"), ACTION='READ', STATUS='OLD', IOSTAT=ios)
      if (ios /= 0) call do_stop('failed to open whichz.txt')
      read(io_unit, fmt=*, iostat=ios) whichz
      if (ios /= 0) call do_stop('failed to read Z from whichz.txt')
      close(io_unit)

      !open(UNIT=io_unit, FILE=trim("include_radiation.txt"), ACTION='READ', STATUS='OLD', IOSTAT=ios)
      !if (ios /= 0) call do_stop('failed to open include_radiation.txt')
      !read(io_unit, fmt=*, iostat=ios) irad
      !if (ios /= 0) call do_stop('failed to read 0 or 1 from include_radiation.txt')
      !close(io_unit)
      
      irad = 1 ! include radiation

      call setup_eos
      
      do ix = 1, num_Xs
         call Make_EoS_Files(whichz, Xs(ix), irad /= 0)
      end do
      
      if (whichz == 0) call Make_EoS_Files(whichz, 1d0, irad /= 0)
      
      call free_helm_table(eos_ht)
      

      contains
      
      
      subroutine do_stop(str)
         character (len=*) :: str
         write(*,*) trim(str)
         stop 1
      end subroutine do_stop

      
      subroutine Make_EoS_Files(Z_in, X_in, include_radiation)
		use helm_opal_scvh_driver
		
      real(dp), intent(in) :: Z_in, X_in
      logical, intent(in) :: include_radiation
      
      character (len=256) :: dir, fname
      integer :: io_unit, i, j, num_logTs, num_logQs, info
         
      real(dp) :: abar,zbar,z53bar,X,Z

      real(dp) :: logQ, logRho, logT, Rho, T
      real(dp) :: logPgas, logE, logS, chiRho, chiT, Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >      mu, free_e, gamma1, gamma3, grad_ad, eta, HELM_fraction
      character (len=64) :: fname_prefix
      
		logical :: helm_only = .false., opal_scvh_only = .false., 
     >			opal_only = .false., scvh_only = .false., search_for_SCVH = .true.
      
      !opal_only = .true.
      scvh_only = .true.
      
      dir = 'data' ! where to put the new data files
      io_unit = 40
      
      Z = Z_in; X = X_in
      call get_composition_info(X, Z, abar, zbar, z53bar)

!..other initialization
		info = 0
      
      if (opal_only) then
         fname_prefix = '/eosOPAL_data/mesa-OPAL_0'
      else if (scvh_only) then
         fname_prefix = '/eosSCVH_data/mesa-SCVH_0'
         del_logT = 0.02d0
         del_logQ = 0.03d0
         logT_min = 2.1d0 
         logT_rad_max = 6.5d0    
         logT_no_rad_max = 6.5d0 
         logQ_min = -2.09d0 
         logQ_max = 5.3d0 
         logRho_min = -9d0
      else
         fname_prefix = '/eosDT_data/mesa-eosDT_0'
      end if

      logT_max = logT_rad_max
      ! NOTE: if you change logT_max,
      ! you should also change logT1 and logT2 in eos_regions_defs.dek
      ! logT1 should = the new logT_max, and logT2 should be about 0.1 smaller.
      ! between logT1 and logT2, the mesa tables are a blend of HELM and OPAL.
      ! so you don't get pure OPAL until logT < logT2.
            
      num_logQs = 1 + int((logQ_max - logQ_min) / del_logQ)
      num_logTs = 1 + int((logT_max - logT_min) / del_logT)
      
      if (X < 0.005d0) then
         write(fname,'(a,a,i1,a)') trim(dir), trim(fname_prefix), floor(100d0*Z + 0.5d0), 'z00x.data'
      else if (X < 1) then
         write(fname,'(a,a,i1,a,i2,a)') trim(dir), trim(fname_prefix), 
     >		floor(100d0*Z + 0.5d0), 'z', floor(100d0*X + 0.5d0), 'x.data'
      else
         fname = trim(dir) // trim(fname_prefix) // '0z100x.data'
      end if
      
      write(*,*) trim(fname)

      open(unit=io_unit,file=trim(fname))
      
      write(io_unit,'(99(a14))') 'version', 'X', 'Z', 'num logTs', 'logT min', 'logT max', 'del logT', 
     1		'num logQs', 'logQ min', 'logQ max', 'del logQ'
      
      write(io_unit,'(i14,2f14.4,2(i10,4x,3(f14.4)))') 
     >      version_number, X, Z, num_logTs, logT_min, logT_max, del_logT,
     >		num_logQs, logQ_min, logQ_max, del_logQ

      do i = 1, num_logQs
         logQ = logQ_min + (i-1) * del_logQ

         write(io_unit,'(/,7x,a)') 'logQ = logRho - 2*logT + 12' 
         write(io_unit,'(2x,f14.6/)') logQ
         write(io_unit,'(a4,1x,3(a9,1x),7(a12,1x),1(a7,1x),1(a11),3(a9,1x),2(a9,1x))') 
     >            'logT', 'logPgas', 'logE', 'logS', 'chiRho', 'chiT', 
     >            'Cp', 'Cv', 'dE_dRho', 'dS_dT', 'dS_dRho', 
     >            'mu', 'log_free_e', 'gamma1', 'gamma3', 'grad_ad', 'eta', 'HELM'
         
         do j = 1, num_logTs
            logT = logT_min + (j-1) * del_logT
            T = 10**logT
            logRho = max(logRho_min,logQ + 2*logT - 12d0)
            Rho = 10**logRho
            info = 0
            call helm_opal_scvh(
     >            helm_only, opal_scvh_only, opal_only, scvh_only, search_for_SCVH,
     >				include_radiation, logT, logRho, T, Rho, abar, zbar, z53bar, X, Z,
     >            logPgas, logE, logS, chiRho, chiT, Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >            mu, free_e, gamma1, gamma3, grad_ad, eta, HELM_fraction, data_dir, info)
				if (info /= 0) then
					write(*,*) 'logT', logT
					write(*,*) 'logRho', logRho
					call do_stop('failed in helm_opal_scvh')
				end if

            write(io_unit,'(f4.2,3(f10.5),7(1pe13.5),1(0pf9.5),4(0pf10.5),2(0pf11.5))') 
     >         logT, logPgas, logE, logS, chiRho, chiT, Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >         mu, log10(max(1d-99,free_e)), gamma1, gamma3, grad_ad, eta, HELM_fraction
         end do
         
         write(io_unit,*)
         
      end do
      
      write(io_unit,*)
      write(io_unit,*)

      close(io_unit)
      
      write(*,*)

      end subroutine Make_EoS_Files
      
      end

