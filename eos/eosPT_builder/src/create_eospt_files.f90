! ***********************************************************************
!
!   Copyright (C) 2009-2019  Bill Paxton & The MESA Team
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

      module create_eosPT_files
      
      use eos_def
      use eos_lib
      use const_def
      use chem_def
      
      implicit none



      logical, parameter :: use_shared_data_dir = .true.   ! MUST BE .true. FOR RELEASE
      !logical, parameter :: use_shared_data_dir = .false.

      integer, parameter :: version_number = 51 ! update this to force rebuilding of caches
      ! update min_version in eosPT_load_tables to force rebuild of data files

      integer :: ix, ios, io_unit, handle
      double precision :: whichz

      ! control what is done by saving the Z in whichz.txt and, if necessary, setting the Xs below
      
      integer, parameter :: num_Xs = 5, num_Zs = 1 ! currently limited to a single Z data file per run
      double precision :: Xs(num_Xs)
      
      
      contains
      
      subroutine do_create_eosPT_files
         integer :: iscvh
         logical :: expanded_scvh
         !open(UNIT=io_unit, FILE=trim("expanded_scvh_coverage.txt"), &
         !     ACTION='READ', STATUS='OLD', IOSTAT=ios)
         !if (ios /= 0) call do_stop('failed to open expanded_scvh_coverage.txt')
         !read(io_unit, fmt=*, iostat=ios) iscvh
         !if (ios /= 0) call do_stop('failed to read 0 or 1 from expanded_scvh_coverage.txt')
         !close(io_unit)
         
         iscvh = 1 ! expanded_scvh
         expanded_scvh = (iscvh /= 0)
      
         Xs(1:num_Xs) = (/ 0.00d0, 0.20d0, 0.40d0, 0.60d0, 0.80d0 /)
         io_unit = 40
         if (expanded_scvh) then
            whichz = 0
         else
            open(UNIT=io_unit, FILE=trim("whichz.txt"), ACTION='READ', STATUS='OLD', IOSTAT=ios)
            if (ios /= 0) call do_stop('failed to open whichz.txt')
            read(io_unit, fmt=*, iostat=ios) whichz
            if (ios /= 0) call do_stop('failed to read Z from whichz.txt')
         end if
         
         do ix = 1, num_Xs
            call make_eosPT_files(whichz, Xs(ix), expanded_scvh)
         end do
         if (whichz == 0) call make_eosPT_files(whichz, 1d0, expanded_scvh)
         
      end subroutine do_create_eosPT_files
      
      
      subroutine make_eosPT_files(Z, X, expanded_scvh)
         use eval_eosPT
         double precision, intent(in) :: Z, X
         logical, intent(in) :: expanded_scvh

         double precision :: &
            abar,zbar, logPgas, Pgas, logW_min, logW_max, del_logW, &
            logRho, logE, logS, logT_min, logT_max, del_logT, &
            logW, logT, T, Rho, E, S, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta, logRho_guess
         integer :: ierr, i, j, num_logWs, num_logTs, call_number
         character (len=256) :: dir, fname, prefix
            
         include 'formats.dek'

         ierr = 0

         dir = 'data' ! where to put the new data files

         call Init_Composition(X, Z, abar, zbar)
         
         if (expanded_scvh) then
            logW_min = -12d0
            logW_max = -3.0d0
            del_logW = 0.01d0 ! 0.02d0
            logT_min = 2.1d0
            logT_max = 7.0d0
            del_logT = 0.01d0 ! 0.02d0
         else
            logW_min = -20d0
            logW_max = -2.9d0
            del_logW = 0.1d0
            logT_min = 2.1d0
            logT_max = 8.20d0
            del_logT = 0.02d0
         end if
            
         num_logWs = 1 + int((logW_max - logW_min) / del_logW)
         num_logTs = 1 + int((logT_max - logT_min) / del_logT)
         
         if (expanded_scvh) then
            prefix =  '/eosPT_SCVH_data/mesa-eosPT_SCVH'
         else
            prefix = '/eosPT_data/mesa-eosPT'
         end if
         if (X < 0.005d0) then
            write(fname,'(a,a,i1,a)') trim(dir), &
               trim(prefix) // '_0', floor(100d0*Z + 0.5d0), 'z00x.data'
         else if (X < 1) then
            write(fname,'(a,a,i1,a,i2,a)') trim(dir), trim(prefix) // '_0',  &
               floor(100d0*Z + 0.5d0), 'z', floor(100d0*X + 0.5d0), 'x.data'
         else
            fname = trim(dir) // trim(prefix) // '_00z100x.data'
         end if
      
         write(*,*) trim(fname)

         open(unit=io_unit, file=trim(fname), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname)
            call mesa_error(__FILE__,__LINE__)
         end if
      
         write(io_unit,'(99(a14))') &
            'version', 'X', 'Z', 'num logTs', 'logT min', 'logT max', 'del logT',  &
            'num logWs', 'logW min', 'logW max', 'del logW'
      
         write(io_unit,'(i14,2f14.4,2(i10,4x,3(f14.4)))') &
               version_number, X, Z, num_logTs, logT_min, logT_max, del_logT, &
               num_logWs, logW_min, logW_max, del_logW
         
         call_number = 0
         do i = 1, num_logWs
            logW = logW_min + (i-1) * del_logW
            
            if (mod(i,10) == 0) write(*,*) 'logW', i

            write(io_unit,'(/,7x,a)') 'logW = logPgas - 4*logT' 
            write(io_unit,'(2x,f14.6/)') logW
            write(io_unit,'(a7,10x,99(a11,7x))') 'logT',  &
               'logRho', 'logE', 'logS', 'chiRho', 'chiT', &
               'Cp', 'Cv', 'dE_dRho', 'dS_dT', 'dS_dRho', &
               'mu', 'log_free_e', 'gamma1', 'gamma3', 'grad_ad', 'eta'
         
            logRho_guess = 0
            do j = 1, num_logTs
               logT = logT_min + (j-1) * del_logT
               T = 10**logT
               logPgas = logW + 4*logT
               Pgas = 10**logPgas
               ierr = 0
               call_number = call_number + 1
               
               call eos_4_Pgas_T( &
                     logW, logT, T, logPgas, Pgas, &
                     abar, zbar, X, Z, logRho_guess, &
                     logRho, logE, logS, chiRho, chiT, &
                     Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
                     mu, log_free_e, gamma1, gamma3, grad_ad, eta, &
                     expanded_scvh, call_number, ierr)
               if (ierr /= 0) then
                  write(*,*) 'logT', logT
                  write(*,*) 'logPgas', logPgas
                  call do_stop('failed in eos_4_PT')
               end if
               
               logRho_guess = logRho

               write(io_unit,'(f11.6,99(e18.9))') logT, &
                     logRho, logE, logS, chiRho, chiT, &
                     Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
                     mu, log_free_e, gamma1, gamma3, grad_ad, eta
                         
            end do
         
            write(io_unit,*)
         
         end do
      
         write(io_unit,*)
         write(io_unit,*)

         close(io_unit)
      
         write(*,*)

      end subroutine make_eosPT_files


      subroutine setup_eos
         use const_lib
         use chem_lib
         use scvh_eval
         use math_lib
         use eval_eosPT, only: eos_handle

         character (len=256) :: my_mesa_dir, data_dir, scvh_data_dir
         integer :: info
         logical, parameter :: use_cache = .true.
         
         info = 0
         
         my_mesa_dir = '../..'         
         call const_init(my_mesa_dir,info)     
         if (info /= 0) then
            write(*,*) 'const_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if        
         
         call math_init()

         if (use_shared_data_dir) then
            data_dir = '../../data' ! test using shared data
         else
            data_dir = '../data' ! test using local data
            
            write(*,*)
            write(*,*) 'TESTING WITH LOCAL DATA'
            write(*,*)
            
         end if
         
         call chem_init('isotopes.data', info)
         if (info /= 0) then
            write(*,*) 'chem_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if

         call eos_init('', '', '', use_cache, info)
         if (info /= 0) then
            write(*,*) 'failed in eos_init'
            call mesa_error(__FILE__,__LINE__)
         end if

         eos_handle = alloc_eos_handle(info)
         if (info /= 0) then
            write(*,*) 'failed in alloc_eos_handle'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         scvh_data_dir = '../eosDT_builder/eos_input_data/'
         call setup_scvh(scvh_data_dir)
      
      end subroutine setup_eos


      subroutine do_stop(str)
         character (len=*) :: str
         write(*,*) trim(str)
         call mesa_error(__FILE__,__LINE__)
      end subroutine do_stop
      
      
      end module create_eosPT_files
