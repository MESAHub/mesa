c ***********************************************************************
!
!   Copyright (C) 2006, 2007, 2008  Bill Paxton, Frank Timmes
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
c ***********************************************************************

      module helm_opal_scvh_driver
      use chem_def
      implicit none
      
      
      character (len=256) :: data_dir
		integer, parameter :: imax = 261, jmax = 101  ! dimensions of our version of helm table
      

		contains


		subroutine setup_eos
		   use helm_alloc
		   use eos_def
		   use const_lib
		   use chem_lib
         use utils_lib
         use math_lib
		
		   integer :: ierr
         character (len=256) :: eos_file_prefix, my_mesa_dir, cache_dir, temp_cache_dir
         logical, parameter :: use_cache = .true.
		      
	      data_dir = 'eos_input_data' ! where to find the input data files
         cache_dir = trim(data_dir) // '/cache'
         temp_cache_dir = '.mesa_temp_cache/cache'
         call mkdir(temp_cache_dir)
                  
         my_mesa_dir = '../..'         
         call const_init(my_mesa_dir,ierr)     
      	if (ierr /= 0) then
      	   write(*,*) 'const_init failed'
      	   stop 1
      	end if        

         call math_init()

         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'chem_init failed'
            stop 1
         end if

         call eos_def_init
         
         ierr = 0
			call alloc_helm_table(eos_ht, imax, jmax, ierr)
			if (ierr /= 0) then
			   write(*,*) 'alloc helm table failed'
			   stop 1
			end if
		
			call read_helm_table(eos_ht, data_dir, cache_dir, temp_cache_dir, use_cache, ierr) ! initialize helm
			if (ierr /= 0) then
			   write(*,*) 'read helm table failed'
			   stop 1
			end if
			
			eos_ht% with_coulomb_corrections = .true.
			
			if (.not. eos_ht% with_coulomb_corrections) write(*,*) 'no coulomb corrections for helm'

		end subroutine setup_eos

      
      subroutine get_composition_info(X_in, Z_in, abar, zbar, z53bar)
         use chem_def, only: ih1, ihe4, ic12, in14, io16, ine20, img24
         use chem_lib, only: basic_composition_info
         double precision, intent(in) :: X_in, Z_in
         double precision, intent(out) :: abar, zbar, z53bar

         integer, parameter :: ionmax = 7
         double precision :: xmass(ionmax), X, Y, Z, z2bar, ye, mass_correction, sumx
         integer :: i, chem_id(ionmax)

			double precision, parameter :: Zfrac_C = 0.173312d0
			double precision, parameter :: Zfrac_N = 0.053177d0
			double precision, parameter :: Zfrac_O = 0.482398d0
			double precision, parameter :: Zfrac_Ne = 0.098675d0
			double precision, parameter :: Zfrac_Mg = 1d0 - (Zfrac_C + Zfrac_N + Zfrac_O + Zfrac_Ne)
         
         X = X_in
         Z = Z_in
         
         Y = 1 - (X+Z)
         
         i = 1
         chem_id(i) = ih1
         xmass(i) = X

         i = i+1
         chem_id(i) = ihe4
         xmass(i) = Y
         
         i = i+1
         chem_id(i) = ic12
         xmass(i) = Z * Zfrac_C
         
         i = i+1
         chem_id(i) = in14
         xmass(i) = Z * Zfrac_N
         
         i = i+1
         chem_id(i) = io16
         xmass(i) = Z * Zfrac_O
         
         i = i+1
         chem_id(i) = ine20
         xmass(i) = Z * Zfrac_Ne
         
         i = i+1
         chem_id(i) = img24
         xmass(i) = Z * Zfrac_Mg
         
         call basic_composition_info( 
     >		ionmax, chem_id, xmass, X, Y, Z, 
     >		abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

      end subroutine get_composition_info


      subroutine helm_opal_scvh(
     >					helm_only, opal_scvh_only, opal_only, scvh_only, search_for_SCVH,
     >					include_radiation, logT, logRho, temp, den, abar, zbar, z53bar, X, Z,
     >               logPgas, logE, logS, chiRho, chiT, 
     >               Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >               mu, free_e, gamma1, gamma3, grad_ad, eta, HELM_fraction,
     >               data_dir, ierr)

!     >     helm_only, opal_scvh_only, opal_only, scvh_only,
!     >		include_radiation,logT,logRho,temp,den,abar,zbar,X,Z,
!     >		pout,dpoutdd,dpoutdt,
!     >		eout,deoutdd,gamma3,
!     >		sout,dsoutdd,dsoutdt,
!     >		xout,dxoutdd,dxoutdt,
!     >		mu_M_out,logNe_out,eta_ele_out,data_dir,ierr)
     
		use eos_def
		use helm
		use opal_scvh_driver
		use utils_lib, only: is_bad_num
		use const_def, only: crad, ln10
		
!..logT and logRho are log10's of temp and den
      implicit none
      save

!..mixes the helmholtz and opal_scvh equation of state

		logical, intent(in) :: helm_only, opal_scvh_only, opal_only, 
     >		scvh_only, search_for_SCVH, include_radiation
      double precision, intent(in) :: logT,logRho,temp,den,abar,zbar,z53bar,X,Z
      double precision, intent(out) ::  
     >               logPgas, logE, logS, chiRho, chiT, 
     >               Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >               mu, free_e, gamma1, gamma3, grad_ad, eta, HELM_fraction
      character (len=*), intent(in) ::  data_dir
		integer, intent(out) :: ierr


!..define the boundary between helm and opal/scvh

      include "eos_regions_defs.dek"


      double precision :: logNe, 
     >   logPgas_helm, logE_helm, logS_helm, chiRho_helm, chiT_helm, 
     >   Cp_helm, Cv_helm, dE_dRho_helm, dS_dT_helm, dS_dRho_helm, 
     >   mu_helm, gamma1_helm, gamma3_helm, grad_ad_helm, logNe_helm, eta_helm, 
     >   logPgas_opalscvh, logE_opalscvh, logS_opalscvh, chiRho_opalscvh, chiT_opalscvh, 
     >   Cp_opalscvh, Cv_opalscvh, dE_dRho_opalscvh, dS_dT_opalscvh, dS_dRho_opalscvh, 
     >   mu_opalscvh, logNe_opalscvh, gamma1_opalscvh, gamma3_opalscvh, grad_ad_opalscvh,  
     >   eta_opalscvh
     
      integer iregion
      double precision :: pa,pb,ea,eb,sa,sb
		double precision :: a, b
		double precision, parameter :: pi = 3.1415926535897932384d0
		logical :: have_called_helm
		double precision, parameter :: helm_min_temp = 1d3

!..some physical constants
      double precision clight
      parameter        (clight  = 2.99792458d10)
      double precision avo
      parameter        (avo  = 6.0221367d23)


!..loading the opal_scvh tables
      integer          ifirst
      data             ifirst/0/


      double precision :: helm_res(num_helm_results), dlnPgas_dlnY, Pgas, Prad
      
      include 'formats'

		ierr = 0
		have_called_helm = .false.
		
		if (helm_only) then
			alfa = 1
		else if (opal_scvh_only .or. opal_only .or. scvh_only) then
			alfa = 0
		else
		   if (include_radiation) then
		      logT1 = logT1_default
		      logT2 = logT2_default
		   else
		      logT1 = logT1_no_rad_default
		      logT2 = logT2_no_rad_default
		   end if
         include 'eos_regions_code.dek'
		end if
		   
      beta = 1-alfa

      if (beta .ne. 0D0) then
         call interpolate_opal_scvh(
     >      opal_only, scvh_only, include_radiation, search_for_SCVH,
     >      logT, logRho, temp, den, abar, zbar, X, Z,
     >      logPgas_opalscvh, logE_opalscvh, logS_opalscvh, chiRho_opalscvh, chiT_opalscvh, 
     >      Cp_opalscvh, Cv_opalscvh, dE_dRho_opalscvh, dS_dT_opalscvh, dS_dRho_opalscvh, 
     >      mu_opalscvh, logNe_opalscvh, gamma1_opalscvh, gamma3_opalscvh, grad_ad_opalscvh,  
     >      eta_opalscvh, dlnPgas_dlnY,
     >      data_dir, ierr)
			
			Pgas = 10d0**logPgas_opalscvh
			Prad = crad*temp**4/3
			call check_results
			if (ierr /= 0) then
            if (scvh_only) then ! try opal
               ierr = 0
               call interpolate_opal_scvh(
     >            .true., .false., include_radiation, search_for_SCVH,
     >            logT, logRho, temp, den, abar, zbar, X, Z,
     >            logPgas_opalscvh, logE_opalscvh, logS_opalscvh, chiRho_opalscvh, chiT_opalscvh, 
     >            Cp_opalscvh, Cv_opalscvh, dE_dRho_opalscvh, dS_dT_opalscvh, dS_dRho_opalscvh, 
     >            mu_opalscvh, logNe_opalscvh, gamma1_opalscvh, gamma3_opalscvh, grad_ad_opalscvh,  
     >            eta_opalscvh, dlnPgas_dlnY,
     >            data_dir, ierr)
			      Pgas = 10d0**logPgas_opalscvh
			      Prad = crad*temp**4/3
			      call check_results
               if (ierr /= 0) then
                  write(*,*) 'bail to OPAL failed for logT logRho logQ', logT, logRho, logRho - 2*logT + 12
               end if
            end if
            if (ierr /= 0) then ! use helm instead
   			   beta = 0
   			   alfa = 1
            end if
			end if
			
      end if
      
      if (alfa .ne. 0D0) then
			call get_helmeos(include_radiation)
			if (ierr /= 0) then
			   write(*,*) 'failed in get_helmeos'
			   return
			end if
      end if
      
      if (alfa .eq. 0d0) then ! no HELM

         logPgas = logPgas_opalscvh
         logE = logE_opalscvh
         logS = logS_opalscvh
         chiRho = chiRho_opalscvh
         chiT = chiT_opalscvh
         Cp = Cp_opalscvh
         Cv = Cv_opalscvh
         dE_dRho = dE_dRho_opalscvh
         dS_dT = dS_dT_opalscvh
         dS_dRho = dS_dRho_opalscvh
         mu = mu_opalscvh
         logNe = logNe_opalscvh
         gamma1 = gamma1_opalscvh
         gamma3 = gamma3_opalscvh
         grad_ad = grad_ad_opalscvh
         eta = eta_opalscvh
         
      else if (beta .eq. 0d0) then ! pure HELM

         logPgas = logPgas_helm
         logE = logE_helm
         logS = logS_helm
         chiRho = chiRho_helm
         chiT = chiT_helm
         Cp = Cp_helm
         Cv = Cv_helm
         dE_dRho = dE_dRho_helm
         dS_dT = dS_dT_helm
         dS_dRho = dS_dRho_helm
         mu = mu_helm
         logNe = logNe_helm
         gamma1 = gamma1_helm
         gamma3 = gamma3_helm
         grad_ad = grad_ad_helm
         eta = eta_helm
      
      else ! combine alfa * helm + beta * opalscvh
      
         if (.false.) then
            logPgas = alfa*logPgas_helm + beta*logPgas_opalscvh
            logE = alfa*logE_helm + beta*logE_opalscvh
            logS = alfa*logS_helm + beta*logS_opalscvh
            chiRho = alfa*chiRho_helm + beta*chiRho_opalscvh
            chiT = alfa*chiT_helm + beta*chiT_opalscvh
            Cp = alfa*Cp_helm + beta*Cp_opalscvh
            Cv = alfa*Cv_helm + beta*Cv_opalscvh
            dE_dRho = alfa*dE_dRho_helm + beta*dE_dRho_opalscvh
            dS_dT = alfa*dS_dT_helm + beta*dS_dT_opalscvh
            dS_dRho = alfa*dS_dRho_helm + beta*dS_dRho_opalscvh
            mu = alfa*mu_helm + beta*mu_opalscvh
            gamma1 = alfa*gamma1_helm + beta*gamma1_opalscvh
            gamma3 = alfa*gamma3_helm + beta*gamma3_opalscvh
            grad_ad = alfa*grad_ad_helm + beta*grad_ad_opalscvh
            logNe = alfa*logNe_helm + beta*logNe_opalscvh
         
         else
         
            call blend(
     >         alfa, beta, den, temp, Prad, 
     >         logPgas_helm, logPgas_opalscvh, 
     >         logS_helm, logS_opalscvh, dS_dT_helm, dS_dT_opalscvh, dS_dRho_helm, dS_dRho_opalscvh,
     >         chiT_helm, chiT_opalscvh, chiRho_helm, chiRho_opalscvh, mu_helm, mu_opalscvh, logNe_helm, logNe_opalscvh,
     >         logE_helm, logE_opalscvh, Cv_helm, Cv_opalscvh, dE_dRho_helm, dE_dRho_opalscvh,
     >         gamma1_helm, gamma1_opalscvh, gamma3_helm, gamma3_opalscvh, grad_ad_helm, grad_ad_opalscvh,
     >         logPgas, logE, logS, chiRho, chiT, Cp, Cv, dE_dRho, dS_dT, dS_dRho,
     >         mu, gamma1, gamma3, grad_ad, logNe)
            
         end if
         
         eta = alfa*eta_helm + beta*eta_opalscvh

      end if
		
      HELM_fraction = alfa
		eta = min(200d0, eta)
		free_e = 10**logNe / (avo * den) ! convert to mean number of free electrons per nucleon     
		
      if (is_bad_num(logPgas) .or. is_bad_num(logE) .or. is_bad_num(logS) .or.
     >      is_bad_num(chiRho) .or. is_bad_num(chiT) .or. is_bad_num(Cp) .or. 
     >      is_bad_num(gamma1) .or. is_bad_num(gamma3) .or. is_bad_num(grad_ad) .or. 
     >      logPgas < -50d0) then
         write(*,1) 'HELM_fraction', HELM_fraction
         write(*,1) 'logT', logT
         write(*,1) 'logRho', logRho
         write(*,1) 'logQ', logRho - 2*logT + 12
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'X', X
         write(*,1) 'Z', Z
         write(*,*)
         write(*,1) 'logPgas', logPgas
         write(*,1) 'logE', logE
         write(*,1) 'logS', logS
         write(*,1) 'chiRho', chiRho
         write(*,1) 'chiT', chiT
         write(*,1) 'Cp', Cp
         write(*,1) 'gamma1', gamma1
         write(*,1) 'gamma3', gamma3
         write(*,1) 'grad_ad', grad_ad
         write(*,*)
   		write(*,1) 'alfa', alfa
   		write(*,1) 'beta', beta
   		write(*,*) 'helm_opal_scvh'
         stop 1
      end if

		
		contains

      subroutine check_results
         include 'formats'
         if (is_bad_num(logPgas_opalscvh) .or. 
     >       is_bad_num(logE_opalscvh) .or. 
     >       is_bad_num(logS_opalscvh) .or. 
     >       is_bad_num(chiRho_opalscvh) .or. 
     >       is_bad_num(chiT_opalscvh) .or. 
     >       is_bad_num(Cp_opalscvh) .or. 
     >       is_bad_num(Cv_opalscvh) .or. 
     >       is_bad_num(dE_dRho_opalscvh) .or. 
     >       is_bad_num(dS_dT_opalscvh) .or. 
     >       is_bad_num(dS_dRho_opalscvh) .or. 
     >       is_bad_num(mu_opalscvh) .or. 
     >       is_bad_num(logNe_opalscvh) .or. 
     >       is_bad_num(gamma1_opalscvh) .or. 
     >       is_bad_num(gamma3_opalscvh) .or. 
     >       is_bad_num(grad_ad_opalscvh) .or. 
     >       is_bad_num(eta_opalscvh) .or. 
     >       Cv_opalscvh <= 0d0 .or. 
     >       logPgas_opalscvh <= -50d0 .or. 
     >       logE_opalscvh <= -50d0 .or. 
     >       logS_opalscvh <= -50d0 .or. 
     >       logPgas_opalscvh <= -50d0 .or. 
     >       gamma1_opalscvh <= 1d0 .or. 
     >       gamma1_opalscvh >= 2d0 .or. 
     >       grad_ad_opalscvh >= 50d0) then
            ierr = -1
            !write(*,*) 'bail to HELM for logT logRho logQ', logT, logRho, logRho - 2*logT + 12
            if (.false.) then
               write(*,1) 'logT', logT
               write(*,1) 'logRho', logRho
               write(*,1) 'logQ', logRho - 2*logT + 12
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'X', X
               write(*,1) 'Z', Z
               write(*,*)
               stop 1
            end if
         end if
      end subroutine check_results

		
		subroutine get_helmeos(include_radiation)
		   logical, intent(in) :: include_radiation
		   
		   logical, parameter :: clip_to_table_boundaries = .false.
		   logical, parameter :: always_skip_elec_pos = .false.
		   logical, parameter :: always_include_elec_pos = .false.
		   double precision, parameter :: logT_ion = 5d0
		   double precision, parameter :: logT_neutral = 4.5d0
         logical :: off_table
		   
		   include 'formats'
		   
			if (have_called_helm) return
			have_called_helm = .true.
     
         call helmeos2(temp, logT, den, logRho, X, abar, zbar,
     >            1d6, 1d3, helm_res, clip_to_table_boundaries, .true., 
     >            always_skip_elec_pos, always_include_elec_pos, 
     >            logT_ion, logT_neutral, off_table, ierr)
			if (ierr /= 0) then
			   write(*,*) 'failed in helmeos2'
			   write(*,1) 'temp', temp
			   write(*,1) 'logT', logT
			   write(*,1) 'den', den
			   write(*,1) 'logRho', logRho
			   write(*,1) 'Z', Z
			   write(*,1) 'X', X
			   write(*,1) 'abar', abar
			   write(*,1) 'zbar', zbar
			   write(*,*) 'clip_to_table_boundaries', clip_to_table_boundaries
			   write(*,*) 'include_radiation', include_radiation
				if (.not. include_radiation) then
				   write(*,*) 'try it with radiation included'
				   ierr = 0
               call helmeos2(temp, logT, den, logRho, X, abar, zbar,
     >            1d6, 1d3, helm_res, clip_to_table_boundaries, .true., 
     >            always_skip_elec_pos, always_include_elec_pos, 
     >            logT_ion, logT_neutral, off_table, ierr)
               write(*,2) 'ierr with radiation', ierr
				end if
				write(*,*) 'stop in get_helmeos'
			   stop 1
			end if

			logPgas_helm = log10(helm_res(h_pgas))
			logE_helm = log10(helm_res(h_etot))
			logS_helm = log10(helm_res(h_stot))
			gamma1_helm = helm_res(h_gam1)
         gamma3_helm = helm_res(h_gam3)
			grad_ad_helm = helm_res(h_nabad)
			dE_dRho_helm = helm_res(h_ded)
			dS_dRho_helm = helm_res(h_dsd)
			dS_dT_helm = helm_res(h_dst)
			chiRho_helm = helm_res(h_chid)
			chiT_helm = helm_res(h_chit)
			Cv_helm = helm_res(h_cv)
			Cp_helm = helm_res(h_cp)

			if (helm_res(h_xne) > 0d0) then
				logNe_helm = log10(helm_res(h_xne)) ! assuming complete ionization
			else
				logNe_helm = -99d0
			end if
			
			mu_helm = abar / (1 + zbar)
			eta_helm = helm_res(h_etaele)
			
			if (is_bad_num(logE_helm)) then
			   write(*,1) 'helm_res(h_etot)', helm_res(h_etot)
				write(*,*) 'stop in get_helmeos'
			   stop 1
			end if
			
		end subroutine get_helmeos

      end subroutine


      end module helm_opal_scvh_driver
