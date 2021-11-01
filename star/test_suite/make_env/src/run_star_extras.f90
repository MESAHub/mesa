! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none
      
      include "test_suite_extras_def.inc"      
            
      
      contains

      include "test_suite_extras.inc"

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         
         if (.not. s% x_logical_ctrl(1)) return
         
         call create_env(id, s, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in create_env')

      end subroutine extras_controls
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         include 'formats'
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step
      
      
      subroutine create_env(id, s, ierr)
         use eos_lib
         use eos_def, only: i_lnfree_e, num_eos_basic_results, num_eos_d_dxa_results
         use chem_lib, only: basic_composition_info
         use utils_lib, only: is_bad
         use atm_lib, only: atm_Teff
         
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k, nz, species, iter, max_iters, &
            i_lnd, i_lnR, i_lnT, i_lum
         character (len=256):: net_name, save_atm_option, &
            save_atm_T_tau_relation, save_atm_T_tau_opacity
         logical :: skip_partials
         real(dp) :: &
            ln_dq1, dq1, dq_factor, sumx, tau_surf, Teff, area, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            T_surf, P_surf, T_m1, P_m1, T_00, P_00, logRho, logT, &
            save_Pextra_factor, gr_factor, d_gr_factor_dlnR00, frac_Type2, &
            kap, dlnkap_dlnd, dlnkap_dlnT, grav, rho_00, vol, &
            P_face, Pgas, Prad, r_00, r_p1, rho1_00, logRho_m1, &
            T_face, d_T_face_dT, dlnP_face, dm_face, dT, &
            gradT, d_gradT_dT, resid, d_resid_dT, &
            T_expected, d_T_expected_dT, tau_base, T4, Teff4, r_phot
         real(dp) :: res(num_eos_basic_results)
         real(dp) :: dres_dlnRho(num_eos_basic_results)
         real(dp) :: dres_dlnT(num_eos_basic_results)
         real(dp), allocatable :: dres_dxa(:,:)
         real(dp), parameter :: LOGRHO_TOL = 1d-11
         real(dp), parameter :: LOGPGAS_TOL = 1d-11
         
         include 'formats'
         
         ierr = 0
         
         nz = s% x_integer_ctrl(1)
         s% nz = nz     
         max_iters = 100
          
         net_name = s% x_character_ctrl(1)
         s% mstar = s% x_ctrl(1)*Msun
         s% xmstar = s% x_ctrl(2)*s% mstar
         s% M_center = s% mstar - s% xmstar
         s% star_mass = s% mstar/Msun

         call star_set_net(id, net_name, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in star_set_net')    
         
         call star_set_var_info(id, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in star_set_var_info')     
         
         call star_set_chem_names(id, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in star_set_chem_names')  
         
         call star_allocate_arrays(id, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in star_allocate_arrays') 
         
         s% m(1) = s% mstar
         s% m_grav(1) = s% mstar
         s% cgrav(1) = standard_cgrav
         s% tau_factor = s% x_ctrl(6)

         species = s% species

         allocate(dres_dxa(num_eos_d_dxa_results, species))
         
         if (s% x_logical_ctrl(2)) then ! R and L in cgs units
            s% r(1) = s% x_ctrl(3)
            s% L(1:nz) = s% x_ctrl(4)
         else ! x_ctrl(3) = Teff; x_ctrl(4) = L/Lsun
            s% L(1:nz) = s% x_ctrl(4)*Lsun
            ! L = 4*pi*R^2*boltz_sigma*Tsurf**4
            tau_base = two_thirds ! just use Eddington for this
            tau_surf = s% tau_factor*tau_base
            Teff = s% x_ctrl(3)
            Teff4 = pow4(Teff)
            T4 = Teff4*0.75d0*(tau_surf + tau_base)
            s% r(1) = sqrt(s% L(1)/(4d0*pi*boltz_sigma*T4))
            r_phot = sqrt(s% L(1)/(4d0*pi*boltz_sigma*Teff4))
            if (abs(r_phot - s% r(1)) > 1d-3*r_phot) then
               write(*,1) 's% r(1)', s% r(1)
               write(*,1) 'r_phot', r_phot
               write(*,1) 's% r(1)/r_phot', s% r(1)/r_phot
            end if
         end if
         
         i_lnd = s% i_lnd
         i_lnR = s% i_lnR
         i_lnT = s% i_lnT
         i_lum = s% i_lum

         s% L_center = s% L(nz)
         s% r_start(1) = s% r(1)        
         
         ln_dq1 = s% x_ctrl(5)*ln10
         dq1 = exp(ln_dq1)
         dq_factor = calc_dq_factor(nz,dq1)
                  
         s% q(1) = 1d0
         s% dq(1) = dq1
         do k=2, nz
            s% q(k) = s% q(k-1) - s% dq(k-1)
            s% dq(k) = s% dq(k-1)*dq_factor
         end do

         call star_normalize_dqs(s% id, nz, s% dq, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in star_normalize_dqs')  

         call star_set_qs(s% id, nz, s% q, s% dq, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in star_set_qs')  
         
         call star_set_m_and_dm(s% id, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in star_set_qs')  
         
         call star_set_dm_bar(s% id, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in star_set_dm_bar') 
         
         call change_to_xa_for_accretion(s% id, 1, nz, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in change_to_xa_for_accretion') 

         do k=1,nz
            s% m_grav(k) = s% m(k)
            call basic_composition_info( &
               species, s% chem_id, s% xa(1:species,k), s% X(k), s% Y(k), s% Z(k), &
               s% abar(k), s% zbar(k), s% z2bar(k), s% z53bar(k), s% ye(k), &
               s% mass_correction(k), sumx)
         end do

         s% tau_base = two_thirds
         tau_surf = s% tau_factor*s% tau_base
         save_atm_option = s% atm_option
         save_atm_T_tau_relation = s% atm_T_tau_relation
         save_atm_T_tau_opacity = s% atm_T_tau_opacity
         save_Pextra_factor = s% Pextra_factor
         s% atm_option = 'T_tau'
         s% atm_T_tau_relation = 'Eddington'
         s% atm_T_tau_opacity = 'iterated'
         s% Pextra_factor = 2
         s% Teff = atm_Teff(s% L(1), s% r(1))
         call get_initial_guess_for_atm(ierr)
         if (ierr /= 0) then
            write(*, *) 'Call get_initial_guess_for_atm failed', k
            stop
         end if

         if (.not. s% x_logical_ctrl(2) .and. &
             abs(Teff - s% Teff) > 1d-5*Teff) then ! report Teff
            write(*,1) 'desired Teff', Teff
            write(*,1) '1st actual s% Teff', s% Teff
            write(*,1) 'old r(1)', s% r(1)
            s% r(1) = s% r(1)*pow2(s% Teff/Teff)
            write(*,1) 'new r(1)', s% r(1)
         end if
         
         ! do cell k=1 and then redo the atm.
         T_m1 = T_surf
         P_m1 = P_surf
         logRho_m1 = logRho
         k = 1
         call do1_cell(ierr)
         if (ierr /= 0) then
            write(*, *) 'Call do1_cell initial failed', k
            stop
         end if
         s% atm_option = save_atm_option
         s% atm_T_tau_relation = save_atm_T_tau_relation
         s% atm_T_tau_opacity = save_atm_T_tau_opacity
         s% Pextra_factor = save_Pextra_factor
         call get_atm(ierr)
         if (ierr /= 0) then
            write(*, *) 'Call get_atm initial failed', k
            stop
         end if

         if (.not. s% x_logical_ctrl(2) .and. &
             abs(Teff - s% Teff) > 1d-5*Teff) then ! report Teff
            write(*,1) 'final actual s% Teff', s% Teff
         end if

         ! update info for surface for final atm
         T_m1 = T_surf
         P_m1 = P_surf
         logRho_m1 = logRho
         
         do k=1,nz
            call do1_cell(ierr)
            if (ierr /= 0) then
               write(*, *) 'Call do1_cell failed', k
               stop
            end if
            T_m1 = T_00
            P_m1 = P_00
            logRho_m1 = logRho
         end do
         
         s% model_number = 0
         s% star_age = 0
                     
         write(*,2) 'start.mod', nz
         call star_write_model(id, 'start.mod', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in star_write_model') 

         write(*,'(A)')
         write(*,*) 'finished create_env'
         write(*,'(A)')
         !stop

         deallocate(dres_dxa)
         
         contains

         
         subroutine get_initial_guess_for_atm(ierr)
            integer, intent(out) :: ierr
            skip_partials = .true.
            s% opacity(1) = 1d-2 ! kap_guess
            call star_get_atm_PT( &
                s% id, tau_surf, s% L(1), s% r(1), s% m(1), s% cgrav(1), skip_partials, &
                s% Teff, &
                lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
                lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
                ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in get_atm_PT')
            T_surf = exp(lnT_surf)
            P_surf = exp(lnP_surf)
            ! get rho_surf
            Prad = one_third*crad*pow4(T_surf)
            Pgas = P_surf - Prad
            logT = log10(T_surf)
            k = 1
            call star_solve_eos_given_PgasT_auto( &
               s% id, k, s% xa(:,k), &
               logT, log10(Pgas), LOGRHO_TOL, LOGPGAS_TOL, &
               logRho, res, dres_dlnRho, dres_dlnT, dres_dxa, &
               ierr)
            if (ierr /= 0) then
               write(*, *) 'Call star_solve_eos_given_PgasT_auto failed', k
               stop
            end if
         end subroutine get_initial_guess_for_atm

         
         subroutine get_atm(ierr)
            integer, intent(out) :: ierr
            logical, parameter :: &
               need_atm_Psurf = .true., need_atm_Tsurf = .true.
            include 'formats'
            ierr = 0
            skip_partials = .true.
            if (s% opacity(1) <= 0d0 .or. is_bad(s% opacity(1))) then
               write(*,1) 's% opacity(1)', s% opacity(1)
               call mesa_error(__FILE__,__LINE__,'run_star_extras get_atm')
            end if
            s% opacity_start(1) = s% opacity(1)
            call star_get_surf_PT( &
               s% id, skip_partials, &
               need_atm_Psurf, need_atm_Tsurf, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'get_atm failed in star_get_surf_PT')
            T_surf = exp(lnT_surf)
            P_surf = exp(lnP_surf)
            ! get rho_surf
            Prad = one_third*crad*pow4(T_surf)
            Pgas = P_surf - Prad
            logT = log10(T_surf)
            k = 1
            call star_solve_eos_given_PgasT( &
               s% id, k, s% xa(:,k), &
               logT, log10(Pgas), logRho, LOGRHO_TOL, LOGPGAS_TOL, &
               logRho, res, dres_dlnRho, dres_dlnT, dres_dxa, &
               ierr)
            if (ierr /= 0) then
               write(*, *) 'Call star_solve_eos_given_PgasT failed in get_atm'
               stop
            end if
         end subroutine get_atm
         
         
         subroutine do1_cell(ierr) ! uses r(k), T_m1, P_m1, logRho_m1
            use eos_def, only: i_grad_ad
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            r_00 = s% r(k)
            s% cgrav(k) = standard_cgrav
            grav = -s% cgrav(k)*s% m(k)/r_00**2
            
            if (k > 1) then
               dm_face = 0.5d0*(s% dm(k) + s% dm(k-1))
            else ! k == 1
               dm_face = 0.5d0*s% dm(k)
            end if
            area = 4d0*pi*r_00**2
            
            P_00 = P_m1 - grav*dm_face/area
            P_face = 0.5d0*(P_m1 + P_00)
            dlnP_face = (P_m1 - P_00)/P_face
            T_00 = T_m1 ! initialize for loop
            do iter = 1,max_iters ! solve implicit eqn for T_00
               Prad = one_third*crad*pow4(T_00)
               Pgas = P_00 - Prad
               logT = log10(T_00)
               
               call star_solve_eos_given_PgasT( &
                  s% id, k, s% xa(:,k), &
                  logT, log10(Pgas), logRho_m1, LOGRHO_TOL, LOGPGAS_TOL, &
                  logRho, res, dres_dlnRho, dres_dlnT, dres_dxa, &
                  ierr)
               if (ierr /= 0) then
                  write(*, *) 'Call star_solve_eos_given_PgasT failed', k
                  stop
               end if
               s% lnd(k) = logRho*ln10
               s% lnT(k) = logT*ln10
               s% rho(k) = exp10(logRho)
               s% T(k) = T_00
               call star_do_eos_for_cell(s% id, k, ierr)
               if (ierr /= 0) then
                  write(*, *) 'Call star_do_eos_for_cell failed', k
                  stop
               end if
               
               s% extra_opacity_factor(k) = 1d0
               call star_do_kap_for_cell(s% id, k, ierr)
               if (ierr /= 0) then
                  write(*, *) 'Call star_do_kap_for_cell failed', k
                  stop
               end if

               ! need these quantities to be set in MLT
               s% lnT_start(k) = s% lnT(k)
               s% csound_face(k) = s% csound(k)
               s% csound_start(k) = s% csound(k)
               s% mlt_gradT_fraction = -1d0
               s% adjust_mlt_gradT_fraction(k) = -1d0
               
               ! skipping use_other_alpha_mlt and other_gradr_factor
               s% alpha_mlt(k) = s% mixing_length_alpha
               s% gradr_factor(k) = 1d0               
               call star_set_mlt_vars(s% id, k, k, ierr)
               if (ierr /= 0) then
                  write(*, *) 'Call set_mlt_vars failed', k
                  stop
               end if
               
               gradT = s% gradT(k)
               d_gradT_dT = s% gradT_ad(k)%d1Array(i_lnT_00)/T_00
               T_face = 0.5d0*(T_m1 + T_00)
               d_T_face_dT = 0.5d0
               T_expected = T_m1 - T_face*gradT*dlnP_face
               d_T_expected_dT = -dlnP_face*(d_T_face_dT*gradT + T_face*d_gradT_dT)
               resid = T_expected - T_00
               d_resid_dT = d_T_expected_dT - 1d0
               dT = -resid/d_resid_dT
               if (abs(dT/T_00) < 1d-4) exit
               if (iter == max_iters) then
                  write(*, *) 'Failed to find T', k
                  stop
               end if
               T_00 = T_00 + dT
            end do
            
            rho_00 = exp10(logRho)
            vol = four_thirds*pi*pow3(r_00) - s% dm(k)/rho_00
            r_p1 = pow(0.75d0*vol/pi, one_third)
            if (k == nz) then
               s% R_center = r_p1
            else
               s% r(k+1) = r_p1
            end if
            
            s% lnR(k) = log(r_00)
            s% xh(i_lnR, k) = s% lnR(k)
            s% lnd(k) = log(rho_00)
            s% xh(i_lnd, k) = s% lnd(k)
            s% lnT(k) = log(T_00)
            s% xh(i_lnT, k) = s% lnT(k)
            s% xh(i_lum, k) = s% L(k)
         
         end subroutine do1_cell
         
         
         real(dp) function dq_f(r, dfdr, lrpar, rpar, lipar, ipar, ierr)
            ! returns with ierr = 0 if was able to evaluate f and df/dx at x
            ! if df/dx not available, it is okay to set it to 0
            use const_def, only: dp
            integer, intent(in) :: lrpar, lipar
            real(dp), intent(in) :: r
            real(dp), intent(out) :: dfdr
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr
            integer :: n
            real(dp) :: r_to_n, a1, Sn
            include 'formats'
            n = ipar(1)
            a1 = rpar(1)
            Sn = 1d0
            r_to_n = pow(r,n)
            dq_f = Sn*(1d0 - r) - a1*(1d0 - r_to_n)
            dfdr = -Sn + a1*n*r_to_n/r
            ipar(2) = ipar(2) + 1
            !write(*,*) ipar(2), 'r, dq_f, dfdr', r, dq_f, dfdr
         end function dq_f
         
         real(dp) function calc_dq_factor(n,dq1) result(dq_factor)
            use num_lib, only: safe_root_with_guess
            integer, intent(in) :: n
            real(dp), intent(in) :: dq1
            real(dp) :: x_guess, dx, x1, x3, y1, y3, epsx, epsy
            integer :: newt_imax, imax
            integer, parameter :: lipar = 2, lrpar = 1
            integer, target :: ipar_array(lipar)
            integer, pointer :: ipar(:) ! (lipar)
            real(dp), target :: rpar_array(lrpar)
            real(dp), pointer :: rpar(:) ! (lrpar)            
            include 'formats'
            ierr = 0
            ipar => ipar_array
            rpar => rpar_array
            ipar(1) = n
            ipar(2) = 0
            rpar(1) = dq1
            x_guess = 1.08d0
            dx = 0d0 ! not used since give x1 and x3
            x1 = 1.01d0
            x3 = 1.05d0
            y1 = arg_not_provided
            y3 = arg_not_provided
            newt_imax = 0 ! 10
            imax = 100
            epsx = 1d-10
            epsy = 1d-10
            dq_factor = safe_root_with_guess( &
               dq_f, x_guess, dx, x1, x3, y1, y3, newt_imax, imax, &
               epsx, epsy, lrpar, rpar, lipar, ipar, ierr)
            !write(*,*) 'dq_factor', dq_factor
            !stop
         end function calc_dq_factor
         
      end subroutine create_env
      

      end module run_star_extras
      
