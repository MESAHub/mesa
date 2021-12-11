! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
      use auto_diff
      use chem_def
      use utils_lib
      use rates_def, only: i_rate

      use interp_1d_def, only: pm_work_size
      use interp_1d_lib, only: interp_pm, interp_values, interp_value
      
      implicit none
      
      include "test_suite_extras_def.inc"

      logical :: dbg = .false.

      logical :: pgstar_flag
      integer :: pgstar_interval, pgstar_file_interval

      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! here are definitions for the indices of all s% xtra, s% lxtra and s% ixtra variables used in the run
      !!!!!!!!!!!!!!!!!!!!!!!!!

      ! lxtra(lx_hydro_on) is true while hydro is on
      integer, parameter :: lx_hydro_on = 1
      ! lxtra(lx_hydro_has_been_on) is true after hydro is turned on for the first time
      integer, parameter :: lx_hydro_has_been_on = 2
      ! lxtra(lx_have_saved_gamma_prof) is true, after a profile is saved when gamma_integral=0
      ! One profile is saved during each hydro phase, this is turned
      ! back to false after a relax
      integer, parameter :: lx_have_saved_gamma_prof = 3
      ! lx_have_reached_gamma_limit is same as before, but it's kept as .true. after the first time
      ! gamma_integral=0. Used to count time from first pulse.
      integer, parameter :: lx_have_reached_gamma_limit = 4
      ! lxtra(lx_he_zams) indicates if star has reached he zams
      integer, parameter :: lx_he_zams = 5

      ! xtra(x_time_start_pulse) contains the time at which a pulse starts
      integer, parameter :: x_time_start_pulse = 1
      ! xtra(x_dyn_time) contains the estimate for the dynamical time when that happens
      integer, parameter :: x_dyn_time = 2
      ! xtra(x_gamma_int_bound) contains gamma_integral in bound regions
      integer, parameter :: x_gamma_int_bound = 3
      ! xtra(x_time_since_first_gamma_zero) contains time in seconds since first time gamma_integral=0
      integer, parameter :: x_time_since_first_gamma_zero = 4
      ! xtra(x_star_age_at_relax) Stores s% star_age at the point the last relax was made
      integer, parameter :: x_star_age_at_relax = 8

      ! ixtra(ix_steps_met_relax_cond) counts the steps during which the conditions to relax a model are met
      integer, parameter :: ix_steps_met_relax_cond = 1
      ! ixtra(ix_num_relaxations) counts the number of times the star has relaxed
      integer, parameter :: ix_num_relaxations = 2
      ! ixtra(ix_steps_since_relax) counts the number of steps since last relax
      integer, parameter :: ix_steps_since_relax = 3
      ! ixtra(ix_steps_since_hydro_on) counts the number of steps since hydro was turned on
      integer, parameter :: ix_steps_since_hydro_on = 4

      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! These variables are loaded up from x_ctrl, x_integer_ctrl and x_logical_ctrl
      ! values specified on inlist_ppisn
      !!!!!!!!!!!!!!!!!!!!!!!!!

      real(dp) :: min_gamma_sub_43_for_hydro
      real(dp) :: max_v_for_pulse, q_for_dyn_ts, num_dyn_ts_for_relax
      real(dp) :: q_for_relax_check, max_v_for_relax, max_machn_for_relax, &
         max_Lneu_for_relax, max_Lnuc_for_relax
      integer :: num_steps_before_relax
      !logical :: remove_extended_layers, in_inlist_pulses
      logical :: in_inlist_pulses
      real(dp) :: max_dt_before_pulse
      real(dp) :: max_Lneu_for_mass_loss
      real(dp) :: delta_lgLnuc_limit, max_Lphoto_for_lgLnuc_limit, max_Lphoto_for_lgLnuc_limit2
      real(dp) :: delta_lgRho_cntr_hard_limit, dt_div_min_dr_div_cs_limit
      real(dp) :: logT_for_v_flag, logLneu_for_v_flag
      logical :: stop_100d_after_pulse
      
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
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         ! we turn pgstar on and off at some points, so we store the original setting
         pgstar_flag = s% job% pgstar_flag
         pgstar_interval = s% pgstar_interval
         pgstar_file_interval = s% Grid2_file_interval

         ! this is optional
         s% other_wind => brott_wind
         s% other_adjust_mdot => my_adjust_mdot
         s% other_before_struct_burn_mix => my_before_struct_burn_mix
         s% other_eval_fp_ft => my_other_eval_fp_ft
         s% other_kap_get => my_other_kap_get

         ! store user provided options from the inlist
         min_gamma_sub_43_for_hydro = s% x_ctrl(1)
         max_v_for_pulse = s% x_ctrl(2)
         q_for_dyn_ts = s% x_ctrl(3)
         num_dyn_ts_for_relax = s% x_ctrl(4)
         q_for_relax_check = s% x_ctrl(5)
         max_v_for_relax = s% x_ctrl(6)
         max_machn_for_relax = s% x_ctrl(7)
         max_Lneu_for_relax = s% x_ctrl(8)
         max_Lnuc_for_relax = s% x_ctrl(9)
         num_steps_before_relax = s% x_integer_ctrl(1)
         in_inlist_pulses = s% x_logical_ctrl(2)
         max_dt_before_pulse = s% x_ctrl(10)
         max_Lneu_for_mass_loss = s% x_ctrl(11)
         delta_lgLnuc_limit = s% x_ctrl(12)
         max_Lphoto_for_lgLnuc_limit = s% x_ctrl(13)
         max_Lphoto_for_lgLnuc_limit2 = s% x_ctrl(14)
         logT_for_v_flag = s% x_ctrl(15)
         logLneu_for_v_flag = s% x_ctrl(16)
         stop_100d_after_pulse = s% x_logical_ctrl(1)

         ! we store the value given in inlist_ppisn and deactivate it at
         ! high T
         delta_lgRho_cntr_hard_limit = s% delta_lgRho_cntr_hard_limit
         ! we also store dt_div_min_dr_div_cs_limit, we keep it at a
         ! high value until the onset of a pulse to prevent unnecesarily
         ! small timesteps before a pulsation
         dt_div_min_dr_div_cs_limit = s% dt_div_min_dr_div_cs_limit

      end subroutine extras_controls

      subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         integer :: h1, he4
         real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
            vink_wind, nieu_wind, hamann_wind
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         L1 = Lsurf
         M1 = Msurf
         R1 = Rsurf
         T1 = Tsurf

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         Xs = s% xa(h1,1)
         Ys = s% xa(he4,1)
         ! Z=0.0142 is Z from Asplund et al. 2009
         Z_div_Z_solar = s% kap_rq% Zbase/0.0142d0
         ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
         Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z_div_Z_solar)))

         vink_wind = 0d0
         nieu_wind = 0d0
         hamann_wind = 0d0
         w = 0

         call eval_Vink_wind(vink_wind)
         call eval_Nieuwenhuijzen_wind(nieu_wind)
         call eval_Hamann_wind(hamann_wind)

         ! use 1/10 hamann
         hamann_wind = hamann_wind/10d0

         if (T1 < Teff_jump) then
            ! low T wind
            w = max(vink_wind, nieu_wind)
         else
            ! high T wind
            alfa = 0d0
            if (Xs > 0.7d0) then
               alfa = 1d0
            else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
               alfa = (Xs - 0.4d0)/0.3d0
            end if
            w = alfa * vink_wind + (1d0-alfa) * hamann_wind
         end if

         ierr = 0

         contains

         subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = (T1 - (Teff_jump - dT)) / (2*dT)
               end if
            end if
            
            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w1 = exp10(logMdot)
            else
               w1 = 0
            end if
            
            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10(L1/Lsun/1d5) &
                  - 1.339d0*log10(M1/Msun/30) &
                  - 1.601d0*log10(vinf_div_vesc/2d0) &
                  + 1.07d0*log10(T1/2d4) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w2 = exp10(logMdot)
            else
               w2 = 0
            end if
            
            w = alfa*w1 + (1 - alfa)*w2
            
         end subroutine eval_Vink_wind

         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 &
                     +1.24d0*log10(L1/Lsun) &
                     +0.16d0*log10(M1/Msun) &
                     +0.81d0*log10(R1/Rsun) &
                     +0.85d0*log10(Z_div_Z_solar)
            w = exp10(log10w)
         end subroutine eval_Nieuwenhuijzen_wind

         subroutine eval_Hamann_wind(w)
            ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -11.95d0 &
                     +1.5d0*log10(L1/Lsun) &
                     -2.85d0*Xs &
                     + 0.85d0*log10(Z_div_Z_solar)
            w = exp10(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind
      
      subroutine my_adjust_mdot(id, ierr)
         use star_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: Lrad_div_Ledd
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% generations > 2) then
            write(*,*) "check mdots", s% mstar_dot, s% mstar_dot_old
            if (abs(s% mstar_dot) > 1.05d0*abs(s% mstar_dot_old)) then
               s% mstar_dot = 1.05d0*s% mstar_dot_old
            else if (abs(s% mstar_dot) < 0.95d0*abs(s% mstar_dot_old)) then
               s% mstar_dot = 0.95d0*s% mstar_dot_old
            end if
         end if
      end subroutine my_adjust_mdot

      subroutine my_other_eval_fp_ft( &
            id, nz, xm, r, rho, aw, ft, fp, r_polar, r_equatorial, report_ierr, ierr)
         use num_lib
         integer, intent(in) :: id
         integer, intent(in) :: nz
         real(dp), intent(in) :: aw(:), r(:), rho(:), xm(:) ! (nz)
         type(auto_diff_real_star_order1), intent(out) :: ft(:), fp(:) ! (nz)
         real(dp), intent(inout) :: r_polar(:), r_equatorial(:) ! (nz)
         logical, intent(in) :: report_ierr
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: j
         real(dp) :: alpha

         type(auto_diff_real_1var_order1) :: A_omega,fp_numerator, ft_numerator, w, w2, w3, w4, w5, w6, lg_one_sub_w4, &
            fp_temp, ft_temp

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

!$OMP PARALLEL DO PRIVATE(j, A_omega, fp_numerator, ft_numerator, fp_temp, ft_temp, w, w2, w3, w4, w5, w6, lg_one_sub_w4) SCHEDULE(dynamic,2)
            do j=1, s% nz
               !Compute fp, ft, re and rp using fits to the Roche geometry of a single star.
               !by this point in the code, w_div_w_crit_roche is set
               w = s% w_div_w_crit_roche(j)
               w%d1val1 = 1d0

               w2 = pow2(w)
               w4 = pow4(w)
               w6 = pow6(w)
               lg_one_sub_w4 = log(1d0-w4)
               A_omega = (1d0-0.1076d0*w4-0.2336d0*w6-0.5583d0*lg_one_sub_w4)
               fp_numerator = (1d0-two_thirds*w2-0.06837d0*w4-0.2495d0*w6)
               ft_numerator = (1d0+0.2185d0*w4-0.1109d0*w6)
               !fits for fp, ft
               fp_temp = fp_numerator/A_omega
               ft_temp = ft_numerator/A_omega
               !re and rp can be derived analytically from w_div_wcrit
               r_equatorial(j) = r(j)*(1d0+w2% val/6d0-0.0002507d0*w4% val+0.06075d0*w6% val)
               r_polar(j) = r_equatorial(j)/(1d0+0.5d0*w2% val)
               ! Be sure they are consistent with r_Phi
               r_equatorial(j) = max(r_equatorial(j),r(j))
               r_polar(j) = min(r_polar(j),r(j))

               fp(j) = 0d0
               ft(j) = 0d0
               fp(j)% val = fp_temp% val
               ft(j)% val = ft_temp% val
               if (s% w_div_wc_flag) then
                  fp(j)% d1Array(i_w_div_wc_00) = fp_temp% d1val1
                  ft(j)% d1Array(i_w_div_wc_00) = ft_temp% d1val1
               end if
            end do
!$OMP END PARALLEL DO

         if (s% u_flag) then
            !make fp and ft 1 in the outer 0.001 mass fraction of the star. softly turn to zero from the outer 0.002
            do j=1, s% nz
               if (s% q(j) > 0.999) then
                  fp(j) = 1d0
                  ft(j) = 1d0
               else if (s% q(j) > 0.998) then
                  alpha = (1d0-(s% q(j)-0.998)/(0.001))
                  fp(j) = fp(j)*alpha + 1d0*(1-alpha)
                  ft(j) = ft(j)*alpha + 1d0*(1-alpha)
               end if
            end do
         end if

      end subroutine my_other_eval_fp_ft


      subroutine my_other_kap_get( &
            id, k, handle, species, chem_id, net_iso, xa, &
            log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)

         use kap_def, only: num_kap_fracs
         use kap_lib
 
         ! INPUT
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! from alloc_kap_handle
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         real(dp), intent(in) :: log10_rho ! density
         real(dp), intent(in) :: log10_T ! temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
            ! eta := electron degeneracy parameter

         ! OUTPUT
         real(dp), intent(out) :: kap_fracs(num_kap_fracs)
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
         real(dp), intent(out) :: dln_kap_dxa(:) ! partial derivative w.r.t. to species
         integer, intent(out) :: ierr ! 0 means AOK.

         type (star_info), pointer :: s
         real(dp) :: velocity
         real(dp) :: radius, logR
         real(dp) :: logT_alt, inv_diff
         real(dp) :: log_kap, alpha

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
                  
         kap = 0; dln_kap_dlnRho = 0; dln_kap_dlnT = 0; dln_kap_dxa = 0
         velocity = 0
         radius = 0

         !if (k==1 .and. s% u_flag .and. .not. is_nan(s% lnR_start(1))) then !very surface cell can go mad, things are more stable if we fix opacity
         !   if (s% xh_start(s% i_u,1)>sqrt(2*s% cgrav(1)*s% m(1)/exp(s% lnR_start(1)))) then
         if (k==1 .and. s% u_flag) then !very surface cell can go mad, things are more stable if we fix opacity
            ! this is to support restarts, as xh_start and r_start are
            ! not properly set when model loads
            if (s% solver_iter > 0) then
               velocity = s% xh_start(s% i_u,1)
               radius = s% r_start(1)
            else
               velocity = s% xh(s% i_u,1)
               radius = s% r(1)
            end if
            if (velocity>sqrt(2*s% cgrav(1)*s% m(1)/radius)) then
               kap = 0.2d0*(1 + s% X(1))
               dln_kap_dlnRho = 0d0
               dln_kap_dlnT = 0d0
               return
            else
               call kap_get( &
                  s% kap_handle, species, chem_id, net_iso, xa, &
                  log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  eta, d_eta_dlnRho, d_eta_dlnT, &
                  kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)
            end if
         else
            call kap_get( &
               s% kap_handle, species, chem_id, net_iso, xa, &
               log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               eta, d_eta_dlnRho, d_eta_dlnT, &
               kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)
         end if


      end subroutine my_other_kap_get

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         if (.not. restart) then
            s% lxtra(lx_hydro_on) = .false.
            s% lxtra(lx_hydro_has_been_on) = .false.
            s% lxtra(lx_have_saved_gamma_prof) = .false.
            s% lxtra(lx_have_reached_gamma_limit) = .false.
            s% lxtra(lx_he_zams) = .false.
            s% xtra(x_time_start_pulse) = -1d0
            s% xtra(x_dyn_time) = -1d0
            s% xtra(x_gamma_int_bound) = -1d0
            s% xtra(x_time_since_first_gamma_zero) = 0d0
            s% ixtra(ix_steps_met_relax_cond) = 0
            s% ixtra(ix_num_relaxations) = 0
            s% ixtra(ix_steps_since_relax) = 0
            s% ixtra(ix_steps_since_hydro_on) = 0
            s% xtra(x_star_age_at_relax) = -1d0
            ! to avoid showing pgstar stuff during initial model creation
            s% pgstar_interval = 100000000
            s% Grid2_file_interval = 100000000
         else ! it is a restart
            if (s% lxtra(lx_hydro_on)) then
               call star_read_controls(id, 'inlist_hydro_on', ierr)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return
               call star_set_u_flag(id, .true., ierr)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return
            else if (s% lxtra(lx_hydro_has_been_on)) then
               call star_read_controls(id, 'inlist_hydro_off', ierr)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return
            end if
            if (s% lxtra(lx_hydro_has_been_on)) then
               call star_read_controls(id, 'inlist_after_first_pulse', ierr)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return
            end if
         end if
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         character (len=strlen) :: test
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr, k
         real(dp) :: max_v
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         

      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 22
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n), v_esc
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, k0
         integer, parameter :: h_relax_count = 1, &
                               h_log_R_098 = 2, &
                               h_log_R_095 = 3, &
                               h_log_R_090 = 4, &
                               h_helium_core_mass = 5, &
                               h_co_core_mass = 6, &
                               h_log_R_100 = 7, &
                               h_log_R_vesc = 8, &
                               h_log_R_vesc_098 = 9, &
                               h_log_R_vesc_095 = 10, &
                               h_log_R_vesc_090 = 11, &
                               h_M_below_vesc = 12, &
                               h_u_flag = 13, &
                               h_U_ejecta = 14, &
                               h_T_ejecta = 15, &
                               h_Omega_ejecta = 16, &
                               h_gamma_integral = 17, &
                               h_max_velocity = 18, &
                               h_min_velocity = 19, &
                               h_yr_since_coll = 20, &
                               h_total_J = 21, &
                               h_total_J_bel_vesc = 22

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(h_relax_count) = "relax_count"
         vals(h_relax_count) = s% ixtra(ix_num_relaxations)

         names(h_log_R_098) = "log_R_098" ! Radius of 98% of the mass of the star
         names(h_log_R_095) = "log_R_095" ! of 95%
         names(h_log_R_090) = "log_R_090" ! and 90%
         names(h_helium_core_mass) = "helium_core_mass"
         names(h_co_core_mass) = "co_core_mass"
         names(h_log_R_100) = "log_R_100" ! Radius of outermost layer
         names(h_log_R_vesc) = "log_R_vesc" ! Radius of layers just below the escape velocity
         names(h_log_R_vesc_098) = "log_R_vesc_098" ! Radius of 98% of the mass just below the esc. vel.
         names(h_log_R_vesc_095) = "log_R_vesc_095" ! of 95%
         names(h_log_R_vesc_090) = "log_R_vesc_090" ! and 90%
         names(h_M_below_vesc) = "M_below_vesc" ! Mass of layers below vesc
         names(h_u_flag) = "u_flag" ! 1 if hydro is on
         names(h_U_ejecta) = "U_ejecta" ! internal energy of layers above vesc
         names(h_T_ejecta) = "T_ejecta" ! kinetic energy of layers above vesc
         names(h_Omega_ejecta) = "Omega_ejecta" ! gravitational energy of layers above vesc
         names(h_gamma_integral) = "gamma_integral" ! integral of gamma1-4/3 in bound layers
         names(h_max_velocity) = "max_velocity"
         names(h_min_velocity) = "min_velocity"
         names(h_yr_since_coll) = "yr_since_coll"
         names(h_total_J) = "total_J"
         names(h_total_J_bel_vesc) = "total_J_bel_vesc"

         vals(h_gamma_integral) = s% xtra(x_gamma_int_bound)

         if (s% u_flag) then
            vals(h_max_velocity) = maxval(s% u(1:s% nz))
         else
            vals(h_max_velocity) = 0
         end if

         if (s% u_flag) then
            vals(h_min_velocity) = minval(s% u(1:s% nz))
         else
            vals(h_min_velocity) = 0
         end if

         vals(h_yr_since_coll) = s% xtra(x_time_since_first_gamma_zero)/secyer
         
         vals(h_log_R_098) = -100
         vals(h_log_R_095) = -100
         vals(h_log_R_090) = -100
         vals(h_helium_core_mass) = -100
         vals(h_co_core_mass) = 0
         vals(h_log_R_100) = log10(s% r(1)/Rsun)
         vals(h_log_R_vesc_098) = -100
         vals(h_log_R_vesc_095) = -100
         vals(h_log_R_vesc_090) = -100
         vals(h_U_ejecta) = 0d0
         vals(h_T_ejecta) = 0d0
         vals(h_Omega_ejecta) = 0d0
         do k=1, s% nz
            if (s% q(k) < 0.98d0 .and. vals(h_log_R_098) == -100) then
               vals(h_log_R_098) = log10(s% r(k)/Rsun)
            else if (s% q(k) < 0.95d0 .and. vals(h_log_R_095) == -100) then
               vals(h_log_R_095) = log10(s% r(k)/Rsun)
            else if (s% q(k) < 0.9d0 .and. vals(h_log_R_090) == -100) then
               vals(h_log_R_090) = log10(s% r(k)/Rsun)
            end if
            if (vals(h_helium_core_mass) < 0d0 .and. s% X(k) < 0.01d0) then
               vals(h_helium_core_mass) = s% m(k)/Msun
            end if
            if (vals(h_co_core_mass) <= 0d0 .and. s% Y(k) < 0.01d0) then
               vals(h_co_core_mass) = s% m(k)/Msun
            end if
         end do

         if (s% u_flag) then
            do k = s% nz, 1, -1
               v_esc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
               if (s% u(k) > v_esc) then
                  exit
               end if
            end do
            if (k==0) then
               k=1
            end if
         else
            k=1
         end if

         if (.not. s% rotation_flag) then
            vals(h_total_J) = 0d0
            vals(h_total_J_bel_vesc) = 0d0
         else
            vals(h_total_J) = dot_product(s% dm_bar(1:s% nz), s% j_rot(1:s% nz))
            vals(h_total_J_bel_vesc) = dot_product(s% dm_bar(k:s% nz), s% j_rot(k:s% nz))
         end if

         vals(h_log_R_vesc) = log10(s% r(k)/Rsun)
         vals(h_M_below_vesc) = s% m(k)/Msun
         ! get internal radii
         do k0=k, s% nz
            if (s% q(k0)/s% q(k) < 0.98d0 .and. vals(h_log_R_vesc_098) == -100) then
               vals(h_log_R_vesc_098) = log10(s% r(k0)/Rsun)
            else if (s% q(k0)/s% q(k) < 0.95d0 .and. vals(h_log_R_vesc_095) == -100) then
               vals(h_log_R_vesc_095) = log10(s% r(k0)/Rsun)
            else if (s% q(k0)/s% q(k) < 0.9d0 .and. vals(h_log_R_vesc_090) == -100) then
               vals(h_log_R_vesc_090) = log10(s% r(k0)/Rsun)
            end if
         end do
         ! get energies
         if (k>1) then
            do k0 = 1, k
               vals(h_U_ejecta) = vals(h_U_ejecta) + s% dm(k0)*s% energy(k0)
               if (s% u_flag) then
                  vals(h_T_ejecta) = vals(h_T_ejecta) + 0.5d0*s% dm(k0)*s% u(k0)*s% u(k0)
               end if
               vals(h_Omega_ejecta) = vals(h_Omega_ejecta) - s% dm_bar(k0)*s% cgrav(k0)*s% m(k0)/s% r(k0)
            end do 
         end if
         
         if (s% u_flag) then
            vals(h_u_flag) = 1d0
         else
            vals(h_u_flag) = 0d0
         end if

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% u_flag) then
            how_many_extra_profile_columns = 8
         else
            how_many_extra_profile_columns = 1
         end if
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n), alpha, J_inside
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         J_inside = 0d0
         if (s% u_flag) then
            names(1) = "vesc"
            names(2) = "v_div_vesc"
            names(3) = "specific_grav_e"
            names(4) = "specific_kin_e"
            names(5) = "specific_thermal_e"
            names(6) = "total_specific_e"
            names(7) = "mlt_vc"
            names(8) = "spin_parameter" 
            
            do k = s% nz, 1, -1
               vals(k,1) = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
               vals(k,2) = s% u(k)/vals(k,1)
               vals(k,3) = -s% cgrav(k)*s% m(k)/s% r(k)
               vals(k,4) = 0.5d0*s% u(k)*s% u(k)
               vals(k,5) = two_thirds*avo*kerg*s% T(k)/(2*s% mu(k)*s% rho(k)) &
                           + crad*pow4(s% T(k))/s% rho(k)
               vals(k,6) = vals(k,3) + vals(k,4) + vals(k,5)
               vals(k,7) = s% mlt_vc(k)
               if (s% rotation_flag) then
                  J_inside = J_inside + s% j_rot(k)*s% dm_bar(k)
                  vals(k,8) = J_inside*clight/(pow2(s% m(k))*standard_cgrav)
               else
                  vals(k,8) = 0d0
               end if
            end do
         else
            names(1) = "spin_parameter" 

            if (s% rotation_flag) then
               do k = s% nz, 1, -1
                  J_inside = J_inside + s% j_rot(k)*s% dm_bar(k)
                  vals(k,1) = J_inside*clight/(pow2(s% m(k))*standard_cgrav)
               end do
            else
               vals(:,1) = 0d0
            end if
         end if
         
         
      end subroutine data_for_extra_profile_columns
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         integer :: k, k0, k1, num_pts, species, model_number, num_trace_history_values
         real(dp) :: v_esc, time, gamma1_integral, integral_norm, tdyn, &
            max_center_cell_dq, avg_v_div_vesc, energy_removed_layers, dt_next, dt, &
            max_years_for_timestep, omega_crit, &
            denergy
         real(dp) :: core_mass, rmax, alfa, log10_r, lburn_div_lsurf
         logical :: just_did_relax
         character (len=200) :: fname
         include 'formats'
         extras_start_step = terminate
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !this is used to ensure we read the right inlist options
         s% use_other_before_struct_burn_mix = .true.

         ! be sure power info is stored
         call star_set_power_info(s)

         if (.not. in_inlist_pulses .and. .not. s% lxtra(lx_he_zams)) then
            lburn_div_lsurf  = abs(s% L_nuc_burn_total*Lsun/s% L(1))
            if (lburn_div_lsurf > 0d0) then
               if((abs(log10(lburn_div_lsurf))) < 0.01 .and. &
                  (s% star_age > 1d3 .or. s% center_he4 < 0.98d0)) then
                  s% use_other_before_struct_burn_mix = .false.
                  call star_relax_uniform_omega(id, 1, s% job% new_omega_div_omega_crit,&
                                                s% job% num_steps_to_relax_rotation, 1d0, ierr)
                  s% use_other_before_struct_burn_mix = .true.
                  if (ierr/=0) then
                     write(*,*) "error when relaxing omega at he ZAMS"
                     stop
                  end if
               s% lxtra(lx_he_zams) = .true.
               end if
            end if
         end if

         ! Ignore energy checks before first time hydro is turned on
         ! otherwise need small steps during core helium burning and it
         ! slows down the test_suite
         if (.not. s% lxtra(lx_hydro_has_been_on)) then
            s% cumulative_energy_error = 0d0
            s% cumulative_extra_heating = 0d0
            if(mod(s%model_number, s%terminal_interval) == 0) then
               write(*,*) &
                  "Setting energy conservation error to zero until hydro is turned on for the first time"
            end if
         end if

         if (pgstar_flag) then
            s% pgstar_interval = pgstar_interval
            s% Grid2_file_interval = pgstar_file_interval
         end if
         just_did_relax = .false.
         if (s% u_flag) then ! get point where v<vesc
            do k = s% nz-1, 1, -1
               v_esc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
               if (s% u(k) > v_esc) then
                  exit
               end if
            end do
            if (k>0 .and. k < s% nz) k = k+1
         else
            k = 0
         end if

         ! s% xtra(x_gamma_int_bound) stores the gamma integral
         ! during pulses integrate gamma only inside bound regions
         ! otherwise compute it for the whole star
         s% xtra(x_gamma_int_bound) = -1d0

         ! can be adjusted below if nearing breakout
         s% profile_interval = 100

         if (s% u_flag .and. k > 0 .and. s% xtra(x_time_start_pulse) > 0d0) then

            ! check energy and average escape velocity of outer layers
            avg_v_div_vesc = 0d0
            energy_removed_layers = 0d0
            do k0 = 1, k
               avg_v_div_vesc = avg_v_div_vesc + s% dm(k0)*s% u(k0)/sqrt(2*s% cgrav(k0)*s% m(k0)/(s% r(k0)))
               energy_removed_layers = energy_removed_layers + &
                  0.5d0*s% dm(k0)*s% u(k0)*s% u(k0) - s% dm(k0)*s% cgrav(k0)*s% m(k0)/s% r(k0) &
                  +s% energy(k0)*s% dm(k0)
            end do
            ! adjust location of boundary to remove by considering also 
            ! material below the escape velocity that has a positive net
            ! total specific energy. 
            if (energy_removed_layers > 0d0) then ! possible to eject material
               if(mod(s%model_number, s%terminal_interval) == 0) then
                  write(*,*) "k, q, energy_removed_layers before adjustment is", k, s% q(k), energy_removed_layers
               end if
               do k0 = k+1, s% nz
                  denergy = &
                     0.5d0*s% u(k0)*s% u(k0) - s% cgrav(k0)*s% m(k0)/s% r(k0) &
                        +s% energy(k0)
                  if (denergy < 0d0) then
                     k = k0-1
                     exit
                  else
                     energy_removed_layers = energy_removed_layers+denergy*s% dm(k0)
                     avg_v_div_vesc = avg_v_div_vesc + &
                        s% dm(k0)*s% u(k0)/sqrt(2*s% cgrav(k0)*s% m(k0)/(s% r(k0)))
                  end if
               end do 
               if(mod(s%model_number, s%terminal_interval) == 0) then
                  write(*,*) "k, q, energy_removed_layers after adjustment is", k, s% q(k), energy_removed_layers
               end if
            end if

            avg_v_div_vesc = avg_v_div_vesc/(s% m(1) - s% m(k))

            ! compute gamma integral in what will remain of the star
            integral_norm = 0d0
            gamma1_integral = 0d0
            do k0=k,s% nz
               integral_norm = integral_norm + (s% Pgas(k0)+s% Prad(k0))*s% dm(k0)/s% rho(k0)
               gamma1_integral = gamma1_integral + &
                  (s% gamma1(k0)-4.d0/3.d0)*(s% Pgas(k0)+s% Prad(k0))*s% dm(k0)/s% rho(k0)
            end do
            gamma1_integral = gamma1_integral/max(1d-99,integral_norm)
            s% xtra(x_gamma_int_bound) = gamma1_integral

            do k1 = 1, s% nz
               if (s% q(k1) < 0.9d0) then
                  !increase profile density near breakout, check that at q=0.9 velocities are
                  ! above 500 km/s, and that at the surface they're below that.
                  ! for the first pulse, increase profile density from the onset of instability
                  ! to breakout
                  if ((s% u(k1)>5d7 .and. s% u(1)<5d7) &
                     .or. (s% ixtra(ix_num_relaxations) == 0 .and. gamma1_integral < 0d0 .and. s% u(1)<5d7)) then
                     s% profile_interval = 10
                  else
                     s% profile_interval = 100
                  end if
                  exit
               end if
            end do

            ! To relax the star after a pulse, we check for a series of conditions that
            ! must apply in layers below q=q_for_relax_check
            do k0 = k, s% nz
               if (s% q(k0) < q_for_relax_check*s% q(k)) then
                  exit
               end if
            end do
 
            if(mod(s%model_number, s%terminal_interval) == 0) then
               write(*,*) 'Layers above q=', s% q(k), 'will be removed'
               write(*,*) 'checking for conditions inside q=', q_for_relax_check, 'of material that will remain'
               write(*,*) 'check time left', &
                  s% xtra(x_time_start_pulse) + s% xtra(x_dyn_time)*num_dyn_ts_for_relax - s% time
               write(*,*) 'max vel inside fraction of cutoff', maxval(abs(s% u(k0:s% nz)))/1e5
               write(*,*) 'max c/cs inside fraction of cutoff', maxval(abs(s% u(k0:s% nz)/s% csound(k0:s% nz)))
               write(*,*) 'average v/vesc outside cutoff', avg_v_div_vesc
               write(*,*) 'Kinetic plus potential energy outside cutoff', energy_removed_layers
               write(*,*) 'mass inside cutoff', k, s% m(k)/Msun
               write(*,*) 'relax counter', s% ixtra(ix_steps_met_relax_cond)
               write(*,*) 'log max v [km/s]=', safe_log10(maxval(s% u(1:s% nz))/1e5)
            end if
            ! If enough time has passed after a pulse, then turn off Riemann hydro
            ! also, wait for the inner q_for_relax_check of what would remain of the star to be below a threshold
            ! in v/cs, and below a threshold in v
            ! Verify also that the conditions to turn on Riemann hydro are not still satisfied
            ! For details on all these options check inlist_project

            ! ignore if s% q(k) < 1d-3, in that case it's very likely a PISN
            if (s% q(k) > 1d-3) then
               if (s% lxtra(lx_hydro_on) .and. s% xtra(x_time_start_pulse) > 0d0 &
                  .and. s% time > s% xtra(x_time_start_pulse) + s% xtra(x_dyn_time)*num_dyn_ts_for_relax &
                  .and. maxval(abs(s% u(k0:s% nz)))/1d5 < max_v_for_relax &
                  .and. maxval(abs(s% u(k0:s% nz)/s% csound(k0:s% nz))) < max_machn_for_relax &
                  .and. gamma1_integral > min_gamma_sub_43_for_hydro &
                  .and. log10(s% power_neutrinos) < max_Lneu_for_relax &
                  .and. log10(s% power_nuc_burn) < max_Lnuc_for_relax &
                  .and. energy_removed_layers > 0d0) then
                  s% ixtra(ix_steps_met_relax_cond) = s% ixtra(ix_steps_met_relax_cond) + 1
               else
                  s% ixtra(ix_steps_met_relax_cond) = 0
               end if
            else
               ! escape velocity reached within a tiny fraction of the
               ! core. Before marking as PISN verify if any cell above 
               ! this is below the escape velocity
               do k0 = k, 1, -1
                  v_esc = sqrt(2*s% cgrav(k0)*s% m(k0)/(s% r(k0)))
                  if (s% u(k0) < v_esc) then
                     write(*,*) "Likely PISN", s% q(k0), s% u(k0)/v_esc
                     exit
                  end if
               end do
               ! If all cells were above the escape velocity, we
               ! mark this as a PISN
               if (k0 == 0) then
                  extras_start_step = terminate
                  write(fname, fmt="(a19)") 'LOGS/pisn_prof.data'
                  call star_write_profile_info(id, fname, ierr)
                  write(*,*) "Entire star is expanding above the escape velocity, PISN!"
                  return
               end if
            end if
            if(s% ixtra(ix_steps_met_relax_cond) >= num_steps_before_relax .and. s% q(k) > 1d-3) then
               write(*,*) "Relaxing model to lower mass!"
               s% ixtra(ix_num_relaxations) = s% ixtra(ix_num_relaxations) + 1
               s% xtra(x_star_age_at_relax) = s% star_age
               
               !save a profile just before relaxation
               write(fname, fmt="(a18,i0.3,a5)") 'LOGS/prerelax_prof', s% ixtra(ix_num_relaxations), '.data'
               call star_write_profile_info(id, fname, ierr)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return

               s% lxtra(lx_hydro_on) = .false.
               s% lxtra(lx_have_saved_gamma_prof) = .false.

               s% ixtra(ix_steps_since_relax) = 0

               write(*,*) "removing cells", k, s% m(k)/Msun

               max_center_cell_dq = s% max_center_cell_dq
               s% max_center_cell_dq = s% dq(s% nz)
               dt = s% dt
               dt_next = s% dt_next

               call star_read_controls(id, 'inlist_hydro_off', ierr)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return

               call star_set_u_flag(id, .false., ierr)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return

               call my_before_struct_burn_mix(s% id, s% dt, extras_start_step)

               ! to avoid showing pgstar stuff during relax
               s% pgstar_interval = 100000000
               s% Grid2_file_interval = 100000000
               s% job% pgstar_flag = .false.

               max_years_for_timestep = s% max_years_for_timestep
               s% max_years_for_timestep = 1d0

               s% delta_lgL_nuc_limit = -1d0
               s% delta_lgL_nuc_hard_limit = -1d0
               s% use_other_before_struct_burn_mix = .false.
               s% timestep_hold = 0

               call star_relax_to_star_cut(s% id, k, .true., .true., .true., ierr)
               if (ierr /= 0) then
                  write(*,*) "error when removing mass through star_relax_to_star_cut", ierr
                  stop
               end if

               s% max_years_for_timestep = max_years_for_timestep
               s% photo_interval = 100
               s% dt_next = min(1d2, dt_next)
               s% dt = min(1d2, dt)

               s% ixtra(ix_steps_met_relax_cond) = 0

               s% max_center_cell_dq = max_center_cell_dq
               if (pgstar_flag) then
                  s% pgstar_interval = pgstar_interval
                  s% Grid2_file_interval = pgstar_file_interval
                  s% job% pgstar_flag = .true.
               end if
               just_did_relax = .true.

               !save a profile and a model just after relaxation
               write(fname, fmt="(a19,i0.3,a5)") 'LOGS/postrelax_prof', s% ixtra(ix_num_relaxations), '.data'
               call star_write_profile_info(id, fname, ierr)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return
               write(fname, fmt="(a20,i0.3,a4)") 'LOGS/postrelax_model', s% ixtra(ix_num_relaxations), '.mod'
               call star_write_model(id, fname, ierr)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return

               s% use_other_before_struct_burn_mix = .true.
            end if
         end if

         ! turn on Riemman hydro if close to instability
         integral_norm = 0d0
         gamma1_integral = 0d0
         do k=1,s% nz
            integral_norm = integral_norm + (s% Pgas(k)+s% Prad(k))*s% dm(k)/s% rho(k)
            gamma1_integral = gamma1_integral + &
               (s% gamma1(k)-4.d0/3.d0)*(s% Pgas(k)+s% Prad(k))*s% dm(k)/s% rho(k)
         end do
         gamma1_integral = gamma1_integral/max(1d-99,integral_norm)
         if (s% xtra(x_gamma_int_bound) == -1d0) then
            ! if xtra(x_gamma_int_bound) is different from -1d0 it means it was computed for the bound material,
            ! use that instead
            s% xtra(x_gamma_int_bound) = gamma1_integral
         end if
         ! Save profile if gamma1_integral becomes negative
         if (gamma1_integral < 0d0 .and. .not. s% lxtra(lx_have_saved_gamma_prof)) then
            s% lxtra(lx_have_saved_gamma_prof) = .true.
            write(fname, fmt="(a16,i0.3,a5)") 'LOGS/gamma1_prof', s% ixtra(ix_num_relaxations)+1, '.data'
            call star_write_profile_info(id, fname, ierr)
            if (dbg) write(*,*) "check ierr", ierr
            if (ierr /= 0) return
         end if
         if(mod(s%model_number, s%terminal_interval) == 0) then
             write(*,*) "check gamma integral", gamma1_integral
         end if

         if (just_did_relax) then
            extras_start_step = keep_going
            return
         end if

         if (s% ixtra(ix_steps_since_relax) > 10 .and. .not. just_did_relax .and. .not. s% lxtra(lx_hydro_on) &
            .and. gamma1_integral < min_gamma_sub_43_for_hydro) then
            write(*,*) "Turning on Riemann hydro!"

            call star_read_controls(id, 'inlist_hydro_on', ierr)
            if (dbg) write(*,*) "check ierr", ierr
            if (ierr /= 0) return
            call star_read_controls(id, 'inlist_after_first_pulse', ierr)
            if (dbg) write(*,*) "check ierr", ierr
            if (ierr /= 0) return
            s% dt_next = min(1d2,s% dt_next)
            s% dt = min(1d2,s% dt)

            !ensure implicit wind is not used
            s% was_in_implicit_wind_limit = .false.

            write(fname, fmt="(a18,i0.3,a5)") 'LOGS/prehydro_prof', s% ixtra(ix_num_relaxations)+1, '.data'
            call star_write_profile_info(id, fname, ierr)
            if (dbg) write(*,*) "check ierr", ierr
            if (ierr /= 0) return
            write(fname, fmt="(a19,i0.3,a4)") 'LOGS/prehydro_model', s% ixtra(ix_num_relaxations)+1, '.mod'
            call star_write_model(id, fname, ierr)
            if (dbg) write(*,*) "check ierr", ierr
            if (ierr /= 0) return

            call star_set_v_flag(id, .false., ierr)!TODO store v
            if (dbg) write(*,*) "check ierr", ierr
            if (ierr /= 0) return

            !s% job% pgstar_flag = .false.
            call star_set_u_flag(id, .true., ierr)
            if (dbg) write(*,*) "check ierr", ierr
            if (ierr /= 0) return
            s% lxtra(lx_hydro_on) = .true.
            s% lxtra(lx_hydro_has_been_on) = .true.
            s% xtra(x_time_start_pulse) = -1
            s% ixtra(ix_steps_since_hydro_on) = 0

            !force initial velocities to zero to prevent issues in outer layers
            s% xh(s% i_u,:) = 0d0
            s% xh_old(s% i_u,:) = 0d0
            s% generations = 1
         end if

         ! check time to relax model
         ! the star is considered to be in a pulse if the velocity is high enough
         ! only check velocity in q_for_relax_check of the mass of the star, this avoids
         ! considering as a pulse cases where the surface goes a bit nuts
         do k0 = 1, s% nz
            if (s% q(k0) < q_for_relax_check) then
               exit
            end if
         end do
         if (s% lxtra(lx_hydro_on) .and. s% xtra(x_time_start_pulse) < 0d0) then
            if ((maxval(abs(s% xh(s% i_u, k0:s% nz)))/1d5 > max_v_for_pulse .or. gamma1_integral < 1d-3)) then
               core_mass = s% m(1)
               !if(s% o_core_mass > 0d0) then
               !   core_mass = s% o_core_mass*Msun
               !else if(s% c_core_mass > 0d0) then
               if(s% co_core_mass > 0d0) then
                  core_mass = s% co_core_mass*Msun
               else if (s% he_core_mass > 0d0) then
                  core_mass = s% he_core_mass*Msun
               end if
               do k=1, s% nz
                  if (s% m(k)/core_mass < q_for_dyn_ts) then
                     exit
                  end if
               end do
               tdyn = 1d0/sqrt(s% m(k)/(4d0/3d0*pi*pow3(s% r(k)))*standard_cgrav)
               s% xtra(x_time_start_pulse) = s% time
               s% xtra(x_dyn_time) = tdyn
               write(*,*) "reached high velocities", maxval(abs(s% xh(s% i_u, k0:s% nz)))/1e5, &
                  "stopping in at least", num_dyn_ts_for_relax*tdyn, "seconds"
            end if
         end if

         !After a relax, wait for ten days before turning on v_flag
         !this avoids the surface post-relax from going crazy
         if((logT_for_v_flag < log10(s% T(s% nz)) .or. logLneu_for_v_flag < safe_log10(s% power_neutrinos)) &
               .and. .not. s% u_flag .and. .not. s% v_flag) then
               write(*,*) "log10 central T has lowered below logT_for_v_flag, turn on v_flag"
               call star_set_v_flag(id, .true., ierr)
               s% dt_next = min(1d2, s% dt_next)
               s% dt = min(1d2, s% dt)
               if (dbg) write(*,*) "check ierr", ierr
               if (ierr /= 0) return
         end if

         !Always call this at the end to ensure we are using the correct
         !inlists
         call my_before_struct_burn_mix(s% id, s% dt, extras_start_step)

         extras_start_step = keep_going
      end function extras_start_step
   
      subroutine my_before_struct_burn_mix(id, dt, res)
         use const_def, only: dp
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: dt
         integer, intent(out) :: res ! keep_going, redo, retry, terminate
         real(dp) :: power_photo, v_esc
         integer :: ierr, k
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !read proper options when hydro is on/off
         !do this to ensure proper behaviour of retries
         if(s% u_flag) then
            call star_read_controls(id, 'inlist_hydro_on', ierr)
            if (s% xtra(x_time_start_pulse) > 0d0) then
               s% max_timestep = 1d99
               do k = s% nz, 1, -1
                  v_esc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
                  if (s% u(k) > 2*v_esc) then
                     exit
                  end if
               end do
            else
               s% max_timestep = max_dt_before_pulse
            end if
         else
            s% max_timestep = 1d99
            call star_read_controls(id, 'inlist_hydro_off', ierr)
         end if

         !ignore L_nuc limit if L_phot is too high or if we just did a relax
         !(ixtra(ix_steps_since_relax) is set to zero right after a relax)
         
         !when L_phot exceeds max_Lphoto_for_lgLnuc_limit, the timestep limit is applied
         !to L_phot instead

         ! be sure power info is stored
         call star_set_power_info(s)

         power_photo = dot_product(s% dm(1:s% nz), s% eps_nuc_categories(iphoto,1:s% nz))/Lsun
         if (safe_log10(abs(power_photo)) > max_Lphoto_for_lgLnuc_limit2) then
            s% delta_lgL_nuc_limit = -1d0
            s% delta_lgL_nuc_hard_limit = -1d0
            s% delta_lgL_power_photo_limit = -1d0
            s% delta_lgL_power_photo_hard_limit = -1d0
         else
            if (s% ixtra(ix_steps_since_relax) == 0 &
                  .or. safe_log10(abs(power_photo)) > max_Lphoto_for_lgLnuc_limit) then
               s% delta_lgL_nuc_limit = -1d0
               s% delta_lgL_nuc_hard_limit = -1d0
            else
               s% delta_lgL_nuc_limit = delta_lgLnuc_limit
               s% delta_lgL_nuc_hard_limit = 2d0*delta_lgLnuc_limit
            end if
            if (safe_log10(abs(power_photo)) > max_Lphoto_for_lgLnuc_limit) then
               s% delta_lgL_power_photo_limit = delta_lgLnuc_limit
               s% delta_lgL_power_photo_hard_limit = 2d0*delta_lgLnuc_limit
            else
               s% delta_lgL_power_photo_limit = -1d0
               s% delta_lgL_power_photo_hard_limit = -1d0
            end if
         end if

         !ignore winds if neutrino luminosity is too high or for a few steps after
         !a relax
         if(s% ixtra(ix_steps_since_relax) < 50 &
               .or. safe_log10(s% power_neutrinos) > max_Lneu_for_mass_loss &
               .or. s% u_flag) then
            s% use_other_wind = .false.
            s% was_in_implicit_wind_limit = .false.
         else
            s% use_other_wind = .true.
         end if

         if (maxval(s% T(1:s% nz)) > 9.8) then
            s% delta_lgRho_cntr_hard_limit = -1d0
         else
            s% delta_lgRho_cntr_hard_limit = delta_lgRho_cntr_hard_limit
         end if

         ! reading inlists can turn this flag off for some reason
         s% use_other_before_struct_burn_mix = .true.

         res = keep_going
      end subroutine my_before_struct_burn_mix
      
      subroutine null_binary_controls(id, binary_id, ierr)
         integer, intent(in) :: id, binary_id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_binary_controls

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         use run_star_support
         integer, intent(in) :: id
         integer :: ierr,k
         real(dp) :: max_vel_inside
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_finish_step = keep_going

         !count time since first collapse
         if (s% lxtra(lx_have_reached_gamma_limit)) then
            s% xtra(x_time_since_first_gamma_zero) = &
               s% xtra(x_time_since_first_gamma_zero) + s% dt 
         end if

         s% ixtra(ix_steps_since_relax) = s% ixtra(ix_steps_since_relax) + 1
         s% ixtra(ix_steps_since_hydro_on) = s% ixtra(ix_steps_since_hydro_on) + 1

         if (s% ixtra(ix_num_relaxations) == 1 .and. stop_100d_after_pulse &
               .and. s% star_age - s% xtra(x_star_age_at_relax) > 100d0/dayyer) then
            !for the test_suite, terminate at the onset of the second pulse
            extras_finish_step = terminate
            s% termination_code = t_xtra1
            termination_code_str(t_xtra1) = "Successful test: evolved 100 days past first relax"
            return
         end if

         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

      end function extras_finish_step
      
      end module run_star_extras
      
