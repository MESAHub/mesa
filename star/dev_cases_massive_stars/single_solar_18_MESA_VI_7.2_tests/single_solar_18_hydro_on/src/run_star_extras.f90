! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
      use chem_def
      
      implicit none

      ! s% xtra(x_old_wind) every step we save the used wind mass loss rate here,
      ! this is used to soften too large changes in wind mass loss rates 
      integer, parameter :: x_old_wind = 1

      ! s% xtra(x_time_thermal_eq) holds time the star has been in thermal equilibrium
      ! used to determine the start of ZAMS
      integer, parameter :: x_time_thermal_eq = 2

      ! s% lxtra(lx_pre_ZAMS) true if the star has not already properly settled in the ZAMS
      integer, parameter :: lx_pre_ZAMS = 1

      ! s% xctrl(x_cool_scaling_factor) holds time the star has been in thermal equilibrium
      integer, parameter :: x_cool_scaling_factor = 10
      ! s% xctrl(x_hot_scaling_factor) holds time the star has been in thermal equilibrium
      integer, parameter :: x_hot_scaling_factor = 11
      ! s% xctrl(x_total_scaling_factor) holds time the star has been in thermal equilibrium
      integer, parameter :: x_total_scaling_factor = 12
      
      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         s% other_wind => brott_wind_modified

      end subroutine extras_controls
      
      subroutine brott_wind_modified(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
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
            vink_wind, nieu_wind, hamann_wind, highT_w, lowT_w, Twindow
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
         ! For wind scaling we use the ratio of iron abundance to the A09 value
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

         lowT_w = max(vink_wind, nieu_wind)

         alfa = 0d0
         if (Xs > 0.7d0) then
            alfa = 1d0
         else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
            alfa = (Xs - 0.4d0)/0.3d0
         end if
         highT_w = alfa * vink_wind + (1d0-alfa) * hamann_wind

         ! have a 10% Teff_jump window to switch from the lowT to the highT wind
         Twindow = Teff_jump*0.10d0
         alfa = 0d0
         if (T1 < Teff_jump - Twindow/2d0) then
            alfa = 1d0
         else if (T1 > Teff_jump - Twindow/2d0 .and. T1 < Teff_jump + Twindow/2d0) then
            alfa = ((Teff_jump + Twindow/2d0)-T1)/Twindow
         end if

         lowT_w = s% x_ctrl(x_cool_scaling_factor)*lowT_w
         highT_w = s% x_ctrl(x_hot_scaling_factor)*highT_w

         w = alfa * lowT_w + (1d0-alfa) * highT_w
         
         if (s% center_h1 < 1e-3) then
            w = s% x_ctrl(x_total_scaling_factor) * w
         end if 

         ! soften change in wind to avoid things going bad
         if (s% xtra(x_old_wind) /= 0) then
            if(abs(w) > abs(s% xtra(x_old_wind))*1.05) then
               w = s% xtra(x_old_wind)*1.05
            else if(abs(w) < abs(s% xtra(x_old_wind))*0.95) then
               w = s% xtra(x_old_wind)*0.95
            end if
         end if
         s% xtra(x_old_wind) = w

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
               w1 = 10**(logMdot)
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
               w2 = 10**(logMdot)
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
                     +0.81d0*log10(R1/Rsun) !&
             !+0.85d0*log10(Z_div_Z_solar) ! we do not apply the Vink Z scaling here
            w = 10**(log10w)
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
            w = 10**(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind_modified
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (.not. restart) then
            s% xtra(x_old_wind) = 0d0
            s% xtra(x_time_thermal_eq) = 0d0
            s% lxtra(lx_pre_ZAMS) = .true.
         end if
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: omega_crit
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0

         ! check if the star is not yet at ZAMS. If that's the case, keep rotation fixed
         ! s% xtra(x_time_thermal_eq) is set in extras_check_model
         if (s% lxtra(lx_pre_ZAMS)) then
            write(*,*) "Check thermal timescale", standard_cgrav*(s% m(1))**2/(s% r(1)*s% L(1))/secyer
            if (s% xtra(x_time_thermal_eq) > standard_cgrav*(s% m(1))**2/(s% r(1)*s% L(1))) then
               write(*,*) "thermal equilibrium for GM^2/RL, normal evolution from here"
               s% lxtra(lx_pre_ZAMS) = .false.
            else if (s% center_h1 < 0.69d0) then
               write(*,*) "Have burned significant hydrogen, normal evolution from here"
               s% lxtra(lx_pre_ZAMS) = .false.
            else if (s% star_age*secyer > 5*standard_cgrav*(s% m(1))**2/(s% r(1)*s% L(1))) then
               write(*,*) "WARNING: no equilibrium found after evolving for 5GM^2/RL"
               write(*,*) "Switching to normal evolution"
               s% lxtra(lx_pre_ZAMS) = .false.
            else
               ! keep rotation fixed
               write(*,*) "Not at ZAMS yet, keeping omega_div_omega_crit fixed"
               omega_crit = star_surface_omega_crit(id,ierr)
               if (ierr /= 0) then
                  write(*,*) "Error in star_surface_omega_crit"
                  return
               end if
               call star_set_uniform_omega(s% id, omega_crit*s% job% new_omega_div_omega_crit, ierr)
               if (ierr /= 0) then
                  write(*,*) "Error in star_set_uniform_omega"
                  return
               end if
            end if
         end if

         !!slowly turn on superad reduction to make it easier to produce pre-MS models
         !if (s% star_age >= 0.01d0) then
         !   s% superad_reduction_diff_grads_limit = 1d-2
         !else
         !   s% superad_reduction_diff_grads_limit = 10**(2-4d0*(s% star_age)/0.01d0)
         !end if

         !write(*,*) "check superad_reduction_diff_grads_limit", s% superad_reduction_diff_grads_limit

         if (s% center_he4 < 1d-3 .and. s% center_c12 < 1d-3) then
             ! use stricter timestep control on evolution of center density
             s% delta_lgRho_cntr_limit = 0.005d0 
             s% delta_lgRho_cntr_hard_limit = 0.01d0 
             s% convergence_ignore_equL_residuals = .true.
         else
             s% delta_lgRho_cntr_limit = 0.02d0 
             s% delta_lgRho_cntr_hard_limit = 0.04d0
             if (s% lxtra(lx_pre_ZAMS)) then
                s% convergence_ignore_equL_residuals = .true.
             else
                s% convergence_ignore_equL_residuals = .false.
             end if
         end if

      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 2
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
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         
         names(1) = 'pre_ZAMS'
         names(2) = 'max_superad_reduction'
         if (s% lxtra(lx_pre_ZAMS)) then
            vals(1) = 1d0
         else
            vals(1) = 0d0
         end if
         vals(2) = maxval(s% superad_reduction_factor(1:s% nz))
         write(*,*) "Check max superad reduction", vals(2)

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 3
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, j
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.
         names(1) = "Lrad_div_Leddd"
         names(2) = "Lrad_div_Leddd2"
         names(3) = "conv_vel_km_s"
         do k=1,s% nz
         vals(k,1) = 4d0*crad/3d0*pow4(s% T(k))/s% Peos(k)*s%gradT(k)
         vals(k,2) = 4d0*crad/3d0*pow4(s% T(k))/s% Peos(k)*s%mlt_gradT(k)
         vals(k,3) = s% conv_vel(k) / 1e5 ! in km/s 
         end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! check for Lnuc ~ L, this is used in star step to determine if we
         ! keep rotation fixed
         if (abs(log10(abs(s% L_nuc_burn_total * Lsun / s% L(1)))) < 0.005) then
            s% xtra(x_time_thermal_eq) = s% xtra(x_time_thermal_eq) + s% dt_old
         else
            s% xtra(x_time_thermal_eq) = 0d0
         end if

         write(*,*) 'Central He abundance:', s% center_he4
         write(*,*) 'Central He abundance:', s% center_c12
      
         ! if (s% center_he4 < 1d-3 .and. s% center_c12 < 1d-2) then
         !    s% need_to_save_profiles_now = .true.
         !    extras_finish_step = terminate
         !    write(*,*) "Terminate due to carbon depletion"
         ! end if
         if (s% center_he4 < 1d-2 .and. s% center_c12 < 1d-2 .and. s% center_o16 < 1d-2) then
             extras_finish_step = terminate
             write(*,*) "Terminate due to oxygen depletion"
         end if

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      end module run_star_extras
      
