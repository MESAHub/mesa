    ! ***********************************************************************
    !
    !   Copyright (C) 2010  Bill Paxton
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

    ! these routines are called by the standard run_star check_model
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
        
        if (s% job% create_pre_main_sequence_model) return
        
        s% other_adjust_mlt_gradT_fraction => other_adjust_mlt_gradT_fraction_Peclet
        s% other_overshooting_scheme => extended_convective_penetration
        
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
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step

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
      

      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model
      
      
    subroutine extras_after_evolve(id, ierr)
        use num_lib
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        real(dp) :: dt
        integer :: k, k1, k2, k3, k4
        real(dp) :: Peclet_number, fraction, gradient
        type (star_info), pointer :: s
        logical :: okay

        include 'formats'

        okay = .true.

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        call test_suite_after_evolve(s, ierr)
        
        if (s% job% create_pre_main_sequence_model) return

        write(*,*)
        k1 = 0
        k2 = 0
        k3 = 0
        k4 = 0
        
        do k = s% nz, 1, -1
           if (s% m(k) > 0.8_dp*Msun .and. k1 == 0) k1 = k
           if (s% m(k) > 0.95_dp*Msun .and. k2 == 0) k2 = k
           if (s% m(k) > 1.1_dp*Msun .and. k3 == 0) k3 = k
           if (s% m(k) > 1.5_dp*Msun .and. k4 == 0) k4 = k
        end do

        call check_int('mixing type at 0.8 Msun', s% mixing_type(k1), convective_mixing, convective_mixing)
        call check_int('mixing type at 0.95 Msun', s% mixing_type(k2), overshoot_mixing, overshoot_mixing)
        call check_int('mixing type at 1.1 Msun', s% mixing_type(k3), overshoot_mixing, overshoot_mixing)
        call check_int('mixing type at 1.5 Msun', s% mixing_type(k4), minimum_mixing, minimum_mixing)

        call calc_Peclet_number(id, k3, Peclet_number, ierr)
        fraction = (safe_log10(Peclet_number)+2.0_dp)/4.0_dp
        gradient = (fraction*s%grada(k3) + (1.0_dp - fraction) * s%gradr(k3))

        write(*,*)
        call check('grada fraction of gradT at 1.1 Msun', fraction, 0.1_dp, 0.9_dp )
        write(*,*)

        call check('gradT at 0.8 Msun', s%gradT(k1), s%grada(k1) - 1.0d-3, s%grada(k1) + 1.0d-3)   ! gradT = grada
        call check('gradT at 0.95 Msun', s%gradT(k2), s%grada(k2) - 1.0d-3, s%grada(k2) + 1.0d-3)   ! gradT = grada
        call check('gradT at 1.1 Msun', s%gradT(k3), gradient - 1.0d-3, gradient + 1.0d-3)   ! gradT = f*grada + (1-f)*gradr
        call check('gradT at 1.5 Msun', s%gradT(k4), s%gradr(k4) - 1.0d-3, s%gradr(k4)+ 1.0d-3)   ! gradT = gradr
        write(*,*)

        call check('Dmix at 0.95 Msun', s%d_mix(k2), s%d_mix(s%conv_bdy_loc(1))*0.95d0, s%d_mix(s%conv_bdy_loc(1))*1.05d0)
        call check('Dmix at 1.1 Msun', s%d_mix(k3), 1.0d4, s%d_mix(s%conv_bdy_loc(1)) * 1.0d-2)

        write(*,*)
        if (okay) write(*,'(a)') 'All values are within tolerances'
        write(*,*)

        contains

        subroutine check(str, val, low, hi)
            real(dp), intent(in) :: val, low, hi
            character (len=*) :: str
            include 'formats'
            if (low <= val .and. val <= hi) then
                write(*,1) trim(str), val, low, hi
            else
                write(*,1) '*** BAD *** ' // trim(str), val, low, hi
                okay = .false.
            end if
        end subroutine check

        subroutine check_int(str, val, low, hi)
            integer, intent(in) :: val, low, hi
            character (len=*) :: str
            include 'formats'
            if (low <= val .and. val <= hi) then
                write(*,11) trim(str), val, low, hi
            else
                write(*,11) '*** BAD *** ' // trim(str), val, low, hi
                okay = .false.
            end if
        end subroutine check_int
      
    end subroutine extras_after_evolve


    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine extended_convective_penetration(id, i, j, k_a, k_b, D, vc, ierr)
        integer, intent(in) :: id, i, j
        integer, intent(out) :: k_a, k_b
        real(dp), intent(out), dimension(:) :: D, vc
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        logical, parameter :: DEBUG = .FALSE.
        real(dp) :: f, f2, f0
        real(dp) :: D0, Delta0
        real(dp) :: w
        real(dp) :: factor
        real(dp) :: r_cb, Hp_cb
        real(dp) :: r_ob, D_ob, vc_ob
        logical  :: outward
        integer  :: dk, k, k_ob
        real(dp) :: r, dr, r_step

        ! Evaluate the overshoot diffusion coefficient D(k_a:k_b) and
        ! mixing velocity vc(k_a:k_b) at the i'th convective boundary,
        ! using the j'th set of overshoot parameters. The overshoot
        ! follows the extended convective penetration scheme description by Mathias
        ! Michielsen, "Probing the shape of the mixing profile and of the thermal
        ! structure at the convective core boundary through asteroseismology",
        ! A&A, 628, 76 (2019)

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        ! Extract parameters
        f = s%overshoot_f(j)
        f0 = s%overshoot_f0(j)
        f2 = s%x_ctrl(j)

        D0 = s%overshoot_D0(j)
        Delta0 = s%overshoot_Delta0(j)

        if (f < 0.0_dp .OR. f0 <= 0.0_dp .OR. f2 < 0.0_dp) then
            write(*,*) 'ERROR: for extended convective penetration, must set f0 > 0, and f and f2 >= 0'
            write(*,*) 'see description of overshooting in star/defaults/control.defaults'
            ierr = -1
            return
        end if

        ! Apply mass limits
        if (s%star_mass < s%overshoot_mass_full_on(j)) then
            if (s%star_mass > s%overshoot_mass_full_off(j)) then
                w = (s%star_mass - s%overshoot_mass_full_off(j)) / &
                  (s%overshoot_mass_full_on(j) - s%overshoot_mass_full_off(j))
                factor = 0.5_dp*(1.0_dp - cospi(w))
                f = f*factor
                f0 = f0*factor
                f2 = f2*factor
            else
                f = 0.0_dp
                f0 = 0.0_dp
                f2 = 0.0_dp
            endif
        endif

        ! Evaluate convective boundary (_cb) parameters
        call star_eval_conv_bdy_r(s, i, r_cb, ierr)
        if (ierr /= 0) return

        call star_eval_conv_bdy_Hp(s, i, Hp_cb, ierr)
        if (ierr /= 0) return

        ! Evaluate overshoot boundary (_ob) parameters
        call star_eval_over_bdy_params(s, i, f0, k_ob, r_ob, D_ob, vc_ob, ierr)
        if (ierr /= 0) return

        ! Loop over cell faces, adding overshoot until D <= overshoot_D_min
        outward = s%top_conv_bdy(i)

        if (outward) then
            k_a = k_ob
            k_b = 1
            dk = -1
        else
            k_a = k_ob+1
            k_b = s%nz
            dk = 1
        endif

        if (f > 0.0_dp) then
            r_step = f*Hp_cb
        else
            r_step = 0.0_dp
        endif

        face_loop : do k = k_a, k_b, dk
            ! Evaluate the extended convective penetration factor
            r = s%r(k)
            if (outward) then
                dr = r - r_ob
            else
                dr = r_ob - r
            endif

            if (dr < r_step .AND. f > 0.0_dp) then  ! step factor
                factor = 1.0_dp
            else
                if ( f2 > 0.0_dp) then                ! exponential factor
                    factor = exp(-2.0_dp*(dr-r_step)/(f2*Hp_cb))
                else
                    factor = 0.0_dp
                endif
            endif

            ! Store the diffusion coefficient and velocity
            D(k) = (D0 + Delta0*D_ob)*factor
            vc(k) = (D0/D_ob + Delta0)*vc_ob*factor

            ! Check for early overshoot completion

            if (D(k) < s%overshoot_D_min) then
                k_b = k
                exit face_loop
            endif

        end do face_loop

        if (DEBUG) then
            write(*,*) 'step exponential overshoot:'
            write(*,*) '  k_a, k_b   =', k_a, k_b
            write(*,*) '  r_a, r_b   =', s%r(k_a), s%r(k_b)
            write(*,*) '  r_ob, r_cb =', r_ob, r_cb
            write(*,*) '  Hp_cb      =', Hp_cb
        end if

    end subroutine extended_convective_penetration

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Adjust temperature gradient in overshoot zone to be adiabatic (Pe>1d2) or radiative (Pe<1d-2) based upon the Peclet number,
    ! with a gradual transition between the two regimes.
    ! Only works if conv_premix = .true. since the last iteration in a time step has NaNs in D_mix if conv_premix = .false.
    subroutine other_adjust_mlt_gradT_fraction_Peclet(id, ierr)
        use utils_lib, only: is_bad
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type(star_info), pointer :: s
        real(dp) :: fraction, Peclet_number       ! f is fraction to compose grad_T = f*grad_ad + (1-f)*grad_rad
        integer :: k
        logical, parameter :: DEBUG = .FALSE.

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        if(is_bad(s% D_mix(1))) return

        if (s%num_conv_boundaries < 1) then ! Is zero at initialisation of the run
           if (DEBUG) then
              write(*,*) 'runstarex_gradT: skip since there are no convective boundaries'
           end if
           return
        endif

        do k= s%nz, 1, -1
            if (s%D_mix(k) <= s% min_D_mix) exit

            call calc_Peclet_number(id, k, Peclet_number, ierr)

            if (Peclet_number >= 100.0_dp) then
                fraction = 1.0_dp
            else if (Peclet_number .le. 0.01_dp) then
                fraction = 0.0_dp
            else
                fraction = (safe_log10(Peclet_number)+2.0_dp)/4.0_dp
            end if

            s% adjust_mlt_gradT_fraction(k) = fraction
        end do

    end subroutine other_adjust_mlt_gradT_fraction_Peclet

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Calculate the Peclet number for the specified cell
    subroutine calc_Peclet_number(id, k, Peclet_nr, ierr)
        integer, intent(in) :: id, k
        real(dp), intent(out) :: Peclet_nr
        integer, intent(out) :: ierr
        type(star_info), pointer :: s
        real(dp) :: conductivity, Hp

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        conductivity = 16.0_dp * boltz_sigma * pow3(s% T(k)) / ( 3.0_dp * s% opacity(k) * pow2(s% rho(k)) * s% cp(k) )
        call evaluate_Hp (s, k, s%r(k), Hp, ierr)
        Peclet_nr = s% conv_vel(k) * Hp * s% mixing_length_alpha / conductivity
    end subroutine calc_Peclet_number

    ! Evaluate the pressure scale height at a given location in the star
    subroutine evaluate_Hp (s, k, r, Hp, ierr)
        type(star_info), pointer :: s
        integer, intent(in)      :: k
        real(dp), intent(in)     :: r
        real(dp), intent(out)    :: Hp
        integer, intent(out)     :: ierr
        real(dp) :: P, rho
        ierr = 0

        P = exp(s%lnP(k))
        rho = exp(s%lnd(k))
        Hp = P/(rho*s%cgrav(k)* (s%M_center + s%xmstar*s%q(k))/(r*r))

    end subroutine evaluate_Hp

    end module run_star_extras
