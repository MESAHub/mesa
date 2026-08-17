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
      use auto_diff

      implicit none

      integer :: time0, time1, clock_rate, num_breakouts
      logical :: currently_in_breakout

      real(dp) :: injection_eps, injection_L, run_total_e


      contains


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
         s% other_energy => energy_dep
      end subroutine extras_controls


      subroutine energy_dep(id, ierr)
         integer, intent(in)  :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, nz
         real(dp) :: r_inject
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         injection_eps = 0d0
         injection_L = 0d0
         if (s% time >= s% x_ctrl(2) .or. s% x_ctrl(1) <= 0d0) then ! finished injection
            if (s% max_timestep /= s% x_ctrl(4)) then
               s% max_timestep = s% x_ctrl(4)
               write(*,1) 'overwrite max_timestep, log10', &
                  s% max_timestep, safe_log10(s% max_timestep)
            end if
            return
         end if
         nz = s% nz
         r_inject = s% x_ctrl(5)
         do k=1,nz-1
            if (s% r(k+1) < r_inject .and. r_inject <= s% r(k)) then
               injection_eps = s% x_ctrl(1)/s% dm(k)
               injection_L = s% x_ctrl(1)
               s% extra_heat(k) = injection_eps
               exit
            end if
         end do
         if (injection_L == 0d0) then
            write(*,*) 'failed to find specified r_inject in model', r_inject
            ierr = -1
            return
         end if
      end subroutine energy_dep


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call system_clock(time0,clock_rate)
         injection_eps = 0d0
         injection_L = 0d0
         if (.not. restart) then
            run_total_e = 0d0
            num_breakouts = 0
            currently_in_breakout = .false.
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
         if (s% x_logical_ctrl(1)) call switch_BCs(s)
         if (s% use_other_energy .and. s% x_ctrl(2) > 0d0) then
            s% max_timestep = s% x_ctrl(3)
            write(*,1) 'overwrite max_timestep, log10', &
               s% max_timestep, safe_log10(s% max_timestep)
         end if
      end subroutine extras_startup


      subroutine switch_BCs(s)
         type (star_info), pointer :: s
         if (s% x_logical_ctrl(2)) then
            s% use_compression_outer_BC = .false.
            s% use_momentum_outer_BC = .false.
            s% use_fixed_vsurf_outer_BC = .true.
         else
            s% use_compression_outer_BC = .true.
            s% use_momentum_outer_BC = .false.
            s% use_fixed_vsurf_outer_BC = .false.
         end if
      end subroutine switch_BCs


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. s% x_logical_ctrl(2)) return
         call system_clock(time1,clock_rate)
         dt = dble(time1 - time0) / clock_rate / 60
         if (.not. currently_in_breakout .and. num_breakouts == 3 .and. &
             abs(s% time - s% max_age_in_seconds) < 1d-6) then
            write(*,*) 'did 3 shock breakouts'
         else
            write(*,*) 'ERROR: failed to do 3 shock breakouts'
         end if
         write(*,'(/,a50,f12.2,99i10/)') 'runtime (minutes), retries, steps', &
            dt, s% num_retries, s% model_number
         ierr = 0
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
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
         if (s% use_other_energy) then
            how_many_extra_history_columns = 34
         else if (s% x_logical_ctrl(1)) then ! doing inlist_final
            how_many_extra_history_columns = 31
         else
            how_many_extra_history_columns = 0
         end if
      end function how_many_extra_history_columns

      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp), pointer :: v(:)
         real(dp) :: r_peak, tol
         real(dp) :: dr_peak, energy_result
         integer :: k_peak, k_front, k_back, k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. (s% x_logical_ctrl(1) .or. s% use_other_energy)) return
         if (s% u_flag) then
            v => s% u
         else
            v => s% v
         end if
         k_peak = maxloc(v(1:s% nz), dim=1)
         r_peak = s% r(k_peak)
         names(1) = 'r_peak_u'
         names(2) = 'r_peak'
         names(3) = 'r_peak_cs'
         names(4) = 'r_peak_u_div_cs'
         names(5) = 'r_peak_log'
         names(6) = 'r_peak_rho'
         names(7) = 'r_peak_T'
         names(8) = 'r_peak_gamma1'
         names(9) = 'r_peak_P'
         names(10) = 'r_peak_entropy'
         names(11) = 'r_peak_Cp'
         names(12) = 'r_peak_Cv'
         names(13) = 'r_peak_L'
         names(14) = 'r_peak_m'
         vals(1) = v(k_peak)
         vals(2) = r_peak
         vals(3) = s% csound(k_peak)
         vals(4) = v(k_peak)/s% csound(k_peak)
         vals(5) = safe_log10(r_peak)
         vals(6) = s% rho(k_peak)
         vals(7) = s% T(k_peak)
         vals(8) = s% gamma1(k_peak)
         vals(9) = s% Peos(k_peak)
         vals(10) = s% entropy(k_peak)*kerg*avo
         vals(11) = s% Cp(k_peak)
         vals(12) = s% Cv(k_peak)
         vals(13) = s% L(k_peak)
         vals(14) = s% m(k_peak)


         tol = s% x_ctrl(21)
         k_back = k_peak
         do while ((v(k_back) > v(k_peak)*tol) .and. (k_back < s% nz))
            k_back = k_back + 1
         end do
         k_front = k_peak
         do while (v(k_front) > v(k_peak)*tol .and. (k_front> 1))
             k_front = k_front - 1
         enddo

         dr_peak = s% r(k_front) - s% r(k_back)
         energy_result = 0d0
         do k = k_front, k_back
            energy_result = energy_result + 0.5*v(k)*v(k)* s% dm(k)
         end do
         names(15) = 'pulse_dr'
         vals(15) = dr_peak
         names(16) = 'pulse_KE'
         vals(16) = energy_result ! = 2 Pi r^2 rho v^2 dr !1/2 v^2 dm

         names(17) = 'pre_shock_u'
         names(18) = 'r_front'
         names(19) = 'k_peak_minus_k_front'
         names(20) = 'pre_shock_cs'
         names(21) = 'pre_shock_m'
         names(22) = 'pre_shock_rho'
         names(23) = 'pre_shock_T'
         names(24) = 'pre_shock_gamma1'
         names(25) = 'pre_shock_P'
         names(26) = 'pre_shock_entropy'
         names(27) = 'pre_shock_Cp'
         names(28) = 'pre_shock_Cv'
         names(29) = 'pre_shock_L'

         names(30) = 'r_front_minus_r_peak'
         names(31) = 'shock_deltaS'

         vals(17) = v(k_front)
         vals(18) = s%r(k_front)
         vals(19) = k_peak - k_front
         vals(20) = s% csound(k_front)
         vals(21) = s% m(k_front)
         vals(22) = s% rho(k_front)
         vals(23) = s% T(k_front)
         vals(24) = s% gamma1(k_front)
         vals(25) = s% Peos(k_front)
         vals(26) = s% entropy(k_front)*kerg*avo
         vals(27) = s% Cp(k_front)
         vals(28) = s% Cv(k_front)
         vals(29) = s% L(k_front)
         vals(30) = s% r(k_front) - s% r(k_peak)  ! This is v_peak/c_sound_shock_front
         vals(31) = s% entropy(k_peak) - s% entropy(k_front)

         if (.not. s% use_other_energy) return
         names(32) = 'heating_L_div_Lsun'
         vals(32) = injection_L/Lsun
         names(33) = 'heating_run_etot'
         vals(33) = run_total_e
         names(34) = 'heating_eps'
         vals(34) = injection_eps

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
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         integer :: k, k0, k1
         real(dp) :: v_esc
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = keep_going
         if (.not. s% x_logical_ctrl(1) .or. .not. s% u_flag) return
         k0 = s% nz+1
         do k = 1, s% nz
            v_esc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
            if (s% q(k) < s% x_ctrl(19)) exit ! only check outer layer
            if (s% u(k) > v_esc) then
               k0 = k
               exit
            end if
         end do
         if (k0 >= s% nz) return
         ! k0 is location below surface with u > escape velocity.
         ! remove inward from k0 where u large enough compared to escape velocity.
         k1 = k0
         do k = k0+1, s% nz
            v_esc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
            if (s% u(k) < s% x_ctrl(17)*v_esc) then ! stop removing here
               k1 = k-1
               exit
            end if
         end do
         if (s% q(k1) > s% x_ctrl(18)) return
         k1 = min(s% nz - 1, max(k1, s% x_integer_ctrl(20)))
         do while (k1 < s% nz)
            if (s% L(k1) > 0d0) exit
            k1 = k1 + 1
         end do
         if (k1 >= s% nz .or. s% L(k1) <= 0d0) return
         write(*,3) 'eject: model, k, delta, M_new, M_old', s% model_number, k1, &
            (s% m(1) - s% m(k1))/Msun, sum(s% dm(1:k1))/Msun
         extras_start_step = terminate
         call star_remove_surface_at_cell_k(s% id, k1, ierr)
         if (ierr /= 0) return
         if (s% x_logical_ctrl(1)) call switch_BCs(s)
         extras_start_step = keep_going
      end function extras_start_step


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         integer :: k_peak
         real(dp) :: r_peak
         real(dp), pointer :: v(:)
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         ! Count injected energy once per accepted step, not for each other_energy call.
         if (s% use_other_energy) run_total_e = run_total_e + s% total_extra_heating
         if (s% x_logical_ctrl(1)) then
            call switch_BCs(s)
            if (currently_in_breakout) then
               if (s% u(1) < s% csound(1)) currently_in_breakout = .false.
            else ! not currently_in_breakout
               if (s% u(1) >= s% csound(1)) then
                  num_breakouts = num_breakouts + 1
                  currently_in_breakout = .true.
               end if
            end if
         end if

         if (s% u_flag) then
            v => s% u
         else
            v => s% v
         end if
         k_peak = maxloc(v(1:s% nz), dim=1)
         r_peak = s% r(k_peak)

         if ((r_peak > s% x_ctrl(22)) .and. (s% x_logical_ctrl(2))) then
            write(*,*) 'terminating because r_peak > r_peak_max'
            extras_finish_step = terminate
            call store_extra_info(s)
            return
         end if
         call store_extra_info(s)
      end function extras_finish_step


      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         call move_flg(currently_in_breakout)
         call move_int(num_breakouts)
         num_ints = i

         i = 0
         ! call move_dbl
         call move_dbl(run_total_e)

         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            call mesa_error(__FILE__,__LINE__)
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras
