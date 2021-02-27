! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module timestep

      use star_private_def
      use utils_lib, only: is_bad
      use const_def
      use chem_def


      implicit none

      private
      public :: timestep_controller, check_change, check_integer_limit

      logical, parameter :: dbg_timestep = .false.
      real(dp) :: max_dt
      integer :: why_Tlim

      contains


      integer function timestep_controller(s, max_timestep)
         ! if don't return keep_going, then set result_reason to say why.
         type (star_info), pointer :: s
         real(dp), intent(in) :: max_timestep

         include 'formats'

         timestep_controller = do_timestep_limits(s, s% dt)
         if (timestep_controller /= keep_going) s% result_reason = timestep_limits

         ! strictly enforce maximum timestep
         max_dt = max_timestep
         if (max_dt <= 0) max_dt = 1d99
         if (s% dt_next > max_dt) then
            s% dt_next = max_dt
            s% why_Tlim = Tlim_max_timestep
         end if

         if (s% timestep_hold > s% model_number .and. s% dt_next > s% dt) then
            s% dt_next = s% dt
            s% why_Tlim = Tlim_timestep_hold
            if (s% report_dt_hard_limit_retries .and. timestep_controller == keep_going) &
               write(*,3) 'timestep_hold > model_number, so no timestep increase', &
                  s% timestep_hold, s% model_number
         end if

         if (is_bad_num(s% dt_next)) then
            write(*, *) 'timestep_controller: dt_next', s% dt_next
            if (s% stop_for_bad_nums) stop 'timestep_controller'
            timestep_controller = terminate
            s% termination_code = t_timestep_controller
            return
         end if

         if (dbg_timestep) then
            write(*,*) 'final result from timestep_controller for model_number', &
               timestep_controller, s% model_number
            write(*,1) 'lg dt/secyer', log10(s% dt/secyer)
            write(*,1) 'lg dt_next/secyer', log10(s% dt_next/secyer)
            write(*,1) 'dt_next/dt', s% dt_next/s% dt
            write(*,*)
            write(*,*)
            stop 'timestep_controller'
         end if

      end function timestep_controller


      integer function do_timestep_limits(s, dt)
         use rates_def, only: i_rate
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt ! timestep just completed

         real(dp) :: dt_limit_ratio(numTlim), dt_limit, limit, order, max_timestep_factor
         integer :: i_limit, nz, ierr
         logical :: skip_hard_limit
         integer :: num_mix_boundaries ! boundaries of regions with mixing_type /= no_mixing
         real(dp), pointer :: mix_bdy_q(:) ! (num_mix_boundaries)
         integer, pointer :: mix_bdy_loc(:) ! (num_mix_boundaries)

         include 'formats'

         if (s% never_skip_hard_limits) then
            skip_hard_limit = .false.
         else
            skip_hard_limit = (s% timestep_hold >= s% model_number) .or. &
               (s% relax_hard_limits_after_retry .and. &
                s% model_number_for_last_retry == s% model_number)
         end if

         nz = s% nz

         ! NOTE: when we get here, complete_model has called the report routine,
         ! so we can use information that it has calculated

         ierr = 0

         num_mix_boundaries = s% num_mix_boundaries
         mix_bdy_q => s% mix_bdy_q
         mix_bdy_loc => s% mix_bdy_loc

         dt_limit_ratio(:) = 0d0

         do_timestep_limits = check_varcontrol_limit( &
            s, dt_limit_ratio(Tlim_struc))
         if (return_now(Tlim_struc)) return

         if (.not. s% doing_first_model_of_run) then
         
            if (s% use_other_timestep_limit) then
               do_timestep_limits = s% other_timestep_limit( &
                  s% id, skip_hard_limit, dt, dt_limit_ratio(Tlim_other_timestep_limit))
               if (return_now(Tlim_other_timestep_limit)) return
            end if

            do_timestep_limits = check_solver_iters_timestep_limit( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_solver))
            if (return_now(Tlim_solver)) return

            do_timestep_limits = check_burn_steps_limit( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_num_burn_steps))
            if (return_now(Tlim_num_burn_steps)) return

            do_timestep_limits = check_diffusion_steps_limit( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_num_diff_solver_steps))
            if (return_now(Tlim_num_diff_solver_steps)) return

            do_timestep_limits = check_diffusion_iters_limit( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_num_diff_solver_iters))
            if (return_now(Tlim_num_diff_solver_iters)) return

            do_timestep_limits = check_dX(s, 0, skip_hard_limit, dt, &
               num_mix_boundaries, mix_bdy_loc, mix_bdy_q, &
               dt_limit_ratio(Tlim_dH), dt_limit_ratio(Tlim_dH_div_H))
            if (return_now(Tlim_dH_div_H)) return

            do_timestep_limits = check_dX(s, 1, skip_hard_limit, dt, &
               num_mix_boundaries, mix_bdy_loc, mix_bdy_q, &
               dt_limit_ratio(Tlim_dHe), dt_limit_ratio(Tlim_dHe_div_He))
            if (return_now(Tlim_dHe_div_He)) return

            do_timestep_limits = check_dX(s, 2, skip_hard_limit, dt, &
               num_mix_boundaries, mix_bdy_loc, mix_bdy_q, &
               dt_limit_ratio(Tlim_dHe3), dt_limit_ratio(Tlim_dHe3_div_He3))
            if (return_now(Tlim_dHe3_div_He3)) return

            do_timestep_limits = check_dX(s, -1, skip_hard_limit, dt, &
               num_mix_boundaries, mix_bdy_loc, mix_bdy_q, &
               dt_limit_ratio(Tlim_dX), dt_limit_ratio(Tlim_dX_div_X))
            if (return_now(Tlim_dX_div_X)) return

            do_timestep_limits = check_dL_div_L( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dL_div_L))
            if (return_now(Tlim_dL_div_L)) return

            do_timestep_limits = check_dlgP_change( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_dlgP))
            if (return_now(Tlim_dlgP)) return

            do_timestep_limits = check_dlgRho_change( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_dlgRho))
            if (return_now(Tlim_dlgRho)) return

            do_timestep_limits = check_dlgT_change( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_dlgT))
            if (return_now(Tlim_dlgT)) return

            do_timestep_limits = check_dlgE_change( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_dlgE))
            if (return_now(Tlim_dlgE)) return

            do_timestep_limits = check_dlgR_change( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_dlgR))
            if (return_now(Tlim_dlgR)) return

            do_timestep_limits = check_lgL_nuc_cat_change( &
               s, num_mix_boundaries, mix_bdy_q, &
               skip_hard_limit, dt_limit_ratio(Tlim_dlgL_nuc_cat))
            if (return_now(Tlim_dlgL_nuc_cat)) return

            do_timestep_limits = check_lgL_H_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgL_H))
            if (return_now(Tlim_dlgL_H)) return

            do_timestep_limits = check_lgL_He_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgL_He))
            if (return_now(Tlim_dlgL_He)) return

            do_timestep_limits = check_lgL_z_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgL_z))
            if (return_now(Tlim_dlgL_z)) return

            do_timestep_limits = check_lgL_power_photo_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgL_power_photo))
            if (return_now(Tlim_dlgL_power_photo)) return

            do_timestep_limits = check_lgL_nuc_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgL_nuc))
            if (return_now(Tlim_dlgL_nuc)) return

            do_timestep_limits = check_dlgTeff_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgTeff))
            if (return_now(Tlim_dlgTeff)) return

            do_timestep_limits = check_dlgRho_cntr_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgRho_cntr))
            if (return_now(Tlim_dlgRho_cntr)) return

            do_timestep_limits = check_dlgT_cntr_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgT_cntr))
            if (return_now(Tlim_dlgT_cntr)) return

            do_timestep_limits = check_dlgT_max_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgT_max))
            if (return_now(Tlim_dlgT_max)) return

            do_timestep_limits = check_dlgT_max_at_high_T_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgT_max_at_high_T))
            if (return_now(Tlim_dlgT_max_at_high_T)) return

            do_timestep_limits = check_dlgP_cntr_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlgP_cntr))
            if (return_now(Tlim_dlgP_cntr)) return

            do_timestep_limits = check_dYe_highT_change( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_delta_Ye_highT))
            if (return_now(Tlim_delta_Ye_highT)) return

            do_timestep_limits = check_dlog_eps_nuc_change( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dlog_eps_nuc))
            if (return_now(Tlim_dlog_eps_nuc)) return

            do_timestep_limits = check_dX_div_X_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_dX_div_X_cntr))
            if (return_now(Tlim_dX_div_X_cntr)) return

            do_timestep_limits = check_lg_XH_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_lg_XH_cntr))
            if (return_now(Tlim_lg_XH_cntr)) return

            do_timestep_limits = check_lg_XHe_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_lg_XHe_cntr))
            if (return_now(Tlim_lg_XHe_cntr)) return

            do_timestep_limits = check_lg_XC_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_lg_XC_cntr))
            if (return_now(Tlim_lg_XC_cntr)) return

            do_timestep_limits = check_lg_XNe_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_lg_XNe_cntr))
            if (return_now(Tlim_lg_XNe_cntr)) return

            do_timestep_limits = check_lg_XO_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_lg_XO_cntr))
            if (return_now(Tlim_lg_XO_cntr)) return

            do_timestep_limits = check_lg_XSi_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_lg_XSi_cntr))
            if (return_now(Tlim_lg_XSi_cntr)) return

            do_timestep_limits = check_XH_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_XH_cntr))
            if (return_now(Tlim_XH_cntr)) return

            do_timestep_limits = check_XHe_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_XHe_cntr))
            if (return_now(Tlim_XHe_cntr)) return

            do_timestep_limits = check_XC_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_XC_cntr))
            if (return_now(Tlim_XC_cntr)) return

            do_timestep_limits = check_XNe_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_XNe_cntr))
            if (return_now(Tlim_XNe_cntr)) return

            do_timestep_limits = check_XO_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_XO_cntr))
            if (return_now(Tlim_XO_cntr)) return

            do_timestep_limits = check_XSi_cntr( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_XSi_cntr))
            if (return_now(Tlim_XSi_cntr)) return

            do_timestep_limits = check_delta_mstar( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dmstar))
            if (return_now(Tlim_dmstar)) return

            do_timestep_limits = check_delta_mdot( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_del_mdot))
            if (return_now(Tlim_del_mdot)) return

            do_timestep_limits = check_adjust_J_q( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_adjust_J_q))
            if (return_now(Tlim_adjust_J_q)) return

            do_timestep_limits = check_delta_lgL( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_lgL))
            if (return_now(Tlim_lgL)) return

            do_timestep_limits = check_dt_div_dt_cell_collapse( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dt_div_dt_cell_collapse))
            if (return_now(Tlim_dt_div_dt_cell_collapse)) return

            do_timestep_limits = check_dt_div_min_dr_div_cs( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dt_div_min_dr_div_cs))
            if (return_now(Tlim_dt_div_min_dr_div_cs)) return

            do_timestep_limits = check_delta_HR( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_delta_HR))
            if (return_now(Tlim_delta_HR)) return

            do_timestep_limits = check_rel_error_in_energy( &
               s, skip_hard_limit, dt_limit_ratio(Tlim_error_in_energy_conservation))
            if (return_now(Tlim_error_in_energy_conservation)) return

            do_timestep_limits = check_dX_nuc_drop( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_dX_nuc_drop))
            if (return_now(Tlim_dX_nuc_drop)) return

         end if

         i_limit = maxloc(dt_limit_ratio(1:numTlim), dim=1)
         
         order = 1
         call filter_dt_next(s, order, dt_limit_ratio(i_limit)) ! sets s% dt_next
         
         if (s% log_max_temperature > s% min_logT_for_max_timestep_factor_at_high_T) then
            max_timestep_factor = s% max_timestep_factor_at_high_T
         else
            max_timestep_factor = s% max_timestep_factor
         end if
         if (max_timestep_factor > 0 .and. s% dt_next > max_timestep_factor*s% dt) then
            s% dt_next = max_timestep_factor*s% dt
            if (i_limit == Tlim_struc) i_limit = Tlim_max_timestep_factor
         end if

         if (s% min_timestep_factor > 0 .and. s% dt_next < s% min_timestep_factor*s% dt) then
            s% dt_next = s% min_timestep_factor*s% dt
            if (i_limit == Tlim_struc) i_limit = Tlim_min_timestep_factor
         end if

         s% why_Tlim = i_limit
         if (i_limit > 0) s% dt_why_count(i_limit) = s% dt_why_count(i_limit) + 1

         contains

         logical function return_now(i_limit)
            integer, intent(in) :: i_limit
            integer :: k
            if (do_timestep_limits == keep_going) then
               return_now = .false.
               return
            end if
            return_now = .true.
            s% why_Tlim = i_limit
            if (i_limit > 0) s% dt_why_retry_count(i_limit) = s% dt_why_retry_count(i_limit) + 1
         end function return_now

      end function do_timestep_limits


      integer function check_integer_limit( &
            s, limit, hard_limit, value, msg, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         integer, intent(in) :: limit, hard_limit, value
         character (len=*), intent(in) :: msg
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         include 'formats'
         if (value > hard_limit .and. hard_limit > 0 .and. (.not. skip_hard_limit)) then
            check_integer_limit = retry
            s% retry_message = trim(msg) // ' hard limit'
            if (s% report_dt_hard_limit_retries) then
               write(*,*) trim(msg) // ' hard limit', hard_limit, value
               write(*,3) trim(msg), s% model_number, value
               write(*,3) trim(msg) // ' hard limit', s% model_number, hard_limit
            end if
            return
         end if
         check_integer_limit = keep_going
         if (value <= 0 .or. limit <= 0) return
         dt_limit_ratio = dble(value)/dble(limit) - 0.05d0
            ! subtract a bit so that allow dt to grow even if value == limit
         if (dt_limit_ratio <= 1d0) dt_limit_ratio = 0
      end function check_integer_limit


      integer function check_burn_steps_limit(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         integer :: max_steps, i
         check_burn_steps_limit = keep_going
         if (.not. s% op_split_burn .or. maxval(s% T_start(1:s%nz)) < s% op_split_burn_min_T) return
         max_steps = maxval(s% burn_num_iters(1:s% nz))
         check_burn_steps_limit = check_integer_limit( &
           s, s% burn_steps_limit, s% burn_steps_hard_limit, max_steps,  &
           'num_burn_solver_steps', skip_hard_limit, dt, dt_limit_ratio)
      end function check_burn_steps_limit


      integer function check_solver_iters_timestep_limit(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         integer :: iters, limit, hard_limit
         include 'formats'
         check_solver_iters_timestep_limit = keep_going
         iters = s% num_solver_iterations
         limit = s% solver_iters_timestep_limit
         hard_limit = -1
         if (s% using_gold_tolerances) then
            limit = s% gold_solver_iters_timestep_limit
         end if
         if (s% used_extra_iter_in_solver_for_accretion) iters = iters - 1
         check_solver_iters_timestep_limit = check_integer_limit( &
           s, limit, hard_limit, &
           iters, 'num_solver_iterations', skip_hard_limit, dt, dt_limit_ratio)
      end function check_solver_iters_timestep_limit


      integer function check_diffusion_steps_limit(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         check_diffusion_steps_limit = keep_going
         if (.not. s% do_element_diffusion) return
         check_diffusion_steps_limit = check_integer_limit( &
           s, s% diffusion_steps_limit, s% diffusion_steps_hard_limit, &
           s% num_diffusion_solver_steps,  &
           'num_diffusion_solver_steps', skip_hard_limit, dt, dt_limit_ratio)
      end function check_diffusion_steps_limit


      integer function check_diffusion_iters_limit(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         check_diffusion_iters_limit = keep_going
         if (.not. s% do_element_diffusion) return
         check_diffusion_iters_limit = check_integer_limit( &
           s, s% diffusion_iters_limit, s% diffusion_iters_hard_limit, &
           s% num_diffusion_solver_iters,  &
           'num_diffusion_solver_iters', skip_hard_limit, dt, dt_limit_ratio)
      end function check_diffusion_iters_limit


      integer function check_dX(s, which, skip_hard_limit, dt, &
            n_mix_bdy, mix_bdy_loc, mix_bdy_q, &
            dX_dt_limit_ratio, dX_div_X_dt_limit_ratio)
         use num_lib, only: binary_search
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         integer, intent(in) :: which, n_mix_bdy, mix_bdy_loc(:)
         real(dp), intent(in) :: dt
         real(dp), intent(in), pointer :: mix_bdy_q(:)
         real(dp), intent(inout) :: dX_dt_limit_ratio, dX_div_X_dt_limit_ratio

         real(dp) :: X, X_old, delta_dX, delta_dX_div_X, max_dX, max_dX_div_X, &
            bdy_dist_dm, max_dX_bdy_dist_dm, max_dX_div_X_bdy_dist_dm, cz_dist_limit
         integer :: j, k, cid, bdy, max_dX_j, max_dX_k, max_dX_div_X_j, max_dX_div_X_k
         real(dp) :: D_mix_cutoff, dX_limit_min_X, dX_limit, dX_hard_limit
         real(dp) :: dX_div_X_limit_min_X, dX_div_X_limit, dX_div_X_hard_limit
         real(dp) :: dX_div_X_at_high_T_limit, dX_div_X_at_high_T_hard_limit
         real(dp) :: dX_div_X_at_high_T_limit_lgT_min
         logical :: decreases_only

         include 'formats'

         check_dX = keep_going

         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return

         if (which == 0) then ! hydrogen
            dX_limit_min_X = s% dH_limit_min_H
            dX_limit = s% dH_limit
            dX_hard_limit = s% dH_hard_limit
            dX_div_X_limit_min_X = s% dH_div_H_limit_min_H
            dX_div_X_limit = s% dH_div_H_limit
            dX_div_X_hard_limit = s% dH_div_H_hard_limit
            decreases_only = s% dH_decreases_only
         else if (which == 1) then ! helium
            dX_limit_min_X = s% dHe_limit_min_He
            dX_limit = s% dHe_limit
            dX_hard_limit = s% dHe_hard_limit
            dX_div_X_limit_min_X = s% dHe_div_He_limit_min_He
            dX_div_X_limit = s% dHe_div_He_limit
            dX_div_X_hard_limit = s% dHe_div_He_hard_limit
            decreases_only = s% dHe_decreases_only
         else if (which == 2) then ! He3
            dX_limit_min_X = s% dHe3_limit_min_He3
            dX_limit = s% dHe3_limit
            dX_hard_limit = s% dHe3_hard_limit
            dX_div_X_limit_min_X = s% dHe3_div_He3_limit_min_He3
            dX_div_X_limit = s% dHe3_div_He3_limit
            dX_div_X_hard_limit = s% dHe3_div_He3_hard_limit
            decreases_only = s% dHe3_decreases_only
         else ! metals
            dX_limit_min_X = s% dX_limit_min_X
            dX_limit = s% dX_limit
            dX_hard_limit = s% dX_hard_limit
            dX_div_X_limit_min_X = s% dX_div_X_limit_min_X
            if (s% log_max_temperature > s% dX_div_X_at_high_T_limit_lgT_min) then
               dX_div_X_limit = s% dX_div_X_at_high_T_limit
               dX_div_X_hard_limit = s% dX_div_X_at_high_T_hard_limit
            else
               dX_div_X_limit = s% dX_div_X_limit
               dX_div_X_hard_limit = s% dX_div_X_hard_limit
            end if
            decreases_only = s% dX_decreases_only
         end if
         
         dX_limit = dX_limit*s% time_delta_coeff
         dX_hard_limit = dX_hard_limit*s% time_delta_coeff

         dX_div_X_limit = dX_div_X_limit*s% time_delta_coeff
         dX_div_X_hard_limit = dX_div_X_hard_limit*s% time_delta_coeff

         if (  dX_limit_min_X >= 1 .and. &
               dX_limit >= 1 .and. &
               dX_hard_limit >= 1 .and. &
               dX_div_X_limit_min_X >= 1 .and. &
               dX_div_X_limit >= 1 .and. &
               dX_div_X_hard_limit >= 1) then
            return
         end if

         max_dX = -1; max_dX_j = -1; max_dX_k = -1
         max_dX_div_X = -1; max_dX_div_X_j = -1; max_dX_div_X_k = -1
         bdy = 0
         max_dX_bdy_dist_dm = 0
         max_dX_div_X_bdy_dist_dm = 0
         cz_dist_limit = s% dX_mix_dist_limit*Msun

         if (s% set_min_D_mix .and. s% ye(s% nz) >= s% min_center_Ye_for_min_D_mix) then
            D_mix_cutoff = s% min_D_mix
         else
            D_mix_cutoff = 0
         end if

         do k = 1, s% nz

            if (s% D_mix(k) > D_mix_cutoff) then
               cycle
            end if
            if (k < s% nz) then
               if (s% D_mix(k+1) > D_mix_cutoff) then
                  cycle
               end if
            end if

            ! find the nearest mixing boundary
            bdy = binary_search(n_mix_bdy, mix_bdy_q, bdy, s% q(k))
            ! don't check cells near a mixing boundary
            if (bdy > 0 .and. bdy < n_mix_bdy) then
               bdy_dist_dm = s% xmstar*abs(s% q(k) - mix_bdy_q(bdy))
               if (bdy_dist_dm < cz_dist_limit) cycle
            else
               bdy_dist_dm = 0
            end if

            do j = 1, s% species

               cid = s% chem_id(j)
               if (which == 0) then ! hydrogen
                  if (cid /= ih1) cycle
               else if (which == 1) then ! helium
                  if (cid /= ihe4) cycle
               else if (which == 2) then ! he3
                  if (cid /= ihe3) cycle
               else ! other
                  if (chem_isos% Z(cid) <= 2) cycle
               end if

               X = s% xa(j,k)
               X_old = s% xa_old(j,k)
               delta_dX = X_old - X ! decrease in abundance

               if ((.not. decreases_only) .and. delta_dX < 0) delta_dX = -delta_dX

               if (X >= dX_limit_min_X) then
                  if ((.not. skip_hard_limit) .and. delta_dX > dX_hard_limit) then
                     check_dX = retry
                     s% why_Tlim = Tlim_dX
                     s% Tlim_dX_species = j
                     s% Tlim_dX_cell = k
                     s% retry_message = 'dX ' // trim(chem_isos% name(s% chem_id(j))) // ' hard limit'
                     s% retry_message_k = k
                     if (s% report_dt_hard_limit_retries) then
                        write(*,2) 'old xa', s% model_number, X_old
                        write(*,2) 'new xa', s% model_number, X
                        write(*,2) 'delta xa', s% model_number, delta_dX
                        write(*,2) 'hard limit delta xa', s% model_number, dX_hard_limit
                     end if
                     return
                  end if
                  if (delta_dX > max_dX) then
                     max_dX = delta_dX
                     max_dX_j = j
                     max_dX_k = k
                     max_dX_bdy_dist_dm = bdy_dist_dm
                  end if
               end if
               if (X >= dX_div_X_limit_min_X) then
                  delta_dX_div_X = delta_dX/X
                  if ((.not. skip_hard_limit) .and. delta_dX_div_X > dX_div_X_hard_limit) then
                     check_dX = retry
                     s% why_Tlim = Tlim_dX_div_X
                     s% Tlim_dX_div_X_species = j
                     s% Tlim_dX_div_X_cell = k            
                     s% retry_message = 'dX_div_X ' // trim(chem_isos% name(s% chem_id(j))) // ' hard limit'
                     s% retry_message_k = k
                     if (s% report_dt_hard_limit_retries) then
                        write(*, '(a30, i5, 99(/,a30,e20.10))') &
                           'delta_dX_div_X ' // trim(chem_isos% name(s% chem_id(j))), &
                           k, 'delta_dX_div_X', delta_dX_div_X, &
                           'dX_div_X_hard_limit', dX_div_X_hard_limit
                        write(*,2) 'old xa', s% model_number, X_old
                        write(*,2) 'new xa', s% model_number, X
                        write(*,2) 'delta_dX_div_X', s% model_number, delta_dX_div_X
                        write(*,2) 'dX_div_X_hard_limit', s% model_number, dX_div_X_hard_limit
                     end if
                     return
                  end if
                  if (delta_dX_div_X > max_dX_div_X) then
                     max_dX_div_X = delta_dX_div_X
                     max_dX_div_X_j = j
                     max_dX_div_X_k = k
                     max_dX_div_X_bdy_dist_dm = bdy_dist_dm
                  end if
               end if
            end do
         end do

         if (dX_limit > 0) then
            dX_dt_limit_ratio = max_dX/dX_limit
            if (dX_dt_limit_ratio <= 1d0) then
               dX_dt_limit_ratio = 0
            else
               j = max_dX_j
               k = max_dX_k
               s% Tlim_dX_species = j
               s% Tlim_dX_cell = k
               write(*, '(a30, i5, 99e20.10)') &
                  'dX ' // trim(chem_isos% name(s% chem_id(j))), &
                  k, max_dX, dX_limit, (s% M_center + s% xmstar*(s% q(k) - s% dq(k)/2))/Msun, &
                  max_dX_bdy_dist_dm/Msun
            end if

         end if

         if (dX_div_X_limit > 0) then
            dX_div_X_dt_limit_ratio = max_dX_div_X/dX_div_X_limit
            if (dX_div_X_dt_limit_ratio <= 1d0) then
               dX_div_X_dt_limit_ratio = 0
            else
               s% Tlim_dX_div_X_species = max_dX_div_X_j
               s% Tlim_dX_div_X_cell = max_dX_div_X_k
               if (which > 2) write(*, '(a30, i5, 99e20.10)') &
                  'limit dt because of large dX_div_X ' // &
                     trim(chem_isos% name(s% chem_id(max_dX_div_X_j))) // &
                     ' k, max, lim, m ', &
                  max_dX_div_X_k, max_dX_div_X, dX_div_X_limit, &
                  max_dX_div_X_bdy_dist_dm/Msun
            end if
         end if

      end function check_dX


      integer function check_dL_div_L(s, skip_hard_limit, dt, dL_div_L_dt_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dL_div_L_dt_ratio

         real(dp) :: L, abs_dL, abs_dL_div_L, max_dL_div_L
         integer :: k, max_dL_div_L_k
         real(dp) :: dL_div_L_limit_min_L, dL_div_L_limit, dL_div_L_hard_limit
         
         include 'formats'

         check_dL_div_L = keep_going

         dL_div_L_limit_min_L = Lsun*s% dL_div_L_limit_min_L
         dL_div_L_limit = s% dL_div_L_limit*s% time_delta_coeff
         dL_div_L_hard_limit = s% dL_div_L_hard_limit*s% time_delta_coeff

         if (dL_div_L_limit_min_L <= 0) return
         if (dL_div_L_limit <= 0 .and. dL_div_L_hard_limit <= 0) return

         max_dL_div_L = -1
         max_dL_div_L_k = -1
         abs_dL_div_L = 0; L=0 ! to quiet gfortran

         do k = 1, s% nz
            L = s% L(k)
            abs_dL = abs(L - s% L_start(k))
            if (L >= dL_div_L_limit_min_L) then
               abs_dL_div_L = abs_dL/L
               if (dL_div_L_hard_limit > 0 .and. (.not. skip_hard_limit) &
                     .and. abs_dL_div_L > dL_div_L_hard_limit) then
                  check_dL_div_L= retry
                  s% retry_message = 'dL_div_L hard limit'
                  s% retry_message_k = k
                  if (s% report_dt_hard_limit_retries) then
                     write(*,2) 'L', L
                     write(*,2) 'L_start', s% L_start(k)
                     write(*,2) 'abs_dL_div_L', abs_dL_div_L
                     write(*,2) 'dL_div_L_hard_limit', dL_div_L_hard_limit
                  end if                  
                  return
               end if
               if (abs_dL_div_L > max_dL_div_L) then
                  max_dL_div_L = abs_dL_div_L
                  max_dL_div_L_k = k
               end if
            end if
         end do

         if (dL_div_L_limit <= 0) return
         dL_div_L_dt_ratio = max_dL_div_L/dL_div_L_limit

      end function check_dL_div_L


      integer function check_change( &
            s, delta_value, lim_in, hard_lim_in, max_k, msg, &
            skip_hard_limit, dt_limit_ratio, relative_excess)
         use const_def, only:ln10
         type (star_info), pointer :: s
         real(dp), intent(in) :: delta_value, lim_in, hard_lim_in
         integer, intent(in) :: max_k
         character (len=*), intent(in) :: msg
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp), intent(out) :: relative_excess
         real(dp) :: abs_change, lim, hard_lim
         include 'formats'
         if (is_bad(delta_value)) then
            write(*,1) trim(msg) // ' delta_value', delta_value
            stop 'check_change'
         end if
         check_change = keep_going
         abs_change = abs(delta_value)
         lim = lim_in*s% time_delta_coeff
         hard_lim = hard_lim_in*s% time_delta_coeff
         if (hard_lim > 0 .and. abs_change > hard_lim .and. (.not. skip_hard_limit)) then
            s% retry_message = trim(msg) // ' hard limit'
            s% retry_message_k = max_k
            check_change = retry
            return
         end if
         if (lim <= 0) return
         relative_excess = (abs_change - lim) / lim
         dt_limit_ratio = abs_change/lim ! 1d0/(s% timestep_dt_factor**relative_excess)
         if (is_bad(dt_limit_ratio)) then
            write(*,1) trim(msg) // ' dt_limit_ratio', dt_limit_ratio, abs_change, lim
            stop 'check_change'
         end if
         if (dt_limit_ratio <= 1d0) dt_limit_ratio = 0
      end function check_change


      subroutine get_dlgP_info(s, i, max_dlnP)
         use const_def, only:ln10
         type (star_info), pointer :: s
         integer, intent(out) :: i
         real(dp), intent(out) :: max_dlnP
         real(dp) :: lim, dlnP
         integer :: k
         include 'formats'
         lim = ln10*s% delta_lgP_limit_min_lgP
         i = 0
         max_dlnP = 0
         do k=1,s% nz
            if (s% lnP(k) < lim) cycle
            dlnP = abs(s% lnP(k) - s% lnP_start(k))
            if (dlnP > max_dlnP) then               
               max_dlnP = dlnP
               i = k
            end if
         end do
      end subroutine get_dlgP_info


      integer function check_dlgP_change(s, skip_hard_limit, dt_limit_ratio)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, max_dlnP
         integer :: i
         include 'formats'
         check_dlgP_change = keep_going
         call get_dlgP_info(s, i, max_dlnP)
         if (i == 0) return
         check_dlgP_change = check_change(s, max_dlnP/ln10, &
            s% delta_lgP_limit, s% delta_lgP_hard_limit, &
            i, 'check_dlgP_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgP_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,3) 'lgP', i, s% lnP(i)/ln10
            write(*,3) 'lgP_old', i, s% lnP_start(i)/ln10
            write(*,3) 'dlgP', i, (s% lnP(i) - s% lnP_start(i))/ln10
            write(*,3) 'hard_limit', i, s% delta_lgP_hard_limit
         end if
      end function check_dlgP_change


      subroutine get_dlgRho_info(s, i, max_dlnRho)
         use const_def, only:ln10
         type (star_info), pointer :: s
         integer, intent(out) :: i
         real(dp), intent(out) :: max_dlnRho
         real(dp) :: lim, dlnRho, max_abs_dlnRho
         integer :: k
         include 'formats'
         lim = ln10*s% delta_lgRho_limit_min_lgRho
         i = 0
         max_abs_dlnRho = 0
         do k=1,s% nz
            if (s% lnd(k) < lim) cycle
            dlnRho = s% lnd(k) - s% lnd_start(k)
            if (abs(dlnRho) > max_abs_dlnRho) then
               max_dlnRho = dlnRho
               max_abs_dlnRho = abs(dlnRho)
               i = k
            end if
         end do
      end subroutine get_dlgRho_info


      integer function check_dlgRho_change(s, skip_hard_limit, dt_limit_ratio)
         ! check max change in log10(density)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, max_dlnd
         integer :: i
         include 'formats'
         check_dlgRho_change = keep_going
         call get_dlgRho_info(s, i, max_dlnd)
         if (i == 0) return
         check_dlgRho_change = check_change(s, abs(max_dlnd)/ln10, &
            s% delta_lgRho_limit, s% delta_lgRho_hard_limit, &
            i, 'check_dlgRho_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgRho_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,3) 'lgRho', i, s% lnd(i)/ln10
            write(*,3) 'lgRho_old', i, s% lnd_start(i)/ln10
            write(*,3) 'dlgRho', i, (s% lnd(i) - s% lnd_start(i))/ln10
            write(*,3) 'hard_limit', i, s% delta_lgRho_hard_limit
         end if
      end function check_dlgRho_change


      subroutine get_dlgT_info(s, i, max_dlnT)
         use const_def, only:ln10
         type (star_info), pointer :: s
         integer, intent(out) :: i
         real(dp), intent(out) :: max_dlnT
         real(dp) :: lim, dlnT, abs_max_dlnT
         integer :: k
         include 'formats'
         lim = ln10*s% delta_lgT_limit_min_lgT
         i = 0
         abs_max_dlnT = 0
         do k=1,s% nz
            if (s% lnT(k) < lim) cycle
            dlnT = s% lnT(k) - s% lnT_start(k)
            if (abs(dlnT) > abs_max_dlnT) then
               max_dlnT = dlnT
               abs_max_dlnT = abs(dlnT)
               i = k
            end if
         end do
      end subroutine get_dlgT_info


      integer function check_dlgT_change(s, skip_hard_limit, dt_limit_ratio)
         ! check max change in log10(temperature)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, max_dlnT
         integer :: i
         include 'formats'
         check_dlgT_change = keep_going
         call get_dlgT_info(s, i, max_dlnT)
         if (i == 0) return
         check_dlgT_change = check_change(s, abs(max_dlnT)/ln10, &
            s% delta_lgT_limit, s% delta_lgT_hard_limit, &
            i, 'check_dlgT_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgT_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,3) 'lgT', i, s% lnT(i)/ln10
            write(*,3) 'lgT_old', i, s% lnT_start(i)/ln10
            write(*,3) 'dlgT', i, (s% lnT(i) - s% lnT_start(i))/ln10
            write(*,3) 'hard_limit', i, s% delta_lgT_hard_limit
         end if
      end function check_dlgT_change


      subroutine get_dlgE_info(s, i, max_dlnE)
         use const_def, only:ln10
         type (star_info), pointer :: s
         integer, intent(out) :: i
         real(dp), intent(out) :: max_dlnE
         real(dp) :: lim, dlnE
         integer :: k
         include 'formats'
         lim = ln10*s% delta_lgE_limit_min_lgE
         i = 0
         max_dlnE = 0
         do k=1,s% nz
            if (s% lnE(k) < lim) cycle
            dlnE = abs(s% energy(k) - s% energy_start(k))/s% energy(k)
            if (dlnE > max_dlnE) then
               max_dlnE = dlnE
               i = k
            end if
         end do
      end subroutine get_dlgE_info


      integer function check_dlgE_change(s, skip_hard_limit, dt_limit_ratio)
         ! check max change in log10(internal energy)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, max_dlnE
         integer :: i
         include 'formats'
         check_dlgE_change = keep_going
         call get_dlgE_info(s, i, max_dlnE)
         if (i == 0) return
         check_dlgE_change = check_change(s, max_dlnE/ln10, &
            s% delta_lgE_limit, s% delta_lgE_hard_limit, &
            i, 'check_dlgE_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgE_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,3) 'lgE', i, safe_log10(s% energy(i))
            write(*,3) 'lgE_old', i, safe_log10(s% energy_start(i))
            write(*,3) 'dlgE', i, (s% energy(i) - s% energy_start(i))/s% energy(i)/ln10
            write(*,3) 'hard_limit', i, s% delta_lgE_hard_limit
         end if
      end function check_dlgE_change


      subroutine get_dlgR_info(s, i, max_dlnR)
         use const_def, only:ln10
         type (star_info), pointer :: s
         integer, intent(out) :: i
         real(dp), intent(out) :: max_dlnR
         real(dp) :: lim, dlnR
         integer :: k
         include 'formats'
         lim = ln10*s% delta_lgR_limit_min_lgR
         i = 0
         max_dlnR = 0
         do k=1,s% nz
            if (s% lnR(k) < lim) cycle
            dlnR = abs(s% lnR(k) - s% lnR_start(k))
            if (dlnR > max_dlnR) then
               max_dlnR = dlnR
               i = k
            end if
         end do
      end subroutine get_dlgR_info


      integer function check_dlgR_change(s, skip_hard_limit, dt_limit_ratio)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, max_dlnR
         integer :: i
         include 'formats'
         check_dlgR_change = keep_going
         call get_dlgR_info(s, i, max_dlnR)
         if (i == 0) return
         check_dlgR_change = check_change(s, max_dlnR/ln10, &
            s% delta_lgR_limit, s% delta_lgR_hard_limit, &
            i, 'check_dlgR_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgR_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,3) 'lgR', i, s% lnR(i)/ln10
            write(*,3) 'lgR_old', i, s% lnR_start(i)/ln10
            write(*,3) 'dlgR', i, (s% lnR(i) - s% lnR_start(i))/ln10
            write(*,3) 'hard_limit', i, s% delta_lgR_hard_limit
         end if
      end function check_dlgR_change


      integer function check_lgL( &
            s, iso_in, msg, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         integer, intent(in) :: iso_in
         character (len=*), intent(in) :: msg
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio

         real(dp) :: &
            new_L, max_other_L, old_L, lim, hard_lim, lgL_min, &
            drop_factor, relative_limit, lgL, max_other_lgL, lgL_old, &
            abs_change, relative_excess
         logical, parameter :: dbg = .false.
         integer :: iso
         include 'formats'
         check_lgL = keep_going
         
         iso = iso_in
         if (iso == iprot) then ! check_lgL_power_photo_change
            if (s% log_max_temperature < s% min_lgT_for_lgL_power_photo_limit) return
            new_L = abs(s% power_photo)
            max_other_L = 0d0
            old_L = abs(s% power_photo_old)
            lim = s% delta_lgL_power_photo_limit
            hard_lim = s% delta_lgL_power_photo_hard_limit
            lgL_min = s% lgL_power_photo_burn_min
            drop_factor = s% lgL_power_photo_drop_factor
            relative_limit = 0d0
         else if (iso == ineut) then ! check_lgL_nuc_change
            if (s% log_max_temperature > s% max_lgT_for_lgL_nuc_limit) return
            new_L = s% power_nuc_burn
            max_other_L = 0d0
            old_L = s% power_nuc_burn_old
            if (s% log_max_temperature > &
                  s% delta_lgL_nuc_at_high_T_limit_lgT_min) then
               lim = s% delta_lgL_nuc_at_high_T_limit
               hard_lim = s% delta_lgL_nuc_at_high_T_hard_limit
            else
               lim = s% delta_lgL_nuc_limit
               hard_lim = s% delta_lgL_nuc_hard_limit
            end if
            lgL_min = s% lgL_nuc_burn_min
            drop_factor = s% lgL_nuc_drop_factor
            relative_limit = 0d0
         else if (iso == ih1) then
            new_L = s% power_h_burn
            max_other_L = max(s% power_he_burn, s% power_z_burn)
            old_L = s% power_h_burn_old
            lim = s% delta_lgL_H_limit
            hard_lim = s% delta_lgL_H_hard_limit
            lgL_min = s% lgL_H_burn_min
            drop_factor = s% lgL_H_drop_factor
            relative_limit = s% lgL_H_burn_relative_limit
         else if (iso == ihe4) then
            new_L = s% power_he_burn
            max_other_L = max(s% power_h_burn, s% power_z_burn)
            old_L = s% power_he_burn_old
            lim = s% delta_lgL_He_limit
            hard_lim = s% delta_lgL_He_hard_limit
            lgL_min = s% lgL_He_burn_min
            drop_factor = s% lgL_He_drop_factor
            relative_limit = s% lgL_He_burn_relative_limit
         else if (iso == isi28) then
            new_L = s% power_z_burn
            max_other_L = max(s% power_h_burn, s% power_he_burn)
            old_L = s% power_z_burn_old
            lim = s% delta_lgL_z_limit
            hard_lim = s% delta_lgL_z_hard_limit
            lgL_min = s% lgL_z_burn_min
            drop_factor = s% lgL_z_drop_factor
            relative_limit = s% lgL_z_burn_relative_limit
         else
            stop 'bad iso arg for check_lgL'
         end if
         
         if (old_L < 0d0) return
         
         lim = lim*s% time_delta_coeff
         hard_lim = hard_lim*s% time_delta_coeff

         if (new_L < old_L) then
            lim = lim*drop_factor
            hard_lim = hard_lim*drop_factor
         end if

         if (dbg) write(*,*)
         if (dbg) write(*,1) trim(msg) // ' new_L', new_L
         if (dbg) write(*,1) 'old_L', old_L
         if (new_L <= 0 .or. old_L <= 0) return

         lgL = safe_log10(new_L)
         if (dbg) write(*,1) 'lgL', lgL
         if (lgL < lgL_min) return

         if (max_other_L > 0) then
            max_other_lgL = safe_log10(max_other_L)
            if (dbg) write(*,1) 'max_other_lgL', max_other_lgL
            if (max_other_lgL - relative_limit > lgL) return
         end if

         lgL_old = safe_log10(old_L)
         if (dbg) write(*,1) 'lgL_old', lgL_old
         abs_change = abs(lgL - lgL_old)
         if (dbg) write(*,1) 'abs_change', abs_change
         if (dbg) write(*,1) 'hard_lim', hard_lim
         if (hard_lim > 0 .and. abs_change > hard_lim .and. (.not. skip_hard_limit)) then
            if (s% report_dt_hard_limit_retries) then
               write(*,1) trim(msg) // ' end', lgL
               write(*,1) trim(msg) // ' start', lgL_old
               write(*,1) trim(msg) // ' delta', lgL - lgL_old
               write(*,1) trim(msg) // ' hard_lim', hard_lim
            end if
            check_lgL = retry
            s% retry_message = 'lgL hard limit'
            return
         end if

         if (dbg) write(*,1) 'lim', lim
         if (lim <= 0) return

         relative_excess = (abs_change - lim) / lim
         if (dbg) write(*,1) 'relative_excess', relative_excess
         dt_limit_ratio = 1d0/pow(s% timestep_dt_factor,relative_excess)
         if (dt_limit_ratio <= 1d0) dt_limit_ratio = 0

      end function check_lgL


      integer function check_lgL_H_change(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         check_lgL_H_change = check_lgL( &
            s, ih1, 'check_lgL_H_change', skip_hard_limit, dt, dt_limit_ratio)
      end function check_lgL_H_change


      integer function check_lgL_He_change(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         check_lgL_He_change = check_lgL( &
            s, ihe4, 'check_lgL_He_change', skip_hard_limit, dt, dt_limit_ratio)
      end function check_lgL_He_change


      integer function check_lgL_z_change(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         check_lgL_z_change = check_lgL( &
            s, isi28, 'check_lgL_z_change', skip_hard_limit, dt, dt_limit_ratio)
      end function check_lgL_z_change


      integer function check_lgL_power_photo_change(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         check_lgL_power_photo_change = check_lgL( &
            s, iprot, 'check_lgL_power_photo_change', skip_hard_limit, dt, dt_limit_ratio)
      end function check_lgL_power_photo_change


      integer function check_lgL_nuc_change(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         include 'formats'
         check_lgL_nuc_change = check_lgL( &
            s, ineut, 'check_lgL_nuc_change', skip_hard_limit, dt, dt_limit_ratio)
      end function check_lgL_nuc_change


      integer function check_lgL_nuc_cat_change( &
            s, n_mix_bdy, mix_bdy_q, skip_hard_limit, dt_limit_ratio)
         use rates_def
         use num_lib, only: binary_search
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         integer, intent(in) :: n_mix_bdy
         real(dp), intent(in), pointer :: mix_bdy_q(:)
         real(dp), intent(inout) :: dt_limit_ratio

         integer :: k, max_j, max_k, bdy
         real(dp) :: max_lgL_diff, relative_excess, max_diff, cat_burn_min, &
            max_luminosity, max_luminosity_start

         include 'formats'

         check_lgL_nuc_cat_change = keep_going

         if (s% delta_lgL_nuc_cat_limit <= 0 .and. s% delta_lgL_nuc_cat_hard_limit <= 0) return

         cat_burn_min = exp10(s% lgL_nuc_cat_burn_min)
         max_diff = 0
         max_j = -1
         max_k = -1
         bdy = -1

         do k = 1, s% nz

            ! find the nearest mixing boundary
            bdy = binary_search(n_mix_bdy, mix_bdy_q, bdy, s% q(k))
            ! don't check cells near a mixing boundary
            if (bdy > 0) then
               if (abs(s% q(k) - mix_bdy_q(bdy)) < s% lgL_nuc_mix_dist_limit) cycle
            end if

            call do1_category(ipp,k)
            call do1_category(icno,k)
            call do1_category(i3alf,k)
            call do1_category(i_burn_c,k)
            call do1_category(i_burn_n,k)
            call do1_category(i_burn_o,k)
            call do1_category(i_burn_ne,k)
            call do1_category(i_burn_na,k)
            call do1_category(i_burn_mg,k)
            call do1_category(i_burn_si,k)
            call do1_category(i_burn_s,k)
            call do1_category(i_burn_ar,k)
            call do1_category(i_burn_ca,k)
            call do1_category(i_burn_ti,k)
            call do1_category(i_burn_cr,k)
            call do1_category(i_burn_fe,k)
            call do1_category(icc,k)
            call do1_category(ico,k)
            call do1_category(ioo,k)

         end do


         if (max_diff <= 0) return

         max_lgL_diff = log10(max_diff/Lsun)
         s% Tlim_dlgL_nuc_category = max_j
         s% Tlim_dlgL_nuc_cell = max_k

         check_lgL_nuc_cat_change = check_change(s, max_lgL_diff, &
            s% delta_lgL_nuc_cat_limit, s% delta_lgL_nuc_cat_hard_limit, &
            max_k, 'check_lgL_nuc_cat_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_lgL_nuc_cat_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'max_luminosity ' // trim(category_name(max_j)), max_luminosity
            write(*,1) 'max_luminosity_start ' // trim(category_name(max_j)), max_luminosity_start
         end if

         contains

         subroutine do1_category(j, k)
            integer, intent(in) :: j, k
            real(dp) :: diff, abs_diff
            if (s% luminosity_by_category(j,k) < cat_burn_min) return
            if (s% luminosity_by_category_start(j,k) < cat_burn_min) return
            diff = s% luminosity_by_category(j,k) - s% luminosity_by_category_start(j,k)
            abs_diff = abs(diff)
            if (abs_diff <= max_diff) return
            max_luminosity = s% luminosity_by_category(j,k)
            max_luminosity_start = s% luminosity_by_category_start(j,k)
            max_diff = abs_diff
            max_j = j
            max_k = k
         end subroutine do1_category

      end function check_lgL_nuc_cat_change


      integer function check_dlgTeff_change(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess
         include 'formats'
         check_dlgTeff_change = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax .or. s% Teff_old <= 0 .or. s% Teff <= 0) return
         check_dlgTeff_change = check_change(s, safe_log10(s% Teff/s% Teff_old), &
            s% delta_lgTeff_limit, s% delta_lgTeff_hard_limit, &
            0, 'check_dlgTeff_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgTeff_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lgTeff', safe_log10(s% Teff)
            write(*,1) 'lgTeff_old', safe_log10(s% Teff_old)
         end if
      end function check_dlgTeff_change


      integer function check_dYe_highT_change( &
            s, skip_hard_limit, dt_limit_ratio) ! check max change in Ye
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, max_diff, ye_diff, T_limit
         integer :: i, k
         include 'formats'
         check_dYe_highT_change = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         i = 0
         max_diff = 0
         T_limit = s% minT_for_highT_Ye_limit
         do k=1, s% nz
            if (s% T(k) < T_limit) cycle
            ye_diff = abs(s% ye(k) - s% ye_start(k))
            if (ye_diff <= max_diff) cycle
            max_diff = ye_diff
            i = k
         end do
         check_dYe_highT_change = check_change(s, max_diff, &
            s% delta_Ye_highT_limit, s% delta_Ye_highT_hard_limit, &
            i, 'check_dYe_highT_change', .false., dt_limit_ratio, relative_excess)
         if (check_dYe_highT_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,2) 'ye', i, s% ye(i)
            write(*,2) 'ye_start', i, s% ye_start(i)
         end if
      end function check_dYe_highT_change


      integer function check_dlgT_max_change(s, skip_hard_limit, dt, dt_limit_ratio)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, change, lnTmax, lnTmax_start
         integer :: lnTmax_k
         include 'formats'
         check_dlgT_max_change = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         if (s% delta_lgT_max_limit_lgT_min < 0d0) return
         if (s% delta_lgT_max_limit_only_after_near_zams) then
            if (s% X(s% nz) > 0.1d0 .and. &
                s% L_nuc_burn_total/s% L_phot < s% Lnuc_div_L_zams_limit ) return
         end if
         lnTmax_k = maxloc(s% lnT(1:s% nz),dim=1)
         lnTmax = s% lnT(lnTmax_k)
         if (lnTmax < s% delta_lgT_max_limit_lgT_min*ln10) return
         lnTmax_start = maxval(s% lnT_start(1:s% nz))
         change = (lnTmax - lnTmax_start)/ln10
         check_dlgT_max_change = check_change(s, change, &
            s% delta_lgT_max_limit, s% delta_lgT_max_hard_limit, &
            lnTmax_k, 'check_dlgT_max_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgT_max_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,2) 'lgT_max', lnTmax_k, lnTmax/ln10
            write(*,2) 'lgT_max_old', lnTmax_k, lnTmax_start/ln10
         end if
      end function check_dlgT_max_change


      integer function check_dlgT_max_at_high_T_change(s, skip_hard_limit, dt, dt_limit_ratio)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, change, lnTmax_at_high_T, lnTmax_at_high_T_start
         integer :: lnTmax_k
         include 'formats'
         check_dlgT_max_at_high_T_change = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         if (s% delta_lgT_max_at_high_T_limit_lgT_min < 0d0) return
         lnTmax_k = maxloc(s% lnT(1:s% nz),dim=1)
         lnTmax_at_high_T = s% lnT(lnTmax_k)
         if (lnTmax_at_high_T < s% delta_lgT_max_at_high_T_limit_lgT_min*ln10) return
         lnTmax_at_high_T_start = maxval(s% lnT_start(1:s% nz))
         change = (lnTmax_at_high_T - lnTmax_at_high_T_start)/ln10
         check_dlgT_max_at_high_T_change = check_change(s, change, &
            s% delta_lgT_max_at_high_T_limit, s% delta_lgT_max_at_high_T_hard_limit, &
            lnTmax_k, 'check_dlgT_max_at_high_T_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgT_max_at_high_T_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,2) 'lgT_max', lnTmax_k, lnTmax_at_high_T/ln10
            write(*,2) 'lgT_max_old', lnTmax_k, lnTmax_at_high_T_start/ln10
         end if
      end function check_dlgT_max_at_high_T_change


      integer function check_dlgT_cntr_change(s, skip_hard_limit, dt, dt_limit_ratio)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, change
         include 'formats'
         check_dlgT_cntr_change = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         if (s% delta_lgT_cntr_limit_only_after_near_zams) then
            if (s% X(s% nz) > 0.1d0 .and. &
                s% L_nuc_burn_total/s% L_phot < s% Lnuc_div_L_zams_limit) return
         end if
         change = (s% lnT(s% nz) - s% lnT_start(s% nz))/ln10
         check_dlgT_cntr_change = check_change(s, change, &
            s% delta_lgT_cntr_limit, s% delta_lgT_cntr_hard_limit, &
            s% nz, 'check_dlgT_cntr_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgT_cntr_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lgT_cntr', s% lnT(s% nz)/ln10
            write(*,1) 'lgT_cntr_old', s% lnT_start(s% nz)/ln10
         end if
      end function check_dlgT_cntr_change


      integer function check_dlgP_cntr_change(s, skip_hard_limit, dt, dt_limit_ratio)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, change
         include 'formats'
         check_dlgP_cntr_change = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         change = (s% lnP(s% nz) - s% lnP_start(s% nz))/ln10
         check_dlgP_cntr_change = check_change(s, change, &
            s% delta_lgP_cntr_limit, s% delta_lgP_cntr_hard_limit, &
            s% nz, 'check_dlgP_cntr_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgP_cntr_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lgP_cntr', s% lnP(s% nz)/ln10
            write(*,1) 'lgP_cntr_old', s% lnP_start(s% nz)/ln10
         end if
      end function check_dlgP_cntr_change


      integer function check_dlgRho_cntr_change(s, skip_hard_limit, dt, dt_limit_ratio)
         use const_def, only:ln10
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, dlgRho_cntr
         integer :: nz
         include 'formats'
         check_dlgRho_cntr_change = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         nz = s% nz
         dlgRho_cntr = (s% lnd(nz) - s% lnd_start(nz))/ln10
         check_dlgRho_cntr_change = check_change(s, dlgRho_cntr, &
            s% delta_lgRho_cntr_limit, s% delta_lgRho_cntr_hard_limit, &
            nz, 'check_dlgRho_cntr_change', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlgRho_cntr_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lgRho_cntr', s% lnd(s% nz)/ln10
            write(*,1) 'lgRho_cntr_old', s% lnd_start(s% nz)/ln10
         end if
      end function check_dlgRho_cntr_change


      integer function check_dlog_eps_nuc_change(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, max_ratio, ratio, delta, &
            limit_ratio, delta_r
         integer :: k, j, nz, k_max
         include 'formats'
         check_dlog_eps_nuc_change = keep_going
         nz = s% nz
         limit_ratio = exp10(s% delta_log_eps_nuc_limit)
         max_ratio = limit_ratio
         k_max = 0
         zoneloop: do k=1,nz
            if (s% eps_nuc_start(k) < 1) cycle zoneloop
            ratio = s% eps_nuc(k)/s% eps_nuc_start(k)
            if (s% mixing_type(k) /= convective_mixing .and. &
                s% mixing_type(min(nz,k+1)) /= convective_mixing .and. &
                ratio > max_ratio) then
               do j = 1, s% num_conv_boundaries
                  delta_r = abs(s% r(s% conv_bdy_loc(j)) - s% r(k))
                  if (delta_r <= s% scale_height(k)) then
                     cycle zoneloop ! skip ones that are too close to convection zone
                  end if
               end do
               max_ratio = ratio
               k_max = k
            end if
         end do zoneloop
         if (k_max > 0) then
            delta = log10(max_ratio)
         else
            delta = 0
         end if
         check_dlog_eps_nuc_change = check_change(s, delta, &
            s% delta_log_eps_nuc_limit, s% delta_log_eps_nuc_hard_limit, &
            k_max, 'check_dlog_eps_nuc_change', &
            skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dlog_eps_nuc_change /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,2) 'log_eps_nuc', k_max, safe_log10(abs(s% eps_nuc(k_max)))
            write(*,2) 'log_eps_nuc_old', k_max, safe_log10(abs(s% eps_nuc_start(k_max)))
         end if
      end function check_dlog_eps_nuc_change


      integer function check_dX_div_X_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, max_abs_dX_div_X, X, dX, abs_dX_div_X
         integer :: j, nz, j_max
         include 'formats'
         check_dX_div_X_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         nz = s% nz
         max_abs_dX_div_X = -1
         do j=1,s% species
            X = s% xa(j,nz)
            if (X > s% delta_dX_div_X_cntr_max) cycle
            if (X < s% delta_dX_div_X_cntr_min) cycle
            if (X <= 0d0) cycle
            dX = X - s% xa_old(j,nz)
            if (s% delta_dX_div_X_drop_only .and. dX > 0) cycle
            abs_dX_div_X = abs(dX/X)
            if (abs_dX_div_X > max_abs_dX_div_X) then
               max_abs_dX_div_X = abs_dX_div_X
               j_max = j
            end if
         end do
         if (max_abs_dX_div_X <= 0d0) return
         check_dX_div_X_cntr = check_change(s, max_abs_dX_div_X, &
            s% delta_dX_div_X_cntr_limit, s% delta_dX_div_X_cntr_hard_limit, &
            0, 'check_dX_div_X_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dX_div_X_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) chem_isos% name(s% chem_id(j_max)) // ' X', s% xa(j_max,nz)
            write(*,1) chem_isos% name(s% chem_id(j_max)) // ' X old', s% xa_old(j_max,nz)
         end if
      end function check_dX_div_X_cntr


      integer function check_lg_XH_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, lg_XH_cntr, lg_XH_cntr_old
         integer :: h1, nz
         include 'formats'
         check_lg_XH_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         h1 = s% net_iso(ih1)
         if (h1 == 0) return
         nz = s% nz
         if (s% xa(h1,nz) < 1d-10) return
         lg_XH_cntr = log10(s% xa(h1,nz))
         if (lg_XH_cntr > s% delta_lg_XH_cntr_max) return
         if (lg_XH_cntr < s% delta_lg_XH_cntr_min) return
         lg_XH_cntr_old = safe_log10(s% xa_old(h1,nz))
         if (s% delta_lg_XH_drop_only .and. lg_XH_cntr >= lg_XH_cntr_old) return
         check_lg_XH_cntr = check_change(s, lg_XH_cntr - lg_XH_cntr_old, &
            s% delta_lg_XH_cntr_limit, s% delta_lg_XH_cntr_hard_limit, &
            nz, 'check_lg_XH_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_lg_XH_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lg_XH_cntr', lg_XH_cntr
            write(*,1) 'lg_XH_cntr_old', lg_XH_cntr_old
            write(*,1) 'delta', lg_XH_cntr - lg_XH_cntr_old
         end if
      end function check_lg_XH_cntr


      integer function check_lg_XHe_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, lg_XHe_cntr, lg_XHe_cntr_old
         integer :: h1, he4, nz
         real(dp) :: xh1, xhe4
         include 'formats'
         check_lg_XHe_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         if (h1 == 0 .or. he4 == 0) return
         nz = s% nz
         xh1 = s% xa(h1,nz)
         xhe4 = s% xa(he4,nz)
         if (xhe4 < max(xh1, 1d-10)) return
         lg_XHe_cntr = log10(xhe4)
         if (lg_XHe_cntr > s% delta_lg_XHe_cntr_max) return
         if (lg_XHe_cntr < s% delta_lg_XHe_cntr_min) return
         lg_XHe_cntr_old = safe_log10(s% xa_old(he4,nz))
         if (s% delta_lg_XHe_drop_only .and. lg_XHe_cntr >= lg_XHe_cntr_old) return
         check_lg_XHe_cntr = check_change(s, lg_XHe_cntr - lg_XHe_cntr_old, &
            s% delta_lg_XHe_cntr_limit, s% delta_lg_XHe_cntr_hard_limit, &
            nz, 'check_lg_XHe_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_lg_XHe_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lg_XHe_cntr', lg_XHe_cntr
            write(*,1) 'lg_XHe_cntr_old', lg_XHe_cntr_old
            write(*,1) 'delta', lg_XHe_cntr - lg_XHe_cntr_old
         end if
      end function check_lg_XHe_cntr


      integer function check_lg_XC_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, lg_XC_cntr, lg_XC_cntr_old
         integer :: h1, he4, c12, nz
         real(dp) :: xh1, xhe4, xc12
         include 'formats'
         check_lg_XC_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         c12 = s% net_iso(ic12)
         if (h1 == 0 .or. he4 == 0 .or. c12 == 0) return
         nz = s% nz
         xh1 = s% xa(h1,nz)
         xhe4 = s% xa(he4,nz)
         xc12 = s% xa(c12,nz)
         if (xc12 < max(xh1, xhe4, 1d-10)) return
         if (s% xa(c12,nz) < 1d-10) return
         lg_XC_cntr = log10(xc12)
         if (lg_XC_cntr > s% delta_lg_XC_cntr_max) return
         if (lg_XC_cntr < s% delta_lg_XC_cntr_min) return
         lg_XC_cntr_old = safe_log10(s% xa_old(c12,nz))
         if (s% delta_lg_XC_drop_only .and. lg_XC_cntr >= lg_XC_cntr_old) return
         check_lg_XC_cntr = check_change(s, lg_XC_cntr - lg_XC_cntr_old, &
            s% delta_lg_XC_cntr_limit, s% delta_lg_XC_cntr_hard_limit, &
            nz, 'check_lg_XC_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_lg_XC_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lg_XC_cntr', lg_XC_cntr
            write(*,1) 'lg_XC_cntr_old', lg_XC_cntr_old
            write(*,1) 'delta', lg_XC_cntr - lg_XC_cntr_old
         end if
      end function check_lg_XC_cntr


      integer function check_lg_XNe_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, lg_XNe_cntr, lg_XNe_cntr_old
         integer :: h1, he4, c12, o16, nz
         real(dp) :: xh1, xhe4, xc12, XNe16
         include 'formats'
         check_lg_XNe_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         c12 = s% net_iso(ic12)
         o16 = s% net_iso(io16)
         if (h1 == 0 .or. he4 == 0 .or. c12 == 0 .or. o16 == 0) return
         nz = s% nz
         xh1 = s% xa(h1,nz)
         xhe4 = s% xa(he4,nz)
         xc12 = s% xa(c12,nz)
         XNe16 = s% xa(o16,nz)
         if (XNe16 < max(xh1, xhe4, xc12, 1d-10)) return
         lg_XNe_cntr = log10(XNe16)
         if (lg_XNe_cntr > s% delta_lg_XNe_cntr_max) return
         if (lg_XNe_cntr < s% delta_lg_XNe_cntr_min) return
         lg_XNe_cntr_old = safe_log10(s% xa_old(o16,nz))
         if (s% delta_lg_XNe_drop_only .and. lg_XNe_cntr >= lg_XNe_cntr_old) return
         check_lg_XNe_cntr = check_change(s, lg_XNe_cntr - lg_XNe_cntr_old, &
            s% delta_lg_XNe_cntr_limit, s% delta_lg_XNe_cntr_hard_limit, &
            nz, 'check_lg_XNe_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_lg_XNe_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lg_XNe_cntr', lg_XNe_cntr
            write(*,1) 'lg_XNe_cntr_old', lg_XNe_cntr_old
            write(*,1) 'delta', lg_XNe_cntr - lg_XNe_cntr_old
         end if
      end function check_lg_XNe_cntr


      integer function check_lg_XO_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, lg_XO_cntr, lg_XO_cntr_old
         integer :: h1, he4, c12, o16, nz
         real(dp) :: xh1, xhe4, xc12, xo16
         include 'formats'
         check_lg_XO_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         c12 = s% net_iso(ic12)
         o16 = s% net_iso(io16)
         if (h1 == 0 .or. he4 == 0 .or. c12 == 0 .or. o16 == 0) return
         nz = s% nz
         xh1 = s% xa(h1,nz)
         xhe4 = s% xa(he4,nz)
         xc12 = s% xa(c12,nz)
         xo16 = s% xa(o16,nz)
         if (xo16 < max(xh1, xhe4, xc12, 1d-10)) return
         lg_XO_cntr = log10(xo16)
         if (lg_XO_cntr > s% delta_lg_XO_cntr_max) return
         if (lg_XO_cntr < s% delta_lg_XO_cntr_min) return
         lg_XO_cntr_old = safe_log10(s% xa_old(o16,nz))
         if (s% delta_lg_XO_drop_only .and. lg_XO_cntr >= lg_XO_cntr_old) return
         check_lg_XO_cntr = check_change(s, lg_XO_cntr - lg_XO_cntr_old, &
            s% delta_lg_XO_cntr_limit, s% delta_lg_XO_cntr_hard_limit, &
            nz, 'check_lg_XO_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_lg_XO_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lg_XO_cntr', lg_XO_cntr
            write(*,1) 'lg_XO_cntr_old', lg_XO_cntr_old
            write(*,1) 'delta', lg_XO_cntr - lg_XO_cntr_old
         end if
      end function check_lg_XO_cntr


      integer function check_lg_XSi_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, lg_XSi_cntr, lg_XSi_cntr_old
         integer :: h1, he4, c12, o16, nz
         real(dp) :: xh1, xhe4, xc12, XSi16
         include 'formats'
         check_lg_XSi_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         c12 = s% net_iso(ic12)
         o16 = s% net_iso(io16)
         if (h1 == 0 .or. he4 == 0 .or. c12 == 0 .or. o16 == 0) return
         nz = s% nz
         xh1 = s% xa(h1,nz)
         xhe4 = s% xa(he4,nz)
         xc12 = s% xa(c12,nz)
         XSi16 = s% xa(o16,nz)
         if (XSi16 < max(xh1, xhe4, xc12, 1d-10)) return
         lg_XSi_cntr = log10(XSi16)
         if (lg_XSi_cntr > s% delta_lg_XSi_cntr_max) return
         if (lg_XSi_cntr < s% delta_lg_XSi_cntr_min) return
         lg_XSi_cntr_old = safe_log10(s% xa_old(o16,nz))
         if (s% delta_lg_XSi_drop_only .and. lg_XSi_cntr >= lg_XSi_cntr_old) return
         check_lg_XSi_cntr = check_change(s, lg_XSi_cntr - lg_XSi_cntr_old, &
            s% delta_lg_XSi_cntr_limit, s% delta_lg_XSi_cntr_hard_limit, &
            nz, 'check_lg_XSi_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_lg_XSi_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'lg_XSi_cntr', lg_XSi_cntr
            write(*,1) 'lg_XSi_cntr_old', lg_XSi_cntr_old
            write(*,1) 'delta', lg_XSi_cntr - lg_XSi_cntr_old
         end if
      end function check_lg_XSi_cntr


      integer function check_XH_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, XH_cntr, XH_cntr_old
         integer :: h1, nz
         include 'formats'
         check_XH_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         h1 = s% net_iso(ih1)
         if (h1 == 0) return
         nz = s% nz
         XH_cntr = s% xa(h1,nz)
         XH_cntr_old = s% xa_old(h1,nz)
         if (s% delta_XH_drop_only .and. XH_cntr >= XH_cntr_old) return
         check_XH_cntr = check_change(s, XH_cntr - XH_cntr_old, &
            s% delta_XH_cntr_limit, s% delta_XH_cntr_hard_limit, &
            nz, 'check_XH_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_XH_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'XH_cntr', XH_cntr
            write(*,1) 'XH_cntr_old', XH_cntr_old
            write(*,1) 'delta', XH_cntr - XH_cntr_old
         end if
      end function check_XH_cntr


      integer function check_XHe_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, XHe_cntr, XHe_cntr_old
         integer :: he4, nz
         include 'formats'
         check_XHe_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         he4 = s% net_iso(ihe4)
         if (he4 == 0) return
         nz = s% nz
         XHe_cntr = s% xa(he4,nz)
         XHe_cntr_old = s% xa_old(he4,nz)
         if (s% delta_XHe_drop_only .and. XHe_cntr >= XHe_cntr_old) return
         check_XHe_cntr = check_change(s, XHe_cntr - XHe_cntr_old, &
            s% delta_XHe_cntr_limit, s% delta_XHe_cntr_hard_limit, &
            nz, 'check_XHe_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_XHe_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'XHe_cntr', XHe_cntr
            write(*,1) 'XHe_cntr_old', XHe_cntr_old
            write(*,1) 'delta', XHe_cntr - XHe_cntr_old
         end if
      end function check_XHe_cntr


      integer function check_XC_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, XC_cntr, XC_cntr_old
         integer :: c12, nz
         include 'formats'
         check_XC_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         c12 = s% net_iso(ic12)
         if (c12 == 0) return
         nz = s% nz
         XC_cntr = s% xa(c12,nz)
         XC_cntr_old = s% xa_old(c12,nz)
         if (s% delta_XC_drop_only .and. XC_cntr >= XC_cntr_old) return
         check_XC_cntr = check_change(s, XC_cntr - XC_cntr_old, &
            s% delta_XC_cntr_limit, s% delta_XC_cntr_hard_limit, &
            nz, 'check_XC_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_XC_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'XC_cntr', XC_cntr
            write(*,1) 'XC_cntr_old', XC_cntr_old
            write(*,1) 'delta', XC_cntr - XC_cntr_old
         end if
      end function check_XC_cntr


      integer function check_XNe_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, XNe_cntr, XNe_cntr_old
         integer :: ne20, nz
         include 'formats'
         check_XNe_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         ne20 = s% net_iso(ine20)
         if (ne20 == 0) return
         nz = s% nz
         XNe_cntr = s% xa(ne20,nz)
         XNe_cntr_old = s% xa_old(ne20,nz)
         if (s% delta_XNe_drop_only .and. XNe_cntr >= XNe_cntr_old) return
         check_XNe_cntr = check_change(s, XNe_cntr - XNe_cntr_old, &
            s% delta_XNe_cntr_limit, s% delta_XNe_cntr_hard_limit, &
            nz, 'check_XNe_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_XNe_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'XNe_cntr', XNe_cntr
            write(*,1) 'XNe_cntr_old', XNe_cntr_old
            write(*,1) 'delta', XNe_cntr - XNe_cntr_old
         end if
      end function check_XNe_cntr


      integer function check_XO_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, XO_cntr, XO_cntr_old
         integer :: o16, nz
         include 'formats'
         check_XO_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         o16 = s% net_iso(io16)
         if (o16 == 0) return
         nz = s% nz
         XO_cntr = s% xa(o16,nz)
         XO_cntr_old = s% xa_old(o16,nz)
         if (s% delta_XO_drop_only .and. XO_cntr >= XO_cntr_old) return
         check_XO_cntr = check_change(s, XO_cntr - XO_cntr_old, &
            s% delta_XO_cntr_limit, s% delta_XO_cntr_hard_limit, &
            nz, 'check_XO_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_XO_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'XO_cntr', XO_cntr
            write(*,1) 'XO_cntr_old', XO_cntr_old
            write(*,1) 'delta', XO_cntr - XO_cntr_old
         end if
      end function check_XO_cntr


      integer function check_XSi_cntr(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess, XSi_cntr, XSi_cntr_old
         integer :: si28, nz
         include 'formats'
         check_XSi_cntr = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return
         si28 = s% net_iso(isi28)
         if (si28 == 0) return
         nz = s% nz
         XSi_cntr = s% xa(si28,nz)
         XSi_cntr_old = s% xa_old(si28,nz)
         if (s% delta_XSi_drop_only .and. XSi_cntr >= XSi_cntr_old) return
         check_XSi_cntr = check_change(s, XSi_cntr - XSi_cntr_old, &
            s% delta_XSi_cntr_limit, s% delta_XSi_cntr_hard_limit, &
            nz, 'check_XSi_cntr', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_XSi_cntr /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'XSi_cntr', XSi_cntr
            write(*,1) 'XSi_cntr_old', XSi_cntr_old
            write(*,1) 'delta', XSi_cntr - XSi_cntr_old
         end if
      end function check_XSi_cntr


      integer function check_delta_mdot(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: mdot, mdot_old, delta_mdot, lim, hard_lim
         check_delta_mdot = keep_going
         mdot = s% mstar_dot
         mdot_old = s% mstar_dot_old
         delta_mdot = abs(mdot - mdot_old)/ &
            (s% delta_mdot_atol*Msun/secyer + &
               s% delta_mdot_rtol*max(abs(mdot),abs(mdot_old)))
         if (delta_mdot == 0) return
         lim = s% delta_mdot_limit*s% time_delta_coeff
         hard_lim = s% delta_mdot_hard_limit*s% time_delta_coeff
         if (hard_lim > 0 .and. (.not. skip_hard_limit)) then
            if (delta_mdot > hard_lim) then
               if (s% report_dt_hard_limit_retries) &
                  write(*, '(a30, f20.10, 99e20.10)') 'delta_mdot_hard_limit', &
                     delta_mdot, hard_lim
               s% retry_message = 'delta_mdot_hard_limit'
               check_delta_mdot = retry
               return
            end if
         end if
         if (lim <= 0) return
         dt_limit_ratio = delta_mdot/lim
         if (dt_limit_ratio <= 1d0) dt_limit_ratio = 0
      end function check_delta_mdot


      integer function check_delta_mstar(s, skip_hard_limit, dt, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: delta_lg_star_mass, lim, hard_lim
         check_delta_mstar = keep_going
         delta_lg_star_mass = abs(log10(s% mstar/s% mstar_old))
         lim = s% delta_lg_star_mass_limit*s% time_delta_coeff
         hard_lim = s% delta_lg_star_mass_hard_limit*s% time_delta_coeff
         if (hard_lim > 0 .and. (.not. skip_hard_limit)) then
            if (delta_lg_star_mass > hard_lim) then
               if (s% report_dt_hard_limit_retries) &
                  write(*, '(a30, f20.10, 99e20.10)') 'delta_lg_star_mass_hard_limit', &
                     delta_lg_star_mass, hard_lim
               check_delta_mstar = retry
               s% retry_message = 'delta_lg_star_mass_hard_limit'
               return
            end if
         end if
         if (lim <= 0) return
         dt_limit_ratio = delta_lg_star_mass/lim
         if (dt_limit_ratio <= 1d0) dt_limit_ratio = 0
      end function check_delta_mstar


      integer function check_adjust_J_q(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: relative_excess
         include 'formats'
         check_adjust_J_q = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         if (.not. (s% rotation_flag .and. s% do_adjust_J_lost .and. s% mstar_dot < 0d0)) return
         ! we care about s% adjust_J_q remaining above a given limit
         ! so we use 1-S% adjust_J_q
         check_adjust_J_q = check_change(s, 1-s% adjust_J_q, &
            1-s% adjust_J_q_limit, &
            1-s% adjust_J_q_hard_limit, &
            0, 'check_adjust_J_q', &
            .false., dt_limit_ratio, relative_excess)
      end function check_adjust_J_q


      integer function check_delta_lgL(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: dlgL, relative_excess
         include 'formats'
         check_delta_lgL = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         if (s% L_surf < s% delta_lgL_limit_L_min .or. s% L_surf_old <= 0d0) then
            dlgL = 0
         else
            dlgL = log10(s% L_surf/s% L_surf_old)
         end if
         if (is_bad(dlgL)) then
            write(*,2) 's% L_surf', s% model_number, s% L_surf
            write(*,2) 's% L_surf_old', s% model_number, s% L_surf_old
            stop 'check_delta_lgL'
         end if
         check_delta_lgL = check_change(s, dlgL, &
            s% delta_lgL_limit, s% delta_lgL_hard_limit, &
            0, 'check_delta_lgL', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_delta_lgL /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'L_surf', s% L_surf
            write(*,1) 'L_surf_old', s% L_surf_old
         end if
      end function check_delta_lgL


      integer function check_delta_HR(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: dlgL, dlgTeff, dHR, relative_excess
         include 'formats'
         check_delta_HR = keep_going
         if (s% L_phot_old <= 0 .or. s% Teff_old <= 0 .or. &
             s% L_phot <= 0 .or. s% Teff <= 0) return
         dlgL = log10(s% L_phot/s% L_phot_old)
         dlgTeff = log10(s% Teff/s% Teff_old)
         dHR = sqrt(pow2(s% delta_HR_ds_L*dlgL) + pow2(s% delta_HR_ds_Teff*dlgTeff))
         if (is_bad(dHR)) then
            write(*,1) 's% L_phot_old', s% L_phot_old
            write(*,1) 's% L_phot', s% L_phot
            write(*,1) 's% Teff_old', s% Teff_old
            write(*,1) 's% Teff', s% Teff
            write(*,1) 'dHR', dHR
            if (s% stop_for_bad_nums) stop 'check_delta_HR'
         end if
         check_delta_HR = check_change(s, dHR, &
            s% delta_HR_limit, s% delta_HR_hard_limit, &
            0, 'check_delta_HR', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_delta_HR /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 's% L_phot_old', s% L_phot_old
            write(*,1) 's% L_phot', s% L_phot
            write(*,1) 's% Teff_old', s% Teff_old
            write(*,1) 's% Teff', s% Teff
            write(*,1) 'dHR', dHR
         end if
      end function check_delta_HR


      integer function check_rel_error_in_energy(s, skip_hard_limit, dt_limit_ratio)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: rel_error, relative_excess
         include 'formats'
         check_rel_error_in_energy = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         rel_error = abs(s% error_in_energy_conservation/s% total_energy_end)
         check_rel_error_in_energy = check_change(s, rel_error, &
            s% limit_for_rel_error_in_energy_conservation, &
            s% hard_limit_for_rel_error_in_energy_conservation, &
            0, 'check_rel_error_in_energy', &
            .false., dt_limit_ratio, relative_excess)
         if (check_rel_error_in_energy /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,1) 'error_in_energy_conservation', s% error_in_energy_conservation
            write(*,1) 'total_energy_end', s% total_energy_end
            write(*,1) 'rel_error', s% error_in_energy_conservation/s% total_energy_end
         end if
      end function check_rel_error_in_energy


      integer function check_dt_div_dt_cell_collapse(s, skip_hard_limit, dt, dt_limit_ratio)
         use star_utils, only: eval_min_cell_collapse_time
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: ratio, dt_timescale, relative_excess
         integer :: min_collapse_k, ierr
         include 'formats'
         check_dt_div_dt_cell_collapse = keep_going
         if (s% doing_relax) return
         dt_timescale = eval_min_cell_collapse_time( &
            s, 2, s% nz, min_collapse_k, ierr)
         if (ierr /= 0) return
         if (dt_timescale < 1d-30) return
         ratio = dt/dt_timescale
         !write(*,2) 'dt dt_cell_collapse ratio', min_collapse_k, dt, dt_timescale, ratio
         check_dt_div_dt_cell_collapse = check_change(s, ratio, &
            s% dt_div_dt_cell_collapse_limit, s% dt_div_dt_cell_collapse_hard_limit, &
            min_collapse_k, 'check_dt_div_dt_cell_collapse', skip_hard_limit, dt_limit_ratio, relative_excess)
         if (check_dt_div_dt_cell_collapse /= keep_going .and. s% report_dt_hard_limit_retries) then
            write(*,2) 'min dt_cell_collapse', min_collapse_k, dt_timescale
         end if
      end function check_dt_div_dt_cell_collapse


      integer function check_dt_div_min_dr_div_cs(s, skip_hard_limit, dt, dt_limit_ratio)
         use star_utils, only: min_dr_div_cs
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio
         real(dp) :: ratio, dt_x, relative_excess
         include 'formats'
         check_dt_div_min_dr_div_cs = keep_going
         dt_limit_ratio = 0d0
         if (s% doing_relax) return
         if (s% dt_div_min_dr_div_cs_limit <= 0d0) return
         dt_x = min_dr_div_cs(s, s% Tlim_dt_div_min_dr_div_cs_cell)
         !write(*,2) 'log min_dr_div_cs, q, m', s% Tlim_dt_div_min_dr_div_cs_cell, &
         !   safe_log10(dt_x), s% q(s% Tlim_dt_div_min_dr_div_cs_cell), s% m(s% Tlim_dt_div_min_dr_div_cs_cell)/Msun
         ratio = dt/dt_x
         check_dt_div_min_dr_div_cs = check_change(s, ratio, &
            s% dt_div_min_dr_div_cs_limit, s% dt_div_min_dr_div_cs_hard_limit, &
            s% Tlim_dt_div_min_dr_div_cs_cell, 'check_dt_div_min_dr_div_cs', &
            skip_hard_limit, dt_limit_ratio, relative_excess)
         if ((check_dt_div_min_dr_div_cs /= keep_going .and. s% report_dt_hard_limit_retries) .or. &
             (ratio > 1d0 .and. s% report_min_dr_div_cs)) then
            write(*,2) 'min_dr_div_cs', s% Tlim_dt_div_min_dr_div_cs_cell, dt_x
         end if
      end function check_dt_div_min_dr_div_cs


      integer function check_dX_nuc_drop(s, skip_hard_limit, dt, dt_limit_ratio)
         use rates_def, only: i_rate
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio

         integer :: k, nz, max_k, max_j
         real(dp), pointer, dimension(:) :: sig
         real(dp) :: max_dx_nuc_drop, X_limit, A_limit, min_dt, limit, hard_limit

         logical, parameter :: dbg = .false.

         include 'formats'

         check_dX_nuc_drop = keep_going
         if (s% mix_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) return

         X_limit = s% dX_nuc_drop_min_X_limit
         A_limit = s% dX_nuc_drop_max_A_limit
         nz = s% nz
         sig => s% sig

         max_dx_nuc_drop = 0
         max_k = 0
         max_j = 0
         do k = 1, nz
            call do1(k)
         end do

         s% Tlim_dXnuc_drop_cell = max_k
         s% Tlim_dXnuc_drop_species = max_j
         
         hard_limit = s% dX_nuc_drop_hard_limit*s% time_delta_coeff
         if (hard_limit > 0 .and. (.not. skip_hard_limit) .and. &
               max_dx_nuc_drop > hard_limit) then
            if (s% report_dt_hard_limit_retries) then
               write(*,2) trim(chem_isos% name(s% chem_id(max_j))), max_k, s% xa(max_j,max_k)
               write(*,2) trim(chem_isos% name(s% chem_id(max_j))) // ' old', max_k, s% xa_old(max_j,max_k)
               write(*,2) 'drop', max_k, max_dx_nuc_drop
            end if
            s% retry_message = 'dX_nuc_drop_hard_limit'
            s% retry_message_k = max_k
            check_dX_nuc_drop = retry
            return
         end if
         
         limit = s% dX_nuc_drop_limit*s% time_delta_coeff
         if (s% log_max_temperature >= 9.45d0 .and. s% dX_nuc_drop_limit_at_high_T > 0) &
            limit = s% dX_nuc_drop_limit_at_high_T
         if (limit <= 0 .or. max_dx_nuc_drop <= 0) return

         if (dt < secyer*s% dX_nuc_drop_min_yrs_for_dt) return
         min_dt = secyer*s% dX_nuc_drop_min_yrs_for_dt
         dt_limit_ratio = min( &
            max_dx_nuc_drop/limit, &
            1d0*dt/min_dt)
         if (dt_limit_ratio <= 1d0) dt_limit_ratio = 0

         s% dX_nuc_drop_max_j = max_j
         s% dX_nuc_drop_max_k = max_k
         s% dX_nuc_drop_max_drop = max_dx_nuc_drop
         
         contains

         subroutine do1(k)
            integer, intent(in) :: k

            integer :: j, jj, ii
            real(dp) :: dx, dx_drop, dm, dt_dm, dx_burning, dx_inflow, dxdt_nuc
            real(dp) ::dx00, dxp1, sig00, sigp1, flux00, fluxp1

            include 'formats'

            dm = s% dm(k)
            dt_dm = dt/dm

            do j=1, s% species

               if (chem_isos% W(s% chem_id(j)) > A_limit) then
                  if (dbg .and. k == 1387) &
                     write(*,2) 'dX_nuc_max_A_limit ' // trim(chem_isos% name(s% chem_id(j))), &
                        k, s% xa(j,k), chem_isos% W(s% chem_id(j)), A_limit
                  cycle
               end if

               if (chem_isos% Z(s% chem_id(j)) <= 2) then
                  cycle ! skip the little guys
               end if

               if (s% xa(j,k) < X_limit) then
                  if (dbg .and. k == 1387) &
                     write(*,2) &
                        'dX_nuc_drop_min_X_limit ' // trim(chem_isos% name(s% chem_id(j))), &
                        k, s% xa(j,k), X_limit
                  cycle
               end if

               dxdt_nuc = s% dxdt_nuc(j,k)
               if (dxdt_nuc >= 0) cycle

               dx_burning = dxdt_nuc*dt
               sig00 = sig(k)

               if (k < s% nz) then
                  sigp1 = sig(k+1)
               else
                  sigp1 = 0
               end if
               
               if (k > 1) then
                  dx00 = s% xa(j,k-1) - s% xa(j,k)
                  flux00 = -sig00*dx00
               else
                  flux00 = 0
               end if

               if (k < s% nz) then
                  dxp1 = s% xa(j,k) - s% xa(j,k+1)
                  fluxp1 = -sigp1*dxp1
               else
                  fluxp1 = 0
               end if

               dx_inflow = max(0d0, fluxp1, -flux00)*dt_dm
               
               dx_drop = -(dx_burning + dx_inflow) ! dx_burning < 0 for drop

               dx = s% xa_old(j,k) - s% xa(j,k) ! the actual drop
               if (dx < dx_drop) dx_drop = dx

               if (dx_drop > max_dx_nuc_drop) then
                  max_dx_nuc_drop = dx_drop
                  max_k = k
                  max_j = j
               end if

            end do

         end subroutine do1


      end function check_dX_nuc_drop


      integer function check_varcontrol_limit(s, dt_limit_ratio)
         type (star_info), pointer :: s
         real(dp), intent(inout) :: dt_limit_ratio

         real(dp) :: varcontrol, vc_target
         integer :: ierr

         include 'formats'

         ierr = 0
         check_varcontrol_limit = keep_going

         varcontrol = eval_varcontrol(s, ierr)
         if (ierr /= 0) then
            check_varcontrol_limit = retry
            s% retry_message = 'varcontrol hard limit'
            if (s% report_ierr) write(*, *) 'check_varcontrol_limit: eval_varcontrol ierr', ierr
            s% result_reason = nonzero_ierr
            return
         end if
         
         if (s% varcontrol_target < s% min_allowed_varcontrol_target) then
            check_varcontrol_limit = terminate
            write(*, *) 'ERROR: timestep varcontrol_target < min_allowed_varcontrol_target'
            s% result_reason = nonzero_ierr
            return
         end if

         vc_target = max(1d-99,s% varcontrol_target)
         dt_limit_ratio = varcontrol/vc_target

         if (s% report_solver_dt_info) then
            write(*, 1) 's% varcontrol_target', s% varcontrol_target
            write(*, 1) 'vc_target', vc_target
            write(*, 1) 'varcontrol', varcontrol
            write(*, 1) 'dt_limit_ratio', dt_limit_ratio
            write(*, *)
         end if

         if (dt_limit_ratio > s% varcontrol_dt_limit_ratio_hard_max) then
            write(*, '(a50, f20.10, 99e20.10)') 'varcontrol_dt_limit_ratio too large', &
               dt_limit_ratio, varcontrol, vc_target
            check_varcontrol_limit = retry
            s% retry_message = 'varcontrol_dt_limit_ratio_hard_max'
            return
         end if

      end function check_varcontrol_limit


      real(dp) function eval_varcontrol(s, ierr) result(varcontrol)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: j, nterms, nvar_hydro, nz, k, kk, iounit, &
            skip1, skip2, skip3, skip4, skip5, skip6, i_alpha_RTI, i_etrb_RSP, i_etrb
         real(dp) :: sumj, sumvar, sumscales, sumterm(s% nvar_total)
         real(dp), pointer :: vc_data(:,:)
         logical :: dbg
         real(dp), parameter :: xscale_min = 1

         vc_data => null()

         include 'formats'

         ierr = 0

         varcontrol = 1d99
         dbg = .false.
         nvar_hydro = s% nvar_hydro
         nz = s% nz

         if (dbg) then
            allocate(vc_data(nvar_hydro,nz),stat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in alloca in eval_varcontrol'
               return
            end if
         end if
         
         if (s% include_L_in_error_est) then
            skip1 = 0
         else
            skip1 = s% i_lum
         end if

         if (s% include_v_in_error_est) then
            skip2 = 0
         else
            skip2 = s% i_v
         end if

         if (s% include_u_in_error_est) then
            skip3 = 0
         else
            skip3 = s% i_u
         end if

         if (s% solver_use_lnR) then
            skip4 = 0
         else
            skip4 = s% i_lnR
         end if

         if (s% solver_use_lnT) then
            skip5 = 0
         else
            skip5 = s% i_lnT
         end if

         if (s% solver_use_lnd) then
            skip6 = 0
         else
            skip6 = s% i_lnd
         end if

         i_alpha_RTI = s% i_alpha_RTI
         i_etrb_RSP = s% i_etrb_RSP
         i_etrb = s% i_etrb

         nterms = 0
         sumvar = 0
         sumscales = 0
         sumterm(:) = 0

         if (.not. associated(s% xh_old)) then
            stop 'not associated xh_old'
         end if

         ! use differences in smoothed old and new to filter out high frequency noise.
         do j = 1, nvar_hydro

            if (j == skip1 .or. &
                j == skip2 .or. &
                j == skip3 .or. &
                j == skip4 .or. &
                j == skip5 .or. &
                j == skip6 .or. &
                j == s% i_ln_cvpv0 .or. &
                j == s% i_j_rot .or. &
                j == s% i_w_div_wc .or. & ! TODO: check why not including this makes restart varcontrol inconsistent
                j == i_alpha_RTI .or. &
                j == i_etrb .or. &
                j == i_etrb_RSP) cycle

            nterms = nterms + nz
            do k = 3, nz-2
               sumj = abs(sum(s% xh(j,k-2:k+2)) - sum(s% xh_old(j,k-2:k+2)))/5
               sumterm(j) = sumterm(j) + sumj
            end do
            sumterm(j) = sumterm(j) + &
               abs((2*s% xh(j,1) + s% xh(j,2)) - (2*s% xh_old(j,1) + s% xh_old(j,2)))/3 + &
               abs((2*s% xh(j,nz) + s% xh(j,nz-1)) - (2*s% xh_old(j,nz) + s% xh_old(j,nz-1)))/3
            k = 2
            sumj = abs(sum(s% xh(j,k-1:k+1)) - sum(s% xh_old(j,k-1:k+1)))/3
            sumterm(j) = sumterm(j) + sumj
            k = nz-1
            sumj = abs(sum(s% xh(j,k-1:k+1)) - sum(s% xh_old(j,k-1:k+1)))/3

            if (j == s% i_lnd) then
               sumterm(j) = sumterm(j)/3 ! Seems to help. from Eggleton.
            end if

            sumvar = sumvar + sumterm(j)
            sumscales = sumscales + max(xscale_min, abs(s% xh_old(j,1)))

         end do

         sumterm(:) = sumterm(:)/sumscales
         sumvar = sumvar/sumscales

         if (dbg) then
            call show_info
            deallocate(vc_data)
            !stop 'debug: timestep'
         end if

         varcontrol = sumvar/nterms

         contains

         subroutine show_info
            character (len=64) :: filename
            real(dp) :: newterm
            write(*,*)
            write(*, *) 'sumvar', sumvar
            write(*, *) 'nterms', nterms
            write(*, *) 'sumvar/nterms', sumvar/nterms
            write(*,*)
            do j=1, nvar_hydro
               if (j == skip1 .or. &
                   j == skip2) cycle
               write(*, '(a40, d26.16)') &
                  'varcontrol fraction for ' // trim(s% nameofvar(j)) // ' = ', sumterm(j)/sumvar
            end do
            write(*,*)
            do j=1, s% species
               if (sumterm(nvar_hydro + j) /= 0) &
                  write(*, '(a40, d26.16)') 'varcontrol fraction for ' // &
                        trim(chem_isos% name(s% chem_id(j))) // ' = ', &
                        sumterm(nvar_hydro + j)/sumvar
            end do
         end subroutine show_info


      end function eval_varcontrol


      subroutine filter_dt_next(s, order, dt_limit_ratio_in)
         ! H211b "low pass" controller.
         ! Soderlind & Wang, J of Computational and Applied Math 185 (2006) 225 – 243.
         use num_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: order, dt_limit_ratio_in

         real(dp) :: ratio, ratio_prev, limtr, dt_limit_ratio_target, &
            dt_limit_ratio, beta1, beta2, alpha2

         include 'formats'

         beta1 = 0.25d0/order
         beta2 = 0.25d0/order
         alpha2 = 0.25d0
         
         dt_limit_ratio = max(1d-10, dt_limit_ratio_in)
         s% dt_limit_ratio = dt_limit_ratio
         dt_limit_ratio_target = 1d0

         if (s% use_dt_low_pass_controller .and. &
               s% dt_limit_ratio_old > 0 .and. s% dt_old > 0) then ! use 2 values to do "low pass" controller
            ratio = limiter(dt_limit_ratio_target/dt_limit_ratio)
            ratio_prev = limiter(dt_limit_ratio_target/s% dt_limit_ratio_old)
            limtr = limiter( &
               pow(ratio,beta1) * pow(ratio_prev,beta2) * pow(s% dt/s% dt_old,-alpha2))
            s% dt_next = s% dt*limtr

            if (s% report_solver_dt_info) then
               write(*,2) 'dt_limit_ratio_target', s% model_number, dt_limit_ratio_target
               write(*,2) 'dt_limit_ratio', s% model_number, dt_limit_ratio
               write(*,2) 's% dt_limit_ratio_old', s% model_number, s% dt_limit_ratio_old
               write(*,2) 'order', s% model_number, order
               write(*,2) 'ratio', s% model_number, ratio
               write(*,2) 'ratio_prev', s% model_number, ratio_prev
               write(*,2) 'limtr', s% model_number, limtr
               write(*,2) 's% dt_next', s% model_number, s% dt_next
               write(*,*)
            end if

         else ! no history available, so fall back to the 1st order controller
            s% dt_next = s% dt*dt_limit_ratio_target/dt_limit_ratio
         end if


         contains


         real(dp) function limiter(x)
            real(dp), intent(in) :: x
            real(dp), parameter :: kappa = 10 ! 2
            ! for x >= 0 and kappa = 2, limiter value is between 0.07 and 4.14
            ! for x >= 0 and kappa = 10, limiter value is between 0.003 and 16.7
            ! for x = 1, limiter = 1
            limiter = 1 + kappa*atan((x-1)/kappa)
         end function limiter


      end subroutine filter_dt_next


      end module timestep


