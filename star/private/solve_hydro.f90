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

      module solve_hydro

      use star_private_def
      use const_def

      implicit none

      private
      public :: set_L_burn_by_category, &
         set_surf_info, set_tol_correction, do_hydro_converge

      integer, parameter :: stencil_neighbors = 1
            ! number of neighbors on each side (e.g., =1 for 3 point stencil)


      contains


      integer function do_hydro_converge( &
            s, itermin, nvar, skip_global_corr_coeff_limit, &
            tol_correction_norm, tol_max_correction, dt_stage)
         ! return keep_going, retry, or terminate
         use mtx_lib
         use mtx_def
         use num_def
         use star_utils, only: start_time, update_time

         type (star_info), pointer :: s
         integer, intent(in) :: itermin, nvar
         logical, intent(in) :: skip_global_corr_coeff_limit
         real(dp), intent(in) :: tol_correction_norm, tol_max_correction, dt_stage

         integer :: ierr, nz, k, mljac, mujac, n, nzmax, &
            hydro_lwork, hydro_liwork, num_jacobians
         logical :: report

         include 'formats'

         if (dt_stage <= 0d0) then
            do_hydro_converge = keep_going
            return
         end if

         do_hydro_converge = terminate

         ierr = 0

         nz = s% nz
         n = nz*nvar
         mljac = (stencil_neighbors+1)*nvar-1
         mujac = mljac
         nzmax = 0

         call work_sizes_for_solver(ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'do_hydro_converge: work_sizes_for_solver failed'
            do_hydro_converge = retry
            s% result_reason = nonzero_ierr
            s% termination_code = t_solve_hydro
            return
         end if

         call alloc_for_solver(ierr)
         if (ierr /= 0) then
            s% termination_code = t_solve_hydro
            s% result_reason = nonzero_ierr
            return
         end if

         report = (s% report_solver_progress .or. s% report_ierr)
         
         s% solver_call_number = s% solver_call_number + 1

         do_hydro_converge = do_hydro_solver( &
            s, itermin, skip_global_corr_coeff_limit, &
            tol_correction_norm, tol_max_correction, dt_stage, &
            report, nz, nvar, s% hydro_work, hydro_lwork, &
            s% hydro_iwork, hydro_liwork)


         contains


         subroutine alloc_for_solver(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            if (.not. associated(s% hydro_iwork)) then
               allocate(s% hydro_iwork(hydro_liwork))
            else if (size(s% hydro_iwork, dim=1) < hydro_liwork) then
               deallocate(s% hydro_iwork)
               allocate(s% hydro_iwork(int(1.3d0*hydro_liwork)+100))
            end if
            if (.not. associated(s% hydro_work)) then
               allocate(s% hydro_work(hydro_lwork))
            else if (size(s% hydro_work, dim=1) < hydro_lwork) then
               deallocate(s% hydro_work)
               allocate(s% hydro_work(int(1.3d0*hydro_lwork)+100))
            end if
         end subroutine alloc_for_solver


         subroutine work_sizes_for_solver(ierr)
            use star_solver, only: get_solver_work_sizes
            use star_solver_15066, only: get_solver_work_sizes_15066
            integer, intent(out) :: ierr
            if (.not. s% conv_vel_flag) then
               call get_solver_work_sizes(s, nvar, nz, hydro_lwork, hydro_liwork, ierr)
            else
               call get_solver_work_sizes_15066(s, nvar, nz, hydro_lwork, hydro_liwork, ierr)
            end if
         end subroutine work_sizes_for_solver


      end function do_hydro_converge


      subroutine set_surf_info(s, nvar) ! set to values at start of step
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         if (s% i_lnd > 0 .and. s% i_lnd <= nvar) s% surf_lnd = s% xh(s% i_lnd,1)
         if (s% i_lnT > 0 .and. s% i_lnT <= nvar) s% surf_lnT = s% xh(s% i_lnT,1)
         if (s% i_lnR > 0 .and. s% i_lnR <= nvar) s% surf_lnR = s% xh(s% i_lnR,1)
         if (s% i_v > 0 .and. s% i_v <= nvar) s% surf_v = s% xh(s% i_v,1)
         s% surf_lnS = s% lnS(1)
         s% num_surf_revisions = 0
      end subroutine set_surf_info


      subroutine set_xh(s,nvar) ! set xh using current structure info
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         integer :: j1, k, nz
         include 'formats'
         nz = s% nz
         do j1 = 1, min(nvar,s% nvar_hydro)
            if (j1 == s% i_lnd .and. s% i_lnd <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% lnd(k)
               end do
            else if (j1 == s% i_lnT .and. s% i_lnT <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% lnT(k)
               end do
            else if (j1 == s% i_lnR .and. s% i_lnR <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% lnR(k)
               end do
            else if (j1 == s% i_lum .and. s% i_lum <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% L(k)
               end do
            else if (j1 == s% i_eturb .and. s% i_eturb <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% Eturb(k)
               end do
            else if (j1 == s% i_v .and. s% i_v <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% v(k)
               end do
            else if (j1 == s% i_u .and. s% i_u <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% u(k)
               end do
            else if (j1 == s% i_alpha_RTI .and. s% i_alpha_RTI <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% alpha_RTI(k)
               end do
            else if (j1 == s% i_ln_cvpv0 .and. s% i_ln_cvpv0 <= nvar) then
               do k = 1, nz
                  ! create a rough first guess using mlt_vc_start and conv_vel when
                  ! mlt_vc is larger than the starting conv_vel
                  if (s% mlt_vc_start(k) > 0d0 .and. s% mlt_vc_start(k) > s% conv_vel(k)) then
                     s% conv_vel(k) = s% conv_vel(k) + &
                        (s% mlt_vc_start(k) -s% conv_vel(k)) * &
                        min(1d0, s% dt*s% mlt_vc_start(k)/(s% scale_height_start(k)*s% mixing_length_alpha))
                  end if
                  s% xh(j1,k) = log(s% conv_vel(k)+s% conv_vel_v0)
               end do
            end if
         end do
      end subroutine set_xh


      subroutine set_tol_correction( &
            s, T_max, tol_correction_norm, tol_max_correction)
         type (star_info), pointer :: s
         real(dp), intent(in) :: T_max
         real(dp), intent(out) :: tol_correction_norm, tol_max_correction
         include 'formats'
         if (T_max >= s% tol_correction_extreme_T_limit) then
            tol_correction_norm = s% tol_correction_norm_extreme_T
            tol_max_correction = s% tol_max_correction_extreme_T
         else if (T_max >= s% tol_correction_high_T_limit) then
            tol_correction_norm = s% tol_correction_norm_high_T
            tol_max_correction = s% tol_max_correction_high_T
         else
            tol_correction_norm = s% tol_correction_norm
            tol_max_correction = s% tol_max_correction
         end if
      end subroutine set_tol_correction


      integer function do_hydro_solver( &
            s, itermin, skip_global_corr_coeff_limit, &
            tol_correction_norm, tol_max_correction, dt_stage, &
            report, nz, nvar, solver_work, solver_lwork, &
            solver_iwork, solver_liwork)
         ! return keep_going, retry, or terminate

         ! when using solver for hydro step,
         ! do not require that functions have been evaluated for starting configuration.
         ! when finish, will have functions evaluated for the final set of primary variables.
         ! for example, the reaction rates will have been computed, so they can be used
         ! as initial values in the following burn and mix.

         use num_def
         use alloc

         type (star_info), pointer :: s
         integer, intent(in) :: itermin, nz, nvar
         logical, intent(in) :: skip_global_corr_coeff_limit, report
         real(dp), intent(in) :: tol_correction_norm, tol_max_correction, dt_stage

         real(dp), pointer :: dx(:,:), dx1(:) ! dx => dx1
         integer, intent(in) :: solver_lwork, solver_liwork
         real(dp), pointer :: solver_work(:) ! (solver_lwork)
         integer, pointer :: solver_iwork(:) ! (solver_liwork)
         logical :: converged
         integer :: i, k, species, ierr, alph, j1, j2, gold_tolerances_level
         real(dp) :: varscale, r003, rp13, dV, frac, maxT

         real(dp), parameter :: xscale_min = 1d-3

         include 'formats'

         species = s% species
         do_hydro_solver = keep_going
         s% using_gold_tolerances = .false.
         gold_tolerances_level = 0
         
         if ((s% use_gold2_tolerances .and. s% steps_before_use_gold2_tolerances < 0) .or. &
             (s% steps_before_use_gold2_tolerances >= 0 .and. &
                s% model_number > s% steps_before_use_gold2_tolerances + s% init_model_number)) then
            s% using_gold_tolerances = .true.
            gold_tolerances_level = 2
         else if ((s% use_gold_tolerances .and. s% steps_before_use_gold_tolerances < 0) .or. &
             (s% steps_before_use_gold_tolerances >= 0 .and. &
                s% model_number > s% steps_before_use_gold_tolerances + s% init_model_number)) then
            if (s% maxT_for_gold_tolerances > 0) then
               maxT = maxval(s% T(1:nz))
            else
               maxT = -1d0
            end if
            if (maxT > s% maxT_for_gold_tolerances) then 
               !write(*,2) 'exceed maxT_for_gold_tolerances', &
               !   s% model_number, maxT, s% maxT_for_gold_tolerances
            else ! okay for maxT, so check if also ok for eosPC_frac
               s% using_gold_tolerances = .true.
               gold_tolerances_level = 1
            end if
         end if

         ! parameters for solver
         solver_iwork(1:num_iwork_params) = 0
         solver_work(1:num_work_params) = 0

         if (s% doing_first_model_of_run) &
            solver_iwork(i_try_really_hard) = 1 ! try_really_hard for 1st model
         solver_iwork(i_itermin) = itermin

         solver_iwork(i_max_iterations_for_jacobian) = 1

         if (s% doing_first_model_of_run) then
            solver_iwork(i_max_tries) = s% max_tries1
         else if (s% retry_cnt > 20) then
            solver_iwork(i_max_tries) = s% max_tries_after_20_retries
         else if (s% retry_cnt > 10) then
            solver_iwork(i_max_tries) = s% max_tries_after_10_retries
         else if (s% retry_cnt > 5) then
            solver_iwork(i_max_tries) = s% max_tries_after_5_retries
         else if (s% retry_cnt > 0) then
            solver_iwork(i_max_tries) = s% max_tries_for_retry
         else
            solver_iwork(i_max_tries) = s% solver_max_tries_before_reject
         end if

         solver_iwork(i_tiny_min_corr_coeff) = s% tiny_corr_coeff_limit
         if (s% report_solver_progress) then
            solver_iwork(i_debug) = 1
         else
            solver_iwork(i_debug) = 0
         end if
         solver_iwork(i_model_number) = s% model_number

         solver_work(r_tol_abs_slope_min) = -1 ! unused
         solver_work(r_tol_corr_resid_product) = -1 ! unused

         solver_work(r_scale_correction_norm) = s% scale_correction_norm
         solver_work(r_corr_param_factor) = s% corr_param_factor
         solver_work(r_scale_max_correction) = s% scale_max_correction
         solver_work(r_corr_norm_jump_limit) = s% corr_norm_jump_limit
         solver_work(r_max_corr_jump_limit) = s% max_corr_jump_limit
         solver_work(r_resid_norm_jump_limit) = s% resid_norm_jump_limit
         solver_work(r_max_resid_jump_limit) = s% max_resid_jump_limit
         solver_work(r_tiny_corr_factor) = s% tiny_corr_factor
         solver_work(r_dt) = dt_stage
         solver_work(r_tol_max_correction) = tol_max_correction
         
         if (gold_tolerances_level == 2) then
            solver_work(r_tol_residual_norm) = s% gold2_tol_residual_norm1
            solver_work(r_tol_max_residual) = s% gold2_tol_max_residual1
            solver_work(r_tol_residual_norm2) = s% gold2_tol_residual_norm2
            solver_work(r_tol_max_residual2) = s% gold2_tol_max_residual2
            solver_work(r_tol_residual_norm3) = s% gold2_tol_residual_norm3
            solver_work(r_tol_max_residual3) = s% gold2_tol_max_residual3
         else if (gold_tolerances_level == 1) then
            solver_work(r_tol_residual_norm) = s% gold_tol_residual_norm1
            solver_work(r_tol_max_residual) = s% gold_tol_max_residual1
            solver_work(r_tol_residual_norm2) = s% gold_tol_residual_norm2
            solver_work(r_tol_max_residual2) = s% gold_tol_max_residual2
            solver_work(r_tol_residual_norm3) = s% gold_tol_residual_norm3
            solver_work(r_tol_max_residual3) = s% gold_tol_max_residual3
         else
            solver_work(r_tol_residual_norm) = s% tol_residual_norm1
            solver_work(r_tol_max_residual) = s% tol_max_residual1
            solver_work(r_tol_residual_norm2) = s% tol_residual_norm2
            solver_work(r_tol_max_residual2) = s% tol_max_residual2
            solver_work(r_tol_residual_norm3) = s% tol_residual_norm3
            solver_work(r_tol_max_residual3) = s% tol_max_residual3
         end if
         
         if (skip_global_corr_coeff_limit) then
            solver_work(r_min_corr_coeff) = 1
         else if (s% conv_vel_flag) then
            solver_work(r_min_corr_coeff) = s% conv_vel_corr_coeff_limit
         else
            solver_work(r_min_corr_coeff) = s% corr_coeff_limit
         end if

         call non_crit_get_work_array(s, dx1, nvar*nz, nvar*nz_alloc_extra, 'solver', ierr)
         if (ierr /= 0) return
         dx(1:nvar,1:nz) => dx1(1:nvar*nz)
         s% solver_dx(1:nvar,1:nz) => dx1(1:nvar*nz)

         call set_xh(s, nvar) ! set xh using current structure info

         do k = 1, nz
            do j1 = 1, min(nvar, s% nvar_hydro)
               dx(j1,k) = s% xh(j1,k) - s% xh_start(j1,k)
            end do
         end do

         if (nvar >= s% i_chem1) then
            do k = 1, nz
               j2 = 1
               do j1 = s% i_chem1, nvar
                  s% xa_sub_xa_start(j2,k) = s% xa(j2,k) - s% xa_start(j2,k)
                  dx(j1,k) = s% xa_sub_xa_start(j2,k)
                  j2 = j2+1
               end do
            end do
         end if

         converged = .false.
         call hydro_solver_step( &
            s, nz, s% nvar_hydro, nvar, dx1, dt_stage, &
            gold_tolerances_level, tol_correction_norm, &
            solver_work, solver_lwork, &
            solver_iwork, solver_liwork, &
            converged, ierr)

         call non_crit_return_work_array(s, dx1, 'solver')
         nullify(s% solver_dx)

         if (ierr /= 0) then
            if (report) then
               write(*, *) 'hydro_solver_step returned ierr', ierr
               write(*, *) 's% model_number', s% model_number
               write(*, *) 'nz', nz
               write(*, *) 's% num_retries', s% num_retries
               write(*, *)
            end if
            do_hydro_solver = retry
            s% result_reason = nonzero_ierr
            s% dt_why_retry_count(Tlim_solver) = &
               s% dt_why_retry_count(Tlim_solver) + 1
            return
         end if

         s% num_solver_iterations = solver_iwork(i_num_jacobians)

         if (converged) then ! sanity checks before accept it
            converged = check_after_converge(s, report, ierr)
            if (converged .and. s% RTI_flag) & ! special checks
               converged = RTI_check_after_converge(s, report, ierr)
         end if

         if (.not. converged) then
            do_hydro_solver = retry
            s% result_reason = hydro_failed_to_converge
            s% dt_why_retry_count(Tlim_solver) = &
               s% dt_why_retry_count(Tlim_solver) + 1
            if (report) then
               write(*,2) 'solver rejected trial model'
               write(*,2) 's% model_number', s% model_number
               write(*,2) 's% solver_call_number', s% solver_call_number
               write(*,2) 'nz', nz
               write(*,2) 's% num_retries', s% num_retries
               write(*,1) 'log dt/secyer', log10(dt_stage/secyer)
               write(*, *)
            end if
            return
         end if

      end function do_hydro_solver


      logical function RTI_check_after_converge(s, report, ierr) result(converged)
         use mesh_adjust, only: set_lnT_for_energy
         use micro, only: do_eos_for_cell
         use chem_def, only: ih1, ihe3, ihe4
         type (star_info), pointer :: s
         logical, intent(in) :: report
         integer, intent(out) :: ierr
         integer :: k, nz
         real(dp) :: old_energy, old_IE, new_IE, old_KE, new_KE, new_u, &
            revised_energy, new_lnT
         include 'formats'
         ierr = 0
         nz = s% nz
         converged = .true.
         !return
         do k=1,nz
            if (k < nz .and. s% alpha_RTI(k) < 1d-10) cycle
            old_energy = s% energy(k)
            old_IE = old_energy*s% dm(k)
            if (s% energy(k) < s% RTI_energy_floor) then
               ! try to take from KE to give to IE
               ! else just bump energy and take hit to energy conservation
               s% energy(k) = s% RTI_energy_floor
               s% lnE(k) = log(s% energy(k))
               call set_lnT_for_energy(s, k, &
                  s% net_iso(ih1), s% net_iso(ihe3), s% net_iso(ihe4), &
                  s% species, s% xa(:,k), &
                  s% rho(k), s% lnd(k)/ln10, s% energy(k), s% lnT(k), &
                  new_lnT, revised_energy, ierr)
               if (ierr /= 0) return ! stop 'do_merge failed in set_lnT_for_energy'
               s% xh(s% i_lnT,k) = new_lnT
               s% lnT(k) = new_lnT
               s% T(k) = exp(new_lnT)
            end if
            new_IE = s% energy(k)*s% dm(k)
            old_KE = 0.5d0*s% dm(k)*s% u(k)*s% u(k)
            new_KE = max(0d0, old_KE + old_IE - new_IE)
            new_u = sqrt(new_KE/(0.5d0*s% dm(k)))
            if (s% u(k) > 0d0) then
               s% u(k) = new_u
            else
               s% u(k) = -new_u
            end if
            s% xh(s% i_u, k) = s% u(k)
         end do
      end function RTI_check_after_converge


      logical function check_after_converge(s, report, ierr) result(converged)
         type (star_info), pointer :: s
         logical, intent(in) :: report
         integer, intent(out) :: ierr
         integer :: k, nz
         include 'formats'
         ierr = 0
         nz = s% nz
         converged = .true.
         if (s% R_center > 0) then
            if (s% R_center > exp(s% lnR(nz))) then
               if (report) &
                  write(*,2) 'volume < 0 in cell nz', nz, &
                     s% R_center - exp(s% lnR(nz)), s% R_center, exp(s% lnR(nz)), &
                     s% dm(nz), s% rho(nz), s% dq(nz)
               converged = .false.
               return
            end if
         end if
         do k=1,nz-1
            if (s% lnR(k) <= s% lnR(k+1)) then
               if (report) write(*,2) 'after hydro, negative cell volume in cell k', &
                     k, s% lnR(k) - s% lnR(k+1), s% lnR(k), s% lnR(k+1), &
                     s% lnR_start(k) - s% lnR_start(k+1), s% lnR_start(k), s% lnR_start(k+1)
               converged = .false.; exit
               stop 'check_after_converge'
            else if (s% gamma_law_hydro <= 0d0) then
               if (s% lnT(k) > ln10*12) then
                  if (report) write(*,2) 'after hydro, logT > 12 in cell k', k, s% lnT(k)
                  converged = .false.!; exit
               else if (s% lnT(k) < ln10) then
                  if (report) write(*,*) 'after hydro, logT < 1 in cell k', k
                  converged = .false.!; exit
               else if (s% lnd(k) > ln10*12) then
                  if (report) write(*,*) 'after hydro, logRho > 12 in cell k', k
                  converged = .false.!; exit
               else if (s% lnd(k) < -ln10*30) then
                  if (report) write(*,*) 'after hydro, logRho < -30 in cell k', k
                  converged = .false.!; exit
               end if
            end if
         end do
      end function check_after_converge


      subroutine hydro_solver_step( &
            s, nz, nvar_hydro, nvar, dx1, dt, &
            gold_tolerances_level, tol_correction_norm, &
            solver_work, solver_lwork, &
            solver_iwork, solver_liwork, &
            converged, ierr)
         use num_def
         use chem_def
         use mtx_lib
         use mtx_def
         use alloc
         use hydro_mtx, only: ipar_id, ipar_first_call, hydro_lipar, &
            rpar_dt, hydro_lrpar

         type (star_info), pointer :: s
         integer, intent(in) :: nz, nvar_hydro, nvar
         real(dp), pointer :: dx1(:)
         real(dp), intent(in) :: dt
         real(dp), intent(in) :: tol_correction_norm
         integer, intent(in) :: gold_tolerances_level
         integer, intent(in) :: solver_lwork, solver_liwork
         real(dp), intent(inout), pointer :: solver_work(:) ! (solver_lwork)
         integer, intent(inout), pointer :: solver_iwork(:) ! (solver_liwork)
         logical, intent(out) :: converged
         integer, intent(out) :: ierr

         integer, parameter :: lipar=hydro_lipar, lrpar=hydro_lrpar
         integer, target :: ipar_target(lipar)
         real(dp), target :: rpar_target(lrpar)
         integer, pointer :: ipar(:)
         real(dp), pointer :: rpar(:)

         integer :: mljac, mujac, i, k, j, matrix_type, neq
         logical :: failure
         real(dp) :: varscale
         real(dp), parameter :: xscale_min = 1
         real(dp), pointer :: dx(:,:)

         real(dp), pointer, dimension(:,:) :: x_scale
         real(dp), pointer, dimension(:) :: x_scale1

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         neq = nvar*nz

         dx(1:nvar,1:nz) => dx1(1:neq)

         if (dbg) write(*, *) 'enter hydro_solver_step'

         s% used_extra_iter_in_solver_for_accretion = .false.

         mljac = 2*nvar-1
         mujac = mljac

         ipar => ipar_target
         ipar(ipar_id) = s% id
         ipar(ipar_first_call) = 1

         rpar => rpar_target
         rpar(rpar_dt) = dt

         call check_sizes(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'check_sizes failed'
            return
         end if

         call non_crit_get_work_array( &
            s, x_scale1, neq, nvar*nz_alloc_extra, 'hydro_solver_step', ierr)
         if (ierr /= 0) return
         x_scale(1:nvar,1:nz) => x_scale1(1:neq)

         do i = 1, nvar
            if (i <= s% nvar_hydro) then
               varscale = maxval(abs(s% xh(i,1:nz)))
               varscale = max(xscale_min, varscale)
            else
               varscale = 1
            end if
            x_scale(i, 1:nz) = varscale
         end do

         if (dbg) write(*, *) 'call solver'
         call newt(ierr)
         if (ierr /= 0 .and. s% report_ierr) then
            write(*,*) 'solver failed for hydro'
         end if

         converged = (ierr == 0) .and. (.not. failure)
         if (converged) then
            do k=1,nz
               do j=1,min(nvar,nvar_hydro)
                  s% xh(j,k) = s% xh_start(j,k) + dx(j,k)
               end do
            end do
            ! s% xa has already been updated by final call to set_solver_vars from solver
         end if

         call non_crit_return_work_array(s, x_scale1, 'hydro_solver_step')


         contains


         subroutine newt(ierr)
            use chem_def
            use hydro_solver_procs
            use rates_def, only: warn_rates_for_high_temp
            use star_utils, only: total_times
            use star_solver, only: solver
            use star_solver_15066, only: solver_15066
            use num_lib, only: default_failed_in_setmatrix, &
               default_set_primaries, default_set_secondaries
            integer, intent(out) :: ierr
            integer :: k, j
            logical :: save_warn_rates_flag
            include 'formats'
            solver_work(r_mtx_time) = 0
            solver_work(r_test_time) = 0
            solver_iwork(i_caller_id) = s% id
            s% doing_solver_iterations = .true.
            save_warn_rates_flag = warn_rates_for_high_temp
            warn_rates_for_high_temp = .false.        
            if (.not. s% conv_vel_flag) then
               call solver( &
                  s, nz, nvar, dx1, &
                  gold_tolerances_level, tol_correction_norm, &
                  x_scale1, s% equ1, &
                  solver_work, solver_lwork, &
                  solver_iwork, solver_liwork, &
                  s% AF1, lrpar, rpar, lipar, ipar, failure, ierr)
            else
               call solver_15066( &
                  s, nz, nvar, dx1, &
                  gold_tolerances_level, tol_correction_norm, &
                  x_scale1, s% equ1, &
                  solver_work, solver_lwork, &
                  solver_iwork, solver_liwork, &
                  s% AF1, lrpar, rpar, lipar, ipar, failure, ierr)
            end if    
            s% doing_solver_iterations = .false.
            warn_rates_for_high_temp = save_warn_rates_flag
         end subroutine newt


      end subroutine hydro_solver_step


      subroutine set_L_burn_by_category(s)
         use rates_def, only: i_rate
         type (star_info), pointer :: s
         integer :: k, j
         real(dp) :: L_burn_by_category(num_categories)
         L_burn_by_category(:) = 0
         do k = s% nz, 1, -1
            do j = 1, num_categories
               L_burn_by_category(j) = &
                  L_burn_by_category(j) + s% dm(k)*s% eps_nuc_categories(j, k)
               s% luminosity_by_category(j,k) = L_burn_by_category(j)
            end do
         end do
      end subroutine set_L_burn_by_category


      end module solve_hydro


