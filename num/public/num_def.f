! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module num_def
      
      use const_def, only: dp
      
      implicit none

      
      ! these are the stiff ode solvers
      integer, parameter :: ros2_solver = 1
      integer, parameter :: rose2_solver = 2
      integer, parameter :: ros3p_solver = 3
      integer, parameter :: ros3pl_solver = 4
      integer, parameter :: rodas3_solver = 5
      integer, parameter :: rodas4_solver = 6
      integer, parameter :: rodasp_solver = 7
      
      integer, parameter :: num_solvers = 7
      
      
      ! operation codes for nelder-mead simplex function calls
      integer, parameter :: simplex_initial = 1 ! 'initial'
      integer, parameter :: simplex_reflect = 2 ! 'reflect'
      integer, parameter :: simplex_expand = 3 ! 'expand'
      integer, parameter :: simplex_inside = 4 ! 'inside'
      integer, parameter :: simplex_outside = 5 ! 'outside'
      integer, parameter :: simplex_random = 6 ! 'random'
      integer, parameter :: simplex_shrink = 7 ! 'shrink'
      

      ! these are the matrix type options
      integer, parameter :: square_matrix_type = 1
      integer, parameter :: banded_matrix_type = 2
      integer, parameter :: block_tridiag_dble_matrix_type = 3
      integer, parameter :: block_tridiag_quad_matrix_type = 4
      
      
      ! parameter indices in newton iwork array
      
         integer, parameter :: i_try_really_hard = 1
            ! if "try_really_hard" is nonzero, it tells newton that the caller 
            ! has no backup in case newton fails to converge, 
            ! so extraordinary efforts are appropriate.
         integer, parameter :: i_itermin = i_try_really_hard+1
            ! min number of iterations in n-r
            ! the default value (for iwork(i_itermin)=0) is 2.
         integer, parameter :: i_max_tries = i_itermin+1
            ! max number of iterations in n-r
            ! the default value (for iwork(i_max_tries)=0) is 50.
         integer, parameter :: i_tiny_min_corr_coeff = i_max_tries+1
            ! give up if get a sequence of this many tiny corrections (only applies when min_corr_coeff < 1).
            ! the default value (for iwork(i_tiny_min_corr_coeff)=0) is 25
         integer, parameter :: i_debug = i_tiny_min_corr_coeff+1
            ! if nonzero, then get lots of debugging output from the solver
         integer, parameter :: i_do_core_dump = i_debug+1
            ! if nonzero, huge amounts of stuff are written to the terminal
         integer, parameter :: i_model_number = i_do_core_dump+1
         integer, parameter :: i_num_solves = i_model_number+1
         integer, parameter :: i_num_jacobians = i_num_solves+1
         
         ! the line search routine considers only variables lsvar_lo:lsvar_hi
         ! if iwork(i_lsvar_lo) = 0, it defaults to 1
         ! if iwork(i_lsvar_hi) = 0, it defaults to the number of variables
         integer, parameter :: i_lsvar_lo = i_num_jacobians+1  ! OBSOLETE
         integer, parameter :: i_lsvar_hi = i_lsvar_lo+1  ! OBSOLETE
                  
         integer, parameter :: i_max_iter_for_enforce_resid_tol = i_lsvar_hi+1
            ! use tol_residual_norm & tol_max_residual
            ! for iters from 1 to max_iter_for_enforce_resid_tol
         
         integer, parameter :: i_max_iter_for_resid_tol2 = i_max_iter_for_enforce_resid_tol+1
            ! use tol_residual_norm2 & tol_max_residual2 
            ! for iters from max_iter_for_enforce_resid_tol+1 to max_iter_for_resid_tol2
         
         integer, parameter :: i_max_iter_for_resid_tol3 = i_max_iter_for_resid_tol2+1
            ! use tol_residual_norm3 & tol_max_residual3 
            ! for iters from max_iter_for_resid_tol2+1 to max_iter_for_resid_tol3
            ! for iters > max_iter_for_resid_tol3, ignore residuals for convergence decision
      
         integer, parameter :: i_max_iterations_for_jacobian = i_max_iter_for_resid_tol3+1
            ! jacobian is always created fresh for 1st iteration.
            ! if this param > 1, then will try to reuse jacobian.
            ! after use jacobian this many times, remake it.
            ! e.g., if = 2, then will make a new jacobian for every other iteration.
      
         integer, parameter :: i_refine_solution = i_max_iterations_for_jacobian+1
            ! if this is non-zero, then allow code to do an extra iteration with same jacobian
            ! to refine an already converged solution.
            ! it will do the refinement only if the correction norm is > 0.1 * tolerance,
            ! and if it hasn't already reached the max allowed number of iterations,
            ! and if the previous iteration was done using a new jacobian (i.e., won't use old J to refine). 
      
         integer, parameter :: i_refine_mtx_solution = i_refine_solution+1
            ! if this is non-zero, then refine the linear algebra solution of A*x = b
            
         integer, parameter :: i_min_for_check_D_norm_converging = i_refine_mtx_solution+1
            ! once iter reaches this limit, start checking to see if D_norm is converging fast enough.
      
         integer, parameter :: i_caller_id = i_min_for_check_D_norm_converging+1

         integer, parameter :: num_iwork_params = i_caller_id
      
      
      ! parameter indices in newton work array
      
         integer, parameter :: r_tol_residual_norm = 1
            ! require avg_residual < tol_residual_norm before consider solution okay
            ! residual_norm is calculated by the size_equ routine.
            ! the default value (for work(r_tol_residual_norm)=0) is 1d99
         integer, parameter :: r_tol_max_residual = r_tol_residual_norm+1
            ! require max_residual < tol_max_residual before consider solution okay
            ! max_residual is calculated by the size_equ routine.
            ! the default value (for work(r_tol_max_residual)=0) is 1d99            
         integer, parameter :: r_tol_residual_norm2 = r_tol_max_residual+1
         integer, parameter :: r_tol_max_residual2 = r_tol_residual_norm2+1
         integer, parameter :: r_tol_residual_norm3 = r_tol_max_residual2+1
         integer, parameter :: r_tol_max_residual3 = r_tol_residual_norm3+1
         integer, parameter :: r_tol_max_correction = r_tol_max_residual3+1
            ! require max correction < tol_max_correction before consider solution okay
            ! max_correction is calculated by the sizeB routine.
            ! the default value (for work(r_tol_max_correction)=0) is 1d99
         integer, parameter :: r_target_corr_factor = r_tol_max_correction+1
            ! if correction_norm > target_corr_factor*corr_norm_previous and not first iteration and
            !     not using a freshly made jacobian, then make a new one and try again
            ! the default value (for work(r_target_corr_factor)=0) is 0.9d0
         integer, parameter :: r_tol_abs_slope_min = r_target_corr_factor+1
            ! if abs(slope) <= tol_abs_slope_min, then accept solution as converged
         integer, parameter :: r_tol_corr_resid_product = r_tol_abs_slope_min+1
            ! if correction_norm*residual_norm <= tol_corr_resid_product, then accept solution as converged
         integer, parameter :: r_scale_correction_norm = r_tol_corr_resid_product+1
            ! if correction_norm > corr_param_factor*scale_correction_norm, then
            !     if not using a new jacobian, make one and retry
            !     else consider such a large correction to be an indication of trouble and give up.
            ! if correction_norm > scale_correction_norm, then
            !     rescale the correction vector by scale_correction_norm/correction_norm
            ! the default value (for work(r_scale_correction_norm)=0) is 2d0
         integer, parameter :: r_corr_param_factor = r_scale_correction_norm+1
            ! used along with scale_correction_norm as described above
            ! the default value (for work(r_corr_param_factor)=0) is 10d0
         integer, parameter :: r_scale_max_correction = r_corr_param_factor+1
            ! if max_correction > scale_max_correction, then
            !     rescale the correction vector by scale_max_correction/max_correction
            ! the default value (for work(r_scale_max_correction)=0) is 1d99
         integer, parameter :: r_corr_norm_jump_limit = r_scale_max_correction+1
            ! if not try_really_hard and not first iteration 
            !        and correction_norm > corr_norm_jump_limit*corr_norm_previous then
            !     if not using a new jacobian, make one and retry
            !     else consider such a large jump in the correction to be an indication of trouble and give up.
            ! the default value (for work(r_corr_norm_jump_limit)=0) is 1d99
         integer, parameter :: r_max_corr_jump_limit = r_corr_norm_jump_limit+1
            ! same as corr_norm_jump_limit, but for the max correction instead.
            ! the default value (for work(r_max_corr_jump_limit)=0) is 1d99
         integer, parameter :: r_resid_norm_jump_limit = r_max_corr_jump_limit+1
            ! if not try_really_hard and not first iteration 
            !        and residual_norm > resid_norm_jump_limit*resid_norm_previous then
            !     if not using a new jacobian, make one and retry
            !     else consider such a large jump in the residual to be an indication of trouble and give up.
            ! the default value (for work(r_resid_norm_jump_limit)=0) is 1d99
         integer, parameter :: r_max_resid_jump_limit = r_resid_norm_jump_limit+1
            ! same as resid_norm_jump_limit, but for the max residual instead.
            ! the default value (for work(r_max_resid_jump_limit)=0) is 1d99
         integer, parameter :: r_min_corr_coeff = r_max_resid_jump_limit+1
            ! lower bound on correction coefficient for globally convergent newton procedure
            ! you can set this to 1 to disable the global convergence parts of the solver
            ! the default value (for work(r_min_corr_coeff)=0) is 1d-3
         integer, parameter :: r_slope_alert_level = r_min_corr_coeff+1
            ! "slope" here indicates the change in the residuals in the direction of the correction, 
            !     and therefore it should be negative.
            ! a positive slope is indication of substantial error in the approximate jacobian
            !     (with an exact jacobian, the slope would never be positive)
            ! when near a root, the slope becomes very tiny and roundoff can give
            !     a tiny positive slope instead of a tiny negative one.
            ! so we ignore these cases and only look for slope > slope_limit
            !     as a signal of needing a better jacobian
            ! the default value (for work(r_slope_alert_level)=0) is 1d0
         integer, parameter :: r_slope_crisis_level = r_slope_alert_level+1
            ! similar to slope_alert_level, but applies when already using a new jacobian.
            ! in this case, if find slope > slope_crisis_level, give up.
            ! the default value (for work(r_slope_crisis_level)=0) is 1d0
         integer, parameter :: r_tiny_corr_factor = r_slope_crisis_level+1
            ! a correction counts as tiny if correction coeff < this factor times the min_corr_coeff 
            ! the default value (for work(r_tiny_corr_factor)=0) is 2d0
         integer, parameter :: r_dt = r_tiny_corr_factor+1
            ! current timestep. (just for debugging output).
         integer, parameter :: r_sparse_non_zero_max_factor = r_dt+1
         
         integer, parameter :: r_mtx_time = r_sparse_non_zero_max_factor+1
            ! if nonzero, then return time spent in mtx routines (in seconds)
         integer, parameter :: r_test_time = r_mtx_time+1
         
         integer, parameter :: r_D_norm_kappa = r_test_time+1
         integer, parameter :: r_D_norm_err_est = r_D_norm_kappa+1
         
         integer, parameter :: num_work_params = r_D_norm_err_est


         ! Interface for integrators
         abstract interface
            real(dp) function integrand(x, args, ierr)
               ! Evalaute the function at point x, with possible extra args and return error in ierr
               use const_def, only: dp
               implicit none
               real(dp), intent(in) :: x
               real(dp), intent(in) :: args(:)
               integer, intent(inout) :: ierr
            end function integrand
         end interface


      end module num_def

