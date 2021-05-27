! ***********************************************************************
!
!   Copyright (C) 2012-2019  The MESA Team
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


      module mod_newton

      use const_def, only: dp, qp
      use math_lib
      use utils_lib, only: is_bad, mesa_error
      use num_def
      use mtx_def
      use mtx_lib, only: band_multiply_xa,
     >   lapack_work_sizes, 
     >   block_multiply_xa,
     >   multiply_xa
      
      implicit none
      
      
      
      contains


      subroutine do_newton_wrapper(
     >   nz, nvar, x, xold,
     >   matrix_type, mljac, mujac,
     >   decsol, decsolblk, 
     >   lrd, rpar_decsol, lid, ipar_decsol, which_decsol,
     >   tol_correction_norm,
     >   set_primaries, set_secondaries, set_xscale, Bdomain, xdomain, eval_equations,
     >   size_equ, sizeb, inspectB,
     >   enter_setmatrix, exit_setmatrix, failed_in_setmatrix, force_another_iteration,
     >   xscale, equ, ldy, nsec, y, work, lwork, iwork, liwork, AF,
     >   lrpar, rpar, lipar, ipar, convergence_failure, ierr)
         
         ! the primary variables
         integer, intent(in) :: nz ! number of zones
         integer, intent(in) :: nvar ! number of variables per zone
         ! the total number of primary variables is neq
         real(dp), pointer, dimension(:) :: x ! =(nvar,nz) 
         ! new vector of primaries
         real(dp), pointer, dimension(:) :: xold ! =(nvar,nz)
         ! old vector of primaries

         ! information about the jacobian matrix
         integer, intent(in) :: matrix_type ! see num_def.f for values
         ! if matrix_type == banded_matrix_type, mljac and mujac give the bandwidths.
         ! if matrix_type /= banded_matrix_type, mljac and mujac are not used.

         integer, intent(in) :: mljac ! number of subdiagonals within the band of the jacobian
         integer, intent(in) :: mujac ! number of superdiagonals
         ! for example if you have a centered 3 zone stencil, 
         ! then you will have mljac and mujac both equal to 2*nvar-1
         ! mljac and mujac are only used for matrix_type == banded matrix type
                  
         ! matrix routines
         ! there are implementations of the matrix routines available in mesa/mtx.
         ! for example, the LAPACK versions are called lapack_dec, lapack_sol, etc.
         interface
#include "mtx_decsol.dek"
#include "mtx_decsolblk_dble.dek"
         end interface
         ! these arrays provide optional extra working storage for the matrix routines.
         ! the implementations in mesa/mtx include routines to determine the sizes.
         ! for example, the LAPACK version is called lapack_work_sizes.
         integer, intent(in) :: lrd, lid
         integer, intent(inout), pointer :: ipar_decsol(:) ! (lid)
         real(dp), intent(inout), pointer :: rpar_decsol(:) ! (lrd)
         integer, intent(in) :: which_decsol

         real(dp), pointer, dimension(:) :: xscale ! =(nvar,nz)
         ! typical values for x.  set by set_xscale.
         real(dp), pointer, dimension(:) :: equ ! =(nvar,nz)
         ! equ(i) has the residual for equation i, i.e., the difference between
         ! the left and right hand sides of the equation.

         ! the secondary variables
            ! the "secondaries" for zone k depend only on the primaries of zone k
            ! and therefore need not be recomputed is the zone k primaries have not been modified.
            ! using this information can significantly accelerate the computation of numerical jacobians.
            ! for stellar evolution, the secondaries include such expensive-to-compute items
            ! as equation of state, 
            ! nuclear reaction rates, and opacities.
         integer, intent(in) :: ldy ! leading dimension of y, >= nz
         integer, intent(in) :: nsec ! number of secondaries per zone
         real(dp), pointer, dimension(:) :: y ! the values. =(ldy,nsec)

         ! work arrays. required sizes provided by the routine newton_work_sizes.
         ! for standard use, set work and iwork to 0 before calling.
         ! NOTE: these arrays contain some optional parameter settings and outputs.
         ! see num_def for details.
         integer, intent(in) :: lwork, liwork
         real(dp), intent(inout), target :: work(:) ! (lwork)
         integer, intent(inout), target :: iwork(:) ! (liwork)
         real(dp), pointer, dimension(:) :: AF ! for factored jacobian
            ! will be allocated or reallocated as necessary.  

         ! convergence criteria
         real(dp), intent(in) :: tol_correction_norm
            ! a trial solution is considered to have converged if
            ! max_correction <= tol_max_correction and
            !
            ! either
            !          (correction_norm <= tol_correction_norm)  
            !    .and. (residual_norm <= tol_residual_norm)
            ! or
            !          (correction_norm*residual_norm <= tol_corr_resid_product)
            !    .and. (abs(slope) <= tol_abs_slope_min)
            !
            ! where "slope" is slope of the line for line search in the newton solver,
            ! and is analogous to the slope of df/dx in a 1D newton root finder.

         ! parameters for caller-supplied routines
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         
         ! output
         logical, intent(out) :: convergence_failure
         integer, intent(out) :: ierr ! 0 means okay.
         
         ! the following routines implement the problem-specific aspects of the newton solver.
         ! see num/include/newton_procs.dek for documentation.
         ! there are default implementations for most of the routines (see below).
         ! the only one without a default is the "eval_equations" routine that computes
         ! the equation residuals for your particular problem.
         interface
#include "newton_procs.dek" 
         end interface
         
         integer :: ldAF, neqns
         real(dp), pointer :: AF_copy(:) ! =(ldAF, neq)
         
         ! for sparse
         integer :: n, nzmax, need_lrd, need_lid, lfil, maxits, iout
         real(dp) :: eps, droptol
               
         real(dp) :: sparse_nzmax_factor
         integer(8) :: test_time0, test_time1, clock_rate
         logical :: do_test_timing
         
         include 'formats'

         do_test_timing = (work(r_test_time) /= 0)
         work(r_test_time) = 0

         ierr = 0
         
         nzmax = 0
         if (which_decsol == lapack) then
            call lapack_work_sizes(n, need_lrd, need_lid)
         else
            write(*,*) 'newton: unknown value for matrix solver option'
            ierr = -1
            return 
         end if
         
         if (need_lrd > lrd .or. need_lid > lid) then
            write(*,*) 'bad lrd or lid for newton'
            write(*,2) 'need_lrd', need_lrd
            write(*,2) '     lrd', lrd
            write(*,2) 'need_lid', need_lid
            write(*,2) '     lid', lid
            ierr = -1
            return
         end if

         neqns = nvar*nz
         if (matrix_type == block_tridiag_dble_matrix_type) then
            ldAF = 3*nvar
         else if (matrix_type == square_matrix_type) then
            ldAF = neqns
         else
            ldAF = 2*mljac+mujac+1
         end if
         
         if (associated(AF)) then
            if (size(AF,dim=1) < ldAF*neqns) then
               deallocate(AF)
               nullify(AF)
            end if
         end if
         
         if (.not. associated(AF)) then
            allocate(AF((ldAF+2)*(neqns+200)), stat=ierr)
            if (ierr /= 0) return
         end if
         
         AF_copy => AF

         if (do_test_timing) call system_clock(test_time0,clock_rate)
         
         call do_newton(
     >      nz, nvar, x, xold, AF_copy, ldAF, neqns, 
     >      matrix_type, mljac, mujac,
     >      decsol, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,
     >      tol_correction_norm,
     >      set_primaries, set_secondaries, set_xscale, Bdomain, xdomain, eval_equations,
     >      size_equ, sizeb, inspectB,
     >      enter_setmatrix, exit_setmatrix, failed_in_setmatrix, force_another_iteration,
     >      xscale, equ, ldy, nsec, y, work, lwork, iwork, liwork,
     >      lrpar, rpar, lipar, ipar, convergence_failure, ierr)

         if (do_test_timing) then
            call system_clock(test_time1,clock_rate)
            work(r_test_time) = work(r_test_time) + dble(test_time1 - test_time0) / clock_rate
         end if
        
         
         contains
         
               
         logical function bad_isize(a,sz,str)
            integer :: a(:)
            integer, intent(in) :: sz
            character (len=*), intent(in) :: str
            bad_isize = (size(a,dim=1) < sz)
            if (.not. bad_isize) return
            ierr = -1
            return
         end function bad_isize
         
      
         logical function bad_size(a,sz,str)
            real(dp) :: a(:)
            integer, intent(in) :: sz
            character (len=*), intent(in) :: str
            bad_size = (size(a,dim=1) < sz)
            if (.not. bad_size) return
            ierr = -1
            return
         end function bad_size
         
      
         logical function bad_size_dble(a,sz,str)
            real(dp) :: a(:)
            integer, intent(in) :: sz
            character (len=*), intent(in) :: str
            bad_size_dble = (size(a,dim=1) < sz)
            if (.not. bad_size_dble) return
            ierr = -1
            return
         end function bad_size_dble
         
      
         logical function bad_sizes(a,sz1,sz2,str)
            real(dp) :: a(:,:)
            integer, intent(in) :: sz1,sz2
            character (len=*), intent(in) :: str
            bad_sizes = (size(a,dim=1) < sz1 .or. size(a,dim=2) < sz2)
            if (.not. bad_sizes) return
            ierr = -1
            return
         end function bad_sizes
         
         
      end subroutine do_newton_wrapper


      subroutine do_newton(
     >   nz, nvar, x1, xold1, AF1, ldAF, neq,
     >   matrix_type, mljac, mujac, 
     >   decsol, decsolblk, 
     >   lrd, rpar_decsol, lid, ipar_decsol,
     >   tol_correction_norm,
     >   set_primaries, set_secondaries, set_xscale, Bdomain, xdomain, eval_equations,
     >   size_equ, sizeb, inspectB,
     >   enter_setmatrix, exit_setmatrix, failed_in_setmatrix, force_another_iteration,
     >   xscale1, equ1, ldy, nsec, y_in1, 
     >   work, lwork, iwork, liwork, 
     >   lrpar, rpar, lipar, ipar, convergence_failure, ierr)

         integer, intent(in) :: nz, nvar, mljac, mujac, ldy, nsec, ldAF, neq
         
         integer, intent(in) :: matrix_type

         real(dp), pointer, dimension(:) :: AF1 ! =(ldAF, neq), neq = neq
         real(dp), pointer, dimension(:) :: x1, xold1, equ1, xscale1 
         real(dp), pointer, dimension(:) :: y_in1 ! the values. =(ldy,nsec)
                           
         ! matrix routines
         interface
#include "mtx_decsol.dek"
#include "mtx_decsolblk_dble.dek"
         end interface
         integer, intent(in) :: lrd, lid
         integer, intent(inout), pointer :: ipar_decsol(:) ! (lid)
         real(dp), intent(inout), pointer :: rpar_decsol(:) ! (lrd)

         ! controls         
         real(dp), intent(in) :: tol_correction_norm

         ! parameters for caller-supplied routines
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)

         ! work arrays
         integer, intent(in) :: lwork, liwork
         real(dp), intent(inout), target :: work(:) ! (lwork)
         integer, intent(inout), target :: iwork(:) ! (liwork)

         ! output
         logical, intent(out) :: convergence_failure
         integer, intent(out) :: ierr
         
         ! procedures
         interface
#include "newton_procs.dek" 
         end interface

         ! info saved in work and iwork
         real(dp), dimension(:,:), pointer :: A, Acopy
         real(dp), dimension(:), pointer :: A1, Acopy1
         
         real(dp), dimension(:,:), pointer :: xsave, dxsave, B, grad_f, B_init
         real(dp), dimension(:), pointer :: xsave1, dxsave1, B1, B_init1, grad_f1
         real(dp), dimension(:,:), pointer ::  rhs
         integer, dimension(:), pointer :: ipiv1
         real(dp), dimension(:,:), pointer :: dx, xgg, dxd, dxdd, xder, equsave
         real(dp), dimension(:,:), pointer :: y1, y2
         
         real(dp), dimension(:), pointer :: lblk1, dblk1, ublk1
         real(dp), dimension(:), pointer :: lblkF1, dblkF1, ublkF1
         integer, dimension(:), pointer :: ipiv_blk1
         
         ! locals
         real(dp)  :: 
     >      coeff, f, slope, residual_norm, max_residual, corr_norm_min, resid_norm_min, correction_factor,
     >      residual_norm_save, corr_norm_min_save, resid_norm_min_save, correction_factor_save,
     >      correction_norm, corr_norm_initial, max_correction, slope_extra,
     >      tol_max_correction, tol_residual_norm, tol_abs_slope_min, tol_corr_resid_product,
     >      min_corr_coeff, tol_max_residual, max_corr_min, max_resid_min
         integer :: iiter, max_tries, ndiag, zone, idiag, tiny_corr_cnt, ldA, i, j, k, info,
     >      last_jac_iter, max_iterations_for_jacobian, force_iter_value,
     >      max_iter_for_enforce_resid_tol, max_iter_for_resid_tol2, max_iter_for_resid_tol3,
     >      caller_id
         integer(8) :: test_time0, test_time1, time0, time1, clock_rate
         character (len=256) :: err_msg
         logical :: first_try, dbg_msg, passed_tol_tests, 
     >      overlay_AF, do_mtx_timing, do_test_timing, doing_extra
         integer, parameter :: num_tol_msgs = 15
         character (len=32) :: tol_msg(num_tol_msgs)
         character (len=64) :: message
         real(dp), pointer, dimension(:) :: p1_1, p1_2

         ! set pointers to 1D data
         real(dp), pointer, dimension(:,:) :: x, xold, equ, xscale ! (nvar,nz)       
         real(dp), pointer, dimension(:,:) :: y ! (ldy,nsec)
         real(dp), pointer, dimension(:,:) :: AF ! (ldAF,neq)
         real(dp), pointer, dimension(:,:,:) :: ublk, dblk, lblk ! (nvar,nvar,nz)
         real(dp), dimension(:,:,:), pointer :: lblkF, dblkF, ublkF ! (nvar,nvar,nz)
         
         include 'formats'
         
         x(1:nvar,1:nz) => x1(1:neq)
         xold(1:nvar,1:nz) => xold1(1:neq)
         equ(1:nvar,1:nz) => equ1(1:neq)
         xscale(1:nvar,1:nz) => xscale1(1:neq)
         AF(1:ldAF,1:neq) => AF1(1:ldAF*neq)
         y(1:ldy,1:nsec) => y_in1(1:ldy*nsec)
         
         do_mtx_timing = (work(r_mtx_time) /= 0)
         work(r_mtx_time) = 0

         tol_msg(1)  = 'avg corr'
         tol_msg(2)  = 'max corr '
         tol_msg(3)  = 'avg+max corr'
         tol_msg(4)  = 'avg resid'
         tol_msg(5)  = 'avg corr+resid'
         tol_msg(6)  = 'max corr, avg resid'
         tol_msg(7)  = 'avg+max corr, avg resid'
         tol_msg(8)  = 'max resid'
         tol_msg(9)  = 'avg corr, max resid'
         tol_msg(10) = 'max corr+resid'
         tol_msg(11) = 'avg+max corr, max resid'
         tol_msg(12) = 'avg+max resid'
         tol_msg(13) = 'avg corr, avg+max resid'
         tol_msg(14) = 'max corr, avg+max resid'
         tol_msg(15) = 'avg+max corr+resid'
 
         ierr = 0
         iiter = 0

         call set_param_defaults
         dbg_msg = (iwork(i_debug) /= 0)
         tol_residual_norm = work(r_tol_residual_norm)
         tol_max_residual = work(r_tol_max_residual)
         tol_max_correction = work(r_tol_max_correction)
         tol_abs_slope_min = work(r_tol_abs_slope_min)
         tol_corr_resid_product = work(r_tol_corr_resid_product)
         min_corr_coeff = work(r_min_corr_coeff)
         
         max_iter_for_enforce_resid_tol = iwork(i_max_iter_for_enforce_resid_tol)
         max_iter_for_resid_tol2 = iwork(i_max_iter_for_resid_tol2)
         max_iter_for_resid_tol3 = iwork(i_max_iter_for_resid_tol3)
         
         caller_id = iwork(i_caller_id)
         
         if (ldy < nz .and. nsec > 0) then
            ierr = -1
            return
         end if
         
         idiag = 1
         if (matrix_type == block_tridiag_dble_matrix_type) then
            ndiag = 3*nvar
         else if (matrix_type == square_matrix_type) then
            ndiag = neq
         else
            idiag = mujac+1
            ndiag = mljac+mujac+1
         end if
            
         ldA = ndiag
         call pointers(ierr)
         if (ierr /= 0) return

         if (iwork(i_do_core_dump) /= 0) then
            call newton_core_dump(x, dx, xold)
            return
         end if
      
         doing_extra = .false.
         passed_tol_tests = .false. ! goes true when pass the tests
         convergence_failure = .false. ! goes true when time to give up
         coeff = 1.
         xscale = 1.
  
         residual_norm=0
         max_residual=0
         corr_norm_min=1d99
         max_corr_min=1d99
         max_resid_min=1d99
         resid_norm_min=1d99
         correction_factor=0
         
         do k=1,nz
            do i=1,nvar
               dx(i,k) = x(i,k) - xold(i,k)
            end do
         end do
         
         call xdomain(iiter, nvar, nz, x, dx, xold, lrpar, rpar, lipar, ipar, ierr)
         if (ierr /= 0) then
            if (dbg_msg)
     >         write(*, *) 'newton failure: xdomain returned ierr', ierr
            convergence_failure = .true.
            return
         end if
         call set_xscale(nvar, nz, xold, xscale, lrpar, rpar, lipar, ipar, ierr) ! set xscale
         if (ierr /= 0) then
            if (dbg_msg)
     >         write(*, *) 'newton failure: set_xscale returned ierr', ierr
            convergence_failure = .true.
            return
         end if
         call setequ(nvar, nz, x, equ, lrpar, rpar, lipar, ipar, ierr)
         if (ierr /= 0) then
            if (dbg_msg)
     >         write(*, *) 'newton failure: setequ returned ierr', ierr
            convergence_failure = .true.
            return
         end if
         call size_equ(
     >      iiter, nvar, nz, equ, residual_norm, max_residual, lrpar, rpar, lipar, ipar, ierr)
         if (ierr /= 0) then
            if (dbg_msg)
     >         write(*, *) 'newton failure: size_equ returned ierr', ierr
            convergence_failure = .true.
            return
         end if

         first_try = .true.
         iiter = 1
         max_tries = abs(iwork(i_max_tries))
         last_jac_iter = 0
         tiny_corr_cnt = 0
         
         if (iwork(i_max_iterations_for_jacobian) == 0) then
            max_iterations_for_jacobian = 1000000
         else
            max_iterations_for_jacobian = iwork(i_max_iterations_for_jacobian)
         end if

         do while (.not. passed_tol_tests)
            
            if (dbg_msg .and. first_try) write(*, *)
                  
            if (iiter >= max_iter_for_enforce_resid_tol) then
               if (iiter >= max_iter_for_resid_tol2) then
                  if (iiter >= max_iter_for_resid_tol3) then ! shut down
                     tol_residual_norm = 1d200
                     tol_max_residual = 1d200
                  else ! >= max_iter_for_resid_tol2 and but < max_iter_for_resid_tol3
                     tol_residual_norm = work(r_tol_residual_norm3)
                     tol_max_residual = work(r_tol_max_residual3)
                  end if
               else ! >= max_iter_for_enforce_resid_tol but < max_iter_for_resid_tol2
                  tol_residual_norm = work(r_tol_residual_norm2)
                  tol_max_residual = work(r_tol_max_residual2)
               end if
            end if

            overlay_AF = (min_corr_coeff == 1) .and. 
     >            (matrix_type == banded_matrix_type .or.            
     >               matrix_type == block_tridiag_dble_matrix_type)
            
            ! NOTE: for banded matrix, the jacobian A is a part of the array AF
            ! AF has extra rows for storing banded LU factored matrix.
            if (overlay_AF) then
               A1 => AF1
               A => AF
               ldA = ldAF
               if (matrix_type == banded_matrix_type) then
                  idiag = mljac+mujac+1
               else if (matrix_type == block_tridiag_dble_matrix_type) then
                  ublk1 => ublkF1
                  dblk1 => dblkF1
                  lblk1 => lblkF1
                  lblk(1:nvar,1:nvar,1:nz) => lblk1(1:nvar*neq)
                  dblk(1:nvar,1:nvar,1:nz) => dblk1(1:nvar*neq)
                  ublk(1:nvar,1:nvar,1:nz) => ublk1(1:nvar*neq)
               else
                  stop 'confusion about matrix_type'
               end if
            else
               idiag = mujac+1
            end if

            call setmatrix(neq, x, dx, xscale, xsave, dxsave, lrpar, rpar, lipar, ipar, ierr)
            if (ierr /= 0) then
               if (any(dx /= 0)) then
                  call write_msg('setmatrix returned ierr /= 0; retry with dx = 0.')
                  do i=1,neq
                     x1(i) = xold1(i)
                  end do
                  dx = 0
                  iiter=iiter+1
                  first_try = .false.
                  cycle
               end if
               call write_msg('setmatrix returned ierr /= 0; dx already = 0, so give up.')
               convergence_failure = .true.; exit
            end if
            iwork(i_num_jacobians) = iwork(i_num_jacobians) + 1
            last_jac_iter = iiter
            
            if (.not. solve_equ()) then ! either singular or horribly ill-conditioned
               write(err_msg, '(a, i5, 3x, a)') 'info', ierr, 'bad_matrix'
               call oops(err_msg)
               exit
            end if
            iwork(i_num_solves) = iwork(i_num_solves) + 1

            ! inform caller about the correction
            call inspectB(iiter, nvar, nz, x, B, xscale, lrpar, rpar, lipar, ipar, ierr)
            if (ierr /= 0) then
               call oops('inspectB returned ierr')
               exit
            end if

            ! compute size of scaled correction B
            call sizeB(iiter, nvar, nz, x, B, xscale, max_correction, correction_norm, 
     >               lrpar, rpar, lipar, ipar, ierr)
            if (ierr /= 0) then
               call oops('correction rejected by sizeB')
               exit
            end if
            correction_norm = abs(correction_norm)
            max_correction = abs(max_correction)
            corr_norm_min = min(correction_norm, corr_norm_min)
            max_corr_min = min(max_correction, max_corr_min)

            if (is_bad(correction_norm) .or. is_bad(max_correction)) then 
               ! bad news -- bogus correction
               call oops('bad result from sizeb -- correction info either NaN or Inf')
               exit
            end if

            if ((correction_norm > work(r_corr_param_factor)*work(r_scale_correction_norm)) .and.
     >            (iwork(i_try_really_hard) == 0)) then
               call oops('avg corr too large')
               exit
            endif
         
            ! shrink the correction if it is too large
            correction_factor = 1
            
            if (correction_norm*correction_factor > work(r_scale_correction_norm)) then
               correction_factor = work(r_scale_correction_norm)/correction_norm
            end if
            
            if (max_correction*correction_factor > work(r_scale_max_correction)) then
               correction_factor = work(r_scale_max_correction)/max_correction
            end if
            
            ! fix B if out of definition domain
            call Bdomain(
     >         iiter, nvar, nz, B, x, xscale, correction_factor, lrpar, rpar, lipar, ipar, ierr)
            if (ierr /= 0) then ! correction cannot be fixed
               call oops('correction rejected by Bdomain')
               exit
            end if

            ! save previous
            residual_norm_save = residual_norm
            corr_norm_min_save = corr_norm_min
            resid_norm_min_save = resid_norm_min
            correction_factor_save = correction_factor

            if (min_corr_coeff < 1) then
               ! compute gradient of f = equ<dot>jacobian
               ! NOTE: NOT jacobian<dot>equ
               if (matrix_type == block_tridiag_dble_matrix_type) then
                  call block_multiply_xa(nvar, nz, lblk1, dblk1, ublk1, equ1, grad_f1)                  
               else if (matrix_type == square_matrix_type) then
                  call multiply_xa(neq, A1, equ1, grad_f1)
               else
                  call band_multiply_xa(neq, mljac, mujac, A1, ldA, equ1, grad_f1)
               end if
            
               slope = eval_slope(nvar, nz, grad_f, B)
               !write(*,*) 'slope', slope
               !if (is_bad(slope)) then
               !   call oops('bad slope value')
               !   exit
               !end if
               if (is_bad(slope) .or. slope > 0) then ! a bad sign
                  ! but give it a chance before give up
                  !write(*,*) 'slope', slope
                  slope = 0
                  min_corr_coeff = 1
               end if
               
            else
            
               slope = 0

            end if
      
            call adjust_correction(
     >         min_corr_coeff, correction_factor, grad_f1, f, slope, coeff,
     >         err_msg, lrpar, rpar, lipar, ipar, ierr)
            if (ierr /= 0) then
               call oops(err_msg)
               exit
            end if
            
            ! coeff is factor by which adjust_correction rescaled the correction vector
            if (coeff > work(r_tiny_corr_factor)*min_corr_coeff) then
               tiny_corr_cnt = 0
            else
               tiny_corr_cnt = tiny_corr_cnt + 1
            end if

            ! check the residuals for the equations
            call size_equ(iiter, nvar, nz, equ, residual_norm, max_residual, lrpar, rpar, lipar, ipar, ierr)
            if (ierr /= 0) then
               call oops('size_equ returned ierr')
               exit
            end if
            if (is_bad(residual_norm)) then
               call oops('residual_norm is a a bad number (NaN or Infinity)')
               exit
            end if
            if (is_bad(max_residual)) then
               call oops('max_residual is a a bad number (NaN or Infinity)')
               exit
            end if
            residual_norm = abs(residual_norm)
            max_residual = abs(max_residual)
            resid_norm_min = min(residual_norm, resid_norm_min)
            max_resid_min = min(max_residual, max_resid_min)
            
            if (max_correction > tol_max_correction*coeff .or. max_residual > tol_max_residual*coeff) then
               passed_tol_tests = .false.
            else
               passed_tol_tests =
     >               (correction_norm <= tol_correction_norm*coeff .and. 
     >                residual_norm <= tol_residual_norm*coeff)
     >          .or.      
     >               (abs(slope) <= tol_abs_slope_min .and. 
     >                correction_norm*residual_norm <= tol_corr_resid_product*coeff*coeff)
            end if
            
            if (.not. passed_tol_tests) then
               if (iiter >= max_tries) then
                  if (dbg_msg) then
                     call get_message
                     message = trim(message) // ' -- give up'
                     call write_msg(message)
                  end if
                  convergence_failure = .true.; exit
               else if (iwork(i_try_really_hard) == 0) then
                  if (coeff < min(min_corr_coeff,correction_factor)) then
                     call oops('coeff too small')
                     exit
                  else if (correction_norm > tol_correction_norm*coeff
     >                  .and. (correction_norm > work(r_corr_norm_jump_limit)*corr_norm_min)
     >                  .and. (.not. first_try)) then
                     call oops('avg corrrection jumped')
                     exit
                  else if (residual_norm > tol_residual_norm*coeff
     >                  .and. (residual_norm > work(r_resid_norm_jump_limit)*resid_norm_min)
     >                  .and. (.not. first_try)) then
                     call oops('avg residual jumped')
                     exit
                  else if (max_correction > tol_max_correction*coeff
     >                  .and. (max_correction > work(r_max_corr_jump_limit)*max_corr_min)
     >                  .and. (.not. first_try)) then
                     call oops('max corrrection jumped')
                     exit
                  else if (residual_norm > tol_residual_norm*coeff
     >                  .and. (max_residual > work(r_max_resid_jump_limit)*max_resid_min)
     >                  .and. (.not. first_try)) then
                     call oops('max residual jumped')
                     exit
                  else if (tiny_corr_cnt >= iwork(i_tiny_min_corr_coeff)
     >                  .and. min_corr_coeff < 1) then
                     call oops('tiny corrections')
                     exit
                  end if
               end if
            end if
            
            if (dbg_msg) then
               if (.not. passed_tol_tests) then
                  call get_message
                  call write_msg(message)
               else if (iiter < iwork(i_itermin)) then     
                  call write_msg('iiter < itermin')
               else
                  call write_msg('okay!')
               end if
            end if
            
            if (passed_tol_tests .and. (iiter+1 < max_tries)) then 
               ! about to declare victory... but may want to do another iteration
               force_iter_value = force_another_iteration(
     >                              iiter, iwork(i_itermin), lrpar, rpar, lipar, ipar)
               if (force_iter_value > 0) then
                  passed_tol_tests = .false. ! force another
                  tiny_corr_cnt = 0 ! reset the counter
                  corr_norm_min = 1d99
                  resid_norm_min = 1d99
                  max_corr_min = 1d99
                  max_resid_min = 1d99
               else if (force_iter_value < 0) then ! failure
                  call oops('force iter')
                  convergence_failure = .true.
                  exit
               end if
            end if

            iiter=iiter+1
            first_try = .false.

         end do
         

         contains
         
         
         subroutine get_message
            include 'formats'
            i = 0
            if (correction_norm > tol_correction_norm*coeff) i = i+1
            if (max_correction > tol_max_correction*coeff) i = i+2
            if (residual_norm > tol_residual_norm*coeff) i = i+4
            if (max_residual > tol_max_residual*coeff) i = i+8
            if (i == 0) then
               message = 'out of tries'
            else
               message = tol_msg(i)
            end if
         end subroutine get_message

         
         subroutine set_param_defaults
         
            if (iwork(i_itermin) == 0) iwork(i_itermin) = 2
            if (iwork(i_max_tries) == 0) iwork(i_max_tries) = 50
            if (iwork(i_tiny_min_corr_coeff) == 0) iwork(i_tiny_min_corr_coeff) = 25
            
            if (work(r_tol_residual_norm)==0) work(r_tol_residual_norm)=1d99
            if (work(r_tol_max_residual)==0) work(r_tol_max_residual)=1d99
            if (work(r_tol_max_correction)==0) work(r_tol_max_correction)=1d99
            if (work(r_target_corr_factor) == 0) work(r_target_corr_factor) = 0.9d0
            if (work(r_scale_correction_norm) == 0) work(r_scale_correction_norm) = 2d0
            if (work(r_corr_param_factor) == 0) work(r_corr_param_factor) = 10d0
            if (work(r_scale_max_correction) == 0) work(r_scale_max_correction) = 1d99
            if (work(r_corr_norm_jump_limit) == 0) work(r_corr_norm_jump_limit) = 1d99
            if (work(r_max_corr_jump_limit) == 0) work(r_max_corr_jump_limit) = 1d99
            if (work(r_resid_norm_jump_limit) == 0) work(r_resid_norm_jump_limit) = 1d99
            if (work(r_max_resid_jump_limit) == 0) work(r_max_resid_jump_limit) = 1d99
            if (work(r_min_corr_coeff) == 0) work(r_min_corr_coeff) = 1d-3
            if (work(r_slope_alert_level) == 0) work(r_slope_alert_level) = 1d0
            if (work(r_slope_crisis_level) == 0) work(r_slope_crisis_level) = 1d0
            if (work(r_tiny_corr_factor) == 0) work(r_tiny_corr_factor) = 2d0

         end subroutine set_param_defaults
         
         
         subroutine oops(msg)
            character (len=*), intent(in) :: msg
            character (len=256) :: full_msg
            full_msg = trim(msg) // ' -- give up'
            call write_msg(full_msg)
            convergence_failure = .true.
         end subroutine oops

      
         subroutine setequ(nvar, nz, x, equ, lrpar, rpar, lipar, ipar, ierr)
            integer, intent(in) :: nvar, nz
            real(dp), pointer :: x(:,:) ! (nvar, nz)
            real(dp), pointer :: equ(:,:) ! (nvar, nz)
            integer, intent(in) :: lrpar, lipar
            real(dp), intent(inout) :: rpar(:) ! (lrpar)
            integer, intent(inout) :: ipar(:) ! (lipar)
            integer, intent(out) :: ierr
            call set_primaries(nvar, nz, x, lrpar, rpar, lipar, ipar, ierr); if (ierr /= 0) return
            call set_secondaries(0, lrpar, rpar, lipar, ipar, ierr); if (ierr /= 0) return
            call eval_equations(iiter, nvar, nz, x, xscale, equ, lrpar, rpar, lipar, ipar, ierr)
            if (ierr /= 0) return
         end subroutine setequ


         subroutine adjust_correction(
     >         min_corr_coeff_in, max_corr_coeff, grad_f, f, slope, coeff, 
     >         err_msg, lrpar, rpar, lipar, ipar, ierr)
            real(dp), intent(in) :: min_corr_coeff_in
            real(dp), intent(in) :: max_corr_coeff
            real(dp), intent(in) :: grad_f(:) ! (neq) ! gradient df/dx at xold
            real(dp), intent(out) :: f ! 1/2 fvec^2. minimize this.
            real(dp), intent(in) :: slope 
            real(dp), intent(out) :: coeff 

            ! the new correction is coeff*xscale*B
            ! with min_corr_coeff <= coeff <= max_corr_coeff
            ! if all goes well, the new x will give an improvement in f
            
            character (len=*), intent(out) :: err_msg
            integer, intent(in) :: lrpar, lipar
            real(dp), intent(inout) :: rpar(:) ! (lrpar)
            integer, intent(inout) :: ipar(:) ! (lipar)
            integer, intent(out) :: ierr
      
            integer :: i, j, k, iter, k_max_corr, i_max_corr
            character (len=256) :: message
            logical :: first_time
            real(dp) :: a1, alam, alam2, alamin, a2, disc, f2, tmp1,
     >         rhs1, rhs2, temp, test, tmplam, max_corr, fold, min_corr_coeff
            real(dp) :: frac, f_target
            logical :: skip_eval_f
     
            real(dp), parameter :: alf = 1d-2 ! ensures sufficient decrease in f

            real(dp), parameter :: alam_factor = 0.2d0
            
            include 'formats'
         
            ierr = 0                  
            coeff = 0
            
            skip_eval_f = (min_corr_coeff_in == 1)
            if (skip_eval_f) then
               f = 0
            else
               do k=1,nz
                  do i=1,nvar
                     xsave(i,k) = x(i,k)
                     dxsave(i,k) = dx(i,k)
                  end do
               end do
               f = eval_f(nvar,nz,equ)
               if (is_bad(f)) then
                  ierr = -1
                  write(err_msg,*) 'adjust_correction failed in eval_f'
                  if (dbg_msg) write(*,*) 'adjust_correction: eval_f(nvar,nz,equ)', eval_f(nvar,nz,equ)
                  return
               end if
            end if
            fold = f
            
            min_corr_coeff = min(min_corr_coeff_in, max_corr_coeff) ! make sure min <= max
            alam = max_corr_coeff
            first_time = .true.
            f2 = 0
            alam2 = 0

         search_loop: do iter = 1, 1000
            
               coeff = max(min_corr_coeff, alam) 
               
               call apply_coeff(nvar, nz, x, xsave, B, xscale, coeff, skip_eval_f)
               do k=1,nz
                  do i=1,nvar
                     dx(i,k) = x(i,k) - xold(i,k)
                  end do
               end do
               call xdomain(iiter, nvar, nz, x, dx, xold, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  write(err_msg,*) 'adjust_correction failed in xdomain'
                  if (dbg_msg) write(*,*) 'adjust_correction failed in xdomain: alam', alam
                  if (alam <= min_corr_coeff) return
                  ierr = 0
                  alam = max(alam*alam_factor, min_corr_coeff)
                  cycle
               end if
               
               call setequ(nvar, nz, x, equ, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  if (alam > min_corr_coeff) then
                     alam = max(alam/10, min_corr_coeff)
                     ierr = 0
                     cycle
                  end if
                  ierr = -1
                  write(err_msg,*) 'adjust_correction failed in setequ'
                  if (dbg_msg) write(*,*) 'adjust_correction: setequ returned ierr', ierr
                  exit search_loop
               end if
               
               if (min_corr_coeff == 1) return
            
               f = eval_f(nvar,nz,equ)
               if (is_bad(f)) then
                  if (alam > min_corr_coeff) then
                     alam = max(alam/10, min_corr_coeff)
                     ierr = 0
                     cycle
                  end if
                  err_msg = 'equ norm is NaN or other bad num'
                  ierr = -1
                  exit search_loop
               end if
               
               f_target = max(fold/2, fold + alf*coeff*slope)
               if (f <= f_target) then
                  return ! sufficient decrease in f
               end if

               if (alam <= min_corr_coeff) then
                  return ! time to give up
               end if

               ! reduce alam and try again
               if (first_time) then
                  tmplam = -slope/(2*(f-fold-slope))
                  first_time = .false.
               else ! have two prior f values to work with
                  rhs1 = f - fold - alam*slope
                  rhs2 = f2 - fold - alam2*slope
                  tmp1 = 1d0/(alam2*alam2*(alam - alam2))
                  a1 = (rhs1 - rhs2)*tmp1
                  a2 = (-alam2*rhs1 + alam*rhs2)*tmp1
                  if (a1 == 0) then
                     tmplam = -slope/(2*a2)
                  else
                     disc = a2*a2-3*a1*slope
                     if (disc < 0) then
                        tmplam = alam*alam_factor
                     else if (a2 <= 0) then
                        tmplam = (-a2+sqrt(disc))/(3*a1)
                     else
                        tmplam = -slope/(a2+sqrt(disc))
                     end if
                  end if
                  if (tmplam > alam*alam_factor) tmplam = alam*alam_factor
               end if
            
               alam2 = alam
               f2 = f
               alam = max(tmplam, alam*alam_factor, min_corr_coeff)
     
            end do search_loop

            do k=1,nz
               do i=1,nvar
                  x(i,k) = xsave(i,k)
                  dx(i,k) = dxsave(i,k)
               end do
            end do
         
         end subroutine adjust_correction
         
         
         subroutine apply_coeff(nvar, nz, x, xsave, B, xscale, coeff, just_use_x)
            integer, intent(in) :: nvar, nz
            real(dp), intent(inout), dimension(:,:) :: x
            real(dp), intent(in), dimension(:,:) :: xsave, B, xscale
            real(dp), intent(in) :: coeff
            logical, intent(in) :: just_use_x
            integer :: i, k
            include 'formats'
            if (just_use_x) then
               !write(*,1) 'apply_coeff just_use_x', coeff
               if (coeff == 1d0) then
                  do k=1,nz
                     do i=1,nvar
                        !if (i==1 .and. x(i,k) + xscale(i,k)*B(i,k) < -1d-10) then
                        !   write(*,3) 'x(i,k)', i, k, x(i,k), xscale(i,k), B(i,k)
                        !   stop 'apply_coeff'
                        !end if
                        x(i,k) = x(i,k) + xscale(i,k)*B(i,k)
                     end do
                  end do
               else
                  do k=1,nz
                     do i=1,nvar
                        x(i,k) = x(i,k) + coeff*xscale(i,k)*B(i,k)
                     end do
                  end do
               end if
               return
            end if
            ! else use xsave instead of x
            if (coeff == 1d0) then
               do k=1,nz
                  do i=1,nvar
                     x(i,k) = xsave(i,k) + xscale(i,k)*B(i,k)
                  end do
               end do
               return
            end if
            do k=1,nz
               do i=1,nvar
                  x(i,k) = xsave(i,k) + coeff*xscale(i,k)*B(i,k)
               end do
            end do
         end subroutine apply_coeff


         logical function solve_equ()    
            integer ::  nrhs, ldafb, ldb, ldx, lda, i, j, n, sprs_nz
            real(dp) :: ferr, berr
            
            include 'formats'

            solve_equ=.true.
            do k=1,nz
               do i=1,nvar
                  b(i,k) = -equ(i,k)
               end do
            end do
            n = nvar*nz

            nrhs=1
            lda=mljac+1+mujac
            ldafb=2*mljac+mujac+1
            ldb=n
            ldx=n
            
            info = 0
            if (do_mtx_timing) call system_clock(time0,clock_rate)            
            call factor_mtx(n, ldafb, sprs_nz)
            if (info == 0) call solve_mtx(n, ldafb, sprs_nz)
            if (do_mtx_timing) then
               call system_clock(time1,clock_rate)
               work(r_mtx_time) = work(r_mtx_time) + dble(time1 - time0) / clock_rate
            end if

            if (info /= 0) then 
               solve_equ=.false.
               b(1:nvar,1:nz)=0
            end if
         
         end function solve_equ
         
         
         subroutine factor_mtx(n, ldafb, sprs_nz)
            integer, intent(in) :: n, ldafb
            integer, intent(out) :: sprs_nz
            integer :: i, j, k, info_dealloc
            include 'formats'
            sprs_nz = 0
            if (matrix_type == block_tridiag_dble_matrix_type) then
               if (.not. overlay_AF) then
                  do k = 1,nvar*neq
                     lblkF1(k) = lblk1(k)
                     dblkF1(k) = dblk1(k)
                     ublkF1(k) = ublk1(k)
                  end do
               end if          
               call decsolblk(
     >                  0, caller_id, nvar, nz, lblkF1, dblkF1, ublkF1, B1, ipiv_blk1,
     >                  lrd, rpar_decsol, lid, ipar_decsol, info)
               if (info /= 0) then
                  call decsolblk(
     >               2, caller_id, nvar, nz, lblkF1, dblkF1, ublkF1, B1, ipiv_blk1, 
     >               lrd, rpar_decsol, lid, ipar_decsol, info_dealloc)  
               end if
            else if (matrix_type == square_matrix_type) then
               do k = 1,n*n
                  AF1(k) = A1(k)
               end do
               call decsol(0, n, n, AF1, n, n, B1, ipiv1, 
     >               lrd, rpar_decsol, lid, ipar_decsol, info)
            else ! banded_matrix_type
               if (.not. overlay_AF) then
                  do j=1,neq
                     do i=1,ldA
                        AF(mljac+i,j) = A(i,j)
                     end do
                  end do
               end if                  
               call decsol(0, n, ldafb, AF1, mljac, mujac, B1, ipiv1, 
     >               lrd, rpar_decsol, lid, ipar_decsol, info)
            end if
         end subroutine factor_mtx
         
         
         subroutine solve_mtx(n, ldafb, sprs_nz)
            integer, intent(in) :: n, ldafb, sprs_nz
            character(1) :: trans
            integer :: info_solve, info_dealloc, i, j
            info = 0; info_solve=0; info_dealloc=0
            trans = 'N'
            if (matrix_type == block_tridiag_dble_matrix_type) then
               call decsolblk(
     >            1, caller_id, nvar, nz, lblkF1, dblkF1, ublkF1, B1, ipiv_blk1,
     >            lrd, rpar_decsol, lid, ipar_decsol, info_solve)
               call decsolblk(
     >            2, caller_id, nvar, nz, lblkF1, dblkF1, ublkF1, B1, ipiv_blk1, 
     >            lrd, rpar_decsol, lid, ipar_decsol, info_dealloc)               
            else if (matrix_type == square_matrix_type) then
               call decsol(
     >            1, n, n, AF1, n, n, B1, ipiv1, 
     >            lrd, rpar_decsol, lid, ipar_decsol, info_solve)
               call decsol(
     >            2, n, n, AF1, n, n, B1, ipiv1, 
     >            lrd, rpar_decsol, lid, ipar_decsol, info_dealloc)               
            else ! banded_matrix_type
               call decsol(
     >            1, n, ldafb, AF1, mljac, mujac, B1, ipiv1, 
     >            lrd, rpar_decsol, lid, ipar_decsol, info_solve)     
               call decsol(
     >            2, n, ldafb, AF1, mljac, mujac, B1, ipiv1, 
     >            lrd, rpar_decsol, lid, ipar_decsol, info_dealloc)
            end if
            if (info_solve /= 0 .or. info_dealloc /= 0) info = -1
         end subroutine solve_mtx
         
         
         logical function do_enter_setmatrix(
     >            neq, x, dx, xscale, lrpar, rpar, lipar, ipar, ierr)
            ! create jacobian by using numerical differences to approximate the partial derivatives
            implicit none
            integer, intent(in) :: neq
            real(dp), pointer, dimension(:,:) :: x, dx, xscale
            integer, intent(in) :: lrpar, lipar
            real(dp), intent(inout) :: rpar(:) ! (lrpar)
            integer, intent(inout) :: ipar(:) ! (lipar)
            integer, intent(out) :: ierr
            logical :: need_solver_to_eval_jacobian
            integer :: i, j, k
            include 'formats'
            need_solver_to_eval_jacobian = .true.
            call enter_setmatrix(iiter, 
     >                  nvar, nz, neq, x, xold, xscale, xder, need_solver_to_eval_jacobian, 
     >                  size(A,dim=1), A1, idiag, lrpar, rpar, lipar, ipar, ierr)
            do_enter_setmatrix = need_solver_to_eval_jacobian
         end function do_enter_setmatrix


         subroutine setmatrix(neq, x, dx, xscale, xsave, dxsave, lrpar, rpar, lipar, ipar, ierr)
            ! create jacobian by using numerical differences to approximate the partial derivatives
            implicit none
            integer, intent(in) :: neq
            real(dp), pointer, dimension(:,:) :: x, dx, xscale, xsave, dxsave
            integer, intent(in) :: lrpar, lipar
            real(dp), intent(inout) :: rpar(:) ! (lrpar)
            integer, intent(inout) :: ipar(:) ! (lipar)
            integer, intent(out) :: ierr

            integer :: i, j, ii, jj, k, kk, ij, ik, ivar, jvar, iz, jz, jzz, ideb, ifin
            integer, dimension(nvar) :: nskip, gskip, dskip
            real(dp) :: dscale, partial
            logical :: need_solver_to_eval_jacobian
            
            include 'formats'

            ierr = 0
            
            need_solver_to_eval_jacobian = do_enter_setmatrix(
     >            neq, x, dx, xscale, lrpar, rpar, lipar, ipar, ierr)     
            if (ierr /= 0) return
            if (.not. need_solver_to_eval_jacobian) return
            
            if (matrix_type == block_tridiag_dble_matrix_type) then
               write(*,'(a)') 'sorry: newton numerical jacobian does ' //
     >               'not support numerical block triangular jacobians.'
               write(*,*) 'requested matrix_type', matrix_type
               write(*,*) 'try using a banded matrix instead'
               ierr = -1
               return
            end if
            
            ! allocate working arrays for numerical jacobian calculation
            allocate(xgg(nvar,nz), dxd(nvar,nz), dxdd(nvar,nz), equsave(nvar,nz), stat=ierr)

            do k=1,nz
               do j=1,nvar
                  xsave(j,k) = x(j,k)
                  dxsave(j,k) = dx(j,k)
                  equsave(j,k) = equ(j,k)
               end do
            end do
            if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y1=y      

            ! some info about the stencil
            ! gskip zones on left
            ! dskip zones on right
            ! nskip zones in total
            gskip=mljac/nvar
            dskip=mujac/nvar
            nskip=1+dskip+gskip

            A=0
            ! loop on variables
            do ivar=1, nvar
               do k=1,nz
                  do j=1,nvar
                     dxd(j,k) = dxsave(j,k)
                  end do
               end do
               do k=1, nz
                  do ii=1, 20 ! may need to increase xder
                     dxd(ivar,k)=dxd(ivar,k)+xder(ivar,k)
                     if (dxd(ivar,k)-dxsave(ivar,k) /= 0) exit
                     xder(ivar,k)=xder(ivar,k)*2
                  end do
               end do
               do k=1,nz
                  do j=1,nvar
                     dx(j,k) = dxd(j,k)
                     x(j,k) = xold(j,k) + dx(j,k)
                     dx(j,k) = x(j,k) - xold(j,k)
                  end do
               end do
               call xdomain(iiter, nvar, nz, x, dx, xold, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y = y1 
                  call cleanup_after_setmatrix
                  call failed_in_setmatrix(0, lrpar, rpar, lipar, ipar, ierr)
                  return
               end if
               do k=1,nz
                  do j=1,nvar
                     dxd(j,k) = dx(j,k)
                  end do
               end do
               call set_primaries(nvar, nz, x, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y = y1
                  call cleanup_after_setmatrix
                  call failed_in_setmatrix(0, lrpar, rpar, lipar, ipar, ierr)
                  return
               end if
               ! compute secondary variables for modified primaries
               call set_secondaries(ivar, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y = y1
                  call cleanup_after_setmatrix
                  call failed_in_setmatrix(0, lrpar, rpar, lipar, ipar, ierr)
                  return
               end if
               if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y2 = y

               ! now use the modified primaries and secondaries to get modified equations
               do kk=0, nskip(ivar)-1

                  if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y = y1
                  do k=1,nz
                     do j=1,nvar
                        x(j,k) = dxsave(j,k) 
                     end do
                  end do
                  ! primaries are changed only on the zones of the comb
                  do k = 1+kk, nz, nskip(ivar)
                     dx(ivar,k) = dxd(ivar,k)
                  end do
                  do k=1,nz
                     do j=1,nvar
                        x(j,k) = xold(j,k) + dx(j,k)
                        dxdd(j,k) = dx(j,k)
                     end do
                  end do
                  call set_primaries(nvar, nz, x, lrpar, rpar, lipar, ipar, ierr)
                  if (ierr /= 0) then
                     if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y = y1
                     call cleanup_after_setmatrix
                     call failed_in_setmatrix(0, lrpar, rpar, lipar, ipar, ierr)
                     return
                  end if
                  
                  if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  then
                     ! note that we can use the previously computed secondaries
                     ! since, by definition, they depend only on the primaries of their own zone.
                     !do j=1+kk, nz, nskip(ivar)
                     !   y(j, 1:nsec)=y2(j, 1:nsec)
                     !enddo
                  !end if
                  
                  ! compute the equations using these primaries and secondaries
                  call eval_equations(iiter, nvar, nz, x, xscale, equ, lrpar, rpar, lipar, ipar, ierr)
                  if (ierr /= 0) then
                     if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y = y1
                     call cleanup_after_setmatrix
                     call failed_in_setmatrix(0, lrpar, rpar, lipar, ipar, ierr)
                     return
                  end if

                  ! compute derivatives
                  do j = ivar+kk*nvar, neq, nvar*nskip(ivar)
                     zone = (j-1)/nvar + 1
                     if (dxdd(ivar,zone) == dxsave(ivar,zone)) then 
                        ! can happen if the xdomain routine changed dx in a bad way.
                        ierr = -1
                        write(*, '(a, i5, 99e20.10)') 
     >                     'failed trying to create numerical derivative for variable ',
     >                     j, dxsave(ivar,zone), xsave(ivar,zone), xder(ivar,zone)
                        if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y = y1 
                        call cleanup_after_setmatrix
                        call failed_in_setmatrix(j, lrpar, rpar, lipar, ipar, ierr)
                        return
                     endif
                     ideb=max(1, (zone-gskip(ivar)-1)*nvar+1)
                     ifin=min(neq, (zone+dskip(ivar))*nvar)
                     ideb=max(ideb, j-mljac)
                     ifin=min(ifin, j+mujac)
                     do i = ideb, ifin
                        ik = (i-1)/nvar + 1
                        ij = i - (ik-1)*nvar
                        partial=xscale(ivar,zone)*
     >                     (equ(ij,ik)-equsave(ij,ik))/(dxdd(ivar,zone)-dxsave(ivar,zone))
                        if (matrix_type == square_matrix_type) then
                           A(i,j)=partial
                        else
                           A(i-j+idiag,j)=partial
                        end if
                     end do
                  end do

                  if (nsec > 0) then ! restore the secondaries that correspond to the unmodified primaries
                     !do j=1+kk, nz, nskip(ivar)       
                     !   y(j, 1:nsec)=y1(j, 1:nsec)
                     !enddo
                  end if

               enddo
         
            enddo

            if (nsec > 0) call mesa_error(__FILE__,__LINE__) !  y = y1 
            call cleanup_after_setmatrix

            call exit_setmatrix(iiter, nvar, nz, neq, 
     >            dx, ldA, A1, idiag, xscale, lrpar, rpar, lipar, ipar, ierr)

         end subroutine setmatrix
         
         
         subroutine cleanup_after_setmatrix
            integer :: i, k
            do k=1,nz
               do i=1,nvar
                  x(i,k) = xsave(i,k)
                  dx(i,k) = dxsave(i,k)
                  equ(i,k) = equsave(i,k)
               end do
            end do
            deallocate(xgg, dxd, dxdd, equsave)
         end subroutine cleanup_after_setmatrix
         
      
         subroutine write_msg(msg)
            real(dp), parameter :: secyer = 3.1558149984d7 ! seconds per year
            character(*)  :: msg
            if (.not. dbg_msg) return
            
    1       format(i6, 2x, i3, 2x, a, f8.4, 6(2x, a, 1x, e10.3), 2x, a, f6.2, 2x, a)            
            write(*,1)
     >         iwork(i_model_number), iiter,
     >         'coeff', coeff, 
     >         'slope', slope, 
     >         'f', f,
     >         'avg resid', residual_norm, 
     >         'max resid', max_residual, 
     >         'avg corr', correction_norm, 
     >         'max corr', max_correction, 
     >         'lg dt/yr', log10(max(1d-99,work(r_dt)/secyer)), 
     >         trim(msg)            
         end subroutine write_msg
      
      
         subroutine newton_core_dump(x, dx, xold)
            real(dp), dimension(:,:) :: x 
            ! new vector of primaries, x = xold+dx
            real(dp), dimension(:,:) :: dx 
            ! increment vector from previous vector of primaries.
            real(dp), dimension(:,:) :: xold 
            ! xold = x-dx.  xold is kept constant; x and dx change.
            integer :: i, j, k
         
    1       format(a20, i16) ! integers
    2       format(a20, 1pe26.16) ! reals
    3       format(a20, i6, 1x, 1pe26.16) ! 1 index reals
    4       format(a20, 2(i6, 1x), 1pe26.16) ! 2 index reals
    5       format(a20, i6, 1x, i16) ! 1 index integers
         
            ! only printout args and things that are carried over from one call to next
            ! e.g., skip work arrays that are written on each call before they are read
         
            write(*, *) 'newton core dump'
            write(*, 1) 'nz', nz
            write(*, 1) 'nvar', nvar
            write(*, 1) 'mljac', mljac
            write(*, 1) 'mujac', mujac
            write(*, 1) 'liwork', liwork
            write(*, 1) 'lwork', lwork
            write(*, 1) 'ldy', ldy
            write(*, 1) 'nsec', nsec
            write(*, 1) 'ldAF', ldAF
            write(*, 1) 'ndiag', ndiag

            write(*, 2) 'tol_correction_norm', tol_correction_norm
         
            do j=1, ndiag
               do k=1, nz
                  write(*, 4) 'A', j, k, A(j, k)
               end do
            end do
         
            do j=1, ldAF
               do k=1, nz
                  write(*, 4) 'AF', j, k, AF(j, k)
               end do
            end do
         
            do k=1, nz
               write(*, 5) 'ipiv1', k, ipiv1(k)
            end do
         
            do k=1, nz
               do j=1, nvar
                  write(*, 4) 'x', j, k, x(j, k)
                  write(*, 4) 'dx', j, k, x(j, k)
                  write(*, 4) 'xold', j, k, x(j, k)
               end do
            end do
         
         end subroutine newton_core_dump


         subroutine pointers(ierr)
            integer, intent(out) :: ierr
      
            integer :: i, j
            character (len=256) :: err_msg

            ierr = 0         
            i = num_work_params+1
            
            A1(1:ndiag*neq) => work(i:i+ndiag*neq-1); i = i+ndiag*neq
            A(1:ndiag,1:neq) => A1(1:ndiag*neq)
            Acopy1 => A1
            Acopy => A
            
            xsave1(1:neq) => work(i:i+neq-1); i = i+neq
            xsave(1:nvar,1:nz) => xsave1(1:neq)
            
            dxsave1(1:neq) => work(i:i+neq-1); i = i+neq
            dxsave(1:nvar,1:nz) => dxsave1(1:neq)
            
            B1 => work(i:i+neq-1); i = i+neq
            B(1:nvar,1:nz) => B1(1:neq)
            
            B_init1 => work(i:i+neq-1); i = i+neq
            B_init(1:nvar,1:nz) => B_init1(1:neq)
            
            grad_f1(1:neq) => work(i:i+neq-1); i = i+neq
            grad_f(1:nvar,1:nz) => grad_f1(1:neq)
            
            rhs(1:nvar,1:nz) => work(i:i+neq-1); i = i+neq
            
            xder(1:nvar,1:nz) => work(i:i+neq-1); i = i+neq
            
            dx(1:nvar,1:nz) => work(i:i+neq-1); i = i+neq
            
            if (nsec > 0) then
               !y1(1:nvar,1:nz) => work(i:i+nsec*neq-1); i = i+nsec*neq
               !y2(1:nvar,1:nz) => work(i:i+nsec*neq-1); i = i+nsec*neq
            else
               nullify(y1)
               nullify(y2)
            end if

            if (i-1 > lwork) then
               ierr = -1
               write(*, 
     >                  '(a, i6, a, 99i6)') 'newton: lwork is too small.  must be at least', i-1,
     >                  '   but is only ', lwork, neq, ndiag, ldAF, nsec
               return
            end if
         
            i = num_iwork_params+1
            ipiv1(1:neq) => iwork(i:i+neq-1); i = i+neq
            if (i-1 > liwork) then
               ierr = -1
               write(*, '(a, i6, a, i6)') 
     >                  'newton: liwork is too small.  must be at least', i, 
     >                  '   but is only ', liwork
               return
            end if
            
            if (matrix_type == block_tridiag_dble_matrix_type) then
     
               ipiv_blk1(1:neq) => ipiv1(1:neq)
               ublk1(1:nvar*neq) => A1(1:nvar*neq)
               dblk1(1:nvar*neq) => A1(1+nvar*neq:2*nvar*neq)
               lblk1(1:nvar*neq) => A1(1+2*nvar*neq:3*nvar*neq)
               lblk(1:nvar,1:nvar,1:nz) => lblk1(1:nvar*neq)
               dblk(1:nvar,1:nvar,1:nz) => dblk1(1:nvar*neq)
               ublk(1:nvar,1:nvar,1:nz) => ublk1(1:nvar*neq)
               
               ! testing
               k = 2*nvar*neq - nvar*nvar
               do i=1,nvar
                  do j=1,nvar
                     dblk(i,j,nz) = i*100+j
                     if (dblk(i,j,nz) /= A1(k + nvar*(j-1) + i)) then
                        call mesa_error(__FILE__,__LINE__)
                     end if
                     !if (dblk1(nvar*nvar*(nz-1)+j) /= dblk(i,j,nz)) then
                     !   call mesa_error(__FILE__,__LINE__)
                     !end if
                  end do
               end do

            end if
               
            if (matrix_type == block_tridiag_dble_matrix_type) then
            
               ublkF1(1:nvar*neq) => AF1(1:nvar*neq)
               dblkF1(1:nvar*neq) => AF1(1+nvar*neq:2*nvar*neq)
               lblkF1(1:nvar*neq) => AF1(1+2*nvar*neq:3*nvar*neq)

               lblkF(1:nvar,1:nvar,1:nz) => lblkF1(1:nvar*neq)
               dblkF(1:nvar,1:nvar,1:nz) => dblkF1(1:nvar*neq)
               ublkF(1:nvar,1:nvar,1:nz) => ublkF1(1:nvar*neq)
               
            end if
         
         end subroutine pointers
         
         
         real(dp) function eval_slope(nvar, nz, grad_f, B)
            integer, intent(in) :: nvar, nz
            real(dp), intent(in), dimension(:,:) :: grad_f, B
            integer :: k, i
            eval_slope = 0
            do i=1,nvar
               eval_slope = eval_slope + dot_product(grad_f(i,1:nz),B(i,1:nz))
            end do
         end function eval_slope
         
         
         real(dp) function eval_f(nvar, nz, equ)
            integer, intent(in) :: nvar, nz
            real(dp), intent(in), dimension(:,:) :: equ
            integer :: k, i
            real*8 :: q
            include 'formats'
            eval_f = 0
            do k = 1, nz
               do i = 1, nvar
                  q = equ(i,k)
                  eval_f = eval_f + q*q
               end do
            end do
            eval_f = eval_f/2
            !write(*,1) 'do_newton: eval_f', eval_f
         end function eval_f


      end subroutine do_newton
      
   
      subroutine get_newton_work_sizes(
     >      mljac, mujac, nvar, nz, nsec, matrix_type, lwork, liwork, ierr)
         integer, intent(in) :: mljac, mujac, nvar, nz, nsec
         integer, intent(in) :: matrix_type
         integer, intent(out) :: lwork, liwork
         integer, intent(out) :: ierr
         
         integer :: i, ndiag, ldAF, neq
         
         include 'formats'

         ierr = 0
         neq = nvar*nz
         
         if (matrix_type == square_matrix_type) then
            ndiag = neq
            ldAF = ndiag
         else if (matrix_type == block_tridiag_dble_matrix_type) then
            ndiag = 3*nvar
            ldAF = ndiag
         else
            ndiag = mljac+mujac+1
            ldAF = mljac+ndiag
         end if
         
         liwork = num_iwork_params + neq     
         lwork = num_work_params + neq*(ndiag + 9 + 2*nsec)
         
      end subroutine get_newton_work_sizes


      end module mod_newton
