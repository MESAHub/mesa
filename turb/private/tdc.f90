! ***********************************************************************
!
!   Copyright (C) 2010-2021  The MESA Team
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


module tdc

use const_def
use num_lib
use utils_lib
use auto_diff
use star_data_def
use tdc_support

implicit none

private
public :: get_TDC_solution

contains


   !> Computes the outputs of time-dependent convection theory following the model specified in
   !! Radek Smolec's thesis [https://users.camk.edu.pl/smolec/phd_smolec.pdf], which in turn
   !! follows the model of Kuhfuss 1986.
   !!
   !! Internally this solves the equation L = L_conv + L_rad.
   !!
   !! @param info tdc_info type storing various quantities that are independent of Y.
   !! @param scale The scale for computing residuals to the luminosity equation (erg/s).
   !! @param Zlb Lower bound on Z.
   !! @param Zub Upper bound on Z.
   !! @param conv_vel The convection speed (cm/s).
   !! @param Y_face The superadiabaticity (dlnT/dlnP - grada, output).
   !! @param tdc_num_iters Number of iterations taken in the TDC solver.
   !! @param ierr Tracks errors (output).
   subroutine get_TDC_solution(info, scale, Zlb, Zub, conv_vel, Y_face, tdc_num_iters, ierr)
      type(tdc_info), intent(in) :: info
      real(dp), intent(in) :: scale
      type(auto_diff_real_tdc), intent(in) :: Zlb, Zub
      type(auto_diff_real_star_order1),intent(out) :: conv_vel, Y_face
      integer, intent(out) :: tdc_num_iters, ierr
      
      logical :: Y_is_positive
      type(auto_diff_real_tdc) :: Af, Y, Y0, Y1, Z0, Z1, radY
      type(auto_diff_real_tdc) :: Q, Q0
      logical :: has_root
      integer :: iter
      integer, parameter :: max_iter = 30
      include 'formats'

      ierr = 0
      if (info%mixing_length_alpha == 0d0 .or. info%dt <= 0d0) then
         call mesa_error(__FILE__,__LINE__,'bad call to TDC get_TDC_solution')
      end if         

      ! Determine the sign of the solution.
      !
      ! If Q(Y=0) is positive then the luminosity is too great to be carried radiatively, so
      ! we'll necessarily have Y > 0.
      !
      ! If Q(Y=0) is negative then the luminosity can be carried by radiation alone, so we'll
      ! necessarily have Y < 0.
      Y = 0d0
      call compute_Q(info, Y, Q0, Af)
      if (Q0 > 0d0) then
         Y_is_positive = .true.
      else
         Y_is_positive = .false.
      end if

      if (info%report) then
         open(unit=4,file='out.data')
         do iter=1,1000
            Z0 = (Zlb + (Zub-Zlb)*(iter-1)/1000)
            Y = set_Y(Y_is_positive,Z0)
            call compute_Q(info, Y, Q, Af)
            write(4,*) Y%val, Q%val
         end do
         write(*,*) 'Wrote Q(Y) to out.data'
      end if

      ! Start down the chain of logic...
      if (Y_is_positive) then
         ! If Y > 0 then Q(Y) is monotone and we can jump straight to the search.
         call bracket_plus_Newton_search(info, scale, Y_is_positive, Zlb, Zub, Y_face, Af, tdc_num_iters, ierr)
         Y = convert(Y_face)
         if (ierr /= 0) return
         if (info%report) write(*,*) 'Y is positive, Y=',Y_face%val
      else
         if (info%report) write(*,*) 'Y is negative.'
         ! If Y < 0 then Q(Y) is not guaranteed to be monotone, so we have to be more careful.
         ! One root we could have is the radiative solution (with Af==0), given by
         radY = (info%L - info%L0 * info%gradL) / info%L0

         ! If A0 == 0 then, because Af(Y) is monotone-increasing in Y, we know that for Y < 0, Af(Y) = 0.
         ! As a result we can directly write down the solution and get just radY.
         if (info%A0 == 0) then
            Y = radY
            if (info%report) write(*,*) 'A0 == 0, Y=',Y%val
         else
            if (info%report) write(*,*) 'A0 > 0.'
            ! Otherwise, we keep going.
            ! We next identify the point where Af(Y) = 0. Call this Y0, corresponding to Z0.
            call Af_bisection_search(info, Zlb, Zub, Z0, Af, ierr)
            if (ierr /= 0) return
            Y0 = set_Y(.false., Z0)
            call compute_Q(info, Y0, Q, Af)
            if (info%report) write(*,*) 'Bisected Af. Y0=',Y0%val,'Af(Y0)=',Af%val

            ! We next need to do a bracket search for where dQdZ = 0 over the interval [Y0,0] (equivalently from Z=lower_bound to Z=Z0).
            call dQdZ_bisection_search(info, Zlb, Z0, Z1, has_root)
            if (has_root) then
               Y1 = set_Y(.false., Z1)
               if (info%report) write(*,*) 'Bisected dQdZ, found root, ',Y1%val
               call compute_Q(info, Y1, Q, Af)
               if (Q < 0) then ! Means there are no roots with Af > 0.
                  if (info%report) write(*,*) 'Root has Q<0, Q=',Q%val,'Y=',radY%val
                  Y = radY
               else
                  if (info%report) write(*,*) 'Root has Q>0. Q(Y1)=',Q%val
                  ! Do a search over [lower_bound, Z1]. If we find a root, that's the root closest to zero so call it done.
                  if (info%report) write(*,*) 'Searching from Y=',-exp(Zlb%val),'to Y=',-exp(Z1%val)
                  call bracket_plus_Newton_search(info, scale, Y_is_positive, Zlb, Z1, Y_face, Af, tdc_num_iters, ierr)
                  Y = convert(Y_face)
                  if (info%report) write(*,*) 'ierr',ierr, tdc_num_iters
                  if (ierr /= 0) then
                     if (info%report) write(*,*) 'No root found. Searching from Y=',-exp(Z1%val),'to Y=',-exp(Z0%val)
                     ! Do a search over [Z1, Z0]. If we find a root, that's the root closest to zero so call it done.
                     ! Note that if we get to this stage there is (mathematically) guaranteed to be a root, modulo precision issues.
                     call bracket_plus_Newton_search(info, scale, Y_is_positive, Z1, Z0, Y_face, Af, tdc_num_iters, ierr)
                     Y = convert(Y_face)
                  end if
                  if (info%report) write(*,*) 'Y=',Y%val
               end if
            else
               if (info%report) write(*,*) 'Bisected dQdZ, no root found.'
               call compute_Q(info, Y0, Q, Af)
               if (Q > 0) then ! Means there's a root in [Y0,0] so we bracket search from [lower_bound,Z0]
                  call bracket_plus_Newton_search(info, scale, Y_is_positive, Zlb, Z0, Y_face, Af, tdc_num_iters, ierr)
                  Y = convert(Y_face)
                  if (info%report) write(*,*) 'Q(Y0) > 0, bisected and found Y=',Y%val
               else ! Means there's no root in [Y0,0] so the only root is radY.
                  if (info%report) write(*,*) 'Q(Y0) < 0, Y=',radY%val
                  Y = radY
               end if
            end if
         end if
      end if


      ! Process Y into the various outputs.
      call compute_Q(info, Y, Q, Af)
      Y_face = unconvert(Y)
      conv_vel = sqrt_2_div_3*unconvert(Af)   

   end subroutine get_TDC_solution

   !> Performs a bracket search for the solution to Q=0 over
   !! a domain in which Q is guaranteed to be monotone. Then
   !! refines the result with a Newton solver to both get a better
   !! solution and endow the solution with partial derivatives.
   !!
   !! @param info tdc_info type storing various quantities that are independent of Y.
   !! @param scale The scale for computing residuals to the luminosity equation (erg/s).
   !! @param Zlb Lower bound on Z.
   !! @param Zub Upper bound on Z.
   !! @param Y_face The superadiabaticity (dlnT/dlnP - grada, output).
   !! @param tdc_num_iters Number of iterations taken in the TDC solver.
   !! @param ierr Tracks errors (output).
   subroutine bracket_plus_Newton_search(info, scale, Y_is_positive, Zlb, Zub, Y_face, Af, tdc_num_iters, ierr)
      type(tdc_info), intent(in) :: info
      logical, intent(in) :: Y_is_positive
      real(dp), intent(in) :: scale
      type(auto_diff_real_tdc), intent(in) :: Zlb, Zub
      type(auto_diff_real_star_order1),intent(out) :: Y_face
      type(auto_diff_real_tdc), intent(out) :: Af
      integer, intent(out) :: tdc_num_iters
      integer, intent(out) :: ierr
      
      type(auto_diff_real_tdc) :: Y, Z, Q, Q_lb, Q_ub, Qc, Z_new, correction, lower_bound_Z, upper_bound_Z
      type(auto_diff_real_tdc) :: dQdZ, Q0
      integer :: iter, line_iter, i
      logical :: converged, have_derivatives, corr_has_derivatives
      real(dp), parameter :: correction_tolerance = 1d-13
      real(dp), parameter :: residual_tolerance = 1d-8
      real(dp), parameter :: alpha_c = (1d0/2d0)*sqrt_2_div_3
      integer, parameter :: max_iter = 200
      integer, parameter :: max_line_search_iter = 5
      include 'formats'

      ierr = 0

      ! We start by bisecting to find a narrow interval around the root.
      lower_bound_Z = Zlb
      upper_bound_Z = Zub

      ! Perform bisection search.
      call Q_bisection_search(info, Y_is_positive, lower_bound_Z, upper_bound_Z, Z, ierr)
      if (ierr /= 0) return

      ! Set up Z from bisection search
      Z%d1val1 = 1d0 ! Set derivative dZ/dZ=1 for Newton iterations.
      if (info%report) write(*,*) 'Z from bisection search', Z%val
      if (info%report) write(*,*) 'lower_bound_Z, upper_bound_Z',lower_bound_Z%val,upper_bound_Z%val

      ! Now we refine the solution with a Newton solve.
      ! This also let's us pick up the derivative of the solution with respect to input parameters.

      ! Initialize starting values for TDC Newton iterations.
      dQdz = 0d0
      converged = .false.
      have_derivatives = .false. ! Tracks if we've done at least one Newton iteration.
                                 ! Need to do this before returning to endow Y with partials
                                 ! with respect to the structure variables.
      do iter = 1, max_iter
         Y = set_Y(Y_is_positive, Z)
         call compute_Q(info, Y, Q, Af)

         if (abs(Q%val)/scale <= residual_tolerance .and. have_derivatives) then
            ! Can't exit on the first iteration, otherwise we have no derivative information.
            if (info%report) write(*,2) 'converged', iter, abs(Q%val)/scale, residual_tolerance
            converged = .true.
            exit
         end if

         ! We use the fact that Q(Y) is monotonic to iteratively refined bounds on Q.
         dQdZ = differentiate_1(Q)
         if (Q > 0d0 .and. dQdZ < 0d0) then
            lower_bound_Z = Z
         else if (Q > 0d0 .and. dQdZ > 0d0) then
            upper_bound_Z = Z
         else if (Q < 0d0 .and. dQdZ < 0d0) then
            upper_bound_Z = Z
         else
            lower_bound_Z = Z
         end if

         if (is_bad(dQdZ%val) .or. abs(dQdZ%val) < 1d-99) then
            ierr = 1
            exit
         end if

         correction = -Q/dQdz
         corr_has_derivatives = .true.

         ! Do a line search to avoid steps that are too big.
         do line_iter=1,max_line_search_iter

            if (abs(correction) < correction_tolerance .and. have_derivatives) then
               ! Can't get much more precision than this.
               converged = .true.
               exit
            end if

            Z_new = Z + correction
            if (corr_has_derivatives) then
               have_derivatives = .true.
            end if

            ! If the correction pushes the solution out of bounds then we know
            ! that was a bad step. Bad steps are still in the same direction, they just
            ! go too far, so we replace that result with one that's halfway to the relevant bound.
            if (Z_new > upper_bound_Z) then
               Z_new = (Z + upper_bound_Z) / 2d0
            else if (Z_new < lower_bound_Z) then
               Z_new = (Z + lower_bound_Z) / 2d0
            end if

            Y = set_Y(Y_is_positive,Z_new)

            call compute_Q(info, Y, Qc, Af)

            if (abs(Qc) < abs(Q)) then
               exit
            else
               correction = 0.5d0 * correction
            end if
         end do

         if (info%report) write(*,3) 'i, li, Z_new, Z, low_bnd, upr_bnd, Q, dQdZ, corr', iter, line_iter, &
            Z_new%val, Z%val, lower_bound_Z%val, upper_bound_Z%val, Q%val, dQdZ%val, correction%val
         Z_new%d1val1 = 1d0 ! Ensures that dZ/dZ = 1.
         Z = Z_new

         Y = set_Y(Y_is_positive,Z)
         if (converged) exit

      end do

      if (.not. converged) then
         ierr = 1
         if (info%report) then
         !$OMP critical (tdc_crit0)
            write(*,*) 'failed get_TDC_solution TDC_iter', &
               iter
            write(*,*) 'scale', scale
            write(*,*) 'Q/scale', Q%val/scale
            write(*,*) 'tolerance', residual_tolerance
            write(*,*) 'dQdZ', dQdZ%val
            write(*,*) 'Y', Y%val
            write(*,*) 'dYdZ', Y%d1val1
            write(*,*) 'exp(Z)', exp(Z%val)
            write(*,*) 'Z', Z%val
            write(*,*) 'Af', Af%val
            write(*,*) 'A0', info%A0%val
            write(*,*) 'c0', info%c0%val
            write(*,*) 'L', info%L%val
            write(*,*) 'L0', info%L0%val
            write(*,*) 'grada', info%grada%val
            write(*,*) 'gradL', info%gradL%val
            write(*,'(A)')
         !$OMP end critical (tdc_crit0)
         end if
         return
      end if

      ! Unpack output
      Y_face = unconvert(Y)
      tdc_num_iters = iter          
   end subroutine bracket_plus_Newton_search

end module tdc
