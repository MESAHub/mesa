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


module tdc_support

use const_def
use num_lib
use utils_lib
use auto_diff
use star_data_def

implicit none

private
public :: set_Y, Q_bisection_search, dQdZ_bisection_search, Af_bisection_search, convert, unconvert, safe_atan, safe_tanh

contains

   !> Stores the information which is required to evaluate TDC-related quantities and which
   !! do not depend on Y.
   !!
   !! @param mixing_length_alpha Mixing length parameter
   !! @param alpha_TDC_DAMP TDC turbulent damping parameter
   !! @param alpha_TDC_DAMPR TDC radiative damping parameter
   !! @param alpha_TDC_PtdVdt TDC coefficient on P_turb*dV/dt. Physically should probably be 1.
   !! @param dt Time-step
   !! @param c0 A proportionality factor for the convective luminosity
   !! @param L luminosity
   !! @param L0 L0 = (Lrad / grad_rad) is the luminosity radiation would carry if dlnT/dlnP = 1.
   !! @param A0 Initial convection speed
   !! @param T Temperature
   !! @param rho Density (g/cm^3)
   !! @param dV 1/rho_face - 1/rho_start_face (change in specific volume at the face)
   !! @param Cp Heat capacity
   !! @param kap Opacity
   !! @param Hp Pressure scale height
   !! @param gradL gradL is the neutrally buoyant dlnT/dlnP (= grad_ad + grad_mu),
   !! @param grada grada is the adiabatic dlnT/dlnP,
   type tdc_info
      real(dp) :: mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt
      type(auto_diff_real_tdc) :: A0, c0, L, L0, gradL, grada
      type(auto_diff_real_star_order1) :: T, rho, dV, Cp, kap, Hp
   end type tdc_info

   !> Y = +- exp(Z)
   !! If Y > 0, Y = exp(Z)
   !! If Y < 0, Y = -exp(Z)
   !!
   !! @param Y_is_positive True if Y > 0, False otherwise.
   !! @param Z log|Y|
   !! @param Y (output)
   type(auto_diff_real_tdc) function set_Y(Y_is_positive, Z) result(Y)
      logical, intent(in) :: Y_is_positive
      type(auto_diff_real_tdc), intent(in) :: Z
      if (Y_is_positive) then
         Y = exp(Z)
      else
         Y = -exp(Z)
      end if
   end function set_Y

   !> This routine performs a bisection search for Q=0 over a domain in Z for which Q is monotone.
   subroutine Q_bisection_search(&
         mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
         Y_is_positive, lower_bound_Z, upper_bound_Z, &
         c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, ierr)
      ! Inputs
      real(dp), intent(in) :: mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt
      logical, intent(in) :: Y_is_positive
      type(auto_diff_real_tdc), intent(in) :: A0, c0, L, L0, gradL, grada
      type(auto_diff_real_star_order1), intent(in) :: T, rho, dV, Cp, kap, Hp

      ! In/Out
      type(auto_diff_real_tdc), intent(inout) :: lower_bound_Z, upper_bound_Z

      ! Outputs
      integer, intent(out) :: ierr

      real(dp), parameter :: bracket_tolerance = 1d0
      integer, parameter :: max_iter = 30

      type(auto_diff_real_tdc) :: Z_new, Y, Af, Q, Q_ub, Q_lb
      integer :: iter

      ierr = 0

      Y = set_Y(Y_is_positive, lower_bound_Z)
      call compute_Q(mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
         Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q_lb, Af)

      Y = set_Y(Y_is_positive, upper_bound_Z)
      call compute_Q(mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
         Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q_ub, Af)

      ! Check to make sure that the lower and upper bounds on Z actually bracket
      ! a solution to Q(Y(Z)) = 0.
      if (Q_lb * Q_ub > 0d0) then
         if (report) then
            write(*,*) 'TDC Error. Initial Z window does not bracket a solution.'
            write(*,*) 'Q(Lower Z)',Q_lb%val
            write(*,*) 'Q(Upper Z)',Q_ub%val
            write(*,*) 'scale', scale
            write(*,*) 'tolerance', residual_tolerance
            write(*,*) 'Y', Y%val
            write(*,*) 'dYdZ', Y%d1val1
            write(*,*) 'exp(Z)', exp(Z%val)
            write(*,*) 'Z', Z%val
            write(*,*) 'A0', A0%val
            write(*,*) 'c0', c0%val
            write(*,*) 'L', L%val
            write(*,*) 'L0', L0%val
            write(*,*) 'grada', grada%val
            write(*,*) 'gradL', gradL%val
            write(*,'(A)')
         end if
         ierr = 1
         return
      end if

      do iter=1,max_iter
         Z_new = (upper_bound_Z + lower_bound_Z) / 2d0
         Y = set_Y(Y_is_positive, Z_new)

         call compute_Q(mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
            Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q, Af)

         if (Q > 0d0 .and. Q_ub > 0d0) then
            upper_bound_Z = Z_new
            Q_ub = Q
         else if (Q > 0d0 .and. Q_lb > 0d0) then
            lower_bound_Z = Z_new
            Q_lb = Q
         else if (Q < 0d0 .and. Q_ub < 0d0) then
            upper_bound_Z = Z_new
            Q_ub = Q
         else if (Q < 0d0 .and. Q_lb < 0d0) then
            lower_bound_Z = Z_new
            Q_lb = Q
         end if

         if (upper_bound_Z - lower_bound_Z < bracket_tolerance) return
      end do

   end subroutine Q_bisection_search

   !> This routine performs a bisection search for dQ/dZ=0 with Y < 0.
   !! The domain is assumed to be restricted to have Af > 0, so that dQ/dZ is
   !! continuous and monotone.
   subroutine dQdZ_bisection_search(&
         mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
         lower_bound_Z, upper_bound_Z, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, has_root)
      ! Inputs
      real(dp), intent(in) :: mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt
      type(auto_diff_real_tdc), intent(in) :: A0, c0, L, L0, gradL, grada
      type(auto_diff_real_star_order1), intent(in) :: T, rho, dV, Cp, kap, Hp

      ! In/Out
      type(auto_diff_real_tdc), intent(inout) :: lower_bound_Z, upper_bound_Z

      ! Outputs
      logical, intent(out) :: has_root

      real(dp), parameter :: bracket_tolerance = 1d-3
      integer, parameter :: max_iter = 30

      type(auto_diff_real_tdc) :: Z_new, Y, Af, Q, dQdZ, dQdZ_lb, dQdZ_ub
      integer :: iter

      ierr = 0

      Y = set_Y(Y_is_positive, lower_bound_Z)
      call compute_Q(mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
         Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q_lb, Af)
      dQdZ_lb = differentiate_1(Q_lb)

      Y = set_Y(Y_is_positive, upper_bound_Z)
      call compute_Q(mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
         Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q_ub, Af)
      dQdZ_ub = differentiate_1(Q_ub)

      ! Check to make sure that the lower and upper bounds on Z actually bracket
      ! a solution to dQ/dZ = 0.
      has_root = .true.
      if (dQdZ_lb * dQdZ_ub > 0d0) then
         if (report) then
            write(*,*) 'TDC Error. Initial Z window does not bracket a solution.'
            write(*,*) 'Q(Lower Z)',Q_lb%val
            write(*,*) 'Q(Upper Z)',Q_ub%val
            write(*,*) 'dQdZ(Lower Z)',dQdZ_lb%val
            write(*,*) 'dQdZ(Upper Z)',dQdZ_ub%val
            write(*,*) 'scale', scale
            write(*,*) 'tolerance', residual_tolerance
            write(*,*) 'Y', Y%val
            write(*,*) 'dYdZ', Y%d1val1
            write(*,*) 'exp(Z)', exp(Z%val)
            write(*,*) 'Z', Z%val
            write(*,*) 'A0', A0%val
            write(*,*) 'c0', c0%val
            write(*,*) 'L', L%val
            write(*,*) 'L0', L0%val
            write(*,*) 'grada', grada%val
            write(*,*) 'gradL', gradL%val
            write(*,'(A)')
         end if
         has_root = .false.
         return
      end if

      do iter=1,max_iter
         Z_new = (upper_bound_Z + lower_bound_Z) / 2d0
         Y = set_Y(.false., Z_new)

         call compute_Q(mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
            Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q, Af)
         dQdZ = differentiate_1(Q)

         ! We only ever call this when Y < 0.
         ! In this regime, dQ/dZ can take on either sign, and has at most one stationary point.


         if (dQdZ > 0d0 .and. dQdZ_ub > 0d0) then
            upper_bound_Z = Z_new
            dQdZ_ub = dQdZ
         else if (dQdZ > 0d0 .and. dQdZ_lb > 0d0) then
            lower_bound_Z = Z_new
            dQdZ_lb = dQdZ
         else if (dQdZ < 0d0 .and. dQdZ_ub < 0d0) then
            upper_bound_Z = Z_new
            dQdZ_ub = dQdZ
         else if (dQdZ < 0d0 .and. dQdZ_lb < 0d0) then
            lower_bound_Z = Z_new
            dQdZ_lb = dQdZ
         end if

         if (upper_bound_Z - lower_bound_Z < bracket_tolerance) return         
      end do

   end subroutine dQdZ_bisection_search

   !> This routine performs a bisection search for the least-negative Y such that Af(Y) = 0.
   !! The search halts after max_iter iterations or when the bisection has converged to an interval
   !! less than Ztol in width.
   subroutine Af_bisection_search(&
         mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
         lower_bound_Z, upper_bound_Z, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Z, Af, ierr)
      ! Inputs
      real(dp), intent(in) :: mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt
      type(auto_diff_real_tdc), intent(in) :: A0, c0, L, L0, gradL, grada
      type(auto_diff_real_star_order1), intent(in) :: T, rho, dV, Cp, kap, Hp

      ! In/Out
      type(auto_diff_real_tdc), intent(inout) :: lower_bound_Z, upper_bound_Z

      ! Outputs
      type(auto_diff_real_tdc), intent(out) :: Z_new, Af
      integer, intent(out) :: ierr

      type(auto_diff_real_tdc) :: Y, Q
      real(dp), parameter :: bracket_tolerance = 1d-3
      integer, parameter :: max_iter = 30

      Y = set_Y(.false., upper_bound_Z)
      call compute_Q(mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
         Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q, Af)
      if (Af > 0) then ! d(Af)/dZ < 0, so if Af(upper_bound_Z) > 0 there's no solution in this interval.
         ierr = 1
         return
      end if

      do iter=1,max_iter
         Z_new = (upper_bound_Z + lower_bound_Z) / 2d0
         Y = set_Y(.false., Z_new)

         call compute_Q(mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, &
            Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q, Af)

         ! Y < 0 so increasing Y means decreasing Z.
         ! d(Af)/dY > 0 so d(Af)/dZ < 0.

         if (Af > 0d0) then ! Means we are at too-low Z.
            lower_bound_Z = Z_new
         else
            upper_bound_Z = Z_new
         end if

         if (upper_bound_Z - lower_bound_Z < bracket_tolerance) return
      end do

   end subroutine Af_bracket_search

   !> Computes the hyperbolic tangent of x in a way that is numerically safe.
   !!
   !! @param x Input
   !! @param z Output
   type(auto_diff_real_tdc) function safe_tanh(x) result(z)
      type(auto_diff_real_tdc), intent(in) :: x

      if (x > 50d0) then
         z = 1d0
      else if (x < -50d0) then
         z = -1d0
      else
         z = tanh(x)
      end if
   end function safe_tanh

   !> Computes the arctangent of y/x in a way that is numerically safe near x=0.
   !!
   !! @param x x coordinate for the arctangent.
   !! @param y y coordinate for the arctangent.
   !! @param z Polar angle z such that tan(z) = y / x.
   type(auto_diff_real_tdc) function safe_atan(x,y) result(z)
      type(auto_diff_real_tdc), intent(in) :: x,y
      type(auto_diff_real_tdc) :: x1, y1
      if (abs(x) < 1d-50) then
         ! x is basically zero, so for ~any non-zero y the ratio y/x is ~infinity.
         ! That means that z = +- pi. We want z to be positive, so we return pi.
         z = pi
      else
         z = atan(y/x)
      end if
   end function safe_atan

   !> The TDC newton solver needs higher-order partial derivatives than
   !! the star newton solver, because the TDC one needs to pass back a result
   !! which itself contains the derivatives that the star solver needs.
   !! These additional derivatives are provided by the auto_diff_real_tdc type.
   !!
   !! This method converts a auto_diff_real_star_order1 variable into a auto_diff_real_tdc,
   !! setting the additional partial derivatives to zero. This 'upgrades' variables storing
   !! stellar structure to a form the TDC solver can use.
   !!
   !! @param K_in, input, an auto_diff_real_star_order1 variable
   !! @param K, output, an auto_diff_real_tdc variable.
   type(auto_diff_real_tdc) function convert(K_in) result(K)
      type(auto_diff_real_star_order1), intent(in) :: K_in
      K%val = K_in%val
      K%d1Array(1:auto_diff_star_num_vars) = K_in%d1Array(1:auto_diff_star_num_vars)
      K%d1val1 = 0d0
      K%d1val1_d1Array(1:auto_diff_star_num_vars) = 0d0
   end function convert

   !> The TDC newton solver needs higher-order partial derivatives than
   !! the star newton solver, because the TDC one needs to pass back a result
   !! which itself contains the derivatives that the star solver needs.
   !! These additional derivatives are provided by the auto_diff_real_tdc type.
   !!
   !! This method converts a auto_diff_real_tdc variable into a auto_diff_real_star_order1,
   !! dropping the additional partial derivatives which (after the TDC solver is done) are
   !! no longer needed. This allows the output of the TDC solver to be passed back to the star solver.
   !!
   !! @param K_in, input, an auto_diff_real_tdc variable
   !! @param K, output, an auto_diff_real_star_order1 variable.      
   type(auto_diff_real_star_order1) function unconvert(K_in) result(K)
      type(auto_diff_real_tdc), intent(in) :: K_in
      K%val = K_in%val
      K%d1Array(1:auto_diff_star_num_vars) = K_in%d1Array(1:auto_diff_star_num_vars)
   end function unconvert

end module tdc_support
