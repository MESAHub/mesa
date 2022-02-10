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
public :: set_Y, Q_bisection_search, dQdZ_bisection_search, Af_bisection_search, &
         convert, unconvert, safe_atan, safe_tanh, tdc_info, &
         eval_Af, eval_xis, compute_Q

   !> Stores the information which is required to evaluate TDC-related quantities and which
   !! do not depend on Y.
   !!
   !! @param report Write debug output if true, not if false.
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
   !! @param Gamma Gamma is the MLT Gamma efficiency parameter, which we evaluate in steady state from MLT.
   type tdc_info
      logical :: report
      real(dp) :: mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt
      type(auto_diff_real_tdc) :: A0, c0, L, L0, gradL, grada
      type(auto_diff_real_star_order1) :: T, rho, dV, Cp, kap, Hp, Gamma
   end type tdc_info

contains

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
   !! Monotonicity is assumed, not verified! This means failures can occur if monotonicity fails.
   !! The search continues until the domain is narrowed to less than a width of bracket_tolerance,
   !! or until more than max_iter iterations have been taken. Because this is just used to get us in
   !! the right ballpark, bracket_tolerance is set quite wide, to 1.
   !! 
   !! There is a check at the start to verify that Q takes on opposite signs on either end of the
   !! domain. This is allows us to bail early if there is no root in the domain.
   !!
   !! @param info tdc_info type storing various quantities that are independent of Y.
   !! @param Y_is_positive True if Y > 0, False otherwise.
   !! @param lower_bound_Z Lower bound on Z. Note that this is an input *and* an output.
   !! @param upper_bound_Z Upper bound on Z. Note that this is an input *and* an output.
   !! @param Z Output: the midpoint of the final bounds.
   !! @param ierr 0 if everything worked, 1 if the domain does not contain a root.
   subroutine Q_bisection_search(info, Y_is_positive, lower_bound_Z, upper_bound_Z, Z, ierr)
      ! Inputs
      type(tdc_info), intent(in) :: info
      logical, intent(in) :: Y_is_positive

      ! Outputs
      type(auto_diff_real_tdc), intent(inout) :: lower_bound_Z, upper_bound_Z
      type(auto_diff_real_tdc), intent(out) :: Z
      integer, intent(out) :: ierr

      ! Parameters
      real(dp), parameter :: bracket_tolerance = 1d0
      integer, parameter :: max_iter = 30

      ! Intermediates
      type(auto_diff_real_tdc) :: Y, Af, Q, Q_ub, Q_lb
      integer :: iter

      ierr = 0

      Y = set_Y(Y_is_positive, lower_bound_Z)
      call compute_Q(info, Y, Q_lb, Af)

      Y = set_Y(Y_is_positive, upper_bound_Z)
      call compute_Q(info, Y, Q_ub, Af)

      ! Check to make sure that the lower and upper bounds on Z actually bracket
      ! a solution to Q(Y(Z)) = 0.
      if (Q_lb * Q_ub > 0d0) then
         if (info%report) then
            write(*,*) 'Q bisection error. Initial Z window does not bracket a solution.'
            write(*,*) 'Q(Lower Z)',Q_lb%val
            write(*,*) 'Q(Upper Z)',Q_ub%val
            write(*,*) 'tolerance', bracket_tolerance
            write(*,*) 'Y', Y%val
            write(*,*) 'dYdZ', Y%d1val1
            write(*,*) 'exp(Z)', exp(Z%val)
            write(*,*) 'Z', Z%val
            write(*,*) 'A0', info%A0%val
            write(*,*) 'c0', info%c0%val
            write(*,*) 'L', info%L%val
            write(*,*) 'L0', info%L0%val
            write(*,*) 'grada', info%grada%val
            write(*,*) 'gradL', info%gradL%val
            write(*,'(A)')
         end if
         ierr = 1
         return
      end if

      do iter=1,max_iter
         Z = (upper_bound_Z + lower_bound_Z) / 2d0
         Y = set_Y(Y_is_positive, Z)

         call compute_Q(info, Y, Q, Af)

         if (Q > 0d0 .and. Q_ub > 0d0) then
            upper_bound_Z = Z
            Q_ub = Q
         else if (Q > 0d0 .and. Q_lb > 0d0) then
            lower_bound_Z = Z
            Q_lb = Q
         else if (Q < 0d0 .and. Q_ub < 0d0) then
            upper_bound_Z = Z
            Q_ub = Q
         else if (Q < 0d0 .and. Q_lb < 0d0) then
            lower_bound_Z = Z
            Q_lb = Q
         end if

         if (upper_bound_Z - lower_bound_Z < bracket_tolerance) then
            Z = (upper_bound_Z + lower_bound_Z) / 2d0
            call compute_Q(info, Y, Q, Af)
            return
         end if
      end do

   end subroutine Q_bisection_search

   !> This routine performs a bisection search for dQ/dZ=0 with Y < 0.
   !! The domain is assumed to be restricted to have Af > 0, so that dQ/dZ is
   !! continuous and monotone. This is checked, and if it fails a MESA ERROR is called
   !! rather than throwing an ierr because if the rest of the TDC logic is correct this should
   !! not happen.
   !!
   !! The search continues until the domain is narrowed to less than a width of bracket_tolerance (1d-4),
   !! or until more than max_iter iterations have been taken.
   !! 
   !! There is a check at the start to verify that dQ/dZ takes on opposite signs on either end of the
   !! domain. This is allows us to bail early if there is no root in the domain.
   !!
   !! @param info tdc_info type storing various quantities that are independent of Y.
   !! @param lower_bound_Z Lower bound on Z. Note that this is an input *and* an output.
   !! @param upper_bound_Z Upper bound on Z. Note that this is an input *and* an output.
   !! @param Z Output: the midpoint of the final bounds.
   !! @param has_root True if a root was found, False otherwise.
   subroutine dQdZ_bisection_search(info, lower_bound_Z_in, upper_bound_Z_in, Z, has_root)
      ! Inputs
      type(tdc_info), intent(in) :: info
      type(auto_diff_real_tdc), intent(in) :: lower_bound_Z_in, upper_bound_Z_in

      ! Outputs
      type(auto_diff_real_tdc), intent(out) :: Z
      logical, intent(out) :: has_root

      ! Parameters
      real(dp), parameter :: bracket_tolerance = 1d-4
      integer, parameter :: max_iter = 50

      ! Intermediates
      type(auto_diff_real_tdc) :: lower_bound_Z, upper_bound_Z
      type(auto_diff_real_tdc) :: Y, Af, Q, Q_lb, Q_ub, dQdZ, dQdZ_lb, dQdZ_ub
      integer :: iter

      ! Set up
      lower_bound_Z = lower_bound_Z_in!lower_bound_Z_in
      lower_bound_Z%d1val1 = 1d0
      upper_bound_Z = upper_bound_Z_in
      upper_bound_Z%d1val1 = 1d0

      ! Check bounds
      Y = set_Y(.false., lower_bound_Z)
      call compute_Q(info, Y, Q_lb, Af)
      if (Af == 0) then
         write(*,*) 'Z_lb, A0, Af', lower_bound_Z%val, info%A0%val, Af%val
         call mesa_error(__FILE__,__LINE__,'bad call to tdc_support dQdZ_bisection_search: Af == 0.')
      end if
      dQdZ_lb = differentiate_1(Q_lb)

      Y = set_Y(.false., upper_bound_Z)
      call compute_Q(info, Y, Q_ub, Af)
      if (Af == 0) then
         write(*,*) 'Z_ub, A0, Af', lower_bound_Z%val, info%A0%val, Af%val
         call mesa_error(__FILE__,__LINE__,'bad call to tdc_support dQdZ_bisection_search: Af == 0.')
      end if
      dQdZ_ub = differentiate_1(Q_ub)

      ! Check to make sure that the lower and upper bounds on Z actually bracket
      ! a solution to dQ/dZ = 0.
      has_root = .true.
      if (dQdZ_lb * dQdZ_ub > 0d0) then
         if (info%report) then
            write(*,*) 'dQdZ bisection error. Initial Z window does not bracket a solution.'
            write(*,*) 'Q(Lower Z)',Q_lb%val
            write(*,*) 'Q(Upper Z)',Q_ub%val
            write(*,*) 'dQdZ(Lower Z)',dQdZ_lb%val
            write(*,*) 'dQdZ(Upper Z)',dQdZ_ub%val
            write(*,*) 'tolerance', bracket_tolerance
            write(*,*) 'Y', Y%val
            write(*,*) 'dYdZ', Y%d1val1
            write(*,*) 'exp(Z)', exp(Z%val)
            write(*,*) 'Z', Z%val
            write(*,*) 'A0', info%A0%val
            write(*,*) 'c0', info%c0%val
            write(*,*) 'L', info%L%val
            write(*,*) 'L0', info%L0%val
            write(*,*) 'grada', info%grada%val
            write(*,*) 'gradL', info%gradL%val
            write(*,'(A)')
         end if
         has_root = .false.
         return
      end if

      ! Bisection search
      do iter=1,max_iter
         Z = (upper_bound_Z + lower_bound_Z) / 2d0
         Z%d1val1 = 1d0
         Y = set_Y(.false., Z)

         call compute_Q(info, Y, Q, Af)
         dQdZ = differentiate_1(Q)

         ! We only ever call this when Y < 0.
         ! In this regime, dQ/dZ can take on either sign, and has at most one stationary point.

         if (info%report) write(*,*) 'Bisecting dQdZ. Z, dQdZ, Z_lb, dQdZ_lb, Z_ub, dQdZ_ub', Z%val, dQdZ%val, lower_bound_Z%val, dQdZ_lb%val, upper_bound_Z%val, dQdZ_ub%val

         if (dQdZ > 0d0 .and. dQdZ_ub > 0d0) then
            upper_bound_Z = Z
            dQdZ_ub = dQdZ
         else if (dQdZ > 0d0 .and. dQdZ_lb > 0d0) then
            lower_bound_Z = Z
            dQdZ_lb = dQdZ
         else if (dQdZ < 0d0 .and. dQdZ_ub < 0d0) then
            upper_bound_Z = Z
            dQdZ_ub = dQdZ
         else if (dQdZ < 0d0 .and. dQdZ_lb < 0d0) then
            lower_bound_Z = Z
            dQdZ_lb = dQdZ
         end if

         if (upper_bound_Z - lower_bound_Z < bracket_tolerance) then
            Z = (upper_bound_Z + lower_bound_Z) / 2d0
            call compute_Q(info, Y, Q, Af)
            return         
         end if
      end do

   end subroutine dQdZ_bisection_search

   !> This routine performs a bisection search for the least-negative Y such that Af(Y) = 0.
   !! Once we find that, we return a Y that is just slightly less negative so that Af > 0.
   !! This is important for our later bisection of dQ/dZ, because we want the upper-Z end of
   !! the domain we bisect to have dQ/dZ < 0, which means it has to capture the fact that d(Af)/dZ < 0,
   !! and so the Z we return has to be on the positive-Af side of the discontinuity in dQdZ.
   !!
   !! Note that monotonicity is assumed, not verified!
   !! Af(Y) is monotonic in Y for negative Y by construction, so this shouldn't be a problem.
   !!
   !! The search continues until the domain is narrowed to less than a width of bracket_tolerance (1d-4),
   !! or until more than max_iter iterations have been taken.
   !! 
   !! There is a check at the start to verify that Af == 0 at the most-negative end of the domain.
   !! This is allows us to bail early if there is no root in the domain.
   !!
   !! @param info tdc_info type storing various quantities that are independent of Y.
   !! @param lower_bound_Z Lower bound on Z. Note that this is an input *and* an output.
   !! @param upper_bound_Z Upper bound on Z. Note that this is an input *and* an output.
   !! @param Z Output: The bound on the positive-Af side of the root (which always is the lower bound).
   !! @param Af Output: Af(Z).
   !! @param ierr 0 if everything worked, 1 if the domain does not contain a root.
   subroutine Af_bisection_search(info, lower_bound_Z_in, upper_bound_Z_in, Z, Af, ierr)
      ! Inputs
      type(tdc_info), intent(in) :: info
      type(auto_diff_real_tdc), intent(in) :: lower_bound_Z_in, upper_bound_Z_in

      ! Outputs
      integer, intent(out) :: ierr
      type(auto_diff_real_tdc), intent(out) :: Z, Af

      ! Parameters
      real(dp), parameter :: bracket_tolerance = 1d-4
      integer, parameter :: max_iter = 50

      ! Intermediates
      type(auto_diff_real_tdc) :: lower_bound_Z, upper_bound_Z
      type(auto_diff_real_tdc) :: Y, Q
      integer :: iter

      ! Set up
      ierr = 0
      lower_bound_Z = lower_bound_Z_in
      upper_bound_Z = upper_bound_Z_in

      Y = set_Y(.false., upper_bound_Z)
      call compute_Q(info, Y, Q, Af)
      if (Af > 0) then ! d(Af)/dZ < 0, so if Af(upper_bound_Z) > 0 there's no solution in this interval.
         ierr = 1
         return
      end if

      Y = set_Y(.false., lower_bound_Z)
      call compute_Q(info, Y, Q, Af)
      if (Af == 0) then
         ! We want to find Z such that Af(Z) is just barely above zero.
         ! We do this by finding Z such that Af(Z) == 0, then backing off to slightly smaller Z.
         ! Because d(Af)/dZ < 0, this gives a Z such that Af(Z) > 0.
         ! Hence, if Af(lower_bound_Z) == 0 then Af = 0 uniformly in this interval and we cannot
         ! return Z such that Af(Z) is just barely non-zero.
         ierr = 2
         return
      end if

      do iter=1,max_iter
         Z = (upper_bound_Z + lower_bound_Z) / 2d0
         Y = set_Y(.false., Z)

         call compute_Q(info, Y, Q, Af)

         ! Y < 0 so increasing Y means decreasing Z.
         ! d(Af)/dY > 0 so d(Af)/dZ < 0.

         if (Af > 0d0) then ! Means we are at too-low Z.
            lower_bound_Z = Z
         else
            upper_bound_Z = Z
         end if

         if (upper_bound_Z - lower_bound_Z < bracket_tolerance) then
            ! We return the lower bound because this is guaranteed to have Af > 0 (just barely).
            ! This is important for our later bisection of dQ/dZ, because we want the upper-Z end of
            ! the domain we bisect to have dQ/dZ < 0, which means it has to capture the fact that d(Af)/dZ < 0.
            Z = lower_bound_Z
            call compute_Q(info, Y, Q, Af)
            return
         end if
      end do

   end subroutine Af_bisection_search

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

   !> Q is the residual in the TDC equation, namely:
   !!
   !! Q = (L - L0 * gradL) - (L0 + c0 * Af) * Y
   !!
   !! @param info tdc_info type storing various quantities that are independent of Y.
   !! @param Y superadiabaticity
   !! @param Q The residual of the above equation (an output).
   !! @param Af The final convection speed (an output).
   subroutine compute_Q(info, Y, Q, Af)
      type(tdc_info), intent(in) :: info
      type(auto_diff_real_tdc), intent(in) :: Y
      type(auto_diff_real_tdc), intent(out) :: Q, Af
      type(auto_diff_real_tdc) :: xi0, xi1, xi2, Y_env

      ! Y = grad-gradL
      ! Gamma=(grad-gradE)/(gradE-gradL)
      ! So
      ! Y_env = grad-gradE = (grad-gradL)*Gamma/(1+Gamma) = Y*Gamma/(1+Gamma)
      ! So overall we just multiply the Y by Gamma/(1+Gamma) to get Y_env.
      !
      ! We only use Y_env /= Y when Y > 0 (i.e. the system is convectively unstable)
      ! because we only have a Gamma from MLT in that case.
      ! so when Y < 0 we just use Y_env = Y.
      if (Y > 0) then
         Y_env = Y * convert(info%Gamma/(1+info%Gamma))
      else
         Y_env = Y
      end if

      ! Y_env sets the acceleration of blobs.
      call eval_xis(info, Y_env, xi0, xi1, xi2)          
      Af = eval_Af(info%dt, info%A0, xi0, xi1, xi2)

      ! Y_env sets the convective flux but not the radiative flux.
      Q = (info%L - info%L0*info%gradL) - info%L0 * Y - info%c0*Af*Y_env

   end subroutine compute_Q

   !! Calculates the coefficients of the TDC velocity equation.
   !! The velocity equation is
   !!
   !! 2 dw/dt = xi0 + w * xi1 + w**2 * xi2 - Lambda
   !!
   !! where Lambda (currently set to zero) captures coupling between cells.
   !!
   !! The coefficients xi0/1/2 are given by
   !!
   !! xi0 = (epsilon_q + S * Y) / w
   !! xi1 = (-D_r + p_turb dV/dt) / w^2    [here V = 1/rho is the volume]
   !! xi2 = -D / w^3
   !!
   !! Note that these terms are evaluated on faces, because the convection speed
   !! is evaluated on faces. As a result all the inputs must either natively
   !! live on faces or be interpolated to faces.
   !!
   !! This function also has some side effects in terms of storing some of the terms
   !! it calculates for plotting purposes.
   !!
   !! @param info tdc_info type storing various quantities that are independent of Y.
   !! @param Y superadiabaticity
   !! @param xi0 Output, the constant term in the convective velocity equation.
   !! @param xi1 Output, the prefactor of the linear term in the convective velocity equation.
   !! @param xi2 Output, the prefactor of the quadratic term in the convective velocity equation.
   subroutine eval_xis(info, Y, xi0, xi1, xi2) 
      ! eval_xis sets up Y with partial wrt Z
      ! so results come back with partials wrt Z
      type(tdc_info), intent(in) :: info
      type(auto_diff_real_tdc), intent(in) :: Y
      type(auto_diff_real_tdc), intent(out) :: xi0, xi1, xi2
      type(auto_diff_real_tdc) :: S0, D0, DR0
      type(auto_diff_real_star_order1) :: gammar_div_alfa, Pt0, dVdt
      real(dp), parameter :: x_ALFAS = (1.d0/2.d0)*sqrt_2_div_3
      real(dp), parameter :: x_CEDE  = (8.d0/3.d0)*sqrt_2_div_3
      real(dp), parameter :: x_ALFAP = 2.d0/3.d0
      real(dp), parameter :: x_GAMMAR = 2.d0*sqrt(3.d0)

      S0 = convert(x_ALFAS*info%mixing_length_alpha*info%Cp*info%T/info%Hp)*info%grada
      S0 = S0*Y
      D0 = convert(info%alpha_TDC_DAMP*x_CEDE/(info%mixing_length_alpha*info%Hp))
      gammar_div_alfa = info%alpha_TDC_DAMPR*x_GAMMAR/(info%mixing_length_alpha*info%Hp)
      DR0 = convert(4d0*boltz_sigma*pow2(gammar_div_alfa)*pow3(info%T)/(pow2(info%rho)*info%Cp*info%kap))
      Pt0 = info%alpha_TDC_PtdVdt*x_ALFAP*info%rho
      dVdt = info%dV/info%dt

      xi0 = S0
      xi1 = -(DR0 + convert(Pt0*dVdt))
      xi2 = -D0
   end subroutine eval_xis

   !! Calculates the solution to the TDC velocity equation.
   !! The velocity equation is
   !!
   !! 2 dw/dt = xi0 + w * xi1 + w**2 * xi2 - Lambda
   !!
   !! where Lambda (currently set to zero) captures coupling between cells.
   !! The xi0/1/2 variables are constants for purposes of solving this equation.
   !!
   !! An important related parameter is J:
   !! 
   !! J^2 = xi1^2 - 4 * xi0 * xi2
   !!
   !! When J^2 > 0 the solution for w is hyperbolic in time.
   !! When J^2 < 0 the solution is trigonometric, behaving like tan(J * t / 4 + c).
   !!
   !! In the trigonometric branch note that once the solution passes through its first
   !! root the solution for all later times is w(t) = 0. This is not a solution to the
   !! convective velocity equation as written above, but *is* a solution to the convective
   !! energy equation (multiply the velocity equation by w), which is the one we actually
   !! want to solve (per Radek Smolec's thesis). We just solve the velocity form
   !! because it's more convenient.
   !!
   !! @param dt Time-step
   !! @param A0 convection speed from the start of the step (cm/s)
   !! @param xi0 The constant term in the convective velocity equation.
   !! @param xi1 The prefactor of the linear term in the convective velocity equation.
   !! @param xi2 The prefactor of the quadratic term in the convective velocity equation.            
   !! @param Af Output, the convection speed at the end of the step (cm/s)
   function eval_Af(dt, A0, xi0, xi1, xi2) result(Af)
      real(dp), intent(in) :: dt    
      type(auto_diff_real_tdc), intent(in) :: A0, xi0, xi1, xi2
      type(auto_diff_real_tdc) :: Af ! output
      type(auto_diff_real_tdc) :: J2, J, Jt, Jt4, num, den, y_for_atan, root, lk 

      J2 = pow2(xi1) - 4d0 * xi0 * xi2

      if (J2 > 0d0) then ! Hyperbolic branch
         J = sqrt(J2)
         Jt = dt * J
         Jt4 = 0.25d0 * Jt
         num = safe_tanh(Jt4) * (2d0 * xi0 + A0 * xi1) + A0 * J
         den = safe_tanh(Jt4) * (xi1 + 2d0 * A0 * xi2) - J
         Af = num / den 
         if (Af < 0d0) then
            Af = -Af
         end if
      else if (J2 < 0d0) then ! Trigonometric branch
         J = sqrt(-J2)
         Jt = dt * J

         ! This branch contains decaying solutions that reach A = 0, at which point
         ! they switch onto the 'zero' branch. So we have to calculate the position of
         ! the first root to check it against dt.
         y_for_atan = xi1 + 2d0 * A0 * xi2
         root = safe_atan(J, xi1) - safe_atan(J, y_for_atan)

         ! The root enters into a tangent, so we can freely shift it by pi and
         ! get another root. We care about the first positive root, and the above prescription
         ! is guaranteed to give an answer between (-2*pi,2*pi) because atan produces an answer in [-pi,pi],
         ! so we add/subtract a multiple of pi to get the root into [0,pi).
         if (root > pi) then
            root = root - pi
         else if (root < -pi) then
            root = root + 2d0*pi
         else if (root < 0d0) then
            root = root + pi
         end if

         if (0.25d0 * Jt < root) then
            num = -xi1 + J * tan(0.25d0 * Jt + atan(y_for_atan / J)) 
            den = 2d0 * xi2
            Af = num / den
         else
            Af = 0d0
         end if
      else ! if (J2 == 0d0) then         
         Af = A0            
      end if

   end function eval_Af

end module tdc_support
