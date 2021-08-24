! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

use star_private_def
use const_def
use num_lib
use utils_lib
use auto_diff_support
use star_utils

implicit none

private
public :: set_TDC, check_if_must_fall_back_to_MLT

contains

   subroutine set_TDC(s, k, &
            mixing_length_alpha, cgrav, m, report, &
            mixing_type, L, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, Y_face, gradT, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(in) :: mixing_length_alpha, cgrav, m
      type(auto_diff_real_star_order1), intent(in) :: &
         L, r, P, T, rho, dV, Cp, opacity, scale_height, gradL, grada
      logical, intent(in) :: report
      type(auto_diff_real_star_order1),intent(out) :: conv_vel, Y_face, gradT
      integer, intent(out) :: mixing_type, ierr
      include 'formats'
      call get_TDC_solution(s, k, &
         mixing_length_alpha, cgrav, m, report, &
         mixing_type, L, r, P, T, rho, dV, Cp, opacity, &
         scale_height, gradL, grada, conv_vel, Y_face, ierr)
      if (ierr /= 0) then
         write(*,2) 'get_TDC_solution failed in set_TDC', k
         write(*,*) 'Repeating call with reporting on.'
         call get_TDC_solution(s, k, &
            mixing_length_alpha, cgrav, m, .true., &
            mixing_type, L, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, Y_face, ierr)
         stop 'get_TDC_solution failed in set_TDC'
      end if
      gradT = Y_face + grada
   end subroutine set_TDC       

   !> Determines if it is safe (physically) to use TDC instead of MLT.
   !!
   !! Currently we only know we have to fall back to MLT in cells that get touched
   !! by adjust_mass, because there the convection speeds at the start of the
   !! step can be badly out of whack.
   !!
   !! @param s star pointer
   !! @param k face index
   !! @param fallback False if we can use TDC, True if we can fall back to MLT.
   logical function check_if_must_fall_back_to_MLT(s, k) result(fallback)
      type (star_info), pointer :: s
      integer, intent(in) :: k

      fallback = .false.
      if (abs(s%mstar_dot) > 1d-99 .and. k < s% k_const_mass) then
         fallback = .true.
      end if
   end function check_if_must_fall_back_to_MLT

   type(auto_diff_real_tdc) function set_Y(Y_is_positive, Z) result(Y)
      logical, intent(in) :: Y_is_positive
      type(auto_diff_real_tdc), intent(in) :: Z
      if (Y_is_positive) then
         Y = exp(Z)
      else
         Y = -exp(Z)
      end if
   end function set_Y

   subroutine TDC_bracket_search(s, k, mixing_length_alpha, Y_is_positive, lower_bound_Z, upper_bound_Z, Q_ub, Q_lb, &
         c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada)
      type (star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(in) :: mixing_length_alpha
      logical, intent(in) :: Y_is_positive
      type(auto_diff_real_star_order1), intent(in) :: &
         c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada
      type(auto_diff_real_tdc), intent(inout) :: lower_bound_Z, upper_bound_Z
      type(auto_diff_real_tdc), intent(inout) :: Q_ub, Q_lb
      type(auto_diff_real_tdc) :: Z_new, Y, Af, Q

      Z_new = (upper_bound_Z + lower_bound_Z) / 2d0
      Y = set_Y(Y_is_positive, Z_new)

      call compute_Q(s, k, mixing_length_alpha, &
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

   end subroutine TDC_bracket_search


   subroutine get_TDC_solution(s, k, &
         mixing_length_alpha, cgrav, m, report, &
         mixing_type, L, r, P, T, rho, dV, Cp, kap, Hp, gradL, grada, cv, Y_face, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(in) :: mixing_length_alpha, cgrav, m
      type(auto_diff_real_star_order1), intent(in) :: &
         L, r, P, T, rho, dV, Cp, kap, Hp, gradL, grada
      logical, intent(in) :: report
      type(auto_diff_real_star_order1),intent(out) :: cv, Y_face
      integer, intent(out) :: mixing_type, ierr
      
      type(auto_diff_real_star_order1) :: A0, c0, L0
      type(auto_diff_real_tdc) :: Af, Y, Z, Q, Q_lb, Q_ub, Qc, Z_new, correction, lower_bound_Z, upper_bound_Z
      type(auto_diff_real_tdc) :: dQdZ, prev_dQdZ
      real(dp) ::  gradT, Lr, Lc, scale
      integer :: iter, line_iter, i
      logical :: converged, Y_is_positive, first_Q_is_positive, have_derivatives, corr_has_derivatives
      real(dp), parameter :: bracket_tolerance = 1d0
      real(dp), parameter :: correction_tolerance = 1d-13
      real(dp), parameter :: residual_tolerance = 1d-8
      real(dp), parameter :: alpha_c = (1d0/2d0)*sqrt_2_div_3
      integer, parameter :: max_iter = 200
      integer, parameter :: max_line_search_iter = 5
      include 'formats'

      ierr = 0
      if (mixing_length_alpha == 0d0 .or. k < 1 .or. s% dt <= 0d0) then
         stop 'bad call to TDC get_TDC_solution'
      end if         

      ! Set up inputs.
      c0 = mixing_length_alpha*alpha_c*rho*T*Cp*4d0*pi*pow2(r)
      L0 = (16d0*pi*crad*clight/3d0)*cgrav*m*pow4(T)/(P*kap) ! assumes QHSE for dP/dm
      if (s% okay_to_set_mlt_vc) then
         A0 = s% mlt_vc_old(k)/sqrt_2_div_3
      else
         A0 = s% mlt_vc(k)/sqrt_2_div_3
      end if

      ! Set scale for judging the solution to Q(Y)=0.
      ! Q has units of a luminosity, so the scale should be a luminosity.
      if (s% solver_iter == 0) then
         scale = max(abs(s% L(k)), 1d-3*maxval(s% L(1:s% nz)))
      else
         scale = max(abs(s% L_start(k)), 1d-3*maxval(s% L_start(1:s% nz)))
      end if

      ! Determine the sign of the solution.
      !
      ! If Q(Y=0) is positive then the luminosity is too great to be carried radiatively, so
      ! we'll necessarily have Y > 0.
      !
      ! If Q(Y=0) is negative then the luminosity can be carried by radiation alone, so we'll
      ! necessarily have Y < 0.
      Y = 0d0
      call compute_Q(s, k, mixing_length_alpha, &
         Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q, Af)
      if (Q > 0d0) then
         Y_is_positive = .true.
      else
         Y_is_positive = .false.
      end if

      ! We start by bisecting to find a narrow interval around the root.
      lower_bound_Z = -200d0
      upper_bound_Z = 100d0 

      Y = set_Y(Y_is_positive, lower_bound_Z)
      call compute_Q(s, k, mixing_length_alpha, &
         Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q_lb, Af)

      Y = set_Y(Y_is_positive, upper_bound_Z)
      call compute_Q(s, k, mixing_length_alpha, &
         Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q_ub, Af)

      ! Check to make sure that the lower and upper bounds on Z actually bracket
      ! a solution to Q(Y(Z)) = 0.
      if (Q_lb * Q_ub > 0d0) then
            write(*,*) 'TDC Error. Initial Z window does not bracket a solution.'
            write(*,*) 'Q(Lower Z)',Q_lb%val
            write(*,*) 'Q(Upper Z)',Q_ub%val
            write(*,2) 'Q(Y=0)', k, Q%val
            write(*,2) 'scale', k, scale
            write(*,2) 'Q/scale', k, Q%val/scale
            write(*,2) 'tolerance', k, residual_tolerance
            write(*,2) 'dQdZ', k, dQdZ%val
            write(*,2) 'Y', k, Y%val
            write(*,2) 'dYdZ', k, Y%d1val1
            write(*,2) 'exp(Z)', k, exp(Z%val)
            write(*,2) 'Z', k, Z%val
            write(*,2) 'Af', k, Af%val
            write(*,2) 'dAfdZ', k, Af%d1val1
            write(*,2) 'A0', k, A0%val
            write(*,2) 'c0', k, c0%val
            write(*,2) 'L', k, L%val
            write(*,2) 'L0', k, L0%val
            write(*,2) 'grada', k, grada%val
            write(*,2) 'gradL', k, gradL%val
            write(*,*)
         ierr = 1
         return
      end if

      ! Perform bisection search.
      do iter = 1, max_iter
         call TDC_bracket_search(s, k, mixing_length_alpha, Y_is_positive, lower_bound_Z, upper_bound_Z, Q_ub, Q_lb, &
                                 c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada)
         if (upper_bound_Z - lower_bound_Z < bracket_tolerance) exit
      end do
      Z = (upper_bound_Z + lower_bound_Z) / 2d0
      Z%d1val1 = 1d0 ! Set derivative dZ/dZ=1 for Newton iterations.
      if (report) write(*,2) 'Z from bracket search', k, Z%val

      ! Now we refine the solution with a Newton solve.
      ! This also let's us pick up the derivative of the solution with respect
      ! to input parameters.

      ! Initialize starting values for TDC Newton iterations.
      dQdz = 0d0
      converged = .false.
      have_derivatives = .false. ! Tracks if we've done at least one un-clipped Newton iteration.
                                 ! Need to do this before returning to endow Y with partials
                                 ! with respect to the structure variables.
      do iter = 1, max_iter
         Y = set_Y(Y_is_positive, Z)
         call compute_Q(s, k, mixing_length_alpha, &
            Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q, Af)

         if (abs(Q%val)/scale <= residual_tolerance .and. have_derivatives) then
            ! Can't exit on the first iteration, otherwise we have no derivative information.
            if (report) write(*,2) 'converged', iter, abs(Q%val)/scale, residual_tolerance
            converged = .true.
            exit
         end if

         if (gradL > 0d0 .and. Y_is_positive) then
            ! We use the fact that Q(Y) is monotonic for Y > 0 to produce iteratively refined bounds on Q.
            if (Q > 0d0) then
               ! Q(Y) is monotonic so this means Z is a lower-bound.
               lower_bound_Z = Z
            else
               ! Q(Y) is monotonic so this means Z is an upper-bound.
               upper_bound_Z = Z
            end if
         end if

         prev_dQdZ = dQdZ
         dQdZ = differentiate_1(Q)
         if (is_bad(dQdZ%val) .or. abs(dQdZ%val) < 1d-99) then
            ierr = 1
            exit
         end if

         if (prev_dQdZ * dQdZ < 0d0) then ! Means we're sitting around a stationary point.
            call TDC_bracket_search(s, k, mixing_length_alpha, Y_is_positive, lower_bound_Z, upper_bound_Z, Q_ub, Q_lb, &
                                    c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada)
            Z_new = (upper_bound_Z + lower_bound_Z) / 2d0
            have_derivatives = .false. ! Bracket search eliminates derivative information.
         else
            correction = -Q/dQdz
            corr_has_derivatives = .true.

            ! Clip steps.
            ! Because Z = log|Y|, large steps in Z
            ! correspond to enormous steps in Y, and
            ! usually indicate that something went wrong.
            if (abs(correction) > 2d0) then
               ! If we end up clipping the correction it loses derivative
               ! information.
               corr_has_derivatives = .false.
            end if
            correction = max(correction, -2d0)
            correction = min(correction, 2d0)

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

               call compute_Q(s, k, mixing_length_alpha, &
               Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Qc, Af)

               if (abs(Qc) < abs(Q)) then
                  exit
               else
                  correction = 0.5d0 * correction
               end if
            end do
         end if

         if (report) write(*,3) 'i, li, Z_new, Z, low_bnd, upr_bnd, Q, dQdZ, pdQdZ, corr', iter, line_iter, &
            Z_new%val, Z%val, lower_bound_Z%val, upper_bound_Z%val, Q%val, dQdZ%val, prev_dQdZ%val, correction%val
         Z_new%d1val1 = 1d0 ! Ensures that dZ/dZ = 1.
         Z = Z_new

         Y = set_Y(Y_is_positive,Z)

      end do

      if (.not. converged) then
         ierr = 1
         if (report .or. s% x_integer_ctrl(19) <= 0) then
         !$OMP critical (tdc_crit0)
            write(*,5) 'failed get_TDC_solution k slvr_iter model TDC_iter', &
               k, s% solver_iter, s% model_number, iter
            write(*,2) 'Q', k, Q%val
            write(*,2) 'scale', k, scale
            write(*,2) 'Q/scale', k, Q%val/scale
            write(*,2) 'tolerance', k, residual_tolerance
            write(*,2) 'dQdZ', k, dQdZ%val
            write(*,2) 'Y', k, Y%val
            write(*,2) 'dYdZ', k, Y%d1val1
            write(*,2) 'exp(Z)', k, exp(Z%val)
            write(*,2) 'Z', k, Z%val
            write(*,2) 'Af', k, Af%val
            write(*,2) 'dAfdZ', k, Af%d1val1
            write(*,2) 'A0', k, A0%val
            write(*,2) 'c0', k, c0%val
            write(*,2) 'L', k, L%val
            write(*,2) 'L0', k, L0%val
            write(*,2) 'grada', k, grada%val
            write(*,2) 'gradL', k, gradL%val
            write(*,*)
         !$OMP end critical (tdc_crit0)
         end if
         return
      end if

      ! Process Y into the various outputs.
      cv = sqrt_2_div_3*unconvert(Af)   
      Y_face = unconvert(Y)
      gradT = Y_face%val + gradL%val
      Lr = L0%val*gradT
      Lc = L%val - Lr
      if (cv > 0d0) then
         mixing_type = convective_mixing
      else
         mixing_type = no_mixing
      end if
      if (k > 0) s% tdc_num_iters(k) = iter          
   end subroutine get_TDC_solution
         

   !> Q is the residual in the TDC equation, namely:
   !!
   !! Q = (L - L0 * gradL) - (L0 + c0 * Af) * Y
   !!
   !! @param s star pointer
   !! @param k face index
   !! @param Y superadiabaticity
   !! @param c0_in A proportionality factor for the convective luminosity
   !! @param L_in luminosity
   !! @param L0_in L0 = (Lrad / grad_rad) is the luminosity radiation would carry if dlnT/dlnP = 1.
   !! @param A0 Initial convection speed
   !! @param T Temperature
   !! @param rho Density (g/cm^3)
   !! @param dV 1/rho_face - 1/rho_start_face (change in specific volume at the face)
   !! @param Cp Heat capacity
   !! @param kap Opacity
   !! @param Hp Pressure scale height
   !! @param gradL_in gradL is the neutrally buoyant dlnT/dlnP (= grad_ad + grad_mu),
   !! @param grada_in grada is the adiabatic dlnT/dlnP,
   !! @param Q The residual of the above equaiton (an output).
   !! @param Af The final convection speed (an output).
   subroutine compute_Q(s, k, mixing_length_alpha, &
         Y, c0_in, L_in, L0_in, A0, T, rho, dV, Cp, kap, Hp, gradL_in, grada_in, Q, Af)
      type (star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(in) :: mixing_length_alpha
      type(auto_diff_real_star_order1), intent(in) :: &
         c0_in, L_in, L0_in, A0, T, rho, dV, Cp, kap, Hp, gradL_in, grada_in
      type(auto_diff_real_tdc), intent(in) :: Y
      type(auto_diff_real_tdc), intent(out) :: Q, Af
      type(auto_diff_real_tdc) :: xi0, xi1, xi2, c0, L0, L, gradL

      call eval_xis(s, k, mixing_length_alpha, &
         Y, T, rho, Cp, dV, kap, Hp, grada_in, xi0, xi1, xi2)          

      Af = eval_Af(s, k, A0, xi0, xi1, xi2)
      L = convert(L_in)
      L0 = convert(L0_in)
      gradL = convert(gradL_in)
      c0 = convert(c0_in)
      Q = (L - L0*gradL) - (L0 + c0*Af)*Y

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
   !! @param s star pointer
   !! @param k face index
   !! @param mixing_length_alpha mixing length parameter
   !! @param Y superadiabaticity
   !! @param T temperature (K)
   !! @param rho density (g/cm^3)
   !! @param Cp heat capacity (erg/g/K)
   !! @param dV 1/rho_face - 1/rho_start_face (change in specific volume at the face)
   !! @param kap opacity (cm^2/g)
   !! @param Hp pressure scale height (cm)
   !! @param grada is the adiabatic dlnT/dlnP.
   !! @param xi0 Output, the constant term in the convective velocity equation.
   !! @param xi1 Output, the prefactor of the linear term in the convective velocity equation.
   !! @param xi2 Output, the prefactor of the quadratic term in the convective velocity equation.
   subroutine eval_xis(s, k, mixing_length_alpha, &
         Y, T, rho, Cp, dV, kap, Hp, grada, xi0, xi1, xi2) 
      ! eval_xis sets up Y with partial wrt Z
      ! so results come back with partials wrt Z
      type (star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(in) :: mixing_length_alpha
      type(auto_diff_real_star_order1), intent(in) :: T, rho, dV, Cp, kap, Hp, grada
      type(auto_diff_real_tdc), intent(in) :: Y
      type(auto_diff_real_tdc), intent(out) :: xi0, xi1, xi2
      type(auto_diff_real_tdc) :: S0, D0, DR0
      type(auto_diff_real_star_order1) :: gammar_div_alfa, Pt0, dVdt
      real(dp), parameter :: x_ALFAS = (1.d0/2.d0)*sqrt_2_div_3
      real(dp), parameter :: x_CEDE  = (8.d0/3.d0)*sqrt_2_div_3
      real(dp), parameter :: x_ALFAP = 2.d0/3.d0
      real(dp), parameter :: x_GAMMAR = 2.d0*sqrt(3.d0)
      include 'formats'

      S0 = convert(x_ALFAS*mixing_length_alpha*Cp*T*grada/Hp)
      S0 = S0*Y
      D0 = convert(s% alpha_TDC_DAMP*x_CEDE/(mixing_length_alpha*Hp))
      if (s% alpha_TDC_DAMPR == 0d0) then
         DR0 = 0d0
      else
         gammar_div_alfa = s% alpha_TDC_DAMPR*x_GAMMAR/(mixing_length_alpha*Hp)
         DR0 = convert(4d0*boltz_sigma*pow2(gammar_div_alfa)*pow3(T)/(pow2(rho)*Cp*kap))
      end if
      Pt0 = s% alpha_TDC_PtdVdt*x_ALFAP*rho
      if (s% dt > 0) then
         dVdt = dV/s% dt
      else
         dVdt = 0d0
      end if
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
   !! This function also has some side effects in terms of storing some of the terms
   !! it calculates for plotting purposes.
   !!
   !! @param s star pointer
   !! @param k face index
   !! @param A0_in convection speed from the start of the step (cm/s)
   !! @param xi0 The constant term in the convective velocity equation.
   !! @param xi1 The prefactor of the linear term in the convective velocity equation.
   !! @param xi2 The prefactor of the quadratic term in the convective velocity equation.            
   !! @param Af Output, the convection speed at the end of the step (cm/s)
   function eval_Af(s, k, A0_in, xi0, xi1, xi2) result(Af)
      type (star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1), intent(in) :: A0_in
      type(auto_diff_real_tdc), intent(in) :: xi0, xi1, xi2
      type(auto_diff_real_tdc) :: Af ! output
      type(auto_diff_real_tdc) :: J2, J, Jt, Jt4, num, den, y_for_atan, root, A0, lk 
      real(dp) :: dt    

      ! Debugging
      logical :: dbg = .false.
      integer :: dbg_k = 1448  

      include 'formats'

      J2 = pow2(xi1) - 4d0 * xi0 * xi2
      dt = s%dt
      A0 = convert(A0_in)

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
      if (k > 0) then ! save for plots
         s% SOURCE(k) = xi0%val*Af%val
         s% DAMPR(k) = -xi1%val*pow2(Af%val)
         s% DAMP(k) = -xi2%val*pow3(Af%val)
         s% COUPL(k) = s% SOURCE(k) - s% DAMP(k) - s% DAMPR(k)
      end if


      if (k == dbg_k .and. dbg) then
         write(*,*) J2%val,J2%d1val1
         write(*,*) Jt%val,Jt%d1val1
         write(*,*) xi0%val,xi0%d1val1
         write(*,*) xi1%val,xi1%d1val1
         write(*,*) xi2%val,xi2%d1val1
         write(*,*) num%val
         write(*,*) den%val
         write(*,*) Af%val,Af%d1val1
         write(*,*) ''
      end if


   end function eval_Af

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

end module tdc