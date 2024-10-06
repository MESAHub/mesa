! ***********************************************************************
!
!   Copyright (C) 2010-2024  The MESA Team
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

module thermohaline

   use const_def, only: dp
   use num_lib
   use math_lib
   use utils_lib
   use auto_diff
   use chem_def, only: chem_isos
   use fingering_modes
   use parasite_model

   implicit none

   ! Derived-type definitions

   type th_results_t
      real(dp) :: K_therm = 0._dp   ! Thermal conductivity
      real(dp) :: K_T = 0._dp       ! Thermal diffusivity
      real(dp) :: K_C = 0._dp       ! Chemical diffusivity
      real(dp) :: nu = 0._dp        ! Viscosity
      real(dp) :: Pr = 0._dp        ! Prandtl number
      real(dp) :: tau = 0._dp       ! Chemical diffusivity ratio
      real(dp) :: R_0 = 0._dp       ! Density ratio
      real(dp) :: r = 0._dp         ! Reduced density ratio
      real(dp) :: H_B = 0._dp       ! Lorentz force coefficient
      real(dp) :: Pm = 0._dp        ! Magnetic Prandtl number
      real(dp) :: D_B = 0._dp       ! Magnetic diffusivity
      real(dp) :: lam_hat = 0._dp   ! Growth rate of fastest-growing fingering mode
      real(dp) :: l2_hat = 0._dp    ! Horizontal wavenumber squared of fastest-growing fingering mode
      real(dp) :: sigma_max = 0._dp ! Growth rate of fastest-growing parasitic mode
      real(dp) :: k_z_max = 0._dp   ! Vertical wavenumber of fastest-growing parasitic mode
      real(dp) :: w = 0._dp         ! Saturation flow speed
      real(dp) :: w_HG19 = 0._dp    ! Saturation flow speed in HG19 treatment
      real(dp) :: w_FRG24 = 0._dp   ! Saturation flow speed in FRG24 treatment
      real(dp) :: nu_C = 0._dp      ! Compositional Nusselt number
      real(dp) :: D_thrm = 0._dp    ! Effective thermohaline mixing diffusivity
   end type th_results_t

   private
   public :: get_thermohaline_results
   public :: set_results_HG19, set_results_FRG24 ! Used by plotter routines
   public :: th_results_t

contains

   !> Computes detialed results for thermohaline mixing --- when the
   !! thermal gradient is stable and the composition gradient is unstable.
   !!
   !! @param thermohaline_option A string specifying which thermohaline prescription to use.
   !! @param grada Adiabatic gradient dlnT/dlnP
   !! @param gradr Radiative temperature gradient dlnT/dlnP, equals the actual gradient because there's no convection
   !! @param N2_T Structure part of brunt squared (excludes composition term)
   !! @param T Temperature
   !! @param opacity opacity
   !! @param rho Density
   !! @param Cp Heat capacity at constant pressure
   !! @param gradL_composition_term dlnMu/dlnP where Mu is the mean molecular weight.
   !! @param iso The index of the species that drives thermohaline mixing.
   !! @param XH1 Mass fraction of H1.
   !! @param eta Magnetic diffusivity.
   !! @param thermohaline_coeff Free parameter multiplying the thermohaline diffusivity.
   !! @param thermohaline_mag_B Magnetic field strength (guass) for HG19 and FRG24 prescriptions.
   !! @param thermohaline_FRG24_safety Safety parameter for choosing approximations in FRG24 prescription.
   !! @param thermohaline_FRG24_nks Number of vertical wavenumbers to search over in FRG24 prescription.
   !! @param thermohaline_FRG24_N Maximal mode index in FRG24 prescription.
   !! @param th_results Output, thermohaline results stored in th_results_t structure
   !! @param ierr Output, error index.
   subroutine get_thermohaline_results(thermohaline_option, &
      grada, gradr, N2_T, T, rho, Cp, opacity, &
      gradL_composition_term, XH1, eta, iso, &
      thermohaline_coeff, thermohaline_mag_B, thermohaline_FRG24_safety, thermohaline_FRG24_nks, thermohaline_FRG24_N, &
      th_results, ierr)

      character(*), intent(in)        :: thermohaline_option
      real(dp), intent(in)            :: grada
      real(dp), intent(in)            :: gradr
      real(dp), intent(in)            :: N2_T
      real(dp), intent(in)            :: T
      real(dp), intent(in)            :: rho
      real(dp), intent(in)            :: Cp
      real(dp), intent(in)            :: opacity
      real(dp), intent(in)            :: gradL_composition_term
      real(dp), intent(in)            :: XH1
      real(dp), intent(in)            :: eta
      integer, intent(in)             :: iso
      real(dp), intent(in)            :: thermohaline_coeff
      real(dp), intent(in)            :: thermohaline_mag_B
      integer, intent(in)             :: thermohaline_FRG24_safety
      integer, intent(in)             :: thermohaline_FRG24_nks
      integer, intent(in)             :: thermohaline_FRG24_N
      type(th_results_t), intent(out) :: th_results
      integer, intent(out)            :: ierr
      
      include 'formats'

      ierr = 0

      ! Calculate common data

      call set_results_coeffs(T, rho, Cp, opacity, iso, XH1, eta, N2_T, thermohaline_mag_B, th_results)
      call set_results_strat(grada, gradr, gradL_composition_term, th_results)

      ! Handle cases where this routine shouldn't have been called in the first place
 
      if (th_results%Pr < 0._dp) then
         write(*, *) 'warning: get_thermohaline_results being called when Pr < 0'
         return
      else if (th_results%r > 1._dp) then
         write(*, *) 'warning: get_thermohaline_results being called when r > 1'
         return
      end if

      ! Dispatch to the various implementations
      
      select case (thermohaline_option)

      case ('Kippenhahn')

         call set_results_KRT80(grada, gradr, gradL_composition_term, rho, Cp, &
            th_results)

      case ('Traxler_Garaud_Stellmach_11')

         call set_results_TGS11(th_results)

      case ('Brown_Garaud_Stellmach_13')
         
         call set_results_BGS13(th_results, ierr)

      case ('Harrington_Garaud_19')

         call set_results_HG19(th_results, ierr)

      case ('Fraser_Reifenstein_Garaud_24')

         call set_results_FRG24(thermohaline_FRG24_safety, thermohaline_FRG24_nks, thermohaline_FRG24_N, &
            th_results, ierr)

      case default

         th_results%D_thrm = 0._dp
         ierr = -1
         write(*,*) 'unknown for MLT thermohaline_option' // trim(thermohaline_option)

      end select

      ! Rescale the thermohaline diffusivity

      th_results%D_thrm = thermohaline_coeff*th_results%D_thrm

   end subroutine get_thermohaline_results

   !****

   subroutine set_results_coeffs(T, rho, Cp, opacity, iso, XH1, eta, N2_T, B0, th_results)

      real(dp), intent(in)              :: T
      real(dp), intent(in)              :: rho
      real(dp), intent(in)              :: Cp
      real(dp), intent(in)              :: opacity
      integer, intent(in)               :: iso
      real(dp), intent(in)              :: XH1
      real(dp), intent(in)              :: eta
      real(dp), intent(in)              :: N2_T
      real(dp), intent(in)              :: B0
      type(th_results_t), intent(inout) :: th_results

      real(dp) :: K_therm, qe4, loglambdah, loglambdacx, loglambdacy, ccx, ccy
      real(dp) :: Bcoeff, chemA, chemZ, acx, acy, nu_mol, nu_rad
      real(dp) :: K_T, K_C, nu
      real(dp) :: d2

      ! Set fluid coefficients in th_results

      K_therm = 4._dp*crad*clight*pow3(T)/(3._dp*opacity*rho) ! thermal conductivity

      K_T = K_therm/(Cp*rho) ! thermal diffusivity (assumes radiatively dominated)
      qe4=pow4(qe)

      ! Log Lambda for pure H (equation 10 from Proffitt Michaud 93)

      loglambdah = -19.26d0 - 0.5d0*log(rho) + 1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+XH1)) 

      ! From Spitzer "Physics of Fully Ionized Gases equation 5-54
      ! Assumes pure H. Still trying to work out what it would be for a mixture. 

      nu_rad = 4d0*crad*pow4(T)/(15d0*clight*opacity*pow2(rho)) ! radiative viscosity
      nu_mol = 0.406d0*sqrt(amu)*pow(boltzm*T,2.5d0)/(qe4*loglambdah*rho) 
      nu = nu_mol + nu_rad   ! total viscosity

      ! The following is from Proffitt & Michaud, 1993.
      ! Their constant B (equation 15)

      Bcoeff = (15.d0/16.d0)*sqrt(2.d0*amu/(5*pi))*pow(boltzm,2.5d0)/qe4

      ! Extract what species drives the thermohaline concvection

      chemA = chem_isos%Z_plus_N(iso)
      chemZ = chem_isos%Z(iso)

      if (chemZ.gt.2) then

         ! This is if the driving chemical is NOT He.
         ! Log Lambda for H-dominant chem mixture (equation 10)

         loglambdacx = loglambdah - log(chemz)  

         ! Log Lambda for He-dominant chem mixture (equation 10)

         loglambdacy = loglambdah - log(2.d0*chemz)

         ! Calculation of C_ij coeffs (equation 12)

         ccx = log(exp(1.2d0*loglambdacx)+1.)/1.2d0
         ccy = log(exp(1.2d0*loglambdacy)+1.)/1.2d0

         ! Reduced masses (I had to guess, from Bahcall & Loeb 1990), with H and He

         acx = (1.d0*chemA)/(1.d0+chemA)
         acy = 4*chemA/(4.d0+chemA)

         ! My formula (see notes) based on Proffitt and Michaud 1993
         K_C = 2*Bcoeff*pow(T,2.5d0)/(sqrt(5._dp)*rho*chemZ*chemZ)/ &
            (XH1*sqrt(acx)*ccx + (1-XH1)*sqrt(acy)*ccy)

      else
         
         ! Log Lambda for H-He mixture (equation 10)

         loglambdah = -19.26d0 - log(2d0) - 0.5d0*log(rho) + &
            1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+XH1)) 

         ! Calculation of C_ij coeffs (equation 12)

         ccy = log(exp(1.2d0*loglambdah)+1d0)/1.2d0

         ! My formula (see notes) based on Proffitt and Michaud 1993

         K_C = (Bcoeff*pow(T,2.5d0)/(rho*ccy))*(3+XH1)/((1+XH1)*(3+5*XH1)*(0.7d0+0.3d0*XH1))
      
      endif

      ! Magnetic stuff

      d2 = pow(K_T*nu/N2_T,0.5_dp) ! width of fingers squared

      ! Store results

      th_results%K_therm = K_therm
      th_results%K_T = K_T
      th_results%K_C = K_C
      th_results%nu = nu
      th_results%Pr = nu/K_T
      th_results%Pm = nu/eta
      th_results%tau = K_C/K_T
      th_results%H_B = B0*B0*d2/(pi4*rho*K_T*K_T)
      th_results%D_B = th_results%Pr/th_results%Pm

   end subroutine set_results_coeffs

   !****

   subroutine set_results_strat(grada, gradr, gradL_composition_term, th_results)

      real(dp), intent(in) :: grada, gradr, gradL_composition_term
      type(th_results_t), intent(inout) :: th_results

      ! Set stratification coefficients in th_results

      th_results%R_0 = (gradr - grada)/gradL_composition_term
      th_results%r = (th_results%R_0 - 1._dp)/(1._dp/th_results%tau - 1._dp)

   end subroutine set_results_strat
      
   !****

   subroutine set_results_KRT80(grada, gradr, gradL_composition_term, rho, Cp, &
      th_results)

      real(dp), intent(in) :: grada, gradr, gradL_composition_term, rho, Cp
      type(th_results_t), intent(inout) :: th_results

      real(dp) :: dgrad

      ! Set components of th_results following Kippenhahn, R.,
      ! Ruschenplatt, G., & Thomas, H.-C. 1980, A&A, 91, 175 (KRT80)

      dgrad = max(1d-40, grada - gradr) ! this seems a bit janky, and is incompatible with set_results_strat

      th_results%D_thrm = -3._dp*th_results%K_therm/(2*rho*Cp)*gradL_composition_term/dgrad

   end subroutine set_results_KRT80

   !****

   subroutine set_results_TGS11(th_results)

      type(th_results_t), intent(inout) :: th_results
      
      ! Set components of th_results following Traxler, Garaud, &
      ! Stellmach, ApJ Letters, 728:L29 (2011). Also see
      ! Denissenkov. ApJ 723:563–579, 2010.

      th_results%D_thrm = 101._dp*sqrt(th_results%K_C*th_results%nu)* &
         exp(-3.6_dp*th_results%r)*pow(1._dp - th_results%r, 1.1_dp) ! eqn. (24)
      
   end subroutine set_results_TGS11

   !****

   subroutine set_results_BGS13(th_results, ierr)

      type(th_results_t), intent(inout) :: th_results
      integer, intent(out) :: ierr
      
      ! Set components of th_results following Brown, Garaud, &
      ! Stellmach, ApJ 768:34 (2013)

      call eval_fastest_fingering(th_results%Pr, th_results%tau, th_results%R_0, th_results%lam_hat, th_results%l2_hat, ierr)
      if (ierr /= 0) return

      th_results%nu_C = Nu_C_brown(th_results%tau, th_results%l2_hat, th_results%lam_hat)
      th_results%D_thrm = th_results%K_C*(th_results%Nu_C - 1._dp)

   end subroutine set_results_BGS13

   !****

   subroutine set_results_HG19(th_results, ierr)

      type(th_results_t), intent(inout) :: th_results
      integer, intent(out) :: ierr

      real(dp), parameter :: K_B = 1.24
      
      ! Set componets of th_results following Harrington & Garaud, ApJ
      ! Letters, 870:L5 (2019; HG19)

      call eval_fastest_fingering(th_results%Pr, th_results%tau, th_results%R_0, th_results%lam_hat, th_results%l2_hat, ierr)
      if (ierr /= 0) return

      ! Solve for w_HG19

      call solve_HG19_eqn32(th_results%H_B, th_results%l2_hat, th_results%lam_hat, th_results%w_HG19, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in solve_HG19_eqn32'
         write(*,*) 'H_B', th_results%H_B
         write(*,*) 'l2_hat', th_results%l2_hat
         write(*,*) 'lam_hat', th_results%lam_hat
         write(*,*) 'w', th_results%w_HG19
         return
      end if

      th_results%w = th_results%w_HG19

      ! Evaluate Nu_C and D_thrm

      th_results%Nu_C = Nu_C(th_results%tau, th_results%w, th_results%lam_hat, th_results%l2_hat, K_B)
      th_results%D_thrm = th_results%K_C*(th_results%Nu_C - 1._dp)

   end subroutine set_results_HG19

   !****
  
   subroutine set_results_FRG24(safety, nks, N, th_results, ierr)

      integer, intent(in) :: safety, nks, N
      type(th_results_t), intent(inout) :: th_results
      integer, intent(out) :: ierr

      real(dp), parameter :: K_B = 0.62_dp

      real(dp), allocatable :: k_z(:)
      integer :: j
      
      ! Set components of th_results following Fraser, Reifenstein, &
      ! Garaud, ApJ 964:184 (2024; FRG24)

      ! Define grid of vertical wavenumbers. This may evolve. Rich is
      ! working on optimization TODO: I don't think I've ever seen
      ! k/lhat > 1 being relevant

      allocate(k_z(2*nks))

      do j = 1, nks
         ! first nks entries are log space from 1e-6 to 0.1 (don't include endpoint)
         k_z(j) = pow(10._dp, -6._dp + (j-1)*(-1._dp + 6._dp)/nks)
         ! last nks entries are linear space from 0.1 to 2
         k_z(j+nks) = 0.1_dp + (j-1)*(2._dp - 0.1_dp)/(nks - 1)
      end do

      ! Solve for w_FRG24

      call eval_parasite_saturation(th_results%Pr, th_results%tau, th_results%R_0, th_results%H_B, th_results%D_B, &
         th_results%lam_hat, th_results%l2_hat, k_z, N, safety, th_results%sigma_max, th_results%k_z_max, th_results%w_FRG24, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in eval_parasite_saturation'
         write(*,*) 'Pr', th_results%Pr
         write(*,*) 'tau', th_results%tau
         write(*,*) 'R_0', th_results%R_0
         write(*,*) 'H_B', th_results%H_B
         write(*,*) 'D_B', th_results%D_B
         write(*,*) 'l2_hat', th_results%l2_hat
         write(*,*) 'lam_hat', th_results%lam_hat
         write(*,*) 'w_HG19', th_results%w_HG19
         return
      end if

      ! For safety = 0, merge with w_HG19

      if (safety == 0) then
         call solve_HG19_eqn32(th_results%H_B, th_results%l2_hat, th_results%lam_hat, th_results%w_HG19, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in solve_HG19_eqn32'
            write(*,*) 'H_B', th_results%H_B
            write(*,*) 'l2_hat', th_results%l2_hat
            write(*,*) 'lam_hat', th_results%lam_hat
            write(*,*) 'w_HG19', th_results%w_HG19
            return
         end if
         th_results%w = MIN(th_results%w_FRG24, th_results%w_HG19)  ! TODO: is this a sane way to go about this?
      else
         th_results%w = th_results%w_FRG24
      end if

      ! Evaluate Nu_C and D_thrm

      th_results%Nu_C = Nu_C(th_results%tau, th_results%w, th_results%lam_hat, th_results%l2_hat, K_B)
      th_results%D_thrm = th_results%K_C*(th_results%Nu_C - 1._dp)

   end subroutine set_results_FRG24

   !****

   function Nu_C_brown(tau, l2_hat, lam_hat) result(Nu_C)

      real(dp), intent(in) :: tau
      real(dp), intent(in) :: l2_hat
      real(dp), intent(in) :: lam_hat
      real(dp)             :: Nu_C

      ! The Nusselt number Nu_C reduces to this simpler form for the
      ! Brown model, Formula (33) from Brown et al, with C = 7.

      Nu_C = 1.d0 + 49.d0*lam_hat*lam_hat/(tau*l2_hat*(lam_hat+tau*l2_hat))

   end function Nu_C_brown

   !****

   function Nu_C(tau, w, lam_hat, l2_hat, KB)

      real(dp), intent(in) :: tau
      real(dp), intent(in) :: w
      real(dp), intent(in) :: lam_hat
      real(dp), intent(in) :: l2_hat
      real(dp), intent(in) :: KB
      real(dp)             :: Nu_C

      ! More general expression for the Nusselt number Nu_C in terms
      ! of w, for Harrington and Fraser models KB = 1.24 for
      ! Harrington, 0.62 for Fraser

      Nu_C =  1d0 + KB * pow2(w) / (tau * (lam_hat + tau * l2_hat))

   end function Nu_C
   
   !****

   ! Solver for HG19's eqn. 32

   subroutine solve_HG19_eqn32(H_B, l2_hat, lam_hat, w, ierr)

      real(dp), intent(in)           :: H_B
      real(dp), intent(in)           :: l2_hat
      real(dp), intent(in)           :: lam_hat
      real(dp), intent(out)          :: w
      integer, intent(out)           :: ierr

      integer, parameter  :: NEWT_IMAX = 100
      integer, parameter  :: IMAX = 0
      real(dp), parameter :: EPSX = 2e-12_dp
      real(dp), parameter :: EPSY = 0._dp
      real(dp), parameter :: CH = 1.66_dp

      real(dp)          :: w0
      real(dp)          :: l_hat
      real(dp), pointer :: rpar(:) => null() ! not used, but needed to pass to safe_root_with_brackets
      integer, pointer  :: ipar(:) => null() ! not used, but needed to pass to safe_root_with_brackets
      integer :: lrpar = 0, lipar = 0        ! not used, but needed to pass to safe_root_with_brackets
      
      l_hat = sqrt(l2_hat)

      w0 = MAX(sqrt(2.0_dp*H_B), 2.0_dp * PI * lam_hat/l_hat)

      w = safe_root_with_guess(root_func_, w0, 0.5_dp*w0, &
           arg_not_provided, arg_not_provided, &
           arg_not_provided, arg_not_provided, &
           NEWT_IMAX, IMAX, EPSX, EPSY, &
           lrpar, rpar, lipar, ipar, &
           ierr)

      return

   contains

      function root_func_(x, df_dx, lrpar, rpar, lipar, ipar, ierr) result(f)

         real(dp), intent(in)             :: x
         real(dp), intent(out)            :: df_dx
         integer, intent(in)              :: lrpar
         real(dp), intent(inout), pointer :: rpar(:)
         integer, intent(in)              :: lipar
         integer, intent(inout), pointer  :: ipar(:)
         integer, intent(out)             :: ierr
         real(dp)                         :: f

         real(dp) :: LHS
         real(dp) :: RHS1
         real(dp) :: RHS2

         ! Evaluate eqn. 32 of HG19

         associate(w => x)

            LHS = 0.5_dp * w**2.0_dp - H_B
            RHS1 = (CH*lam_hat/(0.42_dp*l_hat))**1.5_dp
            RHS2 = sqrt(w)

            f = LHS - (RHS1*RHS2)
            df_dx = w - 0.5_dp * (CH*lam_hat/(0.42_dp*l_hat))**1.5_dp / sqrt(w)         

         end associate

         ierr = 0

      end function root_func_

   end subroutine solve_hg19_eqn32

         
         
         
   ! end subroutine solve_hg19_eqn32

   ! !****

   ! ! RHS function for HG19's eqn. 32 (function and derivatives, suitable for use
   ! ! with safe_root_with_guess)

   ! real(dp) function hg19_eqn32_func(x, df_dx, lrpar, rpar, lipar, ipar, ierr) result(f)

   !    use const_def, only: dp

   !    real(dp), intent(in)             :: x
   !    real(dp), intent(out)            :: df_dx
   !    integer, intent(in)              :: lrpar
   !    real(dp), intent(inout), pointer :: rpar(:)
   !    integer, intent(in)              :: lipar
   !    integer, intent(inout), pointer  :: ipar(:)
   !    integer, intent(out)             :: ierr

   !    real(dp) :: LHS
   !    real(dp) :: RHS1
   !    real(dp) :: RHS2

   !    ! Evaluate eqn. 32 of HG19

   !    if (lrpar /= 5 .OR. lipar /= 0) then
   !       ierr = -1
   !       return
   !    end if

   !    associate( &
   !         w => x, &
   !         l2_hat => rpar(1), &
   !         lhat => rpar(2), &
   !         lam_hat => rpar(3), &
   !         H_B => rpar(4), &
   !         CH => rpar(5))

   !      LHS = 0.5_dp * w**2.0_dp - H_B
   !      RHS1 = (CH*lam_hat/(0.42_dp*lhat))**1.5_dp
   !      RHS2 = sqrt(w)

   !      f = LHS - (RHS1*RHS2)
   !      df_dx = w - 0.5_dp * (CH*lam_hat/(0.42_dp*lhat))**1.5_dp / sqrt(w)         

   !    end associate

   !    ierr = 0

   !    return

   ! end function hg19_eqn32_func

end module thermohaline
