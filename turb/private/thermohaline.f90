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

   use const_def, only: dp, no_mixing, thermohaline_mixing
   use num_lib
   use math_lib
   use chem_def, only: chem_isos
   use fingering_modes
   use parasite_model
   use turb_def

   implicit none

   private
   public :: get_thermohaline_info
   public :: set_info_HG19, set_info_FRG24 ! Used by plotter routines

contains

   !> Computes detialed information for thermohaline mixing --- when the
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
   !! @param gradL_composition_term B = (phi/delta)*dlnMu/dlnP where Mu is the mean molecular weight.
   !! @param iso The index of the species that drives thermohaline mixing.
   !! @param XH1 Mass fraction of H1.
   !! @param eta Magnetic diffusivity.
   !! @param thermohaline_coeff Free parameter multiplying the thermohaline diffusivity.
   !! @param thermohaline_mag_B Magnetic field strength (guass) for HG19 and FRG24 prescriptions.
   !! @param thermohaline_r_ext Reduced density ratio at which to switch to extrapolation for HG19 and FRG24 prescriptions.
   !! @param thermohaline_FRG24_safety Safety parameter for choosing approximations in FRG24 prescription.
   !! @param thermohaline_FRG24_nks Number of vertical wavenumbers to search over in FRG24 prescription.
   !! @param thermohaline_FRG24_N Maximal mode index in FRG24 prescription.
   !! @param th_info Output, detailed thermohaline info stored in th_info_t structure
   !! @param ierr Output, error index.
   subroutine get_thermohaline_info(thermohaline_option, &
      grada, gradr, N2_T, T, rho, Cp, opacity, &
      gradL_composition_term, XH1, eta, iso, &
      thermohaline_coeff, thermohaline_mag_B, thermohaline_r_ext, &
      thermohaline_FRG24_safety, thermohaline_FRG24_nks, thermohaline_FRG24_N, &
      th_info, ierr)

      character(*), intent(in)     :: thermohaline_option
      real(dp), intent(in)         :: grada
      real(dp), intent(in)         :: gradr
      real(dp), intent(in)         :: N2_T
      real(dp), intent(in)         :: T
      real(dp), intent(in)         :: rho
      real(dp), intent(in)         :: Cp
      real(dp), intent(in)         :: opacity
      real(dp), intent(in)         :: gradL_composition_term
      real(dp), intent(in)         :: XH1
      real(dp), intent(in)         :: eta
      integer, intent(in)          :: iso
      real(dp), intent(in)         :: thermohaline_coeff
      real(dp), intent(in)         :: thermohaline_mag_B
      real(dp), intent(in)         :: thermohaline_r_ext
      integer, intent(in)          :: thermohaline_FRG24_safety
      integer, intent(in)          :: thermohaline_FRG24_nks
      integer, intent(in)          :: thermohaline_FRG24_N
      type(th_info_t), intent(out) :: th_info
      integer, intent(out)         :: ierr
      
      include 'formats'

      ierr = 0

      ! Calculate common data

      call set_info_coeffs(T, rho, Cp, opacity, iso, XH1, eta, N2_T, thermohaline_mag_B, th_info)
      call set_info_strat(grada, gradr, gradL_composition_term, thermohaline_r_ext, th_info)

      ! Check for sensible Prandtl number
 
      if (th_info%Pr < 0._dp) then
         write(*, *) 'warning: get_thermohaline_info being called when Pr < 0'
         return
      end if

      ! Check for thermohaline instability based on r or r_prime

      th_info%mixing_type = no_mixing

      select case (thermohaline_option)
      case ('Kippenhahn')
      case ('Harrington_Garaud_19')
         if (th_info%r_prime >= 1._dp) return
      case ('Fraser_Reifenstein_Garaud_24')
         if (th_info%r_prime >= 1._dp) return
      case default
         if (th_info%r >= 1._dp) return
      end select

      th_info%mixing_type = thermohaline_mixing

      ! Dispatch to the various implementations
      
      select case (thermohaline_option)

      case ('Kippenhahn')

         call set_info_KRT80(grada, gradr, gradL_composition_term, rho, Cp, &
            th_info)

      case ('Traxler_Garaud_Stellmach_11')

         call set_info_TGS11(th_info)

      case ('Brown_Garaud_Stellmach_13')
         
         call set_info_BGS13(th_info, ierr)

      case ('Harrington_Garaud_19')

         call set_info_HG19(th_info, ierr)

      case ('Fraser_Reifenstein_Garaud_24')

         call set_info_FRG24(thermohaline_FRG24_safety, thermohaline_FRG24_nks, thermohaline_FRG24_N, &
            th_info, ierr)

      case default

         th_info%D_thrm = 0._dp
         ierr = -1
         write(*,*) 'unknown thermohaline_option' // trim(thermohaline_option)

      end select

      ! Rescale the thermohaline diffusivity

      th_info%D_thrm = thermohaline_coeff*th_info%D_thrm

   end subroutine get_thermohaline_info

   !****

   subroutine set_info_coeffs(T, rho, Cp, opacity, iso, XH1, eta, N2_T, B0, th_info)

      real(dp), intent(in)           :: T
      real(dp), intent(in)           :: rho
      real(dp), intent(in)           :: Cp
      real(dp), intent(in)           :: opacity
      integer, intent(in)            :: iso
      real(dp), intent(in)           :: XH1
      real(dp), intent(in)           :: eta
      real(dp), intent(in)           :: N2_T
      real(dp), intent(in)           :: B0
      type(th_info_t), intent(inout) :: th_info

      real(dp) :: K_therm, qe4, loglambdah, loglambdacx, loglambdacy, ccx, ccy
      real(dp) :: Bcoeff, chemA, chemZ, acx, acy, nu_mol, nu_rad
      real(dp) :: K_T, K_C, nu
      real(dp) :: d2

      ! Set fluid coefficients in th_info

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

      th_info%K_therm = K_therm
      th_info%K_T = K_T
      th_info%K_C = K_C
      th_info%nu = nu
      th_info%Pr = nu/K_T
      th_info%Pm = nu/eta
      th_info%tau = K_C/K_T
      th_info%H_B = B0*B0*d2/(pi4*rho*K_T*K_T)
      th_info%D_B = th_info%Pr/th_info%Pm

   end subroutine set_info_coeffs

   !****

   subroutine set_info_strat(grada, gradr, gradL_composition_term, r_ext, th_info)

      real(dp), intent(in)           :: grada
      real(dp), intent(in)           :: gradr
      real(dp), intent(in)           :: gradL_composition_term
      real(dp), intent(in)           :: r_ext
      type(th_info_t), intent(inout) :: th_info

      ! Set stratification coefficients in th_info

      th_info%R_0 = (gradr - grada)/gradL_composition_term
      th_info%r = (th_info%R_0 - 1._dp)/(1._dp/th_info%tau - 1._dp)

      th_info%r_prime = MIN(th_info%r, r_ext)
      th_info%R_0_prime = th_info%r_prime*(1._dp/th_info%tau - 1._dp) + 1._dp

!      print *,'set strat:',th_info%R_0, th_info%r, th_info%R_0_prime, th_info%r_prime, gradL_composition_term

   end subroutine set_info_strat
      
   !****

   subroutine set_info_KRT80(grada, gradr, gradL_composition_term, rho, Cp, &
      th_info)

      real(dp), intent(in)           :: grada
      real(dp), intent(in)           :: gradr
      real(dp), intent(in)           :: gradL_composition_term
      real(dp), intent(in)           :: rho
      real(dp), intent(in)           :: Cp
      type(th_info_t), intent(inout) :: th_info

      real(dp) :: dgrad

      ! Set components of th_info following Kippenhahn, R.,
      ! Ruschenplatt, G., & Thomas, H.-C. 1980, A&A, 91, 175 (KRT80)

      dgrad = max(1d-40, grada - gradr) ! this seems a bit janky, and is incompatible with set_info_strat

      th_info%D_thrm = -3._dp*th_info%K_therm/(2*rho*Cp)*gradL_composition_term/dgrad

   end subroutine set_info_KRT80

   !****

   subroutine set_info_TGS11(th_info)

      type(th_info_t), intent(inout) :: th_info
      
      ! Set components of th_info following Traxler, Garaud, &
      ! Stellmach, ApJ Letters, 728:L29 (2011). Also see
      ! Denissenkov. ApJ 723:563–579, 2010.

      th_info%D_thrm = 101._dp*sqrt(th_info%K_C*th_info%nu)* &
         exp(-3.6_dp*th_info%r)*pow(1._dp - th_info%r, 1.1_dp) ! eqn. (24)
      
   end subroutine set_info_TGS11

   !****

   subroutine set_info_BGS13(th_info, ierr)

      type(th_info_t), intent(inout) :: th_info
      integer, intent(out)           :: ierr
      
      ! Set components of th_info following Brown, Garaud, &
      ! Stellmach, ApJ 768:34 (2013)

      call eval_fastest_fingering(th_info%Pr, th_info%tau, th_info%R_0, th_info%lam_hat, th_info%l2_hat, ierr)
      if (ierr /= 0) return

      th_info%Nu_C = Nu_C_brown(th_info%tau, th_info%l2_hat, th_info%lam_hat)
      th_info%D_thrm = th_info%K_C*(th_info%Nu_C - 1._dp)

   end subroutine set_info_BGS13

   !****

   subroutine set_info_HG19(th_info, ierr)

      type(th_info_t), intent(inout) :: th_info
      integer, intent(out)           :: ierr

      real(dp), parameter :: K_B = 1.24_dp
      
      ! Set componets of th_info following Harrington & Garaud, ApJ
      ! Letters, 870:L5 (2019; HG19)

      call eval_fastest_fingering(th_info%Pr, th_info%tau, th_info%R_0_prime, th_info%lam_hat, th_info%l2_hat, ierr)
      if (ierr /= 0) return

      ! Solve for w_HG19

      call solve_HG19_eqn32(th_info%H_B, th_info%l2_hat, th_info%lam_hat, th_info%w_HG19, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in solve_HG19_eqn32'
         write(*,*) 'H_B', th_info%H_B
         write(*,*) 'l2_hat', th_info%l2_hat
         write(*,*) 'lam_hat', th_info%lam_hat
         write(*,*) 'w', th_info%w_HG19
         return
      end if

      th_info%w = th_info%w_HG19

      ! Evaluate Nu_C and D_thrm

      th_info%Nu_C = Nu_C(th_info%tau, th_info%w, th_info%lam_hat, th_info%l2_hat, K_B)
!      th_info%D_thrm = th_info%K_C*(th_info%Nu_C - 1._dp)*th_info%R_0_prime/th_info%R_0
      th_info%D_thrm = th_info%K_C*(th_info%Nu_C - 1._dp)*exp(-(th_info%r - th_info%r_prime))

   end subroutine set_info_HG19

   !****
  
   subroutine set_info_FRG24(safety, nks, N, th_info, ierr)

      integer, intent(in)            :: safety
      integer, intent(in)            :: nks
      integer, intent(in)            :: N
      type(th_info_t), intent(inout) :: th_info
      integer, intent(out)           :: ierr

      real(dp), parameter :: K_B = 0.62_dp

      real(dp), allocatable :: k_z(:)
      integer :: j
      
      ! Set components of th_info following Fraser, Reifenstein, &
      ! Garaud, ApJ 964:184 (2024; FRG24)

      call eval_fastest_fingering(th_info%Pr, th_info%tau, th_info%R_0_prime, th_info%lam_hat, th_info%l2_hat, ierr)
      if (ierr /= 0) return

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

      call eval_parasite_saturation(th_info%Pr, th_info%tau, th_info%R_0_prime, th_info%H_B, th_info%D_B, &
         th_info%lam_hat, th_info%l2_hat, k_z, N, safety, th_info%sigma_max, th_info%k_z_max, th_info%w_FRG24, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in eval_parasite_saturation'
         write(*,*) 'Pr', th_info%Pr
         write(*,*) 'tau', th_info%tau
         write(*,*) 'R_0', th_info%R_0
         write(*,*) 'R_0', th_info%R_0_prime
         write(*,*) 'H_B', th_info%H_B
         write(*,*) 'D_B', th_info%D_B
         write(*,*) 'l2_hat', th_info%l2_hat
         write(*,*) 'lam_hat', th_info%lam_hat
         write(*,*) 'w_HG19', th_info%w_HG19
         return
      end if

      ! For safety = 0, merge with w_HG19

      if (safety == 0) then
         call solve_HG19_eqn32(th_info%H_B, th_info%l2_hat, th_info%lam_hat, th_info%w_HG19, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in solve_HG19_eqn32'
            write(*,*) 'H_B', th_info%H_B
            write(*,*) 'l2_hat', th_info%l2_hat
            write(*,*) 'lam_hat', th_info%lam_hat
            write(*,*) 'w_HG19', th_info%w_HG19
            return
         end if
         th_info%w = MIN(th_info%w_FRG24, th_info%w_HG19)  ! TODO: is this a sane way to go about this?
      else
         th_info%w = th_info%w_FRG24
      end if

      ! Evaluate Nu_C and D_thrm

      th_info%Nu_C = Nu_C(th_info%tau, th_info%w, th_info%lam_hat, th_info%l2_hat, K_B)
      th_info%D_thrm = th_info%K_C*(th_info%Nu_C - 1._dp)*th_info%R_0_prime/th_info%R_0

   end subroutine set_info_FRG24

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
