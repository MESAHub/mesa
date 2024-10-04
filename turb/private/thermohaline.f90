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


module thermohaline

use const_def
use num_lib
use utils_lib
use auto_diff
use fingering_modes
use parasite_model

implicit none

! Indices for thrm_extras

integer, parameter :: I_PR = 1
integer, parameter :: I_TAU = 2
integer, parameter :: I_R0 = 3
integer, parameter :: I_DB = 4
integer, parameter :: I_HB = 5
integer, parameter :: I_LAMHAT = 6
integer, parameter :: I_L2HAT = 7
integer, parameter :: I_W = 8
integer, parameter :: I_W_TC = 9
integer, parameter :: I_W_HG19 = 10
integer, parameter :: I_D_THRM = 11

integer, parameter :: N_THRM_EXTRAS = 11

private
public :: get_D_thermohaline, nuC, solve_hg19_eqn32
public :: I_PR, I_TAU, I_R0, I_DB, I_HB, I_LAMHAT, I_L2HAT, &
          I_W, I_W_TC, I_W_HG19, I_D_THRM, N_THRM_EXTRAS

contains

   !> Computes the diffusivity of thermohaline mixing when the
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
   !! @param thermohaline_coeff Free parameter multiplying the thermohaline diffusivity.
   !! @param thermohaline_mag_B Magnetic field strength (guass) for HG19 and FRG24 prescriptions.
   !! @param thermohaline_FRG24_safety Safety parameter for choosing approximations in FRG24 prescription.
   !! @param thermohaline_FRG24_nks Number of vertical wavenumbers in FRG24 prescription.
   !! @param thermohaline_FRG24_res Resolution in FRG24 prescription.
   !! @param thermohaline_FRG24_with_TC Include temperature & composition in FRG24 prescription.
   !! @param D_thrm Output, diffusivity.
   !! @param thrm_extras Output, debugging data (ignored if null)
   !! @param ierr Output, error index.
   subroutine get_D_thermohaline(thermohaline_option, &
         grada, gradr, N2_T, T, opacity, rho, Cp, gradL_composition_term, &
         iso, XH1, thermohaline_coeff, thermohaline_mag_B, &
         thermohaline_FRG24_safety, thermohaline_FRG24_nks, &
         thermohaline_FRG24_res, thermohaline_FRG24_with_TC, &
         eta, D_thrm, thrm_extras, ierr)
      character(len=*) :: thermohaline_option
      real(dp), intent(in) :: &
         grada, gradr, N2_T, T, opacity, rho, Cp, gradL_composition_term, XH1, &
         thermohaline_coeff, thermohaline_mag_B, eta
      integer, intent(in) :: iso, thermohaline_FRG24_safety, thermohaline_FRG24_nks, thermohaline_FRG24_res
      logical, intent(in) :: thermohaline_FRG24_with_TC
      real(dp), intent(out) :: D_thrm
      real(dp), pointer, intent(in) :: thrm_extras(:)
      integer, intent(out) :: ierr
      real(dp) :: dgrad, K_therm, K_T, K_mu, nu, R0, Pr, tau, r_th, Pm, DB
      real(dp) :: l2hat, lamhat, w, w_TC, w_HG19, d2, HB, B0, l2hat_test, lamhat_test, reldiff, reldiff2
      real(dp), allocatable :: ks(:)
      integer :: j, nks, spectral_resolution, safety ! for FRG model
      logical :: withTC
      logical, parameter :: dbg = .true.
      
      include 'formats'     
      dgrad = max(1d-40, grada - gradr) ! positive since Schwarzschild stable
      K_therm = 4d0*crad*clight*pow3(T)/(3d0*opacity*rho) ! thermal conductivity

      if (ASSOCIATED(thrm_extras)) thrm_extras = 0._dp

      if (thermohaline_option == 'Kippenhahn') then
         ! Kippenhahn, R., Ruschenplatt, G., & Thomas, H.-C. 1980, A&A, 91, 175
         D_thrm = -3d0*K_therm/(2*rho*cp)*gradL_composition_term/dgrad
      else if (thermohaline_option == 'Traxler_Garaud_Stellmach_11' .or. &
               thermohaline_option == 'Brown_Garaud_Stellmach_13' .or. &
               thermohaline_option == 'Harrington_Garaud_19' .or. &
               thermohaline_option == 'Fraser_Reifenstein_Garaud_24') then
         call get_diff_coeffs(K_therm, Cp, rho, T, opacity, iso, XH1, K_T, K_mu, nu)
         R0 = (gradr - grada)/gradL_composition_term
         Pr = nu/K_T
         tau = K_mu/K_T

         if (ASSOCIATED(thrm_extras)) then
            thrm_extras(I_PR) = Pr
            thrm_extras(I_TAU) = tau
            thrm_extras(I_R0) = R0
         end if

         r_th = (R0 - 1d0)/(1d0/tau - 1d0)
         if (r_th >= 1d0) then ! stable if R0 >= 1/tau
            D_thrm = 0d0
         else if (Pr < 0d0) then
            ! Bad results from get_diff_coeffs will just result in NaNs from thermohaline options, so skip
            D_thrm = 0d0
         else if (thermohaline_option == 'Traxler_Garaud_Stellmach_11') then 
            ! Traxler, Garaud, & Stellmach, ApJ Letters, 728:L29 (2011).
            ! also see Denissenkov. ApJ 723:563–579, 2010.
            D_thrm = 101d0*sqrt(K_mu*nu)*exp(-3.6d0*r_th)*pow(1d0 - r_th,1.1d0) ! eqn 24
         else if (thermohaline_option == 'Brown_Garaud_Stellmach_13') then
            call gaml2max(pr, tau, R0, lamhat, l2hat, ierr)
            D_thrm = K_mu*(nuC_brown(tau, l2hat, lamhat) - 1d0)

            if (ASSOCIATED(thrm_extras)) then
               thrm_extras(I_LAMHAT) = lamhat
               thrm_extras(I_L2HAT) = l2hat
            endif

         else if (thermohaline_option == 'Harrington_Garaud_19') then
            call gaml2max(pr, tau, R0, lamhat, l2hat, ierr)

            B0 = thermohaline_mag_B
            d2 = pow(K_T*nu/N2_T,0.5d0) ! width of fingers squared
            HB = B0*B0*d2/(pi4*rho*K_T*K_T)

            ! solve for w based on Harrington & Garaud model
            call solve_hg19_eqn32(HB,l2hat,lamhat,w,ierr)
            if(ierr /= 0 .and. dbg) then
               write(*,*) "failed in solve_hg19_eqn32"
               write(*,*) "HB", HB
               write(*,*) "l2hat", l2hat
               write(*,*) "lamhat", lamhat
               write(*,*) "w", w
            end if

            ! KB = 1.24 for Harrington model.
            D_thrm = K_mu*(nuC(tau, w, lamhat, l2hat, 1.24d0) - 1d0)

            if (ASSOCIATED(thrm_extras)) then
               thrm_extras(I_HB) = HB
               thrm_extras(I_LAMHAT) = lamhat
               thrm_extras(I_L2HAT) = l2hat
            endif

         else if (thermohaline_option == 'Fraser_Reifenstein_Garaud_24') then
            call gaml2max(pr, tau, R0, lamhat, l2hat, ierr)

            ! TODO: promote these to inlist options [done]
            nks = thermohaline_FRG24_nks
            spectral_resolution = thermohaline_FRG24_res
            withTC = thermohaline_FRG24_with_TC
            safety = thermohaline_FRG24_safety
            
            B0 = thermohaline_mag_B
            d2 = pow(K_T*nu/N2_T,0.5d0) ! width of fingers squared
            HB = B0*B0*d2/(pi4*rho*K_T*K_T)

            Pm = nu/eta ! Magnetic Prandtl number. TODO: calculate self-consistently from plasma conditions [DONE!]
            DB = Pr/Pm

            ! This may evolve. Rich is working on optimization
            ! TODO: I don't think I've ever seen k/lhat > 1 being relevant
            allocate(ks(nks*2))
            do j = 1,nks
               ! first nks entries are log space from 1e-6 to 0.1 (don't include endpoint)
               ks(j) = pow(10d0,-6d0 + (j-1)*(-1d0 + 6d0)/nks)
               ! last nks entries are linear space from 0.1 to 2
               ks(j+nks) = 0.1d0 + (j-1)*(2d0 - 0.1d0)/(nks - 1)
            end do
            
            ! solve for w based on Fraser model
            if(withTC) then
               w_TC = wf_withTC(pr, tau, R0, HB, DB, ks, spectral_resolution, .false., safety, lamhat, l2hat)
               if (safety == 0) then
                  call solve_hg19_eqn32(HB,l2hat,lamhat,w_HG19,ierr)
                  if(ierr /= 0 .and. dbg) then
                     write(*,*) "failed in solve_hg19_eqn32"
                     write(*,*) "HB", HB
                     write(*,*) "l2hat", l2hat
                     write(*,*) "lamhat", lamhat
                     write(*,*) "w_TC", w_TC
                  end if
                  w = MIN(w_TC, w_HG19)  ! TODO: is this a sane way to go about this?
               end if
            else  ! TODO: we really should just ditch everything that isn't withTC, in hindsight
               w_TC = 0._dp
               w_HG19 = 0._dp
               w = wf(pr, tau, R0, HB, DB, ks, spectral_resolution, 0d0, 0, .false., .false., lamhat, l2hat)
            end if

            ! KB = 0.62 for Fraser model.
            D_thrm = K_mu*(nuC(tau, w, lamhat, l2hat, 0.62d0) - 1d0)

            deallocate(ks)

            if (ASSOCIATED(thrm_extras)) then
               thrm_extras(I_HB) = HB
               thrm_extras(I_DB) = DB
               thrm_extras(I_LAMHAT) = lamhat
               thrm_extras(I_L2HAT) = l2hat
               thrm_extras(I_W) = w
               thrm_extras(I_W_TC) = w_TC
               thrm_extras(I_W_HG19) = w_HG19
            end if

         endif

         if (ASSOCIATED(thrm_extras)) then
            thrm_extras(I_D_THRM) = D_thrm
         endif

      else
         D_thrm = 0
         ierr = -1
         write(*,*) 'unknown for MLT thermohaline_option' // trim(thermohaline_option)
      end if
      D_thrm = thermohaline_coeff*D_thrm
   end subroutine get_D_thermohaline


   subroutine get_diff_coeffs(K_therm, Cp, rho, T, opacity, iso, XH1, kt, kmu, vis)
      use chem_def, only: chem_isos
      real(dp), intent(in) :: K_therm, Cp, rho, T, opacity, XH1
      integer, intent(in) :: iso      
      real(dp), intent(out) :: kt, kmu, vis
      real(dp) :: loglambdah, loglambdacx, loglambdacy, ccx, ccy, qe4
      real(dp) :: Bcoeff, chemA, chemZ, acx, acy, nu_mol, nu_rad      
      real(dp), parameter :: sqrt5 = sqrt(5d0)           
      kt = K_therm/(Cp*rho)       ! thermal diffusivity (assumes radiatively dominated)
      qe4=pow4(qe)

      ! Log Lambda for pure H (equation 10 from Proffitt Michaud 93)
      loglambdah = -19.26d0 - 0.5d0*log(rho) + 1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+XH1)) 
      nu_rad = 4d0*crad*pow4(T)/(15d0*clight*opacity*pow2(rho)) ! radiative viscosity
      nu_mol = 0.406d0*sqrt(amu)*pow(boltzm*T,2.5d0)/(qe4*loglambdah*rho) 
      ! From Spitzer "Physics of Fully Ionized Gases equation 5-54
      ! Assumes pure H. Still trying to work out what it would be for a mixture. 
      vis = nu_mol + nu_rad   ! total viscosity

      ! The following is from Proffitt & Michaud, 1993.
      ! Their constant B (equation 15)
      Bcoeff = (15.d0/16.d0)*sqrt(2.d0*amu/(5*pi))*pow(boltzm,2.5d0)/qe4
      ! Extract what species drives the thermohaline concvection
      chemA = chem_isos%Z_plus_N(iso)
      chemZ = chem_isos%Z(iso)

      if(chemZ.gt.2) then
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
         kmu = 2*Bcoeff*pow(T,2.5d0)/(sqrt5*rho*chemZ*chemZ)/ &
            (XH1*sqrt(acx)*ccx + (1-XH1)*sqrt(acy)*ccy)

      else
         ! Log Lambda for H-He mixture (equation 10)
         loglambdah = -19.26d0 - log(2d0) - 0.5d0*log(rho) + &
            1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+XH1)) 
         ! Calculation of C_ij coeffs (equation 12)
         ccy = log(exp(1.2d0*loglambdah)+1d0)/1.2d0
         ! My formula (see notes) based on Proffitt and Michaud 1993
         kmu = (Bcoeff*pow(T,2.5d0)/(rho*ccy))*(3+XH1)/((1+XH1)*(3+5*XH1)*(0.7d0+0.3d0*XH1))
      
      endif
      ! write(57,*) kt,kmu,vis,chemZ

   end subroutine get_diff_coeffs

   real(dp) function nuC_brown(tau, l2hat, lamhat)
      real(dp), intent(in) :: tau, l2hat, lamhat
      
      ! Nu_C reduces to this simpler form for the Brown model,
      ! Formula (33) from Brown et al, with C = 7.
      nuC_brown = 1.d0 + 49.d0*lamhat*lamhat/(tau*l2hat*(lamhat+tau*l2hat))
      return
   end function nuC_brown

   ! More general expression for Nu_C in terms of w, for Harrington and Fraser models
   ! KB = 1.24 for Harrington, 0.62 for Fraser
   real(dp) function nuC(tau, w, lamhat, l2hat, KB)
      real(dp), intent(in) :: tau, w, lamhat, l2hat, KB
      nuC =  1d0 + KB * pow2(w) / (tau * (lamhat + tau * l2hat))
      return
   end function nuC
   
   !****

   ! Solver for HG19's eqn. 32

   subroutine solve_hg19_eqn32(HB, l2hat, lamhat, w, ierr, CH)

      real(dp), intent(in)           :: HB
      real(dp), intent(in)           :: l2hat
      real(dp), intent(in)           :: lamhat
      real(dp), intent(out)          :: w
      integer, intent(out)           :: ierr
      real(dp), intent(in), optional :: CH

      integer, parameter  :: NEWT_IMAX = 100
      integer, parameter  :: IMAX = 0
      real(dp), parameter :: EPSX = 2e-12_dp
      real(dp), parameter :: EPSY = 0._dp

      real(dp)          :: CH_
      real(dp)          :: w0
      real(dp)          :: lhat
      real(dp), target  :: rpar(5)
      integer,  target  :: ipar(0)
      real(dp), pointer :: rpar_(:)
      integer, pointer  :: ipar_(:)

      if (PRESENT(CH)) then
         CH_ = CH
      else
         CH_ = 1.66_dp
      end if

      lhat = sqrt(l2hat)

      w0 = MAX(sqrt(2.0_dp*HB), 2.0_dp * PI * lamhat/lhat)

      rpar = [l2hat, lhat, lamhat, HB, CH_]

      rpar_ => rpar
      ipar_ => ipar

      w = safe_root_with_guess(hg19_eqn32_func, w0, 0.5*w0, &
           arg_not_provided, arg_not_provided, &
           arg_not_provided, arg_not_provided, &
           NEWT_IMAX, IMAX, EPSX, EPSY, &
           SIZE(rpar_), rpar_, SIZE(ipar_), ipar_, &
           ierr)

      return

   end subroutine solve_hg19_eqn32

   !****

   ! RHS function for HG19's eqn. 32 (function and derivatives, suitable for use
   ! with safe_root_with_guess)

   real(dp) function hg19_eqn32_func(x, df_dx, lrpar, rpar, lipar, ipar, ierr) result(f)

      use const_def, only: dp

      real(dp), intent(in)             :: x
      real(dp), intent(out)            :: df_dx
      integer, intent(in)              :: lrpar
      real(dp), intent(inout), pointer :: rpar(:)
      integer, intent(in)              :: lipar
      integer, intent(inout), pointer  :: ipar(:)
      integer, intent(out)             :: ierr

      real(dp) :: LHS
      real(dp) :: RHS1
      real(dp) :: RHS2

      ! Evaluate eqn. 32 of HG19

      if (lrpar /= 5 .OR. lipar /= 0) then
         ierr = -1
         return
      end if

      associate( &
           w => x, &
           l2hat => rpar(1), &
           lhat => rpar(2), &
           lamhat => rpar(3), &
           HB => rpar(4), &
           CH => rpar(5))

        LHS = 0.5_dp * w**2.0_dp - HB
        RHS1 = (CH*lamhat/(0.42_dp*lhat))**1.5_dp
        RHS2 = sqrt(w)

        f = LHS - (RHS1*RHS2)
        df_dx = w - 0.5_dp * (CH*lamhat/(0.42_dp*lhat))**1.5_dp / sqrt(w)         

      end associate

      ierr = 0

      return

   end function hg19_eqn32_func

end module thermohaline
