! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module skye_coulomb
   use math_lib
   use math_def
   use auto_diff
   use const_def, only: dp, pi, rbohr, qe, amu, me, boltzm, five_thirds, kerg
   use skye_composition_ad
   use skye_coulomb_solid
   use skye_coulomb_liquid

   implicit none

   logical, parameter :: dbg = .false.

   private

   public :: nonideal_corrections
   public :: nonideal_corrections_dxa
   public :: nonideal_corrections_dya

   contains

   !! Computes the non-ideal free energy correction for a Coulomb system.
   !! This is done for both the liquid phase and the solid phase, and the resulting
   !! free energies are then combined in a way that blends from the solid phase when
   !! F_solid < F_liquid to the liquid phase when F_liquid < F_solid, with a blend width
   !! of order kT/abar. The blend is done in a way that propagates derivatives, so that the
   !! result is thermodynamically consistent.
   !!
   !! @param NMIX The number of different elements.
   !! @param XA The mass fractions of the different species.
   !! @param AZion The charges of the different species.
   !! @param ACMI The weight of the different species in amu.
   !! @param min_gamma_for_solid The minimum Gamma_i at which to use the solid free energy fit (below this, extrapolate).
   !! @param max_gamma_for_liquid The maximum Gamma_i at which to use the liquid free energy fit (above this, extrapolate).
   !! @param RHO The density of the system in g/cm^3.
   !! @param temp The temperature of the system in K.
   !! @param xnefer The electron density in 1/cm^3.
   !! @param abar The mean atomic mass number.
   !! @param dF The non-ideal correction to the free energy in erg/g.
   !! @param phase The blended phase. 0 for liquid, 1 for solid, smoothly interpolates in between.
   !! @param latent_ddlnT The latent heat of the smoothed phase transition in lnT (T dS/dlnT)
   !! @param latent_ddlnRho The latent heat of the smoothed phase transition in lnRho (T dS/dlnRho)
   subroutine nonideal_corrections(NMIX,AY,AZion,ACMI, min_gamma_for_solid, max_gamma_for_liquid,&
                                   Skye_solid_mixing_rule, RHO,temp, xnefer, abar, &
                                   dF, latent_ddlnT, latent_ddlnRho,phase, phase_ocp)
      integer, intent(in) :: NMIX
      real(dp), intent(in) :: AZion(:), ACMI(:), abar, AY(:), min_gamma_for_solid, max_gamma_for_liquid
      type(auto_diff_real_2var_order3), intent(in) :: RHO, temp, xnefer
      type(auto_diff_real_2var_order3), intent(out) :: dF, phase, latent_ddlnT, latent_ddlnRho
      type(auto_diff_real_2var_order3), intent(out), optional :: phase_ocp(:)
      character(len=128), intent(in) :: Skye_solid_mixing_rule

      integer :: IX
      integer :: LIQSOL
      real(dp) :: Zmean, Z2mean, Z52, Z53, Z321
      type(auto_diff_real_2var_order3) :: GAME, RS, DENS, Smix, F_phase_independent
      type(auto_diff_real_2var_order3) :: kT, dF_sol, dF_liq
      type(auto_diff_real_2var_order3) :: dF_sol_ocp(NMIX), dF_liq_ocp(NMIX)

      ! Compute various mean charge quantities
      Zmean=0d0
      Z2mean=0d0
      Z52=0d0
      Z53=0d0
      Z321=0d0
      do IX=1,NMIX
         Zmean=Zmean+AY(IX)*AZion(IX)
         Z2mean=Z2mean+AY(IX)*AZion(IX)*AZion(IX)
         Z53=Z53+AY(IX)*pow(AZION(IX), five_thirds)
         Z52=Z52+AY(IX)*pow5(sqrt(AZion(IX)))
         Z321=Z321+AY(IX)*AZion(IX)*pow3(sqrt(AZion(IX)+1.d0))
      end do

      DENS = xnefer * pow3(rbohr)  ! DENS = (electrons per cubic bohr)
      RS = pow(3d0 / (4d0 * PI * DENS),1d0/3d0)  ! r_s - electron density parameter
      GAME = qe * qe / (rbohr * boltzm * temp * RS)  ! electron Coulomb parameter Gamma_e

      ! Electron exchange-correlation free energy
      F_phase_independent = Zmean * EXCOR7(RS,GAME)

      ! Linear mixing entropy
      Smix = linear_mixing_entropy(NMIX, Azion, AY)

      ! Incorporate mixing entropy into phase-independent free energy.
      ! Usually this would be incorporated by subtracting T*Smix (because F = E - TS),
      ! but here F is per ion per kB T and dSmix is per ion per kB, so
      ! F -> F - Smix.
      F_phase_independent = F_phase_independent - Smix


      ! Compute free energy correction for liquid and solid phase.
      LIQSOL = 0
      if (present(phase_ocp)) then
         dF_liq = nonideal_corrections_phase( &
            NMIX,AY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid, &
            Skye_solid_mixing_rule, temp,abar,GAME,RS,LIQSOL, &
            Zmean, Z2mean, Z52, Z53, Z321, dF_liq_ocp)
      else
         dF_liq = nonideal_corrections_phase( &
            NMIX,AY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid, &
            Skye_solid_mixing_rule, temp,abar,GAME,RS,LIQSOL, &
            Zmean, Z2mean, Z52, Z53, Z321)
      end if

      LIQSOL = 1
      if (present(phase_ocp)) then
         dF_sol = nonideal_corrections_phase( &
            NMIX,AY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid, &
            Skye_solid_mixing_rule, temp,abar,GAME,RS,LIQSOL, &
            Zmean, Z2mean, Z52, Z53, Z321, dF_sol_ocp)
      else
         dF_sol = nonideal_corrections_phase( &
            NMIX,AY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid, &
            Skye_solid_mixing_rule, temp,abar,GAME,RS,LIQSOL, &
            Zmean, Z2mean, Z52, Z53, Z321)
      end if

      ! Add electron exchange-correlation energy
      dF_liq = dF_liq + F_phase_independent
      dF_sol = dF_sol + F_phase_independent

      ! Change the units from (free energy per kT per ion) to (erg/g).
      kT = temp * kerg / (abar * amu)
      dF_liq = dF_liq * kT
      dF_sol = dF_sol * kT

      ! Produce a smoothed version of the phase transition to extract the latent heat
      call decide_phase(dF_liq, dF_sol, kT, temp, rho, dF, phase, latent_ddlnT, latent_ddlnRho)

      if (present(phase_ocp)) then
         if (phase% val == 0d0) then
            do IX=1,NMIX
               phase_ocp(IX) = dF_liq_ocp(IX)
            end do
         else if (phase% val == 1d0) then
            do IX=1,NMIX
               phase_ocp(IX) = dF_sol_ocp(IX)
            end do
         end if
      end if

      if (dbg) then
         write(*,*) 'GAME',GAME%val,'Phase', phase%val
      end if

   end subroutine nonideal_corrections


   subroutine nonideal_corrections_dxa( &
         NMIX,AY,dAY,AZion,ACMI, min_gamma_for_solid, max_gamma_for_liquid,&
         Skye_solid_mixing_rule, RHO,temp, xnefer, dxnefer_dxa, abar, dabar_dxa, &
         dF_dxa, latent_ddlnT_dxa, latent_ddlnRho_dxa, phase_dxa, composition_only_in_AY, phase_hint)
      integer, intent(in) :: NMIX
      real(dp), intent(in) :: AZion(:), ACMI(:), abar, AY(:), dAY(:)
      real(dp), intent(in) :: min_gamma_for_solid, max_gamma_for_liquid
      real(dp), intent(in) :: dabar_dxa
      type(auto_diff_real_2var_order3), intent(in) :: RHO, temp, xnefer, dxnefer_dxa
      type(auto_diff_real_2var_order3), intent(out) :: dF_dxa
      type(auto_diff_real_2var_order3), intent(out), optional :: phase_dxa
      type(auto_diff_real_2var_order3), intent(out), optional :: latent_ddlnT_dxa, latent_ddlnRho_dxa
      logical, intent(in), optional :: composition_only_in_AY
      integer, intent(in), optional :: phase_hint
      character(len=128), intent(in) :: Skye_solid_mixing_rule

      integer :: IX, LIQSOL, phase_work
      logical :: only_AY
      type(auto_diff_real_2var_order3) :: zero_ad, abar_ad, dabar_ad, AY_ad, dAY_ad
      type(skye_composition_ad_real) :: rho_c, temp_c, xnefer_c, abar_c, y
      type(skye_composition_ad_real) :: Zmean, Z2mean, Z52, Z53, Z321
      type(skye_composition_ad_real) :: GAME, RS, DENS, Smix, F_phase_independent
      type(skye_composition_ad_real) :: kT, dF_sol, dF_liq, dF_c, phase_c
      type(skye_composition_ad_real) :: latent_ddlnT_c, latent_ddlnRho_c

      only_AY = .false.
      if (present(composition_only_in_AY)) only_AY = composition_only_in_AY
      phase_work = -1
      if (present(phase_hint)) then
         if (phase_hint == 0 .or. phase_hint == 1) phase_work = phase_hint
      end if

      zero_ad = 0d0
      rho_c = make_skye_composition_ad(RHO, zero_ad)
      temp_c = make_skye_composition_ad(temp, zero_ad)
      xnefer_c = make_skye_composition_ad(xnefer, dxnefer_dxa)
      abar_ad = abar
      dabar_ad = dabar_dxa
      abar_c = make_skye_composition_ad(abar_ad, dabar_ad)

      Zmean=0d0
      Z2mean=0d0
      Z52=0d0
      Z53=0d0
      Z321=0d0
      do IX=1,NMIX
         AY_ad = AY(IX)
         dAY_ad = dAY(IX)
         y = make_skye_composition_ad(AY_ad, dAY_ad)
         Zmean=Zmean+y*AZion(IX)
         Z2mean=Z2mean+y*AZion(IX)*AZion(IX)
         Z53=Z53+y*pow(AZION(IX), five_thirds)
         Z52=Z52+y*pow5(sqrt(AZion(IX)))
         Z321=Z321+y*AZion(IX)*pow3(sqrt(AZion(IX)+1.d0))
      end do

      DENS = xnefer_c * pow3(rbohr)
      RS = pow(3d0 / (4d0 * PI * DENS),1d0/3d0)
      GAME = qe * qe / (rbohr * boltzm * temp_c * RS)

      F_phase_independent = Zmean * EXCOR7_dxa(RS,GAME)

      Smix = linear_mixing_entropy_dxa(NMIX, Azion, AY, dAY)

      F_phase_independent = F_phase_independent - Smix

      kT = temp_c * kerg / (abar_c * amu)

      if (phase_work == 0 .or. phase_work == 1) then
         LIQSOL = phase_work
         dF_c = nonideal_corrections_phase_dxa(NMIX,AY,dAY,AZion,ACMI,&
             min_gamma_for_solid, max_gamma_for_liquid, Skye_solid_mixing_rule, &
             temp_c,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321, only_AY)
         dF_c = (dF_c + F_phase_independent) * kT
         if (present(phase_dxa)) phase_c = phase_work
         if (present(latent_ddlnT_dxa)) latent_ddlnT_c = 0d0
         if (present(latent_ddlnRho_dxa)) latent_ddlnRho_c = 0d0
      else
         LIQSOL = 0
         dF_liq = nonideal_corrections_phase_dxa(NMIX,AY,dAY,AZion,ACMI,&
             min_gamma_for_solid, max_gamma_for_liquid, Skye_solid_mixing_rule, &
             temp_c,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321, only_AY)

         LIQSOL = 1
         dF_sol = nonideal_corrections_phase_dxa(NMIX,AY,dAY,AZion,ACMI,&
             min_gamma_for_solid, max_gamma_for_liquid, Skye_solid_mixing_rule, &
             temp_c,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321, only_AY)

         dF_liq = (dF_liq + F_phase_independent) * kT
         dF_sol = (dF_sol + F_phase_independent) * kT

         if (present(latent_ddlnT_dxa) .or. present(latent_ddlnRho_dxa)) then
            call decide_phase_dxa( &
               dF_liq, dF_sol, kT, temp_c, rho_c, dF_c, &
               phase=phase_c, latent_ddlnT=latent_ddlnT_c, &
               latent_ddlnRho=latent_ddlnRho_c)
         else if (present(phase_dxa)) then
            call decide_phase_dxa(dF_liq, dF_sol, kT, temp_c, rho_c, dF_c, phase=phase_c)
         else
            call decide_phase_dxa(dF_liq, dF_sol, kT, temp_c, rho_c, dF_c)
         end if
      end if

      dF_dxa = dF_c% d
      if (present(phase_dxa)) phase_dxa = phase_c% d
      if (present(latent_ddlnT_dxa)) latent_ddlnT_dxa = latent_ddlnT_c% d
      if (present(latent_ddlnRho_dxa)) latent_ddlnRho_dxa = latent_ddlnRho_c% d

   end subroutine nonideal_corrections_dxa


   subroutine nonideal_corrections_dya( &
         NMIX,AY,AZion,ACMI, min_gamma_for_solid, max_gamma_for_liquid,&
         Skye_solid_mixing_rule, RHO,temp, xnefer, abar, dF_dya, phase_hint, phase_ocp)
      integer, intent(in) :: NMIX
      real(dp), intent(in) :: AZion(:), ACMI(:), abar, AY(:)
      real(dp), intent(in) :: min_gamma_for_solid, max_gamma_for_liquid
      type(auto_diff_real_2var_order3), intent(in) :: RHO, temp, xnefer
      type(auto_diff_real_2var_order3), intent(out) :: dF_dya(:)
      integer, intent(in), optional :: phase_hint
      type(auto_diff_real_2var_order3), intent(in), optional :: phase_ocp(:)
      character(len=128), intent(in) :: Skye_solid_mixing_rule

      integer :: IX, LIQSOL, phase_work
      real(dp) :: Zmean, Z2mean, Z52, Z53, Z321
      type(auto_diff_real_2var_order3) :: GAME, RS, DENS, Smix, F_phase_independent
      type(auto_diff_real_2var_order3) :: kT, dF_sol, dF_liq, excor
      type(auto_diff_real_2var_order3) :: dF_liq_dya(NMIX), dF_sol_dya(NMIX)

      phase_work = -1
      if (present(phase_hint)) then
         if (phase_hint == 0 .or. phase_hint == 1) phase_work = phase_hint
      end if

      Zmean=0d0
      Z2mean=0d0
      Z52=0d0
      Z53=0d0
      Z321=0d0
      do IX=1,NMIX
         Zmean=Zmean+AY(IX)*AZion(IX)
         Z2mean=Z2mean+AY(IX)*AZion(IX)*AZion(IX)
         Z53=Z53+AY(IX)*pow(AZION(IX), five_thirds)
         Z52=Z52+AY(IX)*pow5(sqrt(AZion(IX)))
         Z321=Z321+AY(IX)*AZion(IX)*pow3(sqrt(AZion(IX)+1.d0))
      end do

      DENS = xnefer * pow3(rbohr)
      RS = pow(3d0 / (4d0 * PI * DENS),1d0/3d0)
      GAME = qe * qe / (rbohr * boltzm * temp * RS)

      excor = EXCOR7(RS,GAME)

      kT = temp * kerg / (abar * amu)

      if (phase_work == 0 .or. phase_work == 1) then
         LIQSOL = phase_work
         if (present(phase_ocp)) then
            call nonideal_corrections_phase_dya(NMIX,AY,AZion,ACMI,&
                min_gamma_for_solid, max_gamma_for_liquid, Skye_solid_mixing_rule, &
                temp,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321, &
                dF_liq, dF_liq_dya, phase_ocp)
         else
            call nonideal_corrections_phase_dya(NMIX,AY,AZion,ACMI,&
                min_gamma_for_solid, max_gamma_for_liquid, Skye_solid_mixing_rule, &
                temp,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321, &
                dF_liq, dF_liq_dya)
         end if
         do IX=1,NMIX
            dF_dya(IX) = (dF_liq_dya(IX) + AZion(IX)*excor + log(AY(IX)) + 1d0) * kT
         end do
      else
         F_phase_independent = Zmean * excor

         Smix = linear_mixing_entropy(NMIX, Azion, AY)

         F_phase_independent = F_phase_independent - Smix

         LIQSOL = 0
         call nonideal_corrections_phase_dya(NMIX,AY,AZion,ACMI,&
             min_gamma_for_solid, max_gamma_for_liquid, Skye_solid_mixing_rule, &
             temp,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321, dF_liq, dF_liq_dya)

         LIQSOL = 1
         call nonideal_corrections_phase_dya(NMIX,AY,AZion,ACMI,&
             min_gamma_for_solid, max_gamma_for_liquid, Skye_solid_mixing_rule, &
             temp,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321, dF_sol, dF_sol_dya)

         do IX=1,NMIX
            dF_liq_dya(IX) = dF_liq_dya(IX) + AZion(IX)*excor + log(AY(IX)) + 1d0
            dF_sol_dya(IX) = dF_sol_dya(IX) + AZion(IX)*excor + log(AY(IX)) + 1d0
         end do

         dF_liq = (dF_liq + F_phase_independent) * kT
         dF_sol = (dF_sol + F_phase_independent) * kT
         do IX=1,NMIX
            dF_liq_dya(IX) = dF_liq_dya(IX) * kT
            dF_sol_dya(IX) = dF_sol_dya(IX) * kT
         end do

         if (dF_liq < dF_sol) then
            do IX=1,NMIX
               dF_dya(IX) = dF_liq_dya(IX)
            end do
         else
            do IX=1,NMIX
               dF_dya(IX) = dF_sol_dya(IX)
            end do
         end if
      end if

   end subroutine nonideal_corrections_dya


   !> Computes the free energy, phase, and latent heat across the phase transition
   !! between liquid and solid. The latent heat is blurred / smeared out over a finite
   !! width in lnT and lnRho so that the solver doesn't need to deal with a Dirac Delta.
   !! @param dF_liq The free energy of the liquid phase in erg/g.
   !! @param dF_sol The free energy of the solid phase in erg/g.
   !! @param kT The factor kB*T/(abar * amu), the thermal energy per mean ion mass. This is in erg/g.
   !! @param temp The temperature in K.
   !! @param rho The density in g/cm^3.
   !! @param phase The phase (0 for liquid, 1 for solid, in-between to represent the blur).
   !! @param latent_ddlnT Equals T^2 dS/dT for the latent heat. Equivalently, d(latent heat)/dlnT.
   !! @param latent_ddlnRho Equals T Rho dS/dRho for the latent heat. Equivalently, d(latent heat)/dlnRho.
   subroutine decide_phase(dF_liq, dF_sol, kT, temp, rho, dF, phase, latent_ddlnT, latent_ddlnRho)
      type(auto_diff_real_2var_order3), intent(in) :: dF_liq, dF_sol, kT, temp, rho
      type(auto_diff_real_2var_order3), intent(out) :: dF, phase, latent_ddlnT, latent_ddlnRho

      type(auto_diff_real_2var_order3) :: blur, dF_blur, latent_S
      real(dp), parameter :: blur_width = 1d2

      ! Pick the phase with the minimum free energy
      dF = min(dF_liq, dF_sol)

      ! Next, blur the free energy so it's twice differentiable across the phase
      ! transition, and use that to extract the latent heat.
      blur = kT / blur_width

      ! Avoid overflow or underflow
      if (dF_liq < dF_sol - 20d0*blur) then
         phase = 0d0
      else if (dF_sol < dF_liq - 20d0*blur) then
         phase = 1d0
      else
         ! Mix phases with a blur (softmin)
         phase = exp((dF_liq-dF_sol)/blur) / (exp((dF_liq-dF_sol)/blur) + 1d0)
      end if
      dF_blur = dF_liq * (1d0 - phase) + dF_sol * phase

      ! Now subtract off the regular (un-blurred) free energy so we don't
      ! double-count the off-transition dS/dT and dS/dRho in the latent heat.
      dF_blur = dF_blur - dF

      ! Latent entropy
      latent_S = -(differentiate_1(dF_blur))  ! S = -dF/dT

      ! T dS/dlnT = T^2 dS/dT
      latent_ddlnT = differentiate_1(latent_S) *  pow2(temp)

      ! T dS/dlnRho = T Rho dS/dRho
      latent_ddlnRho = temp * rho * differentiate_2(latent_S)

   end subroutine decide_phase


   subroutine decide_phase_dxa( &
         dF_liq, dF_sol, kT, temp, rho, dF, phase, latent_ddlnT, latent_ddlnRho)
      type(skye_composition_ad_real), intent(in) :: dF_liq, dF_sol, kT, temp, rho
      type(skye_composition_ad_real), intent(out) :: dF
      type(skye_composition_ad_real), intent(out), optional :: phase, latent_ddlnT, latent_ddlnRho

      type(skye_composition_ad_real) :: blur, dF_blur, latent_S, phase_work
      logical :: need_latent
      real(dp), parameter :: blur_width = 1d2

      dF = min(dF_liq, dF_sol)
      need_latent = present(latent_ddlnT) .or. present(latent_ddlnRho)
      if (.not. present(phase) .and. .not. need_latent) return

      blur = kT / blur_width

      if (dF_liq < dF_sol - 20d0*blur) then
         phase_work = 0d0
      else if (dF_sol < dF_liq - 20d0*blur) then
         phase_work = 1d0
      else
         phase_work = exp((dF_liq-dF_sol)/blur) / (exp((dF_liq-dF_sol)/blur) + 1d0)
      end if
      if (present(phase)) phase = phase_work
      if (.not. need_latent) return

      dF_blur = dF_liq * (1d0 - phase_work) + dF_sol * phase_work

      dF_blur = dF_blur - dF

      latent_S = -(differentiate_1(dF_blur))

      if (present(latent_ddlnT)) latent_ddlnT = differentiate_1(latent_S) *  pow2(temp)

      if (present(latent_ddlnRho)) latent_ddlnRho = temp * rho * differentiate_2(latent_S)

   end subroutine decide_phase_dxa


   function substitute_composition_binop(F, x, y) result(z)
      type(skye_composition_ad_real), intent(in) :: F, x, y

      type(auto_diff_real_2var_order3) :: zval, zd
      type(skye_composition_ad_real) :: z

      zval = make_binop(x% val, y% val, &
         F% val% val, F% val% d1val1, F% val% d1val2, &
         F% val% d2val1, F% val% d1val1_d1val2, F% val% d2val2, &
         F% val% d3val1, F% val% d2val1_d1val2, &
         F% val% d1val1_d2val2, F% val% d3val2)
      zd = make_binop(x% val, y% val, &
         F% d% val, F% d% d1val1, F% d% d1val2, &
         F% d% d2val1, F% d% d1val1_d1val2, F% d% d2val2, &
         F% d% d3val1, F% d% d2val1_d1val2, &
         F% d% d1val1_d2val2, F% d% d3val2)

      z = make_skye_composition_ad(zval, zd)
   end function substitute_composition_binop

   !> Computes the linear mixing entropy per ion per kB.
   !! @param Nmix The number of species in the mixture.
   !! @param Azion An array of the charges of those species.
   !! @param AY An array of the number fractions of those species.
   !! @param Smix The linear mixing entropy per ion per kB.
   type(auto_diff_real_2var_order3) function linear_mixing_entropy(Nmix, AZion, AY) result(Smix)
      ! Inputs
      integer, intent(in) :: Nmix
      real(dp), intent(in) :: AZion(:), AY(:)

      ! Intermediates and constants
      integer :: i

      Smix = 0d0
      do i=1,Nmix
         Smix = Smix - AY(i)*log(AY(i))
      end do
   end function linear_mixing_entropy


   function linear_mixing_entropy_dxa(Nmix, AZion, AY, dAY) result(Smix)
      ! Inputs
      integer, intent(in) :: Nmix
      real(dp), intent(in) :: AZion(:), AY(:), dAY(:)

      ! Intermediates and constants
      integer :: i
      type(auto_diff_real_2var_order3) :: AY_ad, dAY_ad
      type(skye_composition_ad_real) :: y

      ! Output
      type(skye_composition_ad_real) :: Smix

      Smix = 0d0
      do i=1,Nmix
         AY_ad = AY(i)
         dAY_ad = dAY(i)
         y = make_skye_composition_ad(AY_ad, dAY_ad)
         Smix = Smix - y*log(y)
      end do
   end function linear_mixing_entropy_dxa

   !> Computes the non-ideal one-component plasma corrections to the free energy.
   !! This involves a loop over species to compute the free energy of a
   !! one-component plasma, which is added to the non-ideal mixing free energy.
   !!
   !! @param Nmix The number of species in the mixture.
   !! @param AY An array of the number fractions of those species.
   !! @param Azion An array of the charges in electron charges of those species.
   !! @param ACMI An array of the masses in AMU of those species.
   !! @param min_gamma_for_solid The minimum Gamma_i at which to use the solid free energy fit (below this, extrapolate).
   !! @param max_gamma_for_liquid The maximum Gamma_i at which to use the liquid free energy fit (above this, extrapolate).
   !! @param temp The temperature in K.
   !! @param abar The mean atomic mass number.
   !! @param RS Electron density parameter for component species
   !! @param GAME Electron coupling parameter (Gamma_i)
   !! @param LIQSOL Integer specifying the phase: 0 for liquid, 1 for solid
   !! @param Zmean The mean ion charge (mass fraction weighted)
   !! @param Z2mean The mean squared ion charge (mass fraction weighted)
   !! @param Z53mean The mean of ion charge to the 5/3 power (mass fraction weighted)
   !! @param Z321mean The mean of Z(Z+1)^(3/2), where Z is the ion charge (mass fraction weighted)
   function nonideal_corrections_phase( &
         NMIX,AY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid, &
         Skye_solid_mixing_rule, temp,abar,GAME,RS,LIQSOL, &
         Zmean, Z2mean, Z52, Z53, Z321, phase_ocp) result(dF)
      ! Inputs
      integer, intent(in) :: NMIX
      integer, intent(in) :: LIQSOL
      real(dp), intent(in) :: AZion(:), ACMI(:), abar, AY(:), Zmean, Z2mean, Z52, Z53, Z321, &
                              min_gamma_for_solid, max_gamma_for_liquid
      type(auto_diff_real_2var_order3), intent(in) :: temp, GAME, RS
      type(auto_diff_real_2var_order3), intent(out), optional :: phase_ocp(:)
      character(len=128), intent(in) :: Skye_solid_mixing_rule

      ! Intermediates and constants
      integer :: i
      type(auto_diff_real_2var_order3) :: FMIX, f
      real(dp), parameter :: TINY=1.d-7

      ! Output
      type(auto_diff_real_2var_order3) :: dF

      dF = 0d0
      if (present(phase_ocp)) then
         do i=1,NMIX
            phase_ocp(i) = 0d0
         end do
      end if

      ! Composition loop
      do i=1,nmix
         if (AY(i) > TINY .and. AZion(i) /= 0d0) then  ! skip low-abundance species and neutrons

            ! Add up non-ideal corrections
            f = extrapolate_free_energy(LIQSOL, temp, RS, AZion(i), ACMI(i), min_gamma_for_solid, max_gamma_for_liquid)
            if (LIQSOL == 0) then
               f = f + ocp_liquid_screening_free_energy_correction(AZion(i), ACMI(i), GAME, RS)  ! screening corrections
            else
               f = f + ocp_solid_screening_free_energy_correction(AZion(i), ACMI(i), GAME, RS)  ! screening corrections
            end if
            dF = dF + AY(i) * f
            if (present(phase_ocp)) phase_ocp(i) = f

         end if
      end do

      ! Corrections to the linear mixing rule:
      if (LIQSOL == 0) then  ! liquid phase
         FMIX = liquid_mixing_rule_correction(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321)
      else  ! solid phase (only Madelung contribution) [22.12.12]
         FMIX = solid_mixing_rule_correction(Skye_solid_mixing_rule, NMIX, AY, AZion, GAME)
      end if
      dF = dF + FMIX

   end function nonideal_corrections_phase


   function nonideal_corrections_phase_dxa( &
         NMIX,AY,dAY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid, &
         Skye_solid_mixing_rule, temp,abar,GAME,RS,LIQSOL, &
         Zmean, Z2mean, Z52, Z53, Z321, composition_only_in_AY) result(dF)
      ! Inputs
      integer, intent(in) :: NMIX
      integer, intent(in) :: LIQSOL
      real(dp), intent(in) :: AZion(:), ACMI(:), abar, AY(:), dAY(:), &
                              min_gamma_for_solid, max_gamma_for_liquid
      type(skye_composition_ad_real), intent(in) :: temp, GAME, RS
      type(skye_composition_ad_real), intent(in) :: Zmean, Z2mean, Z52, Z53, Z321
      logical, intent(in), optional :: composition_only_in_AY
      character(len=128), intent(in) :: Skye_solid_mixing_rule

      ! Intermediates and constants
      integer :: i
      type(auto_diff_real_2var_order3) :: AY_ad, dAY_ad
      type(auto_diff_real_2var_order3) :: f_ad
      type(skye_composition_ad_real) :: FMIX, f, y
      logical :: only_AY
      real(dp), parameter :: TINY=1.d-7

      ! Output
      type(skye_composition_ad_real) :: dF

      only_AY = .false.
      if (present(composition_only_in_AY)) only_AY = composition_only_in_AY

      dF = 0d0

      ! Composition loop
      do i=1,nmix
         if (AY(i) > TINY .and. AZion(i) /= 0d0) then

            ! For d/dYA basis vectors, only AY carries the composition perturbation.
            ! Use the ordinary OCP AD path and attach dAY by the product rule.
            if (only_AY) then
               f_ad = extrapolate_free_energy(LIQSOL, temp%val, RS%val, AZion(i), ACMI(i), &
                  min_gamma_for_solid, max_gamma_for_liquid)
               if (LIQSOL == 0) then
                  f_ad = f_ad + ocp_liquid_screening_free_energy_correction( &
                     AZion(i), ACMI(i), GAME%val, RS%val)
               else
                  f_ad = f_ad + ocp_solid_screening_free_energy_correction( &
                     AZion(i), ACMI(i), GAME%val, RS%val)
               end if
               dF%val = dF%val + AY(i)*f_ad
               dF%d = dF%d + dAY(i)*f_ad
            else
               f = extrapolate_free_energy_dxa(LIQSOL, temp, RS, AZion(i), ACMI(i), &
                  min_gamma_for_solid, max_gamma_for_liquid)
               if (LIQSOL == 0) then
                  f = f + ocp_liquid_screening_free_energy_correction_dxa(AZion(i), ACMI(i), GAME, RS)
               else
                  f = f + ocp_solid_screening_free_energy_correction_dxa(AZion(i), ACMI(i), GAME, RS)
               end if
               AY_ad = AY(i)
               dAY_ad = dAY(i)
               y = make_skye_composition_ad(AY_ad, dAY_ad)
               dF = dF + y * f
            end if

         end if
      end do

      ! Corrections to the linear mixing rule:
      if (LIQSOL == 0) then
         FMIX = liquid_mixing_rule_correction_dxa(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321)
      else
         FMIX = solid_mixing_rule_correction_dxa(Skye_solid_mixing_rule, NMIX, AY, dAY, AZion, GAME)
      end if
      dF = dF + FMIX

   end function nonideal_corrections_phase_dxa


   subroutine nonideal_corrections_phase_dya( &
         NMIX,AY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid, &
         Skye_solid_mixing_rule, temp,abar,GAME,RS,LIQSOL, &
         Zmean, Z2mean, Z52, Z53, Z321, dF, dF_dya, phase_ocp)
      integer, intent(in) :: NMIX
      integer, intent(in) :: LIQSOL
      real(dp), intent(in) :: AZion(:), ACMI(:), abar, AY(:), &
                              min_gamma_for_solid, max_gamma_for_liquid
      type(auto_diff_real_2var_order3), intent(in) :: temp, GAME, RS
      real(dp), intent(in) :: Zmean, Z2mean, Z52, Z53, Z321
      type(auto_diff_real_2var_order3), intent(in), optional :: phase_ocp(:)
      character(len=128), intent(in) :: Skye_solid_mixing_rule
      type(auto_diff_real_2var_order3), intent(out) :: dF, dF_dya(:)

      integer :: i
      real(dp) :: dAY(NMIX)
      type(auto_diff_real_2var_order3) :: FMIX, f, zero_ad
      type(auto_diff_real_2var_order3) :: Zmean_ad, Z2mean_ad, Z52_ad, Z53_ad, Z321_ad
      type(auto_diff_real_2var_order3) :: dZmean, dZ2mean, dZ52, dZ53, dZ321
      type(skye_composition_ad_real) :: FMIX_c
      type(skye_composition_ad_real) :: RS_c, GAME_c, Zmean_c, Z2mean_c, Z52_c, Z53_c, Z321_c
      real(dp), parameter :: TINY=1.d-7

      dF = 0d0
      do i=1,NMIX
         dF_dya(i) = 0d0
      end do

      ! Batch the d/dYA path: each OCP leaf is independent of dYA, so evaluate
      ! it once per active species and reuse it for all number-fraction columns.
      do i=1,nmix
         if (AY(i) > TINY .and. AZion(i) /= 0d0) then
            if (present(phase_ocp)) then
               f = phase_ocp(i)
            else
               f = extrapolate_free_energy(LIQSOL, temp, RS, AZion(i), ACMI(i), &
                  min_gamma_for_solid, max_gamma_for_liquid)
               if (LIQSOL == 0) then
                  f = f + ocp_liquid_screening_free_energy_correction( &
                     AZion(i), ACMI(i), GAME, RS)
               else
                  f = f + ocp_solid_screening_free_energy_correction( &
                     AZion(i), ACMI(i), GAME, RS)
               end if
            end if
            dF = dF + AY(i) * f
            dF_dya(i) = f
         end if
      end do

      zero_ad = 0d0
      RS_c = make_skye_composition_ad(RS, zero_ad)
      GAME_c = make_skye_composition_ad(GAME, zero_ad)
      if (LIQSOL == 0) then
         FMIX = liquid_mixing_rule_correction(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321)
         dF = dF + FMIX
         Zmean_ad = Zmean
         Z2mean_ad = Z2mean
         Z52_ad = Z52
         Z53_ad = Z53
         Z321_ad = Z321
         do i=1,NMIX
            dZmean = AZion(i)
            dZ2mean = AZion(i)*AZion(i)
            dZ53 = pow(AZION(i), five_thirds)
            dZ52 = pow5(sqrt(AZion(i)))
            dZ321 = AZion(i)*pow3(sqrt(AZion(i)+1.d0))
            Zmean_c = make_skye_composition_ad(Zmean_ad, dZmean)
            Z2mean_c = make_skye_composition_ad(Z2mean_ad, dZ2mean)
            Z52_c = make_skye_composition_ad(Z52_ad, dZ52)
            Z53_c = make_skye_composition_ad(Z53_ad, dZ53)
            Z321_c = make_skye_composition_ad(Z321_ad, dZ321)
            FMIX_c = liquid_mixing_rule_correction_dxa( &
               RS_c,GAME_c,Zmean_c,Z2mean_c,Z52_c,Z53_c,Z321_c)
            dF_dya(i) = dF_dya(i) + FMIX_c% d
         end do
      else
         FMIX = solid_mixing_rule_correction(Skye_solid_mixing_rule, NMIX, AY, AZion, GAME)
         dF = dF + FMIX
         dAY = 0d0
         do i=1,NMIX
            if (i > 1) dAY(i-1) = 0d0
            dAY(i) = 1d0
            FMIX_c = solid_mixing_rule_correction_dxa( &
               Skye_solid_mixing_rule, NMIX, AY, dAY, AZion, GAME_c)
            dF_dya(i) = dF_dya(i) + FMIX_c% d
         end do
      end if

   end subroutine nonideal_corrections_phase_dya

   !> Extrapolates the free energy of a one-component plasma beyond the boundaries of the data used
   !! to construct the fitting functions. This is done by fixing the probability distribution over states
   !! to its value on the boundary, which then fixes the internal energy and entropy to their boundary values.
   !! The free energy is then just linear in temperature (F = F(boundary) - S(boundary) * (T - T(boundary))),
   !!
   !! @param LIQSOL Integer specifying the phase: 0 for liquid, 1 for solid
   !! @param temp Temperature (K)
   !! @param RS Electron density parameter
   !! @param Zion Charge of the species of interest in electron charges.
   !! @param CMI Mass of the species of interest in AMU.
   !! @param F non-ideal free energy per ion per kT
   function extrapolate_free_energy(LIQSOL, temp, RS, Zion, CMI, min_gamma_for_solid, max_gamma_for_liquid) result(F)
      ! Inputs
      integer, intent(in) :: LIQSOL
      type(auto_diff_real_2var_order3), intent(in) :: temp, RS
      real(dp), intent(in) :: Zion, CMI, min_gamma_for_solid, max_gamma_for_liquid

      ! Intermediates
      real(dp) :: COTPT, gamma_boundary
      type(auto_diff_real_2var_order3) :: temp_boundary, fake_dens, GAMI, TPT, g, tp, dF_dlnT

      ! Output
      type(auto_diff_real_2var_order3) :: F

      GAMI = pre_z(int(Zion))% z5_3 * qe * qe / (rbohr * boltzm * temp * RS)  ! ion Coulomb parameter Gamma_i
      COTPT=sqrt(3d0 * me_in_amu / CMI)/pre_z(int(Zion))%z7_6  ! auxiliary coefficient
      TPT=GAMI*COTPT/sqrt(RS)                   ! T_p/T

      if ((LIQSOL == 0 .and. GAMI < max_gamma_for_liquid) .or. (LIQSOL == 1 .and. GAMI > min_gamma_for_solid)) then
         ! No extrapolation needed
         F = ocp_free_energy(LIQSOL, Zion, CMI, GAMI, TPT)
         if (dbg) then
            write(*,*) 'Species:', Zion, 'LIQSOL', LIQSOL, 'Normal, GAMI:', GAMI%val, 'F', F%val
         end if
      else
         ! Extrapolate past the boundary
         if (dbg) then
            write(*,*) 'Species:', Zion, 'LIQSOL', LIQSOL, 'Extrapolated, GAMI:', GAMI%val
         end if

         ! Identify the boundary
         if (LIQSOL == 0) then
            gamma_boundary = max_gamma_for_liquid
         else
            gamma_boundary = min_gamma_for_solid
         end if

         ! Find the temperature where Gamma = Gamma_boundary
         temp_boundary = temp%val * GAMI%val / gamma_boundary

         ! Make d(temp_boundary)/dT = 1 so we can extract dF/dT at the boundary.
         temp_boundary%d1val1 = 1d0

         ! Compute new (differentiable) Gamma and TPT at the boundary
         g = pre_z(int(Zion))% z5_3 * qe * qe / (rbohr * boltzm * temp_boundary * RS)  ! ion Coulomb parameter Gamma_i
         tp=g*COTPT/sqrt(RS)                  ! T_p/T

         ! Compute boundary free energy
         F = ocp_free_energy(LIQSOL, Zion, CMI, g, tp)

         ! Extract derivative at boundary
         dF_dlnT = differentiate_1(F) * temp_boundary

         ! Now make the substitution T -> T_boundary(rho).
         ! We do this by treating F and dF_dlnT as binary operators of (T, rho),
         ! and letting T = T_boundary(rho).

         ! First, find the temperature where Gamma = Gamma_boundary
         ! Note that because GAMI \propto 1/temp, this ends up being
         ! independent of temperature but dependent on density.
         temp_boundary = temp * GAMI / gamma_boundary

         ! Next, construct a dummy operator
         ! with d(fake_dens)/d(dens) = 1.
         ! The value doesn't matter because only derivatives
         ! of this will be used in the chain rule for the
         ! substitution.
         fake_dens = 0d0
         fake_dens%d1val2 = 1d0

         ! Then make the substitution. This uses the chain rule to replace all dependences on temperature with dependences
         ! on temp_boundary, which in turn depends on rho. Note that the order of the arguments matters here: above we've
         ! made temperature the first independent variable, so the thing we're substituting for it (temp_boundary) has to be
         ! the first argument here.
         dF_dlnT = make_binop(temp_boundary, fake_dens, &
                     dF_dlnT%val, dF_dlnT%d1val1, dF_dlnT%d1val2, dF_dlnT%d2val1, dF_dlnT%d1val1_d1val2, &
                     dF_dlnT%d2val2, dF_dlnT%d3val1, dF_dlnT%d2val1_d1val2, dF_dlnT%d1val1_d2val2, dF_dlnT%d3val2)

         F = make_binop(temp_boundary, fake_dens, F%val, F%d1val1, F%d1val2, F%d2val1, &
                     F%d1val1_d1val2, F%d2val2, F%d3val1,  &
                     F%d2val1_d1val2, F%d1val1_d2val2, F%d3val2)

         ! Extrapolate
         F = F + (1d0 - temp_boundary / temp) * dF_dlnT
      end if

   end function extrapolate_free_energy


   function extrapolate_free_energy_dxa(LIQSOL, temp, RS, Zion, CMI, &
         min_gamma_for_solid, max_gamma_for_liquid) result(F)
      ! Inputs
      integer, intent(in) :: LIQSOL
      type(skye_composition_ad_real), intent(in) :: temp, RS
      real(dp), intent(in) :: Zion, CMI, min_gamma_for_solid, max_gamma_for_liquid

      ! Intermediates
      real(dp) :: COTPT, gamma_boundary
      type(skye_composition_ad_real) :: temp_boundary, fake_dens
      type(skye_composition_ad_real) :: GAMI, TPT, g, tp, dF_dlnT

      ! Output
      type(skye_composition_ad_real) :: F

      GAMI = pre_z(int(Zion))% z5_3 * qe * qe / (rbohr * boltzm * temp * RS)
      COTPT=sqrt(3d0 * me_in_amu / CMI)/pre_z(int(Zion))%z7_6
      TPT=GAMI*COTPT/sqrt(RS)

      if ((LIQSOL == 0 .and. GAMI < max_gamma_for_liquid) .or. &
            (LIQSOL == 1 .and. GAMI > min_gamma_for_solid)) then
         F = ocp_free_energy_dxa(LIQSOL, Zion, CMI, GAMI, TPT)
      else
         if (LIQSOL == 0) then
            gamma_boundary = max_gamma_for_liquid
         else
            gamma_boundary = min_gamma_for_solid
         end if

         temp_boundary = temp * GAMI / gamma_boundary
         temp_boundary% val = temp_boundary% val% val
         temp_boundary% val% d1val1 = 1d0

         g = pre_z(int(Zion))% z5_3 * qe * qe / (rbohr * boltzm * temp_boundary * RS)
         tp=g*COTPT/sqrt(RS)

         F = ocp_free_energy_dxa(LIQSOL, Zion, CMI, g, tp)

         dF_dlnT = differentiate_1(F) * temp_boundary

         temp_boundary = temp * GAMI / gamma_boundary

         fake_dens = 0d0
         fake_dens% val% d1val2 = 1d0

         dF_dlnT = substitute_composition_binop(dF_dlnT, temp_boundary, fake_dens)
         F = substitute_composition_binop(F, temp_boundary, fake_dens)

         F = F + (1d0 - temp_boundary / temp) * dF_dlnT
      end if

   end function extrapolate_free_energy_dxa

   !> Calculates the free energy of a one-component
   !! plasma in the specified phase. This routine is
   !! just responsible for assembling different terms together.
   !! Based on EOSFI8 by Potekhin and Chabrier.
   !!
   !! @param LIQSOL Integer specifying the phase: 0 for liquid, 1 for solid.
   !! @param Zion Charge of the species of interest in electron charges.
   !! @param CMI Mass of the species of interest in AMU.
   !! @param GAMI Ion coupling parameter (Gamma_i)
   !! @param TPT effective T_p/T - ion quantum parameter
   !! @param F non-ideal free energy per ion per kT
   function ocp_free_energy(LIQSOL, Zion, CMI, GAMI, TPT) result(F)
      ! Inputs
      real(dp), intent(in) :: Zion, CMI
      integer, intent(in) :: LIQSOL
      type(auto_diff_real_2var_order3), intent(in) :: GAMI, TPT

      ! Output
      type(auto_diff_real_2var_order3) :: F

      if (LIQSOL == 0) then
         F = classical_ocp_liquid_free_energy(GAMI)                  ! classical ion-ion interaction
         F = F + quantum_ocp_liquid_free_energy_correction(TPT)   ! quantum ion-ion corrections
      else
         F = ocp_solid_harmonic_free_energy(GAMI,TPT)  ! harmonic classical and quantum ion-ion corrections
         F = F + ocp_solid_anharmonic_free_energy(GAMI,TPT)  ! anharmonic classical and quantum ion-ion corrections
      end if

   end function ocp_free_energy


   function ocp_free_energy_dxa(LIQSOL, Zion, CMI, GAMI, TPT) result(F)
      ! Inputs
      real(dp), intent(in) :: Zion, CMI
      integer, intent(in) :: LIQSOL
      type(skye_composition_ad_real), intent(in) :: GAMI, TPT

      ! Output
      type(skye_composition_ad_real) :: F

      if (LIQSOL == 0) then
         F = classical_ocp_liquid_free_energy_dxa(GAMI)
         F = F + quantum_ocp_liquid_free_energy_correction_dxa(TPT)
      else
         F = ocp_solid_harmonic_free_energy_dxa(GAMI,TPT)
         F = F + ocp_solid_anharmonic_free_energy_dxa(GAMI,TPT)
      end if

   end function ocp_free_energy_dxa


   !> Calculates the electron exchange-correlation non-ideal free energy correction.
   !! Based on the results of Tanaka & Ichimaru 85-87.
   !! and the routine EXCOR7 Version 09.06.07 by Potekhin and Chabrier
   !!
   !! @param RS Electron density parameter for component species
   !! @param GAME Electron coupling parameter (Gamma_i)
   !! @param FXC Electron exchange-correlation non-ideal free energy correction per electron per kT
   function EXCOR7(RS,GAME) result(FXC)
      ! Inputs
      type(auto_diff_real_2var_order3), intent(in) :: RS, GAME

      ! Intermediates
      type(auto_diff_real_2var_order3) :: THETA, SQTH, THETA2, THETA3, THETA4, EXP1TH
      type(auto_diff_real_2var_order3) :: CHT1, SHT1, CHT2, SHT2
      type(auto_diff_real_2var_order3) :: T1, T2
      type(auto_diff_real_2var_order3) :: A0, A1, A, B0, B1, B
      type(auto_diff_real_2var_order3) :: C, C3
      type(auto_diff_real_2var_order3) :: D0, D1, D, E0, E1, E
      type(auto_diff_real_2var_order3) :: DISCR, SQGE
      type(auto_diff_real_2var_order3) :: B2, B3, R3, S1, S2, S3, B4, C4
      type(auto_diff_real_2var_order3) :: S4A, S4B, S4C, S4

      ! Output
      type(auto_diff_real_2var_order3) :: FXC

      THETA=0.543d0*RS/GAME  ! non-relativistic degeneracy parameter
      SQTH=sqrt(THETA)
      THETA2=THETA*THETA
      THETA3=THETA2*THETA
      THETA4=THETA3*THETA
      if (THETA>.007d0) then
         CHT1=cosh(1.d0/THETA)
         SHT1=sinh(1.d0/THETA)
         CHT2=cosh(1.d0/SQTH)
         SHT2=sinh(1.d0/SQTH)
         T1=SHT1/CHT1  ! dtanh(1.d0/THETA)
         T2=SHT2/CHT2  ! dtanh(1./sqrt(THETA))
      else
         T1=1.d0
         T2=1.d0
      end if

      A0=0.75d0+3.04363d0*THETA2-0.09227d0*THETA3+1.7035d0*THETA4
      A1=1d0+8.31051d0*THETA2+5.1105d0*THETA4
      A=0.610887d0*A0/A1*T1  ! HF fit of Perrot and Dharma-wardana

      B0=0.341308d0+12.070873d0*THETA2+1.148889d0*THETA4
      B1=1d0+10.495346d0*THETA2+1.326623d0*THETA4
      B=SQTH*T2*B0/B1

      D0=0.614925d0+16.996055d0*THETA2+1.489056d0*THETA4
      D1=1d0+10.10935d0*THETA2+1.22184d0*THETA4
      D=SQTH*T2*D0/D1

      E0=0.539409d0+2.522206d0*THETA2+0.178484d0*THETA4
      E1=1d0+2.555501d0*THETA2+0.146319d0*THETA4
      E=THETA*T1*E0/E1

      EXP1TH=exp(-1.d0/THETA)
      C=(0.872496d0+0.025248d0*EXP1TH)*E

      DISCR=SQRT(4.0d0*E-D*D)

      S1=-C/E*GAME

      B2=B-C*D/E

      SQGE=SQRT(GAME)
      S2=-2.d0/E*B2*SQGE

      R3=E*GAME+D*SQGE+1.0d0

      B3=A-C/E
      C3=(D/E*B2-B3)/E
      S3=C3*log(R3)

      B4=2.d0-D*D/E

      C4=2d0*E*SQGE+D

      S4A=2.0d0/E/DISCR
      S4B=D*B3+B4*B2
      S4C=atan(C4/DISCR)-atan(D/DISCR)
      S4=S4A*S4B*S4C

      FXC=S1+S2+S3+S4

   end function EXCOR7


   function EXCOR7_dxa(RS,GAME) result(FXC)
      ! Inputs
      type(skye_composition_ad_real), intent(in) :: RS, GAME

      ! Intermediates
      type(skye_composition_ad_real) :: THETA, SQTH, THETA2, THETA3, THETA4, EXP1TH
      type(skye_composition_ad_real) :: CHT1, SHT1, CHT2, SHT2
      type(skye_composition_ad_real) :: T1, T2
      type(skye_composition_ad_real) :: A0, A1, A, B0, B1, B
      type(skye_composition_ad_real) :: C, C3
      type(skye_composition_ad_real) :: D0, D1, D, E0, E1, E
      type(skye_composition_ad_real) :: DISCR, SQGE
      type(skye_composition_ad_real) :: B2, B3, R3, S1, S2, S3, B4, C4
      type(skye_composition_ad_real) :: S4A, S4B, S4C, S4

      ! Output
      type(skye_composition_ad_real) :: FXC

      THETA=0.543d0*RS/GAME
      SQTH=sqrt(THETA)
      THETA2=THETA*THETA
      THETA3=THETA2*THETA
      THETA4=THETA3*THETA
      if (THETA>.007d0) then
         CHT1=cosh(1.d0/THETA)
         SHT1=sinh(1.d0/THETA)
         CHT2=cosh(1.d0/SQTH)
         SHT2=sinh(1.d0/SQTH)
         T1=SHT1/CHT1
         T2=SHT2/CHT2
      else
         T1=1.d0
         T2=1.d0
      end if

      A0=0.75d0+3.04363d0*THETA2-0.09227d0*THETA3+1.7035d0*THETA4
      A1=1d0+8.31051d0*THETA2+5.1105d0*THETA4
      A=0.610887d0*A0/A1*T1

      B0=0.341308d0+12.070873d0*THETA2+1.148889d0*THETA4
      B1=1d0+10.495346d0*THETA2+1.326623d0*THETA4
      B=SQTH*T2*B0/B1

      D0=0.614925d0+16.996055d0*THETA2+1.489056d0*THETA4
      D1=1d0+10.10935d0*THETA2+1.22184d0*THETA4
      D=SQTH*T2*D0/D1

      E0=0.539409d0+2.522206d0*THETA2+0.178484d0*THETA4
      E1=1d0+2.555501d0*THETA2+0.146319d0*THETA4
      E=THETA*T1*E0/E1

      EXP1TH=exp(-1.d0/THETA)
      C=(0.872496d0+0.025248d0*EXP1TH)*E

      DISCR=SQRT(4.0d0*E-D*D)

      S1=-C/E*GAME

      B2=B-C*D/E

      SQGE=SQRT(GAME)
      S2=-2.d0/E*B2*SQGE

      R3=E*GAME+D*SQGE+1.0d0

      B3=A-C/E
      C3=(D/E*B2-B3)/E
      S3=C3*log(R3)

      B4=2.d0-D*D/E

      C4=2d0*E*SQGE+D

      S4A=2.0d0/E/DISCR
      S4B=D*B3+B4*B2
      S4C=atan(C4/DISCR)-atan(D/DISCR)
      S4=S4A*S4B*S4C

      FXC=S1+S2+S3+S4

   end function EXCOR7_dxa

   end module skye_coulomb
