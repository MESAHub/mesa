module skye_coulomb
   use math_lib
   use auto_diff
   use const_def, only: dp, PI, rbohr, qe, amu, me
   use skye_coulomb_solid
   use skye_coulomb_liquid

   implicit none

   logical, parameter :: dbg = .false.
   !logical, parameter :: dbg = .true.

   private

   public :: nonideal_corrections
   
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
                                   dF, latent_ddlnT, latent_ddlnRho,phase)
      integer, intent(in) :: NMIX
      real(dp), intent(in) :: AZion(:), ACMI(:), abar, AY(:), min_gamma_for_solid, max_gamma_for_liquid
      type(auto_diff_real_2var_order3), intent(in) :: RHO, temp, xnefer
      type(auto_diff_real_2var_order3), intent(out) :: dF, phase, latent_ddlnT, latent_ddlnRho
      character(len=128), intent(in) :: Skye_solid_mixing_rule

      integer :: IX
      integer :: LIQSOL
      real(dp) :: Zion, Zmean, Z2mean, Z52, Z53, Z321, norm
      type(auto_diff_real_2var_order3) :: GAME, RS, DENS, Smix, F_phase_independent
      type(auto_diff_real_2var_order3) :: kT, dF_sol, dF_liq, latent_S, min_S

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
      enddo

      DENS = xnefer * pow3(rbohr) ! DENS = (electrons per cubic bohr)
      RS = pow(3d0 / (4d0 * PI * DENS),1d0/3d0) ! r_s - electron density parameter
      GAME = qe * qe / (rbohr * boltzm * temp * RS) ! electron Coulomb parameter Gamma_e

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
      dF_liq = nonideal_corrections_phase(NMIX,AY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid,&
          Skye_solid_mixing_rule, temp,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321)

      LIQSOL = 1
      dF_sol = nonideal_corrections_phase(NMIX,AY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid,&
          Skye_solid_mixing_rule, temp,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321)

      ! Add electron exchange-correlation energy
      dF_liq = dF_liq + F_phase_independent
      dF_sol = dF_sol + F_phase_independent

      ! Change the units from (free energy per kT per ion) to (erg/g).
      kT = temp * kerg / (abar * amu)
      dF_liq = dF_liq * kT
      dF_sol = dF_sol * kT

      ! Produce a smoothed version of the phase transition to extract the latent heat
      call decide_phase(dF_liq, dF_sol, kT, temp, rho, dF, phase, latent_ddlnT, latent_ddlnRho)

      if (dbg) then
         write(*,*) 'Phase', phase%val
      end if

   end subroutine nonideal_corrections
   

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
      latent_S = -(differentiate_1(dF_blur)) ! S = -dF/dT

      ! T dS/dlnT = T^2 dS/dT 
      latent_ddlnT = differentiate_1(latent_S) *  pow2(temp)

      ! T dS/dlnRho = T Rho dS/dRho
      latent_ddlnRho = temp * rho * differentiate_2(latent_S)      

   end subroutine decide_phase

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
   function nonideal_corrections_phase(NMIX,AY,AZion,ACMI,min_gamma_for_solid, max_gamma_for_liquid,Skye_solid_mixing_rule,&
                                       temp,abar,GAME,RS,LIQSOL,Zmean, Z2mean, Z52, Z53, Z321) result(dF)
      ! Inputs
      integer, intent(in) :: NMIX
      integer, intent(in) :: LIQSOL
      real(dp), intent(in) :: AZion(:), ACMI(:), abar, AY(:), Zmean, Z2mean, Z52, Z53, Z321, min_gamma_for_solid, max_gamma_for_liquid
      type(auto_diff_real_2var_order3), intent(in) :: temp, GAME, RS
      character(len=128), intent(in) :: Skye_solid_mixing_rule

      ! Intermediates and constants
      integer :: i,j
      type(auto_diff_real_2var_order3) :: FMIX, f
      real(dp), parameter :: TINY=1.d-7

      ! Output
      type(auto_diff_real_2var_order3) :: dF

      dF = 0d0

      ! Composition loop
      do i=1,nmix
         if (AY(i) > TINY .and. AZion(i) /= 0d0) then ! skip low-abundance species and neutrons
            ! Add up non-ideal corrections
            f = extrapolate_free_energy(LIQSOL, temp, RS, AZion(i), ACMI(i), min_gamma_for_solid, max_gamma_for_liquid)
            dF = dF + AY(i) * f

         end if
      enddo

      ! Corrections to the linear mixing rule:
      if (LIQSOL == 0) then ! liquid phase
         FMIX = liquid_mixing_rule_correction(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321)
      else ! solid phase (only Madelung contribution) [22.12.12]
         FMIX = solid_mixing_rule_correction(Skye_solid_mixing_rule, NMIX, AY, AZion, GAME)
      endif
      dF = dF + FMIX

   end function nonideal_corrections_phase

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

      real(dp), parameter :: AUM = amu / me

      ! Output
      type(auto_diff_real_2var_order3) :: F

      GAMI = pow(Zion,5d0/3d0) * qe * qe / (rbohr * boltzm * temp * RS) ! ion Coulomb parameter Gamma_i
      COTPT=sqrt(3d0/AUM/CMI)/pow(Zion,7d0/6d0) ! auxiliary coefficient
      TPT=GAMI/sqrt(RS)*COTPT                   ! T_p/T

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
         g = pow(Zion,5d0/3d0) * qe * qe / (rbohr * boltzm * temp_boundary * RS) ! ion Coulomb parameter Gamma_i
         tp=g/sqrt(RS)*COTPT                   ! T_p/T

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
         dF_dlnT = make_binop(temp_boundary, fake_dens, dF_dlnT%val, dF_dlnT%d1val1, dF_dlnT%d1val2, dF_dlnT%d2val1, dF_dlnT%d1val1_d1val2, &
                     dF_dlnT%d2val2, dF_dlnT%d3val1, dF_dlnT%d2val1_d1val2, dF_dlnT%d1val1_d2val2, dF_dlnT%d3val2)

         F = make_binop(temp_boundary, fake_dens, F%val, F%d1val1, F%d1val2, F%d2val1, &
                     F%d1val1_d1val2, F%d2val2, F%d3val1,  &
                     F%d2val1_d1val2, F%d1val1_d2val2, F%d3val2)

         ! Extrapolate
         F = F + (1d0 - temp_boundary / temp) * dF_dlnT
      end if

   end function extrapolate_free_energy

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

      ! Intermediates
      type(auto_diff_real_2var_order3) :: screening_factor

      ! Output
      type(auto_diff_real_2var_order3) :: F

      screening_factor = pow2(TPT / GAMI) ! Proportional to RS_{ion}^{-1} ~ rho^{1/3}
      screening_factor = pow3(screening_factor / (1d-4 + screening_factor))

      if (LIQSOL == 0) then
         F = classical_ocp_liquid_free_energy(GAMI)                  ! classical ion-ion interaction
         F = F + quantum_ocp_liquid_free_energy_correction(TPT)   ! quantum ion-ion corrections
         F = F + screening_factor * ocp_liquid_screening_free_energy_correction(Zion, CMI*AMU, GAMI, TPT) ! screening corrections
      else     
         F = ocp_solid_harmonic_free_energy(GAMI,TPT) ! harmonic classical and quantum ion-ion corrections
         F = F + ocp_solid_anharmonic_free_energy(GAMI,TPT) ! anharmonic classical and quantum ion-ion corrections
         F = F + screening_factor * ocp_solid_screening_free_energy_correction(Zion, CMI*AMU, GAMI, TPT) ! screening corrections
      endif

   end function ocp_free_energy


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
      type(auto_diff_real_2var_order3) :: C, CDH, CDHH, C3, C3DH, C3DHH
      type(auto_diff_real_2var_order3) :: D0, D1, D, E0, E1, E
      type(auto_diff_real_2var_order3) :: DISCR, SQGE
      type(auto_diff_real_2var_order3) :: B2, B3, R3, S1, S2, S3, B4, C4
      type(auto_diff_real_2var_order3) :: S4A, S4B, S4C, S4

      ! Output
      type(auto_diff_real_2var_order3) :: FXC

      THETA=0.543d0*RS/GAME ! non-relativistic degeneracy parameter
      SQTH=sqrt(THETA)
      THETA2=THETA*THETA
      THETA3=THETA2*THETA
      THETA4=THETA3*THETA
      if (THETA.gt..007d0) then
         CHT1=cosh(1.d0/THETA)
         SHT1=sinh(1.d0/THETA)
         CHT2=cosh(1.d0/SQTH)
         SHT2=sinh(1.d0/SQTH)
         T1=SHT1/CHT1 ! dtanh(1.d0/THETA)
         T2=SHT2/CHT2 ! dtanh(1./sqrt(THETA))
      else
         T1=1.d0
         T2=1.d0
      endif

      A0=0.75d0+3.04363d0*THETA2-0.09227d0*THETA3+1.7035d0*THETA4
      A1=1d0+8.31051d0*THETA2+5.1105d0*THETA4
      A=0.610887d0*A0/A1*T1 ! HF fit of Perrot and Dharma-wardana

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

   end module skye_coulomb
