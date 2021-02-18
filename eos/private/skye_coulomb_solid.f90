module skye_coulomb_solid
   use math_lib
   use auto_diff
   use const_def

   implicit none

   contains

   !> Computes the classical and quantum anharmonic non-ideal
   !! free energy of a one-component plasma in the solid phase
   !! using classical fits due to
   !! Farouki & Hamaguchi'93.
   !!
   !! Quantum corrections to this are due to
   !! Potekhin & Chabrier 2010 (DOI 10.1002/ctpp.201010017)
   !!
   !! For a more recent reference and more mathematical detail see Appendix B.2 of
   !! A. Y. Potekhin and G. Chabrier: Equation of state for magnetized Coulomb plasmas
   !!
   !! Based on the routine ANHARM8 due to Potekhin and Chabrier.
   !!
   !! @param GAMI Ion coupling parameter (Gamma_i)
   !! @param TPT effective T_p/T - ion quantum parameter
   !! @param F free energy per kT per ion
   function ocp_solid_anharmonic_free_energy(GAMI,TPT) result(F)
      type(auto_diff_real_2var_order3), intent(in) :: GAMI,TPT
      type(auto_diff_real_2var_order3) :: F

      real(dp), parameter :: A1 = 10.9d0
      real(dp), parameter :: A2 = 247d0
      real(dp), parameter :: A3 = 1.765d5
      real(dp), parameter :: B1 = 0.12d0 ! coefficient of \eta^2/\Gamma at T=0

      type(auto_diff_real_2var_order3) :: TPT2

      TPT2=TPT*TPT

      F = -(A1 / GAMI + A2 / (2d0 * pow2(GAMI)) + A3 / (3d0 * pow3(GAMI)))
      F = F * exp(-(B1 / A1) * TPT2) ! suppress.factor of classical anharmonicity
      F = F - B1 * TPT2 / GAMI ! Quantum correction
   end function ocp_solid_anharmonic_free_energy

   !> Computes the harmonic non-ideal free energy of a
   !! one-component plasma in the solid phase using fits due to
   !! Baiko, Potekhin, & Yakovlev (2001). Includes both classical
   !! and quantum corrections.
   !!
   !! For a more recent reference and more mathematical detail see Appendix B.2 of
   !! A. Y. Potekhin and G. Chabrier: Equation of state for magnetized Coulomb plasmas
   !!
   !! Based on the routines HLfit12 and FHARM12 due to Potekhin and Chabrier.
   !!
   !! @param GAMI Ion coupling parameter (Gamma_i)
   !! @param TPT effective T_p/T - ion quantum parameter
   !! @param F free energy per kT per ion
   function ocp_solid_harmonic_free_energy(GAMI,TPT_in) result(F)
      ! Inputs
      type(auto_diff_real_2var_order3), intent(in) :: GAMI,TPT_in

      ! Intermediates
      type(auto_diff_real_2var_order3) :: TPT, UP, DN, EA, EB, EG, UP1, UP2, DN1, DN2, E0
      type(auto_diff_real_2var_order3) :: Fth, U0
      
      ! Output
      type(auto_diff_real_2var_order3) :: F

      real(dp), parameter :: CM = .895929256d0 ! Madelung
      real(dp), parameter :: EPS=1d-5

      ! BCC lattice.
      ! Empirically the Potekhin & Chabrier fit to the BCC lattice produces
      ! a lower free energy than their fit to the FCC lattice in all cases of
      ! interest (and no cases of which I am aware), so we can specialize
      ! to the BCC case.
      real(dp), parameter :: CLM=-2.49389d0 ! 3*ln<\omega/\omega_p>
      real(dp), parameter :: U1=0.5113875d0
      real(dp), parameter :: ALPHA=0.265764d0
      real(dp), parameter :: BETA=0.334547d0
      real(dp), parameter :: GAMMA=0.932446d0
      real(dp), parameter :: A1=0.1839d0
      real(dp), parameter :: A2=0.593586d0
      real(dp), parameter :: A3=0.0054814d0
      real(dp), parameter :: A4=5.01813d-4
      real(dp), parameter :: A6=3.9247d-7
      real(dp), parameter :: A8=5.8356d-11
      real(dp), parameter :: B0=261.66d0
      real(dp), parameter :: B2=7.07997d0
      real(dp), parameter :: B4=0.0409484d0
      real(dp), parameter :: B5=0.000397355d0
      real(dp), parameter :: B6=5.11148d-5
      real(dp), parameter :: B7=2.19749d-6
      real(dp), parameter :: C9=0.004757014d0
      real(dp), parameter :: C11=0.0047770935d0
      real(dp), parameter :: B9=A6*C9
      real(dp), parameter :: B11=A8*C11


      TPT = TPT_in

      if (TPT > 1d0/EPS) then ! asymptote of Eq.(13) of BPY'01
         F=-1d0 / (C11*TPT*TPT*TPT)
      else if (TPT < EPS) then ! Eq.(17) of BPY'01
         F=3d0*log(TPT)+CLM-1.5d0*U1*TPT+TPT*TPT/24.d0
      else
         UP=1d0+TPT*(A1+TPT*(A2+TPT*(A3+TPT*(A4+TPT*TPT*(A6+TPT*TPT*A8)))))
         DN=B0+TPT*TPT*(B2+TPT*TPT*(B4+TPT*(B5+TPT*(B6+TPT*(B7+TPT*TPT*(B9+TPT*TPT*B11))))))

         EA=exp(-ALPHA*TPT)
         EB=exp(-BETA*TPT)
         EG=exp(-GAMMA*TPT)

         F=log(1.d0-EA)+log(1.d0-EB)+log(1.d0-EG)-UP/DN ! F_{thermal}/NT
      end if

      U0=-CM*GAMI       ! perfect lattice
      E0=1.5d0*U1*TPT   ! zero-point energy
      Fth=F+E0          ! Thermal component
      F=U0+Fth          ! Total

      F = F -3d0 * log(TPT) + 1.5d0*log(GAMI) + 1.323515d0     ! Subtract out ideal free energy
                                                               ! We use this form because it's what PC fit against
                                                               ! in producing FHARM

   end function ocp_solid_harmonic_free_energy

   !> Calculates a log with cutoffs to prevent over/underflow.
   !!
   !! @param x Input to take the exponential of.
   !! @param ex Output (exp(x) clipped to avoid over/underflow).
   type(auto_diff_real_2var_order3) function safe_exp(x) result(ex)
      type(auto_diff_real_2var_order3), intent(in) :: x
      ex = exp(max(-3d1,min(3d1,x)))
   end function safe_exp

   !> Calculates the electron-ion screening corrections to the free energy
   !! of a one-component plasma in the solid phase using the fits of Potekhin & Chabrier 2013.
   !!
   !! @param Z ion charge
   !! @param mi ion mass in grams
   !! @param g ion interaction parameter
   !! @param TPT effective T_p/T - ion quantum parameter
   !! @param F non-ideal free energy
   function ocp_solid_screening_free_energy_correction(Z, mi, g, TPT) result(F)
         real(dp), intent(in) :: Z, mi
         type(auto_diff_real_2var_order3), intent(in) :: g, TPT

         real(dp) :: s, b1, b2, b3, b4
         real(dp), parameter :: e = 2.718281828459045
         real(dp), parameter :: aTF = 0.00352

         type(auto_diff_real_2var_order3) :: x, f_inf, A, Q, xr, eta, rs, supp
         type(auto_diff_real_2var_order3) :: F

         s = 1d0 / (1d0 + 1d-2 * pow(log(Z), 1.5d0) + 0.097d0 / pow2(Z))
         b1 = 1d0 - 1.1866d0 * pow(Z, -0.267d0) + 0.27d0 / Z
         b2 = 1d0 + (2.25d0 / pow(Z, 1d0/3d0)) * (1d0 + 0.684d0 * pow5(Z) + 0.222d0 * pow6(Z)) / (1d0 + 0.222d0 * pow6(Z))
         b3 = 41.5d0 / (1d0 + log(Z))
         b4 = 0.395d0 * log(Z) + 0.347d0 * pow(Z, -1.5d0)

         rs = (me / mi) * (3d0 * pow2(g / TPT)) * pow(Z, -7d0/3d0)
         xr = 0.014005d0 / rs
         supp = safe_exp(-pow2(0.205d0 * TPT))
         A = (b3 + 17.9d0 * pow2(xr)) / (1d0 + b4 * pow2(xr))
         Q = sqrt(log(1d0 + 1d0/supp)) / sqrt(log(e - (e - 2d0) * supp))
         f_inf = aTF * pow(Z, 2d0/3d0) * b1 * sqrt(1d0 + b2 / pow2(xr))

         F = -f_inf * g * (1d0 + A * pow(Q / g, s))

   end function ocp_solid_screening_free_energy_correction

   !> Computes the correction deltaG to the linear mixing rule for a two-component Coulomb solid mixture.
   !! From Shuji Ogata, Hiroshi Iyetomi, and Setsuo Ichimaru 1993
   !!
   !! Based on simulations done for charge ratio 4/3 <= Rz <= 4 and 163 <= Gamma <= 383, where
   !! Gamma here is the ionic interaction strength measured for the lower-charge species.
   !! deltaG is defined as deltaF / (x1*x2*Gamma), where deltaF is a free energy per kT per ion
   !! and x1 is the abundance of species 1.
   !!
   !! @param x2 Abundance of the higher-charge species
   !! @param Rz Charge ratio of species (> 1 by definition).
   real(dp) function deltaG_Ogata93(x2, Rz) result(dG)
      ! Inputs
      real(dp), intent(in) :: x2
      real(dp), intent(in) :: Rz

      ! Intermediates
      real(dp) :: CR

      CR = 0.05d0 * pow2(Rz - 1d0) / ((1d0 + 0.64d0 * (Rz - 1d0)) * (1d0 + 0.5d0 * pow2(Rz - 1d0)))
      dG = CR / (1 + (sqrt(x2) * (sqrt(x2) - 0.3d0) * (sqrt(x2) - 0.7d0) * (sqrt(x2) - 1d0)) * 27d0 * (Rz - 1d0) / (1d0 + 0.1d0 * (Rz - 1d0)))

   end function deltaG_Ogata93

   !> Computes the correction deltaG to the linear mixing rule for a two-component Coulomb solid mixture.
   !! Originally from PhysRevE.79.016411 (Equation of state of classical Coulomb plasma mixtures)
   !! by Potekhin, Alexander Y. and Chabrier, Gilles and Rogers, Forrest J.
   !! https://link.aps.org/doi/10.1103/PhysRevE.79.016411
   !!
   !! @param x Abundance of the higher-charge species
   !! @param Rz Charge ratio of species (> 1 by definition).
   real(dp) function deltaG_PC13(x2, Rz) result(dG)
      ! Inputs
      real(dp), intent(in) :: x2
      real(dp), intent(in) :: Rz

      ! Intermediates
      real(dp) :: x

      x = x2 / Rz + (1d0 - 1d0 / Rz) * pow(x2, Rz)
      dG = 0.012d0 * ((x*(1d0-x)) / (x2*(1d0-x2))) * (1d0 - 1d0/pow2(Rz)) * (1d0 - x2 + x2 * pow(Rz,5d0/3d0))

   end function deltaG_PC13

   !> Calculates the correction to the linear mixing rule for a Coulomb solid mixture
   !! by extending a two-component deltaG prescription to the multi-component case, using the
   !! prescription of Medin & Cumming 2010.
   !! 
   !! @param n Number of species
   !! @param AY Array of length NMIX holding the masses of species
   !! @param AZion Array of length NMIX holding the charges of species
   !! @param GAME election interaction parameter
   !! @param F mixing free energy correction per ion per kT.
   function solid_mixing_rule_correction(Skye_solid_mixing_rule, n, AY, AZion, GAME) result(F)      
      ! Inputs
      character(len=128), intent(in) :: Skye_solid_mixing_rule
      integer, intent(in) :: n
      real(dp), intent(in) :: AZion(:), AY(:)
      type(auto_diff_real_2var_order3), intent(in) :: GAME

      ! Intermediates
      integer :: i,j, num_unique_charges
      real(dp) :: unique_charges(n), charge_abundances(n)
      logical :: found
      integer :: found_index
      real(dp) :: RZ, aj, dG
      type(auto_diff_real_2var_order3) :: GAMI

      ! Output
      type(auto_diff_real_2var_order3) :: F

      ! Parameters
      real(dp), parameter :: C = 0.012d0
      real(dp), parameter :: eps = 1d-40

      ! Identify and group unique charges
      num_unique_charges = 0
      do i=1,n
         found = .false.

         do j=1,num_unique_charges
            if (unique_charges(j) == AZion(i)) then
               found = .true.
               found_index = j
               exit
            end if
         end do

         if (.not. found) then
            num_unique_charges = num_unique_charges + 1
            unique_charges(num_unique_charges) = AZion(i)
            charge_abundances(num_unique_charges) = AY(i)
         else
            charge_abundances(found_index) = charge_abundances(found_index) + AY(i)
         end if
      end do

      F = 0d0
      do i=1,num_unique_charges
         if (unique_charges(i) == 0d0) cycle
         do j=1,num_unique_charges
            if (unique_charges(j) == 0d0) cycle

            ! Expression needs R > 1.
            ! From PC2013: 'dF_sol = Sum_i Sum_{j>i} ... where the indices are arranged so that Z_j < Z_{j+1}'
            ! So if j > i then Z_j > Z_i, which means R = Z_j / Z_i > 1.
            ! We extend to the case of equality by grouping equal-charge species together, as above.
            if (unique_charges(j) < unique_charges(i)) cycle

            RZ = unique_charges(j)/unique_charges(i) ! Charge ratio
         
            ! max avoids divergence.
            ! The contribution to F scales as abundance_sum^2, so in cases where the max returns eps
            ! we don't care much about the error this incurs.
            aj = charge_abundances(j) / max(eps, charge_abundances(i) + charge_abundances(j))! = x2 / (x1 + x2) in MC10's language
            if (Skye_solid_mixing_rule == 'Ogata') then
               dG = deltaG_Ogata93(aj, RZ)
            else if (Skye_solid_mixing_rule == 'PC') then
               dG = deltaG_PC13(aj, RZ)
            else
               write(*,*) 'Error: Invalid choice for Skye_solid_mixing_rule.'
               stop
            end if

            GAMI=pow(unique_charges(i),5d0/3d0)*GAME
            F = F +  GAMI * (charge_abundances(i) * charge_abundances(j) * dG)
         end do
      end do
   end function solid_mixing_rule_correction


end module skye_coulomb_solid
