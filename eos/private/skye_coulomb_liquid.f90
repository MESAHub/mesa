module skye_coulomb_liquid
   use math_lib
   use auto_diff
   use const_def

   implicit none

   contains

   !> Calculates the free energy of a classical one-component
   !! plasma in the liquid phase using the fitting form of
   !! A.Y.Potekhin and G.Chabrier,Phys.Rev.E62,8554(2000)
   !! fit to the data of DeWitt & Slattery (1999).
   !! This choice ensures consistency with the quantum_ocp_liquid_free_energy_correction
   !! routine, which is based on the fits of Baiko & Yakolev 2019 (who were correcting relative
   !! to the DeWitt & Slattery fits).
   !!
   !! @param g ion degeneracy parameter
   !! @param F non-ideal free energy
   function classical_ocp_liquid_free_energy(g) result(F)
      type(auto_diff_real_2var_order3_array), intent(in) :: g
      type(auto_diff_real_2var_order3_array) :: FA, FB
      type(auto_diff_real_2var_order3_array) :: F

      real(dp), parameter :: SQ32=sqrt(3d0) / 2d0

      real(dp), parameter :: A1=-0.9070d0
      real(dp), parameter :: A2=0.62954d0
      real(dp), parameter :: A3 = -SQ32 - A1 / sqrt(A2)

      real(dp), parameter :: B1 = 4.56d-3
      real(dp), parameter :: B2 = 211.6d0
      real(dp), parameter :: B3 = -1d-4
      real(dp), parameter :: B4 = 4.62d-3
      
      FA = A1 * (sqrt(g * (A2 + g)) - A2 * log(sqrt(g / A2) +  sqrt(1 + g / A2))) &
            + 2d0 * A3 * (sqrt(g) - atan(sqrt(g)))
      FB = B1 * (g - B2 * log(1d0 + g / B2)) + 0.5d0 * B3 * log(1d0 + pow2(g) / B4)
      F = FA + FB 

   end function classical_ocp_liquid_free_energy


   !> Calculates the quantum corrections to the free energy of a 
   !! one-component plasma in the liquid phase using the fits due to 
   !! Baiko & Yakovlev 2019.
   !!
   !! @param TPT effective T_p/T - ion quantum parameter
   !! @param F non-ideal free energy
   function quantum_ocp_liquid_free_energy_correction(TPT) result(F)
         type(auto_diff_real_2var_order3_array), intent(in) :: TPT
         type(auto_diff_real_2var_order3_array) :: eta
         type(auto_diff_real_2var_order3_array) :: F

         real(dp), parameter :: Q1 = 5.994d0
         real(dp), parameter :: Q2 = 70.3d0
         real(dp), parameter :: Q4 = 22.7d0
         real(dp), parameter :: Q3 = (0.25d0 - Q1 / Q2) * Q4

         eta = TPT / sqrt(3d0)

         F = Q1 * eta - Q1 * Q2 * log(1d0 + eta / Q2) + 0.5d0 * Q3 * log(1d0 + pow2(eta) / Q4)
   end function quantum_ocp_liquid_free_energy_correction


   !> Calculates the correction to the linear mixing rule for a Coulomb liquid mixture.
   !! Based on CORMIX Version 02.07.09 by Potekhin and Chabrier
   !!
   !! @param RS Electron density parameter for component species
   !! @param GAME Electron coupling parameter (Gamma_i)
   !! @param Zmean Mean charge of ions.
   !! @param Z2mean Mean of squared charge of ions.
   !! @param Z52 Mean of Z^(5/2) for ions.
   !! @param Z53 Mean of Z^(5/3) for ions.
   !! @param Z321 Mean of Z*(1+Z)^(3/2) for ions.
   !! @param FMIX mixing free energy correction per ion per kT.
   function liquid_mixing_rule_correction(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321) result(FMIX)
      type(auto_diff_real_2var_order3_array), intent(in) :: RS,GAME
      real(dp), intent(in) :: Zmean,Z2mean,Z52,Z53,Z321

      type(auto_diff_real_2var_order3_array) :: GAMImean, Dif0, DifR, DifFDH, D
      type(auto_diff_real_2var_order3_array) :: P3, D0, GP, FMIX0, Q, R, GQ
      
      type(auto_diff_real_2var_order3_array) :: FMIX

      real(dp), parameter :: TINY = 1d-9
      
      GAMImean=GAME*Z53
      if (RS.lt.TINY) then ! OCP
         Dif0=Z52-sqrt(Z2mean*Z2mean*Z2mean/Zmean)
      else
         Dif0=Z321-sqrt((Z2mean+Zmean)*(Z2mean+Zmean)*(Z2mean+Zmean)/Zmean)
      endif
      DifR=Dif0/Z52
      DifFDH=Dif0*GAME*sqrt(GAME/3d0) ! F_DH - F_LM(DH)
      D=Z2mean/(Zmean*Zmean)
      if (abs(D-1.d0).lt.TINY) then ! no correction
         FMIX=0d0
         return
      endif
      P3=pow(D,-0.2d0)
      D0=(2.6d0*DifR+14d0*DifR*DifR*DifR)/(1.d0-P3)
      GP=D0*pow(GAMImean,P3)
      FMIX0=DifFDH/(1d0+GP)

      Q=D*D*0.0117d0
      R=1.5d0/P3-1d0
      GQ=Q*GP
      FMIX=FMIX0/pow(1d0+GQ,R)
   end function liquid_mixing_rule_correction

end module skye_coulomb_liquid

