module skye_thermodynamics
   use math_lib
   use auto_diff

   implicit none

   private

   public :: compute_derived_quantities, pack_for_export, thermodynamics_from_free_energy

   contains

   !> Computes the entropy, internal energy, and pressure of a system
   !! given its free energy, temperature, and density.
   !!
   !! @param F Free energy (erg/g)
   !! @param temp Temperature (K)
   !! @param den Density (g/cm^3)
   !! @param s Entropy (erg/g/K)
   !! @param e Internal energy (erg/g)
   !! @param p Pressure (erg/cm^3)
   subroutine thermodynamics_from_free_energy(F, temp, den, s, e, p)
      type(auto_diff_real_2var_order3), intent(in) :: F
      type(auto_diff_real_2var_order3), intent(in) :: temp, den
      type(auto_diff_real_2var_order3), intent(out) :: s, e, p

      ! Entropy
      ! s = - F/dT
      s = -differentiate_1(F)

      ! Pgas
      ! p = -dF/dV = -(1/V)dF/dlnV = (1/V)dF/dlnRho = Rho dF/dlnRho = Rho^2 dF/dRho
      p = pow2(den) * differentiate_2(F)

      ! Energy
      ! E = F + T * S
      e = F + temp * s

   end subroutine thermodynamics_from_free_energy

   !> Computes various thermodynamic derivatives and related quantities given
   !! the entropy, internal energy, pressure, temperature, and density of a system.
   !!
   !! @param temp Temperature (K)
   !! @param den Density (g/cm^3)
   !! @param s Entropy (erg/g/K)
   !! @param e Internal energy (erg/g)
   !! @param p Pressure (erg/cm^3)
   !! @param cv Specific heat at constant volume (erg/g/K)
   !! @param cp Specific heat at constant pressure (erg/g/K)
   !! @param chit Thermal susceptibility (dlnP/dlnT at constant Rho)
   !! @param chid Volume susceptibility (dlnP/dlnRho at constant T)
   !! @param gam1 First adiabatic index
   !! @param gam2 Second adiabatic index
   !! @param gam3 Third adiabatic index
   !! @param nabad Adiabatic gradient (dlnT/dlnP at constant entropy)
   !! @param cs Sound speed (cm/s)
   subroutine compute_derived_quantities(temp, dens, s, e, p, cv, cp, chit, chid, gam1, gam2, gam3, nabad, cs)
      use const_def, only: clight
      type(auto_diff_real_2var_order3), intent(in) :: temp, dens, s, e, p
      type(auto_diff_real_2var_order3), intent(out) :: cv, cp, chit, chid, gam1, gam2, gam3, nabad, cs

      ! Susceptibilities
      chit = differentiate_1(p) * temp / p
      chid = differentiate_2(p) * dens / p

      ! Specific heat at constant volume
      cv = differentiate_1(e)

      ! Adiabatic indices
      gam3 = 1d0 + (p / dens) * chit / (temp * cv)
      gam1 = chid + (gam3 - 1d0) * chit
      nabad = (gam3 - 1d0) / gam1
      gam2 = 1d0 - nabad

      ! Specific heat at constant pressure
      cp = cv * gam1 / chid

      ! Sound speed
      cs = clight * sqrt(gam1 / (1d0 + (dens / p) * (e + clight**2)))
   end subroutine compute_derived_quantities

   !> Computes thermodynamic quantities from Skye and packs them into the EOS return vectors.
   !!
   !! @param F_ideal_ion Ideal ion free energy (erg/g)
   !! @param F_coul Coulomb corrections to the free energy (erg/g)
   !! @param F_rad Radiation free energy (erg/g)
   !! @param F_ele Ideal electron-positron free energy (erg/g)
   !! @param temp Temperature (K)
   !! @param den Density (g/cm^3)
   !! @param xnefer Electron density (1/cm^3)
   !! @param etaele The electron chemical potential in units of kT.
   !! @param abar The mean atomic mass number.
   !! @param zbar The mean atomic charge number.
   !! @param abar The mean atomic mass number.
   !! @param phase The blended phase. 0 for liquid, 1 for solid, smoothly interpolates in between.
   !! @param latent_ddlnT The latent heat of the smoothed phase transition in lnT (T dS/dlnT)
   !! @param latent_ddlnRho The latent heat of the smoothed phase transition in lnRho (T dS/dlnRho)
   !! @param include_radiation True to include radiation effects, false otherwise.
   !! @param res The EOS return vector.
   !! @param d_dlnRho The derivative of the EOS return vector with respect to lnRho.
   !! @param d_dlnT The derivative of the EOS return vector with respect to lnT.
   subroutine pack_for_export(F_ideal_ion, F_coul, F_rad, F_ele, temp, dens, xnefer, etaele, abar, zbar, &
                                          phase, latent_ddlnT, latent_ddlnRho, res, d_dlnRho, d_dlnT)
      use eos_def
      type(auto_diff_real_2var_order3), intent(in) :: F_ideal_ion, F_coul, F_rad, F_ele, temp, dens, xnefer, etaele
      type(auto_diff_real_2var_order3), intent(in) :: phase, latent_ddlnT, latent_ddlnRho
      real(dp), intent(in) :: abar, zbar

      ! Intermediates
      type(auto_diff_real_2var_order3) :: srad, erad, prad, sgas, egas, pgas, p, e, s
      type(auto_diff_real_2var_order3) :: F_gas, lnS, lnE, lnPgas, mu, lnfree_e
      type(auto_diff_real_2var_order3) :: cv, cp, chit, chid, gam1, gam2, gam3, nabad, cs

      ! Outputs
      real(dp), intent(inout) :: res(nv)
      real(dp), intent(inout) :: d_dlnRho(nv)
      real(dp), intent(inout) :: d_dlnT(nv)

      ! Form the electron-positron-ion gas plus Coulomb corrections
      F_gas = F_ideal_ion + F_coul + F_ele

      ! Compute base thermodynamic quantities
      call thermodynamics_from_free_energy(F_gas, temp, dens, sgas, egas, pgas)
      call thermodynamics_from_free_energy(F_rad, temp, dens, srad, erad, prad)
      p = prad + pgas
      e = erad + egas
      s = srad + sgas

      lnS = log(s)
      lnE = log(e)
      lnPgas = log(pgas)
      
      ! assuming complete ionization
      mu = abar / (1d0 + zbar)
      lnfree_e = log(max(1d-99, xnefer)/(avo*dens))


      call compute_derived_quantities(temp, dens, s, e, p, cv, cp, chit, chid, gam1, gam2, gam3, nabad, cs)

      res(i_lnS) = lnS%val
      res(i_lnE) = lnE%val
      res(i_lnPgas) = lnPgas%val
      res(i_mu) = mu%val
      res(i_grad_ad) = nabad%val
      res(i_chiRho) = chid%val
      res(i_chiT) = chit%val
      res(i_Cp) = cp%val
      res(i_Cv) = cv%val
      res(i_dE_dRho) = e%d1val2
      res(i_dS_dT) = s%d1val1
      res(i_dS_dRho) = s%d1val2
      res(i_lnfree_e) = lnfree_e%val
      res(i_gamma1) = gam1%val
      res(i_gamma3) = gam3%val
      res(i_eta) = etaele%val
      res(i_phase) = phase%val
      res(i_latent_ddlnT) = latent_ddlnT%val
      res(i_latent_ddlnRho) = latent_ddlnRho%val

      d_dlnT(i_lnS) = lnS%d1val1 * temp%val
      d_dlnT(i_lnE) = lnE%d1val1 * temp%val
      d_dlnT(i_lnPgas) = lnPgas%d1val1 * temp%val
      d_dlnT(i_mu) = mu%d1val1 * temp%val
      d_dlnT(i_grad_ad) = nabad%d1val1 * temp%val
      d_dlnT(i_chiRho) = chid%d1val1 * temp%val
      d_dlnT(i_chiT) = chit%d1val1 * temp%val
      d_dlnT(i_Cp) = cp%d1val1 * temp%val
      d_dlnT(i_Cv) = cv%d1val1 * temp%val
      d_dlnT(i_dE_dRho) = e%d1val1_d1val2 * temp%val
      d_dlnT(i_dS_dT) = s%d2val1 * temp%val
      d_dlnT(i_dS_dRho) = s%d1val1_d1val2 * temp%val
      d_dlnT(i_lnfree_e) = lnfree_e%d1val1 * temp%val
      d_dlnT(i_gamma1) = gam1%d1val1 * temp%val
      d_dlnT(i_gamma3) = gam3%d1val1 * temp%val
      d_dlnT(i_eta) = etaele%d1val1 * temp%val
      d_dlnT(i_phase) = phase%d1val1 * temp%val
      d_dlnT(i_latent_ddlnT) = latent_ddlnT%d1val1 * temp%val
      d_dlnT(i_latent_ddlnRho) = latent_ddlnRho%d1val1 * temp%val     

      d_dlnRho(i_lnS) = lnS%d1val2 * dens%val
      d_dlnRho(i_lnE) = lnE%d1val2 * dens%val
      d_dlnRho(i_lnPgas) = lnPgas%d1val2 * dens%val
      d_dlnRho(i_mu) = mu%d1val2 * dens%val
      d_dlnRho(i_grad_ad) = nabad%d1val2 * dens%val
      d_dlnRho(i_chiRho) = chid%d1val2 * dens%val
      d_dlnRho(i_chiT) = chit%d1val2 * dens%val
      d_dlnRho(i_Cp) = cp%d1val2 * dens%val
      d_dlnRho(i_Cv) = cv%d1val2 * dens%val
      d_dlnRho(i_dE_dRho) = e%d2val2 * dens%val
      d_dlnRho(i_dS_dT) = s%d1val1_d1val2 * dens%val
      d_dlnRho(i_dS_dRho) = s%d2val2 * dens%val
      d_dlnRho(i_lnfree_e) = lnfree_e%d1val2 * dens%val
      d_dlnRho(i_gamma1) = gam1%d1val2 * dens%val
      d_dlnRho(i_gamma3) = gam3%d1val2 * dens%val
      d_dlnRho(i_eta) = etaele%d1val2 * dens%val
      d_dlnRho(i_phase) = phase%d1val2 * dens%val
      d_dlnRho(i_latent_ddlnT) = latent_ddlnT%d1val2 * dens%val
      d_dlnRho(i_latent_ddlnRho) = latent_ddlnRho%d1val2 * dens%val

   end subroutine pack_for_export



end module skye_thermodynamics
