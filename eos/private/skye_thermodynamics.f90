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

module skye_thermodynamics
   use math_lib
   use auto_diff

   implicit none

   private

   public :: compute_derived_quantities, pack_for_export, pack_composition_partials, &
      thermodynamics_from_free_energy

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
      cs = clight * sqrt(gam1 / (1d0 + (dens / p) * (e + clight*clight)))
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
                                          phase, latent_ddlnT, latent_ddlnRho, res, d_dlnRho, d_dlnT, ierr)
      use const_def, only: dp, crad, avo
      use eos_def
      type(auto_diff_real_2var_order3), intent(in) :: F_ideal_ion, F_coul, F_rad, F_ele, temp, dens, xnefer, etaele
      type(auto_diff_real_2var_order3), intent(in) :: phase, latent_ddlnT, latent_ddlnRho
      real(dp), intent(in) :: abar, zbar

      integer, intent(out) :: ierr

      ! Intermediates
      type(auto_diff_real_2var_order3) :: srad, erad, prad, sgas, egas, pgas, p, e, s
      type(auto_diff_real_2var_order3) :: F_gas, lnS, lnE, lnPgas, mu, lnfree_e
      type(auto_diff_real_2var_order3) :: cv, cp, chit, chid, gam1, gam2, gam3, nabad, cs

      ! Outputs
      real(dp), intent(inout) :: res(nv)
      real(dp), intent(inout) :: d_dlnRho(nv)
      real(dp), intent(inout) :: d_dlnT(nv)

      ierr = 0

      ! Form the electron-positron-ion gas plus Coulomb corrections
      F_gas = F_ideal_ion + F_coul + F_ele

      ! Compute base thermodynamic quantities
      call thermodynamics_from_free_energy(F_gas, temp, dens, sgas, egas, pgas)

      ! Write the radiation terms explicitly to avoid having rho^2/rho^2 term that
      ! auto_diff gives for prad when calling thermodynamics_from_free_energy on F_rad.
      ! This avoids some subtractions for quantities that should be 0 like chid.
      prad = crad * pow4(temp) / 3d0
      erad = crad * pow4(temp) / dens
      srad = 4d0 * crad * pow3(temp) / (3d0 * dens)

      p = prad + pgas
      e = erad + egas
      s = srad + sgas

      if(s<0 .or. e<0 .or. pgas<0) then
         ierr = -1
         return
      end if

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


   subroutine pack_composition_partials( &
         F_gas_x, F_ideal_ion, F_coul, F_ele, temp, dens, xnefer, abar, zbar, &
         dabar_dxa, dzbar_dxa, phase_x, latent_ddlnT_x, latent_ddlnRho_x, d_dxa, ierr, &
         dxa_rows)
      use const_def, only: dp, crad
      use eos_def
      use eos_composition_partials, only: want_eos_dxa_row
      type(auto_diff_real_2var_order3), intent(in) :: F_gas_x(:)
      type(auto_diff_real_2var_order3), intent(in) :: F_ideal_ion, F_coul, F_ele
      type(auto_diff_real_2var_order3), intent(in) :: temp, dens, xnefer
      real(dp), intent(in) :: abar, zbar, dabar_dxa(:), dzbar_dxa(:)
      type(auto_diff_real_2var_order3), intent(in) :: phase_x(:), latent_ddlnT_x(:), latent_ddlnRho_x(:)
      real(dp), intent(inout) :: d_dxa(:,:)
      integer, intent(out) :: ierr
      integer, intent(in), optional :: dxa_rows(:)

      integer :: j, species
      real(dp) :: dmu
      type(auto_diff_real_2var_order3) :: srad, erad, prad, sgas, egas, pgas, p, e, s
      type(auto_diff_real_2var_order3) :: F_gas, cv, cp, chit, chid, gam1, gam2, gam3, nabad, cs
      type(auto_diff_real_2var_order3) :: s_x, e_x, p_x, cv_x, chiT_x, chiRho_x
      type(auto_diff_real_2var_order3) :: gamma1_x, gamma3_x, grad_ad_x, cp_x
      type(auto_diff_real_2var_order3) :: de_x_drho, ds_x_dT, ds_x_drho
      logical :: need_lnS, need_lnE, need_lnPgas, need_derived
      logical :: need_chiRho, need_chiT, need_Cp, need_Cv, need_dE_dRho
      logical :: need_dS_dT, need_dS_dRho, need_gamma1, need_gamma3, need_grad_ad
      logical :: need_mu, need_lnfree_e, need_eta, need_phase, need_latent_T, need_latent_Rho
      logical :: need_s_x, need_e_x, need_p_x

      ierr = 0
      species = size(d_dxa, dim=2)
      if (present(dxa_rows)) then
         do j = 1, size(dxa_rows)
            if (dxa_rows(j) < 1 .or. dxa_rows(j) > size(d_dxa, dim=1)) cycle
            d_dxa(dxa_rows(j),:) = 0d0
         end do
      else
         d_dxa(:,:) = 0d0
      end if

      need_lnS = want_eos_dxa_row(i_lnS, dxa_rows)
      need_lnE = want_eos_dxa_row(i_lnE, dxa_rows)
      need_lnPgas = want_eos_dxa_row(i_lnPgas, dxa_rows)
      need_chiRho = want_eos_dxa_row(i_chiRho, dxa_rows)
      need_chiT = want_eos_dxa_row(i_chiT, dxa_rows)
      need_Cp = want_eos_dxa_row(i_Cp, dxa_rows)
      need_Cv = want_eos_dxa_row(i_Cv, dxa_rows)
      need_dE_dRho = want_eos_dxa_row(i_dE_dRho, dxa_rows)
      need_dS_dT = want_eos_dxa_row(i_dS_dT, dxa_rows)
      need_dS_dRho = want_eos_dxa_row(i_dS_dRho, dxa_rows)
      need_gamma1 = want_eos_dxa_row(i_gamma1, dxa_rows)
      need_gamma3 = want_eos_dxa_row(i_gamma3, dxa_rows)
      need_grad_ad = want_eos_dxa_row(i_grad_ad, dxa_rows)
      need_mu = want_eos_dxa_row(i_mu, dxa_rows)
      need_lnfree_e = want_eos_dxa_row(i_lnfree_e, dxa_rows)
      need_eta = want_eos_dxa_row(i_eta, dxa_rows)
      need_phase = want_eos_dxa_row(i_phase, dxa_rows)
      need_latent_T = want_eos_dxa_row(i_latent_ddlnT, dxa_rows)
      need_latent_Rho = want_eos_dxa_row(i_latent_ddlnRho, dxa_rows)
      need_derived = need_chiRho .or. need_chiT .or. need_Cp .or. need_Cv .or. &
         need_dE_dRho .or. need_dS_dT .or. need_dS_dRho .or. &
         need_gamma1 .or. need_gamma3 .or. need_grad_ad
      need_s_x = need_lnS .or. need_lnE .or. need_derived
      need_e_x = need_lnE .or. need_derived
      need_p_x = need_lnPgas .or. need_derived

      F_gas = F_ideal_ion + F_coul + F_ele
      call thermodynamics_from_free_energy(F_gas, temp, dens, sgas, egas, pgas)

      prad = crad * pow4(temp) / 3d0
      erad = crad * pow4(temp) / dens
      srad = 4d0 * crad * pow3(temp) / (3d0 * dens)

      p = prad + pgas
      e = erad + egas
      s = srad + sgas

      if(s<0 .or. e<0 .or. pgas<0) then
         ierr = -1
         return
      end if

      if (need_derived) &
         call compute_derived_quantities(temp, dens, s, e, p, cv, cp, chit, chid, gam1, gam2, gam3, nabad, cs)

      do j = 1, species
         if (need_s_x) s_x = -differentiate_1(F_gas_x(j))
         if (need_p_x) p_x = pow2(dens) * differentiate_2(F_gas_x(j))
         if (need_e_x) e_x = F_gas_x(j) + temp*s_x

         if (need_lnS) d_dxa(i_lnS,j) = s_x%val/s%val
         if (need_lnE) d_dxa(i_lnE,j) = e_x%val/e%val
         if (need_lnPgas) d_dxa(i_lnPgas,j) = p_x%val/pgas%val

         if (need_derived) then
            cv_x = differentiate_1(e_x)
            chiT_x = temp*differentiate_1(p_x)/p - chit*p_x/p
            chiRho_x = dens*differentiate_2(p_x)/p - chid*p_x/p
            gamma3_x = (gam3 - 1d0)*(p_x/p + chiT_x/chit - cv_x/cv)
            gamma1_x = chiRho_x + gamma3_x*chit + (gam3 - 1d0)*chiT_x
            grad_ad_x = (gamma3_x*gam1 - (gam3 - 1d0)*gamma1_x)/pow2(gam1)
            cp_x = cp*(cv_x/cv + gamma1_x/gam1 - chiRho_x/chid)
            de_x_drho = differentiate_2(e_x)
            ds_x_dT = differentiate_1(s_x)
            ds_x_drho = differentiate_2(s_x)

            if (need_chiRho) d_dxa(i_chiRho,j) = chiRho_x%val
            if (need_chiT) d_dxa(i_chiT,j) = chiT_x%val
            if (need_Cp) d_dxa(i_Cp,j) = cp_x%val
            if (need_Cv) d_dxa(i_Cv,j) = cv_x%val
            if (need_dE_dRho) d_dxa(i_dE_dRho,j) = de_x_drho%val
            if (need_dS_dT) d_dxa(i_dS_dT,j) = ds_x_dT%val
            if (need_dS_dRho) d_dxa(i_dS_dRho,j) = ds_x_drho%val
            if (need_gamma1) d_dxa(i_gamma1,j) = gamma1_x%val
            if (need_gamma3) d_dxa(i_gamma3,j) = gamma3_x%val
            if (need_grad_ad) d_dxa(i_grad_ad,j) = grad_ad_x%val
         end if

         if (need_mu) then
            dmu = (dabar_dxa(j)*(1d0 + zbar) - abar*dzbar_dxa(j))/pow2(1d0 + zbar)
            d_dxa(i_mu,j) = dmu
         end if
         if (need_lnfree_e .and. xnefer%val > 1d-99 .and. zbar > 0d0 .and. abar > 0d0) then
            d_dxa(i_lnfree_e,j) = dzbar_dxa(j)/zbar - dabar_dxa(j)/abar
         end if

         if (need_eta) d_dxa(i_eta,j) = 0d0
         if (need_phase) d_dxa(i_phase,j) = phase_x(j)%val
         if (need_latent_T) d_dxa(i_latent_ddlnT,j) = latent_ddlnT_x(j)%val
         if (need_latent_Rho) d_dxa(i_latent_ddlnRho,j) = latent_ddlnRho_x(j)%val
      end do

   end subroutine pack_composition_partials



end module skye_thermodynamics
