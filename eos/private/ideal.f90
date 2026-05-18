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

module ideal

   use const_def, only: dp
   use math_lib
   use auto_diff
   use eos_def
   use eos_timing, only: eos_timing_start, eos_timing_record_component

   implicit none

   private
   public :: get_ideal_eos_results, get_ideal_alfa, get_ideal_for_eosdt

   contains

   subroutine get_ideal_alfa( &
            rq, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr)
      use const_def, only: dp
      use eos_blend
      type (EoS_General_Info), pointer :: rq
      real(dp), intent(in) :: logRho, logT, Z, abar, zbar
      real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
      integer, intent(out) :: ierr

      alfa = 0d0
      d_alfa_dlogRho = 0d0
      d_alfa_dlogT = 0d0
   end subroutine get_ideal_alfa

   subroutine get_ideal_for_eosdt( &
         handle, dbg, Z, X, abar, zbar, species, chem_id, net_iso, xa, &
         rho, logRho, T, logT, remaining_fraction, res, d_dlnd, d_dlnT, d_dxa, &
         include_composition_partials, skip, ierr, dxa_rows, &
         dabar_dxa_in, dzbar_dxa_in, dz2bar_dxa_in, dz53bar_dxa_in, &
         dye_dxa_in, dmc_dxa_in)
   integer, intent(in) :: handle
   logical, intent(in) :: dbg
   real(dp), intent(in) :: &
      Z, X, abar, zbar, remaining_fraction
   integer, intent(in) :: species
   integer, pointer :: chem_id(:), net_iso(:)
   real(dp), intent(in) :: xa(:)
   real(dp), intent(in) :: rho, logRho, T, logT
   real(dp), intent(inout), dimension(nv) :: res, d_dlnd, d_dlnT
   real(dp), intent(inout), dimension(nv, species) :: d_dxa
   logical, intent(in) :: include_composition_partials
   logical, intent(out) :: skip
   integer, intent(out) :: ierr
   integer, intent(in), optional :: dxa_rows(:)
   real(dp), intent(in), optional :: &
      dabar_dxa_in(:), dzbar_dxa_in(:), dz2bar_dxa_in(:), &
      dz53bar_dxa_in(:), dye_dxa_in(:), dmc_dxa_in(:)
   type (EoS_General_Info), pointer :: rq
   integer :: time0, clock_rate

   rq => eos_handles(handle)

   if (include_composition_partials) call eos_timing_start(time0, clock_rate)
   call get_ideal_eos_results(rq, Z, X, abar, zbar, rho, logRho, T, logT, species, chem_id, xa, &
                               res, d_dlnd, d_dlnT, d_dxa, ierr, include_composition_partials, dxa_rows, &
                               dabar_dxa_in, dzbar_dxa_in, dz2bar_dxa_in, dz53bar_dxa_in, &
                              dye_dxa_in, dmc_dxa_in)
   skip = .false.

   ! zero all components
   res(i_frac:i_frac+num_eos_frac_results-1) = 0.0d0
   d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0.0d0
   d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0.0d0
   d_dxa(i_frac:i_frac+num_eos_frac_results-1,:) = 0.0d0

   ! mark this one
   res(i_frac_ideal) = 1.0d0
   if (include_composition_partials) call eos_timing_record_component( &
      i_eos_ideal, time0, clock_rate)

   end subroutine get_ideal_for_eosdt

   subroutine get_ideal_eos_results( &
         rq, Z, X, abar, zbar, Rho, logRho, T, logT, &
         species, chem_id, xa, res, d_dlnd, d_dlnT, d_dxa, ierr, include_composition_partials, dxa_rows, &
         dabar_dxa_in, dzbar_dxa_in, dz2bar_dxa_in, dz53bar_dxa_in, dye_dxa_in, dmc_dxa_in)
   type (EoS_General_Info), pointer :: rq
   real(dp), intent(in) :: Z, X, abar, zbar
   real(dp), intent(in) :: Rho, logRho, T, logT
   integer, intent(in) :: species
   integer, pointer :: chem_id(:)
   real(dp), intent(in) :: xa(:)
   integer, intent(out) :: ierr
   real(dp), intent(out), dimension(nv) :: res, d_dlnd, d_dlnT
   real(dp), intent(out), dimension(nv, species) :: d_dxa
   logical, intent(in) :: include_composition_partials
   integer, intent(in), optional :: dxa_rows(:)
   real(dp), intent(in), optional :: &
      dabar_dxa_in(:), dzbar_dxa_in(:), dz2bar_dxa_in(:), &
      dz53bar_dxa_in(:), dye_dxa_in(:), dmc_dxa_in(:)

   ierr = 0

   call ideal_eos( &
      T, Rho, X, abar, zbar, &
      species, chem_id, xa, &
      res, d_dlnd, d_dlnT, d_dxa, ierr, include_composition_partials, dxa_rows, &
      dabar_dxa_in, dzbar_dxa_in, dz2bar_dxa_in, dz53bar_dxa_in, dye_dxa_in, dmc_dxa_in)

   end subroutine get_ideal_eos_results

   subroutine ideal_eos( &
         temp_in, den_in, Xfrac, abar, zbar,  &
         species, chem_id, xa, &
         res, d_dlnd, d_dlnT, d_dxa, ierr, include_composition_partials, dxa_rows, &
         dabar_dxa_in, dzbar_dxa_in, dz2bar_dxa_in, dz53bar_dxa_in, dye_dxa_in, dmc_dxa_in)

      use eos_def
      use const_def, only: dp
      use utils_lib, only: is_bad
      use chem_def, only: chem_isos
      use chem_lib, only: basic_composition_info
      use eos_composition_partials, only: get_eos_composition_partials
      use ion_offset, only: compute_ion_offset
      use skye_ideal
      use skye_thermodynamics
      use auto_diff

      integer :: j
      integer, intent(in) :: species
      integer, pointer :: chem_id(:)
      real(dp), intent(in) :: xa(:)
      real(dp), intent(in) :: temp_in, den_in
      real(dp), intent(in) :: Xfrac, abar, zbar
      logical, intent(in) :: include_composition_partials
      integer, intent(in), optional :: dxa_rows(:)
      real(dp), intent(in), optional :: &
         dabar_dxa_in(:), dzbar_dxa_in(:), dz2bar_dxa_in(:), &
         dz53bar_dxa_in(:), dye_dxa_in(:), dmc_dxa_in(:)
      integer, intent(out) :: ierr
      real(dp), intent(out), dimension(nv) :: res, d_dlnd, d_dlnT
      real(dp), intent(out), dimension(nv, species) :: d_dxa
      real(dp), parameter :: mass_fraction_limit = 1d-10

      integer :: relevant_species, lookup(species)
      type(auto_diff_real_2var_order3) :: temp, den
      real(dp) :: ACMI(species), A(species), ya(species), select_xa(species), norm
      real(dp) :: active_ytot, X_tmp, Y_tmp, Z_tmp, abar_tmp, zbar_tmp, z2bar, z53bar
      real(dp) :: ye, mass_correction, sumx
      real(dp) :: dabar_dxa(species), dzbar_dxa(species), dz2bar_dxa(species)
      real(dp) :: dz53bar_dxa(species), dye_dxa(species), dmc_dxa(species)
      type(auto_diff_real_2var_order3) :: etaele, xnefer, phase, latent_ddlnT, latent_ddlnRho
      type(auto_diff_real_2var_order3) :: F_rad, F_ideal_ion, F_ele, F_coul
      type(auto_diff_real_2var_order3) :: F_ideal_ion_dxa(species), F_gas_dxa(species)
      type(auto_diff_real_2var_order3) :: &
         phase_dxa(species), latent_ddlnT_dxa(species), latent_ddlnRho_dxa(species)

      ! Set up
      ierr = 0
      if (present(dxa_rows)) then
         do j = 1, size(dxa_rows)
            if (dxa_rows(j) < 1 .or. dxa_rows(j) > nv) cycle
            d_dxa(dxa_rows(j),:) = 0d0
         end do
      else
         d_dxa = 0d0
      end if
      lookup = 0
      F_ele = 0d0
      F_coul = 0d0

      ! No electrons, so extreme negative chemical potential
      etaele = -1d99
      xnefer = 1d-20

      ! no latent heat, no phase
      latent_ddlnT = 0d0
      latent_ddlnRho = 0d0
      phase = 0d0

      ! Construct rho,T and partials
      temp = temp_in
      temp%d1val1 = 1d0
      den = den_in
      den%d1val2 = 1d0

      if (include_composition_partials) then
         if (present(dabar_dxa_in) .and. present(dzbar_dxa_in) .and. &
               present(dz2bar_dxa_in) .and. present(dz53bar_dxa_in) .and. &
               present(dye_dxa_in) .and. present(dmc_dxa_in)) then
            dabar_dxa = dabar_dxa_in
            dzbar_dxa = dzbar_dxa_in
            dz2bar_dxa = dz2bar_dxa_in
            dz53bar_dxa = dz53bar_dxa_in
            dye_dxa = dye_dxa_in
            dmc_dxa = dmc_dxa_in
         else
            call basic_composition_info( &
               species, chem_id, xa, X_tmp, Y_tmp, Z_tmp, &
               abar_tmp, zbar_tmp, z2bar, z53bar, ye, mass_correction, sumx)
            call get_eos_composition_partials( &
               species, chem_id, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx, &
               dabar_dxa, dzbar_dxa, dz2bar_dxa, dz53bar_dxa, dye_dxa, dmc_dxa)
         end if
      end if

      ! Count and pack relevant species for Coulomb corrections. Relevant means mass fraction above limit.
      relevant_species = 0
      norm = 0d0
      do j=1,species
         if (xa(j) > mass_fraction_limit) then
            relevant_species = relevant_species + 1
            ACMI(relevant_species) = chem_isos% W(chem_id(j))
            A(relevant_species) = chem_isos% Z_plus_N(chem_id(j))
            select_xa(relevant_species) = xa(j)
            lookup(j) = relevant_species
            norm = norm + xa(j)
         end if
      end do

      active_ytot = 0d0
      do j=1,relevant_species
         active_ytot = active_ytot + select_xa(j)/A(j)
      end do

      ! Normalize
      do j=1,relevant_species
         select_xa(j) = select_xa(j) / norm
      end do

      ! Compute number fractions
      norm = 0d0
      do j=1,relevant_species
         ya(j) = select_xa(j) / A(j)
         norm = norm + ya(j)
      end do
      do j=1,relevant_species
         ya(j) = ya(j) / norm
      end do

      ! Ideal ion free energy, only depends on abar
      F_ideal_ion = compute_F_ideal_ion(temp, den, abar, relevant_species, ACMI, ya)
      if (include_composition_partials) then
         call compute_F_ideal_ion_partials( &
            temp, den, abar, species, relevant_species, lookup, ACMI, A, ya, active_ytot, &
            dabar_dxa, F_ideal_ion_dxa)
      end if

      ! The compute_ion_offset correction should only be applied to ionized matter.
      ! Right now, ideal always assumes neutral, so leave this off unless/until we consider
      ! adding assumption of full ionizatino for high T regions on ideal gas.
      !F_ideal_ion = F_ideal_ion + compute_ion_offset(species, xa, chem_id) ! Offset so ion ground state energy is zero.

      ! Radiation free energy, independent of composition
      F_rad = compute_F_rad(temp, den)


      call  pack_for_export(F_ideal_ion, F_coul, F_rad, F_ele, temp, den, xnefer, etaele, abar, zbar, &
                        phase, latent_ddlnT, latent_ddlnRho, res, d_dlnd, d_dlnT, ierr)
      if(ierr/=0) return

      res(i_mu) = abar  ! ideal assumes neutral matter, whereas pack_for_export assumes ionized matter. So we patch it up here.
      if (include_composition_partials) then
         do j = 1, species
            F_gas_dxa(j) = F_ideal_ion_dxa(j)
            phase_dxa(j) = 0d0
            latent_ddlnT_dxa(j) = 0d0
            latent_ddlnRho_dxa(j) = 0d0
         end do
         call pack_composition_partials( &
            F_gas_dxa, F_ideal_ion, F_coul, F_ele, temp, den, xnefer, abar, zbar, &
            dabar_dxa, dzbar_dxa, phase_dxa, latent_ddlnT_dxa, latent_ddlnRho_dxa, d_dxa, ierr, &
            dxa_rows)
         if(ierr/=0) return
         d_dxa(i_mu,:) = dabar_dxa
         d_dxa(i_lnfree_e,:) = 0d0
      end if

   end subroutine ideal_eos

end module ideal
