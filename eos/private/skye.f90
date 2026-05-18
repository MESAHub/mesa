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

module skye
      use const_def, only: dp, crad, kerg, ln10, mp, avo
      use math_lib
      use auto_diff
      use eos_def
      use eos_timing, only: eos_timing_start, eos_timing_record_component, &
         eos_timing_record_skye, i_skye_dxa_ideal, i_skye_dxa_coul, i_skye_dxa_pack

      implicit none

      logical, parameter :: dbg = .false.


      private
      public :: Get_Skye_EOS_Results, Get_Skye_alfa, Get_Skye_alfa_simple, get_Skye_for_eosdt

      contains

      subroutine Get_Skye_alfa( &
            rq, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr, d_alfa_dabar, d_alfa_dzbar)
         use const_def, only: dp
         use eos_blend
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: logRho, logT, Z, abar, zbar
         real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         integer, intent(out) :: ierr
         real(dp), intent(out), optional :: d_alfa_dabar, d_alfa_dzbar

         logical :: contained
         type(auto_diff_real_2var_order1) :: p(2), blend, dist

         ! Blend parameters
         real(dp) :: skye_blend_width
         integer, parameter :: num_points = 8
         real(dp) :: bounds(8,2)
         real(dp) :: dbounds_dabar(8,2), dbounds_dzbar(8,2)
         real(dp) :: d_dist_dabar, d_dist_dzbar
         type (Helm_Table), pointer :: ht

         ierr = 0
         if (present(d_alfa_dabar)) d_alfa_dabar = 0d0
         if (present(d_alfa_dzbar)) d_alfa_dzbar = 0d0
         ht => eos_ht
         skye_blend_width = 0.1d0
         d_dist_dabar = 0d0
         d_dist_dzbar = 0d0
         dbounds_dabar = 0d0
         dbounds_dzbar = 0d0

         ! Avoid catastrophic loss of precision in HELM tables
         bounds(1,1) = ht% logdlo
         bounds(1,2) = 8.3d0

         ! Rough ionization temperature from Jermyn+2021 Equation 52 (treating denominator as ~1).
         ! We put a lower bound of logT=7.3 to ensure that solar models never use Skye.
         ! This is because the blend even in regions that are 99+% ionized produces noticeable
         ! kinks in the sound speed profile on a scale testable by the observations.
         bounds(2,1) = ht% logdlo
         bounds(2,2) = max(7.3d0,log10(1d5 * pow2(zbar))) + skye_blend_width
         if (zbar > 0d0 .and. bounds(2,2) > 7.3d0 + skye_blend_width) &
            dbounds_dzbar(2,2) = 2d0/(ln10*zbar)

         ! Rough ionization density from Jermyn+2021 Equation 53, dividing by 3 so we get closer to Dragons.
         ! Don't let the density get below rho = 200 g/cc so that solar model stays away from blend into Skye.
         bounds(3,1) = max(2.3d0,log10(abar * pow3(zbar))) + skye_blend_width
         bounds(3,2) = max(7.3d0,log10(1d5 * pow2(zbar))) + skye_blend_width
         if (abar > 0d0 .and. zbar > 0d0 .and. bounds(3,1) > 2.3d0 + skye_blend_width) then
            dbounds_dabar(3,1) = 1d0/(ln10*abar)
            dbounds_dzbar(3,1) = 3d0/(ln10*zbar)
         end if
         dbounds_dzbar(3,2) = dbounds_dzbar(2,2)

         ! HELM low-T bound
         bounds(4,1) = max(2.3d0,log10(abar * pow3(zbar))) + skye_blend_width
         bounds(4,2) = ht% logtlo
         dbounds_dabar(4,1) = dbounds_dabar(3,1)
         dbounds_dzbar(4,1) = dbounds_dzbar(3,1)

         ! Lower-right of (rho,T) plane
         bounds(5,1) = ht% logdhi
         bounds(5,2) = ht% logtlo

         ! Upper-right of (rho,T) plane
         bounds(6,1) = ht% logdhi
         bounds(6,2) = ht% logthi

         ! Avoid catastrophic loss of precision in HELM tables
         bounds(7,1) = 3d0 * ht% logthi + log10(abar * mp * crad / (3d0 * kerg * (zbar + 1d0))) - 6d0
         bounds(7,2) =  ht% logthi
         if (abar > 0d0 .and. zbar > -1d0) then
            dbounds_dabar(7,1) = 1d0/(ln10*abar)
            dbounds_dzbar(7,1) = -1d0/(ln10*(zbar + 1d0))
         end if

         ! Avoid catastrophic loss of precision in HELM tables
         bounds(8,1) = 3d0 * 8.3d0 + log10(abar * mp * crad / (3d0 * kerg * (zbar + 1d0))) - 6d0
         bounds(8,2) = 8.3d0
         dbounds_dabar(8,1) = dbounds_dabar(7,1)
         dbounds_dzbar(8,1) = dbounds_dzbar(7,1)

         ! Set up auto_diff point
         p(1) = logRho
         p(1)%d1val1 = 1d0
         p(2) = logT
         p(2)%d1val2 = 1d0

         contained = is_contained(num_points, bounds, p)
         dist = min_distance_to_polygon(num_points, bounds, p)
         if (present(d_alfa_dabar) .or. present(d_alfa_dzbar)) then
            call min_distance_to_polygon_partials( &
               num_points, bounds, dbounds_dabar, dbounds_dzbar, &
               [logRho, logT], d_dist_dabar, d_dist_dzbar)
         end if

         if (contained) then  ! Make distance negative for points inside the polygon
            dist = -dist
            d_dist_dabar = -d_dist_dabar
            d_dist_dzbar = -d_dist_dzbar
         end if

         if (dist%val > 0d0 .and. dist%val < skye_blend_width) then
            if (present(d_alfa_dabar)) d_alfa_dabar = d_dist_dabar/skye_blend_width
            if (present(d_alfa_dzbar)) d_alfa_dzbar = d_dist_dzbar/skye_blend_width
         end if

         dist = dist / skye_blend_width
         blend = max(dist, 0d0)
         blend = min(blend, 1d0)

         alfa = blend%val
         d_alfa_dlogRho = blend%d1val1
         d_alfa_dlogT = blend%d1val2

         contains

         subroutine min_distance_to_polygon_partials( &
               num_points, coords, dcoords_dabar, dcoords_dzbar, p, &
               d_dabar, d_dzbar)
            integer, intent(in) :: num_points
            real(dp), intent(in) :: coords(num_points,2)
            real(dp), intent(in) :: dcoords_dabar(num_points,2)
            real(dp), intent(in) :: dcoords_dzbar(num_points,2)
            real(dp), intent(in) :: p(2)
            real(dp), intent(out) :: d_dabar, d_dzbar

            integer :: i, i_next
            real(dp) :: line_dist, min_dist, line_dabar, line_dzbar

            min_dist = 1d99
            d_dabar = 0d0
            d_dzbar = 0d0
            do i = 1, num_points
               if (i == num_points) then
                  i_next = 1
               else
                  i_next = i + 1
               end if
               call distance_to_line_segment_partials( &
                  coords(i,1:2), coords(i_next,1:2), &
                  dcoords_dabar(i,1:2), dcoords_dabar(i_next,1:2), &
                  dcoords_dzbar(i,1:2), dcoords_dzbar(i_next,1:2), &
                  p, line_dist, line_dabar, line_dzbar)
               if (line_dist < min_dist) then
                  min_dist = line_dist
                  d_dabar = line_dabar
                  d_dzbar = line_dzbar
               end if
            end do

         end subroutine min_distance_to_polygon_partials


         subroutine distance_to_line_segment_partials( &
               line_start, line_end, dstart_dabar, dend_dabar, &
               dstart_dzbar, dend_dzbar, p, d, d_dabar, d_dzbar)
            real(dp), intent(in) :: line_start(2), line_end(2)
            real(dp), intent(in) :: dstart_dabar(2), dend_dabar(2)
            real(dp), intent(in) :: dstart_dzbar(2), dend_dzbar(2)
            real(dp), intent(in) :: p(2)
            real(dp), intent(out) :: d, d_dabar, d_dzbar

            real(dp), parameter :: eps = 1d-10
            real(dp) :: diff_start(2), diff_end(2), diff_line(2)
            real(dp) :: length_squared, lambda, nearest(2)

            diff_start = p - line_start
            diff_end = p - line_end
            diff_line = line_end - line_start
            length_squared = dot_product(diff_line, diff_line)
            if (length_squared <= 0d0) then
               d = sqrt(dot_product(diff_start, diff_start))
               call endpoint_partial(diff_start, dstart_dabar, d, d_dabar)
               call endpoint_partial(diff_start, dstart_dzbar, d, d_dzbar)
               return
            end if

            lambda = dot_product(diff_start, diff_line)/length_squared
            if (lambda < 0d0) then
               d = sqrt(dot_product(diff_start, diff_start))
               call endpoint_partial(diff_start, dstart_dabar, d, d_dabar)
               call endpoint_partial(diff_start, dstart_dzbar, d, d_dzbar)
            else if (lambda > 1d0) then
               d = sqrt(dot_product(diff_end, diff_end))
               call endpoint_partial(diff_end, dend_dabar, d, d_dabar)
               call endpoint_partial(diff_end, dend_dzbar, d, d_dzbar)
            else
               nearest = line_start + lambda*diff_line - p
               d = sqrt(dot_product(nearest, nearest) + eps*eps)
               call interior_partial( &
                  line_start, line_end, dstart_dabar, dend_dabar, p, &
                  lambda, nearest, d, d_dabar)
               call interior_partial( &
                  line_start, line_end, dstart_dzbar, dend_dzbar, p, &
                  lambda, nearest, d, d_dzbar)
            end if

         end subroutine distance_to_line_segment_partials


         subroutine endpoint_partial(diff, dpoint, d, dd)
            real(dp), intent(in) :: diff(2), dpoint(2), d
            real(dp), intent(out) :: dd

            if (d <= 0d0) then
               dd = 0d0
            else
               dd = -dot_product(diff, dpoint)/d
            end if

         end subroutine endpoint_partial


         subroutine interior_partial( &
               line_start, line_end, dstart, dend, p, lambda, nearest, d, dd)
            real(dp), intent(in) :: line_start(2), line_end(2), dstart(2), dend(2)
            real(dp), intent(in) :: p(2), lambda, nearest(2), d
            real(dp), intent(out) :: dd

            real(dp) :: diff_start(2), diff_line(2), dline(2)
            real(dp) :: length_squared, numerator, dnumerator, dlength_squared
            real(dp) :: dlambda, dnearest(2)

            if (d <= 0d0) then
               dd = 0d0
               return
            end if

            diff_start = p - line_start
            diff_line = line_end - line_start
            dline = dend - dstart
            length_squared = dot_product(diff_line, diff_line)
            numerator = dot_product(diff_start, diff_line)
            dnumerator = -dot_product(dstart, diff_line) + &
               dot_product(diff_start, dline)
            dlength_squared = 2d0*dot_product(diff_line, dline)
            dlambda = (dnumerator*length_squared - numerator*dlength_squared)/ &
               pow2(length_squared)

            dnearest = dstart + dlambda*diff_line + lambda*dline
            dd = dot_product(nearest, dnearest)/d

         end subroutine interior_partial

      end subroutine Get_Skye_alfa


      subroutine Get_Skye_alfa_simple( &
            rq, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr)
         use eos_blend
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: logRho, logT, Z, abar, zbar
         real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         integer, intent(out) :: ierr

         type(auto_diff_real_2var_order1) :: logT_auto, logRho_auto
         type(auto_diff_real_2var_order1) :: blend, blend_logT, blend_logRho

         include 'formats'

         ierr = 0

         ! logRho is val1
         logRho_auto% val = logRho
         logRho_auto% d1val1 = 1d0
         logRho_auto% d1val2 = 0d0

         ! logT is val2
         logT_auto% val = logT
         logT_auto% d1val1 = 0d0
         logT_auto% d1val2 = 1d0

         ! logT blend
         if (logT_auto < rq% logT_min_for_any_Skye) then
            blend_logT = 0d0
         else if (logT_auto <= rq% logT_min_for_all_Skye) then
            blend_logT = (logT_auto - rQ% logT_min_for_any_Skye) / (rq% logT_min_for_all_Skye - rq% logT_min_for_any_Skye)
         else if (logT_auto > rq% logT_min_for_all_Skye) then
            blend_logT = 1d0
         end if


         ! logRho blend
         if (logRho_auto < rq% logRho_min_for_any_Skye) then
            blend_logRho = 0d0
         else if (logRho_auto <= rq% logRho_min_for_all_Skye) then
            blend_logRho = (logRho_auto - rQ% logRho_min_for_any_Skye) / (rq% logRho_min_for_all_Skye - rq% logRho_min_for_any_Skye)
         else if (logRho_auto > rq% logRho_min_for_all_Skye) then
            blend_logRho = 1d0
         end if

         ! combine blends
         blend = (1d0 - blend_logRho) * (1d0 - blend_logT)

         alfa = blend% val
         d_alfa_dlogRho = blend% d1val1
         d_alfa_dlogT = blend% d1val2

      end subroutine get_Skye_alfa_simple


      subroutine get_Skye_for_eosdt( &
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

         skip = .false.
         if (include_composition_partials) call eos_timing_start(time0, clock_rate)
         call Get_Skye_EOS_Results( &
            rq, Z, X, abar, zbar, rho, logRho, T, logT, species, chem_id, xa, &
            res, d_dlnd, d_dlnT, d_dxa, ierr, include_composition_partials, dxa_rows, &
            dabar_dxa_in, dzbar_dxa_in, dz2bar_dxa_in, dz53bar_dxa_in, &
            dye_dxa_in, dmc_dxa_in)
         if (include_composition_partials) call eos_timing_record_component( &
            i_eos_Skye, time0, clock_rate)
         if (ierr /= 0) return

         ! zero all components
         res(i_frac:i_frac+num_eos_frac_results-1) = 0.0d0
         d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0.0d0
         d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0.0d0
         d_dxa(i_frac:i_frac+num_eos_frac_results-1,:) = 0.0d0

         ! mark this one
         res(i_frac_Skye) = 1.0d0

      end subroutine get_Skye_for_eosdt

      subroutine Get_Skye_EOS_Results( &
               rq, Z, X, abar, zbar, Rho, logRho, T, logT, &
               species, chem_id, xa, res, d_dlnd, d_dlnT, d_dxa, ierr, &
               include_composition_partials, dxa_rows, &
               dabar_dxa_in, dzbar_dxa_in, dz2bar_dxa_in, dz53bar_dxa_in, &
               dye_dxa_in, dmc_dxa_in)
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar
         real(dp), intent(in) :: Rho, logRho, T, logT
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)
         real(dp), intent(in) :: xa(:)
         integer, intent(out) :: ierr
         real(dp), intent(out), dimension(nv) :: res, d_dlnd, d_dlnT
         real(dp), intent(out), dimension(nv, species) :: d_dxa
         logical, intent(in), optional :: include_composition_partials
         integer, intent(in), optional :: dxa_rows(:)
         real(dp), intent(in), optional :: &
            dabar_dxa_in(:), dzbar_dxa_in(:), dz2bar_dxa_in(:), &
            dz53bar_dxa_in(:), dye_dxa_in(:), dmc_dxa_in(:)

         real(dp) :: logT_ion, logT_neutral
         logical :: do_composition_partials

         include 'formats'

         ierr = 0
         do_composition_partials = .false.
         if (present(include_composition_partials)) &
            do_composition_partials = include_composition_partials

         call skye_eos( &
            T, Rho, X, abar, zbar, &
            rq%Skye_min_gamma_for_solid, rq%Skye_max_gamma_for_liquid, &
            rq%Skye_solid_mixing_rule, rq%mass_fraction_limit_for_Skye, &
            rq%Skye_use_ion_offsets, &
            species, chem_id, xa, &
            res, d_dlnd, d_dlnT, d_dxa, ierr, do_composition_partials, dxa_rows, &
            dabar_dxa_in, dzbar_dxa_in, dz2bar_dxa_in, dz53bar_dxa_in, &
            dye_dxa_in, dmc_dxa_in)

         if (ierr /= 0) then
            if (dbg) then
               write(*,*) 'failed in Get_Skye_EOS_Results'
               write(*,1) 'T', T
               write(*,1) 'logT', logT
               write(*,1) 'Rho', Rho
               write(*,1) 'logRho', logRho
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'X', X
               call mesa_error(__FILE__,__LINE__,'Get_Skye_EOS_Results')
            end if
            return
         end if

      end subroutine Get_Skye_EOS_Results


      !>..given a temperature temp [K], density den [g/cm**3], and a composition
      !!..this routine returns most of the other
      !!..thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
      !!..specific thermal energy [erg/gr], the entropy [erg/g/K], along with
      !!..their derivatives with respect to temperature, density, abar, and zbar.
      !!..other quantities such the normalized chemical potential eta (plus its
      !!..derivatives), number density of electrons and positron pair (along
      !!..with their derivatives), adiabatic indices, specific heats, and
      !!..relativistically correct sound speed are also returned.
      !!..
      !!..this routine assumes planckian photons, an ideal gas of ions,
      !!..and an electron-positron gas with an arbitrary degree of relativity
      !!..and degeneracy. interpolation in a table of the helmholtz free energy
      !!..is used to return the electron-positron thermodynamic quantities.
      !!..all other derivatives are analytic.
      !!..
      !!..references: cox & giuli chapter 24 ; timmes & swesty apj 1999

      !!..this routine assumes a call to subroutine read_helm_table has
      !!..been performed prior to calling this routine.
         subroutine skye_eos( &
            temp_in, den_in, Xfrac, abar, zbar,  &
            Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
            Skye_solid_mixing_rule, &
            mass_fraction_limit, use_ion_offsets, &
            species, chem_id, xa, &
            res, d_dlnd, d_dlnT, d_dxa, ierr, include_composition_partials, dxa_rows, &
            dabar_dxa_in, dzbar_dxa_in, dz2bar_dxa_in, dz53bar_dxa_in, &
            dye_dxa_in, dmc_dxa_in)

         use eos_def
         use utils_lib, only: is_bad
         use chem_def, only: chem_isos
         use chem_lib, only: basic_composition_info
         use ion_offset, only: compute_ion_offset, compute_ion_offset_partials
         use eos_composition_partials, only: &
            get_eos_composition_partials, get_active_number_fraction_partial, &
            want_eos_dxa_row
         use skye_ideal
         use skye_coulomb
         use skye_thermodynamics
         use auto_diff

         integer :: i, j, row
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: temp_in, den_in, mass_fraction_limit, Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid
         real(dp), intent(in) :: Xfrac, abar, zbar
         logical, intent(in) :: use_ion_offsets
         logical, intent(in) :: include_composition_partials
         integer, intent(in), optional :: dxa_rows(:)
         real(dp), intent(in), optional :: &
            dabar_dxa_in(:), dzbar_dxa_in(:), dz2bar_dxa_in(:), &
            dz53bar_dxa_in(:), dye_dxa_in(:), dmc_dxa_in(:)
         character(len=128), intent(in) :: Skye_solid_mixing_rule
         integer, intent(out) :: ierr
         real(dp), intent(out), dimension(nv) :: res, d_dlnd, d_dlnT
         real(dp), intent(out), dimension(nv, species) :: d_dxa

         integer :: relevant_species, lookup(species), phase_hint
         integer :: time0, clock_rate
         type(auto_diff_real_2var_order3) :: temp, logtemp, den, logden, din
         real(dp) :: AZION(species), ACMI(species), A(species), select_xa(species), ya(species)
         type (Helm_Table), pointer :: ht
         real(dp) :: ytot1, ye, ye_tmp, norm, active_ytot
         real(dp) :: X_tmp, Y_tmp, Z_tmp, abar_tmp, zbar_tmp, z2bar, z53bar, mass_correction, sumx
         real(dp) :: dabar_dxa(species), dzbar_dxa(species), dz2bar_dxa(species)
         real(dp) :: dz53bar_dxa(species), dye_dxa(species), dmc_dxa(species)
         real(dp) :: doffset_dxa(species)
         real(dp) :: dya_dxa_j(species)
         type(auto_diff_real_2var_order3) :: etaele, xnefer, phase, latent_ddlnT, latent_ddlnRho
         type(auto_diff_real_2var_order3) :: F_ion_gas, F_rad, F_ideal_ion, F_coul
         type(auto_diff_real_2var_order3) :: F_ele, F_ele_dye, etaele_dye
         type(auto_diff_real_2var_order3) :: F_coul_dxa, dxnefer_dxa
         type(auto_diff_real_2var_order3) :: F_coul_dye, F_coul_dabar, F_coul_dya(species)
         type(auto_diff_real_2var_order3) :: F_coul_ocp(species)
         type(auto_diff_real_2var_order3) :: phase_x, latent_ddlnT_x, latent_ddlnRho_x
         type(auto_diff_real_2var_order3) :: phase_dye, phase_dabar, phase_dya(species)
         type(auto_diff_real_2var_order3) :: latent_ddlnT_dye, latent_ddlnT_dabar
         type(auto_diff_real_2var_order3) :: latent_ddlnT_dya(species)
         type(auto_diff_real_2var_order3) :: latent_ddlnRho_dye, latent_ddlnRho_dabar
         type(auto_diff_real_2var_order3) :: latent_ddlnRho_dya(species)
         type(auto_diff_real_2var_order3) :: F_ideal_ion_dxa(species), F_gas_dxa(species)
         type(auto_diff_real_2var_order3) :: &
            phase_dxa(species), latent_ddlnT_dxa(species), latent_ddlnRho_dxa(species)
         logical :: need_phase_dxa, need_latent_ddlnT_dxa, need_latent_ddlnRho_dxa

         ht => eos_ht

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

         temp = temp_in
         temp%d1val1 = 1d0
         logtemp = log10(temp)

         den = den_in
         den%d1val2 = 1d0
         logden = log10(den)

         ! HELM table lookup uses din rather than den
         ytot1 = 1.0d0 / abar
         ye = ytot1 * zbar
         din = ye*den

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
                  abar_tmp, zbar_tmp, z2bar, z53bar, ye_tmp, mass_correction, sumx)
               call get_eos_composition_partials( &
                  species, chem_id, abar, zbar, z2bar, z53bar, ye_tmp, mass_correction, sumx, &
                  dabar_dxa, dzbar_dxa, dz2bar_dxa, dz53bar_dxa, dye_dxa, dmc_dxa)
            end if
         end if

         F_rad = 0d0
         F_ion_gas = 0d0
         F_ideal_ion = 0d0
         F_coul = 0d0
         F_ele = 0d0

         ! Radiation free energy, independent of composition
         F_rad = compute_F_rad(temp, den)

         ! Count and pack relevant species for Coulomb corrections. Relevant means mass fraction above limit.
         relevant_species = 0
         norm = 0d0
         do j=1,species
            if (xa(j) > mass_fraction_limit) then
               relevant_species = relevant_species + 1
               AZION(relevant_species) = chem_isos% Z(chem_id(j))
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

         if (use_ion_offsets) then
            F_ideal_ion = F_ideal_ion + compute_ion_offset(species, xa, chem_id)  ! Offset so ion ground state energy is zero.
         end if

         ! Ideal electron-positron thermodynamics (s, e, p)
         ! Derivatives are handled by HELM code, so we don't pass *in* any auto_diff types (just get them as return values).
         call compute_ideal_ele(temp%val, den%val, din%val, logtemp%val, logden%val, zbar, ytot1, ye, ht, &
                               F_ele, F_ele_dye, etaele, etaele_dye, xnefer, ierr)

         xnefer = compute_xne(den, ytot1, zbar)

         ! Normalize mass fractions
         do j=1,relevant_species
            select_xa(j) = select_xa(j) / norm
         end do

         ! Compute non-ideal corrections
         if (include_composition_partials) then
            call nonideal_corrections(relevant_species, ya(1:relevant_species), &
                                      AZION(1:relevant_species), ACMI(1:relevant_species), &
                                      Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                                      Skye_solid_mixing_rule, den, temp, xnefer, abar, &
                                      F_coul, latent_ddlnT, latent_ddlnRho, phase, &
                                      F_coul_ocp(1:relevant_species))
         else
            call nonideal_corrections(relevant_species, ya(1:relevant_species), &
                                      AZION(1:relevant_species), ACMI(1:relevant_species), &
                                      Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                                      Skye_solid_mixing_rule, den, temp, xnefer, abar, &
                                      F_coul, latent_ddlnT, latent_ddlnRho, phase)
         end if

         call  pack_for_export(F_ideal_ion, F_coul, F_rad, F_ele, temp, den, xnefer, etaele, abar, zbar, &
                                 phase, latent_ddlnT, latent_ddlnRho, res, d_dlnd, d_dlnT, ierr)
         if(ierr/=0) return

         if (include_composition_partials) then
            need_phase_dxa = want_eos_dxa_row(i_phase, dxa_rows)
            need_latent_ddlnT_dxa = want_eos_dxa_row(i_latent_ddlnT, dxa_rows)
            need_latent_ddlnRho_dxa = want_eos_dxa_row(i_latent_ddlnRho, dxa_rows)
            phase_hint = -1
            if (.not. (need_phase_dxa .or. need_latent_ddlnT_dxa .or. need_latent_ddlnRho_dxa)) then
               ! Off the phase transition, the base EOS already picked the hard-min
               ! Coulomb branch; reuse it so dxa calls only evaluate that branch.
               if (phase% val == 0d0) phase_hint = 0
               if (phase% val == 1d0) phase_hint = 1
            end if
            call eos_timing_start(time0, clock_rate)
            call compute_F_ideal_ion_partials( &
               temp, den, abar, species, relevant_species, lookup, ACMI, A, ya, active_ytot, &
               dabar_dxa, F_ideal_ion_dxa)

            if (use_ion_offsets) then
               call compute_ion_offset_partials(species, xa, chem_id, doffset_dxa)
               do j = 1, species
                  F_ideal_ion_dxa(j) = F_ideal_ion_dxa(j) + doffset_dxa(j)
               end do
            end if
            call eos_timing_record_skye(i_skye_dxa_ideal, time0, clock_rate)

            call eos_timing_start(time0, clock_rate)
            dya_dxa_j(1:relevant_species) = 0d0
            if (phase_hint /= -1 .and. ye > 0d0) then
               ! Off the phase transition, dYe enters Coulomb terms only through
               ! xnefer, and dabar only through the kT/abar conversion.
               F_coul_dye = differentiate_2(F_coul)*den/ye
               F_coul_dabar = -F_coul/abar
            else if (need_phase_dxa .or. need_latent_ddlnT_dxa .or. need_latent_ddlnRho_dxa) then
               dxnefer_dxa = den*avo
               call nonideal_corrections_dxa( &
                  relevant_species, ya(1:relevant_species), dya_dxa_j(1:relevant_species), &
                  AZION(1:relevant_species), ACMI(1:relevant_species), &
                  Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                  Skye_solid_mixing_rule, den, temp, xnefer, dxnefer_dxa, &
                  abar, 0d0, F_coul_dye, &
                  latent_ddlnT_dxa=latent_ddlnT_dye, latent_ddlnRho_dxa=latent_ddlnRho_dye, &
                  phase_dxa=phase_dye)
               dxnefer_dxa = 0d0
               call nonideal_corrections_dxa( &
                  relevant_species, ya(1:relevant_species), dya_dxa_j(1:relevant_species), &
                  AZION(1:relevant_species), ACMI(1:relevant_species), &
                  Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                  Skye_solid_mixing_rule, den, temp, xnefer, dxnefer_dxa, &
                  abar, 1d0, F_coul_dabar, &
                  latent_ddlnT_dxa=latent_ddlnT_dabar, latent_ddlnRho_dxa=latent_ddlnRho_dabar, &
                  phase_dxa=phase_dabar)
            else
               dxnefer_dxa = den*avo
               call nonideal_corrections_dxa( &
                  relevant_species, ya(1:relevant_species), dya_dxa_j(1:relevant_species), &
                  AZION(1:relevant_species), ACMI(1:relevant_species), &
                  Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                  Skye_solid_mixing_rule, den, temp, xnefer, dxnefer_dxa, &
                  abar, 0d0, F_coul_dye, phase_hint=phase_hint)
               dxnefer_dxa = 0d0
               call nonideal_corrections_dxa( &
                  relevant_species, ya(1:relevant_species), dya_dxa_j(1:relevant_species), &
                  AZION(1:relevant_species), ACMI(1:relevant_species), &
                  Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                  Skye_solid_mixing_rule, den, temp, xnefer, dxnefer_dxa, &
                  abar, 1d0, F_coul_dabar, phase_hint=phase_hint)
            end if

            if (need_phase_dxa .or. need_latent_ddlnT_dxa .or. need_latent_ddlnRho_dxa) then
               do i = 1, relevant_species
                  dya_dxa_j(1:relevant_species) = 0d0
                  dya_dxa_j(i) = 1d0
                  call nonideal_corrections_dxa( &
                     relevant_species, ya(1:relevant_species), dya_dxa_j(1:relevant_species), &
                     AZION(1:relevant_species), ACMI(1:relevant_species), &
                     Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                     Skye_solid_mixing_rule, den, temp, xnefer, dxnefer_dxa, &
                     abar, 0d0, F_coul_dya(i), &
                     latent_ddlnT_dxa=latent_ddlnT_dya(i), &
                     latent_ddlnRho_dxa=latent_ddlnRho_dya(i), &
                     phase_dxa=phase_dya(i), composition_only_in_AY=.true.)
               end do
            else
               ! The star solve needs only lnPgas/lnE rows, so batch the d/dYA
               ! Coulomb basis and avoid recomputing the OCP sum for each species.
               call nonideal_corrections_dya( &
                  relevant_species, ya(1:relevant_species), &
                  AZION(1:relevant_species), ACMI(1:relevant_species), &
                  Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                  Skye_solid_mixing_rule, den, temp, xnefer, abar, &
                  F_coul_dya(1:relevant_species), phase_hint=phase_hint, &
                  phase_ocp=F_coul_ocp(1:relevant_species))
            end if

            do j = 1, species
               call get_active_number_fraction_partial( &
                  relevant_species, lookup(j), A, ya, active_ytot, dya_dxa_j(1:relevant_species))
               ! Preserve exact zero coefficients in the basis expansion.
               F_coul_dxa = 0d0
               if (dye_dxa(j) /= 0d0) F_coul_dxa = F_coul_dxa + dye_dxa(j)*F_coul_dye
               if (dabar_dxa(j) /= 0d0) F_coul_dxa = F_coul_dxa + dabar_dxa(j)*F_coul_dabar
               if (need_phase_dxa) then
                  phase_x = 0d0
                  if (dye_dxa(j) /= 0d0) phase_x = phase_x + dye_dxa(j)*phase_dye
                  if (dabar_dxa(j) /= 0d0) phase_x = phase_x + dabar_dxa(j)*phase_dabar
               end if
               if (need_latent_ddlnT_dxa) then
                  latent_ddlnT_x = 0d0
                  if (dye_dxa(j) /= 0d0) &
                     latent_ddlnT_x = latent_ddlnT_x + dye_dxa(j)*latent_ddlnT_dye
                  if (dabar_dxa(j) /= 0d0) &
                     latent_ddlnT_x = latent_ddlnT_x + dabar_dxa(j)*latent_ddlnT_dabar
               end if
               if (need_latent_ddlnRho_dxa) then
                  latent_ddlnRho_x = 0d0
                  if (dye_dxa(j) /= 0d0) &
                     latent_ddlnRho_x = latent_ddlnRho_x + dye_dxa(j)*latent_ddlnRho_dye
                  if (dabar_dxa(j) /= 0d0) &
                     latent_ddlnRho_x = latent_ddlnRho_x + dabar_dxa(j)*latent_ddlnRho_dabar
               end if
               do i = 1, relevant_species
                  if (dya_dxa_j(i) == 0d0) cycle
                  F_coul_dxa = F_coul_dxa + dya_dxa_j(i)*F_coul_dya(i)
                  if (need_phase_dxa) phase_x = phase_x + dya_dxa_j(i)*phase_dya(i)
                  if (need_latent_ddlnT_dxa) &
                     latent_ddlnT_x = latent_ddlnT_x + dya_dxa_j(i)*latent_ddlnT_dya(i)
                  if (need_latent_ddlnRho_dxa) &
                     latent_ddlnRho_x = latent_ddlnRho_x + dya_dxa_j(i)*latent_ddlnRho_dya(i)
               end do
               F_gas_dxa(j) = F_ideal_ion_dxa(j) + F_coul_dxa
               if (dye_dxa(j) /= 0d0) F_gas_dxa(j) = F_gas_dxa(j) + dye_dxa(j)*F_ele_dye
               if (need_phase_dxa) phase_dxa(j) = phase_x
               if (need_latent_ddlnT_dxa) latent_ddlnT_dxa(j) = latent_ddlnT_x
               if (need_latent_ddlnRho_dxa) latent_ddlnRho_dxa(j) = latent_ddlnRho_x
            end do
            call eos_timing_record_skye(i_skye_dxa_coul, time0, clock_rate)

            call eos_timing_start(time0, clock_rate)
            call pack_composition_partials( &
               F_gas_dxa, F_ideal_ion, F_coul, F_ele, temp, den, xnefer, abar, zbar, &
               dabar_dxa, dzbar_dxa, phase_dxa, latent_ddlnT_dxa, latent_ddlnRho_dxa, d_dxa, ierr, &
               dxa_rows)
            if(ierr/=0) then
               call eos_timing_record_skye(i_skye_dxa_pack, time0, clock_rate)
               return
            end if
            if (want_eos_dxa_row(i_eta, dxa_rows)) then
               do j = 1, species
                  d_dxa(i_eta,j) = etaele_dye%val*dye_dxa(j)
               end do
            end if
            call eos_timing_record_skye(i_skye_dxa_pack, time0, clock_rate)
            if (present(dxa_rows)) then
               do i = 1, size(dxa_rows)
                  row = dxa_rows(i)
                  if (row < 1 .or. row > nv) cycle
                  do j = 1, species
                     if (is_bad(d_dxa(row,j))) then
                        ierr = -1
                        return
                     end if
                  end do
               end do
            else
               do row = 1, nv
                  do j = 1, species
                     if (is_bad(d_dxa(row,j))) then
                        ierr = -1
                        return
                     end if
                  end do
               end do
            end if
         end if

      end subroutine skye_eos


end module skye
