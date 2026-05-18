! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
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

      module eos_composition_partials

      use const_def, only: dp
      use chem_def, only: chem_isos

      implicit none

      private
      public :: get_eos_composition_partials
      public :: get_active_number_fraction_partials
      public :: get_active_number_fraction_partial
      public :: want_eos_dxa_row

      contains

      logical function want_eos_dxa_row(row, dxa_rows)

         integer, intent(in) :: row
         integer, intent(in), optional :: dxa_rows(:)

         integer :: i

         want_eos_dxa_row = .true.
         if (.not. present(dxa_rows)) return

         want_eos_dxa_row = .false.
         do i = 1, size(dxa_rows)
            if (row == dxa_rows(i)) then
               want_eos_dxa_row = .true.
               return
            end if
         end do

      end function want_eos_dxa_row

      subroutine get_eos_composition_partials( &
            species, chem_id, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx, &
            dabar_dxa, dzbar_dxa, dz2bar_dxa, dz53bar_dxa, dye_dxa, dmc_dxa)

         integer, intent(in) :: species
         integer, pointer :: chem_id(:)
         real(dp), intent(in) :: abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
         real(dp), intent(out), dimension(species) :: &
            dabar_dxa, dzbar_dxa, dz2bar_dxa, dz53bar_dxa, dye_dxa, dmc_dxa

         integer :: j, cid
         real(dp) :: a, z, z2, z53, w, abar_over_a_sumx

         if (sumx <= 0d0) then
            dabar_dxa(:) = 0d0
            dzbar_dxa(:) = 0d0
            dz2bar_dxa(:) = 0d0
            dz53bar_dxa(:) = 0d0
            dye_dxa(:) = 0d0
            dmc_dxa(:) = 0d0
            return
         end if

         do j = 1, species
            cid = chem_id(j)
            a = dble(chem_isos% Z_plus_N(cid))
            z = dble(chem_isos% Z(cid))
            z2 = z*z
            z53 = chem_isos% Z53(cid)
            w = chem_isos% W(cid)
            abar_over_a_sumx = abar/(a*sumx)

            dabar_dxa(j) = abar_over_a_sumx*(a - abar)
            dzbar_dxa(j) = abar_over_a_sumx*(z - zbar)
            dz2bar_dxa(j) = abar_over_a_sumx*(z2 - z2bar)
            dz53bar_dxa(j) = abar_over_a_sumx*(z53 - z53bar)
            dye_dxa(j) = (z/a - ye)/sumx
            dmc_dxa(j) = w/a - mass_correction
         end do

      end subroutine get_eos_composition_partials


      subroutine get_active_number_fraction_partials( &
            species, relevant_species, lookup, aion, ya, active_ytot, dya_dxa)

         integer, intent(in) :: species, relevant_species
         integer, intent(in) :: lookup(species)
         real(dp), intent(in) :: aion(relevant_species), ya(relevant_species)
         real(dp), intent(in) :: active_ytot
         real(dp), intent(out) :: dya_dxa(relevant_species,species)

         integer :: i, j, i_active

         dya_dxa = 0d0
         if (active_ytot <= 0d0) return

         do j = 1, species
            i_active = lookup(j)
            if (i_active <= 0) cycle
            do i = 1, relevant_species
               dya_dxa(i,j) = -ya(i)/(active_ytot*aion(i_active))
               if (i == i_active) &
                  dya_dxa(i,j) = dya_dxa(i,j) + 1d0/(aion(i)*active_ytot)
            end do
         end do

      end subroutine get_active_number_fraction_partials


      subroutine get_active_number_fraction_partial( &
            relevant_species, i_active, aion, ya, active_ytot, dya_dxa)

         integer, intent(in) :: relevant_species, i_active
         real(dp), intent(in) :: aion(relevant_species), ya(relevant_species)
         real(dp), intent(in) :: active_ytot
         real(dp), intent(out) :: dya_dxa(relevant_species)

         integer :: i

         dya_dxa = 0d0
         if (active_ytot <= 0d0 .or. i_active <= 0) return

         do i = 1, relevant_species
            dya_dxa(i) = -ya(i)/(active_ytot*aion(i_active))
            if (i == i_active) &
               dya_dxa(i) = dya_dxa(i) + 1d0/(aion(i)*active_ytot)
         end do

      end subroutine get_active_number_fraction_partial


      end module eos_composition_partials
