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

module ion_offset
      use const_def, only: dp, amu, ev2erg
      use math_lib

      implicit none

      logical, parameter :: dbg = .false.

      private
      public :: compute_ion_offset

      contains

      !> Computes the energy required to fully ionize
      !! the specified composition. This is an offset
      !! between Skye/HELM and FreeEOS.
      !!
      !! @param species Number of species
      !! @param xa Mass fractions (vector)
      !! @param offset Energy in erg/g
      real(dp) function compute_ion_offset(species, xa, chem_id) result(offset)
         use chem_def, only: chem_isos
         integer, intent(in) :: species
         real(dp), intent(in) :: xa(species)
         integer, pointer :: chem_id(:)

         ! First 28 ionization energies, same species as used in FreeEOS, obtained from
         ! NIST (https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html)
         ! using the query 'H-Ni' asking for 'Total binding energy'.
         real(dp), parameter :: ionization_table(28) = [13.598434599702d0,79.0051545d0,203.4861711d0,399.14864d0,670.9809d0,&
                                                   1030.1085d0,1486.058d0,2043.8428d0,2715.89d0,3511.696d0,4420.0d0,5451.06d0,6604.95d0,&
                                                   7888.53d0,9305.8d0,10859.7d0,12556.4d0,14400.8d0,16382.0d0,18510.0d0,20788.0d0,&
                                                   23221.0d0,25820.0d0,28582.0d0,31514.0d0,34619.0d0,37899.0d0,41356.0d0]

         integer :: k, Z(species)
         real(dp) :: A(species)

         ! Get basic species info
         do k=1,species
            A(k) = chem_isos% Z_plus_N(chem_id(k))  ! baryon number
            Z(k) = chem_isos% Z(chem_id(k))         ! charge
         end do

         ! Compute offset in eV/baryon
         offset = 0d0
         do k=1,species
            if (Z(k) <= 28 .and. Z(k) >= 1) then
               ! offset = \sum_k E_k * X_k / A_k, here in units of eV/amu
               offset = offset + ionization_table(Z(k)) * xa(k) / A(k)
            end if
         end do

         ! Translate eV/amu into erg/g
         offset = offset * ev2erg / amu

      end function compute_ion_offset

end module ion_offset
