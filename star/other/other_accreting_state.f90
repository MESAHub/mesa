! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

module other_accreting_state

   ! consult star/other/README for general usage instructions
   ! control name: use_other_accreting_state = .true.
   ! procedure pointer: s% other_accreting_state => my_routine

   implicit none

contains

   ! Note that your routine will be called before many star variables have been set.
   ! If you rely on these, you should call the star_set_vars_in_part1 routine from star_lib
   ! to ensure that they are set.
   subroutine null_other_accreting_state(id, total_specific_energy, accretion_pressure, accretion_density, ierr)
      use star_def
      integer, intent(in) :: id
      real(dp), intent(out) :: total_specific_energy, accretion_pressure, accretion_density
      integer, intent(out) :: ierr
      total_specific_energy = 0d0  ! erg/g
      accretion_pressure = 0d0  ! erg/cm^3
      accretion_density = 0d0  ! g/cm^3
      ierr = 0
   end subroutine null_other_accreting_state

end module other_accreting_state

