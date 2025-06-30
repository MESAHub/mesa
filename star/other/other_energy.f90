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

module other_energy

   ! consult star/other/README for general usage instructions
   ! control name: use_other_energy = .true.
   ! procedure pointer: s% other_energy => my_routine

   implicit none

contains

   subroutine default_other_energy(id, ierr)
      use star_def
      use auto_diff
      use const_def, only: Rsun
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      !s% extra_heat(:) = s% extra_power_source
      ! note that extra_heat is type(auto_diff_real_star_order1) so includes partials.
   end subroutine default_other_energy

end module other_energy

