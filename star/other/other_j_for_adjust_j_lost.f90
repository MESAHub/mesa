! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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

module other_j_for_adjust_J_lost

   ! consult star/other/README for general usage instructions
   ! control name: use_other_j_for_adjust_J_lost = .true.
   ! procedure pointer: s% other_j_for_adjust_J_lost => my_routine

   implicit none

contains

   ! your routine will be called after winds and before mass adjustment

   subroutine null_other_j_for_adjust_J_lost(id, starting_j_rot_surf, j_for_mass_loss, ierr)
      use star_def
      integer, intent(in) :: id
      real(dp), intent(in) :: starting_j_rot_surf
      real(dp), intent(out) :: j_for_mass_loss
      integer, intent(out) :: ierr
      write (*, *) 'no implementation for other_j_for_adjust_J_lost'
      ierr = -1
      j_for_mass_loss = -1d99
   end subroutine null_other_j_for_adjust_J_lost

end module other_j_for_adjust_J_lost

