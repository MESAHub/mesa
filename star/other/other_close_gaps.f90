! ***********************************************************************
!
!   Copyright (C) 2023  The MESA Team & Simon Guichandut
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

module other_close_gaps

   ! consult star/other/README for general usage instructions
   ! control name: use_other_close_gaps = .true.
   ! procedure pointer: s% other_close_gaps => my_routine
   ! This also requires the control remove_mixing_glitches = .true.

   implicit none

contains

   subroutine null_other_close_gaps(id, mix_type, min_gap, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(in) :: mix_type
      real(dp), intent(in) :: min_gap
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_other_close_gaps

end module other_close_gaps

