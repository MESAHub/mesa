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

module other_brunt_smoothing

   ! consult star/other/README for general usage instructions
   ! control name: use_other_brunt_smoothing = .true.
   ! procedure pointer: s% other_brunt_smoothing => my_routine

   implicit none

contains

   subroutine null_other_brunt_smoothing(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      !type(star_info), pointer :: s
      !integer :: k
      ierr = 0
   end subroutine null_other_brunt_smoothing

end module other_brunt_smoothing

