! ***********************************************************************
!
!   Copyright (C) 2016  The MESA Team
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

module other_gradr_factor

   ! consult star/other/README for general usage instructions
   ! control name: use_other_gradr_factor = .true.
   ! procedure pointer: s% other_gradr_factor => my_routine

   implicit none

contains

   subroutine null_other_gradr_factor(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      do k = 1, s%nz
         s%gradr_factor(k) = 1d0
      end do
   end subroutine null_other_gradr_factor

end module other_gradr_factor

