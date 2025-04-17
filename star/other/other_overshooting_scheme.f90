! ***********************************************************************
!
!   Copyright (C) 2019  The MESA Team
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

module other_overshooting_scheme

   ! consult star/other/README for general usage instructions
   ! procedure pointer: s% other_overshooting_scheme => my_routine
   ! note that this is enabled by setting s%overshooting_scheme = 'other'
   ! so there is no use_other_overshooting_scheme flag.

   implicit none

contains

   subroutine null_other_overshooting_scheme(id, i, j, k_a, k_b, D, vc, ierr)
      use star_def
      integer, intent(in) :: id, i, j
      integer, intent(out) :: k_a, k_b
      real(dp), intent(out), dimension(:) :: D, vc
      integer, intent(out) :: ierr
      k_a = -1
      k_b = -1
      D = 0d0
      vc = 0d0

      ierr = -1
   end subroutine null_other_overshooting_scheme

end module other_overshooting_scheme
