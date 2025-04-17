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

module other_am_mixing

   ! consult star/other/README for general usage instructions
   ! control name: use_other_am_mixing = .true.
   ! procedure pointer: s% other_am_mixing => my_routine

   implicit none

contains

   ! for extra angular momentum mixing
   ! you can revise vectors s% am_nu_omega(:) and s% am_nu_j(:)

   subroutine null_other_am_mixing(id, ierr)
      use star_def
      use const_def, only: dp
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_other_am_mixing

end module other_am_mixing

