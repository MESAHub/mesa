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

module other_torque

   ! consult star/other/README for general usage instructions
   ! control name: use_other_torque = .true.
   ! procedure pointer: s% other_torque => my_routine

   implicit none

contains

   subroutine default_other_torque(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      ! note that can set extra_omegadot instead of extra_jdot if that is more convenient
      ! set one or the other, not both.  set the one you are not using to 0 as in the following line.
      s%extra_jdot(:) = 0
      s%extra_omegadot(:) = 0
   end subroutine default_other_torque

end module other_torque

