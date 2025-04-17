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

module other_torque_implicit

   ! consult star/other/README for general usage instructions
   ! control name: use_other_torque_implicit = .true.
   ! procedure pointer: s% other_torque_implicit => my_routine

   implicit none

contains

   subroutine default_other_torque_implicit(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      s%extra_jdot(:) = 0
      s%extra_omegadot(:) = 0
      s%d_extra_jdot_domega_m1(:) = 0
      s%d_extra_omegadot_domega_m1(:) = 0
      s%d_extra_jdot_domega_00(:) = 0
      s%d_extra_omegadot_domega_00(:) = 0
      s%d_extra_jdot_domega_p1(:) = 0
      s%d_extra_omegadot_domega_p1(:) = 0
   end subroutine default_other_torque_implicit

end module other_torque_implicit

