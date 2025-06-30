! ***********************************************************************
!
!   Copyright (C) 2018  The MESA Team
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

module other_rsp_linear_analysis

   ! consult star/other/README for general usage instructions
   ! control name: use_other_RSP_linear_analysis = .true.
   ! procedure pointer: s% other_rsp_linear_analysis => my_routine

   use star_def

   implicit none

contains

   subroutine null_other_rsp_linear_analysis(id, restart, ierr)
      use star_def
      integer, intent(in) :: id  ! star id if available; 0 otherwise
      logical, intent(in) :: restart
      integer, intent(out) :: ierr  ! 0 means AOK.

      write (*, *) 'no implementation for other_rsp_linear_analysis'
      ierr = -1

   end subroutine null_other_rsp_linear_analysis

end module other_rsp_linear_analysis

