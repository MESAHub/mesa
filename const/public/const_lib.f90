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

module const_lib

   use const_def

   implicit none

   private
   public :: const_init

contains

   subroutine const_init(mesa_dir_init, ierr)
      character(len=*), intent(in) :: mesa_dir_init
      integer, intent(out) :: ierr
      ierr = 0
      call do_const_init(mesa_dir_init, ierr)
   end subroutine const_init

end module const_lib
