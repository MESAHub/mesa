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

module star_data_lib

   implicit none

contains

   subroutine star_data_init(mesa_dir_init, ierr)
      use star_data_def, only: do_star_def_init
      character(len=*), intent(in) :: mesa_dir_init
      integer, intent(out) :: ierr
      ierr = 0
      call do_star_def_init(mesa_dir_init, ierr)
   end subroutine star_data_init

end module star_data_lib
