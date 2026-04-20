! ***********************************************************************
!
!   Copyright (C) 2024  The MESA Team
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

program run_timing

   use const_lib, only: const_init
   use math_lib, only: math_init
   use time_star_bcyclic, only: time_solvers

   implicit none

   character(len=32) :: my_mesa_dir
   integer :: ierr

   my_mesa_dir = '../..'
   call const_init(my_mesa_dir, ierr)
   if (ierr /= 0) then
      write(*, *) 'const_init failed'
      stop 'const_init failed'
   end if

   call math_init()

   call time_solvers()

end program run_timing
