! ***********************************************************************
!
!   Copyright (C) 2011-2018  Bill Paxton & The MESA Team
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

program test_mtx
   use const_lib, only: const_init
   use math_lib, only: math_init
   use mtx_lib
   use test_mtx_support
   use test_square
   use test_block_tri_dble, only: do_test_block_tri_dble
   use utils_lib, only: mesa_error

   implicit none

   character(len=32) :: my_mesa_dir
   integer :: ierr

   my_mesa_dir = '../..'
   call const_init(my_mesa_dir, ierr)
   if (ierr /= 0) then
      write (*, *) 'const_init failed'
      call mesa_error(__FILE__, __LINE__)
   end if

   call math_init()

   call do_test_square
   call do_test_block_tri_dble
   call test_format_conversion

end program test_mtx
