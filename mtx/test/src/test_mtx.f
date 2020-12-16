! ***********************************************************************
!
!   Copyright (C) 2011-2018  Bill Paxton & The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not,write to the Free Software
!   Foundation,Inc.,59 Temple Place,Suite 330,Boston,MA 02111-1307 USA
!
! ***********************************************************************

      program test_mtx
      use const_lib
      use math_lib, only: math_init
      use mtx_lib
      use test_mtx_support

      use test_square
      use test_square_quad
      
      use test_block_tri_dble, only: do_test_block_tri_dble
      use test_block_tri_quad, only: do_test_block_tri_quad
      
      use utils_lib, only: mesa_error

      implicit none
      
      character (len=32) :: my_mesa_dir
      integer :: ierr, i

      my_mesa_dir = '../..'         
      call const_init(my_mesa_dir,ierr)     
      if (ierr /= 0) then
         write(*,*) 'const_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if        
      
      call math_init()
      
      call do_test_square
      call do_test_square_quad
      call do_test_block_tri_dble
      call do_test_block_tri_quad      
      call test_format_conversion      


      end program

