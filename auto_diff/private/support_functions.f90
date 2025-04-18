! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
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

module support_functions

   use const_def, only: dp
   use math_lib

   implicit none

   interface pow
      module procedure int_real_pow
      module procedure int_int_pow
   end interface pow

   interface log
      module procedure log_int
   end interface log

   interface max
      module procedure max_int_real
      module procedure max_real_int
   end interface max

   interface min
      module procedure min_int_real
      module procedure min_real_int
   end interface min

   contains

   pure real(dp) function sgn(x) result(res)
      real(dp), intent(in) :: x
      if (x < 0d0) then
         res = -1d0
      else if (x == 0d0) then
         res = 0d0
      else
         res = 1d0
      end if
   end function sgn

   pure real(dp) function Heaviside(x) result(res)
      real(dp), intent(in) :: x
      if (x < 0d0) then
         res = 0d0
      else if (x == 0d0) then
         res = 0.5d0
      else
         res = 1d0
      end if
   end function Heaviside

   pure real(dp) function int_real_pow(x, y) result(z)
      integer, intent(in) :: x
      real(dp), intent(in) :: y
      real(dp) :: x_real

      x_real = x
      z = pow(x_real,y)
   end function int_real_pow

   pure real(dp) function int_int_pow(x, y) result(z)
      integer, intent(in) :: x, y

      z = x**y
   end function int_int_pow

   pure real(dp) function log_int(x) result(res)
      integer, intent(in) :: x
      real(dp) :: x_real

      x_real = x
      res = log(x_real)

   end function log_int

   pure real(dp) function max_real_int(x, y) result(z)
      real(dp), intent(in) :: x
      integer, intent(in) :: y
      if (x > y) then
         z = x
      else
         z = y
      end if
   end function max_real_int

   pure real(dp) function max_int_real(x, y) result(z)
      integer, intent(in) :: x
      real(dp), intent(in) :: y
      if (x > y) then
         z = x
      else
         z = y
      end if
   end function max_int_real

   pure real(dp) function min_int_real(x, y) result(z)
      integer, intent(in) :: x
      real(dp), intent(in) :: y
      if (x < y) then
         z = x
      else
         z = y
      end if
   end function min_int_real


   pure real(dp) function min_real_int(x, y) result(z)
      integer, intent(in) :: y
      real(dp), intent(in) :: x
      if (x < y) then
         z = x
      else
         z = y
      end if
   end function min_real_int

end module support_functions