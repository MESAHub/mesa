! ***********************************************************************
!
!   Copyright (C) 2010-2022  The MESA Team
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

module math_def

   use const_def, only: dp

   implicit none

   integer, parameter :: max_precomp_ints = 1000

   type precomp_int
      real(dp) :: z2, z3, z4, z5, z6, z7, z8

      real(dp) :: z1_3, z2_3, z4_3, z5_3, z7_3
      real(dp) :: zm1_3, zm2_3, zm4_3, zm5_3, zm7_3

      real(dp) :: z7_6

      real(dp) :: z1_5, z2_5, z3_5, z4_5

      real(dp) :: z1_2, z3_2, zm1_2, zm1_4, zm3_2

      real(dp) :: logz, sqlogz  ! pow2(log(z))
      real(dp) :: logz_3_2  ! pow(log(Z), 1.5d0)

      real(dp) :: z0p475, zp1_3_2, zm0p267

   end type precomp_int

   type(precomp_int), dimension(max_precomp_ints) :: pre_z  ! Set in the math_lib

end module math_def
