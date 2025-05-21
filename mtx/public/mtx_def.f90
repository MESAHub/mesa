! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton & The MESA Team
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

module mtx_def

   implicit none

   ! matrix solver options
   integer, parameter :: lapack = 1
   integer, parameter :: block_thomas_dble = 2
   integer, parameter :: block_thomas_refine = 4
   integer, parameter :: bcyclic_dble = 5

   ! sparse matrix formats
   integer, parameter :: compressed_column_sparse = 0
   integer, parameter :: compressed_row_sparse = 1

   ! 0 based formats for internal use only
   integer, parameter :: compressed_col_sparse_0_based = 2
   integer, parameter :: compressed_row_sparse_0_based = 3

end module mtx_def

