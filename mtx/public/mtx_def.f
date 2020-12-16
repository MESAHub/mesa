! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module mtx_def
      
      implicit none
      
      ! matrix solver options
      integer, parameter :: lapack = 1
      integer, parameter :: block_thomas_dble = 2
      integer, parameter :: block_thomas_quad = 3
      integer, parameter :: block_thomas_refine = 4
      integer, parameter :: bcyclic_dble = 5

      ! sparse matrix formats
      integer, parameter :: compressed_column_sparse = 0
      integer, parameter :: compressed_row_sparse = 1
      
      ! 0 based formats for internal use only
      integer, parameter :: compressed_col_sparse_0_based = 2
      integer, parameter :: compressed_row_sparse_0_based = 3

      end module mtx_def

