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

      module interp_1d_def

      implicit none


      integer, parameter :: mp_work_size = 14
      integer, parameter :: pm_work_size = 3

      ! these are the limiter options for the monotonicity preserving interpolation

      integer, parameter :: average = 1
      integer, parameter :: quartic = 2
      integer, parameter :: super_bee = 3

      end module interp_1d_def
