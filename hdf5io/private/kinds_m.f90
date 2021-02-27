! Module  : kinds_m
! Purpose : kind type parameter definitions
!
! Copyright 2021 Rich Townsend
!
! This file is part of the ForUM (Fortran Utility Modules)
! package. ForUM is free software: you can redistribute it and/or
! modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, version 3.
!
! ForUM is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module kinds_m

  ! Uses

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: IS = INT32
  integer, parameter :: ID = INT64

  integer, parameter :: RS = REAL32
  integer, parameter :: RD = REAL64

  ! Access specifiers

  private

  public :: IS
  public :: ID
  public :: RS
  public :: RD

end module kinds_m
