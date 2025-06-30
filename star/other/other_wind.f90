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

module other_wind

   ! consult star/other/README for general usage instructions
   ! control name: use_other_wind = .true.
   ! procedure pointer: s% other_wind => my_routine

   ! you can add your own wind routine for use when wind scheme == 'other'

   implicit none

contains

   ! Note that your routine will be called before many star variables have been set.
   ! If you rely on these, you should call the star_set_vars_in_part1 routine from star_lib
   ! to ensure that they are set.
   subroutine null_other_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
      use star_def
      integer, intent(in) :: id
      real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z  ! surface values (cgs)
      ! NOTE: surface is outermost cell. not necessarily at photosphere.
      ! NOTE: don't assume that vars are set at this point.
      ! so if you want values other than those given as args,
      ! you should use values from s% xh(:,:) and s% xa(:,:) only.
      ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
      real(dp), intent(out) :: w  ! wind in units of Msun/year (value is >= 0)
      integer, intent(out) :: ierr
      w = 0
      ierr = 0
   end subroutine null_other_wind

end module other_wind

