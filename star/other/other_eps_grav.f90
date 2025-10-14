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

module other_eps_grav

   ! consult star/other/README for general usage instructions
   ! control name: use_other_eps_grav = .true.
   ! procedure pointer: s% other_eps_grav => my_routine

   implicit none

contains

   subroutine null_other_eps_grav(id, k, dt, ierr)
      use star_def
      integer, intent(in) :: id, k
      real(dp), intent(in) :: dt
      integer, intent(out) :: ierr
      ierr = 0

      ! NOTE: this is called after 1st doing the standard eps_grav calculation.
      ! so if you decide you don't want to do anything special, just return.

      ! NOTE: need to set the auto_diff variable s% eps_grav_ad(k)
      ! this is type(auto_diff_real_star_order1) so includes partials.
      ! in addition to setting the value,
      ! you must also set the partials, and there are lots of them.

   end subroutine null_other_eps_grav

end module other_eps_grav

