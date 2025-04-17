! ***********************************************************************
!
!   Copyright (C) 2019  The MESA Team
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

module other_timestep_limit

   ! consult star/other/README for general usage instructions
   ! control name: use_other_timestep_limit = .true.
   ! procedure pointer: s% other_timestep_limit => my_routine

   implicit none

contains

   integer function null_other_timestep_limit( &
      id, skip_hard_limit, dt, dt_limit_ratio)
      use const_def, only: dp
      use star_def
      integer, intent(in) :: id
      logical, intent(in) :: skip_hard_limit
      real(dp), intent(in) :: dt
      real(dp), intent(inout) :: dt_limit_ratio
      null_other_timestep_limit = keep_going
   end function null_other_timestep_limit

end module other_timestep_limit

