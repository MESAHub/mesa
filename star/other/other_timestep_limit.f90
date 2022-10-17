! ***********************************************************************
!
!   Copyright (C) 2019  The MESA Team
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

module other_timestep_limit
   
   ! consult star/other/README for general usage instructions
   ! control name: use_other_timestep_limit = .true.
   ! procedure pointer: s% other_timestep_limit => my_routine
   
   use star_def
   
   implicit none

contains
   
   integer function null_other_timestep_limit(&
      id, skip_hard_limit, dt, dt_limit_ratio)
      use const_def, only : dp
      integer, intent(in) :: id
      logical, intent(in) :: skip_hard_limit
      real(dp), intent(in) :: dt
      real(dp), intent(inout) :: dt_limit_ratio
      null_other_timestep_limit = keep_going
   end function null_other_timestep_limit

end module other_timestep_limit
      
      
      
      
