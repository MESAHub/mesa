! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   https://mesastar.org/
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

module other_before_struct_burn_mix

   ! consult star/other/README for general usage instructions
   ! control name: use_other_before_struct_burn_mix = .true.
   ! procedure pointer: s% other_before_struct_burn_mix => my_routine

   implicit none

contains

   ! your routine will be called before the standard struct_burn_mix routine

   subroutine null_other_before_struct_burn_mix(id, dt, res)
      use const_def, only: dp
      use star_def
      integer, intent(in) :: id
      real(dp), intent(in) :: dt
      integer, intent(out) :: res  ! keep_going, redo, retry, terminate
      res = keep_going
   end subroutine null_other_before_struct_burn_mix

end module other_before_struct_burn_mix

