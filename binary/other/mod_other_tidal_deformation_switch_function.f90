! ***********************************************************************
!
!   Copyright (C) 2025  Matthias Fabry and the MESA Team
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

module mod_other_tidal_deformation_switch_function

   ! NOTE: remember to set true:
   ! use_other_tidal_deformation_function = .true.

   ! you can add your own routine for use instead of the default one.
   ! the default routine uses the synchronicity parameter to switch from single rotating to tidal corrections.
   ! When the shell is quite synchronous, f_switch -> 1, when not synchronous, f_switch -> 0. It uses a sigmoid around
   ! omega/omega_sync = b% f_sync_switch_from_rot_defor to smoothly vary f_switch from 0 to 1.

   use const_def

   implicit none

   contains

   subroutine null_other_tidal_deformation_switch_function(id, k, omega_in, f_switch, ierr)
      integer, intent(in) :: id, k
      real(dp), intent(in) :: omega_in
      real(dp), intent(out) :: f_switch
      integer, intent(out) :: ierr
      include 'formats'

      write (*, 1) 'no implementation for other_tidal_deformation_function'
      ! must set f_switch
      ierr = -1
   end subroutine null_other_tidal_deformation_switch_function


end module mod_other_tidal_deformation_switch_function

