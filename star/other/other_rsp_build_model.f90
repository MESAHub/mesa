! ***********************************************************************
!
!   Copyright (C) 2018  The MESA Team
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

module other_rsp_build_model
   
   ! consult star/other/README for general usage instructions
   ! control name: use_other_RSP_build_model = .true.
   ! procedure pointer: s% other_rsp_build_model => my_routine
   
   use star_def
   
   implicit none


contains
   
   ! here is a list of what this routine needs to set.
   ! star arrays have been allocated but nothing else.
   ! note that nz is set by RSP_nz.
   
   ! scalars
   ! M_center, L_center, R_center, v_center
   ! star_mass, mstar, xmstar, tau_factor, rsp_period
   
   ! vectors
   ! m, dm, dm_bar, r, Vol, v, T, w, Fr, erad, L
   ! don't need to set xh or xa
   
   subroutine null_other_rsp_build_model(id, ierr)
      integer, intent(in) :: id ! star id if available; 0 otherwise
      integer, intent(out) :: ierr ! 0 means AOK.
      write(*, *) 'no implementation for other_rsp_build_model'
      ierr = -1
   end subroutine null_other_rsp_build_model


end module other_rsp_build_model
      
      
      
      
