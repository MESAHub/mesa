! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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
 
      module other_eps_grav

      ! consult star/other/README for general usage instructions
      ! control name: use_other_eps_grav = .true.
      ! procedure pointer: s% other_eps_grav => my_routine


      use star_def

      implicit none
      
            
      contains
      
      
      subroutine null_other_eps_grav(id, k, dt, ierr)
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
      
      
      
      
