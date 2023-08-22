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
 
      module other_momentum_implicit

      ! consult star/other/README for general usage instructions
      ! control name: use_other_momentum_implicit = .true.
      ! procedure pointer: s% other_momentum_implicit => my_routine


      implicit none
      
            
      contains
      
      
      subroutine default_other_momentum_implicit(id, ierr)
         use star_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         ! s% extra_grav(1:s%nz) = 0 ! this is added to -G*m/r^2
         ! note that extra_grav is type(auto_diff_real_star_order1) so includes partials.
      end subroutine default_other_momentum_implicit


      end module other_momentum_implicit
      
      
      
      
