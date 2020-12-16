! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
 
      module other_diffusion

      ! consult star/other/README for general usage instructions
      ! control name: use_other_diffusion = .true.
      ! procedure pointer: s% other_diffusion => my_routine



      use star_def

      implicit none
      
            
      contains
      
      
      subroutine null_other_diffusion(id, dt, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: dt 
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         !write(*,*) 'null_other_diffusion'
         !write(*,*) 'associated(s% edv)', associated(s% edv)
         s% edv(:,1:s% nz) = 0
         !write(*,*) 'done null_other_diffusion'
      end subroutine null_other_diffusion


      end module other_diffusion
      
      
      
      
