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
 
      module other_astero_freq_corr

      ! consult star/other/README for general usage instructions
      ! control name: use_other_astero_freq_corr = .true.
      ! procedure pointer: s% other_astero_freq_corr => my_routine


      use star_def

      implicit none
      
            
      contains
      
      
      subroutine default_other_astero_freq_corr(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         ! see star/astero/src/astero_support.f, subroutine get_freq_corr.
         ! note that your routine can use both astero_data and astero_support.
      end subroutine default_other_astero_freq_corr


      end module other_astero_freq_corr
      
      
      
      
