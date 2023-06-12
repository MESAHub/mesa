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
 
      module other_am_mixing

      ! consult star/other/README for general usage instructions
      ! control name: use_other_am_mixing = .true.
      ! procedure pointer: s% other_am_mixing => my_routine

      implicit none
      
      contains

      ! for extra angular momentum mixing
      ! you can revise vectors s% am_nu_omega(:) and s% am_nu_j(:)
      
      subroutine null_other_am_mixing(id, ierr)
         use star_def
         use const_def   
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_other_am_mixing
      
      
      end module other_am_mixing
      
      
      
      
