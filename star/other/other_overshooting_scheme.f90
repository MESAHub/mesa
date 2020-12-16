! ***********************************************************************
!
!   Copyright (C) 2019  Bill Paxton and MESA Team
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
 
      module other_overshooting_scheme

      ! consult star/other/README for general usage instructions
      ! procedure pointer: s% other_overshooting_scheme => my_routine
      ! note that this is enabled by setting s%overshooting_scheme = 'other'
      ! so there is no use_other_overshooting_scheme flag.


      use star_def

      implicit none
      
      contains
            
      subroutine null_other_overshooting_scheme(id, i, j, k_a, k_b, D, vc, ierr)
         integer, intent(in) :: id, i, j
         integer, intent(out) :: k_a, k_b
         real(dp), intent(out), dimension(:) :: D, vc
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         k_a = -1
         k_b = -1
         D = 0d0
         vc = 0d0

         ierr = -1
      end subroutine null_other_overshooting_scheme

      end module other_overshooting_scheme
      
      
      
      
