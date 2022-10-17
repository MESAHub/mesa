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

module other_cgrav
   
   ! consult star/other/README for general usage instructions
   ! control name: use_other_cgrav = .true.
   ! procedure pointer: s% other_cgrav => my_routine
   
   ! note that other_cgrav only changes G in the stellar structure
   ! the binary module is unaffected by changes in cgrav
   
   use star_def
   
   implicit none


contains
   
   
   subroutine default_other_cgrav(id, ierr)
      use const_def, only : standard_cgrav
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      s% cgrav(:) = standard_cgrav
   end subroutine default_other_cgrav

end module other_cgrav
      
      
      
      
