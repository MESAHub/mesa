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

module other_rsp_linear_analysis

   ! consult star/other/README for general usage instructions
   ! control name: use_other_RSP_linear_analysis = .true.
   ! procedure pointer: s% other_rsp_linear_analysis => my_routine

   use star_def

   implicit none

contains

   subroutine null_other_rsp_linear_analysis(id, restart, ierr)
      use star_def
      integer, intent(in) :: id  ! star id if available; 0 otherwise
      logical, intent(in) :: restart
      integer, intent(out) :: ierr  ! 0 means AOK.

      write (*, *) 'no implementation for other_rsp_linear_analysis'
      ierr = -1

   end subroutine null_other_rsp_linear_analysis

end module other_rsp_linear_analysis

