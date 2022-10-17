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

module pgstar_astero_plots
   use star_lib
   use star_def
   use star_pgstar
   
   implicit none


contains
   
   
   subroutine astero_pgstar_plots_info(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine astero_pgstar_plots_info
   
   
   subroutine write_plot_to_file(s, p, file_prefix, number, ierr)
      type (star_info), pointer :: s
      type (pgstar_win_file_data), pointer :: p
      character (len = *), intent(in) :: file_prefix
      integer, intent(in) :: number
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine write_plot_to_file


end module pgstar_astero_plots
      
      
      
      
