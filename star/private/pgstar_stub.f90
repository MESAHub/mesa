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

module pgstar
   
   use star_def
   use const_def
   use star_pgstar
   
   implicit none


contains
   
   
   subroutine do_create_file_name(s, dir, prefix, name)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: dir, prefix
      character (len = *), intent(out) :: name
      name = ''
   end subroutine do_create_file_name
   
   
   subroutine do_write_plot_to_file(s, p, filename, ierr)
      type (star_info), pointer :: s
      type (pgstar_win_file_data), pointer :: p
      character (len = *), intent(in) :: filename
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine do_write_plot_to_file
   
   
   subroutine do_show_pgstar_annotations(&
      s, show_annotation1, show_annotation2, show_annotation3)
      type (star_info), pointer :: s
      logical, intent(in) :: &
         show_annotation1, show_annotation2, show_annotation3
   end subroutine do_show_pgstar_annotations
   
   
   subroutine do_start_new_run_for_pgstar(s, ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine do_start_new_run_for_pgstar
   
   
   subroutine do_restart_run_for_pgstar(s, ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine do_restart_run_for_pgstar
   
   
   subroutine do_read_pgstar_controls(s, inlist_fname, ierr)
      type (star_info), pointer :: s
      character(*), intent(in) :: inlist_fname
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine do_read_pgstar_controls
   
   
   subroutine do_pgstar_plots(&
      s, must_write_files, &
      ierr)
      type (star_info), pointer :: s
      logical, intent(in) :: must_write_files
      
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine do_pgstar_plots
   
   
   subroutine do_set_xaxis_bounds(&
      s, xaxis_by, win_xmin_in, win_xmax_in, xmargin, &
      xvec, xmin, xmax, xleft, xright, dx, &
      grid_min, grid_max, npts, ierr)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: xaxis_by
      real, intent(in) :: win_xmin_in, win_xmax_in, xmargin
      real, allocatable, dimension(:) :: xvec
      real, intent(out) :: xmin, xmax, xleft, xright, dx
      integer, intent(out) :: grid_min, grid_max, npts
      integer, intent(out) :: ierr
      xmin = 0; xmax = 0; xleft = 0; xright = 0; dx = 0
      grid_min = 0; grid_max = 0; npts = 0; ierr = 0
   end subroutine do_set_xaxis_bounds
   
   
   subroutine do_show_xaxis_by(s, by, ierr)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: by
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine do_show_xaxis_by
   
   
   subroutine show_box_pgstar(s, str1, str2)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: str1, str2
   end subroutine show_box_pgstar
   
   
   subroutine show_title_pgstar(s, title, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: title
      real, intent(in) :: pad
      optional pad
   end subroutine show_title_pgstar
   
   
   subroutine show_xaxis_label_pgstar(s, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
   end subroutine show_xaxis_label_pgstar
   
   
   subroutine show_left_yaxis_label_pgstar(s, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
   end subroutine show_left_yaxis_label_pgstar
   
   
   subroutine show_right_yaxis_label_pgstar(s, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
   end subroutine show_right_yaxis_label_pgstar
   
   
   subroutine show_left_yaxis_label_pgmtxt_pgstar(&
      s, coord, fjust, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
   end subroutine show_left_yaxis_label_pgmtxt_pgstar
   
   
   subroutine show_right_yaxis_label_pgmtxt_pgstar(&
      s, coord, fjust, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
   end subroutine show_right_yaxis_label_pgmtxt_pgstar
   
   
   subroutine show_model_number_pgstar(s)
      type (star_info), pointer :: s
   end subroutine show_model_number_pgstar
   
   
   subroutine show_age_pgstar(s)
      type (star_info), pointer :: s
   end subroutine show_age_pgstar
   
   subroutine update_pgstar_data(s, ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine update_pgstar_data
   
   subroutine shutdown_pgstar(s)
      type (star_info), pointer :: s
   end subroutine shutdown_pgstar

end module pgstar
      
