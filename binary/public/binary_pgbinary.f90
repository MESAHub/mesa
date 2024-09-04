! ***********************************************************************
!
!   Copyright (C) 2024  Matthias Fabry and the MESA Team
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
!
! ***********************************************************************

module binary_pgbinary

   use const_def
   use star_pgstar

   abstract interface

      subroutine pgbinary_plot_interface(id, device_id, array_ix, ierr)
         integer, intent(in) :: id, device_id, array_ix
         integer, intent(out) :: ierr
      end subroutine pgbinary_plot_interface

      subroutine other_do_plot_in_binary_grid_interface(&
         id, device_id, xleft, xright, ybot, ytop, txt_scale, ierr)
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         integer, intent(out) :: ierr
      end subroutine other_do_plot_in_binary_grid_interface

      subroutine pgbinary_decorator_interface(id, xmin, xmax, ymin, ymax, plot_num, ierr)
         integer, intent(in) :: id
         !Not doubles
         real, intent(in) :: xmin, xmax, ymin, ymax
         integer, intent(in) :: plot_num
         integer, intent(out) :: ierr
      end subroutine pgbinary_decorator_interface

   end interface

   type pgbinary_win_file_data
      integer :: id, array_ix
      character (len = strlen) :: name
      logical :: win_flag, file_flag, do_win, do_file
      integer :: id_win, id_file, file_interval
      real :: win_width, prev_win_width
      real :: win_aspect_ratio, prev_win_aspect_ratio
      real :: file_width, file_aspect_ratio
      character (len = strlen) :: file_dir, file_prefix, most_recent_filename
      character (len = strlen) :: file_dir_for_previous_mkdir
      logical :: have_called_mkdir
      procedure (pgbinary_plot_interface), pointer, nopass :: plot => null()
      ! the following make it possible to use "other" plots in grids
      logical :: okay_to_call_do_plot_in_binary_grid
      procedure (other_do_plot_in_binary_grid_interface), pointer, nopass :: &
         do_plot_in_binary_grid => null()
   end type pgbinary_win_file_data

   integer, parameter :: max_num_Binary_History_Panels = 16
   integer, parameter :: max_num_Summary_Binary_History_Lines = 16
   integer, parameter :: max_num_Binary_Other_plots = 16

   integer, parameter :: pgbinary_array_length = 9

   integer, parameter :: i_Summary_Binary_History = 1
   integer, parameter :: i_Binary_Text_Summary = i_Summary_Binary_History + 1
   integer, parameter :: i_Binary_Hist_Track = i_Binary_Text_Summary + pgbinary_array_length
   integer, parameter :: i_Binary_Hist_Panels = i_Binary_Hist_Track + pgbinary_array_length
   integer, parameter :: i_Binary_Grid = i_Binary_Hist_Panels + pgbinary_array_length
   integer, parameter :: i_Binary_Star1 = i_Binary_Grid + pgbinary_array_length
   integer, parameter :: i_Binary_Star2 = i_Binary_Star1 + 1
   integer, parameter :: i_Binary_Other = i_Binary_Star2 + 1

   integer, parameter :: num_pgbinary_plots = i_Binary_Other + max_num_Binary_Other_plots

   integer, parameter :: max_num_pgbinary_Grid_plots = 10

   integer, parameter :: max_num_pgbinary_trace_Binary_History_values = 20

   type pgbinary_hist_node
      real(dp) :: age
      integer :: step
      real(dp), pointer :: vals(:) => null() ! values of items in history_columns list
      type (pgbinary_hist_node), pointer :: next => null()
      ! list kept in strictly decreasing order by age & step
   end type pgbinary_hist_node

   type pgbinary_controls
      include "pgbinary_controls.inc"

      type (pgbinary_hist_node), pointer :: pgbinary_hist => null()

      procedure(pgbinary_decorator_interface), pointer, nopass :: &
         History_Panels1_pgbinary_decorator => null(), &
         History_Panels2_pgbinary_decorator => null(), &
         History_Panels3_pgbinary_decorator => null(), &
         History_Panels4_pgbinary_decorator => null(), &
         History_Panels5_pgbinary_decorator => null(), &
         History_Panels6_pgbinary_decorator => null(), &
         History_Panels7_pgbinary_decorator => null(), &
         History_Panels8_pgbinary_decorator => null(), &
         History_Panels9_pgbinary_decorator => null(), &
         History_Track1_pgbinary_decorator => null(), &
         History_Track2_pgbinary_decorator => null(), &
         History_Track3_pgbinary_decorator => null(), &
         History_Track4_pgbinary_decorator => null(), &
         History_Track5_pgbinary_decorator => null(), &
         History_Track6_pgbinary_decorator => null(), &
         History_Track7_pgbinary_decorator => null(), &
         History_Track8_pgbinary_decorator => null(), &
         History_Track9_pgbinary_decorator => null(), &
         summary_history_pgbinary_decorator => null()

   end type pgbinary_controls

end module binary_pgbinary
