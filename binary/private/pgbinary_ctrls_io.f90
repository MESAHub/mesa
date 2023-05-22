! ***********************************************************************
!
!   Copyright (C) 2010-2022  The MESA Team, Bill Paxton & Matthias Fabry
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

module pgbinary_ctrls_io

   use const_def
   use binary_private_def

   implicit none

   include "pgbinary_controls.inc"

   namelist /pgbinary/ &

      file_device, &
      file_extension, &
      file_digits, &
      pgbinary_interval, &
      pause, &
      pause_interval, &
      pgbinary_sleep, &
      clear_history, &

      file_white_on_black_flag, &
      win_white_on_black_flag, &
      pgbinary_show_model_number, &
      pgbinary_show_age, &
      pgbinary_show_age_in_seconds, &
      pgbinary_show_age_in_minutes, &
      pgbinary_show_age_in_hours, &
      pgbinary_show_age_in_days, &
      pgbinary_show_age_in_years, &
      pgbinary_show_log_age_in_years, &

      pgbinary_report_writing_files, &

      pgbinary_show_title, &
      pgbinary_title_scale, &
      pgbinary_title_disp, &
      pgbinary_title_coord, &
      pgbinary_title_fjust, &
      pgbinary_title_lw, &

      pgbinary_grid_show_title, &
      pgbinary_grid_title_scale, &
      pgbinary_grid_title_disp, &
      pgbinary_grid_title_coord, &
      pgbinary_grid_title_fjust, &
      pgbinary_grid_title_lw, &

      pgbinary_age_scale, &
      pgbinary_age_disp, &
      pgbinary_age_coord, &
      pgbinary_age_fjust, &
      pgbinary_age_lw, &

      pgbinary_model_scale, &
      pgbinary_model_disp, &
      pgbinary_model_coord, &
      pgbinary_model_fjust, &
      pgbinary_model_lw, &

      pgbinary_xaxis_label_scale, &
      pgbinary_left_yaxis_label_scale, &
      pgbinary_right_yaxis_label_scale, &
      pgbinary_xaxis_label_lw, &
      pgbinary_left_yaxis_label_lw, &
      pgbinary_right_yaxis_label_lw, &
      pgbinary_xaxis_label_disp, &
      pgbinary_left_yaxis_label_disp, &
      pgbinary_right_yaxis_label_disp, &
      pgbinary_num_scale, &
      pgbinary_lw, &
      pgbinary_model_lw, &
      pgbinary_box_lw, &

      Text_Summary_win_flag, &
      Text_Summary_file_flag, &
      Text_Summary_file_interval, &
      Text_Summary_file_dir, &
      Text_Summary_file_prefix, &
      Text_Summary_num_cols, Text_Summary_num_rows, Text_Summary_name, &
      Text_Summary_win_width, &
      Text_Summary_win_aspect_ratio, &
      Text_Summary_file_width, &
      Text_Summary_file_aspect_ratio, &
      Text_Summary_title, Text_Summary_xleft, Text_Summary_xright, &
      Text_Summary_ybot, Text_Summary_ytop, Text_Summary_txt_scale, &

      History_Track_win_flag, &
      History_Track_file_flag, &
      History_Track_file_interval, &
      History_Track_step_min, &
      History_Track_step_max, &
      show_History_Track_target_box, &
      History_Track_n_sigma, &
      History_Track_xtarget, &
      History_Track_xsigma, &
      History_Track_ytarget, &
      History_Track_ysigma, &
      History_Track_xname, &
      History_Track_xaxis_label, &
      History_Track_yname, &
      History_Track_yaxis_label, &
      History_Track_file_dir, &
      History_Track_file_prefix, &
      show_History_Track_annotation1, &
      show_History_Track_annotation2, &
      show_History_Track_annotation3, &
      History_Track_fname, &
      History_Track_reverse_xaxis, &
      History_Track_reverse_yaxis, &
      History_Track_log_xaxis, &
      History_Track_log_yaxis, &
      History_Track_xmin, &
      History_Track_xmax, &
      History_Track_ymin, &
      History_Track_ymax, &
      History_Track_xmargin, &
      History_Track_ymargin, &
      History_Track_dxmin, &
      History_Track_dymin, &
      History_Track_win_width, &
      History_Track_win_aspect_ratio, &
      History_Track_file_width, &
      History_Track_file_aspect_ratio, &
      History_Track_xleft, &
      History_Track_xright, &
      History_Track_ybot, &
      History_Track_ytop, &
      History_Track_txt_scale, &
      History_Track_title, &
      History_Track_use_decorator, &

      Star_History_track_win_flag, &
      Star_History_track_file_flag, &
      Star_History_track_file_interval, &
      Star_History_track_step_min, &
      Star_History_track_step_max, &
      Star_History_track_xname, &
      Star_History_track_xaxis_label, &
      Star_History_track_yname, &
      Star_History_track_yaxis_label, &
      Star_History_track_file_dir, &
      Star_History_track_file_prefix, &
      show_Star_History_track_annotation1, &
      show_Star_History_track_annotation2, &
      show_Star_History_track_annotation3, &
      Star_History_track_fname, &
      Star_History_track_reverse_xaxis, &
      Star_History_track_reverse_yaxis, &
      Star_History_track_log_xaxis, &
      Star_History_track_log_yaxis, &
      Star_History_track_xmin, &
      Star_History_track_xmax, &
      Star_History_track_ymin, &
      Star_History_track_ymax, &
      Star_History_track_xmargin, &
      Star_History_track_ymargin, &
      Star_History_track_dxmin, &
      Star_History_track_dymin, &
      Star_History_track_win_width, &
      Star_History_track_win_aspect_ratio, &
      Star_History_track_file_width, &
      Star_History_track_file_aspect_ratio, &
      Star_History_track_xleft, &
      Star_History_track_xright, &
      Star_History_track_ybot, &
      Star_History_track_ytop, &
      Star_History_track_txt_scale, &
      Star_History_track_title, &
      Star_History_track_use_decorator, &
      
      History_Panels_win_flag, &
      History_Panels_win_width, &
      History_Panels_win_aspect_ratio, &
      History_Panels_xleft, &
      History_Panels_xright, &
      History_Panels_ybot, &
      History_Panels_ytop, &
      History_Panels_txt_scale, &
      History_Panels_title, &
      History_Panels_xmax, &
      History_Panels_xmin, &
      History_Panels_dxmin, &
      History_Panels_max_width, &
      History_Panels_num_panels, &
      History_Panels_xaxis_name, &
      History_Panels_yaxis_name, &
      History_Panels_xaxis_reversed, &
      History_Panels_yaxis_reversed, &
      History_Panels_yaxis_log, &
      History_Panels_ymin, &
      History_Panels_ymax, &
      History_Panels_dymin, &
      History_Panels_other_yaxis_name, &
      History_Panels_other_yaxis_reversed, &
      History_Panels_xaxis_log, &
      History_Panels_other_yaxis_log, &
      History_Panels_other_ymin, &
      History_Panels_other_ymax, &
      History_Panels_other_dymin, &
      History_Panels_points_name, &
      History_Panels_file_flag, &
      History_Panels_file_dir, &
      History_Panels_file_prefix, &
      History_Panels_file_interval, &
      History_Panels_file_width, &
      History_Panels_file_aspect_ratio, &
      History_Panels_xmargin, &
      History_Panels_ymargin, &
      History_Panels_other_ymargin, &
      History_Panels_use_decorator, &

      History_Panel_points_error_bars, &
      History_Panel_points_interval, &
      History_Panel_points_marker, &
      History_Panel_points_ci, &
      History_Panel_points_lw, &
      History_Panel_points_ch, &

      Summary_History_win_flag, &
      Summary_History_file_flag, &
      Summary_History_file_interval, &
      Summary_History_file_dir, &
      Summary_History_file_prefix, &
      Summary_History_scaled_value, &
      Summary_History_xmin, &
      Summary_History_xmax, &
      Summary_History_max_width, &
      Summary_History_win_width, &
      Summary_History_win_aspect_ratio, &
      Summary_History_file_width, &
      Summary_History_file_aspect_ratio, &
      Summary_History_xleft, &
      Summary_History_xright, &
      Summary_History_ybot, &
      Summary_History_ytop, &
      Summary_History_txt_scale, &
      Summary_History_title, &
      Summary_History_name, &
      Summary_History_legend, &
      Summary_History_num_lines, &
      Summary_History_use_decorator, &

      Grid_win_flag, &
      Grid_win_width, &
      Grid_win_aspect_ratio, &
      Grid_xleft, &
      Grid_xright, &
      Grid_ybot, &
      Grid_ytop, &
      Grid_title, &
      Grid_txt_scale_factor, &
      Grid_num_cols, &
      Grid_num_rows, &
      Grid_num_plots, &
      Grid_plot_name, &
      Grid_plot_row, &
      Grid_plot_rowspan, &
      Grid_plot_col, &
      Grid_plot_colspan, &
      Grid_plot_pad_left, &
      Grid_plot_pad_right, &
      Grid_plot_pad_top, &
      Grid_plot_pad_bot, &
      Grid_file_flag, &
      Grid_file_dir, &
      Grid_file_prefix, &
      Grid_file_interval, &
      Grid_file_width, &
      Grid_file_aspect_ratio, &

      Star1_win_flag, &
      Star1_file_flag, &
      Star1_file_interval, &
      Star1_file_dir, &
      Star1_file_prefix, &
      Star1_win_width, &
      Star1_win_aspect_ratio, &
      Star1_xleft, &
      Star1_xright, &
      Star1_ybot, &
      Star1_ytop, &
      Star1_file_width, &
      Star1_file_aspect_ratio, &
      Star1_txt_scale_factor, &
      Star1_title, &
      Star1_plot_name, &
      do_star1_box, &
      star1_box_pad_left, &
      star1_box_pad_right, &
      star1_box_pad_bot, &
      star1_box_pad_top, &

      Star2_win_flag, &
      Star2_file_flag, &
      Star2_file_interval, &
      Star2_file_dir, &
      Star2_file_prefix, &
      Star2_win_width, &
      Star2_win_aspect_ratio, &
      Star2_xleft, &
      Star2_xright, &
      Star2_ybot, &
      Star2_ytop, &
      Star2_file_width, &
      Star2_file_aspect_ratio, &
      Star2_txt_scale_factor, &
      Star2_title, &
      Star2_plot_name, &
      do_star2_box, &
      star2_box_pad_left, &
      star2_box_pad_right, &
      star2_box_pad_bot, &
      star2_box_pad_top, &

      show_mtrans_status, &

      Orbit_win_flag, &
      Orbit_file_flag, &
      Orbit_file_interval, &
      Orbit_file_dir, &
      Orbit_file_prefix, &
      Orbit_title, &
      Orbit_win_width, &
      Orbit_win_aspect_ratio, &
      Orbit_xleft, &
      Orbit_xright, &
      Orbit_ybot, &
      Orbit_ytop, &
      Orbit_file_width, &
      Orbit_file_aspect_ratio, &
      Orbit_txt_scale, &
      Orbit_show_RL, &
      Orbit_show_stars, &

      annotation1_ci, &
      annotation1_ch, &
      annotation1_lw, &
      annotation1_cf, &
      annotation1_text, &
      annotation1_side, &
      annotation1_disp, &
      annotation1_coord, &
      annotation1_fjust, &

      annotation2_ci, &
      annotation2_ch, &
      annotation2_lw, &
      annotation2_cf, &
      annotation2_text, &
      annotation2_side, &
      annotation2_disp, &
      annotation2_coord, &
      annotation2_fjust, &

      annotation3_ci, &
      annotation3_ch, &
      annotation3_lw, &
      annotation3_cf, &
      annotation3_text, &
      annotation3_side, &
      annotation3_disp, &
      annotation3_coord, &
      annotation3_fjust, &

      read_extra_pgbinary_inlist, &
      extra_pgbinary_inlist_name


contains


   subroutine read_pgbinary(b, filename, ierr)
      use binary_private_def
      use utils_lib
      type (binary_info), pointer :: b
      character(*), intent(in) :: filename
      integer, intent(out) :: ierr
      !      character (len = strlen) :: pgbinary_namelist_name
      !      pgbinary_namelist_name = ''
      ierr = 0
      call set_default_pgbinary_controls
      call read_pgbinary_file(b, filename, 1, ierr)
   end subroutine read_pgbinary


   recursive subroutine read_pgbinary_file(b, filename, level, ierr)
      use binary_private_def
      use utils_lib
      character(*), intent(in) :: filename
      type (binary_info), pointer :: b
      integer, intent(in) :: level
      integer, intent(out) :: ierr
      logical, dimension(max_extra_inlists) :: read_extra
      character (len=strlen) :: message
      character (len=strlen), dimension(max_extra_inlists) :: extra
      integer :: unit, i

      ierr = 0

      if (level >= 10) then
         write(*, *) 'ERROR: too many levels of nested extra pgbinary inlist files'
         ierr = -1
         return
      end if
      if (len_trim(filename) > 0) then
         open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', status = 'old', iostat = ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open pgbinary namelist file ', trim(filename)
            return
         end if
         read(unit, nml = pgbinary, iostat = ierr)
         close(unit)
         if (ierr /= 0) then
            write(*, *)
            write(*, *)
            write(*, '(a)') &
               'Failed while trying to read pgbinary namelist file: ' // trim(filename)
            write(*, '(a)') &
               'Perhaps the following runtime error message will help you find the problem.'
            write(*, *)
            open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', status = 'old', iostat = ierr)
            read(unit, nml = pgbinary)
            close(unit)
            return
         end if
      end if

      call store_pgbinary_controls(b, ierr)

 ! recursive calls to read other inlists
         do i=1, max_extra_inlists
            read_extra(i) = read_extra_pgbinary_inlist(i)
            read_extra_pgbinary_inlist(i) = .false.
            extra(i) = extra_pgbinary_inlist_name(i)
            extra_pgbinary_inlist_name(i) = 'undefined'
            
            if (read_extra(i)) then
               call read_pgbinary_file(b, extra(i), level+1, ierr)
               if (ierr /= 0) return
            end if
         end do

   end subroutine read_pgbinary_file


   subroutine store_pgbinary_controls(b, ierr)
      use binary_private_def
      type (binary_info), pointer :: b
      type (pgbinary_controls), pointer :: pg
      integer, intent(out) :: ierr

      ierr = 0
      pg => b% pg

      pg% file_device = file_device
      pg% file_extension = file_extension
      pg% file_digits = file_digits
      pg% pgbinary_interval = pgbinary_interval
      pg% pause = pause
      pg% pause_interval = pause_interval
      pg% pgbinary_sleep = pgbinary_sleep
      pg% clear_history = clear_history

      pg% file_white_on_black_flag = file_white_on_black_flag
      pg% win_white_on_black_flag = win_white_on_black_flag
      pg% pgbinary_show_model_number = pgbinary_show_model_number
      pg% pgbinary_show_age = pgbinary_show_age
      pg% pgbinary_show_age_in_seconds = pgbinary_show_age_in_seconds
      pg% pgbinary_show_age_in_minutes = pgbinary_show_age_in_minutes
      pg% pgbinary_show_age_in_hours = pgbinary_show_age_in_hours
      pg% pgbinary_show_age_in_days = pgbinary_show_age_in_days
      pg% pgbinary_show_age_in_years = pgbinary_show_age_in_years
      pg% pgbinary_show_log_age_in_years = pgbinary_show_log_age_in_years

      pg% pgbinary_report_writing_files = pgbinary_report_writing_files

      pg% pgbinary_show_title = pgbinary_show_title
      pg% pgbinary_title_scale = pgbinary_title_scale
      pg% pgbinary_title_disp = pgbinary_title_disp
      pg% pgbinary_title_coord = pgbinary_title_coord
      pg% pgbinary_title_fjust = pgbinary_title_fjust
      pg% pgbinary_title_lw = pgbinary_title_lw

      pg% pgbinary_grid_show_title = pgbinary_grid_show_title
      pg% pgbinary_grid_title_scale = pgbinary_grid_title_scale
      pg% pgbinary_grid_title_disp = pgbinary_grid_title_disp
      pg% pgbinary_grid_title_coord = pgbinary_grid_title_coord
      pg% pgbinary_grid_title_fjust = pgbinary_grid_title_fjust
      pg% pgbinary_grid_title_lw = pgbinary_grid_title_lw

      pg% pgbinary_age_scale = pgbinary_age_scale
      pg% pgbinary_age_disp = pgbinary_age_disp
      pg% pgbinary_age_coord = pgbinary_age_coord
      pg% pgbinary_age_fjust = pgbinary_age_fjust
      pg% pgbinary_age_lw = pgbinary_age_lw

      pg% pgbinary_model_scale = pgbinary_model_scale
      pg% pgbinary_model_disp = pgbinary_model_disp
      pg% pgbinary_model_coord = pgbinary_model_coord
      pg% pgbinary_model_fjust = pgbinary_model_fjust
      pg% pgbinary_model_lw = pgbinary_model_lw

      pg% pgbinary_xaxis_label_scale = pgbinary_xaxis_label_scale
      pg% pgbinary_left_yaxis_label_scale = pgbinary_left_yaxis_label_scale
      pg% pgbinary_right_yaxis_label_scale = pgbinary_right_yaxis_label_scale
      pg% pgbinary_xaxis_label_lw = pgbinary_xaxis_label_lw
      pg% pgbinary_left_yaxis_label_lw = pgbinary_left_yaxis_label_lw
      pg% pgbinary_right_yaxis_label_lw = pgbinary_right_yaxis_label_lw
      pg% pgbinary_xaxis_label_disp = pgbinary_xaxis_label_disp
      pg% pgbinary_left_yaxis_label_disp = pgbinary_left_yaxis_label_disp
      pg% pgbinary_right_yaxis_label_disp = pgbinary_right_yaxis_label_disp
      pg% pgbinary_num_scale = pgbinary_num_scale
      pg% pgbinary_lw = pgbinary_lw
      pg% pgbinary_model_lw = pgbinary_model_lw
      pg% pgbinary_box_lw = pgbinary_box_lw

      pg% Text_Summary_win_flag = Text_Summary_win_flag
      pg% Text_Summary_file_flag = Text_Summary_file_flag
      pg% Text_Summary_file_interval = Text_Summary_file_interval
      pg% Text_Summary_file_dir = Text_Summary_file_dir
      pg% Text_Summary_file_prefix = Text_Summary_file_prefix
      pg% Text_Summary_num_cols = Text_Summary_num_cols
      pg% Text_Summary_num_rows = Text_Summary_num_rows
      pg% Text_Summary_name = Text_Summary_name
      pg% Text_Summary_win_width = Text_Summary_win_width
      pg% Text_Summary_win_aspect_ratio = Text_Summary_win_aspect_ratio
      pg% Text_Summary_file_width = Text_Summary_file_width
      pg% Text_Summary_file_aspect_ratio = Text_Summary_file_aspect_ratio
      pg% Text_Summary_title = Text_Summary_title
      pg% Text_Summary_xleft = Text_Summary_xleft
      pg% Text_Summary_xright = Text_Summary_xright
      pg% Text_Summary_ybot = Text_Summary_ybot
      pg% Text_Summary_ytop = Text_Summary_ytop
      pg% Text_Summary_txt_scale = Text_Summary_txt_scale

      pg% History_track_win_flag = History_track_win_flag
      pg% History_track_file_flag = History_track_file_flag
      pg% History_track_file_interval = History_track_file_interval
      pg% History_track_step_min = History_track_step_min
      pg% History_track_step_max = History_track_step_max
      pg% show_History_track_target_box = show_History_track_target_box
      pg% History_track_n_sigma = History_track_n_sigma
      pg% History_track_xtarget = History_track_xtarget
      pg% History_track_xsigma = History_track_xsigma
      pg% History_track_ytarget = History_track_ytarget
      pg% History_track_ysigma = History_track_ysigma
      pg% History_track_file_dir = History_track_file_dir
      pg% History_track_file_prefix = History_track_file_prefix
      pg% show_History_track_annotation1 = show_History_track_annotation1
      pg% show_History_track_annotation2 = show_History_track_annotation2
      pg% show_History_track_annotation3 = show_History_track_annotation3
      pg% History_track_fname = History_track_fname
      pg% History_track_xname = History_track_xname
      pg% History_track_xaxis_label = History_track_xaxis_label
      pg% History_track_yname = History_track_yname
      pg% History_track_yaxis_label = History_track_yaxis_label
      pg% History_track_reverse_xaxis = History_track_reverse_xaxis
      pg% History_track_reverse_yaxis = History_track_reverse_yaxis
      pg% History_track_log_xaxis = History_track_log_xaxis
      pg% History_track_log_yaxis = History_track_log_yaxis
      pg% History_track_xmin = History_track_xmin
      pg% History_track_xmax = History_track_xmax
      pg% History_track_ymin = History_track_ymin
      pg% History_track_ymax = History_track_ymax
      pg% History_track_xmargin = History_track_xmargin
      pg% History_track_ymargin = History_track_ymargin
      pg% History_track_dxmin = History_track_dxmin
      pg% History_track_dymin = History_track_dymin
      pg% History_track_win_width = History_track_win_width
      pg% History_track_win_aspect_ratio = History_track_win_aspect_ratio
      pg% History_track_file_width = History_track_file_width
      pg% History_track_file_aspect_ratio = History_track_file_aspect_ratio
      pg% History_track_xleft = History_track_xleft
      pg% History_track_xright = History_track_xright
      pg% History_track_ybot = History_track_ybot
      pg% History_track_ytop = History_track_ytop
      pg% History_track_txt_scale = History_track_txt_scale
      pg% History_track_title = History_track_title
      pg% History_track_use_decorator = History_track_use_decorator
      pg% Star_History_track_win_flag = Star_History_track_win_flag
      pg% Star_History_track_file_flag = Star_History_track_file_flag
      pg% Star_History_track_file_interval = Star_History_track_file_interval
      pg% Star_History_track_step_min = Star_History_track_step_min
      pg% Star_History_track_step_max = Star_History_track_step_max
      pg% Star_History_track_file_dir = Star_History_track_file_dir
      pg% Star_History_track_file_prefix = Star_History_track_file_prefix
      pg% show_Star_History_track_annotation1 = show_Star_History_track_annotation1
      pg% show_Star_History_track_annotation2 = show_Star_History_track_annotation2
      pg% show_Star_History_track_annotation3 = show_Star_History_track_annotation3
      pg% Star_History_track_fname = Star_History_track_fname
      pg% Star_History_track_xname = Star_History_track_xname
      pg% Star_History_track_xaxis_label = Star_History_track_xaxis_label
      pg% Star_History_track_yname = Star_History_track_yname
      pg% Star_History_track_yaxis_label = Star_History_track_yaxis_label
      pg% Star_History_track_reverse_xaxis = Star_History_track_reverse_xaxis
      pg% Star_History_track_reverse_yaxis = Star_History_track_reverse_yaxis
      pg% Star_History_track_log_xaxis = Star_History_track_log_xaxis
      pg% Star_History_track_log_yaxis = Star_History_track_log_yaxis
      pg% Star_History_track_xmin = Star_History_track_xmin
      pg% Star_History_track_xmax = Star_History_track_xmax
      pg% Star_History_track_ymin = Star_History_track_ymin
      pg% Star_History_track_ymax = Star_History_track_ymax
      pg% Star_History_track_xmargin = Star_History_track_xmargin
      pg% Star_History_track_ymargin = Star_History_track_ymargin
      pg% Star_History_track_dxmin = Star_History_track_dxmin
      pg% Star_History_track_dymin = Star_History_track_dymin
      pg% Star_History_track_win_width = Star_History_track_win_width
      pg% Star_History_track_win_aspect_ratio = Star_History_track_win_aspect_ratio
      pg% Star_History_track_file_width = Star_History_track_file_width
      pg% Star_History_track_file_aspect_ratio = Star_History_track_file_aspect_ratio
      pg% Star_History_track_xleft = Star_History_track_xleft
      pg% Star_History_track_xright = Star_History_track_xright
      pg% Star_History_track_ybot = Star_History_track_ybot
      pg% Star_History_track_ytop = Star_History_track_ytop
      pg% Star_History_track_txt_scale = Star_History_track_txt_scale
      pg% Star_History_track_title = Star_History_track_title
      pg% Star_History_track_use_decorator = Star_History_track_use_decorator

      pg% History_Panels_win_flag = History_Panels_win_flag
      pg% History_Panels_win_width = History_Panels_win_width
      pg% History_Panels_win_aspect_ratio = History_Panels_win_aspect_ratio
      pg% History_Panels_xleft = History_Panels_xleft
      pg% History_Panels_xright = History_Panels_xright
      pg% History_Panels_ybot = History_Panels_ybot
      pg% History_Panels_ytop = History_Panels_ytop
      pg% History_Panels_txt_scale = History_Panels_txt_scale
      pg% History_Panels_title = History_Panels_title
      pg% History_Panels_xmax = History_Panels_xmax
      pg% History_Panels_xmin = History_Panels_xmin
      pg% History_Panels_dxmin = History_Panels_dxmin
      pg% History_Panels_max_width = History_Panels_max_width
      pg% History_Panels_num_panels = History_Panels_num_panels
      pg% History_Panels_xaxis_name = History_Panels_xaxis_name
      pg% History_Panels_yaxis_name = History_Panels_yaxis_name
      pg% History_Panels_xaxis_reversed = History_Panels_xaxis_reversed
      pg% History_Panels_yaxis_reversed = History_Panels_yaxis_reversed
      pg% History_Panels_xaxis_log = History_Panels_xaxis_log
      pg% History_Panels_yaxis_log = History_Panels_yaxis_log
      pg% History_Panels_ymin = History_Panels_ymin
      pg% History_Panels_ymax = History_Panels_ymax
      pg% History_Panels_dymin = History_Panels_dymin
      pg% History_Panels_other_yaxis_name = History_Panels_other_yaxis_name
      pg% History_Panels_other_yaxis_reversed = History_Panels_other_yaxis_reversed
      pg% History_Panels_other_yaxis_log = History_Panels_other_yaxis_log
      pg% History_Panels_other_ymin = History_Panels_other_ymin
      pg% History_Panels_other_ymax = History_Panels_other_ymax
      pg% History_Panels_other_dymin = History_Panels_other_dymin
      pg% History_Panels_file_flag = History_Panels_file_flag
      pg% History_Panels_points_name = History_Panels_points_name
      pg% History_Panels_file_dir = History_Panels_file_dir
      pg% History_Panels_file_prefix = History_Panels_file_prefix
      pg% History_Panels_file_interval = History_Panels_file_interval
      pg% History_Panels_file_width = History_Panels_file_width
      pg% History_Panels_file_aspect_ratio = History_Panels_file_aspect_ratio
      pg% History_Panels_xmargin = History_Panels_xmargin
      pg% History_Panels_ymargin = History_Panels_ymargin
      pg% History_Panels_other_ymargin = History_Panels_other_ymargin
      pg% History_Panels_use_decorator = History_Panels_use_decorator

      pg% History_Panel_points_error_bars = History_Panel_points_error_bars
      pg% History_Panel_points_interval = History_Panel_points_interval
      pg% History_Panel_points_marker = History_Panel_points_marker
      pg% History_Panel_points_ci = History_Panel_points_ci
      pg% History_Panel_points_lw = History_Panel_points_lw
      pg% History_Panel_points_ch = History_Panel_points_ch

      pg% Summary_History_win_flag = Summary_History_win_flag
      pg% Summary_History_file_flag = Summary_History_file_flag
      pg% Summary_History_file_interval = Summary_History_file_interval
      pg% Summary_History_file_dir = Summary_History_file_dir
      pg% Summary_History_file_prefix = Summary_History_file_prefix
      pg% Summary_History_scaled_value = Summary_History_scaled_value
      pg% Summary_History_xmin = Summary_History_xmin
      pg% Summary_History_xmax = Summary_History_xmax
      pg% Summary_History_max_width = Summary_History_max_width
      pg% Summary_History_win_width = Summary_History_win_width
      pg% Summary_History_win_aspect_ratio = Summary_History_win_aspect_ratio
      pg% Summary_History_file_width = Summary_History_file_width
      pg% Summary_History_file_aspect_ratio = Summary_History_file_aspect_ratio
      pg% Summary_History_xleft = Summary_History_xleft
      pg% Summary_History_xright = Summary_History_xright
      pg% Summary_History_ybot = Summary_History_ybot
      pg% Summary_History_ytop = Summary_History_ytop
      pg% Summary_History_txt_scale = Summary_History_txt_scale
      pg% Summary_History_title = Summary_History_title
      pg% Summary_History_name = Summary_History_name
      pg% Summary_History_legend = Summary_History_legend
      pg% Summary_History_num_lines = Summary_History_num_lines
      pg% Summary_History_use_decorator = Summary_History_use_decorator

      pg% Grid_win_flag = Grid_win_flag
      pg% Grid_win_width = Grid_win_width
      pg% Grid_win_aspect_ratio = Grid_win_aspect_ratio
      pg% Grid_xleft = Grid_xleft
      pg% Grid_xright = Grid_xright
      pg% Grid_ybot = Grid_ybot
      pg% Grid_ytop = Grid_ytop
      pg% Grid_title = Grid_title
      pg% Grid_txt_scale_factor = Grid_txt_scale_factor
      pg% Grid_num_cols = Grid_num_cols
      pg% Grid_num_rows = Grid_num_rows
      pg% Grid_num_plots = Grid_num_plots
      pg% Grid_plot_name = Grid_plot_name
      pg% Grid_plot_row = Grid_plot_row
      pg% Grid_plot_rowspan = Grid_plot_rowspan
      pg% Grid_plot_col = Grid_plot_col
      pg% Grid_plot_colspan = Grid_plot_colspan
      pg% Grid_plot_pad_left = Grid_plot_pad_left
      pg% Grid_plot_pad_right = Grid_plot_pad_right
      pg% Grid_plot_pad_top = Grid_plot_pad_top
      pg% Grid_plot_pad_bot = Grid_plot_pad_bot
      pg% Grid_file_flag = Grid_file_flag
      pg% Grid_file_dir = Grid_file_dir
      pg% Grid_file_prefix = Grid_file_prefix
      pg% Grid_file_interval = Grid_file_interval
      pg% Grid_file_width = Grid_file_width
      pg% Grid_file_aspect_ratio = Grid_file_aspect_ratio

      pg% Star1_win_flag = Star1_win_flag
      pg% Star1_file_flag = Star1_file_flag
      pg% Star1_file_interval = Star1_file_interval
      pg% Star1_file_dir = Star1_file_dir
      pg% Star1_file_prefix = Star1_file_prefix
      pg% Star1_win_width = Star1_win_width
      pg% Star1_win_aspect_ratio = Star1_win_aspect_ratio
      pg% Star1_xleft = Star1_xleft
      pg% Star1_xright = Star1_xright
      pg% Star1_ybot = Star1_ybot
      pg% Star1_ytop = Star1_ytop
      pg% Star1_file_width = Star1_file_width
      pg% Star1_file_aspect_ratio = Star1_file_aspect_ratio
      pg% Star1_txt_scale_factor = Star1_txt_scale_factor
      pg% Star1_title = Star1_title
      pg% Star1_plot_name = Star1_plot_name
      pg% do_star1_box = do_star1_box
      pg% star1_box_pad_left = star1_box_pad_left
      pg% star1_box_pad_right = star1_box_pad_right
      pg% star1_box_pad_bot = star1_box_pad_bot
      pg% star1_box_pad_top = star1_box_pad_top

      pg% Star2_win_flag = Star2_win_flag
      pg% Star2_file_flag = Star2_file_flag
      pg% Star2_file_interval = Star2_file_interval
      pg% Star2_file_dir = Star2_file_dir
      pg% Star2_file_prefix = Star2_file_prefix
      pg% Star2_win_width = Star2_win_width
      pg% Star2_win_aspect_ratio = Star2_win_aspect_ratio
      pg% Star2_xleft = Star2_xleft
      pg% Star2_xright = Star2_xright
      pg% Star2_ybot = Star2_ybot
      pg% Star2_ytop = Star2_ytop
      pg% Star2_file_width = Star2_file_width
      pg% Star2_file_aspect_ratio = Star2_file_aspect_ratio
      pg% Star2_txt_scale_factor = Star2_txt_scale_factor
      pg% Star2_title = Star2_title
      pg% Star2_plot_name = Star2_plot_name
      pg% do_star2_box = do_star2_box
      pg% star2_box_pad_left = star2_box_pad_left
      pg% star2_box_pad_right = star2_box_pad_right
      pg% star2_box_pad_bot = star2_box_pad_bot
      pg% star2_box_pad_top = star2_box_pad_top

      pg% show_mtrans_status = show_mtrans_status

      pg% Orbit_win_flag = Orbit_win_flag
      pg% Orbit_file_flag = Orbit_file_flag
      pg% Orbit_file_interval = Orbit_file_interval
      pg% Orbit_file_dir = Orbit_file_dir
      pg% Orbit_file_prefix = Orbit_file_prefix
      pg% Orbit_title = Orbit_title
      pg% Orbit_win_width = Orbit_win_width
      pg% Orbit_win_aspect_ratio = Orbit_win_aspect_ratio
      pg% Orbit_xleft = Orbit_xleft
      pg% Orbit_xright = Orbit_xright
      pg% Orbit_ybot = Orbit_ybot
      pg% Orbit_ytop = Orbit_ytop
      pg% Orbit_file_width = Orbit_file_width
      pg% Orbit_file_aspect_ratio = Orbit_file_aspect_ratio
      pg% Orbit_txt_scale = Orbit_txt_scale
      pg% Orbit_show_RL = Orbit_show_RL
      pg% Orbit_show_stars = Orbit_show_stars

      pg% annotation1_ci = annotation1_ci
      pg% annotation1_ch = annotation1_ch
      pg% annotation1_lw = annotation1_lw
      pg% annotation1_cf = annotation1_cf
      pg% annotation1_text = annotation1_text
      pg% annotation1_side = annotation1_side
      pg% annotation1_disp = annotation1_disp
      pg% annotation1_coord = annotation1_coord
      pg% annotation1_fjust = annotation1_fjust

      pg% annotation2_ci = annotation2_ci
      pg% annotation2_ch = annotation2_ch
      pg% annotation2_lw = annotation2_lw
      pg% annotation2_cf = annotation2_cf
      pg% annotation2_text = annotation2_text
      pg% annotation2_side = annotation2_side
      pg% annotation2_disp = annotation2_disp
      pg% annotation2_coord = annotation2_coord
      pg% annotation2_fjust = annotation2_fjust

      pg% annotation3_ci = annotation3_ci
      pg% annotation3_ch = annotation3_ch
      pg% annotation3_lw = annotation3_lw
      pg% annotation3_cf = annotation3_cf
      pg% annotation3_text = annotation3_text
      pg% annotation3_side = annotation3_side
      pg% annotation3_disp = annotation3_disp
      pg% annotation3_coord = annotation3_coord
      pg% annotation3_fjust = annotation3_fjust

      pg% read_extra_pgbinary_inlist = read_extra_pgbinary_inlist
      pg% extra_pgbinary_inlist_name = extra_pgbinary_inlist_name

   end subroutine store_pgbinary_controls


   subroutine set_default_pgbinary_controls

      include 'pgbinary.defaults'

   end subroutine set_default_pgbinary_controls


end module pgbinary_ctrls_io

