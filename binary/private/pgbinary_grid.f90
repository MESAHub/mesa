! ***********************************************************************
!
!   Copyright (C) 2014-2022  Bill Paxton, Matthias Fabry
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
!   Foundation, Inc., 59 Temple Place, Suite 440, Boston, MA 02111-1407 USA
!
! ***********************************************************************

module pgbinary_grid

   use binary_def
   use const_def
   use pgbinary_support

   implicit none


contains


   subroutine grid1_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call Grid_plot(b, id, device_id, &
         b% Grid1_xleft, b% Grid1_xright, &
         b% Grid1_ybot, b% Grid1_ytop, .false., b% Grid1_title, &
         b% Grid1_txt_scale_factor, &
         b% Grid1_num_cols, &
         b% Grid1_num_rows, &
         b% Grid1_num_plots, &
         b% Grid1_plot_name, &
         b% Grid1_plot_row, &
         b% Grid1_plot_rowspan, &
         b% Grid1_plot_col, &
         b% Grid1_plot_colspan, &
         b% Grid1_plot_pad_left, &
         b% Grid1_plot_pad_right, &
         b% Grid1_plot_pad_top, &
         b% Grid1_plot_pad_bot, &
         ierr)
   end subroutine grid1_plot


   subroutine grid2_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call Grid_plot(b, id, device_id, &
         b% Grid2_xleft, b% Grid2_xright, &
         b% Grid2_ybot, b% Grid2_ytop, .false., b% Grid2_title, &
         b% Grid2_txt_scale_factor, &
         b% Grid2_num_cols, &
         b% Grid2_num_rows, &
         b% Grid2_num_plots, &
         b% Grid2_plot_name, &
         b% Grid2_plot_row, &
         b% Grid2_plot_rowspan, &
         b% Grid2_plot_col, &
         b% Grid2_plot_colspan, &
         b% Grid2_plot_pad_left, &
         b% Grid2_plot_pad_right, &
         b% Grid2_plot_pad_top, &
         b% Grid2_plot_pad_bot, &
         ierr)
   end subroutine grid2_plot


   subroutine grid3_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call Grid_plot(b, id, device_id, &
         b% Grid3_xleft, b% Grid3_xright, &
         b% Grid3_ybot, b% Grid3_ytop, .false., b% Grid3_title, &
         b% Grid3_txt_scale_factor, &
         b% Grid3_num_cols, &
         b% Grid3_num_rows, &
         b% Grid3_num_plots, &
         b% Grid3_plot_name, &
         b% Grid3_plot_row, &
         b% Grid3_plot_rowspan, &
         b% Grid3_plot_col, &
         b% Grid3_plot_colspan, &
         b% Grid3_plot_pad_left, &
         b% Grid3_plot_pad_right, &
         b% Grid3_plot_pad_top, &
         b% Grid3_plot_pad_bot, &
         ierr)
   end subroutine grid3_plot


   subroutine grid4_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call Grid_plot(b, id, device_id, &
         b% Grid4_xleft, b% Grid4_xright, &
         b% Grid4_ybot, b% Grid4_ytop, .false., b% Grid4_title, &
         b% Grid4_txt_scale_factor, &
         b% Grid4_num_cols, &
         b% Grid4_num_rows, &
         b% Grid4_num_plots, &
         b% Grid4_plot_name, &
         b% Grid4_plot_row, &
         b% Grid4_plot_rowspan, &
         b% Grid4_plot_col, &
         b% Grid4_plot_colspan, &
         b% Grid4_plot_pad_left, &
         b% Grid4_plot_pad_right, &
         b% Grid4_plot_pad_top, &
         b% Grid4_plot_pad_bot, &
         ierr)
   end subroutine grid4_plot


   subroutine grid5_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call Grid_plot(b, id, device_id, &
         b% Grid5_xleft, b% Grid5_xright, &
         b% Grid5_ybot, b% Grid5_ytop, .false., b% Grid5_title, &
         b% Grid5_txt_scale_factor, &
         b% Grid5_num_cols, &
         b% Grid5_num_rows, &
         b% Grid5_num_plots, &
         b% Grid5_plot_name, &
         b% Grid5_plot_row, &
         b% Grid5_plot_rowspan, &
         b% Grid5_plot_col, &
         b% Grid5_plot_colspan, &
         b% Grid5_plot_pad_left, &
         b% Grid5_plot_pad_right, &
         b% Grid5_plot_pad_top, &
         b% Grid5_plot_pad_bot, &
         ierr)
   end subroutine grid5_plot


   subroutine grid6_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call Grid_plot(b, id, device_id, &
         b% Grid6_xleft, b% Grid6_xright, &
         b% Grid6_ybot, b% Grid6_ytop, .false., b% Grid6_title, &
         b% Grid6_txt_scale_factor, &
         b% Grid6_num_cols, &
         b% Grid6_num_rows, &
         b% Grid6_num_plots, &
         b% Grid6_plot_name, &
         b% Grid6_plot_row, &
         b% Grid6_plot_rowspan, &
         b% Grid6_plot_col, &
         b% Grid6_plot_colspan, &
         b% Grid6_plot_pad_left, &
         b% Grid6_plot_pad_right, &
         b% Grid6_plot_pad_top, &
         b% Grid6_plot_pad_bot, &
         ierr)
   end subroutine grid6_plot


   subroutine grid7_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call Grid_plot(b, id, device_id, &
         b% Grid7_xleft, b% Grid7_xright, &
         b% Grid7_ybot, b% Grid7_ytop, .false., b% Grid7_title, &
         b% Grid7_txt_scale_factor, &
         b% Grid7_num_cols, &
         b% Grid7_num_rows, &
         b% Grid7_num_plots, &
         b% Grid7_plot_name, &
         b% Grid7_plot_row, &
         b% Grid7_plot_rowspan, &
         b% Grid7_plot_col, &
         b% Grid7_plot_colspan, &
         b% Grid7_plot_pad_left, &
         b% Grid7_plot_pad_right, &
         b% Grid7_plot_pad_top, &
         b% Grid7_plot_pad_bot, &
         ierr)
   end subroutine grid7_plot


   subroutine grid8_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call Grid_plot(b, id, device_id, &
         b% Grid8_xleft, b% Grid8_xright, &
         b% Grid8_ybot, b% Grid8_ytop, .false., b% Grid8_title, &
         b% Grid8_txt_scale_factor, &
         b% Grid8_num_cols, &
         b% Grid8_num_rows, &
         b% Grid8_num_plots, &
         b% Grid8_plot_name, &
         b% Grid8_plot_row, &
         b% Grid8_plot_rowspan, &
         b% Grid8_plot_col, &
         b% Grid8_plot_colspan, &
         b% Grid8_plot_pad_left, &
         b% Grid8_plot_pad_right, &
         b% Grid8_plot_pad_top, &
         b% Grid8_plot_pad_bot, &
         ierr)
   end subroutine grid8_plot


   subroutine grid9_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call Grid_plot(b, id, device_id, &
         b% Grid9_xleft, b% Grid9_xright, &
         b% Grid9_ybot, b% Grid9_ytop, .false., b% Grid9_title, &
         b% Grid9_txt_scale_factor, &
         b% Grid9_num_cols, &
         b% Grid9_num_rows, &
         b% Grid9_num_plots, &
         b% Grid9_plot_name, &
         b% Grid9_plot_row, &
         b% Grid9_plot_rowspan, &
         b% Grid9_plot_col, &
         b% Grid9_plot_colspan, &
         b% Grid9_plot_pad_left, &
         b% Grid9_plot_pad_right, &
         b% Grid9_plot_pad_top, &
         b% Grid9_plot_pad_bot, &
         ierr)
   end subroutine grid9_plot


   subroutine Grid_plot(b, id, device_id, &
      Grid_xleft, Grid_xright, &
      Grid_ybot, Grid_ytop, subplot, Grid_title, &
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
      ierr)

      use utils_lib, only : StrLowCase
      use pgbinary_summary_history, only : do_summary_history_plot
      use pgbinary_summary, only : &
         do_Text_Summary1_plot, do_Text_Summary2_plot, do_Text_Summary3_plot, &
         do_Text_Summary4_plot, do_Text_Summary5_plot, do_Text_Summary6_plot, &
         do_Text_Summary7_plot, do_Text_Summary8_plot, do_Text_Summary9_plot
      use pgbinary_history_panels, only : &
         do_History_Panels1_plot, do_History_Panels2_plot, do_History_Panels3_plot, &
         do_History_Panels4_plot, do_History_Panels5_plot, do_History_Panels6_plot, &
         do_History_Panels7_plot, do_History_Panels8_plot, do_History_Panels9_plot
      use pgbinary_hist_track, only : &
         do_History_Track1_plot, do_History_Track2_plot, do_History_Track3_plot, &
         do_History_Track4_plot, do_History_Track5_plot, do_History_Track6_plot, &
         do_History_Track7_plot, do_History_Track8_plot, do_History_Track9_plot
      use pgbinary_star, only : &
         do_Star1_plot, do_Star2_plot

      type (binary_info), pointer :: b
      logical, intent(in) :: subplot
      integer, intent(in) :: id, device_id, &
         Grid_num_cols, &
         Grid_num_rows, &
         Grid_num_plots, &
         Grid_plot_row(:), &
         Grid_plot_rowspan(:), &
         Grid_plot_col(:), &
         Grid_plot_colspan(:)
      real, intent(in) :: &
         Grid_xleft, Grid_xright, &
         Grid_ybot, Grid_ytop, &
         Grid_txt_scale_factor(:), &
         Grid_plot_pad_left(:), &
         Grid_plot_pad_right(:), &
         Grid_plot_pad_top(:), &
         Grid_plot_pad_bot(:)
      character (len = *) :: Grid_title, Grid_plot_name(:)
      integer, intent(out) :: ierr

      integer :: i, j, plot_id
      logical :: found_it
      real :: xleft, xright, ybot, ytop
      real :: row_height, col_width
      type (pgbinary_win_file_data), pointer :: p
      logical, parameter :: grid_subplot = .true.

      include 'formats'

      ierr = 0
      if (Grid_num_plots <= 0 .or. &
         Grid_num_cols <= 0 .or. Grid_num_rows <= 0) return

      col_width = (Grid_xright - Grid_xleft) / Grid_num_cols
      row_height = (Grid_ytop - Grid_ybot) / Grid_num_rows

      if (col_width <= 0d0 .or. row_height <= 0d0) then
         ierr = -1
         write(*, 1) 'Grid: col_width', col_width
         write(*, 1) 'row_height', row_height
         write(*, 1) 'Grid_xleft', Grid_xleft
         write(*, 1) 'Grid_xright', Grid_xright
         write(*, 1) 'Grid_ybot', Grid_ybot
         write(*, 1) 'Grid_ytop', Grid_ytop
         write(*, 2) 'Grid_num_cols', Grid_num_cols
         write(*, 2) 'Grid_num_rows', Grid_num_rows
         return
      end if

      call pgslct(device_id)
      call pgbbuf()
      call pgeras()

      call pgsave
      call pgsvp(Grid_xleft, Grid_xright, Grid_ybot, Grid_ytop)
      if (.not. subplot) then
         call show_model_number_pgbinary(b)
         call show_age_pgbinary(b)
      end if
      call show_grid_title_pgbinary(b, Grid_title)
      call pgunsa

      do i = 1, Grid_num_plots

         if (len_trim(Grid_plot_name(i))==0) exit

         xleft = Grid_xleft + col_width * (Grid_plot_col(i) - 1)
         xright = xleft + col_width * Grid_plot_colspan(i)

         ytop = Grid_ytop - row_height * (Grid_plot_row(i) - 1)
         ybot = ytop - row_height * Grid_plot_rowspan(i)

         xleft = xleft + Grid_plot_pad_left(i)
         xright = xright - Grid_plot_pad_right(i)
         ybot = ybot + Grid_plot_pad_bot(i)
         ytop = ytop - Grid_plot_pad_top(i)

         if (xright <= xleft .or. ytop <= ybot) then
            write(*, 2) 'Bad pgbinary grid spec', i
            write(*, *) 'xright <= xleft', xright <= xleft
            write(*, *) 'ytop <= ybot', ytop <= ybot
            write(*, 2) 'xleft', i, xleft
            write(*, 2) 'xright', i, xright
            write(*, 2) 'ybot', i, ybot
            write(*, 2) 'ytop', i, ytop
            write(*, 2) 'Grid_plot_pad_left(i)', i, Grid_plot_pad_left(i)
            write(*, 2) 'Grid_plot_pad_right(i)', i, Grid_plot_pad_right(i)
            write(*, 2) 'Grid_plot_pad_top(i)', i, Grid_plot_pad_top(i)
            write(*, 2) 'Grid_plot_pad_bot(i)', i, Grid_plot_pad_bot(i)
            write(*, 2) 'col_width', i, col_width
            write(*, 2) 'row_height', i, row_height
            exit
         end if

         call pgsave

         select case(StrLowCase(Grid_plot_name(i)))
         case ('summary_history')
            call do_summary_history_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Summary_History_title, &
               Grid_txt_scale_factor(i) * b% Summary_History_txt_scale, ierr)
         case ('history_panels1')
            call do_History_Panels1_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Panels1_title, &
               Grid_txt_scale_factor(i) * b% History_Panels1_txt_scale, ierr)
         case ('history_panels2')
            call do_History_Panels2_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Panels2_title, &
               Grid_txt_scale_factor(i) * b% History_Panels2_txt_scale, ierr)
         case ('history_panels3')
            call do_History_Panels3_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Panels3_title, &
               Grid_txt_scale_factor(i) * b% History_Panels3_txt_scale, ierr)
         case ('history_panels4')
            call do_History_Panels4_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Panels4_title, &
               Grid_txt_scale_factor(i) * b% History_Panels4_txt_scale, ierr)
         case ('history_panels5')
            call do_History_Panels5_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Panels5_title, &
               Grid_txt_scale_factor(i) * b% History_Panels5_txt_scale, ierr)
         case ('history_panels6')
            call do_History_Panels6_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Panels6_title, &
               Grid_txt_scale_factor(i) * b% History_Panels6_txt_scale, ierr)
         case ('history_panels7')
            call do_History_Panels7_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Panels7_title, &
               Grid_txt_scale_factor(i) * b% History_Panels7_txt_scale, ierr)
         case ('history_panels8')
            call do_History_Panels8_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Panels8_title, &
               Grid_txt_scale_factor(i) * b% History_Panels8_txt_scale, ierr)
         case ('history_panels9')
            call do_History_Panels9_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Panels9_title, &
               Grid_txt_scale_factor(i) * b% History_Panels9_txt_scale, ierr)
         case ('history_track1')
            call do_History_Track1_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Track1_title, &
               Grid_txt_scale_factor(i) * b% History_Track1_txt_scale, ierr)
         case ('history_track2')
            call do_History_Track2_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Track2_title, &
               Grid_txt_scale_factor(i) * b% History_Track2_txt_scale, ierr)
         case ('history_track3')
            call do_History_Track3_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Track3_title, &
               Grid_txt_scale_factor(i) * b% History_Track3_txt_scale, ierr)
         case ('history_track4')
            call do_History_Track4_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Track4_title, &
               Grid_txt_scale_factor(i) * b% History_Track4_txt_scale, ierr)
         case ('history_track5')
            call do_History_Track5_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Track5_title, &
               Grid_txt_scale_factor(i) * b% History_Track5_txt_scale, ierr)
         case ('history_track6')
            call do_History_Track6_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Track6_title, &
               Grid_txt_scale_factor(i) * b% History_Track6_txt_scale, ierr)
         case ('history_track7')
            call do_History_Track7_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Track7_title, &
               Grid_txt_scale_factor(i) * b% History_Track7_txt_scale, ierr)
         case ('history_track8')
            call do_History_Track8_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Track8_title, &
               Grid_txt_scale_factor(i) * b% History_Track8_txt_scale, ierr)
         case ('history_track9')
            call do_History_Track9_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% History_Track9_title, &
               Grid_txt_scale_factor(i) * b% History_Track9_txt_scale, ierr)
         case ('text_summary1')
            call do_Text_Summary1_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Text_Summary1_title, &
               Grid_txt_scale_factor(i) * b% Text_Summary1_txt_scale, ierr)
         case ('text_summary2')
            call do_Text_Summary2_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Text_Summary2_title, &
               Grid_txt_scale_factor(i) * b% Text_Summary2_txt_scale, ierr)
         case ('text_summary3')
            call do_Text_Summary3_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Text_Summary3_title, &
               Grid_txt_scale_factor(i) * b% Text_Summary3_txt_scale, ierr)
         case ('text_summary4')
            call do_Text_Summary4_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Text_Summary4_title, &
               Grid_txt_scale_factor(i) * b% Text_Summary4_txt_scale, ierr)
         case ('text_summary5')
            call do_Text_Summary5_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Text_Summary5_title, &
               Grid_txt_scale_factor(i) * b% Text_Summary5_txt_scale, ierr)
         case ('text_summary6')
            call do_Text_Summary6_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Text_Summary6_title, &
               Grid_txt_scale_factor(i) * b% Text_Summary6_txt_scale, ierr)
         case ('text_summary7')
            call do_Text_Summary7_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Text_Summary7_title, &
               Grid_txt_scale_factor(i) * b% Text_Summary7_txt_scale, ierr)
         case ('text_summary8')
            call do_Text_Summary8_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Text_Summary8_title, &
               Grid_txt_scale_factor(i) * b% Text_Summary8_txt_scale, ierr)
         case ('text_summary9')
            call do_Text_Summary9_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Text_Summary9_title, &
               Grid_txt_scale_factor(i) * b% Text_Summary9_txt_scale, ierr)
         case ('star1')
            call do_Star1_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Star1_Title, &
               Grid_txt_scale_factor(i) * b% Star1_txt_scale_factor, b% Star1_plot_name, ierr)
         case ('star2')
            call do_Star2_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% Star2_Title, &
               Grid_txt_scale_factor(i) * b% Star2_txt_scale_factor, b% Star2_plot_name, ierr)

         case default
            ! check for "other" plot
            found_it = .false.
            do j = 1, max_num_Other_plots
               plot_id = i_Other + j - 1
               p => b% pgbinary_win_file_ptr(plot_id)
               if (p% okay_to_call_do_plot_in_binary_grid .and. &
                  StrLowCase(p% name) == StrLowCase(Grid_plot_name(i))) then
                  call p% do_plot_in_binary_grid(&
                     id, device_id, xleft, xright, ybot, ytop, &
                     Grid_txt_scale_factor(i), ierr)
                  found_it = .true.
                  exit
               end if
            end do

            if (.not. found_it) then

               write(*, *) 'FAILED TO RECOGNIZE NAME FOR GRID PLOT: ' // trim(Grid_plot_name(i))
               write(*, '(a)') &
                  'here are the valid names:', &
                  'Summary_History', &
                  'Text_Summary1,..,9', &
                  'History_Panels1,..,9', &
                  'History_Tracks1,..,9', &
                  'Star1,2'
               write(*, *)

            end if

         end select

         call pgunsa

         if (ierr /= 0) exit

      end do

      call pgebuf()

   end subroutine Grid_plot


end module pgbinary_grid

