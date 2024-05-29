! ***********************************************************************
!
!   Copyright (C) 2014-2022  The MESA Team, Bill Paxton & Matthias Fabry
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


   subroutine grid_plot(id, device_id, array_ix, ierr)
      integer, intent(in) :: id, device_id, array_ix
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call do_Grid_plot(b, id, device_id, &
         b% pg% Grid_xleft(array_ix), b% pg% Grid_xright(array_ix), &
         b% pg% Grid_ybot(array_ix), b% pg% Grid_ytop(array_ix), .false., &
         b% pg% Grid_title(array_ix), &
         b% pg% Grid_txt_scale_factor(array_ix, :), &
         b% pg% Grid_num_cols(array_ix), &
         b% pg% Grid_num_rows(array_ix), &
         b% pg% Grid_num_plots(array_ix), &
         b% pg% Grid_plot_name(array_ix, :), &
         b% pg% Grid_plot_row(array_ix, :), &
         b% pg% Grid_plot_rowspan(array_ix, :), &
         b% pg% Grid_plot_col(array_ix, :), &
         b% pg% Grid_plot_colspan(array_ix, :), &
         b% pg% Grid_plot_pad_left(array_ix, :), &
         b% pg% Grid_plot_pad_right(array_ix, :), &
         b% pg% Grid_plot_pad_top(array_ix, :), &
         b% pg% Grid_plot_pad_bot(array_ix, :), &
         ierr)
   end subroutine grid_plot

   subroutine do_Grid_plot(b, id, device_id, &
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

      use utils_lib, only: StrLowCase
      use pgbinary_summary_history, only: do_summary_history_plot
      use pgbinary_summary, only: do_Text_Summary_plot
      use pgbinary_history_panels, only: do_History_Panels_plot
      use pgbinary_hist_track, only: do_History_Track_plot
      use pgbinary_star, only: do_Star1_plot, do_Star2_plot
      use pgbinary_orbit, only: do_orbit_plot

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
         Grid_plot_pad_left(:), Grid_plot_pad_right(:), &
         Grid_plot_pad_top(:), Grid_plot_pad_bot(:)
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
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% pg% Summary_History_title, &
               Grid_txt_scale_factor(i) * b% pg% Summary_History_txt_scale, ierr)
         case ('history_panels(1)')
            call do_History_Panels_plot(&
               b, id, device_id, 1, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Panels_title(1), &
               Grid_txt_scale_factor(i) * b% pg% History_Panels_txt_scale(1), ierr)
         case ('history_panels(2)')
            call do_History_Panels_plot(&
               b, id, device_id, 2, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Panels_title(2), &
               Grid_txt_scale_factor(i) * b% pg% History_Panels_txt_scale(2), ierr)
         case ('history_panels(3)')
            call do_History_Panels_plot(&
               b, id, device_id, 3, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Panels_title(3), &
               Grid_txt_scale_factor(i) * b% pg% History_Panels_txt_scale(3), ierr)
         case ('history_panels(4)')
            call do_History_Panels_plot(&
               b, id, device_id, 4, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Panels_title(4), &
               Grid_txt_scale_factor(i) * b% pg% History_Panels_txt_scale(4), ierr)
         case ('history_panels(5)')
            call do_History_Panels_plot(&
               b, id, device_id, 5, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Panels_title(5), &
               Grid_txt_scale_factor(i) * b% pg% History_Panels_txt_scale(5), ierr)
         case ('history_panels(6)')
            call do_History_Panels_plot(&
               b, id, device_id, 6, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Panels_title(6), &
               Grid_txt_scale_factor(i) * b% pg% History_Panels_txt_scale(6), ierr)
         case ('history_panels(7)')
            call do_History_Panels_plot(&
               b, id, device_id, 7, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Panels_title(7), &
               Grid_txt_scale_factor(i) * b% pg% History_Panels_txt_scale(7), ierr)
         case ('history_panels(8)')
            call do_History_Panels_plot(&
               b, id, device_id, 8, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Panels_title(8), &
               Grid_txt_scale_factor(i) * b% pg% History_Panels_txt_scale(8), ierr)
         case ('history_panels(9)')
            call do_History_Panels_plot(&
               b, id, device_id, 9, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Panels_title(9), &
               Grid_txt_scale_factor(i) * b% pg% History_Panels_txt_scale(9), ierr)
         case ('history_track(1)')
            call do_History_Track_plot(&
               b, id, device_id, 1, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Track_title(1), &
               Grid_txt_scale_factor(i) * b% pg% History_Track_txt_scale(1), ierr)
         case ('history_track(2)')
            call do_History_Track_plot(&
               b, id, device_id, 2, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Track_title(2), &
               Grid_txt_scale_factor(i) * b% pg% History_Track_txt_scale(2), ierr)
         case ('history_track(3)')
            call do_History_Track_plot(&
               b, id, device_id, 3, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Track_title(3), &
               Grid_txt_scale_factor(i) * b% pg% History_Track_txt_scale(3), ierr)
         case ('history_track(4)')
            call do_History_Track_plot(&
               b, id, device_id, 4, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Track_title(4), &
               Grid_txt_scale_factor(i) * b% pg% History_Track_txt_scale(4), ierr)
         case ('history_track(5)')
            call do_History_Track_plot(&
               b, id, device_id, 5, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Track_title(5), &
               Grid_txt_scale_factor(i) * b% pg% History_Track_txt_scale(5), ierr)
         case ('history_track(6)')
            call do_History_Track_plot(&
               b, id, device_id, 6, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Track_title(6), &
               Grid_txt_scale_factor(i) * b% pg% History_Track_txt_scale(6), ierr)
         case ('history_track(7)')
            call do_History_Track_plot(&
               b, id, device_id, 7, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Track_title(7), &
               Grid_txt_scale_factor(i) * b% pg% History_Track_txt_scale(7), ierr)
         case ('history_track(8)')
            call do_History_Track_plot(&
               b, id, device_id, 8, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Track_title(8), &
               Grid_txt_scale_factor(i) * b% pg% History_Track_txt_scale(8), ierr)
         case ('history_track(9)')
            call do_History_Track_plot(&
               b, id, device_id, 9, xleft, xright, ybot, ytop, grid_subplot, b% pg% History_Track_title(9), &
               Grid_txt_scale_factor(i) * b% pg% History_Track_txt_scale(9), ierr)
         case ('text_summary(1)')
            call do_Text_Summary_plot(&
               b, id, device_id, 1, xleft, xright, ybot, ytop, grid_subplot, b% pg% Text_Summary_title(1), &
               Grid_txt_scale_factor(i) * b% pg% Text_Summary_txt_scale(1), ierr)
         case ('text_summary(2)')
            call do_Text_Summary_plot(&
               b, id, device_id, 2, xleft, xright, ybot, ytop, grid_subplot, b% pg% Text_Summary_title(2), &
               Grid_txt_scale_factor(i) * b% pg% Text_Summary_txt_scale(2), ierr)
         case ('text_summary(3)')
            call do_Text_Summary_plot(&
               b, id, device_id, 3, xleft, xright, ybot, ytop, grid_subplot, b% pg% Text_Summary_title(3), &
               Grid_txt_scale_factor(i) * b% pg% Text_Summary_txt_scale(3), ierr)
         case ('text_summary(4)')
            call do_Text_Summary_plot(&
               b, id, device_id, 4, xleft, xright, ybot, ytop, grid_subplot, b% pg% Text_Summary_title(4), &
               Grid_txt_scale_factor(i) * b% pg% Text_Summary_txt_scale(4), ierr)
         case ('text_summary(5)')
            call do_Text_Summary_plot(&
               b, id, device_id, 5, xleft, xright, ybot, ytop, grid_subplot, b% pg% Text_Summary_title(5), &
               Grid_txt_scale_factor(i) * b% pg% Text_Summary_txt_scale(5), ierr)
         case ('text_summary(6)')
            call do_Text_Summary_plot(&
               b, id, device_id, 6, xleft, xright, ybot, ytop, grid_subplot, b% pg% Text_Summary_title(6), &
               Grid_txt_scale_factor(i) * b% pg% Text_Summary_txt_scale(6), ierr)
         case ('text_summary(7)')
            call do_Text_Summary_plot(&
               b, id, device_id, 7, xleft, xright, ybot, ytop, grid_subplot, b% pg% Text_Summary_title(7), &
               Grid_txt_scale_factor(i) * b% pg% Text_Summary_txt_scale(7), ierr)
         case ('text_summary(8)')
            call do_Text_Summary_plot(&
               b, id, device_id, 8, xleft, xright, ybot, ytop, grid_subplot, b% pg% Text_Summary_title(8), &
               Grid_txt_scale_factor(i) * b% pg% Text_Summary_txt_scale(8), ierr)
         case ('text_summary(9)')
            call do_Text_Summary_plot(&
               b, id, device_id, 9, xleft, xright, ybot, ytop, grid_subplot, b% pg% Text_Summary_title(9), &
               Grid_txt_scale_factor(i) * b% pg% Text_Summary_txt_scale(9), ierr)
         case ('star1')
            call do_Star1_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% pg% Star1_Title, &
               Grid_txt_scale_factor(i) * b% pg% Star1_txt_scale_factor, b% pg% Star1_plot_name, ierr)
         case ('star2')
            call do_Star2_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% pg% Star2_Title, &
               Grid_txt_scale_factor(i) * b% pg% Star2_txt_scale_factor, b% pg% Star2_plot_name, ierr)
         case ('orbit')
            call do_orbit_plot(&
               b, id, device_id, xleft, xright, ybot, ytop, grid_subplot, b% pg% Orbit_title, &
               Grid_txt_scale_factor(i) * b% pg% Orbit_txt_scale_factor, ierr)
         case default
            ! check for "other" plot
            found_it = .false.
            do j = 1, max_num_Other_plots
               plot_id = i_Other + j - 1
               p => b% pg% pgbinary_win_file_ptr(plot_id)
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
                  'Text_Summary1..9', &
                  'History_Panels1..9', &
                  'History_Tracks1..9', &
                  'Star1,2'
               write(*, *)

            end if

         end select

         call pgunsa

         if (ierr /= 0) exit

      end do

      call pgebuf()

   end subroutine do_Grid_plot


end module pgbinary_grid

