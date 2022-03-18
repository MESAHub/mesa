! ***********************************************************************
!
!   Copyright (C) 2014  The MESA Team
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

      module pgstar_grid

      use star_private_def
      use const_def
      use pgstar_support
      use star_pgstar

      implicit none


      contains


      subroutine grid1_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call Grid_plot(s, id, device_id, &
            s% pg% Grid1_xleft, s% pg% Grid1_xright, &
            s% pg% Grid1_ybot, s% pg% Grid1_ytop, .false., s% pg% Grid1_title, &
            s% pg% Grid1_txt_scale_factor, &
            s% pg% Grid1_num_cols, &
            s% pg% Grid1_num_rows, &
            s% pg% Grid1_num_plots, &
            s% pg% Grid1_plot_name, &
            s% pg% Grid1_plot_row, &
            s% pg% Grid1_plot_rowspan, &
            s% pg% Grid1_plot_col, &
            s% pg% Grid1_plot_colspan, &
            s% pg% Grid1_plot_pad_left, &
            s% pg% Grid1_plot_pad_right, &
            s% pg% Grid1_plot_pad_top, &
            s% pg% Grid1_plot_pad_bot, &
            ierr)
      end subroutine grid1_plot


      subroutine grid2_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call Grid_plot(s, id, device_id, &
            s% pg% Grid2_xleft, s% pg% Grid2_xright, &
            s% pg% Grid2_ybot, s% pg% Grid2_ytop, .false., s% pg% Grid2_title, &
            s% pg% Grid2_txt_scale_factor, &
            s% pg% Grid2_num_cols, &
            s% pg% Grid2_num_rows, &
            s% pg% Grid2_num_plots, &
            s% pg% Grid2_plot_name, &
            s% pg% Grid2_plot_row, &
            s% pg% Grid2_plot_rowspan, &
            s% pg% Grid2_plot_col, &
            s% pg% Grid2_plot_colspan, &
            s% pg% Grid2_plot_pad_left, &
            s% pg% Grid2_plot_pad_right, &
            s% pg% Grid2_plot_pad_top, &
            s% pg% Grid2_plot_pad_bot, &
            ierr)
      end subroutine grid2_plot


      subroutine grid3_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call Grid_plot(s, id, device_id, &
            s% pg% Grid3_xleft, s% pg% Grid3_xright, &
            s% pg% Grid3_ybot, s% pg% Grid3_ytop, .false., s% pg% Grid3_title, &
            s% pg% Grid3_txt_scale_factor, &
            s% pg% Grid3_num_cols, &
            s% pg% Grid3_num_rows, &
            s% pg% Grid3_num_plots, &
            s% pg% Grid3_plot_name, &
            s% pg% Grid3_plot_row, &
            s% pg% Grid3_plot_rowspan, &
            s% pg% Grid3_plot_col, &
            s% pg% Grid3_plot_colspan, &
            s% pg% Grid3_plot_pad_left, &
            s% pg% Grid3_plot_pad_right, &
            s% pg% Grid3_plot_pad_top, &
            s% pg% Grid3_plot_pad_bot, &
            ierr)
      end subroutine grid3_plot


      subroutine grid4_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call Grid_plot(s, id, device_id, &
            s% pg% Grid4_xleft, s% pg% Grid4_xright, &
            s% pg% Grid4_ybot, s% pg% Grid4_ytop, .false., s% pg% Grid4_title, &
            s% pg% Grid4_txt_scale_factor, &
            s% pg% Grid4_num_cols, &
            s% pg% Grid4_num_rows, &
            s% pg% Grid4_num_plots, &
            s% pg% Grid4_plot_name, &
            s% pg% Grid4_plot_row, &
            s% pg% Grid4_plot_rowspan, &
            s% pg% Grid4_plot_col, &
            s% pg% Grid4_plot_colspan, &
            s% pg% Grid4_plot_pad_left, &
            s% pg% Grid4_plot_pad_right, &
            s% pg% Grid4_plot_pad_top, &
            s% pg% Grid4_plot_pad_bot, &
            ierr)
      end subroutine grid4_plot


      subroutine grid5_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call Grid_plot(s, id, device_id, &
            s% pg% Grid5_xleft, s% pg% Grid5_xright, &
            s% pg% Grid5_ybot, s% pg% Grid5_ytop, .false., s% pg% Grid5_title, &
            s% pg% Grid5_txt_scale_factor, &
            s% pg% Grid5_num_cols, &
            s% pg% Grid5_num_rows, &
            s% pg% Grid5_num_plots, &
            s% pg% Grid5_plot_name, &
            s% pg% Grid5_plot_row, &
            s% pg% Grid5_plot_rowspan, &
            s% pg% Grid5_plot_col, &
            s% pg% Grid5_plot_colspan, &
            s% pg% Grid5_plot_pad_left, &
            s% pg% Grid5_plot_pad_right, &
            s% pg% Grid5_plot_pad_top, &
            s% pg% Grid5_plot_pad_bot, &
            ierr)
      end subroutine grid5_plot


      subroutine grid6_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call Grid_plot(s, id, device_id, &
            s% pg% Grid6_xleft, s% pg% Grid6_xright, &
            s% pg% Grid6_ybot, s% pg% Grid6_ytop, .false., s% pg% Grid6_title, &
            s% pg% Grid6_txt_scale_factor, &
            s% pg% Grid6_num_cols, &
            s% pg% Grid6_num_rows, &
            s% pg% Grid6_num_plots, &
            s% pg% Grid6_plot_name, &
            s% pg% Grid6_plot_row, &
            s% pg% Grid6_plot_rowspan, &
            s% pg% Grid6_plot_col, &
            s% pg% Grid6_plot_colspan, &
            s% pg% Grid6_plot_pad_left, &
            s% pg% Grid6_plot_pad_right, &
            s% pg% Grid6_plot_pad_top, &
            s% pg% Grid6_plot_pad_bot, &
            ierr)
      end subroutine grid6_plot


      subroutine grid7_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call Grid_plot(s, id, device_id, &
            s% pg% Grid7_xleft, s% pg% Grid7_xright, &
            s% pg% Grid7_ybot, s% pg% Grid7_ytop, .false., s% pg% Grid7_title, &
            s% pg% Grid7_txt_scale_factor, &
            s% pg% Grid7_num_cols, &
            s% pg% Grid7_num_rows, &
            s% pg% Grid7_num_plots, &
            s% pg% Grid7_plot_name, &
            s% pg% Grid7_plot_row, &
            s% pg% Grid7_plot_rowspan, &
            s% pg% Grid7_plot_col, &
            s% pg% Grid7_plot_colspan, &
            s% pg% Grid7_plot_pad_left, &
            s% pg% Grid7_plot_pad_right, &
            s% pg% Grid7_plot_pad_top, &
            s% pg% Grid7_plot_pad_bot, &
            ierr)
      end subroutine grid7_plot


      subroutine grid8_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call Grid_plot(s, id, device_id, &
            s% pg% Grid8_xleft, s% pg% Grid8_xright, &
            s% pg% Grid8_ybot, s% pg% Grid8_ytop, .false., s% pg% Grid8_title, &
            s% pg% Grid8_txt_scale_factor, &
            s% pg% Grid8_num_cols, &
            s% pg% Grid8_num_rows, &
            s% pg% Grid8_num_plots, &
            s% pg% Grid8_plot_name, &
            s% pg% Grid8_plot_row, &
            s% pg% Grid8_plot_rowspan, &
            s% pg% Grid8_plot_col, &
            s% pg% Grid8_plot_colspan, &
            s% pg% Grid8_plot_pad_left, &
            s% pg% Grid8_plot_pad_right, &
            s% pg% Grid8_plot_pad_top, &
            s% pg% Grid8_plot_pad_bot, &
            ierr)
      end subroutine grid8_plot


      subroutine grid9_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call Grid_plot(s, id, device_id, &
            s% pg% Grid9_xleft, s% pg% Grid9_xright, &
            s% pg% Grid9_ybot, s% pg% Grid9_ytop, .false., s% pg% Grid9_title, &
            s% pg% Grid9_txt_scale_factor, &
            s% pg% Grid9_num_cols, &
            s% pg% Grid9_num_rows, &
            s% pg% Grid9_num_plots, &
            s% pg% Grid9_plot_name, &
            s% pg% Grid9_plot_row, &
            s% pg% Grid9_plot_rowspan, &
            s% pg% Grid9_plot_col, &
            s% pg% Grid9_plot_colspan, &
            s% pg% Grid9_plot_pad_left, &
            s% pg% Grid9_plot_pad_right, &
            s% pg% Grid9_plot_pad_top, &
            s% pg% Grid9_plot_pad_bot, &
            ierr)
      end subroutine grid9_plot


      subroutine Grid_plot(s, id, device_id, &
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
         use pgstar_kipp, only: do_Kipp_Plot
         use pgstar_L_R, only: do_L_R_Plot
         use pgstar_L_v, only: do_L_v_Plot
         use pgstar_L_Teff, only: do_L_Teff_Plot
         use pgstar_logL_R, only: do_logL_R_Plot
         use pgstar_logL_v, only: do_logL_v_Plot
         use pgstar_logL_Teff, only: do_logL_Teff_Plot
         use pgstar_r_L, only: do_R_L_Plot
         use pgstar_r_Teff, only: do_R_Teff_Plot
         use pgstar_logg_Teff, only: do_logg_Teff_Plot
         use pgstar_logg_logT, only: do_logg_logT_Plot
         use pgstar_dPg_dnu, only: do_dPg_dnu_Plot
         use pgstar_hr, only: do_HR_Plot
         use pgstar_trho, only: do_TRho_Plot
         use pgstar_tmaxrho, only: do_TmaxRho_Plot
         use pgstar_dynamo, only: do_Dynamo_plot
         use pgstar_mixing_Ds, only: do_Mixing_plot
         use pgstar_trho_profile, only: do_TRho_Profile_plot
         use pgstar_power, only: do_power_plot
         use pgstar_mode_prop, only: do_mode_propagation_plot
         use pgstar_abundance, only: do_abundance_plot
         use pgstar_summary_burn, only: do_summary_burn_plot
         use pgstar_summary_profile, only: do_summary_profile_plot
         use pgstar_summary_history, only: do_summary_history_plot
         use pgstar_network, only : do_network_plot
         use pgstar_production, only : do_production_plot
         use pgstar_summary, only: &
            do_Text_Summary1_plot, do_Text_Summary2_plot, do_Text_Summary3_plot, &
            do_Text_Summary4_plot, do_Text_Summary5_plot, do_Text_Summary6_plot, &
            do_Text_Summary7_plot, do_Text_Summary8_plot, do_Text_Summary9_plot
         use pgstar_profile_panels, only: &
            do_Profile_Panels1_plot, do_Profile_Panels2_plot, do_Profile_Panels3_plot, &
            do_Profile_Panels4_plot, do_Profile_Panels5_plot, do_Profile_Panels6_plot, &
            do_Profile_Panels7_plot, do_Profile_Panels8_plot, do_Profile_Panels9_plot
         use pgstar_history_panels, only: &
            do_History_Panels1_plot, do_History_Panels2_plot, do_History_Panels3_plot, &
            do_History_Panels4_plot, do_History_Panels5_plot, do_History_Panels6_plot, &
            do_History_Panels7_plot, do_History_Panels8_plot, do_History_Panels9_plot
         use pgstar_hist_track, only: &
            do_History_Track1_plot, do_History_Track2_plot, do_History_Track3_plot, &
            do_History_Track4_plot, do_History_Track5_plot, do_History_Track6_plot, &
            do_History_Track7_plot, do_History_Track8_plot, do_History_Track9_plot
         use pgstar_Color_Magnitude, only: &
            do_Color_Magnitude1_plot, do_Color_Magnitude2_plot, do_Color_Magnitude3_plot, &
            do_Color_Magnitude4_plot, do_Color_Magnitude5_plot, do_Color_Magnitude6_plot, &
            do_Color_Magnitude7_plot, do_Color_Magnitude8_plot, do_Color_Magnitude9_plot

         type (star_info), pointer :: s
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
         character (len=*) :: Grid_title, Grid_plot_name(:)
         integer, intent(out) :: ierr

         integer :: i, j, plot_id
         logical :: found_it
         real :: xleft, xright, ybot, ytop
         real :: row_height, col_width
         type (pgstar_win_file_data), pointer :: p
         logical, parameter :: grid_subplot = .true.

         include 'formats'

         ierr = 0
         if (Grid_num_plots <= 0 .or. &
               Grid_num_cols <= 0 .or. Grid_num_rows <= 0) return

         col_width = (Grid_xright - Grid_xleft)/Grid_num_cols
         row_height = (Grid_ytop - Grid_ybot)/Grid_num_rows

         if (col_width <= 0d0 .or. row_height <= 0d0) then
            ierr = -1
            write(*,1) 'Grid: col_width', col_width
            write(*,1) 'row_height', row_height
            write(*,1) 'Grid_xleft', Grid_xleft
            write(*,1) 'Grid_xright', Grid_xright
            write(*,1) 'Grid_ybot', Grid_ybot
            write(*,1) 'Grid_ytop', Grid_ytop
            write(*,2) 'Grid_num_cols', Grid_num_cols
            write(*,2) 'Grid_num_rows', Grid_num_rows
            return
         end if

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call pgsave
         call pgsvp(Grid_xleft, Grid_xright, Grid_ybot, Grid_ytop)
         if (.not. subplot) then
            call show_model_number_pgstar(s)
            call show_age_pgstar(s)
         end if
         call show_grid_title_pgstar(s, Grid_title)
         call pgunsa

         do i=1,Grid_num_plots

            if (len_trim(Grid_plot_name(i))==0) exit

            xleft = Grid_xleft + col_width*(Grid_plot_col(i)-1)
            xright = xleft + col_width*Grid_plot_colspan(i)

            ytop = Grid_ytop - row_height*(Grid_plot_row(i)-1)
            ybot = ytop - row_height*Grid_plot_rowspan(i)

            xleft = xleft + Grid_plot_pad_left(i)
            xright = xright - Grid_plot_pad_right(i)
            ybot = ybot + Grid_plot_pad_bot(i)
            ytop = ytop - Grid_plot_pad_top(i)

            if (xright <= xleft .or. ytop <= ybot) then
               write(*,2) 'Bad pgstar grid spec', i
               write(*,*) 'xright <= xleft', xright <= xleft
               write(*,*) 'ytop <= ybot', ytop <= ybot
               write(*,2) 'xleft', i, xleft
               write(*,2) 'xright', i, xright
               write(*,2) 'ybot', i, ybot
               write(*,2) 'ytop', i, ytop
               write(*,2) 'Grid_plot_pad_left(i)', i, Grid_plot_pad_left(i)
               write(*,2) 'Grid_plot_pad_right(i)', i, Grid_plot_pad_right(i)
               write(*,2) 'Grid_plot_pad_top(i)', i, Grid_plot_pad_top(i)
               write(*,2) 'Grid_plot_pad_bot(i)', i, Grid_plot_pad_bot(i)
               write(*,2) 'col_width', i, col_width
               write(*,2) 'row_height', i, row_height
               exit
            end if

            call pgsave

            select case(StrLowCase(Grid_plot_name(i)))
            case ('abundance')
               call do_abundance_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
                  s% pg% Abundance_title, Grid_txt_scale_factor(i)*s% pg% Abundance_txt_scale, ierr)
            case ('power')
               call do_power_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
                  s% pg% Power_title, Grid_txt_scale_factor(i)*s% pg% Power_txt_scale, ierr)
            case ('mixing')
               call do_Mixing_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
                  s% pg% Mixing_title, Grid_txt_scale_factor(i)*s% pg% Mixing_txt_scale, ierr)
            case ('dynamo')
               call do_Dynamo_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
                  s% pg% Dynamo_title, Grid_txt_scale_factor(i)*s% pg% Dynamo_txt_scale, ierr)
            case ('trho')
               call do_TRho_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
                  s% pg% TRho_title, Grid_txt_scale_factor(i)*s% pg% TRho_txt_scale, ierr)
            case ('tmaxrho')
               call do_TmaxRho_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
                  s% pg% TmaxRho_title, Grid_txt_scale_factor(i)*s% pg% TmaxRho_txt_scale, ierr)
            case ('mode_prop')
               call do_mode_propagation_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Mode_Prop_title, &
                  Grid_txt_scale_factor(i)*s% pg% Mode_Prop_txt_scale, ierr)
            case ('summary_burn')
               call do_summary_burn_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Summary_Burn_title, &
                  Grid_txt_scale_factor(i)*s% pg% Summary_Burn_txt_scale, ierr)
            case ('summary_profile')
               call do_summary_profile_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Summary_Profile_title, &
                  Grid_txt_scale_factor(i)*s% pg% Summary_Profile_txt_scale, ierr)
            case ('summary_history')
               call do_summary_history_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Summary_History_title, &
                  Grid_txt_scale_factor(i)*s% pg% Summary_History_txt_scale, ierr)
            case ('trho_profile')
               call do_TRho_Profile_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% TRho_Profile_title, &
                  Grid_txt_scale_factor(i)*s% pg% TRho_Profile_txt_scale, ierr)
            case ('profile_panels1')
               call do_Profile_Panels1_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels1_title, &
                  Grid_txt_scale_factor(i)*s% pg% Profile_Panels1_txt_scale, ierr)
            case ('profile_panels2')
               call do_Profile_Panels2_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels2_title, &
                  Grid_txt_scale_factor(i)*s% pg% Profile_Panels2_txt_scale, ierr)
            case ('profile_panels3')
               call do_Profile_Panels3_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels3_title, &
                  Grid_txt_scale_factor(i)*s% pg% Profile_Panels3_txt_scale, ierr)
            case ('profile_panels4')
               call do_Profile_Panels4_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels4_title, &
                  Grid_txt_scale_factor(i)*s% pg% Profile_Panels4_txt_scale, ierr)
            case ('profile_panels5')
               call do_Profile_Panels5_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels5_title, &
                  Grid_txt_scale_factor(i)*s% pg% Profile_Panels5_txt_scale, ierr)
            case ('profile_panels6')
               call do_Profile_Panels6_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels6_title, &
                  Grid_txt_scale_factor(i)*s% pg% Profile_Panels6_txt_scale, ierr)
            case ('profile_panels7')
               call do_Profile_Panels7_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels7_title, &
                  Grid_txt_scale_factor(i)*s% pg% Profile_Panels7_txt_scale, ierr)
            case ('profile_panels8')
               call do_Profile_Panels8_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels8_title, &
                  Grid_txt_scale_factor(i)*s% pg% Profile_Panels8_txt_scale, ierr)
            case ('profile_panels9')
               call do_Profile_Panels9_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels9_title, &
                  Grid_txt_scale_factor(i)*s% pg% Profile_Panels9_txt_scale, ierr)
            case ('logg_teff')
               call do_logg_Teff_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% logg_Teff_title, &
                  Grid_txt_scale_factor(i)*s% pg% logg_Teff_txt_scale, ierr)
            case ('logg_logt')
               call do_logg_logT_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% logg_logT_title, &
                  Grid_txt_scale_factor(i)*s% pg% logg_logT_txt_scale, ierr)
            case ('hr')
               call do_HR_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% HR_title, &
                  Grid_txt_scale_factor(i)*s% pg% HR_txt_scale, ierr)
            case ('logl_r')
               call do_logL_R_plot( &
                  s, id, device_id, s% pg% show_logL_photosphere_r, xleft, xright, ybot, ytop, &
                  grid_subplot, s% pg% logL_R_title, &
                  Grid_txt_scale_factor(i)*s% pg% logL_R_txt_scale, ierr)
            case ('logl_v')
               call do_logL_v_plot( &
                  s, id, device_id, s% pg% show_logL_photosphere_v, xleft, xright, ybot, ytop, &
                  grid_subplot, s% pg% logL_v_title, &
                  Grid_txt_scale_factor(i)*s% pg% logL_v_txt_scale, ierr)
            case ('logl_teff')
               call do_logL_Teff_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% logL_Teff_title, &
                  Grid_txt_scale_factor(i)*s% pg% logL_Teff_txt_scale, ierr)
            case ('l_r')
               call do_L_R_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% L_R_title, &
                  Grid_txt_scale_factor(i)*s% pg% L_R_txt_scale, ierr)
            case ('l_v')
               call do_L_v_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% L_v_title, &
                  Grid_txt_scale_factor(i)*s% pg% L_v_txt_scale, ierr)
            case ('l_teff')
               call do_L_Teff_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% L_Teff_title, &
                  Grid_txt_scale_factor(i)*s% pg% L_Teff_txt_scale, ierr)
            case ('r_l')
               call do_R_L_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% R_L_title, &
                  Grid_txt_scale_factor(i)*s% pg% R_L_txt_scale, ierr)
            case ('r_teff')
               call do_R_Teff_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% R_Teff_title, &
                  Grid_txt_scale_factor(i)*s% pg% R_Teff_txt_scale, ierr)
            case ('dpg_dnu')
               call do_dPg_dnu_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% dPg_dnu_title, &
                  Grid_txt_scale_factor(i)*s% pg% dPg_dnu_txt_scale, ierr)
            case ('history_panels1')
               call do_History_Panels1_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels1_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Panels1_txt_scale, ierr)
            case ('history_panels2')
               call do_History_Panels2_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels2_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Panels2_txt_scale, ierr)
            case ('history_panels3')
               call do_History_Panels3_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels3_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Panels3_txt_scale, ierr)
            case ('history_panels4')
               call do_History_Panels4_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels4_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Panels4_txt_scale, ierr)
            case ('history_panels5')
               call do_History_Panels5_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels5_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Panels5_txt_scale, ierr)
            case ('history_panels6')
               call do_History_Panels6_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels6_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Panels6_txt_scale, ierr)
            case ('history_panels7')
               call do_History_Panels7_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels7_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Panels7_txt_scale, ierr)
            case ('history_panels8')
               call do_History_Panels8_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels8_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Panels8_txt_scale, ierr)
            case ('history_panels9')
               call do_History_Panels9_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels9_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Panels9_txt_scale, ierr)
            case ('color_magnitude1')
               call do_Color_Magnitude1_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude1_title, &
                  Grid_txt_scale_factor(i)*s% pg% Color_Magnitude1_txt_scale, ierr)
            case ('color_magnitude2')
               call do_Color_Magnitude2_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude2_title, &
                  Grid_txt_scale_factor(i)*s% pg% Color_Magnitude2_txt_scale, ierr)
            case ('color_magnitude3')
               call do_Color_Magnitude3_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude3_title, &
                  Grid_txt_scale_factor(i)*s% pg% Color_Magnitude3_txt_scale, ierr)
            case ('color_magnitude4')
               call do_Color_Magnitude4_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude4_title, &
                  Grid_txt_scale_factor(i)*s% pg% Color_Magnitude4_txt_scale, ierr)
            case ('color_magnitude5')
               call do_Color_Magnitude5_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude5_title, &
                  Grid_txt_scale_factor(i)*s% pg% Color_Magnitude5_txt_scale, ierr)
            case ('color_magnitude6')
               call do_Color_Magnitude6_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude6_title, &
                  Grid_txt_scale_factor(i)*s% pg% Color_Magnitude6_txt_scale, ierr)
            case ('color_magnitude7')
               call do_Color_Magnitude7_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude7_title, &
                  Grid_txt_scale_factor(i)*s% pg% Color_Magnitude7_txt_scale, ierr)
            case ('color_magnitude8')
               call do_Color_Magnitude8_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude8_title, &
                  Grid_txt_scale_factor(i)*s% pg% Color_Magnitude8_txt_scale, ierr)
            case ('color_magnitude9')
               call do_Color_Magnitude9_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude9_title, &
                  Grid_txt_scale_factor(i)*s% pg% Color_Magnitude9_txt_scale, ierr)
            case ('history_track1')
               call do_History_Track1_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track1_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Track1_txt_scale, ierr)
            case ('history_track2')
               call do_History_Track2_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track2_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Track2_txt_scale, ierr)
            case ('history_track3')
               call do_History_Track3_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track3_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Track3_txt_scale, ierr)
            case ('history_track4')
               call do_History_Track4_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track4_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Track4_txt_scale, ierr)
            case ('history_track5')
               call do_History_Track5_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track5_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Track5_txt_scale, ierr)
            case ('history_track6')
               call do_History_Track6_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track6_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Track6_txt_scale, ierr)
            case ('history_track7')
               call do_History_Track7_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track7_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Track7_txt_scale, ierr)
            case ('history_track8')
               call do_History_Track8_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track8_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Track8_txt_scale, ierr)
            case ('history_track9')
               call do_History_Track9_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track9_title, &
                  Grid_txt_scale_factor(i)*s% pg% History_Track9_txt_scale, ierr)
            case ('kipp')
               call do_Kipp_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Kipp_title, &
                  Grid_txt_scale_factor(i)*s% pg% Kipp_txt_scale, ierr)
            case ('network')
               call do_Network_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Network_title, &
                  Grid_txt_scale_factor(i)*s% pg% Network_txt_scale, ierr)
            case ('production')
               call do_Production_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Production_title, &
                  Grid_txt_scale_factor(i)*s% pg% Production_txt_scale, ierr)
            case ('text_summary1')
               call do_Text_Summary1_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary1_title, &
                  Grid_txt_scale_factor(i)*s% pg% Text_Summary1_txt_scale, s% pg% Text_Summary1_dxval, ierr)
            case ('text_summary2')
               call do_Text_Summary2_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary2_title, &
                  Grid_txt_scale_factor(i)*s% pg% Text_Summary2_txt_scale, s% pg% Text_Summary2_dxval, ierr)
            case ('text_summary3')
               call do_Text_Summary3_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary3_title, &
                  Grid_txt_scale_factor(i)*s% pg% Text_Summary3_txt_scale, s% pg% Text_Summary3_dxval, ierr)
            case ('text_summary4')
               call do_Text_Summary4_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary4_title, &
                  Grid_txt_scale_factor(i)*s% pg% Text_Summary4_txt_scale, s% pg% Text_Summary4_dxval, ierr)
            case ('text_summary5')
               call do_Text_Summary5_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary5_title, &
                  Grid_txt_scale_factor(i)*s% pg% Text_Summary5_txt_scale, s% pg% Text_Summary5_dxval, ierr)
            case ('text_summary6')
               call do_Text_Summary6_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary6_title, &
                  Grid_txt_scale_factor(i)*s% pg% Text_Summary6_txt_scale, s% pg% Text_Summary6_dxval, ierr)
            case ('text_summary7')
               call do_Text_Summary7_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary7_title, &
                  Grid_txt_scale_factor(i)*s% pg% Text_Summary7_txt_scale, s% pg% Text_Summary7_dxval, ierr)
            case ('text_summary8')
               call do_Text_Summary8_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary8_title, &
                  Grid_txt_scale_factor(i)*s% pg% Text_Summary8_txt_scale, s% pg% Text_Summary8_dxval, ierr)
            case ('text_summary9')
               call do_Text_Summary9_plot( &
                  s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary9_title, &
                  Grid_txt_scale_factor(i)*s% pg% Text_Summary9_txt_scale, s% pg% Text_Summary9_dxval, ierr)
            case default
               ! check for "other" plot
               found_it = .false.
               do j = 1, max_num_Other_plots
                  plot_id = i_Other + j - 1
                  p => s% pg% pgstar_win_file_ptr(plot_id)
                  if (p% okay_to_call_do_plot_in_grid .and. &
                        StrLowCase(p% name) == StrLowCase(Grid_plot_name(i))) then
                     call p% do_plot_in_grid( &
                        id, device_id, xleft, xright, ybot, ytop, &
                        Grid_txt_scale_factor(i), ierr)
                     found_it = .true.
                     exit
                  end if
               end do

               if (.not. found_it) then

                  write(*,*) 'FAILED TO RECOGNIZE NAME FOR GRID PLOT: ' // trim(Grid_plot_name(i))
                  write(*,'(a)') &
                     'here are the valid names', &
                     'Kipp', &
                     'HR', &
                     'TRho', &
                     'R_Teff', &
                     'R_L', &
                     'L_Teff', &
                     'L_R', &
                     'L_v', &
                     'logL_Teff', &
                     'logL_R', &
                     'logL_v', &
                     'logg_Teff', &
                     'logg_logT', &
                     'dPg_dnu', &
                     'TRho_Profile', &
                     'Summary_Burn', &
                     'Summary_Profile', &
                     'Summary_History', &
                     'Abundance', &
                     'Network', &
                     'Production', &
                     'Power', &
                     'Dynamo', &
                     'Mixing', &
                     'Mode_Prop', &
                     'Text_Summary1,..,9', &
                     'Profile_Panels1,..,9', &
                     'History_Panels1,..,9', &
                     'History_Tracks1,..,9', &
                     'Color_Magnitude1,..,9', &
                     'and if you are using star/astero', &
                     'Echelle', &
                     'Ratios'
                  write(*,'(A)')

               end if

            end select

            call pgunsa

            if (ierr /= 0) exit

         end do

         call pgebuf()

      end subroutine Grid_plot


      end module pgstar_grid

