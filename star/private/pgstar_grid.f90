! ***********************************************************************
!
!   Copyright (C) 2014-2022  Matthias Fabry & The MESA Team
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

   use star_def
   use pgstar_support
   use star_pgstar

   implicit none


contains


   subroutine grid_plot(id, device_id, array_ix, ierr)
      integer, intent(in) :: id, device_id, array_ix
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return

      call pgslct(device_id)
      call pgbbuf()
      call pgeras()

      call do_grid_plot(s, id, device_id, array_ix, &
         s% pg% Grid_xleft(array_ix), s% pg% Grid_xright(array_ix), &
         s% pg% Grid_ybot(array_ix), s% pg% Grid_ytop(array_ix), .false., &
         s% pg% Grid_title(array_ix), s% pg% Grid_txt_scale_factor(array_ix, :), &
         ierr)

      call pgebuf()
   end subroutine grid_plot

   subroutine do_grid_plot(s, id, device_id, array_ix, &
      winxmin, winxmax, winymin, winymax, subplot, title, &
      txt_scale_factor, ierr)
      integer, intent(in) :: id, device_id, array_ix
      real, intent(in) :: winxmin, winxmax, winymin, winymax
      real, intent(in), dimension(max_num_pgstar_grid_plots) :: txt_scale_factor
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call do_one_Grid_plot(s, id, device_id, array_ix, &
         winxmin, winxmax, winymin, winymax, subplot, title, &
         txt_scale_factor, &
         s% pg% Grid_num_cols(array_ix), &
         s% pg% Grid_num_rows(array_ix), &
         s% pg% Grid_num_plots(array_ix), &
         s% pg% Grid_plot_name(array_ix, :), &
         s% pg% Grid_plot_row(array_ix, :), &
         s% pg% Grid_plot_rowspan(array_ix, :), &
         s% pg% Grid_plot_col(array_ix, :), &
         s% pg% Grid_plot_colspan(array_ix, :), &
         s% pg% Grid_plot_pad_left(array_ix, :), &
         s% pg% Grid_plot_pad_right(array_ix, :), &
         s% pg% Grid_plot_pad_top(array_ix, :), &
         s% pg% Grid_plot_pad_bot(array_ix, :), &
         ierr)
   end subroutine do_grid_plot


   subroutine do_one_Grid_plot(s, id, device_id, array_ix, &
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
      use pgstar_kipp, only : do_Kipp_Plot
      use pgstar_L_R, only : do_L_R_Plot
      use pgstar_L_v, only : do_L_v_Plot
      use pgstar_L_Teff, only : do_L_Teff_Plot
      use pgstar_logL_R, only : do_logL_R_Plot
      use pgstar_logL_v, only : do_logL_v_Plot
      use pgstar_logL_Teff, only : do_logL_Teff_Plot
      use pgstar_r_L, only : do_R_L_Plot
      use pgstar_r_Teff, only : do_R_Teff_Plot
      use pgstar_logg_Teff, only : do_logg_Teff_Plot
      use pgstar_logg_logT, only : do_logg_logT_Plot
      use pgstar_dPg_dnu, only : do_dPg_dnu_Plot
      use pgstar_hr, only : do_HR_Plot
      use pgstar_trho, only : do_TRho_Plot
      use pgstar_tmaxrho, only : do_TmaxRho_Plot
      use pgstar_dynamo, only : do_Dynamo_plot
      use pgstar_mixing_Ds, only : do_Mixing_plot
      use pgstar_trho_profile, only : do_TRho_Profile_plot
      use pgstar_power, only : do_power_plot
      use pgstar_mode_prop, only : do_mode_propagation_plot
      use pgstar_abundance, only : do_abundance_plot
      use pgstar_summary_burn, only : do_summary_burn_plot
      use pgstar_summary_profile, only : do_summary_profile_plot
      use pgstar_summary_history, only : do_summary_history_plot
      use pgstar_network, only : do_network_plot
      use pgstar_production, only : do_production_plot
      use pgstar_summary, only : do_Text_Summary_plot
      use pgstar_profile_panels, only : do_Profile_Panels_plot
      use pgstar_history_panels, only : do_History_Panels_plot
      use pgstar_hist_track, only : do_History_Track_plot
      use pgstar_Color_Magnitude, only : do_Color_Magnitude_plot

      type (star_info), pointer :: s
      logical, intent(in) :: subplot
      integer, intent(in) :: id, device_id, array_ix, &
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
      type (pgstar_win_file_data), pointer :: p
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

      call pgsave
      call pgsvp(Grid_xleft, Grid_xright, Grid_ybot, Grid_ytop)
      if (.not. subplot) then
         call show_model_number_pgstar(s)
         call show_age_pgstar(s)
      end if
      call show_grid_title_pgstar(s, Grid_title)
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
            write(*, 2) 'Bad pgstar grid spec', i
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
         case ('abundance')
            call do_abundance_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
               s% pg% Abundance_title, Grid_txt_scale_factor(i) * s% pg% Abundance_txt_scale, ierr)
         case ('power')
            call do_power_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
               s% pg% Power_title, Grid_txt_scale_factor(i) * s% pg% Power_txt_scale, ierr)
         case ('mixing')
            call do_Mixing_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
               s% pg% Mixing_title, Grid_txt_scale_factor(i) * s% pg% Mixing_txt_scale, ierr)
         case ('dynamo')
            call do_Dynamo_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
               s% pg% Dynamo_title, Grid_txt_scale_factor(i) * s% pg% Dynamo_txt_scale, ierr)
         case ('trho')
            call do_TRho_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
               s% pg% TRho_title, Grid_txt_scale_factor(i) * s% pg% TRho_txt_scale, ierr)
         case ('tmaxrho')
            call do_TmaxRho_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, &
               s% pg% TmaxRho_title, Grid_txt_scale_factor(i) * s% pg% TmaxRho_txt_scale, ierr)
         case ('mode_prop')
            call do_mode_propagation_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Mode_Prop_title, &
               Grid_txt_scale_factor(i) * s% pg% Mode_Prop_txt_scale, ierr)
         case ('summary_burn')
            call do_summary_burn_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Summary_Burn_title, &
               Grid_txt_scale_factor(i) * s% pg% Summary_Burn_txt_scale, ierr)
         case ('summary_profile')
            call do_summary_profile_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Summary_Profile_title, &
               Grid_txt_scale_factor(i) * s% pg% Summary_Profile_txt_scale, ierr)
         case ('summary_history')
            call do_summary_history_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Summary_History_title, &
               Grid_txt_scale_factor(i) * s% pg% Summary_History_txt_scale, ierr)
         case ('trho_profile')
            call do_TRho_Profile_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% TRho_Profile_title, &
               Grid_txt_scale_factor(i) * s% pg% TRho_Profile_txt_scale, ierr)
         case ('profile_panels(1)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels_title(1), &
               Grid_txt_scale_factor(i) * s% pg% Profile_Panels_txt_scale(1), ierr)
         case ('profile_panels(2)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels_title(2), &
               Grid_txt_scale_factor(i) * s% pg% Profile_Panels_txt_scale(2), ierr)
         case ('profile_panels(3)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels_title(3), &
               Grid_txt_scale_factor(i) * s% pg% Profile_Panels_txt_scale(3), ierr)
         case ('profile_panels(4)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels_title(4), &
               Grid_txt_scale_factor(i) * s% pg% Profile_Panels_txt_scale(4), ierr)
         case ('profile_panels(5)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels_title(5), &
               Grid_txt_scale_factor(i) * s% pg% Profile_Panels_txt_scale(5), ierr)
         case ('profile_panels(6)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels_title(6), &
               Grid_txt_scale_factor(i) * s% pg% Profile_Panels_txt_scale(6), ierr)
         case ('profile_panels(7)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels_title(7), &
               Grid_txt_scale_factor(i) * s% pg% Profile_Panels_txt_scale(7), ierr)
         case ('profile_panels(8)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels_title(8), &
               Grid_txt_scale_factor(i) * s% pg% Profile_Panels_txt_scale(8), ierr)
         case ('profile_panels(9)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, grid_subplot, s% pg% Profile_Panels_title(9), &
               Grid_txt_scale_factor(i) * s% pg% Profile_Panels_txt_scale(9), ierr)
         case ('logg_teff')
            call do_logg_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% logg_Teff_title, &
               Grid_txt_scale_factor(i) * s% pg% logg_Teff_txt_scale, ierr)
         case ('logg_logt')
            call do_logg_logT_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% logg_logT_title, &
               Grid_txt_scale_factor(i) * s% pg% logg_logT_txt_scale, ierr)
         case ('hr')
            call do_HR_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% HR_title, &
               Grid_txt_scale_factor(i) * s% pg% HR_txt_scale, ierr)
         case ('logl_r')
            call do_logL_R_plot(&
               s, id, device_id, s% pg% show_logL_photosphere_r, xleft, xright, ybot, ytop, &
               grid_subplot, s% pg% logL_R_title, &
               Grid_txt_scale_factor(i) * s% pg% logL_R_txt_scale, ierr)
         case ('logl_v')
            call do_logL_v_plot(&
               s, id, device_id, s% pg% show_logL_photosphere_v, xleft, xright, ybot, ytop, &
               grid_subplot, s% pg% logL_v_title, &
               Grid_txt_scale_factor(i) * s% pg% logL_v_txt_scale, ierr)
         case ('logl_teff')
            call do_logL_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% logL_Teff_title, &
               Grid_txt_scale_factor(i) * s% pg% logL_Teff_txt_scale, ierr)
         case ('l_r')
            call do_L_R_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% L_R_title, &
               Grid_txt_scale_factor(i) * s% pg% L_R_txt_scale, ierr)
         case ('l_v')
            call do_L_v_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% L_v_title, &
               Grid_txt_scale_factor(i) * s% pg% L_v_txt_scale, ierr)
         case ('l_teff')
            call do_L_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% L_Teff_title, &
               Grid_txt_scale_factor(i) * s% pg% L_Teff_txt_scale, ierr)
         case ('r_l')
            call do_R_L_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% R_L_title, &
               Grid_txt_scale_factor(i) * s% pg% R_L_txt_scale, ierr)
         case ('r_teff')
            call do_R_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% R_Teff_title, &
               Grid_txt_scale_factor(i) * s% pg% R_Teff_txt_scale, ierr)
         case ('dpg_dnu')
            call do_dPg_dnu_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% dPg_dnu_title, &
               Grid_txt_scale_factor(i) * s% pg% dPg_dnu_txt_scale, ierr)
         case ('history_panels(1)')
            call do_History_Panels_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels_title(1), &
               Grid_txt_scale_factor(i) * s% pg% History_Panels_txt_scale(1), ierr)
         case ('history_panels(2)')
            call do_History_Panels_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels_title(2), &
               Grid_txt_scale_factor(i) * s% pg% History_Panels_txt_scale(2), ierr)
         case ('history_panels(3)')
            call do_History_Panels_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels_title(3), &
               Grid_txt_scale_factor(i) * s% pg% History_Panels_txt_scale(3), ierr)
         case ('history_panels(4)')
            call do_History_Panels_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels_title(4), &
               Grid_txt_scale_factor(i) * s% pg% History_Panels_txt_scale(4), ierr)
         case ('history_panels(5)')
            call do_History_Panels_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels_title(5), &
               Grid_txt_scale_factor(i) * s% pg% History_Panels_txt_scale(5), ierr)
         case ('history_panels(6)')
            call do_History_Panels_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels_title(6), &
               Grid_txt_scale_factor(i) * s% pg% History_Panels_txt_scale(6), ierr)
         case ('history_panels(7)')
            call do_History_Panels_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels_title(7), &
               Grid_txt_scale_factor(i) * s% pg% History_Panels_txt_scale(7), ierr)
         case ('history_panels(8)')
            call do_History_Panels_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels_title(8), &
               Grid_txt_scale_factor(i) * s% pg% History_Panels_txt_scale(8), ierr)
         case ('history_panels(9)')
            call do_History_Panels_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Panels_title(9), &
               Grid_txt_scale_factor(i) * s% pg% History_Panels_txt_scale(9), ierr)
         case ('color_magnitude(1)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude_title(1), &
               Grid_txt_scale_factor(i) * s% pg% Color_Magnitude_txt_scale(1), ierr)
         case ('color_magnitude(2)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude_title(2), &
               Grid_txt_scale_factor(i) * s% pg% Color_Magnitude_txt_scale(2), ierr)
         case ('color_magnitude(3)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude_title(3), &
               Grid_txt_scale_factor(i) * s% pg% Color_Magnitude_txt_scale(3), ierr)
         case ('color_magnitude(4)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude_title(4), &
               Grid_txt_scale_factor(i) * s% pg% Color_Magnitude_txt_scale(4), ierr)
         case ('color_magnitude(5)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude_title(5), &
               Grid_txt_scale_factor(i) * s% pg% Color_Magnitude_txt_scale(5), ierr)
         case ('color_magnitude(6)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude_title(6), &
               Grid_txt_scale_factor(i) * s% pg% Color_Magnitude_txt_scale(6), ierr)
         case ('color_magnitude(7)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude_title(7), &
               Grid_txt_scale_factor(i) * s% pg% Color_Magnitude_txt_scale(7), ierr)
         case ('color_magnitude(8)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude_title(8), &
               Grid_txt_scale_factor(i) * s% pg% Color_Magnitude_txt_scale(8), ierr)
         case ('color_magnitude(9)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, grid_subplot, s% pg% Color_Magnitude_title(9), &
               Grid_txt_scale_factor(i) * s% pg% Color_Magnitude_txt_scale(9), ierr)
         case ('history_track(1)')
            call do_History_Track_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track_title(1), &
               Grid_txt_scale_factor(i) * s% pg% History_Track_txt_scale(1), ierr)
         case ('history_track(2)')
            call do_History_Track_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track_title(2), &
               Grid_txt_scale_factor(i) * s% pg% History_Track_txt_scale(2), ierr)
         case ('history_track(3)')
            call do_History_Track_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track_title(3), &
               Grid_txt_scale_factor(i) * s% pg% History_Track_txt_scale(3), ierr)
         case ('history_track(4)')
            call do_History_Track_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track_title(4), &
               Grid_txt_scale_factor(i) * s% pg% History_Track_txt_scale(4), ierr)
         case ('history_track(5)')
            call do_History_Track_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track_title(5), &
               Grid_txt_scale_factor(i) * s% pg% History_Track_txt_scale(5), ierr)
         case ('history_track(6)')
            call do_History_Track_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track_title(6), &
               Grid_txt_scale_factor(i) * s% pg% History_Track_txt_scale(6), ierr)
         case ('history_track(7)')
            call do_History_Track_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track_title(7), &
               Grid_txt_scale_factor(i) * s% pg% History_Track_txt_scale(7), ierr)
         case ('history_track(8)')
            call do_History_Track_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track_title(8), &
               Grid_txt_scale_factor(i) * s% pg% History_Track_txt_scale(8), ierr)
         case ('history_track(9)')
            call do_History_Track_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, grid_subplot, s% pg% History_Track_title(9), &
               Grid_txt_scale_factor(i) * s% pg% History_Track_txt_scale(9), ierr)
         case ('kipp')
            call do_Kipp_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Kipp_title, &
               Grid_txt_scale_factor(i) * s% pg% Kipp_txt_scale, ierr)
         case ('network')
            call do_Network_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Network_title, &
               Grid_txt_scale_factor(i) * s% pg% Network_txt_scale, ierr)
         case ('production')
            call do_Production_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, grid_subplot, s% pg% Production_title, &
               Grid_txt_scale_factor(i) * s% pg% Production_txt_scale, ierr)
         case ('text_summary(1)')
            call do_Text_Summary_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary_title(1), &
               Grid_txt_scale_factor(i) * s% pg% Text_Summary_txt_scale(1), s% pg% Text_Summary_dxval(1), ierr)
         case ('text_summary(2)')
            call do_Text_Summary_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary_title(2), &
               Grid_txt_scale_factor(i) * s% pg% Text_Summary_txt_scale(2), s% pg% Text_Summary_dxval(2), ierr)
         case ('text_summary(3)')
            call do_Text_Summary_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary_title(3), &
               Grid_txt_scale_factor(i) * s% pg% Text_Summary_txt_scale(3), s% pg% Text_Summary_dxval(3), ierr)
         case ('text_summary(4)')
            call do_Text_Summary_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary_title(4), &
               Grid_txt_scale_factor(i) * s% pg% Text_Summary_txt_scale(4), s% pg% Text_Summary_dxval(4), ierr)
         case ('text_summary(5)')
            call do_Text_Summary_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary_title(5), &
               Grid_txt_scale_factor(i) * s% pg% Text_Summary_txt_scale(5), s% pg% Text_Summary_dxval(5), ierr)
         case ('text_summary(6)')
            call do_Text_Summary_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary_title(6), &
               Grid_txt_scale_factor(i) * s% pg% Text_Summary_txt_scale(6), s% pg% Text_Summary_dxval(6), ierr)
         case ('text_summary(7)')
            call do_Text_Summary_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary_title(7), &
               Grid_txt_scale_factor(i) * s% pg% Text_Summary_txt_scale(7), s% pg% Text_Summary_dxval(7), ierr)
         case ('text_summary(8)')
            call do_Text_Summary_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary_title(8), &
               Grid_txt_scale_factor(i) * s% pg% Text_Summary_txt_scale(8), s% pg% Text_Summary_dxval(8), ierr)
         case ('text_summary(9)')
            call do_Text_Summary_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, grid_subplot, s% pg% Text_Summary_title(9), &
               Grid_txt_scale_factor(i) * s% pg% Text_Summary_txt_scale(9), s% pg% Text_Summary_dxval(9), ierr)
         case default
            ! check for "other" plot
            found_it = .false.
            do j = 1, max_num_Other_plots
               plot_id = i_Other + j - 1
               p => s% pg% pgstar_win_file_ptr(plot_id)
               if (p% okay_to_call_do_plot_in_grid .and. &
                  StrLowCase(p% name) == StrLowCase(Grid_plot_name(i))) then
                  call p% do_plot_in_grid(&
                     id, device_id, xleft, xright, ybot, ytop, &
                     Grid_txt_scale_factor(i), ierr)
                  found_it = .true.
                  exit
               end if
            end do

            if (.not. found_it) then

               write(*, *) 'FAILED TO RECOGNIZE NAME FOR GRID PLOT: ' // trim(Grid_plot_name(i))
               write(*, '(a)') &
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
                  'Text_Summary(1,..,9)', &
                  'Profile_Panels(1,..,9)', &
                  'History_Panels(1,..,9)', &
                  'History_Tracks(1,..,9)', &
                  'Color_Magnitude(1,..,9)', &
                  'and if you are using star/astero', &
                  'Echelle', &
                  'Ratios'
               write(*, '(A)')

            end if

         end select

         call pgunsa

         if (ierr /= 0) exit

      end do

   end subroutine do_one_Grid_plot


end module pgstar_grid
