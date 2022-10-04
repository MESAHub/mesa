! ***********************************************************************
!
!   Copyright (C) 2022  Matthias Fabry
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

module pgbinary_star

   use const_def
   use binary_def
   use pgbinary_support

   implicit none

contains


   subroutine Star1_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star1_plot(b, id, device_id, &
         b% Star1_xleft, b% Star1_xright, &
         b% Star1_ybot, b% Star1_ytop, .false., b% Star1_title, &
         b% Star1_txt_scale_factor, b% Star1_plot_name, ierr)
      call pgebuf()
   end subroutine Star1_plot

   subroutine do_Star1_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, &
      subplot, title, txt_scale, plot_name, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title, plot_name
      integer, intent(out) :: ierr
      write(*, *) 'calling do star plot'
      call do_star_plot(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         plot_name, 1, ierr)
   end subroutine do_Star1_plot

   subroutine Star2_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star2_plot(b, id, device_id, &
         b% Star2_xleft, b% Star2_xright, &
         b% Star2_ybot, b% Star2_ytop, .false., b% Star2_title, &
         b% Star2_txt_scale_factor, b% Star2_plot_name, ierr)
      call pgebuf()
   end subroutine Star2_plot

   subroutine do_Star2_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, &
      subplot, title, txt_scale, plot_name, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title, plot_name
      integer, intent(out) :: ierr
      call do_star_plot(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         plot_name, 2, ierr)
   end subroutine do_Star2_plot

   subroutine do_star_plot(b, id, device_id, xleft, xright, &
      ybot, ytop, subplot, Star_title, Star_txt_scale_factor, &
      Star_plot_name, Star_number, ierr)

      use utils_lib, only : StrLowCase
      use star_lib, only : read_pgstar_inlist
      use pgstar, only : update_pgstar_data, update_pgstar_history_file
      use pgstar_grid, only : do_grid1_plot, do_grid2_plot, do_grid3_plot, &
         do_grid4_plot, do_grid5_plot, do_grid6_plot, do_grid7_plot, &
         do_grid8_plot, do_grid9_plot
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
      use pgstar_summary, only : &
         do_Text_Summary1_plot, do_Text_Summary2_plot, do_Text_Summary3_plot, &
         do_Text_Summary4_plot, do_Text_Summary5_plot, do_Text_Summary6_plot, &
         do_Text_Summary7_plot, do_Text_Summary8_plot, do_Text_Summary9_plot
      use pgstar_profile_panels, only : &
         do_Profile_Panels1_plot, do_Profile_Panels2_plot, do_Profile_Panels3_plot, &
         do_Profile_Panels4_plot, do_Profile_Panels5_plot, do_Profile_Panels6_plot, &
         do_Profile_Panels7_plot, do_Profile_Panels8_plot, do_Profile_Panels9_plot
      use pgstar_history_panels, only : &
         do_History_Panels1_plot, do_History_Panels2_plot, do_History_Panels3_plot, &
         do_History_Panels4_plot, do_History_Panels5_plot, do_History_Panels6_plot, &
         do_History_Panels7_plot, do_History_Panels8_plot, do_History_Panels9_plot
      use pgstar_hist_track, only : &
         do_History_Track1_plot, do_History_Track2_plot, do_History_Track3_plot, &
         do_History_Track4_plot, do_History_Track5_plot, do_History_Track6_plot, &
         do_History_Track7_plot, do_History_Track8_plot, do_History_Track9_plot
      use pgstar_Color_Magnitude, only : &
         do_Color_Magnitude1_plot, do_Color_Magnitude2_plot, do_Color_Magnitude3_plot, &
         do_Color_Magnitude4_plot, do_Color_Magnitude5_plot, do_Color_Magnitude6_plot, &
         do_Color_Magnitude7_plot, do_Color_Magnitude8_plot, do_Color_Magnitude9_plot

      type (binary_info), pointer :: b
      logical, intent(in) :: subplot
      integer, intent(in) :: id, device_id, star_number
      real, intent(in) :: xleft, xright, ybot, ytop, Star_txt_scale_factor
      character (len = *), intent(in) :: Star_title, Star_plot_name
      integer, intent(out) :: ierr
      real, dimension(5) :: xs, ys
      logical, parameter :: star_subplot = .true.

      include 'formats'

      call pgsave
      call pgsvp(xleft, xright, ybot, ytop)  ! set viewport
      if (.not. subplot) then  ! do title stuff
         call show_model_number_pgbinary(b)
         call show_age_pgbinary(b)
      end if
      call show_grid_title_pgbinary(b, Star_title)
      call pgunsa

      ierr = 0

      select case(star_number)
      case(1)
         if (b% have_star_1) then
            if (b% do_star1_box) then
               write(*, *) b% star1_box_pad_top
               call pgsvp(xleft + b% star1_box_pad_left, &
                  xright + b% star1_box_pad_right, &
                  ybot + b% star1_box_pad_bot, ytop + b% star1_box_pad_top)
               call draw_rect()
               call pgsvp(xleft, xright, ybot, ytop)
            end if

            call read_pgstar_inlist(b% s1, b% job% inlist_names(1), ierr)
            call update_pgstar_data(b% s1, ierr)
            call plot_case(b% s1, b% star_ids(1))
            call update_pgstar_history_file(b% s1, ierr)
         end if
      case(2)
         if (b% have_star_2) then
            if (b% do_star2_box) then
               call pgsvp(xleft + b% star2_box_pad_left, &
                  xright + b% star2_box_pad_right, &
                  ybot + b% star2_box_pad_bot, ytop + b% star2_box_pad_top)
               call draw_rect()
               call pgsvp(xleft, xright, ybot, ytop)
            end if

            call read_pgstar_inlist(b% s2, b% job% inlist_names(2), ierr)
            call update_pgstar_data(b% s2, ierr)
            call plot_case(b% s2, b% star_ids(2))
            call update_pgstar_history_file(b% s2, ierr)
         end if
      end select

   contains

      subroutine plot_case(s, star_id)
         type (star_info), pointer, intent(in) :: s
         integer, intent(in) :: star_id
         type (pgstar_win_file_data), pointer :: p
         integer :: plot_id, j
         logical :: found_it

         select case(StrLowCase(Star_plot_name))
         case ('abundance')
            call do_abundance_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% Abundance_title, Star_txt_scale_factor * s% Abundance_txt_scale, ierr)
         case ('power')
            call do_power_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% Power_title, Star_txt_scale_factor * s% Power_txt_scale, ierr)
         case ('mixing')
            call do_Mixing_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% Mixing_title, Star_txt_scale_factor * s% Mixing_txt_scale, ierr)
         case ('dynamo')
            call do_Dynamo_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% Dynamo_title, Star_txt_scale_factor * s% Dynamo_txt_scale, ierr)
         case ('trho')
            call do_TRho_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% TRho_title, Star_txt_scale_factor * s% TRho_txt_scale, ierr)
         case ('mode_prop')
            call do_mode_propagation_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Mode_Prop_title, &
               Star_txt_scale_factor * s% Mode_Prop_txt_scale, ierr)
         case ('summary_burn')
            call do_summary_burn_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Summary_Burn_title, &
               Star_txt_scale_factor * s% Summary_Burn_txt_scale, ierr)
         case ('summary_profile')
            call do_summary_profile_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Summary_Profile_title, &
               Star_txt_scale_factor * s% Summary_Profile_txt_scale, ierr)
         case ('summary_history')
            call do_summary_history_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Summary_History_title, &
               Star_txt_scale_factor * s% Summary_History_txt_scale, ierr)
         case ('trho_profile')
            call do_TRho_Profile_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% TRho_Profile_title, &
               Star_txt_scale_factor * s% TRho_Profile_txt_scale, ierr)
         case ('profile_panels1')
            call do_Profile_Panels1_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Profile_Panels1_title, &
               Star_txt_scale_factor * s% Profile_Panels1_txt_scale, ierr)
         case ('profile_panels2')
            call do_Profile_Panels2_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Profile_Panels2_title, &
               Star_txt_scale_factor * s% Profile_Panels2_txt_scale, ierr)
         case ('profile_panels3')
            call do_Profile_Panels3_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Profile_Panels3_title, &
               Star_txt_scale_factor * s% Profile_Panels3_txt_scale, ierr)
         case ('profile_panels4')
            call do_Profile_Panels4_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Profile_Panels4_title, &
               Star_txt_scale_factor * s% Profile_Panels4_txt_scale, ierr)
         case ('profile_panels5')
            call do_Profile_Panels5_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Profile_Panels5_title, &
               Star_txt_scale_factor * s% Profile_Panels5_txt_scale, ierr)
         case ('profile_panels6')
            call do_Profile_Panels6_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Profile_Panels6_title, &
               Star_txt_scale_factor * s% Profile_Panels6_txt_scale, ierr)
         case ('profile_panels7')
            call do_Profile_Panels7_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Profile_Panels7_title, &
               Star_txt_scale_factor * s% Profile_Panels7_txt_scale, ierr)
         case ('profile_panels8')
            call do_Profile_Panels8_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Profile_Panels8_title, &
               Star_txt_scale_factor * s% Profile_Panels8_txt_scale, ierr)
         case ('profile_panels9')
            call do_Profile_Panels9_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Profile_Panels9_title, &
               Star_txt_scale_factor * s% Profile_Panels9_txt_scale, ierr)
         case ('logg_teff')
            call do_logg_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% logg_Teff_title, &
               Star_txt_scale_factor * s% logg_Teff_txt_scale, ierr)
         case ('logg_logt')
            call do_logg_logT_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% logg_logT_title, &
               Star_txt_scale_factor * s% logg_logT_txt_scale, ierr)
         case ('hr')
            call do_HR_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% HR_title, &
               Star_txt_scale_factor * s% HR_txt_scale, ierr)
         case ('logl_r')
            call do_logL_R_plot(&
               s, id, device_id, s% show_logL_photosphere_r, xleft, xright, ybot, ytop, &
               star_subplot, s% logL_R_title, &
               Star_txt_scale_factor * s% logL_R_txt_scale, ierr)
         case ('logl_v')
            call do_logL_v_plot(&
               s, id, device_id, s% show_logL_photosphere_v, xleft, xright, ybot, ytop, &
               star_subplot, s% logL_v_title, &
               Star_txt_scale_factor * s% logL_v_txt_scale, ierr)
         case ('logl_teff')
            call do_logL_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% logL_Teff_title, &
               Star_txt_scale_factor * s% logL_Teff_txt_scale, ierr)
         case ('l_r')
            call do_L_R_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% L_R_title, &
               Star_txt_scale_factor * s% L_R_txt_scale, ierr)
         case ('l_v')
            call do_L_v_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% L_v_title, &
               Star_txt_scale_factor * s% L_v_txt_scale, ierr)
         case ('l_teff')
            call do_L_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% L_Teff_title, &
               Star_txt_scale_factor * s% L_Teff_txt_scale, ierr)
         case ('r_l')
            call do_R_L_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% R_L_title, &
               Star_txt_scale_factor * s% R_L_txt_scale, ierr)
         case ('r_teff')
            call do_R_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% R_Teff_title, &
               Star_txt_scale_factor * s% R_Teff_txt_scale, ierr)
         case ('dpg_dnu')
            call do_dPg_dnu_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% dPg_dnu_title, &
               Star_txt_scale_factor * s% dPg_dnu_txt_scale, ierr)
         case ('history_panels1')
            call do_History_Panels1_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Panels1_title, &
               Star_txt_scale_factor * s% History_Panels1_txt_scale, ierr)
         case ('history_panels2')
            call do_History_Panels2_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Panels2_title, &
               Star_txt_scale_factor * s% History_Panels2_txt_scale, ierr)
         case ('history_panels3')
            call do_History_Panels3_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Panels3_title, &
               Star_txt_scale_factor * s% History_Panels3_txt_scale, ierr)
         case ('history_panels4')
            call do_History_Panels4_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Panels4_title, &
               Star_txt_scale_factor * s% History_Panels4_txt_scale, ierr)
         case ('history_panels5')
            call do_History_Panels5_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Panels5_title, &
               Star_txt_scale_factor * s% History_Panels5_txt_scale, ierr)
         case ('history_panels6')
            call do_History_Panels6_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Panels6_title, &
               Star_txt_scale_factor * s% History_Panels6_txt_scale, ierr)
         case ('history_panels7')
            call do_History_Panels7_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Panels7_title, &
               Star_txt_scale_factor * s% History_Panels7_txt_scale, ierr)
         case ('history_panels8')
            call do_History_Panels8_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Panels8_title, &
               Star_txt_scale_factor * s% History_Panels8_txt_scale, ierr)
         case ('history_panels9')
            call do_History_Panels9_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Panels9_title, &
               Star_txt_scale_factor * s% History_Panels9_txt_scale, ierr)
         case ('color_magnitude1')
            call do_Color_Magnitude1_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Color_Magnitude1_title, &
               Star_txt_scale_factor * s% Color_Magnitude1_txt_scale, ierr)
         case ('color_magnitude2')
            call do_Color_Magnitude2_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Color_Magnitude2_title, &
               Star_txt_scale_factor * s% Color_Magnitude2_txt_scale, ierr)
         case ('color_magnitude3')
            call do_Color_Magnitude3_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Color_Magnitude3_title, &
               Star_txt_scale_factor * s% Color_Magnitude3_txt_scale, ierr)
         case ('color_magnitude4')
            call do_Color_Magnitude4_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Color_Magnitude4_title, &
               Star_txt_scale_factor * s% Color_Magnitude4_txt_scale, ierr)
         case ('color_magnitude5')
            call do_Color_Magnitude5_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Color_Magnitude5_title, &
               Star_txt_scale_factor * s% Color_Magnitude5_txt_scale, ierr)
         case ('color_magnitude6')
            call do_Color_Magnitude6_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Color_Magnitude6_title, &
               Star_txt_scale_factor * s% Color_Magnitude6_txt_scale, ierr)
         case ('color_magnitude7')
            call do_Color_Magnitude7_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Color_Magnitude7_title, &
               Star_txt_scale_factor * s% Color_Magnitude7_txt_scale, ierr)
         case ('color_magnitude8')
            call do_Color_Magnitude8_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Color_Magnitude8_title, &
               Star_txt_scale_factor * s% Color_Magnitude8_txt_scale, ierr)
         case ('color_magnitude9')
            call do_Color_Magnitude9_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Color_Magnitude9_title, &
               Star_txt_scale_factor * s% Color_Magnitude9_txt_scale, ierr)
         case ('history_track1')
            call do_History_Track1_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Track1_title, &
               Star_txt_scale_factor * s% History_Track1_txt_scale, ierr)
         case ('history_track2')
            call do_History_Track2_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Track2_title, &
               Star_txt_scale_factor * s% History_Track2_txt_scale, ierr)
         case ('history_track3')
            call do_History_Track3_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Track3_title, &
               Star_txt_scale_factor * s% History_Track3_txt_scale, ierr)
         case ('history_track4')
            call do_History_Track4_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Track4_title, &
               Star_txt_scale_factor * s% History_Track4_txt_scale, ierr)
         case ('history_track5')
            call do_History_Track5_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Track5_title, &
               Star_txt_scale_factor * s% History_Track5_txt_scale, ierr)
         case ('history_track6')
            call do_History_Track6_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Track6_title, &
               Star_txt_scale_factor * s% History_Track6_txt_scale, ierr)
         case ('history_track7')
            call do_History_Track7_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Track7_title, &
               Star_txt_scale_factor * s% History_Track7_txt_scale, ierr)
         case ('history_track8')
            call do_History_Track8_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Track8_title, &
               Star_txt_scale_factor * s% History_Track8_txt_scale, ierr)
         case ('history_track9')
            call do_History_Track9_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% History_Track9_title, &
               Star_txt_scale_factor * s% History_Track9_txt_scale, ierr)
         case ('kipp')
            call do_Kipp_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Kipp_title, &
               Star_txt_scale_factor * s% Kipp_txt_scale, ierr)
         case ('network')
            call do_Network_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Network_title, &
               Star_txt_scale_factor * s% Network_txt_scale, ierr)
         case ('production')
            call do_Production_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Production_title, &
               Star_txt_scale_factor * s% Production_txt_scale, ierr)
         case ('text_summary1')
            call do_Text_Summary1_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Text_Summary1_title, &
               Star_txt_scale_factor * s% Text_Summary1_txt_scale, ierr)
         case ('text_summary2')
            call do_Text_Summary2_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Text_Summary2_title, &
               Star_txt_scale_factor * s% Text_Summary2_txt_scale, ierr)
         case ('text_summary3')
            call do_Text_Summary3_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Text_Summary3_title, &
               Star_txt_scale_factor * s% Text_Summary3_txt_scale, ierr)
         case ('text_summary4')
            call do_Text_Summary4_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Text_Summary4_title, &
               Star_txt_scale_factor * s% Text_Summary4_txt_scale, ierr)
         case ('text_summary5')
            call do_Text_Summary5_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Text_Summary5_title, &
               Star_txt_scale_factor * s% Text_Summary5_txt_scale, ierr)
         case ('text_summary6')
            call do_Text_Summary6_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Text_Summary6_title, &
               Star_txt_scale_factor * s% Text_Summary6_txt_scale, ierr)
         case ('text_summary7')
            call do_Text_Summary7_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Text_Summary7_title, &
               Star_txt_scale_factor * s% Text_Summary7_txt_scale, ierr)
         case ('text_summary8')
            call do_Text_Summary8_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Text_Summary8_title, &
               Star_txt_scale_factor * s% Text_Summary8_txt_scale, ierr)
         case ('text_summary9')
            call do_Text_Summary9_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Text_Summary9_title, &
               Star_txt_scale_factor * s% Text_Summary9_txt_scale, ierr)
         case('grid1')
            call do_grid1_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Grid1_title, &
               Star_txt_scale_factor * s% Grid1_txt_scale_factor, ierr)
         case('grid2')
            call do_grid2_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Grid2_title, &
               Star_txt_scale_factor * s% Grid2_txt_scale_factor, ierr)
         case('grid3')
            call do_grid3_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Grid3_title, &
               Star_txt_scale_factor * s% Grid3_txt_scale_factor, ierr)
         case('grid4')
            call do_grid4_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Grid4_title, &
               Star_txt_scale_factor * s% Grid4_txt_scale_factor, ierr)
         case('grid5')
            call do_grid5_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Grid5_title, &
               Star_txt_scale_factor * s% Grid5_txt_scale_factor, ierr)
         case('grid6')
            call do_grid6_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Grid6_title, &
               Star_txt_scale_factor * s% Grid6_txt_scale_factor, ierr)
         case('grid7')
            call do_grid7_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Grid7_title, &
               Star_txt_scale_factor * s% Grid7_txt_scale_factor, ierr)
         case('grid8')
            call do_grid8_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Grid8_title, &
               Star_txt_scale_factor * s% Grid8_txt_scale_factor, ierr)
         case('grid9')
            call do_grid9_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% Grid9_title, &
               Star_txt_scale_factor * s% Grid9_txt_scale_factor, ierr)
         case default
            ! check for "other" plot
            found_it = .false.
            do j = 1, max_num_Other_plots
               plot_id = i_Other + j - 1
               p => s% pgstar_win_file_ptr(plot_id)
               if (p% okay_to_call_do_plot_in_grid .and. &
                  StrLowCase(p% name) == StrLowCase(Star_plot_name)) then
                  call p% do_plot_in_grid(&
                     id, device_id, xleft, xright, ybot, ytop, &
                     Star_txt_scale_factor, ierr)
                  found_it = .true.
                  exit
               end if
            end do

            if (.not. found_it) then

               write(*, *) 'FAILED TO RECOGNIZE NAME FOR STAR PLOT: ' &
                  // trim(Star_plot_name)
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
                  'Text_Summary1,..,9', &
                  'Profile_Panels1,..,9', &
                  'History_Panels1,..,9', &
                  'History_Tracks1,..,9', &
                  'Color_Magnitude1,..,9', &
                  'Grid1,..,9', &
                  'and if you are using star/astero', &
                  'Echelle', &
                  'Ratios'
               write(*, *)

            end if

         end select

      end subroutine plot_case

   end subroutine do_star_plot


end module pgbinary_star
