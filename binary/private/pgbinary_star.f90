! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team & Matthias Fabry
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
         b% pg% Star1_xleft, b% pg% Star1_xright, &
         b% pg% Star1_ybot, b% pg% Star1_ytop, .false., b% pg% Star1_title, &
         b% pg% Star1_txt_scale_factor, b% pg% Star1_plot_name, ierr)
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
         b% pg% Star2_xleft, b% pg% Star2_xright, &
         b% pg% Star2_ybot, b% pg% Star2_ytop, .false., b% pg% Star2_title, &
         b% pg% Star2_txt_scale_factor, b% pg% Star2_plot_name, ierr)
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
      use pgstar_grid, only : do_grid_plot
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
      use pgstar_summary, only : do_Text_Summary_plot
      use pgstar_profile_panels, only : do_Profile_Panels_plot
      use pgstar_history_panels, only : do_History_Panels_plot
      use pgstar_hist_track, only : do_History_Track_plot
      use pgstar_Color_Magnitude, only : do_Color_Magnitude_plot

      type (binary_info), pointer :: b
      logical, intent(in) :: subplot
      integer, intent(in) :: id, device_id, star_number
      real, intent(in) :: xleft, xright, ybot, ytop, Star_txt_scale_factor
      character (len = *), intent(in) :: Star_title, Star_plot_name
      integer, intent(out) :: ierr

      character (len = strlen) :: title, status, mass
      real, dimension(5) :: xs, ys
      logical, parameter :: star_subplot = .true.

      include 'formats'

      call pgsave
      call pgsvp(xleft, xright, ybot, ytop)  ! set viewport
      if (.not. subplot) then  ! do title stuff
         call show_model_number_pgbinary(b)
         call show_age_pgbinary(b)
      end if
      if (b% pg% show_mtrans_status) then
         if (Star_title /= '') title = trim(Star_title) // ":"
         status = ' Detached'
         if (b% mtransfer_rate /= 0d0) then
            if (b% d_i == star_number) then
               status = ' Donor'
            else
               status = ' Accretor'
            end if
         end if
         title = trim(title) // status
      end if
      call show_grid_title_pgbinary(b, title)
      call pgunsa

      ierr = 0

      select case(star_number)
      case(1)
         if (b% pg% do_star1_box) then
            call pgsvp(xleft + b% pg% Star1_box_pad_left, &
               xright + b% pg% Star1_box_pad_right, &
               ybot + b% pg% Star1_box_pad_bot, ytop + b% pg% Star1_box_pad_top)
            call draw_rect()
            call pgsvp(xleft, xright, ybot, ytop)
         end if
         if (b% point_mass_i /= 1) then
            call read_pgstar_inlist(b% s1, b% job% inlist_names(1), ierr)
            call update_pgstar_data(b% s1, ierr)
            call plot_case(b% s1, b% star_ids(1))
            call update_pgstar_history_file(b% s1, ierr)
         else
            write(mass, '(f3.2)') b% m(1) / Msun
            call pgmtxt('T', -2.0, 0.5, 0.5, 'Star 1 not simulated')
            call pgmtxt('T', -3.0, 0.5, 0.5, 'point mass of ' // trim(adjustl(mass)) // ' M\d\(2281)')
         end if
      case(2)
         if (b% pg% do_star2_box) then
            call pgsvp(xleft + b% pg% Star2_box_pad_left, &
               xright + b% pg% Star2_box_pad_right, &
               ybot + b% pg% Star2_box_pad_bot, ytop + b% pg% Star2_box_pad_top)
            call draw_rect()
            call pgsvp(xleft, xright, ybot, ytop)
         end if
         if (b% point_mass_i /= 2) then

            call read_pgstar_inlist(b% s2, b% job% inlist_names(2), ierr)
            call update_pgstar_data(b% s2, ierr)
            call plot_case(b% s2, b% star_ids(2))
            call update_pgstar_history_file(b% s2, ierr)
         else
            write(mass, '(f3.2)') b% m(2) / Msun
            call pgmtxt('T', -2.0, 0.5, 0.5, 'Star 2 not simulated')
            call pgmtxt('T', -3.0, 0.5, 0.5, 'point mass of ' // trim(adjustl(mass)) // ' M\d\(2281)')
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
               s% pg% Abundance_title, Star_txt_scale_factor * s% pg% Abundance_txt_scale, ierr)
         case ('power')
            call do_power_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% pg% Power_title, Star_txt_scale_factor * s% pg% Power_txt_scale, ierr)
         case ('mixing')
            call do_Mixing_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% pg% Mixing_title, Star_txt_scale_factor * s% pg% Mixing_txt_scale, ierr)
         case ('dynamo')
            call do_Dynamo_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% pg% Dynamo_title, Star_txt_scale_factor * s% pg% Dynamo_txt_scale, ierr)
         case ('trho')
            call do_TRho_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% pg% TRho_title, Star_txt_scale_factor * s% pg% TRho_txt_scale, ierr)
         case ('tmaxrho')
            call do_TmaxRho_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, &
               s% pg% TmaxRho_title, Star_txt_scale_factor * s% pg% TmaxRho_txt_scale, ierr)
         case ('mode_prop')
            call do_mode_propagation_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% Mode_Prop_title, &
               Star_txt_scale_factor * s% pg% Mode_Prop_txt_scale, ierr)
         case ('summary_burn')
            call do_summary_burn_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% Summary_Burn_title, &
               Star_txt_scale_factor * s% pg% Summary_Burn_txt_scale, ierr)
         case ('summary_profile')
            call do_summary_profile_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% Summary_Profile_title, &
               Star_txt_scale_factor * s% pg% Summary_Profile_txt_scale, ierr)
         case ('summary_history')
            call do_summary_history_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% Summary_History_title, &
               Star_txt_scale_factor * s% pg% Summary_History_txt_scale, ierr)
         case ('trho_profile')
            call do_TRho_Profile_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% TRho_Profile_title, &
               Star_txt_scale_factor * s% pg% TRho_Profile_txt_scale, ierr)
         case ('profile_panels(1)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, star_subplot, s% pg% Profile_Panels_title(1), &
               Star_txt_scale_factor * s% pg% Profile_Panels_txt_scale(1), ierr)
         case ('profile_panels(2)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, star_subplot, s% pg% Profile_Panels_title(2), &
               Star_txt_scale_factor * s% pg% Profile_Panels_txt_scale(2), ierr)
         case ('profile_panels(3)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, star_subplot, s% pg% Profile_Panels_title(3), &
               Star_txt_scale_factor * s% pg% Profile_Panels_txt_scale(3), ierr)
         case ('profile_panels(4)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, star_subplot, s% pg% Profile_Panels_title(4), &
               Star_txt_scale_factor * s% pg% Profile_Panels_txt_scale(4), ierr)
         case ('profile_panels(5)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, star_subplot, s% pg% Profile_Panels_title(5), &
               Star_txt_scale_factor * s% pg% Profile_Panels_txt_scale(5), ierr)
         case ('profile_panels(6)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, star_subplot, s% pg% Profile_Panels_title(6), &
               Star_txt_scale_factor * s% pg% Profile_Panels_txt_scale(6), ierr)
         case ('profile_panels(7)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, star_subplot, s% pg% Profile_Panels_title(7), &
               Star_txt_scale_factor * s% pg% Profile_Panels_txt_scale(7), ierr)
         case ('profile_panels(8)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, star_subplot, s% pg% Profile_Panels_title(8), &
               Star_txt_scale_factor * s% pg% Profile_Panels_txt_scale(8), ierr)
         case ('profile_panels(9)')
            call do_Profile_Panels_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, star_subplot, s% pg% Profile_Panels_title(9), &
               Star_txt_scale_factor * s% pg% Profile_Panels_txt_scale(9), ierr)
         case ('logg_teff')
            call do_logg_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% logg_Teff_title, &
               Star_txt_scale_factor * s% pg% logg_Teff_txt_scale, ierr)
         case ('logg_logt')
            call do_logg_logT_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% logg_logT_title, &
               Star_txt_scale_factor * s% pg% logg_logT_txt_scale, ierr)
         case ('hr')
            call do_HR_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% HR_title, &
               Star_txt_scale_factor * s% pg% HR_txt_scale, ierr)
         case ('logl_r')
            call do_logL_R_plot(&
               s, id, device_id, s% pg% show_logL_photosphere_r, xleft, xright, ybot, ytop, &
               star_subplot, s% pg% logL_R_title, &
               Star_txt_scale_factor * s% pg% logL_R_txt_scale, ierr)
         case ('logl_v')
            call do_logL_v_plot(&
               s, id, device_id, s% pg% show_logL_photosphere_v, xleft, xright, ybot, ytop, &
               star_subplot, s% pg% logL_v_title, &
               Star_txt_scale_factor * s% pg% logL_v_txt_scale, ierr)
         case ('logl_teff')
            call do_logL_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% logL_Teff_title, &
               Star_txt_scale_factor * s% pg% logL_Teff_txt_scale, ierr)
         case ('l_r')
            call do_L_R_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% L_R_title, &
               Star_txt_scale_factor * s% pg% L_R_txt_scale, ierr)
         case ('l_v')
            call do_L_v_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% L_v_title, &
               Star_txt_scale_factor * s% pg% L_v_txt_scale, ierr)
         case ('l_teff')
            call do_L_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% L_Teff_title, &
               Star_txt_scale_factor * s% pg% L_Teff_txt_scale, ierr)
         case ('r_l')
            call do_R_L_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% R_L_title, &
               Star_txt_scale_factor * s% pg% R_L_txt_scale, ierr)
         case ('r_teff')
            call do_R_Teff_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% R_Teff_title, &
               Star_txt_scale_factor * s% pg% R_Teff_txt_scale, ierr)
         case ('dpg_dnu')
            call do_dPg_dnu_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% dPg_dnu_title, &
               Star_txt_scale_factor * s% pg% dPg_dnu_txt_scale, ierr)
         case ('history_panels(1)')
            call do_History_Panels_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Panels_title(1), &
               Star_txt_scale_factor * s% pg% History_Panels_txt_scale(1), ierr)
         case ('history_panels(2)')
            call do_History_Panels_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Panels_title(2), &
               Star_txt_scale_factor * s% pg% History_Panels_txt_scale(2), ierr)
         case ('history_panels(3)')
            call do_History_Panels_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Panels_title(3), &
               Star_txt_scale_factor * s% pg% History_Panels_txt_scale(3), ierr)
         case ('history_panels(4)')
            call do_History_Panels_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Panels_title(4), &
               Star_txt_scale_factor * s% pg% History_Panels_txt_scale(4), ierr)
         case ('history_panels(5)')
            call do_History_Panels_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Panels_title(5), &
               Star_txt_scale_factor * s% pg% History_Panels_txt_scale(5), ierr)
         case ('history_panels(6)')
            call do_History_Panels_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Panels_title(6), &
               Star_txt_scale_factor * s% pg% History_Panels_txt_scale(6), ierr)
         case ('history_panels(7)')
            call do_History_Panels_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Panels_title(7), &
               Star_txt_scale_factor * s% pg% History_Panels_txt_scale(7), ierr)
         case ('history_panels(8)')
            call do_History_Panels_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Panels_title(8), &
               Star_txt_scale_factor * s% pg% History_Panels_txt_scale(8), ierr)
         case ('history_panels(9)')
            call do_History_Panels_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Panels_title(9), &
               Star_txt_scale_factor * s% pg% History_Panels_txt_scale(9), ierr)
         case ('color_magnitude(1)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, star_subplot, s% pg% Color_Magnitude_title(1), &
               Star_txt_scale_factor * s% pg% Color_Magnitude_txt_scale(1), ierr)
         case ('color_magnitude(2)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, star_subplot, s% pg% Color_Magnitude_title(2), &
               Star_txt_scale_factor * s% pg% Color_Magnitude_txt_scale(2), ierr)
         case ('color_magnitude(3)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, star_subplot, s% pg% Color_Magnitude_title(3), &
               Star_txt_scale_factor * s% pg% Color_Magnitude_txt_scale(3), ierr)
         case ('color_magnitude(4)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, star_subplot, s% pg% Color_Magnitude_title(4), &
               Star_txt_scale_factor * s% pg% Color_Magnitude_txt_scale(4), ierr)
         case ('color_magnitude(5)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, star_subplot, s% pg% Color_Magnitude_title(5), &
               Star_txt_scale_factor * s% pg% Color_Magnitude_txt_scale(5), ierr)
         case ('color_magnitude(6)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, star_subplot, s% pg% Color_Magnitude_title(6), &
               Star_txt_scale_factor * s% pg% Color_Magnitude_txt_scale(6), ierr)
         case ('color_magnitude(7)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, star_subplot, s% pg% Color_Magnitude_title(7), &
               Star_txt_scale_factor * s% pg% Color_Magnitude_txt_scale(7), ierr)
         case ('color_magnitude(8)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, star_subplot, s% pg% Color_Magnitude_title(8), &
               Star_txt_scale_factor * s% pg% Color_Magnitude_txt_scale(8), ierr)
         case ('color_magnitude(9)')
            call do_Color_Magnitude_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, star_subplot, s% pg% Color_Magnitude_title(9), &
               Star_txt_scale_factor * s% pg% Color_Magnitude_txt_scale(9), ierr)
         case ('history_track(1)')
            call do_History_Track_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Track_title(1), &
               Star_txt_scale_factor * s% pg% History_Track_txt_scale(1), ierr)
         case ('history_track(2)')
            call do_History_Track_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Track_title(2), &
               Star_txt_scale_factor * s% pg% History_Track_txt_scale(2), ierr)
         case ('history_track(3)')
            call do_History_Track_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Track_title(3), &
               Star_txt_scale_factor * s% pg% History_Track_txt_scale(3), ierr)
         case ('history_track(4)')
            call do_History_Track_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Track_title(4), &
               Star_txt_scale_factor * s% pg% History_Track_txt_scale(4), ierr)
         case ('history_track(5)')
            call do_History_Track_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Track_title(5), &
               Star_txt_scale_factor * s% pg% History_Track_txt_scale(5), ierr)
         case ('history_track(6)')
            call do_History_Track_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Track_title(6), &
               Star_txt_scale_factor * s% pg% History_Track_txt_scale(6), ierr)
         case ('history_track(7)')
            call do_History_Track_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Track_title(7), &
               Star_txt_scale_factor * s% pg% History_Track_txt_scale(7), ierr)
         case ('history_track(8)')
            call do_History_Track_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Track_title(8), &
               Star_txt_scale_factor * s% pg% History_Track_txt_scale(8), ierr)
         case ('history_track(9)')
            call do_History_Track_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, star_subplot, s% pg% History_Track_title(9), &
               Star_txt_scale_factor * s% pg% History_Track_txt_scale(9), ierr)
         case ('kipp')
            call do_Kipp_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% Kipp_title, &
               Star_txt_scale_factor * s% pg% Kipp_txt_scale, ierr)
         case ('network')
            call do_Network_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% Network_title, &
               Star_txt_scale_factor * s% pg% Network_txt_scale, ierr)
         case ('production')
            call do_Production_plot(&
               s, id, device_id, xleft, xright, ybot, ytop, star_subplot, s% pg% Production_title, &
               Star_txt_scale_factor * s% pg% Production_txt_scale, ierr)
         case ('text_summary(1)')
            call do_Text_Summary_plot(&
               s, id, device_id, 1, xleft, xright, ybot, ytop, star_subplot, s% pg% Text_Summary_title(1), &
               Star_txt_scale_factor * s% pg% Text_Summary_txt_scale(1), s% pg% Text_Summary_dxval(1), ierr)
         case ('text_summary(2)')
            call do_Text_Summary_plot(&
               s, id, device_id, 2, xleft, xright, ybot, ytop, star_subplot, s% pg% Text_Summary_title(2), &
               Star_txt_scale_factor * s% pg% Text_Summary_txt_scale(2), s% pg% Text_Summary_dxval(2), ierr)
         case ('text_summary(3)')
            call do_Text_Summary_plot(&
               s, id, device_id, 3, xleft, xright, ybot, ytop, star_subplot, s% pg% Text_Summary_title(3), &
               Star_txt_scale_factor * s% pg% Text_Summary_txt_scale(3), s% pg% Text_Summary_dxval(3), ierr)
         case ('text_summary(4)')
            call do_Text_Summary_plot(&
               s, id, device_id, 4, xleft, xright, ybot, ytop, star_subplot, s% pg% Text_Summary_title(4), &
               Star_txt_scale_factor * s% pg% Text_Summary_txt_scale(4), s% pg% Text_Summary_dxval(4), ierr)
         case ('text_summary(5)')
            call do_Text_Summary_plot(&
               s, id, device_id, 5, xleft, xright, ybot, ytop, star_subplot, s% pg% Text_Summary_title(5), &
               Star_txt_scale_factor * s% pg% Text_Summary_txt_scale(5), s% pg% Text_Summary_dxval(5), ierr)
         case ('text_summary(6)')
            call do_Text_Summary_plot(&
               s, id, device_id, 6, xleft, xright, ybot, ytop, star_subplot, s% pg% Text_Summary_title(6), &
               Star_txt_scale_factor * s% pg% Text_Summary_txt_scale(6), s% pg% Text_Summary_dxval(6), ierr)
         case ('text_summary(7)')
            call do_Text_Summary_plot(&
               s, id, device_id, 7, xleft, xright, ybot, ytop, star_subplot, s% pg% Text_Summary_title(7), &
               Star_txt_scale_factor * s% pg% Text_Summary_txt_scale(7), s% pg% Text_Summary_dxval(7), ierr)
         case ('text_summary(8)')
            call do_Text_Summary_plot(&
               s, id, device_id, 8, xleft, xright, ybot, ytop, star_subplot, s% pg% Text_Summary_title(8), &
               Star_txt_scale_factor * s% pg% Text_Summary_txt_scale(8), s% pg% Text_Summary_dxval(8), ierr)
         case ('text_summary(9)')
            call do_Text_Summary_plot(&
               s, id, device_id, 9, xleft, xright, ybot, ytop, star_subplot, s% pg% Text_Summary_title(9), &
               Star_txt_scale_factor * s% pg% Text_Summary_txt_scale(9), s% pg% Text_Summary_dxval(9), ierr)
         case('grid(1)')
            call do_grid_plot(&
               s, star_id, device_id, 1, xleft, xright, ybot, ytop, star_subplot, s% pg% Grid_title(1), &
               Star_txt_scale_factor * s% pg% Grid_txt_scale_factor(1, :), ierr)
         case('grid(2)')
            call do_grid_plot(&
               s, star_id, device_id, 2, xleft, xright, ybot, ytop, star_subplot, s% pg% Grid_title(2), &
               Star_txt_scale_factor * s% pg% Grid_txt_scale_factor(2, :), ierr)
         case('grid(3)')
            call do_grid_plot(&
               s, star_id, device_id, 3, xleft, xright, ybot, ytop, star_subplot, s% pg% Grid_title(3), &
               Star_txt_scale_factor * s% pg% Grid_txt_scale_factor(3, :), ierr)
         case('grid(4)')
            call do_grid_plot(&
               s, star_id, device_id, 4, xleft, xright, ybot, ytop, star_subplot, s% pg% Grid_title(4), &
               Star_txt_scale_factor * s% pg% Grid_txt_scale_factor(4, :), ierr)
         case('grid(5)')
            call do_grid_plot(&
               s, star_id, device_id, 5, xleft, xright, ybot, ytop, star_subplot, s% pg% Grid_title(5), &
               Star_txt_scale_factor * s% pg% Grid_txt_scale_factor(5, :), ierr)
         case('grid(6)')
            call do_grid_plot(&
               s, star_id, device_id, 6, xleft, xright, ybot, ytop, star_subplot, s% pg% Grid_title(6), &
               Star_txt_scale_factor * s% pg% Grid_txt_scale_factor(6, :), ierr)
         case('grid(7)')
            call do_grid_plot(&
               s, star_id, device_id, 7, xleft, xright, ybot, ytop, star_subplot, s% pg% Grid_title(7), &
               Star_txt_scale_factor * s% pg% Grid_txt_scale_factor(7, :), ierr)
         case('grid(8)')
            call do_grid_plot(&
               s, star_id, device_id, 8, xleft, xright, ybot, ytop, star_subplot, s% pg% Grid_title(8), &
               Star_txt_scale_factor * s% pg% Grid_txt_scale_factor(8, :), ierr)
         case('grid(9)')
            call do_grid_plot(&
               s, star_id, device_id, 9, xleft, xright, ybot, ytop, star_subplot, s% pg% Grid_title(9), &
               Star_txt_scale_factor * s% pg% Grid_txt_scale_factor(9, :), ierr)
         case default
            ! check for "other" plot
            found_it = .false.
            do j = 1, max_num_Other_plots
               plot_id = i_Other + j - 1
               p => s% pg% pgstar_win_file_ptr(plot_id)
               if (p% okay_to_call_do_plot_in_grid .and. &
                  StrLowCase(p% name) == StrLowCase(Star_plot_name)) then
                  call p% do_plot_in_grid(&
                     star_id, device_id, xleft, xright, ybot, ytop, &
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
                  'Text_Summary(1,..,9)', &
                  'Profile_Panels(1,..,9)', &
                  'History_Panels(1,..,9)', &
                  'History_Tracks(1,..,9)', &
                  'Color_Magnitude(1,..,9)', &
                  'Grid(1,..,9)', &
                  'and if you are using star/astero', &
                  'Echelle', &
                  'Ratios'
               write(*, *)

            end if

         end select

      end subroutine plot_case

   end subroutine do_star_plot


end module pgbinary_star
