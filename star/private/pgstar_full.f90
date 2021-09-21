! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

      use star_private_def
      use const_def
      use chem_def, only: category_name
      use rates_def, only: i_rate
      use pgstar_support

      implicit none


      contains


      subroutine do_create_file_name(s, dir, prefix, name)
         use pgstar_support, only: create_file_name
         type (star_info), pointer :: s
         character (len=*), intent(in) :: dir, prefix
         character (len=*), intent(out) :: name
         call create_file_name(s, dir, prefix, name)
      end subroutine do_create_file_name


      subroutine do_write_plot_to_file(s, p, filename, ierr)
         use star_def, only: star_info, pgstar_win_file_data
         use pgstar_support, only: write_plot_to_file
         type (star_info), pointer :: s
         type (pgstar_win_file_data), pointer :: p
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr
         call write_plot_to_file(s, p, filename, ierr)
      end subroutine do_write_plot_to_file


      subroutine do_show_pgstar_annotations( &
            s, show_annotation1, show_annotation2, show_annotation3)
         use pgstar_support, only: show_annotations
         type (star_info), pointer :: s
         logical, intent(in) :: &
            show_annotation1, show_annotation2, show_annotation3
         call show_annotations( &
            s, show_annotation1, show_annotation2, show_annotation3)
      end subroutine do_show_pgstar_annotations


      subroutine do_start_new_run_for_pgstar(s, ierr) ! reset logs
         use utils_lib
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: iounit
         character (len=strlen) :: fname
         logical :: fexist
         ierr = 0
         fname = trim(s% photo_directory) // '/pgstar.dat'
         inquire(file=trim(fname), exist=fexist)
         if (fexist) then
            open(newunit=iounit, file=trim(fname), status='replace', action='write')
            close(iounit)
         end if
         call pgstar_clear(s)
      end subroutine do_start_new_run_for_pgstar


      subroutine do_restart_run_for_pgstar(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         if (.not. s% job% pgstar_flag) return
         call pgstar_clear(s)
         call read_pgstar_data(s,ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in read_pgstar_data'
            ierr = 0
         end if
      end subroutine do_restart_run_for_pgstar


      subroutine do_read_pgstar_controls(s, inlist_fname, ierr)
         use pgstar_ctrls_io, only: read_pgstar
         type (star_info), pointer :: s
         character(*), intent(in) :: inlist_fname
         integer, intent(out) :: ierr
         ierr = 0
         call read_pgstar(s, inlist_fname, ierr)
         if (ierr /= 0) then
            write(*,*) 'PGSTAR failed in reading ' // trim(inlist_fname)
            return
         end if
         if (s% use_other_set_pgstar_controls) &
            call s% other_set_pgstar_controls(s% id)
         call set_win_file_data(s, ierr)
      end subroutine do_read_pgstar_controls


      subroutine set_win_file_data(s, ierr)
         use pgstar_kipp, only: Kipp_Plot
         use pgstar_L_R, only: L_R_Plot
         use pgstar_L_v, only: L_v_Plot
         use pgstar_L_Teff, only: L_Teff_Plot
         use pgstar_logL_R, only: logL_R_Plot
         use pgstar_logL_v, only: logL_v_Plot
         use pgstar_logL_Teff, only: logL_Teff_Plot
         use pgstar_r_L, only: R_L_Plot
         use pgstar_r_Teff, only: R_Teff_Plot
         use pgstar_logg_Teff, only: logg_Teff_Plot
         use pgstar_logg_logT, only: logg_logT_Plot
         use pgstar_dPg_dnu, only: dPg_dnu_Plot
         use pgstar_hr, only: HR_Plot
         use pgstar_trho, only: TRho_Plot
         use pgstar_tmaxrho, only: TmaxRho_Plot
         use pgstar_dynamo, only: Dynamo_plot
         use pgstar_mixing_Ds, only: Mixing_plot
         use pgstar_trho_profile, only: TRho_Profile_plot
         use pgstar_power, only: power_plot
         use pgstar_mode_prop, only: mode_propagation_plot
         use pgstar_abundance, only: abundance_plot
         use pgstar_summary_burn, only: summary_burn_plot
         use pgstar_summary_profile, only: summary_profile_plot
         use pgstar_summary_history, only: summary_history_plot
         use pgstar_grid, only: &
            grid1_plot, grid2_plot, grid3_plot, grid4_plot, &
            grid5_plot, grid6_plot, grid7_plot, grid8_plot, grid9_plot
         use pgstar_summary, only: &
            Text_Summary1_Plot, Text_Summary2_Plot, Text_Summary3_Plot, &
            Text_Summary4_Plot, Text_Summary5_Plot, Text_Summary6_Plot, &
            Text_Summary7_Plot, Text_Summary8_Plot, Text_Summary9_Plot
         use pgstar_profile_panels, only: &
            Profile_Panels1_plot, Profile_Panels2_plot, Profile_Panels3_plot, &
            Profile_Panels4_plot, Profile_Panels5_plot, Profile_Panels6_plot, &
            Profile_Panels7_plot, Profile_Panels8_plot, Profile_Panels9_plot
         use pgstar_history_panels, only: &
            History_Panels1_plot, History_Panels2_plot, History_Panels3_plot, &
            History_Panels4_plot, History_Panels5_plot, History_Panels6_plot, &
            History_Panels7_plot, History_Panels8_plot, History_Panels9_plot
         use pgstar_hist_track, only: &
            History_Track1_plot, History_Track2_plot, History_Track3_plot, &
            History_Track4_plot, History_Track5_plot, History_Track6_plot, &
            History_Track7_plot, History_Track8_plot, History_Track9_plot
         use pgstar_network, only : Network_plot
         use pgstar_production, only : Production_plot
         use pgstar_Color_Magnitude, only: &
            Color_Magnitude1_plot, Color_Magnitude2_plot, Color_Magnitude3_plot, &
            Color_Magnitude4_plot, Color_Magnitude5_plot, Color_Magnitude6_plot, &
            Color_Magnitude7_plot, Color_Magnitude8_plot, Color_Magnitude9_plot
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         type (pgstar_win_file_data), pointer :: p
         integer :: i

         ! store win and file info in records

         p => s% pgstar_win_file_ptr(i_Text_Summary1)
         p% plot => Text_Summary1_Plot
         p% id = i_Text_Summary1
         p% name = 'Text_Summary1'
         p% win_flag = s% Text_Summary1_win_flag
         p% win_width = s% Text_Summary1_win_width
         p% win_aspect_ratio = s% Text_Summary1_win_aspect_ratio
         p% file_flag = s% Text_Summary1_file_flag
         p% file_dir = s% Text_Summary1_file_dir
         p% file_prefix = s% Text_Summary1_file_prefix
         p% file_interval = s% Text_Summary1_file_interval
         p% file_width = s% Text_Summary1_file_width
         p% file_aspect_ratio = s% Text_Summary1_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Text_Summary2)
         p% plot => Text_Summary2_Plot
         p% id = i_Text_Summary2
         p% name = 'Text_Summary2'
         p% win_flag = s% Text_Summary2_win_flag
         p% win_width = s% Text_Summary2_win_width
         p% win_aspect_ratio = s% Text_Summary2_win_aspect_ratio
         p% file_flag = s% Text_Summary2_file_flag
         p% file_dir = s% Text_Summary2_file_dir
         p% file_prefix = s% Text_Summary2_file_prefix
         p% file_interval = s% Text_Summary2_file_interval
         p% file_width = s% Text_Summary2_file_width
         p% file_aspect_ratio = s% Text_Summary2_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Text_Summary3)
         p% plot => Text_Summary3_Plot
         p% id = i_Text_Summary3
         p% name = 'Text_Summary3'
         p% win_flag = s% Text_Summary3_win_flag
         p% win_width = s% Text_Summary3_win_width
         p% win_aspect_ratio = s% Text_Summary3_win_aspect_ratio
         p% file_flag = s% Text_Summary3_file_flag
         p% file_dir = s% Text_Summary3_file_dir
         p% file_prefix = s% Text_Summary3_file_prefix
         p% file_interval = s% Text_Summary3_file_interval
         p% file_width = s% Text_Summary3_file_width
         p% file_aspect_ratio = s% Text_Summary3_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Text_Summary4)
         p% plot => Text_Summary4_Plot
         p% id = i_Text_Summary4
         p% name = 'Text_Summary4'
         p% win_flag = s% Text_Summary4_win_flag
         p% win_width = s% Text_Summary4_win_width
         p% win_aspect_ratio = s% Text_Summary4_win_aspect_ratio
         p% file_flag = s% Text_Summary4_file_flag
         p% file_dir = s% Text_Summary4_file_dir
         p% file_prefix = s% Text_Summary4_file_prefix
         p% file_interval = s% Text_Summary4_file_interval
         p% file_width = s% Text_Summary4_file_width
         p% file_aspect_ratio = s% Text_Summary4_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Text_Summary5)
         p% plot => Text_Summary5_Plot
         p% id = i_Text_Summary5
         p% name = 'Text_Summary5'
         p% win_flag = s% Text_Summary5_win_flag
         p% win_width = s% Text_Summary5_win_width
         p% win_aspect_ratio = s% Text_Summary5_win_aspect_ratio
         p% file_flag = s% Text_Summary5_file_flag
         p% file_dir = s% Text_Summary5_file_dir
         p% file_prefix = s% Text_Summary5_file_prefix
         p% file_interval = s% Text_Summary5_file_interval
         p% file_width = s% Text_Summary5_file_width
         p% file_aspect_ratio = s% Text_Summary5_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Text_Summary6)
         p% plot => Text_Summary6_Plot
         p% id = i_Text_Summary6
         p% name = 'Text_Summary6'
         p% win_flag = s% Text_Summary6_win_flag
         p% win_width = s% Text_Summary6_win_width
         p% win_aspect_ratio = s% Text_Summary6_win_aspect_ratio
         p% file_flag = s% Text_Summary6_file_flag
         p% file_dir = s% Text_Summary6_file_dir
         p% file_prefix = s% Text_Summary6_file_prefix
         p% file_interval = s% Text_Summary6_file_interval
         p% file_width = s% Text_Summary6_file_width
         p% file_aspect_ratio = s% Text_Summary6_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Text_Summary7)
         p% plot => Text_Summary7_Plot
         p% id = i_Text_Summary7
         p% name = 'Text_Summary7'
         p% win_flag = s% Text_Summary7_win_flag
         p% win_width = s% Text_Summary7_win_width
         p% win_aspect_ratio = s% Text_Summary7_win_aspect_ratio
         p% file_flag = s% Text_Summary7_file_flag
         p% file_dir = s% Text_Summary7_file_dir
         p% file_prefix = s% Text_Summary7_file_prefix
         p% file_interval = s% Text_Summary7_file_interval
         p% file_width = s% Text_Summary7_file_width
         p% file_aspect_ratio = s% Text_Summary7_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Text_Summary8)
         p% plot => Text_Summary8_Plot
         p% id = i_Text_Summary8
         p% name = 'Text_Summary8'
         p% win_flag = s% Text_Summary8_win_flag
         p% win_width = s% Text_Summary8_win_width
         p% win_aspect_ratio = s% Text_Summary8_win_aspect_ratio
         p% file_flag = s% Text_Summary8_file_flag
         p% file_dir = s% Text_Summary8_file_dir
         p% file_prefix = s% Text_Summary8_file_prefix
         p% file_interval = s% Text_Summary8_file_interval
         p% file_width = s% Text_Summary8_file_width
         p% file_aspect_ratio = s% Text_Summary8_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Text_Summary9)
         p% plot => Text_Summary9_Plot
         p% id = i_Text_Summary9
         p% name = 'Text_Summary9'
         p% win_flag = s% Text_Summary9_win_flag
         p% win_width = s% Text_Summary9_win_width
         p% win_aspect_ratio = s% Text_Summary9_win_aspect_ratio
         p% file_flag = s% Text_Summary9_file_flag
         p% file_dir = s% Text_Summary9_file_dir
         p% file_prefix = s% Text_Summary9_file_prefix
         p% file_interval = s% Text_Summary9_file_interval
         p% file_width = s% Text_Summary9_file_width
         p% file_aspect_ratio = s% Text_Summary9_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_TRho_Profile)
         p% plot => TRho_Profile_plot
         p% id = i_TRho_Profile
         p% name = 'TRho_Profile'
         p% win_flag = s% TRho_Profile_win_flag
         p% win_width = s% TRho_Profile_win_width
         p% win_aspect_ratio = s% TRho_Profile_win_aspect_ratio
         p% file_flag = s% TRho_Profile_file_flag
         p% file_dir = s% TRho_Profile_file_dir
         p% file_prefix = s% TRho_Profile_file_prefix
         p% file_interval = s% TRho_Profile_file_interval
         p% file_width = s% TRho_Profile_file_width
         p% file_aspect_ratio = s% TRho_Profile_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Profile_Panels1)
         p% plot => Profile_Panels1_plot
         p% id = i_Profile_Panels1
         p% name = 'Profile_Panels1'
         p% win_flag = s% Profile_Panels1_win_flag
         p% win_width = s% Profile_Panels1_win_width
         p% win_aspect_ratio = s% Profile_Panels1_win_aspect_ratio
         p% file_flag = s% Profile_Panels1_file_flag
         p% file_dir = s% Profile_Panels1_file_dir
         p% file_prefix = s% Profile_Panels1_file_prefix
         p% file_interval = s% Profile_Panels1_file_interval
         p% file_width = s% Profile_Panels1_file_width
         p% file_aspect_ratio = s% Profile_Panels1_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Profile_Panels2)
         p% plot => Profile_Panels2_plot
         p% id = i_Profile_Panels2
         p% name = 'Profile_Panels2'
         p% win_flag = s% Profile_Panels2_win_flag
         p% win_width = s% Profile_Panels2_win_width
         p% win_aspect_ratio = s% Profile_Panels2_win_aspect_ratio
         p% file_flag = s% Profile_Panels2_file_flag
         p% file_dir = s% Profile_Panels2_file_dir
         p% file_prefix = s% Profile_Panels2_file_prefix
         p% file_interval = s% Profile_Panels2_file_interval
         p% file_width = s% Profile_Panels2_file_width
         p% file_aspect_ratio = s% Profile_Panels2_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Profile_Panels3)
         p% plot => Profile_Panels3_plot
         p% id = i_Profile_Panels3
         p% name = 'Profile_Panels3'
         p% win_flag = s% Profile_Panels3_win_flag
         p% win_width = s% Profile_Panels3_win_width
         p% win_aspect_ratio = s% Profile_Panels3_win_aspect_ratio
         p% file_flag = s% Profile_Panels3_file_flag
         p% file_dir = s% Profile_Panels3_file_dir
         p% file_prefix = s% Profile_Panels3_file_prefix
         p% file_interval = s% Profile_Panels3_file_interval
         p% file_width = s% Profile_Panels3_file_width
         p% file_aspect_ratio = s% Profile_Panels3_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Profile_Panels4)
         p% plot => Profile_Panels4_plot
         p% id = i_Profile_Panels4
         p% name = 'Profile_Panels4'
         p% win_flag = s% Profile_Panels4_win_flag
         p% win_width = s% Profile_Panels4_win_width
         p% win_aspect_ratio = s% Profile_Panels4_win_aspect_ratio
         p% file_flag = s% Profile_Panels4_file_flag
         p% file_dir = s% Profile_Panels4_file_dir
         p% file_prefix = s% Profile_Panels4_file_prefix
         p% file_interval = s% Profile_Panels4_file_interval
         p% file_width = s% Profile_Panels4_file_width
         p% file_aspect_ratio = s% Profile_Panels4_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Profile_Panels5)
         p% plot => Profile_Panels5_plot
         p% id = i_Profile_Panels5
         p% name = 'Profile_Panels5'
         p% win_flag = s% Profile_Panels5_win_flag
         p% win_width = s% Profile_Panels5_win_width
         p% win_aspect_ratio = s% Profile_Panels5_win_aspect_ratio
         p% file_flag = s% Profile_Panels5_file_flag
         p% file_dir = s% Profile_Panels5_file_dir
         p% file_prefix = s% Profile_Panels5_file_prefix
         p% file_interval = s% Profile_Panels5_file_interval
         p% file_width = s% Profile_Panels5_file_width
         p% file_aspect_ratio = s% Profile_Panels5_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Profile_Panels6)
         p% plot => Profile_Panels6_plot
         p% id = i_Profile_Panels6
         p% name = 'Profile_Panels6'
         p% win_flag = s% Profile_Panels6_win_flag
         p% win_width = s% Profile_Panels6_win_width
         p% win_aspect_ratio = s% Profile_Panels6_win_aspect_ratio
         p% file_flag = s% Profile_Panels6_file_flag
         p% file_dir = s% Profile_Panels6_file_dir
         p% file_prefix = s% Profile_Panels6_file_prefix
         p% file_interval = s% Profile_Panels6_file_interval
         p% file_width = s% Profile_Panels6_file_width
         p% file_aspect_ratio = s% Profile_Panels6_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Profile_Panels7)
         p% plot => Profile_Panels7_plot
         p% id = i_Profile_Panels7
         p% name = 'Profile_Panels7'
         p% win_flag = s% Profile_Panels7_win_flag
         p% win_width = s% Profile_Panels7_win_width
         p% win_aspect_ratio = s% Profile_Panels7_win_aspect_ratio
         p% file_flag = s% Profile_Panels7_file_flag
         p% file_dir = s% Profile_Panels7_file_dir
         p% file_prefix = s% Profile_Panels7_file_prefix
         p% file_interval = s% Profile_Panels7_file_interval
         p% file_width = s% Profile_Panels7_file_width
         p% file_aspect_ratio = s% Profile_Panels7_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Profile_Panels8)
         p% plot => Profile_Panels8_plot
         p% id = i_Profile_Panels8
         p% name = 'Profile_Panels8'
         p% win_flag = s% Profile_Panels8_win_flag
         p% win_width = s% Profile_Panels8_win_width
         p% win_aspect_ratio = s% Profile_Panels8_win_aspect_ratio
         p% file_flag = s% Profile_Panels8_file_flag
         p% file_dir = s% Profile_Panels8_file_dir
         p% file_prefix = s% Profile_Panels8_file_prefix
         p% file_interval = s% Profile_Panels8_file_interval
         p% file_width = s% Profile_Panels8_file_width
         p% file_aspect_ratio = s% Profile_Panels8_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Profile_Panels9)
         p% plot => Profile_Panels9_plot
         p% id = i_Profile_Panels9
         p% name = 'Profile_Panels9'
         p% win_flag = s% Profile_Panels9_win_flag
         p% win_width = s% Profile_Panels9_win_width
         p% win_aspect_ratio = s% Profile_Panels9_win_aspect_ratio
         p% file_flag = s% Profile_Panels9_file_flag
         p% file_dir = s% Profile_Panels9_file_dir
         p% file_prefix = s% Profile_Panels9_file_prefix
         p% file_interval = s% Profile_Panels9_file_interval
         p% file_width = s% Profile_Panels9_file_width
         p% file_aspect_ratio = s% Profile_Panels9_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_logg_Teff)
         p% plot => logg_Teff_Plot
         p% id = i_logg_Teff
         p% name = 'logg_Teff'
         p% win_flag = s% logg_Teff_win_flag
         p% win_width = s% logg_Teff_win_width
         p% win_aspect_ratio = s% logg_Teff_win_aspect_ratio
         p% file_flag = s% logg_Teff_file_flag
         p% file_dir = s% logg_Teff_file_dir
         p% file_prefix = s% logg_Teff_file_prefix
         p% file_interval = s% logg_Teff_file_interval
         p% file_width = s% logg_Teff_file_width
         p% file_aspect_ratio = s% logg_Teff_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_logL_Teff)
         p% plot => logL_Teff_Plot
         p% id = i_logL_Teff
         p% name = 'logL_Teff'
         p% win_flag = s% logL_Teff_win_flag
         p% win_width = s% logL_Teff_win_width
         p% win_aspect_ratio = s% logL_Teff_win_aspect_ratio
         p% file_flag = s% logL_Teff_file_flag
         p% file_dir = s% logL_Teff_file_dir
         p% file_prefix = s% logL_Teff_file_prefix
         p% file_interval = s% logL_Teff_file_interval
         p% file_width = s% logL_Teff_file_width
         p% file_aspect_ratio = s% logL_Teff_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_L_Teff)
         p% plot => L_Teff_Plot
         p% id = i_L_Teff
         p% name = 'L_Teff'
         p% win_flag = s% L_Teff_win_flag
         p% win_width = s% L_Teff_win_width
         p% win_aspect_ratio = s% L_Teff_win_aspect_ratio
         p% file_flag = s% L_Teff_file_flag
         p% file_dir = s% L_Teff_file_dir
         p% file_prefix = s% L_Teff_file_prefix
         p% file_interval = s% L_Teff_file_interval
         p% file_width = s% L_Teff_file_width
         p% file_aspect_ratio = s% L_Teff_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_L_v)
         p% plot => L_v_Plot
         p% id = i_L_v
         p% name = 'L_v'
         p% win_flag = s% L_v_win_flag
         p% win_width = s% L_v_win_width
         p% win_aspect_ratio = s% L_v_win_aspect_ratio
         p% file_flag = s% L_v_file_flag
         p% file_dir = s% L_v_file_dir
         p% file_prefix = s% L_v_file_prefix
         p% file_interval = s% L_v_file_interval
         p% file_width = s% L_v_file_width
         p% file_aspect_ratio = s% L_v_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_L_R)
         p% plot => L_R_Plot
         p% id = i_L_R
         p% name = 'L_R'
         p% win_flag = s% L_R_win_flag
         p% win_width = s% L_R_win_width
         p% win_aspect_ratio = s% L_R_win_aspect_ratio
         p% file_flag = s% L_R_file_flag
         p% file_dir = s% L_R_file_dir
         p% file_prefix = s% L_R_file_prefix
         p% file_interval = s% L_R_file_interval
         p% file_width = s% L_R_file_width
         p% file_aspect_ratio = s% L_R_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_R_Teff)
         p% plot => R_Teff_Plot
         p% id = i_R_Teff
         p% name = 'R_Teff'
         p% win_flag = s% R_Teff_win_flag
         p% win_width = s% R_Teff_win_width
         p% win_aspect_ratio = s% R_Teff_win_aspect_ratio
         p% file_flag = s% R_Teff_file_flag
         p% file_dir = s% R_Teff_file_dir
         p% file_prefix = s% R_Teff_file_prefix
         p% file_interval = s% R_Teff_file_interval
         p% file_width = s% R_Teff_file_width
         p% file_aspect_ratio = s% R_Teff_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_R_L)
         p% plot => R_L_Plot
         p% id = i_R_L
         p% name = 'R_L'
         p% win_flag = s% R_L_win_flag
         p% win_width = s% R_L_win_width
         p% win_aspect_ratio = s% R_L_win_aspect_ratio
         p% file_flag = s% R_L_file_flag
         p% file_dir = s% R_L_file_dir
         p% file_prefix = s% R_L_file_prefix
         p% file_interval = s% R_L_file_interval
         p% file_width = s% R_L_file_width
         p% file_aspect_ratio = s% R_L_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_logg_logT)
         p% plot => logg_logT_Plot
         p% id = i_logg_logT
         p% name = 'logg_logT'
         p% win_flag = s% logg_logT_win_flag
         p% win_width = s% logg_logT_win_width
         p% win_aspect_ratio = s% logg_logT_win_aspect_ratio
         p% file_flag = s% logg_logT_file_flag
         p% file_dir = s% logg_logT_file_dir
         p% file_prefix = s% logg_logT_file_prefix
         p% file_interval = s% logg_logT_file_interval
         p% file_width = s% logg_logT_file_width
         p% file_aspect_ratio = s% logg_logT_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_dPg_dnu)
         p% plot => dPg_dnu_Plot
         p% id = i_dPg_dnu
         p% name = 'dPg_dnu'
         p% win_flag = s% dPg_dnu_win_flag
         p% win_width = s% dPg_dnu_win_width
         p% win_aspect_ratio = s% dPg_dnu_win_aspect_ratio
         p% file_flag = s% dPg_dnu_file_flag
         p% file_dir = s% dPg_dnu_file_dir
         p% file_prefix = s% dPg_dnu_file_prefix
         p% file_interval = s% dPg_dnu_file_interval
         p% file_width = s% dPg_dnu_file_width
         p% file_aspect_ratio = s% dPg_dnu_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_HR)
         p% plot => HR_Plot
         p% id = i_HR
         p% name = 'HR'
         p% win_flag = s% HR_win_flag
         p% win_width = s% HR_win_width
         p% win_aspect_ratio = s% HR_win_aspect_ratio
         p% file_flag = s% HR_file_flag
         p% file_dir = s% HR_file_dir
         p% file_prefix = s% HR_file_prefix
         p% file_interval = s% HR_file_interval
         p% file_width = s% HR_file_width
         p% file_aspect_ratio = s% HR_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_TRho)
         p% plot => TRho_Plot
         p% id = i_TRho
         p% name = 'TRho'
         p% win_flag = s% TRho_win_flag
         p% win_width = s% TRho_win_width
         p% win_aspect_ratio = s% TRho_win_aspect_ratio
         p% file_flag = s% TRho_file_flag
         p% file_dir = s% TRho_file_dir
         p% file_prefix = s% TRho_file_prefix
         p% file_interval = s% TRho_file_interval
         p% file_width = s% TRho_file_width
         p% file_aspect_ratio = s% TRho_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_TmaxRho)
         p% plot => TmaxRho_Plot
         p% id = i_TmaxRho
         p% name = 'TmaxRho'
         p% win_flag = s% TmaxRho_win_flag
         p% win_width = s% TmaxRho_win_width
         p% win_aspect_ratio = s% TmaxRho_win_aspect_ratio
         p% file_flag = s% TmaxRho_file_flag
         p% file_dir = s% TmaxRho_file_dir
         p% file_prefix = s% TmaxRho_file_prefix
         p% file_interval = s% TmaxRho_file_interval
         p% file_width = s% TmaxRho_file_width
         p% file_aspect_ratio = s% TmaxRho_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Dynamo)
         p% plot => Dynamo_plot
         p% id = i_Dynamo
         p% name = 'Dynamo'
         p% win_flag = s% Dynamo_win_flag
         p% win_width = s% Dynamo_win_width
         p% win_aspect_ratio = s% Dynamo_win_aspect_ratio
         p% file_flag = s% Dynamo_file_flag
         p% file_dir = s% Dynamo_file_dir
         p% file_prefix = s% Dynamo_file_prefix
         p% file_interval = s% Dynamo_file_interval
         p% file_width = s% Dynamo_file_width
         p% file_aspect_ratio = s% Dynamo_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Mixing)
         p% plot => Mixing_plot
         p% id = i_Mixing
         p% name = 'Mixing'
         p% win_flag = s% Mixing_win_flag
         p% win_width = s% Mixing_win_width
         p% win_aspect_ratio = s% Mixing_win_aspect_ratio
         p% file_flag = s% Mixing_file_flag
         p% file_dir = s% Mixing_file_dir
         p% file_prefix = s% Mixing_file_prefix
         p% file_interval = s% Mixing_file_interval
         p% file_width = s% Mixing_file_width
         p% file_aspect_ratio = s% Mixing_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Kipp)
         p% plot => Kipp_Plot
         p% id = i_Kipp
         p% name = 'Kipp'
         p% win_flag = s% Kipp_win_flag
         p% win_width = s% Kipp_win_width
         p% win_aspect_ratio = s% Kipp_win_aspect_ratio
         p% file_flag = s% Kipp_file_flag
         p% file_dir = s% Kipp_file_dir
         p% file_prefix = s% Kipp_file_prefix
         p% file_interval = s% Kipp_file_interval
         p% file_width = s% Kipp_file_width
         p% file_aspect_ratio = s% Kipp_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Network)
         p% plot => Network_Plot
         p% id = i_Network
         p% name = 'Network'
         p% win_flag = s% Network_win_flag
         p% win_width = s% Network_win_width
         p% win_aspect_ratio = s% Network_win_aspect_ratio
         p% file_flag = s% Network_file_flag
         p% file_dir = s% Network_file_dir
         p% file_prefix = s% Network_file_prefix
         p% file_interval = s% Network_file_interval
         p% file_width = s% Network_file_width
         p% file_aspect_ratio = s% Network_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Production)
         p% plot => Production_Plot
         p% id = i_Production
         p% name = 'Production'
         p% win_flag = s% Production_win_flag
         p% win_width = s% Production_win_width
         p% win_aspect_ratio = s% Production_win_aspect_ratio
         p% file_flag = s% Production_file_flag
         p% file_dir = s% Production_file_dir
         p% file_prefix = s% Production_file_prefix
         p% file_interval = s% Production_file_interval
         p% file_width = s% Production_file_width
         p% file_aspect_ratio = s% Production_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Panels1)
         p% plot => History_Panels1_plot
         p% id = i_Hist_Panels1
         p% name = 'History_Panels1'
         p% win_flag = s% History_Panels1_win_flag
         p% win_width = s% History_Panels1_win_width
         p% win_aspect_ratio = s% History_Panels1_win_aspect_ratio
         p% file_flag = s% History_Panels1_file_flag
         p% file_dir = s% History_Panels1_file_dir
         p% file_prefix = s% History_Panels1_file_prefix
         p% file_interval = s% History_Panels1_file_interval
         p% file_width = s% History_Panels1_file_width
         p% file_aspect_ratio = s% History_Panels1_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Panels2)
         p% plot => History_Panels2_plot
         p% id = i_Hist_Panels2
         p% name = 'History_Panels2'
         p% win_flag = s% History_Panels2_win_flag
         p% win_width = s% History_Panels2_win_width
         p% win_aspect_ratio = s% History_Panels2_win_aspect_ratio
         p% file_flag = s% History_Panels2_file_flag
         p% file_dir = s% History_Panels2_file_dir
         p% file_prefix = s% History_Panels2_file_prefix
         p% file_interval = s% History_Panels2_file_interval
         p% file_width = s% History_Panels2_file_width
         p% file_aspect_ratio = s% History_Panels2_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Panels3)
         p% plot => History_Panels3_plot
         p% id = i_Hist_Panels3
         p% name = 'History_Panels3'
         p% win_flag = s% History_Panels3_win_flag
         p% win_width = s% History_Panels3_win_width
         p% win_aspect_ratio = s% History_Panels3_win_aspect_ratio
         p% file_flag = s% History_Panels3_file_flag
         p% file_dir = s% History_Panels3_file_dir
         p% file_prefix = s% History_Panels3_file_prefix
         p% file_interval = s% History_Panels3_file_interval
         p% file_width = s% History_Panels3_file_width
         p% file_aspect_ratio = s% History_Panels3_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Panels4)
         p% plot => History_Panels4_plot
         p% id = i_Hist_Panels4
         p% name = 'History_Panels4'
         p% win_flag = s% History_Panels4_win_flag
         p% win_width = s% History_Panels4_win_width
         p% win_aspect_ratio = s% History_Panels4_win_aspect_ratio
         p% file_flag = s% History_Panels4_file_flag
         p% file_dir = s% History_Panels4_file_dir
         p% file_prefix = s% History_Panels4_file_prefix
         p% file_interval = s% History_Panels4_file_interval
         p% file_width = s% History_Panels4_file_width
         p% file_aspect_ratio = s% History_Panels4_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Panels5)
         p% plot => History_Panels5_plot
         p% id = i_Hist_Panels5
         p% name = 'History_Panels5'
         p% win_flag = s% History_Panels5_win_flag
         p% win_width = s% History_Panels5_win_width
         p% win_aspect_ratio = s% History_Panels5_win_aspect_ratio
         p% file_flag = s% History_Panels5_file_flag
         p% file_dir = s% History_Panels5_file_dir
         p% file_prefix = s% History_Panels5_file_prefix
         p% file_interval = s% History_Panels5_file_interval
         p% file_width = s% History_Panels5_file_width
         p% file_aspect_ratio = s% History_Panels5_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Panels6)
         p% plot => History_Panels6_plot
         p% id = i_Hist_Panels6
         p% name = 'History_Panels6'
         p% win_flag = s% History_Panels6_win_flag
         p% win_width = s% History_Panels6_win_width
         p% win_aspect_ratio = s% History_Panels6_win_aspect_ratio
         p% file_flag = s% History_Panels6_file_flag
         p% file_dir = s% History_Panels6_file_dir
         p% file_prefix = s% History_Panels6_file_prefix
         p% file_interval = s% History_Panels6_file_interval
         p% file_width = s% History_Panels6_file_width
         p% file_aspect_ratio = s% History_Panels6_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Panels7)
         p% plot => History_Panels7_plot
         p% id = i_Hist_Panels7
         p% name = 'History_Panels7'
         p% win_flag = s% History_Panels7_win_flag
         p% win_width = s% History_Panels7_win_width
         p% win_aspect_ratio = s% History_Panels7_win_aspect_ratio
         p% file_flag = s% History_Panels7_file_flag
         p% file_dir = s% History_Panels7_file_dir
         p% file_prefix = s% History_Panels7_file_prefix
         p% file_interval = s% History_Panels7_file_interval
         p% file_width = s% History_Panels7_file_width
         p% file_aspect_ratio = s% History_Panels7_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Panels8)
         p% plot => History_Panels8_plot
         p% id = i_Hist_Panels8
         p% name = 'History_Panels8'
         p% win_flag = s% History_Panels8_win_flag
         p% win_width = s% History_Panels8_win_width
         p% win_aspect_ratio = s% History_Panels8_win_aspect_ratio
         p% file_flag = s% History_Panels8_file_flag
         p% file_dir = s% History_Panels8_file_dir
         p% file_prefix = s% History_Panels8_file_prefix
         p% file_interval = s% History_Panels8_file_interval
         p% file_width = s% History_Panels8_file_width
         p% file_aspect_ratio = s% History_Panels8_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Panels9)
         p% plot => History_Panels9_plot
         p% id = i_Hist_Panels9
         p% name = 'History_Panels9'
         p% win_flag = s% History_Panels9_win_flag
         p% win_width = s% History_Panels9_win_width
         p% win_aspect_ratio = s% History_Panels9_win_aspect_ratio
         p% file_flag = s% History_Panels9_file_flag
         p% file_dir = s% History_Panels9_file_dir
         p% file_prefix = s% History_Panels9_file_prefix
         p% file_interval = s% History_Panels9_file_interval
         p% file_width = s% History_Panels9_file_width
         p% file_aspect_ratio = s% History_Panels9_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Track1)
         p% plot => History_Track1_plot
         p% id = i_Hist_Track1
         p% name = 'History_Track1'
         p% win_flag = s% History_Track1_win_flag
         p% win_width = s% History_Track1_win_width
         p% win_aspect_ratio = s% History_Track1_win_aspect_ratio
         p% file_flag = s% History_Track1_file_flag
         p% file_dir = s% History_Track1_file_dir
         p% file_prefix = s% History_Track1_file_prefix
         p% file_interval = s% History_Track1_file_interval
         p% file_width = s% History_Track1_file_width
         p% file_aspect_ratio = s% History_Track1_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Track2)
         p% plot => History_Track2_plot
         p% id = i_Hist_Track2
         p% name = 'History_Track2'
         p% win_flag = s% History_Track2_win_flag
         p% win_width = s% History_Track2_win_width
         p% win_aspect_ratio = s% History_Track2_win_aspect_ratio
         p% file_flag = s% History_Track2_file_flag
         p% file_dir = s% History_Track2_file_dir
         p% file_prefix = s% History_Track2_file_prefix
         p% file_interval = s% History_Track2_file_interval
         p% file_width = s% History_Track2_file_width
         p% file_aspect_ratio = s% History_Track2_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Track3)
         p% plot => History_Track3_plot
         p% id = i_Hist_Track3
         p% name = 'History_Track3'
         p% win_flag = s% History_Track3_win_flag
         p% win_width = s% History_Track3_win_width
         p% win_aspect_ratio = s% History_Track3_win_aspect_ratio
         p% file_flag = s% History_Track3_file_flag
         p% file_dir = s% History_Track3_file_dir
         p% file_prefix = s% History_Track3_file_prefix
         p% file_interval = s% History_Track3_file_interval
         p% file_width = s% History_Track3_file_width
         p% file_aspect_ratio = s% History_Track3_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Track4)
         p% plot => History_Track4_plot
         p% id = i_Hist_Track4
         p% name = 'History_Track4'
         p% win_flag = s% History_Track4_win_flag
         p% win_width = s% History_Track4_win_width
         p% win_aspect_ratio = s% History_Track4_win_aspect_ratio
         p% file_flag = s% History_Track4_file_flag
         p% file_dir = s% History_Track4_file_dir
         p% file_prefix = s% History_Track4_file_prefix
         p% file_interval = s% History_Track4_file_interval
         p% file_width = s% History_Track4_file_width
         p% file_aspect_ratio = s% History_Track4_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Track5)
         p% plot => History_Track5_plot
         p% id = i_Hist_Track5
         p% name = 'History_Track5'
         p% win_flag = s% History_Track5_win_flag
         p% win_width = s% History_Track5_win_width
         p% win_aspect_ratio = s% History_Track5_win_aspect_ratio
         p% file_flag = s% History_Track5_file_flag
         p% file_dir = s% History_Track5_file_dir
         p% file_prefix = s% History_Track5_file_prefix
         p% file_interval = s% History_Track5_file_interval
         p% file_width = s% History_Track5_file_width
         p% file_aspect_ratio = s% History_Track5_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Track6)
         p% plot => History_Track6_plot
         p% id = i_Hist_Track6
         p% name = 'History_Track6'
         p% win_flag = s% History_Track6_win_flag
         p% win_width = s% History_Track6_win_width
         p% win_aspect_ratio = s% History_Track6_win_aspect_ratio
         p% file_flag = s% History_Track6_file_flag
         p% file_dir = s% History_Track6_file_dir
         p% file_prefix = s% History_Track6_file_prefix
         p% file_interval = s% History_Track6_file_interval
         p% file_width = s% History_Track6_file_width
         p% file_aspect_ratio = s% History_Track6_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Track7)
         p% plot => History_Track7_plot
         p% id = i_Hist_Track7
         p% name = 'History_Track7'
         p% win_flag = s% History_Track7_win_flag
         p% win_width = s% History_Track7_win_width
         p% win_aspect_ratio = s% History_Track7_win_aspect_ratio
         p% file_flag = s% History_Track7_file_flag
         p% file_dir = s% History_Track7_file_dir
         p% file_prefix = s% History_Track7_file_prefix
         p% file_interval = s% History_Track7_file_interval
         p% file_width = s% History_Track7_file_width
         p% file_aspect_ratio = s% History_Track7_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Track8)
         p% plot => History_Track8_plot
         p% id = i_Hist_Track8
         p% name = 'History_Track8'
         p% win_flag = s% History_Track8_win_flag
         p% win_width = s% History_Track8_win_width
         p% win_aspect_ratio = s% History_Track8_win_aspect_ratio
         p% file_flag = s% History_Track8_file_flag
         p% file_dir = s% History_Track8_file_dir
         p% file_prefix = s% History_Track8_file_prefix
         p% file_interval = s% History_Track8_file_interval
         p% file_width = s% History_Track8_file_width
         p% file_aspect_ratio = s% History_Track8_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Hist_Track9)
         p% plot => History_Track9_plot
         p% id = i_Hist_Track9
         p% name = 'History_Track9'
         p% win_flag = s% History_Track9_win_flag
         p% win_width = s% History_Track9_win_width
         p% win_aspect_ratio = s% History_Track9_win_aspect_ratio
         p% file_flag = s% History_Track9_file_flag
         p% file_dir = s% History_Track9_file_dir
         p% file_prefix = s% History_Track9_file_prefix
         p% file_interval = s% History_Track9_file_interval
         p% file_width = s% History_Track9_file_width
         p% file_aspect_ratio = s% History_Track9_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Mode_Prop)
         p% plot => mode_propagation_plot
         p% id = i_Mode_Prop
         p% name = 'Mode_Propagation'
         p% win_flag = s% Mode_Prop_win_flag
         p% win_width = s% Mode_Prop_win_width
         p% win_aspect_ratio = s% Mode_Prop_win_aspect_ratio
         p% file_flag = s% Mode_Prop_file_flag
         p% file_dir = s% Mode_Prop_file_dir
         p% file_prefix = s% Mode_Prop_file_prefix
         p% file_interval = s% Mode_Prop_file_interval
         p% file_width = s% Mode_Prop_file_width
         p% file_aspect_ratio = s% Mode_Prop_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Power)
         p% plot => power_plot
         p% id = i_Power
         p% name = 'Power'
         p% win_flag = s% Power_win_flag
         p% win_width = s% Power_win_width
         p% win_aspect_ratio = s% Power_win_aspect_ratio
         p% file_flag = s% Power_file_flag
         p% file_dir = s% Power_file_dir
         p% file_prefix = s% Power_file_prefix
         p% file_interval = s% Power_file_interval
         p% file_width = s% Power_file_width
         p% file_aspect_ratio = s% Power_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Abundance)
         p% plot => abundance_plot
         p% id = i_Abundance
         p% name = 'Abundance'
         p% win_flag = s% Abundance_win_flag
         p% win_width = s% Abundance_win_width
         p% win_aspect_ratio = s% Abundance_win_aspect_ratio
         p% file_flag = s% Abundance_file_flag
         p% file_dir = s% Abundance_file_dir
         p% file_prefix = s% Abundance_file_prefix
         p% file_interval = s% Abundance_file_interval
         p% file_width = s% Abundance_file_width
         p% file_aspect_ratio = s% Abundance_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Summary_Burn)
         p% plot => summary_burn_plot
         p% id = i_Summary_Burn
         p% name = 'Summary_Burn'
         p% win_flag = s% Summary_Burn_win_flag
         p% win_width = s% Summary_Burn_win_width
         p% win_aspect_ratio = s% Summary_Burn_win_aspect_ratio
         p% file_flag = s% Summary_Burn_file_flag
         p% file_dir = s% Summary_Burn_file_dir
         p% file_prefix = s% Summary_Burn_file_prefix
         p% file_interval = s% Summary_Burn_file_interval
         p% file_width = s% Summary_Burn_file_width
         p% file_aspect_ratio = s% Summary_Burn_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Summary_Profile)
         p% plot => summary_profile_plot
         p% id = i_Summary_Profile
         p% name = 'Summary_Profile'
         p% win_flag = s% Summary_Profile_win_flag
         p% win_width = s% Summary_Profile_win_width
         p% win_aspect_ratio = s% Summary_Profile_win_aspect_ratio
         p% file_flag = s% Summary_Profile_file_flag
         p% file_dir = s% Summary_Profile_file_dir
         p% file_prefix = s% Summary_Profile_file_prefix
         p% file_interval = s% Summary_Profile_file_interval
         p% file_width = s% Summary_Profile_file_width
         p% file_aspect_ratio = s% Summary_Profile_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Summary_History)
         p% plot => summary_history_plot
         p% id = i_Summary_History
         p% name = 'Summary_History'
         p% win_flag = s% Summary_History_win_flag
         p% win_width = s% Summary_History_win_width
         p% win_aspect_ratio = s% Summary_History_win_aspect_ratio
         p% file_flag = s% Summary_History_file_flag
         p% file_dir = s% Summary_History_file_dir
         p% file_prefix = s% Summary_History_file_prefix
         p% file_interval = s% Summary_History_file_interval
         p% file_width = s% Summary_History_file_width
         p% file_aspect_ratio = s% Summary_History_file_aspect_ratio


         p => s% pgstar_win_file_ptr(i_Col_Mag1)
         p% plot => Color_Magnitude1_plot
         p% id = i_Col_Mag1
         p% name = 'Color_Magnitude1'
         p% win_flag = s% Color_Magnitude1_win_flag
         p% win_width = s% Color_Magnitude1_win_width
         p% win_aspect_ratio = s% Color_Magnitude1_win_aspect_ratio
         p% file_flag = s% Color_Magnitude1_file_flag
         p% file_dir = s% Color_Magnitude1_file_dir
         p% file_prefix = s% Color_Magnitude1_file_prefix
         p% file_interval = s% Color_Magnitude1_file_interval
         p% file_width = s% Color_Magnitude1_file_width
         p% file_aspect_ratio = s% Color_Magnitude1_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Col_Mag2)
         p% plot => Color_Magnitude2_plot
         p% id = i_Col_Mag2
         p% name = 'Color_Magnitude2'
         p% win_flag = s% Color_Magnitude2_win_flag
         p% win_width = s% Color_Magnitude2_win_width
         p% win_aspect_ratio = s% Color_Magnitude2_win_aspect_ratio
         p% file_flag = s% Color_Magnitude2_file_flag
         p% file_dir = s% Color_Magnitude2_file_dir
         p% file_prefix = s% Color_Magnitude2_file_prefix
         p% file_interval = s% Color_Magnitude2_file_interval
         p% file_width = s% Color_Magnitude2_file_width
         p% file_aspect_ratio = s% Color_Magnitude2_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Col_Mag3)
         p% plot => Color_Magnitude3_plot
         p% id = i_Col_Mag3
         p% name = 'Color_Magnitude3'
         p% win_flag = s% Color_Magnitude3_win_flag
         p% win_width = s% Color_Magnitude3_win_width
         p% win_aspect_ratio = s% Color_Magnitude3_win_aspect_ratio
         p% file_flag = s% Color_Magnitude3_file_flag
         p% file_dir = s% Color_Magnitude3_file_dir
         p% file_prefix = s% Color_Magnitude3_file_prefix
         p% file_interval = s% Color_Magnitude3_file_interval
         p% file_width = s% Color_Magnitude3_file_width
         p% file_aspect_ratio = s% Color_Magnitude3_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Col_Mag4)
         p% plot => Color_Magnitude4_plot
         p% id = i_Col_Mag4
         p% name = 'Color_Magnitude4'
         p% win_flag = s% Color_Magnitude4_win_flag
         p% win_width = s% Color_Magnitude4_win_width
         p% win_aspect_ratio = s% Color_Magnitude4_win_aspect_ratio
         p% file_flag = s% Color_Magnitude4_file_flag
         p% file_dir = s% Color_Magnitude4_file_dir
         p% file_prefix = s% Color_Magnitude4_file_prefix
         p% file_interval = s% Color_Magnitude4_file_interval
         p% file_width = s% Color_Magnitude4_file_width
         p% file_aspect_ratio = s% Color_Magnitude4_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Col_Mag5)
         p% plot => Color_Magnitude5_plot
         p% id = i_Col_Mag5
         p% name = 'Color_Magnitude5'
         p% win_flag = s% Color_Magnitude5_win_flag
         p% win_width = s% Color_Magnitude5_win_width
         p% win_aspect_ratio = s% Color_Magnitude5_win_aspect_ratio
         p% file_flag = s% Color_Magnitude5_file_flag
         p% file_dir = s% Color_Magnitude5_file_dir
         p% file_prefix = s% Color_Magnitude5_file_prefix
         p% file_interval = s% Color_Magnitude5_file_interval
         p% file_width = s% Color_Magnitude5_file_width
         p% file_aspect_ratio = s% Color_Magnitude5_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Col_Mag6)
         p% plot => Color_Magnitude6_plot
         p% id = i_Col_Mag6
         p% name = 'Color_Magnitude6'
         p% win_flag = s% Color_Magnitude6_win_flag
         p% win_width = s% Color_Magnitude6_win_width
         p% win_aspect_ratio = s% Color_Magnitude6_win_aspect_ratio
         p% file_flag = s% Color_Magnitude6_file_flag
         p% file_dir = s% Color_Magnitude6_file_dir
         p% file_prefix = s% Color_Magnitude6_file_prefix
         p% file_interval = s% Color_Magnitude6_file_interval
         p% file_width = s% Color_Magnitude6_file_width
         p% file_aspect_ratio = s% Color_Magnitude6_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Col_Mag7)
         p% plot => Color_Magnitude7_plot
         p% id = i_Col_Mag7
         p% name = 'Color_Magnitude7'
         p% win_flag = s% Color_Magnitude7_win_flag
         p% win_width = s% Color_Magnitude7_win_width
         p% win_aspect_ratio = s% Color_Magnitude7_win_aspect_ratio
         p% file_flag = s% Color_Magnitude7_file_flag
         p% file_dir = s% Color_Magnitude7_file_dir
         p% file_prefix = s% Color_Magnitude7_file_prefix
         p% file_interval = s% Color_Magnitude7_file_interval
         p% file_width = s% Color_Magnitude7_file_width
         p% file_aspect_ratio = s% Color_Magnitude7_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Col_Mag8)
         p% plot => Color_Magnitude8_plot
         p% id = i_Col_Mag8
         p% name = 'Color_Magnitude8'
         p% win_flag = s% Color_Magnitude8_win_flag
         p% win_width = s% Color_Magnitude8_win_width
         p% win_aspect_ratio = s% Color_Magnitude8_win_aspect_ratio
         p% file_flag = s% Color_Magnitude8_file_flag
         p% file_dir = s% Color_Magnitude8_file_dir
         p% file_prefix = s% Color_Magnitude8_file_prefix
         p% file_interval = s% Color_Magnitude8_file_interval
         p% file_width = s% Color_Magnitude8_file_width
         p% file_aspect_ratio = s% Color_Magnitude8_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Col_Mag9)
         p% plot => Color_Magnitude9_plot
         p% id = i_Col_Mag9
         p% name = 'Color_Magnitude9'
         p% win_flag = s% Color_Magnitude9_win_flag
         p% win_width = s% Color_Magnitude9_win_width
         p% win_aspect_ratio = s% Color_Magnitude9_win_aspect_ratio
         p% file_flag = s% Color_Magnitude9_file_flag
         p% file_dir = s% Color_Magnitude9_file_dir
         p% file_prefix = s% Color_Magnitude9_file_prefix
         p% file_interval = s% Color_Magnitude9_file_interval
         p% file_width = s% Color_Magnitude9_file_width
         p% file_aspect_ratio = s% Color_Magnitude9_file_aspect_ratio


         p => s% pgstar_win_file_ptr(i_Grid1)
         p% plot => grid1_plot
         p% id = i_Grid1
         p% name = 'Grid1'
         p% win_flag = s% Grid1_win_flag
         p% win_width = s% Grid1_win_width
         p% win_aspect_ratio = s% Grid1_win_aspect_ratio
         p% file_flag = s% Grid1_file_flag
         p% file_dir = s% Grid1_file_dir
         p% file_prefix = s% Grid1_file_prefix
         p% file_interval = s% Grid1_file_interval
         p% file_width = s% Grid1_file_width
         p% file_aspect_ratio = s% Grid1_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Grid2)
         p% plot => grid2_plot
         p% id = i_Grid2
         p% name = 'Grid2'
         p% win_flag = s% Grid2_win_flag
         p% win_width = s% Grid2_win_width
         p% win_aspect_ratio = s% Grid2_win_aspect_ratio
         p% file_flag = s% Grid2_file_flag
         p% file_dir = s% Grid2_file_dir
         p% file_prefix = s% Grid2_file_prefix
         p% file_interval = s% Grid2_file_interval
         p% file_width = s% Grid2_file_width
         p% file_aspect_ratio = s% Grid2_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Grid3)
         p% plot => grid3_plot
         p% id = i_Grid3
         p% name = 'Grid3'
         p% win_flag = s% Grid3_win_flag
         p% win_width = s% Grid3_win_width
         p% win_aspect_ratio = s% Grid3_win_aspect_ratio
         p% file_flag = s% Grid3_file_flag
         p% file_dir = s% Grid3_file_dir
         p% file_prefix = s% Grid3_file_prefix
         p% file_interval = s% Grid3_file_interval
         p% file_width = s% Grid3_file_width
         p% file_aspect_ratio = s% Grid3_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Grid4)
         p% plot => grid4_plot
         p% id = i_Grid4
         p% name = 'Grid4'
         p% win_flag = s% Grid4_win_flag
         p% win_width = s% Grid4_win_width
         p% win_aspect_ratio = s% Grid4_win_aspect_ratio
         p% file_flag = s% Grid4_file_flag
         p% file_dir = s% Grid4_file_dir
         p% file_prefix = s% Grid4_file_prefix
         p% file_interval = s% Grid4_file_interval
         p% file_width = s% Grid4_file_width
         p% file_aspect_ratio = s% Grid4_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Grid5)
         p% plot => grid5_plot
         p% id = i_Grid5
         p% name = 'Grid5'
         p% win_flag = s% Grid5_win_flag
         p% win_width = s% Grid5_win_width
         p% win_aspect_ratio = s% Grid5_win_aspect_ratio
         p% file_flag = s% Grid5_file_flag
         p% file_dir = s% Grid5_file_dir
         p% file_prefix = s% Grid5_file_prefix
         p% file_interval = s% Grid5_file_interval
         p% file_width = s% Grid5_file_width
         p% file_aspect_ratio = s% Grid5_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Grid6)
         p% plot => grid6_plot
         p% id = i_Grid6
         p% name = 'Grid6'
         p% win_flag = s% Grid6_win_flag
         p% win_width = s% Grid6_win_width
         p% win_aspect_ratio = s% Grid6_win_aspect_ratio
         p% file_flag = s% Grid6_file_flag
         p% file_dir = s% Grid6_file_dir
         p% file_prefix = s% Grid6_file_prefix
         p% file_interval = s% Grid6_file_interval
         p% file_width = s% Grid6_file_width
         p% file_aspect_ratio = s% Grid6_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Grid7)
         p% plot => grid7_plot
         p% id = i_Grid7
         p% name = 'Grid7'
         p% win_flag = s% Grid7_win_flag
         p% win_width = s% Grid7_win_width
         p% win_aspect_ratio = s% Grid7_win_aspect_ratio
         p% file_flag = s% Grid7_file_flag
         p% file_dir = s% Grid7_file_dir
         p% file_prefix = s% Grid7_file_prefix
         p% file_interval = s% Grid7_file_interval
         p% file_width = s% Grid7_file_width
         p% file_aspect_ratio = s% Grid7_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Grid8)
         p% plot => grid8_plot
         p% id = i_Grid8
         p% name = 'Grid8'
         p% win_flag = s% Grid8_win_flag
         p% win_width = s% Grid8_win_width
         p% win_aspect_ratio = s% Grid8_win_aspect_ratio
         p% file_flag = s% Grid8_file_flag
         p% file_dir = s% Grid8_file_dir
         p% file_prefix = s% Grid8_file_prefix
         p% file_interval = s% Grid8_file_interval
         p% file_width = s% Grid8_file_width
         p% file_aspect_ratio = s% Grid8_file_aspect_ratio

         p => s% pgstar_win_file_ptr(i_Grid9)
         p% plot => grid9_plot
         p% id = i_Grid9
         p% name = 'Grid9'
         p% win_flag = s% Grid9_win_flag
         p% win_width = s% Grid9_win_width
         p% win_aspect_ratio = s% Grid9_win_aspect_ratio
         p% file_flag = s% Grid9_file_flag
         p% file_dir = s% Grid9_file_dir
         p% file_prefix = s% Grid9_file_prefix
         p% file_interval = s% Grid9_file_interval
         p% file_width = s% Grid9_file_width
         p% file_aspect_ratio = s% Grid9_file_aspect_ratio

         do i = 1, max_num_Other_plots
            p => s% pgstar_win_file_ptr(i_Other + i - 1)
            p% win_flag = .false.
            p% file_flag = .false.
            p% okay_to_call_do_plot_in_grid = .false.
         end do

         if (s% use_other_pgstar_plots) &
            call s% other_pgstar_plots_info(s% id, ierr)

      end subroutine set_win_file_data


      subroutine do_pgstar_plots( &
            s, must_write_files, &
            ierr)
         type (star_info), pointer :: s
         logical, intent(in) :: must_write_files
         integer, intent(out) :: ierr

         integer :: i
         integer(8) :: time0, time1, clock_rate
         logical :: pause

         include 'formats'

         ierr = 0

         if (s% clear_history) call pgstar_clear(s)

         call update_pgstar_data(s, ierr)
         if (failed('update_pgstar_data')) return
         
         call onScreen_Plots(s, must_write_files, ierr)
         if (failed('onScreen_Plots')) return
         
         call update_pgstar_history_file(s,ierr)
         if (failed('save_text_data')) return
         
         if (s% pause_interval > 0) then
            pause = (mod(s% model_number, s% pause_interval) == 0)
         else
            pause = s% pause
         end if
         
         if (pause .and. s% pgstar_interval > 0) &
            pause = (mod(s% model_number, s% pgstar_interval) == 0)
            
         if (pause) then
            write(*,*)
            write(*,*) 'model_number', s% model_number
            write(*,*) 'PGSTAR: paused -- hit RETURN to continue'
            read(*,*)
         end if

         if (s% pgstar_sleep > 0) then
            time0 = s% system_clock_at_start_of_step
            do
               call system_clock(time1,clock_rate)
               if (dble(time1 - time0)/dble(clock_rate) >= s% pgstar_sleep) exit
            end do
         end if

         !write(*,2) 'PGSTAR: done', s% model_number

         contains

         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed

      end subroutine do_pgstar_plots


      ! PGSTAR driver, called after each timestep
      subroutine onScreen_Plots(s, must_write_files_in, ierr)
         use utils_lib
         use chem_def
         use net_def
         use net_lib, only: get_net_reaction_table
         use rates_def, only: rates_reaction_id_max
         use const_def, only: Msun, Rsun

         type (star_info), pointer :: s
         logical :: must_write_files_in
         integer, intent(out) :: ierr

         integer :: i
         type (pgstar_win_file_data), pointer :: p
         logical, parameter :: dbg = .false.
         real(dp) :: dlgL, dlgTeff, dHR
         logical :: must_write_files,show_plot_now,save_plot_now

         include 'formats'
         ierr = 0

         ! initialize pgstar
         if ( .not. have_initialized_pgstar ) then
            call init_pgstar(s, s% model_number, ierr)
            if (failed('init_pgstar')) return
         end if

         ! request files if sufficient movement in HR diagram
         must_write_files = must_write_files_in

         if (s% delta_HR_limit_for_file_output > 0 .and. &
               s% L_phot_old > 0 .and. s% Teff_old > 0 .and. .not. must_write_files) then
            dlgL = log10(s% L_phot/s% L_phot_old)
            dlgTeff = log10(s% Teff/s% Teff_old)
            dHR = sqrt(pow2(s% delta_HR_ds_L*dlgL) + pow2(s% delta_HR_ds_Teff*dlgTeff))
            sum_dHR_since_last_file_write = sum_dHR_since_last_file_write + dHR
            must_write_files = &
               (sum_dHR_since_last_file_write >= s% delta_HR_limit_for_file_output)
         end if
         if (must_write_files) sum_dHR_since_last_file_write = 0

         show_plot_now=.false.
         if (s% pgstar_interval > 0) then
            if(mod(s% model_number, s% pgstar_interval) == 0) then
               show_plot_now=.true.
            end if
         end if

         ! retrieve extra profile data
         s% num_extra_profile_cols = 0
         if (associated(s% how_many_extra_profile_columns) .and. &
             associated(s% data_for_extra_profile_columns)) then
            i = s% how_many_extra_profile_columns(s% id)
            if (i > 0) then
               if (associated(s% extra_profile_col_names)) &
                  deallocate(s% extra_profile_col_names)
               if (associated(s% extra_profile_col_vals)) &
                  deallocate(s% extra_profile_col_vals)
               allocate(s% extra_profile_col_names(i))
               allocate(s% extra_profile_col_vals(s% nz,i))
               ierr = 0
               call s% data_for_extra_profile_columns( &
                  s% id, i, s% nz, &
                  s% extra_profile_col_names, &
                  s% extra_profile_col_vals, ierr)
               if (ierr == 0) s% num_extra_profile_cols = i
               ierr = 0
            end if
         end if

         ! loop through all plots
         do i = 1, num_pgstar_plots
            p => s% pgstar_win_file_ptr(i)

            if(show_plot_now) then
               ! call to check_window opens device
               call check_window(s, p, ierr)
               if (failed('check_window')) return
   
               ! make the plot (window)
               if (p% do_win) then
                  call p% plot(s% id, p% id_win, ierr)
                  if (failed(p% name)) return
               end if
            end if
   
            save_plot_now=must_write_files
            if (p% file_interval > 0) then
               if(mod(s% model_number, p% file_interval) == 0) then
                  save_plot_now=.true.
               end if
            end if

            if(save_plot_now)then
               ! call to check_file opens device and does mkdir
               call check_file(s, p, ierr)
   
               ! make the plot (file)
               if (p% do_file) then
                  call p% plot(s% id, p% id_file, ierr)
                  if (failed(p% name)) return
                  call pgclos
                  if (s% pgstar_report_writing_files) &
                     write(*,*) trim(p% most_recent_filename)
                  p% id_file = 0
                  p% do_file = .false.
               end if
            end if
         end do

         contains

         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed

      end subroutine onScreen_Plots


      subroutine update_pgstar_history_file(s, ierr)
         use utils_lib
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: iounit, i, n
         character (len=1024) :: fname
         type (pgstar_hist_node), pointer :: pg

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         pg => s% pgstar_hist
         if (.not. associated(pg)) return

         n = s% number_of_history_columns
         fname = trim(s% log_directory) // '/pgstar.dat'

         if (associated(pg% next)) then
            open(newunit=iounit, file=trim(fname), action='write', &
               position='append', form='unformatted', iostat=ierr)
         else
            open(newunit=iounit, file=trim(fname), action='write', &
               status='replace', form='unformatted', iostat=ierr)
            if (ierr == 0) write(iounit) n
         end if
         if (ierr /= 0) then
            write(*,*) 'save_pgstar_data: cannot open new file'
            return
         end if

         if (associated(pg% vals)) then
            if (size(pg% vals,dim=1) >= n) then
               write(iounit) pg% age, pg% step, pg% vals(1:n)
            end if
         end if

         close(iounit)

      end subroutine update_pgstar_history_file


      subroutine read_pgstar_data(s, ierr)
         use utils_lib
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         logical :: fexist
         integer :: iounit, i, n
         character (len=1024) :: fname
         type (pgstar_hist_node), pointer :: pg

         logical, parameter :: dbg = .false.

         include 'formats'
         ierr = 0

         fname = trim(s% log_directory) // '/pgstar.dat'
         inquire(file=trim(fname), exist=fexist)
         if (.not.fexist) then
            if (dbg) write(*,*) 'failed to find ' // trim(fname)
            return
         end if

         open(newunit=iounit, file=trim(fname), action='read', &
                  status='old', iostat=ierr, form='unformatted')
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed to open ' // trim(fname)
            return
         end if

         read(iounit, iostat=ierr) n
         if (ierr == 0) then
            if (s% number_of_history_columns < 0) then
               s% number_of_history_columns = n
            else if (s% number_of_history_columns /= n) then
               ierr = -1
            end if
         end if

         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed read pg_star history ' // trim(fname)
         else
            do ! keep reading until reach end of file so take care of restarts
               allocate(pg)
               allocate(pg% vals(n))
               read(iounit, iostat=ierr) pg% age, pg% step, pg% vals(1:n)
               if (ierr /= 0) then
                  ierr = 0
                  deallocate(pg% vals)
                  deallocate(pg)
                  exit
               end if
               call add_to_pgstar_hist(s, pg)
            end do
         end if

         close(iounit)

      end subroutine read_pgstar_data


      subroutine update_pgstar_data(s, ierr)
         use star_utils, only: eval_csound

         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: num, i
         type (pgstar_hist_node), pointer :: pg

         include 'formats'

         ierr = 0
         allocate(pg)
         pg% step = s% model_number
         pg% age = s% star_age
         num = s% number_of_history_columns
         allocate(pg% vals(num))
         call get_hist_values(num,ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_hist_values'
            return
            stop 'pgstar'
         end if
         call add_to_pgstar_hist(s, pg)

         contains

         subroutine get_hist_values(num,ierr)
            use history, only: do_get_data_for_history_columns
            integer, intent(in) :: num
            integer, intent(out) :: ierr
            integer :: i
            ierr = 0
            if (s% need_to_set_history_names_etc .or. &
                  s% model_number_of_history_values /= s% model_number) then
               call do_get_data_for_history_columns( &
                  s, &
                  ierr)
               if (ierr /= 0) return
            end if
            do i=1,num
               pg% vals(i) = s% history_values(i)
            end do
         end subroutine get_hist_values

      end subroutine update_pgstar_data


      subroutine do_set_xaxis_bounds( &
            s, xaxis_by, win_xmin_in, win_xmax_in, xmargin, &
            xvec, xmin, xmax, xleft, xright, dx, &
            grid_min, grid_max, npts, ierr)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: xaxis_by
         real, intent(in) :: win_xmin_in, win_xmax_in, xmargin
         real, allocatable, dimension(:) :: xvec
         real, intent(out) :: xmin, xmax, xleft, xright, dx
         integer, intent(out) :: grid_min, grid_max, npts
         integer, intent(out) :: ierr
         call set_xaxis_bounds( &
            s, xaxis_by, win_xmin_in, win_xmax_in, .false., xmargin, &
            xvec, xmin, xmax, xleft, xright, dx, &
            grid_min, grid_max, npts, ierr)
      end subroutine do_set_xaxis_bounds


      subroutine do_show_xaxis_by(s,by,ierr)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: by
         integer, intent(out) :: ierr
         call show_xaxis_name(s,by,ierr)
      end subroutine do_show_xaxis_by

      subroutine shutdown_pgstar(s)
         use pgstar_support
         type (star_info), pointer :: s
         
         if(have_initialized_pgstar) then
            call dealloc(hydrogen_burn_logT)
            call dealloc(hydrogen_burn_logRho)
            call dealloc(helium_burn_logT)
            call dealloc(helium_burn_logRho)
            call dealloc(carbon_burn_logT)
            call dealloc(carbon_burn_logRho)
            call dealloc(oxygen_burn_logT)
            call dealloc(oxygen_burn_logRho)
            call dealloc(psi4_logT)
            call dealloc(psi4_logRho)
            call dealloc(elect_data_logT)
            call dealloc(elect_data_logRho)
            call dealloc(gamma_4_thirds_logT)
            call dealloc(gamma_4_thirds_logRho)
            call dealloc(kap_rad_cond_eq_logT)
            call dealloc(kap_rad_cond_eq_logRho)
            call dealloc(opal_clip_logT)
            call dealloc(opal_clip_logRho)
            call dealloc(scvh_clip_logT)
            call dealloc(scvh_clip_logRho)
         end if
         
         
         call pgstar_clear(s)
         
         have_initialized_pgstar = .false.

         contains 
         
         subroutine dealloc(x)
            real, dimension(:), allocatable :: x
            
            if(allocated(x))then
               deallocate(x)
            end if
         
         end subroutine dealloc

      end subroutine shutdown_pgstar

      end module pgstar
