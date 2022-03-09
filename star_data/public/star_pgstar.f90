module star_pgstar

   use const_def
   use chem_def

   implicit none

      ! pgstar data

   abstract interface
      
      subroutine pgstar_plot_interface(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
      end subroutine pgstar_plot_interface

      subroutine other_do_plot_in_grid_interface( &
            id, device_id, xleft, xright, ybot, ytop, txt_scale, ierr)
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         integer, intent(out) :: ierr
      end subroutine other_do_plot_in_grid_interface

      subroutine pgstar_decorator_interface(id, xmin, xmax, ymin, ymax, plot_num, ierr)
         integer, intent(in) :: id
         !Not doubles
         real, intent(in) :: xmin, xmax, ymin, ymax 
         integer, intent(in) :: plot_num
         integer, intent(out) :: ierr
      end subroutine pgstar_decorator_interface

   end interface

   type pgstar_win_file_data
      integer :: id
      character (len=64) :: name
      logical :: win_flag, file_flag, do_win, do_file
      integer :: id_win, id_file, file_interval
      real :: win_width, prev_win_width
      real :: win_aspect_ratio, prev_win_aspect_ratio
      real :: file_width, file_aspect_ratio
      character (len=strlen) :: file_dir, file_prefix, most_recent_filename
      character (len=strlen) :: file_dir_for_previous_mkdir
      logical :: have_called_mkdir
      procedure (pgstar_plot_interface), pointer, nopass :: plot => null()
      ! the following make it possible to use "other" plots in grids
      logical :: okay_to_call_do_plot_in_grid
      procedure (other_do_plot_in_grid_interface), pointer, nopass :: &
         do_plot_in_grid => null()
   end type pgstar_win_file_data

   type pgstar_hist_node
      real(dp) :: age
      integer :: step
      real(dp), pointer :: vals(:) => null() ! values of items in history_columns list
      type (pgstar_hist_node), pointer :: next => null()
         ! list kept in strictly decreasing order by age & step
   end type pgstar_hist_node  


   integer, parameter :: max_Abundance_num_isos_to_show = 1000

   integer, parameter :: max_num_Profile_Panels = 16
   integer, parameter :: max_num_History_Panels = 16
   integer, parameter :: max_num_Color_Magnitude = 16
   integer, parameter :: max_num_Summary_Profile_Lines = 16
   integer, parameter :: max_num_Summary_History_Lines = 16
   integer, parameter :: max_num_Other_plots = 16

   integer, parameter :: i_TRho_Profile = 1
   integer, parameter :: i_logg_logT = i_TRho_Profile + 1
   integer, parameter :: i_logg_Teff = i_logg_logT + 1
   integer, parameter :: i_dPg_dnu = i_logg_Teff + 1
   integer, parameter :: i_L_R = i_dPg_dnu + 1
   integer, parameter :: i_L_v = i_L_R + 1
   integer, parameter :: i_L_Teff = i_L_v + 1
   integer, parameter :: i_R_L = i_L_Teff + 1
   integer, parameter :: i_R_Teff = i_R_L + 1
   integer, parameter :: i_logL_Teff = i_R_Teff + 1
   integer, parameter :: i_HR = i_logL_Teff + 1
   integer, parameter :: i_TRho = i_HR + 1
   integer, parameter :: i_TmaxRho = i_TRho + 1
   integer, parameter :: i_Dynamo = i_TRho + 1
   integer, parameter :: i_Mixing = i_Dynamo + 1
   integer, parameter :: i_rti = i_Mixing + 1
   integer, parameter :: i_Kipp = i_rti + 1
   integer, parameter :: i_Network = i_Kipp + 1
   integer, parameter :: i_Production = i_Network + 1
   integer, parameter :: i_Cntr_Hist = i_Production + 1
   integer, parameter :: i_Surf_Hist = i_Cntr_Hist + 1
   integer, parameter :: i_Mode_Prop = i_Surf_Hist + 1
   integer, parameter :: i_Power = i_Mode_Prop + 1
   integer, parameter :: i_Abundance = i_Power + 1
   integer, parameter :: i_Summary_History = i_Abundance + 1
   integer, parameter :: i_Summary_Burn = i_Summary_History + 1
   integer, parameter :: i_Summary_Profile = i_Summary_Burn + 1

   integer, parameter :: i_Text_Summary1 = i_Summary_Profile + 1
   integer, parameter :: i_Text_Summary2 = i_Text_Summary1 + 1
   integer, parameter :: i_Text_Summary3 = i_Text_Summary2 + 1
   integer, parameter :: i_Text_Summary4 = i_Text_Summary3 + 1
   integer, parameter :: i_Text_Summary5 = i_Text_Summary4 + 1
   integer, parameter :: i_Text_Summary6 = i_Text_Summary5 + 1
   integer, parameter :: i_Text_Summary7 = i_Text_Summary6 + 1
   integer, parameter :: i_Text_Summary8 = i_Text_Summary7 + 1
   integer, parameter :: i_Text_Summary9 = i_Text_Summary8 + 1

   integer, parameter :: i_Profile_Panels1 = i_Text_Summary9 + 1
   integer, parameter :: i_Profile_Panels2 = i_Profile_Panels1 + 1
   integer, parameter :: i_Profile_Panels3 = i_Profile_Panels2 + 1
   integer, parameter :: i_Profile_Panels4 = i_Profile_Panels3 + 1
   integer, parameter :: i_Profile_Panels5 = i_Profile_Panels4 + 1
   integer, parameter :: i_Profile_Panels6 = i_Profile_Panels5 + 1
   integer, parameter :: i_Profile_Panels7 = i_Profile_Panels6 + 1
   integer, parameter :: i_Profile_Panels8 = i_Profile_Panels7 + 1
   integer, parameter :: i_Profile_Panels9 = i_Profile_Panels8 + 1

   integer, parameter :: i_Hist_Track1 = i_Profile_Panels9 + 1
   integer, parameter :: i_Hist_Track2 = i_Hist_Track1 + 1
   integer, parameter :: i_Hist_Track3 = i_Hist_Track2 + 1
   integer, parameter :: i_Hist_Track4 = i_Hist_Track3 + 1
   integer, parameter :: i_Hist_Track5 = i_Hist_Track4 + 1
   integer, parameter :: i_Hist_Track6 = i_Hist_Track5 + 1
   integer, parameter :: i_Hist_Track7 = i_Hist_Track6 + 1
   integer, parameter :: i_Hist_Track8 = i_Hist_Track7 + 1
   integer, parameter :: i_Hist_Track9 = i_Hist_Track8 + 1

   integer, parameter :: i_Hist_Panels1 = i_Hist_Track9 + 1
   integer, parameter :: i_Hist_Panels2 = i_Hist_Panels1 + 1
   integer, parameter :: i_Hist_Panels3 = i_Hist_Panels2 + 1
   integer, parameter :: i_Hist_Panels4 = i_Hist_Panels3 + 1
   integer, parameter :: i_Hist_Panels5 = i_Hist_Panels4 + 1
   integer, parameter :: i_Hist_Panels6 = i_Hist_Panels5 + 1
   integer, parameter :: i_Hist_Panels7 = i_Hist_Panels6 + 1
   integer, parameter :: i_Hist_Panels8 = i_Hist_Panels7 + 1
   integer, parameter :: i_Hist_Panels9 = i_Hist_Panels8 + 1

   integer, parameter :: i_Col_Mag1 = i_Hist_Panels9 + 1
   integer, parameter :: i_Col_Mag2 = i_Col_Mag1 + 1
   integer, parameter :: i_Col_Mag3 = i_Col_Mag2 + 1
   integer, parameter :: i_Col_Mag4 = i_Col_Mag3 + 1
   integer, parameter :: i_Col_Mag5 = i_Col_Mag4 + 1
   integer, parameter :: i_Col_Mag6 = i_Col_Mag5 + 1
   integer, parameter :: i_Col_Mag7 = i_Col_Mag6 + 1
   integer, parameter :: i_Col_Mag8 = i_Col_Mag7 + 1
   integer, parameter :: i_Col_Mag9 = i_Col_Mag8 + 1

   integer, parameter :: i_Grid1 = i_Col_Mag9 + 1
   integer, parameter :: i_Grid2 = i_Grid1 + 1
   integer, parameter :: i_Grid3 = i_Grid2 + 1
   integer, parameter :: i_Grid4 = i_Grid3 + 1
   integer, parameter :: i_Grid5 = i_Grid4 + 1
   integer, parameter :: i_Grid6 = i_Grid5 + 1
   integer, parameter :: i_Grid7 = i_Grid6 + 1
   integer, parameter :: i_Grid8 = i_Grid7 + 1
   integer, parameter :: i_Grid9 = i_Grid8 + 1
   integer, parameter :: i_Other = i_Grid9 + 1

   integer, parameter :: num_pgstar_plots = i_Other + max_num_Other_plots

   integer, parameter :: max_num_pgstar_grid_plots = 10
   integer, parameter :: max_num_rows_Text_Summary = 20
   integer, parameter :: max_num_cols_Text_Summary = 20
   integer, parameter :: max_num_profile_mass_points = 10


   ! some Tioga colors for pgstar
   integer :: clr_Black
   integer :: clr_Blue
   integer :: clr_BrightBlue
   integer :: clr_Goldenrod
   integer :: clr_Lilac
   integer :: clr_Coral
   integer :: clr_FireBrick
   integer :: clr_RoyalPurple
   integer :: clr_Gold
   integer :: clr_Crimson
   integer :: clr_SlateGray
   integer :: clr_Teal
   integer :: clr_LightSteelBlue
   integer :: clr_MediumSlateBlue
   integer :: clr_MediumSpringGreen
   integer :: clr_MediumBlue
   integer :: clr_RoyalBlue
   integer :: clr_LightGray
   integer :: clr_Silver
   integer :: clr_DarkGray
   integer :: clr_Gray
   integer :: clr_LightSkyBlue
   integer :: clr_LightSkyGreen
   integer :: clr_SeaGreen
   integer :: clr_Tan
   integer :: clr_IndianRed
   integer :: clr_LightOliveGreen
   integer :: clr_CadetBlue
   integer :: clr_Beige

   integer, parameter :: max_num_pgstar_trace_history_values = 20


   type pgstar_controls
      include "pgstar_controls.inc"

      procedure(pgstar_decorator_interface), pointer, nopass :: &
            Abundance_pgstar_decorator => null(), &
            Color_Magnitude1_pgstar_decorator => null(), &
            Color_Magnitude2_pgstar_decorator => null(), &
            Color_Magnitude3_pgstar_decorator => null(), &
            Color_Magnitude4_pgstar_decorator => null(), &
            Color_Magnitude5_pgstar_decorator => null(), &
            Color_Magnitude6_pgstar_decorator => null(), &
            Color_Magnitude7_pgstar_decorator => null(), &
            Color_Magnitude8_pgstar_decorator => null(), &
            Color_Magnitude9_pgstar_decorator => null(), &
            dPg_dnu_pgstar_decorator => null(), &
            Dynamo_pgstar_decorator => null(), &
            History_Panels1_pgstar_decorator => null(), &
            History_Panels2_pgstar_decorator => null(), &
            History_Panels3_pgstar_decorator => null(), &
            History_Panels4_pgstar_decorator => null(), &
            History_Panels5_pgstar_decorator => null(), &
            History_Panels6_pgstar_decorator => null(), &
            History_Panels7_pgstar_decorator => null(), &
            History_Panels8_pgstar_decorator => null(), &
            History_Panels9_pgstar_decorator => null(), &
            History_Track1_pgstar_decorator => null(), &
            History_Track2_pgstar_decorator => null(), &
            History_Track3_pgstar_decorator => null(), &
            History_Track4_pgstar_decorator => null(), &
            History_Track5_pgstar_decorator => null(), &
            History_Track6_pgstar_decorator => null(), &
            History_Track7_pgstar_decorator => null(), &
            History_Track8_pgstar_decorator => null(), &
            History_Track9_pgstar_decorator => null(), &
            HR_pgstar_decorator => null(), &
            Kipp_pgstar_decorator => null(), &
            logg_logt_pgstar_decorator => null(), &
            logg_teff_pgstar_decorator => null(), &
            logl_r_pgstar_decorator => null(), &
            logl_teff_pgstar_decorator => null(), &
            logl_v_pgstar_decorator => null(), &
            l_r_pgstar_decorator => null(), &
            l_teff_pgstar_decorator => null(), &
            l_v_pgstar_decorator => null(), &
            mixing_pgstar_decorator => null(), &
            mode_prop_pgstar_decorator => null(), &
            network_pgstar_decorator => null(), &
            power_pgstar_decorator => null(), &
            production_pgstar_decorator => null(), &
            Profile_Panels1_pgstar_decorator => null(), &
            Profile_Panels2_pgstar_decorator => null(), &
            Profile_Panels3_pgstar_decorator => null(), &
            Profile_Panels4_pgstar_decorator => null(), &
            Profile_Panels5_pgstar_decorator => null(), &
            Profile_Panels6_pgstar_decorator => null(), &
            Profile_Panels7_pgstar_decorator => null(), &
            Profile_Panels8_pgstar_decorator => null(), &
            Profile_Panels9_pgstar_decorator => null(), &
            R_L_pgstar_decorator => null(), &
            R_Teff_pgstar_decorator => null(), &
            rti_pgstar_decorator => null(), &
            summary_burn_pgstar_decorator => null(), &
            summary_profile_pgstar_decorator => null(), &
            summary_history_pgstar_decorator => null(), &
            trho_pgstar_decorator => null(), &
            tmaxrho_pgstar_decorator => null(), &
            trho_profile_pgstar_decorator => null()



   end type pgstar_controls

end module star_pgstar
