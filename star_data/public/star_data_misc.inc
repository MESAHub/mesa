! star_data_misc.inc
   ! saved in photos

      integer :: recent_log_header ! model_number when last output header info to terminal
      
      real(dp) :: tau_factor, tau_base
      real(dp) :: force_tau_factor
      real(dp) :: force_Tsurf_factor
      real(dp) :: force_opacity_factor
      
      logical :: using_revised_max_yr_dt
      real(dp) :: revised_max_yr_dt
      real(dp) :: revised_max_yr_dt_old
      
      logical :: astero_using_revised_max_yr_dt
      real(dp) :: astero_revised_max_yr_dt
      real(dp) :: astero_revised_max_yr_dt_old

      integer :: binary_id        
      logical :: include_binary_history_in_log_file
         
      logical :: have_j_rot
      
      logical :: have_previous_conv_vel ! needed for conv_vel variables
   
      ! experimental
      logical :: am_nu_rot_flag! am_nu_rot_flag is true if we are using smoothed am_nu_rot variables
      
