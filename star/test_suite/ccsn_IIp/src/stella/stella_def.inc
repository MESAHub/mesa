      ! declarations for stella
         
      integer :: save_stella_data_for_model_number
      logical :: save_stella_data_when_terminate, stella_show_headers
      character (len=256) :: save_stella_data_filename
      integer :: stella_num_points, stella_nz_extra
      real(dp) :: stella_min_surf_logRho, stella_min_velocity, stella_skip_inner_dm, &
         stella_skip_inner_v_limit, &
         stella_mdot_years_for_wind, stella_mdot_for_wind, stella_v_wind
         
      namelist /stella/ &
         save_stella_data_for_model_number, &
         save_stella_data_when_terminate, &
         save_stella_data_filename, &
         stella_min_surf_logRho, &
         stella_min_velocity, &
         stella_skip_inner_v_limit, &
         stella_skip_inner_dm, &
         stella_num_points, &
         stella_nz_extra, &
         stella_mdot_years_for_wind, &
         stella_mdot_for_wind, &
         stella_v_wind, &
         stella_show_headers
            
      ! end of declarations for stella
