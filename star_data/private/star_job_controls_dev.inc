      logical :: change_TDC_flag, change_initial_TDC_flag, new_TDC_flag
      integer :: change_TDC_flag_at_model_number
      
      logical :: change_conv_vel_flag, change_initial_conv_vel_flag, new_conv_vel_flag
      real(dp) :: logT_for_conv_vel_flag
      
      logical :: use_w_div_wc_flag_with_rotation, use_j_rot_flag_with_rotation
      logical :: change_w_div_wc_flag, change_initial_w_div_wc_flag, new_w_div_wc_flag
      
      logical :: change_j_rot_flag, change_initial_j_rot_flag, new_j_rot_flag
      
      logical :: remove_surface_do_entropy, remove_surface_do_jrot, &
         remove_surface_turn_off_energy_sources_and_sinks, &
         remove_surface_by_relax_to_star_cut


      logical :: load_saved_photo
      character (len=strlen) :: saved_photo_name
      
      integer :: save_photo_number
      character (len=strlen) :: save_photo_filename
      !save_photo_when_terminate = .true.
      