      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! binary_job controls

      logical :: show_binary_log_description_at_start
      character (len=1000) :: binary_history_columns_file
      logical :: warn_binary_extra
      character (len=256), dimension(2) :: inlist_names
      logical :: read_extra_binary_job_inlist1
      character (len=256) :: extra_binary_job_inlist1_name 
      logical :: read_extra_binary_job_inlist2
      character (len=256) :: extra_binary_job_inlist2_name 
      logical :: read_extra_binary_job_inlist3
      character (len=256) :: extra_binary_job_inlist3_name 
      logical :: read_extra_binary_job_inlist4
      character (len=256) :: extra_binary_job_inlist4_name 
      logical :: read_extra_binary_job_inlist5
      character (len=256) :: extra_binary_job_inlist5_name 

      logical :: evolve_both_stars
      logical :: relax_primary_to_th_eq
      real(dp) :: log_Lnuc_div_L_for_relax_primary_to_th_eq
      real(dp) :: min_age_for_relax_primary_to_th_eq
      integer :: max_steps_for_relax_primary_to_th_eq
      logical :: no_history_during_relax_primary_to_th_eq
      logical :: reset_age_for_relax_primary_to_th_eq
      real(dp) :: tsync_for_relax_primary_to_th_eq

      logical :: change_ignore_rlof_flag
      logical :: change_initial_ignore_rlof_flag
      logical :: new_ignore_rlof_flag
      logical :: change_model_twins_flag
      logical :: change_initial_model_twins_flag
      logical :: new_model_twins_flag
      logical :: change_point_mass_i
      logical :: change_initial_point_mass_i
      integer :: new_point_mass_i
      logical :: change_m1
      logical :: change_initial_m1
      real(dp) :: new_m1
      logical :: change_m2
      logical :: change_initial_m2
      real(dp) :: new_m2
      logical :: change_separation_eccentricity
      logical :: change_initial_separation_eccentricity
      logical :: change_period_eccentricity
      logical :: change_initial_period_eccentricity
      real(dp) :: new_separation
      real(dp) :: new_period
      real(dp) :: new_eccentricity
      logical :: pgbinary_flag
