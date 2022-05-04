
! to be included in extras_controls in run_star_extras

! include 'set_star_astero_procs.inc'

! also need to add 
! use astero_def, only: star_astero_procs

   star_astero_procs% set_my_vars => set_my_vars
   star_astero_procs% set_my_param => set_my_param
   star_astero_procs% extras_startup => extras_startup
   star_astero_procs% extras_controls => extras_controls
   star_astero_procs% extras_check_model => extras_check_model
   star_astero_procs% extras_finish_step => extras_finish_step
   star_astero_procs% extras_after_evolve => extras_after_evolve
   star_astero_procs% how_many_extra_history_columns => how_many_extra_history_columns
   star_astero_procs% data_for_extra_history_columns => data_for_extra_history_columns
   star_astero_procs% how_many_extra_profile_columns => how_many_extra_profile_columns
   star_astero_procs% data_for_extra_profile_columns => data_for_extra_profile_columns
   