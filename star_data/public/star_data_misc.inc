! star_data_misc.inc
   ! saved in photos

 

! Q: A lot of names are "log" something -- natural or base 10 log?
! A: I use "ln" for natural log and "log" for base 10.


! Q: How is the "center" of a zone defined? is it the center by mass, or radius?
! A: Extensive variables (R, L, M, velocity) are defined at cell boundaries.
!    Intensive variables (T, rho, composition) are cell averages. 
!    We use the cell average as an estimate of the point value at the cell center by mass.
!   The themo variables are cell averages since they are calculated as functions of the cell averages.
!   Note however that gradr and gradT are not thermo variables.
!   They are defined at cell boundaries and depend on both neighbors.
!   If you are comparing to results from a model that doesn't have staggered variables
!   (i.e., one that defines the intensive and extensive vars at the same place), then
!   you will need to adjust.   For example, to compare our sound speed profile to
!   published results, we interpolate the cell centered csound values to find values
!   at cell boundaries.
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
      ! history
         integer, pointer :: profile_column_spec(:)
         integer :: number_of_history_columns
            ! number_of_history_columns = num_standard_history_columns + num extra history columns
         real(dp), pointer :: history_values(:) ! (number_of_history_columns)
         logical, pointer :: history_value_is_integer(:) ! (number_of_history_columns)
         character (len=maxlen_history_column_name), pointer :: &
            history_names(:) ! (number_of_history_columns)
         type (integer_dict), pointer :: history_names_dict
            ! maps names to index in history_values
         integer :: model_number_of_history_values
         logical :: need_to_set_history_names_etc
         ! the following does not include "extra" columns
         integer, pointer :: history_column_spec(:) ! (num_standard_history_columns)
      
      ! profile
         character (len=strlen) :: &
            model_profile_filename, model_controls_filename, model_data_filename
         character (len=strlen) :: &
            most_recent_profile_filename, most_recent_controls_filename, &
            most_recent_model_data_filename
         
      ! pgstar
         type (pgstar_hist_node), pointer :: pgstar_hist => null()
      
      ! random number generator (marsaglia and zaman). period is 2**43 with
      ! 900 million different sequences. the state of the generator (for restarts)
         integer :: rand_i97, rand_j97
         real(dp) :: rand_u(rand_u_len), rand_c, rand_cd, rand_cm
      
      real(dp) :: cumulative_energy_error, cumulative_extra_heating
      integer :: total_num_solver_iterations, total_num_solver_relax_iterations
      integer :: total_num_solver_calls_made, total_num_solver_relax_calls_made
      integer :: total_num_solver_calls_converged, total_num_solver_relax_calls_converged

      integer :: recent_log_header ! model_number when last output header info to terminal
      
      real(dp) :: tau_factor
      
      integer :: phase_of_evolution
      
      
      logical :: using_revised_max_yr_dt
      real(dp) :: revised_max_yr_dt
      real(dp) :: revised_max_yr_dt_old
      
      logical :: astero_using_revised_max_yr_dt
      real(dp) :: astero_revised_max_yr_dt
      real(dp) :: astero_revised_max_yr_dt_old

      integer :: binary_id        
      logical :: include_binary_history_in_log_file
         
      logical :: have_j_rot
   
      
