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
      
      real(dp) :: tau_factor, tau_base
      real(dp) :: force_tau_factor
      real(dp) :: force_Tsurf_factor
      real(dp) :: force_opacity_factor
      
      
      logical :: using_revised_max_yr_dt
      real(dp) :: revised_max_yr_dt
      
      logical :: astero_using_revised_max_yr_dt
      real(dp) :: astero_revised_max_yr_dt

      integer :: binary_id        
      logical :: include_binary_history_in_log_file
         
      logical :: have_j_rot
      
      logical :: have_previous_conv_vel ! only needed when set conv_vel_flag true
   
      ! experimental "dev" items
         integer :: i_dchem1 ! index of the first delta chem variable (0 means no delta chem variables)
         integer :: equdchem1 ! index of the first delta chem equation (0 means no delta chem equations)
