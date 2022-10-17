integer, parameter :: max_num_stars = 10 ! make as large as needed
integer :: star_ids(max_num_stars)
logical :: star_done(max_num_stars)
integer :: num_stars
logical :: restart_flag
character (len = 256), dimension(max_num_stars) :: &
      inlist_names, restart_names
integer :: which_for_pgstar
real(dp) :: stopping_age
logical :: stop_when_1st_finishes

namelist /multi_stars_job/ &
      num_stars, & ! up to max_num_stars
      inlist_names, &
      restart_flag, & ! .false.   .true. means doing a restart
      restart_names, &
      which_for_pgstar, & ! -1 negative means pgstar for all stars
      stopping_age, & ! 1d99
      stop_when_1st_finishes ! .true.
