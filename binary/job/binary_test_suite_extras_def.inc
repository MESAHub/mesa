integer :: test_suite_time0, test_suite_time1, test_suite_clock_rate
real(dp) :: test_suite_runtime
logical :: test_suite_restart_flag

! testhub information
integer, parameter :: NUM_TESTHUB_EXTRAS = 10
character (len=80) :: testhub_extras_names(NUM_TESTHUB_EXTRAS)
real(dp) :: testhub_extras_vals(NUM_TESTHUB_EXTRAS)
