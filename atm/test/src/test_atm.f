      program test_atm
      use const_def, only: dp
      use test_atm_setup, only: setup
      use test_atm_support, only: do_test_atm, &
         test_verbosely, cgrav, eos_handle, kap_handle
      
      logical :: test_verbosely_in
      real(dp) :: cgrav_in
      integer :: eos_handle_in, kap_handle_in
      
      call setup
      test_verbosely = .true.
      test_verbosely_in = test_verbosely
      cgrav_in = cgrav
      eos_handle_in = eos_handle
      kap_handle_in = kap_handle
      
      call do_test_atm( &
         test_verbosely, cgrav, eos_handle, kap_handle)
         
      end program




