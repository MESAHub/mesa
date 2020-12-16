      program test_atm_quietly
      use const_def, only: dp
      use test_atm_setup, only: setup
      use test_atm_support, only: do_test_atm, &
         test_verbosely, cgrav, eos_handle, kap_handle
      
      call setup
      test_verbosely = .false.
      
      call do_test_atm( &
         test_verbosely, cgrav, eos_handle, kap_handle)
         
      end program




