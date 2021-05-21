      program test_eos
      use eos_support, only: Setup_eos
      use test_eos_support
      use math_lib
      use auto_diff
      use test_eos_blend
      
      implicit none
      
      logical, parameter :: quietly = .false.

      call Setup_eos
      
      !call test_eosDT(1) ! eos_use_FreeEOS
      !call test_eosDT(0) ! old form
      !call test_eosPT(1) ! eos_use_max_SCVH_for_PT
      !call test_eosPT(0) ! old form

      !call test_eosDT(2) ! X derivatives

      call Do_One(quietly)

      call test1_eosPT_for_ck(quietly)
      
      call do_test_eos_blend()

      end   

