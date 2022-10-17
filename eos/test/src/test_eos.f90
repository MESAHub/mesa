program test_eos
   use eos_support, only : Setup_eos
   use test_eos_support
   use math_lib
   use auto_diff
   use test_eos_blend
   
   implicit none
   
   logical, parameter :: quietly = .false.
   
   call Setup_eos
   
   call Do_One(quietly)
   
   call test1_eosPT_for_ck(quietly)
   
   call do_test_eos_blend()

end

