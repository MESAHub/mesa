      program test_eos_quietly
      use test_eos_support
      
      implicit none
      
      logical, parameter :: quietly = .true.
                  
      call Setup_eos

      call Do_One(quietly)
      call test1_eosPT_for_ck(quietly)
      call test1_eosDE_for_ck(quietly)

      end   

