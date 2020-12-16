      module test_binary_mod
      use binary_lib
      implicit none

      contains

      subroutine do_test
         
         write(*,*) 'done'
         
      end subroutine do_test 


      end module test_binary_mod




      program test_binary
      use test_binary_mod
      implicit none
      call do_test
      end program
