      module test_astero_mod
      use astero_lib
      implicit none

      contains

      subroutine do_test
         
         write(*,*) 'done'
         
      end subroutine do_test 


      end module test_astero_mod




      program test_astero
      use test_astero_mod
      implicit none
      call do_test
      end program
