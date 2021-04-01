      module test_xxx_mod
      use xxx_lib
      implicit none

      contains

      subroutine do_test
         
         write(*,*) 'done'
         
      end subroutine do_test 


      end module test_xxx_mod




      program test_xxx
      use test_xxx_mod
      implicit none
      call do_test
      end program
