      module test_gyre_mod
      implicit none

      contains

      subroutine do_test
         integer :: ierr
         character (len=1000) :: line
         ierr = 0
         open(33,file='test_output')
         do 
            read(33,fmt='(a)',iostat=ierr) line
            if (ierr /= 0) exit
            write(*,'(a)') trim(line)
         end do
         close(33)
      end subroutine do_test 


      end module test_gyre_mod




      program test_gyre
      use test_gyre_mod
      implicit none
      call do_test
      end program
