      module test_integrate
            
         use num_def
         use num_lib
         use const_def   

         implicit none

         contains

         subroutine run_test_integrate()


            call test_basic()
            call test_sine()
            call test_exp()

         end subroutine run_test_integrate



         subroutine test_basic
            real(dp) :: xlow=0, xhigh=1
            real(dp) :: expected = 0.5d0
            real(dp) :: res
            integer :: ierr

            res = integrate(linear, xlow, xhigh, (/1d0/), 1d-3,1d-3, 10,ierr)

            write(*,*) 'integrate linear expected ',expected, ' got ', res, ierr

            contains


            real(dp) function linear(x, args, ierr)
               real(dp), intent(in) :: x
               real(dp), intent(in) :: args(:)
               integer, intent(inout) :: ierr
 
               ierr = 0
               linear = x

            end function linear 

         end subroutine test_basic

         subroutine test_sine
            real(dp) :: res
            integer :: ierr

            res = integrate(sine, 0d0, pi, (/1d0/), 1d-8,1d-8, 10, ierr)

            write(*,*) 'integrate sine expected ',2d0, ' got ', res, ierr

            contains


            real(dp) function sine(x, args, ierr)
               real(dp), intent(in) :: x
               real(dp), intent(in) :: args(:)
               integer, intent(inout) :: ierr
 
               ierr = 0
               sine = sin(x)

            end function sine

         end subroutine test_sine


         subroutine test_exp
            real(dp) :: res
            integer :: ierr

            res = integrate(iexp, 0d0, 2d0, (/1d0/), 1d-8,1d-8, 50, ierr)

            write(*,*) 'integrate exp expected ',exp(2d0)-1d0, ' got ', res, ierr

            contains


            real(dp) function iexp(x, args, ierr)
               real(dp), intent(in) :: x
               real(dp), intent(in) :: args(:)
               integer, intent(inout) :: ierr
 
               ierr = 0
               iexp = exp(x)

            end function iexp

         end subroutine test_exp


      end module test_integrate
