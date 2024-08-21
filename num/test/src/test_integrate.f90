      module test_integrate

         use math_lib
         use num_def
         use num_lib
         use const_def

         implicit none

         contains

         subroutine run_test_integrate()


            call test_basic()
            call test_sine()
            call test_exp()
            call test_box()

         end subroutine run_test_integrate



         subroutine test_basic
            real(dp) :: xlow=0, xhigh=1
            real(dp) :: expected = 0.5d0
            real(dp) :: res
            integer :: ierr

            res = integrate(linear, xlow, xhigh, (/1d0/), 1d-3, 1d-3, 10, ierr)

            call check_result('linear', expected, res, ierr)

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

            res = integrate(sine, 0d0, pi, (/1d0/), 1d-8, 1d-8, 10, ierr)

            call check_result('sine', 2d0, res, ierr)

            res = integrate(sine, 0d0, 2*pi, (/1d0/), 1d-8, 1d-8, 10, ierr)

            call check_result('sine', 0d0, res, ierr)

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

            res = integrate(iexp, 0d0, 2d0, (/1d0/), 1d-8, 1d-8, 50, ierr)

            call check_result('exp', exp(2d0)-1d0, res, ierr)

            res = integrate(iexp, 0d0, 10d0, (/1d0/), 1d-8, 1d-8, 50, ierr)

            call check_result('exp', exp(10d0)-1d0, res, ierr)

            contains


            real(dp) function iexp(x, args, ierr)
               real(dp), intent(in) :: x
               real(dp), intent(in) :: args(:)
               integer, intent(inout) :: ierr
 
               ierr = 0
               iexp = exp(x)

            end function iexp

         end subroutine test_exp


         subroutine test_box
            real(dp) :: res
            integer :: ierr

            res = integrate(box, 0d0, 2d0, (/1d0/), 1d-8, 1d-8, 50, ierr)

            call check_result('box', 1d0, res, ierr)

            res = integrate(box, 0.99d0, 1.5d0, (/1d0/), 1d-8, 1d-8, 50, ierr)

            call check_result('box', 0.5d0, res, ierr)

            contains


            real(dp) function box(x, args, ierr)
               real(dp), intent(in) :: x
               real(dp), intent(in) :: args(:)
               integer, intent(inout) :: ierr
 
               ierr = 0
         
               if(x<1) then
                  box = 0d0
               else if(x.ge.1d0 .and. x.le.2d0) then
                  box = 1d0
               else
                  box = 0d0
               end if

            end function box

         end subroutine test_box


         subroutine check_result(name, tgt, val, ierr)
            character(len=*), intent(in) :: name
            real(dp), intent(in) :: tgt, val
            integer, intent(in) :: ierr

            write(*, '(a40, 1pd26.16, a7, 1pd26.16, i4)') &
                 'integrate '// trim(name) // ' expected', &
                 tgt, 'got', val, ierr

         end subroutine check_result

      end module test_integrate
