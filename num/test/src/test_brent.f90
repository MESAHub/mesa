module test_brent

   use num_def
   use num_lib
   use math_lib
   use utils_lib, only: mesa_error
   use const_def, only: dp

   implicit none

   logical, parameter :: dbg = .false.

contains

   subroutine do_test_brent
      write (*, *) 'test brent'

      call test_global_min_all
      call test_local_min_all
      call test_brent_zero

   end subroutine do_test_brent

   subroutine test_global_min_all

      !*****************************************************************************80
      !
      !! TEST_GLOMIN_ALL tests the Brent global minimizer on all test functions.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    12 April 2008
      !
      !  Author:
      !
      !    John Burkardt
      !
      implicit none

      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: e
      real(dp) :: m
      real(dp) :: machep
      real(dp) :: t

      write (*, '(a)') ' '
      write (*, '(a)') 'TEST_GLOMIN_ALL'
      write (*, '(a)') '  Test the Brent GLOMIN routine, which seeks'
      write (*, '(a)') '  a global minimizer of a function F(X)'
      write (*, '(a)') '  in an interval [A,B],'
      write (*, '(a)') '  given some upper bound M for F".'

      machep = epsilon(machep)
      e = sqrt(machep)
      t = sqrt(machep)

      a = 7.0D+00
      b = 9.0D+00
      c = (a + b)/2.0D+00
      m = 0.0D+00

      call test_glomin_one(a, b, c, m, machep, e, t, h_01, &
                           'h_01(x) = 2 - x')

      a = 7.0D+00
      b = 9.0D+00
      c = (a + b)/2.0D+00
      m = 100.0D+00

      call test_glomin_one(a, b, c, m, machep, e, t, h_01, &
                           'h_01(x) = 2 - x')

      a = -1.0D+00
      b = +2.0D+00
      c = (a + b)/2.0D+00
      m = 2.0D+00

      call test_glomin_one(a, b, c, m, machep, e, t, h_02, &
                           'h_02(x) = x * x')

      a = -1.0D+00
      b = +2.0D+00
      c = (a + b)/2.0D+00
      m = 2.1D+00

      call test_glomin_one(a, b, c, m, machep, e, t, h_02, &
                           'h_02(x) = x * x')

      a = -0.5D+00
      b = +2.0D+00
      c = (a + b)/2.0D+00
      m = 14.0D+00

      !call test_glomin_one ( a, b, c, m, machep, e, t, h_03, &
      !  'h_03(x) = x^3 + x^2' )

      a = -0.5D+00
      b = +2.0D+00
      c = (a + b)/2.0D+00
      m = 28.0D+00

      !call test_glomin_one ( a, b, c, m, machep, e, t, h_03, &
      !  'h_03(x) = x^3 + x^2' )

      a = -10.0D+00
      b = +10.0D+00
      c = (a + b)/2.0D+00
      m = 72.0D+00

      call test_glomin_one(a, b, c, m, machep, e, t, h_04, &
                           'h_04(x) = ( x + sin(x) ) * exp(-x*x)')

      a = -10.0D+00
      b = +10.0D+00
      c = (a + b)/2.0D+00
      m = 72.0D+00

      call test_glomin_one(a, b, c, m, machep, e, t, h_05, &
                           'h_05(x) = ( x - sin(x) ) * exp(-x*x)')

   end subroutine test_global_min_all

   subroutine test_glomin_one(a, b, c, m, machep, e, t, f, title)
      real(dp), intent(in) :: a, b, c, m, machep, e, t
      interface
         real(dp) function f(x)
            use const_def, only: dp
            implicit none
            real(dp), intent(in) :: x
         end function f
      end interface

      real(dp) :: fa
      real(dp) :: fb
      real(dp) :: fx
      character(len=*) :: title
      real(dp) :: x
      integer :: max_tries, ierr
      include 'formats'

      max_tries = 10000
      ierr = 0
      fx = brent_global_min(max_tries, a, b, c, m, machep, e, t, f, x, ierr)
      if (ierr /= 0) then
         write (*, *) 'failed in brent_global_min'
         call mesa_error(__FILE__, __LINE__)
      end if
      fa = f(a)
      fb = f(b)

      write (*, '(a)') ' '
      write (*, '(a)') trim(title)
      write (*, 1) 'a,  b', a, b
      write (*, 1) 'fa, fb', fa, fb

   end subroutine test_glomin_one

   real(dp) function h_01(x)
      real(dp), intent(in) :: x
      h_01 = 2.0D+00 - x
   end function h_01

   real(dp) function h_02(x)
      real(dp), intent(in) :: x
      h_02 = x*x
   end function h_02

   real(dp) function h_03(x)
      real(dp), intent(in) :: x
      h_03 = x*x*(x + 1.0D+00)
   end function h_03

   real(dp) function h_04(x)
      real(dp), intent(in) :: x
      h_04 = (x + sin(x))*exp(-x*x)
   end function h_04

   real(dp) function h_05(x)
      real(dp), intent(in) :: x
      h_05 = (x - sin(x))*exp(-x*x)
   end function h_05

   subroutine test_local_min_all

      !*****************************************************************************80
      !
      !! TEST_LOCAL_MIN_ALL tests Brent's local minimizer on all test functions.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    12 April 2008
      !
      !  Author:
      !
      !    John Burkardt
      !
      implicit none

      real(dp) :: a
      real(dp) :: b
      real(dp) :: eps
      real(dp) :: t

      write (*, '(a)') ' '
      write (*, '(a)') 'TEST_LOCAL_MIN_ALL'
      write (*, '(a)') '  Test the Brent LOCAL_MIN routine, which seeks'
      write (*, '(a)') '  a local minimizer of a function F(X)'
      write (*, '(a)') '  in an interval [A,B].'

      eps = 10.0D+00*sqrt(epsilon(eps))
      t = eps

      a = 0.0D+00
      b = 3.141592653589793D+00

      call test_local_min_one(a, b, eps, t, g_01, &
                              'g_01(x) = ( x - 2 ) * ( x - 2 ) + 1')

      a = 0.0D+00
      b = 1.0D+00

      call test_local_min_one(a, b, eps, t, g_02, &
                              'g_02(x) = x * x + exp( - x )')

      a = -2.0D+00
      b = 2.0D+00

      call test_local_min_one(a, b, eps, t, g_03, &
                              'g_03(x) = x^4 + 2x^2 + x + 3')

      a = 0.0001D+00
      b = 1.0D+00

      call test_local_min_one(a, b, eps, t, g_04, &
                              'g_04(x) = exp( x ) + 1 / ( 100 x )')

      a = 0.0002D+00
      b = 2.0D+00

      call test_local_min_one(a, b, eps, t, g_05, &
                              'g_05(x) = exp( x ) - 2x + 1/(100x) - 1/(1000000x^2)')

   end subroutine test_local_min_all

   subroutine test_local_min_one(a, b, eps, t, f, title)

      !*****************************************************************************80
      !
      !! TEST_LOCAL_MIN_ONE tests Brent's local minimizer on one test function.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    12 April 2008
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, real(dp) A, B, the endpoints of the interval.
      !
      !    Input, real(dp) EPS, a positive relative error tolerance.
      !
      !    Input, real(dp) T, a positive absolute error tolerance.
      !
      !    Input, external real(dp) F, the name of a user-supplied
      !    function, of the form "FUNCTION F ( X )", which evaluates the
      !    function whose local minimum is being sought.
      !
      !    Input, character ( LEN = * ) TITLE, a title for the problem.
      !
      implicit none
      real(dp), intent(in) :: a, b, eps, t
      interface
         real(dp) function f(x)
            use const_def, only: dp
            implicit none
            real(dp), intent(in) :: x
         end function f
      end interface
      character(len=*) :: title

      real(dp) :: fa
      real(dp) :: fb
      real(dp) :: fx
      real(dp) :: x
      integer :: max_tries, ierr
      include 'formats'

      max_tries = 10000
      ierr = 0
      fx = brent_local_min(max_tries, a, b, eps, t, f, x, ierr)
      if (ierr /= 0) then
         write (*, *) 'failed in brent_local_min'
         call mesa_error(__FILE__, __LINE__)
      end if
      fa = f(a)
      fb = f(b)

      write (*, '(a)') ' '
      write (*, '(a)') trim(title)
      write (*, 1) 'a,  b', a, b
      write (*, 1) 'fa, fb', fa, fb

   end subroutine test_local_min_one

   real(dp) function g_01(x)
      real(dp), intent(in) :: x
      g_01 = (x - 2.0D+00)*(x - 2.0D+00) + 1.0D+00
   end function g_01

   real(dp) function g_02(x)
      real(dp), intent(in) :: x
      g_02 = x*x + exp(-x)
   end function g_02

   real(dp) function g_03(x)
      real(dp), intent(in) :: x
      g_03 = ((x*x + 2.0D+00)*x + 1.0D+00)*x + 3.0D+00
   end function g_03

   real(dp) function g_04(x)
      real(dp), intent(in) :: x
      g_04 = exp(x) + 0.01D+00/x
   end function g_04

   real(dp) function g_05(x)
      real(dp), intent(in) :: x
      g_05 = exp(x) - 2.0D+00*x + 0.01D+00/x - 0.000001D+00/x/x
   end function g_05

   real(dp) function f_01(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
      integer, intent(in) :: lrpar, lipar
      real(dp), intent(in) :: x
      real(dp), intent(out) :: dfdx
      integer, intent(inout), pointer :: ipar(:)  ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:)  ! (lrpar)
      integer, intent(out) :: ierr
      f_01 = sin(x) - 0.5D+00*x
      ierr = 0
      dfdx = 0
   end function f_01

   real(dp) function f_02(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
      integer, intent(in) :: lrpar, lipar
      real(dp), intent(in) :: x
      real(dp), intent(out) :: dfdx
      integer, intent(inout), pointer :: ipar(:)  ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:)  ! (lrpar)
      integer, intent(out) :: ierr
      f_02 = 2.0D+00*x - exp(-x)
      ierr = 0
      dfdx = 0
   end function f_02

   real(dp) function f_03(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
      integer, intent(in) :: lrpar, lipar
      real(dp), intent(in) :: x
      real(dp), intent(out) :: dfdx
      integer, intent(inout), pointer :: ipar(:)  ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:)  ! (lrpar)
      integer, intent(out) :: ierr
      f_03 = x*exp(-x)
      ierr = 0
      dfdx = 0
   end function f_03

   real(dp) function f_04(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
      integer, intent(in) :: lrpar, lipar
      real(dp), intent(in) :: x
      real(dp), intent(out) :: dfdx
      integer, intent(inout), pointer :: ipar(:)  ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:)  ! (lrpar)
      integer, intent(out) :: ierr
      f_04 = exp(x) - 1.0D+00/100.0D+00/x/x
      ierr = 0
      dfdx = 0
   end function f_04

   real(dp) function f_05(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
      integer, intent(in) :: lrpar, lipar
      real(dp), intent(in) :: x
      real(dp), intent(out) :: dfdx
      integer, intent(inout), pointer :: ipar(:)  ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:)  ! (lrpar)
      integer, intent(out) :: ierr
      f_05 = (x + 3.0D+00)*(x - 1.0D+00)*(x - 1.0D+00)
      ierr = 0
      dfdx = 0
   end function f_05

   subroutine test_brent_zero()

      !*****************************************************************************80
      !
      !! test_brent_zero tests Brent's zero finding routine on all test functions.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    12 April 2008
      !
      !  Author:
      !
      !    John Burkardt
      !
      implicit none

      real(dp) :: a
      real(dp) :: b
      real(dp) :: machep
      real(dp) :: t

      machep = epsilon(machep)
      t = machep

      a = 1.0D+00
      b = 2.0D+00

      call test_zero_one(a, b, machep, t, f_01, &
                         'f_01(x) = sin ( x ) - x / 2')

      a = 0.0D+00
      b = 1.0D+00

      call test_zero_one(a, b, machep, t, f_02, &
                         'f_02(x) = 2 * x - exp( - x )')

      a = -1.0D+00
      b = 0.5D+00

      call test_zero_one(a, b, machep, t, f_03, &
                         'f_03(x) = x * exp( - x )')

      a = 0.0001D+00
      b = 20.0D+00

      call test_zero_one(a, b, machep, t, f_04, &
                         'f_04(x) = exp( x ) - 1 / ( 100 * x * x )')

      a = -5.0D+00
      b = 2.0D+00

      call test_zero_one(a, b, machep, t, f_05, &
                         'f_05(x) = (x+3) * (x-1) * (x-1)')

   end subroutine test_brent_zero

   subroutine test_zero_one(a, b, machep, t, f, title)

      !*****************************************************************************80
      !
      !! TEST_ZERO_ONE tests Brent's zero finding routine on one test function.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    12 April 2008
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, real(dp) A, B, the two endpoints of the change of sign
      !    interval.
      !
      !    Input, real(dp) MACHEP, an estimate for the relative machine
      !    precision.
      !
      !    Input, real(dp) T, a positive error tolerance.
      !
      !    Input, external real(dp) F, the name of a user-supplied
      !    function, of the form "FUNCTION F ( X )", which evaluates the
      !    function whose zero is being sought.
      !
      !    Input, character ( LEN = * ) TITLE, a title for the problem.
      !
      implicit none

      interface
         include 'num_root_fcn.dek'  ! f provides function values
      end interface

      real(dp) :: a
      real(dp) :: b
      real(dp) :: fa
      real(dp) :: fb
      real(dp) :: fz
      real(dp) :: machep
      real(dp) :: t
      character(len=*) :: title
      real(dp) :: z
      real(dp) :: dfdx

      integer, parameter :: lrpar = 0, lipar = 0
      integer :: ierr
      real(dp), target :: rpar_ary(lrpar)
      integer, target :: ipar_ary(lipar)
      real(dp), pointer :: rpar(:)
      integer, pointer :: ipar(:)

      include 'formats'

      rpar => rpar_ary
      ipar => ipar_ary

      ierr = 0
      fa = f(a, dfdx, lrpar, rpar, lipar, ipar, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      fb = f(b, dfdx, lrpar, rpar, lipar, ipar, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      z = brent_safe_zero(a, b, machep, t, 0d0, f, fa, fb, lrpar, rpar, lipar, ipar, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      fz = f(z, dfdx, lrpar, rpar, lipar, ipar, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

      if (abs(fz) < 1d-14) return

      write (*, '(a)') ' '
      write (*, '(a)') trim(title)
      write (*, 1) 'a,  z,  b', a, z, b
      write (*, 1) 'fa, fz, fb', fa, fz, fb
      write (*, '(a)') ' '

   end subroutine test_zero_one

end module test_brent
