      module test_simplex
      
      use num_def
      use num_lib
      use math_lib
      
      implicit none
      
      integer :: num_calls
      
      logical, parameter :: show_details = .false.
      

      contains
      
      
      subroutine do_test_simplex  
         !call test_FR2 ! okay -- escapes from local min
         call test_FR4 ! okay -- escapes from local min
         !call test_FR6 ! okay -- escapes from local min
         call test_ER ! okay
         call test_WD ! okay
         call test_BLE ! okay
         call test_PS ! okay
         call test_TR ! okay
      end subroutine do_test_simplex
      

      subroutine test1_simplex( &
            n, x_first, x_lower, x_upper, simplex, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, fcn, str)
         integer, intent(in) :: n ! number of dimensions
         real(dp), dimension(:) :: x_first, x_lower, x_upper
         real(dp) :: simplex(:,:), centroid_weight_power
         logical, intent(in) :: enforce_bounds, adaptive_random_search
         interface
            include 'num_simplex_fcn.dek'
         end interface
         character (len=*) :: str

         real(dp) :: f(n+1), x_final(n), x_atol, x_rtol, f_final
         integer :: iter_max, fcn_calls_max
         integer :: lrpar, lipar
         logical :: start_from_given_simplex_and_f
         integer, pointer :: ipar(:) ! (lipar)
         real(dp), pointer :: rpar(:) ! (lrpar)
         integer :: num_iters, num_fcn_calls, &
            num_fcn_calls_for_ars, num_accepted_for_ars, ierr
         integer :: seed, i, j, k
         real(dp) :: alpha, beta, gamma, delta
         
         include 'formats'
         
         write(*,*) 'testing NM_simplex with ' // trim(str)
         
         num_calls = 0
         lrpar = 0; lipar = 0
         allocate(rpar(lrpar), ipar(lipar))
         
         x_atol = 1d-10
         x_rtol = 1d-10
         
         iter_max = 1000
         fcn_calls_max = iter_max*10
         seed = 1074698122
         
         start_from_given_simplex_and_f = .false.

         alpha = 1d0
         beta = 2d0
         gamma = 0.5d0
         delta = 0.5d0
              
         call NM_simplex( &
            n, x_lower, x_upper, x_first, x_final, f_final, &
            simplex, f, start_from_given_simplex_and_f, &
            fcn, x_atol, x_rtol, &
            iter_max, fcn_calls_max, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, seed, &
            alpha, beta, gamma, delta, &
            lrpar, rpar, lipar, ipar, &
            num_iters, num_fcn_calls, &
            num_fcn_calls_for_ars, num_accepted_for_ars, ierr)
         
         if (ierr /= 0) then
            write(*,*) 'failed in do_simplex'
         else
            if (f_final < 1d-10) then
               write(*,1) 'found f_final < 1d-10'
            else
               write(*,1) 'failed in do_simplex:  f_final too large', f_final
            end if
            if (show_details) then
               write(*,1) 'f_final', f_final
               write(*,1) 'x_final', x_final(1:n)
               write(*,2) 'num_iters', num_iters
               write(*,2) 'num_fcn_calls', num_fcn_calls
               if (num_fcn_calls_for_ars > 0) &
                  write(*,2) 'num_fcn_calls_for_ars', num_fcn_calls_for_ars
               if (num_accepted_for_ars > 0) &
                  write(*,2) 'num_accepted_for_ars', num_accepted_for_ars
            end if
         end if
         write(*,*)

         deallocate(rpar, ipar)

      end subroutine test1_simplex
      
      
      subroutine test_WD
         ! testing with 4 dimensional Wood function
         integer, parameter :: n = 4
         real(dp), dimension(n) :: x_first, x_lower, x_upper
         real(dp) :: simplex(n,n+1), centroid_weight_power
         logical :: enforce_bounds, adaptive_random_search
         
         x_first(1:n) = (/ -3d0, -1d0, -3d0, -1d0 /)
         x_lower(1:n) = (/ -4d0, -2d0, -4d0, -2d0 /)
         x_upper(1:n) = (/ 2d0, 2d0, 2d0, 2d0 /)

         enforce_bounds = .true.
         adaptive_random_search = .true.
         centroid_weight_power = 1d0
         
         call test1_simplex( &
            n, x_first, x_lower, x_upper, simplex, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, fcn_WD, 'WD')
         
      end subroutine test_WD


      real(dp) function fcn_WD(n, x, lrpar, rpar, lipar, ipar, op_code, ierr)
         use const_def, only: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(:) ! (n)
         integer, intent(in) :: lrpar, lipar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: op_code
         integer, intent(out) :: ierr
         ierr = 0
         fcn_WD = WD(x)
      end function fcn_WD

      
      real(dp) function WD ( x ) ! Wood function
         real(dp), intent(in) :: x(:)
         integer :: i, n
         n = size(x,dim=1)
         WD = 0
         do i = 1, n/4
            WD = WD + &
               100d0*pow2(x(i+1) - x(i)*x(i)) + pow2(1d0 - x(i)) + &
               90d0*pow2(x(i+3) - x(i+2)*x(i+2)) + pow2(1d0 - x(i+2)) + &
               100d0*pow2(x(i+1) + x(i+3) - 2d0) + pow2(x(i+1) - x(i+3))/sqrt(10d0)
         end do
         num_calls = num_calls + 1
      end function WD
      
      
      subroutine test_ER
         ! testing with 4 dimensional extended Rosenbrock
         integer, parameter :: n = 4
         real(dp), dimension(n) :: x_first, x_lower, x_upper
         real(dp) :: simplex(n,n+1), centroid_weight_power
         logical :: enforce_bounds, adaptive_random_search
         
         x_first(1:n) = 3d0
         x_lower(1:n) = -2d0
         x_upper(1:n) = 5d0

         enforce_bounds = .true.
         adaptive_random_search = .true.
         centroid_weight_power = 1d0
         
         call test1_simplex( &
            n, x_first, x_lower, x_upper, simplex, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, fcn_ER, 'ER')
         
      end subroutine test_ER


      real(dp) function fcn_ER(n, x, lrpar, rpar, lipar, ipar, op_code, ierr)
         use const_def, only: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(:) ! (n)
         integer, intent(in) :: lrpar, lipar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: op_code
         integer, intent(out) :: ierr
         ierr = 0
         fcn_ER = ER(x)
      end function fcn_ER

      
      real(dp) function ER ( x ) ! extended Rosenbrock
         real(dp), intent(in) :: x(:)
         integer :: i, n
         n = size(x,dim=1)
         ER = 0
         do i = 1, n-1
            ER = ER + 100d0*pow2(x(i+1) - x(i)*x(i)) + pow2(1D0 - x(i))
         end do
         num_calls = num_calls + 1
      end function ER
      
      
      subroutine test_FR2
         ! testing with 2 dimensional Freudenstein and Roth function
         integer, parameter :: n = 2
         real(dp), dimension(n) :: x_first, x_lower, x_upper
         real(dp) :: simplex(n,n+1), centroid_weight_power
         logical :: enforce_bounds, adaptive_random_search
         
         ! NOTE --- this function has a local minimum
         ! and the starting values lead to the false-minimum first.
         
         ! true min = 0 at (5,4)
         ! false min = 48.98... at (11.41..., -0.8968...)
         ! starting at (0.5, -2) leads to the local min.
         
         x_first(1:n) = (/ 0.5d0, -2d0 /)
         x_lower(1:n) = -2d0
         x_upper(1:n) = 6d0

         enforce_bounds = .false.
         centroid_weight_power = 1d0
         adaptive_random_search = .true.
         
         call test1_simplex( &
            n, x_first, x_lower, x_upper, simplex, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, fcn_FR, 'FR2')
         
      end subroutine test_FR2
      
      
      subroutine test_FR4
         ! testing with 4 dimensional Freudenstein and Roth function
         integer, parameter :: n = 4
         real(dp), dimension(n) :: x_first, x_lower, x_upper
         real(dp) :: simplex(n,n+1), centroid_weight_power
         logical :: enforce_bounds, adaptive_random_search
            
         ! = 0 at (5,4)
         ! = 48.98... at (11.41..., -0.8968...)
         ! starting at (0.5, -2) leads to the bad local min.
         
         x_first(1:n) = (/ 0.5d0, -2d0, 0.5d0, -2d0 /)
         x_lower(1:n) = -2d0
         x_upper(1:n) = 6d0

         enforce_bounds = .false.
         centroid_weight_power = 1d0
         adaptive_random_search = .true.
         
         call test1_simplex( &
            n, x_first, x_lower, x_upper, simplex, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, fcn_FR, 'FR4')
         
      end subroutine test_FR4
      
      
      subroutine test_FR6
         ! testing with 6 dimensional Freudenstein and Roth function
         integer, parameter :: n = 6
         real(dp), dimension(n) :: x_first, x_lower, x_upper
         real(dp) :: simplex(n,n+1), centroid_weight_power
         logical :: enforce_bounds, adaptive_random_search
            
         ! = 0 at (5,4)
         ! = 48.98... at (11.41..., -0.8968...)
         ! starting at (0.5, -2) leads to the bad local min.
         
         x_first(1:n) = (/ 0.5d0, -2d0, 0.5d0, -2d0, 0.5d0, -2d0 /)
         x_lower(1:n) = -2d0
         x_upper(1:n) = 6d0

         enforce_bounds = .false.
         centroid_weight_power = 1d0
         adaptive_random_search = .true.
         
         call test1_simplex( &
            n, x_first, x_lower, x_upper, simplex, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, fcn_FR, 'FR6')
         
      end subroutine test_FR6


      real(dp) function fcn_FR(n, x, lrpar, rpar, lipar, ipar, op_code, ierr)
         use const_def, only: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(:) ! (n)
         integer, intent(in) :: lrpar, lipar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: op_code
         integer, intent(out) :: ierr
         ierr = 0
         fcn_FR = FR(x)
      end function fcn_FR

      
      real(dp) function FR ( x ) ! Freudenstein and Roth function
         real(dp), intent(in) :: x(:)
         integer :: i, n
         n = size(x,dim=1)
         FR = 0
         do i = 1, n/2
            FR = FR + &
               pow2(-13d0 + x(2*i-1) + ((5d0 - x(2*i))*x(2*i) - 2d0)*x(2*i)) + &
               pow2(-29d0 + x(2*i-1) + ((1d0 + x(2*i))*x(2*i) - 14d0)*x(2*i))
         end do
         num_calls = num_calls + 1
      end function FR
      

      subroutine test_BLE
         ! testing with 4 dimensional Beale function
         integer, parameter :: n = 4
         real(dp), dimension(n) :: x_first, x_lower, x_upper
         real(dp) :: simplex(n,n+1), centroid_weight_power
         logical :: enforce_bounds, adaptive_random_search
         
         x_first(1:n) = 3d0
         x_lower(1:n) = -2d0
         x_upper(1:n) = 5d0

         enforce_bounds = .true.
         adaptive_random_search = .true.
         centroid_weight_power = 1d0
         
         call test1_simplex( &
            n, x_first, x_lower, x_upper, simplex, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, fcn_BLE, 'BLE')
         
      end subroutine test_BLE


      real(dp) function fcn_BLE(n, x, lrpar, rpar, lipar, ipar, op_code, ierr)
         use const_def, only: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(:) ! (n)
         integer, intent(in) :: lrpar, lipar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: op_code
         integer, intent(out) :: ierr
         ierr = 0
         fcn_BLE = BLE(x)
      end function fcn_BLE

      
      real(dp) function BLE ( x ) ! Beale function
         real(dp), intent(in) :: x(:)
         integer :: i, n
         n = size(x,dim=1)
         BLE = 0
         do i = 1, n/2
            BLE = BLE + &
               pow2(1.5d0 - x(2*i-1)*(1d0 - x(2*i))) + &
               pow2(2.25d0 - x(2*i-1)*(1d0 - x(2*i)*x(2*i))) + &
               pow2(2.625d0 - x(2*i-1)*(1d0 - x(2*i)*x(2*i)*x(2*i)))
         end do
         num_calls = num_calls + 1
      end function BLE
      
      
      subroutine test_PS
         ! testing with 4 dimensional Powell singular function
         integer, parameter :: n = 4
         real(dp), dimension(n) :: x_first, x_lower, x_upper
         real(dp) :: simplex(n,n+1), centroid_weight_power
         logical :: enforce_bounds, adaptive_random_search
         
         x_first(1:n) = (/ 3d0, -1d0, 0d0, 1d0 /)
         x_lower(1:n) = -1d0
         x_upper(1:n) = 3d0

         enforce_bounds = .true.
         adaptive_random_search = .true.
         centroid_weight_power = 1d0
         
         call test1_simplex( &
            n, x_first, x_lower, x_upper, simplex, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, fcn_PS, 'PS')
         
      end subroutine test_PS


      real(dp) function fcn_PS(n, x, lrpar, rpar, lipar, ipar, op_code, ierr)
         use const_def, only: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(:) ! (n)
         integer, intent(in) :: lrpar, lipar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: op_code
         integer, intent(out) :: ierr
         ierr = 0
         fcn_PS = PS(x)
      end function fcn_PS

      
      real(dp) function PS ( x ) ! Powell singular function
         real(dp), intent(in) :: x(:)
         integer :: i, n
         n = size(x,dim=1)
         PS = 0
         do i = 1, n/4
            PS = PS + &
               pow2(x(i) + 10d0*x(i+1)) + &
               5d0*pow2(x(i+2) - x(i+3)) + &
               pow4(x(i+1) - 2d0*x(i+2)) + &
               10d0*pow4(x(i) - x(i+3))
         end do
         num_calls = num_calls + 1
      end function PS
      
      
      subroutine test_TR
         ! testing with 4 dimensional Trigonometric function
         integer, parameter :: n = 4
         real(dp), dimension(n) :: x_first, x_lower, x_upper
         real(dp) :: simplex(n,n+1), centroid_weight_power
         logical :: enforce_bounds, adaptive_random_search
         
         x_first(1:n) = 0.01d0
         x_lower(1:n) = -0.1d0
         x_upper(1:n) = 0.1d0

         enforce_bounds = .true.
         adaptive_random_search = .true.
         centroid_weight_power = 1d0
         
         call test1_simplex( &
            n, x_first, x_lower, x_upper, simplex, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, fcn_TR, 'TR')
         
      end subroutine test_TR


      real(dp) function fcn_TR(n, x, lrpar, rpar, lipar, ipar, op_code, ierr)
         use const_def, only: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(:) ! (n)
         integer, intent(in) :: lrpar, lipar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: op_code
         integer, intent(out) :: ierr
         ierr = 0
         fcn_TR = TR(x)
      end function fcn_TR

      
      real(dp) function TR ( x ) ! Trigonometric function
         real(dp), intent(in) :: x(:)
         integer :: j, i, n
         real(dp) :: sum_cos
         n = size(x,dim=1)
         sum_cos = 0
         do j = 1, n
            sum_cos = sum_cos + cos(x(j))
         end do
         TR = 0
         do i = 1, n
            TR = TR + pow2(n - sum_cos + i*(1d0 - cos(x(i))) - sin(x(i)))
         end do
         num_calls = num_calls + 1
      end function TR



      end module test_simplex
