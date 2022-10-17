program test_auto_diff
   use math_lib
   use auto_diff
   use const_def
   
   call do_test_auto_diff_1var_order1()
   call do_test_auto_diff_2var_order1()
   call do_test_auto_diff_star_order1()

contains
   
   subroutine header(text)
      character(len = *), intent(in) :: text
      
      write(*, '(a)') ' ----------------------------------------------------------------'
      write(*, '(a)') ' '
      write(*, '(a)') ' ' // text
      write(*, '(a)') ' '
      write(*, '(a)') ' ----------------------------------------------------------------'
   
   end subroutine header
   
   subroutine should_print0(affix, a, z)
      character(len = *), intent(in) :: affix ! to insert ' approximately'
      real(dp), intent(in) :: a, z
      
      write(*, '(2(a),1(1pd26.16),a,99(1pd26.16))') &
         ' Should print', affix, a, '  :  ', z
      write(*, '(a)') ''
   end subroutine should_print0
   
   subroutine should_print1(affix, a, b, z)
      character(len = *), intent(in) :: affix ! to insert ' approximately'
      real(dp), intent(in) :: a, b
      type(auto_diff_real_1var_order1), intent(in) :: z
      
      write(*, '(2(a),2(1pd26.16),a,99(1pd26.16))') &
         ' Should print', affix, a, b, '  :  ', z
      write(*, '(a)') ''
   end subroutine should_print1
   
   
   subroutine should_print2(affix, a, b, c, z)
      character(len = *), intent(in) :: affix ! to insert ' approximately'
      real(dp), intent(in) :: a, b, c
      type(auto_diff_real_2var_order1) :: z
      write(*, '(2(a),3(1pd26.16),a,99(1pd26.16))') &
         ' Should print', affix, a, b, c, '  :  ', z
      write(*, '(a)') ''
   end subroutine should_print2
   
   
   subroutine do_test_auto_diff_star_order1()
      type(auto_diff_real_star_order1) :: x, y, z
      real(dp) :: a, b, c
      integer :: i, j, k
      
      call header('Testing assignment')
      x = 3d0
      x%d1Array(4) = 1d0
      do i = 1, 15
         if (i /= 4) then
            call should_print0('', 0d0, x%d1Array(i))
         else
            call should_print0('', 1d0, x%d1Array(i))
         end if
      end do
      
      call header('Testing unary operators')
      
      write(*, *) 'Test y = x**2'
      x = 3d0
      do i = 1, 15
         x%d1Array(i) = 1d0 * i
      end do
      y = pow2(x)
      do i = 1, 15
         call should_print0('', 6d0 * i, y%d1Array(i))
      end do
      
      write(*, *) 'Test x = x**2'
      x = pow2(x)
      do i = 1, 15
         call should_print0('', 6d0 * i, x%d1Array(i))
      end do
      
      call header('Testing binary operators')
      
      write(*, *) 'Test y = x + 1d0'
      x = 3d0
      do i = 1, 15
         x%d1Array(i) = 1d0 * i
      end do
      y = x + 1d0
      do i = 1, 15
         call should_print0('', 1d0 * i, y%d1Array(i))
      end do
      
      write(*, *) 'Test y = x + x'
      x = 3d0
      do i = 1, 15
         x%d1Array(i) = 1d0 * i
      end do
      y = x + x
      do i = 1, 15
         call should_print0('', 2d0 * i, y%d1Array(i))
      end do
      
      write(*, *) 'Test y = exp(x) * z'
      x = 3d0
      z = 2d0
      do i = 1, 15
         x%d1Array(i) = 1d0 * i
         z%d1Array(i) = 1d0 - i
      end do
      y = exp(x) * z
      do i = 1, 15
         call should_print0('', exp(3d0) * (1d0 * i) * 2d0 + (1d0 - i) * exp(3d0), y%d1Array(i))
      end do
   
   end subroutine do_test_auto_diff_star_order1
   
   subroutine do_test_auto_diff_1var_order1()
      type(auto_diff_real_1var_order1) :: x, y, z
      real(dp) :: a, b, c
      integer :: i, j, k
      
      call header('Testing assignment and comparison')
      
      write(*, *) 'Test real -> auto_diff_real_1var_order1 assignment'
      x = 3.1d0
      call should_print1('', 3.1d0, 0d0, x)
      
      write(*, *) 'Test auto_diff_real_1var_order1 == auto_diff_real_1var_order1'
      y = 3.1d0
      write(*, *) 'Should print T:', (x == y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 > auto_diff_real_1var_order1'
      y = 2.9d0
      write(*, *) 'Should print T:', (x > y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 < auto_diff_real_1var_order1'
      write(*, *) 'Should print F:', (x < y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 >= auto_diff_real_1var_order1'
      write(*, *) 'Should print T:', (x >= y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 <= auto_diff_real_1var_order1'
      write(*, *) 'Should print F:', (x <= y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 >= auto_diff_real_1var_order1'
      y = 3.1d0
      write(*, *) 'Should print T:', (x >= y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 <= auto_diff_real_1var_order1'
      write(*, *) 'Should print T:', (x <= y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 == real(dp)'
      write(*, *) 'Should print T:', (x == 3.1d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 > real(dp)'
      write(*, *) 'Should print T:', (x > 2.9d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 < real(dp)'
      write(*, *) 'Should print F:', (x < 2.9d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 >= real(dp)'
      write(*, *) 'Should print T:', (x >= 2.9d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 <= real(dp)'
      write(*, *) 'Should print F:', (x <= 2.9d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 >= real(dp)'
      write(*, *) 'Should print T:', (x >= 3.1d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 <= real(dp)'
      write(*, *) 'Should print T:', (x <= 3.1d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 == integer'
      x = 3
      write(*, *) 'Should print T:', (x == 3)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 > integer'
      write(*, *) 'Should print T:', (x > 2)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 < integer'
      write(*, *) 'Should print F:', (x < 2)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 >= integer'
      write(*, *) 'Should print T:', (x >= 3)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 <= integer'
      write(*, *) 'Should print F:', (x <= 2)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 >= integer'
      write(*, *) 'Should print T:', (x >= 2)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_1var_order1 <= integer'
      write(*, *) 'Should print T:', (x <= 3)
      write(*, *) ''
      
      call header('Testing binary operators')
      
      write(*, *) 'Test auto_diff_real_1var_order1+auto_diff_real_1var_order1'
      x = 1d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val1 = 1d0
      z = x + y
      call should_print1('', 3d0, 2d0, z)
      
      write(*, *) 'Test auto_diff_real_1var_order1+real(dp)'
      x = 1d0
      x%d1val1 = 1d0
      z = x + 2d0
      call should_print1('', 3d0, 1d0, z)
      
      write(*, *) 'Test real(dp)+auto_diff_real_1var_order1'
      x = 1d0
      x%d1val1 = 1d0
      z = 2d0 + x
      call should_print1('', 3d0, 1d0, z)
      
      write(*, *) 'Test auto_diff_real_1var_order1*auto_diff_real_1var_order1'
      x = 3d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val1 = 1d0
      z = x * y
      call should_print1('', 6d0, 5d0, z)
      
      write(*, *) 'Test auto_diff_real_1var_order1*real(dp)'
      x = 3d0
      x%d1val1 = 1d0
      z = x * 2d0
      call should_print1('', 6d0, 2d0, z)
      
      write(*, *) 'Test real(dp)*auto_diff_real_1var_order1'
      x = 3d0
      x%d1val1 = 1d0
      z = 2d0 * x
      call should_print1('', 6d0, 2d0, z)
      
      write(*, *) 'Test max(auto_diff_real_1var_order1, auto_diff_real_1var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val1 = 2d0
      z = max(x, y)
      call should_print1('', 2d0, 2d0, z)
      
      write(*, *) 'Test max(auto_diff_real_1var_order1, real(dp))'
      x = 1d0
      x%d1val1 = 1d0
      z = max(x, 2d0)
      call should_print1('', 2d0, 0d0, z)
      
      write(*, *) 'Test max(auto_diff_real_1var_order1, real(dp))'
      x = 1d0
      x%d1val1 = 1d0
      z = max(x, 0d0)
      call should_print1('', 1d0, 1d0, z)
      
      write(*, *) 'Test max(real(dp), auto_diff_real_1var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      z = max(2d0, x)
      call should_print1('', 2d0, 0d0, z)
      
      write(*, *) 'Test max(real(dp), auto_diff_real_1var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      z = max(0d0, x)
      call should_print1('', 1d0, 1d0, z)
      
      write(*, *) 'Test dim(auto_diff_real_1var_order1, auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 2d0
      y = 1d0
      y%d1val1 = 1d0
      z = dim(x, y)
      call should_print1('', 1d0, 1d0, z)
      
      write(*, *) 'Test dim(auto_diff_real_1var_order1, auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 2d0
      y = 1d0
      y%d1val1 = 1d0
      z = dim(y, x)
      call should_print1('', 0d0, 0d0, z)
      
      write(*, *) 'Test dim(auto_diff_real_1var_order1, real(dp))'
      x = 2d0
      x%d1val1 = 2d0
      z = dim(x, 1d0)
      call should_print1('', 1d0, 2d0, z)
      
      call header('Testing unary operators')
      
      write(*, *) 'Testing exp(auto_diff_real_1var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = exp(x)
      call should_print1('', exp(1d0), exp(1d0), x)
      
      write(*, *) 'Testing abs(auto_diff_real_1var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = abs(x)
      call should_print1('', 1d0, 1d0, x)
      
      write(*, *) 'Testing abs(auto_diff_real_1var_order1)'
      x = -1d0
      x%d1val1 = 1d0
      x = abs(x)
      call should_print1('', 1d0, -1d0, x)
      
      write(*, *) 'Testing exp(auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = exp(x)
      call should_print1('', exp(2d0), exp(2d0), x)
      
      write(*, *) 'Testing sin(auto_diff_real_1var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = sin(x)
      call should_print1('', sin(pi), cos(pi), x)
      
      write(*, *) 'Testing cos(auto_diff_real_1var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = cos(x)
      call should_print1('', cos(pi), -sin(pi), x)
      
      write(*, *) 'Testing tan(auto_diff_real_1var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = tan(x)
      call should_print1(' approximately', tan(pi), 1 / pow2(cos(pi)), x)
      
      write(*, *) 'Testing sinh(auto_diff_real_1var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = sinh(x)
      call should_print1('', sinh(pi), cosh(pi), x)
      
      write(*, *) 'Testing cosh(auto_diff_real_1var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = cosh(x)
      call should_print1('', cosh(pi), sinh(pi), x)
      
      write(*, *) 'Testing tanh(auto_diff_real_1var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = tanh(x)
      call should_print1(' approximately', tanh(pi), 1 / pow2(cosh(pi)), x)
      
      write(*, *) 'Testing asin(auto_diff_real_1var_order1)'
      x = 0.5d0
      x%d1val1 = 1d0
      x = asin(x)
      call should_print1('', asin(0.5d0), 2d0 / sqrt(3d0), x)
      
      write(*, *) 'Testing acos(auto_diff_real_1var_order1)'
      x = 0.5d0
      x%d1val1 = 1d0
      x = acos(x)
      call should_print1('', acos(0.5d0), -2d0 / sqrt(3d0), x)
      
      write(*, *) 'Testing atan(auto_diff_real_1var_order1)'
      x = 0.5d0
      x%d1val1 = 1d0
      x = atan(x)
      call should_print1(' approximately', atan(0.5d0), 4d0 / 5d0, x)
      
      write(*, *) 'Testing asinh(auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = asinh(x)
      call should_print1(' approximately', asinh(2d0), 1d0 / sqrt(5d0), x)
      
      write(*, *) 'Testing acosh(auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = acosh(x)
      call should_print1(' approximately', acosh(2d0), 1d0 / sqrt(3d0), x)
      
      write(*, *) 'Testing atanh(auto_diff_real_1var_order1)'
      x = 0.5d0
      x%d1val1 = 1d0
      x = atanh(x)
      call should_print1(' approximately', atanh(0.5d0), 4d0 / 3d0, x)
      
      write(*, *) 'Testing unary_minus(auto_diff_real_1var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = -x
      call should_print1('', -1d0, -1d0, x)
      
      write(*, *) 'Testing unary_minus(auto_diff_real_1var_order1)'
      x = -1d0
      x%d1val1 = 1d0
      x = -x
      call should_print1('', 1d0, -1d0, x)
      
      write(*, *) 'Testing log(auto_diff_real_1var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = log(x)
      call should_print1('', 0d0, 1d0, x)
      
      write(*, *) 'Testing safe_log(auto_diff_real_1var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = log(x)
      call should_print1('', 0d0, 1d0, x)
      
      write(*, *) 'Testing log10(auto_diff_real_1var_order1)'
      x = 1d1
      x%d1val1 = 1d0
      x = log10(x)
      call should_print1('', 1d0, 1d0 / (10 * ln10), x)
      
      write(*, *) 'Testing safe_log10(auto_diff_real_1var_order1)'
      x = 1d1
      x%d1val1 = 1d0
      x = safe_log10(x)
      call should_print1('', 1d0, 1d0 / (10 * ln10), x)
      
      write(*, *) 'Testing log(auto_diff_real_1var_order1)'
      x = exp(1d0)
      x%d1val1 = 1d0
      x = log(x)
      call should_print1('', 1d0, exp(-1d0), x)
      
      write(*, *) 'Testing safe_log(auto_diff_real_1var_order1)'
      x = exp(1d0)
      x%d1val1 = 1d0
      x = log(x)
      call should_print1(' approximately', 1d0, exp(-1d0), x)
      
      write(*, *) 'Testing log10(auto_diff_real_1var_order1)'
      x = 1d2
      x%d1val1 = 1d0
      x = log10(x)
      call should_print1(' approximately', 2d0, 1d0 / (100 * ln10), x)
      
      write(*, *) 'Testing safe_log10(auto_diff_real_1var_order1)'
      x = 1d2
      x%d1val1 = 1d0
      x = log10(x)
      call should_print1(' approximately', 2d0, 1d0 / (100 * ln10), x)
      
      write(*, *) 'Testing pow2(auto_diff_real_1var_order1)'
      x = 3d0
      x%d1val1 = 1d0
      x = pow2(x)
      call should_print1(' approximately', 9d0, 6d0, x)
      
      write(*, *) 'Testing pow3(auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow3(x)
      call should_print1('', 8d0, 3 * 4d0, x)
      
      write(*, *) 'Testing pow4(auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow4(x)
      call should_print1('', 16d0, 4 * 8d0, x)
      
      write(*, *) 'Testing pow5(auto_diff_real_1var_order1)'
      x = 3d0
      x%d1val1 = 1d0
      x = pow5(x)
      call should_print1('', 243d0, 405d0, x)
      
      write(*, *) 'Testing pow6(auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow6(x)
      call should_print1('', 64d0, 6 * 32d0, x)
      
      write(*, *) 'Testing pow7(auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow7(x)
      call should_print1('', 128d0, 7 * 64d0, x)
      
      write(*, *) 'Testing pow(auto_diff_real_1var_order1, real(dp))'
      x = 4d0
      x%d1val1 = 1d0
      x = pow(x, 0.5d0)
      call should_print1(' approximately', 2d0, 0.25d0, x)
      
      write(*, *) 'Testing pow(real(dp), auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow(exp(1d0), x)
      call should_print1(' approximately', exp(2d0), exp(2d0), x)
      
      write(*, *) 'Testing pow(auto_diff_real_1var_order1, integer)'
      x = 4d0
      x%d1val1 = 1d0
      x = pow(x, 2d0)
      call should_print1('', 16d0, 8d0, x)
      
      write(*, *) 'Testing pow(integer, auto_diff_real_1var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow(2d0, x)
      call should_print1('', 4d0, 4d0 * log(2d0), x)
   
   end subroutine do_test_auto_diff_1var_order1
   
   
   subroutine do_test_auto_diff_2var_order1()
      type(auto_diff_real_2var_order1) :: x, y, z
      real(dp) :: a, b, c
      integer :: i, j, k
      
      call header('Testing assignment and comparison')
      
      write(*, *) 'Test real -> auto_diff_real_2var_order1 assignment'
      x = 3.1d0
      call should_print2('', 3.1d0, 0d0, 0d0, x)
      
      write(*, *) 'Test auto_diff_real_2var_order1 == auto_diff_real_2var_order1'
      y = 3.1d0
      write(*, *) 'Should print T:', (x == y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 > auto_diff_real_2var_order1'
      y = 2.9d0
      write(*, *) 'Should print T:', (x > y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 < auto_diff_real_2var_order1'
      write(*, *) 'Should print F:', (x < y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 >= auto_diff_real_2var_order1'
      write(*, *) 'Should print T:', (x >= y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 <= auto_diff_real_2var_order1'
      write(*, *) 'Should print F:', (x <= y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 >= auto_diff_real_2var_order1'
      y = 3.1d0
      write(*, *) 'Should print T:', (x >= y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 <= auto_diff_real_2var_order1'
      write(*, *) 'Should print T:', (x <= y)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 == real(dp)'
      write(*, *) 'Should print T:', (x == 3.1d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 > real(dp)'
      write(*, *) 'Should print T:', (x > 2.9d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 < real(dp)'
      write(*, *) 'Should print F:', (x < 2.9d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 >= real(dp)'
      write(*, *) 'Should print T:', (x >= 2.9d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 <= real(dp)'
      write(*, *) 'Should print F:', (x <= 2.9d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 >= real(dp)'
      write(*, *) 'Should print T:', (x >= 3.1d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 <= real(dp)'
      write(*, *) 'Should print T:', (x <= 3.1d0)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 == integer'
      x = 3
      write(*, *) 'Should print T:', (x == 3)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 > integer'
      write(*, *) 'Should print T:', (x > 2)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 < integer'
      write(*, *) 'Should print F:', (x < 2)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 >= integer'
      write(*, *) 'Should print T:', (x >= 3)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 <= integer'
      write(*, *) 'Should print F:', (x <= 2)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 >= integer'
      write(*, *) 'Should print T:', (x >= 2)
      write(*, *) ''
      
      write(*, *) 'Test auto_diff_real_2var_order1 <= integer'
      write(*, *) 'Should print T:', (x <= 3)
      write(*, *) ''
      
      call header('Testing binary operators')
      
      write(*, *) 'Test auto_diff_real_2var_order1+auto_diff_real_2var_order1'
      x = 1d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val1 = 1d0
      z = x + y
      call should_print2('', 3d0, 2d0, 0d0, z)
      
      write(*, *) 'Test auto_diff_real_2var_order1+real(dp)'
      x = 1d0
      x%d1val1 = 1d0
      z = x + 2d0
      call should_print2('', 3d0, 1d0, 0d0, z)
      
      write(*, *) 'Test real(dp)+auto_diff_real_2var_order1'
      x = 1d0
      x%d1val1 = 1d0
      z = 2d0 + x
      call should_print2('', 3d0, 1d0, 0d0, z)
      
      write(*, *) 'Test auto_diff_real_2var_order1+integer'
      x = 1d0
      x%d1val1 = 1d0
      z = x + 2
      call should_print2('', 3d0, 1d0, 0d0, z)
      
      write(*, *) 'Test integer+auto_diff_real_2var_order1'
      x = 1d0
      x%d1val1 = 1d0
      z = 2 + x
      call should_print2('', 3d0, 1d0, 0d0, z)
      
      write(*, *) 'Test auto_diff_real_2var_order1*auto_diff_real_2var_order1'
      x = 3d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val2 = 1d0
      z = x * y
      call should_print2('', 6d0, 2d0, 3d0, z)
      
      write(*, *) 'Test auto_diff_real_2var_order1*real(dp)'
      x = 3d0
      x%d1val1 = 1d0
      z = x * 2d0
      call should_print2('', 6d0, 2d0, 0d0, z)
      
      write(*, *) 'Test auto_diff_real_2var_order1/auto_diff_real_2var_order1'
      x = 3d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val2 = 1d0
      z = x / y
      call should_print2('', 1.5d0, 0.5d0, -0.75d0, z)
      
      write(*, *) 'Test auto_diff_real_2var_order1/real(dp)'
      x = 3d0
      x%d1val1 = 1d0
      z = x / 2d0
      call should_print2('', 1.5d0, 0.5d0, 0d0, z)
      
      write(*, *) 'Test auto_diff_real_2var_order1-auto_diff_real_2var_order1'
      x = 1d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val2 = 1d0
      z = x - y
      call should_print2('', -1d0, 1d0, -1d0, z)
      
      write(*, *) 'Test auto_diff_real_2var_order1-auto_diff_real_2var_order1'
      x = 1d0
      x%d1val1 = 1d0
      y = 5d0
      y%d1val2 = 2d0
      z = x - y
      call should_print2('', -4d0, 1d0, -2d0, z)
      
      write(*, *) 'Test auto_diff_real_2var_order1-real(dp)'
      x = 2d0
      x%d1val1 = 1d0
      z = x - 1d0
      call should_print2('', 1d0, 1d0, 0d0, z)
      
      write(*, *) 'Test real(dp)*auto_diff_real_2var_order1'
      x = 3d0
      x%d1val1 = 1d0
      z = 2d0 * x
      call should_print2('', 6d0, 2d0, 0d0, z)
      
      write(*, *) 'Test max(auto_diff_real_2var_order1, auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val2 = 2d0
      z = max(x, y)
      call should_print2('', 2d0, 0d0, 2d0, z)
      
      write(*, *) 'Test max(auto_diff_real_2var_order1, real(dp))'
      x = 1d0
      x%d1val1 = 1d0
      z = max(x, 2d0)
      call should_print2('', 2d0, 0d0, 0d0, z)
      
      write(*, *) 'Test max(auto_diff_real_2var_order1, real(dp))'
      x = 1d0
      x%d1val1 = 1d0
      z = max(x, 0d0)
      call should_print2('', 1d0, 1d0, 0d0, z)
      
      write(*, *) 'Test max(real(dp), auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      z = max(2d0, x)
      call should_print2('', 2d0, 0d0, 0d0, z)
      
      write(*, *) 'Test max(real(dp), auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      z = max(0d0, x)
      call should_print2('', 1d0, 1d0, 0d0, z)
      
      write(*, *) 'Test min(auto_diff_real_2var_order1, auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val2 = 2d0
      z = min(x, y)
      call should_print2('', 1d0, 1d0, 0d0, z)
      
      write(*, *) 'Test min(auto_diff_real_2var_order1, real(dp))'
      x = 1d0
      x%d1val1 = 1d0
      z = min(x, 2d0)
      call should_print2('', 1d0, 1d0, 0d0, z)
      
      write(*, *) 'Test min(auto_diff_real_2var_order1, real(dp))'
      x = 1d0
      x%d1val1 = 1d0
      z = min(x, 0d0)
      call should_print2('', 0d0, 0d0, 0d0, z)
      
      write(*, *) 'Test min(real(dp), auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      z = min(2d0, x)
      call should_print2('', 1d0, 1d0, 0d0, z)
      
      write(*, *) 'Test min(real(dp), auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      z = min(0d0, x)
      call should_print2('', 0d0, 0d0, 0d0, z)
      write(*, *) 'Test dim(auto_diff_real_2var_order1, auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 2d0
      y = 1d0
      y%d1val2 = 1d0
      z = dim(x, y)
      call should_print2('', 1d0, 2d0, -1d0, z)
      
      write(*, *) 'Test dim(auto_diff_real_2var_order1, auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 2d0
      y = 1d0
      y%d1val2 = 1d0
      z = dim(y, x)
      call should_print2('', 0d0, 0d0, 0d0, z)
      
      write(*, *) 'Test dim(auto_diff_real_2var_order1, real(dp))'
      x = 2d0
      x%d1val1 = 2d0
      z = dim(x, 1d0)
      call should_print2('', 1d0, 2d0, 0d0, z)
      
      write(*, *) 'Testing pow(auto_diff_real_2var_order1, auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      y = 2d0
      y%d1val2 = 1d0
      z = pow(x, y)
      call should_print2(' approximately', 4d0, 2d0 * 2d0, 4 * log(2d0), z)
      
      call header('Testing unary operators')
      
      write(*, *) 'Testing exp(auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = exp(x)
      call should_print2('', exp(1d0), exp(1d0), 0d0, x)
      
      write(*, *) 'Testing exp(auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val2 = 1d0
      x = exp(x)
      call should_print2('', exp(1d0), 0d0, exp(1d0), x)
      
      write(*, *) 'Testing exp(auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = exp(x)
      call should_print2('', exp(2d0), exp(2d0), 0d0, x)
      
      write(*, *) 'Testing abs(auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = abs(x)
      call should_print2('', 1d0, 1d0, 0d0, x)
      
      write(*, *) 'Testing abs(auto_diff_real_2var_order1)'
      x = -1d0
      x%d1val1 = 1d0
      x = abs(x)
      call should_print2('', 1d0, -1d0, 0d0, x)
      
      write(*, *) 'Testing sin(auto_diff_real_2var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = sin(x)
      call should_print2('', sin(pi), cos(pi), 0d0, x)
      
      write(*, *) 'Testing cos(auto_diff_real_2var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = cos(x)
      call should_print2('', cos(pi), -sin(pi), 0d0, x)
      
      write(*, *) 'Testing tan(auto_diff_real_2var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = tan(x)
      call should_print2(' approximately', tan(pi), 1 / pow2(cos(pi)), 0d0, x)
      
      write(*, *) 'Testing sinh(auto_diff_real_2var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = sinh(x)
      call should_print2('', sinh(pi), cosh(pi), 0d0, x)
      
      write(*, *) 'Testing cosh(auto_diff_real_2var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = cosh(x)
      call should_print2('', cosh(pi), sinh(pi), 0d0, x)
      
      write(*, *) 'Testing tanh(auto_diff_real_2var_order1)'
      x = pi
      x%d1val1 = 1d0
      x = tanh(x)
      call should_print2(' approximately', tanh(pi), 1 / pow2(cosh(pi)), 0d0, x)
      
      write(*, *) 'Testing asin(auto_diff_real_2var_order1)'
      x = 0.5d0
      x%d1val1 = 1d0
      x = asin(x)
      call should_print2('', asin(0.5d0), 2d0 / sqrt(3d0), 0d0, x)
      
      write(*, *) 'Testing acos(auto_diff_real_2var_order1)'
      x = 0.5d0
      x%d1val1 = 1d0
      x = acos(x)
      call should_print2('', acos(0.5d0), -2d0 / sqrt(3d0), 0d0, x)
      
      write(*, *) 'Testing atan(auto_diff_real_2var_order1)'
      x = 0.5d0
      x%d1val1 = 1d0
      x = atan(x)
      call should_print2(' approximately', atan(0.5d0), 4d0 / 5d0, 0d0, x)
      
      write(*, *) 'Testing asinh(auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = asinh(x)
      call should_print2(' approximately', asinh(2d0), 1d0 / sqrt(5d0), 0d0, x)
      
      write(*, *) 'Testing acosh(auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = acosh(x)
      call should_print2(' approximately', acosh(2d0), 1d0 / sqrt(3d0), 0d0, x)
      
      write(*, *) 'Testing atanh(auto_diff_real_2var_order1)'
      x = 0.5d0
      x%d1val1 = 1d0
      x = atanh(x)
      call should_print2(' approximately', atanh(0.5d0), 4d0 / 3d0, 0d0, x)
      
      write(*, *) 'Testing unary_minus(auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = -x
      call should_print2('', -1d0, -1d0, 0d0, x)
      
      write(*, *) 'Testing unary_minus(auto_diff_real_2var_order1)'
      x = -1d0
      x%d1val1 = 1d0
      x = -x
      call should_print2('', 1d0, -1d0, 0d0, x)
      
      write(*, *) 'Testing log(auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = log(x)
      call should_print2('', 0d0, 1d0, 0d0, x)
      
      write(*, *) 'Testing safe_log(auto_diff_real_2var_order1)'
      x = 1d0
      x%d1val1 = 1d0
      x = log(x)
      call should_print2('', 0d0, 1d0, 0d0, x)
      
      write(*, *) 'Testing log10(auto_diff_real_2var_order1)'
      x = 1d1
      x%d1val1 = 1d0
      x = log10(x)
      call should_print2('', 1d0, 1d0 / (10 * ln10), 0d0, x)
      
      write(*, *) 'Testing safe_log10(auto_diff_real_2var_order1)'
      x = 1d1
      x%d1val1 = 1d0
      x = safe_log10(x)
      call should_print2('', 1d0, 1d0 / (10 * ln10), 0d0, x)
      
      write(*, *) 'Testing log(auto_diff_real_2var_order1)'
      x = exp(1d0)
      x%d1val1 = 1d0
      x = log(x)
      call should_print2('', 1d0, exp(-1d0), 0d0, x)
      
      write(*, *) 'Testing safe_log(auto_diff_real_2var_order1)'
      x = exp(1d0)
      x%d1val1 = 1d0
      x = log(x)
      call should_print2(' approximately', 1d0, exp(-1d0), 0d0, x)
      
      write(*, *) 'Testing log10(auto_diff_real_2var_order1)'
      x = 1d2
      x%d1val1 = 1d0
      x = log10(x)
      call should_print2(' approximately', 2d0, 1d0 / (100 * ln10), 0d0, x)
      
      write(*, *) 'Testing safe_log10(auto_diff_real_2var_order1)'
      x = 1d2
      x%d1val1 = 1d0
      x = log10(x)
      call should_print2(' approximately', 2d0, 1d0 / (100 * ln10), 0d0, x)
      
      write(*, *) 'Testing pow2(auto_diff_real_2var_order1)'
      x = 3d0
      x%d1val1 = 1d0
      x = pow2(x)
      call should_print2(' approximately', 9d0, 6d0, 0d0, x)
      
      write(*, *) 'Testing pow3(auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow3(x)
      call should_print2('', 8d0, 3 * 4d0, 0d0, x)
      
      write(*, *) 'Testing pow4(auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow4(x)
      call should_print2('', 16d0, 4 * 8d0, 0d0, x)
      
      write(*, *) 'Testing pow5(auto_diff_real_2var_order1)'
      x = 3d0
      x%d1val1 = 1d0
      x = pow5(x)
      call should_print2('', 243d0, 405d0, 0d0, x)
      
      write(*, *) 'Testing pow6(auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow6(x)
      call should_print2('', 64d0, 6 * 32d0, 0d0, x)
      
      write(*, *) 'Testing pow7(auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow7(x)
      call should_print2('', 128d0, 7 * 64d0, 0d0, x)
      
      write(*, *) 'Testing pow(auto_diff_real_2var_order1, real(dp))'
      x = 4d0
      x%d1val1 = 1d0
      x = pow(x, 0.5d0)
      call should_print2(' approximately', 2d0, 0.25d0, 0d0, x)
      
      write(*, *) 'Testing pow(real(dp), auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow(exp(1d0), x)
      call should_print2(' approximately', exp(2d0), exp(2d0), 0d0, x)
      
      write(*, *) 'Testing pow(auto_diff_real_2var_order1, integer)'
      x = 4d0
      x%d1val1 = 1d0
      x = pow(x, 2d0)
      call should_print2('', 16d0, 8d0, 0d0, x)
      
      write(*, *) 'Testing pow(integer, auto_diff_real_2var_order1)'
      x = 2d0
      x%d1val1 = 1d0
      x = pow(2d0, x)
      call should_print2('', 4d0, 4d0 * log(2d0), 0d0, x)
   
   end subroutine do_test_auto_diff_2var_order1

end program test_auto_diff
