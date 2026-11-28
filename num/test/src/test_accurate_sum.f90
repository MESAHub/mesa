module test_accurate_sum

   use const_def, only: dp, qp
   use auto_diff, only: auto_diff_real_star_order1, assignment(=)
   use accurate_sum, only: accurate_real, assignment(=), operator(+)
   use accurate_sum_auto_diff_star_order1

   implicit none

contains

   subroutine run_test_accurate_sum()

      call test_double_precision_addition()
      call test_quad_precision_addition()

   end subroutine run_test_accurate_sum

   subroutine test_double_precision_addition()
      real(dp), parameter :: expected = 1d0
      type(auto_diff_real_star_order1) :: term, res
      type(accurate_auto_diff_real_star_order1) :: sum

      ! The auto_diff version of compensated summation uses double precision.
      term = 1d40
      sum = term
      term = 1d0
      sum = sum + term
      term = -1d40
      sum = sum + term
      res = sum

      write (*, '(a40, 1pd26.16, a7, 1pd26.16)') &
         'compensated addition expected', expected, 'got', res%val

   end subroutine test_double_precision_addition

   subroutine test_quad_precision_addition()
      real(qp), parameter :: expected = 1.0_qp
      real(qp) :: res
      type(accurate_real) :: sum

      sum = 1.e40_qp
      sum = sum + 1.0_qp
      sum = sum + (-1.e40_qp)
      res = sum

      write (*, '(a40, 1pd26.16, a7, 1pd26.16)') &
         'compensated addition qp expected', expected, 'got', res

   end subroutine test_quad_precision_addition

end module test_accurate_sum
