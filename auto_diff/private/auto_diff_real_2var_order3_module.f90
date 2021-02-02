module auto_diff_real_2var_order3_module
      use const_def
      use utils_lib
      use support_functions
      use math_lib
   
      implicit none
      private
   public :: auto_diff_real_2var_order3, &
      assignment(=), &
      operator(.eq.), &
      operator(.ne.), &
      operator(.gt.), &
      operator(.lt.), &
      operator(.le.), &
      operator(.ge.), &
      make_unop, &
      make_binop, &
      operator(-), &
      exp, &
      expm1, &
      exp10, &
      powm1, &
      log, &
      log1p, &
      safe_log, &
      log10, &
      safe_log10, &
      log2, &
      sin, &
      cos, &
      tan, &
      sinpi, &
      cospi, &
      tanpi, &
      sinh, &
      cosh, &
      tanh, &
      asin, &
      acos, &
      atan, &
      asinpi, &
      acospi, &
      atanpi, &
      asinh, &
      acosh, &
      atanh, &
      sqrt, &
      pow2, &
      pow3, &
      pow4, &
      pow5, &
      pow6, &
      pow7, &
      pow8, &
      abs, &
      operator(+), &
      operator(*), &
      operator(/), &
      pow, &
      max, &
      min, &
      dim, &
      differentiate_1, &
      differentiate_2
   type :: auto_diff_real_2var_order3
      real(dp) :: val
      real(dp) :: d1val1
      real(dp) :: d1val2
      real(dp) :: d2val1
      real(dp) :: d1val1_d1val2
      real(dp) :: d2val2
      real(dp) :: d3val1
      real(dp) :: d2val1_d1val2
      real(dp) :: d1val1_d2val2
      real(dp) :: d3val2
   end type auto_diff_real_2var_order3
   
   interface assignment(=)
      module procedure assign_from_self
      module procedure assign_from_real_dp
      module procedure assign_from_int
   end interface assignment(=)
   
   interface operator(.eq.)
      module procedure equal_self
      module procedure equal_auto_diff_real_2var_order3_real_dp
      module procedure equal_real_dp_auto_diff_real_2var_order3
      module procedure equal_auto_diff_real_2var_order3_int
      module procedure equal_int_auto_diff_real_2var_order3
   end interface operator(.eq.)
   
   interface operator(.ne.)
      module procedure neq_self
      module procedure neq_auto_diff_real_2var_order3_real_dp
      module procedure neq_real_dp_auto_diff_real_2var_order3
      module procedure neq_auto_diff_real_2var_order3_int
      module procedure neq_int_auto_diff_real_2var_order3
   end interface operator(.ne.)
   
   interface operator(.gt.)
      module procedure greater_self
      module procedure greater_auto_diff_real_2var_order3_real_dp
      module procedure greater_real_dp_auto_diff_real_2var_order3
      module procedure greater_auto_diff_real_2var_order3_int
      module procedure greater_int_auto_diff_real_2var_order3
   end interface operator(.gt.)
   
   interface operator(.lt.)
      module procedure less_self
      module procedure less_auto_diff_real_2var_order3_real_dp
      module procedure less_real_dp_auto_diff_real_2var_order3
      module procedure less_auto_diff_real_2var_order3_int
      module procedure less_int_auto_diff_real_2var_order3
   end interface operator(.lt.)
   
   interface operator(.le.)
      module procedure leq_self
      module procedure leq_auto_diff_real_2var_order3_real_dp
      module procedure leq_real_dp_auto_diff_real_2var_order3
      module procedure leq_auto_diff_real_2var_order3_int
      module procedure leq_int_auto_diff_real_2var_order3
   end interface operator(.le.)
   
   interface operator(.ge.)
      module procedure geq_self
      module procedure geq_auto_diff_real_2var_order3_real_dp
      module procedure geq_real_dp_auto_diff_real_2var_order3
      module procedure geq_auto_diff_real_2var_order3_int
      module procedure geq_int_auto_diff_real_2var_order3
   end interface operator(.ge.)
   
   interface make_unop
      module procedure make_unary_operator
   end interface make_unop
   
   interface make_binop
      module procedure make_binary_operator
   end interface make_binop
   
   interface operator(-)
      module procedure unary_minus_self
   end interface operator(-)
   
   interface exp
      module procedure exp_self
   end interface exp
   
   interface expm1
      module procedure expm1_self
   end interface expm1
   
   interface exp10
      module procedure exp10_self
   end interface exp10
   
   interface powm1
      module procedure powm1_self
   end interface powm1
   
   interface log
      module procedure log_self
   end interface log
   
   interface log1p
      module procedure log1p_self
   end interface log1p
   
   interface safe_log
      module procedure safe_log_self
   end interface safe_log
   
   interface log10
      module procedure log10_self
   end interface log10
   
   interface safe_log10
      module procedure safe_log10_self
   end interface safe_log10
   
   interface log2
      module procedure log2_self
   end interface log2
   
   interface sin
      module procedure sin_self
   end interface sin
   
   interface cos
      module procedure cos_self
   end interface cos
   
   interface tan
      module procedure tan_self
   end interface tan
   
   interface sinpi
      module procedure sinpi_self
   end interface sinpi
   
   interface cospi
      module procedure cospi_self
   end interface cospi
   
   interface tanpi
      module procedure tanpi_self
   end interface tanpi
   
   interface sinh
      module procedure sinh_self
   end interface sinh
   
   interface cosh
      module procedure cosh_self
   end interface cosh
   
   interface tanh
      module procedure tanh_self
   end interface tanh
   
   interface asin
      module procedure asin_self
   end interface asin
   
   interface acos
      module procedure acos_self
   end interface acos
   
   interface atan
      module procedure atan_self
   end interface atan
   
   interface asinpi
      module procedure asinpi_self
   end interface asinpi
   
   interface acospi
      module procedure acospi_self
   end interface acospi
   
   interface atanpi
      module procedure atanpi_self
   end interface atanpi
   
   interface asinh
      module procedure asinh_self
   end interface asinh
   
   interface acosh
      module procedure acosh_self
   end interface acosh
   
   interface atanh
      module procedure atanh_self
   end interface atanh
   
   interface sqrt
      module procedure sqrt_self
   end interface sqrt
   
   interface pow2
      module procedure pow2_self
   end interface pow2
   
   interface pow3
      module procedure pow3_self
   end interface pow3
   
   interface pow4
      module procedure pow4_self
   end interface pow4
   
   interface pow5
      module procedure pow5_self
   end interface pow5
   
   interface pow6
      module procedure pow6_self
   end interface pow6
   
   interface pow7
      module procedure pow7_self
   end interface pow7
   
   interface pow8
      module procedure pow8_self
   end interface pow8
   
   interface abs
      module procedure abs_self
   end interface abs
   
   interface operator(+)
      module procedure add_self
      module procedure add_self_real
      module procedure add_real_self
      module procedure add_self_int
      module procedure add_int_self
   end interface operator(+)
   
   interface operator(-)
      module procedure sub_self
      module procedure sub_self_real
      module procedure sub_real_self
      module procedure sub_self_int
      module procedure sub_int_self
   end interface operator(-)
   
   interface operator(*)
      module procedure mul_self
      module procedure mul_self_real
      module procedure mul_real_self
      module procedure mul_self_int
      module procedure mul_int_self
   end interface operator(*)
   
   interface operator(/)
      module procedure div_self
      module procedure div_self_real
      module procedure div_real_self
      module procedure div_self_int
      module procedure div_int_self
   end interface operator(/)
   
   interface pow
      module procedure pow_self
      module procedure pow_self_real
      module procedure pow_real_self
      module procedure pow_self_int
      module procedure pow_int_self
   end interface pow
   
   interface max
      module procedure max_self
      module procedure max_self_real
      module procedure max_real_self
      module procedure max_self_int
      module procedure max_int_self
   end interface max
   
   interface min
      module procedure min_self
      module procedure min_self_real
      module procedure min_real_self
      module procedure min_self_int
      module procedure min_int_self
   end interface min
   
   interface dim
      module procedure dim_self
      module procedure dim_self_real
      module procedure dim_real_self
      module procedure dim_self_int
      module procedure dim_int_self
   end interface dim
   
   interface differentiate_1
      module procedure differentiate_auto_diff_real_2var_order3_1
   end interface differentiate_1
   
   interface differentiate_2
      module procedure differentiate_auto_diff_real_2var_order3_2
   end interface differentiate_2
   
   contains

   subroutine assign_from_self(this, other)
      type(auto_diff_real_2var_order3), intent(out) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      this%val = other%val
      this%d1val1 = other%d1val1
      this%d1val2 = other%d1val2
      this%d2val1 = other%d2val1
      this%d1val1_d1val2 = other%d1val1_d1val2
      this%d2val2 = other%d2val2
      this%d3val1 = other%d3val1
      this%d2val1_d1val2 = other%d2val1_d1val2
      this%d1val1_d2val2 = other%d1val1_d2val2
      this%d3val2 = other%d3val2
   end subroutine assign_from_self
   
   subroutine assign_from_real_dp(this, other)
      type(auto_diff_real_2var_order3), intent(out) :: this
      real(dp), intent(in) :: other
      this%val = other
      this%d1val1 = 0_dp
      this%d1val2 = 0_dp
      this%d2val1 = 0_dp
      this%d1val1_d1val2 = 0_dp
      this%d2val2 = 0_dp
      this%d3val1 = 0_dp
      this%d2val1_d1val2 = 0_dp
      this%d1val1_d2val2 = 0_dp
      this%d3val2 = 0_dp
   end subroutine assign_from_real_dp
   
   subroutine assign_from_int(this, other)
      type(auto_diff_real_2var_order3), intent(out) :: this
      integer, intent(in) :: other
      this%val = other
      this%d1val1 = 0_dp
      this%d1val2 = 0_dp
      this%d2val1 = 0_dp
      this%d1val1_d1val2 = 0_dp
      this%d2val2 = 0_dp
      this%d3val1 = 0_dp
      this%d2val1_d1val2 = 0_dp
      this%d1val1_d2val2 = 0_dp
      this%d3val2 = 0_dp
   end subroutine assign_from_int
   
   function equal_self(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this%val .eq. other%val)
   end function equal_self
   
   function equal_auto_diff_real_2var_order3_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .eq. other)
   end function equal_auto_diff_real_2var_order3_real_dp
   
   function equal_real_dp_auto_diff_real_2var_order3(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .eq. other%val)
   end function equal_real_dp_auto_diff_real_2var_order3
   
   function equal_auto_diff_real_2var_order3_int(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .eq. other)
   end function equal_auto_diff_real_2var_order3_int
   
   function equal_int_auto_diff_real_2var_order3(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .eq. other%val)
   end function equal_int_auto_diff_real_2var_order3
   
   function neq_self(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this%val .ne. other%val)
   end function neq_self
   
   function neq_auto_diff_real_2var_order3_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .ne. other)
   end function neq_auto_diff_real_2var_order3_real_dp
   
   function neq_real_dp_auto_diff_real_2var_order3(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .ne. other%val)
   end function neq_real_dp_auto_diff_real_2var_order3
   
   function neq_auto_diff_real_2var_order3_int(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .ne. other)
   end function neq_auto_diff_real_2var_order3_int
   
   function neq_int_auto_diff_real_2var_order3(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .ne. other%val)
   end function neq_int_auto_diff_real_2var_order3
   
   function greater_self(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this%val .gt. other%val)
   end function greater_self
   
   function greater_auto_diff_real_2var_order3_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .gt. other)
   end function greater_auto_diff_real_2var_order3_real_dp
   
   function greater_real_dp_auto_diff_real_2var_order3(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .gt. other%val)
   end function greater_real_dp_auto_diff_real_2var_order3
   
   function greater_auto_diff_real_2var_order3_int(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .gt. other)
   end function greater_auto_diff_real_2var_order3_int
   
   function greater_int_auto_diff_real_2var_order3(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .gt. other%val)
   end function greater_int_auto_diff_real_2var_order3
   
   function less_self(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this%val .lt. other%val)
   end function less_self
   
   function less_auto_diff_real_2var_order3_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .lt. other)
   end function less_auto_diff_real_2var_order3_real_dp
   
   function less_real_dp_auto_diff_real_2var_order3(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .lt. other%val)
   end function less_real_dp_auto_diff_real_2var_order3
   
   function less_auto_diff_real_2var_order3_int(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .lt. other)
   end function less_auto_diff_real_2var_order3_int
   
   function less_int_auto_diff_real_2var_order3(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .lt. other%val)
   end function less_int_auto_diff_real_2var_order3
   
   function leq_self(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this%val .le. other%val)
   end function leq_self
   
   function leq_auto_diff_real_2var_order3_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .le. other)
   end function leq_auto_diff_real_2var_order3_real_dp
   
   function leq_real_dp_auto_diff_real_2var_order3(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .le. other%val)
   end function leq_real_dp_auto_diff_real_2var_order3
   
   function leq_auto_diff_real_2var_order3_int(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .le. other)
   end function leq_auto_diff_real_2var_order3_int
   
   function leq_int_auto_diff_real_2var_order3(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .le. other%val)
   end function leq_int_auto_diff_real_2var_order3
   
   function geq_self(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this%val .ge. other%val)
   end function geq_self
   
   function geq_auto_diff_real_2var_order3_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .ge. other)
   end function geq_auto_diff_real_2var_order3_real_dp
   
   function geq_real_dp_auto_diff_real_2var_order3(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .ge. other%val)
   end function geq_real_dp_auto_diff_real_2var_order3
   
   function geq_auto_diff_real_2var_order3_int(this, other) result(z)
      type(auto_diff_real_2var_order3), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .ge. other)
   end function geq_auto_diff_real_2var_order3_int
   
   function geq_int_auto_diff_real_2var_order3(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3), intent(in) :: other
      logical :: z
      z = (this .ge. other%val)
   end function geq_int_auto_diff_real_2var_order3
   
   function make_unary_operator(x, z_val, z_d1x, z_d2x, z_d3x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      real(dp), intent(in) :: z_val
      real(dp), intent(in) :: z_d1x
      real(dp), intent(in) :: z_d2x
      real(dp), intent(in) :: z_d3x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%d1val1)
      q1 = x%d1val2*z_d2x
      q2 = pow2(x%d1val2)
      q3 = x%d1val1*z_d2x
      q4 = 2.0_dp*x%d1val1_d1val2
      unary%val = z_val
      unary%d1val1 = x%d1val1*z_d1x
      unary%d1val2 = x%d1val2*z_d1x
      unary%d2val1 = q0*z_d2x + x%d2val1*z_d1x
      unary%d1val1_d1val2 = q1*x%d1val1 + x%d1val1_d1val2*z_d1x
      unary%d2val2 = q2*z_d2x + x%d2val2*z_d1x
      unary%d3val1 = 3.0_dp*q3*x%d2val1 + x%d3val1*z_d1x + z_d3x*pow3(x%d1val1)
      unary%d2val1_d1val2 = q0*x%d1val2*z_d3x + q1*x%d2val1 + q3*q4 + x%d2val1_d1val2*z_d1x
      unary%d1val1_d2val2 = q1*q4 + q2*x%d1val1*z_d3x + q3*x%d2val2 + x%d1val1_d2val2*z_d1x
      unary%d3val2 = 3.0_dp*q1*x%d2val2 + x%d3val2*z_d1x + z_d3x*pow3(x%d1val2)
   end function make_unary_operator
   
   function make_binary_operator(x, y, z_val, z_d1x, z_d1y, z_d2x, z_d1x_d1y, z_d2y, z_d3x, z_d2x_d1y, z_d1x_d2y, z_d3y) result(binary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3), intent(in) :: y
      real(dp), intent(in) :: z_val
      real(dp), intent(in) :: z_d1x
      real(dp), intent(in) :: z_d1y
      real(dp), intent(in) :: z_d2x
      real(dp), intent(in) :: z_d1x_d1y
      real(dp), intent(in) :: z_d2y
      real(dp), intent(in) :: z_d3x
      real(dp), intent(in) :: z_d2x_d1y
      real(dp), intent(in) :: z_d1x_d2y
      real(dp), intent(in) :: z_d3y
      type(auto_diff_real_2var_order3) :: binary
      real(dp) :: q26
      real(dp) :: q25
      real(dp) :: q24
      real(dp) :: q23
      real(dp) :: q22
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%d1val1)
      q1 = pow2(y%d1val1)
      q2 = x%d1val1*z_d1x_d1y
      q3 = 2.0_dp*q2
      q4 = x%d1val2*z_d2x
      q5 = y%d1val2*z_d1x_d1y
      q6 = x%d1val2*z_d1x_d1y
      q7 = y%d1val2*z_d2y
      q8 = pow2(x%d1val2)
      q9 = pow2(y%d1val2)
      q10 = 2.0_dp*q5
      q11 = x%d1val1*z_d2x
      q12 = 3.0_dp*x%d2val1
      q13 = 3.0_dp*y%d2val1
      q14 = y%d1val1*z_d1x_d1y
      q15 = y%d1val1*z_d2y
      q16 = q1*z_d1x_d2y
      q17 = q0*z_d2x_d1y
      q18 = 2.0_dp*x%d1val1_d1val2
      q19 = 2.0_dp*y%d1val1_d1val2
      q20 = 2.0_dp*x%d1val1*y%d1val1
      q21 = x%d1val2*z_d2x_d1y
      q22 = y%d1val2*z_d1x_d2y
      q23 = q9*z_d1x_d2y
      q24 = q8*z_d2x_d1y
      q25 = 3.0_dp*x%d2val2
      q26 = 3.0_dp*y%d2val2
      binary%val = z_val
      binary%d1val1 = x%d1val1*z_d1x + y%d1val1*z_d1y
      binary%d1val2 = x%d1val2*z_d1x + y%d1val2*z_d1y
      binary%d2val1 = q0*z_d2x + q1*z_d2y + q3*y%d1val1 + x%d2val1*z_d1x + y%d2val1*z_d1y
      binary%d1val1_d1val2 = q4*x%d1val1 + q5*x%d1val1 + q6*y%d1val1 + q7*y%d1val1 + x%d1val1_d1val2*z_d1x + y%d1val1_d1val2*z_d1y
      binary%d2val2 = q10*x%d1val2 + q8*z_d2x + q9*z_d2y + x%d2val2*z_d1x + y%d2val2*z_d1y
      binary%d3val1 = 3.0_dp*q16*x%d1val1 + 3.0_dp*q17*y%d1val1 + q11*q12 + q12*q14 + q13*q15 + q13*q2 + x%d3val1*z_d1x + y%d3val1*z_d1y + z_d3x*pow3(x%d1val1) + z_d3y*pow3(y%d1val1)
      binary%d2val1_d1val2 = q0*x%d1val2*z_d3x + q1*y%d1val2*z_d3y + q11*q18 + q14*q18 + q15*q19 + q16*x%d1val2 + q17*y%d1val2 + q20*q21 + q20*q22 + q3*y%d1val1_d1val2 + q4*x%d2val1 + q5*x%d2val1 + q6*y%d2val1 + q7*y%d2val1 + x%d2val1_d1val2*z_d1x + y%d2val1_d1val2*z_d1y
      binary%d1val1_d2val2 = 2.0_dp*q21*x%d1val1*y%d1val2 + 2.0_dp*q22*x%d1val2*y%d1val1 + q10*x%d1val1_d1val2 + q11*x%d2val2 + q14*x%d2val2 + q15*y%d2val2 + q18*q4 + q19*q6 + q19*q7 + q2*y%d2val2 + q23*x%d1val1 + q24*y%d1val1 + q8*x%d1val1*z_d3x + q9*y%d1val1*z_d3y + x%d1val1_d2val2*z_d1x + y%d1val1_d2val2*z_d1y
      binary%d3val2 = 3.0_dp*q23*x%d1val2 + 3.0_dp*q24*y%d1val2 + q25*q4 + q25*q5 + q26*q6 + q26*q7 + x%d3val2*z_d1x + y%d3val2*z_d1y + z_d3x*pow3(x%d1val2) + z_d3y*pow3(y%d1val2)
   end function make_binary_operator
   
   function unary_minus_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      unary%val = -x%val
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d2val2 = -x%d2val2
      unary%d3val1 = -x%d3val1
      unary%d2val1_d1val2 = -x%d2val1_d1val2
      unary%d1val1_d2val2 = -x%d1val1_d2val2
      unary%d3val2 = -x%d3val2
   end function unary_minus_self
   
   function exp_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = exp(x%val)
      q1 = x%d2val1 + pow2(x%d1val1)
      q2 = x%d2val2 + pow2(x%d1val2)
      q3 = 2.0_dp*x%d1val1_d1val2
      unary%val = q0
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*q1
      unary%d1val1_d1val2 = q0*(x%d1val1*x%d1val2 + x%d1val1_d1val2)
      unary%d2val2 = q0*q2
      unary%d3val1 = q0*(3.0_dp*x%d1val1*x%d2val1 + x%d3val1 + pow3(x%d1val1))
      unary%d2val1_d1val2 = q0*(q1*x%d1val2 + q3*x%d1val1 + x%d2val1_d1val2)
      unary%d1val1_d2val2 = q0*(q2*x%d1val1 + q3*x%d1val2 + x%d1val1_d2val2)
      unary%d3val2 = q0*(3.0_dp*x%d1val2*x%d2val2 + x%d3val2 + pow3(x%d1val2))
   end function exp_self
   
   function expm1_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = exp(x%val)
      q1 = x%d2val1 + pow2(x%d1val1)
      q2 = x%d2val2 + pow2(x%d1val2)
      q3 = 2.0_dp*x%d1val1_d1val2
      unary%val = expm1(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*q1
      unary%d1val1_d1val2 = q0*(x%d1val1*x%d1val2 + x%d1val1_d1val2)
      unary%d2val2 = q0*q2
      unary%d3val1 = q0*(3.0_dp*x%d1val1*x%d2val1 + x%d3val1 + pow3(x%d1val1))
      unary%d2val1_d1val2 = q0*(q1*x%d1val2 + q3*x%d1val1 + x%d2val1_d1val2)
      unary%d1val1_d2val2 = q0*(q2*x%d1val1 + q3*x%d1val2 + x%d1val1_d2val2)
      unary%d3val2 = q0*(3.0_dp*x%d1val2*x%d2val2 + x%d3val2 + pow3(x%d1val2))
   end function expm1_self
   
   function exp10_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow(10.0_dp, x%val)
      q1 = ln10
      q2 = q0*q1
      q3 = q1*pow2(x%d1val1) + x%d2val1
      q4 = q1*x%d1val2
      q5 = q1*pow2(x%d1val2) + x%d2val2
      q6 = q1*x%d1val1
      q7 = pow2(q1)
      q8 = 2.0_dp*x%d1val1_d1val2
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q2*(q4*x%d1val1 + x%d1val1_d1val2)
      unary%d2val2 = q2*q5
      unary%d3val1 = q2*(3.0_dp*q6*x%d2val1 + q7*pow3(x%d1val1) + x%d3val1)
      unary%d2val1_d1val2 = q2*(q3*q4 + q6*q8 + x%d2val1_d1val2)
      unary%d1val1_d2val2 = q2*(q4*q8 + q5*q6 + x%d1val1_d2val2)
      unary%d3val2 = q2*(3.0_dp*q4*x%d2val2 + q7*pow3(x%d1val2) + x%d3val2)
   end function exp10_self
   
   function powm1_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = powm1(q0)
      q2 = powm1(pow3(x%val))
      q3 = x%d2val1*x%val
      q4 = pow2(x%d1val1)
      q5 = 2.0_dp*x%d1val1
      q6 = x%d1val1_d1val2*x%val
      q7 = x%d2val2*x%val
      q8 = -q7
      q9 = pow2(x%d1val2)
      q10 = powm1(pow4(x%val))
      q11 = 4.0_dp*q6
      q12 = 6.0_dp*x%d1val2
      unary%val = powm1(x%val)
      unary%d1val1 = -q1*x%d1val1
      unary%d1val2 = -q1*x%d1val2
      unary%d2val1 = q2*(2.0_dp*q4 - q3)
      unary%d1val1_d1val2 = q2*(q5*x%d1val2 - q6)
      unary%d2val2 = q2*(2.0_dp*q9 + q8)
      unary%d3val1 = q10*(-6.0_dp*pow3(x%d1val1) + 6.0_dp*q3*x%d1val1 - q0*x%d3val1)
      unary%d2val1_d1val2 = q10*(2.0_dp*q3*x%d1val2 - q0*x%d2val1_d1val2 + q11*x%d1val1 - q12*q4)
      unary%d1val1_d2val2 = q10*(-q0*x%d1val1_d2val2 + q11*x%d1val2 - q5*(3.0_dp*q9 + q8))
      unary%d3val2 = q10*(-6.0_dp*pow3(x%d1val2) - q0*x%d3val2 + q12*q7)
   end function powm1_self
   
   function log_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(x%val)
      q1 = pow2(x%val)
      q2 = powm1(q1)
      q3 = x%d2val1*x%val
      q4 = pow2(x%d1val1)
      q5 = x%d1val1_d1val2*x%val
      q6 = x%d2val2*x%val
      q7 = pow2(x%d1val2)
      q8 = powm1(pow3(x%val))
      q9 = 2.0_dp*x%d1val2
      unary%val = log(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*(q3 - q4)
      unary%d1val1_d1val2 = q2*(q5 - x%d1val1*x%d1val2)
      unary%d2val2 = q2*(q6 - q7)
      unary%d3val1 = q8*(-3.0_dp*q3*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q1*x%d3val1)
      unary%d2val1_d1val2 = q8*(-2.0_dp*q5*x%d1val1 + q1*x%d2val1_d1val2 - q3*x%d1val2 + q4*q9)
      unary%d1val1_d2val2 = q8*(q1*x%d1val1_d2val2 - q5*q9 + x%d1val1*(2.0_dp*q7 - q6))
      unary%d3val2 = q8*(-3.0_dp*q6*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q1*x%d3val2)
   end function log_self
   
   function log1p_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val + 1
      q1 = powm1(q0)
      q2 = pow2(q0)
      q3 = powm1(q2)
      q4 = pow2(x%d1val1)
      q5 = q0*x%d2val1
      q6 = q0*x%d1val1_d1val2
      q7 = pow2(x%d1val2)
      q8 = q0*x%d2val2
      q9 = powm1(pow3(q0))
      q10 = 2.0_dp*q6
      unary%val = log1p(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q3*(-q4 + q5)
      unary%d1val1_d1val2 = q3*(q6 - x%d1val1*x%d1val2)
      unary%d2val2 = q3*(-q7 + q8)
      unary%d3val1 = q9*(-3.0_dp*q5*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q2*x%d3val1)
      unary%d2val1_d1val2 = q9*(-q10*x%d1val1 + q2*x%d2val1_d1val2 + q4*x%d1val2 + x%d1val2*(q4 - q5))
      unary%d1val1_d2val2 = q9*(-q10*x%d1val2 + q2*x%d1val1_d2val2 + x%d1val1*(2.0_dp*q7 - q8))
      unary%d3val2 = q9*(-3.0_dp*q8*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q2*x%d3val2)
   end function log1p_self
   
   function safe_log_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(x%val)
      q1 = pow2(x%val)
      q2 = powm1(q1)
      q3 = x%d2val1*x%val
      q4 = pow2(x%d1val1)
      q5 = x%d1val1_d1val2*x%val
      q6 = x%d2val2*x%val
      q7 = pow2(x%d1val2)
      q8 = powm1(pow3(x%val))
      q9 = 2.0_dp*x%d1val2
      unary%val = safe_log(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*(q3 - q4)
      unary%d1val1_d1val2 = q2*(q5 - x%d1val1*x%d1val2)
      unary%d2val2 = q2*(q6 - q7)
      unary%d3val1 = q8*(-3.0_dp*q3*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q1*x%d3val1)
      unary%d2val1_d1val2 = q8*(-2.0_dp*q5*x%d1val1 + q1*x%d2val1_d1val2 - q3*x%d1val2 + q4*q9)
      unary%d1val1_d2val2 = q8*(q1*x%d1val1_d2val2 - q5*q9 + x%d1val1*(2.0_dp*q7 - q6))
      unary%d3val2 = q8*(-3.0_dp*q6*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q1*x%d3val2)
   end function safe_log_self
   
   function log10_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(ln10)
      q1 = q0*powm1(x%val)
      q2 = x%d2val1*x%val
      q3 = pow2(x%d1val1)
      q4 = pow2(x%val)
      q5 = q0*powm1(q4)
      q6 = x%d1val1_d1val2*x%val
      q7 = x%d2val2*x%val
      q8 = pow2(x%d1val2)
      q9 = q0*powm1(pow3(x%val))
      q10 = 2.0_dp*x%d1val2
      unary%val = q0*log(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q5*(q2 - q3)
      unary%d1val1_d1val2 = q5*(q6 - x%d1val1*x%d1val2)
      unary%d2val2 = q5*(q7 - q8)
      unary%d3val1 = q9*(-3.0_dp*q2*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q4*x%d3val1)
      unary%d2val1_d1val2 = q9*(-2.0_dp*q6*x%d1val1 + q10*q3 - q2*x%d1val2 + q4*x%d2val1_d1val2)
      unary%d1val1_d2val2 = q9*(-q10*q6 + q4*x%d1val1_d2val2 + x%d1val1*(2.0_dp*q8 - q7))
      unary%d3val2 = q9*(-3.0_dp*q7*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q4*x%d3val2)
   end function log10_self
   
   function safe_log10_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(ln10)
      q1 = q0*powm1(x%val)
      q2 = x%d2val1*x%val
      q3 = pow2(x%d1val1)
      q4 = pow2(x%val)
      q5 = q0*powm1(q4)
      q6 = x%d1val1_d1val2*x%val
      q7 = x%d2val2*x%val
      q8 = pow2(x%d1val2)
      q9 = q0*powm1(pow3(x%val))
      q10 = 2.0_dp*x%d1val2
      unary%val = q0*safe_log(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q5*(q2 - q3)
      unary%d1val1_d1val2 = q5*(q6 - x%d1val1*x%d1val2)
      unary%d2val2 = q5*(q7 - q8)
      unary%d3val1 = q9*(-3.0_dp*q2*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q4*x%d3val1)
      unary%d2val1_d1val2 = q9*(-2.0_dp*q6*x%d1val1 + q10*q3 - q2*x%d1val2 + q4*x%d2val1_d1val2)
      unary%d1val1_d2val2 = q9*(-q10*q6 + q4*x%d1val1_d2val2 + x%d1val1*(2.0_dp*q8 - q7))
      unary%d3val2 = q9*(-3.0_dp*q7*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q4*x%d3val2)
   end function safe_log10_self
   
   function log2_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(log(2.0_dp))
      q1 = q0*powm1(x%val)
      q2 = x%d2val1*x%val
      q3 = pow2(x%d1val1)
      q4 = pow2(x%val)
      q5 = q0*powm1(q4)
      q6 = x%d1val1_d1val2*x%val
      q7 = x%d2val2*x%val
      q8 = pow2(x%d1val2)
      q9 = q0*powm1(pow3(x%val))
      q10 = 2.0_dp*x%d1val2
      unary%val = q0*log(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q5*(q2 - q3)
      unary%d1val1_d1val2 = q5*(q6 - x%d1val1*x%d1val2)
      unary%d2val2 = q5*(q7 - q8)
      unary%d3val1 = q9*(-3.0_dp*q2*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q4*x%d3val1)
      unary%d2val1_d1val2 = q9*(-2.0_dp*q6*x%d1val1 + q10*q3 - q2*x%d1val2 + q4*x%d2val1_d1val2)
      unary%d1val1_d2val2 = q9*(-q10*q6 + q4*x%d1val1_d2val2 + x%d1val1*(2.0_dp*q8 - q7))
      unary%d3val2 = q9*(-3.0_dp*q7*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q4*x%d3val2)
   end function log2_self
   
   function sin_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = sin(x%val)
      q1 = cos(x%val)
      q2 = q1*x%d1val2
      q3 = pow2(x%d1val1)
      q4 = q0*x%d1val2
      q5 = pow2(x%d1val2)
      q6 = q0*x%d1val1
      q7 = 2.0_dp*x%d1val1_d1val2
      q8 = q0*x%d2val2
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q2
      unary%d2val1 = -q0*q3 + q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2 - q4*x%d1val1
      unary%d2val2 = -q0*q5 + q1*x%d2val2
      unary%d3val1 = -3.0_dp*q6*x%d2val1 + q1*x%d3val1 - q1*pow3(x%d1val1)
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2 - q2*q3 - q4*x%d2val1 - q6*q7
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2 - q4*q7 - x%d1val1*(q1*q5 + q8)
      unary%d3val2 = -3.0_dp*q8*x%d1val2 + q1*x%d3val2 - q1*pow3(x%d1val2)
   end function sin_self
   
   function cos_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = cos(x%val)
      q1 = sin(x%val)
      q2 = q1*x%d1val2
      q3 = pow2(x%d1val1)
      q4 = q0*x%d1val2
      q5 = pow2(x%d1val2)
      q6 = q0*x%d1val1
      q7 = 2.0_dp*x%d1val1_d1val2
      q8 = q0*x%d2val2
      unary%val = q0
      unary%d1val1 = -q1*x%d1val1
      unary%d1val2 = -q2
      unary%d2val1 = -q0*q3 - q1*x%d2val1
      unary%d1val1_d1val2 = -q1*x%d1val1_d1val2 - q4*x%d1val1
      unary%d2val2 = -q0*q5 - q1*x%d2val2
      unary%d3val1 = -3.0_dp*q6*x%d2val1 - q1*x%d3val1 + q1*pow3(x%d1val1)
      unary%d2val1_d1val2 = -q1*x%d2val1_d1val2 + q2*q3 - q4*x%d2val1 - q6*q7
      unary%d1val1_d2val2 = -q1*x%d1val1_d2val2 - q4*q7 + x%d1val1*(q1*q5 - q8)
      unary%d3val2 = -3.0_dp*q8*x%d1val2 - q1*x%d3val2 + q1*pow3(x%d1val2)
   end function cos_self
   
   function tan_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = tan(x%val)
      q1 = pow2(q0)
      q2 = q1 + 1
      q3 = powm1(pow2(cos(x%val)))
      q4 = pow2(x%d1val1)
      q5 = 2.0_dp*q0
      q6 = q4*q5 + x%d2val1
      q7 = q5*x%d1val2
      q8 = pow2(x%d1val2)
      q9 = q0*x%d1val1
      q10 = pow3(x%d1val1)
      q11 = 2.0_dp*q3
      q12 = 4.0_dp*q1
      q13 = 4.0_dp*x%d1val1_d1val2
      q14 = q0*x%d2val2
      q15 = pow3(x%d1val2)
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q3*q6
      unary%d1val1_d1val2 = q3*(q7*x%d1val1 + x%d1val1_d1val2)
      unary%d2val2 = q3*(q5*q8 + x%d2val2)
      unary%d3val1 = q3*(6.0_dp*q9*x%d2val1 + q10*q11 + q10*q12 + x%d3val1)
      unary%d2val1_d1val2 = q3*(q11*q4*x%d1val2 + q13*q9 + q6*q7 + x%d2val1_d1val2)
      unary%d1val1_d2val2 = q3*(2.0_dp*x%d1val1*(2.0_dp*q1*q8 + q14 + q3*q8) + q0*q13*x%d1val2 + x%d1val1_d2val2)
      unary%d3val2 = q3*(6.0_dp*q14*x%d1val2 + q11*q15 + q12*q15 + x%d3val2)
   end function tan_self
   
   function sinpi_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pi*x%val
      q1 = sin(q0)
      q2 = cos(q0)
      q3 = pi*q2
      q4 = pow2(x%d1val1)
      q5 = pi*q1
      q6 = q5*x%d1val2
      q7 = pow2(x%d1val2)
      q8 = q5*x%d1val1
      q9 = q2*pow2(pi)
      q10 = 2.0_dp*x%d1val1_d1val2
      q11 = q1*x%d2val2
      unary%val = q1
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d2val1 = pi*(q2*x%d2val1 - q4*q5)
      unary%d1val1_d1val2 = pi*(q2*x%d1val1_d1val2 - q6*x%d1val1)
      unary%d2val2 = pi*(q2*x%d2val2 - q5*q7)
      unary%d3val1 = pi*(-3.0_dp*q8*x%d2val1 + q2*x%d3val1 - q9*pow3(x%d1val1))
      unary%d2val1_d1val2 = pi*(-q10*q8 + q2*x%d2val1_d1val2 - q4*q9*x%d1val2 - q6*x%d2val1)
      unary%d1val1_d2val2 = -pi*(pi*x%d1val1*(q11 + q3*q7) + q10*q6 - q2*x%d1val1_d2val2)
      unary%d3val2 = pi*(-3.0_dp*pi*q11*x%d1val2 + q2*x%d3val2 - q9*pow3(x%d1val2))
   end function sinpi_self
   
   function cospi_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pi*x%val
      q1 = cos(q0)
      q2 = sin(q0)
      q3 = pi*q2
      q4 = pow2(x%d1val1)
      q5 = pi*q1
      q6 = q5*x%d1val2
      q7 = pow2(x%d1val2)
      q8 = q5*x%d1val1
      q9 = q2*pow2(pi)
      q10 = 2.0_dp*x%d1val1_d1val2
      q11 = q1*x%d2val2
      unary%val = q1
      unary%d1val1 = -q3*x%d1val1
      unary%d1val2 = -q3*x%d1val2
      unary%d2val1 = -pi*(q2*x%d2val1 + q4*q5)
      unary%d1val1_d1val2 = -pi*(q2*x%d1val1_d1val2 + q6*x%d1val1)
      unary%d2val2 = -pi*(q2*x%d2val2 + q5*q7)
      unary%d3val1 = pi*(-3.0_dp*q8*x%d2val1 - q2*x%d3val1 + q9*pow3(x%d1val1))
      unary%d2val1_d1val2 = pi*(-q10*q8 - q2*x%d2val1_d1val2 + q4*q9*x%d1val2 - q6*x%d2val1)
      unary%d1val1_d2val2 = -pi*(-pi*x%d1val1*(-q11 + q3*q7) + q10*q6 + q2*x%d1val1_d2val2)
      unary%d3val2 = pi*(-3.0_dp*pi*q11*x%d1val2 - q2*x%d3val2 + q9*pow3(x%d1val2))
   end function cospi_self
   
   function tanpi_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pi*x%val
      q1 = tan(q0)
      q2 = pow2(q1)
      q3 = q2 + 1
      q4 = pi*q3
      q5 = pi*q1
      q6 = 2.0_dp*pow2(x%d1val1)
      q7 = q5*q6 + x%d2val1
      q8 = 2.0_dp*q5
      q9 = q8*x%d1val2
      q10 = powm1(pow2(cos(q0)))
      q11 = pi*q10
      q12 = pow2(x%d1val2)
      q13 = 6.0_dp*q5
      q14 = pow2(pi)
      q15 = q14*pow3(x%d1val1)
      q16 = 4.0_dp*q2
      q17 = 2.0_dp*q3
      q18 = pi*x%d1val1
      q19 = 4.0_dp*q1
      q20 = q19*x%d1val1_d1val2
      q21 = q14*pow3(x%d1val2)
      unary%val = q1
      unary%d1val1 = q4*x%d1val1
      unary%d1val2 = q4*x%d1val2
      unary%d2val1 = q4*q7
      unary%d1val1_d1val2 = q11*(q9*x%d1val1 + x%d1val1_d1val2)
      unary%d2val2 = q4*(q12*q8 + x%d2val2)
      unary%d3val1 = q4*(q13*x%d1val1*x%d2val1 + q15*q16 + q15*q17 + x%d3val1)
      unary%d2val1_d1val2 = q11*(q10*q14*q6*x%d1val2 + q18*q20 + q7*q9 + x%d2val1_d1val2)
      unary%d1val1_d2val2 = q4*(0.5_dp*q18*(4.0_dp*q12*q4 + 8.0_dp*pi*q12*q2 + q19*x%d2val2) + pi*q20*x%d1val2 + x%d1val1_d2val2)
      unary%d3val2 = q4*(q13*x%d1val2*x%d2val2 + q16*q21 + q17*q21 + x%d3val2)
   end function tanpi_self
   
   function sinh_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = sinh(x%val)
      q1 = cosh(x%val)
      q2 = q1*x%d1val2
      q3 = pow2(x%d1val1)
      q4 = q0*x%d1val2
      q5 = pow2(x%d1val2)
      q6 = q0*x%d1val1
      q7 = 2.0_dp*x%d1val1_d1val2
      q8 = q0*x%d2val2
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q2
      unary%d2val1 = q0*q3 + q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2 + q4*x%d1val1
      unary%d2val2 = q0*q5 + q1*x%d2val2
      unary%d3val1 = 3.0_dp*q6*x%d2val1 + q1*x%d3val1 + q1*pow3(x%d1val1)
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2 + q2*q3 + q4*x%d2val1 + q6*q7
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2 + q4*q7 + x%d1val1*(q1*q5 + q8)
      unary%d3val2 = 3.0_dp*q8*x%d1val2 + q1*x%d3val2 + q1*pow3(x%d1val2)
   end function sinh_self
   
   function cosh_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = cosh(x%val)
      q1 = sinh(x%val)
      q2 = q1*x%d1val2
      q3 = pow2(x%d1val1)
      q4 = q0*x%d1val2
      q5 = pow2(x%d1val2)
      q6 = q0*x%d1val1
      q7 = 2.0_dp*x%d1val1_d1val2
      q8 = q0*x%d2val2
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q2
      unary%d2val1 = q0*q3 + q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2 + q4*x%d1val1
      unary%d2val2 = q0*q5 + q1*x%d2val2
      unary%d3val1 = 3.0_dp*q6*x%d2val1 + q1*x%d3val1 + q1*pow3(x%d1val1)
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2 + q2*q3 + q4*x%d2val1 + q6*q7
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2 + q4*q7 + x%d1val1*(q1*q5 + q8)
      unary%d3val2 = 3.0_dp*q8*x%d1val2 + q1*x%d3val2 + q1*pow3(x%d1val2)
   end function cosh_self
   
   function tanh_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = tanh(x%val)
      q1 = pow2(q0)
      q2 = q1 - 1
      q3 = powm1(pow2(cosh(x%val)))
      q4 = pow2(x%d1val1)
      q5 = 2.0_dp*q0
      q6 = q4*q5 - x%d2val1
      q7 = q5*x%d1val2
      q8 = pow2(x%d1val2)
      q9 = pow3(x%d1val1)
      q10 = q0*x%d1val1
      q11 = 6.0_dp*q1
      q12 = 4.0_dp*x%d1val1_d1val2
      q13 = q8*x%d1val1
      q14 = q0*x%d1val2
      q15 = pow3(x%d1val2)
      unary%val = q0
      unary%d1val1 = -q2*x%d1val1
      unary%d1val2 = -q2*x%d1val2
      unary%d2val1 = -q3*q6
      unary%d1val1_d1val2 = -q3*(q7*x%d1val1 - x%d1val1_d1val2)
      unary%d2val2 = -q3*(q5*q8 - x%d2val2)
      unary%d3val1 = q3*(-2.0_dp*q9 - 6.0_dp*q10*x%d2val1 + q11*q9 + x%d3val1)
      unary%d2val1_d1val2 = -q3*(2.0_dp*q3*q4*x%d1val2 + q10*q12 - q6*q7 - x%d2val1_d1val2)
      unary%d1val1_d2val2 = -q3*(2.0_dp*q13 - q11*q13 + q12*q14 + q5*x%d1val1*x%d2val2 - x%d1val1_d2val2)
      unary%d3val2 = q3*(-2.0_dp*q15 - 6.0_dp*q14*x%d2val2 + q11*q15 + x%d3val2)
   end function tanh_self
   
   function asin_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = 1 - q0
      q2 = powm1(sqrt(q1))
      q3 = powm1(pow3(sqrt(q1)))
      q4 = pow2(x%d1val1)
      q5 = q0 - 1
      q6 = q4*x%val - q5*x%d2val1
      q7 = x%d1val1*x%d1val2
      q8 = pow2(x%d1val2)
      q9 = powm1(pow5(sqrt(q1)))
      q10 = 3.0_dp*q0
      q11 = pow2(q5)
      q12 = q5*x%d1val1
      q13 = 2.0_dp*x%d1val1_d1val2*x%val
      q14 = q5*x%d1val2
      q15 = x%d2val2*x%val
      unary%val = asin(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q3*q6
      unary%d1val1_d1val2 = q3*(q1*x%d1val1_d1val2 + q7*x%val)
      unary%d2val2 = q3*(-q5*x%d2val2 + q8*x%val)
      unary%d3val1 = q9*(q10*pow3(x%d1val1) + q11*x%d3val1 - q12*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = (q1*q6*x%d1val2*x%val + q5*(-2.0_dp*q0*q4*x%d1val2 - q11*x%d2val1_d1val2 + q12*(q13 + q7)))*powm1(pow7(sqrt(q1)))
      unary%d1val1_d2val2 = q9*(q11*x%d1val1_d2val2 - q13*q14 - x%d1val1*(-q10*q8 + q5*(q15 + q8)))
      unary%d3val2 = q9*(q10*pow3(x%d1val2) + q11*x%d3val2 - q14*(3.0_dp*q15 + q8))
   end function asin_self
   
   function acos_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = 1 - q0
      q2 = powm1(sqrt(q1))
      q3 = powm1(pow3(sqrt(q1)))
      q4 = pow2(x%d1val1)
      q5 = q4*x%val
      q6 = q0 - 1
      q7 = q6*x%d2val1
      q8 = x%d1val1*x%d1val2
      q9 = q6*x%d1val1_d1val2
      q10 = pow2(x%d1val2)
      q11 = powm1(pow7(sqrt(q1)))
      q12 = pow3(q6)
      q13 = 3.0_dp*q0
      q14 = q1*q13
      q15 = pow2(q6)
      q16 = 2.0_dp*x%val
      q17 = x%d2val2*x%val
      unary%val = acos(x%val)
      unary%d1val1 = -q2*x%d1val1
      unary%d1val2 = -q2*x%d1val2
      unary%d2val1 = q3*(-q5 + q7)
      unary%d1val1_d1val2 = -q3*(q8*x%val - q9)
      unary%d2val2 = q3*(-q10*x%val + q6*x%d2val2)
      unary%d3val1 = q11*(q12*x%d3val1 - q14*pow3(x%d1val1) - q15*x%d1val1*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = q11*(-q1*x%d1val2*x%val*(q5 - q7) + q6*(2.0_dp*q0*q4*x%d1val2 + q15*x%d2val1_d1val2 - q6*x%d1val1*(q16*x%d1val1_d1val2 + q8)))
      unary%d1val1_d2val2 = (-q15*x%d1val1_d2val2 + q16*q9*x%d1val2 + x%d1val1*(-q10*q13 + q6*(q10 + q17)))*powm1(pow5(sqrt(q1)))
      unary%d3val2 = q11*(q12*x%d3val2 - q14*pow3(x%d1val2) - q15*x%d1val2*(3.0_dp*q17 + q10))
   end function acos_self
   
   function atan_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = q0 + 1
      q2 = powm1(q1)
      q3 = pow2(q1)
      q4 = powm1(q3)
      q5 = pow2(x%d1val1)
      q6 = 2.0_dp*x%val
      q7 = q5*q6
      q8 = q1*x%d2val1
      q9 = x%d1val1*x%d1val2
      q10 = q1*x%d1val1_d1val2
      q11 = pow2(x%d1val2)
      q12 = powm1(pow3(q1))
      q13 = 8.0_dp*q0
      q14 = 2.0_dp*x%d1val1
      q15 = q1*q14
      q16 = 4.0_dp*q0
      q17 = x%d2val2*x%val
      unary%val = atan(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q4*(-q7 + q8)
      unary%d1val1_d1val2 = q4*(q10 - q6*q9)
      unary%d2val2 = q4*(q1*x%d2val2 - q11*q6)
      unary%d3val1 = q12*(q13*pow3(x%d1val1) - q15*(3.0_dp*x%d2val1*x%val + q5) + q3*x%d3val1)
      unary%d2val1_d1val2 = q12*(-q15*(q6*x%d1val1_d1val2 + q9) + q16*q5*x%d1val2 + q3*x%d2val1_d1val2 + q6*x%d1val2*(q7 - q8))
      unary%d1val1_d2val2 = q12*(-4.0_dp*q10*x%d1val2*x%val - q14*(q1*(q11 + q17) - q11*q16) + q3*x%d1val1_d2val2)
      unary%d3val2 = q12*(-2.0_dp*q1*x%d1val2*(3.0_dp*q17 + q11) + q13*pow3(x%d1val2) + q3*x%d3val2)
   end function atan_self
   
   function asinpi_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(pi)
      q1 = pow2(x%val)
      q2 = 1 - q1
      q3 = q0*powm1(sqrt(q2))
      q4 = pow2(x%d1val1)
      q5 = q1 - 1
      q6 = q0*powm1(pow3(sqrt(q2)))
      q7 = x%d1val2*x%val
      q8 = pow2(x%d1val2)
      q9 = 3.0_dp*q1
      q10 = pow2(q5)
      q11 = q0*powm1(pow5(sqrt(q2)))
      q12 = pow4(x%val)
      q13 = 2.0_dp*q1
      q14 = q4*x%d1val2
      q15 = 2.0_dp*x%d1val1_d1val2
      q16 = q15*x%d1val1
      q17 = pow3(x%val)
      q18 = x%d2val2*x%val
      unary%val = q0*asin(x%val)
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d2val1 = q6*(q4*x%val - q5*x%d2val1)
      unary%d1val1_d1val2 = q6*(q2*x%d1val1_d1val2 + q7*x%d1val1)
      unary%d2val2 = q6*(-q5*x%d2val2 + q8*x%val)
      unary%d3val1 = q11*(q10*x%d3val1 - q5*x%d1val1*(3.0_dp*x%d2val1*x%val + q4) + q9*pow3(x%d1val1))
      unary%d2val1_d1val2 = q3*(q12*x%d2val1_d1val2 + q13*q14 - q13*x%d2val1_d1val2 + q14 - q16*q17 + q16*x%val - q17*x%d1val2*x%d2val1 + q7*x%d2val1 + x%d2val1_d1val2)*powm1(q12 - q13 + 1)
      unary%d1val1_d2val2 = q11*(q10*x%d1val1_d2val2 - q15*q5*q7 - x%d1val1*(q5*(q18 + q8) - q8*q9))
      unary%d3val2 = q11*(q10*x%d3val2 - q5*x%d1val2*(3.0_dp*q18 + q8) + q9*pow3(x%d1val2))
   end function asinpi_self
   
   function acospi_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(pi)
      q1 = pow2(x%val)
      q2 = 1 - q1
      q3 = q0*powm1(sqrt(q2))
      q4 = pow2(x%d1val1)
      q5 = q1 - 1
      q6 = q0*powm1(pow3(sqrt(q2)))
      q7 = x%d1val2*x%val
      q8 = q5*x%d1val1_d1val2
      q9 = pow2(x%d1val2)
      q10 = pow3(q5)
      q11 = 3.0_dp*q1
      q12 = q1*(3.0_dp - q11)
      q13 = pow2(q5)
      q14 = q0*powm1(pow7(sqrt(q2)))
      q15 = pow4(x%val)
      q16 = 2.0_dp*q1
      q17 = 2.0_dp*x%d1val1*x%d1val1_d1val2
      q18 = q4*x%d1val2
      q19 = pow3(x%val)
      q20 = x%d2val2*x%val
      unary%val = q0*acos(x%val)
      unary%d1val1 = -q3*x%d1val1
      unary%d1val2 = -q3*x%d1val2
      unary%d2val1 = q6*(-q4*x%val + q5*x%d2val1)
      unary%d1val1_d1val2 = -q6*(q7*x%d1val1 - q8)
      unary%d2val2 = q6*(q5*x%d2val2 - q9*x%val)
      unary%d3val1 = q14*(q10*x%d3val1 - q12*pow3(x%d1val1) - q13*x%d1val1*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = q3*(-q15*x%d2val1_d1val2 - q16*q18 + q16*x%d2val1_d1val2 + q17*q19 - q17*x%val - q18 + q19*x%d1val2*x%d2val1 - q7*x%d2val1 - x%d2val1_d1val2)*powm1(q15 - q16 + 1)
      unary%d1val1_d2val2 = q0*(2.0_dp*q7*q8 - q13*x%d1val1_d2val2 + x%d1val1*(-q11*q9 + q5*(q20 + q9)))*powm1(pow5(sqrt(q2)))
      unary%d3val2 = q14*(q10*x%d3val2 - q12*pow3(x%d1val2) - q13*x%d1val2*(3.0_dp*q20 + q9))
   end function acospi_self
   
   function atanpi_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q24
      real(dp) :: q23
      real(dp) :: q22
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = pi*q0
      q2 = powm1(pi + q1)
      q3 = pow2(x%d1val1)
      q4 = 2.0_dp*x%val
      q5 = pow4(x%val)
      q6 = pi*q5
      q7 = powm1(2.0_dp*q1 + pi + q6)
      q8 = q4*x%d1val2
      q9 = pow2(x%d1val2)
      q10 = powm1(3.0_dp*q1 + 3.0_dp*q6 + pi*pow6(x%val) + pi)
      q11 = pow3(x%d1val1)
      q12 = x%d1val1*x%val
      q13 = 6.0_dp*x%d2val1
      q14 = 2.0_dp*q0
      q15 = pow3(x%val)
      q16 = q15*x%d1val1
      q17 = 6.0_dp*q0
      q18 = 4.0_dp*x%d1val1_d1val2
      q19 = 2.0_dp*x%d1val2
      q20 = q18*x%d1val2
      q21 = 2.0_dp*x%d1val1
      q22 = q15*x%d2val2
      q23 = pow3(x%d1val2)
      q24 = 6.0_dp*x%d1val2
      unary%val = powm1(pi)*atan(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q7*(q0*x%d2val1 - q3*q4 + x%d2val1)
      unary%d1val1_d1val2 = q7*(q0*x%d1val1_d1val2 - q8*x%d1val1 + x%d1val1_d1val2)
      unary%d2val2 = q7*(q0*x%d2val2 - q4*q9 + x%d2val2)
      unary%d3val1 = q10*(-2.0_dp*q11 + q11*q17 - q12*q13 - q13*q16 + q14*x%d3val1 + q5*x%d3val1 + x%d3val1)
      unary%d2val1_d1val2 = q10*(-q12*q18 + q14*x%d2val1_d1val2 - q15*q19*x%d2val1 - q16*q18 + q17*q3*x%d1val2 - q19*q3 + q5*x%d2val1_d1val2 - q8*x%d2val1 + x%d2val1_d1val2)
      unary%d1val1_d2val2 = q10*(q14*x%d1val1_d2val2 - q15*q20 + q17*q9*x%d1val1 - q20*x%val - q21*q22 - q21*q9 - q4*x%d1val1*x%d2val2 + q5*x%d1val1_d2val2 + x%d1val1_d2val2)
      unary%d3val2 = q10*(-2.0_dp*q23 + q14*x%d3val2 + q17*q23 - q22*q24 - q24*x%d2val2*x%val + q5*x%d3val2 + x%d3val2)
   end function atanpi_self
   
   function asinh_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = q0 + 1
      q2 = powm1(sqrt(q1))
      q3 = powm1(pow3(sqrt(q1)))
      q4 = pow2(x%d1val1)
      q5 = q4*x%val
      q6 = q1*x%d2val1
      q7 = x%d1val1*x%d1val2
      q8 = q1*x%d1val1_d1val2
      q9 = pow2(x%d1val2)
      q10 = powm1(pow5(sqrt(q1)))
      q11 = 3.0_dp*q0
      q12 = pow2(q1)
      q13 = q1*x%d1val1
      q14 = x%d1val2*x%val
      q15 = x%d2val2*x%val
      unary%val = asinh(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q3*(-q5 + q6)
      unary%d1val1_d1val2 = q3*(-q7*x%val + q8)
      unary%d2val2 = q3*(q1*x%d2val2 - q9*x%val)
      unary%d3val1 = q10*(q11*pow3(x%d1val1) + q12*x%d3val1 - q13*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = q10*(2.0_dp*q0*q4*x%d1val2 + q12*x%d2val1_d1val2 - q13*(2.0_dp*x%d1val1_d1val2*x%val + q7) + q14*(q5 - q6))
      unary%d1val1_d2val2 = q10*(-2.0_dp*q14*q8 + q12*x%d1val1_d2val2 - x%d1val1*(q1*(q15 + q9) - q11*q9))
      unary%d3val2 = q10*(-q1*x%d1val2*(3.0_dp*q15 + q9) + q11*pow3(x%d1val2) + q12*x%d3val2)
   end function asinh_self
   
   function acosh_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = q0 - 1
      q2 = powm1(sqrt(q1))
      q3 = powm1(pow3(sqrt(q1)))
      q4 = pow2(x%d1val1)
      q5 = q4*x%val
      q6 = q1*x%d2val1
      q7 = x%d1val1*x%d1val2
      q8 = q1*x%d1val1_d1val2
      q9 = pow2(x%d1val2)
      q10 = powm1(pow5(sqrt(q1)))
      q11 = 3.0_dp*q0
      q12 = pow2(q1)
      q13 = q1*x%d1val1
      q14 = x%d1val2*x%val
      q15 = x%d2val2*x%val
      unary%val = acosh(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q3*(-q5 + q6)
      unary%d1val1_d1val2 = q3*(-q7*x%val + q8)
      unary%d2val2 = q3*(q1*x%d2val2 - q9*x%val)
      unary%d3val1 = q10*(q11*pow3(x%d1val1) + q12*x%d3val1 - q13*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = q10*(2.0_dp*q0*q4*x%d1val2 + q12*x%d2val1_d1val2 - q13*(2.0_dp*x%d1val1_d1val2*x%val + q7) + q14*(q5 - q6))
      unary%d1val1_d2val2 = q10*(-2.0_dp*q14*q8 + q12*x%d1val1_d2val2 - x%d1val1*(q1*(q15 + q9) - q11*q9))
      unary%d3val2 = q10*(-q1*x%d1val2*(3.0_dp*q15 + q9) + q11*pow3(x%d1val2) + q12*x%d3val2)
   end function acosh_self
   
   function atanh_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = q0 - 1
      q2 = powm1(q1)
      q3 = pow2(q1)
      q4 = powm1(q3)
      q5 = pow2(x%d1val1)
      q6 = 2.0_dp*x%val
      q7 = -q1*x%d2val1 + q5*q6
      q8 = x%d1val1*x%d1val2
      q9 = q1*x%d1val1_d1val2
      q10 = pow2(x%d1val2)
      q11 = powm1(pow4(q1))
      q12 = pow3(q1)
      q13 = 8.0_dp*q0*(1 - q0)
      q14 = 2.0_dp*x%d1val1
      q15 = powm1(q12)
      q16 = 4.0_dp*q0
      q17 = x%d2val2*x%val
      unary%val = atanh(x%val)
      unary%d1val1 = -q2*x%d1val1
      unary%d1val2 = -q2*x%d1val2
      unary%d2val1 = q4*q7
      unary%d1val1_d1val2 = q4*(q6*q8 - q9)
      unary%d2val2 = q4*(-q1*x%d2val2 + q10*q6)
      unary%d3val1 = q11*(-q12*x%d3val1 + q13*pow3(x%d1val1) + q14*q3*(3.0_dp*x%d2val1*x%val + q5))
      unary%d2val1_d1val2 = q15*(q1*q14*(q6*x%d1val1_d1val2 + q8) - q16*q5*x%d1val2 - q3*x%d2val1_d1val2 - q6*q7*x%d1val2)
      unary%d1val1_d2val2 = q15*(4.0_dp*q9*x%d1val2*x%val + q14*(q1*(q10 + q17) - q10*q16) - q3*x%d1val1_d2val2)
      unary%d3val2 = q11*(2.0_dp*q3*x%d1val2*(3.0_dp*q17 + q10) - q12*x%d3val2 + q13*pow3(x%d1val2))
   end function atanh_self
   
   function sqrt_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = sqrt(x%val)
      q1 = 0.5_dp*powm1(q0)
      q2 = 2.0_dp*x%val
      q3 = q2*x%d2val1
      q4 = pow2(x%d1val1)
      q5 = 0.25_dp*powm1(pow3(sqrt(x%val)))
      q6 = q2*x%d2val2
      q7 = pow2(x%d1val2)
      q8 = x%d1val1*x%val
      q9 = 4.0_dp*pow2(x%val)
      q10 = 0.125_dp*powm1(pow5(sqrt(x%val)))
      q11 = 4.0_dp*x%d1val1_d1val2
      q12 = x%d1val2*x%val
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q5*(q3 - q4)
      unary%d1val1_d1val2 = q5*(q2*x%d1val1_d1val2 - x%d1val1*x%d1val2)
      unary%d2val2 = q5*(q6 - q7)
      unary%d3val1 = q10*(-6.0_dp*q8*x%d2val1 + 3.0_dp*pow3(x%d1val1) + q9*x%d3val1)
      unary%d2val1_d1val2 = q10*(3.0_dp*q4*x%d1val2 - q11*q8 - q3*x%d1val2 + q9*x%d2val1_d1val2)
      unary%d1val1_d2val2 = q10*(-q11*q12 + q9*x%d1val1_d2val2 + x%d1val1*(3.0_dp*q7 - q6))
      unary%d3val2 = q10*(-6.0_dp*q12*x%d2val2 + 3.0_dp*pow3(x%d1val2) + q9*x%d3val2)
   end function sqrt_self
   
   function pow2_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = 2.0_dp*x%val
      q1 = 2.0_dp*x%d1val2
      q2 = 4.0_dp*x%d1val1_d1val2
      unary%val = pow2(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = 2.0_dp*pow2(x%d1val1) + q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q1*x%d1val1
      unary%d2val2 = 2.0_dp*pow2(x%d1val2) + q0*x%d2val2
      unary%d3val1 = 6.0_dp*x%d1val1*x%d2val1 + q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2 + q1*x%d2val1 + q2*x%d1val1
      unary%d1val1_d2val2 = 2.0_dp*x%d1val1*x%d2val2 + q0*x%d1val1_d2val2 + q2*x%d1val2
      unary%d3val2 = 6.0_dp*x%d1val2*x%d2val2 + q0*x%d3val2
   end function pow2_self
   
   function pow3_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = 3.0_dp*pow2(x%val)
      q1 = x%d2val1*x%val
      q2 = 2.0_dp*pow2(x%d1val1) + q1
      q3 = 3.0_dp*x%val
      q4 = x%d1val1_d1val2*x%val
      q5 = x%d2val2*x%val
      q6 = pow2(x%d1val2)
      unary%val = pow3(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(2.0_dp*x%d1val1*x%d1val2 + q4)
      unary%d2val2 = q3*(2.0_dp*q6 + q5)
      unary%d3val1 = 18.0_dp*q1*x%d1val1 + 6.0_dp*pow3(x%d1val1) + q0*x%d3val1
      unary%d2val1_d1val2 = 3.0_dp*q2*x%d1val2 + q3*(4.0_dp*x%d1val1*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val)
      unary%d1val1_d2val2 = 12.0_dp*q4*x%d1val2 + 6.0_dp*x%d1val1*(q5 + q6) + q0*x%d1val1_d2val2
      unary%d3val2 = 18.0_dp*q5*x%d1val2 + 6.0_dp*pow3(x%d1val2) + q0*x%d3val2
   end function pow3_self
   
   function pow4_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = 4.0_dp*pow3(x%val)
      q1 = x%d2val1*x%val
      q2 = 3.0_dp*pow2(x%d1val1) + q1
      q3 = pow2(x%val)
      q4 = 4.0_dp*q3
      q5 = x%d1val1_d1val2*x%val
      q6 = 3.0_dp*x%d1val1
      q7 = x%d2val2*x%val
      q8 = pow2(x%d1val2)
      q9 = 4.0_dp*x%val
      unary%val = pow4(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*q4
      unary%d1val1_d1val2 = q4*(q5 + q6*x%d1val2)
      unary%d2val2 = q4*(3.0_dp*q8 + q7)
      unary%d3val1 = q9*(6.0_dp*pow3(x%d1val1) + 9.0_dp*q1*x%d1val1 + q3*x%d3val1)
      unary%d2val1_d1val2 = q9*(2.0_dp*q2*x%d1val2 + x%val*(6.0_dp*x%d1val1*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d1val1_d2val2 = q9*(6.0_dp*q5*x%d1val2 + q3*x%d1val1_d2val2 + q6*(2.0_dp*q8 + q7))
      unary%d3val2 = q9*(6.0_dp*pow3(x%d1val2) + 9.0_dp*q7*x%d1val2 + q3*x%d3val2)
   end function pow4_self
   
   function pow5_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = 5.0_dp*pow4(x%val)
      q1 = x%d2val1*x%val
      q2 = 4.0_dp*pow2(x%d1val1) + q1
      q3 = 5.0_dp*pow3(x%val)
      q4 = x%d1val1_d1val2*x%val
      q5 = 4.0_dp*x%d1val1
      q6 = x%d2val2*x%val
      q7 = pow2(x%d1val2)
      q8 = pow2(x%val)
      q9 = 5.0_dp*q8
      unary%val = pow5(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(q4 + q5*x%d1val2)
      unary%d2val2 = q3*(4.0_dp*q7 + q6)
      unary%d3val1 = q9*(12.0_dp*q1*x%d1val1 + 12.0_dp*pow3(x%d1val1) + q8*x%d3val1)
      unary%d2val1_d1val2 = q9*(3.0_dp*q2*x%d1val2 + x%val*(8.0_dp*x%d1val1*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d1val1_d2val2 = q9*(8.0_dp*q4*x%d1val2 + q5*(3.0_dp*q7 + q6) + q8*x%d1val1_d2val2)
      unary%d3val2 = q9*(12.0_dp*q6*x%d1val2 + 12.0_dp*pow3(x%d1val2) + q8*x%d3val2)
   end function pow5_self
   
   function pow6_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = 6.0_dp*pow5(x%val)
      q1 = x%d2val1*x%val
      q2 = 5.0_dp*pow2(x%d1val1) + q1
      q3 = 6.0_dp*pow4(x%val)
      q4 = x%d1val1_d1val2*x%val
      q5 = 5.0_dp*x%d1val1
      q6 = x%d2val2*x%val
      q7 = pow2(x%d1val2)
      q8 = pow2(x%val)
      q9 = 6.0_dp*pow3(x%val)
      unary%val = pow6(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(q4 + q5*x%d1val2)
      unary%d2val2 = q3*(5.0_dp*q7 + q6)
      unary%d3val1 = q9*(15.0_dp*q1*x%d1val1 + 20.0_dp*pow3(x%d1val1) + q8*x%d3val1)
      unary%d2val1_d1val2 = q9*(4.0_dp*q2*x%d1val2 + x%val*(10.0_dp*x%d1val1*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d1val1_d2val2 = q9*(10.0_dp*q4*x%d1val2 + q5*(4.0_dp*q7 + q6) + q8*x%d1val1_d2val2)
      unary%d3val2 = q9*(15.0_dp*q6*x%d1val2 + 20.0_dp*pow3(x%d1val2) + q8*x%d3val2)
   end function pow6_self
   
   function pow7_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = 7.0_dp*pow6(x%val)
      q1 = x%d2val1*x%val
      q2 = 6.0_dp*pow2(x%d1val1) + q1
      q3 = 7.0_dp*pow5(x%val)
      q4 = x%d1val1_d1val2*x%val
      q5 = 6.0_dp*x%d1val1
      q6 = x%d2val2*x%val
      q7 = pow2(x%d1val2)
      q8 = pow2(x%val)
      q9 = 7.0_dp*pow4(x%val)
      unary%val = pow7(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(q4 + q5*x%d1val2)
      unary%d2val2 = q3*(6.0_dp*q7 + q6)
      unary%d3val1 = q9*(18.0_dp*q1*x%d1val1 + 30.0_dp*pow3(x%d1val1) + q8*x%d3val1)
      unary%d2val1_d1val2 = q9*(5.0_dp*q2*x%d1val2 + x%val*(12.0_dp*x%d1val1*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d1val1_d2val2 = q9*(12.0_dp*q4*x%d1val2 + q5*(5.0_dp*q7 + q6) + q8*x%d1val1_d2val2)
      unary%d3val2 = q9*(18.0_dp*q6*x%d1val2 + 30.0_dp*pow3(x%d1val2) + q8*x%d3val2)
   end function pow7_self
   
   function pow8_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = 8.0_dp*pow7(x%val)
      q1 = x%d2val1*x%val
      q2 = 7.0_dp*pow2(x%d1val1) + q1
      q3 = 8.0_dp*pow6(x%val)
      q4 = x%d1val1_d1val2*x%val
      q5 = 7.0_dp*x%d1val1
      q6 = x%d2val2*x%val
      q7 = pow2(x%d1val2)
      q8 = pow2(x%val)
      q9 = 8.0_dp*pow5(x%val)
      unary%val = pow8(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(q4 + q5*x%d1val2)
      unary%d2val2 = q3*(7.0_dp*q7 + q6)
      unary%d3val1 = q9*(21.0_dp*q1*x%d1val1 + 42.0_dp*pow3(x%d1val1) + q8*x%d3val1)
      unary%d2val1_d1val2 = q9*(6.0_dp*q2*x%d1val2 + x%val*(14.0_dp*x%d1val1*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d1val1_d2val2 = q9*(14.0_dp*q4*x%d1val2 + q5*(6.0_dp*q7 + q6) + q8*x%d1val1_d2val2)
      unary%d3val2 = q9*(21.0_dp*q6*x%d1val2 + 42.0_dp*pow3(x%d1val2) + q8*x%d3val2)
   end function pow8_self
   
   function abs_self(x) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q0
      q0 = sgn(x%val)
      unary%val = Abs(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function abs_self
   
   function add_self(x, y) result(binary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3), intent(in) :: y
      type(auto_diff_real_2var_order3) :: binary
      binary%val = x%val + y%val
      binary%d1val1 = x%d1val1 + y%d1val1
      binary%d1val2 = x%d1val2 + y%d1val2
      binary%d2val1 = x%d2val1 + y%d2val1
      binary%d1val1_d1val2 = x%d1val1_d1val2 + y%d1val1_d1val2
      binary%d2val2 = x%d2val2 + y%d2val2
      binary%d3val1 = x%d3val1 + y%d3val1
      binary%d2val1_d1val2 = x%d2val1_d1val2 + y%d2val1_d1val2
      binary%d1val1_d2val2 = x%d1val1_d2val2 + y%d1val1_d2val2
      binary%d3val2 = x%d3val2 + y%d3val2
   end function add_self
   
   function add_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      unary%val = x%val + y
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d3val2 = x%d3val2
   end function add_self_real
   
   function add_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      unary%val = x%val + z
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d3val2 = x%d3val2
   end function add_real_self
   
   function add_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val + y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d3val2 = x%d3val2
   end function add_self_int
   
   function add_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = x%val + y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d3val2 = x%d3val2
   end function add_int_self
   
   function sub_self(x, y) result(binary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3), intent(in) :: y
      type(auto_diff_real_2var_order3) :: binary
      binary%val = x%val - y%val
      binary%d1val1 = x%d1val1 - y%d1val1
      binary%d1val2 = x%d1val2 - y%d1val2
      binary%d2val1 = x%d2val1 - y%d2val1
      binary%d1val1_d1val2 = x%d1val1_d1val2 - y%d1val1_d1val2
      binary%d2val2 = x%d2val2 - y%d2val2
      binary%d3val1 = x%d3val1 - y%d3val1
      binary%d2val1_d1val2 = x%d2val1_d1val2 - y%d2val1_d1val2
      binary%d1val1_d2val2 = x%d1val1_d2val2 - y%d1val1_d2val2
      binary%d3val2 = x%d3val2 - y%d3val2
   end function sub_self
   
   function sub_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      unary%val = x%val - y
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d3val2 = x%d3val2
   end function sub_self_real
   
   function sub_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      unary%val = -x%val + z
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d2val2 = -x%d2val2
      unary%d3val1 = -x%d3val1
      unary%d2val1_d1val2 = -x%d2val1_d1val2
      unary%d1val1_d2val2 = -x%d1val1_d2val2
      unary%d3val2 = -x%d3val2
   end function sub_real_self
   
   function sub_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val - y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d3val2 = x%d3val2
   end function sub_self_int
   
   function sub_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = -x%val + y_dp
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d2val2 = -x%d2val2
      unary%d3val1 = -x%d3val1
      unary%d2val1_d1val2 = -x%d2val1_d1val2
      unary%d1val1_d2val2 = -x%d1val1_d2val2
      unary%d3val2 = -x%d3val2
   end function sub_int_self
   
   function mul_self(x, y) result(binary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3), intent(in) :: y
      type(auto_diff_real_2var_order3) :: binary
      real(dp) :: q1
      real(dp) :: q0
      q0 = 2.0_dp*x%d1val1
      q1 = 2.0_dp*y%d1val2
      binary%val = x%val*y%val
      binary%d1val1 = x%d1val1*y%val + x%val*y%d1val1
      binary%d1val2 = x%d1val2*y%val + x%val*y%d1val2
      binary%d2val1 = q0*y%d1val1 + x%d2val1*y%val + x%val*y%d2val1
      binary%d1val1_d1val2 = x%d1val1*y%d1val2 + x%d1val1_d1val2*y%val + x%d1val2*y%d1val1 + x%val*y%d1val1_d1val2
      binary%d2val2 = q1*x%d1val2 + x%d2val2*y%val + x%val*y%d2val2
      binary%d3val1 = 3.0_dp*x%d1val1*y%d2val1 + 3.0_dp*x%d2val1*y%d1val1 + x%d3val1*y%val + x%val*y%d3val1
      binary%d2val1_d1val2 = 2.0_dp*x%d1val1_d1val2*y%d1val1 + q0*y%d1val1_d1val2 + x%d1val2*y%d2val1 + x%d2val1*y%d1val2 + x%d2val1_d1val2*y%val + x%val*y%d2val1_d1val2
      binary%d1val1_d2val2 = 2.0_dp*x%d1val2*y%d1val1_d1val2 + q1*x%d1val1_d1val2 + x%d1val1*y%d2val2 + x%d1val1_d2val2*y%val + x%d2val2*y%d1val1 + x%val*y%d1val1_d2val2
      binary%d3val2 = 3.0_dp*x%d1val2*y%d2val2 + 3.0_dp*x%d2val2*y%d1val2 + x%d3val2*y%val + x%val*y%d3val2
   end function mul_self
   
   function mul_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      unary%val = x%val*y
      unary%d1val1 = x%d1val1*y
      unary%d1val2 = x%d1val2*y
      unary%d2val1 = x%d2val1*y
      unary%d1val1_d1val2 = x%d1val1_d1val2*y
      unary%d2val2 = x%d2val2*y
      unary%d3val1 = x%d3val1*y
      unary%d2val1_d1val2 = x%d2val1_d1val2*y
      unary%d1val1_d2val2 = x%d1val1_d2val2*y
      unary%d3val2 = x%d3val2*y
   end function mul_self_real
   
   function mul_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      unary%val = x%val*z
      unary%d1val1 = x%d1val1*z
      unary%d1val2 = x%d1val2*z
      unary%d2val1 = x%d2val1*z
      unary%d1val1_d1val2 = x%d1val1_d1val2*z
      unary%d2val2 = x%d2val2*z
      unary%d3val1 = x%d3val1*z
      unary%d2val1_d1val2 = x%d2val1_d1val2*z
      unary%d1val1_d2val2 = x%d1val1_d2val2*z
      unary%d3val2 = x%d3val2*z
   end function mul_real_self
   
   function mul_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val*y_dp
      unary%d1val1 = x%d1val1*y_dp
      unary%d1val2 = x%d1val2*y_dp
      unary%d2val1 = x%d2val1*y_dp
      unary%d1val1_d1val2 = x%d1val1_d1val2*y_dp
      unary%d2val2 = x%d2val2*y_dp
      unary%d3val1 = x%d3val1*y_dp
      unary%d2val1_d1val2 = x%d2val1_d1val2*y_dp
      unary%d1val1_d2val2 = x%d1val1_d2val2*y_dp
      unary%d3val2 = x%d3val2*y_dp
   end function mul_self_int
   
   function mul_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = x%val*y_dp
      unary%d1val1 = x%d1val1*y_dp
      unary%d1val2 = x%d1val2*y_dp
      unary%d2val1 = x%d2val1*y_dp
      unary%d1val1_d1val2 = x%d1val1_d1val2*y_dp
      unary%d2val2 = x%d2val2*y_dp
      unary%d3val1 = x%d3val1*y_dp
      unary%d2val1_d1val2 = x%d2val1_d1val2*y_dp
      unary%d1val1_d2val2 = x%d1val1_d2val2*y_dp
      unary%d3val2 = x%d3val2*y_dp
   end function mul_int_self
   
   function div_self(x, y) result(binary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3), intent(in) :: y
      type(auto_diff_real_2var_order3) :: binary
      real(dp) :: q22
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(y%val)
      q1 = powm1(q0)
      q2 = x%d1val1*y%val
      q3 = x%val*y%d1val1
      q4 = x%d1val2*y%val
      q5 = x%val*y%d1val2
      q6 = pow3(y%val)
      q7 = powm1(q6)
      q8 = q0*x%d2val1
      q9 = 2.0_dp*y%d1val1
      q10 = y%d2val1*y%val
      q11 = pow2(y%d1val1)
      q12 = 2.0_dp*q11
      q13 = -q10 + q12
      q14 = q0*x%d1val1_d1val2
      q15 = 2.0_dp*y%d1val2
      q16 = x%d1val2*y%d1val1
      q17 = q0*x%d2val2
      q18 = y%d2val2*y%val
      q19 = pow2(y%d1val2)
      q20 = 2.0_dp*q19 - q18
      q21 = powm1(pow4(y%val))
      q22 = 2.0_dp*y%d1val1_d1val2
      binary%val = x%val*powm1(y%val)
      binary%d1val1 = q1*(q2 - q3)
      binary%d1val2 = q1*(q4 - q5)
      binary%d2val1 = q7*(q13*x%val - q2*q9 + q8)
      binary%d1val1_d1val2 = q7*(q14 + q15*q3 - y%val*(q16 + x%d1val1*y%d1val2 + x%val*y%d1val1_d1val2))
      binary%d2val2 = q7*(-q15*q4 + q17 + q20*x%val)
      binary%d3val1 = q21*(-3.0_dp*q8*y%d1val1 + 3.0_dp*q13*q2 + q6*x%d3val1 - x%val*(-6.0_dp*q10*y%d1val1 + 6.0_dp*pow3(y%d1val1) + q0*y%d3val1))
      binary%d2val1_d1val2 = q21*(-6.0_dp*q11*q5 + 2.0_dp*q10*q5 + 4.0_dp*q2*y%d1val1*y%d1val2 + 4.0_dp*q3*y%d1val1_d1val2*y%val - q0*q22*x%d1val1 - q0*x%d1val2*y%d2val1 - q0*x%val*y%d2val1_d1val2 + q12*q4 - q14*q9 + q6*x%d2val1_d1val2 - q8*y%d1val2)
      binary%d1val1_d2val2 = q21*(-6.0_dp*q19*q3 + 2.0_dp*y%val*(q15*q16 + q19*x%d1val1 + q22*q5 + q3*y%d2val2) - q0*(q15*x%d1val1_d1val2 + q22*x%d1val2 + x%d1val1*y%d2val2 + x%d2val2*y%d1val1 + x%val*y%d1val1_d2val2) + q6*x%d1val1_d2val2)
      binary%d3val2 = q21*(-3.0_dp*q17*y%d1val2 + 3.0_dp*q20*q4 + q6*x%d3val2 - x%val*(-6.0_dp*q18*y%d1val2 + 6.0_dp*pow3(y%d1val2) + q0*y%d3val2))
   end function div_self
   
   function div_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q0
      q0 = powm1(y)
      unary%val = q0*x%val
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function div_self_real
   
   function div_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = z*powm1(q0)
      q2 = x%d2val1*x%val
      q3 = pow2(x%d1val1)
      q4 = z*powm1(pow3(x%val))
      q5 = 2.0_dp*x%d1val1
      q6 = x%d1val1_d1val2*x%val
      q7 = x%d2val2*x%val
      q8 = -q7
      q9 = pow2(x%d1val2)
      q10 = z*powm1(pow4(x%val))
      q11 = 4.0_dp*q6
      q12 = 6.0_dp*x%d1val2
      unary%val = z*powm1(x%val)
      unary%d1val1 = -q1*x%d1val1
      unary%d1val2 = -q1*x%d1val2
      unary%d2val1 = q4*(2.0_dp*q3 - q2)
      unary%d1val1_d1val2 = q4*(q5*x%d1val2 - q6)
      unary%d2val2 = q4*(2.0_dp*q9 + q8)
      unary%d3val1 = q10*(-6.0_dp*pow3(x%d1val1) + 6.0_dp*q2*x%d1val1 - q0*x%d3val1)
      unary%d2val1_d1val2 = q10*(2.0_dp*q2*x%d1val2 - q0*x%d2val1_d1val2 + q11*x%d1val1 - q12*q3)
      unary%d1val1_d2val2 = -q10*(q0*x%d1val1_d2val2 - q11*x%d1val2 + q5*(3.0_dp*q9 + q8))
      unary%d3val2 = q10*(-6.0_dp*pow3(x%d1val2) - q0*x%d3val2 + q12*q7)
   end function div_real_self
   
   function div_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = y
      q0 = powm1(y_dp)
      unary%val = q0*x%val
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function div_self_int
   
   function div_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      y_dp = z
      q0 = pow2(x%val)
      q1 = y_dp*powm1(q0)
      q2 = x%d2val1*x%val
      q3 = pow2(x%d1val1)
      q4 = y_dp*powm1(pow3(x%val))
      q5 = 2.0_dp*x%d1val1
      q6 = x%d1val1_d1val2*x%val
      q7 = x%d2val2*x%val
      q8 = -q7
      q9 = pow2(x%d1val2)
      q10 = y_dp*powm1(pow4(x%val))
      q11 = 4.0_dp*q6
      q12 = 6.0_dp*x%d1val2
      unary%val = y_dp*powm1(x%val)
      unary%d1val1 = -q1*x%d1val1
      unary%d1val2 = -q1*x%d1val2
      unary%d2val1 = q4*(2.0_dp*q3 - q2)
      unary%d1val1_d1val2 = q4*(q5*x%d1val2 - q6)
      unary%d2val2 = q4*(2.0_dp*q9 + q8)
      unary%d3val1 = q10*(-6.0_dp*pow3(x%d1val1) + 6.0_dp*q2*x%d1val1 - q0*x%d3val1)
      unary%d2val1_d1val2 = q10*(2.0_dp*q2*x%d1val2 - q0*x%d2val1_d1val2 + q11*x%d1val1 - q12*q3)
      unary%d1val1_d2val2 = -q10*(q0*x%d1val1_d2val2 - q11*x%d1val2 + q5*(3.0_dp*q9 + q8))
      unary%d3val2 = q10*(-6.0_dp*pow3(x%d1val2) - q0*x%d3val2 + q12*q7)
   end function div_int_self
   
   function pow_self(x, y) result(binary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3), intent(in) :: y
      type(auto_diff_real_2var_order3) :: binary
      real(dp) :: q30
      real(dp) :: q29
      real(dp) :: q28
      real(dp) :: q27
      real(dp) :: q26
      real(dp) :: q25
      real(dp) :: q24
      real(dp) :: q23
      real(dp) :: q22
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow(x%val, y%val - 1)
      q1 = x%d1val1*y%val
      q2 = log(x%val)
      q3 = q2*x%val
      q4 = q1 + q3*y%d1val1
      q5 = x%d1val2*y%val
      q6 = q3*y%d1val2 + q5
      q7 = pow(x%val, -2.0_dp + y%val)
      q8 = pow2(x%d1val1)
      q9 = pow2(x%val)
      q10 = q2*q9
      q11 = x%d2val1*y%val
      q12 = x%d1val1*y%d1val1
      q13 = q10*y%d2val1 - q8*y%val + x%val*(2.0_dp*q12 + q11)
      q14 = q13 + pow2(q4)
      q15 = x%d1val1*y%d1val2
      q16 = x%d1val2*y%d1val1
      q17 = -q1*x%d1val2 + q10*y%d1val1_d1val2 + x%val*(q15 + q16 + x%d1val1_d1val2*y%val)
      q18 = pow2(x%d1val2)
      q19 = x%d2val2*y%val
      q20 = x%d1val2*y%d1val2
      q21 = q10*y%d2val2 - q18*y%val + x%val*(2.0_dp*q20 + q19)
      q22 = q21 + pow2(q6)
      q23 = pow(x%val, -3.0_dp + y%val)
      q24 = 2.0_dp*y%val
      q25 = q2*pow3(x%val)
      q26 = 3.0_dp*x%d1val1
      q27 = 2.0_dp*y%d1val1_d1val2
      q28 = 2.0_dp*x%d1val1_d1val2
      q29 = 2.0_dp*q17
      q30 = 3.0_dp*x%d1val2
      binary%val = pow(x%val, y%val)
      binary%d1val1 = q0*q4
      binary%d1val2 = q0*q6
      binary%d2val1 = q14*q7
      binary%d1val1_d1val2 = (q17 + q4*q6)*pow(x%val, 3.0_dp + y%val)*powm1(pow5(x%val))
      binary%d2val2 = q22*q7
      binary%d3val1 = q23*(3.0_dp*q13*q4 + q24*pow3(x%d1val1) + q25*y%d3val1 - q26*x%val*(q11 + q12) + q9*(3.0_dp*x%d2val1*y%d1val1 + q26*y%d2val1 + x%d3val1*y%val) + pow3(q4))
      binary%d2val1_d1val2 = (2.0_dp*q5*q8 + q14*q6 + q25*y%d2val1_d1val2 + q29*q4 + q9*(q27*x%d1val1 + q28*y%d1val1 + x%d1val2*y%d2val1 + x%d2val1*y%d1val2 + x%d2val1_d1val2*y%val) - x%val*(2.0_dp*q16*x%d1val1 + q1*q28 + q5*x%d2val1 + q8*y%d1val2))*pow(x%val, 10.0_dp + y%val)*powm1(pow(x%val, 13.0_dp))
      binary%d1val1_d2val2 = q23*(2.0_dp*q1*q18 + q22*q4 + q25*y%d1val1_d2val2 + q29*q6 + q9*(q27*x%d1val2 + q28*y%d1val2 + x%d1val1*y%d2val2 + x%d1val1_d2val2*y%val + x%d2val2*y%d1val1) - x%val*(2.0_dp*q15*x%d1val2 + q1*x%d2val2 + q18*y%d1val1 + q28*q5))
      binary%d3val2 = q23*(3.0_dp*q21*q6 + q24*pow3(x%d1val2) + q25*y%d3val2 - q30*x%val*(q19 + q20) + q9*(3.0_dp*x%d2val2*y%d1val2 + q30*y%d2val2 + x%d3val2*y%val) + pow3(q6))
   end function pow_self
   
   function pow_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = y - 1
      q1 = y*pow(x%val, q0)
      q2 = x%d2val1*x%val
      q3 = pow2(x%d1val1)
      q4 = q3*y
      q5 = pow(x%val, -2.0_dp + y)
      q6 = q5*y
      q7 = x%d1val1_d1val2*x%val
      q8 = x%d1val1*x%d1val2
      q9 = x%d2val2*x%val
      q10 = pow2(x%d1val2)
      q11 = q10*y
      q12 = y*(-q10 + q11 + q9)
      q13 = pow2(x%val)
      q14 = pow2(y)
      q15 = -3.0_dp*y + 2.0_dp + q14
      q16 = y*pow(x%val, -3.0_dp + y)
      q17 = 2.0_dp*q7
      q18 = q17*x%d1val1
      q19 = q2*x%d1val2
      q20 = q3*x%d1val2
      q21 = 3.0_dp*x%d1val2
      unary%val = pow(x%val, y)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q6*(q2 - q3 + q4)
      unary%d1val1_d1val2 = q6*(q7 + q8*y - q8)
      unary%d2val2 = q12*q5
      unary%d3val1 = q16*(3.0_dp*q0*q2*x%d1val1 + q13*x%d3val1 + q15*pow3(x%d1val1))
      unary%d2val1_d1val2 = q16*(2.0_dp*q20 + q13*x%d2val1_d1val2 + q14*q20 + q18*y - q18 + q19*y - q19 - q21*q4)
      unary%d1val1_d2val2 = q16*(q0*q17*x%d1val2 + q13*x%d1val1_d2val2 - x%d1val1*(-2.0_dp*q10 + 2.0_dp*q11 - q12 + q9))
      unary%d3val2 = q16*(q0*q21*q9 + q13*x%d3val2 + q15*pow3(x%d1val2))
   end function pow_self_real
   
   function pow_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow(z, x%val)
      q1 = log(z)
      q2 = q0*q1
      q3 = q1*pow2(x%d1val1) + x%d2val1
      q4 = q1*x%d1val2
      q5 = q1*pow2(x%d1val2) + x%d2val2
      q6 = q1*x%d1val1
      q7 = pow2(q1)
      q8 = 2.0_dp*x%d1val1_d1val2
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q2*(q4*x%d1val1 + x%d1val1_d1val2)
      unary%d2val2 = q2*q5
      unary%d3val1 = q2*(3.0_dp*q6*x%d2val1 + q7*pow3(x%d1val1) + x%d3val1)
      unary%d2val1_d1val2 = q2*(q3*q4 + q6*q8 + x%d2val1_d1val2)
      unary%d1val1_d2val2 = q2*(q4*q8 + q5*q6 + x%d1val1_d2val2)
      unary%d3val2 = q2*(3.0_dp*q4*x%d2val2 + q7*pow3(x%d1val2) + x%d3val2)
   end function pow_real_self
   
   function pow_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      y_dp = y
      q0 = y_dp - 1
      q1 = y_dp*pow(x%val, q0)
      q2 = x%d2val1*x%val
      q3 = pow2(x%d1val1)
      q4 = q3*y_dp
      q5 = pow(x%val, -2.0_dp + y_dp)
      q6 = q5*y_dp
      q7 = x%d1val1_d1val2*x%val
      q8 = x%d1val1*x%d1val2
      q9 = x%d2val2*x%val
      q10 = pow2(x%d1val2)
      q11 = q10*y_dp
      q12 = y_dp*(-q10 + q11 + q9)
      q13 = pow2(x%val)
      q14 = pow2(y_dp)
      q15 = -3.0_dp*y_dp + 2.0_dp + q14
      q16 = y_dp*pow(x%val, -3.0_dp + y_dp)
      q17 = 2.0_dp*q7
      q18 = q17*x%d1val1
      q19 = q2*x%d1val2
      q20 = q3*x%d1val2
      q21 = 3.0_dp*x%d1val2
      unary%val = pow(x%val, y_dp)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q6*(q2 - q3 + q4)
      unary%d1val1_d1val2 = q6*(q7 + q8*y_dp - q8)
      unary%d2val2 = q12*q5
      unary%d3val1 = q16*(3.0_dp*q0*q2*x%d1val1 + q13*x%d3val1 + q15*pow3(x%d1val1))
      unary%d2val1_d1val2 = q16*(2.0_dp*q20 + q13*x%d2val1_d1val2 + q14*q20 + q18*y_dp - q18 + q19*y_dp - q19 - q21*q4)
      unary%d1val1_d2val2 = q16*(q0*q17*x%d1val2 + q13*x%d1val1_d2val2 - x%d1val1*(-2.0_dp*q10 + 2.0_dp*q11 - q12 + q9))
      unary%d3val2 = q16*(q0*q21*q9 + q13*x%d3val2 + q15*pow3(x%d1val2))
   end function pow_self_int
   
   function pow_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      y_dp = z
      q0 = pow(y_dp, x%val)
      q1 = log(y_dp)
      q2 = q0*q1
      q3 = q1*pow2(x%d1val1) + x%d2val1
      q4 = q1*x%d1val2
      q5 = q1*pow2(x%d1val2) + x%d2val2
      q6 = q1*x%d1val1
      q7 = pow2(q1)
      q8 = 2.0_dp*x%d1val1_d1val2
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q2*(q4*x%d1val1 + x%d1val1_d1val2)
      unary%d2val2 = q2*q5
      unary%d3val1 = q2*(3.0_dp*q6*x%d2val1 + q7*pow3(x%d1val1) + x%d3val1)
      unary%d2val1_d1val2 = q2*(q3*q4 + q6*q8 + x%d2val1_d1val2)
      unary%d1val1_d2val2 = q2*(q4*q8 + q5*q6 + x%d1val1_d2val2)
      unary%d3val2 = q2*(3.0_dp*q4*x%d2val2 + q7*pow3(x%d1val2) + x%d3val2)
   end function pow_int_self
   
   function max_self(x, y) result(binary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3), intent(in) :: y
      type(auto_diff_real_2var_order3) :: binary
      real(dp) :: q1
      real(dp) :: q0
      q0 = Heaviside(x%val - y%val)
      q1 = Heaviside(-x%val + y%val)
      binary%val = Max(x%val, y%val)
      binary%d1val1 = q0*x%d1val1 + q1*y%d1val1
      binary%d1val2 = q0*x%d1val2 + q1*y%d1val2
      binary%d2val1 = q0*x%d2val1 + q1*y%d2val1
      binary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q1*y%d1val1_d1val2
      binary%d2val2 = q0*x%d2val2 + q1*y%d2val2
      binary%d3val1 = q0*x%d3val1 + q1*y%d3val1
      binary%d2val1_d1val2 = q0*x%d2val1_d1val2 + q1*y%d2val1_d1val2
      binary%d1val1_d2val2 = q0*x%d1val1_d2val2 + q1*y%d1val1_d2val2
      binary%d3val2 = q0*x%d3val2 + q1*y%d3val2
   end function max_self
   
   function max_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q0
      q0 = Heaviside(x%val - y)
      unary%val = Max(x%val, y)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function max_self_real
   
   function max_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q0
      q0 = Heaviside(x%val - z)
      unary%val = Max(x%val, z)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function max_real_self
   
   function max_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = y
      q0 = Heaviside(x%val - y_dp)
      unary%val = Max(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function max_self_int
   
   function max_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = z
      q0 = Heaviside(x%val - y_dp)
      unary%val = Max(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function max_int_self
   
   function min_self(x, y) result(binary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3), intent(in) :: y
      type(auto_diff_real_2var_order3) :: binary
      real(dp) :: q1
      real(dp) :: q0
      q0 = Heaviside(-x%val + y%val)
      q1 = Heaviside(x%val - y%val)
      binary%val = Min(x%val, y%val)
      binary%d1val1 = q0*x%d1val1 + q1*y%d1val1
      binary%d1val2 = q0*x%d1val2 + q1*y%d1val2
      binary%d2val1 = q0*x%d2val1 + q1*y%d2val1
      binary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q1*y%d1val1_d1val2
      binary%d2val2 = q0*x%d2val2 + q1*y%d2val2
      binary%d3val1 = q0*x%d3val1 + q1*y%d3val1
      binary%d2val1_d1val2 = q0*x%d2val1_d1val2 + q1*y%d2val1_d1val2
      binary%d1val1_d2val2 = q0*x%d1val1_d2val2 + q1*y%d1val1_d2val2
      binary%d3val2 = q0*x%d3val2 + q1*y%d3val2
   end function min_self
   
   function min_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q0
      q0 = Heaviside(-x%val + y)
      unary%val = Min(x%val, y)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function min_self_real
   
   function min_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q0
      q0 = Heaviside(-x%val + z)
      unary%val = Min(x%val, z)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function min_real_self
   
   function min_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = y
      q0 = Heaviside(-x%val + y_dp)
      unary%val = Min(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function min_self_int
   
   function min_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = z
      q0 = Heaviside(-x%val + y_dp)
      unary%val = Min(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d3val2 = q0*x%d3val2
   end function min_int_self
   
   function dim_self(x, y) result(binary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3), intent(in) :: y
      type(auto_diff_real_2var_order3) :: binary
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val - y%val
      q1 = 0.5_dp*sgn(q0)
      binary%val = -0.5_dp*y%val + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      binary%d1val1 = -0.5_dp*y%d1val1 + 0.5_dp*x%d1val1 + q1*(x%d1val1 - y%d1val1)
      binary%d1val2 = -0.5_dp*y%d1val2 + 0.5_dp*x%d1val2 + q1*(x%d1val2 - y%d1val2)
      binary%d2val1 = -0.5_dp*y%d2val1 + 0.5_dp*x%d2val1 + q1*(x%d2val1 - y%d2val1)
      binary%d1val1_d1val2 = -0.5_dp*y%d1val1_d1val2 + 0.5_dp*x%d1val1_d1val2 + q1*(x%d1val1_d1val2 - y%d1val1_d1val2)
      binary%d2val2 = -0.5_dp*y%d2val2 + 0.5_dp*x%d2val2 + q1*(x%d2val2 - y%d2val2)
      binary%d3val1 = -0.5_dp*y%d3val1 + 0.5_dp*x%d3val1 + q1*(x%d3val1 - y%d3val1)
      binary%d2val1_d1val2 = -0.5_dp*y%d2val1_d1val2 + 0.5_dp*x%d2val1_d1val2 + q1*(x%d2val1_d1val2 - y%d2val1_d1val2)
      binary%d1val1_d2val2 = -0.5_dp*y%d1val1_d2val2 + 0.5_dp*x%d1val1_d2val2 + q1*(x%d1val1_d2val2 - y%d1val1_d2val2)
      binary%d3val2 = -0.5_dp*y%d3val2 + 0.5_dp*x%d3val2 + q1*(x%d3val2 - y%d3val2)
   end function dim_self
   
   function dim_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val - y
      q1 = sgn(q0)
      q2 = 0.5_dp*q1 + 0.5_dp
      q3 = 0.5_dp*x%d1val1_d1val2
      q4 = 0.5_dp*x%d2val1_d1val2
      q5 = 0.5_dp*x%d1val1_d2val2
      unary%val = -0.5_dp*y + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*x%d2val1
      unary%d1val1_d1val2 = q1*q3 + q3
      unary%d2val2 = q2*x%d2val2
      unary%d3val1 = q2*x%d3val1
      unary%d2val1_d1val2 = q1*q4 + q4
      unary%d1val1_d2val2 = q1*q5 + q5
      unary%d3val2 = q2*x%d3val2
   end function dim_self_real
   
   function dim_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val - z
      q1 = sgn(q0)
      q2 = -0.5_dp + 0.5_dp*q1
      q3 = 0.5_dp*x%d1val1_d1val2
      q4 = 0.5_dp*x%d2val1_d1val2
      q5 = 0.5_dp*x%d1val1_d2val2
      unary%val = -0.5_dp*x%val + 0.5_dp*z + 0.5_dp*Abs(q0)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*x%d2val1
      unary%d1val1_d1val2 = q1*q3 - q3
      unary%d2val2 = q2*x%d2val2
      unary%d3val1 = q2*x%d3val1
      unary%d2val1_d1val2 = q1*q4 - q4
      unary%d1val1_d2val2 = q1*q5 - q5
      unary%d3val2 = q2*x%d3val2
   end function dim_real_self
   
   function dim_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      y_dp = y
      q0 = x%val - y_dp
      q1 = sgn(q0)
      q2 = 0.5_dp*q1 + 0.5_dp
      q3 = 0.5_dp*x%d1val1_d1val2
      q4 = 0.5_dp*x%d2val1_d1val2
      q5 = 0.5_dp*x%d1val1_d2val2
      unary%val = -0.5_dp*y_dp + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*x%d2val1
      unary%d1val1_d1val2 = q1*q3 + q3
      unary%d2val2 = q2*x%d2val2
      unary%d3val1 = q2*x%d3val1
      unary%d2val1_d1val2 = q1*q4 + q4
      unary%d1val1_d2val2 = q1*q5 + q5
      unary%d3val2 = q2*x%d3val2
   end function dim_self_int
   
   function dim_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3), intent(in) :: x
      type(auto_diff_real_2var_order3) :: unary
      real(dp) :: y_dp
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      y_dp = z
      q0 = x%val - y_dp
      q1 = sgn(q0)
      q2 = -0.5_dp + 0.5_dp*q1
      q3 = 0.5_dp*x%d1val1_d1val2
      q4 = 0.5_dp*x%d2val1_d1val2
      q5 = 0.5_dp*x%d1val1_d2val2
      unary%val = -0.5_dp*x%val + 0.5_dp*y_dp + 0.5_dp*Abs(q0)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*x%d2val1
      unary%d1val1_d1val2 = q1*q3 - q3
      unary%d2val2 = q2*x%d2val2
      unary%d3val1 = q2*x%d3val1
      unary%d2val1_d1val2 = q1*q4 - q4
      unary%d1val1_d2val2 = q1*q5 - q5
      unary%d3val2 = q2*x%d3val2
   end function dim_int_self
   
   function differentiate_auto_diff_real_2var_order3_1(this) result(derivative)
      type(auto_diff_real_2var_order3), intent(in) :: this
      type(auto_diff_real_2var_order3) :: derivative
      derivative%val = this%d1val1
      derivative%d1val1 = this%d2val1
      derivative%d1val2 = this%d1val1_d1val2
      derivative%d2val1 = this%d3val1
      derivative%d1val1_d1val2 = this%d2val1_d1val2
      derivative%d2val2 = this%d1val1_d2val2
      derivative%d3val1 = 0_dp
      derivative%d2val1_d1val2 = 0_dp
      derivative%d1val1_d2val2 = 0_dp
      derivative%d3val2 = 0_dp
   end function differentiate_auto_diff_real_2var_order3_1
   
   function differentiate_auto_diff_real_2var_order3_2(this) result(derivative)
      type(auto_diff_real_2var_order3), intent(in) :: this
      type(auto_diff_real_2var_order3) :: derivative
      derivative%val = this%d1val2
      derivative%d1val1 = this%d1val1_d1val2
      derivative%d1val2 = this%d2val2
      derivative%d2val1 = this%d2val1_d1val2
      derivative%d1val1_d1val2 = this%d1val1_d2val2
      derivative%d2val2 = this%d3val2
      derivative%d3val1 = 0_dp
      derivative%d2val1_d1val2 = 0_dp
      derivative%d1val1_d2val2 = 0_dp
      derivative%d3val2 = 0_dp
   end function differentiate_auto_diff_real_2var_order3_2
   
end module auto_diff_real_2var_order3_module