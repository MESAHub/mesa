module auto_diff_real_2var_order3_array_module
      use const_def
      use utils_lib
      use support_functions
      use math_lib
   
      implicit none
      private
   public :: auto_diff_real_2var_order3_array, &
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
      operator(**), &
      max, &
      min, &
      dim, &
      pow, &
      differentiate_1, &
      differentiate_2
   type :: auto_diff_real_2var_order3_array
      real(dp) :: val
      real(dp) :: d1val1
      real(dp) :: d1val2
      real(dp) :: d1Array(5)
      real(dp) :: d2val1
      real(dp) :: d1val1_d1val2
      real(dp) :: d1val1_d1Array(5)
      real(dp) :: d2val2
      real(dp) :: d1val2_d1Array(5)
      real(dp) :: d3val1
      real(dp) :: d2val1_d1val2
      real(dp) :: d2val1_d1Array(5)
      real(dp) :: d1val1_d2val2
      real(dp) :: d1val1_d1val2_d1Array(5)
      real(dp) :: d3val2
      real(dp) :: d2val2_d1Array(5)
   end type auto_diff_real_2var_order3_array
   
   interface assignment(=)
      module procedure assign_from_self
      module procedure assign_from_real_dp
      module procedure assign_from_int
   end interface assignment(=)
   
   interface operator(.eq.)
      module procedure equal_self
      module procedure equal_auto_diff_real_2var_order3_array_real_dp
      module procedure equal_real_dp_auto_diff_real_2var_order3_array
      module procedure equal_auto_diff_real_2var_order3_array_int
      module procedure equal_int_auto_diff_real_2var_order3_array
   end interface operator(.eq.)
   
   interface operator(.ne.)
      module procedure neq_self
      module procedure neq_auto_diff_real_2var_order3_array_real_dp
      module procedure neq_real_dp_auto_diff_real_2var_order3_array
      module procedure neq_auto_diff_real_2var_order3_array_int
      module procedure neq_int_auto_diff_real_2var_order3_array
   end interface operator(.ne.)
   
   interface operator(.gt.)
      module procedure greater_self
      module procedure greater_auto_diff_real_2var_order3_array_real_dp
      module procedure greater_real_dp_auto_diff_real_2var_order3_array
      module procedure greater_auto_diff_real_2var_order3_array_int
      module procedure greater_int_auto_diff_real_2var_order3_array
   end interface operator(.gt.)
   
   interface operator(.lt.)
      module procedure less_self
      module procedure less_auto_diff_real_2var_order3_array_real_dp
      module procedure less_real_dp_auto_diff_real_2var_order3_array
      module procedure less_auto_diff_real_2var_order3_array_int
      module procedure less_int_auto_diff_real_2var_order3_array
   end interface operator(.lt.)
   
   interface operator(.le.)
      module procedure leq_self
      module procedure leq_auto_diff_real_2var_order3_array_real_dp
      module procedure leq_real_dp_auto_diff_real_2var_order3_array
      module procedure leq_auto_diff_real_2var_order3_array_int
      module procedure leq_int_auto_diff_real_2var_order3_array
   end interface operator(.le.)
   
   interface operator(.ge.)
      module procedure geq_self
      module procedure geq_auto_diff_real_2var_order3_array_real_dp
      module procedure geq_real_dp_auto_diff_real_2var_order3_array
      module procedure geq_auto_diff_real_2var_order3_array_int
      module procedure geq_int_auto_diff_real_2var_order3_array
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
   
   interface operator(**)
      module procedure pow_self
      module procedure pow_self_real
      module procedure pow_real_self
      module procedure pow_self_int
      module procedure pow_int_self
   end interface operator(**)
   
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
   
   interface pow
      module procedure pow_self
      module procedure pow_self_real
      module procedure pow_real_self
      module procedure pow_self_int
      module procedure pow_int_self
   end interface pow
   
   interface differentiate_1
      module procedure differentiate_auto_diff_real_2var_order3_array_1
   end interface differentiate_1
   
   interface differentiate_2
      module procedure differentiate_auto_diff_real_2var_order3_array_2
   end interface differentiate_2
   
   contains

   subroutine assign_from_self(this, other)
      type(auto_diff_real_2var_order3_array), intent(out) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      this%val = other%val
      this%d1val1 = other%d1val1
      this%d1val2 = other%d1val2
      this%d1Array = other%d1Array
      this%d2val1 = other%d2val1
      this%d1val1_d1val2 = other%d1val1_d1val2
      this%d1val1_d1Array = other%d1val1_d1Array
      this%d2val2 = other%d2val2
      this%d1val2_d1Array = other%d1val2_d1Array
      this%d3val1 = other%d3val1
      this%d2val1_d1val2 = other%d2val1_d1val2
      this%d2val1_d1Array = other%d2val1_d1Array
      this%d1val1_d2val2 = other%d1val1_d2val2
      this%d1val1_d1val2_d1Array = other%d1val1_d1val2_d1Array
      this%d3val2 = other%d3val2
      this%d2val2_d1Array = other%d2val2_d1Array
   end subroutine assign_from_self
   
   subroutine assign_from_real_dp(this, other)
      type(auto_diff_real_2var_order3_array), intent(out) :: this
      real(dp), intent(in) :: other
      this%val = other
      this%d1val1 = 0_dp
      this%d1val2 = 0_dp
      this%d1Array = 0_dp
      this%d2val1 = 0_dp
      this%d1val1_d1val2 = 0_dp
      this%d1val1_d1Array = 0_dp
      this%d2val2 = 0_dp
      this%d1val2_d1Array = 0_dp
      this%d3val1 = 0_dp
      this%d2val1_d1val2 = 0_dp
      this%d2val1_d1Array = 0_dp
      this%d1val1_d2val2 = 0_dp
      this%d1val1_d1val2_d1Array = 0_dp
      this%d3val2 = 0_dp
      this%d2val2_d1Array = 0_dp
   end subroutine assign_from_real_dp
   
   subroutine assign_from_int(this, other)
      type(auto_diff_real_2var_order3_array), intent(out) :: this
      integer, intent(in) :: other
      this%val = other
      this%d1val1 = 0_dp
      this%d1val2 = 0_dp
      this%d1Array = 0_dp
      this%d2val1 = 0_dp
      this%d1val1_d1val2 = 0_dp
      this%d1val1_d1Array = 0_dp
      this%d2val2 = 0_dp
      this%d1val2_d1Array = 0_dp
      this%d3val1 = 0_dp
      this%d2val1_d1val2 = 0_dp
      this%d2val1_d1Array = 0_dp
      this%d1val1_d2val2 = 0_dp
      this%d1val1_d1val2_d1Array = 0_dp
      this%d3val2 = 0_dp
      this%d2val2_d1Array = 0_dp
   end subroutine assign_from_int
   
   function equal_self(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this%val .eq. other%val)
   end function equal_self
   
   function equal_auto_diff_real_2var_order3_array_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .eq. other)
   end function equal_auto_diff_real_2var_order3_array_real_dp
   
   function equal_real_dp_auto_diff_real_2var_order3_array(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .eq. other%val)
   end function equal_real_dp_auto_diff_real_2var_order3_array
   
   function equal_auto_diff_real_2var_order3_array_int(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .eq. other)
   end function equal_auto_diff_real_2var_order3_array_int
   
   function equal_int_auto_diff_real_2var_order3_array(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .eq. other%val)
   end function equal_int_auto_diff_real_2var_order3_array
   
   function neq_self(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this%val .ne. other%val)
   end function neq_self
   
   function neq_auto_diff_real_2var_order3_array_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .ne. other)
   end function neq_auto_diff_real_2var_order3_array_real_dp
   
   function neq_real_dp_auto_diff_real_2var_order3_array(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .ne. other%val)
   end function neq_real_dp_auto_diff_real_2var_order3_array
   
   function neq_auto_diff_real_2var_order3_array_int(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .ne. other)
   end function neq_auto_diff_real_2var_order3_array_int
   
   function neq_int_auto_diff_real_2var_order3_array(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .ne. other%val)
   end function neq_int_auto_diff_real_2var_order3_array
   
   function greater_self(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this%val .gt. other%val)
   end function greater_self
   
   function greater_auto_diff_real_2var_order3_array_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .gt. other)
   end function greater_auto_diff_real_2var_order3_array_real_dp
   
   function greater_real_dp_auto_diff_real_2var_order3_array(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .gt. other%val)
   end function greater_real_dp_auto_diff_real_2var_order3_array
   
   function greater_auto_diff_real_2var_order3_array_int(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .gt. other)
   end function greater_auto_diff_real_2var_order3_array_int
   
   function greater_int_auto_diff_real_2var_order3_array(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .gt. other%val)
   end function greater_int_auto_diff_real_2var_order3_array
   
   function less_self(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this%val .lt. other%val)
   end function less_self
   
   function less_auto_diff_real_2var_order3_array_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .lt. other)
   end function less_auto_diff_real_2var_order3_array_real_dp
   
   function less_real_dp_auto_diff_real_2var_order3_array(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .lt. other%val)
   end function less_real_dp_auto_diff_real_2var_order3_array
   
   function less_auto_diff_real_2var_order3_array_int(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .lt. other)
   end function less_auto_diff_real_2var_order3_array_int
   
   function less_int_auto_diff_real_2var_order3_array(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .lt. other%val)
   end function less_int_auto_diff_real_2var_order3_array
   
   function leq_self(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this%val .le. other%val)
   end function leq_self
   
   function leq_auto_diff_real_2var_order3_array_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .le. other)
   end function leq_auto_diff_real_2var_order3_array_real_dp
   
   function leq_real_dp_auto_diff_real_2var_order3_array(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .le. other%val)
   end function leq_real_dp_auto_diff_real_2var_order3_array
   
   function leq_auto_diff_real_2var_order3_array_int(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .le. other)
   end function leq_auto_diff_real_2var_order3_array_int
   
   function leq_int_auto_diff_real_2var_order3_array(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .le. other%val)
   end function leq_int_auto_diff_real_2var_order3_array
   
   function geq_self(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this%val .ge. other%val)
   end function geq_self
   
   function geq_auto_diff_real_2var_order3_array_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .ge. other)
   end function geq_auto_diff_real_2var_order3_array_real_dp
   
   function geq_real_dp_auto_diff_real_2var_order3_array(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .ge. other%val)
   end function geq_real_dp_auto_diff_real_2var_order3_array
   
   function geq_auto_diff_real_2var_order3_array_int(this, other) result(z)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .ge. other)
   end function geq_auto_diff_real_2var_order3_array_int
   
   function geq_int_auto_diff_real_2var_order3_array(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_array), intent(in) :: other
      logical :: z
      z = (this .ge. other%val)
   end function geq_int_auto_diff_real_2var_order3_array
   
   function make_unary_operator(x, z_val, z_d1x, z_d2x, z_d3x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      real(dp), intent(in) :: z_val
      real(dp), intent(in) :: z_d1x
      real(dp), intent(in) :: z_d2x
      real(dp), intent(in) :: z_d3x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6(1:5)
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%d1val1)
      q1 = x%d1val1*z_d2x
      q2 = pow2(x%d1val2)
      q3 = x%d1val2*z_d2x
      q4 = 2.0_dp*q1
      q5 = q0*z_d3x
      q6(1:5) = x%d1Array(1:5)*z_d2x
      q7 = 2.0_dp*q3
      q8 = q2*z_d3x
      unary%val = z_val
      unary%d1val1 = x%d1val1*z_d1x
      unary%d1val2 = x%d1val2*z_d1x
      unary%d1Array(1:5) = x%d1Array(1:5)*z_d1x
      unary%d2val1 = q0*z_d2x + x%d2val1*z_d1x
      unary%d1val1_d1val2 = q1*x%d1val2 + x%d1val1_d1val2*z_d1x
      unary%d1val1_d1Array(1:5) = q1*x%d1Array(1:5) + x%d1val1_d1Array(1:5)*z_d1x
      unary%d2val2 = q2*z_d2x + x%d2val2*z_d1x
      unary%d1val2_d1Array(1:5) = q3*x%d1Array(1:5) + x%d1val2_d1Array(1:5)*z_d1x
      unary%d3val1 = 3.0_dp*q1*x%d2val1 + x%d3val1*z_d1x + z_d3x*pow3(x%d1val1)
      unary%d2val1_d1val2 = q3*x%d2val1 + q4*x%d1val1_d1val2 + q5*x%d1val2 + x%d2val1_d1val2*z_d1x
      unary%d2val1_d1Array(1:5) = q4*x%d1val1_d1Array(1:5) + q5*x%d1Array(1:5) + q6*x%d2val1 + x%d2val1_d1Array(1:5)*z_d1x
      unary%d1val1_d2val2 = q1*x%d2val2 + q7*x%d1val1_d1val2 + q8*x%d1val1 + x%d1val1_d2val2*z_d1x
      unary%d1val1_d1val2_d1Array(1:5) = q1*x%d1val2_d1Array(1:5) + q3*x%d1val1_d1Array(1:5) + q6*x%d1val1_d1val2 + x%d1Array(1:5)*x%d1val1*x%d1val2*z_d3x + x%d1val1_d1val2_d1Array(1:5)*z_d1x
      unary%d3val2 = 3.0_dp*q3*x%d2val2 + x%d3val2*z_d1x + z_d3x*pow3(x%d1val2)
      unary%d2val2_d1Array(1:5) = q6*x%d2val2 + q7*x%d1val2_d1Array(1:5) + q8*x%d1Array(1:5) + x%d2val2_d1Array(1:5)*z_d1x
   end function make_unary_operator
   
   function make_binary_operator(x, y, z_val, z_d1x, z_d1y, z_d2x, z_d1x_d1y, z_d2y, z_d3x, z_d2x_d1y, z_d1x_d2y, z_d3y) result(binary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array), intent(in) :: y
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
      type(auto_diff_real_2var_order3_array) :: binary
      real(dp) :: q45(1:5)
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42(1:5)
      real(dp) :: q41(1:5)
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33(1:5)
      real(dp) :: q32(1:5)
      real(dp) :: q31(1:5)
      real(dp) :: q30(1:5)
      real(dp) :: q29(1:5)
      real(dp) :: q28(1:5)
      real(dp) :: q27(1:5)
      real(dp) :: q26(1:5)
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
      q4 = x%d1val1*z_d2x
      q5 = y%d1val1*z_d1x_d1y
      q6 = y%d1val1*z_d2y
      q7 = pow2(x%d1val2)
      q8 = pow2(y%d1val2)
      q9 = y%d1val2*z_d1x_d1y
      q10 = 2.0_dp*q9
      q11 = x%d1val2*z_d2x
      q12 = x%d1val2*z_d1x_d1y
      q13 = y%d1val2*z_d2y
      q14 = 3.0_dp*x%d2val1
      q15 = 3.0_dp*y%d2val1
      q16 = q1*z_d1x_d2y
      q17 = q0*z_d2x_d1y
      q18 = 2.0_dp*x%d1val1_d1val2
      q19 = 2.0_dp*q6
      q20 = q0*z_d3x
      q21 = q1*z_d3y
      q22 = x%d1val2*z_d2x_d1y
      q23 = 2.0_dp*y%d1val1
      q24 = q23*x%d1val1
      q25 = y%d1val2*z_d1x_d2y
      q26(1:5) = x%d1Array(1:5)*z_d2x
      q27(1:5) = x%d1Array(1:5)*z_d1x_d1y
      q28(1:5) = y%d1Array(1:5)*z_d1x_d1y
      q29(1:5) = y%d1Array(1:5)*z_d2y
      q30(1:5) = 2.0_dp*x%d1val1_d1Array(1:5)
      q31(1:5) = x%d1Array(1:5)*x%d1val1
      q32(1:5) = q31*z_d2x_d1y
      q33(1:5) = x%d1val1*y%d1Array(1:5)
      q34 = 2.0_dp*y%d1val1_d1val2
      q35 = q7*z_d3x
      q36 = q8*z_d1x_d2y
      q37 = q7*z_d2x_d1y
      q38 = q8*z_d3y
      q39 = 2.0_dp*q22*y%d1val2
      q40 = q25*x%d1val2
      q41(1:5) = x%d1Array(1:5)*y%d1val1
      q42(1:5) = y%d1Array(1:5)*y%d1val1
      q43 = 3.0_dp*x%d2val2
      q44 = 3.0_dp*y%d2val2
      q45(1:5) = 2.0_dp*y%d1val2_d1Array(1:5)
      binary%val = z_val
      binary%d1val1 = x%d1val1*z_d1x + y%d1val1*z_d1y
      binary%d1val2 = x%d1val2*z_d1x + y%d1val2*z_d1y
      binary%d1Array(1:5) = x%d1Array(1:5)*z_d1x + y%d1Array(1:5)*z_d1y
      binary%d2val1 = q0*z_d2x + q1*z_d2y + q3*y%d1val1 + x%d2val1*z_d1x + y%d2val1*z_d1y
      binary%d1val1_d1val2 = q2*y%d1val2 + q4*x%d1val2 + q5*x%d1val2 + q6*y%d1val2 + x%d1val1_d1val2*z_d1x + y%d1val1_d1val2*z_d1y
      binary%d1val1_d1Array(1:5) = q2*y%d1Array(1:5) + q4*x%d1Array(1:5) + q5*x%d1Array(1:5) + q6*y%d1Array(1:5) + x%d1val1_d1Array(1:5)*z_d1x + y%d1val1_d1Array(1:5)*z_d1y
      binary%d2val2 = q10*x%d1val2 + q7*z_d2x + q8*z_d2y + x%d2val2*z_d1x + y%d2val2*z_d1y
      binary%d1val2_d1Array(1:5) = q11*x%d1Array(1:5) + q12*y%d1Array(1:5) + q13*y%d1Array(1:5) + q9*x%d1Array(1:5) + x%d1val2_d1Array(1:5)*z_d1x + y%d1val2_d1Array(1:5)*z_d1y
      binary%d3val1 = 3.0_dp*q16*x%d1val1 + 3.0_dp*q17*y%d1val1 + q14*q4 + q14*q5 + q15*q2 + q15*q6 + x%d3val1*z_d1x + y%d3val1*z_d1y + z_d3x*pow3(x%d1val1) + z_d3y*pow3(y%d1val1)
      binary%d2val1_d1val2 = q11*x%d2val1 + q12*y%d2val1 + q13*y%d2val1 + q16*x%d1val2 + q17*y%d1val2 + q18*q4 + q18*q5 + q19*y%d1val1_d1val2 + q20*x%d1val2 + q21*y%d1val2 + q22*q24 + q24*q25 + q3*y%d1val1_d1val2 + q9*x%d2val1 + x%d2val1_d1val2*z_d1x + y%d2val1_d1val2*z_d1y
      binary%d2val1_d1Array(1:5) = q16*x%d1Array(1:5) + q17*y%d1Array(1:5) + q19*y%d1val1_d1Array(1:5) + q20*x%d1Array(1:5) + q21*y%d1Array(1:5) + q23*q32 + q23*q33*z_d1x_d2y + q26*x%d2val1 + q27*y%d2val1 + q28*x%d2val1 + q29*y%d2val1 + q3*y%d1val1_d1Array(1:5) + q30*q4 + q30*q5 + x%d2val1_d1Array(1:5)*z_d1x + y%d2val1_d1Array(1:5)*z_d1y
      binary%d1val1_d2val2 = q10*x%d1val1_d1val2 + q11*q18 + q12*q34 + q13*q34 + q2*y%d2val2 + q23*q40 + q35*x%d1val1 + q36*x%d1val1 + q37*y%d1val1 + q38*y%d1val1 + q39*x%d1val1 + q4*x%d2val2 + q5*x%d2val2 + q6*y%d2val2 + x%d1val1_d2val2*z_d1x + y%d1val1_d2val2*z_d1y
      binary%d1val1_d1val2_d1Array(1:5) = q11*x%d1val1_d1Array(1:5) + q12*y%d1val1_d1Array(1:5) + q13*y%d1val1_d1Array(1:5) + q2*y%d1val2_d1Array(1:5) + q22*q33 + q22*q41 + q25*q33 + q25*q41 + q26*x%d1val1_d1val2 + q27*y%d1val1_d1val2 + q28*x%d1val1_d1val2 + q29*y%d1val1_d1val2 + q31*x%d1val2*z_d3x + q32*y%d1val2 + q4*x%d1val2_d1Array(1:5) + q42*x%d1val2*z_d1x_d2y + q42*y%d1val2*z_d3y + q5*x%d1val2_d1Array(1:5) + q6*y%d1val2_d1Array(1:5) + q9*x%d1val1_d1Array(1:5) + x%d1val1_d1val2_d1Array(1:5)*z_d1x + y%d1val1_d1val2_d1Array(1:5)*z_d1y
      binary%d3val2 = 3.0_dp*q36*x%d1val2 + 3.0_dp*q37*y%d1val2 + q11*q43 + q12*q44 + q13*q44 + q43*q9 + x%d3val2*z_d1x + y%d3val2*z_d1y + z_d3x*pow3(x%d1val2) + z_d3y*pow3(y%d1val2)
      binary%d2val2_d1Array(1:5) = 2.0_dp*q11*x%d1val2_d1Array(1:5) + 2.0_dp*q40*y%d1Array(1:5) + q10*x%d1val2_d1Array(1:5) + q12*q45 + q13*q45 + q26*x%d2val2 + q27*y%d2val2 + q28*x%d2val2 + q29*y%d2val2 + q35*x%d1Array(1:5) + q36*x%d1Array(1:5) + q37*y%d1Array(1:5) + q38*y%d1Array(1:5) + q39*x%d1Array(1:5) + x%d2val2_d1Array(1:5)*z_d1x + y%d2val2_d1Array(1:5)*z_d1y
   end function make_binary_operator
   
   function unary_minus_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      unary%val = -x%val
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d1Array(1:5) = -x%d1Array(1:5)
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = -x%d1val1_d1Array(1:5)
      unary%d2val2 = -x%d2val2
      unary%d1val2_d1Array(1:5) = -x%d1val2_d1Array(1:5)
      unary%d3val1 = -x%d3val1
      unary%d2val1_d1val2 = -x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = -x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = -x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = -x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = -x%d3val2
      unary%d2val2_d1Array(1:5) = -x%d2val2_d1Array(1:5)
   end function unary_minus_self
   
   function exp_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = exp(x%val)
      q1 = x%d2val1 + pow2(x%d1val1)
      q2 = x%d1val1*x%d1val2
      q3 = x%d2val2 + pow2(x%d1val2)
      q4 = 2.0_dp*x%d1val1
      q5 = 2.0_dp*x%d1val2
      unary%val = q0
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*q1
      unary%d1val1_d1val2 = q0*(q2 + x%d1val1_d1val2)
      unary%d1val1_d1Array(1:5) = q0*(x%d1Array(1:5)*x%d1val1 + x%d1val1_d1Array(1:5))
      unary%d2val2 = q0*q3
      unary%d1val2_d1Array(1:5) = q0*(x%d1Array(1:5)*x%d1val2 + x%d1val2_d1Array(1:5))
      unary%d3val1 = q0*(3.0_dp*x%d1val1*x%d2val1 + x%d3val1 + pow3(x%d1val1))
      unary%d2val1_d1val2 = q0*(q1*x%d1val2 + q4*x%d1val1_d1val2 + x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q0*(q1*x%d1Array(1:5) + q4*x%d1val1_d1Array(1:5) + x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q0*(q3*x%d1val1 + q5*x%d1val1_d1val2 + x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = q0*(q2*x%d1Array(1:5) + x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2 + x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q0*(3.0_dp*x%d1val2*x%d2val2 + x%d3val2 + pow3(x%d1val2))
      unary%d2val2_d1Array(1:5) = q0*(q3*x%d1Array(1:5) + q5*x%d1val2_d1Array(1:5) + x%d2val2_d1Array(1:5))
   end function exp_self
   
   function expm1_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = exp(x%val)
      q1 = x%d2val1 + pow2(x%d1val1)
      q2 = x%d1val1*x%d1val2
      q3 = x%d2val2 + pow2(x%d1val2)
      q4 = 2.0_dp*x%d1val1
      q5 = 2.0_dp*x%d1val2
      unary%val = expm1(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*q1
      unary%d1val1_d1val2 = q0*(q2 + x%d1val1_d1val2)
      unary%d1val1_d1Array(1:5) = q0*(x%d1Array(1:5)*x%d1val1 + x%d1val1_d1Array(1:5))
      unary%d2val2 = q0*q3
      unary%d1val2_d1Array(1:5) = q0*(x%d1Array(1:5)*x%d1val2 + x%d1val2_d1Array(1:5))
      unary%d3val1 = q0*(3.0_dp*x%d1val1*x%d2val1 + x%d3val1 + pow3(x%d1val1))
      unary%d2val1_d1val2 = q0*(q1*x%d1val2 + q4*x%d1val1_d1val2 + x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q0*(q1*x%d1Array(1:5) + q4*x%d1val1_d1Array(1:5) + x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q0*(q3*x%d1val1 + q5*x%d1val1_d1val2 + x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = q0*(q2*x%d1Array(1:5) + x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2 + x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q0*(3.0_dp*x%d1val2*x%d2val2 + x%d3val2 + pow3(x%d1val2))
      unary%d2val2_d1Array(1:5) = q0*(q3*x%d1Array(1:5) + q5*x%d1val2_d1Array(1:5) + x%d2val2_d1Array(1:5))
   end function expm1_self
   
   function exp10_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6(1:5)
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow(10.0_dp, x%val)
      q1 = ln10
      q2 = q0*q1
      q3 = pow2(x%d1val1)
      q4 = q1*x%d1val1
      q5 = pow2(x%d1val2)
      q6(1:5) = x%d1Array(1:5)*x%d1val2
      q7 = pow2(q1)
      q8 = 2.0_dp*x%d1val1
      q9 = x%d2val1 + log(pow(10.0_dp, q3))
      q10 = 2.0_dp*x%d1val2
      q11 = x%d2val2 + log(pow(10.0_dp, q5))
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1Array(1:5) = q2*x%d1Array(1:5)
      unary%d2val1 = q2*(q1*q3 + x%d2val1)
      unary%d1val1_d1val2 = q2*(q4*x%d1val2 + x%d1val1_d1val2)
      unary%d1val1_d1Array(1:5) = q2*(q4*x%d1Array(1:5) + x%d1val1_d1Array(1:5))
      unary%d2val2 = q2*(q1*q5 + x%d2val2)
      unary%d1val2_d1Array(1:5) = q2*(q1*q6 + x%d1val2_d1Array(1:5))
      unary%d3val1 = q2*(3.0_dp*q4*x%d2val1 + q7*pow3(x%d1val1) + x%d3val1)
      unary%d2val1_d1val2 = log(pow(10.0_dp, q0*(x%d2val1_d1val2 + log(pow(10.0_dp, q8*x%d1val1_d1val2 + q9*x%d1val2)))))
      unary%d2val1_d1Array(1:5) = log(pow(10.0_dp, q0*(x%d2val1_d1Array(1:5) + log(pow(10.0_dp, q8*x%d1val1_d1Array(1:5) + q9*x%d1Array(1:5))))))
      unary%d1val1_d2val2 = log(pow(10.0_dp, q0*(x%d1val1_d2val2 + log(pow(10.0_dp, q10*x%d1val1_d1val2 + q11*x%d1val1)))))
      unary%d1val1_d1val2_d1Array(1:5) = log(pow(10.0_dp, q0*(q6*q7*x%d1val1 + x%d1val1_d1val2_d1Array(1:5) + log(pow(10.0_dp, x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2)))))
      unary%d3val2 = q2*(3.0_dp*q1*x%d1val2*x%d2val2 + q7*pow3(x%d1val2) + x%d3val2)
      unary%d2val2_d1Array(1:5) = log(pow(10.0_dp, q0*(x%d2val2_d1Array(1:5) + log(pow(10.0_dp, q10*x%d1val2_d1Array(1:5) + q11*x%d1Array(1:5))))))
   end function exp10_self
   
   function powm1_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q18(1:5)
      real(dp) :: q17
      real(dp) :: q16(1:5)
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12(1:5)
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7(1:5)
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
      q7(1:5) = x%d1val1_d1Array(1:5)*x%val
      q8 = x%d2val2*x%val
      q9 = -q8
      q10 = pow2(x%d1val2)
      q11 = 2.0_dp*x%d1val2
      q12(1:5) = x%d1val2_d1Array(1:5)*x%val
      q13 = powm1(pow4(x%val))
      q14 = 4.0_dp*x%d1val1
      q15 = 6.0_dp*q4
      q16(1:5) = 2.0_dp*x%d1Array(1:5)
      q17 = 4.0_dp*x%d1val2
      q18(1:5) = 6.0_dp*x%d1Array(1:5)
      unary%val = powm1(x%val)
      unary%d1val1 = -q1*x%d1val1
      unary%d1val2 = -q1*x%d1val2
      unary%d1Array(1:5) = -q1*x%d1Array(1:5)
      unary%d2val1 = q2*(2.0_dp*q4 - q3)
      unary%d1val1_d1val2 = q2*(q5*x%d1val2 - q6)
      unary%d1val1_d1Array(1:5) = q2*(q5*x%d1Array(1:5) - q7)
      unary%d2val2 = q2*(2.0_dp*q10 + q9)
      unary%d1val2_d1Array(1:5) = q2*(q11*x%d1Array(1:5) - q12)
      unary%d3val1 = q13*(-6.0_dp*pow3(x%d1val1) + 6.0_dp*q3*x%d1val1 - q0*x%d3val1)
      unary%d2val1_d1val2 = q13*(-q0*x%d2val1_d1val2 + q11*q3 + q14*q6 - q15*x%d1val2)
      unary%d2val1_d1Array(1:5) = q13*(-q0*x%d2val1_d1Array(1:5) + q14*q7 - q15*x%d1Array(1:5) + q16*q3)
      unary%d1val1_d2val2 = q13*(-q0*x%d1val1_d2val2 + q17*q6 - q5*(3.0_dp*q10 + q9))
      unary%d1val1_d1val2_d1Array(1:5) = q13*(2.0_dp*x%val*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2) - q0*x%d1val1_d1val2_d1Array(1:5) - q18*x%d1val1*x%d1val2)
      unary%d3val2 = q13*(-6.0_dp*pow3(x%d1val2) + 6.0_dp*q8*x%d1val2 - q0*x%d3val2)
      unary%d2val2_d1Array(1:5) = q13*(-q0*x%d2val2_d1Array(1:5) - q10*q18 + q12*q17 + q16*q8)
   end function powm1_self
   
   function log_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7(1:5)
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
      q6 = x%d1val1*x%d1val2
      q7(1:5) = x%d1val1_d1Array(1:5)*x%val
      q8 = x%d2val2*x%val
      q9 = pow2(x%d1val2)
      q10(1:5) = x%d1val2_d1Array(1:5)*x%val
      q11 = powm1(pow3(x%val))
      q12 = 2.0_dp*x%d1val1
      q13 = 2.0_dp*q4
      q14 = 2.0_dp*x%d1val2
      q15 = 2.0_dp*q9
      unary%val = log(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q2*(q3 - q4)
      unary%d1val1_d1val2 = q2*(q5 - q6)
      unary%d1val1_d1Array(1:5) = q2*(q7 - x%d1Array(1:5)*x%d1val1)
      unary%d2val2 = q2*(q8 - q9)
      unary%d1val2_d1Array(1:5) = q2*(q10 - x%d1Array(1:5)*x%d1val2)
      unary%d3val1 = q11*(-3.0_dp*q3*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q1*x%d3val1)
      unary%d2val1_d1val2 = q11*(q1*x%d2val1_d1val2 - q12*q5 + q13*x%d1val2 - q3*x%d1val2)
      unary%d2val1_d1Array(1:5) = q11*(q1*x%d2val1_d1Array(1:5) - q12*q7 + q13*x%d1Array(1:5) - q3*x%d1Array(1:5))
      unary%d1val1_d2val2 = q11*(q1*x%d1val1_d2val2 - q14*q5 + x%d1val1*(q15 - q8))
      unary%d1val1_d1val2_d1Array(1:5) = q11*(2.0_dp*q6*x%d1Array(1:5) + q1*x%d1val1_d1val2_d1Array(1:5) - x%val*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2))
      unary%d3val2 = q11*(-3.0_dp*q8*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q1*x%d3val2)
      unary%d2val2_d1Array(1:5) = q11*(q1*x%d2val2_d1Array(1:5) - q10*q14 + q15*x%d1Array(1:5) - q8*x%d1Array(1:5))
   end function log_self
   
   function log1p_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8(1:5)
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
      q6 = x%d1val1*x%d1val2
      q7 = q0*x%d1val1_d1val2
      q8(1:5) = q0*x%d1val1_d1Array(1:5)
      q9 = pow2(x%d1val2)
      q10 = q0*x%d2val2
      q11(1:5) = q0*x%d1val2_d1Array(1:5)
      q12 = powm1(pow3(q0))
      q13 = 2.0_dp*x%d1val1
      q14 = q4 - q5
      q15 = 2.0_dp*x%d1val2
      q16 = -q10
      unary%val = log1p(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q3*(-q4 + q5)
      unary%d1val1_d1val2 = q3*(-q6 + q7)
      unary%d1val1_d1Array(1:5) = q3*(q8 - x%d1Array(1:5)*x%d1val1)
      unary%d2val2 = q3*(q10 - q9)
      unary%d1val2_d1Array(1:5) = q3*(q11 - x%d1Array(1:5)*x%d1val2)
      unary%d3val1 = q12*(-3.0_dp*q5*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q2*x%d3val1)
      unary%d2val1_d1val2 = q12*(-q13*q7 + q14*x%d1val2 + q2*x%d2val1_d1val2 + q4*x%d1val2)
      unary%d2val1_d1Array(1:5) = q12*(-q13*q8 + q14*x%d1Array(1:5) + q2*x%d2val1_d1Array(1:5) + q4*x%d1Array(1:5))
      unary%d1val1_d2val2 = q12*(-q15*q7 + q2*x%d1val1_d2val2 + x%d1val1*(2.0_dp*q9 + q16))
      unary%d1val1_d1val2_d1Array(1:5) = q12*(2.0_dp*q6*x%d1Array(1:5) - q0*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2) + q2*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q12*(-3.0_dp*q10*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q2*x%d3val2)
      unary%d2val2_d1Array(1:5) = q12*(-q11*q15 + q2*x%d2val2_d1Array(1:5) + q9*x%d1Array(1:5) + x%d1Array(1:5)*(q16 + q9))
   end function log1p_self
   
   function safe_log_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7(1:5)
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
      q6 = x%d1val1*x%d1val2
      q7(1:5) = x%d1val1_d1Array(1:5)*x%val
      q8 = x%d2val2*x%val
      q9 = pow2(x%d1val2)
      q10(1:5) = x%d1val2_d1Array(1:5)*x%val
      q11 = powm1(pow3(x%val))
      q12 = 2.0_dp*x%d1val1
      q13 = 2.0_dp*q4
      q14 = 2.0_dp*x%d1val2
      q15 = 2.0_dp*q9
      unary%val = safe_log(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q2*(q3 - q4)
      unary%d1val1_d1val2 = q2*(q5 - q6)
      unary%d1val1_d1Array(1:5) = q2*(q7 - x%d1Array(1:5)*x%d1val1)
      unary%d2val2 = q2*(q8 - q9)
      unary%d1val2_d1Array(1:5) = q2*(q10 - x%d1Array(1:5)*x%d1val2)
      unary%d3val1 = q11*(-3.0_dp*q3*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q1*x%d3val1)
      unary%d2val1_d1val2 = q11*(q1*x%d2val1_d1val2 - q12*q5 + q13*x%d1val2 - q3*x%d1val2)
      unary%d2val1_d1Array(1:5) = q11*(q1*x%d2val1_d1Array(1:5) - q12*q7 + q13*x%d1Array(1:5) - q3*x%d1Array(1:5))
      unary%d1val1_d2val2 = q11*(q1*x%d1val1_d2val2 - q14*q5 + x%d1val1*(q15 - q8))
      unary%d1val1_d1val2_d1Array(1:5) = q11*(2.0_dp*q6*x%d1Array(1:5) + q1*x%d1val1_d1val2_d1Array(1:5) - x%val*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2))
      unary%d3val2 = q11*(-3.0_dp*q8*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q1*x%d3val2)
      unary%d2val2_d1Array(1:5) = q11*(q1*x%d2val2_d1Array(1:5) - q10*q14 + q15*x%d1Array(1:5) - q8*x%d1Array(1:5))
   end function safe_log_self
   
   function log10_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8(1:5)
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
      q7 = x%d1val1*x%d1val2
      q8(1:5) = x%d1val1_d1Array(1:5)*x%val
      q9 = x%d2val2*x%val
      q10 = pow2(x%d1val2)
      q11(1:5) = x%d1val2_d1Array(1:5)*x%val
      q12 = q0*powm1(pow3(x%val))
      q13 = 2.0_dp*x%d1val1
      q14 = 2.0_dp*q3
      q15 = 2.0_dp*x%d1val2
      q16 = 2.0_dp*q10
      unary%val = q0*log(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q5*(q2 - q3)
      unary%d1val1_d1val2 = q5*(q6 - q7)
      unary%d1val1_d1Array(1:5) = q5*(q8 - x%d1Array(1:5)*x%d1val1)
      unary%d2val2 = q5*(-q10 + q9)
      unary%d1val2_d1Array(1:5) = q5*(q11 - x%d1Array(1:5)*x%d1val2)
      unary%d3val1 = q12*(-3.0_dp*q2*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q4*x%d3val1)
      unary%d2val1_d1val2 = q12*(-q13*q6 + q14*x%d1val2 - q2*x%d1val2 + q4*x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q12*(-q13*q8 + q14*x%d1Array(1:5) - q2*x%d1Array(1:5) + q4*x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q12*(-q15*q6 + q4*x%d1val1_d2val2 + x%d1val1*(q16 - q9))
      unary%d1val1_d1val2_d1Array(1:5) = q12*(2.0_dp*q7*x%d1Array(1:5) + q4*x%d1val1_d1val2_d1Array(1:5) - x%val*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2))
      unary%d3val2 = q12*(-3.0_dp*q9*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q4*x%d3val2)
      unary%d2val2_d1Array(1:5) = q12*(-q11*q15 + q16*x%d1Array(1:5) + q4*x%d2val2_d1Array(1:5) - q9*x%d1Array(1:5))
   end function log10_self
   
   function safe_log10_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8(1:5)
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
      q7 = x%d1val1*x%d1val2
      q8(1:5) = x%d1val1_d1Array(1:5)*x%val
      q9 = x%d2val2*x%val
      q10 = pow2(x%d1val2)
      q11(1:5) = x%d1val2_d1Array(1:5)*x%val
      q12 = q0*powm1(pow3(x%val))
      q13 = 2.0_dp*x%d1val1
      q14 = 2.0_dp*q3
      q15 = 2.0_dp*x%d1val2
      q16 = 2.0_dp*q10
      unary%val = q0*safe_log(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q5*(q2 - q3)
      unary%d1val1_d1val2 = q5*(q6 - q7)
      unary%d1val1_d1Array(1:5) = q5*(q8 - x%d1Array(1:5)*x%d1val1)
      unary%d2val2 = q5*(-q10 + q9)
      unary%d1val2_d1Array(1:5) = q5*(q11 - x%d1Array(1:5)*x%d1val2)
      unary%d3val1 = q12*(-3.0_dp*q2*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q4*x%d3val1)
      unary%d2val1_d1val2 = q12*(-q13*q6 + q14*x%d1val2 - q2*x%d1val2 + q4*x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q12*(-q13*q8 + q14*x%d1Array(1:5) - q2*x%d1Array(1:5) + q4*x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q12*(-q15*q6 + q4*x%d1val1_d2val2 + x%d1val1*(q16 - q9))
      unary%d1val1_d1val2_d1Array(1:5) = q12*(2.0_dp*q7*x%d1Array(1:5) + q4*x%d1val1_d1val2_d1Array(1:5) - x%val*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2))
      unary%d3val2 = q12*(-3.0_dp*q9*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q4*x%d3val2)
      unary%d2val2_d1Array(1:5) = q12*(-q11*q15 + q16*x%d1Array(1:5) + q4*x%d2val2_d1Array(1:5) - q9*x%d1Array(1:5))
   end function safe_log10_self
   
   function log2_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8(1:5)
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
      q7 = x%d1val1*x%d1val2
      q8(1:5) = x%d1val1_d1Array(1:5)*x%val
      q9 = x%d2val2*x%val
      q10 = pow2(x%d1val2)
      q11(1:5) = x%d1val2_d1Array(1:5)*x%val
      q12 = q0*powm1(pow3(x%val))
      q13 = 2.0_dp*x%d1val1
      q14 = 2.0_dp*q3
      q15 = 2.0_dp*x%d1val2
      q16 = 2.0_dp*q10
      unary%val = q0*log(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q5*(q2 - q3)
      unary%d1val1_d1val2 = q5*(q6 - q7)
      unary%d1val1_d1Array(1:5) = q5*(q8 - x%d1Array(1:5)*x%d1val1)
      unary%d2val2 = q5*(-q10 + q9)
      unary%d1val2_d1Array(1:5) = q5*(q11 - x%d1Array(1:5)*x%d1val2)
      unary%d3val1 = q12*(-3.0_dp*q2*x%d1val1 + 2.0_dp*pow3(x%d1val1) + q4*x%d3val1)
      unary%d2val1_d1val2 = q12*(-q13*q6 + q14*x%d1val2 - q2*x%d1val2 + q4*x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q12*(-q13*q8 + q14*x%d1Array(1:5) - q2*x%d1Array(1:5) + q4*x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q12*(-q15*q6 + q4*x%d1val1_d2val2 + x%d1val1*(q16 - q9))
      unary%d1val1_d1val2_d1Array(1:5) = q12*(2.0_dp*q7*x%d1Array(1:5) + q4*x%d1val1_d1val2_d1Array(1:5) - x%val*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2))
      unary%d3val2 = q12*(-3.0_dp*q9*x%d1val2 + 2.0_dp*pow3(x%d1val2) + q4*x%d3val2)
      unary%d2val2_d1Array(1:5) = q12*(-q11*q15 + q16*x%d1Array(1:5) + q4*x%d2val2_d1Array(1:5) - q9*x%d1Array(1:5))
   end function log2_self
   
   function sin_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4(1:5)
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = sin(x%val)
      q1 = cos(x%val)
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4(1:5) = q1*x%d1Array(1:5)
      q5 = pow2(x%d1val1)
      q6 = q0*x%d1val1
      q7 = pow2(x%d1val2)
      q8 = q0*x%d1val2
      q9 = 2.0_dp*q6
      q10(1:5) = q0*x%d1Array(1:5)
      q11 = 2.0_dp*q8
      q12 = q0*x%d2val2
      unary%val = q0
      unary%d1val1 = q2
      unary%d1val2 = q3
      unary%d1Array(1:5) = q4
      unary%d2val1 = -q0*q5 + q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2 - q6*x%d1val2
      unary%d1val1_d1Array(1:5) = q1*x%d1val1_d1Array(1:5) - q6*x%d1Array(1:5)
      unary%d2val2 = -q0*q7 + q1*x%d2val2
      unary%d1val2_d1Array(1:5) = q1*x%d1val2_d1Array(1:5) - q8*x%d1Array(1:5)
      unary%d3val1 = -3.0_dp*q6*x%d2val1 + q1*x%d3val1 - q1*pow3(x%d1val1)
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2 - q3*q5 - q8*x%d2val1 - q9*x%d1val1_d1val2
      unary%d2val1_d1Array(1:5) = q1*x%d2val1_d1Array(1:5) - q10*x%d2val1 - q4*q5 - q9*x%d1val1_d1Array(1:5)
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2 - q11*x%d1val1_d1val2 - x%d1val1*(q1*q7 + q12)
      unary%d1val1_d1val2_d1Array(1:5) = q1*x%d1val1_d1val2_d1Array(1:5) - q10*x%d1val1_d1val2 - q2*x%d1Array(1:5)*x%d1val2 - q6*x%d1val2_d1Array(1:5) - q8*x%d1val1_d1Array(1:5)
      unary%d3val2 = -3.0_dp*q12*x%d1val2 + q1*x%d3val2 - q1*pow3(x%d1val2)
      unary%d2val2_d1Array(1:5) = q1*x%d2val2_d1Array(1:5) - q11*x%d1val2_d1Array(1:5) - q12*x%d1Array(1:5) - q4*q7
   end function sin_self
   
   function cos_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4(1:5)
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = cos(x%val)
      q1 = sin(x%val)
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4(1:5) = q1*x%d1Array(1:5)
      q5 = pow2(x%d1val1)
      q6 = q0*x%d1val1
      q7 = pow2(x%d1val2)
      q8 = q0*x%d1val2
      q9 = 2.0_dp*q6
      q10(1:5) = q0*x%d1Array(1:5)
      q11 = 2.0_dp*q8
      q12 = q0*x%d2val2
      unary%val = q0
      unary%d1val1 = -q2
      unary%d1val2 = -q3
      unary%d1Array(1:5) = -q4
      unary%d2val1 = -q0*q5 - q1*x%d2val1
      unary%d1val1_d1val2 = -q1*x%d1val1_d1val2 - q6*x%d1val2
      unary%d1val1_d1Array(1:5) = -q1*x%d1val1_d1Array(1:5) - q6*x%d1Array(1:5)
      unary%d2val2 = -q0*q7 - q1*x%d2val2
      unary%d1val2_d1Array(1:5) = -q1*x%d1val2_d1Array(1:5) - q8*x%d1Array(1:5)
      unary%d3val1 = -3.0_dp*q6*x%d2val1 - q1*x%d3val1 + q1*pow3(x%d1val1)
      unary%d2val1_d1val2 = -q1*x%d2val1_d1val2 + q3*q5 - q8*x%d2val1 - q9*x%d1val1_d1val2
      unary%d2val1_d1Array(1:5) = -q1*x%d2val1_d1Array(1:5) - q10*x%d2val1 + q4*q5 - q9*x%d1val1_d1Array(1:5)
      unary%d1val1_d2val2 = -q1*x%d1val1_d2val2 - q11*x%d1val1_d1val2 + x%d1val1*(q1*q7 - q12)
      unary%d1val1_d1val2_d1Array(1:5) = -q1*x%d1val1_d1val2_d1Array(1:5) - q10*x%d1val1_d1val2 + q2*x%d1Array(1:5)*x%d1val2 - q6*x%d1val2_d1Array(1:5) - q8*x%d1val1_d1Array(1:5)
      unary%d3val2 = -3.0_dp*q12*x%d1val2 - q1*x%d3val2 + q1*pow3(x%d1val2)
      unary%d2val2_d1Array(1:5) = -q1*x%d2val2_d1Array(1:5) - q11*x%d1val2_d1Array(1:5) - q12*x%d1Array(1:5) + q4*q7
   end function cos_self
   
   function tan_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q24
      real(dp) :: q23(1:5)
      real(dp) :: q22
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19(1:5)
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
      real(dp) :: q4(1:5)
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = tan(x%val)
      q1 = powm1(pow2(cos(x%val)))
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4(1:5) = q1*x%d1Array(1:5)
      q5 = pow2(x%d1val1)
      q6 = 2.0_dp*q0
      q7 = q5*q6 + x%d2val1
      q8 = q6*x%d1val1
      q9 = pow2(x%d1val2)
      q10 = q6*q9 + x%d2val2
      q11 = q6*x%d1val2
      q12 = q0*x%d1val1
      q13 = pow3(x%d1val1)
      q14 = 2.0_dp*q1
      q15 = pow2(q0)
      q16 = 4.0_dp*q15
      q17 = 4.0_dp*q12
      q18 = 2.0_dp*q5
      q19(1:5) = q6*x%d1Array(1:5)
      q20 = 4.0_dp*q0*x%d1val2
      q21 = q0*x%d2val2
      q22 = 2.0_dp*q9
      q23(1:5) = x%d1Array(1:5)*x%d1val2
      q24 = pow3(x%d1val2)
      unary%val = q0
      unary%d1val1 = q2
      unary%d1val2 = q3
      unary%d1Array(1:5) = q4
      unary%d2val1 = q1*q7
      unary%d1val1_d1val2 = q1*(q8*x%d1val2 + x%d1val1_d1val2)
      unary%d1val1_d1Array(1:5) = q1*(q8*x%d1Array(1:5) + x%d1val1_d1Array(1:5))
      unary%d2val2 = q1*q10
      unary%d1val2_d1Array(1:5) = q1*(q11*x%d1Array(1:5) + x%d1val2_d1Array(1:5))
      unary%d3val1 = q1*(6.0_dp*q12*x%d2val1 + q13*q14 + q13*q16 + x%d3val1)
      unary%d2val1_d1val2 = q1*(q11*q7 + q17*x%d1val1_d1val2 + q18*q3 + x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q1*(q17*x%d1val1_d1Array(1:5) + q18*q4 + q19*q7 + x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q1*(2.0_dp*x%d1val1*(q1*q9 + q15*q22 + q21) + q20*x%d1val1_d1val2 + x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = q1*(2.0_dp*q2*q23 + q11*x%d1val1_d1Array(1:5) + q16*q23*x%d1val1 + q19*x%d1val1_d1val2 + q8*x%d1val2_d1Array(1:5) + x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q1*(6.0_dp*q21*x%d1val2 + q14*q24 + q16*q24 + x%d3val2)
      unary%d2val2_d1Array(1:5) = q1*(q10*q19 + q20*x%d1val2_d1Array(1:5) + q22*q4 + x%d2val2_d1Array(1:5))
   end function tan_self
   
   function sinpi_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15(1:5)
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12(1:5)
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
      q6 = q5*x%d1val1
      q7 = pow2(x%d1val2)
      q8 = q5*x%d1val2
      q9 = q2*pow2(pi)
      q10 = 2.0_dp*q6
      q11 = q4*q9
      q12(1:5) = q5*x%d1Array(1:5)
      q13 = 2.0_dp*q8
      q14 = q1*x%d2val2
      q15(1:5) = q9*x%d1Array(1:5)
      q16 = pi*q14
      unary%val = q1
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d1Array(1:5) = q3*x%d1Array(1:5)
      unary%d2val1 = pi*(q2*x%d2val1 - q4*q5)
      unary%d1val1_d1val2 = pi*(q2*x%d1val1_d1val2 - q6*x%d1val2)
      unary%d1val1_d1Array(1:5) = pi*(q2*x%d1val1_d1Array(1:5) - q6*x%d1Array(1:5))
      unary%d2val2 = pi*(q2*x%d2val2 - q5*q7)
      unary%d1val2_d1Array(1:5) = pi*(q2*x%d1val2_d1Array(1:5) - q8*x%d1Array(1:5))
      unary%d3val1 = pi*(-3.0_dp*q6*x%d2val1 + q2*x%d3val1 - q9*pow3(x%d1val1))
      unary%d2val1_d1val2 = pi*(-q10*x%d1val1_d1val2 - q11*x%d1val2 + q2*x%d2val1_d1val2 - q8*x%d2val1)
      unary%d2val1_d1Array(1:5) = pi*(-q10*x%d1val1_d1Array(1:5) - q11*x%d1Array(1:5) - q12*x%d2val1 + q2*x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = -pi*(pi*x%d1val1*(q14 + q3*q7) + q13*x%d1val1_d1val2 - q2*x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = pi*(-q12*x%d1val1_d1val2 - q15*x%d1val1*x%d1val2 + q2*x%d1val1_d1val2_d1Array(1:5) - q6*x%d1val2_d1Array(1:5) - q8*x%d1val1_d1Array(1:5))
      unary%d3val2 = pi*(-3.0_dp*q16*x%d1val2 + q2*x%d3val2 - q9*pow3(x%d1val2))
      unary%d2val2_d1Array(1:5) = pi*(-q13*x%d1val2_d1Array(1:5) - q15*q7 - q16*x%d1Array(1:5) + q2*x%d2val2_d1Array(1:5))
   end function sinpi_self
   
   function cospi_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15(1:5)
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12(1:5)
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
      q6 = q5*x%d1val1
      q7 = pow2(x%d1val2)
      q8 = q5*x%d1val2
      q9 = q2*pow2(pi)
      q10 = 2.0_dp*q6
      q11 = q4*q9
      q12(1:5) = q5*x%d1Array(1:5)
      q13 = 2.0_dp*q8
      q14 = q1*x%d2val2
      q15(1:5) = q9*x%d1Array(1:5)
      q16 = pi*q14
      unary%val = q1
      unary%d1val1 = -q3*x%d1val1
      unary%d1val2 = -q3*x%d1val2
      unary%d1Array(1:5) = -q3*x%d1Array(1:5)
      unary%d2val1 = -pi*(q2*x%d2val1 + q4*q5)
      unary%d1val1_d1val2 = -pi*(q2*x%d1val1_d1val2 + q6*x%d1val2)
      unary%d1val1_d1Array(1:5) = -pi*(q2*x%d1val1_d1Array(1:5) + q6*x%d1Array(1:5))
      unary%d2val2 = -pi*(q2*x%d2val2 + q5*q7)
      unary%d1val2_d1Array(1:5) = -pi*(q2*x%d1val2_d1Array(1:5) + q8*x%d1Array(1:5))
      unary%d3val1 = pi*(-3.0_dp*q6*x%d2val1 - q2*x%d3val1 + q9*pow3(x%d1val1))
      unary%d2val1_d1val2 = pi*(-q10*x%d1val1_d1val2 + q11*x%d1val2 - q2*x%d2val1_d1val2 - q8*x%d2val1)
      unary%d2val1_d1Array(1:5) = pi*(-q10*x%d1val1_d1Array(1:5) + q11*x%d1Array(1:5) - q12*x%d2val1 - q2*x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = -pi*(-pi*x%d1val1*(-q14 + q3*q7) + q13*x%d1val1_d1val2 + q2*x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = pi*(-q12*x%d1val1_d1val2 + q15*x%d1val1*x%d1val2 - q2*x%d1val1_d1val2_d1Array(1:5) - q6*x%d1val2_d1Array(1:5) - q8*x%d1val1_d1Array(1:5))
      unary%d3val2 = pi*(-3.0_dp*q16*x%d1val2 - q2*x%d3val2 + q9*pow3(x%d1val2))
      unary%d2val2_d1Array(1:5) = pi*(-q13*x%d1val2_d1Array(1:5) + q15*q7 - q16*x%d1Array(1:5) - q2*x%d2val2_d1Array(1:5))
   end function cospi_self
   
   function tanpi_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q27
      real(dp) :: q26
      real(dp) :: q25(1:5)
      real(dp) :: q24
      real(dp) :: q23
      real(dp) :: q22(1:5)
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
      q2 = powm1(pow2(cos(q0)))
      q3 = pi*q2
      q4 = pow2(x%d1val1)
      q5 = 2.0_dp*pi
      q6 = q1*q5
      q7 = q4*q6 + x%d2val1
      q8 = q6*x%d1val1
      q9 = pow2(x%d1val2)
      q10 = q6*q9 + x%d2val2
      q11 = q6*x%d1val2
      q12 = pi*q1
      q13 = q12*x%d1val1
      q14 = pow2(pi)
      q15 = q14*pow3(x%d1val1)
      q16 = 2.0_dp*q2
      q17 = pow2(q1)
      q18 = 4.0_dp*q17
      q19 = 4.0_dp*q13
      q20 = q14*q16
      q21 = q20*q4
      q22(1:5) = q6*x%d1Array(1:5)
      q23 = 4.0_dp*q12*x%d1val2
      q24 = q1*x%d2val2
      q25(1:5) = q20*x%d1Array(1:5)
      q26 = x%d1val1*x%d1val2
      q27 = q14*pow3(x%d1val2)
      unary%val = q1
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d1Array(1:5) = q3*x%d1Array(1:5)
      unary%d2val1 = q3*q7
      unary%d1val1_d1val2 = q3*(q8*x%d1val2 + x%d1val1_d1val2)
      unary%d1val1_d1Array(1:5) = q3*(q8*x%d1Array(1:5) + x%d1val1_d1Array(1:5))
      unary%d2val2 = q10*q3
      unary%d1val2_d1Array(1:5) = q3*(q11*x%d1Array(1:5) + x%d1val2_d1Array(1:5))
      unary%d3val1 = q3*(6.0_dp*q13*x%d2val1 + q15*q16 + q15*q18 + x%d3val1)
      unary%d2val1_d1val2 = q3*(q11*q7 + q19*x%d1val1_d1val2 + q21*x%d1val2 + x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q3*(q19*x%d1val1_d1Array(1:5) + q21*x%d1Array(1:5) + q22*q7 + x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q3*(q23*x%d1val1_d1val2 + q5*x%d1val1*(q17*q5*q9 + q24 + q3*q9) + x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = q3*(q11*x%d1val1_d1Array(1:5) + q14*q18*q26*x%d1Array(1:5) + q22*x%d1val1_d1val2 + q25*q26 + q8*x%d1val2_d1Array(1:5) + x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q3*(6.0_dp*pi*q24*x%d1val2 + q16*q27 + q18*q27 + x%d3val2)
      unary%d2val2_d1Array(1:5) = q3*(q10*q22 + q23*x%d1val2_d1Array(1:5) + q25*q9 + x%d2val2_d1Array(1:5))
   end function tanpi_self
   
   function sinh_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4(1:5)
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = sinh(x%val)
      q1 = cosh(x%val)
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4(1:5) = q1*x%d1Array(1:5)
      q5 = pow2(x%d1val1)
      q6 = q0*x%d1val1
      q7 = pow2(x%d1val2)
      q8 = q0*x%d1val2
      q9 = 2.0_dp*q6
      q10(1:5) = q0*x%d1Array(1:5)
      q11 = 2.0_dp*q8
      q12 = q0*x%d2val2
      unary%val = q0
      unary%d1val1 = q2
      unary%d1val2 = q3
      unary%d1Array(1:5) = q4
      unary%d2val1 = q0*q5 + q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2 + q6*x%d1val2
      unary%d1val1_d1Array(1:5) = q1*x%d1val1_d1Array(1:5) + q6*x%d1Array(1:5)
      unary%d2val2 = q0*q7 + q1*x%d2val2
      unary%d1val2_d1Array(1:5) = q1*x%d1val2_d1Array(1:5) + q8*x%d1Array(1:5)
      unary%d3val1 = 3.0_dp*q6*x%d2val1 + q1*x%d3val1 + q1*pow3(x%d1val1)
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2 + q3*q5 + q8*x%d2val1 + q9*x%d1val1_d1val2
      unary%d2val1_d1Array(1:5) = q1*x%d2val1_d1Array(1:5) + q10*x%d2val1 + q4*q5 + q9*x%d1val1_d1Array(1:5)
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2 + q11*x%d1val1_d1val2 + x%d1val1*(q1*q7 + q12)
      unary%d1val1_d1val2_d1Array(1:5) = q1*x%d1val1_d1val2_d1Array(1:5) + q10*x%d1val1_d1val2 + q2*x%d1Array(1:5)*x%d1val2 + q6*x%d1val2_d1Array(1:5) + q8*x%d1val1_d1Array(1:5)
      unary%d3val2 = 3.0_dp*q12*x%d1val2 + q1*x%d3val2 + q1*pow3(x%d1val2)
      unary%d2val2_d1Array(1:5) = q1*x%d2val2_d1Array(1:5) + q11*x%d1val2_d1Array(1:5) + q12*x%d1Array(1:5) + q4*q7
   end function sinh_self
   
   function cosh_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4(1:5)
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = cosh(x%val)
      q1 = sinh(x%val)
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4(1:5) = q1*x%d1Array(1:5)
      q5 = pow2(x%d1val1)
      q6 = q0*x%d1val1
      q7 = pow2(x%d1val2)
      q8 = q0*x%d1val2
      q9 = 2.0_dp*q6
      q10(1:5) = q0*x%d1Array(1:5)
      q11 = 2.0_dp*q8
      q12 = q0*x%d2val2
      unary%val = q0
      unary%d1val1 = q2
      unary%d1val2 = q3
      unary%d1Array(1:5) = q4
      unary%d2val1 = q0*q5 + q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2 + q6*x%d1val2
      unary%d1val1_d1Array(1:5) = q1*x%d1val1_d1Array(1:5) + q6*x%d1Array(1:5)
      unary%d2val2 = q0*q7 + q1*x%d2val2
      unary%d1val2_d1Array(1:5) = q1*x%d1val2_d1Array(1:5) + q8*x%d1Array(1:5)
      unary%d3val1 = 3.0_dp*q6*x%d2val1 + q1*x%d3val1 + q1*pow3(x%d1val1)
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2 + q3*q5 + q8*x%d2val1 + q9*x%d1val1_d1val2
      unary%d2val1_d1Array(1:5) = q1*x%d2val1_d1Array(1:5) + q10*x%d2val1 + q4*q5 + q9*x%d1val1_d1Array(1:5)
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2 + q11*x%d1val1_d1val2 + x%d1val1*(q1*q7 + q12)
      unary%d1val1_d1val2_d1Array(1:5) = q1*x%d1val1_d1val2_d1Array(1:5) + q10*x%d1val1_d1val2 + q2*x%d1Array(1:5)*x%d1val2 + q6*x%d1val2_d1Array(1:5) + q8*x%d1val1_d1Array(1:5)
      unary%d3val2 = 3.0_dp*q12*x%d1val2 + q1*x%d3val2 + q1*pow3(x%d1val2)
      unary%d2val2_d1Array(1:5) = q1*x%d2val2_d1Array(1:5) + q11*x%d1val2_d1Array(1:5) + q12*x%d1Array(1:5) + q4*q7
   end function cosh_self
   
   function tanh_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q22
      real(dp) :: q21(1:5)
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17(1:5)
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
      q0 = tanh(x%val)
      q1 = powm1(pow2(cosh(x%val)))
      q2 = pow2(q0)
      q3 = q2 - 1
      q4 = pow2(x%d1val1)
      q5 = 2.0_dp*q0
      q6 = q4*q5 - x%d2val1
      q7 = q5*x%d1val1
      q8 = pow2(x%d1val2)
      q9 = q5*q8 - x%d2val2
      q10 = q5*x%d1val2
      q11 = q0*x%d1val1
      q12 = pow3(x%d1val1)
      q13 = 4.0_dp*q2
      q14 = 2.0_dp*q3
      q15 = 4.0_dp*q11
      q16 = q14*q4
      q17(1:5) = q5*x%d1Array(1:5)
      q18 = 4.0_dp*q0*x%d1val2
      q19 = q0*x%d2val2
      q20 = q3*q8
      q21(1:5) = x%d1Array(1:5)*x%d1val1*x%d1val2
      q22 = pow3(x%d1val2)
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q3*q6
      unary%d1val1_d1val2 = q3*(q7*x%d1val2 - x%d1val1_d1val2)
      unary%d1val1_d1Array(1:5) = q3*(q7*x%d1Array(1:5) - x%d1val1_d1Array(1:5))
      unary%d2val2 = q3*q9
      unary%d1val2_d1Array(1:5) = q3*(q10*x%d1Array(1:5) - x%d1val2_d1Array(1:5))
      unary%d3val1 = -q3*(-6.0_dp*q11*x%d2val1 + q12*q13 + q12*q14 + x%d3val1)
      unary%d2val1_d1val2 = -q3*(q10*q6 - q15*x%d1val1_d1val2 + q16*x%d1val2 + x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = -q3*(-q15*x%d1val1_d1Array(1:5) + q16*x%d1Array(1:5) + q17*q6 + x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = -q3*(2.0_dp*x%d1val1*(2.0_dp*q2*q8 - q19 + q20) - q18*x%d1val1_d1val2 + x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = -q3*(-q10*x%d1val1_d1Array(1:5) + q13*q21 + q14*q21 - q17*x%d1val1_d1val2 - q7*x%d1val2_d1Array(1:5) + x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = -q3*(-6.0_dp*q19*x%d1val2 + q13*q22 + q14*q22 + x%d3val2)
      unary%d2val2_d1Array(1:5) = -q3*(2.0_dp*q20*x%d1Array(1:5) + q17*q9 - q18*x%d1val2_d1Array(1:5) + x%d2val2_d1Array(1:5))
   end function tanh_self
   
   function asin_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q26(1:5)
      real(dp) :: q25
      real(dp) :: q24
      real(dp) :: q23(1:5)
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
      real(dp) :: q11(1:5)
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8(1:5)
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
      q8(1:5) = x%d1Array(1:5)*x%d1val1
      q9 = pow2(x%d1val2)
      q10 = -q5*x%d2val2 + q9*x%val
      q11(1:5) = x%d1Array(1:5)*x%d1val2
      q12 = powm1(pow5(sqrt(q1)))
      q13 = 3.0_dp*q0
      q14 = pow2(q5)
      q15 = q5*x%d1val1
      q16 = powm1(pow7(sqrt(q1)))
      q17 = x%d1val2*x%val
      q18 = q1*q6
      q19 = 2.0_dp*q0
      q20 = q19*q4
      q21 = 2.0_dp*x%val
      q22 = q21*x%d1val1_d1val2
      q23(1:5) = x%d1Array(1:5)*x%val
      q24 = q5*x%d1val2
      q25 = x%d2val2*x%val
      q26(1:5) = q7*x%d1Array(1:5)
      unary%val = asin(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1Array(1:5) = q2*x%d1Array(1:5)
      unary%d2val1 = q3*q6
      unary%d1val1_d1val2 = q3*(q1*x%d1val1_d1val2 + q7*x%val)
      unary%d1val1_d1Array(1:5) = q3*(q1*x%d1val1_d1Array(1:5) + q8*x%val)
      unary%d2val2 = q10*q3
      unary%d1val2_d1Array(1:5) = q3*(q1*x%d1val2_d1Array(1:5) + q11*x%val)
      unary%d3val1 = q12*(q13*pow3(x%d1val1) + q14*x%d3val1 - q15*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = q16*(q17*q18 + q5*(-q14*x%d2val1_d1val2 + q15*(q22 + q7) - q20*x%d1val2))
      unary%d2val1_d1Array(1:5) = q16*(q18*q23 + q5*(-q14*x%d2val1_d1Array(1:5) + q15*(q21*x%d1val1_d1Array(1:5) + q8) - q20*x%d1Array(1:5)))
      unary%d1val1_d2val2 = q12*(q14*x%d1val1_d2val2 - q22*q24 - x%d1val1*(-q13*q9 + q5*(q25 + q9)))
      unary%d1val1_d1val2_d1Array(1:5) = q12*(q1*(q17*x%d1val1_d1Array(1:5) + q23*x%d1val1_d1val2 + q26 + x%d1val1*x%d1val2_d1Array(1:5)*x%val) + q13*q26 + x%d1val1_d1val2_d1Array(1:5)*pow2(q1))
      unary%d3val2 = q12*(q13*pow3(x%d1val2) + q14*x%d3val2 - q24*(3.0_dp*q25 + q9))
      unary%d2val2_d1Array(1:5) = q16*(q1*q10*q23 + q5*(-q14*x%d2val2_d1Array(1:5) - q19*q9*x%d1Array(1:5) + q24*(q11 + q21*x%d1val2_d1Array(1:5))))
   end function asin_self
   
   function acos_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q29(1:5)
      real(dp) :: q28
      real(dp) :: q27
      real(dp) :: q26(1:5)
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
      real(dp) :: q14(1:5)
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
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
      q10(1:5) = x%d1Array(1:5)*x%d1val1
      q11 = pow2(x%d1val2)
      q12 = q11*x%val
      q13 = q6*x%d2val2
      q14(1:5) = x%d1Array(1:5)*x%d1val2
      q15 = powm1(pow7(sqrt(q1)))
      q16 = pow3(q6)
      q17 = 3.0_dp*q0
      q18 = q1*q17
      q19 = pow2(q6)
      q20 = x%d1val2*x%val
      q21 = q1*(q5 - q7)
      q22 = 2.0_dp*q0
      q23 = q22*q4
      q24 = 2.0_dp*x%val
      q25 = q6*x%d1val1
      q26(1:5) = x%d1Array(1:5)*x%val
      q27 = powm1(pow5(sqrt(q1)))
      q28 = x%d2val2*x%val
      q29(1:5) = q8*x%d1Array(1:5)
      unary%val = acos(x%val)
      unary%d1val1 = -q2*x%d1val1
      unary%d1val2 = -q2*x%d1val2
      unary%d1Array(1:5) = -q2*x%d1Array(1:5)
      unary%d2val1 = q3*(-q5 + q7)
      unary%d1val1_d1val2 = -q3*(q8*x%val - q9)
      unary%d1val1_d1Array(1:5) = -q3*(q10*x%val - q6*x%d1val1_d1Array(1:5))
      unary%d2val2 = q3*(-q12 + q13)
      unary%d1val2_d1Array(1:5) = -q3*(q14*x%val - q6*x%d1val2_d1Array(1:5))
      unary%d3val1 = q15*(q16*x%d3val1 - q18*pow3(x%d1val1) - q19*x%d1val1*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = q15*(-q20*q21 + q6*(q19*x%d2val1_d1val2 + q23*x%d1val2 - q25*(q24*x%d1val1_d1val2 + q8)))
      unary%d2val1_d1Array(1:5) = q15*(-q21*q26 + q6*(q19*x%d2val1_d1Array(1:5) + q23*x%d1Array(1:5) - q25*(q10 + q24*x%d1val1_d1Array(1:5))))
      unary%d1val1_d2val2 = q27*(-q19*x%d1val1_d2val2 + q24*q9*x%d1val2 + x%d1val1*(-q11*q17 + q6*(q11 + q28)))
      unary%d1val1_d1val2_d1Array(1:5) = -q27*(q17*q29 + q19*x%d1val1_d1val2_d1Array(1:5) - q6*(q20*x%d1val1_d1Array(1:5) + q26*x%d1val1_d1val2 + q29 + x%d1val1*x%d1val2_d1Array(1:5)*x%val))
      unary%d3val2 = q15*(q16*x%d3val2 - q18*pow3(x%d1val2) - q19*x%d1val2*(3.0_dp*q28 + q11))
      unary%d2val2_d1Array(1:5) = q15*(-q1*q26*(q12 - q13) + q6*(q11*q22*x%d1Array(1:5) + q19*x%d2val2_d1Array(1:5) - q6*x%d1val2*(q14 + q24*x%d1val2_d1Array(1:5))))
   end function acos_self
   
   function atan_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q28
      real(dp) :: q27
      real(dp) :: q26(1:5)
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
      real(dp) :: q15(1:5)
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
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
      q11(1:5) = x%d1Array(1:5)*x%d1val1
      q12 = pow2(x%d1val2)
      q13 = q12*q6
      q14 = q1*x%d2val2
      q15(1:5) = x%d1Array(1:5)*x%d1val2
      q16 = powm1(pow3(q1))
      q17 = 8.0_dp*q0
      q18 = 2.0_dp*x%d1val1
      q19 = q1*q18
      q20 = 4.0_dp*q0
      q21 = q20*q5
      q22 = q6*(q7 - q8)
      q23 = x%d1val2*x%val
      q24 = q12*q20
      q25 = x%d2val2*x%val
      q26(1:5) = q9*x%d1Array(1:5)
      q27 = 2.0_dp*q1
      q28 = q27*x%d1val2
      unary%val = atan(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1Array(1:5) = q2*x%d1Array(1:5)
      unary%d2val1 = q4*(-q7 + q8)
      unary%d1val1_d1val2 = q4*(q10 - q6*q9)
      unary%d1val1_d1Array(1:5) = q4*(q1*x%d1val1_d1Array(1:5) - q11*q6)
      unary%d2val2 = q4*(-q13 + q14)
      unary%d1val2_d1Array(1:5) = q4*(q1*x%d1val2_d1Array(1:5) - q15*q6)
      unary%d3val1 = q16*(q17*pow3(x%d1val1) - q19*(3.0_dp*x%d2val1*x%val + q5) + q3*x%d3val1)
      unary%d2val1_d1val2 = q16*(-q19*(q6*x%d1val1_d1val2 + q9) + q21*x%d1val2 + q22*x%d1val2 + q3*x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q16*(-q19*(q11 + q6*x%d1val1_d1Array(1:5)) + q21*x%d1Array(1:5) + q22*x%d1Array(1:5) + q3*x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q16*(-4.0_dp*q10*q23 - q18*(q1*(q12 + q25) - q24) + q3*x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = q16*(q17*q26 - q27*(q23*x%d1val1_d1Array(1:5) + q26 + x%d1Array(1:5)*x%d1val1_d1val2*x%val + x%d1val1*x%d1val2_d1Array(1:5)*x%val) + q3*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q16*(q17*pow3(x%d1val2) - q28*(3.0_dp*q25 + q12) + q3*x%d3val2)
      unary%d2val2_d1Array(1:5) = q16*(q24*x%d1Array(1:5) - q28*(q15 + q6*x%d1val2_d1Array(1:5)) + q3*x%d2val2_d1Array(1:5) + q6*x%d1Array(1:5)*(q13 - q14))
   end function atan_self
   
   function asinpi_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q27(1:5)
      real(dp) :: q26(1:5)
      real(dp) :: q25(1:5)
      real(dp) :: q24
      real(dp) :: q23(1:5)
      real(dp) :: q22(1:5)
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
      q0 = powm1(pi)
      q1 = pow2(x%val)
      q2 = 1 - q1
      q3 = q0*powm1(sqrt(q2))
      q4 = pow2(x%d1val1)
      q5 = q1 - 1
      q6 = q0*powm1(pow3(sqrt(q2)))
      q7 = x%d1val1*x%val
      q8 = pow2(x%d1val2)
      q9 = x%d1val2*x%val
      q10 = 3.0_dp*q1
      q11 = pow2(q5)
      q12 = x%d2val1*x%val
      q13 = q0*powm1(pow5(sqrt(q2)))
      q14 = q4*x%d1val2
      q15 = pow4(x%val)
      q16 = 2.0_dp*x%d1val1_d1val2
      q17 = 2.0_dp*q1
      q18 = pow3(x%val)
      q19 = q18*x%d1val1
      q20 = q18*x%d2val1
      q21 = q3*powm1(q15 - q17 + 1)
      q22(1:5) = q4*x%d1Array(1:5)
      q23(1:5) = 2.0_dp*x%d1val1_d1Array(1:5)
      q24 = x%d2val2*x%val
      q25(1:5) = x%d1Array(1:5)*x%d1val1*x%d1val2
      q26(1:5) = q8*x%d1Array(1:5)
      q27(1:5) = 2.0_dp*x%d1val2_d1Array(1:5)
      unary%val = q0*asin(x%val)
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d1Array(1:5) = q3*x%d1Array(1:5)
      unary%d2val1 = q6*(q4*x%val - q5*x%d2val1)
      unary%d1val1_d1val2 = q6*(q2*x%d1val1_d1val2 + q7*x%d1val2)
      unary%d1val1_d1Array(1:5) = q6*(q2*x%d1val1_d1Array(1:5) + q7*x%d1Array(1:5))
      unary%d2val2 = q6*(-q5*x%d2val2 + q8*x%val)
      unary%d1val2_d1Array(1:5) = q6*(q2*x%d1val2_d1Array(1:5) + q9*x%d1Array(1:5))
      unary%d3val1 = q13*(q10*pow3(x%d1val1) + q11*x%d3val1 - q5*x%d1val1*(3.0_dp*q12 + q4))
      unary%d2val1_d1val2 = q21*(q14*q17 + q14 + q15*x%d2val1_d1val2 - q16*q19 + q16*q7 - q17*x%d2val1_d1val2 - q20*x%d1val2 + q9*x%d2val1 + x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q21*(q12*x%d1Array(1:5) + q15*x%d2val1_d1Array(1:5) + q17*q22 - q17*x%d2val1_d1Array(1:5) - q19*q23 - q20*x%d1Array(1:5) + q22 + q23*q7 + x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q13*(q11*x%d1val1_d2val2 - q16*q5*q9 - x%d1val1*(-q10*q8 + q5*(q24 + q8)))
      unary%d1val1_d1val2_d1Array(1:5) = q13*(q10*q25 + q2*(q25 + q7*x%d1val2_d1Array(1:5) + q9*x%d1val1_d1Array(1:5) + x%d1Array(1:5)*x%d1val1_d1val2*x%val) + x%d1val1_d1val2_d1Array(1:5)*pow2(q2))
      unary%d3val2 = q13*(q10*pow3(x%d1val2) + q11*x%d3val2 - q5*x%d1val2*(3.0_dp*q24 + q8))
      unary%d2val2_d1Array(1:5) = q21*(q15*x%d2val2_d1Array(1:5) + q17*q26 - q17*x%d2val2_d1Array(1:5) - q18*q27*x%d1val2 - q18*x%d1Array(1:5)*x%d2val2 + q24*x%d1Array(1:5) + q26 + q27*q9 + x%d2val2_d1Array(1:5))
   end function asinpi_self
   
   function acospi_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q30(1:5)
      real(dp) :: q29
      real(dp) :: q28
      real(dp) :: q27(1:5)
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
      real(dp) :: q15(1:5)
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
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
      q5 = q4*x%val
      q6 = q1 - 1
      q7 = q6*x%d2val1
      q8 = q0*powm1(pow3(sqrt(q2)))
      q9 = x%d1val1*x%d1val2
      q10 = q6*x%d1val1_d1val2
      q11(1:5) = x%d1Array(1:5)*x%d1val1
      q12 = pow2(x%d1val2)
      q13 = q12*x%val
      q14 = q6*x%d2val2
      q15(1:5) = x%d1Array(1:5)*x%d1val2
      q16 = pow3(q6)
      q17 = 3.0_dp*q1
      q18 = q1*(3.0_dp - q17)
      q19 = pow2(q6)
      q20 = q0*powm1(pow7(sqrt(q2)))
      q21 = x%d1val2*x%val
      q22 = q6*(q5 - q7)
      q23 = 2.0_dp*q1
      q24 = q23*q4
      q25 = 2.0_dp*x%val
      q26 = q6*x%d1val1
      q27(1:5) = x%d1Array(1:5)*x%val
      q28 = x%d2val2*x%val
      q29 = q0*powm1(pow5(sqrt(q2)))
      q30(1:5) = q9*x%d1Array(1:5)
      unary%val = q0*acos(x%val)
      unary%d1val1 = -q3*x%d1val1
      unary%d1val2 = -q3*x%d1val2
      unary%d1Array(1:5) = -q3*x%d1Array(1:5)
      unary%d2val1 = q8*(-q5 + q7)
      unary%d1val1_d1val2 = -q8*(-q10 + q9*x%val)
      unary%d1val1_d1Array(1:5) = -q8*(q11*x%val - q6*x%d1val1_d1Array(1:5))
      unary%d2val2 = q8*(-q13 + q14)
      unary%d1val2_d1Array(1:5) = -q8*(q15*x%val - q6*x%d1val2_d1Array(1:5))
      unary%d3val1 = q20*(q16*x%d3val1 - q18*pow3(x%d1val1) - q19*x%d1val1*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = q20*(q2*(-q19*x%d2val1_d1val2 - q24*x%d1val2 + q26*(q25*x%d1val1_d1val2 + q9)) + q21*q22)
      unary%d2val1_d1Array(1:5) = q20*(q2*(-q19*x%d2val1_d1Array(1:5) - q24*x%d1Array(1:5) + q26*(q11 + q25*x%d1val1_d1Array(1:5))) + q22*q27)
      unary%d1val1_d2val2 = q29*(q10*q25*x%d1val2 - q19*x%d1val1_d2val2 + x%d1val1*(-q12*q17 + q6*(q12 + q28)))
      unary%d1val1_d1val2_d1Array(1:5) = -q29*(q17*q30 + q19*x%d1val1_d1val2_d1Array(1:5) - q6*(q21*x%d1val1_d1Array(1:5) + q27*x%d1val1_d1val2 + q30 + x%d1val1*x%d1val2_d1Array(1:5)*x%val))
      unary%d3val2 = q20*(q16*x%d3val2 - q18*pow3(x%d1val2) - q19*x%d1val2*(3.0_dp*q28 + q12))
      unary%d2val2_d1Array(1:5) = q20*(q2*(-q12*q23*x%d1Array(1:5) - q19*x%d2val2_d1Array(1:5) + q6*x%d1val2*(q15 + q25*x%d1val2_d1Array(1:5))) + q27*q6*(q13 - q14))
   end function acospi_self
   
   function atanpi_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q29
      real(dp) :: q28
      real(dp) :: q27(1:5)
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
      real(dp) :: q16(1:5)
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12(1:5)
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
      q2 = q1 + 1
      q3 = q0*powm1(q2)
      q4 = pow2(x%d1val1)
      q5 = 2.0_dp*x%val
      q6 = q4*q5
      q7 = q2*x%d2val1
      q8 = pow2(q2)
      q9 = q0*powm1(q8)
      q10 = x%d1val1*x%d1val2
      q11 = q2*x%d1val1_d1val2
      q12(1:5) = x%d1Array(1:5)*x%d1val1
      q13 = pow2(x%d1val2)
      q14 = q13*q5
      q15 = q2*x%d2val2
      q16(1:5) = x%d1Array(1:5)*x%d1val2
      q17 = 8.0_dp*q1
      q18 = 2.0_dp*x%d1val1
      q19 = q18*q2
      q20 = q0*powm1(pow3(q2))
      q21 = 4.0_dp*q1
      q22 = q21*q4
      q23 = q5*(q6 - q7)
      q24 = x%d1val2*x%val
      q25 = q13*q21
      q26 = x%d2val2*x%val
      q27(1:5) = q10*x%d1Array(1:5)
      q28 = 2.0_dp*q2
      q29 = q28*x%d1val2
      unary%val = q0*atan(x%val)
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d1Array(1:5) = q3*x%d1Array(1:5)
      unary%d2val1 = q9*(-q6 + q7)
      unary%d1val1_d1val2 = q9*(-q10*q5 + q11)
      unary%d1val1_d1Array(1:5) = q9*(-q12*q5 + q2*x%d1val1_d1Array(1:5))
      unary%d2val2 = q9*(-q14 + q15)
      unary%d1val2_d1Array(1:5) = q9*(-q16*q5 + q2*x%d1val2_d1Array(1:5))
      unary%d3val1 = q20*(q17*pow3(x%d1val1) - q19*(3.0_dp*x%d2val1*x%val + q4) + q8*x%d3val1)
      unary%d2val1_d1val2 = q20*(-q19*(q10 + q5*x%d1val1_d1val2) + q22*x%d1val2 + q23*x%d1val2 + q8*x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q20*(-q19*(q12 + q5*x%d1val1_d1Array(1:5)) + q22*x%d1Array(1:5) + q23*x%d1Array(1:5) + q8*x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q20*(-4.0_dp*q11*q24 - q18*(q2*(q13 + q26) - q25) + q8*x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = q20*(q17*q27 - q28*(q24*x%d1val1_d1Array(1:5) + q27 + x%d1Array(1:5)*x%d1val1_d1val2*x%val + x%d1val1*x%d1val2_d1Array(1:5)*x%val) + q8*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q20*(q17*pow3(x%d1val2) - q29*(3.0_dp*q26 + q13) + q8*x%d3val2)
      unary%d2val2_d1Array(1:5) = q20*(q25*x%d1Array(1:5) - q29*(q16 + q5*x%d1val2_d1Array(1:5)) + q5*x%d1Array(1:5)*(q14 - q15) + q8*x%d2val2_d1Array(1:5))
   end function atanpi_self
   
   function asinh_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q25
      real(dp) :: q24(1:5)
      real(dp) :: q23(1:5)
      real(dp) :: q22
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13(1:5)
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9(1:5)
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
      q9(1:5) = x%d1Array(1:5)*x%d1val1
      q10 = pow2(x%d1val2)
      q11 = q10*x%val
      q12 = q1*x%d2val2
      q13(1:5) = x%d1Array(1:5)*x%d1val2
      q14 = powm1(pow5(sqrt(q1)))
      q15 = 3.0_dp*q0
      q16 = pow2(q1)
      q17 = q1*x%d1val1
      q18 = 2.0_dp*q0
      q19 = q18*q4
      q20 = 2.0_dp*x%val
      q21 = x%val*(q5 - q6)
      q22 = x%d2val2*x%val
      q23(1:5) = q7*x%d1Array(1:5)
      q24(1:5) = x%d1Array(1:5)*x%val
      q25 = q1*x%d1val2
      unary%val = asinh(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1Array(1:5) = q2*x%d1Array(1:5)
      unary%d2val1 = q3*(-q5 + q6)
      unary%d1val1_d1val2 = q3*(-q7*x%val + q8)
      unary%d1val1_d1Array(1:5) = q3*(q1*x%d1val1_d1Array(1:5) - q9*x%val)
      unary%d2val2 = q3*(-q11 + q12)
      unary%d1val2_d1Array(1:5) = q3*(q1*x%d1val2_d1Array(1:5) - q13*x%val)
      unary%d3val1 = q14*(q15*pow3(x%d1val1) + q16*x%d3val1 - q17*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = q14*(q16*x%d2val1_d1val2 - q17*(q20*x%d1val1_d1val2 + q7) + q19*x%d1val2 + q21*x%d1val2)
      unary%d2val1_d1Array(1:5) = q14*(q16*x%d2val1_d1Array(1:5) - q17*(q20*x%d1val1_d1Array(1:5) + q9) + q19*x%d1Array(1:5) + q21*x%d1Array(1:5))
      unary%d1val1_d2val2 = q14*(q16*x%d1val1_d2val2 - q20*q8*x%d1val2 - x%d1val1*(q1*(q10 + q22) - q10*q15))
      unary%d1val1_d1val2_d1Array(1:5) = q14*(-q1*(q23 + q24*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5)*x%val + x%d1val1_d1Array(1:5)*x%d1val2*x%val) + q15*q23 + q16*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q14*(q15*pow3(x%d1val2) + q16*x%d3val2 - q25*(3.0_dp*q22 + q10))
      unary%d2val2_d1Array(1:5) = q14*(q10*q18*x%d1Array(1:5) + q16*x%d2val2_d1Array(1:5) + q24*(q11 - q12) - q25*(q13 + q20*x%d1val2_d1Array(1:5)))
   end function asinh_self
   
   function acosh_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q25
      real(dp) :: q24(1:5)
      real(dp) :: q23(1:5)
      real(dp) :: q22
      real(dp) :: q21
      real(dp) :: q20
      real(dp) :: q19
      real(dp) :: q18
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13(1:5)
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9(1:5)
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
      q9(1:5) = x%d1Array(1:5)*x%d1val1
      q10 = pow2(x%d1val2)
      q11 = q10*x%val
      q12 = q1*x%d2val2
      q13(1:5) = x%d1Array(1:5)*x%d1val2
      q14 = powm1(pow5(sqrt(q1)))
      q15 = 3.0_dp*q0
      q16 = pow2(q1)
      q17 = q1*x%d1val1
      q18 = 2.0_dp*q0
      q19 = q18*q4
      q20 = 2.0_dp*x%val
      q21 = x%val*(q5 - q6)
      q22 = x%d2val2*x%val
      q23(1:5) = q7*x%d1Array(1:5)
      q24(1:5) = x%d1Array(1:5)*x%val
      q25 = q1*x%d1val2
      unary%val = acosh(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1Array(1:5) = q2*x%d1Array(1:5)
      unary%d2val1 = q3*(-q5 + q6)
      unary%d1val1_d1val2 = q3*(-q7*x%val + q8)
      unary%d1val1_d1Array(1:5) = q3*(q1*x%d1val1_d1Array(1:5) - q9*x%val)
      unary%d2val2 = q3*(-q11 + q12)
      unary%d1val2_d1Array(1:5) = q3*(q1*x%d1val2_d1Array(1:5) - q13*x%val)
      unary%d3val1 = q14*(q15*pow3(x%d1val1) + q16*x%d3val1 - q17*(3.0_dp*x%d2val1*x%val + q4))
      unary%d2val1_d1val2 = q14*(q16*x%d2val1_d1val2 - q17*(q20*x%d1val1_d1val2 + q7) + q19*x%d1val2 + q21*x%d1val2)
      unary%d2val1_d1Array(1:5) = q14*(q16*x%d2val1_d1Array(1:5) - q17*(q20*x%d1val1_d1Array(1:5) + q9) + q19*x%d1Array(1:5) + q21*x%d1Array(1:5))
      unary%d1val1_d2val2 = q14*(q16*x%d1val1_d2val2 - q20*q8*x%d1val2 - x%d1val1*(q1*(q10 + q22) - q10*q15))
      unary%d1val1_d1val2_d1Array(1:5) = q14*(-q1*(q23 + q24*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5)*x%val + x%d1val1_d1Array(1:5)*x%d1val2*x%val) + q15*q23 + q16*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q14*(q15*pow3(x%d1val2) + q16*x%d3val2 - q25*(3.0_dp*q22 + q10))
      unary%d2val2_d1Array(1:5) = q14*(q10*q18*x%d1Array(1:5) + q16*x%d2val2_d1Array(1:5) + q24*(q11 - q12) - q25*(q13 + q20*x%d1val2_d1Array(1:5)))
   end function acosh_self
   
   function atanh_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q28
      real(dp) :: q27(1:5)
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
      real(dp) :: q13(1:5)
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
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
      q10(1:5) = x%d1Array(1:5)*x%d1val1
      q11 = pow2(x%d1val2)
      q12 = -q1*x%d2val2 + q11*q6
      q13(1:5) = x%d1Array(1:5)*x%d1val2
      q14 = powm1(pow4(q1))
      q15 = pow3(q1)
      q16 = 8.0_dp*q0
      q17 = q16*(1 - q0)
      q18 = 2.0_dp*x%d1val1
      q19 = powm1(q15)
      q20 = 4.0_dp*q0
      q21 = q20*q5
      q22 = q1*q18
      q23 = q6*q7
      q24 = x%d1val2*x%val
      q25 = q11*q20
      q26 = x%d2val2*x%val
      q27(1:5) = q8*x%d1Array(1:5)
      q28 = 2.0_dp*q1
      unary%val = atanh(x%val)
      unary%d1val1 = -q2*x%d1val1
      unary%d1val2 = -q2*x%d1val2
      unary%d1Array(1:5) = -q2*x%d1Array(1:5)
      unary%d2val1 = q4*q7
      unary%d1val1_d1val2 = q4*(q6*q8 - q9)
      unary%d1val1_d1Array(1:5) = q4*(-q1*x%d1val1_d1Array(1:5) + q10*q6)
      unary%d2val2 = q12*q4
      unary%d1val2_d1Array(1:5) = q4*(-q1*x%d1val2_d1Array(1:5) + q13*q6)
      unary%d3val1 = q14*(-q15*x%d3val1 + q17*pow3(x%d1val1) + q18*q3*(3.0_dp*x%d2val1*x%val + q5))
      unary%d2val1_d1val2 = q19*(-q21*x%d1val2 + q22*(q6*x%d1val1_d1val2 + q8) - q23*x%d1val2 - q3*x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q19*(-q21*x%d1Array(1:5) + q22*(q10 + q6*x%d1val1_d1Array(1:5)) - q23*x%d1Array(1:5) - q3*x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q19*(4.0_dp*q24*q9 + q18*(q1*(q11 + q26) - q25) - q3*x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = q19*(-q16*q27 + q28*(q24*x%d1val1_d1Array(1:5) + q27 + x%d1Array(1:5)*x%d1val1_d1val2*x%val + x%d1val1*x%d1val2_d1Array(1:5)*x%val) - q3*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q14*(2.0_dp*q3*x%d1val2*(3.0_dp*q26 + q11) - q15*x%d3val2 + q17*pow3(x%d1val2))
      unary%d2val2_d1Array(1:5) = q19*(-q12*q6*x%d1Array(1:5) - q25*x%d1Array(1:5) + q28*x%d1val2*(q13 + q6*x%d1val2_d1Array(1:5)) - q3*x%d2val2_d1Array(1:5))
   end function atanh_self
   
   function sqrt_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
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
      q0 = sqrt(x%val)
      q1 = 0.5_dp*powm1(q0)
      q2 = 2.0_dp*x%val
      q3 = q2*x%d2val1
      q4 = pow2(x%d1val1)
      q5 = 0.25_dp*powm1(pow3(sqrt(x%val)))
      q6 = x%d1val1*x%d1val2
      q7 = q2*x%d2val2
      q8 = pow2(x%d1val2)
      q9 = x%d1val1*x%val
      q10 = 4.0_dp*pow2(x%val)
      q11 = 0.125_dp*powm1(pow5(sqrt(x%val)))
      q12 = 4.0_dp*q9
      q13 = 3.0_dp*q4
      q14 = x%d1val2*x%val
      q15 = 4.0_dp*q14
      q16 = 3.0_dp*q8
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q5*(q3 - q4)
      unary%d1val1_d1val2 = q5*(q2*x%d1val1_d1val2 - q6)
      unary%d1val1_d1Array(1:5) = q5*(q2*x%d1val1_d1Array(1:5) - x%d1Array(1:5)*x%d1val1)
      unary%d2val2 = q5*(q7 - q8)
      unary%d1val2_d1Array(1:5) = q5*(q2*x%d1val2_d1Array(1:5) - x%d1Array(1:5)*x%d1val2)
      unary%d3val1 = q11*(-6.0_dp*q9*x%d2val1 + 3.0_dp*pow3(x%d1val1) + q10*x%d3val1)
      unary%d2val1_d1val2 = q11*(q10*x%d2val1_d1val2 - q12*x%d1val1_d1val2 + q13*x%d1val2 - q3*x%d1val2)
      unary%d2val1_d1Array(1:5) = q11*(q10*x%d2val1_d1Array(1:5) - q12*x%d1val1_d1Array(1:5) + q13*x%d1Array(1:5) - q3*x%d1Array(1:5))
      unary%d1val1_d2val2 = q11*(q10*x%d1val1_d2val2 - q15*x%d1val1_d1val2 + x%d1val1*(q16 - q7))
      unary%d1val1_d1val2_d1Array(1:5) = q11*(3.0_dp*q6*x%d1Array(1:5) + q10*x%d1val1_d1val2_d1Array(1:5) - q2*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2))
      unary%d3val2 = q11*(-6.0_dp*q14*x%d2val2 + 3.0_dp*pow3(x%d1val2) + q10*x%d3val2)
      unary%d2val2_d1Array(1:5) = q11*(q10*x%d2val2_d1Array(1:5) - q15*x%d1val2_d1Array(1:5) + q16*x%d1Array(1:5) - q7*x%d1Array(1:5))
   end function sqrt_self
   
   function pow2_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q5
      real(dp) :: q4(1:5)
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = 2.0_dp*x%val
      q1 = 2.0_dp*x%d1val1
      q2 = 2.0_dp*x%d1val2
      q3 = 4.0_dp*x%d1val1
      q4(1:5) = 2.0_dp*x%d1Array(1:5)
      q5 = 4.0_dp*x%d1val2
      unary%val = pow2(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = 2.0_dp*pow2(x%d1val1) + q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q1*x%d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5) + q1*x%d1Array(1:5)
      unary%d2val2 = 2.0_dp*pow2(x%d1val2) + q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5) + q2*x%d1Array(1:5)
      unary%d3val1 = 6.0_dp*x%d1val1*x%d2val1 + q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2 + q2*x%d2val1 + q3*x%d1val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5) + q3*x%d1val1_d1Array(1:5) + q4*x%d2val1
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2 + q1*x%d2val2 + q5*x%d1val1_d1val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5) + q1*x%d1val2_d1Array(1:5) + q2*x%d1val1_d1Array(1:5) + q4*x%d1val1_d1val2
      unary%d3val2 = 6.0_dp*x%d1val2*x%d2val2 + q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5) + q4*x%d2val2 + q5*x%d1val2_d1Array(1:5)
   end function pow2_self
   
   function pow3_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6(1:5)
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
      q5 = 2.0_dp*x%d1val1
      q6(1:5) = x%d1val1_d1Array(1:5)*x%val
      q7 = x%d2val2*x%val
      q8 = pow2(x%d1val2)
      q9 = 2.0_dp*q8 + q7
      q10(1:5) = x%d1val2_d1Array(1:5)*x%val
      q11(1:5) = x%d1Array(1:5)*x%d1val2
      q12 = 3.0_dp*q2
      q13 = 4.0_dp*x%d1val1
      q14 = 6.0_dp*x%d1val1
      unary%val = pow3(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(q4 + q5*x%d1val2)
      unary%d1val1_d1Array(1:5) = q3*(q5*x%d1Array(1:5) + q6)
      unary%d2val2 = q3*q9
      unary%d1val2_d1Array(1:5) = q3*(2.0_dp*q11 + q10)
      unary%d3val1 = 18.0_dp*q1*x%d1val1 + 6.0_dp*pow3(x%d1val1) + q0*x%d3val1
      unary%d2val1_d1val2 = q12*x%d1val2 + q3*(q13*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val)
      unary%d2val1_d1Array(1:5) = q12*x%d1Array(1:5) + q3*(q13*x%d1val1_d1Array(1:5) + x%d1Array(1:5)*x%d2val1 + x%d2val1_d1Array(1:5)*x%val)
      unary%d1val1_d2val2 = 12.0_dp*q4*x%d1val2 + q0*x%d1val1_d2val2 + q14*(q7 + q8)
      unary%d1val1_d1val2_d1Array(1:5) = 6.0_dp*q4*x%d1Array(1:5) + 6.0_dp*q6*x%d1val2 + q0*x%d1val1_d1val2_d1Array(1:5) + q10*q14 + q11*q14
      unary%d3val2 = 18.0_dp*q7*x%d1val2 + 6.0_dp*pow3(x%d1val2) + q0*x%d3val2
      unary%d2val2_d1Array(1:5) = 3.0_dp*q9*x%d1Array(1:5) + q3*(4.0_dp*x%d1val2*x%d1val2_d1Array(1:5) + x%d1Array(1:5)*x%d2val2 + x%d2val2_d1Array(1:5)*x%val)
   end function pow3_self
   
   function pow4_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12(1:5)
      real(dp) :: q11(1:5)
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7(1:5)
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
      q7(1:5) = x%d1val1_d1Array(1:5)*x%val
      q8 = x%d2val2*x%val
      q9 = pow2(x%d1val2)
      q10 = 3.0_dp*q9 + q8
      q11(1:5) = x%d1val2_d1Array(1:5)*x%val
      q12(1:5) = 3.0_dp*x%d1Array(1:5)
      q13 = 4.0_dp*x%val
      q14 = 2.0_dp*q2
      q15 = 6.0_dp*x%d1val1
      q16 = 6.0_dp*x%d1val2
      unary%val = pow4(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q2*q4
      unary%d1val1_d1val2 = q4*(q5 + q6*x%d1val2)
      unary%d1val1_d1Array(1:5) = q4*(q6*x%d1Array(1:5) + q7)
      unary%d2val2 = q10*q4
      unary%d1val2_d1Array(1:5) = q4*(q11 + q12*x%d1val2)
      unary%d3val1 = q13*(6.0_dp*pow3(x%d1val1) + 9.0_dp*q1*x%d1val1 + q3*x%d3val1)
      unary%d2val1_d1val2 = q13*(q14*x%d1val2 + x%val*(q15*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d2val1_d1Array(1:5) = q13*(q14*x%d1Array(1:5) + x%val*(q15*x%d1val1_d1Array(1:5) + x%d1Array(1:5)*x%d2val1 + x%d2val1_d1Array(1:5)*x%val))
      unary%d1val1_d2val2 = q13*(q16*q5 + q3*x%d1val1_d2val2 + q6*(2.0_dp*q9 + q8))
      unary%d1val1_d1val2_d1Array(1:5) = q13*(3.0_dp*q7*x%d1val2 + q11*q6 + q12*q5 + q15*x%d1Array(1:5)*x%d1val2 + q3*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q13*(6.0_dp*pow3(x%d1val2) + 9.0_dp*q8*x%d1val2 + q3*x%d3val2)
      unary%d2val2_d1Array(1:5) = q13*(2.0_dp*q10*x%d1Array(1:5) + x%val*(q16*x%d1val2_d1Array(1:5) + x%d1Array(1:5)*x%d2val2 + x%d2val2_d1Array(1:5)*x%val))
   end function pow4_self
   
   function pow5_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q17
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6(1:5)
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
      q6(1:5) = x%d1val1_d1Array(1:5)*x%val
      q7 = x%d2val2*x%val
      q8 = pow2(x%d1val2)
      q9 = 4.0_dp*q8 + q7
      q10(1:5) = x%d1val2_d1Array(1:5)*x%val
      q11(1:5) = 4.0_dp*x%d1Array(1:5)
      q12 = pow2(x%val)
      q13 = 12.0_dp*x%d1val1
      q14 = 5.0_dp*q12
      q15 = 3.0_dp*q2
      q16 = 8.0_dp*x%d1val1
      q17 = 8.0_dp*x%d1val2
      unary%val = pow5(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(q4 + q5*x%d1val2)
      unary%d1val1_d1Array(1:5) = q3*(q5*x%d1Array(1:5) + q6)
      unary%d2val2 = q3*q9
      unary%d1val2_d1Array(1:5) = q3*(q10 + q11*x%d1val2)
      unary%d3val1 = q14*(12.0_dp*pow3(x%d1val1) + q1*q13 + q12*x%d3val1)
      unary%d2val1_d1val2 = q14*(q15*x%d1val2 + x%val*(q16*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d2val1_d1Array(1:5) = q14*(q15*x%d1Array(1:5) + x%val*(q16*x%d1val1_d1Array(1:5) + x%d1Array(1:5)*x%d2val1 + x%d2val1_d1Array(1:5)*x%val))
      unary%d1val1_d2val2 = q14*(q12*x%d1val1_d2val2 + q17*q4 + q5*(3.0_dp*q8 + q7))
      unary%d1val1_d1val2_d1Array(1:5) = q14*(4.0_dp*q6*x%d1val2 + q10*q5 + q11*q4 + q12*x%d1val1_d1val2_d1Array(1:5) + q13*x%d1Array(1:5)*x%d1val2)
      unary%d3val2 = q14*(12.0_dp*q7*x%d1val2 + 12.0_dp*pow3(x%d1val2) + q12*x%d3val2)
      unary%d2val2_d1Array(1:5) = q14*(3.0_dp*q9*x%d1Array(1:5) + x%val*(q17*x%d1val2_d1Array(1:5) + x%d1Array(1:5)*x%d2val2 + x%d2val2_d1Array(1:5)*x%val))
   end function pow5_self
   
   function pow6_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6(1:5)
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
      q6(1:5) = x%d1val1_d1Array(1:5)*x%val
      q7 = x%d2val2*x%val
      q8 = pow2(x%d1val2)
      q9 = 5.0_dp*q8 + q7
      q10(1:5) = x%d1val2_d1Array(1:5)*x%val
      q11(1:5) = 5.0_dp*x%d1Array(1:5)
      q12 = pow2(x%val)
      q13 = 6.0_dp*pow3(x%val)
      q14 = 4.0_dp*q2
      q15 = 10.0_dp*x%d1val1
      q16 = 10.0_dp*x%d1val2
      unary%val = pow6(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(q4 + q5*x%d1val2)
      unary%d1val1_d1Array(1:5) = q3*(q5*x%d1Array(1:5) + q6)
      unary%d2val2 = q3*q9
      unary%d1val2_d1Array(1:5) = q3*(q10 + q11*x%d1val2)
      unary%d3val1 = q13*(15.0_dp*q1*x%d1val1 + 20.0_dp*pow3(x%d1val1) + q12*x%d3val1)
      unary%d2val1_d1val2 = q13*(q14*x%d1val2 + x%val*(q15*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d2val1_d1Array(1:5) = q13*(q14*x%d1Array(1:5) + x%val*(q15*x%d1val1_d1Array(1:5) + x%d1Array(1:5)*x%d2val1 + x%d2val1_d1Array(1:5)*x%val))
      unary%d1val1_d2val2 = q13*(q12*x%d1val1_d2val2 + q16*q4 + q5*(4.0_dp*q8 + q7))
      unary%d1val1_d1val2_d1Array(1:5) = q13*(20.0_dp*x%d1Array(1:5)*x%d1val1*x%d1val2 + 5.0_dp*q6*x%d1val2 + q10*q5 + q11*q4 + q12*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q13*(15.0_dp*q7*x%d1val2 + 20.0_dp*pow3(x%d1val2) + q12*x%d3val2)
      unary%d2val2_d1Array(1:5) = q13*(4.0_dp*q9*x%d1Array(1:5) + x%val*(q16*x%d1val2_d1Array(1:5) + x%d1Array(1:5)*x%d2val2 + x%d2val2_d1Array(1:5)*x%val))
   end function pow6_self
   
   function pow7_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6(1:5)
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
      q6(1:5) = x%d1val1_d1Array(1:5)*x%val
      q7 = x%d2val2*x%val
      q8 = pow2(x%d1val2)
      q9 = 6.0_dp*q8 + q7
      q10(1:5) = x%d1val2_d1Array(1:5)*x%val
      q11(1:5) = 6.0_dp*x%d1Array(1:5)
      q12 = pow2(x%val)
      q13 = 7.0_dp*pow4(x%val)
      q14 = 5.0_dp*q2
      q15 = 12.0_dp*x%d1val1
      q16 = 12.0_dp*x%d1val2
      unary%val = pow7(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(q4 + q5*x%d1val2)
      unary%d1val1_d1Array(1:5) = q3*(q5*x%d1Array(1:5) + q6)
      unary%d2val2 = q3*q9
      unary%d1val2_d1Array(1:5) = q3*(q10 + q11*x%d1val2)
      unary%d3val1 = q13*(18.0_dp*q1*x%d1val1 + 30.0_dp*pow3(x%d1val1) + q12*x%d3val1)
      unary%d2val1_d1val2 = q13*(q14*x%d1val2 + x%val*(q15*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d2val1_d1Array(1:5) = q13*(q14*x%d1Array(1:5) + x%val*(q15*x%d1val1_d1Array(1:5) + x%d1Array(1:5)*x%d2val1 + x%d2val1_d1Array(1:5)*x%val))
      unary%d1val1_d2val2 = q13*(q12*x%d1val1_d2val2 + q16*q4 + q5*(5.0_dp*q8 + q7))
      unary%d1val1_d1val2_d1Array(1:5) = q13*(30.0_dp*x%d1Array(1:5)*x%d1val1*x%d1val2 + 6.0_dp*q6*x%d1val2 + q10*q5 + q11*q4 + q12*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q13*(18.0_dp*q7*x%d1val2 + 30.0_dp*pow3(x%d1val2) + q12*x%d3val2)
      unary%d2val2_d1Array(1:5) = q13*(5.0_dp*q9*x%d1Array(1:5) + x%val*(q16*x%d1val2_d1Array(1:5) + x%d1Array(1:5)*x%d2val2 + x%d2val2_d1Array(1:5)*x%val))
   end function pow7_self
   
   function pow8_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q16
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11(1:5)
      real(dp) :: q10(1:5)
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7
      real(dp) :: q6(1:5)
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
      q6(1:5) = x%d1val1_d1Array(1:5)*x%val
      q7 = x%d2val2*x%val
      q8 = pow2(x%d1val2)
      q9 = 7.0_dp*q8 + q7
      q10(1:5) = x%d1val2_d1Array(1:5)*x%val
      q11(1:5) = 7.0_dp*x%d1Array(1:5)
      q12 = pow2(x%val)
      q13 = 8.0_dp*pow5(x%val)
      q14 = 6.0_dp*q2
      q15 = 14.0_dp*x%d1val1
      q16 = 14.0_dp*x%d1val2
      unary%val = pow8(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q3*(q4 + q5*x%d1val2)
      unary%d1val1_d1Array(1:5) = q3*(q5*x%d1Array(1:5) + q6)
      unary%d2val2 = q3*q9
      unary%d1val2_d1Array(1:5) = q3*(q10 + q11*x%d1val2)
      unary%d3val1 = q13*(21.0_dp*q1*x%d1val1 + 42.0_dp*pow3(x%d1val1) + q12*x%d3val1)
      unary%d2val1_d1val2 = q13*(q14*x%d1val2 + x%val*(q15*x%d1val1_d1val2 + x%d1val2*x%d2val1 + x%d2val1_d1val2*x%val))
      unary%d2val1_d1Array(1:5) = q13*(q14*x%d1Array(1:5) + x%val*(q15*x%d1val1_d1Array(1:5) + x%d1Array(1:5)*x%d2val1 + x%d2val1_d1Array(1:5)*x%val))
      unary%d1val1_d2val2 = q13*(q12*x%d1val1_d2val2 + q16*q4 + q5*(6.0_dp*q8 + q7))
      unary%d1val1_d1val2_d1Array(1:5) = q13*(42.0_dp*x%d1Array(1:5)*x%d1val1*x%d1val2 + 7.0_dp*q6*x%d1val2 + q10*q5 + q11*q4 + q12*x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q13*(21.0_dp*q7*x%d1val2 + 42.0_dp*pow3(x%d1val2) + q12*x%d3val2)
      unary%d2val2_d1Array(1:5) = q13*(6.0_dp*q9*x%d1Array(1:5) + x%val*(q16*x%d1val2_d1Array(1:5) + x%d1Array(1:5)*x%d2val2 + x%d2val2_d1Array(1:5)*x%val))
   end function pow8_self
   
   function abs_self(x) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q0
      q0 = sgn(x%val)
      unary%val = Abs(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function abs_self
   
   function add_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: binary
      binary%val = x%val + y%val
      binary%d1val1 = x%d1val1 + y%d1val1
      binary%d1val2 = x%d1val2 + y%d1val2
      binary%d1Array(1:5) = x%d1Array(1:5) + y%d1Array(1:5)
      binary%d2val1 = x%d2val1 + y%d2val1
      binary%d1val1_d1val2 = x%d1val1_d1val2 + y%d1val1_d1val2
      binary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5) + y%d1val1_d1Array(1:5)
      binary%d2val2 = x%d2val2 + y%d2val2
      binary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5) + y%d1val2_d1Array(1:5)
      binary%d3val1 = x%d3val1 + y%d3val1
      binary%d2val1_d1val2 = x%d2val1_d1val2 + y%d2val1_d1val2
      binary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5) + y%d2val1_d1Array(1:5)
      binary%d1val1_d2val2 = x%d1val1_d2val2 + y%d1val1_d2val2
      binary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5) + y%d1val1_d1val2_d1Array(1:5)
      binary%d3val2 = x%d3val2 + y%d3val2
      binary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5) + y%d2val2_d1Array(1:5)
   end function add_self
   
   function add_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      unary%val = x%val + y
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1Array(1:5) = x%d1Array(1:5)
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)
      unary%d2val2 = x%d2val2
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = x%d3val2
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)
   end function add_self_real
   
   function add_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      unary%val = x%val + z
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1Array(1:5) = x%d1Array(1:5)
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)
      unary%d2val2 = x%d2val2
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = x%d3val2
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)
   end function add_real_self
   
   function add_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val + y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1Array(1:5) = x%d1Array(1:5)
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)
      unary%d2val2 = x%d2val2
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = x%d3val2
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)
   end function add_self_int
   
   function add_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = x%val + y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1Array(1:5) = x%d1Array(1:5)
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)
      unary%d2val2 = x%d2val2
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = x%d3val2
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)
   end function add_int_self
   
   function sub_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: binary
      binary%val = x%val - y%val
      binary%d1val1 = x%d1val1 - y%d1val1
      binary%d1val2 = x%d1val2 - y%d1val2
      binary%d1Array(1:5) = x%d1Array(1:5) - y%d1Array(1:5)
      binary%d2val1 = x%d2val1 - y%d2val1
      binary%d1val1_d1val2 = x%d1val1_d1val2 - y%d1val1_d1val2
      binary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5) - y%d1val1_d1Array(1:5)
      binary%d2val2 = x%d2val2 - y%d2val2
      binary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5) - y%d1val2_d1Array(1:5)
      binary%d3val1 = x%d3val1 - y%d3val1
      binary%d2val1_d1val2 = x%d2val1_d1val2 - y%d2val1_d1val2
      binary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5) - y%d2val1_d1Array(1:5)
      binary%d1val1_d2val2 = x%d1val1_d2val2 - y%d1val1_d2val2
      binary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5) - y%d1val1_d1val2_d1Array(1:5)
      binary%d3val2 = x%d3val2 - y%d3val2
      binary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5) - y%d2val2_d1Array(1:5)
   end function sub_self
   
   function sub_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      unary%val = x%val - y
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1Array(1:5) = x%d1Array(1:5)
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)
      unary%d2val2 = x%d2val2
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = x%d3val2
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)
   end function sub_self_real
   
   function sub_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      unary%val = -x%val + z
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d1Array(1:5) = -x%d1Array(1:5)
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = -x%d1val1_d1Array(1:5)
      unary%d2val2 = -x%d2val2
      unary%d1val2_d1Array(1:5) = -x%d1val2_d1Array(1:5)
      unary%d3val1 = -x%d3val1
      unary%d2val1_d1val2 = -x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = -x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = -x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = -x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = -x%d3val2
      unary%d2val2_d1Array(1:5) = -x%d2val2_d1Array(1:5)
   end function sub_real_self
   
   function sub_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val - y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1Array(1:5) = x%d1Array(1:5)
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)
      unary%d2val2 = x%d2val2
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = x%d3val2
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)
   end function sub_self_int
   
   function sub_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = -x%val + y_dp
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d1Array(1:5) = -x%d1Array(1:5)
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = -x%d1val1_d1Array(1:5)
      unary%d2val2 = -x%d2val2
      unary%d1val2_d1Array(1:5) = -x%d1val2_d1Array(1:5)
      unary%d3val1 = -x%d3val1
      unary%d2val1_d1val2 = -x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = -x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = -x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = -x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = -x%d3val2
      unary%d2val2_d1Array(1:5) = -x%d2val2_d1Array(1:5)
   end function sub_int_self
   
   function mul_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: binary
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = 2.0_dp*x%d1val1
      q1 = 2.0_dp*y%d1val2
      q2 = 2.0_dp*y%d1val1
      q3 = 2.0_dp*x%d1val2
      binary%val = x%val*y%val
      binary%d1val1 = x%d1val1*y%val + x%val*y%d1val1
      binary%d1val2 = x%d1val2*y%val + x%val*y%d1val2
      binary%d1Array(1:5) = x%d1Array(1:5)*y%val + x%val*y%d1Array(1:5)
      binary%d2val1 = q0*y%d1val1 + x%d2val1*y%val + x%val*y%d2val1
      binary%d1val1_d1val2 = x%d1val1*y%d1val2 + x%d1val1_d1val2*y%val + x%d1val2*y%d1val1 + x%val*y%d1val1_d1val2
      binary%d1val1_d1Array(1:5) = x%d1Array(1:5)*y%d1val1 + x%d1val1*y%d1Array(1:5) + x%d1val1_d1Array(1:5)*y%val + x%val*y%d1val1_d1Array(1:5)
      binary%d2val2 = q1*x%d1val2 + x%d2val2*y%val + x%val*y%d2val2
      binary%d1val2_d1Array(1:5) = x%d1Array(1:5)*y%d1val2 + x%d1val2*y%d1Array(1:5) + x%d1val2_d1Array(1:5)*y%val + x%val*y%d1val2_d1Array(1:5)
      binary%d3val1 = 3.0_dp*x%d1val1*y%d2val1 + 3.0_dp*x%d2val1*y%d1val1 + x%d3val1*y%val + x%val*y%d3val1
      binary%d2val1_d1val2 = q0*y%d1val1_d1val2 + q2*x%d1val1_d1val2 + x%d1val2*y%d2val1 + x%d2val1*y%d1val2 + x%d2val1_d1val2*y%val + x%val*y%d2val1_d1val2
      binary%d2val1_d1Array(1:5) = q0*y%d1val1_d1Array(1:5) + q2*x%d1val1_d1Array(1:5) + x%d1Array(1:5)*y%d2val1 + x%d2val1*y%d1Array(1:5) + x%d2val1_d1Array(1:5)*y%val + x%val*y%d2val1_d1Array(1:5)
      binary%d1val1_d2val2 = q1*x%d1val1_d1val2 + q3*y%d1val1_d1val2 + x%d1val1*y%d2val2 + x%d1val1_d2val2*y%val + x%d2val2*y%d1val1 + x%val*y%d1val1_d2val2
      binary%d1val1_d1val2_d1Array(1:5) = x%d1Array(1:5)*y%d1val1_d1val2 + x%d1val1*y%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*y%d1val2 + x%d1val1_d1val2*y%d1Array(1:5) + x%d1val1_d1val2_d1Array(1:5)*y%val + x%d1val2*y%d1val1_d1Array(1:5) + x%d1val2_d1Array(1:5)*y%d1val1 + x%val*y%d1val1_d1val2_d1Array(1:5)
      binary%d3val2 = 3.0_dp*x%d1val2*y%d2val2 + 3.0_dp*x%d2val2*y%d1val2 + x%d3val2*y%val + x%val*y%d3val2
      binary%d2val2_d1Array(1:5) = q1*x%d1val2_d1Array(1:5) + q3*y%d1val2_d1Array(1:5) + x%d1Array(1:5)*y%d2val2 + x%d2val2*y%d1Array(1:5) + x%d2val2_d1Array(1:5)*y%val + x%val*y%d2val2_d1Array(1:5)
   end function mul_self
   
   function mul_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      unary%val = x%val*y
      unary%d1val1 = x%d1val1*y
      unary%d1val2 = x%d1val2*y
      unary%d1Array(1:5) = x%d1Array(1:5)*y
      unary%d2val1 = x%d2val1*y
      unary%d1val1_d1val2 = x%d1val1_d1val2*y
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)*y
      unary%d2val2 = x%d2val2*y
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)*y
      unary%d3val1 = x%d3val1*y
      unary%d2val1_d1val2 = x%d2val1_d1val2*y
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)*y
      unary%d1val1_d2val2 = x%d1val1_d2val2*y
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)*y
      unary%d3val2 = x%d3val2*y
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)*y
   end function mul_self_real
   
   function mul_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      unary%val = x%val*z
      unary%d1val1 = x%d1val1*z
      unary%d1val2 = x%d1val2*z
      unary%d1Array(1:5) = x%d1Array(1:5)*z
      unary%d2val1 = x%d2val1*z
      unary%d1val1_d1val2 = x%d1val1_d1val2*z
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)*z
      unary%d2val2 = x%d2val2*z
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)*z
      unary%d3val1 = x%d3val1*z
      unary%d2val1_d1val2 = x%d2val1_d1val2*z
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)*z
      unary%d1val1_d2val2 = x%d1val1_d2val2*z
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)*z
      unary%d3val2 = x%d3val2*z
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)*z
   end function mul_real_self
   
   function mul_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val*y_dp
      unary%d1val1 = x%d1val1*y_dp
      unary%d1val2 = x%d1val2*y_dp
      unary%d1Array(1:5) = x%d1Array(1:5)*y_dp
      unary%d2val1 = x%d2val1*y_dp
      unary%d1val1_d1val2 = x%d1val1_d1val2*y_dp
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)*y_dp
      unary%d2val2 = x%d2val2*y_dp
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)*y_dp
      unary%d3val1 = x%d3val1*y_dp
      unary%d2val1_d1val2 = x%d2val1_d1val2*y_dp
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)*y_dp
      unary%d1val1_d2val2 = x%d1val1_d2val2*y_dp
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)*y_dp
      unary%d3val2 = x%d3val2*y_dp
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)*y_dp
   end function mul_self_int
   
   function mul_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = x%val*y_dp
      unary%d1val1 = x%d1val1*y_dp
      unary%d1val2 = x%d1val2*y_dp
      unary%d1Array(1:5) = x%d1Array(1:5)*y_dp
      unary%d2val1 = x%d2val1*y_dp
      unary%d1val1_d1val2 = x%d1val1_d1val2*y_dp
      unary%d1val1_d1Array(1:5) = x%d1val1_d1Array(1:5)*y_dp
      unary%d2val2 = x%d2val2*y_dp
      unary%d1val2_d1Array(1:5) = x%d1val2_d1Array(1:5)*y_dp
      unary%d3val1 = x%d3val1*y_dp
      unary%d2val1_d1val2 = x%d2val1_d1val2*y_dp
      unary%d2val1_d1Array(1:5) = x%d2val1_d1Array(1:5)*y_dp
      unary%d1val1_d2val2 = x%d1val1_d2val2*y_dp
      unary%d1val1_d1val2_d1Array(1:5) = x%d1val1_d1val2_d1Array(1:5)*y_dp
      unary%d3val2 = x%d3val2*y_dp
      unary%d2val2_d1Array(1:5) = x%d2val2_d1Array(1:5)*y_dp
   end function mul_int_self
   
   function div_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: binary
      real(dp) :: q43
      real(dp) :: q42(1:5)
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33(1:5)
      real(dp) :: q32
      real(dp) :: q31(1:5)
      real(dp) :: q30
      real(dp) :: q29
      real(dp) :: q28
      real(dp) :: q27(1:5)
      real(dp) :: q26
      real(dp) :: q25(1:5)
      real(dp) :: q24(1:5)
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
      real(dp) :: q9(1:5)
      real(dp) :: q8(1:5)
      real(dp) :: q7
      real(dp) :: q6
      real(dp) :: q5
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(y%val)
      q1 = q0*x%val
      q2 = q0*x%d1val1
      q3 = powm1(pow2(y%val))
      q4 = q3*x%val
      q5 = q4*y%d1val1
      q6 = q0*x%d1val2
      q7 = q4*y%d1val2
      q8(1:5) = q0*x%d1Array(1:5)
      q9(1:5) = q4*y%d1Array(1:5)
      q10 = 2.0_dp*q2
      q11 = pow2(y%d1val1)
      q12 = 2.0_dp*q0
      q13 = -q11*q12 + y%d2val1
      q14 = -q1*q13 - q10*y%d1val1 + x%d2val1
      q15 = q0*x%d1val1_d1val2
      q16 = q3*x%d1val1
      q17 = q16*y%d1val2
      q18 = q3*y%d1val1
      q19 = q18*x%d1val2
      q20 = 2.0_dp*y%d1val2
      q21 = powm1(pow3(y%val))
      q22 = q21*x%val
      q23 = q20*q22
      q24(1:5) = q0*x%d1val1_d1Array(1:5)
      q25(1:5) = q16*y%d1Array(1:5)
      q26 = 2.0_dp*y%d1val1
      q27(1:5) = q26*y%d1Array(1:5)
      q28 = pow2(y%d1val2)
      q29 = -q12*q28 + y%d2val2
      q30 = -q1*q29 - q20*q6 + x%d2val2
      q31(1:5) = q0*x%d1val2_d1Array(1:5)
      q32 = q3*y%d1val2
      q33(1:5) = q3*y%d1Array(1:5)
      q34 = q0*y%d1val1
      q35 = 6.0_dp*q3
      q36 = 4.0_dp*q34
      q37 = 2.0_dp*q11
      q38 = 2.0_dp*q6
      q39 = 2.0_dp*q28
      q40 = 6.0_dp*y%d1val1
      q41 = q20*q21
      q42(1:5) = 2.0_dp*y%d1Array(1:5)
      q43 = q0*y%d1val2
      binary%val = q1
      binary%d1val1 = q2 - q5
      binary%d1val2 = q6 - q7
      binary%d1Array(1:5) = q8 - q9
      binary%d2val1 = q0*q14
      binary%d1val1_d1val2 = q15 - q17 - q19 + q23*y%d1val1 - q4*y%d1val1_d1val2
      binary%d1val1_d1Array(1:5) = -q18*x%d1Array(1:5) + q22*q27 + q24 - q25 - q4*y%d1val1_d1Array(1:5)
      binary%d2val2 = q0*q30
      binary%d1val2_d1Array(1:5) = q23*y%d1Array(1:5) + q31 - q32*x%d1Array(1:5) - q33*x%d1val2 - q4*y%d1val2_d1Array(1:5)
      binary%d3val1 = q0*(-3.0_dp*q13*q2 - 3.0_dp*q34*x%d2val1 - q1*(-6.0_dp*q34*y%d2val1 + q35*pow3(y%d1val1) + y%d3val1) + x%d3val1)
      binary%d2val1_d1val2 = q0*(-q1*(q32*q37 - q36*y%d1val1_d1val2 + y%d2val1_d1val2) - q10*y%d1val1_d1val2 - q13*q6 + q13*q7 - q15*q26 + q17*q26 + x%d2val1_d1val2) - q14*q32
      binary%d2val1_d1Array(1:5) = q0*(-q1*(q33*q37 - q36*y%d1val1_d1Array(1:5) + y%d2val1_d1Array(1:5)) - q10*y%d1val1_d1Array(1:5) - q13*q8 + q13*q9 - q24*q26 + q25*q26 + x%d2val1_d1Array(1:5)) - q14*q33
      binary%d1val1_d2val2 = q0*(2.0_dp*q5*y%d2val2 + 4.0_dp*q19*y%d1val2 + 4.0_dp*q7*y%d1val1_d1val2 - q1*y%d1val1_d2val2 - q15*q20 + q16*q39 - q2*y%d2val2 - q22*q28*q40 - q34*x%d2val2 - q38*y%d1val1_d1val2 + x%d1val1_d2val2)
      binary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5) - q16*y%d1val2_d1Array(1:5) - q18*x%d1val2_d1Array(1:5) + q21*q27*x%d1val2 + q22*q26*y%d1val2_d1Array(1:5) + q22*q42*y%d1val1_d1val2 + q23*y%d1val1_d1Array(1:5) - q3*x%d1Array(1:5)*y%d1val1_d1val2 - q3*x%d1val2*y%d1val1_d1Array(1:5) - q32*x%d1val1_d1Array(1:5) - q33*x%d1val1_d1val2 - q4*y%d1val1_d1val2_d1Array(1:5) - q40*x%val*y%d1Array(1:5)*y%d1val2*powm1(pow4(y%val)) + q41*x%d1Array(1:5)*y%d1val1 + q41*x%d1val1*y%d1Array(1:5)
      binary%d3val2 = q0*(-3.0_dp*q29*q6 - 3.0_dp*q43*x%d2val2 - q1*(-6.0_dp*q43*y%d2val2 + q35*pow3(y%d1val2) + y%d3val2) + x%d3val2)
      binary%d2val2_d1Array(1:5) = q0*(-q1*(-4.0_dp*q43*y%d1val2_d1Array(1:5) + q33*q39 + y%d2val2_d1Array(1:5)) - q20*q31 - q29*q8 + q29*q9 + q32*q42*x%d1val2 - q38*y%d1val2_d1Array(1:5) + x%d2val2_d1Array(1:5)) - q30*q33
   end function div_self
   
   function div_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q0
      q0 = powm1(y)
      unary%val = q0*x%val
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function div_self_real
   
   function div_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q18(1:5)
      real(dp) :: q17
      real(dp) :: q16(1:5)
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12(1:5)
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7(1:5)
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
      q7(1:5) = x%d1val1_d1Array(1:5)*x%val
      q8 = x%d2val2*x%val
      q9 = -q8
      q10 = pow2(x%d1val2)
      q11 = 2.0_dp*x%d1val2
      q12(1:5) = x%d1val2_d1Array(1:5)*x%val
      q13 = z*powm1(pow4(x%val))
      q14 = 4.0_dp*x%d1val1
      q15 = 6.0_dp*q3
      q16(1:5) = 2.0_dp*x%d1Array(1:5)
      q17 = 4.0_dp*x%d1val2
      q18(1:5) = 6.0_dp*x%d1Array(1:5)
      unary%val = z*powm1(x%val)
      unary%d1val1 = -q1*x%d1val1
      unary%d1val2 = -q1*x%d1val2
      unary%d1Array(1:5) = -q1*x%d1Array(1:5)
      unary%d2val1 = q4*(2.0_dp*q3 - q2)
      unary%d1val1_d1val2 = q4*(q5*x%d1val2 - q6)
      unary%d1val1_d1Array(1:5) = q4*(q5*x%d1Array(1:5) - q7)
      unary%d2val2 = q4*(2.0_dp*q10 + q9)
      unary%d1val2_d1Array(1:5) = q4*(q11*x%d1Array(1:5) - q12)
      unary%d3val1 = -q13*(-6.0_dp*q2*x%d1val1 + 6.0_dp*pow3(x%d1val1) + q0*x%d3val1)
      unary%d2val1_d1val2 = -q13*(q0*x%d2val1_d1val2 - q11*q2 - q14*q6 + q15*x%d1val2)
      unary%d2val1_d1Array(1:5) = -q13*(q0*x%d2val1_d1Array(1:5) - q14*q7 + q15*x%d1Array(1:5) - q16*q2)
      unary%d1val1_d2val2 = -q13*(q0*x%d1val1_d2val2 - q17*q6 + q5*(3.0_dp*q10 + q9))
      unary%d1val1_d1val2_d1Array(1:5) = q13*(2.0_dp*x%val*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2) - q0*x%d1val1_d1val2_d1Array(1:5) - q18*x%d1val1*x%d1val2)
      unary%d3val2 = -q13*(-6.0_dp*q8*x%d1val2 + 6.0_dp*pow3(x%d1val2) + q0*x%d3val2)
      unary%d2val2_d1Array(1:5) = -q13*(q0*x%d2val2_d1Array(1:5) + q10*q18 - q12*q17 - q16*q8)
   end function div_real_self
   
   function div_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = y
      q0 = powm1(y_dp)
      unary%val = q0*x%val
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function div_self_int
   
   function div_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q18(1:5)
      real(dp) :: q17
      real(dp) :: q16(1:5)
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12(1:5)
      real(dp) :: q11
      real(dp) :: q10
      real(dp) :: q9
      real(dp) :: q8
      real(dp) :: q7(1:5)
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
      q7(1:5) = x%d1val1_d1Array(1:5)*x%val
      q8 = x%d2val2*x%val
      q9 = -q8
      q10 = pow2(x%d1val2)
      q11 = 2.0_dp*x%d1val2
      q12(1:5) = x%d1val2_d1Array(1:5)*x%val
      q13 = y_dp*powm1(pow4(x%val))
      q14 = 4.0_dp*x%d1val1
      q15 = 6.0_dp*q3
      q16(1:5) = 2.0_dp*x%d1Array(1:5)
      q17 = 4.0_dp*x%d1val2
      q18(1:5) = 6.0_dp*x%d1Array(1:5)
      unary%val = y_dp*powm1(x%val)
      unary%d1val1 = -q1*x%d1val1
      unary%d1val2 = -q1*x%d1val2
      unary%d1Array(1:5) = -q1*x%d1Array(1:5)
      unary%d2val1 = q4*(2.0_dp*q3 - q2)
      unary%d1val1_d1val2 = q4*(q5*x%d1val2 - q6)
      unary%d1val1_d1Array(1:5) = q4*(q5*x%d1Array(1:5) - q7)
      unary%d2val2 = q4*(2.0_dp*q10 + q9)
      unary%d1val2_d1Array(1:5) = q4*(q11*x%d1Array(1:5) - q12)
      unary%d3val1 = -q13*(-6.0_dp*q2*x%d1val1 + 6.0_dp*pow3(x%d1val1) + q0*x%d3val1)
      unary%d2val1_d1val2 = -q13*(q0*x%d2val1_d1val2 - q11*q2 - q14*q6 + q15*x%d1val2)
      unary%d2val1_d1Array(1:5) = -q13*(q0*x%d2val1_d1Array(1:5) - q14*q7 + q15*x%d1Array(1:5) - q16*q2)
      unary%d1val1_d2val2 = -q13*(q0*x%d1val1_d2val2 - q17*q6 + q5*(3.0_dp*q10 + q9))
      unary%d1val1_d1val2_d1Array(1:5) = q13*(2.0_dp*x%val*(x%d1Array(1:5)*x%d1val1_d1val2 + x%d1val1*x%d1val2_d1Array(1:5) + x%d1val1_d1Array(1:5)*x%d1val2) - q0*x%d1val1_d1val2_d1Array(1:5) - q18*x%d1val1*x%d1val2)
      unary%d3val2 = -q13*(-6.0_dp*q8*x%d1val2 + 6.0_dp*pow3(x%d1val2) + q0*x%d3val2)
      unary%d2val2_d1Array(1:5) = -q13*(q0*x%d2val2_d1Array(1:5) + q10*q18 - q12*q17 - q16*q8)
   end function div_int_self
   
   function pow_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: binary
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69(1:5)
      real(dp) :: q68(1:5)
      real(dp) :: q67(1:5)
      real(dp) :: q66(1:5)
      real(dp) :: q65(1:5)
      real(dp) :: q64
      real(dp) :: q63
      real(dp) :: q62
      real(dp) :: q61
      real(dp) :: q60
      real(dp) :: q59
      real(dp) :: q58
      real(dp) :: q57
      real(dp) :: q56
      real(dp) :: q55
      real(dp) :: q54(1:5)
      real(dp) :: q53(1:5)
      real(dp) :: q52
      real(dp) :: q51(1:5)
      real(dp) :: q50(1:5)
      real(dp) :: q49
      real(dp) :: q48(1:5)
      real(dp) :: q47(1:5)
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37(1:5)
      real(dp) :: q36(1:5)
      real(dp) :: q35(1:5)
      real(dp) :: q34(1:5)
      real(dp) :: q33(1:5)
      real(dp) :: q32(1:5)
      real(dp) :: q31
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
      real(dp) :: q12(1:5)
      real(dp) :: q11(1:5)
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
      q0 = pow(x%val, y%val)
      q1 = log(x%val)
      q2 = q1*y%d1val1
      q3 = powm1(x%val)
      q4 = q3*y%val
      q5 = q4*x%d1val1
      q6 = q0*(q2 + q5)
      q7 = q1*y%d1val2
      q8 = q4*x%d1val2
      q9 = q7 + q8
      q10 = q0*q9
      q11(1:5) = q1*y%d1Array(1:5) + q4*x%d1Array(1:5)
      q12(1:5) = q0*q11
      q13 = q1*y%d2val1
      q14 = q4*x%d2val1
      q15 = q3*x%d1val1
      q16 = q15*y%d1val1
      q17 = pow2(x%d1val1)
      q18 = powm1(pow2(x%val))
      q19 = q18*y%val
      q20 = q17*q19
      q21 = 2.0_dp*q2 + 2.0_dp*q5
      q22 = 0.25_dp*pow2(q21) + 2.0_dp*q16 + q13 + q14 - q20
      q23 = q6*q9
      q24 = q1*y%d1val1_d1val2
      q25 = q15*y%d1val2
      q26 = q4*x%d1val1_d1val2
      q27 = q3*y%d1val1
      q28 = q27*x%d1val2
      q29 = q19*x%d1val1
      q30 = q29*x%d1val2
      q31 = q24 + q25 + q26 + q28 - q30
      q32(1:5) = q1*y%d1val1_d1Array(1:5)
      q33(1:5) = q27*x%d1Array(1:5)
      q34(1:5) = q15*y%d1Array(1:5)
      q35(1:5) = q4*x%d1val1_d1Array(1:5)
      q36(1:5) = q29*x%d1Array(1:5)
      q37(1:5) = q32 + q33 + q34 + q35 - q36
      q38 = q1*y%d2val2
      q39 = q4*x%d2val2
      q40 = q3*y%d1val2
      q41 = q40*x%d1val2
      q42 = pow2(x%d1val2)
      q43 = q19*q42
      q44 = 2.0_dp*q7 + 2.0_dp*q8
      q45 = pow2(q44)
      q46 = 0.25_dp*q45 + 2.0_dp*q41 + q38 + q39 - q43
      q47(1:5) = q1*y%d1val2_d1Array(1:5)
      q48(1:5) = q40*x%d1Array(1:5)
      q49 = q3*x%d1val2
      q50(1:5) = q49*y%d1Array(1:5)
      q51(1:5) = q4*x%d1val2_d1Array(1:5)
      q52 = q19*x%d1val2
      q53(1:5) = q52*x%d1Array(1:5)
      q54(1:5) = q47 + q48 + q50 + q51 - q53
      q55 = 3.0_dp*x%d2val1
      q56 = q17*q18
      q57 = 2.0_dp*y%val*powm1(pow3(x%val))
      q58 = 2.0_dp*q15
      q59 = 2.0_dp*q27
      q60 = 2.0_dp*q29
      q61 = x%d1val2*y%d1val1
      q62 = q18*x%d1val1
      q63 = q17*q57
      q64 = 0.25_dp*q21
      q65(1:5) = q3*x%d1Array(1:5)
      q66(1:5) = q3*y%d1Array(1:5)
      q67(1:5) = q18*x%d1Array(1:5)
      q68(1:5) = q67*x%d1val1
      q69(1:5) = q19*x%d1Array(1:5)
      q70 = 2.0_dp*q40
      q71 = 2.0_dp*q49
      q72 = q62*x%d1val2
      q73 = 2.0_dp*y%d1val2
      q74 = 2.0_dp*q52
      q75 = q18*q42
      q76 = q42*q57
      q77 = 3.0_dp*x%d2val2
      binary%val = q0
      binary%d1val1 = q6
      binary%d1val2 = q10
      binary%d1Array(1:5) = q12
      binary%d2val1 = q0*q22
      binary%d1val1_d1val2 = q0*q31 + q23
      binary%d1val1_d1Array(1:5) = q0*q37 + q11*q6
      binary%d2val2 = q0*q46
      binary%d1val2_d1Array(1:5) = q0*q54 + q10*q11
      binary%d3val1 = q0*(-3.0_dp*q56*y%d1val1 + 0.125_dp*pow3(q21) + 0.75_dp*q21*(-2.0_dp*q20 + 2.0_dp*q13 + 2.0_dp*q14 + 4.0_dp*q16) + 3.0_dp*q15*y%d2val1 + q1*y%d3val1 + q27*q55 - q29*q55 + q4*x%d3val1 + q57*pow3(x%d1val1))
      binary%d2val1_d1val2 = q0*(-2.0_dp*q61*q62 + q1*y%d2val1_d1val2 + q4*x%d2val1_d1val2 + q40*x%d2val1 + q49*y%d2val1 - q52*x%d2val1 - q56*y%d1val2 + q58*y%d1val1_d1val2 + q59*x%d1val1_d1val2 - q60*x%d1val1_d1val2 + q63*x%d1val2 + q64*(-4.0_dp*q30 + 4.0_dp*q24 + 4.0_dp*q25 + 4.0_dp*q26 + 4.0_dp*q28)) + q10*q22
      binary%d2val1_d1Array(1:5) = q0*(-2.0_dp*q68*y%d1val1 + q1*y%d2val1_d1Array(1:5) + q4*x%d2val1_d1Array(1:5) - q56*y%d1Array(1:5) + q58*y%d1val1_d1Array(1:5) + q59*x%d1val1_d1Array(1:5) - q60*x%d1val1_d1Array(1:5) + q63*x%d1Array(1:5) + q64*(-4.0_dp*q36 + 4.0_dp*q32 + 4.0_dp*q33 + 4.0_dp*q34 + 4.0_dp*q35) + q65*y%d2val1 + q66*x%d2val1 - q69*x%d2val1) + q12*q22
      binary%d1val1_d2val2 = q0*(0.125_dp*q21*(-4.0_dp*q43 + 4.0_dp*q38 + 4.0_dp*q39 + 8.0_dp*q41 + q45) + 0.5_dp*q44*(-2.0_dp*q30 + 2.0_dp*q24 + 2.0_dp*q25 + 2.0_dp*q26 + 2.0_dp*q28) + q1*y%d1val1_d2val2 + q15*y%d2val2 + q27*x%d2val2 - q29*x%d2val2 + q4*x%d1val1_d2val2 + q70*x%d1val1_d1val2 + q71*y%d1val1_d1val2 - q72*q73 - q74*x%d1val1_d1val2 - q75*y%d1val1 + q76*x%d1val1)
      binary%d1val1_d1val2_d1Array(1:5) = q0*(q1*y%d1val1_d1val2_d1Array(1:5) + q15*y%d1val2_d1Array(1:5) + q27*x%d1val2_d1Array(1:5) - q29*x%d1val2_d1Array(1:5) + q4*x%d1val1_d1val2_d1Array(1:5) + q40*x%d1val1_d1Array(1:5) + q49*y%d1val1_d1Array(1:5) - q52*x%d1val1_d1Array(1:5) + q57*x%d1Array(1:5)*x%d1val1*x%d1val2 - q61*q67 + q65*y%d1val1_d1val2 + q66*x%d1val1_d1val2 - q68*y%d1val2 - q69*x%d1val1_d1val2 - q72*y%d1Array(1:5)) + q10*q37 + q11*q23 + q12*q31 + q54*q6
      binary%d3val2 = q0*(-3.0_dp*q75*y%d1val2 + 0.125_dp*pow3(q44) + 0.75_dp*q44*(-2.0_dp*q43 + 2.0_dp*q38 + 2.0_dp*q39 + 4.0_dp*q41) + 3.0_dp*q49*y%d2val2 + q1*y%d3val2 + q4*x%d3val2 + q40*q77 - q52*q77 + q57*pow3(x%d1val2))
      binary%d2val2_d1Array(1:5) = q0*(0.25_dp*q44*(-4.0_dp*q53 + 4.0_dp*q47 + 4.0_dp*q48 + 4.0_dp*q50 + 4.0_dp*q51) + q1*y%d2val2_d1Array(1:5) + q4*x%d2val2_d1Array(1:5) + q65*y%d2val2 + q66*x%d2val2 - q67*q73*x%d1val2 - q69*x%d2val2 + q70*x%d1val2_d1Array(1:5) + q71*y%d1val2_d1Array(1:5) - q74*x%d1val2_d1Array(1:5) - q75*y%d1Array(1:5) + q76*x%d1Array(1:5)) + q12*q46
   end function pow_self
   
   function pow_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q33(1:5)
      real(dp) :: q32(1:5)
      real(dp) :: q31(1:5)
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
      real(dp) :: q17(1:5)
      real(dp) :: q16(1:5)
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
      real(dp) :: q9(1:5)
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
      q4 = q2 + q3*y - q3
      q5 = pow(x%val, -2.0_dp + y)
      q6 = q5*y
      q7 = x%d1val1_d1val2*x%val
      q8 = x%d1val1*x%d1val2
      q9(1:5) = x%d1val1_d1Array(1:5)*x%val
      q10(1:5) = x%d1Array(1:5)*x%d1val1
      q11 = x%d2val2*x%val
      q12 = pow2(x%d1val2)
      q13 = q12*y
      q14 = q11 - q12 + q13
      q15 = q14*y
      q16(1:5) = x%d1val2_d1Array(1:5)*x%val
      q17(1:5) = x%d1Array(1:5)*x%d1val2
      q18 = pow2(x%val)
      q19 = q0*x%d1val1
      q20 = -3.0_dp*y + 2.0_dp + pow2(y)
      q21 = y*pow(x%val, -3.0_dp + y)
      q22 = 2.0_dp*q19
      q23 = 1 - y
      q24 = q23*q3
      q25 = q0*q4
      q26 = pow(x%val, 3.0_dp + y)
      q27 = y*powm1(pow6(x%val))
      q28 = q26*q27
      q29 = q0*x%d1val2
      q30 = 2.0_dp*q29
      q31(1:5) = x%d1Array(1:5)*x%d1val1_d1val2
      q32(1:5) = x%d1val1*x%d1val2_d1Array(1:5)
      q33(1:5) = x%d1val1_d1Array(1:5)*x%d1val2
      unary%val = pow(x%val, y)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q4*q6
      unary%d1val1_d1val2 = q6*(q7 + q8*y - q8)
      unary%d1val1_d1Array(1:5) = q6*(q10*y - q10 + q9)
      unary%d2val2 = q15*q5
      unary%d1val2_d1Array(1:5) = q6*(q16 + q17*y - q17)
      unary%d3val1 = q21*(3.0_dp*q19*q2 + q18*x%d3val1 + q20*pow3(x%d1val1))
      unary%d2val1_d1val2 = q28*(q18*x%d2val1_d1val2 + q22*q7 + q24*x%d1val2 + q25*x%d1val2)
      unary%d2val1_d1Array(1:5) = q28*(q18*x%d2val1_d1Array(1:5) + q22*q9 + q24*x%d1Array(1:5) + q25*x%d1Array(1:5))
      unary%d1val1_d2val2 = q21*(q18*x%d1val1_d2val2 + q30*q7 - x%d1val1*(-2.0_dp*q12 + 2.0_dp*q13 + q11 - q15))
      unary%d1val1_d1val2_d1Array(1:5) = q27*(q20*q26*q8*x%d1Array(1:5) + x%d1val1_d1val2_d1Array(1:5)*pow(x%val, 5.0_dp + y) + (q31*y - q31 + q32*y - q32 + q33*y - q33)*pow(x%val, 4.0_dp + y))
      unary%d3val2 = q21*(3.0_dp*q11*q29 + q18*x%d3val2 + q20*pow3(x%d1val2))
      unary%d2val2_d1Array(1:5) = q28*(q0*q14*x%d1Array(1:5) + q12*q23*x%d1Array(1:5) + q16*q30 + q18*x%d2val2_d1Array(1:5))
   end function pow_self_real
   
   function pow_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q10
      real(dp) :: q9(1:5)
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
      q4 = q1*x%d1val1
      q5 = q1*pow2(x%d1val2) + x%d2val2
      q6 = q1*x%d1val2
      q7 = pow2(q1)
      q8 = 2.0_dp*q4
      q9(1:5) = q1*x%d1Array(1:5)
      q10 = 2.0_dp*q6
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1Array(1:5) = q2*x%d1Array(1:5)
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q2*(q4*x%d1val2 + x%d1val1_d1val2)
      unary%d1val1_d1Array(1:5) = q2*(q4*x%d1Array(1:5) + x%d1val1_d1Array(1:5))
      unary%d2val2 = q2*q5
      unary%d1val2_d1Array(1:5) = q2*(q6*x%d1Array(1:5) + x%d1val2_d1Array(1:5))
      unary%d3val1 = q2*(3.0_dp*q4*x%d2val1 + q7*pow3(x%d1val1) + x%d3val1)
      unary%d2val1_d1val2 = q2*(q3*q6 + q8*x%d1val1_d1val2 + x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q2*(q3*q9 + q8*x%d1val1_d1Array(1:5) + x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q2*(q10*x%d1val1_d1val2 + q4*q5 + x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = q2*(q4*x%d1val2_d1Array(1:5) + q6*x%d1val1_d1Array(1:5) + q7*x%d1Array(1:5)*x%d1val1*x%d1val2 + q9*x%d1val1_d1val2 + x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q2*(3.0_dp*q6*x%d2val2 + q7*pow3(x%d1val2) + x%d3val2)
      unary%d2val2_d1Array(1:5) = q2*(q10*x%d1val2_d1Array(1:5) + q5*q9 + x%d2val2_d1Array(1:5))
   end function pow_real_self
   
   function pow_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q33(1:5)
      real(dp) :: q32(1:5)
      real(dp) :: q31(1:5)
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
      real(dp) :: q17(1:5)
      real(dp) :: q16(1:5)
      real(dp) :: q15
      real(dp) :: q14
      real(dp) :: q13
      real(dp) :: q12
      real(dp) :: q11
      real(dp) :: q10(1:5)
      real(dp) :: q9(1:5)
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
      q4 = q2 + q3*y_dp - q3
      q5 = pow(x%val, -2.0_dp + y_dp)
      q6 = q5*y_dp
      q7 = x%d1val1_d1val2*x%val
      q8 = x%d1val1*x%d1val2
      q9(1:5) = x%d1val1_d1Array(1:5)*x%val
      q10(1:5) = x%d1Array(1:5)*x%d1val1
      q11 = x%d2val2*x%val
      q12 = pow2(x%d1val2)
      q13 = q12*y_dp
      q14 = q11 - q12 + q13
      q15 = q14*y_dp
      q16(1:5) = x%d1val2_d1Array(1:5)*x%val
      q17(1:5) = x%d1Array(1:5)*x%d1val2
      q18 = pow2(x%val)
      q19 = q0*x%d1val1
      q20 = -3.0_dp*y_dp + 2.0_dp + pow2(y_dp)
      q21 = y_dp*pow(x%val, -3.0_dp + y_dp)
      q22 = 2.0_dp*q19
      q23 = 1 - y_dp
      q24 = q23*q3
      q25 = q0*q4
      q26 = pow(x%val, 3.0_dp + y_dp)
      q27 = y_dp*powm1(pow6(x%val))
      q28 = q26*q27
      q29 = q0*x%d1val2
      q30 = 2.0_dp*q29
      q31(1:5) = x%d1Array(1:5)*x%d1val1_d1val2
      q32(1:5) = x%d1val1*x%d1val2_d1Array(1:5)
      q33(1:5) = x%d1val1_d1Array(1:5)*x%d1val2
      unary%val = pow(x%val, y_dp)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q4*q6
      unary%d1val1_d1val2 = q6*(q7 + q8*y_dp - q8)
      unary%d1val1_d1Array(1:5) = q6*(q10*y_dp - q10 + q9)
      unary%d2val2 = q15*q5
      unary%d1val2_d1Array(1:5) = q6*(q16 + q17*y_dp - q17)
      unary%d3val1 = q21*(3.0_dp*q19*q2 + q18*x%d3val1 + q20*pow3(x%d1val1))
      unary%d2val1_d1val2 = q28*(q18*x%d2val1_d1val2 + q22*q7 + q24*x%d1val2 + q25*x%d1val2)
      unary%d2val1_d1Array(1:5) = q28*(q18*x%d2val1_d1Array(1:5) + q22*q9 + q24*x%d1Array(1:5) + q25*x%d1Array(1:5))
      unary%d1val1_d2val2 = q21*(q18*x%d1val1_d2val2 + q30*q7 - x%d1val1*(-2.0_dp*q12 + 2.0_dp*q13 + q11 - q15))
      unary%d1val1_d1val2_d1Array(1:5) = q27*(q20*q26*q8*x%d1Array(1:5) + x%d1val1_d1val2_d1Array(1:5)*pow(x%val, 5.0_dp + y_dp) + (q31*y_dp - q31 + q32*y_dp - q32 + q33*y_dp - q33)*pow(x%val, 4.0_dp + y_dp))
      unary%d3val2 = q21*(3.0_dp*q11*q29 + q18*x%d3val2 + q20*pow3(x%d1val2))
      unary%d2val2_d1Array(1:5) = q28*(q0*q14*x%d1Array(1:5) + q12*q23*x%d1Array(1:5) + q16*q30 + q18*x%d2val2_d1Array(1:5))
   end function pow_self_int
   
   function pow_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q10
      real(dp) :: q9(1:5)
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
      q4 = q1*x%d1val1
      q5 = q1*pow2(x%d1val2) + x%d2val2
      q6 = q1*x%d1val2
      q7 = pow2(q1)
      q8 = 2.0_dp*q4
      q9(1:5) = q1*x%d1Array(1:5)
      q10 = 2.0_dp*q6
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1Array(1:5) = q2*x%d1Array(1:5)
      unary%d2val1 = q2*q3
      unary%d1val1_d1val2 = q2*(q4*x%d1val2 + x%d1val1_d1val2)
      unary%d1val1_d1Array(1:5) = q2*(q4*x%d1Array(1:5) + x%d1val1_d1Array(1:5))
      unary%d2val2 = q2*q5
      unary%d1val2_d1Array(1:5) = q2*(q6*x%d1Array(1:5) + x%d1val2_d1Array(1:5))
      unary%d3val1 = q2*(3.0_dp*q4*x%d2val1 + q7*pow3(x%d1val1) + x%d3val1)
      unary%d2val1_d1val2 = q2*(q3*q6 + q8*x%d1val1_d1val2 + x%d2val1_d1val2)
      unary%d2val1_d1Array(1:5) = q2*(q3*q9 + q8*x%d1val1_d1Array(1:5) + x%d2val1_d1Array(1:5))
      unary%d1val1_d2val2 = q2*(q10*x%d1val1_d1val2 + q4*q5 + x%d1val1_d2val2)
      unary%d1val1_d1val2_d1Array(1:5) = q2*(q4*x%d1val2_d1Array(1:5) + q6*x%d1val1_d1Array(1:5) + q7*x%d1Array(1:5)*x%d1val1*x%d1val2 + q9*x%d1val1_d1val2 + x%d1val1_d1val2_d1Array(1:5))
      unary%d3val2 = q2*(3.0_dp*q6*x%d2val2 + q7*pow3(x%d1val2) + x%d3val2)
      unary%d2val2_d1Array(1:5) = q2*(q10*x%d1val2_d1Array(1:5) + q5*q9 + x%d2val2_d1Array(1:5))
   end function pow_int_self
   
   function max_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: binary
      real(dp) :: q1
      real(dp) :: q0
      q0 = Heaviside(x%val - y%val)
      q1 = Heaviside(-x%val + y%val)
      binary%val = Max(x%val, y%val)
      binary%d1val1 = q0*x%d1val1 + q1*y%d1val1
      binary%d1val2 = q0*x%d1val2 + q1*y%d1val2
      binary%d1Array(1:5) = q0*x%d1Array(1:5) + q1*y%d1Array(1:5)
      binary%d2val1 = q0*x%d2val1 + q1*y%d2val1
      binary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q1*y%d1val1_d1val2
      binary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5) + q1*y%d1val1_d1Array(1:5)
      binary%d2val2 = q0*x%d2val2 + q1*y%d2val2
      binary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5) + q1*y%d1val2_d1Array(1:5)
      binary%d3val1 = q0*x%d3val1 + q1*y%d3val1
      binary%d2val1_d1val2 = q0*x%d2val1_d1val2 + q1*y%d2val1_d1val2
      binary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5) + q1*y%d2val1_d1Array(1:5)
      binary%d1val1_d2val2 = q0*x%d1val1_d2val2 + q1*y%d1val1_d2val2
      binary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5) + q1*y%d1val1_d1val2_d1Array(1:5)
      binary%d3val2 = q0*x%d3val2 + q1*y%d3val2
      binary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5) + q1*y%d2val2_d1Array(1:5)
   end function max_self
   
   function max_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q0
      q0 = Heaviside(x%val - y)
      unary%val = Max(x%val, y)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function max_self_real
   
   function max_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q0
      q0 = Heaviside(x%val - z)
      unary%val = Max(x%val, z)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function max_real_self
   
   function max_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = y
      q0 = Heaviside(x%val - y_dp)
      unary%val = Max(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function max_self_int
   
   function max_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = z
      q0 = Heaviside(x%val - y_dp)
      unary%val = Max(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function max_int_self
   
   function min_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: binary
      real(dp) :: q1
      real(dp) :: q0
      q0 = Heaviside(-x%val + y%val)
      q1 = Heaviside(x%val - y%val)
      binary%val = Min(x%val, y%val)
      binary%d1val1 = q0*x%d1val1 + q1*y%d1val1
      binary%d1val2 = q0*x%d1val2 + q1*y%d1val2
      binary%d1Array(1:5) = q0*x%d1Array(1:5) + q1*y%d1Array(1:5)
      binary%d2val1 = q0*x%d2val1 + q1*y%d2val1
      binary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q1*y%d1val1_d1val2
      binary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5) + q1*y%d1val1_d1Array(1:5)
      binary%d2val2 = q0*x%d2val2 + q1*y%d2val2
      binary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5) + q1*y%d1val2_d1Array(1:5)
      binary%d3val1 = q0*x%d3val1 + q1*y%d3val1
      binary%d2val1_d1val2 = q0*x%d2val1_d1val2 + q1*y%d2val1_d1val2
      binary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5) + q1*y%d2val1_d1Array(1:5)
      binary%d1val1_d2val2 = q0*x%d1val1_d2val2 + q1*y%d1val1_d2val2
      binary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5) + q1*y%d1val1_d1val2_d1Array(1:5)
      binary%d3val2 = q0*x%d3val2 + q1*y%d3val2
      binary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5) + q1*y%d2val2_d1Array(1:5)
   end function min_self
   
   function min_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q0
      q0 = Heaviside(-x%val + y)
      unary%val = Min(x%val, y)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function min_self_real
   
   function min_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q0
      q0 = Heaviside(-x%val + z)
      unary%val = Min(x%val, z)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function min_real_self
   
   function min_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = y
      q0 = Heaviside(-x%val + y_dp)
      unary%val = Min(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function min_self_int
   
   function min_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = z
      q0 = Heaviside(-x%val + y_dp)
      unary%val = Min(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1Array(1:5) = q0*x%d1Array(1:5)
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q0*x%d1val1_d1Array(1:5)
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1Array(1:5) = q0*x%d1val2_d1Array(1:5)
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q0*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q0*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1Array(1:5) = q0*x%d2val2_d1Array(1:5)
   end function min_int_self
   
   function dim_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: binary
      real(dp) :: q9(1:5)
      real(dp) :: q8(1:5)
      real(dp) :: q7(1:5)
      real(dp) :: q6
      real(dp) :: q5(1:5)
      real(dp) :: q4(1:5)
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val - y%val
      q1 = sgn(q0)
      q2 = 0.5_dp*q1
      q3 = -0.5_dp*y%d1val1_d1val2 + 0.5_dp*x%d1val1_d1val2
      q4(1:5) = -0.5_dp*y%d1val1_d1Array(1:5) + 0.5_dp*x%d1val1_d1Array(1:5)
      q5(1:5) = -0.5_dp*y%d1val2_d1Array(1:5) + 0.5_dp*x%d1val2_d1Array(1:5)
      q6 = -0.5_dp*y%d2val1_d1val2 + 0.5_dp*x%d2val1_d1val2
      q7(1:5) = -0.5_dp*y%d2val1_d1Array(1:5) + 0.5_dp*x%d2val1_d1Array(1:5)
      q8(1:5) = -0.5_dp*y%d1val1_d1val2_d1Array(1:5) + 0.5_dp*x%d1val1_d1val2_d1Array(1:5)
      q9(1:5) = -0.5_dp*y%d2val2_d1Array(1:5) + 0.5_dp*x%d2val2_d1Array(1:5)
      binary%val = -0.5_dp*y%val + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      binary%d1val1 = -0.5_dp*y%d1val1 + 0.5_dp*x%d1val1 + q2*(x%d1val1 - y%d1val1)
      binary%d1val2 = -0.5_dp*y%d1val2 + 0.5_dp*x%d1val2 + q2*(x%d1val2 - y%d1val2)
      binary%d1Array(1:5) = -0.5_dp*y%d1Array(1:5) + 0.5_dp*x%d1Array(1:5) + q2*(x%d1Array(1:5) - y%d1Array(1:5))
      binary%d2val1 = -0.5_dp*y%d2val1 + 0.5_dp*x%d2val1 + q2*(x%d2val1 - y%d2val1)
      binary%d1val1_d1val2 = q1*q3 + q3
      binary%d1val1_d1Array(1:5) = q1*q4 + q4
      binary%d2val2 = -0.5_dp*y%d2val2 + 0.5_dp*x%d2val2 + q2*(x%d2val2 - y%d2val2)
      binary%d1val2_d1Array(1:5) = q1*q5 + q5
      binary%d3val1 = -0.5_dp*y%d3val1 + 0.5_dp*x%d3val1 + q2*(x%d3val1 - y%d3val1)
      binary%d2val1_d1val2 = q1*q6 + q6
      binary%d2val1_d1Array(1:5) = q1*q7 + q7
      binary%d1val1_d2val2 = -0.5_dp*y%d1val1_d2val2 + 0.5_dp*x%d1val1_d2val2 + q2*(x%d1val1_d2val2 - y%d1val1_d2val2)
      binary%d1val1_d1val2_d1Array(1:5) = q1*q8 + q8
      binary%d3val2 = -0.5_dp*y%d3val2 + 0.5_dp*x%d3val2 + q2*(x%d3val2 - y%d3val2)
      binary%d2val2_d1Array(1:5) = q1*q9 + q9
   end function dim_self
   
   function dim_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val - y
      q1 = 0.5_dp*sgn(q0) + 0.5_dp
      unary%val = -0.5_dp*y + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q1*x%d1val1_d1Array(1:5)
      unary%d2val2 = q1*x%d2val2
      unary%d1val2_d1Array(1:5) = q1*x%d1val2_d1Array(1:5)
      unary%d3val1 = q1*x%d3val1
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q1*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q1*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q1*x%d3val2
      unary%d2val2_d1Array(1:5) = q1*x%d2val2_d1Array(1:5)
   end function dim_self_real
   
   function dim_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val - z
      q1 = -0.5_dp + 0.5_dp*sgn(q0)
      unary%val = -0.5_dp*x%val + 0.5_dp*z + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q1*x%d1val1_d1Array(1:5)
      unary%d2val2 = q1*x%d2val2
      unary%d1val2_d1Array(1:5) = q1*x%d1val2_d1Array(1:5)
      unary%d3val1 = q1*x%d3val1
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q1*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q1*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q1*x%d3val2
      unary%d2val2_d1Array(1:5) = q1*x%d2val2_d1Array(1:5)
   end function dim_real_self
   
   function dim_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q1
      real(dp) :: q0
      y_dp = y
      q0 = x%val - y_dp
      q1 = 0.5_dp*sgn(q0) + 0.5_dp
      unary%val = -0.5_dp*y_dp + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q1*x%d1val1_d1Array(1:5)
      unary%d2val2 = q1*x%d2val2
      unary%d1val2_d1Array(1:5) = q1*x%d1val2_d1Array(1:5)
      unary%d3val1 = q1*x%d3val1
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q1*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q1*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q1*x%d3val2
      unary%d2val2_d1Array(1:5) = q1*x%d2val2_d1Array(1:5)
   end function dim_self_int
   
   function dim_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_array), intent(in) :: x
      type(auto_diff_real_2var_order3_array) :: unary
      real(dp) :: y_dp
      real(dp) :: q1
      real(dp) :: q0
      y_dp = z
      q0 = x%val - y_dp
      q1 = -0.5_dp + 0.5_dp*sgn(q0)
      unary%val = -0.5_dp*x%val + 0.5_dp*y_dp + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1Array(1:5) = q1*x%d1Array(1:5)
      unary%d2val1 = q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2
      unary%d1val1_d1Array(1:5) = q1*x%d1val1_d1Array(1:5)
      unary%d2val2 = q1*x%d2val2
      unary%d1val2_d1Array(1:5) = q1*x%d1val2_d1Array(1:5)
      unary%d3val1 = q1*x%d3val1
      unary%d2val1_d1val2 = q1*x%d2val1_d1val2
      unary%d2val1_d1Array(1:5) = q1*x%d2val1_d1Array(1:5)
      unary%d1val1_d2val2 = q1*x%d1val1_d2val2
      unary%d1val1_d1val2_d1Array(1:5) = q1*x%d1val1_d1val2_d1Array(1:5)
      unary%d3val2 = q1*x%d3val2
      unary%d2val2_d1Array(1:5) = q1*x%d2val2_d1Array(1:5)
   end function dim_int_self
   
   function differentiate_auto_diff_real_2var_order3_array_1(this) result(derivative)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      type(auto_diff_real_2var_order3_array) :: derivative
      derivative%val = this%d1val1
      derivative%d1val1 = this%d2val1
      derivative%d1val2 = this%d1val1_d1val2
      derivative%d1Array = this%d1val1_d1Array
      derivative%d2val1 = this%d3val1
      derivative%d1val1_d1val2 = this%d2val1_d1val2
      derivative%d1val1_d1Array = this%d2val1_d1Array
      derivative%d2val2 = this%d1val1_d2val2
      derivative%d1val2_d1Array = this%d1val1_d1val2_d1Array
      derivative%d3val1 = 0_dp
      derivative%d2val1_d1val2 = 0_dp
      derivative%d2val1_d1Array = 0_dp
      derivative%d1val1_d2val2 = 0_dp
      derivative%d1val1_d1val2_d1Array = 0_dp
      derivative%d3val2 = 0_dp
      derivative%d2val2_d1Array = 0_dp
   end function differentiate_auto_diff_real_2var_order3_array_1
   
   function differentiate_auto_diff_real_2var_order3_array_2(this) result(derivative)
      type(auto_diff_real_2var_order3_array), intent(in) :: this
      type(auto_diff_real_2var_order3_array) :: derivative
      derivative%val = this%d1val2
      derivative%d1val1 = this%d1val1_d1val2
      derivative%d1val2 = this%d2val2
      derivative%d1Array = this%d1val2_d1Array
      derivative%d2val1 = this%d2val1_d1val2
      derivative%d1val1_d1val2 = this%d1val1_d2val2
      derivative%d1val1_d1Array = this%d1val1_d1val2_d1Array
      derivative%d2val2 = this%d3val2
      derivative%d1val2_d1Array = this%d2val2_d1Array
      derivative%d3val1 = 0_dp
      derivative%d2val1_d1val2 = 0_dp
      derivative%d2val1_d1Array = 0_dp
      derivative%d1val1_d2val2 = 0_dp
      derivative%d1val1_d1val2_d1Array = 0_dp
      derivative%d3val2 = 0_dp
      derivative%d2val2_d1Array = 0_dp
   end function differentiate_auto_diff_real_2var_order3_array_2
   
end module auto_diff_real_2var_order3_array_module