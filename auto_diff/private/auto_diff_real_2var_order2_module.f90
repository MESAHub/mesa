module auto_diff_real_2var_order2_module
      use const_def
      use utils_lib
      use support_functions
      use math_lib
   
      implicit none
      private
   public :: auto_diff_real_2var_order2, &
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
   type :: auto_diff_real_2var_order2
      real(dp) :: val
      real(dp) :: d1val1
      real(dp) :: d1val2
      real(dp) :: d2val1
      real(dp) :: d1val1_d1val2
      real(dp) :: d2val2
   end type auto_diff_real_2var_order2
   
   interface assignment(=)
      module procedure assign_from_self
      module procedure assign_from_real_dp
      module procedure assign_from_int
   end interface assignment(=)
   
   interface operator(.eq.)
      module procedure equal_self
      module procedure equal_auto_diff_real_2var_order2_real_dp
      module procedure equal_real_dp_auto_diff_real_2var_order2
      module procedure equal_auto_diff_real_2var_order2_int
      module procedure equal_int_auto_diff_real_2var_order2
   end interface operator(.eq.)
   
   interface operator(.ne.)
      module procedure neq_self
      module procedure neq_auto_diff_real_2var_order2_real_dp
      module procedure neq_real_dp_auto_diff_real_2var_order2
      module procedure neq_auto_diff_real_2var_order2_int
      module procedure neq_int_auto_diff_real_2var_order2
   end interface operator(.ne.)
   
   interface operator(.gt.)
      module procedure greater_self
      module procedure greater_auto_diff_real_2var_order2_real_dp
      module procedure greater_real_dp_auto_diff_real_2var_order2
      module procedure greater_auto_diff_real_2var_order2_int
      module procedure greater_int_auto_diff_real_2var_order2
   end interface operator(.gt.)
   
   interface operator(.lt.)
      module procedure less_self
      module procedure less_auto_diff_real_2var_order2_real_dp
      module procedure less_real_dp_auto_diff_real_2var_order2
      module procedure less_auto_diff_real_2var_order2_int
      module procedure less_int_auto_diff_real_2var_order2
   end interface operator(.lt.)
   
   interface operator(.le.)
      module procedure leq_self
      module procedure leq_auto_diff_real_2var_order2_real_dp
      module procedure leq_real_dp_auto_diff_real_2var_order2
      module procedure leq_auto_diff_real_2var_order2_int
      module procedure leq_int_auto_diff_real_2var_order2
   end interface operator(.le.)
   
   interface operator(.ge.)
      module procedure geq_self
      module procedure geq_auto_diff_real_2var_order2_real_dp
      module procedure geq_real_dp_auto_diff_real_2var_order2
      module procedure geq_auto_diff_real_2var_order2_int
      module procedure geq_int_auto_diff_real_2var_order2
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
      module procedure differentiate_auto_diff_real_2var_order2_1
   end interface differentiate_1
   
   interface differentiate_2
      module procedure differentiate_auto_diff_real_2var_order2_2
   end interface differentiate_2
   
   contains

   subroutine assign_from_self(this, other)
      type(auto_diff_real_2var_order2), intent(out) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      this%val = other%val
      this%d1val1 = other%d1val1
      this%d1val2 = other%d1val2
      this%d2val1 = other%d2val1
      this%d1val1_d1val2 = other%d1val1_d1val2
      this%d2val2 = other%d2val2
   end subroutine assign_from_self
   
   subroutine assign_from_real_dp(this, other)
      type(auto_diff_real_2var_order2), intent(out) :: this
      real(dp), intent(in) :: other
      this%val = other
      this%d1val1 = 0_dp
      this%d1val2 = 0_dp
      this%d2val1 = 0_dp
      this%d1val1_d1val2 = 0_dp
      this%d2val2 = 0_dp
   end subroutine assign_from_real_dp
   
   subroutine assign_from_int(this, other)
      type(auto_diff_real_2var_order2), intent(out) :: this
      integer, intent(in) :: other
      this%val = other
      this%d1val1 = 0_dp
      this%d1val2 = 0_dp
      this%d2val1 = 0_dp
      this%d1val1_d1val2 = 0_dp
      this%d2val2 = 0_dp
   end subroutine assign_from_int
   
   function equal_self(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this%val .eq. other%val)
   end function equal_self
   
   function equal_auto_diff_real_2var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .eq. other)
   end function equal_auto_diff_real_2var_order2_real_dp
   
   function equal_real_dp_auto_diff_real_2var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .eq. other%val)
   end function equal_real_dp_auto_diff_real_2var_order2
   
   function equal_auto_diff_real_2var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .eq. other)
   end function equal_auto_diff_real_2var_order2_int
   
   function equal_int_auto_diff_real_2var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .eq. other%val)
   end function equal_int_auto_diff_real_2var_order2
   
   function neq_self(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this%val .ne. other%val)
   end function neq_self
   
   function neq_auto_diff_real_2var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .ne. other)
   end function neq_auto_diff_real_2var_order2_real_dp
   
   function neq_real_dp_auto_diff_real_2var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .ne. other%val)
   end function neq_real_dp_auto_diff_real_2var_order2
   
   function neq_auto_diff_real_2var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .ne. other)
   end function neq_auto_diff_real_2var_order2_int
   
   function neq_int_auto_diff_real_2var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .ne. other%val)
   end function neq_int_auto_diff_real_2var_order2
   
   function greater_self(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this%val .gt. other%val)
   end function greater_self
   
   function greater_auto_diff_real_2var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .gt. other)
   end function greater_auto_diff_real_2var_order2_real_dp
   
   function greater_real_dp_auto_diff_real_2var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .gt. other%val)
   end function greater_real_dp_auto_diff_real_2var_order2
   
   function greater_auto_diff_real_2var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .gt. other)
   end function greater_auto_diff_real_2var_order2_int
   
   function greater_int_auto_diff_real_2var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .gt. other%val)
   end function greater_int_auto_diff_real_2var_order2
   
   function less_self(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this%val .lt. other%val)
   end function less_self
   
   function less_auto_diff_real_2var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .lt. other)
   end function less_auto_diff_real_2var_order2_real_dp
   
   function less_real_dp_auto_diff_real_2var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .lt. other%val)
   end function less_real_dp_auto_diff_real_2var_order2
   
   function less_auto_diff_real_2var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .lt. other)
   end function less_auto_diff_real_2var_order2_int
   
   function less_int_auto_diff_real_2var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .lt. other%val)
   end function less_int_auto_diff_real_2var_order2
   
   function leq_self(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this%val .le. other%val)
   end function leq_self
   
   function leq_auto_diff_real_2var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .le. other)
   end function leq_auto_diff_real_2var_order2_real_dp
   
   function leq_real_dp_auto_diff_real_2var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .le. other%val)
   end function leq_real_dp_auto_diff_real_2var_order2
   
   function leq_auto_diff_real_2var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .le. other)
   end function leq_auto_diff_real_2var_order2_int
   
   function leq_int_auto_diff_real_2var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .le. other%val)
   end function leq_int_auto_diff_real_2var_order2
   
   function geq_self(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this%val .ge. other%val)
   end function geq_self
   
   function geq_auto_diff_real_2var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .ge. other)
   end function geq_auto_diff_real_2var_order2_real_dp
   
   function geq_real_dp_auto_diff_real_2var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .ge. other%val)
   end function geq_real_dp_auto_diff_real_2var_order2
   
   function geq_auto_diff_real_2var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .ge. other)
   end function geq_auto_diff_real_2var_order2_int
   
   function geq_int_auto_diff_real_2var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order2), intent(in) :: other
      logical :: z
      z = (this .ge. other%val)
   end function geq_int_auto_diff_real_2var_order2
   
   function make_unary_operator(x, z_val, z_d1x, z_d2x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      real(dp), intent(in) :: z_val
      real(dp), intent(in) :: z_d1x
      real(dp), intent(in) :: z_d2x
      type(auto_diff_real_2var_order2) :: unary
      unary%val = z_val
      unary%d1val1 = x%d1val1*z_d1x
      unary%d1val2 = x%d1val2*z_d1x
      unary%d2val1 = x%d2val1*z_d1x + z_d2x*pow2(x%d1val1)
      unary%d1val1_d1val2 = x%d1val1*x%d1val2*z_d2x + x%d1val1_d1val2*z_d1x
      unary%d2val2 = x%d2val2*z_d1x + z_d2x*pow2(x%d1val2)
   end function make_unary_operator
   
   function make_binary_operator(x, y, z_val, z_d1x, z_d1y, z_d2x, z_d1x_d1y, z_d2y) result(binary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2), intent(in) :: y
      real(dp), intent(in) :: z_val
      real(dp), intent(in) :: z_d1x
      real(dp), intent(in) :: z_d1y
      real(dp), intent(in) :: z_d2x
      real(dp), intent(in) :: z_d1x_d1y
      real(dp), intent(in) :: z_d2y
      type(auto_diff_real_2var_order2) :: binary
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%d1val1*z_d1x_d1y
      q1 = x%d1val2*z_d1x_d1y
      binary%val = z_val
      binary%d1val1 = x%d1val1*z_d1x + y%d1val1*z_d1y
      binary%d1val2 = x%d1val2*z_d1x + y%d1val2*z_d1y
      binary%d2val1 = 2.0_dp*q0*y%d1val1 + x%d2val1*z_d1x + y%d2val1*z_d1y + z_d2x*pow2(x%d1val1) + z_d2y*pow2(y%d1val1)
      binary%d1val1_d1val2 = q0*y%d1val2 + q1*y%d1val1 + x%d1val1*x%d1val2*z_d2x + x%d1val1_d1val2*z_d1x + y%d1val1*y%d1val2*z_d2y + y%d1val1_d1val2*z_d1y
      binary%d2val2 = 2.0_dp*q1*y%d1val2 + x%d2val2*z_d1x + y%d2val2*z_d1y + z_d2x*pow2(x%d1val2) + z_d2y*pow2(y%d1val2)
   end function make_binary_operator
   
   function unary_minus_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      unary%val = -x%val
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d2val2 = -x%d2val2
   end function unary_minus_self
   
   function exp_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q0
      q0 = exp(x%val)
      unary%val = q0
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*(x%d2val1 + pow2(x%d1val1))
      unary%d1val1_d1val2 = q0*(x%d1val1*x%d1val2 + x%d1val1_d1val2)
      unary%d2val2 = q0*(x%d2val2 + pow2(x%d1val2))
   end function exp_self
   
   function expm1_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q0
      q0 = exp(x%val)
      unary%val = expm1(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*(x%d2val1 + pow2(x%d1val1))
      unary%d1val1_d1val2 = q0*(x%d1val1*x%d1val2 + x%d1val1_d1val2)
      unary%d2val2 = q0*(x%d2val2 + pow2(x%d1val2))
   end function expm1_self
   
   function exp10_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow(10.0_dp, x%val)
      q1 = ln10
      q2 = q0*q1
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*(q1*pow2(x%d1val1) + x%d2val1)
      unary%d1val1_d1val2 = q2*(q1*x%d1val1*x%d1val2 + x%d1val1_d1val2)
      unary%d2val2 = q2*(q1*pow2(x%d1val2) + x%d2val2)
   end function exp10_self
   
   function powm1_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(pow2(x%val))
      q1 = powm1(pow3(x%val))
      unary%val = powm1(x%val)
      unary%d1val1 = -q0*x%d1val1
      unary%d1val2 = -q0*x%d1val2
      unary%d2val1 = q1*(2.0_dp*pow2(x%d1val1) - x%d2val1*x%val)
      unary%d1val1_d1val2 = q1*(2.0_dp*x%d1val1*x%d1val2 - x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(2.0_dp*pow2(x%d1val2) - x%d2val2*x%val)
   end function powm1_self
   
   function log_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(x%val)
      q1 = powm1(pow2(x%val))
      unary%val = log(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q1*(x%d2val1*x%val - pow2(x%d1val1))
      unary%d1val1_d1val2 = q1*(-x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(x%d2val2*x%val - pow2(x%d1val2))
   end function log_self
   
   function log1p_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val + 1
      q1 = powm1(q0)
      q2 = powm1(pow2(q0))
      unary%val = log1p(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q2*(q0*x%d2val1 - pow2(x%d1val1))
      unary%d1val1_d1val2 = q2*(q0*x%d1val1_d1val2 - x%d1val1*x%d1val2)
      unary%d2val2 = q2*(q0*x%d2val2 - pow2(x%d1val2))
   end function log1p_self
   
   function safe_log_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(x%val)
      q1 = powm1(pow2(x%val))
      unary%val = safe_log(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q1*(x%d2val1*x%val - pow2(x%d1val1))
      unary%d1val1_d1val2 = q1*(-x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(x%d2val2*x%val - pow2(x%d1val2))
   end function safe_log_self
   
   function log10_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(ln10)
      q1 = q0*powm1(x%val)
      q2 = q0*powm1(pow2(x%val))
      unary%val = q0*log(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q2*(x%d2val1*x%val - pow2(x%d1val1))
      unary%d1val1_d1val2 = q2*(-x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q2*(x%d2val2*x%val - pow2(x%d1val2))
   end function log10_self
   
   function safe_log10_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(ln10)
      q1 = q0*powm1(x%val)
      q2 = q0*powm1(pow2(x%val))
      unary%val = q0*safe_log(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q2*(x%d2val1*x%val - pow2(x%d1val1))
      unary%d1val1_d1val2 = q2*(-x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q2*(x%d2val2*x%val - pow2(x%d1val2))
   end function safe_log10_self
   
   function log2_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(log(2.0_dp))
      q1 = q0*powm1(x%val)
      q2 = q0*powm1(pow2(x%val))
      unary%val = q0*log(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q2*(x%d2val1*x%val - pow2(x%d1val1))
      unary%d1val1_d1val2 = q2*(-x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q2*(x%d2val2*x%val - pow2(x%d1val2))
   end function log2_self
   
   function sin_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = sin(x%val)
      q1 = cos(x%val)
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = -q0*pow2(x%d1val1) + q1*x%d2val1
      unary%d1val1_d1val2 = -q0*x%d1val1*x%d1val2 + q1*x%d1val1_d1val2
      unary%d2val2 = -q0*pow2(x%d1val2) + q1*x%d2val2
   end function sin_self
   
   function cos_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = cos(x%val)
      q1 = sin(x%val)
      unary%val = q0
      unary%d1val1 = -q1*x%d1val1
      unary%d1val2 = -q1*x%d1val2
      unary%d2val1 = -q0*pow2(x%d1val1) - q1*x%d2val1
      unary%d1val1_d1val2 = -q0*x%d1val1*x%d1val2 - q1*x%d1val1_d1val2
      unary%d2val2 = -q0*pow2(x%d1val2) - q1*x%d2val2
   end function cos_self
   
   function tan_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = tan(x%val)
      q1 = powm1(pow2(cos(x%val)))
      q2 = 2.0_dp*q0
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q1*(q2*pow2(x%d1val1) + x%d2val1)
      unary%d1val1_d1val2 = q1*(q2*x%d1val1*x%d1val2 + x%d1val1_d1val2)
      unary%d2val2 = q1*(q2*pow2(x%d1val2) + x%d2val2)
   end function tan_self
   
   function sinpi_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pi*x%val
      q1 = sin(q0)
      q2 = cos(q0)
      q3 = pi*q2
      q4 = pi*q1
      unary%val = q1
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d2val1 = pi*(q2*x%d2val1 - q4*pow2(x%d1val1))
      unary%d1val1_d1val2 = pi*(q2*x%d1val1_d1val2 - q4*x%d1val1*x%d1val2)
      unary%d2val2 = pi*(q2*x%d2val2 - q4*pow2(x%d1val2))
   end function sinpi_self
   
   function cospi_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pi*x%val
      q1 = cos(q0)
      q2 = sin(q0)
      q3 = pi*q2
      q4 = pi*q1
      unary%val = q1
      unary%d1val1 = -q3*x%d1val1
      unary%d1val2 = -q3*x%d1val2
      unary%d2val1 = -pi*(q2*x%d2val1 + q4*pow2(x%d1val1))
      unary%d1val1_d1val2 = -pi*(q2*x%d1val1_d1val2 + q4*x%d1val1*x%d1val2)
      unary%d2val2 = -pi*(q2*x%d2val2 + q4*pow2(x%d1val2))
   end function cospi_self
   
   function tanpi_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pi*x%val
      q1 = tan(q0)
      q2 = pi*powm1(pow2(cos(q0)))
      q3 = 2.0_dp*pi*q1
      unary%val = q1
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*(q3*pow2(x%d1val1) + x%d2val1)
      unary%d1val1_d1val2 = q2*(q3*x%d1val1*x%d1val2 + x%d1val1_d1val2)
      unary%d2val2 = q2*(q3*pow2(x%d1val2) + x%d2val2)
   end function tanpi_self
   
   function sinh_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = sinh(x%val)
      q1 = cosh(x%val)
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q0*pow2(x%d1val1) + q1*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1*x%d1val2 + q1*x%d1val1_d1val2
      unary%d2val2 = q0*pow2(x%d1val2) + q1*x%d2val2
   end function sinh_self
   
   function cosh_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = cosh(x%val)
      q1 = sinh(x%val)
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q0*pow2(x%d1val1) + q1*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1*x%d1val2 + q1*x%d1val1_d1val2
      unary%d2val2 = q0*pow2(x%d1val2) + q1*x%d2val2
   end function cosh_self
   
   function tanh_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = tanh(x%val)
      q1 = powm1(pow2(cosh(x%val)))
      q2 = pow2(q0) - 1
      q3 = 2.0_dp*q0
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q2*(q3*pow2(x%d1val1) - x%d2val1)
      unary%d1val1_d1val2 = q2*(q3*x%d1val1*x%d1val2 - x%d1val1_d1val2)
      unary%d2val2 = q2*(q3*pow2(x%d1val2) - x%d2val2)
   end function tanh_self
   
   function asin_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = 1 - q0
      q2 = powm1(sqrt(q1))
      q3 = powm1(pow3(sqrt(q1)))
      q4 = q0 - 1
      unary%val = asin(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q3*(-q4*x%d2val1 + x%val*pow2(x%d1val1))
      unary%d1val1_d1val2 = q3*(q1*x%d1val1_d1val2 + x%d1val1*x%d1val2*x%val)
      unary%d2val2 = q3*(-q4*x%d2val2 + x%val*pow2(x%d1val2))
   end function asin_self
   
   function acos_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val)
      q1 = 1 - q0
      q2 = powm1(sqrt(q1))
      q3 = powm1(pow3(sqrt(q1)))
      q4 = q0 - 1
      unary%val = acos(x%val)
      unary%d1val1 = -q2*x%d1val1
      unary%d1val2 = -q2*x%d1val2
      unary%d2val1 = q3*(q4*x%d2val1 - x%val*pow2(x%d1val1))
      unary%d1val1_d1val2 = -q3*(-q4*x%d1val1_d1val2 + x%d1val1*x%d1val2*x%val)
      unary%d2val2 = q3*(q4*x%d2val2 - x%val*pow2(x%d1val2))
   end function acos_self
   
   function atan_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val) + 1
      q1 = powm1(q0)
      q2 = powm1(pow2(q0))
      q3 = 2.0_dp*x%val
      unary%val = atan(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q2*(q0*x%d2val1 - q3*pow2(x%d1val1))
      unary%d1val1_d1val2 = q2*(q0*x%d1val1_d1val2 - q3*x%d1val1*x%d1val2)
      unary%d2val2 = q2*(q0*x%d2val2 - q3*pow2(x%d1val2))
   end function atan_self
   
   function asinpi_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
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
      q4 = q1 - 1
      q5 = q0*powm1(pow3(sqrt(q2)))
      unary%val = q0*asin(x%val)
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d2val1 = q5*(-q4*x%d2val1 + x%val*pow2(x%d1val1))
      unary%d1val1_d1val2 = q5*(q2*x%d1val1_d1val2 + x%d1val1*x%d1val2*x%val)
      unary%d2val2 = q5*(-q4*x%d2val2 + x%val*pow2(x%d1val2))
   end function asinpi_self
   
   function acospi_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
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
      q4 = q1 - 1
      q5 = q0*powm1(pow3(sqrt(q2)))
      unary%val = q0*acos(x%val)
      unary%d1val1 = -q3*x%d1val1
      unary%d1val2 = -q3*x%d1val2
      unary%d2val1 = q5*(q4*x%d2val1 - x%val*pow2(x%d1val1))
      unary%d1val1_d1val2 = -q5*(-q4*x%d1val1_d1val2 + x%d1val1*x%d1val2*x%val)
      unary%d2val2 = q5*(q4*x%d2val2 - x%val*pow2(x%d1val2))
   end function acospi_self
   
   function atanpi_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = powm1(pi)
      q1 = pow2(x%val) + 1
      q2 = q0*powm1(q1)
      q3 = 2.0_dp*x%val
      q4 = q0*powm1(pow2(q1))
      unary%val = q0*atan(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q4*(q1*x%d2val1 - q3*pow2(x%d1val1))
      unary%d1val1_d1val2 = q4*(q1*x%d1val1_d1val2 - q3*x%d1val1*x%d1val2)
      unary%d2val2 = q4*(q1*x%d2val2 - q3*pow2(x%d1val2))
   end function atanpi_self
   
   function asinh_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val) + 1
      q1 = powm1(sqrt(q0))
      q2 = powm1(pow3(sqrt(q0)))
      unary%val = asinh(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q2*(q0*x%d2val1 - x%val*pow2(x%d1val1))
      unary%d1val1_d1val2 = q2*(q0*x%d1val1_d1val2 - x%d1val1*x%d1val2*x%val)
      unary%d2val2 = q2*(q0*x%d2val2 - x%val*pow2(x%d1val2))
   end function asinh_self
   
   function acosh_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val) - 1
      q1 = powm1(sqrt(q0))
      q2 = powm1(pow3(sqrt(q0)))
      unary%val = acosh(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q2*(q0*x%d2val1 - x%val*pow2(x%d1val1))
      unary%d1val1_d1val2 = q2*(q0*x%d1val1_d1val2 - x%d1val1*x%d1val2*x%val)
      unary%d2val2 = q2*(q0*x%d2val2 - x%val*pow2(x%d1val2))
   end function acosh_self
   
   function atanh_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow2(x%val) - 1
      q1 = powm1(q0)
      q2 = powm1(pow2(q0))
      q3 = 2.0_dp*x%val
      unary%val = atanh(x%val)
      unary%d1val1 = -q1*x%d1val1
      unary%d1val2 = -q1*x%d1val2
      unary%d2val1 = q2*(-q0*x%d2val1 + q3*pow2(x%d1val1))
      unary%d1val1_d1val2 = q2*(-q0*x%d1val1_d1val2 + q3*x%d1val1*x%d1val2)
      unary%d2val2 = q2*(-q0*x%d2val2 + q3*pow2(x%d1val2))
   end function atanh_self
   
   function sqrt_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = sqrt(x%val)
      q1 = 0.5_dp*powm1(q0)
      q2 = 2.0_dp*x%val
      q3 = 0.25_dp*powm1(pow3(sqrt(x%val)))
      unary%val = q0
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q3*(q2*x%d2val1 - pow2(x%d1val1))
      unary%d1val1_d1val2 = q3*(q2*x%d1val1_d1val2 - x%d1val1*x%d1val2)
      unary%d2val2 = q3*(q2*x%d2val2 - pow2(x%d1val2))
   end function sqrt_self
   
   function pow2_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q0
      q0 = 2.0_dp*x%val
      unary%val = pow2(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = 2.0_dp*pow2(x%d1val1) + q0*x%d2val1
      unary%d1val1_d1val2 = 2.0_dp*x%d1val1*x%d1val2 + q0*x%d1val1_d1val2
      unary%d2val2 = 2.0_dp*pow2(x%d1val2) + q0*x%d2val2
   end function pow2_self
   
   function pow3_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = 3.0_dp*pow2(x%val)
      q1 = 3.0_dp*x%val
      unary%val = pow3(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q1*(2.0_dp*pow2(x%d1val1) + x%d2val1*x%val)
      unary%d1val1_d1val2 = q1*(2.0_dp*x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(2.0_dp*pow2(x%d1val2) + x%d2val2*x%val)
   end function pow3_self
   
   function pow4_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = 4.0_dp*pow3(x%val)
      q1 = 4.0_dp*pow2(x%val)
      unary%val = pow4(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q1*(3.0_dp*pow2(x%d1val1) + x%d2val1*x%val)
      unary%d1val1_d1val2 = q1*(3.0_dp*x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(3.0_dp*pow2(x%d1val2) + x%d2val2*x%val)
   end function pow4_self
   
   function pow5_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = 5.0_dp*pow4(x%val)
      q1 = 5.0_dp*pow3(x%val)
      unary%val = pow5(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q1*(4.0_dp*pow2(x%d1val1) + x%d2val1*x%val)
      unary%d1val1_d1val2 = q1*(4.0_dp*x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(4.0_dp*pow2(x%d1val2) + x%d2val2*x%val)
   end function pow5_self
   
   function pow6_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = 6.0_dp*pow5(x%val)
      q1 = 6.0_dp*pow4(x%val)
      unary%val = pow6(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q1*(5.0_dp*pow2(x%d1val1) + x%d2val1*x%val)
      unary%d1val1_d1val2 = q1*(5.0_dp*x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(5.0_dp*pow2(x%d1val2) + x%d2val2*x%val)
   end function pow6_self
   
   function pow7_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = 7.0_dp*pow6(x%val)
      q1 = 7.0_dp*pow5(x%val)
      unary%val = pow7(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q1*(6.0_dp*pow2(x%d1val1) + x%d2val1*x%val)
      unary%d1val1_d1val2 = q1*(6.0_dp*x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(6.0_dp*pow2(x%d1val2) + x%d2val2*x%val)
   end function pow7_self
   
   function pow8_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = 8.0_dp*pow7(x%val)
      q1 = 8.0_dp*pow6(x%val)
      unary%val = pow8(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q1*(7.0_dp*pow2(x%d1val1) + x%d2val1*x%val)
      unary%d1val1_d1val2 = q1*(7.0_dp*x%d1val1*x%d1val2 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(7.0_dp*pow2(x%d1val2) + x%d2val2*x%val)
   end function pow8_self
   
   function abs_self(x) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q0
      q0 = sgn(x%val)
      unary%val = Abs(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
   end function abs_self
   
   function add_self(x, y) result(binary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2), intent(in) :: y
      type(auto_diff_real_2var_order2) :: binary
      binary%val = x%val + y%val
      binary%d1val1 = x%d1val1 + y%d1val1
      binary%d1val2 = x%d1val2 + y%d1val2
      binary%d2val1 = x%d2val1 + y%d2val1
      binary%d1val1_d1val2 = x%d1val1_d1val2 + y%d1val1_d1val2
      binary%d2val2 = x%d2val2 + y%d2val2
   end function add_self
   
   function add_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      unary%val = x%val + y
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
   end function add_self_real
   
   function add_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      unary%val = x%val + z
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
   end function add_real_self
   
   function add_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val + y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
   end function add_self_int
   
   function add_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = x%val + y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
   end function add_int_self
   
   function sub_self(x, y) result(binary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2), intent(in) :: y
      type(auto_diff_real_2var_order2) :: binary
      binary%val = x%val - y%val
      binary%d1val1 = x%d1val1 - y%d1val1
      binary%d1val2 = x%d1val2 - y%d1val2
      binary%d2val1 = x%d2val1 - y%d2val1
      binary%d1val1_d1val2 = x%d1val1_d1val2 - y%d1val1_d1val2
      binary%d2val2 = x%d2val2 - y%d2val2
   end function sub_self
   
   function sub_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      unary%val = x%val - y
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
   end function sub_self_real
   
   function sub_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      unary%val = -x%val + z
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d2val2 = -x%d2val2
   end function sub_real_self
   
   function sub_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val - y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d2val2 = x%d2val2
   end function sub_self_int
   
   function sub_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = -x%val + y_dp
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d2val2 = -x%d2val2
   end function sub_int_self
   
   function mul_self(x, y) result(binary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2), intent(in) :: y
      type(auto_diff_real_2var_order2) :: binary
      binary%val = x%val*y%val
      binary%d1val1 = x%d1val1*y%val + x%val*y%d1val1
      binary%d1val2 = x%d1val2*y%val + x%val*y%d1val2
      binary%d2val1 = 2.0_dp*x%d1val1*y%d1val1 + x%d2val1*y%val + x%val*y%d2val1
      binary%d1val1_d1val2 = x%d1val1*y%d1val2 + x%d1val1_d1val2*y%val + x%d1val2*y%d1val1 + x%val*y%d1val1_d1val2
      binary%d2val2 = 2.0_dp*x%d1val2*y%d1val2 + x%d2val2*y%val + x%val*y%d2val2
   end function mul_self
   
   function mul_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      unary%val = x%val*y
      unary%d1val1 = x%d1val1*y
      unary%d1val2 = x%d1val2*y
      unary%d2val1 = x%d2val1*y
      unary%d1val1_d1val2 = x%d1val1_d1val2*y
      unary%d2val2 = x%d2val2*y
   end function mul_self_real
   
   function mul_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      unary%val = x%val*z
      unary%d1val1 = x%d1val1*z
      unary%d1val2 = x%d1val2*z
      unary%d2val1 = x%d2val1*z
      unary%d1val1_d1val2 = x%d1val1_d1val2*z
      unary%d2val2 = x%d2val2*z
   end function mul_real_self
   
   function mul_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val*y_dp
      unary%d1val1 = x%d1val1*y_dp
      unary%d1val2 = x%d1val2*y_dp
      unary%d2val1 = x%d2val1*y_dp
      unary%d1val1_d1val2 = x%d1val1_d1val2*y_dp
      unary%d2val2 = x%d2val2*y_dp
   end function mul_self_int
   
   function mul_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = x%val*y_dp
      unary%d1val1 = x%d1val1*y_dp
      unary%d1val2 = x%d1val2*y_dp
      unary%d2val1 = x%d2val1*y_dp
      unary%d1val1_d1val2 = x%d1val1_d1val2*y_dp
      unary%d2val2 = x%d2val2*y_dp
   end function mul_int_self
   
   function div_self(x, y) result(binary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2), intent(in) :: y
      type(auto_diff_real_2var_order2) :: binary
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
      q5 = q0*x%d1val2
      q6 = 2.0_dp*y%d1val1
      q7 = 2.0_dp*q0
      binary%val = q1
      binary%d1val1 = q2 - q4*y%d1val1
      binary%d1val2 = -q4*y%d1val2 + q5
      binary%d2val1 = q0*(-q1*(-q7*pow2(y%d1val1) + y%d2val1) - q2*q6 + x%d2val1)
      binary%d1val1_d1val2 = q0*x%d1val1_d1val2 - q3*x%d1val1*y%d1val2 - q3*x%d1val2*y%d1val1 - q4*y%d1val1_d1val2 + q6*x%val*y%d1val2*powm1(pow3(y%val))
      binary%d2val2 = q0*(-2.0_dp*q5*y%d1val2 - q1*(-q7*pow2(y%d1val2) + y%d2val2) + x%d2val2)
   end function div_self
   
   function div_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q0
      q0 = powm1(y)
      unary%val = q0*x%val
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
   end function div_self_real
   
   function div_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = z*powm1(pow2(x%val))
      q1 = z*powm1(pow3(x%val))
      unary%val = z*powm1(x%val)
      unary%d1val1 = -q0*x%d1val1
      unary%d1val2 = -q0*x%d1val2
      unary%d2val1 = q1*(2.0_dp*pow2(x%d1val1) - x%d2val1*x%val)
      unary%d1val1_d1val2 = q1*(2.0_dp*x%d1val1*x%d1val2 - x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(2.0_dp*pow2(x%d1val2) - x%d2val2*x%val)
   end function div_real_self
   
   function div_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
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
   end function div_self_int
   
   function div_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q1
      real(dp) :: q0
      y_dp = z
      q0 = y_dp*powm1(pow2(x%val))
      q1 = y_dp*powm1(pow3(x%val))
      unary%val = y_dp*powm1(x%val)
      unary%d1val1 = -q0*x%d1val1
      unary%d1val2 = -q0*x%d1val2
      unary%d2val1 = q1*(2.0_dp*pow2(x%d1val1) - x%d2val1*x%val)
      unary%d1val1_d1val2 = q1*(2.0_dp*x%d1val1*x%d1val2 - x%d1val1_d1val2*x%val)
      unary%d2val2 = q1*(2.0_dp*pow2(x%d1val2) - x%d2val2*x%val)
   end function div_int_self
   
   function pow_self(x, y) result(binary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2), intent(in) :: y
      type(auto_diff_real_2var_order2) :: binary
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
      q2 = powm1(x%val)
      q3 = q2*y%val
      q4 = q1*y%d1val1 + q3*x%d1val1
      q5 = q0*q4
      q6 = q1*y%d1val2 + q3*x%d1val2
      q7 = q2*x%d1val1
      q8 = y%val*powm1(pow2(x%val))
      q9 = q2*x%d1val2
      binary%val = q0
      binary%d1val1 = q5
      binary%d1val2 = q0*q6
      binary%d2val1 = q0*(2.0_dp*q7*y%d1val1 + q1*y%d2val1 + q3*x%d2val1 - q8*pow2(x%d1val1) + pow2(q4))
      binary%d1val1_d1val2 = q0*(q1*y%d1val1_d1val2 + q3*x%d1val1_d1val2 + q7*y%d1val2 - q8*x%d1val1*x%d1val2 + q9*y%d1val1) + q5*q6
      binary%d2val2 = q0*(2.0_dp*q9*y%d1val2 + q1*y%d2val2 + q3*x%d2val2 - q8*pow2(x%d1val2) + pow2(q6))
   end function pow_self
   
   function pow_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = y*pow(x%val, y - 1)
      q1 = pow2(x%d1val1)
      q2 = y*pow(x%val, -2.0_dp + y)
      q3 = x%d1val1*x%d1val2
      q4 = pow2(x%d1val2)
      unary%val = pow(x%val, y)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*(q1*y - q1 + x%d2val1*x%val)
      unary%d1val1_d1val2 = q2*(q3*y - q3 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q2*(q4*y - q4 + x%d2val2*x%val)
   end function pow_self_real
   
   function pow_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = pow(z, x%val)
      q1 = log(z)
      q2 = q0*q1
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*(q1*pow2(x%d1val1) + x%d2val1)
      unary%d1val1_d1val2 = q2*(q1*x%d1val1*x%d1val2 + x%d1val1_d1val2)
      unary%d2val2 = q2*(q1*pow2(x%d1val2) + x%d2val2)
   end function pow_real_self
   
   function pow_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q4
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      y_dp = y
      q0 = y_dp*pow(x%val, y_dp - 1)
      q1 = pow2(x%d1val1)
      q2 = y_dp*pow(x%val, -2.0_dp + y_dp)
      q3 = x%d1val1*x%d1val2
      q4 = pow2(x%d1val2)
      unary%val = pow(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q2*(q1*y_dp - q1 + x%d2val1*x%val)
      unary%d1val1_d1val2 = q2*(q3*y_dp - q3 + x%d1val1_d1val2*x%val)
      unary%d2val2 = q2*(q4*y_dp - q4 + x%d2val2*x%val)
   end function pow_self_int
   
   function pow_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      y_dp = z
      q0 = pow(y_dp, x%val)
      q1 = log(y_dp)
      q2 = q0*q1
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d2val1 = q2*(q1*pow2(x%d1val1) + x%d2val1)
      unary%d1val1_d1val2 = q2*(q1*x%d1val1*x%d1val2 + x%d1val1_d1val2)
      unary%d2val2 = q2*(q1*pow2(x%d1val2) + x%d2val2)
   end function pow_int_self
   
   function max_self(x, y) result(binary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2), intent(in) :: y
      type(auto_diff_real_2var_order2) :: binary
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
   end function max_self
   
   function max_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q0
      q0 = Heaviside(x%val - y)
      unary%val = Max(x%val, y)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
   end function max_self_real
   
   function max_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q0
      q0 = Heaviside(x%val - z)
      unary%val = Max(x%val, z)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
   end function max_real_self
   
   function max_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
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
   end function max_self_int
   
   function max_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
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
   end function max_int_self
   
   function min_self(x, y) result(binary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2), intent(in) :: y
      type(auto_diff_real_2var_order2) :: binary
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
   end function min_self
   
   function min_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q0
      q0 = Heaviside(-x%val + y)
      unary%val = Min(x%val, y)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
   end function min_self_real
   
   function min_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q0
      q0 = Heaviside(-x%val + z)
      unary%val = Min(x%val, z)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d2val2 = q0*x%d2val2
   end function min_real_self
   
   function min_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
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
   end function min_self_int
   
   function min_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
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
   end function min_int_self
   
   function dim_self(x, y) result(binary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2), intent(in) :: y
      type(auto_diff_real_2var_order2) :: binary
      real(dp) :: q3
      real(dp) :: q2
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val - y%val
      q1 = sgn(q0)
      q2 = 0.5_dp*q1
      q3 = -0.5_dp*y%d1val1_d1val2 + 0.5_dp*x%d1val1_d1val2
      binary%val = -0.5_dp*y%val + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      binary%d1val1 = -0.5_dp*y%d1val1 + 0.5_dp*x%d1val1 + q2*(x%d1val1 - y%d1val1)
      binary%d1val2 = -0.5_dp*y%d1val2 + 0.5_dp*x%d1val2 + q2*(x%d1val2 - y%d1val2)
      binary%d2val1 = -0.5_dp*y%d2val1 + 0.5_dp*x%d2val1 + q2*(x%d2val1 - y%d2val1)
      binary%d1val1_d1val2 = q1*q3 + q3
      binary%d2val2 = -0.5_dp*y%d2val2 + 0.5_dp*x%d2val2 + q2*(x%d2val2 - y%d2val2)
   end function dim_self
   
   function dim_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val - y
      q1 = 0.5_dp*sgn(q0) + 0.5_dp
      unary%val = -0.5_dp*y + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2
      unary%d2val2 = q1*x%d2val2
   end function dim_self_real
   
   function dim_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: q1
      real(dp) :: q0
      q0 = x%val - z
      q1 = -0.5_dp + 0.5_dp*sgn(q0)
      unary%val = -0.5_dp*x%val + 0.5_dp*z + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2
      unary%d2val2 = q1*x%d2val2
   end function dim_real_self
   
   function dim_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q1
      real(dp) :: q0
      y_dp = y
      q0 = x%val - y_dp
      q1 = 0.5_dp*sgn(q0) + 0.5_dp
      unary%val = -0.5_dp*y_dp + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2
      unary%d2val2 = q1*x%d2val2
   end function dim_self_int
   
   function dim_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order2), intent(in) :: x
      type(auto_diff_real_2var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q1
      real(dp) :: q0
      y_dp = z
      q0 = x%val - y_dp
      q1 = -0.5_dp + 0.5_dp*sgn(q0)
      unary%val = -0.5_dp*x%val + 0.5_dp*y_dp + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d2val1 = q1*x%d2val1
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2
      unary%d2val2 = q1*x%d2val2
   end function dim_int_self
   
   function differentiate_auto_diff_real_2var_order2_1(this) result(derivative)
      type(auto_diff_real_2var_order2), intent(in) :: this
      type(auto_diff_real_2var_order2) :: derivative
      derivative%val = this%d1val1
      derivative%d1val1 = this%d2val1
      derivative%d1val2 = this%d1val1_d1val2
      derivative%d2val1 = 0_dp
      derivative%d1val1_d1val2 = 0_dp
      derivative%d2val2 = 0_dp
   end function differentiate_auto_diff_real_2var_order2_1
   
   function differentiate_auto_diff_real_2var_order2_2(this) result(derivative)
      type(auto_diff_real_2var_order2), intent(in) :: this
      type(auto_diff_real_2var_order2) :: derivative
      derivative%val = this%d1val2
      derivative%d1val1 = this%d1val1_d1val2
      derivative%d1val2 = this%d2val2
      derivative%d2val1 = 0_dp
      derivative%d1val1_d1val2 = 0_dp
      derivative%d2val2 = 0_dp
   end function differentiate_auto_diff_real_2var_order2_2
   
end module auto_diff_real_2var_order2_module