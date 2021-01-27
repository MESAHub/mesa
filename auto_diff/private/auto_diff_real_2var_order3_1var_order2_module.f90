module auto_diff_real_2var_order3_1var_order2_module
      use const_def
      use utils_lib
      use support_functions
      use math_lib
   
      implicit none
      private
   public :: auto_diff_real_2var_order3_1var_order2, &
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
      differentiate_2, &
      differentiate_3
   type :: auto_diff_real_2var_order3_1var_order2
      real(dp) :: val
      real(dp) :: d1val1
      real(dp) :: d1val2
      real(dp) :: d1val3
      real(dp) :: d2val1
      real(dp) :: d1val1_d1val2
      real(dp) :: d1val1_d1val3
      real(dp) :: d2val2
      real(dp) :: d1val2_d1val3
      real(dp) :: d2val3
      real(dp) :: d3val1
      real(dp) :: d2val1_d1val2
      real(dp) :: d2val1_d1val3
      real(dp) :: d1val1_d2val2
      real(dp) :: d1val1_d1val2_d1val3
      real(dp) :: d1val1_d2val3
      real(dp) :: d3val2
      real(dp) :: d2val2_d1val3
      real(dp) :: d1val2_d2val3
      real(dp) :: d3val1_d1val3
      real(dp) :: d2val1_d1val2_d1val3
      real(dp) :: d2val1_d2val3
      real(dp) :: d1val1_d2val2_d1val3
      real(dp) :: d1val1_d1val2_d2val3
      real(dp) :: d3val2_d1val3
      real(dp) :: d2val2_d2val3
      real(dp) :: d3val1_d2val3
      real(dp) :: d2val1_d1val2_d2val3
      real(dp) :: d1val1_d2val2_d2val3
      real(dp) :: d3val2_d2val3
   end type auto_diff_real_2var_order3_1var_order2
   
   interface assignment(=)
      module procedure assign_from_self
      module procedure assign_from_real_dp
      module procedure assign_from_int
   end interface assignment(=)
   
   interface operator(.eq.)
      module procedure equal_self
      module procedure equal_auto_diff_real_2var_order3_1var_order2_real_dp
      module procedure equal_real_dp_auto_diff_real_2var_order3_1var_order2
      module procedure equal_auto_diff_real_2var_order3_1var_order2_int
      module procedure equal_int_auto_diff_real_2var_order3_1var_order2
   end interface operator(.eq.)
   
   interface operator(.ne.)
      module procedure neq_self
      module procedure neq_auto_diff_real_2var_order3_1var_order2_real_dp
      module procedure neq_real_dp_auto_diff_real_2var_order3_1var_order2
      module procedure neq_auto_diff_real_2var_order3_1var_order2_int
      module procedure neq_int_auto_diff_real_2var_order3_1var_order2
   end interface operator(.ne.)
   
   interface operator(.gt.)
      module procedure greater_self
      module procedure greater_auto_diff_real_2var_order3_1var_order2_real_dp
      module procedure greater_real_dp_auto_diff_real_2var_order3_1var_order2
      module procedure greater_auto_diff_real_2var_order3_1var_order2_int
      module procedure greater_int_auto_diff_real_2var_order3_1var_order2
   end interface operator(.gt.)
   
   interface operator(.lt.)
      module procedure less_self
      module procedure less_auto_diff_real_2var_order3_1var_order2_real_dp
      module procedure less_real_dp_auto_diff_real_2var_order3_1var_order2
      module procedure less_auto_diff_real_2var_order3_1var_order2_int
      module procedure less_int_auto_diff_real_2var_order3_1var_order2
   end interface operator(.lt.)
   
   interface operator(.le.)
      module procedure leq_self
      module procedure leq_auto_diff_real_2var_order3_1var_order2_real_dp
      module procedure leq_real_dp_auto_diff_real_2var_order3_1var_order2
      module procedure leq_auto_diff_real_2var_order3_1var_order2_int
      module procedure leq_int_auto_diff_real_2var_order3_1var_order2
   end interface operator(.le.)
   
   interface operator(.ge.)
      module procedure geq_self
      module procedure geq_auto_diff_real_2var_order3_1var_order2_real_dp
      module procedure geq_real_dp_auto_diff_real_2var_order3_1var_order2
      module procedure geq_auto_diff_real_2var_order3_1var_order2_int
      module procedure geq_int_auto_diff_real_2var_order3_1var_order2
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
      module procedure differentiate_auto_diff_real_2var_order3_1var_order2_1
   end interface differentiate_1
   
   interface differentiate_2
      module procedure differentiate_auto_diff_real_2var_order3_1var_order2_2
   end interface differentiate_2
   
   interface differentiate_3
      module procedure differentiate_auto_diff_real_2var_order3_1var_order2_3
   end interface differentiate_3
   
   contains

   subroutine assign_from_self(this, other)
      type(auto_diff_real_2var_order3_1var_order2), intent(out) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      this%val = other%val
      this%d1val1 = other%d1val1
      this%d1val2 = other%d1val2
      this%d1val3 = other%d1val3
      this%d2val1 = other%d2val1
      this%d1val1_d1val2 = other%d1val1_d1val2
      this%d1val1_d1val3 = other%d1val1_d1val3
      this%d2val2 = other%d2val2
      this%d1val2_d1val3 = other%d1val2_d1val3
      this%d2val3 = other%d2val3
      this%d3val1 = other%d3val1
      this%d2val1_d1val2 = other%d2val1_d1val2
      this%d2val1_d1val3 = other%d2val1_d1val3
      this%d1val1_d2val2 = other%d1val1_d2val2
      this%d1val1_d1val2_d1val3 = other%d1val1_d1val2_d1val3
      this%d1val1_d2val3 = other%d1val1_d2val3
      this%d3val2 = other%d3val2
      this%d2val2_d1val3 = other%d2val2_d1val3
      this%d1val2_d2val3 = other%d1val2_d2val3
      this%d3val1_d1val3 = other%d3val1_d1val3
      this%d2val1_d1val2_d1val3 = other%d2val1_d1val2_d1val3
      this%d2val1_d2val3 = other%d2val1_d2val3
      this%d1val1_d2val2_d1val3 = other%d1val1_d2val2_d1val3
      this%d1val1_d1val2_d2val3 = other%d1val1_d1val2_d2val3
      this%d3val2_d1val3 = other%d3val2_d1val3
      this%d2val2_d2val3 = other%d2val2_d2val3
      this%d3val1_d2val3 = other%d3val1_d2val3
      this%d2val1_d1val2_d2val3 = other%d2val1_d1val2_d2val3
      this%d1val1_d2val2_d2val3 = other%d1val1_d2val2_d2val3
      this%d3val2_d2val3 = other%d3val2_d2val3
   end subroutine assign_from_self
   
   subroutine assign_from_real_dp(this, other)
      type(auto_diff_real_2var_order3_1var_order2), intent(out) :: this
      real(dp), intent(in) :: other
      this%val = other
      this%d1val1 = 0_dp
      this%d1val2 = 0_dp
      this%d1val3 = 0_dp
      this%d2val1 = 0_dp
      this%d1val1_d1val2 = 0_dp
      this%d1val1_d1val3 = 0_dp
      this%d2val2 = 0_dp
      this%d1val2_d1val3 = 0_dp
      this%d2val3 = 0_dp
      this%d3val1 = 0_dp
      this%d2val1_d1val2 = 0_dp
      this%d2val1_d1val3 = 0_dp
      this%d1val1_d2val2 = 0_dp
      this%d1val1_d1val2_d1val3 = 0_dp
      this%d1val1_d2val3 = 0_dp
      this%d3val2 = 0_dp
      this%d2val2_d1val3 = 0_dp
      this%d1val2_d2val3 = 0_dp
      this%d3val1_d1val3 = 0_dp
      this%d2val1_d1val2_d1val3 = 0_dp
      this%d2val1_d2val3 = 0_dp
      this%d1val1_d2val2_d1val3 = 0_dp
      this%d1val1_d1val2_d2val3 = 0_dp
      this%d3val2_d1val3 = 0_dp
      this%d2val2_d2val3 = 0_dp
      this%d3val1_d2val3 = 0_dp
      this%d2val1_d1val2_d2val3 = 0_dp
      this%d1val1_d2val2_d2val3 = 0_dp
      this%d3val2_d2val3 = 0_dp
   end subroutine assign_from_real_dp
   
   subroutine assign_from_int(this, other)
      type(auto_diff_real_2var_order3_1var_order2), intent(out) :: this
      integer, intent(in) :: other
      this%val = other
      this%d1val1 = 0_dp
      this%d1val2 = 0_dp
      this%d1val3 = 0_dp
      this%d2val1 = 0_dp
      this%d1val1_d1val2 = 0_dp
      this%d1val1_d1val3 = 0_dp
      this%d2val2 = 0_dp
      this%d1val2_d1val3 = 0_dp
      this%d2val3 = 0_dp
      this%d3val1 = 0_dp
      this%d2val1_d1val2 = 0_dp
      this%d2val1_d1val3 = 0_dp
      this%d1val1_d2val2 = 0_dp
      this%d1val1_d1val2_d1val3 = 0_dp
      this%d1val1_d2val3 = 0_dp
      this%d3val2 = 0_dp
      this%d2val2_d1val3 = 0_dp
      this%d1val2_d2val3 = 0_dp
      this%d3val1_d1val3 = 0_dp
      this%d2val1_d1val2_d1val3 = 0_dp
      this%d2val1_d2val3 = 0_dp
      this%d1val1_d2val2_d1val3 = 0_dp
      this%d1val1_d1val2_d2val3 = 0_dp
      this%d3val2_d1val3 = 0_dp
      this%d2val2_d2val3 = 0_dp
      this%d3val1_d2val3 = 0_dp
      this%d2val1_d1val2_d2val3 = 0_dp
      this%d1val1_d2val2_d2val3 = 0_dp
      this%d3val2_d2val3 = 0_dp
   end subroutine assign_from_int
   
   function equal_self(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this%val .eq. other%val)
   end function equal_self
   
   function equal_auto_diff_real_2var_order3_1var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .eq. other)
   end function equal_auto_diff_real_2var_order3_1var_order2_real_dp
   
   function equal_real_dp_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .eq. other%val)
   end function equal_real_dp_auto_diff_real_2var_order3_1var_order2
   
   function equal_auto_diff_real_2var_order3_1var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .eq. other)
   end function equal_auto_diff_real_2var_order3_1var_order2_int
   
   function equal_int_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .eq. other%val)
   end function equal_int_auto_diff_real_2var_order3_1var_order2
   
   function neq_self(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this%val .ne. other%val)
   end function neq_self
   
   function neq_auto_diff_real_2var_order3_1var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .ne. other)
   end function neq_auto_diff_real_2var_order3_1var_order2_real_dp
   
   function neq_real_dp_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .ne. other%val)
   end function neq_real_dp_auto_diff_real_2var_order3_1var_order2
   
   function neq_auto_diff_real_2var_order3_1var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .ne. other)
   end function neq_auto_diff_real_2var_order3_1var_order2_int
   
   function neq_int_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .ne. other%val)
   end function neq_int_auto_diff_real_2var_order3_1var_order2
   
   function greater_self(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this%val .gt. other%val)
   end function greater_self
   
   function greater_auto_diff_real_2var_order3_1var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .gt. other)
   end function greater_auto_diff_real_2var_order3_1var_order2_real_dp
   
   function greater_real_dp_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .gt. other%val)
   end function greater_real_dp_auto_diff_real_2var_order3_1var_order2
   
   function greater_auto_diff_real_2var_order3_1var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .gt. other)
   end function greater_auto_diff_real_2var_order3_1var_order2_int
   
   function greater_int_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .gt. other%val)
   end function greater_int_auto_diff_real_2var_order3_1var_order2
   
   function less_self(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this%val .lt. other%val)
   end function less_self
   
   function less_auto_diff_real_2var_order3_1var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .lt. other)
   end function less_auto_diff_real_2var_order3_1var_order2_real_dp
   
   function less_real_dp_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .lt. other%val)
   end function less_real_dp_auto_diff_real_2var_order3_1var_order2
   
   function less_auto_diff_real_2var_order3_1var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .lt. other)
   end function less_auto_diff_real_2var_order3_1var_order2_int
   
   function less_int_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .lt. other%val)
   end function less_int_auto_diff_real_2var_order3_1var_order2
   
   function leq_self(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this%val .le. other%val)
   end function leq_self
   
   function leq_auto_diff_real_2var_order3_1var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .le. other)
   end function leq_auto_diff_real_2var_order3_1var_order2_real_dp
   
   function leq_real_dp_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .le. other%val)
   end function leq_real_dp_auto_diff_real_2var_order3_1var_order2
   
   function leq_auto_diff_real_2var_order3_1var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .le. other)
   end function leq_auto_diff_real_2var_order3_1var_order2_int
   
   function leq_int_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .le. other%val)
   end function leq_int_auto_diff_real_2var_order3_1var_order2
   
   function geq_self(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this%val .ge. other%val)
   end function geq_self
   
   function geq_auto_diff_real_2var_order3_1var_order2_real_dp(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      real(dp), intent(in) :: other
      logical :: z
      z = (this%val .ge. other)
   end function geq_auto_diff_real_2var_order3_1var_order2_real_dp
   
   function geq_real_dp_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      real(dp), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .ge. other%val)
   end function geq_real_dp_auto_diff_real_2var_order3_1var_order2
   
   function geq_auto_diff_real_2var_order3_1var_order2_int(this, other) result(z)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      integer, intent(in) :: other
      logical :: z
      z = (this%val .ge. other)
   end function geq_auto_diff_real_2var_order3_1var_order2_int
   
   function geq_int_auto_diff_real_2var_order3_1var_order2(this, other) result(z)
      integer, intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: other
      logical :: z
      z = (this .ge. other%val)
   end function geq_int_auto_diff_real_2var_order3_1var_order2
   
   function make_unary_operator(x, z_val, z_d1x, z_d2x, z_d3x, z_d4x, z_d5x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      real(dp), intent(in) :: z_val
      real(dp), intent(in) :: z_d1x
      real(dp), intent(in) :: z_d2x
      real(dp), intent(in) :: z_d3x
      real(dp), intent(in) :: z_d4x
      real(dp), intent(in) :: z_d5x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q1 = x%d1val1*z_d2x
      q2 = pow2(x%d1val2)
      q3 = x%d1val2*z_d2x
      q4 = pow2(x%d1val3)
      q5 = pow3(x%d1val1)
      q6 = 3.0_dp*q1
      q7 = 2.0_dp*q1
      q8 = q0*z_d3x
      q9 = x%d1val3*z_d2x
      q10 = 2.0_dp*q3
      q11 = x%d1val1*z_d3x
      q12 = x%d1val2*x%d1val3
      q13 = 2.0_dp*q9
      q14 = pow3(x%d1val2)
      q15 = 3.0_dp*q3
      q16 = q2*z_d3x
      q17 = q4*z_d3x
      q18 = x%d2val1*z_d2x
      q19 = 3.0_dp*x%d1val1_d1val3
      q20 = x%d1val3*z_d4x
      q21 = q11*x%d1val3
      q22 = 3.0_dp*x%d2val1
      q23 = 2.0_dp*z_d2x
      q24 = q23*x%d1val1_d1val2
      q25 = x%d2val1*z_d3x
      q26 = 2.0_dp*q11
      q27 = q26*x%d1val3
      q28 = q26*x%d1val2
      q29 = q0*x%d1val2
      q30 = pow2(x%d1val1_d1val3)
      q31 = 4.0_dp*x%d1val1_d1val3
      q32 = q4*z_d4x
      q33 = x%d2val2*z_d2x
      q34 = q11*x%d2val2
      q35 = q12*z_d3x
      q36 = 2.0_dp*x%d1val1_d1val2
      q37 = q2*x%d1val1
      q38 = x%d2val3*z_d2x
      q39 = q23*x%d1val2_d1val3
      q40 = x%d1val2*x%d2val3
      q41 = 2.0_dp*q35
      q42 = q32*x%d1val2
      q43 = 3.0_dp*x%d1val2_d1val3
      q44 = 3.0_dp*x%d2val2
      q45 = pow2(x%d1val2_d1val3)
      q46 = 4.0_dp*x%d1val2_d1val3
      q47 = x%d1val1_d1val3*z_d2x
      q48 = 3.0_dp*x%d1val1_d2val3
      q49 = x%d2val3*z_d4x
      q50 = 6.0_dp*q11
      q51 = 6.0_dp*x%d1val1_d1val3
      q52 = q25*x%d1val3
      q53 = q4*z_d5x
      q54 = q32*x%d1val1
      q55 = q0*q20
      q56 = 4.0_dp*x%d1val1_d1val2_d1val3
      q57 = q31*x%d1val2_d1val3
      q58 = x%d1val3*z_d3x
      q59 = q58*x%d1val1_d1val2
      q60 = 2.0_dp*x%d1val2_d1val3
      q61 = x%d1val2*z_d3x
      q62 = q20*x%d1val1*x%d1val2
      q63 = x%d1val2_d1val3*z_d2x
      q64 = q40*z_d3x
      q65 = 2.0_dp*x%d1val1_d1val3
      q66 = q58*x%d2val2
      q67 = q2*q20
      q68 = 6.0_dp*x%d2val2_d1val3
      q69 = 3.0_dp*x%d1val2_d2val3
      q70 = 6.0_dp*x%d1val2_d1val3
      unary%val = z_val
      unary%d1val1 = x%d1val1*z_d1x
      unary%d1val2 = x%d1val2*z_d1x
      unary%d1val3 = x%d1val3*z_d1x
      unary%d2val1 = q0*z_d2x + x%d2val1*z_d1x
      unary%d1val1_d1val2 = q1*x%d1val2 + x%d1val1_d1val2*z_d1x
      unary%d1val1_d1val3 = q1*x%d1val3 + x%d1val1_d1val3*z_d1x
      unary%d2val2 = q2*z_d2x + x%d2val2*z_d1x
      unary%d1val2_d1val3 = q3*x%d1val3 + x%d1val2_d1val3*z_d1x
      unary%d2val3 = q4*z_d2x + x%d2val3*z_d1x
      unary%d3val1 = q5*z_d3x + q6*x%d2val1 + x%d3val1*z_d1x
      unary%d2val1_d1val2 = q3*x%d2val1 + q7*x%d1val1_d1val2 + q8*x%d1val2 + x%d2val1_d1val2*z_d1x
      unary%d2val1_d1val3 = q7*x%d1val1_d1val3 + q8*x%d1val3 + q9*x%d2val1 + x%d2val1_d1val3*z_d1x
      unary%d1val1_d2val2 = q1*x%d2val2 + q10*x%d1val1_d1val2 + q11*q2 + x%d1val1_d2val2*z_d1x
      unary%d1val1_d1val2_d1val3 = q1*x%d1val2_d1val3 + q11*q12 + q3*x%d1val1_d1val3 + q9*x%d1val1_d1val2 + x%d1val1_d1val2_d1val3*z_d1x
      unary%d1val1_d2val3 = q1*x%d2val3 + q11*q4 + q13*x%d1val1_d1val3 + x%d1val1_d2val3*z_d1x
      unary%d3val2 = q14*z_d3x + q15*x%d2val2 + x%d3val2*z_d1x
      unary%d2val2_d1val3 = q10*x%d1val2_d1val3 + q16*x%d1val3 + q9*x%d2val2 + x%d2val2_d1val3*z_d1x
      unary%d1val2_d2val3 = q13*x%d1val2_d1val3 + q17*x%d1val2 + q3*x%d2val3 + x%d1val2_d2val3*z_d1x
      unary%d3val1_d1val3 = q18*q19 + q19*q8 + q20*q5 + q21*q22 + q6*x%d2val1_d1val3 + q9*x%d3val1 + x%d3val1_d1val3*z_d1x
      unary%d2val1_d1val2_d1val3 = q12*q25 + q18*x%d1val2_d1val3 + q20*q29 + q24*x%d1val1_d1val3 + q27*x%d1val1_d1val2 + q28*x%d1val1_d1val3 + q3*x%d2val1_d1val3 + q7*x%d1val1_d1val2_d1val3 + q8*x%d1val2_d1val3 + q9*x%d2val1_d1val2 + x%d2val1_d1val2_d1val3*z_d1x
      unary%d2val1_d2val3 = q0*q32 + q13*x%d2val1_d1val3 + q17*x%d2val1 + q18*x%d2val3 + q21*q31 + q23*q30 + q7*x%d1val1_d2val3 + q8*x%d2val3 + x%d2val1_d2val3*z_d1x
      unary%d1val1_d2val2_d1val3 = q1*x%d2val2_d1val3 + q10*x%d1val1_d1val2_d1val3 + q16*x%d1val1_d1val3 + q20*q37 + q24*x%d1val2_d1val3 + q28*x%d1val2_d1val3 + q33*x%d1val1_d1val3 + q34*x%d1val3 + q35*q36 + q9*x%d1val1_d2val2 + x%d1val1_d2val2_d1val3*z_d1x
      unary%d1val1_d1val2_d2val3 = q1*x%d1val2_d2val3 + q11*q40 + q13*x%d1val1_d1val2_d1val3 + q17*x%d1val1_d1val2 + q27*x%d1val2_d1val3 + q3*x%d1val1_d2val3 + q38*x%d1val1_d1val2 + q39*x%d1val1_d1val3 + q41*x%d1val1_d1val3 + q42*x%d1val1 + x%d1val1_d1val2_d2val3*z_d1x
      unary%d3val2_d1val3 = q14*q20 + q15*x%d2val2_d1val3 + q16*q43 + q33*q43 + q35*q44 + q9*x%d3val2 + x%d3val2_d1val3*z_d1x
      unary%d2val2_d2val3 = q10*x%d1val2_d2val3 + q13*x%d2val2_d1val3 + q16*x%d2val3 + q17*x%d2val2 + q2*q32 + q23*q45 + q33*x%d2val3 + q35*q46 + x%d2val2_d2val3*z_d1x
      unary%d3val1_d2val3 = 6.0_dp*q47*x%d2val1_d1val3 + q11*q22*x%d2val3 + q13*x%d3val1_d1val3 + q17*x%d3val1 + q18*q48 + q22*q54 + q30*q50 + q38*x%d3val1 + q48*q8 + q49*q5 + q5*q53 + q50*x%d1val3*x%d2val1_d1val3 + q51*q52 + q51*q55 + q6*x%d2val1_d2val3 + x%d3val1_d2val3*z_d1x
      unary%d2val1_d1val2_d2val3 = 2.0_dp*q30*q61 + q11*q57 + q13*x%d2val1_d1val2_d1val3 + q17*x%d2val1_d1val2 + q18*x%d1val2_d2val3 + q21*q56 + q24*x%d1val1_d2val3 + q25*q40 + q26*x%d1val1_d1val2*x%d2val3 + q28*x%d1val1_d2val3 + q29*q49 + q29*q53 + q3*x%d2val1_d2val3 + q31*q59 + q31*q62 + q36*q54 + q38*x%d2val1_d1val2 + q39*x%d2val1_d1val3 + q41*x%d2val1_d1val3 + q42*x%d2val1 + q47*q56 + q52*q60 + q55*q60 + q7*x%d1val1_d1val2_d2val3 + q8*x%d1val2_d2val3 + x%d2val1_d1val2_d2val3*z_d1x
      unary%d1val1_d2val2_d2val3 = q1*x%d2val2_d2val3 + q10*x%d1val1_d1val2_d2val3 + q13*x%d1val1_d2val2_d1val3 + q16*x%d1val1_d2val3 + q17*x%d1val1_d2val2 + q23*x%d1val1_d1val3*x%d2val2_d1val3 + q24*x%d1val2_d2val3 + q26*q45 + q27*x%d2val2_d1val3 + q28*x%d1val2_d2val3 + q33*x%d1val1_d2val3 + q34*x%d2val3 + q35*q56 + q36*q42 + q36*q64 + q37*q49 + q37*q53 + q38*x%d1val1_d2val2 + q46*q59 + q46*q62 + q54*x%d2val2 + q56*q63 + q57*q61 + q65*q66 + q65*q67 + x%d1val1_d2val2_d2val3*z_d1x
      unary%d3val2_d2val3 = 6.0_dp*q45*q61 + q13*x%d3val2_d1val3 + q14*q49 + q14*q53 + q15*x%d2val2_d2val3 + q16*q69 + q17*x%d3val2 + q33*q69 + q35*q68 + q38*x%d3val2 + q42*q44 + q44*q64 + q63*q68 + q66*q70 + q67*q70 + x%d3val2_d2val3*z_d1x
   end function make_unary_operator
   
   function make_binary_operator(x, y, z_val, z_d1x, z_d1y, z_d2x, z_d1x_d1y, z_d2y, z_d3x, z_d2x_d1y, z_d1x_d2y, z_d3y, z_d4x, z_d3x_d1y, z_d2x_d2y, z_d1x_d3y, z_d4y, z_d5x, z_d4x_d1y, z_d3x_d2y, z_d2x_d3y, z_d1x_d4y, z_d5y) result(binary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: y
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
      real(dp), intent(in) :: z_d4x
      real(dp), intent(in) :: z_d3x_d1y
      real(dp), intent(in) :: z_d2x_d2y
      real(dp), intent(in) :: z_d1x_d3y
      real(dp), intent(in) :: z_d4y
      real(dp), intent(in) :: z_d5x
      real(dp), intent(in) :: z_d4x_d1y
      real(dp), intent(in) :: z_d3x_d2y
      real(dp), intent(in) :: z_d2x_d3y
      real(dp), intent(in) :: z_d1x_d4y
      real(dp), intent(in) :: z_d5y
      type(auto_diff_real_2var_order3_1var_order2) :: binary
      real(dp) :: q441
      real(dp) :: q440
      real(dp) :: q439
      real(dp) :: q438
      real(dp) :: q437
      real(dp) :: q436
      real(dp) :: q435
      real(dp) :: q434
      real(dp) :: q433
      real(dp) :: q432
      real(dp) :: q431
      real(dp) :: q430
      real(dp) :: q429
      real(dp) :: q428
      real(dp) :: q427
      real(dp) :: q426
      real(dp) :: q425
      real(dp) :: q424
      real(dp) :: q423
      real(dp) :: q422
      real(dp) :: q421
      real(dp) :: q420
      real(dp) :: q419
      real(dp) :: q418
      real(dp) :: q417
      real(dp) :: q416
      real(dp) :: q415
      real(dp) :: q414
      real(dp) :: q413
      real(dp) :: q412
      real(dp) :: q411
      real(dp) :: q410
      real(dp) :: q409
      real(dp) :: q408
      real(dp) :: q407
      real(dp) :: q406
      real(dp) :: q405
      real(dp) :: q404
      real(dp) :: q403
      real(dp) :: q402
      real(dp) :: q401
      real(dp) :: q400
      real(dp) :: q399
      real(dp) :: q398
      real(dp) :: q397
      real(dp) :: q396
      real(dp) :: q395
      real(dp) :: q394
      real(dp) :: q393
      real(dp) :: q392
      real(dp) :: q391
      real(dp) :: q390
      real(dp) :: q389
      real(dp) :: q388
      real(dp) :: q387
      real(dp) :: q386
      real(dp) :: q385
      real(dp) :: q384
      real(dp) :: q383
      real(dp) :: q382
      real(dp) :: q381
      real(dp) :: q380
      real(dp) :: q379
      real(dp) :: q378
      real(dp) :: q377
      real(dp) :: q376
      real(dp) :: q375
      real(dp) :: q374
      real(dp) :: q373
      real(dp) :: q372
      real(dp) :: q371
      real(dp) :: q370
      real(dp) :: q369
      real(dp) :: q368
      real(dp) :: q367
      real(dp) :: q366
      real(dp) :: q365
      real(dp) :: q364
      real(dp) :: q363
      real(dp) :: q362
      real(dp) :: q361
      real(dp) :: q360
      real(dp) :: q359
      real(dp) :: q358
      real(dp) :: q357
      real(dp) :: q356
      real(dp) :: q355
      real(dp) :: q354
      real(dp) :: q353
      real(dp) :: q352
      real(dp) :: q351
      real(dp) :: q350
      real(dp) :: q349
      real(dp) :: q348
      real(dp) :: q347
      real(dp) :: q346
      real(dp) :: q345
      real(dp) :: q344
      real(dp) :: q343
      real(dp) :: q342
      real(dp) :: q341
      real(dp) :: q340
      real(dp) :: q339
      real(dp) :: q338
      real(dp) :: q337
      real(dp) :: q336
      real(dp) :: q335
      real(dp) :: q334
      real(dp) :: q333
      real(dp) :: q332
      real(dp) :: q331
      real(dp) :: q330
      real(dp) :: q329
      real(dp) :: q328
      real(dp) :: q327
      real(dp) :: q326
      real(dp) :: q325
      real(dp) :: q324
      real(dp) :: q323
      real(dp) :: q322
      real(dp) :: q321
      real(dp) :: q320
      real(dp) :: q319
      real(dp) :: q318
      real(dp) :: q317
      real(dp) :: q316
      real(dp) :: q315
      real(dp) :: q314
      real(dp) :: q313
      real(dp) :: q312
      real(dp) :: q311
      real(dp) :: q310
      real(dp) :: q309
      real(dp) :: q308
      real(dp) :: q307
      real(dp) :: q306
      real(dp) :: q305
      real(dp) :: q304
      real(dp) :: q303
      real(dp) :: q302
      real(dp) :: q301
      real(dp) :: q300
      real(dp) :: q299
      real(dp) :: q298
      real(dp) :: q297
      real(dp) :: q296
      real(dp) :: q295
      real(dp) :: q294
      real(dp) :: q293
      real(dp) :: q292
      real(dp) :: q291
      real(dp) :: q290
      real(dp) :: q289
      real(dp) :: q288
      real(dp) :: q287
      real(dp) :: q286
      real(dp) :: q285
      real(dp) :: q284
      real(dp) :: q283
      real(dp) :: q282
      real(dp) :: q281
      real(dp) :: q280
      real(dp) :: q279
      real(dp) :: q278
      real(dp) :: q277
      real(dp) :: q276
      real(dp) :: q275
      real(dp) :: q274
      real(dp) :: q273
      real(dp) :: q272
      real(dp) :: q271
      real(dp) :: q270
      real(dp) :: q269
      real(dp) :: q268
      real(dp) :: q267
      real(dp) :: q266
      real(dp) :: q265
      real(dp) :: q264
      real(dp) :: q263
      real(dp) :: q262
      real(dp) :: q261
      real(dp) :: q260
      real(dp) :: q259
      real(dp) :: q258
      real(dp) :: q257
      real(dp) :: q256
      real(dp) :: q255
      real(dp) :: q254
      real(dp) :: q253
      real(dp) :: q252
      real(dp) :: q251
      real(dp) :: q250
      real(dp) :: q249
      real(dp) :: q248
      real(dp) :: q247
      real(dp) :: q246
      real(dp) :: q245
      real(dp) :: q244
      real(dp) :: q243
      real(dp) :: q242
      real(dp) :: q241
      real(dp) :: q240
      real(dp) :: q239
      real(dp) :: q238
      real(dp) :: q237
      real(dp) :: q236
      real(dp) :: q235
      real(dp) :: q234
      real(dp) :: q233
      real(dp) :: q232
      real(dp) :: q231
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q9 = x%d1val2*z_d1x_d1y
      q10 = 2.0_dp*q9
      q11 = x%d1val2*z_d2x
      q12 = y%d1val2*z_d1x_d1y
      q13 = y%d1val2*z_d2y
      q14 = pow2(x%d1val3)
      q15 = pow2(y%d1val3)
      q16 = x%d1val3*z_d1x_d1y
      q17 = 2.0_dp*q16
      q18 = pow3(x%d1val1)
      q19 = pow3(y%d1val1)
      q20 = 3.0_dp*x%d2val1
      q21 = 3.0_dp*y%d2val1
      q22 = q1*z_d1x_d2y
      q23 = 3.0_dp*q22
      q24 = q0*z_d2x_d1y
      q25 = 3.0_dp*q24
      q26 = 2.0_dp*x%d1val1_d1val2
      q27 = 2.0_dp*q6
      q28 = q0*z_d3x
      q29 = q1*z_d3y
      q30 = x%d1val1*x%d1val2
      q31 = y%d1val1*z_d2x_d1y
      q32 = 2.0_dp*q31
      q33 = y%d1val1*y%d1val2
      q34 = x%d1val1*z_d1x_d2y
      q35 = 2.0_dp*q34
      q36 = x%d1val3*z_d2x
      q37 = y%d1val3*z_d1x_d1y
      q38 = y%d1val3*z_d2y
      q39 = 2.0_dp*x%d1val1_d1val3
      q40 = x%d1val1*x%d1val3
      q41 = y%d1val1*y%d1val3
      q42 = 2.0_dp*q13
      q43 = x%d1val1*z_d3x
      q44 = y%d1val1*z_d3y
      q45 = x%d1val1*z_d2x_d1y
      q46 = 2.0_dp*x%d1val2
      q47 = q46*y%d1val2
      q48 = y%d1val1*z_d1x_d2y
      q49 = x%d1val2*x%d1val3
      q50 = x%d1val2*y%d1val3
      q51 = x%d1val3*y%d1val2
      q52 = y%d1val2*y%d1val3
      q53 = 2.0_dp*q38
      q54 = 2.0_dp*y%d1val3
      q55 = q45*x%d1val3
      q56 = q48*x%d1val3
      q57 = pow3(x%d1val2)
      q58 = pow3(y%d1val2)
      q59 = 3.0_dp*x%d2val2
      q60 = 3.0_dp*y%d2val2
      q61 = q8*z_d1x_d2y
      q62 = 3.0_dp*q61
      q63 = q7*z_d2x_d1y
      q64 = 3.0_dp*q63
      q65 = 2.0_dp*x%d1val2_d1val3
      q66 = q7*z_d3x
      q67 = q8*z_d3y
      q68 = 2.0_dp*y%d1val2
      q69 = q68*z_d2x_d1y
      q70 = q50*z_d1x_d2y
      q71 = q14*z_d3x
      q72 = q15*z_d1x_d2y
      q73 = q14*z_d2x_d1y
      q74 = q15*z_d3y
      q75 = q49*z_d2x_d1y
      q76 = q51*z_d1x_d2y
      q77 = 3.0_dp*x%d2val1_d1val3
      q78 = 3.0_dp*y%d2val1_d1val3
      q79 = x%d2val1*z_d2x
      q80 = 3.0_dp*x%d1val1_d1val3
      q81 = y%d2val1*z_d1x_d1y
      q82 = x%d2val1*z_d1x_d1y
      q83 = 3.0_dp*y%d1val1_d1val3
      q84 = y%d2val1*z_d2y
      q85 = x%d1val3*z_d4x
      q86 = y%d1val3*z_d3x_d1y
      q87 = x%d1val3*z_d1x_d3y
      q88 = y%d1val3*z_d4y
      q89 = 6.0_dp*x%d1val1_d1val3
      q90 = q31*x%d1val1
      q91 = q20*x%d1val3
      q92 = q20*y%d1val3
      q93 = q34*y%d1val1_d1val3
      q94 = 6.0_dp*y%d1val1
      q95 = q21*y%d1val3
      q96 = q1*z_d2x_d2y
      q97 = 3.0_dp*q96
      q98 = 3.0_dp*x%d1val1
      q99 = q1*z_d1x_d3y
      q100 = q99*y%d1val3
      q101 = 3.0_dp*y%d1val1
      q102 = q0*z_d3x_d1y
      q103 = q102*x%d1val3
      q104 = q0*z_d2x_d2y
      q105 = 2.0_dp*x%d1val1_d1val2_d1val3
      q106 = 2.0_dp*z_d2x
      q107 = q106*x%d1val1_d1val2
      q108 = y%d1val1_d1val3*z_d1x_d1y
      q109 = y%d1val1_d1val2*z_d1x_d1y
      q110 = 2.0_dp*z_d2y
      q111 = q110*y%d1val1_d1val2
      q112 = x%d2val1*z_d3x
      q113 = y%d2val1*z_d2x_d1y
      q114 = x%d2val1*z_d2x_d1y
      q115 = y%d2val1*z_d1x_d2y
      q116 = x%d2val1*z_d1x_d2y
      q117 = y%d2val1*z_d3y
      q118 = q26*x%d1val3
      q119 = q26*y%d1val3
      q120 = q39*x%d1val2
      q121 = q39*y%d1val2
      q122 = q46*y%d1val1_d1val3
      q123 = 2.0_dp*y%d1val1_d1val2
      q124 = q35*y%d1val1
      q125 = q35*y%d1val3
      q126 = q54*y%d1val1_d1val2
      q127 = q44*y%d1val1_d1val3
      q128 = q0*x%d1val2
      q129 = q1*y%d1val2
      q130 = 2.0_dp*y%d1val1
      q131 = q130*x%d1val1
      q132 = q49*z_d3x_d1y
      q133 = q131*z_d2x_d2y
      q134 = q131*z_d1x_d3y
      q135 = 2.0_dp*x%d1val1_d2val3
      q136 = 4.0_dp*x%d1val1_d1val3
      q137 = 2.0_dp*x%d2val1_d1val3
      q138 = pow2(x%d1val1_d1val3)
      q139 = pow2(y%d1val1_d1val3)
      q140 = q136*x%d1val3
      q141 = q136*y%d1val3
      q142 = 4.0_dp*y%d1val1_d1val3
      q143 = x%d1val1*x%d2val3
      q144 = 4.0_dp*y%d1val3
      q145 = q54*x%d1val3
      q146 = x%d1val3*z_d2x_d2y
      q147 = q146*x%d1val1
      q148 = 4.0_dp*q41
      q149 = q14*z_d4x
      q150 = q15*z_d2x_d2y
      q151 = q14*z_d2x_d2y
      q152 = q15*z_d4y
      q153 = q14*z_d3x_d1y
      q154 = q153*x%d1val1
      q155 = q15*z_d1x_d3y
      q156 = q155*x%d1val1
      q157 = q0*q86
      q158 = 2.0_dp*x%d1val3
      q159 = q54*q87
      q160 = x%d2val2*z_d2x
      q161 = y%d2val2*z_d1x_d1y
      q162 = x%d2val2*z_d1x_d1y
      q163 = y%d2val2*z_d2y
      q164 = y%d1val2_d1val3*z_d1x_d1y
      q165 = x%d1val3*x%d2val2
      q166 = x%d1val3*y%d2val2
      q167 = x%d2val2*y%d1val3
      q168 = y%d1val3*y%d2val2
      q169 = q65*x%d1val2
      q170 = q46*y%d1val2_d1val3
      q171 = q65*y%d1val2
      q172 = q68*y%d1val2_d1val3
      q173 = q49*z_d3x
      q174 = q26*z_d2x_d1y
      q175 = q52*z_d1x_d2y
      q176 = y%d1val2*z_d2x_d1y
      q177 = q46*z_d1x_d2y
      q178 = y%d1val1_d1val3*y%d1val2
      q179 = q52*z_d3y
      q180 = q7*x%d1val1
      q181 = q8*x%d1val1
      q182 = q181*z_d1x_d3y
      q183 = q7*y%d1val1
      q184 = q183*z_d3x_d1y
      q185 = q183*z_d2x_d2y
      q186 = q8*y%d1val1
      q187 = x%d1val1*y%d1val2
      q188 = 2.0_dp*q187
      q189 = q50*z_d2x_d2y
      q190 = q49*z_d2x_d2y
      q191 = 2.0_dp*q33
      q192 = q50*z_d1x_d3y
      q193 = x%d2val3*z_d2x
      q194 = y%d2val3*z_d1x_d1y
      q195 = x%d2val3*z_d1x_d1y
      q196 = y%d2val3*z_d2y
      q197 = q106*x%d1val2_d1val3
      q198 = q110*y%d1val2_d1val3
      q199 = x%d1val2*x%d2val3
      q200 = x%d1val2*y%d2val3
      q201 = x%d2val3*y%d1val2
      q202 = y%d1val2*y%d2val3
      q203 = q65*x%d1val3
      q204 = q65*y%d1val3
      q205 = 2.0_dp*y%d1val2_d1val3
      q206 = y%d1val3*z_d2x_d1y
      q207 = q39*z_d2x_d1y
      q208 = 2.0_dp*y%d1val1_d1val3
      q209 = x%d1val3*z_d1x_d2y
      q210 = q54*y%d1val2_d1val3
      q211 = x%d1val2*y%d1val1
      q212 = q49*q86
      q213 = 2.0_dp*q212
      q214 = q54*z_d2x_d2y
      q215 = q214*q51
      q216 = 2.0_dp*q41
      q217 = q52*q87
      q218 = 3.0_dp*x%d2val2_d1val3
      q219 = 3.0_dp*y%d2val2_d1val3
      q220 = 3.0_dp*x%d1val2_d1val3
      q221 = 3.0_dp*y%d1val2_d1val3
      q222 = 6.0_dp*x%d1val2
      q223 = q176*x%d1val2_d1val3
      q224 = q59*z_d2x_d1y
      q225 = q222*z_d1x_d2y
      q226 = y%d1val2*y%d1val2_d1val3
      q227 = 3.0_dp*q8
      q228 = q227*z_d1x_d3y
      q229 = 3.0_dp*q7
      q230 = q51*z_d3x_d1y
      q231 = q229*z_d2x_d2y
      q232 = 2.0_dp*x%d1val2_d2val3
      q233 = 4.0_dp*x%d1val2_d1val3
      q234 = 2.0_dp*x%d2val2_d1val3
      q235 = pow2(x%d1val2_d1val3)
      q236 = pow2(y%d1val2_d1val3)
      q237 = q233*z_d2x_d1y
      q238 = 4.0_dp*y%d1val2_d1val3
      q239 = q200*z_d1x_d2y
      q240 = q165*z_d2x_d1y
      q241 = q166*z_d1x_d2y
      q242 = 4.0_dp*q52
      q243 = q7*q86
      q244 = 3.0_dp*x%d2val1_d2val3
      q245 = 3.0_dp*y%d2val1_d2val3
      q246 = y%d2val1_d1val3*z_d1x_d1y
      q247 = 3.0_dp*x%d1val1_d2val3
      q248 = 2.0_dp*x%d3val1_d1val3
      q249 = 3.0_dp*y%d1val1_d2val3
      q250 = 6.0_dp*x%d2val1_d1val3
      q251 = y%d1val1_d1val3*z_d2y
      q252 = 6.0_dp*y%d2val1_d1val3
      q253 = x%d2val3*z_d4x
      q254 = y%d2val3*z_d3x_d1y
      q255 = x%d2val3*z_d1x_d3y
      q256 = y%d2val3*z_d4y
      q257 = q45*y%d1val1_d1val3
      q258 = 12.0_dp*x%d1val1_d1val3
      q259 = q250*x%d1val3
      q260 = q20*x%d2val3
      q261 = q20*y%d2val3
      q262 = q250*y%d1val3
      q263 = q21*x%d2val3
      q264 = q34*y%d1val1_d2val3
      q265 = q252*y%d1val3
      q266 = q21*y%d2val3
      q267 = q89*x%d1val3
      q268 = q89*y%d1val3
      q269 = q48*y%d1val1_d1val3
      q270 = 6.0_dp*y%d1val1_d1val3
      q271 = q270*x%d1val3
      q272 = q145*z_d2x_d1y
      q273 = q145*z_d1x_d2y
      q274 = q270*y%d1val3
      q275 = 6.0_dp*q138
      q276 = 6.0_dp*q139
      q277 = q40*y%d1val1*z_d3x_d1y
      q278 = q41*x%d1val1
      q279 = q278*z_d2x_d2y
      q280 = 6.0_dp*q40
      q281 = 12.0_dp*y%d1val1_d1val3
      q282 = q147*y%d1val1
      q283 = 6.0_dp*y%d2val1
      q284 = q278*z_d1x_d3y
      q285 = 6.0_dp*q41
      q286 = q14*z_d5x
      q287 = q15*z_d3x_d2y
      q288 = q14*z_d2x_d3y
      q289 = q15*z_d5y
      q290 = q20*x%d1val1
      q291 = q0*q85
      q292 = q145*z_d4x_d1y
      q293 = q1*q270
      q294 = q145*z_d1x_d4y
      q295 = q20*y%d1val1
      q296 = q21*y%d1val1
      q297 = q1*z_d2x_d3y
      q298 = q0*z_d3x_d2y
      q299 = q14*z_d3x_d2y
      q300 = q1*q98
      q301 = q15*z_d1x_d4y
      q302 = q14*z_d4x_d1y
      q303 = q0*q101
      q304 = q15*z_d2x_d3y
      q305 = 2.0_dp*x%d1val1_d1val2_d2val3
      q306 = q26*z_d1x_d1y
      q307 = x%d1val1_d1val2_d1val3*z_d2x
      q308 = 4.0_dp*x%d1val1_d1val2_d1val3
      q309 = y%d1val1_d1val2_d1val3*z_d1x_d1y
      q310 = 2.0_dp*x%d2val1_d1val2_d1val3
      q311 = 4.0_dp*y%d1val1_d1val2_d1val3
      q312 = q26*x%d2val3
      q313 = q26*y%d2val3
      q314 = q308*x%d1val3
      q315 = q308*y%d1val3
      q316 = q136*x%d1val2_d1val3
      q317 = q136*y%d1val2_d1val3
      q318 = q135*x%d1val2
      q319 = q135*y%d1val2
      q320 = q46*y%d1val1_d2val3
      q321 = q123*x%d2val3
      q322 = q311*y%d1val3
      q323 = x%d1val1_d1val2*z_d3x
      q324 = x%d1val1_d1val2*z_d2x_d1y
      q325 = q324*x%d1val3
      q326 = q142*y%d1val3
      q327 = x%d1val1_d1val2*z_d1x_d2y
      q328 = x%d1val2*z_d2x_d1y
      q329 = y%d1val1_d1val2*z_d1x_d2y
      q330 = q136*z_d1x_d2y
      q331 = 2.0_dp*y%d2val1_d1val3
      q332 = q137*z_d2x_d1y
      q333 = q158*y%d1val2_d1val3
      q334 = q209*y%d1val1_d1val2
      q335 = q44*q68
      q336 = y%d1val1_d1val2*z_d3y
      q337 = x%d1val2*x%d2val1
      q338 = x%d1val2*y%d2val1
      q339 = x%d2val1*y%d1val2
      q340 = y%d1val2*y%d2val1
      q341 = 4.0_dp*x%d1val1_d1val2
      q342 = q136*q30
      q343 = q136*x%d1val1
      q344 = q52*z_d2x_d2y
      q345 = q142*x%d1val1
      q346 = q199*z_d3x_d1y
      q347 = q144*y%d1val1_d1val2
      q348 = q51*z_d2x_d2y
      q349 = q52*z_d1x_d3y
      q350 = q136*y%d1val1
      q351 = q142*y%d1val1
      q352 = q142*q33
      q353 = q26*x%d1val1
      q354 = q26*y%d1val1
      q355 = q1*q205
      q356 = q123*y%d1val1
      q357 = q148*x%d1val1
      q358 = q49*z_d3x_d2y
      q359 = q0*y%d1val2
      q360 = q1*x%d1val2
      q361 = q130*q30
      q362 = q130*q187
      q363 = q49*q54
      q364 = q51*q54
      q365 = y%d2val2_d1val3*z_d1x_d1y
      q366 = 2.0_dp*x%d1val1_d2val2_d1val3
      q367 = y%d1val2_d1val3*z_d2y
      q368 = x%d2val2*x%d2val3
      q369 = x%d2val2*y%d2val3
      q370 = x%d2val3*y%d2val2
      q371 = y%d2val2*y%d2val3
      q372 = q232*x%d1val2
      q373 = q46*y%d1val2_d2val3
      q374 = q233*y%d1val2_d1val3
      q375 = q232*y%d1val2
      q376 = q234*x%d1val3
      q377 = 2.0_dp*y%d2val2_d1val3
      q378 = q234*y%d1val3
      q379 = q199*z_d3x
      q380 = q202*z_d1x_d2y
      q381 = q144*y%d1val2_d1val3
      q382 = q308*z_d2x_d1y
      q383 = q165*z_d3x
      q384 = q168*z_d1x_d2y
      q385 = q199*z_d2x_d1y
      q386 = q167*z_d1x_d2y
      q387 = q201*z_d1x_d2y
      q388 = q202*z_d3y
      q389 = q168*z_d3y
      q390 = 2.0_dp*q235
      q391 = 2.0_dp*q236
      q392 = x%d1val1*x%d2val2
      q393 = x%d1val1*y%d2val2
      q394 = x%d2val2*y%d1val1
      q395 = y%d1val1*y%d2val2
      q396 = q233*q30
      q397 = q238*x%d1val1
      q398 = q233*x%d1val1
      q399 = q136*y%d1val2
      q400 = q233*y%d1val1
      q401 = q238*y%d1val1
      q402 = q142*y%d1val2
      q403 = q199*z_d2x_d2y
      q404 = q238*q33
      q405 = q26*x%d1val2
      q406 = q26*y%d1val2
      q407 = q7*q85
      q408 = q39*q8
      q409 = y%d1val3*z_d1x_d3y
      q410 = q46*y%d1val1_d1val2
      q411 = q7*y%d1val1_d1val3
      q412 = q208*q8
      q413 = q68*y%d1val1_d1val2
      q414 = q49*z_d2x_d3y
      q415 = q30*q68
      q416 = q211*q68
      q417 = 3.0_dp*x%d2val2_d2val3
      q418 = 3.0_dp*y%d2val2_d2val3
      q419 = 6.0_dp*x%d1val2_d1val3
      q420 = 3.0_dp*x%d1val2_d2val3
      q421 = 2.0_dp*x%d3val2_d1val3
      q422 = 3.0_dp*y%d1val2_d2val3
      q423 = 6.0_dp*x%d2val2_d1val3
      q424 = 6.0_dp*y%d2val2_d1val3
      q425 = 12.0_dp*x%d1val2_d1val3
      q426 = q423*z_d2x_d1y
      q427 = q419*z_d2x_d1y
      q428 = 6.0_dp*y%d1val2_d1val3
      q429 = 6.0_dp*y%d1val2
      q430 = q425*y%d1val2
      q431 = 12.0_dp*q226
      q432 = q59*x%d1val2
      q433 = q60*x%d1val2
      q434 = q428*q7
      q435 = q419*q8
      q436 = q428*q8
      q437 = q59*y%d1val2
      q438 = q60*y%d1val2
      q439 = 6.0_dp*y%d1val3
      q440 = q227*x%d1val2
      q441 = q229*y%d1val2
      binary%val = z_val
      binary%d1val1 = x%d1val1*z_d1x + y%d1val1*z_d1y
      binary%d1val2 = x%d1val2*z_d1x + y%d1val2*z_d1y
      binary%d1val3 = x%d1val3*z_d1x + y%d1val3*z_d1y
      binary%d2val1 = q0*z_d2x + q1*z_d2y + q3*y%d1val1 + x%d2val1*z_d1x + y%d2val1*z_d1y
      binary%d1val1_d1val2 = q2*y%d1val2 + q4*x%d1val2 + q5*x%d1val2 + q6*y%d1val2 + x%d1val1_d1val2*z_d1x + y%d1val1_d1val2*z_d1y
      binary%d1val1_d1val3 = q2*y%d1val3 + q4*x%d1val3 + q5*x%d1val3 + q6*y%d1val3 + x%d1val1_d1val3*z_d1x + y%d1val1_d1val3*z_d1y
      binary%d2val2 = q10*y%d1val2 + q7*z_d2x + q8*z_d2y + x%d2val2*z_d1x + y%d2val2*z_d1y
      binary%d1val2_d1val3 = q11*x%d1val3 + q12*x%d1val3 + q13*y%d1val3 + q9*y%d1val3 + x%d1val2_d1val3*z_d1x + y%d1val2_d1val3*z_d1y
      binary%d2val3 = q14*z_d2x + q15*z_d2y + q17*y%d1val3 + x%d2val3*z_d1x + y%d2val3*z_d1y
      binary%d3val1 = q18*z_d3x + q19*z_d3y + q2*q21 + q20*q4 + q20*q5 + q21*q6 + q23*x%d1val1 + q25*y%d1val1 + x%d3val1*z_d1x + y%d3val1*z_d1y
      binary%d2val1_d1val2 = q11*x%d2val1 + q12*x%d2val1 + q13*y%d2val1 + q22*x%d1val2 + q24*y%d1val2 + q26*q4 + q26*q5 + q27*y%d1val1_d1val2 + q28*x%d1val2 + q29*y%d1val2 + q3*y%d1val1_d1val2 + q30*q32 + q33*q35 + q9*y%d2val1 + x%d2val1_d1val2*z_d1x + y%d2val1_d1val2*z_d1y
      binary%d2val1_d1val3 = q16*y%d2val1 + q22*x%d1val3 + q24*y%d1val3 + q27*y%d1val1_d1val3 + q28*x%d1val3 + q29*y%d1val3 + q3*y%d1val1_d1val3 + q32*q40 + q35*q41 + q36*x%d2val1 + q37*x%d2val1 + q38*y%d2val1 + q39*q4 + q39*q5 + x%d2val1_d1val3*z_d1x + y%d2val1_d1val3*z_d1y
      binary%d1val1_d2val2 = q10*y%d1val1_d1val2 + q11*q26 + q12*q26 + q2*y%d2val2 + q31*q7 + q34*q8 + q4*x%d2val2 + q42*y%d1val1_d1val2 + q43*q7 + q44*q8 + q45*q47 + q47*q48 + q5*x%d2val2 + q6*y%d2val2 + x%d1val1_d2val2*z_d1x + y%d1val1_d2val2*z_d1y
      binary%d1val1_d1val2_d1val3 = q11*x%d1val1_d1val3 + q12*x%d1val1_d1val3 + q13*y%d1val1_d1val3 + q16*y%d1val1_d1val2 + q2*y%d1val2_d1val3 + q31*q49 + q34*q52 + q36*x%d1val1_d1val2 + q37*x%d1val1_d1val2 + q38*y%d1val1_d1val2 + q4*x%d1val2_d1val3 + q43*q49 + q44*q52 + q45*q50 + q45*q51 + q48*q50 + q48*q51 + q5*x%d1val2_d1val3 + q6*y%d1val2_d1val3 + q9*y%d1val1_d1val3 + x%d1val1_d1val2_d1val3*z_d1x + y%d1val1_d1val2_d1val3*z_d1y
      binary%d1val1_d2val3 = q14*q31 + q14*q43 + q15*q34 + q15*q44 + q17*y%d1val1_d1val3 + q2*y%d2val3 + q36*q39 + q37*q39 + q4*x%d2val3 + q5*x%d2val3 + q53*y%d1val1_d1val3 + q54*q55 + q54*q56 + q6*y%d2val3 + x%d1val1_d2val3*z_d1x + y%d1val1_d2val3*z_d1y
      binary%d3val2 = q11*q59 + q12*q59 + q13*q60 + q57*z_d3x + q58*z_d3y + q60*q9 + q62*x%d1val2 + q64*y%d1val2 + x%d3val2*z_d1x + y%d3val2*z_d1y
      binary%d2val2_d1val3 = q10*y%d1val2_d1val3 + q11*q65 + q12*q65 + q16*y%d2val2 + q36*x%d2val2 + q37*x%d2val2 + q38*y%d2val2 + q42*y%d1val2_d1val3 + q49*q69 + q61*x%d1val3 + q63*y%d1val3 + q66*x%d1val3 + q67*y%d1val3 + q68*q70 + x%d2val2_d1val3*z_d1x + y%d2val2_d1val3*z_d1y
      binary%d1val2_d2val3 = q11*x%d2val3 + q12*x%d2val3 + q13*y%d2val3 + q17*y%d1val2_d1val3 + q36*q65 + q37*q65 + q53*y%d1val2_d1val3 + q54*q75 + q54*q76 + q71*x%d1val2 + q72*x%d1val2 + q73*y%d1val2 + q74*y%d1val2 + q9*y%d2val3 + x%d1val2_d2val3*z_d1x + y%d1val2_d2val3*z_d1y
      binary%d3val1_d1val3 = 3.0_dp*q104*q41 + q100*q98 + q101*q103 + q16*y%d3val1 + q18*q85 + q18*q86 + q19*q87 + q19*q88 + q2*q78 + q21*q55 + q21*q56 + q23*x%d1val1_d1val3 + q25*y%d1val1_d1val3 + q28*q80 + q29*q83 + q31*q91 + q34*q95 + q36*x%d3val1 + q37*x%d3val1 + q38*y%d3val1 + q4*q77 + q40*q97 + q43*q91 + q44*q95 + q45*q92 + q48*q92 + q5*q77 + q6*q78 + q79*q80 + q80*q81 + q82*q83 + q83*q84 + q89*q90 + q93*q94 + x%d3val1_d1val3*z_d1x + y%d3val1_d1val3*z_d1y
      binary%d2val1_d1val2_d1val3 = q102*q51 + q104*q52 + q105*q4 + q105*q5 + q107*x%d1val1_d1val3 + q108*q26 + q109*q39 + q11*x%d2val1_d1val3 + q111*y%d1val1_d1val3 + q112*q49 + q113*q49 + q114*q50 + q114*q51 + q115*q50 + q115*q51 + q116*q52 + q117*q52 + q118*q31 + q118*q43 + q119*q45 + q119*q48 + q12*x%d2val1_d1val3 + q120*q31 + q120*q43 + q121*q45 + q121*q48 + q122*q45 + q122*q48 + q123*q55 + q123*q56 + q124*y%d1val2_d1val3 + q125*y%d1val1_d1val2 + q126*q44 + q127*q68 + q128*q85 + q128*q86 + q129*q87 + q129*q88 + q13*y%d2val1_d1val3 + q131*q132 + q133*q50 + q133*q51 + q134*q52 + q16*y%d2val1_d1val2 + q22*x%d1val2_d1val3 + q24*y%d1val2_d1val3 + q27*y%d1val1_d1val2_d1val3 + q28*x%d1val2_d1val3 + q29*y%d1val2_d1val3 + q3*y%d1val1_d1val2_d1val3 + q36*x%d2val1_d1val2 + q37*x%d2val1_d1val2 + q38*y%d2val1_d1val2 + q49*q96 + q50*q99 + q65*q90 + q68*q93 + q79*x%d1val2_d1val3 + q81*x%d1val2_d1val3 + q82*y%d1val2_d1val3 + q84*y%d1val2_d1val3 + q9*y%d2val1_d1val3 + x%d2val1_d1val2_d1val3*z_d1x + y%d2val1_d1val2_d1val3*z_d1y
      binary%d2val1_d2val3 = q0*q149 + q0*q150 + q1*q151 + q1*q152 + q1*q159 + q106*q138 + q108*q136 + q110*q139 + q114*q145 + q115*q145 + q124*y%d2val3 + q127*q144 + q130*q154 + q130*q156 + q135*q4 + q135*q5 + q137*q36 + q137*q37 + q140*q31 + q140*q43 + q141*q45 + q141*q48 + q142*q55 + q142*q56 + q143*q32 + q144*q93 + q147*q148 + q157*q158 + q17*y%d2val1_d1val3 + q22*x%d2val3 + q24*y%d2val3 + q27*y%d1val1_d2val3 + q28*x%d2val3 + q29*y%d2val3 + q3*y%d1val1_d2val3 + q53*y%d2val1_d1val3 + q71*x%d2val1 + q72*x%d2val1 + q73*y%d2val1 + q74*y%d2val1 + q79*x%d2val3 + q81*x%d2val3 + q82*y%d2val3 + q84*y%d2val3 + x%d2val1_d2val3*z_d1x + y%d2val1_d2val3*z_d1y
      binary%d1val1_d2val2_d1val3 = q10*y%d1val1_d1val2_d1val3 + q105*q11 + q105*q12 + q107*x%d1val2_d1val3 + q109*q65 + q111*y%d1val2_d1val3 + q120*q176 + q123*q179 + q123*q70 + q123*q75 + q123*q76 + q132*q188 + q146*q181 + q16*y%d1val1_d2val2 + q160*x%d1val1_d1val3 + q161*x%d1val1_d1val3 + q162*y%d1val1_d1val3 + q163*y%d1val1_d1val3 + q164*q26 + q165*q31 + q165*q43 + q166*q45 + q166*q48 + q167*q45 + q167*q48 + q168*q34 + q168*q44 + q169*q31 + q169*q43 + q170*q45 + q170*q48 + q171*q45 + q171*q48 + q172*q34 + q172*q44 + q173*q26 + q174*q50 + q174*q51 + q175*q26 + q177*q178 + q180*q85 + q180*q86 + q182*y%d1val3 + q184*x%d1val3 + q185*y%d1val3 + q186*q87 + q186*q88 + q188*q189 + q190*q191 + q191*q192 + q2*y%d2val2_d1val3 + q36*x%d1val1_d2val2 + q37*x%d1val1_d2val2 + q38*y%d1val1_d2val2 + q4*x%d2val2_d1val3 + q42*y%d1val1_d1val2_d1val3 + q5*x%d2val2_d1val3 + q6*y%d2val2_d1val3 + q61*x%d1val1_d1val3 + q63*y%d1val1_d1val3 + q66*x%d1val1_d1val3 + q67*y%d1val1_d1val3 + x%d1val1_d2val2_d1val3*z_d1x + y%d1val1_d2val2_d1val3*z_d1y
      binary%d1val1_d1val2_d2val3 = q105*q36 + q105*q37 + q108*q65 + q11*x%d1val1_d2val3 + q118*q206 + q12*x%d1val1_d2val3 + q125*y%d1val2_d1val3 + q126*q209 + q13*y%d1val1_d2val3 + q130*q217 + q149*q30 + q150*q30 + q151*q33 + q152*q33 + q153*q187 + q153*q211 + q155*q187 + q155*q211 + q164*q39 + q17*y%d1val1_d1val2_d1val3 + q173*q39 + q175*q39 + q179*q208 + q190*q216 + q193*x%d1val1_d1val2 + q194*x%d1val1_d1val2 + q195*y%d1val1_d1val2 + q196*y%d1val1_d1val2 + q197*x%d1val1_d1val3 + q198*y%d1val1_d1val3 + q199*q31 + q199*q43 + q2*y%d1val2_d2val3 + q200*q45 + q200*q48 + q201*q45 + q201*q48 + q202*q34 + q202*q44 + q203*q31 + q203*q43 + q204*q45 + q204*q48 + q205*q55 + q205*q56 + q207*q50 + q207*q51 + q208*q70 + q208*q75 + q208*q76 + q210*q44 + q213*x%d1val1 + q215*x%d1val1 + q4*x%d1val2_d2val3 + q5*x%d1val2_d2val3 + q53*y%d1val1_d1val2_d1val3 + q6*y%d1val2_d2val3 + q71*x%d1val1_d1val2 + q72*x%d1val1_d1val2 + q73*y%d1val1_d1val2 + q74*y%d1val1_d1val2 + q9*y%d1val1_d2val3 + x%d1val1_d1val2_d2val3*z_d1x + y%d1val1_d1val2_d2val3*z_d1y
      binary%d3val2_d1val3 = q11*q218 + q12*q218 + q13*q219 + q16*y%d3val2 + q160*q220 + q161*q220 + q162*q221 + q163*q221 + q173*q59 + q175*q59 + q179*q60 + q190*q227 + q219*q9 + q220*q66 + q221*q67 + q222*q223 + q224*q50 + q224*q51 + q225*q226 + q228*q50 + q229*q230 + q231*q52 + q36*x%d3val2 + q37*x%d3val2 + q38*y%d3val2 + q57*q85 + q57*q86 + q58*q87 + q58*q88 + q60*q70 + q60*q75 + q60*q76 + q62*x%d1val2_d1val3 + q64*y%d1val2_d1val3 + x%d3val2_d1val3*z_d1x + y%d3val2_d1val3*z_d1y
      binary%d2val2_d2val3 = q10*y%d1val2_d2val3 + q106*q235 + q11*q232 + q110*q236 + q12*q232 + q149*q7 + q150*q7 + q151*q8 + q152*q8 + q153*q47 + q155*q47 + q158*q243 + q159*q8 + q160*x%d2val3 + q161*x%d2val3 + q162*y%d2val3 + q163*y%d2val3 + q164*q233 + q17*y%d2val2_d1val3 + q173*q233 + q175*q233 + q179*q238 + q190*q242 + q199*q69 + q234*q36 + q234*q37 + q237*q50 + q237*q51 + q238*q70 + q238*q75 + q238*q76 + q239*q68 + q240*q54 + q241*q54 + q42*y%d1val2_d2val3 + q53*y%d2val2_d1val3 + q61*x%d2val3 + q63*y%d2val3 + q66*x%d2val3 + q67*y%d2val3 + q71*x%d2val2 + q72*x%d2val2 + q73*y%d2val2 + q74*y%d2val2 + x%d2val2_d2val3*z_d1x + y%d2val2_d2val3*z_d1y
      binary%d3val1_d2val3 = 6.0_dp*q90*x%d1val1_d2val3 + q100*q89 + q101*q102*x%d2val3 + q101*q104*y%d2val3 + q103*q270 + q104*q274 + q108*q250 + q112*q267 + q113*q267 + q114*q268 + q114*q271 + q115*q268 + q115*q271 + q116*q274 + q117*q274 + q143*q97 + q146*q285*x%d2val1 + q147*q283*y%d1val3 + q149*q290 + q150*q290 + q151*q296 + q152*q296 + q153*q295 + q154*q21 + q155*q295 + q156*q21 + q157*q89 + q17*y%d3val1_d1val3 + q18*q253 + q18*q254 + q18*q286 + q18*q287 + q18*q292 + q19*q255 + q19*q256 + q19*q288 + q19*q289 + q19*q294 + q193*x%d3val1 + q194*x%d3val1 + q195*y%d3val1 + q196*y%d3val1 + q2*q245 + q23*x%d1val1_d2val3 + q244*q4 + q244*q5 + q245*q6 + q246*q89 + q247*q28 + q247*q79 + q247*q81 + q248*q36 + q248*q37 + q249*q29 + q249*q82 + q249*q84 + q25*y%d1val1_d2val3 + q251*q252 + q252*q55 + q252*q56 + q257*q258 + q258*q269 + q258*q277 + q258*q279 + q259*q31 + q259*q43 + q260*q31 + q260*q43 + q261*q45 + q261*q48 + q262*q45 + q262*q48 + q263*q45 + q263*q48 + q264*q94 + q265*q34 + q265*q44 + q266*q34 + q266*q44 + q267*q96 + q272*x%d3val1 + q273*y%d3val1 + q275*q31 + q275*q43 + q276*q34 + q276*q44 + q280*q297*y%d1val3 + q280*q86*x%d2val1 + q281*q282 + q281*q284 + q283*q41*q87 + q285*q298*x%d1val3 + q291*q89 + q293*q87 + q293*q88 + q299*q300 + q300*q301 + q302*q303 + q303*q304 + q53*y%d3val1_d1val3 + q71*x%d3val1 + q72*x%d3val1 + q73*y%d3val1 + q74*y%d3val1 + q89*x%d2val1_d1val3*z_d2x + q98*q99*y%d2val3 + x%d3val1_d2val3*z_d1x + y%d3val1_d2val3*z_d1y
      binary%d2val1_d1val2_d2val3 = 2.0_dp*q217*y%d2val1 + q0*q363*z_d4x_d1y + q1*q364*z_d1x_d4y + q100*q65 + q102*q201 + q103*q205 + q104*q202 + q104*q210 + q107*x%d1val1_d2val3 + q108*q308 + q109*q135 + q11*x%d2val1_d2val3 + q111*y%d1val1_d2val3 + q112*q199 + q112*q203 + q113*q199 + q113*q203 + q114*q200 + q114*q201 + q114*q204 + q114*q333 + q115*q200 + q115*q201 + q115*q204 + q115*q333 + q116*q202 + q116*q210 + q117*q202 + q117*q210 + q12*x%d2val1_d2val3 + q123*q154 + q123*q156 + q123*q44*y%d2val3 + q124*y%d1val2_d2val3 + q127*q238 + q128*q253 + q128*q254 + q128*q286 + q128*q287 + q129*q255 + q129*q256 + q129*q288 + q129*q289 + q13*y%d2val1_d2val3 + q131*q346 + q132*q345 + q132*q350 + q133*q200 + q133*q201 + q134*q202 + q136*q307 + q136*q309 + q136*q328*y%d1val1_d1val3 + q137*q164 + q137*q173 + q137*q175 + q138*q46*z_d3x + q138*q69 + q139*q177 + q139*q68*z_d3y + q140*q323 + q140*y%d1val1_d1val2*z_d2x_d1y + q141*q324 + q141*q329 + q142*q325 + q142*q334 + q146*q148*x%d1val1_d1val2 + q147*q347 + q148*q87*y%d1val1_d1val2 + q149*q337 + q149*q353 + q150*q337 + q150*q353 + q151*q340 + q151*q356 + q152*q340 + q152*q356 + q153*q338 + q153*q339 + q153*q354 + q155*q338 + q155*q339 + q155*q354 + q157*q65 + q17*y%d2val1_d1val2_d1val3 + q178*q330 + q179*q331 + q189*q345 + q189*q350 + q190*q351 + q190*q54*y%d2val1 + q192*q351 + q193*x%d2val1_d1val2 + q194*x%d2val1_d1val2 + q195*y%d2val1_d1val2 + q196*y%d2val1_d1val2 + q197*x%d2val1_d1val3 + q198*y%d2val1_d1val3 + q199*q96 + q200*q99 + q203*q96 + q213*x%d2val1 + q215*x%d2val1 + q22*x%d1val2_d2val3 + q230*q343 + q232*q90 + q233*q257 + q233*q269 + q233*q277 + q233*q279 + q238*q282 + q238*q284 + q238*q93 + q24*y%d1val2_d2val3 + q246*q65 + q251*q311 + q264*q68 + q27*y%d1val1_d1val2_d2val3 + q272*x%d2val1_d1val2 + q273*y%d2val1_d1val2 + q28*x%d1val2_d2val3 + q29*y%d1val2_d2val3 + q291*q65 + q297*q363 + q298*q364 + q299*q360 + q299*q362 + q3*y%d1val1_d1val2_d2val3 + q301*q360 + q301*q362 + q302*q359 + q302*q361 + q304*q359 + q304*q361 + q305*q4 + q305*q5 + q306*y%d1val1_d2val3 + q31*q312 + q31*q314 + q31*q316 + q31*q318 + q310*q36 + q310*q37 + q311*q55 + q311*q56 + q312*q43 + q313*q45 + q313*q48 + q314*q43 + q315*q45 + q315*q48 + q316*q43 + q317*q45 + q317*q48 + q318*q43 + q319*q45 + q319*q48 + q320*q45 + q320*q48 + q321*q45 + q321*q48 + q322*q34 + q322*q44 + q326*q327 + q326*q336 + q331*q70 + q331*q75 + q331*q76 + q332*q50 + q332*q51 + q335*y%d1val1_d2val3 + q341*q40*q86 + q342*q85 + q342*q86 + q343*q344 + q345*q348 + q345*q349 + q348*q350 + q349*q350 + q35*y%d1val1_d1val2*y%d2val3 + q352*q87 + q352*q88 + q355*q87 + q355*q88 + q357*q358 + q357*q51*z_d2x_d3y + q53*y%d2val1_d1val2_d1val3 + q71*x%d2val1_d1val2 + q72*x%d2val1_d1val2 + q73*y%d2val1_d1val2 + q74*y%d2val1_d1val2 + q79*x%d1val2_d2val3 + q81*x%d1val2_d2val3 + q82*y%d1val2_d2val3 + q84*y%d1val2_d2val3 + q9*y%d2val1_d2val3 + x%d2val1_d1val2_d2val3*z_d1x + y%d2val1_d1val2_d2val3*z_d1y
      binary%d1val1_d2val2_d2val3 = 2.0_dp*q165*q86*x%d1val1 + 4.0_dp*q217*y%d1val1_d1val2 + q10*y%d1val1_d1val2_d2val3 + q106*x%d1val1_d1val3*x%d2val2_d1val3 + q107*x%d1val2_d2val3 + q108*q234 + q109*q232 + q11*q305 + q110*y%d1val1_d1val3*y%d2val2_d1val3 + q111*y%d1val2_d2val3 + q12*q305 + q123*q239 + q123*q385 + q123*q387 + q123*q388 + q125*y%d2val2_d1val3 + q130*q168*q87 + q132*q397 + q132*q399 + q132*q400 + q136*q223 + q142*q226*z_d3y + q142*x%d1val2*y%d1val2_d1val3*z_d1x_d2y + q144*q348*x%d1val1_d1val2 + q145*q181*z_d2x_d3y + q145*q183*z_d3x_d2y + q146*q408 + q149*q392 + q149*q405 + q150*q392 + q150*q405 + q151*q395 + q151*q413 + q152*q395 + q152*q413 + q153*q393 + q153*q394 + q153*q406 + q153*q410 + q155*q393 + q155*q394 + q155*q406 + q155*q410 + q158*q411*z_d3x_d1y + q160*x%d1val1_d2val3 + q161*x%d1val1_d2val3 + q162*y%d1val1_d2val3 + q163*y%d1val1_d2val3 + q164*q308 + q165*q216*z_d2x_d2y + q166*q207 + q166*q214*x%d1val1 + q167*q207 + q17*y%d1val1_d2val2_d1val3 + q173*q308 + q174*q200 + q174*q201 + q175*q308 + q176*q318 + q177*y%d1val1_d2val3*y%d1val2 + q178*q233*z_d1x_d2y + q179*q311 + q180*q253 + q180*q254 + q180*q286 + q180*q287 + q180*q292 + q181*q299 + q181*q301 + q181*x%d2val3*z_d2x_d2y + q182*y%d2val3 + q183*q302 + q183*q304 + q184*x%d2val3 + q185*y%d2val3 + q186*q255 + q186*q256 + q186*q288 + q186*q289 + q186*q294 + q188*q200*z_d2x_d2y + q188*q346 + q189*q397 + q189*q399 + q189*q400 + q190*q347 + q190*q401 + q190*q402 + q191*q200*z_d1x_d3y + q191*q403 + q192*q401 + q192*q402 + q193*x%d1val1_d2val2 + q194*x%d1val1_d2val2 + q195*y%d1val1_d2val2 + q196*y%d1val1_d2val2 + q2*y%d2val2_d2val3 + q206*q233*x%d1val1_d1val2 + q208*q240 + q208*q241 + q208*q386 + q208*q389 + q212*q341 + q214*q411 + q226*q330 + q230*q398 + q233*q307 + q233*q309 + q233*q323*x%d1val3 + q233*q329*y%d1val3 + q237*x%d1val2*y%d1val1_d1val3 + q237*x%d1val3*y%d1val1_d1val2 + q238*q325 + q238*q334 + q242*q358*x%d1val1 + q242*q414*y%d1val1 + q243*q39 + q26*q379 + q26*q380 + q272*x%d1val1_d2val2 + q273*y%d1val1_d2val2 + q299*q416 + q301*q416 + q302*q415 + q304*q415 + q306*y%d1val2_d2val3 + q31*q368 + q31*q372 + q31*q376 + q31*q390 + q311*q367 + q311*q70 + q311*q75 + q311*q76 + q316*x%d1val2*z_d3x + q317*q328 + q327*q381 + q335*y%d1val2_d2val3 + q336*q381 + q34*q371 + q34*q391 + q34*q68*y%d1val2_d2val3 + q344*q398 + q348*q397 + q348*q400 + q349*q397 + q349*q400 + q36*q366 + q365*q39 + q366*q37 + q368*q43 + q369*q45 + q369*q48 + q370*q45 + q370*q48 + q371*q44 + q372*q43 + q373*q45 + q373*q48 + q374*q45 + q374*q48 + q375*q45 + q375*q48 + q376*q43 + q377*q55 + q377*q56 + q378*q45 + q378*q48 + q382*q50 + q382*q51 + q383*q39 + q384*q39 + q39*q407 + q390*q43 + q391*q44 + q396*q85 + q396*q86 + q4*x%d2val2_d2val3 + q404*q87 + q404*q88 + q408*q409 + q412*q87 + q412*q88 + q42*y%d1val1_d1val2_d2val3 + q44*q54*y%d2val2_d1val3 + q5*x%d2val2_d2val3 + q53*y%d1val1_d2val2_d1val3 + q6*y%d2val2_d2val3 + q61*x%d1val1_d2val3 + q63*y%d1val1_d2val3 + q66*x%d1val1_d2val3 + q67*y%d1val1_d2val3 + q71*x%d1val1_d2val2 + q72*x%d1val1_d2val2 + q73*y%d1val1_d2val2 + q74*y%d1val1_d2val2 + x%d1val1_d2val2_d2val3*z_d1x + y%d1val1_d2val2_d2val3*z_d1y
      binary%d3val2_d2val3 = 6.0_dp*q167*q348 + 6.0_dp*q168*q190 + 6.0_dp*q212*x%d2val2 + 6.0_dp*q217*y%d2val2 + q11*q417 + q12*q417 + q13*q418 + q132*q430 + q146*q435 + q149*q432 + q150*q432 + q151*q438 + q152*q438 + q153*q433 + q153*q437 + q155*q433 + q155*q437 + q160*q420 + q161*q420 + q162*q422 + q163*q422 + q164*q423 + q166*q427 + q167*q427 + q17*y%d3val2_d1val3 + q173*q423 + q175*q423 + q176*q222*x%d1val2_d2val3 + q179*q424 + q189*q430 + q190*q431 + q192*q431 + q193*x%d3val2 + q194*x%d3val2 + q195*y%d3val2 + q196*y%d3val2 + q200*q224 + q200*q228 + q201*q224 + q201*q229*z_d3x_d1y + q202*q231 + q222*q235*z_d3x + q225*q236 + q225*y%d1val2*y%d1val2_d2val3 + q226*q425*z_d1x_d2y + q227*q403 + q235*q429*z_d2x_d1y + q236*q429*z_d3y + q239*q60 + q240*q428 + q241*q428 + q243*q419 + q253*q57 + q254*q57 + q255*q58 + q256*q58 + q272*x%d3val2 + q273*y%d3val2 + q286*q57 + q287*q57 + q288*q58 + q289*q58 + q292*q57 + q294*q58 + q299*q440 + q301*q440 + q302*q441 + q304*q441 + q328*q425*y%d1val2_d1val3 + q36*q421 + q365*q419 + q367*q424 + q37*q421 + q379*q59 + q380*q59 + q383*q419 + q384*q419 + q385*q60 + q386*q428 + q387*q60 + q388*q60 + q389*q428 + q407*q419 + q409*q435 + q414*q439*q8 + q418*q9 + q419*x%d2val2_d1val3*z_d2x + q420*q66 + q422*q67 + q424*q70 + q424*q75 + q424*q76 + q426*q50 + q426*q51 + q434*x%d1val3*z_d3x_d1y + q434*y%d1val3*z_d2x_d2y + q436*q87 + q436*q88 + q439*q51*q7*z_d3x_d2y + q53*y%d3val2_d1val3 + q62*x%d1val2_d2val3 + q64*y%d1val2_d2val3 + q71*x%d3val2 + q72*x%d3val2 + q73*y%d3val2 + q74*y%d3val2 + x%d3val2_d2val3*z_d1x + y%d3val2_d2val3*z_d1y
   end function make_binary_operator
   
   function unary_minus_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      unary%val = -x%val
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d1val3 = -x%d1val3
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d1val1_d1val3 = -x%d1val1_d1val3
      unary%d2val2 = -x%d2val2
      unary%d1val2_d1val3 = -x%d1val2_d1val3
      unary%d2val3 = -x%d2val3
      unary%d3val1 = -x%d3val1
      unary%d2val1_d1val2 = -x%d2val1_d1val2
      unary%d2val1_d1val3 = -x%d2val1_d1val3
      unary%d1val1_d2val2 = -x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = -x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = -x%d1val1_d2val3
      unary%d3val2 = -x%d3val2
      unary%d2val2_d1val3 = -x%d2val2_d1val3
      unary%d1val2_d2val3 = -x%d1val2_d2val3
      unary%d3val1_d1val3 = -x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = -x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = -x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = -x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = -x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = -x%d3val2_d1val3
      unary%d2val2_d2val3 = -x%d2val2_d2val3
      unary%d3val1_d2val3 = -x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = -x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = -x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = -x%d3val2_d2val3
   end function unary_minus_self
   
   function exp_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = exp(x%val)
      q1 = q0*x%d1val1
      q2 = q0*x%d1val2
      q3 = q0*x%d1val3
      q4 = pow2(x%d1val1)
      q5 = q4 + x%d2val1
      q6 = q0*q5
      q7 = q1*x%d1val2
      q8 = pow2(x%d1val2)
      q9 = q8 + x%d2val2
      q10 = q2*x%d1val3
      q11 = pow2(x%d1val3)
      q12 = q11 + x%d2val3
      q13 = pow3(x%d1val1)
      q14 = 3.0_dp*x%d1val1
      q15 = q13 + q14*x%d2val1 + x%d3val1
      q16 = 2.0_dp*x%d1val1
      q17 = q16*x%d1val1_d1val2 + x%d2val1_d1val2
      q18 = q16*x%d1val1_d1val3 + x%d2val1_d1val3
      q19 = 2.0_dp*x%d1val2
      q20 = 16.0_dp*q8 + 16.0_dp*x%d2val2
      q21 = 0.0625_dp*x%d1val1
      q22 = q19*x%d1val1_d1val2 + q20*q21 + x%d1val1_d2val2
      q23 = 2.0_dp*x%d1val3
      q24 = pow3(x%d1val2)
      q25 = 3.0_dp*x%d1val2
      q26 = q24 + q25*x%d2val2 + x%d3val2
      q27 = q19*x%d1val2_d1val3 + x%d2val2_d1val3
      q28 = 3.0_dp*x%d1val1_d1val3
      q29 = 2.0_dp*x%d1val1_d1val2
      q30 = q16*x%d1val1_d1val2_d1val3 + q29*x%d1val1_d1val3 + x%d2val1_d1val2_d1val3
      q31 = 4.0_dp*q4 + 4.0_dp*x%d2val1
      q32 = 0.25_dp*q12
      q33 = pow2(x%d1val1_d1val3)
      q34 = 2.0_dp*q33 + q16*x%d1val1_d2val3 + x%d2val1_d2val3
      q35 = 2.0_dp*x%d1val2_d1val3
      q36 = x%d1val1*x%d1val2
      q37 = q19*x%d1val3
      q38 = 3.0_dp*x%d1val2_d1val3
      q39 = 4.0_dp*q8 + 4.0_dp*x%d2val2
      q40 = pow2(x%d1val2_d1val3)
      q41 = 2.0_dp*q40 + q19*x%d1val2_d2val3 + x%d2val2_d2val3
      q42 = 3.0_dp*x%d1val1_d2val3
      q43 = 0.125_dp*q12
      q44 = 12.0_dp*x%d1val1_d1val3
      q45 = 0.5_dp*x%d1val3
      q46 = 4.0_dp*x%d1val1_d1val2_d1val3
      q47 = 4.0_dp*x%d1val1
      q48 = 0.5_dp*q47*x%d1val1_d1val2 + x%d2val1_d1val2
      q49 = 0.25_dp*q31
      q50 = q49*x%d1val2
      q51 = 3.0_dp*x%d1val2_d2val3
      q52 = 12.0_dp*x%d1val2_d1val3
      unary%val = q0
      unary%d1val1 = q1
      unary%d1val2 = q2
      unary%d1val3 = q3
      unary%d2val1 = q6
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q7
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q1*x%d1val3
      unary%d2val2 = q0*q9
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q10
      unary%d2val3 = q0*q12
      unary%d3val1 = q0*q15
      unary%d2val1_d1val2 = q0*q17 + q2*q5
      unary%d2val1_d1val3 = q0*q18 + q3*q5
      unary%d1val1_d2val2 = q0*q22
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3 + q1*x%d1val2_d1val3 + q2*x%d1val1_d1val3 + q3*x%d1val1_d1val2 + q7*x%d1val3
      unary%d1val1_d2val3 = q0*(q12*x%d1val1 + q23*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q0*q26
      unary%d2val2_d1val3 = q0*q27 + q3*q9
      unary%d1val2_d2val3 = q0*(q12*x%d1val2 + q23*x%d1val2_d1val3 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q0*(q14*x%d2val1_d1val3 + q28*q4 + q28*x%d2val1 + x%d3val1_d1val3) + q15*q3
      unary%d2val1_d1val2_d1val3 = q0*q30 + q10*q5 + q17*q3 + q18*q2 + q6*x%d1val2_d1val3
      unary%d2val1_d2val3 = q0*(q18*q23 + q31*q32 + q34)
      unary%d1val1_d2val2_d1val3 = q0*(0.0625_dp*q20*x%d1val1_d1val3 + q19*x%d1val1_d1val2_d1val3 + q21*(16.0_dp*x%d2val2_d1val3 + 32.0_dp*x%d1val2*x%d1val2_d1val3) + q29*x%d1val2_d1val3 + x%d1val1_d2val2_d1val3) + q22*q3
      unary%d1val1_d1val2_d2val3 = q0*(q11*q36 + q11*x%d1val1_d1val2 + q16*x%d1val2_d1val3*x%d1val3 + q23*x%d1val1_d1val2_d1val3 + q35*x%d1val1_d1val3 + q36*x%d2val3 + q37*x%d1val1_d1val3 + x%d1val1*x%d1val2_d2val3 + x%d1val1_d1val2*x%d2val3 + x%d1val1_d1val2_d2val3 + x%d1val1_d2val3*x%d1val2)
      unary%d3val2_d1val3 = q0*(q25*x%d2val2_d1val3 + q38*q8 + q38*x%d2val2 + x%d3val2_d1val3) + q26*q3
      unary%d2val2_d2val3 = q0*(q23*q27 + q32*q39 + q41)
      unary%d3val1_d2val3 = q0*(6.0_dp*q33*x%d1val1 + 6.0_dp*x%d1val1_d1val3*x%d2val1_d1val3 + q14*x%d2val1_d2val3 + q4*q42 + q42*x%d2val1 + q43*(24.0_dp*x%d1val1*x%d2val1 + 8.0_dp*q13 + 8.0_dp*x%d3val1) + q45*(12.0_dp*x%d1val1*x%d2val1_d1val3 + 4.0_dp*x%d3val1_d1val3 + q4*q44 + q44*x%d2val1) + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q0*(q11*q48 + q11*q50 + q16*x%d1val1_d1val2_d2val3 + q18*q35 + q18*q37 + q23*q30 + q29*x%d1val1_d2val3 + q31*q45*x%d1val2_d1val3 + q34*x%d1val2 + q46*x%d1val1_d1val3 + q48*x%d2val3 + q49*x%d1val2_d2val3 + q50*x%d2val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q0*(0.25_dp*q39*x%d1val1_d2val3 + 2.0_dp*q27*x%d1val1_d1val3 + q19*x%d1val1_d1val2_d2val3 + q29*x%d1val2_d2val3 + q41*x%d1val1 + q43*(16.0_dp*x%d1val1_d1val2*x%d1val2 + 8.0_dp*x%d1val1_d2val2 + q16*q39) + q45*(4.0_dp*x%d1val1_d2val2_d1val3 + 8.0_dp*x%d1val1_d1val2*x%d1val2_d1val3 + 8.0_dp*x%d1val1_d1val2_d1val3*x%d1val2 + q27*q47 + q39*x%d1val1_d1val3) + q46*x%d1val2_d1val3 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q0*(6.0_dp*q40*x%d1val2 + 6.0_dp*x%d1val2_d1val3*x%d2val2_d1val3 + q25*x%d2val2_d2val3 + q43*(24.0_dp*x%d1val2*x%d2val2 + 8.0_dp*q24 + 8.0_dp*x%d3val2) + q45*(12.0_dp*x%d1val2*x%d2val2_d1val3 + 4.0_dp*x%d3val2_d1val3 + q52*q8 + q52*x%d2val2) + q51*q8 + q51*x%d2val2 + x%d3val2_d2val3)
   end function exp_self
   
   function expm1_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = exp(x%val)
      q1 = q0*x%d1val1
      q2 = q0*x%d1val2
      q3 = q0*x%d1val3
      q4 = pow2(x%d1val1)
      q5 = q4 + x%d2val1
      q6 = q0*q5
      q7 = q1*x%d1val2
      q8 = pow2(x%d1val2)
      q9 = q8 + x%d2val2
      q10 = q2*x%d1val3
      q11 = pow2(x%d1val3)
      q12 = q11 + x%d2val3
      q13 = pow3(x%d1val1)
      q14 = 3.0_dp*x%d1val1
      q15 = q13 + q14*x%d2val1 + x%d3val1
      q16 = 2.0_dp*x%d1val1
      q17 = q16*x%d1val1_d1val2 + x%d2val1_d1val2
      q18 = q16*x%d1val1_d1val3 + x%d2val1_d1val3
      q19 = 2.0_dp*x%d1val2
      q20 = 16.0_dp*q8 + 16.0_dp*x%d2val2
      q21 = 0.0625_dp*x%d1val1
      q22 = q19*x%d1val1_d1val2 + q20*q21 + x%d1val1_d2val2
      q23 = 2.0_dp*x%d1val3
      q24 = pow3(x%d1val2)
      q25 = 3.0_dp*x%d1val2
      q26 = q24 + q25*x%d2val2 + x%d3val2
      q27 = q19*x%d1val2_d1val3 + x%d2val2_d1val3
      q28 = 3.0_dp*x%d1val1_d1val3
      q29 = 2.0_dp*x%d1val1_d1val2
      q30 = q16*x%d1val1_d1val2_d1val3 + q29*x%d1val1_d1val3 + x%d2val1_d1val2_d1val3
      q31 = 4.0_dp*q4 + 4.0_dp*x%d2val1
      q32 = 0.25_dp*q12
      q33 = pow2(x%d1val1_d1val3)
      q34 = 2.0_dp*q33 + q16*x%d1val1_d2val3 + x%d2val1_d2val3
      q35 = 2.0_dp*x%d1val2_d1val3
      q36 = x%d1val1*x%d1val2
      q37 = q19*x%d1val3
      q38 = 3.0_dp*x%d1val2_d1val3
      q39 = 4.0_dp*q8 + 4.0_dp*x%d2val2
      q40 = pow2(x%d1val2_d1val3)
      q41 = 2.0_dp*q40 + q19*x%d1val2_d2val3 + x%d2val2_d2val3
      q42 = 3.0_dp*x%d1val1_d2val3
      q43 = 0.125_dp*q12
      q44 = 12.0_dp*x%d1val1_d1val3
      q45 = 0.5_dp*x%d1val3
      q46 = 4.0_dp*x%d1val1_d1val2_d1val3
      q47 = 4.0_dp*x%d1val1
      q48 = 0.5_dp*q47*x%d1val1_d1val2 + x%d2val1_d1val2
      q49 = 0.25_dp*q31
      q50 = q49*x%d1val2
      q51 = 3.0_dp*x%d1val2_d2val3
      q52 = 12.0_dp*x%d1val2_d1val3
      unary%val = expm1(x%val)
      unary%d1val1 = q1
      unary%d1val2 = q2
      unary%d1val3 = q3
      unary%d2val1 = q6
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q7
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q1*x%d1val3
      unary%d2val2 = q0*q9
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q10
      unary%d2val3 = q0*q12
      unary%d3val1 = q0*q15
      unary%d2val1_d1val2 = q0*q17 + q2*q5
      unary%d2val1_d1val3 = q0*q18 + q3*q5
      unary%d1val1_d2val2 = q0*q22
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3 + q1*x%d1val2_d1val3 + q2*x%d1val1_d1val3 + q3*x%d1val1_d1val2 + q7*x%d1val3
      unary%d1val1_d2val3 = q0*(q12*x%d1val1 + q23*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q0*q26
      unary%d2val2_d1val3 = q0*q27 + q3*q9
      unary%d1val2_d2val3 = q0*(q12*x%d1val2 + q23*x%d1val2_d1val3 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q0*(q14*x%d2val1_d1val3 + q28*q4 + q28*x%d2val1 + x%d3val1_d1val3) + q15*q3
      unary%d2val1_d1val2_d1val3 = q0*q30 + q10*q5 + q17*q3 + q18*q2 + q6*x%d1val2_d1val3
      unary%d2val1_d2val3 = q0*(q18*q23 + q31*q32 + q34)
      unary%d1val1_d2val2_d1val3 = q0*(0.0625_dp*q20*x%d1val1_d1val3 + q19*x%d1val1_d1val2_d1val3 + q21*(16.0_dp*x%d2val2_d1val3 + 32.0_dp*x%d1val2*x%d1val2_d1val3) + q29*x%d1val2_d1val3 + x%d1val1_d2val2_d1val3) + q22*q3
      unary%d1val1_d1val2_d2val3 = q0*(q11*q36 + q11*x%d1val1_d1val2 + q16*x%d1val2_d1val3*x%d1val3 + q23*x%d1val1_d1val2_d1val3 + q35*x%d1val1_d1val3 + q36*x%d2val3 + q37*x%d1val1_d1val3 + x%d1val1*x%d1val2_d2val3 + x%d1val1_d1val2*x%d2val3 + x%d1val1_d1val2_d2val3 + x%d1val1_d2val3*x%d1val2)
      unary%d3val2_d1val3 = q0*(q25*x%d2val2_d1val3 + q38*q8 + q38*x%d2val2 + x%d3val2_d1val3) + q26*q3
      unary%d2val2_d2val3 = q0*(q23*q27 + q32*q39 + q41)
      unary%d3val1_d2val3 = q0*(6.0_dp*q33*x%d1val1 + 6.0_dp*x%d1val1_d1val3*x%d2val1_d1val3 + q14*x%d2val1_d2val3 + q4*q42 + q42*x%d2val1 + q43*(24.0_dp*x%d1val1*x%d2val1 + 8.0_dp*q13 + 8.0_dp*x%d3val1) + q45*(12.0_dp*x%d1val1*x%d2val1_d1val3 + 4.0_dp*x%d3val1_d1val3 + q4*q44 + q44*x%d2val1) + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q0*(q11*q48 + q11*q50 + q16*x%d1val1_d1val2_d2val3 + q18*q35 + q18*q37 + q23*q30 + q29*x%d1val1_d2val3 + q31*q45*x%d1val2_d1val3 + q34*x%d1val2 + q46*x%d1val1_d1val3 + q48*x%d2val3 + q49*x%d1val2_d2val3 + q50*x%d2val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q0*(0.25_dp*q39*x%d1val1_d2val3 + 2.0_dp*q27*x%d1val1_d1val3 + q19*x%d1val1_d1val2_d2val3 + q29*x%d1val2_d2val3 + q41*x%d1val1 + q43*(16.0_dp*x%d1val1_d1val2*x%d1val2 + 8.0_dp*x%d1val1_d2val2 + q16*q39) + q45*(4.0_dp*x%d1val1_d2val2_d1val3 + 8.0_dp*x%d1val1_d1val2*x%d1val2_d1val3 + 8.0_dp*x%d1val1_d1val2_d1val3*x%d1val2 + q27*q47 + q39*x%d1val1_d1val3) + q46*x%d1val2_d1val3 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q0*(6.0_dp*q40*x%d1val2 + 6.0_dp*x%d1val2_d1val3*x%d2val2_d1val3 + q25*x%d2val2_d2val3 + q43*(24.0_dp*x%d1val2*x%d2val2 + 8.0_dp*q24 + 8.0_dp*x%d3val2) + q45*(12.0_dp*x%d1val2*x%d2val2_d1val3 + 4.0_dp*x%d3val2_d1val3 + q52*q8 + q52*x%d2val2) + q51*q8 + q51*x%d2val2 + x%d3val2_d2val3)
   end function expm1_self
   
   function exp10_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = pow(10.0_dp, x%val)
      q1 = ln10
      q2 = q0*q1
      q3 = pow2(x%d1val1)
      q4 = q1*q3
      q5 = q4 + x%d2val1
      q6 = pow2(q1)
      q7 = q0*q6
      q8 = x%d1val1*x%d1val2
      q9 = q7*x%d1val3
      q10 = pow2(x%d1val2)
      q11 = q1*q10
      q12 = q11 + x%d2val2
      q13 = pow2(x%d1val3)
      q14 = q1*q13 + x%d2val3
      q15 = q1*x%d1val1
      q16 = 3.0_dp*q15
      q17 = q6*pow3(x%d1val1)
      q18 = q16*x%d2val1 + q17 + x%d3val1
      q19 = 2.0_dp*q1
      q20 = q19*x%d1val1
      q21 = q20*x%d1val1_d1val2 + x%d2val1_d1val2
      q22 = q7*x%d1val2
      q23 = q20*x%d1val1_d1val3 + x%d2val1_d1val3
      q24 = q19*x%d1val2
      q25 = 16.0_dp*q11 + 16.0_dp*x%d2val2
      q26 = 0.0625_dp*q15
      q27 = q24*x%d1val1_d1val2 + q25*q26 + x%d1val1_d2val2
      q28 = q7*x%d1val2_d1val3
      q29 = q6*x%d1val1_d1val2
      q30 = q0*x%d1val3
      q31 = pow3(q1)
      q32 = q31*q8
      q33 = q19*x%d1val3
      q34 = q1*q14
      q35 = q1*x%d1val2
      q36 = 3.0_dp*q35
      q37 = q6*pow3(x%d1val2)
      q38 = q36*x%d2val2 + q37 + x%d3val2
      q39 = q24*x%d1val2_d1val3 + x%d2val2_d1val3
      q40 = q1*x%d1val1_d1val3
      q41 = q40*x%d2val1
      q42 = 3.0_dp*q6
      q43 = q3*x%d1val1_d1val3
      q44 = q31*x%d1val2
      q45 = q19*x%d1val1_d1val2
      q46 = q20*x%d1val1_d1val2_d1val3 + q45*x%d1val1_d1val3 + x%d2val1_d1val2_d1val3
      q47 = 4.0_dp*q4 + 4.0_dp*x%d2val1
      q48 = 0.25_dp*q34
      q49 = pow2(x%d1val1_d1val3)
      q50 = q19*q49 + q20*x%d1val1_d2val3 + x%d2val1_d2val3
      q51 = q1*x%d2val3
      q52 = q19*x%d1val2_d1val3
      q53 = q6*x%d2val3
      q54 = q6*x%d1val2_d1val3
      q55 = 2.0_dp*x%d1val3
      q56 = q55*q6*x%d1val2
      q57 = 3.0_dp*q1
      q58 = x%d1val2_d1val3*x%d2val2
      q59 = q10*x%d1val2_d1val3
      q60 = 4.0_dp*q11 + 4.0_dp*x%d2val2
      q61 = pow2(x%d1val2_d1val3)
      q62 = q19*q61 + q24*x%d1val2_d2val3 + x%d2val2_d2val3
      q63 = 6.0_dp*q6
      q64 = 12.0_dp*q6
      q65 = 0.5_dp*x%d1val3
      q66 = q1*q65
      q67 = 0.125_dp*q34
      q68 = 4.0_dp*x%d1val1_d1val2_d1val3
      q69 = 4.0_dp*q15
      q70 = 0.5_dp*q69*x%d1val1_d1val2 + x%d2val1_d1val2
      q71 = 0.25_dp*q1
      q72 = 0.25_dp*q47
      q73 = q1*x%d1val2_d1val3
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q2*q5
      unary%d1val1_d1val2 = q2*x%d1val1_d1val2 + q7*q8
      unary%d1val1_d1val3 = q2*x%d1val1_d1val3 + q9*x%d1val1
      unary%d2val2 = q12*q2
      unary%d1val2_d1val3 = q2*x%d1val2_d1val3 + q9*x%d1val2
      unary%d2val3 = q14*q2
      unary%d3val1 = q18*q2
      unary%d2val1_d1val2 = q2*q21 + q22*q5
      unary%d2val1_d1val3 = q2*q23 + q5*q9
      unary%d1val1_d2val2 = q2*q27
      unary%d1val1_d1val2_d1val3 = q2*x%d1val1_d1val2_d1val3 + q22*x%d1val1_d1val3 + q28*x%d1val1 + q29*q30 + q30*q32
      unary%d1val1_d2val3 = q2*(q33*x%d1val1_d1val3 + q34*x%d1val1 + x%d1val1_d2val3)
      unary%d3val2 = q2*q38
      unary%d2val2_d1val3 = q12*q9 + q2*q39
      unary%d1val2_d2val3 = q2*(q33*x%d1val2_d1val3 + q34*x%d1val2 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q18*q9 + q2*(3.0_dp*q41 + q16*x%d2val1_d1val3 + q42*q43 + x%d3val1_d1val3)
      unary%d2val1_d1val2_d1val3 = q2*q46 + q21*q9 + q22*q23 + q28*q5 + q30*q44*q5
      unary%d2val1_d2val3 = q2*(q23*q33 + q47*q48 + q50)
      unary%d1val1_d2val2_d1val3 = q2*(0.0625_dp*q25*q40 + q24*x%d1val1_d1val2_d1val3 + q26*(16.0_dp*x%d2val2_d1val3 + 32.0_dp*q35*x%d1val2_d1val3) + q45*x%d1val2_d1val3 + x%d1val1_d2val2_d1val3) + q27*q9
      unary%d1val1_d1val2_d2val3 = q2*(q13*q29 + q13*q32 + q15*x%d1val2_d2val3 + q33*x%d1val1_d1val2_d1val3 + q35*x%d1val1_d2val3 + q51*x%d1val1_d1val2 + q52*x%d1val1_d1val3 + q53*q8 + q54*q55*x%d1val1 + q56*x%d1val1_d1val3 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q2*(q36*x%d2val2_d1val3 + q42*q59 + q57*q58 + x%d3val2_d1val3) + q38*q9
      unary%d2val2_d2val3 = q2*(q33*q39 + q48*q60 + q62)
      unary%d3val1_d2val3 = q2*(6.0_dp*q40*x%d2val1_d1val3 + q16*x%d2val1_d2val3 + q3*q42*x%d1val1_d2val3 + q49*q63*x%d1val1 + q57*x%d1val1_d2val3*x%d2val1 + q66*(12.0_dp*q15*x%d2val1_d1val3 + 12.0_dp*q41 + 4.0_dp*x%d3val1_d1val3 + q43*q64) + q67*(24.0_dp*q15*x%d2val1 + 8.0_dp*q17 + 8.0_dp*x%d3val1) + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(q13*q44*q72 + q13*q6*q70 + q20*x%d1val1_d1val2_d2val3 + q23*q52 + q23*q56 + q33*q46 + q35*q50 + q40*q68 + q45*x%d1val1_d2val3 + q47*q54*q65 + q47*q71*x%d1val2_d2val3 + q51*q70 + q53*q72*x%d1val2 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(q15*q62 + q19*q39*x%d1val1_d1val3 + q24*x%d1val1_d1val2_d2val3 + q45*x%d1val2_d2val3 + q60*q71*x%d1val1_d2val3 + q66*(4.0_dp*x%d1val1_d2val2_d1val3 + 8.0_dp*q35*x%d1val1_d1val2_d1val3 + 8.0_dp*q73*x%d1val1_d1val2 + q39*q69 + q40*q60) + q67*(16.0_dp*q35*x%d1val1_d1val2 + 8.0_dp*x%d1val1_d2val2 + q20*q60) + q68*q73 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(6.0_dp*q73*x%d2val2_d1val3 + q10*q42*x%d1val2_d2val3 + q36*x%d2val2_d2val3 + q57*x%d1val2_d2val3*x%d2val2 + q61*q63*x%d1val2 + q66*(12.0_dp*q1*q58 + 12.0_dp*q35*x%d2val2_d1val3 + 4.0_dp*x%d3val2_d1val3 + q59*q64) + q67*(24.0_dp*q35*x%d2val2 + 8.0_dp*q37 + 8.0_dp*x%d3val2) + x%d3val2_d2val3)
   end function exp10_self
   
   function powm1_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = powm1(x%val)
      q1 = powm1(pow2(x%val))
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4 = q1*x%d1val3
      q5 = 2.0_dp*x%d2val1
      q6 = pow2(x%d1val1)
      q7 = 4.0_dp*q0
      q8 = q6*q7
      q9 = -q5 + q8
      q10 = 0.5_dp*q1
      q11 = q1*x%d1val1_d1val2
      q12 = powm1(pow3(x%val))
      q13 = q12*x%d1val2
      q14 = 2.0_dp*x%d1val1
      q15 = q1*x%d1val1_d1val3
      q16 = q12*x%d1val3
      q17 = 2.0_dp*x%d2val2
      q18 = pow2(x%d1val2)
      q19 = q18*q7
      q20 = -q17 + q19
      q21 = q1*x%d1val2_d1val3
      q22 = 2.0_dp*q16
      q23 = pow2(x%d1val3)
      q24 = q0*q23
      q25 = 2.0_dp*x%d3val1
      q26 = q0*x%d1val1
      q27 = 12.0_dp*q26
      q28 = q27*x%d2val1
      q29 = pow3(x%d1val1)
      q30 = 12.0_dp*q1
      q31 = q29*q30
      q32 = -q25 + q28 - q31
      q33 = q12*q9
      q34 = 2.0_dp*x%d2val1_d1val2
      q35 = 8.0_dp*q26
      q36 = q35*x%d1val1_d1val2
      q37 = 4.0_dp*q6
      q38 = q3*q37
      q39 = -q34 + q36 - q38
      q40 = -2.0_dp*x%d2val1_d1val3 + q35*x%d1val1_d1val3 - q37*q4
      q41 = 2.0_dp*x%d1val1_d2val2
      q42 = q0*x%d1val2
      q43 = 8.0_dp*x%d1val1_d1val2
      q44 = q42*q43
      q45 = q0*q18
      q46 = -12.0_dp*q45 + 4.0_dp*x%d2val2
      q47 = q26*q46 - q41 + q44
      q48 = q12*q14
      q49 = 2.0_dp*x%d1val1_d1val3
      q50 = 6.0_dp*x%d1val2
      q51 = powm1(pow4(x%val))
      q52 = q51*x%d1val3
      q53 = q7*x%d1val3
      q54 = -3.0_dp*q24 + x%d2val3
      q55 = 2.0_dp*q26
      q56 = 2.0_dp*x%d3val2
      q57 = 12.0_dp*q42
      q58 = q57*x%d2val2
      q59 = pow3(x%d1val2)
      q60 = q30*q59
      q61 = -q56 + q58 - q60
      q62 = q42*x%d1val2_d1val3
      q63 = q18*q4
      q64 = q0*q54
      q65 = 2.0_dp*x%d1val2
      q66 = q0*x%d1val1_d1val3
      q67 = q66*x%d2val1
      q68 = q2*x%d1val3
      q69 = 12.0_dp*q68
      q70 = q15*q6
      q71 = 24.0_dp*q16
      q72 = 8.0_dp*x%d1val2
      q73 = q2*x%d1val1_d1val3
      q74 = q16*q6
      q75 = q5 - q8
      q76 = q7*x%d1val1
      q77 = 2.0_dp*q6
      q78 = q4*q77 - q76*x%d1val1_d1val3 + x%d2val1_d1val3
      q79 = pow2(x%d1val1_d1val3)
      q80 = q1*x%d2val3
      q81 = q12*q23
      q82 = -8.0_dp*q68*x%d1val1_d1val3 + q37*q81 + q7*q79 + q76*x%d1val1_d2val3 - q77*q80 - x%d2val1_d2val3
      q83 = q0*x%d1val2_d1val3
      q84 = 8.0_dp*x%d1val1_d1val2_d1val3
      q85 = q3*x%d1val3
      q86 = q0*x%d2val3
      q87 = q7*x%d1val2_d1val3
      q88 = q0*x%d1val1_d2val3
      q89 = q2*x%d2val3
      q90 = 12.0_dp*q85
      q91 = q13*q23
      q92 = 12.0_dp*q83
      q93 = q18*q21
      q94 = -x%d2val2_d2val3
      q95 = q7*x%d1val2_d2val3
      q96 = pow2(x%d1val2_d1val3)
      q97 = 8.0_dp*x%d1val2_d1val3
      q98 = q18*q80
      q99 = q18*q81
      q100 = 6.0_dp*q26
      q101 = 12.0_dp*x%d2val1_d1val3
      q102 = 6.0_dp*x%d2val1
      q103 = 12.0_dp*x%d2val1
      q104 = q4*x%d1val1_d1val3
      q105 = 18.0_dp*q1
      q106 = 12.0_dp*q29
      q107 = q12*x%d2val3
      q108 = q81*x%d1val1
      q109 = q23*q51
      q110 = 36.0_dp*q109
      q111 = q7*x%d1val1_d1val2
      q112 = 4.0_dp*x%d1val1_d1val2
      q113 = 4.0_dp*x%d1val2
      q114 = 4.0_dp*q3
      q115 = q16*x%d1val2
      q116 = q0*x%d1val2_d2val3
      q117 = q3*x%d2val3
      q118 = q4*x%d1val2_d1val3
      q119 = 12.0_dp*q91
      q120 = q34 - q36 + q38
      q121 = q7*x%d1val2
      q122 = q114*x%d1val1_d1val2
      q123 = -6.0_dp*q45 + q17
      q124 = 6.0_dp*q0
      q125 = q124*x%d1val2
      q126 = 3.0_dp*q63 - q125*x%d1val2_d1val3 + x%d2val2_d1val3
      q127 = 6.0_dp*q116
      q128 = 6.0_dp*x%d2val2
      q129 = x%d1val2_d1val3*x%d2val2
      q130 = 12.0_dp*q59
      unary%val = q0
      unary%d1val1 = -q2
      unary%d1val2 = -q3
      unary%d1val3 = -q4
      unary%d2val1 = q10*q9
      unary%d1val1_d1val2 = -q11 + q13*q14
      unary%d1val1_d1val3 = q14*q16 - q15
      unary%d2val2 = q10*q20
      unary%d1val2_d1val3 = -q21 + q22*x%d1val2
      unary%d2val3 = q1*(2.0_dp*q24 - x%d2val3)
      unary%d3val1 = q10*q32
      unary%d2val1_d1val2 = q10*q39 - q33*x%d1val2
      unary%d2val1_d1val3 = q10*q40 - q33*x%d1val3
      unary%d1val1_d2val2 = q10*q47
      unary%d1val1_d1val2_d1val3 = -q1*x%d1val1_d1val2_d1val3 + q13*q49 + q22*x%d1val1_d1val2 + q48*x%d1val2_d1val3 - q50*q52*x%d1val1
      unary%d1val1_d2val3 = q1*(q53*x%d1val1_d1val3 + q54*q55 - x%d1val1_d2val3)
      unary%d3val2 = q10*q61
      unary%d2val2_d1val3 = q10*(-2.0_dp*x%d2val2_d1val3 - 4.0_dp*q63 + 8.0_dp*q62) - q16*q20
      unary%d1val2_d2val3 = q1*(q53*x%d1val2_d1val3 + q64*q65 - x%d1val2_d2val3)
      unary%d3val1_d1val3 = q10*(-2.0_dp*x%d3val1_d1val3 - 36.0_dp*q70 + 12.0_dp*q67 + q27*x%d2val1_d1val3 + q29*q71 - q69*x%d2val1) - q16*q32
      unary%d2val1_d1val2_d1val3 = 3.0_dp*q52*q9*x%d1val2 + q10*(-2.0_dp*x%d2val1_d1val2_d1val3 - q21*q37 + q35*x%d1val1_d1val2_d1val3 + q43*q66 - q43*q68 - q72*q73 + q72*q74) - q13*q40 - q16*q39 - q33*x%d1val2_d1val3
      unary%d2val1_d2val3 = q1*(q53*q78 + q64*q75 + q82)
      unary%d1val1_d2val2_d1val3 = q10*(-2.0_dp*x%d1val1_d2val2_d1val3 + q26*(-24.0_dp*q62 + 12.0_dp*q63 + 4.0_dp*x%d2val2_d1val3) + q42*q84 + q43*q83 - q43*q85 + q46*q66 - q46*q68) - q16*q47
      unary%d1val1_d1val2_d2val3 = q1*(-6.0_dp*q11*q23 + 2.0_dp*q86*x%d1val1_d1val2 + 24.0_dp*q91*x%d1val1 - q50*q89 + q53*x%d1val1_d1val2_d1val3 + q55*x%d1val2_d2val3 + q65*q88 - q69*x%d1val2_d1val3 + q87*x%d1val1_d1val3 - q90*x%d1val1_d1val3 - x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q10*(-2.0_dp*x%d3val2_d1val3 - 36.0_dp*q93 + q57*x%d2val2_d1val3 + q59*q71 - q90*x%d2val2 + q92*x%d2val2) - q16*q61
      unary%d2val2_d2val3 = q1*(-2.0_dp*q98 + 4.0_dp*q99 + q53*(2.0_dp*q63 - q87*x%d1val2 + x%d2val2_d1val3) + q64*(q17 - q19) + q7*q96 - q85*q97 + q94 + q95*x%d1val2)
      unary%d3val1_d2val3 = q1*(-36.0_dp*q2*q79 + 72.0_dp*q74*x%d1val1_d1val3 + q100*x%d2val1_d2val3 + q101*q66 - q101*q68 + q102*q88 - q102*q89 - q103*q104 + q103*q108 - q105*q6*x%d1val1_d2val3 + q106*q107 - q110*q29 + q53*(-6.0_dp*q67 + 18.0_dp*q70 - q100*x%d2val1_d1val3 + q102*q68 - q106*q16 + x%d3val1_d1val3) + q64*(q25 - q28 + q31) - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q1*(-12.0_dp*q109*q6*x%d1val2 - 2.0_dp*q42*q82 - 3.0_dp*q1*q120*q23 - 3.0_dp*q117*q75 - 6.0_dp*q118*q75 + 16.0_dp*q115*x%d1val1*x%d1val1_d1val3 - q1*q77*x%d1val2_d2val3 - q104*q43 + q108*q43 + q111*x%d1val1_d2val3 - q112*q89 - q113*q2*x%d1val1_d2val3 - q114*q79 + q116*q75 + q119*q75 + q120*q86 + q13*q37*x%d2val3 + q53*(-q111*x%d1val1_d1val3 + q112*q68 + q113*q73 - q115*q37 + q21*q77 - q76*x%d1val1_d1val2_d1val3 + x%d2val1_d1val2_d1val3) + q66*q84 - q68*q84 - q73*q97 + q74*q97 + q76*x%d1val1_d1val2_d2val3 + q78*q87 - q78*q90 - x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q1*(-4.0_dp*q126*q68 - q118*q43 + q121*x%d1val1_d1val2_d2val3 - q122*x%d2val3 + q123*q23*q48 - q123*q4*q49 + q123*q88 - q123*q89 + q126*q7*x%d1val1_d1val3 + q43*q91 + q53*(-q121*x%d1val1_d1val2_d1val3 + q122*x%d1val3 - q123*q66 + q123*q68 - q126*q55 - q87*x%d1val1_d1val2 + x%d1val1_d2val2_d1val3) - q55*(-3.0_dp*q98 + 6.0_dp*q99 + q124*q96 + q127*x%d1val2 - q90*x%d1val2_d1val3 + q94) + q64*(-q123*q55 + q41 - q44) + q83*q84 - q84*q85 + q95*x%d1val1_d1val2 - x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q1*(-12.0_dp*q129*q4 - 36.0_dp*q3*q96 + 72.0_dp*q16*q18*x%d1val2_d1val3 - q105*q18*x%d1val2_d2val3 + q107*q130 - q110*q59 - q117*q128 + q119*x%d2val2 + q125*x%d2val2_d2val3 + q127*x%d2val2 + q53*(18.0_dp*q93 - q124*q129 - q125*x%d2val2_d1val3 + q128*q85 - q130*q16 + x%d3val2_d1val3) + q64*(q56 - q58 + q60) - q90*x%d2val2_d1val3 + q92*x%d2val2_d1val3 - x%d3val2_d2val3)
   end function powm1_self
   
   function log_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = powm1(x%val)
      q1 = q0*x%d1val1
      q2 = q0*x%d1val2
      q3 = q0*x%d1val3
      q4 = pow2(x%d1val1)
      q5 = 4.0_dp*q0
      q6 = 4.0_dp*x%d2val1 - q4*q5
      q7 = 0.25_dp*q0
      q8 = q0*x%d1val1_d1val2
      q9 = powm1(pow2(x%val))
      q10 = q9*x%d1val1
      q11 = q10*x%d1val2
      q12 = q0*x%d1val1_d1val3
      q13 = q10*x%d1val3
      q14 = pow2(x%d1val2)
      q15 = -q14*q5
      q16 = 4.0_dp*x%d2val2 + q15
      q17 = q0*x%d1val2_d1val3
      q18 = q9*x%d1val3
      q19 = q18*x%d1val2
      q20 = pow2(x%d1val3)
      q21 = q0*q20
      q22 = 6.0_dp*q1
      q23 = pow3(x%d1val1)
      q24 = 4.0_dp*q9
      q25 = 2.0_dp*x%d3val1 - q22*x%d2val1 + q23*q24
      q26 = 0.5_dp*q0
      q27 = q25*q26
      q28 = q9*x%d1val2
      q29 = 0.25_dp*q6
      q30 = 8.0_dp*q1
      q31 = q24*q4
      q32 = 4.0_dp*x%d2val1_d1val2 - q30*x%d1val1_d1val2 + q31*x%d1val2
      q33 = q18*q4
      q34 = 4.0_dp*q33 + 4.0_dp*x%d2val1_d1val3 - q30*x%d1val1_d1val3
      q35 = 4.0_dp*q2
      q36 = 2.0_dp*x%d2val2
      q37 = q15 + q36
      q38 = 2.0_dp*x%d1val1_d2val2 - q1*q37 - q35*x%d1val1_d1val2
      q39 = q26*q38
      q40 = q10*x%d1val2_d1val3
      q41 = q18*x%d1val1_d1val2
      q42 = powm1(pow3(x%val))
      q43 = q42*x%d1val2
      q44 = q43*x%d1val3
      q45 = 2.0_dp*q3
      q46 = -2.0_dp*q21 + x%d2val3
      q47 = 6.0_dp*q2
      q48 = pow3(x%d1val2)
      q49 = 2.0_dp*x%d3val2 + q24*q48 - q47*x%d2val2
      q50 = q26*q49
      q51 = 0.25_dp*q18
      q52 = q2*x%d1val2_d1val3
      q53 = q14*q18
      q54 = -8.0_dp*q52 + 4.0_dp*q53
      q55 = 0.5_dp*q18
      q56 = 6.0_dp*q12
      q57 = 6.0_dp*x%d2val1
      q58 = q4*q9
      q59 = q58*x%d1val1_d1val3
      q60 = q23*q42
      q61 = 8.0_dp*x%d1val3
      q62 = q9*x%d1val2_d1val3
      q63 = q8*x%d1val1_d1val3
      q64 = q11*x%d1val1_d1val3
      q65 = q43*q61
      q66 = 2.0_dp*q0
      q67 = 2.0_dp*x%d2val1 - q4*q66
      q68 = q26*q46
      q69 = 2.0_dp*q1
      q70 = q33 - q69*x%d1val1_d1val3 + x%d2val1_d1val3
      q71 = pow2(x%d1val1_d1val3)
      q72 = 4.0_dp*q13
      q73 = q9*x%d2val3
      q74 = 2.0_dp*q4
      q75 = q20*q42
      q76 = -q4*q73 + q66*q71 + q69*x%d1val1_d2val3 - q72*x%d1val1_d1val3 + q74*q75 - x%d2val1_d2val3
      q77 = q8*x%d1val2_d1val3
      q78 = 4.0_dp*x%d1val1_d1val2
      q79 = q12*q37
      q80 = q13*q37
      q81 = 2.0_dp*q12
      q82 = 4.0_dp*q19
      q83 = q20*q9
      q84 = 2.0_dp*x%d1val1_d1val2
      q85 = q75*x%d1val1
      q86 = 6.0_dp*x%d1val2
      q87 = 6.0_dp*q17
      q88 = 6.0_dp*q19
      q89 = q14*q62
      q90 = q42*q48
      q91 = -x%d2val2_d2val3
      q92 = 2.0_dp*x%d1val2_d2val3
      q93 = pow2(x%d1val2_d1val3)
      q94 = q14*q73
      q95 = q14*q75
      q96 = 3.0_dp*q1
      q97 = 3.0_dp*x%d2val1
      q98 = q10*x%d2val3
      q99 = q18*x%d1val1_d1val3
      q100 = 4.0_dp*q60
      q101 = q4*q42*x%d1val3
      q102 = q20*powm1(pow4(x%val))
      q103 = 12.0_dp*q102
      q104 = 2.0_dp*x%d1val1_d2val3
      q105 = 4.0_dp*x%d1val1_d1val2_d1val3
      q106 = 4.0_dp*x%d1val1_d1val3
      q107 = 4.0_dp*x%d1val2_d1val3
      q108 = q28*x%d2val3
      q109 = q18*x%d1val2_d1val3
      q110 = q75*x%d1val2
      q111 = 2.0_dp*x%d2val1_d1val2 - q1*q78 + q28*q74
      q112 = 2.0_dp*q2
      q113 = 2.0_dp*q53 - q35*x%d1val2_d1val3 + x%d2val2_d1val3
      q114 = 3.0_dp*q2
      q115 = 3.0_dp*x%d2val2
      q116 = 4.0_dp*q90
      unary%val = log(x%val)
      unary%d1val1 = q1
      unary%d1val2 = q2
      unary%d1val3 = q3
      unary%d2val1 = q6*q7
      unary%d1val1_d1val2 = -q11 + q8
      unary%d1val1_d1val3 = q12 - q13
      unary%d2val2 = q16*q7
      unary%d1val2_d1val3 = q17 - q19
      unary%d2val3 = q0*(-q21 + x%d2val3)
      unary%d3val1 = q27
      unary%d2val1_d1val2 = -q28*q29 + q32*q7
      unary%d2val1_d1val3 = -q18*q29 + q34*q7
      unary%d1val1_d2val2 = q39
      unary%d1val1_d1val2_d1val3 = 2.0_dp*q44*x%d1val1 + q0*x%d1val1_d1val2_d1val3 - q28*x%d1val1_d1val3 - q40 - q41
      unary%d1val1_d2val3 = q0*(-q1*q46 - q45*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q50
      unary%d2val2_d1val3 = -q16*q51 + q7*(4.0_dp*x%d2val2_d1val3 + q54)
      unary%d1val2_d2val3 = q0*(-q2*q46 - q45*x%d1val2_d1val3 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q25*q55 + q26*(12.0_dp*q59 + 2.0_dp*x%d3val1_d1val3 + q13*q57 - q22*x%d2val1_d1val3 - q56*x%d2val1 - q60*q61)
      unary%d2val1_d1val2_d1val3 = -0.25_dp*q28*q34 + 0.5_dp*q44*q6 - q29*q62 - q32*q51 + q7*(-8.0_dp*q63 + 4.0_dp*x%d2val1_d1val2_d1val3 + 8.0_dp*q13*x%d1val1_d1val2 + 8.0_dp*q64 - q30*x%d1val1_d1val2_d1val3 + q31*x%d1val2_d1val3 - q4*q65)
      unary%d2val1_d2val3 = -q0*(q45*q70 + q67*q68 + q76)
      unary%d1val1_d2val2_d1val3 = q26*(-4.0_dp*q77 + 2.0_dp*x%d1val1_d2val2_d1val3 - q1*(2.0_dp*x%d2val2_d1val3 + q54) + q19*q78 - q35*x%d1val1_d1val2_d1val3 - q79 + q80) - q38*q55
      unary%d1val1_d1val2_d2val3 = q0*(2.0_dp*q11*x%d2val3 - q1*x%d1val2_d2val3 - q2*x%d1val1_d2val3 - q45*x%d1val1_d1val2_d1val3 + q72*x%d1val2_d1val3 - q8*x%d2val3 - q81*x%d1val2_d1val3 + q82*x%d1val1_d1val3 + q83*q84 - q85*q86 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q26*(12.0_dp*q89 + 2.0_dp*x%d3val2_d1val3 - q47*x%d2val2_d1val3 - q61*q90 - q87*x%d2val2 + q88*x%d2val2) - q49*q55
      unary%d2val2_d2val3 = -q0*(2.0_dp*q95 + q2*q92 + q45*(-2.0_dp*q52 + q53 + x%d2val2_d1val3) + q66*q93 + q68*(-q14*q66 + q36) - q82*x%d1val2_d1val3 + q91 - q94)
      unary%d3val1_d2val3 = q0*(-24.0_dp*q101*x%d1val1_d1val3 + 12.0_dp*q10*q71 + 6.0_dp*q13*x%d2val1_d1val3 + 6.0_dp*q58*x%d1val1_d2val3 - q0*q97*x%d1val1_d2val3 - q100*x%d2val3 + q103*q23 - q27*q46 - q45*(6.0_dp*q59 - q100*x%d1val3 - q12*q97 + q13*q97 - q96*x%d2val1_d1val3 + x%d3val1_d1val3) - q56*x%d2val1_d1val3 - q57*q85 + q57*q99 - q96*x%d2val1_d2val3 + q97*q98 + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q0*(-2.0_dp*q17*q70 - 3.0_dp*q110*q67 + 2.0_dp*q109*q67 + 2.0_dp*q28*q71 - q101*q107 + q102*q4*q86 + q104*q11 - q104*q8 - q105*q12 + q105*q13 + q106*q40 + q106*q41 + q108*q67 - q111*q26*x%d2val3 + q111*q83 + q2*q76 - q26*q67*x%d1val2_d2val3 - q43*q74*x%d2val3 - q45*(-2.0_dp*q63 + 2.0_dp*q64 + q13*q84 - q44*q74 + q58*x%d1val2_d1val3 - q69*x%d1val1_d1val2_d1val3 + x%d2val1_d1val2_d1val3) + q58*x%d1val2_d2val3 - q65*x%d1val1*x%d1val1_d1val3 - q69*x%d1val1_d1val2_d2val3 + q70*q82 - q78*q85 + q84*q98 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q0*(0.5_dp*q37*q98 + 2.0_dp*q113*q13 + q1*(-2.0_dp*q94 - 8.0_dp*q19*x%d1val2_d1val3 + 4.0_dp*q95 + q35*x%d1val2_d2val3 + q5*q93 + q91) - q105*q17 + q105*q19 + q107*q41 + q108*q84 - q110*q78 - q112*x%d1val1_d1val2_d2val3 - q113*q81 - q26*q37*x%d1val1_d2val3 - q37*q85 + q37*q99 - q39*q46 - q45*(-0.5_dp*q79 - 2.0_dp*q77 + 0.5_dp*q80 - q1*q113 - q112*x%d1val1_d1val2_d1val3 + q19*q84 + x%d1val1_d2val2_d1val3) - q8*q92 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q0*(-24.0_dp*q14*q42*x%d1val2_d1val3*x%d1val3 + 12.0_dp*q28*q93 + 6.0_dp*q109*x%d2val2 + 6.0_dp*q14*q9*x%d1val2_d2val3 - q0*q115*x%d1val2_d2val3 + q103*q48 + q108*q115 - q114*x%d2val2_d2val3 - q116*x%d2val3 - q45*(6.0_dp*q89 - q114*x%d2val2_d1val3 - q115*q17 + q115*q19 - q116*x%d1val3 + x%d3val2_d1val3) - q46*q50 - q75*q86*x%d2val2 - q87*x%d2val2_d1val3 + q88*x%d2val2_d1val3 + x%d3val2_d2val3)
   end function log_self
   
   function log1p_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = x%val + 1
      q1 = powm1(q0)
      q2 = pow2(x%d1val1)
      q3 = 12.0_dp*x%val + 12.0_dp
      q4 = powm1(q3)
      q5 = 48.0_dp*q4
      q6 = 4.0_dp*x%d2val1 - q2*q5
      q7 = 3.0_dp*q4
      q8 = powm1(pow2(q0))
      q9 = q8*x%d1val1
      q10 = pow2(x%d1val2)
      q11 = -q10*q5
      q12 = 4.0_dp*x%d2val2 + q11
      q13 = q8*x%d1val3
      q14 = pow2(x%d1val3)
      q15 = 2.0_dp*x%val + 2.0_dp
      q16 = powm1(q15)
      q17 = 2.0_dp*q16
      q18 = 2.0_dp*x%d3val1
      q19 = q4*x%d1val1
      q20 = 72.0_dp*q19
      q21 = pow3(x%d1val1)
      q22 = powm1(pow2(q3))
      q23 = 576.0_dp*q22
      q24 = q18 - q20*x%d2val1 + q21*q23
      q25 = 6.0_dp*q4
      q26 = 36.0_dp*q22
      q27 = q6*x%d1val2
      q28 = 96.0_dp*q19
      q29 = q2*q23
      q30 = 4.0_dp*x%d2val1_d1val2 - q28*x%d1val1_d1val2 + q29*x%d1val2
      q31 = q26*x%d1val3
      q32 = 4.0_dp*x%d2val1_d1val3 - q28*x%d1val1_d1val3 + q29*x%d1val3
      q33 = 2.0_dp*x%d1val1_d2val2
      q34 = q5*x%d1val1_d1val2
      q35 = 2.0_dp*x%d2val2
      q36 = q11 + q35
      q37 = 12.0_dp*q19
      q38 = q33 - q34*x%d1val2 - q36*q37
      q39 = x%d1val1_d1val3*x%d1val2
      q40 = x%d1val2*x%d1val3
      q41 = 2.0_dp*x%d1val1
      q42 = 4.0_dp*q16
      q43 = q42*x%d1val3
      q44 = -q14*q42 + x%d2val3
      q45 = q17*q44
      q46 = 2.0_dp*x%d3val2
      q47 = q4*x%d1val2
      q48 = 72.0_dp*q47
      q49 = pow3(x%d1val2)
      q50 = q23*q49 + q46 - q48*x%d2val2
      q51 = q10*x%d1val3
      q52 = -96.0_dp*q47*x%d1val2_d1val3 + q23*q51
      q53 = q22*x%d1val3
      q54 = 72.0_dp*q53
      q55 = q4*x%d1val1_d1val3
      q56 = x%d1val1*x%d2val1
      q57 = 864.0_dp*q53
      q58 = 1728.0_dp*q22
      q59 = q2*x%d1val1_d1val3
      q60 = q21*x%d1val3
      q61 = powm1(pow3(q3))
      q62 = 13824.0_dp*q61
      q63 = 1152.0_dp*x%d1val1
      q64 = q2*q40
      q65 = 2.0_dp*x%d2val1 - q2*q42
      q66 = q16*q44
      q67 = q42*x%d1val1
      q68 = powm1(pow2(q15))
      q69 = 4.0_dp*q68
      q70 = q2*q69
      q71 = -q67*x%d1val1_d1val3 + q70*x%d1val3 + x%d2val1_d1val3
      q72 = pow2(x%d1val1_d1val3)
      q73 = 16.0_dp*q68
      q74 = q73*x%d1val3
      q75 = q74*x%d1val1
      q76 = q2*x%d2val3
      q77 = powm1(pow3(q15))
      q78 = q14*q77
      q79 = 16.0_dp*q78
      q80 = q2*q79 + q42*q72 + q67*x%d1val1_d2val3 - q69*q76 - q75*x%d1val1_d1val3 - x%d2val1_d2val3
      q81 = x%d1val1_d1val2_d1val3*x%d1val2
      q82 = q17*x%d1val1
      q83 = x%d1val1_d1val2*x%d2val3
      q84 = q42*x%d1val2_d1val3
      q85 = q17*x%d1val2
      q86 = 8.0_dp*q68
      q87 = q86*x%d1val2
      q88 = q87*x%d1val1
      q89 = q86*x%d1val1_d1val2
      q90 = 48.0_dp*q78
      q91 = x%d1val2_d1val3*x%d2val2
      q92 = x%d1val2*x%d2val2
      q93 = q10*x%d1val2_d1val3
      q94 = q49*x%d1val3
      q95 = -x%d2val2_d2val3
      q96 = q42*x%d1val2_d2val3
      q97 = pow2(x%d1val2_d1val3)
      q98 = x%d1val2*x%d1val2_d1val3
      q99 = q10*x%d2val3
      q100 = 6.0_dp*q16
      q101 = q100*x%d1val1
      q102 = 12.0_dp*q16
      q103 = q100*x%d2val1
      q104 = 24.0_dp*q68
      q105 = q104*x%d1val3
      q106 = 12.0_dp*q56
      q107 = q68*x%d2val3
      q108 = 48.0_dp*q68
      q109 = 32.0_dp*q77
      q110 = q109*x%d2val3
      q111 = 192.0_dp*q77
      q112 = q14*powm1(pow4(q15))
      q113 = 192.0_dp*q112
      q114 = q68*x%d1val3
      q115 = q42*x%d1val1_d1val2
      q116 = 8.0_dp*q16
      q117 = q116*x%d1val1_d1val2_d1val3
      q118 = q86*x%d1val1
      q119 = q74*x%d1val1_d1val2
      q120 = x%d1val1*x%d1val3
      q121 = 32.0_dp*q78
      q122 = q121*x%d1val1_d1val2
      q123 = 16.0_dp*q77
      q124 = x%d1val2_d1val3*x%d1val3
      q125 = q16*x%d1val2_d2val3
      q126 = q65*x%d1val2
      q127 = q116*x%d1val1_d1val2
      q128 = 2.0_dp*x%d2val1_d1val2 - q127*x%d1val1 + q2*q87
      q129 = -q10*q116 + q35
      q130 = q129*q16
      q131 = q129*q41
      q132 = -q116*q98 + q51*q86 + x%d2val2_d1val3
      q133 = q100*x%d1val2
      q134 = 12.0_dp*q92
      unary%val = log1p(x%val)
      unary%d1val1 = q1*x%d1val1
      unary%d1val2 = q1*x%d1val2
      unary%d1val3 = q1*x%d1val3
      unary%d2val1 = q6*q7
      unary%d1val1_d1val2 = q1*x%d1val1_d1val2 - q9*x%d1val2
      unary%d1val1_d1val3 = q1*x%d1val1_d1val3 - q9*x%d1val3
      unary%d2val2 = q12*q7
      unary%d1val2_d1val3 = q1*x%d1val2_d1val3 - q13*x%d1val2
      unary%d2val3 = q17*(-q14*q17 + x%d2val3)
      unary%d3val1 = q24*q25
      unary%d2val1_d1val2 = -q26*q27 + q30*q7
      unary%d2val1_d1val3 = -q31*q6 + q32*q7
      unary%d1val1_d2val2 = q25*q38
      unary%d1val1_d1val2_d1val3 = q1*x%d1val1_d1val2_d1val3 - q13*x%d1val1_d1val2 - q39*q8 + q40*q41*powm1(pow3(q0)) - q9*x%d1val2_d1val3
      unary%d1val1_d2val3 = q17*(-q43*x%d1val1_d1val3 - q45*x%d1val1 + x%d1val1_d2val3)
      unary%d3val2 = q25*q50
      unary%d2val2_d1val3 = -q12*q31 + q7*(4.0_dp*x%d2val2_d1val3 + q52)
      unary%d1val2_d2val3 = q17*(-q43*x%d1val2_d1val3 - q45*x%d1val2 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q24*q54 + q25*(-72.0_dp*q55*x%d2val1 + 2.0_dp*x%d3val1_d1val3 - q20*x%d2val1_d1val3 + q56*q57 + q58*q59 - q60*q62)
      unary%d2val1_d1val2_d1val3 = 864.0_dp*q27*q61*x%d1val3 - q26*q32*x%d1val2 - q26*q6*x%d1val2_d1val3 - q30*q31 + q7*(-96.0_dp*q55*x%d1val1_d1val2 + 4.0_dp*x%d2val1_d1val2_d1val3 + q22*q39*q63 - q28*x%d1val1_d1val2_d1val3 + q29*x%d1val2_d1val3 + q53*q63*x%d1val1_d1val2 - q62*q64)
      unary%d2val1_d2val3 = -q17*(q43*q71 + q65*q66 + q80)
      unary%d1val1_d2val2_d1val3 = q25*(-12.0_dp*q36*q55 + 144.0_dp*q36*q53*x%d1val1 + 2.0_dp*x%d1val1_d2val2_d1val3 + q23*q40*x%d1val1_d1val2 - q34*x%d1val2_d1val3 - q37*(2.0_dp*x%d2val2_d1val3 + q52) - q5*q81) - q38*q54
      unary%d1val1_d1val2_d2val3 = q17*(q14*q89 - q17*q83 + q39*q74 - q43*x%d1val1_d1val2_d1val3 + q75*x%d1val2_d1val3 - q82*x%d1val2_d2val3 - q84*x%d1val1_d1val3 - q85*x%d1val1_d2val3 + q88*x%d2val3 - q90*x%d1val1*x%d1val2 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q25*(-72.0_dp*q4*q91 + 2.0_dp*x%d3val2_d1val3 - q48*x%d2val2_d1val3 + q57*q92 + q58*q93 - q62*q94) - q50*q54
      unary%d2val2_d2val3 = -q17*(q10*q79 + q42*q97 + q43*(q51*q69 - q84*x%d1val2 + x%d2val2_d1val3) + q66*(-q10*q42 + q35) - q69*q99 - q74*q98 + q95 + q96*x%d1val2)
      unary%d3val1_d2val3 = q17*(-q101*x%d2val1_d2val3 - q102*x%d1val1_d1val3*x%d2val1_d1val3 - q103*x%d1val1_d2val3 + q104*q2*x%d1val1_d2val3 + q105*x%d1val1*x%d2val1_d1val3 + q105*x%d1val1_d1val3*x%d2val1 + q106*q107 + q108*q72*x%d1val1 - q110*q21 - q111*q59*x%d1val3 + q113*q21 - q43*(-q101*x%d2val1_d1val3 - q103*x%d1val1_d1val3 + q104*q59 + q106*q114 - q109*q60 + x%d3val1_d1val3) - q56*q90 - q66*(-q102*q56 + q18 + q21*q73) + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q17*(-24.0_dp*q126*q78 - 64.0_dp*q120*q39*q77 + 96.0_dp*q112*q2*x%d1val2 - q109*q124*q2 - q115*x%d1val1_d2val3 - q117*x%d1val1_d1val3 + q118*q83 + q119*x%d1val1_d1val3 - q122*x%d1val1 - q123*q76*x%d1val2 + q124*q65*q86 - q125*q65 + q126*q69*x%d2val3 + q128*q14*q69 - q128*q16*x%d2val3 + q40*q71*q73 - q43*(-q115*x%d1val1_d1val3 + q118*q39 + q120*q89 - q123*q64 - q67*x%d1val1_d1val2_d1val3 + q70*x%d1val2_d1val3 + x%d2val1_d1val2_d1val3) - q67*x%d1val1_d1val2_d2val3 + q70*x%d1val2_d2val3 - q71*q84 + q72*q87 + q73*x%d1val1*x%d1val1_d1val3*x%d1val2_d1val3 + q75*x%d1val1_d1val2_d1val3 + q80*q85 + q88*x%d1val1_d2val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q17*(-8.0_dp*q129*q78*x%d1val1 + q107*q131 - q117*x%d1val2_d1val3 + q119*x%d1val2_d1val3 + q120*q132*q86 - q122*x%d1val2 + q129*q69*x%d1val1_d1val3*x%d1val3 - q130*x%d1val1_d2val3 - q132*q42*x%d1val1_d1val3 - q42*x%d1val1_d1val2_d2val3*x%d1val2 - q43*(q114*q131 - q130*x%d1val1_d1val3 - q132*q82 + q40*q89 - q42*q81 - q84*x%d1val1_d1val2 + x%d1val1_d2val2_d1val3) - q66*(-q127*x%d1val2 - q129*q82 + q33) + q74*q81 + q82*(-32.0_dp*q114*q98 + 8.0_dp*q125*x%d1val2 + q10*q121 + q116*q97 - q86*q99 + q95) + q83*q87 - q96*x%d1val1_d1val2 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q17*(-6.0_dp*q125*x%d2val2 + q10*q104*x%d1val2_d2val3 - q102*x%d1val2_d1val3*x%d2val2_d1val3 + q104*q40*x%d2val2_d1val3 + q105*q91 + q107*q134 + q108*q97*x%d1val2 - q110*q49 - q111*q51*x%d1val2_d1val3 + q113*q49 - q133*x%d2val2_d2val3 - q43*(-q100*q91 + q104*q93 - q109*q94 + q114*q134 - q133*x%d2val2_d1val3 + x%d3val2_d1val3) - q66*(-q102*q92 + q46 + q49*q73) - q90*q92 + x%d3val2_d2val3)
   end function log1p_self
   
   function safe_log_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = powm1(x%val)
      q1 = q0*x%d1val1
      q2 = q0*x%d1val2
      q3 = q0*x%d1val3
      q4 = pow2(x%d1val1)
      q5 = 4.0_dp*q0
      q6 = 4.0_dp*x%d2val1 - q4*q5
      q7 = 0.25_dp*q0
      q8 = q0*x%d1val1_d1val2
      q9 = powm1(pow2(x%val))
      q10 = q9*x%d1val1
      q11 = q10*x%d1val2
      q12 = q0*x%d1val1_d1val3
      q13 = q10*x%d1val3
      q14 = pow2(x%d1val2)
      q15 = -q14*q5
      q16 = 4.0_dp*x%d2val2 + q15
      q17 = q0*x%d1val2_d1val3
      q18 = q9*x%d1val3
      q19 = q18*x%d1val2
      q20 = pow2(x%d1val3)
      q21 = q0*q20
      q22 = 6.0_dp*q1
      q23 = pow3(x%d1val1)
      q24 = 4.0_dp*q9
      q25 = 2.0_dp*x%d3val1 - q22*x%d2val1 + q23*q24
      q26 = 0.5_dp*q0
      q27 = q25*q26
      q28 = q9*x%d1val2
      q29 = 0.25_dp*q6
      q30 = 8.0_dp*q1
      q31 = q24*q4
      q32 = 4.0_dp*x%d2val1_d1val2 - q30*x%d1val1_d1val2 + q31*x%d1val2
      q33 = q18*q4
      q34 = 4.0_dp*q33 + 4.0_dp*x%d2val1_d1val3 - q30*x%d1val1_d1val3
      q35 = 4.0_dp*q2
      q36 = 2.0_dp*x%d2val2
      q37 = q15 + q36
      q38 = 2.0_dp*x%d1val1_d2val2 - q1*q37 - q35*x%d1val1_d1val2
      q39 = q26*q38
      q40 = q10*x%d1val2_d1val3
      q41 = q18*x%d1val1_d1val2
      q42 = powm1(pow3(x%val))
      q43 = q42*x%d1val2
      q44 = q43*x%d1val3
      q45 = 2.0_dp*q3
      q46 = -2.0_dp*q21 + x%d2val3
      q47 = 6.0_dp*q2
      q48 = pow3(x%d1val2)
      q49 = 2.0_dp*x%d3val2 + q24*q48 - q47*x%d2val2
      q50 = q26*q49
      q51 = 0.25_dp*q18
      q52 = q2*x%d1val2_d1val3
      q53 = q14*q18
      q54 = -8.0_dp*q52 + 4.0_dp*q53
      q55 = 0.5_dp*q18
      q56 = 6.0_dp*q12
      q57 = 6.0_dp*x%d2val1
      q58 = q4*q9
      q59 = q58*x%d1val1_d1val3
      q60 = q23*q42
      q61 = 8.0_dp*x%d1val3
      q62 = q9*x%d1val2_d1val3
      q63 = q8*x%d1val1_d1val3
      q64 = q11*x%d1val1_d1val3
      q65 = q43*q61
      q66 = 2.0_dp*q0
      q67 = 2.0_dp*x%d2val1 - q4*q66
      q68 = q26*q46
      q69 = 2.0_dp*q1
      q70 = q33 - q69*x%d1val1_d1val3 + x%d2val1_d1val3
      q71 = pow2(x%d1val1_d1val3)
      q72 = 4.0_dp*q13
      q73 = q9*x%d2val3
      q74 = 2.0_dp*q4
      q75 = q20*q42
      q76 = -q4*q73 + q66*q71 + q69*x%d1val1_d2val3 - q72*x%d1val1_d1val3 + q74*q75 - x%d2val1_d2val3
      q77 = q8*x%d1val2_d1val3
      q78 = 4.0_dp*x%d1val1_d1val2
      q79 = q12*q37
      q80 = q13*q37
      q81 = 2.0_dp*q12
      q82 = 4.0_dp*q19
      q83 = q20*q9
      q84 = 2.0_dp*x%d1val1_d1val2
      q85 = q75*x%d1val1
      q86 = 6.0_dp*x%d1val2
      q87 = 6.0_dp*q17
      q88 = 6.0_dp*q19
      q89 = q14*q62
      q90 = q42*q48
      q91 = -x%d2val2_d2val3
      q92 = 2.0_dp*x%d1val2_d2val3
      q93 = pow2(x%d1val2_d1val3)
      q94 = q14*q73
      q95 = q14*q75
      q96 = 3.0_dp*q1
      q97 = 3.0_dp*x%d2val1
      q98 = q10*x%d2val3
      q99 = q18*x%d1val1_d1val3
      q100 = 4.0_dp*q60
      q101 = q4*q42*x%d1val3
      q102 = q20*powm1(pow4(x%val))
      q103 = 12.0_dp*q102
      q104 = 2.0_dp*x%d1val1_d2val3
      q105 = 4.0_dp*x%d1val1_d1val2_d1val3
      q106 = 4.0_dp*x%d1val1_d1val3
      q107 = 4.0_dp*x%d1val2_d1val3
      q108 = q28*x%d2val3
      q109 = q18*x%d1val2_d1val3
      q110 = q75*x%d1val2
      q111 = 2.0_dp*x%d2val1_d1val2 - q1*q78 + q28*q74
      q112 = 2.0_dp*q2
      q113 = 2.0_dp*q53 - q35*x%d1val2_d1val3 + x%d2val2_d1val3
      q114 = 3.0_dp*q2
      q115 = 3.0_dp*x%d2val2
      q116 = 4.0_dp*q90
      unary%val = safe_log(x%val)
      unary%d1val1 = q1
      unary%d1val2 = q2
      unary%d1val3 = q3
      unary%d2val1 = q6*q7
      unary%d1val1_d1val2 = -q11 + q8
      unary%d1val1_d1val3 = q12 - q13
      unary%d2val2 = q16*q7
      unary%d1val2_d1val3 = q17 - q19
      unary%d2val3 = q0*(-q21 + x%d2val3)
      unary%d3val1 = q27
      unary%d2val1_d1val2 = -q28*q29 + q32*q7
      unary%d2val1_d1val3 = -q18*q29 + q34*q7
      unary%d1val1_d2val2 = q39
      unary%d1val1_d1val2_d1val3 = 2.0_dp*q44*x%d1val1 + q0*x%d1val1_d1val2_d1val3 - q28*x%d1val1_d1val3 - q40 - q41
      unary%d1val1_d2val3 = q0*(-q1*q46 - q45*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q50
      unary%d2val2_d1val3 = -q16*q51 + q7*(4.0_dp*x%d2val2_d1val3 + q54)
      unary%d1val2_d2val3 = q0*(-q2*q46 - q45*x%d1val2_d1val3 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q25*q55 + q26*(12.0_dp*q59 + 2.0_dp*x%d3val1_d1val3 + q13*q57 - q22*x%d2val1_d1val3 - q56*x%d2val1 - q60*q61)
      unary%d2val1_d1val2_d1val3 = -0.25_dp*q28*q34 + 0.5_dp*q44*q6 - q29*q62 - q32*q51 + q7*(-8.0_dp*q63 + 4.0_dp*x%d2val1_d1val2_d1val3 + 8.0_dp*q13*x%d1val1_d1val2 + 8.0_dp*q64 - q30*x%d1val1_d1val2_d1val3 + q31*x%d1val2_d1val3 - q4*q65)
      unary%d2val1_d2val3 = -q0*(q45*q70 + q67*q68 + q76)
      unary%d1val1_d2val2_d1val3 = q26*(-4.0_dp*q77 + 2.0_dp*x%d1val1_d2val2_d1val3 - q1*(2.0_dp*x%d2val2_d1val3 + q54) + q19*q78 - q35*x%d1val1_d1val2_d1val3 - q79 + q80) - q38*q55
      unary%d1val1_d1val2_d2val3 = q0*(2.0_dp*q11*x%d2val3 - q1*x%d1val2_d2val3 - q2*x%d1val1_d2val3 - q45*x%d1val1_d1val2_d1val3 + q72*x%d1val2_d1val3 - q8*x%d2val3 - q81*x%d1val2_d1val3 + q82*x%d1val1_d1val3 + q83*q84 - q85*q86 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q26*(12.0_dp*q89 + 2.0_dp*x%d3val2_d1val3 - q47*x%d2val2_d1val3 - q61*q90 - q87*x%d2val2 + q88*x%d2val2) - q49*q55
      unary%d2val2_d2val3 = -q0*(2.0_dp*q95 + q2*q92 + q45*(-2.0_dp*q52 + q53 + x%d2val2_d1val3) + q66*q93 + q68*(-q14*q66 + q36) - q82*x%d1val2_d1val3 + q91 - q94)
      unary%d3val1_d2val3 = q0*(-24.0_dp*q101*x%d1val1_d1val3 + 12.0_dp*q10*q71 + 6.0_dp*q13*x%d2val1_d1val3 + 6.0_dp*q58*x%d1val1_d2val3 - q0*q97*x%d1val1_d2val3 - q100*x%d2val3 + q103*q23 - q27*q46 - q45*(6.0_dp*q59 - q100*x%d1val3 - q12*q97 + q13*q97 - q96*x%d2val1_d1val3 + x%d3val1_d1val3) - q56*x%d2val1_d1val3 - q57*q85 + q57*q99 - q96*x%d2val1_d2val3 + q97*q98 + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q0*(-2.0_dp*q17*q70 - 3.0_dp*q110*q67 + 2.0_dp*q109*q67 + 2.0_dp*q28*q71 - q101*q107 + q102*q4*q86 + q104*q11 - q104*q8 - q105*q12 + q105*q13 + q106*q40 + q106*q41 + q108*q67 - q111*q26*x%d2val3 + q111*q83 + q2*q76 - q26*q67*x%d1val2_d2val3 - q43*q74*x%d2val3 - q45*(-2.0_dp*q63 + 2.0_dp*q64 + q13*q84 - q44*q74 + q58*x%d1val2_d1val3 - q69*x%d1val1_d1val2_d1val3 + x%d2val1_d1val2_d1val3) + q58*x%d1val2_d2val3 - q65*x%d1val1*x%d1val1_d1val3 - q69*x%d1val1_d1val2_d2val3 + q70*q82 - q78*q85 + q84*q98 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q0*(0.5_dp*q37*q98 + 2.0_dp*q113*q13 + q1*(-2.0_dp*q94 - 8.0_dp*q19*x%d1val2_d1val3 + 4.0_dp*q95 + q35*x%d1val2_d2val3 + q5*q93 + q91) - q105*q17 + q105*q19 + q107*q41 + q108*q84 - q110*q78 - q112*x%d1val1_d1val2_d2val3 - q113*q81 - q26*q37*x%d1val1_d2val3 - q37*q85 + q37*q99 - q39*q46 - q45*(-0.5_dp*q79 - 2.0_dp*q77 + 0.5_dp*q80 - q1*q113 - q112*x%d1val1_d1val2_d1val3 + q19*q84 + x%d1val1_d2val2_d1val3) - q8*q92 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q0*(-24.0_dp*q14*q42*x%d1val2_d1val3*x%d1val3 + 12.0_dp*q28*q93 + 6.0_dp*q109*x%d2val2 + 6.0_dp*q14*q9*x%d1val2_d2val3 - q0*q115*x%d1val2_d2val3 + q103*q48 + q108*q115 - q114*x%d2val2_d2val3 - q116*x%d2val3 - q45*(6.0_dp*q89 - q114*x%d2val2_d1val3 - q115*q17 + q115*q19 - q116*x%d1val3 + x%d3val2_d1val3) - q46*q50 - q75*q86*x%d2val2 - q87*x%d2val2_d1val3 + q88*x%d2val2_d1val3 + x%d3val2_d2val3)
   end function safe_log_self
   
   function log10_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = powm1(ln10)
      q1 = powm1(x%val)
      q2 = q0*q1
      q3 = pow2(x%d1val1)
      q4 = 4.0_dp*q1
      q5 = 4.0_dp*x%d2val1 - q3*q4
      q6 = 0.25_dp*q2
      q7 = powm1(pow2(x%val))
      q8 = q0*q7
      q9 = x%d1val1*x%d1val2
      q10 = q8*x%d1val3
      q11 = pow2(x%d1val2)
      q12 = -q11*q4
      q13 = 4.0_dp*x%d2val2 + q12
      q14 = pow2(x%d1val3)
      q15 = q1*q14
      q16 = q1*x%d1val1
      q17 = 6.0_dp*q16
      q18 = pow3(x%d1val1)
      q19 = 4.0_dp*q7
      q20 = 2.0_dp*x%d3val1 - q17*x%d2val1 + q18*q19
      q21 = 0.5_dp*q2
      q22 = q8*x%d1val2
      q23 = 0.25_dp*q5
      q24 = 8.0_dp*q16
      q25 = q19*q3
      q26 = 4.0_dp*x%d2val1_d1val2 - q24*x%d1val1_d1val2 + q25*x%d1val2
      q27 = 4.0_dp*x%d2val1_d1val3 - q24*x%d1val1_d1val3 + q25*x%d1val3
      q28 = q4*x%d1val1_d1val2
      q29 = 2.0_dp*x%d2val2
      q30 = q12 + q29
      q31 = 2.0_dp*x%d1val1_d2val2 - q16*q30 - q28*x%d1val2
      q32 = q8*x%d1val2_d1val3
      q33 = 2.0_dp*x%d1val3
      q34 = powm1(pow3(x%val))
      q35 = q0*q34
      q36 = 2.0_dp*q1
      q37 = q36*x%d1val3
      q38 = -2.0_dp*q15 + x%d2val3
      q39 = q1*x%d1val2
      q40 = 6.0_dp*q39
      q41 = pow3(x%d1val2)
      q42 = 2.0_dp*x%d3val2 + q19*q41 - q40*x%d2val2
      q43 = 0.25_dp*q10
      q44 = q39*x%d1val2_d1val3
      q45 = q19*x%d1val3
      q46 = -8.0_dp*q44 + q11*q45
      q47 = 0.5_dp*q10
      q48 = q1*x%d1val1_d1val3
      q49 = 6.0_dp*q48
      q50 = q7*x%d1val1
      q51 = 6.0_dp*x%d1val3
      q52 = q50*q51
      q53 = q3*q7
      q54 = q53*x%d1val1_d1val3
      q55 = q18*q34
      q56 = 8.0_dp*x%d1val3
      q57 = q7*x%d1val1_d1val3
      q58 = q57*q9
      q59 = q3*q34
      q60 = 2.0_dp*x%d2val1 - q3*q36
      q61 = 0.5_dp*q1
      q62 = q38*q61
      q63 = 2.0_dp*q16
      q64 = q53*x%d1val3 - q63*x%d1val1_d1val3 + x%d2val1_d1val3
      q65 = pow2(x%d1val1_d1val3)
      q66 = q45*x%d1val1
      q67 = q14*q34
      q68 = 2.0_dp*q67
      q69 = q3*q68 + q36*q65 - q53*x%d2val3 + q63*x%d1val1_d2val3 - q66*x%d1val1_d1val3 - x%d2val1_d2val3
      q70 = q4*x%d1val1_d1val2_d1val3
      q71 = q45*x%d1val2
      q72 = q30*q48
      q73 = q30*x%d1val3
      q74 = q50*q73
      q75 = q1*x%d2val3
      q76 = q36*x%d1val2_d1val3
      q77 = 2.0_dp*q7
      q78 = q77*q9
      q79 = q14*q7
      q80 = 2.0_dp*x%d1val1_d1val2
      q81 = 6.0_dp*q67
      q82 = q1*x%d2val2
      q83 = 6.0_dp*x%d1val2_d1val3
      q84 = q7*x%d1val2
      q85 = q84*x%d2val2
      q86 = q11*q7
      q87 = q34*q41
      q88 = -x%d2val2_d2val3
      q89 = 2.0_dp*q39
      q90 = pow2(x%d1val2_d1val3)
      q91 = q86*x%d2val3
      q92 = q86*x%d1val3
      q93 = 3.0_dp*q16
      q94 = 3.0_dp*x%d2val1
      q95 = q50*q94
      q96 = 4.0_dp*q55
      q97 = q67*x%d1val1
      q98 = q59*x%d1val3
      q99 = q14*powm1(pow4(x%val))
      q100 = 12.0_dp*q99
      q101 = q36*x%d1val1_d1val2
      q102 = q50*q80
      q103 = q45*x%d1val1_d1val2
      q104 = 2.0_dp*x%d1val2
      q105 = 4.0_dp*q67
      q106 = q105*x%d1val1_d1val2
      q107 = q84*x%d2val3
      q108 = x%d1val2_d1val3*x%d1val3
      q109 = 2.0_dp*x%d2val1_d1val2 + q104*q53 - q28*x%d1val1
      q110 = q4*x%d1val2
      q111 = 2.0_dp*q92 - q110*x%d1val2_d1val3 + x%d2val2_d1val3
      q112 = 3.0_dp*q39
      q113 = 3.0_dp*q82
      q114 = 4.0_dp*q87
      unary%val = q0*log(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q5*q6
      unary%d1val1_d1val2 = q2*x%d1val1_d1val2 - q8*q9
      unary%d1val1_d1val3 = -q10*x%d1val1 + q2*x%d1val1_d1val3
      unary%d2val2 = q13*q6
      unary%d1val2_d1val3 = -q10*x%d1val2 + q2*x%d1val2_d1val3
      unary%d2val3 = q2*(-q15 + x%d2val3)
      unary%d3val1 = q20*q21
      unary%d2val1_d1val2 = -q22*q23 + q26*q6
      unary%d2val1_d1val3 = -q10*q23 + q27*q6
      unary%d1val1_d2val2 = q21*q31
      unary%d1val1_d1val2_d1val3 = -q10*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3 - q22*x%d1val1_d1val3 - q32*x%d1val1 + q33*q35*q9
      unary%d1val1_d2val3 = q2*(-q16*q38 - q37*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q21*q42
      unary%d2val2_d1val3 = -q13*q43 + q6*(4.0_dp*x%d2val2_d1val3 + q46)
      unary%d1val2_d2val3 = q2*(-q37*x%d1val2_d1val3 - q38*q39 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q20*q47 + q21*(12.0_dp*q54 + 2.0_dp*x%d3val1_d1val3 - q17*x%d2val1_d1val3 - q49*x%d2val1 + q52*x%d2val1 - q55*q56)
      unary%d2val1_d1val2_d1val3 = -0.25_dp*q22*q27 + 0.5_dp*q35*q5*x%d1val2*x%d1val3 - q23*q32 - q26*q43 + q6*(-8.0_dp*q48*x%d1val1_d1val2 + 4.0_dp*x%d2val1_d1val2_d1val3 + 8.0_dp*q58 - q24*x%d1val1_d1val2_d1val3 + q25*x%d1val2_d1val3 + q50*q56*x%d1val1_d1val2 - q56*q59*x%d1val2)
      unary%d2val1_d2val3 = -q2*(q37*q64 + q60*q62 + q69)
      unary%d1val1_d2val2_d1val3 = q21*(2.0_dp*x%d1val1_d2val2_d1val3 - q16*(2.0_dp*x%d2val2_d1val3 + q46) - q28*x%d1val2_d1val3 - q70*x%d1val2 + q71*x%d1val1_d1val2 - q72 + q74) - q31*q47
      unary%d1val1_d1val2_d2val3 = q2*(-q16*x%d1val2_d2val3 - q37*x%d1val1_d1val2_d1val3 - q39*x%d1val1_d2val3 + q66*x%d1val2_d1val3 + q71*x%d1val1_d1val3 - q75*x%d1val1_d1val2 - q76*x%d1val1_d1val3 + q78*x%d2val3 + q79*q80 - q81*q9 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q21*(12.0_dp*q86*x%d1val2_d1val3 + 2.0_dp*x%d3val2_d1val3 - q40*x%d2val2_d1val3 + q51*q85 - q56*q87 - q82*q83) - q42*q47
      unary%d2val2_d2val3 = -q2*(q11*q68 + q36*q90 + q37*(-2.0_dp*q44 + q92 + x%d2val2_d1val3) + q62*(-q11*q36 + q29) - q71*x%d1val2_d1val3 + q88 + q89*x%d1val2_d2val3 - q91)
      unary%d3val1_d2val3 = q2*(-24.0_dp*q98*x%d1val1_d1val3 - 6.0_dp*q97*x%d2val1 + 12.0_dp*q50*q65 + 6.0_dp*q53*x%d1val1_d2val3 - q1*q94*x%d1val1_d2val3 + q100*q18 - q20*q62 - q37*(6.0_dp*q54 - q48*q94 - q93*x%d2val1_d1val3 + q95*x%d1val3 - q96*x%d1val3 + x%d3val1_d1val3) - q49*x%d2val1_d1val3 + q51*q57*x%d2val1 + q52*x%d2val1_d1val3 - q93*x%d2val1_d2val3 + q95*x%d2val3 - q96*x%d2val3 + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(-0.5_dp*q109*q75 - 3.0_dp*q60*q67*x%d1val2 - 4.0_dp*q98*x%d1val2_d1val3 + 6.0_dp*q3*q99*x%d1val2 - q101*x%d1val1_d2val3 + q102*x%d2val3 + q103*x%d1val1_d1val3 - q104*q59*x%d2val3 + q104*q65*q7 - q106*x%d1val1 + q107*q60 + q108*q60*q77 + q109*q79 + q19*x%d1val1*x%d1val1_d1val3*x%d1val2_d1val3 - q34*q56*q9*x%d1val1_d1val3 - q37*(2.0_dp*q58 - q101*x%d1val1_d1val3 + q102*x%d1val3 - q104*q98 + q53*x%d1val2_d1val3 - q63*x%d1val1_d1val2_d1val3 + x%d2val1_d1val2_d1val3) + q39*q69 + q53*x%d1val2_d2val3 - q60*q61*x%d1val2_d2val3 - q63*x%d1val1_d1val2_d2val3 + q64*q71 - q64*q76 + q66*x%d1val1_d1val2_d1val3 - q70*x%d1val1_d1val3 + q78*x%d1val1_d2val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(0.5_dp*q30*q50*x%d2val3 - q101*x%d1val2_d2val3 + q103*x%d1val2_d1val3 - q106*x%d1val2 + q107*q80 + q111*q33*q50 - q111*q36*x%d1val1_d1val3 + q16*(-2.0_dp*q91 + q105*q11 + q110*x%d1val2_d2val3 + q4*q90 - q56*q84*x%d1val2_d1val3 + q88) - q30*q61*x%d1val1_d2val3 - q30*q97 - q31*q62 - q37*(-0.5_dp*q72 + 0.5_dp*q74 - q111*q16 - q76*x%d1val1_d1val2 + q80*q84*x%d1val3 - q89*x%d1val1_d1val2_d1val3 + x%d1val1_d2val2_d1val3) + q57*q73 - q70*x%d1val2_d1val3 + q71*x%d1val1_d1val2_d1val3 - q89*x%d1val1_d1val2_d2val3 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(-24.0_dp*q108*q11*q34 + 12.0_dp*q84*q90 + 3.0_dp*q107*x%d2val2 + 6.0_dp*q86*x%d1val2_d2val3 - q1*q83*x%d2val2_d1val3 + q100*q41 - q112*x%d2val2_d2val3 - q113*x%d1val2_d2val3 - q114*x%d2val3 - q37*(3.0_dp*q85*x%d1val3 - q112*x%d2val2_d1val3 - q113*x%d1val2_d1val3 - q114*x%d1val3 + q83*q86 + x%d3val2_d1val3) - q42*q62 + q51*q84*x%d2val2_d1val3 + q7*q83*x%d1val3*x%d2val2 - q81*x%d1val2*x%d2val2 + x%d3val2_d2val3)
   end function log10_self
   
   function safe_log10_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = powm1(ln10)
      q1 = powm1(x%val)
      q2 = q0*q1
      q3 = pow2(x%d1val1)
      q4 = 4.0_dp*q1
      q5 = 4.0_dp*x%d2val1 - q3*q4
      q6 = 0.25_dp*q2
      q7 = powm1(pow2(x%val))
      q8 = q0*q7
      q9 = x%d1val1*x%d1val2
      q10 = q8*x%d1val3
      q11 = pow2(x%d1val2)
      q12 = -q11*q4
      q13 = 4.0_dp*x%d2val2 + q12
      q14 = pow2(x%d1val3)
      q15 = q1*q14
      q16 = q1*x%d1val1
      q17 = 6.0_dp*q16
      q18 = pow3(x%d1val1)
      q19 = 4.0_dp*q7
      q20 = 2.0_dp*x%d3val1 - q17*x%d2val1 + q18*q19
      q21 = 0.5_dp*q2
      q22 = q8*x%d1val2
      q23 = 0.25_dp*q5
      q24 = 8.0_dp*q16
      q25 = q19*q3
      q26 = 4.0_dp*x%d2val1_d1val2 - q24*x%d1val1_d1val2 + q25*x%d1val2
      q27 = 4.0_dp*x%d2val1_d1val3 - q24*x%d1val1_d1val3 + q25*x%d1val3
      q28 = q4*x%d1val1_d1val2
      q29 = 2.0_dp*x%d2val2
      q30 = q12 + q29
      q31 = 2.0_dp*x%d1val1_d2val2 - q16*q30 - q28*x%d1val2
      q32 = q8*x%d1val2_d1val3
      q33 = 2.0_dp*x%d1val3
      q34 = powm1(pow3(x%val))
      q35 = q0*q34
      q36 = 2.0_dp*q1
      q37 = q36*x%d1val3
      q38 = -2.0_dp*q15 + x%d2val3
      q39 = q1*x%d1val2
      q40 = 6.0_dp*q39
      q41 = pow3(x%d1val2)
      q42 = 2.0_dp*x%d3val2 + q19*q41 - q40*x%d2val2
      q43 = 0.25_dp*q10
      q44 = q39*x%d1val2_d1val3
      q45 = q19*x%d1val3
      q46 = -8.0_dp*q44 + q11*q45
      q47 = 0.5_dp*q10
      q48 = q1*x%d1val1_d1val3
      q49 = 6.0_dp*q48
      q50 = q7*x%d1val1
      q51 = 6.0_dp*x%d1val3
      q52 = q50*q51
      q53 = q3*q7
      q54 = q53*x%d1val1_d1val3
      q55 = q18*q34
      q56 = 8.0_dp*x%d1val3
      q57 = q7*x%d1val1_d1val3
      q58 = q57*q9
      q59 = q3*q34
      q60 = 2.0_dp*x%d2val1 - q3*q36
      q61 = 0.5_dp*q1
      q62 = q38*q61
      q63 = 2.0_dp*q16
      q64 = q53*x%d1val3 - q63*x%d1val1_d1val3 + x%d2val1_d1val3
      q65 = pow2(x%d1val1_d1val3)
      q66 = q45*x%d1val1
      q67 = q14*q34
      q68 = 2.0_dp*q67
      q69 = q3*q68 + q36*q65 - q53*x%d2val3 + q63*x%d1val1_d2val3 - q66*x%d1val1_d1val3 - x%d2val1_d2val3
      q70 = q4*x%d1val1_d1val2_d1val3
      q71 = q45*x%d1val2
      q72 = q30*q48
      q73 = q30*x%d1val3
      q74 = q50*q73
      q75 = q1*x%d2val3
      q76 = q36*x%d1val2_d1val3
      q77 = 2.0_dp*q7
      q78 = q77*q9
      q79 = q14*q7
      q80 = 2.0_dp*x%d1val1_d1val2
      q81 = 6.0_dp*q67
      q82 = q1*x%d2val2
      q83 = 6.0_dp*x%d1val2_d1val3
      q84 = q7*x%d1val2
      q85 = q84*x%d2val2
      q86 = q11*q7
      q87 = q34*q41
      q88 = -x%d2val2_d2val3
      q89 = 2.0_dp*q39
      q90 = pow2(x%d1val2_d1val3)
      q91 = q86*x%d2val3
      q92 = q86*x%d1val3
      q93 = 3.0_dp*q16
      q94 = 3.0_dp*x%d2val1
      q95 = q50*q94
      q96 = 4.0_dp*q55
      q97 = q67*x%d1val1
      q98 = q59*x%d1val3
      q99 = q14*powm1(pow4(x%val))
      q100 = 12.0_dp*q99
      q101 = q36*x%d1val1_d1val2
      q102 = q50*q80
      q103 = q45*x%d1val1_d1val2
      q104 = 2.0_dp*x%d1val2
      q105 = 4.0_dp*q67
      q106 = q105*x%d1val1_d1val2
      q107 = q84*x%d2val3
      q108 = x%d1val2_d1val3*x%d1val3
      q109 = 2.0_dp*x%d2val1_d1val2 + q104*q53 - q28*x%d1val1
      q110 = q4*x%d1val2
      q111 = 2.0_dp*q92 - q110*x%d1val2_d1val3 + x%d2val2_d1val3
      q112 = 3.0_dp*q39
      q113 = 3.0_dp*q82
      q114 = 4.0_dp*q87
      unary%val = q0*safe_log(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q5*q6
      unary%d1val1_d1val2 = q2*x%d1val1_d1val2 - q8*q9
      unary%d1val1_d1val3 = -q10*x%d1val1 + q2*x%d1val1_d1val3
      unary%d2val2 = q13*q6
      unary%d1val2_d1val3 = -q10*x%d1val2 + q2*x%d1val2_d1val3
      unary%d2val3 = q2*(-q15 + x%d2val3)
      unary%d3val1 = q20*q21
      unary%d2val1_d1val2 = -q22*q23 + q26*q6
      unary%d2val1_d1val3 = -q10*q23 + q27*q6
      unary%d1val1_d2val2 = q21*q31
      unary%d1val1_d1val2_d1val3 = -q10*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3 - q22*x%d1val1_d1val3 - q32*x%d1val1 + q33*q35*q9
      unary%d1val1_d2val3 = q2*(-q16*q38 - q37*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q21*q42
      unary%d2val2_d1val3 = -q13*q43 + q6*(4.0_dp*x%d2val2_d1val3 + q46)
      unary%d1val2_d2val3 = q2*(-q37*x%d1val2_d1val3 - q38*q39 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q20*q47 + q21*(12.0_dp*q54 + 2.0_dp*x%d3val1_d1val3 - q17*x%d2val1_d1val3 - q49*x%d2val1 + q52*x%d2val1 - q55*q56)
      unary%d2val1_d1val2_d1val3 = -0.25_dp*q22*q27 + 0.5_dp*q35*q5*x%d1val2*x%d1val3 - q23*q32 - q26*q43 + q6*(-8.0_dp*q48*x%d1val1_d1val2 + 4.0_dp*x%d2val1_d1val2_d1val3 + 8.0_dp*q58 - q24*x%d1val1_d1val2_d1val3 + q25*x%d1val2_d1val3 + q50*q56*x%d1val1_d1val2 - q56*q59*x%d1val2)
      unary%d2val1_d2val3 = -q2*(q37*q64 + q60*q62 + q69)
      unary%d1val1_d2val2_d1val3 = q21*(2.0_dp*x%d1val1_d2val2_d1val3 - q16*(2.0_dp*x%d2val2_d1val3 + q46) - q28*x%d1val2_d1val3 - q70*x%d1val2 + q71*x%d1val1_d1val2 - q72 + q74) - q31*q47
      unary%d1val1_d1val2_d2val3 = q2*(-q16*x%d1val2_d2val3 - q37*x%d1val1_d1val2_d1val3 - q39*x%d1val1_d2val3 + q66*x%d1val2_d1val3 + q71*x%d1val1_d1val3 - q75*x%d1val1_d1val2 - q76*x%d1val1_d1val3 + q78*x%d2val3 + q79*q80 - q81*q9 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q21*(12.0_dp*q86*x%d1val2_d1val3 + 2.0_dp*x%d3val2_d1val3 - q40*x%d2val2_d1val3 + q51*q85 - q56*q87 - q82*q83) - q42*q47
      unary%d2val2_d2val3 = -q2*(q11*q68 + q36*q90 + q37*(-2.0_dp*q44 + q92 + x%d2val2_d1val3) + q62*(-q11*q36 + q29) - q71*x%d1val2_d1val3 + q88 + q89*x%d1val2_d2val3 - q91)
      unary%d3val1_d2val3 = q2*(-24.0_dp*q98*x%d1val1_d1val3 - 6.0_dp*q97*x%d2val1 + 12.0_dp*q50*q65 + 6.0_dp*q53*x%d1val1_d2val3 - q1*q94*x%d1val1_d2val3 + q100*q18 - q20*q62 - q37*(6.0_dp*q54 - q48*q94 - q93*x%d2val1_d1val3 + q95*x%d1val3 - q96*x%d1val3 + x%d3val1_d1val3) - q49*x%d2val1_d1val3 + q51*q57*x%d2val1 + q52*x%d2val1_d1val3 - q93*x%d2val1_d2val3 + q95*x%d2val3 - q96*x%d2val3 + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(-0.5_dp*q109*q75 - 3.0_dp*q60*q67*x%d1val2 - 4.0_dp*q98*x%d1val2_d1val3 + 6.0_dp*q3*q99*x%d1val2 - q101*x%d1val1_d2val3 + q102*x%d2val3 + q103*x%d1val1_d1val3 - q104*q59*x%d2val3 + q104*q65*q7 - q106*x%d1val1 + q107*q60 + q108*q60*q77 + q109*q79 + q19*x%d1val1*x%d1val1_d1val3*x%d1val2_d1val3 - q34*q56*q9*x%d1val1_d1val3 - q37*(2.0_dp*q58 - q101*x%d1val1_d1val3 + q102*x%d1val3 - q104*q98 + q53*x%d1val2_d1val3 - q63*x%d1val1_d1val2_d1val3 + x%d2val1_d1val2_d1val3) + q39*q69 + q53*x%d1val2_d2val3 - q60*q61*x%d1val2_d2val3 - q63*x%d1val1_d1val2_d2val3 + q64*q71 - q64*q76 + q66*x%d1val1_d1val2_d1val3 - q70*x%d1val1_d1val3 + q78*x%d1val1_d2val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(0.5_dp*q30*q50*x%d2val3 - q101*x%d1val2_d2val3 + q103*x%d1val2_d1val3 - q106*x%d1val2 + q107*q80 + q111*q33*q50 - q111*q36*x%d1val1_d1val3 + q16*(-2.0_dp*q91 + q105*q11 + q110*x%d1val2_d2val3 + q4*q90 - q56*q84*x%d1val2_d1val3 + q88) - q30*q61*x%d1val1_d2val3 - q30*q97 - q31*q62 - q37*(-0.5_dp*q72 + 0.5_dp*q74 - q111*q16 - q76*x%d1val1_d1val2 + q80*q84*x%d1val3 - q89*x%d1val1_d1val2_d1val3 + x%d1val1_d2val2_d1val3) + q57*q73 - q70*x%d1val2_d1val3 + q71*x%d1val1_d1val2_d1val3 - q89*x%d1val1_d1val2_d2val3 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(-24.0_dp*q108*q11*q34 + 12.0_dp*q84*q90 + 3.0_dp*q107*x%d2val2 + 6.0_dp*q86*x%d1val2_d2val3 - q1*q83*x%d2val2_d1val3 + q100*q41 - q112*x%d2val2_d2val3 - q113*x%d1val2_d2val3 - q114*x%d2val3 - q37*(3.0_dp*q85*x%d1val3 - q112*x%d2val2_d1val3 - q113*x%d1val2_d1val3 - q114*x%d1val3 + q83*q86 + x%d3val2_d1val3) - q42*q62 + q51*q84*x%d2val2_d1val3 + q7*q83*x%d1val3*x%d2val2 - q81*x%d1val2*x%d2val2 + x%d3val2_d2val3)
   end function safe_log10_self
   
   function log2_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = powm1(log(2.0_dp))
      q1 = powm1(x%val)
      q2 = q0*q1
      q3 = pow2(x%d1val1)
      q4 = 4.0_dp*q1
      q5 = 4.0_dp*x%d2val1 - q3*q4
      q6 = 0.25_dp*q2
      q7 = powm1(pow2(x%val))
      q8 = q0*q7
      q9 = x%d1val1*x%d1val2
      q10 = q8*x%d1val3
      q11 = pow2(x%d1val2)
      q12 = -q11*q4
      q13 = 4.0_dp*x%d2val2 + q12
      q14 = pow2(x%d1val3)
      q15 = q1*q14
      q16 = q1*x%d1val1
      q17 = 6.0_dp*q16
      q18 = pow3(x%d1val1)
      q19 = 4.0_dp*q7
      q20 = 2.0_dp*x%d3val1 - q17*x%d2val1 + q18*q19
      q21 = 0.5_dp*q2
      q22 = q8*x%d1val2
      q23 = 0.25_dp*q5
      q24 = 8.0_dp*q16
      q25 = q19*q3
      q26 = 4.0_dp*x%d2val1_d1val2 - q24*x%d1val1_d1val2 + q25*x%d1val2
      q27 = 4.0_dp*x%d2val1_d1val3 - q24*x%d1val1_d1val3 + q25*x%d1val3
      q28 = q4*x%d1val1_d1val2
      q29 = 2.0_dp*x%d2val2
      q30 = q12 + q29
      q31 = 2.0_dp*x%d1val1_d2val2 - q16*q30 - q28*x%d1val2
      q32 = q8*x%d1val2_d1val3
      q33 = 2.0_dp*x%d1val3
      q34 = powm1(pow3(x%val))
      q35 = q0*q34
      q36 = 2.0_dp*q1
      q37 = q36*x%d1val3
      q38 = -2.0_dp*q15 + x%d2val3
      q39 = q1*x%d1val2
      q40 = 6.0_dp*q39
      q41 = pow3(x%d1val2)
      q42 = 2.0_dp*x%d3val2 + q19*q41 - q40*x%d2val2
      q43 = 0.25_dp*q10
      q44 = q39*x%d1val2_d1val3
      q45 = q19*x%d1val3
      q46 = -8.0_dp*q44 + q11*q45
      q47 = 0.5_dp*q10
      q48 = q1*x%d1val1_d1val3
      q49 = 6.0_dp*q48
      q50 = q7*x%d1val1
      q51 = 6.0_dp*x%d1val3
      q52 = q50*q51
      q53 = q3*q7
      q54 = q53*x%d1val1_d1val3
      q55 = q18*q34
      q56 = 8.0_dp*x%d1val3
      q57 = q7*x%d1val1_d1val3
      q58 = q57*q9
      q59 = q3*q34
      q60 = 2.0_dp*x%d2val1 - q3*q36
      q61 = 0.5_dp*q1
      q62 = q38*q61
      q63 = 2.0_dp*q16
      q64 = q53*x%d1val3 - q63*x%d1val1_d1val3 + x%d2val1_d1val3
      q65 = pow2(x%d1val1_d1val3)
      q66 = q45*x%d1val1
      q67 = q14*q34
      q68 = 2.0_dp*q67
      q69 = q3*q68 + q36*q65 - q53*x%d2val3 + q63*x%d1val1_d2val3 - q66*x%d1val1_d1val3 - x%d2val1_d2val3
      q70 = q4*x%d1val1_d1val2_d1val3
      q71 = q45*x%d1val2
      q72 = q30*q48
      q73 = q30*x%d1val3
      q74 = q50*q73
      q75 = q1*x%d2val3
      q76 = q36*x%d1val2_d1val3
      q77 = 2.0_dp*q7
      q78 = q77*q9
      q79 = q14*q7
      q80 = 2.0_dp*x%d1val1_d1val2
      q81 = 6.0_dp*q67
      q82 = q1*x%d2val2
      q83 = 6.0_dp*x%d1val2_d1val3
      q84 = q7*x%d1val2
      q85 = q84*x%d2val2
      q86 = q11*q7
      q87 = q34*q41
      q88 = -x%d2val2_d2val3
      q89 = 2.0_dp*q39
      q90 = pow2(x%d1val2_d1val3)
      q91 = q86*x%d2val3
      q92 = q86*x%d1val3
      q93 = 3.0_dp*q16
      q94 = 3.0_dp*x%d2val1
      q95 = q50*q94
      q96 = 4.0_dp*q55
      q97 = q67*x%d1val1
      q98 = q59*x%d1val3
      q99 = q14*powm1(pow4(x%val))
      q100 = 12.0_dp*q99
      q101 = q36*x%d1val1_d1val2
      q102 = q50*q80
      q103 = q45*x%d1val1_d1val2
      q104 = 2.0_dp*x%d1val2
      q105 = 4.0_dp*q67
      q106 = q105*x%d1val1_d1val2
      q107 = q84*x%d2val3
      q108 = x%d1val2_d1val3*x%d1val3
      q109 = 2.0_dp*x%d2val1_d1val2 + q104*q53 - q28*x%d1val1
      q110 = q4*x%d1val2
      q111 = 2.0_dp*q92 - q110*x%d1val2_d1val3 + x%d2val2_d1val3
      q112 = 3.0_dp*q39
      q113 = 3.0_dp*q82
      q114 = 4.0_dp*q87
      unary%val = q0*log(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q5*q6
      unary%d1val1_d1val2 = q2*x%d1val1_d1val2 - q8*q9
      unary%d1val1_d1val3 = -q10*x%d1val1 + q2*x%d1val1_d1val3
      unary%d2val2 = q13*q6
      unary%d1val2_d1val3 = -q10*x%d1val2 + q2*x%d1val2_d1val3
      unary%d2val3 = q2*(-q15 + x%d2val3)
      unary%d3val1 = q20*q21
      unary%d2val1_d1val2 = -q22*q23 + q26*q6
      unary%d2val1_d1val3 = -q10*q23 + q27*q6
      unary%d1val1_d2val2 = q21*q31
      unary%d1val1_d1val2_d1val3 = -q10*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3 - q22*x%d1val1_d1val3 - q32*x%d1val1 + q33*q35*q9
      unary%d1val1_d2val3 = q2*(-q16*q38 - q37*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q21*q42
      unary%d2val2_d1val3 = -q13*q43 + q6*(4.0_dp*x%d2val2_d1val3 + q46)
      unary%d1val2_d2val3 = q2*(-q37*x%d1val2_d1val3 - q38*q39 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q20*q47 + q21*(12.0_dp*q54 + 2.0_dp*x%d3val1_d1val3 - q17*x%d2val1_d1val3 - q49*x%d2val1 + q52*x%d2val1 - q55*q56)
      unary%d2val1_d1val2_d1val3 = -0.25_dp*q22*q27 + 0.5_dp*q35*q5*x%d1val2*x%d1val3 - q23*q32 - q26*q43 + q6*(-8.0_dp*q48*x%d1val1_d1val2 + 4.0_dp*x%d2val1_d1val2_d1val3 + 8.0_dp*q58 - q24*x%d1val1_d1val2_d1val3 + q25*x%d1val2_d1val3 + q50*q56*x%d1val1_d1val2 - q56*q59*x%d1val2)
      unary%d2val1_d2val3 = -q2*(q37*q64 + q60*q62 + q69)
      unary%d1val1_d2val2_d1val3 = q21*(2.0_dp*x%d1val1_d2val2_d1val3 - q16*(2.0_dp*x%d2val2_d1val3 + q46) - q28*x%d1val2_d1val3 - q70*x%d1val2 + q71*x%d1val1_d1val2 - q72 + q74) - q31*q47
      unary%d1val1_d1val2_d2val3 = q2*(-q16*x%d1val2_d2val3 - q37*x%d1val1_d1val2_d1val3 - q39*x%d1val1_d2val3 + q66*x%d1val2_d1val3 + q71*x%d1val1_d1val3 - q75*x%d1val1_d1val2 - q76*x%d1val1_d1val3 + q78*x%d2val3 + q79*q80 - q81*q9 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q21*(12.0_dp*q86*x%d1val2_d1val3 + 2.0_dp*x%d3val2_d1val3 - q40*x%d2val2_d1val3 + q51*q85 - q56*q87 - q82*q83) - q42*q47
      unary%d2val2_d2val3 = -q2*(q11*q68 + q36*q90 + q37*(-2.0_dp*q44 + q92 + x%d2val2_d1val3) + q62*(-q11*q36 + q29) - q71*x%d1val2_d1val3 + q88 + q89*x%d1val2_d2val3 - q91)
      unary%d3val1_d2val3 = q2*(-24.0_dp*q98*x%d1val1_d1val3 - 6.0_dp*q97*x%d2val1 + 12.0_dp*q50*q65 + 6.0_dp*q53*x%d1val1_d2val3 - q1*q94*x%d1val1_d2val3 + q100*q18 - q20*q62 - q37*(6.0_dp*q54 - q48*q94 - q93*x%d2val1_d1val3 + q95*x%d1val3 - q96*x%d1val3 + x%d3val1_d1val3) - q49*x%d2val1_d1val3 + q51*q57*x%d2val1 + q52*x%d2val1_d1val3 - q93*x%d2val1_d2val3 + q95*x%d2val3 - q96*x%d2val3 + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(-0.5_dp*q109*q75 - 3.0_dp*q60*q67*x%d1val2 - 4.0_dp*q98*x%d1val2_d1val3 + 6.0_dp*q3*q99*x%d1val2 - q101*x%d1val1_d2val3 + q102*x%d2val3 + q103*x%d1val1_d1val3 - q104*q59*x%d2val3 + q104*q65*q7 - q106*x%d1val1 + q107*q60 + q108*q60*q77 + q109*q79 + q19*x%d1val1*x%d1val1_d1val3*x%d1val2_d1val3 - q34*q56*q9*x%d1val1_d1val3 - q37*(2.0_dp*q58 - q101*x%d1val1_d1val3 + q102*x%d1val3 - q104*q98 + q53*x%d1val2_d1val3 - q63*x%d1val1_d1val2_d1val3 + x%d2val1_d1val2_d1val3) + q39*q69 + q53*x%d1val2_d2val3 - q60*q61*x%d1val2_d2val3 - q63*x%d1val1_d1val2_d2val3 + q64*q71 - q64*q76 + q66*x%d1val1_d1val2_d1val3 - q70*x%d1val1_d1val3 + q78*x%d1val1_d2val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(0.5_dp*q30*q50*x%d2val3 - q101*x%d1val2_d2val3 + q103*x%d1val2_d1val3 - q106*x%d1val2 + q107*q80 + q111*q33*q50 - q111*q36*x%d1val1_d1val3 + q16*(-2.0_dp*q91 + q105*q11 + q110*x%d1val2_d2val3 + q4*q90 - q56*q84*x%d1val2_d1val3 + q88) - q30*q61*x%d1val1_d2val3 - q30*q97 - q31*q62 - q37*(-0.5_dp*q72 + 0.5_dp*q74 - q111*q16 - q76*x%d1val1_d1val2 + q80*q84*x%d1val3 - q89*x%d1val1_d1val2_d1val3 + x%d1val1_d2val2_d1val3) + q57*q73 - q70*x%d1val2_d1val3 + q71*x%d1val1_d1val2_d1val3 - q89*x%d1val1_d1val2_d2val3 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(-24.0_dp*q108*q11*q34 + 12.0_dp*q84*q90 + 3.0_dp*q107*x%d2val2 + 6.0_dp*q86*x%d1val2_d2val3 - q1*q83*x%d2val2_d1val3 + q100*q41 - q112*x%d2val2_d2val3 - q113*x%d1val2_d2val3 - q114*x%d2val3 - q37*(3.0_dp*q85*x%d1val3 - q112*x%d2val2_d1val3 - q113*x%d1val2_d1val3 - q114*x%d1val3 + q83*q86 + x%d3val2_d1val3) - q42*q62 + q51*q84*x%d2val2_d1val3 + q7*q83*x%d1val3*x%d2val2 - q81*x%d1val2*x%d2val2 + x%d3val2_d2val3)
   end function log2_self
   
   function sin_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = sin(x%val)
      q1 = cos(x%val)
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4 = q1*x%d1val3
      q5 = q1*x%d2val1
      q6 = pow2(x%d1val1)
      q7 = q0*q6
      q8 = q1*x%d1val1_d1val2
      q9 = q0*x%d1val1
      q10 = q9*x%d1val2
      q11 = q1*x%d1val1_d1val3
      q12 = q1*x%d2val2
      q13 = pow2(x%d1val2)
      q14 = q0*q13
      q15 = q1*x%d1val2_d1val3
      q16 = q0*x%d1val2
      q17 = q16*x%d1val3
      q18 = q1*x%d2val3
      q19 = pow2(x%d1val3)
      q20 = q0*q19
      q21 = q1*x%d3val1
      q22 = 3.0_dp*q9
      q23 = pow3(x%d1val1)
      q24 = q1*q23
      q25 = q1*x%d2val1_d1val2
      q26 = 2.0_dp*q9
      q27 = q3*q6
      q28 = q0*x%d1val3
      q29 = q1*x%d1val1_d2val2
      q30 = 2.0_dp*q16
      q31 = q0*x%d2val2
      q32 = q1*q13
      q33 = x%d1val2*x%d1val3
      q34 = q1*x%d1val1_d2val3
      q35 = 2.0_dp*q28
      q36 = q0*x%d2val3
      q37 = q1*q19
      q38 = q36 + q37
      q39 = 3.0_dp*q16
      q40 = pow3(x%d1val2)
      q41 = q1*x%d1val2_d2val3
      q42 = q0*x%d2val1
      q43 = 3.0_dp*q42
      q44 = q2*x%d1val3
      q45 = 3.0_dp*x%d2val1
      q46 = 3.0_dp*q6
      q47 = 2.0_dp*q0
      q48 = q47*x%d1val1_d1val2
      q49 = 2.0_dp*q44
      q50 = 2.0_dp*q2
      q51 = q50*x%d1val2
      q52 = q3*x%d2val1
      q53 = pow2(x%d1val1_d1val3)
      q54 = 4.0_dp*x%d1val1_d1val3
      q55 = 2.0_dp*q3
      q56 = q55*x%d1val3
      q57 = q0*x%d2val2_d1val3
      q58 = q3*x%d1val2_d1val3
      q59 = q4*x%d2val2
      q60 = q14*x%d1val3
      q61 = q47*x%d1val2_d1val3
      q62 = x%d1val2*x%d2val3
      q63 = q20*x%d1val2
      q64 = 3.0_dp*q31
      q65 = q3*x%d1val3
      q66 = 3.0_dp*x%d2val2
      q67 = 3.0_dp*q13
      q68 = pow2(x%d1val2_d1val3)
      q69 = 6.0_dp*x%d1val1_d1val3
      q70 = 6.0_dp*q2
      q71 = q4*x%d2val1
      q72 = q20*x%d1val1
      q73 = q7*x%d1val3
      q74 = x%d1val1_d1val2*x%d2val3
      q75 = 4.0_dp*x%d1val1_d1val2_d1val3
      q76 = q4*x%d1val1_d1val2
      q77 = 2.0_dp*x%d1val2_d1val3
      q78 = 2.0_dp*x%d1val1_d1val2
      q79 = 4.0_dp*q0
      q80 = 4.0_dp*x%d2val3
      q81 = 6.0_dp*x%d1val2_d1val3
      q82 = 6.0_dp*q3
      unary%val = q0
      unary%d1val1 = q2
      unary%d1val2 = q3
      unary%d1val3 = q4
      unary%d2val1 = q5 - q7
      unary%d1val1_d1val2 = -q10 + q8
      unary%d1val1_d1val3 = q11 - q9*x%d1val3
      unary%d2val2 = q12 - q14
      unary%d1val2_d1val3 = q15 - q17
      unary%d2val3 = q18 - q20
      unary%d3val1 = q21 - q22*x%d2val1 - q24
      unary%d2val1_d1val2 = -q16*x%d2val1 + q25 - q26*x%d1val1_d1val2 - q27
      unary%d2val1_d1val3 = q1*x%d2val1_d1val3 - q26*x%d1val1_d1val3 - q28*x%d2val1 - q4*q6
      unary%d1val1_d2val2 = -0.0625_dp*x%d1val1*(16.0_dp*q31 + 16.0_dp*q32) + q29 - q30*x%d1val1_d1val2
      unary%d1val1_d1val2_d1val3 = q1*x%d1val1_d1val2_d1val3 - q16*x%d1val1_d1val3 - q2*q33 - q28*x%d1val1_d1val2 - q9*x%d1val2_d1val3
      unary%d1val1_d2val3 = q34 - q35*x%d1val1_d1val3 - q38*x%d1val1
      unary%d3val2 = -q1*q40 + q1*x%d3val2 - q39*x%d2val2
      unary%d2val2_d1val3 = q1*x%d2val2_d1val3 - q13*q4 - q28*x%d2val2 - q30*x%d1val2_d1val3
      unary%d1val2_d2val3 = -q35*x%d1val2_d1val3 - q38*x%d1val2 + q41
      unary%d3val1_d1val3 = q1*x%d3val1_d1val3 - q11*q46 - q22*x%d2val1_d1val3 + q23*q28 - q28*x%d3val1 - q43*x%d1val1_d1val3 - q44*q45
      unary%d2val1_d1val2_d1val3 = q1*x%d2val1_d1val2_d1val3 - q15*q6 - q16*x%d2val1_d1val3 - q26*x%d1val1_d1val2_d1val3 - q28*x%d2val1_d1val2 + q33*q7 - q42*x%d1val2_d1val3 - q48*x%d1val1_d1val3 - q49*x%d1val1_d1val2 - q51*x%d1val1_d1val3 - q52*x%d1val3
      unary%d2val1_d2val3 = q1*x%d2val1_d2val3 - q18*q6 - q19*q5 + q19*q7 - q26*x%d1val1_d2val3 - q35*x%d2val1_d1val3 - q36*x%d2val1 - q44*q54 - q47*q53
      unary%d1val1_d2val2_d1val3 = 4.0_dp*x%d1val1*(-0.25_dp*q57 - 0.25_dp*q59 - 0.5_dp*q58 + 0.25_dp*q60) + q1*x%d1val1_d2val2_d1val3 - q28*x%d1val1_d2val2 - q30*x%d1val1_d1val2_d1val3 - q48*x%d1val2_d1val3 + q54*(-0.25_dp*q31 - 0.25_dp*q32) - q56*x%d1val1_d1val2
      unary%d1val1_d1val2_d2val3 = q1*x%d1val1_d1val2_d2val3 - q16*x%d1val1_d2val3 - q19*q8 - q2*q62 - q35*x%d1val1_d1val2_d1val3 - q36*x%d1val1_d1val2 - q49*x%d1val2_d1val3 - q56*x%d1val1_d1val3 - q61*x%d1val1_d1val3 + q63*x%d1val1 - q9*x%d1val2_d2val3
      unary%d3val2_d1val3 = q1*x%d3val2_d1val3 - q15*q67 + q28*q40 - q28*x%d3val2 - q39*x%d2val2_d1val3 - q64*x%d1val2_d1val3 - q65*q66
      unary%d2val2_d2val3 = -4.0_dp*q58*x%d1val3 + q1*x%d2val2_d2val3 - q12*q19 - q13*q18 + q14*q19 - q30*x%d1val2_d2val3 - q35*x%d2val2_d1val3 - q36*x%d2val2 - q47*q68
      unary%d3val1_d2val3 = -q0*q69*x%d2val1_d1val3 + q1*x%d3val1_d2val3 - q19*q21 + q19*q24 - q2*q45*x%d2val3 - q22*x%d2val1_d2val3 + q23*q36 - q34*q46 - q35*x%d3val1_d1val3 - q36*x%d3val1 - q43*x%d1val1_d2val3 + q45*q72 - q53*q70 - q69*q71 + q69*q73 - q70*x%d1val3*x%d2val1_d1val3
      unary%d2val1_d1val2_d2val3 = -q0*q54*x%d1val1_d1val2_d1val3 + q1*x%d2val1_d1val2_d2val3 + q10*q54*x%d1val3 - q16*x%d2val1_d2val3 - q19*q25 + q19*q27 - q2*q54*x%d1val2_d1val3 - q26*x%d1val1_d1val2_d2val3 - q35*x%d2val1_d1val2_d1val3 - q36*x%d2val1_d1val2 - q41*q6 - q42*x%d1val2_d2val3 - q44*q75 - q48*x%d1val1_d2val3 - q50*q74 - q51*x%d1val1_d2val3 - q52*x%d2val3 - q53*q55 - q54*q76 - q56*x%d2val1_d1val3 - q61*x%d2val1_d1val3 + q62*q7 + q63*x%d2val1 - q71*q77 + q72*q78 + q73*q77
      unary%d1val1_d2val2_d2val3 = -0.25_dp*x%d1val1*(-16.0_dp*q17*x%d1val2_d1val3 - 4.0_dp*q13*q37 - 4.0_dp*q20*x%d2val2 + 8.0_dp*q1*q68 + 8.0_dp*q3*x%d1val2_d2val3 + 8.0_dp*q4*x%d2val2_d1val3 + q12*q80 - q14*q80 + q79*x%d2val2_d2val3) - 0.25_dp*x%d1val1_d2val3*(4.0_dp*q31 + 4.0_dp*q32) - 0.5_dp*x%d1val1_d1val3*(-4.0_dp*q60 + 4.0_dp*q57 + 4.0_dp*q59 + 8.0_dp*q58) - 4.0_dp*q76*x%d1val2_d1val3 + q1*x%d1val1_d2val2_d2val3 - q19*q29 - q30*x%d1val1_d1val2_d2val3 - q35*x%d1val1_d2val2_d1val3 - q36*x%d1val1_d2val2 - q48*x%d1val2_d2val3 - q55*q74 + q63*q78 - q65*q75 - q79*x%d1val1_d1val2_d1val3*x%d1val2_d1val3
      unary%d3val2_d2val3 = q1*x%d3val2_d2val3 - q3*q66*x%d2val3 - q35*x%d3val2_d1val3 + q36*q40 - q36*x%d3val2 + q37*q40 - q37*x%d3val2 - q39*x%d2val2_d2val3 - q41*q67 - q57*q81 - q59*q81 + q60*q81 + q63*q66 - q64*x%d1val2_d2val3 - q68*q82 - q82*x%d1val3*x%d2val2_d1val3
   end function sin_self
   
   function cos_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = cos(x%val)
      q1 = sin(x%val)
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4 = q1*x%d1val3
      q5 = q1*x%d2val1
      q6 = pow2(x%d1val1)
      q7 = q0*q6
      q8 = q1*x%d1val1_d1val2
      q9 = q0*x%d1val1
      q10 = q9*x%d1val2
      q11 = q1*x%d1val1_d1val3
      q12 = q1*x%d2val2
      q13 = pow2(x%d1val2)
      q14 = q0*q13
      q15 = q1*x%d1val2_d1val3
      q16 = q0*x%d1val2
      q17 = q16*x%d1val3
      q18 = q1*x%d2val3
      q19 = pow2(x%d1val3)
      q20 = q0*q19
      q21 = q1*x%d3val1
      q22 = pow3(x%d1val1)
      q23 = q1*q22
      q24 = 3.0_dp*q9
      q25 = q1*x%d2val1_d1val2
      q26 = 2.0_dp*q9
      q27 = q3*q6
      q28 = q0*x%d1val3
      q29 = q1*x%d1val1_d2val2
      q30 = 2.0_dp*q16
      q31 = q0*x%d2val2
      q32 = q1*q13
      q33 = x%d1val2*x%d1val3
      q34 = q1*x%d1val1_d2val3
      q35 = 2.0_dp*q28
      q36 = q0*x%d2val3
      q37 = q1*q19
      q38 = q36 - q37
      q39 = pow3(x%d1val2)
      q40 = 3.0_dp*q16
      q41 = q1*x%d1val2_d2val3
      q42 = q0*x%d2val1
      q43 = 3.0_dp*q42
      q44 = q2*x%d1val3
      q45 = 3.0_dp*x%d2val1
      q46 = 3.0_dp*q6
      q47 = 2.0_dp*q0
      q48 = q47*x%d1val1_d1val2
      q49 = q3*x%d2val1
      q50 = 2.0_dp*q44
      q51 = 2.0_dp*q2
      q52 = q51*x%d1val2
      q53 = pow2(x%d1val1_d1val3)
      q54 = 4.0_dp*x%d1val1_d1val3
      q55 = 2.0_dp*q3
      q56 = q55*x%d1val3
      q57 = q0*x%d2val2_d1val3
      q58 = q3*x%d1val2_d1val3
      q59 = q4*x%d2val2
      q60 = q14*x%d1val3
      q61 = q47*x%d1val2_d1val3
      q62 = x%d1val2*x%d2val3
      q63 = q20*x%d1val2
      q64 = 3.0_dp*q31
      q65 = q3*x%d1val3
      q66 = 3.0_dp*x%d2val2
      q67 = 3.0_dp*q13
      q68 = pow2(x%d1val2_d1val3)
      q69 = 6.0_dp*x%d2val1_d1val3
      q70 = 6.0_dp*x%d1val1_d1val3
      q71 = q4*x%d2val1
      q72 = q20*x%d1val1
      q73 = q7*x%d1val3
      q74 = x%d1val1_d1val2*x%d2val3
      q75 = 4.0_dp*x%d1val1_d1val2_d1val3
      q76 = q4*x%d1val1_d1val2
      q77 = 2.0_dp*x%d1val2_d1val3
      q78 = 2.0_dp*x%d1val1_d1val2
      q79 = 4.0_dp*q0
      q80 = 4.0_dp*x%d2val3
      q81 = 6.0_dp*x%d1val2_d1val3
      q82 = 6.0_dp*q3
      unary%val = q0
      unary%d1val1 = -q2
      unary%d1val2 = -q3
      unary%d1val3 = -q4
      unary%d2val1 = -q5 - q7
      unary%d1val1_d1val2 = -q10 - q8
      unary%d1val1_d1val3 = -q11 - q9*x%d1val3
      unary%d2val2 = -q12 - q14
      unary%d1val2_d1val3 = -q15 - q17
      unary%d2val3 = -q18 - q20
      unary%d3val1 = -q21 + q23 - q24*x%d2val1
      unary%d2val1_d1val2 = -q16*x%d2val1 - q25 - q26*x%d1val1_d1val2 + q27
      unary%d2val1_d1val3 = -q1*x%d2val1_d1val3 - q26*x%d1val1_d1val3 - q28*x%d2val1 + q4*q6
      unary%d1val1_d2val2 = -0.0625_dp*x%d1val1*(-16.0_dp*q32 + 16.0_dp*q31) - q29 - q30*x%d1val1_d1val2
      unary%d1val1_d1val2_d1val3 = -q1*x%d1val1_d1val2_d1val3 - q16*x%d1val1_d1val3 + q2*q33 - q28*x%d1val1_d1val2 - q9*x%d1val2_d1val3
      unary%d1val1_d2val3 = -q34 - q35*x%d1val1_d1val3 - q38*x%d1val1
      unary%d3val2 = q1*q39 - q1*x%d3val2 - q40*x%d2val2
      unary%d2val2_d1val3 = -q1*x%d2val2_d1val3 + q13*q4 - q28*x%d2val2 - q30*x%d1val2_d1val3
      unary%d1val2_d2val3 = -q35*x%d1val2_d1val3 - q38*x%d1val2 - q41
      unary%d3val1_d1val3 = -q1*x%d3val1_d1val3 + q11*q46 + q22*q28 - q24*x%d2val1_d1val3 - q28*x%d3val1 - q43*x%d1val1_d1val3 + q44*q45
      unary%d2val1_d1val2_d1val3 = -q1*x%d2val1_d1val2_d1val3 + q15*q6 - q16*x%d2val1_d1val3 - q26*x%d1val1_d1val2_d1val3 - q28*x%d2val1_d1val2 + q33*q7 - q42*x%d1val2_d1val3 - q48*x%d1val1_d1val3 + q49*x%d1val3 + q50*x%d1val1_d1val2 + q52*x%d1val1_d1val3
      unary%d2val1_d2val3 = -q1*x%d2val1_d2val3 + q18*q6 + q19*q5 + q19*q7 - q26*x%d1val1_d2val3 - q35*x%d2val1_d1val3 - q36*x%d2val1 + q44*q54 - q47*q53
      unary%d1val1_d2val2_d1val3 = 4.0_dp*x%d1val1*(-0.25_dp*q57 + 0.25_dp*q59 + 0.25_dp*q60 + 0.5_dp*q58) - q1*x%d1val1_d2val2_d1val3 - q28*x%d1val1_d2val2 - q30*x%d1val1_d1val2_d1val3 - q48*x%d1val2_d1val3 + q54*(-0.25_dp*q31 + 0.25_dp*q32) + q56*x%d1val1_d1val2
      unary%d1val1_d1val2_d2val3 = -q1*x%d1val1_d1val2_d2val3 - q16*x%d1val1_d2val3 + q19*q8 + q2*q62 - q35*x%d1val1_d1val2_d1val3 - q36*x%d1val1_d1val2 + q50*x%d1val2_d1val3 + q56*x%d1val1_d1val3 - q61*x%d1val1_d1val3 + q63*x%d1val1 - q9*x%d1val2_d2val3
      unary%d3val2_d1val3 = -q1*x%d3val2_d1val3 + q15*q67 + q28*q39 - q28*x%d3val2 - q40*x%d2val2_d1val3 - q64*x%d1val2_d1val3 + q65*q66
      unary%d2val2_d2val3 = 4.0_dp*q58*x%d1val3 - q1*x%d2val2_d2val3 + q12*q19 + q13*q18 + q14*q19 - q30*x%d1val2_d2val3 - q35*x%d2val2_d1val3 - q36*x%d2val2 - q47*q68
      unary%d3val1_d2val3 = 6.0_dp*q2*q53 - q0*q69*x%d1val1_d1val3 - q1*x%d3val1_d2val3 + q19*q21 - q19*q23 + q2*q45*x%d2val3 + q22*q36 - q24*x%d2val1_d2val3 + q34*q46 - q35*x%d3val1_d1val3 - q36*x%d3val1 - q43*x%d1val1_d2val3 + q44*q69 + q45*q72 + q70*q71 + q70*q73
      unary%d2val1_d1val2_d2val3 = -q0*q54*x%d1val1_d1val2_d1val3 - q1*x%d2val1_d1val2_d2val3 + q10*q54*x%d1val3 - q16*x%d2val1_d2val3 + q19*q25 - q19*q27 + q2*q54*x%d1val2_d1val3 - q26*x%d1val1_d1val2_d2val3 - q35*x%d2val1_d1val2_d1val3 - q36*x%d2val1_d1val2 + q41*q6 - q42*x%d1val2_d2val3 + q44*q75 - q48*x%d1val1_d2val3 + q49*x%d2val3 + q51*q74 + q52*x%d1val1_d2val3 + q53*q55 + q54*q76 + q56*x%d2val1_d1val3 - q61*x%d2val1_d1val3 + q62*q7 + q63*x%d2val1 + q71*q77 + q72*q78 + q73*q77
      unary%d1val1_d2val2_d2val3 = 0.25_dp*x%d1val1*(-4.0_dp*q13*q37 + 16.0_dp*q17*x%d1val2_d1val3 + 4.0_dp*q20*x%d2val2 + 8.0_dp*q1*q68 + 8.0_dp*q3*x%d1val2_d2val3 + 8.0_dp*q4*x%d2val2_d1val3 + q12*q80 + q14*q80 - q79*x%d2val2_d2val3) + 0.25_dp*x%d1val1_d2val3*(-4.0_dp*q31 + 4.0_dp*q32) + 0.5_dp*x%d1val1_d1val3*(-4.0_dp*q57 + 4.0_dp*q59 + 4.0_dp*q60 + 8.0_dp*q58) + 4.0_dp*q76*x%d1val2_d1val3 - q1*x%d1val1_d2val2_d2val3 + q19*q29 - q30*x%d1val1_d1val2_d2val3 - q35*x%d1val1_d2val2_d1val3 - q36*x%d1val1_d2val2 - q48*x%d1val2_d2val3 + q55*q74 + q63*q78 + q65*q75 - q79*x%d1val1_d1val2_d1val3*x%d1val2_d1val3
      unary%d3val2_d2val3 = -q1*x%d3val2_d2val3 + q3*q66*x%d2val3 - q35*x%d3val2_d1val3 + q36*q39 - q36*x%d3val2 - q37*q39 + q37*x%d3val2 - q40*x%d2val2_d2val3 + q41*q67 - q57*q81 + q59*q81 + q60*q81 + q63*q66 - q64*x%d1val2_d2val3 + q68*q82 + q82*x%d1val3*x%d2val2_d1val3
   end function cos_self
   
   function tan_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q3 = q2*x%d1val1
      q4 = q2*x%d1val2
      q5 = q2*x%d1val3
      q6 = 0.125_dp*q1 + 0.125_dp
      q7 = pow2(x%d1val1)
      q8 = 16.0_dp*q0
      q9 = q7*q8
      q10 = 8.0_dp*x%d2val1 + q9
      q11 = q2*x%d1val1_d1val2
      q12 = 2.0_dp*q0
      q13 = q12*x%d1val2
      q14 = q2*x%d1val1_d1val3
      q15 = q12*x%d1val3
      q16 = pow2(x%d1val2)
      q17 = 8.0_dp*x%d2val2 + q16*q8
      q18 = q2*x%d1val2_d1val3
      q19 = pow2(x%d1val3)
      q20 = 0.03125_dp*q1 + 0.03125_dp
      q21 = q0*x%d1val1
      q22 = 192.0_dp*x%d2val1
      q23 = pow3(x%d1val1)
      q24 = 128.0_dp*q1
      q25 = 64.0_dp*q2
      q26 = 32.0_dp*x%d3val1 + q21*q22 + q23*q24 + q23*q25
      q27 = 0.25_dp*q0
      q28 = q10*q4
      q29 = 32.0_dp*q21
      q30 = 16.0_dp*q7
      q31 = 8.0_dp*x%d2val1_d1val2 + q29*x%d1val1_d1val2 + q30*q4
      q32 = q27*q5
      q33 = 8.0_dp*x%d2val1_d1val3 + q29*x%d1val1_d1val3 + q30*q5
      q34 = q0*x%d1val2
      q35 = 128.0_dp*x%d1val1_d1val2
      q36 = q1*q16
      q37 = q16*q2
      q38 = 16.0_dp*q37 + 32.0_dp*q36 + q8*x%d2val2
      q39 = 4.0_dp*x%d1val1
      q40 = 32.0_dp*x%d1val1_d2val2 + q34*q35 + q38*q39
      q41 = 2.0_dp*q1
      q42 = 2.0_dp + q41
      q43 = q34*x%d1val1_d1val3
      q44 = q3*x%d1val3
      q45 = q44*x%d1val2
      q46 = 4.0_dp*q0
      q47 = q46*x%d1val3
      q48 = q0*x%d2val3
      q49 = q19*q41
      q50 = q19*q2
      q51 = q48 + q49 + q50
      q52 = 2.0_dp*q51
      q53 = 192.0_dp*q34
      q54 = pow3(x%d1val2)
      q55 = 32.0_dp*x%d3val2 + q24*q54 + q25*q54 + q53*x%d2val2
      q56 = 32.0_dp*x%d1val2_d1val3
      q57 = 16.0_dp*q5
      q58 = q0*q5
      q59 = 0.0625_dp*q58
      q60 = q21*x%d2val1_d1val3
      q61 = q0*x%d1val1_d1val3
      q62 = q1*q7
      q63 = q62*x%d1val1_d1val3
      q64 = q14*q7
      q65 = 384.0_dp*q58
      q66 = q0*x%d1val2_d1val3
      q67 = 0.25_dp*q1 + 0.25_dp
      q68 = q10*q67
      q69 = 0.5_dp*q1
      q70 = q4*x%d1val3
      q71 = 32.0_dp*x%d1val1_d1val2
      q72 = q3*x%d1val2
      q73 = q72*x%d1val1_d1val3
      q74 = q7*q70
      q75 = 0.5_dp + q69
      q76 = 2.0_dp*x%d2val1 + q46*q7
      q77 = 8.0_dp*q0
      q78 = q77*x%d1val1
      q79 = 4.0_dp*q5
      q80 = 2.0_dp*x%d2val1_d1val3 + q7*q79 + q78*x%d1val1_d1val3
      q81 = pow2(x%d1val1_d1val3)
      q82 = 16.0_dp*q44
      q83 = 4.0_dp*q2
      q84 = q7*q83
      q85 = q50*q77
      q86 = 2.0_dp*x%d2val1_d2val3 + q7*q85 + q77*q81 + q78*x%d1val1_d2val3 + q82*x%d1val1_d1val3 + q84*x%d2val3
      q87 = 4.0_dp*x%d1val1_d1val3
      q88 = x%d1val2*x%d1val2_d1val3
      q89 = q0*x%d1val2_d2val3
      q90 = q77*x%d1val1_d1val2_d1val3
      q91 = q77*x%d1val1_d1val3
      q92 = x%d1val1_d2val3*x%d1val2
      q93 = 8.0_dp*q1
      q94 = x%d1val1*x%d1val2
      q95 = 16.0_dp*q1
      q96 = x%d1val2_d1val3*x%d1val3
      q97 = x%d1val1_d1val3*x%d1val3
      q98 = q19*x%d1val1_d1val2
      q99 = q19*pow3(q0)
      q100 = 4.0_dp*x%d2val3
      q101 = 8.0_dp*q44
      q102 = 8.0_dp*q4
      q103 = q19*q3
      q104 = q66*x%d2val2
      q105 = q36*x%d1val2_d1val3
      q106 = q16*q18
      q107 = 2.0_dp*x%d2val2_d2val3
      q108 = x%d1val2*x%d1val2_d2val3
      q109 = pow2(x%d1val2_d1val3)
      q110 = 16.0_dp*x%d1val2_d1val3
      q111 = 4.0_dp*q37
      q112 = 2.0_dp*x%d2val2_d1val3
      q113 = 2.0_dp*x%d2val2
      q114 = 12.0_dp*x%d2val1
      q115 = 48.0_dp*q1
      q116 = 24.0_dp*q3
      q117 = q3*x%d2val3
      q118 = 24.0_dp*x%d2val1
      q119 = 24.0_dp*q23
      q120 = q2*q48
      q121 = 144.0_dp*q5
      q122 = q19*pow2(q2)
      q123 = q115*q50
      q124 = 8.0_dp*q2
      q125 = q77*x%d1val1_d1val2
      q126 = q8*x%d1val1_d1val2_d1val3
      q127 = q57*x%d1val1_d1val2
      q128 = q8*q98
      q129 = q102*q7
      q130 = q5*x%d1val2_d1val3
      q131 = q76*x%d1val2
      q132 = q19*q4
      q133 = q4*x%d2val3
      q134 = q8*x%d1val1_d1val2
      q135 = 4.0_dp*x%d2val1_d1val2 + q129 + q134*x%d1val1
      q136 = q80*x%d1val3
      q137 = 4.0_dp*q4
      q138 = q102*x%d1val1_d1val2
      q139 = q46*x%d2val2
      q140 = 8.0_dp*q36 + q111 + q139
      q141 = 12.0_dp*q16
      q142 = q0*q112 + q113*q5 + q137*x%d1val2_d1val3 + q141*q58 + q88*q93
      q143 = 2.0_dp*x%d1val1
      q144 = 12.0_dp*q37
      q145 = 24.0_dp*q36
      q146 = 12.0_dp*q34
      q147 = 12.0_dp*x%d2val2
      q148 = 24.0_dp*q4
      q149 = 24.0_dp*x%d2val2
      q150 = 24.0_dp*q54
      unary%val = q0
      unary%d1val1 = q3
      unary%d1val2 = q4
      unary%d1val3 = q5
      unary%d2val1 = q10*q6
      unary%d1val1_d1val2 = q11 + q13*q3
      unary%d1val1_d1val3 = q14 + q15*q3
      unary%d2val2 = q17*q6
      unary%d1val2_d1val3 = q15*q4 + q18
      unary%d2val3 = q2*(q12*q19 + x%d2val3)
      unary%d3val1 = q20*q26
      unary%d2val1_d1val2 = q27*q28 + q31*q6
      unary%d2val1_d1val3 = q10*q32 + q33*q6
      unary%d1val1_d2val2 = q20*q40
      unary%d1val1_d1val2_d1val3 = 4.0_dp*q1*q45 + q12*q5*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3 + q21*q42*x%d1val2_d1val3 + q42*q43 + q42*q45
      unary%d1val1_d2val3 = q2*(q47*x%d1val1_d1val3 + q52*x%d1val1 + x%d1val1_d2val3)
      unary%d3val2 = q20*q55
      unary%d2val2_d1val3 = q17*q32 + q6*(8.0_dp*x%d2val2_d1val3 + q16*q57 + q34*q56)
      unary%d1val2_d2val3 = q2*(q47*x%d1val2_d1val3 + q52*x%d1val2 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q20*(192.0_dp*q60 + 192.0_dp*q64 + 32.0_dp*x%d3val1_d1val3 + 384.0_dp*q63 + q22*q44 + q22*q61 + q23*q65) + q26*q59
      unary%d2val1_d1val2_d1val3 = q28*q69*x%d1val3 + q31*q32 + q33*q34*q67 + q6*(32.0_dp*q0*q74 + 32.0_dp*q73 + 8.0_dp*x%d2val1_d1val2_d1val3 + q18*q30 + q29*x%d1val1_d1val2_d1val3 + q44*q71 + q61*q71) + q66*q68 + q68*q70
      unary%d2val1_d2val3 = q75*(q47*q80 + q52*q76 + q86)
      unary%d1val1_d2val2_d1val3 = q20*(128.0_dp*q34*x%d1val1_d1val2_d1val3 + 32.0_dp*x%d1val1_d2val2_d1val3 + q35*q66 + q35*q70 + q38*q87 + q39*(64.0_dp*q1*q88 + 96.0_dp*q16*q58 + q4*q56 + q57*x%d2val2 + q8*x%d2val2_d1val3)) + q40*q59
      unary%d1val1_d1val2_d2val3 = q75*(16.0_dp*q94*q99 + 2.0_dp*x%d1val1_d1val2_d2val3 + 32.0_dp*q103*q34 + 4.0_dp*q11*q19 + 4.0_dp*q48*x%d1val1_d1val2 + q100*q72 + q101*x%d1val2_d1val3 + q102*q97 + q39*q89 + q46*q92 + q90*x%d1val3 + q91*x%d1val2_d1val3 + q93*q94*x%d2val3 + q93*q98 + q95*q96*x%d1val1 + q95*q97*x%d1val2)
      unary%d3val2_d1val3 = q20*(192.0_dp*q104 + 192.0_dp*q106 + 192.0_dp*q70*x%d2val2 + 32.0_dp*x%d3val2_d1val3 + 384.0_dp*q105 + q53*x%d2val2_d1val3 + q54*q65) + q55*q59
      unary%d2val2_d2val3 = q75*(q107 + q108*q77 + q109*q77 + q110*q70 + q111*x%d2val3 + q16*q85 + q47*(q112 + q16*q79 + q77*q88) + q52*(q113 + q16*q46))
      unary%d3val1_d2val3 = q75*(12.0_dp*q2*q7*x%d1val1_d2val3 + 12.0_dp*q21*x%d2val1_d2val3 + 2.0_dp*x%d3val1_d2val3 + 24.0_dp*q61*x%d2val1_d1val3 + 24.0_dp*q62*x%d1val1_d2val3 + q0*q103*q118 + q0*q114*x%d1val1_d2val3 + q114*q117 + q115*q81*x%d1val1 + q116*q81 + q116*x%d1val3*x%d2val1_d1val3 + q118*q5*x%d1val1_d1val3 + q119*q120 + q119*q122 + q121*q61*q7 + q123*q23 + q47*(12.0_dp*q60 + 12.0_dp*q64 + 2.0_dp*x%d3val1_d1val3 + 24.0_dp*q63 + q114*q44 + q114*q61 + q119*q58) + q51*(4.0_dp*x%d3val1 + q118*q21 + q124*q23 + q23*q95))
      unary%d2val1_d1val2_d2val3 = q75*(2.0_dp*q133*q76 + 2.0_dp*x%d2val1_d1val2_d2val3 + 32.0_dp*q43*q44 + 8.0_dp*q117*x%d1val1_d1val2 + 8.0_dp*q122*q7*x%d1val2 + 8.0_dp*q131*q99 + 8.0_dp*q3*q92 + q1*q100*q131 + q102*q81 + q110*q3*x%d1val1_d1val3 + q12*q76*x%d1val2_d2val3 + q125*x%d1val1_d2val3 + q126*x%d1val1_d1val3 + q127*x%d1val1_d1val3 + q128*q3 + q129*q48 + q13*q86 + q130*q9 + q132*q7*q95 + q132*q76*q8 + q135*q48 + q135*q49 + q135*q50 + q136*q137 + q136*q93*x%d1val2 + q46*q80*x%d1val2_d1val3 + q47*(2.0_dp*x%d2val1_d1val2_d1val3 + 4.0_dp*q18*q7 + 8.0_dp*q73 + q101*x%d1val1_d1val2 + q74*q77 + q78*x%d1val1_d1val2_d1val3 + q91*x%d1val1_d1val2) + q76*q79*x%d1val2_d1val3 + q76*q93*q96 + q78*x%d1val1_d1val2_d2val3 + q82*x%d1val1_d1val2_d1val3 + q84*x%d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q75*(16.0_dp*q70*x%d1val1_d1val2_d1val3 + 2.0_dp*x%d1val1_d2val2_d2val3 + q125*x%d1val2_d2val3 + q126*x%d1val2_d1val3 + q127*x%d1val2_d1val3 + q128*q4 + q138*x%d2val3 + q140*x%d1val1_d2val3 + q142*q87 + q143*(48.0_dp*q66*q70 + q0*q107 + q108*q93 + q109*q83 + q109*q93 + q113*q2*x%d2val3 + q122*q141 + q137*x%d1val2_d2val3 + q139*q50 + q144*q48 + q145*q50 + q79*x%d2val2_d1val3) + q47*(2.0_dp*x%d1val1_d2val2_d1val3 + q125*x%d1val2_d1val3 + q138*x%d1val3 + q140*x%d1val1_d1val3 + q142*q143 + q90*x%d1val2) + q51*(4.0_dp*x%d1val1_d2val2 + q134*x%d1val2 + q140*q143) + q77*x%d1val1_d1val2_d2val3*x%d1val2)
      unary%d3val2_d2val3 = q75*(2.0_dp*x%d3val2_d2val3 + 24.0_dp*q66*x%d2val2_d1val3 + q0*q132*q149 + q109*q115*x%d1val2 + q109*q148 + q120*q150 + q121*q16*q66 + q122*q150 + q123*q54 + q130*q149 + q133*q147 + q144*x%d1val2_d2val3 + q145*x%d1val2_d2val3 + q146*x%d2val2_d2val3 + q147*q89 + q148*x%d1val3*x%d2val2_d1val3 + q47*(12.0_dp*q104 + 12.0_dp*q106 + 2.0_dp*x%d3val2_d1val3 + 24.0_dp*q105 + q146*x%d2val2_d1val3 + q147*q70 + q150*q58) + q51*(4.0_dp*x%d3val2 + q124*q54 + q149*q34 + q54*q95))
   end function tan_self
   
   function sinpi_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q1 = sin(q0)
      q2 = cos(q0)
      q3 = pi*q2
      q4 = q3*x%d1val2
      q5 = q3*x%d1val3
      q6 = 16.0_dp*q2
      q7 = q6*x%d2val1
      q8 = pow2(x%d1val1)
      q9 = pi*q1
      q10 = 16.0_dp*q9
      q11 = 0.0625_dp*pi
      q12 = pow2(pi)
      q13 = q1*q12
      q14 = x%d1val1*x%d1val2
      q15 = q13*x%d1val3
      q16 = pow2(x%d1val2)
      q17 = q2*x%d2val3
      q18 = pow2(x%d1val3)
      q19 = 64.0_dp*q2
      q20 = q9*x%d1val1
      q21 = 192.0_dp*q20
      q22 = pow3(x%d1val1)
      q23 = q12*q19
      q24 = 0.015625_dp*pi
      q25 = 32.0_dp*q20
      q26 = q10*x%d2val1
      q27 = q12*q6
      q28 = q8*x%d1val2
      q29 = q27*x%d1val3
      q30 = q9*x%d1val2
      q31 = 128.0_dp*x%d1val1_d1val2
      q32 = 16.0_dp*q1
      q33 = q16*q3
      q34 = 16.0_dp*q33 + q32*x%d2val2
      q35 = pi*x%d1val1
      q36 = 4.0_dp*q35
      q37 = x%d1val1*x%d1val2_d1val3
      q38 = x%d1val1_d1val3*x%d1val2
      q39 = pow3(pi)
      q40 = q39*x%d1val3
      q41 = q14*q40
      q42 = q2*x%d1val1_d2val3
      q43 = q9*x%d1val1_d1val3
      q44 = 2.0_dp*x%d1val3
      q45 = q1*x%d2val3
      q46 = pi*(q18*q3 + q45)
      q47 = 192.0_dp*q30
      q48 = pow3(x%d1val2)
      q49 = 32.0_dp*x%d1val2_d1val3
      q50 = q10*x%d1val3
      q51 = q2*x%d1val2_d2val3
      q52 = q9*x%d1val2_d1val3
      q53 = 192.0_dp*x%d2val1
      q54 = 64.0_dp*q9*x%d1val3
      q55 = q12*q2
      q56 = q55*x%d1val3
      q57 = q56*x%d1val1
      q58 = 192.0_dp*q55
      q59 = q8*x%d1val1_d1val3
      q60 = q1*q40
      q61 = 64.0_dp*q60
      q62 = 32.0_dp*x%d1val1_d1val2
      q63 = q10*x%d1val2
      q64 = q12*x%d1val2
      q65 = q64*x%d1val3
      q66 = q8*x%d1val2_d1val3
      q67 = 4.0_dp*q2
      q68 = 8.0_dp*q9
      q69 = q68*x%d1val1
      q70 = q68*x%d1val3
      q71 = 4.0_dp*pi
      q72 = q45*q71
      q73 = pow2(x%d1val1_d1val3)
      q74 = q29*x%d1val1
      q75 = q12*q17
      q76 = 4.0_dp*q75
      q77 = q12*q18
      q78 = q67*q77
      q79 = 4.0_dp*q1
      q80 = q18*q39
      q81 = q79*q80
      q82 = 0.25_dp*pi
      q83 = q56*x%d1val2
      q84 = q71*x%d1val1_d1val3
      q85 = q5*x%d2val2
      q86 = q12*q16
      q87 = q86*x%d1val3
      q88 = pi*q45
      q89 = q9*x%d1val1_d1val2_d1val3
      q90 = q44*q55
      q91 = q18*x%d1val1_d1val2
      q92 = q1*q80
      q93 = x%d1val3*x%d2val2
      q94 = q16*x%d1val2_d1val3
      q95 = pow2(x%d1val2_d1val3)
      q96 = 8.0_dp*q2
      q97 = 48.0_dp*x%d2val1_d1val3
      q98 = q9*x%d2val1
      q99 = 8.0_dp*q88
      q100 = 24.0_dp*x%d1val1*x%d2val1
      q101 = 48.0_dp*q55
      q102 = x%d1val3*x%d2val1
      q103 = q12*q8
      q104 = q77*q96
      q105 = 8.0_dp*q39
      q106 = q105*q45
      q107 = 48.0_dp*q60
      q108 = q18*pow4(pi)
      q109 = q108*q96
      q110 = 0.125_dp*pi
      q111 = q75*x%d1val1_d1val2
      q112 = q64*q96
      q113 = x%d1val2*x%d2val1
      q114 = 4.0_dp*q45
      q115 = q79*x%d2val2
      q116 = 24.0_dp*x%d2val2
      q117 = q101*x%d1val2
      q118 = q116*x%d1val2
      unary%val = q1
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q4
      unary%d1val3 = q5
      unary%d2val1 = q11*(-q10*q8 + q7)
      unary%d1val1_d1val2 = -q13*q14 + q3*x%d1val1_d1val2
      unary%d1val1_d1val3 = -q15*x%d1val1 + q3*x%d1val1_d1val3
      unary%d2val2 = q11*(-q10*q16 + q6*x%d2val2)
      unary%d1val2_d1val3 = -q15*x%d1val2 + q3*x%d1val2_d1val3
      unary%d2val3 = pi*(q17 - q18*q9)
      unary%d3val1 = q24*(q19*x%d3val1 - q21*x%d2val1 - q22*q23)
      unary%d2val1_d1val2 = q11*(-q25*x%d1val1_d1val2 - q26*x%d1val2 - q27*q28 + q6*x%d2val1_d1val2)
      unary%d2val1_d1val3 = q11*(-q25*x%d1val1_d1val3 - q26*x%d1val3 - q29*q8 + q6*x%d2val1_d1val3)
      unary%d1val1_d2val2 = q24*(q19*x%d1val1_d2val2 - q30*q31 - q34*q36)
      unary%d1val1_d1val2_d1val3 = -q13*q37 - q13*q38 - q15*x%d1val1_d1val2 - q2*q41 + q3*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = pi*(q42 - q43*q44 - q46*x%d1val1)
      unary%d3val2 = q24*(q19*x%d3val2 - q23*q48 - q47*x%d2val2)
      unary%d2val2_d1val3 = q11*(-q16*q29 - q30*q49 - q50*x%d2val2 + q6*x%d2val2_d1val3)
      unary%d1val2_d2val3 = pi*(-q44*q52 - q46*x%d1val2 + q51)
      unary%d3val1_d1val3 = q24*(q19*x%d3val1_d1val3 - q21*x%d2val1_d1val3 + q22*q61 - q43*q53 - q53*q57 - q54*x%d3val1 - q58*q59)
      unary%d2val1_d1val2_d1val3 = q11*(-32.0_dp*q14*q55*x%d1val1_d1val3 - q25*x%d1val1_d1val2_d1val3 - q26*x%d1val2_d1val3 - q27*q66 + q28*q32*q40 - q43*q62 - q50*x%d2val1_d1val2 - q57*q62 + q6*x%d2val1_d1val2_d1val3 - q63*x%d2val1_d1val3 - q65*q7)
      unary%d2val1_d2val3 = -q82*(-q67*x%d2val1_d2val3 + q68*q73 + q69*x%d1val1_d2val3 + q70*x%d2val1_d1val3 + q72*x%d2val1 + q74*x%d1val1_d1val3 + q76*q8 + q78*x%d2val1 - q8*q81)
      unary%d1val1_d2val2_d1val3 = q24*(-128.0_dp*q30*x%d1val1_d1val2_d1val3 + q19*x%d1val1_d2val2_d1val3 - q31*q52 - q31*q83 - q34*q84 - q36*(16.0_dp*q85 - q32*q87 + q32*x%d2val2_d1val3 + q4*q49) - q54*x%d1val1_d2val2)
      unary%d1val1_d1val2_d2val3 = pi*(-2.0_dp*q43*x%d1val2_d1val3 - q14*q75 + q14*q92 + q2*x%d1val1_d1val2_d2val3 - q20*x%d1val2_d2val3 - q30*x%d1val1_d2val3 - q37*q90 - q38*q90 - q44*q89 - q55*q91 - q88*x%d1val1_d1val2)
      unary%d3val2_d1val3 = q24*(-192.0_dp*q52*x%d2val2 + q19*x%d3val2_d1val3 - q47*x%d2val2_d1val3 + q48*q61 - q54*x%d3val2 - q58*q93*x%d1val2 - q58*q94)
      unary%d2val2_d2val3 = -q82*(q16*q76 - q16*q81 + q29*x%d1val2*x%d1val2_d1val3 - q67*x%d2val2_d2val3 + q68*q95 + q68*x%d1val2*x%d1val2_d2val3 + q70*x%d2val2_d1val3 + q72*x%d2val2 + q78*x%d2val2)
      unary%d3val1_d2val3 = -q110*(24.0_dp*q103*q42 + 24.0_dp*q20*x%d2val1_d2val3 + 24.0_dp*q98*x%d1val1_d2val3 + q100*q75 - q100*q92 + q101*q102*x%d1val1_d1val3 + q101*q73*x%d1val1 + q104*x%d3val1 - q106*q22 - q107*q59 - q109*q22 + q43*q97 + q50*x%d3val1_d1val3 + q57*q97 - q96*x%d3val1_d2val3 + q99*x%d3val1)
      unary%d2val1_d1val2_d2val3 = -q82*(-8.0_dp*q60*q66 + 4.0_dp*q103*q51 + 4.0_dp*q30*x%d2val1_d2val3 + 4.0_dp*q98*x%d1val2_d2val3 + 8.0_dp*q111*x%d1val1 + 8.0_dp*q12*q14*q42 - q1*q105*q91*x%d1val1 + q10*x%d1val1_d1val2_d1val3*x%d1val1_d1val3 + q102*q12*q96*x%d1val2_d1val3 - q108*q28*q67 + q112*q73 + q112*x%d1val3*x%d2val1_d1val3 + q113*q76 - q113*q81 - q114*q28*q39 + q27*q37*x%d1val1_d1val3 + q29*x%d1val1_d1val2*x%d1val1_d1val3 - q32*q41*x%d1val1_d1val3 - q67*x%d2val1_d1val2_d2val3 + q68*x%d1val1_d1val2*x%d1val1_d2val3 + q68*x%d1val2_d1val3*x%d2val1_d1val3 + q69*x%d1val1_d1val2_d2val3 + q70*x%d2val1_d1val2_d1val3 + q72*x%d2val1_d1val2 + q74*x%d1val1_d1val2_d1val3 + q78*x%d2val1_d1val2)
      unary%d1val1_d2val2_d2val3 = -q110*(16.0_dp*q111*x%d1val2 + 2.0_dp*pi*x%d1val1_d2val3*(4.0_dp*q33 + q115) + 2.0_dp*q35*(8.0_dp*pi*q51*x%d1val2 + 8.0_dp*q3*q95 + 8.0_dp*q5*x%d2val2_d1val3 - q114*q86 - q115*q77 - q16*q67*q80 + q17*q71*x%d2val2 - q32*q65*x%d1val2_d1val3 + q79*x%d2val2_d2val3) + 32.0_dp*q83*x%d1val1_d1val2_d1val3 + q10*x%d1val1_d1val2*x%d1val2_d2val3 + q104*x%d1val1_d2val2 - q32*q39*q91*x%d1val2 + q49*q56*x%d1val1_d1val2 + q49*q89 + q50*x%d1val1_d2val2_d1val3 + q63*x%d1val1_d1val2_d2val3 + q84*(4.0_dp*q85 + 8.0_dp*q4*x%d1val2_d1val3 - q79*q87 + q79*x%d2val2_d1val3) - q96*x%d1val1_d2val2_d2val3 + q99*x%d1val1_d2val2)
      unary%d3val2_d2val3 = -q110*(24.0_dp*q30*x%d2val2_d2val3 + 24.0_dp*q51*q86 + 48.0_dp*q52*x%d2val2_d1val3 + q101*q93*x%d1val2_d1val3 + q104*x%d3val2 - q106*q48 - q107*q94 - q109*q48 + q116*q9*x%d1val2_d2val3 + q117*q95 + q117*x%d1val3*x%d2val2_d1val3 + q118*q75 - q118*q92 + q50*x%d3val2_d1val3 - q96*x%d3val2_d2val3 + q99*x%d3val2)
   end function sinpi_self
   
   function cospi_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q1 = cos(q0)
      q2 = sin(q0)
      q3 = pi*q2
      q4 = q3*x%d1val2
      q5 = q3*x%d1val3
      q6 = 16.0_dp*q2
      q7 = q6*x%d2val1
      q8 = pow2(x%d1val1)
      q9 = pi*q1
      q10 = 16.0_dp*q9
      q11 = 0.0625_dp*pi
      q12 = pow2(pi)
      q13 = q1*q12
      q14 = x%d1val1*x%d1val2
      q15 = q13*x%d1val3
      q16 = pow2(x%d1val2)
      q17 = q2*x%d2val3
      q18 = pow2(x%d1val3)
      q19 = 64.0_dp*q2
      q20 = q9*x%d1val1
      q21 = 192.0_dp*q20
      q22 = pow3(x%d1val1)
      q23 = q12*q19
      q24 = 0.015625_dp*pi
      q25 = 32.0_dp*q20
      q26 = q10*x%d2val1
      q27 = q12*q6
      q28 = q8*x%d1val2
      q29 = q27*x%d1val3
      q30 = q9*x%d1val2
      q31 = 128.0_dp*x%d1val1_d1val2
      q32 = 16.0_dp*q1
      q33 = q16*q3
      q34 = 16.0_dp*q33 - q32*x%d2val2
      q35 = pi*x%d1val1
      q36 = 4.0_dp*q35
      q37 = x%d1val1*x%d1val2_d1val3
      q38 = x%d1val1_d1val3*x%d1val2
      q39 = pow3(pi)
      q40 = q39*x%d1val3
      q41 = q14*q40
      q42 = q2*x%d1val1_d2val3
      q43 = q9*x%d1val1_d1val3
      q44 = 2.0_dp*x%d1val3
      q45 = q1*x%d2val3
      q46 = pi*(q18*q3 - q45)
      q47 = 192.0_dp*q30
      q48 = pow3(x%d1val2)
      q49 = 32.0_dp*x%d1val2_d1val3
      q50 = q10*x%d1val3
      q51 = q2*x%d1val2_d2val3
      q52 = q9*x%d1val2_d1val3
      q53 = 192.0_dp*x%d2val1
      q54 = 64.0_dp*q9*x%d1val3
      q55 = q12*q2
      q56 = q55*x%d1val3
      q57 = q56*x%d1val1
      q58 = 192.0_dp*q55
      q59 = q8*x%d1val1_d1val3
      q60 = q1*q40
      q61 = 64.0_dp*q60
      q62 = 32.0_dp*x%d1val1_d1val2
      q63 = q10*x%d1val2
      q64 = q12*x%d1val2
      q65 = q64*x%d1val3
      q66 = q8*x%d1val2_d1val3
      q67 = 4.0_dp*q2
      q68 = 8.0_dp*q9
      q69 = q68*x%d1val1
      q70 = q68*x%d1val3
      q71 = 4.0_dp*pi
      q72 = q45*q71
      q73 = pow2(x%d1val1_d1val3)
      q74 = q29*x%d1val1
      q75 = q12*q17
      q76 = 4.0_dp*q75
      q77 = q12*q18
      q78 = q67*q77
      q79 = 4.0_dp*q1
      q80 = q18*q39
      q81 = q79*q80
      q82 = 0.25_dp*pi
      q83 = q56*x%d1val2
      q84 = q71*x%d1val1_d1val3
      q85 = q5*x%d2val2
      q86 = q12*q16
      q87 = q86*x%d1val3
      q88 = pi*q45
      q89 = q9*x%d1val1_d1val2_d1val3
      q90 = q18*x%d1val1_d1val2
      q91 = q44*q55
      q92 = q1*q80
      q93 = x%d1val3*x%d2val2
      q94 = q16*x%d1val2_d1val3
      q95 = pow2(x%d1val2_d1val3)
      q96 = 8.0_dp*q2
      q97 = 48.0_dp*x%d2val1_d1val3
      q98 = q9*x%d2val1
      q99 = 8.0_dp*q88
      q100 = q75*x%d1val1
      q101 = 24.0_dp*x%d2val1
      q102 = 48.0_dp*q55
      q103 = x%d1val3*x%d2val1
      q104 = q12*q8
      q105 = q77*q96
      q106 = q39*q45
      q107 = 8.0_dp*q106
      q108 = 48.0_dp*q60
      q109 = q18*pow4(pi)
      q110 = q109*q96
      q111 = 0.125_dp*pi
      q112 = q64*q96
      q113 = x%d1val2*x%d2val1
      q114 = q39*q90
      q115 = q75*x%d1val2
      q116 = q79*x%d2val2
      q117 = 24.0_dp*x%d2val2
      q118 = q102*x%d1val2
      unary%val = q1
      unary%d1val1 = -q3*x%d1val1
      unary%d1val2 = -q4
      unary%d1val3 = -q5
      unary%d2val1 = -q11*(q10*q8 + q7)
      unary%d1val1_d1val2 = -q13*q14 - q3*x%d1val1_d1val2
      unary%d1val1_d1val3 = -q15*x%d1val1 - q3*x%d1val1_d1val3
      unary%d2val2 = -q11*(q10*q16 + q6*x%d2val2)
      unary%d1val2_d1val3 = -q15*x%d1val2 - q3*x%d1val2_d1val3
      unary%d2val3 = -pi*(q17 + q18*q9)
      unary%d3val1 = q24*(-q19*x%d3val1 - q21*x%d2val1 + q22*q23)
      unary%d2val1_d1val2 = -q11*(q25*x%d1val1_d1val2 + q26*x%d1val2 - q27*q28 + q6*x%d2val1_d1val2)
      unary%d2val1_d1val3 = -q11*(q25*x%d1val1_d1val3 + q26*x%d1val3 - q29*q8 + q6*x%d2val1_d1val3)
      unary%d1val1_d2val2 = q24*(-q19*x%d1val1_d2val2 - q30*q31 + q34*q36)
      unary%d1val1_d1val2_d1val3 = -q13*q37 - q13*q38 - q15*x%d1val1_d1val2 + q2*q41 - q3*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = pi*(-q42 - q43*q44 + q46*x%d1val1)
      unary%d3val2 = q24*(-q19*x%d3val2 + q23*q48 - q47*x%d2val2)
      unary%d2val2_d1val3 = -q11*(-q16*q29 + q30*q49 + q50*x%d2val2 + q6*x%d2val2_d1val3)
      unary%d1val2_d2val3 = pi*(-q44*q52 + q46*x%d1val2 - q51)
      unary%d3val1_d1val3 = q24*(-q19*x%d3val1_d1val3 - q21*x%d2val1_d1val3 + q22*q61 - q43*q53 + q53*q57 - q54*x%d3val1 + q58*q59)
      unary%d2val1_d1val2_d1val3 = -q11*(-32.0_dp*q14*q55*x%d1val1_d1val3 + q25*x%d1val1_d1val2_d1val3 + q26*x%d1val2_d1val3 - q27*q66 - q28*q32*q40 + q43*q62 + q50*x%d2val1_d1val2 - q57*q62 + q6*x%d2val1_d1val2_d1val3 + q63*x%d2val1_d1val3 - q65*q7)
      unary%d2val1_d2val3 = -q82*(q67*x%d2val1_d2val3 + q68*q73 + q69*x%d1val1_d2val3 + q70*x%d2val1_d1val3 + q72*x%d2val1 - q74*x%d1val1_d1val3 - q76*q8 - q78*x%d2val1 - q8*q81)
      unary%d1val1_d2val2_d1val3 = q24*(-128.0_dp*q30*x%d1val1_d1val2_d1val3 - q19*x%d1val1_d2val2_d1val3 - q31*q52 + q31*q83 + q34*q84 + q36*(16.0_dp*q85 + q32*q87 - q32*x%d2val2_d1val3 + q4*q49) - q54*x%d1val1_d2val2)
      unary%d1val1_d1val2_d2val3 = pi*(-2.0_dp*q43*x%d1val2_d1val3 + q14*q75 + q14*q92 - q2*x%d1val1_d1val2_d2val3 - q20*x%d1val2_d2val3 - q30*x%d1val1_d2val3 + q37*q91 + q38*q91 - q44*q89 + q55*q90 - q88*x%d1val1_d1val2)
      unary%d3val2_d1val3 = q24*(-192.0_dp*q52*x%d2val2 - q19*x%d3val2_d1val3 - q47*x%d2val2_d1val3 + q48*q61 - q54*x%d3val2 + q58*q93*x%d1val2 + q58*q94)
      unary%d2val2_d2val3 = -q82*(-q16*q76 - q16*q81 - q29*x%d1val2*x%d1val2_d1val3 + q67*x%d2val2_d2val3 + q68*q95 + q68*x%d1val2*x%d1val2_d2val3 + q70*x%d2val2_d1val3 + q72*x%d2val2 - q78*x%d2val2)
      unary%d3val1_d2val3 = q111*(-24.0_dp*q20*x%d2val1_d2val3 - 24.0_dp*q98*x%d1val1_d2val3 + 24.0_dp*q104*q42 + q100*q101 + q101*q92*x%d1val1 + q102*q103*x%d1val1_d1val3 + q102*q73*x%d1val1 + q105*x%d3val1 + q107*q22 + q108*q59 - q110*q22 - q43*q97 - q50*x%d3val1_d1val3 + q57*q97 - q96*x%d3val1_d2val3 - q99*x%d3val1)
      unary%d2val1_d1val2_d2val3 = q82*(-4.0_dp*q30*x%d2val1_d2val3 - 4.0_dp*q98*x%d1val2_d2val3 + 4.0_dp*q104*q51 + 4.0_dp*q106*q28 + 8.0_dp*q1*q114*x%d1val1 + 8.0_dp*q100*x%d1val1_d1val2 + 8.0_dp*q12*q14*q42 + 8.0_dp*q60*q66 - q10*x%d1val1_d1val2_d1val3*x%d1val1_d1val3 + q103*q12*q96*x%d1val2_d1val3 - q109*q28*q67 + q112*q73 + q112*x%d1val3*x%d2val1_d1val3 + q113*q76 + q113*q81 + q27*q37*x%d1val1_d1val3 + q29*x%d1val1_d1val2*x%d1val1_d1val3 + q32*q41*x%d1val1_d1val3 - q67*x%d2val1_d1val2_d2val3 - q68*x%d1val1_d1val2*x%d1val1_d2val3 - q68*x%d1val2_d1val3*x%d2val1_d1val3 - q69*x%d1val1_d1val2_d2val3 - q70*x%d2val1_d1val2_d1val3 - q72*x%d2val1_d1val2 + q74*x%d1val1_d1val2_d1val3 + q78*x%d2val1_d1val2)
      unary%d1val1_d2val2_d2val3 = q111*(16.0_dp*q115*x%d1val1_d1val2 + 2.0_dp*pi*x%d1val1_d2val3*(4.0_dp*q33 - q116) + 2.0_dp*q35*(4.0_dp*q45*q86 + 8.0_dp*pi*q51*x%d1val2 + 8.0_dp*q3*q95 + 8.0_dp*q5*x%d2val2_d1val3 + q116*q77 - q16*q67*q80 + q17*q71*x%d2val2 + q32*q65*x%d1val2_d1val3 - q79*x%d2val2_d2val3) + 32.0_dp*q83*x%d1val1_d1val2_d1val3 - q10*x%d1val1_d1val2*x%d1val2_d2val3 + q105*x%d1val1_d2val2 + q114*q32*x%d1val2 + q49*q56*x%d1val1_d1val2 - q49*q89 - q50*x%d1val1_d2val2_d1val3 - q63*x%d1val1_d1val2_d2val3 + q84*(4.0_dp*q85 + 8.0_dp*q4*x%d1val2_d1val3 + q79*q87 - q79*x%d2val2_d1val3) - q96*x%d1val1_d2val2_d2val3 - q99*x%d1val1_d2val2)
      unary%d3val2_d2val3 = q111*(-24.0_dp*q30*x%d2val2_d2val3 - 48.0_dp*q52*x%d2val2_d1val3 + 24.0_dp*q51*q86 + q102*q93*x%d1val2_d1val3 + q105*x%d3val2 + q107*q48 + q108*q94 - q110*q48 + q115*q117 - q117*q9*x%d1val2_d2val3 + q117*q92*x%d1val2 + q118*q95 + q118*x%d1val3*x%d2val2_d1val3 - q50*x%d3val2_d1val3 - q96*x%d3val2_d2val3 - q99*x%d3val2)
   end function cospi_self
   
   function tanpi_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = tan(pi*x%val)
      q1 = pow2(q0)
      q2 = q1 + 1
      q3 = pi*q2
      q4 = q3*x%d1val2
      q5 = q3*x%d1val3
      q6 = pow2(x%d1val1)
      q7 = 2.0_dp*q0
      q8 = pi*q7
      q9 = q6*q8 + x%d2val1
      q10 = pow2(pi)
      q11 = q10*q2
      q12 = q11*x%d1val2
      q13 = q7*x%d1val1
      q14 = q11*x%d1val3
      q15 = pow2(x%d1val2)
      q16 = q15*q8 + x%d2val2
      q17 = q12*q7
      q18 = pow2(x%d1val3)
      q19 = pi*x%d1val1
      q20 = q0*q19
      q21 = 6.0_dp*x%d2val1
      q22 = pow3(x%d1val1)
      q23 = q10*q22
      q24 = 4.0_dp*q1
      q25 = 2.0_dp*q2
      q26 = q20*q21 + q23*q24 + q23*q25 + x%d3val1
      q27 = 4.0_dp*q0
      q28 = q19*q27
      q29 = q10*q25
      q30 = q6*x%d1val2
      q31 = q28*x%d1val1_d1val2 + q29*q30 + x%d2val1_d1val2
      q32 = q14*q7
      q33 = q29*x%d1val3
      q34 = q28*x%d1val1_d1val3 + q33*q6 + x%d2val1_d1val3
      q35 = pi*q27
      q36 = q35*x%d1val2
      q37 = 16.0_dp*q0
      q38 = pi*q1
      q39 = q15*q38
      q40 = q15*q3
      q41 = 16.0_dp*q40 + 32.0_dp*q39 + q37*x%d2val2
      q42 = 0.125_dp*q19
      q43 = q36*x%d1val1_d1val2 + q41*q42 + x%d1val1_d2val2
      q44 = x%d1val1*x%d1val3
      q45 = pow2(q2)
      q46 = pow3(pi)
      q47 = q45*q46
      q48 = 2.0_dp*x%d1val2
      q49 = q47*q48
      q50 = q11*x%d1val2_d1val3
      q51 = q2*q46
      q52 = q51*x%d1val2
      q53 = q24*q52
      q54 = pi*x%d1val1_d1val3
      q55 = q27*q54
      q56 = q0*x%d2val3
      q57 = pi*(2.0_dp*q18*q38 + q18*q3 + q56)
      q58 = 2.0_dp*q57
      q59 = pi*q0
      q60 = q59*x%d1val2
      q61 = 6.0_dp*q60
      q62 = pow3(x%d1val2)
      q63 = q10*q62
      q64 = q24*q63 + q25*q63 + q61*x%d2val2 + x%d3val2
      q65 = q35*x%d1val2_d1val3
      q66 = q20*x%d2val1_d1val3
      q67 = q0*q54
      q68 = q6*x%d1val1_d1val3
      q69 = q1*q10
      q70 = 12.0_dp*q69
      q71 = q21*x%d1val1
      q72 = 6.0_dp*q11
      q73 = q0*x%d1val3
      q74 = 12.0_dp*q22
      q75 = q51*q74
      q76 = q9*x%d1val3
      q77 = 4.0_dp*q11
      q78 = q77*x%d1val3
      q79 = q78*x%d1val1_d1val2
      q80 = q77*x%d1val2
      q81 = x%d1val1*x%d1val1_d1val3
      q82 = q29*q6
      q83 = q27*q51
      q84 = pow2(x%d1val1_d1val3)
      q85 = 8.0_dp*q2
      q86 = q10*q85
      q87 = q44*q86
      q88 = q18*q6
      q89 = 8.0_dp*q20
      q90 = q6*x%d1val3
      q91 = 2.0_dp*x%d2val1_d1val3 + q77*q90 + q89*x%d1val1_d1val3
      q92 = q8*x%d1val3
      q93 = 2.0_dp*x%d2val1 + q35*q6
      q94 = x%d1val2*x%d1val2_d1val3
      q95 = q4*x%d1val2_d1val3
      q96 = q0*q14
      q97 = pi*x%d1val1_d1val2
      q98 = x%d1val1_d1val2_d1val3*x%d1val3
      q99 = x%d1val1_d2val3*x%d1val2
      q100 = q10*q24
      q101 = x%d1val1*x%d1val2
      q102 = 8.0_dp*q69
      q103 = x%d1val1_d1val3*x%d1val2
      q104 = q18*x%d1val1_d1val2
      q105 = q29*x%d2val3
      q106 = q18*q46
      q107 = q106*pow3(q0)
      q108 = q37*x%d1val2
      q109 = q108*q51
      q110 = q18*x%d1val1
      q111 = 6.0_dp*x%d2val2
      q112 = q59*x%d1val2_d1val3
      q113 = q15*x%d1val2_d1val3
      q114 = x%d1val2*x%d1val3
      q115 = q114*x%d2val2
      q116 = 12.0_dp*q62
      q117 = q51*q73
      q118 = pow2(x%d1val2_d1val3)
      q119 = q86*x%d1val3
      q120 = q15*q18
      q121 = 2.0_dp*x%d2val2_d1val3
      q122 = 8.0_dp*q59
      q123 = 2.0_dp*x%d2val2
      q124 = 12.0_dp*q67
      q125 = 24.0_dp*q69
      q126 = q84*x%d1val1
      q127 = q6*x%d1val1_d2val3
      q128 = 12.0_dp*q11
      q129 = q128*q44
      q130 = q128*x%d1val3
      q131 = q18*pow4(pi)
      q132 = q131*q45
      q133 = 72.0_dp*q117
      q134 = 24.0_dp*q22
      q135 = q1*q131
      q136 = q135*q2
      q137 = 16.0_dp*q1
      q138 = 0.5_dp*q57
      q139 = q35*x%d1val1_d1val2
      q140 = 8.0_dp*q67
      q141 = q77*x%d1val1
      q142 = x%d1val1_d1val2*x%d2val3
      q143 = q86*q98
      q144 = q81*q86
      q145 = q119*x%d1val1_d1val2
      q146 = q59*x%d1val2_d2val3
      q147 = 4.0_dp*q30
      q148 = q46*q85
      q149 = q0*q148
      q150 = q104*q149
      q151 = q51*q56
      q152 = q93*x%d2val3
      q153 = q93*x%d1val2_d1val3
      q154 = q93*x%d1val2
      q155 = 4.0_dp*x%d2val1_d1val2 + q19*q37*x%d1val1_d1val2 + q30*q86
      q156 = 0.5_dp*pi
      q157 = q155*q18
      q158 = q6*q77
      q159 = q122*x%d1val2_d1val3
      q160 = q27*x%d2val2
      q161 = 4.0_dp*q40 + 8.0_dp*q39 + q160
      q162 = 8.0_dp*q38
      q163 = q128*q15
      q164 = 4.0_dp*q95 + q0*q121 + q123*q5 + q162*q94 + q163*q73
      q165 = 2.0_dp*q19
      q166 = 12.0_dp*x%d2val2_d1val3
      q167 = q118*x%d1val2
      q168 = q15*x%d1val2_d2val3
      q169 = 12.0_dp*x%d2val2
      q170 = 24.0_dp*q62
      unary%val = q0
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q4
      unary%d1val3 = q5
      unary%d2val1 = q3*q9
      unary%d1val1_d1val2 = q12*q13 + q3*x%d1val1_d1val2
      unary%d1val1_d1val3 = q13*q14 + q3*x%d1val1_d1val3
      unary%d2val2 = q16*q3
      unary%d1val2_d1val3 = q17*x%d1val3 + q3*x%d1val2_d1val3
      unary%d2val3 = q3*(q18*q8 + x%d2val3)
      unary%d3val1 = q26*q3
      unary%d2val1_d1val2 = q17*q9 + q3*q31
      unary%d2val1_d1val3 = q3*q34 + q32*q9
      unary%d1val1_d2val2 = q3*q43
      unary%d1val1_d1val2_d1val3 = q13*q50 + q17*x%d1val1_d1val3 + q3*x%d1val1_d1val2_d1val3 + q32*x%d1val1_d1val2 + q44*q49 + q44*q53
      unary%d1val1_d2val3 = q3*(q55*x%d1val3 + q58*x%d1val1 + x%d1val1_d2val3)
      unary%d3val2 = q3*q64
      unary%d2val2_d1val3 = q16*q32 + q3*(q15*q33 + q36*x%d1val2_d1val3 + x%d2val2_d1val3)
      unary%d1val2_d2val3 = q3*(q58*x%d1val2 + q65*x%d1val3 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q26*q32 + q3*(6.0_dp*q66 + q14*q71 + q21*q67 + q68*q70 + q68*q72 + q73*q75 + x%d3val1_d1val3)
      unary%d2val1_d1val2_d1val3 = q17*q34 + q3*(q28*x%d1val1_d1val2_d1val3 + q30*q83*x%d1val3 + q55*x%d1val1_d1val2 + q79*x%d1val1 + q80*q81 + q82*x%d1val2_d1val3 + x%d2val1_d1val2_d1val3) + q31*q32 + q49*q76 + q50*q7*q9 + q53*q76
      unary%d2val1_d2val3 = q3*(q28*x%d1val1_d2val3 + q35*q84 + q57*q93 + q82*x%d2val3 + q83*q88 + q87*x%d1val1_d1val3 + q91*q92 + x%d2val1_d2val3)
      unary%d1val1_d2val2_d1val3 = q3*(0.125_dp*q41*q54 + q36*x%d1val1_d1val2_d1val3 + q42*(16.0_dp*q5*x%d2val2 + 32.0_dp*q95 + 64.0_dp*q38*q94 + 96.0_dp*q15*q96 + q37*x%d2val2_d1val3) + q65*x%d1val1_d1val2 + q79*x%d1val2 + x%d1val1_d2val2_d1val3) + q32*q43
      unary%d1val1_d1val2_d2val3 = q3*(2.0_dp*q56*q97 + 8.0_dp*q101*q107 + q100*q101*x%d2val3 + q100*q104 + q101*q105 + q102*q103*x%d1val3 + q102*q44*x%d1val2_d1val3 + q103*q78 + q104*q29 + q109*q110 + q19*q7*x%d1val2_d2val3 + q35*q98 + q55*x%d1val2_d1val3 + q78*x%d1val1*x%d1val2_d1val3 + q8*q99 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q3*(q111*q112 + q113*q70 + q113*q72 + q115*q72 + q116*q117 + q61*x%d2val2_d1val3 + x%d3val2_d1val3) + q32*q64
      unary%d2val2_d2val3 = q3*(q105*q15 + q118*q35 + q119*q94 + q120*q83 + q36*x%d1val2_d2val3 + q57*(q123 + q15*q35) + q92*(q121 + q122*q94 + q15*q78) + x%d2val2_d2val3)
      unary%d3val1_d2val3 = q3*(12.0_dp*q0*q110*q51*x%d2val1 + 6.0_dp*q20*x%d2val1_d2val3 + q11*q71*x%d2val3 + q124*x%d2val1_d1val3 + q125*q126 + q126*q128 + q127*q70 + q127*q72 + q129*x%d2val1_d1val3 + q130*x%d1val1_d1val3*x%d2val1 + q132*q74 + q133*q68 + q134*q136 + q138*(24.0_dp*q20*x%d2val1 + 4.0_dp*x%d3val1 + q137*q23 + q23*q85) + q21*q59*x%d1val1_d2val3 + q56*q75 + q92*(12.0_dp*q66 + 2.0_dp*x%d3val1_d1val3 + q117*q134 + q124*x%d2val1 + q125*q68 + q128*q68 + q129*x%d2val1) + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q3*(0.5_dp*q11*q157 + 4.0_dp*q107*q154 + q0*q106*q154*q85 + q100*q114*q91 + q100*q153*x%d1val3 + q109*q44*x%d1val1_d1val3 + q12*q152 + q132*q147 + q135*q30*q85 + q139*x%d1val1_d2val3 + q140*x%d1val1_d1val2_d1val3 + q141*q142 + q141*q99 + q143*x%d1val1 + q144*x%d1val2_d1val3 + q145*x%d1val1_d1val3 + q146*q93 + q147*q151 + q149*q90*x%d1val2_d1val3 + q150*x%d1val1 + q152*q48*q69 + q153*q33 + q155*q156*q56 + q157*q69 + q28*x%d1val1_d1val2_d2val3 + q33*q91*x%d1val2 + q60*(16.0_dp*q14*q81 + 2.0_dp*x%d2val1_d2val3 + q122*q84 + q149*q88 + q158*x%d2val3 + q89*x%d1val1_d2val3) + q8*q91*x%d1val2_d1val3 + q80*q84 + q82*x%d1val2_d2val3 + q92*(2.0_dp*x%d2val1_d1val2_d1val3 + q140*x%d1val1_d1val2 + q144*x%d1val2 + q148*q30*q73 + q158*x%d1val2_d1val3 + q87*x%d1val1_d1val2 + q89*x%d1val1_d1val2_d1val3) + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q3*(2.0_dp*q164*q54 + q138*(4.0_dp*x%d1val1_d2val2 + q108*q97 + q161*q165) + q139*x%d1val2_d2val3 + q142*q80 + q143*x%d1val2 + q145*x%d1val2_d1val3 + q150*x%d1val2 + q156*q161*x%d1val1_d2val3 + q159*x%d1val1_d1val2_d1val3 + q19*(12.0_dp*q120*q47 + 24.0_dp*q1*q120*q51 + 4.0_dp*q118*q3 + 4.0_dp*q4*x%d1val2_d2val3 + 4.0_dp*q5*x%d2val2_d1val3 + 48.0_dp*q94*q96 + q11*q160*q18 + q118*q162 + q123*q3*x%d2val3 + q162*x%d1val2*x%d1val2_d2val3 + q163*q56 + q7*x%d2val2_d2val3) + q36*x%d1val1_d1val2_d2val3 + q92*(2.0_dp*x%d1val1_d2val2_d1val3 + q114*q86*x%d1val1_d1val2 + q122*x%d1val1_d1val2_d1val3*x%d1val2 + q159*x%d1val1_d1val2 + q161*q54 + q164*q165) + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q3*(q0*q169*q18*q52 + q111*q146 + q112*q166 + q113*q133 + q114*q128*x%d2val2_d1val3 + q116*q132 + q116*q151 + q125*q167 + q128*q167 + q130*x%d1val2_d1val3*x%d2val2 + q136*q170 + q138*(24.0_dp*q60*x%d2val2 + 4.0_dp*x%d3val2 + q137*q63 + q63*q85) + q168*q70 + q168*q72 + q61*x%d2val2_d2val3 + q72*x%d1val2*x%d2val2*x%d2val3 + q92*(2.0_dp*x%d3val2_d1val3 + q112*q169 + q113*q125 + q113*q128 + q115*q128 + q117*q170 + q166*q60) + x%d3val2_d2val3)
   end function tanpi_self
   
   function sinh_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = sinh(x%val)
      q1 = cosh(x%val)
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4 = q1*x%d1val3
      q5 = q1*x%d2val1
      q6 = pow2(x%d1val1)
      q7 = q0*q6
      q8 = q1*x%d1val1_d1val2
      q9 = q0*x%d1val1
      q10 = q9*x%d1val2
      q11 = q1*x%d1val1_d1val3
      q12 = q1*x%d2val2
      q13 = pow2(x%d1val2)
      q14 = q0*q13
      q15 = q1*x%d1val2_d1val3
      q16 = q0*x%d1val2
      q17 = q16*x%d1val3
      q18 = q1*x%d2val3
      q19 = pow2(x%d1val3)
      q20 = q0*q19
      q21 = q1*x%d3val1
      q22 = pow3(x%d1val1)
      q23 = q1*q22
      q24 = 3.0_dp*q9
      q25 = q1*x%d2val1_d1val2
      q26 = 2.0_dp*q9
      q27 = q3*q6
      q28 = q0*x%d1val3
      q29 = q1*x%d1val1_d2val2
      q30 = 2.0_dp*q16
      q31 = q0*x%d2val2
      q32 = q1*q13
      q33 = x%d1val2*x%d1val3
      q34 = q1*x%d1val1_d2val3
      q35 = 2.0_dp*q28
      q36 = q0*x%d2val3
      q37 = q1*q19
      q38 = q36 + q37
      q39 = pow3(x%d1val2)
      q40 = 3.0_dp*q16
      q41 = q1*x%d1val2_d2val3
      q42 = q0*x%d2val1
      q43 = 3.0_dp*q42
      q44 = q2*x%d1val3
      q45 = 3.0_dp*x%d2val1
      q46 = 3.0_dp*q6
      q47 = 2.0_dp*q0
      q48 = q47*x%d1val1_d1val2
      q49 = q3*x%d2val1
      q50 = 2.0_dp*q44
      q51 = 2.0_dp*q2
      q52 = q51*x%d1val2
      q53 = pow2(x%d1val1_d1val3)
      q54 = 4.0_dp*x%d1val1_d1val3
      q55 = 2.0_dp*q3
      q56 = q55*x%d1val3
      q57 = q0*x%d2val2_d1val3
      q58 = q3*x%d1val2_d1val3
      q59 = q4*x%d2val2
      q60 = q14*x%d1val3
      q61 = q47*x%d1val2_d1val3
      q62 = x%d1val2*x%d2val3
      q63 = q20*x%d1val2
      q64 = 3.0_dp*q31
      q65 = q3*x%d1val3
      q66 = 3.0_dp*x%d2val2
      q67 = 3.0_dp*q13
      q68 = pow2(x%d1val2_d1val3)
      q69 = 6.0_dp*x%d1val1_d1val3
      q70 = 6.0_dp*q2
      q71 = q4*x%d2val1
      q72 = q20*x%d1val1
      q73 = q7*x%d1val3
      q74 = x%d1val1_d1val2*x%d2val3
      q75 = 4.0_dp*x%d1val1_d1val2_d1val3
      q76 = q4*x%d1val1_d1val2
      q77 = 2.0_dp*x%d1val2_d1val3
      q78 = 2.0_dp*x%d1val1_d1val2
      q79 = 4.0_dp*q0
      q80 = 4.0_dp*x%d2val3
      q81 = 6.0_dp*x%d1val2_d1val3
      q82 = 6.0_dp*q3
      unary%val = q0
      unary%d1val1 = q2
      unary%d1val2 = q3
      unary%d1val3 = q4
      unary%d2val1 = q5 + q7
      unary%d1val1_d1val2 = q10 + q8
      unary%d1val1_d1val3 = q11 + q9*x%d1val3
      unary%d2val2 = q12 + q14
      unary%d1val2_d1val3 = q15 + q17
      unary%d2val3 = q18 + q20
      unary%d3val1 = q21 + q23 + q24*x%d2val1
      unary%d2val1_d1val2 = q16*x%d2val1 + q25 + q26*x%d1val1_d1val2 + q27
      unary%d2val1_d1val3 = q1*x%d2val1_d1val3 + q26*x%d1val1_d1val3 + q28*x%d2val1 + q4*q6
      unary%d1val1_d2val2 = 0.0625_dp*x%d1val1*(16.0_dp*q31 + 16.0_dp*q32) + q29 + q30*x%d1val1_d1val2
      unary%d1val1_d1val2_d1val3 = q1*x%d1val1_d1val2_d1val3 + q16*x%d1val1_d1val3 + q2*q33 + q28*x%d1val1_d1val2 + q9*x%d1val2_d1val3
      unary%d1val1_d2val3 = q34 + q35*x%d1val1_d1val3 + q38*x%d1val1
      unary%d3val2 = q1*q39 + q1*x%d3val2 + q40*x%d2val2
      unary%d2val2_d1val3 = q1*x%d2val2_d1val3 + q13*q4 + q28*x%d2val2 + q30*x%d1val2_d1val3
      unary%d1val2_d2val3 = q35*x%d1val2_d1val3 + q38*x%d1val2 + q41
      unary%d3val1_d1val3 = q1*x%d3val1_d1val3 + q11*q46 + q22*q28 + q24*x%d2val1_d1val3 + q28*x%d3val1 + q43*x%d1val1_d1val3 + q44*q45
      unary%d2val1_d1val2_d1val3 = q1*x%d2val1_d1val2_d1val3 + q15*q6 + q16*x%d2val1_d1val3 + q26*x%d1val1_d1val2_d1val3 + q28*x%d2val1_d1val2 + q33*q7 + q42*x%d1val2_d1val3 + q48*x%d1val1_d1val3 + q49*x%d1val3 + q50*x%d1val1_d1val2 + q52*x%d1val1_d1val3
      unary%d2val1_d2val3 = q1*x%d2val1_d2val3 + q18*q6 + q19*q5 + q19*q7 + q26*x%d1val1_d2val3 + q35*x%d2val1_d1val3 + q36*x%d2val1 + q44*q54 + q47*q53
      unary%d1val1_d2val2_d1val3 = 4.0_dp*x%d1val1*(0.25_dp*q57 + 0.25_dp*q59 + 0.25_dp*q60 + 0.5_dp*q58) + q1*x%d1val1_d2val2_d1val3 + q28*x%d1val1_d2val2 + q30*x%d1val1_d1val2_d1val3 + q48*x%d1val2_d1val3 + q54*(0.25_dp*q31 + 0.25_dp*q32) + q56*x%d1val1_d1val2
      unary%d1val1_d1val2_d2val3 = q1*x%d1val1_d1val2_d2val3 + q16*x%d1val1_d2val3 + q19*q8 + q2*q62 + q35*x%d1val1_d1val2_d1val3 + q36*x%d1val1_d1val2 + q50*x%d1val2_d1val3 + q56*x%d1val1_d1val3 + q61*x%d1val1_d1val3 + q63*x%d1val1 + q9*x%d1val2_d2val3
      unary%d3val2_d1val3 = q1*x%d3val2_d1val3 + q15*q67 + q28*q39 + q28*x%d3val2 + q40*x%d2val2_d1val3 + q64*x%d1val2_d1val3 + q65*q66
      unary%d2val2_d2val3 = 4.0_dp*q58*x%d1val3 + q1*x%d2val2_d2val3 + q12*q19 + q13*q18 + q14*q19 + q30*x%d1val2_d2val3 + q35*x%d2val2_d1val3 + q36*x%d2val2 + q47*q68
      unary%d3val1_d2val3 = q0*q69*x%d2val1_d1val3 + q1*x%d3val1_d2val3 + q19*q21 + q19*q23 + q2*q45*x%d2val3 + q22*q36 + q24*x%d2val1_d2val3 + q34*q46 + q35*x%d3val1_d1val3 + q36*x%d3val1 + q43*x%d1val1_d2val3 + q45*q72 + q53*q70 + q69*q71 + q69*q73 + q70*x%d1val3*x%d2val1_d1val3
      unary%d2val1_d1val2_d2val3 = q0*q54*x%d1val1_d1val2_d1val3 + q1*x%d2val1_d1val2_d2val3 + q10*q54*x%d1val3 + q16*x%d2val1_d2val3 + q19*q25 + q19*q27 + q2*q54*x%d1val2_d1val3 + q26*x%d1val1_d1val2_d2val3 + q35*x%d2val1_d1val2_d1val3 + q36*x%d2val1_d1val2 + q41*q6 + q42*x%d1val2_d2val3 + q44*q75 + q48*x%d1val1_d2val3 + q49*x%d2val3 + q51*q74 + q52*x%d1val1_d2val3 + q53*q55 + q54*q76 + q56*x%d2val1_d1val3 + q61*x%d2val1_d1val3 + q62*q7 + q63*x%d2val1 + q71*q77 + q72*q78 + q73*q77
      unary%d1val1_d2val2_d2val3 = 0.25_dp*x%d1val1*(16.0_dp*q17*x%d1val2_d1val3 + 4.0_dp*q13*q37 + 4.0_dp*q20*x%d2val2 + 8.0_dp*q1*q68 + 8.0_dp*q3*x%d1val2_d2val3 + 8.0_dp*q4*x%d2val2_d1val3 + q12*q80 + q14*q80 + q79*x%d2val2_d2val3) + 0.25_dp*x%d1val1_d2val3*(4.0_dp*q31 + 4.0_dp*q32) + 0.5_dp*x%d1val1_d1val3*(4.0_dp*q57 + 4.0_dp*q59 + 4.0_dp*q60 + 8.0_dp*q58) + 4.0_dp*q76*x%d1val2_d1val3 + q1*x%d1val1_d2val2_d2val3 + q19*q29 + q30*x%d1val1_d1val2_d2val3 + q35*x%d1val1_d2val2_d1val3 + q36*x%d1val1_d2val2 + q48*x%d1val2_d2val3 + q55*q74 + q63*q78 + q65*q75 + q79*x%d1val1_d1val2_d1val3*x%d1val2_d1val3
      unary%d3val2_d2val3 = q1*x%d3val2_d2val3 + q3*q66*x%d2val3 + q35*x%d3val2_d1val3 + q36*q39 + q36*x%d3val2 + q37*q39 + q37*x%d3val2 + q40*x%d2val2_d2val3 + q41*q67 + q57*q81 + q59*q81 + q60*q81 + q63*q66 + q64*x%d1val2_d2val3 + q68*q82 + q82*x%d1val3*x%d2val2_d1val3
   end function sinh_self
   
   function cosh_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = cosh(x%val)
      q1 = sinh(x%val)
      q2 = q1*x%d1val1
      q3 = q1*x%d1val2
      q4 = q1*x%d1val3
      q5 = q1*x%d2val1
      q6 = pow2(x%d1val1)
      q7 = q0*q6
      q8 = q1*x%d1val1_d1val2
      q9 = q0*x%d1val1
      q10 = q9*x%d1val2
      q11 = q1*x%d1val1_d1val3
      q12 = q1*x%d2val2
      q13 = pow2(x%d1val2)
      q14 = q0*q13
      q15 = q1*x%d1val2_d1val3
      q16 = q0*x%d1val2
      q17 = q16*x%d1val3
      q18 = q1*x%d2val3
      q19 = pow2(x%d1val3)
      q20 = q0*q19
      q21 = q1*x%d3val1
      q22 = pow3(x%d1val1)
      q23 = q1*q22
      q24 = 3.0_dp*q9
      q25 = q1*x%d2val1_d1val2
      q26 = 2.0_dp*q9
      q27 = q3*q6
      q28 = q0*x%d1val3
      q29 = q1*x%d1val1_d2val2
      q30 = 2.0_dp*q16
      q31 = q0*x%d2val2
      q32 = q1*q13
      q33 = x%d1val2*x%d1val3
      q34 = q1*x%d1val1_d2val3
      q35 = 2.0_dp*q28
      q36 = q0*x%d2val3
      q37 = q1*q19
      q38 = q36 + q37
      q39 = pow3(x%d1val2)
      q40 = 3.0_dp*q16
      q41 = q1*x%d1val2_d2val3
      q42 = q0*x%d2val1
      q43 = 3.0_dp*q42
      q44 = q2*x%d1val3
      q45 = 3.0_dp*x%d2val1
      q46 = 3.0_dp*q6
      q47 = 2.0_dp*q0
      q48 = q47*x%d1val1_d1val2
      q49 = q3*x%d2val1
      q50 = 2.0_dp*q44
      q51 = 2.0_dp*q2
      q52 = q51*x%d1val2
      q53 = pow2(x%d1val1_d1val3)
      q54 = 4.0_dp*x%d1val1_d1val3
      q55 = 2.0_dp*q3
      q56 = q55*x%d1val3
      q57 = q0*x%d2val2_d1val3
      q58 = q3*x%d1val2_d1val3
      q59 = q4*x%d2val2
      q60 = q14*x%d1val3
      q61 = q47*x%d1val2_d1val3
      q62 = x%d1val2*x%d2val3
      q63 = q20*x%d1val2
      q64 = 3.0_dp*q31
      q65 = q3*x%d1val3
      q66 = 3.0_dp*x%d2val2
      q67 = 3.0_dp*q13
      q68 = pow2(x%d1val2_d1val3)
      q69 = 6.0_dp*x%d1val1_d1val3
      q70 = 6.0_dp*q2
      q71 = q4*x%d2val1
      q72 = q20*x%d1val1
      q73 = q7*x%d1val3
      q74 = x%d1val1_d1val2*x%d2val3
      q75 = 4.0_dp*x%d1val1_d1val2_d1val3
      q76 = q4*x%d1val1_d1val2
      q77 = 2.0_dp*x%d1val2_d1val3
      q78 = 2.0_dp*x%d1val1_d1val2
      q79 = 4.0_dp*q0
      q80 = 4.0_dp*x%d2val3
      q81 = 6.0_dp*x%d1val2_d1val3
      q82 = 6.0_dp*q3
      unary%val = q0
      unary%d1val1 = q2
      unary%d1val2 = q3
      unary%d1val3 = q4
      unary%d2val1 = q5 + q7
      unary%d1val1_d1val2 = q10 + q8
      unary%d1val1_d1val3 = q11 + q9*x%d1val3
      unary%d2val2 = q12 + q14
      unary%d1val2_d1val3 = q15 + q17
      unary%d2val3 = q18 + q20
      unary%d3val1 = q21 + q23 + q24*x%d2val1
      unary%d2val1_d1val2 = q16*x%d2val1 + q25 + q26*x%d1val1_d1val2 + q27
      unary%d2val1_d1val3 = q1*x%d2val1_d1val3 + q26*x%d1val1_d1val3 + q28*x%d2val1 + q4*q6
      unary%d1val1_d2val2 = 0.0625_dp*x%d1val1*(16.0_dp*q31 + 16.0_dp*q32) + q29 + q30*x%d1val1_d1val2
      unary%d1val1_d1val2_d1val3 = q1*x%d1val1_d1val2_d1val3 + q16*x%d1val1_d1val3 + q2*q33 + q28*x%d1val1_d1val2 + q9*x%d1val2_d1val3
      unary%d1val1_d2val3 = q34 + q35*x%d1val1_d1val3 + q38*x%d1val1
      unary%d3val2 = q1*q39 + q1*x%d3val2 + q40*x%d2val2
      unary%d2val2_d1val3 = q1*x%d2val2_d1val3 + q13*q4 + q28*x%d2val2 + q30*x%d1val2_d1val3
      unary%d1val2_d2val3 = q35*x%d1val2_d1val3 + q38*x%d1val2 + q41
      unary%d3val1_d1val3 = q1*x%d3val1_d1val3 + q11*q46 + q22*q28 + q24*x%d2val1_d1val3 + q28*x%d3val1 + q43*x%d1val1_d1val3 + q44*q45
      unary%d2val1_d1val2_d1val3 = q1*x%d2val1_d1val2_d1val3 + q15*q6 + q16*x%d2val1_d1val3 + q26*x%d1val1_d1val2_d1val3 + q28*x%d2val1_d1val2 + q33*q7 + q42*x%d1val2_d1val3 + q48*x%d1val1_d1val3 + q49*x%d1val3 + q50*x%d1val1_d1val2 + q52*x%d1val1_d1val3
      unary%d2val1_d2val3 = q1*x%d2val1_d2val3 + q18*q6 + q19*q5 + q19*q7 + q26*x%d1val1_d2val3 + q35*x%d2val1_d1val3 + q36*x%d2val1 + q44*q54 + q47*q53
      unary%d1val1_d2val2_d1val3 = 4.0_dp*x%d1val1*(0.25_dp*q57 + 0.25_dp*q59 + 0.25_dp*q60 + 0.5_dp*q58) + q1*x%d1val1_d2val2_d1val3 + q28*x%d1val1_d2val2 + q30*x%d1val1_d1val2_d1val3 + q48*x%d1val2_d1val3 + q54*(0.25_dp*q31 + 0.25_dp*q32) + q56*x%d1val1_d1val2
      unary%d1val1_d1val2_d2val3 = q1*x%d1val1_d1val2_d2val3 + q16*x%d1val1_d2val3 + q19*q8 + q2*q62 + q35*x%d1val1_d1val2_d1val3 + q36*x%d1val1_d1val2 + q50*x%d1val2_d1val3 + q56*x%d1val1_d1val3 + q61*x%d1val1_d1val3 + q63*x%d1val1 + q9*x%d1val2_d2val3
      unary%d3val2_d1val3 = q1*x%d3val2_d1val3 + q15*q67 + q28*q39 + q28*x%d3val2 + q40*x%d2val2_d1val3 + q64*x%d1val2_d1val3 + q65*q66
      unary%d2val2_d2val3 = 4.0_dp*q58*x%d1val3 + q1*x%d2val2_d2val3 + q12*q19 + q13*q18 + q14*q19 + q30*x%d1val2_d2val3 + q35*x%d2val2_d1val3 + q36*x%d2val2 + q47*q68
      unary%d3val1_d2val3 = q0*q69*x%d2val1_d1val3 + q1*x%d3val1_d2val3 + q19*q21 + q19*q23 + q2*q45*x%d2val3 + q22*q36 + q24*x%d2val1_d2val3 + q34*q46 + q35*x%d3val1_d1val3 + q36*x%d3val1 + q43*x%d1val1_d2val3 + q45*q72 + q53*q70 + q69*q71 + q69*q73 + q70*x%d1val3*x%d2val1_d1val3
      unary%d2val1_d1val2_d2val3 = q0*q54*x%d1val1_d1val2_d1val3 + q1*x%d2val1_d1val2_d2val3 + q10*q54*x%d1val3 + q16*x%d2val1_d2val3 + q19*q25 + q19*q27 + q2*q54*x%d1val2_d1val3 + q26*x%d1val1_d1val2_d2val3 + q35*x%d2val1_d1val2_d1val3 + q36*x%d2val1_d1val2 + q41*q6 + q42*x%d1val2_d2val3 + q44*q75 + q48*x%d1val1_d2val3 + q49*x%d2val3 + q51*q74 + q52*x%d1val1_d2val3 + q53*q55 + q54*q76 + q56*x%d2val1_d1val3 + q61*x%d2val1_d1val3 + q62*q7 + q63*x%d2val1 + q71*q77 + q72*q78 + q73*q77
      unary%d1val1_d2val2_d2val3 = 0.25_dp*x%d1val1*(16.0_dp*q17*x%d1val2_d1val3 + 4.0_dp*q13*q37 + 4.0_dp*q20*x%d2val2 + 8.0_dp*q1*q68 + 8.0_dp*q3*x%d1val2_d2val3 + 8.0_dp*q4*x%d2val2_d1val3 + q12*q80 + q14*q80 + q79*x%d2val2_d2val3) + 0.25_dp*x%d1val1_d2val3*(4.0_dp*q31 + 4.0_dp*q32) + 0.5_dp*x%d1val1_d1val3*(4.0_dp*q57 + 4.0_dp*q59 + 4.0_dp*q60 + 8.0_dp*q58) + 4.0_dp*q76*x%d1val2_d1val3 + q1*x%d1val1_d2val2_d2val3 + q19*q29 + q30*x%d1val1_d1val2_d2val3 + q35*x%d1val1_d2val2_d1val3 + q36*x%d1val1_d2val2 + q48*x%d1val2_d2val3 + q55*q74 + q63*q78 + q65*q75 + q79*x%d1val1_d1val2_d1val3*x%d1val2_d1val3
      unary%d3val2_d2val3 = q1*x%d3val2_d2val3 + q3*q66*x%d2val3 + q35*x%d3val2_d1val3 + q36*q39 + q36*x%d3val2 + q37*q39 + q37*x%d3val2 + q40*x%d2val2_d2val3 + q41*q67 + q57*q81 + q59*q81 + q60*q81 + q63*q66 + q64*x%d1val2_d2val3 + q68*q82 + q82*x%d1val3*x%d2val2_d1val3
   end function cosh_self
   
   function tanh_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q2 = 1 - q1
      q3 = q2*x%d1val1
      q4 = q2*x%d1val2
      q5 = q2*x%d1val3
      q6 = -0.125_dp + 0.125_dp*q1
      q7 = pow2(x%d1val1)
      q8 = 16.0_dp*q0
      q9 = q7*q8
      q10 = -8.0_dp*x%d2val1 + q9
      q11 = 2.0_dp*q0
      q12 = q11*x%d1val2
      q13 = q11*x%d1val3
      q14 = pow2(x%d1val2)
      q15 = -8.0_dp*x%d2val2 + q14*q8
      q16 = q2*x%d1val2_d1val3
      q17 = q1 - 1
      q18 = pow2(x%d1val3)
      q19 = -0.03125_dp + 0.03125_dp*q1
      q20 = q0*x%d1val1
      q21 = 192.0_dp*x%d2val1
      q22 = pow3(x%d1val1)
      q23 = 128.0_dp*q1
      q24 = 64.0_dp*q17
      q25 = -32.0_dp*x%d3val1 + q20*q21 - q22*q23 - q22*q24
      q26 = 0.25_dp*q0
      q27 = q10*q4
      q28 = 32.0_dp*q20
      q29 = 16.0_dp*q7
      q30 = -8.0_dp*x%d2val1_d1val2 + q28*x%d1val1_d1val2 + q29*q4
      q31 = q26*q5
      q32 = q28*x%d1val1_d1val3
      q33 = -8.0_dp*x%d2val1_d1val3 + q29*q5 + q32
      q34 = q0*x%d1val2
      q35 = 128.0_dp*x%d1val1_d1val2
      q36 = q1*q14
      q37 = q14*q17
      q38 = 16.0_dp*q37 + 32.0_dp*q36 - q8*x%d2val2
      q39 = 4.0_dp*x%d1val1
      q40 = -32.0_dp*x%d1val1_d2val2 + q34*q35 - q38*q39
      q41 = 2.0_dp*q1
      q42 = -2.0_dp + q41
      q43 = q3*x%d1val3
      q44 = q43*x%d1val2
      q45 = 4.0_dp*q1
      q46 = 4.0_dp*q0
      q47 = q46*x%d1val3
      q48 = q0*x%d2val3
      q49 = q18*q41
      q50 = q17*q18
      q51 = -q48 + q49 + q50
      q52 = 2.0_dp*q51
      q53 = 192.0_dp*q34
      q54 = pow3(x%d1val2)
      q55 = -32.0_dp*x%d3val2 - q23*q54 - q24*q54 + q53*x%d2val2
      q56 = 32.0_dp*x%d1val2_d1val3
      q57 = 16.0_dp*q5
      q58 = q0*q5
      q59 = 0.0625_dp*q58
      q60 = q20*x%d2val1_d1val3
      q61 = q0*x%d1val1_d1val3
      q62 = q7*x%d1val1_d1val3
      q63 = q1*q62
      q64 = q17*q62
      q65 = 384.0_dp*q58
      q66 = q0*x%d1val2_d1val3
      q67 = -0.25_dp*q1 + 0.25_dp
      q68 = q10*q67
      q69 = 0.5_dp*q1
      q70 = q4*x%d1val3
      q71 = 32.0_dp*x%d1val1_d1val2
      q72 = x%d1val1_d1val3*x%d1val2
      q73 = -0.5_dp + q69
      q74 = 2.0_dp*x%d2val1 - q46*q7
      q75 = 8.0_dp*q0
      q76 = q75*x%d1val1
      q77 = 4.0_dp*q17
      q78 = q7*q77
      q79 = 2.0_dp*x%d2val1_d1val3 - q76*x%d1val1_d1val3 + q78*x%d1val3
      q80 = pow2(x%d1val1_d1val3)
      q81 = q17*x%d1val1_d1val3
      q82 = x%d1val1*x%d1val3
      q83 = 16.0_dp*q82
      q84 = q50*q75
      q85 = -2.0_dp*x%d2val1_d2val3 + q7*q84 + q75*q80 + q76*x%d1val1_d2val3 - q78*x%d2val3 - q81*q83
      q86 = 4.0_dp*x%d1val1_d1val3
      q87 = x%d1val2*x%d1val2_d1val3
      q88 = q17*x%d1val2
      q89 = q0*x%d1val2_d2val3
      q90 = 4.0_dp*x%d1val1_d1val2
      q91 = q75*x%d1val1_d1val2_d1val3
      q92 = q75*x%d1val1_d1val3
      q93 = x%d1val1_d2val3*x%d1val2
      q94 = 8.0_dp*q1
      q95 = x%d1val1*x%d1val2
      q96 = 16.0_dp*q1
      q97 = q82*x%d1val2_d1val3
      q98 = q18*pow3(q0)
      q99 = q88*x%d2val3
      q100 = 8.0_dp*q17
      q101 = q100*x%d1val2
      q102 = q101*x%d1val1_d1val3
      q103 = q50*x%d1val2
      q104 = q66*x%d2val2
      q105 = q36*x%d1val2_d1val3
      q106 = q37*x%d1val2_d1val3
      q107 = 2.0_dp*x%d2val2_d2val3
      q108 = x%d1val2*x%d1val2_d2val3
      q109 = pow2(x%d1val2_d1val3)
      q110 = 16.0_dp*x%d1val2_d1val3
      q111 = q88*x%d1val3
      q112 = 4.0_dp*q37
      q113 = 2.0_dp*x%d2val2_d1val3
      q114 = 2.0_dp*x%d2val2
      q115 = 12.0_dp*x%d2val1
      q116 = q80*x%d1val1
      q117 = 48.0_dp*q1
      q118 = q7*x%d1val1_d2val3
      q119 = 24.0_dp*q17
      q120 = 12.0_dp*q17
      q121 = q120*x%d2val1
      q122 = x%d1val1*x%d2val3
      q123 = 24.0_dp*x%d2val1
      q124 = q81*x%d1val3
      q125 = q119*q22
      q126 = q123*q20
      q127 = q0*x%d1val3
      q128 = q18*pow2(q17)
      q129 = 24.0_dp*q128
      q130 = q117*q50
      q131 = q75*x%d1val1_d1val2
      q132 = q8*x%d1val1_d1val2_d1val3
      q133 = q100*x%d1val1_d1val2
      q134 = q8*x%d1val1_d1val2
      q135 = q134*x%d1val1
      q136 = q101*q7
      q137 = x%d1val2_d1val3*x%d1val3
      q138 = q137*q17
      q139 = 8.0_dp*x%d1val2
      q140 = q74*x%d1val2
      q141 = q94*x%d1val3
      q142 = q74*x%d1val2_d1val3
      q143 = q77*x%d1val3
      q144 = 4.0_dp*x%d2val1_d1val2 - q135 + q136
      q145 = q79*x%d1val2
      q146 = x%d1val1_d1val2*x%d1val3
      q147 = q134*x%d1val2
      q148 = q46*x%d2val2
      q149 = 8.0_dp*q36 + q112 - q148
      q150 = q114*q17
      q151 = 12.0_dp*q37
      q152 = -q0*q113 - q127*q151 + q150*x%d1val3 + q77*q87 + q87*q94
      q153 = 2.0_dp*x%d1val1
      q154 = 24.0_dp*q36
      q155 = 12.0_dp*q34
      q156 = 12.0_dp*x%d2val2
      q157 = 24.0_dp*q88
      q158 = q119*q54
      q159 = 24.0_dp*q34*x%d2val2
      unary%val = q0
      unary%d1val1 = q3
      unary%d1val2 = q4
      unary%d1val3 = q5
      unary%d2val1 = q10*q6
      unary%d1val1_d1val2 = -q12*q3 + q2*x%d1val1_d1val2
      unary%d1val1_d1val3 = -q13*q3 + q2*x%d1val1_d1val3
      unary%d2val2 = q15*q6
      unary%d1val2_d1val3 = -q13*q4 + q16
      unary%d2val3 = q17*(q11*q18 - x%d2val3)
      unary%d3val1 = q19*q25
      unary%d2val1_d1val2 = q26*q27 + q30*q6
      unary%d2val1_d1val3 = q10*q31 + q33*q6
      unary%d1val1_d2val2 = q19*q40
      unary%d1val1_d1val2_d1val3 = -q11*q5*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3 + q20*q42*x%d1val2_d1val3 + q34*q42*x%d1val1_d1val3 + q42*q44 + q44*q45
      unary%d1val1_d2val3 = q17*(q47*x%d1val1_d1val3 - q52*x%d1val1 - x%d1val1_d2val3)
      unary%d3val2 = q19*q55
      unary%d2val2_d1val3 = q15*q31 + q6*(-8.0_dp*x%d2val2_d1val3 + q14*q57 + q34*q56)
      unary%d1val2_d2val3 = q17*(q47*x%d1val2_d1val3 - q52*x%d1val2 - x%d1val2_d2val3)
      unary%d3val1_d1val3 = q19*(-192.0_dp*q64 - 32.0_dp*x%d3val1_d1val3 - 384.0_dp*q63 + 192.0_dp*q60 + q21*q43 + q21*q61 - q22*q65) + q25*q59
      unary%d2val1_d1val2_d1val3 = -q27*q69*x%d1val3 + q30*q31 + q33*q34*q67 + q6*(-32.0_dp*q0*q7*q70 - 8.0_dp*x%d2val1_d1val2_d1val3 + 32.0_dp*q3*q72 + q16*q29 + q28*x%d1val1_d1val2_d1val3 + q43*q71 + q61*q71) + q66*q68 + q68*q70
      unary%d2val1_d2val3 = q73*(q47*q79 - q52*q74 + q85)
      unary%d1val1_d2val2_d1val3 = q19*(-32.0_dp*x%d1val1_d2val2_d1val3 + 128.0_dp*q34*x%d1val1_d1val2_d1val3 + q35*q66 + q35*q70 - q38*q86 - q39*(64.0_dp*q1*q87 + 96.0_dp*q14*q58 + q56*q88 - q57*x%d2val2 - q8*x%d2val2_d1val3)) + q40*q59
      unary%d1val1_d1val2_d2val3 = q73*(-2.0_dp*x%d1val1_d1val2_d2val3 + 16.0_dp*q95*q98 - q100*q97 - q102*x%d1val3 + q103*q28 - q18*q94*x%d1val1_d1val2 + q39*q89 - q39*q99 + q46*q93 + q48*q90 - q50*q90 - q72*q96*x%d1val3 + q91*x%d1val3 + q92*x%d1val2_d1val3 - q94*q95*x%d2val3 - q96*q97)
      unary%d3val2_d1val3 = q19*(-192.0_dp*q106 - 32.0_dp*x%d3val2_d1val3 - 384.0_dp*q105 + 192.0_dp*q104 + 192.0_dp*q70*x%d2val2 + q53*x%d2val2_d1val3 - q54*q65) + q55*q59
      unary%d2val2_d2val3 = q73*(-q107 + q108*q75 + q109*q75 - q110*q111 - q112*x%d2val3 + q14*q84 + q47*(q112*x%d1val3 + q113 - q75*q87) - q52*(q114 - q14*q46))
      unary%d3val1_d2val3 = q73*(-2.0_dp*x%d3val1_d2val3 - 24.0_dp*q1*q118 + 12.0_dp*q20*x%d2val1_d2val3 + 144.0_dp*q127*q64 + 24.0_dp*q61*x%d2val1_d1val3 + q0*q115*x%d1val1_d2val3 - q116*q117 - q116*q119 - q118*q120 - q119*q82*x%d2val1_d1val3 - q121*q122 - q123*q124 + q125*q48 + q126*q50 - q129*q22 - q130*q22 + q47*(-12.0_dp*q60 + 12.0_dp*q64 + 2.0_dp*x%d3val1_d1val3 + 24.0_dp*q63 - q115*q61 + q121*q82 - q125*q127) - q51*(4.0_dp*x%d3val1 + q100*q22 - q126 + q22*q96))
      unary%d2val1_d1val2_d2val3 = q73*(-16.0_dp*q124*x%d1val1_d1val2 - 2.0_dp*q74*q99 - 2.0_dp*x%d2val1_d1val2_d2val3 - q100*q93*x%d1val1 - q101*q80 - q103*q7*q96 + q11*q74*x%d1val2_d2val3 - q110*q81*x%d1val1 + q111*q32 - q12*q85 - q122*q133 - q128*q139*q7 + q131*x%d1val1_d2val3 + q132*x%d1val1_d1val3 + q135*q50 + q136*q48 + q138*q9 + q139*q74*q98 - q140*q45*x%d2val3 + q140*q50*q8 - q141*q142 - q141*q145 - q142*q143 - q143*q145 + q144*q48 - q144*q49 - q144*q50 - q17*q83*x%d1val1_d1val2_d1val3 + q46*q79*x%d1val2_d1val3 + q47*(2.0_dp*x%d2val1_d1val2_d1val3 + q102*x%d1val1 - q111*q7*q75 + q133*q82 - q76*x%d1val1_d1val2_d1val3 + q78*x%d1val2_d1val3 - q92*x%d1val1_d1val2) + q76*x%d1val1_d1val2_d2val3 - q78*x%d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q73*(-16.0_dp*q111*x%d1val1_d1val2_d1val3 - 2.0_dp*x%d1val1_d2val2_d2val3 - q101*x%d1val1_d1val2*x%d2val3 - q110*q146*q17 + q131*x%d1val2_d2val3 + q132*x%d1val2_d1val3 + q147*q50 - q149*x%d1val1_d2val3 - q152*q86 - q153*(-48.0_dp*q138*q34 + 12.0_dp*q128*q14 - q0*q107 + q108*q77 + q108*q94 + q109*q77 + q109*q94 + q143*x%d2val2_d1val3 - q148*q50 + q150*x%d2val3 - q151*q48 + q154*q50) + q47*(2.0_dp*x%d1val1_d2val2_d1val3 + q101*q146 - q131*x%d1val2_d1val3 + q149*x%d1val1_d1val3 + q152*q153 - q91*x%d1val2) - q51*(4.0_dp*x%d1val1_d2val2 - q147 + q149*q153) + q75*x%d1val1_d1val2_d2val3*x%d1val2)
      unary%d3val2_d2val3 = q73*(-2.0_dp*x%d3val2_d2val3 + 144.0_dp*q37*q66*x%d1val3 + 24.0_dp*q66*x%d2val2_d1val3 - q109*q117*x%d1val2 - q109*q157 - q119*q137*x%d2val2 - q129*q54 - q130*q54 - q151*x%d1val2_d2val3 - q154*x%d1val2_d2val3 + q155*x%d2val2_d2val3 + q156*q89 - q156*q99 - q157*x%d1val3*x%d2val2_d1val3 + q158*q48 + q159*q50 + q47*(-12.0_dp*q104 + 12.0_dp*q106 + 2.0_dp*x%d3val2_d1val3 + 24.0_dp*q105 + q111*q156 - q127*q158 - q155*x%d2val2_d1val3) - q51*(4.0_dp*x%d3val2 + q100*q54 - q159 + q54*q96))
   end function tanh_self
   
   function asin_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q237
      real(dp) :: q236
      real(dp) :: q235
      real(dp) :: q234
      real(dp) :: q233
      real(dp) :: q232
      real(dp) :: q231
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q3 = pow2(x%d1val1)
      q4 = powm1(q1)
      q5 = q4*x%val
      q6 = q3*q5 + x%d2val1
      q7 = powm1(pow3(sqrt(q1)))
      q8 = q7*x%val
      q9 = x%d1val1*x%d1val2
      q10 = q8*x%d1val3
      q11 = pow2(x%d1val2)
      q12 = q11*q5 + x%d2val2
      q13 = pow2(x%d1val3)
      q14 = pow3(x%d1val1)
      q15 = q5*x%d1val1
      q16 = 3.0_dp*x%d2val1
      q17 = q0*q14
      q18 = powm1(pow2(q1))
      q19 = 3.0_dp*q18
      q20 = q14*q4 + q15*q16 + q17*q19 + x%d3val1
      q21 = q8*x%d1val2
      q22 = q3*q4
      q23 = 2.0_dp*q15
      q24 = 2.0_dp*x%d1val2
      q25 = q0*q3
      q26 = q18*q25
      q27 = q22*x%d1val2 + q23*x%d1val1_d1val2 + q24*q26 + x%d2val1_d1val2
      q28 = 2.0_dp*x%d1val3
      q29 = q22*x%d1val3 + q23*x%d1val1_d1val3 + q26*q28 + x%d2val1_d1val3
      q30 = q5*x%d1val1_d1val2
      q31 = 48.0_dp*x%d2val2
      q32 = 48.0_dp*q11
      q33 = 144.0_dp*q0
      q34 = -144.0_dp + q33
      q35 = powm1(q34)
      q36 = q0*q11
      q37 = -20736.0_dp*q35*q36 + q31*x%val + q32
      q38 = q4*x%d1val1
      q39 = 0.020833333333333332_dp*q38
      q40 = q24*q30 + q37*q39 + x%d1val1_d2val2
      q41 = q7*x%d1val3
      q42 = q8*x%d1val2_d1val3
      q43 = q0*x%d1val3
      q44 = 3.0_dp*q43*powm1(pow5(sqrt(q1)))
      q45 = q28*q5
      q46 = 2.0_dp*x%val
      q47 = 4.0_dp*q0
      q48 = -4.0_dp + q47
      q49 = powm1(q48)
      q50 = 24.0_dp*q49
      q51 = q0*q13
      q52 = 2.0_dp*q13 + q46*x%d2val3 - q50*q51
      q53 = 0.5_dp*q52
      q54 = pow3(x%d1val2)
      q55 = q5*x%d1val2
      q56 = 3.0_dp*x%d2val2
      q57 = q0*q54
      q58 = q19*q57 + q4*q54 + q55*q56 + x%d3val2
      q59 = q11*q4
      q60 = x%d1val2*x%d1val2_d1val3
      q61 = 2.0_dp*q5
      q62 = q18*q36
      q63 = x%d1val2_d1val3*x%d1val3
      q64 = q4*x%d1val2
      q65 = 3.0_dp*x%d1val1_d1val3
      q66 = q14*x%d1val3
      q67 = q18*x%val
      q68 = 8.0_dp*q67
      q69 = 6.0_dp*q18
      q70 = x%d1val3*x%d2val1
      q71 = q0*x%d1val1
      q72 = q70*q71
      q73 = pow3(x%val)
      q74 = q66*q73
      q75 = powm1(pow3(q1))
      q76 = 12.0_dp*q75
      q77 = q6*x%d1val2
      q78 = q28*x%d1val1_d1val2
      q79 = q4*x%d1val1_d1val3
      q80 = 2.0_dp*q30
      q81 = x%d1val1*x%d1val1_d1val2
      q82 = q18*q47
      q83 = q82*x%d1val3
      q84 = q9*x%d1val1_d1val3
      q85 = q3*x%d1val2
      q86 = q67*x%d1val3
      q87 = q3*q73
      q88 = q87*x%d1val3
      q89 = q75*x%d1val2
      q90 = x%d1val1*x%d1val1_d1val3
      q91 = 16.0_dp*q49
      q92 = q91*x%d1val3
      q93 = q90*q92
      q94 = 8.0_dp*q49
      q95 = q94*x%val
      q96 = q95*x%d1val1
      q97 = q96*x%d1val1_d2val3
      q98 = q3*q49
      q99 = 4.0_dp*q98
      q100 = q99*x%d2val3
      q101 = pow2(x%d1val1_d1val3)
      q102 = q101*q95
      q103 = powm1(pow2(q48))
      q104 = 128.0_dp*q103
      q105 = q104*q43
      q106 = q105*q90
      q107 = q103*q13
      q108 = q107*x%val
      q109 = 96.0_dp*q108
      q110 = q109*q3
      q111 = 32.0_dp*q103
      q112 = q111*q25
      q113 = q112*x%d2val3
      q114 = powm1(pow3(q48))
      q115 = q114*q13
      q116 = 512.0_dp*q115
      q117 = q116*q87
      q118 = 2.0_dp*x%d2val1 - q3*q95
      q119 = q118*q4
      q120 = 0.25_dp*q52
      q121 = 4.0_dp*x%d1val3
      q122 = q112*x%d1val3 - q121*q98 - q90*q95 + x%d2val1_d1val3
      q123 = x%d1val1_d1val2*x%d1val2
      q124 = x%d2val2_d1val3*x%val
      q125 = 41472.0_dp*q35
      q126 = q0*q60
      q127 = q11*x%d1val3
      q128 = q127*q73
      q129 = q9*x%d2val3
      q130 = q13*q4
      q131 = q28*q64
      q132 = q61*x%d1val2_d1val3
      q133 = q13*q9
      q134 = q0*q19
      q135 = q43*x%d1val2
      q136 = q135*q69
      q137 = 3.0_dp*x%d1val2_d1val3
      q138 = q54*x%d1val3
      q139 = q138*q73
      q140 = x%d1val2*x%d1val2_d2val3
      q141 = 4.0_dp*q11
      q142 = q141*q49
      q143 = pow2(x%d1val2_d1val3)
      q144 = q104*x%d1val3
      q145 = q111*q36
      q146 = q11*q73
      q147 = 2.0_dp*x%d2val2
      q148 = q120*q4
      q149 = q50*x%d1val1
      q150 = x%d1val3*x%d2val1_d1val3
      q151 = 12.0_dp*q49
      q152 = q151*x%d1val1
      q153 = x%d2val1*x%d2val3
      q154 = q50*x%d1val1_d1val3
      q155 = x%d2val1_d1val3*x%val
      q156 = x%d2val1*x%val
      q157 = q151*q156
      q158 = 12.0_dp*q98
      q159 = q14*x%d2val3
      q160 = q104*x%val
      q161 = 128.0_dp*q107
      q162 = 192.0_dp*q103
      q163 = x%d1val1*x%d2val1
      q164 = 96.0_dp*q103
      q165 = 768.0_dp*q103
      q166 = x%d1val3*x%val
      q167 = q166*x%d1val1_d1val3
      q168 = q0*q162
      q169 = 288.0_dp*q103
      q170 = q103*q33
      q171 = q170*q3
      q172 = 768.0_dp*q114
      q173 = q172*q73
      q174 = q115*q73
      q175 = 1536.0_dp*q174
      q176 = q114*q88
      q177 = pow4(x%val)
      q178 = q13*q177*powm1(pow4(q48))
      q179 = 18432.0_dp*q178
      q180 = q81*q94
      q181 = q92*x%d1val1_d1val2_d1val3
      q182 = q90*x%d1val2_d1val3
      q183 = q9*q94
      q184 = x%d1val1_d1val2*x%d1val1_d1val3
      q185 = q95*x%d1val1_d1val2
      q186 = q91*x%val
      q187 = q186*x%d1val1_d1val2_d1val3
      q188 = 4.0_dp*x%d1val2_d2val3
      q189 = q94*x%d1val2
      q190 = 192.0_dp*q108
      q191 = 64.0_dp*q103
      q192 = q191*q71
      q193 = x%d1val1_d1val2*x%d2val3
      q194 = q0*q104
      q195 = q0*q191
      q196 = q164*q85
      q197 = x%d2val3*x%val
      q198 = q63*x%val
      q199 = q191*x%d1val2
      q200 = q0*q199
      q201 = q73*x%d1val3
      q202 = 1024.0_dp*q174
      q203 = 512.0_dp*x%d1val2
      q204 = q114*q87
      q205 = q114*q51
      q206 = 0.5_dp*q118
      q207 = q118*x%d1val2
      q208 = 1.5_dp*q18
      q209 = q0*x%d2val3
      q210 = 2.0_dp*x%d2val1_d1val2 - q186*q81 - q189*q3 + q199*q25
      q211 = 0.5_dp*q210
      q212 = x%d1val1_d1val2*x%d1val3
      q213 = q189*x%d1val1_d1val2
      q214 = q63*x%d1val1_d1val2
      q215 = q95*x%d1val2
      q216 = x%d2val2*x%val
      q217 = q32*q49
      q218 = 4.0_dp*q216 - q0*q217 + q141
      q219 = q218*q49
      q220 = q218*x%d1val1
      q221 = 8.0_dp*q220
      q222 = q103*q221
      q223 = 48.0_dp*q49
      q224 = q217*x%val
      q225 = 4.0_dp*q60 - q126*q223 + q128*q162 + q147*x%d1val3 - q224*x%d1val3 + q46*x%d2val2_d1val3
      q226 = q225*q49
      q227 = 2.0_dp*x%d1val1
      q228 = q0*q143
      q229 = q50*x%d1val2
      q230 = x%d1val3*x%d2val2_d1val3
      q231 = q151*x%d1val2
      q232 = q231*x%d2val2
      q233 = q63*x%d2val2
      q234 = q11*x%d1val2_d2val3
      q235 = q54*x%d2val3
      q236 = q164*x%d2val2
      q237 = q151*x%d1val2_d1val3
      unary%val = asin(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q2*q6
      unary%d1val1_d1val2 = q2*x%d1val1_d1val2 + q8*q9
      unary%d1val1_d1val3 = q10*x%d1val1 + q2*x%d1val1_d1val3
      unary%d2val2 = q12*q2
      unary%d1val2_d1val3 = q10*x%d1val2 + q2*x%d1val2_d1val3
      unary%d2val3 = q2*(q13*q5 + x%d2val3)
      unary%d3val1 = q2*q20
      unary%d2val1_d1val2 = q2*q27 + q21*q6
      unary%d2val1_d1val3 = q10*q6 + q2*q29
      unary%d1val1_d2val2 = q2*q40
      unary%d1val1_d1val2_d1val3 = q10*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3 + q21*x%d1val1_d1val3 + q41*q9 + q42*x%d1val1 + q44*q9
      unary%d1val1_d2val3 = q2*(q38*q53 + q45*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q2*q58
      unary%d2val2_d1val3 = q10*q12 + q2*(q28*q62 + q59*x%d1val3 + q60*q61 + x%d2val2_d1val3)
      unary%d1val2_d2val3 = q2*(q53*q64 + q61*q63 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q10*q20 + q2*(3.0_dp*q15*x%d2val1_d1val3 + 9.0_dp*q26*x%d1val1_d1val3 + q16*q38*x%d1val3 + q22*q65 + q5*q65*x%d2val1 + q66*q68 + q69*q72 + q74*q76 + x%d3val1_d1val3)
      unary%d2val1_d1val2_d1val3 = q10*q27 + q2*(2.0_dp*q26*x%d1val2_d1val3 + 2.0_dp*q79*q9 + 6.0_dp*q85*q86 + 8.0_dp*q88*q89 + q22*x%d1val2_d1val3 + q23*x%d1val1_d1val2_d1val3 + q38*q78 + q80*x%d1val1_d1val3 + q81*q83 + q82*q84 + x%d2val1_d1val2_d1val3) + q21*q29 + q41*q77 + q42*q6 + q44*q77
      unary%d2val1_d2val3 = q2*(-q100 - q102 + q106 + q110 + q113 - q117 + q119*q120 + q122*q45 - q93 - q97 + x%d2val1_d2val3)
      unary%d1val1_d2val2_d1val3 = q10*q40 + q2*(0.020833333333333332_dp*q37*q79 + 0.041666666666666664_dp*q37*q86*x%d1val1 + 2.0_dp*q55*x%d1val1_d1val2_d1val3 + q123*q83 + q39*(48.0_dp*q124 + 5971968.0_dp*q128*powm1(pow2(q34)) + 96.0_dp*q60 - q125*q126 - q125*q127*x%val + q31*x%d1val3) + q64*q78 + q80*x%d1val2_d1val3 + x%d1val1_d2val2_d1val3)
      unary%d1val1_d1val2_d2val3 = q2*(15.0_dp*q133*q73*q75 + 2.0_dp*q38*q63 + 9.0_dp*q133*q67 + q129*q134 + q129*q4 + q130*x%d1val1_d1val2 + q131*x%d1val1_d1val3 + q132*x%d1val1_d1val3 + q136*x%d1val1_d1val3 + q15*x%d1val2_d2val3 + q19*q51*x%d1val1_d1val2 + q30*x%d2val3 + q45*x%d1val1_d1val2_d1val3 + q55*x%d1val1_d2val3 + q63*q69*q71 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q10*q58 + q2*(3.0_dp*q55*x%d2val2_d1val3 + 9.0_dp*q62*x%d1val2_d1val3 + q136*x%d2val2 + q137*q5*x%d2val2 + q137*q59 + q138*q68 + q139*q76 + q56*q64*x%d1val3 + x%d3val2_d1val3)
      unary%d2val2_d2val3 = q2*(q109*q11 - q116*q146 + q126*q144 - q140*q95 - q142*x%d2val3 - q143*q95 + q145*x%d2val3 + q148*(-q11*q95 + q147) + q45*(-q142*x%d1val3 + q145*x%d1val3 - q60*q95 + x%d2val2_d1val3) - q60*q92 + x%d2val2_d2val3)
      unary%d3val1_d2val3 = q2*(-4352.0_dp*q115*q17 - 4608.0_dp*q176*x%d1val1_d1val3 + 288.0_dp*q108*q163 - q101*q149 + q101*q169*q71 + q14*q161 + q14*q179 + q148*(2.0_dp*x%d3val1 - q14*q94 - q149*q156 + q164*q17) - q149*q150 + q150*q162*q71 - q152*q153 - q152*x%d2val1_d2val3*x%val + q153*q164*q71 - q154*q155 - q154*q70 - q157*x%d1val1_d2val3 - q158*x%d1val1_d2val3 + q159*q160 - q159*q173 - q163*q175 + q165*q167*q3 + q168*q70*x%d1val1_d1val3 + q171*x%d1val1_d2val3 + q45*(-q152*q155 - q152*q70 - q157*x%d1val1_d1val3 - q158*x%d1val1_d1val3 + q160*q66 + q164*q72 + q171*x%d1val1_d1val3 - q172*q74 + x%d3val1_d1val3) + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(-1024.0_dp*q204*q63 - 2048.0_dp*q114*q201*q84 - 3072.0_dp*q205*q85 + 12288.0_dp*q178*q85 + 384.0_dp*q103*q166*q84 + 4.5_dp*q13*q207*q67 + 7.5_dp*q118*q13*q73*q89 + 96.0_dp*q107*q85 - q101*q189 + q101*q200 + q105*q184 + q112*x%d1val2_d2val3 + q118*q134*q63 + q119*q63 + q122*q131 + q122*q132 + q122*q136 + q130*q211 + q144*q71*x%d1val1_d1val2_d1val3 + q162*q198*q3 - q180*x%d2val3 - q181*x%d1val1 + q182*q194 - q182*q91 - q183*x%d1val1_d2val3 - q184*q92 - q185*x%d1val1_d2val3 - q187*x%d1val1_d1val3 - q188*q98 + q190*q81 + q192*q193 + q195*q9*x%d1val1_d2val3 + q196*q197 - q202*q81 - q203*q204*x%d2val3 + q206*q5*x%d1val2_d2val3 + q206*q64*x%d2val3 + q207*q208*q209 + q208*q210*q51 + q211*q5*x%d2val3 + q45*(q112*x%d1val2_d1val3 + q166*q196 - q176*q203 - q180*x%d1val3 - q183*x%d1val1_d1val3 - q184*q95 + q192*q212 + q195*q84 - q96*x%d1val1_d1val2_d1val3 - q99*x%d1val2_d1val3 + x%d2val1_d1val2_d1val3) - q55*(q100 + q102 - q106 - q110 - q113 + q117 + q93 + q97 - x%d2val1_d2val3) - q96*x%d1val1_d1val2_d2val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(-128.0_dp*q205*q220 - 4.0_dp*q226*x%d1val1_d1val3 + 16.0_dp*q103*q167*q218 + q105*x%d1val1_d1val2_d1val3*x%d1val2 + q107*q221 + q111*q166*q225*x%d1val1 + q123*q190 - q123*q202 + q148*(2.0_dp*x%d1val1_d2val2 - q123*q186 - q219*q227) - q181*x%d1val2 - q185*x%d1val2_d2val3 - q187*x%d1val2_d1val3 + q193*q200 + q194*q214 + q197*q222 - q213*x%d2val3 - q214*q91 - q215*x%d1val1_d1val2_d2val3 - q219*x%d1val1_d2val3 - q227*q49*(-192.0_dp*q166*q49*q60 - 3072.0_dp*q11*q115*q177 + 4.0_dp*q143 + 960.0_dp*q107*q36 - q0*q140*q223 + q121*x%d2val2_d1val3 - q13*q217 + q146*q162*x%d2val3 + q147*x%d2val3 + q165*q201*q60 + q188*x%d1val2 - q223*q228 - q224*x%d2val3 + q46*x%d2val2_d2val3) - q45*(-q166*q222 + q185*x%d1val2_d1val3 - q200*q212 + q213*x%d1val3 + q215*x%d1val1_d1val2_d1val3 + q219*x%d1val1_d1val3 + q226*q227 - x%d1val1_d2val2_d1val3) + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(-4352.0_dp*q205*q54 - 4608.0_dp*q114*q146*q63 + 288.0_dp*q107*q216*x%d1val2 + q11*q165*q198 - q124*q50*x%d1val2_d1val3 - q143*q229 + q148*(2.0_dp*x%d3val2 + q164*q57 - q216*q229 - q54*q94) - q151*q216*x%d1val2_d2val3 - q151*q234 + q160*q235 + q161*q54 + q168*q230*x%d1val2 + q168*q233 + q169*q228*x%d1val2 + q170*q234 - q173*q235 - q175*x%d1val2*x%d2val2 + q179*q54 + q209*q236*x%d1val2 - q229*q230 - q231*x%d2val2_d2val3*x%val - q232*x%d2val3 - q233*q50 + q45*(q11*q170*x%d1val2_d1val3 - q11*q237 - q124*q231 + q135*q236 + q138*q160 - q139*q172 - q216*q237 - q232*x%d1val3 + x%d3val2_d1val3) + x%d3val2_d2val3)
   end function asin_self
   
   function acos_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q3 = pow2(x%d1val1)
      q4 = powm1(q1)
      q5 = q4*x%val
      q6 = q3*q5 + x%d2val1
      q7 = powm1(pow3(sqrt(q1)))
      q8 = q7*x%val
      q9 = x%d1val1*x%d1val2
      q10 = q8*x%d1val3
      q11 = pow2(x%d1val2)
      q12 = q11*q5 + x%d2val2
      q13 = pow2(x%d1val3)
      q14 = pow3(x%d1val1)
      q15 = q5*x%d1val1
      q16 = 3.0_dp*x%d2val1
      q17 = q0*q14
      q18 = powm1(pow2(q1))
      q19 = 3.0_dp*q18
      q20 = q14*q4 + q15*q16 + q17*q19 + x%d3val1
      q21 = q8*x%d1val2
      q22 = q3*q4
      q23 = 2.0_dp*q15
      q24 = 2.0_dp*x%d1val2
      q25 = q0*q3
      q26 = q18*q25
      q27 = q22*x%d1val2 + q23*x%d1val1_d1val2 + q24*q26 + x%d2val1_d1val2
      q28 = 2.0_dp*x%d1val3
      q29 = q22*x%d1val3 + q23*x%d1val1_d1val3 + q26*q28 + x%d2val1_d1val3
      q30 = q5*x%d1val1_d1val2
      q31 = 48.0_dp*x%d2val2
      q32 = 48.0_dp*q11
      q33 = 144.0_dp*q0
      q34 = -144.0_dp + q33
      q35 = powm1(q34)
      q36 = q0*q11
      q37 = -20736.0_dp*q35*q36 + q31*x%val + q32
      q38 = q4*x%d1val1
      q39 = 0.020833333333333332_dp*q38
      q40 = q24*q30 + q37*q39 + x%d1val1_d2val2
      q41 = q7*x%d1val3
      q42 = q8*x%d1val2_d1val3
      q43 = q0*x%d1val3
      q44 = 3.0_dp*q43*powm1(pow5(sqrt(q1)))
      q45 = q28*q5
      q46 = 2.0_dp*x%val
      q47 = 4.0_dp*q0
      q48 = -4.0_dp + q47
      q49 = powm1(q48)
      q50 = 24.0_dp*q49
      q51 = q0*q13
      q52 = 2.0_dp*q13 + q46*x%d2val3 - q50*q51
      q53 = 0.5_dp*q52
      q54 = pow3(x%d1val2)
      q55 = q5*x%d1val2
      q56 = 3.0_dp*x%d2val2
      q57 = q0*q54
      q58 = q19*q57 + q4*q54 + q55*q56 + x%d3val2
      q59 = q11*q4
      q60 = x%d1val2*x%d1val2_d1val3
      q61 = 2.0_dp*q5
      q62 = q18*q36
      q63 = x%d1val2_d1val3*x%d1val3
      q64 = q4*x%d1val2
      q65 = 3.0_dp*x%d1val1_d1val3
      q66 = q14*x%d1val3
      q67 = q18*x%val
      q68 = 8.0_dp*q67
      q69 = 6.0_dp*q18
      q70 = q0*x%d1val1
      q71 = x%d1val3*x%d2val1
      q72 = q70*q71
      q73 = pow3(x%val)
      q74 = q66*q73
      q75 = powm1(pow3(q1))
      q76 = 12.0_dp*q75
      q77 = q6*x%d1val2
      q78 = q28*x%d1val1_d1val2
      q79 = q4*x%d1val1_d1val3
      q80 = 2.0_dp*q30
      q81 = x%d1val1*x%d1val1_d1val2
      q82 = q18*q47
      q83 = q82*x%d1val3
      q84 = q9*x%d1val1_d1val3
      q85 = q3*x%d1val2
      q86 = q85*x%d1val3
      q87 = q73*q75
      q88 = 8.0_dp*q49
      q89 = q88*x%val
      q90 = 2.0_dp*x%d2val1 - q3*q89
      q91 = q4*q90
      q92 = 0.25_dp*q52
      q93 = q89*x%d1val1
      q94 = 4.0_dp*x%d1val3
      q95 = q3*q49
      q96 = powm1(pow2(q48))
      q97 = 32.0_dp*q96
      q98 = q25*q97
      q99 = -q93*x%d1val1_d1val3 - q94*q95 + q98*x%d1val3 + x%d2val1_d1val3
      q100 = q49*x%d1val1_d1val3
      q101 = 16.0_dp*x%d1val3
      q102 = q101*x%d1val1
      q103 = 4.0_dp*q95
      q104 = pow2(x%d1val1_d1val3)
      q105 = 128.0_dp*q96
      q106 = q105*x%d1val3
      q107 = q106*q70
      q108 = q3*x%val
      q109 = q13*q96
      q110 = 96.0_dp*q109
      q111 = q13*q3
      q112 = powm1(pow3(q48))
      q113 = q112*q73
      q114 = 512.0_dp*q113
      q115 = q100*q102 + q103*x%d2val3 + q104*q89 - q107*x%d1val1_d1val3 - q108*q110 + q111*q114 + q93*x%d1val1_d2val3 - q98*x%d2val3 - x%d2val1_d2val3
      q116 = x%d1val1_d1val2*x%d1val2
      q117 = x%d1val1*x%d1val3
      q118 = x%d2val2_d1val3*x%val
      q119 = 41472.0_dp*q35
      q120 = q0*q60
      q121 = q11*x%d1val3
      q122 = q121*q73
      q123 = q9*x%d2val3
      q124 = q13*q4
      q125 = q28*q64
      q126 = q61*x%d1val2_d1val3
      q127 = q13*q9
      q128 = q0*q19
      q129 = q43*x%d1val1_d1val3
      q130 = q69*x%d1val2
      q131 = 3.0_dp*x%d1val2_d1val3
      q132 = q54*x%d1val3
      q133 = q130*q43
      q134 = q101*q49
      q135 = x%d1val2*x%d1val2_d2val3
      q136 = 4.0_dp*q11
      q137 = q136*q49
      q138 = pow2(x%d1val2_d1val3)
      q139 = q11*x%val
      q140 = q36*q97
      q141 = q11*q13
      q142 = 2.0_dp*x%d2val2
      q143 = q4*q92
      q144 = q50*x%d1val1
      q145 = x%d1val3*x%d2val1_d1val3
      q146 = 12.0_dp*q49
      q147 = q146*x%d1val1
      q148 = x%d2val1*x%d2val3
      q149 = q147*x%val
      q150 = q50*x%d1val1_d1val3
      q151 = x%d2val1*x%val
      q152 = 12.0_dp*q95
      q153 = q14*x%d2val3
      q154 = q105*x%val
      q155 = 128.0_dp*q109
      q156 = 192.0_dp*q96
      q157 = q109*x%d1val1
      q158 = 96.0_dp*q96
      q159 = x%d1val1_d1val3*x%d1val3
      q160 = 768.0_dp*q96
      q161 = q0*q156
      q162 = 288.0_dp*q96
      q163 = q33*q96
      q164 = q163*q3
      q165 = 768.0_dp*q113
      q166 = q113*q13
      q167 = 1536.0_dp*q166
      q168 = q113*q3
      q169 = 4352.0_dp*q112
      q170 = pow4(x%val)
      q171 = q170*powm1(pow4(q48))
      q172 = 18432.0_dp*q13*q171
      q173 = q81*q88
      q174 = q88*q9
      q175 = q89*x%d1val1_d1val2
      q176 = 16.0_dp*x%val
      q177 = 4.0_dp*x%d1val2_d2val3
      q178 = q88*x%d1val2
      q179 = q96*x%val
      q180 = q84*x%d1val3
      q181 = 192.0_dp*x%val
      q182 = q109*q181
      q183 = 64.0_dp*q96
      q184 = q183*q70
      q185 = x%d1val1_d1val2*x%d2val3
      q186 = q0*q183
      q187 = q158*x%d1val2
      q188 = q108*q187
      q189 = q183*x%d1val2
      q190 = q0*q189
      q191 = 1024.0_dp*q166
      q192 = 3072.0_dp*q112
      q193 = 0.5_dp*q90
      q194 = q90*x%d1val2
      q195 = q13*q194
      q196 = 1.5_dp*q18
      q197 = q0*x%d2val3
      q198 = q176*q49
      q199 = 2.0_dp*x%d2val1_d1val2 - q178*q3 + q189*q25 - q198*q81
      q200 = 0.5_dp*q199
      q201 = x%d1val1_d1val2*x%d1val3
      q202 = q178*x%d1val1_d1val2
      q203 = q63*x%d1val1_d1val2
      q204 = x%d1val1_d1val2_d1val3*x%d1val2
      q205 = q89*x%d1val2
      q206 = x%d2val2*x%val
      q207 = q32*q49
      q208 = 4.0_dp*q206 - q0*q207 + q136
      q209 = q208*q49
      q210 = 8.0_dp*q208
      q211 = q210*q96
      q212 = 48.0_dp*q49
      q213 = q207*x%val
      q214 = 4.0_dp*q60 - q120*q212 + q122*q156 + q142*x%d1val3 - q213*x%d1val3 + q46*x%d2val2_d1val3
      q215 = q117*x%val
      q216 = 2.0_dp*x%d1val1
      q217 = q60*x%d1val3
      q218 = q0*q138
      q219 = q216*q49
      q220 = q50*x%d1val2
      q221 = x%d1val3*x%d2val2_d1val3
      q222 = q146*x%d1val2
      q223 = q222*x%d2val2
      q224 = q63*x%d2val2
      q225 = q11*x%d1val2_d2val3
      q226 = q54*x%d2val3
      q227 = q187*x%d2val2
      q228 = q146*x%d1val2_d1val3
      unary%val = acos(x%val)
      unary%d1val1 = -q2*x%d1val1
      unary%d1val2 = -q2*x%d1val2
      unary%d1val3 = -q2*x%d1val3
      unary%d2val1 = -q2*q6
      unary%d1val1_d1val2 = -q2*x%d1val1_d1val2 - q8*q9
      unary%d1val1_d1val3 = -q10*x%d1val1 - q2*x%d1val1_d1val3
      unary%d2val2 = -q12*q2
      unary%d1val2_d1val3 = -q10*x%d1val2 - q2*x%d1val2_d1val3
      unary%d2val3 = -q2*(q13*q5 + x%d2val3)
      unary%d3val1 = -q2*q20
      unary%d2val1_d1val2 = -q2*q27 - q21*q6
      unary%d2val1_d1val3 = -q10*q6 - q2*q29
      unary%d1val1_d2val2 = -q2*q40
      unary%d1val1_d1val2_d1val3 = -q10*x%d1val1_d1val2 - q2*x%d1val1_d1val2_d1val3 - q21*x%d1val1_d1val3 - q41*q9 - q42*x%d1val1 - q44*q9
      unary%d1val1_d2val3 = -q2*(q38*q53 + q45*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = -q2*q58
      unary%d2val2_d1val3 = -q10*q12 - q2*(q28*q62 + q59*x%d1val3 + q60*q61 + x%d2val2_d1val3)
      unary%d1val2_d2val3 = -q2*(q53*q64 + q61*q63 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q10*q20 - q2*(3.0_dp*q15*x%d2val1_d1val3 + 9.0_dp*q26*x%d1val1_d1val3 + q16*q38*x%d1val3 + q22*q65 + q5*q65*x%d2val1 + q66*q68 + q69*q72 + q74*q76 + x%d3val1_d1val3)
      unary%d2val1_d1val2_d1val3 = -q10*q27 - q2*(2.0_dp*q26*x%d1val2_d1val3 + 2.0_dp*q79*q9 + 6.0_dp*q67*q86 + 8.0_dp*q86*q87 + q22*x%d1val2_d1val3 + q23*x%d1val1_d1val2_d1val3 + q38*q78 + q80*x%d1val1_d1val3 + q81*q83 + q82*q84 + x%d2val1_d1val2_d1val3) - q21*q29 - q41*q77 - q42*q6 - q44*q77
      unary%d2val1_d2val3 = q2*(q115 - q45*q99 - q91*q92)
      unary%d1val1_d2val2_d1val3 = -q10*q40 - q2*(0.020833333333333332_dp*q37*q79 + 0.041666666666666664_dp*q117*q37*q67 + 2.0_dp*q55*x%d1val1_d1val2_d1val3 + q116*q83 + q39*(48.0_dp*q118 + 5971968.0_dp*q122*powm1(pow2(q34)) + 96.0_dp*q60 - q119*q120 - q119*q121*x%val + q31*x%d1val3) + q64*q78 + q80*x%d1val2_d1val3 + x%d1val1_d2val2_d1val3)
      unary%d1val1_d1val2_d2val3 = -q2*(15.0_dp*q127*q87 + 2.0_dp*q38*q63 + 9.0_dp*q127*q67 + q123*q128 + q123*q4 + q124*x%d1val1_d1val2 + q125*x%d1val1_d1val3 + q126*x%d1val1_d1val3 + q129*q130 + q15*x%d1val2_d2val3 + q19*q51*x%d1val1_d1val2 + q30*x%d2val3 + q45*x%d1val1_d1val2_d1val3 + q55*x%d1val1_d2val3 + q63*q69*q70 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = -q10*q58 - q2*(3.0_dp*q55*x%d2val2_d1val3 + 9.0_dp*q62*x%d1val2_d1val3 + q131*q5*x%d2val2 + q131*q59 + q132*q68 + q132*q73*q76 + q133*x%d2val2 + q56*q64*x%d1val3 + x%d3val2_d1val3)
      unary%d2val2_d2val3 = q2*(-q106*q120 - q110*q139 + q114*q141 + q134*q60 + q135*q89 + q137*x%d2val3 + q138*q89 - q140*x%d2val3 - q143*(-q11*q89 + q142) - q45*(-q137*x%d1val3 + q140*x%d1val3 - q60*q89 + x%d2val2_d1val3) - x%d2val2_d2val3)
      unary%d3val1_d2val3 = q2*(-288.0_dp*q151*q157 + 4608.0_dp*q159*q168 + q104*q144 - q104*q162*q70 - q108*q159*q160 + q13*q169*q17 - q14*q155 - q14*q172 - q143*(2.0_dp*x%d3val1 - q14*q88 - q144*q151 + q158*q17) + q144*q145 - q145*q156*q70 + q146*q151*x%d1val1_d2val3 + q147*q148 - q148*q158*q70 + q149*x%d2val1_d2val3 + q150*q71 + q150*x%d2val1_d1val3*x%val + q152*x%d1val1_d2val3 - q153*q154 + q153*q165 - q161*q71*x%d1val1_d1val3 - q164*x%d1val1_d2val3 + q167*x%d1val1*x%d2val1 - q45*(-12.0_dp*q100*q151 - 768.0_dp*q112*q74 - q147*q71 - q149*x%d2val1_d1val3 - q152*x%d1val1_d1val3 + q154*q66 + q158*q72 + q164*x%d1val1_d1val3 + x%d3val1_d1val3) - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(-12288.0_dp*q111*q171*x%d1val2 - 384.0_dp*q179*q180 - 4.5_dp*q195*q67 - 7.5_dp*q195*q87 + 1024.0_dp*q168*q63 + 16.0_dp*q100*x%d1val1*x%d1val2_d1val3 + 2048.0_dp*q113*q180 + q100*q101*x%d1val1_d1val2 + q100*q176*x%d1val1_d1val2_d1val3 + q102*q49*x%d1val1_d1val2_d1val3 + q104*q178 - q104*q190 - q105*q129*x%d1val1_d1val2 - q105*q70*x%d1val1_d1val3*x%d1val2_d1val3 - q107*x%d1val1_d1val2_d1val3 - q108*q156*q63 - q110*q85 + q114*q85*x%d2val3 + q115*q55 - q124*q200 - q125*q99 - q126*q99 - q128*q63*q90 - q133*q99 + q173*x%d2val3 + q174*x%d1val1_d2val3 + q175*x%d1val1_d2val3 + q177*q95 - q182*q81 - q184*q185 - q186*q9*x%d1val1_d2val3 - q188*x%d2val3 + q191*q81 + q192*q51*q85 - q193*q5*x%d1val2_d2val3 - q193*q64*x%d2val3 - q194*q196*q197 - q196*q199*q51 - q200*q5*x%d2val3 - q45*(-q103*x%d1val2_d1val3 - q114*q86 - q173*x%d1val3 - q174*x%d1val1_d1val3 - q175*x%d1val1_d1val3 + q184*q201 + q186*q84 + q188*x%d1val3 - q93*x%d1val1_d1val2_d1val3 + q98*x%d1val2_d1val3 + x%d2val1_d1val2_d1val3) - q63*q91 + q93*x%d1val1_d1val2_d2val3 - q98*x%d1val2_d2val3 - x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(128.0_dp*q112*q208*q51*x%d1val1 + 16.0_dp*q203*q49 + 4.0_dp*q100*q214 - q0*q105*q203 - q101*q179*q208*x%d1val1_d1val3 - q105*q204*q43 - q116*q182 + q116*q191 + q134*q204 - q143*(2.0_dp*x%d1val1_d2val2 - q116*q198 - q209*q216) - q157*q210 + q175*x%d1val2_d2val3 - q185*q190 + q198*x%d1val1_d1val2_d1val3*x%d1val2_d1val3 + q202*x%d2val3 + q205*x%d1val1_d1val2_d2val3 + q209*x%d1val1_d2val3 - q211*x%d1val1*x%d2val3*x%val - q214*q215*q97 + q219*(4.0_dp*q138 + 960.0_dp*q109*q36 - q0*q135*q212 + q11*q156*q73*x%d2val3 - q13*q207 - q141*q170*q192 + q142*x%d2val3 + q160*q217*q73 + q177*x%d1val2 - q181*q217*q49 - q212*q218 - q213*x%d2val3 + q46*x%d2val2_d2val3 + q94*x%d2val2_d1val3) + q45*(q175*x%d1val2_d1val3 - q190*q201 + q202*x%d1val3 + q205*x%d1val1_d1val2_d1val3 + q209*x%d1val1_d1val3 - q211*q215 + q214*q219 - x%d1val1_d2val2_d1val3) - x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(-288.0_dp*q109*q206*x%d1val2 + 4608.0_dp*q11*q113*q63 + q118*q50*x%d1val2_d1val3 + q138*q220 - q139*q160*q63 - q143*(2.0_dp*x%d3val2 + q158*q57 - q206*q220 - q54*q88) + q146*q206*x%d1val2_d2val3 + q146*q225 - q154*q226 - q155*q54 - q161*q221*x%d1val2 - q161*q224 - q162*q218*x%d1val2 - q163*q225 + q165*q226 + q167*x%d1val2*x%d2val2 + q169*q51*q54 - q172*q54 - q197*q227 + q220*q221 + q222*x%d2val2_d2val3*x%val + q223*x%d2val3 + q224*q50 - q45*(q11*q163*x%d1val2_d1val3 - q11*q228 - q118*q222 + q132*q154 - q132*q165 - q206*q228 - q223*x%d1val3 + q227*q43 + x%d3val2_d1val3) - x%d3val2_d2val3)
   end function acos_self
   
   function atan_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q232
      real(dp) :: q231
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q3 = 2.0_dp*x%d2val1
      q4 = pow2(x%d1val1)
      q5 = 144.0_dp*q0 + 144.0_dp
      q6 = powm1(q5)
      q7 = 576.0_dp*q6
      q8 = q7*x%val
      q9 = q3 - q4*q8
      q10 = 72.0_dp*q6
      q11 = powm1(pow2(q1))
      q12 = 2.0_dp*x%val
      q13 = q11*q12
      q14 = x%d1val1*x%d1val2
      q15 = q13*x%d1val3
      q16 = 2.0_dp*x%d2val2
      q17 = pow2(x%d1val2)
      q18 = q16 - q17*q8
      q19 = pow2(x%d1val3)
      q20 = 4.0_dp*q0 + 4.0_dp
      q21 = powm1(q20)
      q22 = q19*q21
      q23 = 8.0_dp*q22
      q24 = 4.0_dp*q21
      q25 = 2.0_dp*x%d3val1
      q26 = pow3(x%d1val1)
      q27 = q6*x%d1val1
      q28 = q27*x%val
      q29 = 1728.0_dp*x%d2val1
      q30 = q0*q26
      q31 = powm1(pow2(q5))
      q32 = 331776.0_dp*q31
      q33 = q25 - q26*q7 - q28*q29 + q30*q32
      q34 = q31*x%d1val2
      q35 = 20736.0_dp*x%val
      q36 = q35*q9
      q37 = 2.0_dp*x%d2val1_d1val2
      q38 = 1152.0_dp*q28
      q39 = q4*q7
      q40 = q0*q4
      q41 = 165888.0_dp*q40
      q42 = q34*q41 + q37 - q38*x%d1val1_d1val2 - q39*x%d1val2
      q43 = q31*x%d1val3
      q44 = 2.0_dp*x%d2val1_d1val3 - q38*x%d1val1_d1val3 - q39*x%d1val3 + q41*q43
      q45 = 2.0_dp*x%d1val1_d2val2
      q46 = q6*x%d1val2
      q47 = 1152.0_dp*x%d1val1_d1val2
      q48 = q47*x%val
      q49 = 48.0_dp*x%d2val2
      q50 = q49*x%val
      q51 = q17*q6
      q52 = -27648.0_dp*q0*q51 + 48.0_dp*q17 + q50
      q53 = 12.0_dp*q27
      q54 = q45 - q46*q48 - q52*q53
      q55 = q14*x%d1val3
      q56 = x%d1val1*x%d1val2_d1val3
      q57 = x%d1val1_d1val3*x%d1val2
      q58 = 16.0_dp*q21
      q59 = q58*x%val
      q60 = q59*x%d1val3
      q61 = -32.0_dp*q0*q22 + 2.0_dp*q19 + q12*x%d2val3
      q62 = q24*q61
      q63 = 2.0_dp*x%d3val2
      q64 = pow3(x%d1val2)
      q65 = x%d2val2*x%val
      q66 = 1728.0_dp*q46
      q67 = q0*q64
      q68 = q32*q67 + q63 - q64*q7 - q65*q66
      q69 = q35*q43
      q70 = 2.0_dp*x%d2val2_d1val3
      q71 = x%d1val2*x%d1val2_d1val3
      q72 = q6*q71
      q73 = 1152.0_dp*x%val
      q74 = q17*x%d1val3
      q75 = q0*q31
      q76 = q27*x%d1val3
      q77 = q6*x%d1val1_d1val3
      q78 = 1728.0_dp*q77
      q79 = x%d2val1*x%val
      q80 = q26*x%val
      q81 = 829440.0_dp*q43
      q82 = q0*x%d1val1
      q83 = q4*x%d1val1_d1val3
      q84 = 995328.0_dp*q75
      q85 = pow3(x%val)
      q86 = q26*x%d1val3
      q87 = q85*q86
      q88 = powm1(pow3(q5))
      q89 = 191102976.0_dp*q88
      q90 = q34*x%d1val3
      q91 = q31*x%d1val2_d1val3
      q92 = q88*x%d1val3
      q93 = q0*x%d1val2
      q94 = x%d1val1*x%d1val1_d1val2
      q95 = q0*q32
      q96 = q95*x%d1val3
      q97 = q14*x%d1val1_d1val3
      q98 = q4*x%val
      q99 = 497664.0_dp*q90
      q100 = q4*x%d1val2
      q101 = q3 - q4*q59
      q102 = 2.0_dp*q21
      q103 = q102*q61
      q104 = q59*x%d1val1
      q105 = 8.0_dp*q21
      q106 = q105*q4
      q107 = powm1(pow2(q20))
      q108 = 64.0_dp*q107
      q109 = q108*q40
      q110 = -q104*x%d1val1_d1val3 - q106*x%d1val3 + q109*x%d1val3 + x%d2val1_d1val3
      q111 = x%d1val1_d1val3*x%d1val3
      q112 = 32.0_dp*q21
      q113 = q112*x%d1val1
      q114 = q4*x%d2val3
      q115 = pow2(x%d1val1_d1val3)
      q116 = 256.0_dp*q107
      q117 = q116*q82
      q118 = q107*q19
      q119 = 192.0_dp*q118
      q120 = q0*q108
      q121 = q19*q4
      q122 = powm1(pow3(q20))
      q123 = q122*q85
      q124 = 1024.0_dp*q123
      q125 = q104*x%d1val1_d2val3 + q105*q114 + q111*q113 - q111*q117 - q114*q120 + q115*q59 - q119*q98 + q121*q124 - x%d2val1_d2val3
      q126 = x%d1val1_d1val2*x%d1val2
      q127 = q49*x%d1val3
      q128 = x%d2val2_d1val3*x%val
      q129 = 48.0_dp*q128
      q130 = x%d1val3*x%val
      q131 = q74*q85
      q132 = q105*x%d2val3
      q133 = q58*x%d1val3
      q134 = q105*x%val
      q135 = x%d1val1_d1val2*x%val
      q136 = q58*x%d1val2
      q137 = q59*x%d1val2_d1val3
      q138 = q134*x%d1val2
      q139 = 384.0_dp*q118
      q140 = q0*q107
      q141 = q140*x%d2val3
      q142 = 128.0_dp*q14
      q143 = q0*q116
      q144 = q143*x%d1val3
      q145 = q0*q118
      q146 = 128.0_dp*x%d1val1_d1val2
      q147 = q123*q19
      q148 = 3072.0_dp*q147
      q149 = x%d1val3*x%d2val2
      q150 = 1728.0_dp*x%d1val2_d1val3
      q151 = q64*x%val
      q152 = q17*x%d1val2_d1val3
      q153 = q64*x%d1val3
      q154 = q112*x%d1val3
      q155 = x%d1val2*x%d1val2_d2val3
      q156 = pow2(x%d1val2_d1val3)
      q157 = q119*x%val
      q158 = 64.0_dp*q17
      q159 = q158*x%d2val3
      q160 = q17*q19
      q161 = 48.0_dp*q21
      q162 = q161*x%d1val1
      q163 = x%d1val3*x%d2val1_d1val3
      q164 = 24.0_dp*q21
      q165 = q164*x%d1val1
      q166 = q165*x%d2val1
      q167 = q165*x%val
      q168 = q111*x%d2val1
      q169 = x%d1val1_d1val3*x%val
      q170 = q164*q79
      q171 = q115*x%d1val1
      q172 = q4*x%d1val1_d2val3
      q173 = 320.0_dp*q107
      q174 = q173*x%d2val3
      q175 = 320.0_dp*q118
      q176 = 384.0_dp*q140
      q177 = q118*x%d1val1
      q178 = x%d1val1*x%d2val1
      q179 = q107*q130
      q180 = 768.0_dp*q140
      q181 = 2048.0_dp*q123
      q182 = q181*x%d2val3
      q183 = 12288.0_dp*q123
      q184 = q122*q19
      q185 = 11264.0_dp*q184
      q186 = pow4(x%val)
      q187 = q186*powm1(pow4(q20))
      q188 = 49152.0_dp*q187*q19
      q189 = 192.0_dp*x%d1val3
      q190 = x%d1val1_d1val2_d1val3*x%d1val3
      q191 = q56*x%d1val1_d1val3
      q192 = q136*x%d1val1
      q193 = q111*x%d1val1_d1val2
      q194 = q59*x%d1val1_d1val2
      q195 = q112*x%d1val1_d1val2_d1val3
      q196 = q146*x%d1val1
      q197 = q107*x%d1val2
      q198 = x%d1val2_d1val3*x%d1val3
      q199 = 128.0_dp*x%d1val2
      q200 = q181*q19
      q201 = q101*q24
      q202 = x%d1val2*x%d2val3
      q203 = q101*q198
      q204 = q101*x%d1val2
      q205 = 128.0_dp*q140
      q206 = q136*x%d1val3
      q207 = q107*q199*q40 - q113*q135 - q136*q4 + q37
      q208 = q207*q24
      q209 = x%d2val3*x%val
      q210 = q116*q93
      q211 = q140*x%d1val3
      q212 = x%d1val1_d1val2_d1val3*x%d1val2
      q213 = x%d1val2_d1val3*x%val
      q214 = q135*x%d1val2
      q215 = q0*q21
      q216 = 4.0_dp*q17 + 4.0_dp*q65 - q158*q215
      q217 = q102*q216
      q218 = 16.0_dp*q216
      q219 = q218*x%d1val1
      q220 = 64.0_dp*q215
      q221 = q21*x%val
      q222 = -64.0_dp*q221*q74 + 4.0_dp*q71 + q116*q131 + q16*x%d1val3 - q220*q71 + q70*x%val
      q223 = q24*x%d1val1
      q224 = x%d1val3*x%d2val2_d1val3
      q225 = q224*x%d1val2
      q226 = q202*x%d2val2
      q227 = q164*x%d1val2
      q228 = q21*x%d1val2_d1val3
      q229 = q164*q65
      q230 = q156*x%d1val2
      q231 = q164*q17
      q232 = 192.0_dp*q140
      unary%val = atan(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q10*q9
      unary%d1val1_d1val2 = -q13*q14 + q2*x%d1val1_d1val2
      unary%d1val1_d1val3 = -q15*x%d1val1 + q2*x%d1val1_d1val3
      unary%d2val2 = q10*q18
      unary%d1val2_d1val3 = -q15*x%d1val2 + q2*x%d1val2_d1val3
      unary%d2val3 = q24*(-q23*x%val + x%d2val3)
      unary%d3val1 = q10*q33
      unary%d2val1_d1val2 = q10*q42 - q34*q36
      unary%d2val1_d1val3 = q10*q44 - q36*q43
      unary%d1val1_d2val2 = q10*q54
      unary%d1val1_d1val2_d1val3 = -2.0_dp*q11*q55 + 8.0_dp*q0*q55*powm1(pow3(q1)) - q13*q56 - q13*q57 - q15*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q24*(-q60*x%d1val1_d1val3 - q62*x%d1val1 + x%d1val1_d2val3)
      unary%d3val2 = q10*q68
      unary%d2val2_d1val3 = q10*(165888.0_dp*q74*q75 - q7*q74 + q70 - q72*q73) - q18*q69
      unary%d1val2_d2val3 = q24*(-q60*x%d1val2_d1val3 - q62*x%d1val2 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q10*(-1728.0_dp*q28*x%d2val1_d1val3 + 2.0_dp*x%d3val1_d1val3 + 497664.0_dp*q43*q82*x%d2val1 - q29*q76 - q4*q78 - q78*q79 + q80*q81 + q83*q84 - q87*q89) - q33*q69
      unary%d2val1_d1val2_d1val3 = -20736.0_dp*q9*q90 + 11943936.0_dp*q9*q92*q93 + q10*(-1152.0_dp*q27*q57 - 95551488.0_dp*q100*q85*q92 + 2.0_dp*x%d2val1_d1val2_d1val3 - q38*x%d1val1_d1val2_d1val3 - q39*x%d1val2_d1val3 + q41*q91 - q47*q76 - q48*q77 + q94*q96 + q95*q97 + q98*q99) - q34*q35*q44 - q36*q91 - q42*q69
      unary%d2val1_d2val3 = -q24*(q101*q103 + q110*q60 + q125)
      unary%d1val1_d2val2_d1val3 = q10*(-12.0_dp*q52*q77 + 2.0_dp*x%d1val1_d2val2_d1val3 + 3456.0_dp*q43*q52*x%d1val1*x%val + q126*q96 - q46*q47*x%d1val3 - q46*q73*x%d1val1_d1val2_d1val3 - q48*q6*x%d1val2_d1val3 - q53*(-55296.0_dp*q0*q72 - 55296.0_dp*q130*q51 + 7962624.0_dp*q131*q31 + 96.0_dp*q71 + q127 + q129)) - q54*q69
      unary%d1val1_d1val2_d2val3 = q24*(-q111*q136 - q132*q135 - q132*q14 - q133*q56 - q134*x%d1val1*x%d1val2_d2val3 - q137*x%d1val1_d1val3 - q138*x%d1val1_d2val3 + q139*q14*x%val - q14*q148 + q141*q142 + q144*q56 + q144*q57 + q145*q146 - q23*x%d1val1_d1val2 - q60*x%d1val1_d1val2_d1val3 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q10*(2.0_dp*x%d3val2_d1val3 + q0*q99*x%d2val2 - q128*q66 - q149*q66 - q150*q51 - q150*q6*q65 + q151*q81 + q152*q84 - q153*q85*q89) - q68*q69
      unary%d2val2_d2val3 = -q24*(q103*(q16 - q17*q59) + q124*q160 + q132*q17 - q140*q159 - q144*q71 + q154*q71 + q155*q59 + q156*q59 - q157*q17 + q60*(-q105*q74 + q120*q74 - q59*q71 + x%d2val2_d1val3) - x%d2val2_d2val3)
      unary%d3val1_d2val3 = -q24*(-192.0_dp*q141*q178 - 1920.0_dp*q179*q83 - 576.0_dp*q177*q79 + q103*(q116*q30 - q162*q79 + q25 - q26*q58) + q148*q178 + q161*q168 + q161*q169*x%d2val1_d1val3 + q161*q171 + q162*q163 - q163*q176*x%d1val1 + q164*q172 + q166*x%d2val3 + q167*x%d2val1_d2val3 - q168*q176 + q170*x%d1val1_d2val3 - q171*q180 - q172*q176 - q174*q80 - q175*q26 + q182*q26 + q183*q83*x%d1val3 + q185*q30 - q188*q26 + q60*(-2048.0_dp*q122*q87 + q140*q178*q189 - q164*q83 - q166*x%d1val3 - q167*x%d2val1_d1val3 - q170*x%d1val1_d1val3 + q173*q86*x%val + q176*q83 + x%d3val1_d1val3) - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q24*(-1536.0_dp*q147*q204 - 4096.0_dp*q123*q55*x%d1val1_d1val3 - 6144.0_dp*q184*q40*x%d1val2 + 192.0_dp*q114*q197*x%val + 24576.0_dp*q121*q187*x%d1val2 + 384.0_dp*q107*q198*q98 + 384.0_dp*q135*q177 + 64.0_dp*q145*q207 + 768.0_dp*q107*q169*q55 + q100*q119 + q101*q120*q202 - q104*x%d1val1_d1val2_d2val3 - q105*q203 - q106*x%d1val2_d2val3 + q109*x%d1val2_d2val3 - q110*q137 - q110*q206 + q110*q210*x%d1val3 - q112*q191 - q112*q193 - q113*q190 - q114*q124*x%d1val2 - q115*q136 + q115*q140*q199 + q117*q190 + q125*q138 + q140*q142*x%d1val1_d2val3 + q141*q196 + q143*q191 + q143*q193 + q157*q204 - q169*q195 - q181*q198*q4 - q19*q208 - q192*x%d1val1_d2val3 - q194*x%d1val1_d2val3 - q200*q94 - q201*q202 - q201*x%d1val2_d2val3*x%val + q203*q205 - q208*q209 - q58*q94*x%d2val3 - q60*(-q100*q124*x%d1val3 - q104*x%d1val1_d1val2_d1val3 - q106*x%d1val2_d1val3 + q109*x%d1val2_d1val3 - q133*q94 + q189*q197*q98 - q192*x%d1val1_d1val3 - q194*x%d1val1_d1val3 + q196*q211 + q205*q97 + x%d2val1_d1val2_d1val3) + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q24*(-256.0_dp*q184*q216*q82 + 32.0_dp*q107*q111*q216*x%val - q103*(-q112*q214 - q216*q223 + q45) - q105*q222*x%d1val1_d1val3 + q107*q209*q219 + q108*q130*q222*x%d1val1 - q126*q200 - q136*x%d1val1_d1val2*x%d2val3 + q139*q214 + q140*q146*q202 + q143*q198*x%d1val1_d1val2 - q154*q212 - q154*x%d1val1_d1val2*x%d1val2_d1val3 + q177*q218 + q190*q210 - q194*x%d1val2_d2val3 - q195*q213 - q217*x%d1val1_d2val3 - q223*(-256.0_dp*q130*q21*q71 - 4096.0_dp*q122*q160*q186 + 1024.0_dp*q107*q71*q85*x%d1val3 + 1280.0_dp*q145*q17 + 4.0_dp*q155 + 4.0_dp*q156 + 4.0_dp*q224 + q116*q17*q85*x%d2val3 + q12*x%d2val2_d2val3 - q155*q220 - q156*q220 - q158*q22 - q159*q221 + q16*x%d2val3) - q59*x%d1val1_d1val2_d2val3*x%d1val2 + q60*(q137*x%d1val1_d1val2 - q146*q211*x%d1val2 - q179*q219 + q206*x%d1val1_d1val2 + q212*q59 + q217*x%d1val1_d1val3 + q222*q223 - x%d1val1_d2val2_d1val3) + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = -q24*(-1920.0_dp*q107*q213*q74 - 576.0_dp*q118*q65*x%d1val2 + q103*(q116*q67 - q21*q50*x%d1val2 - q58*q64 + q63) + q127*q228 + q129*q228 + q148*x%d1val2*x%d2val2 - q149*q176*x%d1val2_d1val3 - q151*q174 + q161*q225 + q161*q230 + q164*q226 - q17*q176*x%d1val2_d2val3 - q175*q64 - q176*q225 - q180*q230 + q182*q64 + q183*q74*x%d1val2_d1val3 + q185*q67 - q188*q64 - q226*q232 + q227*x%d2val2_d2val3*x%val + q229*x%d1val2_d2val3 + q231*x%d1val2_d2val3 + q60*(-q128*q227 + q130*q173*q64 - q149*q227 + q149*q232*x%d1val2 + q152*q176 - q153*q181 - q229*x%d1val2_d1val3 - q231*x%d1val2_d1val3 + x%d3val2_d1val3) - x%d3val2_d2val3)
   end function atan_self
   
   function asinpi_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q233
      real(dp) :: q232
      real(dp) :: q231
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q5 = 144.0_dp*q1
      q6 = -144.0_dp + q5
      q7 = powm1(q6)
      q8 = 144.0_dp*q7
      q9 = q8*x%val
      q10 = -q4*q9 + x%d2val1
      q11 = x%d1val1*x%val
      q12 = powm1(pow3(sqrt(q2)))
      q13 = q0*x%d1val2
      q14 = q12*q13
      q15 = q0*q12
      q16 = q11*x%d1val3
      q17 = pow2(x%d1val2)
      q18 = -q17*q9 + x%d2val2
      q19 = q14*x%val
      q20 = 4.0_dp*x%val
      q21 = pow2(x%d1val3)
      q22 = 4.0_dp*q1
      q23 = -4.0_dp + q22
      q24 = powm1(q23)
      q25 = q21*q24
      q26 = pow3(x%d1val1)
      q27 = 432.0_dp*q7
      q28 = q27*x%d2val1
      q29 = q1*q26
      q30 = powm1(pow2(q6))
      q31 = 62208.0_dp*q30
      q32 = -q11*q28 - q26*q8 + q29*q31 + x%d3val1
      q33 = 288.0_dp*q7
      q34 = q11*x%d1val1_d1val2
      q35 = q4*q8
      q36 = q4*x%d1val2
      q37 = q1*q30
      q38 = 41472.0_dp*q37
      q39 = -q33*q34 - q35*x%d1val2 + q36*q38 + x%d2val1_d1val2
      q40 = x%d1val3*x%val
      q41 = q15*q40
      q42 = q11*q33
      q43 = q38*q4
      q44 = -q35*x%d1val3 - q42*x%d1val1_d1val3 + q43*x%d1val3 + x%d2val1_d1val3
      q45 = 2.0_dp*x%val
      q46 = powm1(q2)
      q47 = q46*x%d1val1_d1val2
      q48 = q47*x%d1val2
      q49 = 48.0_dp*x%d2val2
      q50 = 48.0_dp*q17
      q51 = q17*q7
      q52 = -20736.0_dp*q1*q51 + q49*x%val + q50
      q53 = q46*x%d1val1
      q54 = 0.020833333333333332_dp*q53
      q55 = q45*q48 + q52*q54 + x%d1val1_d2val2
      q56 = x%d1val1*x%d1val3
      q57 = q15*x%d1val2_d1val3
      q58 = 3.0_dp*q1*q13*powm1(pow5(sqrt(q2)))
      q59 = q46*x%d1val1_d1val3
      q60 = q45*x%d1val3
      q61 = -24.0_dp*q1*q25 + 2.0_dp*q21 + q45*x%d2val3
      q62 = 0.5_dp*q61
      q63 = pow3(x%d1val2)
      q64 = q27*x%d1val2
      q65 = x%d2val2*x%val
      q66 = q1*q63
      q67 = q31*q66 - q63*q8 - q64*q65 + x%d3val2
      q68 = x%d1val2*x%d1val2_d1val3
      q69 = q33*x%val
      q70 = q17*x%d1val3
      q71 = x%d1val2_d1val3*x%d1val3
      q72 = q45*q46
      q73 = q46*x%d1val2
      q74 = x%d1val1_d1val3*x%val
      q75 = q4*x%d1val1_d1val3
      q76 = q26*q40
      q77 = 165888.0_dp*q30
      q78 = 124416.0_dp*q37
      q79 = q56*x%d2val1
      q80 = 186624.0_dp*q37
      q81 = pow3(x%val)
      q82 = q26*q81
      q83 = q82*x%d1val3
      q84 = powm1(pow3(q6))
      q85 = 35831808.0_dp*q84
      q86 = q10*x%d1val3
      q87 = q56*x%d1val1_d1val2
      q88 = x%d1val1*x%d1val1_d1val3
      q89 = q88*x%d1val2
      q90 = x%d1val1_d1val2*x%d1val1_d1val3
      q91 = 82944.0_dp*q37
      q92 = q36*q40
      q93 = q81*x%d1val3
      q94 = q36*q93
      q95 = x%d1val1_d1val3*x%d1val3
      q96 = q24*x%d1val1
      q97 = 16.0_dp*q96
      q98 = q95*q97
      q99 = 8.0_dp*q24
      q100 = q99*x%val
      q101 = q100*x%d1val1
      q102 = q101*x%d1val1_d2val3
      q103 = q24*q4
      q104 = 4.0_dp*q103
      q105 = q104*x%d2val3
      q106 = pow2(x%d1val1_d1val3)
      q107 = q100*q106
      q108 = powm1(pow2(q23))
      q109 = 128.0_dp*q108
      q110 = q1*q109
      q111 = q56*x%d1val1_d1val3
      q112 = q110*q111
      q113 = q108*q21
      q114 = 96.0_dp*q113
      q115 = q4*x%val
      q116 = q114*q115
      q117 = q108*q4
      q118 = 32.0_dp*q117
      q119 = q1*x%d2val3
      q120 = q118*q119
      q121 = powm1(pow3(q23))
      q122 = q121*q21
      q123 = q122*q81
      q124 = 512.0_dp*q123
      q125 = q124*q4
      q126 = 2.0_dp*x%d2val1 - q100*q4
      q127 = q126*q46
      q128 = 0.25_dp*q61
      q129 = 4.0_dp*x%d1val3
      q130 = q1*x%d1val3
      q131 = -q100*q88 - q103*q129 + q118*q130 + x%d2val1_d1val3
      q132 = q46*q60
      q133 = 2.0_dp*x%d1val3
      q134 = q45*x%d1val2_d1val3
      q135 = x%d1val1_d1val2*x%d1val3
      q136 = powm1(pow2(q2))
      q137 = q136*x%d1val2
      q138 = x%d2val2_d1val3*x%val
      q139 = q70*q81
      q140 = x%d1val2*x%d2val3
      q141 = x%d1val2_d2val3*x%val
      q142 = x%d2val3*x%val
      q143 = q73*x%val
      q144 = q137*q21
      q145 = q1*q136
      q146 = 3.0_dp*q145
      q147 = 6.0_dp*q145
      q148 = q81*x%d1val2
      q149 = q148*q21*powm1(pow3(q2))
      q150 = x%d1val3*x%d2val2
      q151 = q65*x%d1val2_d1val3
      q152 = q40*q63
      q153 = q17*x%d1val2_d1val3
      q154 = q63*q81
      q155 = q154*x%d1val3
      q156 = 16.0_dp*x%d1val3
      q157 = q24*q68
      q158 = x%d1val2*x%d1val2_d2val3
      q159 = 4.0_dp*q17
      q160 = q159*q24
      q161 = pow2(x%d1val2_d1val3)
      q162 = q109*q130
      q163 = 32.0_dp*q108
      q164 = 2.0_dp*x%d2val2
      q165 = q128*q46
      q166 = 24.0_dp*q96
      q167 = 12.0_dp*q96
      q168 = q24*x%d1val1_d1val3
      q169 = 24.0_dp*q168
      q170 = x%d1val3*x%d2val1
      q171 = x%d2val1_d1val3*x%val
      q172 = 12.0_dp*x%d1val1_d2val3
      q173 = x%d2val1*x%val
      q174 = q109*q142
      q175 = 128.0_dp*q113
      q176 = q1*q108
      q177 = 192.0_dp*q176
      q178 = 96.0_dp*q108
      q179 = x%d1val1*x%d2val1
      q180 = 768.0_dp*q108
      q181 = q106*q176
      q182 = 768.0_dp*q121
      q183 = q182*x%d2val3
      q184 = 1536.0_dp*q123
      q185 = 4608.0_dp*q121
      q186 = 4352.0_dp*q122
      q187 = pow4(x%val)
      q188 = q187*q21*powm1(pow4(q23))
      q189 = 18432.0_dp*q188
      q190 = 12.0_dp*q24
      q191 = q1*q178
      q192 = q108*q5
      q193 = q99*x%d1val1_d1val2
      q194 = q99*x%d1val2
      q195 = x%d1val1*x%d1val1_d2val3
      q196 = q100*x%d1val1_d1val2
      q197 = 16.0_dp*x%d1val1_d1val2_d1val3*x%val
      q198 = 4.0_dp*x%d1val2_d2val3
      q199 = q108*q16
      q200 = 192.0_dp*q113
      q201 = x%d1val1*x%d1val1_d1val2
      q202 = 64.0_dp*q176
      q203 = q1*q118
      q204 = 64.0_dp*x%d1val2
      q205 = 1024.0_dp*q123
      q206 = 512.0_dp*q121
      q207 = q4*q81
      q208 = 3072.0_dp*q122
      q209 = 1.5_dp*q145
      q210 = x%d1val1_d1val2*x%val
      q211 = 2.0_dp*x%d2val1_d1val2 - q194*q4 + q202*q36 - q210*q97
      q212 = 0.5_dp*q46
      q213 = q21*q211
      q214 = q24*q71
      q215 = q24*x%d1val2
      q216 = q24*x%d1val2_d1val3
      q217 = q100*x%d1val2
      q218 = q1*q24
      q219 = q159 + q20*x%d2val2 - q218*q50
      q220 = q219*q24
      q221 = 8.0_dp*q219
      q222 = q221*x%d1val1
      q223 = 48.0_dp*q218
      q224 = q24*q50
      q225 = 192.0_dp*q108
      q226 = 4.0_dp*q68 + q139*q225 + q164*x%d1val3 - q223*q68 - q224*q40 + q45*x%d2val2_d1val3
      q227 = 2.0_dp*q96
      q228 = 24.0_dp*q215
      q229 = q190*x%d2val2
      q230 = q190*x%d1val2
      q231 = q17*q190
      q232 = 288.0_dp*x%d1val2
      q233 = x%d1val2*x%d2val2
      unary%val = q0*asin(x%val)
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d1val3 = q3*x%d1val3
      unary%d2val1 = q10*q3
      unary%d1val1_d1val2 = q11*q14 + q3*x%d1val1_d1val2
      unary%d1val1_d1val3 = q15*q16 + q3*x%d1val1_d1val3
      unary%d2val2 = q18*q3
      unary%d1val2_d1val3 = q19*x%d1val3 + q3*x%d1val2_d1val3
      unary%d2val3 = q3*(-q20*q25 + x%d2val3)
      unary%d3val1 = q3*q32
      unary%d2val1_d1val2 = q10*q19 + q3*q39
      unary%d2val1_d1val3 = q10*q41 + q3*q44
      unary%d1val1_d2val2 = q3*q55
      unary%d1val1_d1val2_d1val3 = q11*q57 + q14*q56 + q19*x%d1val1_d1val3 + q3*x%d1val1_d1val2_d1val3 + q41*x%d1val1_d1val2 + q56*q58
      unary%d1val1_d2val3 = q3*(q53*q62 + q59*q60 + x%d1val1_d2val3)
      unary%d3val2 = q3*q67
      unary%d2val2_d1val3 = q18*q41 + q3*(q38*q70 - q68*q69 - q70*q8 + x%d2val2_d1val3)
      unary%d1val2_d2val3 = q3*(q62*q73 + q71*q72 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q3*(-q11*q27*x%d2val1_d1val3 - q27*q75 - q28*q56 - q28*q74 + q75*q80 + q76*q77 + q78*q79 - q83*q85 + x%d3val1_d1val3) + q32*q41
      unary%d2val1_d1val2_d1val3 = q10*q57*x%val + q14*q86 + q19*q44 + q3*(-23887872.0_dp*q84*q94 + 124416.0_dp*q30*q92 - q33*q87 - q33*q89 - q35*x%d1val2_d1val3 - q42*x%d1val1_d1val2_d1val3 + q43*x%d1val2_d1val3 - q69*q90 + q87*q91 + q89*q91 + x%d2val1_d1val2_d1val3) + q39*q41 + q58*q86
      unary%d2val1_d2val3 = q3*(-q102 - q105 - q107 + q112 + q116 + q120 - q125 + q127*q128 + q131*q132 - q98 + x%d2val1_d2val3)
      unary%d1val1_d2val2_d1val3 = q3*(0.020833333333333332_dp*q52*q59 + 0.041666666666666664_dp*q136*q16*q52 + q133*q48 + q134*q47 + q135*q137*q22 + q45*q73*x%d1val1_d1val2_d1val3 + q54*(-41472.0_dp*q1*q68*q7 - 41472.0_dp*q40*q51 + 48.0_dp*q138 + 5971968.0_dp*q139*q30 + 96.0_dp*q68 + q49*x%d1val3) + x%d1val1_d2val2_d1val3) + q41*q55
      unary%d1val1_d1val2_d2val3 = q3*(15.0_dp*q149*x%d1val1 + 2.0_dp*q53*q71 + 2.0_dp*q73*q95 + 9.0_dp*q11*q144 + q132*x%d1val1_d1val2_d1val3 + q134*q59 + q140*q146*x%d1val1 + q140*q53 + q141*q53 + q142*q47 + q143*x%d1val1_d2val3 + q146*q21*x%d1val1_d1val2 + q147*q71*x%d1val1 + q147*q95*x%d1val2 + q21*q47 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q3*(-432.0_dp*q51*x%d1val2_d1val3 - q138*q64 - q150*q64 + q150*q78*x%d1val2 - q151*q27 + q152*q77 + q153*q80 - q155*q85 + x%d3val2_d1val3) + q41*q67
      unary%d2val2_d2val3 = q3*(-q100*q158 - q100*q161 + q114*q17*x%val + q119*q163*q17 - q124*q17 + q132*(q1*q163*q70 - q100*q68 - q160*x%d1val3 + x%d2val2_d1val3) - q156*q157 - q160*x%d2val3 + q162*q68 + q165*(-q100*q17 + q164) + x%d2val2_d2val3)
      unary%d3val1_d2val3 = q3*(288.0_dp*q11*q113*x%d2val1 + 288.0_dp*q181*x%d1val1 - q103*q172 - q106*q166 + q117*q5*x%d1val1_d2val3 + q119*q178*q179 + q132*(q109*q76 - q167*q170 - q167*q171 - q182*q83 - q190*q74*x%d2val1 - q190*q75 + q191*q79 + q192*q75 + x%d3val1_d1val3) + q165*(2.0_dp*x%d3val1 - q166*q173 + q178*q29 - q26*q99) - q166*x%d1val3*x%d2val1_d1val3 - q167*x%d2val1*x%d2val3 - q167*x%d2val1_d2val3*x%val - q169*q170 - q169*q171 - q172*q173*q24 + q174*q26 + q175*q26 + q177*q56*x%d2val1_d1val3 + q177*q95*x%d2val1 - q179*q184 + q180*q40*q75 - q183*q82 - q185*q75*q93 - q186*q29 + q189*q26 + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q3*(-1024.0_dp*q121*q207*q71 - 2048.0_dp*q111*q121*q148 + 0.5_dp*q126*q73*x%d2val3 + 0.5_dp*q127*q141 + 12288.0_dp*q188*q36 + 192.0_dp*q117*q71*x%val + 384.0_dp*q199*x%d1val1_d1val3*x%d1val2 + 4.5_dp*q126*q144*x%val + 6.0_dp*q130*q131*q137 + 64.0_dp*q108*q119*q201 + 7.5_dp*q126*q149 - q1*q208*q36 - q101*x%d1val1_d1val2_d2val3 - q103*q198 - q106*q194 + q110*q56*x%d1val1_d1val2_d1val3 + q110*q88*x%d1val2_d1val3 + q114*q36 + q115*q140*q178 + q126*q140*q209 + q126*q146*q71 + q127*q71 + q131*q133*q73 + q131*q72*x%d1val2_d1val3 + q132*(-q100*q90 - q101*x%d1val1_d1val2_d1val3 - q104*x%d1val2_d1val3 + q178*q92 - q194*q88 + q202*q87 + q202*q89 + q203*x%d1val2_d1val3 - q206*q94 - q87*q99 + x%d2val1_d1val2_d1val3) - q140*q206*q207 + q142*q211*q212 - q143*(q102 + q105 + q107 - q112 - q116 - q120 + q125 + q98 - x%d2val1_d2val3) - q156*q168*x%d1val1_d1val2 + q162*q90 - q168*q197 + q181*q204 - q193*x%d1val1*x%d2val3 - q194*q195 + q195*q202*x%d1val2 - q196*x%d1val1_d2val3 + q200*q34 - q201*q205 + q203*x%d1val2_d2val3 + q209*q213 + q212*q213 - q97*x%d1val1_d1val2_d1val3*x%d1val3 - q97*x%d1val1_d1val3*x%d1val2_d1val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q3*(-128.0_dp*q1*q122*q219*x%d1val1 - 16.0_dp*q214*x%d1val1_d1val2 - 4.0_dp*q168*q226 + 16.0_dp*q108*q219*q40*x%d1val1_d1val3 + q108*q142*q222 + q110*q71*x%d1val1_d1val2 + q113*q222 - q132*(-q108*q130*q204*x%d1val1_d1val2 + q135*q194 + q196*x%d1val2_d1val3 - q199*q221 + q217*x%d1val1_d1val2_d1val3 + q220*x%d1val1_d1val3 + q226*q227 - x%d1val1_d2val2_d1val3) - q140*q193 + q140*q202*x%d1val1_d1val2 - q156*q215*x%d1val1_d1val2_d1val3 + q16*q163*q226 + q162*x%d1val1_d1val2_d1val3*x%d1val2 + q165*(-16.0_dp*q210*q215 - 2.0_dp*q220*x%d1val1 + 2.0_dp*x%d1val1_d2val2) - q196*x%d1val2_d2val3 - q197*q216 + q200*q210*x%d1val2 - q205*x%d1val1_d1val2*x%d1val2 - q217*x%d1val1_d1val2_d2val3 - q220*x%d1val1_d2val3 - q227*(-192.0_dp*q157*q40 + 4.0_dp*q161 + 960.0_dp*q1*q113*q17 + q129*x%d2val2_d1val3 - q142*q224 - q158*q223 - q161*q223 + q164*x%d2val3 - q17*q187*q208 + q17*q225*q81*x%d2val3 + q180*q68*q93 + q198*x%d1val2 - q25*q50 + q45*x%d2val2_d2val3) + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q3*(-24.0_dp*q138*q216 - 24.0_dp*q214*x%d2val2 + q113*q232*q65 + q130*q225*x%d1val2*x%d2val2_d1val3 + q132*(q109*q152 + q130*q178*q233 - q138*q230 - q150*q230 - q151*q190 + q153*q192 - q155*q182 - q231*x%d1val2_d1val3 + x%d3val2_d1val3) - q139*q185*x%d1val2_d1val3 + q140*q191*x%d2val2 - q140*q229 - q141*q229 - q154*q183 + q161*q176*q232 - q161*q228 + q165*(2.0_dp*x%d3val2 + q178*q66 - q228*q65 - q63*q99) + q17*q192*x%d1val2_d2val3 + q174*q63 + q175*q63 + q177*q71*x%d2val2 + q180*q70*x%d1val2_d1val3*x%val - q184*q233 - q186*q66 + q189*q63 - q228*x%d1val3*x%d2val2_d1val3 - q230*x%d2val2_d2val3*x%val - q231*x%d1val2_d2val3 + x%d3val2_d2val3)
   end function asinpi_self
   
   function acospi_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q5 = 144.0_dp*q1
      q6 = -144.0_dp + q5
      q7 = powm1(q6)
      q8 = 144.0_dp*q7
      q9 = q8*x%val
      q10 = -q4*q9 + x%d2val1
      q11 = x%d1val1*x%val
      q12 = powm1(pow3(sqrt(q2)))
      q13 = q0*x%d1val2
      q14 = q12*q13
      q15 = q0*q12
      q16 = q11*x%d1val3
      q17 = pow2(x%d1val2)
      q18 = -q17*q9 + x%d2val2
      q19 = q14*x%val
      q20 = 4.0_dp*x%val
      q21 = pow2(x%d1val3)
      q22 = 4.0_dp*q1
      q23 = -4.0_dp + q22
      q24 = powm1(q23)
      q25 = q21*q24
      q26 = pow3(x%d1val1)
      q27 = 432.0_dp*q7
      q28 = q27*x%d2val1
      q29 = q1*q26
      q30 = powm1(pow2(q6))
      q31 = 62208.0_dp*q30
      q32 = -q11*q28 - q26*q8 + q29*q31 + x%d3val1
      q33 = 288.0_dp*q7
      q34 = q11*x%d1val1_d1val2
      q35 = q4*q8
      q36 = q4*x%d1val2
      q37 = q1*q30
      q38 = 41472.0_dp*q37
      q39 = -q33*q34 - q35*x%d1val2 + q36*q38 + x%d2val1_d1val2
      q40 = x%d1val3*x%val
      q41 = q15*q40
      q42 = q11*q33
      q43 = q38*q4
      q44 = -q35*x%d1val3 - q42*x%d1val1_d1val3 + q43*x%d1val3 + x%d2val1_d1val3
      q45 = 2.0_dp*x%val
      q46 = powm1(q2)
      q47 = q46*x%d1val1_d1val2
      q48 = q47*x%d1val2
      q49 = 48.0_dp*x%d2val2
      q50 = 48.0_dp*q17
      q51 = q17*q7
      q52 = -20736.0_dp*q1*q51 + q49*x%val + q50
      q53 = q46*x%d1val1
      q54 = 0.020833333333333332_dp*q53
      q55 = q45*q48 + q52*q54 + x%d1val1_d2val2
      q56 = x%d1val1*x%d1val3
      q57 = q15*x%d1val2_d1val3
      q58 = 3.0_dp*q1*q13*powm1(pow5(sqrt(q2)))
      q59 = q46*x%d1val1_d1val3
      q60 = q45*x%d1val3
      q61 = -24.0_dp*q1*q25 + 2.0_dp*q21 + q45*x%d2val3
      q62 = 0.5_dp*q61
      q63 = pow3(x%d1val2)
      q64 = q27*x%d1val2
      q65 = x%d2val2*x%val
      q66 = q1*q63
      q67 = q31*q66 - q63*q8 - q64*q65 + x%d3val2
      q68 = x%d1val2*x%d1val2_d1val3
      q69 = q33*x%val
      q70 = q17*x%d1val3
      q71 = x%d1val2_d1val3*x%d1val3
      q72 = q45*q46
      q73 = q46*x%d1val2
      q74 = x%d1val1_d1val3*x%val
      q75 = q4*x%d1val1_d1val3
      q76 = q26*q40
      q77 = 165888.0_dp*q30
      q78 = 124416.0_dp*q37
      q79 = q56*x%d2val1
      q80 = 186624.0_dp*q37
      q81 = pow3(x%val)
      q82 = q26*q81
      q83 = q82*x%d1val3
      q84 = powm1(pow3(q6))
      q85 = 35831808.0_dp*q84
      q86 = q10*x%d1val3
      q87 = q56*x%d1val1_d1val2
      q88 = x%d1val1*x%d1val1_d1val3
      q89 = q88*x%d1val2
      q90 = x%d1val1_d1val2*x%d1val1_d1val3
      q91 = 82944.0_dp*q37
      q92 = q36*q40
      q93 = q81*x%d1val3
      q94 = q36*q93
      q95 = 8.0_dp*q24
      q96 = q95*x%val
      q97 = 2.0_dp*x%d2val1 - q4*q96
      q98 = q46*q97
      q99 = 0.25_dp*q61
      q100 = 4.0_dp*x%d1val3
      q101 = q24*q4
      q102 = powm1(pow2(q23))
      q103 = q1*q102
      q104 = 32.0_dp*q103
      q105 = q104*q4
      q106 = -q100*q101 + q105*x%d1val3 - q88*q96 + x%d2val1_d1val3
      q107 = q46*q60
      q108 = x%d1val1_d1val3*x%d1val3
      q109 = q24*x%d1val1
      q110 = 16.0_dp*q109
      q111 = q96*x%d1val1
      q112 = 4.0_dp*q101
      q113 = pow2(x%d1val1_d1val3)
      q114 = 128.0_dp*q103
      q115 = q56*x%d1val1_d1val3
      q116 = q102*q21
      q117 = 96.0_dp*q116
      q118 = q4*x%val
      q119 = powm1(pow3(q23))
      q120 = q119*q21
      q121 = q120*q81
      q122 = 512.0_dp*q121
      q123 = -q105*x%d2val3 + q108*q110 + q111*x%d1val1_d2val3 + q112*x%d2val3 + q113*q96 - q114*q115 - q117*q118 + q122*q4 - x%d2val1_d2val3
      q124 = 2.0_dp*x%d1val3
      q125 = q45*x%d1val2_d1val3
      q126 = x%d1val1_d1val2*x%d1val3
      q127 = powm1(pow2(q2))
      q128 = q127*x%d1val2
      q129 = x%d2val2_d1val3*x%val
      q130 = q70*q81
      q131 = x%d1val2*x%d2val3
      q132 = x%d1val2_d2val3*x%val
      q133 = x%d2val3*x%val
      q134 = q73*x%val
      q135 = q128*q21
      q136 = q1*q127
      q137 = 3.0_dp*q136
      q138 = 6.0_dp*q136
      q139 = q138*x%d1val2
      q140 = q81*x%d1val2
      q141 = q140*q21*powm1(pow3(q2))
      q142 = x%d1val3*x%d2val2
      q143 = q65*x%d1val2_d1val3
      q144 = q40*q63
      q145 = q142*x%d1val2
      q146 = q17*x%d1val2_d1val3
      q147 = q63*q81
      q148 = q147*x%d1val3
      q149 = 16.0_dp*x%d1val3
      q150 = q24*q68
      q151 = x%d1val2*x%d1val2_d2val3
      q152 = 4.0_dp*q17
      q153 = q152*q24
      q154 = pow2(x%d1val2_d1val3)
      q155 = q114*x%d1val3
      q156 = q17*x%d2val3
      q157 = 2.0_dp*x%d2val2
      q158 = q46*q99
      q159 = 24.0_dp*q109
      q160 = 12.0_dp*q109
      q161 = x%d2val1*x%d2val3
      q162 = q160*x%val
      q163 = x%d1val3*x%d2val1
      q164 = q24*x%d1val1_d1val3
      q165 = 24.0_dp*q164
      q166 = 12.0_dp*x%d1val1_d2val3
      q167 = x%d2val1*x%val
      q168 = 128.0_dp*q102
      q169 = q133*q168
      q170 = 128.0_dp*q116
      q171 = 192.0_dp*q103
      q172 = 288.0_dp*q116
      q173 = 96.0_dp*q102
      q174 = q1*q173
      q175 = 768.0_dp*q102
      q176 = q103*q113
      q177 = q102*q5
      q178 = 768.0_dp*q119
      q179 = q178*x%d2val3
      q180 = q121*x%d1val1
      q181 = 4608.0_dp*q119
      q182 = 4352.0_dp*q120
      q183 = pow4(x%val)
      q184 = q183*q21*powm1(pow4(q23))
      q185 = 18432.0_dp*q184
      q186 = 12.0_dp*q24
      q187 = q95*x%d1val1_d1val2
      q188 = x%d1val1*x%d2val3
      q189 = x%d1val1_d1val2_d1val3*x%d1val3
      q190 = q95*x%d1val2
      q191 = x%d1val1*x%d1val1_d2val3
      q192 = q96*x%d1val1_d1val2
      q193 = 16.0_dp*x%d1val1_d1val2_d1val3*x%val
      q194 = 4.0_dp*x%d1val2_d2val3
      q195 = q102*q16
      q196 = 192.0_dp*q116
      q197 = 64.0_dp*q103
      q198 = q197*x%d1val1_d1val2
      q199 = q197*x%d1val2
      q200 = 192.0_dp*q102
      q201 = 1024.0_dp*x%d1val1_d1val2
      q202 = 512.0_dp*q119
      q203 = q4*q81
      q204 = 3072.0_dp*q120
      q205 = 1.5_dp*q136
      q206 = x%d1val1_d1val2*x%val
      q207 = 2.0_dp*x%d2val1_d1val2 - q110*q206 - q190*q4 + q197*q36
      q208 = 0.5_dp*q46
      q209 = q207*q21
      q210 = q24*q71
      q211 = q24*x%d1val2
      q212 = q24*x%d1val2_d1val3
      q213 = q96*x%d1val2
      q214 = q103*x%d1val2
      q215 = q121*x%d1val2
      q216 = q1*q24
      q217 = q152 + q20*x%d2val2 - q216*q50
      q218 = q217*q24
      q219 = 8.0_dp*q217
      q220 = q219*x%d1val1
      q221 = 48.0_dp*q216
      q222 = q24*q50
      q223 = 4.0_dp*q68 + q130*q200 + q157*x%d1val3 - q221*q68 - q222*q40 + q45*x%d2val2_d1val3
      q224 = 2.0_dp*q109
      q225 = 24.0_dp*q211
      q226 = x%d1val3*x%d2val2_d1val3
      q227 = q186*x%d2val2
      q228 = q186*x%d1val2
      q229 = q17*q186
      unary%val = q0*acos(x%val)
      unary%d1val1 = -q3*x%d1val1
      unary%d1val2 = -q3*x%d1val2
      unary%d1val3 = -q3*x%d1val3
      unary%d2val1 = -q10*q3
      unary%d1val1_d1val2 = -q11*q14 - q3*x%d1val1_d1val2
      unary%d1val1_d1val3 = -q15*q16 - q3*x%d1val1_d1val3
      unary%d2val2 = -q18*q3
      unary%d1val2_d1val3 = -q19*x%d1val3 - q3*x%d1val2_d1val3
      unary%d2val3 = -q3*(-q20*q25 + x%d2val3)
      unary%d3val1 = -q3*q32
      unary%d2val1_d1val2 = -q10*q19 - q3*q39
      unary%d2val1_d1val3 = -q10*q41 - q3*q44
      unary%d1val1_d2val2 = -q3*q55
      unary%d1val1_d1val2_d1val3 = -q11*q57 - q14*q56 - q19*x%d1val1_d1val3 - q3*x%d1val1_d1val2_d1val3 - q41*x%d1val1_d1val2 - q56*q58
      unary%d1val1_d2val3 = -q3*(q53*q62 + q59*q60 + x%d1val1_d2val3)
      unary%d3val2 = -q3*q67
      unary%d2val2_d1val3 = -q18*q41 - q3*(q38*q70 - q68*q69 - q70*q8 + x%d2val2_d1val3)
      unary%d1val2_d2val3 = -q3*(q62*q73 + q71*q72 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q3*(-q11*q27*x%d2val1_d1val3 - q27*q75 - q28*q56 - q28*q74 + q75*q80 + q76*q77 + q78*q79 - q83*q85 + x%d3val1_d1val3) - q32*q41
      unary%d2val1_d1val2_d1val3 = -q10*q57*x%val - q14*q86 - q19*q44 - q3*(-23887872.0_dp*q84*q94 + 124416.0_dp*q30*q92 - q33*q87 - q33*q89 - q35*x%d1val2_d1val3 - q42*x%d1val1_d1val2_d1val3 + q43*x%d1val2_d1val3 - q69*q90 + q87*q91 + q89*q91 + x%d2val1_d1val2_d1val3) - q39*q41 - q58*q86
      unary%d2val1_d2val3 = q3*(-q106*q107 + q123 - q98*q99)
      unary%d1val1_d2val2_d1val3 = -q3*(0.020833333333333332_dp*q52*q59 + 0.041666666666666664_dp*q127*q16*q52 + q124*q48 + q125*q47 + q126*q128*q22 + q45*q73*x%d1val1_d1val2_d1val3 + q54*(-41472.0_dp*q1*q68*q7 - 41472.0_dp*q40*q51 + 48.0_dp*q129 + 5971968.0_dp*q130*q30 + 96.0_dp*q68 + q49*x%d1val3) + x%d1val1_d2val2_d1val3) - q41*q55
      unary%d1val1_d1val2_d2val3 = -q3*(15.0_dp*q141*x%d1val1 + 2.0_dp*q108*q73 + 2.0_dp*q53*q71 + 9.0_dp*q11*q135 + q107*x%d1val1_d1val2_d1val3 + q108*q139 + q125*q59 + q131*q137*x%d1val1 + q131*q53 + q132*q53 + q133*q47 + q134*x%d1val1_d2val3 + q137*q21*x%d1val1_d1val2 + q138*q71*x%d1val1 + q21*q47 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = -q3*(-432.0_dp*q51*x%d1val2_d1val3 - q129*q64 - q142*q64 - q143*q27 + q144*q77 + q145*q78 + q146*q80 - q148*q85 + x%d3val2_d1val3) - q41*q67
      unary%d2val2_d2val3 = q3*(-q104*q156 - q107*(q104*q70 - q153*x%d1val3 - q68*q96 + x%d2val2_d1val3) - q117*q17*x%val + q122*q17 + q149*q150 + q151*q96 + q153*x%d2val3 + q154*q96 - q155*q68 - q158*(q157 - q17*q96) - x%d2val2_d2val3)
      unary%d3val1_d2val3 = q3*(-288.0_dp*q176*x%d1val1 + 1536.0_dp*q180*x%d2val1 + q101*q166 - q107*(-q160*q163 - q162*x%d2val1_d1val3 + q168*q76 + q174*q79 + q177*q75 - q178*q83 - q186*q74*x%d2val1 - q186*q75 + x%d3val1_d1val3) - q108*q171*x%d2val1 - q11*q172*x%d2val1 + q113*q159 - q158*(2.0_dp*x%d3val1 - q159*q167 + q173*q29 - q26*q95) + q159*x%d1val3*x%d2val1_d1val3 + q160*q161 - q161*q174*x%d1val1 + q162*x%d2val1_d2val3 + q163*q165 + q165*x%d2val1_d1val3*x%val + q166*q167*q24 - q169*q26 - q170*q26 - q171*q56*x%d2val1_d1val3 - q175*q40*q75 - q177*q4*x%d1val1_d2val3 + q179*q82 + q181*q75*q93 + q182*q29 - q185*q26 - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q3*(-0.5_dp*q132*q98 - 0.5_dp*q73*q97*x%d2val3 - 12288.0_dp*q184*q36 - 384.0_dp*q195*x%d1val1_d1val3*x%d1val2 - 4.5_dp*q135*q97*x%val - 64.0_dp*q176*x%d1val2 - 7.5_dp*q141*q97 + 1024.0_dp*q119*q203*q71 + 2048.0_dp*q115*q119*q140 + q1*q204*q36 + q101*q194 - q105*x%d1val2_d2val3 - q106*q124*q73 - q106*q139*x%d1val3 - q106*q72*x%d1val2_d1val3 - q107*(q105*x%d1val2_d1val3 - q111*x%d1val1_d1val2_d1val3 - q112*x%d1val2_d1val3 + q173*q92 - q190*q88 + q197*q87 + q197*q89 - q202*q94 - q87*q95 - q90*q96 + x%d2val1_d1val2_d1val3) + q110*q189 + q110*x%d1val1_d1val3*x%d1val2_d1val3 + q111*x%d1val1_d1val2_d2val3 + q113*q190 - q114*q56*x%d1val1_d1val2_d1val3 - q114*q88*x%d1val2_d1val3 - q117*q36 - q118*q131*q173 - q118*q200*q71 + q123*q134 + q131*q202*q203 - q131*q205*q97 - q133*q207*q208 - q137*q71*q97 + q149*q164*x%d1val1_d1val2 - q155*q90 + q164*q193 + q180*q201 + q187*q188 - q188*q198 + q190*q191 - q191*q199 + q192*x%d1val1_d2val3 - q196*q34 - q205*q209 - q208*q209 - q71*q98 - x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q3*(-128.0_dp*q189*q214 - 16.0_dp*q102*q217*q40*x%d1val1_d1val3 - 32.0_dp*q195*q223 + 128.0_dp*q1*q120*q217*x%d1val1 + 16.0_dp*q210*x%d1val1_d1val2 + 4.0_dp*q164*q223 - q102*q133*q220 + q107*(q126*q190 - q126*q199 + q192*x%d1val2_d1val3 - q195*q219 + q213*x%d1val1_d1val2_d1val3 + q218*x%d1val1_d1val3 + q223*q224 - x%d1val1_d2val2_d1val3) - q114*q71*x%d1val1_d1val2 - q116*q220 + q131*q187 - q131*q198 + q149*q211*x%d1val1_d1val2_d1val3 - q158*(-16.0_dp*q206*q211 - 2.0_dp*q218*x%d1val1 + 2.0_dp*x%d1val1_d2val2) + q192*x%d1val2_d2val3 + q193*q212 - q196*q206*x%d1val2 + q201*q215 + q213*x%d1val1_d1val2_d2val3 + q218*x%d1val1_d2val3 + q224*(-192.0_dp*q150*q40 + 4.0_dp*q154 + 960.0_dp*q1*q116*q17 + q100*x%d2val2_d1val3 - q133*q222 - q151*q221 - q154*q221 + q156*q200*q81 + q157*x%d2val3 - q17*q183*q204 + q175*q68*q93 + q194*x%d1val2 - q25*q50 + q45*x%d2val2_d2val3) - x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q3*(-192.0_dp*q214*q226 - 288.0_dp*q154*q214 + 1536.0_dp*q215*x%d2val2 + 24.0_dp*q129*q212 + 24.0_dp*q210*x%d2val2 - q107*(-q129*q228 - q142*q228 - q143*q186 + q144*q168 + q145*q174 + q146*q177 - q148*q178 - q229*x%d1val2_d1val3 + x%d3val2_d1val3) + q130*q181*x%d1val2_d1val3 - q131*q174*x%d2val2 + q131*q227 + q132*q227 + q147*q179 + q154*q225 - q158*(2.0_dp*x%d3val2 + q173*q66 - q225*q65 - q63*q95) - q169*q63 - q17*q177*x%d1val2_d2val3 - q170*q63 - q171*q71*x%d2val2 - q172*q65*x%d1val2 - q175*q70*x%d1val2_d1val3*x%val + q182*q66 - q185*q63 + q225*q226 + q228*x%d2val2_d2val3*x%val + q229*x%d1val2_d2val3 - x%d3val2_d2val3)
   end function acospi_self
   
   function atanpi_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q233
      real(dp) :: q232
      real(dp) :: q231
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q2 = q1 + 1
      q3 = q0*powm1(q2)
      q4 = 2.0_dp*x%d2val1
      q5 = pow2(x%d1val1)
      q6 = 144.0_dp*q1 + 144.0_dp
      q7 = powm1(q6)
      q8 = 576.0_dp*q7
      q9 = q8*x%val
      q10 = q0*(q4 - q5*q9)
      q11 = 72.0_dp*q7
      q12 = x%d1val1*x%d1val2
      q13 = 2.0_dp*x%val
      q14 = q0*powm1(pow2(q2))
      q15 = q13*q14
      q16 = x%d1val1*x%d1val3
      q17 = 2.0_dp*x%d2val2
      q18 = pow2(x%d1val2)
      q19 = q17 - q18*q9
      q20 = q0*q11
      q21 = x%d1val2*x%d1val3
      q22 = pow2(x%d1val3)
      q23 = 4.0_dp*q1 + 4.0_dp
      q24 = powm1(q23)
      q25 = q22*q24
      q26 = 8.0_dp*q25
      q27 = 4.0_dp*q24
      q28 = q0*q27
      q29 = 2.0_dp*x%d3val1
      q30 = pow3(x%d1val1)
      q31 = q7*x%d1val1
      q32 = q31*x%val
      q33 = 1728.0_dp*x%d2val1
      q34 = q1*q30
      q35 = powm1(pow2(q6))
      q36 = 331776.0_dp*q35
      q37 = q29 - q30*q8 - q32*q33 + q34*q36
      q38 = q35*x%d1val2
      q39 = 20736.0_dp*x%val
      q40 = q10*q39
      q41 = 2.0_dp*x%d2val1_d1val2
      q42 = 1152.0_dp*q32
      q43 = q5*q8
      q44 = q1*q5
      q45 = 165888.0_dp*q44
      q46 = q38*q45 + q41 - q42*x%d1val1_d1val2 - q43*x%d1val2
      q47 = q35*x%d1val3
      q48 = 2.0_dp*x%d2val1_d1val3 - q42*x%d1val1_d1val3 - q43*x%d1val3 + q45*q47
      q49 = 2.0_dp*x%d1val1_d2val2
      q50 = q7*x%d1val2
      q51 = 1152.0_dp*x%d1val1_d1val2
      q52 = q51*x%val
      q53 = 48.0_dp*x%d2val2
      q54 = q53*x%val
      q55 = q18*q7
      q56 = -27648.0_dp*q1*q55 + 48.0_dp*q18 + q54
      q57 = 12.0_dp*q31
      q58 = q49 - q50*q52 - q56*q57
      q59 = q16*x%d1val2
      q60 = x%d1val1*x%d1val2_d1val3
      q61 = x%d1val1_d1val2*x%d1val3
      q62 = x%d1val1_d1val3*x%d1val2
      q63 = 16.0_dp*q24
      q64 = q63*x%val
      q65 = q64*x%d1val3
      q66 = -32.0_dp*q1*q25 + 2.0_dp*q22 + q13*x%d2val3
      q67 = q27*x%d1val1
      q68 = 2.0_dp*x%d3val2
      q69 = pow3(x%d1val2)
      q70 = x%d2val2*x%val
      q71 = 1728.0_dp*q50
      q72 = q1*q69
      q73 = q36*q72 + q68 - q69*q8 - q70*q71
      q74 = q0*q39
      q75 = q47*q74
      q76 = 2.0_dp*x%d2val2_d1val3
      q77 = x%d1val2*x%d1val2_d1val3
      q78 = q7*q77
      q79 = 1152.0_dp*x%val
      q80 = q18*x%d1val3
      q81 = q35*q80
      q82 = 4.0_dp*x%d1val2
      q83 = q24*q66
      q84 = q31*x%d1val3
      q85 = q7*x%d1val1_d1val3
      q86 = 1728.0_dp*q85
      q87 = x%d2val1*x%val
      q88 = q30*x%val
      q89 = 829440.0_dp*q47
      q90 = q1*x%d1val1
      q91 = 995328.0_dp*q35
      q92 = q44*x%d1val1_d1val3
      q93 = pow3(x%val)
      q94 = q30*x%d1val3
      q95 = q93*q94
      q96 = powm1(pow3(q6))
      q97 = 191102976.0_dp*q96
      q98 = q38*x%d1val3
      q99 = q35*x%d1val2_d1val3
      q100 = q1*q16
      q101 = q5*x%val
      q102 = q5*x%d1val2
      q103 = q102*x%d1val3
      q104 = q4 - q5*q64
      q105 = 2.0_dp*q83
      q106 = q64*x%d1val1
      q107 = 8.0_dp*q24
      q108 = q107*q5
      q109 = powm1(pow2(q23))
      q110 = 64.0_dp*q109
      q111 = q110*q44
      q112 = -q106*x%d1val1_d1val3 - q108*x%d1val3 + q111*x%d1val3 + x%d2val1_d1val3
      q113 = q24*x%d1val1_d1val3
      q114 = 32.0_dp*q113
      q115 = q5*x%d2val3
      q116 = pow2(x%d1val1_d1val3)
      q117 = 256.0_dp*q109
      q118 = q100*q117
      q119 = q109*q22
      q120 = 192.0_dp*q119
      q121 = q1*q110
      q122 = q22*q5
      q123 = powm1(pow3(q23))
      q124 = q123*q93
      q125 = 1024.0_dp*q124
      q126 = -q101*q120 + q106*x%d1val1_d2val3 + q107*q115 + q114*q16 - q115*q121 + q116*q64 - q118*x%d1val1_d1val3 + q122*q125 - x%d2val1_d2val3
      q127 = q1*x%d1val2
      q128 = q53*x%d1val3
      q129 = x%d2val2_d1val3*x%val
      q130 = 48.0_dp*q129
      q131 = x%d1val3*x%val
      q132 = q107*x%d2val3
      q133 = q16*q63
      q134 = q107*x%val
      q135 = x%d1val1_d1val2*x%val
      q136 = q63*x%d1val2
      q137 = x%d1val1_d1val3*x%d1val3
      q138 = q64*x%d1val2_d1val3
      q139 = q134*x%d1val2
      q140 = 384.0_dp*q119
      q141 = q1*q109
      q142 = 128.0_dp*q141
      q143 = q12*q142
      q144 = q1*q117
      q145 = q1*q119
      q146 = 128.0_dp*x%d1val1_d1val2
      q147 = q124*q22
      q148 = 3072.0_dp*q147
      q149 = x%d1val3*x%d2val2
      q150 = 1728.0_dp*x%d1val2_d1val3
      q151 = q69*x%val
      q152 = q18*x%d1val2_d1val3
      q153 = q1*q152
      q154 = q69*x%d1val3
      q155 = 32.0_dp*q24
      q156 = q77*x%d1val3
      q157 = x%d1val2*x%d1val2_d2val3
      q158 = pow2(x%d1val2_d1val3)
      q159 = q120*x%val
      q160 = 64.0_dp*q18
      q161 = q160*x%d2val3
      q162 = q18*q22
      q163 = 48.0_dp*q24
      q164 = q16*x%d2val1_d1val3
      q165 = 24.0_dp*q24
      q166 = q165*x%d2val1
      q167 = x%d1val1*x%d2val3
      q168 = q165*x%d1val1
      q169 = q137*x%d2val1
      q170 = x%d2val1_d1val3*x%val
      q171 = q24*x%d1val1_d2val3
      q172 = 24.0_dp*q87
      q173 = q116*x%d1val1
      q174 = 24.0_dp*q5
      q175 = 320.0_dp*q109
      q176 = q175*x%d2val3
      q177 = 320.0_dp*q119
      q178 = 384.0_dp*q141
      q179 = q87*x%d1val1
      q180 = 576.0_dp*q119
      q181 = 192.0_dp*x%d2val1
      q182 = q141*q167
      q183 = q101*q109
      q184 = 768.0_dp*q141
      q185 = 384.0_dp*q109
      q186 = 2048.0_dp*q124
      q187 = q186*x%d2val3
      q188 = 12288.0_dp*q124
      q189 = q123*q22
      q190 = 11264.0_dp*q189
      q191 = pow4(x%val)
      q192 = q191*powm1(pow4(q23))
      q193 = 49152.0_dp*q192*q22
      q194 = q141*q16
      q195 = q155*x%d1val1_d1val2_d1val3
      q196 = q136*x%d1val1
      q197 = q64*x%d1val1_d1val2
      q198 = q109*x%val
      q199 = q16*q62
      q200 = q119*x%d1val1
      q201 = q115*x%d1val2
      q202 = x%d1val2_d1val3*x%d1val3
      q203 = 128.0_dp*x%d1val2
      q204 = q186*q22*x%d1val1_d1val2
      q205 = q104*x%d2val3
      q206 = q104*q202
      q207 = q27*x%val
      q208 = q104*x%d1val2
      q209 = q112*x%d1val3
      q210 = q135*q155
      q211 = q109*q203*q44 - q136*q5 - q210*x%d1val1 + q41
      q212 = q61*x%d1val2_d1val3
      q213 = x%d1val2_d1val3*x%val
      q214 = q64*x%d1val2
      q215 = q141*q146
      q216 = q1*q24
      q217 = 4.0_dp*q18 + 4.0_dp*q70 - q160*q216
      q218 = 2.0_dp*q217
      q219 = 16.0_dp*q217
      q220 = q198*q219
      q221 = 64.0_dp*q216
      q222 = q24*x%val
      q223 = q117*q93
      q224 = -64.0_dp*q222*q80 + 4.0_dp*q77 + q17*x%d1val3 - q221*q77 + q223*q80 + q76*x%val
      q225 = x%d1val3*x%d2val2_d1val3
      q226 = q225*x%d1val2
      q227 = q165*x%d1val2
      q228 = x%d2val2*x%d2val3
      q229 = q24*x%d1val2_d1val3
      q230 = q165*q70
      q231 = q158*x%d1val2
      q232 = q18*x%d1val2_d2val3
      q233 = 192.0_dp*q141*x%d1val2
      unary%val = q0*atan(x%val)
      unary%d1val1 = q3*x%d1val1
      unary%d1val2 = q3*x%d1val2
      unary%d1val3 = q3*x%d1val3
      unary%d2val1 = q10*q11
      unary%d1val1_d1val2 = -q12*q15 + q3*x%d1val1_d1val2
      unary%d1val1_d1val3 = -q15*q16 + q3*x%d1val1_d1val3
      unary%d2val2 = q19*q20
      unary%d1val2_d1val3 = -q15*q21 + q3*x%d1val2_d1val3
      unary%d2val3 = q28*(-q26*x%val + x%d2val3)
      unary%d3val1 = q20*q37
      unary%d2val1_d1val2 = q20*q46 - q38*q40
      unary%d2val1_d1val3 = q20*q48 - q40*q47
      unary%d1val1_d2val2 = q20*q58
      unary%d1val1_d1val2_d1val3 = -2.0_dp*q14*q59 + 8.0_dp*q0*q1*q59*powm1(pow3(q2)) - q15*q60 - q15*q61 - q15*q62 + q3*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q28*(-q65*x%d1val1_d1val3 - q66*q67 + x%d1val1_d2val3)
      unary%d3val2 = q20*q73
      unary%d2val2_d1val3 = -q19*q75 + q20*(165888.0_dp*q1*q81 + q76 - q78*q79 - q8*q80)
      unary%d1val2_d2val3 = q28*(-q65*x%d1val2_d1val3 - q82*q83 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q20*(-1728.0_dp*q32*x%d2val1_d1val3 + 2.0_dp*x%d3val1_d1val3 + 497664.0_dp*q47*q90*x%d2val1 - q33*q84 - q5*q86 - q86*q87 + q88*q89 + q91*q92 - q95*q97) - q37*q75
      unary%d2val1_d1val2_d1val3 = -20736.0_dp*q10*q98 + 11943936.0_dp*q1*q10*q21*q96 + q20*(-1152.0_dp*q31*q62 - 95551488.0_dp*q103*q93*q96 + 2.0_dp*x%d2val1_d1val2_d1val3 + 497664.0_dp*q101*q98 + q100*q36*x%d1val1_d1val2 + q36*q62*q90 - q42*x%d1val1_d1val2_d1val3 - q43*x%d1val2_d1val3 + q45*q99 - q51*q84 - q52*q85) - q38*q48*q74 - q40*q99 - q46*q75
      unary%d2val1_d2val3 = -q28*(q104*q105 + q112*q65 + q126)
      unary%d1val1_d2val2_d1val3 = q20*(-12.0_dp*q56*q85 + 2.0_dp*x%d1val1_d2val2_d1val3 + 3456.0_dp*q47*q56*x%d1val1*x%val + q127*q36*q61 - q50*q51*x%d1val3 - q50*q79*x%d1val1_d1val2_d1val3 - q52*q7*x%d1val2_d1val3 - q57*(-55296.0_dp*q1*q78 - 55296.0_dp*q131*q55 + 7962624.0_dp*q81*q93 + 96.0_dp*q77 + q128 + q130)) - q58*q75
      unary%d1val1_d1val2_d2val3 = q28*(q118*x%d1val2_d1val3 - q12*q132 + q12*q140*x%val - q12*q148 - q132*q135 - q133*x%d1val2_d1val3 - q134*x%d1val1*x%d1val2_d2val3 - q136*q137 - q138*x%d1val1_d1val3 - q139*x%d1val1_d2val3 + q143*x%d2val3 + q144*q62*x%d1val3 + q145*q146 - q26*x%d1val1_d1val2 - q65*x%d1val1_d1val2_d1val3 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q20*(2.0_dp*x%d3val2_d1val3 + 497664.0_dp*q1*q149*q38 - q129*q71 - q149*q71 - q150*q55 - q150*q7*q70 + q151*q89 + q153*q91 - q154*q93*q97) - q73*q75
      unary%d2val2_d2val3 = -q28*(q105*(q17 - q18*q64) + q125*q162 + q132*q18 - q141*q161 - q144*q156 + q155*q156 + q157*q64 + q158*q64 - q159*q18 + q65*(-q107*q80 + q121*q80 - q64*q77 + x%d2val2_d1val3) - x%d2val2_d2val3)
      unary%d3val1_d2val3 = -q28*(-1920.0_dp*q137*q183 + q105*(q117*q34 - q163*q179 + q29 - q30*q63) + q137*q188*q5 + q148*x%d1val1*x%d2val1 + q163*q164 + q163*q169 + q163*q170*x%d1val1_d1val3 + q163*q173 - q164*q178 + q166*q167 + q168*x%d2val1_d2val3*x%val - q169*q178 + q171*q172 + q171*q174 - q173*q184 - q176*q88 - q177*q30 - q179*q180 - q181*q182 - q185*q44*x%d1val1_d2val3 + q187*q30 + q190*q34 - q193*q30 + q65*(-2048.0_dp*q123*q95 - q113*q172 - q113*q174 - q16*q166 - q168*q170 + q175*q94*x%val + q181*q194 + q185*q92 + x%d3val1_d1val3) - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q28*(-1536.0_dp*q147*q208 - 4.0_dp*q211*q25 - 4096.0_dp*q124*q199 - 6144.0_dp*q189*q44*x%d1val2 + 192.0_dp*q198*q201 + 24576.0_dp*q122*q192*x%d1val2 + 384.0_dp*q135*q200 + 64.0_dp*q145*q211 + 768.0_dp*q198*q199 + q101*q185*q202 + q102*q120 - q104*q207*x%d1val2_d2val3 - q106*x%d1val1_d1val2_d2val3 - q107*q206 - q108*x%d1val2_d2val3 + q111*x%d1val2_d2val3 - q112*q138 - q114*q60 - q114*q61 - q114*x%d1val1_d1val2_d1val3*x%val - q116*q136 + q116*q141*q203 + q117*q127*q209 + q118*x%d1val1_d1val2_d1val3 + q121*q205*x%d1val2 - q125*q201 + q126*q139 - q136*q209 + q137*q144*x%d1val1_d1val2 + q142*q206 + q143*x%d1val1_d2val3 + q144*q60*x%d1val1_d1val3 + q146*q182 + q159*q208 - q16*q195 - q167*q63*x%d1val1_d1val2 - q186*q202*q5 - q196*x%d1val1_d2val3 - q197*x%d1val1_d2val3 - q204*x%d1val1 - q205*q24*q82 - q207*q211*x%d2val3 - q65*(192.0_dp*q183*q21 - q103*q125 - q106*x%d1val1_d1val2_d1val3 - q108*x%d1val2_d1val3 + q111*x%d1val2_d1val3 - q133*x%d1val1_d1val2 + q142*q62*x%d1val1 + q146*q194 - q196*x%d1val1_d1val3 - q197*x%d1val1_d1val3 + x%d2val1_d1val2_d1val3) + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q28*(-256.0_dp*q189*q217*q90 + 32.0_dp*q109*q131*q217*x%d1val1_d1val3 - q105*(-q210*x%d1val2 - q217*q67 + q49) - q107*q224*x%d1val1_d1val3 + q110*q16*q224*x%val + q135*q140*x%d1val2 - q136*x%d1val1_d1val2*x%d2val3 + q144*q21*x%d1val1_d1val2_d1val3 + q144*q212 - q155*q212 + q167*q220 - q171*q218 - q195*q21 - q195*q213 - q197*x%d1val2_d2val3 + q200*q219 - q204*x%d1val2 - q214*x%d1val1_d1val2_d2val3 + q215*x%d1val2*x%d2val3 + q65*(q113*q218 + q136*q61 + q138*x%d1val1_d1val2 - q16*q220 - q21*q215 + q214*x%d1val1_d1val2_d1val3 + q224*q67 - x%d1val1_d2val2_d1val3) - q67*(-256.0_dp*q131*q24*q77 - 4096.0_dp*q123*q162*q191 + 1024.0_dp*q109*q156*q93 + 1280.0_dp*q145*q18 + 4.0_dp*q158 + 4.0_dp*q225 + q13*x%d2val2_d2val3 - q157*q221 - q158*q221 - q160*q25 - q161*q222 + q17*x%d2val3 + q18*q223*x%d2val3 + q82*x%d1val2_d2val3) + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = -q28*(-1920.0_dp*q109*q213*q80 + q105*(q117*q72 - q24*q54*x%d1val2 - q63*q69 + q68) + q128*q229 + q130*q229 + q148*x%d1val2*x%d2val2 - q149*q178*x%d1val2_d1val3 - q151*q176 + q163*q226 + q163*q231 + q165*q232 - q177*q69 - q178*q226 - q178*q232 - q180*q70*x%d1val2 - q184*q231 + q187*q69 + q188*q80*x%d1val2_d1val3 + q190*q72 - q193*q69 + q227*q228 + q227*x%d2val2_d2val3*x%val - q228*q233 + q230*x%d1val2_d2val3 + q65*(-q129*q227 + q131*q175*q69 - q149*q227 + q149*q233 - q152*q165 + q153*q185 - q154*q186 - q230*x%d1val2_d1val3 + x%d3val2_d1val3) - x%d3val2_d2val3)
   end function atanpi_self
   
   function asinh_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q3 = pow2(x%d1val1)
      q4 = powm1(q1)
      q5 = q4*x%val
      q6 = -q3*q5 + x%d2val1
      q7 = powm1(pow3(sqrt(q1)))
      q8 = q7*x%val
      q9 = x%d1val1*x%d1val2
      q10 = q8*x%d1val3
      q11 = pow2(x%d1val2)
      q12 = -q11*q5 + x%d2val2
      q13 = pow2(x%d1val3)
      q14 = pow3(x%d1val1)
      q15 = q5*x%d1val1
      q16 = 3.0_dp*x%d2val1
      q17 = q0*q14
      q18 = powm1(pow2(q1))
      q19 = 3.0_dp*q18
      q20 = -q14*q4 - q15*q16 + q17*q19 + x%d3val1
      q21 = q8*x%d1val2
      q22 = 2.0_dp*q15
      q23 = q3*q4
      q24 = 2.0_dp*x%d1val2
      q25 = q0*q3
      q26 = q18*q25
      q27 = -q22*x%d1val1_d1val2 - q23*x%d1val2 + q24*q26 + x%d2val1_d1val2
      q28 = 2.0_dp*x%d1val3
      q29 = -q22*x%d1val1_d1val3 - q23*x%d1val3 + q26*q28 + x%d2val1_d1val3
      q30 = q5*x%d1val1_d1val2
      q31 = 48.0_dp*x%d2val2
      q32 = 48.0_dp*q11
      q33 = 144.0_dp*q0
      q34 = 144.0_dp + q33
      q35 = powm1(q34)
      q36 = q0*q11
      q37 = -20736.0_dp*q35*q36 + q31*x%val + q32
      q38 = q4*x%d1val1
      q39 = 0.020833333333333332_dp*q38
      q40 = -q24*q30 - q37*q39 + x%d1val1_d2val2
      q41 = q7*x%d1val3
      q42 = q8*x%d1val2_d1val3
      q43 = q0*x%d1val3
      q44 = 3.0_dp*q43*powm1(pow5(sqrt(q1)))
      q45 = q28*q5
      q46 = 2.0_dp*x%val
      q47 = 4.0_dp*q0
      q48 = 4.0_dp + q47
      q49 = powm1(q48)
      q50 = 24.0_dp*q49
      q51 = q0*q13
      q52 = 2.0_dp*q13 + q46*x%d2val3 - q50*q51
      q53 = 0.5_dp*q52
      q54 = pow3(x%d1val2)
      q55 = q5*x%d1val2
      q56 = 3.0_dp*x%d2val2
      q57 = q0*q54
      q58 = q19*q57 - q4*q54 - q55*q56 + x%d3val2
      q59 = x%d1val2*x%d1val2_d1val3
      q60 = 2.0_dp*q5
      q61 = q11*q4
      q62 = q18*q36
      q63 = x%d1val2_d1val3*x%d1val3
      q64 = q4*x%d1val2
      q65 = 3.0_dp*x%d1val1_d1val3
      q66 = q14*x%d1val3
      q67 = q66*x%val
      q68 = 8.0_dp*q18
      q69 = x%d1val1*x%d1val3
      q70 = q18*q69
      q71 = q0*x%d2val1
      q72 = pow3(x%val)
      q73 = q72*powm1(pow3(q1))
      q74 = 12.0_dp*q73
      q75 = q6*x%d1val2
      q76 = q28*x%d1val1_d1val2
      q77 = q4*x%d1val1_d1val3
      q78 = 2.0_dp*q30
      q79 = q18*q47
      q80 = q9*x%d1val1_d1val3
      q81 = q3*x%val
      q82 = q81*x%d1val3
      q83 = 6.0_dp*q18
      q84 = q83*x%d1val2
      q85 = q3*x%d1val2
      q86 = q85*x%d1val3
      q87 = 8.0_dp*q49
      q88 = q87*x%val
      q89 = 2.0_dp*x%d2val1 - q3*q88
      q90 = q4*q89
      q91 = 0.25_dp*q52
      q92 = q88*x%d1val1
      q93 = 4.0_dp*x%d1val3
      q94 = q3*q49
      q95 = powm1(pow2(q48))
      q96 = 32.0_dp*q95
      q97 = q25*q96
      q98 = -q92*x%d1val1_d1val3 - q93*q94 + q97*x%d1val3 + x%d2val1_d1val3
      q99 = 16.0_dp*q49
      q100 = q99*x%d1val1_d1val3
      q101 = 4.0_dp*q94
      q102 = pow2(x%d1val1_d1val3)
      q103 = q0*q69
      q104 = 128.0_dp*q95
      q105 = q104*x%d1val1_d1val3
      q106 = q13*q95
      q107 = 96.0_dp*q106
      q108 = q13*q3
      q109 = powm1(pow3(q48))
      q110 = q109*q72
      q111 = 512.0_dp*q110
      q112 = q100*q69 + q101*x%d2val3 + q102*q88 - q103*q105 - q107*q81 + q108*q111 + q92*x%d1val1_d2val3 - q97*x%d2val3 - x%d2val1_d2val3
      q113 = x%d1val1_d1val2*x%d1val3
      q114 = x%d2val2_d1val3*x%val
      q115 = 41472.0_dp*q35
      q116 = q0*q59
      q117 = q11*x%d1val3
      q118 = q117*q72
      q119 = q9*x%d2val3
      q120 = q28*q64
      q121 = q60*x%d1val2_d1val3
      q122 = q13*q4
      q123 = q13*q9
      q124 = q18*x%val
      q125 = q0*q19
      q126 = q0*x%d1val1
      q127 = q43*q84
      q128 = 3.0_dp*x%d1val2_d1val3
      q129 = q54*x%d1val3
      q130 = q129*x%val
      q131 = q99*x%d1val3
      q132 = x%d1val2*x%d1val2_d2val3
      q133 = 4.0_dp*q11
      q134 = q133*q49
      q135 = pow2(x%d1val2_d1val3)
      q136 = q11*x%val
      q137 = q36*q96
      q138 = q11*q13
      q139 = 2.0_dp*x%d2val2
      q140 = q4*q91
      q141 = q50*x%d1val1
      q142 = 12.0_dp*q49
      q143 = q142*x%d2val1
      q144 = x%d1val1*x%d2val3
      q145 = q142*x%d1val1
      q146 = q50*x%d1val1_d1val3
      q147 = x%d1val3*x%d2val1
      q148 = x%d2val1_d1val3*x%val
      q149 = x%d2val1*x%val
      q150 = q142*q149
      q151 = 12.0_dp*q94
      q152 = q14*x%d2val3
      q153 = q104*x%val
      q154 = 128.0_dp*q106
      q155 = 192.0_dp*q95
      q156 = q106*x%d1val1
      q157 = 96.0_dp*q95
      q158 = q157*x%d2val3
      q159 = 768.0_dp*q95
      q160 = q0*q155
      q161 = 288.0_dp*q95
      q162 = q33*q95
      q163 = q162*q3
      q164 = 768.0_dp*q110
      q165 = q110*q13
      q166 = q165*x%d1val1
      q167 = q110*q3
      q168 = x%d1val1_d1val3*x%d1val3
      q169 = 4352.0_dp*q109
      q170 = pow4(x%val)
      q171 = q170*powm1(pow4(q48))
      q172 = 18432.0_dp*q13*q171
      q173 = q87*x%d1val1_d1val2
      q174 = q99*x%d1val1_d1val2_d1val3
      q175 = q87*q9
      q176 = q88*x%d1val1_d1val2
      q177 = 4.0_dp*x%d1val2_d2val3
      q178 = q87*x%d1val2
      q179 = q95*x%val
      q180 = q80*x%d1val3
      q181 = 192.0_dp*x%val
      q182 = q181*x%d1val1_d1val2
      q183 = 64.0_dp*q95
      q184 = x%d1val1_d1val2*x%d2val3
      q185 = q0*q183
      q186 = q158*x%d1val2
      q187 = q183*x%d1val2
      q188 = q0*q187
      q189 = 1024.0_dp*x%d1val1_d1val2
      q190 = 3072.0_dp*q109
      q191 = 0.5_dp*q89
      q192 = q89*x%d1val2
      q193 = q13*q192
      q194 = 1.5_dp*q18
      q195 = q99*x%d1val1_d1val2
      q196 = q195*x%val
      q197 = 2.0_dp*x%d2val1_d1val2 - q178*q3 + q187*q25 - q196*x%d1val1
      q198 = 0.5_dp*q197
      q199 = q157*x%d1val2
      q200 = x%d1val1_d1val2_d1val3*x%d1val2
      q201 = q106*x%d1val2
      q202 = q165*x%d1val2
      q203 = x%d2val2*x%val
      q204 = q32*q49
      q205 = 4.0_dp*q203 - q0*q204 + q133
      q206 = q205*q49
      q207 = 8.0_dp*q205
      q208 = q179*q207
      q209 = 48.0_dp*q49
      q210 = q204*x%val
      q211 = 4.0_dp*q59 - q116*q209 + q118*q155 + q139*x%d1val3 - q210*x%d1val3 + q46*x%d2val2_d1val3
      q212 = q211*q49
      q213 = 2.0_dp*x%d1val1
      q214 = q59*x%d1val3
      q215 = q0*q135
      q216 = q50*x%d1val2
      q217 = x%d1val3*x%d2val2_d1val3
      q218 = q142*x%d1val2
      q219 = q218*x%d2val2
      q220 = q63*x%d2val2
      q221 = q11*x%d1val2_d2val3
      q222 = q54*x%d2val3
      q223 = q142*x%d1val2_d1val3
      unary%val = asinh(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q2*q6
      unary%d1val1_d1val2 = q2*x%d1val1_d1val2 - q8*q9
      unary%d1val1_d1val3 = -q10*x%d1val1 + q2*x%d1val1_d1val3
      unary%d2val2 = q12*q2
      unary%d1val2_d1val3 = -q10*x%d1val2 + q2*x%d1val2_d1val3
      unary%d2val3 = q2*(-q13*q5 + x%d2val3)
      unary%d3val1 = q2*q20
      unary%d2val1_d1val2 = q2*q27 - q21*q6
      unary%d2val1_d1val3 = -q10*q6 + q2*q29
      unary%d1val1_d2val2 = q2*q40
      unary%d1val1_d1val2_d1val3 = -q10*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3 - q21*x%d1val1_d1val3 - q41*q9 - q42*x%d1val1 + q44*q9
      unary%d1val1_d2val3 = q2*(-q38*q53 - q45*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q2*q58
      unary%d2val2_d1val3 = -q10*q12 + q2*(q28*q62 - q59*q60 - q61*x%d1val3 + x%d2val2_d1val3)
      unary%d1val2_d2val3 = q2*(-q53*q64 - q60*q63 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q10*q20 + q2*(-3.0_dp*q15*x%d2val1_d1val3 + 6.0_dp*q70*q71 + 9.0_dp*q26*x%d1val1_d1val3 - q16*q38*x%d1val3 - q23*q65 - q5*q65*x%d2val1 - q66*q74 + q67*q68 + x%d3val1_d1val3)
      unary%d2val1_d1val2_d1val3 = -q10*q27 + q2*(-2.0_dp*q77*q9 - 8.0_dp*q73*q86 + 2.0_dp*q26*x%d1val2_d1val3 - q22*x%d1val1_d1val2_d1val3 - q23*x%d1val2_d1val3 - q38*q76 + q47*q70*x%d1val1_d1val2 - q78*x%d1val1_d1val3 + q79*q80 + q82*q84 + x%d2val1_d1val2_d1val3) - q21*q29 - q41*q75 - q42*q6 + q44*q75
      unary%d2val1_d2val3 = -q2*(q112 + q45*q98 + q90*q91)
      unary%d1val1_d2val2_d1val3 = -q10*q40 + q2*(-0.020833333333333332_dp*q37*q77 - 2.0_dp*q55*x%d1val1_d1val2_d1val3 + 0.041666666666666664_dp*q37*q70*x%val + q113*q79*x%d1val2 - q39*(48.0_dp*q114 + 5971968.0_dp*q118*powm1(pow2(q34)) + 96.0_dp*q59 - q115*q116 - q115*q117*x%val + q31*x%d1val3) - q64*q76 - q78*x%d1val2_d1val3 + x%d1val1_d2val2_d1val3)
      unary%d1val1_d1val2_d2val3 = q2*(-15.0_dp*q123*q73 - 2.0_dp*q38*q63 + 9.0_dp*q123*q124 + q119*q125 - q119*q4 - q120*x%d1val1_d1val3 - q121*x%d1val1_d1val3 - q122*x%d1val1_d1val2 + q126*q63*q83 + q127*x%d1val1_d1val3 - q15*x%d1val2_d2val3 + q19*q51*x%d1val1_d1val2 - q30*x%d2val3 - q45*x%d1val1_d1val2_d1val3 - q55*x%d1val1_d2val3 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = -q10*q58 + q2*(-3.0_dp*q55*x%d2val2_d1val3 + 9.0_dp*q62*x%d1val2_d1val3 + q127*x%d2val2 - q128*q5*x%d2val2 - q128*q61 - q129*q74 + q130*q68 - q56*q64*x%d1val3 + x%d3val2_d1val3)
      unary%d2val2_d2val3 = -q2*(-q104*q116*x%d1val3 - q107*q136 + q111*q138 + q131*q59 + q132*q88 + q134*x%d2val3 + q135*q88 - q137*x%d2val3 + q140*(-q11*q88 + q139) + q45*(-q134*x%d1val3 + q137*x%d1val3 - q59*q88 + x%d2val2_d1val3) - x%d2val2_d2val3)
      unary%d3val1_d2val3 = -q2*(-288.0_dp*q149*q156 + 1536.0_dp*q166*x%d2val1 + 4608.0_dp*q167*q168 - q102*q126*q161 + q102*q141 - q103*q155*x%d2val1_d1val3 - q126*q158*x%d2val1 + q13*q169*q17 - q14*q154 - q14*q172 + q140*(2.0_dp*x%d3val1 - q14*q87 - q141*q149 + q157*q17) + q141*x%d1val3*x%d2val1_d1val3 + q143*q144 + q145*x%d2val1_d2val3*x%val + q146*q147 + q146*q148 - q147*q160*x%d1val1_d1val3 + q150*x%d1val1_d2val3 + q151*x%d1val1_d2val3 - q152*q153 + q152*q164 - q159*q82*x%d1val1_d1val3 - q163*x%d1val1_d2val3 + q45*(q104*q67 - q143*q69 - q145*q148 - q150*x%d1val1_d1val3 - q151*x%d1val1_d1val3 + q157*q69*q71 + q163*x%d1val1_d1val3 - q164*q66 + x%d3val1_d1val3) - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(-1024.0_dp*q167*q63 - 2048.0_dp*q110*q180 - 7.5_dp*q193*q73 + 12288.0_dp*q108*q171*x%d1val2 + 384.0_dp*q179*q180 + 4.5_dp*q124*q193 + q0*q105*q113 + q0*q192*q194*x%d2val3 - q100*q113 - q100*x%d1val1*x%d1val2_d1val3 - q100*x%d1val1_d1val2_d1val3*x%val - q102*q178 + q102*q188 + q103*q104*x%d1val1_d1val2_d1val3 + q105*q126*x%d1val2_d1val3 + q107*q85 - q111*q85*x%d2val3 + q112*q55 - q120*q98 - q121*q98 - q122*q198 + q125*q63*q89 + q126*q183*q184 + q127*q98 - q144*q173 + q155*q63*q81 + q156*q182 - q166*q189 - q174*q69 - q175*x%d1val1_d2val3 - q176*x%d1val1_d2val3 - q177*q94 + q185*q9*x%d1val1_d2val3 + q186*q81 - q190*q51*q85 - q191*q5*x%d1val2_d2val3 - q191*q64*x%d2val3 + q194*q197*q51 - q198*q5*x%d2val3 - q45*(-q101*x%d1val2_d1val3 + q103*q183*x%d1val1_d1val2 - q111*q86 - q173*q69 - q175*x%d1val1_d1val3 - q176*x%d1val1_d1val3 + q185*q80 + q199*q82 - q92*x%d1val1_d1val2_d1val3 + q97*x%d1val2_d1val3 + x%d2val1_d1val2_d1val3) - q63*q90 - q92*x%d1val1_d1val2_d2val3 + q97*x%d1val2_d2val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(-128.0_dp*q109*q205*q51*x%d1val1 - 4.0_dp*q212*x%d1val1_d1val3 + 16.0_dp*q168*q179*q205 + q0*q104*q63*x%d1val1_d1val2 + q104*q200*q43 - q131*q200 - q140*(2.0_dp*x%d1val1_d2val2 - q196*x%d1val2 - q206*q213) + q144*q208 + q156*q207 - q174*x%d1val2_d1val3*x%val - q176*x%d1val2_d2val3 - q178*q184 + q182*q201 + q184*q188 - q189*q202 - q195*q63 - q206*x%d1val1_d2val3 + q211*q69*q96*x%val - q213*q49*(4.0_dp*q135 + 960.0_dp*q106*q36 - q0*q132*q209 + q11*q155*q72*x%d2val3 - q13*q204 - q138*q170*q190 + q139*x%d2val3 + q159*q214*q72 + q177*x%d1val2 - q181*q214*q49 - q209*q215 - q210*x%d2val3 + q46*x%d2val2_d2val3 + q93*x%d2val2_d1val3) + q45*(q113*q178 - q113*q188 + q176*x%d1val2_d1val3 + q200*q88 + q206*x%d1val1_d1val3 - q208*q69 + q212*q213 - x%d1val1_d2val2_d1val3) - q88*x%d1val1_d1val2_d2val3*x%d1val2 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = -q2*(-288.0_dp*q201*q203 + 1536.0_dp*q202*x%d2val2 + 4608.0_dp*q11*q110*q63 - q0*q186*x%d2val2 + q114*q50*x%d1val2_d1val3 + q135*q216 - q136*q159*q63 + q140*(2.0_dp*x%d3val2 + q157*q57 - q203*q216 - q54*q87) + q142*q203*x%d1val2_d2val3 + q142*q221 - q153*q222 - q154*q54 - q160*q217*x%d1val2 - q160*q220 - q161*q215*x%d1val2 - q162*q221 + q164*q222 + q169*q51*q54 - q172*q54 + q216*q217 + q218*x%d2val2_d2val3*x%val + q219*x%d2val3 + q220*q50 + q45*(q104*q130 + q11*q162*x%d1val2_d1val3 - q11*q223 - q114*q218 - q129*q164 + q199*q43*x%d2val2 - q203*q223 - q219*x%d1val3 + x%d3val2_d1val3) - x%d3val2_d2val3)
   end function asinh_self
   
   function acosh_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q3 = pow2(x%d1val1)
      q4 = powm1(q1)
      q5 = q4*x%val
      q6 = -q3*q5 + x%d2val1
      q7 = powm1(pow3(sqrt(q1)))
      q8 = q7*x%val
      q9 = x%d1val1*x%d1val2
      q10 = q8*x%d1val3
      q11 = pow2(x%d1val2)
      q12 = -q11*q5 + x%d2val2
      q13 = pow2(x%d1val3)
      q14 = pow3(x%d1val1)
      q15 = q5*x%d1val1
      q16 = 3.0_dp*x%d2val1
      q17 = q0*q14
      q18 = powm1(pow2(q1))
      q19 = 3.0_dp*q18
      q20 = -q14*q4 - q15*q16 + q17*q19 + x%d3val1
      q21 = q8*x%d1val2
      q22 = 2.0_dp*q15
      q23 = q3*q4
      q24 = 2.0_dp*x%d1val2
      q25 = q0*q3
      q26 = q18*q25
      q27 = -q22*x%d1val1_d1val2 - q23*x%d1val2 + q24*q26 + x%d2val1_d1val2
      q28 = 2.0_dp*x%d1val3
      q29 = -q22*x%d1val1_d1val3 - q23*x%d1val3 + q26*q28 + x%d2val1_d1val3
      q30 = q5*x%d1val1_d1val2
      q31 = 48.0_dp*x%d2val2
      q32 = 48.0_dp*q11
      q33 = 144.0_dp*q0
      q34 = -144.0_dp + q33
      q35 = powm1(q34)
      q36 = q0*q11
      q37 = -20736.0_dp*q35*q36 + q31*x%val + q32
      q38 = q4*x%d1val1
      q39 = 0.020833333333333332_dp*q38
      q40 = -q24*q30 - q37*q39 + x%d1val1_d2val2
      q41 = q7*x%d1val3
      q42 = q8*x%d1val2_d1val3
      q43 = q0*x%d1val3
      q44 = 3.0_dp*q43*powm1(pow5(sqrt(q1)))
      q45 = q28*q5
      q46 = 2.0_dp*x%val
      q47 = 4.0_dp*q0
      q48 = -4.0_dp + q47
      q49 = powm1(q48)
      q50 = 24.0_dp*q49
      q51 = q0*q13
      q52 = 2.0_dp*q13 + q46*x%d2val3 - q50*q51
      q53 = 0.5_dp*q52
      q54 = pow3(x%d1val2)
      q55 = q5*x%d1val2
      q56 = 3.0_dp*x%d2val2
      q57 = q0*q54
      q58 = q19*q57 - q4*q54 - q55*q56 + x%d3val2
      q59 = x%d1val2*x%d1val2_d1val3
      q60 = 2.0_dp*q5
      q61 = q11*q4
      q62 = q18*q36
      q63 = x%d1val2_d1val3*x%d1val3
      q64 = q4*x%d1val2
      q65 = 3.0_dp*x%d1val1_d1val3
      q66 = q14*x%d1val3
      q67 = q66*x%val
      q68 = 8.0_dp*q18
      q69 = x%d1val1*x%d1val3
      q70 = q18*q69
      q71 = q0*x%d2val1
      q72 = pow3(x%val)
      q73 = q72*powm1(pow3(q1))
      q74 = 12.0_dp*q73
      q75 = q6*x%d1val2
      q76 = q28*x%d1val1_d1val2
      q77 = q4*x%d1val1_d1val3
      q78 = 2.0_dp*q30
      q79 = q18*q47
      q80 = q9*x%d1val1_d1val3
      q81 = q3*x%val
      q82 = q81*x%d1val3
      q83 = 6.0_dp*q18
      q84 = q83*x%d1val2
      q85 = q3*x%d1val2
      q86 = q85*x%d1val3
      q87 = 8.0_dp*q49
      q88 = q87*x%val
      q89 = 2.0_dp*x%d2val1 - q3*q88
      q90 = q4*q89
      q91 = 0.25_dp*q52
      q92 = q88*x%d1val1
      q93 = 4.0_dp*x%d1val3
      q94 = q3*q49
      q95 = powm1(pow2(q48))
      q96 = 32.0_dp*q95
      q97 = q25*q96
      q98 = -q92*x%d1val1_d1val3 - q93*q94 + q97*x%d1val3 + x%d2val1_d1val3
      q99 = 16.0_dp*q49
      q100 = q99*x%d1val1_d1val3
      q101 = 4.0_dp*q94
      q102 = pow2(x%d1val1_d1val3)
      q103 = q0*q69
      q104 = 128.0_dp*q95
      q105 = q104*x%d1val1_d1val3
      q106 = q13*q95
      q107 = 96.0_dp*q106
      q108 = q13*q3
      q109 = powm1(pow3(q48))
      q110 = q109*q72
      q111 = 512.0_dp*q110
      q112 = q100*q69 + q101*x%d2val3 + q102*q88 - q103*q105 - q107*q81 + q108*q111 + q92*x%d1val1_d2val3 - q97*x%d2val3 - x%d2val1_d2val3
      q113 = x%d1val1_d1val2*x%d1val3
      q114 = x%d2val2_d1val3*x%val
      q115 = 41472.0_dp*q35
      q116 = q0*q59
      q117 = q11*x%d1val3
      q118 = q117*q72
      q119 = q9*x%d2val3
      q120 = q28*q64
      q121 = q60*x%d1val2_d1val3
      q122 = q13*q4
      q123 = q13*q9
      q124 = q18*x%val
      q125 = q0*q19
      q126 = q0*x%d1val1
      q127 = q43*q84
      q128 = 3.0_dp*x%d1val2_d1val3
      q129 = q54*x%d1val3
      q130 = q129*x%val
      q131 = q99*x%d1val3
      q132 = x%d1val2*x%d1val2_d2val3
      q133 = 4.0_dp*q11
      q134 = q133*q49
      q135 = pow2(x%d1val2_d1val3)
      q136 = q11*x%val
      q137 = q36*q96
      q138 = q11*q13
      q139 = 2.0_dp*x%d2val2
      q140 = q4*q91
      q141 = q50*x%d1val1
      q142 = 12.0_dp*q49
      q143 = q142*x%d2val1
      q144 = x%d1val1*x%d2val3
      q145 = q142*x%d1val1
      q146 = q50*x%d1val1_d1val3
      q147 = x%d1val3*x%d2val1
      q148 = x%d2val1_d1val3*x%val
      q149 = x%d2val1*x%val
      q150 = q142*q149
      q151 = 12.0_dp*q94
      q152 = q14*x%d2val3
      q153 = q104*x%val
      q154 = 128.0_dp*q106
      q155 = 192.0_dp*q95
      q156 = q106*x%d1val1
      q157 = 96.0_dp*q95
      q158 = q157*x%d2val3
      q159 = 768.0_dp*q95
      q160 = q0*q155
      q161 = 288.0_dp*q95
      q162 = q33*q95
      q163 = q162*q3
      q164 = 768.0_dp*q110
      q165 = q110*q13
      q166 = q165*x%d1val1
      q167 = q110*q3
      q168 = x%d1val1_d1val3*x%d1val3
      q169 = 4352.0_dp*q109
      q170 = pow4(x%val)
      q171 = q170*powm1(pow4(q48))
      q172 = 18432.0_dp*q13*q171
      q173 = q87*x%d1val1_d1val2
      q174 = q99*x%d1val1_d1val2_d1val3
      q175 = q87*q9
      q176 = q88*x%d1val1_d1val2
      q177 = 4.0_dp*x%d1val2_d2val3
      q178 = q87*x%d1val2
      q179 = q95*x%val
      q180 = q80*x%d1val3
      q181 = 192.0_dp*x%val
      q182 = q181*x%d1val1_d1val2
      q183 = 64.0_dp*q95
      q184 = x%d1val1_d1val2*x%d2val3
      q185 = q0*q183
      q186 = q158*x%d1val2
      q187 = q183*x%d1val2
      q188 = q0*q187
      q189 = 1024.0_dp*x%d1val1_d1val2
      q190 = 3072.0_dp*q109
      q191 = 0.5_dp*q89
      q192 = q89*x%d1val2
      q193 = q13*q192
      q194 = 1.5_dp*q18
      q195 = q99*x%d1val1_d1val2
      q196 = q195*x%val
      q197 = 2.0_dp*x%d2val1_d1val2 - q178*q3 + q187*q25 - q196*x%d1val1
      q198 = 0.5_dp*q197
      q199 = q157*x%d1val2
      q200 = x%d1val1_d1val2_d1val3*x%d1val2
      q201 = q106*x%d1val2
      q202 = q165*x%d1val2
      q203 = x%d2val2*x%val
      q204 = q32*q49
      q205 = 4.0_dp*q203 - q0*q204 + q133
      q206 = q205*q49
      q207 = 8.0_dp*q205
      q208 = q179*q207
      q209 = 48.0_dp*q49
      q210 = q204*x%val
      q211 = 4.0_dp*q59 - q116*q209 + q118*q155 + q139*x%d1val3 - q210*x%d1val3 + q46*x%d2val2_d1val3
      q212 = q211*q49
      q213 = 2.0_dp*x%d1val1
      q214 = q59*x%d1val3
      q215 = q0*q135
      q216 = q50*x%d1val2
      q217 = x%d1val3*x%d2val2_d1val3
      q218 = q142*x%d1val2
      q219 = q218*x%d2val2
      q220 = q63*x%d2val2
      q221 = q11*x%d1val2_d2val3
      q222 = q54*x%d2val3
      q223 = q142*x%d1val2_d1val3
      unary%val = acosh(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q2*q6
      unary%d1val1_d1val2 = q2*x%d1val1_d1val2 - q8*q9
      unary%d1val1_d1val3 = -q10*x%d1val1 + q2*x%d1val1_d1val3
      unary%d2val2 = q12*q2
      unary%d1val2_d1val3 = -q10*x%d1val2 + q2*x%d1val2_d1val3
      unary%d2val3 = q2*(-q13*q5 + x%d2val3)
      unary%d3val1 = q2*q20
      unary%d2val1_d1val2 = q2*q27 - q21*q6
      unary%d2val1_d1val3 = -q10*q6 + q2*q29
      unary%d1val1_d2val2 = q2*q40
      unary%d1val1_d1val2_d1val3 = -q10*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3 - q21*x%d1val1_d1val3 - q41*q9 - q42*x%d1val1 + q44*q9
      unary%d1val1_d2val3 = q2*(-q38*q53 - q45*x%d1val1_d1val3 + x%d1val1_d2val3)
      unary%d3val2 = q2*q58
      unary%d2val2_d1val3 = -q10*q12 + q2*(q28*q62 - q59*q60 - q61*x%d1val3 + x%d2val2_d1val3)
      unary%d1val2_d2val3 = q2*(-q53*q64 - q60*q63 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = -q10*q20 + q2*(-3.0_dp*q15*x%d2val1_d1val3 + 6.0_dp*q70*q71 + 9.0_dp*q26*x%d1val1_d1val3 - q16*q38*x%d1val3 - q23*q65 - q5*q65*x%d2val1 - q66*q74 + q67*q68 + x%d3val1_d1val3)
      unary%d2val1_d1val2_d1val3 = -q10*q27 + q2*(-2.0_dp*q77*q9 - 8.0_dp*q73*q86 + 2.0_dp*q26*x%d1val2_d1val3 - q22*x%d1val1_d1val2_d1val3 - q23*x%d1val2_d1val3 - q38*q76 + q47*q70*x%d1val1_d1val2 - q78*x%d1val1_d1val3 + q79*q80 + q82*q84 + x%d2val1_d1val2_d1val3) - q21*q29 - q41*q75 - q42*q6 + q44*q75
      unary%d2val1_d2val3 = -q2*(q112 + q45*q98 + q90*q91)
      unary%d1val1_d2val2_d1val3 = -q10*q40 + q2*(-0.020833333333333332_dp*q37*q77 - 2.0_dp*q55*x%d1val1_d1val2_d1val3 + 0.041666666666666664_dp*q37*q70*x%val + q113*q79*x%d1val2 - q39*(48.0_dp*q114 + 5971968.0_dp*q118*powm1(pow2(q34)) + 96.0_dp*q59 - q115*q116 - q115*q117*x%val + q31*x%d1val3) - q64*q76 - q78*x%d1val2_d1val3 + x%d1val1_d2val2_d1val3)
      unary%d1val1_d1val2_d2val3 = q2*(-15.0_dp*q123*q73 - 2.0_dp*q38*q63 + 9.0_dp*q123*q124 + q119*q125 - q119*q4 - q120*x%d1val1_d1val3 - q121*x%d1val1_d1val3 - q122*x%d1val1_d1val2 + q126*q63*q83 + q127*x%d1val1_d1val3 - q15*x%d1val2_d2val3 + q19*q51*x%d1val1_d1val2 - q30*x%d2val3 - q45*x%d1val1_d1val2_d1val3 - q55*x%d1val1_d2val3 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = -q10*q58 + q2*(-3.0_dp*q55*x%d2val2_d1val3 + 9.0_dp*q62*x%d1val2_d1val3 + q127*x%d2val2 - q128*q5*x%d2val2 - q128*q61 - q129*q74 + q130*q68 - q56*q64*x%d1val3 + x%d3val2_d1val3)
      unary%d2val2_d2val3 = -q2*(-q104*q116*x%d1val3 - q107*q136 + q111*q138 + q131*q59 + q132*q88 + q134*x%d2val3 + q135*q88 - q137*x%d2val3 + q140*(-q11*q88 + q139) + q45*(-q134*x%d1val3 + q137*x%d1val3 - q59*q88 + x%d2val2_d1val3) - x%d2val2_d2val3)
      unary%d3val1_d2val3 = -q2*(-288.0_dp*q149*q156 + 1536.0_dp*q166*x%d2val1 + 4608.0_dp*q167*q168 - q102*q126*q161 + q102*q141 - q103*q155*x%d2val1_d1val3 - q126*q158*x%d2val1 + q13*q169*q17 - q14*q154 - q14*q172 + q140*(2.0_dp*x%d3val1 - q14*q87 - q141*q149 + q157*q17) + q141*x%d1val3*x%d2val1_d1val3 + q143*q144 + q145*x%d2val1_d2val3*x%val + q146*q147 + q146*q148 - q147*q160*x%d1val1_d1val3 + q150*x%d1val1_d2val3 + q151*x%d1val1_d2val3 - q152*q153 + q152*q164 - q159*q82*x%d1val1_d1val3 - q163*x%d1val1_d2val3 + q45*(q104*q67 - q143*q69 - q145*q148 - q150*x%d1val1_d1val3 - q151*x%d1val1_d1val3 + q157*q69*q71 + q163*x%d1val1_d1val3 - q164*q66 + x%d3val1_d1val3) - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(-1024.0_dp*q167*q63 - 2048.0_dp*q110*q180 - 7.5_dp*q193*q73 + 12288.0_dp*q108*q171*x%d1val2 + 384.0_dp*q179*q180 + 4.5_dp*q124*q193 + q0*q105*q113 + q0*q192*q194*x%d2val3 - q100*q113 - q100*x%d1val1*x%d1val2_d1val3 - q100*x%d1val1_d1val2_d1val3*x%val - q102*q178 + q102*q188 + q103*q104*x%d1val1_d1val2_d1val3 + q105*q126*x%d1val2_d1val3 + q107*q85 - q111*q85*x%d2val3 + q112*q55 - q120*q98 - q121*q98 - q122*q198 + q125*q63*q89 + q126*q183*q184 + q127*q98 - q144*q173 + q155*q63*q81 + q156*q182 - q166*q189 - q174*q69 - q175*x%d1val1_d2val3 - q176*x%d1val1_d2val3 - q177*q94 + q185*q9*x%d1val1_d2val3 + q186*q81 - q190*q51*q85 - q191*q5*x%d1val2_d2val3 - q191*q64*x%d2val3 + q194*q197*q51 - q198*q5*x%d2val3 - q45*(-q101*x%d1val2_d1val3 + q103*q183*x%d1val1_d1val2 - q111*q86 - q173*q69 - q175*x%d1val1_d1val3 - q176*x%d1val1_d1val3 + q185*q80 + q199*q82 - q92*x%d1val1_d1val2_d1val3 + q97*x%d1val2_d1val3 + x%d2val1_d1val2_d1val3) - q63*q90 - q92*x%d1val1_d1val2_d2val3 + q97*x%d1val2_d2val3 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(-128.0_dp*q109*q205*q51*x%d1val1 - 4.0_dp*q212*x%d1val1_d1val3 + 16.0_dp*q168*q179*q205 + q0*q104*q63*x%d1val1_d1val2 + q104*q200*q43 - q131*q200 - q140*(2.0_dp*x%d1val1_d2val2 - q196*x%d1val2 - q206*q213) + q144*q208 + q156*q207 - q174*x%d1val2_d1val3*x%val - q176*x%d1val2_d2val3 - q178*q184 + q182*q201 + q184*q188 - q189*q202 - q195*q63 - q206*x%d1val1_d2val3 + q211*q69*q96*x%val - q213*q49*(4.0_dp*q135 + 960.0_dp*q106*q36 - q0*q132*q209 + q11*q155*q72*x%d2val3 - q13*q204 - q138*q170*q190 + q139*x%d2val3 + q159*q214*q72 + q177*x%d1val2 - q181*q214*q49 - q209*q215 - q210*x%d2val3 + q46*x%d2val2_d2val3 + q93*x%d2val2_d1val3) + q45*(q113*q178 - q113*q188 + q176*x%d1val2_d1val3 + q200*q88 + q206*x%d1val1_d1val3 - q208*q69 + q212*q213 - x%d1val1_d2val2_d1val3) - q88*x%d1val1_d1val2_d2val3*x%d1val2 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = -q2*(-288.0_dp*q201*q203 + 1536.0_dp*q202*x%d2val2 + 4608.0_dp*q11*q110*q63 - q0*q186*x%d2val2 + q114*q50*x%d1val2_d1val3 + q135*q216 - q136*q159*q63 + q140*(2.0_dp*x%d3val2 + q157*q57 - q203*q216 - q54*q87) + q142*q203*x%d1val2_d2val3 + q142*q221 - q153*q222 - q154*q54 - q160*q217*x%d1val2 - q160*q220 - q161*q215*x%d1val2 - q162*q221 + q164*q222 + q169*q51*q54 - q172*q54 + q216*q217 + q218*x%d2val2_d2val3*x%val + q219*x%d2val3 + q220*q50 + q45*(q104*q130 + q11*q162*x%d1val2_d1val3 - q11*q223 - q114*q218 - q129*q164 + q199*q43*x%d2val2 - q203*q223 - q219*x%d1val3 + x%d3val2_d1val3) - x%d3val2_d2val3)
   end function acosh_self
   
   function atanh_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q2 = powm1(q1)
      q3 = 2.0_dp*x%d2val1
      q4 = pow2(x%d1val1)
      q5 = -144.0_dp + 144.0_dp*q0
      q6 = powm1(q5)
      q7 = 576.0_dp*q6
      q8 = q7*x%val
      q9 = -q3 + q4*q8
      q10 = 72.0_dp*q6
      q11 = powm1(pow2(q1))
      q12 = 2.0_dp*x%val
      q13 = q11*q12
      q14 = x%d1val1*x%d1val2
      q15 = q13*x%d1val3
      q16 = 2.0_dp*x%d2val2
      q17 = pow2(x%d1val2)
      q18 = -q16 + q17*q8
      q19 = pow2(x%d1val3)
      q20 = -4.0_dp + 4.0_dp*q0
      q21 = powm1(q20)
      q22 = q19*q21
      q23 = 8.0_dp*q22
      q24 = 4.0_dp*q21
      q25 = 2.0_dp*x%d3val1
      q26 = pow3(x%d1val1)
      q27 = q6*x%d1val1
      q28 = q27*x%val
      q29 = 1728.0_dp*x%d2val1
      q30 = q0*q26
      q31 = powm1(pow2(q5))
      q32 = 331776.0_dp*q31
      q33 = -q25 + q26*q7 + q28*q29 - q30*q32
      q34 = q31*x%d1val2
      q35 = 20736.0_dp*x%val
      q36 = q35*q9
      q37 = 2.0_dp*x%d2val1_d1val2
      q38 = 1152.0_dp*q28
      q39 = q4*q7
      q40 = q0*q4
      q41 = 165888.0_dp*q40
      q42 = -q34*q41 - q37 + q38*x%d1val1_d1val2 + q39*x%d1val2
      q43 = q31*x%d1val3
      q44 = -2.0_dp*x%d2val1_d1val3 + q38*x%d1val1_d1val3 + q39*x%d1val3 - q41*q43
      q45 = 2.0_dp*x%d1val1_d2val2
      q46 = q6*x%d1val2
      q47 = 1152.0_dp*x%d1val1_d1val2
      q48 = q47*x%val
      q49 = 48.0_dp*x%d2val2
      q50 = q49*x%val
      q51 = q17*q6
      q52 = -27648.0_dp*q0*q51 + 48.0_dp*q17 + q50
      q53 = 12.0_dp*q27
      q54 = -q45 + q46*q48 + q52*q53
      q55 = q14*x%d1val3
      q56 = x%d1val1*x%d1val2_d1val3
      q57 = x%d1val1_d1val3*x%d1val2
      q58 = 16.0_dp*q21
      q59 = q58*x%val
      q60 = q59*x%d1val3
      q61 = -32.0_dp*q0*q22 + 2.0_dp*q19 + q12*x%d2val3
      q62 = q24*q61
      q63 = 2.0_dp*x%d3val2
      q64 = pow3(x%d1val2)
      q65 = x%d2val2*x%val
      q66 = 1728.0_dp*q46
      q67 = q0*q64
      q68 = -q32*q67 - q63 + q64*q7 + q65*q66
      q69 = q35*q43
      q70 = 2.0_dp*x%d2val2_d1val3
      q71 = x%d1val2*x%d1val2_d1val3
      q72 = q6*q71
      q73 = 1152.0_dp*x%val
      q74 = q17*x%d1val3
      q75 = q0*q31
      q76 = q27*x%d1val3
      q77 = q6*x%d1val1_d1val3
      q78 = 1728.0_dp*q77
      q79 = x%d2val1*x%val
      q80 = q26*x%val
      q81 = 829440.0_dp*q43
      q82 = q0*x%d1val1
      q83 = q4*x%d1val1_d1val3
      q84 = 995328.0_dp*q75
      q85 = pow3(x%val)
      q86 = q26*q85
      q87 = x%d1val3*powm1(pow3(q5))
      q88 = 191102976.0_dp*q87
      q89 = q34*x%d1val3
      q90 = q31*x%d1val2_d1val3
      q91 = q87*x%d1val2
      q92 = q0*q32
      q93 = x%d1val1*x%d1val1_d1val2
      q94 = q93*x%d1val3
      q95 = q14*x%d1val1_d1val3
      q96 = q4*x%val
      q97 = 497664.0_dp*q89
      q98 = q4*q85
      q99 = q3 - q4*q59
      q100 = 2.0_dp*q21
      q101 = q100*q61
      q102 = q59*x%d1val1
      q103 = 8.0_dp*q21
      q104 = q103*q4
      q105 = powm1(pow2(q20))
      q106 = 64.0_dp*q105
      q107 = q106*q40
      q108 = -q102*x%d1val1_d1val3 - q104*x%d1val3 + q107*x%d1val3 + x%d2val1_d1val3
      q109 = x%d1val1_d1val3*x%d1val3
      q110 = 32.0_dp*q21
      q111 = q110*x%d1val1
      q112 = q4*x%d2val3
      q113 = pow2(x%d1val1_d1val3)
      q114 = 256.0_dp*q105
      q115 = q114*q82
      q116 = q105*q19
      q117 = 192.0_dp*q116
      q118 = q0*q106
      q119 = powm1(pow3(q20))
      q120 = q119*q19
      q121 = 1024.0_dp*q120
      q122 = q102*x%d1val1_d2val3 + q103*q112 + q109*q111 - q109*q115 - q112*q118 + q113*q59 - q117*q96 + q121*q98 - x%d2val1_d2val3
      q123 = x%d1val2*x%d1val3
      q124 = x%d1val1*x%val
      q125 = q49*x%d1val3
      q126 = x%d2val2_d1val3*x%val
      q127 = 48.0_dp*q126
      q128 = x%d1val3*x%val
      q129 = q74*q85
      q130 = q103*x%d2val3
      q131 = q58*x%d1val3
      q132 = q103*x%val
      q133 = x%d1val1_d1val2*x%val
      q134 = q58*x%d1val2
      q135 = q59*x%d1val2_d1val3
      q136 = q132*x%d1val2
      q137 = 384.0_dp*x%val
      q138 = 128.0_dp*q0
      q139 = q105*q138
      q140 = q139*q14
      q141 = q0*q114
      q142 = q141*x%d1val3
      q143 = q116*x%d1val1_d1val2
      q144 = q120*q85
      q145 = 3072.0_dp*q144
      q146 = x%d1val3*x%d2val2
      q147 = 1728.0_dp*x%d1val2_d1val3
      q148 = q64*x%val
      q149 = q17*x%d1val2_d1val3
      q150 = q64*q85
      q151 = q110*x%d1val3
      q152 = x%d1val2*x%d1val2_d2val3
      q153 = pow2(x%d1val2_d1val3)
      q154 = q117*x%val
      q155 = q0*q105
      q156 = 64.0_dp*q17
      q157 = q156*x%d2val3
      q158 = q17*q85
      q159 = 48.0_dp*q21
      q160 = q159*x%d1val1
      q161 = x%d1val3*x%d2val1_d1val3
      q162 = 24.0_dp*q21
      q163 = q162*x%d1val1
      q164 = q163*x%d2val1
      q165 = q163*x%val
      q166 = q109*x%d2val1
      q167 = x%d1val1_d1val3*x%val
      q168 = q162*q79
      q169 = q113*x%d1val1
      q170 = q4*x%d1val1_d2val3
      q171 = 320.0_dp*q105
      q172 = q171*q80
      q173 = 320.0_dp*q116
      q174 = 384.0_dp*q155
      q175 = q116*x%d1val1
      q176 = x%d1val1*x%d2val1
      q177 = q155*q176
      q178 = q105*q128
      q179 = 768.0_dp*q155
      q180 = 2048.0_dp*q119
      q181 = q180*q86
      q182 = q119*q85
      q183 = 11264.0_dp*q120
      q184 = pow4(x%val)
      q185 = q184*q19*powm1(pow4(q20))
      q186 = 49152.0_dp*q185
      q187 = 192.0_dp*x%d1val3
      q188 = q93*x%d2val3
      q189 = x%d1val1_d1val2_d1val3*x%d1val3
      q190 = q56*x%d1val1_d1val3
      q191 = q134*x%d1val1
      q192 = q109*x%d1val1_d1val2
      q193 = q59*x%d1val1_d1val2
      q194 = q110*x%d1val1_d1val2_d1val3
      q195 = q105*x%d1val2
      q196 = x%d1val2_d1val3*x%d1val3
      q197 = q4*x%d1val2
      q198 = q138*q195
      q199 = 2048.0_dp*q144
      q200 = q24*q99
      q201 = x%d1val2*x%d2val3
      q202 = q196*q99
      q203 = q99*x%d1val2
      q204 = q134*x%d1val3
      q205 = 128.0_dp*q195*q40 - q111*q133 - q134*q4 + q37
      q206 = q205*q24
      q207 = x%d2val3*x%val
      q208 = q0*q116
      q209 = x%d1val1_d1val2*x%d2val3
      q210 = x%d1val2_d1val3*x%val
      q211 = q59*x%d1val2
      q212 = q0*q21
      q213 = 4.0_dp*q17 + 4.0_dp*q65 - q156*q212
      q214 = q100*q213
      q215 = 16.0_dp*q213
      q216 = q215*x%d1val1
      q217 = 64.0_dp*q212
      q218 = q21*x%val
      q219 = -64.0_dp*q218*q74 + 4.0_dp*q71 + q114*q129 + q16*x%d1val3 - q217*q71 + q70*x%val
      q220 = q24*x%d1val1
      q221 = x%d1val3*x%d2val2_d1val3
      q222 = q221*x%d1val2
      q223 = q201*x%d2val2
      q224 = q162*x%d1val2
      q225 = q21*x%d1val2_d1val3
      q226 = q162*q65
      q227 = q153*x%d1val2
      q228 = q162*q17
      q229 = 192.0_dp*q155
      q230 = q150*q180
      unary%val = atanh(x%val)
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q10*q9
      unary%d1val1_d1val2 = q13*q14 + q2*x%d1val1_d1val2
      unary%d1val1_d1val3 = q15*x%d1val1 + q2*x%d1val1_d1val3
      unary%d2val2 = q10*q18
      unary%d1val2_d1val3 = q15*x%d1val2 + q2*x%d1val2_d1val3
      unary%d2val3 = q24*(q23*x%val - x%d2val3)
      unary%d3val1 = q10*q33
      unary%d2val1_d1val2 = q10*q42 - q34*q36
      unary%d2val1_d1val3 = q10*q44 - q36*q43
      unary%d1val1_d2val2 = q10*q54
      unary%d1val1_d1val2_d1val3 = 2.0_dp*q11*q55 + 8.0_dp*q0*q55*powm1(pow3(q1)) + q13*q56 + q13*q57 + q15*x%d1val1_d1val2 + q2*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q24*(q60*x%d1val1_d1val3 + q62*x%d1val1 - x%d1val1_d2val3)
      unary%d3val2 = q10*q68
      unary%d2val2_d1val3 = q10*(-165888.0_dp*q74*q75 + q7*q74 - q70 + q72*q73) - q18*q69
      unary%d1val2_d2val3 = q24*(q60*x%d1val2_d1val3 + q62*x%d1val2 - x%d1val2_d2val3)
      unary%d3val1_d1val3 = q10*(-2.0_dp*x%d3val1_d1val3 - 497664.0_dp*q43*q82*x%d2val1 + 1728.0_dp*q28*x%d2val1_d1val3 + q29*q76 + q4*q78 + q78*q79 - q80*q81 - q83*q84 + q86*q88) - q33*q69
      unary%d2val1_d1val2_d1val3 = -20736.0_dp*q89*q9 + 11943936.0_dp*q0*q9*q91 + q10*(-2.0_dp*x%d2val1_d1val2_d1val3 + 1152.0_dp*q27*q57 + 95551488.0_dp*q91*q98 + q38*x%d1val1_d1val2_d1val3 + q39*x%d1val2_d1val3 - q41*q90 + q47*q76 + q48*q77 - q92*q94 - q92*q95 - q96*q97) - q34*q35*q44 - q36*q90 - q42*q69
      unary%d2val1_d2val3 = q24*(q101*q99 + q108*q60 + q122)
      unary%d1val1_d2val2_d1val3 = q10*(-2.0_dp*x%d1val1_d2val2_d1val3 - 3456.0_dp*q124*q43*q52 + 12.0_dp*q52*q77 - q123*q92*x%d1val1_d1val2 + q46*q47*x%d1val3 + q46*q73*x%d1val1_d1val2_d1val3 + q48*q6*x%d1val2_d1val3 + q53*(-55296.0_dp*q0*q72 - 55296.0_dp*q128*q51 + 7962624.0_dp*q129*q31 + 96.0_dp*q71 + q125 + q127)) - q54*q69
      unary%d1val1_d1val2_d2val3 = q24*(q109*q134 - q116*q137*q14 + q130*q133 + q130*q14 + q131*q56 + q132*x%d1val1*x%d1val2_d2val3 + q135*x%d1val1_d1val3 + q136*x%d1val1_d2val3 - q138*q143 + q14*q145 - q140*x%d2val3 - q142*q56 - q142*q57 + q23*x%d1val1_d1val2 + q60*x%d1val1_d1val2_d1val3 - x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q10*(-2.0_dp*x%d3val2_d1val3 - q0*q97*x%d2val2 + q126*q66 + q146*q66 + q147*q51 + q147*q6*q65 - q148*q81 - q149*q84 + q150*q88) - q68*q69
      unary%d2val2_d2val3 = q24*(q101*(q16 - q17*q59) + q121*q158 + q130*q17 - q142*q71 + q151*q71 + q152*q59 + q153*q59 - q154*q17 - q155*q157 + q60*(-q103*q74 + q118*q74 - q59*q71 + x%d2val2_d1val3) - x%d2val2_d2val3)
      unary%d3val1_d2val3 = q24*(-192.0_dp*q177*x%d2val3 - 1920.0_dp*q178*q83 - 576.0_dp*q175*q79 + 12288.0_dp*q182*q83*x%d1val3 + q101*(q114*q30 - q160*q79 + q25 - q26*q58) + q145*q176 + q159*q166 + q159*q167*x%d2val1_d1val3 + q159*q169 + q160*q161 - q161*q174*x%d1val1 + q162*q170 + q164*x%d2val3 + q165*x%d2val1_d2val3 - q166*q174 + q168*x%d1val1_d2val3 - q169*q179 - q170*q174 - q172*x%d2val3 - q173*q26 + q181*x%d2val3 + q183*q30 - q186*q26 + q60*(-q162*q83 - q164*x%d1val3 - q165*x%d2val1_d1val3 - q168*x%d1val1_d1val3 + q172*x%d1val3 + q174*q83 + q177*q187 - q181*x%d1val3 + x%d3val1_d1val3) - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q24*(-192.0_dp*q112*q195*x%val - 24576.0_dp*q185*q197 - 384.0_dp*q105*q196*q96 - 384.0_dp*q124*q143 - 64.0_dp*q205*q208 - 768.0_dp*q105*q167*q55 + 1024.0_dp*q112*q182*x%d1val2 + 1536.0_dp*q144*q203 + 4096.0_dp*q182*q55*x%d1val1_d1val3 + 6144.0_dp*q120*q40*x%d1val2 + q102*x%d1val1_d1val2_d2val3 + q103*q202 + q104*x%d1val2_d2val3 - q107*x%d1val2_d2val3 - q108*q123*q141 + q108*q135 + q108*q204 + q110*q190 + q110*q192 + q111*q189 + q113*q134 - q113*q198 - q115*q189 - q117*q197 - q118*q201*q99 - q122*q136 - q139*q188 - q139*q202 - q140*x%d1val1_d2val3 - q141*q190 - q141*q192 - q154*q203 + q167*q194 + q180*q196*q98 + q188*q58 + q19*q206 + q191*x%d1val1_d2val3 + q193*x%d1val1_d2val3 + q199*q93 + q200*q201 + q200*x%d1val2_d2val3*x%val + q206*q207 + q60*(-1024.0_dp*q119*q123*q98 - q102*x%d1val1_d1val2_d1val3 - q104*x%d1val2_d1val3 + q107*x%d1val2_d1val3 - q131*q93 + q139*q94 + q139*q95 + q187*q195*q96 - q191*x%d1val1_d1val3 - q193*x%d1val1_d1val3 + x%d2val1_d1val2_d1val3) - x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q24*(-32.0_dp*q105*q109*q213*x%val + 256.0_dp*q120*q213*q82 + q101*(-q110*q133*x%d1val2 - q213*q220 + q45) + q103*q219*x%d1val1_d1val3 - q105*q207*q216 - q106*q128*q219*x%d1val1 + q134*q209 - q137*q143*x%d1val2 - q141*q189*x%d1val2 - q141*q196*x%d1val1_d1val2 + q151*x%d1val1_d1val2*x%d1val2_d1val3 + q151*x%d1val1_d1val2_d1val3*x%d1val2 - q175*q215 + q193*x%d1val2_d2val3 + q194*q210 - q198*q209 + q199*x%d1val1_d1val2*x%d1val2 + q211*x%d1val1_d1val2_d2val3 + q214*x%d1val1_d2val3 + q220*(-256.0_dp*q128*q21*q71 - 4096.0_dp*q120*q17*q184 + 1024.0_dp*q105*q71*q85*x%d1val3 + 1280.0_dp*q17*q208 + 4.0_dp*q152 + 4.0_dp*q153 + 4.0_dp*q221 + q114*q158*x%d2val3 + q12*x%d2val2_d2val3 - q152*q217 - q153*q217 - q156*q22 - q157*q218 + q16*x%d2val3) - q60*(q135*x%d1val1_d1val2 - q178*q216 - q198*x%d1val1_d1val2*x%d1val3 + q204*x%d1val1_d1val2 + q211*x%d1val1_d1val2_d1val3 + q214*x%d1val1_d1val3 + q219*q220 - x%d1val1_d2val2_d1val3) - x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q24*(-1920.0_dp*q105*q210*q74 - 576.0_dp*q116*q65*x%d1val2 + 12288.0_dp*q119*q129*x%d1val2_d1val3 + q101*(q114*q67 - q21*q50*x%d1val2 - q58*q64 + q63) + q125*q225 + q127*q225 + q145*x%d1val2*x%d2val2 - q146*q174*x%d1val2_d1val3 - q148*q171*x%d2val3 + q159*q222 + q159*q227 + q162*q223 - q17*q174*x%d1val2_d2val3 - q173*q64 - q174*q222 - q179*q227 + q183*q67 - q186*q64 - q223*q229 + q224*x%d2val2_d2val3*x%val + q226*x%d1val2_d2val3 + q228*x%d1val2_d2val3 + q230*x%d2val3 + q60*(-q126*q224 + q128*q171*q64 - q146*q224 + q146*q229*x%d1val2 + q149*q174 - q226*x%d1val2_d1val3 - q228*x%d1val2_d1val3 - q230*x%d1val3 + x%d3val2_d1val3) - x%d3val2_d2val3)
   end function atanh_self
   
   function sqrt_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q1 = powm1(q0)
      q2 = 0.5_dp*q1
      q3 = pow2(x%d1val1)
      q4 = powm1(x%val)
      q5 = 16.0_dp*q4
      q6 = 32.0_dp*x%d2val1 - q3*q5
      q7 = 0.015625_dp*q1
      q8 = powm1(pow3(sqrt(x%val)))
      q9 = q8*x%d1val2
      q10 = 0.25_dp*x%d1val1
      q11 = q8*x%d1val3
      q12 = pow2(x%d1val2)
      q13 = 32.0_dp*x%d2val2 - q12*q5
      q14 = 0.25_dp*q9
      q15 = 2.0_dp*x%d2val3
      q16 = pow2(x%d1val3)
      q17 = q16*q4
      q18 = q4*x%d1val1
      q19 = 384.0_dp*q18
      q20 = pow3(x%d1val1)
      q21 = powm1(pow2(x%val))
      q22 = 192.0_dp*q21
      q23 = 256.0_dp*x%d3val1 - q19*x%d2val1 + q20*q22
      q24 = 0.001953125_dp*q1
      q25 = 0.0078125_dp*q6
      q26 = 32.0_dp*q18
      q27 = q21*q3
      q28 = q27*x%d1val2
      q29 = 16.0_dp*q28 + 32.0_dp*x%d2val1_d1val2 - q26*x%d1val1_d1val2
      q30 = q21*x%d1val3
      q31 = q3*q30
      q32 = 16.0_dp*q31 + 32.0_dp*x%d2val1_d1val3 - q26*x%d1val1_d1val3
      q33 = q4*x%d1val2
      q34 = 64.0_dp*x%d1val1_d1val2
      q35 = q12*q4
      q36 = -12.0_dp*q35 + 8.0_dp*x%d2val2
      q37 = 4.0_dp*q18
      q38 = 64.0_dp*x%d1val1_d2val2 - q33*q34 - q36*q37
      q39 = 0.0078125_dp*q1
      q40 = q8*x%d1val2_d1val3
      q41 = x%d1val2*x%d1val3
      q42 = q41*powm1(pow5(sqrt(x%val)))
      q43 = q4*x%d1val1_d1val3
      q44 = 4.0_dp*x%d1val3
      q45 = -1.5_dp*q17 + x%d2val3
      q46 = 2.0_dp*q4
      q47 = q45*q46
      q48 = 0.125_dp*q1
      q49 = 384.0_dp*q33
      q50 = pow3(x%d1val2)
      q51 = 256.0_dp*x%d3val2 + q22*q50 - q49*x%d2val2
      q52 = 0.0078125_dp*q11
      q53 = q33*x%d1val2_d1val3
      q54 = q12*q30
      q55 = q4*q44
      q56 = 0.0009765625_dp*q11
      q57 = 384.0_dp*x%d2val1
      q58 = q30*x%d1val1
      q59 = q27*x%d1val1_d1val3
      q60 = powm1(pow3(x%val))
      q61 = q20*q60
      q62 = 384.0_dp*x%d1val3
      q63 = 32.0_dp*q43
      q64 = q58*x%d1val1_d1val2
      q65 = q21*x%d1val2
      q66 = q65*x%d1val1
      q67 = q66*x%d1val1_d1val3
      q68 = q27*x%d1val2_d1val3
      q69 = q3*q60
      q70 = 8.0_dp*q4
      q71 = x%d1val1*x%d1val1_d2val3
      q72 = pow2(x%d1val1_d1val3)
      q73 = q58*x%d1val1_d1val3
      q74 = q27*x%d2val3
      q75 = q16*q60
      q76 = 8.0_dp*q75
      q77 = q46*x%d1val1
      q78 = 2.0_dp*x%d2val1_d1val3 + q31 - q77*x%d1val1_d1val3
      q79 = 4.0_dp*x%d2val1 - q3*q46
      q80 = q4*q45
      q81 = 0.0625_dp*q1
      q82 = q4*x%d1val2_d1val3
      q83 = q33*x%d1val1_d1val2_d1val3
      q84 = q30*x%d1val2
      q85 = 4.0_dp*q36
      q86 = q5*x%d1val2_d2val3
      q87 = x%d1val1_d1val2*x%d2val3
      q88 = 32.0_dp*x%d1val1_d1val2_d1val3
      q89 = q4*x%d1val3
      q90 = q5*x%d1val2
      q91 = 24.0_dp*x%d2val3
      q92 = 48.0_dp*q84
      q93 = q16*q21
      q94 = q75*x%d1val1
      q95 = 384.0_dp*x%d2val2
      q96 = q12*q21
      q97 = q96*x%d1val2_d1val3
      q98 = q50*q60
      q99 = q70*x%d1val2_d2val3
      q100 = pow2(x%d1val2_d1val3)
      q101 = x%d1val2*x%d1val2_d1val3
      q102 = q101*q30
      q103 = q96*x%d2val3
      q104 = 2.0_dp*x%d2val2_d1val3
      q105 = 4.0_dp*x%d2val2
      q106 = 48.0_dp*x%d2val1_d1val3
      q107 = q4*x%d1val1_d2val3
      q108 = 24.0_dp*x%d2val1
      q109 = q21*x%d1val1
      q110 = q109*x%d2val3
      q111 = x%d1val1_d1val3*x%d2val1
      q112 = 48.0_dp*q30
      q113 = q21*q72
      q114 = q69*x%d1val3
      q115 = q16*powm1(pow4(x%val))
      q116 = 72.0_dp*q115
      q117 = 6.0_dp*q21
      q118 = 6.0_dp*q4
      q119 = 6.0_dp*x%d1val3
      q120 = 0.03125_dp*q1
      q121 = 128.0_dp*x%d1val1_d1val2_d1val3
      q122 = 128.0_dp*x%d1val2_d1val3
      q123 = 64.0_dp*x%d1val2
      q124 = x%d1val1_d1val2*x%d1val1_d1val3
      q125 = q41*q60
      q126 = q30*x%d1val2_d1val3
      q127 = q75*x%d1val2
      q128 = 2.0_dp*q28 + 4.0_dp*x%d2val1_d1val2 - q37*x%d1val1_d1val2
      q129 = 2.0_dp*q3
      q130 = x%d1val1_d1val2*x%d1val2
      q131 = -6.0_dp*q35 + q105
      q132 = 4.0_dp*q30
      q133 = 3.0_dp*q54 - q101*q118 + q104
      q134 = q118*x%d1val2
      q135 = 0.5_dp*q131
      q136 = 48.0_dp*x%d2val2_d1val3
      q137 = x%d1val2_d1val3*x%d2val2
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q6*q7
      unary%d1val1_d1val2 = -q10*q9 + q2*x%d1val1_d1val2
      unary%d1val1_d1val3 = -q10*q11 + q2*x%d1val1_d1val3
      unary%d2val2 = q13*q7
      unary%d1val2_d1val3 = -q14*x%d1val3 + q2*x%d1val2_d1val3
      unary%d2val3 = 0.25_dp*q1*(q15 - q17)
      unary%d3val1 = q23*q24
      unary%d2val1_d1val2 = -q25*q9 + q29*q7
      unary%d2val1_d1val3 = -q11*q25 + q32*q7
      unary%d1val1_d2val2 = q38*q39
      unary%d1val1_d1val2_d1val3 = -0.25_dp*q11*x%d1val1_d1val2 + 0.375_dp*q42*x%d1val1 - q10*q40 - q14*x%d1val1_d1val3 + q2*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q48*(4.0_dp*x%d1val1_d2val3 - q43*q44 - q47*x%d1val1)
      unary%d3val2 = q24*q51
      unary%d2val2_d1val3 = -q13*q52 + q7*(-32.0_dp*q53 + 16.0_dp*q54 + 32.0_dp*x%d2val2_d1val3)
      unary%d1val2_d2val3 = q48*(4.0_dp*x%d1val2_d2val3 - q47*x%d1val2 - q55*x%d1val2_d1val3)
      unary%d3val1_d1val3 = -q23*q56 + q24*(256.0_dp*x%d3val1_d1val3 + 576.0_dp*q59 - q19*x%d2val1_d1val3 - q43*q57 + q57*q58 - q61*q62)
      unary%d2val1_d1val2_d1val3 = -0.0078125_dp*q32*q9 + 0.01171875_dp*q42*q6 - q25*q40 - q29*q52 + q7*(-32.0_dp*q41*q69 + 16.0_dp*q68 + 32.0_dp*q64 + 32.0_dp*q67 + 32.0_dp*x%d2val1_d1val2_d1val3 - q26*x%d1val1_d1val2_d1val3 - q63*x%d1val1_d1val2)
      unary%d2val1_d2val3 = -q81*(-16.0_dp*q73 - 4.0_dp*q74 - 8.0_dp*x%d2val1_d2val3 + q3*q76 + q55*q78 + q70*q71 + q70*q72 + q79*q80)
      unary%d1val1_d2val2_d1val3 = -0.00390625_dp*q11*q38 + q39*(-64.0_dp*q83 + 64.0_dp*x%d1val1_d2val2_d1val3 - q34*q82 + q34*q84 - q37*(-24.0_dp*q53 + 12.0_dp*q54 + 8.0_dp*x%d2val2_d1val3) - q43*q85 + q58*q85)
      unary%d1val1_d1val2_d2val3 = q7*(-60.0_dp*q94*x%d1val2 + 24.0_dp*q93*x%d1val1_d1val2 + 32.0_dp*x%d1val1_d1val2_d2val3 + 48.0_dp*q58*x%d1val2_d1val3 - q5*q87 - q63*x%d1val2_d1val3 + q66*q91 - q86*x%d1val1 - q88*q89 - q90*x%d1val1_d2val3 + q92*x%d1val1_d1val3)
      unary%d3val2_d1val3 = q24*(256.0_dp*x%d3val2_d1val3 + 576.0_dp*q97 - q49*x%d2val2_d1val3 - q62*q98 - q82*q95 + q84*q95) - q51*q56
      unary%d2val2_d2val3 = -q81*(-16.0_dp*q102 - 4.0_dp*q103 - 8.0_dp*x%d2val2_d2val3 + q100*q70 + q12*q76 + q55*(-q101*q46 + q104 + q54) + q80*(-2.0_dp*q35 + q105) + q99*x%d1val2)
      unary%d3val1_d2val3 = q120*(-144.0_dp*q114*x%d1val1_d1val3 - 24.0_dp*q18*x%d2val1_d2val3 - 48.0_dp*q94*x%d2val1 + 16.0_dp*x%d3val1_d2val3 + 36.0_dp*q27*x%d1val1_d2val3 + 72.0_dp*q113*x%d1val1 - q106*q43 + q106*q58 - q107*q108 + q108*q110 + q111*q112 + q116*q20 - q55*(-6.0_dp*q18*x%d2val1_d1val3 + 4.0_dp*x%d3val1_d1val3 + 6.0_dp*q58*x%d2val1 + 9.0_dp*q59 - q111*q118 - q119*q61) - q61*q91 - q80*(-12.0_dp*q18*x%d2val1 + 8.0_dp*x%d3val1 + q117*q20))
      unary%d2val1_d1val2_d2val3 = q39*(-128.0_dp*q94*x%d1val1_d1val2 - 256.0_dp*q125*x%d1val1*x%d1val1_d1val3 - 30.0_dp*q127*q79 - 32.0_dp*q78*q82 - 32.0_dp*q89*(2.0_dp*q64 + 2.0_dp*q67 + 2.0_dp*x%d2val1_d1val2_d1val3 - q124*q46 - q125*q129 + q68 - q77*x%d1val1_d1val2_d1val3) - 64.0_dp*q18*x%d1val1_d1val2_d2val3 + 12.0_dp*q128*q93 + 12.0_dp*q65*q79*x%d2val3 + 128.0_dp*q124*q30 + 192.0_dp*q115*q3*x%d1val2 + 24.0_dp*q126*q79 + 32.0_dp*q27*x%d1val2_d2val3 + 64.0_dp*x%d2val1_d1val2_d2val3 - q107*q34 + q109*q122*x%d1val1_d1val3 + q110*q34 + q113*q123 - q114*q122 - q121*q43 + q121*q58 + q123*q21*q71 - q123*q69*x%d2val3 - q128*q70*x%d2val3 + q78*q92 - q79*q99 + q90*(-2.0_dp*x%d2val1_d2val3 - 4.0_dp*q73 + q129*q75 + q46*q72 - q74 + q77*x%d1val1_d2val3))
      unary%d1val1_d2val2_d2val3 = q120*(-32.0_dp*q130*q75 - 4.0_dp*q131*q94 + 16.0_dp*q65*q87 + 16.0_dp*x%d1val1_d2val2_d2val3 + 32.0_dp*q126*x%d1val1_d1val2 + 8.0_dp*q133*q58 + q109*q131*q15 + q131*q132*x%d1val1_d1val3 - q131*q46*x%d1val1_d2val3 - q133*q70*x%d1val1_d1val3 + q37*(-12.0_dp*q102 - 2.0_dp*x%d2val2_d2val3 - 3.0_dp*q103 + 6.0_dp*q12*q75 + q100*q118 + q134*x%d1val2_d2val3) - q55*(-4.0_dp*q82*x%d1val1_d1val2 - 4.0_dp*q83 + 4.0_dp*x%d1val1_d2val2_d1val3 + q130*q132 - q133*q18 - q135*q43 + q135*q58) - q80*(8.0_dp*x%d1val1_d2val2 - q130*q70 - q131*q18) - q82*q88 + q84*q88 - q86*x%d1val1_d1val2 - q90*x%d1val1_d1val2_d2val3)
      unary%d3val2_d2val3 = q120*(-144.0_dp*q12*q60*x%d1val2_d1val3*x%d1val3 - 24.0_dp*q33*x%d2val2_d2val3 - 24.0_dp*q4*x%d1val2_d2val3*x%d2val2 - 48.0_dp*q127*x%d2val2 + 16.0_dp*x%d3val2_d2val3 + 36.0_dp*q96*x%d1val2_d2val3 + 72.0_dp*q100*q65 + q112*q137 + q116*q50 - q136*q82 + q136*q84 - q55*(4.0_dp*x%d3val2_d1val3 + 6.0_dp*q84*x%d2val2 + 9.0_dp*q97 - q118*q137 - q119*q98 - q134*x%d2val2_d1val3) + q65*q91*x%d2val2 - q80*(-12.0_dp*q33*x%d2val2 + 8.0_dp*x%d3val2 + q117*q50) - q91*q98)
   end function sqrt_self
   
   function pow2_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
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
      q0 = 2.0_dp*x%val
      q1 = 2.0_dp*x%d1val1
      q2 = 2.0_dp*x%d1val2
      q3 = 6.0_dp*x%d1val1
      q4 = 4.0_dp*x%d1val1
      q5 = 2.0_dp*x%d1val3
      q6 = 4.0_dp*x%d1val2
      q7 = 4.0_dp*x%d1val3
      q8 = 6.0_dp*x%d1val2
      q9 = 6.0_dp*x%d2val1
      q10 = 4.0_dp*x%d1val1_d1val2
      q11 = 2.0_dp*x%d2val1
      q12 = 2.0_dp*x%d2val2
      q13 = 2.0_dp*x%d2val3
      q14 = 4.0_dp*x%d1val2_d1val3
      q15 = 6.0_dp*x%d2val2
      q16 = 8.0_dp*x%d1val1_d1val2_d1val3
      unary%val = pow2(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = 2.0_dp*pow2(x%d1val1) + q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q1*x%d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q1*x%d1val3
      unary%d2val2 = 2.0_dp*pow2(x%d1val2) + q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q2*x%d1val3
      unary%d2val3 = 2.0_dp*pow2(x%d1val3) + q0*x%d2val3
      unary%d3val1 = q0*x%d3val1 + q3*x%d2val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2 + q2*x%d2val1 + q4*x%d1val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3 + q4*x%d1val1_d1val3 + q5*x%d2val1
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2 + q1*x%d2val2 + q6*x%d1val1_d1val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3 + q1*x%d1val2_d1val3 + q2*x%d1val1_d1val3 + q5*x%d1val1_d1val2
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3 + q1*x%d2val3 + q7*x%d1val1_d1val3
      unary%d3val2 = q0*x%d3val2 + q8*x%d2val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3 + q5*x%d2val2 + q6*x%d1val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3 + q2*x%d2val3 + q7*x%d1val2_d1val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3 + q3*x%d2val1_d1val3 + q5*x%d3val1 + q9*x%d1val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3 + q10*x%d1val1_d1val3 + q11*x%d1val2_d1val3 + q2*x%d2val1_d1val3 + q4*x%d1val1_d1val2_d1val3 + q5*x%d2val1_d1val2
      unary%d2val1_d2val3 = 4.0_dp*pow2(x%d1val1_d1val3) + q0*x%d2val1_d2val3 + q11*x%d2val3 + q4*x%d1val1_d2val3 + q7*x%d2val1_d1val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3 + q1*x%d2val2_d1val3 + q10*x%d1val2_d1val3 + q12*x%d1val1_d1val3 + q5*x%d1val1_d2val2 + q6*x%d1val1_d1val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3 + q1*x%d1val2_d2val3 + q13*x%d1val1_d1val2 + q14*x%d1val1_d1val3 + q2*x%d1val1_d2val3 + q7*x%d1val1_d1val2_d1val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3 + q15*x%d1val2_d1val3 + q5*x%d3val2 + q8*x%d2val2_d1val3
      unary%d2val2_d2val3 = 4.0_dp*pow2(x%d1val2_d1val3) + q0*x%d2val2_d2val3 + q12*x%d2val3 + q6*x%d1val2_d2val3 + q7*x%d2val2_d1val3
      unary%d3val1_d2val3 = 12.0_dp*x%d1val1_d1val3*x%d2val1_d1val3 + q0*x%d3val1_d2val3 + q13*x%d3val1 + q3*x%d2val1_d2val3 + q7*x%d3val1_d1val3 + q9*x%d1val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3 + q10*x%d1val1_d2val3 + q11*x%d1val2_d2val3 + q13*x%d2val1_d1val2 + q14*x%d2val1_d1val3 + q16*x%d1val1_d1val3 + q2*x%d2val1_d2val3 + q4*x%d1val1_d1val2_d2val3 + q7*x%d2val1_d1val2_d1val3
      unary%d1val1_d2val2_d2val3 = 4.0_dp*x%d1val1_d1val3*x%d2val2_d1val3 + q0*x%d1val1_d2val2_d2val3 + q1*x%d2val2_d2val3 + q10*x%d1val2_d2val3 + q12*x%d1val1_d2val3 + q13*x%d1val1_d2val2 + q16*x%d1val2_d1val3 + q6*x%d1val1_d1val2_d2val3 + q7*x%d1val1_d2val2_d1val3
      unary%d3val2_d2val3 = 12.0_dp*x%d1val2_d1val3*x%d2val2_d1val3 + q0*x%d3val2_d2val3 + q13*x%d3val2 + q15*x%d1val2_d2val3 + q7*x%d3val2_d1val3 + q8*x%d2val2_d2val3
   end function pow2_self
   
   function pow3_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = 3.0_dp*pow2(x%val)
      q1 = 24.0_dp*x%val
      q2 = pow2(x%d1val1)
      q3 = 0.125_dp*x%val
      q4 = x%d1val1*x%val
      q5 = 6.0_dp*q4
      q6 = pow2(x%d1val2)
      q7 = 48.0_dp*q6
      q8 = x%d1val2*x%val
      q9 = 6.0_dp*x%d1val3
      q10 = 2.0_dp*x%val
      q11 = q10*x%d2val3
      q12 = pow2(x%d1val3)
      q13 = 1.5_dp*x%val
      q14 = 18.0_dp*x%d2val1
      q15 = 72.0_dp*x%d2val1
      q16 = 0.041666666666666664_dp*q15*x%val + 6.0_dp*q2
      q17 = 288.0_dp*x%d1val1
      q18 = 72.0_dp*x%val
      q19 = q15*x%d1val2 + q17*x%d1val1_d1val2 + q18*x%d2val1_d1val2
      q20 = 0.041666666666666664_dp*x%val
      q21 = q15*x%d1val3 + q17*x%d1val1_d1val3 + q18*x%d2val1_d1val3
      q22 = 12.0_dp*x%d1val1_d1val2
      q23 = x%d2val2*x%val
      q24 = x%d1val1*x%d1val2
      q25 = q9*x%val
      q26 = 6.0_dp*x%d1val1_d1val3
      q27 = x%d1val1_d1val3*x%val
      q28 = 12.0_dp*x%d1val3
      q29 = 3.0_dp*q11 + 6.0_dp*q12
      q30 = 18.0_dp*x%d1val2
      q31 = 0.041666666666666664_dp*x%d1val3
      q32 = x%d1val2*x%d1val2_d1val3
      q33 = 72.0_dp*x%d1val3
      q34 = x%d1val2_d1val3*x%val
      q35 = q14*x%d1val1
      q36 = 18.0_dp*q4
      q37 = 18.0_dp*q2
      q38 = x%d1val1_d1val2*x%d1val1_d1val3
      q39 = x%d1val2*x%d2val1_d1val3
      q40 = 4.0_dp*x%d2val1
      q41 = 8.0_dp*q2 + q40*x%val
      q42 = 0.75_dp*x%d2val3
      q43 = 8.0_dp*x%d1val1
      q44 = 2.0_dp*x%d2val1
      q45 = q10*x%d2val1_d1val3 + q43*x%d1val1_d1val3 + q44*x%d1val3
      q46 = 3.0_dp*x%d1val3
      q47 = 4.0_dp*x%d1val3
      q48 = pow2(x%d1val1_d1val3)
      q49 = 8.0_dp*q48 + q10*x%d2val1_d2val3 + q43*x%d1val1_d2val3 + q44*x%d2val3 + q47*x%d2val1_d1val3
      q50 = q22*x%d1val2
      q51 = 12.0_dp*q8
      q52 = x%d1val3*x%d2val2
      q53 = 6.0_dp*x%d2val3
      q54 = q53*x%val
      q55 = q28*x%val
      q56 = 6.0_dp*q12
      q57 = 18.0_dp*q8
      q58 = 18.0_dp*q23
      q59 = 18.0_dp*q6
      q60 = 4.0_dp*q23
      q61 = 2.0_dp*q52 + q10*x%d2val2_d1val3
      q62 = pow2(x%d1val2_d1val3)
      q63 = x%d1val2*x%d1val2_d2val3
      q64 = 2.0_dp*x%d2val3
      q65 = q10*x%d2val2_d2val3 + q47*x%d2val2_d1val3 + q64*x%d2val2
      q66 = 36.0_dp*x%d1val1
      q67 = 24.0_dp*x%d1val3
      q68 = 36.0_dp*x%d1val2
      unary%val = pow3(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q3*(48.0_dp*q2 + q1*x%d2val1)
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q5*x%d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q5*x%d1val3
      unary%d2val2 = q3*(q1*x%d2val2 + q7)
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q8*q9
      unary%d2val3 = q13*(4.0_dp*q12 + q11)
      unary%d3val1 = 6.0_dp*pow3(x%d1val1) + q0*x%d3val1 + q14*q4
      unary%d2val1_d1val2 = q16*x%d1val2 + q19*q20
      unary%d2val1_d1val3 = q16*x%d1val3 + q20*q21
      unary%d1val1_d2val2 = 0.125_dp*x%d1val1*(48.0_dp*q23 + q7) + q0*x%d1val1_d2val2 + q22*q8
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3 + q24*q9 + q25*x%d1val1_d1val2 + q26*q8 + q5*x%d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3 + q27*q28 + q29*x%d1val1
      unary%d3val2 = 6.0_dp*pow3(x%d1val2) + q0*x%d3val2 + q23*q30
      unary%d2val2_d1val3 = q20*(288.0_dp*q32 + q18*x%d2val2_d1val3 + q33*x%d2val2) + q31*(144.0_dp*q6 + q18*x%d2val2)
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3 + q28*q34 + q29*x%d1val2
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3 + q14*q27 + q25*x%d3val1 + q35*x%d1val3 + q36*x%d2val1_d1val3 + q37*x%d1val1_d1val3
      unary%d2val1_d1val2_d1val3 = 0.041666666666666664_dp*q21*x%d1val2 + q16*x%d1val2_d1val3 + q19*q31 + q20*(288.0_dp*q38 + 72.0_dp*q39 + q15*x%d1val2_d1val3 + q17*x%d1val1_d1val2_d1val3 + q18*x%d2val1_d1val2_d1val3 + q33*x%d2val1_d1val2)
      unary%d2val1_d2val3 = q13*q49 + q41*q42 + q45*q46
      unary%d1val1_d2val2_d1val3 = 4.0_dp*x%d1val1*(1.5_dp*q52 + 3.0_dp*q32 + q13*x%d2val2_d1val3) + 4.0_dp*x%d1val1_d1val3*(1.5_dp*q6 + q13*x%d2val2) + q0*x%d1val1_d2val2_d1val3 + q22*q34 + q25*x%d1val1_d2val2 + q50*x%d1val3 + q51*x%d1val1_d1val2_d1val3
      unary%d1val1_d1val2_d2val3 = 12.0_dp*q27*x%d1val2_d1val3 + 6.0_dp*q8*x%d1val1_d2val3 + q0*x%d1val1_d1val2_d2val3 + q24*q53 + q28*x%d1val1*x%d1val2_d1val3 + q28*x%d1val1_d1val3*x%d1val2 + q5*x%d1val2_d2val3 + q54*x%d1val1_d1val2 + q55*x%d1val1_d1val2_d1val3 + q56*x%d1val1_d1val2
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3 + q25*x%d3val2 + q30*q52 + q57*x%d2val2_d1val3 + q58*x%d1val2_d1val3 + q59*x%d1val2_d1val3
      unary%d2val2_d2val3 = q13*(8.0_dp*q62 + 8.0_dp*q63 + q65) + q42*(8.0_dp*q6 + q60) + q46*(8.0_dp*q32 + q61)
      unary%d3val1_d2val3 = 36.0_dp*q27*x%d2val1_d1val3 + 36.0_dp*x%d1val1_d1val3*x%d1val3*x%d2val1 + q0*x%d3val1_d2val3 + q14*x%d1val1_d2val3*x%val + q35*x%d2val3 + q36*x%d2val1_d2val3 + q37*x%d1val1_d2val3 + q48*q66 + q54*x%d3val1 + q55*x%d3val1_d1val3 + q56*x%d3val1 + q66*x%d1val3*x%d2val1_d1val3
      unary%d2val1_d1val2_d2val3 = 0.75_dp*q41*x%d1val2_d2val3 + 1.5_dp*q49*x%d1val2 + 3.0_dp*q45*x%d1val2_d1val3 + q13*(16.0_dp*x%d1val1_d1val2_d1val3*x%d1val1_d1val3 + 2.0_dp*x%d1val2*x%d2val1_d2val3 + 4.0_dp*x%d1val2_d1val3*x%d2val1_d1val3 + 8.0_dp*x%d1val1_d1val2*x%d1val1_d2val3 + q10*x%d2val1_d1val2_d2val3 + q43*x%d1val1_d1val2_d2val3 + q44*x%d1val2_d2val3 + q47*x%d2val1_d1val2_d1val3 + q64*x%d2val1_d1val2) + q42*(16.0_dp*x%d1val1*x%d1val1_d1val2 + 4.0_dp*x%d2val1_d1val2*x%val + q40*x%d1val2) + q46*(2.0_dp*q39 + 2.0_dp*x%d1val3*x%d2val1_d1val2 + 8.0_dp*q38 + q10*x%d2val1_d1val2_d1val3 + q43*x%d1val1_d1val2_d1val3 + q44*x%d1val2_d1val3)
      unary%d1val1_d2val2_d2val3 = 1.5_dp*x%d1val1_d2val3*(4.0_dp*q6 + q60) + 3.0_dp*x%d1val1*(4.0_dp*q62 + 4.0_dp*q63 + q65) + q0*x%d1val1_d2val2_d2val3 + q1*x%d1val1_d1val2_d1val3*x%d1val2_d1val3 + q22*x%d1val2_d2val3*x%val + q26*(4.0_dp*q32 + q61) + q50*x%d2val3 + q51*x%d1val1_d1val2_d2val3 + q54*x%d1val1_d2val2 + q55*x%d1val1_d2val2_d1val3 + q56*x%d1val1_d2val2 + q67*x%d1val1_d1val2*x%d1val2_d1val3 + q67*x%d1val1_d1val2_d1val3*x%d1val2
      unary%d3val2_d2val3 = 36.0_dp*q34*x%d2val2_d1val3 + 36.0_dp*q52*x%d1val2_d1val3 + q0*x%d3val2_d2val3 + q30*x%d2val2*x%d2val3 + q54*x%d3val2 + q55*x%d3val2_d1val3 + q56*x%d3val2 + q57*x%d2val2_d2val3 + q58*x%d1val2_d2val3 + q59*x%d1val2_d2val3 + q62*q68 + q68*x%d1val3*x%d2val2_d1val3
   end function pow3_self
   
   function pow4_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = 4.0_dp*pow3(x%val)
      q1 = 48.0_dp*x%val
      q2 = pow2(x%d1val1)
      q3 = 144.0_dp*q2 + q1*x%d2val1
      q4 = pow2(x%val)
      q5 = 0.08333333333333333_dp*q4
      q6 = 12.0_dp*x%d1val1
      q7 = q4*q6
      q8 = pow2(x%d1val2)
      q9 = 144.0_dp*q8 + q1*x%d2val2
      q10 = q4*x%d1val3
      q11 = 12.0_dp*x%d1val2
      q12 = 2.0_dp*x%d2val3
      q13 = q12*x%val
      q14 = pow2(x%d1val3)
      q15 = 2.0_dp*q4
      q16 = pow3(x%d1val1)
      q17 = 2592.0_dp*x%d1val1
      q18 = x%d2val1*x%val
      q19 = 288.0_dp*q4
      q20 = 1728.0_dp*q16 + q17*q18 + q19*x%d3val1
      q21 = 0.013888888888888888_dp*x%val
      q22 = x%d1val2*x%val
      q23 = 0.16666666666666666_dp*q3
      q24 = 288.0_dp*x%d1val1
      q25 = 48.0_dp*x%d2val1
      q26 = q1*x%d2val1_d1val2 + q24*x%d1val1_d1val2 + q25*x%d1val2
      q27 = x%d1val3*x%val
      q28 = 1728.0_dp*x%d1val1_d1val2
      q29 = x%d2val2*x%val
      q30 = 36.0_dp*x%d1val1
      q31 = q19*x%d1val1_d2val2 + q22*q28 + q30*(24.0_dp*q29 + 48.0_dp*q8)
      q32 = 24.0_dp*x%d1val1
      q33 = q22*q32
      q34 = 12.0_dp*x%d1val1_d1val2
      q35 = 12.0_dp*x%d1val1_d1val3
      q36 = 24.0_dp*x%d1val1_d1val3
      q37 = 4.0_dp*q4
      q38 = 24.0_dp*q14 + 6.0_dp*q13
      q39 = pow3(x%d1val2)
      q40 = 2592.0_dp*x%d1val2
      q41 = 1728.0_dp*q39 + q19*x%d3val2 + q29*q40
      q42 = 0.16666666666666666_dp*q27
      q43 = x%d1val2*x%d1val2_d1val3
      q44 = 48.0_dp*x%d1val3
      q45 = 24.0_dp*x%d1val2_d1val3
      q46 = 0.013888888888888888_dp*x%d1val3
      q47 = x%d1val3*x%d2val1
      q48 = x%d2val1_d1val3*x%val
      q49 = q18*x%d1val1_d1val3
      q50 = 576.0_dp*q27
      q51 = q2*x%d1val1_d1val3
      q52 = 0.027777777777777776_dp*x%d2val1
      q53 = 24.0_dp*q2 + 288.0_dp*q52*x%val
      q54 = x%d1val2*x%d1val3
      q55 = x%d1val2_d1val3*x%val
      q56 = x%d1val1*x%d1val1_d1val3
      q57 = x%d1val1_d1val2*x%d1val1_d1val3
      q58 = x%d1val2*x%d2val1_d1val3
      q59 = 2.0_dp*x%d1val3
      q60 = 12.0_dp*q56 + 2.0_dp*q48 + q59*x%d2val1
      q61 = 8.0_dp*q27
      q62 = 4.0_dp*x%d2val1
      q63 = 12.0_dp*q2 + q62*x%val
      q64 = 2.0_dp*q14
      q65 = q13 + q64
      q66 = 4.0_dp*x%d1val3
      q67 = 2.0_dp*x%val
      q68 = pow2(x%d1val1_d1val3)
      q69 = 12.0_dp*q68 + q12*x%d2val1 + q6*x%d1val1_d2val3 + q66*x%d2val1_d1val3 + q67*x%d2val1_d2val3
      q70 = q22*x%d1val1_d1val2_d1val3
      q71 = x%d1val3*x%d2val2
      q72 = x%d2val2_d1val3*x%val
      q73 = q1*x%d1val2_d1val3
      q74 = 24.0_dp*x%d1val1_d1val2
      q75 = q74*x%val
      q76 = q29*x%d1val2_d1val3
      q77 = q8*x%d1val2_d1val3
      q78 = q59*x%d2val2 + q67*x%d2val2_d1val3
      q79 = 4.0_dp*q29
      q80 = pow2(x%d1val2_d1val3)
      q81 = x%d1val2*x%d1val2_d2val3
      q82 = q12*x%d2val2 + q66*x%d2val2_d1val3 + q67*x%d2val2_d2val3
      q83 = 72.0_dp*x%d1val1
      q84 = 8.0_dp*x%d3val1
      q85 = 0.5_dp*x%d2val3
      q86 = 72.0_dp*x%d1val1_d1val3
      q87 = 16.0_dp*q27
      q88 = x%d2val3*x%val
      q89 = 8.0_dp*q60
      q90 = 4.0_dp*x%val
      q91 = q32*x%d1val1_d1val2 + q62*x%d1val2 + q90*x%d2val1_d1val2
      q92 = 2.0_dp*x%d2val1
      q93 = 8.0_dp*x%d1val1_d2val2
      q94 = 8.0_dp*q8 + q79
      q95 = 6.0_dp*x%d1val1
      q96 = 3.0_dp*q94
      q97 = 8.0_dp*q43 + q78
      q98 = x%d1val2*x%d2val3
      q99 = 8.0_dp*x%d3val2
      q100 = 36.0_dp*x%d1val2
      q101 = 72.0_dp*x%d1val2_d1val3
      unary%val = pow4(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q3*q5
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q7*x%d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q7*x%d1val3
      unary%d2val2 = q5*q9
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q10*q11
      unary%d2val3 = q15*(6.0_dp*q14 + q13)
      unary%d3val1 = q20*q21
      unary%d2val1_d1val2 = q22*q23 + q26*q5
      unary%d2val1_d1val3 = q23*q27 + q5*(q1*x%d2val1_d1val3 + q24*x%d1val1_d1val3 + q25*x%d1val3)
      unary%d1val1_d2val2 = q21*q31
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3 + q10*q34 + q33*x%d1val3 + q35*q4*x%d1val2 + q7*x%d1val2_d1val3
      unary%d1val1_d2val3 = x%val*(q27*q36 + q37*x%d1val1_d2val3 + q38*x%d1val1)
      unary%d3val2 = q21*q41
      unary%d2val2_d1val3 = q42*q9 + q5*(288.0_dp*q43 + q1*x%d2val2_d1val3 + q44*x%d2val2)
      unary%d1val2_d2val3 = x%val*(q27*q45 + q37*x%d1val2_d2val3 + q38*x%d1val2)
      unary%d3val1_d1val3 = q20*q46 + q21*(2592.0_dp*q49 + 5184.0_dp*q51 + q17*q47 + q17*q48 + q19*x%d3val1_d1val3 + q50*x%d3val1)
      unary%d2val1_d1val2_d1val3 = 288.0_dp*q22*(0.027777777777777776_dp*q48 + 0.16666666666666666_dp*q56 + q52*x%d1val3) + q26*q42 + q5*(288.0_dp*q57 + 48.0_dp*q58 + q1*x%d2val1_d1val2_d1val3 + q24*x%d1val1_d1val2_d1val3 + q25*x%d1val2_d1val3 + q44*x%d2val1_d1val2) + q53*q54 + q53*q55
      unary%d2val1_d2val3 = q15*q69 + q60*q61 + q63*q65
      unary%d1val1_d2val2_d1val3 = q21*(1728.0_dp*q70 + 4.0_dp*x%d1val1*(216.0_dp*q71 + 216.0_dp*q72 + 864.0_dp*q43) + 4.0_dp*x%d1val1_d1val3*(216.0_dp*q29 + 432.0_dp*q8) + q19*x%d1val1_d2val2_d1val3 + q28*q54 + q28*q55 + q50*x%d1val1_d2val2) + q31*q46
      unary%d1val1_d1val2_d2val3 = 24.0_dp*q10*x%d1val1_d1val2_d1val3 + q0*x%d1val1_d1val2_d2val3 + q1*q54*x%d1val1_d1val3 + q11*q4*x%d1val1_d2val3 + q14*q32*x%d1val2 + q14*q75 + q33*x%d2val3 + q34*q4*x%d2val3 + q36*q4*x%d1val2_d1val3 + q7*x%d1val2_d2val3 + q73*x%d1val1*x%d1val3
      unary%d3val2_d1val3 = q21*(2592.0_dp*q76 + 5184.0_dp*q77 + q19*x%d3val2_d1val3 + q40*q71 + q40*q72 + q50*x%d3val2) + q41*q46
      unary%d2val2_d2val3 = q15*(12.0_dp*q80 + 12.0_dp*q81 + q82) + q61*(12.0_dp*q43 + q78) + q65*(12.0_dp*q8 + q79)
      unary%d3val1_d2val3 = q59*(36.0_dp*q49 + 72.0_dp*q51 + q27*q84 + q30*q47 + q30*q48 + q37*x%d3val1_d1val3) + q85*(48.0_dp*q16 + q18*q83 + q4*q84) + x%val*(144.0_dp*q68*x%d1val1 + 36.0_dp*q18*x%d1val1_d2val3 + 72.0_dp*q2*x%d1val1_d2val3 + q14*q84 + q30*x%d2val1*x%d2val3 + q30*x%d2val1_d2val3*x%val + q37*x%d3val1_d2val3 + q47*q86 + q48*q86 + q83*x%d1val3*x%d2val1_d1val3 + q84*q88 + q87*x%d3val1_d1val3)
      unary%d2val1_d1val2_d2val3 = q12*q63*x%d1val2 + q13*q91 + q15*(2.0_dp*x%d1val2*x%d2val1_d2val3 + 4.0_dp*x%d1val2_d1val3*x%d2val1_d1val3 + q12*x%d2val1_d1val2 + q34*x%d1val1_d2val3 + q36*x%d1val1_d1val2_d1val3 + q6*x%d1val1_d1val2_d2val3 + q66*x%d2val1_d1val2_d1val3 + q67*x%d2val1_d1val2_d2val3 + q92*x%d1val2_d2val3) + q54*q89 + q55*q89 + q61*(12.0_dp*q57 + 2.0_dp*q58 + q59*x%d2val1_d1val2 + q6*x%d1val1_d1val2_d1val3 + q67*x%d2val1_d1val2_d1val3 + q92*x%d1val2_d1val3) + q63*q66*x%d1val2_d1val3 + q63*q67*x%d1val2_d2val3 + q64*q91 + q69*q90*x%d1val2
      unary%d1val1_d2val2_d2val3 = q59*(24.0_dp*q70 + q27*q93 + q37*x%d1val1_d2val2_d1val3 + q45*x%d1val1_d1val2*x%val + q54*q74 + q95*q97 + q96*x%d1val1_d1val3) + q85*(q1*x%d1val1_d1val2*x%d1val2 + q4*q93 + q94*q95) + x%val*(24.0_dp*q22*x%d1val1_d1val2_d2val3 + q14*q93 + q35*q97 + q37*x%d1val1_d2val2_d2val3 + q44*x%d1val1_d1val2*x%d1val2_d1val3 + q44*x%d1val1_d1val2_d1val3*x%d1val2 + q73*x%d1val1_d1val2_d1val3 + q74*q98 + q75*x%d1val2_d2val3 + q87*x%d1val1_d2val2_d1val3 + q88*q93 + q95*(8.0_dp*q80 + 8.0_dp*q81 + q82) + q96*x%d1val1_d2val3)
      unary%d3val2_d2val3 = q59*(36.0_dp*q76 + 72.0_dp*q77 + q100*q71 + q100*q72 + q27*q99 + q37*x%d3val2_d1val3) + q85*(48.0_dp*q39 + 72.0_dp*q29*x%d1val2 + q4*q99) + x%val*(144.0_dp*q80*x%d1val2 + 36.0_dp*q22*x%d2val2_d2val3 + 36.0_dp*q29*x%d1val2_d2val3 + 36.0_dp*q98*x%d2val2 + 72.0_dp*q54*x%d2val2_d1val3 + 72.0_dp*q8*x%d1val2_d2val3 + q101*q71 + q101*q72 + q14*q99 + q37*x%d3val2_d2val3 + q87*x%d3val2_d1val3 + q88*q99)
   end function pow4_self
   
   function pow5_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = 5.0_dp*pow4(x%val)
      q1 = 24.0_dp*x%val
      q2 = pow2(x%d1val1)
      q3 = 96.0_dp*q2 + q1*x%d2val1
      q4 = pow3(x%val)
      q5 = 0.20833333333333334_dp*q4
      q6 = 20.0_dp*q4
      q7 = x%d1val1*x%d1val2
      q8 = q6*x%d1val3
      q9 = pow2(x%d1val2)
      q10 = 96.0_dp*q9 + q1*x%d2val2
      q11 = 2.0_dp*x%d2val3
      q12 = q11*x%val
      q13 = pow2(x%d1val3)
      q14 = 8.0_dp*q13
      q15 = 2.5_dp*q4
      q16 = pow3(x%d1val1)
      q17 = x%d1val1*x%d2val1
      q18 = 3456.0_dp*q17
      q19 = pow2(x%val)
      q20 = 288.0_dp*q19
      q21 = 3456.0_dp*q16 + q18*x%val + q20*x%d3val1
      q22 = 0.017361111111111112_dp*q19
      q23 = 0.625_dp*q19
      q24 = q3*x%d1val2
      q25 = 192.0_dp*x%d1val1
      q26 = 24.0_dp*x%d2val1
      q27 = q1*x%d2val1_d1val2 + q25*x%d1val1_d1val2 + q26*x%d1val2
      q28 = q23*x%d1val3
      q29 = q1*x%d2val1_d1val3 + q25*x%d1val1_d1val3 + q26*x%d1val3
      q30 = x%d1val1_d1val2*x%val
      q31 = 2304.0_dp*x%d1val2
      q32 = x%d2val2*x%val
      q33 = 48.0_dp*q32
      q34 = 144.0_dp*q9
      q35 = 24.0_dp*x%d1val1*(q33 + q34) + q20*x%d1val1_d2val2 + q30*q31
      q36 = x%d1val1_d1val3*x%d1val2
      q37 = 2.0_dp*x%val
      q38 = q37*x%d1val3
      q39 = 0.25_dp*q19
      q40 = 0.5_dp*q12 + 3.0_dp*q13
      q41 = 20.0_dp*q19
      q42 = pow3(x%d1val2)
      q43 = 3456.0_dp*x%d1val2
      q44 = 3456.0_dp*q42 + q20*x%d3val2 + q32*q43
      q45 = x%d1val2*x%d1val2_d1val3
      q46 = 24.0_dp*x%d1val3
      q47 = x%d1val3*x%val
      q48 = 0.034722222222222224_dp*q47
      q49 = 3456.0_dp*x%val
      q50 = x%d1val1*x%d2val1_d1val3
      q51 = x%d1val1_d1val3*x%d2val1
      q52 = 576.0_dp*q47
      q53 = q2*x%d1val1_d1val3
      q54 = x%d1val1_d1val2*x%d1val1_d1val3
      q55 = x%d1val2*x%d2val1_d1val3
      q56 = 16.0_dp*x%d1val1
      q57 = 2.0_dp*x%d1val3
      q58 = q37*x%d2val1_d1val3 + q56*x%d1val1_d1val3 + q57*x%d2val1
      q59 = q1*x%d1val3
      q60 = 4.0_dp*x%d2val1
      q61 = 16.0_dp*q2 + q60*x%val
      q62 = 12.0_dp*q13 + 3.0_dp*q12
      q63 = 4.0_dp*x%d1val3
      q64 = pow2(x%d1val1_d1val3)
      q65 = 16.0_dp*q64 + q11*x%d2val1 + q37*x%d2val1_d2val3 + q56*x%d1val1_d2val3 + q63*x%d2val1_d1val3
      q66 = 4.0_dp*q19
      q67 = 0.625_dp*x%val
      q68 = x%d1val1_d1val2*x%d1val3
      q69 = q30*x%d1val2_d1val3
      q70 = x%d1val1_d1val2_d1val3*x%val
      q71 = 4.0_dp*x%d1val1_d1val3
      q72 = x%d1val3*x%d2val2
      q73 = x%d2val2_d1val3*x%val
      q74 = x%d2val3*x%val
      q75 = q47*x%d1val2_d1val3
      q76 = q13*x%d1val2
      q77 = 32.0_dp*x%d1val1
      q78 = q19*x%d1val2_d2val3
      q79 = 32.0_dp*q19
      q80 = x%d1val1_d1val2*x%d2val3
      q81 = 64.0_dp*q19
      q82 = x%d1val1_d1val2_d1val3*x%d1val3
      q83 = q37*x%d2val2_d1val3 + q57*x%d2val2
      q84 = 4.0_dp*q32
      q85 = pow2(x%d1val2_d1val3)
      q86 = x%d1val2*x%d1val2_d2val3
      q87 = q11*x%d2val2 + q37*x%d2val2_d2val3 + q63*x%d2val2_d1val3
      q88 = 96.0_dp*x%val
      q89 = 8.0_dp*q19
      q90 = 0.625_dp*q12 + 1.25_dp*q13
      q91 = 48.0_dp*q17
      q92 = 48.0_dp*x%val
      q93 = 8.0_dp*x%d3val1
      q94 = 5.0_dp*q47
      q95 = 96.0_dp*x%d1val3
      q96 = 16.0_dp*q47
      q97 = 1.25_dp*q19
      q98 = 7.5_dp*q61
      q99 = 15.0_dp*q19
      q100 = 4.0_dp*x%d2val1_d1val2*x%val + q60*x%d1val2 + q77*x%d1val1_d1val2
      q101 = 2.0_dp*x%d2val1
      q102 = 32.0_dp*x%d1val1_d1val2_d1val3
      q103 = 64.0_dp*x%d1val2
      q104 = 12.0_dp*q9 + q84
      q105 = 8.0_dp*x%d1val1
      q106 = 32.0_dp*x%d1val2
      q107 = 8.0_dp*x%d1val1_d2val2
      q108 = 12.0_dp*q45 + q83
      q109 = 64.0_dp*x%d1val2_d1val3
      q110 = 48.0_dp*x%d1val2
      q111 = 8.0_dp*x%d3val2
      q112 = 96.0_dp*x%d1val2_d1val3
      unary%val = pow5(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q3*q5
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q6*q7
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q8*x%d1val1
      unary%d2val2 = q10*q5
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q8*x%d1val2
      unary%d2val3 = q15*(q12 + q14)
      unary%d3val1 = q21*q22
      unary%d2val1_d1val2 = q23*q24 + q27*q5
      unary%d2val1_d1val3 = q28*q3 + q29*q5
      unary%d1val1_d2val2 = q22*q35
      unary%d1val1_d1val2_d1val3 = 60.0_dp*q19*q7*x%d1val3 + q0*x%d1val1_d1val2_d1val3 + q36*q6 + q6*x%d1val1*x%d1val2_d1val3 + q8*x%d1val1_d1val2
      unary%d1val1_d2val3 = q41*(q38*x%d1val1_d1val3 + q39*x%d1val1_d2val3 + q40*x%d1val1)
      unary%d3val2 = q22*q44
      unary%d2val2_d1val3 = q10*q28 + q5*(192.0_dp*q45 + q1*x%d2val2_d1val3 + q46*x%d2val2)
      unary%d1val2_d2val3 = q41*(q38*x%d1val2_d1val3 + q39*x%d1val2_d2val3 + q40*x%d1val2)
      unary%d3val1_d1val3 = q21*q48 + q22*(10368.0_dp*q53 + q18*x%d1val3 + q20*x%d3val1_d1val3 + q49*q50 + q49*q51 + q52*x%d3val1)
      unary%d2val1_d1val2_d1val3 = 1.25_dp*q24*q47 + q23*q29*x%d1val2 + q23*q3*x%d1val2_d1val3 + q27*q28 + q5*(192.0_dp*q54 + 24.0_dp*q55 + q1*x%d2val1_d1val2_d1val3 + q25*x%d1val1_d1val2_d1val3 + q26*x%d1val2_d1val3 + q46*x%d2val1_d1val2)
      unary%d2val1_d2val3 = q67*(q58*q59 + q61*q62 + q65*q66)
      unary%d1val1_d2val2_d1val3 = q22*(2304.0_dp*q69 + 4.0_dp*x%d1val1*(1728.0_dp*q45 + 288.0_dp*q72 + 288.0_dp*q73) + q20*x%d1val1_d2val2_d1val3 + q31*q68 + q31*q70 + q52*x%d1val1_d2val2 + q71*(288.0_dp*q32 + 864.0_dp*q9)) + q35*q48
      unary%d1val1_d1val2_d2val3 = q67*(192.0_dp*q36*q47 + 8.0_dp*q4*x%d1val1_d1val2_d2val3 + 96.0_dp*q13*q30 + 96.0_dp*q7*q74 + q25*q75 + q25*q76 + q77*q78 + q79*q80 + q79*x%d1val1_d2val3*x%d1val2 + q81*q82 + q81*x%d1val1_d1val3*x%d1val2_d1val3)
      unary%d3val2_d1val3 = q22*(10368.0_dp*q9*x%d1val2_d1val3 + 3456.0_dp*q32*x%d1val2_d1val3 + q20*x%d3val2_d1val3 + q43*q72 + q43*q73 + q52*x%d3val2) + q44*q48
      unary%d2val2_d2val3 = q67*(q59*(16.0_dp*q45 + q83) + q62*(16.0_dp*q9 + q84) + q66*(16.0_dp*q85 + 16.0_dp*q86 + q87))
      unary%d3val1_d2val3 = q90*(96.0_dp*q16 + q17*q88 + q89*x%d3val1) + q94*(144.0_dp*q53 + q47*q93 + q50*q92 + q51*q92 + q66*x%d3val1_d1val3 + q91*x%d1val3) + q97*(144.0_dp*q2*x%d1val1_d2val3 + 288.0_dp*q64*x%d1val1 + q14*x%d3val1 + q50*q95 + q51*q95 + q66*x%d3val1_d2val3 + q74*q93 + q88*x%d1val1_d1val3*x%d2val1_d1val3 + q91*x%d2val3 + q92*x%d1val1*x%d2val1_d2val3 + q92*x%d1val1_d2val3*x%d2val1 + q96*x%d3val1_d1val3)
      unary%d2val1_d1val2_d2val3 = 15.0_dp*q61*q75 + 3.75_dp*q100*q19*x%d2val3 + 3.75_dp*q61*q78 + 30.0_dp*q47*q58*x%d1val2 + 7.5_dp*q100*q13*x%val + 7.5_dp*q19*q65*x%d1val2 + q15*(16.0_dp*x%d1val1_d1val2*x%d1val1_d2val3 + 2.0_dp*x%d1val2*x%d2val1_d2val3 + 4.0_dp*x%d1val2_d1val3*x%d2val1_d1val3 + q101*x%d1val2_d2val3 + q102*x%d1val1_d1val3 + q11*x%d2val1_d1val2 + q37*x%d2val1_d1val2_d2val3 + q56*x%d1val1_d1val2_d2val3 + q63*x%d2val1_d1val2_d1val3) + q58*q99*x%d1val2_d1val3 + q74*q98*x%d1val2 + q76*q98 + q99*x%d1val3*(16.0_dp*q54 + 2.0_dp*q55 + q101*x%d1val2_d1val3 + q37*x%d2val1_d1val2_d1val3 + q56*x%d1val1_d1val2_d1val3 + q57*x%d2val1_d1val2)
      unary%d1val1_d2val2_d2val3 = q90*(q103*q30 + q104*q105 + q89*x%d1val1_d2val2) + q94*(32.0_dp*q69 + q102*x%d1val2*x%val + q104*q71 + q105*q108 + q106*q68 + q107*q47 + q66*x%d1val1_d2val2_d1val3) + q97*(16.0_dp*q108*x%d1val1_d1val3 + 32.0_dp*q30*x%d1val2_d2val3 + 4.0_dp*q104*x%d1val1_d2val3 + q103*q82 + q105*(12.0_dp*q85 + 12.0_dp*q86 + q87) + q106*q80 + q106*x%d1val1_d1val2_d2val3*x%val + q107*q74 + q109*q68 + q109*q70 + q14*x%d1val1_d2val2 + q66*x%d1val1_d2val2_d2val3 + q96*x%d1val1_d2val2_d1val3)
      unary%d3val2_d2val3 = q90*(96.0_dp*q32*x%d1val2 + 96.0_dp*q42 + q89*x%d3val2) + q94*(q110*q72 + q110*q73 + q111*q47 + q33*x%d1val2_d1val3 + q34*x%d1val2_d1val3 + q66*x%d3val2_d1val3) + q97*(288.0_dp*q85*x%d1val2 + q110*x%d2val2*x%d2val3 + q111*q74 + q112*q72 + q112*q73 + q14*x%d3val2 + q33*x%d1val2_d2val3 + q34*x%d1val2_d2val3 + q66*x%d3val2_d2val3 + q92*x%d1val2*x%d2val2_d2val3 + q95*x%d1val2*x%d2val2_d1val3 + q96*x%d3val2_d1val3)
   end function pow5_self
   
   function pow6_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = 6.0_dp*pow5(x%val)
      q1 = 48.0_dp*x%val
      q2 = pow2(x%d1val1)
      q3 = 240.0_dp*q2
      q4 = q1*x%d2val1 + q3
      q5 = pow4(x%val)
      q6 = 0.125_dp*q5
      q7 = 30.0_dp*q5
      q8 = x%d1val1*x%d1val2
      q9 = q7*x%d1val3
      q10 = pow2(x%d1val2)
      q11 = 240.0_dp*q10
      q12 = q1*x%d2val2 + q11
      q13 = 2.0_dp*x%d2val3
      q14 = q13*x%val
      q15 = pow2(x%d1val3)
      q16 = 10.0_dp*q15
      q17 = pow3(x%d1val1)
      q18 = x%d2val1*x%val
      q19 = 4320.0_dp*x%d1val1
      q20 = pow2(x%val)
      q21 = 288.0_dp*q20
      q22 = 5760.0_dp*q17 + q18*q19 + q21*x%d3val1
      q23 = pow3(x%val)
      q24 = 0.020833333333333332_dp*q23
      q25 = 0.5_dp*q23
      q26 = q25*q4
      q27 = 480.0_dp*x%d1val1
      q28 = 48.0_dp*x%d2val1
      q29 = q1*x%d2val1_d1val2 + q27*x%d1val1_d1val2 + q28*x%d1val2
      q30 = q1*x%d2val1_d1val3 + q27*x%d1val1_d1val3 + q28*x%d1val3
      q31 = x%d1val1_d1val2*x%val
      q32 = 2880.0_dp*x%d1val2
      q33 = x%d2val2*x%val
      q34 = 60.0_dp*x%d1val1
      q35 = q21*x%d1val1_d2val2 + q31*q32 + q34*(24.0_dp*q33 + 96.0_dp*q10)
      q36 = 40.0_dp*x%d1val1_d1val3
      q37 = x%d1val3*x%val
      q38 = 4.0_dp*q20
      q39 = 8.0_dp*q15
      q40 = 10.0_dp*q14 + 10.0_dp*q39
      q41 = 1.5_dp*q23
      q42 = pow3(x%d1val2)
      q43 = 4320.0_dp*x%d1val2
      q44 = 5760.0_dp*q42 + q21*x%d3val2 + q33*q43
      q45 = q25*x%d1val3
      q46 = x%d1val2*x%d1val2_d1val3
      q47 = 48.0_dp*x%d1val3
      q48 = q37*x%d1val2_d1val3
      q49 = q20*x%d1val3
      q50 = 0.0625_dp*q49
      q51 = x%d1val3*x%d2val1
      q52 = x%d2val1_d1val3*x%val
      q53 = q18*x%d1val1_d1val3
      q54 = 576.0_dp*q37
      q55 = x%d1val1_d1val2*x%d1val1_d1val3
      q56 = x%d1val2*x%d2val1_d1val3
      q57 = 20.0_dp*x%d1val1
      q58 = 2.0_dp*x%d1val3
      q59 = 2.0_dp*x%val
      q60 = q57*x%d1val1_d1val3 + q58*x%d2val1 + q59*x%d2val1_d1val3
      q61 = q59*x%d1val3
      q62 = 4.0_dp*x%d2val1
      q63 = 20.0_dp*q2 + q62*x%val
      q64 = 0.25_dp*q14 + 1.5_dp*q15
      q65 = 4.0_dp*x%d1val3
      q66 = pow2(x%d1val1_d1val3)
      q67 = 20.0_dp*q66 + q13*x%d2val1 + q57*x%d1val1_d2val3 + q59*x%d2val1_d2val3 + q65*x%d2val1_d1val3
      q68 = 0.25_dp*q20
      q69 = 12.0_dp*q20
      q70 = x%d1val1_d1val2*x%d1val3
      q71 = x%d1val1_d1val2_d1val3*x%val
      q72 = x%d1val3*x%d2val2
      q73 = x%d2val2_d1val3*x%val
      q74 = 10.0_dp*x%d1val1
      q75 = x%d2val3*x%val
      q76 = 20.0_dp*x%d1val1_d1val3
      q77 = q37*x%d1val2
      q78 = 2.5_dp*q20
      q79 = q20*x%d1val2_d1val3
      q80 = q33*x%d1val2_d1val3
      q81 = q58*x%d2val2 + q59*x%d2val2_d1val3
      q82 = 4.0_dp*q33
      q83 = pow2(x%d1val2_d1val3)
      q84 = x%d1val2*x%d1val2_d2val3
      q85 = q13*x%d2val2 + q59*x%d2val2_d2val3 + q65*x%d2val2_d1val3
      q86 = 120.0_dp*x%d1val1
      q87 = 8.0_dp*q20
      q88 = 12.0_dp*q15 + 3.0_dp*q14
      q89 = 8.0_dp*x%d3val1
      q90 = 24.0_dp*q37
      q91 = 120.0_dp*x%d1val1_d1val3
      q92 = 16.0_dp*q37
      q93 = 0.375_dp*x%val
      q94 = q63*x%d1val2
      q95 = 16.0_dp*q20
      q96 = 40.0_dp*x%d1val1_d1val2
      q97 = 4.0_dp*x%d2val1_d1val2*x%val + q62*x%d1val2 + q96*x%d1val1
      q98 = 2.0_dp*x%d2val1
      q99 = 80.0_dp*x%d1val2
      q100 = 16.0_dp*q10 + q82
      q101 = q96*x%d1val2
      q102 = q96*x%val
      q103 = 40.0_dp*x%d1val2
      q104 = 8.0_dp*x%d1val1_d2val2
      q105 = 5.0_dp*q100
      q106 = 16.0_dp*q46 + q81
      q107 = 80.0_dp*x%d1val2_d1val3
      q108 = 120.0_dp*x%d1val2
      q109 = 60.0_dp*x%d1val2
      q110 = 8.0_dp*x%d3val2
      q111 = 120.0_dp*x%d1val2_d1val3
      unary%val = pow6(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q4*q6
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q7*q8
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q9*x%d1val1
      unary%d2val2 = q12*q6
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q9*x%d1val2
      unary%d2val3 = 3.0_dp*q5*(q14 + q16)
      unary%d3val1 = q22*q24
      unary%d2val1_d1val2 = q26*x%d1val2 + q29*q6
      unary%d2val1_d1val3 = q26*x%d1val3 + q30*q6
      unary%d1val1_d2val2 = q24*q35
      unary%d1val1_d1val2_d1val3 = 120.0_dp*q23*q8*x%d1val3 + q0*x%d1val1_d1val2_d1val3 + q7*x%d1val1*x%d1val2_d1val3 + q7*x%d1val1_d1val3*x%d1val2 + q9*x%d1val1_d1val2
      unary%d1val1_d2val3 = q41*(q36*q37 + q38*x%d1val1_d2val3 + q40*x%d1val1)
      unary%d3val2 = q24*q44
      unary%d2val2_d1val3 = q12*q45 + q6*(480.0_dp*q46 + q1*x%d2val2_d1val3 + q47*x%d2val2)
      unary%d1val2_d2val3 = q41*(40.0_dp*q48 + q38*x%d1val2_d2val3 + q40*x%d1val2)
      unary%d3val1_d1val3 = q22*q50 + q24*(17280.0_dp*q2*x%d1val1_d1val3 + 4320.0_dp*q53 + q19*q51 + q19*q52 + q21*x%d3val1_d1val3 + q54*x%d3val1)
      unary%d2val1_d1val2_d1val3 = 1.5_dp*q4*q49*x%d1val2 + q25*q30*x%d1val2 + q26*x%d1val2_d1val3 + q29*q45 + q6*(48.0_dp*q56 + 480.0_dp*q55 + q1*x%d2val1_d1val2_d1val3 + q27*x%d1val1_d1val2_d1val3 + q28*x%d1val2_d1val3 + q47*x%d2val1_d1val2)
      unary%d2val1_d2val3 = q69*(q60*q61 + q63*q64 + q67*q68)
      unary%d1val1_d2val2_d1val3 = q24*(2880.0_dp*q31*x%d1val2_d1val3 + 4.0_dp*x%d1val1*(2880.0_dp*q46 + 360.0_dp*q72 + 360.0_dp*q73) + 4.0_dp*x%d1val1_d1val3*(1440.0_dp*q10 + 360.0_dp*q33) + q21*x%d1val1_d2val2_d1val3 + q32*q70 + q32*q71 + q54*x%d1val1_d2val2) + q35*q50
      unary%d1val1_d1val2_d2val3 = q69*(30.0_dp*q15*q8 + 5.0_dp*q49*x%d1val1_d1val2_d1val3 + 5.0_dp*q79*x%d1val1_d1val3 + q16*q31 + q25*x%d1val1_d1val2_d2val3 + q48*q57 + q74*q75*x%d1val2 + q76*q77 + q78*x%d1val1*x%d1val2_d2val3 + q78*x%d1val1_d1val2*x%d2val3 + q78*x%d1val1_d2val3*x%d1val2)
      unary%d3val2_d1val3 = q24*(17280.0_dp*q10*x%d1val2_d1val3 + 4320.0_dp*q80 + q21*x%d3val2_d1val3 + q43*q72 + q43*q73 + q54*x%d3val2) + q44*q50
      unary%d2val2_d2val3 = q69*(q61*(20.0_dp*q46 + q81) + q64*(20.0_dp*q10 + q82) + q68*(20.0_dp*q83 + 20.0_dp*q84 + q85))
      unary%d3val1_d2val3 = q93*(q38*(60.0_dp*q18*x%d1val1_d2val3 + q27*q66 + q3*x%d1val1_d2val3 + q34*x%d2val1*x%d2val3 + q34*x%d2val1_d2val3*x%val + q38*x%d3val1_d2val3 + q39*x%d3val1 + q51*q91 + q52*q91 + q75*q89 + q86*x%d1val3*x%d2val1_d1val3 + q92*x%d3val1_d1val3) + q88*(160.0_dp*q17 + q18*q86 + q87*x%d3val1) + q90*(60.0_dp*q53 + q3*x%d1val1_d1val3 + q34*q51 + q34*q52 + q37*q89 + q38*x%d3val1_d1val3))
      unary%d2val1_d1val2_d2val3 = q93*(192.0_dp*q60*q77 + 32.0_dp*q20*q67*x%d1val2 + 64.0_dp*q49*(2.0_dp*q56 + 20.0_dp*q55 + q57*x%d1val1_d1val2_d1val3 + q58*x%d2val1_d1val2 + q59*x%d2val1_d1val2_d1val3 + q98*x%d1val2_d1val3) + 64.0_dp*q60*q79 + 8.0_dp*q23*(2.0_dp*x%d1val2*x%d2val1_d2val3 + 20.0_dp*x%d1val1_d1val2*x%d1val1_d2val3 + 4.0_dp*x%d1val2_d1val3*x%d2val1_d1val3 + q13*x%d2val1_d1val2 + q36*x%d1val1_d1val2_d1val3 + q57*x%d1val1_d1val2_d2val3 + q59*x%d2val1_d1val2_d2val3 + q65*x%d2val1_d1val2_d1val3 + q98*x%d1val2_d2val3) + 96.0_dp*q15*q94 + 96.0_dp*q48*q63 + q1*q15*q97 + q1*q94*x%d2val3 + q63*q95*x%d1val2_d2val3 + q95*q97*x%d2val3)
      unary%d1val1_d2val2_d2val3 = q93*(q38*(q101*x%d2val3 + q102*x%d1val2_d2val3 + q103*x%d1val1_d1val2_d2val3*x%val + q104*q75 + q105*x%d1val1_d2val3 + q106*q76 + q107*q70 + q107*q71 + q38*x%d1val1_d2val2_d2val3 + q39*x%d1val1_d2val2 + q74*(16.0_dp*q83 + 16.0_dp*q84 + q85) + q92*x%d1val1_d2val2_d1val3 + q99*x%d1val1_d1val2_d1val3*x%d1val3) + q88*(q100*q74 + q31*q99 + q87*x%d1val1_d2val2) + q90*(q101*x%d1val3 + q102*x%d1val2_d1val3 + q103*q71 + q104*q37 + q105*x%d1val1_d1val3 + q106*q74 + q38*x%d1val1_d2val2_d1val3))
      unary%d3val2_d2val3 = q93*(q38*(480.0_dp*q83*x%d1val2 + 60.0_dp*q33*x%d1val2_d2val3 + q108*x%d1val3*x%d2val2_d1val3 + q109*x%d2val2*x%d2val3 + q109*x%d2val2_d2val3*x%val + q11*x%d1val2_d2val3 + q110*q75 + q111*q72 + q111*q73 + q38*x%d3val2_d2val3 + q39*x%d3val2 + q92*x%d3val2_d1val3) + q88*(160.0_dp*q42 + q108*q33 + q87*x%d3val2) + q90*(60.0_dp*q80 + q109*q72 + q109*q73 + q11*x%d1val2_d1val3 + q110*q37 + q38*x%d3val2_d1val3))
   end function pow6_self
   
   function pow7_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = 7.0_dp*pow6(x%val)
      q1 = 24.0_dp*x%val
      q2 = pow2(x%d1val1)
      q3 = 144.0_dp*q2 + q1*x%d2val1
      q4 = pow5(x%val)
      q5 = 0.2916666666666667_dp*q4
      q6 = 42.0_dp*q4
      q7 = x%d1val1*x%d1val2
      q8 = q6*x%d1val3
      q9 = pow2(x%d1val2)
      q10 = 144.0_dp*q9 + q1*x%d2val2
      q11 = 2.0_dp*x%d2val3
      q12 = q11*x%val
      q13 = pow2(x%d1val3)
      q14 = pow3(x%d1val1)
      q15 = x%d1val1*x%d2val1
      q16 = 5184.0_dp*q15
      q17 = pow2(x%val)
      q18 = 288.0_dp*q17
      q19 = 8640.0_dp*q14 + q16*x%val + q18*x%d3val1
      q20 = pow4(x%val)
      q21 = 0.024305555555555556_dp*q20
      q22 = 1.4583333333333333_dp*q20
      q23 = q3*x%d1val2
      q24 = 288.0_dp*x%d1val1
      q25 = 24.0_dp*x%d2val1
      q26 = q1*x%d2val1_d1val2 + q24*x%d1val1_d1val2 + q25*x%d1val2
      q27 = q22*x%d1val3
      q28 = q1*x%d2val1_d1val3 + q24*x%d1val1_d1val3 + q25*x%d1val3
      q29 = x%d1val1_d1val2*x%val
      q30 = 3456.0_dp*x%d1val2
      q31 = x%d2val2*x%val
      q32 = 36.0_dp*x%d1val1*(240.0_dp*q9 + 48.0_dp*q31) + q18*x%d1val1_d2val2 + q29*q30
      q33 = x%d1val1*x%d1val2_d1val3
      q34 = x%d1val1_d1val3*x%d1val2
      q35 = 48.0_dp*x%d1val1_d1val3
      q36 = x%d1val3*x%val
      q37 = 4.0_dp*q17
      q38 = 12.0_dp*q12 + 120.0_dp*q13
      q39 = 1.75_dp*q20
      q40 = pow3(x%d1val2)
      q41 = 5184.0_dp*x%d1val2
      q42 = 8640.0_dp*q40 + q18*x%d3val2 + q31*q41
      q43 = x%d1val2*x%d1val2_d1val3
      q44 = 24.0_dp*x%d1val3
      q45 = 48.0_dp*x%d1val2_d1val3
      q46 = pow3(x%val)
      q47 = q46*x%d1val3
      q48 = 0.09722222222222222_dp*q47
      q49 = 5184.0_dp*x%val
      q50 = x%d1val1*x%d2val1_d1val3
      q51 = x%d1val1_d1val3*x%d2val1
      q52 = 576.0_dp*q36
      q53 = q2*x%d1val1_d1val3
      q54 = x%d1val1_d1val2*x%d1val1_d1val3
      q55 = 24.0_dp*x%d1val2
      q56 = 24.0_dp*x%d1val1
      q57 = 2.0_dp*x%d1val3
      q58 = 2.0_dp*x%val
      q59 = q56*x%d1val1_d1val3 + q57*x%d2val1 + q58*x%d2val1_d1val3
      q60 = 40.0_dp*q36
      q61 = 4.0_dp*x%d2val1
      q62 = 24.0_dp*q2 + q61*x%val
      q63 = 8.0_dp*q13
      q64 = 5.0_dp*q12 + 5.0_dp*q63
      q65 = 4.0_dp*x%d1val3
      q66 = pow2(x%d1val1_d1val3)
      q67 = 24.0_dp*q66 + q11*x%d2val1 + q56*x%d1val1_d2val3 + q58*x%d2val1_d2val3 + q65*x%d2val1_d1val3
      q68 = 0.875_dp*q46
      q69 = x%d1val1_d1val2*x%d1val3
      q70 = x%d1val1_d1val2_d1val3*x%val
      q71 = x%d1val3*x%d2val2
      q72 = x%d2val2_d1val3*x%val
      q73 = 8.0_dp*q46
      q74 = x%d2val3*x%val
      q75 = 480.0_dp*q36
      q76 = 48.0_dp*x%d1val1
      q77 = 240.0_dp*q13
      q78 = 48.0_dp*q17
      q79 = x%d1val1_d1val2*x%d2val3
      q80 = 96.0_dp*q17
      q81 = x%d1val1_d1val2_d1val3*x%d1val3
      q82 = q31*x%d1val2_d1val3
      q83 = q9*x%d1val2_d1val3
      q84 = q57*x%d2val2 + q58*x%d2val2_d1val3
      q85 = 4.0_dp*q31
      q86 = pow2(x%d1val2_d1val3)
      q87 = q11*x%d2val2 + q58*x%d2val2_d2val3 + q65*x%d2val2_d1val3
      q88 = 144.0_dp*x%val
      q89 = 8.0_dp*q17
      q90 = 24.0_dp*q13 + 4.0_dp*q12
      q91 = 72.0_dp*q15
      q92 = 72.0_dp*x%val
      q93 = 8.0_dp*x%d3val1
      q94 = 32.0_dp*q36
      q95 = 144.0_dp*x%d1val3
      q96 = 16.0_dp*q36
      q97 = 0.4375_dp*q17
      q98 = q62*x%d1val2
      q99 = 20.0_dp*x%d1val2_d2val3
      q100 = 80.0_dp*q17
      q101 = 4.0_dp*x%d2val1_d1val2*x%val + q61*x%d1val2 + q76*x%d1val1_d1val2
      q102 = 2.0_dp*x%d1val2
      q103 = 2.0_dp*x%d2val1
      q104 = 96.0_dp*x%d1val2
      q105 = 20.0_dp*q9 + q85
      q106 = 12.0_dp*x%d1val1
      q107 = 48.0_dp*x%d1val2
      q108 = 8.0_dp*x%d1val1_d2val2
      q109 = 6.0_dp*q105
      q110 = 20.0_dp*q43 + q84
      q111 = 96.0_dp*x%d1val2_d1val3
      q112 = 72.0_dp*x%d1val2
      q113 = 8.0_dp*x%d3val2
      q114 = 144.0_dp*x%d1val2_d1val3
      unary%val = pow7(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q3*q5
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q6*q7
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q8*x%d1val1
      unary%d2val2 = q10*q5
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q8*x%d1val2
      unary%d2val3 = 3.5_dp*q4*(12.0_dp*q13 + q12)
      unary%d3val1 = q19*q21
      unary%d2val1_d1val2 = q22*q23 + q26*q5
      unary%d2val1_d1val3 = q27*q3 + q28*q5
      unary%d1val1_d2val2 = q21*q32
      unary%d1val1_d1val2_d1val3 = 210.0_dp*q20*q7*x%d1val3 + q0*x%d1val1_d1val2_d1val3 + q33*q6 + q34*q6 + q8*x%d1val1_d1val2
      unary%d1val1_d2val3 = q39*(q35*q36 + q37*x%d1val1_d2val3 + q38*x%d1val1)
      unary%d3val2 = q21*q42
      unary%d2val2_d1val3 = q10*q27 + q5*(288.0_dp*q43 + q1*x%d2val2_d1val3 + q44*x%d2val2)
      unary%d1val2_d2val3 = q39*(q36*q45 + q37*x%d1val2_d2val3 + q38*x%d1val2)
      unary%d3val1_d1val3 = q19*q48 + q21*(25920.0_dp*q53 + q16*x%d1val3 + q18*x%d3val1_d1val3 + q49*q50 + q49*q51 + q52*x%d3val1)
      unary%d2val1_d1val2_d1val3 = 5.833333333333333_dp*q23*q47 + q22*q28*x%d1val2 + q22*q3*x%d1val2_d1val3 + q26*q27 + q5*(288.0_dp*q54 + q1*x%d2val1_d1val2_d1val3 + q24*x%d1val1_d1val2_d1val3 + q25*x%d1val2_d1val3 + q44*x%d2val1_d1val2 + q55*x%d2val1_d1val3)
      unary%d2val1_d2val3 = q68*(q37*q67 + q59*q60 + q62*q64)
      unary%d1val1_d2val2_d1val3 = q21*(3456.0_dp*q29*x%d1val2_d1val3 + 4.0_dp*x%d1val1*(432.0_dp*q71 + 432.0_dp*q72 + 4320.0_dp*q43) + 4.0_dp*x%d1val1_d1val3*(2160.0_dp*q9 + 432.0_dp*q31) + q18*x%d1val1_d2val2_d1val3 + q30*q69 + q30*q70 + q52*x%d1val1_d2val2) + q32*q48
      unary%d1val1_d1val2_d2val3 = q68*(240.0_dp*q7*q74 + 960.0_dp*q13*q7 + q17*q76*x%d1val2_d2val3 + q29*q77 + q33*q75 + q34*q75 + q73*x%d1val1_d1val2_d2val3 + q78*q79 + q78*x%d1val1_d2val3*x%d1val2 + q80*q81 + q80*x%d1val1_d1val3*x%d1val2_d1val3)
      unary%d3val2_d1val3 = q21*(25920.0_dp*q83 + 5184.0_dp*q82 + q18*x%d3val2_d1val3 + q41*q71 + q41*q72 + q52*x%d3val2) + q42*q48
      unary%d2val2_d2val3 = q68*(q37*(24.0_dp*q86 + q55*x%d1val2_d2val3 + q87) + q60*(24.0_dp*q43 + q84) + q64*(24.0_dp*q9 + q85))
      unary%d3val1_d2val3 = q97*(q37*(360.0_dp*q2*x%d1val1_d2val3 + 720.0_dp*q66*x%d1val1 + q37*x%d3val1_d2val3 + q50*q95 + q51*q95 + q63*x%d3val1 + q74*q93 + q88*x%d1val1_d1val3*x%d2val1_d1val3 + q91*x%d2val3 + q92*x%d1val1*x%d2val1_d2val3 + q92*x%d1val1_d2val3*x%d2val1 + q96*x%d3val1_d1val3) + q90*(240.0_dp*q14 + q15*q88 + q89*x%d3val1) + q94*(360.0_dp*q53 + q36*q93 + q37*x%d3val1_d1val3 + q50*q92 + q51*q92 + q91*x%d1val3))
      unary%d2val1_d1val2_d2val3 = q97*(160.0_dp*q36*q62*x%d1val2_d1val3 + 20.0_dp*q101*q17*x%d2val3 + 320.0_dp*q36*q59*x%d1val2 + 40.0_dp*q17*q67*x%d1val2 + 80.0_dp*q101*q13*x%val + 80.0_dp*q74*q98 + q100*q59*x%d1val2_d1val3 + q100*x%d1val3*(24.0_dp*q54 + q102*x%d2val1_d1val3 + q103*x%d1val2_d1val3 + q56*x%d1val1_d1val2_d1val3 + q57*x%d2val1_d1val2 + q58*x%d2val1_d1val2_d1val3) + q17*q62*q99 + q73*(24.0_dp*x%d1val1_d1val2*x%d1val1_d2val3 + 4.0_dp*x%d1val2_d1val3*x%d2val1_d1val3 + q102*x%d2val1_d2val3 + q103*x%d1val2_d2val3 + q11*x%d2val1_d1val2 + q35*x%d1val1_d1val2_d1val3 + q56*x%d1val1_d1val2_d2val3 + q58*x%d2val1_d1val2_d2val3 + q65*x%d2val1_d1val2_d1val3) + q77*q98)
      unary%d1val1_d2val2_d2val3 = q97*(q37*(24.0_dp*q110*x%d1val1_d1val3 + 48.0_dp*q29*x%d1val2_d2val3 + q104*q81 + q106*(20.0_dp*q86 + q87 + q99*x%d1val2) + q107*q79 + q107*x%d1val1_d1val2_d2val3*x%val + q108*q74 + q109*x%d1val1_d2val3 + q111*q69 + q111*q70 + q37*x%d1val1_d2val2_d2val3 + q63*x%d1val1_d2val2 + q96*x%d1val1_d2val2_d1val3) + q90*(q104*q29 + q105*q106 + q89*x%d1val1_d2val2) + q94*(q106*q110 + q107*q69 + q107*q70 + q108*q36 + q109*x%d1val1_d1val3 + q29*q45 + q37*x%d1val1_d2val2_d1val3))
      unary%d3val2_d2val3 = q97*(q37*(360.0_dp*q9*x%d1val2_d2val3 + 72.0_dp*q31*x%d1val2_d2val3 + 720.0_dp*q86*x%d1val2 + q112*x%d2val2*x%d2val3 + q113*q74 + q114*q71 + q114*q72 + q37*x%d3val2_d2val3 + q63*x%d3val2 + q92*x%d1val2*x%d2val2_d2val3 + q95*x%d1val2*x%d2val2_d1val3 + q96*x%d3val2_d1val3) + q90*(144.0_dp*q31*x%d1val2 + 240.0_dp*q40 + q89*x%d3val2) + q94*(360.0_dp*q83 + 72.0_dp*q82 + q112*q71 + q112*q72 + q113*q36 + q37*x%d3val2_d1val3))
   end function pow7_self
   
   function pow8_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = 8.0_dp*pow7(x%val)
      q1 = 48.0_dp*x%val
      q2 = pow2(x%d1val1)
      q3 = 336.0_dp*q2 + q1*x%d2val1
      q4 = pow6(x%val)
      q5 = 0.16666666666666666_dp*q4
      q6 = 56.0_dp*x%d1val1
      q7 = q4*q6
      q8 = pow2(x%d1val2)
      q9 = 336.0_dp*q8 + q1*x%d2val2
      q10 = q4*x%d1val2
      q11 = 56.0_dp*x%d1val3
      q12 = 2.0_dp*x%d2val3
      q13 = q12*x%val
      q14 = pow2(x%d1val3)
      q15 = pow3(x%d1val1)
      q16 = x%d2val1*x%val
      q17 = 6048.0_dp*x%d1val1
      q18 = pow2(x%val)
      q19 = 288.0_dp*q18
      q20 = 12096.0_dp*q15 + q16*q17 + q19*x%d3val1
      q21 = pow5(x%val)
      q22 = 0.027777777777777776_dp*q21
      q23 = q21*q3
      q24 = 672.0_dp*x%d1val1
      q25 = 48.0_dp*x%d2val1
      q26 = q1*x%d2val1_d1val2 + q24*x%d1val1_d1val2 + q25*x%d1val2
      q27 = q1*x%d2val1_d1val3 + q24*x%d1val1_d1val3 + q25*x%d1val3
      q28 = x%d1val1_d1val2*x%val
      q29 = 4032.0_dp*x%d1val2
      q30 = x%d2val2*x%val
      q31 = 84.0_dp*x%d1val1
      q32 = q19*x%d1val1_d2val2 + q28*q29 + q31*(144.0_dp*q8 + 24.0_dp*q30)
      q33 = q11*x%d1val1_d1val2
      q34 = 56.0_dp*x%d1val1_d1val3
      q35 = q21*x%d1val3
      q36 = 336.0_dp*x%d1val1*x%d1val2
      q37 = x%d1val3*x%val
      q38 = 4.0_dp*q18
      q39 = 14.0_dp*q13 + 168.0_dp*q14
      q40 = 2.0_dp*q21
      q41 = pow3(x%d1val2)
      q42 = 6048.0_dp*x%d1val2
      q43 = 12096.0_dp*q41 + q19*x%d3val2 + q30*q42
      q44 = x%d1val2*x%d1val2_d1val3
      q45 = 48.0_dp*x%d1val3
      q46 = x%d1val2_d1val3*x%val
      q47 = pow4(x%val)
      q48 = q47*x%d1val3
      q49 = 0.1388888888888889_dp*q48
      q50 = x%d1val3*x%d2val1
      q51 = x%d2val1_d1val3*x%val
      q52 = q16*x%d1val1_d1val3
      q53 = 576.0_dp*q37
      q54 = q2*x%d1val1_d1val3
      q55 = x%d1val1_d1val2*x%d1val1_d1val3
      q56 = x%d1val2*x%d2val1_d1val3
      q57 = 28.0_dp*x%d1val1
      q58 = 2.0_dp*x%d1val3
      q59 = 2.0_dp*x%val
      q60 = q57*x%d1val1_d1val3 + q58*x%d2val1 + q59*x%d2val1_d1val3
      q61 = q1*x%d1val3
      q62 = 4.0_dp*x%d2val1
      q63 = 28.0_dp*q2 + q62*x%val
      q64 = 6.0_dp*q13 + 60.0_dp*q14
      q65 = 4.0_dp*x%d1val3
      q66 = pow2(x%d1val1_d1val3)
      q67 = 28.0_dp*q66 + q12*x%d2val1 + q57*x%d1val1_d2val3 + q59*x%d2val1_d2val3 + q65*x%d2val1_d1val3
      q68 = x%d1val1_d1val2*x%d1val3
      q69 = q28*x%d1val2_d1val3
      q70 = x%d1val1_d1val2_d1val3*x%val
      q71 = x%d1val3*x%d2val2
      q72 = x%d2val2_d1val3*x%val
      q73 = pow3(x%val)
      q74 = 8.0_dp*q73
      q75 = x%d2val3*x%val
      q76 = q37*x%d1val2_d1val3
      q77 = q37*x%d1val2
      q78 = q14*x%d1val2
      q79 = q18*x%d1val2_d2val3
      q80 = 56.0_dp*q18
      q81 = x%d1val1_d1val2*x%d2val3
      q82 = 112.0_dp*q18
      q83 = x%d1val1_d1val2_d1val3*x%d1val3
      q84 = q30*x%d1val2_d1val3
      q85 = q8*x%d1val2_d1val3
      q86 = q58*x%d2val2 + q59*x%d2val2_d1val3
      q87 = 4.0_dp*q30
      q88 = pow2(x%d1val2_d1val3)
      q89 = x%d1val2*x%d1val2_d2val3
      q90 = q12*x%d2val2 + q59*x%d2val2_d2val3 + q65*x%d2val2_d1val3
      q91 = 168.0_dp*x%d1val1
      q92 = 8.0_dp*q18
      q93 = 8.0_dp*q14
      q94 = 5.0_dp*q13 + 5.0_dp*q93
      q95 = 8.0_dp*x%d3val1
      q96 = 40.0_dp*q37
      q97 = 168.0_dp*x%d1val1_d1val3
      q98 = 16.0_dp*q37
      q99 = 0.5_dp*q73
      q100 = 6.0_dp*q18
      q101 = 4.0_dp*x%d2val1_d1val2*x%val + q6*x%d1val1_d1val2 + q62*x%d1val2
      q102 = 2.0_dp*x%d2val1
      q103 = 112.0_dp*x%d1val2
      q104 = 24.0_dp*q8 + q87
      q105 = 14.0_dp*x%d1val1
      q106 = 56.0_dp*x%d1val2
      q107 = 8.0_dp*x%d1val1_d2val2
      q108 = 7.0_dp*q104
      q109 = 24.0_dp*q44 + q86
      q110 = 168.0_dp*x%d1val2
      q111 = 84.0_dp*x%d1val2
      q112 = 8.0_dp*x%d3val2
      q113 = 168.0_dp*x%d1val2_d1val3
      unary%val = pow8(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q3*q5
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q7*x%d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q7*x%d1val3
      unary%d2val2 = q5*q9
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q10*q11
      unary%d2val3 = 4.0_dp*q4*(14.0_dp*q14 + q13)
      unary%d3val1 = q20*q22
      unary%d2val1_d1val2 = q23*x%d1val2 + q26*q5
      unary%d2val1_d1val3 = q23*x%d1val3 + q27*q5
      unary%d1val1_d2val2 = q22*q32
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3 + q10*q34 + q33*q4 + q35*q36 + q7*x%d1val2_d1val3
      unary%d1val1_d2val3 = q40*(q34*q37 + q38*x%d1val1_d2val3 + q39*x%d1val1)
      unary%d3val2 = q22*q43
      unary%d2val2_d1val3 = q35*q9 + q5*(672.0_dp*q44 + q1*x%d2val2_d1val3 + q45*x%d2val2)
      unary%d1val2_d2val3 = q40*(q11*q46 + q38*x%d1val2_d2val3 + q39*x%d1val2)
      unary%d3val1_d1val3 = q20*q49 + q22*(36288.0_dp*q54 + 6048.0_dp*q52 + q17*q50 + q17*q51 + q19*x%d3val1_d1val3 + q53*x%d3val1)
      unary%d2val1_d1val2_d1val3 = 5.0_dp*q3*q48*x%d1val2 + q21*q27*x%d1val2 + q23*x%d1val2_d1val3 + q26*q35 + q5*(48.0_dp*q56 + 672.0_dp*q55 + q1*x%d2val1_d1val2_d1val3 + q24*x%d1val1_d1val2_d1val3 + q25*x%d1val2_d1val3 + q45*x%d2val1_d1val2)
      unary%d2val1_d2val3 = q47*(q38*q67 + q60*q61 + q63*q64)
      unary%d1val1_d2val2_d1val3 = q22*(4.0_dp*x%d1val1*(504.0_dp*q71 + 504.0_dp*q72 + 6048.0_dp*q44) + 4.0_dp*x%d1val1_d1val3*(3024.0_dp*q8 + 504.0_dp*q30) + 4032.0_dp*q69 + q19*x%d1val1_d2val2_d1val3 + q29*q68 + q29*q70 + q53*x%d1val1_d2val2) + q32*q49
      unary%d1val1_d1val2_d2val3 = q47*(1680.0_dp*q78*x%d1val1 + 336.0_dp*q14*q28 + 672.0_dp*q77*x%d1val1_d1val3 + q24*q76 + q36*q75 + q6*q79 + q74*x%d1val1_d1val2_d2val3 + q80*q81 + q80*x%d1val1_d2val3*x%d1val2 + q82*q83 + q82*x%d1val1_d1val3*x%d1val2_d1val3)
      unary%d3val2_d1val3 = q22*(36288.0_dp*q85 + 6048.0_dp*q84 + q19*x%d3val2_d1val3 + q42*q71 + q42*q72 + q53*x%d3val2) + q43*q49
      unary%d2val2_d2val3 = q47*(q38*(28.0_dp*q88 + 28.0_dp*q89 + q90) + q61*(28.0_dp*q44 + q86) + q64*(28.0_dp*q8 + q87))
      unary%d3val1_d2val3 = q99*(q38*(1008.0_dp*q66*x%d1val1 + 504.0_dp*q2*x%d1val1_d2val3 + 84.0_dp*q16*x%d1val1_d2val3 + q31*x%d2val1*x%d2val3 + q31*x%d2val1_d2val3*x%val + q38*x%d3val1_d2val3 + q50*q97 + q51*q97 + q75*q95 + q91*x%d1val3*x%d2val1_d1val3 + q93*x%d3val1 + q98*x%d3val1_d1val3) + q94*(336.0_dp*q15 + q16*q91 + q92*x%d3val1) + q96*(504.0_dp*q54 + 84.0_dp*q52 + q31*q50 + q31*q51 + q37*q95 + q38*x%d3val1_d1val3))
      unary%d2val1_d1val2_d2val3 = q74*(1.5_dp*q101*q18*x%d2val3 + 1.5_dp*q63*q79 + 15.0_dp*q63*q76 + 3.0_dp*q18*q67*x%d1val2 + 30.0_dp*q60*q77 + 30.0_dp*q63*q78 + 7.5_dp*q101*q14*x%val + 7.5_dp*q63*q75*x%d1val2 + q100*q60*x%d1val2_d1val3 + q100*x%d1val3*(2.0_dp*q56 + 28.0_dp*q55 + q102*x%d1val2_d1val3 + q57*x%d1val1_d1val2_d1val3 + q58*x%d2val1_d1val2 + q59*x%d2val1_d1val2_d1val3) + q99*(2.0_dp*x%d1val2*x%d2val1_d2val3 + 28.0_dp*x%d1val1_d1val2*x%d1val1_d2val3 + 4.0_dp*x%d1val2_d1val3*x%d2val1_d1val3 + q102*x%d1val2_d2val3 + q12*x%d2val1_d1val2 + q34*x%d1val1_d1val2_d1val3 + q57*x%d1val1_d1val2_d2val3 + q59*x%d2val1_d1val2_d2val3 + q65*x%d2val1_d1val2_d1val3))
      unary%d1val1_d2val2_d2val3 = q99*(q38*(112.0_dp*q46*x%d1val1_d1val2_d1val3 + 112.0_dp*q68*x%d1val2_d1val3 + 28.0_dp*q109*x%d1val1_d1val3 + 56.0_dp*q28*x%d1val2_d2val3 + q103*q83 + q105*(24.0_dp*q88 + 24.0_dp*q89 + q90) + q106*q81 + q106*x%d1val1_d1val2_d2val3*x%val + q107*q75 + q108*x%d1val1_d2val3 + q38*x%d1val1_d2val2_d2val3 + q93*x%d1val1_d2val2 + q98*x%d1val1_d2val2_d1val3) + q94*(q103*q28 + q104*q105 + q92*x%d1val1_d2val2) + q96*(56.0_dp*q69 + q105*q109 + q106*q70 + q107*q37 + q108*x%d1val1_d1val3 + q33*x%d1val2 + q38*x%d1val1_d2val2_d1val3))
      unary%d3val2_d2val3 = q99*(q38*(1008.0_dp*q88*x%d1val2 + 504.0_dp*q8*x%d1val2_d2val3 + 84.0_dp*q30*x%d1val2_d2val3 + q110*x%d1val3*x%d2val2_d1val3 + q111*x%d2val2*x%d2val3 + q111*x%d2val2_d2val3*x%val + q112*q75 + q113*q71 + q113*q72 + q38*x%d3val2_d2val3 + q93*x%d3val2 + q98*x%d3val2_d1val3) + q94*(336.0_dp*q41 + q110*q30 + q92*x%d3val2) + q96*(504.0_dp*q85 + 84.0_dp*q84 + q111*q71 + q111*q72 + q112*q37 + q38*x%d3val2_d1val3))
   end function pow8_self
   
   function abs_self(x) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q0
      q0 = sgn(x%val)
      unary%val = Abs(x%val)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function abs_self
   
   function add_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: binary
      binary%val = x%val + y%val
      binary%d1val1 = x%d1val1 + y%d1val1
      binary%d1val2 = x%d1val2 + y%d1val2
      binary%d1val3 = x%d1val3 + y%d1val3
      binary%d2val1 = x%d2val1 + y%d2val1
      binary%d1val1_d1val2 = x%d1val1_d1val2 + y%d1val1_d1val2
      binary%d1val1_d1val3 = x%d1val1_d1val3 + y%d1val1_d1val3
      binary%d2val2 = x%d2val2 + y%d2val2
      binary%d1val2_d1val3 = x%d1val2_d1val3 + y%d1val2_d1val3
      binary%d2val3 = x%d2val3 + y%d2val3
      binary%d3val1 = x%d3val1 + y%d3val1
      binary%d2val1_d1val2 = x%d2val1_d1val2 + y%d2val1_d1val2
      binary%d2val1_d1val3 = x%d2val1_d1val3 + y%d2val1_d1val3
      binary%d1val1_d2val2 = x%d1val1_d2val2 + y%d1val1_d2val2
      binary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3 + y%d1val1_d1val2_d1val3
      binary%d1val1_d2val3 = x%d1val1_d2val3 + y%d1val1_d2val3
      binary%d3val2 = x%d3val2 + y%d3val2
      binary%d2val2_d1val3 = x%d2val2_d1val3 + y%d2val2_d1val3
      binary%d1val2_d2val3 = x%d1val2_d2val3 + y%d1val2_d2val3
      binary%d3val1_d1val3 = x%d3val1_d1val3 + y%d3val1_d1val3
      binary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3 + y%d2val1_d1val2_d1val3
      binary%d2val1_d2val3 = x%d2val1_d2val3 + y%d2val1_d2val3
      binary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3 + y%d1val1_d2val2_d1val3
      binary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3 + y%d1val1_d1val2_d2val3
      binary%d3val2_d1val3 = x%d3val2_d1val3 + y%d3val2_d1val3
      binary%d2val2_d2val3 = x%d2val2_d2val3 + y%d2val2_d2val3
      binary%d3val1_d2val3 = x%d3val1_d2val3 + y%d3val1_d2val3
      binary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3 + y%d2val1_d1val2_d2val3
      binary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3 + y%d1val1_d2val2_d2val3
      binary%d3val2_d2val3 = x%d3val2_d2val3 + y%d3val2_d2val3
   end function add_self
   
   function add_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      unary%val = x%val + y
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1val3 = x%d1val3
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1val3 = x%d1val1_d1val3
      unary%d2val2 = x%d2val2
      unary%d1val2_d1val3 = x%d1val2_d1val3
      unary%d2val3 = x%d2val3
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1val3 = x%d2val1_d1val3
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = x%d1val1_d2val3
      unary%d3val2 = x%d3val2
      unary%d2val2_d1val3 = x%d2val2_d1val3
      unary%d1val2_d2val3 = x%d1val2_d2val3
      unary%d3val1_d1val3 = x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = x%d3val2_d1val3
      unary%d2val2_d2val3 = x%d2val2_d2val3
      unary%d3val1_d2val3 = x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = x%d3val2_d2val3
   end function add_self_real
   
   function add_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      unary%val = x%val + z
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1val3 = x%d1val3
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1val3 = x%d1val1_d1val3
      unary%d2val2 = x%d2val2
      unary%d1val2_d1val3 = x%d1val2_d1val3
      unary%d2val3 = x%d2val3
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1val3 = x%d2val1_d1val3
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = x%d1val1_d2val3
      unary%d3val2 = x%d3val2
      unary%d2val2_d1val3 = x%d2val2_d1val3
      unary%d1val2_d2val3 = x%d1val2_d2val3
      unary%d3val1_d1val3 = x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = x%d3val2_d1val3
      unary%d2val2_d2val3 = x%d2val2_d2val3
      unary%d3val1_d2val3 = x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = x%d3val2_d2val3
   end function add_real_self
   
   function add_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val + y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1val3 = x%d1val3
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1val3 = x%d1val1_d1val3
      unary%d2val2 = x%d2val2
      unary%d1val2_d1val3 = x%d1val2_d1val3
      unary%d2val3 = x%d2val3
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1val3 = x%d2val1_d1val3
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = x%d1val1_d2val3
      unary%d3val2 = x%d3val2
      unary%d2val2_d1val3 = x%d2val2_d1val3
      unary%d1val2_d2val3 = x%d1val2_d2val3
      unary%d3val1_d1val3 = x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = x%d3val2_d1val3
      unary%d2val2_d2val3 = x%d2val2_d2val3
      unary%d3val1_d2val3 = x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = x%d3val2_d2val3
   end function add_self_int
   
   function add_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = x%val + y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1val3 = x%d1val3
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1val3 = x%d1val1_d1val3
      unary%d2val2 = x%d2val2
      unary%d1val2_d1val3 = x%d1val2_d1val3
      unary%d2val3 = x%d2val3
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1val3 = x%d2val1_d1val3
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = x%d1val1_d2val3
      unary%d3val2 = x%d3val2
      unary%d2val2_d1val3 = x%d2val2_d1val3
      unary%d1val2_d2val3 = x%d1val2_d2val3
      unary%d3val1_d1val3 = x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = x%d3val2_d1val3
      unary%d2val2_d2val3 = x%d2val2_d2val3
      unary%d3val1_d2val3 = x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = x%d3val2_d2val3
   end function add_int_self
   
   function sub_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: binary
      binary%val = x%val - y%val
      binary%d1val1 = x%d1val1 - y%d1val1
      binary%d1val2 = x%d1val2 - y%d1val2
      binary%d1val3 = x%d1val3 - y%d1val3
      binary%d2val1 = x%d2val1 - y%d2val1
      binary%d1val1_d1val2 = x%d1val1_d1val2 - y%d1val1_d1val2
      binary%d1val1_d1val3 = x%d1val1_d1val3 - y%d1val1_d1val3
      binary%d2val2 = x%d2val2 - y%d2val2
      binary%d1val2_d1val3 = x%d1val2_d1val3 - y%d1val2_d1val3
      binary%d2val3 = x%d2val3 - y%d2val3
      binary%d3val1 = x%d3val1 - y%d3val1
      binary%d2val1_d1val2 = x%d2val1_d1val2 - y%d2val1_d1val2
      binary%d2val1_d1val3 = x%d2val1_d1val3 - y%d2val1_d1val3
      binary%d1val1_d2val2 = x%d1val1_d2val2 - y%d1val1_d2val2
      binary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3 - y%d1val1_d1val2_d1val3
      binary%d1val1_d2val3 = x%d1val1_d2val3 - y%d1val1_d2val3
      binary%d3val2 = x%d3val2 - y%d3val2
      binary%d2val2_d1val3 = x%d2val2_d1val3 - y%d2val2_d1val3
      binary%d1val2_d2val3 = x%d1val2_d2val3 - y%d1val2_d2val3
      binary%d3val1_d1val3 = x%d3val1_d1val3 - y%d3val1_d1val3
      binary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3 - y%d2val1_d1val2_d1val3
      binary%d2val1_d2val3 = x%d2val1_d2val3 - y%d2val1_d2val3
      binary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3 - y%d1val1_d2val2_d1val3
      binary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3 - y%d1val1_d1val2_d2val3
      binary%d3val2_d1val3 = x%d3val2_d1val3 - y%d3val2_d1val3
      binary%d2val2_d2val3 = x%d2val2_d2val3 - y%d2val2_d2val3
      binary%d3val1_d2val3 = x%d3val1_d2val3 - y%d3val1_d2val3
      binary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3 - y%d2val1_d1val2_d2val3
      binary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3 - y%d1val1_d2val2_d2val3
      binary%d3val2_d2val3 = x%d3val2_d2val3 - y%d3val2_d2val3
   end function sub_self
   
   function sub_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      unary%val = x%val - y
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1val3 = x%d1val3
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1val3 = x%d1val1_d1val3
      unary%d2val2 = x%d2val2
      unary%d1val2_d1val3 = x%d1val2_d1val3
      unary%d2val3 = x%d2val3
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1val3 = x%d2val1_d1val3
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = x%d1val1_d2val3
      unary%d3val2 = x%d3val2
      unary%d2val2_d1val3 = x%d2val2_d1val3
      unary%d1val2_d2val3 = x%d1val2_d2val3
      unary%d3val1_d1val3 = x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = x%d3val2_d1val3
      unary%d2val2_d2val3 = x%d2val2_d2val3
      unary%d3val1_d2val3 = x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = x%d3val2_d2val3
   end function sub_self_real
   
   function sub_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      unary%val = -x%val + z
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d1val3 = -x%d1val3
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d1val1_d1val3 = -x%d1val1_d1val3
      unary%d2val2 = -x%d2val2
      unary%d1val2_d1val3 = -x%d1val2_d1val3
      unary%d2val3 = -x%d2val3
      unary%d3val1 = -x%d3val1
      unary%d2val1_d1val2 = -x%d2val1_d1val2
      unary%d2val1_d1val3 = -x%d2val1_d1val3
      unary%d1val1_d2val2 = -x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = -x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = -x%d1val1_d2val3
      unary%d3val2 = -x%d3val2
      unary%d2val2_d1val3 = -x%d2val2_d1val3
      unary%d1val2_d2val3 = -x%d1val2_d2val3
      unary%d3val1_d1val3 = -x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = -x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = -x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = -x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = -x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = -x%d3val2_d1val3
      unary%d2val2_d2val3 = -x%d2val2_d2val3
      unary%d3val1_d2val3 = -x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = -x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = -x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = -x%d3val2_d2val3
   end function sub_real_self
   
   function sub_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val - y_dp
      unary%d1val1 = x%d1val1
      unary%d1val2 = x%d1val2
      unary%d1val3 = x%d1val3
      unary%d2val1 = x%d2val1
      unary%d1val1_d1val2 = x%d1val1_d1val2
      unary%d1val1_d1val3 = x%d1val1_d1val3
      unary%d2val2 = x%d2val2
      unary%d1val2_d1val3 = x%d1val2_d1val3
      unary%d2val3 = x%d2val3
      unary%d3val1 = x%d3val1
      unary%d2val1_d1val2 = x%d2val1_d1val2
      unary%d2val1_d1val3 = x%d2val1_d1val3
      unary%d1val1_d2val2 = x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = x%d1val1_d2val3
      unary%d3val2 = x%d3val2
      unary%d2val2_d1val3 = x%d2val2_d1val3
      unary%d1val2_d2val3 = x%d1val2_d2val3
      unary%d3val1_d1val3 = x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = x%d3val2_d1val3
      unary%d2val2_d2val3 = x%d2val2_d2val3
      unary%d3val1_d2val3 = x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = x%d3val2_d2val3
   end function sub_self_int
   
   function sub_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = -x%val + y_dp
      unary%d1val1 = -x%d1val1
      unary%d1val2 = -x%d1val2
      unary%d1val3 = -x%d1val3
      unary%d2val1 = -x%d2val1
      unary%d1val1_d1val2 = -x%d1val1_d1val2
      unary%d1val1_d1val3 = -x%d1val1_d1val3
      unary%d2val2 = -x%d2val2
      unary%d1val2_d1val3 = -x%d1val2_d1val3
      unary%d2val3 = -x%d2val3
      unary%d3val1 = -x%d3val1
      unary%d2val1_d1val2 = -x%d2val1_d1val2
      unary%d2val1_d1val3 = -x%d2val1_d1val3
      unary%d1val1_d2val2 = -x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = -x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = -x%d1val1_d2val3
      unary%d3val2 = -x%d3val2
      unary%d2val2_d1val3 = -x%d2val2_d1val3
      unary%d1val2_d2val3 = -x%d1val2_d2val3
      unary%d3val1_d1val3 = -x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = -x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = -x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = -x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = -x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = -x%d3val2_d1val3
      unary%d2val2_d2val3 = -x%d2val2_d2val3
      unary%d3val1_d2val3 = -x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = -x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = -x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = -x%d3val2_d2val3
   end function sub_int_self
   
   function mul_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: binary
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
      q0 = 2.0_dp*x%d1val1
      q1 = 2.0_dp*y%d1val2
      q2 = 2.0_dp*y%d1val3
      q3 = 3.0_dp*x%d1val1
      q4 = 3.0_dp*x%d2val1
      q5 = 2.0_dp*y%d1val1
      q6 = 2.0_dp*x%d1val2
      q7 = 2.0_dp*x%d1val3
      q8 = 3.0_dp*x%d1val2
      q9 = 3.0_dp*x%d2val2
      q10 = 3.0_dp*y%d2val1
      q11 = 3.0_dp*y%d1val1
      q12 = 2.0_dp*x%d1val1_d1val2
      q13 = 2.0_dp*y%d1val1_d1val2
      q14 = 4.0_dp*y%d1val1_d1val3
      q15 = 2.0_dp*y%d1val2_d1val3
      q16 = 2.0_dp*x%d1val2_d1val3
      q17 = 3.0_dp*y%d2val2
      q18 = 3.0_dp*y%d1val2
      q19 = 4.0_dp*y%d1val2_d1val3
      q20 = 4.0_dp*y%d1val1_d1val2_d1val3
      binary%val = x%val*y%val
      binary%d1val1 = x%d1val1*y%val + x%val*y%d1val1
      binary%d1val2 = x%d1val2*y%val + x%val*y%d1val2
      binary%d1val3 = x%d1val3*y%val + x%val*y%d1val3
      binary%d2val1 = q0*y%d1val1 + x%d2val1*y%val + x%val*y%d2val1
      binary%d1val1_d1val2 = x%d1val1*y%d1val2 + x%d1val1_d1val2*y%val + x%d1val2*y%d1val1 + x%val*y%d1val1_d1val2
      binary%d1val1_d1val3 = x%d1val1*y%d1val3 + x%d1val1_d1val3*y%val + x%d1val3*y%d1val1 + x%val*y%d1val1_d1val3
      binary%d2val2 = q1*x%d1val2 + x%d2val2*y%val + x%val*y%d2val2
      binary%d1val2_d1val3 = x%d1val2*y%d1val3 + x%d1val2_d1val3*y%val + x%d1val3*y%d1val2 + x%val*y%d1val2_d1val3
      binary%d2val3 = q2*x%d1val3 + x%d2val3*y%val + x%val*y%d2val3
      binary%d3val1 = q3*y%d2val1 + q4*y%d1val1 + x%d3val1*y%val + x%val*y%d3val1
      binary%d2val1_d1val2 = q0*y%d1val1_d1val2 + q5*x%d1val1_d1val2 + x%d1val2*y%d2val1 + x%d2val1*y%d1val2 + x%d2val1_d1val2*y%val + x%val*y%d2val1_d1val2
      binary%d2val1_d1val3 = q0*y%d1val1_d1val3 + q5*x%d1val1_d1val3 + x%d1val3*y%d2val1 + x%d2val1*y%d1val3 + x%d2val1_d1val3*y%val + x%val*y%d2val1_d1val3
      binary%d1val1_d2val2 = q1*x%d1val1_d1val2 + q6*y%d1val1_d1val2 + x%d1val1*y%d2val2 + x%d1val1_d2val2*y%val + x%d2val2*y%d1val1 + x%val*y%d1val1_d2val2
      binary%d1val1_d1val2_d1val3 = x%d1val1*y%d1val2_d1val3 + x%d1val1_d1val2*y%d1val3 + x%d1val1_d1val2_d1val3*y%val + x%d1val1_d1val3*y%d1val2 + x%d1val2*y%d1val1_d1val3 + x%d1val2_d1val3*y%d1val1 + x%d1val3*y%d1val1_d1val2 + x%val*y%d1val1_d1val2_d1val3
      binary%d1val1_d2val3 = q2*x%d1val1_d1val3 + q7*y%d1val1_d1val3 + x%d1val1*y%d2val3 + x%d1val1_d2val3*y%val + x%d2val3*y%d1val1 + x%val*y%d1val1_d2val3
      binary%d3val2 = q8*y%d2val2 + q9*y%d1val2 + x%d3val2*y%val + x%val*y%d3val2
      binary%d2val2_d1val3 = q1*x%d1val2_d1val3 + q6*y%d1val2_d1val3 + x%d1val3*y%d2val2 + x%d2val2*y%d1val3 + x%d2val2_d1val3*y%val + x%val*y%d2val2_d1val3
      binary%d1val2_d2val3 = q2*x%d1val2_d1val3 + q7*y%d1val2_d1val3 + x%d1val2*y%d2val3 + x%d1val2_d2val3*y%val + x%d2val3*y%d1val2 + x%val*y%d1val2_d2val3
      binary%d3val1_d1val3 = q10*x%d1val1_d1val3 + q11*x%d2val1_d1val3 + q3*y%d2val1_d1val3 + q4*y%d1val1_d1val3 + x%d1val3*y%d3val1 + x%d3val1*y%d1val3 + x%d3val1_d1val3*y%val + x%val*y%d3val1_d1val3
      binary%d2val1_d1val2_d1val3 = q0*y%d1val1_d1val2_d1val3 + q12*y%d1val1_d1val3 + q13*x%d1val1_d1val3 + q5*x%d1val1_d1val2_d1val3 + x%d1val2*y%d2val1_d1val3 + x%d1val2_d1val3*y%d2val1 + x%d1val3*y%d2val1_d1val2 + x%d2val1*y%d1val2_d1val3 + x%d2val1_d1val2*y%d1val3 + x%d2val1_d1val2_d1val3*y%val + x%d2val1_d1val3*y%d1val2 + x%val*y%d2val1_d1val2_d1val3
      binary%d2val1_d2val3 = q0*y%d1val1_d2val3 + q14*x%d1val1_d1val3 + q2*x%d2val1_d1val3 + q5*x%d1val1_d2val3 + q7*y%d2val1_d1val3 + x%d2val1*y%d2val3 + x%d2val1_d2val3*y%val + x%d2val3*y%d2val1 + x%val*y%d2val1_d2val3
      binary%d1val1_d2val2_d1val3 = q1*x%d1val1_d1val2_d1val3 + q12*y%d1val2_d1val3 + q13*x%d1val2_d1val3 + q6*y%d1val1_d1val2_d1val3 + x%d1val1*y%d2val2_d1val3 + x%d1val1_d1val3*y%d2val2 + x%d1val1_d2val2*y%d1val3 + x%d1val1_d2val2_d1val3*y%val + x%d1val3*y%d1val1_d2val2 + x%d2val2*y%d1val1_d1val3 + x%d2val2_d1val3*y%d1val1 + x%val*y%d1val1_d2val2_d1val3
      binary%d1val1_d1val2_d2val3 = q15*x%d1val1_d1val3 + q16*y%d1val1_d1val3 + q2*x%d1val1_d1val2_d1val3 + q7*y%d1val1_d1val2_d1val3 + x%d1val1*y%d1val2_d2val3 + x%d1val1_d1val2*y%d2val3 + x%d1val1_d1val2_d2val3*y%val + x%d1val1_d2val3*y%d1val2 + x%d1val2*y%d1val1_d2val3 + x%d1val2_d2val3*y%d1val1 + x%d2val3*y%d1val1_d1val2 + x%val*y%d1val1_d1val2_d2val3
      binary%d3val2_d1val3 = q17*x%d1val2_d1val3 + q18*x%d2val2_d1val3 + q8*y%d2val2_d1val3 + q9*y%d1val2_d1val3 + x%d1val3*y%d3val2 + x%d3val2*y%d1val3 + x%d3val2_d1val3*y%val + x%val*y%d3val2_d1val3
      binary%d2val2_d2val3 = q1*x%d1val2_d2val3 + q19*x%d1val2_d1val3 + q2*x%d2val2_d1val3 + q6*y%d1val2_d2val3 + q7*y%d2val2_d1val3 + x%d2val2*y%d2val3 + x%d2val2_d2val3*y%val + x%d2val3*y%d2val2 + x%val*y%d2val2_d2val3
      binary%d3val1_d2val3 = 6.0_dp*x%d1val1_d1val3*y%d2val1_d1val3 + 6.0_dp*x%d2val1_d1val3*y%d1val1_d1val3 + q10*x%d1val1_d2val3 + q11*x%d2val1_d2val3 + q2*x%d3val1_d1val3 + q3*y%d2val1_d2val3 + q4*y%d1val1_d2val3 + q7*y%d3val1_d1val3 + x%d2val3*y%d3val1 + x%d3val1*y%d2val3 + x%d3val1_d2val3*y%val + x%val*y%d3val1_d2val3
      binary%d2val1_d1val2_d2val3 = q0*y%d1val1_d1val2_d2val3 + q12*y%d1val1_d2val3 + q13*x%d1val1_d2val3 + q14*x%d1val1_d1val2_d1val3 + q15*x%d2val1_d1val3 + q16*y%d2val1_d1val3 + q2*x%d2val1_d1val2_d1val3 + q20*x%d1val1_d1val3 + q5*x%d1val1_d1val2_d2val3 + q7*y%d2val1_d1val2_d1val3 + x%d1val2*y%d2val1_d2val3 + x%d1val2_d2val3*y%d2val1 + x%d2val1*y%d1val2_d2val3 + x%d2val1_d1val2*y%d2val3 + x%d2val1_d1val2_d2val3*y%val + x%d2val1_d2val3*y%d1val2 + x%d2val3*y%d2val1_d1val2 + x%val*y%d2val1_d1val2_d2val3
      binary%d1val1_d2val2_d2val3 = 2.0_dp*x%d1val1_d1val3*y%d2val2_d1val3 + 2.0_dp*x%d2val2_d1val3*y%d1val1_d1val3 + q1*x%d1val1_d1val2_d2val3 + q12*y%d1val2_d2val3 + q13*x%d1val2_d2val3 + q19*x%d1val1_d1val2_d1val3 + q2*x%d1val1_d2val2_d1val3 + q20*x%d1val2_d1val3 + q6*y%d1val1_d1val2_d2val3 + q7*y%d1val1_d2val2_d1val3 + x%d1val1*y%d2val2_d2val3 + x%d1val1_d2val2*y%d2val3 + x%d1val1_d2val2_d2val3*y%val + x%d1val1_d2val3*y%d2val2 + x%d2val2*y%d1val1_d2val3 + x%d2val2_d2val3*y%d1val1 + x%d2val3*y%d1val1_d2val2 + x%val*y%d1val1_d2val2_d2val3
      binary%d3val2_d2val3 = 6.0_dp*x%d1val2_d1val3*y%d2val2_d1val3 + 6.0_dp*x%d2val2_d1val3*y%d1val2_d1val3 + q17*x%d1val2_d2val3 + q18*x%d2val2_d2val3 + q2*x%d3val2_d1val3 + q7*y%d3val2_d1val3 + q8*y%d2val2_d2val3 + q9*y%d1val2_d2val3 + x%d2val3*y%d3val2 + x%d3val2*y%d2val3 + x%d3val2_d2val3*y%val + x%val*y%d3val2_d2val3
   end function mul_self
   
   function mul_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      unary%val = x%val*y
      unary%d1val1 = x%d1val1*y
      unary%d1val2 = x%d1val2*y
      unary%d1val3 = x%d1val3*y
      unary%d2val1 = x%d2val1*y
      unary%d1val1_d1val2 = x%d1val1_d1val2*y
      unary%d1val1_d1val3 = x%d1val1_d1val3*y
      unary%d2val2 = x%d2val2*y
      unary%d1val2_d1val3 = x%d1val2_d1val3*y
      unary%d2val3 = x%d2val3*y
      unary%d3val1 = x%d3val1*y
      unary%d2val1_d1val2 = x%d2val1_d1val2*y
      unary%d2val1_d1val3 = x%d2val1_d1val3*y
      unary%d1val1_d2val2 = x%d1val1_d2val2*y
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3*y
      unary%d1val1_d2val3 = x%d1val1_d2val3*y
      unary%d3val2 = x%d3val2*y
      unary%d2val2_d1val3 = x%d2val2_d1val3*y
      unary%d1val2_d2val3 = x%d1val2_d2val3*y
      unary%d3val1_d1val3 = x%d3val1_d1val3*y
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3*y
      unary%d2val1_d2val3 = x%d2val1_d2val3*y
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3*y
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3*y
      unary%d3val2_d1val3 = x%d3val2_d1val3*y
      unary%d2val2_d2val3 = x%d2val2_d2val3*y
      unary%d3val1_d2val3 = x%d3val1_d2val3*y
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3*y
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3*y
      unary%d3val2_d2val3 = x%d3val2_d2val3*y
   end function mul_self_real
   
   function mul_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      unary%val = x%val*z
      unary%d1val1 = x%d1val1*z
      unary%d1val2 = x%d1val2*z
      unary%d1val3 = x%d1val3*z
      unary%d2val1 = x%d2val1*z
      unary%d1val1_d1val2 = x%d1val1_d1val2*z
      unary%d1val1_d1val3 = x%d1val1_d1val3*z
      unary%d2val2 = x%d2val2*z
      unary%d1val2_d1val3 = x%d1val2_d1val3*z
      unary%d2val3 = x%d2val3*z
      unary%d3val1 = x%d3val1*z
      unary%d2val1_d1val2 = x%d2val1_d1val2*z
      unary%d2val1_d1val3 = x%d2val1_d1val3*z
      unary%d1val1_d2val2 = x%d1val1_d2val2*z
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3*z
      unary%d1val1_d2val3 = x%d1val1_d2val3*z
      unary%d3val2 = x%d3val2*z
      unary%d2val2_d1val3 = x%d2val2_d1val3*z
      unary%d1val2_d2val3 = x%d1val2_d2val3*z
      unary%d3val1_d1val3 = x%d3val1_d1val3*z
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3*z
      unary%d2val1_d2val3 = x%d2val1_d2val3*z
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3*z
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3*z
      unary%d3val2_d1val3 = x%d3val2_d1val3*z
      unary%d2val2_d2val3 = x%d2val2_d2val3*z
      unary%d3val1_d2val3 = x%d3val1_d2val3*z
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3*z
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3*z
      unary%d3val2_d2val3 = x%d3val2_d2val3*z
   end function mul_real_self
   
   function mul_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      y_dp = y
      unary%val = x%val*y_dp
      unary%d1val1 = x%d1val1*y_dp
      unary%d1val2 = x%d1val2*y_dp
      unary%d1val3 = x%d1val3*y_dp
      unary%d2val1 = x%d2val1*y_dp
      unary%d1val1_d1val2 = x%d1val1_d1val2*y_dp
      unary%d1val1_d1val3 = x%d1val1_d1val3*y_dp
      unary%d2val2 = x%d2val2*y_dp
      unary%d1val2_d1val3 = x%d1val2_d1val3*y_dp
      unary%d2val3 = x%d2val3*y_dp
      unary%d3val1 = x%d3val1*y_dp
      unary%d2val1_d1val2 = x%d2val1_d1val2*y_dp
      unary%d2val1_d1val3 = x%d2val1_d1val3*y_dp
      unary%d1val1_d2val2 = x%d1val1_d2val2*y_dp
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3*y_dp
      unary%d1val1_d2val3 = x%d1val1_d2val3*y_dp
      unary%d3val2 = x%d3val2*y_dp
      unary%d2val2_d1val3 = x%d2val2_d1val3*y_dp
      unary%d1val2_d2val3 = x%d1val2_d2val3*y_dp
      unary%d3val1_d1val3 = x%d3val1_d1val3*y_dp
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3*y_dp
      unary%d2val1_d2val3 = x%d2val1_d2val3*y_dp
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3*y_dp
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3*y_dp
      unary%d3val2_d1val3 = x%d3val2_d1val3*y_dp
      unary%d2val2_d2val3 = x%d2val2_d2val3*y_dp
      unary%d3val1_d2val3 = x%d3val1_d2val3*y_dp
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3*y_dp
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3*y_dp
      unary%d3val2_d2val3 = x%d3val2_d2val3*y_dp
   end function mul_self_int
   
   function mul_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      y_dp = z
      unary%val = x%val*y_dp
      unary%d1val1 = x%d1val1*y_dp
      unary%d1val2 = x%d1val2*y_dp
      unary%d1val3 = x%d1val3*y_dp
      unary%d2val1 = x%d2val1*y_dp
      unary%d1val1_d1val2 = x%d1val1_d1val2*y_dp
      unary%d1val1_d1val3 = x%d1val1_d1val3*y_dp
      unary%d2val2 = x%d2val2*y_dp
      unary%d1val2_d1val3 = x%d1val2_d1val3*y_dp
      unary%d2val3 = x%d2val3*y_dp
      unary%d3val1 = x%d3val1*y_dp
      unary%d2val1_d1val2 = x%d2val1_d1val2*y_dp
      unary%d2val1_d1val3 = x%d2val1_d1val3*y_dp
      unary%d1val1_d2val2 = x%d1val1_d2val2*y_dp
      unary%d1val1_d1val2_d1val3 = x%d1val1_d1val2_d1val3*y_dp
      unary%d1val1_d2val3 = x%d1val1_d2val3*y_dp
      unary%d3val2 = x%d3val2*y_dp
      unary%d2val2_d1val3 = x%d2val2_d1val3*y_dp
      unary%d1val2_d2val3 = x%d1val2_d2val3*y_dp
      unary%d3val1_d1val3 = x%d3val1_d1val3*y_dp
      unary%d2val1_d1val2_d1val3 = x%d2val1_d1val2_d1val3*y_dp
      unary%d2val1_d2val3 = x%d2val1_d2val3*y_dp
      unary%d1val1_d2val2_d1val3 = x%d1val1_d2val2_d1val3*y_dp
      unary%d1val1_d1val2_d2val3 = x%d1val1_d1val2_d2val3*y_dp
      unary%d3val2_d1val3 = x%d3val2_d1val3*y_dp
      unary%d2val2_d2val3 = x%d2val2_d2val3*y_dp
      unary%d3val1_d2val3 = x%d3val1_d2val3*y_dp
      unary%d2val1_d1val2_d2val3 = x%d2val1_d1val2_d2val3*y_dp
      unary%d1val1_d2val2_d2val3 = x%d1val1_d2val2_d2val3*y_dp
      unary%d3val2_d2val3 = x%d3val2_d2val3*y_dp
   end function mul_int_self
   
   function div_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: binary
      real(dp) :: q285
      real(dp) :: q284
      real(dp) :: q283
      real(dp) :: q282
      real(dp) :: q281
      real(dp) :: q280
      real(dp) :: q279
      real(dp) :: q278
      real(dp) :: q277
      real(dp) :: q276
      real(dp) :: q275
      real(dp) :: q274
      real(dp) :: q273
      real(dp) :: q272
      real(dp) :: q271
      real(dp) :: q270
      real(dp) :: q269
      real(dp) :: q268
      real(dp) :: q267
      real(dp) :: q266
      real(dp) :: q265
      real(dp) :: q264
      real(dp) :: q263
      real(dp) :: q262
      real(dp) :: q261
      real(dp) :: q260
      real(dp) :: q259
      real(dp) :: q258
      real(dp) :: q257
      real(dp) :: q256
      real(dp) :: q255
      real(dp) :: q254
      real(dp) :: q253
      real(dp) :: q252
      real(dp) :: q251
      real(dp) :: q250
      real(dp) :: q249
      real(dp) :: q248
      real(dp) :: q247
      real(dp) :: q246
      real(dp) :: q245
      real(dp) :: q244
      real(dp) :: q243
      real(dp) :: q242
      real(dp) :: q241
      real(dp) :: q240
      real(dp) :: q239
      real(dp) :: q238
      real(dp) :: q237
      real(dp) :: q236
      real(dp) :: q235
      real(dp) :: q234
      real(dp) :: q233
      real(dp) :: q232
      real(dp) :: q231
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = powm1(y%val)
      q1 = q0*x%val
      q2 = q0*x%d1val1
      q3 = powm1(pow2(y%val))
      q4 = q3*x%val
      q5 = q4*y%d1val1
      q6 = q0*x%d1val2
      q7 = q4*y%d1val2
      q8 = q0*x%d1val3
      q9 = q4*y%d1val3
      q10 = 4.0_dp*q2
      q11 = pow2(y%d1val1)
      q12 = 4.0_dp*q0
      q13 = 2.0_dp*y%d2val1 - q11*q12
      q14 = 2.0_dp*x%d2val1 - q1*q13 - q10*y%d1val1
      q15 = 0.5_dp*q0
      q16 = q14*q15
      q17 = q0*x%d1val1_d1val2
      q18 = q3*x%d1val1
      q19 = q18*y%d1val2
      q20 = q3*y%d1val1
      q21 = q20*x%d1val2
      q22 = q4*y%d1val1_d1val2
      q23 = 2.0_dp*y%d1val1
      q24 = powm1(pow3(y%val))
      q25 = q24*y%d1val2
      q26 = q25*x%val
      q27 = q0*x%d1val1_d1val3
      q28 = q18*y%d1val3
      q29 = q20*x%d1val3
      q30 = q4*y%d1val1_d1val3
      q31 = 2.0_dp*y%d1val3
      q32 = q24*y%d1val1
      q33 = q32*x%val
      q34 = 2.0_dp*x%d2val2
      q35 = 4.0_dp*y%d1val2
      q36 = 2.0_dp*y%d2val2
      q37 = pow2(y%d1val2)
      q38 = -q12*q37 + q36
      q39 = -q1*q38 + q34 - q35*q6
      q40 = q15*q39
      q41 = q0*x%d1val2_d1val3
      q42 = q3*y%d1val3
      q43 = q42*x%d1val2
      q44 = q3*y%d1val2
      q45 = q44*x%d1val3
      q46 = q4*y%d1val2_d1val3
      q47 = q26*q31
      q48 = pow2(y%d1val3)
      q49 = 2.0_dp*q48
      q50 = -q0*q49 + y%d2val3
      q51 = q0*y%d1val1
      q52 = 6.0_dp*x%d2val1
      q53 = 3.0_dp*q2
      q54 = 12.0_dp*q51
      q55 = pow3(y%d1val1)
      q56 = 12.0_dp*q3
      q57 = 2.0_dp*y%d3val1 - q54*y%d2val1 + q55*q56
      q58 = 2.0_dp*x%d3val1 - q1*q57 - q13*q53 - q51*q52
      q59 = q15*q58
      q60 = 0.5_dp*q14
      q61 = 4.0_dp*y%d1val1
      q62 = 8.0_dp*q51
      q63 = 4.0_dp*q11
      q64 = 2.0_dp*y%d2val1_d1val2 + q44*q63 - q62*y%d1val1_d1val2
      q65 = 2.0_dp*x%d2val1_d1val2 - q1*q64 - q10*y%d1val1_d1val2 - q13*q6 + q13*q7 - q17*q61 + q19*q61
      q66 = q15*q65
      q67 = q13*q8
      q68 = q13*q9
      q69 = 2.0_dp*y%d2val1_d1val3 + q42*q63 - q62*y%d1val1_d1val3
      q70 = 2.0_dp*x%d2val1_d1val3 - q1*q69 - q10*y%d1val1_d1val3 - q27*q61 + q28*q61 - q67 + q68
      q71 = 4.0_dp*q6
      q72 = 2.0_dp*y%d1val1_d2val2
      q73 = 8.0_dp*y%d1val2
      q74 = 4.0_dp*y%d2val2
      q75 = 8.0_dp*y%d1val1_d1val2
      q76 = 4.0_dp*q37
      q77 = 12.0_dp*q37
      q78 = 2.0_dp*x%d1val1_d2val2 - q1*q72 - q17*q35 + q18*q76 - q2*q36 + q21*q73 - q33*q77 - q34*q51 + q5*q74 + q7*q75 - q71*y%d1val1_d1val2
      q79 = q15*q78
      q80 = q0*x%d1val1_d1val2_d1val3
      q81 = q18*y%d1val2_d1val3
      q82 = q42*x%d1val1_d1val2
      q83 = q44*x%d1val1_d1val3
      q84 = q3*x%d1val2
      q85 = q84*y%d1val1_d1val3
      q86 = q20*x%d1val2_d1val3
      q87 = q3*x%d1val3
      q88 = q87*y%d1val1_d1val2
      q89 = q25*q31
      q90 = q25*x%d1val3
      q91 = 2.0_dp*y%d1val2_d1val3
      q92 = q24*x%val
      q93 = q31*q92
      q94 = 2.0_dp*y%d1val1_d1val3
      q95 = 6.0_dp*y%d1val2
      q96 = x%val*y%d1val1
      q97 = q95*q96
      q98 = powm1(pow4(y%val))
      q99 = q98*y%d1val3
      q100 = 2.0_dp*q8
      q101 = q0*x%d2val3
      q102 = 4.0_dp*y%d1val3
      q103 = 2.0_dp*y%d2val3
      q104 = 4.0_dp*y%d1val1_d1val3
      q105 = q24*q48
      q106 = q105*x%val
      q107 = 6.0_dp*q106
      q108 = q0*y%d1val2
      q109 = 6.0_dp*q108
      q110 = 3.0_dp*q6
      q111 = 12.0_dp*q108
      q112 = pow3(y%d1val2)
      q113 = 2.0_dp*y%d3val2 - q111*y%d2val2 + q112*q56
      q114 = 2.0_dp*x%d3val2 - q1*q113 - q109*x%d2val2 - q110*q38
      q115 = q114*q15
      q116 = 0.5_dp*q42
      q117 = 2.0_dp*x%d2val2_d1val3
      q118 = q38*q8
      q119 = q38*q9
      q120 = 2.0_dp*y%d2val2_d1val3
      q121 = 8.0_dp*y%d1val2_d1val3
      q122 = -q108*q121 + q120 + q42*q76
      q123 = q42*x%d1val3
      q124 = 4.0_dp*y%d1val2_d1val3
      q125 = q0*y%d1val1_d1val3
      q126 = 6.0_dp*q51
      q127 = q20*y%d1val3
      q128 = q13*q27
      q129 = 3.0_dp*q13
      q130 = q57*q8
      q131 = q57*q9
      q132 = 12.0_dp*y%d2val1
      q133 = q11*q3
      q134 = 36.0_dp*y%d1val1_d1val3
      q135 = q24*y%d1val3
      q136 = 24.0_dp*q135
      q137 = q3*y%d1val2_d1val3
      q138 = q135*y%d1val2
      q139 = 0.5_dp*q44
      q140 = 4.0_dp*y%d1val1_d1val2
      q141 = q102*q20
      q142 = q35*x%d1val1_d1val3
      q143 = 8.0_dp*q138
      q144 = q143*x%d1val1
      q145 = q13*q41
      q146 = q13*q43
      q147 = q13*q45
      q148 = q13*q46
      q149 = q64*q8
      q150 = q64*q9
      q151 = 8.0_dp*q125
      q152 = q20*y%d1val1_d1val3
      q153 = 2.0_dp*q2
      q154 = 2.0_dp*q27
      q155 = q12*y%d1val1
      q156 = 2.0_dp*q42
      q157 = q11*q156 - q155*y%d1val1_d1val3 + y%d2val1_d1val3
      q158 = -0.5_dp*q67 + 0.5_dp*q68 - q1*q157 - q153*y%d1val1_d1val3 - q154*y%d1val1 + q23*q28 + x%d2val1_d1val3
      q159 = q0*q31
      q160 = 2.0_dp*q51
      q161 = q103*q18
      q162 = q105*x%d1val1
      q163 = q15*x%d2val3
      q164 = q4*y%d2val3
      q165 = 0.5_dp*q164
      q166 = 2.0_dp*q9
      q167 = pow2(y%d1val1_d1val3)
      q168 = 8.0_dp*y%d1val1_d1val3
      q169 = -q103*q133 + q105*q63 + q12*q167 - q127*q168 + q155*y%d1val1_d2val3 - y%d2val1_d2val3
      q170 = q1*q169 - q100*q157 - q104*q27 + q104*q28 - q106*q13 + q123*q13 - q13*q163 + q13*q165 + q141*x%d1val1_d1val3 - q153*y%d1val1_d2val3 + q157*q166 - q160*x%d1val1_d2val3 + q161*y%d1val1 - q162*q61 + x%d2val1_d2val3
      q171 = q8*y%d1val1_d2val2
      q172 = q1*y%d1val1_d2val2_d1val3
      q173 = q44*x%d1val2
      q174 = 4.0_dp*y%d2val2_d1val3
      q175 = 8.0_dp*y%d1val1_d1val2_d1val3
      q176 = q3*x%d1val1_d1val3
      q177 = 16.0_dp*q138
      q178 = x%d1val2*y%d1val1
      q179 = q26*y%d1val1
      q180 = 24.0_dp*y%d1val2_d1val3
      q181 = q135*q96
      q182 = 8.0_dp*y%d2val2
      q183 = q177*y%d1val1_d1val2
      q184 = q135*x%d1val1
      q185 = 8.0_dp*q37
      q186 = q32*x%d1val3
      q187 = q92*y%d1val1_d1val3
      q188 = q37*q96
      q189 = q188*q99
      q190 = 2.0_dp*q41
      q191 = 2.0_dp*y%d1val2
      q192 = q20*x%d2val3
      q193 = 2.0_dp*q5
      q194 = 4.0_dp*y%d1val1_d1val2_d1val3
      q195 = 2.0_dp*y%d1val1_d2val3
      q196 = 12.0_dp*q138
      q197 = x%d1val3*y%d1val1
      q198 = q24*y%d2val3
      q199 = 12.0_dp*y%d1val2_d1val3
      q200 = q105*x%d1val2
      q201 = q48*q98
      q202 = q201*y%d1val2
      q203 = 24.0_dp*q202
      q204 = q0*x%d2val2
      q205 = 6.0_dp*y%d1val2_d1val3
      q206 = q42*q95
      q207 = q38*q41
      q208 = 3.0_dp*q38
      q209 = q113*q8
      q210 = q113*q9
      q211 = q0*y%d2val2
      q212 = q42*y%d1val2
      q213 = 12.0_dp*y%d2val2
      q214 = q137*q37
      q215 = 2.0_dp*q6
      q216 = 2.0_dp*q108
      q217 = q42*x%d1val2_d1val3
      q218 = q12*y%d1val2
      q219 = q156*q37 - q218*y%d1val2_d1val3 + y%d2val2_d1val3
      q220 = pow2(y%d1val2_d1val3)
      q221 = q3*q37
      q222 = -q103*q221 + q105*q76 + q12*q220 - q121*q212 + q218*y%d1val2_d2val3 - y%d2val2_d2val3
      q223 = q0*y%d1val1_d2val3
      q224 = 3.0_dp*x%d2val1
      q225 = 6.0_dp*x%d2val1_d1val3
      q226 = 3.0_dp*q51
      q227 = q20*y%d2val3
      q228 = q42*y%d1val1_d1val3
      q229 = q105*y%d1val1
      q230 = q0*x%d1val1_d2val3
      q231 = 1.5_dp*q13
      q232 = q18*y%d2val3
      q233 = q42*x%d1val1_d1val3
      q234 = 6.0_dp*q157
      q235 = 6.0_dp*y%d2val1
      q236 = 18.0_dp*q133
      q237 = 12.0_dp*q55
      q238 = -q125*q235 - q126*y%d2val1_d1val3 + q127*q235 - q135*q237 + q236*y%d1val1_d1val3 + y%d3val1_d1val3
      q239 = 12.0_dp*y%d2val1_d1val3
      q240 = q11*q135
      q241 = 36.0_dp*q201
      q242 = 2.0_dp*q17
      q243 = 2.0_dp*y%d1val1_d1val2
      q244 = q20*x%d1val1_d1val2
      q245 = q20*x%d1val1_d1val3
      q246 = 2.0_dp*x%d1val1_d2val3
      q247 = q198*q35
      q248 = x%d1val1*y%d1val1
      q249 = q184*y%d1val1
      q250 = q105*x%d1val1_d1val2
      q251 = 12.0_dp*q202
      q252 = q13*y%d1val2_d1val3
      q253 = 0.5_dp*q13
      q254 = q84*y%d2val3
      q255 = q13*x%val
      q256 = 2.0_dp*q157
      q257 = q135*x%val
      q258 = q44*y%d2val3
      q259 = q105*y%d1val2
      q260 = q12*y%d1val1_d1val2
      q261 = 2.0_dp*q133
      q262 = q127*q140 + q152*q35 - q155*y%d1val1_d1val2_d1val3 - q240*q35 - q260*y%d1val1_d1val3 + q261*y%d1val2_d1val3 + y%d2val1_d1val2_d1val3
      q263 = q35*q42
      q264 = q20*q35
      q265 = 4.0_dp*q44
      q266 = q177*y%d1val1_d1val3
      q267 = 2.0_dp*q80
      q268 = q0*x%d1val2_d2val3
      q269 = q20*x%d2val2
      q270 = 4.0_dp*y%d1val2_d2val3
      q271 = q198*q73
      q272 = q135*y%d1val2_d1val3
      q273 = 16.0_dp*q272
      q274 = x%val*y%d1val1_d1val2
      q275 = q257*y%d2val2
      q276 = 6.0_dp*q37
      q277 = q37*q99
      q278 = 72.0_dp*q37
      q279 = 3.0_dp*q204
      q280 = 3.0_dp*q108
      q281 = 3.0_dp*x%d2val2
      q282 = 1.5_dp*q38
      q283 = 6.0_dp*q219
      q284 = 12.0_dp*q112
      q285 = 18.0_dp*q214 - q109*y%d2val2_d1val3 - q135*q284 - q205*q211 + q206*y%d2val2 + y%d3val2_d1val3
      binary%val = q1
      binary%d1val1 = q2 - q5
      binary%d1val2 = q6 - q7
      binary%d1val3 = q8 - q9
      binary%d2val1 = q16
      binary%d1val1_d1val2 = q17 - q19 - q21 - q22 + q23*q26
      binary%d1val1_d1val3 = q27 - q28 - q29 - q30 + q31*q33
      binary%d2val2 = q40
      binary%d1val2_d1val3 = q41 - q43 - q45 - q46 + q47
      binary%d2val3 = q0*(-q1*q50 - q31*q8 + x%d2val3)
      binary%d3val1 = q59
      binary%d2val1_d1val2 = -q44*q60 + q66
      binary%d2val1_d1val3 = q15*q70 - q42*q60
      binary%d1val1_d2val2 = q79
      binary%d1val1_d1val2_d1val3 = q23*q90 + q26*q94 + q31*q32*x%d1val2 + q33*q91 - q4*y%d1val1_d1val2_d1val3 + q80 - q81 - q82 - q83 - q85 - q86 - q88 + q89*x%d1val1 + q93*y%d1val1_d1val2 - q97*q99
      binary%d1val1_d2val3 = q0*(-q1*y%d1val1_d2val3 - q100*y%d1val1_d1val3 - q101*y%d1val1 + q102*q29 + q103*q5 + q104*q9 - q107*y%d1val1 + q18*q49 - q2*y%d2val3 - q27*q31 + x%d1val1_d2val3)
      binary%d3val2 = q115
      binary%d2val2_d1val3 = -q116*q39 + q15*(-q1*q122 + q117 - q118 + q119 - q35*q41 + q35*q43 - q71*y%d1val2_d1val3)
      binary%d1val2_d2val3 = q0*(-q1*y%d1val2_d2val3 - q100*y%d1val2_d1val3 - q101*y%d1val2 + q103*q7 - q107*y%d1val2 + q123*q35 + q124*q9 - q31*q41 + q49*q84 - q6*y%d2val3 + x%d1val2_d2val3)
      binary%d3val1_d1val3 = -q116*q58 + q15*(-3.0_dp*q128 + 2.0_dp*x%d3val1_d1val3 - q1*(2.0_dp*y%d3val1_d1val3 - q125*q132 + q127*q132 + q133*q134 - q136*q55 - q54*y%d2val1_d1val3) - q125*q52 - q126*x%d2val1_d1val3 + q127*q52 + q129*q28 - q130 + q131 - q53*q69)
      binary%d2val1_d1val2_d1val3 = -q116*q65 - q137*q60 + q138*q14 - q139*q70 + q15*(2.0_dp*x%d2val1_d1val2_d1val3 - q1*(2.0_dp*y%d2val1_d1val2_d1val3 - q11*q143 + q127*q75 + q137*q63 - q151*y%d1val1_d1val2 + q152*q73 - q62*y%d1val1_d1val2_d1val3) - q10*y%d1val1_d1val2_d1val3 - q104*q17 + q104*q19 - q13*q47 - q140*q27 + q140*q28 + q141*x%d1val1_d1val2 + q142*q20 - q144*y%d1val1 - q145 + q146 + q147 + q148 - q149 + q150 - q6*q69 - q61*q80 + q61*q81 + q69*q7)
      binary%d2val1_d2val3 = q0*(-q158*q159 - q16*q50 + q170)
      binary%d1val1_d2val2_d1val3 = -q116*q78 + q15*(-2.0_dp*q171 - 2.0_dp*q172 + 2.0_dp*x%d1val1_d2val2_d1val3 + 36.0_dp*q189 - q117*q51 - q120*q2 + q121*q19 + q121*q21 + q121*q22 - q124*q17 - q125*q34 + q127*q34 - q140*q41 + q140*q43 + q168*q173 + q174*q5 + q175*q7 + q176*q76 - q177*q178 - q179*q180 - q181*q182 - q183*x%val - q184*q185 - q186*q77 - q187*q77 - q27*q36 + q28*q36 + q29*q74 + q30*q74 - q35*q80 + q35*q82 + q45*q75 - q71*y%d1val1_d1val2_d1val3 + q72*q9 + q73*q86)
      binary%d1val1_d1val2_d2val3 = q0*(-6.0_dp*q200*y%d1val1 - q1*y%d1val1_d1val2_d2val3 - q100*y%d1val1_d1val2_d1val3 - q101*y%d1val1_d1val2 + q102*q86 + q103*q19 + q103*q21 + q103*q22 + q104*q43 + q104*q45 - q107*y%d1val1_d1val2 - q108*x%d1val1_d2val3 + q123*q140 + q124*q28 + q124*q29 + q124*q30 + q142*q42 - q154*y%d1val2_d1val3 - q162*q95 - q17*y%d2val3 - q181*q199 - q190*y%d1val1_d1val3 + q191*q192 + q193*y%d1val2_d2val3 + q194*q9 + q195*q7 - q196*q197 - q196*x%val*y%d1val1_d1val3 - q198*q97 - q2*y%d1val2_d2val3 + q203*q96 + q3*q49*x%d1val1_d1val2 - q31*q80 - q51*x%d1val2_d2val3 - q6*y%d1val1_d2val3 + x%d1val1_d1val2_d2val3)
      binary%d3val2_d1val3 = -q114*q116 + q15*(-3.0_dp*q207 + 2.0_dp*x%d3val2_d1val3 - q1*(2.0_dp*y%d3val2_d1val3 + 36.0_dp*q214 - q111*y%d2val2_d1val3 - q112*q136 - q199*q211 + q212*q213) - q109*x%d2val2_d1val3 - q110*q122 - q204*q205 + q206*x%d2val2 + q208*q43 - q209 + q210)
      binary%d2val2_d2val3 = q0*(q1*q222 - q100*q219 + q103*q173 - q106*q38 + q123*q38 - q124*q41 + q124*q43 - q159*(-0.5_dp*q118 + 0.5_dp*q119 - q1*q219 - q190*y%d1val2 + q191*q43 - q215*y%d1val2_d1val3 + x%d2val2_d1val3) - q163*q38 + q165*q38 + q166*q219 - q200*q35 - q215*y%d1val2_d2val3 - q216*x%d1val2_d2val3 + q217*q35 - q40*q50 + x%d2val2_d2val3)
      binary%d3val1_d2val3 = q0*(-q1*(-72.0_dp*q240*y%d1val1_d1val3 + 36.0_dp*q167*q20 - q125*q239 - q126*y%d2val1_d2val3 + q127*q239 + q132*q228 - q132*q229 - q198*q237 - q223*q235 + q227*q235 + q236*y%d1val1_d2val3 + q241*q55 + y%d3val1_d2val3) - q100*q238 - q106*q57 + q123*q57 - q125*q225 + q127*q225 - q129*q162 + q129*q233 + q159*(-0.5_dp*q131 + 0.5_dp*q130 + 1.5_dp*q128 + q1*q238 + q125*q224 - q127*q224 + q157*q53 + q226*x%d2val1_d1val3 - q231*q28 - x%d3val1_d1val3) - q163*q57 + q165*q57 + q166*q238 + q169*q53 - q223*q224 + q224*q227 - q226*x%d2val1_d2val3 + q228*q52 - q229*q52 - q230*q231 + q231*q232 - q234*q27 + q234*q28 - q50*q59 + x%d3val1_d2val3)
      binary%d2val1_d1val2_d2val3 = q0*(-3.0_dp*q14*q259 - q0*q158*q91 - q1*(-q11*q247 + q11*q251 + q121*q152 - q121*q240 + q127*q175 + q140*q227 - q151*y%d1val1_d1val2_d1val3 - q155*y%d1val1_d1val2_d2val3 + q167*q265 + q228*q75 - q229*q75 - q260*y%d1val1_d2val3 + q261*y%d1val2_d2val3 + q264*y%d1val1_d2val3 - q266*y%d1val1 + y%d2val1_d1val2_d2val3) - q100*q262 + q103*q244 - q104*q80 + q104*q81 + q104*q82 + q104*q83 - q106*q64 - q108*q170 - q121*q249 + q123*q64 + q124*q245 + q129*q202*x%val + q13*q139*x%d2val3 - q13*q15*x%d1val2_d2val3 - q13*q200 + q13*q217 - q13*q89*x%d1val3 + q14*q156*y%d1val2_d1val3 + q14*q258 - q140*q162 + q140*q233 + q141*x%d1val1_d1val2_d1val3 - q143*x%d1val1_d1val3*y%d1val1 - q144*y%d1val1_d1val3 - q153*y%d1val1_d1val2_d2val3 - q157*q190 - q157*q257*q35 + q158*q263 - q159*(-0.5_dp*q145 - 0.5_dp*q149 + 0.5_dp*q146 + 0.5_dp*q147 + 0.5_dp*q148 + 0.5_dp*q150 - q1*q262 - q138*q255 - q153*y%d1val1_d1val2_d1val3 - q154*y%d1val1_d1val2 - q157*q6 + q157*q7 + q19*q94 + q191*q245 + q23*q81 - q242*y%d1val1_d1val3 + q243*q28 + q244*q31 - q249*q35 - q267*y%d1val1 + x%d2val1_d1val2_d1val3) - q16*y%d1val2_d2val3 - q160*x%d1val1_d1val2_d2val3 + q161*y%d1val1_d1val2 - q163*q64 + q165*q64 + q166*q262 + q169*q6 - q169*q7 + q18*q23*y%d1val2_d2val3 + q19*q195 - q194*q27 + q194*q28 - q198*q255*y%d1val2 + q20*q246*y%d1val2 - q230*q243 - q242*y%d1val1_d2val3 - q247*q248 + q248*q251 - q250*q61 + q252*q87 - q252*q93 + q253*q254 + q253*q4*y%d1val2_d2val3 + q256*q43 + q256*q45 + q256*q46 + q3*q48*q65 - q66*y%d2val3 + x%d2val1_d1val2_d2val3)
      binary%d1val1_d2val2_d2val3 = q0*(-12.0_dp*q179*y%d1val2_d2val3 - 12.0_dp*q220*q33 - 8.0_dp*q181*y%d2val2_d1val3 + 18.0_dp*q188*q98*y%d2val3 + 36.0_dp*q197*q277 + 4.0_dp*q18*q220 + 4.0_dp*q7*y%d1val1_d1val2_d2val3 + 72.0_dp*q96*q99*y%d1val2*y%d1val2_d1val3 - q1*y%d1val1_d2val2_d2val3 - q100*y%d1val1_d2val2_d1val3 - q101*y%d1val1_d2val2 + q103*q44*x%d1val1_d1val2 + q103*q84*y%d1val1_d1val2 + q104*q87*y%d2val2 - q117*q125 + q117*q127 - q120*q27 + q120*q28 + q121*q83 + q121*q85 + q121*q86 + q121*q88 + q123*q72 - q124*q80 + q124*q82 + q134*q277*x%val - q135*q182*q197 - q135*q185*x%d1val1_d1val3 - q140*q200 + q140*q217 + q140*q44*x%d2val3 - q159*(18.0_dp*q189 + 2.0_dp*q176*q37 + q104*q173 + q120*q5 + q124*q19 + q124*q21 + q124*q22 - q138*q75*x%val + q140*q45 - q143*q178 - q171 - q172 - q179*q199 - q184*q76 - q186*q276 - q187*q276 - q190*y%d1val1_d1val2 + q191*q82 + q194*q7 - q2*y%d2val2_d1val3 - q204*y%d1val1_d1val3 - q215*y%d1val1_d1val2_d1val3 - q242*y%d1val2_d1val3 + q243*q43 - q267*y%d1val2 + q269*y%d1val3 - q27*y%d2val2 - q275*q61 + q28*y%d2val2 + q29*q36 + q30*q36 + q35*q86 - q51*x%d2val2_d1val3 + q9*y%d1val1_d2val2 + x%d1val1_d2val2_d1val3) - q162*q36 + q164*y%d1val1_d2val2 + q166*y%d1val1_d2val2_d1val3 - q168*q275 + q168*q44*x%d1val2_d1val3 + q174*q29 + q174*q30 + q175*q45 + q175*q46 - q177*x%d1val1*y%d1val2_d1val3 - q177*x%d1val2_d1val3*y%d1val1 - q177*x%val*y%d1val1_d1val2_d1val3 + q178*q203 - q178*q271 - q178*q273 - q180*q26*y%d1val1_d1val3 - q180*q90*y%d1val1 - q183*x%d1val3 + q19*q270 + q192*q36 + q193*y%d2val2_d2val3 - q194*q41 + q194*q43 - q198*q61*x%val*y%d2val2 - q198*q76*x%d1val1 - q2*y%d2val2_d2val3 + q201*q213*q96 + q201*q77*x%d1val1 + q203*q274 - q204*y%d1val1_d2val3 + q21*q270 - q211*x%d1val1_d2val3 - q215*y%d1val1_d1val2_d2val3 - q216*x%d1val1_d1val2_d2val3 + q22*q270 + q221*q246 + q228*q34 - q229*q34 + q232*y%d2val2 + q233*q36 - q24*q77*x%d1val3*y%d1val1_d1val3 - q242*y%d1val2_d2val3 - q243*q268 - q250*q35 + q263*x%d1val1_d1val2_d1val3 + q264*x%d1val2_d2val3 + q265*x%d1val2*y%d1val1_d2val3 - q266*x%d1val2 + q269*y%d2val3 - q271*q274 - q273*q274 - q276*q32*x%d2val3 - q276*q92*y%d1val1_d2val3 - q278*q48*q96*powm1(pow5(y%val)) + q36*q4*y%d1val1_d2val3 - q49*q92*y%d1val1_d2val2 - q50*q79 - q51*x%d2val2_d2val3 + x%d1val1_d2val2_d2val3)
      binary%d3val2_d2val3 = q0*(-q0*q205*x%d2val2_d1val3 - q1*(-6.0_dp*q211*y%d1val2_d2val3 + 12.0_dp*q212*y%d2val2_d1val3 + 18.0_dp*q221*y%d1val2_d2val3 + 36.0_dp*q220*q44 + 6.0_dp*q258*y%d2val2 - q0*q199*y%d2val2_d1val3 - q109*y%d2val2_d2val3 + q112*q241 - q198*q284 + q199*q42*y%d2val2 - q213*q259 - q272*q278 + y%d3val2_d2val3) - q100*q285 - q105*q95*x%d2val2 - q106*q113 + q110*q222 + q113*q123 - q113*q163 + q113*q165 - q115*q50 + q159*(-0.5_dp*q210 + 0.5_dp*q209 + 1.5_dp*q207 + q1*q285 + q110*q219 - q212*q281 + q279*y%d1val2_d1val3 + q280*x%d2val2_d1val3 - q282*q43 - x%d3val2_d1val3) + q166*q285 - q200*q208 + q205*q42*x%d2val2 + q206*x%d2val2_d1val3 + q208*q217 + q254*q282 + q258*q281 - q268*q282 - q279*y%d1val2_d2val3 - q280*x%d2val2_d2val3 - q283*q41 + q283*q43 + x%d3val2_d2val3)
   end function div_self
   
   function div_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q0
      q0 = powm1(y)
      unary%val = q0*x%val
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function div_self_real
   
   function div_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = powm1(x%val)
      q1 = powm1(pow2(x%val))
      q2 = q1*z
      q3 = pow2(x%d1val1)
      q4 = 4.0_dp*q0
      q5 = 2.0_dp*x%d2val1 - q3*q4
      q6 = 0.5_dp*q2
      q7 = powm1(pow3(x%val))
      q8 = q7*z
      q9 = 2.0_dp*x%d1val2
      q10 = q8*q9
      q11 = q8*x%d1val3
      q12 = 2.0_dp*x%d1val1
      q13 = 2.0_dp*x%d2val2
      q14 = pow2(x%d1val2)
      q15 = q13 - q14*q4
      q16 = pow2(x%d1val3)
      q17 = q0*q16
      q18 = q0*x%d1val1
      q19 = 12.0_dp*q18
      q20 = pow3(x%d1val1)
      q21 = 12.0_dp*q1
      q22 = 2.0_dp*x%d3val1 - q19*x%d2val1 + q20*q21
      q23 = q5*x%d1val2
      q24 = 8.0_dp*q18
      q25 = 4.0_dp*q1
      q26 = q3*x%d1val2
      q27 = 2.0_dp*x%d2val1_d1val2 - q24*x%d1val1_d1val2 + q25*q26
      q28 = q25*x%d1val3
      q29 = 2.0_dp*x%d2val1_d1val3 - q24*x%d1val1_d1val3 + q28*q3
      q30 = 2.0_dp*x%d1val1_d2val2
      q31 = q0*x%d1val2
      q32 = 8.0_dp*x%d1val1_d1val2
      q33 = q31*q32
      q34 = q0*q14
      q35 = -12.0_dp*q34 + 4.0_dp*x%d2val2
      q36 = q18*q35 - q30 + q33
      q37 = q8*x%d1val2_d1val3
      q38 = 2.0_dp*x%d1val1_d1val2
      q39 = 6.0_dp*x%d1val2
      q40 = powm1(pow4(x%val))
      q41 = q40*x%d1val3*z
      q42 = q4*x%d1val3
      q43 = -3.0_dp*q17 + x%d2val3
      q44 = 2.0_dp*q18
      q45 = 12.0_dp*q31
      q46 = pow3(x%d1val2)
      q47 = 2.0_dp*x%d3val2 + q21*q46 - q45*x%d2val2
      q48 = q31*x%d1val2_d1val3
      q49 = q0*q43
      q50 = q0*x%d1val1_d1val3
      q51 = q50*x%d2val1
      q52 = q21*x%d1val3
      q53 = x%d1val1*x%d2val1
      q54 = q1*q3
      q55 = q54*x%d1val1_d1val3
      q56 = q20*q7
      q57 = 24.0_dp*x%d1val3
      q58 = 3.0_dp*q23
      q59 = q1*x%d1val1
      q60 = q59*x%d1val3
      q61 = 8.0_dp*x%d1val1_d1val3
      q62 = q59*q61
      q63 = q3*x%d1val2_d1val3
      q64 = 8.0_dp*x%d1val3
      q65 = q26*q7
      q66 = q4*x%d1val1
      q67 = 2.0_dp*q54
      q68 = -q66*x%d1val1_d1val3 + q67*x%d1val3 + x%d2val1_d1val3
      q69 = pow2(x%d1val1_d1val3)
      q70 = q16*q7
      q71 = 4.0_dp*q70
      q72 = q3*q71 + q4*q69 - q60*q61 + q66*x%d1val1_d2val3 - q67*x%d2val3 - x%d2val1_d2val3
      q73 = q0*x%d1val2_d1val3
      q74 = 8.0_dp*x%d1val1_d1val2_d1val3
      q75 = q1*x%d1val2
      q76 = q32*x%d1val3
      q77 = q14*x%d1val3
      q78 = q0*x%d2val3
      q79 = q4*x%d1val2_d1val3
      q80 = q0*x%d1val1_d2val3
      q81 = q59*x%d2val3
      q82 = q52*x%d1val1
      q83 = q52*x%d1val2
      q84 = q1*q16
      q85 = q70*x%d1val1
      q86 = 12.0_dp*q73
      q87 = q1*q14
      q88 = q87*x%d1val2_d1val3
      q89 = q46*q7
      q90 = -x%d2val2_d2val3
      q91 = q4*x%d1val2_d2val3
      q92 = pow2(x%d1val2_d1val3)
      q93 = q87*x%d2val3
      q94 = q1*q77
      q95 = 6.0_dp*q18
      q96 = 6.0_dp*x%d2val1
      q97 = 12.0_dp*q56
      q98 = 12.0_dp*q70
      q99 = x%d1val1_d1val3*x%d1val3
      q100 = 72.0_dp*q7
      q101 = q16*q40
      q102 = 36.0_dp*q101
      q103 = q4*x%d1val1_d1val2
      q104 = x%d1val1*x%d1val1_d1val2
      q105 = q25*x%d1val2
      q106 = q105*x%d1val1
      q107 = q1*q76
      q108 = 4.0_dp*q65
      q109 = q0*x%d1val2_d2val3
      q110 = 6.0_dp*x%d1val3
      q111 = q4*x%d1val2
      q112 = -6.0_dp*q34 + q13
      q113 = 6.0_dp*q0
      q114 = q113*x%d1val2
      q115 = 3.0_dp*q94 - q114*x%d1val2_d1val3 + x%d2val2_d1val3
      q116 = 6.0_dp*q109
      q117 = q75*x%d2val2
      q118 = x%d1val2_d1val3*x%d2val2
      q119 = 12.0_dp*q89
      unary%val = q0*z
      unary%d1val1 = -q2*x%d1val1
      unary%d1val2 = -q2*x%d1val2
      unary%d1val3 = -q2*x%d1val3
      unary%d2val1 = -q5*q6
      unary%d1val1_d1val2 = q10*x%d1val1 - q2*x%d1val1_d1val2
      unary%d1val1_d1val3 = q11*q12 - q2*x%d1val1_d1val3
      unary%d2val2 = -q15*q6
      unary%d1val2_d1val3 = q11*q9 - q2*x%d1val2_d1val3
      unary%d2val3 = -q2*(-2.0_dp*q17 + x%d2val3)
      unary%d3val1 = -q22*q6
      unary%d2val1_d1val2 = q23*q8 - q27*q6
      unary%d2val1_d1val3 = q11*q5 - q29*q6
      unary%d1val1_d2val2 = q36*q6
      unary%d1val1_d1val2_d1val3 = q10*x%d1val1_d1val3 + q11*q38 + q12*q37 - q2*x%d1val1_d1val2_d1val3 - q39*q41*x%d1val1
      unary%d1val1_d2val3 = q2*(q42*x%d1val1_d1val3 + q43*q44 - x%d1val1_d2val3)
      unary%d3val2 = -q47*q6
      unary%d2val2_d1val3 = q11*q15 - q6*(-8.0_dp*q48 + 2.0_dp*x%d2val2_d1val3 + q14*q28)
      unary%d1val2_d2val3 = q2*(q42*x%d1val2_d1val3 + q49*q9 - x%d1val2_d2val3)
      unary%d3val1_d1val3 = q11*q22 - q6*(-12.0_dp*q51 + 2.0_dp*x%d3val1_d1val3 + 36.0_dp*q55 - q19*x%d2val1_d1val3 + q52*q53 - q56*q57)
      unary%d2val1_d1val2_d1val3 = q11*q27 + q29*q8*x%d1val2 + q37*q5 - q41*q58 - q6*(2.0_dp*x%d2val1_d1val2_d1val3 - q24*x%d1val1_d1val2_d1val3 + q25*q63 - q32*q50 + q32*q60 + q62*x%d1val2 - q64*q65)
      unary%d2val1_d2val3 = q2*(q42*q68 + q49*q5 + q72)
      unary%d1val1_d2val2_d1val3 = -q11*q36 + q6*(-2.0_dp*x%d1val1_d2val2_d1val3 + q18*(-24.0_dp*q48 + 4.0_dp*x%d2val2_d1val3 + q21*q77) + q31*q74 + q32*q73 + q35*q50 - q35*q60 - q75*q76)
      unary%d1val1_d1val2_d2val3 = q2*(-6.0_dp*q84*x%d1val1_d1val2 + 24.0_dp*q85*x%d1val2 + q38*q78 - q39*q81 + q42*x%d1val1_d1val2_d1val3 + q44*x%d1val2_d2val3 + q79*x%d1val1_d1val3 + q80*q9 - q82*x%d1val2_d1val3 - q83*x%d1val1_d1val3 - x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q11*q47 - q6*(2.0_dp*x%d3val2_d1val3 + 36.0_dp*q88 - q45*x%d2val2_d1val3 - q57*q89 + q83*x%d2val2 - q86*x%d2val2)
      unary%d2val2_d2val3 = q2*(-2.0_dp*q93 + q14*q71 + q15*q49 + q4*q92 + q42*(2.0_dp*q94 - q79*x%d1val2 + x%d2val2_d1val3) - q64*q75*x%d1val2_d1val3 + q90 + q91*x%d1val2)
      unary%d3val1_d2val3 = q2*(-18.0_dp*q54*x%d1val1_d2val3 - 36.0_dp*q59*q69 + 12.0_dp*q50*x%d2val1_d1val3 + q100*q3*q99 - q102*q20 + q22*q49 + q42*(-6.0_dp*q51 + 18.0_dp*q55 + q60*q96 - q95*x%d2val1_d1val3 - q97*x%d1val3 + x%d3val1_d1val3) - q52*x%d1val1_d1val3*x%d2val1 + q53*q98 + q80*q96 - q81*q96 - q82*x%d2val1_d1val3 + q95*x%d2val1_d2val3 + q97*x%d2val3 - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(-12.0_dp*q101*q26 - 2.0_dp*q31*q72 - 3.0_dp*q27*q84 + 16.0_dp*q7*q99*x%d1val1*x%d1val2 - q1*q110*q5*x%d1val2_d1val3 - q1*q58*x%d2val3 + q103*x%d1val1_d2val3 - q104*q25*x%d2val3 - q105*q69 - q106*x%d1val1_d2val3 - q107*x%d1val1_d1val3 + q108*x%d2val3 + q109*q5 + q23*q98 + q27*q78 + q32*q85 + q42*(-q103*x%d1val1_d1val3 + q104*q28 + q106*x%d1val1_d1val3 - q108*x%d1val3 - q66*x%d1val1_d1val2_d1val3 + q67*x%d1val2_d1val3 + x%d2val1_d1val2_d1val3) + q50*q74 - q60*q74 - q62*x%d1val2_d1val3 + q63*q64*q7 + q66*x%d1val1_d1val2_d2val3 - q67*x%d1val2_d2val3 + q68*q79 - q68*q83 - x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(-2.0_dp*q1*q112*q99 - q105*x%d1val1_d1val2*x%d2val3 - q107*x%d1val2_d1val3 + q111*x%d1val1_d1val2_d2val3 + q112*q12*q70 + q112*q80 - q112*q81 - q115*q28*x%d1val1 + q115*q4*x%d1val1_d1val3 + q32*q70*x%d1val2 + q42*(-q111*x%d1val1_d1val2_d1val3 - q112*q50 + q112*q60 - q115*q44 + q28*x%d1val1_d1val2*x%d1val2 - q79*x%d1val1_d1val2 + x%d1val1_d2val2_d1val3) - q44*(-3.0_dp*q93 + 6.0_dp*q14*q70 + q113*q92 + q116*x%d1val2 - q83*x%d1val2_d1val3 + q90) + q49*(-q112*q44 + q30 - q33) + q73*q74 - q74*q75*x%d1val3 + q91*x%d1val1_d1val2 - x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(-18.0_dp*q87*x%d1val2_d2val3 - 36.0_dp*q75*q92 - 6.0_dp*q117*x%d2val3 + q100*q77*x%d1val2_d1val3 - q102*q46 + q114*x%d2val2_d2val3 + q116*x%d2val2 - q118*q52 + q119*x%d2val3 + q42*(18.0_dp*q88 + q110*q117 - q113*q118 - q114*x%d2val2_d1val3 - q119*x%d1val3 + x%d3val2_d1val3) + q47*q49 - q83*x%d2val2_d1val3 + q86*x%d2val2_d1val3 + q98*x%d1val2*x%d2val2 - x%d3val2_d2val3)
   end function div_real_self
   
   function div_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = y
      q0 = powm1(y_dp)
      unary%val = q0*x%val
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function div_self_int
   
   function div_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = powm1(x%val)
      q1 = powm1(pow2(x%val))
      q2 = q1*y_dp
      q3 = pow2(x%d1val1)
      q4 = 4.0_dp*q0
      q5 = 2.0_dp*x%d2val1 - q3*q4
      q6 = 0.5_dp*q2
      q7 = powm1(pow3(x%val))
      q8 = q7*y_dp
      q9 = 2.0_dp*x%d1val2
      q10 = q8*q9
      q11 = q8*x%d1val3
      q12 = 2.0_dp*x%d1val1
      q13 = 2.0_dp*x%d2val2
      q14 = pow2(x%d1val2)
      q15 = q13 - q14*q4
      q16 = pow2(x%d1val3)
      q17 = q0*q16
      q18 = q0*x%d1val1
      q19 = 12.0_dp*q18
      q20 = pow3(x%d1val1)
      q21 = 12.0_dp*q1
      q22 = 2.0_dp*x%d3val1 - q19*x%d2val1 + q20*q21
      q23 = q5*x%d1val2
      q24 = 8.0_dp*q18
      q25 = 4.0_dp*q1
      q26 = q3*x%d1val2
      q27 = 2.0_dp*x%d2val1_d1val2 - q24*x%d1val1_d1val2 + q25*q26
      q28 = q25*x%d1val3
      q29 = 2.0_dp*x%d2val1_d1val3 - q24*x%d1val1_d1val3 + q28*q3
      q30 = 2.0_dp*x%d1val1_d2val2
      q31 = q0*x%d1val2
      q32 = 8.0_dp*x%d1val1_d1val2
      q33 = q31*q32
      q34 = q0*q14
      q35 = -12.0_dp*q34 + 4.0_dp*x%d2val2
      q36 = q18*q35 - q30 + q33
      q37 = q8*x%d1val2_d1val3
      q38 = 2.0_dp*x%d1val1_d1val2
      q39 = 6.0_dp*x%d1val2
      q40 = powm1(pow4(x%val))
      q41 = q40*x%d1val3*y_dp
      q42 = q4*x%d1val3
      q43 = -3.0_dp*q17 + x%d2val3
      q44 = 2.0_dp*q18
      q45 = 12.0_dp*q31
      q46 = pow3(x%d1val2)
      q47 = 2.0_dp*x%d3val2 + q21*q46 - q45*x%d2val2
      q48 = q31*x%d1val2_d1val3
      q49 = q0*q43
      q50 = q0*x%d1val1_d1val3
      q51 = q50*x%d2val1
      q52 = q21*x%d1val3
      q53 = x%d1val1*x%d2val1
      q54 = q1*q3
      q55 = q54*x%d1val1_d1val3
      q56 = q20*q7
      q57 = 24.0_dp*x%d1val3
      q58 = 3.0_dp*q23
      q59 = q1*x%d1val1
      q60 = q59*x%d1val3
      q61 = 8.0_dp*x%d1val1_d1val3
      q62 = q59*q61
      q63 = q3*x%d1val2_d1val3
      q64 = 8.0_dp*x%d1val3
      q65 = q26*q7
      q66 = q4*x%d1val1
      q67 = 2.0_dp*q54
      q68 = -q66*x%d1val1_d1val3 + q67*x%d1val3 + x%d2val1_d1val3
      q69 = pow2(x%d1val1_d1val3)
      q70 = q16*q7
      q71 = 4.0_dp*q70
      q72 = q3*q71 + q4*q69 - q60*q61 + q66*x%d1val1_d2val3 - q67*x%d2val3 - x%d2val1_d2val3
      q73 = q0*x%d1val2_d1val3
      q74 = 8.0_dp*x%d1val1_d1val2_d1val3
      q75 = q1*x%d1val2
      q76 = q32*x%d1val3
      q77 = q14*x%d1val3
      q78 = q0*x%d2val3
      q79 = q4*x%d1val2_d1val3
      q80 = q0*x%d1val1_d2val3
      q81 = q59*x%d2val3
      q82 = q52*x%d1val1
      q83 = q52*x%d1val2
      q84 = q1*q16
      q85 = q70*x%d1val1
      q86 = 12.0_dp*q73
      q87 = q1*q14
      q88 = q87*x%d1val2_d1val3
      q89 = q46*q7
      q90 = -x%d2val2_d2val3
      q91 = q4*x%d1val2_d2val3
      q92 = pow2(x%d1val2_d1val3)
      q93 = q87*x%d2val3
      q94 = q1*q77
      q95 = 6.0_dp*q18
      q96 = 6.0_dp*x%d2val1
      q97 = 12.0_dp*q56
      q98 = 12.0_dp*q70
      q99 = x%d1val1_d1val3*x%d1val3
      q100 = 72.0_dp*q7
      q101 = q16*q40
      q102 = 36.0_dp*q101
      q103 = q4*x%d1val1_d1val2
      q104 = x%d1val1*x%d1val1_d1val2
      q105 = q25*x%d1val2
      q106 = q105*x%d1val1
      q107 = q1*q76
      q108 = 4.0_dp*q65
      q109 = q0*x%d1val2_d2val3
      q110 = 6.0_dp*x%d1val3
      q111 = q4*x%d1val2
      q112 = -6.0_dp*q34 + q13
      q113 = 6.0_dp*q0
      q114 = q113*x%d1val2
      q115 = 3.0_dp*q94 - q114*x%d1val2_d1val3 + x%d2val2_d1val3
      q116 = 6.0_dp*q109
      q117 = q75*x%d2val2
      q118 = x%d1val2_d1val3*x%d2val2
      q119 = 12.0_dp*q89
      unary%val = q0*y_dp
      unary%d1val1 = -q2*x%d1val1
      unary%d1val2 = -q2*x%d1val2
      unary%d1val3 = -q2*x%d1val3
      unary%d2val1 = -q5*q6
      unary%d1val1_d1val2 = q10*x%d1val1 - q2*x%d1val1_d1val2
      unary%d1val1_d1val3 = q11*q12 - q2*x%d1val1_d1val3
      unary%d2val2 = -q15*q6
      unary%d1val2_d1val3 = q11*q9 - q2*x%d1val2_d1val3
      unary%d2val3 = -q2*(-2.0_dp*q17 + x%d2val3)
      unary%d3val1 = -q22*q6
      unary%d2val1_d1val2 = q23*q8 - q27*q6
      unary%d2val1_d1val3 = q11*q5 - q29*q6
      unary%d1val1_d2val2 = q36*q6
      unary%d1val1_d1val2_d1val3 = q10*x%d1val1_d1val3 + q11*q38 + q12*q37 - q2*x%d1val1_d1val2_d1val3 - q39*q41*x%d1val1
      unary%d1val1_d2val3 = q2*(q42*x%d1val1_d1val3 + q43*q44 - x%d1val1_d2val3)
      unary%d3val2 = -q47*q6
      unary%d2val2_d1val3 = q11*q15 - q6*(-8.0_dp*q48 + 2.0_dp*x%d2val2_d1val3 + q14*q28)
      unary%d1val2_d2val3 = q2*(q42*x%d1val2_d1val3 + q49*q9 - x%d1val2_d2val3)
      unary%d3val1_d1val3 = q11*q22 - q6*(-12.0_dp*q51 + 2.0_dp*x%d3val1_d1val3 + 36.0_dp*q55 - q19*x%d2val1_d1val3 + q52*q53 - q56*q57)
      unary%d2val1_d1val2_d1val3 = q11*q27 + q29*q8*x%d1val2 + q37*q5 - q41*q58 - q6*(2.0_dp*x%d2val1_d1val2_d1val3 - q24*x%d1val1_d1val2_d1val3 + q25*q63 - q32*q50 + q32*q60 + q62*x%d1val2 - q64*q65)
      unary%d2val1_d2val3 = q2*(q42*q68 + q49*q5 + q72)
      unary%d1val1_d2val2_d1val3 = -q11*q36 + q6*(-2.0_dp*x%d1val1_d2val2_d1val3 + q18*(-24.0_dp*q48 + 4.0_dp*x%d2val2_d1val3 + q21*q77) + q31*q74 + q32*q73 + q35*q50 - q35*q60 - q75*q76)
      unary%d1val1_d1val2_d2val3 = q2*(-6.0_dp*q84*x%d1val1_d1val2 + 24.0_dp*q85*x%d1val2 + q38*q78 - q39*q81 + q42*x%d1val1_d1val2_d1val3 + q44*x%d1val2_d2val3 + q79*x%d1val1_d1val3 + q80*q9 - q82*x%d1val2_d1val3 - q83*x%d1val1_d1val3 - x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q11*q47 - q6*(2.0_dp*x%d3val2_d1val3 + 36.0_dp*q88 - q45*x%d2val2_d1val3 - q57*q89 + q83*x%d2val2 - q86*x%d2val2)
      unary%d2val2_d2val3 = q2*(-2.0_dp*q93 + q14*q71 + q15*q49 + q4*q92 + q42*(2.0_dp*q94 - q79*x%d1val2 + x%d2val2_d1val3) - q64*q75*x%d1val2_d1val3 + q90 + q91*x%d1val2)
      unary%d3val1_d2val3 = q2*(-18.0_dp*q54*x%d1val1_d2val3 - 36.0_dp*q59*q69 + 12.0_dp*q50*x%d2val1_d1val3 + q100*q3*q99 - q102*q20 + q22*q49 + q42*(-6.0_dp*q51 + 18.0_dp*q55 + q60*q96 - q95*x%d2val1_d1val3 - q97*x%d1val3 + x%d3val1_d1val3) - q52*x%d1val1_d1val3*x%d2val1 + q53*q98 + q80*q96 - q81*q96 - q82*x%d2val1_d1val3 + q95*x%d2val1_d2val3 + q97*x%d2val3 - x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(-12.0_dp*q101*q26 - 2.0_dp*q31*q72 - 3.0_dp*q27*q84 + 16.0_dp*q7*q99*x%d1val1*x%d1val2 - q1*q110*q5*x%d1val2_d1val3 - q1*q58*x%d2val3 + q103*x%d1val1_d2val3 - q104*q25*x%d2val3 - q105*q69 - q106*x%d1val1_d2val3 - q107*x%d1val1_d1val3 + q108*x%d2val3 + q109*q5 + q23*q98 + q27*q78 + q32*q85 + q42*(-q103*x%d1val1_d1val3 + q104*q28 + q106*x%d1val1_d1val3 - q108*x%d1val3 - q66*x%d1val1_d1val2_d1val3 + q67*x%d1val2_d1val3 + x%d2val1_d1val2_d1val3) + q50*q74 - q60*q74 - q62*x%d1val2_d1val3 + q63*q64*q7 + q66*x%d1val1_d1val2_d2val3 - q67*x%d1val2_d2val3 + q68*q79 - q68*q83 - x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(-2.0_dp*q1*q112*q99 - q105*x%d1val1_d1val2*x%d2val3 - q107*x%d1val2_d1val3 + q111*x%d1val1_d1val2_d2val3 + q112*q12*q70 + q112*q80 - q112*q81 - q115*q28*x%d1val1 + q115*q4*x%d1val1_d1val3 + q32*q70*x%d1val2 + q42*(-q111*x%d1val1_d1val2_d1val3 - q112*q50 + q112*q60 - q115*q44 + q28*x%d1val1_d1val2*x%d1val2 - q79*x%d1val1_d1val2 + x%d1val1_d2val2_d1val3) - q44*(-3.0_dp*q93 + 6.0_dp*q14*q70 + q113*q92 + q116*x%d1val2 - q83*x%d1val2_d1val3 + q90) + q49*(-q112*q44 + q30 - q33) + q73*q74 - q74*q75*x%d1val3 + q91*x%d1val1_d1val2 - x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(-18.0_dp*q87*x%d1val2_d2val3 - 36.0_dp*q75*q92 - 6.0_dp*q117*x%d2val3 + q100*q77*x%d1val2_d1val3 - q102*q46 + q114*x%d2val2_d2val3 + q116*x%d2val2 - q118*q52 + q119*x%d2val3 + q42*(18.0_dp*q88 + q110*q117 - q113*q118 - q114*x%d2val2_d1val3 - q119*x%d1val3 + x%d3val2_d1val3) + q47*q49 - q83*x%d2val2_d1val3 + q86*x%d2val2_d1val3 + q98*x%d1val2*x%d2val2 - x%d3val2_d2val3)
   end function div_int_self
   
   function pow_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: binary
      real(dp) :: q493
      real(dp) :: q492
      real(dp) :: q491
      real(dp) :: q490
      real(dp) :: q489
      real(dp) :: q488
      real(dp) :: q487
      real(dp) :: q486
      real(dp) :: q485
      real(dp) :: q484
      real(dp) :: q483
      real(dp) :: q482
      real(dp) :: q481
      real(dp) :: q480
      real(dp) :: q479
      real(dp) :: q478
      real(dp) :: q477
      real(dp) :: q476
      real(dp) :: q475
      real(dp) :: q474
      real(dp) :: q473
      real(dp) :: q472
      real(dp) :: q471
      real(dp) :: q470
      real(dp) :: q469
      real(dp) :: q468
      real(dp) :: q467
      real(dp) :: q466
      real(dp) :: q465
      real(dp) :: q464
      real(dp) :: q463
      real(dp) :: q462
      real(dp) :: q461
      real(dp) :: q460
      real(dp) :: q459
      real(dp) :: q458
      real(dp) :: q457
      real(dp) :: q456
      real(dp) :: q455
      real(dp) :: q454
      real(dp) :: q453
      real(dp) :: q452
      real(dp) :: q451
      real(dp) :: q450
      real(dp) :: q449
      real(dp) :: q448
      real(dp) :: q447
      real(dp) :: q446
      real(dp) :: q445
      real(dp) :: q444
      real(dp) :: q443
      real(dp) :: q442
      real(dp) :: q441
      real(dp) :: q440
      real(dp) :: q439
      real(dp) :: q438
      real(dp) :: q437
      real(dp) :: q436
      real(dp) :: q435
      real(dp) :: q434
      real(dp) :: q433
      real(dp) :: q432
      real(dp) :: q431
      real(dp) :: q430
      real(dp) :: q429
      real(dp) :: q428
      real(dp) :: q427
      real(dp) :: q426
      real(dp) :: q425
      real(dp) :: q424
      real(dp) :: q423
      real(dp) :: q422
      real(dp) :: q421
      real(dp) :: q420
      real(dp) :: q419
      real(dp) :: q418
      real(dp) :: q417
      real(dp) :: q416
      real(dp) :: q415
      real(dp) :: q414
      real(dp) :: q413
      real(dp) :: q412
      real(dp) :: q411
      real(dp) :: q410
      real(dp) :: q409
      real(dp) :: q408
      real(dp) :: q407
      real(dp) :: q406
      real(dp) :: q405
      real(dp) :: q404
      real(dp) :: q403
      real(dp) :: q402
      real(dp) :: q401
      real(dp) :: q400
      real(dp) :: q399
      real(dp) :: q398
      real(dp) :: q397
      real(dp) :: q396
      real(dp) :: q395
      real(dp) :: q394
      real(dp) :: q393
      real(dp) :: q392
      real(dp) :: q391
      real(dp) :: q390
      real(dp) :: q389
      real(dp) :: q388
      real(dp) :: q387
      real(dp) :: q386
      real(dp) :: q385
      real(dp) :: q384
      real(dp) :: q383
      real(dp) :: q382
      real(dp) :: q381
      real(dp) :: q380
      real(dp) :: q379
      real(dp) :: q378
      real(dp) :: q377
      real(dp) :: q376
      real(dp) :: q375
      real(dp) :: q374
      real(dp) :: q373
      real(dp) :: q372
      real(dp) :: q371
      real(dp) :: q370
      real(dp) :: q369
      real(dp) :: q368
      real(dp) :: q367
      real(dp) :: q366
      real(dp) :: q365
      real(dp) :: q364
      real(dp) :: q363
      real(dp) :: q362
      real(dp) :: q361
      real(dp) :: q360
      real(dp) :: q359
      real(dp) :: q358
      real(dp) :: q357
      real(dp) :: q356
      real(dp) :: q355
      real(dp) :: q354
      real(dp) :: q353
      real(dp) :: q352
      real(dp) :: q351
      real(dp) :: q350
      real(dp) :: q349
      real(dp) :: q348
      real(dp) :: q347
      real(dp) :: q346
      real(dp) :: q345
      real(dp) :: q344
      real(dp) :: q343
      real(dp) :: q342
      real(dp) :: q341
      real(dp) :: q340
      real(dp) :: q339
      real(dp) :: q338
      real(dp) :: q337
      real(dp) :: q336
      real(dp) :: q335
      real(dp) :: q334
      real(dp) :: q333
      real(dp) :: q332
      real(dp) :: q331
      real(dp) :: q330
      real(dp) :: q329
      real(dp) :: q328
      real(dp) :: q327
      real(dp) :: q326
      real(dp) :: q325
      real(dp) :: q324
      real(dp) :: q323
      real(dp) :: q322
      real(dp) :: q321
      real(dp) :: q320
      real(dp) :: q319
      real(dp) :: q318
      real(dp) :: q317
      real(dp) :: q316
      real(dp) :: q315
      real(dp) :: q314
      real(dp) :: q313
      real(dp) :: q312
      real(dp) :: q311
      real(dp) :: q310
      real(dp) :: q309
      real(dp) :: q308
      real(dp) :: q307
      real(dp) :: q306
      real(dp) :: q305
      real(dp) :: q304
      real(dp) :: q303
      real(dp) :: q302
      real(dp) :: q301
      real(dp) :: q300
      real(dp) :: q299
      real(dp) :: q298
      real(dp) :: q297
      real(dp) :: q296
      real(dp) :: q295
      real(dp) :: q294
      real(dp) :: q293
      real(dp) :: q292
      real(dp) :: q291
      real(dp) :: q290
      real(dp) :: q289
      real(dp) :: q288
      real(dp) :: q287
      real(dp) :: q286
      real(dp) :: q285
      real(dp) :: q284
      real(dp) :: q283
      real(dp) :: q282
      real(dp) :: q281
      real(dp) :: q280
      real(dp) :: q279
      real(dp) :: q278
      real(dp) :: q277
      real(dp) :: q276
      real(dp) :: q275
      real(dp) :: q274
      real(dp) :: q273
      real(dp) :: q272
      real(dp) :: q271
      real(dp) :: q270
      real(dp) :: q269
      real(dp) :: q268
      real(dp) :: q267
      real(dp) :: q266
      real(dp) :: q265
      real(dp) :: q264
      real(dp) :: q263
      real(dp) :: q262
      real(dp) :: q261
      real(dp) :: q260
      real(dp) :: q259
      real(dp) :: q258
      real(dp) :: q257
      real(dp) :: q256
      real(dp) :: q255
      real(dp) :: q254
      real(dp) :: q253
      real(dp) :: q252
      real(dp) :: q251
      real(dp) :: q250
      real(dp) :: q249
      real(dp) :: q248
      real(dp) :: q247
      real(dp) :: q246
      real(dp) :: q245
      real(dp) :: q244
      real(dp) :: q243
      real(dp) :: q242
      real(dp) :: q241
      real(dp) :: q240
      real(dp) :: q239
      real(dp) :: q238
      real(dp) :: q237
      real(dp) :: q236
      real(dp) :: q235
      real(dp) :: q234
      real(dp) :: q233
      real(dp) :: q232
      real(dp) :: q231
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q11 = q1*y%d1val3 + q4*x%d1val3
      q12 = q0*q11
      q13 = 16.0_dp*q1
      q14 = q3*x%d1val1
      q15 = q14*y%d1val1
      q16 = 16.0_dp*q4
      q17 = powm1(pow2(x%val))
      q18 = pow2(x%d1val1)
      q19 = q17*q18
      q20 = 16.0_dp*y%val
      q21 = 4.0_dp*q2 + 4.0_dp*q5
      q22 = pow2(q21)
      q23 = 32.0_dp*q15 + q13*y%d2val1 + q16*x%d2val1 - q19*q20 + q22
      q24 = 0.0625_dp*q0
      q25 = q6*q9
      q26 = q1*y%d1val1_d1val2
      q27 = q14*y%d1val2
      q28 = q4*x%d1val1_d1val2
      q29 = q3*y%d1val1
      q30 = q29*x%d1val2
      q31 = q17*x%d1val1
      q32 = x%d1val2*y%val
      q33 = q31*q32
      q34 = q26 + q27 + q28 + q30 - q33
      q35 = q1*y%d1val1_d1val3
      q36 = q14*y%d1val3
      q37 = q4*x%d1val1_d1val3
      q38 = q29*x%d1val3
      q39 = q31*y%val
      q40 = q39*x%d1val3
      q41 = q35 + q36 + q37 + q38 - q40
      q42 = q3*y%d1val2
      q43 = q42*x%d1val2
      q44 = pow2(x%d1val2)
      q45 = q17*q44
      q46 = 4.0_dp*q7 + 4.0_dp*q8
      q47 = pow2(q46)
      q48 = 32.0_dp*q43 + q13*y%d2val2 + q16*x%d2val2 - q20*q45 + q47
      q49 = q10*q11
      q50 = q1*y%d1val2_d1val3
      q51 = q3*y%d1val3
      q52 = q51*x%d1val2
      q53 = q4*x%d1val2_d1val3
      q54 = q42*x%d1val3
      q55 = q17*x%d1val3
      q56 = q32*q55
      q57 = q50 + q52 + q53 + q54 - q56
      q58 = q0*q57
      q59 = pow2(q11)
      q60 = 2.0_dp*q51
      q61 = pow2(x%d1val3)
      q62 = q17*q61
      q63 = q1*y%d2val3 + q4*x%d2val3 + q60*x%d1val3 - q62*y%val
      q64 = q59 + q63
      q65 = 64.0_dp*q1
      q66 = 192.0_dp*q14
      q67 = 192.0_dp*q29
      q68 = 64.0_dp*q4
      q69 = 192.0_dp*x%d2val1
      q70 = 192.0_dp*q19
      q71 = powm1(pow3(x%val))
      q72 = 128.0_dp*q71
      q73 = pow3(x%d1val1)
      q74 = q73*y%val
      q75 = 4.0_dp*q1
      q76 = 4.0_dp*q4
      q77 = 4.0_dp*y%val
      q78 = 8.0_dp*q15 - q19*q77 + q75*y%d2val1 + q76*x%d2val1
      q79 = 12.0_dp*q21*q78 - q39*q69 + q65*y%d3val1 + q66*y%d2val1 + q67*x%d2val1 + q68*x%d3val1 - q70*y%d1val1 + q72*q74 + pow3(q21)
      q80 = 0.015625_dp*q0
      q81 = 0.0625_dp*q23
      q82 = 32.0_dp*q14
      q83 = 32.0_dp*q29
      q84 = q3*x%d1val2
      q85 = q84*y%d2val1
      q86 = q42*x%d2val1
      q87 = 32.0_dp*q39
      q88 = 32.0_dp*y%d1val1
      q89 = q31*x%d1val2
      q90 = q17*x%d2val1
      q91 = q20*x%d1val2
      q92 = q19*y%d1val2
      q93 = 32.0_dp*q71
      q94 = q18*q32
      q95 = -8.0_dp*q33 + 8.0_dp*q26 + 8.0_dp*q27 + 8.0_dp*q28 + 8.0_dp*q30
      q96 = -16.0_dp*q92 + 16.0_dp*q85 + 16.0_dp*q86 + q13*y%d2val1_d1val2 + q16*x%d2val1_d1val2 + q21*q95 + q82*y%d1val1_d1val2 + q83*x%d1val1_d1val2 - q87*x%d1val1_d1val2 - q88*q89 - q90*q91 + q93*q94
      q97 = q1*y%d2val1_d1val3
      q98 = q3*x%d1val3
      q99 = q98*y%d2val1
      q100 = q51*x%d2val1
      q101 = q4*x%d2val1_d1val3
      q102 = q31*x%d1val3
      q103 = q55*x%d2val1
      q104 = q19*y%d1val3
      q105 = q18*y%val
      q106 = q93*x%d1val3
      q107 = -16.0_dp*q104 + 16.0_dp*q100 + 16.0_dp*q101 + 16.0_dp*q97 + 16.0_dp*q99 - q102*q88 - q103*q20 + q105*q106 + q21*(-8.0_dp*q40 + 8.0_dp*q35 + 8.0_dp*q36 + 8.0_dp*q37 + 8.0_dp*q38) + q82*y%d1val1_d1val3 + q83*x%d1val1_d1val3 - q87*x%d1val1_d1val3
      q108 = q1*y%d1val1_d2val2
      q109 = q14*y%d2val2
      q110 = q4*x%d1val1_d2val2
      q111 = q29*x%d2val2
      q112 = 2.0_dp*x%d1val1_d1val2
      q113 = 2.0_dp*q84
      q114 = 2.0_dp*q89
      q115 = q39*x%d2val2
      q116 = q17*q32
      q117 = q45*y%d1val1
      q118 = 2.0_dp*y%val
      q119 = q118*q71
      q120 = q44*x%d1val1
      q121 = -4.0_dp*q33 + 4.0_dp*q26 + 4.0_dp*q27 + 4.0_dp*q28 + 4.0_dp*q30
      q122 = 0.015625_dp*q21*q48 + 0.125_dp*q121*q46 + q108 + q109 + q110 + q111 - q112*q116 + q112*q42 + q113*y%d1val1_d1val2 - q114*y%d1val2 - q115 - q117 + q119*q120
      q123 = q1*y%d1val1_d1val2_d1val3
      q124 = q14*y%d1val2_d1val3
      q125 = q51*x%d1val1_d1val2
      q126 = q4*x%d1val1_d1val2_d1val3
      q127 = q42*x%d1val1_d1val3
      q128 = q84*y%d1val1_d1val3
      q129 = q29*x%d1val2_d1val3
      q130 = q98*y%d1val1_d1val2
      q131 = q31*y%d1val3
      q132 = q131*x%d1val2
      q133 = q39*x%d1val2_d1val3
      q134 = q102*y%d1val2
      q135 = q55*y%val
      q136 = q135*x%d1val1_d1val2
      q137 = q116*x%d1val1_d1val3
      q138 = q55*y%d1val1
      q139 = q138*x%d1val2
      q140 = 2.0_dp*x%d1val3
      q141 = q71*x%d1val1
      q142 = q141*q32
      q143 = q123 + q124 + q125 + q126 + q127 + q128 + q129 + q130 - q132 - q133 - q134 - q136 - q137 - q139 + q140*q142
      q144 = 2.0_dp*q2 + 2.0_dp*q5
      q145 = 0.5_dp*q64
      q146 = 2.0_dp*q11
      q147 = 2.0_dp*q98
      q148 = 2.0_dp*q135
      q149 = q119*q61
      q150 = q1*y%d1val1_d2val3 - q131*q140 + q14*y%d2val3 + q147*y%d1val1_d1val3 - q148*x%d1val1_d1val3 + q149*x%d1val1 + q29*x%d2val3 - q39*x%d2val3 + q4*x%d1val1_d2val3 + q60*x%d1val1_d1val3 - q62*y%d1val1
      q151 = 192.0_dp*q84
      q152 = 192.0_dp*q42
      q153 = 192.0_dp*x%d2val2
      q154 = 192.0_dp*q45
      q155 = pow3(x%d1val2)
      q156 = q155*q72
      q157 = 8.0_dp*q43 - q45*q77 + q75*y%d2val2 + q76*x%d2val2
      q158 = 12.0_dp*q157*q46 - q116*q153 + q151*y%d2val2 + q152*x%d2val2 - q154*y%d1val2 + q156*y%val + q65*y%d3val2 + q68*x%d3val2 + pow3(q46)
      q159 = 0.0625_dp*q12
      q160 = q1*y%d2val2_d1val3
      q161 = q84*y%d1val2_d1val3
      q162 = q42*x%d1val2_d1val3
      q163 = q98*y%d2val2
      q164 = q51*x%d2val2
      q165 = q4*x%d2val2_d1val3
      q166 = q116*x%d1val2_d1val3
      q167 = q55*x%d1val2
      q168 = 32.0_dp*y%d1val2
      q169 = q20*q55
      q170 = q45*y%d1val3
      q171 = q44*y%val
      q172 = -16.0_dp*q170 - 32.0_dp*q166 + 16.0_dp*q160 + 16.0_dp*q163 + 16.0_dp*q164 + 16.0_dp*q165 + 32.0_dp*q161 + 32.0_dp*q162 + q106*q171 - q167*q168 - q169*x%d2val2 + q46*(-8.0_dp*q56 + 8.0_dp*q50 + 8.0_dp*q52 + 8.0_dp*q53 + 8.0_dp*q54)
      q173 = 2.0_dp*q7 + 2.0_dp*q8
      q174 = 2.0_dp*y%d1val3
      q175 = q61*q71
      q176 = 2.0_dp*q175
      q177 = q1*y%d1val2_d2val3 - q116*x%d2val3 + q147*y%d1val2_d1val3 - q148*x%d1val2_d1val3 - q167*q174 + q176*q32 + q4*x%d1val2_d2val3 + q42*x%d2val3 + q60*x%d1val2_d1val3 - q62*y%d1val2 + q84*y%d2val3
      q178 = 0.015625_dp*q12
      q179 = q3*y%d2val1
      q180 = q179*x%d1val1_d1val3
      q181 = 64.0_dp*q98
      q182 = q3*x%d2val1
      q183 = q182*y%d1val1_d1val3
      q184 = 64.0_dp*q51
      q185 = 384.0_dp*x%d1val1_d1val3
      q186 = q31*y%d1val1
      q187 = q102*y%d2val1
      q188 = q17*y%val
      q189 = q188*x%d2val1
      q190 = q189*x%d1val1_d1val3
      q191 = 64.0_dp*q135
      q192 = q73*y%d1val3
      q193 = 384.0_dp*x%d1val3
      q194 = q193*y%val
      q195 = q141*x%d2val1
      q196 = q105*q71
      q197 = q18*q71
      q198 = q197*y%d1val1
      q199 = powm1(pow4(x%val))
      q200 = q199*q74
      q201 = q14*y%d1val1_d1val3
      q202 = q29*x%d1val1_d1val3
      q203 = q39*x%d1val1_d1val3
      q204 = q102*y%d1val1
      q205 = 8.0_dp*q71
      q206 = q205*x%d1val3
      q207 = q1*y%d2val1_d1val2_d1val3
      q208 = 32.0_dp*q3
      q209 = x%d1val1_d1val2*y%d1val1_d1val3
      q210 = x%d1val1_d1val3*y%d1val1_d1val2
      q211 = q84*y%d2val1_d1val3
      q212 = q179*x%d1val2_d1val3
      q213 = q98*y%d2val1_d1val2
      q214 = q182*y%d1val2_d1val3
      q215 = q51*x%d2val1_d1val2
      q216 = q4*x%d2val1_d1val2_d1val3
      q217 = q42*x%d2val1_d1val3
      q218 = 32.0_dp*x%d1val1_d1val2
      q219 = q31*x%d1val1_d1val3
      q220 = q31*x%d1val2_d1val3
      q221 = q102*y%d1val1_d1val2
      q222 = q188*x%d1val1_d1val3
      q223 = q17*y%d1val1
      q224 = x%d1val1_d1val3*x%d1val2
      q225 = q167*y%d2val1
      q226 = x%d1val2*y%d1val3
      q227 = q226*q90
      q228 = q17*q91
      q229 = q90*x%d1val2_d1val3
      q230 = q103*y%d1val2
      q231 = q19*y%d1val2_d1val3
      q232 = x%d1val1*y%val
      q233 = q232*x%d1val1_d1val2
      q234 = q71*x%d1val3
      q235 = 64.0_dp*q234
      q236 = x%d1val1*y%d1val1
      q237 = q236*x%d1val2
      q238 = q32*x%d2val1
      q239 = q18*q226
      q240 = q105*x%d1val2_d1val3
      q241 = q18*y%d1val2
      q242 = 96.0_dp*x%d1val3
      q243 = q199*q242
      q244 = q141*x%d1val3
      q245 = q119*q18
      q246 = -2.0_dp*q203 - 2.0_dp*q204 + 2.0_dp*q201 + 2.0_dp*q202 + q100 + q101 - q104 - q135*x%d2val1 + q245*x%d1val3 + q97 + q99
      q247 = q144*q41 + q246
      q248 = pow2(q144)
      q249 = q248 + q78
      q250 = 0.25_dp*q64
      q251 = pow2(q41)
      q252 = q3*y%d2val3
      q253 = q3*x%d2val3
      q254 = 2.0_dp*q14
      q255 = q3*x%d1val1_d1val3
      q256 = 4.0_dp*y%d1val1_d1val3
      q257 = 2.0_dp*q29
      q258 = 4.0_dp*x%d1val1_d1val3
      q259 = 2.0_dp*q39
      q260 = q31*x%d2val3
      q261 = 2.0_dp*q260
      q262 = 2.0_dp*x%d2val1_d1val3
      q263 = q188*x%d2val3
      q264 = pow2(x%d1val1_d1val3)
      q265 = 2.0_dp*q188
      q266 = q206*q232
      q267 = 4.0_dp*q175
      q268 = q197*y%d1val3
      q269 = 4.0_dp*x%d1val3
      q270 = q199*q61
      q271 = 6.0_dp*q270
      q272 = q1*y%d2val1_d2val3 - q102*q256 - q103*q174 - q105*q271 - q131*q258 - q135*q262 - q138*q258 + q147*y%d2val1_d1val3 + q149*x%d2val1 - q19*y%d2val3 + q236*q267 + q245*x%d2val3 + q252*x%d2val1 + q253*y%d2val1 + q254*y%d1val1_d2val3 + q255*q256 + q257*x%d1val1_d2val3 - q259*x%d1val1_d2val3 - q261*y%d1val1 - q263*x%d2val1 - q264*q265 + q266*x%d1val1_d1val3 + q268*q269 + q4*x%d2val1_d2val3 + q60*x%d2val1_d1val3 - q62*y%d2val1
      q273 = 2.0_dp*q251 + q144*q150 + q272
      q274 = q1*y%d1val1_d2val2_d1val3
      q275 = q14*y%d2val2_d1val3
      q276 = q3*y%d2val2
      q277 = q276*x%d1val1_d1val3
      q278 = q51*x%d1val1_d2val2
      q279 = q4*x%d1val1_d2val2_d1val3
      q280 = q98*y%d1val1_d2val2
      q281 = q3*x%d2val2
      q282 = q281*y%d1val1_d1val3
      q283 = q29*x%d2val2_d1val3
      q284 = q3*y%d1val2_d1val3
      q285 = 2.0_dp*q42
      q286 = 2.0_dp*q3
      q287 = x%d1val2_d1val3*y%d1val1_d1val2
      q288 = 2.0_dp*y%d1val2
      q289 = q102*y%d2val2
      q290 = q131*x%d2val2
      q291 = q39*x%d2val2_d1val3
      q292 = q17*q226
      q293 = q112*q188
      q294 = q112*q55
      q295 = 2.0_dp*q116
      q296 = q17*y%d1val2
      q297 = 2.0_dp*q296
      q298 = q188*x%d2val2
      q299 = 2.0_dp*q223
      q300 = x%d1val2*x%d1val2_d1val3
      q301 = 2.0_dp*q167
      q302 = q138*x%d2val2
      q303 = q45*y%d1val1_d1val3
      q304 = 4.0_dp*x%d1val2_d1val3
      q305 = x%d1val2*y%d1val2
      q306 = q141*q269
      q307 = x%d1val1*x%d1val3
      q308 = q119*x%d2val2
      q309 = 4.0_dp*x%d1val1_d1val2
      q310 = q234*q32
      q311 = q120*q71
      q312 = q119*q44
      q313 = q44*q71
      q314 = q313*y%d1val1
      q315 = q120*y%val
      q316 = q199*x%d1val3
      q317 = 6.0_dp*q316
      q318 = q55*q77
      q319 = q32*x%d1val1
      q320 = -2.0_dp*q33 + 2.0_dp*q26 + 2.0_dp*q27 + 2.0_dp*q28 + 2.0_dp*q30
      q321 = 0.5_dp*q320
      q322 = 0.5_dp*q144
      q323 = 0.5_dp*q173
      q324 = 2.0_dp*q57
      q325 = q11*q57
      q326 = q11*q173
      q327 = q173*q59
      q328 = 0.25_dp*q144
      q329 = q173*q63
      q330 = 2.0_dp*q255
      q331 = q286*x%d1val2_d1val3
      q332 = x%d1val2*y%d2val3
      q333 = 2.0_dp*x%d1val2_d1val3
      q334 = 2.0_dp*y%d1val2_d1val3
      q335 = x%d2val3*y%d1val2
      q336 = 2.0_dp*x%d1val1_d1val3
      q337 = q265*x%d1val2_d1val3
      q338 = q55*x%d1val1_d1val3
      q339 = x%d1val2*x%d2val3
      q340 = 2.0_dp*x%d2val3
      q341 = q77*x%d1val2_d1val3
      q342 = q176*y%d1val2
      q343 = q176*x%d1val2
      q344 = q1*y%d1val1_d1val2_d2val3 - q102*q334 + q112*q175*y%val - q116*x%d1val1_d2val3 - q131*q333 - q138*q333 + q14*y%d1val2_d2val3 + q142*q340 + q147*y%d1val1_d1val2_d1val3 - q148*x%d1val1_d1val2_d1val3 - q223*q339 + q226*q306 + q244*q341 + q252*x%d1val1_d1val2 + q253*y%d1val1_d1val2 + q258*q310 - q263*x%d1val1_d1val2 - q271*q319 - q288*q338 + q29*x%d1val2_d2val3 - q292*q336 - q294*y%d1val3 - q301*y%d1val1_d1val3 - q31*q332 - q31*q335 + q330*y%d1val2_d1val3 + q331*y%d1val1_d1val3 - q337*x%d1val1_d1val3 + q342*x%d1val1 + q343*y%d1val1 - q39*x%d1val2_d2val3 + q4*x%d1val1_d1val2_d2val3 + q42*x%d1val1_d2val3 + q60*x%d1val1_d1val2_d1val3 - q62*y%d1val1_d1val2 + q84*y%d1val1_d2val3
      q345 = q276*x%d1val2_d1val3
      q346 = q281*y%d1val2_d1val3
      q347 = q296*q300
      q348 = q167*y%d2val2
      q349 = q298*x%d1val2_d1val3
      q350 = q55*y%d1val2
      q351 = q32*q71
      q352 = q351*x%d2val2
      q353 = q171*q71
      q354 = q353*x%d1val2_d1val3
      q355 = q313*y%d1val2
      q356 = q155*q199
      q357 = 8.0_dp*y%d1val2
      q358 = -2.0_dp*q166 + 2.0_dp*q162 + q113*y%d1val2_d1val3 - q135*x%d2val2 + q160 + q163 + q164 + q165 - q167*q288 - q170 + q312*x%d1val3
      q359 = q173*q57 + q358
      q360 = pow2(q173)
      q361 = q157 + q360
      q362 = pow2(q57)
      q363 = 4.0_dp*q284
      q364 = 4.0_dp*y%d1val2_d1val3
      q365 = 2.0_dp*q335
      q366 = q174*q55
      q367 = pow2(x%d1val2_d1val3)
      q368 = q206*x%d1val2_d1val3
      q369 = q313*y%d1val3
      q370 = q1*y%d2val2_d2val3 + q113*y%d1val2_d2val3 + q147*y%d2val2_d1val3 - q148*x%d2val2_d1val3 + q149*x%d2val2 - q167*q364 - q17*q365*x%d1val2 - q171*q271 + q252*x%d2val2 + q253*y%d2val2 - q263*x%d2val2 - q265*q367 + q267*q305 + q269*q369 + q285*x%d1val2_d2val3 - q292*q304 - q295*x%d1val2_d2val3 - q304*q350 + q312*x%d2val3 + q32*q368 + q363*x%d1val2_d1val3 - q366*x%d2val2 + q4*x%d2val2_d2val3 - q45*y%d2val3 + q60*x%d2val2_d1val3 - q62*y%d2val2
      q371 = 2.0_dp*q362 + q173*q177 + q370
      q372 = 8.0_dp*q1
      q373 = 24.0_dp*q14
      q374 = 48.0_dp*y%d2val1_d1val3
      q375 = 24.0_dp*x%d1val1_d2val3
      q376 = 16.0_dp*q98
      q377 = q3*y%d1val1_d1val3
      q378 = 48.0_dp*x%d2val1_d1val3
      q379 = 24.0_dp*q29
      q380 = 8.0_dp*y%d3val1
      q381 = 8.0_dp*x%d3val1
      q382 = 16.0_dp*q51
      q383 = 8.0_dp*q4
      q384 = 96.0_dp*y%d1val1_d1val3
      q385 = 48.0_dp*x%d1val1_d2val3
      q386 = q31*y%d2val3
      q387 = 24.0_dp*x%d2val1
      q388 = 48.0_dp*y%d2val1
      q389 = x%d1val1_d1val3*y%d1val3
      q390 = 16.0_dp*q55*y%d1val3
      q391 = q223*x%d2val3
      q392 = 24.0_dp*q19
      q393 = q71*q73
      q394 = 16.0_dp*y%d2val3
      q395 = q234*x%d1val1_d1val3
      q396 = q242*y%d1val3
      q397 = q242*q71
      q398 = 48.0_dp*x%d2val3
      q399 = q232*x%d2val1
      q400 = x%d1val1_d1val3*y%val
      q401 = q264*q71
      q402 = q175*x%d1val1
      q403 = q175*y%d1val1
      q404 = q175*q20
      q405 = 144.0_dp*q270
      q406 = q316*x%d1val1_d1val3
      q407 = q61*powm1(pow5(x%val))
      q408 = 192.0_dp*q407
      q409 = 2.0_dp*q1
      q410 = 2.0_dp*q4
      q411 = 4.0_dp*q15 - q118*q19 + q409*y%d2val1 + q410*x%d2val1
      q412 = 6.0_dp*q144
      q413 = 4.0_dp*q98
      q414 = 12.0_dp*x%d2val1_d1val3
      q415 = 4.0_dp*q51
      q416 = 24.0_dp*x%d1val1_d1val3
      q417 = 12.0_dp*x%d2val1
      q418 = 24.0_dp*x%d1val3
      q419 = 4.0_dp*q11
      q420 = 0.125_dp*q0
      q421 = q112*q3
      q422 = 4.0_dp*y%d1val1_d1val2_d1val3
      q423 = q286*y%d1val1_d1val2
      q424 = 4.0_dp*x%d1val1_d1val2_d1val3
      q425 = q288*q31
      q426 = q17*q309
      q427 = 4.0_dp*q55
      q428 = q17*x%d1val1_d1val3
      q429 = q428*q77
      q430 = q17*q224
      q431 = q258*x%d1val2_d1val3
      q432 = x%d1val1_d2val3*x%d1val2
      q433 = q17*q339
      q434 = q205*y%d1val3
      q435 = q77*x%d1val1
      q436 = q435*x%d1val1_d1val2
      q437 = q205*x%d1val1_d1val3
      q438 = q226*x%d1val1
      q439 = q437*x%d1val2_d1val3
      q440 = q206*x%d1val1
      q441 = 4.0_dp*q71
      q442 = q319*q441
      q443 = q440*x%d1val2
      q444 = q206*x%d1val1_d1val2
      q445 = q206*q224
      q446 = q269*q71
      q447 = q238*q71
      q448 = 4.0_dp*y%d1val1_d1val2
      q449 = 2.0_dp*q197
      q450 = 12.0_dp*q270
      q451 = 12.0_dp*q316
      q452 = 6.0_dp*q199*x%d2val3
      q453 = 24.0_dp*q407
      q454 = 0.125_dp*q249
      q455 = 4.0_dp*q143
      q456 = -0.5_dp*q295*x%d2val1 - 0.5_dp*q309*q39 - 2.0_dp*q89*y%d1val1 + 0.5_dp*q14*q448 + 0.5_dp*q144*q320 + 0.5_dp*q29*q309 + 0.5_dp*q409*y%d2val1_d1val2 + 0.5_dp*q410*x%d2val1_d1val2 + 0.5_dp*q441*q94 + q85 + q86 - q92
      q457 = q3*x%d1val2_d1val3
      q458 = 2.0_dp*x%d2val2_d1val3
      q459 = q112*q17
      q460 = x%d1val2_d1val3*y%d1val3
      q461 = q55*y%d1val2_d1val3
      q462 = q55*y%d2val2
      q463 = q17*q300
      q464 = x%d1val2*x%d1val2_d2val3
      q465 = q367*q71
      q466 = q175*y%d1val2
      q467 = q175*x%d1val2
      q468 = q199*q418
      q469 = q44*y%d1val1
      q470 = 2.0_dp*q359
      q471 = 4.0_dp*q320
      q472 = 8.0_dp*x%d1val1_d1val2
      q473 = 8.0_dp*x%d1val1_d1val2_d1val3
      q474 = q188*x%d1val2_d1val3
      q475 = 24.0_dp*q84
      q476 = 48.0_dp*y%d2val2_d1val3
      q477 = 24.0_dp*x%d1val2_d2val3
      q478 = 48.0_dp*x%d2val2_d1val3
      q479 = 24.0_dp*q42
      q480 = 8.0_dp*y%d3val2
      q481 = 8.0_dp*x%d3val2
      q482 = 48.0_dp*q296
      q483 = 24.0_dp*x%d2val2
      q484 = q17*q483
      q485 = 48.0_dp*x%d2val2
      q486 = 24.0_dp*q45
      q487 = q155*q71
      q488 = q397*x%d2val2
      q489 = q356*y%val
      q490 = 4.0_dp*q43 - q118*q45 + q409*y%d2val2 + q410*x%d2val2
      q491 = 6.0_dp*q173
      q492 = 12.0_dp*x%d2val2_d1val3
      q493 = 12.0_dp*x%d2val2
      binary%val = q0
      binary%d1val1 = q6
      binary%d1val2 = q10
      binary%d1val3 = q12
      binary%d2val1 = q23*q24
      binary%d1val1_d1val2 = q0*q34 + q25
      binary%d1val1_d1val3 = q0*q41 + q11*q6
      binary%d2val2 = q24*q48
      binary%d1val2_d1val3 = q49 + q58
      binary%d2val3 = q0*q64
      binary%d3val1 = q79*q80
      binary%d2val1_d1val2 = q10*q81 + q24*q96
      binary%d2val1_d1val3 = q107*q24 + q12*q81
      binary%d1val1_d2val2 = q0*q122
      binary%d1val1_d1val2_d1val3 = q0*q143 + q10*q41 + q11*q25 + q12*q34 + q57*q6
      binary%d1val1_d2val3 = q0*(q144*q145 + q146*q41 + q150)
      binary%d3val2 = q158*q80
      binary%d2val2_d1val3 = q159*q48 + q172*q24
      binary%d1val2_d2val3 = q0*(q145*q173 + q146*q57 + q177)
      binary%d3val1_d1val3 = q178*q79 + q80*(-192.0_dp*q187 - 192.0_dp*q190 - 192.0_dp*q39*x%d2val1_d1val3 + 192.0_dp*q180 + 192.0_dp*q183 - q131*q69 - q138*q69 + q181*y%d3val1 + q184*x%d3val1 - q185*q186 + q185*q196 - q191*x%d3val1 + q192*q72 + q193*q198 - q193*q200 + q194*q195 + q22*(-12.0_dp*q40 + 12.0_dp*q35 + 12.0_dp*q36 + 12.0_dp*q37 + 12.0_dp*q38) + q65*y%d3val1_d1val3 + q66*y%d2val1_d1val3 + q67*x%d2val1_d1val3 + q68*x%d3val1_d1val3 - q70*y%d1val1_d1val3 + q78*(-48.0_dp*q40 + 48.0_dp*q35 + 48.0_dp*q36 + 48.0_dp*q37 + 48.0_dp*q38) + (48.0_dp*q2 + 48.0_dp*q5)*(-4.0_dp*q104 - 8.0_dp*q203 - 8.0_dp*q204 + 4.0_dp*q100 + 4.0_dp*q101 + 4.0_dp*q97 + 4.0_dp*q99 + 8.0_dp*q201 + 8.0_dp*q202 - q103*q77 + q105*q206))
      binary%d2val1_d1val2_d1val3 = 0.0625_dp*q10*q107 + q159*q96 + q24*(-16.0_dp*q225 - 16.0_dp*q227 - 16.0_dp*q230 - 16.0_dp*q231 - 32.0_dp*q221 - 32.0_dp*q223*q224 - 32.0_dp*q89*y%d1val1_d1val3 + 16.0_dp*q207 + 16.0_dp*q211 + 16.0_dp*q212 + 16.0_dp*q213 + 16.0_dp*q214 + 16.0_dp*q215 + 16.0_dp*q216 + 16.0_dp*q217 + 64.0_dp*q142*x%d1val1_d1val3 + q106*q238 + q106*q241 - q131*q218 - q138*q218 - q168*q219 - q169*x%d2val1_d1val2 - q20*q229 + q208*q209 + q208*q210 + q21*(-8.0_dp*q132 - 8.0_dp*q133 - 8.0_dp*q134 - 8.0_dp*q136 - 8.0_dp*q137 - 8.0_dp*q139 + 8.0_dp*q123 + 8.0_dp*q124 + 8.0_dp*q125 + 8.0_dp*q126 + 8.0_dp*q127 + 8.0_dp*q128 + 8.0_dp*q129 + 8.0_dp*q130 + q244*q91) - q218*q222 - q220*q88 - q228*x%d2val1_d1val3 + q233*q235 + q235*q237 + q239*q93 + q240*q93 - q243*q94 + q82*y%d1val1_d1val2_d1val3 + q83*x%d1val1_d1val2_d1val3 - q87*x%d1val1_d1val2_d1val3 + q95*(-4.0_dp*q40 + 4.0_dp*q35 + 4.0_dp*q36 + 4.0_dp*q37 + 4.0_dp*q38)) + q49*q81 + q58*q81
      binary%d2val1_d2val3 = q0*(q146*q247 + q249*q250 + q273)
      binary%d1val1_d2val2_d1val3 = q0*(q112*q284 - q112*q292 + q113*y%d1val1_d1val2_d1val3 - q114*y%d1val2_d1val3 + q121*(-0.5_dp*q56 + 0.5_dp*q50 + 0.5_dp*q52 + 0.5_dp*q53 + 0.5_dp*q54) - q135*x%d1val1_d2val2 + q140*q314 + q142*q304 + q172*(0.0625_dp*q2 + 0.0625_dp*q5) + q174*q311 - q220*q288 - q224*q297 + q274 + q275 + q277 + q278 + q279 + q280 + q282 + q283 + q285*x%d1val1_d1val2_d1val3 + q286*q287 - q289 - q290 - q291 - q293*x%d1val2_d1val3 - q294*y%d1val2 - q295*x%d1val1_d1val2_d1val3 - q298*x%d1val1_d1val3 - q299*q300 - q301*y%d1val1_d1val2 - q302 - q303 + q305*q306 + q307*q308 + q309*q310 + q312*x%d1val1_d1val3 - q315*q317 + q48*(-0.0625_dp*q40 + 0.0625_dp*q35 + 0.0625_dp*q36 + 0.0625_dp*q37 + 0.0625_dp*q38) + (0.5_dp*q7 + 0.5_dp*q8)*(-4.0_dp*q132 - 4.0_dp*q133 - 4.0_dp*q134 - 4.0_dp*q137 - 4.0_dp*q139 + 4.0_dp*q123 + 4.0_dp*q124 + 4.0_dp*q125 + 4.0_dp*q126 + 4.0_dp*q127 + 4.0_dp*q128 + 4.0_dp*q129 + 4.0_dp*q130 + q206*q319 - q318*x%d1val1_d1val2)) + q12*q122
      binary%d1val1_d1val2_d2val3 = q0*(q143*q146 + q144*q325 + q150*q323 + q177*q322 + q321*q59 + q321*q63 + q324*q41 + q326*q41 + q327*q328 + q328*q329 + q344)
      binary%d3val2_d1val3 = q158*q178 + q80*(-192.0_dp*q116*x%d2val2_d1val3 - 192.0_dp*q348 - 192.0_dp*q349 - 384.0_dp*q347 + 192.0_dp*q345 + 192.0_dp*q346 + 384.0_dp*q354 + q151*y%d2val2_d1val3 + q152*x%d2val2_d1val3 - q153*q292 - q153*q350 - q154*y%d1val2_d1val3 + q156*y%d1val3 + q157*(-48.0_dp*q56 + 48.0_dp*q50 + 48.0_dp*q52 + 48.0_dp*q53 + 48.0_dp*q54) + q181*y%d3val2 + q184*x%d3val2 - q191*x%d3val2 + q193*q352 + q193*q355 - q194*q356 + q47*(-12.0_dp*q56 + 12.0_dp*q50 + 12.0_dp*q52 + 12.0_dp*q53 + 12.0_dp*q54) + q65*y%d3val2_d1val3 + q68*x%d3val2_d1val3 + (48.0_dp*q7 + 48.0_dp*q8)*(-4.0_dp*q170 - 8.0_dp*q166 + 4.0_dp*q160 + 4.0_dp*q163 + 4.0_dp*q164 + 4.0_dp*q165 + 8.0_dp*q161 + 8.0_dp*q162 - q167*q357 + q171*q206 - q318*x%d2val2))
      binary%d2val2_d2val3 = q0*(q146*q359 + q250*q361 + q371)
      binary%d3val1_d2val3 = q420*(-24.0_dp*q260*y%d2val1 - 24.0_dp*q39*x%d2val1_d2val3 - 288.0_dp*q105*q406 - 48.0_dp*q103*y%d1val1_d1val3 - 48.0_dp*q223*q264 - 48.0_dp*q389*q90 + 12.0_dp*q144*q272 + 12.0_dp*q150*q411 + 192.0_dp*q236*q395 + 24.0_dp*q144*q251 + 24.0_dp*q182*y%d1val1_d2val3 + 48.0_dp*q246*q41 + 48.0_dp*q403*x%d2val1 + 6.0_dp*q150*q248 + 96.0_dp*q232*q401 + 96.0_dp*q268*x%d1val1_d1val3 - q102*q374 - q131*q378 - q138*q378 - q169*x%d3val1_d1val3 + q179*q375 - q18*q405*y%d1val1 - q186*q385 - q189*q375 - q192*q243 + q195*q396 + q196*q385 + q197*q384*x%d1val3 + q198*q398 - q200*q398 - q219*q384 - q222*q378 + q232*q397*x%d2val1_d1val3 + q252*q381 + q253*q380 + q255*q374 - q263*q381 - q338*q388 + q372*y%d3val1_d2val3 + q373*y%d2val1_d2val3 + q376*y%d3val1_d1val3 + q377*q378 + q379*x%d2val1_d2val3 - q380*q62 + q382*x%d3val1_d1val3 + q383*x%d3val1_d2val3 - q386*q387 - q387*q391 + q388*q402 - q390*x%d3val1 - q392*y%d1val1_d2val3 + q393*q394 + q397*q400*x%d2val1 + q398*q399*q71 - q399*q405 + q404*x%d3val1 + q408*q74 + q419*(-12.0_dp*q187 - 12.0_dp*q19*y%d1val1_d1val3 - 12.0_dp*q190 + 12.0_dp*q14*y%d2val1_d1val3 + 12.0_dp*q180 + 12.0_dp*q183 + 3.0_dp*q248*q41 + 6.0_dp*q41*q411 - q131*q417 - q138*q417 - q186*q416 + q192*q205 + q196*q416 + q198*q418 - q200*q418 + q232*q234*q387 + q246*q412 + q29*q414 - q318*x%d3val1 - q39*q414 + q413*y%d3val1 + q415*x%d3val1 + q75*y%d3val1_d1val3 + q76*x%d3val1_d1val3) + q64*(q20*q393 + q372*y%d3val1 + q373*y%d2val1 + q379*x%d2val1 + q381*q4 - q387*q39 - q392*y%d1val1 + q411*q412 + pow3(q144)))
      binary%d2val1_d1val2_d2val3 = q0*(-2.0_dp*q186*x%d1val2_d2val3 + 0.25_dp*q177*q249 + 0.5_dp*q249*q325 + 4.0_dp*q32*q401 + q1*y%d2val1_d1val2_d2val3 - q102*q422 - q103*q334 - q112*q386 - q112*q391 - q114*y%d1val1_d2val3 - q116*x%d2val1_d2val3 - q131*q424 - q138*q424 + q144*q344 + q146*(-2.0_dp*q221 - q112*q131 - q112*q138 - q112*q222 + q112*q377 - q114*y%d1val1_d1val3 - q116*x%d2val1_d1val3 - q135*x%d2val1_d1val2 + q140*q197*y%d1val2 + q140*q447 + q142*q258 + q143*q144 - q186*q333 - q189*x%d1val2_d1val3 + q207 + q210*q286 + q211 + q212 + q213 + q214 + q215 + q216 + q217 - q219*q288 - q224*q299 - q225 + q226*q449 - q227 - q230 - q231 + q234*q436 + q237*q446 + q245*x%d1val2_d1val3 + q254*y%d1val1_d1val2_d1val3 + q257*x%d1val1_d1val2_d1val3 - q259*x%d1val1_d1val2_d1val3 - q317*q94 + q320*q41) + q147*y%d2val1_d1val2_d1val3 - q148*x%d2val1_d1val2_d1val3 + q149*x%d2val1_d1val2 + q150*q320 - q174*q229 + q179*x%d1val2_d2val3 + q182*y%d1val2_d2val3 - q189*x%d1val2_d2val3 - q19*y%d1val2_d2val3 + q197*q269*y%d1val2_d1val3 + q197*q365 - q209*q427 - q210*q427 - q219*q364 - q220*q256 - q223*q431 + q226*q446*x%d2val1 + q232*q439 - q233*q450 + q234*q341*x%d2val1 + q236*q339*q441 + q236*q368 - q237*q450 - q238*q271 - q239*q451 - q240*q451 - q241*q271 + q245*x%d1val2_d2val3 + q247*q324 + q247*q326 + q252*x%d2val1_d1val2 + q253*y%d2val1_d1val2 + q254*y%d1val1_d1val2_d2val3 + q255*q422 + q256*q3*x%d1val1_d1val2_d1val3 - q256*q430 + q257*x%d1val1_d1val2_d2val3 - q259*x%d1val1_d1val2_d2val3 - q261*y%d1val1_d1val2 + q262*q284 - q262*q292 - q262*q350 - q263*x%d2val1_d1val2 - q264*q297 + q266*x%d1val1_d1val2_d1val3 + q268*q304 + q269*q351*x%d2val1_d1val3 + q273*q323 - q293*x%d1val1_d2val3 - q299*q432 - q301*y%d2val1_d1val3 + q307*q434*x%d1val1_d1val2 + q309*q403 - q316*q319*q416 + q327*q454 + q329*q454 + q331*y%d2val1_d1val3 + q332*q449 - q332*q90 - q333*q55*y%d2val1 - q335*q90 - q337*x%d2val1_d1val3 + q340*q447 + q342*x%d2val1 + q343*y%d2val1 - q366*x%d2val1_d1val2 - q389*q426 + q4*x%d2val1_d1val2_d2val3 + q400*q444 + q402*q448 + q41*q455 + q42*x%d2val1_d2val3 + q421*y%d1val1_d2val3 + q423*x%d1val1_d2val3 - q425*x%d1val1_d2val3 - q429*x%d1val1_d1val2_d1val3 - q433*y%d2val1 + q436*q71*x%d2val3 + q437*q438 + q440*x%d1val1_d1val3*y%d1val2 + q442*x%d1val1_d2val3 + q443*y%d1val1_d1val3 + q445*y%d1val1 - q452*q94 + q453*q94 + q456*q59 + q456*q63 + q60*x%d2val1_d1val2_d1val3 - q62*y%d2val1_d1val2 + q84*y%d2val1_d2val3)
      binary%d1val1_d2val2_d2val3 = q0*(-12.0_dp*q171*q406 - 2.0_dp*q102*y%d2val2_d1val3 - 2.0_dp*q433*y%d1val1_d1val2 - 2.0_dp*q55*x%d2val2*y%d1val1_d1val3 + 0.125_dp*q64*(-16.0_dp*q89*y%d1val2 - 8.0_dp*q115 - 8.0_dp*q117 + 16.0_dp*q42*x%d1val1_d1val2 + 16.0_dp*q84*y%d1val1_d1val2 + 8.0_dp*q108 + 8.0_dp*q109 + 8.0_dp*q110 + 8.0_dp*q111 + q144*q361 + q173*q471 + q20*q311 - q228*x%d1val1_d1val2) + 0.25_dp*q150*q361 + 0.5_dp*q11*(-4.0_dp*q289 - 4.0_dp*q290 - 4.0_dp*q291 - 4.0_dp*q302 - 4.0_dp*q303 - 8.0_dp*q167*y%d1val1_d1val2 - 8.0_dp*q223*q300 - 8.0_dp*q224*q296 - 8.0_dp*q89*y%d1val2_d1val3 + 16.0_dp*q244*q305 + 4.0_dp*q274 + 4.0_dp*q275 + 4.0_dp*q277 + 4.0_dp*q278 + 4.0_dp*q279 + 4.0_dp*q280 + 4.0_dp*q282 + 4.0_dp*q283 + 8.0_dp*q287*q3 + 8.0_dp*q84*y%d1val1_d1val2_d1val3 - q116*q473 + q120*q434 + q141*q91*x%d1val2_d1val3 + q144*q470 + q171*q437 + q173*q455 + q206*q469 - q220*q357 + q234*q91*x%d1val1_d1val2 + q266*x%d2val2 + q284*q472 - q292*q472 - q315*q468 - q318*x%d1val1_d2val2 - q350*q472 + q361*q41 + q42*q473 - q429*x%d2val2 + q471*q57 - q472*q474) + 2.0_dp*q311*y%d2val3 + q1*y%d1val1_d2val2_d2val3 + q113*y%d1val1_d1val2_d2val3 - q114*y%d1val2_d2val3 - q120*q451*y%d1val3 - q131*q458 - q138*q458 + q14*y%d2val2_d2val3 + q147*y%d1val1_d2val2_d1val3 - q148*x%d1val1_d2val2_d1val3 + q149*x%d1val1_d2val2 - q167*q422 - q17*q341*x%d1val1_d1val2_d1val3 + q173*q344 - q174*q428*x%d2val2 + q176*x%d1val1*y%d2val2 + q176*x%d2val2*y%d1val1 + q177*q320 + q205*q438*x%d1val2_d1val3 + q206*q300*y%d1val1 + q206*q32*x%d1val1_d1val2_d1val3 - q220*q364 + q226*q444 - q232*q271*x%d2val2 + q234*q435*x%d2val2_d1val3 + q252*x%d1val1_d2val2 + q253*y%d1val1_d2val2 + q256*q313*x%d1val3 - q256*q463 + q258*q369 - q260*y%d2val2 - q263*x%d1val1_d2val2 - q265*x%d1val1_d1val3*x%d2val2_d1val3 - q271*q469 + q276*x%d1val1_d2val3 + q281*y%d1val1_d2val3 + q285*x%d1val1_d1val2_d2val3 + q286*x%d2val2_d1val3*y%d1val1_d1val3 - q287*q427 + q29*x%d2val2_d2val3 - q292*q424 - q293*x%d1val2_d2val3 - q295*x%d1val1_d1val2_d2val3 - q296*q431 - q297*q432 - q298*x%d1val1_d2val3 - q299*q367 - q299*q464 - q305*q450*x%d1val1 + q306*x%d2val2*y%d1val3 + q308*x%d1val1*x%d2val3 + q309*q351*x%d2val3 - q309*q461 + q309*q466 + q312*x%d1val1_d2val3 + q314*q340 - q315*q452 + q315*q453 - q319*q468*x%d1val2_d1val3 + q32*q439 - q32*q450*x%d1val1_d1val2 + q322*q371 + q330*y%d2val2_d1val3 - q332*q459 + q335*q441*x%d1val1*x%d1val2 - q335*q459 - q336*q462 - q350*q424 + q363*x%d1val1_d1val2_d1val3 - q364*q430 - q366*x%d1val1_d2val2 + q368*x%d1val1*y%d1val2 + q368*x%d1val1_d1val2*y%val - q386*x%d2val2 - q39*x%d2val2_d2val3 - q391*x%d2val2 + q395*q77*x%d2val2 + q4*x%d1val1_d2val2_d2val3 + q41*q470 + q421*y%d1val2_d2val3 + q422*q457 + q423*x%d1val2_d2val3 - q425*x%d1val2_d2val3 - q426*q460 + q435*q465 + q442*x%d1val2_d2val3 + q443*y%d1val2_d1val3 + q445*y%d1val2 + q448*q467 - q45*y%d1val1_d2val3 + q455*q57 + q60*x%d1val1_d2val2_d1val3 - q62*y%d1val1_d2val2)
      binary%d3val2_d2val3 = q420*(-24.0_dp*q116*x%d2val2_d2val3 - 24.0_dp*q433*y%d2val2 - 288.0_dp*q171*q316*x%d1val2_d1val3 - 48.0_dp*q462*x%d1val2_d1val3 - 96.0_dp*q463*y%d1val2_d1val3 + 12.0_dp*q173*q370 + 12.0_dp*q177*q490 + 192.0_dp*q234*q300*y%d1val2 + 24.0_dp*q173*q362 + 24.0_dp*q281*y%d1val2_d2val3 + 48.0_dp*q313*q335 + 48.0_dp*q353*x%d1val2_d2val3 + 48.0_dp*q358*q57 + 48.0_dp*q467*y%d2val2 + 6.0_dp*q177*q360 + 96.0_dp*q32*q465 + 96.0_dp*q369*x%d1val2_d1val3 + q155*q408*y%val - q167*q476 - q169*x%d3val2_d1val3 - q17*q460*q485 + q226*q488 + q242*q313*y%d1val2_d1val3 + q242*q351*x%d2val2_d1val3 + q252*q481 + q253*q480 - q263*q481 + q276*q477 + q284*q478 - q292*q478 - q298*q477 - q32*q405*x%d2val2 - q332*q484 - q335*q484 - q350*q478 + q352*q398 - q356*q396 - q367*q482 + q372*y%d3val2_d2val3 + q376*y%d3val2_d1val3 + q382*x%d3val2_d1val3 + q383*x%d3val2_d2val3 - q390*x%d3val2 + q394*q487 - q398*q489 + q404*x%d3val2 - q405*q44*y%d1val2 + q419*(-12.0_dp*q348 - 12.0_dp*q349 - 12.0_dp*q45*y%d1val2_d1val3 - 24.0_dp*q347 + 12.0_dp*q345 + 12.0_dp*q346 + 12.0_dp*q84*y%d2val2_d1val3 + 24.0_dp*q354 + 3.0_dp*q360*q57 + 6.0_dp*q490*q57 - q116*q492 + q155*q434 - q292*q493 - q318*x%d3val2 - q350*q493 + q352*q418 + q355*q418 + q358*q491 + q413*y%d3val2 + q415*x%d3val2 - q418*q489 + q42*q492 + q75*y%d3val2_d1val3 + q76*x%d3val2_d1val3) + q457*q476 - q461*q485 - q464*q482 + q466*q485 - q474*q478 + q475*y%d2val2_d2val3 + q479*x%d2val2_d2val3 - q480*q62 - q486*y%d1val2_d2val3 + q488*x%d1val2_d1val3*y%val + q64*(-q116*q483 + q20*q487 + q372*y%d3val2 + q383*x%d3val2 + q475*y%d2val2 + q479*x%d2val2 - q486*y%d1val2 + q490*q491 + pow3(q173)))
   end function pow_self
   
   function pow_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q319
      real(dp) :: q318
      real(dp) :: q317
      real(dp) :: q316
      real(dp) :: q315
      real(dp) :: q314
      real(dp) :: q313
      real(dp) :: q312
      real(dp) :: q311
      real(dp) :: q310
      real(dp) :: q309
      real(dp) :: q308
      real(dp) :: q307
      real(dp) :: q306
      real(dp) :: q305
      real(dp) :: q304
      real(dp) :: q303
      real(dp) :: q302
      real(dp) :: q301
      real(dp) :: q300
      real(dp) :: q299
      real(dp) :: q298
      real(dp) :: q297
      real(dp) :: q296
      real(dp) :: q295
      real(dp) :: q294
      real(dp) :: q293
      real(dp) :: q292
      real(dp) :: q291
      real(dp) :: q290
      real(dp) :: q289
      real(dp) :: q288
      real(dp) :: q287
      real(dp) :: q286
      real(dp) :: q285
      real(dp) :: q284
      real(dp) :: q283
      real(dp) :: q282
      real(dp) :: q281
      real(dp) :: q280
      real(dp) :: q279
      real(dp) :: q278
      real(dp) :: q277
      real(dp) :: q276
      real(dp) :: q275
      real(dp) :: q274
      real(dp) :: q273
      real(dp) :: q272
      real(dp) :: q271
      real(dp) :: q270
      real(dp) :: q269
      real(dp) :: q268
      real(dp) :: q267
      real(dp) :: q266
      real(dp) :: q265
      real(dp) :: q264
      real(dp) :: q263
      real(dp) :: q262
      real(dp) :: q261
      real(dp) :: q260
      real(dp) :: q259
      real(dp) :: q258
      real(dp) :: q257
      real(dp) :: q256
      real(dp) :: q255
      real(dp) :: q254
      real(dp) :: q253
      real(dp) :: q252
      real(dp) :: q251
      real(dp) :: q250
      real(dp) :: q249
      real(dp) :: q248
      real(dp) :: q247
      real(dp) :: q246
      real(dp) :: q245
      real(dp) :: q244
      real(dp) :: q243
      real(dp) :: q242
      real(dp) :: q241
      real(dp) :: q240
      real(dp) :: q239
      real(dp) :: q238
      real(dp) :: q237
      real(dp) :: q236
      real(dp) :: q235
      real(dp) :: q234
      real(dp) :: q233
      real(dp) :: q232
      real(dp) :: q231
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = pow(x%val, y)
      q1 = powm1(x%val)
      q2 = q1*x%d1val1
      q3 = q0*y
      q4 = q1*x%d1val2
      q5 = q1*q3
      q6 = pow2(x%d1val1)
      q7 = 4.0_dp*q1
      q8 = q6*q7
      q9 = 4.0_dp*x%d2val1 + q8*y - q8
      q10 = 0.25_dp*q5
      q11 = q1*x%d1val1_d1val2
      q12 = powm1(pow2(x%val))
      q13 = q12*x%d1val2
      q14 = q13*x%d1val1
      q15 = pow2(y)
      q16 = q0*q15
      q17 = q12*x%d1val3
      q18 = q17*x%d1val1
      q19 = pow2(x%d1val2)
      q20 = q19*q7
      q21 = -q20
      q22 = 4.0_dp*x%d2val2 + q20*y + q21
      q23 = q13*x%d1val3
      q24 = pow2(x%d1val3)
      q25 = q1*q24
      q26 = q25*y - q25 + x%d2val3
      q27 = 12.0_dp*q2
      q28 = q27*x%d2val1
      q29 = pow3(x%d1val1)
      q30 = 8.0_dp*q12
      q31 = 12.0_dp*q12
      q32 = q29*y
      q33 = 4.0_dp*q12
      q34 = q29*q33
      q35 = 4.0_dp*x%d3val1 + q15*q34 + q28*y - q28 + q29*q30 - q31*q32
      q36 = 0.25_dp*q9
      q37 = q13*q36
      q38 = 8.0_dp*q2
      q39 = q38*x%d1val1_d1val2
      q40 = q12*q6
      q41 = q40*x%d1val2
      q42 = 4.0_dp*q41
      q43 = 4.0_dp*x%d2val1_d1val2 + q39*y - q39 - q42*y + q42
      q44 = q17*q36
      q45 = q38*x%d1val1_d1val3
      q46 = q40*x%d1val3
      q47 = 4.0_dp*q46
      q48 = 4.0_dp*x%d2val1_d1val3 + q45*y - q45 - q47*y + q47
      q49 = 8.0_dp*q11
      q50 = q49*x%d1val2
      q51 = q30*x%d1val1
      q52 = q19*y
      q53 = 2.0_dp*x%d2val2
      q54 = q21 + q53
      q55 = 2.0_dp*q1
      q56 = q55*x%d1val1
      q57 = q2*y
      q58 = 4.0_dp*x%d1val1_d2val2 + q22*q57 + q50*y - q50 - q51*q52 - q54*q56
      q59 = x%d1val1*x%d1val2_d1val3
      q60 = q12*q59
      q61 = q17*x%d1val1_d1val2
      q62 = q13*x%d1val1_d1val3
      q63 = powm1(pow3(x%val))
      q64 = q63*x%d1val1
      q65 = q64*x%d1val3
      q66 = q65*x%d1val2
      q67 = x%d1val2*pow3(y)
      q68 = q0*q67
      q69 = q55*x%d1val3
      q70 = q69*x%d1val1_d1val3
      q71 = q12*q24
      q72 = q71*y
      q73 = 2.0_dp*q72
      q74 = q1*(-2.0_dp*q25 + x%d2val3)
      q75 = 12.0_dp*q4
      q76 = q75*x%d2val2
      q77 = pow3(x%d1val2)
      q78 = q77*y
      q79 = q33*q77
      q80 = 4.0_dp*x%d3val2 + q15*q79 + q30*q77 - q31*q78 + q76*y - q76
      q81 = 0.25_dp*q17
      q82 = q22*q81
      q83 = q12*q19
      q84 = q83*x%d1val3
      q85 = 4.0_dp*q84
      q86 = 8.0_dp*q4
      q87 = q86*x%d1val2_d1val3
      q88 = q85 - q87
      q89 = 4.0_dp*x%d2val2_d1val3 - q85*y + q87*y + q88
      q90 = q69*x%d1val2_d1val3
      q91 = q26*y
      q92 = q35*q81
      q93 = q27*x%d2val1_d1val3
      q94 = q1*x%d1val1_d1val3
      q95 = q94*x%d2val1
      q96 = 12.0_dp*q95
      q97 = q31*x%d1val1
      q98 = q97*x%d1val3*x%d2val1
      q99 = q40*x%d1val1_d1val3
      q100 = q29*q63
      q101 = 16.0_dp*x%d1val3
      q102 = 36.0_dp*y
      q103 = q32*q63
      q104 = 24.0_dp*x%d1val3
      q105 = 12.0_dp*q15
      q106 = 8.0_dp*x%d1val3
      q107 = q100*q15
      q108 = q12*q36*x%d1val2_d1val3
      q109 = q63*x%d1val3
      q110 = q109*q9*x%d1val2
      q111 = 0.25_dp*q13*q48
      q112 = q43*q81
      q113 = q38*x%d1val1_d1val2_d1val3
      q114 = q49*x%d1val1_d1val3
      q115 = x%d1val1_d1val2*x%d1val3
      q116 = q115*q51
      q117 = x%d1val1_d1val3*x%d1val2
      q118 = q117*q51
      q119 = q40*x%d1val2_d1val3
      q120 = 4.0_dp*q119
      q121 = q6*q63
      q122 = q121*x%d1val2
      q123 = q106*q122
      q124 = q56*x%d1val1_d1val3
      q125 = q124*y - q124 - q46*y + q46 + x%d2val1_d1val3
      q126 = q125*q69
      q127 = q55*q6
      q128 = 2.0_dp*x%d2val1 + q127*y - q127
      q129 = 0.5_dp*q74
      q130 = 0.5_dp*q1
      q131 = q128*q130
      q132 = q56*x%d1val1_d2val3
      q133 = pow2(x%d1val1_d1val3)
      q134 = q133*q55
      q135 = q40*x%d2val3
      q136 = q33*x%d1val3
      q137 = x%d1val1*x%d1val1_d1val3
      q138 = q136*q137
      q139 = q24*q63
      q140 = 2.0_dp*q139
      q141 = q140*q6
      q142 = q132*y - q132 + q134*y - q134 - q135*y + q135 - q138*y + q138 + q141*y - q141 + x%d2val1_d2val3
      q143 = q58*q81
      q144 = q49*x%d1val2_d1val3
      q145 = q86*x%d1val1_d1val2_d1val3
      q146 = q30*x%d1val2
      q147 = q115*q146
      q148 = x%d1val2_d1val3*y
      q149 = q52*x%d1val1_d1val3
      q150 = q55*x%d1val1_d1val3
      q151 = q54*x%d1val1
      q152 = 2.0_dp*q12
      q153 = q152*x%d1val3
      q154 = q22*y
      q155 = q2*x%d1val2_d2val3
      q156 = q11*x%d2val3
      q157 = q69*x%d1val1_d1val2_d1val3
      q158 = q150*x%d1val2_d1val3
      q159 = q4*x%d1val1_d2val3
      q160 = q152*x%d1val2
      q161 = x%d1val1*x%d2val3
      q162 = q71*x%d1val1_d1val2
      q163 = q13*x%d2val3
      q164 = x%d1val1*y
      q165 = q59*y
      q166 = 6.0_dp*q12
      q167 = q166*x%d1val3
      q168 = q139*x%d1val1
      q169 = 6.0_dp*x%d1val2
      q170 = q168*q169
      q171 = q15*x%d1val1
      q172 = q15*q160*x%d1val3
      q173 = q80*q81
      q174 = q75*x%d2val2_d1val3
      q175 = q1*x%d1val2_d1val3
      q176 = q175*x%d2val2
      q177 = 12.0_dp*q176
      q178 = q31*x%d1val2
      q179 = x%d1val3*x%d2val2
      q180 = q178*q179
      q181 = q83*x%d1val2_d1val3
      q182 = q63*q77
      q183 = q63*q78
      q184 = q15*q182
      q185 = q55*x%d1val2
      q186 = q185*x%d1val2_d1val3
      q187 = q186*y - q186 - q84*y + q84 + x%d2val2_d1val3
      q188 = q187*q69
      q189 = q19*q55
      q190 = q189*y - q189 + q53
      q191 = q130*q91
      q192 = q185*x%d1val2_d2val3
      q193 = pow2(x%d1val2_d1val3)
      q194 = q193*q55
      q195 = q83*x%d2val3
      q196 = q136*x%d1val2
      q197 = q196*x%d1val2_d1val3
      q198 = q140*q19
      q199 = q192*y - q192 + q194*y - q194 - q195*y + q195 - q197*y + q197 + q198*y - q198 + x%d2val2_d2val3
      q200 = 3.0_dp*q2
      q201 = 6.0_dp*x%d2val1_d1val3
      q202 = q201*q94
      q203 = q1*x%d1val1_d2val3
      q204 = 3.0_dp*x%d2val1
      q205 = q203*q204
      q206 = q18*q201
      q207 = q12*q161
      q208 = q204*q207
      q209 = 3.0_dp*q57
      q210 = 6.0_dp*x%d2val1
      q211 = q17*x%d1val1_d1val3
      q212 = q210*q211
      q213 = q40*x%d1val1_d2val3
      q214 = 4.0_dp*q100
      q215 = q168*q210
      q216 = q121*x%d1val1_d1val3
      q217 = 9.0_dp*y
      q218 = 6.0_dp*q103
      q219 = q24*powm1(pow4(x%val))
      q220 = 12.0_dp*q219
      q221 = q216*x%d1val3
      q222 = 3.0_dp*q15
      q223 = 18.0_dp*q219
      q224 = 2.0_dp*q107
      q225 = q15*q29
      q226 = 6.0_dp*q219
      q227 = 2.0_dp*x%d3val1 + q152*q225 - q166*q32 - q2*q210 + q210*q57 + q34
      q228 = 3.0_dp*q95
      q229 = q18*q204
      q230 = q69*(6.0_dp*q99 - q200*x%d2val1_d1val3 + q209*x%d2val1_d1val3 - q214*x%d1val3 - q217*q99 + q218*x%d1val3 + q222*q99 - q224*x%d1val3 + q228*y - q228 - q229*y + q229 + x%d3val1_d1val3)
      q231 = q56*x%d1val1_d1val2_d2val3
      q232 = q55*x%d1val1_d1val2
      q233 = q232*x%d1val1_d2val3
      q234 = q7*x%d1val1_d1val2_d1val3
      q235 = q234*x%d1val1_d1val3
      q236 = q40*x%d1val2_d2val3
      q237 = q152*q161*x%d1val1_d1val2
      q238 = q136*x%d1val1_d1val2_d1val3
      q239 = q33*x%d1val1_d1val3
      q240 = q239*q59
      q241 = q160*x%d1val1_d2val3
      q242 = q115*q239
      q243 = q133*q160
      q244 = q106*q117
      q245 = 4.0_dp*q139
      q246 = q245*x%d1val1_d1val2
      q247 = 2.0_dp*q122
      q248 = q247*x%d2val3
      q249 = x%d1val2_d1val3*x%d1val3
      q250 = 4.0_dp*q121*q249
      q251 = q169*q219*q6
      q252 = q131*x%d1val2_d2val3
      q253 = q128*q163
      q254 = q128*x%d1val2_d1val3
      q255 = 1.5_dp*y
      q256 = 3.0_dp*q128
      q257 = q139*x%d1val2
      q258 = q128*q257
      q259 = 0.5_dp*q15
      q260 = q125*q55*x%d1val2_d1val3
      q261 = q166*x%d1val2
      q262 = q261*x%d1val3
      q263 = q7*x%d1val1_d1val2
      q264 = q263*x%d1val1
      q265 = 2.0_dp*q41
      q266 = 2.0_dp*x%d2val1_d1val2 + q264*y - q264 - q265*y + q265
      q267 = q130*q266*x%d2val3
      q268 = q266*q71
      q269 = q142*q4
      q270 = q56*x%d1val1_d1val2_d1val3
      q271 = q150*x%d1val1_d1val2
      q272 = q115*q152
      q273 = q137*q160
      q274 = q247*x%d1val3
      q275 = q69*(-q119*y + q119 - q164*q272 + q270*y - q270 + q271*y - q271 + q272*x%d1val1 - q273*y + q273 + q274*y - q274 + x%d2val1_d1val2_d1val3)
      q276 = q232*x%d1val2_d2val3
      q277 = q234*x%d1val2_d1val3
      q278 = q185*x%d1val1_d1val2_d2val3
      q279 = q160*x%d1val1_d1val2*x%d2val3
      q280 = q115*q33*x%d1val2_d1val3
      q281 = q196*x%d1val1_d1val2_d1val3
      q282 = q164*q33
      q283 = q246*x%d1val2
      q284 = 2.0_dp*q83*y
      q285 = q165*x%d1val2
      q286 = 4.0_dp*q52
      q287 = 0.5_dp*q54
      q288 = q7*x%d1val2
      q289 = 2.0_dp*q84 - q288*x%d1val2_d1val3 + x%d2val2_d1val3
      q290 = q190*y
      q291 = 0.5_dp*q290
      q292 = q263*x%d1val2
      q293 = -4.0_dp*q164*q83 + 2.0_dp*x%d1val1_d2val2 + q190*q57 - q2*q54 + q292*y - q292
      q294 = q232*x%d1val2_d1val3
      q295 = q185*x%d1val1_d1val2_d1val3
      q296 = q115*q160
      q297 = q69*(q18*q287 - q18*q291 + q187*q57 - q2*q289 - q284*x%d1val1_d1val3 - q285*q33 + q286*q65 - q287*q94 + q291*q94 + q294*y - q294 + q295*y - q295 - q296*y + q296 + x%d1val1_d2val2_d1val3)
      q298 = 3.0_dp*q4
      q299 = q298*x%d2val2_d2val3
      q300 = 6.0_dp*q175*x%d2val2_d1val3
      q301 = 3.0_dp*x%d2val2
      q302 = q1*q301*x%d1val2_d2val3
      q303 = q262*x%d2val2_d1val3
      q304 = q163*q301
      q305 = q166*q179*x%d1val2_d1val3
      q306 = q83*x%d1val2_d2val3
      q307 = 4.0_dp*q182
      q308 = 6.0_dp*x%d2val2
      q309 = q257*q308
      q310 = q19*q63
      q311 = 6.0_dp*q183
      q312 = 2.0_dp*q184
      q313 = q15*q77
      q314 = q308*q4
      q315 = 2.0_dp*x%d3val2 + q152*q313 - q166*q78 + q314*y - q314 + q79
      q316 = q298*x%d2val2_d1val3
      q317 = 3.0_dp*q176
      q318 = q23*q301
      q319 = q69*(6.0_dp*q181 - q181*q217 + q181*q222 - q307*x%d1val3 + q311*x%d1val3 - q312*x%d1val3 + q316*y - q316 + q317*y - q317 - q318*y + q318 + x%d3val2_d1val3)
      unary%val = q0
      unary%d1val1 = q2*q3
      unary%d1val2 = q3*q4
      unary%d1val3 = q5*x%d1val3
      unary%d2val1 = q10*q9
      unary%d1val1_d1val2 = q11*q3 + q14*q16 - q14*q3
      unary%d1val1_d1val3 = q16*q18 - q18*q3 + q5*x%d1val1_d1val3
      unary%d2val2 = q10*q22
      unary%d1val2_d1val3 = q16*q23 - q23*q3 + q5*x%d1val2_d1val3
      unary%d2val3 = q26*q5
      unary%d3val1 = q10*q35
      unary%d2val1_d1val2 = q10*q43 + q16*q37 - q3*q37
      unary%d2val1_d1val3 = q10*q48 + q16*q44 - q3*q44
      unary%d1val1_d2val2 = q10*q58
      unary%d1val1_d1val2_d1val3 = -3.0_dp*q16*q66 + 2.0_dp*q3*q66 + q16*q60 + q16*q61 + q16*q62 - q3*q60 - q3*q61 - q3*q62 + q5*x%d1val1_d1val2_d1val3 + q65*q68
      unary%d1val1_d2val3 = q5*(q26*q57 + q70*y - q70 - q73*x%d1val1 - q74*x%d1val1 + x%d1val1_d2val3)
      unary%d3val2 = q10*q80
      unary%d2val2_d1val3 = q10*q89 + q16*q82 - q3*q82
      unary%d1val2_d2val3 = q5*(q4*q91 - q73*x%d1val2 - q74*x%d1val2 + q90*y - q90 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q10*(24.0_dp*q99 + 4.0_dp*x%d3val1_d1val3 - q100*q101 - q102*q99 + q103*q104 + q105*q99 - q106*q107 + q93*y - q93 + q96*y - q96 - q98*y + q98) + q16*q92 - q3*q92
      unary%d2val1_d1val2_d1val3 = -0.75_dp*q110*q16 + 0.5_dp*q110*q3 + q10*(4.0_dp*x%d2val1_d1val2_d1val3 + q113*y - q113 + q114*y - q114 - q116*y + q116 - q118*y + q118 - q120*y + q120 + q123*y - q123) + q108*q16 - q108*q3 + q109*q36*q68 + q111*q16 - q111*q3 + q112*q16 - q112*q3
      unary%d2val1_d2val3 = q5*(q126*y - q126 - q128*q129 - q128*q72 + q131*q91 + q142)
      unary%d1val1_d2val2_d1val3 = q10*(-16.0_dp*q14*q148 + 4.0_dp*x%d1val1_d2val2_d1val3 + q101*q52*q64 + q144*y - q144 + q145*y - q145 - q147*y + q147 - q149*q30 - q150*q54 + q151*q153 - q154*q18 + q154*q94 - q56*(2.0_dp*x%d2val2_d1val3 + q88) + q57*q89) + q143*q16 - q143*q3
      unary%d1val1_d1val2_d2val3 = q5*(-3.0_dp*q163*q164 - 3.0_dp*q72*x%d1val1_d1val2 + 11.0_dp*q168*x%d1val2*y + 2.0_dp*q162 + q117*q136 - q117*q167*y + q136*q59 + q15*q153*q59 + q15*q162 - q15*q170 + q155*y - q155 + q156*y - q156 + q157*y - q157 + q158*y - q158 + q159*y - q159 + q160*q161 + q163*q171 - q165*q167 + q168*q67 - q170 + q172*x%d1val1_d1val3 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q10*(24.0_dp*q181 + 4.0_dp*x%d3val2_d1val3 - q101*q182 - q102*q181 + q104*q183 + q105*q181 - q106*q184 + q174*y - q174 + q177*y - q177 - q180*y + q180) + q16*q173 - q173*q3
      unary%d2val2_d2val3 = q5*(-q129*q190 + q188*y - q188 + q190*q191 - q190*q72 + q199)
      unary%d3val1_d2val3 = q5*(-18.0_dp*q12*q133*q164 + 6.0_dp*q213 + q102*q221 - q104*q216 - q105*q221 - q129*q227 + q133*q166*q171 + q133*q97 + q191*q227 - q200*x%d2val1_d2val3 + q202*y - q202 + q205*y - q205 - q206*y + q206 - q208*y + q208 + q209*x%d2val1_d2val3 - q212*y + q212 - q213*q217 + q213*q222 - q214*x%d2val3 + q215*y - q215 + q218*x%d2val3 + q220*q29 - q223*q32 - q224*x%d2val3 + q225*q226 - q227*q72 + q230*y - q230 + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q5*(0.5_dp*q128*q139*q67 + 5.5_dp*q258*y + q125*q172 + q125*q196 - q125*q262*y - q148*q17*q256 + q15*q17*q254 + q153*q254 - q164*q238 - q164*q241 + q164*q244*q63 + q164*q246 - q222*q258 + q231*y - q231 + q233*y - q233 + q235*y - q235 - q236*y + q236 - q237*y + q237 + q238*x%d1val1 - q240*y + q240 + q241*x%d1val1 - q242*y + q242 - q243*y + q243 - q244*q64 - q246*x%d1val1 + q248*y - q248 + q250*y - q250 - q251*y + q251 + q252*y - q252 - q253*q255 + q253*q259 + q253 - q255*q268 - q256*q257 + q259*q268 + q260*y - q260 + q267*y - q267 + q268 + q269*y - q269 + q275*y - q275 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q5*(q101*q285*q63 + q106*q149*q63 - q117*q148*q30 - q129*q293 - q139*q151 + q150*q187*y - q150*q289 - q153*q164*q187 + q153*q289*x%d1val1 + q161*q286*q63 + q168*q290 + q191*q293 - q193*q282 + q199*q57 + q2*(-2.0_dp*q195 - q146*q249 + q19*q245 + q193*q7 + q288*x%d1val2_d2val3 - x%d2val2_d2val3) - q203*q287 + q203*q291 + q207*q287 - q207*q291 - q211*q290 + q211*q54 - q220*q52*x%d1val1 + q276*y - q276 + q277*y - q277 + q278*y - q278 - q279*y + q279 - q280*y + q280 - q281*y + q281 - q282*x%d1val2*x%d1val2_d2val3 + q283*y - q283 - q284*x%d1val1_d2val3 - q293*q72 + q297*y - q297 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q5*(-18.0_dp*q13*q193*y + 36.0_dp*q249*q52*q63 + 6.0_dp*q306 - q104*q310*x%d1val2_d1val3 - q105*q249*q310 - q129*q315 + q15*q193*q261 + q178*q193 + q191*q315 - q217*q306 + q220*q77 + q222*q306 - q223*q78 + q226*q313 + q299*y - q299 + q300*y - q300 + q302*y - q302 - q303*y + q303 - q304*y + q304 - q305*y + q305 - q307*x%d2val3 + q309*y - q309 + q311*x%d2val3 - q312*x%d2val3 - q315*q72 + q319*y - q319 + x%d3val2_d2val3)
   end function pow_self_real
   
   function pow_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = pow(z, x%val)
      q1 = log(z)
      q2 = q0*q1
      q3 = pow2(x%d1val1)
      q4 = q1*q3
      q5 = q4 + x%d2val1
      q6 = pow2(q1)
      q7 = q0*q6
      q8 = x%d1val1*x%d1val2
      q9 = q7*x%d1val3
      q10 = pow2(x%d1val2)
      q11 = q1*q10
      q12 = q11 + x%d2val2
      q13 = pow2(x%d1val3)
      q14 = q1*q13 + x%d2val3
      q15 = q1*x%d1val1
      q16 = 3.0_dp*q15
      q17 = q6*pow3(x%d1val1)
      q18 = q16*x%d2val1 + q17 + x%d3val1
      q19 = 2.0_dp*q1
      q20 = q19*x%d1val1
      q21 = q20*x%d1val1_d1val2 + x%d2val1_d1val2
      q22 = q7*x%d1val2
      q23 = q20*x%d1val1_d1val3 + x%d2val1_d1val3
      q24 = q19*x%d1val2
      q25 = 16.0_dp*q11 + 16.0_dp*x%d2val2
      q26 = 0.0625_dp*q15
      q27 = q24*x%d1val1_d1val2 + q25*q26 + x%d1val1_d2val2
      q28 = q7*x%d1val2_d1val3
      q29 = q6*x%d1val1_d1val2
      q30 = q0*x%d1val3
      q31 = pow3(q1)
      q32 = q31*q8
      q33 = q19*x%d1val3
      q34 = q1*q14
      q35 = q1*x%d1val2
      q36 = 3.0_dp*q35
      q37 = q6*pow3(x%d1val2)
      q38 = q36*x%d2val2 + q37 + x%d3val2
      q39 = q24*x%d1val2_d1val3 + x%d2val2_d1val3
      q40 = q1*x%d1val1_d1val3
      q41 = q40*x%d2val1
      q42 = 3.0_dp*q6
      q43 = q3*x%d1val1_d1val3
      q44 = q31*x%d1val2
      q45 = q19*x%d1val1_d1val2
      q46 = q20*x%d1val1_d1val2_d1val3 + q45*x%d1val1_d1val3 + x%d2val1_d1val2_d1val3
      q47 = 4.0_dp*q4 + 4.0_dp*x%d2val1
      q48 = 0.25_dp*q34
      q49 = pow2(x%d1val1_d1val3)
      q50 = q19*q49 + q20*x%d1val1_d2val3 + x%d2val1_d2val3
      q51 = q1*x%d2val3
      q52 = q19*x%d1val2_d1val3
      q53 = q6*x%d2val3
      q54 = q6*x%d1val2_d1val3
      q55 = 2.0_dp*x%d1val3
      q56 = q55*q6*x%d1val2
      q57 = 3.0_dp*q1
      q58 = x%d1val2_d1val3*x%d2val2
      q59 = q10*x%d1val2_d1val3
      q60 = 4.0_dp*q11 + 4.0_dp*x%d2val2
      q61 = pow2(x%d1val2_d1val3)
      q62 = q19*q61 + q24*x%d1val2_d2val3 + x%d2val2_d2val3
      q63 = 6.0_dp*q6
      q64 = 12.0_dp*q6
      q65 = 0.5_dp*x%d1val3
      q66 = q1*q65
      q67 = 0.125_dp*q34
      q68 = 4.0_dp*x%d1val1_d1val2_d1val3
      q69 = 4.0_dp*q15
      q70 = 0.5_dp*q69*x%d1val1_d1val2 + x%d2val1_d1val2
      q71 = 0.25_dp*q1
      q72 = 0.25_dp*q47
      q73 = q1*x%d1val2_d1val3
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q2*q5
      unary%d1val1_d1val2 = q2*x%d1val1_d1val2 + q7*q8
      unary%d1val1_d1val3 = q2*x%d1val1_d1val3 + q9*x%d1val1
      unary%d2val2 = q12*q2
      unary%d1val2_d1val3 = q2*x%d1val2_d1val3 + q9*x%d1val2
      unary%d2val3 = q14*q2
      unary%d3val1 = q18*q2
      unary%d2val1_d1val2 = q2*q21 + q22*q5
      unary%d2val1_d1val3 = q2*q23 + q5*q9
      unary%d1val1_d2val2 = q2*q27
      unary%d1val1_d1val2_d1val3 = q2*x%d1val1_d1val2_d1val3 + q22*x%d1val1_d1val3 + q28*x%d1val1 + q29*q30 + q30*q32
      unary%d1val1_d2val3 = q2*(q33*x%d1val1_d1val3 + q34*x%d1val1 + x%d1val1_d2val3)
      unary%d3val2 = q2*q38
      unary%d2val2_d1val3 = q12*q9 + q2*q39
      unary%d1val2_d2val3 = q2*(q33*x%d1val2_d1val3 + q34*x%d1val2 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q18*q9 + q2*(3.0_dp*q41 + q16*x%d2val1_d1val3 + q42*q43 + x%d3val1_d1val3)
      unary%d2val1_d1val2_d1val3 = q2*q46 + q21*q9 + q22*q23 + q28*q5 + q30*q44*q5
      unary%d2val1_d2val3 = q2*(q23*q33 + q47*q48 + q50)
      unary%d1val1_d2val2_d1val3 = q2*(0.0625_dp*q25*q40 + q24*x%d1val1_d1val2_d1val3 + q26*(16.0_dp*x%d2val2_d1val3 + 32.0_dp*q35*x%d1val2_d1val3) + q45*x%d1val2_d1val3 + x%d1val1_d2val2_d1val3) + q27*q9
      unary%d1val1_d1val2_d2val3 = q2*(q13*q29 + q13*q32 + q15*x%d1val2_d2val3 + q33*x%d1val1_d1val2_d1val3 + q35*x%d1val1_d2val3 + q51*x%d1val1_d1val2 + q52*x%d1val1_d1val3 + q53*q8 + q54*q55*x%d1val1 + q56*x%d1val1_d1val3 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q2*(q36*x%d2val2_d1val3 + q42*q59 + q57*q58 + x%d3val2_d1val3) + q38*q9
      unary%d2val2_d2val3 = q2*(q33*q39 + q48*q60 + q62)
      unary%d3val1_d2val3 = q2*(6.0_dp*q40*x%d2val1_d1val3 + q16*x%d2val1_d2val3 + q3*q42*x%d1val1_d2val3 + q49*q63*x%d1val1 + q57*x%d1val1_d2val3*x%d2val1 + q66*(12.0_dp*q15*x%d2val1_d1val3 + 12.0_dp*q41 + 4.0_dp*x%d3val1_d1val3 + q43*q64) + q67*(24.0_dp*q15*x%d2val1 + 8.0_dp*q17 + 8.0_dp*x%d3val1) + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(q13*q44*q72 + q13*q6*q70 + q20*x%d1val1_d1val2_d2val3 + q23*q52 + q23*q56 + q33*q46 + q35*q50 + q40*q68 + q45*x%d1val1_d2val3 + q47*q54*q65 + q47*q71*x%d1val2_d2val3 + q51*q70 + q53*q72*x%d1val2 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(q15*q62 + q19*q39*x%d1val1_d1val3 + q24*x%d1val1_d1val2_d2val3 + q45*x%d1val2_d2val3 + q60*q71*x%d1val1_d2val3 + q66*(4.0_dp*x%d1val1_d2val2_d1val3 + 8.0_dp*q35*x%d1val1_d1val2_d1val3 + 8.0_dp*q73*x%d1val1_d1val2 + q39*q69 + q40*q60) + q67*(16.0_dp*q35*x%d1val1_d1val2 + 8.0_dp*x%d1val1_d2val2 + q20*q60) + q68*q73 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(6.0_dp*q73*x%d2val2_d1val3 + q10*q42*x%d1val2_d2val3 + q36*x%d2val2_d2val3 + q57*x%d1val2_d2val3*x%d2val2 + q61*q63*x%d1val2 + q66*(12.0_dp*q1*q58 + 12.0_dp*q35*x%d2val2_d1val3 + 4.0_dp*x%d3val2_d1val3 + q59*q64) + q67*(24.0_dp*q35*x%d2val2 + 8.0_dp*q37 + 8.0_dp*x%d3val2) + x%d3val2_d2val3)
   end function pow_real_self
   
   function pow_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q319
      real(dp) :: q318
      real(dp) :: q317
      real(dp) :: q316
      real(dp) :: q315
      real(dp) :: q314
      real(dp) :: q313
      real(dp) :: q312
      real(dp) :: q311
      real(dp) :: q310
      real(dp) :: q309
      real(dp) :: q308
      real(dp) :: q307
      real(dp) :: q306
      real(dp) :: q305
      real(dp) :: q304
      real(dp) :: q303
      real(dp) :: q302
      real(dp) :: q301
      real(dp) :: q300
      real(dp) :: q299
      real(dp) :: q298
      real(dp) :: q297
      real(dp) :: q296
      real(dp) :: q295
      real(dp) :: q294
      real(dp) :: q293
      real(dp) :: q292
      real(dp) :: q291
      real(dp) :: q290
      real(dp) :: q289
      real(dp) :: q288
      real(dp) :: q287
      real(dp) :: q286
      real(dp) :: q285
      real(dp) :: q284
      real(dp) :: q283
      real(dp) :: q282
      real(dp) :: q281
      real(dp) :: q280
      real(dp) :: q279
      real(dp) :: q278
      real(dp) :: q277
      real(dp) :: q276
      real(dp) :: q275
      real(dp) :: q274
      real(dp) :: q273
      real(dp) :: q272
      real(dp) :: q271
      real(dp) :: q270
      real(dp) :: q269
      real(dp) :: q268
      real(dp) :: q267
      real(dp) :: q266
      real(dp) :: q265
      real(dp) :: q264
      real(dp) :: q263
      real(dp) :: q262
      real(dp) :: q261
      real(dp) :: q260
      real(dp) :: q259
      real(dp) :: q258
      real(dp) :: q257
      real(dp) :: q256
      real(dp) :: q255
      real(dp) :: q254
      real(dp) :: q253
      real(dp) :: q252
      real(dp) :: q251
      real(dp) :: q250
      real(dp) :: q249
      real(dp) :: q248
      real(dp) :: q247
      real(dp) :: q246
      real(dp) :: q245
      real(dp) :: q244
      real(dp) :: q243
      real(dp) :: q242
      real(dp) :: q241
      real(dp) :: q240
      real(dp) :: q239
      real(dp) :: q238
      real(dp) :: q237
      real(dp) :: q236
      real(dp) :: q235
      real(dp) :: q234
      real(dp) :: q233
      real(dp) :: q232
      real(dp) :: q231
      real(dp) :: q230
      real(dp) :: q229
      real(dp) :: q228
      real(dp) :: q227
      real(dp) :: q226
      real(dp) :: q225
      real(dp) :: q224
      real(dp) :: q223
      real(dp) :: q222
      real(dp) :: q221
      real(dp) :: q220
      real(dp) :: q219
      real(dp) :: q218
      real(dp) :: q217
      real(dp) :: q216
      real(dp) :: q215
      real(dp) :: q214
      real(dp) :: q213
      real(dp) :: q212
      real(dp) :: q211
      real(dp) :: q210
      real(dp) :: q209
      real(dp) :: q208
      real(dp) :: q207
      real(dp) :: q206
      real(dp) :: q205
      real(dp) :: q204
      real(dp) :: q203
      real(dp) :: q202
      real(dp) :: q201
      real(dp) :: q200
      real(dp) :: q199
      real(dp) :: q198
      real(dp) :: q197
      real(dp) :: q196
      real(dp) :: q195
      real(dp) :: q194
      real(dp) :: q193
      real(dp) :: q192
      real(dp) :: q191
      real(dp) :: q190
      real(dp) :: q189
      real(dp) :: q188
      real(dp) :: q187
      real(dp) :: q186
      real(dp) :: q185
      real(dp) :: q184
      real(dp) :: q183
      real(dp) :: q182
      real(dp) :: q181
      real(dp) :: q180
      real(dp) :: q179
      real(dp) :: q178
      real(dp) :: q177
      real(dp) :: q176
      real(dp) :: q175
      real(dp) :: q174
      real(dp) :: q173
      real(dp) :: q172
      real(dp) :: q171
      real(dp) :: q170
      real(dp) :: q169
      real(dp) :: q168
      real(dp) :: q167
      real(dp) :: q166
      real(dp) :: q165
      real(dp) :: q164
      real(dp) :: q163
      real(dp) :: q162
      real(dp) :: q161
      real(dp) :: q160
      real(dp) :: q159
      real(dp) :: q158
      real(dp) :: q157
      real(dp) :: q156
      real(dp) :: q155
      real(dp) :: q154
      real(dp) :: q153
      real(dp) :: q152
      real(dp) :: q151
      real(dp) :: q150
      real(dp) :: q149
      real(dp) :: q148
      real(dp) :: q147
      real(dp) :: q146
      real(dp) :: q145
      real(dp) :: q144
      real(dp) :: q143
      real(dp) :: q142
      real(dp) :: q141
      real(dp) :: q140
      real(dp) :: q139
      real(dp) :: q138
      real(dp) :: q137
      real(dp) :: q136
      real(dp) :: q135
      real(dp) :: q134
      real(dp) :: q133
      real(dp) :: q132
      real(dp) :: q131
      real(dp) :: q130
      real(dp) :: q129
      real(dp) :: q128
      real(dp) :: q127
      real(dp) :: q126
      real(dp) :: q125
      real(dp) :: q124
      real(dp) :: q123
      real(dp) :: q122
      real(dp) :: q121
      real(dp) :: q120
      real(dp) :: q119
      real(dp) :: q118
      real(dp) :: q117
      real(dp) :: q116
      real(dp) :: q115
      real(dp) :: q114
      real(dp) :: q113
      real(dp) :: q112
      real(dp) :: q111
      real(dp) :: q110
      real(dp) :: q109
      real(dp) :: q108
      real(dp) :: q107
      real(dp) :: q106
      real(dp) :: q105
      real(dp) :: q104
      real(dp) :: q103
      real(dp) :: q102
      real(dp) :: q101
      real(dp) :: q100
      real(dp) :: q99
      real(dp) :: q98
      real(dp) :: q97
      real(dp) :: q96
      real(dp) :: q95
      real(dp) :: q94
      real(dp) :: q93
      real(dp) :: q92
      real(dp) :: q91
      real(dp) :: q90
      real(dp) :: q89
      real(dp) :: q88
      real(dp) :: q87
      real(dp) :: q86
      real(dp) :: q85
      real(dp) :: q84
      real(dp) :: q83
      real(dp) :: q82
      real(dp) :: q81
      real(dp) :: q80
      real(dp) :: q79
      real(dp) :: q78
      real(dp) :: q77
      real(dp) :: q76
      real(dp) :: q75
      real(dp) :: q74
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = pow(x%val, y_dp)
      q1 = powm1(x%val)
      q2 = q1*x%d1val1
      q3 = q0*y_dp
      q4 = q1*x%d1val2
      q5 = q1*q3
      q6 = pow2(x%d1val1)
      q7 = 4.0_dp*q1
      q8 = q6*q7
      q9 = 4.0_dp*x%d2val1 + q8*y_dp - q8
      q10 = 0.25_dp*q5
      q11 = q1*x%d1val1_d1val2
      q12 = powm1(pow2(x%val))
      q13 = q12*x%d1val2
      q14 = q13*x%d1val1
      q15 = pow2(y_dp)
      q16 = q0*q15
      q17 = q12*x%d1val3
      q18 = q17*x%d1val1
      q19 = pow2(x%d1val2)
      q20 = q19*q7
      q21 = -q20
      q22 = 4.0_dp*x%d2val2 + q20*y_dp + q21
      q23 = q13*x%d1val3
      q24 = pow2(x%d1val3)
      q25 = q1*q24
      q26 = q25*y_dp - q25 + x%d2val3
      q27 = 12.0_dp*q2
      q28 = q27*x%d2val1
      q29 = pow3(x%d1val1)
      q30 = 8.0_dp*q12
      q31 = 12.0_dp*q12
      q32 = q29*y_dp
      q33 = 4.0_dp*q12
      q34 = q29*q33
      q35 = 4.0_dp*x%d3val1 + q15*q34 + q28*y_dp - q28 + q29*q30 - q31*q32
      q36 = 0.25_dp*q9
      q37 = q13*q36
      q38 = 8.0_dp*q2
      q39 = q38*x%d1val1_d1val2
      q40 = q12*q6
      q41 = q40*x%d1val2
      q42 = 4.0_dp*q41
      q43 = 4.0_dp*x%d2val1_d1val2 + q39*y_dp - q39 - q42*y_dp + q42
      q44 = q17*q36
      q45 = q38*x%d1val1_d1val3
      q46 = q40*x%d1val3
      q47 = 4.0_dp*q46
      q48 = 4.0_dp*x%d2val1_d1val3 + q45*y_dp - q45 - q47*y_dp + q47
      q49 = 8.0_dp*q11
      q50 = q49*x%d1val2
      q51 = q30*x%d1val1
      q52 = q19*y_dp
      q53 = 2.0_dp*x%d2val2
      q54 = q21 + q53
      q55 = 2.0_dp*q1
      q56 = q55*x%d1val1
      q57 = q2*y_dp
      q58 = 4.0_dp*x%d1val1_d2val2 + q22*q57 + q50*y_dp - q50 - q51*q52 - q54*q56
      q59 = x%d1val1*x%d1val2_d1val3
      q60 = q12*q59
      q61 = q17*x%d1val1_d1val2
      q62 = q13*x%d1val1_d1val3
      q63 = powm1(pow3(x%val))
      q64 = q63*x%d1val1
      q65 = q64*x%d1val3
      q66 = q65*x%d1val2
      q67 = x%d1val2*pow3(y_dp)
      q68 = q0*q67
      q69 = q55*x%d1val3
      q70 = q69*x%d1val1_d1val3
      q71 = q12*q24
      q72 = q71*y_dp
      q73 = 2.0_dp*q72
      q74 = q1*(-2.0_dp*q25 + x%d2val3)
      q75 = 12.0_dp*q4
      q76 = q75*x%d2val2
      q77 = pow3(x%d1val2)
      q78 = q77*y_dp
      q79 = q33*q77
      q80 = 4.0_dp*x%d3val2 + q15*q79 + q30*q77 - q31*q78 + q76*y_dp - q76
      q81 = 0.25_dp*q17
      q82 = q22*q81
      q83 = q12*q19
      q84 = q83*x%d1val3
      q85 = 4.0_dp*q84
      q86 = 8.0_dp*q4
      q87 = q86*x%d1val2_d1val3
      q88 = q85 - q87
      q89 = 4.0_dp*x%d2val2_d1val3 - q85*y_dp + q87*y_dp + q88
      q90 = q69*x%d1val2_d1val3
      q91 = q26*y_dp
      q92 = q35*q81
      q93 = q27*x%d2val1_d1val3
      q94 = q1*x%d1val1_d1val3
      q95 = q94*x%d2val1
      q96 = 12.0_dp*q95
      q97 = q31*x%d1val1
      q98 = q97*x%d1val3*x%d2val1
      q99 = q40*x%d1val1_d1val3
      q100 = q29*q63
      q101 = 16.0_dp*x%d1val3
      q102 = 36.0_dp*y_dp
      q103 = q32*q63
      q104 = 24.0_dp*x%d1val3
      q105 = 12.0_dp*q15
      q106 = 8.0_dp*x%d1val3
      q107 = q100*q15
      q108 = q12*q36*x%d1val2_d1val3
      q109 = q63*x%d1val3
      q110 = q109*q9*x%d1val2
      q111 = 0.25_dp*q13*q48
      q112 = q43*q81
      q113 = q38*x%d1val1_d1val2_d1val3
      q114 = q49*x%d1val1_d1val3
      q115 = x%d1val1_d1val2*x%d1val3
      q116 = q115*q51
      q117 = x%d1val1_d1val3*x%d1val2
      q118 = q117*q51
      q119 = q40*x%d1val2_d1val3
      q120 = 4.0_dp*q119
      q121 = q6*q63
      q122 = q121*x%d1val2
      q123 = q106*q122
      q124 = q56*x%d1val1_d1val3
      q125 = q124*y_dp - q124 - q46*y_dp + q46 + x%d2val1_d1val3
      q126 = q125*q69
      q127 = q55*q6
      q128 = 2.0_dp*x%d2val1 + q127*y_dp - q127
      q129 = 0.5_dp*q74
      q130 = 0.5_dp*q1
      q131 = q128*q130
      q132 = q56*x%d1val1_d2val3
      q133 = pow2(x%d1val1_d1val3)
      q134 = q133*q55
      q135 = q40*x%d2val3
      q136 = q33*x%d1val3
      q137 = x%d1val1*x%d1val1_d1val3
      q138 = q136*q137
      q139 = q24*q63
      q140 = 2.0_dp*q139
      q141 = q140*q6
      q142 = q132*y_dp - q132 + q134*y_dp - q134 - q135*y_dp + q135 - q138*y_dp + q138 + q141*y_dp - q141 + x%d2val1_d2val3
      q143 = q58*q81
      q144 = q49*x%d1val2_d1val3
      q145 = q86*x%d1val1_d1val2_d1val3
      q146 = q30*x%d1val2
      q147 = q115*q146
      q148 = x%d1val2_d1val3*y_dp
      q149 = q52*x%d1val1_d1val3
      q150 = q55*x%d1val1_d1val3
      q151 = q54*x%d1val1
      q152 = 2.0_dp*q12
      q153 = q152*x%d1val3
      q154 = q22*y_dp
      q155 = q2*x%d1val2_d2val3
      q156 = q11*x%d2val3
      q157 = q69*x%d1val1_d1val2_d1val3
      q158 = q150*x%d1val2_d1val3
      q159 = q4*x%d1val1_d2val3
      q160 = q152*x%d1val2
      q161 = x%d1val1*x%d2val3
      q162 = q71*x%d1val1_d1val2
      q163 = q13*x%d2val3
      q164 = x%d1val1*y_dp
      q165 = q59*y_dp
      q166 = 6.0_dp*q12
      q167 = q166*x%d1val3
      q168 = q139*x%d1val1
      q169 = 6.0_dp*x%d1val2
      q170 = q168*q169
      q171 = q15*x%d1val1
      q172 = q15*q160*x%d1val3
      q173 = q80*q81
      q174 = q75*x%d2val2_d1val3
      q175 = q1*x%d1val2_d1val3
      q176 = q175*x%d2val2
      q177 = 12.0_dp*q176
      q178 = q31*x%d1val2
      q179 = x%d1val3*x%d2val2
      q180 = q178*q179
      q181 = q83*x%d1val2_d1val3
      q182 = q63*q77
      q183 = q63*q78
      q184 = q15*q182
      q185 = q55*x%d1val2
      q186 = q185*x%d1val2_d1val3
      q187 = q186*y_dp - q186 - q84*y_dp + q84 + x%d2val2_d1val3
      q188 = q187*q69
      q189 = q19*q55
      q190 = q189*y_dp - q189 + q53
      q191 = q130*q91
      q192 = q185*x%d1val2_d2val3
      q193 = pow2(x%d1val2_d1val3)
      q194 = q193*q55
      q195 = q83*x%d2val3
      q196 = q136*x%d1val2
      q197 = q196*x%d1val2_d1val3
      q198 = q140*q19
      q199 = q192*y_dp - q192 + q194*y_dp - q194 - q195*y_dp + q195 - q197*y_dp + q197 + q198*y_dp - q198 + x%d2val2_d2val3
      q200 = 3.0_dp*q2
      q201 = 6.0_dp*x%d2val1_d1val3
      q202 = q201*q94
      q203 = q1*x%d1val1_d2val3
      q204 = 3.0_dp*x%d2val1
      q205 = q203*q204
      q206 = q18*q201
      q207 = q12*q161
      q208 = q204*q207
      q209 = 3.0_dp*q57
      q210 = 6.0_dp*x%d2val1
      q211 = q17*x%d1val1_d1val3
      q212 = q210*q211
      q213 = q40*x%d1val1_d2val3
      q214 = 4.0_dp*q100
      q215 = q168*q210
      q216 = q121*x%d1val1_d1val3
      q217 = 9.0_dp*y_dp
      q218 = 6.0_dp*q103
      q219 = q24*powm1(pow4(x%val))
      q220 = 12.0_dp*q219
      q221 = q216*x%d1val3
      q222 = 3.0_dp*q15
      q223 = 18.0_dp*q219
      q224 = 2.0_dp*q107
      q225 = q15*q29
      q226 = 6.0_dp*q219
      q227 = 2.0_dp*x%d3val1 + q152*q225 - q166*q32 - q2*q210 + q210*q57 + q34
      q228 = 3.0_dp*q95
      q229 = q18*q204
      q230 = q69*(6.0_dp*q99 - q200*x%d2val1_d1val3 + q209*x%d2val1_d1val3 - q214*x%d1val3 - q217*q99 + q218*x%d1val3 + q222*q99 - q224*x%d1val3 + q228*y_dp - q228 - q229*y_dp + q229 + x%d3val1_d1val3)
      q231 = q56*x%d1val1_d1val2_d2val3
      q232 = q55*x%d1val1_d1val2
      q233 = q232*x%d1val1_d2val3
      q234 = q7*x%d1val1_d1val2_d1val3
      q235 = q234*x%d1val1_d1val3
      q236 = q40*x%d1val2_d2val3
      q237 = q152*q161*x%d1val1_d1val2
      q238 = q136*x%d1val1_d1val2_d1val3
      q239 = q33*x%d1val1_d1val3
      q240 = q239*q59
      q241 = q160*x%d1val1_d2val3
      q242 = q115*q239
      q243 = q133*q160
      q244 = q106*q117
      q245 = 4.0_dp*q139
      q246 = q245*x%d1val1_d1val2
      q247 = 2.0_dp*q122
      q248 = q247*x%d2val3
      q249 = x%d1val2_d1val3*x%d1val3
      q250 = 4.0_dp*q121*q249
      q251 = q169*q219*q6
      q252 = q131*x%d1val2_d2val3
      q253 = q128*q163
      q254 = q128*x%d1val2_d1val3
      q255 = 1.5_dp*y_dp
      q256 = 3.0_dp*q128
      q257 = q139*x%d1val2
      q258 = q128*q257
      q259 = 0.5_dp*q15
      q260 = q125*q55*x%d1val2_d1val3
      q261 = q166*x%d1val2
      q262 = q261*x%d1val3
      q263 = q7*x%d1val1_d1val2
      q264 = q263*x%d1val1
      q265 = 2.0_dp*q41
      q266 = 2.0_dp*x%d2val1_d1val2 + q264*y_dp - q264 - q265*y_dp + q265
      q267 = q130*q266*x%d2val3
      q268 = q266*q71
      q269 = q142*q4
      q270 = q56*x%d1val1_d1val2_d1val3
      q271 = q150*x%d1val1_d1val2
      q272 = q115*q152
      q273 = q137*q160
      q274 = q247*x%d1val3
      q275 = q69*(-q119*y_dp + q119 - q164*q272 + q270*y_dp - q270 + q271*y_dp - q271 + q272*x%d1val1 - q273*y_dp + q273 + q274*y_dp - q274 + x%d2val1_d1val2_d1val3)
      q276 = q232*x%d1val2_d2val3
      q277 = q234*x%d1val2_d1val3
      q278 = q185*x%d1val1_d1val2_d2val3
      q279 = q160*x%d1val1_d1val2*x%d2val3
      q280 = q115*q33*x%d1val2_d1val3
      q281 = q196*x%d1val1_d1val2_d1val3
      q282 = q164*q33
      q283 = q246*x%d1val2
      q284 = 2.0_dp*q83*y_dp
      q285 = q165*x%d1val2
      q286 = 4.0_dp*q52
      q287 = 0.5_dp*q54
      q288 = q7*x%d1val2
      q289 = 2.0_dp*q84 - q288*x%d1val2_d1val3 + x%d2val2_d1val3
      q290 = q190*y_dp
      q291 = 0.5_dp*q290
      q292 = q263*x%d1val2
      q293 = -4.0_dp*q164*q83 + 2.0_dp*x%d1val1_d2val2 + q190*q57 - q2*q54 + q292*y_dp - q292
      q294 = q232*x%d1val2_d1val3
      q295 = q185*x%d1val1_d1val2_d1val3
      q296 = q115*q160
      q297 = q69*(q18*q287 - q18*q291 + q187*q57 - q2*q289 - q284*x%d1val1_d1val3 - q285*q33 + q286*q65 - q287*q94 + q291*q94 + q294*y_dp - q294 + q295*y_dp - q295 - q296*y_dp + q296 + x%d1val1_d2val2_d1val3)
      q298 = 3.0_dp*q4
      q299 = q298*x%d2val2_d2val3
      q300 = 6.0_dp*q175*x%d2val2_d1val3
      q301 = 3.0_dp*x%d2val2
      q302 = q1*q301*x%d1val2_d2val3
      q303 = q262*x%d2val2_d1val3
      q304 = q163*q301
      q305 = q166*q179*x%d1val2_d1val3
      q306 = q83*x%d1val2_d2val3
      q307 = 4.0_dp*q182
      q308 = 6.0_dp*x%d2val2
      q309 = q257*q308
      q310 = q19*q63
      q311 = 6.0_dp*q183
      q312 = 2.0_dp*q184
      q313 = q15*q77
      q314 = q308*q4
      q315 = 2.0_dp*x%d3val2 + q152*q313 - q166*q78 + q314*y_dp - q314 + q79
      q316 = q298*x%d2val2_d1val3
      q317 = 3.0_dp*q176
      q318 = q23*q301
      q319 = q69*(6.0_dp*q181 - q181*q217 + q181*q222 - q307*x%d1val3 + q311*x%d1val3 - q312*x%d1val3 + q316*y_dp - q316 + q317*y_dp - q317 - q318*y_dp + q318 + x%d3val2_d1val3)
      unary%val = q0
      unary%d1val1 = q2*q3
      unary%d1val2 = q3*q4
      unary%d1val3 = q5*x%d1val3
      unary%d2val1 = q10*q9
      unary%d1val1_d1val2 = q11*q3 + q14*q16 - q14*q3
      unary%d1val1_d1val3 = q16*q18 - q18*q3 + q5*x%d1val1_d1val3
      unary%d2val2 = q10*q22
      unary%d1val2_d1val3 = q16*q23 - q23*q3 + q5*x%d1val2_d1val3
      unary%d2val3 = q26*q5
      unary%d3val1 = q10*q35
      unary%d2val1_d1val2 = q10*q43 + q16*q37 - q3*q37
      unary%d2val1_d1val3 = q10*q48 + q16*q44 - q3*q44
      unary%d1val1_d2val2 = q10*q58
      unary%d1val1_d1val2_d1val3 = -3.0_dp*q16*q66 + 2.0_dp*q3*q66 + q16*q60 + q16*q61 + q16*q62 - q3*q60 - q3*q61 - q3*q62 + q5*x%d1val1_d1val2_d1val3 + q65*q68
      unary%d1val1_d2val3 = q5*(q26*q57 + q70*y_dp - q70 - q73*x%d1val1 - q74*x%d1val1 + x%d1val1_d2val3)
      unary%d3val2 = q10*q80
      unary%d2val2_d1val3 = q10*q89 + q16*q82 - q3*q82
      unary%d1val2_d2val3 = q5*(q4*q91 - q73*x%d1val2 - q74*x%d1val2 + q90*y_dp - q90 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q10*(24.0_dp*q99 + 4.0_dp*x%d3val1_d1val3 - q100*q101 - q102*q99 + q103*q104 + q105*q99 - q106*q107 + q93*y_dp - q93 + q96*y_dp - q96 - q98*y_dp + q98) + q16*q92 - q3*q92
      unary%d2val1_d1val2_d1val3 = -0.75_dp*q110*q16 + 0.5_dp*q110*q3 + q10*(4.0_dp*x%d2val1_d1val2_d1val3 + q113*y_dp - q113 + q114*y_dp - q114 - q116*y_dp + q116 - q118*y_dp + q118 - q120*y_dp + q120 + q123*y_dp - q123) + q108*q16 - q108*q3 + q109*q36*q68 + q111*q16 - q111*q3 + q112*q16 - q112*q3
      unary%d2val1_d2val3 = q5*(q126*y_dp - q126 - q128*q129 - q128*q72 + q131*q91 + q142)
      unary%d1val1_d2val2_d1val3 = q10*(-16.0_dp*q14*q148 + 4.0_dp*x%d1val1_d2val2_d1val3 + q101*q52*q64 + q144*y_dp - q144 + q145*y_dp - q145 - q147*y_dp + q147 - q149*q30 - q150*q54 + q151*q153 - q154*q18 + q154*q94 - q56*(2.0_dp*x%d2val2_d1val3 + q88) + q57*q89) + q143*q16 - q143*q3
      unary%d1val1_d1val2_d2val3 = q5*(-3.0_dp*q163*q164 - 3.0_dp*q72*x%d1val1_d1val2 + 11.0_dp*q168*x%d1val2*y_dp + 2.0_dp*q162 + q117*q136 - q117*q167*y_dp + q136*q59 + q15*q153*q59 + q15*q162 - q15*q170 + q155*y_dp - q155 + q156*y_dp - q156 + q157*y_dp - q157 + q158*y_dp - q158 + q159*y_dp - q159 + q160*q161 + q163*q171 - q165*q167 + q168*q67 - q170 + q172*x%d1val1_d1val3 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q10*(24.0_dp*q181 + 4.0_dp*x%d3val2_d1val3 - q101*q182 - q102*q181 + q104*q183 + q105*q181 - q106*q184 + q174*y_dp - q174 + q177*y_dp - q177 - q180*y_dp + q180) + q16*q173 - q173*q3
      unary%d2val2_d2val3 = q5*(-q129*q190 + q188*y_dp - q188 + q190*q191 - q190*q72 + q199)
      unary%d3val1_d2val3 = q5*(-18.0_dp*q12*q133*q164 + 6.0_dp*q213 + q102*q221 - q104*q216 - q105*q221 - q129*q227 + q133*q166*q171 + q133*q97 + q191*q227 - q200*x%d2val1_d2val3 + q202*y_dp - q202 + q205*y_dp - q205 - q206*y_dp + q206 - q208*y_dp + q208 + q209*x%d2val1_d2val3 - q212*y_dp + q212 - q213*q217 + q213*q222 - q214*x%d2val3 + q215*y_dp - q215 + q218*x%d2val3 + q220*q29 - q223*q32 - q224*x%d2val3 + q225*q226 - q227*q72 + q230*y_dp - q230 + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q5*(0.5_dp*q128*q139*q67 + 5.5_dp*q258*y_dp + q125*q172 + q125*q196 - q125*q262*y_dp - q148*q17*q256 + q15*q17*q254 + q153*q254 - q164*q238 - q164*q241 + q164*q244*q63 + q164*q246 - q222*q258 + q231*y_dp - q231 + q233*y_dp - q233 + q235*y_dp - q235 - q236*y_dp + q236 - q237*y_dp + q237 + q238*x%d1val1 - q240*y_dp + q240 + q241*x%d1val1 - q242*y_dp + q242 - q243*y_dp + q243 - q244*q64 - q246*x%d1val1 + q248*y_dp - q248 + q250*y_dp - q250 - q251*y_dp + q251 + q252*y_dp - q252 - q253*q255 + q253*q259 + q253 - q255*q268 - q256*q257 + q259*q268 + q260*y_dp - q260 + q267*y_dp - q267 + q268 + q269*y_dp - q269 + q275*y_dp - q275 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q5*(q101*q285*q63 + q106*q149*q63 - q117*q148*q30 - q129*q293 - q139*q151 + q150*q187*y_dp - q150*q289 - q153*q164*q187 + q153*q289*x%d1val1 + q161*q286*q63 + q168*q290 + q191*q293 - q193*q282 + q199*q57 + q2*(-2.0_dp*q195 - q146*q249 + q19*q245 + q193*q7 + q288*x%d1val2_d2val3 - x%d2val2_d2val3) - q203*q287 + q203*q291 + q207*q287 - q207*q291 - q211*q290 + q211*q54 - q220*q52*x%d1val1 + q276*y_dp - q276 + q277*y_dp - q277 + q278*y_dp - q278 - q279*y_dp + q279 - q280*y_dp + q280 - q281*y_dp + q281 - q282*x%d1val2*x%d1val2_d2val3 + q283*y_dp - q283 - q284*x%d1val1_d2val3 - q293*q72 + q297*y_dp - q297 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q5*(-18.0_dp*q13*q193*y_dp + 36.0_dp*q249*q52*q63 + 6.0_dp*q306 - q104*q310*x%d1val2_d1val3 - q105*q249*q310 - q129*q315 + q15*q193*q261 + q178*q193 + q191*q315 - q217*q306 + q220*q77 + q222*q306 - q223*q78 + q226*q313 + q299*y_dp - q299 + q300*y_dp - q300 + q302*y_dp - q302 - q303*y_dp + q303 - q304*y_dp + q304 - q305*y_dp + q305 - q307*x%d2val3 + q309*y_dp - q309 + q311*x%d2val3 - q312*x%d2val3 - q315*q72 + q319*y_dp - q319 + x%d3val2_d2val3)
   end function pow_self_int
   
   function pow_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q73
      real(dp) :: q72
      real(dp) :: q71
      real(dp) :: q70
      real(dp) :: q69
      real(dp) :: q68
      real(dp) :: q67
      real(dp) :: q66
      real(dp) :: q65
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
      real(dp) :: q54
      real(dp) :: q53
      real(dp) :: q52
      real(dp) :: q51
      real(dp) :: q50
      real(dp) :: q49
      real(dp) :: q48
      real(dp) :: q47
      real(dp) :: q46
      real(dp) :: q45
      real(dp) :: q44
      real(dp) :: q43
      real(dp) :: q42
      real(dp) :: q41
      real(dp) :: q40
      real(dp) :: q39
      real(dp) :: q38
      real(dp) :: q37
      real(dp) :: q36
      real(dp) :: q35
      real(dp) :: q34
      real(dp) :: q33
      real(dp) :: q32
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
      q0 = pow(y_dp, x%val)
      q1 = log(y_dp)
      q2 = q0*q1
      q3 = pow2(x%d1val1)
      q4 = q1*q3
      q5 = q4 + x%d2val1
      q6 = pow2(q1)
      q7 = q0*q6
      q8 = x%d1val1*x%d1val2
      q9 = q7*x%d1val3
      q10 = pow2(x%d1val2)
      q11 = q1*q10
      q12 = q11 + x%d2val2
      q13 = pow2(x%d1val3)
      q14 = q1*q13 + x%d2val3
      q15 = q1*x%d1val1
      q16 = 3.0_dp*q15
      q17 = q6*pow3(x%d1val1)
      q18 = q16*x%d2val1 + q17 + x%d3val1
      q19 = 2.0_dp*q1
      q20 = q19*x%d1val1
      q21 = q20*x%d1val1_d1val2 + x%d2val1_d1val2
      q22 = q7*x%d1val2
      q23 = q20*x%d1val1_d1val3 + x%d2val1_d1val3
      q24 = q19*x%d1val2
      q25 = 16.0_dp*q11 + 16.0_dp*x%d2val2
      q26 = 0.0625_dp*q15
      q27 = q24*x%d1val1_d1val2 + q25*q26 + x%d1val1_d2val2
      q28 = q7*x%d1val2_d1val3
      q29 = q6*x%d1val1_d1val2
      q30 = q0*x%d1val3
      q31 = pow3(q1)
      q32 = q31*q8
      q33 = q19*x%d1val3
      q34 = q1*q14
      q35 = q1*x%d1val2
      q36 = 3.0_dp*q35
      q37 = q6*pow3(x%d1val2)
      q38 = q36*x%d2val2 + q37 + x%d3val2
      q39 = q24*x%d1val2_d1val3 + x%d2val2_d1val3
      q40 = q1*x%d1val1_d1val3
      q41 = q40*x%d2val1
      q42 = 3.0_dp*q6
      q43 = q3*x%d1val1_d1val3
      q44 = q31*x%d1val2
      q45 = q19*x%d1val1_d1val2
      q46 = q20*x%d1val1_d1val2_d1val3 + q45*x%d1val1_d1val3 + x%d2val1_d1val2_d1val3
      q47 = 4.0_dp*q4 + 4.0_dp*x%d2val1
      q48 = 0.25_dp*q34
      q49 = pow2(x%d1val1_d1val3)
      q50 = q19*q49 + q20*x%d1val1_d2val3 + x%d2val1_d2val3
      q51 = q1*x%d2val3
      q52 = q19*x%d1val2_d1val3
      q53 = q6*x%d2val3
      q54 = q6*x%d1val2_d1val3
      q55 = 2.0_dp*x%d1val3
      q56 = q55*q6*x%d1val2
      q57 = 3.0_dp*q1
      q58 = x%d1val2_d1val3*x%d2val2
      q59 = q10*x%d1val2_d1val3
      q60 = 4.0_dp*q11 + 4.0_dp*x%d2val2
      q61 = pow2(x%d1val2_d1val3)
      q62 = q19*q61 + q24*x%d1val2_d2val3 + x%d2val2_d2val3
      q63 = 6.0_dp*q6
      q64 = 12.0_dp*q6
      q65 = 0.5_dp*x%d1val3
      q66 = q1*q65
      q67 = 0.125_dp*q34
      q68 = 4.0_dp*x%d1val1_d1val2_d1val3
      q69 = 4.0_dp*q15
      q70 = 0.5_dp*q69*x%d1val1_d1val2 + x%d2val1_d1val2
      q71 = 0.25_dp*q1
      q72 = 0.25_dp*q47
      q73 = q1*x%d1val2_d1val3
      unary%val = q0
      unary%d1val1 = q2*x%d1val1
      unary%d1val2 = q2*x%d1val2
      unary%d1val3 = q2*x%d1val3
      unary%d2val1 = q2*q5
      unary%d1val1_d1val2 = q2*x%d1val1_d1val2 + q7*q8
      unary%d1val1_d1val3 = q2*x%d1val1_d1val3 + q9*x%d1val1
      unary%d2val2 = q12*q2
      unary%d1val2_d1val3 = q2*x%d1val2_d1val3 + q9*x%d1val2
      unary%d2val3 = q14*q2
      unary%d3val1 = q18*q2
      unary%d2val1_d1val2 = q2*q21 + q22*q5
      unary%d2val1_d1val3 = q2*q23 + q5*q9
      unary%d1val1_d2val2 = q2*q27
      unary%d1val1_d1val2_d1val3 = q2*x%d1val1_d1val2_d1val3 + q22*x%d1val1_d1val3 + q28*x%d1val1 + q29*q30 + q30*q32
      unary%d1val1_d2val3 = q2*(q33*x%d1val1_d1val3 + q34*x%d1val1 + x%d1val1_d2val3)
      unary%d3val2 = q2*q38
      unary%d2val2_d1val3 = q12*q9 + q2*q39
      unary%d1val2_d2val3 = q2*(q33*x%d1val2_d1val3 + q34*x%d1val2 + x%d1val2_d2val3)
      unary%d3val1_d1val3 = q18*q9 + q2*(3.0_dp*q41 + q16*x%d2val1_d1val3 + q42*q43 + x%d3val1_d1val3)
      unary%d2val1_d1val2_d1val3 = q2*q46 + q21*q9 + q22*q23 + q28*q5 + q30*q44*q5
      unary%d2val1_d2val3 = q2*(q23*q33 + q47*q48 + q50)
      unary%d1val1_d2val2_d1val3 = q2*(0.0625_dp*q25*q40 + q24*x%d1val1_d1val2_d1val3 + q26*(16.0_dp*x%d2val2_d1val3 + 32.0_dp*q35*x%d1val2_d1val3) + q45*x%d1val2_d1val3 + x%d1val1_d2val2_d1val3) + q27*q9
      unary%d1val1_d1val2_d2val3 = q2*(q13*q29 + q13*q32 + q15*x%d1val2_d2val3 + q33*x%d1val1_d1val2_d1val3 + q35*x%d1val1_d2val3 + q51*x%d1val1_d1val2 + q52*x%d1val1_d1val3 + q53*q8 + q54*q55*x%d1val1 + q56*x%d1val1_d1val3 + x%d1val1_d1val2_d2val3)
      unary%d3val2_d1val3 = q2*(q36*x%d2val2_d1val3 + q42*q59 + q57*q58 + x%d3val2_d1val3) + q38*q9
      unary%d2val2_d2val3 = q2*(q33*q39 + q48*q60 + q62)
      unary%d3val1_d2val3 = q2*(6.0_dp*q40*x%d2val1_d1val3 + q16*x%d2val1_d2val3 + q3*q42*x%d1val1_d2val3 + q49*q63*x%d1val1 + q57*x%d1val1_d2val3*x%d2val1 + q66*(12.0_dp*q15*x%d2val1_d1val3 + 12.0_dp*q41 + 4.0_dp*x%d3val1_d1val3 + q43*q64) + q67*(24.0_dp*q15*x%d2val1 + 8.0_dp*q17 + 8.0_dp*x%d3val1) + x%d3val1_d2val3)
      unary%d2val1_d1val2_d2val3 = q2*(q13*q44*q72 + q13*q6*q70 + q20*x%d1val1_d1val2_d2val3 + q23*q52 + q23*q56 + q33*q46 + q35*q50 + q40*q68 + q45*x%d1val1_d2val3 + q47*q54*q65 + q47*q71*x%d1val2_d2val3 + q51*q70 + q53*q72*x%d1val2 + x%d2val1_d1val2_d2val3)
      unary%d1val1_d2val2_d2val3 = q2*(q15*q62 + q19*q39*x%d1val1_d1val3 + q24*x%d1val1_d1val2_d2val3 + q45*x%d1val2_d2val3 + q60*q71*x%d1val1_d2val3 + q66*(4.0_dp*x%d1val1_d2val2_d1val3 + 8.0_dp*q35*x%d1val1_d1val2_d1val3 + 8.0_dp*q73*x%d1val1_d1val2 + q39*q69 + q40*q60) + q67*(16.0_dp*q35*x%d1val1_d1val2 + 8.0_dp*x%d1val1_d2val2 + q20*q60) + q68*q73 + x%d1val1_d2val2_d2val3)
      unary%d3val2_d2val3 = q2*(6.0_dp*q73*x%d2val2_d1val3 + q10*q42*x%d1val2_d2val3 + q36*x%d2val2_d2val3 + q57*x%d1val2_d2val3*x%d2val2 + q61*q63*x%d1val2 + q66*(12.0_dp*q1*q58 + 12.0_dp*q35*x%d2val2_d1val3 + 4.0_dp*x%d3val2_d1val3 + q59*q64) + q67*(24.0_dp*q35*x%d2val2 + 8.0_dp*q37 + 8.0_dp*x%d3val2) + x%d3val2_d2val3)
   end function pow_int_self
   
   function max_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: binary
      real(dp) :: q1
      real(dp) :: q0
      q0 = Heaviside(x%val - y%val)
      q1 = Heaviside(-x%val + y%val)
      binary%val = Max(x%val, y%val)
      binary%d1val1 = q0*x%d1val1 + q1*y%d1val1
      binary%d1val2 = q0*x%d1val2 + q1*y%d1val2
      binary%d1val3 = q0*x%d1val3 + q1*y%d1val3
      binary%d2val1 = q0*x%d2val1 + q1*y%d2val1
      binary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q1*y%d1val1_d1val2
      binary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q1*y%d1val1_d1val3
      binary%d2val2 = q0*x%d2val2 + q1*y%d2val2
      binary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q1*y%d1val2_d1val3
      binary%d2val3 = q0*x%d2val3 + q1*y%d2val3
      binary%d3val1 = q0*x%d3val1 + q1*y%d3val1
      binary%d2val1_d1val2 = q0*x%d2val1_d1val2 + q1*y%d2val1_d1val2
      binary%d2val1_d1val3 = q0*x%d2val1_d1val3 + q1*y%d2val1_d1val3
      binary%d1val1_d2val2 = q0*x%d1val1_d2val2 + q1*y%d1val1_d2val2
      binary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3 + q1*y%d1val1_d1val2_d1val3
      binary%d1val1_d2val3 = q0*x%d1val1_d2val3 + q1*y%d1val1_d2val3
      binary%d3val2 = q0*x%d3val2 + q1*y%d3val2
      binary%d2val2_d1val3 = q0*x%d2val2_d1val3 + q1*y%d2val2_d1val3
      binary%d1val2_d2val3 = q0*x%d1val2_d2val3 + q1*y%d1val2_d2val3
      binary%d3val1_d1val3 = q0*x%d3val1_d1val3 + q1*y%d3val1_d1val3
      binary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3 + q1*y%d2val1_d1val2_d1val3
      binary%d2val1_d2val3 = q0*x%d2val1_d2val3 + q1*y%d2val1_d2val3
      binary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3 + q1*y%d1val1_d2val2_d1val3
      binary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3 + q1*y%d1val1_d1val2_d2val3
      binary%d3val2_d1val3 = q0*x%d3val2_d1val3 + q1*y%d3val2_d1val3
      binary%d2val2_d2val3 = q0*x%d2val2_d2val3 + q1*y%d2val2_d2val3
      binary%d3val1_d2val3 = q0*x%d3val1_d2val3 + q1*y%d3val1_d2val3
      binary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3 + q1*y%d2val1_d1val2_d2val3
      binary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3 + q1*y%d1val1_d2val2_d2val3
      binary%d3val2_d2val3 = q0*x%d3val2_d2val3 + q1*y%d3val2_d2val3
   end function max_self
   
   function max_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q0
      q0 = Heaviside(x%val - y)
      unary%val = Max(x%val, y)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function max_self_real
   
   function max_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q0
      q0 = Heaviside(x%val - z)
      unary%val = Max(x%val, z)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function max_real_self
   
   function max_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = y
      q0 = Heaviside(x%val - y_dp)
      unary%val = Max(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function max_self_int
   
   function max_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = z
      q0 = Heaviside(x%val - y_dp)
      unary%val = Max(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function max_int_self
   
   function min_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: binary
      real(dp) :: q1
      real(dp) :: q0
      q0 = Heaviside(-x%val + y%val)
      q1 = Heaviside(x%val - y%val)
      binary%val = Min(x%val, y%val)
      binary%d1val1 = q0*x%d1val1 + q1*y%d1val1
      binary%d1val2 = q0*x%d1val2 + q1*y%d1val2
      binary%d1val3 = q0*x%d1val3 + q1*y%d1val3
      binary%d2val1 = q0*x%d2val1 + q1*y%d2val1
      binary%d1val1_d1val2 = q0*x%d1val1_d1val2 + q1*y%d1val1_d1val2
      binary%d1val1_d1val3 = q0*x%d1val1_d1val3 + q1*y%d1val1_d1val3
      binary%d2val2 = q0*x%d2val2 + q1*y%d2val2
      binary%d1val2_d1val3 = q0*x%d1val2_d1val3 + q1*y%d1val2_d1val3
      binary%d2val3 = q0*x%d2val3 + q1*y%d2val3
      binary%d3val1 = q0*x%d3val1 + q1*y%d3val1
      binary%d2val1_d1val2 = q0*x%d2val1_d1val2 + q1*y%d2val1_d1val2
      binary%d2val1_d1val3 = q0*x%d2val1_d1val3 + q1*y%d2val1_d1val3
      binary%d1val1_d2val2 = q0*x%d1val1_d2val2 + q1*y%d1val1_d2val2
      binary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3 + q1*y%d1val1_d1val2_d1val3
      binary%d1val1_d2val3 = q0*x%d1val1_d2val3 + q1*y%d1val1_d2val3
      binary%d3val2 = q0*x%d3val2 + q1*y%d3val2
      binary%d2val2_d1val3 = q0*x%d2val2_d1val3 + q1*y%d2val2_d1val3
      binary%d1val2_d2val3 = q0*x%d1val2_d2val3 + q1*y%d1val2_d2val3
      binary%d3val1_d1val3 = q0*x%d3val1_d1val3 + q1*y%d3val1_d1val3
      binary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3 + q1*y%d2val1_d1val2_d1val3
      binary%d2val1_d2val3 = q0*x%d2val1_d2val3 + q1*y%d2val1_d2val3
      binary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3 + q1*y%d1val1_d2val2_d1val3
      binary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3 + q1*y%d1val1_d1val2_d2val3
      binary%d3val2_d1val3 = q0*x%d3val2_d1val3 + q1*y%d3val2_d1val3
      binary%d2val2_d2val3 = q0*x%d2val2_d2val3 + q1*y%d2val2_d2val3
      binary%d3val1_d2val3 = q0*x%d3val1_d2val3 + q1*y%d3val1_d2val3
      binary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3 + q1*y%d2val1_d1val2_d2val3
      binary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3 + q1*y%d1val1_d2val2_d2val3
      binary%d3val2_d2val3 = q0*x%d3val2_d2val3 + q1*y%d3val2_d2val3
   end function min_self
   
   function min_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q0
      q0 = Heaviside(-x%val + y)
      unary%val = Min(x%val, y)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function min_self_real
   
   function min_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: q0
      q0 = Heaviside(-x%val + z)
      unary%val = Min(x%val, z)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function min_real_self
   
   function min_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = y
      q0 = Heaviside(-x%val + y_dp)
      unary%val = Min(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function min_self_int
   
   function min_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
      real(dp) :: q0
      y_dp = z
      q0 = Heaviside(-x%val + y_dp)
      unary%val = Min(x%val, y_dp)
      unary%d1val1 = q0*x%d1val1
      unary%d1val2 = q0*x%d1val2
      unary%d1val3 = q0*x%d1val3
      unary%d2val1 = q0*x%d2val1
      unary%d1val1_d1val2 = q0*x%d1val1_d1val2
      unary%d1val1_d1val3 = q0*x%d1val1_d1val3
      unary%d2val2 = q0*x%d2val2
      unary%d1val2_d1val3 = q0*x%d1val2_d1val3
      unary%d2val3 = q0*x%d2val3
      unary%d3val1 = q0*x%d3val1
      unary%d2val1_d1val2 = q0*x%d2val1_d1val2
      unary%d2val1_d1val3 = q0*x%d2val1_d1val3
      unary%d1val1_d2val2 = q0*x%d1val1_d2val2
      unary%d1val1_d1val2_d1val3 = q0*x%d1val1_d1val2_d1val3
      unary%d1val1_d2val3 = q0*x%d1val1_d2val3
      unary%d3val2 = q0*x%d3val2
      unary%d2val2_d1val3 = q0*x%d2val2_d1val3
      unary%d1val2_d2val3 = q0*x%d1val2_d2val3
      unary%d3val1_d1val3 = q0*x%d3val1_d1val3
      unary%d2val1_d1val2_d1val3 = q0*x%d2val1_d1val2_d1val3
      unary%d2val1_d2val3 = q0*x%d2val1_d2val3
      unary%d1val1_d2val2_d1val3 = q0*x%d1val1_d2val2_d1val3
      unary%d1val1_d1val2_d2val3 = q0*x%d1val1_d1val2_d2val3
      unary%d3val2_d1val3 = q0*x%d3val2_d1val3
      unary%d2val2_d2val3 = q0*x%d2val2_d2val3
      unary%d3val1_d2val3 = q0*x%d3val1_d2val3
      unary%d2val1_d1val2_d2val3 = q0*x%d2val1_d1val2_d2val3
      unary%d1val1_d2val2_d2val3 = q0*x%d1val1_d2val2_d2val3
      unary%d3val2_d2val3 = q0*x%d3val2_d2val3
   end function min_int_self
   
   function dim_self(x, y) result(binary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: binary
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
      q0 = x%val - y%val
      q1 = sgn(q0)
      q2 = 0.5_dp*q1
      q3 = -0.5_dp*y%d2val1 + 0.5_dp*x%d2val1
      q4 = -0.5_dp*y%d1val1_d1val2 + 0.5_dp*x%d1val1_d1val2
      q5 = -0.5_dp*y%d1val1_d1val3 + 0.5_dp*x%d1val1_d1val3
      q6 = -0.5_dp*y%d2val2 + 0.5_dp*x%d2val2
      q7 = -0.5_dp*y%d1val2_d1val3 + 0.5_dp*x%d1val2_d1val3
      q8 = -0.5_dp*y%d3val1 + 0.5_dp*x%d3val1
      q9 = -0.5_dp*y%d2val1_d1val2 + 0.5_dp*x%d2val1_d1val2
      q10 = -0.5_dp*y%d2val1_d1val3 + 0.5_dp*x%d2val1_d1val3
      q11 = -0.5_dp*y%d1val1_d2val2 + 0.5_dp*x%d1val1_d2val2
      q12 = -0.5_dp*y%d1val1_d1val2_d1val3 + 0.5_dp*x%d1val1_d1val2_d1val3
      q13 = -0.5_dp*y%d3val2 + 0.5_dp*x%d3val2
      q14 = -0.5_dp*y%d2val2_d1val3 + 0.5_dp*x%d2val2_d1val3
      q15 = -0.5_dp*y%d3val1_d1val3 + 0.5_dp*x%d3val1_d1val3
      q16 = -0.5_dp*y%d2val1_d1val2_d1val3 + 0.5_dp*x%d2val1_d1val2_d1val3
      q17 = -0.5_dp*y%d1val1_d2val2_d1val3 + 0.5_dp*x%d1val1_d2val2_d1val3
      q18 = -0.5_dp*y%d3val2_d1val3 + 0.5_dp*x%d3val2_d1val3
      binary%val = -0.5_dp*y%val + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      binary%d1val1 = -0.5_dp*y%d1val1 + 0.5_dp*x%d1val1 + q2*(x%d1val1 - y%d1val1)
      binary%d1val2 = -0.5_dp*y%d1val2 + 0.5_dp*x%d1val2 + q2*(x%d1val2 - y%d1val2)
      binary%d1val3 = -0.5_dp*y%d1val3 + 0.5_dp*x%d1val3 + q2*(x%d1val3 - y%d1val3)
      binary%d2val1 = q1*q3 + q3
      binary%d1val1_d1val2 = q1*q4 + q4
      binary%d1val1_d1val3 = q1*q5 + q5
      binary%d2val2 = q1*q6 + q6
      binary%d1val2_d1val3 = q1*q7 + q7
      binary%d2val3 = -0.5_dp*y%d2val3 + 0.5_dp*x%d2val3 + q2*(x%d2val3 - y%d2val3)
      binary%d3val1 = q1*q8 + q8
      binary%d2val1_d1val2 = q1*q9 + q9
      binary%d2val1_d1val3 = q1*q10 + q10
      binary%d1val1_d2val2 = q1*q11 + q11
      binary%d1val1_d1val2_d1val3 = q1*q12 + q12
      binary%d1val1_d2val3 = -0.5_dp*y%d1val1_d2val3 + 0.5_dp*x%d1val1_d2val3 + q2*(x%d1val1_d2val3 - y%d1val1_d2val3)
      binary%d3val2 = q1*q13 + q13
      binary%d2val2_d1val3 = q1*q14 + q14
      binary%d1val2_d2val3 = -0.5_dp*y%d1val2_d2val3 + 0.5_dp*x%d1val2_d2val3 + q2*(x%d1val2_d2val3 - y%d1val2_d2val3)
      binary%d3val1_d1val3 = q1*q15 + q15
      binary%d2val1_d1val2_d1val3 = q1*q16 + q16
      binary%d2val1_d2val3 = -0.5_dp*y%d2val1_d2val3 + 0.5_dp*x%d2val1_d2val3 + q2*(x%d2val1_d2val3 - y%d2val1_d2val3)
      binary%d1val1_d2val2_d1val3 = q1*q17 + q17
      binary%d1val1_d1val2_d2val3 = -0.5_dp*y%d1val1_d1val2_d2val3 + 0.5_dp*x%d1val1_d1val2_d2val3 + q2*(x%d1val1_d1val2_d2val3 - y%d1val1_d1val2_d2val3)
      binary%d3val2_d1val3 = q1*q18 + q18
      binary%d2val2_d2val3 = -0.5_dp*y%d2val2_d2val3 + 0.5_dp*x%d2val2_d2val3 + q2*(x%d2val2_d2val3 - y%d2val2_d2val3)
      binary%d3val1_d2val3 = -0.5_dp*y%d3val1_d2val3 + 0.5_dp*x%d3val1_d2val3 + q2*(x%d3val1_d2val3 - y%d3val1_d2val3)
      binary%d2val1_d1val2_d2val3 = -0.5_dp*y%d2val1_d1val2_d2val3 + 0.5_dp*x%d2val1_d1val2_d2val3 + q2*(x%d2val1_d1val2_d2val3 - y%d2val1_d1val2_d2val3)
      binary%d1val1_d2val2_d2val3 = -0.5_dp*y%d1val1_d2val2_d2val3 + 0.5_dp*x%d1val1_d2val2_d2val3 + q2*(x%d1val1_d2val2_d2val3 - y%d1val1_d2val2_d2val3)
      binary%d3val2_d2val3 = -0.5_dp*y%d3val2_d2val3 + 0.5_dp*x%d3val2_d2val3 + q2*(x%d3val2_d2val3 - y%d3val2_d2val3)
   end function dim_self
   
   function dim_self_real(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      real(dp), intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
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
      q0 = x%val - y
      q1 = 0.5_dp*x%d1val1
      q2 = sgn(q0)
      q3 = 0.5_dp*x%d1val2
      q4 = 0.5_dp*x%d1val3
      q5 = 0.5_dp*x%d2val1
      q6 = 0.5_dp*x%d1val1_d1val2
      q7 = 0.5_dp*x%d1val1_d1val3
      q8 = 0.5_dp*x%d2val2
      q9 = 0.5_dp*x%d1val2_d1val3
      q10 = 0.5_dp*x%d2val3
      q11 = 0.5_dp*x%d3val1
      q12 = 0.5_dp*x%d2val1_d1val2
      q13 = 0.5_dp*x%d2val1_d1val3
      q14 = 0.5_dp*x%d1val1_d2val2
      q15 = 0.5_dp*x%d1val1_d1val2_d1val3
      q16 = 0.5_dp*x%d1val1_d2val3
      q17 = 0.5_dp*x%d3val2
      q18 = 0.5_dp*x%d2val2_d1val3
      q19 = 0.5_dp*x%d1val2_d2val3
      q20 = 0.5_dp*x%d3val1_d1val3
      q21 = 0.5_dp*x%d2val1_d1val2_d1val3
      q22 = 0.5_dp*x%d2val1_d2val3
      q23 = 0.5_dp*x%d1val1_d2val2_d1val3
      q24 = 0.5_dp*x%d1val1_d1val2_d2val3
      q25 = 0.5_dp*x%d3val2_d1val3
      q26 = 0.5_dp*x%d2val2_d2val3
      q27 = 0.5_dp*x%d3val1_d2val3
      q28 = 0.5_dp*x%d2val1_d1val2_d2val3
      q29 = 0.5_dp*x%d1val1_d2val2_d2val3
      q30 = 0.5_dp*x%d3val2_d2val3
      unary%val = -0.5_dp*y + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*q2 + q1
      unary%d1val2 = q2*q3 + q3
      unary%d1val3 = q2*q4 + q4
      unary%d2val1 = q2*q5 + q5
      unary%d1val1_d1val2 = q2*q6 + q6
      unary%d1val1_d1val3 = q2*q7 + q7
      unary%d2val2 = q2*q8 + q8
      unary%d1val2_d1val3 = q2*q9 + q9
      unary%d2val3 = q10*q2 + q10
      unary%d3val1 = q11*q2 + q11
      unary%d2val1_d1val2 = q12*q2 + q12
      unary%d2val1_d1val3 = q13*q2 + q13
      unary%d1val1_d2val2 = q14*q2 + q14
      unary%d1val1_d1val2_d1val3 = q15*q2 + q15
      unary%d1val1_d2val3 = q16*q2 + q16
      unary%d3val2 = q17*q2 + q17
      unary%d2val2_d1val3 = q18*q2 + q18
      unary%d1val2_d2val3 = q19*q2 + q19
      unary%d3val1_d1val3 = q2*q20 + q20
      unary%d2val1_d1val2_d1val3 = q2*q21 + q21
      unary%d2val1_d2val3 = q2*q22 + q22
      unary%d1val1_d2val2_d1val3 = q2*q23 + q23
      unary%d1val1_d1val2_d2val3 = q2*q24 + q24
      unary%d3val2_d1val3 = q2*q25 + q25
      unary%d2val2_d2val3 = q2*q26 + q26
      unary%d3val1_d2val3 = q2*q27 + q27
      unary%d2val1_d1val2_d2val3 = q2*q28 + q28
      unary%d1val1_d2val2_d2val3 = q2*q29 + q29
      unary%d3val2_d2val3 = q2*q30 + q30
   end function dim_self_real
   
   function dim_real_self(z, x) result(unary)
      real(dp), intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
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
      q0 = x%val - z
      q1 = 0.5_dp*x%d1val1
      q2 = sgn(q0)
      q3 = 0.5_dp*x%d1val2
      q4 = 0.5_dp*x%d1val3
      q5 = 0.5_dp*x%d2val1
      q6 = 0.5_dp*x%d1val1_d1val2
      q7 = 0.5_dp*x%d1val1_d1val3
      q8 = 0.5_dp*x%d2val2
      q9 = 0.5_dp*x%d1val2_d1val3
      q10 = 0.5_dp*x%d2val3
      q11 = 0.5_dp*x%d3val1
      q12 = 0.5_dp*x%d2val1_d1val2
      q13 = 0.5_dp*x%d2val1_d1val3
      q14 = 0.5_dp*x%d1val1_d2val2
      q15 = 0.5_dp*x%d1val1_d1val2_d1val3
      q16 = 0.5_dp*x%d1val1_d2val3
      q17 = 0.5_dp*x%d3val2
      q18 = 0.5_dp*x%d2val2_d1val3
      q19 = 0.5_dp*x%d1val2_d2val3
      q20 = 0.5_dp*x%d3val1_d1val3
      q21 = 0.5_dp*x%d2val1_d1val2_d1val3
      q22 = 0.5_dp*x%d2val1_d2val3
      q23 = 0.5_dp*x%d1val1_d2val2_d1val3
      q24 = 0.5_dp*x%d1val1_d1val2_d2val3
      q25 = 0.5_dp*x%d3val2_d1val3
      q26 = 0.5_dp*x%d2val2_d2val3
      q27 = 0.5_dp*x%d3val1_d2val3
      q28 = 0.5_dp*x%d2val1_d1val2_d2val3
      q29 = 0.5_dp*x%d1val1_d2val2_d2val3
      q30 = 0.5_dp*x%d3val2_d2val3
      unary%val = -0.5_dp*x%val + 0.5_dp*z + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*q2 - q1
      unary%d1val2 = q2*q3 - q3
      unary%d1val3 = q2*q4 - q4
      unary%d2val1 = q2*q5 - q5
      unary%d1val1_d1val2 = q2*q6 - q6
      unary%d1val1_d1val3 = q2*q7 - q7
      unary%d2val2 = q2*q8 - q8
      unary%d1val2_d1val3 = q2*q9 - q9
      unary%d2val3 = q10*q2 - q10
      unary%d3val1 = q11*q2 - q11
      unary%d2val1_d1val2 = q12*q2 - q12
      unary%d2val1_d1val3 = q13*q2 - q13
      unary%d1val1_d2val2 = q14*q2 - q14
      unary%d1val1_d1val2_d1val3 = q15*q2 - q15
      unary%d1val1_d2val3 = q16*q2 - q16
      unary%d3val2 = q17*q2 - q17
      unary%d2val2_d1val3 = q18*q2 - q18
      unary%d1val2_d2val3 = q19*q2 - q19
      unary%d3val1_d1val3 = q2*q20 - q20
      unary%d2val1_d1val2_d1val3 = q2*q21 - q21
      unary%d2val1_d2val3 = q2*q22 - q22
      unary%d1val1_d2val2_d1val3 = q2*q23 - q23
      unary%d1val1_d1val2_d2val3 = q2*q24 - q24
      unary%d3val2_d1val3 = q2*q25 - q25
      unary%d2val2_d2val3 = q2*q26 - q26
      unary%d3val1_d2val3 = q2*q27 - q27
      unary%d2val1_d1val2_d2val3 = q2*q28 - q28
      unary%d1val1_d2val2_d2val3 = q2*q29 - q29
      unary%d3val2_d2val3 = q2*q30 - q30
   end function dim_real_self
   
   function dim_self_int(x, y) result(unary)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      integer, intent(in) :: y
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
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
      y_dp = y
      q0 = x%val - y_dp
      q1 = 0.5_dp*x%d1val1
      q2 = sgn(q0)
      q3 = 0.5_dp*x%d1val2
      q4 = 0.5_dp*x%d1val3
      q5 = 0.5_dp*x%d2val1
      q6 = 0.5_dp*x%d1val1_d1val2
      q7 = 0.5_dp*x%d1val1_d1val3
      q8 = 0.5_dp*x%d2val2
      q9 = 0.5_dp*x%d1val2_d1val3
      q10 = 0.5_dp*x%d2val3
      q11 = 0.5_dp*x%d3val1
      q12 = 0.5_dp*x%d2val1_d1val2
      q13 = 0.5_dp*x%d2val1_d1val3
      q14 = 0.5_dp*x%d1val1_d2val2
      q15 = 0.5_dp*x%d1val1_d1val2_d1val3
      q16 = 0.5_dp*x%d1val1_d2val3
      q17 = 0.5_dp*x%d3val2
      q18 = 0.5_dp*x%d2val2_d1val3
      q19 = 0.5_dp*x%d1val2_d2val3
      q20 = 0.5_dp*x%d3val1_d1val3
      q21 = 0.5_dp*x%d2val1_d1val2_d1val3
      q22 = 0.5_dp*x%d2val1_d2val3
      q23 = 0.5_dp*x%d1val1_d2val2_d1val3
      q24 = 0.5_dp*x%d1val1_d1val2_d2val3
      q25 = 0.5_dp*x%d3val2_d1val3
      q26 = 0.5_dp*x%d2val2_d2val3
      q27 = 0.5_dp*x%d3val1_d2val3
      q28 = 0.5_dp*x%d2val1_d1val2_d2val3
      q29 = 0.5_dp*x%d1val1_d2val2_d2val3
      q30 = 0.5_dp*x%d3val2_d2val3
      unary%val = -0.5_dp*y_dp + 0.5_dp*x%val + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*q2 + q1
      unary%d1val2 = q2*q3 + q3
      unary%d1val3 = q2*q4 + q4
      unary%d2val1 = q2*q5 + q5
      unary%d1val1_d1val2 = q2*q6 + q6
      unary%d1val1_d1val3 = q2*q7 + q7
      unary%d2val2 = q2*q8 + q8
      unary%d1val2_d1val3 = q2*q9 + q9
      unary%d2val3 = q10*q2 + q10
      unary%d3val1 = q11*q2 + q11
      unary%d2val1_d1val2 = q12*q2 + q12
      unary%d2val1_d1val3 = q13*q2 + q13
      unary%d1val1_d2val2 = q14*q2 + q14
      unary%d1val1_d1val2_d1val3 = q15*q2 + q15
      unary%d1val1_d2val3 = q16*q2 + q16
      unary%d3val2 = q17*q2 + q17
      unary%d2val2_d1val3 = q18*q2 + q18
      unary%d1val2_d2val3 = q19*q2 + q19
      unary%d3val1_d1val3 = q2*q20 + q20
      unary%d2val1_d1val2_d1val3 = q2*q21 + q21
      unary%d2val1_d2val3 = q2*q22 + q22
      unary%d1val1_d2val2_d1val3 = q2*q23 + q23
      unary%d1val1_d1val2_d2val3 = q2*q24 + q24
      unary%d3val2_d1val3 = q2*q25 + q25
      unary%d2val2_d2val3 = q2*q26 + q26
      unary%d3val1_d2val3 = q2*q27 + q27
      unary%d2val1_d1val2_d2val3 = q2*q28 + q28
      unary%d1val1_d2val2_d2val3 = q2*q29 + q29
      unary%d3val2_d2val3 = q2*q30 + q30
   end function dim_self_int
   
   function dim_int_self(z, x) result(unary)
      integer, intent(in) :: z
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: x
      type(auto_diff_real_2var_order3_1var_order2) :: unary
      real(dp) :: y_dp
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
      y_dp = z
      q0 = x%val - y_dp
      q1 = 0.5_dp*x%d1val1
      q2 = sgn(q0)
      q3 = 0.5_dp*x%d1val2
      q4 = 0.5_dp*x%d1val3
      q5 = 0.5_dp*x%d2val1
      q6 = 0.5_dp*x%d1val1_d1val2
      q7 = 0.5_dp*x%d1val1_d1val3
      q8 = 0.5_dp*x%d2val2
      q9 = 0.5_dp*x%d1val2_d1val3
      q10 = 0.5_dp*x%d2val3
      q11 = 0.5_dp*x%d3val1
      q12 = 0.5_dp*x%d2val1_d1val2
      q13 = 0.5_dp*x%d2val1_d1val3
      q14 = 0.5_dp*x%d1val1_d2val2
      q15 = 0.5_dp*x%d1val1_d1val2_d1val3
      q16 = 0.5_dp*x%d1val1_d2val3
      q17 = 0.5_dp*x%d3val2
      q18 = 0.5_dp*x%d2val2_d1val3
      q19 = 0.5_dp*x%d1val2_d2val3
      q20 = 0.5_dp*x%d3val1_d1val3
      q21 = 0.5_dp*x%d2val1_d1val2_d1val3
      q22 = 0.5_dp*x%d2val1_d2val3
      q23 = 0.5_dp*x%d1val1_d2val2_d1val3
      q24 = 0.5_dp*x%d1val1_d1val2_d2val3
      q25 = 0.5_dp*x%d3val2_d1val3
      q26 = 0.5_dp*x%d2val2_d2val3
      q27 = 0.5_dp*x%d3val1_d2val3
      q28 = 0.5_dp*x%d2val1_d1val2_d2val3
      q29 = 0.5_dp*x%d1val1_d2val2_d2val3
      q30 = 0.5_dp*x%d3val2_d2val3
      unary%val = -0.5_dp*x%val + 0.5_dp*y_dp + 0.5_dp*Abs(q0)
      unary%d1val1 = q1*q2 - q1
      unary%d1val2 = q2*q3 - q3
      unary%d1val3 = q2*q4 - q4
      unary%d2val1 = q2*q5 - q5
      unary%d1val1_d1val2 = q2*q6 - q6
      unary%d1val1_d1val3 = q2*q7 - q7
      unary%d2val2 = q2*q8 - q8
      unary%d1val2_d1val3 = q2*q9 - q9
      unary%d2val3 = q10*q2 - q10
      unary%d3val1 = q11*q2 - q11
      unary%d2val1_d1val2 = q12*q2 - q12
      unary%d2val1_d1val3 = q13*q2 - q13
      unary%d1val1_d2val2 = q14*q2 - q14
      unary%d1val1_d1val2_d1val3 = q15*q2 - q15
      unary%d1val1_d2val3 = q16*q2 - q16
      unary%d3val2 = q17*q2 - q17
      unary%d2val2_d1val3 = q18*q2 - q18
      unary%d1val2_d2val3 = q19*q2 - q19
      unary%d3val1_d1val3 = q2*q20 - q20
      unary%d2val1_d1val2_d1val3 = q2*q21 - q21
      unary%d2val1_d2val3 = q2*q22 - q22
      unary%d1val1_d2val2_d1val3 = q2*q23 - q23
      unary%d1val1_d1val2_d2val3 = q2*q24 - q24
      unary%d3val2_d1val3 = q2*q25 - q25
      unary%d2val2_d2val3 = q2*q26 - q26
      unary%d3val1_d2val3 = q2*q27 - q27
      unary%d2val1_d1val2_d2val3 = q2*q28 - q28
      unary%d1val1_d2val2_d2val3 = q2*q29 - q29
      unary%d3val2_d2val3 = q2*q30 - q30
   end function dim_int_self
   
   function differentiate_auto_diff_real_2var_order3_1var_order2_1(this) result(derivative)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2) :: derivative
      derivative%val = this%d1val1
      derivative%d1val1 = this%d2val1
      derivative%d1val2 = this%d1val1_d1val2
      derivative%d1val3 = this%d1val1_d1val3
      derivative%d2val1 = this%d3val1
      derivative%d1val1_d1val2 = this%d2val1_d1val2
      derivative%d1val1_d1val3 = this%d2val1_d1val3
      derivative%d2val2 = this%d1val1_d2val2
      derivative%d1val2_d1val3 = this%d1val1_d1val2_d1val3
      derivative%d2val3 = this%d1val1_d2val3
      derivative%d3val1 = 0_dp
      derivative%d2val1_d1val2 = 0_dp
      derivative%d2val1_d1val3 = this%d3val1_d1val3
      derivative%d1val1_d2val2 = 0_dp
      derivative%d1val1_d1val2_d1val3 = this%d2val1_d1val2_d1val3
      derivative%d1val1_d2val3 = this%d2val1_d2val3
      derivative%d3val2 = 0_dp
      derivative%d2val2_d1val3 = this%d1val1_d2val2_d1val3
      derivative%d1val2_d2val3 = this%d1val1_d1val2_d2val3
      derivative%d3val1_d1val3 = 0_dp
      derivative%d2val1_d1val2_d1val3 = 0_dp
      derivative%d2val1_d2val3 = this%d3val1_d2val3
      derivative%d1val1_d2val2_d1val3 = 0_dp
      derivative%d1val1_d1val2_d2val3 = this%d2val1_d1val2_d2val3
      derivative%d3val2_d1val3 = 0_dp
      derivative%d2val2_d2val3 = this%d1val1_d2val2_d2val3
      derivative%d3val1_d2val3 = 0_dp
      derivative%d2val1_d1val2_d2val3 = 0_dp
      derivative%d1val1_d2val2_d2val3 = 0_dp
      derivative%d3val2_d2val3 = 0_dp
   end function differentiate_auto_diff_real_2var_order3_1var_order2_1
   
   function differentiate_auto_diff_real_2var_order3_1var_order2_2(this) result(derivative)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2) :: derivative
      derivative%val = this%d1val2
      derivative%d1val1 = this%d1val1_d1val2
      derivative%d1val2 = this%d2val2
      derivative%d1val3 = this%d1val2_d1val3
      derivative%d2val1 = this%d2val1_d1val2
      derivative%d1val1_d1val2 = this%d1val1_d2val2
      derivative%d1val1_d1val3 = this%d1val1_d1val2_d1val3
      derivative%d2val2 = this%d3val2
      derivative%d1val2_d1val3 = this%d2val2_d1val3
      derivative%d2val3 = this%d1val2_d2val3
      derivative%d3val1 = 0_dp
      derivative%d2val1_d1val2 = 0_dp
      derivative%d2val1_d1val3 = this%d2val1_d1val2_d1val3
      derivative%d1val1_d2val2 = 0_dp
      derivative%d1val1_d1val2_d1val3 = this%d1val1_d2val2_d1val3
      derivative%d1val1_d2val3 = this%d1val1_d1val2_d2val3
      derivative%d3val2 = 0_dp
      derivative%d2val2_d1val3 = this%d3val2_d1val3
      derivative%d1val2_d2val3 = this%d2val2_d2val3
      derivative%d3val1_d1val3 = 0_dp
      derivative%d2val1_d1val2_d1val3 = 0_dp
      derivative%d2val1_d2val3 = this%d2val1_d1val2_d2val3
      derivative%d1val1_d2val2_d1val3 = 0_dp
      derivative%d1val1_d1val2_d2val3 = this%d1val1_d2val2_d2val3
      derivative%d3val2_d1val3 = 0_dp
      derivative%d2val2_d2val3 = this%d3val2_d2val3
      derivative%d3val1_d2val3 = 0_dp
      derivative%d2val1_d1val2_d2val3 = 0_dp
      derivative%d1val1_d2val2_d2val3 = 0_dp
      derivative%d3val2_d2val3 = 0_dp
   end function differentiate_auto_diff_real_2var_order3_1var_order2_2
   
   function differentiate_auto_diff_real_2var_order3_1var_order2_3(this) result(derivative)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: this
      type(auto_diff_real_2var_order3_1var_order2) :: derivative
      derivative%val = this%d1val3
      derivative%d1val1 = this%d1val1_d1val3
      derivative%d1val2 = this%d1val2_d1val3
      derivative%d1val3 = this%d2val3
      derivative%d2val1 = this%d2val1_d1val3
      derivative%d1val1_d1val2 = this%d1val1_d1val2_d1val3
      derivative%d1val1_d1val3 = this%d1val1_d2val3
      derivative%d2val2 = this%d2val2_d1val3
      derivative%d1val2_d1val3 = this%d1val2_d2val3
      derivative%d2val3 = 0_dp
      derivative%d3val1 = this%d3val1_d1val3
      derivative%d2val1_d1val2 = this%d2val1_d1val2_d1val3
      derivative%d2val1_d1val3 = this%d2val1_d2val3
      derivative%d1val1_d2val2 = this%d1val1_d2val2_d1val3
      derivative%d1val1_d1val2_d1val3 = this%d1val1_d1val2_d2val3
      derivative%d1val1_d2val3 = 0_dp
      derivative%d3val2 = this%d3val2_d1val3
      derivative%d2val2_d1val3 = this%d2val2_d2val3
      derivative%d1val2_d2val3 = 0_dp
      derivative%d3val1_d1val3 = this%d3val1_d2val3
      derivative%d2val1_d1val2_d1val3 = this%d2val1_d1val2_d2val3
      derivative%d2val1_d2val3 = 0_dp
      derivative%d1val1_d2val2_d1val3 = this%d1val1_d2val2_d2val3
      derivative%d1val1_d1val2_d2val3 = 0_dp
      derivative%d3val2_d1val3 = this%d3val2_d2val3
      derivative%d2val2_d2val3 = 0_dp
      derivative%d3val1_d2val3 = 0_dp
      derivative%d2val1_d1val2_d2val3 = 0_dp
      derivative%d1val1_d2val2_d2val3 = 0_dp
      derivative%d3val2_d2val3 = 0_dp
   end function differentiate_auto_diff_real_2var_order3_1var_order2_3
   
end module auto_diff_real_2var_order3_1var_order2_module