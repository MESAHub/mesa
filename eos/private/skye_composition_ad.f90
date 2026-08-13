! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module skye_composition_ad

      use const_def, only: dp
      use auto_diff

      implicit none

      private
      public :: skye_composition_ad_real
      public :: make_skye_composition_ad
      public :: assignment(=)
      public :: operator(+), operator(-), operator(*), operator(/)
      public :: operator(.eq.), operator(.ne.), operator(.gt.)
      public :: operator(.lt.), operator(.ge.), operator(.le.)
      public :: sqrt, log, exp, atan, tanh, cosh, sinh, abs
      public :: pow, pow2, pow3, pow4, pow5, min, max
      public :: differentiate_1, differentiate_2

      type :: skye_composition_ad_real
         type(auto_diff_real_2var_order3) :: val
         type(auto_diff_real_2var_order3) :: d
      end type skye_composition_ad_real

      interface assignment(=)
         module procedure assign_from_real
         module procedure assign_from_int
         module procedure assign_from_auto_diff
      end interface assignment(=)

      interface operator(+)
         module procedure add_self
         module procedure add_self_real
         module procedure add_real_self
         module procedure add_self_int
         module procedure add_int_self
         module procedure add_self_ad
         module procedure add_ad_self
      end interface operator(+)

      interface operator(-)
         module procedure unary_minus_self
         module procedure sub_self
         module procedure sub_self_real
         module procedure sub_real_self
         module procedure sub_self_int
         module procedure sub_int_self
         module procedure sub_self_ad
         module procedure sub_ad_self
      end interface operator(-)

      interface operator(*)
         module procedure multiply_self
         module procedure multiply_self_real
         module procedure multiply_real_self
         module procedure multiply_self_int
         module procedure multiply_int_self
         module procedure multiply_self_ad
         module procedure multiply_ad_self
      end interface operator(*)

      interface operator(/)
         module procedure divide_self
         module procedure divide_self_real
         module procedure divide_real_self
         module procedure divide_self_int
         module procedure divide_int_self
         module procedure divide_self_ad
         module procedure divide_ad_self
      end interface operator(/)

      interface operator(.eq.)
         module procedure equal_self
         module procedure equal_self_real
         module procedure equal_real_self
      end interface operator(.eq.)

      interface operator(.ne.)
         module procedure neq_self
         module procedure neq_self_real
         module procedure neq_real_self
      end interface operator(.ne.)

      interface operator(.gt.)
         module procedure gt_self
         module procedure gt_self_real
         module procedure gt_real_self
      end interface operator(.gt.)

      interface operator(.lt.)
         module procedure lt_self
         module procedure lt_self_real
         module procedure lt_real_self
      end interface operator(.lt.)

      interface operator(.ge.)
         module procedure ge_self
         module procedure ge_self_real
         module procedure ge_real_self
      end interface operator(.ge.)

      interface operator(.le.)
         module procedure le_self
         module procedure le_self_real
         module procedure le_real_self
      end interface operator(.le.)

      interface sqrt
         module procedure sqrt_self
      end interface sqrt

      interface log
         module procedure log_self
      end interface log

      interface exp
         module procedure exp_self
      end interface exp

      interface atan
         module procedure atan_self
      end interface atan

      interface tanh
         module procedure tanh_self
      end interface tanh

      interface cosh
         module procedure cosh_self
      end interface cosh

      interface sinh
         module procedure sinh_self
      end interface sinh

      interface abs
         module procedure abs_self
      end interface abs

      interface pow
         module procedure pow_self_real
         module procedure pow_self
         module procedure pow_real_self
      end interface pow

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

      interface min
         module procedure min_self
         module procedure min_self_real
         module procedure min_real_self
      end interface min

      interface max
         module procedure max_self
         module procedure max_self_real
         module procedure max_real_self
      end interface max

      interface differentiate_1
         module procedure differentiate_1_self
      end interface differentiate_1

      interface differentiate_2
         module procedure differentiate_2_self
      end interface differentiate_2

      contains

      function make_skye_composition_ad(val, d) result(x)
         type(auto_diff_real_2var_order3), intent(in) :: val, d
         type(skye_composition_ad_real) :: x
         x% val = val
         x% d = d
      end function make_skye_composition_ad


      subroutine assign_from_real(this, other)
         type(skye_composition_ad_real), intent(out) :: this
         real(dp), intent(in) :: other
         this% val = other
         this% d = 0d0
      end subroutine assign_from_real


      subroutine assign_from_int(this, other)
         type(skye_composition_ad_real), intent(out) :: this
         integer, intent(in) :: other
         this% val = other
         this% d = 0d0
      end subroutine assign_from_int


      subroutine assign_from_auto_diff(this, other)
         type(skye_composition_ad_real), intent(out) :: this
         type(auto_diff_real_2var_order3), intent(in) :: other
         this% val = other
         this% d = 0d0
      end subroutine assign_from_auto_diff


      function add_self(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x, y
         type(skye_composition_ad_real) :: z
         z% val = x% val + y% val
         z% d = x% d + y% d
      end function add_self


      function add_self_real(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x% val + y
         z% d = x% d
      end function add_self_real


      function add_real_self(x, y) result(z)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x + y% val
         z% d = y% d
      end function add_real_self


      function add_self_int(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         integer, intent(in) :: y
         type(skye_composition_ad_real) :: z
         z = x + dble(y)
      end function add_self_int


      function add_int_self(x, y) result(z)
         integer, intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z = dble(x) + y
      end function add_int_self


      function add_self_ad(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(auto_diff_real_2var_order3), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x% val + y
         z% d = x% d
      end function add_self_ad


      function add_ad_self(x, y) result(z)
         type(auto_diff_real_2var_order3), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x + y% val
         z% d = y% d
      end function add_ad_self


      function unary_minus_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = -x% val
         z% d = -x% d
      end function unary_minus_self


      function sub_self(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x, y
         type(skye_composition_ad_real) :: z
         z% val = x% val - y% val
         z% d = x% d - y% d
      end function sub_self


      function sub_self_real(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x% val - y
         z% d = x% d
      end function sub_self_real


      function sub_real_self(x, y) result(z)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x - y% val
         z% d = -y% d
      end function sub_real_self


      function sub_self_int(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         integer, intent(in) :: y
         type(skye_composition_ad_real) :: z
         z = x - dble(y)
      end function sub_self_int


      function sub_int_self(x, y) result(z)
         integer, intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z = dble(x) - y
      end function sub_int_self


      function sub_self_ad(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(auto_diff_real_2var_order3), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x% val - y
         z% d = x% d
      end function sub_self_ad


      function sub_ad_self(x, y) result(z)
         type(auto_diff_real_2var_order3), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x - y% val
         z% d = -y% d
      end function sub_ad_self


      function multiply_self(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x, y
         type(skye_composition_ad_real) :: z
         z% val = x% val*y% val
         z% d = x% d*y% val + x% val*y% d
      end function multiply_self


      function multiply_self_real(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x% val*y
         z% d = x% d*y
      end function multiply_self_real


      function multiply_real_self(x, y) result(z)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x*y% val
         z% d = x*y% d
      end function multiply_real_self


      function multiply_self_int(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         integer, intent(in) :: y
         type(skye_composition_ad_real) :: z
         z = x*dble(y)
      end function multiply_self_int


      function multiply_int_self(x, y) result(z)
         integer, intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z = dble(x)*y
      end function multiply_int_self


      function multiply_self_ad(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(auto_diff_real_2var_order3), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x% val*y
         z% d = x% d*y
      end function multiply_self_ad


      function multiply_ad_self(x, y) result(z)
         type(auto_diff_real_2var_order3), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x*y% val
         z% d = x*y% d
      end function multiply_ad_self


      function divide_self(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x, y
         type(skye_composition_ad_real) :: z
         z% val = x% val/y% val
         z% d = (x% d*y% val - x% val*y% d)/pow2(y% val)
      end function divide_self


      function divide_self_real(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x% val/y
         z% d = x% d/y
      end function divide_self_real


      function divide_real_self(x, y) result(z)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x/y% val
         z% d = -x*y% d/pow2(y% val)
      end function divide_real_self


      function divide_self_int(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         integer, intent(in) :: y
         type(skye_composition_ad_real) :: z
         z = x/dble(y)
      end function divide_self_int


      function divide_int_self(x, y) result(z)
         integer, intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z = dble(x)/y
      end function divide_int_self


      function divide_self_ad(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(auto_diff_real_2var_order3), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x% val/y
         z% d = x% d/y
      end function divide_self_ad


      function divide_ad_self(x, y) result(z)
         type(auto_diff_real_2var_order3), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         z% val = x/y% val
         z% d = -x*y% d/pow2(y% val)
      end function divide_ad_self


      logical function equal_self(x, y)
         type(skye_composition_ad_real), intent(in) :: x, y
         equal_self = x% val == y% val
      end function equal_self


      logical function equal_self_real(x, y)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         equal_self_real = x% val == y
      end function equal_self_real


      logical function equal_real_self(x, y)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         equal_real_self = x == y% val
      end function equal_real_self


      logical function neq_self(x, y)
         type(skye_composition_ad_real), intent(in) :: x, y
         neq_self = x% val /= y% val
      end function neq_self


      logical function neq_self_real(x, y)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         neq_self_real = x% val /= y
      end function neq_self_real


      logical function neq_real_self(x, y)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         neq_real_self = x /= y% val
      end function neq_real_self


      logical function gt_self(x, y)
         type(skye_composition_ad_real), intent(in) :: x, y
         gt_self = x% val > y% val
      end function gt_self


      logical function gt_self_real(x, y)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         gt_self_real = x% val > y
      end function gt_self_real


      logical function gt_real_self(x, y)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         gt_real_self = x > y% val
      end function gt_real_self


      logical function lt_self(x, y)
         type(skye_composition_ad_real), intent(in) :: x, y
         lt_self = x% val < y% val
      end function lt_self


      logical function lt_self_real(x, y)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         lt_self_real = x% val < y
      end function lt_self_real


      logical function lt_real_self(x, y)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         lt_real_self = x < y% val
      end function lt_real_self


      logical function ge_self(x, y)
         type(skye_composition_ad_real), intent(in) :: x, y
         ge_self = x% val >= y% val
      end function ge_self


      logical function ge_self_real(x, y)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         ge_self_real = x% val >= y
      end function ge_self_real


      logical function ge_real_self(x, y)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         ge_real_self = x >= y% val
      end function ge_real_self


      logical function le_self(x, y)
         type(skye_composition_ad_real), intent(in) :: x, y
         le_self = x% val <= y% val
      end function le_self


      logical function le_self_real(x, y)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         le_self_real = x% val <= y
      end function le_self_real


      logical function le_real_self(x, y)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         le_real_self = x <= y% val
      end function le_real_self


      function sqrt_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = sqrt(x% val)
         z% d = 0.5d0*x% d/z% val
      end function sqrt_self


      function log_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = log(x% val)
         z% d = x% d/x% val
      end function log_self


      function exp_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = exp(x% val)
         z% d = z% val*x% d
      end function exp_self


      function atan_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = atan(x% val)
         z% d = x% d/(1d0 + pow2(x% val))
      end function atan_self


      function tanh_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = tanh(x% val)
         z% d = (1d0 - pow2(z% val))*x% d
      end function tanh_self


      function cosh_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = cosh(x% val)
         z% d = sinh(x% val)*x% d
      end function cosh_self


      function sinh_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = sinh(x% val)
         z% d = cosh(x% val)*x% d
      end function sinh_self


      function abs_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         if (x% val >= 0d0) then
            z = x
         else
            z = -x
         end if
      end function abs_self


      function pow_self_real(x, a) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: a
         type(skye_composition_ad_real) :: z
         z% val = pow(x% val, a)
         z% d = a*pow(x% val, a - 1d0)*x% d
      end function pow_self_real


      function pow_self(x, a) result(z)
         type(skye_composition_ad_real), intent(in) :: x, a
         type(skye_composition_ad_real) :: z
         z% val = pow(x% val, a% val)
         z% d = z% val*(a% d*log(x% val) + a% val*x% d/x% val)
      end function pow_self


      function pow_real_self(x, a) result(z)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: a
         type(skye_composition_ad_real) :: z
         z% val = pow(x, a% val)
         z% d = z% val*log(x)*a% d
      end function pow_real_self


      function pow2_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = pow2(x% val)
         z% d = 2d0*x% val*x% d
      end function pow2_self


      function pow3_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = pow3(x% val)
         z% d = 3d0*pow2(x% val)*x% d
      end function pow3_self


      function pow4_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = pow4(x% val)
         z% d = 4d0*pow3(x% val)*x% d
      end function pow4_self


      function pow5_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = pow5(x% val)
         z% d = 5d0*pow4(x% val)*x% d
      end function pow5_self


      function min_self(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x, y
         type(skye_composition_ad_real) :: z
         if (x% val <= y% val) then
            z = x
         else
            z = y
         end if
      end function min_self


      function min_self_real(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         type(skye_composition_ad_real) :: z
         if (x% val <= y) then
            z = x
         else
            z = y
         end if
      end function min_self_real


      function min_real_self(x, y) result(z)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         if (x <= y% val) then
            z = x
         else
            z = y
         end if
      end function min_real_self


      function max_self(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x, y
         type(skye_composition_ad_real) :: z
         if (x% val >= y% val) then
            z = x
         else
            z = y
         end if
      end function max_self


      function max_self_real(x, y) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         real(dp), intent(in) :: y
         type(skye_composition_ad_real) :: z
         if (x% val >= y) then
            z = x
         else
            z = y
         end if
      end function max_self_real


      function max_real_self(x, y) result(z)
         real(dp), intent(in) :: x
         type(skye_composition_ad_real), intent(in) :: y
         type(skye_composition_ad_real) :: z
         if (x >= y% val) then
            z = x
         else
            z = y
         end if
      end function max_real_self


      function differentiate_1_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = differentiate_1(x% val)
         z% d = differentiate_1(x% d)
      end function differentiate_1_self


      function differentiate_2_self(x) result(z)
         type(skye_composition_ad_real), intent(in) :: x
         type(skye_composition_ad_real) :: z
         z% val = differentiate_2(x% val)
         z% d = differentiate_2(x% d)
      end function differentiate_2_self

      end module skye_composition_ad
