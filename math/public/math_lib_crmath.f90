! ***********************************************************************
!
!   Copyright (C) 2010-2020  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module math_lib

  ! Uses

  use const_lib, only: dp

  use math_pown
  use math_io
  use utils_lib
  use math_def

  use crmath

  use IEEE_ARITHMETIC

  ! No implicit typing

  implicit none

  ! Parameter definitions

  character(LEN=16), parameter :: MATH_BACKEND = 'CRMATH'

  ! Module variables

  real(dp), save :: ln10_m

  ! Interfaces

  interface safe_sqrt
     module procedure safe_sqrt_
  end interface safe_sqrt

  interface safe_log
     module procedure safe_log_
  end interface safe_log

  interface safe_log10
     module procedure safe_log10_
  end interface safe_log10

  interface exp10
     module procedure exp10_
  end interface exp10

  interface pow
     module procedure pow_i_
     module procedure pow_r_
  end interface pow

  ! Access specifiers

  private

  public :: MATH_BACKEND

  public :: math_init
  public :: safe_sqrt
  public :: log
  public :: safe_log
  public :: log10
  public :: safe_log10
  public :: log1p
  public :: log2
  public :: exp
  public :: exp10
  public :: expm1
  public :: powm1
  public :: pow
  public :: pow2
  public :: pow3
  public :: pow4
  public :: pow5
  public :: pow6
  public :: pow7
  public :: pow8
  public :: cos
  public :: sin
  public :: tan
  public :: cospi
  public :: sinpi
  public :: tanpi
  public :: acos
  public :: asin
  public :: atan
  public :: acospi
  public :: asinpi
  public :: atanpi
  public :: cosh
  public :: sinh
  public :: tanh
  public :: acosh
  public :: asinh
  public :: atanh
  public :: str_to_vector
  public :: str_to_double
  public :: double_to_str

  ! Procedures

contains

  subroutine math_init ()

    call crmath_init()

    ln10_m = log(10._dp)

    call precompute_some_zs()

  end subroutine math_init

  !****

  elemental function safe_sqrt_ (x) result (sqrt_x)

    real(dp), intent(in) :: x
    real(dp)             :: sqrt_x

    sqrt_x = SQRT(MAX(x, 0._dp))

  end function safe_sqrt_

  !****

  elemental function safe_log_ (x) result (log_x)

    real(dp), intent(in) :: x
    real(dp)             :: log_x

    if (.NOT. IEEE_IS_FINITE(x)) then

       log_x = -99._dp

    else

       log_x = log(MAX(1.E-99_dp, x))

    end if

  end function safe_log_

  !****

  elemental function safe_log10_ (x) result (log10_x)

    real(dp), intent(in) :: x
    real(dp)             :: log10_x

    if (.NOT. IEEE_IS_FINITE(x)) then

       log10_x = -99._dp

    else

       log10_x = log10(MAX(1.E-99_dp, x))

    end if

  end function safe_log10_

  !****

  elemental function exp10_ (x) result (exp10_x)

    real(dp), intent(in) :: x
    real(dp)             :: exp10_x

    integer :: ix
    integer :: i

    ix = FLOOR(x)

    if (x == ix) then ! integer power of 10

       exp10_x = 1._dp

       do i = 1, ABS(ix)
          exp10_x = exp10_x*10._dp
       end do

       if (ix < 0) exp10_x = 1._dp/exp10_x

    else

       exp10_x = exp(x*ln10_m)

    endif

  end function exp10_

  !****

  elemental function pow_i_ (x, iy) result (pow_x)

    real(dp), intent(in) :: x
    integer, intent(in)  :: iy
    real(dp)             :: pow_x

    integer :: i

    if (x == 0._dp) then

       pow_x = 0._dp

    else

       pow_x = 1._dp

       do i = 1, ABS(iy)
          pow_x = pow_x*x
       end do

       if (iy < 0) pow_x = 1._dp/pow_x

    endif

  end function pow_i_

  !****

  elemental function pow_r_ (x, y) result (pow_x)

    real(dp), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp)             :: pow_x

    integer :: iy
    integer :: i

    if (x == 0._dp) then

       pow_x = 0._dp

    else

       iy = floor(y)

       if (y == iy .AND. ABS(iy) < 100) then ! integer power of x

          pow_x = 1._dp

          do i = 1, ABS(iy)
             pow_x = pow_x*x
          end do

          if (iy < 0) pow_x = 1._dp/pow_x

       else

          pow_x = exp(log(x)*y)

       end if

    end if

  end function pow_r_


  include 'precompute_zs.inc'


end module math_lib
