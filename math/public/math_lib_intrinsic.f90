! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

  use const_lib, only: dp, PI

  use math_io
  use math_pown
  use math_def

  use IEEE_ARITHMETIC

  ! No implicit typing

  implicit none

  ! Parameter definitions

  character(LEN=16), parameter :: MATH_BACKEND = 'INTRINSIC'

  ! Interfaces

  ! Generic interfaces

  interface safe_sqrt
     module procedure safe_sqrt_
  end interface safe_sqrt

  interface safe_log
     module procedure safe_log_
  end interface safe_log

  interface safe_log10
     module procedure safe_log10_
  end interface safe_log10

  interface log1p
     module procedure log1p_
  end interface log1p

  interface log2
     module procedure log2_
  end interface log2

  interface exp10
     module procedure exp10_
  end interface exp10

  interface expm1
     module procedure expm1_
  end interface expm1

  interface pow
     module procedure pow_i_
     module procedure pow_r_
  end interface pow

  interface cospi
     module procedure cospi_
  end interface cospi

  interface sinpi
     module procedure sinpi_
  end interface sinpi

  interface tanpi
     module procedure tanpi_
  end interface tanpi

  interface acospi
     module procedure acospi_
  end interface acospi

  interface asinpi
     module procedure asinpi_
  end interface asinpi

  interface atanpi
     module procedure atanpi_
  end interface atanpi

  ! Module variables

  real(dp), save :: ln10_m

  ! Access specifiers

  private

  public :: MATH_BACKEND

  public :: math_init
  public :: safe_sqrt
  public :: safe_log
  public :: safe_log10
  public :: log1p
  public :: log2
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
  public :: cospi
  public :: sinpi
  public :: tanpi
  public :: acospi
  public :: asinpi
  public :: atanpi
  public :: str_to_vector
  public :: str_to_double
  public :: double_to_str

  ! Procedures

contains

  subroutine math_init ()

    ln10_m = LOG(10._dp)

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

       log_x = LOG(MAX(1.E-99_dp, x))

    end if

  end function safe_log_

  !****

  elemental function safe_log10_ (x) result (log10_x)

    real(dp), intent(in) :: x
    real(dp)             :: log10_x

    if (.NOT. IEEE_IS_FINITE(x)) then

       log10_x = -99._dp

    else

       log10_x = LOG10(MAX(1.E-99_dp, x))

    end if

  end function safe_log10_

  !****

  elemental function log1p_ (x) result (log1p_x)

    real(dp), intent(in) :: x
    real(dp)             :: log1p_x

    log1p_x = LOG(1._dp + x)

  end function log1p_

  !****

  elemental function log2_ (x) result (log2_x)

    real(dp), intent(in) :: x
    real(dp)             :: log2_x

    log2_x = LOG(x)/LOG(2._dp)

  end function log2_

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

       exp10_x = EXP(x*ln10_m)

    endif

  end function exp10_

  !****

  elemental function expm1_ (x) result (expm1_x)

    real(dp), intent(in) :: x
    real(dp)             :: expm1_x

    expm1_x = EXP(x) - 1._dp

  end function expm1_

  !****

  elemental function pow_i_ (x, iy) result (pow_x)

    real(dp), intent(in) :: x
    integer, intent(in)  :: iy
    real(dp)             :: pow_x

    pow_x = x**iy
    
  end function pow_i_

  !****

  elemental function pow_r_ (x, y) result (pow_x)

    real(dp), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp)             :: pow_x

    pow_x = x**y
    
  end function pow_r_

  !****

  elemental function cospi_ (x) result (cospi_x)

    real(dp), intent(in) :: x
    real(dp)             :: cospi_x

    cospi_x = COS(x*PI)

  end function cospi_

  !****

  elemental function sinpi_ (x) result (sinpi_x)

    real(dp), intent(in) :: x
    real(dp)             :: sinpi_x

    sinpi_x = SIN(x*PI)

  end function sinpi_

  !****

  elemental function tanpi_ (x) result (tanpi_x)

    real(dp), intent(in) :: x
    real(dp)             :: tanpi_x

    tanpi_x = TAN(x*PI)

  end function tanpi_

  !****

  elemental function acospi_ (x) result (acospi_x)

    real(dp), intent(in) :: x
    real(dp)             :: acospi_x

    acospi_x = ACOS(x)/PI

  end function acospi_

  !****

  elemental function asinpi_ (x) result (asinpi_x)

    real(dp), intent(in) :: x
    real(dp)             :: asinpi_x

    asinpi_x = ASIN(x)/PI

  end function asinpi_

  !****

  elemental function atanpi_ (x) result (atanpi_x)

    real(dp), intent(in) :: x
    real(dp)             :: atanpi_x

    atanpi_x = ATAN(x)/PI

  end function atanpi_

  include 'precompute_zs.inc'


end module math_lib
