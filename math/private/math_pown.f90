! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

module math_pown

  use const_def, only: dp

  implicit none

  interface powm1
     module procedure powm1_
  end interface powm1

  interface pow2
     module procedure pow2_
  end interface pow2

  interface pow3
     module procedure pow3_
  end interface pow3

  interface pow4
     module procedure pow4_
  end interface pow4

  interface pow5
     module procedure pow5_
  end interface pow5

  interface pow6
     module procedure pow6_
  end interface pow6

  interface pow7
     module procedure pow7_
  end interface pow7

  interface pow8
     module procedure pow8_
  end interface pow8

  private

  public :: powm1, pow2, pow3, pow4, pow5, pow6, pow7, pow8

contains

  elemental function powm1_ (x) result (powm1_x)

    real(dp), intent(in) :: x
    real(dp)             :: powm1_x

    powm1_x = 1.0_dp / x

  end function powm1_


  elemental function pow2_ (x) result (pow2_x)

    real(dp), intent(in) :: x
    real(dp)             :: pow2_x

    pow2_x = x*x

  end function pow2_


  elemental function pow3_ (x) result (pow3_x)

    real(dp), intent(in) :: x
    real(dp)             :: pow3_x

    pow3_x = x*x*x

  end function pow3_


  elemental function pow4_ (x) result (pow4_x)

    real(dp), intent(in) :: x
    real(dp)             :: pow4_x

    pow4_x = x*x*x*x

  end function pow4_


  elemental function pow5_ (x) result (pow5_x)

    real(dp), intent(in) :: x
    real(dp)             :: pow5_x

    pow5_x = x*x*x*x*x

  end function pow5_


  elemental function pow6_ (x) result (pow6_x)

    real(dp), intent(in) :: x
    real(dp)             :: pow6_x

    pow6_x = x*x*x*x*x*x

  end function pow6_


  elemental function pow7_ (x) result (pow7_x)

    real(dp), intent(in) :: x
    real(dp)             :: pow7_x

    pow7_x = x*x*x*x*x*x*x

  end function pow7_


  elemental function pow8_ (x) result (pow8_x)

    real(dp), intent(in) :: x
    real(dp)             :: pow8_x

    pow8_x = x*x*x*x*x*x*x*x

  end function pow8_

end module math_pown
