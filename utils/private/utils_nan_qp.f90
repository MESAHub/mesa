! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

module utils_nan_qp

  ! Uses

  use const_def, only : qp

  use ISO_FORTRAN_ENV
  use ISO_C_BINDING

  ! No implicit typing
      
  implicit none

  ! Parameters

  integer, parameter :: FRAC_BITS_128_H = 48
  integer, parameter :: FRAC_BITS_128_L = 64
  integer, parameter :: EXPN_BITS_128_H = 15

  integer(INT64), parameter :: QNAN_128_H = INT(z'7fff100000000000', INT64)
  integer(INT64), parameter :: QNAN_128_L = INT(z'0000000000000001', INT64)
  integer(INT64), parameter :: SNAN_128_H = INT(z'7fff000000000000', INT64)
  integer(INT64), parameter :: SNAN_128_L = INT(z'0000000000000001', INT64)

  ! Interfaces

  interface is_nan
     module procedure is_nan_qp
  end interface is_nan

  interface is_inf
     module procedure is_inf_qp
  end interface is_inf

  interface is_bad
     module procedure is_bad_qp
  end interface is_bad

  interface set_nan
     module procedure set_nan_qp_0d
     module procedure set_nan_qp_1d
     module procedure set_nan_qp_2d
     module procedure set_nan_qp_3d
     module procedure set_nan_qp_4d
  end interface set_nan

  ! Access specifiers

  private

  public :: is_nan
  public :: is_inf
  public :: is_bad
  public :: set_nan

  ! Procedures
      
contains

  elemental function is_nan_qp (x, signal) result (is_nan)

    real(qp), target, intent(in)  :: x
    logical, optional, intent(in) :: signal
    logical                       :: is_nan

    integer(INT64) :: ix(2)
    integer(INT64) :: frac_l
    integer(INT64) :: frac_h
    integer(INT64) :: expn
    integer(INT64) :: sign

    ! Convert x to integer

    ix = TRANSFER(x, ix)

    ! Split out IEEE fields
    
    frac_l = IBITS(ix(1), 0, FRAC_BITS_128_L)
    frac_h = IBITS(ix(2), 0, FRAC_BITS_128_H)
    expn = IBITS(ix(2), FRAC_BITS_128_H, EXPN_BITS_128_H)
    sign = IBITS(ix(2), FRAC_BITS_128_H+EXPN_BITS_128_H, 1)

    ! Check for NaN

    is_nan = expn == MASKR(EXPN_BITS_128_H, INT64) .AND. &
             (frac_h /= 0_INT64 .OR. frac_l /= 0_INT64)

    if (PRESENT(signal)) then
       is_nan = is_nan .AND. (BTEST(frac_h, FRAC_BITS_128_H-1) .EQV. .NOT. signal)
    endif

    ! Finish

    return

  end function is_nan_qp

  !****

  elemental function is_inf_qp (x) result (is_inf)

    real(qp), target, intent(in)  :: x
    logical                       :: is_inf

    integer(INT64) :: ix(2)
    integer(INT64) :: frac_l
    integer(INT64) :: frac_h
    integer(INT64) :: expn
    integer(INT64) :: sign

    ! Convert x to integer

    ix = TRANSFER(x, ix)

    frac_l = IBITS(ix(1), 0, FRAC_BITS_128_L)
    frac_h = IBITS(ix(2), 0, FRAC_BITS_128_H)
    expn = IBITS(ix(2), FRAC_BITS_128_H, EXPN_BITS_128_H)
    sign = IBITS(ix(2), FRAC_BITS_128_H+EXPN_BITS_128_H, 1)

    ! Check for infinity

    is_inf = expn == MASKR(EXPN_BITS_128_H, INT64) .AND. &
             (frac_l == 0_INT64 .AND. frac_h == 0_INT64)

    ! Finish

    return

  end function is_inf_qp

  !****

  elemental function is_bad_qp (x) result (is_bad)

    real(qp), intent(in) :: x
    logical              :: is_bad

    ! Check for NaN or infinity

    is_bad = is_nan(x) .OR. is_inf(x)

    ! Finish

    return

  end function is_bad_qp

  !****
      
  subroutine set_nan_qp_0d (x, signal)

    real(qp), target, intent(out) :: x
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix(:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, [2])

    if (PRESENT(signal)) then
       if (signal) then
          ix(1) = SNAN_128_L
          ix(2) = SNAN_128_H
       else
          ix(1) = QNAN_128_L
          ix(2) = QNAN_128_H
       endif
    else
       ix(1) = SNAN_128_L
       ix(2) = SNAN_128_H
    endif

    ! Finish

    return

  end subroutine set_nan_qp_0d

  !****

  subroutine set_nan_qp_1d (x, signal)

    real(qp), target, intent(out) :: x(:)
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix(:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, [2,SHAPE(x)])

    if (PRESENT(signal)) then
       if (signal) then
          ix(1,:) = SNAN_128_L
          ix(2,:) = SNAN_128_H
       else
          ix(1,:) = QNAN_128_L
          ix(2,:) = QNAN_128_H
       endif
    else
       ix(1,:) = SNAN_128_L
       ix(2,:) = SNAN_128_H
    endif

    ! Finish

    return

  end subroutine set_nan_qp_1d

  !****

  subroutine set_nan_qp_2d (x, signal)

    real(qp), target, intent(out) :: x(:,:)
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix(:,:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, [2,SHAPE(x)])

    if (PRESENT(signal)) then
       if (signal) then
          ix(1,:,:) = SNAN_128_L
          ix(2,:,:) = SNAN_128_H
       else
          ix(1,:,:) = QNAN_128_L
          ix(2,:,:) = QNAN_128_H
       endif
    else
       ix(1,:,:) = SNAN_128_L
       ix(2,:,:) = SNAN_128_H
    endif

    ! Finish

    return

  end subroutine set_nan_qp_2d

  !****

  subroutine set_nan_qp_3d (x, signal)

    real(qp), target, intent(out) :: x(:,:,:)
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix(:,:,:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, [2,SHAPE(x)])

    if (PRESENT(signal)) then
       if (signal) then
          ix(1,:,:,:) = SNAN_128_L
          ix(2,:,:,:) = SNAN_128_H
       else
          ix(1,:,:,:) = QNAN_128_L
          ix(2,:,:,:) = QNAN_128_H
       endif
    else
       ix(1,:,:,:) = SNAN_128_L
       ix(2,:,:,:) = SNAN_128_H
    endif

    ! Finish

    return

  end subroutine set_nan_qp_3d

  !****

  subroutine set_nan_qp_4d (x, signal)

    real(qp), target, intent(out) :: x(:,:,:,:)
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix(:,:,:,:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, [2,SHAPE(x)])

    if (PRESENT(signal)) then
       if (signal) then
          ix(1,:,:,:,:) = SNAN_128_L
          ix(2,:,:,:,:) = SNAN_128_H
       else
          ix(1,:,:,:,:) = QNAN_128_L
          ix(2,:,:,:,:) = QNAN_128_H
       endif
    else
       ix(1,:,:,:,:) = SNAN_128_L
       ix(2,:,:,:,:) = SNAN_128_H
    endif

    ! Finish

    return

  end subroutine set_nan_qp_4d

end module utils_nan_qp
