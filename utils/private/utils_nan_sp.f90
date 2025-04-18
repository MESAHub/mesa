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

module utils_nan_sp

  use const_def, only : sp

  use ISO_FORTRAN_ENV
  use ISO_C_BINDING

  implicit none

  integer, parameter :: FRAC_BITS_32 = 23
  integer, parameter :: EXPN_BITS_32 = 8

  integer(INT32), parameter :: QNAN_32 = INT(z'7fc00001', INT32)
  integer(INT32), parameter :: SNAN_32 = INT(z'7f800001', INT32)

  interface is_nan
     module procedure is_nan_sp
  end interface is_nan

  interface is_inf
     module procedure is_inf_sp
  end interface is_inf

  interface is_bad
     module procedure is_bad_sp
  end interface is_bad

  interface set_nan
     module procedure set_nan_sp_0d
     module procedure set_nan_sp_1d
     module procedure set_nan_sp_2d
     module procedure set_nan_sp_3d
     module procedure set_nan_sp_4d
  end interface set_nan

  private

  public :: is_nan
  public :: is_inf
  public :: is_bad
  public :: set_nan

contains

  elemental function is_nan_sp (x, signal) result (is_nan)

    real(sp), target, intent(in)  :: x
    logical, optional, intent(in) :: signal
    logical                       :: is_nan

    integer(INT32) :: ix
    integer(INT32) :: frac
    integer(INT32) :: expn
    integer(INT32) :: sign

    ! Convert x to integer

    ix = TRANSFER(x, ix)

    ! Split out IEEE fields

    frac = IBITS(ix, 0, FRAC_BITS_32)
    expn = IBITS(ix, FRAC_BITS_32, EXPN_BITS_32)
    sign = IBITS(ix, FRAC_BITS_32+EXPN_BITS_32, 1)

    ! Check for NaN

    is_nan = expn == MASKR(EXPN_BITS_32, INT32) .AND. frac /= 0_INT32

    if (PRESENT(signal)) then
       is_nan = is_nan .AND. (BTEST(frac, FRAC_BITS_32-1) .EQV. .NOT. signal)
    end if


    return

  end function is_nan_sp


  elemental function is_inf_sp (x) result (is_inf)

    real(sp), target, intent(in)  :: x
    logical                       :: is_inf

    integer(INT32) :: ix
    integer(INT32) :: frac
    integer(INT32) :: expn
    integer(INT32) :: sign

    ! Convert x to integer

    ix = TRANSFER(x, ix)

    ! Split out IEEE fields

    frac = IBITS(ix, 0, FRAC_BITS_32)
    expn = IBITS(ix, FRAC_BITS_32, EXPN_BITS_32)
    sign = IBITS(ix, FRAC_BITS_32+EXPN_BITS_32, 1)

    ! Check for infinity

    is_inf = expn == MASKR(EXPN_BITS_32, INT32) .AND. frac == 0_INT32


    return

  end function is_inf_sp


  elemental function is_bad_sp (x) result (is_bad)

    real(sp), intent(in) :: x
    logical              :: is_bad

    ! Check for NaN or infinity

    is_bad = is_nan(x) .OR. is_inf(x)


    return

  end function is_bad_sp


  subroutine set_nan_sp_0d (x, signal)

    real(sp), target, intent(out) :: x
    logical, optional, intent(in) :: signal

    integer(INT32), pointer :: ix

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix)

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_32
       else
          ix = QNAN_32
       end if
    else
       ix = SNAN_32
    end if

  end subroutine set_nan_sp_0d


  subroutine set_nan_sp_1d (x, signal)

    real(sp), target, intent(out) :: x(:)
    logical, optional, intent(in) :: signal

    integer(INT32), pointer :: ix(:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, SHAPE(x))

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_32
       else
          ix = QNAN_32
       end if
    else
       ix = SNAN_32
    end if


    return

  end subroutine set_nan_sp_1d


  subroutine set_nan_sp_2d (x, signal)

    real(sp), target, intent(out) :: x(:,:)
    logical, optional, intent(in) :: signal

    integer(INT32), pointer :: ix(:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, SHAPE(x))

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_32
       else
          ix = QNAN_32
       end if
    else
       ix = SNAN_32
    end if


    return

  end subroutine set_nan_sp_2d


  subroutine set_nan_sp_3d (x, signal)

    real(sp), target, intent(out) :: x(:,:,:)
    logical, optional, intent(in) :: signal

    integer(INT32), pointer :: ix(:,:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, SHAPE(x))

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_32
       else
          ix = QNAN_32
       end if
    else
       ix = SNAN_32
    end if


    return

  end subroutine set_nan_sp_3d


  subroutine set_nan_sp_4d (x, signal)

    real(sp), target, intent(out) :: x(:,:,:,:)
    logical, optional, intent(in) :: signal

    integer(INT32), pointer :: ix(:,:,:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, SHAPE(x))

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_32
       else
          ix = QNAN_32
       end if
    else
       ix = SNAN_32
    end if


    return

  end subroutine set_nan_sp_4d

end module utils_nan_sp
