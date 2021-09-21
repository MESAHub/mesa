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

module utils_nan_dp

  ! Uses

  use const_def, only : dp

  use ISO_FORTRAN_ENV
  use ISO_C_BINDING

  ! No implicit typing
      
  implicit none

  ! Parameters

  integer, parameter :: FRAC_BITS_64 = 52
  integer, parameter :: EXPN_BITS_64 = 11

  integer(INT64), parameter :: QNAN_64 = INT(z'7ff8000000000001', INT64)
  integer(INT64), parameter :: SNAN_64 = INT(z'7ff0000000000001', INT64)

  ! Interfaces

  interface is_nan
     module procedure is_nan_dp
  end interface is_nan

  interface is_inf
     module procedure is_inf_dp
  end interface is_inf

  interface is_bad
     module procedure is_bad_dp
  end interface is_bad

  interface set_nan
     module procedure set_nan_dp_0d
     module procedure set_nan_dp_1d
     module procedure set_nan_dp_2d
     module procedure set_nan_dp_3d
     module procedure set_nan_dp_4d
  end interface set_nan

  ! Access specifiers

  private

  public :: is_nan
  public :: is_inf
  public :: is_bad
  public :: set_nan

  ! Procedures
      
contains

  elemental function is_nan_dp (x, signal) result (is_nan)

    real(dp), target, intent(in)  :: x
    logical, optional, intent(in) :: signal
    logical                       :: is_nan

    integer(INT64) :: ix
    integer(INT64) :: frac
    integer(INT64) :: expn
    integer(INT64) :: sign

    ! Convert x to integer

    ix = TRANSFER(x, ix)

    ! Split out IEEE fields
    
    frac = IBITS(ix, 0, FRAC_BITS_64)
    expn = IBITS(ix, FRAC_BITS_64, EXPN_BITS_64)
    sign = IBITS(ix, FRAC_BITS_64+EXPN_BITS_64, 1)

    ! Check for NaN

    is_nan = expn == MASKR(EXPN_BITS_64, INT64) .AND. frac /= 0_INT64

    if (PRESENT(signal)) then
       is_nan = is_nan .AND. (BTEST(frac, FRAC_BITS_64-1) .EQV. .NOT. signal)
    endif

    ! Finish

    return

  end function is_nan_dp

  !****

  elemental function is_inf_dp (x) result (is_inf)

    real(dp), target, intent(in)  :: x
    logical                       :: is_inf

    integer(INT64) :: ix
    integer(INT64) :: frac
    integer(INT64) :: expn
    integer(INT64) :: sign

    ! Convert x to integer

    ix = TRANSFER(x, ix)

    ! Split out IEEE fields
    
    frac = IBITS(ix, 0, FRAC_BITS_64)
    expn = IBITS(ix, FRAC_BITS_64, EXPN_BITS_64)
    sign = IBITS(ix, FRAC_BITS_64+EXPN_BITS_64, 1)

    ! Check for infinity

    is_inf = expn == MASKR(EXPN_BITS_64, INT64) .AND. frac == 0_INT64

    ! Finish

    return

  end function is_inf_dp

  !****

  elemental function is_bad_dp (x) result (is_bad)

    real(dp), intent(in) :: x
    logical              :: is_bad

    ! Check for NaN or infinity

    is_bad = is_nan(x) .OR. is_inf(x)

    ! Finish

    return

  end function is_bad_dp

  !****
      
  subroutine set_nan_dp_0d (x, signal)

    real(dp), target, intent(out) :: x
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix)

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_64
       else
          ix = QNAN_64
       endif
    else
       ix = SNAN_64
    endif

    ! Finish

    return

  end subroutine set_nan_dp_0d

  !****

  subroutine set_nan_dp_1d (x, signal)

    real(dp), target, intent(out) :: x(:)
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix(:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, SHAPE(x))

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_64
       else
          ix = QNAN_64
       endif
    else
       ix = SNAN_64
    endif

    ! Finish

    return

  end subroutine set_nan_dp_1d

  !****

  subroutine set_nan_dp_2d (x, signal)

    real(dp), target, intent(out) :: x(:,:)
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix(:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, SHAPE(x))

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_64
       else
          ix = QNAN_64
       endif
    else
       ix = SNAN_64
    endif

    ! Finish

    return

  end subroutine set_nan_dp_2d

  !****

  subroutine set_nan_dp_3d (x, signal)

    real(dp), target, intent(out) :: x(:,:,:)
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix(:,:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, SHAPE(x))

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_64
       else
          ix = QNAN_64
       endif
    else
       ix = SNAN_64
    endif

    ! Finish

    return

  end subroutine set_nan_dp_3d

  !****

  subroutine set_nan_dp_4d (x, signal)

    real(dp), target, intent(out) :: x(:,:,:,:)
    logical, optional, intent(in) :: signal

    integer(INT64), pointer :: ix(:,:,:,:)

    ! Convert x to a signaling or quiet NaN

    call C_F_POINTER(C_LOC(x), ix, SHAPE(x))

    if (PRESENT(signal)) then
       if (signal) then
          ix = SNAN_64
       else
          ix = QNAN_64
       endif
    else
       ix = SNAN_64
    endif

    ! Finish

    return

  end subroutine set_nan_dp_4d

end module utils_nan_dp
