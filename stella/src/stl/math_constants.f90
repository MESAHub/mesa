!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/mathconstants [1.0] *
!!
!!   NAME
!!     math_constants
!!
!!   FUNCTION
!!     Definition of mathematical constants and functions.
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     Adapted for use in CP2K (JGH)
!!     JGH (16-06-2002) : Added Gamma functions
!!     JGH (10-08-2004) : Added Euler constant (gamma)
!!     Baklanov (13-07-2006) : change module name 
!!   SOURCE
!******************************************************************************

module math_constants

  use kinds,                           only: dp

  implicit none

  private
  public ::  p_pi, p_pio2, p_twopi, p_fourpi, p_root2, p_rootpi, p_oorootpi, &
       p_zero, p_one, p_half, p_degree, p_radians
  public :: deq,dgt,dne,dle,dlt,dge
  public :: dNeqZero, dEqZero, dLtZero, dGtZero
 
  real(kind=dp), public :: one=1.0_dp,  eps=2*epsilon(one)

  real (kind=dp), parameter :: p_pi = 3.14159265358979323846264338_dp
  real (kind=dp), parameter :: p_pio2 = 1.57079632679489661923132169_dp
  real (kind=dp), parameter :: p_twopi = 6.28318530717958647692528677_dp
  real (kind=dp), parameter :: p_fourpi = 12.56637061435917295385057353_dp
  real (kind=dp), parameter :: p_root2 = 1.41421356237309504880168872_dp
  real (kind=dp), parameter :: p_rootpi = 1.77245387556702677717522568_dp
  real (kind=dp), parameter :: p_oorootpi = 1.0_dp/p_rootpi

  real (kind=dp), parameter :: p_euler=0.577215664901532860606512_dp

  real (kind=dp), parameter :: p_zero = 0.0_dp, p_one = 1.0_dp, p_half = 0.5_dp
  real (kind=dp), parameter :: p_degree = 180.0_dp/p_pi, p_radians = p_one/p_degree

contains

   function deq(x,y) result (test)
    real(kind=dp), intent(in) :: x,y
    logical :: test
    test = abs(x-y) <= max(abs(x),abs(y))*eps
  end function deq

   function dgt(x,y) result (test)
    real(kind=dp), intent(in) :: x,y
    logical :: test
    test = (x-y) > max(abs(x),abs(y))*eps
  end function dgt

    function dne(x,y) result (test)
    real(kind=dp), intent(in) :: x,y
    logical :: test
    test = .not. deq(x,y)
  end function dne

    function dle(x,y) result (test)
    real(kind=dp), intent(in) :: x,y
    logical :: test
    test = .not. dgt(x,y)
  end function dle

    function dlt(x,y) result (test)
    real(kind=dp), intent(in) :: x,y
    logical :: test
    test = dle(x,y) .and. dne(x,y)
  end function dlt

    function dge(x,y) result (test)
    real(kind=dp), intent(in) :: x,y
    logical :: test
    test = dgt(x,y) .or. deq(x,y)
  end function dge
  
! compare with 0  
   logical function dGtZero(x)
    real(kind=dp), intent(in) :: x
    dGtZero =  dgt(x, 0._dp)
  end function dGtZero
 
   logical function dLtZero(x)
    real(kind=dp), intent(in) :: x
    dLtZero =  dlt(x, 0._dp)
  end function dLtZero
 
   logical function dEqZero(x)
    real(kind=dp), intent(in) :: x
write(*,*) ' eps=',eps
    dEqZero =  deq(x, 0._dp)
  end function dEqZero
  
   logical function dNeqZero(x)
    real(kind=dp), intent(in) :: x
    dNeqZero =  dne(x, 0._dp)
  end function dNeqZero
  
  
end module math_constants

!******************************************************************************
