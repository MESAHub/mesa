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
!
! ***********************************************************************

module atm_utils

  ! Uses

  use const_def
  use math_lib

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: E2_NPAIRS = 571

  ! Module variables

  real(dp), target, save  :: E2_x(E2_NPAIRS)
  real(dp), save          :: E2_pairs(2*E2_NPAIRS)
  real(dp), target, save  :: E2_f_ary(4*E2_NPAIRS)
  real(dp), pointer, save :: E2_f1(:), E2_f(:,:)
  logical, save           :: have_E2_interpolant = .false.

  ! Access specifiers

  private

  public :: init
  public :: shutdown
  public :: eval_Teff_g
  public :: eval_Paczynski_gradr
  public :: eval_E2

  ! Procedures

contains

  subroutine init(use_cache, ierr)

    use table_atm, only: table_atm_init

    logical, intent(in) :: use_cache
    integer, intent(out) :: ierr

    ierr = 0

    E2_f1 => E2_f_ary
    E2_f(1:4,1:E2_NPAIRS) => E2_f1(1:4*E2_NPAIRS)

    call table_atm_init(use_cache, ierr)

    include 'e2_pairs.dek'

  end subroutine init

  !****

  subroutine shutdown()

    use table_atm, only : table_atm_shutdown

    call table_atm_shutdown()

  end subroutine shutdown

  !****

  subroutine eval_Teff_g(L, R, M, cgrav, Teff, g)

    real(dp), intent(in)  :: L
    real(dp), intent(in)  :: R
    real(dp), intent(in)  :: M
    real(dp), intent(in)  :: cgrav
    real(dp), intent(out) :: Teff
    real(dp), intent(out) :: g

    ! Evaluate the effective temperature and surface gravity

    Teff = pow(L/(4._dp*pi*R*R*boltz_sigma), 0.25_dp)
   
    g = cgrav * M / (R*R)

  end subroutine eval_Teff_g

  !****

  function eval_Paczynski_gradr( &
       T, P, rho, tau, kap, L, M, R, cgrav) result (gradr)

    use eos_lib, only: radiation_pressure
    
    real(dp), intent(in) :: T
    real(dp), intent(in) :: P
    real(dp), intent(in) :: rho
    real(dp), intent(in) :: tau
    real(dp), intent(in) :: kap
    real(dp), intent(in) :: L
    real(dp), intent(in) :: R
    real(dp), intent(in) :: M
    real(dp), intent(in) :: cgrav
    real(dp)             :: gradr
    
    real(dp) :: Prad
    real(dp) :: dilution_factor
    real(dp) :: s
    real(dp) :: f

    ! Evaluate the radiative temperature gradient, using expressions
    ! from Paczynski (1969, Act Ast, 19, 1)

    Prad = radiation_pressure(T)

    gradr = P*kap*L / (16._dp*pi*clight*M*cgrav*Prad)

    if (tau < 2._dp/3._dp) then ! Eqn. 8

       ! Eqn. 15

       s = (2._dp*crad*T*T*T*SQRT(R))/(3._dp*cgrav*M*rho)*pow(L/(8._dp*pi*boltz_sigma), 0.25_dp)

       ! Eqn. 8

       f = 1._dp - 1.5_dp*tau

       dilution_factor = (1._dp + f*s*(4._dp*pi*cgrav*clight*M)/(kap*L))/(1._dp + f*s)

       gradr = gradr*dilution_factor

    end if

    ! Finish

    return

  end function eval_Paczynski_gradr
       
  !****

  subroutine eval_E2(x, E2, dE2_dx, ierr)

    use interp_1d_lib
    use interp_1d_def

    real(dp), intent(in)  :: x
    real(dp), intent(out) :: E2
    real(dp), intent(out) :: dE2_dx
    integer, intent(out)  :: ierr

    real(dp) :: val
    real(dp) :: slope

    ierr = 0

    ! Evaluate the E2 exponential integral

    if (.not. have_E2_interpolant) then
       call create_E2_interpolant(ierr)
       if (ierr /= 0) return
    end if

    call interp_value_and_slope(E2_x, E2_NPAIRS, E2_f1, x, val, slope, ierr)
    if (ierr /= 0) then
       write(*,*) 'call to interp_value_and_slope failed in eval_E2'
       return
    end if

    ! val = log10[E2]
    E2 = exp10(val)
    dE2_dx = slope*ln10*E2

  end subroutine eval_E2
  
  !****
       
  subroutine create_E2_interpolant(ierr)

    use interp_1d_lib
    use interp_1d_def

    use atm_def

    integer, intent(out) :: ierr

    integer, parameter :: NWORK = pm_work_size

    real(dp), target  :: work_ary(E2_NPAIRS*NWORK)
    real(dp), pointer :: work(:)
    integer           :: i

    ierr = 0

    ! Set up an interpolating spline for the E2 exponential integral

    work => work_ary

    do i= 1, E2_NPAIRS
       E2_x(i) = E2_pairs(2*i-1)
       E2_f(1,i) = E2_pairs(2*i)
    end do

    call interp_pm(E2_x, E2_NPAIRS, E2_f1, nwork, work, 'atm_utils', ierr)
    if (ierr /= 0) then
       write(*,*) 'call to interp_pm failed in create_E2_interpolant'
    end if

    have_E2_interpolant = .true.

  end subroutine create_E2_interpolant

end module atm_utils
