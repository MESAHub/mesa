! ***********************************************************************
!
!   Copyright (C) 2013-2021  Josiah Schwab & The MESA Team
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

! In general, this routine follows the notation in Schwab et al. (2015)

module eval_psi

  use const_def, only: dp
  use math_lib

  implicit none

contains

  subroutine fd_integral(dk, y, F)
    ! wrap eos_fermi_dirac_integral for auto_diff
    use auto_diff
    use eos_lib, only: eos_fermi_dirac_integral
    real(dp), intent(in) :: dk
    type(auto_diff_real_2var_order1), intent(in) :: y
    type(auto_diff_real_2var_order1), intent(out) :: F
    real(dp) :: Fval, dF_dy, dF_dz
    call eos_fermi_dirac_integral(dk, y% val, 0.0_dp, Fval, dF_dy, dF_dz)
    F = Fval
    F% d1val1 = dF_dy * y% d1val1
    F% d1val2 = dF_dy * y% d1val2
  end subroutine fd_integral


  subroutine do_psi_Iec_and_Jec(beta, zeta, eta, I, J)

    use auto_diff
    ! calculate the phase space integral for electron emission (beta-decay)


    ! auto_diff variables have
    ! var1: lnT
    ! var2: lnRho

    type(auto_diff_real_2var_order1), intent(in) :: beta  ! mec2 / kT
    type(auto_diff_real_2var_order1), intent(in) :: zeta  ! Q_n / kT
    type(auto_diff_real_2var_order1), intent(in) :: eta   ! chemical potential / kT

    type(auto_diff_real_2var_order1), intent(out) :: I, J   ! phase space integral

    type(auto_diff_real_2var_order1) :: y

    type(auto_diff_real_2var_order1) :: F2, F3, F4, F5
    type(auto_diff_real_2var_order1) :: c2, c3

    ! check that assumptions are met
    if (zeta>-beta) stop "ECAPTURE:  zeta > -beta"

    y = zeta+eta

    call fd_integral(5.0_dp, y, F5)
    call fd_integral(4.0_dp, y, F4)
    call fd_integral(3.0_dp, y, F3)
    call fd_integral(2.0_dp, y, F2)

    c3 = -2.0_dp*zeta
    c2 = zeta*zeta

    I = F4 + c3*F3 + c2*F2
    J = F5 + c3*F4 + c2*F3

    I = I / pow5(beta)
    J = J / pow6(beta)

    return

  end subroutine do_psi_Iec_and_Jec


  subroutine do_psi_Iee_and_Jee(beta, zeta, eta, I, J)

    use auto_diff

    ! calculate the phase space integral for electron emission (beta-decay)


    ! auto_diff variables have
    ! var1: lnT
    ! var2: lnRho

    type(auto_diff_real_2var_order1), intent(in) :: beta  ! mec2 / kT
    type(auto_diff_real_2var_order1), intent(in) :: zeta  ! Q_n / kT
    type(auto_diff_real_2var_order1), intent(in) :: eta   ! chemical potential / kT

    type(auto_diff_real_2var_order1), intent(out) :: I, J   ! phase space integral

    type(auto_diff_real_2var_order1) :: y

    type(auto_diff_real_2var_order1) :: F0, F1, F2, F3, F4, F5
    type(auto_diff_real_2var_order1) :: c0, c1, c2, c3, c4

    ! check that assumptions are met
    if (zeta<beta) stop "ECAPTURE:  zeta < beta"


    y = zeta-eta

    call fd_integral(5.0_dp, y, F5)
    call fd_integral(4.0_dp, y, F4)
    call fd_integral(3.0_dp, y, F3)
    call fd_integral(2.0_dp, y, F2)

    c3 = -2.0_dp*zeta
    c2 = zeta*zeta

    I = F4 + c3*F3 + c2*F2
    J = F5 + c3*F4 + c2*F3


    y = beta-eta

    call fd_integral(5.0_dp, y, F5)
    call fd_integral(4.0_dp, y, F4)
    call fd_integral(3.0_dp, y, F3)
    call fd_integral(2.0_dp, y, F2)
    call fd_integral(1.0_dp, y, F1)
    call fd_integral(0.0_dp, y, F0)

    c3 = (2.0_dp*zeta-4.0_dp*beta)
    c2 = (6.0_dp*beta*beta - 6.0_dp*beta*zeta + zeta*zeta)
    c1 = -2.0_dp*beta*(2.0_dp*Beta*beta - 3.0_dp*beta*zeta + zeta*zeta)
    c0 = beta*beta*(beta-zeta)*(beta-zeta)

    I = I - (F4 + c3*F3 + c2*F2 + c1*F1 + c0*F0)

    c4 = (3.0_dp*zeta-5.0_dp*beta)
    c3 = (10.0_dp*beta*beta-12.0_dp*beta*zeta+3.0_dp*zeta*zeta)
    c2 = (zeta*zeta*zeta - 9.0_dp*beta*zeta*zeta + 18.0_dp*beta*beta*zeta - 10.0_dp*beta*beta*beta)
    c1 = beta*(5.0_dp*beta - 2.0_dp*zeta)*(beta-zeta)*(beta-zeta)
    c0 = -beta*beta*pow3(beta-zeta)

    J = J - (F5 + c4*F4 + c3*F3 + c2*F2 + c1*F1 + c0*F0)


    I = I / pow5(beta)
    J = J / pow6(beta)

    return

  end subroutine do_psi_Iee_and_Jee

end module eval_psi
