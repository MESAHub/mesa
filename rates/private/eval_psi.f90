! ***********************************************************************
!
!   Copyright (C) 2013-2019  Josiah Schwab, Bill Paxton & The MESA Team
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

! In general, this routine follows the notation in Schwab et al. (2013)

module eval_psi

  use const_def, only: dp
  use eos_lib, only: eos_fermi_dirac_integral
  use math_lib

  implicit none

contains

  subroutine do_psi_Iec_and_Jec(beta, zeta, eta, deta_dlnT, deta_dlnRho, &
       I, dI_dlnT, dI_dlnRho, J, dJ_dlnT, dJ_dlnRho)

    ! calulate the phase space integral for electron emission (beta-decay)

    implicit none

    real(dp), intent(in) :: beta  ! mec2 / kT
    real(dp), intent(in) :: zeta  ! Q_n / kT
    real(dp), intent(in) :: eta   ! chemical potential / kT
    real(dp), intent(in) :: deta_dlnT, deta_dlnRho ! and derivs

    real(dp), intent(out) :: I, J   ! phase space integral
    real(dp), intent(out) :: dI_dlnT, dI_dlnRho ! and derivatives
    real(dp), intent(out) :: dJ_dlnT, dJ_dlnRho ! and derivatives

    real(dp) :: betai, eta_e_F, eta_e_L

    real(dp) :: x,y,z
    real(dp) :: dy_dlnT, dy_dlnRho

    real(dp) :: F2, F3, F4, F5
    real(dp) :: dF2_dy, dF3_dy, dF4_dy, dF5_dy
    real(dp) :: dF_dz

    real(dp) :: bm5, bm6
    real(dp) :: c2, c3
    real(dp) :: dc2_dlnT, dc3_dlnT

    ! check that assumptions are met
    if (zeta.gt.-beta) stop "ECAPTURE:  zeta > -beta"

    y = zeta+eta
    z = 0.0_dp

    call eos_fermi_dirac_integral(5.0_dp, y, z, F5, dF5_dy, dF_dz)
    call eos_fermi_dirac_integral(4.0_dp, y, z, F4, dF4_dy, dF_dz)
    call eos_fermi_dirac_integral(3.0_dp, y, z, F3, dF3_dy, dF_dz)
    call eos_fermi_dirac_integral(2.0_dp, y, z, F2, dF2_dy, dF_dz)

    c3 = -2.0_dp*zeta
    c2 = zeta*zeta

    I = F4 + c3*F3 + c2*F2
    J = F5 + c3*F4 + c2*F3

    ! derivatives
    dy_dlnT = -zeta + deta_dlnT
    dy_dlnRho = deta_dlnRho
    
    dc3_dlnT = -c3
    dc2_dlnT = -2.0_dp*c2

    ! dc?_dlnRho = 0

    dI_dlnT = (dF4_dy + c3*dF3_dy + c2*dF2_dy) * dy_dlnT + & 
         (dc3_dlnT * F3 + dc2_dlnT * F2)
    dJ_dlnT = (dF5_dy + c3*dF4_dy + c2*dF3_dy) * dy_dlnT + & 
         (dc3_dlnT * F4 + dc2_dlnT * F3)

    dI_dlnRho = (dF4_dy + c3*dF3_dy + c2*dF2_dy) * dy_dlnRho
    dJ_dlnRho = (dF5_dy + c3*dF4_dy + c2*dF3_dy) * dy_dlnRho

    bm5 = 1d0/pow5(beta)
    bm6 = bm5/beta

    ! put in the powers of beta
    I = bm5 * I
    J = bm6 * J

    dI_dlnT = 5.0_dp*I + bm5 * dI_dlnT
    dJ_dlnT = 6.0_dp*J + bm6 * dJ_dlnT

    dI_dlnRho = bm5 * dI_dlnRho
    dJ_dlnRho = bm6 * dJ_dlnRho

    return

  end subroutine do_psi_Iec_and_Jec


  subroutine do_psi_Iee_and_Jee(beta, zeta, eta, deta_dlnT, deta_dlnRho, &
       I, dI_dlnT, dI_dlnRho, J, dJ_dlnT, dJ_dlnRho)

    ! calulate the phase space integral for electron emission (beta-decay)

    implicit none

    real(dp), intent(in) :: beta  ! mec2 / kT
    real(dp), intent(in) :: zeta  ! Q_n / kT
    real(dp), intent(in) :: eta   ! chemical potential / kT
    real(dp), intent(in) :: deta_dlnT, deta_dlnRho ! and derivs

    real(dp), intent(out) :: I, J   ! phase space integral
    real(dp), intent(out) :: dI_dlnT, dI_dlnRho ! and derivatives
    real(dp), intent(out) :: dJ_dlnT, dJ_dlnRho ! and derivatives

    real(dp) :: betai, eta_e_F, eta_e_L

    real(dp) :: x,y,z
    real(dp) :: dy_dlnT, dy_dlnRho

    real(dp) :: F0, F1, F2, F3, F4, F5
    real(dp) :: dF0_dy, dF1_dy, dF2_dy, dF3_dy, dF4_dy, dF5_dy
    real(dp) :: dF_dz

    real(dp) :: bm5, bm6
    real(dp) :: c0, c1, c2, c3, c4
    real(dp) :: dc0_dlnT, dc1_dlnT, dc2_dlnT, dc3_dlnT, dc4_dlnT

    ! check that assumptions are met
    if (zeta.lt.beta) stop "ECAPTURE:  zeta < beta"

    y = zeta-eta
    z = 0.0_dp

    call eos_fermi_dirac_integral(5.0_dp, y, z, F5, dF5_dy, dF_dz)
    call eos_fermi_dirac_integral(4.0_dp, y, z, F4, dF4_dy, dF_dz)
    call eos_fermi_dirac_integral(3.0_dp, y, z, F3, dF3_dy, dF_dz)
    call eos_fermi_dirac_integral(2.0_dp, y, z, F2, dF2_dy, dF_dz)

    c3 = -2.0_dp*zeta
    c2 = zeta*zeta

    I = F4 + c3*F3 + c2*F2
    J = F5 + c3*F4 + c2*F3

    ! derivatives
    dy_dlnT = -zeta - deta_dlnT
    dy_dlnRho = -deta_dlnRho
    
    dc3_dlnT = -c3
    dc2_dlnT = -2.0_dp*c2

    ! dc?_dlnRho = 0

    dI_dlnT = (dF4_dy + c3*dF3_dy + c2*dF2_dy) * dy_dlnT + & 
         (dc3_dlnT * F3 + dc2_dlnT * F2)
    dJ_dlnT = (dF5_dy + c3*dF4_dy + c2*dF3_dy) * dy_dlnT + & 
         (dc3_dlnT * F4 + dc2_dlnT * F3)

    dI_dlnRho = (dF4_dy + c3*dF3_dy + c2*dF2_dy) * dy_dlnRho
    dJ_dlnRho = (dF5_dy + c3*dF4_dy + c2*dF3_dy) * dy_dlnRho

    ! evalute the fermi-dirac functions
    y = beta-eta
    z = 0.0_dp

    call eos_fermi_dirac_integral(5.0_dp, y, z, F5, dF5_dy, dF_dz)
    call eos_fermi_dirac_integral(4.0_dp, y, z, F4, dF4_dy, dF_dz)
    call eos_fermi_dirac_integral(3.0_dp, y, z, F3, dF3_dy, dF_dz)
    call eos_fermi_dirac_integral(2.0_dp, y, z, F2, dF2_dy, dF_dz)
    call eos_fermi_dirac_integral(1.0_dp, y, z, F1, dF1_dy, dF_dz)
    call eos_fermi_dirac_integral(0.0_dp, y, z, F0, dF0_dy, dF_dz)

    c3 = (2.0_dp*zeta-4.0_dp*beta)
    c2 = (6.0_dp*beta*beta - 6.0_dp*beta*zeta + zeta*zeta)
    c1 = -2.0_dp*beta*(2.0_dp*Beta*beta - 3.0_dp*beta*zeta + zeta*zeta)
    c0 = beta*beta*(beta-zeta)*(beta-zeta)

    I = I - (F4 + c3*F3 + c2*F2 + c1*F1 + c0*F0)

    ! derivatives
    dy_dlnT = -beta - deta_dlnT
    dy_dlnRho = deta_dlnRho

    dc3_dlnT = -c3
    dc2_dlnT = -2.0_dp*c2
    dc1_dlnT = -3.0_dp*c1
    dc0_dlnT = -4.0_dp*c0

    !dc?_dlnRho = 0

    dI_dlnT = dI_dlnT - (&
         (dF4_dy + c3*dF3_dy + c2*dF2_dy + c1*dF1_dy + c0*dF0_dy) * dy_dlnT + & 
         (dc3_dlnT * F3 + dc2_dlnT * F2 + dc1_dlnT * F1 + dc0_dlnT * F0))

    dI_dlnRho = dI_dlnRho - &
         (dF4_dy + c3*dF3_dy + c2*dF2_dy + c1*dF1_dy + c0*dF0_dy) * dy_dlnRho

    c4 = (3.0_dp*zeta-5.0_dp*beta)
    c3 = (10.0_dp*beta*beta-12.0_dp*beta*zeta+3.0_dp*zeta*zeta)
    c2 = (zeta*zeta*zeta - 9.0_dp*beta*zeta*zeta + 18.0_dp*beta*beta*zeta - 10.0_dp*beta*beta*beta)
    c1 = beta*(5.0_dp*beta - 2.0_dp*zeta)*(beta-zeta)*(beta-zeta)
    c0 = -beta*beta*pow3(beta-zeta)

    dc4_dlnT = -c4
    dc3_dlnT = -2.0_dp*c3
    dc2_dlnT = -3.0_dp*c2
    dc1_dlnT = -4.0_dp*c1
    dc0_dlnT = -5.0_dp*c0

    !dc?_dlnRho = 0

    J = J - (F5 + c4*F4 + c3*F3 + c2*F2 + c1*F1 + c0*F0)

    dJ_dlnT = dJ_dlnT - (&
         (dF5_dy + c4*dF4_dy + c3*dF3_dy + c2*dF2_dy + c1*dF1_dy + c0*dF0_dy) * dy_dlnT + & 
         (dc4_dlnT * F4 + dc3_dlnT * F3 + dc2_dlnT * F2 + dc1_dlnT * F1 + dc0_dlnT * F0))

    dJ_dlnRho = dJ_dlnRho - &
         (dF5_dy + c4*dF4_dy + c3*dF3_dy + c2*dF2_dy + c1*dF1_dy + c0*dF0_dy) * dy_dlnRho

    ! precompute needed powers of beta
    bm5 = 1d0/pow5(beta)
    bm6 = bm5/beta

    ! put in the powers of beta
    I = bm5 * I
    J = bm6 * J

    dI_dlnT = 5.0_dp*I + bm5 * dI_dlnT
    dJ_dlnT = 6.0_dp*J + bm6 * dJ_dlnT

    dI_dlnRho = bm5 * dI_dlnRho
    dJ_dlnRho = bm6 * dJ_dlnRho

    return

  end subroutine do_psi_Iee_and_Jee

end module eval_psi
