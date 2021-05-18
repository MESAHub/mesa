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

module atm_lib

  ! Uses
  
  use const_def, only: dp

  use atm_utils, only: &
       atm_init => init, &
       atm_shutdown => shutdown

  use atm_T_tau_uniform, only: &
       atm_eval_T_tau_uniform => eval_T_tau_uniform, &
       atm_build_T_tau_uniform => build_T_tau_uniform
  use atm_T_tau_varying, only: &
       atm_eval_T_tau_varying => eval_T_tau_varying, &
       atm_build_T_tau_varying => build_T_tau_varying
  use atm_T_tau_relations, only: &
       atm_get_T_tau_base => get_T_tau_base, &
       atm_eval_T_tau_dq_dtau => eval_T_tau_dq_dtau

  use atm_table, only: &
       atm_eval_table => eval_table, &
       atm_get_table_alfa_beta => get_table_alfa_beta, &
       atm_get_table_base => get_table_base

  use atm_irradiated, only: &
       atm_eval_irradiated => eval_irradiated
  
  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: atm_init
  public :: atm_shutdown

  public :: atm_eval_T_tau_uniform
  public :: atm_build_T_tau_uniform
  public :: atm_eval_T_tau_varying
  public :: atm_build_T_tau_varying
  public :: atm_get_T_tau_base
  public :: atm_eval_T_tau_dq_dtau

  public :: atm_eval_table
  public :: atm_get_table_alfa_beta
  public :: atm_get_table_base

  public :: atm_eval_irradiated

  public :: atm_Teff
  public :: atm_L
  public :: atm_black_body_T

contains

  ! Crufty utility routines

  real(dp) function atm_Teff(L, R)
    use const_def, only: pi, boltz_sigma
    real(dp), intent(in) :: L, R
    atm_Teff = atm_black_body_T(L, R)
  end function atm_Teff

  real(dp) function atm_L(Teff, R)
    use const_def, only: pi, boltz_sigma
    real(dp), intent(in) :: Teff, R
    atm_L = 4d0*pi*R*R*boltz_sigma*Teff*Teff*Teff*Teff
  end function atm_L

  real(dp) function atm_black_body_T(L, R)
    use math_lib, only: pow
    use const_def, only: pi, boltz_sigma
    real(dp), intent(in) :: L, R
    atm_black_body_T = pow(L / (4d0*pi*R*R*boltz_sigma), 0.25d0)
  end function atm_black_body_T
 
end module atm_lib

