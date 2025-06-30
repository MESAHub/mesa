! ***********************************************************************
!
!   Copyright (C) 2014-2021  Josiah Schwab & The MESA Team
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

module coulomb

  use const_def, only: dp, one_third
  use rates_def
  use math_lib
  use utils_lib, only: is_bad
  use auto_diff

  implicit none

contains

  subroutine do_coulomb_set_context( &
       cc, temp_in, den_in, logT_in, logRho_in, zbar, abar, z2bar)
    type (Coulomb_Info), pointer :: cc
    real(dp), intent(in) ::  &
       temp_in, den_in, logT_in, logRho_in, zbar, abar, z2bar
    real(dp) :: ye
    type(auto_diff_real_2var_order1) :: temp, den

    include 'formats'

    ! auto_diff variables have
    ! var1: lnT
    ! var2: lnRho

    temp = temp_in
    temp% d1val1 = temp_in
    temp% d1val2 = 0d0

    den = den_in
    den% d1val1 = 0d0
    den% d1val2 = den_in

    cc% temp  = temp_in
    cc% den   = den_in
    cc% logT  = logT_in
    cc% logRho = logRho_in
    cc% zbar  = zbar
    cc% abar  = abar
    cc% z2bar = z2bar

    ye = zbar / abar

    ! calculate key plasma parameters
    cc% gamma_e = 2.275d5 * pow(ye * den, one_third) / temp
    cc% rs = 1.388_dp * pow(ye * den, -one_third)

  end subroutine do_coulomb_set_context

end module coulomb

