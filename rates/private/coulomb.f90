! ***********************************************************************
!
!   Copyright (C) 2014-2021  Josiah Schwab & The MESA Team
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

module coulomb

  use const_def
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

    ! calculate key plasma paramters
    cc% gamma_e = 2.275d5 * pow(ye * den, one_third) / temp
    cc% rs = 1.388_dp * pow(ye * den, -one_third)

  end subroutine do_coulomb_set_context

end module coulomb

