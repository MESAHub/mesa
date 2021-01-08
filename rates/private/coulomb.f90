! ***********************************************************************
!
!   Copyright (C) 2014-2019  Josiah Schwab, Bill Paxton & The MESA Team
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

  implicit none

contains

  subroutine do_coulomb_set_context( &
       cc, temp, den, logT, logRho, zbar, abar, z2bar,  &
       num_isos, y, iso_z52)
    type (Coulomb_Info), pointer :: cc
    integer, intent(in) :: num_isos
    real(dp), intent(in) ::  &
         temp, den, logT, logRho, zbar, abar, z2bar, &
         y(:), iso_z52(:) ! Z to the power of 5/2

    real(dp), parameter :: x13   = 1.0d0/3.0d0 
    real(dp), parameter :: x52   = 5.0d0/2.0d0
    real(dp) :: qq
    integer :: j
    
    include 'formats'

    cc% temp  = temp
    cc% den   = den
    cc% logT  = logT
    cc% logRho = logRho
    cc% zbar  = zbar
    cc% abar  = abar
    cc% z2bar = z2bar

    ! get the info that depends only on temp, den, and overall composition         

    cc% abari    = 1.0d0/abar
    cc% rr       = den * cc% abari
    cc% tempi    = 1.0d0/temp
    cc% dtempi   = -cc% tempi * cc% tempi
    cc% deni     = 1.0d0/den
    cc% ye       = cc% zbar * cc% abari

    ! calculate the other powers <z^?> that we need
    qq = 0d0
    do j=1,num_isos
      qq = qq + iso_z52(j) * y(j)
      if (is_bad(qq)) then
         write(*,2) 'qq', j, qq, iso_z52(j), x52, y(j)
         stop 'do_coulomb_set_context'
      end if
    end do
    cc% z52bar = abar * qq

    ! calculate powers of zbar
    cc% zbar13   = pow(zbar,x13)

    ! calculate more complicated, but useful quantities
    cc% gamma_e = 2.275d5 * pow(cc% ye * cc% den, x13) * cc% tempi
    cc% rs = 1.388_dp * pow(cc% ye * cc% den, -x13)

  end subroutine do_coulomb_set_context

end module coulomb

