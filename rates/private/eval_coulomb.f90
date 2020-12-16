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

module eval_coulomb

  use const_def, only: dp
  use math_lib

  implicit none

  integer, parameter :: None = 0

  ! select the option for the ion chemical potential
  integer, parameter :: DGC1973 = 1
  integer, parameter :: I1993 = 2
  integer, parameter :: PCR2009 = 3

  ! select the option for the electron screening
  integer, parameter :: ThomasFermi = 1
  integer, parameter :: Itoh2002 = 2

contains

  function do_mui_coulomb(cc, Z) result(mu)

    use rates_def, only: Coulomb_Info, which_mui_coulomb

    implicit none

    type(Coulomb_Info), intent(in) :: cc
    real(dp), intent(in) :: Z ! nuclear charge

    real(dp) :: mu ! the chemical potential in units of kT

    select case (which_mui_coulomb)
    case (None)
       mu = 0
    case (DGC1973)
       mu = mui_coulomb_DGC1973(cc, Z)
    case (I1993)
       mu = mui_coulomb_I1993(cc, Z)
    case (PCR2009)
       mu = mui_coulomb_PCR2009(cc, Z)
    case default
       stop "STOP (eval_coulomb.f: do_mui_coulomb): incorrect value of which_mui_coulomb"
    end select

    return

  end function do_mui_coulomb


  function mui_coulomb_DGC1973(cc, Z) result(mu)
    use math_lib
    use rates_def, only: Coulomb_Info

    implicit none

    type(Coulomb_Info), intent(in) :: cc
    real(dp), intent(in) :: Z ! nuclear charge

    real(dp) :: mu ! the chemical potential in units of kT

    ! calulate the chemical potential of an ion

    real(dp) :: gamma
    real(dp) :: zr, zr_m1o3

    ! define parameters
    real(dp), parameter :: c0 = 0.9d0
    real(dp), parameter :: c1 = 0.2843d0
    real(dp), parameter :: c2 = -0.054d0
    real(dp), parameter :: d0 = -9d0/ 16d0
    real(dp), parameter :: d1 = 0.460d0

    logical :: debug = .false.

    ! from Dewitt, H. E., Graboske, H. C., & Cooper, M. S. 1973, ApJ, 181, 439
    ! via Couch & Loumos ApJ, vol. 194, Dec. 1, 1974, pt. 1, p. 385-392

    ! coulomb coupling parameter
    gamma = pow(Z,2d0/3d0) * cc% zbar * cc% gamma_e

    ! it appears that Gutierrez et al. use something like
    ! despite citing the above references
    ! gamma = pow(cc% zbar,5d0/3d0) * cc% gamma_e

    ! ratios for convenience
    zr = Z/cc% zbar
    zr_m1o3 = pow(zr,-1d0/3d0) ! zr to the minus 1/3 power

    ! evaluate mu(Z); C&L eq. (11)
    mu = -zr * (gamma * (c0 + (c1  + c2 * zr_m1o3) * zr_m1o3) + &
         (d0 + d1 * zr_m1o3))

    return

  end function mui_coulomb_DGC1973


  function mui_coulomb_I1993(cc, Z) result(mu)
    use math_lib
    use rates_def, only: Coulomb_Info

    implicit none

    type(Coulomb_Info), intent(in) :: cc
    real(dp), intent(in) :: Z ! nuclear charge

    real(dp) :: mu ! the chemical potential in units of kT

    ! calulate the chemical potential of an ion

    ! form from W.L. Slattery, G.D. Doolen, H.E. DeWitt, Phys. Rev. A 26 (1982) 2255.
    ! values from Ichimaru, S. 1993, Reviews of Modern Physics, 65, 255
    ! via Juodagalvis et al. (2010)

    real(dp) :: gamma, gamma14

    ! define parameters
    real(dp), parameter :: a = -0.898004d0
    real(dp), parameter :: b = 0.96786d0
    real(dp), parameter :: c = 0.220703d0
    real(dp), parameter :: d = -0.86097d0
    real(dp), parameter :: e = -2.52692d0

    logical :: debug = .false.

    ! coulomb coupling parameter for species with charge Z
    gamma = cc% gamma_e * pow(Z,5d0/3d0)

    ! ratios for convenience
    gamma14 = pow(gamma, 1d0/4d0)

    mu = a * gamma + 4d0 * b * gamma14  - 4d0 * c / gamma14 + d * log(gamma) + e

    return

  end function mui_coulomb_I1993


  function mui_coulomb_PCR2009(cc, Z) result(mu)
    use math_lib
    use rates_def, only: Coulomb_Info

    ! calulate the chemical potential of an ion

    implicit none

    type(Coulomb_Info), intent(in) :: cc
    real(dp), intent(in) :: Z ! nuclear charge

    real(dp) :: mu ! the chemical potential in units of kT

    real(dp) :: a, b
    real(dp) :: alpha, beta, gamma, delta

    logical :: debug = .false.

    ! expressions taken from 
    ! [CP98]: Chabrier, G., \& Potekhin, A.~Y.\ 1998, \pre, 58, 4941
    ! [PC00]: Potekhin, A.~Y., \& Chabrier, G.\ 2000, \pre, 62, 8554 
    ! [P+09a]: Potekhin, A.~Y., Chabrier, G., \& Rogers, F.~J.\ 2009, \pre, 79, 016411
    ! [P+09b]: Potekhin, A.~Y., Chabrier, G., Chugunov, A.~I., Dewitt, H.~E., \& Rogers, F.~J.\ 2009, \pre, 80, 047401

    ! coulomb coupling parameter for species with charge Z
    gamma = cc% gamma_e * pow(Z, 5d0/3d0)

    ! first, calculate the shift from the linear mixing relation (LMR)
    ! see equation 2 from [P+09b]

    mu = fii(gamma) + fie(cc% rs, cc% gamma_e, Z)

    ! no need to do LMR corrections for now

    ! ! next, calculate the correction to the linear mixing relation (LMR)

    ! ! equations 10a, 11 & 12 from [P+09b]
    ! alpha = pow(zbar, 2d0/5d0) / pow(z2bar, 1d0/5d0)
    ! beta = 3d0 / (2d0 * alpha) - 1d0
    ! delta = 1d0 - pow(zbar2, 3d0/2d0) / (sqrt(zbar) * z52bar)
    ! a = (2.6d0 * delta + 14 * pow(delta, 3d0)) / (1d0 - alpha)
    ! b = 0.0117d0 * pow(z2bar/(zbar*zbar), 2d0) * a

    ! ! equation 9 from [P+09b]
    ! deltaf = (gamma_e * z52bar / sqrt(3d0)) * & 
    !      (delta / ((1d0 + a * pow(gamma, alpha)) * pow(1d0 + b * pow(gamma, alpha), beta)))



    return

  contains

    function fii(gamma) result(f)

      implicit none

      real(dp), intent(in) :: gamma
      real(dp) :: f

      ! define parameters
      real(dp), parameter :: A1 = -0.907347d0
      real(dp), parameter :: A2 = 0.62849d0
      real(dp) :: A3
      real(dp), parameter :: B1 = 0.0045d0
      real(dp), parameter :: B2 = 170d0
      real(dp), parameter :: B3 = -8.4d-5
      real(dp), parameter :: B4 = 0.0037d0

      A3 = -sqrt(3d0)/2d0 - A1/sqrt(A2)

      ! [CP00] Eq. (28)
      f = A1 * (sqrt(gamma * (A2 + gamma)) - A2 * log(sqrt(gamma/A2) + sqrt(1+gamma/A2))) + &
          2 * A3 * (sqrt(gamma) - atan(sqrt(gamma)))

      return 

    end function fii

    function fie(rs, gamma_e, Z) result(f)

      implicit none

      real(dp), intent(in) :: rs
      real(dp), intent(in) :: gamma_e
      real(dp), intent(in) :: Z

      real(dp) :: f

      real(dp) :: logZ
      real(dp) :: c_DH, c_inf, c_TF, g1, g2, h1, h2
      real(dp) :: a, b, x, nu

      ! [CP00], between eq 31 and 32
      logZ = log(Z)
      a = 1.11_dp * pow(Z, 0.475_dp)
      b = 0.2_dp + 0.078_dp * logZ * logZ
      nu = 1.16_dp + 0.08_dp * logZ

      ! [CP00], eq 30 and eq 31
      c_DH = Z / sqrt(3d0) * (pow(1+Z,3d0/2d0) - 1d0 - pow(Z,3d0/2d0))
      c_inf = 0.2513_dp
      c_TF = c_inf * pow(Z,7d0/3d0) * (1d0 - pow(Z,-1d0/3d0) + 0.2d0 * pow(Z,-1d0/2d0))

      ! [CP00], eq 32 and eq 33
      g1 = 1 + (0.78d0 / (21d0 + gamma_e * pow(Z/rs, 3d0))) * sqrt(gamma_e/Z)
      g2 = 1 + (Z-1d0)/9d0 * (1d0 + 1d0/(0.001d0*Z*Z + 2*gamma_e)) * pow(rs,3d0) / (1d0 + 6d0 * rs * rs)

      ! [CP00] eq 4
      x = 0.014005_dp / rs

      ! [CP98], below eq 33
      h1 = 1d0 / (1d0 + pow(x/sqrt(1+x*x),6d0) * pow(Z,-1d0/3d0))
      h2 = 1d0 / sqrt(1+x*x)

      f = - gamma_e * (c_DH * sqrt(gamma_e) + c_TF * a * pow(gamma_e, nu) * g1 * h1) / & 
                      (1d0 + (b * sqrt(gamma_e) + a * g2 * pow(gamma_e, nu)/rs) * h2)


    end function fie

  end function mui_coulomb_PCR2009



  function do_Vs_coulomb(cc, Z) result(Vs)

    use rates_def, only: Coulomb_Info, which_Vs_coulomb

    implicit none

    type(Coulomb_Info), intent(in) :: cc
    real(dp), intent(in) :: Z ! nuclear charge

    real(dp) :: Vs ! the screening potential in units of the fermi energy

    select case (which_Vs_coulomb)
    case (None)
       Vs = 0
    case (ThomasFermi)
       Vs = Vs_coulomb_TF(cc, Z)
    case (Itoh2002)
       Vs = Vs_coulomb_Itoh(cc, Z)
    case default
       stop "STOP (eval_coulomb.f: do_Vs_coulomb): incorrect value of which_Vs_coulomb"
    end select

    return

  end function do_Vs_coulomb


  function Vs_coulomb_TF(cc, Z) result(Vs)
    use const_def, only : fine, pi
    use math_lib
    use rates_def, only: Coulomb_Info

    implicit none

    type(Coulomb_Info), intent(in) :: cc
    real(dp), intent(in) :: Z ! nuclear charge

    real(dp) :: Vs ! the screening potential in units of the fermi energy

    ! this uses the thomas-fermi approximation, in the relativistic limit
    Vs = 2d0 * Z * fine * sqrt(fine/pi)

    return

  end function Vs_coulomb_TF


  function Vs_coulomb_Itoh(cc, Z) result(Vs)
    use math_lib
    use rates_def, only: Coulomb_Info

    implicit none

    type(Coulomb_Info), intent(in) :: cc
    real(dp), intent(in) :: Z ! nuclear charge

    real(dp) :: Vs ! the screening potential in units of the fermi energy

    ! code from Itoh, N., Tomizawa, N., Tamamura, M., Wanajo, S., & Nozawa, S. 2002, ApJ, 579, 380 

    integer :: i
    real(dp) :: rs, rs0, s, fj
    real(dp), dimension(11), parameter :: c(0:10) = (/ &
        0.450861D-01, 0.113078D-02, 0.312104D-02, 0.864302D-03, &
        0.157214D-01, 0.816962D-01, 0.784921D-01,-0.680863D-01, &
       -0.979967D-01, 0.204907D-01, 0.366713D-01 /)

    rs = cc% rs

    rs0=0
    if(rs.lt.0.00001d0)then
       rs0=1d0-rs
       rs=0.00001d0
    endif
    s=(log10(rs)+3d0)/2d0
    fj=0
    if(s.ne.0)then
       do i=0,10
          fj=c(i)*pow(s,dble(i))+fj
       enddo
    else
       fj=c(0)
    endif
    if(rs0.ne.0)then
       rs=1d0-rs0
    endif

    Vs = 1.46d-2 * Z * fj

    return

  end function Vs_coulomb_Itoh



end module eval_coulomb
