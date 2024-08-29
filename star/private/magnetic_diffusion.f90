! ***********************************************************************
!
!   Copyright (C) 2010-2021  The MESA Team
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


      module magnetic_diffusion

      use const_def
      use num_lib
      use math_lib
      use utils_lib

      implicit none

      private
      public :: calc_sige, calc_eta

      contains

      !> Compute the magnetic diffusivity from the electric conductivity.
      !! @param sig The electrical conductivity (1/s).
      !! @param eta The magnetic diffusivity (output, cm^2/s).
      real(dp) function calc_eta(sig) result(eta)
         real(dp), intent(in) :: sig

         eta = (clight * clight / (4d0 * pi)) /sig
      end function calc_eta

      !> Computes the electrical conductivity following
      !! S.-C. YOON Oct. 10, 2003.
      !!
      !! @param abar The mean atomic mass number.
      !! @param zbar The mean atomic charge.
      !! @param rho The density (g/cm^3).
      !! @param T The temperature (K).
      !! @param Cp The specific heat at constant pressure (erg/g/K).
      !! @param kap_cond The electronic thermal opacity (cm^2/g).
      !! @param opacity The opacity (cm^2/g).
      !! @param sig The electrical conductivity (output, 1/s).
      real(dp) function calc_sige(abar, zbar, rho, T, Cp, kap_cond, opacity) result(sig)
         real(dp), intent(in) :: abar, zbar, rho, T, Cp, kap_cond, opacity
         real(dp) :: gamma, xlg, alpha, ffff, xxx, xsig1, xsig2, xsig3

         gamma = 0.2275d0*pow2(zbar) * pow(rho * 1.d-6 / abar, one_third)*1.d8/T
         xlg = log10(gamma)

         alpha = 16d0 * boltz_sigma * pow3(T) / (3d0 * opacity * pow2(rho) * Cp)

         if (xlg < -1.5d0) then
            sig = sige1(zbar,T,gamma)
         else if (xlg >= -1.5d0 .and. xlg <= 0d0) then
            xxx = (xlg + 0.75d0)*4d0/3d0
            ffff = 0.25d0*(2d0-3d0*xxx + pow3(xxx))
            xsig1 = sige1(zbar,T,gamma)
            xsig2 = sige2(T,rho,kap_cond)
            sig = (1d0-ffff)*xsig2 + ffff*xsig1
         else if (xlg > 0d0 .and. xlg < 0.5d0) then
            xsig2 = sige2(T,rho,kap_cond)
            sig = xsig2
         else if (xlg >= 0.5d0 .and. xlg < 1d0) then
            xxx = (xlg-0.75d0)*4d0
            ffff = 0.25d0*(2d0-3d0*xxx + pow3(xxx))
            xsig2 = sige2(T,rho,kap_cond)
            xsig3 = sige3(zbar,T,gamma)
            sig = (1d0-ffff)*xsig3 + ffff*xsig2
         else
            sig = sige3(zbar,T,gamma)
         endif

      end function calc_sige

      !> Computes one regime of the electrical conductivity.
      !! Written by S.-C. Yoon, Oct. 10, 2003
      !! See also Spitzer 1962 and Wendell et al. 1987, ApJ 313:284
      !! @param Z species charge
      !! @param T Temperature (K)
      !! @param xgamma The ion coupling strength (dimensionless).
      !! @param sige1 The electrical conductivity (1/s).
      real(dp) function sige1(z,t,xgamma)
         real(dp), intent(in) :: z, t, xgamma
         real(dp) :: etan, xlambda,f
         if (t >= 4.2d5) then
            f = sqrt(4.2d5/t)
         else
            f = 1.d0
         end if
         xlambda = sqrt(3d0*z*z*z)*pow(xgamma,-1.5d0)*f + 1d0
         etan = 3.d11*z*log(xlambda)*pow(t,-1.5d0)             ! magnetic diffusivity
         etan = etan/(1.d0-1.20487d0*exp(-1.0576d0*pow(z,0.347044d0))) ! correction: gammae
         sige1 = clight*clight/(pi4*etan)                    ! sigma = c^2/(4pi*eta)
      end function sige1

      !> Computes one regime of the electrical conductivity using the conductive opacity.
      !! Written by S.-C. Yoon, Oct. 10, 2003
      !! See Wendell et al. 1987, ApJ 313:284
      !! @param T Temperature (K)
      !! @param rho Temperature (g/cm^3)
      !! @param kap_cond The electronic thermal opacity (cm^2/g).
      !! @param sige2 The electrical conductivity (1/s).
      real(dp) function sige2(T,rho,kap_cond)
         real(dp), intent(in) :: t,rho,kap_cond
         sige2 = 1.11d9*T*T/(rho*kap_cond)
      end function sige2

      !> Computes the electrical conductivity in degenerate matter.
      !! Written by S.-C. Yoon, Oct. 10, 2003
      !! See Nandkumar & Pethick (1984)
      !! @param Z species charge
      !! @param T Temperature (K)
      !! @param xgamma The ion coupling strength (dimensionless).
      !! @param sige3 The electrical conductivity (1/s).
      real(dp) function sige3(z,t,xgamma)
         real(dp), intent(in) :: z, t, xgamma
         real(dp) :: rme, rm23, ctmp, xi
         rme = 8.5646d-23*t*t*t*xgamma*xgamma*xgamma/pow5(z)  ! rme = rho6/mue
         rm23 = pow(rme,2d0/3d0)
         ctmp = 1d0 + 1.018d0*rm23
         xi= sqrt(3.14159d0/3.)*log(z)/3.d0 + 2.d0*log(1.32d0+2.33d0/sqrt(xgamma))/3.d0-0.484d0*rm23/ctmp
         sige3 = 8.630d21*rme/(z*ctmp*xi)
      end function sige3


      end module magnetic_diffusion