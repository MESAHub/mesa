C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: tau_calc.f 352 2006-04-13 02:12:47Z airwin $
C
C       For the latest version of this source code, please contact
C       Alan W. Irwin
C       Department of Physics and Astronomy
C       University of Victoria,
C       Box 3055
C       Victoria, B.C., Canada
C       V8W 3P6
C       e-mail: irwin@beluga.phys.uvic.ca.
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C
C       End of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C*******************************************************************************
      subroutine tau_calc(sum0, sum2, t, thetaxne,
     &  tau, dtau, d2tau)
C       calculate tau(x) (see MDH Paper II, eqs. 11 and 12)
C       and its first derivatives dtau(4) wrt the 4 input parameters,
C       and its second derivatives d2tau(4,4) [lower triangle filled in
C       from symmetry arguments] wrt to the 4 input parameters:
C       sum0 = sum_ions nion, where nion is the number per unit volume,
C         and the sum taken over positive ions.
C       sum2 = ne thetae + sum_ions Z^2 nion,
C         ne is the electron number density, thetae is the degeneracy
C         correction (see PTEH or notes), Z is the charge on the ions.
C       t = temperature in K
C       thetaxne = (n_e k T/P_e)*n_e, where n_e is the number density of
C         free electrons and P_e is associated pressure.
      implicit none
      double precision sum0, sum2, t, thetaxne, tau, dtau(4), d2tau(4,4),
     &  lnx, dlnx(4), d2lnx(4), x, dtaulnx, d2taulnx
      integer i,j
      call lnx_calc(sum0, sum2, t, thetaxne, lnx, dlnx, d2lnx)
      x = exp(lnx)
      if(x.lt.0.1d0) then
C         tau(x) = (3/x^3)(x^3/3-x^4/4+x^5/5-...)
        tau = 1.d0 - 3.d0*x*(1.d0/4.d0 - x*(1.d0/5.d0 - x*(1.d0/6.d0 -
     &    x*(1.d0/7.d0 - x*(1.d0/8.d0 - x*(1.d0/9.d0 - x*(1.d0/10.d0 -
     &    x*(1.d0/11.d0 - x*(1.d0/12.d0 - x*(1.d0/13.d0 -
     &    x*(1.d0/14.d0 - x*(1.d0/15.d0 - x*(1.d0/16.d0 -
     &    x*1.d0/17.d0)))))))))))))
C         dtaulnx = d tau(x)/d ln x
        dtaulnx = - 3.d0*x*(1.d0/4.d0 - x*(2.d0/5.d0 - x*(3.d0/6.d0 -
     &    x*(4.d0/7.d0 - x*(5.d0/8.d0 - x*(6.d0/9.d0 - x*(7.d0/10.d0 -
     &    x*(8.d0/11.d0 - x*(9.d0/12.d0 - x*(10.d0/13.d0 -
     &    x*(11.d0/14.d0 - x*(12.d0/15.d0 - x*(13.d0/16.d0 -
     &    x*14.d0/17.d0)))))))))))))
C         d2taulnx = d 2 tau(x)/d 2(ln x)
        d2taulnx = - 3.d0*x*(1.d0/4.d0 - x*(4.d0/5.d0 - x*(9.d0/6.d0 -
     &    x*(16.d0/7.d0 - x*(25.d0/8.d0 -
     &    x*(36.d0/9.d0 - x*(49.d0/10.d0 -
     &    x*(64.d0/11.d0 - x*(81.d0/12.d0 - x*(100.d0/13.d0 -
     &    x*(121.d0/14.d0 - x*(144.d0/15.d0 - x*(169.d0/16.d0 -
     &    x*196.d0/17.d0)))))))))))))
      else
C         tau(x) = (3/x^3)*(ln(1+x) -x + x^2/2)
        tau = 3.d0*(log(1.d0 + x) - x*(1.d0 - 0.5d0*x))/(x*x*x)
C         dtaulnx = d tau(x)/d ln x
        dtaulnx = 3.d0*(-tau + 1.d0/(1.d0 + x))
C         d2taulnx = d 2 tau(x)/d 2(ln x)
        d2taulnx = 3.d0*(-dtaulnx - x/((1.d0 + x)*(1.d0 + x)))
      endif
      do i = 1,4
        dtau(i) = dtaulnx*dlnx(i)
        d2tau(i,i) = d2taulnx*dlnx(i)*dlnx(i) + dtaulnx*d2lnx(i)
        do j = i+1,4
          d2tau(j,i) = d2taulnx*dlnx(j)*dlnx(i)
          d2tau(i,j) = d2tau(j,i)
        enddo
      enddo
      end
