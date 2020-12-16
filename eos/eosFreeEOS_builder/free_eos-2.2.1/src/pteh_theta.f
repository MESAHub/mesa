C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: pteh_theta.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine pteh_theta(x, t, rhostar,
     &  theta_e, theta_ef, theta_et, theta_eff, theta_eft, theta_ett,
     &  thetane, thetanef, thetanet,
     &  thetat, thetatf, thetatt)
C       routine for calculating theta_e and its derivatives from approximate
C       relation given in pteh eq. (29).  this calculation only replaces
C       the ordinary fermi-dirac calculation under very special circumstances
C       when attempting exact mimicry of the pteh approach.
C       x = number density of electrons/N_a = cd*re
C       t = temperature
C       rhostar(9), in order of undefined (1); partial ln n_e wrt ln f (2);
C         ln t (3); ff(4); ft(5); tt(6); fff(7); fft(8); ftt(9).
      implicit none
      double precision x, t, rhostar(9),
     &  theta_e, theta_ef, theta_et, theta_eff, theta_eft, theta_ett,
     &  thetane, thetanef, thetanet,
     &  thetat, thetatf, thetatt
      double precision c5, c6,
     &  onethird,
     &  y, z, arg1, arg2, darg2,
     &  f, df, d2f, d3f, dpsidf, dpsidff, dpsidfff
      parameter(c5 = 1.07654d0)
      parameter(c6 = 0.61315d0)
      parameter(onethird = 1.d0/3.d0)
C       copy value used by pteh for H2 ionization in K^-1
C       this is a fairly inaccurate value, but worth propagating since
C       the pteh pressure ionization approximation
C       was determined using this value.
      double precision hion
      parameter(hion=11605.0d0*13.6d0)
      save
C       psi(z) = 2 sqrt(f) where
C       f = c5 z  (1 + c6 z)^{1/3})
C       z = x y^{3/2}
C       y  = 13.60 (ev)/kt = 13.60 * 11605.0/t = hion/t
      y = hion/t
      z = x*y**1.5d0
      arg1 = c5*z
      darg2 = c6*z
      arg2 = 1.d0 + darg2
      f = arg1*arg2**onethird
C       derivatives ln f wrt to ln z
      df = 1.d0 + onethird*darg2/arg2
      d2f = onethird*darg2*(1.d0 - darg2/arg2)/arg2
      d3f = d2f + onethird*darg2*darg2*(-2.d0 + 2.d0*darg2/arg2)/
     &  arg2/arg2
C       derivatives of psi wrt fl based on usual relation for psi(f).
      dpsidf = sqrt(1.d0+f)
      dpsidff = 0.5d0*f*dpsidf/(1.d0+f)
      dpsidfff = dpsidff*(1.d0 - 0.5d0*f/(1.d0+f))
C       theta_e = partial ln ne(psi, t)/wrt psi
C       implicit differentiation: ln ne = const + ln z(f) + 3/2 ln t
C       also note z is a function of f alone so theta_e is a function of f
C       alone:
C       d ln z/d ln f = 1/df
C       d2 ln z/... = - d2f/df^3
C       d3 ln z/... = - d3f/df^4 + 3 d2f^2/df^5
C       theta_e = (d ln z/d ln f )/ (d psi/ d ln f)
      theta_e = 1.d0/(df*dpsidf)
C       derivatives wrt local f, t.  Note the t derivatives are zero
C       because theta_e depends on z(f) and f.
      theta_ef = (-d2f/(df*df) - dpsidff/dpsidf)*theta_e
      theta_eff = ((-d3f + 3.d0*d2f*d2f/df)/df**3 +
     &  2.d0*d2f/(df*df)*dpsidff/dpsidf -
     &  (dpsidfff-2.d0*dpsidff*dpsidff/dpsidf)/dpsidf)*theta_e
C       n.b. partial theta_e(ne,T)/partial ln ne =
C         theta_ef*df*partial ln z(ne,T)/ partial ln ne, where
C       ln z = const + ln ne - 3/2 ln t.
      thetane = theta_ef*df
C       derivatives wrt local f, t.  Note the t derivatives are zero
C       because thetane depends on z(f) and f.
      thetanef = theta_eff*df + theta_ef*d2f/df
C       n.b. partial theta_e(ne, T)/partial ln T = 
C         theta_ef*df*partial ln z(ne,T)/ partial ln T
      thetat = theta_ef*df*(-1.5d0)
C       derivatives wrt local f, t.  Note the t derivatives are zero
C       because thetat depends on z(f) and f.
      thetatf = theta_eff*df*(-1.5d0) + theta_ef*d2f/df*(-1.5d0) 
C       local approximate f is a function of z which is a function of 
C       ne(f,t) and t.
C       do appropriate transformation in reverse order so always dealing
C       with untransformed quantities on rhs.
      thetatt = thetatf*df*(rhostar(3)-1.5d0)
      thetatf = thetatf*df*rhostar(2)
      thetanet = thetanef*df*(rhostar(3)-1.5d0)
      thetanef = thetanef*df*rhostar(2)
      theta_ett = theta_eff*(df*(rhostar(3)-1.5d0))*
     &  (df*(rhostar(3)-1.5d0)) +
     &  theta_ef*(d2f*(rhostar(3)-1.5d0)*(rhostar(3)-1.5d0) +
     &  df*rhostar(6))
      theta_eft = theta_eff*(df*(rhostar(3)-1.5d0))*(df*rhostar(2)) +
     &  theta_ef*(d2f*(rhostar(3)-1.5d0)*rhostar(2) + df*rhostar(5))
      theta_eff = theta_eff*(df*rhostar(2))*(df*rhostar(2)) +
     &  theta_ef*(d2f*rhostar(2)*rhostar(2) + df*rhostar(4))
      theta_et = theta_ef*df*(rhostar(3)-1.5d0)
      theta_ef = theta_ef*df*rhostar(2)
      end
