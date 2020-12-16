C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: pteh_pi.f 370 2006-11-29 23:57:39Z airwin $
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
      subroutine pteh_pi(ifmodified,
     &  re, ref, ret, t,
     &  dve, dvef, dvet)
C       routine for modelling PTEH pressure ionization
C      input quantities:
C      ifmodified >0 means_use_modified PTEH pressure ionization model
C      re, ref, ret are the scaled EFF number density of electrons and
C       the derivative of ln re wrt ln f and ln t.
C      t is the temperature.
C      output quantities:
C      dve, dvef, and dvet are the change in ln equilibrium constant due
C      to the electron contribution from the PTEH pressure-ionization
C      free energy model and ln f and ln t derivatives of same.
      implicit none
      include 'constants.h'
      integer ifmodified
      double precision
     &  re, ref, ret, t,
     &  dve, dvef, dvet,
     &  full_sum1, rho, rf, rt, ne, nef, net,
     &  ppi, ppif, ppit, ppir, spi, spif, spit, upi
      double precision c1pi, c2pi, c3pi, c4pi, c5pi, c6pi,
     &  onethird, onethirdm1, onethirdm2, 
     &  xpi, ypi, zpi, arg1pi, arg2pi,
     &  psipi, dpsipi, d2psipi, arg3pi,
     &  da, dax, day, daxx, daxy, dayy,
     &  daf, dat, daxf, daxt, dayf, dayt, 
     &  da0, da0x, da0y, da0xx, da0xy, da0yy,
     &  da0f, da0t, da0xf, da0xt, da0yf, da0yt,
     &  c1pi_original, c1pi_modified,
     &  c2pi_original, c2pi_modified,
     &  c3pi_original, c3pi_modified,
     &  c4pi_original, c4pi_modified
C       pressure ionization setup for chemical potential
C       pteh constants
      parameter(c1pi_original = 3.d0)
      parameter(c2pi_original = 0.25d0)
      parameter(c3pi_original = 2.d0)
      parameter(c4pi_original = 0.03d0)
C       fit to Best EOS (with H2 and H2+) for Z = 0
      parameter(c1pi_modified = 15.d0)
      parameter(c2pi_modified = 0.25d0)
      parameter(c3pi_modified = 2.d0)
      parameter(c4pi_modified = 0.03d0)
      parameter(c5pi = 1.07654d0)
      parameter(c6pi = 0.61315d0)
      parameter(onethird = 1.d0/3.d0)
      parameter(onethirdm1 = onethird-1.d0)
      parameter(onethirdm2 = onethird-2.d0)
C       copy value used by pteh for H2 ionization in K^-1
C       this is a fairly inaccurate value, but worth propagating since
C       the pteh pressure ionization approximation
C       was determined using this value.
      double precision hion
      parameter(hion=11605.0d0*13.6d0)
      double precision free_pi, free_pif, free_pir
      save
C       free energy per unit volume
C       fpi = -k T [n_e g(x,y) - n_e0 g(x0,y)], where
C       n_e is the number density of electrons,
C       n_e0 is the total number density of electrons whether bound
C       into atoms and molecules or not.
C       g(x,y) = exp(-(c1pi/x)**c2pi)*
C         (y + psiprime(z) + c3pi*ln(1 + x/c4pi)
C       psiprime(z) = 2 sqrt(c5pi z  (1 + c6pi z)^{1/3})
C       z = x y^{3/2}
C       x = n_e*H (in PTEH paper H *might* be multiplied by 1.008, but
C         we_use_simpler definition here that would be equivalent
C         to changing c1pi and c4pi by factor of 1.008.  n.b. checked
C         onno pols statef, and he uses our definition.)
C       x0 = n_e0 * mass (grams) of the hydrogen atom
C       y  = 13.60 (ev)/kt = 13.60 * 11605.0/t = hion/t
C       chemical potential = mu_e = partial fpi/partial n_e
C       mu_e = -k T (da + dax), where
C       da = g(x,y), dax = partial g(x,y)/partial ln x
      if(ifmodified.gt.0) then
        c1pi = c1pi_modified
        c2pi = c2pi_modified
        c3pi = c3pi_modified
        c4pi = c4pi_modified
      else
        c1pi = c1pi_original
        c2pi = c2pi_original
        c3pi = c3pi_original
        c4pi = c4pi_original
      endif
      xpi = cd*re
      arg3pi = -(c1pi/xpi)**c2pi
C       underflow protection
      if(arg3pi.lt.-70.d0) then
        da = 0.d0
        dax = 0.d0
        day = 0.d0
        daxx = 0.d0
        daxy = 0.d0
        dayy = 0.d0
        daf = 0.d0
        dat = 0.d0
        daxf = 0.d0
        daxt = 0.d0
        dayf = 0.d0
        dayt = 0.d0
      else
        arg3pi = exp(arg3pi)
        ypi = hion/t
        zpi = xpi*ypi**1.5d0
        arg1pi = c5pi*zpi
        arg2pi = 1.d0 + c6pi*zpi
        psipi = 2.d0*sqrt(arg1pi*arg2pi**onethird)
C         derivatives wrt to zpi
        dpsipi = 2.d0/psipi*c5pi*(
     &    arg2pi**onethird + zpi*c6pi*arg2pi**onethirdm1*onethird)
        d2psipi = -dpsipi*dpsipi/psipi + 2.d0/psipi*c5pi*(
     &    2.d0*c6pi*arg2pi**onethirdm1*onethird +
     &    zpi*c6pi*c6pi*arg2pi**onethirdm2*onethirdm1*onethird)
C         da = g(x,y)
        da = arg3pi*(ypi + psipi + c3pi*log(1.d0 + xpi/c4pi))
C         dax = partial g(x,y)/partial ln x
        dax = da*((c1pi/xpi)**(c2pi)*c2pi) +
     &    arg3pi*(dpsipi*zpi + c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi)
C         day = partial g(x,y)/partial ln y
        day = arg3pi*(ypi + 1.5d0*dpsipi*zpi)
C         first calculate partial dax/partial ln x, ln y
        daxx = dax*((c1pi/xpi)**(c2pi)*c2pi) - da*
     &    ((c1pi/xpi)**(c2pi)*c2pi*c2pi) +
     &    arg3pi*((dpsipi*zpi + c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi)*
     &    ((c1pi/xpi)**(c2pi)*c2pi) +
     &    d2psipi*zpi*zpi + dpsipi*zpi +
     &    c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi*
     &    (1.d0 - 1.d0/(1.d0 + xpi/c4pi)*xpi/c4pi))
        daxy = day*((c1pi/xpi)**(c2pi)*c2pi) +
     &    arg3pi*1.5d0*(d2psipi*zpi + dpsipi)*zpi
        dayy = arg3pi*(ypi + 1.5d0*1.5d0*(d2psipi*zpi + dpsipi)*zpi)
        daf = dax*ref
        dat = dax*ret - day
C         transform from ln x, ln y derivatives to ln f, ln t derivatives.
        daxf = daxx*ref
        dayf = daxy*ref
        daxt = daxx*ret - daxy
        dayt = daxy*ret - dayy
      endif
      dve = (da + dax)
      dvef = (daf + daxf)
      dvet = (dat + daxt)
      return
      entry pteh_pi_pressure(
     &  full_sum1, rho, rf, rt, t, ne, nef, net,
     &  ppi, ppif, ppit, ppir)
C       pteh pressure ionization
C       calculate ppi, ppif, and ppir, the component of the pressure
C       due to PTEH pressure ionization plus its fl and rl derivatives for
C       fixed tl.
C       calculation for n_e0 = (rho/H * full_sum1) not n_e
      xpi = full_sum1*rho
      arg3pi = -(c1pi/xpi)**c2pi
C       underflow protection
      if(arg3pi.lt.-70.d0) then
        da0 = 0.d0
        da0x = 0.d0
        da0y = 0.d0
        da0xx = 0.d0
        da0xy = 0.d0
        da0xt = 0.d0
      else
        arg3pi = exp(arg3pi)
        ypi = hion/t
        zpi = xpi*ypi**1.5d0
        arg1pi = c5pi*zpi
        arg2pi = 1.d0 + c6pi*zpi
        psipi = 2.d0*sqrt(arg1pi*arg2pi**onethird)
        dpsipi = 2.d0/psipi*c5pi*(
     &    arg2pi**onethird + zpi*c6pi*arg2pi**onethirdm1*onethird)
        d2psipi = -dpsipi*dpsipi/psipi + 2.d0/psipi*c5pi*(
     &    2.d0*c6pi*arg2pi**onethirdm1*onethird +
     &    zpi*c6pi*c6pi*arg2pi**onethirdm2*onethirdm1*onethird)
C         da0 = g(x0,y)
        da0 = arg3pi*(ypi + psipi + c3pi*log(1.d0 + xpi/c4pi))
C         da0x = partial g(x0,y)/partial ln x0
        da0x = da0*((c1pi/xpi)**(c2pi)*c2pi) +
     &    arg3pi*(dpsipi*zpi + c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi)
C         da0y = partial g(x0,y)/partial ln y
        da0y = arg3pi*(ypi + 1.5d0*dpsipi*zpi)
C         partial da0x/partial ln x0
        da0xx = da0x*((c1pi/xpi)**(c2pi)*c2pi) - da0*
     &    ((c1pi/xpi)**(c2pi)*c2pi*c2pi) +
     &    arg3pi*((dpsipi*zpi + c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi)*
     &    ((c1pi/xpi)**(c2pi)*c2pi) +
     &    d2psipi*zpi*zpi + dpsipi*zpi +
     &    c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi*
     &    (1.d0 - 1.d0/(1.d0 + xpi/c4pi)*xpi/c4pi))
        da0xy = da0y*((c1pi/xpi)**(c2pi)*c2pi) +
     &    arg3pi*1.5d0*(d2psipi*zpi + dpsipi)*zpi
        da0xt = da0xx*rt - da0xy
      endif
C       free energy per unit volume
C       fpi = -k T [n_e g(x,y) - n_e0 g(x0,y)], where
C       n_e is the number density of electrons,
C       n_e0 is the total number density of electrons whether bound
C       into atoms and molecules or not.
C       delta P = - partial (fpi(T, N_e, V) V)/partial V,
C         = -fpi + n_e mu_e  + n_e0 mu_e0
C         = -kT (n_e dax - n_e0 da0x)
C       where da0x = partial g(x0, y)/partial ln x0
      ppi = -cr*rho*t*(ne*dax - full_sum1*da0x)
      ppif = ppi*rf -
     &  cr*rho*t*(nef*dax + ne*daxf - full_sum1*da0xx*rf)
      ppit = ppi*(rt + 1.d0) -
     &  cr*rho*t*(net*dax + ne*daxt - full_sum1*da0xt)
C      partial wrt ln rho for fixed fl and tl where
C      ne is inversely proportional to rho, and da(x,y) is independent of
C      rho because x is proportional to n_e(f, t) and y depends on t.
C      simplify expression below to get rid of cancellation loss (due
C      to rho*ne being independent of rho).
C      ppir = ppi -
C     &  cr*rho*t*(-ne*dax - full_sum1*da0xx)
      ppir = cr*rho*t*full_sum1*(da0x + da0xx)
      return
      entry pteh_pi_free(full_sum1, rho, rf, t, ne, nef,
     &  free_pi, free_pif, free_pir)
C       pteh pressure ionization
C       calculate free_pi = upi -t*spi
      xpi = full_sum1*rho
      arg3pi = -(c1pi/xpi)**c2pi
C       underflow protection
      if(arg3pi.lt.-70.d0) then
        da0 = 0.d0
        da0x = 0.d0
      else
        arg3pi = exp(arg3pi)
        ypi = hion/t
        zpi = xpi*ypi**1.5d0
        arg1pi = c5pi*zpi
        arg2pi = 1.d0 + c6pi*zpi
        psipi = 2.d0*sqrt(arg1pi*arg2pi**onethird)
        dpsipi = 2.d0/psipi*c5pi*(
     &    arg2pi**onethird + zpi*c6pi*arg2pi**onethirdm1*onethird)
C         da0 = g(x0,y)
        da0 = arg3pi*(ypi + psipi + c3pi*log(1.d0 + xpi/c4pi))
C        da0x = partial g(x0,y)/partial ln x0
        da0x = da0*((c1pi/xpi)**(c2pi)*c2pi) +
     &    arg3pi*(dpsipi*zpi + c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi)
      endif
C      free energy per unit volume
C      fpi = -k T [n_e g(x,y) - n_e0 g(x0,y)], where
C      n_e is the number density of electrons,
C      n_e0 is the total number density of electrons whether bound
C      into atoms and molecules or not.
C      N.B. ne = n_e/(Navogadro*rho) and full_sum1 = sum_i eps_i charge_i
C      is the value of ne (NOT n_e) for full ionization.
C      hence, free_pi is per unit mass.
      free_pi = -cr*t*(ne*da - full_sum1*da0)
      free_pif = -cr*t*(nef*da + ne*daf - full_sum1*da0x*rf)
C      partial of free energy/unit mass wrt ln rho at fixed f, t, where
C      ne is inversely proportional to rho, and da(x,y) is independent of
C      rho because x is proportional to n_e(f, t) and y depends on t.
      free_pir = -cr*t*(-ne*da - full_sum1*da0x)
      return
      entry pteh_pi_end(
     &  full_sum1, rho, rf, rt, t, ne, nef, net,
     &  ppi, ppif, ppit, spi, spif, spit, upi)
C       pteh pressure ionization
C       calculate ppi, spi, ln f, and ln t derivatives, and upi.
C       see commentary above when calculating chemical potential
C       calculation for n_e0 = (rho/H * full_sum1) not n_e
      xpi = full_sum1*rho
      arg3pi = -(c1pi/xpi)**c2pi
C       underflow protection
      if(arg3pi.lt.-70.d0) then
        da0 = 0.d0
        da0x = 0.d0
        da0y = 0.d0
        da0xx = 0.d0
        da0xy = 0.d0
        da0yy = 0.d0
        da0f = 0.d0
        da0t = 0.d0
        da0xf = 0.d0
        da0xt = 0.d0
        da0yf = 0.d0
        da0yt = 0.d0
      else
        arg3pi = exp(arg3pi)
        ypi = hion/t
        zpi = xpi*ypi**1.5d0
        arg1pi = c5pi*zpi
        arg2pi = 1.d0 + c6pi*zpi
        psipi = 2.d0*sqrt(arg1pi*arg2pi**onethird)
        dpsipi = 2.d0/psipi*c5pi*(
     &    arg2pi**onethird + zpi*c6pi*arg2pi**onethirdm1*onethird)
        d2psipi = -dpsipi*dpsipi/psipi + 2.d0/psipi*c5pi*(
     &    2.d0*c6pi*arg2pi**onethirdm1*onethird +
     &    zpi*c6pi*c6pi*arg2pi**onethirdm2*onethirdm1*onethird)
C         da0 = g(x0,y)
        da0 = arg3pi*(ypi + psipi + c3pi*log(1.d0 + xpi/c4pi))
C         da0x = partial g(x0,y)/partial ln x0
        da0x = da0*((c1pi/xpi)**(c2pi)*c2pi) +
     &    arg3pi*(dpsipi*zpi + c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi)
C         da0y = partial g(x0,y)/partial ln y
        da0y = arg3pi*(ypi + 1.5d0*dpsipi*zpi)
C         first calculate partial da0x/partial ln x0, ln y
        da0xx = da0x*((c1pi/xpi)**(c2pi)*c2pi) - da0*
     &    ((c1pi/xpi)**(c2pi)*c2pi*c2pi) +
     &    arg3pi*((dpsipi*zpi + c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi)*
     &    ((c1pi/xpi)**(c2pi)*c2pi) +
     &    d2psipi*zpi*zpi + dpsipi*zpi +
     &    c3pi/(1.d0 + xpi/c4pi)*xpi/c4pi*
     &    (1.d0 - 1.d0/(1.d0 + xpi/c4pi)*xpi/c4pi))
        da0xy = da0y*((c1pi/xpi)**(c2pi)*c2pi) +
     &    arg3pi*1.5d0*(d2psipi*zpi + dpsipi)*zpi
        da0yy = arg3pi*(ypi + 1.5d0*1.5d0*(d2psipi*zpi + dpsipi)*zpi)
        da0f = da0x*rf
        da0t = da0x*rt - da0y
C         transform from ln x, ln y derivatives to ln f, ln t derivatives.
        da0xf = da0xx*rf
        da0yf = da0xy*rf
        da0xt = da0xx*rt - da0xy
        da0yt = da0xy*rt - da0yy
      endif
C       free energy per unit volume
C       fpi = -k T [n_e g(x,y) - n_e0 g(x0,y)], where
C       n_e is the number density of electrons,
C       n_e0 is the total number density of electrons whether bound
C       into atoms and molecules or not.
C       delta P = - partial (fpi(T, N_e, V) V)/partial V,
C         = -fpi + n_e mu_e  + n_e0 mu_e0
C         = -kT (n_e dax - n_e0 da0x)
C       where da0x = partial g(x0, y)/partial ln x0
      ppi = -cr*rho*t*(ne*dax - full_sum1*da0x)
      ppif = ppi*rf -
     &  cr*rho*t*(nef*dax + ne*daxf - full_sum1*da0xf)
      ppit = ppi*(rt + 1.d0) -
     &  cr*rho*t*(net*dax + ne*daxt - full_sum1*da0xt)
C       delta s (per unit mass) = - partial fpi(T, n_e)/partial T/rho
C         = -fpi/(rho T) - k/rho * 
C       [n_e partial g(x,y)/partial ln y - n_e0 partial g(x0,y)/partial ln y]
      spi = cr*(ne*da - full_sum1*da0 - (ne*day - full_sum1*da0y))
      spif = cr*(nef*da + ne*daf - full_sum1*da0f -
     &  (nef*day + ne*dayf - full_sum1*da0yf))
      spit = cr*(net*da + ne*dat - full_sum1*da0t -
     &  (net*day + ne*dayt - full_sum1*da0yt))
C       the energy per unit *mass* = fpi/rho + T delta s
      upi = -cr*t*(ne*day - full_sum1*da0y)
      end
