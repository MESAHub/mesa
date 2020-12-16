C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: coulomb.f 839 2008-07-07 18:59:04Z airwin $
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
      subroutine coulomb(ifcoulomb, if_dc, sum0, sum2, t, x,
     &  fcoulomb, fcoulomb0, fcoulomb2, fcoulombt, fcoulombx,
     &  fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0x,
     &  fcoulomb22, fcoulomb2t, fcoulomb2x, 
     &  fcoulombtt, fcoulombtx, fcoulombxx, lambda, gamma_e)
C       routine to calculate EOS quantitities relevant to the
C       Coulomb interaction using the PTEH approximation (ifcoulomb = 3, 4
C       pteh theta calculated outside, otherwise fermi-dirac theta_e calculated
C       outside) or DH smoothly joined by cubic to modified OCP
C       (preferred, ifcoulomb = 5)
C       or the Debye-Huckel limit (ifcoulomb = 1) or that limit corrected
C       by tau(x) (ifcoulomb = 2) for the Coulomb
C       free energy calculated for the partial ionization case.
C       input variables:
C       ifcoulomb controls whether Debye-Huckel (1), Debye-Huckel with tau(x)
C         correction (2), PTEH approximation (3 or 4 or 9), 
C         DH smoothly joined using cubic with modified OCP (preferred) (5)
C          or alternative version of cubic with different end points (6),
C         DH for Gamma <1 and OCP for Gamma >= 1 (for illustrative purposes)
C         (ifcoulomb = 7 or 8, with different theta_e values calculated by
C          calling routine)
C         ifcoulomb = 3 means theta_e by pteh approximation outside
C         ifcoulomb != 3 means theta_e by usual fermi-dirac formula.
C         ifcouloumb <= 4 means sum0 calculated outside with
C         PTEH definition otherwise calculated outside with DeWitt definition.
C         N.B. this distinction is meaningless for ifcoulomb = 0, 1, or 2.
C       if_dc = 0 means apply no diffraction correction to Lambda
C       if_dc = 1 means apply diffraction correction to Lambda
C       n.b. tau(x) correction and diffraction correction are not
C         simultaneously allowed.
C       sum0 = sum_ions nion, where nion is the number per unit volume,
C         and the sum taken over positive ions when PTEH definition is
C        used outside.  Otherwise n_e thetae has been added outside
C      (DeWitt definition).
C       sum2 = n_e thetae + sum_ions Z^2 nion,
C       n_e is the electron number density, thetae is the degeneracy
C       correction (see PTEH or notes), Z is the charge on the ions.
C       t = temperature in K
C       if if_dc = 0, then
C         x = thetax*n_e = (n_e k T/P_e)*n_e, where n_e is the 
C         number density of free electrons and P_e is associated pressure.
C       n.b. thetax*n_e is only relevant when the tau(x) correction is 
C         being applied to the Debye-Huckle limit, i.e., ifcoulomb = 2
C       if if_dc = 1, then
C         x = n_e thetae or possibly n_e
C       output variables:
C       fcoulomb = free energy per unit volume as a function of
C         sum0, sum2, t, and x.
C       fcoulomb0 = partial fcoulomb(sum0, sum2, t, x) wrt sum0
C       fcoulomb2 = partial fcoulomb(sum0, sum2, t, x) wrt sum2
C       fcoulombt = partial fcoulomb(sum0, sum2, t, x) wrt t
C       fcoulombx = partial fcoulomb(sum0, sum2, t, x) wrt x
C       fcoulomb00, etc., higher order derivatives
C       lambda = plasma interaction parameter (defined by PTEH and consistent 
C       with the simpler CG 15.60 criterion for single component plasmas.)
C       gamma_e ~ gamma_e_i of DeWitt, paper II, with a correction
C       that is only good (currently) for 1 >> gamma_e > Lambda.
      implicit none
      include 'constants.h'
      include 'pi_fit.h'
      double precision sum0, sum2, t, x, fcoulomb, 
     &  fcoulomb0, fcoulomb2, fcoulombt, fcoulombx,
     &  fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0x, 
     &  fcoulomb22, fcoulomb2t, fcoulomb2x, 
     &  fcoulombtt, fcoulombtx, fcoulombxx,
     &  lambda, lambda_const, gamma_e, gamma_e_const, gamma, zeta2,
     &  dc_const, dc_const1, dc_lambda,
!  1 dc_lambda0, 
     &  dc_lambda2, dc_lambdat, dc_lambdax,
!  1 dc_lambda00, dc_lambda02, dc_lambda0t, dc_lambda0x, 
     &  dc_lambda22, dc_lambda2t, dc_lambda2x, 
     &  dc_lambdatt, dc_lambdatx, dc_lambdaxx,
!  1 lambda0, 
     &  lambda2, lambdat, lambdax,
!  1 lambda00, lambda02, lambda0t, lambda0x, 
     &  lambda22, lambda2t, lambda2x, 
     &  lambdatt, lambdatx, lambdaxx,
     &  acon, dcon, gcon, 
     &  dh, ddh, ddh2, ocp, docp, docp2,
     &  gc, dgc, dgc2, tau, dtau(4), d2tau(4,4), ln10,
     &  xvalue, xdh, xmocp
      integer ifcoulomb, ifcoulomb_old, if_dc, ifstart
C      Must be different from 5 or 6.
      data ifcoulomb_old/-1/
      data ifstart/1/
C       PTEH values
      parameter(acon = 0.89752d0)
C      n.b. changing 0.208 to 0.5 gives good LMS fit, but bad solar
C      fit (at least by solar standards)
C      pteh values give even worse solar fit so went to dh for
C      solar conditions (log10(gamma) < -0.4) smoothly joined by cubic
C      to modified ocp result.
      parameter(dcon = 0.208d0*acon)
      parameter(gcon = -1.d0/0.768d0)
C       fitting coefficients from 1995, DeWitt, Slattery, and Chabrier
      double precision a, b, c, s, f1, d, dscacon, dscbcon, dscscon,
     &  dscccon, dscdcon, dscmccon, dscmdcon,
     &  dx, acap, bcap, ccap, dcap, ylo, yhi, dylo2, dyhi, dyhi2
      parameter(a = -0.899126d0)
      parameter(b = 0.60712d0)
      parameter(c = -0.27998d0)
      parameter(s = 0.321308d0)
      parameter(f1 = -0.436484d0)
C       n.b. f1 is from HNC calculation
      parameter(d = f1 -(a +b/s))
      parameter(dscacon = -a)
      parameter(dscbcon = -b/s)
      parameter(dscscon = s)
      parameter(dscccon = -c)
      parameter(dscdcon = -d)
      save
      if(ifcoulomb.eq.2.and.if_dc.eq.1)
     &  stop 'coulomb: ifcoulomb = 2 and if_dc = 1 not allowed'
      if(ifstart.eq.1) then
        ifstart = 0
        ln10 = log(10.d0)
        lambda_const = 2.d0*echarge*echarge*echarge*
     &    sqrt(pi/boltzmann)/boltzmann
C        _use_eq. (36) of DeWitt, 1966 (J. Math. Phys. 7, 616), but
C        _use_different definition of gamma_e, assume m_e/m_i << 1
C         so that gamma_e_e = sqrt(2) gamma_e, gamma_e_i ~ gamma_e.
C         (2 pi k/avogadro/h^2)^(-0.5) 
C           = (alpha2/avogadro^5)**(1/6)*avogdro**(3/6)
C           = alpha2**(1/6)/avogadro**(1/3)
        gamma_e_const = -3.d0*sqrt(pi)/8.d0
        dc_const = gamma_e_const*alpha2**(1.d0/6.d0)/
     &    avogadro**(1.d0/3.d0)/sqrt(4.d0*pi*electron_mass)*
     &    (4.d0*pi*lambda_const)**(1.d0/3.d0)
C         include lambda as part of definition (see below)
        dc_const = dc_const*lambda_const
        dc_const1 = 1.d0/sqrt(2.d0) - 1.d0
      endif
C       PTEH, eq. 19
      zeta2 = sum2/sum0
      lambda = lambda_const*zeta2*sqrt(sum2/(t*t*t))
!  lambda0 = -lambda/sum0
      lambda2 = 1.5d0*lambda/sum2
      lambdat = -1.5d0*lambda/t
!  lambda00 = -2.d0*lambda0/sum0
!  lambda02 = -lambda2/sum0
!  lambda0t = -lambdat/sum0
      lambda22 = 0.5d0*lambda2/sum2
      lambda2t = -1.5d0*lambda2/t
      lambdatt = -2.5d0*lambdat/t  
      if(if_dc.eq.1) then
C        I have long since disabled this option because of its known
C        invalidity for cooler temperatures and higher densities right
C        where the Coulomb interaction has its strongest effect on EOS
C        results.  Note, the implementation of this option was a quick
C        hack before I decided to drop it.  Thus, if you ever want to try
C        this option again, the implementation should be carefully checked
C        against the original literature on the diffraction correction.
C        Also, the implementation should be checked for correct
C        partial derivatives and thermodynamic consistency.
        stop 'coulomb: the diffraction correction is disabled.'
C         delta lambda due to diffraction correction
C           = dc_const*lambda*(sum2 + dc_const1*x)*x/(T * sum2^1.5d0)
C           = dc_const*lambda_const*(sum2 + dc_const1*x)/(t^2.5)*(x/sum0)
        dc_lambda = dc_const*(sum2 + dc_const1*x)/t**2.5d0*(x/sum0)
C         by definition of gamma_e
        gamma_e = dc_lambda/lambda/gamma_e_const
!    dc_lambda0 = -dc_lambda/sum0
        dc_lambda2 = dc_const/t**2.5d0*(x/sum0)
        dc_lambdat = -2.5d0*dc_lambda/t
        dc_lambdax = dc_lambda/x + dc_const*dc_const1/t**2.5d0*(x/sum0)
!    dc_lambda00 = -2.d0*dc_lambda0/sum0
!    dc_lambda02 = -dc_lambda2/sum0
!    dc_lambda0t = -2.5d0*dc_lambda0/t
!    dc_lambda0x = -dc_lambdax/sum0
        dc_lambda22 = 0.d0
        dc_lambda2t = -2.5d0*dc_lambda2/t
        dc_lambda2x = dc_lambda2/x
        dc_lambdatt = -3.5d0*dc_lambdat/t
        dc_lambdatx = -2.5d0*dc_lambdax/t
        dc_lambdaxx = 2.d0*dc_const*dc_const1/t**2.5d0/sum0
        lambda = lambda + dc_lambda
!    lambda0 = lambda0 + dc_lambda0
        lambda2 = lambda2 + dc_lambda2
        lambdat = lambdat + dc_lambdat
        lambdax = dc_lambdax
!    lambda00 = lambda00 + dc_lambda00
!    lambda02 = lambda02 + dc_lambda02
!    lambda0t = lambda0t + dc_lambda0t
!    lambda0x = dc_lambda0x
        lambda22 = lambda22 + dc_lambda22
        lambda2t = lambda2t + dc_lambda2t
        lambda2x = dc_lambda2x
        lambdatt = lambdatt + dc_lambdatt
        lambdatx = dc_lambdatx
        lambdaxx = dc_lambdaxx
      else
        lambdax = 0.d0
!    lambda0x = 0.d0
        lambda2x = 0.d0
        lambdatx = 0.d0
        lambdaxx = 0.d0
      endif
C       just before PTEH eq. 23
      gamma = (lambda*lambda/3.d0)**(1.d0/3.d0)
      xvalue = log(gamma)
      if(ifcoulomb.ge.3.and.ifcoulomb.le.9) then
C         reduces to PTEH eq. 25 for original parameters
C         dh = debye-huckel limit
        dh = gamma**1.5d0/sqrt(3.d0)
        ddh = 1.5d0*dh/gamma
        ddh2 = 0.5d0*ddh/gamma
        if(ifcoulomb.eq.5.or.ifcoulomb.eq.6) then
          if(ifcoulomb.eq.5) then
C            (For solar conditions, Gamma < 0.367 or log10(gamma) < -0.435
C            for proper Lambda and Gamma normalization.)
C            Final attempt at correction to DH so that Coulomb treatment is
C            essentially DH for solar conditions, and for larger Gammas there
C            is a rapid, but smooth transition so this results corrects DH to
C            the OCP result.  In this latter case our mean
C            Gamma is defined in such a way that we obtain a good
C            approximation to the multi-component plasma results (see paper).
            xdh = xdh10*ln10
            xmocp = xmocp10*ln10
          elseif(ifcoulomb.eq.6) then
C            alternative variation to show how sensitive results are
C            to approximation for intermediate Gamma values.
            xdh = -1.d0*ln10
            xmocp = 0.d0
          endif
          if(ifcoulomb.ne.ifcoulomb_old) then
            ifcoulomb_old = ifcoulomb
C            find modified coefficients dscmccon, dscmdcon so that modified
C            OCP result is second-order continuous (via a cubic polynomial)
C            with the DH result.
            call coulomb_adjust(xdh, xmocp,
     &        dscacon, dscbcon, dscscon, dscccon, dscdcon,
     &        dscmccon, dscmdcon)
            dx = xmocp - xdh
            ylo = 1.5d0*xdh - 0.5d0*log(3.d0)
C            dylo = 1.5d0
            dylo2 = 0.d0
            yhi = dscacon*exp(xmocp) +
     &        (dscbcon*exp(dscscon*xmocp) + dscmccon*xmocp +
     &        dscmdcon)
            dyhi = dscacon*exp(xmocp) + 
     &        (dscscon*dscbcon*exp(dscscon*xmocp) +
     &        dscmccon)
            dyhi2 = dscacon*exp(xmocp) +
     &        dscscon*dscscon*dscbcon*exp(dscscon*xmocp)
C            transform to ln yhi.
C            First do commented form using untransformed variables on RHS
C            yhi = log(yhi)
C            dyhi = dyhi/yhi
C            dyhi2 = dyhi2/yhi - dyhi*dyhi/(yhi*yhi)
            dyhi = dyhi/yhi
            dyhi2 = dyhi2/yhi - dyhi*dyhi
            yhi = log(yhi)
          endif
          if(xvalue.le.xdh) then
C            These old expressions work out to the D-H limit once
C            the derivatives transformed to Lambda derivatives, but
C            that transformation has significance loss for dgc2.
C            gc = dh
C            dgc = ddh
C            dgc2 = ddh2
C            Debye-Huckel limit with derivatives already transformed
C            to Lambda derivatives to eliminate significance loss.
            gc = lambda/3.d0
            dgc = 1.d0/3.d0
            dgc2 = 0.d0
          elseif(xvalue.lt.xmocp) then
C            cubic which by design is second-order continuous with
C            ln(DH) at xdh, and second-order continuous with ln(modified
C            ocp) at xmocp, where coulomb_adjust above makes sure
C            of that continuity by adjusting the coefficients that
C            define the modiefied ocp result.
C            acap, bcap, ccap, and dcap from numerical recipes 3.3.2, 3.3.4
            acap = (xmocp-xvalue)/dx
            bcap = (xvalue-xdh)/dx
C            ccap = acap*(acap*acap-1.d0)*dx*dx/6.d0
            ccap = (xmocp-xvalue)*(acap*acap-1.d0)*dx/6.d0
C            dcap = bcap*(bcap*bcap-1.d0)*dx*dx/6.d0
            dcap = (xvalue-xdh)*(bcap*bcap-1.d0)*dx/6.d0
C            interpolated cubic, derivative, and second derivative
C            from numerical recipes 3.3.3, 3.3.5, 3.3.6
            gc = acap*ylo + bcap*yhi + ccap*dylo2 + dcap*dyhi2
            dgc = (yhi-ylo)/dx -
     &        (3.d0*acap*acap-1.d0)*dx*dylo2/6.d0 +
     &        (3.d0*bcap*bcap-1.d0)*dx*dyhi2/6.d0
            dgc2 = acap*dylo2 + bcap*dyhi2
C            *****convert from ln gc to gc.
            gc = exp(gc)
C            commented out versions have RHS expressed in untransformed
C            variables.
C            dgc2 = exp(gc)*(dgc*dgc + dgc2)
C            the dgc on the RHS below is untransformed so this must be done
C            before dgc is transformed below.  The gc on the RHS below
C            is transformed so must be done after the gc transform above.
            dgc2 = gc*(dgc*dgc + dgc2)
C            dgc = exp(gc)*dgc
C            must be done after gc transform above.
            dgc = gc*dgc
C            *****convert from X = ln gamma independent variable to gamma.
C            second derivative done fir so can_use_untransformed
C            quantities on RHS.
            dgc2 = (dgc2-dgc)/(gamma*gamma)
            dgc = dgc/gamma
          else
C           MODIFIED one-component plasma (OCP) approximation taken from
C           DeWitt et al and modified to be second-order continuous
C           (via a cubic polynomial) with DH result.
            gc = dscacon*gamma +
     &        (dscbcon*gamma**dscscon + dscmccon*xvalue +
     &        dscmdcon)
            dgc = dscacon + 
     &        (dscscon*dscbcon*gamma**(-1.d0+dscscon) +
     &        dscmccon/gamma)
            dgc2 = 
     &        ((-1.d0+dscscon)*dscscon*dscbcon*
     &        gamma**(-2.d0+dscscon) -
     &        dscmccon/(gamma*gamma))
          endif
        elseif(ifcoulomb.eq.7.or.ifcoulomb.eq.8) then
C         _Use_same limit as indicated by extraordinarily good DH solar fit.
          if(log10(gamma).lt.-0.4d0) then
            dgc2 = ddh2
            dgc = ddh
            gc = dh
          else
C           one-component plasma (OCP) approximation taken from
C           DeWitt et al.
          gc = dscacon*gamma +
     &      (dscbcon*gamma**dscscon + dscccon*xvalue +
     &      dscdcon)
          dgc = dscacon + 
     &      (dscscon*dscbcon*gamma**(-1.d0+dscscon) +
     &      dscccon/gamma)
          dgc2 = 
     &      ((-1.d0+dscscon)*dscscon*dscbcon*
     &      gamma**(-2.d0+dscscon) -
     &      dscccon/(gamma*gamma))
          endif
        elseif(ifcoulomb.eq.3.or.ifcoulomb.eq.4.or.ifcoulomb.eq.9)
     &      then
C           ocp = simplified one-component plasma approximation taken from
C           pteh (works in conjunction with smooth transition from dh to
C           provide what is actually quite good approximation for ocp
C           above gamma = 1.
          ocp = acon*gamma + dcon
          docp = acon
          docp2 = 0.d0
C         smooth transition between dh and ocp
          gc = (dh**gcon + ocp**gcon)**(1.d0/gcon)
          dgc = (ddh*dh**(-1.d0+gcon) + docp*ocp**(-1.d0+gcon))*
     &      (dh**gcon + ocp**gcon)**(-1.d0 + 1.d0/gcon)
          dgc2 = (ddh2*dh**(-1.d0+gcon) + docp2*ocp**(-1.d0+gcon) +
     &      ddh*ddh*(-1.d0+gcon)*dh**(-2.d0+gcon) +
     &      docp*docp*(-1.d0+gcon)*ocp**(-2.d0+gcon))*
     &      (dh**gcon + ocp**gcon)**(-1.d0 + 1.d0/gcon) +
     &      (1.d0 - gcon)*
     &      (ddh*dh**(-1.d0+gcon) + docp*ocp**(-1.d0+gcon))*
     &      (ddh*dh**(-1.d0+gcon) + docp*ocp**(-1.d0+gcon))*
     &      (dh**gcon + ocp**gcon)**(-2.d0 + 1.d0/gcon)
        else
          stop 'coulomb: logic error'
        endif
C        convert derivatives to lambda except for case where that has
C        been done already.
        if(.not.((ifcoulomb.eq.5.or.ifcoulomb.eq.6).and.
     &      xvalue.le.xdh)) then
          dgc = dgc*(2.d0/3.d0)*gamma/lambda
          dgc2 = dgc2*(4.d0/9.d0)*(gamma/lambda)*(gamma/lambda) -
     &      dgc/(3.d0*lambda)
        endif
      elseif(ifcoulomb.eq.1.or.ifcoulomb.eq.2) then
C         Debye-Huckel limit
        gc = lambda/3.d0
        dgc = 1.d0/3.d0
        dgc2 = 0.d0
      else
        stop 'coulomb: ifcoulomb must be in range from 1-9'
      endif
      fcoulomb = -sum0*boltzmann*t*gc
      if(ifcoulomb.eq.1.or.ifcoulomb.eq.2) then
C         in this case gc proportional lambda proportional to sum0^-1 and
C         fcoulomb is completely independent of sum0.
        fcoulomb0 = 0.d0
      else
        fcoulomb0 = -boltzmann*t*(gc - dgc*lambda)
      endif
      fcoulomb2 = fcoulomb*dgc/gc*lambda2
      fcoulombt = fcoulomb*(1.d0/t + dgc/gc*lambdat)
      fcoulombx = fcoulomb*dgc/gc*lambdax
C      _use_the fact that lambda proportional to sum0^-1 
C       (diffraction corrected or not)
      fcoulomb00 = -boltzmann*t*dgc2*lambda*lambda/sum0
      fcoulomb02 = -boltzmann*t*(-dgc2*lambda2*lambda)
      fcoulomb0t = fcoulomb0/t - boltzmann*t*(-dgc2*lambdat*lambda)
      fcoulomb0x = -boltzmann*t*(-dgc2*lambdax*lambda)
      fcoulomb22 = fcoulomb2*(dgc2*lambda2/dgc + lambda22/lambda2)
      fcoulomb2t = fcoulomb2*
     &  (1.d0/t + dgc2*lambdat/dgc + lambda2t/lambda2)
      fcoulomb2x = fcoulomb2*(dgc2*lambdax/dgc + lambda2x/lambda2)
      fcoulombtt = -sum0*boltzmann*(2.d0*dgc*lambdat +
     &  t*(dgc2*lambdat*lambdat + dgc*lambdatt))
      fcoulombtx = -sum0*boltzmann*(dgc*lambdax +
     &  t*(dgc2*lambdat*lambdax + dgc*lambdatx))
      fcoulombxx = -sum0*boltzmann*t*
     &  (dgc2*lambdax*lambdax + dgc*lambdaxx)
      if(ifcoulomb.eq.2) then
C         prior to this, ifcoulomb = 1 and 2 logic is identical.
C         with this flag on correct DH results by tau(x) correction
C         also note that prior x derivatives are zero because if_dc
C         must be 0 if ifcoulomb = 2.
        call tau_calc(sum0, sum2, t, x, tau, dtau, d2tau)
C         second partials first to_use_untransformed first partials
C         and fcoulomb
        fcoulomb00 = fcoulomb*d2tau(1,1)
        fcoulomb02 = fcoulomb2*dtau(1) + fcoulomb*d2tau(1,2)
        fcoulomb0t = fcoulombt*dtau(1) + fcoulomb*d2tau(1,3)
        fcoulomb0x = fcoulomb*d2tau(1,4)
        fcoulomb22 = fcoulomb22*tau + 2.d0*fcoulomb2*dtau(2) +
     &    fcoulomb*d2tau(2,2)
        fcoulomb2t = fcoulomb2t*tau + fcoulomb2*dtau(3) +
     &    fcoulombt*dtau(2) + fcoulomb*d2tau(2,3)
        fcoulomb2x = fcoulomb2*dtau(4) + fcoulomb*d2tau(2,4)
        fcoulombtt = fcoulombtt*tau + 2.d0*fcoulombt*dtau(3) +
     &    fcoulomb*d2tau(3,3)
        fcoulombtx = fcoulombt*dtau(4) + fcoulomb*d2tau(3,4)
        fcoulombxx = fcoulomb*d2tau(4,4)
C         first partials next to_use_untransformed fcoulomb
C         original fcoulomb0 and fcoulombx are zero in DH limit.
        fcoulomb0 = fcoulomb*dtau(1)
        fcoulomb2 = fcoulomb2*tau + fcoulomb*dtau(2)
        fcoulombt = fcoulombt*tau + fcoulomb*dtau(3)
        fcoulombx = fcoulomb*dtau(4)
        fcoulomb = fcoulomb*tau
      endif
      end
