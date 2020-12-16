C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: master_coulomb.f 370 2006-11-29 23:57:39Z airwin $
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
      subroutine master_coulomb(ifnr_in, rhostar, f,
     &  sum0, sum0ne, sum0f, sum0t, sum2, sum2ne, sum2f, sum2t,
     &  n_e, t, p_e, pstar, lambda, gamma_e,
     &  ifcoulomb, if_dc, if_pteh,
     &  dve, dvef, dvet, dve0, dve2,
     &  dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22)
C       input quantities:
C       ifnr_in = 0, simple iteration, calculate dve, dv0, dv2, and all f and
C         t derivatives with no other variable fixed.
C       ifnr_in = 1, NR iteration, calculate dve, dv0, dv2, and all sum0 and
C         sum2 derivatives.
C       ifnr_in = 2, calculate dve, dv0, dv2, and all f and
C         t derivatives *keeping* sum0 and sum2 fixed.
C       ifnr_in = 3, is combination of ifnr_in = 1 and 2.
C       rhostar(9), 1=re proportional to n_e,
C       the rest are partials of ln n_e wrt fl and tl,
C       2 = f, 3 = t, 4 = ff, 5 = ft, 6 = tt, 7 = fff, 8 = fft 9 = ftt
C       sum0 = sum_ions nion, where nion is the number per unit volume,
C         and the sum taken over positive ions, but excluding n_e.
C       sum2 = sum_ions Z^2 nion,
C         Z is the charge on the ions.
C       outside this routine, sum0 and sum2 are calculated as a function
C         of n_e *or* f and T, but for programming convenience sum0
C         and sum2 are viewed as functions of n_e, f, and T.
C         The derivatives calculated outside are:
C       sum0ne = partial sum0(n_e, f, T) wrt n_e
C       sum0f = partial sum0(n_e, f, T) wrt fl
C       sum0t = partial sum0(n_e, f, T) wrt tl
C       sum2ne = partial sum2(n_e, f, T) wrt n_e
C       sum2f = partial sum2(n_e, f, T) wrt fl
C       sum2t = partial sum2(n_e, f, T) wrt tl
C       Since the actual functional dependence is n_e *or* f, T, it is the
C         calling routine's responsibility to calculate one set of 
C         derivatives and set the the other set to zero as appropriate.
C       n_e = electron number density per unit volume
C       t = temperature (Kelvin)
C       p_e = electron pressure (cgs)
C       pstar(9), 1=proportional to p_e,
C       the rest are partials of ln p_e wrt fl and tl,
C       2 = f, 3 = t, 41 = ff, 5 = ft, 6 = tt, 7 = fff, 8 = fft 9 = ftt
C       ifcoulomb = 0 ignore the effects of the Coulomb interaction.
C       ifcoulomb = 1 treat Coulomb interaction in the Debye-Huckel
C         approximation.
C       ifcoulomb = 2 treat Coulomb interaction in the Debye-Huckel 
C         approximation corrected by tau(x).
C       ifcoulomb = 3 treat Coulomb interaction in the PTEH approximation 
C         with pteh theta_e.
C       ifcoulomb = 4 treat Coulomb interaction in the PTEH approximation
C         with fermi-dirac theta_e.
C       ifcoulomb = 9 treat Coulomb interaction in the PTEH approximation
C         with fermi-dirac theta_e and DeWitt definition of lambda
C         (using sum0a = sum0 + ne*theta_e).
C       ifcoulomb = 5 DH smoothly connected to modified OCP and DeWitt
C         definition of lambda.
C       ifcouloumb = 6 same as 5 with alternative smooth connection
C       ifcoulomb = 7, DH (Gamma < 1) or OCP using new DeWitt lambda.
C       ifcoulomb = 8, same as 7 with theta_e = 0.
C       if_dc = 0, ignore effects of the diffraction correction
C       if_dc = 1, apply diffraction correction to Lambda
C       if_pteh = 0, sum0 and sum2 are functions of n_i only
C       if_pteh = 1, sum0 and sum2 are functions of n_e only (PTEH
C         approximation)
C       output quantities:
C       lambda = plasma interaction parameter (defined by PTEH *or* defined
C       by DeWitt if ifcoulomb > 4 and consistent
C       with the simpler CG 15.60 criterion for single component plasmas.)
C       n.b., if_dc = 1 means this is the diffraction-corrected quantity.
C       gamma_e =  DeWitt paper II gamma_e = diffraction correction parameter.
C       dve, dvef, dvet, dve0, dve2:
C       dv0, dv0f, dv0t, dv00, dv02:
C       dv2, dv2f, dv2t, dv22:
C       N.B. the f, t, 0, and 2 suffixes for dve, dv0, and dv2 type of
C       output variables refer to fl, tl, sum0, and sum2 partial derivatives.
C       dve is the change in the electron component of the equilibrium
C         constant.
C       dv0 is the change in the sum0 component of the equilibrium
C         constant.
C       dv2 is the change in the sum2 component of the equilibrium
C         constant.
      implicit none
      include 'constants.h'
      integer ifnr_in, ifnr, ifcoulomb, if_dc, if_pteh
      double precision
     &  rhostar(9), f, dpsidf, dpsidff, dpsidfff,
     &  sum0, sum0ne, sum0f, sum0t, sum2, sum2ne, sum2f, sum2t,
     &  n_e, t, p_e, pstar(9), lambda, gamma_e,
     &  dve, dvef, dvet, dve0, dve2,
     &  dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22
      double precision free_coulomb, free_coulombf, free_coulomb0, free_coulomb2
      double precision
     &  theta_e, theta_ef, theta_et, theta_eff, theta_eft, theta_ett,
     &  thetane, thetanef, thetanet, thetat, thetatf, thetatt,
     &  theta_pow,
     &  sum0a, sum0ane, sum0anef, sum0anet,
     &  sum0atl, sum0atlf, sum0atlt, sum0af,
     &  sum2a, sum2ane, sum2anef, sum2anet,
     &  sum2atl, sum2atlf, sum2atlt, sum2af,
     &  dmucoulomb, dmucoulombf, dmucoulombt,
     &  dmucoulomb0, dmucoulomb2,
     &  dvcoulomb, dvcoulombf, dvcoulombt,
     &  fcoulomb, fcoulomb0, fcoulomb2, fcoulombt, fcoulombx,
     &  fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0x, 
     &  fcoulomb22, fcoulomb2t, fcoulomb2x, fcoulombtt, fcoulombtx,
     &  fcoulombxx, fcft0f, fcft2f, fcft0t, fcft2t, fcfttf, fcfttt,
     &  dpcoulomb, dpcoulombt, dpcoulombf, dpcoulomb0, dpcoulomb2,
     &  dscoulomb, dscoulombf, dscoulombt, ducoulomb,
     &  dscoulomb_local,
     &  thetax, fcftxf, fcftxt, x, xf, xt, xe, xef, xet,
     &  xfixedt, xfixedtf, xfixedtt
      logical ifnr13, ifnr23
      save
      if(ifcoulomb.ge.1) then
        if(ifnr_in.gt.3)
     &    stop 'master_coulomb: bad ifnr_in'
C        must zero NR derivative calculations when if_pteh = 1
        if(if_pteh.eq.1) then
          ifnr = min(ifnr_in,0)
          dve0 = 0.d0
          dve2 = 0.d0
          dv00 = 0.d0
          dv02 = 0.d0
          dv22 = 0.d0
        else
          ifnr = ifnr_in
        endif
        ifnr13 = ifnr.eq.1.or.ifnr.eq.3
        ifnr23 = ifnr.eq.2.or.ifnr.eq.3
C        don't want tau(x) correction and diffraction correction
C        simultaneously
        if(ifcoulomb.eq.2.and.if_dc.eq.1) stop
     &    'master_coulomb: bad combination of ifcoulomb and if_dc'
C        set up coulomb interaction treatment
C        partial psi/d ln f and higher derivatives
        dpsidf = sqrt(1.d0+f)
        dpsidff = 0.5d0*f*dpsidf/(1.d0+f)
        dpsidfff = dpsidff*(1.d0 - 0.5d0*f/(1.d0+f))
C        Decide power of original theta_e
        if(ifcoulomb.gt.5) then
C          Increase this to get more rapid quenching of electronic
C          component as the degeneracy is increased.
C          Currently experimental, but at some point this will also
C          be used for standard Coulomb treament (ifcoulomb = 5).
!          theta_pow = 1.3d0
          theta_pow = 1.d0
        else
C          Traditional value for MDH, etc.
          theta_pow = 1.d0
        endif
C         partial ln n_e(T, psi)/partial psi
        if(ifcoulomb.ne.3.and.ifcoulomb.ne.8) then
          theta_e = rhostar(2)/dpsidf
          theta_ef = (rhostar(4)-rhostar(2)*dpsidff/dpsidf)/dpsidf
          theta_et = rhostar(5)/dpsidf
          theta_eff = (rhostar(7)-(2.d0*rhostar(4)*dpsidff +
     &      rhostar(2)*(dpsidfff-2.d0*dpsidff/dpsidf*dpsidff))/
     &      dpsidf)/dpsidf
          theta_eft = (rhostar(8)-rhostar(5)*dpsidff/dpsidf)/dpsidf
          theta_ett = rhostar(9)/dpsidf
          if(theta_pow.ne.1.d0) then
C            transform to power of original theta_e
            theta_eff =
     &        theta_eff*theta_pow*theta_e**(theta_pow-1.d0) +
     &        theta_ef*theta_ef*theta_pow*(theta_pow-1.d0)*
     &        theta_e**(theta_pow-2.d0)
            theta_eft =
     &        theta_eft*theta_pow*theta_e**(theta_pow-1.d0) +
     &        theta_ef*theta_et*theta_pow*(theta_pow-1.d0)*
     &        theta_e**(theta_pow-2.d0)
            theta_ett =
     &        theta_ett*theta_pow*theta_e**(theta_pow-1.d0) +
     &        theta_et*theta_et*theta_pow*(theta_pow-1.d0)*
     &        theta_e**(theta_pow-2.d0)
            theta_ef = theta_ef*theta_pow*theta_e**(theta_pow-1.d0)
            theta_et = theta_et*theta_pow*theta_e**(theta_pow-1.d0)
            theta_e = theta_e**theta_pow
          endif
C          n.b. partial theta_e(f(ne,T), T)/partial ln ne =
C          theta_ef/rhostar(2)
          thetane = theta_ef/rhostar(2)
          thetanef = (theta_eff - theta_ef*rhostar(4)/
     &      rhostar(2))/rhostar(2)
          thetanet = (theta_eft - theta_ef*rhostar(5)/
     &      rhostar(2))/rhostar(2)
C          n.b. partial theta_e(ne, T)/partial ln T =
C          theta_et - theta_ef*rhostar(3)/rhostar(2)
          thetat = theta_et - theta_ef*rhostar(3)/rhostar(2)
          thetatf = theta_eft - theta_eff*rhostar(3)/rhostar(2) -
     &      theta_ef*(rhostar(5) - rhostar(3)*rhostar(4)/rhostar(2))/
     &      rhostar(2)
          thetatt = theta_ett - theta_eft*rhostar(3)/rhostar(2) -
     &      theta_ef*(rhostar(6) - rhostar(3)*rhostar(5)/rhostar(2))/
     &      rhostar(2)
        elseif(ifcoulomb.eq.3) then
C         _use_simplified pteh approximation for theta_e and
C          all subsequent derivatives based on their eq. 29.
          call pteh_theta(cd*rhostar(1), t, rhostar,
     &      theta_e, theta_ef, theta_et,
     &      theta_eff, theta_eft, theta_ett,
     &      thetane, thetanef, thetanet,
     &      thetat, thetatf, thetatt)
        elseif(ifcoulomb.eq.8) then
          theta_e = 0.d0
          theta_ef = 0.d0
          theta_et = 0.d0
          theta_eff = 0.d0
          theta_eft = 0.d0
          theta_ett = 0.d0
          thetane = 0.d0
          thetanef = 0.d0
          thetanet = 0.d0
          thetat = 0.d0
          thetatf = 0.d0
          thetatf = 0.d0
        else
          stop 'master_coulomb: should not happen'
        endif
        if(ifcoulomb.gt.4) then
C          DeWitt definition
          sum0a = sum0 + theta_e*n_e
          sum0af = theta_ef*n_e + theta_e*n_e*rhostar(2)
        else
C          PTEH definition
          sum0a = sum0
          sum0af = 0.d0
        endif
        sum2a = sum2 + theta_e*n_e
        sum2af = theta_ef*n_e + theta_e*n_e*rhostar(2)
C        sum0af and sum2af are the partials of sum?a wrt fl for fixed
C        auxiliary variables.  Ordinarily sum0 and sum2 *are* auxiliary
C        variables (or are transformed later to a different set of
C        auxiliary variables) and we are done.  However, for special case
C        of if_pteh = 1, sum0 and sum2 are not auxiliary variables, and
C        we must do more.
        if(if_pteh.eq.1) then
          sum0af = sum0af + sum0ne*rhostar(2)
          sum2af = sum2af + sum2ne*rhostar(2)
        endif
C        n.b. x and its derivatives and the partials of fcoulomb
C        wrt x returned by coulomb below *only* relevant when
C        if_dc.eq.1 or ifcoulomb.eq.2.
C        n.b.
C        xf = partial x(f,T)/partial fl
C        xt = partial x(f,T)/partial ft
C        xe = partial x(n_e,T)/partial n_e
C        xef = partial xe(f,T)/partial fl
C        xet = partial xe(f,T)/partial tl
        if(if_dc.eq.1) then
          x = theta_e*n_e
          xf = theta_ef*n_e + x*rhostar(2)
          xt = theta_et*n_e + x*rhostar(3)
          xe = theta_e + thetane
          xef = theta_ef + thetanef
          xet = theta_et + thetanet
        elseif(ifcoulomb.eq.2) then
C          thetax = (n_e k T/p_e) ~ 1
          thetax = (n_e*boltzmann*t/p_e)
          x = thetax*n_e
          xf = x*(2.d0*rhostar(2) - pstar(2))
          xt = x*(2.d0*rhostar(3) - pstar(3) + 1.d0)
          xe = thetax*(2.d0 - pstar(2)/rhostar(2))
          xef = thetax*((2.d0 - pstar(2)/rhostar(2))*
     &      (rhostar(2)-pstar(2)) +
     &      (-pstar(4) + pstar(2)*rhostar(4)/rhostar(2))/rhostar(2))
          xet = thetax*((2.d0 - pstar(2)/rhostar(2))*
     &      (rhostar(3) - pstar(3) + 1.d0) +
     &      (-pstar(5) + pstar(2)*rhostar(5)/rhostar(2))/rhostar(2))
        else
C          should define since used below (although x derivatives
C          set to zero by coulomb).
          x = 0.d0
          xf = 0.d0
          xt = 0.d0
          xe = 0.d0
          xef = 0.d0
          xet = 0.d0
        endif
        call coulomb(ifcoulomb, if_dc, sum0a, sum2a, t, x,
     &    fcoulomb, fcoulomb0, fcoulomb2, fcoulombt, fcoulombx,
     &    fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0x,
     &    fcoulomb22, fcoulomb2t, fcoulomb2x,
     &    fcoulombtt, fcoulombtx, fcoulombxx, lambda, gamma_e)
C        fcoulomb (free energy per unit volume) is a function of all 
C        n_i (where i ranges over positive ions), n_e, T.  
C        Calculate additions to chemical potential of electrons,
C        partial (V*fcoulomb(n_i, n_e, T)) wrt n_e * partial n_e wrt (V n_e)
C        = partial fcoulomb(sum0a(n_i,n_e), sum2a(n_i,n_e), n_e, T) 
C          wrt n_e
C        n.b. sum0 and sum2 are functions of n_e alone or are 
C          independent of n_e.  (depending on ifcoulomb for
C          sum0) sum0a and sum2a have the extra factor
C          theta_e(f(n_e,T),T)*n_e
        if(ifcoulomb.gt.4) then
C          DeWitt definition
C          sum0a = sum0 + theta_e*n_e
C          calculate partial sum0a(n_e, n_i, T) wrt n_e
C          n.b., when n_i is fixed, sum0 is fixed except in pteh
C          approximation when it varies with n_e.
          sum0ane = sum0ne + theta_e + thetane
C          calculate partial sum0a(n_e, n_i, T) wrt ln T
C          (ignore sum0 dependence on t, since fixed n_i means sum0 is fixed
C          or depends on fixed n_e in initial approximation)
          sum0atl = n_e*thetat
        else
C          PTEH definition
C          sum0a = sum0
          sum0ane = sum0ne
          sum0atl = 0.d0
        endif
C        calculate partial sum2a(n_e, n_i, T) wrt n_e
C        n.b., when n_i is fixed, sum2 is fixed except in pteh
C        approximation when it varies with n_e.
        sum2ane = sum2ne + theta_e + thetane
C        calculate partial sum2a(n_e, n_i, T) wrt ln T
C        (ignore sum2 dependence on t, since fixed n_i means sum2 is fixed
C        or depends on fixed n_e in initial approximation)
        sum2atl = n_e*thetat
C        partial fcoulomb(sum0, sum2, t, ne) wrt ne
        dmucoulomb = fcoulomb0*sum0ane + fcoulomb2*sum2ane +
     &    fcoulombx*xe
C        partials of dmucoulomb wrt sum0 and sum2
C        note that sum[02]ane and xe are completely independent of sum[02].
        dmucoulomb0 = fcoulomb00*sum0ane + fcoulomb02*sum2ane +
     &    fcoulomb0x*xe
        dmucoulomb2 = fcoulomb02*sum0ane + fcoulomb22*sum2ane +
     &    fcoulomb2x*xe
C        change in equilibrium constant related to negative chemical
C        potential/kT.
        dvcoulomb = -dmucoulomb/(boltzmann*t)
        if(ifnr.eq.0.or.ifnr23) then
C          following partial derivatives of fcoulomb0, fcoulomb2, and 
C          fcoulombx are wrt ln f, ln t holding nothing else fixed
C          if(ifnr.eq.0) or holding input sum0 and sum2 fixed (ifnr23)
C          f derivatives first....
          fcft0f = 
     &      n_e*rhostar(2)*(fcoulomb00*sum0ane +
     &      fcoulomb02*sum2ane) + fcoulomb0x*xe*n_e*rhostar(2)
          fcft2f = 
     &      n_e*rhostar(2)*(fcoulomb02*sum0ane +
     &      fcoulomb22*sum2ane) + fcoulomb2x*xe*n_e*rhostar(2)
          fcftxf = 
     &      n_e*rhostar(2)*(fcoulomb0x*sum0ane +
     &      fcoulomb2x*sum2ane) + fcoulombxx*xe*n_e*rhostar(2)
C          we recognize that mixed second partial derivatives of (input) 
C          sum0 and sum2 wrt ne, f = 0.
          if(ifcoulomb.gt.4) then
C            DeWitt definition
C            sum0a = sum0 + theta_e*n_e
            sum0anef = theta_ef + thetanef
          else
C            PTEH definition
C            sum0a = sum0
            sum0anef = 0.d0
          endif
          sum2anef = theta_ef + thetanef
C          now t derivatives....
          fcft0t =
     &      t*fcoulomb0t +
     &      fcoulomb00*sum0atl + fcoulomb02*sum2atl +
     &      n_e*rhostar(3)*(fcoulomb00*sum0ane + fcoulomb02*sum2ane) +
     &      fcoulomb0x*xt
          fcft2t =
     &      t*fcoulomb2t +
     &      fcoulomb02*sum0atl + fcoulomb22*sum2atl +
     &      n_e*rhostar(3)*(fcoulomb02*sum0ane + fcoulomb22*sum2ane) +
     &      fcoulomb2x*xt
          fcftxt =
     &      t*fcoulombtx +
     &      fcoulomb0x*sum0atl + fcoulomb2x*sum2atl +
     &      n_e*rhostar(3)*(fcoulomb0x*sum0ane + fcoulomb2x*sum2ane) +
     &      fcoulombxx*xt
C          we recognize that the second mixed partial derivatives of 
C          (input) sum0 and sum2 wrt ne, t = 0.
          if(ifcoulomb.gt.4) then
C            DeWitt definition
C            sum0a = sum0 + theta_e*n_e
            sum0anet = theta_et + thetanet
          else
C            PTEH definition
C            sum0a = sum0
            sum0anet = 0.d0
          endif
          sum2anet = theta_et + thetanet
          if(ifnr.eq.0) then
C            account for sum0, sum2 dependence on f, t
            fcft0f = fcft0f + fcoulomb00*sum0f + fcoulomb02*sum2f
            fcft2f = fcft2f + fcoulomb02*sum0f + fcoulomb22*sum2f
            fcftxf = fcftxf + fcoulomb0x*sum0f + fcoulomb2x*sum2f
            fcft0t = fcft0t + fcoulomb00*sum0t + fcoulomb02*sum2t
            fcft2t = fcft2t + fcoulomb02*sum0t + fcoulomb22*sum2t
            fcftxt = fcftxt + fcoulomb0x*sum0t + fcoulomb2x*sum2t
          endif
C         dmucoulomb = fcoulomb0*sum0ane + fcoulomb2*sum2ane +
C     &     fcoulombx*xe
          dmucoulombf = fcft0f*sum0ane + fcft2f*sum2ane +
     &      fcoulomb0*sum0anef + fcoulomb2*sum2anef +
     &      fcftxf*xe + fcoulombx*xef
          dvcoulombf = -dmucoulombf/(boltzmann*t)
          dmucoulombt = fcft0t*sum0ane + fcft2t*sum2ane +
     &      fcoulomb0*sum0anet + fcoulomb2*sum2anet +
     &      fcftxt*xe + fcoulombx*xet
          dvcoulombt = (dmucoulomb - dmucoulombt)/(boltzmann*t)
        endif
C        calculate negative of electron component (dvcoulomb)
C        appropriately....
        dve = dvcoulomb
        if(ifnr.eq.0) then
C          derivative wrt f, t holding no other variable fixed.
          dvef = dvcoulombf
          dvet = dvcoulombt
        elseif(ifnr23) then
C          derivative wrt f, t holding sum0 and sum2 fixed.
          dvef = dvcoulombf
          dvet = dvcoulombt
        endif
        if(ifnr13) then
C          derivative wrt sum0, sum2, holding f, t fixed.
          dve0 = -dmucoulomb0/(boltzmann*t)
          dve2 = -dmucoulomb2/(boltzmann*t)
        endif
C        remaining equilibrium constant changes are absolute,
C        *not* relative to next lower ion.
        if(if_pteh.eq.1) then
C          in this approximation, sum0 and sum2 are independent
C          of n_i.
          dv0 = 0.d0
          dv0f = 0.d0
          dv0t = 0.d0
          dv2 = 0.d0
          dv2f = 0.d0
          dv2t = 0.d0
        else
C          calculate dv0 = negative of chemical potential/kT due to
C          sum0 term = -partial f(sum0, sum2, T, x(f,T))/partial sum0 *
C            partial sum0/partial n_i /kT
          dv0 = -fcoulomb0/(boltzmann*t)
C          calculate dv2 = negative of chemical potential/kT/Z^2 due to
C          sum2 term = -partial f(sum0, sum2, T, x(f,T))/partial sum2 *
C            partial sum2/partial n_i /kT
          dv2 = -fcoulomb2/(boltzmann*t)
          if(ifnr.eq.0) then
C            derivative wrt f, t holding no other variable fixed.
            dv0f = -fcft0f/(boltzmann*t)
            dv0t = -fcft0t/(boltzmann*t) -dv0
            dv2f = -fcft2f/(boltzmann*t)
            dv2t = -fcft2t/(boltzmann*t) -dv2
          elseif(ifnr23) then
C            derivative wrt f, t holding sum0 and sum2 fixed.
            dv0f = -fcft0f/(boltzmann*t)
            dv0t = -fcft0t/(boltzmann*t) -dv0
            dv2f = -fcft2f/(boltzmann*t)
            dv2t = -fcft2t/(boltzmann*t) -dv2
          endif
          if(ifnr13) then
C            derivative wrt sum0, sum2, holding f, t fixed.
            dv00 = -fcoulomb00/(boltzmann*t)
            dv02 = -fcoulomb02/(boltzmann*t)
            dv22 = -fcoulomb22/(boltzmann*t)
          endif
        endif
      else
        dve = 0.d0
        dvef = 0.d0
        dvet = 0.d0
        dv0 = 0.d0
        dv0f = 0.d0
        dv0t = 0.d0
        dv2 = 0.d0
        dv2f = 0.d0
        dv2t = 0.d0
        lambda = 0.d0
      endif
      return
      entry master_coulomb_pressure(
     &  ifcoulomb, if_dc, if_pteh, sum0, sum2,
     &  t, n_e, rhostar, pstar,
     &  dpcoulomb, dpcoulombf, dpcoulombt, dpcoulomb0, dpcoulomb2)
      if(ifcoulomb.ge.1) then
C        don't want tau(x) correction and diffraction correction
C        simultaneously
        if(ifcoulomb.eq.2.and.if_dc.eq.1) stop
     &    'master_coulomb_pressure: bad ifcoulomb or if_dc'
C        calculate change to thermodynamic quantities 
C        n.b. fcoulomb is per unit volume so is dscoulomb, ducoulomb
C        n.b. fcoulomb is considered to be a function of n_e, n_i, T
C        because sum0a(n_e, n_i), sum2a(n_e, n_i, T), and x(n_e, T),
C        where x = thetax*n_e (ifcoulomb=2) or theta_e*n_e (if_dc = 1)
C        partial -fcoulomb(n_i, n_e, t) wrt t
C        n.b.
C          xfixedt = partial x(ne,T)/partial T
C          also note that 
C          partial fl(ne, T)/partial ln T = -rhostar(3)/rhostar(2)
        if(if_dc.eq.1) then
          xfixedt = (-theta_ef*rhostar(3)/rhostar(2))*n_e/t
        elseif(ifcoulomb.eq.2) then
          xfixedt = x*(1.d0 + pstar(2)*rhostar(3)/rhostar(2) -
     &      pstar(3))/t
        else
          xfixedt = 0.d0
        endif
C        non-returned version needed for this entry
        dscoulomb_local = -fcoulombt -
     &    fcoulomb0*sum0atl/t - fcoulomb2*sum2atl/t - 
     &    fcoulombx*xfixedt
C        dpcoulomb = - partial(V*fcoulomb(n_e, n_i, T)) wrt V
C          = -fcoulomb + n_e*partial fcoulomb wrt n_e + sum_i n_i*
C          partial fcoulomb wrt n_i
C          = -fcoulomb + n_e*dmucoulomb + sum0*partial fcoulomb
C          wrt sum0a + sum2*partial fcoulomb wrt sum2a
C        n.b. the latter two terms are dropped if sum0 and sum2 
C          depend only on n_e.  This occurs if
C          the PTEH approximation being used for sum0, sum2 
        if(if_pteh.eq.1) then
          dpcoulomb = -fcoulomb + n_e*dmucoulomb
C          dpcoulomb(f,t)/ ln f and ln t:  to do these derivatives
C          fcoulomb (for the fully ionized case) is considered 
C          to be a function of (sum0(n_e), sum2(n_e)+
C          theta_e(n_e,t)*n_e, t, x(n_e,T)) = function(n_e, T),
C          completely independent of n_i.
C          also recall dmucoulomb = partial fcoulomb(n_e, T) wrt n_e
C          in the fully ionized case.
C          There is a useful cancellation of terms when using 
C          this viewpoint....
          dpcoulombf = n_e*dmucoulombf
C          dscoulomb = partial -fcoulomb(n_e, t) wrt t
          dpcoulombt = t*dscoulomb_local + n_e*dmucoulombt
          dpcoulomb0 = 0.d0
          dpcoulomb2 = 0.d0
        else
          dpcoulomb = -fcoulomb + n_e*dmucoulomb +
     &      sum0*fcoulomb0 + sum2*fcoulomb2
C          dpcoulomb(f,t)/ ln f and ln t:  to do these derivatives
C          fcoulomb is considered a function of (sum0(f,t), sum2(f,t)+
C          theta_e(n_e,t)*n_e, t, x(n_e,T)).  In other words, the previous
C          n_i dependence is replaced with sum0 and sum2 dependence.
C          There is a useful cancellation of terms when using 
C          this viewpoint....
          dpcoulombf = n_e*dmucoulombf +
     &      sum0*fcft0f + sum2*fcft2f
C          dscoulomb = partial -fcoulomb(sum0, sum2, n_e, t) wrt t
          dpcoulombt = t*dscoulomb_local + n_e*dmucoulombt +
     &      sum0*fcft0t + sum2*fcft2t
          dpcoulomb0 = n_e*dmucoulomb0 +
     &      sum0*fcoulomb00 + sum2*fcoulomb02
          dpcoulomb2 = n_e*dmucoulomb2 +
     &      sum0*fcoulomb02 + sum2*fcoulomb22
        endif
      else
        dpcoulomb = 0.d0
        dpcoulombf = 0.d0
        dpcoulombt = 0.d0
        dpcoulomb0 = 0.d0
        dpcoulomb2 = 0.d0
      endif
      return
      entry master_coulomb_free(if_pteh,
     &  free_coulomb, free_coulombf, free_coulomb0, free_coulomb2)
C      return free energy per unit volume and its derivatives with
C      respect to fl and the auxiliary variables, sum0, and sum2.
      free_coulomb = fcoulomb
      free_coulombf = fcoulomb0*sum0af + fcoulomb2*sum2af +
     &  fcoulombx*xf
C      n.b. The fcoulomb0 and fcoulomb2 derivatives are actually with
C      respect to sum0a and sum2a with x = thetax*n_e fixed (if there is a
C      tau correction), but with f, and t fixed, thetax*n_e is fixed
C      and also the partials of sum0a with respect to sum0 and sum2a wrt
C      sum2 are both unity.
      if(if_pteh.eq.1) then
C        if the pteh Coulomb approximation is used, then sum0 and
C        sum2 are just functions of f and t.
        free_coulomb0 = 0.d0
        free_coulomb2 = 0.d0
      else
        free_coulomb0 = fcoulomb0
        free_coulomb2 = fcoulomb2
      endif
      return
      entry master_coulomb_end(rhostar,
     &  sum0, sum0ne, sum0f, sum0t,
     &  sum2, sum2ne, sum2f, sum2t,
     &  n_e, t, pstar,
     &  ifcoulomb, if_dc, if_pteh,
     &  dpcoulomb, dpcoulombf, dpcoulombt,
     &  dscoulomb, dscoulombf, dscoulombt, ducoulomb)
      if(ifcoulomb.ge.1) then
C        don't want tau(x) correction and diffraction correction
C        simultaneously
        if(ifcoulomb.eq.2.and.if_dc.eq.1) stop
     &    'master_coulomb_end: bad combination of ifcoulomb and if_dc'
C        calculate change to thermodynamic quantities 
C        n.b. fcoulomb is per unit volume so is dscoulomb, ducoulomb
C        n.b. fcoulomb is considered to be a function of n_e, n_i, T
C        because sum0a(n_e, n_i), sum2a(n_e, n_i, T), and x(n_e, T),
C        where x = thetax*n_e (ifcoulomb=2) or theta_e*n_e (if_dc = 1)
C        partial -fcoulomb(n_i, n_e, t) wrt t
C        n.b.
C          xfixedt = partial x(ne,T)/partial T
C          also note that 
C          partial fl(ne, T)/partial ln T = -rhostar(3)/rhostar(2)
        if(if_dc.eq.1) then
          xfixedt = (-theta_ef*rhostar(3)/rhostar(2))*n_e/t
          xfixedtf = xfixedt*(theta_eff/theta_ef +
     &      rhostar(5)/rhostar(3) - rhostar(4)/rhostar(2) +
     &      rhostar(2))
          xfixedtt = xfixedt*(theta_eft/theta_ef +
     &      rhostar(6)/rhostar(3) - rhostar(5)/rhostar(2) +
     &      rhostar(3) - 1.d0)
        elseif(ifcoulomb.eq.2) then
          xfixedt = x*(1.d0 + pstar(2)*rhostar(3)/rhostar(2) -
     &      pstar(3))/t
          xfixedtf = x*((2.d0*rhostar(2) - pstar(2))*
     &      (1.d0 + pstar(2)*rhostar(3)/rhostar(2) - pstar(3)) +
     &      (pstar(4)*rhostar(3) + pstar(2)*
     &      (rhostar(5)-rhostar(3)*rhostar(4)/rhostar(2)))/
     &      rhostar(2) - pstar(5))/t
          xfixedtt = x*((2.d0*rhostar(3) - pstar(3))*
     &      (1.d0 + pstar(2)*rhostar(3)/rhostar(2) - pstar(3)) +
     &      (pstar(5)*rhostar(3) + pstar(2)*
     &      (rhostar(6)-rhostar(3)*rhostar(5)/rhostar(2)))/
     &      rhostar(2) - pstar(6))/t
        else
          xfixedt = 0.d0
          xfixedtf = 0.d0
          xfixedtt = 0.d0
        endif
        dscoulomb = -fcoulombt -
     &    fcoulomb0*sum0atl/t - fcoulomb2*sum2atl/t - 
     &    fcoulombx*xfixedt
C        partial fcoulombt(f,t) wrt ln f
        fcfttf = 
     &    fcoulomb0t*sum0f + fcoulomb2t*sum2f +
     &    n_e*rhostar(2)*(fcoulomb0t*sum0ane +
     &    fcoulomb2t*sum2ane) + fcoulombtx*xf
C        partial fcoulombt(f,t) wrt ln t
        fcfttt = 
     &    fcoulomb0t*sum0t + fcoulomb2t*sum2t +
     &    n_e*rhostar(3)*(fcoulomb0t*sum0ane +
     &    fcoulomb2t*sum2ane) + fcoulomb0t*sum0atl +
     &    fcoulomb2t*sum2atl + fcoulombtt*t + fcoulombtx*xt
        if(ifcoulomb.gt.4) then
C          DeWitt definition
C          sum0a = sum0 + theta_e*n_e
C          partial sum0atl(f,t) wrt ln f
          sum0atlf = sum0atl*rhostar(2) + n_e*thetatf
C          partial sum0atl(f,t) wrt ln t
          sum0atlt = sum0atl*rhostar(3) + n_e*thetatt
        else
C          PTEH definition
C          sum0a = sum0
          sum0atlf = 0.d0
          sum0atlt = 0.d0
        endif
C        partial sum2atl(f,t) wrt ln f
        sum2atlf = sum2atl*rhostar(2) + n_e*thetatf
C        partial sum2atl(f,t) wrt ln t
        sum2atlt = sum2atl*rhostar(3) + n_e*thetatt
        dscoulombf = -fcfttf -
     &    fcft0f*sum0atl/t - fcoulomb0*sum0atlf/t -
     &    fcft2f*sum2atl/t - fcoulomb2*sum2atlf/t -
     &    fcftxf*xfixedt - fcoulombx*xfixedtf
        dscoulombt = -fcfttt -
     &    fcft0t*sum0atl/t - fcft2t*sum2atl/t -
     &    fcoulomb0*(sum0atlt/t - sum0atl/t) -
     &    fcoulomb2*(sum2atlt/t - sum2atl/t) -
     &    fcftxt*xfixedt - fcoulombx*xfixedtt
C         from u = f + Ts
        ducoulomb = fcoulomb + t*dscoulomb
C        dpcoulomb = - partial(V*fcoulomb(n_e, n_i, T)) wrt V
C          = -fcoulomb + n_e*partial fcoulomb wrt n_e + sum_i n_i*
C          partial fcoulomb wrt n_i
C          = -fcoulomb + n_e*dmucoulomb + sum0*partial fcoulomb
C          wrt sum0a + sum2*partial fcoulomb wrt sum2a
C        n.b. the latter two terms are dropped if sum0 and sum2 
C          depend only on n_e.  This occurs if
C          the PTEH approximation being used for sum0, sum2 
        if(if_pteh.eq.1) then
          dpcoulomb = -fcoulomb + n_e*dmucoulomb
C          dpcoulomb(f,t)/ ln f and ln t:  to do these derivatives
C          fcoulomb (for the fully ionized case) is considered 
C          to be a function of (sum0(n_e), sum2(n_e)+
C          theta_e(n_e,t)*n_e, t, x(n_e,T)) = function(n_e, T),
C          completely independent of n_i.
C          also recall dmucoulomb = partial fcoulomb(n_e, T) wrt n_e
C          in the fully ionized case.
C          There is a useful cancellation of terms when using 
C          this viewpoint....
          dpcoulombf = n_e*dmucoulombf
C          dscoulomb = partial -fcoulomb(n_e, t) wrt t
          dpcoulombt = t*dscoulomb + n_e*dmucoulombt
C          temporary check confirming definition of dmucoulomb
!          dpcoulomb = fcoulomb
!          dpcoulombf = dmucoulomb*n_e*rhostar(2)
!          dpcoulombt = dmucoulomb*n_e*rhostar(3) - t*dscoulomb
        else
          dpcoulomb = -fcoulomb + n_e*dmucoulomb +
     &      sum0*fcoulomb0 + sum2*fcoulomb2
C          dpcoulomb(f,t)/ ln f and ln t:  to do these derivatives
C          fcoulomb is considered a function of (sum0(f,t), sum2(f,t)+
C          theta_e(n_e,t)*n_e, t, x(n_e,T)).  In other words, the previous
C          n_i dependence is replaced with sum0 and sum2 dependence.
C          There is a useful cancellation of terms when using 
C          this viewpoint....
          dpcoulombf = n_e*dmucoulombf +
     &      sum0*fcft0f + sum2*fcft2f
C          dscoulomb = partial -fcoulomb(sum0, sum2, n_e, t) wrt t
          dpcoulombt = t*dscoulomb + n_e*dmucoulombt +
     &      sum0*fcft0t + sum2*fcft2t
        endif
      else
        dpcoulomb = 0.d0
        dpcoulombf = 0.d0
        dpcoulombt = 0.d0
        dscoulomb = 0.d0
        dscoulombf = 0.d0
        dscoulombt = 0.d0
        ducoulomb = 0.d0
      endif
      end
