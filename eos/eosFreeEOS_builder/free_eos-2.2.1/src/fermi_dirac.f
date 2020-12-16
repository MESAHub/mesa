C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fermi_dirac.f 370 2006-11-29 23:57:39Z airwin $
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
      subroutine fermi_dirac(fl, tcl, rhostar, pstar, sstar, ustar,
     &  mstar, morder_in)
C       subroutine to calculate *scaled* n_e, p_e, s_e, and u_e and the
C       derivative of the *ln* of these quantities wrt fl, and tcl.  value
C       and the fl and tcl derivatives are stored
C       in the arrays rhostar(9), pstar(9), sstar(3), and ustar(3).
C       The six extra derivatives of ln rhostar and pstar are in 
C       the order flfl, fltcl, tcltcl, flflfl, flfltcl, fltcltcl.
C       fl is ln f and tcl is ln tc, where tc = kt/mc^2 equiv beta. (see Eggleton, 
C       Faulkner and Flannery paper and vdb notes).
C       subroutine written by Alan W. Irwin, October, 1989.
C       scaling:
C       n_e (number density of electrons) = 8 pi/lambda_c^3 rhostar
C        ==> rhostar = sqrt(2) beta^3/2 [F(1/2 + beta F(3/2)] (C.G. 24.98)
C       p_e (electron pressure) = 8 pi/lambda_c^3 *m c^2 pstar
C        ==> pstar = (2/3)sqrt(2) beta^5/2 [F(3/2 + (1/2)beta F(5/2)] (C.G. 24.99)
C       s_e (entropy per unit volume of the free electrons) CG. 24.76b
C         = (u_e + p_e)/T - psi n_e k
C         = n_e k sstar
C       note u_e (see later) is the internal energy per unit volume of the 
C         electrons.
C       internally we calculate:
C       qstar, where the convention for qstar in this programme is
C       (1+f)/g times the qstar defined in the paper.
C       Therefore, the following transformations apply using this convention:
C       sstar = qstar/(rhostar*sqrt(1+f)) + 2 sqrt(1+f) - psi
C         = qstar/(rhostar*sqrt(1+f)) + ln((1+sqrt(1+f))^2/f)
C       finally we calculate
C       ustar = sstar + psi - pstar sqrt(1+f)/(g rhostar)
C       this is related to u_e via
C       u_e = n_e k T  ustar = n_e m c^2 beta ustar
C        ==> ustar = [F(3/2 + beta F(5/2)]/[F(1/2 + beta F(3/2)] (C.G. 24.100)
C         = s_e T + psi n_e k T - pstar mc^2 n_e/rhostar
C         = s_e T + psi n_e k T - p_e (checks with previous).
C       n.b. the ustar in this routine is different from the ustar in
C       the vdb notes.
C       more relation to CG:  
C       sstar = (u_e + p_e)/n_e k T - psi
C         = ustar + pstar/(rhostar beta) - psi
C         = ustar + (2/3) [F(3/2 + (1/2)beta F(5/2)]/[F(1/2 + beta F(3/2)] - psi
C         = [(5/3)F(3/2 + (4/3)beta F(5/2)]/[F(1/2 + beta F(3/2)] - psi
C       morder_in = 3, 5, or 8  (or 13, 15, or 18) means 3rd, 5th, or 8th order
C       fermi_dirac integral approximation following EFF fit (or modified
C       version of EFF fit which reduces to Cody-Thacher approximation for
C       low relativistic correction).
C       morder_in = -3, -5, or -8 (or -13, -15, or -18) means_use_above
C       approximations in non-relativistic limit.
C       morder_in = -1 or 1 uses Cody-Thacher approximation directly.
C       morder_in = 21 calculate Fermi-Dirac integrals with slow, but precise
C       (~1.d-9 relative errors) numerical integration.
C       morder_in = -21 is same as 21 in non-relativistic limit.
C       morder_in = 23 means original 3rd order eff result
C       morder_in = -23 is same as 23 in non-relativistic limit.
      implicit none
      integer mstar, morder_in, morder, norder
      double precision rhostar(mstar+6),pstar(mstar+6),sstar(mstar), ustar(mstar),
     &  fl, tcl, rhostarcon, pstarcon
      integer mderiv, mordermax
      parameter (mderiv = 2,mordermax = 8)
      integer maxfd_direct
      parameter (maxfd_direct = 9)
C       n.b ccoeff in order of rhocoeff, pcoeff, qcoeff
C       n.b. dimensions on ccoeff and plowcoeff large enough for
C         morder = 8 call to fermi_dirac_coeff.
      double precision ccoeff(3*(mordermax+1)*(mordermax+1)),
     &  sum(mderiv+1+7),
     &  fd(0:4), fermi_dirac_ct,
     &  fd_direct(maxfd_direct,3)
      integer morder_old, i, ideriv, nderiv, ioffset, moffset
      double precision f, wf, tc, g, vf, vg, uf, duf, duf2, ug, dug, dug2,
     &  fdf, psi, qconst, pconst, sumadd, power, powerufm1, poweruf
C       must be invalid value!
      data morder_old/0/
      save
      f = exp(fl)
C       d psi d fl
      wf = sqrt(1.d0+f)
      if(f.gt.1.d0) then
        psi = 2.d0*wf + log((wf-1.d0)/(wf+1.d0))
      else
        psi = fl + 2.d0*(wf - log(wf+1.d0))
      endif
C      a.k.a. beta = kt/mc^2
      tc = exp(tcl)
      morder = abs(morder_in)
      if(morder.eq.1) then
C         this branch taken independent of all else to_use_Cody-Thacher
C         approximations to Fermi-Dirac integral in *non* relativistic
C         limit.
        do ideriv = 0, 4
C           F(3/2) and 4 derivatives
          fd(ideriv) = fermi_dirac_ct(psi,ideriv,0)
        enddo
C        N.B. F 1/2 = (2/3) d F 3/2 d psi so that is where the factor
C        of (2/3) comes from.
        rhostarcon = sqrt(2.d0) * exp(1.5d0*tcl) * (2.d0/3.d0)
        rhostar(1) = rhostarcon*fd(1)
C        ln derivatives in order function, f, t, ff, ft, tt, fff, fft, ftt
        rhostar(2) = fd(2)/fd(1)*wf
        rhostar(3) = 1.5d0
        rhostar(4) = (-fd(2)*fd(2)*(1.d0+f)/fd(1) +
     &    fd(3)*(1.d0+f) + 0.5d0*fd(2)*f/wf)/fd(1)
        rhostar(5) = 0.d0
        rhostar(6) = 0.d0
        rhostar(7) = -rhostar(4)*(fd(2)/fd(1))*wf + (
     &    -fd(2)*fd(2)*f/fd(1) -
     &    2.d0*wf*fd(2)*fd(3)*(1.d0+f)/fd(1) +
     &    wf*fd(2)*(fd(2)/fd(1))*(fd(2)/fd(1))*(1.d0+f) +
     &    fd(4)*(1.d0+f)*wf + 1.5d0*fd(3)*f +
     &    0.5d0*fd(2)*(f/wf)*(1.d0 - 0.5d0*f/(1.d0+f)))/fd(1)
        rhostar(8) = 0.d0
        rhostar(9) = 0.d0
        pstarcon = sqrt(2.d0) * exp(2.5d0*tcl) * (2.d0/3.d0)
        pstar(1) = pstarcon*fd(0)
        pstar(2) = fd(1)/fd(0)*wf
        pstar(3) = 2.5d0
        pstar(4) = (-fd(1)*fd(1)*(1.d0+f)/fd(0) +
     &    fd(2)*(1.d0+f) + 0.5d0*fd(1)*f/wf)/fd(0)
        pstar(5) = 0.d0
        pstar(6) = 0.d0
        pstar(7) = -pstar(4)*(fd(1)/fd(0))*wf + (
     &    -fd(1)*fd(1)*f/fd(0) -
     &    2.d0*wf*fd(1)*fd(2)*(1.d0+f)/fd(0) +
     &    wf*fd(1)*(fd(1)/fd(0))*(fd(1)/fd(0))*(1.d0+f) +
     &    fd(3)*(1.d0+f)*wf + 1.5d0*fd(2)*f +
     &    0.5d0*fd(1)*(f/wf)*(1.d0 - 0.5d0*f/(1.d0+f)))/fd(0)
        pstar(8) = 0.d0
        pstar(9) = 0.d0
        ustar(1) = 1.5d0*fd(0)/fd(1)
        ustar(2) = wf*(fd(1)/fd(0) - fd(2)/fd(1))
        ustar(3) = 0.d0
C        significance loss for large degeneracy to be fixed later
C        if ever. (only loses 3 figures at psi = 100)
        sstar(1) = 2.5d0*fd(0)/fd(1) - psi
        sstar(2) = wf*(1.5d0 - 2.5d0*(fd(0)/fd(1))*(fd(2)/fd(1)))/
     &    sstar(1)
        sstar(3) = 0.d0
        return
      elseif(abs(morder_in).eq.21) then
C        branch for direct, near-exact, but slow calculation of Fermi-Dirac
C        integrals using numerical integration.
        if(morder_in.eq.-21) then
          call fermi_dirac_direct(0.5d0, psi, 0.d0, maxfd_direct-1,
     &      fd_direct(1,1))
          call fermi_dirac_direct(1.5d0, psi, 0.d0, maxfd_direct-1,
     &      fd_direct(1,2))
C          derivatives in order function, f, t, ff, ft, tt, fff, fft, ftt
          rhostar(1) = fd_direct(1,1)
          rhostar(2) = fd_direct(2,1)
          rhostar(4) = fd_direct(4,1)
          rhostar(7) = fd_direct(7,1)
C          Zero beta derivatives.
          rhostar(3) = 0.d0
          rhostar(5) = 0.d0
          rhostar(6) = 0.d0
          rhostar(8) = 0.d0
          rhostar(9) = 0.d0
          rhostar(10) = 0.d0
          pstar(1) = fd_direct(1,2)
          pstar(2) = fd_direct(2,2)
          pstar(4) = fd_direct(4,2)
          pstar(7) = fd_direct(7,2)
C          Zero beta derivatives.
          pstar(3) = 0.d0
          pstar(5) = 0.d0
          pstar(6) = 0.d0
          pstar(8) = 0.d0
          pstar(9) = 0.d0
          pstar(10) = 0.d0
          ustar(1) = fd_direct(1,2)
          ustar(2) = fd_direct(2,2)
C          Zero beta derivatives.
          ustar(3) = 0.d0
        else
          call fermi_dirac_direct(0.5d0, psi, tc, maxfd_direct-1,
     &      fd_direct(1,1))
          call fermi_dirac_direct(1.5d0, psi, tc, maxfd_direct-1,
     &      fd_direct(1,2))
          call fermi_dirac_direct(2.5d0, psi, tc, maxfd_direct-1,
     &      fd_direct(1,3))
C          derivatives in order function, f, t, ff, ft, tt, fff, fft, ftt
          do i = 1, maxfd_direct
            rhostar(i) = fd_direct(i,1) + tc*fd_direct(i,2)
            pstar(i) = fd_direct(i,2) + 0.5d0*tc*fd_direct(i,3)
            if(i.le.3) then
              ustar(i) = fd_direct(i,2) + tc*fd_direct(i,3)
            endif
          enddo
C          Transform to total beta=tc derivative of above rhostar, pstar.
          rhostar(3) = rhostar(3) + fd_direct(1,2)
          rhostar(5) = rhostar(5) + fd_direct(2,2)
          rhostar(6) = rhostar(6) + 2.d0*fd_direct(3,2)
          rhostar(8) = rhostar(8) + fd_direct(4,2)
          rhostar(9) = rhostar(9) + 2.d0*fd_direct(5,2)
          pstar(3) = pstar(3) + 0.5d0*fd_direct(1,3)
          pstar(5) = pstar(5) + 0.5d0*fd_direct(2,3)
          pstar(6) = pstar(6) + fd_direct(3,3)
          pstar(8) = pstar(8) + 0.5d0*fd_direct(4,3)
          pstar(9) = pstar(9) + fd_direct(5,3)
          ustar(3) = ustar(3) + fd_direct(1,3)
C          Transform to tcl = ln beta derivatives from tc=beta derivatives
          rhostar(9) = rhostar(9)*tc*tc + rhostar(5)*tc
          rhostar(8) = rhostar(8)*tc
          rhostar(6) = rhostar(6)*tc*tc + rhostar(3)*tc
          rhostar(5) = rhostar(5)*tc
          rhostar(3) = rhostar(3)*tc
          pstar(9) = pstar(9)*tc*tc + pstar(5)*tc
          pstar(8) = pstar(8)*tc
          pstar(6) = pstar(6)*tc*tc + pstar(3)*tc
          pstar(5) = pstar(5)*tc
          pstar(3) = pstar(3)*tc
          ustar(3) = ustar(3)*tc
        endif
C        Transform to fl derivatives from psi derivatives.  wf = dpsi/dfl
        rhostar(9) = rhostar(9)*wf
        rhostar(8) = rhostar(8)*wf*wf + rhostar(5)*(0.5d0*f/wf)
        rhostar(7) = rhostar(7)*wf*wf*wf + rhostar(4)*(1.5d0*f) +
     &    rhostar(2)*0.5d0*f/wf*(1.d0 - 0.5d0*(f/wf)/wf)
        rhostar(5) = rhostar(5)*wf
        rhostar(4) = rhostar(4)*wf*wf + rhostar(2)*(0.5d0*f/wf)
        rhostar(2) = rhostar(2)*wf
        pstar(9) = pstar(9)*wf
        pstar(8) = pstar(8)*wf*wf + pstar(5)*(0.5d0*f/wf)
        pstar(7) = pstar(7)*wf*wf*wf + pstar(4)*(1.5d0*f) +
     &    pstar(2)*0.5d0*f/wf*(1.d0 - 0.5d0*(f/wf)/wf)
        pstar(5) = pstar(5)*wf
        pstar(4) = pstar(4)*wf*wf + pstar(2)*(0.5d0*f/wf)
        pstar(2) = pstar(2)*wf
        ustar(2) = ustar(2)*wf
C        Transform to derivatives of ln of function, but do not transform
C        function itself (which is the convention used in this routine).
        rhostar(9) = rhostar(9)/rhostar(1) -
     &    (rhostar(6)/rhostar(1))*(rhostar(2)/rhostar(1)) -
     &    2.d0*(rhostar(5)/rhostar(1))*(rhostar(3)/rhostar(1)) +
     &    2.d0*(rhostar(3)/rhostar(1))*(rhostar(3)/rhostar(1))*
     &    (rhostar(2)/rhostar(1))
        rhostar(8) = rhostar(8)/rhostar(1) -
     &    (rhostar(4)/rhostar(1))*(rhostar(3)/rhostar(1)) -
     &    2.d0*(rhostar(5)/rhostar(1))*(rhostar(2)/rhostar(1)) +
     &    2.d0*(rhostar(3)/rhostar(1))*(rhostar(2)/rhostar(1))*
     &    (rhostar(2)/rhostar(1))
        rhostar(7) = rhostar(7)/rhostar(1) -
     &    3.d0*(rhostar(4)/rhostar(1))*(rhostar(2)/rhostar(1)) +
     &    2.d0*(rhostar(2)/rhostar(1))*(rhostar(2)/rhostar(1))*
     &    (rhostar(2)/rhostar(1))
        rhostar(6) = rhostar(6)/rhostar(1) -
     &    (rhostar(3)/rhostar(1))*(rhostar(3)/rhostar(1))
        rhostar(5) = rhostar(5)/rhostar(1) -
     &    (rhostar(2)/rhostar(1))*(rhostar(3)/rhostar(1))
        rhostar(4) = rhostar(4)/rhostar(1) -
     &    (rhostar(2)/rhostar(1))*(rhostar(2)/rhostar(1))
        rhostar(3) = rhostar(3)/rhostar(1)
        rhostar(2) = rhostar(2)/rhostar(1)
        pstar(9) = pstar(9)/pstar(1) -
     &    (pstar(6)/pstar(1))*(pstar(2)/pstar(1)) -
     &    2.d0*(pstar(5)/pstar(1))*(pstar(3)/pstar(1)) +
     &    2.d0*(pstar(3)/pstar(1))*(pstar(3)/pstar(1))*
     &    (pstar(2)/pstar(1))
        pstar(8) = pstar(8)/pstar(1) -
     &    (pstar(4)/pstar(1))*(pstar(3)/pstar(1)) -
     &    2.d0*(pstar(5)/pstar(1))*(pstar(2)/pstar(1)) +
     &    2.d0*(pstar(3)/pstar(1))*(pstar(2)/pstar(1))*
     &    (pstar(2)/pstar(1))
        pstar(7) = pstar(7)/pstar(1) -
     &    3.d0*(pstar(4)/pstar(1))*(pstar(2)/pstar(1)) +
     &    2.d0*(pstar(2)/pstar(1))*(pstar(2)/pstar(1))*
     &    (pstar(2)/pstar(1))
        pstar(6) = pstar(6)/pstar(1) -
     &    (pstar(3)/pstar(1))*(pstar(3)/pstar(1))
        pstar(5) = pstar(5)/pstar(1) -
     &    (pstar(2)/pstar(1))*(pstar(3)/pstar(1))
        pstar(4) = pstar(4)/pstar(1) -
     &    (pstar(2)/pstar(1))*(pstar(2)/pstar(1))
        pstar(3) = pstar(3)/pstar(1)
        pstar(2) = pstar(2)/pstar(1)
        ustar(3) = ustar(3)/ustar(1)
        ustar(2) = ustar(2)/ustar(1)
C        N.B. at this stage we have the following:
C        rhostar = F 1/2 + beta F 3/2
C        pstar = F 3/2 + 0.5 beta F 5/2
C        ustar = F 3/2 + beta F 5/2
C         plus partial derivatives of the ln of these quantities.
C
C        Transform to final forms:
C        ustar to ratio of present ustar and rhostar.
        ustar(1) = ustar(1)/rhostar(1)
        ustar(2) = ustar(2) - rhostar(2)
        ustar(3) = ustar(3) - rhostar(3)
C        Multiply by factor to convert to rhostar (or add 1.5 tcl)
C        to ln rhostar for partial-derivative purposes).
C        N.B. no factor of 2/3 in this version since we are dealing with
C        FD integrals rather than their derivatives as above for the
C        non-releativistic Cody-Thacher case.
        rhostarcon = sqrt(2.d0) * exp(1.5d0*tcl)
        rhostar(1) = rhostarcon*rhostar(1)
        rhostar(3) = rhostar(3) + 1.5d0
C        Multiply by factor to convert to pstar (or add 2.5 tcl)
C        to ln pstar for partial-derivative purposes).
        pstarcon = sqrt(2.d0) * exp(2.5d0*tcl) * (2.d0/3.d0)
        pstar(1) = pstarcon*pstar(1)
        pstar(3) = pstar(3) + 2.5d0
C       finally we calculate
C       sstar = ustar - psi + pstar sqrt(1+f)/(g rhostar)
        sstar(1) = ustar(1) - psi + pstar(1)/(tc*rhostar(1))
        sstar(2) = (ustar(1)*ustar(2) - wf +
     &    (pstar(1)/(tc*rhostar(1)))*(pstar(2)-rhostar(2)))/
     &    sstar(1)
        sstar(3) = (ustar(1)*ustar(3) +
     &    (pstar(1)/(tc*rhostar(1)))*(pstar(3)-rhostar(3)-1.d0))/
     &    sstar(1)
        return
      endif
C      at this stage have dealt with mord_in = -1, 1, -21, or 21 result.
C      now do rest of possibilities (which include special case of
C      morder_in = -23 or 23).
C      N.B. morder = abs(morder_in) (see above).
      if(mstar.ne.mderiv+1) stop 'invalid mstar input to fermi_dirac'
      if(abs(morder_in).gt.20) then
        morder = morder - 20
      elseif(abs(morder_in).gt.10) then
        morder = morder - 10
      endif
      if(abs(morder_in).ne.morder_old) then
        morder_old = abs(morder_in)
C         fill coefficients with thermodynamically consistent set.
        if(abs(morder_in).gt.20) then
          call fermi_dirac_original_coeff(ccoeff,morder)
        elseif(abs(morder_in).gt.10) then
          call fermi_dirac_minusnr_coeff(ccoeff,morder)
C          Precalculate factors to convert from fd to non-relativistic
C          pstar and rhostar.
          pstarcon = sqrt(2.d0) * (2.d0/3.d0)
        else
          call fermi_dirac_coeff(ccoeff,morder)
        endif
      endif
      if((.not.abs(morder_in).gt.20).and.abs(morder_in).gt.10) then
C        Collect zero-temperature Cody-Thacher results that were subtracted
C        in the fit used to generate the fermi_dirac_minusnr_coeff results.
        do ideriv = 0, 4
C          F(3/2) and 4 derivatives
          fd(ideriv) = fermi_dirac_ct(psi,ideriv,0)
        enddo
      endif
C       f = exp(fl)
C       wf = sqrt(1.d0+f)
      vf = 1.d0/(1.d0+f)
C      partial of ln (1+f) wrt ln f
      uf = f*vf
C      partial uf wrt ln f.  N.B. 1-uf = 1/(1+f)
      duf = uf*vf
      duf2 = uf*(1.d0-f)*vf*vf
      if(morder_in.gt.0) then
        g = tc*wf
        vg = 1.d0/(1.d0+g)
        fdf = g*g+g
        fdf = uf*fdf*sqrt(fdf)*(vf*vg)**morder
C      partial of ln (1+g) wrt ln g
        ug = g*vg
C        partial ug wrt ln g
        dug = ug*vg
        dug2 = ug*(1.d0-g)*vg*vg
      else
C         negative morder_in flags NR limit.
        g = 0.d0
        vg = 1.d0
C         d ln(1+g)/dln g in limit where 1+g replaced by 1.
        ug = 0.d0
        fdf = uf*(tc*wf)*sqrt(tc*wf)*vf**morder
        dug = 0.d0
        dug2 = 0.d0
      endif
      if((.not.abs(morder_in).gt.20).and.abs(morder_in).gt.10) then
C        When non-relativistic limit is treated separately, lowest order
C        (in g) coefficients are zero.  Ignore these coefficients in
C        sum and (later) multiply series results by an additional factor
C        of g.
        norder = morder - 1
        moffset = 1 + (morder+1)
      else
        norder = morder
        moffset = 1
      endif
      do i = 1,3
        if(i.le.2) then
          nderiv = 9
        else
          nderiv = 2
        endif
        ioffset = (i-1)*(morder+1)*(morder+1) + moffset
C        Note for 10< abs(morder_in)<= 20, the zero-order g coefficients are
C        zero so that norder and moffset are manipulated so that
C        effsum_calc just calculates polynomial multiplier for g.
        call effsum_calc(f, g, morder, norder, ccoeff(ioffset),
     &    nderiv, sum)
        if((.not.abs(morder_in).gt.20).and.abs(morder_in).gt.10) then
C          Multiply sum by g and adjust sum derivatives accordingly since
C          for this case effsum_calc (with the appropriate norder and
C          moffset) just calculates polynomial multiplier of g.
          if(i.le.2) then
C            ggg
            sum(10) = g*(sum(1) + 3.d0*sum(3) + 3.d0*sum(6) + sum(10))
C            fgg
            sum(9) = g*(sum(2) + 2.d0*sum(5) + sum(9))
C            ffg
            sum(8) = g*(sum(4) + sum(8))
C            fff
            sum(7) = g*sum(7)
C            gg
            sum(6) = g*(sum(1) + 2.d0*sum(3) + sum(6))
C            fg
            sum(5) = g*(sum(2) + sum(5))
C            ff
            sum(4) = g*sum(4)
          endif
C          g
          sum(3) = g*(sum(1) + sum(3))
C          f
          sum(2) = g*sum(2)
          sum(1) = g*sum(1)
C          Add in zero-temperature Cody-Thacher results that were subtracted
C          in the fit used to generate the fermi_dirac_minusnr_coeff results.
          if(i.eq.1) then
            power = dble(morder)+0.25d0
            powerufm1 = power*uf - 1.d0
            sumadd = pstarcon*(1.d0+g)*fd(1)*((1.d0+f)**power)/f
            sum(1) = sum(1) + sumadd
            sum(2) = sum(2) + sumadd*(fd(2)*wf/fd(1) + powerufm1)
            sum(3) = sum(3) + ug*sumadd
            sum(4) = sum(4) + sumadd*(
     &        (fd(3)*wf + 0.5d0*fd(2)*uf)*wf/fd(1) + power*duf +
     &        powerufm1*(2.d0*fd(2)*wf/fd(1) + powerufm1))
            sum(5) = sum(5) + ug*sumadd*(fd(2)*wf/fd(1) + powerufm1)
            sum(6) = sum(6) + ug*sumadd
C            sum(7) = sum(7) + sumadd*(
C     &        ((fd(4)*wf + fd(3)*uf)*wf + 0.5d0*fd(2)*duf)*wf/fd(1) +
C     &        (fd(3)*wf + 0.5d0*fd(2)*uf)*wf*
C     &        (0.5d0*uf - fd(2)*wf/fd(1))/fd(1) + power*duf2 +
C     &        power*duf*(2.d0*fd(2)*wf/fd(1) + 2.d0*powerufm1) +
C     &        powerufm1*2.d0*(
C     &        fd(3)*wf + 0.5d0*fd(2)*uf -
C     &        fd(2)*fd(2)*wf/fd(1))*wf/fd(1) + (
C     &        (fd(3)*wf + 0.5d0*fd(2)*uf)*wf/fd(1) + power*duf +
C     &        powerufm1*(2.d0*fd(2)*wf/fd(1) + powerufm1))*(
C     &        fd(2)*wf/fd(1) + powerufm1))
C            commented out expression above is straight derivative which has
C            been tested to be correct.  Expression below has terms
C            consolidated from above.  This is also tested to be correct,
C            but it a lot less understandable than expression above.
            sum(7) = sum(7) + sumadd*(
     &        (
     &        (fd(4)*wf + 1.5d0*fd(3)*uf)*wf + 0.5d0*fd(2)*duf
     &        )*wf/fd(1) +
     &        (0.5d0*fd(2)*uf)*wf*(0.5d0*uf)/fd(1) +
     &        power*duf2 + 
     &        power*duf*(3.d0*fd(2)*wf/fd(1)) +
     &        powerufm1*(
     &        3.d0*(fd(3)*wf + 0.5d0*fd(2)*uf)*wf/fd(1) +
     &        3.d0*power*duf + powerufm1*(
     &        3.d0*fd(2)*wf/fd(1) + powerufm1))
     &        )
            sum(8) = sum(8) + ug*sumadd*(
     &        (fd(3)*wf + 0.5d0*fd(2)*uf)*wf/fd(1) +
     &        power*duf +
     &        powerufm1*(2.d0*fd(2)*wf/fd(1) +
     &        powerufm1))
            sum(9) = sum(9) + ug*sumadd*(fd(2)*wf/fd(1) + powerufm1)
            sum(10) = sum(10) + ug*sumadd
C            add in second term
            power = dble(morder)-1.25d0
            poweruf = power*uf
            sumadd = 0.5d0*(2.5d0 - dble(morder))*pstarcon*
     &        g*fd(0)*((1.d0+f)**power)
            sum(1) = sum(1) + sumadd
            sum(2) = sum(2) + sumadd*(fd(1)*wf/fd(0) + poweruf)
            sum(3) = sum(3) + sumadd
            sum(4) = sum(4) + sumadd*(
     &        (fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) + power*duf +
     &        poweruf*(2.d0*fd(1)*wf/fd(0) + poweruf))
            sum(5) = sum(5) + sumadd*(fd(1)*wf/fd(0) + poweruf)
            sum(6) = sum(6) + sumadd
C            sum(7) = sum(7) + sumadd*(
C     &        ((fd(3)*wf + fd(2)*uf)*wf + 0.5d0*fd(1)*duf)*wf/fd(0) +
C     &        (fd(2)*wf + 0.5d0*fd(1)*uf)*wf*
C     &        (0.5d0*uf - fd(1)*wf/fd(0))/fd(0) + power*duf2 +
C     &        power*duf*(2.d0*fd(1)*wf/fd(0) + 2.d0*poweruf) +
C     &        poweruf*2.d0*(
C     &        fd(2)*wf + 0.5d0*fd(1)*uf -
C     &        fd(1)*fd(1)*wf/fd(0))*wf/fd(0) + (
C     &        (fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) + power*duf +
C     &        poweruf*(2.d0*fd(1)*wf/fd(0) + poweruf))*(
C     &        fd(1)*wf/fd(0) + poweruf))
C            commented out expression above is straight derivative which has
C            been tested to be correct.  Expression below has terms
C            consolidated from above.  This is also tested to be correct,
C            but it a lot less understandable than expression above.
            sum(7) = sum(7) + sumadd*(
     &        (
     &        (fd(3)*wf + 1.5d0*fd(2)*uf)*wf + 0.5d0*fd(1)*duf
     &        )*wf/fd(0) +
     &        (0.5d0*fd(1)*uf)*wf*(0.5d0*uf)/fd(0) +
     &        power*duf2 + 
     &        power*duf*(3.d0*fd(1)*wf/fd(0)) +
     &        poweruf*(
     &        3.d0*(fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) +
     &        3.d0*power*duf + poweruf*(
     &        3.d0*fd(1)*wf/fd(0) + poweruf))
     &        )
            sum(8) = sum(8) + sumadd*(
     &        (fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) + power*duf +
     &        poweruf*(2.d0*fd(1)*wf/fd(0) + poweruf))
            sum(9) = sum(9) + sumadd*(fd(1)*wf/fd(0) + poweruf)
            sum(10) = sum(10) + sumadd
          elseif(i.eq.2) then
            power = dble(morder)-0.25d0
            powerufm1 = power*uf - 1.d0
            sumadd = pstarcon*(1.d0+g)*fd(0)*((1.d0+f)**power)/f
            sum(1) = sum(1) + sumadd
            sum(2) = sum(2) + sumadd*(fd(1)*wf/fd(0) + powerufm1)
            sum(3) = sum(3) + ug*sumadd
            sum(4) = sum(4) + sumadd*(
     &        (fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) + power*duf +
     &        powerufm1*(2.d0*fd(1)*wf/fd(0) + powerufm1))
            sum(5) = sum(5) + ug*sumadd*(fd(1)*wf/fd(0) + powerufm1)
            sum(6) = sum(6) + ug*sumadd
C            sum(7) = sum(7) + sumadd*(
C     &        ((fd(3)*wf + fd(2)*uf)*wf + 0.5d0*fd(1)*duf)*wf/fd(0) +
C     &        (fd(2)*wf + 0.5d0*fd(1)*uf)*wf*
C     &        (0.5d0*uf - fd(1)*wf/fd(0))/fd(0) + power*duf2 +
C     &        power*duf*(2.d0*fd(1)*wf/fd(0) + 2.d0*powerufm1) +
C     &        powerufm1*2.d0*(
C     &        fd(2)*wf + 0.5d0*fd(1)*uf -
C     &        fd(1)*fd(1)*wf/fd(0))*wf/fd(0) + (
C     &        (fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) + power*duf +
C     &        powerufm1*(2.d0*fd(1)*wf/fd(0) + powerufm1))*(
C     &        fd(1)*wf/fd(0) + powerufm1))
C            commented out expression above is straight derivative which has
C            been tested to be correct.  Expression below has terms
C            consolidated from above.  This is also tested to be correct,
C            but it a lot less understandable than expression above.
            sum(7) = sum(7) + sumadd*(
     &        (
     &        (fd(3)*wf + 1.5d0*fd(2)*uf)*wf + 0.5d0*fd(1)*duf
     &        )*wf/fd(0) +
     &        (0.5d0*fd(1)*uf)*wf*(0.5d0*uf)/fd(0) +
     &        power*duf2 + 
     &        power*duf*(3.d0*fd(1)*wf/fd(0)) +
     &        powerufm1*(
     &        3.d0*(fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) +
     &        3.d0*power*duf + powerufm1*(
     &        3.d0*fd(1)*wf/fd(0) + powerufm1))
     &        )
            sum(8) = sum(8) + ug*sumadd*(
     &        (fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) +
     &        power*duf +
     &        powerufm1*(2.d0*fd(1)*wf/fd(0) +
     &        powerufm1))
            sum(9) = sum(9) + ug*sumadd*(fd(1)*wf/fd(0) + powerufm1)
            sum(10) = sum(10) + ug*sumadd
          elseif(i.eq.3) then
C            add in first term.
            power = dble(morder)+0.75d0
            powerufm1 = power*uf - 1.d0
            sumadd = 2.5d0*pstarcon*
     &        (1.d0+g)*fd(0)*((1.d0+f)**power)/f
            sum(1) = sum(1) + sumadd
            sum(2) = sum(2) + sumadd*(fd(1)*wf/fd(0) + powerufm1)
            sum(3) = sum(3) + ug*sumadd
C            add in second term.
            power = dble(morder)+1.25d0
            powerufm1 = power*uf - 1.d0
            sumadd = -2.d0*pstarcon*
     &        (1.d0+g)*fd(1)*((1.d0+f)**power)/f
            sum(1) = sum(1) + sumadd
            sum(2) = sum(2) + sumadd*(fd(2)*wf/fd(1) + powerufm1)
            sum(3) = sum(3) + ug*sumadd
C            add in third term.
            power = dble(morder)-0.25d0
            powerufm1 = power*uf - 1.d0
            sumadd = (2.5d0 - dble(morder))*pstarcon*
     &        g*fd(0)*((1.d0+f)**power)/f
            sum(1) = sum(1) + sumadd
            sum(2) = sum(2) + sumadd*(fd(1)*wf/fd(0) + powerufm1)
            sum(3) = sum(3) + sumadd
          endif
        endif
C        Transform to derivatives of ln of function, but do not transform
C        function itself (which is the convention used in this routine).
        if(i.le.2) then
          sum(10) = sum(10)/sum(1) -
     &      3.d0*(sum(6)/sum(1))*(sum(3)/sum(1)) +
     &      2.d0*(sum(3)/sum(1))*(sum(3)/sum(1))*
     &      (sum(3)/sum(1))
          sum(9) = sum(9)/sum(1) -
     &      (sum(6)/sum(1))*(sum(2)/sum(1)) -
     &      2.d0*(sum(5)/sum(1))*(sum(3)/sum(1)) +
     &      2.d0*(sum(3)/sum(1))*(sum(3)/sum(1))*
     &      (sum(2)/sum(1))
          sum(8) = sum(8)/sum(1) -
     &      (sum(4)/sum(1))*(sum(3)/sum(1)) -
     &      2.d0*(sum(5)/sum(1))*(sum(2)/sum(1)) +
     &      2.d0*(sum(3)/sum(1))*(sum(2)/sum(1))*
     &      (sum(2)/sum(1))
          sum(7) = sum(7)/sum(1) -
     &      3.d0*(sum(4)/sum(1))*(sum(2)/sum(1)) +
     &      2.d0*(sum(2)/sum(1))*(sum(2)/sum(1))*
     &      (sum(2)/sum(1))
          sum(6) = sum(6)/sum(1) -
     &      (sum(3)/sum(1))*(sum(3)/sum(1))
          sum(5) = sum(5)/sum(1) -
     &      (sum(2)/sum(1))*(sum(3)/sum(1))
          sum(4) = sum(4)/sum(1) -
     &      (sum(2)/sum(1))*(sum(2)/sum(1))
        endif
        sum(3) = sum(3)/sum(1)
        sum(2) = sum(2)/sum(1)
C        Multiply sum by fdf where
C        fdf = f/(1+f) (g*(1+g))^3/2 ((1+f)*(1+g))^(-morder)
        if(i.le.2) then
          sum(10) = sum(10) + (1.5d0-dble(morder))*dug2
          sum(7) = sum(7) - dble(morder+1)*duf2
          sum(6) = sum(6) + (1.5d0-dble(morder))*dug
          sum(4) = sum(4) - dble(morder+1)*duf
        endif
        sum(3) = sum(3) + 1.5d0 + (1.5d0-dble(morder))*ug
        sum(2) = sum(2) + 1.d0 - dble(morder+1)*uf
        sum(1) = fdf*sum(1)
C        Transform to tcl derivative from ln g 
C        derivative recalling that dlng/dlntc = 1, and dlng/dlnf = 0.5*uf.
        if(i.le.2) then
          sum(4) = sum(4) +
     &      0.5d0*(duf*sum(3) + uf*sum(5)) +
     &      0.5d0*uf*(sum(5) + 0.5d0*uf*sum(6))
          sum(7) = sum(7) +
     &      0.5d0*(duf2*sum(3) + 2.d0*duf*sum(5) + uf*sum(8)) +
     &      0.5d0*(duf*sum(5) + uf*sum(8) +
     &      0.5d0*uf*(2.d0*duf*sum(6) + uf*sum(9))) +
     &      0.5d0*uf*(sum(8) +
     &      0.5d0*(duf*sum(6) + uf*sum(9)) +
     &      0.5d0*uf*(sum(9) + 0.5d0*uf*sum(10)))
C          do sum(5) later than sum(7) so that rhs of sum(7) is undisturbed.
          sum(5) = sum(5) + 0.5d0*uf*sum(6)
          sum(8) = sum(8) + 
     &      0.5d0*(duf*sum(6) + uf*sum(9)) +
     &      0.5d0*uf*(sum(9) + 0.5d0*uf*sum(10))
          sum(9) = sum(9) + 0.5d0*uf*sum(10)
        endif
        sum(2) = sum(2) + 0.5d0*uf*sum(3)
        if(i.eq.1) then
          do ideriv = 1,9
            rhostar(ideriv) = sum(ideriv)
          enddo
        elseif(i.eq.2) then
          do ideriv = 1,9
            pstar(ideriv) = sum(ideriv)
          enddo
        elseif(i.eq.3) then
          do ideriv = 1,3
            sstar(ideriv) = sum(ideriv)
          enddo
        endif
      enddo
C       psi = fl+2.d0*(wf-log(1.d0+wf))
C       transform to new functions of fl, tcl.
C       N.B. up to now sstar has stored qstar and its derivatives, transform
C       to derivative of ln se wrt ln f and ln tc.
C       N.B. qstar in this programme is defined with different convention
C       than paper, see initial commentary.
      qconst = sstar(1)/(rhostar(1)*wf)
      sstar(1) = qconst+log((1.d0+wf)*(1.d0+wf)/f)
      sstar(2) = (qconst*(sstar(2)-rhostar(2)-0.5d0*uf)-1.d0/wf)/
     &  sstar(1)
      sstar(3) = (qconst*(sstar(3)-rhostar(3)))/sstar(1)
      pconst = pstar(1)*wf/rhostar(1)
      ustar(1) = sstar(1)+psi-pconst
      ustar(2) = (sstar(1)*sstar(2)+wf-
     &  pconst*(pstar(2)+0.5d0*uf-rhostar(2)))/
     &  ustar(1)
      ustar(3) = (sstar(1)*sstar(3)-pconst*(pstar(3)-rhostar(3)))/
     &  ustar(1)
C       transform pstar = pstar*g
      pstar(1) = (tc*wf)*pstar(1)
      pstar(2) = pstar(2) + 0.5d0*uf
      pstar(3) = pstar(3) + 1.d0
      pstar(4) = pstar(4) + 0.5d0*duf
      pstar(7) = pstar(7) + 0.5d0*duf2
      end
