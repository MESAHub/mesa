C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fermi_dirac_exchange.f 370 2006-11-29 23:57:39Z airwin $
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
      subroutine fermi_dirac_exchange(fl, tl,
     &  rhostar, pstar, sstar, ustar,
     &  nstar, morder,
     &  ifexchange_in, dve, dvef, dvet)
C       calculate transformed fermi_dirac integrals
C       rhostar, pstar, sstar, ustar as a function of fl and tl.
C       also calculate exchange if ifexchange_in != 0.
C      input data:
C      if abs(ifexchange_in) > 100 then_use_linear transform approximation
C        for exchange treatement.  Otherwise,_use_numerical transform as
C        described in research note.
C      ifexchange = mod(ifexchange_in,100):
C      One other wrinkle on deciding ifexchange is that large degeneracy
C      approximations (mod(ifexchange,10) = 2) only allowed above psi_lim.

C      ifexchange:
C      To understand ifexchange, must summarize free-energy model of fex used
C      in research note.  In general from kapusta relation,
C      fex is proportional to I - J + 2^1.5 pi^2/3 beta^2 K, where
C      K and J and known integrals which are functions of psi and beta,
C      and I = K^2. (see Kovetz et al 1972, ApJ 174, 109, hereafter KLVH)
C       0  < ifexchange < 10 --> KLVH treatment (i.e., drop K term).
C      10 <= ifexchange < 20 --> Kapusta treatment (i.e., retain K term).
C      20 <= ifexchange < 30 --> special test case,_use_I term alone.
C      30 <= ifexchange < 40 --> special test case,_use_J term alone.
C      40 <= ifexchange < 50 --> special test case,_use_K term alone.
C      N.B. if ifexchange is negative, then_use_non-relativistic limit
C       of corresponding positive ifexchange option.
C      ifexchange details:
C      ifexchange = 1 --> G(psi) + weak relativistic correction from KLVH
C      ifexchange = 2 --> I-J degenerate expression from KLVH (with corrected
C        sign error on a2 and high numerical precision a2 and a3, see
C        paper IV)
C      ifexchange = 11 or 12 K term added (both series from CG).
C      ifexchange = 21 or 22 I term alone (both series from KLVH)
C      ifexchange = 31 or 32 -J term alone (both series from KLVH).
C      ifexchange = 41 or 42 K term alone.
C      mod(ifexchange,10) = 4 is lowest order fit of J, K
C      mod(ifexchange,10) = 5 is next higher order fit of J, K
C      mod(ifexchange,10) = 6 is highest order fit of J, K
      implicit none
      include 'constants.h'
      integer ifexchange_in, morder, nstar
      double precision fl, tl,
     &  dve, dvef, dvet,
     &  pex, pext, pexf,
     &  sex, sexf, sext, uex
      double precision
     &  rhostar(nstar+6), pstar(nstar+6), sstar(nstar), ustar(nstar),
     &  n_e, t, tcl,
     &  fex, fexf, fext, fexf2, fexft, fext2,
     &  muex, muexf, muext
C       variables for full-order approach
      double precision psi_lim, psi, dpsidf, dpsidf2,
     &  psiprime, dpsiprimedf, dpsiprimedf2,
     &  dflprime, flprime, n_eprime,
     &  rhostarprime(9), pstarprime(9), sstarprime(3), ustarprime(3),
     &  fzero, fzerox, fzerof, fzerot, fzerox2, fzeroxf, fzeroxt,
     &  fzerof2, fzeroft, fzerot2,
     &  flprimef, flprimet, flprimef2, flprimeft, flprimet2,
     &  p_e, p_eprime, 
     &  fexprime, fexprimef, fexprimet,
     &  fexprimef2, fexprimeft, fexprimet2,
     &  fex_psi, fex_psif, fex_psit,
     &  fex_psif2, fex_psift, fex_psit2,
     &  muex2, muex2f, muex2t,
     &  fex2, fex2f, fex2t, fex2f2, fex2ft, fex2t2,
     &  sex2, free_ex, free_exf
      integer icount, ifexchange, iforder
      logical*4 ifcontinue
      save
      if(ifexchange_in.ne.0) then
C        iforder = 1 means exchange using first-order transformation from gcpf
C        iforder = 2 means exchange using numerical transformation from gcpf
        if(abs(ifexchange_in).ge.100) then
          iforder = 1
        else
          iforder = 2
        endif
        ifexchange = mod(ifexchange_in,100)
C        validity limit for Kovetz et al strong degeneracy series.
        psi_lim = max(3.d0,0.5d0/(ct*exp(tl)))
      else
C        signal that exchange not wanted.
        iforder = 0
      endif
C       always calculate fermi_dirac (and return values) for fl.
C       later may require primed values for non-linear transform
C       approach (iforder = 2) which correspond to flprime.
      tcl = tl + log(ct)
      call fermi_dirac(fl, tcl,
     &  rhostar, pstar, sstar, ustar,
     &  nstar, morder)
      if(iforder.eq.1) then
C         first-order transformation from gcpf form
C        large degeneracy approximation.
        if(mod(ifexchange,10).eq.2) then
          dpsidf = sqrt(1.d0+exp(fl))
          psi = fl + 2.d0*(dpsidf - log(1.d0 + dpsidf))
          if(psi.lt.psi_lim) then
            ifexchange = ifexchange -1
          endif
        endif
        call exchange_gcpf(
     &    ifexchange, fl, tl,
     &    psi, dpsidf, dpsidf2,
     &    fex, fexf, fext, fexf2, fexft, fext2,
     &    fex_psi, fex_psif, fex_psit,
     &    fex_psif2, fex_psift, fex_psit2)
C         number density of free electrons
        n_e = c_e*rhostar(1)
        t = exp(tl)
C         electron exchange chemical potential (divided by kt)
C           = partial fex(T, n_e)/partial n_e/k T
        muex = fexf/(n_e*rhostar(2)*boltzmann*t)
        muexf = muex*(fexf2/fexf - rhostar(2) - rhostar(4)/rhostar(2))
        muext = muex*(fexft/fexf - rhostar(3) - rhostar(5)/rhostar(2)
     &    - 1.d0)
C         dve is *negative* chemical potential of electron/kT
        dve = -muex
        dvef = -muexf
        dvet = -muext
      elseif(iforder.eq.2) then
C         full order approach (ignoring differential Coulomb effects)
C         start with non-relativistic degeneracy treatment if considering
C         option of strong degeneracy treatment
        if(mod(ifexchange,10).eq.2) then
          ifexchange = ifexchange -1
        endif
        call exchange_gcpf(
     &    ifexchange, fl, tl,
     &    psi, dpsidf, dpsidf2,
     &    fex, fexf, fext, fexf2, fexft, fext2,
     &    fex_psi, fex_psif, fex_psit,
     &    fex_psif2, fex_psift, fex_psit2)
C         number density of free electrons
        n_e = c_e*rhostar(1)
        t = exp(tl)
C        needed for exchange_end and exchange_free entries later.
        p_e = cpe*pstar(1)
C         solve for fl' such that
C         n_e(fl') - n_e(fl) - (1/kT) partial fex(psi(fl'),T)/partial psi = 0
C         where fex is the first-order free energy/volume = -kT ln Z_x
C         first approximation for delta fl = fl' - fl (see below)
C         function (see below) evaluated at fl' = fl
        fzero = - fex_psi/(boltzmann*t)
C         first derivative wrt x = flprime (see below) evaluated at fl' = fl
        fzerox = n_e*rhostar(2)
     &    - fex_psif/(boltzmann*t)
C         second derivative wrt x = flprime (see below) evaluated at fl' = fl
        fzerox2 = 
     &    n_e*(rhostar(2)*rhostar(2) + rhostar(4))
     &    - fex_psif2/(boltzmann*t)
C         first linear estimate
        dflprime = -fzero/fzerox
C         could solve quadratic, but_use_instead linear estimate of dflprime
C         to estimate quadratic correction and store temporarily in fzerox2
C         variable.
        fzerox2 = -0.5d0*dflprime*(dflprime*fzerox2/fzerox)
C         make second-orer correction if significantly less than first-order
C         estimate.
        if(abs(fzerox2).lt.0.5d0*abs(dflprime))
     &    dflprime = dflprime + fzerox2
        flprime = fl
        ifcontinue = .true.
        icount = 0
        do while(ifcontinue)
C           go one iteration beyond 1.d-7 which should give ~1.d-14 accuracy
C           for quadratic convergence
          icount = icount + 1
!          write(0,*) 'icount, flprime, dflprime = ',
!     &      icount, flprime, dflprime
          ifcontinue = abs(dflprime).gt.1.d-7.and.icount.le.20
C           fzero is monotonic in flprime so assured convergence so long
C           as don't change flprime too wildly.  So far experience is
C           that essentially no limits on change are fine so long as
C           don't cross discontinuity at fl_limit
          flprime = flprime + max(-1.d1,min(1.d1,dflprime))
          call fermi_dirac(flprime, tcl,
     &      rhostarprime, pstarprime, sstarprime, ustarprime,
     &      nstar, morder)
          call exchange_gcpf(
     &      ifexchange, flprime, tl,
     &      psiprime, dpsiprimedf, dpsiprimedf2,
     &      fexprime, fexprimef, fexprimet,
     &      fexprimef2, fexprimeft, fexprimet2,
     &      fex_psi, fex_psif, fex_psit,
     &      fex_psif2, fex_psift, fex_psit2)
          n_eprime = c_e*rhostarprime(1)
C           n_e(fl') - n_e(fl) - (1/kT) partial fex(psi(fl'),T)/partial psi = 0
C           where fex is the first-order free energy/volume = -kT ln Z_x
          fzero = n_eprime - n_e - fex_psi/(boltzmann*t)
C           derivative wrt x = flprime
          fzerox = n_eprime*rhostarprime(2)
     &      - fex_psif/(boltzmann*t)
          dflprime = -fzero/fzerox
        enddo
        if(icount.gt.500) then
          write(0,*) 'fl, tl = ', fl, tl
          write(0,*) 'flprime, dflprime = ', flprime, dflprime
          stop
     &      'fermi_dirac_exchange: first iteration failed to converge'
        endif
        if(mod(mod(ifexchange_in,100),10).eq.2
     &    .and.psiprime.ge.psi_lim) then
          ifexchange = mod(ifexchange_in,100)
C           if considering the possibility of a strong degeneracy treatment
C             and if result for non-relativistic degeneracy is in appropriate
C             regime with psiprime.ge.psi_lim
C             then_use_strong degeneracy treatment
        call exchange_gcpf(
     &    ifexchange, fl, tl,
     &    psi, dpsidf, dpsidf2,
     &    fex, fexf, fext, fexf2, fexft, fext2,
     &    fex_psi, fex_psif, fex_psit,
     &    fex_psif2, fex_psift, fex_psit2)
C         number density of free electrons
        n_e = c_e*rhostar(1)
        t = exp(tl)
C         solve for fl' such that
C         n_e(fl') - n_e(fl) - (1/kT) partial fex(psi(fl'),T)/partial psi = 0
C         where fex is the first-order free energy/volume = -kT ln Z_x
C         first approximation for delta fl = fl' - fl (see below)
C         function (see below) evaluated at fl' = fl
        fzero = - fex_psi/(boltzmann*t)
C         first derivative wrt x = flprime (see below) evaluated at fl' = fl
        fzerox = n_e*rhostar(2)
     &    - fex_psif/(boltzmann*t)
C         second derivative wrt x = flprime (see below) evaluated at fl' = fl
        fzerox2 = 
     &    n_e*(rhostar(2)*rhostar(2) + rhostar(4))
     &    - fex_psif2/(boltzmann*t)
C         first linear estimate
        dflprime = -fzero/fzerox
C         could solve quadratic, but_use_instead linear estimate of dflprime
C         to estimate quadratic correction and store temporarily in fzerox2
C         variable.
        fzerox2 = -0.5d0*dflprime*(dflprime*fzerox2/fzerox)
C         make second-orer correction if significantly less than first-order
C         estimate.
        if(abs(fzerox2).lt.0.5d0*abs(dflprime))
     &    dflprime = dflprime + fzerox2
        flprime = fl
        ifcontinue = .true.
        icount = 0
        do while(ifcontinue)
C           go one iteration beyond 1.d-7 which should give ~1.d-14 accuracy
C           for quadratic convergence
          icount = icount + 1
!          write(0,*) 'icount, flprime, dflprime = ',
!     &      icount, flprime, dflprime
          ifcontinue = abs(dflprime).gt.1.d-7.and.icount.le.20
C           fzero is monotonic in flprime so assured convergence so long
C           as don't change flprime too wildly.  So far experience is
C           that essentially no limits on change are fine so long as
C           don't cross discontinuity at fl_limit
          flprime = flprime + max(-1.d1,min(1.d1,dflprime))
          call fermi_dirac(flprime, tcl,
     &      rhostarprime, pstarprime, sstarprime, ustarprime,
     &      nstar, morder)
          call exchange_gcpf(
     &      ifexchange, flprime, tl,
     &      psiprime, dpsiprimedf, dpsiprimedf2,
     &      fexprime, fexprimef, fexprimet,
     &      fexprimef2, fexprimeft, fexprimet2,
     &      fex_psi, fex_psif, fex_psit,
     &      fex_psif2, fex_psift, fex_psit2)
          n_eprime = c_e*rhostarprime(1)
C           n_e(fl') - n_e(fl) - (1/kT) partial fex(psi(fl'),T)/partial psi = 0
C           where fex is the first-order free energy/volume = -kT ln Z_x
          fzero = n_eprime - n_e - fex_psi/(boltzmann*t)
C           derivative wrt x = flprime
          fzerox = n_eprime*rhostarprime(2)
     &      - fex_psif/(boltzmann*t)
          dflprime = -fzero/fzerox
        enddo
        if(icount.gt.20) then
          write(0,*) 'fl, tl = ', fl, tl
          write(0,*) 'flprime, dflprime = ', flprime, dflprime
          stop
     &      'fermi_dirac_exchange: second iteration failed to converge'
        endif
        endif
C         flprime (and therefore all primed quantities) is an implicit
C         function of fl and tl.  Find the derivatives from the chain
C         rule for implicit functions.
        fzerof = -n_e*rhostar(2)
        fzerot = (n_eprime*rhostarprime(3) - n_e*rhostar(3))
     &    - (fex_psit - fex_psi)/(boltzmann*t)
        flprimef = -fzerof/fzerox
        flprimet = -fzerot/fzerox
C         second partial derivatives:
        fzerox2 = 
     &    n_eprime*(rhostarprime(2)*rhostarprime(2) + rhostarprime(4))
     &    - fex_psif2/(boltzmann*t)
        fzeroxf = 0.d0
        fzeroxt = 
     &    n_eprime*(rhostarprime(2)*rhostarprime(3) + rhostarprime(5))
     &    - (fex_psift - fex_psif)/(boltzmann*t)
        fzerof2 = -n_e*(rhostar(2)*rhostar(2) + rhostar(4))
        fzeroft = -n_e*(rhostar(2)*rhostar(3) + rhostar(5))
        fzerot2 =
     &    n_eprime*(rhostarprime(3)*rhostarprime(3) + rhostarprime(6))
     &    - n_e*(rhostar(3)*rhostar(3) + rhostar(6))
     &    - (fex_psit2 - 2.d0*fex_psit + fex_psi)
     &    /(boltzmann*t)
        flprimef2 = -(
     &    fzerox2*flprimef*flprimef +
     &    fzeroxf*flprimef + fzeroxf*flprimef +
     &    fzerof2
     &    )/fzerox
        flprimeft = -(
     &    fzerox2*flprimef*flprimet +
     &    fzeroxf*flprimet + fzeroxt*flprimef +
     &    fzeroft
     &    )/fzerox
        flprimet2 = -(
     &    fzerox2*flprimet*flprimet +
     &    fzeroxt*flprimet + fzeroxt*flprimet +
     &    fzerot2
     &    )/fzerox
C         free energy/volume =
C         = fex(psi') + kT(psi'-psi) n_e -kT (ln Z_e(psi')-ln Z_e(psi))/V
C         = fex(psi') + kT(psi'-psi) n_e +
C             (f_e(psi')-kt ne' psi') - (f_e(psi)-kt ne psi)
C         = fex(psi') + kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
C         = fex(psi') + kT(n_e-n_e') psi' + (f_e(psi') - f_e(psi))
C         where n_e = n_e(psi) and n_e' = n_e(psi')
C         chemical potential/kt from fex(psi')
        muex = fexprimef*flprimef/(n_e*rhostar(2)*boltzmann*t)
        muexf = muex*(
     &    fexprimef2*flprimef/fexprimef
     &    + flprimef2/flprimef - rhostar(2) - rhostar(4)/rhostar(2))
        muext = muex*(
     &    (fexprimef2*flprimet + fexprimeft)/fexprimef
     &    - 1.d0 +
     &    flprimeft/flprimef - rhostar(3) - rhostar(5)/rhostar(2))
C         chemical potential/kt from rest of free energy:
C         n.b. by definition we have
C         partial f_e(psi')/n_e' = kT psi'
C         partial f_e(psi)/n_e = kT psi
C         so there is much cancellation
C         muex2 = psiprime - psi + (ne-ne') partial psi'/partial ne
C         first all but psiprime-psi
        muex2 = -dpsiprimedf*flprimef/(n_e*rhostar(2))
        muex2f = muex2*(
     &    (dpsiprimedf2/dpsiprimedf)*flprimef +
     &    flprimef2/flprimef - rhostar(2) - rhostar(4)/rhostar(2))
        muex2t = muex2*(
     &    (dpsiprimedf2/dpsiprimedf)*flprimet +
     &    flprimeft/flprimef - rhostar(3) - rhostar(5)/rhostar(2))
        muex2f = (n_eprime*rhostarprime(2)*flprimef - n_e*rhostar(2))*
     &    muex2 + (n_eprime - n_e)*muex2f
        muex2t = (n_eprime*(rhostarprime(2)*flprimet + rhostarprime(3))
     &    - n_e*rhostar(3))*muex2 + (n_eprime - n_e)*muex2t
        muex2 = (n_eprime - n_e)*muex2
C         now do psiprime - psi part
        muex2 = muex2 + psiprime - psi
        muex2f = muex2f + dpsiprimedf*flprimef - dpsidf
        muex2t = muex2t + dpsiprimedf*flprimet
C         dve is *negative* chemical potential of electron/kT
        dve = - muex - muex2
        dvef = - muexf - muex2f
        dvet = - muext - muex2t
      elseif(iforder.eq.0) then
        dve = 0.d0
        dvef = 0.d0
        dvet = 0.d0
      else
        stop 'fermi_dirac_exchange: should not happen'
      endif
      return
      entry exchange_pressure(rhostar, pstar, nstar,
     &  pex, pexf, pext)
      if(iforder.eq.1) then
C         electron exchange pressure
C           = - fex + ne (partial fex(T,ne)/partial ne)
        pex = -fex + n_e*boltzmann*t*muex
        pexf = -fexf + n_e*boltzmann*t*(rhostar(2)*muex + muexf)
        pext = -fext + n_e*boltzmann*t*((rhostar(3)+1.d0)*muex + muext)
      elseif(iforder.eq.2) then
C         full order approach (ignoring differential Coulomb effects)
C         free energy/volume =
C         = fex(psi') + kT(psi'-psi) n_e -kT (ln Z_e(psi')-ln Z_e(psi))/V
C         = fex(psi') + kT(psi'-psi) n_e +
C             (f_e(psi')-kt ne' psi') - (f_e(psi)-kt ne psi)
C         = fex(psi') + kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
C         = fex(psi') + kT(n_e-n_e') psi' + (f_e(psi') - f_e(psi))
C         where n_e = n_e(psi) and n_e' = n_e(psi')
C         first do just electron exchange pressure from fex(psi')
C         pex  = - fex' + ne (partial fex'(T,ne)/partial ne)
        pex = - fexprime + n_e*boltzmann*t*muex
        pexf = -fexprimef*flprimef
     &    + n_e*boltzmann*t*(rhostar(2)*muex + muexf)
        pext = -(fexprimef*flprimet + fexprimet)
     &    + n_e*boltzmann*t*((rhostar(3)+1.d0)*muex + muext)
C         remainder of exchange free energy in full-order approximation
C         fex2  = kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
        p_eprime = cpe*pstarprime(1)
        fex2 = boltzmann*t*(psiprime - psi)*n_e - (p_eprime - p_e)
        fex2f = boltzmann*t*((dpsiprimedf*flprimef - dpsidf)*n_e +
     &    (psiprime - psi)*n_e*rhostar(2)) -
     &    (p_eprime*pstarprime(2)*flprimef - p_e*pstar(2))
        fex2t = boltzmann*t*(
     &    (psiprime - psi + dpsiprimedf*flprimet)*n_e +
     &    (psiprime - psi)*n_e*rhostar(3)) -
     &    (p_eprime*(pstarprime(2)*flprimet + pstarprime(3)) -
     &    p_e*pstar(3))
C           electron exchange pressure
C             = - fex2 + ne (partial fex2(T,ne)/partial ne)
        pex = pex - fex2 + n_e*boltzmann*t*muex2
        pexf = pexf - fex2f + n_e*boltzmann*t*
     &    (rhostar(2)*muex2 + muex2f)
        pext = pext - fex2t + n_e*boltzmann*t*
     &    ((rhostar(3) + 1.d0)*muex2 + muex2t)
      elseif(iforder.eq.0) then
        pex = 0.d0
        pexf = 0.d0
        pext = 0.d0
      else
        stop 'exchange_pressure: should not happen'
      endif
      return
      entry exchange_free(rhostar, pstar, nstar,
     &  free_ex, free_exf)
      if(iforder.eq.1) then
        free_ex = fex
        free_exf = fexf
      elseif(iforder.eq.2) then
C         remainder of exchange free energy in full-order approximation
C         fex2  = kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
        p_eprime = cpe*pstarprime(1)
        fex2 = boltzmann*t*(psiprime - psi)*n_e - (p_eprime - p_e)
        fex2f = boltzmann*t*((dpsiprimedf*flprimef - dpsidf)*n_e +
     &    (psiprime - psi)*n_e*rhostar(2)) -
     &    (p_eprime*pstarprime(2)*flprimef - p_e*pstar(2))
        free_ex = fexprime + fex2
        free_exf = fexprimef*flprimef + fex2f
      elseif(iforder.eq.0) then
        free_ex = 0.d0
        free_exf = 0.d0
      else
        stop 'exchange_free: should not happen'
      endif
      return
      entry exchange_end(rhostar, pstar, nstar,
     &  pex, pext, pexf,
     &  sex, sexf, sext, uex)
      if(iforder.eq.1) then
C         electron exchange pressure
C           = - fex + ne (partial fex(T,ne)/partial ne)
        pex = -fex + n_e*boltzmann*t*muex
        pexf = -fexf + n_e*boltzmann*t*(rhostar(2)*muex + muexf)
        pext = -fext + n_e*boltzmann*t*((rhostar(3)+1.d0)*muex + muext)
C         electron exchange entropy per unit volume
C           = - partial fex(T,n_e)/partial T
C           = - (partial fex(tl,fl)/partial tl +
C                partial fex(tl,fl)/partial fl *
C                partial fl(tl,n_e)/partial tl)/T, where
C         partial fl(tl,n_e)/partial tl = -rhostar(3)/rhostar(2)
        sex = -(fext - fexf*rhostar(3)/rhostar(2))/t
        sexf = -(fexft
     &    - (fexf2/fexf + rhostar(5)/rhostar(3)
     &    - rhostar(4)/rhostar(2))
     &    *fexf*rhostar(3)/rhostar(2))/t
        sext = - sex - (fext2
     &    - (fexft/fexf + rhostar(6)/rhostar(3)
     &    - rhostar(5)/rhostar(2))
     &    *fexf*rhostar(3)/rhostar(2))/t
        uex = fex + t*sex
      elseif(iforder.eq.2) then
C         full order approach (ignoring differential Coulomb effects)
C         free energy/volume =
C         = fex(psi') + kT(psi'-psi) n_e -kT (ln Z_e(psi')-ln Z_e(psi))/V
C         = fex(psi') + kT(psi'-psi) n_e +
C             (f_e(psi')-kt ne' psi') - (f_e(psi)-kt ne psi)
C         = fex(psi') + kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
C         = fex(psi') + kT(n_e-n_e') psi' + (f_e(psi') - f_e(psi))
C         where n_e = n_e(psi) and n_e' = n_e(psi')
C         first do just electron exchange pressure from fex(psi')
C         pex  = - fex' + ne (partial fex'(T,ne)/partial ne)
        pex = - fexprime + n_e*boltzmann*t*muex
        pexf = -fexprimef*flprimef
     &    + n_e*boltzmann*t*(rhostar(2)*muex + muexf)
        pext = -(fexprimef*flprimet + fexprimet)
     &    + n_e*boltzmann*t*((rhostar(3)+1.d0)*muex + muext)
C         electron exchange entropy per unit volume
C           = - partial fex'(T,n_e)/partial T
C           = - (partial fex'(tl,fl)/partial tl +
C                partial fex'(tl,fl)/partial fl *
C                partial fl(tl,n_e)/partial tl)/T, where
C         partial fl(tl,n_e)/partial tl = -rhostar(3)/rhostar(2)
        sex = -(fexprimef*flprimet + fexprimet
     &    - fexprimef*flprimef*rhostar(3)/rhostar(2))/t
        sexf = -(fexprimef2*flprimef*flprimet + fexprimef*flprimeft
     &    + fexprimeft*flprimef
     &    - (fexprimef2*flprimef/fexprimef + flprimef2/flprimef
     &    + rhostar(5)/rhostar(3)
     &    - rhostar(4)/rhostar(2))
     &    *fexprimef*flprimef*rhostar(3)/rhostar(2))/t
        sext = - sex - (
     &    (fexprimef2*flprimet + fexprimeft)*flprimet
     &    + fexprimef*flprimet2
     &    + fexprimeft*flprimet + fexprimet2
     &    - ((fexprimef2*flprimet + fexprimeft)/fexprimef
     &    + flprimeft/flprimef
     &    + rhostar(6)/rhostar(3)
     &    - rhostar(5)/rhostar(2))
     &    *fexprimef*flprimef*rhostar(3)/rhostar(2))/t
        uex = fexprime + t*sex
C         remainder of exchange free energy in full-order approximation
C         fex2  = kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
        p_eprime = cpe*pstarprime(1)
        fex2 = boltzmann*t*(psiprime - psi)*n_e - (p_eprime - p_e)
        fex2f = boltzmann*t*((dpsiprimedf*flprimef - dpsidf)*n_e +
     &    (psiprime - psi)*n_e*rhostar(2)) -
     &    (p_eprime*pstarprime(2)*flprimef - p_e*pstar(2))
        fex2t = boltzmann*t*(
     &    (psiprime - psi + dpsiprimedf*flprimet)*n_e +
     &    (psiprime - psi)*n_e*rhostar(3)) -
     &    (p_eprime*(pstarprime(2)*flprimet + pstarprime(3)) -
     &    p_e*pstar(3))
        fex2f2 = boltzmann*t*((dpsiprimedf2*flprimef*flprimef +
     &    dpsiprimedf*flprimef2 - dpsidf2)*n_e +
     &    2.d0*(dpsiprimedf*flprimef - dpsidf)*n_e*rhostar(2) +
     &    (psiprime - psi)*n_e*(rhostar(2)*rhostar(2) + rhostar(4))) -
     &    (p_eprime*(pstarprime(2)*pstarprime(2) + pstarprime(4))*
     &    flprimef*flprimef +
     &    p_eprime*pstarprime(2)*flprimef2 -
     &    p_e*(pstar(2)*pstar(2) + pstar(4)))
        fex2ft = boltzmann*t*(
     &    (dpsiprimedf*flprimef - dpsidf +
     &    dpsiprimedf2*flprimef*flprimet + dpsiprimedf*flprimeft)*n_e +
     &    (psiprime - psi + dpsiprimedf*flprimet)*n_e*rhostar(2) +
     &    (dpsiprimedf*flprimef - dpsidf)*n_e*rhostar(3) +
     &    (psiprime - psi)*n_e*(rhostar(2)*rhostar(3) + rhostar(5))) -
     &    (p_eprime*pstarprime(2)*flprimef*
     &    (pstarprime(2)*flprimet + pstarprime(3)) +
     &    p_eprime*(pstarprime(4)*flprimef*flprimet +
     &    pstarprime(2)*flprimeft + pstarprime(5)*flprimef) -
     &    p_e*(pstar(2)*pstar(3) + pstar(5)))
        fex2t2 = boltzmann*t*(
     &    (psiprime - psi + dpsiprimedf*flprimet)*n_e +
     &    (psiprime - psi)*n_e*rhostar(3) +
     &    (dpsiprimedf*flprimet +
     &    dpsiprimedf2*flprimet*flprimet + dpsiprimedf*flprimet2)*n_e +
     &    (psiprime - psi + dpsiprimedf*flprimet)*n_e*rhostar(3) +
     &    (dpsiprimedf*flprimet)*n_e*rhostar(3) +
     &    (psiprime - psi)*n_e*(rhostar(3)*rhostar(3) + rhostar(6))) -
     &    (p_eprime*(pstarprime(2)*flprimet + pstarprime(3))*
     &    (pstarprime(2)*flprimet + pstarprime(3)) +
     &    p_eprime*(pstarprime(4)*flprimet*flprimet +
     &    pstarprime(5)*flprimet +
     &    pstarprime(2)*flprimet2 +
     &    pstarprime(5)*flprimet + pstarprime(6)) -
     &    p_e*(pstar(3)*pstar(3) + pstar(6)))
C           electron exchange pressure
C             = - fex2 + ne (partial fex2(T,ne)/partial ne)
        pex = pex - fex2 + n_e*boltzmann*t*muex2
        pexf = pexf - fex2f + n_e*boltzmann*t*
     &    (rhostar(2)*muex2 + muex2f)
        pext = pext - fex2t + n_e*boltzmann*t*
     &    ((rhostar(3) + 1.d0)*muex2 + muex2t)
C       electron exchange entropy per unit volume
C         = - partial fex2(T,n_e)/partial T
C         = - (partial fex2(tl,fl)/partial tl +
C              partial fex2(tl,fl)/partial fl *
C              partial fl(tl,n_e)/partial tl)/T, where
C       partial fl(tl,n_e)/partial tl = -rhostar(3)/rhostar(2)
        sex2 = -(fex2t - fex2f*rhostar(3)/rhostar(2))/t
        sex = sex + sex2
        sexf = sexf - (fex2ft -(fex2f2*rhostar(3) +
     &    fex2f*(rhostar(5) - rhostar(3)*rhostar(4)/rhostar(2)))/
     &    rhostar(2))/t
        sext = sext - sex2 -
     &    (fex2t2 -(fex2ft*rhostar(3) +
     &    fex2f*(rhostar(6) - rhostar(3)*rhostar(5)/rhostar(2)))/
     &    rhostar(2))/t
C           electron exchange internal energy per unit volume
C           from u = f + T s
        uex = uex + fex2 + t*sex2
      elseif(iforder.eq.0) then
        pex = 0.d0
        pexf = 0.d0
        pext = 0.d0
        sex = 0.d0
        sexf = 0.d0
        sext = 0.d0
        uex = 0.d0
      else
        stop 'exchange_end: should not happen'
      endif
      end
