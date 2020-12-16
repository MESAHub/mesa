C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: exchange_gcpf.f 370 2006-11-29 23:57:39Z airwin $
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
      subroutine exchange_gcpf(
     &  ifexchange, fl, tl,
     &  psi, dpsidf, dpsidf2,
     &  fex, fexf, fext,
     &  fexf2, fexft, fext2,
     &  fex_psi, fex_psif, fex_psit,
     &  fex_psif2, fex_psift, fex_psit2)

C      Calculate psi equivalent to CG eta and first and second
C       derivatives wrt fl.

C      Also, calculate
C      free energy/V in linear inversion approximation. this free
C      energy is given by fex = -kT ln Z_x/V, where ln Z_x is the grand
C      canonical partition function (gcpf) exchange term described in the
C      research note.  fex is used either to get the free energy/V in the
C      linear inversion approximation or else is used by the calling routine
C      to solve the full inversion problem numerically.

C      Also calculate higher order partial derivatives of fex wrt
C       fl, tl, fl2, fltl, tl2.

C      Also calculate fex_psi = partial fex(t,psi) wrt psi and higher order
C       partial derivatives of fex_psi wrt to fl, tl, fl2, fltl, tl2.
C
C      input data:
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
C      fl = ln f where f is EFF degeneracy parameter
C      tl = ln T.
C      output data:
C      psi = CG eta degeneracy parameter and derivatives wrt fl.
C      fex (explained above) and various partial derivates.
      implicit none
      include 'constants.h'
      integer ifexchange
      double precision fl, tl, 
     &  fex, fexf, fext,
     &  fexf2, fexft, fext2,
     &  fex_psi, fex_psif, fex_psit,
     &  fex_psif2, fex_psift, fex_psit2,
     &  con_exchange, fex0,
     &  fex0f, fex0g, fex0ff, fex0fg, fex0gg,
     &  fex0fff, fex0ffg, fex0fgg, fex0ggg,
     &  fex0log, fex0log2, d1fex0log2, d2fex0log2, d3fex0log2,
     &  f, t, psi,
     &  fpsi, dfpsi, d2fpsi, d3fpsi,
     &  dfpsidf, dfpsidf2, ddfpsidf, ddfpsidf2, 
     &  dpsidf, dpsidf2, dpsidf3
C       variables needed for relativistic corrections
      double precision con_nr1_ratio, con_d1_ratio,
     &  beta, fermi_dirac_ct,
     &  fd32(6), fd12(5),
     &  fex1, fex1f, fex1t,
     &  fex1f2, fex1ft, fex1t2,
     &  fex1_psi, fex1_psif, fex1_psit,
     &  fex1_psif2, fex1_psift, fex1_psit2,
     &  fex2, fex2f, fex2t,
     &  fex2f2, fex2ft, fex2t2,
     &  fex2_psi, fex2_psif, fex2_psit,
     &  fex2_psif2, fex2_psift, fex2_psit2
      integer ideriv, maxderiv
      parameter (maxderiv=4)
      double precision fd(maxderiv), dfd(maxderiv-1)
      integer iffirst
      data iffirst/1/
C       large degeneracy case variables with parameters taken from KLVH et
C       al. paper.
      double precision a1, a2, con2, con4,
     &  z, x, xz, xz2, xz3,
     &  phi, phiz, 
     &  lnalpha, lnalphaz, lnalphaz2, lnalphaz3
!      parameter (a1 = 0.449d0)  !eq. 21
C       this coefficient is negative of that in paper, but gives vastly
C       improved agreement with G(psi) asympototic series
!      parameter (a2 = 0.504d0) !eq. 22 (n.b. negated!)
!      parameter (a1 = 0.44904d0) !Lamb and Graziani (1997, private communication)
!      parameter (a2 = 0.50475d0) !Lamb and Graziani (1997, private communication)
C      These two values the result of running
C      ./exchange_a1_direct
C      A1 =         0.4490427560
C      ./exchange_a2_direct
C      A2 =         0.5047526561
C      These results are derived using similar numerical integration used
C      for other exchange integral testing which should be good to
C      nominal precision of 10^{-9}.  Actually, various tests with
C      convergence criteria indicates above results are good to 1 digit
C      in last place.
      parameter (a1 = 0.4490427560d0)
      parameter (a2 = 0.5047526561d0)
      parameter (con2 = pi*pi/12.d0)  !eq. 8
      parameter (con4 = 7.d0*pi*pi*pi*pi/720.d0)  !eq.8
!      double precision con6
!      parameter (con6 = 31.d0*pi*pi*pi*pi*pi*pi/30240.d0)  !CG 24.182
C      variables needed for EFF-style summations and derivatives.
      integer mordermax, maxc, maxkind,
     &  ikind, ioffset, nderiv
      parameter (mordermax = 9)
      parameter (maxc = (mordermax+1)*(mordermax+1))
      parameter (maxkind = 6)
      parameter (nderiv = 9)
      double precision ccoeff(maxc, maxkind)
      integer mforder(maxkind), mgorder(maxkind), n_max
      double precision sum(1+nderiv)
      double precision
     &  fexi, fexif, fexit,
     &  fexif2, fexift, fexit2,
     &  fexi_psi, fexi_psif, fexi_psit,
     &  fexi_psif2, fexi_psift, fexi_psit2,
     &  fexj, fexjf, fexjt,
     &  fexjf2, fexjft, fexjt2,
     &  fexj_psi, fexj_psif, fexj_psit,
     &  fexj_psif2, fexj_psift, fexj_psit2,
     &  fexk, fexkf, fexkt,
     &  fexkf2, fexkft, fexkt2,
     &  fexk_psi, fexk_psif, fexk_psit,
     &  fexk_psif2, fexk_psift, fexk_psit2,
     &  fexl, fexlf, fexlt,
     &  fexlf2, fexlft, fexlt2,
     &  fexl_psi, fexl_psif, fexl_psit,
     &  fexl_psif2, fexl_psift, fexl_psit2
      double precision wf, vf, uf, duf, duf2, g, vg, ug, fdf, dug, dug2,
     &  power
      double precision acon, bcon, ccon, dcon, econ
      integer index_exchange
      logical iffexi, iffexj, iffexk
      data index_exchange /0/
      save
      if(iffirst.eq.1) then
        iffirst = 0
C         m_e * e * k/h^2
C         = (m_e*avogadro) * e * k^2/(h^2 * (k*avogadro))
        con_exchange = electron_mass*echarge*clight*clight/
     &    (c2*c2*cr)
C         in non-relativistic degenerate limit,
C         free energy/volume = con_exchange*t*t*G(psi)
C         thus, coefficient of G(psi) is con_exchange*t*t
        con_exchange = -8.d0*pi*con_exchange*con_exchange
C         coefficient of F_{1/2} is con_nr1_ratio*beta^1.5 times coefficient of
C         G(psi).  Also, negative of ratio between I-J and K integral.
        con_nr1_ratio = -2.d0*pi*pi*sqrt(2.d0)/3.d0
C         coefficient of I-J is con_d1_ratio/t^2 times coefficient of G(psi)
C         where con_d1_ratio = -0.5*(m_e*c^2/k)^2.  This also works out from
C         eq. 6a of KLVH.
        con_d1_ratio = -0.5d0*(electron_mass*clight*clight/cr)**2
      endif
      t = exp(tl)
      beta = ct*t
      f = exp(fl)
C       partial psi/d ln f and higher derivatives
      dpsidf = sqrt(1.d0+f)
C       dpsidf2 = 0.5d0*f*dpsidf/(1.d0+f)
      dpsidf2 = 0.5d0*f/dpsidf
      dpsidf3 = dpsidf2*(1.d0 - dpsidf2/dpsidf)
      psi = fl + 2.d0*(dpsidf - log(1.d0 + dpsidf))
      if(mod(abs(ifexchange),10).eq.1) then
C        required for most cases so don't bother with logic to exclude
        fd32(1) = fermi_dirac_ct(psi, 0, 0)
        do ideriv = 1,4
          fd32(ideriv+1) = fermi_dirac_ct(psi, ideriv, 0)
          fd12(ideriv) = (2.d0/3.d0)*fd32(ideriv+1)
        enddo
        if(abs(ifexchange).eq.1.or.abs(ifexchange).eq.11.or.
     &      abs(ifexchange).eq.31) then
C          These cases require non-relativistic limit of J
          call f_psi(psi, fpsi, dfpsi, d2fpsi, d3fpsi)
          dfpsidf = dfpsi*dpsidf
          dfpsidf2 = d2fpsi*dpsidf*dpsidf + dfpsi*dpsidf2
          ddfpsidf = d2fpsi*dpsidf
          ddfpsidf2 = d3fpsi*dpsidf*dpsidf + d2fpsi*dpsidf2
C          definition of fps derived from definitions of fp12 and exct (see their
C          commentary):
C          fpsi(psi) = 
C          integral from -inf to psi of square of derivative of fp12(psi') d psi'
C          fp12(psi) is C&G F(1/2, psi)/Gamma(3/2), DeWitt's, script
C          capital I(1/2) and Rogers and DeWitt's script capital F(1/2).
C          Now, Gamma(3/2)^2 = pi/4, and d F(1/2, psi)/d psi = 1/2 F(-1/2, psi).
C          Thus, fpsi(psi) = G(psi)/pi where G(psi) is defined in paper.
C
C          electron exchange free energy per unit volume.
C          see DeWitt 1969 in "Low Luminosity Stars" ed. S. Kumar, eq. 17
C          for the canonical ensemble equation we are programming.
C          see also Rogers and DeWitt, 1973, Phys. Rev. A8, 1061 eq. 45
C          for the equivalent equation for the corresponding *grand*
C          canonical ensemble equation.
C          extra factor of pi needed because fpsi = G(psi)/pi
          fex0 = con_exchange*pi*t*t
          fex = fex0*fpsi
          fexf = fex0*dfpsidf
          fext = fex*2.d0
          fexf2 = fex0*dfpsidf2
          fexft = fexf*2.d0
          fext2 = fext*2.d0
          fex_psi = fex0*dfpsi
          fex_psif = fex0*ddfpsidf
          fex_psit = fex_psi*2.d0
          fex_psif2 = fex0*ddfpsidf2
          fex_psift = fex_psif*2.d0
          fex_psit2 = fex_psit*2.d0
        elseif(ifexchange.eq.-21) then
C          Term from Non-relativistic I = K^2  integral
C          extra factor is -4 beta^3/(2 beta^2)
          fex0 = con_exchange*t*t*(-2.d0*beta)
          fex = fex0*fd12(1)*fd12(1)
          fexf = 2.d0*fex0*fd12(1)*fd12(2)*dpsidf
          fext = 3.d0*fex 
          fexf2 = 2.d0*fex0*(
     &      (fd12(2)*fd12(2) + fd12(1)*fd12(3))*dpsidf*dpsidf +
     &      fd12(1)*fd12(2)*dpsidf2)
          fexft = 3.d0*fexf
          fext2 = 3.d0*fext
          fex_psi = 2.d0*fex0*fd12(1)*fd12(2)
          fex_psif = 2.d0*fex0*(
     &      fd12(2)*fd12(2) + fd12(1)*fd12(3))*dpsidf
          fex_psit = 3.d0*fex_psi
          fex_psif2 = 2.d0*fex0*(
     &      (3.d0*fd12(2)*fd12(3) + fd12(1)*fd12(4))*dpsidf*dpsidf +
     &      (fd12(2)*fd12(2) + fd12(1)*fd12(3))*dpsidf2)
          fex_psift = 3.d0*fex_psif
          fex_psit2 = 3.d0*fex_psit
        elseif(ifexchange.eq.-41) then
C          Term from non-relativistic K integral
C          extra factor is con_nr1_ratio*beta^2*(2 beta^1.5/(2 beta^2)
          fex0 = con_exchange*t*t*(con_nr1_ratio*beta**1.5d0)
          fex = fex0*fd12(1)
          fexf = fex0*fd12(2)*dpsidf
          fext = 3.5d0*fex 
          fexf2 = fex0*(fd12(3)*dpsidf*dpsidf + fd12(2)*dpsidf2)
          fexft = 3.5d0*fexf
          fext2 = 3.5d0*fext
          fex_psi = fex0*fd12(2)
          fex_psif = fex0*fd12(3)*dpsidf
          fex_psit = 3.5d0*fex_psi
          fex_psif2 = fex0*(fd12(4)*dpsidf*dpsidf + fd12(3)*dpsidf2)
          fex_psift = 3.5d0*fex_psif
          fex_psit2 = 3.5d0*fex_psit
        endif
C        N.B. nothing calculated in above blocks
C        for ifexchange = 21 or ifexchange = 41 cases.
        if(0.lt.ifexchange.and.ifexchange.lt.40) then
C          ifexchange = 1, 11, 21, 31
C          second Kapusta term is equivalent to KLVH so_use_their
C          weakly relativistic series.
C          free energy per unit volume in non-relativisitic case is
C          con_exchange*t*t*G(psi), and here we add on F_{1/2}F_{1/2} and
C          F_{1/2}F_{3/2} terms from KLVH eqs. 38, 36, and 36.  acon is
C          the overall multiplier compared to NR case, while bcon is
C          the ratio of the second to first term for relativistic part of
C          J-I, -I, and J.
          if(ifexchange.eq.1.or.ifexchange.eq.11) then
C            relativistic part of J-I
            acon = -1.5d0
            bcon = 0.75d0
          elseif(ifexchange.eq.21) then
C            -I
            acon = -2.d0
            bcon = 0.5d0
          elseif(ifexchange.eq.31) then
C            J
            acon = 0.5d0
            bcon = -0.25d0
          else
            stop 'exchange_gcpf: logic error 1'
          endif
          fex0 = acon*con_exchange*t*t*beta
C           n.b. multiply by fd12(1) later
          fex1 = fex0*(fd12(1) + bcon*beta*fd32(1))
          fex1f = fex0*(fd12(2) + bcon*beta*fd32(2))*dpsidf
          fex1t = 3.d0*fex1 + fex0*bcon*beta*fd32(1)
          fex1f2 = fex0
     &      *((fd12(3) + bcon*beta*fd32(3))*dpsidf*dpsidf
     &      + (fd12(2) + bcon*beta*fd32(2))*dpsidf2)
          fex1ft = 3.d0*fex1f + fex0*bcon*beta*fd32(2)*dpsidf
          fex1t2 = 3.d0*fex1t + 4.d0*fex0*bcon*beta*fd32(1)
          fex1_psi = fex0*(fd12(2) + bcon*beta*fd32(2))
          fex1_psif = fex0*(fd12(3) + bcon*beta*fd32(3))*dpsidf
          fex1_psit = 3.d0*fex1_psi + fex0*bcon*beta*fd32(2)
          fex1_psif2 = fex0
     &      *((fd12(4) + bcon*beta*fd32(4))*dpsidf*dpsidf
     &      + (fd12(3) + bcon*beta*fd32(3))*dpsidf2)
          fex1_psift = 3.d0*fex1_psif
     &      + fex0*bcon*beta*fd32(3)*dpsidf
          fex1_psit2 = 3.d0*fex1_psit
     &      + 4.d0*fex0*bcon*beta*fd32(2)
C           to get final result must multiply by fd12(1)
C           n.b. the order of assignment statements is important here
C           definition of fex_psi is partial of fex1 wrt psi
C           thus from multiplication transformation we have
C           fex1_psi = fex1_psi*fd12(1) + fex1*fd12(2)
C           **with the rhs evaluated before any multiplication transformation
C           is made**.
C           apply multiplication transformation to fex1_psi
          fex1_psif2 = fex1_psif2*fd12(1)
     &      + 2.d0*fex1_psif*fd12(2)*dpsidf
     &      + fex1_psi*(fd12(3)*dpsidf*dpsidf + fd12(2)*dpsidf2)
     &      + fex1f2*fd12(2) + 2.d0*fex1f*fd12(3)*dpsidf
     &      + fex1*(fd12(4)*dpsidf*dpsidf + fd12(3)*dpsidf2)
          fex1_psift = fex1_psift*fd12(1) + fex1_psit*fd12(2)*dpsidf
     &      + fex1ft*fd12(2) + fex1t*fd12(3)*dpsidf
          fex1_psit2 = fex1_psit2*fd12(1)
     &      + fex1t2*fd12(2)
          fex1_psif = fex1_psif*fd12(1) + fex1_psi*fd12(2)*dpsidf
     &      + fex1f*fd12(2) + fex1*fd12(3)*dpsidf
          fex1_psit = fex1_psit*fd12(1)
     &      + fex1t*fd12(2)
          fex1_psi = fex1_psi*fd12(1)
     &      + fex1*fd12(2)
C           apply multiplication transformation to fex1
          fex1f2 = fex1f2*fd12(1) + 2.d0*fex1f*fd12(2)*dpsidf
     &      + fex1*(fd12(3)*dpsidf*dpsidf + fd12(2)*dpsidf2)
          fex1ft = fex1ft*fd12(1) + fex1t*fd12(2)*dpsidf
          fex1t2 = fex1t2*fd12(1)
          fex1f = fex1f*fd12(1) + fex1*fd12(2)*dpsidf
          fex1t = fex1t*fd12(1)
          fex1 = fex1*fd12(1)
          if(ifexchange.ne.21) then
C            Add in previously calculated NR J term for J-I, and J cases.
            fex = fex + fex1
            fexf = fexf + fex1f
            fext = fext + fex1t
            fexf2 = fexf2 + fex1f2
            fexft = fexft + fex1ft
            fext2 = fext2 + fex1t2
            fex_psi = fex_psi + fex1_psi
            fex_psif = fex_psif + fex1_psif
            fex_psit = fex_psit + fex1_psit
            fex_psif2 = fex_psif2 + fex1_psif2
            fex_psift = fex_psift + fex1_psift
            fex_psit2 = fex_psit2 + fex1_psit2
          else
C            For -I case, no NR term previously calculated.
            fex = fex1
            fexf = fex1f
            fext = fex1t
            fexf2 = fex1f2
            fexft = fex1ft
            fext2 = fex1t2
            fex_psi = fex1_psi
            fex_psif = fex1_psif
            fex_psit = fex1_psit
            fex_psif2 = fex1_psif2
            fex_psift = fex1_psift
            fex_psit2 = fex1_psit2
          endif
        endif
        if(ifexchange.eq.11.or.ifexchange.eq.41) then
C          - K term.  CG, 24.289
          fex0 = con_exchange*con_nr1_ratio*t*t*beta**1.5d0
          fex2 = fex0*(fd12(1) + 0.25d0*beta*fd32(1))
          fex2f = fex0*(fd12(2) + 0.25d0*beta*fd32(2))*dpsidf
          fex2t = 3.5d0*fex2 + fex0*0.25d0*beta*fd32(1)
          fex2f2 = fex0
     &      *((fd12(3) + 0.25d0*beta*fd32(3))*dpsidf*dpsidf
     &      + (fd12(2) + 0.25d0*beta*fd32(2))*dpsidf2)
          fex2ft = 3.5d0*fex2f + fex0*0.25d0*beta*fd32(2)*dpsidf
          fex2t2 = 3.5d0*fex2t + 4.5d0*fex0*0.25d0*beta*fd32(1)
          fex2_psi = fex0*(fd12(2) + 0.25d0*beta*fd32(2))
          fex2_psif = fex0*(fd12(3) + 0.25d0*beta*fd32(3))*dpsidf
          fex2_psit = 3.5d0*fex2_psi + fex0*0.25d0*beta*fd32(2)
          fex2_psif2 = fex0
     &      *((fd12(4) + 0.25d0*beta*fd32(4))*dpsidf*dpsidf
     &      + (fd12(3) + 0.25d0*beta*fd32(3))*dpsidf2)
          fex2_psift = 3.5d0*fex2_psif
     &      + fex0*0.25d0*beta*fd32(3)*dpsidf
          fex2_psit2 = 3.5d0*fex2_psit
     &      + 4.5d0*fex0*0.25d0*beta*fd32(2)
          if(ifexchange.eq.11) then
C            Add in - K term to previously calculated J-I term.
            fex = fex + fex2
            fexf = fexf + fex2f
            fext = fext + fex2t
            fexf2 = fexf2 + fex2f2
            fexft = fexft + fex2ft
            fext2 = fext2 + fex2t2
            fex_psi = fex_psi + fex2_psi
            fex_psif = fex_psif + fex2_psif
            fex_psit = fex_psit + fex2_psit
            fex_psif2 = fex_psif2 + fex2_psif2
            fex_psift = fex_psift + fex2_psift
            fex_psit2 = fex_psit2 + fex2_psit2
          elseif(ifexchange.eq.41) then
C            K term alone with no reference to previously
C            calculated NR component.
            fex = fex2
            fexf = fex2f
            fext = fex2t
            fexf2 = fex2f2
            fexft = fex2ft
            fext2 = fex2t2
            fex_psi = fex2_psi
            fex_psif = fex2_psif
            fex_psit = fex2_psit
            fex_psif2 = fex2_psif2
            fex_psift = fex2_psift
            fex_psit2 = fex2_psit2
          else
            stop 'exchange_gcpf: logic error 2'
          endif
        endif
      elseif(0.lt.ifexchange.and.mod(ifexchange,10).eq.2) then
C        ifexchange = 2, 12, 22, 32, 42
C        strong degeneracy
        if(psi.lt.3.d0)
     &    stop 'exchange_gcpf: psi too low for degenerate series'
        z = psi*beta
        if(z.lt.0.01d0) then
          write(0,*) 'exchange_gcpf: psi beta so low that '//
     &      'significance loss occurs for degenerate series.'
          write(0,*) 'exchange_gcpf: use weakly relativistic '//
     &      'series instead.'
          stop 'exchange_gcpf: psi beta too low'
        endif
        phi = 1.d0 + z
        phiz = 1.d0
        x = sqrt(2.d0*z + z*z)
        xz = phi/x
C        xz2 = 1.d0/x - phi*xz/x^2 = (x^2 - phi^2)/x^3 = -1.d0/x^3
        xz2 = -1.d0/(x*x*x)
        xz3 = (-3.d0*xz2/x)*xz
        lnalpha = x + phi
        lnalphaz = 1.d0 + xz
        lnalphaz2 = xz2
        lnalphaz3 = xz3
C        order is important here.
        lnalphaz3 = (lnalphaz3 - 3.d0*lnalphaz*lnalphaz2/lnalpha
     &    + 2.d0*lnalphaz*lnalphaz*lnalphaz/
     &    (lnalpha*lnalpha))/lnalpha
        lnalphaz2 = (lnalphaz2 - lnalphaz*lnalphaz/lnalpha)/lnalpha
        lnalphaz = lnalphaz/lnalpha
        lnalpha = log(lnalpha)
        if(ifexchange.lt.40) then
C          coefficient of I - J is
C          con_d1_ratio/t^2 * coefficient of G(psi)
          fex0 = con_d1_ratio*con_exchange
          if(ifexchange.eq.2.or.ifexchange.eq.12) then
C            I_1 - J_1 term from KLVH eq. 24.
            acon = 1.5d0
            bcon = -3.0d0
            ccon = 1.5d0
            dcon = 0.5d0
          elseif(ifexchange.eq.22) then
C            I_1 term from KLVH eq. 10.
            acon = 0.5d0
            bcon = -1.0d0
            ccon = 0.5d0
            dcon = 0.5d0
          elseif(ifexchange.eq.32) then
C            -J_1 term from KLVH eq. 17.
            acon = 1.0d0
            bcon = -2.0d0
            ccon = 1.0d0
            dcon = 0.0d0
          else
            stop 'exchange_gcpf: logic error 3'
          endif
          fex = fex0*(acon*lnalpha*lnalpha
     &      + x*(bcon*phi*lnalpha
     &      + x*(ccon + dcon*x*x)))
C          first store z derivatives in fex_psi, fex_psif, etc.
          fex_psi = fex0*(acon*2.d0*lnalpha*lnalphaz
     &      + x*bcon*(phiz*lnalpha + phi*lnalphaz)
     &      + xz*(bcon*phi*lnalpha
     &      + x*(2.d0*ccon + 4.d0*dcon*x*x)))
          fex_psif = fex0*(
     &      acon*2.d0*(lnalphaz*lnalphaz + lnalpha*lnalphaz2)
     &      + x*bcon*(2.d0*phiz*lnalphaz + phi*lnalphaz2)
     &      + 2.d0*xz*bcon*(phiz*lnalpha + phi*lnalphaz)
     &      + xz2*(bcon*phi*lnalpha
     &      + x*(2.d0*ccon + 4.d0*dcon*x*x))
     &      + xz*xz*(2.d0*ccon + 12.d0*dcon*x*x))
          fex_psif2 = fex0*(
     &      acon*2.d0*(3.d0*lnalphaz*lnalphaz2 + lnalpha*lnalphaz3)
     &      + x*bcon*(3.d0*phiz*lnalphaz2 + phi*lnalphaz3)
     &      + 3.d0*xz*bcon*(2.d0*phiz*lnalphaz + phi*lnalphaz2)
     &      + 3.d0*xz2*bcon*(phiz*lnalpha + phi*lnalphaz)
     &      + xz3*(bcon*phi*lnalpha
     &      + x*(2.d0*ccon + 4.d0*dcon*x*x))
     &      + 3.d0*xz*xz2*(2.d0*ccon + 12.d0*dcon*x*x)
     &      + xz*xz*xz*24.d0*dcon*x)
C          convert z derivatives to fl and tl derivatives
          fexf = fex_psi*beta*dpsidf
          fext = fex_psi*z
          fexf2 = fex_psif*beta*beta*dpsidf*dpsidf +
     &      fex_psi*beta*dpsidf2
          fexft = fex_psif*z*beta*dpsidf + fex_psi*beta*dpsidf
          fext2 = fex_psif*z*z + fex_psi*z
C          order of assignment statements is important
          fex_psit2 = fex_psif2*z*z*beta + 3.d0*fex_psif*z*beta
     &      + fex_psi*beta
          fex_psift = fex_psif2*z*beta*beta*dpsidf
     &      + 2.d0*fex_psif*beta*beta*dpsidf 
          fex_psif2 = fex_psif2*beta*beta*beta*dpsidf*dpsidf
     &      + fex_psif*beta*beta*dpsidf2
          fex_psit = fex_psif*z*beta + fex_psi*beta
          fex_psif = fex_psif*beta*beta*dpsidf
          fex_psi = fex_psi*beta
C          note the first term above
C          is the lowest order term in psi which for small x reduces to the
C          equivalent of the first term in the asymptotic expression for
C          G(psi) = 2 psi^2 + ...
C          now do beta^2 term (which brings in psi^0 terms)
          fex0 = fex0*4.d0*beta*beta
          if(ifexchange.eq.2.or.ifexchange.eq.12) then
C            I_2 - J_2 term from KLVH eq. 24.
            acon = -1.0d0
            bcon = 4.0d0
            ccon = 1.0d0
            dcon = -3.0d0
          elseif(ifexchange.eq.22) then
C            I_2 term from KLVH eq. 10.
            acon = 0.0d0
            bcon = 0.0d0
            ccon = 1.0d0
            dcon = -1.0d0
          elseif(ifexchange.eq.32) then
C            -J_2 term from KLVH eq. 20 and 22.
            acon = -1.0d0
            bcon = 4.0d0
            ccon = 0.0d0
            dcon = -2.0d0
          else
            stop 'exchange_gcpf: logic error 4'
          endif
C          I_2-J_2, KLVH eq. 24. excluding ln beta/2) term.
          fex1 = fex0*(
     &      + acon*(2.d0*a1 + a2) + con2*(
     &      + bcon*log(x) + ccon*(1.d0 + x*x) + dcon*phi*lnalpha/x))
C          first store z derivatives in fex1_psi, fex1_psif, etc.
          fex1_psi = fex0*(
     &      con2*(xz*(bcon/x + 2.d0*ccon*x - dcon*phi*lnalpha/(x*x))
     &      + dcon*(phiz*lnalpha + phi*lnalphaz)/x
     &      ))
          fex1_psif = fex0*(
     &      con2*(xz2*(bcon/x + 2.d0*ccon*x - dcon*phi*lnalpha/(x*x))
     &      + xz*xz*(-bcon/(x*x) + 2.d0*ccon
     &      + 2.d0*dcon*phi*lnalpha/(x*x*x))
     &      - 2.d0*dcon*xz*(phiz*lnalpha + phi*lnalphaz)/(x*x)
     &      + dcon*(2.d0*phiz*lnalphaz + phi*lnalphaz2)/x
     &      ))
          fex1_psif2 = fex0*(
     &      con2*(xz3*(bcon/x + 2.d0*ccon*x - dcon*phi*lnalpha/(x*x))
     &      + 3.d0*xz*xz2*(-bcon/(x*x) + 2.d0*ccon
     &      + 2.d0*dcon*phi*lnalpha/(x*x*x))
     &      - 3.d0*dcon*xz2*(phiz*lnalpha + phi*lnalphaz)/(x*x)
     &      + xz*xz*xz*(+2.d0*bcon/(x*x*x)
     &      - 6.d0*dcon*phi*lnalpha/(x*x*x*x))
     &      + 6.d0*dcon*xz*xz*(phiz*lnalpha + phi*lnalphaz)/(x*x*x)
     &      - 3.d0*dcon*xz*(2.d0*phiz*lnalphaz + phi*lnalphaz2)/(x*x)
     &      + dcon*(3.d0*phiz*lnalphaz2 + phi*lnalphaz3)/x
     &      ))
C          convert z derivatives to fl and tl derivatives
          fex1f = fex1_psi*beta*dpsidf
          if(ifexchange.eq.2.or.ifexchange.eq.12.or.
     &        ifexchange.eq.32) then
C            log(beta/2) term from -J_2, KLVH eq. 20.
            fex1 = fex1 + fex0*con2*(-2.d0*log(beta/2.d0))
            fex1t = fex1_psi*z + 2.d0*fex1 - fex0*2.d0*con2
          elseif(ifexchange.eq.22) then
            fex1t = fex1_psi*z + 2.d0*fex1
          else
            stop 'exchange_gcpf: logic error 5'
          endif
          fex1f2 = fex1_psif*beta*beta*dpsidf*dpsidf
     &      + fex1_psi*beta*dpsidf2
          fex1ft = fex1_psif*z*beta*dpsidf
     &      + 3.d0*fex1_psi*beta*dpsidf
          if(ifexchange.eq.2.or.ifexchange.eq.12.or.
     &        ifexchange.eq.32) then
            fex1t2 = fex1_psif*z*z + 3.d0*fex1_psi*z
     &        + 2.d0*fex1t - 4.d0*fex0*con2
          else
            fex1t2 = fex1_psif*z*z + 3.d0*fex1_psi*z
     &        + 2.d0*fex1t
          endif
C          order of assignment statements is important
          fex1_psit2 = fex1_psif2*z*z*beta + 7.d0*fex1_psif*z*beta
     &      + 9.d0*fex1_psi*beta
          fex1_psift = fex1_psif2*z*beta*beta*dpsidf
     &      + 4.d0*fex1_psif*beta*beta*dpsidf 
          fex1_psif2 = fex1_psif2*beta*beta*beta*dpsidf*dpsidf
     &    + fex1_psif*beta*beta*dpsidf2
          fex1_psit = fex1_psif*z*beta + 3.d0*fex1_psi*beta
          fex1_psif = fex1_psif*beta*beta*dpsidf
          fex1_psi = fex1_psi*beta
          fex = fex + fex1
          fexf = fexf + fex1f
          fext = fext + fex1t
          fexf2 = fexf2 + fex1f2
          fexft = fexft + fex1ft
          fext2 = fext2 + fex1t2
          fex_psi = fex_psi + fex1_psi
          fex_psif = fex_psif + fex1_psif
          fex_psit = fex_psit + fex1_psit
          fex_psif2 = fex_psif2 + fex1_psif2
          fex_psift = fex_psift + fex1_psift
          fex_psit2 = fex_psit2 + fex1_psit2
C          now do beta^4 term (which brings in psi^{-2} terms)
          fex0 = fex0*beta*beta
          if(ifexchange.eq.2.or.ifexchange.eq.12) then
C            I_3 - J_3 term from KLVH eq. 24.
            acon = 2.0d0
            bcon = 1.0d0
            ccon = 2.0d0
            dcon = 1.0d0
            econ = 3.0d0
          elseif(ifexchange.eq.22) then
C            I_3 term from KLVH eq. 10.
            acon = 2.0d0
            bcon = 0.0d0
            ccon = -1.0d0
            dcon = -1.0d0
            econ = 1.0d0
          elseif(ifexchange.eq.32) then
C            -J_3 term from KLVH eq. 20 and 22.
            acon = 0.0d0
            bcon = 1.0d0
            ccon = 3.0d0
            dcon = 2.0d0
            econ = 2.0d0
          else
            stop 'exchange_gcpf: logic error 6'
          endif
          fex2 = fex0*(
     &      con2*con2*(acon + (acon + bcon/(x*x))/(x*x))
     &      - 3.d0*con4/(x*x*x*x)*
     &      (ccon + dcon*x*x + econ*phi*lnalpha/x))
C          first store z derivatives in fex2_psi, fex2_psif, etc.
          fex2_psi = fex0*(
     &      - 2.d0*con2*con2*xz*(acon + 2.d0*bcon/(x*x))/(x*x*x)
     &      + 3.d0*con4*xz/(x*x*x*x*x)
     &      *(4.d0*ccon + 2.d0*dcon*x*x + 5.d0*econ*phi*lnalpha/x)
     &      - 3.d0*con4/(x*x*x*x)
     &      *(econ*(phiz*lnalpha + phi*lnalphaz)/x)
     &      )
          fex2_psif = fex0*(
     &      - 2.d0*con2*con2*(xz2*(acon + 2.d0*bcon/(x*x))/(x*x*x)
     &      + xz*xz*(-3.d0*acon - 10.d0*bcon/(x*x))/(x*x*x*x))
     &      + 3.d0*con4*xz2/(x*x*x*x*x)
     &      *(4.d0*ccon + 2.d0*dcon*x*x + 5.d0*econ*phi*lnalpha/x)
     &      + 3.d0*con4*xz*xz/(x*x*x*x*x*x)
     &      *(-20.d0*ccon - 6.d0*dcon*x*x - 30.d0*econ*phi*lnalpha/x)
     &      + 3.d0*con4*xz/(x*x*x*x*x)
     &      *(10.d0*econ*(phiz*lnalpha + phi*lnalphaz)/x)
     &      + 3.d0*con4/(x*x*x*x)
     &      *(-econ*(2.d0*phiz*lnalphaz + phi*lnalphaz2)/x)
     &      )
          fex2_psif2 = fex0*(
     &      - 2.d0*con2*con2*(xz3*(acon + 2.d0*bcon/(x*x))/(x*x*x)
     &      + 3.d0*xz*xz2*(-3.d0*acon - 10.d0*bcon/(x*x))/(x*x*x*x)
     &      + xz*xz*xz*(12.d0*acon + 60.d0*bcon/(x*x))/(x*x*x*x*x))
     &      + 3.d0*con4*xz3/(x*x*x*x*x)
     &      *(4.d0*ccon + 2.d0*dcon*x*x + 5.d0*econ*phi*lnalpha/x)
     &      + 9.d0*con4*xz*xz2/(x*x*x*x*x*x)
     &      *(-20.d0*ccon - 6.d0*dcon*x*x - 30.d0*econ*phi*lnalpha/x)
     &      + 3.d0*con4*xz2/(x*x*x*x*x)
     &      *(15.d0*econ*(phiz*lnalpha + phi*lnalphaz)/x)
     &      + 3.d0*con4*xz*xz/(x*x*x*x*x*x)
     &      *(-90.d0*econ*(phiz*lnalpha + phi*lnalphaz)/x)
     &      + 3.d0*con4*xz*xz*xz/(x*x*x*x*x*x*x)
     &      *(+120.d0*ccon +24.d0*dcon*x*x + 210.d0*econ*phi*lnalpha/x)
     &      + 3.d0*con4*xz/(x*x*x*x*x)
     &      *(15.d0*econ*(2.d0*phiz*lnalphaz + phi*lnalphaz2)/x)
     &      + 3.d0*con4/(x*x*x*x)
     &      *(-econ*(3.d0*phiz*lnalphaz2 + phi*lnalphaz3)/x)
     &      )
C          convert z derivatives to fl and tl derivatives
          fex2f = fex2_psi*beta*dpsidf
          fex2t = fex2_psi*z + 4.d0*fex2
          fex2f2 = fex2_psif*beta*beta*dpsidf*dpsidf
     &      + fex2_psi*beta*dpsidf2
          fex2ft = fex2_psif*z*beta*dpsidf
     &      + 5.d0*fex2_psi*beta*dpsidf
          fex2t2 = fex2_psif*z*z + 5.d0*fex2_psi*z
     &      + 4.d0*fex2t
C          order of assignment statements is important
          fex2_psit2 = fex2_psif2*z*z*beta + 11.d0*fex2_psif*z*beta
     &      + 25.d0*fex2_psi*beta
          fex2_psift = fex2_psif2*z*beta*beta*dpsidf
     &      + 6.d0*fex2_psif*beta*beta*dpsidf 
          fex2_psif2 = fex2_psif2*beta*beta*beta*dpsidf*dpsidf
     &      + fex2_psif*beta*beta*dpsidf2
          fex2_psit = fex2_psif*z*beta + 5.d0*fex2_psi*beta
          fex2_psif = fex2_psif*beta*beta*dpsidf
          fex2_psi = fex2_psi*beta
          fex = fex + fex2
          fexf = fexf + fex2f
          fext = fext + fex2t
          fexf2 = fexf2 + fex2f2
          fexft = fexft + fex2ft
          fext2 = fext2 + fex2t2
          fex_psi = fex_psi + fex2_psi
          fex_psif = fex_psif + fex2_psif
          fex_psit = fex_psit + fex2_psit
          fex_psif2 = fex_psif2 + fex2_psif2
          fex_psift = fex_psift + fex2_psift
          fex_psit2 = fex_psit2 + fex2_psit2
        endif
        if(ifexchange.eq.12.or.ifexchange.eq.42) then
C          add in first Kapusta term (with no KLVH equivalent)
C          expression for strong degeneracy
C          coefficient of F_{1/2} is
C          fex0 = con_exchange*con_nr1_ratio*t*t*beta**1.5d0
C          but F_{1/2} = (2 beta)^{-3/2}(x*phi + ...), hence coefficient
C          of x*phi is ...
          fex0 = con_exchange*con_nr1_ratio*t*t/(2.d0**1.5d0)
          fex1 = fex0*(x*phi - lnalpha)
C          temporary to show how to deal with severe significance loss
C          in the degenerate series for small x.  But comment out because
C          I don't want to deal with this consistently for all series
C          and all partial derivatives, and
C          the weakly relativistic series can be used instead for
C          the large degeneracy case.  In fact, because of the significance
C          loss problem, I stop this code above now if x < 1.d-1 and an
C          attempt is made to_use_the degenerate series.
!          if(x.lt.1.d-1) then
C            CG. eq. 24.191.
!            fex1 = fex0*2.d0*x**3*(
!     &        1.d0/3.d0 - x*x*(
!     &        1.d0/10.d0 - x*x*(
!     &        3.d0/56.d0 - x*x*(
!     &        5.d0/144.d0 - x*x*(
!     &        35.d0/1408.d0)))))
!          endif
C          first store z derivatives in fex1_psi, fex1_psif, etc.
          fex1_psi = fex0*(xz*phi + x*phiz - lnalphaz)
          fex1_psif = fex0*(xz2*phi + 2.d0*xz*phiz - lnalphaz2)
          fex1_psif2 = fex0*(xz3*phi + 3.d0*xz2*phiz - lnalphaz3)
C          convert z derivatives to fl and tl derivatives
C          copied all conversions from above because that fex0 term was also
C          just proportional to t**2 (excluding con2 term).
          fex1f = fex1_psi*beta*dpsidf
          fex1t = fex1_psi*z + 2.d0*fex1
          fex1f2 = fex1_psif*beta*beta*dpsidf*dpsidf
     &      + fex1_psi*beta*dpsidf2
          fex1ft = fex1_psif*z*beta*dpsidf
     &      + 3.d0*fex1_psi*beta*dpsidf
          fex1t2 = fex1_psif*z*z + 3.d0*fex1_psi*z
     &      + 2.d0*fex1t
C          order of assignment statements is important
          fex1_psit2 = fex1_psif2*z*z*beta + 7.d0*fex1_psif*z*beta
     &      + 9.d0*fex1_psi*beta
          fex1_psift = fex1_psif2*z*beta*beta*dpsidf
     &      + 4.d0*fex1_psif*beta*beta*dpsidf 
          fex1_psif2 = fex1_psif2*beta*beta*beta*dpsidf*dpsidf
     &      + fex1_psif*beta*beta*dpsidf2
          fex1_psit = fex1_psif*z*beta + 3.d0*fex1_psi*beta
          fex1_psif = fex1_psif*beta*beta*dpsidf
          fex1_psi = fex1_psi*beta
C          next higher term in beta
          fex0 = fex0*4.d0*con2*beta*beta
          fex2 = fex0*phi/x
C          first store z derivatives in fex1_psi, fex1_psif, etc.
          fex2_psi = fex0*(phiz - phi*xz/x)/x
          fex2_psif = fex0*(-2.d0*phiz*xz - phi*xz2
     &      + 2.d0*phi*xz*xz/x)/(x*x)
          fex2_psif2 = fex0*(3.d0*phiz*(-xz2 + 2.d0*xz*xz/x)
     &      + phi*(-xz3 + 6.d0*xz*(xz2 - xz*xz/x)/x))/(x*x)
C          convert z derivatives to fl and tl derivatives
          fex2f = fex2_psi*beta*dpsidf
          fex2t = fex2_psi*z + 4.d0*fex2
          fex2f2 = fex2_psif*beta*beta*dpsidf*dpsidf
     &      + fex2_psi*beta*dpsidf2
          fex2ft = fex2_psif*z*beta*dpsidf
     &      + 5.d0*fex2_psi*beta*dpsidf
          fex2t2 = fex2_psif*z*z + 5.d0*fex2_psi*z
     &      + 4.d0*fex2t
C          order of assignment statements is important
          fex2_psit2 = fex2_psif2*z*z*beta + 11.d0*fex2_psif*z*beta
     &      + 25.d0*fex2_psi*beta
          fex2_psift = fex2_psif2*z*beta*beta*dpsidf
     &      + 6.d0*fex2_psif*beta*beta*dpsidf 
          fex2_psif2 = fex2_psif2*beta*beta*beta*dpsidf*dpsidf
     &      + fex2_psif*beta*beta*dpsidf2
          fex2_psit = fex2_psif*z*beta + 5.d0*fex2_psi*beta
          fex2_psif = fex2_psif*beta*beta*dpsidf
          fex2_psi = fex2_psi*beta
C          add in lower order term
          fex2 = fex2 + fex1
          fex2f = fex2f + fex1f
          fex2t = fex2t + fex1t
          fex2f2 = fex2f2 + fex1f2
          fex2ft = fex2ft + fex1ft
          fex2t2 = fex2t2 + fex1t2
          fex2_psi = fex2_psi + fex1_psi
          fex2_psif = fex2_psif + fex1_psif
          fex2_psit = fex2_psit + fex1_psit
          fex2_psif2 = fex2_psif2 + fex1_psif2
          fex2_psift = fex2_psift + fex1_psift
          fex2_psit2 = fex2_psit2 + fex1_psit2
C          temporary.  Add more K terms to show how to increase precision
C          at large degeneracy at expense of precision for low degeneracy
C          for low psi beta.  However, this low psi beta region is handled
C          by weakly relativistic series in any case because of
C          significance loss in some of the expressions so leave all this
C          commented out now that I have proved the point with some plots.
!          fex0 = fex0/con2*beta*beta
!          fex2 = fex2 + fex0*(3.d0*con4*phi/x**5 +
!     &      15.d0*con6*beta*beta*phi*(7.d0 + 4.d0*x*x)/x**9)
          if(ifexchange.eq.12) then
            fex = fex + fex2
            fexf = fexf + fex2f
            fext = fext + fex2t
            fexf2 = fexf2 + fex2f2
            fexft = fexft + fex2ft
            fext2 = fext2 + fex2t2
            fex_psi = fex_psi + fex2_psi
            fex_psif = fex_psif + fex2_psif
            fex_psit = fex_psit + fex2_psit
            fex_psif2 = fex_psif2 + fex2_psif2
            fex_psift = fex_psift + fex2_psift
            fex_psit2 = fex_psit2 + fex2_psit2
          elseif(ifexchange.eq.42) then
C            test mode for K term alone.
            fex = fex2
            fexf = fex2f
            fext = fex2t
            fexf2 = fex2f2
            fexft = fex2ft
            fext2 = fex2t2
            fex_psi = fex2_psi
            fex_psif = fex2_psif
            fex_psit = fex2_psit
            fex_psif2 = fex2_psif2
            fex_psift = fex2_psift
            fex_psit2 = fex2_psit2
          else
            stop 'exchange_gcpf: logic error 7'
          endif
        endif
      elseif(4.le.abs(mod(ifexchange,10)).and.
     &    abs(mod(ifexchange,10)).le.6) then
C        Must provide fexi in these cases.
        iffexi = (0.le.abs(ifexchange).and.abs(ifexchange).lt.30)
C        Must provide fexj in all cases except fexi alone or fexk alone.
        iffexj = .not.(
     &      (20.le.abs(ifexchange).and.abs(ifexchange).lt.30).or.
     &    (40.le.abs(ifexchange).and.abs(ifexchange).lt.50))
C        Must provide fexk in these cases.
        iffexk = (10.le.abs(ifexchange).and.abs(ifexchange).lt.20).or.
     &    (40.le.abs(ifexchange).and.abs(ifexchange).lt.50)

        if(index_exchange.ne.mod(abs(ifexchange),10)-3) then
          index_exchange = mod(abs(ifexchange),10)-3
C          fill ccoeff, mforder, and mgorder if that hasn't been done before.
          call exchange_coeff(index_exchange, maxc, maxkind,
     &      ccoeff, mforder, mgorder)
        endif
        fexi = 0.d0
        fexif = 0.d0
        fexit = 0.d0
        fexif2 = 0.d0
        fexift = 0.d0
        fexit2 = 0.d0
        fexi_psi = 0.d0
        fexi_psif = 0.d0
        fexi_psit = 0.d0
        fexi_psif2 = 0.d0
        fexi_psift = 0.d0
        fexi_psit2 = 0.d0
        fexj = 0.d0
        fexjf = 0.d0
        fexjt = 0.d0
        fexjf2 = 0.d0
        fexjft = 0.d0
        fexjt2 = 0.d0
        fexj_psi = 0.d0
        fexj_psif = 0.d0
        fexj_psit = 0.d0
        fexj_psif2 = 0.d0
        fexj_psift = 0.d0
        fexj_psit2 = 0.d0
        fexk = 0.d0
        fexkf = 0.d0
        fexkt = 0.d0
        fexkf2 = 0.d0
        fexkft = 0.d0
        fexkt2 = 0.d0
        fexk_psi = 0.d0
        fexk_psif = 0.d0
        fexk_psit = 0.d0
        fexk_psif2 = 0.d0
        fexk_psift = 0.d0
        fexk_psit2 = 0.d0
        wf = sqrt(1.d0+f)
        vf = 1.d0/(1.d0+f)
        uf = f*vf
C      partial uf wrt ln f.  N.B. 1-uf = 1/(1+f)
        duf = uf/(1.d0 + f)
        duf2 = uf*(1.d0-f)/((1.d0 + f)*(1.d0 + f))
        if(mod(ifexchange,10).gt.0) then
          g = beta*wf
        else
          g = 0.d0
          if(iffexj) then
C            if fexj is provided it contains the NR component.
            iffexi = .false.
            iffexk = .false.
          endif
        endif
        vg = 1.d0/(1.d0+g)
        ug = g*vg
C        partial ug wrt ln f
        dug = ug/(1.d0+g)
        dug2 = ug*(1.d0-g)/((1.d0 + g)*(1.d0 + g))

        if(iffexj) then
C          Calculate approximation to J integral.
C          From Paper IV, eqs. 13 and 32
C          NR component = 2 beta^2 G(eta)/(1+g)^{N_max} =
C             2 pi beta^2 fpsi/(1+g)^{N_max} =
C             2 pi g^2/(1+f) fpsi/(1+g)^{N_max},
C          where N_max is the maximum of N+1, N'+1, N''+3, N'''+1.
          call f_psi(psi, fpsi, dfpsi, d2fpsi, d3fpsi)
          dfpsidf = dfpsi*dpsidf
          dfpsidf2 = d2fpsi*dpsidf*dpsidf + dfpsi*dpsidf2
          ddfpsidf = d2fpsi*dpsidf
          ddfpsidf2 = d3fpsi*dpsidf*dpsidf + d2fpsi*dpsidf2
          n_max = max(mgorder(1)+1, mgorder(2)+1,
     &      mgorder(3)+3, mgorder(4)+1)
C          From Paper IV, eqs. 13 and 32 and G(eta) = pi*fpsi
          fex0 = 2.d0*pi*beta*beta*vg**n_max
          fexj = fex0*fpsi
          fexjf = fex0*dfpsidf +
     &      fexj*(-0.5d0*n_max)*ug*uf
          fexjt = fexj*(2.d0 - n_max*ug)
          fexjf2 = fex0*dfpsidf2 +
     &      fex0*dfpsidf*(-0.5d0*n_max)*ug*uf +
     &      fexjf*(-0.5d0*n_max)*ug*uf +
     &      fexj*(-0.5d0*n_max)*(dug*uf*0.5d0*uf + ug*duf)
          fexjft = fexjf*(2.d0 - n_max*ug) +
     &      fexj*(-0.5d0*n_max)*(dug*uf)
          fexjt2 = fexjt*(2.d0 - n_max*ug) - fexj*n_max*dug
          fexj_psi = fex0*dfpsi
          fexj_psif = fexj_psi*(-0.5d0*n_max)*ug*uf +
     &      fex0*ddfpsidf
          fexj_psit = fexj_psi*(2.d0 - n_max*ug)
          fexj_psif2 = fexj_psi*(-0.5d0*n_max)*
     &      (dug*uf*0.5d0*uf + ug*duf) +
     &      fexj_psif*(-0.5d0*n_max)*ug*uf +
     &      fex0*ddfpsidf*(-0.5d0*n_max)*ug*uf +
     &      fex0*ddfpsidf2
          fexj_psift = fexj_psif*(2.d0 - n_max*ug) +
     &      fexj_psi*(-n_max*dug)*0.5d0*uf
          fexj_psit2 = fexj_psit*(2.d0 - n_max*ug) +
     &      fexj_psi*(-n_max*dug)
          fexj_psi = fexj_psi +
     &      fexj*(-0.5d0*n_max)*ug*uf/dpsidf
          fexj_psif = fexj_psif +
     &      fexjf*(-0.5d0*n_max)*ug*uf/dpsidf +
     &      fexj*(-0.5d0*n_max)*
     &      (dug*uf*0.5d0*uf + ug*duf - ug*uf*dpsidf2/dpsidf)/dpsidf
          fexj_psit = fexj_psit +
     &      fexjt*(-0.5d0*n_max)*ug*uf/dpsidf +
     &      fexj*(-0.5d0*n_max)*dug*uf/dpsidf
          fexj_psif2 = fexj_psif2 +
     &      fexjf2*(-0.5d0*n_max)*ug*uf/dpsidf +
     &      2.d0*fexjf*(-0.5d0*n_max)*
     &      (dug*uf*0.5d0*uf + ug*duf - ug*uf*dpsidf2/dpsidf)/dpsidf +
     &      fexj*(-0.5d0*n_max)*
     &      (dug2*uf*0.5d0*uf*0.5d0*uf + 1.5d0*dug*uf*duf + ug*duf2 -
     &      2.d0*(dug*uf*0.5d0*uf + ug*duf)*dpsidf2/dpsidf -
     &      ug*uf*(dpsidf3 - 2.d0*dpsidf2*dpsidf2/dpsidf)/dpsidf
     &      )/dpsidf
          fexj_psift = fexj_psift +
     &      fexjft*(-0.5d0*n_max)*ug*uf/dpsidf +
     &      fexjt*(-0.5d0*n_max)*
     &      (dug*uf*0.5d0*uf + ug*duf - ug*uf*dpsidf2/dpsidf)/dpsidf +
     &      fexjf*(-0.5d0*n_max)*dug*uf/dpsidf +
     &      fexj*(-0.5d0*n_max)*
     &      (dug2*uf*0.5d0*uf + dug*duf - dug*uf*dpsidf2/dpsidf
     &      )/dpsidf
          fexj_psit2 = fexj_psit2 +
     &      fexjt2*(-0.5d0*n_max)*ug*uf/dpsidf +
     &      2.d0*fexjt*(-0.5d0*n_max)*dug*uf/dpsidf +
     &      fexj*(-0.5d0*n_max)*dug2*uf/dpsidf
C          fdf = (f/(1+f))^2 g^2 (g/(1+g))
          fdf = uf*uf*g*g*ug
          do ikind = 1,4
            if(mforder(ikind).ge.0.and.mgorder(ikind).ge.0) then
C              Note, for negative mod(ifexchange,10), g is zero.
C              Nevertheless, for programming
C              simplicity effsum_calc grinds through entire sum.
              call effsum_calc(f, g, mforder(ikind), mgorder(ikind),
     &          ccoeff(1,ikind), nderiv, sum)
C              from Paper IV, equation 32
              fex0 = fdf*vf**mforder(ikind)*vg**mgorder(ikind)
C              ln fex0 = 2 ln f - 2 ln 1+f + 3 ln g - ln 1+g
C              -mforder ln 1+f -mgorder ln 1+g
              fex0f = fex0*(2.d0 - (2.d0 + mforder(ikind))*uf)
              fex0g = fex0*(3.d0 - (1.d0 + mgorder(ikind))*ug)
              fex0ff = fex0f*(2.d0 - (2.d0 + mforder(ikind))*uf) +
     &          fex0*(-(2.d0 + mforder(ikind)))*duf
              fex0fg = fex0g*(2.d0 - (2.d0 + mforder(ikind))*uf)
              fex0gg = fex0g*(3.d0 - (1.d0 + mgorder(ikind))*ug) +
     &          fex0*(-(1.d0 + mgorder(ikind)))*dug
              fex0fff = fex0ff*(2.d0 - (2.d0 + mforder(ikind))*uf) +
     &          2.d0*(fex0f*(-(2.d0 + mforder(ikind)))*duf) +
     &          fex0*(-(2.d0 + mforder(ikind)))*duf2
              fex0ffg = fex0fg*(2.d0 - (2.d0 + mforder(ikind))*uf) +
     &          fex0g*(-(2.d0 + mforder(ikind))*duf)
              fex0fgg = fex0gg*(2.d0 - (2.d0 + mforder(ikind))*uf)
              fex0ggg = fex0gg*(3.d0 - (1.d0 + mgorder(ikind))*ug) +
     &          2.d0*(fex0g*(-(1.d0 + mgorder(ikind)))*dug) +
     &          fex0*(-(1.d0 + mgorder(ikind)))*dug2
              if(ikind.eq.2) then
                if(g.gt.1.d-2) then
                  fex0log = log(1.d0+g)
                else
                  fex0log = g*(1.d0 - g*(1.d0/2.d0 - g*(1.d0/3.d0 -
     &              g*(1.d0/4.d0 - g*(1.d0/5.d0 - g*(1.d0/6.d0 -
     &              g*(1.d0/7.d0 - g*(1.d0/8.d0))))))))
                endif
                fex0ggg = fex0log*fex0ggg + 3.d0*ug*fex0gg +
     &            3.d0*dug*fex0g + dug2*fex0
                fex0fgg = fex0log*fex0fgg + 2.d0*ug*fex0fg +
     &            dug*fex0f
                fex0ffg = fex0log*fex0ffg + ug*fex0ff
                fex0fff = fex0log*fex0fff
                fex0gg = fex0log*fex0gg + 2.d0*ug*fex0g + dug*fex0
                fex0fg = fex0log*fex0fg + ug*fex0f
                fex0ff = fex0log*fex0ff
                fex0g = fex0log*fex0g + ug*fex0
                fex0f = fex0log*fex0f
                fex0 = fex0log*fex0
              elseif(ikind.eq.3) then
C                N.B. this branch normally not used, but I have tested
C                it anyhow with exchange_test and suitable swapping
C                of ikind=2 and 3.  There was substantial significance
C                loss, but I think all derivatives were okay, but should
C                retest if ever seriously_use_this branch.
                if(g.gt.1.d-2) then
                  fex0log = log(1.d0+g)
                else
                  fex0log = g*(1.d0 - g*(1.d0/2.d0 - g*(1.d0/3.d0 -
     &              g*(1.d0/4.d0 - g*(1.d0/5.d0 - g*(1.d0/6.d0 -
     &              g*(1.d0/7.d0 - g*(1.d0/8.d0))))))))
                endif
                fex0log2 = -fex0log*fex0log*vg*vg
                d1fex0log2 = -2.d0*(
     &            ug*(fex0log*vg*vg + fex0log2)
     &            )
                d2fex0log2 = -2.d0*(
     &            dug*(fex0log*vg*vg + fex0log2) +
     &            ug*(ug*vg*vg - 2.d0*fex0log*vg*vg*ug + d1fex0log2)
     &            )
                d3fex0log2 = -2.d0*(
     &            dug2*(fex0log*vg*vg + fex0log2) +
     &            2.d0*dug*(ug*vg*vg - 2.d0*fex0log*vg*vg*ug +
     &            d1fex0log2) +
     &            ug*(dug*vg*vg - 4.d0*ug*vg*vg*ug -
     &            2.d0*fex0log*(-2.d0*vg*vg*ug*ug + vg*vg*dug) +
     &            d2fex0log2)
     &            )

                fex0ggg = d3fex0log2*fex0 + 3.d0*d2fex0log2*fex0g +
     &            3.d0*d1fex0log2*fex0gg + fex0log2*fex0ggg
                fex0fgg = d2fex0log2*fex0f + 2.d0*d1fex0log2*fex0fg +
     &            fex0log2*fex0fgg
                fex0ffg = d1fex0log2*fex0ff + fex0log2*fex0ffg
                fex0fff = fex0log2*fex0fff
                fex0gg = d2fex0log2*fex0 + 2.d0*d1fex0log2*fex0g +
     &            fex0log2*fex0gg
                fex0fg = d1fex0log2*fex0f + fex0log2*fex0fg
                fex0ff = fex0log2*fex0ff
                fex0g = d1fex0log2*fex0 + fex0log2*fex0g
                fex0f = fex0log2*fex0f
                fex0 = fex0log2*fex0
              elseif(ikind.eq.4) then
                if(f.gt.1.d-2) then
                  fex0log = log(1.d0+f)
                else
                  fex0log = f*(1.d0 - f*(1.d0/2.d0 - f*(1.d0/3.d0 -
     &              f*(1.d0/4.d0 - f*(1.d0/5.d0 - f*(1.d0/6.d0 -
     &              f*(1.d0/7.d0 - f*(1.d0/8.d0))))))))
                endif
                fex0log2 = -fex0log/(1.d0+f)
                d1fex0log2 =
     &            fex0log2*uf*(1.d0/fex0log - 1.d0)
                d2fex0log2 =
     &            d1fex0log2*uf*(1.d0/fex0log - 1.d0) +
     &            fex0log2*(duf*(1.d0/fex0log - 1.d0) -
     &            uf*uf/(fex0log*fex0log))
                d3fex0log2 =
     &            d2fex0log2*uf*(1.d0/fex0log - 1.d0) +
     &            2.d0*d1fex0log2*(duf*(1.d0/fex0log - 1.d0) -
     &            uf*uf/(fex0log*fex0log)) +
     &            fex0log2*(duf2*(1.d0/fex0log - 1.d0) -
     &            3.d0*duf*uf/(fex0log*fex0log) +
     &            2.d0*uf*uf*uf/(fex0log*fex0log*fex0log))

                fex0ggg = fex0log2*fex0ggg
                fex0fgg = d1fex0log2*fex0gg + fex0log2*fex0fgg
                fex0ffg = d2fex0log2*fex0g + 2.d0*d1fex0log2*fex0fg +
     &            fex0log2*fex0ffg
                fex0fff = d3fex0log2*fex0 + 3.d0*d2fex0log2*fex0f +
     &            3.d0*d1fex0log2*fex0ff + fex0log2*fex0fff
                fex0gg = fex0log2*fex0gg
                fex0fg = d1fex0log2*fex0g + fex0log2*fex0fg
                fex0ff = d2fex0log2*fex0 + 2.d0*d1fex0log2*fex0f +
     &            fex0log2*fex0ff
                fex0g = fex0log2*fex0g
                fex0f = d1fex0log2*fex0 + fex0log2*fex0f
                fex0 = fex0log2*fex0
              endif
C              reverse order to preserve values on RHS.
C              ggg
              sum(10) = fex0ggg*sum(1) + 3.d0*fex0gg*sum(3) +
     &          3.d0*fex0g*sum(6) + fex0*sum(10)
C              fgg
              sum(9) = fex0fgg*sum(1) + fex0gg*sum(2) +
     &          2.d0*(fex0fg*sum(3) + fex0g*sum(5)) +
     &          fex0f*sum(6) + fex0*sum(9)
C              ffg
              sum(8) = fex0ffg*sum(1) + fex0ff*sum(3) +
     &          2.d0*(fex0fg*sum(2) + fex0f*sum(5)) +
     &          fex0g*sum(4) + fex0*sum(8)
C              fff
              sum(7) = fex0fff*sum(1) + 3.d0*fex0ff*sum(2) +
     &          3.d0*fex0f*sum(4) + fex0*sum(7)
C              gg
              sum(6) = fex0gg*sum(1) + 2.d0*fex0g*sum(3) + fex0*sum(6)
C              fg
              sum(5) = fex0fg*sum(1) + fex0g*sum(2) + fex0f*sum(3) +
     &          fex0*sum(5)
C              ff
              sum(4) = fex0ff*sum(1) + 2.d0*fex0f*sum(2) + fex0*sum(4)
C              g
              sum(3) = fex0g*sum(1) + fex0*sum(3)
C              f
              sum(2) = fex0f*sum(1) + fex0*sum(2)
              sum(1) = fex0*sum(1)
C        Transform to tcl derivative from ln g
C        derivative recalling that dlng/dlntc = 1, and dlng/dlnf = 0.5*uf.
              sum(4) = sum(4) +
     &          0.5d0*(duf*sum(3) + uf*sum(5)) +
     &          0.5d0*uf*(sum(5) + 0.5d0*uf*sum(6))
              sum(7) = sum(7) +
     &          0.5d0*(duf2*sum(3) + 2.d0*duf*sum(5) + uf*sum(8)) +
     &          0.5d0*(duf*sum(5) + uf*sum(8) +
     &          0.5d0*uf*(2.d0*duf*sum(6) + uf*sum(9))) +
     &          0.5d0*uf*(sum(8) +
     &          0.5d0*(duf*sum(6) + uf*sum(9)) +
     &          0.5d0*uf*(sum(9) + 0.5d0*uf*sum(10)))
C          do sum(5) later than sum(7) so that rhs of sum(7) is undisturbed.
              sum(5) = sum(5) + 0.5d0*uf*sum(6)
              sum(8) = sum(8) +
     &          0.5d0*(duf*sum(6) + uf*sum(9)) +
     &          0.5d0*uf*(sum(9) + 0.5d0*uf*sum(10))
              sum(9) = sum(9) + 0.5d0*uf*sum(10)
              sum(2) = sum(2) + 0.5d0*uf*sum(3)
              fexj = fexj + sum(1)
              fexjf = fexjf + sum(2)
              fexjt = fexjt + sum(3)
              fexjf2 = fexjf2 + sum(4)
              fexjft = fexjft + sum(5)
              fexjt2 = fexjt2 + sum(6)
              fexj_psi = fexj_psi + sum(2)/dpsidf
              fexj_psif = fexj_psif +
     &          (sum(4) - sum(2)*dpsidf2/dpsidf)/dpsidf
              fexj_psit = fexj_psit + sum(5)/dpsidf
              fexj_psif2 = fexj_psif2 +
     &          (sum(7) - (2.d0*sum(4)*dpsidf2 + sum(2)*dpsidf3 -
     &          2.d0*sum(2)*dpsidf2*dpsidf2/dpsidf)/dpsidf)/dpsidf
              fexj_psift = fexj_psift +
     &          (sum(8) - sum(5)*dpsidf2/dpsidf)/dpsidf
              fexj_psit2 = fexj_psit2 + sum(9)/dpsidf
            endif
          enddo
        endif
        if(iffexi.or.iffexk) then
C          Must calculate I=K^2 integral and/or K integral.
C          Calculate approximation to K integral.
C          From Paper IV, eqs. 17 and 37
C          NR component = 2 beta^{3/2} fd(1)/(1+g)^{N + 0.5}
          do ideriv = 1, maxderiv
C            F(1/2) and 3 derivatives
            fd(ideriv) = (2.d0/3.d0)*fermi_dirac_ct(psi,ideriv,0)
          enddo
C          Order is important.
          dfd(1) = fd(2)
          dfd(2) = fd(3)*dpsidf
          dfd(3) = fd(4)*dpsidf*dpsidf + fd(3)*dpsidf2
          fd(3) = fd(3)*dpsidf*dpsidf + fd(2)*dpsidf2
          fd(2) = fd(2)*dpsidf
          power = max(mgorder(5), mgorder(6)) + 0.5d0
          fex0 = 2.d0*beta*sqrt(beta)*vg**power

          fexl = fex0*fd(1)
          fexlf = fex0*fd(2) + fexl*(-0.5d0*power)*ug*uf
          fexlt = fexl*(1.5d0 - power*ug)
          fexlf2 = fex0*fd(3) +
     &      fex0*fd(2)*(-0.5d0*power)*ug*uf +
     &      fexlf*(-0.5d0*power)*ug*uf +
     &      fexl*(-0.5d0*power)*(dug*uf*0.5d0*uf + ug*duf)
          fexlft = fexlf*(1.5d0 - power*ug) +
     &      fexl*(-0.5d0*power)*(dug*uf)
          fexlt2 = fexlt*(1.5d0 - power*ug) - fexl*power*dug

          fexl_psi = fex0*dfd(1)
          fexl_psif = fex0*dfd(2) + fexl_psi*(-0.5d0*power)*ug*uf
          fexl_psit = fexl_psi*(1.5d0 - power*ug)
          fexl_psif2 = fex0*dfd(3) +
     &      fex0*dfd(2)*(-0.5d0*power)*ug*uf +
     &      fexl_psif*(-0.5d0*power)*ug*uf +
     &      fexl_psi*(-0.5d0*power)*(dug*uf*0.5d0*uf + ug*duf)
          fexl_psift = fexl_psif*(1.5d0 - power*ug) +
     &      fexl_psi*(-0.5d0*power)*(dug*uf)
          fexl_psit2 = fexl_psit*(1.5d0 - power*ug) -
     &      fexl_psi*power*dug

          fexl_psi = fexl_psi +
     &      fexl*(-0.5d0*power)*ug*uf/dpsidf
          fexl_psif = fexl_psif + (
     &      fexlf*(-0.5d0*power)*ug*uf +
     &      fexl*(-0.5d0*power)*(
     &      dug*uf*0.5d0*uf + ug*duf - ug*uf*dpsidf2/dpsidf)
     &      )/dpsidf
          fexl_psit = fexl_psit +
     &      fexlt*(-0.5d0*power)*ug*uf/dpsidf +
     &      fexl*(-0.5d0*power)*dug*uf/dpsidf
          fexl_psif2 = fexl_psif2 + (
     &      fexlf2*(-0.5d0*power)*ug*uf +
     &      fexlf*(-0.5d0*power)*(
     &      dug*uf*uf + 2.d0*ug*duf - ug*uf*dpsidf2/dpsidf) +
     &      fexl*(-0.5d0*power)*(
     &      dug2*uf*0.5d0*uf*0.5d0*uf + 1.5d0*dug*uf*duf + ug*duf2 - (
     &      dug*uf*dpsidf2*0.5d0*uf + ug*duf*dpsidf2 + ug*uf*dpsidf3 -
     &      ug*uf*dpsidf2*dpsidf2/dpsidf
     &      )/dpsidf
     &      ) - dpsidf2*(
     &      fexlf*(-0.5d0*power)*ug*uf +
     &      fexl*(-0.5d0*power)*(
     &      dug*uf*0.5d0*uf + ug*duf - ug*uf*dpsidf2/dpsidf)
     &      )/dpsidf
     &      )/dpsidf
          fexl_psift = fexl_psift + (
     &      fexlft*(-0.5d0*power)*ug*uf +
     &      fexlf*(-0.5d0*power)*dug*uf +
     &      fexlt*(-0.5d0*power)*(
     &      dug*uf*0.5d0*uf + ug*duf - ug*uf*dpsidf2/dpsidf) +
     &      fexl*(-0.5d0*power)*(
     &      dug2*uf*0.5d0*uf + dug*duf - dug*uf*dpsidf2/dpsidf)
     &      )/dpsidf
          fexl_psit2 = fexl_psit2 +
     &      fexlt2*(-0.5d0*power)*ug*uf/dpsidf +
     &      2.d0*fexlt*(-0.5d0*power)*dug*uf/dpsidf +
     &      fexl*(-0.5d0*power)*dug2*uf/dpsidf

C          fdf = f/(1+f) g^2 sqrt(g/(1+g))
          fdf = uf*g*g*sqrt(ug)
          do ikind = 5, 6
            if(mforder(ikind).ge.0.and.mgorder(ikind).ge.0) then
C              Note, for negative mod(ifexchange,10), g is zero.
C              Nevertheless, for programming
C              simplicity effsum_calc grinds through entire sum.
              call effsum_calc(f, g, mforder(ikind), mgorder(ikind),
     &          ccoeff(1,ikind), nderiv, sum)
C              from Paper IV, equation 37
              fex0 = fdf*vf**mforder(ikind)*vg**mgorder(ikind)
C              ln fex0 = ln f + 2.5 ln g 
C              -(mforder+1) ln 1+f -(mgorder+0.5) ln 1+g
              fex0f = fex0*(1.d0 - (1.d0 + mforder(ikind))*uf)
              fex0g = fex0*(2.5d0 - (0.5d0 + mgorder(ikind))*ug)
              fex0ff = fex0f*(1.d0 - (1.d0 + mforder(ikind))*uf) +
     &          fex0*(-(1.d0 + mforder(ikind)))*duf
              fex0fg = fex0g*(1.d0 - (1.d0 + mforder(ikind))*uf)
              fex0gg = fex0g*(2.5d0 - (0.5d0 + mgorder(ikind))*ug) +
     &          fex0*(-(0.5d0 + mgorder(ikind)))*dug
              fex0fff = fex0ff*(1.d0 - (1.d0 + mforder(ikind))*uf) +
     &          2.d0*(fex0f*(-(1.d0 + mforder(ikind)))*duf) +
     &          fex0*(-(1.d0 + mforder(ikind)))*duf2
              fex0ffg = fex0fg*(1.d0 - (1.d0 + mforder(ikind))*uf) +
     &          fex0g*(-(1.d0 + mforder(ikind))*duf)
              fex0fgg = fex0gg*(1.d0 - (1.d0 + mforder(ikind))*uf)
              fex0ggg = fex0gg*(2.5d0 - (0.5d0 + mgorder(ikind))*ug) +
     &          2.d0*(fex0g*(-(0.5d0 + mgorder(ikind)))*dug) +
     &          fex0*(-(0.5d0 + mgorder(ikind)))*dug2

              if(ikind.eq.6) then
C                N.B. there seems to be some significance loss here
C                that I cannot track down when looked at in isolation.
C                But it doesn't matter when combined in the do ikind loop.
                if(g.gt.1.d-2) then
                  fex0log = log(1.d0+g)
                else
                  fex0log = g*(1.d0 - g*(1.d0/2.d0 - g*(1.d0/3.d0 -
     &              g*(1.d0/4.d0 - g*(1.d0/5.d0 - g*(1.d0/6.d0 -
     &              g*(1.d0/7.d0 - g*(1.d0/8.d0))))))))
                endif
                fex0log2 = -fex0log*vg*vg
                d1fex0log2 =
     &            -ug*(vg*vg + 2.d0*fex0log2)
                d2fex0log2 =
     &            -dug*(vg*vg + 2.d0*fex0log2) -
     &            2.d0*ug*(-vg*vg*ug + d1fex0log2)
                d3fex0log2 =
     &            -dug2*(vg*vg + 2.d0*fex0log2) -
     &            4.d0*dug*(-vg*vg*ug + d1fex0log2) -
     &            2.d0*ug*(2.d0*vg*vg*ug*ug - vg*vg*dug + d2fex0log2)

                fex0ggg = d3fex0log2*fex0 + 3.d0*d2fex0log2*fex0g +
     &            3.d0*d1fex0log2*fex0gg + fex0log2*fex0ggg
                fex0fgg = d2fex0log2*fex0f + 2.d0*d1fex0log2*fex0fg +
     &            fex0log2*fex0fgg
                fex0ffg = d1fex0log2*fex0ff + fex0log2*fex0ffg
                fex0fff = fex0log2*fex0fff
                fex0gg = d2fex0log2*fex0 + 2.d0*d1fex0log2*fex0g +
     &            fex0log2*fex0gg
                fex0fg = d1fex0log2*fex0f + fex0log2*fex0fg
                fex0ff = fex0log2*fex0ff
                fex0g = d1fex0log2*fex0 + fex0log2*fex0g
                fex0f = fex0log2*fex0f
                fex0 = fex0log2*fex0
              endif
C              reverse order to preserve values on RHS.
C              ggg
              sum(10) = fex0ggg*sum(1) + 3.d0*fex0gg*sum(3) +
     &          3.d0*fex0g*sum(6) + fex0*sum(10)
C              fgg
              sum(9) = fex0fgg*sum(1) + fex0gg*sum(2) +
     &          2.d0*(fex0fg*sum(3) + fex0g*sum(5)) +
     &          fex0f*sum(6) + fex0*sum(9)
C              ffg
              sum(8) = fex0ffg*sum(1) + fex0ff*sum(3) +
     &          2.d0*(fex0fg*sum(2) + fex0f*sum(5)) +
     &          fex0g*sum(4) + fex0*sum(8)
C              fff
              sum(7) = fex0fff*sum(1) + 3.d0*fex0ff*sum(2) +
     &          3.d0*fex0f*sum(4) + fex0*sum(7)
C              gg
              sum(6) = fex0gg*sum(1) + 2.d0*fex0g*sum(3) + fex0*sum(6)
C              fg
              sum(5) = fex0fg*sum(1) + fex0g*sum(2) + fex0f*sum(3) +
     &          fex0*sum(5)
C              ff
              sum(4) = fex0ff*sum(1) + 2.d0*fex0f*sum(2) + fex0*sum(4)
C              g
              sum(3) = fex0g*sum(1) + fex0*sum(3)
C              f
              sum(2) = fex0f*sum(1) + fex0*sum(2)
              sum(1) = fex0*sum(1)
C        Transform to tcl derivative from ln g
C        derivative recalling that dlng/dlntc = 1, and dlng/dlnf = 0.5*uf.
              sum(4) = sum(4) +
     &          0.5d0*(duf*sum(3) + uf*sum(5)) +
     &          0.5d0*uf*(sum(5) + 0.5d0*uf*sum(6))
              sum(7) = sum(7) +
     &          0.5d0*(duf2*sum(3) + 2.d0*duf*sum(5) + uf*sum(8)) +
     &          0.5d0*(duf*sum(5) + uf*sum(8) +
     &          0.5d0*uf*(2.d0*duf*sum(6) + uf*sum(9))) +
     &          0.5d0*uf*(sum(8) +
     &          0.5d0*(duf*sum(6) + uf*sum(9)) +
     &          0.5d0*uf*(sum(9) + 0.5d0*uf*sum(10)))
C          do sum(5) later than sum(7) so that rhs of sum(7) is undisturbed.
              sum(5) = sum(5) + 0.5d0*uf*sum(6)
              sum(8) = sum(8) +
     &          0.5d0*(duf*sum(6) + uf*sum(9)) +
     &          0.5d0*uf*(sum(9) + 0.5d0*uf*sum(10))
              sum(9) = sum(9) + 0.5d0*uf*sum(10)
              sum(2) = sum(2) + 0.5d0*uf*sum(3)
              fexl = fexl + sum(1)
              fexlf = fexlf + sum(2)
              fexlt = fexlt + sum(3)
              fexlf2 = fexlf2 + sum(4)
              fexlft = fexlft + sum(5)
              fexlt2 = fexlt2 + sum(6)
              fexl_psi = fexl_psi + sum(2)/dpsidf
              fexl_psif = fexl_psif +
     &          (sum(4) - sum(2)*dpsidf2/dpsidf)/dpsidf
              fexl_psit = fexl_psit + sum(5)/dpsidf
              fexl_psif2 = fexl_psif2 +
     &          (sum(7) - (2.d0*sum(4)*dpsidf2 + sum(2)*dpsidf3 -
     &          2.d0*sum(2)*dpsidf2*dpsidf2/dpsidf)/dpsidf)/dpsidf
              fexl_psift = fexl_psift +
     &          (sum(8) - sum(5)*dpsidf2/dpsidf)/dpsidf
              fexl_psit2 = fexl_psit2 + sum(9)/dpsidf
            endif
          enddo
C          at this point K integral is stored in fexl.
C          Now,transform fexl depending on iffexi and iffexk values.
          if(iffexi) then
            fexi = fexl*fexl
            fexif = 2.d0*fexl*fexlf
            fexit = 2.d0*fexl*fexlt
            fexif2 = 2.d0*(fexlf*fexlf + fexl*fexlf2)
            fexift = 2.d0*(fexlt*fexlf + fexl*fexlft)
            fexit2 = 2.d0*(fexlt*fexlt + fexl*fexlt2)
            fexi_psi = 2.d0*fexl*fexl_psi
            fexi_psif = 2.d0*(fexlf*fexl_psi + fexl*fexl_psif)
            fexi_psit = 2.d0*(fexlt*fexl_psi + fexl*fexl_psit)
            fexi_psif2 = 2.d0*(fexlf2*fexl_psi + 2.d0*fexlf*fexl_psif +
     &        fexl*fexl_psif2)
            fexi_psift = 2.d0*(fexlft*fexl_psi + fexlf*fexl_psit +
     &        fexlt*fexl_psif + fexl*fexl_psift)
            fexi_psit2 = 2.d0*(fexlt2*fexl_psi + 2.d0*fexlt*fexl_psit +
     &        fexl*fexl_psit2)
          endif
          if(iffexk) then
            fex0 = -con_nr1_ratio*beta*beta
            fexk = fex0*fexl
            fexkf = fex0*fexlf
            fexkt = fex0*(2.d0*fexl + fexlt)
            fexkf2 = fex0*fexlf2
            fexkft = fex0*(2.d0*fexlf + fexlft)
            fexkt2 = fex0*(4.d0*fexl + 4.d0*fexlt + fexlt2)
            fexk_psi = fex0*fexl_psi
            fexk_psif = fex0*fexl_psif
            fexk_psit = fex0*(2.d0*fexl_psi + fexl_psit)
            fexk_psif2 = fex0*fexl_psif2
            fexk_psift = fex0*(2.d0*fexl_psif + fexl_psift)
            fexk_psit2 =
     &        fex0*(4.d0*fexl_psi + 4.d0*fexl_psit + fexl_psit2)
          endif
        endif
C        from Paper IV equation 21:
C        fex = fex0*(K^2 - J - con_nr1_ratio*beta^2 K)
C        and fex0 = 4 pi(e m^2 c^2/h^2)^2 when converted to cgs and
C        multiplied by -kT/V in accordance with fex definition.
        fex0 = con_d1_ratio*con_exchange
        fex = fex0*(fexi - fexj + fexk)
        fexf = fex0*(fexif - fexjf + fexkf)
        fext = fex0*(fexit - fexjt + fexkt)
        fexf2 = fex0*(fexif2 - fexjf2 + fexkf2)
        fexft = fex0*(fexift - fexjft + fexkft)
        fext2 = fex0*(fexit2 - fexjt2 + fexkt2)
        fex_psi = fex0*(fexi_psi - fexj_psi + fexk_psi)
        fex_psif = fex0*(fexi_psif - fexj_psif + fexk_psif)
        fex_psit = fex0*(fexi_psit - fexj_psit + fexk_psit)
        fex_psif2 = fex0*(fexi_psif2 - fexj_psif2 + fexk_psif2)
        fex_psift = fex0*(fexi_psift - fexj_psift + fexk_psift)
        fex_psit2 = fex0*(fexi_psit2 - fexj_psit2 + fexk_psit2)
      else
        stop 'exchange_gcpf: invalid ifexchange'
      endif
      end
