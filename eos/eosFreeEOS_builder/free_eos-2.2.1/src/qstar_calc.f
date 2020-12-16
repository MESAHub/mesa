C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: qstar_calc.f 627 2007-07-19 01:25:19Z airwin $
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
      subroutine qstar_calc(ifpi_fit,
     &  ifpl, ifmhd, ifneutral, ifapprox,
     &  eps_factor, nmin, nmin_max, nmax, iz, tl, x, nx,
     &  qstar, qstart, qstarx, qstart2, qstartx, qstarx2)
C       calculate internal partition function (referred to the ionized
C       state) and derivatives using Rydberg energy levels.
C       
C       The partition function sum is taken for the range
C       of minimum principal quantum numbers (given by nmin to nmin_max)
C       to infinity of 2 n^2 exp(a/n^2) times
C       planck_larkin occupation probability (ifpl true) times
C       mhd occupation probability (ifmhd true).
C       The planck-larkin occupation probability is given by
C         (1 - exp(-a)*(1+a)).
C       The *ln* of the mhd occupation probability is given by
C       -[b(1) + b(2)*n^2*(1+g(n)) + b(3)*n^4*(1+g(n))^2 +
C       b(4)*n^6*(1+g(n))^3 + b(5)*n^(15/2)*(1+h(n))].
C       g(n) = 1/(2n).
C       h(n) = (16/(3*n*K_n))^(3/2) - 1
C         = (16/(3*n))^(3/2) - 1 [for n <= 3]
C         = ((n+1)/n)^3*((n^2 + n + 1/2)/(n^2 + 7/6*n))^(3/2) - 1 [for n >=3]
C       a = c2*R*iz*2/t is calculated internally and
C       the b vector is also calculated internally from the x vector and iz.
C       input quantities:
C       ifpi_fit = 2,_use_best fit to Saumon table
C       ifpi_fit = 1,_use_best fit to opal table + extensions
C       ifpi_fit = 0,_use_best fit to original MDH table.
C       ifpl (logical) controls whether to_use_planck-larkin
C         occupation probability.
C       ifmhd (logical) controls whether to_use_mhd occupation probability.
C       ifneutral (logical) controls when b(1) through b(4) are employed
C         in the occupation probability calculation.
C       ifapprox (logical) controls whether to approximate sum (preferred
C         because much quicker and only 1.d-4 errors) or do actual sum
C         (more exact to test the approximations).
C       eps_factor is a factor used to help terminate the principal
C         quantum number sum. eps_factor = exp(-c2*min(chi)/t), where
C         min(chi) = the minimum ionization potential for all species
C         with the same a value (i.e. charge).  This definition means
C         that at the limit summand contributes approximately eps (see
C         parameter below) to ln(1+q_excited/q_ground).  With exponential
C         cutoff (MHD), this insures maximum cutoff errors are of order eps.
C         if not MHD (i.e., only Planck-Larkin occupation probability) the
C         remainder of the sum is n*summand.  With nmax of order 300 000
C         this means total relative error is of order 1.d-4 which is still
C         fine.
C       nmin to nmin_max is the range of minimum principal quantum numbers
C         required.
C       nmax is the maximum principal quantum number.  In Planck-Larkin case
C         obtain relative errors of 1.d-5 for nmin of 3 if nmax = 300 000.
C       iz is the charge number of the *ion* of the species.
C       tl is ln T.
C       x(nx), where nx = 5 are the auxiliary variables which
C         help to determine b and the mhd occupation probability.
C       n.b. index is *not* in auxiliary variable order and it is the
C         calling routine's responsibility to put x into the following order
C         (and also to reorder the derivatives appropriately):
C         using the notation of the paper on the RHS we have
C         x(i) = sigma^PI_{4-i} for i = 1,4 and
C         x(5) = sigma^PI_Z.
C       output quantities:
C       qstar(nmin_max), qstart(nmin_max), qstarx(nx,nmin_max)
C         qstart2(nmin_max), qstartx(nx,nmin_max), qstarx2(nx,nx,nmin_max) are
C         the partition function plus partial derivatives
C         wrt *tl* and x.
      implicit none
      include 'constants.h'  
      include 'pi_fit.h'  
      integer ifpi_fit,
     &  nmin, nmin_max, nmax, iz, nx, nb, nmax_reached, n, ib, jb
      parameter (nb = 5)
      double precision eps_factor, tl, a, x(nx), b(nb),
     &  qstar(nmin_max), qstart(nmin_max), qstarx(nx, nmin_max),
     &  qstart2(nmin_max), qstartx(nx, nmin_max),
     &  qstarx2(nx,nx,nmin_max)
      logical ifpl, ifmhd, ifneutral, ifapprox
      double precision rnuconst, rionconst3, const, bconst(nb)
      logical iffirst
      data iffirst/.true./
      save
      if(nx.ne.nb) stop 'qstar_calc: bad nx value'
      if(iz.lt.1) stop 'qstar_calc: bad iz value'
      if(iffirst) then
        iffirst = .false.
C         n.b. rnu = rnuconst*n**2*(1+g(n))
C        At the moment, we only calculate "neutral" radii for neutrals
C        and H2+, but because of that latter we must include the Z=iz factor
C        (see p. 117 of Condon and Shortly, "The Theory of Atomic Spectra").
        rnuconst = bohr/dble(iz)
C         n.b. rion^3 = rionconst3*iz**(-9/2)*n**7.5*(1+h(n))
        rionconst3 = 16.d0*(sqrt(3.d0/16.d0)*
     &    echarge*echarge/ergspercmm1/rydberg)**3
        const = 4.d0*pi/3.d0
        bconst(1) = const
      endif
      if(ifpi_fit.eq.0) then
        bconst(2) = 3.d0*const*
     &    exp(pi_fitx_neutral_ln_original/3.d0)*rnuconst
        bconst(3) = bconst(2)*
     &    exp(pi_fitx_neutral_ln_original/3.d0)*rnuconst
        bconst(4) = bconst(3)*
     &    exp(pi_fitx_neutral_ln_original/3.d0)*rnuconst/3.d0
        bconst(5) = const*
     &    exp(pi_fitx_ion_ln_original)*rionconst3/dble(iz)**4.5d0
      elseif(ifpi_fit.eq.1) then
        bconst(2) = 3.d0*const*
     &    exp(pi_fitx_neutral_ln/3.d0)*rnuconst
        bconst(3) = bconst(2)*
     &    exp(pi_fitx_neutral_ln/3.d0)*rnuconst
        bconst(4) = bconst(3)*
     &    exp(pi_fitx_neutral_ln/3.d0)*rnuconst/3.d0
        bconst(5) = const*
     &    exp(pi_fitx_ion_ln)*rionconst3/dble(iz)**4.5d0
      elseif(ifpi_fit.eq.2) then
        bconst(2) = 3.d0*const*
     &    exp(pi_fitx_neutral_ln_saumon/3.d0)*rnuconst
        bconst(3) = bconst(2)*
     &    exp(pi_fitx_neutral_ln_saumon/3.d0)*rnuconst
        bconst(4) = bconst(3)*
     &    exp(pi_fitx_neutral_ln_saumon/3.d0)*rnuconst/3.d0
        bconst(5) = const*
     &    exp(pi_fitx_ion_ln_saumon)*rionconst3/dble(iz)**4.5d0
      else
        stop 'qstar_calc: bad ifpi_fit value'
      endif
      a = c2*rydberg*dble(iz*iz)*exp(-tl)
C       convert rydberg(inf) to rydberg (finite) assuming:
      if(iz.eq.1) then
C       iz=1 dominated by hydrogen,
        a = a/(1.d0 + electron_mass/h_mass)
        bconst(5) = bconst(5)*(1.d0 + electron_mass/h_mass)**3
      elseif(iz.eq.2) then
C       iz =2 dominated by He,
        a = a/(1.d0 + electron_mass/(4.d0*h_mass))
        bconst(5) = bconst(5)*(1.d0 + electron_mass/(4.d0*h_mass))**3
      elseif(iz.gt.2) then
C       and iz >2 dominated by N.
        a = a/(1.d0 + electron_mass/(14.d0*h_mass))
        bconst(5) = bconst(5)*(1.d0 + electron_mass/(14.d0*h_mass))**3
      else
        stop 'qstar_calc: bad iz'
      endif
      if(ifmhd) then
        if(ifneutral) then
          do ib = 1,nb-1
            b(ib) = bconst(ib)*x(ib)
          enddo
        endif
        b(nb) = bconst(nb)*x(nb)
      endif
      if(ifapprox) then
        call qryd_approx(ifpl, ifmhd, ifneutral, eps_factor,
     &    nmin, nmin_max, nmax, a, b, nb,
     &    qstar, qstart, qstarx, qstart2, qstartx, qstarx2,
     &    nmax_reached)
      else
        call qryd_calc(ifpl, ifmhd, ifneutral, eps_factor,
     &    nmin, nmin_max, nmax, a, b, nb,
     &    qstar, qstart, qstarx, qstart2, qstartx, qstarx2,
     &    nmax_reached)
      endif
      do n = nmin, nmin_max
C         convert from b derivative to x derivative.
        if(ifmhd) then
          if(ifneutral) then
            do ib = 1,nb-1
              qstarx(ib,n) = qstarx(ib,n)*bconst(ib)
              do jb = ib, nb
                qstarx2(jb,ib,n) = qstarx2(jb,ib,n)*
     &            bconst(jb)*bconst(ib)
              enddo
C               convert from a derivatives to tl derivatives.
C               a = c2*rydberg*dble(iz)*exp(-tl)
              qstartx(ib,n) = -a*qstartx(ib,n)*bconst(ib)
            enddo
          endif
          qstarx(nb,n) = qstarx(nb,n)*bconst(nb)
          qstarx2(nb,nb,n) = qstarx2(nb,nb,n)*bconst(nb)*bconst(nb)
C           convert from a derivatives to tl derivatives.
C           a = c2*rydberg*dble(iz)*exp(-tl)
          qstartx(nb,n) = -a*qstartx(nb,n)*bconst(nb)
        endif
C         convert from a derivatives to tl derivatives.
C         a = c2*rydberg*dble(iz)*exp(-tl)
        qstart2(n) = a*(qstart(n) + a*qstart2(n))
        qstart(n) = -a*qstart(n)
      enddo
      end
