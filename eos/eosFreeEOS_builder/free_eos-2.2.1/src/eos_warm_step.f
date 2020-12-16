C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: eos_warm_step.f 822 2008-06-24 23:13:34Z airwin $
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
      subroutine eos_warm_step(
     &  allow_log, allow_log_neg, aux_old, partial_aux,
     &  ifrad, debug_dauxdv, itemp1, itemp2, delta1, delta2,
     &  ddv_aux, sumion0, sumion2,
     &  lambda, gamma_e,
     &  ifnr, nux, nuy, nuz, mion_end,
     &  n_partial_ions, n_partial_aux,
     &  sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &  partial_ions, f, eta, wf, t, n_e, pstar, 
     &  dve_exchange, dve_exchangef, dve_exchanget,
     &  full_sum0, full_sum1, full_sum2, charge, charge2,
     &  dv_pl, dv_plt,
     &  ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &  ifsame_zero_abundances,
     &  ifexcited, ifsame_under,
     &  naux, inv_aux,
     &  inv_ion, max_index,
     &  partial_elements, n_partial_elements, ion_end,
     &  ifionized, if_pteh, if_mc, ifreducedmass,
     &  ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,
     &  ifpl, ifmodified, ifh2, ifh2plus,
     &  izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &  eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &  r_ion3, r_neutral, nelements, 
     &  ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &  rhostar,
     &  ne, nef, net, sion, sionf, siont, uion,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het,
     &  h_dv, hd_dv, he_dv,
     &  pnorad, pnoradf, pnoradt, pnorad_daux,
     &  fion, fionf, fion_dv,
     &  sum0, sum0f, sum0t, sum0_dv,
     &  sum2, sum2f, sum2t, sum2_dv,
     &  extrasum, extrasumf, extrasumt, extrasum_dv,
     &  nextrasum, maxnextrasum,
     &  sumpl0, sumpl0f, sumpl0_dv,
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  rl, rf, rt, r_dv,
     &  h2, h2f, h2t, h2_dv, h2plus, h2plusf, h2plust, h2plus_dv,
     &  xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
C         eos_calc calls ionize then fixes up the result consistent with
C         formation of hydrogen molecules.
C         the purpose of the two routines is to calculate ionization fractions
C         and weighted sums over those quantities.  the weighted sums are
C         ultimately used in an iterative way to calculate chemical potentials
C         according to some free energy model.  The appropriate difference in 
C         chemical potentials are then combined to form the dv quantities, 
C         the change in the equilibrium constant of an ion relative to the
C         un-ionized reference state and the free electron
C         dv = (partial F/partial nref - partial F/partial nion - ion*partial F/partial ne)/kT.
C         These dv values (held in an array) are required input to eos_calc.
C       
C         The free energy model:
C         it is the responsibility of the calling programme to 
C         fill in the dv quantities for the various ions in a manner
C         consistent with a particular free energy model.  It should be noted that
C         at the startup of the iteration, the free-energy terms corresponding
C         to pressure ionization and the Coulomb effect are approximated as
C         functions of only ne, so the dv array elements are quite similar
C         to each other.  (The electron exchange term is always just a function
C         of ne.) Later as more complicated models are used for the Coulomb
C         effect and pressure ionization, the dv array elements can become
C         more varied.
C         input quantitites:
C         nux, nuy, nuz are number/volume of maximum possible ionization
C         electrons (divided by rho*NA) for hydrogen, helium, and metals.
C         naux is the total number of auxiliary variables.
C         ifexcited > 0 means_use_excited states (must have Planck-Larkin or
C           ifpi = 3 or 4).
C            0 < ifexcited < 10 means_use_approximation to explicit summation
C           10 < ifexcited < 20 means_use_explicit summation
C           mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules)
C             and helium.
C           mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
C           mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
C         ifsame_under = logical variable passed to ionize to control whether
C           to_use_the same underflow limits or recalculate them.
C           used for removing small discontinuities in eos when attempting
C           to obtain converged solution.
C         inv_ion(nions+2) maps ion index to contiguous ion index used in NR 
C           iteration dot products.
C         max_index is the maximum contiguous index used in NR iteration
C         partial_elements(n_partial_elements+2) index of elements treated as 
C           partially ionized consistent with ifelement.
C         ion_end(nelements) keeps track of largest ion index for each element.
C         ifionized = 0 means all elements treated as partially ionized
C         ifionized = 1 means trace metals fully ionized.
C         ifionized = 2 means all elements treated as fully ionized.
C         if_pteh controls whether pteh approximation used for Coulomb sums (1)
C           or whether_use_detailed Coulomb sum (0).
C         if_mc controls whether special approximation used for metal part
C           of Coulomb sum for the partially ionized metals.
C         if_mc = 1,_use_approximation
C         if_mc = 0, don't_use_approximation
C         ifreducedmass = 1 (use reduced mass in equilibrium constant)
C         ifreducedmass = 0 (use electron mass in equilibrium constant)
C         ifsame_abundances = 1, if call made with same abundances as previous
C           it is the responsibility of the calling programme to set
C           this flag where appropriate.  ifsame_abundances = 0 means
C           slower execution of ionize, but is completely safe against
C           abundance changes.
C         ifmtrace just to keep track of changes in ifelement.
C         iatomic_number(nelements), the atomic number of each element.
C         ifpi = 0, 1, 2, 3, 4 if no, pteh, geff, mhd, or Saumon-style
C           pressure ionization
C         ifpl = 0, 1 if there is no or if there is planck-larkin
C         ifmodified > 0 or not affects ifpi_fit inside excitation_sum.
C         ifh2 non-zero implies h2 is included.
C         ifh2plus non-zero implies h2plus is included.
C         izlo is lowest core charge (1 for H, He, 2 for He+, etc.)
C           for all species.
C         izhi is highest core charge for all species.
C         bmin(izhi) is the minimum bion for all species with a given core charge.
C         nmin(izhi) is the minimum excited principal quantum number for
C           all species with a given core charge.
C         nmin_max(izhi) is the largest minimum excited principal
C           quantum number for all species with a given core charge.
C         nmin_species(nions+2) is the minimum excited principal quantum number
C           organized by species.
C         nmax is maximum principal quantum number included in sum
C           (to be compatible with opal which used nmax = 4 rather than infinity).
C           if(nmax > 300 000) then treated as infinity in qryd_approx.
C           otherwise nmax is meant to be used with qryd_calc only (i.e.,
C           case for mhd approximations not programmed).
C         eps(nelements) is the ratio of mass abundance/atomic mass
C         tl = log(t)
C         tc2 = c2/t
C         bi(nions_in+2) ionization potentials (cm**-1) in ion order:
C           h, he, he+, c, c+, etc.  H2 and H2+ are last two.
C         n.b. the code logic now demands there are iatomic_number ionized
C           states + 1 neutral state for each element.  Thus, the
C           ionization potentials must include all ions, not just
C           first ions (as for trace metals in old code).
C         h2diss, h2 dissociation energy
C         plop(nions+2), plopt, plopt2 planck-larkin occupation probabilities
C           + derivatives.
C         r_ion3(nions+2), *cube of the*
C           effective radii (in ion order but going from neutral
C           to next to bare ion for each species) of MDH interaction between
C           all but bare nucleii species and ionized species.  last 2 are
C           H2 and H2+ (used externally)
C         r_neutral(nelements+2) effective radii for MDH neutral-neutral 
C           interactions.  Last two (used externally) are for for H2 and H2+ (the only ionic
C           species in the MHD model with a non-zero hard-sphere radius).
C         ifelement(nelements) = 1 if element treated as partially ionized
C         (depending on abundance, and trace metal treatment), 0 if
C         element treated as fully ionized or has no abundance.
C         dvzero(nelements) is a zero point shift that is added to dv in ionize,
C           but which is zero for hydrogen.
C         dv(nions) change in equilibrium constant (explained above)
C         dvf(nions) = d dv/d ln f
C         dvt(nions) = d dv/d ln t
C         nion(nions), charge on ion in ion order (must be same order as bi)
C           e.g., for H+, He+, He++, etc.
C         re, ref, ret: fermi-dirac integral and derivatives
C         output quantities:
C         ne, nef, net is nu(e) = n(positive ions)/(Navogadro*rho) and its 
C           derivatives calculated from all elements.
C         sion, sionf, siont is the ideal entropy/R per unit mass and its 
C           derivatives with respect to lnf and ln t.
C         n.b. if ionized, then the contribution from all elements is
C         calculated in ionize.f.  If partially ionized, then ionize separates
C         the hydrogen from the rest of the components.  If molecular, the
C         hydrogen component from ionize is completely ignored, and
C         recalculated in this routine.
C         n.b.  sion returns a component of entropy/R per unit mass,
C         sion = the sum over all non-electron species of
C         -nu_i*[-5/2 - 3/2 ln T + ln(alpha) - ln Na 
C           + ln (n_i/[A_i^{3/2} Q_i]) - dln Q_i/d ln T],
C         where nu_i = n_i/(Na rho), Na is the Avogadro number,
C         A_i is the atomic weight,
C         alpha =  (2 pi k/[Na h^2])^(-3/2) Na
C         (see documentation of alpha^2 in constants.h), and
C         Q_i is the internal ideal partition function
C         of the non-Rydberg states.  (Currently, we calculate this
C         partition function by the statistical weights of the ground states
C         of helium and the combined lower states (roughly
C         approximated) of each of the metals.
C         helium is subsequently corrected for detailed excitation
C         of the non-Rydberg states, and this crude approximation for the metals
C         is not currently corrected.  Thus, in all *current* monatomic
C         cases ln Q_i is a constant, and dln Q_i/d ln T is zero, but this
C         will change for the metals eventually.
C         currently, hydrogen is treated exactly for all cases (full ionization,
C         partial ionization, partial molecular formation).
C         From the equilibrium constant approach and the monatomic species
C         treated in this subroutine (molecules treated outside) we have
C         -nu_i ln (n_i/[A_i^{3/2} Q_i]) =
C         -nu_i * [ln (n_neutral/[A_neutral^{3/2} Q_neutral) +
C         (- chi_i/kT + dv_i)]
C         if we sum this term over all species
C         of an element without molecules we obtain
C         s_element = - eps * ln (eps*n_neutral/sum(n))
C         - eps ln (alpha/(A_neutral^{3/2} Q_neutral)) -
C         - sum over all species of the element of nu_i*(-chi/kT + dv(i))
C         where we have ignored the term
C         eps * (-5/2 - 3/2 ln T + ln rho)
C         (taking into account the first 4 terms above).
C         If hydrogen molecules are included,
C         the result is the same except for the addition of the
C         d ln Q_i/d ln T term (which will also appear for the metals eventually)
C         and the 5/2 factor is multiplied by sum over all species which is
C         corrected to eps by subtracting 5/2 (nu(H2) + nu(H2+)) from sion below.
C         note a final correction of the s zero point occurs
C         in awieos_detailed.f which puts back the ignored terms for both
C         the case of molecules and no molecules.
C         uion is the ideal internal energy (cm^-1 per unit mass divided by
C           avogadro).
C         h, hf, ht is n(H+) and its derivatives with respect to ln f and ln t.
C         hd, hdf, hdt is n(He+) and its derivatives with respect to 
C           ln f and ln t.
C         he, hef, het is n(He++) and its derivatives with respect to 
C           ln f and ln t.
C         the following are returned only if .not.( ifpteh.eq.1.or.if_mc.eq.1)
C         sum0 and f,t,dv derivatives, sum over positive charge number densities
C           with uniform weights.
C         sum2 and f,t,dv derivatives, sum over positive charge number densities
C           weighted by charge^2.
C         the following are returned only if ifpi = 3 or 4.
C         extrasum(maxnextrasum) weighted sums over n(i).
C           for iextrasum = 1,nextrasum-2, sum is only over
C           neutral species + H2+ (the only species with
C           non-zero radii according to the MHD model) and
C           weight is r_neutral^{iextrasum-1}.
C           for iextrasum = nextrasum-1, sum is over all ionized species including
C           bare nucleii, but excluding free electrons, weight is Z^1.5.
C           for iextrasum = nextrasum, sum is over all species excluding bare nucleii
C           and free electrons, the weight is rion^3.
C         extrasumf(maxnextrasum) = partial of extrasum/partial ln f
C         extrasumt(maxnextrasum) = partial of extrasum/partial ln t
C         the following are returned only if ifpl = 1
C         sumpl1 and sumpl2 = weighted sums over non-H nu(i) = n(i)/(rho/H).
C           for sumpl1 sum is over Planck-Larkin occupation
C           probabilities + d ln w/d ln T
C           for sumpl2, sum is over Planck-Larkin d ln w/d ln T
C         sumpl1f and sumpl1t = derivatives of first sum wrt lnf and lnt.
C         rl, rf, rt are the ln mass density and derivatives
C         h2, h2f, h2t are n(H2) and derivatives.
C         h2plus, h2plusf, h2plust are n(H2+) and derivatives.
C         h2plus_dv is the h2plus derivatives wrt dv.  n.b. this vector only
C           returned if molecular hydrogen is calculated.  The calling
C           routine (awieos_detailed) uses it only if 
C           if_mc.eq.1.and.ifh2plus.gt.0.
C         the following are returned only if ifexcited.gt.0 and
C           ifpi.eq.3.or4.
C         xextrasum(4) is the *negative* sum nuvar/(1 + qratio)*
C         partial qratio/partial extrasum(k), for k = 1, 2, 3, and nextrasum-1.
C         xextrasumf(4), xextrasumt(4), xextrasum_dv(nions+2,4) are the
C         fl and tl derivatives (ifnr.eq.0) or dv derivatives (ifnr.eq.1)
      implicit none
      include 'constants.h'
      integer nions, nelements, nion(nions+2)
      integer ion_end(nelements+2)
      integer ifreducedmass, ifmodified, ifexcited,
     &  nmax
      double precision bi(nions+2), h2diss, t, tl,
     &  eps(nelements), eta
      integer nstar
      parameter (nstar = 9)
      double precision rhostar(nstar), pstar(nstar),
     &  ne, nef, net
      integer maxfjs_aux
      parameter (maxfjs_aux = 4)
      double precision dvzero(nelements), 
     &  dv(nions+2), dvf(nions+2), dvt(nions+2),
     &  dve, dvef, dvet, dve0, dve2, 
     &  dve_pi, dve_pif, dve_pit, dve_pi_aux(maxfjs_aux),
     &  dve_coulomb, dve_coulombf, dve_coulombt,
     &  dve_exchange, dve_exchangef, dve_exchanget,
     &  dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22,
     &  plop(nions+2), plopt(nions+2), plopt2(nions+2),
     &  dv_pl(nions+2), dv_plt(nions+2)
      integer maxcorecharge
C      maximum core charge for non-bare ion
      parameter(maxcorecharge = 28)
      integer nmin(maxcorecharge), nmin_max(maxcorecharge),
     &  izlo, izhi
      integer nmin_species(nions+2)
      integer iatomic_number(nelements)
      integer ifnr,
     &  ifh2, ifh2plus, ifmtrace,
     &  if_pteh, ifcoulomb_mod, if_mc, if_dc,
     &  ifionized, ion, index_ion, max_index,
     &  index, maxnextrasum,
     &  nextrasum, ifpi_local, ifpl, nxextrasum,
     &  ifsame_abundances,
     &  ifsame_zero_abundances,
     &  ion0
      parameter (nxextrasum = 4)
      double precision lambda, gamma_e,
     &  full_sum0, full_sum1, full_sum2,
     &  sum0, sum0ne, sum0f, sum0t,
     &  sum0_aux(maxfjs_aux+1), sum2_aux(maxfjs_aux+1),
     &  sum2, sum2ne, sum2f, sum2t,
     &  r_ion3(nions+2), r_neutral(nelements+2),
     &  charge(nions+2), charge2(nions+2),
     &  extrasum(maxnextrasum), extrasumf(maxnextrasum), 
     &  extrasumt(maxnextrasum),
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  tc2
      double precision f, wf, n_e,
     &  sion, sionf, siont, uion,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het, rt, rf,
     &  h2, h2f, h2t, h2_dv(nions+2),
     &  h2plus, h2plusf, h2plust, h2plus_dv(nions+2)
      integer ifdv(nions+2), ifelement(nelements),
     &  partial_ions(nions+2), inv_ion(nions+2),
     &  n_partial_elements, partial_elements(nelements+2)
      double precision rl
C      iextraoff is the number of non-extrasum auxiliary variables
      integer iextraoff, naux
      parameter (iextraoff = 7)
      integer 
     &  inv_aux(naux), ielement
      double precision 
     &  h_dv(nions+2), hd_dv(nions+2), he_dv(nions+2),
     &  r_dv(nions+2), sum0_dv(nions+2), sum2_dv(nions+2),
     &  extrasum_dv(nions+2, maxnextrasum),
     &  xextrasum(nxextrasum), xextrasumf(nxextrasum),
     &  xextrasumt(nxextrasum), xextrasum_dv(nions+2,nxextrasum),
     &  bmin(maxcorecharge)
      logical 
     &  ifsame_under, ifdvzero, ifnuform
      logical allow_log(naux), allow_log_neg(naux)
      double precision aux_old(naux)
      integer partial_aux(naux)
      integer index_aux, n_partial_aux
      double precision ddv_aux(nions+2,naux), rho, rho_new,
     &  dvh, dvhf, dvht, dvh_aux(maxfjs_aux),
     &  dvhe1, dvhe1f, dvhe1t, dvhe1_aux(maxfjs_aux),
     &  dvhe2, dvhe2f, dvhe2t, dvhe2_aux(maxfjs_aux),
     &  sumion0, sumion2,
     &  dsum0, dsum2,
     &  sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &  nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt,
     &  nux, nuy, nuz
      integer
     &  inv_ion_index, ifjs_aux, iaux, jndex_aux,
     &  mion_end, n_partial_ions
C      variables associated with pressure and free energy calculation:
      integer nions_local, naux_local, maxnextrasum_local
      parameter(nions_local = 295)
      parameter(naux_local = 21)
C      Note, this value must be greater than or equal to 4 as well to store
C      fjs_pi_free and corresponding pressure quantities.
      parameter(maxnextrasum_local = 9)
      double precision
     &  pe_cgs,
     &  p0, pion, pionf, piont, pion_dv(nions_local+2),
     &  pex, pexf, pext,
     &  ppi, ppif, ppit, ppir, ppi_dv(nions_local+2),
     &  ppi2_daux(maxnextrasum_local), ppi_daux(naux_local),
     &  pexcited, pexcitedf, pexcitedt, pexcited_dv(nions_local+2),
     &  dpcoulomb, dpcoulombf, dpcoulombt, dpcoulomb0, dpcoulomb2,
     &  pnorad_dv(nions_local+2),
     &  pnorad, pnoradf, pnoradt, pnorad_daux(n_partial_aux+2),
!     &  free_rad, free_e, free_ef,
     &  free_pi, free_pif,! free_pit,
!     &  free_pir, free_pi_dv(nions_local+2),
     &  free_pi2_daux(maxnextrasum_local),
!     &  free_pi_daux(naux_local),
     &  free_excited, free_excitedf, free_excited_dv(nions_local+2),
!     &  free_coulomb, free_coulombf, free_coulomb0, free_coulomb2,
!     &  free_ex, free_exf,
     &  ni, nif, ni_dv(nions_local+2),
!     &  free_pl, free_plf, free_pl_dv(nions_local+2),
     &  fion, fionf, fion_dv(nions+2),
!     &  fion_daux(n_partial_aux+2),
     &  sumpl0, sumpl0f, sumpl0_dv(nions+2)
!     &  free, freef, free_dv(nions_local+2),
!     &  free_daux(n_partial_aux+2)
      logical ifnr03, ifnr13
C      for debugging
      logical debug_dauxdv
      integer itemp1, itemp2
      double precision delta1, delta2
      integer ifrad
      save
C      sanity checks.
      if(nextrasum.gt.maxnextrasum_local)
     &  stop 'eos_warm_step: maxnextrasum_local too small'
      if(.not.(ifnr.eq.0.or.ifnr.eq.1.or.ifnr.eq.3))
     &  stop 'eos_warm_step: ifnr must be 0, 1, or 3'
      if(ifnr.eq.0) stop 'eos_warm_step: ifnr = 0 is disabled'
      if(nions.ne.nions_local)
     &  stop 'eos_warm_step: nions must be equal to nions_local'
      if(naux.ne.naux_local)
     &  stop 'eos_warm_step: naux must be equal to naux_local'
      if(n_partial_aux.gt.naux_local)
     &  stop 'eos_warm_step: n_partial_aux too large'
      ifnr03 = ifnr.eq.0.or.ifnr.eq.3
      ifnr13 = ifnr.eq.1.or.ifnr.eq.3
      rho = exp(rl)
      do index = 1, n_partial_elements
        ielement = partial_elements(index)
        dvzero(ielement) = 0.d0
      enddo
      do index_ion = 1,max_index
        ion = partial_ions(index_ion)
        dv(ion) = 0.d0
        if(ifnr03) then
          dvf(ion) = 0.d0
          dvt(ion) = 0.d0
        endif
      enddo
      if(ifnr13) then
        do index_aux = 1,n_partial_aux
          do index_ion = 1,max_index
            ddv_aux(index_ion, index_aux) = 0.d0
          enddo
        enddo
      endif
      if(ifpi_local.eq.0) then
        dve_pi = 0.d0
        if(ifnr03) then
          dve_pif = 0.d0
          dve_pit = 0.d0
        endif
      elseif(ifpi_local.eq.1) then
        call pteh_pi(ifmodified,
     &    rhostar(1), rhostar(2), rhostar(3), t,
     &    dve_pi, dve_pif, dve_pit)
      elseif(ifpi_local.eq.2) then
C       _use_FJS pressure ionization
        nx = nux*rho*avogadro
        ny = nuy*rho*avogadro
        nz = nuz*rho*avogadro
        if(ifnr.eq.0) then
          nxf = nx*rf
          nxt = nx*rt
          nyf = ny*rf
          nyt = ny*rt
          nzf = nz*rf
          nzt = nz*rt
        elseif(ifnr.eq.3) then
C          For ifnr.eq.3 calculate all f and t derivatives assuming input
C          auxiliary variables (e.g., rho for ifpi_local.eq.2) are fixed.
          nxf = 0.d0
          nxt = 0.d0
          nyf = 0.d0
          nyt = 0.d0
          nzf = 0.d0
          nzt = 0.d0
        endif
        call fjs_pi(ifnr, maxfjs_aux, ifmodified,
     &    t, nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt,
     &    h, hf, ht, hd, hdf, hdt, he, hef, het,
     &    n_e, n_e*rhostar(2), n_e*rhostar(3),
     &    dvh, dvhf, dvht, dvh_aux,
     &    dvhe1, dvhe1f, dvhe1t, dvhe1_aux,
     &    dvhe2, dvhe2f, dvhe2t, dvhe2_aux,
     &    dve_pi, dve_pif, dve_pit, dve_pi_aux)
C        this call must occur before call to eos_calc (which
C        calculates new values of rho (used to calculate nx, ny, and nz),
C        h, hd, and he.
C        n.b. ppi and all its derivatives are equivalent to
C        free_pi and all its derivatives for the fjs form of
C        pressure ionization.
        call fjs_pi_free(nx, ny, nz, h, hd, he,
     &    t, n_e, n_e*rhostar(2), n_e*rhostar(3),
     &    ppi, ppif, ppit, ppi2_daux)
C        update dv keeping in mind that calculated dvh etc. have
C        dve_pi included.  the dve_pi quantity
C        must be subtracted so that further free_eos_detailed logic
C        which adds dve_pi works properly.
        if(eps(1).gt.0.d0) then
          dv(1) = dv(1) + dvh - dve_pi
          if(ifnr03) then
            dvf(1) = dvf(1) + dvhf - dve_pif
            dvt(1) = dvt(1) + dvht - dve_pit
          endif
          if(ifnr13) then
            inv_ion_index = inv_ion(1)
            do ifjs_aux = 1,maxfjs_aux
              iaux = inv_aux(ifjs_aux)
              if(iaux.gt.0) then
                ddv_aux(inv_ion_index,iaux) =
     &            ddv_aux(inv_ion_index,iaux) +
     &            dvh_aux(ifjs_aux) - dve_pi_aux(ifjs_aux)
              endif
            enddo
          endif
        endif
        if(eps(2).gt.0.d0) then
          if(ifdvzero) then
            dv(2) = dv(2) + dvhe1 - dve_pi
            dv(3) = dv(3) + dvhe1 + dvhe2 - 2.d0*dve_pi
          else
            dvzero(2) = dvzero(2) + dvhe1 + dvhe2 - 2.d0*dve_pi
            dv(2) = dv(2) - dvhe2 + dve_pi
          endif
          if(ifnr03) then
            dvf(2) = dvf(2) + dvhe1f - dve_pif
            dvt(2) = dvt(2) + dvhe1t - dve_pit
            dvf(3) = dvf(3) + dvhe1f + dvhe2f - 2.d0*dve_pif
            dvt(3) = dvt(3) + dvhe1t + dvhe2t - 2.d0*dve_pit
          endif
          if(ifnr13) then
            inv_ion_index = inv_ion(2)
            do ifjs_aux = 1,maxfjs_aux
              iaux = inv_aux(ifjs_aux)
              if(iaux.gt.0) then
                ddv_aux(inv_ion_index,iaux) =
     &            ddv_aux(inv_ion_index,iaux) +
     &            dvhe1_aux(ifjs_aux) - dve_pi_aux(ifjs_aux)
              endif
            enddo
            inv_ion_index = inv_ion(3)
            do ifjs_aux = 1,maxfjs_aux
              iaux = inv_aux(ifjs_aux)
              if(iaux.gt.0) then
                ddv_aux(inv_ion_index,iaux) =
     &            ddv_aux(inv_ion_index,iaux) +
     &            dvhe1_aux(ifjs_aux) + dvhe2_aux(ifjs_aux) -
     &            2.d0*dve_pi_aux(ifjs_aux)
              endif
            enddo
          endif
        endif
C        no change in H2 equilibrium constant because fjs_pi is independent
C        n(H) or n(H2).  Change in H2+ equilibrium constant
C        via dve_pi which is handled later as part of dve.
      elseif(ifpi_local.eq.3.or.ifpi_local.eq.4) then
C       _use_mdh-like pressure ionization and dissociation
        call mdh_pi(
     &    ifdvzero, ifnr, inv_aux, iextraoff, naux, inv_ion,
     &    partial_elements, n_partial_elements,
     &    ion_end, mion_end,
     &    ifmodified,
     &    r_ion3, nion, nions, r_neutral, nelements,
     &    extrasum, extrasumf, extrasumt, nextrasum,
     &    ifdv, dvzero, dv, dvf, dvt, ddv_aux)
C        this call must occur before call to eos_calc (which
C        updates extrasum to new values.)
        call mdh_pi_pressure_free(t, extrasum, nextrasum,
     &    ppi, ppif, ppit, ppi2_daux,
     &    free_pi, free_pif, free_pi2_daux)
C        MDH pressure ionization has no N_e dependence
        dve_pi = 0.d0
        if(ifnr03) then
          dve_pif = 0.d0
          dve_pit = 0.d0
        endif
      endif
C      calculate change in equilibrium constant due to electron degeneracy
      dve = dve_pi - eta
      if(ifnr03) then
        dvef = dve_pif - wf
        dvet = dve_pit
      endif
C      n.b. for master_coulomb to work properly, sum0 and sum2
C      are considered to be functions of ne, f, t.
C      for pteh sum approximation, the sums depend on ne, and
C      there is no explicit f, t dependence.
      if(if_pteh.eq.1) then
C        PTEH (full ionization) approximation to sum0, sum2
        sum0ne = full_sum0/full_sum1
        sum0 = sum0ne*n_e
        sum2ne = full_sum2/full_sum1
        sum2 = sum2ne*n_e
        sum0f = 0.d0
        sum0t = 0.d0
        sum2f = 0.d0
        sum2t = 0.d0
      else
C        derivative wrt ne is zero (this must be
C        asserted otherwise master_coulomb doesn't work
C        correctly)
        sum0ne = 0.d0
        sum2ne = 0.d0
      endif
      if(if_mc.eq.1) then
C        in just this case
C        ionize and eos_calc have only calculated fully
C        ionized part/(rho*avogadro).
        sumion0 = sum0_mc*rho*avogadro
        sumion2 = sum2_mc*rho*avogadro
        sum0 = sumion0 + h2plus + h*(1.d0+hcon_mc) +
     &    hd + he
        sum2 = sumion2 + h2plus + h*(1.d0+hcon_mc) +
     &    hd + he*(4.d0+hecon_mc)
        if(ifnr.eq.0) then
C          N.B. this code depends on auxiliary variables that are only
C          defined for ifpi_local = 2 case.  So it only works because
C          if_mc.eq.1 is correlated with that case.  Check this.
          if(ifpi_local.ne.2)
     &      stop 'eos_warm_step: logical screwup wrt ifpi_local'
          sum0f = sumion0*rf + h2plusf + hf*(1.d0+hcon_mc) +
     &      hdf + hef
          sum0t = sumion0*rt + h2plust + ht*(1.d0+hcon_mc) +
     &      hdt + het
          sum2f = sumion2*rf + h2plusf + hf*(1.d0+hcon_mc) +
     &      hdf + hef*(4.d0+hecon_mc)
          sum2t = sumion2*rt + h2plust + ht*(1.d0+hcon_mc) +
     &      hdt + het*(4.d0+hecon_mc)
        elseif(ifnr.eq.3) then
C          For ifnr.eq.3 calculate all f and t derivatives assuming input
C          auxiliary variables are fixed.
          sum0f = 0.d0
          sum0t = 0.d0
          sum2f = 0.d0
          sum2t = 0.d0
        endif
        if(ifnr13) then
          sum0_aux(1) = (1.d0+hcon_mc)
          sum2_aux(1) = (1.d0+hcon_mc)
          sum0_aux(2) = 1.d0
          sum2_aux(2) = 1.d0
          sum0_aux(3) = 1.d0
          sum2_aux(3) = (4.d0+hecon_mc)
          sum0_aux(4) = sumion0
          sum2_aux(4) = sumion2
          sum0_aux(5) = 1.d0
          sum2_aux(5) = 1.d0
        endif
      endif
      call master_coulomb(ifnr,
     &  rhostar, f,
     &  sum0, sum0ne, sum0f, sum0t, sum2, sum2ne, sum2f, sum2t,
     &  n_e, t, cpe*pstar(1), pstar, lambda, gamma_e,
     &  ifcoulomb_mod, if_dc, if_pteh,
     &  dve_coulomb, dve_coulombf, dve_coulombt, dve0, dve2,
     &  dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22)
C      N.B. this routine must be called _before_ sum0 and sum2 are changed
C      in value.
      call master_coulomb_pressure(
     &  ifcoulomb_mod, if_dc, if_pteh, sum0, sum2,
     &  t, n_e, rhostar, pstar,
     &  dpcoulomb, dpcoulombf, dpcoulombt, dpcoulomb0, dpcoulomb2)
      dve = dve + dve_coulomb
      if(ifnr03) then
        dvef = dvef + dve_coulombf
        dvet = dvet + dve_coulombt
      endif
C      add in pre-calculated exchange effects
      dve = dve + dve_exchange
      if(ifnr03) then
        dvef = dvef + dve_exchangef
        dvet = dvet + dve_exchanget
      endif
      if(if_mc.ne.1) then
        index_aux = inv_aux(iextraoff-1)
        jndex_aux = inv_aux(iextraoff)
        do index = 1, n_partial_elements
          ielement = partial_elements(index)
          if(ifdvzero.or.ielement.eq.1) then
          else
            dvzero(ielement) = dvzero(ielement) +
     &        dv0 + charge(ion_end(ielement))*
     &        (dve + charge(ion_end(ielement))*dv2)
          endif
          if(ielement.gt.1) then
            ion0 = ion_end(ielement-1)
          else
            ion0 = 0
          endif
          do ion = ion0+1,ion_end(ielement)
            if(ifdvzero.or.ielement.eq.1) then
              dv(ion) = dv(ion) +
     &          dv0 + charge(ion)*(dve + charge(ion)*dv2)
            else
              dv(ion) = dv(ion) +
     &          dble(nion(ion)-nion(ion_end(ielement)))*
     &          (dve + dble(nion(ion)+nion(ion_end(ielement)))*dv2)
            endif
            if(ifnr03) then
              dvf(ion) = dvf(ion) +
     &          dv0f + charge(ion)*(dvef + charge(ion)*dv2f)
              dvt(ion) = dvt(ion) +
     &          dv0t + charge(ion)*(dvet + charge(ion)*dv2t)
            endif
            if(ifnr13.and.index_aux.ne.0) then
              index_ion = inv_ion(ion)
              ddv_aux(index_ion,index_aux) =
     &          ddv_aux(index_ion,index_aux) +
     &          dv00 + charge(ion)*(dve0 + charge(ion)*dv02)
              ddv_aux(index_ion,jndex_aux) =
     &          ddv_aux(index_ion,jndex_aux) +
     &          dv02 + charge(ion)*(dve2 + charge(ion)*dv22)
            endif
          enddo
        enddo
C        H2+
        if(ifdv(nions+2).eq.1) then
          dv(nions+2) = dv(nions+2) + dv0 + dve + dv2
          if(ifnr03) then
            dvf(nions+2) = dvf(nions+2) + dv0f + dvef + dv2f
            dvt(nions+2) = dvt(nions+2) + dv0t + dvet + dv2t
          endif
          if(ifnr13.and.index_aux.ne.0) then
C            max_index is correct pointer when ifdv(nions+2).eq.1
            ddv_aux(max_index,index_aux) =
     &        ddv_aux(max_index,index_aux) +
     &        dv00 + dve0 + dv02
            ddv_aux(max_index,jndex_aux) =
     &        ddv_aux(max_index,jndex_aux) +
     &        dv02 + dve2 + dv22
          endif
        endif
      else
C        if_mc *is* 1
        do index = 1, n_partial_elements
          ielement = partial_elements(index)
          if(ifdvzero.or.ielement.eq.1) then
          elseif(ielement.eq.2) then
            dvzero(ielement) = dvzero(ielement) +
     &        dv0 + charge(ion_end(ielement))*
     &        (dve + charge(ion_end(ielement))*dv2)
          else
            dvzero(ielement) = dvzero(ielement) +
     &        charge(ion_end(ielement))*dve
          endif
          if(ielement.gt.1) then
            ion0 = ion_end(ielement-1)
          else
            ion0 = 0
          endif
          do ion = ion0+1,ion_end(ielement)
            if(ifnr13) then
              if(ion.ge.4) then
                dsum0 = charge(ion)*dve0
                dsum2 = charge(ion)*dve2
              elseif(ion.eq.1) then
                dsum0 = charge(ion)*dve0 + (1.d0 + hcon_mc)*
     &            (dv00 + charge2(ion)*dv02)
                dsum2 = charge(ion)*dve2 + (1.d0 + hcon_mc)*
     &            (dv02 + charge2(ion)*dv22)
              elseif(ion.eq.2) then
                dsum0 = 
     &            dv00 + charge(ion)*(dve0 + charge(ion)*dv02)
                dsum2 =
     &            dv02 + charge(ion)*(dve2 + charge(ion)*dv22)
              elseif(ion.eq.3) then
                dsum0 = 
     &            charge(ion)*dve0 + dv00 +
     &            (hecon_mc + charge2(ion))*dv02
                dsum2 = 
     &            charge(ion)*dve2 + dv02 +
     &            (hecon_mc + charge2(ion))*dv22
              endif
              index_ion = inv_ion(ion)
              do ifjs_aux = 1, maxfjs_aux+1
                index_aux = inv_aux(ifjs_aux)
                if(index_aux.gt.0)
     &            ddv_aux(index_ion,index_aux) =
     &            ddv_aux(index_ion,index_aux) +
     &            dsum0*sum0_aux(ifjs_aux) + dsum2*sum2_aux(ifjs_aux)
              enddo
            endif
            if(ion.ge.4) then
              if(ifdvzero) then
                dv(ion) = dv(ion) + charge(ion)*dve
              else
                dv(ion) = dv(ion) +
     &            dble(nion(ion)-nion(ion_end(ielement)))*dve
              endif
              if(ifnr03) then
                dvf(ion) = dvf(ion) + charge(ion)*dvef
                dvt(ion) = dvt(ion) + charge(ion)*dvet
              endif
            elseif(ion.eq.1) then
              dv(ion) = dv(ion) +
     &          charge(ion)*dve + (1.d0 + hcon_mc)*
     &          (dv0 + charge2(ion)*dv2)
              if(ifnr03) then
                dvf(ion) = dvf(ion) +
     &            charge(ion)*dvef + (1.d0 + hcon_mc)*
     &            (dv0f + charge2(ion)*dv2f)
                dvt(ion) = dvt(ion) +
     &            charge(ion)*dvet + (1.d0 + hcon_mc)*
     &            (dv0t + charge2(ion)*dv2t)
              endif
            elseif(ion.eq.2) then
              if(ifdvzero) then
                dv(ion) = dv(ion) +
     &            dv0 + charge(ion)*(dve + charge(ion)*dv2)
              else
                dv(ion) = dv(ion) +
     &            dble(nion(ion)-nion(ion_end(ielement)))*
     &            (dve + dble(nion(ion)+nion(ion_end(ielement)))*dv2)
              endif
              if(ifnr03) then
                dvf(ion) = dvf(ion) +
     &            dv0f + charge(ion)*(dvef + charge(ion)*dv2f)
                dvt(ion) = dvt(ion) +
     &            dv0t + charge(ion)*(dvet + charge(ion)*dv2t)
              endif
            elseif(ion.eq.3) then
              if(ifdvzero) then
                dv(ion) = dv(ion) + charge(ion)*dve + dv0 +
     &            (hecon_mc + charge2(ion))*dv2
              else
                dv(ion) = dv(ion) + hecon_mc*dv2
              endif
              if(ifnr03) then
                dvf(ion) = dvf(ion) + charge(ion)*dvef + dv0f +
     &            (hecon_mc + charge2(ion))*dv2f
                dvt(ion) = dvt(ion) + charge(ion)*dvet + dv0t +
     &            (hecon_mc + charge2(ion))*dv2t
              endif
            endif
          enddo
        enddo
C        H2+
        if(ifdv(nions+2).eq.1) then
          dv(nions+2) = dv(nions+2) + dv0 + dve + dv2
          if(ifnr03) then
            dvf(nions+2) = dvf(nions+2) + dv0f + dvef + dv2f
            dvt(nions+2) = dvt(nions+2) + dv0t + dvet + dv2t
          endif
          if(ifnr13) then
            dsum0 = dv00 + dve0 + dv02
            dsum2 = dv02 + dve2 + dv22
            do ifjs_aux = 1, maxfjs_aux+1
              index_aux = inv_aux(ifjs_aux)
C              max_index is correct pointer when ifdv(nions+2).eq.1
              if(index_aux.gt.0)
     &          ddv_aux(max_index,index_aux) =
     &          ddv_aux(max_index,index_aux) +
     &          dsum0*sum0_aux(ifjs_aux) + dsum2*sum2_aux(ifjs_aux)
            enddo
          endif
        endif
      endif
      if(ifpi_local.eq.2.and.ifnr13) then
C        in all monatomic cases added term of charge(ion)*dve to dv
C        for ifpi_local == 2, there is a pressure ionization
C        component to dve.
        do index_ion = 1,n_partial_ions
          ion = partial_ions(index_ion)
          do ifjs_aux = 1,maxfjs_aux
            index_aux = inv_aux(ifjs_aux)
            if(index_aux.gt.0) then
              ddv_aux(index_ion,index_aux) =
     &          ddv_aux(index_ion,index_aux) +
     &          charge(ion)*dve_pi_aux(ifjs_aux)
            endif
          enddo
        enddo
C        for H2+ added term of dve to dv(nions+2)
        if(ifdv(nions+2).eq.1) then
          do ifjs_aux = 1,maxfjs_aux
            index_aux = inv_aux(ifjs_aux)
            if(index_aux.gt.0) then
C              max_index is correct pointer when ifdv(nions+2).eq.1
              ddv_aux(max_index,index_aux) =
     &          ddv_aux(max_index,index_aux) +
     &          dve_pi_aux(ifjs_aux)
            endif
          enddo
        endif
      endif
      if(ifexcited.gt.0)
     &  call excitation_pi(ifexcited,
     &  ifpl, ifpi_local, ifmodified, ifnr, ifsame_zero_abundances,
     &  inv_aux, iextraoff, naux, inv_ion, ifh2, ifh2plus,
     &  partial_elements, n_partial_elements, ion_end, 
     &  tl, izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &  bi, plop, plopt, plopt2,
     &  r_ion3, nion, nions, r_neutral, nelements,
     &  extrasum, extrasumf, extrasumt, nextrasum,
     &  xextrasum, xextrasumf, xextrasumt,
     &  dv, dvf, dvt, ddv_aux)
      if(debug_dauxdv) then
C          specify change in dv as a function of fl and tl
        dv(partial_ions(itemp1)) = dv(partial_ions(itemp1)) +
     &    delta1
        dv(partial_ions(itemp2)) = dv(partial_ions(itemp2)) +
     &    delta2
      endif
C      add in Planck-Larkin occupation probability effect
      if(ifpl.eq.1) then
        do index_ion = 1,max_index
          ion = partial_ions(index_ion)
          dv(ion) = dv(ion) + dv_pl(ion)
          if(ifnr03) then
            dvt(ion) = dvt(ion) + dv_plt(ion)
          endif
        enddo
      endif
C      eos_warm_step always returns n form rather than nu = n/(rho*avagadro)
C      form for h2, h2plus, h, hd, he, and extrasum.
      ifnuform = .false.
      call eos_calc(ifexcited, ifsame_zero_abundances,
     &  ifnuform, ifsame_under, ifnr, inv_ion, max_index,
     &  partial_elements, n_partial_elements, ion_end,
     &  ifionized, if_pteh, if_mc, ifreducedmass, 
     &  ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,
     &  ifpl, ifmodified, ifh2, ifh2plus,
     &  izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &  eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &  r_ion3, r_neutral, nelements, 
     &  ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &  rhostar(1), rhostar(2), rhostar(3),
     &  ne, nef, net, sion, sionf, siont, uion,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het,
     &  h_dv, hd_dv, he_dv, fion, fionf, fion_dv,
     &  sum0, sum0f, sum0t, sum0_dv,
     &  sum2, sum2f, sum2t, sum2_dv,
     &  extrasum, extrasumf, extrasumt, extrasum_dv,
     &  nextrasum, maxnextrasum,
     &  sumpl0, sumpl0f, sumpl0_dv,
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  rl, rf, rt, r_dv,
     &  h2, h2f, h2t, h2_dv, h2plus, h2plusf, h2plust, h2plus_dv,
     &  xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
      rho_new = exp(rl)
C      Calculate pressure and its aux and fl derivatives.
C      free electrons component
      pe_cgs = cpe*pstar(1)
C      free ions component
      p0 = cr*rho_new*t
      pion = full_sum0*p0 - (h2+h2plus)*boltzmann*t
      pionf = full_sum0*p0*rf - (h2f+h2plusf)*boltzmann*t
      piont = full_sum0*p0*(1.d0+rt) -
     &  (h2+h2plus+h2t+h2plust)*boltzmann*t
      do index = 1,max_index
        pion_dv(index) = full_sum0*p0*r_dv(index) -
     &    (h2_dv(index) + h2plus_dv(index))*boltzmann*t
      enddo
C      exchange component
      call exchange_pressure(rhostar, pstar, nstar, pex, pexf, pext)
C      pressure ionization component
      do index = 1,max_index
        ppi_dv(index) = 0.d0
      enddo
      do index_aux = 1,n_partial_aux
        ppi_daux(index_aux) = 0.d0
      enddo
      if(ifpi_local.eq.0) then
C        there is no free energy term from pressure ionization.
        ppi = 0.d0
        ppif = 0.d0
        ppit = 0.d0
      elseif(ifpi_local.eq.1) then
        call pteh_pi_pressure(
     &    full_sum1, rho_new, rf, rt, t, ne, nef, net,
     &    ppi, ppif, ppit, ppir)
        do index = 1,max_index
          ppi_dv(index) = ppir*r_dv(index)
        enddo
      elseif(ifpi_local.eq.2) then
C        ppi, ppif, ppit, and ppi2_daux defined as a result of
C        the call to fjs_pi_free done earlier in this routine.
        do index = 1,4
          index_aux = inv_aux(index)
          if(index_aux.gt.0)
     &      ppi_daux(index_aux) = ppi2_daux(index)
        enddo
      elseif(ifpi_local.eq.3.or.ifpi_local.eq.4) then
C        ppi, ppif, ppit, and ppi2_daux defined as a result of
C        the call to mdh_pi_pressure_free done earlier in this routine.
        do index = 1,nextrasum
          index_aux = inv_aux(index+iextraoff)
          if(index_aux.gt.0)
     &      ppi_daux(index_aux) = ppi2_daux(index)
        enddo
      endif
      if(ifexcited.gt.0) then
        call excitation_pi_pressure_free(
     &    t, rho_new, rf, rt, r_dv, max_index,
     &    pexcited, pexcitedf, pexcitedt, pexcited_dv,
     &    free_excited, free_excitedf, free_excited_dv)
      else
        pexcited = 0.d0
        pexcitedf = 0.d0
        free_excited = 0.d0
        free_excitedf = 0.d0
        do index = 1, max_index
          pexcited_dv(index) = 0.d0
!          free_excited_dv(index) = 0.d0
        enddo
      endif
C      N.B. there are no pressure terms for Planck-Larkin occupation
C      probability.
      pnorad = pe_cgs + pion + pex + ppi + pexcited
      pnoradf = pe_cgs*pstar(2) + pionf + pexf + ppif + pexcitedf
      pnoradt = pe_cgs*pstar(3) + piont + pext + ppit + pexcitedt
      do index = 1,max_index
        pnorad_dv(index) = pion_dv(index) + ppi_dv(index) +
     &    pexcited_dv(index)
      enddo
      do index_aux = 1,n_partial_aux
        pnorad_daux(index_aux) = ppi_daux(index_aux)
        do index = 1,max_index
          pnorad_daux(index_aux) = pnorad_daux(index_aux) +
     &      pnorad_dv(index)*ddv_aux(index,index_aux)
        enddo
      enddo
      pnorad = pnorad + dpcoulomb
      pnoradf = pnoradf + dpcoulombf
      pnoradt = pnoradt + dpcoulombt
      if(if_pteh.ne.1) then
        if(if_mc.eq.1) then
C          sum0 and sum2 calculated according to following formulas:
C          sum0 = sumion0 + h2plus + h*(1.d0+hcon_mc) +
C     &      hd + he
C          sum2 = sumion2 + h2plus + h*(1.d0+hcon_mc) +
C     &      hd + he*(4.d0+hecon_mc)
C          where h, hd, he, and h2plus are the first, second, third
C          and fifth auxiliary variables,
C          hcon_mc and hecon_mc are constants, and
C          sumion0 and sumion2 are constants times exp(rl), where rl
C          is the 4th auxiliary variable.
C          partial wrt h.
          index_aux = inv_aux(1)
          if(index_aux.gt.0)
     &      pnorad_daux(index_aux) = pnorad_daux(index_aux) +
     &      (1.d0+hcon_mc)*(dpcoulomb0+dpcoulomb2)
C          partial wrt hd.
          index_aux = inv_aux(2)
          if(index_aux.gt.0)
     &      pnorad_daux(index_aux) = pnorad_daux(index_aux) +
     &      (dpcoulomb0+dpcoulomb2)
C          partial wrt he.
          index_aux = inv_aux(3)
          if(index_aux.gt.0)
     &      pnorad_daux(index_aux) = pnorad_daux(index_aux) +
     &      dpcoulomb0 + (4.d0+hecon_mc)*dpcoulomb2
C          partial wrt (old) rl.
          index_aux = inv_aux(4)
          if(index_aux.gt.0)
     &      pnorad_daux(index_aux) = pnorad_daux(index_aux) +
     &      sumion0*dpcoulomb0 + sumion2*dpcoulomb2
C          partial wrt h2plus.
          index_aux = inv_aux(5)
          if(index_aux.gt.0)
     &      pnorad_daux(index_aux) = pnorad_daux(index_aux) +
     &      (dpcoulomb0+dpcoulomb2)
        else
C          sum0
          index_aux = inv_aux(6)
          if(index_aux.gt.0)
     &      pnorad_daux(index_aux) = pnorad_daux(index_aux) +
     &      dpcoulomb0
C        sum2
          index_aux = inv_aux(7)
          if(index_aux.gt.0)
     &      pnorad_daux(index_aux) = pnorad_daux(index_aux) +
     &      dpcoulomb2
        endif
      endif
C      Calculate the free-energy/unit mass
C      which is equivalent to u - Ts, where u is the energy per unit
C      mass and s is the entropy per unit mass.
C
C      First, collect the parts which are in per unit volume units.
!      if(ifrad.ge.1) then
!        free_rad = -prad_const*t**4
!      else
!        free_rad = 0.d0
!      endif
!      free_e = -cpe*pstar(1) + eta*n_e*boltzmann*t
!      free_ef = -cpe*pstar(1)*pstar(2) +
!     &  (wf + eta*rhostar(2))*n_e*boltzmann*t
!      call master_coulomb_free(if_pteh,
!     &  free_coulomb, free_coulombf, free_coulomb0, free_coulomb2)
!      call exchange_free(rhostar, pstar, nstar, free_ex, free_exf)
C************************************************
C      pressure ionization
C      deal with free_pi_dv and free_pi_daux separately for now and transform
C      free_pi_dv into equivalent free_pi_daux form later.
!      do index = 1,max_index
!        free_pi_dv(index) = 0.d0
!      enddo
!      do index_aux = 1,n_partial_aux
!        free_pi_daux(index_aux) = 0.d0
!      enddo
!      if(ifpi_local.eq.0) then
!C        there is no free energy term from pressure ionization.
!        free_pi = 0.d0
!        free_pif = 0.d0
!      elseif(ifpi_local.eq.1) then
!C        Note.  Unlike other free_pi which are per volume, this one is
!C        per mass to be accounted for later.
!        call pteh_pi_free(full_sum1, rho_new, rf, t, ne, nef,
!     &    free_pi, free_pif, free_pir)
!        do index = 1,max_index
!          free_pi_dv(index) = free_pir*r_dv(index)
!        enddo
!      elseif(ifpi_local.eq.2) then
!C        free_pi, free_pif, and free_pi2_daux defined as a result of
!C        the call to fjs_pi_free done earlier in this routine.
!        do index = 1,4
!          index_aux = inv_aux(index)
!          if(index_aux.gt.0)
!     &      free_pi_daux(index_aux) = free_pi2_daux(index)
!        enddo
!      elseif(ifpi_local.eq.3.or.ifpi_local.eq.4) then
!C        free_pi, free_pif, and free_pi2_daux defined as a result of
!C        the call to mdh_pi_pressure_free done earlier in this routine.
!        do index = 1,nextrasum
!          index_aux = inv_aux(index+iextraoff)
!          if(index_aux.gt.0)
!     &      free_pi_daux(index_aux) = free_pi2_daux(index)
!        enddo
!      endif
!      if(ifpi_local.ge.2) then
!C        for this case convert from free_pi per unit volume to per unit
!C        mass.
!        free_pif = (free_pif - free_pi*rf)/rho_new
!        free_pi = free_pi/rho_new
!C        do this transformation from POV that free_pi a function of both
!C        dv and (old) aux (actually it is one or the other with the unused
!C        derivatives are zeroed), while rho_new is a function of dv alone.
!        do index = 1,max_index
!          free_pi_dv(index) = (free_pi_dv(index) -
!     &      free_pi*r_dv(index))/rho_new
!        enddo
!        do index_aux = 1,n_partial_aux
!          free_pi_daux(index_aux) = free_pi_daux(index_aux)/rho_new
!        enddo
!      endif
!C      Second, collect the parts which are naturally in per unit mass
!C      units.
!C************************************************
!C      free_excited and derivatives available from call to
!C      excitation_pi_pressure_free done earlier in this routine
!C************************************************
!C      Planck-Larkin component/unit mass from free_eos_detailed
!C      correction to spi and upi and above formula.
!      if(ifpl.eq.1) then
!        free_pl = -cr*t*sumpl0
!        free_plf = -cr*t*sumpl0f
!        do index = 1,max_index
!          free_pl_dv(index) = -cr*t*sumpl0_dv(index)
!        enddo
!      else
!        free_pl = 0.d0
!        free_plf = 0.d0
!        do index = 1,max_index
!          free_pl_dv(index) = 0.d0
!        enddo
!      endif
C      Final transformation of fion:
C      add in (1.d0 + 1.5d0*tl - rl)) term summed over all species,
C      multiply by RT to obtain free-energy in ergs, and put in
C      hydrogen species energy zero point shift of h2diss for the
C      two diatomic species and h2diss/2 for the monatomic species.
      ni = full_sum0 - (h2+h2plus)/(rho_new*avogadro)
      nif = ((h2+h2plus)*rf - (h2f+h2plusf))/(rho_new*avogadro)
      fion = cr*t*(fion - ni*(1.d0 + 1.5d0*tl - rl)) +
     &  0.5d0*(c2*cr)*h2diss*eps(1)
      fionf = cr*t*(fionf - nif*(1.d0 + 1.5d0*tl - rl) + ni*rf)
      do index = 1,max_index
        ni_dv(index) = ((h2+h2plus)*r_dv(index) -
     &    (h2_dv(index)+h2plus_dv(index)))/(rho_new*avogadro)
        fion_dv(index) = cr*t*(fion_dv(index) -
     &    ni_dv(index)*(1.d0 + 1.5d0*tl - rl) + ni*r_dv(index))
      enddo
!      free = free_pi +
!     &    (free_rad + free_e + free_coulomb + free_ex)/rho_new +
!     &    free_excited + free_pl + fion
!        freef = free_pif +
!     &    (free_ef + free_coulombf + free_exf)/rho_new -
!     &    (free_rad + free_e + free_coulomb + free_ex)*rf/rho_new +
!     &    free_excitedf + free_plf + fionf
!      do index = 1,max_index
!        free_dv(index) = free_pi_dv(index) -
!     &    (free_rad + free_e + free_coulomb + free_ex)*
!     &    r_dv(index)/rho_new +
!     &    free_excited_dv(index) + free_pl_dv(index) + fion_dv(index)
!      enddo
!      do index_aux = 1,n_partial_aux
!        free_daux(index_aux) = free_pi_daux(index_aux)
!        do index = 1,max_index
!          free_daux(index_aux) = free_daux(index_aux) +
!     &      free_dv(index)*ddv_aux(index,index_aux)
!        enddo
!      enddo
!      if(if_pteh.ne.1) then
!        if(if_mc.eq.1) then
!C          sum0 and sum2 calculated according to following formulas:
!C          sum0 = sumion0 + h2plus + h*(1.d0+hcon_mc) +
!C     &      hd + he
!C          sum2 = sumion2 + h2plus + h*(1.d0+hcon_mc) +
!C     &      hd + he*(4.d0+hecon_mc)
!C          where h, hd, he, and h2plus are the first, second, third
!C          and fifth auxiliary variables,
!C          hcon_mc and hecon_mc are constants, and
!C          sumion0 and sumion2 are constants times exp(rl), where rl
!C          is the 4th auxiliary variable.
!C          partial wrt h.
!          index_aux = inv_aux(1)
!          if(index_aux.gt.0)
!     &      free_daux(index_aux) = free_daux(index_aux) +
!     &      (1.d0+hcon_mc)*(free_coulomb0+free_coulomb2)/rho_new
!C          partial wrt hd.
!          index_aux = inv_aux(2)
!          if(index_aux.gt.0)
!     &      free_daux(index_aux) = free_daux(index_aux) +
!     &      (free_coulomb0+free_coulomb2)/rho_new
!C          partial wrt he.
!          index_aux = inv_aux(3)
!          if(index_aux.gt.0)
!     &      free_daux(index_aux) = free_daux(index_aux) +
!     &      (free_coulomb0 + (4.d0+hecon_mc)*free_coulomb2)/rho_new
!C          partial wrt (old) rl.
!          index_aux = inv_aux(4)
!          if(index_aux.gt.0)
!     &      free_daux(index_aux) = free_daux(index_aux) +
!     &      (sumion0*free_coulomb0 + sumion2*free_coulomb2)/rho_new
!C          partial wrt h2plus.
!          index_aux = inv_aux(5)
!          if(index_aux.gt.0)
!     &      free_daux(index_aux) = free_daux(index_aux) +
!     &      (free_coulomb0+free_coulomb2)/rho_new
!        else
!C          sum0
!          index_aux = inv_aux(6)
!          if(index_aux.gt.0)
!     &      free_daux(index_aux) = free_daux(index_aux) +
!     &      free_coulomb0/rho_new
!C          sum2
!          index_aux = inv_aux(7)
!          if(index_aux.gt.0)
!     &      free_daux(index_aux) = free_daux(index_aux) +
!     &      free_coulomb2/rho_new
!        endif
!      endif
      do index_aux = 1,n_partial_aux
C        take derivative with respect to log(auxold) or log(-auxold)
        iaux = partial_aux(index_aux)
        if(allow_log(index_aux).or.
     &      allow_log_neg(index_aux)) then
          pnorad_daux(index_aux) =
     &      pnorad_daux(index_aux)*aux_old(iaux)
        endif
      enddo
      pnorad_daux(n_partial_aux+1) = pnoradf
      pnorad_daux(n_partial_aux+2) = pnoradt
!      !temporary
!      if(.false.) then
!      do index_aux = 1,n_partial_aux
!C        take derivative with respect to log(auxold) or log(-auxold)
!        iaux = partial_aux(index_aux)
!        if(allow_log(index_aux).or.
!     &      allow_log_neg(index_aux)) then
!          free_daux(index_aux) =
!     &      free_daux(index_aux)*aux_old(iaux)
!        endif
!      enddo
!      endif
!      free_daux(n_partial_aux+1) = freef
!      free_daux(n_partial_aux+2) = 0.d0
      !temporary (note quantity output is scaled the same as eos_bfgs result
      ! for scaled free energy).
!      write(0,'(a,/,(1p2d25.15))')
!     &  'eos_warm_step: fl, tl, pnorad, (scaled) free = ',
!     &  log(f), tl, pnorad, free/(cr*t*full_sum0)
!      pe_cgs = pe_cgs/pnorad
!      dpcoulomb = dpcoulomb/pnorad
!      pex = pex/pnorad
!      ppi = ppi/pnorad
!      pexcited = pexcited/pnorad
!      pion = pion/pnorad
!      free_rad = free_rad/rho_new/free
!      free_e = free_e/rho_new/free
!      free_coulomb = free_coulomb/rho_new/free
!      free_ex = free_ex/rho_new/free
!      free_pi = free_pi/free
!      free_excited = free_excited/free
!      free_pl = free_pl/free
!      fion = fion/free
!      write(0,*)
!     &  'fractions of e, coulomb, ex, pi, excited, '//
!     &  'and ion components of pnorad = '
!      write(0,*) pe_cgs, dpcoulomb, pex, ppi,
!     &  pexcited, pion
!      write(0,*)
!     &  'total e, coulomb, ex, pi, excited, '//
!     &  ' and ion components of pnorad = '
!      write(0,*) pnorad*pe_cgs, pnorad*dpcoulomb,
!     &  pnorad*pex, pnorad*ppi, pnorad*pexcited,
!     &  pnorad*fion
!      write(0,*)
!     &  'fractions of rad, e, coulomb, ex, pi, excited, '//
!     &  'pl, and ion components of free = '
!      write(0,*) free_rad, free_e, free_coulomb, free_ex, free_pi,
!     &  free_excited, free_pl, fion
!      write(0,*)
!     &  'total rad, e, coulomb, ex, pi, excited, '//
!     &  'pl, and ion components of free = '
!      write(0,'(1p4d25.15)')
!     &  free*free_rad, free*free_e, free*free_coulomb,
!     &  free*free_ex, free*free_pi, free*free_excited,
!     &  free*free_pl, free*fion
!      fion = free*fion
!      call free_non_ideal_calc(
!     &  maxnextrasum,
!     &  ifpi_local, ifnr, ifmodified,
!     &  ifdvzero, inv_aux, naux, partial_ions, inv_ion,
!     &  partial_elements, n_partial_elements,
!     &  ion_end, mion_end,
!     &  ifdv, ifrad,
!     &  r_ion3, nion, nions, r_neutral, nelements,
!     &  tl, f, wf, eta, rhostar, pstar,
!     &  nux, nuy, nuz,
!     &  if_pteh, full_sum0, full_sum1, full_sum2,
!     &  if_mc, hcon_mc, hecon_mc, sum0_mc, sum2_mc,
!     &  ifcoulomb_mod, if_dc,
!     &  ifexcited, ifpl, ifsame_zero_abundances, ifh2, ifh2plus,
!     &  izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax, bi,
!     &  plop, plopt, plopt2,
!     &  sumpl0, sumpl0f,
!     &  h, hf, ht,
!     &  hd, hdf, hdt,
!     &  he, hef, het,
!     &  rl, rf, rt,
!     &  h2plus, h2plusf, h2plust,
!     &  sum0, sum0f, sum0t,
!     &  sum2, sum2f, sum2t,
!     &  extrasum, extrasumf, extrasumt, nextrasum,
!     &  xextrasum, xextrasumf, xextrasumt,
!     &  free_pi, free_pif, free_pit,
!     &  free_rad,
!     &  free_e,
!     &  free_coulomb,
!     &  free_ex,
!     &  free_excited,
!     &  free_pl
!     &  )
!      free = free_pi +
!     &    free_rad + free_e + free_coulomb + free_ex +
!     &    free_excited + free_pl + fion
!      write(0,'(a,/,(1p2d25.15))')
!     &  'eos_warm_step2: fl, tl, pnorad, (scaled) free = ',
!     &  log(f), tl, pnorad, free/(full_sum0*cr*t)
!      write(0,*)
!     &  'total rad, e, coulomb, ex, pi, excited, '//
!     &  'pl, and ion components of free = '
!      write(0,'(1p4d25.15)')
!     &  free_rad, free_e, free_coulomb,
!     &  free_ex, free_pi, free_excited,
!     &  free_pl, fion
      end
