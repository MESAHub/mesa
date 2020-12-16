C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: eos_calc.f 831 2008-06-30 00:22:11Z airwin $
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
      subroutine eos_calc(ifexcited, ifsame_zero_abundances,
     &  ifnuform, ifsame_under, ifnr,
     &  inv_ion, max_index,
     &  partial_elements, n_partial_elements, ion_end,
     &  ifionized, if_pteh, if_mc, ifreducedmass,
     &  ifsame_abundances, ifmtrace, iatomic_number, ifpi, ifpl,
     &  ifmodified, ifh2, ifh2plus,
     &  izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &  eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &  r_ion3, r_neutral, nelements,
     &  ifelement, dvzero, dv, dvf, dvt, nion, nions, re, ref, ret,
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
C         ifexcited > 0 means_use_excited states (must have Planck-Larkin or
C           ifpi = 3 or 4).
C            0 < ifexcited < 10 means_use_approximation to explicit summation
C           10 < ifexcited < 20 means_use_explicit summation
C           mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules)
C             and helium.
C           mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
C           mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
C         ifnuform: usually set to false so that h2, h2plus, h, hd, he, and
C           extrasum returned in standard number density form.
C           however, on last call to eos_calc for a given fl, tl, this
C           quantity should be set to true so that these quantities
C           returned in nu = n/(rho*avogadro) form.  This reduces
C           the significance loss in the entropy due to pressure ionization.
C         ifsame_under = logical variable passed to ionize to control whether
C           to_use_the same underflow limits or recalculate them.
C           used for removing small discontinuities in eos when attempting
C           to obtain converged solution.
C         ifnr = 0
C           calculate derivatives of ne, sion, sumpl1, h2,
C           h, hd, he, sum0, sum2, extrasum, r, h2plus, and xextrasum
C           wrt f, t and sumpl0 wrt f with no other variables fixed, i.e., use
C           input dvf and dvt and chain rule.
C         N.B. chain of calling routines depend on the assertion that ifnr = 0
C         means no *_dv variables should be read or written.
C         ifnr = 1
C           calculate derivatives of
C           h, hd, he, sum0, sum2, extrasum, r, h2plus, and xextrasum
C           n.b. list of variables is subset of those listed for ifnr = 0
C           because we only need dv derivatives for variables used to calculate
C           (output) auxiliary variables.
C         ifnr = 2 (should not occur)
C         ifnr = 3 same as combination of ifnr = 0 and ifnr = 1.
C           Note this is quite different in detail from ifnr = 3 interpretation
C           for many other routines, but the general motivation is the same
C           for all ifnr = 3 results; calculate both f and t derivatives and
C           other derivatives.
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
C         re, ref, ret: ne_star fermi-dirac integral and derivatives
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
C         sumpl0, sumpl1 and sumpl2 = weighted sums over
C         non-H nu(i) = n(i)/(rho/H).
C           for sumpl0 sum is over Planck-Larkin occupation probabilities.
C           for sumpl1 sum is over Planck-Larkin occupation probabilities +
C           d ln w/d ln T
C           for sumpl2, sum is over Planck-Larkin d ln w/d ln T
C         sumpl0f, sumpl0_dv = derivative of sumpl0 wrt ln f and dv.
C         sumpl1f and sumpl1t = derivatives of sumpl1 wrt lnf and lnt.
C         rl, rf, rt are the ln mass density and derivatives
C         h2, h2f, h2t, h2_dv are n(H2) and derivatives.
C         h2plus, h2plusf, h2plust are n(H2+) and derivatives.
C         h2plus_dv is the h2plus derivatives wrt dv.  n.b. this vector only
C           returned if molecular hydrogen is calculated.  The calling
C           routine uses it only if if_mc.eq.1.and.ifh2plus.gt.0.
C         the following are returned only if ifexcited.gt.0 and
C           ifpi.eq.3.or4.
C         xextrasum(4) is the *negative* sum nuvar/(1 + qratio)*
C         partial qratio/partial extrasum(k), for k = 1, 2, 3, and nextrasum-1.
C         xextrasumf(4), xextrasumt(4), xextrasum_dv(nions+2,4) are the
C         fl and tl derivatives (ifnr.eq.0) or dv derivatives (ifnr.eq.1)
      implicit none
C      required so that can store nuh2 and nuh2plus for eos_sum_calc and
C      eos_free_calc.
      include 'aux_scale.h'
      include 'constants.h'
      include 'nuvar.h'
      include 'statistical_weights.h'
      logical ifnuform, ifsame_under, ifabh, ifabg,
     &  ifaah2plus, ifaah2
      integer ifexcited, ifsame_zero_abundances, ifnr,
     &  n_partial_elements, partial_elements(n_partial_elements+2),
     &  ifionized, if_pteh, if_mc, ifreducedmass, 
     &  ifsame_abundances, ifmtrace, iatomic_number, ifpi, ifpl,
     &  ifmodified, ifh2, ifh2plus,
     &  izlo, izhi,  nmin(izhi), nmin_max(izhi),
     &  nions, nmin_species(nions+2), nmax,
     &  nion(nions), inv_ion(nions+2), max_index, index,
     &  index1, index2, index3, indexh, indexh2, indexh2plus,
     &  nelements, ion_end(nelements), ifelement(nelements),
     &  nextrasum, maxnextrasum, nions_local, nxextrasum
      parameter(nxextrasum = 4)
      parameter(nions_local = 295)
      double precision  eps(nelements), tl, tc2, bi(nions+2), h2diss,
     &  plop(nions+2), plopt(nions+2), plopt2(nions+2), bmin(izhi),
     &  dvzero(nelements), dv(nions+2), dvf(nions+2), dvt(nions+2),
     &  re, ref, ret,
     &  ne, nef, net, sion, sionf, siont, uion,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het,
     &  h_dv(nions+2), hd_dv(nions+2), he_dv(nions+2),
     &  fionh, fionhf, fionh_dv(nions_local+2),
     &  fion, fionf, fion_dv(nions+2),
     &  sum0, sum0f, sum0t, sum0_dv(nions+2),
     &  sum2, sum2f, sum2t, sum2_dv(nions+2),
     &  r_ion3(nions+2), r_neutral(nelements+2),
     &  extrasum(maxnextrasum), extrasumf(maxnextrasum),
     &  extrasumt(maxnextrasum), extrasum_dv(nions+2, maxnextrasum),
     &  lextrasum(9), lextrasumf(9),
     &  lextrasumt(9), lextrasum_dv(nions_local+2, 9),
     &  xextrasum(nxextrasum), xextrasumf(nxextrasum),
     &  xextrasumt(nxextrasum), xextrasum_dv(nions+2,nxextrasum),
     &  sumpl0, sumpl0f, sumpl0_dv(nions+2),
     &  sumpl1, sumpl1f, sumpl1t,
     &  sumpl2, sharg, rl, rf, rt,
     &  r_dv(nions+2), h2, h2f, h2t, h2plus, h2plusf, h2plust
      double precision  g, gf, gt, g_dv(nions_local+2),
     &  lognug, nug, nugf, nugt, nug_dv(nions_local+2),
     &  nuh, nuhf, nuht, 
     &  nuh2, nuh2f, nuh2t, nuh2_dv(nions_local+2),
     &  nuh2plus, nuh2plusf, nuh2plust, nuh2plus_dv(nions_local+2),
     &  h2_dv(nions+2), h2plus_dv(nions+2),
     &  sh, shf, sht, hne, hnef, hnet, hne_dv(nions_local+2),
     &  en, rmue, rho,
     &  x, xf, xt, x_dv(nions_local+2),
     &  tlold, qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt,
     &  hequil, hequilf, hequilt, hequil_dv,
     &  h2equil, h2equilf, h2equilt, h2equilt0, h2equil_dv,
     &  h2plusequil, h2plusequilf, h2plusequilt, h2plusequil_dv, 
     &  logalpha,
     &  logmass_h2plus, logmass_h2, logmass_h, logmass_hplus,
     &  logqtl_const_h2plus, logqtl_const_h2,
     &  logqtl_const_h, logqtl_const_hplus,
     &  h2_log, h2plus_log,
     &  ablog, ablogf, ablogt, ablog_dv(nions_local+2),
     &  abh, abhf, abht, abh_dv, 
     &  abg, abgf, abgt, abg_dv, abexp,
     &  abdlog, abdlogf, abdlogt, abdlog_dv(nions_local+2),
     &  abdh, abdhf, abdht, abdh_dv, 
     &  abdg, abdgf, abdgt, abdg_dv,
     &  ac, acf, act, 
     &  aalog, aalogf, aalogt, aalog_dv,
     &  aah2plus, aah2plusf, aah2plust, aah2plus_dv, 
     &  aah2, aah2f, aah2t, aah2_dv, aaexp,
     &  cprime, cprimef, cprimet, cprime_dv,
     &  xprimelog, xprimelogf, xprimelogt,
     &  xprimelog_dv(nions_local+2),
     &  rpower, rpowerh2, rpowerh2plus, constant_sh
      integer iextrasum, jextrasum, iffirst
      data iffirst/1/
C      something ridiculous
      data tlold/1.d-300/
      logical ifnr03, ifnr13
      integer ifreducedmassold
C      need invalid value
      data ifreducedmassold/-1/
C      Must be consistent with atomic_mass(1) in ionize to produce consistent
C      results with and without molecules.
      double precision atomic_mass_h
      data atomic_mass_h/1.007825035d0/
C      variables associated with finding the dominant H species.
      double precision maxh, second_maxh, lnerrcrit
C      lnerrcrit is set so that ln(1 + exp(-lnerrcrit)) ~
C      exp(-lnerrcrit) ~ 1.d-17
      parameter(lnerrcrit = 39.d0)
      integer index_maxh
      logical ifpi34
      save
C      sanity checking:
      if(ifnr.lt.0.or.ifnr.eq.2.or.ifnr.gt.3) stop
     &  'eos_calc: ifnr must be 0, 1, or 3'
      if(nions.ne.nions_local) stop 'eos_calc: bad nions'
C      ifnuform should be true only on last call to eos_calc which coincides
C      with ifnr = 0.
      if(ifnuform.and.ifnr.ne.0) 
     &  stop 'eos_calc: bad combination of ifnuform and ifnr'
      if(ifionized.eq.2.and.ifexcited.gt.0)
     &  stop 'eos_calc: bad combination of ifionized and ifexcited'
      if(iffirst.eq.1.and.ifsame_abundances.eq.1)
     &  stop 'eos_calc: invalid ifsame_abundances on first call'
C      n.b. the atomic masses and the ground statistical weight for H and H+
C      used below *must* be identical with the values used in ionize.
C      Otherwise, you will introduce a discontinuity when switching from
C      negligible molecular hydrogen to a calculation without molecular
C      hydrogen.
      if(iffirst.eq.1) then
        iffirst = 0
C        alpha =  (2 pi k/[Na h^2])^(-3/2) Na
C        (see documentation of alpha^2 in constants.h).
        logalpha = 0.5d0*log(alpha2)
        logmass_h = log(atomic_mass_h)
        logmass_h2 = log(2.d0*atomic_mass_h)
        h2_log = 1.5d0*(logmass_h2 - 2.d0*logmass_h)
      endif
      if(ifreducedmass.ne.ifreducedmassold) then
        ifreducedmassold = ifreducedmass
        if(ifreducedmass.eq.1) then
          logmass_hplus = log(atomic_mass_h - electron_mass)
!          logmass_h2plus = log(2.d0*atomic_mass_h - electron_mass)
        else
          logmass_hplus = log(atomic_mass_h)
!          logmass_h2plus = log(2.d0*atomic_mass_h)
        endif
C        to be consistent with the slightly lame free-energy model
C        underlying the entropy and energy calculation which assumes
C        the mass of H2 and H2+ are equal to each other.
        logmass_h2plus = logmass_h2
        h2plus_log = 1.5d0*(logmass_h2plus - logmass_h2)
C        drop partition function from all molecular values since that not
C        constant with temperature.
        logqtl_const_h2plus = 1.5d0*logmass_h2plus - logalpha
        logqtl_const_h2 = 1.5d0*logmass_h2 - logalpha
        logqtl_const_h =
     &    log(dble(iqneutral(1))) + 1.5d0*logmass_h - logalpha
        logqtl_const_hplus =
     &    log(dble(iqion(1))) + 1.5d0*logmass_hplus - logalpha
      endif
C      calculate constant component of hydrogen entropy to be used
C      in molecular case only (see below).
      if(ifsame_abundances.ne.1)
     &  constant_sh = eps(1)*(1.5d0*logmass_h +
     &  log(dble(iqneutral(1))) - logalpha)
      ifnr03 = ifnr.eq.0.or.ifnr.eq.3
      ifnr13 = ifnr.eq.1.or.ifnr.eq.3
      ifpi34 = ifpi.eq.3.or.ifpi.eq.4
      call ionize(ifexcited, ifsame_under, ifnr, inv_ion,
     &  partial_elements, n_partial_elements, ion_end,
     &  ifionized, if_pteh, if_mc, ifreducedmass, 
     &  ifsame_abundances, ifmtrace, iatomic_number, ifpi, ifpl,
     &  eps, tc2, bi, plop, plopt, plopt2,
     &  r_ion3, r_neutral, nelements, 
     &  ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &  hne, hnef, hnet, hne_dv, sion, sionf, siont, uion,
     &  g, gf, gt, sh, shf, sht,
     &  hequil, hequilf, hequilt,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het,
     &  h_dv, hd_dv, he_dv,
     &  fionh, fionhf, fionh_dv, fion, fionf, fion_dv,
     &  sum0, sum0f, sum0t, sum0_dv,
     &  sum2, sum2f, sum2t, sum2_dv,
     &  extrasum, extrasumf, extrasumt, extrasum_dv,
     &  nextrasum, maxnextrasum,
     &  sumpl0, sumpl0f, sumpl0_dv,
     &  sumpl1, sumpl1f, sumpl1t, sumpl2)
      if(ifionized.ne.2.and.ifnr13) then
C        only:
C        h_dv(inv_ion(1)), 
C        hd_dv(inv_ion(2)), hd_dv(inv_ion(3)), and
C        he_dv(inv_ion(2)), he_dv(inv_ion(3)) defined by
C        ionize, but further programming below skips the
C        undefined values for these arrays so there is no
C        problem.
C        on the other hand hne_dv, sum0_dv, sum2_dv, extrasum_dv, fion_dv,
C        and sumpl0_dv are undefined by ionize for index = inv_ion(nions+1)
C        and inv_ion(nions+1).  
C        This block zeroes these values if required.
        do index = nions+1, nions+2
          index1 = inv_ion(index)
          if(index1.gt.0) then
            hne_dv(index1) = 0.d0
            fion_dv(index1) = 0.d0
            if(ifpl.eq.1)
     &        sumpl0_dv(index1) = 0.d0
            if(.not.(if_mc.eq.1.or.if_pteh.eq.1)) then
              sum0_dv(index1) = 0.d0
              sum2_dv(index1) = 0.d0
            endif
            if(ifpi34) then
              do iextrasum = 1, nextrasum
                extrasum_dv(index1,iextrasum) = 0.d0
              enddo
            endif
          endif
        enddo
      endif
C      rho/mu_e = n_e H = n_e/N_A = cd*re = rmue
      rmue = cd*re
      if(ifionized.eq.2) then
C        everything fully ionized.
C        ignore ifnr (i.e. forget dv derivatives) since no NR iteration
C        required with full ionization.
C        n.b. for full ionization, ionize includes hydrogen effects
C        in hne, sion, uion, and fion.
        ne = hne
        nef = 0.d0
        net = 0.d0
        en = 1.d0/ne
        rho = rmue*en
        rl = log(rho)
        rt = ret-en*net
        rf = ref-en*nef
        if(ifnuform) then
C          calculate nu forms of h, he
          h = eps(1)
          hf = 0.d0
          ht = 0.d0
          he = eps(2)
          hef = 0.d0
          het = 0.d0
        else
C          should never reach this block of code because with
C          full ionization in awieos_detailed should have
C          straight through option which should yield 
C          ifnuform = .true.
          stop 'eos_calc: bad ifnuform setting for full ionization'
C          calculate number densities
!          h = eps(1)*rho*avogadro
!          hf = eps(1)*rf
!          ht = eps(1)*rt
!          he = eps(2)*rho*avogadro
!          hef = eps(2)*rf
!          het = eps(2)*rt
        endif
        hd = 0.d0
        hdf = 0.d0
        hdt = 0.d0
C        sum0, sum2, extrasum, sumpl0, sumpl1, sumpl2
C        not needed.
        h2 = 0.d0
        h2f = 0.d0
        h2t = 0.d0
        h2plus = 0.d0
        h2plusf = 0.d0
        h2plust = 0.d0
C        h2plus_dv not needed.
      elseif(eps(1).le.0.d0) then
C        partially ionized but no hydrogen (or molecular hydrogen).
        ne = hne
C        sion and derivatives, fion and derivative, and uion
C        passed back as is from ionize because no h component.
        en = 1.d0/ne
C        rho/mu_e = n_e H = cd*re = rmue
        rho = rmue*en
        rl = log(rho)
        h2 = 0.d0
        h2plus = 0.d0
        if(ifnr03) then
          nef = hnef
          net = hnet
          rf = ref-en*nef
          rt = ret-en*net
          h2f = 0.d0
          h2t = 0.d0
          h2plusf = 0.d0
          h2plust = 0.d0
        endif
        if(ifnr13) then
          do index = 1, max_index
            r_dv(index) = -en*hne_dv(index)
          enddo
        endif
C        h and derivatives irrelevant if zero abundance,
C        but force to be zero because previous
C        calls may have had non-zero abundance.
        h = 0.d0
        if(ifpi.eq.2.or.if_mc.eq.1) then
          if(ifnr03) then
            hf = 0.d0
            ht = 0.d0
          endif
        endif
C        convert hd, he, and derivatives to number densities.
        if(eps(2).gt.0.d0) then
C          if ifnuform is true then ifnr is zero and leave hd, he
C          and their derivatives in the nu form returned by
C          ionize.
          if(.not.ifnuform) then
            hd = hd*rho*avogadro
            he = he*rho*avogadro
            if(ifpi.eq.2.or.if_mc.eq.1) then
              if(ifnr03) then
                hdf = hdf*rho*avogadro + rf*hd
                hdt = hdt*rho*avogadro + rt*hd
                hef = hef*rho*avogadro + rf*he
                het = het*rho*avogadro + rt*he
              endif
              if(ifnr13) then
                index2 = inv_ion(2)
                index3 = inv_ion(3)
                do index = 1, max_index
                  if(index.eq.index2.or.index.eq.index3) then
                    hd_dv(index) = hd_dv(index)*rho*avogadro +
     &                hd*r_dv(index)
                    he_dv(index) = he_dv(index)*rho*avogadro +
     &                he*r_dv(index)
                  else
                    hd_dv(index) = hd*r_dv(index)
                    he_dv(index) = he*r_dv(index)
                  endif
                enddo
              endif
            endif
          endif
        else
          hd = 0.d0
          he = 0.d0
          if(ifpi.eq.2.or.if_mc.eq.1) then
            if(ifnr03) then
              hdf = 0.d0
              hef = 0.d0
              hdt = 0.d0
              het = 0.d0
            endif
          endif
        endif
        if(.not.(if_mc.eq.1.or.if_pteh.eq.1)) then
C          convert Coulomb sums to number densities.
          sum0 = sum0*(rho*avogadro)
          sum2 = sum2*(rho*avogadro)
          if(ifnr03) then
            sum0f = sum0f*(rho*avogadro) + sum0*rf
            sum0t = sum0t*(rho*avogadro) + sum0*rt
            sum2f = sum2f*(rho*avogadro) + sum2*rf
            sum2t = sum2t*(rho*avogadro) + sum2*rt
          endif
          if(ifnr13) then
            do index = 1, max_index
              sum0_dv(index) = sum0_dv(index)*(rho*avogadro) +
     &          sum0*r_dv(index)
              sum2_dv(index) = sum2_dv(index)*(rho*avogadro) +
     &          sum2*r_dv(index)
            enddo
          endif
        endif
        if(.not.ifnuform.and.ifpi34) then
C          if ifnuform is true then leave
C          extrasum and its derivative in the nu form returned by
C          ionize.
          do iextrasum = 1, nextrasum
C            convert to n form from nu = n/(rho*avogadro) form
            extrasum(iextrasum) = extrasum(iextrasum)*(rho*avogadro)
            if(ifnr03) then
              extrasumf(iextrasum) =
     &          extrasumf(iextrasum)*(rho*avogadro) +
     &          extrasum(iextrasum)*rf
              extrasumt(iextrasum) =
     &          extrasumt(iextrasum)*(rho*avogadro) +
     &          extrasum(iextrasum)*rt
            endif
            if(ifnr13) then
              do index = 1, max_index
                extrasum_dv(index,iextrasum) =
     &            extrasum_dv(index,iextrasum)*(rho*avogadro) +
     &            extrasum(iextrasum)*r_dv(index)
              enddo
            endif
          enddo
        endif
C        return sumpl0, sumpl1, derivatives, and sumpl2 as is from ionize
C        since no h component and keep in nu = n/(rho*avogadro) form.
C        h2plus_dv not needed.
      elseif(ifh2.eq.0) then
C        partially ionized including monatomic forms of hydrogen
C        keep ne, sion, sion derivatives, uion, 
C          sumpl0, sumpl0f, sumpl1, sumpl1 derivatives, 
C          and sumpl2 in nu = n/(rho*avogadro) form.
        ne = hne + h
        uion = uion + bi(1)*h
        sion = sion + sh
        fion = fion + fionh
        if(ifpl.eq.1) then
          sumpl0 = sumpl0 + g*plop(1)
          sumpl1 = sumpl1 + g*(plop(1) + plopt(1))
          sumpl2 = sumpl2 + g*plopt(1)
          if(ifnr03) then
            sumpl0f = sumpl0f + gf*plop(1)
            sumpl1f = sumpl1f + gf*(plop(1) + plopt(1))
            sumpl1t = sumpl1t + gt*(plop(1) + plopt(1)) +
     &        g*(plopt(1) + plopt2(1))
          endif
          if(ifnr13) then
C            as returned from ionize, g+h = eps(1), hence g_dv = -h_dv.
            do index = 1, max_index
              sumpl0_dv(index) = sumpl0_dv(index) -
     &          h_dv(index)*plop(1)
            enddo
          endif
        endif
        en = 1.d0/ne
C        rho/mu_e = n_e/avogadro = cd*re = rmue
        rho = rmue*en
        rl = log(rho)
        if(ifnuform) then
C          will need these for later.
          nug = g
          nugf = gf
          nugt = gt
          nuh = h
          nuhf = hf
          nuht = ht
        endif
        nuvar(1,1) = g
        nuvar(2,1) = h
C        store nuh2 and nuh2plus in the 3rd and 4th index
        nuvar(3,1) = 0.d0
        nuvar(4,1) = 0.d0
C        convert g, h, and derivatives to number densities.
        g = g*rho*avogadro
        h = h*rho*avogadro
        h2 = 0.d0
        h2plus = 0.d0
        if(ifnr03) then
          nef = hnef+hf
          net = hnet+ht
          sionf = sionf + shf
          siont = siont + sht
          fionf = fionf + fionhf
          rf = ref-en*nef
          rt = ret-en*net
          gf = gf*rho*avogadro + rf*g
          gt = gt*rho*avogadro + rt*g
          hf = hf*rho*avogadro + rf*h
          ht = ht*rho*avogadro + rt*h
          h2f = 0.d0
          h2t = 0.d0
          h2plusf = 0.d0
          h2plust = 0.d0
        endif
        if(ifnr13) then
          indexh = inv_ion(1)
          do index = 1, max_index
            if(index.eq.indexh) then
              r_dv(index) = -en*(hne_dv(index) + h_dv(index))
C              before conversion g+h = eps(1), hence g_dv = -h_dv.
              g_dv(index) = -h_dv(index)*rho*avogadro +
     &          g*r_dv(index)
C              the following expression has significance loss so do
C              another way.
!              h_dv(index) = h_dv(index)*rho*avogadro +
!     &          h*r_dv(index)
              h_dv(index) = rho*avogadro*en*h_dv(index)*hne -
     &          en*h*hne_dv(index)
              fion_dv(index) = fion_dv(index) + fionh_dv(index)
            else
              r_dv(index) = -en*hne_dv(index)
              g_dv(index) = g*r_dv(index)
              h_dv(index) = h*r_dv(index)
            endif
          enddo
        endif
C        convert hd, he, and derivatives to number densities.
        if(eps(2).gt.0.d0) then
C          if ifnuform is true then ifnr is zero and leave hd, he
C          and their derivatives in the nu form returned by
C          ionize.
          if(.not.ifnuform) then
            hd = hd*rho*avogadro
            he = he*rho*avogadro
            if(ifpi.eq.2.or.if_mc.eq.1) then
              if(ifnr03) then
                hdf = hdf*rho*avogadro + rf*hd
                hdt = hdt*rho*avogadro + rt*hd
                hef = hef*rho*avogadro + rf*he
                het = het*rho*avogadro + rt*he
              endif
              if(ifnr13) then
                index2 = inv_ion(2)
                index3 = inv_ion(3)
                do index = 1, max_index
                  if(index.eq.index2.or.index.eq.index3) then
                    hd_dv(index) = hd_dv(index)*rho*avogadro +
     &                hd*r_dv(index)
                    he_dv(index) = he_dv(index)*rho*avogadro +
     &                he*r_dv(index)
                  else
                    hd_dv(index) = hd*r_dv(index)
                    he_dv(index) = he*r_dv(index)
                  endif
                enddo
              endif
            endif
          endif
        else
          hd = 0.d0
          he = 0.d0
          if(ifpi.eq.2.or.if_mc.eq.1) then
            if(ifnr03) then
              hdf = 0.d0
              hef = 0.d0
              hdf = 0.d0
              het = 0.d0
            endif
          endif
        endif
        if(.not.(if_mc.eq.1.or.if_pteh.eq.1)) then
C          convert Coulomb sums to number densities and add in
C          effect of monatomic hydrogen species.
          sum0 = sum0*(rho*avogadro)
          sum2 = sum2*(rho*avogadro)
          if(ifnr03) then
            sum0f = sum0f*(rho*avogadro) + sum0*rf + hf*sum0_scale
            sum0t = sum0t*(rho*avogadro) + sum0*rt + ht*sum0_scale
            sum2f = sum2f*(rho*avogadro) + sum2*rf + hf*sum2_scale
            sum2t = sum2t*(rho*avogadro) + sum2*rt + ht*sum2_scale
          endif
          if(ifnr13) then
            do index = 1, max_index
              sum0_dv(index) = sum0_dv(index)*(rho*avogadro) +
     &          sum0*r_dv(index) + h_dv(index)*sum0_scale
              sum2_dv(index) = sum2_dv(index)*(rho*avogadro) +
     &          sum2*r_dv(index) + h_dv(index)*sum2_scale
            enddo
          endif
C          add in monatomic hydrogen species.
          sum0 = sum0 + h*sum0_scale
          sum2 = sum2 + h*sum2_scale
        endif
        if(ifnuform) then
C          replace n by nu form for g, h, and derivatives
          g = nug
          h = nuh
          gf = nugf
          gt = nugt
          hf = nuhf
          ht = nuht
        endif
        if(ifpi34) then
C          if ifnuform is true then leave
C          extrasum and its derivative in the nu form returned by
C          ionize.
          if(.not.ifnuform) then
            do iextrasum = 1, nextrasum
C              convert to n form from nu = n/(rho*avogadro) form
              extrasum(iextrasum) = extrasum(iextrasum)*(rho*avogadro)
              if(ifnr03) then
                extrasumf(iextrasum) =
     &            extrasumf(iextrasum)*(rho*avogadro) +
     &            extrasum(iextrasum)*rf
                extrasumt(iextrasum) =
     &            extrasumt(iextrasum)*(rho*avogadro) +
     &            extrasum(iextrasum)*rt
              endif
              if(ifnr13) then
                do index = 1, max_index
                  extrasum_dv(index,iextrasum) =
     &              extrasum_dv(index,iextrasum)*(rho*avogadro) +
     &              extrasum(iextrasum)*r_dv(index)
                enddo
              endif
            enddo
          endif
C          add in effect of neutral monatomic H in n (or nu
C          form if ifnuform is .true.).
          rpower = 1.d0
          do iextrasum = 1,nextrasum-2
            extrasum(iextrasum) = extrasum(iextrasum) +
     &        g*extrasum_scale(iextrasum)*rpower
            if(ifnr03) then
              extrasumf(iextrasum) = extrasumf(iextrasum) +
     &          gf*extrasum_scale(iextrasum)*rpower
              extrasumt(iextrasum) = extrasumt(iextrasum) +
     &          gt*extrasum_scale(iextrasum)*rpower
            endif
            if(ifnr13) then
              do index = 1, max_index
                extrasum_dv(index,iextrasum) = 
     &            extrasum_dv(index,iextrasum) +
     &            g_dv(index)*extrasum_scale(iextrasum)*rpower
              enddo
            endif
            rpower = rpower*r_neutral(1)
          enddo
C          sum over ions with weight of Z^1.5 which is unity for H+.
          extrasum(nextrasum-1) = extrasum(nextrasum-1) +
     &      h*extrasum_scale(nextrasum-1)
C          sum over all neutral and ionized states except bare nuclei
          extrasum(nextrasum) = extrasum(nextrasum) +
     &      g*extrasum_scale(nextrasum)*r_ion3(1)
          if(ifnr03) then
            extrasumf(nextrasum-1) = extrasumf(nextrasum-1) +
     &        hf*extrasum_scale(nextrasum-1)
            extrasumt(nextrasum-1) = extrasumt(nextrasum-1) +
     &        ht*extrasum_scale(nextrasum-1)
            extrasumf(nextrasum) =
     &        extrasumf(nextrasum) +
     &        gf*extrasum_scale(nextrasum)*r_ion3(1)
            extrasumt(nextrasum) =
     &        extrasumt(nextrasum) +
     &        gt*extrasum_scale(nextrasum)*r_ion3(1)
          endif
          if(ifnr13) then
            do index = 1, max_index
              extrasum_dv(index,nextrasum-1) = 
     &          extrasum_dv(index,nextrasum-1) +
     &          h_dv(index)*extrasum_scale(nextrasum-1)
              extrasum_dv(index,nextrasum) = 
     &          extrasum_dv(index,nextrasum) +
     &          g_dv(index)*extrasum_scale(nextrasum)*r_ion3(1)
            enddo
          endif
        endif
C        h2plus_dv not needed.
      else
C        partially ionized including hydrogen monatomics and molecules.
        if(tl.ne.tlold) then
          tlold = tl
          call molecular_hydrogen(ifh2, ifh2plus, tl,
     &      qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt)
        endif
        if(ifnr13) then
C          required by hequil
          indexh = inv_ion(1)
C          required by h2equil *and* h2plusequil
          indexh2 = inv_ion(nions+1)
C          required by h2plusequil
          indexh2plus = inv_ion(nions+2)
        endif
C        n(H+)/n(H) = exp(hequil), where hequil  returned from ionize
C        n(H2)*avogadro/(n(H)*n(H)) = exp(h2equil)
        if(ifnr13) then
          hequil_dv = 1.d0
        endif
C        n.b. keep h2equil in ln form for under/over flow reasons.
        h2equil = logalpha - log(4.d0) + qh2 -
     &    1.5d0*tl + h2diss*tc2 + h2_log + dv(nions+1)
        if(ifnr03) then
C          d h2equil/ d ln f
          h2equilf = dvf(nions+1)
C          d h2equil/ d ln T (excluding dv effects)
          h2equilt0 = qh2t - 1.5d0 - h2diss*tc2
C          d h2equil/ d ln T
          h2equilt = h2equilt0 + dvt(nions+1)
        endif
        if(ifnr13) then
          h2equil_dv = 1.d0
        endif
        if(ifh2plus.gt.0) then
C          n(H2+)*avogadro/(n(H)*n(H)) = exp(h2plusequil)
          h2plusequil = h2equil + qh2plus - qh2 -
     &      bi(nions+1)*tc2 + dv(nions+2)
          if(ifnr03) then
            h2plusequilf = h2equilf + dvf(nions+2)
            h2plusequilt = h2equilt +
     &        qh2plust - qh2t + bi(nions+1)*tc2 + dvt(nions+2)
          endif
          if(ifnr13) then
C            this derivative for *both* indexh2 and indexh2plus
            h2plusequil_dv = 1.d0
          endif
        endif
C        solve charge equation which is quadratic in x = n(H)/avogadro,
C        i.e., aa*x^2 + ab*x = ac
C        this equation derived from
C        n(H2+)/avogadro + n(H+)/avogadro + hne*rho = n_e/avogadro
C          = cd*re,
C        where
C        hne*rho is the sum of non-hydrogen positive charges per
C          unit volume divided by avogadro.
C        rho = (2 n(H2) + 2 n(H2+) + n(H)(1+e))/(eps(h)*avogadro)
C        aa = exp(h2plusequil)*(1.d0 + 2.d0*hne/eps(1)) +
C          2.d0*hne*exp(h2equil)/eps(1)
C        ab = exp(hequil)*(1.d0 + hne/eps(1)) + hne/eps(1)
C        rho/mu_e = n_e H = cd*re
C        ac = cd*re
C        must_use_care when scaling this equation.  aa, ab, and ac
C        are all positive.  transform to equation in
C        xprime = x/xmax <= 1, where xmax = ac/ab and
C        x = xmax is the solution if no molecular formation (i.e. aa=0).
C        ab = exp(ablog),
C        ablog = log(exp(abh) + exp(abg)), where
C        abh = hequil + log(1.d0 + hne/eps(1))
C        abg = log(hne/eps(1))
C        to avoid significance loss for large hequil define "d" analogs
C        of ab, abh, and abg where exp(hequil) is divided out.
C        abd = (1.d0 + hne/eps(1)) + hne/eps(1)/exp(hequil)
C        abdlog = log(exp(abdh) + exp(abdg))
C        abdh = log(1.d0 + hne/eps(1))
C        abdg = log(hne/eps(1)) - hequil
        
C        form ln quantities first:
        abh = hequil + log(1.d0 + hne/eps(1))
        abdh = log(1.d0 + hne/eps(1))
        if(hne.gt.0.d0) then
          abg = log(hne/eps(1))
          abdg = log(hne/eps(1)) - hequil
        endif
C        exp(-46) ~ 1.d-20
        ifabh = hne.le.0.d0.or.
     &    (hne.gt.0.d0.and.abg.lt.abh-46.d0)
        ifabg = hne.gt.0.d0.and.abh.lt.abg-46.d0
        if(ifabh) then
          ablog = abh
          abdlog = abdh
        elseif(ifabg) then
          ablog = abg
          abdlog = abdg
        elseif(abh.gt.abg) then
          abexp = exp(abg-abh)
          ablog = abh + log(1.d0 + abexp)
          abdlog = abdh + log(1.d0 + abexp)
        else
          abexp = exp(abh-abg)
          ablog = abg + log(1.d0 + abexp)
          abdlog = abdg + log(1.d0 + abexp)
        endif
        ac = cd*re
        if(ifnr03) then
C          abh = hequil + log(1.d0 + hne/eps(1))
          abhf = hequilf + (hnef/eps(1))/(1.d0 + hne/eps(1))
          abht = hequilt + (hnet/eps(1))/(1.d0 + hne/eps(1))
          abdhf = (hnef/eps(1))/(1.d0 + hne/eps(1))
          abdht = (hnet/eps(1))/(1.d0 + hne/eps(1))
          if(hne.gt.0.d0) then
C            abg = log(hne/eps(1))
            abgf = hnef/hne
            abgt = hnet/hne
            abdgf = hnef/hne - hequilf
            abdgt = hnet/hne - hequilt
          endif
          if(ifabh) then
C            ablog = abh
            ablogf = abhf
            ablogt = abht
            abdlogf = abdhf
            abdlogt = abdht
          elseif(ifabg) then
C            ablog = abg
            ablogf = abgf
            ablogt = abgt
            abdlogf = abdgf
            abdlogt = abdgt
          elseif(abh.gt.abg) then
C            abexp = exp(abg-abh)
C            ablog = abh + log(1.d0 + abexp)
            ablogf = abhf + (abgf-abhf)*abexp/(1.d0 + abexp)
            ablogt = abht + (abgt-abht)*abexp/(1.d0 + abexp)
            abdlogf = abdhf + (abgf-abhf)*abexp/(1.d0 + abexp)
            abdlogt = abdht + (abgt-abht)*abexp/(1.d0 + abexp)
          else
C            abexp = exp(abh-abg)
C            ablog = abg + log(1.d0 + abexp)
            ablogf = abgf + (abhf-abgf)*abexp/(1.d0 + abexp)
            ablogt = abgt + (abht-abgt)*abexp/(1.d0 + abexp)
            abdlogf = abdgf + (abhf-abgf)*abexp/(1.d0 + abexp)
            abdlogt = abdgt + (abht-abgt)*abexp/(1.d0 + abexp)
          endif
          acf = ac*ref
          act = ac*ret
        endif
        if(ifnr13) then
          do index = 1,max_index
C            abh = hequil + log(1.d0 + hne/eps(1))
            abh_dv = (hne_dv(index)/eps(1))/(1.d0 + hne/eps(1))
            abdh_dv = abh_dv
            if(index.eq.indexh) abh_dv = abh_dv + hequil_dv
            if(hne.gt.0.d0) then
C              abg = log(hne/eps(1))
              abg_dv = hne_dv(index)/hne
              abdg_dv = abg_dv
              if(index.eq.indexh) abdg_dv = abdg_dv - hequil_dv
            endif
            if(ifabh) then
C              ablog = abh
              ablog_dv(index) = abh_dv
              abdlog_dv(index) = abdh_dv
            elseif(ifabg) then
C              ablog = abg
              ablog_dv(index) = abg_dv
              abdlog_dv(index) = abdg_dv
            elseif(abh.gt.abg) then
C              abexp = exp(abg-abh)
C              ablog = abh + log(1.d0 + abexp)
              ablog_dv(index) = abh_dv + (abg_dv-abh_dv)*abexp/
     &          (1.d0 + abexp)
              abdlog_dv(index) = abdh_dv + (abg_dv-abh_dv)*abexp/
     &          (1.d0 + abexp)
            else
C              abexp = exp(abh-abg)
C              ablog = abg + log(1.d0 + abexp)
              ablog_dv(index) = abg_dv + (abh_dv-abg_dv)*abexp/
     &          (1.d0 + abexp)
              abdlog_dv(index) = abdg_dv + (abh_dv-abg_dv)*abexp/
     &          (1.d0 + abexp)
            endif
          enddo
        endif
        if(ifh2plus.gt.0.or.hne.gt.0.d0) then
C          the equation to solve for xprimelog = log(xprime) is
C          cprime*xprime^2 + xprime = 1.d0,
C          where cprime = (aa/ac)*xmax^2 = aa*ac/(ab*ab)
C          ab = exp(ablog),
C          and aa = exp(aalog),
C          aalog = log(exp(aah2plus) + exp(aah2)), where
C          aah2plus = h2plusequil + ln(1.d0 + 2.d0*hne/eps(1))
C          aah2 = h2equil + ln(2.d0*hne/eps(1))
C          form ln quantities first:
          if(ifh2plus.gt.0)
     &      aah2plus = h2plusequil + log(1.d0 + 2.d0*hne/eps(1))
          if(hne.gt.0.d0)
     &      aah2 = h2equil + log(2.d0*hne/eps(1))
C          exp(-46) ~ 1.d-20
          ifaah2plus = hne.le.0.d0.or.
     &      (ifh2plus.gt.0.and.hne.gt.0.d0.and.aah2.lt.aah2plus-46.d0)
          ifaah2 = ifh2plus.le.0.or.
     &      (ifh2plus.gt.0.and.hne.gt.0.d0.and.aah2plus.lt.aah2-46.d0)
          if(ifaah2plus) then
            aalog = aah2plus
          elseif(ifaah2) then
            aalog = aah2
          elseif(aah2plus.gt.aah2) then
            aaexp = exp(aah2-aah2plus)
            aalog = aah2plus + log(1.d0 + aaexp)
          else
            aaexp = exp(aah2plus-aah2)
            aalog = aah2 + log(1.d0 + aaexp)
          endif
C          form cprime in ln form to start
          cprime = aalog + log(ac) - 2.d0*ablog
          if(cprime.le.-575.d0) then
C            limit corresponds to about 1.d-250
C            cprime = 0 implies xprime = 1.
            xprimelog = 0.d0
            if(ifnr03) then
              xprimelogf = 0.d0
              xprimelogt = 0.d0
            endif
            if(ifnr13) then
              do index = 1,max_index
                xprimelog_dv(index) = 0.d0
              enddo
            endif
          elseif(cprime.le.92.d0) then
C            limit corresponds to about 1.d40
            cprime = exp(cprime)
C            solve cprime*xprime^2 + xprime = 1.d0
C            for xprimelog = log(xprime) and derivatives
C            temporary_use_which is immediately superseded.
            xprimelog = sqrt(1.d0 + 4.d0*cprime)
C            for now xprimelog actually carries xprime.
            xprimelog = 2.d0/(1.d0 + xprimelog)
            if(ifnr03) then
              if(ifh2plus.gt.0) then
C                aah2plus = h2plusequil + log(1.d0 + 2.d0*hne/eps(1))
                aah2plusf = h2plusequilf +
     &            (2.d0*hnef/eps(1))/(1.d0 + 2.d0*hne/eps(1))
                aah2plust = h2plusequilt +
     &            (2.d0*hnet/eps(1))/(1.d0 + 2.d0*hne/eps(1))
              endif
              if(hne.gt.0.d0) then
C                aah2 = h2equil + log(2.d0*hne/eps(1))
                aah2f = h2equilf + hnef/hne
                aah2t = h2equilt + hnet/hne
              endif
              if(ifaah2plus) then
C                aalog = aah2plus
                aalogf = aah2plusf
                aalogt = aah2plust
              elseif(ifaah2) then
C                aalog = aah2
                aalogf = aah2f
                aalogt = aah2t
              elseif(aah2plus.gt.aah2) then
C                aaexp = exp(aah2-aah2plus)
C                aalog = aah2plus + log(1.d0 + aaexp)
                aalogf = aah2plusf + (aah2f-aah2plusf)*
     &            (aaexp/(1.d0+aaexp))
                aalogt = aah2plust + (aah2t-aah2plust)*
     &            (aaexp/(1.d0+aaexp))
              else
C                aaexp = exp(aah2plus-aah2)
C                aalog = aah2 + log(1.d0 + aaexp)
                aalogf = aah2f + (aah2plusf-aah2f)*
     &            (aaexp/(1.d0+aaexp))
                aalogt = aah2t + (aah2plust-aah2t)*
     &            (aaexp/(1.d0+aaexp))
              endif
C              cprime = exp(aalog + log(ac) - 2.d0*ablog)
              cprimef = aalogf + acf/ac - 2.d0*ablogf
              cprimet = aalogt + act/ac - 2.d0*ablogt
              cprimef = cprime*cprimef
              cprimet = cprime*cprimet
C              implicit differentiation of cprime*xprime^2 + xprime = 1.d0
C              xprimelog currently carries xprime,
C              but derivatives are of the log variable.
              xprimelogf = -cprimef*xprimelog/
     &          (2.d0*cprime*xprimelog + 1.d0)
              xprimelogt = -cprimet*xprimelog/
     &          (2.d0*cprime*xprimelog + 1.d0)
            endif
            if(ifnr13) then
              do index = 1,max_index
                if(ifh2plus.gt.0) then
C                  aah2plus = h2plusequil + log(1.d0 + 2.d0*hne/eps(1))
                  aah2plus_dv = 
     &              (2.d0*hne_dv(index)/eps(1))/
     &              (1.d0 + 2.d0*hne/eps(1))
                  if(index.eq.indexh2.or.index.eq.indexh2plus)
     &              aah2plus_dv = aah2plus_dv + h2plusequil_dv
                endif
                if(hne.gt.0.d0) then
C                  aah2 = h2equil + log(2.d0*hne/eps(1))
                  aah2_dv = hne_dv(index)/hne
                  if(index.eq.indexh2)
     &              aah2_dv = aah2_dv + h2equil_dv
                endif
                if(ifaah2plus) then
C                  aalog = aah2plus
                  aalog_dv = aah2plus_dv
                elseif(ifaah2) then
C                  aalog = aah2
                  aalog_dv = aah2_dv
                elseif(aah2plus.gt.aah2) then
C                  aaexp = exp(aah2-aah2plus)
C                  aalog = aah2plus + log(1.d0 + aaexp)
                  aalog_dv = aah2plus_dv +
     &              (aah2_dv-aah2plus_dv)*(aaexp/(1.d0+aaexp))
                else
C                  aaexp = exp(aah2plus-aah2)
C                  aalog = aah2 + log(1.d0 + aaexp)
                  aalog_dv = aah2_dv +
     &              (aah2plus_dv-aah2_dv)*(aaexp/(1.d0+aaexp))
                endif
C                cprime = exp(aalog + log(ac) - 2.d0*ablog)
                cprime_dv = aalog_dv - 2.d0*ablog_dv(index)
                cprime_dv = cprime*cprime_dv
C                implicit differentiation of cprime*xprime^2 + xprime = 1.d0
C                xprimelog currently carries xprime,
C                but derivatives are of the log variable.
                xprimelog_dv(index) = -cprime_dv*xprimelog/
     &            (2.d0*cprime*xprimelog + 1.d0)
              enddo
            endif
            xprimelog = log(xprimelog)
          else
C            cprime is greater than about 1.d40.
C            solve cprime*xprime^2 = 1, but variable cprime
C            actually carries log of that quantity so
C            really solve: cprime + 2.d0*xprimelog = 0
            xprimelog = -0.5d0*cprime
            if(ifnr03) then
              if(ifh2plus.gt.0) then
C                aah2plus = h2plusequil + log(1.d0 + 2.d0*hne/eps(1))
                aah2plusf = h2plusequilf +
     &            (2.d0*hnef/eps(1))/(1.d0 + 2.d0*hne/eps(1))
                aah2plust = h2plusequilt +
     &            (2.d0*hnet/eps(1))/(1.d0 + 2.d0*hne/eps(1))
              endif
              if(hne.gt.0.d0) then
C                aah2 = h2equil + log(2.d0*hne/eps(1))
                aah2f = h2equilf + hnef/hne
                aah2t = h2equilt + hnet/hne
              endif
              if(ifaah2plus) then
C                aalog = aah2plus
                aalogf = aah2plusf
                aalogt = aah2plust
              elseif(ifaah2) then
C                aalog = aah2
                aalogf = aah2f
                aalogt = aah2t
              elseif(aah2plus.gt.aah2) then
C                aaexp = exp(aah2-aah2plus)
C                aalog = aah2plus + log(1.d0 + aaexp)
                aalogf = aah2plusf + (aah2f-aah2plusf)*
     &            (aaexp/(1.d0+aaexp))
                aalogt = aah2plust + (aah2t-aah2plust)*
     &            (aaexp/(1.d0+aaexp))
              else
C                aaexp = exp(aah2plus-aah2)
C                aalog = aah2 + log(1.d0 + aaexp)
                aalogf = aah2f + (aah2plusf-aah2f)*
     &            (aaexp/(1.d0+aaexp))
                aalogt = aah2t + (aah2plust-aah2t)*
     &            (aaexp/(1.d0+aaexp))
              endif
C              cprime = aalog + log(ac) - 2.d0*ablog
              cprimef = aalogf + acf/ac - 2.d0*ablogf
              cprimet = aalogt + act/ac - 2.d0*ablogt
              xprimelogf = -0.5d0*cprimef
              xprimelogt = -0.5d0*cprimet
            endif
            if(ifnr13) then
              do index = 1,max_index
                if(ifh2plus.gt.0) then
C                  aah2plus = h2plusequil + log(1.d0 + 2.d0*hne/eps(1))
                  aah2plus_dv = 
     &              (2.d0*hne_dv(index)/eps(1))/
     &              (1.d0 + 2.d0*hne/eps(1))
                  if(index.eq.indexh2.or.index.eq.indexh2plus)
     &              aah2plus_dv = aah2plus_dv + h2plusequil_dv
                endif
                if(hne.gt.0.d0) then
C                  aah2 = h2equil + log(2.d0*hne/eps(1))
                  aah2_dv = hne_dv(index)/hne
                  if(index.eq.indexh2)
     &              aah2_dv = aah2_dv + h2equil_dv
                endif
                if(ifaah2plus) then
C                  aalog = aah2plus
                  aalog_dv = aah2plus_dv
                elseif(ifaah2) then
C                  aalog = aah2
                  aalog_dv = aah2_dv
                elseif(aah2plus.gt.aah2) then
C                  aaexp = exp(aah2-aah2plus)
C                  aalog = aah2plus + log(1.d0 + aaexp)
                  aalog_dv = aah2plus_dv +
     &              (aah2_dv-aah2plus_dv)*(aaexp/(1.d0+aaexp))
                else
C                  aaexp = exp(aah2plus-aah2)
C                  aalog = aah2 + log(1.d0 + aaexp)
                  aalog_dv = aah2_dv +
     &              (aah2plus_dv-aah2_dv)*(aaexp/(1.d0+aaexp))
                endif
C                cprime = aalog + log(ac) - 2.d0*ablog
                cprime_dv = aalog_dv - 2.d0*ablog_dv(index)
                xprimelog_dv(index) = -0.5d0*cprime_dv
              enddo
            endif
          endif
        else
C          cprime = 0 implies xprime = 1.
          xprimelog = 0.d0
          if(ifnr03) then
            xprimelogf = 0.d0
            xprimelogt = 0.d0
          endif
          if(ifnr13) then
            do index = 1,max_index
              xprimelog_dv(index) = 0.d0
            enddo
          endif
        endif
C        at this point have calculated xprimelog = log(x/xmax)
C        and derivatives.
C        transform to x  = log(xprime*xmax) = log(xprime*ac/ab) where
C        x is log(n(neutral H)/avogadro).
        x = xprimelog + log(ac) - ablog
        if(ifnr03) then
          xf = xprimelogf + acf/ac - ablogf
          xt = xprimelogt + act/ac - ablogt
        endif
        if(ifnr13) then
          do index = 1,max_index
            x_dv(index) = xprimelog_dv(index) - ablog_dv(index)
          enddo
        endif
C        h2 is n(H2)/avogadro
C        first in ln form:
        h2 = 2.d0*x + h2equil
        if(ifnr03) then
          h2f = 2.d0*xf + h2equilf
          h2t = 2.d0*xt + h2equilt
        endif
        if(ifnr13) then
          do index = 1,max_index
            h2_dv(index) = 2.d0*x_dv(index)
            if(index.eq.indexh2) h2_dv(index) = h2_dv(index) +
     &        h2equil_dv
          enddo
        endif
        if(ifh2plus.gt.0) then
C          h2plus is n(H2+)/avogadro
C          first in ln form:
          h2plus = 2.d0*x + h2plusequil
          if(ifnr03) then
            h2plusf = 2.d0*xf + h2plusequilf
            h2plust = 2.d0*xt + h2plusequilt
          endif
          if(ifnr13) then
            do index = 1,max_index
              h2plus_dv(index) = 2.d0*x_dv(index)
              if(index.eq.indexh2.or.index.eq.indexh2plus)
     &          h2plus_dv(index) = h2plus_dv(index) + h2plusequil_dv
            enddo
          endif
        endif
C        h is n(H+)/avogadro
C        first in ln form:
        if(hequil.gt.10.d0) then
C          large pressure ionization can force hequil to large values.
C          avoid significance loss for this case by using the abd form
C          of variables.
C          x = xprimelog + log(ac) - ablog
C          h = x + hequil
          h = xprimelog + log(ac) - abdlog
          if(ifnr03) then
            hf = xprimelogf + acf/ac - abdlogf
            ht = xprimelogt + act/ac - abdlogt
          endif
          if(ifnr13) then
            do index = 1,max_index
              h_dv(index) = xprimelog_dv(index) - abdlog_dv(index)
            enddo
          endif
        else
          h = x + hequil
          if(ifnr03) then
            hf = xf + hequilf
            ht = xt + hequilt
          endif
          if(ifnr13) then
            do index = 1,max_index
              h_dv(index) = x_dv(index)
              if(index.eq.indexh) h_dv(index) = h_dv(index) +
     &          hequil_dv
            enddo
          endif
        endif
C        N.B. x is already ln(n(H)/avogadro)

C        Find whether any of h2, h2plus, h, or x is dominant so can calculate
C        rl without fear of over/underflow where
C        rl = log(2.d0*exp(h2) + 2.d0*exp(h2plus) + exp(h) + exp(x)) -
C     &    log(eps(1))
C        "dominant" means one of x, h, h2, or h2plus is lnerrcrit
C        larger than any of the rest so that the above log
C        factor can be approximated by the maximum component with error
C        less than ln(1 + exp(-lnerrcrit)) ~ exp(-lnerrcrit).
        maxh = h2
        index_maxh = 1
        second_maxh = maxh - lnerrcrit - 1.d0
        if(ifh2plus.gt.0) then
          if(h2plus.gt.maxh) then
            index_maxh = 2
            second_maxh = maxh
            maxh = h2plus
          else
            second_maxh = max(second_maxh, h2plus)
          endif
        endif
        if(h.gt.maxh) then
          index_maxh = 3
          second_maxh = maxh
          maxh = h
        else
          second_maxh = max(second_maxh, h)
        endif
        if(x.gt.maxh) then
          index_maxh = 4
          second_maxh = maxh
          maxh = x
        else
          second_maxh = max(second_maxh, x)
        endif
        if(maxh - second_maxh.lt.lnerrcrit) index_maxh = 0
          
        if(index_maxh.eq.1) then
          rl = maxh - log(eps(1)/2.d0)
          if(ifnr03) then
            rf = h2f
            rt = h2t
          endif
          if(ifnr13) then
            do index = 1,max_index
              r_dv(index) = h2_dv(index)
            enddo
          endif
        elseif(index_maxh.eq.2) then
          if(ifh2plus.le.0) stop 'eos_calc: internal ifh2plus error'
          rl = maxh - log(eps(1)/2.d0)
          if(ifnr03) then
            rf = h2plusf
            rt = h2plust
          endif
          if(ifnr13) then
            do index = 1,max_index
              r_dv(index) = h2plus_dv(index)
            enddo
          endif
        elseif(index_maxh.eq.3) then
          rl = maxh - log(eps(1))
          if(ifnr03) then
            rf = hf
            rt = ht
          endif
          if(ifnr13) then
            do index = 1,max_index
              r_dv(index) = h_dv(index)
            enddo
          endif
        elseif(index_maxh.eq.4) then
          rl = maxh - log(eps(1))
          if(ifnr03) then
            rf = xf
            rt = xt
          endif
          if(ifnr13) then
            do index = 1,max_index
              r_dv(index) = x_dv(index)
            enddo
          endif
        endif

C        transform h2 (-575.d0 corresponds to underflow limit of 1.d-250)
C        from ln n(H2)/avogadro) form to n(H2)
        if(h2.gt.-575.d0) then
          h2 = avogadro*exp(h2)
          if(ifnr03) then
            h2f = h2*h2f
            h2t = h2*h2t
          endif
          if(ifnr13) then
            do index = 1,max_index
              h2_dv(index) = h2*h2_dv(index)
            enddo
          endif
        else
          h2 = 0.d0
          if(ifnr03) then
            h2f = 0.d0
            h2t = 0.d0
          endif
          if(ifnr13) then
            do index = 1,max_index
              h2_dv(index) = 0.d0
            enddo
          endif
        endif
C        transform h2plus (-575.d0 corresponds to underflow limit of 1.d-250)
C        from ln n(H2+)/avogadro) form to n(H2+)
        if(ifh2plus.gt.0.and.h2plus.gt.-575.d0) then
          h2plus = avogadro*exp(h2plus)
          if(ifnr03) then
            h2plusf = h2plus*h2plusf
            h2plust = h2plus*h2plust
          endif
          if(ifnr13) then
            do index = 1,max_index
              h2plus_dv(index) = h2plus*h2plus_dv(index)
            enddo
          endif
        else
          h2plus = 0.d0
          if(ifnr03) then
            h2plusf = 0.d0
            h2plust = 0.d0
          endif
          if(ifnr13) then
            do index = 1,max_index
              h2plus_dv(index) = 0.d0
            enddo
          endif
        endif
C        transform h (-575.d0 corresponds to underflow limit of 1.d-250)
C        from ln n(H+)/avogadro) form to n(H+)
        if(h.gt.-575.d0) then
          h = avogadro*exp(h)
          if(ifnr03) then
            hf = h*hf
            ht = h*ht
          endif
          if(ifnr13) then
            do index = 1,max_index
              h_dv(index) = h*h_dv(index)
            enddo
          endif
        else
          h = 0.d0
          if(ifnr03) then
            hf = 0.d0
            ht = 0.d0
          endif
          if(ifnr13) then
            do index = 1,max_index
              h_dv(index) = 0.d0
            enddo
          endif
        endif
C        transform x (-575.d0 corresponds to underflow limit of 1.d-250)
C        from ln n(H)/avogadro) form to g = n(H)
C        n.b. x used later so don't mess with it.
        if(x.gt.-575.d0) then
          g = avogadro*exp(x)
          if(ifnr03) then
            gf = g*xf
            gt = g*xt
          endif
          if(ifnr13) then
            do index = 1,max_index
              g_dv(index) = g*x_dv(index)
            enddo
          endif
        else
          g = 0.d0
          if(ifnr03) then
            gf = 0.d0
            gt = 0.d0
          endif
          if(ifnr13) then
            do index = 1,max_index
              g_dv(index) = 0.d0
            enddo
          endif
        endif
        if(index_maxh.eq.0) then
C          calculate rho*eps(1)*avogadro from first abundance constraint 
C          equation.
          rho = g + h + 2.d0*h2 + 2.d0*h2plus
          if(ifnr03) then
C            d ln rho(fl, tl, dv(fl, tl)) wrt fl, tl
            rf = (gf+hf+2.d0*h2f+2.d0*h2plusf)/rho
            rt = (gt+ht+2.d0*h2t+2.d0*h2plust)/rho
          endif
          if(ifnr13) then
C            d ln rho(fl, tl, dv) wrt dv
            do index = 1,max_index
              r_dv(index) = (g_dv(index) + h_dv(index) +
     &          2.d0*h2_dv(index)+2.d0*h2plus_dv(index))/rho
            enddo
          endif
C          transform to actual rho from scaled rho (doesn't affect
C          ln derivatives).
          rho = rho/(eps(1)*avogadro)
          rl = log(rho)
        else
          rho = exp(rl)
        endif
C        calculate nu = n/rho/avogadro forms
C        including ne = number density of electrons/(rho*avogadro)
        nug = g/(rho*avogadro)
C        n.b. x = log(g/avogadro)
        lognug = x - rl
        nuh = h/(rho*avogadro)
        nuh2 = h2/(rho*avogadro)
        nuh2plus = h2plus/(rho*avogadro)
        nuvar(1,1) = nug
        nuvar(2,1) = nuh
C        store nuh2 and nuh2plus in the 3rd and 4th index
        nuvar(3,1) = nuh2
        nuvar(4,1) = nuh2plus
        ne = hne + nuh + nuh2plus
        if(ifnr03) then
C        n.b. these derivatives may have some significance loss 
C        (2 nuh2plus + 2 nuh2 + nug + nuh = eps(1)), but
C        so far haven't noticed any problem
          nugf = gf/(rho*avogadro) - nug*rf
          nugt = gt/(rho*avogadro) - nug*rt
          nuhf = hf/(rho*avogadro) - nuh*rf
          nuht = ht/(rho*avogadro) - nuh*rt
          nuh2f = h2f/(rho*avogadro) - nuh2*rf
          nuh2t = h2t/(rho*avogadro) - nuh2*rt
          nuh2plusf = h2plusf/(rho*avogadro) - nuh2plus*rf
          nuh2plust = h2plust/(rho*avogadro) - nuh2plus*rt
          nef = hnef + nuhf + nuh2plusf
          net = hnet + nuht + nuh2plust
        endif
        if(ifnr13.and.(ifpl.eq.1.or.ifexcited.gt.0)) then
C          these quantities only needed for Planck-Larkin or ifexcited case
          do index = 1,max_index
            nug_dv(index) = g_dv(index)/(rho*avogadro) -
     &        nug*r_dv(index)
            nuh2_dv(index) = h2_dv(index)/(rho*avogadro) -
     &        nuh2*r_dv(index)
            nuh2plus_dv(index) = h2plus_dv(index)/(rho*avogadro) -
     &        nuh2plus*r_dv(index)
          enddo
        endif
        en = 1.d0/ne
C        convert hd, he, and derivatives to number densities.
        if(eps(2).gt.0.d0) then
C          if ifnuform is true then ifnr is zero and leave hd, he
C          and their derivatives in the nu form returned by
C          ionize.
          if(.not.ifnuform) then
            hd = hd*rho*avogadro
            he = he*rho*avogadro
            if(ifpi.eq.2.or.if_mc.eq.1) then
              if(ifnr03) then
                hdf = hdf*rho*avogadro + rf*hd
                hdt = hdt*rho*avogadro + rt*hd
                hef = hef*rho*avogadro + rf*he
                het = het*rho*avogadro + rt*he
              endif
              if(ifnr13) then
                index2 = inv_ion(2)
                index3 = inv_ion(3)
                do index = 1, max_index
                  if(index.eq.index2.or.index.eq.index3) then
                    hd_dv(index) = hd_dv(index)*rho*avogadro +
     &                hd*r_dv(index)
                    he_dv(index) = he_dv(index)*rho*avogadro +
     &                he*r_dv(index)
                  else
                    hd_dv(index) = hd*r_dv(index)
                    he_dv(index) = he*r_dv(index)
                  endif
                enddo
              endif
            endif
          endif
        else
          hd = 0.d0
          he = 0.d0
          if(ifpi.eq.2.or.if_mc.eq.1) then
            if(ifnr03) then
              hdf = 0.d0
              hef = 0.d0
              hdf = 0.d0
              het = 0.d0
            endif
          endif
        endif
C        these expressions derived from vdb notes. they_use_the equations:
C        2 nu(H2+) + 2 nu(H2) + nu(H) + nu(H+) = eps(1),
C        where the fortran variables h2plus, h2, g, and h
C        carry n(H2+), n(H2), n(H), n(H+) or the n/(avagadro*rho)
C        equivalents when ifnuform is true.  It follows that
C        alpha nu(H+)/(A(H+)^3/2 Q(H+)) = 
C          alpha nu(H)/(A(H)^3/2 Q(H)) *
C          exp(-tc2*Hion + dv(1))
C        alpha nu(H2)/(A(H2)^3/2 Q(H2)) = 
C          rho T^{-3/2} [alpha nu(H)/(A(H)^3/2 Q(H))]^2 *
C          exp(tc2*hdiss + dv(nions+1))
C        alpha nu(H2+)/(A(H2+)^3/2 Q(H2+)) = 
C          alpha nu(H2)/(A(H2)^3/2 Q(H2)) *
C          exp(-tc2*H2ion + dv(nions+2))
C        where alpha = avogadro*(2 pi k/(avogadro*h^2)^{-3/2}
C        (see commentary in constants.h), and 
C        Q(i) exp(-tc2*E(i)) A(i)^{3/2}/alpha is the 
C        total partition per unit volume divided by avogadro.
C        Furthermore, we have taken advantage of the fact that
C        the individual rho, T dependence of n(H2+) and n(H2)
C        partially cancels the overall rho, T dependence, see
C        expression for s in awieos_detailed.f.
C        note that h2equilt0 = qh2t - 1.5d0 - h2diss*tc2.  Thus, must
C        subtract 1 from this quantity in sion so that can correct by
C        full_sum0 in awieos.f for the 5/2 term.  The equivalent ln T
C        and ln rho terms are taken care of by the above equilibrium
C        constant relations.
        sharg = bi(1)*tc2 - dv(1)
        sion = sion + constant_sh  - eps(1)*lognug +
     &    nuh*sharg +
     &    nuh2*(h2equilt0-dv(nions+1)-1.d0) +
     &    nuh2plus*(h2equilt0-1.d0+qh2plust-qh2t+bi(nions+1)*tc2 -
     &    dv(nions+1) - dv(nions+2))
        sionf = sionf + nuhf*sharg + nuh2f*(h2equilt0-dv(nions+1)) +
     &    nuh2plusf*(h2equilt0+qh2plust-qh2t+bi(nions+1)*tc2 -
     &    dv(nions+1) - dv(nions+2))
        siont = siont + nuht*sharg + nuh2t*(h2equilt0-dv(nions+1)) +
     &    nuh2*(qh2t+qh2tt) +
     &    nuh2plust*(h2equilt0+qh2plust-qh2t+bi(nions+1)*tc2 -
     &    dv(nions+1) - dv(nions+2)) +
     &    nuh2plus*(qh2plust+qh2plustt)
C        correct ideal internal energy cm^-1 per unit mass
C        per avogadro for hydrogen number densities
        uion = uion + nuh*bi(1) + nuh2*(-h2diss + qh2t/tc2) +
     &    nuh2plus*(bi(nions+1) - h2diss + qh2plust/tc2)
C        add in ideal free-energy terms due to all hydrogen species
C        from first principles (MHD II, F1 + F2 term transformed to
C        free-energy per unit mass and ignoring constant
C        [1 + 1.5d0 ln T - ln rho] term).
        if(nuh2plus.gt.0.d0) then
          fion = fion +
     &      nuh2plus*(log(nuh2plus) + tc2*(bi(nions+1)-h2diss) -
     &      qh2plus - logqtl_const_h2plus)
          if(ifnr03)
     &      fionf = fionf +
     &      nuh2plusf*(log(nuh2plus) + tc2*(bi(nions+1)-h2diss) -
     &      qh2plus - logqtl_const_h2plus + 1.d0)
          if(ifnr13) then
C            n.b., nu form is n/(rho*avogadro) so
C            nu_dv = (n_dv/(rho*avogadro) - nu*r_dv)
            do index = 1, max_index
              fion_dv(index) = fion_dv(index) +
     &          (h2plus_dv(index)/(rho*avogadro) -
     &          nuh2plus*r_dv(index))*
     &          (log(nuh2plus) + tc2*(bi(nions+1)-h2diss) -
     &          qh2plus - logqtl_const_h2plus + 1.d0)
            enddo
          endif
        endif
        if(nuh2.gt.0.d0) then
          fion = fion +
     &      nuh2*(log(nuh2) + tc2*(-h2diss) -
     &      qh2 - logqtl_const_h2)
          if(ifnr03)
     &      fionf = fionf +
     &      nuh2f*(log(nuh2) + tc2*(-h2diss) -
     &      qh2 - logqtl_const_h2 + 1.d0)
          if(ifnr13) then
C            n.b., nu form is n/(rho*avogadro) so
C            nu_dv = (n_dv/(rho*avogadro) - nu*r_dv)
            do index = 1, max_index
              fion_dv(index) = fion_dv(index) +
     &          (h2_dv(index)/(rho*avogadro) - nuh2*r_dv(index))*
     &          (log(nuh2) + tc2*(-h2diss) -
     &          qh2 - logqtl_const_h2 + 1.d0)
            enddo
          endif
        endif
        if(nug.gt.0.d0) then
          fion = fion +
     &      nug*(lognug - logqtl_const_h)
          if(ifnr03)
     &      fionf = fionf +
     &      nugf*(lognug - logqtl_const_h + 1.d0)
          if(ifnr13) then
C            n.b., nu form is n/(rho*avogadro) so
C            nu_dv = (n_dv/(rho*avogadro) - nu*r_dv)
            do index = 1, max_index
              fion_dv(index) = fion_dv(index) +
     &          (g_dv(index)/(rho*avogadro) - nug*r_dv(index))*
     &          (lognug - logqtl_const_h + 1.d0)
            enddo
          endif
        endif
        if(nuh.gt.0.d0) then
          fion = fion +
     &      nuh*(log(nuh) + tc2*bi(1) - logqtl_const_hplus)
          if(ifnr03)
     &      fionf = fionf +
     &      nuhf*(log(nuh) + tc2*bi(1) - logqtl_const_hplus + 1.d0)
          if(ifnr13) then
C            n.b., nu form is n/(rho*avogadro) so
C            nu_dv = (n_dv/(rho*avogadro) - nu*r_dv)
            do index = 1, max_index
              fion_dv(index) = fion_dv(index) +
     &          (h_dv(index)/(rho*avogadro) - nuh*r_dv(index))*
     &          (log(nuh) + tc2*bi(1) - logqtl_const_hplus + 1.d0)
            enddo
          endif
        endif

        if(.not.(if_mc.eq.1.or.if_pteh.eq.1)) then
C          convert Coulomb sums to number densities and add in
C          effect of hydrogen species.
          sum0 = sum0*(rho*avogadro)
          sum2 = sum2*(rho*avogadro)
          if(ifnr03) then
            sum0f = sum0f*(rho*avogadro) + sum0*rf + sum0_scale*
     &        (hf + h2plusf)
            sum0t = sum0t*(rho*avogadro) + sum0*rt + sum0_scale*
     &        (ht + h2plust)
            sum2f = sum2f*(rho*avogadro) + sum2*rf + sum2_scale*
     &        (hf + h2plusf)
            sum2t = sum2t*(rho*avogadro) + sum2*rt + sum2_scale*
     &        (ht + h2plust)
          endif
          if(ifnr13) then
            do index = 1, max_index
              sum0_dv(index) = sum0_dv(index)*(rho*avogadro) +
     &          sum0*r_dv(index) + sum0_scale*
     &          (h_dv(index) + h2plus_dv(index))
              sum2_dv(index) = sum2_dv(index)*(rho*avogadro) +
     &          sum2*r_dv(index) + sum2_scale*
     &          (h_dv(index) + h2plus_dv(index))
            enddo
          endif
C          add in hydrogen species.
          sum0 = sum0 + sum0_scale*(h + h2plus)
          sum2 = sum2 + sum2_scale*(h + h2plus)
        endif
        if(ifpl.eq.1) then
C          correct sumpl0, sumpl1 and sumpl2 for hydrogen nu values.
          sumpl0 = sumpl0 +
     &      nug*plop(1) +
     &      nuh2*plop(nions+1) +
     &      nuh2plus*plop(nions+2)
          sumpl1 = sumpl1 + nug*(plop(1) + plopt(1)) +
     &      nuh2*(plop(nions+1) + plopt(nions+1)) +
     &      nuh2plus*(plop(nions+2) + plopt(nions+2))
          sumpl2 = sumpl2 + nug*plopt(1) +
     &      nuh2*plopt(nions+1) + nuh2plus*plopt(nions+2)
          if(ifnr03) then
            sumpl0f = sumpl0f +
     &        nugf*plop(1) +
     &        nuh2f*plop(nions+1) +
     &        nuh2plusf*plop(nions+2)
            sumpl1f = sumpl1f + nugf*(plop(1) + plopt(1)) +
     &        nuh2f*(plop(nions+1) + plopt(nions+1)) +
     &        nuh2plusf*(plop(nions+2) + plopt(nions+2))
            sumpl1t = sumpl1t + nugt*(plop(1) + plopt(1)) +
     &        nuh2t*(plop(nions+1) + plopt(nions+1)) +
     &        nuh2plust*(plop(nions+2) + plopt(nions+2)) +
     &        nug*(plopt(1) + plopt2(1)) +
     &        nuh2*(plopt(nions+1) + plopt2(nions+1)) +
     &        nuh2plus*(plopt(nions+2) + plopt2(nions+2))
          endif
          if(ifnr13) then
            do index = 1, max_index
              sumpl0_dv(index) = sumpl0_dv(index) +
     &          nug_dv(index)*plop(1) +
     &          nuh2_dv(index)*plop(nions+1) +
     &          nuh2plus_dv(index)*plop(nions+2)
            enddo
          endif
        endif
        if(ifnuform) then
C          replace n by nu form for h2plus, h2, g, h,
C          and derivatives
          h2plus = nuh2plus
          h2 = nuh2
          g = nug
          h = nuh
          h2plusf = nuh2plusf
          h2plust = nuh2plust
          h2f = nuh2f
          h2t = nuh2t
          gf = nugf
          gt = nugt
          hf = nuhf
          ht = nuht
        endif
        if(ifpi34) then
C          if ifnuform is true then leave
C          extrasum and its derivative in the nu form returned by
C          ionize.
          if(.not.ifnuform) then
            do iextrasum = 1, nextrasum
C              convert to n form from nu = n/(rho*avogadro) form
              extrasum(iextrasum) = extrasum(iextrasum)*
     &          (rho*avogadro)
              if(ifnr03) then
                extrasumf(iextrasum) =
     &            extrasumf(iextrasum)*(rho*avogadro) +
     &            extrasum(iextrasum)*rf
                extrasumt(iextrasum) =
     &            extrasumt(iextrasum)*(rho*avogadro) +
     &            extrasum(iextrasum)*rt
              endif
              if(ifnr13) then
                do index = 1, max_index
                  extrasum_dv(index,iextrasum) =
     &              extrasum_dv(index,iextrasum)*(rho*avogadro) +
     &              extrasum(iextrasum)*r_dv(index)
                enddo
              endif
            enddo
          endif
C          correct extrasum for hydrogen number densities (or
C          nu values if ifnuform is .true.).
          rpower = 1.d0
          rpowerh2 = 1.d0
          rpowerh2plus = 1.d0
          do iextrasum = 1,nextrasum-2
            extrasum(iextrasum) = extrasum(iextrasum) +
     &        (h2*rpowerh2 + h2plus*rpowerh2plus + g*rpower)*
     &        extrasum_scale(iextrasum)
            if(ifnr03) then
              extrasumf(iextrasum) = extrasumf(iextrasum) +
     &          (h2f*rpowerh2 + h2plusf*rpowerh2plus + gf*rpower)*
     &          extrasum_scale(iextrasum)
              extrasumt(iextrasum) = extrasumt(iextrasum) +
     &          (h2t*rpowerh2 + h2plust*rpowerh2plus + gt*rpower)*
     &          extrasum_scale(iextrasum)
            endif
            if(ifnr13) then
              do index = 1, max_index
                extrasum_dv(index,iextrasum) =
     &            extrasum_dv(index,iextrasum) +
     &            (h2_dv(index)*rpowerh2 +
     &            h2plus_dv(index)*rpowerh2plus + g_dv(index)*rpower)*
     &            extrasum_scale(iextrasum)
              enddo
            endif
            rpower = rpower*r_neutral(1)
            rpowerh2 = rpowerh2*r_neutral(nelements+1)
            rpowerh2plus = rpowerh2plus*r_neutral(nelements+2)
          enddo
C          sum over ions with weight of Z^1.5 which is unity for H2+ and H+.
          extrasum(nextrasum-1) = extrasum(nextrasum-1) +
     &      (h2plus + h)*extrasum_scale(nextrasum-1)
C          sum over all neutral and ionized states except bare nuclei
          extrasum(nextrasum) = extrasum(nextrasum) +
     &      (h2plus*r_ion3(nions+2) +
     &      h2*r_ion3(nions+1) + g*r_ion3(1))*
     &      extrasum_scale(nextrasum)
          if(ifnr03) then
            extrasumf(nextrasum-1) =
     &        extrasumf(nextrasum-1) +
     &        (h2plusf + hf)*extrasum_scale(nextrasum-1)
            extrasumt(nextrasum-1) =
     &        extrasumt(nextrasum-1) +
     &        (h2plust + ht)*extrasum_scale(nextrasum-1)
            extrasumf(nextrasum) =
     &        extrasumf(nextrasum) +
     &        (h2plusf*r_ion3(nions+2) +
     &        h2f*r_ion3(nions+1) + gf*r_ion3(1))*
     &        extrasum_scale(nextrasum)
            extrasumt(nextrasum) =
     &        extrasumt(nextrasum) +
     &        (h2plust*r_ion3(nions+2) +
     &        h2t*r_ion3(nions+1) + gt*r_ion3(1))*
     &        extrasum_scale(nextrasum)
          endif
          if(ifnr13) then
            do index = 1, max_index
              extrasum_dv(index,nextrasum-1) = 
     &          extrasum_dv(index,nextrasum-1) +
     &          (h2plus_dv(index) + h_dv(index))*
     &          extrasum_scale(nextrasum-1)
              extrasum_dv(index,nextrasum) =
     &          extrasum_dv(index,nextrasum) +
     &          (h2plus_dv(index)*r_ion3(nions+2) +
     &          h2_dv(index)*r_ion3(nions+1) +
     &          g_dv(index)*r_ion3(1))*
     &          extrasum_scale(nextrasum)
            enddo
          endif
        endif
      endif
      if(ifexcited.gt.0) then
        if(ifpi34) then
          if(ifnuform) then
            do jextrasum = 1, 4
              if(jextrasum.lt.4) then
                iextrasum = jextrasum
              else
                iextrasum = nextrasum-1
              endif
C              convert to n form from nu = n/(rho*avogadro) form
              lextrasum(iextrasum) =
     &          (extrasum(iextrasum)/extrasum_scale(iextrasum))*
     &          (rho*avogadro)
              if(ifnr03) then
                lextrasumf(iextrasum) =
     &            (extrasumf(iextrasum)/extrasum_scale(iextrasum))*
     &            (rho*avogadro) +
     &            lextrasum(iextrasum)*rf
                lextrasumt(iextrasum) =
     &            (extrasumt(iextrasum)/extrasum_scale(iextrasum))*
     &            (rho*avogadro) +
     &            lextrasum(iextrasum)*rt
              endif
              if(ifnr13) then
                do index = 1, max_index
                  lextrasum_dv(index,iextrasum) =
     &              (extrasum_dv(index,iextrasum)/
     &              extrasum_scale(iextrasum))*
     &              (rho*avogadro) +
     &              lextrasum(iextrasum)*r_dv(index)
                enddo
              endif
            enddo
          else
            do jextrasum = 1, 4
              if(jextrasum.lt.4) then
                iextrasum = jextrasum
              else
                iextrasum = nextrasum-1
              endif
C              copy n form 
              lextrasum(iextrasum) =
     &          extrasum(iextrasum)/extrasum_scale(iextrasum)
              if(ifnr03) then
                lextrasumf(iextrasum) =
     &            extrasumf(iextrasum)/extrasum_scale(iextrasum)
                lextrasumt(iextrasum) =
     &            extrasumt(iextrasum)/extrasum_scale(iextrasum)
              endif
              if(ifnr13) then
                do index = 1, max_index
                  lextrasum_dv(index,iextrasum) =
     &              extrasum_dv(index,iextrasum)/
     &              extrasum_scale(iextrasum)
                enddo
              endif
            enddo
          endif
        endif
        call excitation_sum(ifexcited, ifsame_zero_abundances,
     &    ifpl, ifpi, ifmodified, ifnr, inv_ion, ifh2, ifh2plus,
     &    partial_elements, n_partial_elements, ion_end,
     &    tl, izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &    bi, plop, plopt, plopt2,
     &    r_ion3, nion, nions, r_neutral, nelements,
     &    lextrasum, lextrasumf, lextrasumt, lextrasum_dv, nextrasum,
     &    max_index, nug, nugf, nugt, nug_dv,
     &    nuh2, nuh2f, nuh2t, nuh2_dv,
     &    nuh2plus, nuh2plusf, nuh2plust, nuh2plus_dv,
     &    xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
        if(.not.ifnuform.and.ifpi34) then
          do iextrasum = 1,4
            xextrasum(iextrasum) = xextrasum(iextrasum)*(rho*avogadro)
            if(ifnr03) then
              xextrasumf(iextrasum) = xextrasumf(iextrasum)*
     &          (rho*avogadro) + xextrasum(iextrasum)*rf
              xextrasumt(iextrasum) = xextrasumt(iextrasum)*
     &          (rho*avogadro) + xextrasum(iextrasum)*rt
            endif
            if(ifnr13) then
              do index = 1, max_index
                xextrasum_dv(index, iextrasum) =
     &            xextrasum_dv(index, iextrasum)*(rho*avogadro) +
     &            xextrasum(iextrasum)*r_dv(index)
              enddo
            endif
          enddo
        endif
      endif
C      unscale all calculated scaled auxiliary variables and their
C      derivatives and do smooth transition to zero for all variables
C      and their derivatives if the magnitude of the auxiliary variable
C      is less than aux_underflow.
      sum0 = sum0/sum0_scale
      if(sum0.lt.aux_underflow) sum0 = 0.d0
      sum2 = sum2/sum2_scale
      if(sum2.lt.aux_underflow) sum2 = 0.d0
      if(ifnr03) then
        if(sum0.gt.0.d0) then
          sum0f = sum0f/sum0_scale
          sum0t = sum0t/sum0_scale
        else
          sum0f = 0.d0
          sum0t = 0.d0
        endif
        if(sum2.gt.0.d0) then
          sum2f = sum2f/sum2_scale
          sum2t = sum2t/sum2_scale
        else
          sum2f = 0.d0
          sum2t = 0.d0
        endif
      endif
      if(ifnr13) then
        do index = 1, max_index
          if(sum0.gt.0.d0) then
            sum0_dv(index) = sum0_dv(index)/sum0_scale
          else
            sum0_dv(index) = 0.d0
          endif
          if(sum2.gt.0.d0) then
            sum2_dv(index) = sum2_dv(index)/sum2_scale
          else
            sum2_dv(index) = 0.d0
          endif
        enddo
      endif
      if(ifpi34) then
        do iextrasum = 1, nextrasum
          extrasum(iextrasum) = extrasum(iextrasum)/
     &      extrasum_scale(iextrasum)
          if(extrasum(iextrasum).lt.aux_underflow)
     &      extrasum(iextrasum) = 0.d0
          if(ifnr03) then
            if(extrasum(iextrasum).gt.0.d0) then
              extrasumf(iextrasum) = extrasumf(iextrasum)/
     &          extrasum_scale(iextrasum)
              extrasumt(iextrasum) = extrasumt(iextrasum)/
     &          extrasum_scale(iextrasum)
            else
              extrasumf(iextrasum) = 0.d0
              extrasumt(iextrasum) = 0.d0
            endif
          endif
          if(ifnr13) then
            if(extrasum(iextrasum).gt.0.d0) then
              do index = 1, max_index
                extrasum_dv(index,iextrasum) =
     &            extrasum_dv(index,iextrasum)/
     &            extrasum_scale(iextrasum)
              enddo
            else
              do index = 1, max_index
                extrasum_dv(index,iextrasum) = 0.d0
              enddo
            endif
          endif
        enddo
        if(ifexcited.gt.0) then
          do iextrasum = 1,4
            xextrasum(iextrasum) = xextrasum(iextrasum)/
     &        xextrasum_scale(iextrasum)
            if(abs(xextrasum(iextrasum)).lt.aux_underflow)
     &        xextrasum(iextrasum) = 0.d0
            if(ifnr03) then
              if(abs(xextrasum(iextrasum)).gt.0.d0) then
                xextrasumf(iextrasum) = xextrasumf(iextrasum)/
     &            xextrasum_scale(iextrasum)
                xextrasumt(iextrasum) = xextrasumt(iextrasum)/
     &            xextrasum_scale(iextrasum)
              else
                xextrasumf(iextrasum) = 0.d0
                xextrasumt(iextrasum) = 0.d0
              endif
            endif
            if(ifnr13) then
              if(abs(xextrasum(iextrasum)).gt.0.d0) then
                do index = 1, max_index
                  xextrasum_dv(index,iextrasum) =
     &              xextrasum_dv(index,iextrasum)/
     &              xextrasum_scale(iextrasum)
                enddo
              else
                do index = 1, max_index
                  xextrasum_dv(index,iextrasum) = 0.d0
                enddo
              endif
            endif
          enddo
        endif
      endif
      end
