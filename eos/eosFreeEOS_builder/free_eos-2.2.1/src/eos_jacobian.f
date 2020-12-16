C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: eos_jacobian.f 842 2008-07-07 21:18:51Z airwin $
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
      subroutine eos_jacobian(
     &  allow_log, allow_log_neg, partial_aux,
     &  ifrad,
     &  debug_dauxdv, itemp1, itemp2, delta1, delta2,
     &  debug_dvdaux, debug_jacobian,
     &  degeneracy, pressure, density, energy, enthalpy, entropy,
     &  iaux1, iaux2,
     &  match_variable, match_variable_save, kif, fl, tl_save,
     &  aux_old, aux, auxf, auxt, daux_dv,
     &  njacobian, njacobian_final,
     &  faux, jacobian, p, pr,
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
     &  pnorad, pnoradf, pnoradt, pnorad_daux,
     &  free, free_daux,
     &  nextrasum, maxnextrasum,
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  h2, h2f, h2t, h2_dv)
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
C         the following are returned only if .not.( if_pteh.eq.1.or.if_mc.eq.1)
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
C     &  dve, dvef, dvet, dve0, dve2, 
C     &  dve_pi, dve_pif, dve_pit, dve_pi_aux(maxfjs_aux),
C     &  dve_coulomb, dve_coulombf, dve_coulombt,
     &  dve_exchange, dve_exchangef, dve_exchanget,
C     &  dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22,
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
     &  ifionized, max_index,
     &  index, maxnextrasum,
     &  nextrasum, ifpi_local, ifpl, nxextrasum,
     &  ifsame_abundances,
     &  ifsame_zero_abundances
      parameter (nxextrasum = 4)
      double precision lambda, gamma_e,
     &  full_sum0, full_sum1, full_sum2,
C     &  sum0_aux(maxfjs_aux+1), sum2_aux(maxfjs_aux+1),
     &  r_ion3(nions+2), r_neutral(nelements+2),
     &  charge(nions+2), charge2(nions+2),
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  tc2
      double precision f, wf, n_e,
     &  sion, sionf, siont, uion,
     &  h2, h2f, h2t, h2_dv(nions+2)
      integer ifdv(nions+2), ifelement(nelements),
     &  partial_ions(nions+2), inv_ion(nions+2),
     &  n_partial_elements, partial_elements(nelements+2)
C      iextraoff is the number of non-extrasum auxiliary variables
      integer iextraoff, naux
      parameter (iextraoff = 7)
      integer 
     &  inv_aux(naux)
      double precision bmin(maxcorecharge)
      logical 
     &  ifsame_under, ifdvzero
      logical allow_log(naux), allow_log_neg(naux),
     &  ifcsum, ifextrasum, ifxextrasum, if_dv
      double precision aux_old(naux)
      integer partial_aux(naux)
      integer index_aux, n_partial_aux
      double precision ddv_aux(nions+2,naux),
C     &  rho, rho_new,
C     &  dvh, dvhf, dvht, dvh_aux(maxfjs_aux),
C     &  dvhe1, dvhe1f, dvhe1t, dvhe1_aux(maxfjs_aux),
C     &  dvhe2, dvhe2f, dvhe2t, dvhe2_aux(maxfjs_aux),
     &  sumion0, sumion2,
C     &  dsum0, dsum2,
     &  sum0_mc, sum2_mc, hcon_mc, hecon_mc,
C     &  nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt,
     &  nux, nuy, nuz
      integer
     &  iaux, jndex_aux,
     &  mion_end, n_partial_ions
C      variables associated with pressure and free energy calculation:
      integer nions_local, naux_local, maxnextrasum_local,
     &  nxextrasum_local
      parameter(nions_local = 295)
      parameter(naux_local = 21)
C      Note, this value must be greater than or equal to 4 as well to store
C      fjs_pi_free and corresponding pressure quantities.
      parameter(maxnextrasum_local = 9)
      parameter(nxextrasum_local = 4)
      double precision
C     &  pe_cgs,
C     &  p0, pion, pionf, piont, pion_dv(nions_local+2),
C     &  pex, pexf, pext,
C     &  ppi, ppif, ppit, ppir, ppi_dv(nions_local+2),
C     &  ppi2_daux(maxnextrasum_local), ppi_daux(naux_local),
C     &  pexcited, pexcitedf, pexcitedt, pexcited_dv(nions_local+2),
C     &  dpcoulomb, dpcoulombf, dpcoulombt, dpcoulomb0, dpcoulomb2,
C     &  pnorad_dv(nions_local+2),
     &  pnorad, pnoradf, pnoradt, pnorad_daux(n_partial_aux+2),
     &  fion, fionf, !fion_daux(naux_local+2),
     &  fion_dv(nions_local+2),! free_dv(nions_local+2),
     &  free, free_daux(n_partial_aux+2),
     &  h, hf, ht, h_dv(nions_local+2),
     &  hd, hdf, hdt, hd_dv(nions_local+2),
     &  he, hef, het, he_dv(nions_local+2),
     &  rl, rt, rf, r_dv(nions_local+2),
     &  h2plus, h2plusf, h2plust, h2plus_dv(nions_local+2),
     &  sum0, sum0f, sum0t, sum0_dv(nions_local+2),
     &  sum2, sum2f, sum2t, sum2_dv(nions_local+2),
     &  sumpl0, sumpl0f, sumpl0_dv(nions_local+2),
     &  extrasum(maxnextrasum_local),
     &  extrasumf(maxnextrasum_local), 
     &  extrasumt(maxnextrasum_local),
     &  extrasum_dv(nions_local+2, maxnextrasum_local),
     &  xextrasum(nxextrasum_local),
     &  xextrasumf(nxextrasum_local),
     &  xextrasumt(nxextrasum_local),
     &  xextrasum_dv(nions_local+2,nxextrasum_local)
C      for debugging
      logical debug_dauxdv, debug_dvdaux, debug_jacobian
      integer itemp1, itemp2,
     &  jtemp1, jtemp2, jtemp3, jtemp4, jtemp5, jtemp6
      double precision delta1, delta2
      double precision degeneracy(3), pressure(3), density(3), energy(3),
     &  enthalpy(3), entropy(3)
      integer iaux1, iaux2, kif, njacobian, njacobian_final
      double precision match_variable, match_variable_save, fl, tl_save,
     &  aux(naux), auxf(naux), auxt(naux), daux_dv(nions+2,naux)
      integer ifrad
      double precision faux(naux), jacobian(naux,naux),
     &  daux_daux(naux_local, naux_local),
     &  p, pr, deriv_lnp_factor
      integer jaux
      logical debug_output
      save
C      sanity checks.
      if(nextrasum.gt.maxnextrasum_local)
     &  stop 'eos_jacobian: maxnextrasum_local too small'
      if(.not.(ifnr.eq.0.or.ifnr.eq.1.or.ifnr.eq.3))
     &  stop 'eos_jacobian: ifnr must be 0, 1, or 3'
      if(ifnr.eq.0) stop 'eos_jacobian: ifnr = 0 is disabled'
      if(nions.ne.nions_local)
     &  stop 'eos_jacobian: nions must be equal to nions_local'
      if(naux.ne.naux_local)
     &  stop 'eos_jacobian: naux must be equal to naux_local'
      if(n_partial_aux.gt.naux_local)
     &  stop 'eos_jacobian: n_partial_aux too large'
      if(nxextrasum.gt.nxextrasum_local)
     &  stop 'eos_jacobian: nxextrasum too large'

      ifcsum = .not.(if_mc.eq.1.or.if_pteh.eq.1)
      ifextrasum = ifpi_local.eq.3.or.ifpi_local.eq.4
      ifxextrasum = ifextrasum .and. (ifexcited.gt.0)
      if_dv = ifnr.eq.1.or.ifnr.eq.3
      
      call aux_to_traditional(
     &  nions, max_index, naux, inv_aux,
     &  aux_old,
     &  h,
     &  hd,
     &  he,
     &  rl,
     &  h2plus,
     &  sum0,
     &  sum2,
     &  iextraoff, maxnextrasum, nxextrasum,
     &  extrasum,
     &  xextrasum)

      call eos_warm_step(
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
     &  ifsame_abundances, ifmtrace,
     &  iatomic_number, ifpi_local,
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
     &  h2, h2f, h2t, h2_dv,
     &  h2plus, h2plusf, h2plust, h2plus_dv,
     &  xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
      call traditional_to_aux(
     &  if_mc, ifcsum, ifextrasum, ifxextrasum, if_dv,
     &  h, hf, ht, h_dv,
     &  hd, hdf, hdt, hd_dv,
     &  he, hef, het, he_dv,
     &  rl, rf, rt, r_dv,
     &  h2plus, h2plusf, h2plust, h2plus_dv,
     &  sum0, sum0f, sum0t, sum0_dv,
     &  sum2, sum2f, sum2t, sum2_dv,
     &  iextraoff, maxnextrasum, nxextrasum,
     &  extrasum, extrasumf, extrasumt, extrasum_dv,
     &  xextrasum, xextrasumf, xextrasumt, xextrasum_dv,
     &  nions, max_index, naux,
     &  aux, auxf, auxt, daux_dv)
      
      if(debug_dvdaux) then
C        partial derivative of dv wrt auxiliary variables
C        indices correspond to H+, He+, He++, C+, and H2 and H2+.
        jtemp1 = 1
        jtemp2 = 2
        jtemp3 = 3
        jtemp4 = 4
        jtemp5 = max_index-1
        jtemp6 = max_index
        degeneracy(1) = dv(partial_ions(jtemp1))
        if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
          degeneracy(2) =
     &      ddv_aux(jtemp1,itemp1)*aux_old(iaux1)
        elseif(itemp1.eq.n_partial_aux+1) then
          degeneracy(2) = dvf(partial_ions(jtemp1))
        elseif(itemp1.eq.n_partial_aux+2) then
          degeneracy(2) = dvt(partial_ions(jtemp1))
        endif
        if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
          degeneracy(3) =
     &      ddv_aux(jtemp1,itemp2)*aux_old(iaux2)
        elseif(itemp2.eq.n_partial_aux+1) then
          degeneracy(3) = dvf(partial_ions(jtemp1))
        elseif(itemp2.eq.n_partial_aux+2) then
          degeneracy(3) = dvt(partial_ions(jtemp1))
        endif
        pressure(1) = dv(partial_ions(jtemp2))
        if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
          pressure(2) =
     &      ddv_aux(jtemp2,itemp1)*aux_old(iaux1)
        elseif(itemp1.eq.n_partial_aux+1) then
          pressure(2) = dvf(partial_ions(jtemp2))
        elseif(itemp1.eq.n_partial_aux+2) then
          pressure(2) = dvt(partial_ions(jtemp2))
        endif
        if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
          pressure(3) =
     &      ddv_aux(jtemp2,itemp2)*aux_old(iaux2)
        elseif(itemp2.eq.n_partial_aux+1) then
          pressure(3) = dvf(partial_ions(jtemp2))
        elseif(itemp2.eq.n_partial_aux+2) then
          pressure(3) = dvt(partial_ions(jtemp2))
        endif
        density(1) = dv(partial_ions(jtemp3))
        if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
          density(2) =
     &      ddv_aux(jtemp3,itemp1)*aux_old(iaux1)
        elseif(itemp1.eq.n_partial_aux+1) then
          density(2) = dvf(partial_ions(jtemp3))
        elseif(itemp1.eq.n_partial_aux+2) then
          density(2) = dvt(partial_ions(jtemp3))
        endif
        if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
          density(3) =
     &      ddv_aux(jtemp3,itemp2)*aux_old(iaux2)
        elseif(itemp2.eq.n_partial_aux+1) then
          density(3) = dvf(partial_ions(jtemp3))
        elseif(itemp2.eq.n_partial_aux+2) then
          density(3) = dvt(partial_ions(jtemp3))
        endif
        energy(1) = dv(partial_ions(jtemp4))
        if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
          energy(2) =
     &      ddv_aux(jtemp4,itemp1)*aux_old(iaux1)
        elseif(itemp1.eq.n_partial_aux+1) then
          energy(2) = dvf(partial_ions(jtemp4))
        elseif(itemp1.eq.n_partial_aux+2) then
          energy(2) = dvt(partial_ions(jtemp4))
        endif
        if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
          energy(3) =
     &      ddv_aux(jtemp4,itemp2)*aux_old(iaux2)
        elseif(itemp2.eq.n_partial_aux+1) then
          energy(3) = dvf(partial_ions(jtemp4))
        elseif(itemp2.eq.n_partial_aux+2) then
          energy(3) = dvt(partial_ions(jtemp4))
        endif
        enthalpy(1) = dv(partial_ions(jtemp5))
        if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
          enthalpy(2) =
     &      ddv_aux(jtemp5,itemp1)*aux_old(iaux1)
        elseif(itemp1.eq.n_partial_aux+1) then
          enthalpy(2) = dvf(partial_ions(jtemp5))
        elseif(itemp1.eq.n_partial_aux+2) then
          enthalpy(2) = dvt(partial_ions(jtemp5))
        endif
        if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
          enthalpy(3) =
     &      ddv_aux(jtemp5,itemp2)*aux_old(iaux2)
        elseif(itemp2.eq.n_partial_aux+1) then
          enthalpy(3) = dvf(partial_ions(jtemp5))
        elseif(itemp2.eq.n_partial_aux+2) then
          enthalpy(3) = dvt(partial_ions(jtemp5))
        endif
        entropy(1) = dv(partial_ions(jtemp6))
        if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
          entropy(2) =
     &      ddv_aux(jtemp6,itemp1)*aux_old(iaux1)
        elseif(itemp1.eq.n_partial_aux+1) then
          entropy(2) = dvf(partial_ions(jtemp6))
        elseif(itemp1.eq.n_partial_aux+2) then
          entropy(2) = dvt(partial_ions(jtemp6))
        endif
        if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
          entropy(3) =
     &      ddv_aux(jtemp6,itemp2)*aux_old(iaux2)
        elseif(itemp2.eq.n_partial_aux+1) then
          entropy(3) = dvf(partial_ions(jtemp6))
        elseif(itemp2.eq.n_partial_aux+2) then
          entropy(3) = dvt(partial_ions(jtemp6))
        endif
        match_variable = match_variable_save
        if(kif.eq.0) fl = match_variable_save
        tl = tl_save
        return
      endif
      if(if_mc.eq.1) then
C        in just this case
C        ionize and eos_calc have only calculated fully
C        ionized part/(rho*avogadro).
C        save it.
        sum0_mc = sum0
        sum2_mc = sum2
      endif
      if(debug_dauxdv) then
C        test of daux_dv derivatives
C        Choose jtemp1 ==> jtemp6 to correspond to auxiliary variables
C        actually used for particular free-energy model being tested.
C        reminder, there are only 15 of these variables for EOS1 in
C        order sum0, sum2, 7 neutral sums, 2 ionized sums, and 4 xextra
C        sums.
        if(.false.) then
          jtemp1 = 1
          jtemp2 = 2
          jtemp3 = 3
          jtemp4 = 4
          jtemp5 = 5
          jtemp6 = 6
        elseif(.false.) then
          jtemp1 = 7
          jtemp2 = 8
          jtemp3 = 9
          jtemp4 = 10
          jtemp5 = 11
          jtemp6 = 12
        elseif(.true.) then
          jtemp1 = 10
          jtemp2 = 11
          jtemp3 = 12
C          prior to this are repeats.
          jtemp4 = 13
          jtemp5 = 14
          jtemp6 = 15
        endif
        degeneracy(1) = aux(partial_aux(jtemp1))
        degeneracy(2) = daux_dv(itemp1,partial_aux(jtemp1))
        degeneracy(3) = daux_dv(itemp2,partial_aux(jtemp1))
        pressure(1) = aux(partial_aux(jtemp2))
        pressure(2) = daux_dv(itemp1,partial_aux(jtemp2))
        pressure(3) = daux_dv(itemp2,partial_aux(jtemp2))
        density(1) = aux(partial_aux(jtemp3))
        density(2) = daux_dv(itemp1,partial_aux(jtemp3))
        density(3) = daux_dv(itemp2,partial_aux(jtemp3))
        energy(1) = aux(partial_aux(jtemp4))
        energy(2) = daux_dv(itemp1,partial_aux(jtemp4))
        energy(3) = daux_dv(itemp2,partial_aux(jtemp4))
        enthalpy(1) = aux(partial_aux(jtemp5))
        enthalpy(2) = daux_dv(itemp1,partial_aux(jtemp5))
        enthalpy(3) = daux_dv(itemp2,partial_aux(jtemp5))
        entropy(1) = aux(partial_aux(jtemp6))
        entropy(2) = daux_dv(itemp1,partial_aux(jtemp6))
        entropy(3) = daux_dv(itemp2,partial_aux(jtemp6))
        return
      endif

C      I previously had logic here to_use_njacobian = n_partial_aux
C      (i.e., the fixed fl solution) for the poorly converged case and
C      only_use_njacobian = njacobian_final (i.e., the combined
C      solution for auxiliary variables and fl) for final convergence.
C      I now always_use_the combined solution subject to scaling of
C      the raw solution vector.
      njacobian = njacobian_final
C      change allow_log and allow_log_neg to refer to both aux and aux_old
      do jndex_aux = 1,n_partial_aux
        jaux = partial_aux(jndex_aux)
C        n.b. 4th auxiliary variable (rl = log(rho)) is already in 
C        log form.
C        N.B. "allow_log" connotes "allow taking log of auxiliary variable
C        delivered by eos_warm_step"
        allow_log(jndex_aux) = jaux.ne.4.and.
     &    aux(jaux).gt.0.d0.and.
     &    aux_old(jaux).gt.0.d0
C        sign of xextrasum depends on extrasum.  it's okay if negative,
C        but subsequent logic must adjust for this.
C        N.B. "allow_log_neg" connotes "allow taking log of *negative* of
C        auxiliary variable delivered by eos_warm_step"
        allow_log_neg(jndex_aux) = jaux.gt.maxnextrasum.and.
     &    aux(jaux).lt.0.d0.and.
     &    aux_old(jaux).lt.0.d0
      enddo
C      form RHS and Jacobian
      do jndex_aux = 1,n_partial_aux
        jaux = partial_aux(jndex_aux)
C        RHS that will be transformed to change in aux from aux_old.
        if(allow_log(jndex_aux)) then
          faux(jndex_aux) = log(aux(jaux))
        elseif(allow_log_neg(jndex_aux)) then
          faux(jndex_aux) = log(-aux(jaux))
        else
          faux(jndex_aux) = aux(jaux)
        endif
        do index_aux = 1,n_partial_aux
C          negative partial aux(jaux) wrt aux_old(iaux)
          jacobian(jndex_aux,index_aux) = 0.d0
C          n.b. ddv_aux(index,index_aux)  has compact index for 
C          index and index_aux
C          n.b. daux_dv(index,jaux)  has compact index for 
C          index and *uncompact*index for jaux
          do index = 1,max_index
            jacobian(jndex_aux,index_aux) = 
     &        jacobian(jndex_aux,index_aux) -
     &        daux_dv(index,jaux)*ddv_aux(index,index_aux)
          enddo
C          partial aux(jaux) wrt aux_old(iaux)
          daux_daux(jndex_aux,index_aux) =
     &      -jacobian(jndex_aux,index_aux)
C          transform to derivative of faux = log(aux) or log(-aux)
          if(allow_log(jndex_aux).or.
     &        allow_log_neg(jndex_aux)) then
            jacobian(jndex_aux,index_aux) =
     &        jacobian(jndex_aux,index_aux)/aux(jaux)
          else
            jacobian(jndex_aux,index_aux) =
     &        jacobian(jndex_aux,index_aux)
          endif
C          transform to derivative with respect to log(auxold) or
C          log(-auxold)
          iaux = partial_aux(index_aux)
          if(allow_log(index_aux).or.
     &        allow_log_neg(index_aux)) then
            jacobian(jndex_aux,index_aux) =
     &        jacobian(jndex_aux,index_aux)*aux_old(iaux)
          endif
        enddo
        if(allow_log(jndex_aux)) then
          faux(jndex_aux) = faux(jndex_aux) - log(aux_old(jaux))
        elseif(allow_log_neg(jndex_aux)) then
          faux(jndex_aux) = faux(jndex_aux) - log(-aux_old(jaux))
        else
          faux(jndex_aux) = faux(jndex_aux) - aux_old(jaux)
        endif
C        take negative derivative wrt to log(auxold), log(-auxold)
C        or auxold whichever is appropriate
        jacobian(jndex_aux,jndex_aux) = 
     &    jacobian(jndex_aux,jndex_aux) + 1.d0
C        negative partial of aux(new) wrt fl with aux(old) fixed.
        if(njacobian.gt.n_partial_aux) then
          jacobian(jndex_aux,njacobian) = -auxf(jaux)
C          transform to derivative of faux = log(aux) or log(-aux)
          if(allow_log(jndex_aux).or.
     &        allow_log_neg(jndex_aux)) then
            if(
     &          (jaux.ne.4.and.
     &          (aux(jaux).eq.0.d0.or.aux_old(jaux).eq.0.d0)))
     &          then
C              This branch should not be taken because of above definition
C              of allow_log variables referring to both aux_old and aux.
              stop 'eos_jacobian: bad logic'
            else
              jacobian(jndex_aux,njacobian) =
     &          jacobian(jndex_aux,njacobian)/aux(jaux)
            endif
          else
            jacobian(jndex_aux,njacobian) =
     &        jacobian(jndex_aux,njacobian)
          endif
        endif
      enddo
      if(njacobian.gt.n_partial_aux) then
        if(kif.eq.1) then
          p = pnorad + pr
C          Additional RHS function that must be zeroed
          if(ifrad.gt.1) then
            if(pnorad.le.0.d0)
     &        stop 'eos_jacobian: negative pnorad calculated'
            faux(njacobian) = match_variable - log(pnorad)
            deriv_lnp_factor = 1.d0/pnorad
          else
            if(p.le.0.d0)
     &        stop 'eos_jacobian: negative p calculated'
            faux(njacobian) = match_variable - log(p)
            deriv_lnp_factor = 1.d0/p
          endif
C          negative partial derivative of faux wrt fl.
          jacobian(njacobian,njacobian) =
     &      deriv_lnp_factor*pnoradf
C          negative partial derivative of faux wrt old
C          auxiliary variables.
          do index_aux = 1,n_partial_aux
C            pnorad_daux already with respect to log(auxold) or log(-auxold)
C            if appropriate, see eos_warm_start.
            jacobian(njacobian,index_aux) = pnorad_daux(index_aux)
C            convert to ln(pnorad) or ln(p) derivatives
            jacobian(njacobian,index_aux) = deriv_lnp_factor*
     &        jacobian(njacobian,index_aux)
          enddo
        elseif(kif.eq.2) then
C          Additional RHS function that must be zeroed
          faux(njacobian) = match_variable - rl
C          negative partial derivative of faux wrt fl.
          jacobian(njacobian,njacobian) = rf
C          negative partial derivative of faux wrt old
C          auxiliary variables.
          do index_aux = 1,n_partial_aux
            jacobian(njacobian,index_aux) = 0.d0
            do index = 1,max_index
              jacobian(njacobian,index_aux) = 
     &          jacobian(njacobian,index_aux) +
     &          r_dv(index)*ddv_aux(index,index_aux)
            enddo
C            take derivative with respect to log(auxold) or log(-auxold)
            iaux = partial_aux(index_aux)
            if(allow_log(index_aux).or.
     &          allow_log_neg(index_aux)) then
              jacobian(njacobian,index_aux) =
     &          jacobian(njacobian,index_aux)*aux_old(iaux)
            endif
          enddo
        else
          stop 'eos_jacobian: bad internal logic'
        endif
      endif
      call free_non_ideal_calc(
     &  maxnextrasum,
     &  ifpi_local, ifnr, ifmodified,
     &  ifdvzero, inv_aux, naux, partial_ions, inv_ion,
     &  partial_elements, n_partial_elements,
     &  ion_end, mion_end,
     &  ifdv, ifrad,
     &  r_ion3, nion, nions, r_neutral, nelements,
     &  tl, f, wf, eta, rhostar, pstar,
     &  nux, nuy, nuz,
     &  if_pteh, full_sum0, full_sum1, full_sum2,
     &  if_mc, hcon_mc, hecon_mc, sum0_mc, sum2_mc,
     &  ifcoulomb_mod, if_dc,
     &  ifexcited, ifpl, ifsame_zero_abundances, ifh2, ifh2plus,
     &  izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax, bi,
     &  plop, plopt, plopt2,
     &  sumpl0, sumpl0f, sumpl0_dv,
     &  max_index, n_partial_aux,
     &  allow_log, allow_log_neg, aux_old, partial_aux,
     &  ddv_aux,
     &  aux, auxf, auxt, daux_dv, daux_daux, nextrasum,
     &  fion, fionf, fion_dv,
     &  ne, nef,
     &  free, free_daux)
      if(debug_jacobian) then
C        Choose jtemp1 ==> jtemp6 to correspond to auxiliary variables
C        actually used for particular free-energy model being tested.
C        reminder, there are only 15 of these variables for EOS1 in
C        order sum0, sum2, 7 neutral sums, 2 ionized sums, and 4 xextra
C        sums.
        if(.true.) then
          jtemp1 = 1
          jtemp2 = 2
          jtemp3 = 1
          jtemp4 = 2
          jtemp5 = 3
          jtemp6 = 4
        elseif(.false.) then
          jtemp1 = 7
          jtemp2 = 8
          jtemp3 = 9
          jtemp4 = 10
          jtemp5 = 11
          jtemp6 = 12
        elseif(.true.) then
C         11 and 12 repeat to get good alignment
          jtemp1 = 11
          jtemp2 = 12
          jtemp3 = 13
          jtemp4 = njacobian_final-2
          jtemp5 = njacobian_final-1
          jtemp6 = njacobian_final
        endif
        degeneracy(1) = faux(jtemp1)
        degeneracy(2) = -jacobian(jtemp1,itemp1)
        degeneracy(3) = -jacobian(jtemp1,itemp2)
        pressure(1) = faux(jtemp2)
        pressure(2) = -jacobian(jtemp2,itemp1)
        pressure(3) = -jacobian(jtemp2,itemp2)
        degeneracy(1) = pnorad
        degeneracy(2) = pnorad_daux(itemp1)
        degeneracy(3) = pnorad_daux(itemp2)
        pressure(1) = free
        pressure(2) = free_daux(itemp1)
        pressure(3) = free_daux(itemp2)
        density(1) = faux(jtemp3)
        density(2) = -jacobian(jtemp3,itemp1)
        density(3) = -jacobian(jtemp3,itemp2)
        energy(1) = faux(jtemp4)
        energy(2) = -jacobian(jtemp4,itemp1)
        energy(3) = -jacobian(jtemp4,itemp2)
        enthalpy(1) = faux(jtemp5)
        enthalpy(2) = -jacobian(jtemp5,itemp1)
        enthalpy(3) = -jacobian(jtemp5,itemp2)
        entropy(1) = faux(jtemp6)
        entropy(2) = -jacobian(jtemp6,itemp1)
        entropy(3) = -jacobian(jtemp6,itemp2)
        return
      endif
      end
