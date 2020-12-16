C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: eos_free_calc.f 822 2008-06-24 23:13:34Z airwin $
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
      subroutine eos_free_calc(
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
     &  partial_ions, t, morder, ifexchange_in,
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
     &  ne, nef, net, sion, sionf, siont, uion,
     &  pnorad, pnoradf, pnoradt, pnorad_daux,
     &  free, free_daux,
     &  nextrasum, maxnextrasum,
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  h2, h2f, h2t, h2_dv, info)
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
C         info is a status code which is set to non-zero values if
C           there is a convergence problem in eos_free_calc.
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
     &  eps(nelements)
      double precision ne, nef, net
      integer maxfjs_aux
      parameter (maxfjs_aux = 4)
      double precision dvzero(nelements), 
     &  dv(nions+2), dvf(nions+2), dvt(nions+2),
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
     &  maxnextrasum,
     &  nextrasum, ifpi_local, ifpl, nxextrasum,
     &  ifsame_abundances,
     &  ifsame_zero_abundances
      parameter (nxextrasum = 4)
      double precision lambda, gamma_e,
     &  full_sum0, full_sum1, full_sum2,
     &  r_ion3(nions+2), r_neutral(nelements+2),
     &  charge(nions+2), charge2(nions+2),
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  tc2
      double precision
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
      logical allow_log(naux), allow_log_neg(naux)
      double precision aux_old(naux)
      integer partial_aux(naux)
      integer n_partial_aux
      double precision ddv_aux(nions+2,naux),
     &  sumion0, sumion2,
     &  sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &  nux, nuy, nuz
      integer
     &  mion_end, n_partial_ions
C      variables associated with pressure and free energy calculation:
      integer nions_local, naux_local, maxnextrasum_local,
     &  nxextrasum_local
      parameter(nions_local = 295)
      parameter(naux_local = 21)
      double precision rx(naux_local)
C      Note, this value must be greater than or equal to 4 as well to store
C      fjs_pi_free and corresponding pressure quantities.
      parameter(maxnextrasum_local = 9)
      parameter(nxextrasum_local = 4)
      double precision
     &  pnorad, pnoradf, pnoradt, pnorad_daux(n_partial_aux+2),
     &  free, free_daux(n_partial_aux+2)
C      for debugging
      logical debug_dauxdv, debug_dvdaux, debug_jacobian
      integer itemp1, itemp2
      double precision delta1, delta2
      double precision degeneracy(3), pressure(3), density(3), energy(3),
     &  enthalpy(3), entropy(3)
      integer iaux1, iaux2, kif, njacobian, njacobian_final
      double precision match_variable, match_variable_save, fl, tl_save,
     &  aux(naux), auxf(naux), auxt(naux), daux_dv(nions+2,naux)
      integer ifrad
      double precision faux(naux), jacobian(naux,naux), p, pr
C      variables needed for fl iteration
      double precision flold, paaplus, paaminus, flplus, flminus, paa,
     &  pac, pab, lnrho_5
      integer lter, iflast
      integer ltermax
C      maximum number of allowed ln f iterations
      parameter (ltermax = 200)
      double precision rhostar(9), pstar(9), sstar(3), ustar(3), 
     &  dve_exchange, dve_exchangef, dve_exchanget,
     &  f, wf, eta, n_e, rmue, fl_restore
      integer morder, ifexchange_in
      integer iaux, index_aux, index
      integer info
      save
C      default good status value returned.
      info = 0
      !temporary
!      write(0,*) 'entering eos_free_calc'
C      sanity checks.
      if(nextrasum.gt.maxnextrasum_local)
     &  stop 'eos_free_calc: maxnextrasum_local too small'
      if(.not.(ifnr.eq.0.or.ifnr.eq.1.or.ifnr.eq.3))
     &  stop 'eos_free_calc: ifnr must be 0, 1, or 3'
      if(ifnr.eq.0) stop 'eos_free_calc: ifnr = 0 is disabled'
      if(nions.ne.nions_local)
     &  stop 'eos_free_calc: nions must be equal to nions_local'
      if(naux.ne.naux_local)
     &  stop 'eos_free_calc: naux must be equal to naux_local'
      if(n_partial_aux.gt.naux_local)
     &  stop 'eos_free_calc: n_partial_aux too large'
      if(nxextrasum.gt.nxextrasum_local)
     &  stop 'eos_free_calc: nxextrasum too large'
      if(kif.ne.2) stop 'eos_free_calc: kif must be 2'
      !temporary
!      write(0,*) 'fixed aux_old for fl iteration'
!      write(0,'(5(1pd14.5,2x))')
!     &  (aux_old(partial_aux(index_aux)),
!     &  index_aux = 1, n_partial_aux)

      fl_restore = fl
C      iterate on fl until calculated density converges to match_variable.
C      Initialization for loop convergence criteria
      flold = 1.d30
      paaplus = 1.d30
      paaminus = -1.d30
C      mark undefined by ridiculous values
      flplus = -1.d30
      flminus = 1.d30
      lter = 0
      iflast = 0
      paa = 1.d0  !assure at least twice through loop.
      do while(iflast.ne.1)
        lter = lter + 1
C        Newton-Raphson iteration is quadratic, so obtain
C        machine precision (within significance loss noise of say
C        1.d-14) if do one more iteration after 10^-7 convergence.
        if(lter.ge.ltermax.or.abs(fl-flold).le.1.d-7) iflast = 1
C        find themodynamically consistent set of fermi-dirac integrals and
C        put the values and their derivatives into local versions
C        of rhostar, pstar, sstar, and ustar
C        and their f and t derivatives.  Also calculate exchange (when
C        ifexchange_in > 0) and its effects on dv.
        call fermi_dirac_exchange(fl, tl,
     &    rhostar, pstar, sstar, ustar, 3, morder, ifexchange_in,
     &    dve_exchange, dve_exchangef, dve_exchanget)
        f = exp(fl)
        wf = sqrt(1.d0 + f)
        eta = fl+2.d0*(wf-log(1.d0+wf))
C        number density of free electrons and electron pressure
        n_e = c_e*rhostar(1)
C        rho/mu_e = n_e H = cd*re
        rmue = cd*rhostar(1)
        call eos_jacobian(
     &    allow_log, allow_log_neg, partial_aux,
     &    ifrad, debug_dauxdv, itemp1, itemp2, delta1, delta2,
     &    debug_dvdaux, debug_jacobian,
     &    degeneracy, pressure, density, energy, enthalpy, entropy,
     &    iaux1, iaux2,
     &    match_variable, match_variable_save, kif, fl, tl_save,
     &    aux_old, aux, auxf, auxt, daux_dv,
     &    njacobian, njacobian_final,
     &    faux, jacobian, p, pr,
     &    ddv_aux, sumion0, sumion2,
     &    lambda, gamma_e,
     &    ifnr, nux, nuy, nuz, mion_end,
     &    n_partial_ions, n_partial_aux,
     &    sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &    partial_ions, f, eta, wf, t, n_e, pstar, 
     &    dve_exchange, dve_exchangef, dve_exchanget,
     &    full_sum0, full_sum1, full_sum2, charge, charge2,
     &    dv_pl, dv_plt,
     &    ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &    ifsame_zero_abundances,
     &    ifexcited, ifsame_under,
     &    naux, inv_aux,
     &    inv_ion, max_index,
     &    partial_elements, n_partial_elements, ion_end,
     &    ifionized, if_pteh, if_mc, ifreducedmass,
     &    ifsame_abundances, ifmtrace,
     &    iatomic_number, ifpi_local,
     &    ifpl, ifmodified, ifh2, ifh2plus,
     &    izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &    eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &    r_ion3, r_neutral, nelements, 
     &    ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &    rhostar,
     &    ne, nef, net, sion, sionf, siont, uion,
     &    pnorad, pnoradf, pnoradt, pnorad_daux,
     &    free, free_daux,
     &    nextrasum, maxnextrasum,
     &    sumpl1, sumpl1f, sumpl1t, sumpl2,
     &    h2, h2f, h2t, h2_dv)
        paa = faux(njacobian)
        pac = paa/jacobian(njacobian,njacobian)
        flold = fl       
C        Apply limits to potential fl change.
        if(fl + pac .lt. -10.d0) then
C         Cannot get into much trouble for derived  fl < -10.
          pab = min(20.d0,max(-20.d0,pac))
        else
          lnrho_5 = aux(4) - 1.5d0*(tl - log(1.d5))
          if(tl.gt.log(1.d6).or.
     &        lnrho_5 + auxf(4)*min(10.d0,pac).lt.log(1.d-3)) then
            pab = min(10.d0,max(-10.d0,pac))
          elseif(
     &        lnrho_5 + auxf(4)*min(3.d0,pac).lt.log(1.d-2)) then
            pab = min(3.d0,max(-3.d0,pac))
          elseif(
     &        lnrho_5 + auxf(4)*min(1.d0,pac).lt.log(1.d-1)) then
            pab = min(1.d0,max(-1.d0,pac))
          else
            pab = min(0.5d0,max(-0.5d0,pac))
          endif
        endif
        if(paa.ge.0.d0) then
          paaplus = paa
          flplus = fl
        else
          paaminus = paa
          flminus = fl
        endif
        if(paaminus.gt.-1.d30.and.paaplus.lt.1.d30) then
C          have bracket already!
C          make sure step would keep within bracket
          if((flminus.le.fl+pab.and.fl+pab.le.flplus).or.
     &        (flminus.ge.fl+pab.and.fl+pab.ge.flplus)) then
            fl = fl + pab
          else
            fl = 0.5d0*(flminus+flplus)
          endif
        else
C          if no bracket, yet, then normal step is only allowed for
C          when rf = auxf(4) = jacobian(njacobian,njacobian) is positive.
C          that is local paa is following global trend that paa
C          and rl generally increases with fl.
C          if rf negative with no bracket, 
C          then move out of region in the direction 
C          which should produce (global) bracket.
          if(auxf(4).lt.0.d0) then
            if(fl.eq.flminus) then
              pab = 0.5d0
            else
              pab = -0.5d0
            endif
          endif
          fl = fl + pab
        endif
        !temporary
!        write(0,'(a,/,i5,1p5d25.15)')
!     &    'lter, flold, paa, pac, pab, fl =',
!     &    lter, flold, paa, pac, pab, fl
C        End of iteration on fl to match match_variable with
C        eos_jacobian result.
      enddo
      if(lter.ge.ltermax) then
!        write(0,*) 'match_variable/ln(10), fl, tl/ln(10) ='
!        write(0,'(1p5d25.15)') match_variable/log(10.d0), fl,
!     &    tl/log(10.d0)
!        write(0,*) 'eps ='
!        write(0,'(1p5d25.15)') eps
!        write(0,*) 'ERROR: probably a multi-valued EOS '//
!     &    'for the eos_free_calc case'
!        write(0,*) 'flplus, paa(flplus), flminus, paa(flminus) ='
!        write(0,'(1p5d25.15)') flplus, paaplus, flminus, paaminus
!        write(0,'(a,/,i5,1p5d25.15)')
!     &    'lter, flold, paa, pac, pab, fl =',
!     &    lter, flold, paa, pac, pab, fl
!        write(0,*) 'eos_free_calc: too many fl iterations'
        info = 11
        return
      endif
C      implicitly eliminate fl dependence from gradient of free.
C      first calculate partial log rho wrt old auxiliary variables.
      do index_aux = 1,n_partial_aux
        rx(index_aux) = 0.d0
        do index = 1,max_index
          rx(index_aux) = rx(index_aux) +
     &      daux_dv(index,4)*ddv_aux(index,index_aux)
        enddo
C        take derivative with respect to log(auxold) or log(-auxold)
        iaux = partial_aux(index_aux)
        if(allow_log(index_aux).or.
     &      allow_log_neg(index_aux)) then
          rx(index_aux) = rx(index_aux)*aux_old(iaux)
        endif
      enddo
C      then eliminate fl dependence based on chain rule and rules of implicit
C      differentiation.
      do index_aux = 1,n_partial_aux
        free_daux(index_aux) = free_daux(index_aux) -
     &    (free_daux(n_partial_aux+1)/auxf(4))*rx(index_aux)
      enddo
      !temporary
!      write(0,'(a,/,(1p5d25.15))')
!     &  'eos_free_calc: aux_old= ',
!     &  (aux_old(partial_aux(index_aux)), index_aux = 1,n_partial_aux)
!      write(0,'(a,/,1pd25.15,/(1p5d25.15))')
!     &  'eos_free_calc: ln(rho) and its gradient (including fl) =',
!     &  aux(4), (rx(index_aux),
!     &  index_aux = 1,n_partial_aux), auxf(4)
!      
!      t = exp(tl)
!      write(0,'(a,/,(1p2d25.15))')
!     &  'eos_free_calc: fl, tl, pnorad, (scaled) free = ',
!     &  log(f), tl, pnorad, free/(full_sum0*cr*t)
!      write(0,'(a,/(1p5d25.15))')
!     &  'eos_free_calc: (scaled) free gradient =',
!     &  (free_daux(index_aux)/(full_sum0*cr*t),
!     &  index_aux = 1,n_partial_aux+1)
!      
      !temporary restore input
!      fl = fl_restore
!      do index_aux = 1,n_partial_aux
!        iaux = partial_aux(index_aux)
!        aux(iaux) = aux_old(iaux)
!      enddo
!      !temporary
!      write(0,*) 'exiting eos_free_calc'
      end
