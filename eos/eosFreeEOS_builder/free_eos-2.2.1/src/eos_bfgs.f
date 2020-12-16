C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: eos_bfgs.f 822 2008-06-24 23:13:34Z airwin $
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
      subroutine eos_bfgs(
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
      logical allow_log_bfgs(naux_local),
     &  allow_log_neg_bfgs(naux_local)
      integer transform_bfgs(naux_local), index_bfgs
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
      integer morder, ifexchange_in
C      variables to help with restoring initial values.
      double precision fl_restore, aux_restore(naux_local)
C      variables to store acceptable points found by BFGS technique.
      double precision fl_acceptable, aux_acceptable(naux_local)
C      variables associated with the bfgs minimization.
      integer nbfgs
      double precision x(naux_local), xi(naux_local), dx(naux_local),
     &  func, gradient(naux_local), x_original(naux_local)
      double precision fold, step_ratio
      double precision fbar, epsilon, df_criterion, gradient_criterion
C      Setting fbar to a large negative number essentially disables the
C      Fletcher bracket-limiting test which should be fine for well-behaved
C      functions.  This parameter is in absolute units of func, but
C      func = free/(full_sum0*RT) is expected to stay within a few orders
C      of magnitude of unity.
      parameter(fbar = -1.d200)
C      epsilon is a delta f expected from roundoff error to add
C      some robustness in the presence of round-off error.  df_criterion
C      is the function value convergence parameter for the bfgs iteration.
C      gradient_criterion is the maximum absolute value of all gradient
C      components before can exit from bfgs iteration.
C      these three parameters are in absolute units of func, but
C      func = free/(full_sum0*RT) is expected to stay within a few orders
C      of magnitude of unity.
C      These values may be fairly critical so some experimentation may be
C      required for a large range of physical conditions to find the best
C      values.
      parameter(epsilon = 1.d-14)
      parameter(df_criterion = 1.d-14)
      parameter(gradient_criterion = 1.d-12)
      character*100 status, name
      integer linesearch_count, function_count, gradient_count,
     &  max_linesearch_count
C      each line search represents something like ~3 function and gradient
C      calculations and therefore represents something like ~10 calls to
C      eos_calc (due to the embedded fl iteration).  Thus, be careful that
C      you do not specify too large an integer here since the number of line
C      searches may often be the best termination criterion for the bfgs
C      technique.
      parameter(max_linesearch_count = 20)
      integer index_aux, iaux
      double precision max_gradient_component
      integer info
      logical debug_output
      data debug_output/.false./
      integer lnblnk
      save
      if(debug_output) then
        write(0,'(a,/,1p2d25.15)') 'entering eos_bfgs for '//
     &    'log rho, log T = ',
     &    match_variable/log(10.d0), tl/log(10.d0)
      endif
C      sanity checks.
      if(nextrasum.gt.maxnextrasum_local)
     &  stop 'eos_bfgs: maxnextrasum_local too small'
      if(.not.(ifnr.eq.0.or.ifnr.eq.1.or.ifnr.eq.3))
     &  stop 'eos_bfgs: ifnr must be 0, 1, or 3'
      if(ifnr.eq.0) stop 'eos_bfgs: ifnr = 0 is disabled'
      if(nions.ne.nions_local)
     &  stop 'eos_bfgs: nions must be equal to nions_local'
      if(naux.ne.naux_local)
     &  stop 'eos_bfgs: naux must be equal to naux_local'
      if(n_partial_aux.gt.naux_local)
     &  stop 'eos_bfgs: n_partial_aux too large'
      if(nxextrasum.gt.nxextrasum_local)
     &  stop 'eos_bfgs: nxextrasum too large'
      if(kif.ne.2) stop 'eos_bfgs: kif must be 2'

C      save values for later restoration.
      fl_restore = fl
      do index_aux = 1,n_partial_aux
        iaux = partial_aux(index_aux)
        aux_restore(iaux) = aux_old(iaux)
      enddo

C      starting solution:
      index_bfgs = 0
      do index_aux = 1,n_partial_aux
        iaux = partial_aux(index_aux)
C        ignore certain auxiliary variables for BFGS minimization.
        if(allow_log(index_aux).or.allow_log_neg(index_aux).or.
     &      iaux.eq.4) then
          index_bfgs = index_bfgs + 1
          transform_bfgs(index_bfgs) = iaux
C          eos_free calls eos_jacobian which can change allow_log and
C          allow_log_neg depending on whether aux agrees with aux_old
C          in sign.  Save allow_log_bfgs and allow_log_neg_bfgs here which
C          refer strictly to aux_old and which are indexed differently
C          so that the relationship between xi and aux_old remains consistent.
          allow_log_bfgs(index_bfgs) = allow_log(index_aux)
          allow_log_neg_bfgs(index_bfgs) = allow_log_neg(index_aux)
          if(allow_log_bfgs(index_bfgs)) then
            x(index_bfgs) = log(aux_old(iaux))
          elseif(allow_log_neg_bfgs(index_bfgs)) then
            x(index_bfgs) = log(-aux_old(iaux))
          elseif(iaux.eq.4) then
            x(index_bfgs) = aux_old(iaux)
          else
            stop 'eos_bfgs: logic error 1'
          endif
          x_original(index_bfgs) = x(index_bfgs)
        endif
      enddo
      nbfgs = index_bfgs
      status = 'start'
      linesearch_count = 0
      function_count = 0
      gradient_count = 0
      max_gradient_component = 2.d0*gradient_criterion
C      force at least two iterations.
      func = 1.d300
      fold = 2.d0*func
      do while(linesearch_count.le.max_linesearch_count.and.
     &    fold.gt.func+df_criterion.and.
     &    max_gradient_component.gt.gradient_criterion.and.
     &    status(:5).ne.'error')
        fold = func
C        the combination of this bfgs_iterate and the completion of
C        the following do while loop corresponds to the complete
C        initialization of the bfgs technique (for linesearch_count.eq.0)
C        or one line search followed by a bfgs correction in direction
C        (linesearch_count.gt.0)
        call bfgs_iterate(status, name, fbar, epsilon,
     &    nbfgs, x, xi, func, gradient, step_ratio, dx)
        do while(.not.
     &      (status(:5).eq.'error'.or.status(:8).eq.'complete'))
C          a gradient calculation always asked for right after a function
C          calculation done at the same xi.
          if(status(:8).eq.'function'.or.status(:4).eq.'both') then
            do index_bfgs = 1,nbfgs
              iaux = transform_bfgs(index_bfgs)
              if(allow_log_bfgs(index_bfgs)) then
                aux_old(iaux) = exp(xi(index_bfgs))
              elseif(allow_log_neg_bfgs(index_bfgs)) then
                aux_old(iaux) = -exp(xi(index_bfgs))
              else
                aux_old(iaux) = xi(index_bfgs)
              endif
            enddo
            if(.false..and.debug_output) then
              write(0,'(a,/,(1p5d25.15))')
     &          'eos_bfgs: xi= ',
     &          (xi(index_bfgs), index_bfgs = 1,nbfgs)
              write(0,'(a,/,(1p5d25.15))')
     &          'eos_bfgs: aux_old= ',
     &          (aux_old(partial_aux(index_aux)),
     &          index_aux = 1,n_partial_aux)
            endif
            call eos_free_calc(
     &        allow_log, allow_log_neg, partial_aux,
     &        ifrad,
     &        debug_dauxdv, itemp1, itemp2, delta1, delta2,
     &        debug_dvdaux, debug_jacobian,
     &        degeneracy, pressure, density,
     &        energy, enthalpy, entropy,
     &        iaux1, iaux2,
     &        match_variable, match_variable_save, kif, fl, tl_save,
     &        aux_old, aux, auxf, auxt, daux_dv,
     &        njacobian, njacobian_final,
     &        faux, jacobian, p, pr,
     &        ddv_aux, sumion0, sumion2,
     &        lambda, gamma_e,
     &        ifnr, nux, nuy, nuz, mion_end,
     &        n_partial_ions, n_partial_aux,
     &        sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &        partial_ions, t, morder, ifexchange_in,
     &        full_sum0, full_sum1, full_sum2, charge, charge2,
     &        dv_pl, dv_plt,
     &        ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &        ifsame_zero_abundances,
     &        ifexcited, ifsame_under,
     &        naux, inv_aux,
     &        inv_ion, max_index,
     &        partial_elements, n_partial_elements, ion_end,
     &        ifionized, if_pteh, if_mc, ifreducedmass,
     &        ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,
     &        ifpl, ifmodified, ifh2, ifh2plus,
     &        izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &        eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &        r_ion3, r_neutral, nelements, 
     &        ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &        ne, nef, net, sion, sionf, siont, uion,
     &        pnorad, pnoradf, pnoradt, pnorad_daux,
     &        free, free_daux,
     &        nextrasum, maxnextrasum,
     &        sumpl1, sumpl1f, sumpl1t, sumpl2,
     &        h2, h2f, h2t, h2_dv, info)
            if(info.ne.0) return
C            Calculate func and gradient for the requested xi.
            func = free/(full_sum0*cr*t)
            do index_bfgs = 1,nbfgs
              iaux = transform_bfgs(index_bfgs)
              index_aux = inv_aux(iaux)
              if(index_aux.le.0) stop 'eos_bfgs: logic error 2'
              gradient(index_bfgs) =
     &          free_daux(index_aux)/(full_sum0*cr*t)
C              correct xi gradient for differences between bfgs
C              allow_log variables used to compute relationship between
C              xi and aux_old and allow_log variables modified inside
C              eos_jacobian to calculate free_daux.
              if(allow_log_bfgs(index_bfgs)) then
                if(.not.(allow_log(index_aux).or.
     &            allow_log_neg(index_aux)))
     &            gradient(index_bfgs) = gradient(index_bfgs)*
     &            aux_old(iaux)
              elseif(allow_log_neg_bfgs(index_bfgs)) then
                if(.not.(allow_log(index_aux).or.
     &            allow_log_neg(index_aux)))
     &            gradient(index_bfgs) = gradient(index_bfgs)*
     &            aux_old(iaux)
              else
                if(allow_log(index_aux).or.
     &            allow_log_neg(index_aux)) then
C                  n.b. for this case aux_old must be non-zero (see allow_log
C                  and allow_log_neg logic in eos_jacobian).
                  gradient(index_bfgs) = gradient(index_bfgs)/
     &              aux_old(iaux)
                endif
              endif
            enddo
          endif
          if(status(:4).eq.'both') then
            function_count = function_count + 1
            gradient_count = gradient_count + 1
            if(.false..and.debug_output) then
              write(0,'(a,/,(1p5d25.15))') 'xi = ',
     &          (xi(index_aux), index_aux = 1,nbfgs)
              write(0,'(a,/,(1p5d25.15))') 'func = ', func
              write(0,'(a,/,(1p5d25.15))') 'gradient = ',
     &          (gradient(index_bfgs), index_bfgs = 1,nbfgs)
            endif
          elseif(status(:8).eq.'function') then
            function_count = function_count + 1
            if(.false..and.debug_output) then
              write(0,'(a,/,(1p5d25.15))') 'xi = ',
     &          (xi(index_bfgs), index_bfgs = 1,nbfgs)
              write(0,'(a,/,(1p5d25.15))') 'func = ', func
            endif
          elseif(status(:8).eq.'gradient') then
            gradient_count = gradient_count + 1
            if(.false..and.debug_output) then
              write(0,'(a,/,(1p5d25.15))') 'xi = ',
     &          (xi(index_bfgs), index_bfgs = 1,nbfgs)
              write(0,'(a,/,(1p5d25.15))') 'gradient = ',
     &          (gradient(index_bfgs), index_bfgs = 1,nbfgs)
            endif
          else
            stop 'eos_bfgs: logic error 3'
          endif
          call bfgs_iterate(status, name, fbar, epsilon,
     &      nbfgs, x, xi, func, gradient, step_ratio, dx)
        enddo
        if(status(:8).eq.'complete') then
C          save latest acceptable fl and aux_old (signalled by 'complete'
C          status code) found by BFGS technique.
C          n.b. this logic only works because the last evaluation done
C          by the bfgs line search is always an acceptable point.
          fl_acceptable = fl
          do index_aux = 1,n_partial_aux
            iaux = partial_aux(index_aux)
            aux_acceptable(iaux) = aux_old(iaux)
          enddo
          if(debug_output) then
            write(0,'(a,/,3i5)')
     &        'linesearch_count, function_count, gradient_count =',
     &        linesearch_count, function_count, gradient_count
            write(0,'(a,/,(1p5d25.15))') 'x = ',
     &        (x(index_bfgs), index_bfgs = 1,nbfgs)
            write(0,'(a,/,(1p5d25.15))') 'x - x_original = ',
     &        (x(index_bfgs)-x_original(index_bfgs),
     &        index_bfgs = 1,nbfgs)
            write(0,'(a,/,(1p5d25.15))') 'fl = ', fl
            write(0,'(a,/,(1p5d25.15))') 'fl - fl_original = ',
     &        fl-fl_restore
            if(linesearch_count.gt.0) then
              write(0,'(a,/,(1p5d25.15))') 'dx = ',
     &          (dx(index_bfgs), index_bfgs = 1,nbfgs)
              write(0,'(a,/,1pd25.15)') 'step_ratio = ', step_ratio
            else
              write(0,'(a,/,(1p5d25.15))') 'dx = ',
     &          (0.d0, index_bfgs = 1,nbfgs)
              write(0,'(a,/,1pd25.15)') 'step_ratio = ', 0.d0
            endif
            write(0,'(a,/,1pd25.15)') 'func = ', func
          endif
          max_gradient_component = 0.d0
          do index_bfgs = 1, nbfgs
            max_gradient_component =
     &        max(max_gradient_component, abs(gradient(index_bfgs)))
          enddo
          if(debug_output) then
            write(0,'(a,/,1pd25.15)') 'func - fold = ', func-fold
            write(0,'(a,/,1pd25.15)')
     &        'maximum absolute value of gradient component = ',
     &        max_gradient_component
            write(0,'(a,/,(1p5d25.15))') 'gradient = ',
     &        (gradient(index_bfgs), index_bfgs = 1,nbfgs)
          endif
C          increment at end of loop since linesearch = 0 corresponds to
C          status = 'start' which does initialization calculations at
C          initial x.
          linesearch_count = linesearch_count + 1
        !else
          !write(0,*) 'eos_bfgs: bad status = ', status(:lnblnk(status))
        endif
      enddo

      !return last acceptable point to calling routine.
      fl = fl_acceptable
      do index_aux = 1,n_partial_aux
        iaux = partial_aux(index_aux)
        aux_old(iaux) = aux_acceptable(iaux)
      enddo

      if(debug_output) then
        write(0,*) 'exiting eos_bfgs'
      endif
      end
