C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: eos_cold_start.f 842 2008-07-07 21:18:51Z airwin $
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
      subroutine eos_cold_start(
     &  ifsame_under, lambda, gamma_e,
     &  partial_ions, f, eta, wf, t, n_e, pstar, 
     &  dve_exchange, dve_exchangef, dve_exchanget,
     &  full_sum0, full_sum1, full_sum2, charge, dv_pl, dv_plt,
     &  ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &  ifsame_zero_abundances, ifexcited, ifnuform,
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
     &  nextrasum,
     &  ne, nef, net, sion, sionf, siont, uion,
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  h2, h2f, h2t,
     &  aux, auxf, auxt)

C      Determine the EOS (including all "output" auxiliary variables and
C      their derivatives) for a simplified form of free-energy model that
C      does not depend on "input" auxiliary variables. This routine is called
C      for both the case where the simplified EOS is the final result that is
C      wanted and the case where the "output" auxiliary variables calculated
C      for the simplified EOS are used for a starting approximation for a
C      more realistic EOS.
C
C      This routine first calculates dv and its partial derivatives wrt fl
C      and tl where dv is the change in the equilibrium constant of an ion
C      relative to the un-ionized reference state and the free electron dv =
C      (partial F/partial nref - partial F/partial nion - ion*partial
C      F/partial ne)/kT. N.B. with the simplified free-energy model, dv is
C      independent of "input" auxiliary variables.  After the dv, dvf, and
C      dvt values (held in arrays) are determined, eos_calc is called to
C      calculate the ionization/molecular equilibrium and all "output"
C      auxiliary variables and their fl and tl derivatives
         
C         input quantitites:
C         naux is the total number of auxiliary variables.
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
C         the following are returned only if ifexcited.gt.0 and
C           ifpi.eq.3.or4.
C         xextrasum(4) is the *negative* sum nuvar/(1 + qratio)*
C         partial qratio/partial extrasum(k), for k = 1, 2, 3, and nextrasum-1.
C         xextrasumf(4), xextrasumt(4) are the
C         fl and tl derivatives
      implicit none
      include 'constants.h'
      integer nions, nelements, nion(nions+2)
      integer ion_end(nelements+2)
      integer ifreducedmass, ifmodified, ifexcited,
     &  nmax
      double precision bi(nions+2), h2diss, t, tl,
     &  eps(nelements), eta
      double precision rhostar(9), pstar(9),
     &  ne, nef, net
      double precision dvzero(nelements), 
     &  dv(nions+2), dvf(nions+2), dvt(nions+2),
     &  dve, dvef, dvet, dve0, dve2, 
     &  dve_pi, dve_pif, dve_pit,
     &  dve_coulomb, dve_coulombf, dve_coulombt,
     &  dve_exchange, dve_exchangef, dve_exchanget,
     &  dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22,
     &  plop(nions+2), plopt(nions+2), plopt2(nions+2),
     &  dv_pl(nions+2), dv_plt(nions+2)
      integer maxcorecharge
C       maximum core charge for non-bare ion
      parameter(maxcorecharge = 28)
      integer nmin(maxcorecharge), nmin_max(maxcorecharge),
     & izlo, izhi
      integer nmin_species(nions+2)
      integer iatomic_number(nelements)
      integer ifnr,
     &  ifh2, ifh2plus, ifmtrace,
     &  if_pteh, ifcoulomb_mod, if_mc, if_dc,
     &  ifionized, ion, index_ion, max_index,
     &  index, maxnextrasum_local,
     &  nextrasum, ifpi_local, ifpl, nxextrasum_local,
     &  ifsame_abundances,
     &  ifsame_zero_abundances,
     &  ion0
      parameter (nxextrasum_local = 4)
      parameter (maxnextrasum_local = 9)
      double precision lambda, gamma_e,
     &  full_sum0, full_sum1, full_sum2,
     &  sum0, sum0ne, sum0f, sum0t,
     &  sum2, sum2ne, sum2f, sum2t,
     &  r_ion3(nions+2), r_neutral(nelements+2),
     &  charge(nions+2),
     &  extrasum(maxnextrasum_local), extrasumf(maxnextrasum_local), 
     &  extrasumt(maxnextrasum_local),
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  tc2
      integer nions_local
      parameter (nions_local = 295)
      double precision f, wf, n_e,
     &  sion, sionf, siont, uion,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het, rt, rf,
     &  h2, h2f, h2t, 
     &  h2plus, h2plusf, h2plust,
     &  sumpl0, sumpl0f
      integer ifdv(nions+2), ifelement(nelements),
     &  partial_ions(nions+2), inv_ion(nions+2),
     &  n_partial_elements, partial_elements(nelements+2)
      double precision rl
C       iextraoff is the number of non-extrasum auxiliary variables
      integer iextraoff, naux
      parameter (iextraoff = 7)
      integer 
     &  inv_aux(naux), ielement
      double precision 
     &  fion, fionf,
     &  xextrasum(nxextrasum_local), xextrasumf(nxextrasum_local),
     &  xextrasumt(nxextrasum_local),
     &  bmin(maxcorecharge)
      logical 
     &  ifsame_under, ifdvzero, ifnuform
C      For ifnr = 0 calculations, this array is never touched below so don't
C      bother with dimensions.
      double precision dummy
      logical ifcsum, ifextrasum, ifxextrasum, if_dv
      double precision aux(naux), auxf(naux), auxt(naux)
      save
      ifnr = 0
      do index = 1, n_partial_elements
        ielement = partial_elements(index)
        dvzero(ielement) = 0.d0
      enddo
      do index_ion = 1,max_index
        ion = partial_ions(index_ion)
        dv(ion) = 0.d0
        dvf(ion) = 0.d0
        dvt(ion) = 0.d0
      enddo
      if(ifpi_local.eq.0) then
        dve_pi = 0.d0
        dve_pif = 0.d0
        dve_pit = 0.d0
      else
        call pteh_pi(ifmodified,
     &    rhostar(1), rhostar(2), rhostar(3), t,
     &    dve_pi,
     &    dve_pif, dve_pit)
      endif
C      calculate change in equilibrium constant due to electron degeneracy
      dve = dve_pi - eta
      dvef = dve_pif - wf
      dvet = dve_pit
C      PTEH (full ionization) approximation to sum0, sum2
      sum0ne = full_sum0/full_sum1
      sum0 = sum0ne*n_e
      sum2ne = full_sum2/full_sum1
      sum2 = sum2ne*n_e
      sum0f = 0.d0
      sum0t = 0.d0
      sum2f = 0.d0
      sum2t = 0.d0
C      n.b. the 1 in the argument list is to_use_if_pteh = 1 inside
C      master_coulomb consistent with the above approximation.  This yields a
C      consistent calculation of dv and derivatives according to this
C      free-energy model approximation. However, later in the eos_calc
C      argument list the actual if_pteh value is used so that the output
C      sum0, sum2, and derivatives are calculated when needed for further
C      calculations without this cold-start approximation.
      call master_coulomb(ifnr,
     &  rhostar, f,
     &  sum0, sum0ne, sum0f, sum0t, sum2, sum2ne, sum2f, sum2t,
     &  n_e, t, cpe*pstar(1), pstar, lambda, gamma_e,
     &  ifcoulomb_mod, if_dc, 1,
     &  dve_coulomb, dve_coulombf, dve_coulombt, dve0, dve2,
     &  dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22)
      dve = dve + dve_coulomb
      dvef = dvef + dve_coulombf
      dvet = dvet + dve_coulombt
C      add in pre-calculated exchange effects
      dve = dve + dve_exchange
      dvef = dvef + dve_exchangef
      dvet = dvet + dve_exchanget
      do index = 1, n_partial_elements
        ielement = partial_elements(index)
        if(ifdvzero.or.ielement.eq.1) then
        else
          dvzero(ielement) = dvzero(ielement) +
     &      dv0 + charge(ion_end(ielement))*
     &      (dve + charge(ion_end(ielement))*dv2)
        endif
        if(ielement.gt.1) then
          ion0 = ion_end(ielement-1)
        else
          ion0 = 0
        endif
        do ion = ion0+1,ion_end(ielement)
          if(ifdvzero.or.ielement.eq.1) then
            dv(ion) = dv(ion) +
     &        dv0 + charge(ion)*(dve + charge(ion)*dv2)
          else
            dv(ion) = dv(ion) +
     &        dble(nion(ion)-nion(ion_end(ielement)))*
     &        (dve + dble(nion(ion)+nion(ion_end(ielement)))*dv2)
          endif
          dvf(ion) = dvf(ion) +
     &      dv0f + charge(ion)*(dvef + charge(ion)*dv2f)
          dvt(ion) = dvt(ion) +
     &      dv0t + charge(ion)*(dvet + charge(ion)*dv2t)
        enddo
      enddo
C      H2+
      if(ifdv(nions+2).eq.1) then
        dv(nions+2) = dv(nions+2) + dv0 + dve + dv2
        dvf(nions+2) = dvf(nions+2) + dv0f + dvef + dv2f
        dvt(nions+2) = dvt(nions+2) + dv0t + dvet + dv2t
      endif
C      add in Planck-Larkin occupation probability effect
      if(ifpl.eq.1) then
        do index_ion = 1,max_index
          ion = partial_ions(index_ion)
          dv(ion) = dv(ion) + dv_pl(ion)
          dvt(ion) = dvt(ion) + dv_plt(ion)
        enddo
      endif
C      apply Planck-Larkin excitation for pteh pressure ionization
C      N.B. note because ifpi = 1 for this call, extrasum, extrasumf,
C      extrasumt, xextrasum, xextrasumf, xextrasumt, and ddv_aux are
C      not read or written inside excitation_pi so all these
C      arguments can be replaced with dummy.
      if(ifexcited.gt.0.and.ifpl.eq.1)
     &  call excitation_pi(ifexcited,
     &  ifpl, 1, ifmodified, ifnr, ifsame_zero_abundances,
     &  inv_aux, iextraoff, naux, inv_ion, ifh2, ifh2plus,
     &  partial_elements, n_partial_elements, ion_end, 
     &  tl, izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &  bi, plop, plopt, plopt2,
     &  r_ion3, nion, nions, r_neutral, nelements,
     &  dummy, dummy, dummy, nextrasum,
     &  dummy, dummy, dummy,
     &  dv, dvf, dvt, dummy)
C      n.b. since ifnr = 0 (see above) all *_dv variables in this
C      call to eos_calc can be replaced by dummy since no *_dv variables
C      are read or written by eos_calc for ifnr = 0.
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
     &  dummy, dummy, dummy,
     &  fion, fionf, dummy,
     &  sum0, sum0f, sum0t, dummy,
     &  sum2, sum2f, sum2t, dummy,
     &  extrasum, extrasumf, extrasumt, dummy,
     &  nextrasum, maxnextrasum_local,
     &  sumpl0, sumpl0f, dummy,
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  rl, rf, rt, dummy,
     &  h2, h2f, h2t, dummy,
     &  h2plus, h2plusf, h2plust, dummy,
     &  xextrasum, xextrasumf, xextrasumt, dummy)

      ifcsum = .not.(if_mc.eq.1.or.if_pteh.eq.1)
      ifextrasum = ifpi_local.eq.3.or.ifpi_local.eq.4
      ifxextrasum = ifextrasum .and. (ifexcited.gt.0)
      if_dv = ifnr.eq.1.or.ifnr.eq.3
      call traditional_to_aux(
     &  if_mc, ifcsum, ifextrasum, ifxextrasum, if_dv,
     &  h, hf, ht, dummy,
     &  hd, hdf, hdt, dummy,
     &  he, hef, het, dummy,
     &  rl, rf, rt, dummy,
     &  h2plus, h2plusf, h2plust, dummy,
     &  sum0, sum0f, sum0t, dummy,
     &  sum2, sum2f, sum2t, dummy,
     &  iextraoff, maxnextrasum_local, nxextrasum_local,
     &  extrasum, extrasumf, extrasumt, dummy,
     &  xextrasum, xextrasumf, xextrasumt, dummy,
     &  nions, max_index, naux,
     &  aux, auxf, auxt, dummy)
      end
