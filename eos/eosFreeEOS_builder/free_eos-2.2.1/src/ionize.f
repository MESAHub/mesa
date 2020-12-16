C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: ionize.f 825 2008-06-28 20:56:54Z airwin $
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
      subroutine ionize(ifexcited, ifsame_under, ifnr, inv_ion,
     &  partial_elements, n_partial_elements, ion_end,
     &  ifionized, if_pteh, if_mc, ifreducedmass,
     &  ifsame_abundances, ifmtrace, iatomic_number, ifpi, ifpl,
     &  eps, tc2, bi, plop, plopt, plopt2,
     &  r_ion3, r_neutral, nelements_in, 
     &  ifelement, dvzero, dv, dvf, dvt, nion, nions_in,
     &  hne, hnef, hnet, hne_dv, s, sf, st, u, g, gf, gt, sh, shf, sht,
     &  hequil, hequilf, hequilt, h, hf, ht, hd, hdf, hdt, he, hef, het,
     &  h_dv, hd_dv, he_dv,
     &  fionh, fionhf, fionh_dv, fion, fionf, fion_dv,
     &  sum0, sum0f, sum0t, sum0_dv,
     &  sum2, sum2f, sum2t, sum2_dv,
     &  extrasum, extrasumf, extrasumt, extrasum_dv,
     &  nextrasum, maxnextrasum,
     &  sumpl0, sumpl0f, sumpl0_dv,
     &  sumpl1, sumpl1f, sumpl1t, sumpl2)
C       the purpose of this routine is to calculate ionization fractions
C       and weighted sums over those quantities.  the weighted sums are
C       ultimately used in an iterative way to calculate chemical potentials
C       according to some free energy model.  The appropriate difference in 
C       chemical potentials are then combined to form the dv quantities, 
C       the change in the equilibrium constant of an ion relative to the
C       un-ionized reference state and the free electron
C       dv = (partial F/partial nref - partial F/partial nion -
C       ion*partial F/partial ne)/kT.
C       These dv values (held in an array) are required input to ionize.
C
C       The free energy model:
C       it is the responsibility of the calling programme to get this
C       calculated correctly and fill in the dv quantities for the
C       various ions in a consistent manner.  It should be noted that
C       at the startup of the iteration, the free-energy terms corresponding
C       to pressure ionization and the Coulomb effect are approximated as
C       functions of only ne, so the dv array elements are quite similar
C       to each other.  (The electron exchange term is always just a function
C       of ne.) Later as more complicated models are used for the Coulomb
C       effect and pressure ionization, the dv array elements can become
C       more varied.
C       input quantitites:
C       ifexcited > 0 means_use_excited states (must have Planck-Larkin or
C         ifpi = 3 or 4).
C          0 < ifexcited < 10 means_use_approximation to explicit summation
C         10 < ifexcited < 20 means_use_explicit summation
C         mod(ifexcited,10) = 1 means just apply to hydrogen (without
C           molecules) and helium.
C         mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
C         mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
C       ifsame_under = .true. means_use_same underflow zeroing for each
C         ion as in previous call.  This option is useful for removing
C         small discontinuities caused by variations in the underflow
C         zeroing which sometimes foil the last stages
C         of convergence.
C       ifsame_under = .false. means calculate underflow zeroing.
C       ifnr = 0
C         calculate derivatives of nuvar (delivered through common block),
C         hne, sum0, sum2, s, extrasum, h, sh, hd, he, g, hequil, and sumpl1 
C         wrt f, t and sumpl0 wrt f with no other variables fixed, i.e., use
C         input dvf and dvt and chain rule.
C       N.B. chain of calling routines depend on the assertion that ifnr = 0
C       means no *_dv variables should be read or written.
C       ifnr = 1
C         calculate derivatives of nuvar (delivered through common block),
C         hne, sum0, sum2, extrasum, h, hd, and he wrt dv with f, t fixed.
C         n.b. list of variables is subset of those listed for ifnr = 0
C         because we only need dv derivatives for variables used to calculate
C         (output) auxiliary variables.
C       ifnr = 2 (should not occur)
C       ifnr = 3 same as combination of ifnr = 0 and ifnr = 1.
C         Note this is quite different in detail from ifnr = 3 interpretation
C         for many other routines, but the general motivation is the same
C         for all ifnr = 3 results; calculate both f and t derivatives and
C         other derivatives.
C       inv_ion(nions_in+2) maps ion index to contiguous ion index used in NR 
C         iteration dot products.
C       partial_elements(n_partial_elements+2) index of elements treated as 
C         partially ionized consistent with ifelement.
C       ion_end(nelements_in) keeps track of largest ion index for
C       each element.
C       ifionized = 0 means all elements treated as partially ionized
C       ifionized = 1 means trace metals fully ionized.
C       ifionized = 2 means all elements treated as fully ionized.
C       if_pteh controls whether pteh approximation used for Coulomb sums (1)
C         or whether_use_detailed Coulomb sum (0).
C       if_mc controls whether special approximation used for metal part
C         of Coulomb sum for the partially ionized metals.
C       if_mc = 1,_use_approximation
C       if_mc = 0, don't_use_approximation
C       ifreducedmass = 1 (use reduced mass in equilibrium constant)
C       ifreducedmass = 0 (use electron mass in equilibrium constant)
C       ifsame_abundances = 1, if call made with same abundances as previous
C         it is the responsibility of the calling programme to set
C         this flag where appropriate.  ifsame_abundances = 0 means
C         slower execution of ionize, but is completely safe against
C         abundance changes.
C       ifmtrace just to keep track of changes in ifelement.
C       n.b. if ifionized is different from previous call or 
C         ifsame_abundances = 0, or ifmtrace different from previous
C         call, then local variables ionized_hne,
C         ionized_sum0, ionized_sum2, and ionized_extrasum(2) are
C         recalculated.
C       iatomic_number(nelements_in), the atomic number of each element.
C       n.b. the code logic now demands there are iatomic_number ionized
C         states + 1 neutral state for each element.  Thus, the
C         ionization potentials must include all ions, not just
C         first ions (as for trace metals in old code).
C       eps(nelements_in) is an array of neps = nelements = 20 values
C         of relative abundance by weight divided by the appropriate atomic
C         weight. We_use_the atomic weight scale where (un-ionized) C(12) has
C         a weight of 12.00000000....  All weights are for the un-ionized
C         element. The eps value for an element should be the sum of the
C         individual isotopic eps values for that element.  The eps array
C         refers to the elements in the following order:
C         H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni
C       tc2 = c2/t
C       bi(nions_in+2) ionization potentials (cm**-1) in ion order:
C         h, he, he+, c, c+, etc.  H2 and H2+ are last two.
C       plop(nions_in+2), plopt, plopt2, planck-larkin occupation probability
C         and derivatives in ion order (starting with neutral in each case
C         and avoiding bare nuclei as in bi.
C       r_ion3(nions_in+2), *the cube of the*
C         effective radii (in ion order but going from neutral
C         to next to bare ion for each species) of MDH interaction between
C         all but bare nucleii species and ionized species.  last 2 are
C         H2 and H2+ (used externally)
C       r_neutral(nelements+2) effective radii for MDH neutral-neutral 
C         interactions.   Last two (used externally) are H2 and H2+ (the
C         only ionic species in the MHD model with a non-zero hard-sphere
C         radius).
C       ifelement(nelements_in) = 1 if element treated as partially ionized
C        (depending on abundance, and trace metal treatment), 0 if
C        element treated as fully ionized or has no abundance.
C       dvzero(nelements_in) is a zero point shift (zero for hydrogen) that
C         must be added to dv in all uses below.  Sometimes
C         (high density, high ionization) this quantity can be
C         larger than 1.d9 so if we avoid using it in differences
C         of dv quanties we can gain many significant digits.
C       dv(nions_in) change in equilibrium constant (explained above)
C       dvf(nions_in) = d dv/d ln f
C       dvt(nions_in) = d dv/d ln t
C       nion(nions_in), charge on ion in ion order (must be same order as bi)
C       e.g., for H+, He+, He++, etc.
C       output quantities:
C       hne, hnef, hnet, hne_dv is the accumulated nu(e) = n_e/(Navogadro*rho)
C         from the metals and helium
C         (and hydrogen if treated as fully ionized)
C         and its derivatives wrt ln f, and ln t, dv
C       s, sf, st is the accumulated entropy/R per unit mass and its 
C         derivatives with respect to lnf and ln t.
C       n.b. if not completely ionized, then the returned value of s excludes
C       the hydrogen component and sh (see below) returns the hydrogen
C       component (which is ignored from this routine and recalculated outside
C       when molecular formation is important).
C       if completely ionized s contains both the hydrogen and non
C       hydrogen components of the entropy.
C       n.b.  s returns a component of entropy/R per unit mass,
C       s = the sum over all non-electron species of
C       -nu_i*[-5/2 - 3/2 ln T + ln(alpha) - ln Na 
C         + ln (n_i/[A_i^{3/2} Q_i]) - dln Q_i/d ln T],
C       where nu_i = n_i/(Na rho), Na is the Avogadro number,
C       A_i is the atomic weight,
C       alpha =  (2 pi k/[Na h^2])^(-3/2) Na
C       (see documentation of alpha^2 in constants.h), and
C       Q_i is the internal ideal partition function
C       of the non-Rydberg states.  (Currently, we calculate this
C       partition function by the statistical weights of the ground states
C       of monatomic hydrogen and helium and the combined lower states (roughly
C       approximated) of each of the metals.  This is a superb approximation
C       for hydrogen, helium is subsequently corrected for detailed excitation
C       of the non-Rydberg states, and this crude approximation for the metals
C       is not currently corrected.  Thus, in all *current* monatomic
C       cases ln Q_i is a constant, and dln Q_i/d ln T is zero, but this
C       will change for the metals eventually.
C       From the equilibrium constant approach and the monatomic species
C       treated in this subroutine (molecules treated outside) we have
C       -nu_i ln (n_i/[A_i^{3/2} Q_i]) =
C       -nu_i * [ln (n_neutral/[A_neutral^{3/2} Q_neutral) +
C       (- chi_i/kT + dv_i)]
C       if we sum this term over all species
C       of an element without molecules (true for all species treated here
C       since molecular hydrogen treated specially outside when important)
C       we obtain
C       s_element = - eps * ln (eps*n_neutral/sum(n))
C       - eps ln (alpha/(A_neutral^{3/2} Q_neutral)) -
C       - sum over all species of the element of nu_i*(-chi/kT + dv(i))
C       where we have ignored the term
C       eps * (-5/2 - 3/2 ln T + ln rho)
C       (taking into account the first 4 terms above)
C       until a final correction of the s zero point
C       in free_eos_detailed.
C       u*avogadro, is the accumulated internal energy (cm^-1) per unit mass.
C       g, gf, gt is the h neutral fraction (n(h)/(n(h)+n(h+))) times eps(1)
C         and its derivatives with respect to ln f and ln t.
C         if no molecular formation of H2 and H2+, then
C         g =  n(H)/(rho*NA)
C       sh, shf, sht is the hydrogen contribution to the entropy and its
C         derivatives.  for the full ionization case, these terms are 
C         not calculated but instead are included in s and derivatives.
C         for the case where molecular formation is important, these
C         terms are not valid and are ignored.
C       hequil, hequilf, hequilt are the h equilibrium constant
C         ln [n(H+)/n(H)] and derivatives (only used for the case 
C         of h2 and h2+ formation).
C       h, hf, ht, h_dv is the h ionization fraction (n(h+)/(n(h)+n(h+)))
C         times eps(1) and its derivatives with respect to ln f, ln t,
C         and dv.
C         if no molecular formation of H2 and H2+, then
C         h =  n(H+)/(rho*NA)
C       hd, hdf, hdt, hd_dv = n(He+)/(rho*NA) and its 
C         derivatives with respect to ln f, ln t, and dv.
C       he, hef, het, he_dv = n(He++)/(rho*NA) and its 
C         derivatives with respect to ln f, ln t, and dv.
C         the following are returned only if ifpi = 3 or 4.
C       fionh, fionhf, fionh_dv = hydrogen ion free_energy/(rho*kT*NA)
C                     (= tc2*uh - sh at equilibrium)
C         and its derivatives with respect to ln f and dv.
C       fion, fionf, fion_dv = ion free_energy/(rho*kT*NA)
C                  (= tc2*u - s at equilibrium)
C         and its derivatives with respect to ln f and dv.
C       extrasum(maxnextrasum) weighted sums over non-H nu(i) = n(i)/(rho/H).
C         for iextrasum = 1,nextrasum-2, sum is only over neutral species and
C         weight is r_neutral^{iextrasum-1} for iextrasum = nextrasum-1, sum is
C         over all ionized species including bare nucleii, but excluding free
C         electrons, weight is Z^1.5. for iextrasum = nextrasum, sum is over
C         all species excluding bare nucleii and free electrons, the weight
C         is rion^3.
C       extrasumf(maxnextrasum) = partial of extrasum/partial ln f
C       extrasumt(maxnextrasum) = partial of extrasum/partial ln t
C       the following are returned only if ifpl = 1
C       sumpl0, sumpl1 and sumpl2 = weighted sums over
C         non-H nu(i) = n(i)/(rho/H).
C         for sumpl0 sum is over Planck-Larkin occupation probabilities.
C         for sumpl1 sum is over Planck-Larkin occupation probabilities +
C         d ln w/d ln T
C         for sumpl2, sum is over Planck-Larkin d ln w/d ln T
C       sumpl0f, sumpl0_dv = derivatives of sumpl0 wrt ln f and dv.
C       sumpl1f and sumpl1t = derivatives of sumpl1 wrt lnf and lnt.
      implicit none
      include 'aux_scale.h'
      include 'nuvar.h'
      integer nions
      parameter(nions = 295)
      integer nelements
      parameter(nelements = 20)  !number of elements considered.
      logical 
     &  ifsame_under, ifneutral_zero(nelements), ifion_zero(nions)
      integer ifexcited, ifnr,
     &  n_partial_elements, partial_elements(n_partial_elements+2),
     &  ifionized, if_pteh, if_mc, ifreducedmass, ifpi, ifpl,
     &  nelements_in, ifelement(nelements_in),
     &  nions_in, inv_ion(nions_in+2), 
     &  ion_end(nelements_in), nextrasum, maxnextrasum
      double precision ionized_hne, ionized_sum0, ionized_sum2,
     &  ionized_extrasum(2),
     &  ionized_u, ionized_s, eps(nelements_in), tc2,
     &  dvzero(nelements_in),
     &  dv(nions_in), dvf(nions_in), dvt(nions_in),
     &  hne, hnef, hnet, hne_dv(nions_in),
     &  s, sf, st, u, g, gf, gt, sh, shf, sht,
     &  hequil, hequilf, hequilt,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het,
     &  h_dv(nions_in), hd_dv(nions_in), he_dv(nions_in),
     &  fionh, fionhf, fionh_dv(nions_in),
     &  fion, fionf, fion_dv(nions_in),
     &  sum0, sum0f, sum0t, sum0_dv(nions_in),
     &  sum2, sum2f, sum2t, sum2_dv(nions_in),
     &  r_ion3(nions_in+2), r_neutral(nelements_in+2),
     &  extrasum(maxnextrasum), extrasumf(maxnextrasum),
     &  extrasumt(maxnextrasum), extrasum_dv(nions_in+2, nions_in),
     &  sumpl0, sumpl0f, sumpl0_dv(nions_in),
     &  sumpl1, sumpl1f, sumpl1t, sumpl2
C       ionization potentials in cm-1 in ion order for the ions,
C       H+, 
C       He+, He++,
C       Li, Be, B missing from list
C       C+ through C6+,
C       N+ through N7+,
C       O+ through O8+,
C       F missing from list
C       Ne+ through Ne10+,
C       Na+ through Na11+,
C       Mg+ through Mg12+,
C       Al+ through Al13+,
C       Si+ through Si14+,
C       P+ through P15+,
C       S+ through S16+,
C       Cl+ through Cl17+,
C       A+ through A18+,
C       K missing from list
C       Ca+ through Ca20+,
C       Sc missing from list
C       Ti+ through Ti22+,
C       V missing from list
C       Cr+ through Cr24+,
C       Mn+ through Mn25+,
C       Fe+ through Fe26+,
C       Co missing from list
C       Ni+ through Ni28+,
      double precision bi(nions_in), bi_ref(nions),
     &  ce0(nions), ce(nions),
     &  logqtl_const0(nelements), logqtl_const(nelements+nions+2),
     &  plop(nions_in), plopt(nions_in), plopt2(nions_in)
      include 'statistical_weights.h'
      include 'constants.h'
      integer  nion(nions_in)  !ionization state in ion order.
C      atomic masses of most abundant isotope.  this is most consistent way
C      to treat isotopes in equilibrium constants (aside from having
C      separate abundance constraint equations for each isotope).  These
C      numbers *should not be confused* with the mean atomic weights
C      required for abundance and density calculations.  The element order
C      is H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni. We_use_the
C      atomic weight scale where (un-ionized) C(12) has a weight of
C      12.00000000....  All weights are for the un-ionized element.
      double precision atomic_mass(nelements)
      data atomic_mass/
     &  1.007825035d0, 4.00260324d0,
     &  12.0000000d0, 14.003074002d0, 15.99491463d0,
     &  19.9924356d0, 22.9897677d0, 23.9850423d0,
     &  26.9815386d0, 27.9769271d0, 30.9737620d0,
     &  31.972070698d0, 34.968852728d0, 39.9623837d0,
     &  39.9625906d0, 47.9479473d0, 51.9405098d0,
     &  54.9380471d0, 55.9349393d0, 57.9353462d0/
      integer maxionstage
      parameter(maxionstage = 28)  !treat up to 28 ions of one element.
      double precision fract(maxionstage), fractf(maxionstage),
     &  fractt(maxionstage), fractneutral, vz_dv(maxionstage),
     &  rpower_dv(maxionstage), logalpha, constant_sh, constant_s
      integer ifsame_abundances,
     &  ifmtrace, iatomic_number(nelements_in),
     &  ifmtrace_old, ifionized_old, ifpi_old, if_pteh_old
      data ifmtrace_old, ifionized_old, ifpi_old,
     &  if_pteh_old/4*-10000/
      integer iffirst
      data iffirst/1/
      integer ion, ielement, index_f, jndex_f, ion0, ion_index,
     &  inv_ion_index, index_max, index_maxnb, iextrasum,
     &  index_element
      double precision fractmax, sum_full, logsum0, logsum0f,
     &  logsum0_dv(maxionstage), arg,
     &  sum, sumf, sumt, sumnb,
     &  vz, vzf, vzt, sarg, sarg0, epsilon, rpower, rpowerf,
     &  rpowert
C      exp(2*arglim) is the maximum ratio of largest to least ionization
C      component.
C      n.b. typical overflow limit ~ 1.d305, but we back off
C      by more than a factor of two because subsequent calculations use
C      the square of these factors, and we want to leave some room for
C      bad scaling.
C      exp(2*2*arglim) = exp(2*2*144) ~ 1.d250 which gives ~50 orders
C      of magnitude of room for bad scaling.  Thus, 144 is the maximum
C      arglim.
      double precision arglim
      parameter (arglim = 144.d0)
C      144.d0 is not gross overkill.  The whole EOS is expressed in terms of
C      ne so *some* ionized species must be present.
C      Some of the auxiliary variables depend
C      only on neutrals, only on ions, or only on non-bare ions.  Thus,
C      near full neutralization or near full ionization can cause 
C      discontinuous behaviour whenever some species is zeroed because
C      of underflow concerns.  These discontinuities shouldn't matter
C      for converged EOS results, but the discontinuities will mess
C      up convergence criteria if the discontinuity occurs near the
C      converged solution.  Even with a criterion of 144.d0, we still
C      found this behaviour so ifsame_under logic was introduced to
C      avoid these discontinuities (see ifsame_under logic comments above).
      integer ifreducedmassold
C      need invalid value
      data ifreducedmassold/-1/
      logical ifpi34, ifnr03, ifnr13, ifcsum
      logical iffion
      parameter(iffion = .true.)
      save
C      sanity checks:
      if(ifnr.lt.0.or.ifnr.eq.2.or.ifnr.gt.3) stop
     &  'ionize: ifnr must be 0, 1, or 3'
      if(if_pteh.eq.1.and.if_mc.eq.1) stop
     &  'ionize: if_pteh and if_mc cannot simultaneously be non-zero'
      if(nions_in.ne.nions.or.nions_in.ne.nions_stat)
     &  stop 'ionize: invalid nions_in'
      if(nelements_in.ne.nelements.or.nelements_in.ne.nelements_stat)
     &  stop 'ionize: invalid nelements_in'
      if(ifpi.lt.0.or.ifpi.gt.4)
     &  stop 'ionize: invalid ifpi'
      if(.not.((ifpi.ne.3.and.ifpi.ne.4.and.nextrasum.eq.0).or.
     &  ((ifpi.eq.3.or.ifpi.eq.4).and.
     &  (nextrasum.eq.6.or.nextrasum.eq.9))))
     &  stop 'ionize: invalid nextrasum'
      if(iffirst.eq.1.and.ifsame_abundances.eq.1)
     &  stop 'ionize: invalid ifsame_abundances on first call'
      if(iffirst.eq.1) then
        iffirst = 0
C        alpha =  (2 pi k/[Na h^2])^(-3/2) Na
C        (see documentation of alpha^2 in constants.h).
        logalpha = 0.5d0*log(alpha2)
C         refer energies to neutral species.
        do ion = 1,nions
          bi_ref(ion) = bi(ion)
          if(nion(ion).gt.1)
     &      bi_ref(ion) = bi_ref(ion) + bi_ref(ion-1)
        enddo
        ielement = 0
        do ion = 1,nions
          if(nion(ion).eq.1) ielement = ielement + 1
          ce0(ion) = log(dble(iqion(ion))/dble(iqneutral(ielement)))
        enddo
        if(ifsame_under)
     &    stop 'ionize: ifsame_under must be .false. on first call'
      endif
C      useful combinations and flags in following tests:
      ifpi34 = ifpi.eq.3.or.ifpi.eq.4
      ifnr03 = ifnr.eq.0.or.ifnr.eq.3
      ifnr13 = ifnr.eq.1.or.ifnr.eq.3
      ifcsum = .not.(if_mc.eq.1.or.if_pteh.eq.1)
      if(ifreducedmass.ne.ifreducedmassold) then
        ifreducedmassold = ifreducedmass
        ielement = 0
        do ion = 1,nions
          ce(ion) = ce0(ion)
          if(nion(ion).eq.1) then
            ielement = ielement + 1
C            temperature-independent part of log product of translational and
C            internal partition functions/avogadro
C            logqtl_const0 is the part that is constant for each element
            logqtl_const0(ielement) =
     &        1.5d0*log(atomic_mass(ielement)) - logalpha
            logqtl_const(ielement) = log(dble(iqneutral(ielement))) +
     &        logqtl_const0(ielement)
          endif
          logqtl_const(nelements+ion) = log(dble(iqion(ion))) +
     &      logqtl_const0(ielement)
          if(ifreducedmass.eq.1) then
            ce(ion) = ce(ion) +
     &      1.5d0*log((atomic_mass(ielement) -
     &        dble(nion(ion))*electron_mass)/atomic_mass(ielement))
            logqtl_const(nelements+ion) =
     &        logqtl_const(nelements+ion) + 
     &      1.5d0*log((atomic_mass(ielement) -
     &        dble(nion(ion))*electron_mass)/atomic_mass(ielement))
          endif
        enddo
      endif
C      calculate constant component of entropy.  Keep hydrogen separate
C      for now, but add in with constant_s (see below) in full ionization
C      case.
      if(ifsame_abundances.ne.1) then
        constant_sh = eps(1)*(1.5d0*log(atomic_mass(1)) +
     &    log(dble(iqneutral(1))) - logalpha)
        constant_s = 0.d0
        do ielement = 2, nelements
          constant_s = constant_s +
     &      eps(ielement)*(1.5d0*log(atomic_mass(ielement)) +
     &      log(dble(iqneutral(ielement))) - logalpha)
        enddo
      endif
C      ionized_hne is the contribution of fully ionized elements 
C      (depending on ifelement and ifionized) to hne.
C      ionized_u is the contribution of fully ionized elements 
C      (depending on ifelement and ifionized) to u.
C      ionized_s is the contribution of fully ionized elements 
C      (depending on ifelement and ifionized) to s
C      ionized_sum0 is the contribution of fully ionized elements 
C      (depending on ifelement and ifionized) to sum0.
C      ionized_sum2 is the contribution of fully ionized elements 
C      (depending on ifelement and ifionized) to sum2.
C      ionized_extrasum(2) is the contribution of fully ionized elements 
C      (depending on ifelement and ifionized) to the last two elements of 
C      extrasum.
      if(ifionized.ne.ifionized_old.or.ifpi.ne.ifpi_old.or.
     &    if_pteh.ne.if_pteh_old.or.ifsame_abundances.ne.1.or.
     &    ifmtrace.ne.ifmtrace_old) then
C        if different ifionized flag, or if the pressure ionization
C        has changed, or if the Coulomb approximation has changed,
C        or if the abundances have changed,
C        or ifmtrace (ifelement) has changed,
C        then recalculate full-ionized sum approximations.
        ifionized_old = ifionized
        ifpi_old = ifpi
        if_pteh_old = if_pteh
        ifmtrace_old = ifmtrace
        ionized_hne = 0.d0
        ionized_u = 0.d0
        ionized_s = 0.d0
        ionized_sum0 = 0.d0
        ionized_sum2 = 0.d0
        if(ifpi34) then
          ionized_extrasum(1) = 0.d0
          ionized_extrasum(2) = 0.d0
        endif
        index_f = 0
        do ielement = 1, nelements
C          index_f should point to last index of ion for this
C          ielement.
          index_f = index_f + iatomic_number(ielement)
C          if element treated as fully ionized (or zero abundance)
C          ifelement(ielement) = 0.  For ifionized = 2, overides
C          ifelement and all elements are treated as fully ionized.
          if(eps(ielement).gt.0.d0.and.
     &        (ifionized.eq.2.or.ifelement(ielement).eq.0)) then
            ionized_hne = ionized_hne +
     &        dble(iatomic_number(ielement))*eps(ielement)
            ionized_u = ionized_u + eps(ielement)*bi_ref(index_f)
C            this expression is limit of expressions below
C            when all lower ions may be finite but negligible
C            relative to bare nucleus. n.b. much cancellation!
            ionized_s = ionized_s + eps(ielement)*
     &        (ce(index_f) - log(eps(ielement)))
            if(if_pteh.ne.1) then
              ionized_sum0 = ionized_sum0 +  eps(ielement)
              ionized_sum2 = ionized_sum2 +
     &          dble(iatomic_number(ielement))**2*eps(ielement)
            endif
            if(ifpi34) then
              ionized_extrasum(1) = ionized_extrasum(1) +
     &          dble(iatomic_number(ielement))**1.5d0*eps(ielement)
C              n.b. there should be *no* contribution
C              of bare nuclei to last extrasum
            endif
          endif
        enddo
      endif
      if(ifnr03) then
        hnef = 0.d0
        hnet = 0.d0
        sum0f = 0.d0
        sum0t = 0.d0
        sum2f = 0.d0
        sum2t = 0.d0
        sf = 0.d0
        st = 0.d0
        fionf = 0.d0
        fionhf = 0.d0
      endif
      if(ifpi34) then
        do iextrasum = 1,nextrasum
          extrasum(iextrasum) = 0.d0
          if(ifnr03) then
            extrasumf(iextrasum) = 0.d0
            extrasumt(iextrasum) = 0.d0
          endif
        enddo
      endif
      if(ifpl.eq.1) then
C        n.b. for fully ionized elements, there is no contribution
C        to Planck-Larkin sums.
        sumpl0 = 0.d0
        sumpl1 = 0.d0
        sumpl2 = 0.d0
        if(ifnr03) then
          sumpl0f = 0.d0
          sumpl1f = 0.d0
          sumpl1t = 0.d0
        endif
      endif
      if(ifionized.eq.2) then
C        everything (H, He, metals) fully ionized and ifnr irrelevant
        hne = ionized_hne
        u = ionized_u
        s = ionized_s + constant_s + constant_sh
C        verified as complete ionization limit for detailed expression below
C        for fion.
        fion = tc2*u - s
        sum0 = ionized_sum0*sum0_scale
        sum2 = ionized_sum2*sum2_scale
        if(ifpi34) then
          extrasum(nextrasum-1) = ionized_extrasum(1)*
     &      extrasum_scale(nextrasum-1)
          extrasum(nextrasum) = ionized_extrasum(2)*
     &      extrasum_scale(nextrasum)
        endif
      else
        if(ifionized.eq.1) then
C          trace metals fully ionized, hydrogen effect split off
          hne = ionized_hne
          u = ionized_u
          s = ionized_s + constant_s
C          verified as limit for detailed expression below
C          for fion of completely ionized trace metals.
          fion = tc2*u - s
          fionh = 0.d0
          sum0 = ionized_sum0*sum0_scale
          sum2 = ionized_sum2*sum2_scale
          if(ifpi34) then
            extrasum(nextrasum-1) = ionized_extrasum(1)*
     &      extrasum_scale(nextrasum-1)
            extrasum(nextrasum) = ionized_extrasum(2)*
     &      extrasum_scale(nextrasum)
          endif
        else
C          ifionized not 2 or 1 so everything treated as partially
C          ionized, hydrogen effect split off
          hne = 0.d0
          u = 0.d0
          s = constant_s
C          fion = tc2*u - s
          fion = 0.d0
          fionh = 0.d0
          sum0 = 0.d0
          sum2 = 0.d0
C          the following commands are redundant so comment out
C          if(ifpi34) then
C             extrasum(nextrasum-1) = 0.d0
C             extrasum(nextrasum) = 0.d0
C          endif
        endif
        nuvar_nelements = n_partial_elements
        do index_element = 1, n_partial_elements
          ielement = partial_elements(index_element)
          epsilon = eps(ielement)
          if(ielement.gt.1) then
            ion0 = ion_end(ielement-1)
          else
            ion0 = 0
          endif
C          ln fraction of neutral relative to neutral.
C          to get previous zero point must add dvzero to
C          fractmax and fract, but such zero point shifts
C          don't matter so long as they are constant for
C          an element.
          fractmax = -dvzero(ielement)
          index_max = 0
C          n.b., the index_maxnb logic below does not get triggered for the
C          hydrogen case (iatomic_number(ielement) = 1) or used later on,
C          but that is an accident waiting to happen (in case the subsequent
C          code is changed) so always initialize index_maxnb here.
          index_maxnb = 0
          do index_f = 1, iatomic_number(ielement)
            ion = ion0+index_f
            fract(index_f) = ce(ion) + dv(ion) - bi_ref(ion)*tc2
            if(fract(index_f).ge.fractmax) then
              index_max = index_f
              fractmax = fract(index_f)
            endif
C            this is the maximum component excluding the
C            bare nucleus (used only for ifpi.eq.3 or 4)
            if(index_f.eq.iatomic_number(ielement)-1)
     &        index_maxnb = index_max
            if(ifnr03) then
              fractf(index_f) = dvf(ion)
              fractt(index_f) = dvt(ion) + bi_ref(ion)*tc2
            endif
          enddo
C          calculate ionization fractions and their derivatives.
C          normalization near overflow limit exp(arglim) for maximum accuracy
          fractmax = fractmax - arglim
C          evaluate neutral contribution to sum where 
C          fraction = 1, ln fract = 0.
          arg = -fractmax - dvzero(ielement)
C          normalized so that maximum fraction is exp(arglim) so 
C          that exp(-arglim) or less can safely be ignored.
C          n.b. note ifsame_under = .false. check on first call above.
          if(ifsame_under) then
            if(ifneutral_zero(ielement)) then
              fractneutral = 0.d0
            else
              fractneutral = exp(arg)
            endif
          else
            if(arg.le.-arglim) then
              fractneutral = 0.d0
              ifneutral_zero(ielement) = .true.
            else
              fractneutral = exp(arg)
              ifneutral_zero(ielement) = .false.
            endif
          endif
C          note that sum skips maximum component to solve significance
C          loss problems.
          if(index_max.ne.0) then
            sum = fractneutral
          else
            sum = 0.d0
          endif
          if(ifpi34) sumnb = 0.d0
          logsum0 = arg
          if(ifnr03) then
            sumf = 0.d0
            sumt = 0.d0
            logsum0f = 0.d0
          endif
          do jndex_f = 1,iatomic_number(ielement)
C            difference in zero point doesn't matter here
C            and may reduce significance loss.
            arg = fract(jndex_f) - fractmax
C            normalized so that maximum fraction is exp(arglim) so 
C            that exp(-arglim) or less can safely be ignored.
C          n.b. note ifsame_under = .false. check on first call above.
            if(ifsame_under) then
              if(ifion_zero(ion0+jndex_f)) then
                fract(jndex_f) = 0.d0
              else
                fract(jndex_f) = exp(arg)
              endif
            else
              if(arg.le.-arglim) then
                fract(jndex_f) = 0.d0
                ifion_zero(ion0+jndex_f) = .true.
              else
                fract(jndex_f) = exp(arg)
                ifion_zero(ion0+jndex_f) = .false.
              endif
            endif
C            sum, sumf, sumt skips maximum component to solve
C            significance loss problems.
            if(jndex_f.ne.index_max) then
              sum = sum + fract(jndex_f)
              if(ifnr03) then
                sumf = sumf + fract(jndex_f)*fractf(jndex_f)
                sumt = sumt + fract(jndex_f)*fractt(jndex_f)
              endif
            endif
            if(ifpi34.and.
     &        jndex_f.lt.iatomic_number(ielement))
     &        sumnb = sumnb + fract(jndex_f)
          enddo
          if(index_max.ne.0) then
            sum_full = sum + fract(index_max)
          else
            sum_full = sum + fractneutral
          endif
C          n.b. -logsum0 = log(sum_full/fract_neut)
C          = log(fract_max/fract_neut) + log(1+sum/fract_max)
          logsum0 = logsum0 - log(sum_full)
          if(ifnr03) then
            if(index_max.ne.0) then
              logsum0f = logsum0f -
     &          (sumf + fract(index_max)*fractf(index_max))/sum_full
            else
              logsum0f = logsum0f - sumf/sum_full
            endif
C            convert to derivatives of ln(sum_full) with missing maximum
C            component (except where the maximum is the neutral component
C            in which case the derivative of that component is zero
C            in any case.)
            sumf = sumf/sum_full
            sumt = sumt/sum_full
          endif
          if(ifnr13) then
            if(ifpi34) then
              do iextrasum = 1,nextrasum
                do index_f = ion0+1,ion0+iatomic_number(ielement)
                  inv_ion_index = inv_ion(index_f)
                  extrasum_dv(inv_ion_index,iextrasum) = 0.d0
                enddo
              enddo
            endif
            do jndex_f = ion0+1,ion0+iatomic_number(ielement)
              inv_ion_index = inv_ion(jndex_f)
              if(ifcsum) then
                sum2_dv(inv_ion_index) = 0.d0
              endif
C              h, hd, he derivatives ignored outside except
C              for relevant dv indices *always* updated below.
C              h_dv(inv_ion_index) = 0.d0
C              hd_dv(inv_ion_index) = 0.d0
C              he_dv(inv_ion_index) = 0.d0
              hne_dv(inv_ion_index) = 0.d0
              fion_dv(inv_ion_index) = 0.d0
              fionh_dv(inv_ion_index) = 0.d0
              if(ifpl.eq.1) sumpl0_dv(inv_ion_index) = 0.d0
              logsum0_dv(jndex_f-ion0) = -fract(jndex_f-ion0)/
     &          sum_full
            enddo
          endif
          rpower = fractneutral*epsilon/sum_full
          if(rpower.gt.0.d0.and.iffion) then
            if(ielement.eq.1) then
              fionh = fionh + rpower*
     &          (log(rpower) - logqtl_const(ielement))
            else
              fion = fion + rpower*
     &          (log(rpower) - logqtl_const(ielement))
            endif
          endif
C          in all cases and all indices save n/(rho*avogadro) in nuvar
C          this is needed by excitation_sum
C          and eos_sum_calc and eos_free_calc.
          nuvar(1,index_element) = rpower
          if(ielement.eq.1.or.
     &        (rpower.gt.0.d0.and.
     &        (iffion.or.ifpl.eq.1.or.ifpi34))) then
            if(ifnr03) then
              if(index_max.ne.0) then
C                sumf, sumt skip maximum component.
                rpowerf = -rpower*(sumf +
     &            fract(index_max)*fractf(index_max)/sum_full)
                rpowert = -rpower*(sumt +
     &            fract(index_max)*fractt(index_max)/sum_full)
              else
C                fractf, fractt = 0, and sumf, sumt complete, for
C                index_max = 0.
                rpowerf = -rpower*sumf
                rpowert = -rpower*sumt
              endif
            endif
            if(ifnr13.and.(iffion.or.ifpl.eq.1.or.ifpi34)) then
              do index_f = 1, iatomic_number(ielement)
C                n.b. d *ln* rpower/d dv
                rpower_dv(index_f) = -fract(index_f)/sum_full
              enddo
            endif
            if(rpower.gt.0.d0.and.iffion) then
              if(ielement.eq.1) then
                if(ifnr03) then
                  fionhf = fionhf + rpowerf*
     &              (log(rpower) - logqtl_const(ielement) + 1.d0) 
                endif
                if(ifnr13) then
                  do index_f = 1, iatomic_number(ielement)
                    fionh_dv(index_f) = fionh_dv(index_f) +
     &                rpower*rpower_dv(index_f)*
     &                (log(rpower) - logqtl_const(ielement) + 1.d0)
                  enddo
                endif
              else
                if(ifnr03) then
                  fionf = fionf + rpowerf*
     &              (log(rpower) - logqtl_const(ielement) + 1.d0) 
                endif
                if(ifnr13) then
                  do index_f = 1, iatomic_number(ielement)
                    inv_ion_index = inv_ion(index_f+ion0)
                    fion_dv(inv_ion_index) = fion_dv(inv_ion_index) +
     &                rpower*rpower_dv(index_f)*
     &                (log(rpower) - logqtl_const(ielement) + 1.d0)
                  enddo
                endif
              endif
            endif
            
            if(ifexcited.gt.0) then
              if(ifnr03) then
                nuvarf(1,index_element) = rpowerf
                nuvart(1,index_element) = rpowert
              endif
              if(ifnr13) then
                do index_f = 1, iatomic_number(ielement)
                  nuvar_dv(index_f,1,index_element) =
     &              rpower*rpower_dv(index_f)
                enddo
              endif
            endif
            if(ielement.eq.1) then
C              save special quantities for hydrogen.
              g = rpower
              hequil = ce(1) + dv(1) - bi_ref(1)*tc2
              if(ifnr03) then
                gf = rpowerf
                gt = rpowert
                hequilf = dvf(1)
                hequilt = dvt(1) + bi_ref(1)*tc2
              endif
            else
              if(ifpl.eq.1) then
C                Planck-Larkin occupation probability sums
                sumpl0 = sumpl0 + rpower*plop(ion0+1)
                sumpl1 = sumpl1 + rpower*(plop(ion0+1)+plopt(ion0+1))
                sumpl2 = sumpl2 + rpower*plopt(ion0+1)
                if(ifnr03) then
                  sumpl0f = sumpl0f +
     &              rpowerf*plop(ion0+1)
                  sumpl1f = sumpl1f +
     &              rpowerf*(plop(ion0+1)+plopt(ion0+1))
                  sumpl1t = sumpl1t +
     &              rpowert*(plop(ion0+1)+plopt(ion0+1)) +
     &              rpower*(plopt(ion0+1)+plopt2(ion0+1))
                endif
                if(ifnr13) then
                  do index_f = 1,iatomic_number(ielement)
                    inv_ion_index = inv_ion(index_f+ion0)
                    sumpl0_dv(inv_ion_index) = 
     &                sumpl0_dv(inv_ion_index) +
     &                rpower_dv(index_f)*rpower*plop(ion0+1)
                  enddo
                endif
              endif
              if(ifpi34) then
                extrasum(nextrasum) = extrasum(nextrasum) +
     &            rpower*extrasum_scale(nextrasum)*r_ion3(ion0+1)
                if(ifnr03) then
                  extrasumf(nextrasum) = extrasumf(nextrasum) +
     &              rpowerf*extrasum_scale(nextrasum)*r_ion3(ion0+1)
                  extrasumt(nextrasum) = extrasumt(nextrasum) +
     &              rpowert*extrasum_scale(nextrasum)*r_ion3(ion0+1)
                endif
                if(ifnr13) then
                  do index_f = 1,iatomic_number(ielement)
                    inv_ion_index = inv_ion(index_f+ion0)
                    extrasum_dv(inv_ion_index,nextrasum) = 
     &                extrasum_dv(inv_ion_index,nextrasum) +
     &                rpower_dv(index_f)*rpower*
     &                extrasum_scale(nextrasum)*r_ion3(ion0+1)
                  enddo
                endif
                do iextrasum = 1,nextrasum-2
                  extrasum(iextrasum) = extrasum(iextrasum) +
     &              extrasum_scale(iextrasum)*rpower
                  if(ifnr03) then
                    extrasumf(iextrasum) = extrasumf(iextrasum) +
     &                extrasum_scale(iextrasum)*rpowerf
                    extrasumt(iextrasum) = extrasumt(iextrasum) +
     &                extrasum_scale(iextrasum)*rpowert
                    rpowerf = rpowerf*r_neutral(ielement)
                    rpowert = rpowert*r_neutral(ielement)
                  endif
                  if(ifnr13) then
                    do index_f = 1,iatomic_number(ielement)
                      inv_ion_index = inv_ion(index_f+ion0)
                      extrasum_dv(inv_ion_index,iextrasum) = 
     &                  extrasum_dv(inv_ion_index,iextrasum) +
     &                  rpower_dv(index_f)*
     &                  extrasum_scale(iextrasum)*rpower
                    enddo
                  endif
                  rpower = rpower*r_neutral(ielement)
                enddo
              endif
            endif
          else
C            nuvar(1,index_element) always calculated (see above).
            if(ifexcited.gt.0) then
              if(ifnr03) then
                nuvarf(1,index_element) = 0.d0
                nuvart(1,index_element) = 0.d0
              endif
              if(ifnr13) then
                do index_f = 1, iatomic_number(ielement)
                  nuvar_dv(index_f,1,index_element) = 0.d0
                enddo
              endif
            endif
          endif
          if(index_max.ne.0) then
            sarg0 = bi_ref(ion0+index_max)*tc2-dv(ion0+index_max)
          endif
          nuvar_atomic_number(index_element) =
     &      iatomic_number(ielement)
          nuvar_index_element(index_element) = ielement
          do jndex_f = 1,iatomic_number(ielement)
            if(ielement.eq.1.or.
     &          (fract(jndex_f).gt.0.d0)) then
C              sum_full is the sum of fract values for this element so
C              fract(jndex_f)/sum_full is the fraction of the number
C              density of this element = the fraction of
C              (rho*Navogadro)*epsilon in this particular form.  Thus,
C              vz = number density/(rho*Navogadro) of this particular form
C              of the element, and sum of vz over all forms of the element
C              is equal to epsilon.
              vz = fract(jndex_f)*epsilon/sum_full
              if(ifnr03) then
C                sumf and sumt skip the maximum component
C                expressions below have negligible significance loss.
                if(index_max.eq.jndex_f) then
                  vzf = vz*(fractf(jndex_f)*sum/sum_full-sumf)
                  vzt = vz*(fractt(jndex_f)*sum/sum_full-sumt)
                elseif(index_max.ne.0) then
                  vzf = vz*(fractf(jndex_f) - sumf -
     &              fract(index_max)*fractf(index_max)/sum_full)
                  vzt = vz*(fractt(jndex_f) - sumt -
     &              fract(index_max)*fractt(index_max)/sum_full)
                else
C                  for index_max = 0, sumf and sumt are complete.
                  vzf = vz*(fractf(jndex_f)-sumf)
                  vzt = vz*(fractt(jndex_f)-sumt)
                endif
              endif
              if(ifnr13) then
                do index_f = 1, iatomic_number(ielement)
                  if(index_f.ne.jndex_f) then
                    vz_dv(index_f) = -vz*fract(index_f)/sum_full
                  endif
                enddo
                if(index_max.eq.jndex_f) then
                  vz_dv(jndex_f) = vz*sum/sum_full
                else
                  vz_dv(jndex_f) = vz*
     &              (1.d0 - fract(jndex_f)/sum_full)
                endif
              endif
              if(iffion.and.vz.gt.0.d0) then
                if(ielement.eq.1) then
                  fionh = fionh + vz*
     &              (log(vz) - logqtl_const(nelements+ion0+jndex_f) +
     &              tc2*bi_ref(ion0+jndex_f))
                  if(ifnr03)
     &              fionhf = fionhf + vzf*
     &              (log(vz) - logqtl_const(nelements+ion0+jndex_f) +
     &              tc2*bi_ref(ion0+jndex_f) + 1.d0)
                  if(ifnr13) then
                    do index_f = 1, iatomic_number(ielement)
                      fionh_dv(index_f) = fionh_dv(index_f) +
     &                  vz_dv(index_f)*
     &                  (log(vz) -
     &                  logqtl_const(nelements+ion0+jndex_f) +
     &                  tc2*bi_ref(ion0+jndex_f) + 1.d0)
                    enddo
                  endif
                else
                  fion = fion + vz*
     &              (log(vz) - logqtl_const(nelements+ion0+jndex_f) +
     &              tc2*bi_ref(ion0+jndex_f))
                  if(ifnr03)
     &              fionf = fionf + vzf*
     &              (log(vz) - logqtl_const(nelements+ion0+jndex_f) +
     &              tc2*bi_ref(ion0+jndex_f) + 1.d0)
                  if(ifnr13) then
                    do index_f = 1, iatomic_number(ielement)
                      inv_ion_index = inv_ion(index_f+ion0)
                      fion_dv(inv_ion_index) =
     &                  fion_dv(inv_ion_index) +
     &                  vz_dv(index_f)*
     &                  (log(vz) -
     &                  logqtl_const(nelements+ion0+jndex_f) +
     &                  tc2*bi_ref(ion0+jndex_f) + 1.d0)
                    enddo
                  endif
                endif
              endif
C              in all cases and all indices save vz in nuvar
C              this is needed by excitation_sum (excluding bare ion)
C              and eos_sum_calc and eos_free_calc (including bare ion).
              nuvar(jndex_f+1,index_element) = vz
C              save derivatives as well for excitation and for all but
C              bare ion index range.
              if(ifexcited.gt.0.and.
     &            jndex_f.lt.iatomic_number(ielement)) then
                if(ifnr03) then
                  nuvarf(jndex_f+1,index_element) = vzf
                  nuvart(jndex_f+1,index_element) = vzt
                endif
                if(ifnr13) then
                  do index_f = 1, iatomic_number(ielement)
                    nuvar_dv(index_f,jndex_f+1,index_element) =
     &                vz_dv(index_f)
                  enddo
                endif
              endif
              if(ielement.eq.1) then
C                save special quantities for hydrogen.
C                note, jndex_f is 1.
                sarg = bi_ref(ion0+1)*tc2-dv(ion0+1)
                if(index_max.eq.0) then
C                 _use_regular expression.
                  sh = vz*sarg - epsilon*(log(epsilon) + logsum0)
                else
C                  when h is maximum component of sum_full, then
C                  cancellations occur which cause significance loss
C                  unless you_use_this revised expression.
C                  n.b. -logsum0 = log(sum_full/fract_neut)
C                    = log(fract_max/fract_neut) + log(1+sum/fract_max)
                  sh = epsilon*(-log(epsilon) + ce(ion0+index_max) +
     &              log(1.d0+sum/fract(index_max)) - (sum/sum_full)*
     &              (bi_ref(ion0+index_max)*tc2-dv(ion0+index_max)))
                endif
                sh = sh + constant_sh
                h = vz
                if(ifnr03) then
                  shf = vzf*sarg
                  sht = vzt*sarg
                  hf = vzf
                  ht = vzt
                endif
                if(ifnr13) then
                  inv_ion_index = inv_ion(ion0+1)
                  h_dv(inv_ion_index) = vz_dv(1)
C                 hydrogen contribution to sum0 is zero inside
C                 this routine.
C                 sum0_dv not pre-zeroed before do loop over fract
                  if(ifcsum) then
                    inv_ion_index = inv_ion(ion0+1)
                    sum0_dv(inv_ion_index) = 0.d0
                  endif
                endif
              else
C                save special quanties for helium.
                if(ielement.eq.2.and.jndex_f.eq.1) then
                  hd = vz
                  if(ifnr03) then
                    hdf = vzf
                    hdt = vzt
                  endif
                  if(ifnr13) then
                    inv_ion_index = inv_ion(2)
                    hd_dv(inv_ion_index) = vz_dv(1)
                    inv_ion_index = inv_ion(3)
                    hd_dv(inv_ion_index) = vz_dv(2)
                  endif
                elseif(ielement.eq.2.and.jndex_f.eq.2) then
                  he = vz
                  if(ifnr03) then
                    hef = vzf
                    het = vzt
                  endif
                  if(ifnr13) then
                    inv_ion_index = inv_ion(2)
                    he_dv(inv_ion_index) = vz_dv(1)
                    inv_ion_index = inv_ion(3)
                    he_dv(inv_ion_index) = vz_dv(2)
                  endif
                endif
                ion_index = ion0+jndex_f
                u = u + bi_ref(ion_index)*vz
C                sarg should have dvzero(ielement) subtracted, but
C                this split off to avoid significance loss
C                for s, sf, and st
C                for same reason subtract off sarg of maximum
C                component to be added later.
                if(index_max.eq.0) then
                  sarg = bi_ref(ion_index)*tc2-dv(ion_index)
                elseif(jndex_f.ne.index_max) then
                  sarg = bi_ref(ion_index)*tc2-dv(ion_index) - sarg0
                else
C                  sarg = bi_ref(ion_index)*tc2-dv(ion_index) - sarg0
                  sarg = 0.d0
                endif
C                ignore entropy term due to maximum fraction.  this term
C                will be added later in a way which avoids significance
C                loss.  the derivatives are unaffected by this significance
C                loss so they are done without the complications.
C                n.b. we also ignore *complete* sum over vz*(sarg0-dvzero)
C                which will be added later.
                if(jndex_f.ne.index_max) then
                  s = s + vz*sarg
                endif
                if(ifcsum) then
C                  accumulate number density of positive ions/
C                  (rho/H = rho*Navogadro)
C                  do the following sum analytically to avoid
C                  significance loss. also subtract out constant
C                  component of sum2.
C                  sum0 = sum0 + vz 
C                  accumulate number density of positive ions *
C                  charge^2/(rho/H = rho*Navogadro)
                  sum2 = sum2 +
     &              dble(jndex_f*jndex_f-index_max*index_max)*
     &              sum2_scale*vz
                  if(ifnr03) then
C                    see above comment.
C                    sum0f = sum0f + vzf
C                    sum0t = sum0t + vzt
                    sum2f = sum2f +
     &                dble(jndex_f*jndex_f-index_max*index_max)*
     &                sum2_scale*vzf
                    sum2t = sum2t +
     &                dble(jndex_f*jndex_f-index_max*index_max)*
     &                sum2_scale*vzt
                  endif
                endif
C                accumulate number density of positive charges/
C                (rho/H = rho*Navogadro)
C                add in constant component later except for 
C                index_max = 0 case.
                hne = hne + dble(jndex_f-index_max)*vz
                if(ifnr03) then
                  hnef = hnef + dble(jndex_f-index_max)*vzf
                  hnet = hnet + dble(jndex_f-index_max)*vzt
                endif
                if((ifpl.eq.1).and.
     &              jndex_f.lt.iatomic_number(ielement)) then
C                  Planck-Larkin occupation probability sums
C                  add into sum if not bare nucleus.
                  sumpl0 = sumpl0 + vz*plop(ion_index+1)
                  sumpl1 = sumpl1 + vz*(plop(ion_index+1) +
     &              plopt(ion_index+1))
                  sumpl2 = sumpl2 + vz*plopt(ion_index+1)
                  if(ifnr03) then
                    sumpl0f = sumpl0f + vzf*plop(ion_index+1)
                    sumpl1f = sumpl1f +
     &                vzf*(plop(ion_index+1)+plopt(ion_index+1))
                    sumpl1t = sumpl1t +
     &                vzt*(plop(ion_index+1)+plopt(ion_index+1)) +
     &                vz*(plopt(ion_index+1)+plopt2(ion_index+1))
                  endif
                  if(ifnr13) then
                    do index_f = 1,iatomic_number(ielement)
                      inv_ion_index = inv_ion(index_f+ion0)
                      sumpl0_dv(inv_ion_index) =
     &                  sumpl0_dv(inv_ion_index) +
     &                  vz_dv(index_f)*plop(ion_index+1)
                    enddo
                  endif
                endif
                if(ifpi34) then
                  rpower = dble(jndex_f)**1.5d0
                  extrasum(nextrasum-1) = extrasum(nextrasum-1) +
     &              rpower*extrasum_scale(nextrasum-1)*vz
                  if(ifnr03) then
                    extrasumf(nextrasum-1) = extrasumf(nextrasum-1) +
     &                rpower*extrasum_scale(nextrasum-1)*vzf
                    extrasumt(nextrasum-1) = extrasumt(nextrasum-1) +
     &                rpower*extrasum_scale(nextrasum-1)*vzt
                  endif
                  if(ifnr13) then
                    do index_f = 1,iatomic_number(ielement)
                      inv_ion_index = inv_ion(index_f+ion0)
                      extrasum_dv(inv_ion_index,nextrasum-1) = 
     &                  extrasum_dv(inv_ion_index,nextrasum-1) +
     &                  rpower*
     &                  extrasum_scale(nextrasum-1)*vz_dv(index_f)
                    enddo
                  endif
C                  the ion radii indices are offset by one, and the bare
C                  nucleii should be skipped.
C                  n.b. the if statement means that hydrogen is always
C                  skipped, and eos_calc must correct for that.
                  if(jndex_f.lt.iatomic_number(ielement)) then
                    if(index_maxnb.eq.0) then
                      rpower = r_ion3(ion_index+1)
                    elseif(jndex_f.ne.index_maxnb) then
                      rpower = r_ion3(ion_index+1) -
     &                  r_ion3(ion0+index_maxnb+1)
                    else
C                      offset value.
                      rpower = r_ion3(ion0+index_maxnb+1)
                    endif
                    if(jndex_f.ne.index_maxnb) then
C                      n.b. this if block includes
C                      the case where index_maxnb = 0,
C                      and no offset value is needed or calculated.
                      extrasum(nextrasum) = extrasum(nextrasum) +
     &                  rpower*extrasum_scale(nextrasum)*vz
                      if(ifnr03) then
                        extrasumf(nextrasum) =
     &                    extrasumf(nextrasum) +
     &                    rpower*extrasum_scale(nextrasum)*vzf
                        extrasumt(nextrasum) =
     &                    extrasumt(nextrasum) +
     &                    rpower*extrasum_scale(nextrasum)*vzt
                      endif
                      if(ifnr13) then
                        do index_f = 1,iatomic_number(ielement)
                          inv_ion_index = inv_ion(index_f+ion0)
                          extrasum_dv(inv_ion_index,nextrasum) = 
     &                      extrasum_dv(inv_ion_index,nextrasum) +
     &                      rpower*
     &                      extrasum_scale(nextrasum)*vz_dv(index_f)
                        enddo
                      endif
                    else
C                      offset maximum component is zero.
C                      maximum component index
C                      encountered once (if at all) per element.
C                      for this index add in
C                      sum over all ions of offset.
                      extrasum(nextrasum) = extrasum(nextrasum) +
     &                  rpower*epsilon*sumnb*
     &                  (extrasum_scale(nextrasum)/sum_full)
C                      to obtain derivatives take sumnb/sum_full =
C                      1-(fractneutral+fract(iatomic_number(ielement)))/
C                      sum_full
                      if(ifnr03) then
C                        sumf and sumt skip the maximum component
C                        expressions below have negligible significance loss.
                        if(index_max.eq.iatomic_number(ielement))
     &                      then
                          extrasumf(nextrasum) =
     &                      extrasumf(nextrasum) -
     &                      rpower*epsilon*
     &                      (extrasum_scale(nextrasum)/sum_full)*
     &                      (fract(iatomic_number(ielement))*
     &                      (fractf(iatomic_number(ielement))*
     &                      sum/sum_full-sumf) - fractneutral*
     &                      (sumf + fract(index_max)*
     &                      fractf(index_max)/sum_full))
                          extrasumt(nextrasum) =
     &                      extrasumt(nextrasum) -
     &                      rpower*epsilon*
     &                      (extrasum_scale(nextrasum)/sum_full)*
     &                      (fract(iatomic_number(ielement))*
     &                      (fractt(iatomic_number(ielement))*
     &                      sum/sum_full-sumt) - fractneutral*
     &                      (sumt + fract(index_max)*
     &                      fractt(index_max)/sum_full))
                        else
                          extrasumf(nextrasum) =
     &                      extrasumf(nextrasum) -
     &                      rpower*epsilon*
     &                      (extrasum_scale(nextrasum)/sum_full)*
     &                      (fract(iatomic_number(ielement))*
     &                      fractf(iatomic_number(ielement)) -
     &                      (fractneutral +
     &                      fract(iatomic_number(ielement)))*
     &                      (sumf + fract(index_max)*
     &                      fractf(index_max)/sum_full))
                          extrasumt(nextrasum) =
     &                      extrasumt(nextrasum) -
     &                      rpower*epsilon*
     &                      (extrasum_scale(nextrasum)/sum_full)*
     &                      (fract(iatomic_number(ielement))*
     &                      fractt(iatomic_number(ielement)) -
     &                      (fractneutral +
     &                      fract(iatomic_number(ielement)))*
     &                      (sumt + fract(index_max)*
     &                      fractt(index_max)/sum_full))
                        endif
                      endif
                      if(ifnr13) then
                        do index_f = 1, iatomic_number(ielement)-1
                          inv_ion_index = inv_ion(index_f+ion0)
                          extrasum_dv(inv_ion_index,nextrasum) = 
     &                      extrasum_dv(inv_ion_index,nextrasum) +
     &                      rpower*epsilon*
     &                      (extrasum_scale(nextrasum)/sum_full)*
     &                      (fractneutral +
     &                      fract(iatomic_number(ielement)))*
     &                      (fract(index_f)/sum_full)
                        enddo
                        inv_ion_index =
     &                    inv_ion(iatomic_number(ielement)+ion0)
                        extrasum_dv(inv_ion_index,nextrasum) = 
     &                    extrasum_dv(inv_ion_index,nextrasum) -
     &                    rpower*epsilon*sumnb*
     &                    (extrasum_scale(nextrasum)/sum_full)*
     &                    (fract(iatomic_number(ielement))/sum_full)
                      endif
                    endif
                  endif
                endif
C                naively, this expression doesn't seem right, but have
C                used constancy of sum of vz over all species of helium and
C                metals to derive this expression.
                if(ifnr03) then
C                  skip sarg0 - dvzero(ielement) term which will
C                  be added in later in a way that reduces
C                  significance loss
                  sf = sf + vzf*sarg
                  st = st + vzt*sarg
                endif
                if(ifnr13) then
C                  to avoid significance loss_use_special
C                  pre-summed form for sum0_dv.
                  if(ifcsum) then
                    inv_ion_index = inv_ion(ion0+jndex_f)
                    sum0_dv(inv_ion_index) =
     &                vz*fractneutral*(sum0_scale/sum_full)
                    sum2_dv(inv_ion_index) =
     &                sum2_dv(inv_ion_index) +
     &                dble(index_max*index_max)*
     &                vz*fractneutral*(sum2_scale/sum_full)
                  endif
                  do index_f = 1, iatomic_number(ielement)
                    inv_ion_index = inv_ion(ion0+index_f)
                    if(ifcsum) then
C                      avoid the following commented out expression
C                      because of bad significance loss.
!                      sum0_dv(inv_ion_index) = 
!     &                  sum0_dv(inv_ion_index) + vz_dv(index_f)
                      sum2_dv(inv_ion_index) =
     &                  sum2_dv(inv_ion_index) +
     &                  dble(jndex_f*jndex_f-index_max*index_max)*
     &                  sum2_scale*vz_dv(index_f)
                    endif
                    hne_dv(inv_ion_index) = hne_dv(inv_ion_index) +
     &                dble(jndex_f-index_max)*vz_dv(index_f)
                  enddo
C                _use_sum rule for vz_dv for constant part to
C                 avoid significance loss.
                  inv_ion_index = inv_ion(ion0+jndex_f)
                  hne_dv(inv_ion_index) = hne_dv(inv_ion_index) +
     &              dble(index_max)*vz*fractneutral/sum_full
                endif
              endif
            else
C              n.b. ielement.ne.1.and.fract(jndex_f).le.0.d0
C              save special quantities for helium.
              if(ielement.eq.2.and.jndex_f.eq.1) then
                hd = 0.d0
                if(ifnr03) then
                  hdf = 0.d0
                  hdt = 0.d0
                endif
                if(ifnr13) then
                  inv_ion_index = inv_ion(2)
                  hd_dv(inv_ion_index) = 0.d0
                  inv_ion_index = inv_ion(3)
                  hd_dv(inv_ion_index) = 0.d0
                endif
              elseif(ielement.eq.2.and.jndex_f.eq.2) then
                he = 0.d0
                if(ifnr03) then
                  hef = 0.d0
                  het = 0.d0
                endif
                if(ifnr13) then
                  inv_ion_index = inv_ion(2)
                  he_dv(inv_ion_index) = 0.d0
                  inv_ion_index = inv_ion(3)
                  he_dv(inv_ion_index) = 0.d0
                endif
              endif
              if(ifnr13) then
C                sum0_dv not pre-zeroed before do loop over fract
                if(ifcsum) then
                  inv_ion_index = inv_ion(ion0+jndex_f)
                  sum0_dv(inv_ion_index) = 0.d0
                endif
              endif
              nuvar(jndex_f+1,index_element) = 0.d0
              if(ifexcited.gt.0.and.
     &            jndex_f.lt.iatomic_number(ielement)) then
                if(ifnr03) then
                  nuvarf(jndex_f+1,index_element) = 0.d0
                  nuvart(jndex_f+1,index_element) = 0.d0
                endif
                if(ifnr13) then
                  do index_f = 1, iatomic_number(ielement)
                    nuvar_dv(index_f,jndex_f+1,index_element) = 0.d0
                  enddo
                endif
              endif
            endif
C         end of jndex_f loop
          enddo
C          have used constancy of sum of vz over all ionization species 
C          to derive this expression.
C          n.b. there are no corresponding additions to sf and st, see above.
          if(ielement.gt.1.and.epsilon.gt.0.d0) then
            if(index_max.eq.0) then
C              maximum term not dropped from entropy sum above so
C             _use_regular expression.
              s = s - epsilon*(log(epsilon) + logsum0)
C              correct for dvzero term that is missing from *complete*
C              sum over dv.
              s = s - epsilon*(sum/sum_full)*dvzero(ielement)
              if(ifcsum)
     &          sum0 = sum0 + epsilon*sum*(sum0_scale/sum_full)
            else
C              maximum term was dropped from entropy sum above so
C             _use_regular expression corrected by this term in
C              a way that avoids significance loss.
C              n.b. -logsum0 = log(sum_full/fract_neut)
C                = log(fract_max/fract_neut) + log(1+sum/fract_max)
              s = s + (epsilon*(-log(epsilon) + ce(ion0+index_max) +
     &          log(1.d0+sum/fract(index_max)) - (sum/sum_full)*
     &          (bi_ref(ion0+index_max)*tc2 - dv(ion0+index_max))))
C              correct for dvzero term that is missing from *complete*
C              sum over dv.  the -epsilon*dvzero term has already been
C              analytically subtracted without significance loss
              s = s + epsilon*(fractneutral/sum_full)*dvzero(ielement)
C              correct for sarg0 term that must be added for all vz
C              except maximum term
              s = s + epsilon*((sum - fractneutral)/sum_full)*sarg0
              hne = hne +
     &          epsilon*(1.d0 - fractneutral/sum_full)*dble(index_max)
              if(ifcsum) then
                sum0 = sum0 + epsilon*sum0_scale*
     &            (1.d0 - fractneutral/sum_full)
                sum2 = sum2 + epsilon*sum2_scale*
     &            (1.d0 - fractneutral/sum_full)*
     &            dble(index_max*index_max)
              endif
            endif
            if(ifnr03) then
C              add in sum over vzf and vzt times -dvzero(ielement) term.
C              these sum derivatives are the negative of the vzneutral
C              derivatives which I derive from the vzf and vzt expressions
C              above with vzneutral = fractneutral*epsilon/sum_full
C              replacing vz.
C              note from outside logic that dvzero is non-zero only
C              when fractneutral is small.
              if(index_max.eq.0) then
                sf = sf +
     &            sumf*(fractneutral*epsilon/sum_full)*
     &            (-dvzero(ielement))
                st = st +
     &            sumt*(fractneutral*epsilon/sum_full)*
     &            (-dvzero(ielement))
                if(ifcsum) then
                  sum0f = sum0f + sumf*epsilon*
     &              fractneutral*(sum0_scale/sum_full)
                  sum0t = sum0t + sumt*epsilon*
     &              fractneutral*(sum0_scale/sum_full)
                endif
              else
                sf = sf +
     &            (sumf + fract(index_max)*
     &            fractf(index_max)/sum_full)*
     &            (fractneutral*epsilon/sum_full)*
     &            (sarg0-dvzero(ielement))
                st = st +
     &            (sumt + fract(index_max)*
     &            fractt(index_max)/sum_full)*
     &            (fractneutral*epsilon/sum_full)*
     &            (sarg0-dvzero(ielement))
                hnef = hnef +
     &            (sumf + fract(index_max)*
     &            fractf(index_max)/sum_full)*
     &            (fractneutral*epsilon/sum_full)*dble(index_max)
                hnet = hnet +
     &            (sumt + fract(index_max)*
     &            fractt(index_max)/sum_full)*
     &            (fractneutral*epsilon/sum_full)*dble(index_max)
                if(ifcsum) then
                  sum0f = sum0f + (sumf + fract(index_max)*
     &              fractf(index_max)/sum_full)*epsilon*
     &              fractneutral*(sum0_scale/sum_full)
                  sum0t = sum0t + (sumt + fract(index_max)*
     &              fractt(index_max)/sum_full)*epsilon*
     &              fractneutral*(sum0_scale/sum_full)
                  sum2f = sum2f + (sumf + fract(index_max)*
     &              fractf(index_max)/sum_full)*epsilon*
     &              fractneutral*(sum2_scale/sum_full)*
     &              dble(index_max*index_max)
                  sum2t = sum2t + (sumt + fract(index_max)*
     &              fractt(index_max)/sum_full)*epsilon*
     &              fractneutral*(sum2_scale/sum_full)*
     &              dble(index_max*index_max)
                endif
              endif
            endif
C          end of if(ielement.gt.1.and.epsilon.gt.0.d0) then
          endif
C        end of do index_element = 1, n_partial_elements
        enddo
      endif
      end
