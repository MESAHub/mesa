C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: excitation_pi.f 626 2007-07-19 00:03:57Z airwin $
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
      subroutine excitation_pi(ifexcited,
     &  ifpl, ifpi, ifmodified, ifnr, ifsame_zero_abundances,
     &  inv_aux, iextraoff, naux, inv_ion, ifh2, ifh2plus,
     &  partial_elements, n_partial_elements, ion_end,
     &  tl, izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &  bion, plop, plopt, plopt2,
     &  r_ion3, nion, nions, r_neutral, nelements,
     &  extrasum, extrasumf, extrasumt, nextrasum,
     &  xextrasum, xextrasumf, xextrasumt,
     &  dv, dvf, dvt, ddv_aux)
C       the purpose of this entry is to calculate the chemical potential/kT
C       for each species (neutral and non-bare ions) for the free
C       energy associated with excited states.  Afterward, calculate the
C       corresponding change in the equilibrium constants dv(nions+2).
C       The excitation free energy (per unit volume) model is
C       f = -kT sum n(species) delta ln Z, where
C       delta ln Z = ln(1 + gion*exp(-c2 eion/T)*qstar/(g*plop*piop)),
C       where:
C       gion is the internal partition function of the ion of the species
C       g is the internal partition function of the species
C       qstar is the excited state partition function (relative to the
C         energy of the ion) in the hydrogenic approximation corrected for
C         the Planck-Larkin and MHD occupation probabilities
C         = sum (from n=nmin to nmax) 2 n^2 plop(n)*piop(n)* exp[arg(n)];
C       arg(n) = c2 R iz^2/(n^2 T);
C       R is Rydbergs constant
C       plop = is the ground state Planck-Larkin occupation probability
C       plop(n) = 1 - exp(-arg)*(1 + arg) (see notes for qstar_calc);
C       piop is the ground state MHD occupation probability;
C       piop(n) is the same thing for the excited states (see notes for
C         qstar_calc).
C       input quantities:
C       ifexcited > 0 means_use_excited states (must have Planck-Larkin or
C         ifpi = 3 or 4).
C          0 < ifexcited < 10 means_use_approximation to explicit summation
C         10 < ifexcited < 20 means_use_explicit summation
C         mod(ifexcited,10) = 1 means just apply to hydrogen (without
C           molecules) and helium.
C         mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
C         mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
C       ifpl = 1, worry about Planck-Larkin occupation probability, otherwise
C         set to unity;
C       ifpi = 3 or 4, worry about MHD occupation probability, otherwise
C         set to unity;
C       ifmodified > 0 (or not) affects ifpi_fit
C       ifnr = 0, calculate fl and tl derivatives of output using chain rule
C         and auxf and auxt
C       ifnr = 1, calculate aux derivatives of output with fl, tl fixed
C       ifnr = 2, calculate fl and tl derivatives of output with aux fixed.
C       ifnr = 3 is combination of ifnr = 1 and ifnr = 2.
C       ifsame_zero_abundances = 1 means izlo, izhi, bmin, nmin, nmin_max
C         have not changed since last call.
C       inv_aux(naux) is pointer to actual auxiliary variable index
C       iextraoff (= 7 ) is the index offset between
C         extrasum (*not xextrasum*) and aux.
C       naux is number of active auxiliary variables
C       inv_ion(nions+2) maps ion index to contiguous ion index used in NR 
C         iteration dot products.
C       ifh2 = 0  no h2
C       ifh2 = 1  vdb h2
C       ifh2 = 2  st h2
C       ifh2 = 3  irwin h2 (recommended)
C       ifh2 = 4  pteh h2
C       ifh2plus = 0  no h2plus
C       ifh2plus = 1  st h2plus
C       ifh2plus = 2  irwin h2plus (recommended)
C       partial_elements(n_partial_elements+2) index of elements treated as 
C         partially ionized consistent with ifelement.
C       ion_end(nelements) keeps track of largest ion index for each element
C       tl = log(t).
C       izlo is lowest core charge (1 for H, He, 2 for He+, etc.)
C         for all species.
C       izhi is highest core charge for all species.
C       bmin(izhi) is the minimum bion for all species with a given
C         core charge.
C       nmin(izhi) is the minimum excited principal quantum number for
C         all species with a given core charge.
C       nmin_max(izhi) is the largest minimum excited principal
C         quantum number for all species with a given core charge.
C       nmin_species(nions+2) is the minimum excited principal quantum number
C         organized by species.
C       nmax is maximum principal quantum number included in sum
C         (to be compatible with opal which used nmax = 4 rather than
C         infinity).
C         if(nmax > 300 000) then treated as infinity in qryd_approx.
C         otherwise nmax is meant to be used with qryd_calc only (i.e.,
C         case for mhd approximations not programmed).
C       bion(nions+2) = ionization energy of next higher ion relative to
C         species that is being calculated.  nions+1 refers to H2,
C         nions+2 refers to H2+.
C       plop(nions+2), plopt(nions+2), plopt2(nions+2) = *ln* of ground
C         state Planck-Larkin occupation probabilty and first and
C         second tl derivatives
C       r_ion3(nions+2) is *the cube of the*
C         effective radii (in ion order but going from neutral
C         to next to bare ion for each species) of MDH interaction between
C         all but bare nucleii species and ionized species.  last 2 are
C         H2 and H2+ (note, only need r_ion3 of particular species)
C       nion(nions+2), charge on ion in ion order (must be same order as bi)
C       e.g., for H+, He+, He++, etc.
C       r_neutral(nelements+2) effective radii for MDH neutral-neutral 
C         interactions.  Last two are for for H2 and H2+ (the only ionic
C         species in the MHD model with a non-zero hard-sphere radius).
C       extrasum(nextrasum = 9) weighted sums over n(i)
C         for iextrasum = 1,nextrasum-2,
C         sum is only over neutral species (and H2+) and
C         weight is r_neutral^{iextrasum-1}
C         for iextrasum = nextrasum-1, sum is over all ionized species
C         including bare nucleii, but excluding free electrons,
C         weight is Z^1.5.
C         for iextrasum = nextrasum, sum is over all species
C         excluding bare nucleii and free electrons, the weight is rion^3.
C       extrasumf(nextrasum) = partial of extrasum/partial ln f
C       extrasumt(nextrasum) = partial of extrasum/partial ln t
C       xextrasum(nxextrasum=4) is the *negative* sum nuvar/(1 + qratio)*
C         partial qratio/partial extrasum(k), for k = 1, 2, 3, and nextrasum-1.
C       xextrasumf(nxextrasum) = partial of xextrasum/partial ln f
C       xextrasumt(nxextrasum) = partial of xextrasum/partial ln t
C       modified arrays:
C       dv(nions+2), dvf(nions+2), dvt(nions+2), ddv_aux(nions+2, naux)
C         equilibrium constants and derivatives (depending on ifnr).
      implicit none
      include 'constants.h'
      include 'statistical_weights.h'
C      specific neff, ell, and weight data for helium 1 (supplied by Dappen)
      include 'helium1.h'
      include 'excitation_block.h'
      integer ifexcited, ifpl, ifpi, ifmodified,
     &  ifnr, ifsame_zero_abundances,
     &  naux, inv_aux(naux), iextraoff, ifh2, ifh2plus,
     &  n_partial_elements, partial_elements(n_partial_elements+2),
     &  nelements, ion_end(nelements),
     &  izlo, izhi, nmin(izhi), nmin_max(izhi), nmax,
     &  nions, nion(nions+2), nmin_species(nions+2),
     &  nextrasum, inv_ion(nions+2), inv_ion_index,
C     next line for temporary test
!     &  ion_index, ion_index1,
     &  maxnextrasum, nxextrasum
C      number of xextrasum variables
      parameter(nxextrasum = 4)
      parameter (maxnextrasum = 9)
      double precision tl, bmin(izhi), bion(nions+2),
     &  plop(nions+2), plopt(nions+2), plopt2(nions+2),
     &  r_ion3(nions+2), r_neutral(nelements+2),
     &  extrasum(nextrasum), extrasumf(nextrasum),
     &  extrasumt(nextrasum),
     &  xextrasum(nxextrasum), xextrasumf(nxextrasum),
     &  xextrasumt(nxextrasum), 
     &  dv(nions+2), dvf(nions+2),
     &  dvt(nions+2), ddv_aux(nions+2, naux)
      logical*4 ifmhd_logical, ifpl_logical, ifapprox, 
     &  ifneutral
      integer iz, izqstar, index, 
     &  ielement, ion_start, ion, nmin_s, ifpi_fit
      integer n, ix, min_nmin_max, 
     &  nions_local, nmin_rydberg_he1
C      checked later to agree with calling routine
      parameter(nions_local = 295)
C      minimum value of principal quantum number where approximations are
C      used
      parameter(min_nmin_max = 3)
C      minimum rydberg level for neutral helium (lower levels treated with
C      exact neff, weight, and ell.
C      4 is indistinguishable from 10 for solar comparisons.
      parameter(nmin_rydberg_he1 = 4)
      double precision exparg, exparg_lim, exp_lim, ground_state_lim,
     &  eps_factor, eps_factor_lim,
     &  mu(nions_local+2), muf(nions_local+2), mut(nions_local+2),
     &  mu_aux(2*nxextrasum, nions_local+2),
     &  occ_const
C      corresponds to ~1.d-20
      parameter (exparg_lim = -500.d0)
C      corresponds to ~1.d-300
      parameter (eps_factor_lim = -690.d0)
C      This seems like a reasonable limit on the ln (ground state occupation
C      probability) for when the ground state is worth correcting for
C      excitation.
      parameter(ground_state_lim = -50.d0)
C      multiply qratio and its derivatives by
C      qratio_scale = exp(exparg_shift)
C      in order to get more dynamic range.
C      qratio_us is the unscaled quantity and onepq = 1 + qratio_us.
      double precision exparg_shift, qratio_scale, qratio_us, onepq
      parameter(exparg_shift = 600.d0)
      parameter (occ_const = -4.d0*pi/3.d0)
      double precision 
     &  rfactor, ln_ground_occ
      integer index_x(nx), index_aux
      data index_x /4, 3, 2, 1, 0/
      logical ifnr13, ifnr23
      integer iffirst
      data iffirst/1/
      save
      if(iffirst.eq.1) then
        iffirst = 0
C        Something non-astrophysical
        tl_old_excitation = -1.d30
C        this required so that x_old_excitation doesn't need to be
C        initialized
        ifmhd_logical_old_excitation = .false.
C        The rest of the initial values of the "old_excitation" variables
C        don't matter, but initialize them anyway.
        ifpi_fit_old_excitation = 0
        ifh2_old_excitation = 0
        ifh2plus_old_excitation = 0
        ifpl_logical_old_excitation = .false.
        ifapprox_old_excitation = .false.
        exp_lim = exp(exparg_lim)
        qratio_scale = exp(exparg_shift)
      endif
      ifnr13 = ifnr.eq.1.or.ifnr.eq.3
      ifnr23 = ifnr.eq.2.or.ifnr.eq.3
      c2t = c2*exp(-tl)
      if(ifexcited.le.0) stop 'excitation_pi: invalid ifexcited'
      if(nelements.ne.nelements_stat)
     &  stop 'excitation_pi: invalid nelements'
      if(nions.ne.nions_stat.or.nions.ne.nions_local)
     &  stop 'excitation_pi: invalid nions'
      if(izhi+1.gt.max_izhi) stop 'excitation_pi: izhi too large'
      ifpl_logical = ifpl.eq.1
      ifmhd_logical = ifpi.eq.3.or.ifpi.eq.4
      ifapprox = ifexcited.lt.10
      if(.not.(ifpl_logical.or.ifmhd_logical))
     &  stop 'excitation_pi: invalid combination of ifpl and ifpi'
C      sort out what fitting factors will be applied to pressure ionization.
      if(ifpi.eq.4) then
C        factors to fit Saumon results.
        ifpi_fit = 2
      elseif(ifmodified.gt.0) then
C        factors to fit opal results.
        ifpi_fit = 1
      else
C        unity factors to mimic mdh results as closely as possible.
        ifpi_fit = 0
      endif
C      The following if block is identical to code in excitation_sum
C      (except for routine identification on stop statement message).
C      The idea is all the qh2, qh2plus, and qstar values will be taken
C      from previous calculations (either excitation_sum or here)
C      if nothing has been changed.

C      ifdiff_x_excitation is true if ifmhd_logical is true and
C      the old version not or if both true and the old version has
C      different x values.
C      also initialize x and x_old_excitation if needed.
      if(ifmhd_logical) then
        index_x(nx) = nextrasum-1
        do ix = 1,nx
          x(ix) = extrasum(index_x(ix))
        enddo
        if(ifmhd_logical_old_excitation) then
          ix = 1
          do while(ix.lt.nx.and.x(ix).eq.x_old_excitation(ix))
            ix = ix + 1
          enddo
          ifdiff_x_excitation = x(ix).ne.x_old_excitation(ix)
        else
          ifdiff_x_excitation = .true.
        endif
        do ix = 1,nx
          x_old_excitation(ix) = x(ix)
        enddo
      else
        ifdiff_x_excitation = .false.
      endif
      if(ifdiff_x_excitation.or.
     &    ifpi_fit.ne.ifpi_fit_old_excitation.or.
     &    (ifpl_logical.neqv.ifpl_logical_old_excitation).or.
     &    (ifmhd_logical.neqv.ifmhd_logical_old_excitation).or.
     &    (ifapprox.neqv.ifapprox_old_excitation).or.
     &    tl.ne.tl_old_excitation.or.
     &    ifsame_zero_abundances.ne.1) then
        ifpi_fit_old_excitation = ifpi_fit
        ifpl_logical_old_excitation = ifpl_logical
        ifmhd_logical_old_excitation = ifmhd_logical
        ifapprox_old_excitation = ifapprox
        if(
     &    (
     &    ifh2.ne.ifh2_old_excitation.or.
     &    ifh2plus.ne.ifh2plus_old_excitation.or.
     &    tl.ne.tl_old_excitation
     &    ).and.
     &    (ifh2.gt.0.and.ifh2plus.gt.0.and.mod(ifexcited,10).gt.1))
     &    call molecular_hydrogen(ifh2, ifh2plus, tl,
     &    qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2)
        ifh2_old_excitation = ifh2
        ifh2plus_old_excitation = ifh2plus
        tl_old_excitation = tl
        do iz = izlo, izhi
          if(nmin_max(iz).gt.max_nmin_max)
     &      stop 'excitation_pi: nmin_max too large'
C          eps_factor determines the limit on the excited-state
C          partition function sum.  (Include summands beyond nmin where
C          summand > eps/eps_factor where the error due to the
C          cutoff is roughly eps = 1.d-10.)
C          N.B. the limit on ln(eps_factor) avoids underflows
C          and divide by zeros in the summand test.  Also, for such near-zero
C          eps_factor values, the partition function sum never passes the
C          summand test so is just given by the minimum nmin value.
C          N.B. for ground-state occupation probabilities near unity,
C          eps_factor is always larger than epsarg below so one might
C          be tempted to avoid the call to qstar_calc altogether if
C          eps_factor less than exparg_lim.  However, that logic does
C          not work for small ground-state occupation probabilities so
C          _use_simple floor logic on eps_factor and always call
C          qstar_calc.
          eps_factor = max(eps_factor_lim,-c2t*bmin(iz))
!          if(eps_factor.gt.exparg_lim) then
            eps_factor = exp(eps_factor)
            ifneutral = iz.eq.1
            call qstar_calc(ifpi_fit,
     &        ifpl_logical, ifmhd_logical, ifneutral, ifapprox,
     &        eps_factor, nmin(iz),
     &        max(min_nmin_max,nmin_max(iz)),
     &        nmax,
     &        iz, tl, x, nx,
     &        qstar(1,iz), qstart(1,iz), qstarx(1,1,iz),
     &        qstart2(1,iz), qstartx(1,1,iz), qstarx2(1,1,1,iz))
!C            zero results (including derivatives below) to provide
!C            consistent small value zeroing rather than inconsistent
!C            underflow zeroing.
!            do n = nmin(iz), nmin_max(iz)
!              if(qstar(n,iz).lt.exp_lim) qstar(n,iz) = 0.d0
!            enddo
            if(iz.eq.2.and.ifh2plus.gt.0) then
              ifneutral = .true.
C              special values of qstar and friends calculated including
C              "neutral" occupation probability for H2+.
              call qstar_calc(ifpi_fit,
     &          ifpl_logical, ifmhd_logical, ifneutral, ifapprox,
     &          eps_factor, nmin(iz),
     &          max(min_nmin_max,nmin_max(iz)),
     &          nmax,
     &          iz, tl, x, nx,
     &          qstar(1,izhi+1), qstart(1,izhi+1),
     &          qstarx(1,1,izhi+1),
     &          qstart2(1,izhi+1), qstartx(1,1,izhi+1),
     &          qstarx2(1,1,1,izhi+1))
!C              zero results (including derivatives below) to provide
!C              consistent small value zeroing rather than inconsistent
!C              underflow zeroing.
!              do n = nmin(iz), nmin_max(iz)
!                if(qstar(n,izhi+1).lt.exp_lim) qstar(n,izhi+1) = 0.d0
!              enddo
            endif
!          else
!            do n = nmin(iz), nmin_max(iz)
!              qstar(n,iz) = 0.d0
!            enddo
!            if(iz.eq.2.and.ifh2plus.gt.0) then
!              do n = nmin(iz), nmin_max(iz)
!                qstar(n,izhi+1) = 0.d0
!              enddo
!            endif
!          endif
        enddo
      endif !test for previous calculation of qh2, qh2plus, and qstar
C      do loop over atomic species only
      do index = 1, n_partial_elements
        ielement = partial_elements(index)
        if(ielement.gt.1) then
          ion_start = ion_end(ielement-1) + 1
        else
          ion_start = 1
        endif
C        calculate chemical potential/kt of neutral through next to bare ion
        do ion = ion_start, ion_end(ielement)
          iz = nion(ion)
          exparg = exparg_shift - c2t*bion(ion) - plop(ion)
C          correct for MHD ground state occupation probability.
          if(ifmhd_logical) then
            if(iz.eq.1) then
              ln_ground_occ = occ_const*(
     &          x(1) + r_neutral(ielement)*(
     &          x(2)*3.d0 + r_neutral(ielement)*(
     &          x(3)*3.d0 + r_neutral(ielement)*(
     &          x(4)))))
            else
              ln_ground_occ = 0.d0
            endif
            ln_ground_occ = ln_ground_occ +
     &        occ_const*r_ion3(ion)*x(5)
C            ignore excitation correction if ground state wiped out by
C            pressure ionization in any case.
            if(ln_ground_occ.gt.ground_state_lim) then
              exparg = exparg - ln_ground_occ
            else
C             signal to ignore this species
              exparg = -1000.d0!exparg_lim - 1.d0
            endif
          endif
!          if(exparg.gt.exparg_lim) then
            exparg = exp(exparg)
!          else
!            exparg = 0.d0
!          endif
          if(ielement.le.2.or.mod(ifexcited,10).eq.3) then
            nmin_s = nmin_species(ion)
            if(nmin_s.lt.nmin(iz).or.
     &        nmin_s.gt.nmin_max(iz))
     &        stop 'excitation_pi: invalid nmin or nmin_max'
            if(iz.eq.1) then
C              deal with neutral monatomic species of element
C              ratio of excited to ground state partition functions
              if(ifhe1_special.and.ion.eq.2.and.ifmhd_logical) then
C                special for neutral helium
C                eps_factor determines the limit on the excited-state
C                partition function sum.  (Include summands beyond nmin where
C                summand > eps/eps_factor where the error due to the
C                cutoff is roughly eps = 1.d-10.)
C                N.B. the limit on ln(eps_factor) avoids underflows
C                and divide by zeros in the summand test.  Also, for
C                such near-zero eps_factor values, the partition
C                function sum never passes the
C                summand test so is just given by the minimum nmin value.
C                N.B. for ground-state occupation probabilities near unity,
C                eps_factor is the same as epsarg above so one might
C                be tempted to avoid the call to qstar_calc altogether if
C                eps_factor less than exparg_lim.  However, that logic does
C                not work for small ground-state occupation probabilities so
C                _use_simple floor logic on eps_factor and always call
C                qstar_calc.
                eps_factor = max(eps_factor_lim,-c2t*bion(ion))
!                if(eps_factor.gt.exparg_lim) then
                  eps_factor = exp(eps_factor)
                  rfactor = 1.d0/
     &              (1.d0 + electron_mass/(4.d0*h_mass))
                  call qmhd_calc(ifpi_fit, dble(iqion(ion)),
     &              ifpl_logical, ifapprox, eps_factor,
     &              nmin_rydberg_he1, nmax,
     &              tl, iz, rfactor, nlevels_helium1,
     &              weight_helium1, neff_helium1,
     &              ell_helium1, nx, x,
     &              qmhd_he1, qmhd_he1t, qmhd_he1x,
     &              qmhd_he1t2, qmhd_he1tx, qmhd_he1x2)
!C                  zero results (including derivatives below) to provide
!C                  consistent small value zeroing rather than inconsistent
!C                  underflow zeroing.
!                  if(qmhd_he1.lt.exp_lim) qmhd_he1 = 0.d0
                  mu(ion) = exparg*qmhd_he1/
     &              dble(iqneutral(ielement))
!                else
!                  mu(ion) = 0.d0
!                endif
              else
                mu(ion) = exparg*qstar(nmin_s,iz)*
     &            dble(iqion(ion))/dble(iqneutral(ielement))
              endif
            else
              mu(ion) = exparg*qstar(nmin_s,iz)*
     &          dble(iqion(ion))/dble(iqion(ion-1))
            endif
          else
            mu(ion) = 0.d0
          endif
          if(mu(ion).gt.exp_lim) then
            qratio_us = mu(ion)/qratio_scale
            onepq = 1.d0 + qratio_us
            if(ifnr.eq.0.or.ifnr23) then
C              calculate fl, tl derivatives with extrasum fixed
              muf(ion) = 0.d0
              if(ifhe1_special.and.ion.eq.2.and.ifmhd_logical) then
C                special for neutral helium
                mut(ion) = mu(ion)*(c2t*bion(ion) - plopt(ion) +
     &            qmhd_he1t/qmhd_he1)
              else
                mut(ion) = mu(ion)*(c2t*bion(ion) - plopt(ion) +
     &            qstart(nmin_s,iz)/qstar(nmin_s,iz))
              endif
            endif
            if(ifmhd_logical) then
              if(iz.eq.1) then
                if(ifnr.eq.0) then
                  muf(ion) = muf(ion) - 
     &              occ_const*mu(ion)*r_neutral(ielement)*(
     &              3.d0*extrasumf(3) + r_neutral(ielement)*(
     &              3.d0*extrasumf(2) + r_neutral(ielement)*(
     &              extrasumf(1))))
                  mut(ion) = mut(ion) -
     &              occ_const*mu(ion)*r_neutral(ielement)*(
     &              3.d0*extrasumt(3) + r_neutral(ielement)*(
     &              3.d0*extrasumt(2) + r_neutral(ielement)*(
     &              extrasumt(1))))
                  if(ifhe1_special.and.ion.eq.2) then
C                    special for neutral helium
                    muf(ion) = muf(ion) +
     &                mu(ion)/qmhd_he1*(
     &                qmhd_he1x(4)*extrasumf(1) +
     &                qmhd_he1x(3)*extrasumf(2) +
     &                qmhd_he1x(2)*extrasumf(3))
                    mut(ion) = mut(ion) +
     &                mu(ion)/qmhd_he1*(
     &                qmhd_he1x(4)*extrasumt(1) +
     &                qmhd_he1x(3)*extrasumt(2) +
     &                qmhd_he1x(2)*extrasumt(3))
                  else
                    muf(ion) = muf(ion) +
     &                mu(ion)/qstar(nmin_s,iz)*(
     &                qstarx(4,nmin_s,iz)*extrasumf(1) +
     &                qstarx(3,nmin_s,iz)*extrasumf(2) +
     &                qstarx(2,nmin_s,iz)*extrasumf(3))
                    mut(ion) = mut(ion) +
     &                mu(ion)/qstar(nmin_s,iz)*(
     &                qstarx(4,nmin_s,iz)*extrasumt(1) +
     &                qstarx(3,nmin_s,iz)*extrasumt(2) +
     &                qstarx(2,nmin_s,iz)*extrasumt(3))
                  endif
                elseif(ifnr13) then 
                  mu_aux(3,ion) =
     &              -occ_const*mu(ion)*r_neutral(ielement)*3.d0
                  mu_aux(2,ion) = mu_aux(3,ion)*r_neutral(ielement)
                  mu_aux(1,ion) = mu_aux(2,ion)*r_neutral(ielement)/
     &              3.d0
                  if(ifhe1_special.and.ion.eq.2) then
C                    special for neutral helium
                    mu_aux(3,ion) = mu_aux(3,ion) + mu(ion)*
     &                qmhd_he1x(2)/qmhd_he1
                    mu_aux(2,ion) = mu_aux(2,ion) + mu(ion)*
     &                qmhd_he1x(3)/qmhd_he1
                    mu_aux(1,ion) = mu_aux(1,ion) + mu(ion)*
     &                qmhd_he1x(4)/qmhd_he1
                  else
                    mu_aux(3,ion) = mu_aux(3,ion) + mu(ion)*
     &                qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)
                    mu_aux(2,ion) = mu_aux(2,ion) + mu(ion)*
     &                qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)
                    mu_aux(1,ion) = mu_aux(1,ion) + mu(ion)*
     &                qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)
                  endif
                endif
              endif
              if(ifnr.eq.0) then
                muf(ion) = muf(ion) -occ_const*mu(ion)*
     &            r_ion3(ion)*extrasumf(nextrasum-1)
                mut(ion) = mut(ion) - occ_const*mu(ion)*
     &            r_ion3(ion)*extrasumt(nextrasum-1)
                if(ifhe1_special.and.ion.eq.2) then
C                  special for neutral helium
                  muf(ion) = muf(ion) +
     &              mu(ion)/qmhd_he1*(
     &              qmhd_he1x(5)*extrasumf(nextrasum-1))
                  mut(ion) = mut(ion) +
     &              mu(ion)/qmhd_he1*(
     &              qmhd_he1x(5)*extrasumt(nextrasum-1))
                else
                  muf(ion) = muf(ion) +
     &              mu(ion)/qstar(nmin_s,iz)*(
     &              qstarx(5,nmin_s,iz)*extrasumf(nextrasum-1))
                  mut(ion) = mut(ion) +
     &              mu(ion)/qstar(nmin_s,iz)*(
     &              qstarx(5,nmin_s,iz)*extrasumt(nextrasum-1))
                endif
              elseif(ifnr13) then 
                mu_aux(4,ion) = -occ_const*mu(ion)*r_ion3(ion)
                if(ifhe1_special.and.ion.eq.2) then
C                  special for neutral helium
                  mu_aux(4,ion) = mu_aux(4,ion) +
     &              mu(ion)*qmhd_he1x(5)/qmhd_he1
                else
                  mu_aux(4,ion) = mu_aux(4,ion) +
     &              mu(ion)*qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)
                endif
              endif
            endif
            if(ifnr.eq.0.or.ifnr23) then
C              transform to chemical potential/kt
              muf(ion) = -muf(ion)/onepq
              mut(ion) = -mut(ion)/onepq
            endif
            if(ifnr13.and.ifmhd_logical) then
              if(iz.eq.1) then
                mu_aux(1,ion) = -mu_aux(1,ion)/onepq
                mu_aux(2,ion) = -mu_aux(2,ion)/onepq
                mu_aux(3,ion) = -mu_aux(3,ion)/onepq
              endif
              mu_aux(4,ion) = -mu_aux(4,ion)/onepq
            endif
            if(qratio_us.gt.1.d-3) then
              mu(ion) = -qratio_scale*log(onepq)
            else
C              alternating series so relative error is less than first
C              missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
              mu(ion) = -mu(ion)*
     &          (1.d0      - qratio_us*
     &          (1.d0/2.d0 - qratio_us*
     &          (1.d0/3.d0 - qratio_us*
     &          (1.d0/4.d0 - qratio_us*
     &          (1.d0/5.d0)))))
            endif
          else
            mu(ion) = 0.d0
            if(ifnr.eq.0.or.ifnr23) then
              muf(ion) = 0.d0
              mut(ion) = 0.d0
            endif
            if(ifnr13.and.ifmhd_logical) then
              if(iz.eq.1) then
                mu_aux(1,ion) = 0.d0
                mu_aux(2,ion) = 0.d0
                mu_aux(3,ion) = 0.d0
              endif
              mu_aux(4,ion) = 0.d0
            endif
          endif !if(mu(ion).gt.exp_lim) then
          if(ifmhd_logical) then
            if(iz.eq.1) then
              mu(ion) = mu(ion) + qratio_scale*(
     &          xextrasum(1) +
     &          xextrasum(2)*r_neutral(ielement) +
     &          xextrasum(3)*r_neutral(ielement)*r_neutral(ielement))
              if(ifnr.eq.0) then
                muf(ion) = muf(ion) + qratio_scale*(
     &            xextrasumf(1) +
     &            xextrasumf(2)*r_neutral(ielement) +
     &            xextrasumf(3)*r_neutral(ielement)*
     &            r_neutral(ielement))
                mut(ion) = mut(ion) + qratio_scale*(
     &            xextrasumt(1) +
     &            xextrasumt(2)*r_neutral(ielement) +
     &            xextrasumt(3)*r_neutral(ielement)*
     &            r_neutral(ielement))
              elseif(ifnr13) then
                mu_aux(5,ion) = qratio_scale
                mu_aux(6,ion) = qratio_scale*r_neutral(ielement)
                mu_aux(7,ion) = qratio_scale*r_neutral(ielement)*
     &            r_neutral(ielement)
              endif
            else
              mu(ion) = mu(ion) + qratio_scale*
     &          xextrasum(4)*dble(iz-1)**1.5d0
              if(ifnr.eq.0) then
                muf(ion) = muf(ion) + qratio_scale*
     &            xextrasumf(4)*dble(iz-1)**1.5d0
                mut(ion) = mut(ion) + qratio_scale*
     &            xextrasumt(4)*dble(iz-1)**1.5d0
              elseif(ifnr13) then
                mu_aux(8,ion) = qratio_scale*dble(iz-1)**1.5d0
              endif
            endif
          endif
C          temporary test of derivatives
!          dv(ion) = mu(ion)
!          inv_ion_index = inv_ion(ion)
!          if(iz.eq.1) then
!            index_aux = inv_aux(iextraoff+1)
!            ddv_aux(inv_ion_index,index_aux) = mu_aux(1,ion)
!            index_aux = inv_aux(iextraoff+2)
!            ddv_aux(inv_ion_index,index_aux) = mu_aux(2,ion)
!            index_aux = inv_aux(iextraoff+3)
!            ddv_aux(inv_ion_index,index_aux) = mu_aux(3,ion)
!          else
!            index_aux = inv_aux(iextraoff+1)
!            ddv_aux(inv_ion_index,index_aux) = 0.d0
!            index_aux = inv_aux(iextraoff+2)
!            ddv_aux(inv_ion_index,index_aux) = 0.d0
!            index_aux = inv_aux(iextraoff+3)
!            ddv_aux(inv_ion_index,index_aux) = 0.d0
!          endif
!          index_aux = inv_aux(iextraoff+nextrasum-1)
!          ddv_aux(inv_ion_index,index_aux) = mu_aux(4,ion)
        enddo
!        return
C        calculate change in equilibrium constant of monatomic ions
C        relative to neutral monatomic.
C        n.b. ion now refers to equilibrium constant of h+, he+, he++, etc.,
C        i.e. offset by 1 from previous ion meaning.
        do ion = ion_start, ion_end(ielement)
          if(ion.lt.ion_end(ielement)) then
            dv(ion) = dv(ion) + (mu(ion_start) - mu(ion+1))/
     &        qratio_scale
            if(ifnr.eq.0.or.ifnr23) then
              dvf(ion) = dvf(ion) + (muf(ion_start) - muf(ion+1))/
     &          qratio_scale
              dvt(ion) = dvt(ion) + (mut(ion_start) - mut(ion+1))/
     &          qratio_scale
            endif
            if(ifnr13.and.ifmhd_logical) then
              inv_ion_index = inv_ion(ion)
              index_aux = inv_aux(iextraoff+1)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(1,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+2)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(2,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+3)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(3,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+nextrasum-1)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          (mu_aux(4,ion_start) - mu_aux(4,ion+1))/
     &          qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+1)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(5,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+2)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(6,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+3)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(7,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+4)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) -
     &          mu_aux(8,ion+1)/qratio_scale
            endif
          else
            dv(ion) = dv(ion) + mu(ion_start)/qratio_scale
            if(ifmhd_logical) then
C              n.b. bare ion has a chemical potential of xextrasum(4)*iz^{3/2}
              iz = nion(ion)
              dv(ion) = dv(ion) - xextrasum(4)*dble(iz)**1.5d0
              if(ifnr.eq.0) then
                dvf(ion) = dvf(ion) - xextrasumf(4)*dble(iz)**1.5d0
                dvt(ion) = dvt(ion) - xextrasumt(4)*dble(iz)**1.5d0
              endif
            endif
            if(ifnr.eq.0.or.ifnr23) then
              dvf(ion) = dvf(ion) + muf(ion_start)/qratio_scale
              dvt(ion) = dvt(ion) + mut(ion_start)/qratio_scale
            endif
            if(ifnr13.and.ifmhd_logical) then
              inv_ion_index = inv_ion(ion)
              index_aux = inv_aux(iextraoff+1)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(1,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+2)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(2,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+3)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(3,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+nextrasum-1)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(4,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+1)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(5,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+2)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(6,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+3)
              ddv_aux(inv_ion_index,index_aux) =
     &          ddv_aux(inv_ion_index,index_aux) +
     &          mu_aux(7,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+4)
              ddv_aux(inv_ion_index,index_aux) = 
     &          ddv_aux(inv_ion_index,index_aux) -
     &          dble(iz)**1.5d0
            endif
          endif
        enddo
      enddo
      if(partial_elements(1).eq.1) then
C        do all the molecular stuff below only if hydrogen has a non-zero
C        abundance
C        now do H2 equilibrium constant change due to neutral monatomic H.
        if(ifh2.gt.0) then
C          calculate change in H2 equilibrium constant 
C          relative to neutral monatomic H
          ion = nions+1
          ion_start = 1
          dv(ion) = dv(ion) + 2.d0*mu(ion_start)/qratio_scale
          if(ifnr.eq.0.or.ifnr23) then
            dvf(ion) = dvf(ion) + 2.d0*muf(ion_start)/qratio_scale
            dvt(ion) = dvt(ion) + 2.d0*mut(ion_start)/qratio_scale
          endif
          if(ifnr13.and.ifmhd_logical) then
            inv_ion_index = inv_ion(ion)
            index_aux = inv_aux(iextraoff+1)
            ddv_aux(inv_ion_index,index_aux) =
     &        ddv_aux(inv_ion_index,index_aux) +
     &        2.d0*mu_aux(1,ion_start)/qratio_scale
            index_aux = inv_aux(iextraoff+2)
            ddv_aux(inv_ion_index,index_aux) =
     &        ddv_aux(inv_ion_index,index_aux) +
     &        2.d0*mu_aux(2,ion_start)/qratio_scale
            index_aux = inv_aux(iextraoff+3)
            ddv_aux(inv_ion_index,index_aux) =
     &        ddv_aux(inv_ion_index,index_aux) +
     &        2.d0*mu_aux(3,ion_start)/qratio_scale
            index_aux = inv_aux(iextraoff+nextrasum-1)
            ddv_aux(inv_ion_index,index_aux) =
     &        ddv_aux(inv_ion_index,index_aux) +
     &        2.d0*mu_aux(4,ion_start)/qratio_scale
            index_aux = inv_aux(iextraoff+maxnextrasum+1)
            ddv_aux(inv_ion_index,index_aux) =
     &        ddv_aux(inv_ion_index,index_aux) +
     &        2.d0*mu_aux(5,ion_start)/qratio_scale
            index_aux = inv_aux(iextraoff+maxnextrasum+2)
            ddv_aux(inv_ion_index,index_aux) =
     &        ddv_aux(inv_ion_index,index_aux) +
     &        2.d0*mu_aux(6,ion_start)/qratio_scale
            index_aux = inv_aux(iextraoff+maxnextrasum+3)
            ddv_aux(inv_ion_index,index_aux) =
     &        ddv_aux(inv_ion_index,index_aux) +
     &        2.d0*mu_aux(7,ion_start)/qratio_scale
          endif
        endif
C        now do H2 and H2plus chemical potentials/kT
C        n.b. ifh2, ifh2plus, and ifexcited control whether *xextrasum* has
C        a molecular component.  However, *only* ifh2 and ifh2plus control
C        whether *extrasum* has a molecular component.
        if((ifh2.gt.0.and.ifh2plus.gt.0.and.
     &      mod(ifexcited,10).gt.1).or.ifmhd_logical) then
          do iz = 1,2
            if((iz.eq.1.and.ifh2.gt.0).or.
     &          (iz.eq.2.and.ifh2plus.gt.0)) then
              if(iz.eq.1) then
                ielement = nelements + 1
                ion = nions+1
                izqstar = iz
              else
                ielement = nelements + 2
                ion = nions+2
C                special values of qstar and friends calculated including
C                "neutral" occupation probability.
                izqstar = izhi+1
              endif
              if(ifh2.gt.0.and.ifh2plus.gt.0.and.
     &            mod(ifexcited,10).gt.1) then
                if(iz.eq.1) then
                  exparg = exparg_shift -
     &              c2t*bion(ion) - plop(ion) + qh2plus - qh2
                else
C                  n.b. excited state core is two protons with unity
C                  statistical weight (following usual convention that
C                  nuclear spin statistical weights are divided out).
                  exparg = exparg_shift -
     &              c2t*bion(ion) - plop(ion) - qh2plus
                endif
C                correct for MHD ground state occupation probability.
                if(ifmhd_logical) then
                  ln_ground_occ = occ_const*(
     &              x(1) + r_neutral(ielement)*(
     &              x(2)*3.d0 + r_neutral(ielement)*(
     &              x(3)*3.d0 + r_neutral(ielement)*(
     &              x(4)))))
                  ln_ground_occ = ln_ground_occ +
     &              occ_const*r_ion3(ion)*x(5)
C                  ignore excitation correction if ground state wiped out by
C                  pressure ionization in any case.
                  if(ln_ground_occ.gt.ground_state_lim) then
                    exparg = exparg - ln_ground_occ
                  else
C                    signal to ignore this species
                    exparg = -1000.d0!exparg_lim - 1.d0
                  endif
                endif
!                if(exparg.gt.exparg_lim) then
                  exparg = exp(exparg)
!                else
!                  exparg = 0.d0
!                endif
                nmin_s = nmin_species(ion)
                if(nmin_s.lt.nmin(iz).or.
     &            nmin_s.gt.nmin_max(iz))
     &            stop 'excitation_pi: invalid nmin or nmin_max'
                mu(ion) = exparg*qstar(nmin_s,izqstar)
              else
                mu(ion) = 0.d0
              endif
              if(mu(ion).gt.exp_lim) then
                qratio_us = mu(ion)/qratio_scale
                onepq = 1.d0 + qratio_us
                if(ifnr.eq.0.or.ifnr23) then
C                  calculate fl, tl derivatives with extrasum fixed
                  muf(ion) = 0.d0
                  mut(ion) = mu(ion)*(c2t*bion(ion) - plopt(ion) +
     &              qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar))
                  if(iz.eq.1) then
                    mut(ion) = mut(ion) + mu(ion)*(qh2plust - qh2t)
                  else
                    mut(ion) = mut(ion) + mu(ion)*(-qh2plust)
                  endif
                endif
                if(ifmhd_logical) then
                  if(ifnr.eq.0) then
                    muf(ion) = muf(ion) -
     &                occ_const*mu(ion)*r_neutral(ielement)*(
     &                3.d0*extrasumf(3) + r_neutral(ielement)*(
     &                3.d0*extrasumf(2) + r_neutral(ielement)*(
     &                extrasumf(1)))) +
     &                mu(ion)/qstar(nmin_s,izqstar)*(
     &                qstarx(4,nmin_s,izqstar)*extrasumf(1) +
     &                qstarx(3,nmin_s,izqstar)*extrasumf(2) +
     &                qstarx(2,nmin_s,izqstar)*extrasumf(3))
                    mut(ion) = mut(ion) -
     &                occ_const*mu(ion)*r_neutral(ielement)*(
     &                3.d0*extrasumt(3) + r_neutral(ielement)*(
     &                3.d0*extrasumt(2) + r_neutral(ielement)*(
     &                extrasumt(1)))) +
     &                mu(ion)/qstar(nmin_s,izqstar)*(
     &                qstarx(4,nmin_s,izqstar)*extrasumt(1) +
     &                qstarx(3,nmin_s,izqstar)*extrasumt(2) +
     &                qstarx(2,nmin_s,izqstar)*extrasumt(3))
                  elseif(ifnr13) then 
                    mu_aux(3,ion) =
     &                -occ_const*mu(ion)*r_neutral(ielement)*3.d0
                    mu_aux(2,ion) = mu_aux(3,ion)*
     &                r_neutral(ielement)
                    mu_aux(1,ion) = mu_aux(2,ion)*
     &                r_neutral(ielement)/3.d0
                    mu_aux(3,ion) = mu_aux(3,ion) + mu(ion)*
     &                qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)
                    mu_aux(2,ion) = mu_aux(2,ion) + mu(ion)*
     &                qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)
                    mu_aux(1,ion) = mu_aux(1,ion) + mu(ion)*
     &                qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)
                  endif
                  if(ifnr.eq.0) then
                    muf(ion) = muf(ion) -occ_const*mu(ion)*
     &                r_ion3(ion)*extrasumf(nextrasum-1) +
     &                mu(ion)/qstar(nmin_s,izqstar)*(
     &                qstarx(5,nmin_s,izqstar)*
     &                extrasumf(nextrasum-1))
                    mut(ion) = mut(ion) - occ_const*mu(ion)*
     &                r_ion3(ion)*extrasumt(nextrasum-1) +
     &                mu(ion)/qstar(nmin_s,izqstar)*(
     &                qstarx(5,nmin_s,izqstar)*
     &                extrasumt(nextrasum-1))
                  elseif(ifnr13) then 
                    mu_aux(4,ion) = -occ_const*mu(ion)*r_ion3(ion) +
     &                mu(ion)*qstarx(5,nmin_s,izqstar)/
     &                qstar(nmin_s,izqstar)
                  endif
                endif
                if(ifnr.eq.0.or.ifnr23) then
C                  transform to chemical potential/kt
                  muf(ion) = -muf(ion)/onepq
                  mut(ion) = -mut(ion)/onepq
                endif
                if(ifnr13.and.ifmhd_logical) then
                  mu_aux(1,ion) = -mu_aux(1,ion)/onepq
                  mu_aux(2,ion) = -mu_aux(2,ion)/onepq
                  mu_aux(3,ion) = -mu_aux(3,ion)/onepq
                  mu_aux(4,ion) = -mu_aux(4,ion)/onepq
                endif
                if(qratio_us.gt.1.d-3) then
                  mu(ion) = -qratio_scale*log(onepq)
                else
C                  alternating series so relative error is less than first
C                  missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
                  mu(ion) = -mu(ion)*
     &              (1.d0      - qratio_us*
     &              (1.d0/2.d0 - qratio_us*
     &              (1.d0/3.d0 - qratio_us*
     &              (1.d0/4.d0 - qratio_us*
     &              (1.d0/5.d0)))))
                endif
              else
                mu(ion) = 0.d0
                if(ifnr.eq.0.or.ifnr23) then
                  muf(ion) = 0.d0
                  mut(ion) = 0.d0
                endif
                if(ifnr13.and.ifmhd_logical) then
                  if(iz.eq.1) then
                    mu_aux(1,ion) = 0.d0
                    mu_aux(2,ion) = 0.d0
                    mu_aux(3,ion) = 0.d0
                  endif
                  mu_aux(4,ion) = 0.d0
                endif !if(mu(ion).gt.exp_lim) then
              endif
              if(ifmhd_logical) then
                mu(ion) = mu(ion) + qratio_scale*(
     &            xextrasum(1) +
     &            xextrasum(2)*r_neutral(ielement) +
     &            xextrasum(3)*r_neutral(ielement)*
     &            r_neutral(ielement))
                if(ifnr.eq.0) then
                  muf(ion) = muf(ion) + qratio_scale*(
     &              xextrasumf(1) +
     &              xextrasumf(2)*r_neutral(ielement) +
     &              xextrasumf(3)*r_neutral(ielement)*
     &              r_neutral(ielement))
                  mut(ion) = mut(ion) + qratio_scale*(
     &              xextrasumt(1) +
     &              xextrasumt(2)*r_neutral(ielement) +
     &              xextrasumt(3)*r_neutral(ielement)*
     &              r_neutral(ielement))
                elseif(ifnr13) then
                  mu_aux(5,ion) = qratio_scale
                  mu_aux(6,ion) = qratio_scale*r_neutral(ielement)
                  mu_aux(7,ion) = qratio_scale*r_neutral(ielement)*
     &              r_neutral(ielement)
                endif
                if(iz.eq.2) then
                  mu(ion) = mu(ion) + qratio_scale*xextrasum(4)*
     &              dble(iz-1)**1.5d0
                  if(ifnr.eq.0) then
                    muf(ion) = muf(ion) + qratio_scale*xextrasumf(4)*
     &                dble(iz-1)**1.5d0
                    mut(ion) = mut(ion) + qratio_scale*xextrasumt(4)*
     &                dble(iz-1)**1.5d0
                  elseif(ifnr13) then
                    mu_aux(8,ion) = qratio_scale*dble(iz-1)**1.5d0
                  endif
                endif
              endif
              if(iz.eq.1) then
C                calculate change in H2 equilibrium constant 
C                due to chemical potential of H2 (effect of chemical potential
C                of neutral H already applied).
                dv(ion) = dv(ion) - mu(ion)/qratio_scale
                if(ifnr.eq.0.or.ifnr23) then
                  dvf(ion) = dvf(ion) - muf(ion)/qratio_scale
                  dvt(ion) = dvt(ion) - mut(ion)/qratio_scale
                endif
                if(ifnr13.and.ifmhd_logical) then
                  inv_ion_index = inv_ion(ion)
                  index_aux = inv_aux(iextraoff+1)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) -
     &              mu_aux(1,ion)/qratio_scale
                  index_aux = inv_aux(iextraoff+2)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) -
     &              mu_aux(2,ion)/qratio_scale
                  index_aux = inv_aux(iextraoff+3)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) -
     &              mu_aux(3,ion)/qratio_scale
                  index_aux = inv_aux(iextraoff+nextrasum-1)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) -
     &              mu_aux(4,ion)/qratio_scale
                  index_aux = inv_aux(iextraoff+maxnextrasum+1)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) -
     &              mu_aux(5,ion)/qratio_scale
                  index_aux = inv_aux(iextraoff+maxnextrasum+2)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) -
     &              mu_aux(6,ion)/qratio_scale
                  index_aux = inv_aux(iextraoff+maxnextrasum+3)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) -
     &              mu_aux(7,ion)/qratio_scale
                endif
              else
C                calculate change in H2+ equilibrium constant 
C                relative to H2 and e-.
                ion_start = nions + 1
                dv(ion) = dv(ion) + (mu(ion_start) - mu(ion))/
     &            qratio_scale
                if(ifnr.eq.0.or.ifnr23) then
                  dvf(ion) = dvf(ion) + (muf(ion_start) - muf(ion))/
     &              qratio_scale
                  dvt(ion) = dvt(ion) + (mut(ion_start) - mut(ion))/
     &              qratio_scale
                endif
                if(ifnr13.and.ifmhd_logical) then
                  inv_ion_index = inv_ion(ion)
                  index_aux = inv_aux(iextraoff+1)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) +
     &              (mu_aux(1,ion_start) - mu_aux(1,ion))/
     &              qratio_scale
                  index_aux = inv_aux(iextraoff+2)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) +
     &              (mu_aux(2,ion_start) - mu_aux(2,ion))/
     &              qratio_scale
                  index_aux = inv_aux(iextraoff+3)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) +
     &              (mu_aux(3,ion_start) - mu_aux(3,ion))/
     &              qratio_scale
                  index_aux = inv_aux(iextraoff+nextrasum-1)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) +
     &              (mu_aux(4,ion_start) - mu_aux(4,ion))/
     &              qratio_scale
                  index_aux = inv_aux(iextraoff+maxnextrasum+1)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) +
     &              (mu_aux(5,ion_start) - mu_aux(5,ion))/
     &              qratio_scale
                  index_aux = inv_aux(iextraoff+maxnextrasum+2)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) +
     &              (mu_aux(6,ion_start) - mu_aux(6,ion))/
     &              qratio_scale
                  index_aux = inv_aux(iextraoff+maxnextrasum+3)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) +
     &              (mu_aux(7,ion_start) - mu_aux(7,ion))/
     &              qratio_scale
                  index_aux = inv_aux(iextraoff+maxnextrasum+4)
                  ddv_aux(inv_ion_index,index_aux) =
     &              ddv_aux(inv_ion_index,index_aux) -
     &              mu_aux(8,ion)/qratio_scale
                endif
              endif
            endif
          enddo
        endif
      endif
      end
