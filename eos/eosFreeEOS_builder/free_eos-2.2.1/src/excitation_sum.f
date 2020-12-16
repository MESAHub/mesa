C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: excitation_sum.f 827 2008-06-29 01:32:24Z airwin $
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
      subroutine excitation_sum(ifexcited, ifsame_zero_abundances,
     &  ifpl, ifpi, ifmodified, ifnr, inv_ion, ifh2, ifh2plus,
     &  partial_elements, n_partial_elements, ion_end,
     &  tl, izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &  bion, plop, plopt, plopt2,
     &  r_ion3, nion, nions, r_neutral, nelements,
     &  extrasum, extrasumf, extrasumt, extrasum_dv, nextrasum,
     &  max_index, nug, nugf, nugt, nug_dv,
     &  nuh2, nuh2f, nuh2t, nuh2_dv,
     &  nuh2plus, nuh2plusf, nuh2plust, nuh2plus_dv,
     &  xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
C       the purpose of this subroutine is to calculate auxiliary variables
C       and their derivatives that are associated with the free energy of
C       excited states.
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
C         mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules)
C           and helium.
C         mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
C         mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
C       ifpl = 1, calculate Planck-Larkin occupation probability, otherwise
C         set to unity;
C       ifpi = 3 or 4, calculate MHD occupation probability, otherwise
C         set to unity;
C       ifmodified > 0 (or not) affects ifpi_fit
C       ifnr = 0
C         calculate derivatives of psum and usum (delivered through common
C         block) and xextrasum wrt f, t with no other variables fixed,
C         i.e.,_use_input f and t derivatives and chain rule.
C       N.B. chain of calling routines depend on the assertion that ifnr = 0
C       means no *_dv variables should be read or written.
C       ifnr = 1
C         calculate derivatives of xextrasum wrt dv with f, t fixed.
C         n.b. this variable is subset of the list for ifnr = 0
C         because we only need dv derivatives for variables used to calculate
C         (output) auxiliary variables.
C       ifnr = 2 (should not occur)
C       ifnr = 3 same as combination of ifnr = 0 and ifnr = 1.
C         Note this is quite different in detail from ifnr = 3 interpretation
C         for many other routines, but the general motivation is the same
C         for all ifnr = 3 results; calculate both f and t derivatives and
C         other derivatives.
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
C       bmin(izhi) is the minimum bion for all species with a given core charge.
C       nmin(izhi) is the minimum excited principal quantum number for
C         all species with a given core charge.
C       nmin_max(izhi) is the largest minimum excited principal
C         quantum number for all species with a given core charge.
C       nmin_species(nions+2) is the minimum excited principal quantum number
C         organized by species.
C       nmax is maximum principal quantum number included in sum
C         (to be compatible with opal which used nmax = 4 rather than infinity).
C         if(nmax > 300 000) then treated as infinity in qryd_approx.
C         otherwise nmax is meant to be used with qryd_calc only (i.e.,
C         case for mhd approximations not programmed).
C       bion(nions+2) = ionization energy of next higher ion relative to species
C         that is being calculated.  nions+1 refers to H2, nions+2 refers to H2+.
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
C         interactions.   Last two are H2 and H2+ (the only ionic
C         species in the MHD model with a non-zero hard-sphere radius).
C       extrasum(nextrasum = 9) weighted sums over n(i)
C         for iextrasum = 1,nextrasum-2,
C         sum is only over neutral species (and H2+) and
C         weight is r_neutral^{iextrasum-1}
C         for iextrasum = nextrasum-1, sum is over all ionized species including
C         bare nucleii, but excluding free electrons, weight is Z^1.5.
C         for iextrasum = nextrasum, sum is over all species excluding bare nucleii
C         and free electrons, the weight is rion^3.
C       extrasumf(nextrasum) = partial of extrasum/partial ln f
C       extrasumt(nextrasum) = partial of extrasum/partial ln t
C       extrasum_dv(nions+2,nextrasum) = partial derivatives wrt dv
C       max_index
C       nug
C       nuh2
C       nuh2plus
C       xextrasum(4) is the *negative* sum nuvar/(1 + qratio)*
C       partial qratio/partial extrasum(k), k = 1, 2, 3, and nextrasum-1.
      implicit none
      include 'aux_scale.h'
      include 'constants.h'
      include 'nuvar.h'
      include 'statistical_weights.h'
C      specific neff, ell, and weight data for helium 1 (supplied by Dappen)
      include 'helium1.h'
      include 'excitation_block.h'
      integer ifexcited, ifsame_zero_abundances,
     &  ifpl, ifpi, ifmodified,
     &  ifnr, ifh2, ifh2plus,
     &  n_partial_elements, partial_elements(n_partial_elements+2),
     &  nelements, ion_end(nelements),
     &  izlo, izhi, nmin(izhi), nmin_max(izhi), nmax,
     &  nions, nion(nions+2), nmin_species(nions+2),
     &  nextrasum, inv_ion(nions+2), inv_ion_index,
     &  ion_index, ion_index1, nxextrasum
C      number of xextrasum variables
      parameter(nxextrasum = 4)
      double precision tl, bmin(izhi), bion(nions+2),
     &  plop(nions+2), plopt(nions+2), plopt2(nions+2),
     &  r_ion3(nions+2), r_neutral(nelements+2),
     &  extrasum(nextrasum), extrasumf(nextrasum),
     &  extrasumt(nextrasum), extrasum_dv(nions+2,nextrasum),
     &  nug, nugf, nugt, nug_dv(nions+2),
     &  nuh2, nuh2f, nuh2t, nuh2_dv(nions+2),
     &  nuh2plus, nuh2plusf, nuh2plust, nuh2plus_dv(nions+2),
     &  xextrasum(nxextrasum), xextrasumf(nxextrasum),
     &  xextrasumt(nxextrasum), xextrasum_dv(nions+2,nxextrasum)
      logical*4 ifpl_logical, ifmhd_logical, ifneutral, ifapprox
      integer iz, izqstar, index, max_index,
     &  ielement, ion_start, ion, nmin_s, ifpi_fit
      integer n, ix, min_nmin_max,
     &  nions_local, nmin_rydberg_he1
C      checked later to agree with calling routine
      parameter(nions_local = 295)
C      minimum value of principal quantum number where approximations are
C        used
      parameter(min_nmin_max = 3)
C      minimum rydberg level for neutral helium (lower levels treated with
C      exact neff, weight, and ell.
C      4 is indistinguishable from 10 for solar comparisons.
      parameter(nmin_rydberg_he1 = 4)
      double precision exparg, exparg_lim, exp_lim, ground_state_lim,
     &  eps_factor, eps_factor_lim, rfactor,
     &  qratio, qratiot, qratio_dx(nxextrasum),
     &  qratiot2, qratiot_dx(nxextrasum),
     &  qratio_dx2(nxextrasum,nxextrasum),
     &  qratiov, qratiovt, qratiov_dx(nxextrasum),
     &  occ_const, 
     &  dqratiof, dqratiotf, dqratio_dxf(nxextrasum),
     &  dqratiot, dqratiott, dqratio_dxt(nxextrasum),
     &  dqratio_dv(nions_local+2),
     &  dqratio_dxdv(nions_local+2,nxextrasum),
     &  dqratiovf, dqratiovt, dqratiov_dv(nions_local+2)
C      corresponds to ~1.d-20
      parameter (exparg_lim = -500.d0)
C      corresponds to ~1.d-300
      parameter (eps_factor_lim = -690.d0)
C      This seems like a reasonable limit on the ln (ground state occupation
C      probability) for when the ground state is worth correcting for
C      excitation.
      parameter(ground_state_lim = -50.d0)
      parameter (occ_const = -4.d0*pi/3.d0)
C      multiply qratio and its derivatives by
C      qratio_scale = exp(exparg_shift)
C      in order to get more dynamic range.
C      qratio_us is the unscaled quantity, onepq = 1 + qratio_us,
C      and onepqs = onepq*qratio_scale
      double precision exparg_shift, qratio_scale, qratio_us, onepq, onepqs
      parameter(exparg_shift = 600.d0)
      double precision ln_ground_occ
      integer index_x(nx)
      data index_x /4, 3, 2, 1, 0/
      logical ifnr03, ifnr13
      integer iffirst
      data iffirst/1/
      save
      if(iffirst.eq.1) then
C        go through this just in case this routine is called without
C        excitation_pi called first.  In normal case this means we
C        have an extra evaluation of qh2, qh2plus, and qstar on first
C        call, but from then on, no extra calls at all.
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
C      sanity checking:
      if(nions.ne.nions_local.or.nions.ne.nions_excitation)
     &  stop 'excitation_sum: inconsistent nions values'
      if(ifnr.eq.2.or.ifnr.gt.3) stop
     &  'excitation_sum: ifnr must be negative, 0, 1, or 3'
      if(ifexcited.le.0) stop 'excitation_sum: invalid ifexcited'
      ifnr03 = ifnr.eq.0.or.ifnr.eq.3
      ifnr13 = ifnr.eq.1.or.ifnr.eq.3
      ifapprox = ifexcited.lt.10
      ifpl_logical = ifpl.eq.1
      ifmhd_logical = ifpi.eq.3.or.ifpi.eq.4
      if(.not.(ifpl_logical.or.ifmhd_logical))
     &  stop 'excitation_sum: invalid combination of ifpl and ifpi'
      c2t = c2*exp(-tl)
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
C      The following if block is identical to code in excitation_pi
C      (except for routine identification on stop statement message).
C      The idea is all the qh2, qh2plus, and qstar values will be taken
C      from previous calculations (either excitation_pi or here)
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
      if( ifdiff_x_excitation.or.
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
     &      stop 'excitation_sum: nmin_max too large'
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
          eps_factor = max(eps_factor_lim, -c2t*bmin(iz))
!          if(eps_factor.gt.eps_factor_lim) then
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
      psum = 0.d0
      ssum = 0.d0
      usum = 0.d0
C      free_sum = ssum - usum proportional to negative of free energy
C      associated with excitation.
      free_sum = 0.d0
      if(ifmhd_logical) then
        xextrasum(1) = 0.d0
        xextrasum(2) = 0.d0
        xextrasum(3) = 0.d0
        xextrasum(4) = 0.d0
      endif
      if(ifnr03) then
        psumf = 0.d0
        psumt = 0.d0
        ssumf = 0.d0
        ssumt = 0.d0
        free_sumf = 0.d0
        if(ifmhd_logical) then
          xextrasumf(1) = 0.d0
          xextrasumf(2) = 0.d0
          xextrasumf(3) = 0.d0
          xextrasumf(4) = 0.d0
          xextrasumt(1) = 0.d0
          xextrasumt(2) = 0.d0
          xextrasumt(3) = 0.d0
          xextrasumt(4) = 0.d0
        endif
      endif
      if(ifnr13) then
        do index = 1, max_index
          psum_dv(index) = 0.d0
          free_sum_dv(index) = 0.d0
        enddo
      endif
      if(ifmhd_logical.and.ifnr13) then
        do index = 1, max_index
          xextrasum_dv(index,1) = 0.d0
          xextrasum_dv(index,2) = 0.d0
          xextrasum_dv(index,3) = 0.d0
          xextrasum_dv(index,4) = 0.d0
        enddo
      endif
C      this do loop only over atoms only.
      do index = 1, n_partial_elements
        ielement = partial_elements(index)
        if((ielement.gt.1.or.ifh2.eq.0).and.
     &      (ielement.le.2.or.mod(ifexcited,10).eq.3)) then
          if(ielement.gt.1) then
            ion_start = ion_end(ielement-1) + 1
          else
            ion_start = 1
          endif
C          calculate delta ln Z = ln(1 + qratio), where
C          qratio is the ratio of excited to ground
C          electronic state partition function
C          n.b. index goes from neutral to next to bare ion
          do ion = ion_start, ion_end(ielement)
            iz = nion(ion)
            nmin_s = nmin_species(ion)
            if(nmin_s.lt.nmin(iz).or.
     &        nmin_s.gt.nmin_max(iz))
     &        stop 'excitation_sum: invalid nmin or nmin_max'
            exparg = exparg_shift - c2t*bion(ion) - plop(ion)
C            correct for MHD ground state occupation probability.
            if(ifmhd_logical) then
              if(iz.eq.1) then
                ln_ground_occ = occ_const*(
     &            x(1) + r_neutral(ielement)*(
     &            x(2)*3.d0 + r_neutral(ielement)*(
     &            x(3)*3.d0 + r_neutral(ielement)*(
     &            x(4)))))
              else
                ln_ground_occ = 0.d0
              endif
              ln_ground_occ = ln_ground_occ +
     &          occ_const*r_ion3(ion)*x(5)
C              ignore excitation correction if ground state wiped out by
C              pressure ionization in any case.
              if(ln_ground_occ.gt.ground_state_lim) then
                exparg = exparg - ln_ground_occ
              else
C                signal to ignore this species
                exparg = -1000.d0!exparg_lim - 1.d0
              endif
            endif
!            if(exparg.gt.exparg_lim) then
              exparg = exp(exparg)
!            else
!              exparg = 0.d0
!            endif
C            qratio is first ratio of excited to ground partition function,
C            but then is transformed to ln(1 + qratio).  Similarly, there
C            is a subsequent transformation of all derivatives.
C            N.B. important convention on partial derivative variable names:
C            names starting with "qratio" are partial derivatives assuming
C            that qratio is a function of tl and x.
C            names starting with "dqratio" are the *change* to the partial
C            derivative caused by x being a function of fl, tl, dv.
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
                eps_factor = max(eps_factor_lim, -c2t*bion(ion))
!                if(eps_factor.gt.eps_factor_lim) then
                  eps_factor = exp(eps_factor)
                  rfactor = 1.d0/(1.d0 + electron_mass/(4.d0*h_mass))
                  call qmhd_calc(ifpi_fit, dble(iqion(ion)),
     &              ifpl_logical, ifapprox, eps_factor,
     &              nmin_rydberg_he1, nmax,
     &              tl, iz, rfactor, nlevels_helium1,
     &              weight_helium1, neff_helium1, ell_helium1, nx, x,
     &              qmhd_he1, qmhd_he1t, qmhd_he1x,
     &              qmhd_he1t2, qmhd_he1tx, qmhd_he1x2)
C                  special helium case.
!C                  zero results (including derivatives below) to provide
!C                  consistent small value zeroing rather than inconsistent
!C                  underflow zeroing.
!                  if(qmhd_he1.lt.exp_lim) qmhd_he1 = 0.d0
                  qratio = exparg*qmhd_he1/
     &              dble(iqneutral(ielement))
!                else
!                  qratio = 0.d0
!                endif
              else
C                neutral atomic case
                qratio = exparg*qstar(nmin_s,iz)*
     &            dble(iqion(ion))/dble(iqneutral(ielement))
              endif
            else
C              ionized atomic case
              qratio = exparg*qstar(nmin_s,iz)*
     &          dble(iqion(ion))/dble(iqion(ion-1))
            endif
            if(qratio.gt.exp_lim) then
              qratio_us = qratio/qratio_scale
              onepq = 1.d0 + qratio_us
              onepqs = onepq*qratio_scale
C              qratiot and qratiot2 are the first and second partial wrt tl
              if(ifhe1_special.and.ion.eq.2.and.ifmhd_logical) then
C                special for neutral helium
                qratiot = qratio*(c2t*bion(ion) - plopt(ion) +
     &            qmhd_he1t/qmhd_he1)
                qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) +
     &            qmhd_he1t/qmhd_he1) +
     &            qratio*(-c2t*bion(ion) - plopt2(ion) +
     &            qmhd_he1t2/qmhd_he1 -
     &            (qmhd_he1t/qmhd_he1)*
     &            (qmhd_he1t/qmhd_he1))
              else
                qratiot = qratio*(c2t*bion(ion) - plopt(ion) +
     &            qstart(nmin_s,iz)/qstar(nmin_s,iz))
                qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) +
     &            qstart(nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(-c2t*bion(ion) - plopt2(ion) +
     &            qstart2(nmin_s,iz)/qstar(nmin_s,iz) -
     &            (qstart(nmin_s,iz)/qstar(nmin_s,iz))*
     &            (qstart(nmin_s,iz)/qstar(nmin_s,iz)))
              endif
C              transform second order quantities to derivative of
C              ln(1 + qratio)
              qratiot2 = (qratiot2 - qratiot*
     &          (qratiot/onepqs))/onepq
              if(ifmhd_logical) then
C                qratio_dx(k), qratiot_dx(k), qratio_dx2(k,l) are 
C                partials wrt tl, extrasum(k) and extrasum(l) following
C                special convention for k or l = 4.
                if(ifhe1_special.and.ion.eq.2) then
C                  special for neutral helium
                  qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +
     &              qmhd_he1x(5)/qmhd_he1)
                  qratiot_dx(4) = qratiot*
     &              (-occ_const*r_ion3(ion) +
     &              qmhd_he1x(5)/qmhd_he1) +
     &              qratio*(qmhd_he1tx(5) -
     &              qmhd_he1x(5)*qmhd_he1t/
     &              qmhd_he1)/qmhd_he1
                  qratio_dx2(4,4) = qratio_dx(4)*
     &              (-occ_const*r_ion3(ion) +
     &              qmhd_he1x(5)/qmhd_he1) +
     &              qratio*(qmhd_he1x2(5,5) - qmhd_he1x(5)*
     &              qmhd_he1x(5)/qmhd_he1)/
     &              qmhd_he1
                else
                  qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +
     &              qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))
                  qratiot_dx(4) = qratiot*
     &              (-occ_const*r_ion3(ion) +
     &              qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +
     &              qratio*(qstartx(5,nmin_s,iz) -
     &              qstarx(5,nmin_s,iz)*qstart(nmin_s,iz)/
     &              qstar(nmin_s,iz))/qstar(nmin_s,iz)
                  qratio_dx2(4,4) = qratio_dx(4)*
     &              (-occ_const*r_ion3(ion) +
     &              qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +
     &              qratio*(qstarx2(5,5,nmin_s,iz) -
     &              qstarx(5,nmin_s,iz)*
     &              qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))/
     &              qstar(nmin_s,iz)
                endif
C                partial wrt ln V, tl, and extrasum(k) following special
C                k convention.
                qratiov = -qratio_dx(4)*extrasum(nextrasum-1)
                qratiovt = -qratiot_dx(4)*extrasum(nextrasum-1)
                qratiov_dx(4) = -qratio_dx2(4,4)*
     &            extrasum(nextrasum-1) -
     &            qratio_dx(4)
C                transform second order quantities to derivative of
C                ln(1 + qratio)
                qratiot_dx(4) = (qratiot_dx(4) - qratiot*
     &            (qratio_dx(4)/onepqs))/onepq
                qratio_dx2(4,4) = (qratio_dx2(4,4) - qratio_dx(4)*
     &            (qratio_dx(4)/onepqs))/onepq
                if(iz.eq.1) then
C                  N.B. x(1) (or extrasum(4)) dependence divides out of qratio
C                  n.b. indices are reordered here so that
C                  qratio_dx(k) refers to derivative wrt extrasum(k) for
C                  k = 1, 2, 3, and qratio_dx(4) refers to the derivative
C                  wrt extrasum(nextrasum-1).
                  if(ifhe1_special.and.ion.eq.2) then
C                    special for neutral helium
                    qratio_dx(1) = qratio*(
     &                -occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qmhd_he1x(4)/qmhd_he1)
                    qratio_dx(2) = qratio*(
     &                -occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qmhd_he1x(3)/qmhd_he1)
                    qratio_dx(3) = qratio*(
     &                -occ_const*r_neutral(ielement)*
     &                3.d0 +
     &                qmhd_he1x(2)/qmhd_he1)
                    qratiot_dx(1) = qratiot*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qmhd_he1x(4)/qmhd_he1) +
     &                qratio*(qmhd_he1tx(4) -
     &                qmhd_he1x(4)*qmhd_he1t/
     &                qmhd_he1)/qmhd_he1
                    qratiot_dx(2) = qratiot*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qmhd_he1x(3)/qmhd_he1) +
     &                qratio*(qmhd_he1tx(3) -
     &                qmhd_he1x(3)*qmhd_he1t/
     &                qmhd_he1)/qmhd_he1
                    qratiot_dx(3) = qratiot*
     &                (-occ_const*r_neutral(ielement)*
     &                3.d0 +
     &                qmhd_he1x(2)/qmhd_he1) +
     &                qratio*(qmhd_he1tx(2) -
     &                qmhd_he1x(2)*qmhd_he1t/
     &                qmhd_he1)/qmhd_he1
                    qratio_dx2(1,1) = qratio_dx(1)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qmhd_he1x(4)/qmhd_he1) +
     &                qratio*(qmhd_he1x2(4,4) -
     &                qmhd_he1x(4)*
     &                qmhd_he1x(4)/qmhd_he1)/
     &                qmhd_he1
                    qratio_dx2(2,1) = qratio_dx(2)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qmhd_he1x(4)/qmhd_he1) +
     &                qratio*(qmhd_he1x2(4,3) -
     &                qmhd_he1x(4)*
     &                qmhd_he1x(3)/qmhd_he1)/
     &                qmhd_he1
                    qratio_dx2(3,1) = qratio_dx(3)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qmhd_he1x(4)/qmhd_he1) +
     &                qratio*(qmhd_he1x2(4,2) -
     &                qmhd_he1x(4)*
     &                qmhd_he1x(2)/qmhd_he1)/
     &                qmhd_he1
                    qratio_dx2(4,1) = qratio_dx(4)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qmhd_he1x(4)/qmhd_he1) +
     &                qratio*(qmhd_he1x2(5,4) -
     &                qmhd_he1x(5)*
     &                qmhd_he1x(4)/qmhd_he1)/
     &                qmhd_he1
                    qratio_dx2(2,2) = qratio_dx(2)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qmhd_he1x(3)/qmhd_he1) +
     &                qratio*(qmhd_he1x2(3,3) -
     &                qmhd_he1x(3)*
     &                qmhd_he1x(3)/qmhd_he1)/
     &                qmhd_he1
                    qratio_dx2(3,2) = qratio_dx(3)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qmhd_he1x(3)/qmhd_he1) +
     &                qratio*(qmhd_he1x2(3,2) -
     &                qmhd_he1x(3)*
     &                qmhd_he1x(2)/qmhd_he1)/
     &                qmhd_he1
                    qratio_dx2(4,2) = qratio_dx(4)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qmhd_he1x(3)/qmhd_he1) +
     &                qratio*(qmhd_he1x2(5,3) -
     &                qmhd_he1x(5)*
     &                qmhd_he1x(3)/qmhd_he1)/
     &                qmhd_he1
                    qratio_dx2(3,3) = qratio_dx(3)*
     &                (-occ_const*r_neutral(ielement)*
     &                3.d0 +
     &                qmhd_he1x(2)/qmhd_he1) +
     &                qratio*(qmhd_he1x2(2,2) -
     &                qmhd_he1x(2)*
     &                qmhd_he1x(2)/qmhd_he1)/
     &                qmhd_he1
                    qratio_dx2(4,3) = qratio_dx(4)*
     &                (-occ_const*r_neutral(ielement)*
     &                3.d0 +
     &                qmhd_he1x(2)/qmhd_he1) +
     &                qratio*(qmhd_he1x2(5,2) -
     &                qmhd_he1x(5)*
     &                qmhd_he1x(2)/qmhd_he1)/
     &                qmhd_he1
                  else
                    qratio_dx(1) = qratio*(
     &                -occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))
                    qratio_dx(2) = qratio*(
     &                -occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))
                    qratio_dx(3) = qratio*(
     &                -occ_const*r_neutral(ielement)*
     &                3.d0 +
     &                qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))
                    qratiot_dx(1) = qratiot*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstartx(4,nmin_s,iz) -
     &                qstarx(4,nmin_s,iz)*qstart(nmin_s,iz)/
     &                qstar(nmin_s,iz))/qstar(nmin_s,iz)
                    qratiot_dx(2) = qratiot*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstartx(3,nmin_s,iz) -
     &                qstarx(3,nmin_s,iz)*qstart(nmin_s,iz)/
     &                qstar(nmin_s,iz))/qstar(nmin_s,iz)
                    qratiot_dx(3) = qratiot*
     &                (-occ_const*r_neutral(ielement)*
     &                3.d0 +
     &                qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstartx(2,nmin_s,iz) -
     &                qstarx(2,nmin_s,iz)*qstart(nmin_s,iz)/
     &                qstar(nmin_s,iz))/qstar(nmin_s,iz)
                    qratio_dx2(1,1) = qratio_dx(1)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstarx2(4,4,nmin_s,iz) -
     &                qstarx(4,nmin_s,iz)*
     &                qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))/
     &                qstar(nmin_s,iz)
                    qratio_dx2(2,1) = qratio_dx(2)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstarx2(4,3,nmin_s,iz) -
     &                qstarx(4,nmin_s,iz)*
     &                qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))/
     &                qstar(nmin_s,iz)
                    qratio_dx2(3,1) = qratio_dx(3)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstarx2(4,2,nmin_s,iz) -
     &                qstarx(4,nmin_s,iz)*
     &                qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &                qstar(nmin_s,iz)
                    qratio_dx2(4,1) = qratio_dx(4)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*r_neutral(ielement) +
     &                qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstarx2(5,4,nmin_s,iz) -
     &                qstarx(5,nmin_s,iz)*
     &                qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))/
     &                qstar(nmin_s,iz)
                    qratio_dx2(2,2) = qratio_dx(2)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstarx2(3,3,nmin_s,iz) -
     &                qstarx(3,nmin_s,iz)*
     &                qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))/
     &                qstar(nmin_s,iz)
                    qratio_dx2(3,2) = qratio_dx(3)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstarx2(3,2,nmin_s,iz) -
     &                qstarx(3,nmin_s,iz)*
     &                qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &                qstar(nmin_s,iz)
                    qratio_dx2(4,2) = qratio_dx(4)*
     &                (-occ_const*r_neutral(ielement)*
     &                r_neutral(ielement)*3.d0 +
     &                qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstarx2(5,3,nmin_s,iz) -
     &                qstarx(5,nmin_s,iz)*
     &                qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))/
     &                qstar(nmin_s,iz)
                    qratio_dx2(3,3) = qratio_dx(3)*
     &                (-occ_const*r_neutral(ielement)*
     &                3.d0 +
     &                qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstarx2(2,2,nmin_s,iz) -
     &                qstarx(2,nmin_s,iz)*
     &                qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &                qstar(nmin_s,iz)
                    qratio_dx2(4,3) = qratio_dx(4)*
     &                (-occ_const*r_neutral(ielement)*
     &                3.d0 +
     &                qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +
     &                qratio*(qstarx2(5,2,nmin_s,iz) -
     &                qstarx(5,nmin_s,iz)*
     &                qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &                qstar(nmin_s,iz)
                  endif !if(ifhe1_special.and.ion.eq.2) ...
C                  derivative of qratio wrt ln V
                  qratiov = qratiov -
     &              qratio_dx(1)*extrasum(1) -
     &              qratio_dx(2)*extrasum(2) -
     &              qratio_dx(3)*extrasum(3)
                  qratiovt = qratiovt -
     &              qratiot_dx(1)*extrasum(1) -
     &              qratiot_dx(2)*extrasum(2) -
     &              qratiot_dx(3)*extrasum(3)
                  qratiov_dx(1) = -
     &              qratio_dx2(1,1)*extrasum(1) -
     &              qratio_dx2(2,1)*extrasum(2) -
     &              qratio_dx2(3,1)*extrasum(3) -
     &              qratio_dx2(4,1)*extrasum(nextrasum-1) -
     &              qratio_dx(1)
                  qratiov_dx(2) = -
     &              qratio_dx2(2,1)*extrasum(1) -
     &              qratio_dx2(2,2)*extrasum(2) -
     &              qratio_dx2(3,2)*extrasum(3) -
     &              qratio_dx2(4,2)*extrasum(nextrasum-1) -
     &              qratio_dx(2)
                  qratiov_dx(3) = -
     &              qratio_dx2(3,1)*extrasum(1) -
     &              qratio_dx2(3,2)*extrasum(2) -
     &              qratio_dx2(3,3)*extrasum(3) -
     &              qratio_dx2(4,3)*extrasum(nextrasum-1) -
     &              qratio_dx(3)
                  qratiov_dx(4) = qratiov_dx(4) -
     &              qratio_dx2(4,1)*extrasum(1) -
     &              qratio_dx2(4,2)*extrasum(2) -
     &              qratio_dx2(4,3)*extrasum(3)
C                  transform second order quantities to derivative of
C                  ln(1 + qratio)
                  qratiot_dx(1) = (qratiot_dx(1) - qratiot*
     &              (qratio_dx(1)/onepqs))/onepq
                  qratiot_dx(2) = (qratiot_dx(2) - qratiot*
     &              (qratio_dx(2)/onepqs))/onepq
                  qratiot_dx(3) = (qratiot_dx(3) - qratiot*
     &              (qratio_dx(3)/onepqs))/onepq
                  qratio_dx2(1,1) = (qratio_dx2(1,1) - qratio_dx(1)*
     &              (qratio_dx(1)/onepqs))/onepq
                  qratio_dx2(2,1) = (qratio_dx2(2,1) - qratio_dx(2)*
     &              (qratio_dx(1)/onepqs))/onepq
                  qratio_dx2(3,1) = (qratio_dx2(3,1) - qratio_dx(3)*
     &              (qratio_dx(1)/onepqs))/onepq
                  qratio_dx2(4,1) = (qratio_dx2(4,1) - qratio_dx(4)*
     &              (qratio_dx(1)/onepqs))/onepq
                  qratio_dx2(2,2) = (qratio_dx2(2,2) - qratio_dx(2)*
     &              (qratio_dx(2)/onepqs))/onepq
                  qratio_dx2(3,2) = (qratio_dx2(3,2) - qratio_dx(3)*
     &              (qratio_dx(2)/onepqs))/onepq
                  qratio_dx2(4,2) = (qratio_dx2(4,2) - qratio_dx(4)*
     &              (qratio_dx(2)/onepqs))/onepq
                  qratio_dx2(3,3) = (qratio_dx2(3,3) - qratio_dx(3)*
     &              (qratio_dx(3)/onepqs))/onepq
                  qratio_dx2(4,3) = (qratio_dx2(4,3) - qratio_dx(4)*
     &              (qratio_dx(3)/onepqs))/onepq
                  qratiov_dx(1) = (qratiov_dx(1) - qratiov*
     &              (qratio_dx(1)/onepqs))/onepq
                  qratiov_dx(2) = (qratiov_dx(2) - qratiov*
     &              (qratio_dx(2)/onepqs))/onepq
                  qratiov_dx(3) = (qratiov_dx(3) - qratiov*
     &              (qratio_dx(3)/onepqs))/onepq
C                  transform first order quantities to derivative of
C                  ln(1 + qratio)
                  qratio_dx(1) = qratio_dx(1)/onepq
                  qratio_dx(2) = qratio_dx(2)/onepq
                  qratio_dx(3) = qratio_dx(3)/onepq
                endif
C                transform second order quantities to derivative of
C                ln(1 + qratio)
                qratiovt = (qratiovt - qratiov*
     &            (qratiot/onepqs))/onepq
                qratiov_dx(4) = (qratiov_dx(4) - qratiov*
     &            (qratio_dx(4)/onepqs))/onepq
C                transform first order quantities to derivative of
C                ln(1 + qratio)
                qratio_dx(4) = qratio_dx(4)/onepq
                qratiov = qratiov/onepq
              endif
C              transform first order quantities to derivative of
C              ln(1 + qratio)
              qratiot = qratiot/onepq
              if(qratio_us.gt.1.d-3) then
                qratio = qratio_scale*log(onepq)
              else
C                alternating series so relative error is less than first
C                missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
                qratio = qratio*
     &            (1.d0      - qratio_us*
     &            (1.d0/2.d0 - qratio_us*
     &            (1.d0/3.d0 - qratio_us*
     &            (1.d0/4.d0 - qratio_us*
     &            (1.d0/5.d0)))))
              endif
              if(ifmhd_logical) then
C                calculate derivative *correction* due to x dependence on fl
C                and tl.
                if(ifnr03) then
                  dqratiof = qratio_dx(4)*extrasumf(nextrasum-1)
                  dqratiot = qratio_dx(4)*extrasumt(nextrasum-1)
                  dqratiotf = qratiot_dx(4)*extrasumf(nextrasum-1)
                  dqratiott = qratiot_dx(4)*extrasumt(nextrasum-1)
                  dqratio_dxf(4) = qratio_dx2(4,4)*
     &              extrasumf(nextrasum-1)
                  dqratio_dxt(4) = qratio_dx2(4,4)*
     &              extrasumt(nextrasum-1)
                  dqratiovf = qratiov_dx(4)*extrasumf(nextrasum-1)
                  dqratiovt = qratiov_dx(4)*extrasumt(nextrasum-1)
                endif
                if(ifnr13) then
                  do inv_ion_index = 1, max_index
                    dqratio_dv(inv_ion_index) = qratio_dx(4)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratiov_dv(inv_ion_index) = qratiov_dx(4)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,4) = qratio_dx2(4,4)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                  enddo
                endif
                if(iz.eq.1) then
                  if(ifnr03) then
                    dqratiof = dqratiof +
     &                qratio_dx(1)*extrasumf(1) +
     &                qratio_dx(2)*extrasumf(2) +
     &                qratio_dx(3)*extrasumf(3)
                    dqratiot = dqratiot +
     &                qratio_dx(1)*extrasumt(1) +
     &                qratio_dx(2)*extrasumt(2) +
     &                qratio_dx(3)*extrasumt(3)
                    dqratiotf = dqratiotf +
     &                qratiot_dx(1)*extrasumf(1) +
     &                qratiot_dx(2)*extrasumf(2) +
     &                qratiot_dx(3)*extrasumf(3)
                    dqratiott = dqratiott +
     &                qratiot_dx(1)*extrasumt(1) +
     &                qratiot_dx(2)*extrasumt(2) +
     &                qratiot_dx(3)*extrasumt(3)
                    dqratio_dxf(1) =
     &                qratio_dx2(1,1)*extrasumf(1) +
     &                qratio_dx2(2,1)*extrasumf(2) +
     &                qratio_dx2(3,1)*extrasumf(3) +
     &                qratio_dx2(4,1)*extrasumf(nextrasum-1)
                    dqratio_dxt(1) =
     &                qratio_dx2(1,1)*extrasumt(1) +
     &                qratio_dx2(2,1)*extrasumt(2) +
     &                qratio_dx2(3,1)*extrasumt(3) +
     &                qratio_dx2(4,1)*extrasumt(nextrasum-1)
                    dqratio_dxf(2) =
     &                qratio_dx2(2,1)*extrasumf(1) +
     &                qratio_dx2(2,2)*extrasumf(2) +
     &                qratio_dx2(3,2)*extrasumf(3) +
     &                qratio_dx2(4,2)*extrasumf(nextrasum-1)
                    dqratio_dxt(2) =
     &                qratio_dx2(2,1)*extrasumt(1) +
     &                qratio_dx2(2,2)*extrasumt(2) +
     &                qratio_dx2(3,2)*extrasumt(3) +
     &                qratio_dx2(4,2)*extrasumt(nextrasum-1)
                    dqratio_dxf(3) =
     &                qratio_dx2(3,1)*extrasumf(1) +
     &                qratio_dx2(3,2)*extrasumf(2) +
     &                qratio_dx2(3,3)*extrasumf(3) +
     &                qratio_dx2(4,3)*extrasumf(nextrasum-1)
                    dqratio_dxt(3) =
     &                qratio_dx2(3,1)*extrasumt(1) +
     &                qratio_dx2(3,2)*extrasumt(2) +
     &                qratio_dx2(3,3)*extrasumt(3) +
     &                qratio_dx2(4,3)*extrasumt(nextrasum-1)
                    dqratio_dxf(4) = dqratio_dxf(4) +
     &                qratio_dx2(4,1)*extrasumf(1) +
     &                qratio_dx2(4,2)*extrasumf(2) +
     &                qratio_dx2(4,3)*extrasumf(3)
                    dqratio_dxt(4) = dqratio_dxt(4) +
     &                qratio_dx2(4,1)*extrasumt(1) +
     &                qratio_dx2(4,2)*extrasumt(2) +
     &                qratio_dx2(4,3)*extrasumt(3)
                    dqratiovf = dqratiovf +
     &                qratiov_dx(1)*extrasumf(1) +
     &                qratiov_dx(2)*extrasumf(2) +
     &                qratiov_dx(3)*extrasumf(3)
                    dqratiovt = dqratiovt +
     &                qratiov_dx(1)*extrasumt(1) +
     &                qratiov_dx(2)*extrasumt(2) +
     &                qratiov_dx(3)*extrasumt(3)
                  endif
                  if(ifnr13) then
                    do inv_ion_index = 1, max_index
                      dqratio_dv(inv_ion_index) =
     &                  dqratio_dv(inv_ion_index) +
     &                  qratio_dx(1)*extrasum_dv(inv_ion_index,1) +
     &                  qratio_dx(2)*extrasum_dv(inv_ion_index,2) +
     &                  qratio_dx(3)*extrasum_dv(inv_ion_index,3)
                      dqratiov_dv(inv_ion_index) =
     &                  dqratiov_dv(inv_ion_index) +
     &                  qratiov_dx(1)*extrasum_dv(inv_ion_index,1) +
     &                  qratiov_dx(2)*extrasum_dv(inv_ion_index,2) +
     &                  qratiov_dx(3)*extrasum_dv(inv_ion_index,3)
                      dqratio_dxdv(inv_ion_index,1) =
     &                  qratio_dx2(1,1)*
     &                  extrasum_dv(inv_ion_index,1) +
     &                  qratio_dx2(2,1)*
     &                  extrasum_dv(inv_ion_index,2) +
     &                  qratio_dx2(3,1)*
     &                  extrasum_dv(inv_ion_index,3) +
     &                  qratio_dx2(4,1)*
     &                  extrasum_dv(inv_ion_index,nextrasum-1)
                      dqratio_dxdv(inv_ion_index,2) =
     &                  qratio_dx2(2,1)*
     &                  extrasum_dv(inv_ion_index,1) +
     &                  qratio_dx2(2,2)*
     &                  extrasum_dv(inv_ion_index,2) +
     &                  qratio_dx2(3,2)*
     &                  extrasum_dv(inv_ion_index,3) +
     &                  qratio_dx2(4,2)*
     &                  extrasum_dv(inv_ion_index,nextrasum-1)
                      dqratio_dxdv(inv_ion_index,3) =
     &                  qratio_dx2(3,1)*
     &                  extrasum_dv(inv_ion_index,1) +
     &                  qratio_dx2(3,2)*
     &                  extrasum_dv(inv_ion_index,2) +
     &                  qratio_dx2(3,3)*
     &                  extrasum_dv(inv_ion_index,3) +
     &                  qratio_dx2(4,3)*
     &                  extrasum_dv(inv_ion_index,nextrasum-1)
                      dqratio_dxdv(inv_ion_index,4) =
     &                  dqratio_dxdv(inv_ion_index,4) +
     &                  qratio_dx2(4,1)*
     &                  extrasum_dv(inv_ion_index,1) +
     &                  qratio_dx2(4,2)*
     &                  extrasum_dv(inv_ion_index,2) +
     &                  qratio_dx2(4,3)*
     &                  extrasum_dv(inv_ion_index,3)
                    enddo
                  endif ! if(ifnr03) then.
                endif ! if(iz.eq.1) then ...
              endif ! if(ifmhd_logical) then.....
              call exsum_component_add(
     &          ifmhd_logical, iz.eq.1, ifnr03, ifnr13,.false.,
     &          nxextrasum, nions, maxionstage_nuvar,
     &          ion_start, ion_end(ielement),
     &          1, max_index,
     &          inv_ion,
     &          nuvar(iz,index), nuvarf(iz,index),
     &          nuvart(iz,index), nuvar_dv(1, iz,index),
     &          xextrasum_scale, qratio_scale,
     &          qratio, qratiot, qratio_dx,
     &          qratiot2, qratiot_dx,
     &          qratio_dx2,
     &          qratiov, qratiovt, qratiov_dx,
     &          dqratiof, dqratiotf, dqratio_dxf,
     &          dqratiot, dqratiott, dqratio_dxt,
     &          dqratio_dv,
     &          dqratio_dxdv,
     &          dqratiovf, dqratiovt, dqratiov_dv,
     &          psum, psumf, psumt, psum_dv,
     &          ssum, ssumf, ssumt,
     &          usum,
     &          free_sum, free_sumf, free_sum_dv,
     &          xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
            endif ! if(qratio.gt.exp_lim) then
          enddo ! do ion = ion_start, i... neutral to next to bare.
        endif ! if((ielement.gt
      enddo ! do index = 1,... for atoms only
      if(partial_elements(1).eq.1) then
C         do all the molecular stuff below only if hydrogen has a non-zero
C         abundance
C         now do H2 equilibrium constant change due to neutral monatomic H.
C         now do neutral monatomic H
C        (which skipped previously if molecular formation).
        if(ifh2.gt.0) then
          ielement = 1
          iz = 1
          ion = 1
          nmin_s = nmin_species(ion)
          if(nmin_s.lt.nmin(iz).or.
     &      nmin_s.gt.nmin_max(iz))
     &      stop 'excitation_sum: invalid nmin or nmin_max'
          exparg = exparg_shift -
     &      c2t*bion(ion) - plop(ion)
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
C              signal to ignore this species
              exparg = -1000.d0!exparg_lim - 1.d0
            endif
          endif
!          if(exparg.gt.exparg_lim) then
            exparg = exp(exparg)
!          else
!            exparg = 0.d0
!          endif
C          qratio is first ratio of excited to ground partition function,
C          but then is transformed to ln(1 + qratio).  Similarly, there
C          is a subsequent transformation of all derivatives.
C          N.B. important convention on partial derivative variable names:
C          names starting with "qratio" are partial derivatives assuming
C          that qratio is a function of tl and x.
C          names starting with "dqratio" are the *change* to the partial
C          derivative caused by x being a function of fl, tl, dv.
          if(iz.eq.1) then
C            neutral monatomic H as part of H2 calculation
            qratio = exparg*qstar(nmin_s,iz)*
     &        dble(iqion(ion))/dble(iqneutral(ielement))
          else
C            unused since iz set to 1 above.
            qratio = exparg*qstar(nmin_s,iz)*
     &        dble(iqion(ion))/dble(iqion(ion-1))
          endif
          if(qratio.gt.exp_lim) then
            qratio_us = qratio/qratio_scale
            onepq = 1.d0 + qratio_us
            onepqs = onepq*qratio_scale
C            qratiot and qratiot2 are the first and second partial wrt tl
            qratiot = qratio*(c2t*bion(ion) - plopt(ion) +
     &        qstart(nmin_s,iz)/qstar(nmin_s,iz))
            qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) +
     &        qstart(nmin_s,iz)/qstar(nmin_s,iz)) +
     &        qratio*(-c2t*bion(ion) - plopt2(ion) +
     &        qstart2(nmin_s,iz)/qstar(nmin_s,iz) -
     &        (qstart(nmin_s,iz)/qstar(nmin_s,iz))*
     &        (qstart(nmin_s,iz)/qstar(nmin_s,iz)))
C            transform second order quantities to derivative of
C            ln(1 + qratio)
            qratiot2 = (qratiot2 - qratiot*
     &        (qratiot/onepqs))/onepq
            if(ifmhd_logical) then
C              qratio_dx(k), qratiot_dx(k), qratio_dx2(k,l) are 
C              partials wrt tl, extrasum(k) and extrasum(l) following
C              special convention for k or l = 4.
              qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +
     &          qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))
              qratiot_dx(4) = qratiot*
     &          (-occ_const*r_ion3(ion) +
     &          qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +
     &          qratio*(qstartx(5,nmin_s,iz) -
     &          qstarx(5,nmin_s,iz)*qstart(nmin_s,iz)/
     &          qstar(nmin_s,iz))/qstar(nmin_s,iz)
              qratio_dx2(4,4) = qratio_dx(4)*
     &          (-occ_const*r_ion3(ion) +
     &          qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +
     &          qratio*(qstarx2(5,5,nmin_s,iz) -
     &          qstarx(5,nmin_s,iz)*
     &          qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))/
     &          qstar(nmin_s,iz)
C              derivative of qratio wrt ln V
              qratiov = -qratio_dx(4)*extrasum(nextrasum-1)
              qratiovt = -qratiot_dx(4)*extrasum(nextrasum-1)
              qratiov_dx(4) =
     &          -qratio_dx2(4,4)*extrasum(nextrasum-1) -
     &          qratio_dx(4)
C              transform second order quantities to derivative of
C              ln(1 + qratio)
              qratiot_dx(4) = (qratiot_dx(4) - qratiot*
     &          (qratio_dx(4)/onepqs))/onepq
              qratio_dx2(4,4) = (qratio_dx2(4,4) - qratio_dx(4)*
     &          (qratio_dx(4)/onepqs))/onepq
              if(iz.eq.1) then
C                N.B. x(1) (or extrasum(4)) dependence divides out of qratio
C                n.b. indices are reordered here so that
C                qratio_dx(k) refers to derivative wrt extrasum(k) for
C                k = 1, 2, 3, and qratio_dx(4) refers to the derivative
C                wrt extrasum(nextrasum-1).
                qratio_dx(1) = qratio*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))
                qratio_dx(2) = qratio*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))
                qratio_dx(3) = qratio*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))
                qratiot_dx(1) = qratiot*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstartx(4,nmin_s,iz) -
     &            qstarx(4,nmin_s,iz)*qstart(nmin_s,iz)/
     &            qstar(nmin_s,iz))/qstar(nmin_s,iz)
                qratiot_dx(2) = qratiot*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstartx(3,nmin_s,iz) -
     &            qstarx(3,nmin_s,iz)*qstart(nmin_s,iz)/
     &            qstar(nmin_s,iz))/qstar(nmin_s,iz)
                qratiot_dx(3) = qratiot*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstartx(2,nmin_s,iz) -
     &            qstarx(2,nmin_s,iz)*qstart(nmin_s,iz)/
     &            qstar(nmin_s,iz))/qstar(nmin_s,iz)
                qratio_dx2(1,1) = qratio_dx(1)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(4,4,nmin_s,iz) -
     &            qstarx(4,nmin_s,iz)*
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(2,1) = qratio_dx(2)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(4,3,nmin_s,iz) -
     &            qstarx(4,nmin_s,iz)*
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(3,1) = qratio_dx(3)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(4,2,nmin_s,iz) -
     &            qstarx(4,nmin_s,iz)*
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(4,1) = qratio_dx(4)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(5,4,nmin_s,iz) -
     &            qstarx(5,nmin_s,iz)*
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(2,2) = qratio_dx(2)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(3,3,nmin_s,iz) -
     &            qstarx(3,nmin_s,iz)*
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(3,2) = qratio_dx(3)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(3,2,nmin_s,iz) -
     &            qstarx(3,nmin_s,iz)*
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(4,2) = qratio_dx(4)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(5,3,nmin_s,iz) -
     &            qstarx(5,nmin_s,iz)*
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(3,3) = qratio_dx(3)*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(2,2,nmin_s,iz) -
     &            qstarx(2,nmin_s,iz)*
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(4,3) = qratio_dx(4)*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(5,2,nmin_s,iz) -
     &            qstarx(5,nmin_s,iz)*
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
C                derivative of qratio wrt ln V
                qratiov = qratiov -
     &            qratio_dx(1)*extrasum(1) -
     &            qratio_dx(2)*extrasum(2) -
     &            qratio_dx(3)*extrasum(3)
                qratiovt = qratiovt -
     &            qratiot_dx(1)*extrasum(1) -
     &            qratiot_dx(2)*extrasum(2) -
     &            qratiot_dx(3)*extrasum(3)
                qratiov_dx(1) = -
     &            qratio_dx2(1,1)*extrasum(1) -
     &            qratio_dx2(2,1)*extrasum(2) -
     &            qratio_dx2(3,1)*extrasum(3) -
     &            qratio_dx2(4,1)*extrasum(nextrasum-1) -
     &            qratio_dx(1)
                qratiov_dx(2) = -
     &            qratio_dx2(2,1)*extrasum(1) -
     &            qratio_dx2(2,2)*extrasum(2) -
     &            qratio_dx2(3,2)*extrasum(3) -
     &            qratio_dx2(4,2)*extrasum(nextrasum-1) -
     &            qratio_dx(2)
                qratiov_dx(3) = -
     &            qratio_dx2(3,1)*extrasum(1) -
     &            qratio_dx2(3,2)*extrasum(2) -
     &            qratio_dx2(3,3)*extrasum(3) -
     &            qratio_dx2(4,3)*extrasum(nextrasum-1) -
     &            qratio_dx(3)
                qratiov_dx(4) = qratiov_dx(4) -
     &            qratio_dx2(4,1)*extrasum(1) -
     &            qratio_dx2(4,2)*extrasum(2) -
     &            qratio_dx2(4,3)*extrasum(3)
C                transform second order quantities to derivative of
C                ln(1 + qratio)
                qratiot_dx(1) = (qratiot_dx(1) - qratiot*
     &            (qratio_dx(1)/onepqs))/onepq
                qratiot_dx(2) = (qratiot_dx(2) - qratiot*
     &            (qratio_dx(2)/onepqs))/onepq
                qratiot_dx(3) = (qratiot_dx(3) - qratiot*
     &            (qratio_dx(3)/onepqs))/onepq
                qratio_dx2(1,1) = (qratio_dx2(1,1) - qratio_dx(1)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(2,1) = (qratio_dx2(2,1) - qratio_dx(2)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(3,1) = (qratio_dx2(3,1) - qratio_dx(3)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(4,1) = (qratio_dx2(4,1) - qratio_dx(4)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(2,2) = (qratio_dx2(2,2) - qratio_dx(2)*
     &            (qratio_dx(2)/onepqs))/onepq
                qratio_dx2(3,2) = (qratio_dx2(3,2) - qratio_dx(3)*
     &            (qratio_dx(2)/onepqs))/onepq
                qratio_dx2(4,2) = (qratio_dx2(4,2) - qratio_dx(4)*
     &            (qratio_dx(2)/onepqs))/onepq
                qratio_dx2(3,3) = (qratio_dx2(3,3) - qratio_dx(3)*
     &            (qratio_dx(3)/onepqs))/onepq
                qratio_dx2(4,3) = (qratio_dx2(4,3) - qratio_dx(4)*
     &            (qratio_dx(3)/onepqs))/onepq
                qratiov_dx(1) = (qratiov_dx(1) - qratiov*
     &            (qratio_dx(1)/onepqs))/onepq
                qratiov_dx(2) = (qratiov_dx(2) - qratiov*
     &            (qratio_dx(2)/onepqs))/onepq
                qratiov_dx(3) = (qratiov_dx(3) - qratiov*
     &            (qratio_dx(3)/onepqs))/onepq
C                transform first order quantities to derivative of
C                ln(1 + qratio)
                qratio_dx(1) = qratio_dx(1)/onepq
                qratio_dx(2) = qratio_dx(2)/onepq
                qratio_dx(3) = qratio_dx(3)/onepq
              endif ! if(iz.eq.1) then....
C              transform second order quantities to derivative of
C              ln(1 + qratio)
              qratiovt = (qratiovt - qratiov*
     &          (qratiot/onepqs))/onepq
              qratiov_dx(4) = (qratiov_dx(4) - qratiov*
     &          (qratio_dx(4)/onepqs))/onepq
C              transform first order quantities to derivative of
C              ln(1 + qratio)
              qratio_dx(4) = qratio_dx(4)/onepq
              qratiov = qratiov/onepq
            endif ! if(ifmhd_logical) then....
C            transform first order quantities to derivative of
C            ln(1 + qratio)
            qratiot = qratiot/onepq
            if(qratio_us.gt.1.d-3) then
              qratio = qratio_scale*log(onepq)
            else
C              alternating series so relative error is less than first
C              missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
              qratio = qratio*
     &          (1.d0      - qratio_us*
     &          (1.d0/2.d0 - qratio_us*
     &          (1.d0/3.d0 - qratio_us*
     &          (1.d0/4.d0 - qratio_us*
     &          (1.d0/5.d0)))))
            endif
            if(ifmhd_logical) then
C              calculate derivative *correction* due to x dependence on fl
C              and tl.
              if(ifnr03) then
                dqratiof = qratio_dx(4)*extrasumf(nextrasum-1)
                dqratiot = qratio_dx(4)*extrasumt(nextrasum-1)
                dqratiotf = qratiot_dx(4)*extrasumf(nextrasum-1)
                dqratiott = qratiot_dx(4)*extrasumt(nextrasum-1)
                dqratio_dxf(4) = qratio_dx2(4,4)*
     &            extrasumf(nextrasum-1)
                dqratio_dxt(4) = qratio_dx2(4,4)*
     &            extrasumt(nextrasum-1)
                dqratiovf = qratiov_dx(4)*extrasumf(nextrasum-1)
                dqratiovt = qratiov_dx(4)*extrasumt(nextrasum-1)
              endif
              if(ifnr13) then
                do inv_ion_index = 1, max_index
                  dqratio_dv(inv_ion_index) = qratio_dx(4)*
     &              extrasum_dv(inv_ion_index,nextrasum-1)
                  dqratiov_dv(inv_ion_index) = qratiov_dx(4)*
     &              extrasum_dv(inv_ion_index,nextrasum-1)
                  dqratio_dxdv(inv_ion_index,4) = qratio_dx2(4,4)*
     &              extrasum_dv(inv_ion_index,nextrasum-1)
                enddo
              endif
              if(iz.eq.1) then
                if(ifnr03) then
                  dqratiof = dqratiof +
     &              qratio_dx(1)*extrasumf(1) +
     &              qratio_dx(2)*extrasumf(2) +
     &              qratio_dx(3)*extrasumf(3)
                  dqratiot = dqratiot +
     &              qratio_dx(1)*extrasumt(1) +
     &              qratio_dx(2)*extrasumt(2) +
     &              qratio_dx(3)*extrasumt(3)
                  dqratiotf = dqratiotf +
     &              qratiot_dx(1)*extrasumf(1) +
     &              qratiot_dx(2)*extrasumf(2) +
     &              qratiot_dx(3)*extrasumf(3)
                  dqratiott = dqratiott +
     &              qratiot_dx(1)*extrasumt(1) +
     &              qratiot_dx(2)*extrasumt(2) +
     &              qratiot_dx(3)*extrasumt(3)
                  dqratio_dxf(1) =
     &              qratio_dx2(1,1)*extrasumf(1) +
     &              qratio_dx2(2,1)*extrasumf(2) +
     &              qratio_dx2(3,1)*extrasumf(3) +
     &              qratio_dx2(4,1)*extrasumf(nextrasum-1)
                  dqratio_dxt(1) =
     &              qratio_dx2(1,1)*extrasumt(1) +
     &              qratio_dx2(2,1)*extrasumt(2) +
     &              qratio_dx2(3,1)*extrasumt(3) +
     &              qratio_dx2(4,1)*extrasumt(nextrasum-1)
                  dqratio_dxf(2) =
     &              qratio_dx2(2,1)*extrasumf(1) +
     &              qratio_dx2(2,2)*extrasumf(2) +
     &              qratio_dx2(3,2)*extrasumf(3) +
     &              qratio_dx2(4,2)*extrasumf(nextrasum-1)
                  dqratio_dxt(2) =
     &              qratio_dx2(2,1)*extrasumt(1) +
     &              qratio_dx2(2,2)*extrasumt(2) +
     &              qratio_dx2(3,2)*extrasumt(3) +
     &              qratio_dx2(4,2)*extrasumt(nextrasum-1)
                  dqratio_dxf(3) =
     &              qratio_dx2(3,1)*extrasumf(1) +
     &              qratio_dx2(3,2)*extrasumf(2) +
     &              qratio_dx2(3,3)*extrasumf(3) +
     &              qratio_dx2(4,3)*extrasumf(nextrasum-1)
                  dqratio_dxt(3) =
     &              qratio_dx2(3,1)*extrasumt(1) +
     &              qratio_dx2(3,2)*extrasumt(2) +
     &              qratio_dx2(3,3)*extrasumt(3) +
     &              qratio_dx2(4,3)*extrasumt(nextrasum-1)
                  dqratio_dxf(4) = dqratio_dxf(4) +
     &              qratio_dx2(4,1)*extrasumf(1) +
     &              qratio_dx2(4,2)*extrasumf(2) +
     &              qratio_dx2(4,3)*extrasumf(3)
                  dqratio_dxt(4) = dqratio_dxt(4) +
     &              qratio_dx2(4,1)*extrasumt(1) +
     &              qratio_dx2(4,2)*extrasumt(2) +
     &              qratio_dx2(4,3)*extrasumt(3)
                  dqratiovf = dqratiovf +
     &              qratiov_dx(1)*extrasumf(1) +
     &              qratiov_dx(2)*extrasumf(2) +
     &              qratiov_dx(3)*extrasumf(3)
                  dqratiovt = dqratiovt +
     &              qratiov_dx(1)*extrasumt(1) +
     &              qratiov_dx(2)*extrasumt(2) +
     &              qratiov_dx(3)*extrasumt(3)
                endif
                if(ifnr13) then
                  do inv_ion_index = 1, max_index
                    dqratio_dv(inv_ion_index) =
     &                dqratio_dv(inv_ion_index) +
     &                qratio_dx(1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx(2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx(3)*extrasum_dv(inv_ion_index,3)
                    dqratiov_dv(inv_ion_index) =
     &                dqratiov_dv(inv_ion_index) +
     &                qratiov_dx(1)*extrasum_dv(inv_ion_index,1) +
     &                qratiov_dx(2)*extrasum_dv(inv_ion_index,2) +
     &                qratiov_dx(3)*extrasum_dv(inv_ion_index,3)
                    dqratio_dxdv(inv_ion_index,1) =
     &                qratio_dx2(1,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(2,1)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(3,1)*extrasum_dv(inv_ion_index,3) +
     &                qratio_dx2(4,1)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,2) =
     &                qratio_dx2(2,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(2,2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(3,2)*extrasum_dv(inv_ion_index,3) +
     &                qratio_dx2(4,2)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,3) =
     &                qratio_dx2(3,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(3,2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(3,3)*extrasum_dv(inv_ion_index,3) +
     &                qratio_dx2(4,3)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,4) =
     &                dqratio_dxdv(inv_ion_index,4) +
     &                qratio_dx2(4,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(4,2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(4,3)*extrasum_dv(inv_ion_index,3)
                  enddo
                endif ! if(ifnr03) then...
              endif ! if(iz.eq.1) then...
            endif ! if(ifmhd_logical)...
C            monatomic H component of excitation sums.
            call exsum_component_add(
     &        ifmhd_logical, iz.eq.1, ifnr03, ifnr13,.true.,
     &        nxextrasum, nions, max_index,
     &        1, max_index,
     &        1, max_index,
     &        inv_ion,
     &        nug, nugf,
     &        nugt, nug_dv,
     &        xextrasum_scale, qratio_scale,
     &        qratio, qratiot, qratio_dx,
     &        qratiot2, qratiot_dx,
     &        qratio_dx2,
     &        qratiov, qratiovt, qratiov_dx,
     &        dqratiof, dqratiotf, dqratio_dxf,
     &        dqratiot, dqratiott, dqratio_dxt,
     &        dqratio_dv,
     &        dqratio_dxdv,
     &        dqratiovf, dqratiovt, dqratiov_dv,
     &        psum, psumf, psumt, psum_dv,
     &        ssum, ssumf, ssumt,
     &        usum,
     &        free_sum, free_sumf, free_sum_dv,
     &        xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
          endif ! if(qratio.gt.exp_lim) then
        endif ! if(ifh2.gt.0) then... neutral monatomic H change to H2
C        now do H2 change to H2 equilibrium constant
        if(ifh2.gt.0.and.ifh2plus.gt.0.and.
     &      mod(ifexcited,10).gt.1) then
          ielement = nelements + 1
          iz = 1
          ion = nions + 1
          nmin_s = nmin_species(ion)
          if(nmin_s.lt.nmin(iz).or.
     &      nmin_s.gt.nmin_max(iz))
     &      stop 'excitation_sum: invalid nmin or nmin_max'
          exparg = exparg_shift -
     &      c2t*bion(ion) - plop(ion) + qh2plus - qh2
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
C              signal to ignore this species
              exparg = -1000.d0!exparg_lim - 1.d0
            endif
          endif
!          if(exparg.gt.exparg_lim) then
            exparg = exp(exparg)
!          else
!            exparg = 0.d0
!          endif
C          qratio is first ratio of excited to ground partition function,
C          but then is transformed to ln(1 + qratio).  Similarly, there
C          is a subsequent transformation of all derivatives.
C          N.B. important convention on partial derivative variable names:
C          names starting with "qratio" are partial derivatives assuming
C          that qratio is a function of tl and x.
C          names starting with "dqratio" are the *change* to the partial
C          derivative caused by x being a function of fl, tl, dv.
C          ratio of excited to ground state partition functions
C          h2
          qratio = exparg*qstar(nmin_s,iz)
          if(qratio.gt.exp_lim) then
            qratio_us = qratio/qratio_scale
            onepq = 1.d0 + qratio_us
            onepqs = onepq*qratio_scale
C            qratiot and qratiot2 are the first and second partial wrt tl
            qratiot = qratio*(c2t*bion(ion) - plopt(ion) +
     &        qh2plust - qh2t +
     &        qstart(nmin_s,iz)/qstar(nmin_s,iz))
            qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) +
     &        qh2plust - qh2t +
     &        qstart(nmin_s,iz)/qstar(nmin_s,iz)) +
     &        qratio*(-c2t*bion(ion) - plopt2(ion) +
     &        qh2plust2 - qh2t2 +
     &        qstart2(nmin_s,iz)/qstar(nmin_s,iz) -
     &        (qstart(nmin_s,iz)/qstar(nmin_s,iz))*
     &        (qstart(nmin_s,iz)/qstar(nmin_s,iz)))
C            transform second order quantities to derivative of
C            ln(1 + qratio)
            qratiot2 = (qratiot2 - qratiot*
     &        (qratiot/onepqs))/onepq
            if(ifmhd_logical) then
C              qratio_dx(k), qratiot_dx(k), qratio_dx2(k,l) are 
C              partials wrt tl, extrasum(k) and extrasum(l) following
C              special convention for k or l = 4.
              qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +
     &          qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))
              qratiot_dx(4) = qratiot*
     &          (-occ_const*r_ion3(ion) +
     &          qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +
     &          qratio*(qstartx(5,nmin_s,iz) -
     &          qstarx(5,nmin_s,iz)*qstart(nmin_s,iz)/
     &          qstar(nmin_s,iz))/qstar(nmin_s,iz)
              qratio_dx2(4,4) = qratio_dx(4)*
     &          (-occ_const*r_ion3(ion) +
     &          qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +
     &          qratio*(qstarx2(5,5,nmin_s,iz) -
     &          qstarx(5,nmin_s,iz)*
     &          qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))/
     &          qstar(nmin_s,iz)
C              derivative of qratio wrt ln V
              qratiov = -qratio_dx(4)*extrasum(nextrasum-1)
              qratiovt = -qratiot_dx(4)*extrasum(nextrasum-1)
              qratiov_dx(4) =
     &          -qratio_dx2(4,4)*extrasum(nextrasum-1) -
     &          qratio_dx(4)
C              transform second order quantities to derivative of
C              ln(1 + qratio)
              qratiot_dx(4) = (qratiot_dx(4) - qratiot*
     &          (qratio_dx(4)/onepqs))/onepq
              qratio_dx2(4,4) = (qratio_dx2(4,4) - qratio_dx(4)*
     &          (qratio_dx(4)/onepqs))/onepq
              if(iz.eq.1) then
C                N.B. x(1) (or extrasum(4)) dependence divides out of qratio
C                n.b. indices are reordered here so that
C                qratio_dx(k) refers to derivative wrt extrasum(k) for
C                k = 1, 2, 3, and qratio_dx(4) refers to the derivative
C                wrt extrasum(nextrasum-1).
                qratio_dx(1) = qratio*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))
                qratio_dx(2) = qratio*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))
                qratio_dx(3) = qratio*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))
                qratiot_dx(1) = qratiot*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstartx(4,nmin_s,iz) -
     &            qstarx(4,nmin_s,iz)*qstart(nmin_s,iz)/
     &            qstar(nmin_s,iz))/qstar(nmin_s,iz)
                qratiot_dx(2) = qratiot*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstartx(3,nmin_s,iz) -
     &            qstarx(3,nmin_s,iz)*qstart(nmin_s,iz)/
     &            qstar(nmin_s,iz))/qstar(nmin_s,iz)
                qratiot_dx(3) = qratiot*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstartx(2,nmin_s,iz) -
     &            qstarx(2,nmin_s,iz)*qstart(nmin_s,iz)/
     &            qstar(nmin_s,iz))/qstar(nmin_s,iz)
                qratio_dx2(1,1) = qratio_dx(1)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(4,4,nmin_s,iz) -
     &            qstarx(4,nmin_s,iz)*
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(2,1) = qratio_dx(2)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(4,3,nmin_s,iz) -
     &            qstarx(4,nmin_s,iz)*
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(3,1) = qratio_dx(3)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(4,2,nmin_s,iz) -
     &            qstarx(4,nmin_s,iz)*
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(4,1) = qratio_dx(4)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(5,4,nmin_s,iz) -
     &            qstarx(5,nmin_s,iz)*
     &            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(2,2) = qratio_dx(2)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(3,3,nmin_s,iz) -
     &            qstarx(3,nmin_s,iz)*
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(3,2) = qratio_dx(3)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(3,2,nmin_s,iz) -
     &            qstarx(3,nmin_s,iz)*
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(4,2) = qratio_dx(4)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(5,3,nmin_s,iz) -
     &            qstarx(5,nmin_s,iz)*
     &            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(3,3) = qratio_dx(3)*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(2,2,nmin_s,iz) -
     &            qstarx(2,nmin_s,iz)*
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
                qratio_dx2(4,3) = qratio_dx(4)*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +
     &            qratio*(qstarx2(5,2,nmin_s,iz) -
     &            qstarx(5,nmin_s,iz)*
     &            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))/
     &            qstar(nmin_s,iz)
C                derivative of qratio wrt ln V
                qratiov = qratiov -
     &            qratio_dx(1)*extrasum(1) -
     &            qratio_dx(2)*extrasum(2) -
     &            qratio_dx(3)*extrasum(3)
                qratiovt = qratiovt -
     &            qratiot_dx(1)*extrasum(1) -
     &            qratiot_dx(2)*extrasum(2) -
     &            qratiot_dx(3)*extrasum(3)
                qratiov_dx(1) = -
     &            qratio_dx2(1,1)*extrasum(1) -
     &            qratio_dx2(2,1)*extrasum(2) -
     &            qratio_dx2(3,1)*extrasum(3) -
     &            qratio_dx2(4,1)*extrasum(nextrasum-1) -
     &            qratio_dx(1)
                qratiov_dx(2) = -
     &            qratio_dx2(2,1)*extrasum(1) -
     &            qratio_dx2(2,2)*extrasum(2) -
     &            qratio_dx2(3,2)*extrasum(3) -
     &            qratio_dx2(4,2)*extrasum(nextrasum-1) -
     &            qratio_dx(2)
                qratiov_dx(3) = -
     &            qratio_dx2(3,1)*extrasum(1) -
     &            qratio_dx2(3,2)*extrasum(2) -
     &            qratio_dx2(3,3)*extrasum(3) -
     &            qratio_dx2(4,3)*extrasum(nextrasum-1) -
     &            qratio_dx(3)
                qratiov_dx(4) = qratiov_dx(4) -
     &            qratio_dx2(4,1)*extrasum(1) -
     &            qratio_dx2(4,2)*extrasum(2) -
     &            qratio_dx2(4,3)*extrasum(3)
C                transform second order quantities to derivative of
C                ln(1 + qratio)
                qratiot_dx(1) = (qratiot_dx(1) - qratiot*
     &            (qratio_dx(1)/onepqs))/onepq
                qratiot_dx(2) = (qratiot_dx(2) - qratiot*
     &            (qratio_dx(2)/onepqs))/onepq
                qratiot_dx(3) = (qratiot_dx(3) - qratiot*
     &            (qratio_dx(3)/onepqs))/onepq
                qratio_dx2(1,1) = (qratio_dx2(1,1) - qratio_dx(1)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(2,1) = (qratio_dx2(2,1) - qratio_dx(2)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(3,1) = (qratio_dx2(3,1) - qratio_dx(3)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(4,1) = (qratio_dx2(4,1) - qratio_dx(4)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(2,2) = (qratio_dx2(2,2) - qratio_dx(2)*
     &            (qratio_dx(2)/onepqs))/onepq
                qratio_dx2(3,2) = (qratio_dx2(3,2) - qratio_dx(3)*
     &            (qratio_dx(2)/onepqs))/onepq
                qratio_dx2(4,2) = (qratio_dx2(4,2) - qratio_dx(4)*
     &            (qratio_dx(2)/onepqs))/onepq
                qratio_dx2(3,3) = (qratio_dx2(3,3) - qratio_dx(3)*
     &            (qratio_dx(3)/onepqs))/onepq
                qratio_dx2(4,3) = (qratio_dx2(4,3) - qratio_dx(4)*
     &            (qratio_dx(3)/onepqs))/onepq
                qratiov_dx(1) = (qratiov_dx(1) - qratiov*
     &            (qratio_dx(1)/onepqs))/onepq
                qratiov_dx(2) = (qratiov_dx(2) - qratiov*
     &            (qratio_dx(2)/onepqs))/onepq
                qratiov_dx(3) = (qratiov_dx(3) - qratiov*
     &            (qratio_dx(3)/onepqs))/onepq
C                transform first order quantities to derivative of
C                ln(1 + qratio)
                qratio_dx(1) = qratio_dx(1)/onepq
                qratio_dx(2) = qratio_dx(2)/onepq
                qratio_dx(3) = qratio_dx(3)/onepq
              endif
C              transform second order quantities to derivative of
C              ln(1 + qratio)
              qratiovt = (qratiovt - qratiov*
     &          (qratiot/onepqs))/onepq
              qratiov_dx(4) = (qratiov_dx(4) - qratiov*
     &          (qratio_dx(4)/onepqs))/onepq
C              transform first order quantities to derivative of
C              ln(1 + qratio)
              qratio_dx(4) = qratio_dx(4)/onepq
              qratiov = qratiov/onepq
            endif
C            transform first order quantities to derivative of
C            ln(1 + qratio)
            qratiot = qratiot/onepq
            if(qratio_us.gt.1.d-3) then
              qratio = qratio_scale*log(onepq)
            else
C              alternating series so relative error is less than first
C              missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
              qratio = qratio*
     &          (1.d0      - qratio_us*
     &          (1.d0/2.d0 - qratio_us*
     &          (1.d0/3.d0 - qratio_us*
     &          (1.d0/4.d0 - qratio_us*
     &          (1.d0/5.d0)))))
            endif
            if(ifmhd_logical) then
C              calculate derivative *correction* due to x dependence on fl
C              and tl.
              if(ifnr03) then
                dqratiof = qratio_dx(4)*extrasumf(nextrasum-1)
                dqratiot = qratio_dx(4)*extrasumt(nextrasum-1)
                dqratiotf = qratiot_dx(4)*extrasumf(nextrasum-1)
                dqratiott = qratiot_dx(4)*extrasumt(nextrasum-1)
                dqratio_dxf(4) = qratio_dx2(4,4)*
     &            extrasumf(nextrasum-1)
                dqratio_dxt(4) = qratio_dx2(4,4)*
     &            extrasumt(nextrasum-1)
                dqratiovf = qratiov_dx(4)*extrasumf(nextrasum-1)
                dqratiovt = qratiov_dx(4)*extrasumt(nextrasum-1)
              endif
              if(ifnr13) then
                do inv_ion_index = 1, max_index
                  dqratio_dv(inv_ion_index) = qratio_dx(4)*
     &              extrasum_dv(inv_ion_index,nextrasum-1)
                  dqratiov_dv(inv_ion_index) = qratiov_dx(4)*
     &              extrasum_dv(inv_ion_index,nextrasum-1)
                  dqratio_dxdv(inv_ion_index,4) = qratio_dx2(4,4)*
     &              extrasum_dv(inv_ion_index,nextrasum-1)
                enddo
              endif
              if(iz.eq.1) then
                if(ifnr03) then
                  dqratiof = dqratiof +
     &              qratio_dx(1)*extrasumf(1) +
     &              qratio_dx(2)*extrasumf(2) +
     &              qratio_dx(3)*extrasumf(3)
                  dqratiot = dqratiot +
     &              qratio_dx(1)*extrasumt(1) +
     &              qratio_dx(2)*extrasumt(2) +
     &              qratio_dx(3)*extrasumt(3)
                  dqratiotf = dqratiotf +
     &              qratiot_dx(1)*extrasumf(1) +
     &              qratiot_dx(2)*extrasumf(2) +
     &              qratiot_dx(3)*extrasumf(3)
                  dqratiott = dqratiott +
     &              qratiot_dx(1)*extrasumt(1) +
     &              qratiot_dx(2)*extrasumt(2) +
     &              qratiot_dx(3)*extrasumt(3)
                  dqratio_dxf(1) =
     &              qratio_dx2(1,1)*extrasumf(1) +
     &              qratio_dx2(2,1)*extrasumf(2) +
     &              qratio_dx2(3,1)*extrasumf(3) +
     &              qratio_dx2(4,1)*extrasumf(nextrasum-1)
                  dqratio_dxt(1) =
     &              qratio_dx2(1,1)*extrasumt(1) +
     &              qratio_dx2(2,1)*extrasumt(2) +
     &              qratio_dx2(3,1)*extrasumt(3) +
     &              qratio_dx2(4,1)*extrasumt(nextrasum-1)
                  dqratio_dxf(2) =
     &              qratio_dx2(2,1)*extrasumf(1) +
     &              qratio_dx2(2,2)*extrasumf(2) +
     &              qratio_dx2(3,2)*extrasumf(3) +
     &              qratio_dx2(4,2)*extrasumf(nextrasum-1)
                  dqratio_dxt(2) =
     &              qratio_dx2(2,1)*extrasumt(1) +
     &              qratio_dx2(2,2)*extrasumt(2) +
     &              qratio_dx2(3,2)*extrasumt(3) +
     &              qratio_dx2(4,2)*extrasumt(nextrasum-1)
                  dqratio_dxf(3) =
     &              qratio_dx2(3,1)*extrasumf(1) +
     &              qratio_dx2(3,2)*extrasumf(2) +
     &              qratio_dx2(3,3)*extrasumf(3) +
     &              qratio_dx2(4,3)*extrasumf(nextrasum-1)
                  dqratio_dxt(3) =
     &              qratio_dx2(3,1)*extrasumt(1) +
     &              qratio_dx2(3,2)*extrasumt(2) +
     &              qratio_dx2(3,3)*extrasumt(3) +
     &              qratio_dx2(4,3)*extrasumt(nextrasum-1)
                  dqratio_dxf(4) = dqratio_dxf(4) +
     &              qratio_dx2(4,1)*extrasumf(1) +
     &              qratio_dx2(4,2)*extrasumf(2) +
     &              qratio_dx2(4,3)*extrasumf(3)
                  dqratio_dxt(4) = dqratio_dxt(4) +
     &              qratio_dx2(4,1)*extrasumt(1) +
     &              qratio_dx2(4,2)*extrasumt(2) +
     &              qratio_dx2(4,3)*extrasumt(3)
                  dqratiovf = dqratiovf +
     &              qratiov_dx(1)*extrasumf(1) +
     &              qratiov_dx(2)*extrasumf(2) +
     &              qratiov_dx(3)*extrasumf(3)
                  dqratiovt = dqratiovt +
     &              qratiov_dx(1)*extrasumt(1) +
     &              qratiov_dx(2)*extrasumt(2) +
     &              qratiov_dx(3)*extrasumt(3)
                endif
                if(ifnr13) then
                  do inv_ion_index = 1, max_index
                    dqratio_dv(inv_ion_index) =
     &                dqratio_dv(inv_ion_index) +
     &                qratio_dx(1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx(2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx(3)*extrasum_dv(inv_ion_index,3)
                    dqratiov_dv(inv_ion_index) =
     &                dqratiov_dv(inv_ion_index) +
     &                qratiov_dx(1)*extrasum_dv(inv_ion_index,1) +
     &                qratiov_dx(2)*extrasum_dv(inv_ion_index,2) +
     &                qratiov_dx(3)*extrasum_dv(inv_ion_index,3)
                    dqratio_dxdv(inv_ion_index,1) =
     &                qratio_dx2(1,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(2,1)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(3,1)*extrasum_dv(inv_ion_index,3) +
     &                qratio_dx2(4,1)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,2) =
     &                qratio_dx2(2,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(2,2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(3,2)*extrasum_dv(inv_ion_index,3) +
     &                qratio_dx2(4,2)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,3) =
     &                qratio_dx2(3,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(3,2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(3,3)*extrasum_dv(inv_ion_index,3) +
     &                qratio_dx2(4,3)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,4) =
     &                dqratio_dxdv(inv_ion_index,4) +
     &                qratio_dx2(4,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(4,2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(4,3)*extrasum_dv(inv_ion_index,3)
                  enddo
                endif
              endif
            endif
C            H2 component of excitation sums.
            call exsum_component_add(
     &        ifmhd_logical, iz.eq.1, ifnr03, ifnr13,.true.,
     &        nxextrasum, nions, max_index,
     &        1, max_index,
     &        1, max_index,
     &        inv_ion,
     &        nuh2, nuh2f,
     &        nuh2t, nuh2_dv,
     &        xextrasum_scale, qratio_scale,
     &        qratio, qratiot, qratio_dx,
     &        qratiot2, qratiot_dx,
     &        qratio_dx2,
     &        qratiov, qratiovt, qratiov_dx,
     &        dqratiof, dqratiotf, dqratio_dxf,
     &        dqratiot, dqratiott, dqratio_dxt,
     &        dqratio_dv,
     &        dqratio_dxdv,
     &        dqratiovf, dqratiovt, dqratiov_dv,
     &        psum, psumf, psumt, psum_dv,
     &        ssum, ssumf, ssumt,
     &        usum,
     &        free_sum, free_sumf, free_sum_dv,
     &        xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
          endif !if(qratio.gt.exp_lim) then
C         finished H2 change to H2 equilibrium constant.
        endif
C        now do H2+
        if(ifh2.gt.0.and.ifh2plus.gt.0.and.
     &      mod(ifexcited,10).gt.1) then
          ielement = nelements + 2
          iz = 2
          izqstar = izhi + 1
          ion = nions + 2
          nmin_s = nmin_species(ion)
          if(nmin_s.lt.nmin(iz).or.
     &      nmin_s.gt.nmin_max(iz))
     &      stop 'excitation_sum: invalid nmin or nmin_max'
C          n.b. excited state core is two protons with unity statistical weight
C          (following usual convention that nuclear spin statistical weights
C          are divided out).
          exparg = exparg_shift -
     &      c2t*bion(ion) - plop(ion) - qh2plus
C          correct for MHD ground state occupation probability.
          if(ifmhd_logical) then
C            H2+ part of "neutral" list in this case.
            if(iz.eq.2) then
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
C              signal to ignore this species
              exparg = -1000.d0!exparg_lim - 1.d0
            endif
          endif
!          if(exparg.gt.exparg_lim) then
            exparg = exp(exparg)
!          else
!            exparg = 0.d0
!          endif
C          qratio is first ratio of excited to ground partition function,
C          but then is transformed to ln(1 + qratio).  Similarly, there
C          is a subsequent transformation of all derivatives.
C          N.B. important convention on partial derivative variable names:
C          names starting with "qratio" are partial derivatives assuming
C          that qratio is a function of tl and x.
C          names starting with "dqratio" are the *change* to the partial
C          derivative caused by x being a function of fl, tl, dv.
C            ratio of excited to ground state partition functions
C          h2+
          qratio = exparg*qstar(nmin_s,izqstar)
          if(qratio.gt.exp_lim) then
            qratio_us = qratio/qratio_scale
            onepq = 1.d0 + qratio_us
            onepqs = onepq*qratio_scale
C            qratiot and qratiot2 are the first and second partial wrt tl
            qratiot = qratio*(c2t*bion(ion) - plopt(ion) -
     &        qh2plust +
     &        qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar))
            qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) -
     &        qh2plust +
     &        qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &        qratio*(-c2t*bion(ion) - plopt2(ion) -
     &        qh2plust2 +
     &        qstart2(nmin_s,izqstar)/qstar(nmin_s,izqstar) -
     &        (qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar))*
     &        (qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar)))
C            transform second order quantities to derivative of
C            ln(1 + qratio)
            qratiot2 = (qratiot2 - qratiot*
     &        (qratiot/onepqs))/onepq
            if(ifmhd_logical) then
C              qratio_dx(k), qratiot_dx(k), qratio_dx2(k,l) are 
C              partials wrt tl, extrasum(k) and extrasum(l) following
C              special convention for k or l = 4.
              qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +
     &          qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))
              qratiot_dx(4) = qratiot*
     &          (-occ_const*r_ion3(ion) +
     &          qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &          qratio*(qstartx(5,nmin_s,izqstar) -
     &          qstarx(5,nmin_s,izqstar)*qstart(nmin_s,izqstar)/
     &          qstar(nmin_s,izqstar))/qstar(nmin_s,izqstar)
              qratio_dx2(4,4) = qratio_dx(4)*
     &          (-occ_const*r_ion3(ion) +
     &          qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &          qratio*(qstarx2(5,5,nmin_s,izqstar) -
     &          qstarx(5,nmin_s,izqstar)*
     &          qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &          qstar(nmin_s,izqstar)
C              derivative of qratio wrt ln V
              qratiov = -qratio_dx(4)*extrasum(nextrasum-1)
              qratiovt = -qratiot_dx(4)*extrasum(nextrasum-1)
              qratiov_dx(4) =
     &          -qratio_dx2(4,4)*extrasum(nextrasum-1) -
     &          qratio_dx(4)
C              transform second order quantities to derivative of
C              ln(1 + qratio)
              qratiot_dx(4) = (qratiot_dx(4) - qratiot*
     &          (qratio_dx(4)/onepqs))/onepq
              qratio_dx2(4,4) = (qratio_dx2(4,4) - qratio_dx(4)*
     &          (qratio_dx(4)/onepqs))/onepq
              if(iz.eq.2) then
C                N.B. x(1) (or extrasum(4)) dependence divides out of qratio
C                n.b. indices are reordered here so that
C                qratio_dx(k) refers to derivative wrt extrasum(k) for
C                k = 1, 2, 3, and qratio_dx(4) refers to the derivative
C                wrt extrasum(nextrasum-1).
                qratio_dx(1) =
     &            qratio*(-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                qratio_dx(2) =
     &            qratio*(-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                qratio_dx(3) =
     &            qratio*(-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                qratiot_dx(1) = qratiot*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstartx(4,nmin_s,izqstar) -
     &            qstarx(4,nmin_s,izqstar)*qstart(nmin_s,izqstar)/
     &            qstar(nmin_s,izqstar))/qstar(nmin_s,izqstar)
                qratiot_dx(2) = qratiot*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstartx(3,nmin_s,izqstar) -
     &            qstarx(3,nmin_s,izqstar)*qstart(nmin_s,izqstar)/
     &            qstar(nmin_s,izqstar))/qstar(nmin_s,izqstar)
                qratiot_dx(3) = qratiot*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstartx(2,nmin_s,izqstar) -
     &            qstarx(2,nmin_s,izqstar)*qstart(nmin_s,izqstar)/
     &            qstar(nmin_s,izqstar))/qstar(nmin_s,izqstar)
                qratio_dx2(1,1) = qratio_dx(1)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstarx2(4,4,nmin_s,izqstar) -
     &            qstarx(4,nmin_s,izqstar)*
     &            qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &            qstar(nmin_s,izqstar)
                qratio_dx2(2,1) = qratio_dx(2)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstarx2(4,3,nmin_s,izqstar) -
     &            qstarx(4,nmin_s,izqstar)*
     &            qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &            qstar(nmin_s,izqstar)
                qratio_dx2(3,1) = qratio_dx(3)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstarx2(4,2,nmin_s,izqstar) -
     &            qstarx(4,nmin_s,izqstar)*
     &            qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &            qstar(nmin_s,izqstar)
                qratio_dx2(4,1) = qratio_dx(4)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*r_neutral(ielement) +
     &            qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstarx2(5,4,nmin_s,izqstar) -
     &            qstarx(5,nmin_s,izqstar)*
     &            qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &            qstar(nmin_s,izqstar)
                qratio_dx2(2,2) = qratio_dx(2)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstarx2(3,3,nmin_s,izqstar) -
     &            qstarx(3,nmin_s,izqstar)*
     &            qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &            qstar(nmin_s,izqstar)
                qratio_dx2(3,2) = qratio_dx(3)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstarx2(3,2,nmin_s,izqstar) -
     &            qstarx(3,nmin_s,izqstar)*
     &            qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &            qstar(nmin_s,izqstar)
                qratio_dx2(4,2) = qratio_dx(4)*
     &            (-occ_const*r_neutral(ielement)*
     &            r_neutral(ielement)*3.d0 +
     &            qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstarx2(5,3,nmin_s,izqstar) -
     &            qstarx(5,nmin_s,izqstar)*
     &            qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &            qstar(nmin_s,izqstar)
                qratio_dx2(3,3) = qratio_dx(3)*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstarx2(2,2,nmin_s,izqstar) -
     &            qstarx(2,nmin_s,izqstar)*
     &            qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &            qstar(nmin_s,izqstar)
                qratio_dx2(4,3) = qratio_dx(4)*
     &            (-occ_const*r_neutral(ielement)*
     &            3.d0 +
     &            qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +
     &            qratio*(qstarx2(5,2,nmin_s,izqstar) -
     &            qstarx(5,nmin_s,izqstar)*
     &            qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar))/
     &            qstar(nmin_s,izqstar)
C                derivative of qratio wrt ln V
                qratiov = qratiov -
     &            qratio_dx(1)*extrasum(1) -
     &            qratio_dx(2)*extrasum(2) -
     &            qratio_dx(3)*extrasum(3)
                qratiovt = qratiovt -
     &            qratiot_dx(1)*extrasum(1) -
     &            qratiot_dx(2)*extrasum(2) -
     &            qratiot_dx(3)*extrasum(3)
                qratiov_dx(1) = -
     &            qratio_dx2(1,1)*extrasum(1) -
     &            qratio_dx2(2,1)*extrasum(2) -
     &            qratio_dx2(3,1)*extrasum(3) -
     &            qratio_dx2(4,1)*extrasum(nextrasum-1) -
     &            qratio_dx(1)
                qratiov_dx(2) = -
     &            qratio_dx2(2,1)*extrasum(1) -
     &            qratio_dx2(2,2)*extrasum(2) -
     &            qratio_dx2(3,2)*extrasum(3) -
     &            qratio_dx2(4,2)*extrasum(nextrasum-1) -
     &            qratio_dx(2)
                qratiov_dx(3) = -
     &            qratio_dx2(3,1)*extrasum(1) -
     &            qratio_dx2(3,2)*extrasum(2) -
     &            qratio_dx2(3,3)*extrasum(3) -
     &            qratio_dx2(4,3)*extrasum(nextrasum-1) -
     &            qratio_dx(3)
                qratiov_dx(4) = qratiov_dx(4) -
     &            qratio_dx2(4,1)*extrasum(1) -
     &            qratio_dx2(4,2)*extrasum(2) -
     &            qratio_dx2(4,3)*extrasum(3)
C                transform second order quantities to derivative of
C                ln(1 + qratio)
                qratiot_dx(1) = (qratiot_dx(1) - qratiot*
     &            (qratio_dx(1)/onepqs))/onepq
                qratiot_dx(2) = (qratiot_dx(2) - qratiot*
     &            (qratio_dx(2)/onepqs))/onepq
                qratiot_dx(3) = (qratiot_dx(3) - qratiot*
     &            (qratio_dx(3)/onepqs))/onepq
                qratio_dx2(1,1) = (qratio_dx2(1,1) - qratio_dx(1)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(2,1) = (qratio_dx2(2,1) - qratio_dx(2)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(3,1) = (qratio_dx2(3,1) - qratio_dx(3)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(4,1) = (qratio_dx2(4,1) - qratio_dx(4)*
     &            (qratio_dx(1)/onepqs))/onepq
                qratio_dx2(2,2) = (qratio_dx2(2,2) - qratio_dx(2)*
     &            (qratio_dx(2)/onepqs))/onepq
                qratio_dx2(3,2) = (qratio_dx2(3,2) - qratio_dx(3)*
     &            (qratio_dx(2)/onepqs))/onepq
                qratio_dx2(4,2) = (qratio_dx2(4,2) - qratio_dx(4)*
     &            (qratio_dx(2)/onepqs))/onepq
                qratio_dx2(3,3) = (qratio_dx2(3,3) - qratio_dx(3)*
     &            (qratio_dx(3)/onepqs))/onepq
                qratio_dx2(4,3) = (qratio_dx2(4,3) - qratio_dx(4)*
     &            (qratio_dx(3)/onepqs))/onepq
                qratiov_dx(1) = (qratiov_dx(1) - qratiov*
     &            (qratio_dx(1)/onepqs))/onepq
                qratiov_dx(2) = (qratiov_dx(2) - qratiov*
     &            (qratio_dx(2)/onepqs))/onepq
                qratiov_dx(3) = (qratiov_dx(3) - qratiov*
     &            (qratio_dx(3)/onepqs))/onepq
C                transform first order quantities to derivative of
C                ln(1 + qratio)
                qratio_dx(1) = qratio_dx(1)/onepq
                qratio_dx(2) = qratio_dx(2)/onepq
                qratio_dx(3) = qratio_dx(3)/onepq
              endif ! if(iz.eq.2) then...
C              transform second order quantities to derivative of
C              ln(1 + qratio)
              qratiovt = (qratiovt - qratiov*
     &          (qratiot/onepqs))/onepq
              qratiov_dx(4) = (qratiov_dx(4) - qratiov*
     &          (qratio_dx(4)/onepqs))/onepq
C              transform first order quantities to derivative of
C              ln(1 + qratio)
              qratio_dx(4) = qratio_dx(4)/onepq
              qratiov = qratiov/onepq
            endif ! if(ifmhd_logical) then
C            transform first order quantities to derivative of
C            ln(1 + qratio)
            qratiot = qratiot/onepq
            if(qratio_us.gt.1.d-3) then
              qratio = qratio_scale*log(onepq)
            else
C              alternating series so relative error is less than first
C              missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
              qratio = qratio*
     &          (1.d0      - qratio_us*
     &          (1.d0/2.d0 - qratio_us*
     &          (1.d0/3.d0 - qratio_us*
     &          (1.d0/4.d0 - qratio_us*
     &          (1.d0/5.d0)))))
            endif
            if(ifmhd_logical) then
C              calculate derivative *correction* due to x dependence on fl
C              and tl.
              if(ifnr03) then
                dqratiof = qratio_dx(4)*extrasumf(nextrasum-1)
                dqratiot = qratio_dx(4)*extrasumt(nextrasum-1)
                dqratiotf = qratiot_dx(4)*extrasumf(nextrasum-1)
                dqratiott = qratiot_dx(4)*extrasumt(nextrasum-1)
                dqratio_dxf(4) = qratio_dx2(4,4)*
     &            extrasumf(nextrasum-1)
                dqratio_dxt(4) = qratio_dx2(4,4)*
     &            extrasumt(nextrasum-1)
                dqratiovf = qratiov_dx(4)*extrasumf(nextrasum-1)
                dqratiovt = qratiov_dx(4)*extrasumt(nextrasum-1)
              endif
              if(ifnr13) then
                do inv_ion_index = 1, max_index
                  dqratio_dv(inv_ion_index) = qratio_dx(4)*
     &              extrasum_dv(inv_ion_index,nextrasum-1)
                  dqratiov_dv(inv_ion_index) = qratiov_dx(4)*
     &              extrasum_dv(inv_ion_index,nextrasum-1)
                  dqratio_dxdv(inv_ion_index,4) = qratio_dx2(4,4)*
     &              extrasum_dv(inv_ion_index,nextrasum-1)
                enddo
              endif
              if(iz.eq.2) then
                if(ifnr03) then
                  dqratiof = dqratiof +
     &              qratio_dx(1)*extrasumf(1) +
     &              qratio_dx(2)*extrasumf(2) +
     &              qratio_dx(3)*extrasumf(3)
                  dqratiot = dqratiot +
     &              qratio_dx(1)*extrasumt(1) +
     &              qratio_dx(2)*extrasumt(2) +
     &              qratio_dx(3)*extrasumt(3)
                  dqratiotf = dqratiotf +
     &              qratiot_dx(1)*extrasumf(1) +
     &              qratiot_dx(2)*extrasumf(2) +
     &              qratiot_dx(3)*extrasumf(3)
                  dqratiott = dqratiott +
     &              qratiot_dx(1)*extrasumt(1) +
     &              qratiot_dx(2)*extrasumt(2) +
     &              qratiot_dx(3)*extrasumt(3)
                  dqratio_dxf(1) =
     &              qratio_dx2(1,1)*extrasumf(1) +
     &              qratio_dx2(2,1)*extrasumf(2) +
     &              qratio_dx2(3,1)*extrasumf(3) +
     &              qratio_dx2(4,1)*extrasumf(nextrasum-1)
                  dqratio_dxt(1) =
     &              qratio_dx2(1,1)*extrasumt(1) +
     &              qratio_dx2(2,1)*extrasumt(2) +
     &              qratio_dx2(3,1)*extrasumt(3) +
     &              qratio_dx2(4,1)*extrasumt(nextrasum-1)
                  dqratio_dxf(2) =
     &              qratio_dx2(2,1)*extrasumf(1) +
     &              qratio_dx2(2,2)*extrasumf(2) +
     &              qratio_dx2(3,2)*extrasumf(3) +
     &              qratio_dx2(4,2)*extrasumf(nextrasum-1)
                  dqratio_dxt(2) =
     &              qratio_dx2(2,1)*extrasumt(1) +
     &              qratio_dx2(2,2)*extrasumt(2) +
     &              qratio_dx2(3,2)*extrasumt(3) +
     &              qratio_dx2(4,2)*extrasumt(nextrasum-1)
                  dqratio_dxf(3) =
     &              qratio_dx2(3,1)*extrasumf(1) +
     &              qratio_dx2(3,2)*extrasumf(2) +
     &              qratio_dx2(3,3)*extrasumf(3) +
     &              qratio_dx2(4,3)*extrasumf(nextrasum-1)
                  dqratio_dxt(3) =
     &              qratio_dx2(3,1)*extrasumt(1) +
     &              qratio_dx2(3,2)*extrasumt(2) +
     &              qratio_dx2(3,3)*extrasumt(3) +
     &              qratio_dx2(4,3)*extrasumt(nextrasum-1)
                  dqratio_dxf(4) = dqratio_dxf(4) +
     &              qratio_dx2(4,1)*extrasumf(1) +
     &              qratio_dx2(4,2)*extrasumf(2) +
     &              qratio_dx2(4,3)*extrasumf(3)
                  dqratio_dxt(4) = dqratio_dxt(4) +
     &              qratio_dx2(4,1)*extrasumt(1) +
     &              qratio_dx2(4,2)*extrasumt(2) +
     &              qratio_dx2(4,3)*extrasumt(3)
                  dqratiovf = dqratiovf +
     &              qratiov_dx(1)*extrasumf(1) +
     &              qratiov_dx(2)*extrasumf(2) +
     &              qratiov_dx(3)*extrasumf(3)
                  dqratiovt = dqratiovt +
     &              qratiov_dx(1)*extrasumt(1) +
     &              qratiov_dx(2)*extrasumt(2) +
     &              qratiov_dx(3)*extrasumt(3)
                endif
                if(ifnr13) then
                  do inv_ion_index = 1, max_index
                    dqratio_dv(inv_ion_index) =
     &                dqratio_dv(inv_ion_index) +
     &                qratio_dx(1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx(2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx(3)*extrasum_dv(inv_ion_index,3)
                    dqratiov_dv(inv_ion_index) =
     &                dqratiov_dv(inv_ion_index) +
     &                qratiov_dx(1)*extrasum_dv(inv_ion_index,1) +
     &                qratiov_dx(2)*extrasum_dv(inv_ion_index,2) +
     &                qratiov_dx(3)*extrasum_dv(inv_ion_index,3)
                    dqratio_dxdv(inv_ion_index,1) =
     &                qratio_dx2(1,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(2,1)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(3,1)*extrasum_dv(inv_ion_index,3) +
     &                qratio_dx2(4,1)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,2) =
     &                qratio_dx2(2,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(2,2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(3,2)*extrasum_dv(inv_ion_index,3) +
     &                qratio_dx2(4,2)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,3) =
     &                qratio_dx2(3,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(3,2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(3,3)*extrasum_dv(inv_ion_index,3) +
     &                qratio_dx2(4,3)*
     &                extrasum_dv(inv_ion_index,nextrasum-1)
                    dqratio_dxdv(inv_ion_index,4) =
     &                dqratio_dxdv(inv_ion_index,4) +
     &                qratio_dx2(4,1)*extrasum_dv(inv_ion_index,1) +
     &                qratio_dx2(4,2)*extrasum_dv(inv_ion_index,2) +
     &                qratio_dx2(4,3)*extrasum_dv(inv_ion_index,3)
                  enddo
                endif ! if(ifnr03) then
              endif ! if(iz.eq.2) then
            endif ! if(ifmhd_logical) then
C            H2plus component of excitation sums.
            call exsum_component_add(
     &        ifmhd_logical, iz.eq.2, ifnr03, ifnr13,.true.,
     &        nxextrasum, nions, max_index,
     &        1, max_index,
     &        1, max_index,
     &        inv_ion,
     &        nuh2plus, nuh2plusf,
     &        nuh2plust, nuh2plus_dv,
     &        xextrasum_scale, qratio_scale,
     &        qratio, qratiot, qratio_dx,
     &        qratiot2, qratiot_dx,
     &        qratio_dx2,
     &        qratiov, qratiovt, qratiov_dx,
     &        dqratiof, dqratiotf, dqratio_dxf,
     &        dqratiot, dqratiott, dqratio_dxt,
     &        dqratio_dv,
     &        dqratio_dxdv,
     &        dqratiovf, dqratiovt, dqratiov_dv,
     &        psum, psumf, psumt, psum_dv,
     &        ssum, ssumf, ssumt,
     &        usum,
     &        free_sum, free_sumf, free_sum_dv,
     &        xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
          endif ! if(qratio.gt.exp_lim) then
        endif ! if(ifh2.gt.0.and.ifh2plus.gt.0.and.... now do H2+
      endif ! if(partial_element.... non-zero H abundance
      end
      subroutine exsum_component_add(
     &  ifmhd_logical, iz_logical, ifnr03, ifnr13, ifsimple_index,
     &  nxextrasum, nions, max_dvindex,
     &  ion_start1, ion_end1,
     &  ion_start2, ion_end2,
     &  inv_ion, 
     &  nuvar, nuvarf,
     &  nuvart, nuvar_dv,
     &  xextrasum_scale, qratio_scale,
     &  qratio, qratiot, qratio_dx,
     &  qratiot2, qratiot_dx,
     &  qratio_dx2,
     &  qratiov, qratiovt, qratiov_dx,
     &  dqratiof, dqratiotf, dqratio_dxf,
     &  dqratiot, dqratiott, dqratio_dxt,
     &  dqratio_dv,
     &  dqratio_dxdv,
     &  dqratiovf, dqratiovt, dqratiov_dv,
     &  psum, psumf, psumt, psum_dv,
     &  ssum, ssumf, ssumt,
     &  usum,
     &  free_sum, free_sumf, free_sum_dv,
     &  xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
C      add components to each excitation sum and its derivatives.
      implicit none
C      input variables:
      integer nxextrasum, nions, max_dvindex
      integer ion_start1, ion_end1, ion_start2, ion_end2,
     &  inv_ion(nions+2)
      logical ifmhd_logical, iz_logical, ifnr03, ifnr13,
     &  ifsimple_index
      double precision nuvar, nuvarf, nuvart, nuvar_dv(max_dvindex),
     &  xextrasum_scale(nxextrasum), qratio_scale,
     &  qratio, qratiot, qratio_dx(nxextrasum),
     &  qratiot2, qratiot_dx(nxextrasum),
     &  qratio_dx2(nxextrasum,nxextrasum),
     &  qratiov, qratiovt, qratiov_dx(nxextrasum),
     &  dqratiof, dqratiotf, dqratio_dxf(nxextrasum),
     &  dqratiot, dqratiott, dqratio_dxt(nxextrasum),
     &  dqratio_dv(nions+2),
     &  dqratio_dxdv(nions+2,nxextrasum),
     &  dqratiovf, dqratiovt, dqratiov_dv(nions+2)
C      output variables:
      double precision
     &  psum, psumf, psumt, psum_dv(nions+2),
     &  ssum, ssumf, ssumt,
     &  usum,
     &  free_sum, free_sumf, free_sum_dv(nions+2),
     &  xextrasum(nxextrasum), xextrasumf(nxextrasum),
     &  xextrasumt(nxextrasum), xextrasum_dv(nions+2,nxextrasum)
C      local variables:
      integer ion_index, ion_index1, inv_ion_index 
      double precision
     &  nuvarsum_scale, nuvarfsum_scale,
     &  nuvartsum_scale, nuvar_dvsum_scale,
     &  nuvar_scale(4), nuvarf_scale(4),
     &  nuvart_scale(4), nuvar_dv_scale(4)
      nuvarsum_scale = nuvar/qratio_scale
      nuvar_scale(1) = nuvar*(xextrasum_scale(1)/qratio_scale)
      nuvar_scale(2) = nuvar*(xextrasum_scale(2)/qratio_scale)
      nuvar_scale(3) = nuvar*(xextrasum_scale(3)/qratio_scale)
      nuvar_scale(4) = nuvar*(xextrasum_scale(4)/qratio_scale)
      nuvarfsum_scale = nuvarf/qratio_scale
      nuvarf_scale(1) = nuvarf*(xextrasum_scale(1)/qratio_scale)
      nuvarf_scale(2) = nuvarf*(xextrasum_scale(2)/qratio_scale)
      nuvarf_scale(3) = nuvarf*(xextrasum_scale(3)/qratio_scale)
      nuvarf_scale(4) = nuvarf*(xextrasum_scale(4)/qratio_scale)
      nuvartsum_scale = nuvart/qratio_scale
      nuvart_scale(1) = nuvart*(xextrasum_scale(1)/qratio_scale)
      nuvart_scale(2) = nuvart*(xextrasum_scale(2)/qratio_scale)
      nuvart_scale(3) = nuvart*(xextrasum_scale(3)/qratio_scale)
      nuvart_scale(4) = nuvart*(xextrasum_scale(4)/qratio_scale)
C      psum = sum n(species)/(rho*avogadro) 
C        partial delta ln Z wrt ln V
C      psumf and psumt are derivatives (including dependence of n/rho
C      on fl and tl) wrt fl and tl of psum.
C      ssum = sum n(species)/(rho*avogadro) (delta ln Z +
C        partial delta ln Z wrt ln t)
C      ssumf and ssumt are derivatives (including dependence of n/rho
C      on fl and tl) wrt fl and tl of ssum.
C      usum = sum n(species)/(rho*avogadro) 
C      partial delta ln Z wrt ln t
      ssum = ssum + nuvarsum_scale*(qratio + qratiot)
      usum = usum + nuvarsum_scale*qratiot
      free_sum = free_sum + nuvarsum_scale*qratio
      if(ifmhd_logical) then
        psum = psum + nuvarsum_scale*qratiov
        if(iz_logical) then
          xextrasum(1) = xextrasum(1) -
     &      nuvar_scale(1)*qratio_dx(1)
          xextrasum(2) = xextrasum(2) -
     &      nuvar_scale(2)*qratio_dx(2)
          xextrasum(3) = xextrasum(3) -
     &      nuvar_scale(3)*qratio_dx(3)
        endif
        xextrasum(4) = xextrasum(4) -
     &    nuvar_scale(4)*qratio_dx(4)
      endif
      if(ifnr03) then
        ssumf = ssumf + nuvarfsum_scale*(qratio + qratiot)
        ssumt = ssumt + nuvartsum_scale*(qratio + qratiot) +
     &    nuvarsum_scale*(qratiot + qratiot2)
        free_sumf = free_sumf + nuvarfsum_scale*qratio
        if(ifmhd_logical) then
          psumf = psumf + nuvarfsum_scale*qratiov
          psumt = psumt + nuvartsum_scale*qratiov +
     &      nuvarsum_scale*qratiovt
          if(iz_logical) then
            xextrasumf(1) = xextrasumf(1) -
     &        nuvarf_scale(1)*qratio_dx(1) -
     &        nuvar_scale(1)*dqratio_dxf(1)
            xextrasumf(2) = xextrasumf(2) -
     &        nuvarf_scale(2)*qratio_dx(2) -
     &        nuvar_scale(2)*dqratio_dxf(2)
            xextrasumf(3) = xextrasumf(3) -
     &        nuvarf_scale(3)*qratio_dx(3) -
     &        nuvar_scale(3)*dqratio_dxf(3)
            xextrasumt(1) = xextrasumt(1) -
     &        nuvart_scale(1)*qratio_dx(1) -
     &        nuvar_scale(1)*(qratiot_dx(1) +
     &        dqratio_dxt(1))
            xextrasumt(2) = xextrasumt(2) -
     &        nuvart_scale(2)*qratio_dx(2) -
     &        nuvar_scale(2)*(qratiot_dx(2) +
     &        dqratio_dxt(2))
            xextrasumt(3) = xextrasumt(3) -
     &        nuvart_scale(3)*qratio_dx(3) -
     &        nuvar_scale(3)*(qratiot_dx(3) +
     &        dqratio_dxt(3))
          endif
          xextrasumf(4) = xextrasumf(4) -
     &      nuvarf_scale(4)*qratio_dx(4) -
     &      nuvar_scale(4)*dqratio_dxf(4)
          xextrasumt(4) = xextrasumt(4) -
     &      nuvart_scale(4)*qratio_dx(4) -
     &      nuvar_scale(4)*(qratiot_dx(4) +
     &      dqratio_dxt(4))
          psumf = psumf + nuvarsum_scale*dqratiovf
          psumt = psumt + nuvarsum_scale*dqratiovt
          ssumf = ssumf + nuvarsum_scale*
     &      (dqratiof + dqratiotf)
          ssumt = ssumt + nuvarsum_scale*
     &      (dqratiot + dqratiott)
          free_sumf = free_sumf + nuvarsum_scale*dqratiof
        endif
      endif
      if(ifnr13) then
        do ion_index = ion_start1, ion_end1
          if(ifsimple_index) then
            inv_ion_index = ion_index
            ion_index1 = inv_ion_index
          else
            inv_ion_index = inv_ion(ion_index)
            ion_index1 = ion_index - ion_start1 + 1
          endif
          nuvar_dvsum_scale = nuvar_dv(ion_index1)/qratio_scale
          free_sum_dv(inv_ion_index) =
     &      free_sum_dv(inv_ion_index) + nuvar_dvsum_scale*qratio
        enddo
      endif
      if(ifmhd_logical.and.ifnr13) then
        do ion_index = ion_start1, ion_end1
          if(ifsimple_index) then
            inv_ion_index = ion_index
            ion_index1 = inv_ion_index
          else
            inv_ion_index = inv_ion(ion_index)
            ion_index1 = ion_index - ion_start1 + 1
          endif
          nuvar_dvsum_scale = nuvar_dv(ion_index1)/qratio_scale
          nuvar_dv_scale(1) = nuvar_dv(ion_index1)*
     &      (xextrasum_scale(1)/qratio_scale)
          nuvar_dv_scale(2) = nuvar_dv(ion_index1)*
     &      (xextrasum_scale(2)/qratio_scale)
          nuvar_dv_scale(3) = nuvar_dv(ion_index1)*
     &      (xextrasum_scale(3)/qratio_scale)
          nuvar_dv_scale(4) = nuvar_dv(ion_index1)*
     &      (xextrasum_scale(4)/qratio_scale)
          psum_dv(inv_ion_index) = psum_dv(inv_ion_index) +
     &      nuvar_dvsum_scale*qratiov
          if(iz_logical) then
            xextrasum_dv(inv_ion_index,1) =
     &        xextrasum_dv(inv_ion_index,1) -
     &        nuvar_dv_scale(1)*qratio_dx(1)
            xextrasum_dv(inv_ion_index,2) =
     &        xextrasum_dv(inv_ion_index,2) -
     &        nuvar_dv_scale(2)*qratio_dx(2)
            xextrasum_dv(inv_ion_index,3) =
     &        xextrasum_dv(inv_ion_index,3) -
     &        nuvar_dv_scale(3)*qratio_dx(3)
          endif
          xextrasum_dv(inv_ion_index,4) =
     &      xextrasum_dv(inv_ion_index,4) -
     &      nuvar_dv_scale(4)*qratio_dx(4)
        enddo
        do inv_ion_index = ion_start2, ion_end2
          free_sum_dv(inv_ion_index) =
     &      free_sum_dv(inv_ion_index) +
     &      nuvarsum_scale*dqratio_dv(inv_ion_index)
          psum_dv(inv_ion_index) =
     &      psum_dv(inv_ion_index) +
     &      nuvarsum_scale*dqratiov_dv(inv_ion_index)
          if(iz_logical) then
            xextrasum_dv(inv_ion_index,1) =
     &        xextrasum_dv(inv_ion_index,1) -
     &        nuvar_scale(1)*dqratio_dxdv(inv_ion_index,1)
            xextrasum_dv(inv_ion_index,2) =
     &        xextrasum_dv(inv_ion_index,2) -
     &        nuvar_scale(2)*dqratio_dxdv(inv_ion_index,2)
            xextrasum_dv(inv_ion_index,3) =
     &        xextrasum_dv(inv_ion_index,3) -
     &        nuvar_scale(3)*dqratio_dxdv(inv_ion_index,3)
          endif
          xextrasum_dv(inv_ion_index,4) =
     &      xextrasum_dv(inv_ion_index,4) -
     &      nuvar_scale(4)*dqratio_dxdv(inv_ion_index,4)
        enddo
      endif
      end
