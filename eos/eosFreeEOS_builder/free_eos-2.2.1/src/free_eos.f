C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: free_eos.f 830 2008-06-29 20:34:26Z airwin $
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
      subroutine free_eos(ifoption, ifmodified,
     &  ifion, kif_in, eps, neps, match_variable, tl, fl,
     &  t, rho, rl, p, pl, cf, cp, qf, qp, sf, st, grada, rtp,
     &  rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,
     &  gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e,
     &  degeneracy, pressure, density, energy, enthalpy, entropy,
     &  iteration_count)
C*******input flags:
C       ifoption = 1 (pteh style)
C       ifoption = 2 (fjs style)
C       ifoption = 3 (mdh style)
C       ifoption = 4 (variant of mdh style to fit SCVH)
C       ifoption = 5 (no explicit pteh, fjs, or mdh pressure ionization,
C              just_use_Planck-Larkin)
C       ifmodified = 0 (original formulation suitable for stellar interiors
C         calculation in particular style)
C       ifmodified = 1 (modified version suitable for stellar interiors calculation.
C       ifmodified = 2 (modify style for best fit of opal tables + extension)
C       ifmodified = many other values for a number of test cases below.
C       ifion is a flag that controls the way that ionization is done.  In
C         general, the lower ifion, the slower the code, and the more
C         ionization details that are calculated.
C       ifion = -2 sets lw=4 and ifmtrace = 0, which implies all 295 ionization
C         stages of the 20 elements are treated in detail.
C       ifion = -1 sets lw=3 and ifmtrace = 0, which implies minor metals
C         are treated as fully ionized while H, He, and the major
C         metals are treated in detail.
C       ifion = 0 (recommended) sets lw = 3, and ifmtrace = 0 below T = 1.d6
C         and ifmtrace = 1 (both major and minor metals treated as
C         fully ionized) above T = 1.d6.
C       ifion = 1 sets lw = 3, and ifmtrace = 1, which implies all major
C         and minor metals are always treated as fully ionized
C       ifion = 2 sets lw = 2, and ifmtrace = 1, which implies all elements
C         are treated as fully ionized.
C       n.b. major metals controlled by array iftracemetal in free_eos_detailed.
C         currently list includes C, N, O, Ne, Mg, Si, S, Fe
C       kif_in = 0, ln f = match_variable and tl are independent variables.
C       kif_in = 1, ln p = match_variable and tl are independent variables.
C       kif_in = -1 signal kif = 1, and ifrad = 2 (see calculated flags below).
C       kif_in = 2, ln rho = match_variable and tl are independent variables.
C       n.b. match_variable and tl are *always* the independent variables,
C       and it is the calling routines responsibility to place the correct
C       value consistent with kif_in in match_variable for each call.
C       the fl, pl, and rl values in the argument list are used *strictly for
C       output*
C       n.b. if kif_in is not equal to 0, then much more computer time is
C       required to compute the EOS because an fl iteration must
C       be used to to match the match_variable.  However, the initial guess
C       for fl is improved by a Taylor series approach in this case to
C       reduce these fl iterations to a minimum.
C*******input real quantities which are not flags:
C       eps(nelements_in) is an array of neps = nelements = 20 values of relative
C         abundance by weight divided by the appropriate atomic weight.
C         We_use_the atomic weight scale where (un-ionized) C(12)
C         has a weight of 12.00000000....  All weights are for the 
C         un-ionized element. The eps value for an element
C         should be the sum of the individual isotopic eps values for
C         that element.  The eps array refers to the
C         elements in the following order:
C         H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni
C       match_variable is matched by iterative adjustment of fl when kif > 0.
C       for kif = 0, match_variable is interpreted as fl
C       tl = ln t
C*******calculated flags used in free_eos_detailed call:
C       kif is ordinarily equal to kif_in, but the kif_in = -1 case is
C         used to signal kif = 1, and ifrad = 2. 
C       ifh2 = 0  no h2
C       ifh2 = 1  vdb h2
C       ifh2 = 2  st h2
C       ifh2 = 3  irwin h2
C       ifh2plus = 0  no h2plus
C       ifh2plus = 1  st h2plus
C       ifh2plus = 2  irwin h2plus
C       morder = 3, 5, or 8  (or 13, 15, or 18) means 3rd, 5th, or 8th order
C       fermi_dirac integral approximation following EFF fit (or modified
C       version of EFF fit which reduces to Cody-Thacher approximation for
C       low relativistic correction).
C       morder = -3, -5, or -8 (or -13, -15, or -18) means_use_above
C       approximations in non-relativistic limit.
C       morder = 1 uses Cody-Thacher approximation directly.
C       morder = 21 calculate Fermi-Dirac integrals with slow, but precise
C       (~1.d-9 relative errors) numerical integration.
C       morder = -21 is same as 21 in non-relativistic limit.
C       morder = 23 means original 3rd order eff result
C       morder = -23 is same as 23 in non-relativistic limit.
C      if abs(ifexchange_in) > 100 then_use_linear transform approximation
C        for exchange treatement.  Otherwise,_use_numerical transform as
C        described in Paper IV.
C      ifexchange = mod(ifexchange_in,100):
C      One other wrinkle on deciding ifexchange is that large degeneracy
C      approximations (mod(ifexchange,10) = 2) only allowed above psi_lim.

C      ifexchange:
C      To understand ifexchange, must summarize free-energy model of fex used
C      in research note.  In general from kapusta relation,
C      fex is proportional to I - J + 2^1.5 pi^2/3 beta^2 K, where
C      K and J and known integrals which are functions of psi and beta,
C      and I = K^2. (see Kovetz et al 1972, ApJ 174, 109, hereafter KLVH)
C       0  < ifexchange < 10 --> KLVH treatment (i.e., drop K term).
C      10 <= ifexchange < 20 --> Kapusta treatment (i.e., retain K term).
C      20 <= ifexchange < 30 --> special test case,_use_I term alone.
C      30 <= ifexchange < 40 --> special test case,_use_J term alone.
C      40 <= ifexchange < 50 --> special test case,_use_K term alone.
C      N.B. if ifexchange is negative, then_use_non-relativistic limit
C       of corresponding positive ifexchange option.
C      ifexchange details:
C      ifexchange = 1 --> G(psi) + weak relativistic correction from KLVH
C      ifexchange = 2 --> I-J degenerate expression from KLVH (with corrected
C        sign error on a2 and high numerical precision a2 and a3, see
C        paper IV)
C      ifexchange = 11 or 12 K term added (both series from CG).
C      ifexchange = 21 or 22 I term alone (both series from KLVH)
C      ifexchange = 31 or 32 -J term alone (both series from KLVH).
C      ifexchange = 41 or 42 K term alone.
C      mod(ifexchange,10) = 4 is lowest order fit of J, K
C      mod(ifexchange,10) = 5 is next higher order fit of J, K
C      mod(ifexchange,10) = 6 is highest order fit of J, K

C       ifcoulomb > 9, do diffraction correction
C       ifcoulomb > 0 means do metal Coulomb contribution to sum0 and sum2 exactly
C       -10 < ifcoulomb < 0 means do metal Coulomb contribution to sum0 and 
C         sum2 using the "metal Coulomb" approximation.
C         n.b. this approximation is internally replaced by PTEH 
C         approximation (next line) if either eps(1) or eps(2) are zero.
C       ifcoulomb < -9 means_use_PTEH approximation for sum0 and sum2.
C         n.b. if this approximation is combined with the PTEH 
C           approximation for pressure ionization, then the 
C           routine is much faster because no ionization fraction
C           iterations are required.
C       mod(|ifcoulomb|,10) = 0 ignore Coulomb interaction.
C       mod(|ifcoulomb|,10) = 1_use_Debye-Huckel Coulomb approximation.
C       mod(|ifcoulomb|,10) = 2_use_Debye-Huckel Coulomb approximation with tau(x) correction
C       mod(|ifcoulomb|,10) = 3_use_PTEH Coulomb approximation with their theta_e
C       mod(|ifcoulomb|,10) = 4_use_PTEH Coulomb approximation with fermi-dirac theta_e
C       mod(|ifcoulomb|,10) = 9 same as 4 with DeWitt definition of lambda
C        (using sum0a = sum0 + ne*theta_e).
C       mod(|ifcoulomb|,10) = 5_use_DH smoothly connected to modified OCP
C         with DeWitt definition of lambda.
C       mod(|ifcoulomb|,10) = 6 same as 5 with alternative smooth connection
C       mod(|ifcoulomb|,10) = 7 DH (Gamma < 1) or OCP using new DeWitt lambda.
C       mod(|ifcoulomb|,10) = 8 same as 7 with theta_e = 0.
C       
C       ifpi contains meaning of two flags:
C       ifpi > 0 means_use_Planck-Larkin occupation probability, otherwise not.
C       remaining meaning in absolute value of ifpi
C       |ifpi| = 0,_use_no pressure ionization
C       |ifpi| = 1,_use_pteh pressure ionization
C       |ifpi| = 2,_use_fjs pressure ionization
C       |ifpi| = 3,_use_MDH pressure ionization
C       |ifpi| = 4,_use_Saumon-like variation of MDH pressure ionization
C       |ifpi| > 4 same as zero, i.e.,_use_no pressure ionization.
C       ifrad = 0, no radiation pressure,
C         input match_variable is consistent (ln P excluding radiation
C         pressure) for kif = 1.
C       ifrad = 1, radiation pressure included,
C         input match_variable is consistent (ln P including radiation
C         pressure) for kif = 1.
C       ifrad = 2, radiation pressure included,
C         input match_variable is ln(ptotal-prad) for kif = 1, but
C         all output quantities are calculated with radiation included.
C         this feature is used to reduce significance loss in regions which
C         are dominated by radiation pressure.
C         for kif = 2, the input match_variable is ln rho as per normal, but
C         the output pressure(1) is ln(ptotal-prad).
C         n.b. this latter case is only used for some tables which have
C         rho and T as the independent variable and all quantities including
C         pressure derivatives calculated with radiation pressure *except for*
C         the pressure itself.  Also note for this latter case
C         (ifrad = 2, kif =2) that pressure(1) = ln (ptotal-prad) will be 
C         inconsistent with pressure(2) and pressure(3) which will be
C         partials of ln ptotal wrt ln rho and ln t.
C       lw < 3 every element treated as fully ionized.
C       lw = 3 trace metals treated as fully ionized.
C       lw > 3 very slow option with all elements treated as partially ionized.
C       iftc =1 only used for thermodynamic consistency tests on entropy(2)
C       (in which case entropy(2) returned via rtp).
C*******output quantities:
C       fl = EFF degeneracy parameter.  For kif=0 this is merely determined
C         from fl = match_variable.  For kif>0 this is determined
C         from a Taylor series approach as an initial guess, then
C         refined through iteration inside free_eos_detailed.
C       t = temperature
C       rho = density
C       rl = ln rho
C       p = pressure (total if ifrad = 1, prad subtracted if ifrad = 0)
C       pl = ln p
C       cf = cp (- partial ln rho(T,P)/partial T)^{1/2}
C       cp = specific heat at constant pressure
C       qp and qf are legacy parameters which were used to calculate an
C         approximate gravitational energy generation rate when
C         abundances were changing.  These parameters have now been
C         disabled (set to 1.d300) since the exact gravitational energy
C         generation rate for the case of changing abundances can be
C         calculated using the precepts in Strittmatter, P.A.,
C         Faulkner, J.. Robertson, J.W., and Faulkner, D.J. 1970,
C         "A Question of Entropy", Ap.J. 161, 369-373.
C         http://adsbit.harvard.edu/cgi-bin/ (join URL to next line)
C         nph-iarticle_query?bibcode=1970ApJ...161..369S
C         From eq. 10a of that paper the preferred rate form is eps_grav =
C         -(eps_nuclear - d L/dm - eps_neutrino) = - du/dt - P d (1/rho)/dt,
C         where p, rho, and energy(1) (=u) are variables already returned by
C         free_eos.  That paper also gives an alternative form (eq. 10b) of
C         the rate equation which is equivalent to Cox and Giuli eq.
C         17.75''', but that form is deprecated since it will be invalid for
C         discontinuous abundances and is much more complicated than eq.
C         10a.
C       sf = partial entropy(fl, tl)/partial fl
C       st = partial entropy(fl, tl)/partial tl
C       grada = the adiabatic temperature gradient
C       rtp = - partial ln rho(T,P) wrt ln T.  n.b. note negative sign so
C         quantity should be positive for normal EOS.
C       rmue = rho/mu_e, where mu_e is the mean molecular 
C         weight per electron.
C       fh2 = n(H+)/n(all H)
C       fhe2 = n(He+)/n(all He)
C       fhe3 = n(He++)/n(all He)
C       xmu1 = 1/mu, where mu is the mean molecular
C         weight per particle.
C       xmu3 = 1/mu_e, where mu is the mean molecular
C         weight per electron.
C       eta = degeneracy parameter (Cox and Guili eta), related
C       to EFF ln f = fl by
C        wf = sqrt(1.d0 + f) = d eta/d fl
C        eta = fl+2.d0*(wf-log(1.d0+wf))
C       gamma1-gamma3 are the gammas as defined in Cox and Guili
C       h2rat = 2*n(H2)/n(all H)
C       h2plusrat = 2*n(H2+)/n(all H)
C       lambda = Coulomb interaction parameter (two definitions depending on
C         ifcoulomb).
C       gamma_e = Coulomb diffraction parameter.
C       iteration_count (if positive) =
C         the total number of ionization fraction loops completed
C         to get the solution
C       iteration_count (if negative) =
C         negative of info status variable returned by free_eos_detailed
C         (only if that info is non-zero, i.e., it signals an abormal end).
C       degeneracy, pressure, density, energy, enthalpy, and entropy
C         are all 3-vectors, with the second component being the derivative
C         of the first component wrt match_variable (except for the case where
C         kif_in = -1), and the third component being the derivative of
C         the first component wrt tl.
C       definitions:
C         degeneracy(1) is EFF degeneracy parameter ln f defined above.
C         pressure(1) = ln pressure (except for kif=2, ifrad=2).
C         density(1) = ln density.
C         energy(1) = internal energy per unit mass.
C         enthalpy(1) = enthalpy per unit mass = energy(1) + p/rho.
C         entropy(1) = entropy per unit mass.
C*******internal quantities:
C       fm = partial fl(tl, match_variable)/partial match_variable (kif>0)
C       ft = partial fl(tl, match_variable)/partial tl (kif>0)
C*******
      implicit none
      integer ifoption, ifmodified, ifion,
     &  ifh2, ifh2plus, morder, ifexchange_in, 
     &  ifmtrace, 
     &  ifcoulomb, ifpi, ifrad, lw, ifexcited, nmax,
     &  ifreducedmass, iftc,
     &  kif_in, kif, neps, iteration_count, info
      double precision
     &  eps(neps), match_variable, fl, tl, fl_old, fm, ft,
     &  t, rho, rl, p, pl, cf, cp, qf, qp, sf, st, grada, rtp,
     &  rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,
     &  gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e,
     &  degeneracy(3), pressure(3), density(3), energy(3), enthalpy(3),
     &  entropy(3),  match_variable_old, tl_old, dm, dt, tllim
C       must be ridiculous values
      data match_variable_old, tl_old, fl_old/3*1.d30/
C       must be non-zero
      data fm, ft/2*1.d0/
      integer ncall, ncall_abrupt, ncall_start
      data ncall, ncall_abrupt, ncall_start /2*0,1/
      logical ifcheck_abrupt
      data ifcheck_abrupt/.true./
C      This variable should be changed to .true. if any of the debug
C      options are being tried in free_eos_detailed.
      logical debug_any
      parameter(debug_any = .false.)
      save
C       takes care of kif_in = -1 case.
      kif = iabs(kif_in)
C       temperature limit for a number of options
      tllim = log(1.d6)
C       starting (or final when kif=0) value for fl:
      if(kif.eq.0) then
        fl = match_variable
      elseif(kif.eq.1.or.kif.eq.2) then
        dm = match_variable - match_variable_old
        dt = tl - tl_old
C        make sure change in match_variable and tl is not too large
C        for first-order Taylor series.  This is my guess of
C        the larger limit allowed, but if you run into trouble
C        you can always reduce this limit at the expense of starting
C        with a safe fl initial value which then will require a
C        large number of iterations to refine.
C        the following values are slightly larger than the opal
C        grid spacing: delta log10 rho = .25d0, delta log10 t ~ 0.1d0.
C        n.b. maximum delta tl, fl larger than warm start criterion
C        in free_eos_detailed, but cold start in free_eos_detailed
C        is reliable and more efficient than cold start for "safe" fl
C        which then must be iterated a lot more times.
        if(abs(dm).le.0.6d0.and.abs(dt).le.0.25d0) then
          fl = fl_old + fm*dm + ft*dt
        elseif(abs(fm*dm + ft*dt).le.20.d0) then
C           if not a grotesquely large step, then decrementing
C           by 20 is a safe option.
          fl = max(-100.d0, fl_old - 20.d0)
        else
C           this is a safe, but inefficient starting value
          fl = -100.d0
        endif
C        must always specify same starting fl if doing derivative tests
C        in free_eos_detailed.
        if(debug_any) fl = -100.d0
      else
        stop 'free_eos: bad kif value'
      endif
      if((2.le.ifoption.and.ifoption.le.4).and.ifcheck_abrupt) then
C         check for abrupt changes in variables and warn if this
C         occurs too often when auxiliary variable iterations
C         are required
        ncall = ncall + 1
        if(abs(match_variable-match_variable_old).gt.0.25d0.or.
     &      abs(tl-tl_old).gt.0.05d0) then
          ncall_abrupt = ncall_abrupt + 1
        elseif(ncall_start.eq.1) then
C          Start real counting after initial cold start with lots of
C          abrupt changes.
          ncall = 0
          ncall_abrupt = 0
          ncall_start = 0
        endif
C        sample first 100 real calls after initial cold start to see
C        what fraction have excessively abrupt changes.
        if(ncall.eq.100) then
          if(ncall_abrupt.ge.10) then
            ifcheck_abrupt = .false.
            write(0,*)
     &        'free_eos warning: too many cold starts of the eos'
            write(0,*) 'smaller step sizes corresponding to'
            write(0,*)
     &        '|delta match_variable| < 0.25 and |delta tl| < 0.05'
            write(0,*) 'will make the eos run *much* faster per call'
          endif
          ncall = 0
          ncall_abrupt = 0
        endif
      endif
      tl_old = tl
      match_variable_old = match_variable
C       set ifmtrace and lw according to ifion
      if(ifion.eq.-2) then
C         all ionization stages of all elements treated in detail.
        lw = 4
        ifmtrace = 0
      elseif(ifion.eq.-1) then
C         minor metals approximated as fully ionized.
        lw = 3
        ifmtrace = 0
      elseif(ifion.eq.0) then
        lw = 3
        if(tl.lt.tllim) then
C           only minor metals approximated as fully ionized.
          ifmtrace = 0
        else
C           minor and major metals approximated as fully ionized.
          ifmtrace = 1
        endif
      elseif(ifion.eq.1) then
C         minor and major metals approximated as fully ionized.
        lw = 3
        ifmtrace = 1
      elseif(ifion.eq.2) then
C         all elements approximated as fully ionized.
        lw = 2
        ifmtrace = 1
      else
        stop 'free_eos: bad ifion flag'
      endif
C       default behaviour unless otherwise specified below
      ifh2 = 3
      ifh2plus = 2
      morder = 13
      if(ifmodified.le.0) then
        ifexchange_in = 0
      else
C        Non-linear transform, Kapusta term, lowest-order EFF-style
C        approximation for J and K integrals.
        ifexchange_in = 14
      endif
C       default treatment of radiation pressure.
      if(kif_in.eq.-1) then
        ifrad = 2
      else
        ifrad = 1
      endif
      if(lw.le.2.or.(ifmodified.le.0)) then
        ifexcited = 0
      else
C         excitation approximation works well here
        ifexcited = 2
        nmax = 300 001
      endif
C       usually iftc = 0 (entropy(2) done with straight derivatives.  there
C       is some evidence that where entropy derivative is done with 
C       thermodynamic relations, significance loss is more severe (especially
C       in calculation of Q = rtp), but the difference is usually at the 1.d-13
C       level or better.
      iftc = 0
C       in all but one case,_use_reduced mass of electron in equilibrium
C         constant equation.
      ifreducedmass = 1
      if(ifoption.eq.1) then
C         pteh style
C         default ifcouloumb 
        if(tl.ge.tllim) then
          ifcoulomb = -15
        else
          ifcoulomb = 5
        endif
        if(ifmodified.eq.0) then
C           original version for interior calculations.
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
        elseif(ifmodified.eq.1) then
C           modified version with enhancements for interior calculations.
          ifpi = 1
        elseif(ifmodified.eq.2) then
C           modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 1
          ifrad = 0
C          non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
            ifexcited = 2
            nmax = 4
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.102) then
C           modified version to fit old opal table
          ifh2plus = 0
          morder = 1
C           same as 2 except for Debye-Huckel
          ifcoulomb = 1
          ifpi = 1
          ifrad = 0
C          non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
            ifexcited = 2
            nmax = 4
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.3) then
C           modified version to fit old opal table extension
          ifh2 = 0
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 1
          ifrad = 0
C          non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
            ifexcited = 2
            nmax = 4
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.4) then
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with no exchange
          ifexchange_in = 0
        elseif(ifmodified.eq.5) then
C          non-linear transform of G(eta) exchange
          ifexchange_in = -1
C           modified version with enhancements for interior calculations.
          ifpi = 1
        elseif(ifmodified.eq.6) then
C          non-linear, Kovetz et al, no strong degeneracy approximation.
C           G(eta) + relativistic Kapusta II equiv Kovetz terms
          ifexchange_in = 1
C           modified version with enhancements for interior calculations.
          ifpi = 1
        elseif(ifmodified.eq.7) then
C          non-linear, Kovetz et al, possible strong degeneracy approximation.
          ifexchange_in = 2
C           modified version with enhancements for interior calculations.
          ifpi = 1
        elseif(ifmodified.eq.8) then
C           non-linear, Kapusta, no strong degeneracy approximation
          ifexchange_in = 11
C           modified version with enhancements for interior calculations.
          ifpi = 1
        elseif(ifmodified.eq.9) then
C          linear, Kapusta, possible strong degeneracy approximation.
          ifexchange_in = 112
C           modified version with enhancements for interior calculations.
          ifpi = 1
        elseif(ifmodified.eq.11) then
C           modified version with enhancements for interior calculations.
C           but *always* used PTEH Coulomb sum approximation
          ifcoulomb = -15
          ifpi = 1
        elseif(ifmodified.eq.12) then
C           modified version with enhancements for interior calculations.
C           but_use_case b (metal Coulomb sum approximation)
          ifcoulomb = -5
          ifpi = 1
        elseif(ifmodified.eq.13) then
C           modified version with enhancements for interior calculations.
C           but_use_no Coulomb sum approximation (recommended for work
C           of high precision)
          ifcoulomb = 5
          ifpi = 1
        elseif(ifmodified.eq.15) then
C           Cody-Thacher
C           modified version with enhancements for interior calculations.
          morder = 1
          ifpi = 1
        elseif(ifmodified.eq.16) then
C           non-relativistic
C           modified version with enhancements for interior calculations.
          morder = -3
          ifpi = 1
        elseif(ifmodified.eq.17) then
C           non-relativistic
C           modified version with enhancements for interior calculations.
          morder = -5
          ifpi = 1
        elseif(ifmodified.eq.18) then
C           non-relativistic
C           modified version with enhancements for interior calculations.
          morder = -8
          ifpi = 1
        elseif(ifmodified.eq.20) then
C           drop h2 and h2+
C           modified version with enhancements for interior calculations.
          ifh2 = 0
          ifh2plus = 0
          ifpi = 1
        elseif(ifmodified.eq.21) then
C           qh2
C           modified version with enhancements for interior calculations.
          ifh2 = 4
          ifpi = 1
        elseif(ifmodified.eq.22) then
C           drop h2+
C           modified version with enhancements for interior calculations.
          ifh2plus = 0
          ifpi = 1
        elseif(ifmodified.eq.24) then
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with no Coulomb
          ifcoulomb = 0
        elseif(ifmodified.eq.25) then
C           no excitation or Planck-Larkin
C           modified version with enhancements for interior calculations.
          ifpi = -1
          ifexcited = 0
        elseif(ifmodified.eq.26) then
C           no excitation with Planck-Larkin
C           modified version with enhancements for interior calculations.
          ifpi = 1
          ifexcited = 0
        elseif(ifmodified.eq.27) then
C           Debye-Hueckel
C           modified version with enhancements for interior calculations.
          ifcoulomb = 1
          ifpi = 1
        elseif(ifmodified.eq.28) then
C           Debye-Hueckel with tau correction
C           modified version with enhancements for interior calculations.
          ifcoulomb = 2
          ifpi = 1
        elseif(ifmodified.eq.29) then
C           PTEH original Coulomb
C           modified version with enhancements for interior calculations.
          ifcoulomb = 3
          ifpi = 1
        elseif(ifmodified.eq.30) then
C           thermodynamic consistency
          iftc = 0
C           modified version with enhancements for interior calculations.
          ifpi = 1
        elseif(ifmodified.eq.31) then
C           thermodynamic consistency
          iftc = 1
C           modified version with enhancements for interior calculations.
          ifpi = 1
        elseif(ifmodified.eq.40) then
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with "exact" excitation calculation for H, He, and molecules
          if(lw.gt.2) then
            ifexcited = 12
            nmax = 10
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.41) then
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with approximate excitation calculation for H, He only
C           with default "infinite" nmax
          if(lw.gt.2) then
            ifexcited = 1
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.42) then
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with approximation excitation calculation for H, He
C           + molecules + metals.
C           with default "infinite" nmax
          if(lw.gt.2) then
            ifexcited = 3
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.43) then
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with approximate excitation calculation for H, He, and molecules
C           with special nmax
          if(lw.gt.2) then
            ifexcited = 2
            nmax = 4
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.44) then
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with approximate excitation calculation for H, He, and molecules
C           with special nmax
          if(lw.gt.2) then
            ifexcited = 2
            nmax = 10
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.45) then
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with approximate excitation calculation for H, He, and molecules
C           with special nmax
          if(lw.gt.2) then
            ifexcited = 2
            nmax = 50
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.101) then
C           ****EOS4****
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           *always* used PTEH Coulomb sum approximation
C           (causes substantial errors for LMS, but much quicker)
          ifcoulomb = -15
C           linearly transformed relativistic exchange
C           (causes substantial errors for LMS, but quicker)
C           n.b. found didn't matter that much
C          linear, Kapusta, possible strong degeneracy approximation.
!          ifexchange_in = 112
C          _use_lower-order approximation for fermi-dirac integrals?
C           The maximum errors are similar to the higher-order approximation,
C           (1.d-3 in ln P), but the larger errors are more widespread than
C           the morder = 5 case.  The fermi-dirac overhead doesn't matter
C           much for other forms of the EOS, but for the quick form it might
C           be useful to eliminate this extra source of overhead.
C           n.b. turned out to matter very little
!          morder = 3
        elseif(ifmodified.eq.103) then
C           emulation of Sweigart EOS
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with no Coulomb and no exchange
          ifcoulomb = 0
          ifexchange_in = 0
        elseif(ifmodified.eq.104) then
C           emulation of Sweigart EOS
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with no Coulomb and no exchange
          ifcoulomb = 0
          ifexchange_in = 0
C           and no radiation!
          ifrad = 0
        elseif(ifmodified.eq.105) then
C           emulation of Sweigart EOS
C           modified version with enhancements for interior calculations.
          ifpi = 1
C           but with no Coulomb and no exchange
          ifcoulomb = 0
          ifexchange_in = 0
C           and Cody-Thacher F-D integrals
          morder = 1
        elseif(ifmodified.eq.110) then
C           EOS3 with no radiation pressure.
          ifpi = 1
          ifrad = 0
        elseif(ifmodified.eq.111) then
C           EOS4 with no radiation pressure.
          ifpi = 1
          ifrad = 0
          ifcoulomb = -15
        elseif(ifmodified.eq.-1) then
C           original version for comparison with pteh table.
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
          ifrad = 0
        elseif(ifmodified.eq.-101) then
C           original pteh with radiation pressure added (by default)
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
        elseif(ifmodified.eq.-6) then
C           modified version with enhancements for interior calculations.
C           ifmodified = 1 defaults
C          non-linear, Kapusta, possible strong degeneracy approximation.
          ifexchange_in = 12
C           but with original (via ifmodified, ifpi) pteh 
C           pressure ionization without planck-larkin occupation probability
C           or excitation
          ifexcited = 0
          ifpi = -1
        elseif(ifmodified.eq.-7) then
C           modified version with enhancements for interior calculations.
C           ifmodified = 1 defaults
C          non-linear, Kapusta, possible strong degeneracy approximation.
          ifexchange_in = 12
C           but with original (via ifmodified, ifpi) geff
C           pressure ionization without planck-larkin occupation probability
C           or excitation
          ifexcited = 0
          ifpi = -2
        elseif(ifmodified.eq.-8) then
C           modified version with enhancements for interior calculations.
C           ifmodified = 1 defaults
C          non-linear, Kapusta, possible strong degeneracy approximation.
          ifexchange_in = 12
C           but with original (via ifmodified, ifpi) mhd
C           pressure ionization without planck-larkin occupation probability
C           or excitation
          ifexcited = 0
          ifpi = -3
        elseif(ifmodified.eq.-30) then
C           thermodynamic consistency
          iftc = 0
C           original version for interior calculations.
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
        elseif(ifmodified.eq.-31) then
C           thermodynamic consistency
          iftc = 1
C           original version for interior calculations.
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
        else
          stop 'free_eos: bad ifmodified for ifoption=1'
        endif
      elseif(ifoption.eq.2) then
C         fjs style
        if(ifmodified.eq.0) then
C           original version for interior calculations
          ifreducedmass = 0
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
        elseif(ifmodified.eq.1) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = -5
          ifpi = 2
        elseif(ifmodified.eq.2) then
C           modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ifcoulomb = -5
          ifpi = 2
          ifrad = 0
C          non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
            ifexcited = 2
            nmax = 4
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.3) then
C           modified version to fit old opal table extension
          ifh2 = 0
          ifh2plus = 0
          morder = 1
          ifcoulomb = -5
          ifpi = 2
          ifrad = 0
C           non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
            ifexcited = 1
            nmax = 4
          else
            ifexcited = 0
          endif
        elseif(ifmodified.eq.4) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = -5
          ifpi = 2
        elseif(ifmodified.eq.26) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = -5
          ifpi = 2
C           special with no excitation
          ifexcited = 0
        elseif(ifmodified.eq.-1) then
C           original version for comparison with geff code results
          ifreducedmass = 0
C           but with different H2 t limit (for strict mimicry)
          if(exp(tl).gt.1.d5) then
            ifh2 = 0
          endif
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
C           version of geff code I have ignores radiation pressure
C           so mimic this behaviour with the current code.
          ifrad = 0
        elseif(ifmodified.eq.-101) then
C           original geff with radiation pressure added
          ifreducedmass = 0
C           but with different H2 t limit (for strict mimicry)
          if(exp(tl).gt.1.d5) then
            ifh2 = 0
          endif
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
C           version of geff code I have ignores radiation pressure, but
C           put it in anyway (by default) for this option of the current code.
        elseif(ifmodified.eq.-2) then
C           original version for comparison with sireff code results
C           g(eta) exchange with linear transform
          ifexchange_in = -101
          ifreducedmass = 0
C           but with different H2 t limit (for strict mimicry)
          if(exp(tl).gt.1.d5) then
            ifh2 = 0
          endif
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
C           version of sireff code I have ignores radiation pressure
C           so mimic this behaviour with the current code.
          ifrad = 0
        elseif(ifmodified.eq.-201) then
C           Same as -2 above except for including radiation pressure.
C           g(eta) exchange with linear transform
          ifexchange_in = -101
          ifreducedmass = 0
C           but with different H2 t limit (for strict mimicry)
          if(exp(tl).gt.1.d5) then
            ifh2 = 0
          endif
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
C           version of sireff code I have ignores radiation pressure
C           but put it in anyway (by default) for this option of the current code.
        elseif(ifmodified.eq.-30) then
C           thermodynamic consistency
          iftc = 0
C           original version for interior calculations
          ifreducedmass = 0
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
        elseif(ifmodified.eq.-31) then
C           thermodynamic consistency
          iftc = 1
C           original version for interior calculations
          ifreducedmass = 0
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
        elseif(ifmodified.eq.30) then
C           thermodynamic consistency
          iftc = 0
C           modified version with enhancements for interior calculations.
          ifcoulomb = -5
          ifpi = 2
        elseif(ifmodified.eq.31) then
C           thermodynamic consistency
          iftc = 1
C           modified version with enhancements for interior calculations.
          ifcoulomb = -5
          ifpi = 2
        else
          stop 'free_eos: bad ifmodified for ifoption=2'
        endif
      elseif(ifoption.eq.3) then
C         mdh style
        if(lw.gt.2.and.ifmodified.gt.0) then
C          n.b._use_truncated sum because approximations don't work
C          very well for mhd case.
C          n.b. _use_excited states for H, He, and molecules but
C          not metals (for now).
          ifexcited = 12
          nmax = 10
        elseif(lw.gt.2.and.ifmodified.le.0) then
C           n.b._use_truncated sum because approximations don't work very well for mhd case
C           original mode doesn't have molecular electronic excitation
          ifexcited = 11
          nmax = 10
        else
          ifexcited = 0
        endif
        if(ifmodified.eq.0) then
C           original version for interior calculations.
C           we know morder=1 (Cody-Thacher Fermi-Dirac without relativistic
C           correction) is a bad assumption for red giant cores, but the
C           point is to demonstrate this for the original mhd mode.
          morder = 1
          ifcoulomb = 2
          ifpi = -3
        elseif(ifmodified.eq.1) then
C           EOS1 modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
        elseif(ifmodified.eq.11) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
C           for special table without radiation pressure
          ifrad = 0
        elseif(ifmodified.eq.2) then
C           modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 3
          ifrad = 0
C           non-linear transform of G(eta) exchange
          ifexchange_in = -1
C           follow opal table
          nmax = 4
        elseif(ifmodified.eq.102) then
C           modified version to fit old opal table
          ifh2plus = 0
          morder = 1
C           same as 2 except_use_Debye-Huckel
          ifcoulomb = 1
          ifpi = 3
          ifrad = 0
C           non-linear transform of G(eta) exchange
          ifexchange_in = -1
C           follow opal table
          nmax = 4
        elseif(ifmodified.eq.103) then
C           modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 3
C          _use_default ifrad behaviour (only change from 2)
!          ifrad = 0
C           non-linear transform of G(eta) exchange
          ifexchange_in = -1
C           follow opal table
          nmax = 4
        elseif(ifmodified.eq.20) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
C           drop h2 and h2+
          ifh2 = 0
          ifh2plus = 0
        elseif(ifmodified.ge.201.and.ifmodified.le.231) then
C          Same as ifmodified = 1 except for ifred = 0 to match EOS2005 table
          ifcoulomb = 5
          ifpi = 3
          ifrad = 0
C          only difference for 201 <= ifmodified < 210
          if(ifmodified.eq.201) then
C            non-linear, Kapusta, possible strong degeneracy approximation.
            ifexchange_in = 12
          elseif(ifmodified.eq.202) then
C            alternative Coulomb
            ifcoulomb = 6
C            non-linear, Kapusta, possible strong degeneracy approximation.
            ifexchange_in = 12
          elseif(ifmodified.eq.203) then
C            no exchange
            ifexchange_in = 0
          elseif(ifmodified.eq.204) then
C            no Coulomb
            ifcoulomb = 0
C            non-linear, Kapusta, possible strong degeneracy approximation.
            ifexchange_in = 12
          elseif(ifmodified.eq.205) then
C            no Coulomb or exchange
            ifcoulomb = 0
            ifexchange_in = 0
          elseif(ifmodified.gt.210) then
C            Standard Coulomb and exchange treatment
C            but change Fermi-Dirac morder through all variations.
            if(ifmodified.eq.211) then
              morder = 1
            elseif(ifmodified.eq.212) then
              morder = -3
            elseif(ifmodified.eq.213) then
              morder = -5
            elseif(ifmodified.eq.214) then
              morder = -8
            elseif(ifmodified.eq.215) then
              morder = -21
            elseif(ifmodified.eq.216) then
              morder = -23
            elseif(ifmodified.eq.217) then
              morder = 3
            elseif(ifmodified.eq.218) then
              morder = 5
            elseif(ifmodified.eq.219) then
              morder = 8
            elseif(ifmodified.eq.220) then
              morder = 13
            elseif(ifmodified.eq.221) then
              morder = 15
            elseif(ifmodified.eq.222) then
              morder = 18
            elseif(ifmodified.eq.223) then
              morder = 21
            elseif(ifmodified.eq.224) then
              morder = 23
            endif
          endif
        elseif(301.le.ifmodified.and.ifmodified.le.320) then
C          This series to be compared with ifmodified.eq.1 (ifcoulomb=5, ifpi=3)
          ifcoulomb = 5
          ifpi = 3
          if(ifmodified.eq.301) then
C            nocoulomb
            ifcoulomb = 0
          elseif(ifmodified.eq.302) then
C            dh
            ifcoulomb = 1
          elseif(ifmodified.eq.303) then
C            dh tau
            ifcoulomb = 2
          elseif(ifmodified.eq.304) then
C            g_pteh(gamma) with everything else the same as ifcoulomb=5
            ifcoulomb = 9
          elseif(ifmodified.eq.311) then
C            no exchange
            ifexchange_in = 0
          elseif(ifmodified.eq.312) then
C            normal exchange (14) with linear inversion approximation
            ifexchange_in = 114
          elseif(ifmodified.eq.313) then
C            normal exchange (14) with non-relativistic approximation
            ifexchange_in = -14
          elseif(ifmodified.eq.314) then
C            series approximation(s) for exchange
            ifexchange_in = 12
          elseif(ifmodified.eq.315) then
C            higher-order approximation for exchange
            ifexchange_in = 15
          elseif(ifmodified.eq.316) then
C            highest-order approximation for exchange
            ifexchange_in = 16
          endif
        elseif(ifmodified.eq.3) then
C           modified version to fit old opal table extension
          ifh2 = 0
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 3
          ifrad = 0
C           non-linear transform of G(eta) exchange
          ifexchange_in = -1
C           follow opal table
          nmax = 4
        elseif(ifmodified.eq.4) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
C           but with no exchange
          ifexchange_in = 0
        elseif(ifmodified.eq.24) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
C           but with no Coulomb
          ifcoulomb = 0
        elseif(ifmodified.eq.26) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
C           but with no excitation
          ifexcited = 0
        elseif(ifmodified.eq.30.or.ifmodified.eq.1030) then
C           thermodynamic consistency
          iftc = 0
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          if(ifmodified.eq.1030) morder = 23
        elseif(ifmodified.eq.31.or.ifmodified.eq.1031) then
C           thermodynamic consistency
          iftc = 1
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          if(ifmodified.eq.1031) morder = 23
        elseif(ifmodified.eq.43) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
C           but with special nmax
          nmax = 4
        elseif(ifmodified.eq.44) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
C           but with special nmax
          nmax = 25
        elseif(ifmodified.eq.45) then
C           modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
C           but with special nmax
          nmax = 50
        elseif(ifmodified.eq.-1) then
C           original version for interior calculations.
C           we know morder=1 (Cody-Thacher Fermi-Dirac without relativistic
C           correction) is a bad assumption for red giant cores, but the
C           point is to demonstrate this for the original mhd mode.
          morder = 1
          ifcoulomb = 2
          ifpi = -3
C           but_use_no H2+
          ifh2plus = 0
        elseif(ifmodified.eq.-30) then
C           thermodynamic consistency
          iftc = 0
C           original version for interior calculations.
          morder = 1
          ifcoulomb = 2
          ifpi = -3
        elseif(ifmodified.eq.-31) then
C           thermodynamic consistency
          iftc = 1
C           original version for interior calculations.
          morder = 1
          ifcoulomb = 2
          ifpi = -3
        else
          stop 'free_eos: bad ifmodified for ifoption=3'
        endif
      elseif(ifoption.eq.4) then
        if(ifmodified.eq.2) then
C           modified version to fit opal table for Saumon style
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 4
          ifrad = 0
C           non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(tl.lt.tllim) then
            ifexcited = 12
            nmax = 4
          else
            lw = 2
            ifexcited = 0
          endif
        elseif(ifmodified.eq.3) then
C           modified version to fit He Saumon tables
          ifh2 = 0
          ifh2plus = 0
          morder = 1
          ifexchange_in = 0
C           He table calculated with SOCP for log10 rho > 0.5, otherwise DH
C           theory.  For He keep comparison less than log10 rho = 0.5 and use
C           DH theory.
          ifcoulomb = 1
          ifpi = -4
          ifrad = 0
          ifexcited = 0
          if(tl.ge.tllim) then
            lw = 2
          endif
        elseif(ifmodified.eq.4
     &      .or.ifmodified.eq.30.or.ifmodified.eq.31) then
C           modified version to fit H Saumon tables
          ifh2plus = 0
          morder = 1
C          linear, Kapusta, possible strong degeneracy approximation.
          ifexchange_in = 112
          ifcoulomb = 5
          ifpi = -4
          ifrad = 0
          if(tl.ge.tllim) then
            lw = 2
            ifexcited = 0
          else
            ifexcited = 12
            nmax = 100
          endif
          if(ifmodified.eq.31) iftc = 1
        else
          stop 'free_eos: bad ifmodified for ifoption=4'
        endif
      elseif(ifoption.eq.5) then
C         ifpi = 5 (Planck-Larkin pressure ionization, only) options
        if(ifmodified.eq.1) then
C           same as standard (ifoption.eq.1.and.ifmodified.eq.1) except
C          _use_Planck-Larkin only.
          if(tl.ge.tllim) then
            ifcoulomb = -15
C             all elements approximated as fully ionized.
            lw = 2
            ifmtrace = 1
            ifexcited = 0
            ifpi = 0
          else
            ifcoulomb = 5
            ifpi = 5
          endif
        elseif(ifmodified.eq.2) then
          morder = 1
          ifcoulomb = 5
          ifpi = 5
          ifrad = 0
        else
          stop 'free_eos: bad ifmodified for ifoption=5'
        endif
      else
        stop 'free_eos: bad ifoption'
      endif
      if(tl.ge.tllim) then
C         turn off molecules for T > 1.d6
C         the partition functions are approximated above T = 1.d5
C         by a Taylor series approach which assures continuity for
C         second-order thermodynamic quanties (e.g. grada).
C         This approximation is well behaved (at least), and the
C         errors in it don't matter very much because n(H2) and n(H2+)
C         are so small for these temperatures (for the densities 
C         associated with stars).  Eventually the Taylor series would
C         underflow (negative curvature parabola in ln Q vs ln T in
C         all cases) so that is why we turn it off above t=10^6.
C         (also code is more efficient).  We have also tried turning
C         off above T = 10^5, but leaves tiny but noticable discontinuity.
        ifh2 = 0
        ifh2plus = 0
      endif
      call free_eos_detailed(
     &  ifh2, ifh2plus, morder, ifexchange_in, 
     &  ifmtrace, ifcoulomb, ifpi, ifrad, lw, ifexcited, nmax,
     &  ifreducedmass, iftc, ifmodified, kif,
     &  eps, neps, match_variable, fl, tl, fm, ft,
     &  t, rho, rl, p, pl, cf, cp, sf, st, grada, rtp,
     &  rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,
     &  gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e,
     &  degeneracy, pressure, density, energy, enthalpy, entropy,
     &  iteration_count, info)
C      disable qp and qf (see comment above)
      qp = 1.d300
      qf = 1.d300
C      signal calling routine if abnormal ending for free_eos_detailed.
      if(info.gt.0) iteration_count = -info
      fl_old = fl
C       special call to test thermodynamic consistency
      if(abs(ifmodified).eq.30.or.abs(ifmodified).eq.31.or.
     &    abs(ifmodified).eq.1030.or.abs(ifmodified).eq.1031) then
        if(entropy(2).eq.0.d0) then
          pl = 0.d0
        else
C           this coefficient is normally negative
          pl = log(abs(entropy(2)))
        endif
      endif
      end
