C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: free_eos_detailed.f 841 2008-07-07 19:12:20Z airwin $
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
      subroutine free_eos_detailed(
     &  ifh2, ifh2plus, morder, ifexchange_in, 
     &  ifmtrace, ifcoulomb, ifpi, ifrad, lw, ifexcited, nmax,
     &  ifreducedmass, iftc, ifmodified, kif,
     &  eps, neps, match_variable, fl, tl, fm, ft,
     &  t, rho, rlout, p, pl, cf, cp, sf, st, grada, rtp,
     &  rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,
     &  gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e,
     &  degeneracy, pressure, density, energy, enthalpy, entropy,
     &  iteration_count, info)
C       input quantities:
C       ifh2 = 0  no h2
C       ifh2 = 1  vdb h2
C       ifh2 = 2  st h2
C       ifh2 = 3  irwin h2
C       ifh2plus = 0  no h2plus
C       ifh2plus = 1  st h2plus
C       ifh2plus = 2  irwin h2plus
C       morder_in = 3, 5, or 8  (or 13, 15, or 18) means 3rd, 5th, or 8th order
C       fermi_dirac integral approximation following EFF fit (or modified
C       version of EFF fit which reduces to Cody-Thacher approximation for
C       low relativistic correction).
C       morder_in = -3, -5, or -8 (or -13, -15, or -18) means_use_above
C       approximations in non-relativistic limit.
C       morder_in = 1 uses Cody-Thacher approximation directly.
C       morder_in = 21 calculate Fermi-Dirac integrals with slow, but precise
C       (~1.d-9 relative errors) numerical integration.
C       morder_in = -21 is same as 21 in non-relativistic limit.
C       morder_in = 23 means original 3rd order eff result
C       morder_in = -23 is same as 23 in non-relativistic limit.
C      if abs(ifexchange_in) > 100 then_use_linear transform approximation
C        for exchange treatement.  Otherwise,_use_numerical transform as
C        described in research note.
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
      
C       n.b. important metals controlled by array iftracemetal below.
C         currently list includes C, N, O, Ne, Mg, Si, S, Fe
C       ifmtrace = 0 treat important metals like H and He.
C       ifmtrace = 1 treat everything but H, He like trace metals; e.g.,
C         partially ionized (lw.gt.3) or fully ionized (lw.le.3).
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
C       lw > 3 very slow option with all elements treated as partially ionized
C         using all stages of ionization.
C       ifexcited > 0 means_use_excited states (must have Planck-Larkin or
C         |ifpi| = 3 or 4).
C          0 < ifexcited < 10 means_use_approximation to explicit summation
C         10 < ifexcited < 20 means_use_explicit summation
C         mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules)
C           and helium.
C         mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
C         mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
C       nmax is maximum principal quantum number included in excited state
C         sum for explicit summation.  Also, same interpretation when approximate
C         summation used for pure Planck-Larkin case with special signal that
C         nmax > 300 000 means_use_approximate for infinite Planck-Larkin sum
C         rather than approximation for finite Planck-Larkin sum up to nmax.
C         However, if approximate and mhd (Planck-Larkin or not),
C         then nmax not programmed.
C       ifreducedmass = 1 (use reduced mass in equilibrium constant)
C       ifreducedmass = 0 (use electron mass in equilibrium constant)
C       iftc is only currently meaningful with kif = 1 or 2.
C       iftc = 1 (use thermodynamic consistency arguments to obtain appropriate
C         derivatives of entropy).
C       iftc = 0 (use direct analytical derivatives of entropy.))
      
C       ifmodified <= 0 (original form of pressure ionization for |ifpi| = 1,2,3,4
C       ifmodified > 0 (modified form of pressure ionization)
C       kif = 0, ln f and ln t are independent variables.
C       kif = 1, ln p = match_variable and ln t are independent variables.
C       kif = 2, ln rho = match_variable and ln t are independent variables.
C       n.b. kif > 0 are much slower options because ln f must be iterated to
C         match match_variable
C       eps(nelements_in) is an array of neps = nelements = 20 values of relative
C         abundance by weight divided by the appropriate atomic weight.
C         We_use_the atomic weight scale where (un-ionized) C(12)
C         has a weight of 12.00000000....  All weights are for the 
C         un-ionized element. The eps value for an element
C         should be the sum of the individual isotopic eps values for
C         that element.  The eps array refers to the
C         elements in the following order:
C         H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni
C       match_variable is matched by iterative adjustment of fl when kif > 0
C       fl = ln of EFF degeneracy parameter related to eta by equation below.
C         fl is starting value for iterative adjustment upon input
C         and iteratively adjusted for output when kif > 0
C       tl = ln t
C       output quantities:
C       fm = partial fl(tl, match_variable)/ partial match_variable (kif > 0).
C       ft = partial fl(tl, match_variable)/ partial tl  (kif > 0).
C       n.b. when ifrad = 2, kif = 1, these derivatives are with respect to
C       ln P_gas.  All other output derivatives for this combination of flags
C       are taken w.r.t. ln pressure.
C       t = temperature
C       rho = density
C       rlout = ln rho
C       p = pressure (total if ifrad = 1 or 2, excluding radiation if ifrad = 0)
C       pl = ln p
C       cf = cp (- partial ln rho(T,P)/partial T)^{1/2}
C       cp = specific heat at constant pressure
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
C       iteration_count = the total number of ionization fraction loops completed
C         to get the solution
C       info = status code.  IMPORTANT: anything other than zero means
C         there was an abnormal ending to the FreeEOS calculation and no
C         returned quantities should be considered to be reliable.
C       degeneracy, pressure, density, energy, enthalpy, and entropy
C         are all 3-vectors, with the second component being the derivative
C         of the first component wrt match_variable (except for the case where
C         ifrad = 2), and the third component being the derivative of the first
C         component wrt tl.
C       definitions:
C         degeneracy(1) is EFF degeneracy parameter ln f defined above.
C         pressure(1) = ln pressure (except for kif=2, ifrad=2).
C         density(1) = ln density.
C         energy(1) = internal energy per unit mass.
C         enthalpy(1) = enthalpy per unit mass = energy(1) + p/rho.
C         entropy(1) = entropy per unit mass.
C       N.B. qf, qp used to be returned to help calculate the
C         approximate gravitational energy generation rate for the case
C         when abundances are changing, but there is a much better way
C         to calculate the exact gravitational energy generation rate
C         (see comment in free_eos.f) so qf and qp have been dropped from
C         the free_eos_detailed argument list.
      implicit none
      include 'aux_scale.h'
      include 'constants.h'
      include 'ionization.h'
      double precision pr_ratio, ppie_ratio, pc_ratio, pex_ratio
      common /diagnostics/pr_ratio, ppie_ratio, pc_ratio, pex_ratio
      integer kif, ifreducedmass, iftc, ifmodified, lw, ifexcited,
     &  nmax, neps, neps_local, ifpi_fit, ifpi_fit_old
      data ifpi_fit_old/-1000/
      parameter(neps_local = 20)
      double precision match_variable, match_variablef, match_variablet,
     &  fm, ft
      double precision t, tl,
     &  rho, rlout, p, pl, cf, cp, qf, qp, sf, st, grada,
     &  fl, rmue, fh2, fhe2, fhe3, xmu1, xmu3,
     &  eps(neps), eps_old(neps_local), eta
C      required for all debugging tests.
C      dvdaux refers to the partial of dv wrt old auxiliary variables
C      and fl and tl.
C      dauxdv refers to the partial of the new auxiliary variables wrt dv
C      and fl and tl.
C      jacobian refers to the partial of the zeroed functions (either
C      logarithmic form of aux(old) - aux(new) or calculated match_variable
C      - input match_variable) wrt old auxiliary variables and fl (and ft).
C      n.b. MUST be zero to start to force ifsame_abundances and
C      ifsame_zero_abundances to be 0 on first entry into routine.
      data eps_old/20*0.d0/
      double precision degeneracy(3), pressure(3), density(3), energy(3), 
     &  enthalpy(3), entropy(3)
      double precision flminus, flplus, paaminus, paaplus
      double precision rhostar(9), pstar(9), sstar(3), ustar(3), 
     &  dni(2), dne(2), ne, ni, nef, net
      double precision re, ref, ret, pe, pef, pet, se, sef, set, ue, uef, uet,
     &  pe_cgs
      equivalence
     &  (rhostar(1),re), (rhostar(2),ref), (rhostar(3),ret),
     &  (pstar(1),pe), (pstar(2),pef), (pstar(3),pet),
     &  (sstar(1),se), (sstar(2),sef), (sstar(3),set),
     &  (ustar(1),ue), (ustar(2),uef), (ustar(3),uet)
      integer ifprintrtp, ifprintcp, ifprintgrada
      data ifprintrtp/1/  !one warning on negative rtp
      data ifprintcp/1/  !one warning on negative cp
      data ifprintgrada/1/  !one warning on negative grada
      data dni/0.75d0, 0.166667d0/
      data dne/0.5d0, 0.0d0/
      integer ltermax
C      maximum number of allowed ln f iterations
      parameter (ltermax = 200)
      integer nelements
C      these two parameters must agree exactly with corresponding parameters
C      in the ionize subroutine.
      parameter (nelements = 20)
      integer ion_end(nelements+2), mion_end, number_electrons
      integer maxfjs_aux
      parameter (maxfjs_aux = 4)
      double precision dvzero(nelements), 
     &  dv(nions+2), dvf(nions+2), dvt(nions+2),
     &  dve_exchange, dve_exchangef, dve_exchanget,
     &  plop(nions+2), plopt(nions+2), plopt2(nions+2),
     &  dv_pl(nions+2), dv_plt(nions+2)
      integer maxcorecharge
C      maximum core charge for non-bare ion
      parameter(maxcorecharge = 28)
      integer nmin(maxcorecharge), nmin_max(maxcorecharge),
     &  iz, izlo, izhi, ion_start
      integer nmin_species(nions+2)
      integer iatomic_number(nelements)
      data iatomic_number/1,2,6,7,8,10,11,12,13,14,15,16,17,18,20,
     &  22,24,25,26,28/
      integer ifnr,
     &  ifh2, ifh2plus, morder, ifexchange_in, ifmtrace, ifcoulomb,
     &  if_pteh, ifcoulomb_mod, if_mc, if_dc, i, lter, iteration_count,
     &  info,
     &  ifionized, ioncount, ion, max_index, maxioncount,
     &  index, maxnextrasum,
     &  nextrasum, ifpi, ifpi_local, ifpl, ifrad, nxextrasum,
     &  ifsame_abundances, ifnear_abundances,
     &  ifsame_zero_abundances,
     &  iflast, lw_old, ifexcited_old, nmax_old,
     &  ifh2_old, ifh2plus_old, morder_old, ifmtrace_old, ifcoulomb_old,
     &  ifmtrace_olda, lw_olda, if_mc_olda, ifh2_olda, ifh2plus_olda
      parameter (maxioncount=2000)
      data lw_old, ifexcited_old, nmax_old,
     &  ifh2_old, ifh2plus_old, morder_old, ifmtrace_old, ifcoulomb_old,
     &  ifmtrace_olda, lw_olda, if_mc_olda, ifh2_olda, ifh2plus_olda
     &  /13*-1000/
C      mdh-like treatment has a maximum of 9 extra parameters
      parameter (maxnextrasum = 9)
      parameter (nxextrasum = 4)
      double precision gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, 
     &  gamma_e,
     &  full_sum0, full_sum1, full_sum2,
     &  sum0, sum0a, sumion0, sum0_mc, sum0ne, sum0f, sum0t,
     &  sum2, sum2a, sumion2, sum2_mc, sum2ne, sum2f, sum2t,
     &  r_ion3(nions+2), r_neutral(nelements+2),
     &  charge(nions+2), charge2(nions+2),
     &  extrasum(maxnextrasum), extrasumf(maxnextrasum), 
     &  extrasumt(maxnextrasum),
     &  sumpl1, sumpl1f, sumpl1t, sumpl2,
     &  pr, tc2
      double precision paa, pab, pac_nr, pac, lnrho_5, f, wf, n_e,
     &  hcon_mc, hecon_mc,
     &  dpcoulomb, dpcoulombt, dpcoulombf,
     &  dscoulomb, dscoulombf, dscoulombt, ducoulomb,
     &  ne_old, lambda_old, s,
     &  sion, sionf, siont, u, uion,
     &  h, hf, ht, hd, hdf, hdt, he, hef, het, rt, rf,
     &  h2, h2f, h2t, h2_dv(nions+2),
     &  h2plus, h2plusf, h2plust, h2plus_dv(nions+2),
     &  p0, pion, pionf, piont, pf, sr, pt, rpt, rtp,
     &  vx, vz, gamma, 
     &  pex, pexf, pext, sex, sexf, sext, uex,
     &  ppi, ppif, ppit, spi, spif, spit, upi,
     &  pexcited, pexcitedf, pexcitedt,
     &  sexcited, sexcitedf, sexcitedt, uexcited,
     &  nux, nuy, nuz,
     &  pnorad, pnoradf, pnoradt,
     &  match_variableold, tlold, flold
C      must be wild values.
      data match_variableold, tlold, flold/3*1.d30/
      integer iftracemetal(nelements)
C      iftracemetal controls what elements are treated as fully ionized
C      when lw = 3.
C      H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni
      data iftracemetal/
C      MDH dela residual ridge < 0.02
!      1 0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
C      MDH dela residual ridge < 0.01 if Fe not trace
!      1 0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1/
C      MDH dela residual ridge < 0.001 if Mg, Si, and S not trace
     &  0,0,0,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,0,1/
      integer ifdv(nions+2), ion_stop, ifelement(nelements),
     &  n_partial_ions, partial_ions(nions+2), inv_ion(nions+2),
     &  n_partial_elements, partial_elements(nelements+2)
C      set up equivalenced auxiliary variables that must be iterated to 
C      consistency.
      double precision rl, rl_old
C      iextraoff is the number of non-extrasum auxiliary variables
C      naux is the total number of auxiliary variables.
      integer iextraoff, naux
      parameter (iextraoff = 7)
C      one extra in case degeneracy parameter is included as one
C      of the auxiliary variables.
      parameter (naux = iextraoff + maxnextrasum + nxextrasum + 1)
      integer ifaux(naux), n_partial_aux, njacobian,
     &  njacobian_final, partial_aux(naux), inv_aux(naux),
     &  index_aux, iaux, jndex_aux, jaux,
     &  ielement, index_max,
     &  ioncountzero, ifzerocount
      double precision aux(naux), auxf(naux), auxt(naux),
     &  aux_old(naux), aux_restore(naux), fl_restore,
     &  faux(naux), faux_nr(naux), 
     &  zerolim(naux),
     &  eps_aux, faux_limit, faux_limit_small_start,
     &  faux_limit_small, faux_limit_large, faux_scale,
     &  maxfaux, maxfaux_diag, maxfaux_diag_old,
     &  daux_dv(nions+2,naux), ddv_aux(nions+2,naux),
     &  h_dv(nions+2), hd_dv(nions+2), he_dv(nions+2),
     &  r_dv(nions+2), sum0_dv(nions+2), sum2_dv(nions+2),
     &  extrasum_dv(nions+2, maxnextrasum),
     &  xextrasum(nxextrasum), xextrasumf(nxextrasum),
     &  xextrasumt(nxextrasum), xextrasum_dv(nions+2,nxextrasum),
     &  rhs(naux), rhs_save(naux),
     &  jacobian(naux,naux), jacobian_save(naux, naux),
     &  pnorad_daux(naux+1),
     &  bmin(maxcorecharge)
      logical allow_log(naux), allow_log_neg(naux),
     &  ifnoaux_iteration, ifwarm, any_ifsimple,
     &  ifsame_under, ifdvzero, ifnuform
      double precision simple_lambda_max, simple_lambda_min,
     &  simple_lambda_ratio, simple_lambda(naux)
C      Set up fast attack (immediate raise of simple lambda to
C      simple_lambda_max if simple criteria are satisfied) and slow
C      decay (multiply simple_lambda by simple_lambda_ratio whenever the
C      criteria are not met).
      parameter(simple_lambda_max = 1.d3)
C      relatively slow decay (6 iterations to zero simple_lambda)
      parameter(simple_lambda_ratio = 1.d-1)
C      floor value below which simple_lambda is set to zero.
      parameter(simple_lambda_min = 1.1d-2)
C      lapack variables
      double precision lu_lapack(naux,naux),
     &  row_lapack(naux), col_lapack(naux), sol_lapack(naux,2),
     &  rcond_lapack, ferr_lapack(2), berr_lapack(2),
     &  work_lapack(4*naux), rhs_lapack(naux,2)
      integer ipiv_lapack(naux), iwork_lapack(naux), info_lapack
      character*1 equed_lapack
      parameter (faux_limit_small_start = 2.d0)
C      used 50 with dgeco, and 50 seems to be okay with lapack svdc.
C      14 of these steps would be close to the over or under flow limit.
      parameter (faux_limit_large = 50.d0)
C      new control logic:
C      limit used to help control when diag solution used.  That diag solution
C      only used if the (log) change in auxiliary variable is greater than
C      this value and diag solution smaller than NR solution.
      double precision faux_limit_diag
      parameter(faux_limit_diag = 0.5d0)
C      limit on smallest value of faux_limit_small (which is ordinarily
C      decreased when maxfaux changes sign)
      double precision faux_limit_min
C      minimum faux_limit
      parameter(faux_limit_min = 1.d-10)
C      this equivalence organization is *important* for later logic on
C      deciding ifaux.  The order is h, hd, he, rl, h2plus, sum0, sum2, 
C      extrasum(maxnextrasum), xextrasum(nxextrasum)
      equivalence
     &  (aux(1),h),(aux(2),hd),(aux(3),he),
     &  (auxf(1),hf),(auxf(2),hdf),(auxf(3),hef),
     &  (auxt(1),ht),(auxt(2),hdt),(auxt(3),het),
     &  (aux(4),rl),(aux(5),h2plus),
     &  (auxf(4),rf),(auxf(5),h2plusf),
     &  (auxt(4),rt),(auxt(5),h2plust),
     &  (aux(6),sum0),(aux(7),sum2),
     &  (auxf(6),sum0f),(auxf(7),sum2f),
     &  (auxt(6),sum0t),(auxt(7),sum2t),
     &  (aux(iextraoff+1),extrasum(1)),
     &  (auxf(iextraoff+1),extrasumf(1)),
     &  (auxt(iextraoff+1),extrasumt(1)),
     &  (aux(iextraoff+maxnextrasum+1),xextrasum(1)),
     &  (auxf(iextraoff+maxnextrasum+1),xextrasumf(1)),
     &  (auxt(iextraoff+maxnextrasum+1),xextrasumt(1))
C      equivalence daux_dv and h_dv, hd_dv, he_dv, r_dv, h2plus, extrasum_dv
C      note aux index of daux_dv cannot be indexed because these arrays
C      used separately inside eos_calc so that aux position must be fixed.
C      Note, however, that in all _dv arrays, the dv index *is* indexed.
      equivalence
     &  (daux_dv(1,1), h_dv),
     &  (daux_dv(1,2), hd_dv),
     &  (daux_dv(1,3), he_dv), 
     &  (daux_dv(1,4), r_dv),
     &  (daux_dv(1,5), h2plus_dv),
     &  (daux_dv(1,6), sum0_dv), 
     &  (daux_dv(1,7), sum2_dv), 
     &  (daux_dv(1,iextraoff+1), extrasum_dv),
     &  (daux_dv(1,iextraoff+maxnextrasum+1), xextrasum_dv)
      integer iffirst
      data iffirst/1/
      integer ifmodified_old
C      silly value
      data ifmodified_old/-1000/
      !temporary
C      variables associated with free energy calculation:
      double precision
     &  free_rad, free_e, free_ion, free_pi, free_excited,
     &  free_coulomb,
     &  free_ex, free_pl
      double precision free, free_fp, free_old, free_daux(naux+1)
      logical eos_bfgs_called
C      variables required for debugging tests.
      logical debug_dvdaux, debug_dauxdv,
     &  debug_jacobian, debug_bfgs, debug_any
C      normally, all of these are false, but set one (and only one) of
C      these to be true to do some debugging of the dvdaux, dauxdv, or
C      jacobian partial derivatives internally calculated by
C      free_eos_detailed.
      parameter (debug_dvdaux = .false.)
      parameter (debug_dauxdv = .false.)
      parameter (debug_jacobian = .false.)
      parameter (debug_bfgs = .false.)
      parameter (debug_any = debug_dvdaux.or.debug_dauxdv.or.
     &  debug_jacobian.or.debug_bfgs)
      integer itemp1, itemp2, iaux1, iaux2
      double precision match_variable_in, match_variable_save,
     &  tl_in, tl_save, delta1, delta2,
     &  x(nions+2+nelements)!, y, g(nions+2+nelements)
      integer debug_count
      data debug_count/0/
C      end of variables required for debugging tests.
C      variables required for determining ifmajor
      logical ifmajor(naux)
      integer ijacobian, jjacobian
      double precision row_norm, major_crit
      parameter(major_crit = 1.d-3)
C      end of variables required for determining ifmajor
      integer isimple_bfgs
      logical debug_output
      parameter(debug_output = .false.)
      save
C      default good status
      info = 0
C      approximation for excited states has derivative error and/or
C      significance loss for combination of mhd and pl occupation
C      probability.  Direct summation works quite well with nmax ~ 100
C      when mhd occupation probabilities are active so drop excited
C      approximation in all cases but pl-excitation (where it works
C      very well and saves enormous amounts of time).
      if(ifexcited.gt.0.and.ifexcited.lt.10.and.(
     &  abs(ifpi).eq.3.or.abs(ifpi).eq.4))
     &  stop 'free_eos_detailed: ifexcited approximation not debugged'
      if(ifexcited.gt.0.and..not.
     &  (abs(ifpi).eq.3.or.abs(ifpi).eq.4.or.ifpi.gt.0))
     &  stop 'free_eos_detailed: bad ifexcited and ifpi combination'
C      NR auxiliary variable iteration.  Go one more step than
C      this criterion so should guarantee near machine precision
C      for the expected (and often tested) quadratic convergence.
      eps_aux = 1.d-7
      if(neps.ne.nelements.or.neps.ne.neps_local)
     &  stop 'free_eos_detailed: bad neps value'
      if(debug_any) then
C        this logic dependent on free_eos_test logic of original call
C        plus 8 stanzas of +/- delta1 and +/- delta2 for a total
C        of 33 calls to free_eos (and free_eos_detailed) per grid point.
        if(mod(debug_count,33).eq.0) then
          match_variable_save = match_variable
          tl_save = tl
        endif
        debug_count = debug_count + 1
      endif
      if(iffirst.eq.1) then
        iffirst = 0
C        aux_scale initialization:
C        default is a scale factor of unity
        do iaux = 1, naux_scale
          aux_scale(iaux) = 1.d0
        enddo
        if(.true.) then
C        certain aux_scale factors are non-unity because of the equivalance
C        statements between aux_scale and the following *_scale variables:
        sum0_scale = 4.d200
        sum2_scale = 1.d200
        extrasum_scale(1) = 1.d200
        do iaux = 2, maxnextrasum_scale-2
          extrasum_scale(iaux) = 1.d8*extrasum_scale(iaux-1)
        enddo
        extrasum_scale(maxnextrasum_scale-1) = extrasum_scale(1)
        extrasum_scale(maxnextrasum_scale) = extrasum_scale(1)
        !else
        xextrasum_scale(1) = 1.01d200
        xextrasum_scale(2) = 1.02d200
        xextrasum_scale(3) = 1.03d200
        xextrasum_scale(4) = 1.04d200
        endif
        if(debug_output) then
          write(0,*) 'aux_scale ='
          write(0,'(5(1pd14.5,2x))') aux_scale
        endif

C        For Planck-Larkin and effective radius need ionization 
C        potential of H2+.  This is the
C        dissociation energy of H2+ (which dissociates to H + H+) + the
C        ionization potential of monatomic H (= bi(1)).  
C        But the dissociation energy
C        of H2+ = h2diss - h2ion + bi(1)
        bi(nions+2) = h2diss + 2.d0*bi(1) - bi(nions+1)
!        temporary change to get energy scale onto the opal system.
C        Forrest Rogers believes he used 1/k = 11605.4 (where k expressed
C        in ev/K) and he also used 1 Rydberg = 13.6058ev
C        put Helium on Rydberg scale
!        bi(3) = bi(3)/bi(1)
!        bi(2) = bi(2)/bi(1)
C        1 Rydberg in cm^-1 under opal energy system
!        bi(1) = 11605.4d0*13.6058d0/c2
!        bi(2) = bi(2)*bi(1)
!        bi(3) = bi(3)*bi(1)
        ion_end(1) = iatomic_number(1)
        do ielement = 2, nelements
          ion_end(ielement) = ion_end(ielement-1) +
     &      iatomic_number(ielement)
        enddo
        do ion = 1,nions+2
          charge(ion) = dble(nion(ion))
          charge2(ion) = charge(ion)*charge(ion)
        enddo
        ion_start = 1
        do ielement = 1, nelements
C          for each element, and hydrogen molecules calculate nmin
C          nmin is 2 for 1 - 2-electron systems, 
C          nmin is 3 for 3 - 10-electron systems, 
C          nmin is 4 for 11 - 28-electron systems
          do ion = ion_start, ion_end(ielement)
            number_electrons = ion_end(ielement) - ion + 1
            if(number_electrons.le.2) then
              nmin_species(ion) = 2
            elseif(number_electrons.le.10) then
              nmin_species(ion) = 3
            else
              nmin_species(ion) = 4
            endif
          enddo
          ion_start = ion_end(ielement) + 1
        enddo
C        H2 is a 2-electron system
        nmin_species(nions+1) = 2
C        H2+ is a 1-electron system
        nmin_species(nions+2) = 2
      endif
      if(debug_any) then
C        keep track of current and original values
        match_variable_in = match_variable
        match_variable = match_variable_save
        if(kif.eq.0) fl = match_variable_save
        delta1 = 1.d0*(match_variable_in - match_variable_save)
        tl_in = tl
        tl = tl_save
        delta2 = 1.d0*(tl_in - tl_save)
      endif
C      sort out what fitting factors will be applied to pressure ionization.
      if(abs(ifpi).eq.4) then
C        factors to fit Saumon results.
        ifpi_fit = 2
      elseif(ifmodified.gt.0) then
C        factors to fit opal results.
        ifpi_fit = 1
      else
C        unity factors to mimic mdh results as closely as possible.
        ifpi_fit = 0
      endif
      if(ifpi_fit.ne.ifpi_fit_old.and.
     &    (abs(ifpi).eq.3.or.abs(ifpi).eq.4)) then
        ifpi_fit_old = ifpi_fit
C        calculate effective MDH radii for anything-ion interactions
C        (r_ion) and neutral-neutral (r_neutral) interactions.
C        n.b. "neutral" here refers to neutral species and H2+, the
C        only species considered to have non-zero radii in MHD model.
        call effective_radius(ifpi_fit,
     &    bi, nion, nions, r_ion3, r_neutral, nelements)
      endif
C      above 10^5 K, hydrogen is ionized and helium is usually at least
C      first ionized.  In these circumstances it is best to_use_the
C      bare ion as the dv pressure ionization zero rather than the neutral
C      to avoid significance loss.
      if(tl.le.log(1.d5)) then
        ifdvzero = .true.
      else
        ifdvzero = .false.
      endif
      if(ifcoulomb.gt.9) then
        if_dc = 1
        stop 'free_eos_detailed: diffraction correction disabled'
      else
        if_dc = 0
      endif
      if(ifpi.gt.0) then
        ifpl = 1
      else
        ifpl = 0
      endif
      ifpi_local = abs(ifpi)
      if(ifpi_local.gt.4) ifpi_local = 0
      if(lw.lt.3) then
C        full ionization approximation...
C       _use_pteh (e.g., full ionization) approximation for Coulomb sums.
        if_pteh = 1
        if_mc = 0
C       _use_no pressure ionization, since all forms of pressure
C        ionization have zero effect on pressure and entropy for full
C        ionization
        ifpi_local = 0
C        Planck-Larkin reduces to zero for full ionization.
        ifpl = 0
      elseif(ifcoulomb.ge.0) then
        if_mc = 0
        if_pteh = 0
      elseif(ifcoulomb.gt.-10.and.
     &    ((eps(1).gt.0.d0.and.eps(2).gt.0.d0).or.
     &    (lw.eq.3.and.ifmtrace.eq.1))) then
C        non-zero H, and He or metals fully ionized
        if_mc = 1
        if_pteh = 0
      else
        if_mc = 0
        if_pteh = 1
      endif
      ifnoaux_iteration = (ifcoulomb.eq.0.or.if_pteh.eq.1).and.
     &  ifpi_local.le.1
      ifcoulomb_mod = mod(abs(ifcoulomb),10)
      if(ifpi_local.eq.3.or.ifpi_local.eq.4) then
        nextrasum = 9
!        if(ifmodified.le.0) nextrasum = 9
      else
        nextrasum = 0
      endif
C      each entry to free_eos_detailed could potentially have a different
C      set of abundances or a change from zero to non-zero abundance.  
C      check this.
      ifsame_abundances = 1
      ifnear_abundances = 1
      ifsame_zero_abundances = 1
      do i = 1,neps
        if(eps(i).ne.eps_old(i)) then
          if(eps(i).eq.0.d0.or.eps_old(i).eq.0.d0) then
            ifsame_zero_abundances = 0
            ifnear_abundances = 0
          elseif(abs(eps(i)-eps_old(i)).gt.0.05d0*abs(eps(i))) then
            ifnear_abundances = 0
          endif
          ifsame_abundances = 0
          eps_old(i) = eps(i)
        endif
      enddo
      if(ifsame_zero_abundances.ne.1.or.ifmtrace.ne.ifmtrace_old.or.
     &    lw.ne.lw_old.or.ifexcited.ne.ifexcited_old.or.
     &    ifh2.ne.ifh2_olda.or.ifh2plus.ne.ifh2plus_olda) then
        ifh2_olda = ifh2
        ifh2plus_olda = ifh2plus
        ion_stop = 0
        n_partial_ions = 0
        n_partial_elements = 0
C        requires ridiculous values.
        izlo = 1000
        izhi = 0
        do i = 1,nelements
          if(eps(i).gt.0.d0.and.
     &        ((iftracemetal(i).ne.1.and.ifmtrace.ne.1).or.
     &        i.le.2.or.lw.ge.4)) then
C            accumulate data for non-zero eps if any of following conditions are true:
C            1) element not treated as trace metal.
C            2) element is H or He.
C            3) all elements (including trace metals) treated as partially ionized
            ifelement(i) = 1
            n_partial_elements = n_partial_elements + 1
            partial_elements(n_partial_elements) = i
            do ion = ion_stop+1, ion_stop+iatomic_number(i)
C              ifdv(ion) = 1
              n_partial_ions = n_partial_ions + 1
              partial_ions(n_partial_ions) = ion
              inv_ion(ion) = n_partial_ions
              if(i.le.2.or.mod(ifexcited,10).eq.3) then
                iz = ion - ion_stop
                if(iz.lt.izlo.or.iz.gt.izhi) then
C                  new iz value
                  izlo = min(izlo,iz)
                  izhi = max(izhi,iz)
                  if(iz.gt.maxcorecharge)
     &              stop 'free_eos_detailed: internal logic failure'
                  bmin(iz) = bi(ion)
                  nmin(iz) = nmin_species(ion)
                  nmin_max(iz) = nmin_species(ion)
                else
                  bmin(iz) = min(bmin(iz), bi(ion))
                  nmin(iz) = min(nmin(iz), nmin_species(ion))
                  nmin_max(iz) = max(nmin_max(iz), nmin_species(ion))
                endif
              endif
            enddo
          else
            ifelement(i) = 0
            do ion = ion_stop+1, ion_stop+iatomic_number(i)
C              ifdv(ion) = 0
              inv_ion(ion) = 0
            enddo
          endif
          ion_stop = ion_stop + iatomic_number(i)
        enddo
        if(n_partial_elements.eq.0)
     &    stop 'free_eos_detailed: ionized pure metals not implemented'
        max_index = n_partial_ions
        mion_end = partial_elements(n_partial_elements)
C        set ifdv for H2 and H2+
        if(eps(1).gt.0.d0) then
          if(ifh2.gt.0) then
            mion_end = nelements + 1
            ion_end(mion_end) = ion_end(mion_end-1) + 1
            if(ion_end(mion_end).ne.nions+1)
     &        stop 'free_eos_detailed: internal logic failure'
            ifdv(nions+1) = 1
            partial_elements(n_partial_elements+1) = nelements + 1
            max_index = max_index + 1
            partial_ions(n_partial_ions+1) = nions+1
            inv_ion(nions+1) = n_partial_ions+1
C            H2+ required for excited H2.
            if(ifh2plus.gt.0.and.mod(ifexcited,10).gt.1) then
              iz = 1
              ion = nions + 1
              if(iz.lt.izlo.or.iz.gt.izhi) then
C                new iz value
                izlo = min(izlo, iz)
                izhi = max(izhi,iz)
                if(iz.gt.maxcorecharge)
     &            stop 'free_eos_detailed: internal logic failure'
                bmin(iz) = bi(ion)
                nmin(iz) = nmin_species(ion)
                nmin_max(iz) = nmin_species(ion)
              else
                bmin(iz) = min(bmin(iz), bi(ion))
                nmin(iz) = min(nmin(iz), nmin_species(ion))
                nmin_max(iz) = max(nmin_max(iz), nmin_species(ion))
              endif
            endif
          else
            ifdv(nions+1) = 0
            inv_ion(nions+1) = 0
          endif
          if(ifh2plus.gt.0) then
            mion_end = nelements + 2
            ion_end(mion_end) = ion_end(mion_end-1) + 1
            if(ion_end(mion_end).ne.nions+2)
     &        stop 'free_eos_detailed: internal logic failure'
            ifdv(nions+2) = 1
            partial_elements(n_partial_elements+2) = nelements + 2
            max_index = max_index + 1
            partial_ions(n_partial_ions+2) = nions+2
            inv_ion(nions+2) = n_partial_ions+2
            if(mod(ifexcited,10).gt.1) then
              iz = 2
              ion = nions + 2
              if(iz.lt.izlo.or.iz.gt.izhi) then
C                new iz value
                izlo = min(izlo, iz)
                izhi = max(izhi,iz)
                if(iz.gt.maxcorecharge)
     &            stop 'free_eos_detailed: internal logic failure'
                bmin(iz) = bi(ion)
                nmin(iz) = nmin_species(ion)
                nmin_max(iz) = nmin_species(ion)
              else
                bmin(iz) = min(bmin(iz), bi(ion))
                nmin(iz) = min(nmin(iz), nmin_species(ion))
                nmin_max(iz) = max(nmin_max(iz), nmin_species(ion))
              endif
            endif
          else
            ifdv(nions+2) = 0
            inv_ion(nions+2) = 0
          endif
        else
          ifdv(nions+1) = 0
          ifdv(nions+2) = 0
          inv_ion(nions+1) = 0
          inv_ion(nions+2) = 0
        endif
      endif
C      n.b. the conditions on ifmtrace and lw arise because partial_elements
C      may change due to these flags.
      if(ifsame_abundances.ne.1.or.ifmtrace.ne.ifmtrace_olda.or.
     &    lw.ne.lw_olda.or.if_mc.ne.if_mc_olda) then
C        these are distingushed from other "old" flags with "a" suffix
        ifmtrace_olda = ifmtrace
        lw_olda = lw
        if_mc_olda = if_mc
C        calculate total number of electrons available to be freed.
C        n.b. done once each call in case of abundance changes.
C        full_sum0 is full ionization approximation to sum0/rho*NA,
C        where sum0 is the sum over positive ion number densities,
C        rho is the density, and NA is Avogadro's number.
        full_sum0 = 0.d0
C        full_sum2 is full ionization approximation to (sum2-ne*thetae)/rho*NA
C        where sum2-ne*thetae = sum over positive ion number densities times
C        the charge squared on those ions.
        full_sum2 = 0.d0
C        nux, nuy, nuz are number/volume of maximum possible ionization 
C        electrons (divided by rho*NA) for hydrogen, helium, and metals.
        nux = eps(1)
        nuy = 2.d0*eps(2)
        nuz = 0.d0
        do i = 1,nelements
          if(eps(i).gt.0.d0) then
            full_sum0 = full_sum0 + eps(i)
            full_sum2 = full_sum2 + eps(i)*
     &        dble(iatomic_number(i)*iatomic_number(i))
            if(i.gt.2) nuz = nuz + eps(i)*dble(iatomic_number(i))
          endif
        enddo
C        define to ridiculous values as kludgey check that these values
C        never used unless if_mc is 1 (in which case, good values defined).
        hcon_mc = 1.d70
        hecon_mc = 1.d70
        if(if_mc.eq.1) then
C          n.b. other flags supersede if_mc.  only accumulate metal
C          Coulomb constants for metals that are allowed to be
C          partially ionized.  Coulomb effect of fully ionized metals 
C          accounted for with different code.
C          either H and He non-zero abundance or metals fully ionized
C          if latter case, then hcon_mc and hecon_mc will be zero.
          hcon_mc = 0.d0
          hecon_mc = 0.d0
          do index = 1, n_partial_elements
            ielement = partial_elements(index)
            if(ielement.gt.2) then
              hcon_mc = hcon_mc + eps(ielement)
              hecon_mc = hecon_mc + eps(ielement)*
     &          dble(iatomic_number(ielement)*
     &          iatomic_number(ielement))
            endif
          enddo
          hecon_mc = hecon_mc - hcon_mc
          if(eps(2).gt.0.d0) then
            hecon_mc = hecon_mc/eps(2)
          elseif(hecon_mc.ne.0.d0) then
            stop 'free_eos_detailed: hecon_mc logic screwup'
          endif
          if(eps(1).gt.0.d0) then
            hcon_mc = hcon_mc/eps(1)
          elseif(hcon_mc.ne.0.d0) then
            stop 'free_eos_detailed: hcon_mc logic screwup'
          endif
        endif
C        full_sum1 is number/volume of maximum possible ionization electrons
C        divided by rho*NA.
        full_sum1 = nux + nuy + nuz
      endif
C      set up auxiliary variable logic
C      h, hd, he
      if(if_mc.eq.1.or.ifpi_local.eq.2) then
        if(eps(1).gt.0.d0) then
          ifaux(1) = 1
        else
          ifaux(1) = 0
        endif
        if(eps(2).gt.0.d0) then
          ifaux(2) = 1
          ifaux(3) = 1
        else
          ifaux(2) = 0
          ifaux(3) = 0
        endif
C        this somewhat redundant when kif=2 since an auxiliary variable is
C        equal to the match_variable, but it works without any special
C        programming so I won't try to reduce number of auxiliary variables
C        in this special case.
        ifaux(4) = 1
      else
        ifaux(1) = 0
        ifaux(2) = 0
        ifaux(3) = 0
        ifaux(4) = 0
      endif
      if(if_mc.eq.1.and.ifh2plus.gt.0) then
C        need h2plus for if_mc approximation for Coulomb sums.
        ifaux(5) = 1
      else
        ifaux(5) = 0
      endif
C      sum0, sum2
      if(if_pteh.eq.1.or.if_mc.eq.1) then
        ifaux(6) = 0
        ifaux(7) = 0
      else
        ifaux(6) = 1
        ifaux(7) = 1
      endif
C      extrasum
      do iaux = iextraoff+1,iextraoff+nextrasum
        ifaux(iaux) = 1
      enddo
      do iaux = iextraoff+1+nextrasum,iextraoff+maxnextrasum
        ifaux(iaux) = 0
      enddo
C      xextrasum
      if((ifpi_local.eq.3.or.ifpi_local.eq.4).and.ifexcited.gt.0) then
        do iaux = iextraoff+maxnextrasum+1,iextraoff+maxnextrasum+4
          ifaux(iaux) = 1
        enddo
      else
        do iaux = iextraoff+maxnextrasum+1,iextraoff+maxnextrasum+4
          ifaux(iaux) = 0
        enddo
      endif
C      fl is nominally the nauxth auxiliary variable, but it is treated
C      in a different manner than the other auxiliary variables.
      ifaux(naux) = 0
      n_partial_aux = 0
      do iaux = 1,naux
        if(ifaux(iaux).eq.1) then
          n_partial_aux = n_partial_aux + 1
          partial_aux(n_partial_aux) = iaux
          inv_aux(iaux) = n_partial_aux
        else
          inv_aux(iaux) = 0
        endif
      enddo
      if(kif.eq.0) then
        njacobian_final = n_partial_aux
      else
        njacobian_final = n_partial_aux + 1
      endif
      t = exp(tl)
      tc2 = c2/t
C      calculate planck-larkin occupation probabilities and equilibrium
C      constant changes.
      if(ifpl.eq.1)
     &  call pl_prepare(partial_ions, n_partial_ions, max_index,
     &  nions, nion, ifdv,
     &  bi, tc2, plop, plopt, plopt2, dv_pl, dv_plt)
      if(ifrad.ge.1) then
        pr = prad_const*t**4
      else
        pr = 0.d0
      endif
C      criterion for doing warm start with previous auxiliary
C      variables.
C      Cannot get into much trouble if both fl and flold < -10 so
C      relax delta fl criterion in that case.
      
      if(ifnear_abundances.eq.1.and.abs(tlold-tl).le.0.05001d0.and.
     &    (abs(flold-fl).le.0.50001d0.or.
     &    (max(flold,fl).lt.-10.d0.and.abs(flold-fl).le.20.d0)).and.
     &    lw.eq.lw_old.and.
     &    ifexcited.eq.ifexcited_old.and.
     &    nmax.eq.nmax_old.and.
     &    ifh2.eq.ifh2_old.and.
     &    ifh2plus.eq.ifh2plus_old.and.
     &    morder.eq.morder_old.and.
     &    ifmtrace.eq.ifmtrace_old.and.
     &    ifcoulomb.eq.ifcoulomb_old.and.
     &    ifmodified.eq.ifmodified_old) then
        ifwarm = .true.
      else
        ifwarm = .false.
      endif
      if(debug_any) then
C        always_use_cold start when testing derivatives
        ifwarm = .false.
      endif
      if(debug_output.and..not.ifwarm) then
        write(0,*) 'ifwarm is .false.'
        write(0,*) 'tlold, tl, abs(tlold-tl) = ',
     &    tlold, tl, abs(tlold-tl)
        write(0,*) 'flold, fl, abs(flold-fl) = ',
     &    flold, fl, abs(flold-fl)
        write(0,*) 'lw_old, lw = ', lw_old, lw
        write(0,*) 'ifexcited_old, ifexcited = ',
     &    ifexcited_old, ifexcited
        write(0,*) 'nmax_old, nmax = ', nmax_old, nmax
        write(0,*) 'ifh2_old, ifh2 = ', ifh2_old, ifh2
        write(0,*) 'ifh2plus_old, ifh2plus = ',
     &    ifh2plus_old, ifh2plus
        write(0,*) 'morder_old, morder = ', morder_old, morder
        write(0,*) 'ifmtrace_old, ifmtrace = ',
     &    ifmtrace_old, ifmtrace
        write(0,*) 'ifcoulomb_old, ifcoulomb = ',
     &    ifcoulomb_old, ifcoulomb
        write(0,*) 'ifmodified_old, ifmodified = ',
     &    ifmodified_old, ifmodified
      endif
C      Save flags for the next time the above test is done.
C      flold and tlold saved later.
      lw_old = lw
      ifexcited_old = ifexcited
      nmax_old = nmax
      ifh2_old = ifh2
      ifh2plus_old = ifh2plus
      morder_old = morder
      ifmtrace_old = ifmtrace
      ifcoulomb_old = ifcoulomb
      ifmodified_old = ifmodified
      match_variableold = match_variable
      iteration_count = 0
      if(ifnoaux_iteration.or..not.ifwarm) then
C        if no auxiliary variable iteration required or cold start
C        then iterate on fl to make it consistent with match_variable.
C        Initialization for loop convergence criteria
        paaplus = 1.d30
        paaminus = -1.d30
C        mark undefined by ridiculous values
        flplus = -1.d30
        flminus = 1.d30
        lter = 0
        iflast = 0
        paa = 1.d0  !assure at least twice (once if kif=0) through loop.
        pab = 1.d0
        do while(iflast.ne.1)
          lter = lter + 1
C          Newton-Raphson iteration is quadratic, so obtain
C          machine precision (within significance loss noise of say
C          1.d-14) if do one more iteration after 10^-7 convergence.
          if(kif.eq.0.or.lter.ge.ltermax.or.
     &      abs(fl-flold).le.1.d-7) iflast = 1
C          find themodynamically consistent set of fermi-dirac integrals and
C          put the values and their derivatives into rhostar, pstar,
C          sstar, and ustar which are equivalenced to re, pe, se, and ue
C          and their f and t derivatives.  Also calculate exchange (when
C          ifexchange_in > 0) and its effects on dv.
          call fermi_dirac_exchange(fl, tl,
     &      rhostar, pstar, sstar, ustar, 3, morder, ifexchange_in,
     &      dve_exchange, dve_exchangef, dve_exchanget)
          f = exp(fl)
          wf = sqrt(1.d0 + f)
          eta = fl+2.d0*(wf-log(1.d0+wf))
C          number density of free electrons and electron pressure
          n_e = c_e*re
C          rho/mu_e = n_e H = cd*re
          rmue = cd*re
C          set ifionized.
          if(lw.gt.3) then
C            partial ionization with all stages of ionization
C            of all elements.
            ifionized = 0
          elseif(lw.eq.3) then
C            full ionization of all trace metals.
            ifionized = 1
          else
C            full ionization of all elements.
            ifionized = 2
          endif
C          All calls of eos_cold_start should calculate nu = n/(rho*avogadro)
C          form of h2, h2plus, h, hd, he, and extrasum for the
C          ifnoaux_iteration case or kif=0 case.  This allows the pressure
C          and its derivatives to be calculated properly during the
C          eos_cold_start iteration on fl for kif=1, affects nothing
C          else that is used during the iterations, and gives the right
C          result for entropy, etc., after the iterations are completed for
C          the ifnoaux_iteration case.  (For the
C          .not.ifnoaux_iteration.and.kif.eq.1 case, a final call to
C          eos_cold_start is required with ifnuform = .false., see below).
          if(ifnoaux_iteration.or.kif.eq.1) then
            ifnuform = .true.
          else
            ifnuform = .false.
          endif
          iteration_count = iteration_count + 1
C          N.B. The actual eos_cold_start EOS calculation is done adopting
C          no or PTEH pressure ionization and adopting the PTEH approximation
C          for the Coulomb sums.  However, the sums resulting from that
C          EOS calculation are calculated depending on the values of
C          if_pteh, if_mc, and ifpi_local which depend on the overall
C          free-energy model, not the special model used for the eos_cold_start
C          calculation.
C          if near fl convergence, make sure to fix underflow criterion
C          in ionize called by eos_calc.
C          note paa = 1. (and pab = 1.) on first time through so always false
C          on first time
          if(kif.gt.0.and.max(abs(paa),abs(pab)).le.1.d-2) then
            ifsame_under = .true.
          else
            ifsame_under = .false.
          endif
          call eos_cold_start(
     &      ifsame_under, lambda, gamma_e,
     &      partial_ions, f, eta, wf, t, n_e, pstar, 
     &      dve_exchange, dve_exchangef, dve_exchanget,
     &      full_sum0, full_sum1, full_sum2, charge, dv_pl, dv_plt,
     &      ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &      ifsame_zero_abundances, ifexcited, ifnuform,
     &      naux, inv_aux,
     &      inv_ion, max_index,
     &      partial_elements, n_partial_elements, ion_end,
     &      ifionized, if_pteh, if_mc, ifreducedmass,
     &      ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,
     &      ifpl, ifmodified, ifh2, ifh2plus,
     &      izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &      eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &      r_ion3, r_neutral, nelements, 
     &      ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &      rhostar,
     &      nextrasum,
     &      ne, nef, net, sion, sionf, siont, uion,
     &      sumpl1, sumpl1f, sumpl1t, sumpl2,
     &      h2, h2f, h2t,
     &      aux, auxf, auxt)
C          abundances don't change within auxiliary *AND*
C          match_variable loops.
          ifsame_abundances = 1
          ifnear_abundances = 1
C          update rho to be consistent with calculated rl
          rho = exp(rl)
          if(if_mc.eq.1) then
C            in just this case
C            ionize and eos_calc have only calculated fully
C            ionized part/(rho*avogadro).
C            save it.
            sum0_mc = sum0
            sum2_mc = sum2
          endif
          tlold = tl
          flold = fl
          if(kif.eq.0) then
          else
            if(kif.eq.1) then
C              there is the possibility of pressure ionization and 
C              Planck-Larkin terms
              if(ifpi_local.eq.0) then
C                there are no p, s, u terms from pressure ionization.
                ppi = 0.d0
                ppif = 0.d0
              else
                call pteh_pi_end(
     &            full_sum1, rho, rf, rt, t, ne, nef, net,
     &            ppi, ppif, ppit, spi, spif, spit, upi)
              endif
              if(ifexcited.gt.0.and.ifpl.eq.1.and.(
     &            ifpi_local.eq.3.or.ifpi_local.eq.4)) then
C                For the above combination of flags, the excitation_pi
C                code evaluates non-zero values of psum and its
C                derivatives so that excitation_pi_end
C                returns non-zero values of pexcited and its derivatives.
C                n.b. for .not.ifwarm, the simplified free-energy model
C                used to calculate the equilibrium constants, dv, is different
C                from the full free-energy model used to calculate psum,
C                etc., in excitation_sum.
                call excitation_pi_end(t, rho, rf, rt,
     &            pexcited, pexcitedf, pexcitedt,
     &            sexcited, sexcitedf, sexcitedt, uexcited)
              else
C                otherwise, excitation_pi_end would return zero for
C                pexcited and its derivatives so save some time by not
C                calling it.
                pexcited = 0.d0
                pexcitedf = 0.d0
              endif
C              All the Coulomb stuff is calculated using if_pteh = 1 for this
C              case.
C              PTEH (full ionization) approximation to sum0, sum2
C              reassert this approximation because eos_calc
C              calls ionize which messes a bit with sum0 and sum2
              sum0ne = full_sum0/full_sum1
              sum0 = sum0ne*n_e
              sum0f = 0.d0
              sum0t = 0.d0
              sum2ne = full_sum2/full_sum1
              sum2 = sum2ne*n_e
              sum2f = 0.d0
              sum2t = 0.d0
              call master_coulomb_end(rhostar,
     &          sum0, sum0ne, sum0f, sum0t,
     &          sum2, sum2ne, sum2f, sum2t,
     &          n_e, t, pstar,
     &          ifcoulomb_mod, if_dc, 1,
     &          dpcoulomb, dpcoulombf, dpcoulombt,
     &          dscoulomb, dscoulombf, dscoulombt, ducoulomb)
              call exchange_end(rhostar, pstar, 3,
     &          pex, pext, pexf,
     &          sex, sexf, sext, uex)
              p0 = cr*rho*t
              ni = full_sum0 - (h2+h2plus)
              pion = ni*p0
              pionf = ni*p0*rf - (h2f+h2plusf)*p0
              pe_cgs = cpe*pe
              pnorad = pe_cgs + pion + ppi + pexcited +
     &          dpcoulomb + pex
              p = pnorad + pr
              if(p.le.0.d0) then
C                write(0,*) 'match_variable/ln(10), fl, tl/ln(10) ='
C                write(0,'(1p5d25.15)') match_variable/log(10.d0), fl,
C     &            tl/log(10.d0)
C                write(0,*) 'eps ='
C                write(0,'(1p5d25.15)') eps
C                write(0,*) 'free_eos_detailed: (1) '//
C     &            'negative p calculated'
                info = 1
                return
              endif
C              d ln p(f,t)/d ln f
              pf = (pe_cgs*pef + pionf + ppif + pexcitedf +
     &          dpcoulombf + pexf)/p
              if(ifrad.gt.1) then
                paa = log(pnorad) - match_variable
C                the p/pnorad factor converts ln ptotal derivative
C                to ln pnorad derivative.
                pnoradf = pf*p/pnorad
                pac = -paa/(pnoradf)
              else
                paa = log(p) - match_variable
                pac = -paa/pf
              endif
            elseif(kif.eq.2) then
              paa = log(rho) - match_variable
              pac = -paa/rf
            else
              stop 'free_eos_detailed: kif must be 0, 1, or 2'
            endif
C            Apply limits to potential fl change.
            if(fl + pac .lt. -10.d0) then
C             Cannot get into much trouble for derived  fl < -10.
              pab = min(20.d0,max(-20.d0,pac))
            else
              lnrho_5 = log(rho) - 1.5d0*(tl - log(1.d5))
              if(tl.gt.log(1.d6).or.
     &            lnrho_5 + rf*min(10.d0,pac).lt.log(1.d-3)) then
                pab = min(10.d0,max(-10.d0,pac))
              elseif(ifnoaux_iteration.or.
     &            lnrho_5 + rf*min(3.d0,pac).lt.log(1.d-2)) then
                pab = min(3.d0,max(-3.d0,pac))
              elseif(
     &            lnrho_5 + rf*min(1.d0,pac).lt.log(1.d-1)) then
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
C              have bracket already!
C              make sure step would keep within bracket
              if((flminus.le.fl+pab.and.fl+pab.le.flplus).or.
     &            (flminus.ge.fl+pab.and.fl+pab.ge.flplus)) then
                fl = fl + pab
              else
                fl = 0.5d0*(flminus+flplus)
              endif
            else
C              if no bracket, yet, then normal step is only allowed for
C              when pf (or rf) is positive.  that is local paa
C              is following global trend that paa
C              and pl (or rl) generally increases with fl.
C              if pf (or rf) negative with no bracket, 
C              then move out of region in the direction 
C              which should produce (global) bracket.
              if((kif.eq.1.and.pf.lt.0.d0).or.
     &            (kif.eq.2.and.rf.lt.0.d0)) then
                if(fl.eq.flminus) then
                  pab = 0.5d0
                else
                  pab = -0.5d0
                endif
              endif
              fl = fl + pab
            endif
          endif
          if(debug_output) then
            write(0,'(a,/,i5,1p5d25.15)')
     &        'lter, flold, paa, pac, pab, fl =',
     &        lter, flold, paa, pac, pab, fl
          endif
C          End of iteration on fl to match match_variable with
C          eos_cold_start
        enddo
        if(kif.gt.0.and.lter.ge.ltermax) then
C          write(0,*) 'match_variable/ln(10), fl, tl/ln(10) ='
C          write(0,'(1p5d25.15)') match_variable/log(10.d0), fl,
C     &      tl/log(10.d0)
C          write(0,*) 'eps ='
C          write(0,'(1p5d25.15)') eps
C          write(0,*) 'ERROR: probably a multi-valued EOS '//
C     &      'for the eos_cold_start case'
C          write(0,*) 'flplus, paa(flplus), flminus, paa(flminus) ='
C          write(0,'(1p5d25.15)') flplus, paaplus, flminus, paaminus
C          write(0,'(a,/,i5,1p5d25.15)')
C     &      'lter, flold, paa, pac, pab, fl =',
C     &      lter, flold, paa, pac, pab, fl
C          write(0,*) 'free_eos_detailed: too many fl iterations'
          info = 2
          return
        endif
C      end of if(ifnoaux_iteration.or..not.ifwarm) then branch
      endif
      if(.not.ifnoaux_iteration) then
        lw_old = lw
        ifexcited_old = ifexcited
        nmax_old = nmax
        ifh2_old = ifh2
        ifh2plus_old = ifh2plus
        morder_old = morder
        ifmtrace_old = ifmtrace
        ifcoulomb_old = ifcoulomb
        ifmodified_old = ifmodified
        match_variableold = match_variable
C        set ifionized.
        if(lw.gt.3) then
C          partial ionization with all stages of ionization
C          of all elements.
          ifionized = 0
        elseif(lw.eq.3) then
C          full ionization of all trace metals.
          ifionized = 1
        else
C          full ionization of all elements.
          ifionized = 2
        endif
        if(ifwarm) then
C          warm start:
          
C          must update rho from local variable (equivalenced to
C          an auxiliary variable) because
C          in the argument list it is commented as output only
C          (and the value may not be saved between calls by the
C          calling programme, bug found by jcd.)
C          rf and rt are okay because they are saved local variables
          rho = exp(rl)
C          For warm start, previous calculation of eos was done with
C          eos_tqft which produces nu form of many auxiliary variables.
          do index_aux = 1,n_partial_aux
            iaux = partial_aux(index_aux)
            if(iaux.ne.4.and.iaux.ne.6.and.iaux.ne.7) then
C              for all but rl, sum0, and sum2 ... convert
C              to n form.
              aux(iaux) = aux(iaux)*rho*avogadro
              auxf(iaux) = auxf(iaux)*rho*avogadro + aux(iaux)*rf
              auxt(iaux) = auxt(iaux)*rho*avogadro + aux(iaux)*rt
            endif
C            Logarithmic Taylor series unless auxiliary variable is zero.
            if(iaux.eq.4) then
C              already in log form.
              faux(index_aux) =
     &          (fl-flold)*auxf(iaux) + (tl-tlold)*auxt(iaux)
              if(abs(faux(index_aux)).le.faux_limit_small_start)
     &          aux(iaux) = aux(iaux) + faux(index_aux)
            elseif(aux(iaux).gt.0.d0) then
              faux(index_aux) =
     &          ((fl-flold)*auxf(iaux) + (tl-tlold)*auxt(iaux))/
     &          aux(iaux)
              if(abs(faux(index_aux)).le.faux_limit_small_start)
     &          aux(iaux) = exp(log(aux(iaux)) + faux(index_aux))
            elseif(iaux.gt.maxnextrasum.and.aux(iaux).lt.0.d0) then
              faux(index_aux) =
     &          ((fl-flold)*auxf(iaux) + (tl-tlold)*auxt(iaux))/
     &          aux(iaux)
              if(abs(faux(index_aux)).le.faux_limit_small_start)
     &          aux(iaux) = -exp(log(-aux(iaux)) + faux(index_aux))
            else
              aux(iaux) = 0.d0
            endif
          enddo
        elseif(kif.eq.1) then
C          for cold start must recalculate eos_cold_start with correct
C          ifnuform = .false. for this special case.  See ifnuform
C          shenanigans for kif=1 above.
C          N.B. here ifwarm is false so ifsame_under will be identical to
C          the value used for the above call to eos_cold_start.
          ifnuform = .false.
          iteration_count = iteration_count + 1
          call eos_cold_start(
     &      ifsame_under, lambda, gamma_e,
     &      partial_ions, f, eta, wf, t, n_e, pstar, 
     &      dve_exchange, dve_exchangef, dve_exchanget,
     &      full_sum0, full_sum1, full_sum2, charge, dv_pl, dv_plt,
     &      ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &      ifsame_zero_abundances, ifexcited, ifnuform,
     &      naux, inv_aux,
     &      inv_ion, max_index,
     &      partial_elements, n_partial_elements, ion_end,
     &      ifionized, if_pteh, if_mc, ifreducedmass,
     &      ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,
     &      ifpl, ifmodified, ifh2, ifh2plus,
     &      izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &      eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &      r_ion3, r_neutral, nelements, 
     &      ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &      rhostar,
     &      nextrasum,
     &      ne, nef, net, sion, sionf, siont, uion,
     &      sumpl1, sumpl1f, sumpl1t, sumpl2,
     &      h2, h2f, h2t,
     &      aux, auxf, auxt)
        endif
        tlold = tl
        flold = fl
        do index_aux = 1,n_partial_aux
          iaux = partial_aux(index_aux)
C          initially only allow small magnitude above aux_underflow
C          for zero to non-zero auxiliary variable change.
          zerolim(iaux) = 1.d10*aux_underflow
C          overall underflow trap (exp(-575) = 10^{-249.7})
!          if(abs(aux(iaux)).lt.1.d-250)
!     &      aux(iaux) = 0.d0
        enddo
C        save values for later_use_to start bfgs iteration.
        fl_restore = fl
        do index_aux = 1,n_partial_aux
          iaux = partial_aux(index_aux)
          aux_restore(iaux) = aux(iaux)
        enddo
C        this factor must be unity in order to judge
C        the initial relative change in ne.
        faux_scale = 1.d0
C        initialize logic for counting number of zero trys.
        ioncountzero = -1
        ifzerocount = 0
        faux_limit_small = faux_limit_small_start
        maxfaux = 1.d0    !assure at least two iterations
        maxfaux_diag = 1.d0 !needed to initialize maxfaux_diag_old
        ioncount = 0
        iflast = 0
C        some value that is too large to be astrophysically possible
C        for a free energy per unit mass.  This forces at least two
C        calls to eos_jacobian below before evaluating (for kif = 2)
C        whether free is decreasing or not.
        free_old = 1.d300
C        This allows the possibility of one eos_bfg call in the do loop
C        to settle things down for kif = 2, cold start.
        eos_bfgs_called = .false.
!        eos_bfgs_called = .true.
C        if not straight through (iflast.ne.1) then either warm
C        start (with all auxiliary variables defined from  previous call)
C        or else cold start (with all auxiliary variables defined by
C        above call to eos_cold_start using pteh_pi and pteh sum approximation).

C        For first iteration assume NR solution unless simple iteration
C        solution criteria are fulfilled.
        any_ifsimple = .false.
        do index_aux = 1, njacobian_final
          simple_lambda(index_aux) = 0.d0
        enddo
        do while(iflast.ne.1.and.ioncount.lt.maxioncount)
C          find themodynamically consistent set of fermi-dirac integrals and
C          put the values and their derivatives into rhostar, pstar,
C          sstar, and ustar which are equivalenced to re, pe, se, and ue
C          and their f and t derivatives.  Also calculate exchange (when
C          ifexchange_in > 0) and its effects on dv.
          call fermi_dirac_exchange(fl, tl,
     &      rhostar, pstar, sstar, ustar, 3, morder,
     &      ifexchange_in, dve_exchange, dve_exchangef, dve_exchanget)
          f = exp(fl)
          wf = sqrt(1.d0 + f)
          eta = fl+2.d0*(wf-log(1.d0+wf))
C          number density of free electrons and electron pressure
          n_e = c_e*re
C          rho/mu_e = n_e H = cd*re
          rmue = cd*re
          ioncount = ioncount + 1
          maxfaux_diag_old = maxfaux_diag
C          either cold or warm start may require huge ln changes in 
C          auxiliary variables.  However, these huge ln changes
C          are often only required of auxiliary variables that
C          are approaching zero and which have little practical
C          effect on the equilibrium constants.  If we_use_a
C          over-cautious approach when limiting the maximum
C          ln change, then a large number of auxiliary variable
C          iterations will be required.  However, with our present
C          approach of allowing large ln changes, the maximum
C          number of iterations is usually less than 10.
C          n.b._use_one extra iteration after satisfy ending criteria.
          if(.not.any_ifsimple.and.
     &      ((ifwarm.and.ioncount.ge.maxioncount).or.
     &      ioncount.ge.maxioncount.or.
     &      abs(maxfaux).le.eps_aux))
     &      iflast = 1
C          if near final stages of auxiliary variable convergence,
C          make sure to fix underflow criterion in ionize called by eos_calc.
!          if(ioncount.ge.2.and.
!     &        max(abs(rl-rl_old),
!     &        abs(ne-ne_old)/ne,
!     &        abs(lambda-lambda_old)/max(1.d-15,lambda))/
!     &        faux_scale.le.1.d-2) then
C          N.B. maxfaux = 1. on first time through this loop so ifsame_under
C          always false on first time.
          if(ioncount.ge.2.and.
     &        abs(maxfaux).le.1.d-3) then
            ifsame_under = .true.
          else
            ifsame_under = .false.
          endif
          if(.false..and.debug_output) then
            write(0,*) 'aux'
            write(0,*) (aux(partial_aux(index_aux)),
     &        index_aux = 1, n_partial_aux)
          endif
          if(debug_dvdaux.or.debug_jacobian) then
C            specify change in auxiliary variables as a function of fl and tl
C            leave dv zero point as is.
            ifdvzero = .true.
C            Choose itemp1 and itemp2 to correspond to auxiliary variables
C            actually used for particular free-energy model being tested.
C            For example, EOS1 has 15 auxiliary variables in order
C            sum0, sum2, 7 neutral auxiliary variables, 2 ion auxiliary
C            variables, and 4 xextrasum auxiliary variables.
            if(.false.) then
              itemp1 = 1
              itemp2 = 2
            elseif(.false.) then
              itemp1 = 3
              itemp2 = 4
            elseif(.false.) then
              itemp1 = 5
              itemp2 = 6
            elseif(.false.) then
              itemp1 = 7
              itemp2 = 8
            elseif(.false.) then
C              repeat 8 for alignment of like results.  8 and 9 are
C              last two derivatives wrt neutral sums.
              itemp1 = 8
              itemp2 = 9
            elseif(.false.) then
              itemp1 = 10
              itemp2 = 11
            elseif(.false.) then
              itemp1 = 12
              itemp2 = 13
            elseif(.false.) then
              itemp1 = 14
              itemp2 = 15
            elseif(.true.) then
              itemp1 = n_partial_aux + 1
              itemp2 = n_partial_aux + 2
            endif
            if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
              iaux1 = partial_aux(itemp1)
              if(iaux1.eq.4) then
                aux(iaux1) = aux(iaux1) + delta1 !density
                rho = exp(rl)
              else
                aux(iaux1) = aux(iaux1)*(1.d0 + delta1)
              endif
            elseif(itemp1.eq.n_partial_aux + 1) then
              fl = fl + delta1
            elseif(itemp1.eq.n_partial_aux + 2) then
              tl = tl + delta1
            else
              stop 'free_eos_detailed: bad itemp1'
            endif
            if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
              iaux2 = partial_aux(itemp2)
              if(iaux2.eq.4) then
                aux(iaux2) = aux(iaux2) + delta2 !density
                rho = exp(rl)
              else
                aux(iaux2) = aux(iaux2)*(1.d0 + delta2)
              endif
            elseif(itemp2.eq.n_partial_aux + 1) then
              fl = fl + delta2
            elseif(itemp2.eq.n_partial_aux + 2) then
              tl = tl + delta2
            else
              stop 'free_eos_detailed: bad itemp2'
            endif
C            just in case there was a change in tl
            t = exp(tl)
            tc2 = c2/t
C            calculate planck-larkin occupation probabilities and equilibrium
C            constant changes.
            if(ifpl.eq.1)
     &        call pl_prepare(partial_ions, n_partial_ions, max_index,
     &        nions, nion, ifdv,
     &        bi, tc2, plop, plopt, plopt2, dv_pl, dv_plt)
            if(ifrad.ge.1) then
              pr = prad_const*t**4
            else
              pr = 0.d0
            endif
C            just in case there was a change in fl or tl.
            call fermi_dirac_exchange(fl, tl,
     &        rhostar, pstar, sstar, ustar, 3, morder,
     &        ifexchange_in, dve_exchange, dve_exchangef,
     &        dve_exchanget)
            f = exp(fl)
            wf = sqrt(1.d0 + f)
            eta = fl+2.d0*(wf-log(1.d0+wf))
C            number density of free electrons and electron pressure
            n_e = c_e*re
C            rho/mu_e = n_e H = cd*re
            rmue = cd*re
          elseif(debug_dauxdv.or.debug_bfgs) then
C            choose first two dv or x values to vary
C            partial derivative of auxiliary variables wrt dv
C            or partial derivative of y wrt x.
C            1 ==> 4 indices correspond to H+, He+, He++, C+, while
C            max_index-1 ==> max_index correspond to H2 and H2+.
            if(.true.) then
              itemp1 = 2
              itemp2 = 3
            elseif(.true.) then
              itemp1 = 3
              itemp2 = 4
            elseif(.true.) then
              itemp1 = max_index-1
              itemp2 = max_index
            endif
C            the next part of this initial step for this debugging mode is
C            done inside eos_jacobian after consistent dv values are formed.
          endif
          rl_old = rl
          ne_old = ne
          lambda_old = lambda
          do index_aux = 1,n_partial_aux
            iaux = partial_aux(index_aux)
            aux_old(iaux) = aux(iaux)
          enddo
          flold = fl
          if(kif.eq.0.and.iflast.ne.1) then
C            No need to calculate fl and tl derivatives until last
C            iteration for kif = 0
            ifnr = 1
          else
C            Calculate combined auxiliary variable and fl (and ft)
C            partial derivatives.
            ifnr = 3
          endif
          if(.true..and.
     &        kif.eq.2.and.ioncount.eq.100) then
C            write(0,*)
C     &        'start one-time emergency eos_bfgs per NR cycle'
            if(.true.) then
            fl = fl_restore
            do index_aux = 1,n_partial_aux
              iaux = partial_aux(index_aux)
              aux_old(iaux) = aux_restore(iaux)
            enddo
            do index_aux = 1,n_partial_aux
              iaux = partial_aux(index_aux)
C              n.b. 4th auxiliary variable (rl = log(rho)) is already in 
C              log form.
C              N.B. "allow_log" connotes "allow taking log of auxiliary variable
C              delivered by eos_jacobian"
              allow_log(index_aux) = iaux.ne.4.and.
     &          aux_old(iaux).gt.0.d0
C              sign of xextrasum depends on extrasum.  it's okay if negative,
C              but subsequent logic must adjust for this.
C              N.B. "allow_log_neg" connotes "allow taking log of *negative* of
C              auxiliary variable delivered by eos_jacobian"
              allow_log_neg(index_aux) = iaux.gt.maxnextrasum.and.
     &          aux_old(iaux).lt.0.d0
            enddo
            isimple_bfgs = 1
C            ordinarily do not do any preliminary simple iterations at fixed
C            fl prior to using the bfgs technique.
            do while(isimple_bfgs.le.0)
            !temporary.  
            !do while(isimple_bfgs.le.5)
              call eos_jacobian(
     &          allow_log, allow_log_neg, partial_aux,
     &          ifrad, debug_dauxdv, itemp1, itemp2, delta1, delta2,
     &          debug_dvdaux, debug_jacobian,
     &          degeneracy, pressure, density,
     &          energy, enthalpy, entropy,
     &          iaux1, iaux2,
     &          match_variable, match_variable_save, kif, fl, tl_save,
     &          aux_old, aux, auxf, auxt, daux_dv,
     &          njacobian, njacobian_final,
     &          rhs, jacobian, p, pr,
     &          ddv_aux, sumion0, sumion2,
     &          lambda, gamma_e,
     &          ifnr, nux, nuy, nuz, mion_end,
     &          n_partial_ions, n_partial_aux,
     &          sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &          partial_ions, f, eta, wf, t, n_e, pstar, 
     &          dve_exchange, dve_exchangef, dve_exchanget,
     &          full_sum0, full_sum1, full_sum2, charge, charge2,
     &          dv_pl, dv_plt,
     &          ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &          ifsame_zero_abundances,
     &          ifexcited, ifsame_under,
     &          naux, inv_aux,
     &          inv_ion, max_index,
     &          partial_elements, n_partial_elements, ion_end,
     &          ifionized, if_pteh, if_mc, ifreducedmass,
     &          ifsame_abundances, ifmtrace,
     &          iatomic_number, ifpi_local,
     &          ifpl, ifmodified, ifh2, ifh2plus,
     &          izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &          eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &          r_ion3, r_neutral, nelements, 
     &          ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &          rhostar,
     &          ne, nef, net, sion, sionf, siont, uion,
     &          pnorad, pnoradf, pnoradt, pnorad_daux,
     &          free_fp, free_daux,
     &          nextrasum, maxnextrasum,
     &          sumpl1, sumpl1f, sumpl1t, sumpl2,
     &          h2, h2f, h2t, h2_dv)
              if(debug_output) then
                write(0,*) 'pre-bfgs isimple_bfgs = ', isimple_bfgs
                write(0,*)
     &            'old pre-bfgs value of aux for frozen fl'
                write(0,'(5(1pd14.5,2x))')
     &            (aux_old(partial_aux(index_aux)),
     &            index_aux = 1, n_partial_aux)
                write(0,*)
     &            'new pre-bfgs value of aux for frozen fl'
                write(0,'(5(1pd14.5,2x))')
     &            (aux(partial_aux(index_aux)),
     &            index_aux = 1, n_partial_aux)
              endif
              do index_aux = 1,n_partial_aux
                iaux = partial_aux(index_aux)
                aux_old(iaux) = aux(iaux)
              enddo
              do index_aux = 1,n_partial_aux
                iaux = partial_aux(index_aux)
C                n.b. 4th auxiliary variable (rl = log(rho)) is already in 
C                log form.
C                N.B. "allow_log" connotes "allow taking log of auxiliary variable
C                delivered by eos_jacobian"
                allow_log(index_aux) = iaux.ne.4.and.
     &            aux_old(iaux).gt.0.d0
C                sign of xextrasum depends on extrasum.  it's okay if negative,
C                but subsequent logic must adjust for this.
C                N.B. "allow_log_neg" connotes "allow taking log of *negative* of
C                auxiliary variable delivered by eos_jacobian"
                allow_log_neg(index_aux) = iaux.gt.maxnextrasum.and.
     &            aux_old(iaux).lt.0.d0
              enddo
              isimple_bfgs = isimple_bfgs + 1
            enddo
            endif
            call eos_bfgs(
     &        allow_log, allow_log_neg, partial_aux,
     &        ifrad, debug_dauxdv, itemp1, itemp2, delta1, delta2,
     &        debug_dvdaux, debug_jacobian,
     &        degeneracy, pressure, density, energy, enthalpy, entropy,
     &        iaux1, iaux2,
     &        match_variable, match_variable_save, kif, fl, tl_save,
     &        aux_old, aux, auxf, auxt, daux_dv,
     &        njacobian, njacobian_final,
     &        rhs, jacobian, p, pr,
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
     &        ifsame_abundances, ifmtrace,
     &        iatomic_number, ifpi_local,
     &        ifpl, ifmodified, ifh2, ifh2plus,
     &        izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &        eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &        r_ion3, r_neutral, nelements, 
     &        ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &        ne, nef, net, sion, sionf, siont, uion,
     &        pnorad, pnoradf, pnoradt, pnorad_daux,
     &        free_fp, free_daux,
     &        nextrasum, maxnextrasum,
     &        sumpl1, sumpl1f, sumpl1t, sumpl2,
     &        h2, h2f, h2t, h2_dv, info)
            if(info.ne.0) return
            call fermi_dirac_exchange(fl, tl,
     &        rhostar, pstar, sstar, ustar, 3, morder,
     &        ifexchange_in, dve_exchange, dve_exchangef,
     &        dve_exchanget)
            flold = fl
            f = exp(fl)
            wf = sqrt(1.d0 + f)
            eta = fl+2.d0*(wf-log(1.d0+wf))
C            number density of free electrons and electron pressure
            n_e = c_e*re
C            rho/mu_e = n_e H = cd*re
            rmue = cd*re
            isimple_bfgs = 1
C            ordinarily do not do any post-BFGS simple iterations at fixed
C            fl prior to using the NR technique.
            do while(isimple_bfgs.le.0)
            !temporary.  
            !do while(isimple_bfgs.le.5)
              call eos_jacobian(
     &          allow_log, allow_log_neg, partial_aux,
     &          ifrad, debug_dauxdv, itemp1, itemp2, delta1, delta2,
     &          debug_dvdaux, debug_jacobian,
     &          degeneracy, pressure, density,
     &          energy, enthalpy, entropy,
     &          iaux1, iaux2,
     &          match_variable, match_variable_save, kif, fl, tl_save,
     &          aux_old, aux, auxf, auxt, daux_dv,
     &          njacobian, njacobian_final,
     &          rhs, jacobian, p, pr,
     &          ddv_aux, sumion0, sumion2,
     &          lambda, gamma_e,
     &          ifnr, nux, nuy, nuz, mion_end,
     &          n_partial_ions, n_partial_aux,
     &          sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &          partial_ions, f, eta, wf, t, n_e, pstar, 
     &          dve_exchange, dve_exchangef, dve_exchanget,
     &          full_sum0, full_sum1, full_sum2, charge, charge2,
     &          dv_pl, dv_plt,
     &          ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &          ifsame_zero_abundances,
     &          ifexcited, ifsame_under,
     &          naux, inv_aux,
     &          inv_ion, max_index,
     &          partial_elements, n_partial_elements, ion_end,
     &          ifionized, if_pteh, if_mc, ifreducedmass,
     &          ifsame_abundances, ifmtrace,
     &          iatomic_number, ifpi_local,
     &          ifpl, ifmodified, ifh2, ifh2plus,
     &          izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &          eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &          r_ion3, r_neutral, nelements, 
     &          ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &          rhostar,
     &          ne, nef, net, sion, sionf, siont, uion,
     &          pnorad, pnoradf, pnoradt, pnorad_daux,
     &          free_fp, free_daux,
     &          nextrasum, maxnextrasum,
     &          sumpl1, sumpl1f, sumpl1t, sumpl2,
     &          h2, h2f, h2t, h2_dv)
              if(debug_output) then
                write(0,*) 'post-bfgs isimple_bfgs = ', isimple_bfgs
                write(0,*)
     &            'old post-bfgs value of aux for frozen fl'
                write(0,'(5(1pd14.5,2x))')
     &            (aux_old(partial_aux(index_aux)),
     &            index_aux = 1, n_partial_aux)
                write(0,*)
     &            'new post-bfgs value of aux for frozen fl'
                write(0,'(5(1pd14.5,2x))')
     &            (aux(partial_aux(index_aux)),
     &            index_aux = 1, n_partial_aux)
              endif
              do index_aux = 1,n_partial_aux
                iaux = partial_aux(index_aux)
                aux_old(iaux) = aux(iaux)
              enddo
              do index_aux = 1,n_partial_aux
                iaux = partial_aux(index_aux)
C                n.b. 4th auxiliary variable (rl = log(rho)) is already in 
C                log form.
C                N.B. "allow_log" connotes "allow taking log of auxiliary variable
C                delivered by eos_jacobian"
                allow_log(index_aux) = iaux.ne.4.and.
     &            aux_old(iaux).gt.0.d0
C                sign of xextrasum depends on extrasum.  it's okay if negative,
C                but subsequent logic must adjust for this.
C                N.B. "allow_log_neg" connotes "allow taking log of *negative* of
C                auxiliary variable delivered by eos_jacobian"
                allow_log_neg(index_aux) = iaux.gt.maxnextrasum.and.
     &            aux_old(iaux).lt.0.d0
              enddo
              isimple_bfgs = isimple_bfgs + 1
            enddo
C            restore initial size of faux_limit_small.
            faux_limit_small = faux_limit_small_start
C            For first post-bfgs iteration assume NR solution unless
C            simple iteration solution criteria are fulfilled.
            any_ifsimple = .false.
            do index_aux = 1, njacobian_final
              simple_lambda(index_aux) = 0.d0
            enddo
C            force at least two NR iterations after bfgs minimization.
            iflast = 0
C            write(0,*)
C     &        'end one-time emergency eos_bfgs per NR cycle'
C            only call preliminary eos_bfgs once per Newton-Raphson
C            iteration cycle.
            eos_bfgs_called = .true.
            if(debug_bfgs) then
C              y and g are no longer available to this routine so comment
C              out these test results unless and until they are needed again.
!              degeneracy(1) = y
!              degeneracy(2) = g(itemp1)
!              degeneracy(3) = g(itemp2)
!              pressure(1) = y
!              pressure(2) = g(itemp1)
!              pressure(3) = g(itemp2)
!              density(1) = y
!              density(2) = g(itemp1)
!              density(3) = g(itemp2)
!              energy(1) = y
!              energy(2) = g(itemp1)
!              energy(3) = g(itemp2)
!              enthalpy(1) = y
!              enthalpy(2) = g(itemp1)
!              enthalpy(3) = g(itemp2)
!              entropy(1) = y
!              entropy(2) = g(itemp1)
!              entropy(3) = g(itemp2)
              match_variable = match_variable_save
              if(kif.eq.0) fl = match_variable_save
              tl = tl_save
              return
            endif
          endif
          do index_aux = 1,n_partial_aux
            iaux = partial_aux(index_aux)
C            n.b. 4th auxiliary variable (rl = log(rho)) is already in 
C            log form.
C            N.B. "allow_log" connotes "allow taking log of auxiliary variable
C            delivered by eos_jacobian"
            allow_log(index_aux) = iaux.ne.4.and.
     &        aux_old(iaux).gt.0.d0
C            sign of xextrasum depends on extrasum.  it's okay if negative,
C            but subsequent logic must adjust for this.
C            N.B. "allow_log_neg" connotes "allow taking log of *negative* of
C            auxiliary variable delivered by eos_jacobian"
            allow_log_neg(index_aux) = iaux.gt.maxnextrasum.and.
     &        aux_old(iaux).lt.0.d0
          enddo
          iteration_count = iteration_count + 1
          call eos_jacobian(
     &      allow_log, allow_log_neg, partial_aux,
     &      ifrad, debug_dauxdv, itemp1, itemp2, delta1, delta2,
     &      debug_dvdaux, debug_jacobian,
     &      degeneracy, pressure, density, energy, enthalpy, entropy,
     &      iaux1, iaux2,
     &      match_variable, match_variable_save, kif, fl, tl_save,
     &      aux_old, aux, auxf, auxt, daux_dv,
     &      njacobian, njacobian_final,
     &      rhs, jacobian, p, pr,
     &      ddv_aux, sumion0, sumion2,
     &      lambda, gamma_e,
     &      ifnr, nux, nuy, nuz, mion_end,
     &      n_partial_ions, n_partial_aux,
     &      sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &      partial_ions, f, eta, wf, t, n_e, pstar, 
     &      dve_exchange, dve_exchangef, dve_exchanget,
     &      full_sum0, full_sum1, full_sum2, charge, charge2,
     &      dv_pl, dv_plt,
     &      ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &      ifsame_zero_abundances,
     &      ifexcited, ifsame_under,
     &      naux, inv_aux,
     &      inv_ion, max_index,
     &      partial_elements, n_partial_elements, ion_end,
     &      ifionized, if_pteh, if_mc, ifreducedmass,
     &      ifsame_abundances, ifmtrace,
     &      iatomic_number, ifpi_local,
     &      ifpl, ifmodified, ifh2, ifh2plus,
     &      izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &      eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &      r_ion3, r_neutral, nelements, 
     &      ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &      rhostar,
     &      ne, nef, net, sion, sionf, siont, uion,
     &      pnorad, pnoradf, pnoradt, pnorad_daux,
     &      free_fp, free_daux,
     &      nextrasum, maxnextrasum,
     &      sumpl1, sumpl1f, sumpl1t, sumpl2,
     &      h2, h2f, h2t, h2_dv)
          if(debug_jacobian) then
            match_variable = match_variable_save
            if(kif.eq.0) fl = match_variable_save
            tl = tl_save
            return
          endif

C        if aux or aux_old is zero (aside from the log density), this
C        is the result of underflow.  this problem should
C        be isolated from the rest of the solution as much as possible
C        by zeroing all rows and columns of the Jacobian that are
C        concerned with the auxiliary variable that is zero while
C        the diagonal element is set to unity.
C        (the assumption here is that one must be careful with auxiliary
C        variables that are the result of an error.  We do the most
C        conservative thing which is a simple iteration.  This should leave
C        the remainder of the solution largely unaffected since *usually*
C        a transition to or from the underflow condition means the associated
C        non-zero variable is small.)
C          n.b. this transformation has already been effectively done
C          so comment out.
!        do jndex_aux = 1, n_partial_aux
!          jaux = partial_aux(jndex_aux)
!          do index_aux = 1, n_partial_aux
!            iaux = partial_aux(index_aux)
!C            zero rows *or* columns....
!            if(
!     &          ((jaux.ne.4.and.
!     &          (aux(jaux).eq.0.d0.or.aux_old(jaux).eq.0.d0)).or.
!     &          (iaux.ne.4.and.
!     &          (aux(iaux).eq.0.d0.or.aux_old(iaux).eq.0.d0)))
!     &          ) then
!              if(iaux.eq.jaux) then
!                jacobian(jndex_aux,index_aux) = 1.d0
!              else
!                jacobian(jndex_aux,index_aux) = 0.d0
!              endif
!            endif
!          enddo
!          if(njacobian.gt.n_partial_aux) then
!            if(jaux.ne.4.and.
!     &          (aux(jaux).eq.0.d0.or.aux_old(jaux).eq.0.d0)) then
!              jacobian(jndex_aux,njacobian) = 0.d0
!              jacobian(njacobian,jndex_aux) = 0.d0
!            endif
!          endif
!        enddo
C          Jacobian(i,j) is negative partial ith equation wrt jth auxiliary
C          variable in (usually) log form.  Decide on whether an auxiliary
C          variable is of major importance based on whether any off-diagonal
C          element of its Jacobian column is greater than major_crit * row
C          norm.
          do jjacobian = 1,njacobian
            ifmajor(jjacobian) = .false.
          enddo
          do ijacobian = 1,njacobian
            rhs_save(ijacobian) = rhs(ijacobian)
            row_norm = 0.d0
            do jjacobian = 1,njacobian
              jacobian_save(ijacobian,jjacobian) =
     &          jacobian(ijacobian,jjacobian)
              row_norm = max(row_norm,
     &          abs(jacobian(ijacobian,jjacobian)))
            enddo
            do jjacobian = 1,njacobian
              if(jjacobian.ne.ijacobian)
     &          ifmajor(jjacobian) = ifmajor(jjacobian).or.
     &          (abs(jacobian(ijacobian,jjacobian)).gt.
     &          major_crit*row_norm)
            enddo
          enddo
          if(.true..and.debug_output) then
            write(0,*) 'njacobian, naux = ', njacobian, naux
            do jndex_aux = 1,njacobian
              if(jndex_aux.le.n_partial_aux) then
                jaux = partial_aux(jndex_aux)
                write(0,*) 'jndex_aux, allow_log, allow_log_neg, '//
     &            'aux_old, aux, rhs ='
                write(0,'(i5,2l5,1p3d20.10)')
     &            jndex_aux, allow_log(jndex_aux),
     &            allow_log_neg(jndex_aux),
     &            aux_old(jaux), aux(jaux),
     &            rhs(jndex_aux)
              else
                write(0,*) 'jndex_aux, rhs ='
                write(0,'(i5,50x,1pd20.10)')
     &            jndex_aux, rhs(jndex_aux)
              endif
              write(0,*) 'jacobian(jndex_aux, index_aux) ='
              write(0,'(5(1pd14.5,2x))')
     &          (jacobian(jndex_aux, index_aux),
     &          index_aux = 1, njacobian)
            enddo
          endif
C          LU factorization solution from lapack.
C          n.b. scaling is done inside if necessary and jacobian and rhs *may*
C          be modified by these scaling factors.
C          the returned sol_lapack is the solution of the original unscaled system
C          of equations.
          call dgesvx('E', 'N',
     &      njacobian, 1, jacobian, naux,
     &      lu_lapack, naux, ipiv_lapack, equed_lapack,
     &      row_lapack, col_lapack, rhs, naux,
     &      sol_lapack, naux,
     &      rcond_lapack, ferr_lapack, berr_lapack,
     &      work_lapack, iwork_lapack, info_lapack)
          if(info_lapack.ne.0) then
C            write(0,*) 'match_variable/ln(10), fl, tl/ln(10) ='
C            write(0,'(1p5d25.15)') match_variable/log(10.d0), fl,
C     &        tl/log(10.d0)
C            write(0,*) 'eps ='
C            write(0,'(1p5d25.15)') eps
C            write(0,*) 'ioncount, ne, ne_old ='
C            write(0,'(i5,1p5d25.15)') ioncount, ne, ne_old
C            write(0,'(a,/,i5,1p5d25.15)')
C     &        'info_lapack, rcond_lapack = ',
C     &        info_lapack, rcond_lapack
C            write(0 ,*) 'free_eos_detailed: (1) no dgesvx solution'
            info = 3
            return
          endif
!          if(debug_output.or.rcond_lapack.lt.1.d-10)
!     &      write(0,'(a,/,i5,1p5d25.15)')
!     &      'ioncount, rcond_lapack = ',
!     &      ioncount, rcond_lapack
          maxfaux = 0.d0
          do index_aux = 1, n_partial_aux
C            find solution
            faux_nr(index_aux) = sol_lapack(index_aux,1)
            iaux = partial_aux(index_aux)
C            update allow_log and allow_log_neg for scaled results.
C            n.b. 4th auxiliary variable (rl = log(rho)) is already in 
C            log form.
C            N.B. "allow_log" connotes "allow taking log of auxiliary variable
C            delivered by eos_jacobian"
            allow_log(index_aux) =
     &        iaux.ne.4.and.
     &        aux_old(iaux).gt.0.d0.and.
     &        aux(iaux).gt.0.d0
C            sign of xextrasum depends on extrasum.  it's okay if negative,
C            but subsequent logic must adjust for this.
C            N.B. "allow_log_neg" connotes "allow taking log of *negative* of
C            auxiliary variable delivered by eos_jacobian"
            allow_log_neg(index_aux) =
     &        iaux.gt.maxnextrasum.and.
     &        aux_old(iaux).lt.0.d0.and.
     &        aux(iaux).lt.0.d0
C            4th auxiliary variable is already in log form so the if
C            statement selects all log variables.
            if((iaux.eq.4.or.allow_log(index_aux).or.
     &          allow_log_neg(index_aux)).and.
     &          abs(faux_nr(index_aux)).gt.abs(maxfaux)) then
              index_max = index_aux
              maxfaux = faux_nr(index_aux)
            endif
          enddo
          if(njacobian.gt.n_partial_aux) then
C            raw change in fl
            pac_nr = sol_lapack(njacobian,1)
            if(abs(pac_nr).gt.abs(maxfaux)) then
              index_max = njacobian
              maxfaux = pac_nr
            endif
          endif
          if(debug_output) then
            if(njacobian.gt.n_partial_aux) then
              write(0,*) 'raw NR solution including fl'
              write(0,'(5(1pd14.5,l2))')
     &          (faux_nr(index_aux), ifmajor(index_aux),
     &          index_aux = 1, n_partial_aux),
     &          pac_nr, ifmajor(njacobian)
            else
              write(0,*) 'raw NR solution with frozen fl'
              write(0,'(5(1pd14.5,l2))')
     &          (faux_nr(index_aux), ifmajor(index_aux),
     &          index_aux = 1, n_partial_aux)
            endif
          endif
          !temporary put no limit on change by disabling following do loop
          !do index_aux = 1, 0
          any_ifsimple = .false.
          do index_aux = 1, n_partial_aux
C            default is the NR solution.
            faux(index_aux) = faux_nr(index_aux)
C            to control changes far from the solution it is usually
C            best to_use_a lower-order solution, e.g., the simple
C            iteration solution where auxnew = result delivered by
C            eos_calc.  We only substitute the lower-order solution if
C            the NR solution is larger than both faux_limit_diag *and* the
C            lower-order solution.
            iaux = partial_aux(index_aux)
            if(allow_log(index_aux).or.
     &          allow_log_neg(index_aux)) then
              if(abs(faux_nr(index_aux)).gt.max(faux_limit_diag,
     &            abs(log(aux(iaux)/aux_old(iaux))))) then
                any_ifsimple = .true.
                simple_lambda(index_aux) = simple_lambda_max
              else
                simple_lambda(index_aux) = simple_lambda_ratio*
     &            simple_lambda(index_aux)
                if(simple_lambda(index_aux).ge.simple_lambda_min.and.
     &              abs(maxfaux).ge.simple_lambda_min) then
                  any_ifsimple = .true.
                else
                  simple_lambda(index_aux) = 0.d0
                endif
              endif
            else
C              no log transformation allowed.
C              iaux = 4 (already logarithmic) or
C              aux or aux_old was zero or had opposite signs.
              if(abs(faux_nr(index_aux)).gt.max(faux_limit_diag,
     &          abs(aux(iaux)-aux_old(iaux)))) then
                any_ifsimple = .true.
                simple_lambda(index_aux) = simple_lambda_max
              else
                simple_lambda(index_aux) = simple_lambda_ratio*
     &            simple_lambda(index_aux)
                if(simple_lambda(index_aux).ge.simple_lambda_min.and.
     &              abs(maxfaux).ge.simple_lambda_min) then
                  any_ifsimple = .true.
                else
                  simple_lambda(index_aux) = 0.d0
                endif
              endif
            endif
          enddo
          if(njacobian.gt.n_partial_aux) then
C            default value.
            pac = pac_nr
            simple_lambda(njacobian) = 0.d0
          endif
          if(any_ifsimple) then
C            multiply RHS and diagonal by (1+lambda).  For large lambda this
C            has similar effect to zeroing off-diagonal elements.
            do ijacobian = 1,njacobian
              rhs(ijacobian) = (1.d0 + simple_lambda(ijacobian))*
     &          rhs_save(ijacobian)
              do jjacobian = 1, njacobian
                jacobian(ijacobian,jjacobian) =
     &            jacobian_save(ijacobian,jjacobian)
              enddo
              jacobian(ijacobian,ijacobian) =
     &          jacobian(ijacobian,ijacobian) + 
     &          simple_lambda(ijacobian)
            enddo
C            LU factorization solution from lapack.
C            n.b. scaling is done inside if necessary and jacobian and rhs
C            *may* be modified by these scaling factors.
C            the returned sol_lapack is the solution of the original
C            unscaled system of equations.
            call dgesvx('E', 'N',
     &        njacobian, 1, jacobian, naux,
     &        lu_lapack, naux, ipiv_lapack, equed_lapack,
     &        row_lapack, col_lapack, rhs, naux,
     &        sol_lapack, naux,
     &        rcond_lapack, ferr_lapack, berr_lapack,
     &        work_lapack, iwork_lapack, info_lapack)
            if(info_lapack.ne.0) then
C              write(0,*) 'match_variable/ln(10), fl, tl/ln(10) ='
C              write(0,'(1p5d25.15)') match_variable/log(10.d0), fl,
C     &          tl/log(10.d0)
C              write(0,*) 'eps ='
C              write(0,'(1p5d25.15)') eps
C              write(0,*) 'ioncount, ne, ne_old ='
C              write(0,'(i5,1p5d25.15)') ioncount, ne, ne_old
C              write(0,'(a,/,i5,1p5d25.15)')
C     &          'info_lapack, rcond_lapack = ',
C     &          info_lapack, rcond_lapack
C              write(0,*) 'free_eos_detailed: (2) no dgesvx solution'
              info = 4
              return
            endif
!            if(debug_output.or.rcond_lapack.lt.1.d-10)
!     &        write(0,'(a,/,i5,1p5d25.15)')
!     &        'ioncount, rcond_lapack = ',
!     &        ioncount, rcond_lapack
            do ijacobian = 1, njacobian
              if(ijacobian.gt.n_partial_aux) then
                pac = sol_lapack(ijacobian,1)
              else
                faux(ijacobian) = sol_lapack(ijacobian,1)
              endif
            enddo
          endif
          if(debug_output) then
            write(0,*) 'simple_lambda'
            write(0,'(5(1pd14.5,l2))')
     &        (simple_lambda(index_aux),
     &        simple_lambda(index_aux).gt.0.d0,
     &        index_aux = 1, njacobian)
            if(njacobian.gt.n_partial_aux) then
              write(0,*)
     &          'Mixed NR and simple-iteration solution including fl'
              write(0,'(5(1pd14.5,l2))')
     &          (faux(index_aux),
     &          simple_lambda(index_aux).gt.0.d0,
     &          index_aux = 1, n_partial_aux),
     &          pac, simple_lambda(njacobian).gt.0.d0
            else
              write(0,*) 'Mixed NR and simple-iteration solution '//
     &          'with frozen fl'
              write(0,'(5(1pd14.5,l2))')
     &          (faux(index_aux),
     &          simple_lambda(index_aux).gt.0.d0,
     &          index_aux = 1, n_partial_aux)
            endif
          endif
          maxfaux_diag = 0.d0
          do index_aux = 1, n_partial_aux
C            4th auxiliary variable is already in log form so the if
C            statement selects all log variables and only if major
C            auxiliary variable.
            iaux = partial_aux(index_aux)
            if(ifmajor(index_aux).and.
     &          (iaux.eq.4.or.allow_log(index_aux).or.
     &          allow_log_neg(index_aux)).and.
     &          abs(faux(index_aux)).gt.abs(maxfaux_diag)) then
C              n.b. maxfaux_diag perturbed by diagonal component and
C              used for scaling, but maxfaux used for iteration control
              maxfaux_diag = faux(index_aux)
            endif
          enddo
          if(njacobian.gt.n_partial_aux.and.
     &      ifmajor(njacobian).and.
     &      abs(pac).gt.abs(maxfaux_diag))
     &      maxfaux_diag = pac

C          if change in sign, halve faux_limit_small (unless would force
C          below faux_limit_min)
          if(maxfaux_diag*maxfaux_diag_old.lt.0.d0.and.
     &      ioncount.ge.2)
     &      faux_limit_small =
     &      max(faux_limit_min,0.5d0*faux_limit_small)
            
C          only worry about maximum changes if there would have been
C          a substantial change in ne proportional to n_e/rho for the
C          previous *unscaled* solution.  ne is an overall
C          measure of ionization fractions.  (We also control rho for those
C          cases [molecular formation important] where it is more independent
C          of ne.) Also, lambda.  This logic is meant to 
C          take care of the case where there are large ln changes in an
C          auxiliary variable that have little or no effect on overall
C          ionization fractions because the auxiliary variable is
C          approaching zero.
          if(max(abs(rl-rl_old),
     &        abs(ne-ne_old)/ne,
     &        abs(lambda-lambda_old)/max(1.d-15,lambda))/
     &        faux_scale.le.1.d-2) then
            faux_limit = faux_limit_large
            faux_limit_small = faux_limit_small_start
          else
            faux_limit = faux_limit_small
          endif
          faux_scale = min(1.d0,
     &      faux_limit/max(1.d-16,abs(maxfaux_diag)))
          if(njacobian.gt.n_partial_aux) then
C            scale so that maximum change in fl is always <= 0.01
            faux_scale = min(faux_scale,
     &        min(faux_limit,0.01d0)/max(1.d-16,abs(pac)))
          endif
C          cannot get into much trouble if current fl and NR-predicted fl
C          are less than -10 so turn off limits to NR change in that case.
          if(njacobian.gt.n_partial_aux) then
            if(max(fl,fl+pac).lt.-10.d0) faux_scale = 1.d0
          else
C            fixed fl.
            if(fl.lt.-10.d0) faux_scale = 1.d0
          endif
          !temporary (no limit on NR change)
          !faux_scale = 1.d0
          do index_aux = 1, n_partial_aux
            iaux = partial_aux(index_aux)
C            do nothing when eos_calc returns a zero aux (index.ne.4)
C            because solution may add some noise.
            if(allow_log(index_aux).or.
     &          allow_log_neg(index_aux).or.
     &          iaux.eq.4) then
C              these comments pertain to each of branches below:
C              keep direction of faux vector, but scale it so no
C              log component exceeds faux_limit.
              faux(index_aux) = faux_scale*faux(index_aux)
              if(allow_log(index_aux)) then
                aux(iaux) = log(aux_old(iaux)) + faux(index_aux)
C                underflow trap
!                if(aux(iaux).gt.-575.d0) then
                  aux(iaux) = exp(aux(iaux))
!                else
!                  aux(iaux) = 0.d0
!                endif
              elseif(allow_log_neg(index_aux)) then
                aux(iaux) = log(-aux_old(iaux)) + faux(index_aux)
C                underflow trap
!                if(aux(iaux).gt.-575.d0) then
                  aux(iaux) = -exp(aux(iaux))
!                else
!                  aux(iaux) = 0.d0
!                endif
              elseif(iaux.eq.4) then
                aux(iaux) = aux_old(iaux) + faux(index_aux)
C                n.b. last two branches are old code that is avoided
C                by above outer if statement.
              elseif(sign(1.d0,aux(iaux)).eq.
     &            sign(1.d0, aux_old(iaux) + faux(index_aux))) then
C                by logic of code this branch occurs only if sign change
C                or change from zero to non-zero
C                only_use_Newton-Raphson value if its new sign agrees with
C                sign of value delivered by eos_calc.
                aux(iaux) = aux_old(iaux) + faux(index_aux)
              elseif(aux_old(iaux).ne.0.d0) then
C                scale change which involves sign change that disagrees with
C                NR sign
                aux(iaux) = aux_old(iaux) +
     &            faux_scale*(aux(iaux)-aux_old(iaux))
              endif
            endif
C            n.b. above logic drops through (accepts aux delivered by
C            eos_calc) in case aux_old or aux is zero or
C            opposite signs (aside from iaux.eq.4 which is logarithmic).

C            overall underflow trap (exp(-575) = 10^{-249.7})
!            if(abs(aux(iaux)).lt.1.d-250)
!     &        aux(iaux) = 0.d0
C            if trying to go from non-zero to zero or vice versa
C           _use_special limits on aux
            if(iaux.ne.4.and.
     &          (aux_old(iaux).ne.aux(iaux).and.
     &          (aux_old(iaux).eq.0.d0.or.
     &          aux(iaux).eq.0.d0))) then
C            continue with iteration unless zerolim is zero
              if(zerolim(iaux).gt.0.d0) iflast = 0
              if(aux_old(iaux).ne.0.d0) then
C                if trying to go from non-zero to zero...
                aux(iaux) = aux_old(iaux)/exp(faux_limit)
!                if(abs(aux(iaux)).lt.1.d-250) aux(iaux) = 0.d0
C                define maximum magnitude for subsequent zero to non-zero
                zerolim(iaux) = abs(aux(iaux))
C                mark infinite decrease in log auxiliary variable.
                maxfaux = -50.d0
                index_max = index_aux
C                if this occurred on previous iteration, increase
C                ifzerocount and allow to go to zero if 5 in a row
                if(ioncount-1.eq.ioncountzero) then
C                  zero occurred on previous iteration
C                  n.b. this branch only once per iteration
                  ifzerocount = ifzerocount + 1
                elseif(ioncount.ne.ioncountzero) then
C                  zero did not occur for previous iteration.
C                  n.b. this branch only once per iteration
                  ifzerocount = 1
                endif
                ioncountzero = ioncount
                if(mod(ifzerocount,5).eq.0) aux(iaux) = 0.d0
              else
C                if going from zero to non-zero, don't go above
C                previous magnitude found on non-zero to zero step
C                or if no previous such step (i.e., initial parameter was
C                zero), then only allow magnitude of 1.d-240 (see
C                initialization of zerolim above).
                if(abs(aux(iaux)).gt.zerolim(iaux))
     &            aux(iaux) = sign(zerolim(iaux), aux(iaux))
C                mark infinite increase in log auxiliary variable.
                maxfaux = 50.d0
                index_max = index_aux
              endif
            endif
          enddo
          if(njacobian.gt.n_partial_aux) then
C            scale fl change just like raw NR changes were
C            scaled above.
            pab = faux_scale*pac
            fl = fl + pab
          endif
C          rho must be updated also.
          rho = exp(rl)
          if(debug_output) then
            if(njacobian.gt.n_partial_aux) then
              write(0,*)
     &          'old value of aux including fl'
              write(0,'(5(1pd14.5,l2))')
     &          (aux_old(partial_aux(index_aux)),
     &          ifmajor(index_aux),
     &          index_aux = 1, n_partial_aux),
     &          flold, ifmajor(njacobian)
              write(0,*)
     &          'new value of aux including fl'
              write(0,'(5(1pd14.5,l2))')
     &          (aux(partial_aux(index_aux)),
     &          ifmajor(index_aux),
     &          index_aux = 1, n_partial_aux),
     &          fl, ifmajor(njacobian)
            else
              write(0,*)
     &          'old value of aux for frozen fl'
              write(0,'(5(1pd14.5,l2))')
     &          (aux_old(partial_aux(index_aux)),
     &          ifmajor(index_aux),
     &          index_aux = 1, n_partial_aux)
              write(0,*)
     &          'new value of aux for frozen fl'
              write(0,'(5(1pd14.5,l2))')
     &          (aux(partial_aux(index_aux)),
     &          ifmajor(index_aux),
     &          index_aux = 1, n_partial_aux)
            endif
          endif
          if(ioncount.ge.100.or.debug_output) then
C            write(0,*) 'ioncount, index_max, maxfaux, '//
C     &        'maxfaux_diag, faux_limit, faux_scale'
C            write(0,*) ioncount, index_max, maxfaux, maxfaux_diag,
C     &        faux_limit, faux_scale
          endif
        enddo  !end of iteration for consistent auxiliary variables and fl
        if((ifwarm.and.ioncount.ge.maxioncount).or.
     &      ioncount.ge.maxioncount) then
C          write(0,*) 'match_variable/ln(10), fl, tl/ln(10) ='
C          write(0,'(1p5d25.15)') match_variable/log(10.d0), fl,
C     &      tl/log(10.d0)
C          write(0,*) 'eps ='
C          write(0,'(1p5d25.15)') eps
C          write(0,*) 'ioncount, ne, ne_old ='
C          write(0,'(i5,1p5d25.15)') ioncount, ne, ne_old
C          write(0,*) 'maxfaux ='
C          write(0,'(1pd25.15)') maxfaux
C          write(0,*) 'auxiliary variable iteration did not converge'
          info = 5
          return
        endif
        !temporary just to get comparisons, but messes up derivatives.
        if(.false.) then
          write(0,*) 'start final eos_bfgs per NR cycle'
          call eos_bfgs(
     &      allow_log, allow_log_neg, partial_aux,
     &      ifrad, debug_dauxdv, itemp1, itemp2, delta1, delta2,
     &      debug_dvdaux, debug_jacobian,
     &      degeneracy, pressure, density, energy, enthalpy, entropy,
     &      iaux1, iaux2,
     &      match_variable, match_variable_save, kif, fl, tl_save,
     &      aux_old, aux, auxf, auxt, daux_dv,
     &      njacobian, njacobian_final,
     &      rhs, jacobian, p, pr,
     &      ddv_aux, sumion0, sumion2,
     &      lambda, gamma_e,
     &      ifnr, nux, nuy, nuz, mion_end,
     &      n_partial_ions, n_partial_aux,
     &      sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &      partial_ions, t, morder, ifexchange_in,
     &      full_sum0, full_sum1, full_sum2, charge, charge2,
     &      dv_pl, dv_plt,
     &      ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &      ifsame_zero_abundances,
     &      ifexcited, ifsame_under,
     &      naux, inv_aux,
     &      inv_ion, max_index,
     &      partial_elements, n_partial_elements, ion_end,
     &      ifionized, if_pteh, if_mc, ifreducedmass,
     &      ifsame_abundances, ifmtrace,
     &      iatomic_number, ifpi_local,
     &      ifpl, ifmodified, ifh2, ifh2plus,
     &      izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &      eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &      r_ion3, r_neutral, nelements, 
     &      ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &      ne, nef, net, sion, sionf, siont, uion,
     &      pnorad, pnoradf, pnoradt, pnorad_daux,
     &      free_fp, free_daux,
     &      nextrasum, maxnextrasum,
     &      sumpl1, sumpl1f, sumpl1t, sumpl2,
     &      h2, h2f, h2t, h2_dv, info)
          if(info.ne.0) return
          write(0,*) 'end final eos_bfgs per NR cycle'
        endif
C        save aux so that can restore after eos_calc below.
C        (It is possible that eos_calc could move off of
C        NR-converged value by propagating numerical
C        noise.)
        do index_aux = 1,n_partial_aux
          iaux = partial_aux(index_aux)
          aux_old(iaux) = aux(iaux)
        enddo
C        from partial of RHS wrt match_variable and tl_use_chain rule
C        and pre-factored Jacobian to calculate partials of 
C        converged auxiliary variables and fl wrt match_variable and tl.
        do index_aux = 1, n_partial_aux
          iaux = partial_aux(index_aux)
          if(allow_log(index_aux).or.
     &        allow_log_neg(index_aux)) then
C            transform to log form of faux
            if(njacobian.gt.n_partial_aux) then
              rhs_lapack(index_aux,1) = 0.d0
            else
              rhs_lapack(index_aux,1) = auxf(iaux)/aux_old(iaux)
            endif
            rhs_lapack(index_aux,2) = auxt(iaux)/aux_old(iaux)
          else
            if(njacobian.gt.n_partial_aux) then
              rhs_lapack(index_aux,1) = 0.d0
            else
              rhs_lapack(index_aux,1) = auxf(iaux)
            endif
            rhs_lapack(index_aux,2) = auxt(iaux)
          endif
        enddo
        if(njacobian.gt.n_partial_aux) then
          if(kif.eq.1) then
            if(ifrad.gt.1) then
C              For kif = 1, RHS = match_variable - log(pnorad)
C              partial of RHS wrt tl
              rhs_lapack(index_aux,2) = -pnoradt/pnorad
            else
C              For kif = 1, RHS = match_variable - log(p)
              rhs_lapack(index_aux,2) = -(4.d0*pr + pnoradt)/p
            endif
          elseif(kif.eq.2) then
C            For kif = 2, RHS = match_variable - rl
            rhs_lapack(index_aux,2) = -rt
          endif
C          partial of RHS wrt to match_variable.
          rhs_lapack(index_aux,1) = 1.d0
        endif
C        LU factorization solution from lapack with jacobian already
C        scaled (potentially) and factored to lu_lapack
C        n.b. equed_lapack keeps track of what scaling occurred on
C        previous call that factored (scaled) jacobian and this
C        scaling if any is applied to rhs_lapack
C        the returned sol_lapack is the solution of the original
C        unscaled system of equations.
        call dgesvx('F', 'N',
     &    njacobian, 2, jacobian, naux,
     &    lu_lapack, naux, ipiv_lapack, equed_lapack,
     &    row_lapack, col_lapack, rhs_lapack, naux,
     &    sol_lapack, naux,
     &    rcond_lapack, ferr_lapack, berr_lapack,
     &    work_lapack, iwork_lapack, info_lapack)
        if(info_lapack.ne.0.or.rcond_lapack.lt.1.d-10) then
C          write(0,*) 'match_variable/ln(10), fl, tl/ln(10) ='
C          write(0,'(1p5d25.15)') match_variable/log(10.d0), fl,
C     &      tl/log(10.d0)
C          write(0,*) 'eps ='
C          write(0,'(1p5d25.15)') eps
C          write(0,*) 'ioncount, ne, ne_old ='
C          write(0,'(i5,1p5d25.15)') ioncount, ne, ne_old
C          write(0,'(a,/,i5,1p5d25.15)')
C     &      'info_lapack, rcond_lapack = ',
C     &      info_lapack, rcond_lapack
C          write(0,*) 'free_eos_detailed: (3) no dgesvx solution'
          info = 6
          return
        endif
!        if(debug_output.or.rcond_lapack.lt.1.d-10)
!     &    write(0,'(a,/,i5,1p5d25.15)')
!     &    'ioncount, rcond_lapack = ',
!     &    ioncount, rcond_lapack
        do index_aux = 1, n_partial_aux
          iaux = partial_aux(index_aux)
          if(allow_log(index_aux).or.
     &        allow_log_neg(index_aux)) then
C            transform from log auxiliary variable derivatives
            auxf(iaux) = sol_lapack(index_aux,1)*aux_old(iaux)
            auxt(iaux) = sol_lapack(index_aux,2)*aux_old(iaux)
          else
            auxf(iaux) = sol_lapack(index_aux,1)
            auxt(iaux) = sol_lapack(index_aux,2)
          endif
        enddo
        if(njacobian.gt.n_partial_aux) then
C          transform independent variable from match_variable to fl.
C          this could introduce significance loss when we transform
C          back again below, but this does not appear to be a problem.
C          sol_lapack(njacobian,1) is the partial of
C          fl(match_variable,tl) wrt match_variable
C          sol_lapack(njacobian,2) is the partial of
C          fl(match_variable,tl) wrt tl.
C          partial of calculated match_variable(fl,tl) wrt fl
          match_variablef = 1.d0/sol_lapack(njacobian,1)
C          partial of calculated match_variable(fl,tl) wrt tl
          match_variablet =
     &      -match_variablef*sol_lapack(njacobian,2)
          do index_aux = 1, n_partial_aux
            iaux = partial_aux(index_aux)
C            auxf and auxt contain partial of aux wrt match_variable and tl.
C            so transform appropriately to fl and tl derivatives.
            auxt(iaux) = auxf(iaux)*match_variablet + auxt(iaux)
            auxf(iaux) = auxf(iaux)*match_variablef
          enddo
        endif
        if(debug_output) then
          write(0,*) 'auxf = '
          write(0,'(5(1pd14.5,2x))')
     &      (auxf(partial_aux(index_aux)),
     &      index_aux = 1, n_partial_aux)
          write(0,*) 'auxt = '
          write(0,'(5(1pd14.5,2x))')
     &      (auxt(partial_aux(index_aux)),
     &      index_aux = 1, n_partial_aux)
        endif
C*************************************************************
C        now that all auxiliary variable derivatives have
C        been determined from NR procedure, call whole
C        eos procedure again to determine fl and tl derivatives
C        of all thermodynamic quantities (tq).
C        n.b.  eos_tqft returns nu form of many auxiliary variables
C        since this form is required by the entropy calculation.
        iteration_count = iteration_count + 1
        call eos_tqft(
     &    ddv_aux, sumion0, sumion2,
     &    lambda, gamma_e,
     &    nux, nuy, nuz, mion_end,
     &    n_partial_ions, n_partial_aux,
     &    sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &    partial_ions, f, eta, wf, t, n_e, pstar, 
     &    dve_exchange, dve_exchangef, dve_exchanget,
     &    full_sum0, full_sum1, full_sum2, charge, charge2,
     &    dv_pl, dv_plt,
     &    ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &    ifsame_zero_abundances, ifexcited,
     &    naux, inv_aux,
     &    inv_ion, max_index,
     &    partial_elements, n_partial_elements, ion_end,
     &    ifionized, if_pteh, if_mc, ifreducedmass,
     &    ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,
     &    ifpl, ifmodified, ifh2, ifh2plus,
     &    izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &    eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,
     &    r_ion3, r_neutral, nelements, 
     &    ifelement, dvzero, dv, dvf, dvt, nion, nions,
     &    rhostar,
     &    ne, nef, net, sion, sionf, siont, uion,
     &    h, hf, ht, hd, hdf, hdt, he, hef, het,
     &    h_dv, hd_dv, he_dv,
     &    sum0, sum0f, sum0t, sum0_dv,
     &    sum2, sum2f, sum2t, sum2_dv,
     &    extrasum, extrasumf, extrasumt, extrasum_dv,
     &    nextrasum, maxnextrasum,
     &    sumpl1, sumpl1f, sumpl1t, sumpl2,
     &    rl, rf, rt, r_dv,
     &    h2, h2f, h2t, h2plus, h2plusf, h2plust, h2plus_dv,
     &    xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
C        do not restore aux, since part of it has been put
C        in nu form by eos_calc.
C*************************************************************
C        end of if block to calculate derivatives of all
C        thermodynamic quantities wrt fl and tl using the chain rule
C        for NR-converged results.
      endif !if(.not.ifnoaux_iteration)
C      there is the possibility of pressure ionization and 
C      Planck-Larkin terms.
C      n.b. For kif = 1, some of the pressure-related stuff (but not all of
C      it in the ifnoaux_iteration case) is calculated previously.  But just
C      repeat those calculations here rather than worrying about special logic.
      if(ifpi_local.eq.0) then
C        there are no p, s, u terms from pressure ionization.
        ppi = 0.d0
        ppif = 0.d0
        ppit = 0.d0
        spi = 0.d0
        spif = 0.d0
        spit = 0.d0
        upi = 0.d0
      elseif(ifpi_local.eq.1) then
        call pteh_pi_end(
     &    full_sum1, rho, rf, rt, t, ne, nef, net,
     &    ppi, ppif, ppit, spi, spif, spit, upi)
      elseif(ifpi_local.eq.2) then
C        n.b. call to eos_tqft produces nu = n/(rho*avogadro)
C        form of h, hd, and he.  Also_use_nu form
C        of nx, ny, nz, ne, etc.
        call fjs_pi_end(
     &    t, rho, rf, rt,
     &    nux, 0.d0, 0.d0, nuy, 0.d0, 0.d0, nuz, 0.d0, 0.d0,
     &    h, hf, ht, hd, hdf, hdt, he, hef, het, ne, nef, net,
     &    ppi, ppif, ppit, spi, spif, spit, upi)
      elseif(ifpi_local.eq.3.or.ifpi_local.eq.4) then
C        n.b. call to eos_tqft produces nu = n/(rho*avogadro)
C        form of extrasum.
        call mdh_pi_end(
     &    t, rho, rf, rt,
     &    r_ion3, nion, nions, r_neutral, nelements,
     &    extrasum, extrasumf, extrasumt, nextrasum,
     &    ppi, ppif, ppit, spi, spif, spit, upi)
      endif
      if(ifexcited.gt.0) then
        call excitation_pi_end(t, rho, rf, rt,
     &    pexcited, pexcitedf, pexcitedt,
     &    sexcited, sexcitedf, sexcitedt, uexcited)
      else
        pexcited = 0.d0
        pexcitedf = 0.d0
        pexcitedt = 0.d0
        sexcited = 0.d0
        sexcitedf = 0.d0
        sexcitedt = 0.d0
        uexcited = 0.d0
      endif
      if(ifpl.eq.1) then
C        no pressure terms for Planck-Larkin occupation probability
C        but there are entropy and energy terms.
        spi = spi + cr*sumpl1
        spif = spif + cr*sumpl1f
        spit = spit + cr*sumpl1t
        upi = upi + cr*t*sumpl2
      endif
      if(if_pteh.eq.1) then
C        PTEH (full ionization) approximation to sum0, sum2
C        reassert this approximation because eos_calc
C        calls ionize which messes a bit with sum0 and sum2
        sum0ne = full_sum0/full_sum1
        sum0 = sum0ne*n_e
        sum0f = 0.d0
        sum0t = 0.d0
        sum2ne = full_sum2/full_sum1
        sum2 = sum2ne*n_e
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
C        sumion0 = sum0_mc*rho*avogadro
C        sumion2 = sum2_mc*rho*avogadro
C        n.b. call to eos_tqft produces nu = n/(rho*avogadro)
C        form of h2plus, h, hd, he.
        sum0a = (h2plus + h*(1.d0+hcon_mc) +
     &    hd +  he)*(rho*avogadro)
        sum2a = (h2plus + h*(1.d0+hcon_mc) +
     &    hd +  he*(4.d0+hecon_mc))*(rho*avogadro)
        sum0 = sumion0 + sum0a
        sum2 = sumion2 + sum2a 
        sum0f = sumion0*rf + (h2plusf + hf*(1.d0+hcon_mc) +
     &    hdf + hef)*(rho*avogadro) + sum0a*rf
        sum0t = sumion0*rt + (h2plust + ht*(1.d0+hcon_mc) +
     &    hdt + het)*(rho*avogadro) + sum0a*rt
        sum2f = sumion2*rf + (h2plusf + hf*(1.d0+hcon_mc) +
     &    hdf + hef*(4.d0+hecon_mc))*(rho*avogadro) + sum2a*rf
        sum2t = sumion2*rt + (h2plust + ht*(1.d0+hcon_mc) +
     &    hdt + het*(4.d0+hecon_mc))*(rho*avogadro) + sum2a*rt
      endif
      call master_coulomb_end(rhostar,
     &  sum0, sum0ne, sum0f, sum0t,
     &  sum2, sum2ne, sum2f, sum2t,
     &  n_e, t, pstar,
     &  ifcoulomb_mod, if_dc, if_pteh,
     &  dpcoulomb, dpcoulombf, dpcoulombt,
     &  dscoulomb, dscoulombf, dscoulombt, ducoulomb)
      call exchange_end(rhostar, pstar, 3,
     &  pex, pext, pexf,
     &  sex, sexf, sext, uex)
      p0 = cr*rho*t
      ni = full_sum0 - (h2+h2plus)
      pion = ni*p0
      pionf = ni*p0*rf - (h2f+h2plusf)*p0
      pe_cgs = cpe*pe
      pnorad = pe_cgs + pion + ppi + pexcited + dpcoulomb + pex
      p = pnorad + pr
      if(p.le.0.d0) then
C        write(0,*) 'match_variable/ln(10), fl, tl/ln(10) ='
C        write(0,'(1p5d25.15)') match_variable/log(10.d0), fl,
C     &    tl/log(10.d0)
C        write(0,*) 'eps ='
C        write(0,'(1p5d25.15)') eps
C        write(0,*) 'ioncount, ne, ne_old ='
C        write(0,'(i5,1p5d25.15)') ioncount, ne, ne_old
C        write(0,*) 'free_eos_detailed: (2) negative p calculated'
        info = 7
        return
      endif
C      d ln p(f,t)/d ln f
      pf = (pe_cgs*pef + pionf + ppif + pexcitedf +
     &  dpcoulombf + pexf)/p
      pl = log(p)
C      pr_ratio, ppie_ratio, pc_ratio, and pex_ratio  calculated only for
C      diagnostic and teaching purposes.
      pr_ratio = pr/p
C      Combine pressure-ionization correction to lower states and entire
C      effect of pressure-ionization corrected Rydberg states since the
C      two effects tend to offset each other and should be treated together
C      (or not).
      ppie_ratio = (ppi + pexcited)/p
      pc_ratio = dpcoulomb/p
      pex_ratio = pex/p
      
C      pion = ni*p0
      piont = ni*p0*(1.d0 + rt) - (h2t+h2plust)*p0
C      d ln p(f,T)/d ln t
      pnoradt = (pe_cgs*pet + piont + ppit + pexcitedt +
     &  dpcoulombt + pext)/pnorad
      pt = (pnorad*pnoradt + 4.d0*pr)/p
C      rpt = d ln rho(T,P)/ d ln P = 1/chi_rho
      rpt = rf/pf
C      rtp = - d ln rho(T,P)/ d ln T  n.b. *negative* sign = chi_t/chi_rho
      rtp = rpt*pt - rt
      rlout = rl
      degeneracy(1) = fl
      if(kif.eq.0) then
        degeneracy(2) = 1.d0
        degeneracy(3) = 0.d0
      elseif(kif.eq.1) then
        fm = 1.d0/pf
        ft = -pt/pf
        degeneracy(2) = fm
        degeneracy(3) = ft
        if(ifrad.eq.2) then
C          kif=1, and ifrad = 2 returns most derivatives wrt ln p, but in the
C          special case of fm and ft return derivatives wrt ln pnorad
C          since fm and ft used outside in calling free_eos programme
C          for Taylor series.
          fm = 1.d0/pnoradf
          ft = -pnoradt/pnoradf
        endif
      elseif(kif.eq.2) then
        fm = 1.d0/rf
        ft = -rt/rf
        degeneracy(2) = fm
        degeneracy(3) = ft
      endif
      pressure(1) = pl
      if(kif.eq.0) then
        pressure(2) = pf
        pressure(3) = pt
C        special tests
C        pressure(1) = pe_cgs
C        pressure(2) = pef*pe_cgs
C        pressure(3) = pet*pe_cgs
C        pressure(1) = pion
C        pressure(2) = pionf
C        pressure(3) = piont
C        pressure(1) = pexcited
C        pressure(2) = pexcitedf
C        pressure(3) = pexcitedt
C        index_aux = 12
C        iaux = partial_aux(index_aux)
C        pressure(1) = aux(iaux)
C        pressure(2) = auxf(iaux)
C        pressure(3) = auxt(iaux)
C        pressure(1) = dpcoulomb
C        pressure(2) = dpcoulombf
C        pressure(3) = dpcoulombt
C        pressure(1) = pex
C        pressure(2) = pexf
C        pressure(3) = pext
      elseif(kif.eq.1) then
        pressure(2) = 1.d0
        pressure(3) = 0.d0
      elseif(kif.eq.2) then
        if(ifrad.gt.1) then
C          n.b. this makes derivatives of this quantity 
C          inconsistent, but nevertheless some tables have
C          all quantities calculated with radiation
C          pressure except for pg = p - pr.
          pressure(1) = log(pnorad)
        else
          pressure(1) = log(p)
        endif
        pressure(2) = 1.d0/rpt
        pressure(3) = rtp/rpt
      endif
      density(1) = rl
      if(kif.eq.0) then
        density(2) = rf
        density(3) = rt
C        temporary check
C        density(1) = ne
C        density(2) = nef
C        density(3) = net
      elseif(kif.eq.1) then
        density(2) = rpt
        density(3) = -rtp
      elseif(kif.eq.2) then
        density(2) = 1.d0
        density(3) = 0.d0
      endif
      sr = 4.d0*pr/p0
C      the third term corrects for the dropped ideal terms in
C      ionize and eos_calc.  (See the notes for those routines).
      s = ne*se + sr + full_sum0*(2.5d0 + 1.5d0*tl - rl) +  sion
      s = cr*s + spi + sexcited + dscoulomb/rho + sex/rho
      sf = nef*se + ne*se*sef - rf*sr - ni*rf + sionf
      sf = cr*sf + spif + sexcitedf +
     &  (dscoulombf + sexf - (dscoulomb+sex)*rf)/rho
      st = net*se + ne*se*set + (3.d0-rt)*sr + ni*(1.5d0-rt) + siont
      st = cr*st + spit + sexcitedt +
     &  (dscoulombt + sext - (dscoulomb+sex)*rt)/rho
C      by definition of d S(T,P)/d ln T and by chain rule
      cp = st - sf*(pt/pf)
      entropy(1) = s
      if(kif.eq.0) then
        if(iftc.eq.1) then
C          this derivatives from thermodynamic consistency arguments.
C          partial s(p,t) wrt ln p = (p/(rho*t)*d ln rho(P,T)/d ln T
          entropy(2) = -pf*(p/(rho*t))*rtp
        else
C          straight derivative by chain rule
          entropy(2) = sf
        endif
        entropy(3) = st
C        special tests
C        entropy(1) = cr*ne*se
C        entropy(2) = cr*(nef*se + ne*se*sef)
C        entropy(3) = cr*(net*se + ne*se*set)
C        ideal part of entropy
C        entropy(1) = cr*(eps(19)*(2.5d0 + 1.5d0*tl - rl) +  sion)
C        entropy(2) = cr*(- eps(19)*rf + sionf)
C        entropy(3) = cr*(eps(19)*(1.5d0-rt) + siont)
C        entropy(1) = cr*(full_sum0*(2.5d0 + 1.5d0*tl - rl) +  sion)
C        entropy(2) = cr*(-ni*rf + sionf)
C        entropy(3) = cr*(ni*(1.5d0-rt) + siont)
C        entropy(1) = sexcited
C        entropy(2) = sexcitedf
C        entropy(3) = sexcitedt
C        entropy(1) = dv(3) + dvzero(2)
C        entropy(2) = dvf(3)
C        entropy(3) = dvt(3)
C        index_aux = 13
C        iaux = partial_aux(index_aux)
C        entropy(1) = aux(iaux)
C        entropy(2) = auxf(iaux)
C        entropy(3) = auxt(iaux)
C        entropy(1) = dscoulomb/rho
C        entropy(2) = dscoulombf/rho - dscoulomb*rf/rho
C        entropy(3) = dscoulombt/rho - dscoulomb*rt/rho
C        entropy(1) = sex/rho
C        entropy(2) = sexf/rho - sex*rf/rho
C        entropy(3) = sext/rho - sex*rt/rho
      elseif(kif.eq.1) then
        if(iftc.eq.1) then
C          this derivatives from thermodynamic consistency arguments.
C          partial s(p,t) wrt ln p = (p/(rho*t)*d ln rho(P,T)/d ln T
          entropy(2) = -(p/(rho*t))*rtp
        else
C          straight derivative by chain rule
          entropy(2) = sf/pf
        endif
C        partial s(p,t) wrt ln t = C_p by definition (see above).
        entropy(3) = cp
      elseif(kif.eq.2) then
        if(iftc.eq.1) then
C          this derivatives from thermodynamic consistency arguments.
C          partial s(rho,t) wrt ln rho 
          entropy(2) = -(p/(rho*t))*rtp/rpt
        else
C          straight derivative by chain rule
          entropy(2) = sf/rf
        endif
C        straight derivatives by chain rule
C        partial s(rho,t) wrt ln t 
        entropy(3) = st - sf*rt/rf
C        this derivative from thermodynamic consistency arguments.
C        however, it is straightforward to derive this from straight
C        derivatives and thermodynamic consistency for entropy(2),
C        so it is not an independent test.
C        entropy(3) = cp - rtp*rtp*p/rpt/rho/t
      endif
C      vdb expression with metal, h2, and h2+, and coulomb contribution added.
C      note energy of electrons/volume = n_e k T ue, but n_e = ne*rho/H, thus,
C      energy of electrons/volume = rho * ne R T ue, and
C      energy/mass = ne R T ue.
      u = cr*t*(ue*ne+1.5d0*ni+0.75d0*sr) + (c2*cr)*uion +
     &  upi + uexcited + (ducoulomb+uex)/rho
C      for hydrogen shift energy zero from monatomic to diatomic ground state,
C      that is add h2diss/2 to monatomics and h2diss to diatomics
C      n.b. this zero point shift does not affect pressure, entropy, or
C      chemical potentials, but does keep ideal u positive.
      u = u + 0.5d0*(c2*cr)*h2diss*eps(1)
      !temporary
      free_rad = -cr*t*sr/4.d0
      free_e = cr*t*(ue*ne-ne*se)
      free_coulomb = (ducoulomb -t*dscoulomb)/rho
      free_ex = (uex -t*sex)/rho
      free_pl = cr*t*(sumpl2 - sumpl1)
      free_pi = upi -t*spi - free_pl
      free_excited = uexcited -t*sexcited
      free_ion = cr*t*(1.5d0*ni - full_sum0*
     &  (2.5d0 + 1.5d0*tl - rl) - sion) + (c2*cr)*uion +
     &  0.5d0*(c2*cr)*h2diss*eps(1)
      free =
     &  free_rad + free_e + free_coulomb + free_ex +
     &  free_pi +
     &  free_excited + free_pl + free_ion
      if(debug_output) then
        write(0,'(a,/,(1p2d25.15))')
     &    'free_eos_detailed: fl, tl, pnorad, (scaled) free = ',
     &    fl, tl, pnorad, free/(full_sum0*cr*t)
      endif
      free_rad = free_rad/free
      free_e = free_e/free
      free_coulomb = free_coulomb/free
      free_ex = free_ex/free
      free_pi = free_pi/free
      free_excited = free_excited/free
      free_pl = free_pl/free
      free_ion = free_ion/free
!      write(0,*)
!     &  'fractions of rad, e, coulomb, ex, pi, excited, '//
!     &  'pl, and ion components of free = '
!      write(0,*) free_rad, free_e, free_coulomb, free_ex, free_pi,
!     &  free_excited, free_pl, free_ion
!      write(0,*) 'total rad, e, coulomb, ex, pi, star, '//
!     &  'pl, and ion components of free = '
!      write(0,'(1p4d25.15)')
!     &  free*free_rad, free*free_e, free*free_coulomb,
!     &  free*free_ex, free*free_pi, free*free_excited,
!     &  free*free_pl, free*free_ion
C      adiabiatic gradient d log T/d log P
      grada = -sf/(cp*pf)
C      do in this way to avoid overflow problems
      vz = (rt/cp)*(sf/pf)
      vx = (rf/cp)*(st/pf)
      gamma1 = 1.d0/(vx-vz)
      if(kif.eq.0) then
        energy(1) = u
C        these derivatives from thermodynamic consistency arguments.
        energy(2) = pf*(rpt - rtp)*p/rho
C        it is tedious but straightforward to show this expression has
C        large significance loss in the radiation dominated case so
C        don't depend on it in that case.
        energy(3) = (pt*(rpt - rtp)*p/rho - rtp*p/rho) + cp*t
C        special tests
C        energy(1) = h2
C        energy(2) = h2f
C        energy(3) = h2t
C        energy(1) = dv(1) + dvzero(1)
C        energy(2) = dvf(1)
C        energy(3) = dvt(1)
C        index_aux = 14
C        iaux = partial_aux(index_aux)
C        energy(1) = aux(iaux)
C        energy(2) = auxf(iaux)
C        energy(3) = auxt(iaux)
      elseif(kif.eq.1) then
        energy(1) = u
C        these derivatives from thermodynamic consistency arguments.
        energy(2) = (rpt - rtp)*p/rho
        energy(3) = cp*t - rtp*p/rho
      elseif(kif.eq.2) then
        energy(1) = u
C        these derivatives from thermodynamic consistency arguments.
        energy(2) = (rpt - rtp)*p/rho/rpt
        energy(3) = cp*t - rtp*rtp*p/rpt/rho
      endif
      if(kif.eq.0) then
        enthalpy(1) = u + p/rho
C        these derivatives from thermodynamic consistency arguments.
        enthalpy(2) = pf*(1.d0-rtp)*p/rho
C        it is tedious but straightforward to show this expression has
C        large significance loss in the radiation dominated case so
C        don't depend on it in that case.
        enthalpy(3) = pt*(1.d0-rtp)*p/rho + cp*t
C        special tests
C        enthalpy(1) = dve + eta
C        enthalpy(2) = dvef + wf
C        enthalpy(3) = dvet
C        enthalpy(1) = dv(2) + dvzero(2)
C        enthalpy(2) = dvf(2)
C        enthalpy(3) = dvt(2)
C        index_aux = 15
C        iaux = partial_aux(index_aux)
C        enthalpy(1) = aux(iaux)
C        enthalpy(2) = auxf(iaux)
C        enthalpy(3) = auxt(iaux)
C        enthalpy(1) = h2plus
C        enthalpy(2) = h2plusf
C        enthalpy(3) = h2plust
C        enthalpy(1) = dve_pi
C        enthalpy(2) = dve_pif
C        enthalpy(3) = dve_pit
      elseif(kif.eq.1) then
        enthalpy(1) = u + p/rho
C        these derivatives from thermodynamic consistency arguments.
        enthalpy(2) = (1.d0-rtp)*p/rho
        enthalpy(3) = cp*t
!C        temporary special test of free derivatives.
!        if(ifnoaux_iteration) then
!C          in this case eos_jacobian not called so must combine results
!C          for u and s to calculate free energy.
!          enthalpy(1) = u -t*s
!        else
!!          enthalpy(1) = free_fp
!          enthalpy(1) = free
!        endif
!        enthalpy(2) = (p/rho)*density(2)
!        enthalpy(3) = -s*t + (p/rho)*density(3)
      elseif(kif.eq.2) then
        enthalpy(1) = u + p/rho
C        these derivatives from thermodynamic consistency arguments.
        enthalpy(2) = (1.d0-rtp)*p/rho/rpt
        enthalpy(3) = cp*t - (rtp-1.d0)*rtp*p/rpt/rho
!C        temporary special test of free derivatives.
!        if(ifnoaux_iteration) then
!C          in this case eos_jacobian not called so must combine results
!C          for u and s to calculate free energy.
!          enthalpy(1) = u -t*s
!        else
!!          enthalpy(1) = free_fp
!          enthalpy(1) = free
!        endif
!        enthalpy(2) = p/(rho)
!        enthalpy(3) = -s*t
      endif
      if(cp.lt.0.d0.and.ifprintcp.eq.1) then
        ifprintcp = 0
C        write(0,*)
C     &    'free_eos_detailed: negative C_p '//
C     &    'occurred at least once.'
C        write(0,*) 'watch out for convection implications.'
      endif
      if(grada.lt.0.d0.and.ifprintgrada.eq.1) then
        ifprintgrada = 0
C        write(0,*)
C     &    'free_eos_detailed: negative grada '//
C     &    'occurred at least once.'
C        write(0,*) 'watch out for convection implications.'
      endif
      if(rtp.lt.0.d0) then
        if(ifprintrtp.eq.1) then
          ifprintrtp = 0
C          write(0,*)
C     &      'free_eos_detailed: negative rtp  '//
C     &      'occurred at least once.'
C          write(0,*) 'watch out for convection implications.'
        endif
        cf = 0.d0
      else
        cf = cp*sqrt(rtp)
      endif
      gamma = gamma1*rf/pf
      gamma2 = 1.d0/(1.d0-grada)
      gamma3 = 1.d0+gamma1*grada
      xmu1 = ni+ne
      xmu3 = ne
      i = 2
      if(eps(1).gt.0.d0) i = 1
      vx = dne(i)*ni/ne-dni(i)
      qf = cr*((se+eta+sr/ne)*dne(i)+2.5d0*dni(i)+vx)
      qp = qf+vx*sf*p0/(p*pf)
      if(eps(1).gt.0.d0) then
        fh2 = h/eps(1)      !fraction of H+
        h2rat = 2.d0*h2/eps(1)    !fraction of H2
        h2plusrat = 2.d0*h2plus/eps(1)  !fraction of H2+
      else
        fh2 = 0.d0
        h2rat = 0.d0
        h2plusrat = 0.d0
      endif
      if(eps(2).gt.0.d0) then
        fhe2 = hd/eps(2)    !fraction of He+
        fhe3 = he/eps(2)    !fraction of He++
      else
        fhe2 = 0.d0
        fhe3 = 0.d0
      endif
      end
