C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: free_non_ideal_calc.f 821 2008-06-24 23:12:33Z airwin $
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
      subroutine free_non_ideal_calc(
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
      implicit none
      include 'constants.h'
      integer nions, nelements, nion(nions+2)
      integer ion_end(nelements+2)
      integer ifmodified, ifexcited,
     &  nmax
      double precision bi(nions+2), t, tl,
     &  eta
      integer nstar
      parameter (nstar = 9)
      double precision rhostar(nstar), pstar(nstar),
     &  ne, nef
      integer maxfjs_aux
      parameter (maxfjs_aux = 4)
      double precision dvzero(nelements), 
     &  dve0, dve2, 
     &  dve_pi, dve_pif, dve_pit, dve_pi_aux(maxfjs_aux),
     &  dve_coulomb, dve_coulombf, dve_coulombt,
     &  dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22,
     &  plop(nions+2), plopt(nions+2), plopt2(nions+2)
      integer maxcorecharge
C      maximum core charge for non-bare ion
      parameter(maxcorecharge = 28)
      integer nmin(maxcorecharge), nmin_max(maxcorecharge),
     &  izlo, izhi
      integer nmin_species(nions+2)
      integer ifnr,
     &  ifh2, ifh2plus,
     &  if_pteh, ifcoulomb_mod, if_mc, if_dc,
     &  max_index,
     &  index, maxnextrasum,
     &  nextrasum, ifpi_local, ifpl, nxextrasum,
     &  ifsame_zero_abundances
      parameter (nxextrasum = 4)
      double precision lambda, gamma_e,
     &  full_sum0, full_sum1, full_sum2,
     &  sum0, sum0ne, sum0f, sum0t,
     &  sum0_aux(maxfjs_aux+1), sum2_aux(maxfjs_aux+1),
     &  sum2, sum2ne, sum2f, sum2t,
     &  r_ion3(nions+2), r_neutral(nelements+2)
      double precision f, wf, n_e
      integer ifdv(nions+2),
     &  partial_ions(nions+2), inv_ion(nions+2),
     &  n_partial_elements, partial_elements(nelements+2)
C      iextraoff is the number of non-extrasum auxiliary variables
      integer iextraoff, naux
      parameter (iextraoff = 7)
      logical allow_log(naux), allow_log_neg(naux)
      double precision aux_old(naux)
      integer partial_aux(naux)
      integer iaux
      integer 
     &  inv_aux(naux)
      double precision aux(naux), auxf(naux), auxt(naux),
     &  daux_dv(nions+2, naux), daux_daux(naux, naux)
      double precision bmin(maxcorecharge)
      logical 
     &  ifdvzero
      integer index_aux, jindex_aux, n_partial_aux
      double precision ddv_aux(nions+2,naux), rho,
     &  dvh, dvhf, dvht, dvh_aux(maxfjs_aux),
     &  dvhe1, dvhe1f, dvhe1t, dvhe1_aux(maxfjs_aux),
     &  dvhe2, dvhe2f, dvhe2t, dvhe2_aux(maxfjs_aux),
     &  sumion0, sumion2,
     &  sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &  nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt,
     &  nux, nuy, nuz
      integer
     &  mion_end
C      variables associated with pressure and free energy calculation:
      integer nions_local, naux_local, maxnextrasum_local
      parameter(nions_local = 295)
      parameter(naux_local = 21)
C      Note, this value must be greater than or equal to 4 as well to store
C      fjs_pi_free and corresponding pressure quantities.
      parameter(maxnextrasum_local = 9)
      double precision
     &  dv(nions_local+2), dvf(nions_local+2), dvt(nions_local+2),
     &  ppi, ppif, ppit,
     &  ppi2_daux(maxnextrasum_local), 
     &  pexcited, pexcitedf, pexcitedt, pexcited_dv(nions_local+2),
!     &  free_rad, free_e, free_ef,
     &  free_pi, free_pif, free_pit,
     &  free_pir, !free_pi_dv(nions_local+2),
     &  free_pi2_daux(maxnextrasum_local),
!     &  free_pi_daux(naux_local),
     &  free_excited, free_excitedf, free_excited_dv(nions_local+2),
     &  free_coulomb, free_coulombf, free_coulomb0, free_coulomb2,
     &  free_ex, free_exf,
!     &  free_pl, free_plf, free_pl_dv(nions_local+2),
     &  sumpl0, sumpl0f, sumpl0_dv(nions_local+2)
      logical ifnr03, ifnr13
      integer ifrad
      integer ielement, index_ion, ion
      double precision fion, fionf, fion_dv(nions+2)
      double precision free, freef, free_dv(nions_local+2),
     &  free_dauxnew(naux_local), free_daux(n_partial_aux+2)
      double precision freev, freevf, freev_dv(nions_local+2),
     &  freev_dauxnew(naux_local)
      double precision ddv_aux_dummy(nions_local+2,naux)
      !temporary
      double precision free_rad, free_e, free_pl
      save
C      sanity checks.
      if(nextrasum.gt.maxnextrasum_local)
     &  stop 'free_non_ideal_calc: maxnextrasum_local too small'
      if(.not.(ifnr.eq.0.or.ifnr.eq.1.or.ifnr.eq.3))
     &  stop 'free_non_ideal_calc: ifnr must be 0, 1, or 3'
      if(ifnr.eq.0) stop 'free_non_ideal_calc: ifnr = 0 is disabled'
      if(nions.ne.nions_local)
     &  stop 'free_non_ideal_calc: nions must be equal to nions_local'
      if(naux.ne.naux_local)
     &  stop 'free_non_ideal_calc: naux must be equal to naux_local'
      if(n_partial_aux.gt.naux_local)
     &  stop 'free_non_ideal_calc: n_partial_aux too large'
      ifnr03 = ifnr.eq.0.or.ifnr.eq.3
      ifnr13 = ifnr.eq.1.or.ifnr.eq.3
      rho = exp(aux(4))
      t = exp(tl)
      n_e = c_e*rhostar(1)

      do index = 1, n_partial_elements
        ielement = partial_elements(index)
        dvzero(ielement) = 0.d0
      enddo
      do index_ion = 1,max_index
        ion = partial_ions(index_ion)
        dv(ion) = 0.d0
        if(ifnr03) then
          dvf(ion) = 0.d0
          dvt(ion) = 0.d0
        endif
      enddo
      if(ifnr13) then
        do index_aux = 1,n_partial_aux
          do index_ion = 1,max_index
            ddv_aux_dummy(index_ion, index_aux) = 0.d0
          enddo
        enddo
      endif

C      Free energy components are organized initially into free (per mass)
C      and freev (per volume) which are a function of f, dv, and auxnew.
C      Later will combine free and freev and transform dv and aux derivatives
C      to aux (old) derivatives.
      free = 0.d0
      freev = 0.d0
      freef = 0.d0
      freevf = 0.d0
      do index_ion = 1,max_index
        free_dv(index_ion) = 0.d0
        freev_dv(index_ion) = 0.d0
      enddo
      do index_aux = 1,n_partial_aux
        free_dauxnew(index_aux) = 0.d0
        freev_dauxnew(index_aux) = 0.d0
      enddo
      
C!!!!  radiation (per volume) component.
      if(ifrad.ge.1) then
        freev = freev - prad_const*t**4
        free_rad = (- prad_const*t**4)/rho !for debug output
      endif

C!!!!  ideal (per mass) component from all particles other than free electrons.
      free = free + fion
      freef = freef + fionf
      do index_ion = 1,max_index
        free_dv(index_ion) = free_dv(index_ion) + fion_dv(index_ion)
      enddo

C!!!! ideal (per volume) component due to free electons.

      freev = freev - cpe*pstar(1) + eta*n_e*boltzmann*t
      freevf = freevf - cpe*pstar(1)*pstar(2) +
     &  (wf + eta*rhostar(2))*n_e*boltzmann*t
      free_e = (- cpe*pstar(1) + eta*n_e*boltzmann*t)/rho !for debug output

C!!!!  Coulomb (per volume) component.

C      First initialize variables for master_coulomb call.

C      n.b. for master_coulomb to work properly, sum0 and sum2
C      are considered to be functions of ne, f, t.
C      for pteh sum approximation, the sums depend on ne, and
C      there is no explicit f, t dependence.
      if(if_pteh.eq.1) then
C        PTEH (full ionization) approximation to sum0, sum2
        sum0ne = full_sum0/full_sum1
        sum0 = sum0ne*n_e
        sum2ne = full_sum2/full_sum1
        sum2 = sum2ne*n_e
        sum0f = 0.d0
        sum0t = 0.d0
        sum2f = 0.d0
        sum2t = 0.d0
      else
        sum0 = aux(6)
        sum0f = auxf(6)
        sum0t = auxt(6)
        sum2 = aux(7)
        sum2f = auxf(7)
        sum2t = auxt(7)
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
        sumion0 = sum0_mc*rho*avogadro
        sumion2 = sum2_mc*rho*avogadro
        sum0 = sumion0 + aux(5) +
     &    aux(1)*(1.d0+hcon_mc) + aux(2) + aux(3)
        sum2 = sumion2 + aux(5) +
     &    aux(1)*(1.d0+hcon_mc) + aux(2) + aux(3)*(4.d0+hecon_mc)
        if(ifnr.eq.0) then
C          N.B. this code depends on auxiliary variables that are only
C          defined for ifpi_local = 2 case.  So it only works because
C          if_mc.eq.1 is correlated with that case.  Check this.
          if(ifpi_local.ne.2)
     &      stop 'free_non_ideal_calc: logical screwup wrt ifpi_local'
          sum0f = sumion0*auxf(4) + auxf(5) +
     &      auxf(1)*(1.d0+hcon_mc) + auxf(2) + auxf(3)
          sum0t = sumion0*auxt(4) + auxt(5) +
     &      auxt(1)*(1.d0+hcon_mc) + auxt(2) + auxt(3)
          sum2f = sumion2*auxf(4) + auxf(5) +
     &      auxf(1)*(1.d0+hcon_mc) +  auxf(2) + auxf(3)*(4.d0+hecon_mc)
          sum2t = sumion2*auxt(4) + auxt(5) +
     &      auxt(1)*(1.d0+hcon_mc) +  auxt(2) + auxt(3)*(4.d0+hecon_mc)
        elseif(ifnr.eq.3) then
C          For ifnr.eq.3 calculate all f and t derivatives assuming input
C          auxiliary variables are fixed.
          sum0f = 0.d0
          sum0t = 0.d0
          sum2f = 0.d0
          sum2t = 0.d0
        endif
        if(ifnr13) then
          sum0_aux(1) = (1.d0+hcon_mc)
          sum2_aux(1) = (1.d0+hcon_mc)
          sum0_aux(2) = 1.d0
          sum2_aux(2) = 1.d0
          sum0_aux(3) = 1.d0
          sum2_aux(3) = (4.d0+hecon_mc)
          sum0_aux(4) = sumion0
          sum2_aux(4) = sumion2
          sum0_aux(5) = 1.d0
          sum2_aux(5) = 1.d0
        endif
      endif
      call master_coulomb(ifnr,
     &  rhostar, f,
     &  sum0, sum0ne, sum0f, sum0t,
     &  sum2, sum2ne, sum2f, sum2t,
     &  n_e, t, cpe*pstar(1), pstar, lambda, gamma_e,
     &  ifcoulomb_mod, if_dc, if_pteh,
     &  dve_coulomb, dve_coulombf, dve_coulombt, dve0, dve2,
     &  dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22)
      call master_coulomb_free(if_pteh,
     &  free_coulomb, free_coulombf, free_coulomb0, free_coulomb2)

      freev = freev + free_coulomb
      free_coulomb = free_coulomb/rho !for debug  output
      freevf = freevf + free_coulombf
      if(if_pteh.ne.1) then
        if(if_mc.eq.1) then
C          sum0 and sum2 calculated according to following formulas:
C          sum0 = sumion0 + h2plus + h*(1.d0+hcon_mc) +
C     &      hd + he
C          sum2 = sumion2 + h2plus + h*(1.d0+hcon_mc) +
C     &      hd + he*(4.d0+hecon_mc)
C          where h, hd, he, and h2plus are the first, second, third
C          and fifth auxiliary variables,
C          hcon_mc and hecon_mc are constants, and
C          sumion0 and sumion2 are constants times exp(rl), where rl
C          is the 4th auxiliary variable.
C          partial wrt h.
          index_aux = inv_aux(1)
          if(index_aux.gt.0)
     &      freev_dauxnew(index_aux) = freev_dauxnew(index_aux) +
     &      (1.d0+hcon_mc)*(free_coulomb0+free_coulomb2)
C          partial wrt hd.
          index_aux = inv_aux(2)
          if(index_aux.gt.0)
     &      freev_dauxnew(index_aux) = freev_dauxnew(index_aux) +
     &      (free_coulomb0+free_coulomb2)
C          partial wrt he.
          index_aux = inv_aux(3)
          if(index_aux.gt.0)
     &      freev_dauxnew(index_aux) = freev_dauxnew(index_aux) +
     &      (free_coulomb0 + (4.d0+hecon_mc)*free_coulomb2)
C          partial wrt (old) rl.
          index_aux = inv_aux(4)
          if(index_aux.gt.0)
     &      freev_dauxnew(index_aux) = freev_dauxnew(index_aux) +
     &      (sumion0*free_coulomb0 + sumion2*free_coulomb2)
C          partial wrt h2plus.
          index_aux = inv_aux(5)
          if(index_aux.gt.0)
     &      freev_dauxnew(index_aux) = freev_dauxnew(index_aux) +
     &      (free_coulomb0+free_coulomb2)
        else
C          sum0
          index_aux = inv_aux(6)
          if(index_aux.gt.0)
     &      freev_dauxnew(index_aux) = freev_dauxnew(index_aux) +
     &      free_coulomb0
C          sum2
          index_aux = inv_aux(7)
          if(index_aux.gt.0)
     &      freev_dauxnew(index_aux) = freev_dauxnew(index_aux) +
     &      free_coulomb2
        endif
      endif

C!!!!  Exchange (per volume) component
      call exchange_free(rhostar, pstar, nstar, free_ex, free_exf)
      freev = freev + free_ex
      free_ex = free_ex/rho !for debug output
      freevf = freevf + free_exf

C!!!!  Pressure ionization component (per unit mass for ifpi_local.lt.2,
C      per unit volume, otherwise.)
      if(ifpi_local.eq.1) then
        call pteh_pi(ifmodified,
     &    rhostar(1), rhostar(2), rhostar(3), t,
     &    dve_pi, dve_pif, dve_pit)
C        Note.  Unlike other free_pi which are per volume, this one is
C        per mass.
        call pteh_pi_free(full_sum1, rho, auxf(4), t, ne, nef,
     &    free_pi, free_pif, free_pir)
        free = free + free_pi
        freef = freef + free_pif
C        partial wrt ln rho.
        index_aux = inv_aux(4)
        if(index_aux.gt.0)
     &    free_dauxnew(index_aux) = free_dauxnew(index_aux) +
     &    free_pir
      elseif(ifpi_local.eq.2) then
C       _use_FJS pressure ionization
        nx = nux*rho*avogadro
        ny = nuy*rho*avogadro
        nz = nuz*rho*avogadro
        if(ifnr.eq.0) then
          nxf = nx*auxf(4)
          nxt = nx*auxt(4)
          nyf = ny*auxf(4)
          nyt = ny*auxt(4)
          nzf = nz*auxf(4)
          nzt = nz*auxt(4)
        elseif(ifnr.eq.3) then
C          For ifnr.eq.3 calculate all f and t derivatives assuming input
C          auxiliary variables (e.g., rho for ifpi_local.eq.2) are fixed.
          nxf = 0.d0
          nxt = 0.d0
          nyf = 0.d0
          nyt = 0.d0
          nzf = 0.d0
          nzt = 0.d0
        endif
        call fjs_pi(ifnr, maxfjs_aux, ifmodified,
     &    t, nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt,
!     &    h, hf, ht, hd, hdf, hdt, he, hef, het,
     &    aux(1), auxf(1), auxt(1),
     &    aux(2), auxf(2), auxt(2),
     &    aux(3), auxf(3), auxt(3),
     &    n_e, n_e*rhostar(2), n_e*rhostar(3),
     &    dvh, dvhf, dvht, dvh_aux,
     &    dvhe1, dvhe1f, dvhe1t, dvhe1_aux,
     &    dvhe2, dvhe2f, dvhe2t, dvhe2_aux,
     &    dve_pi, dve_pif, dve_pit, dve_pi_aux)
        call fjs_pi_free(nx, ny, nz,! h, hd, he,
     &    aux(1), aux(2), aux(3),
     &    t, n_e, n_e*rhostar(2), n_e*rhostar(3),
     &    free_pi, free_pif, free_pit, free_pi2_daux)
        freev = freev + free_pi
        free_pi = free_pi/rho !for debug output
        freevf = freevf + free_pif
        do index = 1,4
          index_aux = inv_aux(index)
          if(index_aux.gt.0)
     &      freev_dauxnew(index_aux) = free_pi2_daux(index)
        enddo
      elseif(ifpi_local.eq.3.or.ifpi_local.eq.4) then
C       _use_mdh-like pressure ionization and dissociation
        call mdh_pi(
     &    ifdvzero, ifnr, inv_aux, iextraoff, naux, inv_ion,
     &    partial_elements, n_partial_elements,
     &    ion_end, mion_end,
     &    ifmodified,
     &    r_ion3, nion, nions, r_neutral, nelements,
     &    aux(1+iextraoff), auxf(1+iextraoff), auxt(1+iextraoff),
     &    nextrasum,
     &    ifdv, dvzero, dv, dvf, dvt, ddv_aux_dummy)
        call mdh_pi_pressure_free(t, aux(1+iextraoff), nextrasum,
     &    ppi, ppif, ppit, ppi2_daux,
     &    free_pi, free_pif, free_pi2_daux)
        freev = freev + free_pi
        free_pi = free_pi/rho !for debug output
        freevf = freevf + free_pif
        do index = 1,nextrasum
          index_aux = inv_aux(index+iextraoff)
          if(index_aux.gt.0)
     &      freev_dauxnew(index_aux) = freev_dauxnew(index_aux) +
     &      free_pi2_daux(index)
        enddo
      endif

C!!!!  excitation (per mass) component
      if(ifexcited.gt.0) then
        call excitation_pi(ifexcited,
     &    ifpl, ifpi_local, ifmodified, ifnr, ifsame_zero_abundances,
     &    inv_aux, iextraoff, naux, inv_ion, ifh2, ifh2plus,
     &    partial_elements, n_partial_elements, ion_end, 
     &    tl, izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,
     &    bi, plop, plopt, plopt2,
     &    r_ion3, nion, nions, r_neutral, nelements,
     &    aux(1+iextraoff), auxf(1+iextraoff), auxt(1+iextraoff),
     &    nextrasum,
     &    aux(1+iextraoff+maxnextrasum),
     &    auxf(1+iextraoff+maxnextrasum),
     &    auxt(1+iextraoff+maxnextrasum),
     &    dv, dvf, dvt, ddv_aux_dummy)
        call excitation_pi_pressure_free(
     &    t, rho, auxf(4), auxt(4), daux_dv(1,4), max_index,
     &    pexcited, pexcitedf, pexcitedt, pexcited_dv,
     &    free_excited, free_excitedf, free_excited_dv)
        free = free + free_excited
        freef = freef + free_excitedf
        do index_ion = 1,max_index
          free_dv(index_ion) =
     &      free_dv(index_ion) + free_excited_dv(index_ion)
        enddo
      endif

C!!!!  Planck-Larkin (per mass) component
      if(ifpl.eq.1) then
        free = free - cr*t*sumpl0
        free_pl = - cr*t*sumpl0 !for debug output
        freef = freef - cr*t*sumpl0f
        do index_ion = 1,max_index
          free_dv(index_ion) =
     &      free_dv(index_ion) - cr*t*sumpl0_dv(index_ion)
        enddo
      endif

C      combine free and freev.
      free = free + freev/rho
      freef = freef + (freevf - freev*auxf(4))/rho
      do index_ion = 1,max_index
        free_dv(index_ion) = free_dv(index_ion) +
     &    (freev_dv(index_ion) - freev*daux_dv(index_ion,4))/rho
      enddo
      do index_aux = 1,n_partial_aux
        free_dauxnew(index_aux) = free_dauxnew(index_aux) +
     &    freev_dauxnew(index_aux)/rho
      enddo

C      transform dv partials and dauxnew partials to dauxold partials 
      do index_aux = 1,n_partial_aux
        free_daux(index_aux) = 0.d0
        do index_ion = 1,max_index
          free_daux(index_aux) = free_daux(index_aux) +
     &      free_dv(index_ion)*ddv_aux(index_ion,index_aux)
        enddo
        iaux = partial_aux(index_aux)
        freef = freef + free_dauxnew(index_aux)*auxf(iaux)
        do jindex_aux = 1,n_partial_aux
          free_daux(index_aux) = free_daux(index_aux) +
     &      free_dauxnew(jindex_aux)*daux_daux(jindex_aux, index_aux)
        enddo
      enddo
      do index_aux = 1,n_partial_aux
C        take derivative with respect to either log(auxold) or log(-auxold)
        iaux = partial_aux(index_aux)
        if(allow_log(index_aux).or.
     &      allow_log_neg(index_aux)) then
          free_daux(index_aux) =
     &      free_daux(index_aux)*aux_old(iaux)
        endif
      enddo
      free_daux(n_partial_aux+1) = freef
      free_daux(n_partial_aux+2) = 0.d0
      
!      write(0,'(a,/,(1p2d25.15))')
!     &  'free_non_ideal_calc: fl, tl, pnorad, (scaled) free = ',
!     &  log(f), tl, 0.d0 , free/(full_sum0*cr*t)
!      write(0,*)
!     &  'total rad, e, coulomb, ex, pi, excited, '//
!     &  'pl, and ion components of free = '
!      write(0,'(1p4d25.15)')
!     &  free_rad, free_e, free_coulomb,
!     &  free_ex, free_pi, free_excited,
!     &  free_pl, fion
      end
