C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: mdh_pi.f 627 2007-07-19 01:25:19Z airwin $
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
      subroutine mdh_pi(ifdvzero, ifnr, inv_aux,
     &  iextraoff, naux, inv_ion,
     &  partial_elements, n_partial_elements,
     &  ion_end, mion_end,
     &  ifmodified,
     &  r_ion3, nion, nions, r_neutral, nelements,
     &  extrasum, extrasumf, extrasumt, nextrasum,
     &  ifdv, dvzero, dv, dvf, dvt, ddv_aux)
C       the purpose of this routine is to calculate the change in equilibrium
C       constant (relative to ground state monatomic and free electrons)
C       from the appropriate combinations of pressure ionization
C       free energy/kT.
C       input:
C       ifdvzero = .true. implies dvzero (the quantity added to
C         all dv in ionize) is zero.  this is the low temperature
C         option.  For high temperatures, ifdvzero is false,
C         and dvzero is a zero point shift that renders the
C         pressure ionization contribution to dv of all bare ions zero.
C         this option gives the smallest significance loss near
C         full ionization and high densities.
C       ifnr = 0, calculate fl and tl derivatives of output using chain rule
C         and auxf and auxt
C       ifnr = 1, calculate aux derivatives of output with fl, tl fixed
C       ifnr = 2, calculate fl and tl derivatives of output with aux fixed.
C       ifnr = 3 is combination of ifnr = 1 and ifnr = 2.
C       inv_aux(naux) is pointer to actual auxiliary variable index
C       iextraoff is index offset between extrasum and aux.
C       naux is number of active auxiliary variables
C       inv_ion(nions+2) maps ion index to contiguous ion index used in NR 
C         iteration dot products.
C       partial_elements(n_partial_elements+2) index of elements treated as 
C         partially ionized consistent with ifelement.
C       ion_end(mion_end) keeps track of largest ion index for each element
C         and largest ion index for H2 (ielement = nelements+1 when
C         H2 is treated in the EOS and H2+ (ielement = nelements+2 when
C         H2+ is treated in the EOS.)
C       ifmodified > 0 means quad terms modified.
C       r_ion3(nions+2), *the cube of the*
C         effective radii (in ion order but going from neutral
C         to next to bare ion for each species) of MDH interaction between
C         all but bare nucleii species and ionized species.  last 2 are
C         H2 and H2+
C       nion(nions), charge on ion in ion order (must be same order as bi)
C       e.g., for H+, He+, He++, etc.
C       r_neutral(nelements+2) effective radii for MDH neutral-neutral 
C         interactions.  Last two are H2 and H2+ (the only ionic
C         species in the MHD model with a non-zero hard-sphere radius).
C       extrasum(nextrasum = 9) weighted sums over n(i)
C         for iextrasum = 1,nextrasum-2,
C         sum is only over neutral species (and H2+) and
C         weight is r_neutral^{iextrasum-1}
C         for iextrasum = nextrasum-1, sum is over all ionized species
C         including bare nucleii, but excluding free electrons, weight is
C         Z^1.5.
C         for iextrasum = nextrasum, sum is over all species excluding
C         bare nucleii and free electrons, the weight is rion^3.
C       extrasumf(nextrasum) = partial of extrasum/partial ln f
C       extrasumt(nextrasum) = partial of extrasum/partial ln t
C       ifdv(nions+2) = 0 means it is safe to exclude this ion from the
C       calculations
C       output:
C       dvzero(nelements) is the quantity added to dv in ionize.  it
C         is a useful zero point shift to avoid significance loss.
C       dv(nions+2) = change in equilibrium constant in ion order with last
C         two H2 and H2+.
C       dvf(nions+2) and dvt(nions+2) the ln f and ln t derivatives.
      implicit none
      logical ifdvzero
      include 'constants.h'
      include 'pi_fit.h'
      integer nions, ifnr, naux, inv_aux(naux), iextraoff,
     &  inv_ion(nions+2),
     &  n_partial_elements, partial_elements(n_partial_elements+2),
     &  mion_end, ion_end(mion_end),
     &  ifmodified, nion(nions), nelements, ion,
     &  ielement, nextrasum, ifdv(nions+2),
     &  index, max_indexe, ion_start,
     &  iextrasum, maxextrasum, index_inv_ion, index_aux
      parameter(maxextrasum = 9)
      double precision r_ion3(nions+2), r_neutral(nelements+2),
     &  extrasum(nextrasum),
     &  extrasumf(nextrasum), extrasumt(nextrasum),
     &  dvzero(nelements), 
     &  dv(nions+2), dvf(nions+2),
     &  dvt(nions+2), ddv_aux(nions+2, naux),
     &  muneutral, muneutralf,
     &  muneutralt, muneutral_aux(maxextrasum-1),
     &  muneutralh, muneutralhf, muneutralht,
     &  muneutralh_aux(maxextrasum-1),
     &  muneutralh2, muneutralh2f, muneutralh2t, 
     &  muneutralh2_aux(maxextrasum-1),
     &  muneutralh2plus, muneutralh2plusf, muneutralh2plust, 
     &  muneutralh2plus_aux(maxextrasum-1),
     &  quad,
     &  t, rho, rf, rt, squad, squadf, squadt,
     &  spi, spif, spit, upi
      double precision fdt, fdt_daux(maxextrasum),
     &  fquaddt, fquaddt_daux(maxextrasum),
     &  ppi, ppif, ppit, ppi2_daux(nextrasum),
     &  free_pi, free_pif, free_pi2_daux(nextrasum)
C      must zero so that later logic involving these arrays works properly.
      data fdt_daux, fquaddt_daux/18*0.d0/
      logical ifnr13
      save
C      calculate change in equilibrium constant from 
C      MDH-like free energy/unit volume model given by:
C      f = kt*4 pi/3 *[
C        2(3 alpha1*alpha2 + alpha0*alpha3) +
C        gamma*beta] + f'
C      f' is an ad hoc extra term designed (MDH paper II, Appendix B) to
C        force pressure ionization for low temperatures.  It has the
C        right functional form for the next higher order interaction.
C        if one expands out the various powers of r in the MDH expression, then
C      f' = quad*kt*(4 pi/3)^2 *[
C        3 alpha0*alpha3*alpha3 +
C        30 alpha1*alpha2*alpha3 +
C        9 alpha2*alpha2*alpha2 +
C        6 alpha0*alpha2*alpha4 +
C        9 alpha1*alpha1*alpha4 +
C        6 alpha0*alpha1*alpha5 +
C        1 alpha0*alpha0*alpha6]
C        
C      alphak = sum(over all neutrals + H2+) n(i)* r_neutral(i)^k
C      beta = sum(over all ions including bare nucleii, but excluding ne)
C        n(i)*nion(i)^1.5.
C      gamma = sum (all neutral and ionized species except ne, and bare nucleii) of:
C        n_i * rion(i)^3.
C      n_i and rion are in *shifted* ion order, e.g., h, he, he+, ...,
C        where each element starts at neutral and ends with next to
C        bare nucleii.
C      n.b. extrasum(1-->nextrasum-2) = alpha0-->alpha{nextrasum-3}. 
C        extrasum(nextrasum-1) = beta, and
C        extrasum(nextrasum) = gamma.
C      n.b. MDH paper II introduced quad term to force pressure ionization for
C        low T, high densities.  The quad term corresponds to alpha of
C        the MDH paper which was set to 10 for their computations.
      ifnr13 = ifnr.eq.1.or.ifnr.eq.3
      if(ifmodified.le.0) then
C        original mdh value.
        quad = quad_original
      else
        quad = quad_modified
      endif
C      sanity check:
      if(nextrasum.gt.maxextrasum)
     &  stop 'mdh_pi: nextrasum too large'
      max_indexe = n_partial_elements
      if(ifdv(nions+1).eq.1) max_indexe = max_indexe + 1
      if(ifdv(nions+2).eq.1) max_indexe = max_indexe + 1
      if(partial_elements(max_indexe).ne.mion_end) then
        write(0,*) 'nelements, max_indexe, partial_elements('//
     &    'max_indexe), mion_end = '
        write(0,*) nelements, max_indexe,
     &    partial_elements(max_indexe), mion_end
        stop 'mdh_pi: above quantities inconsistent ==> logic error'
      endif
C      n.b. ifnr=2 component is zero since fixed extrasum.  this routine
C      therefore should do nothing to dvf and dvt for ifnr = 2.
      do index = 1, max_indexe
        ielement = partial_elements(index)
        if(ielement.gt.1) then
          ion_start = ion_end(ielement-1) + 1
        else
          ion_start = 1
        endif
        do ion = ion_start, ion_end(ielement)
          if(ion.gt.nions.or.nion(min(nions,ion)).eq.1) then
C            ielement points to index of neutral monatomic or nelement+1
C            for H2 or nelement+2 for H2+.
C            Calculate chemical potential of neutral /kt.
C            neutral-neutral (where "neutral" includes any species with
C            non-zero radius (all neutrals plus H2+ according to MHD model).
            muneutral = (4.d0*pi/3.d0)*2.d0*(
     &        extrasum(4) +
     &        r_neutral(ielement)*(3.d0*extrasum(3) +
     &        r_neutral(ielement)*(3.d0*extrasum(2) +
     &        r_neutral(ielement)*(extrasum(1)))))
C            neutral-ion
            muneutral = muneutral +
     &        (4.d0*pi/3.d0)*(r_ion3(ion)*extrasum(nextrasum-1))
            if(ifnr.eq.0) then
C              calculate fl, tl derivatives assuming extrasum
C              is a function of fl and tl with specified
C              derivatives.
              muneutralf = (4.d0*pi/3.d0)*2.d0*(
     &          extrasumf(4) +
     &          r_neutral(ielement)*(3.d0*extrasumf(3) +
     &          r_neutral(ielement)*(3.d0*extrasumf(2) +
     &          r_neutral(ielement)*(extrasumf(1)))))
              muneutralt = (4.d0*pi/3.d0)*2.d0*(
     &          extrasumt(4) +
     &          r_neutral(ielement)*(3.d0*extrasumt(3) +
     &          r_neutral(ielement)*(3.d0*extrasumt(2) +
     &          r_neutral(ielement)*(extrasumt(1)))))
              muneutralf = muneutralf +
     &          (4.d0*pi/3.d0)*(r_ion3(ion)*extrasumf(nextrasum-1))
              muneutralt = muneutralt +
     &          (4.d0*pi/3.d0)*(r_ion3(ion)*extrasumt(nextrasum-1))
            elseif(ifnr13) then
              muneutral_aux(7) = 0.d0
              muneutral_aux(6) = 0.d0
              muneutral_aux(5) = 0.d0
              muneutral_aux(4) = (4.d0*pi/3.d0)*2.d0
              muneutral_aux(3) =
     &          muneutral_aux(4)*3.d0*r_neutral(ielement)
              muneutral_aux(2) =
     &          muneutral_aux(3)*r_neutral(ielement)
              muneutral_aux(1) =
     &          muneutral_aux(2)*r_neutral(ielement)/3.d0
              muneutral_aux(nextrasum-1) = (4.d0*pi/3.d0)*r_ion3(ion)
            endif
C            ad hoc neutral-neutral-neutral
            if(quad.gt.0.d0) then
              muneutral = muneutral +
     &          quad*(16.d0*pi*pi/9.d0)*(
     &          3.d0*extrasum(4)*extrasum(4) +
     &          6.d0*extrasum(3)*extrasum(5) +
     &          6.d0*extrasum(2)*extrasum(6) +
     &          2.d0*extrasum(1)*extrasum(7) +
     &          r_neutral(ielement)*(
     &          30.d0*extrasum(3)*extrasum(4) +
     &          18.d0*extrasum(2)*extrasum(5) +
     &          6.d0*extrasum(1)*extrasum(6) +
     &          r_neutral(ielement)*(
     &          27.d0*extrasum(3)*extrasum(3) +
     &          30.d0*extrasum(2)*extrasum(4) +
     &          6.d0*extrasum(1)*extrasum(5) +
     &          r_neutral(ielement)*(
     &          30.d0*extrasum(2)*extrasum(3) +
     &          6.d0*extrasum(1)*extrasum(4) +
     &          r_neutral(ielement)*(
     &          9.d0*extrasum(2)*extrasum(2) +
     &          6.d0*extrasum(1)*extrasum(3) +
     &          r_neutral(ielement)*(
     &          6.d0*extrasum(1)*extrasum(2) +
     &          r_neutral(ielement)*(
     &          1.d0*extrasum(1)*extrasum(1))))))))
              if(ifnr.eq.0) then
C                calculate fl, tl derivatives assuming extrasum
C                is a function of fl and tl with specified
C                derivatives.
                muneutralf = muneutralf +
     &            quad*(16.d0*pi*pi/9.d0)*(
     &            3.d0*(extrasumf(4)*extrasum(4) +
     &            extrasum(4)*extrasumf(4)) +
     &            6.d0*(extrasumf(3)*extrasum(5) +
     &            extrasum(3)*extrasumf(5)) +
     &            6.d0*(extrasumf(2)*extrasum(6) +
     &            extrasum(2)*extrasumf(6)) +
     &            2.d0*(extrasumf(1)*extrasum(7) +
     &            extrasum(1)*extrasumf(7)) +
     &            r_neutral(ielement)*(
     &            30.d0*(extrasumf(3)*extrasum(4) +
     &            extrasum(3)*extrasumf(4)) +
     &            18.d0*(extrasumf(2)*extrasum(5) +
     &            extrasum(2)*extrasumf(5)) +
     &            6.d0*(extrasumf(1)*extrasum(6) +
     &            extrasum(1)*extrasumf(6)) +
     &            r_neutral(ielement)*(
     &            27.d0*(extrasumf(3)*extrasum(3) +
     &            extrasum(3)*extrasumf(3)) +
     &            30.d0*(extrasumf(2)*extrasum(4) +
     &            extrasum(2)*extrasumf(4)) +
     &            6.d0*(extrasumf(1)*extrasum(5) +
     &            extrasum(1)*extrasumf(5)) +
     &            r_neutral(ielement)*(
     &            30.d0*(extrasumf(2)*extrasum(3) +
     &            extrasum(2)*extrasumf(3)) +
     &            6.d0*(extrasumf(1)*extrasum(4) +
     &            extrasum(1)*extrasumf(4)) +
     &            r_neutral(ielement)*(
     &            9.d0*(extrasumf(2)*extrasum(2) +
     &            extrasum(2)*extrasumf(2)) +
     &            6.d0*(extrasumf(1)*extrasum(3) +
     &            extrasum(1)*extrasumf(3)) +
     &            r_neutral(ielement)*(
     &            6.d0*(extrasumf(1)*extrasum(2) +
     &            extrasum(1)*extrasumf(2)) +
     &            r_neutral(ielement)*(
     &            1.d0*(extrasumf(1)*extrasum(1) +
     &            extrasum(1)*extrasumf(1)))))))))
                muneutralt = muneutralt +
     &            quad*(16.d0*pi*pi/9.d0)*(
     &            3.d0*(extrasumt(4)*extrasum(4) +
     &            extrasum(4)*extrasumt(4)) +
     &            6.d0*(extrasumt(3)*extrasum(5) +
     &            extrasum(3)*extrasumt(5)) +
     &            6.d0*(extrasumt(2)*extrasum(6) +
     &            extrasum(2)*extrasumt(6)) +
     &            2.d0*(extrasumt(1)*extrasum(7) +
     &            extrasum(1)*extrasumt(7)) +
     &            r_neutral(ielement)*(
     &            30.d0*(extrasumt(3)*extrasum(4) +
     &            extrasum(3)*extrasumt(4)) +
     &            18.d0*(extrasumt(2)*extrasum(5) +
     &            extrasum(2)*extrasumt(5)) +
     &            6.d0*(extrasumt(1)*extrasum(6) +
     &            extrasum(1)*extrasumt(6)) +
     &            r_neutral(ielement)*(
     &            27.d0*(extrasumt(3)*extrasum(3) +
     &            extrasum(3)*extrasumt(3)) +
     &            30.d0*(extrasumt(2)*extrasum(4) +
     &            extrasum(2)*extrasumt(4)) +
     &            6.d0*(extrasumt(1)*extrasum(5) +
     &            extrasum(1)*extrasumt(5)) +
     &            r_neutral(ielement)*(
     &            30.d0*(extrasumt(2)*extrasum(3) +
     &            extrasum(2)*extrasumt(3)) +
     &            6.d0*(extrasumt(1)*extrasum(4) +
     &            extrasum(1)*extrasumt(4)) +
     &            r_neutral(ielement)*(
     &            9.d0*(extrasumt(2)*extrasum(2) +
     &            extrasum(2)*extrasumt(2)) +
     &            6.d0*(extrasumt(1)*extrasum(3) +
     &            extrasum(1)*extrasumt(3)) +
     &            r_neutral(ielement)*(
     &            6.d0*(extrasumt(1)*extrasum(2) +
     &            extrasum(1)*extrasumt(2)) +
     &            r_neutral(ielement)*(
     &            1.d0*(extrasumt(1)*extrasum(1) +
     &            extrasum(1)*extrasumt(1)))))))))
              elseif(ifnr13) then
                muneutral_aux(1) = muneutral_aux(1) +
     &            quad*(16.d0*pi*pi/9.d0)*(
     &            2.d0*extrasum(7) + r_neutral(ielement)*(
     &            6.d0*extrasum(6) + r_neutral(ielement)*(
     &            6.d0*extrasum(5) + r_neutral(ielement)*(
     &            6.d0*extrasum(4) + r_neutral(ielement)*(
     &            6.d0*extrasum(3) + r_neutral(ielement)*(
     &            6.d0*extrasum(2) + r_neutral(ielement)*(
     &            2.d0*extrasum(1) )))))))
                muneutral_aux(2) = muneutral_aux(2) +
     &            quad*(16.d0*pi*pi/9.d0)*(
     &            6.d0*extrasum(6) + r_neutral(ielement)*(
     &            18.d0*extrasum(5) + r_neutral(ielement)*(
     &            30.d0*extrasum(4) + r_neutral(ielement)*(
     &            30.d0*extrasum(3) + r_neutral(ielement)*(
     &            18.d0*extrasum(2) + r_neutral(ielement)*(
     &            6.d0*extrasum(1) ))))))
                muneutral_aux(3) = muneutral_aux(3) +
     &            quad*(16.d0*pi*pi/9.d0)*(
     &            6.d0*extrasum(5) + r_neutral(ielement)*(
     &            30.d0*extrasum(4) + r_neutral(ielement)*(
     &            54.d0*extrasum(3) + r_neutral(ielement)*(
     &            30.d0*extrasum(2) + r_neutral(ielement)*(
     &            6.d0*extrasum(1) )))))
                muneutral_aux(4) = muneutral_aux(4) +
     &            quad*(16.d0*pi*pi/9.d0)*(
     &            6.d0*extrasum(4) + r_neutral(ielement)*(
     &            30.d0*extrasum(3) + r_neutral(ielement)*(
     &            30.d0*extrasum(2) + r_neutral(ielement)*(
     &            6.d0*extrasum(1) ))))
                muneutral_aux(5) = muneutral_aux(5) +
     &            quad*(16.d0*pi*pi/9.d0)*(
     &            6.d0*extrasum(3) + r_neutral(ielement)*(
     &            18.d0*extrasum(2) + r_neutral(ielement)*(
     &            6.d0*extrasum(1) )))
                muneutral_aux(6) = muneutral_aux(6) +
     &            quad*(16.d0*pi*pi/9.d0)*(
     &            6.d0*extrasum(2) + r_neutral(ielement)*(
     &            6.d0*extrasum(1) ))
                muneutral_aux(7) = muneutral_aux(7) +
     &            quad*(16.d0*pi*pi/9.d0)*(
     &            2.d0*extrasum(1) )
              endif !elseif(ifnr13)....
            endif !if(quad.gt.0.d....
            if(ion.le.nions) then
              if(ifdvzero.or.ielement.eq.1) then
C                significance loss not affected for
C                hydrogen with just one neutral and
C                one ionized state so don't bother
C                with zero point shift.
              else
C                dvzero is a zero point shift that
C                is added to dv in ionize so must
C                be subtracted from dv here.  for this
C                case_use_a zero point shift that
C                renders the bare ion pressure ionization
C                dv contribution zero.
                dvzero(ielement) = dvzero(ielement) + muneutral - 
     &            (4.d0*pi/3.d0)*(dble(nion(ion_end(ielement)))
     &            **1.5d0*
     &            extrasum(nextrasum))
              endif
            endif
          endif !if(ion.gt.nions.or.nion.....
          if(ion.le.nions) then
            if(ifdvzero.or.ielement.eq.1) then
              dv(ion) = dv(ion) + muneutral
C              add portion of negative of chemical potential of ion species/kt
              dv(ion) = dv(ion) -
     &          (4.d0*pi/3.d0)*(dble(nion(ion))**1.5d0*
     &          extrasum(nextrasum))
            else
C              subtract out dvzero without causing
C              significance loss.
C              n.b. i^1.5-j^1.5 = (i^3-j^3)/(i^1.5+j^1.5)
              dv(ion) = dv(ion) - (4.d0*pi/3.d0)*extrasum(nextrasum)*
     &          dble(nion(ion)**3-nion(ion_end(ielement))**3)/
     &          (dble(nion(ion))**1.5d0+
     &          dble(nion(ion_end(ielement)))**1.5d0)
            endif
            if(ifnr.eq.0) then
              dvf(ion) = dvf(ion) + muneutralf
              dvt(ion) = dvt(ion) + muneutralt
              dvf(ion) = dvf(ion) -
     &          (4.d0*pi/3.d0)*(dble(nion(ion))**1.5d0*
     &          extrasumf(nextrasum))
              dvt(ion) = dvt(ion) -
     &          (4.d0*pi/3.d0)*(dble(nion(ion))**1.5d0*
     &          extrasumt(nextrasum))
            elseif(ifnr13) then
              index_inv_ion = inv_ion(ion)
              do iextrasum = 1, nextrasum-1
                index_aux = inv_aux(iextraoff+iextrasum)
                ddv_aux(index_inv_ion,index_aux) =
     &            ddv_aux(index_inv_ion,index_aux) +
     &            muneutral_aux(iextrasum)
              enddo
              index_aux = inv_aux(iextraoff+nextrasum)
              ddv_aux(index_inv_ion,index_aux) =
     &          ddv_aux(index_inv_ion,index_aux) -
     &          (4.d0*pi/3.d0)*dble(nion(ion))**1.5d0
            endif
            if(ion.lt.nions.and.nion(min(nions,ion+1)).gt.1) then
C              if not last ion (usually bare nucleus, but when
C              not bare nucleus this is done consistently) add
C              in relevant additional component.
C              add negative of chemical potential of species/kt
              dv(ion) = dv(ion) -
     &          (4.d0*pi/3.d0)*(r_ion3(ion+1)*
     &          extrasum(nextrasum-1))
              if(ifnr.eq.0) then
                dvf(ion) = dvf(ion) -
     &            (4.d0*pi/3.d0)*(r_ion3(ion+1)*
     &            extrasumf(nextrasum-1))
                dvt(ion) = dvt(ion) -
     &            (4.d0*pi/3.d0)*(r_ion3(ion+1)*
     &            extrasumt(nextrasum-1))
              elseif(ifnr13) then
                index_aux = inv_aux(iextraoff+nextrasum-1)
                ddv_aux(index_inv_ion,index_aux) =
     &            ddv_aux(index_inv_ion,index_aux) -
     &            (4.d0*pi/3.d0)*r_ion3(ion+1)
              endif
            endif
          endif !if(ion.le......
          if(ion.eq.1) then
C            save muneutralh for later H2, H2+ calculation.
            muneutralh = muneutral
            if(ifnr.eq.0) then
              muneutralhf = muneutralf
              muneutralht = muneutralt
            elseif(ifnr13) then
              do iextrasum = 1, nextrasum-1
                muneutralh_aux(iextrasum) = muneutral_aux(iextrasum)
              enddo
            endif
          endif
          if(ion.eq.nions+1) then
C            save muneutralh2 for later H2, H2+ calculation.
            muneutralh2 = muneutral
            if(ifnr.eq.0) then
              muneutralh2f = muneutralf
              muneutralh2t = muneutralt
            elseif(ifnr13) then
              do iextrasum = 1, nextrasum-1
                muneutralh2_aux(iextrasum) = muneutral_aux(iextrasum)
              enddo
            endif
          endif
          if(ion.eq.nions+2) then
C            save muneutralh2plus for later H2+ calculation.
            muneutralh2plus = muneutral
            if(ifnr.eq.0) then
              muneutralh2plusf = muneutralf
              muneutralh2plust = muneutralt
            elseif(ifnr13) then
              do iextrasum = 1, nextrasum-1
                muneutralh2plus_aux(iextrasum) =
     &            muneutral_aux(iextrasum)
              enddo
            endif
          endif
        enddo !do ion = ion_start,....
      enddo !do index = 1,.......
C      special handling of H2, and H2+.  At this point:
C      muneutralh holds partial F partial n(H)/kt,
C      muneutralh2 holds partial F partial n(H2)/kt
C      muneutralh2plus holds partial F partial n(H2+)/kt
      if(ifdv(nions+1).eq.1) then
C        by definition of H2 equilibrium constant relative to neutral monatomic
        dv(nions+1) = dv(nions+1) + 2.d0*muneutralh - muneutralh2
        if(ifnr.eq.0) then
          dvf(nions+1) = dvf(nions+1) + 2.d0*muneutralhf -
     &      muneutralh2f
          dvt(nions+1) = dvt(nions+1) + 2.d0*muneutralht -
     &      muneutralh2t
        elseif(ifnr13) then
          index_inv_ion = inv_ion(nions+1)
          do iextrasum = 1, nextrasum-1
            index_aux = inv_aux(iextraoff+iextrasum)
            ddv_aux(index_inv_ion,index_aux) =
     &        ddv_aux(index_inv_ion,index_aux) +
     &        2.d0*muneutralh_aux(iextrasum) -
     &        muneutralh2_aux(iextrasum)
          enddo
        endif
      endif
      if(ifdv(nions+2).eq.1) then
C        by definition of H2+ equilibrium constant relative to H2 and e-.
C        n.b. the extrasum(nextrasum-1) component is already taken care
C        of in muneutralh2plus
        dv(nions+2) = dv(nions+2) + muneutralh2 -
     &    muneutralh2plus -
     &    (4.d0*pi/3.d0)*extrasum(nextrasum)
        if(ifnr.eq.0) then
          dvf(nions+2) = dvf(nions+2) + muneutralh2f -
     &      muneutralh2plusf -
     &      (4.d0*pi/3.d0)*extrasumf(nextrasum)
          dvt(nions+2) = dvt(nions+2) + muneutralh2t -
     &      muneutralh2plust -
     &      (4.d0*pi/3.d0)*extrasumt(nextrasum)
        elseif(ifnr13) then
          index_inv_ion = inv_ion(nions+2)
          do iextrasum = 1, nextrasum-1
            index_aux = inv_aux(iextraoff+iextrasum)
            ddv_aux(index_inv_ion,index_aux) =
     &        ddv_aux(index_inv_ion,index_aux) +
     &        muneutralh2_aux(iextrasum) -
     &        muneutralh2plus_aux(iextrasum)
          enddo
          index_aux = inv_aux(iextraoff+nextrasum)
          ddv_aux(index_inv_ion,index_aux) =
     &      ddv_aux(index_inv_ion,index_aux) -
     &      (4.d0*pi/3.d0)
        endif
      endif
      return
      entry mdh_pi_pressure_free(t, extrasum, nextrasum,
     &  ppi, ppif, ppit, ppi2_daux,
     &  free_pi, free_pif, free_pi2_daux)
C      sanity check:
      if(nextrasum.gt.maxextrasum)
     &  stop 'mdh_pi_pressure_free: nextrasum too large'
C      MDH-like free energy/unit volume model given by:
C      f = kt*4 pi/3 *[
C        2(3 alpha1*alpha2 + alpha0*alpha3) +
C        gamma*beta] + f'
C      f' is an ad hoc extra term designed (MDH paper II, Appendix B) to
C        force pressure ionization for low temperatures.  It has the
C        right functional form for the next higher order interaction.
C        if one expands out the various powers of r in the MDH expression, then
C      f' = quad*kt*(4 pi/3)^2 *[
C        3 alpha0*alpha3*alpha3 +
C        30 alpha1*alpha2*alpha3 +
C        9 alpha2*alpha2*alpha2 +
C        6 alpha0*alpha2*alpha4 +
C        9 alpha1*alpha1*alpha4 +
C        6 alpha0*alpha1*alpha5 +
C        1 alpha0*alpha0*alpha6]
C      Calculate f and fquad divided by t.
C      n.b. extrasum is in number per unit volume form
      fdt = boltzmann*(4.d0*pi/3.d0)*(
     &  2.d0*(3.d0*extrasum(2)*extrasum(3) +
     &  extrasum(1)*extrasum(4)) +
     &  extrasum(nextrasum-1)*extrasum(nextrasum))
      fdt_daux(1) = boltzmann*(4.d0*pi/3.d0)*
     &  (2.d0*extrasum(4))
      fdt_daux(2) = boltzmann*(4.d0*pi/3.d0)*
     &  (6.d0*extrasum(3))
      fdt_daux(3) = boltzmann*(4.d0*pi/3.d0)*
     &  (6.d0*extrasum(2))
      fdt_daux(4) = boltzmann*(4.d0*pi/3.d0)*
     &  (2.d0*extrasum(1))
      fdt_daux(nextrasum-1) = boltzmann*(4.d0*pi/3.d0)*
     &  extrasum(nextrasum)
      fdt_daux(nextrasum) = boltzmann*(4.d0*pi/3.d0)*
     &  extrasum(nextrasum-1)
      fquaddt = quad*boltzmann*(16.d0*pi*pi/9.d0)*(
     &  3.d0*extrasum(1)*extrasum(4)*extrasum(4) +
     &  30.d0*extrasum(2)*extrasum(3)*extrasum(4) +
     &  9.d0*extrasum(3)*extrasum(3)*extrasum(3) +
     &  6.d0*extrasum(1)*extrasum(3)*extrasum(5) +
     &  9.d0*extrasum(2)*extrasum(2)*extrasum(5) +
     &  6.d0*extrasum(1)*extrasum(2)*extrasum(6) +
     &  1.d0*extrasum(1)*extrasum(1)*extrasum(7))
      fquaddt_daux(1) = boltzmann*(4.d0*pi/3.d0)*
     &  (
     &  quad*(4.d0*pi/3.d0)*(
     &  3.d0*extrasum(4)*extrasum(4) +
     &  6.d0*extrasum(3)*extrasum(5) +
     &  6.d0*extrasum(2)*extrasum(6) +
     &  2.d0*extrasum(1)*extrasum(7)
     &  ))
      fquaddt_daux(2) = boltzmann*(4.d0*pi/3.d0)*
     &  (
     &  quad*(4.d0*pi/3.d0)*(
     &  30.d0*extrasum(3)*extrasum(4) +
     &  18.d0*extrasum(2)*extrasum(5) +
     &  6.d0*extrasum(1)*extrasum(6)
     &  ))
      fquaddt_daux(3) = boltzmann*(4.d0*pi/3.d0)*
     &  (
     &  quad*(4.d0*pi/3.d0)*(
     &  30.d0*extrasum(2)*extrasum(4) +
     &  27.d0*extrasum(3)*extrasum(3) +
     &  6.d0*extrasum(1)*extrasum(5)
     &  ))
      fquaddt_daux(4) = boltzmann*(4.d0*pi/3.d0)*
     &  (
     &  quad*(4.d0*pi/3.d0)*(
     &  6.d0*extrasum(1)*extrasum(4) +
     &  30.d0*extrasum(2)*extrasum(3) +
     &  6.d0*extrasum(1)*extrasum(5)
     &  ))
      fquaddt_daux(5) = boltzmann*(4.d0*pi/3.d0)*
     &  (
     &  quad*(4.d0*pi/3.d0)*(
     &  6.d0*extrasum(1)*extrasum(3) +
     &  9.d0*extrasum(2)*extrasum(2)
     &  ))
      fquaddt_daux(6) = boltzmann*(4.d0*pi/3.d0)*
     &  (
     &  quad*(4.d0*pi/3.d0)*(
     &  6.d0*extrasum(1)*extrasum(2)
     &  ))
      fquaddt_daux(7) = boltzmann*(4.d0*pi/3.d0)*
     &  (
     &  quad*(4.d0*pi/3.d0)*(
     &  1.d0*extrasum(1)*extrasum(1)
     &  ))
C      F(T,V,N) = t*V*fdt + t*V*fquaddt.
C      additional V factors divide to convert extrasum products to N form.
C      Therefore, first term is proportional to V^{-1} and
C      second term is proportional to V^{-2}
C      P = -partial F(T,V,N)/partial V
      ppi = t*(fdt + 2.d0*fquaddt)
      ppif = 0.d0
      ppit = ppi
      free_pi = t*(fdt + fquaddt)
      free_pif = 0.d0
C      n.b. some parts of these arrays are initialized to zero by
C      data statement earlier in routine.
      do iextrasum = 1, nextrasum
        ppi2_daux(iextrasum) = t*(fdt_daux(iextrasum) +
     &    2.d0*fquaddt_daux(iextrasum))
        free_pi2_daux(iextrasum) = t*(fdt_daux(iextrasum) +
     &    fquaddt_daux(iextrasum))
      enddo
      return
      entry mdh_pi_end(t, rho, rf, rt,
     &  r_ion3, nion, nions, r_neutral, nelements,
     &  extrasum, extrasumf, extrasumt, nextrasum,
     &  ppi, ppif, ppit, spi, spif, spit, upi)
C      compute remaining quantities having found ionization balance
C      MDH-like free energy/unit volume model given by:
C      f = kt*4 pi/3 *[
C        2(3 alpha1*alpha2 + alpha0*alpha3) +
C        gamma*beta] + f'
C      f' is an ad hoc extra term designed (MDH paper II, Appendix B) to
C        force pressure ionization for low temperatures.  It has the
C        right functional form for the next higher order interaction.
C        if one expands out the various powers of r in the MDH expression, then
C      f' = quad*kt*(4 pi/3)^2 *[
C        3 alpha0*alpha3*alpha3 +
C        30 alpha1*alpha2*alpha3 +
C        9 alpha2*alpha2*alpha2 +
C        6 alpha0*alpha2*alpha4 +
C        9 alpha1*alpha1*alpha4 +
C        6 alpha0*alpha1*alpha5 +
C        1 alpha0*alpha0*alpha6]
C      spi = s per unit mass = - partial fV/partial T/m = -f/(rho*T)
C      correct for quad part later.
C      n.b. extrasum is in nu = n/(rho*avogadro) form for mdh_pi_end call
      spi = -cr*(4.d0*pi/3.d0)*(
     &  2.d0*(3.d0*extrasum(2)*extrasum(3) +
     &  extrasum(1)*extrasum(4)) +
     &  extrasum(nextrasum-1)*extrasum(nextrasum))*(rho*avogadro)
      spif = -cr*(4.d0*pi/3.d0)*(
     &  2.d0*(3.d0*(extrasumf(2)*extrasum(3) +
     &  extrasum(2)*extrasumf(3)) +
     &  extrasumf(1)*extrasum(4) + extrasum(1)*extrasumf(4)) +
     &  extrasumf(nextrasum-1)*extrasum(nextrasum) +
     &  extrasum(nextrasum-1)*extrasumf(nextrasum))*(rho*avogadro) +
     &  spi*rf
      spit = -cr*(4.d0*pi/3.d0)*(
     &  2.d0*(3.d0*(extrasumt(2)*extrasum(3) +
     &  extrasum(2)*extrasumt(3)) +
     &  extrasumt(1)*extrasum(4) + extrasum(1)*extrasumt(4)) +
     &  extrasumt(nextrasum-1)*extrasum(nextrasum) +
     &  extrasum(nextrasum-1)*extrasumt(nextrasum))*(rho*avogadro) +
     &  spi*rt
C      p = -partial (fV)/partial V
C      fV proportional to V^-1 and ppi = f = -spi*rho*t
      ppi = -spi*rho*t
      ppif = -(spif+spi*rf)*rho*t
      ppit = -(spit+spi*(1.d0+rt))*rho*t
C      upi = internal energy/unit mass = fpi/rho + t*spi
C        = -spi*t + t*spi = 0
      upi = 0.d0
      if(quad.le.0.d0) return
C      s per unit mass = - partial fquad V/partial T/m = -fquad/(rho*t)
C      n.b. extrasum is in nu = n/(rho*avogadro) form for mdh_pi_end call
      squad = -quad*cr*(16.d0*pi*pi/9.d0)*(
     &  3.d0*extrasum(1)*extrasum(4)*extrasum(4) +
     &  30.d0*extrasum(2)*extrasum(3)*extrasum(4) +
     &  9.d0*extrasum(3)*extrasum(3)*extrasum(3) +
     &  6.d0*extrasum(1)*extrasum(3)*extrasum(5) +
     &  9.d0*extrasum(2)*extrasum(2)*extrasum(5) +
     &  6.d0*extrasum(1)*extrasum(2)*extrasum(6) +
     &  1.d0*extrasum(1)*extrasum(1)*extrasum(7))*
     &  (rho*avogadro)*(rho*avogadro)
      squadf = -quad*cr*(16.d0*pi*pi/9.d0)*(
     &  3.d0*extrasumf(1)*extrasum(4)*extrasum(4) +
     &  3.d0*extrasum(1)*extrasumf(4)*extrasum(4) +
     &  3.d0*extrasum(1)*extrasum(4)*extrasumf(4) +
     &  30.d0*extrasumf(2)*extrasum(3)*extrasum(4) +
     &  30.d0*extrasum(2)*extrasumf(3)*extrasum(4) +
     &  30.d0*extrasum(2)*extrasum(3)*extrasumf(4) +
     &  9.d0*extrasumf(3)*extrasum(3)*extrasum(3) +
     &  9.d0*extrasum(3)*extrasumf(3)*extrasum(3) +
     &  9.d0*extrasum(3)*extrasum(3)*extrasumf(3) +
     &  6.d0*extrasumf(1)*extrasum(3)*extrasum(5) +
     &  6.d0*extrasum(1)*extrasumf(3)*extrasum(5) +
     &  6.d0*extrasum(1)*extrasum(3)*extrasumf(5) +
     &  9.d0*extrasumf(2)*extrasum(2)*extrasum(5) +
     &  9.d0*extrasum(2)*extrasumf(2)*extrasum(5) +
     &  9.d0*extrasum(2)*extrasum(2)*extrasumf(5) +
     &  6.d0*extrasumf(1)*extrasum(2)*extrasum(6) +
     &  6.d0*extrasum(1)*extrasumf(2)*extrasum(6) +
     &  6.d0*extrasum(1)*extrasum(2)*extrasumf(6) +
     &  1.d0*extrasumf(1)*extrasum(1)*extrasum(7) +
     &  1.d0*extrasum(1)*extrasumf(1)*extrasum(7) +
     &  1.d0*extrasum(1)*extrasum(1)*extrasumf(7))*
     &  (rho*avogadro)*(rho*avogadro) + 2.d0*squad*rf
      squadt = -quad*cr*(16.d0*pi*pi/9.d0)*(
     &  3.d0*extrasumt(1)*extrasum(4)*extrasum(4) +
     &  3.d0*extrasum(1)*extrasumt(4)*extrasum(4) +
     &  3.d0*extrasum(1)*extrasum(4)*extrasumt(4) +
     &  30.d0*extrasumt(2)*extrasum(3)*extrasum(4) +
     &  30.d0*extrasum(2)*extrasumt(3)*extrasum(4) +
     &  30.d0*extrasum(2)*extrasum(3)*extrasumt(4) +
     &  9.d0*extrasumt(3)*extrasum(3)*extrasum(3) +
     &  9.d0*extrasum(3)*extrasumt(3)*extrasum(3) +
     &  9.d0*extrasum(3)*extrasum(3)*extrasumt(3) +
     &  6.d0*extrasumt(1)*extrasum(3)*extrasum(5) +
     &  6.d0*extrasum(1)*extrasumt(3)*extrasum(5) +
     &  6.d0*extrasum(1)*extrasum(3)*extrasumt(5) +
     &  9.d0*extrasumt(2)*extrasum(2)*extrasum(5) +
     &  9.d0*extrasum(2)*extrasumt(2)*extrasum(5) +
     &  9.d0*extrasum(2)*extrasum(2)*extrasumt(5) +
     &  6.d0*extrasumt(1)*extrasum(2)*extrasum(6) +
     &  6.d0*extrasum(1)*extrasumt(2)*extrasum(6) +
     &  6.d0*extrasum(1)*extrasum(2)*extrasumt(6) +
     &  1.d0*extrasumt(1)*extrasum(1)*extrasum(7) +
     &  1.d0*extrasum(1)*extrasumt(1)*extrasum(7) +
     &  1.d0*extrasum(1)*extrasum(1)*extrasumt(7))*
     &  (rho*avogadro)*(rho*avogadro) + 2.d0*squad*rt
      spi = spi + squad
      spif = spif + squadf
      spit = spit + squadt
C      n.b. fquadV is proportional to V^-2, so 
C        pquad = 2*fquad = -2*squad*rho*t
      ppi = ppi - 2.d0*squad*rho*t
      ppif = ppif - 2.d0*(squadf+squad*rf)*rho*t
      ppit = ppit - 2.d0*(squadt+squad*(1.d0+rt))*rho*t
C      upi = internal energy/unit mass = fquad/rho + t*squad
C        = pquad/(2*rho) + t*squad = -squad*t + t*squad = 0
C      upi = upi + 0.d0
      end
