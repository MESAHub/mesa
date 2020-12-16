C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: pl_prepare.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine pl_prepare(partial_ions, n_partial_ions, max_index,
     &  nions, nion,
     &  ifdv, bi, tc2, plop, plopt, plopt2, dv_pl, dv_plt)
C       calculate planck-larkin occupation probability and consequent
C       partial_ions(n_partial_ions+2) is an integer index corresponding to when
C         ifdv(ion) is 1.
C       changes in equilibrium constants as a function of tc2.
C       nions is the number ions (which is also the number of neutrals +
C         the number of non-bare ions).  nions+1 refers to H2, 
C         nions+2 refers to H2+.
C       nion(nions+2) is the charge of the ion associated with each species
C         neutral = 1, once ionized = 2, etc.
C       ifdv(nions+2) is in species order (H, He, He+, C, C+,...,C+++++, etc.)
C         with entries for each neutral and non-bare ion.  
C         The calling programmes
C         responsibility is to set all ifdv associated with an element
C         as 1 or zero.  If 1, then planck-larkin occupation probabilities
C         and associated changes to the equilibrium constant are
C         calculated for the neutral and non-bare ions of the element.
C         if zero, then all species of the element are ignored.
C       bi(nions+2) is the ionization potential in cm^-1 of each species
C         in species order.
C       tc2 = hc/kt
C       plop(nions+2) = ln planck-larkin occupation probabilities in species
C         order.
C       plopt(nions+2) = partial plop wrt ln t.
C       plopt2(nions+2) = partial plopt wrt ln t.
C       dv_pl(nions+2) = change in equilibrium constant in species order.
C       dv_plt(nions+2) = partial dv_pl wrt ln t
      implicit none
      integer n_partial_ions, partial_ions(n_partial_ions+2), 
     &  nions, nion(nions+2), ifdv(nions+2), ion, ion_ref, index,
     &  max_index
      double precision bi(nions+2), tc2,
     &  plop(nions+2), plopt(nions+2), plopt2(nions+2),
     &  dv_pl(nions+2), dv_plt(nions+2), arg_pl, expmarg
C       calculate log of planck-larkin occupation probabilities to be used
C       to calculate change in equilibrium constants and also later to
C       calculate planck-larkin sums in ionize.
C       n.b. ion index here refers to neutral + all non-bare ions including
C       h molecules
      do index = 1, max_index
        ion = partial_ions(index)
        arg_pl = tc2*bi(ion)
        if(arg_pl.gt.50.d0) then
          plop(ion) = 0.d0
          plopt(ion) = 0.d0
          plopt2(ion) = 0.d0
        else
          expmarg = exp(-arg_pl)
          plop(ion) = 1.d0 - expmarg*(1.d0 + arg_pl)
          plopt(ion) = -expmarg*arg_pl*arg_pl
          plopt2(ion) = expmarg*arg_pl*arg_pl*(2.d0-arg_pl)
          plopt2(ion) = (plopt2(ion) - plopt(ion)*
     &      plopt(ion)/plop(ion))/plop(ion)
          plopt(ion) = plopt(ion)/plop(ion)
          plop(ion) = log(plop(ion))
        endif
      enddo
C       calculate change in equilibrium constant = mu/kt(ref) - mu/kt(ion)
C       where (n.b.) ion now refers to first ion, second ion,.... bare nucleus
C       of each species.
      do index = 1, n_partial_ions
        ion = partial_ions(index)
C         n.b. chemical potential/kt = *negative* plop
C        This logic depends on first nion encountered for each element
C        referring to neutral monatomic of that species.  This
C        condition is fulfilled (see ionization.h).
        if(nion(ion).eq.1) ion_ref = ion
C         add effect of chemical potential of reference species
        dv_pl(ion) = -plop(ion_ref)
        dv_plt(ion) = -plopt(ion_ref)
        if(ion.ne.nions.and.nion(min(nions,ion+1)).ne.1) then
C           if not bare nucleus subtract effect of ion
C           chemical potential.  ion now has different interpretation
C           hence must add 1 to index on RHS.
          dv_pl(ion) = dv_pl(ion) + plop(ion+1)
          dv_plt(ion) = dv_plt(ion) + plopt(ion+1)
        endif
      enddo
      if(ifdv(nions+1).eq.1) then
C         H2 relative to atomic H
C         n.b. chemical potential/kt = *negative* plop
        dv_pl(nions+1) = -2.d0*plop(1) + plop(nions+1)
        dv_plt(nions+1) = -2.d0*plopt(1) + plopt(nions+1)
      else
C         plop(nions+1) used in awieos_detailed without ifdv check 
C         so should define even though multiplied by zero when ifdv
C         is zero.
        plop(nions+1) = 0.d0
        plopt(nions+1) = 0.d0
        plopt2(nions+1) = 0.d0
      endif
      if(ifdv(nions+2).eq.1) then
C         H2+ relative to *H2*
C         n.b. chemical potential/kt = *negative* plop
        dv_pl(nions+2) = -plop(nions+1) + plop(nions+2)
        dv_plt(nions+2) = -plopt(nions+1) + plopt(nions+2)
      else
C         plop(nions+2) used in awieos_detailed without ifdv check 
C         so should define even though multiplied by zero when ifdv
C         is zero.
        plop(nions+2) = 0.d0
        plopt(nions+2) = 0.d0
        plopt2(nions+2) = 0.d0
      endif
      end
