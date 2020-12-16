C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: effective_radius.f 627 2007-07-19 01:25:19Z airwin $
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
      subroutine effective_radius(ifpi_fit,
     &  bi, nion, nions, r_ion3, r_neutral, nelements)
C       calculate effective radii
C         to be used in an MDH-style pressure ionization
C       ifpi_fit = 2,_use_best fit to Saumon table
C       ifpi_fit = 1,_use_best fit to opal table + extensions
C       ifpi_fit = 0,_use_best fit to original MDH table.
C       bi(nions+2) are the ionization potentials in cm^-1
C         (last two H2 and H2+).
C       nion(nions+2) are the ionic charge of the parent ion
C         (last two H2 and H2+).
C       r_ion3(nions+2) is the *cube of the*
C         effective radius (cm) of the ion-ion products.  aside from
C         a proportionality fitting factor this is the MDH value for the
C         ground state assuming K_n = unity.  The values are ordered
C         the same as bi, but the reference is to the lower ionization
C         state, with the bare nucleii skipped, e.g, h, he, he+, c-c+++++, etc.
C         the nions+1 and nions+2 values are for H2 and H2+
C       r_neutral(nelements+2) are the radii (cm) of the ground state of
C         the "neutral" species.  (These are actually species with
C         non-zero hard-sphere radii which happen to all be neutral in the
C         MHD model except for H2+.)  Aside from a proportionality fitting
C         factor these are the MDH values, i.e., aside from some
C         adopted values for H2, H2+, H, and He, the remaining values are
C         calculated from r = pi_fit*3 n^2 a0/2 Z, (hydrogenic values), where
C         n is the principal quantum number calculated from chi = Z^2 R/n^2.
C         The values are in order of the elements encountered for
C         bi, with the H2 and H2+ ground state radii appended as the
C         nelements+1 and nelements+2 values.
      implicit none
      include 'constants.h'
      include 'pi_fit.h'
      integer ifpi_fit, nions, nelements, nion(nions+2)
      double precision bi(nions+2), r_ion3(nions+2), 
     &    r_neutral(nelements+2)
      integer nelements_local
      parameter (nelements_local = 20)
      integer* 4 ielement, ions
      double precision rconst, principal2
      save
C      sanity checks
      if(nions.ne.nions_pi_fit) stop 'effective_radius: bad input nions'
C       calculate ion perturbation radii for all monatomic species
C       (excluding bare nucleii) and for H2 and H2+.  Assume K_n = 1
      rconst = 16.d0**(1.d0/3.d0)*echarge*echarge/ergspercmm1
      if(ifpi_fit.eq.0) then
        do ions = 1, nions+2
          r_ion3(ions) = (exp(pi_fit_ion_ln_original(ions)/3.d0)*
     &      sqrt(dble(nion(ions)))*rconst/bi(ions))**3
        enddo
      elseif(ifpi_fit.eq.1) then
        do ions = 1, nions+2
          r_ion3(ions) = (exp(pi_fit_ion_ln(ions)/3.d0)*
     &      sqrt(dble(nion(ions)))*rconst/bi(ions))**3
        enddo
      elseif(ifpi_fit.eq.2) then
        do ions = 1, nions+2
          r_ion3(ions) = (exp(pi_fit_ion_ln_saumon(ions)/3.d0)*
     &      sqrt(dble(nion(ions)))*rconst/bi(ions))**3
        enddo
      else
        stop 'effective_radius: bad ifpi_fit value input'
      endif
C       calculate hydrogenic radii for all neutrals
      ielement = 0
      do ions = 1, nions
      if(nion(ions).eq.1) then
C         only if first ion
        ielement = ielement + 1
        if(ielement.gt.nelements)
     &    stop 'effective_radius: inconsistent input'
        principal2 = rydberg/bi(ions)
        r_neutral(ielement) = principal2*bohr*
     &    (1.d0 + 0.5d0/sqrt(principal2))
      endif
      enddo
C      more sanity checks
      if(ielement.ne.nelements)
     &  stop 'effective_radius: inconsistent input'
      if(nelements.ne.nelements_pi_fit)
     &  stop 'effective_radius: bad nelements value'
C       adopt MHD values (multiplied by fitting factor later) for
C       H, He, H2, and H2+ radii
C       r_neutral(1) = bohr (see Table 3 of MHD II paper for hydrogen
C      species radii and text of same paper Section b) other elements
C      for adopted helium ground-state radius).
      r_neutral(1) = 0.529d-8
      r_neutral(2) = 0.5d-8
      r_neutral(nelements+1) = 1.45d-8
      r_neutral(nelements+2) = 1.56d-8
      if(ifpi_fit.eq.0) then
        do ielement = 1,nelements+2
          r_neutral(ielement) =
     &      exp(pi_fit_neutral_ln_original(ielement)/3.d0)*
     &      r_neutral(ielement)
        enddo
      elseif(ifpi_fit.eq.1) then
        do ielement = 1,nelements+2
          r_neutral(ielement) =
     &      exp(pi_fit_neutral_ln(ielement)/3.d0)*
     &      r_neutral(ielement)
        enddo
      elseif(ifpi_fit.eq.2) then
        do ielement = 1,nelements+2
          r_neutral(ielement) =
     &      exp(pi_fit_neutral_ln_saumon(ielement)/3.d0)*
     &      r_neutral(ielement)
        enddo
      endif
      end
