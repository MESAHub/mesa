C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: excitation_pi_end.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine excitation_pi_end(t, rho, rf, rt,
     &  pexcited, pexcitedf, pexcitedt,
     &  sexcited, sexcitedf, sexcitedt, uexcited)
      implicit none
      include 'constants.h'
      include 'excitation_block.h'
      integer index, max_index
      double precision t, rho, rf, rt, r_dv(max_index),
     &  pexcited, pexcitedf, pexcitedt,
     &  pexcited_dv(max_index),
     &  sexcited, sexcitedf, sexcitedt, uexcited
      double precision free_excited, free_excitedf,
     &  free_excited_dv(max_index)
C       free energy  = -kTV sum n(species) delta ln Z
C       psum is kept in n/(rho*avogardro) form.
      pexcited = cr*t*psum*rho
      pexcitedf = cr*t*psumf*rho + pexcited*rf
      pexcitedt = cr*t*psumt*rho + pexcited*(rt+1.d0)
C      ssum is kept in n/(rho*avogardro) form above so sexcited
C      is per unit mass.
      sexcited = cr*ssum
      sexcitedf = cr*ssumf
      sexcitedt = cr*ssumt
C       E = -T^2 partial (F/T) wrt T -->
C       energy/unit mass = R T sum (n/(rho*avogadro))
C         partial delta ln Z wrt ln t
      uexcited = cr*t*usum
      return
      entry excitation_pi_pressure_free(
     &  t, rho, rf, rt, r_dv, max_index,
     &  pexcited, pexcitedf, pexcitedt, pexcited_dv,
     &  free_excited, free_excitedf, free_excited_dv)
      pexcited = cr*t*psum*rho
      pexcitedf = cr*t*psumf*rho + pexcited*rf
      pexcitedt = cr*t*psumt*rho + pexcited*(rt+1.d0)
C      free_sum kept in n/(rho*avogardro) form so free_excited
C      is per unit mass.
      free_excited = -cr*t*free_sum
      free_excitedf = -cr*t*free_sumf
      do index = 1, max_index
        pexcited_dv(index) = cr*t*rho*
     &    (psum_dv(index) + psum*r_dv(index))
        free_excited_dv(index) = -cr*t*free_sum_dv(index)
      enddo
      end
