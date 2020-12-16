C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: lnx_calc.f 370 2006-11-29 23:57:39Z airwin $
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
      subroutine lnx_calc(sum0, sum2, t, thetaxne, lnx, dlnx, d2lnx)
C       calculate lnx = ln(x) (mdh Paper II, eq.12) 
C       and its first derivatives dlnx(4) wrt the 4 input parameters,
C       and its second derivatives d2lnx(4) wrt to the input parameters
C       (note, mixed partials are zero).
C       sum0 = sum_ions nion, where nion is the number per unit volume,
C         and the sum taken over positive ions.
C       sum2 = ne thetae + sum_ions Z^2 nion,
C         ne is the electron number density, thetae is the degeneracy
C         correction (see PTEH or notes), Z is the charge on the ions.
C       t = temperature in K
C       thetaxne = (n_e k T/P_e)*n_e, where n_e is the number density of
C         free electrons and P_e is associated pressure.
      implicit none
      include 'constants.h'
      double precision sum0, sum2, t, thetaxne, xlcon, lnx, dlnx(4), d2lnx(4)
      integer ifstart
      data ifstart/1/
      save
      if(ifstart.eq.1) then
        ifstart = 0
        xlcon = log((4.d0/3.d0)*dsqrt(pi)*echarge*echarge*echarge/
     &    boltzmann**(1.5d0))
      endif
C       lnx = ln(x) from mdh eq. 12 and assuming
C         thetaxne [equiv (n_e kT/P_e)*n_e]
C           = 3/2 F_1/2(eta)/F_3/2(eta) sum n_alpha Z,
C       The equality follows in the non-relativistic limit by definition.
C       However, n_e and P_e are defined also for the partially relativistic
C       case, so this equality is a generalization (as yet unjustified) of the 
C       MDH equation for the partially relativistic case.  This shouldn't
C       matter because in practice the tau(x) term is only used when trying
C       to exactly mimic MDH, and in order to do that you must_use_the
C       non-relativistic limit in any case.
      lnx = xlcon - 1.5d0*log(t) + log(thetaxne/sum0) + 0.5d0*log(sum2)
      dlnx(1) = -1.d0/sum0
      dlnx(2) = 0.5d0/sum2
      dlnx(3) = -1.5d0/t
      dlnx(4) = 1.d0/thetaxne
      d2lnx(1) = 1.d0/(sum0*sum0)
      d2lnx(2) = -0.5d0/(sum2*sum2)
      d2lnx(3) = 1.5d0/(t*t)
      d2lnx(4) = -1.d0/(thetaxne*thetaxne)
      end
