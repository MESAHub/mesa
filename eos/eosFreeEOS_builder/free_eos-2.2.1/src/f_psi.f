C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: f_psi.f 816 2008-06-24 18:43:01Z airwin $
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
C      Note:
C      The original version of this routine is in the public domain.  It was
C      written in the late 60's or early 70's by Darwin Harwood for the US
C      government, and it was brought to the attention of AWI by Forrest J.
C      Rogers and Fritz J. Swenson. AWI relicensed the routine under the GPL
C      and put in many changes to convert to structured programming and
C      double precision, evaluate partial derivatives, and insert commentary.
C
      subroutine f_psi(psi, fpsi, dfpsi, d2fpsi, d3fpsi)
C       forrest's electron exchange integral subroutine.
C       fpsi and its first through third derivatives wrt psi, where
C       psi is the Cox and Giuli degeneracy parameter eta.
C
C       definition is derived from definitions of fp12 and exct (see their
C       commentary):
C       fpsi(psi) = 
C       integral from -inf to psi of square of derivative of fp12(psi') d psi'
C       fp12(psi) is C&G F(1/2, psi)/Gamma(3/2).
C
      implicit none
      double precision psi, fpsi, dfpsi, d2fpsi, d3fpsi,
     &  f1, df1, d2f1, d3f1,
     &  f2, df2, d2f2, d3f2
      call fp12_calc(psi, f1, df1, d2f1, d3f1)
C       df2, etc., are derivatives of f2 wrt f1
      call exct_calc(f1, f2, df2, d2f2, d3f2)
      fpsi = f1*f1*f2
      dfpsi = 2.d0*f1*df1*f2 + f1*f1*df2*df1
      d2fpsi = 2.d0*(df1*df1 + f1*d2f1)*f2 +
     &  4.d0*f1*df1*df2*df1 +
     &  f1*f1*(d2f2*df1*df1 + df2*d2f1)
      d3fpsi = 2.d0*(3.d0*df1*d2f1 + f1*d3f1)*f2 +
     &  6.d0*(df1*df1 + f1*d2f1)*df2*df1 +
     &  6.d0*f1*df1*(d2f2*df1*df1 + df2*d2f1) +
     &  f1*f1*(d3f2*df1*df1*df1 + 3.d0*d2f2*df1*d2f1 + df2*d3f1)
      end
