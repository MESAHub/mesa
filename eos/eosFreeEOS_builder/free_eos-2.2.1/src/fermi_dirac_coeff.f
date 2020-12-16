C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fermi_dirac_coeff.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine fermi_dirac_coeff(ccoeff,morder)
C      subroutine to calculate coefficients for thermodynamically
C      consistent set of fermi-dirac integrals as in Eggleton et al.
C      1973, A&A 23,325.
C      have choice of 3rd, 5th, or 8th order thermodynamically
C      consistent approximations based on 2nd, 4th, and 7th order
C      fits to Pstar(f,g).

C      The 2nd and 4th order PMN coefficients have been completely redone,
C      and the 7th order PMN coefficients are entirely new.  All sets of
C      coefficients have been determined using eff_fit in my utilities
C      programmes with results stored in
C      utils/eff_results/eff_fit_rounded.out.
C      The 7-figure rounding should introduce relative errors of 1
C      part in 10^7 at most which is far below the fitting errors.
C      The P00 coefficient is frozen at the correct limiting value for
C      low degeneracy and degree of relativity, but most importantly the
C      fitting ranges are larger than in the EFF paper and
C      sampling is 4 times as dense both in f and g which
C      makes the fitting solution more reliable at high order.

      implicit none
      integer morder
C       n.b. ccoeff in order of rho, p, q.
      double precision ccoeff(0:morder,0:morder,3),
     &  p2coeff(0:2,0:2),p4coeff(0:4,0:4), p7coeff(0:7,0:7)

C     Coefficients for mf, mg, m =     3    3    9
C     with maximum relative residual =       0.0026100
      data p2coeff/
     &        2.3152021D0,        4.4330377D0,        2.1328286D0,
     &        5.5248695D0,        9.1657045D0,        3.3446714D0,
     &        3.6937637D0,        5.1661867D0,        1.3337625D0/
C     Coefficients for mf, mg, m =     5    5   25
C     with maximum relative residual =       0.0002298
      data p4coeff/
     &        2.3152021D0,        9.0753953D0,       13.3831949D0,
     &        8.6413847D0,        2.1332672D0,       10.1304904D0,
     &       38.2546252D0,       53.9451641D0,       33.2470532D0,
     &        7.6185569D0,       17.0830091D0,       61.9683343D0,
     &       83.0286801D0,       48.4339385D0,       10.1606269D0,
     &       12.9308973D0,       45.1991895D0,       57.6815749D0,
     &       31.6865562D0,        6.0005977D0,        3.6944723D0,
     &       12.5145460D0,       15.2592897D0,        7.9029114D0,
     &        1.3333618D0/
C     Coefficients for mf, mg, m =     8    8   64
C     with maximum relative residual =       0.0000043
      data p7coeff/
     &        2.3152021D0,       16.0295559D0,       47.5274166D0,
     &       78.3124529D0,       77.1994138D0,       45.8264429D0,
     &       15.0221403D0,        2.1333340D0,       17.0746141D0,
     &      116.7580333D0,      341.6114955D0,      554.9601789D0,
     &      538.7109181D0,      314.3226338D0,      101.1345314D0,
     &       14.0190658D0,       54.4247243D0,      367.3316300D0,
     &     1059.9499091D0,     1695.4593318D0,     1619.2183133D0,
     &      926.9175258D0,      292.1527640D0,       39.4155629D0,
     &       96.8843740D0,      645.6945133D0,     1836.9357406D0,
     &     2894.2091419D0,     2715.9112713D0,     1524.0543152D0,
     &      469.6015071D0,       61.4708699D0,      103.8623404D0,
     &      683.8011759D0,     1919.6727945D0,     2976.9551039D0,
     &     2746.6774581D0,     1508.6850884D0,      453.7692063D0,
     &       57.4337361D0,       66.9610005D0,      435.9460086D0,
     &     1208.1387561D0,     1846.0983187D0,     1673.9427291D0,
     &      899.8225245D0,      263.7719776D0,       32.1655595D0,
     &       24.0145174D0,      154.7428212D0,      423.8039766D0,
     &      638.2887657D0,      569.3358498D0,      299.3380353D0,
     &       85.4389821D0,       10.0000291D0,        3.6945284D0,
     &       23.5876179D0,       63.9072286D0,       94.9578338D0,
     &       83.3911378D0,       42.8892993D0,       11.9134229D0,
     &        1.3333330D0/
      save
      if(morder.eq.3) then
        call fermi_dirac_recursion(p2coeff, ccoeff, morder)
      elseif(morder.eq.5) then
        call fermi_dirac_recursion(p4coeff, ccoeff, morder)
      elseif(morder.eq.8) then
        call fermi_dirac_recursion(p7coeff, ccoeff, morder)
      else
        stop 'invalid morder argument to fermi_dirac_coeff'
      endif
      end
