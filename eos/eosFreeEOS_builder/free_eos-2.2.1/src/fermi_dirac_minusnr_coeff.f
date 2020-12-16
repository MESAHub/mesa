C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fermi_dirac_minusnr_coeff.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine fermi_dirac_minusnr_coeff(ccoeff ,morder)
C      subroutine to calculate coefficients for thermodynamically
C      consistent set of fermi-dirac integrals as in Eggleton et al.
C      1973, A&A 23,325, but with non-relativistic limit subtracted. There
C      is choice of 3rd, 5th, or 8th order thermodynamically
C      consistent approximations based on 2nd, 4th, and 7th order
C      fits to Pstar(f,g).
C
C      All sets of
C      coefficients have been determined using eff_cdfit in my utilities
C      programmes with results stored in
C      utils/eff_results/eff_cdfit_rounded.out.
C      The 7-figure rounding should introduce relative errors of 1
C      part in 10^7 at most which is far below the fitting errors.
C      The low-order g coefficients were frozen at zero consistent
C      with subtracting off the non-relativistic result.
C
C      When these coefficients are used to calculate an f, g multinomial
C      approximation and the results combined with a precise approximation
C      of the non-relativistic limit (such as the Cody-Thacher approximation),
C      the results are better in that important limit than the pure EFF-style
C      multinomial approximation for the full Fermi-Dirac integrals.
      implicit none
      integer morder
C       n.b. ccoeff in order of rho, p, q.
      double precision ccoeff(0:morder,0:morder,3),
     &  p2coeff(0:2,0:2),p4coeff(0:4,0:4), p7coeff(0:7,0:7)

C     Coefficients for mf, mg, m =     3    3    9
C     with maximum relative residual =       0.0027119
      data p2coeff/
     &        0.0000000D0,        0.0000000D0,        0.0000000D0,
     &        5.5051486D0,        9.1814586D0,        3.3509116D0,
     &        3.6945576D0,        5.1656731D0,        1.3336029D0/
C     Coefficients for mf, mg, m =     5    5   25
C     with maximum relative residual =       0.0001772
      data p4coeff/
     &        0.0000000D0,        0.0000000D0,        0.0000000D0,
     &        0.0000000D0,        0.0000000D0,       10.1292483D0,
     &       38.2616936D0,       53.9322152D0,       33.2523591D0,
     &        7.6188636D0,       17.0862052D0,       61.9601529D0,
     &       83.0404668D0,       48.4264382D0,       10.1598500D0,
     &       12.9298674D0,       45.2012106D0,       57.6788896D0,
     &       31.6885718D0,        6.0008022D0,        3.6944790D0,
     &       12.5145345D0,       15.2593039D0,        7.9029001D0,
     &        1.3333607D0/
C     Coefficients for mf, mg, m =     8    8   64
C     with maximum relative residual =       0.0000045
      data p7coeff/
     &        0.0000000D0,        0.0000000D0,        0.0000000D0,
     &        0.0000000D0,        0.0000000D0,        0.0000000D0,
     &        0.0000000D0,        0.0000000D0,       17.0746137D0,
     &      116.7565267D0,      341.6206921D0,      554.9387132D0,
     &      538.7259249D0,      314.3185834D0,      101.1346265D0,
     &       14.0190497D0,       54.4249544D0,      367.3451784D0,
     &     1059.8746088D0,     1695.6501699D0,     1619.0830291D0,
     &      926.9575547D0,      292.1522041D0,       39.4157661D0,
     &       96.8834221D0,      645.6603616D0,     1837.1092878D0,
     &     2893.7518429D0,     2716.2321541D0,     1523.9549922D0,
     &      469.6022201D0,       61.4703123D0,      103.8633443D0,
     &      683.8297286D0,     1919.5366558D0,     2977.3218243D0,
     &     2746.4224691D0,     1508.7658062D0,      453.7688812D0,
     &       57.4341957D0,       66.9606940D0,      435.9384048D0,
     &     1208.1733455D0,     1846.0041046D0,     1674.0076566D0,
     &      899.8017708D0,      263.7720349D0,       32.1654442D0,
     &       24.0145386D0,      154.7433006D0,      423.8018745D0,
     &      638.2945138D0,      569.3319258D0,      299.3392925D0,
     &       85.4389791D0,       10.0000358D0,        3.6945283D0,
     &       23.5876169D0,       63.9072327D0,       94.9578227D0,
     &       83.3911452D0,       42.8892970D0,       11.9134229D0,
     &        1.3333330D0/
      save
      if(morder.eq.3) then
        call fermi_dirac_recursion(p2coeff, ccoeff, morder)
      elseif(morder.eq.5) then
        call fermi_dirac_recursion(p4coeff, ccoeff, morder)
      elseif(morder.eq.8) then
        call fermi_dirac_recursion(p7coeff, ccoeff, morder)
      else
        stop 'invalid morder argument to fermi_dirac_minusnr_coeff'
      endif
      end
