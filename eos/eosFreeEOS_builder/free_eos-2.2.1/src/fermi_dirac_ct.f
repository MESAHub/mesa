C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fermi_dirac_ct.f 370 2006-11-29 23:57:39Z airwin $
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
      function fermi_dirac_ct(eta,nderiv,ifsimple)
C       calculate nderiv derivative of fermi-dirac integral of order 3/2.
C       derived from cody and thacher 1967 math. comp. 21, 30.
C       n.b. see errata for table iic 1967 math. comp. 21, 525!!!!!!!!!
C       according to c-t paper, worst precision of these approximations is one
C       part in 10^(8.58).  Errors in derivatives will be worse (probably
C       something like an order of magnitude per order of derivative).
C       This routine also tested against
C       eggleton, faulkner, and flannery, generalized
C       fermi-dirac integrals and also tables in back of cox and guili vol. 2.
C       if(ifsimple = 1)_use_low degeneracy approximation.
      implicit none
      include 'constants.h'
      double precision fermi_dirac_ct,eta,p,q,rootpi,expeta,y,ratio1,top,
     &  poly_sum,bot,
     &  deriv_ratio,exp4,const4,deriv_product,ratio4,ratio5,ratio3,
     &  ratio2
      integer maxderiv,nderiv,ifsimple,j,i,ifstart,itab,ideriv,ix
      parameter (maxderiv = 5)
      dimension p(0:4,0:maxderiv,3), q(0:4,0:maxderiv,3)
C       from cody and thacher tables iic, iiic, ivc.
      data ((p(i,0,j),q(i,0,j),i=0,4),j=1,3)/
Cerrata  1  -2.34996 39854 06 d-01,  1.05000 00000 00 d 00,
     &    -2.34996 39854 06 d-01,  1.00000 00000 00 d 00,
     &    -2.92737 36375 47 d-01,  1.60859 71091 46 d 00,
     &    -9.88309 75887 38 d-02,  8.27528 95308 80 d-01,
     &    -8.25138 63795 51 d-03,  1.52232 23828 50 d-01,
     &    -1.87438 41532 23 d-05,  7.69512 04750 64 d-03,
     &    1.15302 13402    d+00,  1.00000 00000    d+00,
     &    1.05915 58972    d+00,  3.73489 53841    d-02,
     &    4.68988 03095    d-01,  2.32484 58137    d-02,
     &    1.18829 08784    d-01, -1.37667 70874    d-03,
     &    1.94387 55787    d-02,  4.64663 92781    d-05,
     &    2.46740 02368 4  d 00,  1.00000 00000 0  d 00,
     &    2.19167 58236 8  d 02,  8.91125 14061 9  d 01,
     &    1.23829 37907 5  d 04,  5.04575 66966 7  d 03,
     &    2.20667 72496 8  d 05,  9.09075 94630 4  d 04,
     &    8.49442 92003 4  d 05,  3.89960 91564 1  d 05/
      data ifstart/1/
      dimension top(0:maxderiv), bot(0:maxderiv),
     &  ratio1(0:maxderiv), ratio2(0:maxderiv), ratio3(0:maxderiv),
     &  ratio4(0:maxderiv), ratio5(0:maxderiv)
C       save rootpi
      save
      if(nderiv.lt.0.or.nderiv.gt.maxderiv)
     &  stop 'illegal nderiv argument to fermi_dirac_ct'
      if(ifstart.eq.1) then
        ifstart = 0
        rootpi = dsqrt(pi)
C         calculate derivative coefficients.
        do itab = 1,3
        do ideriv = 1,maxderiv
        do ix = 0,4-ideriv
          p(ix,ideriv,itab) = dble(ix+1)*p(ix+1,ideriv-1,itab) 
          q(ix,ideriv,itab) = dble(ix+1)*q(ix+1,ideriv-1,itab) 
        enddo
        enddo
        enddo
      endif
      if(ifsimple.eq.1) then
        expeta = dexp(eta)
        fermi_dirac_ct = 0.75d0*rootpi*expeta
      elseif(eta.le.1.d0) then
C        _use_cody and thacher table iic
        itab = 1
        y = dexp(eta)
        do ideriv = 0,nderiv
          if(ideriv.eq.0) then
            ratio1(ideriv) = y
          elseif(ideriv.eq.1) then
            ratio1(ideriv) = 1.d0
          else
            ratio1(ideriv) = 0.d0
          endif
          top(ideriv) = poly_sum(y,p(0,ideriv,itab),5-ideriv)
          bot(ideriv) = poly_sum(y,q(0,ideriv,itab),5-ideriv)
          ratio2(ideriv) = deriv_ratio(top,bot,ideriv)
          ratio3(ideriv) = deriv_product(ratio1,ratio2,ideriv)
          if(ideriv.eq.0) ratio3(0) = ratio3(0) + 0.75d0*rootpi
          ratio5(ideriv) = deriv_product(ratio1,ratio3,ideriv)
        enddo
C         ratio1 is y and its x derivatives.
        ratio1(0) = y
        do ideriv = 1,nderiv
          ratio1(ideriv) = y
        enddo
C         transform from y derivatives to x derivatives.
        if(nderiv.eq.0) then
          fermi_dirac_ct = ratio5(0)
        elseif(nderiv.eq.1) then
          fermi_dirac_ct = ratio5(1)*ratio1(1)
        elseif(nderiv.eq.2) then
          fermi_dirac_ct =
     &      ratio5(2)*ratio1(1)*ratio1(1) + ratio5(1)*ratio1(2)
        elseif(nderiv.eq.3) then
          fermi_dirac_ct = ratio5(3)*ratio1(1)*ratio1(1)*ratio1(1) +
     &      3.d0*ratio5(2)*ratio1(1)*ratio1(2) + ratio5(1)*ratio1(3)
        elseif(nderiv.eq.4) then
          fermi_dirac_ct =
     &      ratio5(4)*ratio1(1)*ratio1(1)*ratio1(1)*ratio1(1) +
     &      6.d0*ratio5(3)*ratio1(1)*ratio1(1)*ratio1(2) +
     &      4.d0*ratio5(2)*ratio1(1)*ratio1(3) +
     &      3.d0*ratio5(2)*ratio1(2)*ratio1(2) + ratio5(1)*ratio1(4)
        elseif(nderiv.eq.5) then
          fermi_dirac_ct =
     &      ratio5(5)*
     &      ratio1(1)*ratio1(1)*ratio1(1)*ratio1(1)*ratio1(1) +
     &      10.d0*ratio5(4)*ratio1(1)*ratio1(1)*ratio1(1)*ratio1(2) +
     &      10.d0*ratio5(3)*ratio1(1)*ratio1(1)*ratio1(3) +
     &      15.d0*ratio5(3)*ratio1(1)*ratio1(2)*ratio1(2) +
     &      5.d0*ratio5(2)*ratio1(1)*ratio1(4) +
     &      10.d0*ratio5(2)*ratio1(2)*ratio1(3) + ratio5(1)*ratio1(5)
        endif
      elseif(eta.le.4.d0) then
C        _use_cody and thacher table iiic
        itab = 2
        y = eta
        do ideriv = 0,nderiv
          top(ideriv) = poly_sum(y,p(0,ideriv,itab),5-ideriv)
          bot(ideriv) = poly_sum(y,q(0,ideriv,itab),5-ideriv)
        enddo
        fermi_dirac_ct = deriv_ratio(top,bot,nderiv)
      else
C        _use_cody and thacher table ivc
        itab = 3
        y = 1.d0/(eta*eta)
        do ideriv = 0,nderiv
          if(ideriv.eq.0) then
            ratio1(ideriv) = y
            exp4 = -1.25d0
            const4 = 1.d0
          elseif(ideriv.eq.1) then
            ratio1(ideriv) = 1.d0
            const4 = const4*exp4
            exp4 = exp4 - 1.d0
          else
            ratio1(ideriv) = 0.d0
            const4 = const4*exp4
            exp4 = exp4 - 1.d0
          endif
          ratio4(ideriv) = const4*y**exp4
          top(ideriv) = poly_sum(y,p(0,ideriv,itab),5-ideriv)
          bot(ideriv) = poly_sum(y,q(0,ideriv,itab),5-ideriv)
          ratio2(ideriv) = deriv_ratio(top,bot,ideriv)
          ratio3(ideriv) = deriv_product(ratio1,ratio2,ideriv)
          if(ideriv.eq.0) ratio3(0) = ratio3(0) + 0.4d0
          ratio5(ideriv) = deriv_product(ratio4,ratio3,ideriv)
        enddo
C         ratio1 is y and its x derivatives.
        ratio1(0) = y
        do ideriv = 1,nderiv
          ratio1(ideriv) = -dble(ideriv+1)*ratio1(ideriv-1)/eta
        enddo
C         transform from y derivatives to x derivatives.
        if(nderiv.eq.0) then
          fermi_dirac_ct = ratio5(0)
        elseif(nderiv.eq.1) then
          fermi_dirac_ct = ratio5(1)*ratio1(1)
        elseif(nderiv.eq.2) then
          fermi_dirac_ct =
     &      ratio5(2)*ratio1(1)*ratio1(1) + ratio5(1)*ratio1(2)
        elseif(nderiv.eq.3) then
          fermi_dirac_ct = ratio5(3)*ratio1(1)*ratio1(1)*ratio1(1) +
     &      3.d0*ratio5(2)*ratio1(1)*ratio1(2) + ratio5(1)*ratio1(3)
        elseif(nderiv.eq.4) then
          fermi_dirac_ct =
     &      ratio5(4)*ratio1(1)*ratio1(1)*ratio1(1)*ratio1(1) +
     &      6.d0*ratio5(3)*ratio1(1)*ratio1(1)*ratio1(2) +
     &      4.d0*ratio5(2)*ratio1(1)*ratio1(3) +
     &      3.d0*ratio5(2)*ratio1(2)*ratio1(2) + ratio5(1)*ratio1(4)
        elseif(nderiv.eq.5) then
          fermi_dirac_ct =
     &      ratio5(5)*
     &      ratio1(1)*ratio1(1)*ratio1(1)*ratio1(1)*ratio1(1) +
     &      10.d0*ratio5(4)*ratio1(1)*ratio1(1)*ratio1(1)*ratio1(2) +
     &      10.d0*ratio5(3)*ratio1(1)*ratio1(1)*ratio1(3) +
     &      15.d0*ratio5(3)*ratio1(1)*ratio1(2)*ratio1(2) +
     &      5.d0*ratio5(2)*ratio1(1)*ratio1(4) +
     &      10.d0*ratio5(2)*ratio1(2)*ratio1(3) + ratio5(1)*ratio1(5)
        endif
      endif
      end
