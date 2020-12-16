C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: plsum_approx.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine plsum_approx(nlow, a, sum, dsum, dsum2)
C       calculate approximation to planck-larkin sum where
C       sum = approximation to sum from nlow to infinity of
C       2 n^2 (exp(a/n^2) -(1+a/n^2))
C       dsum is the derivative of sum wrt to a
C       dsum2 is the second derivative of sum wrt to a
      implicit none
      integer nlow, i
      double precision a, sum, dsum, dsum2
      double precision exp_max
C       this is a standard value for most of the approximations
C       which is used to limit their range of applicability
      parameter(exp_max = 1.d0)
      integer nblock, nc
      parameter(nblock = 29)
      parameter(nc = 174)
      integer istart(nblock), istop(nblock)
      double precision coeff(nc)
      data istart/
     &  1,
     &  7,
     &  13,
     &  19,
     &  25,
     &  31,
     &  37,
     &  43,
     &  49,
     &  55,
     &  61,
     &  67,
     &  73,
     &  79,
     &  85,
     &  91,
     &  97,
     &  103,
     &  109,
     &  115,
     &  121,
     &  127,
     &  133,
     &  139,
     &  145,
     &  151,
     &  157,
     &  163,
     &  169/
      data istop/
     &  6,
     &  12,
     &  18,
     &  24,
     &  30,
     &  36,
     &  42,
     &  48,
     &  54,
     &  60,
     &  66,
     &  72,
     &  78,
     &  84,
     &  90,
     &  96,
     &  102,
     &  108,
     &  114,
     &  120,
     &  126,
     &  132,
     &  138,
     &  144,
     &  150,
     &  156,
     &  162,
     &  168,
     &  174/
      data coeff/
     &  0.64493379440D0,
     &  0.027442692300D0,
     &  0.0014429412090D0,
     &  0.000069399073200D0,
     &  0.0000023183767810D0,
     &  0.00000016231026170D0,
     &  0.39493394420D0,
     &  0.0066080674370D0,
     &  0.00014296613790D0,
     &  0.0000029087927880D0,
     &  0.000000042230264560D0,
     &  0.0000000012771898410D0,
     &  0.28382288500D0,
     &  0.0024926229090D0,
     &  0.000028822472370D0,
     &  0.00000031735823540D0,
     &  0.0000000025425126030D0,
     &  4.1909882160D-11,
     &  0.22132290890D0,
     &  0.0011904792530D0,
     &  0.0000085047618610D0,
     &  0.000000058209565780D0,
     &  0.00000000029397810190D0,
     &  3.0126577460D-12,
     &  0.18132292190D0,
     &  0.00065712383250D0,
     &  0.0000031780560560D0,
     &  0.000000014773016170D0,
     &  5.1199129740D-11,
     &  3.5536844680D-13,
     &  0.15354515200D0,
     &  0.00039991247700D0,
     &  0.0000013940112180D0,
     &  0.0000000046789925720D0,
     &  1.1801540130D-11,
     &  5.8924487680D-14,
     &  0.13313699390D0,
     &  0.00026107658150D0,
     &  0.00000068646906340D0,
     &  0.0000000017399165930D0,
     &  3.3342434700D-12,
     &  1.2521406030D-14,
     &  0.11751199740D0,
     &  0.00017969370620D0,
     &  0.00000036891145780D0,
     &  0.00000000073055796580D0,
     &  1.0991856350D-12,
     &  3.2128765430D-15,
     &  0.10516632100D0,
     &  0.00012888686620D0,
     &  0.00000021226333330D0,
     &  0.00000000033733901520D0,
     &  4.0894862060D-13,
     &  9.5590406600D-16,
     &  0.095166323060D0,
     &  0.000095552545750D0,
     &  0.00000012901127490D0,
     &  1.6813362860D-10,
     &  1.6769544050D-13,
     &  3.2040707030D-16,
     &  0.086901861790D0,
     &  0.000072784782960D0,
     &  0.000000082016323360D0,
     &  8.9222806310D-11,
     &  7.4488930220D-14,
     &  1.1845535660D-16,
     &  0.079957418560D0,
     &  0.000056709238600D0,
     &  0.000000054133913460D0,
     &  4.9893950530D-11,
     &  3.5374324490D-14,
     &  4.7534713320D-17,
     &  0.074040259790D0,
     &  0.000045038001250D0,
     &  0.000000036884805220D0,
     &  2.9168052560D-11,
     &  1.7779035170D-14,
     &  2.0449386150D-17,
     &  0.068938219800D0,
     &  0.000036360836240D0,
     &  0.000000025827066500D0,
     &  1.7714305280D-11,
     &  9.3816209760D-15,
     &  9.3391113500D-18,
     &  0.064493776040D0,
     &  0.000029776307570D0,
     &  0.000000018517463670D0,
     &  1.1120075040D-11,
     &  5.1642625690D-15,
     &  4.4922784960D-18,
     &  0.060587526620D0,
     &  0.000024689918680D0,
     &  0.000000013554656220D0,
     &  7.1857933550D-12,
     &  2.9500334770D-15,
     &  2.2613901900D-18,
     &  0.057127319510D0,
     &  0.000020698809590D0,
     &  0.000000010105132760D0,
     &  4.7638023940D-12,
     &  1.7412365120D-15,
     &  1.1850307570D-18,
     &  0.054040900200D0,
     &  0.000017523404920D0,
     &  0.0000000076570711610D0,
     &  3.2308497310D-12,
     &  1.0581184140D-15,
     &  6.4358508690D-19,
     &  0.051270817480D0,
     &  0.000014965558120D0,
     &  0.0000000058872066350D0,
     &  2.2362913510D-12,
     &  6.5998646280D-16,
     &  3.6089599170D-19,
     &  0.048770817810D0,
     &  0.000012882176620D0,
     &  0.0000000045861820630D0,
     &  1.5765494200D-12,
     &  4.2143940870D-16,
     &  2.0829390140D-19,
     &  0.046503244410D0,
     &  0.000011168174040D0,
     &  0.0000000036153297030D0,
     &  1.1300522580D-12,
     &  2.7489649090D-16,
     &  1.2339678830D-19,
     &  0.044437128980D0,
     &  0.0000097451972920D0,
     &  0.0000000028809239860D0,
     &  8.2233337340D-13,
     &  1.8281106310D-16,
     &  7.4858107840D-20,
     &  0.042546770050D0,
     &  0.0000085540181700D0,
     &  0.0000000023184427720D0,
     &  6.0671858680D-13,
     &  1.2373952120D-16,
     &  4.6407612160D-20,
     &  0.040810659150D0,
     &  0.0000075493021010D0,
     &  0.0000000018827200040D0,
     &  4.5333472780D-13,
     &  8.5123757470D-17,
     &  2.9347665050D-20,
     &  0.039210659350D0,
     &  0.0000066959501480D0,
     &  0.0000000015416526840D0,
     &  3.4269265320D-13,
     &  5.9438683070D-17,
     &  1.8901821500D-20,
     &  0.037731369590D0,
     &  0.0000059665014670D0,
     &  0.0000000012721008190D0,
     &  2.6185260980D-13,
     &  4.2079208910D-17,
     &  1.2381417260D-20,
     &  0.036359627640D0,
     &  0.0000053392625460D0,
     &  0.0000000010571683890D0,
     &  2.0208332240D-13,
     &  3.0172045180D-17,
     &  8.2381910580D-21,
     &  0.035084117580D0,
     &  0.0000047969422520D0,
     &  0.00000000088437055360D0,
     &  1.5740462100D-13,
     &  2.1892057470D-17,
     &  5.5616564710D-21,
     &  0.033895057080D0,
     &  0.0000043256439000D0,
     &  0.00000000074437948750D0,
     &  1.2366354190D-13,
     &  1.6060516410D-17,
     &  3.8058628800D-21/
      if(nlow.lt.2.or.nlow.gt.nblock+1)
     &  stop 'plsum_approx: bad nlow value'
C       limit of approximation
      if(a.gt.exp_max*dble(nlow)*dble(nlow))
     &  stop 'plsum_approx: invalid a value'
C       derivatives evaluated using NR p. 137 synthetic division trick.
      sum = coeff(istop(nlow-1))
      dsum = 0.d0
      dsum2 = 0.d0
      do i = istop(nlow-1)-1, istart(nlow-1), -1
        dsum2 = dsum2*a + dsum
        dsum = dsum*a + sum
        sum = sum*a + coeff(i)
      enddo
C       multiply power series by a^2
      dsum2 = dsum2*a + dsum
      dsum = dsum*a + sum
      sum = sum*a
      dsum2 = dsum2*a + dsum
      dsum = dsum*a + sum
      sum = sum*a
C       multiply nth derivative by n factorial
      dsum2 = 2.d0*dsum2
      end
