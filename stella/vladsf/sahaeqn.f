c------------------------------------------------------------------------
c            s u b r o u t i n e     s a h a e q n
c------------------------------------------------------------------------
c This routine computes the ionization fraction of an element by solving
c the Saha-Boltzmann equation. This routine contains data on parition
c functions, mostly in the form of ground configuration statistical weights,
c and ionization potentials.

      subroutine sahaeqn(tempflag, kurucz, press, temp, eden, z,
     .     maxion, firstion, fraction, partition)
      implicit real*8 (a-h, o-z)
      intent(in) tempflag, kurucz, press, temp, eden, z, maxion
      intent(out) firstion, fraction, partition

c Input parameters:

c (integer) TEMPFLAG determines whether, for ionization stages greater than
c  the first six, whether the partition function is taken to equal the
c  statistical weight of the ground term (a cold gas), or the statistical
c  weight of the ground configuration (a ``hot'' gas). 
c   TEMPFLAG = 0 => assume gas is cold.
c   TEMPFLAG = 1 => assume gas is hot.
c (logical) KURUCZ determines whether or not the first six ionization stage 
c  partition functions should be taken from PFSAHA.
c   KURUCZ = T => take first six ionization stage partition from PFSAHA.
c (real*8) PRESS must contain the pressure if KURUCZ is set to true.
c (real*8) TEMP is the temperature.
c (real*8) EDEN is the electron density.
c (integer) Z is the atomic number.
c (integer) MAXION is the maximum number of ionization stages to return
c  in FRACTION. For element Z, the routine will return the fraction of ions 
c  of the MAXION most abundant and successive ionization stages, starting
c  with ionization stages FISRTION.
c (integer) FIRSTION is the ionization stage of the first ion whose fraction
c  is returned in FRACTION. FIRSTION >= 1, which denotes the neutral atom.
c (real*8) PARTITION is used to return the assumed parition function
c for each of the maxion ions.

      integer tempflag, z, maxion, firstion
      logical kurucz
      real*8 press, temp, eden, fraction(0:maxion), partition(0:maxion)
      save init
      real*8 frac(0:100), partit(0:100)

      real*8 jg(465), jc(465), ionpot(465)

      real*8 jgH(1), jcH(1), ionpotH(1)
      real*8 ionpotHminus(1)
      equivalence(jg(1), jgH(1))
      equivalence(jc(1), jcH(1))
      equivalence(ionpot(1), ionpotH(1))
      real*8 jgHe(2), jcHe(2), ionpotHe(2)
      equivalence(jg(2), jgHe(1))
      equivalence(jc(2), jcHe(1))
      equivalence(ionpot(2), ionpotHe(1))
      real*8 jgLi(3), jcLi(3), ionpotLi(3)
      equivalence(jg(4), jgLi(1))
      equivalence(jc(4), jcLi(1))
      equivalence(ionpot(4), ionpotLi(1))
      real*8 jgBe(4), jcBe(4), ionpotBe(4)
      equivalence(jg(7), jgBe(1))
      equivalence(jc(7), jcBe(1))
      equivalence(ionpot(7), ionpotBe(1))
      real*8 jgB(5), jcB(5), ionpotB(5)
      equivalence(jg(11), jgB(1))
      equivalence(jc(11), jcB(1))
      equivalence(ionpot(11), ionpotB(1))
      real*8 jgC(6), jcC(6), ionpotC(6)
      equivalence(jg(16), jgC(1))
      equivalence(jc(16), jcC(1))
      equivalence(ionpot(16), ionpotC(1))
      real*8 jgN(7), jcN(7), ionpotN(7)
      equivalence(jg(22), jgN(1))
      equivalence(jc(22), jcN(1))
      equivalence(ionpot(22), ionpotN(1))
      real*8 jgO(8), jcO(8), ionpotO(8)
      equivalence(jg(29), jgO(1))
      equivalence(jc(29), jcO(1))
      equivalence(ionpot(29), ionpotO(1))
      real*8 jgF(9), jcF(9), ionpotF(9)
      equivalence(jg(37), jgF(1))
      equivalence(jc(37), jcF(1))
      equivalence(ionpot(37), ionpotF(1))
      real*8 jgNe(10), jcNe(10), ionpotNe(10)
      equivalence(jg(46), jgNe(1))
      equivalence(jc(46), jcNe(1))
      equivalence(ionpot(46), ionpotNe(1))
      real*8 jgNa(11), jcNa(11), ionpotNa(11)
      equivalence(jg(56), jgNa(1))
      equivalence(jc(56), jcNa(1))
      equivalence(ionpot(56), ionpotNa(1))
      real*8 jgMg(12), jcMg(12), ionpotMg(12)
      equivalence(jg(67), jgMg(1))
      equivalence(jc(67), jcMg(1))
      equivalence(ionpot(67), ionpotMg(1))
      real*8 jgAl(13), jcAl(13), ionpotAl(13)
      equivalence(jg(79), jgAl(1))
      equivalence(jc(79), jcAl(1))
      equivalence(ionpot(79), ionpotAl(1))
      real*8 jgSi(14), jcSi(14), ionpotSi(14)
      equivalence(jg(92), jgSi(1))
      equivalence(jc(92), jcSi(1))
      equivalence(ionpot(92), ionpotSi(1))
      real*8 jgP(15), jcP(15), ionpotP(15)
      equivalence(jg(106), jgP(1))
      equivalence(jc(106), jcP(1))
      equivalence(ionpot(106), ionpotP(1))
      real*8 jgS(16), jcS(16), ionpotS(16)
      equivalence(jg(121), jgS(1))
      equivalence(jc(121), jcS(1))
      equivalence(ionpot(121), ionpotS(1))
      real*8 jgCl(17), jcCl(17), ionpotCl(17)
      equivalence(jg(137), jgCl(1))
      equivalence(jc(137), jcCl(1))
      equivalence(ionpot(137), ionpotCl(1))
      real*8 jgAr(18), jcAr(18), ionpotAr(18)
      equivalence(jg(154), jgAr(1))
      equivalence(jc(154), jcAr(1))
      equivalence(ionpot(154), ionpotAr(1))
      real*8 jgK(19), jcK(19), ionpotK(19)
      equivalence(jg(172), jgK(1))
      equivalence(jc(172), jcK(1))
      equivalence(ionpot(172), ionpotK(1))
      real*8 jgCa(20), jcCa(20), ionpotCa(20)
      equivalence(jg(191), jgCa(1))
      equivalence(jc(191), jcCa(1))
      equivalence(ionpot(191), ionpotCa(1))
      real*8 jgSc(21), jcSc(21), ionpotSc(21)
      equivalence(jg(211), jgSc(1))
      equivalence(jc(211), jcSc(1))
      equivalence(ionpot(211), ionpotSc(1))
      real*8 jgTi(22), jcTi(22), ionpotTi(22)
      equivalence(jg(232), jgTi(1))
      equivalence(jc(232), jcTi(1))
      equivalence(ionpot(232), ionpotTi(1))
      real*8 jgV(23), jcV(23), ionpotV(23)
      equivalence(jg(254), jgV(1))
      equivalence(jc(254), jcV(1))
      equivalence(ionpot(254), ionpotV(1))
      real*8 jgCr(24), jcCr(24), ionpotCr(24)
      equivalence(jg(277), jgCr(1))
      equivalence(jc(277), jcCr(1))
      equivalence(ionpot(277), ionpotCr(1))
      real*8 jgMn(25), jcMn(25), ionpotMn(25)
      equivalence(jg(301), jgMn(1))
      equivalence(jc(301), jcMn(1))
      equivalence(ionpot(301), ionpotMn(1))
      real*8 jgFe(26), jcFe(26), ionpotFe(26)
      equivalence(jg(326), jgFe(1))
      equivalence(jc(326), jcFe(1))
      equivalence(ionpot(326), ionpotFe(1))
      real*8 jgCo(27), jcCo(27), ionpotCo(27)
      equivalence(jg(352), jgCo(1))
      equivalence(jc(352), jcCo(1))
      equivalence(ionpot(352), ionpotCo(1))
      real*8 jgNi(28), jcNi(28), ionpotNi(28)
      equivalence(jg(379), jgNi(1))
      equivalence(jc(379), jcNi(1))
      equivalence(ionpot(379), ionpotNi(1))
      real*8 jgCu(29), jcCu(29), ionpotCu(29)
      equivalence(jg(407), jgCu(1))
      equivalence(jc(407), jcCu(1))
      equivalence(ionpot(407), ionpotCu(1))
      real*8 jgZn(30), jcZn(30), ionpotZn(30)
      equivalence(jg(436), jgZn(1))
      equivalence(jc(436), jcZn(1))
      equivalence(ionpot(436), ionpotZn(1))
      
      data jgH/0.5/
      data jcH/0.5/
      data ionpotH/13.595/
      data ionpotHminus/0.754546971/  ! IonpotH * 2 * 0.0277509
      data jgHe/0.0, 0.5/
      data jcHe/0.0, 0.5/
      data ionpotHe/24.580, 54.4/
      data jgLi/0.5, 0.0, 0.5/
      data jcLi/0.5, 0.0, 0.5/
      data ionpotLi/5.390, 75.6193, 122.420/
      data jgBe/0.0, 0.5, 0.0, 0.5/
      data jcBe/0.0, 0.5, 0.0, 0.5/
      data ionpotBe/9.320, 18.206, 153.850, 217.657/
      data jgB/2.5, 0.0, 0.5, 0.0, 0.5/
      data jcB/2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotB/8.296, 25.149, 37.920, 259.298, 340.127/
      data jgC/4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcC/7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotC/11.264, 24.376, 47.864, 64.476, 391.986, 489.840/
      data jgN/1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcN/9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotN/14.540, 29.605, 47.426, 77.450, 97.863, 551.925,
     .     666.83/
      data jgO/4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcO/7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotO/13.614, 35.146, 54.934, 77.394, 113.873, 138.080,
     .     739.114, 871.12/
      data jgF/2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcF/2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotF/17.42, 34.98, 62.646, 87.23, 114.214, 157.117,
     .     185.139, 953.60, 1103.1/
      data jgNe/0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcNe/0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotNe/21.559, 41.07, 64.0, 97.16, 126.4, 157.91,
     .     207.3, 239.1, 1103.1, 1362.2/
      data jgNa/0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5,
     .     0.0, 0.5/
      data jcNa/0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5,
     .     0.0, 0.5/
      data ionpotNa/5.138, 47.29, 71.65, 98.9, 138.4, 172.1, 208.5,
     .     264.2, 299.9, 1465.1, 1648.7/
      data jgMg/0.0, 0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5,
     .     0.0, 0.5/
      data jcMg/0.0, 0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5,
     .     0.0, 0.5/
      data ionpotMg/7.644, 15.0, 80.1, 109.3, 141.3, 186.5, 225.0,
     .     265.9, 328.0, 367.4, 1761.6, 1962.6/
      data jgAl/2.5, 0.0, 0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0,
     .     0.5, 0.0, 0.5/
      data jcAl/2.5, 0.0, 0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0,
     .     0.5, 0.0, 0.5/
      data ionpotAl/5.984, 18.8, 28.4, 120.0, 153.8, 190.5, 241.4,
     .     284.6, 330.2, 398.6, 441.9, 2086.0, 2304.1/
      data jgSi/4.0, 2.5, 0.0, 0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5,
     .     0.0, 0.5, 0.0, 0.5/
      data jcSi/7.0, 2.5, 0.0, 0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5,
     .     0.0, 0.5, 0.0, 0.5/
      data ionpotSi/8.15, 16.3, 33.5, 45.1, 166.8, 205.2, 246.5,
     .     303.1, 351.1, 401.4, 476.1, 523.3, 2437.7, 2673.1/
      data jgP/1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 2.5, 4.0, 1.5, 4.0,
     .     2.5, 0.0, 0.5, 0.0, 0.5/
      data jcP/9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 2.5, 7.0, 9.5, 7.0,
     .     2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotP/10.49, 19.7, 30.2, 51.5, 65.0, 220.5, 263.4,
     .     309.3, 371.7, 424.4, 479.5, 560.5, 611.6, 2816.9, 3069.8/
      data jgS/4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 2.5, 4.0, 1.5,
     .     4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcS/7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 2.5, 7.0, 9.5,
     .     7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotS/10.36, 23.4, 35.0, 47.3, 72.7, 88.0, 281.0, 328.4,
     .     379.0, 447.1, 504.6, 564.5, 651.9, 706.8, 3223.8, 3494.0/
      data jgCl/2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 2.5, 4.0,
     .     1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcCl/2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 2.5, 7.0,
     .     9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotCl/12.97, 23.8, 39.9, 53.5, 67.6, 97.0, 114.2,
     .     348.5, 400.3, 455.6, 529.3, 591.7, 656.4, 750.1, 809.0,
     .     3658.0, 3946.0/
      data jgAr/0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 2.5,
     .     4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcAr/0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 2.5,
     .     7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotAr/15.76, 27.6, 40.9, 59.7, 75.2, 91.2, 124.5,
     .     143.4, 422.8, 479.1, 539.0, 618.5, 685.7, 755.3, 855.3,
     .     918.0, 4121.0, 4426.0/
      data jgK/0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0,
     .     2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcK/0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0,
     .     2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotK/4.339, 31.7, 45.8, 61.1, 87.2, 100.0, 117.8,
     .     155.1, 175.8, 503.9, 564.7, 629.4, 714.6, 786.6, 861.1,
     .     968.0, 1034.0, 4611.0, 4934.0/
      data jgCa/0.0, 0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5,
     .     0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcCa/0.0, 0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5,
     .     0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotCa/6.111, 11.87, 51.2, 67.3, 84.5, 108.8, 128.0,
     .     147.6, 188.8, 211.3, 591.9, 657.2, 726.6, 817.6, 894.5,
     .     974.0, 1087.0, 1157.0, 5129.0, 5470.0/
      data jgSc/4.5, 7.0, 4.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0,
     .     0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcSc/4.5, 9.5, 4.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0,
     .     0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotSc/6.56, 12.89, 24.75, 73.9, 91.9, 111.0, 138.0,
     .     159.0, 180.5, 225.6, 249.8, 686.8, 756.7, 830.8, 927.5,
     .     1009.0, 1094.0, 1213.0, 1288.0, 5675.0, 6034.0/
      data jgTi/10.0, 13.5, 10.0, 4.5, 0.0, 2.5, 4.0, 1.5, 4.0,
     .     2.5, 0.0, 0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5,
     .     0.0, 0.5/
      data jcTi/22.0, 44.5, 22.0, 4.5, 0.0, 2.5, 7.0, 9.5, 7.0,
     .     2.5, 0.0, 0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5,
     .     0.0, 0.5/
      data ionpotTi/6.83, 13.63, 28.14, 43.24, 99.9, 119.7, 140.7,
     .     170.4, 193.2, 216.5, 265.5, 291.6, 788.6, 863.1, 941.9,
     .     1044.0, 1131.0, 1221.0, 1346.0, 1425.0, 6249.0, 6626.0/
      data jgV/13.5, 12.0, 13.5, 10.0, 4.5, 0.0, 2.5, 4.0, 1.5,
     .     4.0, 2.5, 0.0, 0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0,
     .     0.5, 0.0, 0.5/
      data jcV/59.5, 104.5, 59.5, 22.0, 4.5, 0.0, 2.5, 7.0, 9.5,
     .     7.0, 2.5, 0.0, 0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0,
     .     0.5, 0.0, 0.5/
      data ionpotV/6.74, 14.2, 29.7, 48.0, 65.2, 129.0, 150.6, 173.4,
     .     205.8, 230.5, 255.7, 308.6, 336.4, 897.3, 976.0, 1060.0,
     .     1168.0, 1260.0, 1355.0, 1486.0, 1569.0, 6851.0, 7246.0/
      data jgCr/3.0, 2.5, 12.0, 13.5, 10.0, 4.5, 0.0, 2.5, 4.0,
     .     1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5,
     .     0.0, 0.5, 0.0, 0.5/
      data jcCr/251.5, 125.5, 104.5, 59.5, 22.0, 4.5, 0.0, 2.5,
     .     7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 2.5, 7.0, 9.5, 7.0,
     .     2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotCr/6.77, 16.49, 30.96, 49.16, 69.46, 90.64, 161.1,
     .     184.7, 209.3, 244.4, 270.8, 298.0, 354.8, 384.4, 1013.0,
     .     1097.0, 1185.0, 1299.0, 1396.0, 1496.0, 1634.0, 1721.0,
     .     7482.0, 7895.0/
      data jgMn/2.5, 3.0, 2.5, 12.0, 13.5, 10.0, 4.5, 0.0, 2.5,
     .     4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 2.5, 4.0, 1.5, 4.0,
     .     2.5, 0.0, 0.5, 0.0, 0.5/
      data jcMn/125.5, 251.5, 125.5, 104.5, 59.5, 22.0, 4.5, 0.0,
     .     2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 2.5, 7.0, 9.5,
     .     7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotMn/7.43, 15.64, 33.67, 51.2, 72.4, 95.6, 119.2,
     .     195.6, 221.8, 248.3, 286.0, 314.4, 343.6, 404.1, 436.0,
     .     1136.0, 1224.0, 1317.0, 1437.0, 1539.0, 1644.0, 1788.0,
     .     1880.0, 8141.0, 8572.0/
      data jgFe/12.0, 14.5, 12.0, 2.5, 12.0, 13.5, 10.0, 4.5, 0.0,
     .     2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 2.5, 4.0, 1.5,
     .     4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcFe/104.5, 209.5, 104.5, 125.5, 104.5, 59.5, 22.0, 4.5,
     .     0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 2.5, 7.0,
     .     9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotFe/7.90, 16.188, 30.652, 54.8, 75.0, 99.1, 124.98,
     .     151.061, 234.9, 262.0, 290.3, 330.8, 361.0, 392.0, 457.0,
     .     490.0, 1265.0, 1358.0, 1456.0, 1582.0, 1689.0, 1799.0,
     .     1950.0, 2045.0, 8828.0, 9278.0/
      data jgCo/13.5, 10.0, 13.5, 12.0, 2.5, 12.0, 13.5, 10.0, 4.5,
     .     0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 2.5, 4.0,
     .     1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcCo/59.5, 22.0, 59.5, 104.5, 125.5, 104.5, 59.5, 22.0,
     .     4.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 2.5,
     .     7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotCo/7.86, 17.083, 33.50, 51.3, 79.5, 102.0, 128.9,
     .     157.8, 186.13, 276.4, 305.0, 336.0, 379.0, 411.0, 444.0,
     .     512.0, 547.0, 1402.0, 1500.0, 1603.0, 1735.0, 1846.0,
     .     1962.0, 2119.0, 2218.0, 9544.0, 10012.2/
      data jgNi/10.0, 4.5, 10.0, 13.5, 12.0, 2.5, 12.0, 13.5, 10.0,
     .     4.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 2.5,
     .     4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcNi/22.0, 4.5, 22.0, 59.5, 104.5, 125.5, 104.5, 59.5,
     .     22.0, 4.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0,
     .     2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotNi/7.6375, 18.169, 35.19, 54.9, 76.06, 108.0, 133.0,
     .     162.0, 193.0, 224.6, 321.0, 352.0, 384.0, 430.0, 464.0,
     .     499.0, 571.0, 608.0, 1546.0, 1648.0, 1756.0, 1894.0,
     .     2011.0, 2131.0, 2295.0, 2399.0, 10288.8, 10775.48/
      data jgCu/0.5, 0.0, 4.5, 10.0, 13.5, 12.0, 2.5, 12.0, 13.5,
     .     10.0, 4.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0,
     .     2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcCu/0.5, 0.0, 4.5, 22.0, 59.5, 104.5, 125.5, 104.5,
     .     59.5, 22.0, 4.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5,
     .     0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotCu/7.724, 20.29, 36.83, 55.2, 79.9, 103., 139.,
     .     166., 199., 232., 266., 368.8, 401.0, 435.0, 484.0, 520.0,
     .     557.0, 633.0, 672.0, 1697.0, 1804.0, 1916.0, 2060.0,
     .     2182.0, 2308.0, 2478.0, 2573.0, 11036.0, 11561.0/
      data jgZn/0.0, 0.5, 0.0, 4.5, 10.0, 13.5, 12.0, 2.5, 12.0,
     .     13.5, 10.0, 4.5, 0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5,
     .     0.0, 2.5, 4.0, 1.5, 4.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data jcZn/0.0, 0.5, 0.0, 4.5, 22.0, 59.5, 104.5, 125.5,
     .     104.5, 59.5, 22.0, 4.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0,
     .     0.5, 0.0, 2.5, 7.0, 9.5, 7.0, 2.5, 0.0, 0.5, 0.0, 0.5/
      data ionpotZn/9.391, 17.96, 39.722, 59.4, 82.6, 108.0, 134.0,
     .     174.0, 203.0, 238.0, 274.0, 310.8, 419.7, 454.0, 490.0,
     .     542.0, 579.0, 619.0, 698.0, 739.0, 1854.0, 1966.0, 2084.0,
     .     2234.0, 2361.0, 2493.0, 2652.0, 2754.0, 11741.0, 12370.0/


      data pi/3.141592653589793d+00/, fourpi/12.5637061d+00/
      data bc/1.380626d-16/, h/6.626205d-27/, hev/4.1357d-15/
      data c/2.997925d+10/, elecxsec/6.6524d-25/
      data evtoerg/1.6022d-12/, a/7.56464d-15/
      data stefbltz/5.66956d-05/, hmass/1.67352d-24/
      data esu/4.80298d-10/, emass/9.1091d-28/
      data srpi/1.772453851d+00/, bcev/8.617064d-05/

      data init/1/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (init .eq. 1) then
         sahacoef = 2.d+00 * ((emass / h) * (bc / h) * 2.d+00 * pi)**1.5
     ~
         init = 0
      end if
C add bakl  error here -- sahacoef !=0 bakl 4.828724981518439E+015
C add bakl  sahacoef = 2.d+00 * ((emass / h) * (bc / h) * 2.d+00 * pi)**1.5
	 sahacoef = 4.828724981518439d+015

      saha = sahacoef * temp**1.5 / eden
      bctev = bcev * temp
      bctevinv = 1.d+00 / bctev
      smallnum = dmach(2)

      do i = 1, min(maxion, z+1)
         fraction(i) = 0.d+00
         partition(i) = 1.d+00
      end do

      if (z .gt. 30) then
         firstion = 1
         call pfsaha(z, min(maxion, 6), 12, temp,
     .        press, eden, fraction(1:6))
         call pfsaha(z, min(maxion, 6), 13, temp,
     .        press, eden, partition(1:6))
         return
      end if

      frac(1) = 1.d+00
      ipos = (z * (z - 1)) / 2

      if (kurucz) then
         call pfsaha(z, 1, 3, temp, press, eden, partit(1))
      else if (tempflag .eq. 0) then
         partit(1) = 2.d+00 * jg(ipos+1) + 1.d+00
      else
         partit(1) = 2.d+00 * jc(ipos+1) + 1.d+00
      end if

      do i = 2, z + 1
         if (i .lt. (z + 1)) then
            if (kurucz .and. (i .lt. 7)) then
               call pfsaha(z, i, 3, temp, press, eden, partit(i))
            else if (tempflag .eq. 0) then
               partit(i) = 2.d+00 * jg(ipos+i) + 1.d+00
            else
               partit(i) = 2.d+00 * jc(ipos+i) + 1.d+00
            end if
         else
            partit(i) = 1.d+00
         end if

         frac(i) = saha * dexp(-ionpot(ipos+i-1) * bctevinv) *
     .        partit(i) / partit(i-1)
c      if(z.le.3)write(0,*)' saha bctevinv partit i i-1=', saha,
c     + bctevinv, partit(i),partit(i-1)
!         if(power .lt. 1.d0/sqrt(smallnum))then
!           power = power * frac(i)
!         else
!           alogpow=min(log(power)+log(frac(i)),-log(smallnum))
!           power = exp(alogpow)
!         endif
         if (frac(i) .lt. smallnum) go to 10
      end do

 10   icutoff = i - 1

c      if(z.le.3)write(0,*)' i  fract=', frac(i),i
      imax = 1

      do i = 2, icutoff
         if (frac(i) .gt. 1.d+00) imax = i
      end do

      sum = 1.d+00
      power = 1.d+00

      do i = imax + 1, icutoff
         power = power * frac(i)
         frac(i) = power
         sum = sum + frac(i)
      end do

      do i = 1, imax - 1
         frac(i) = 1.d+00 / frac(i+1)
      end do

      power = 1.d+00

      do k = 1, imax - 1
         i = imax - k
         power = power * frac(i)
         frac(i) = power
         sum = sum + frac(i)
      end do

      frac(imax) = 1.d+00 / sum

      do i = 1, imax - 1
         frac(i) = frac(i) * frac(imax)
      end do

      do i = imax + 1, icutoff
         frac(i) = frac(i) * frac(imax)
      end do

      if (icutoff .gt. maxion) then
         imax = 1

         do i = 2, icutoff
            if (frac(i) .gt. frac(imax)) imax = i
         end do

         i1 = imax
         i2 = imax

         do k = 2, min(maxion, (z+1))
            if (i1 .gt. 1) then
               if (i2 .lt. (z+1)) then
                  if (frac(i1-1) .gt. frac(i2+1)) then
                     i1 = i1 - 1
                  else
                     i2 = i2 + 1
                  end if
               else
                  i1 = i1 - 1
               end if
            else
               i2 = i2 + 1
            end if
         end do

         firstion = i1
      else
         firstion = 1
      end if

      do i = firstion, min(firstion + maxion - 1, icutoff)
         fraction(i - firstion + 1) = frac(i)
         partition(i - firstion + 1) = partit(i)
      end do

!  For negative hydrogen ion      
      if( z .eq. 1 ) then
       fraction(0) = dexp( ionpotHminus(1) * bctevinv ) / saha * 0.5d0
      endif
      return

      end
