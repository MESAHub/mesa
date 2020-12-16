C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: exchange_coeff.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine exchange_coeff(index_exchange, maxc, maxkind,
     &  ccoeff, mforder, mgorder)
C      subroutine to deliver J and K exchange integral coefficients
C      calculated from EFF-style fit with non-relativistic limit
C      subtracted as described in Paper IV.
C      index_exchange gives a choice of approximations
C      to the J and K=sqrt(I) integrals which represent a variety of
C      compromises between accuracy and speed (number of coefficients).
C
C      All sets of
C      coefficients have been determined using j_fit and sqrti_fit
C       in my utilities programmes with results stored in
C      utils/eff_results/j_fit_fortran.out_new23 and
C      utils/eff_results/sqrti_fit_fortran.out18
C      The 7-figure rounding should introduce relative errors of 1
C      part in 10^7 at most which is far below the fitting errors.
C
      implicit none
      integer maxc, maxkind, index_exchange,
     &  mforder(maxkind), mgorder(maxkind),
     &  j1mf, j1mg, j1mflg, j1mglg, j1mflg2, j1mglg2, j1mflf, j1mglf,
     &  j2mf, j2mg, j2mflg, j2mglg, j2mflg2, j2mglg2, j2mflf, j2mglf,
     &  j3mf, j3mg, j3mflg, j3mglg, j3mflg2, j3mglg2, j3mflf, j3mglf,
     &  k1mf, k1mg, k1mflg, k1mglg,
     &  k2mf, k2mg, k2mflg, k2mglg,
     &  k3mf, k3mg, k3mflg, k3mglg,
     &  i, j, m, ioffset
      integer j1dim, j2dim, j3dim
      parameter (j1dim=3*3+3*3+0*0+3*3)
      parameter (j2dim=3*3+5*5+0*0+5*5)
      parameter (j3dim=4*4+7*7+0*0+7*7)
      integer k1dim, k2dim, k3dim
      parameter (k1dim=3*3+1)
      parameter (k2dim=5*5+2)
      parameter (k3dim=8*8+2)
C       n.b. ccoeff in order of J, Jlng, Jlng2, Jlnf, K, and Klng.
      double precision ccoeff(maxc,maxkind),
     &  j1coeff(j1dim),j2coeff(j2dim), j3coeff(j3dim),
     &  k1coeff(k1dim),k2coeff(k2dim), k3coeff(k3dim)

C     Coefficients for mf, mg, mflg, mglg, mflg2, mglg2, mflf, mglf = 
C         3    3    3    3    0    0    3    3
C     with maximum relative residual =       0.0039994
      data j1mf, j1mg, j1mflg, j1mglg, j1mflg2, j1mglg2, j1mflf, j1mglf
     &/    3,    3,    3,    3,    0,    0,    3,    3/
      data j1coeff/
     &       34.8358213D0,       86.1918396D0,       51.5419186D0,
     &        6.1265167D0,       25.7316016D0,       29.3774454D0,
     &        0.8600275D0,        5.6471122D0,        7.0535536D0,
     &       33.2562696D0,       66.7030494D0,       28.4278181D0,
     &       31.5537757D0,       65.6524569D0,       30.6499630D0,
     &        6.8170355D0,       14.8834241D0,        8.0074649D0,
     &        7.6128464D0,       16.9077576D0,        7.1078540D0,
     &        2.7809264D0,        8.2810716D0,        9.0118836D0,
     &        2.4421025D0,        6.7537105D0,        1.1701926D0/

C     Coefficients for mf, mg, mflg, mglg, mflg2, mglg2, mflf, mglf = 
C         3    3    5    5    0    0    5    5
C     with maximum relative residual =       0.0005652
      data j2mf, j2mg, j2mflg, j2mglg, j2mflg2, j2mglg2, j2mflf, j2mglf
     &/    3,    3,    5,    5,    0,    0,    5,    5/
      data j2coeff/
     &       56.2805336D0,      155.6222406D0,       83.5532182D0,
     &        1.4600744D0,      -38.3239264D0,       28.8234751D0,
     &        0.8040612D0,       10.8897054D0,        7.0883559D0,
     &        5.6437776D0,      174.7225632D0,     1174.0556985D0,
     &       69.7350999D0,      -19.1101716D0,       72.0785970D0,
     &      700.7413357D0,     1553.6277210D0,      461.3841453D0,
     &       32.3236302D0,      105.3785812D0,      598.6273693D0,
     &      854.7360263D0,      494.9380271D0,       87.5200915D0,
     &       46.7695962D0,      198.3020460D0,      291.6944973D0,
     &      200.3802976D0,       47.0467342D0,        6.8234346D0,
     &       28.4677241D0,       43.6666247D0,       31.6890797D0,
     &        8.0006636D0,       28.6410004D0,      108.8331950D0,
     &      144.3375035D0,       94.8689197D0,       15.7800513D0,
     &      165.6488584D0,     1460.9695944D0,     2059.8613336D0,
     &      649.1216512D0,       43.7140906D0,      274.3616284D0,
     &      848.1550036D0,      932.1227056D0,      513.4043189D0,
     &       40.9849434D0,       -7.8112104D0,      -24.4941518D0,
     &      -37.9448618D0,      -49.7713793D0,       10.9901399D0,
     &        7.6299746D0,       23.9006008D0,       32.1617806D0,
     &       31.8665472D0,        3.4726043D0/

C     Coefficients for mf, mg, mflg, mglg, mflg2, mglg2, mflf, mglf = 
C         4    4    7    7    0    0    7    7
C     with maximum relative residual =       0.0000561
      data j3mf, j3mg, j3mflg, j3mglg, j3mflg2, j3mglg2, j3mflf, j3mglf
     &/    4,    4,    7,    7,    0,    0,    7,    7/
      data j3coeff/
     &       77.7223322D0,      329.9690751D0,      339.1739022D0,
     &      115.5562291D0,       15.1130448D0,    -1711.5422468D0,
     &      778.5679695D0,       74.2173982D0,       -2.6865083D0,
     &     1199.9753651D0,       54.6086181D0,       31.4786125D0,
     &        0.7914812D0,      -73.9603649D0,       14.3235218D0,
     &        7.0902001D0,       -5.4000121D0,     5695.7957148D0,
     &    11216.4052787D0,    -5004.0576957D0,    -2978.6276336D0,
     &     -990.8463971D0,      -61.0768342D0,      146.7136763D0,
     &    -7092.1165083D0,    28977.5973434D0,   -17432.3188130D0,
     &    -6011.6709017D0,    -1902.4091931D0,       -3.5689832D0,
     &      411.6203402D0,   -13101.9904352D0,    27389.3394169D0,
     &    -4881.7422341D0,    -4714.8750809D0,    -1278.8279881D0,
     &      243.7667578D0,      441.1299352D0,    -2447.2072200D0,
     &    13360.4808816D0,     7561.9107837D0,      776.1294503D0,
     &      974.1793642D0,      349.9756271D0,      232.4649446D0,
     &     1014.6845310D0,     4202.0496258D0,     4667.0296065D0,
     &     2627.0527294D0,     1238.7406226D0,      213.4913981D0,
     &       61.5101678D0,      368.7097096D0,      953.5922014D0,
     &     1274.5396464D0,      942.4763969D0,      392.2869249D0,
     &       63.9909089D0,        6.8247249D0,       42.0150160D0,
     &      107.4877403D0,      149.2296305D0,      113.1130650D0,
     &       48.0123563D0,        8.0000291D0,       76.9669717D0,
     &      403.8465365D0,      877.1286674D0,     1009.2265882D0,
     &      625.3080827D0,      215.1119554D0,       23.7703342D0,
     &     4196.5960409D0,    13947.6028263D0,     6386.7560942D0,
     &    -8509.1033779D0,    -5068.4236377D0,     -116.6120752D0,
     &       60.5395737D0,   -14690.7086472D0,   -10481.5778168D0,
     &   -10134.5858478D0,   -31981.8034156D0,   -14050.4153207D0,
     &      220.8791890D0,      138.3380767D0,   -11483.9761233D0,
     &    -8157.4031117D0,     4139.5945610D0,    -2456.0463932D0,
     &    -4921.9776737D0,     -980.5436208D0,       49.8463981D0,
     &      273.9069680D0,     7797.1164005D0,    19163.2886021D0,
     &    15475.1207915D0,     5080.5195211D0,     1423.8520116D0,
     &      107.3457476D0,      925.7506712D0,     4392.2634040D0,
     &     7721.2581060D0,     6649.2247613D0,     2808.9762429D0,
     &      365.3071337D0,       20.1232319D0,      -78.6975448D0,
     &     -351.4769520D0,     -595.0931958D0,     -468.7739531D0,
     &     -177.7836871D0,        4.6112748D0,        3.3403911D0/

C     Coefficients for mf, mg, mflg, mglg =     3    3    1    1
C     with maximum relative residual =       0.0030394
      data k1mf, k1mg, k1mflg, k1mglg
     &/    3,    3,    1,    1/
      data k1coeff/
     &        9.4150106D0,       19.3949179D0,       10.5602307D0,
     &        9.0399638D0,       18.4548660D0,        9.8746980D0,
     &        2.6135398D0,        5.2690297D0,        2.8270724D0,
     &        0.2219891D0/
C     Coefficients for mf, mg, mflg, mglg =     5    5    2    1
C     with maximum relative residual =       0.0002268
      data k2mf, k2mg, k2mflg, k2mglg
     &/    5,    5,    2,    1/
      data k2coeff/
     &       15.9616640D0,       65.2137462D0,      100.0779497D0,
     &       68.9805226D0,       18.1027203D0,       31.6264744D0,
     &      128.5287114D0,      196.2665661D0,      134.3128503D0,
     &       35.0222537D0,       30.8574217D0,      124.7941936D0,
     &      189.7875295D0,      128.9954761D0,       33.5069500D0,
     &       14.3505923D0,       58.0368244D0,       88.3645023D0,
     &       59.8503499D0,       15.5482948D0,        2.6124722D0,
     &       10.5422097D0,       16.0646752D0,       10.8120877D0,
     &        2.8284087D0,        0.7085843D0,        0.4384647D0/
C     Coefficients for mf, mg, mflg, mglg =     8    8    2    1
C     with maximum relative residual =       0.0000309
      data k3mf, k3mg, k3mflg, k3mglg
     &/    8,    8,    2,    1/
      data k3coeff/
     &       25.7842747D0,      182.7491698D0,      556.4181655D0,
     &      938.5331301D0,      960.7519007D0,      583.5117994D0,
     &      200.4860761D0,       29.4157259D0,       89.6532756D0,
     &      634.2399245D0,     1927.2484177D0,     3244.9585385D0,
     &     3313.2266232D0,     2009.3901546D0,      688.3557251D0,
     &      100.8001430D0,      178.2531213D0,     1258.5795649D0,
     &     3816.1015364D0,     6413.7471264D0,     6532.0560019D0,
     &     3952.7213124D0,     1350.3664972D0,      197.3155508D0,
     &      220.2391122D0,     1552.0338194D0,     4696.8562020D0,
     &     7879.1826389D0,     8002.6234013D0,     4839.3366977D0,
     &     1647.8294240D0,      240.4593766D0,      171.7202945D0,
     &     1209.2423616D0,     3655.4138590D0,     6129.2942522D0,
     &     6214.1077832D0,     3753.8298340D0,     1275.4438080D0,
     &      185.9307018D0,       82.5677241D0,      580.9794546D0,
     &     1754.9056233D0,     2940.4105339D0,     2977.6103202D0,
     &     1798.9814958D0,      610.3316895D0,       89.0447118D0,
     &       22.2022652D0,      156.3493557D0,      472.5034538D0,
     &      792.5785975D0,      802.1268594D0,      485.3881563D0,
     &      164.4720854D0,       24.0404251D0,        2.6124292D0,
     &       18.3854785D0,       55.5290929D0,       93.0966512D0,
     &       94.1489718D0,       56.9861334D0,       19.2907139D0,
     &        2.8284359D0,        1.0175489D0,        0.5914964D0/
      save
C      checks
      if(j1dim.ne.j1mf*j1mg+j1mflg*j1mglg+j1mflg2*j1mglg2+j1mflf*j1mglf)
     &  stop 'exchange_coeff: inconsistent j1 dimension data'
      if(j2dim.ne.j2mf*j2mg+j2mflg*j2mglg+j2mflg2*j2mglg2+j2mflf*j2mglf)
     &  stop 'exchange_coeff: inconsistent j2 dimension data'
      if(j3dim.ne.j3mf*j3mg+j3mflg*j3mglg+j3mflg2*j3mglg2+j3mflf*j3mglf)
     &  stop 'exchange_coeff: inconsistent j3 dimension data'
      if(k1dim.ne.k1mf*k1mg+k1mflg*k1mglg)
     &  stop 'exchange_coeff: inconsistent k1 dimension data'
      if(k2dim.ne.k2mf*k2mg+k2mflg*k2mglg)
     &  stop 'exchange_coeff: inconsistent k2 dimension data'
      if(k3dim.ne.k3mf*k3mg+k3mflg*k3mglg)
     &  stop 'exchange_coeff: inconsistent k3 dimension data'
      do j = 1, maxkind
        if(index_exchange.eq.1) then
          if(j.eq.1) then
            mforder(j) = j1mf-1
            mgorder(j) = j1mg-1
            m = j1mf*j1mg
            ioffset = 0
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j1coeff(i+ioffset)
            enddo
          elseif(j.eq.2) then
            mforder(j) = j1mflg-1
            mgorder(j) = j1mglg-1
            m = j1mflg*j1mglg
            ioffset = j1mf*j1mg
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j1coeff(i+ioffset)
            enddo
          elseif(j.eq.3) then
            mforder(j) = j1mflg2-1
            mgorder(j) = j1mglg2-1
            m = j1mflg2*j1mglg2
            ioffset = j1mf*j1mg+j1mflg*j1mglg
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j1coeff(i+ioffset)
            enddo
          elseif(j.eq.4) then
            mforder(j) = j1mflf-1
            mgorder(j) = j1mglf-1
            m = j1mflf*j1mglf
            ioffset = j1mf*j1mg+j1mflg*j1mglg+j1mflg2*j1mglg2
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j1coeff(i+ioffset)
            enddo
          elseif(j.eq.5) then
            mforder(j) = k1mf-1
            mgorder(j) = k1mg-1
            m = k1mf*k1mg
            ioffset = 0
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = k1coeff(i+ioffset)
            enddo
          elseif(j.eq.6) then
            mforder(j) = k1mflg-1
            mgorder(j) = k1mglg-1
            m = k1mflg*k1mglg
            ioffset = k1mf*k1mg
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = k1coeff(i+ioffset)
            enddo
          endif
        elseif(index_exchange.eq.2) then
          if(j.eq.1) then
            mforder(j) = j2mf-1
            mgorder(j) = j2mg-1
            m = j2mf*j2mg
            ioffset = 0
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j2coeff(i+ioffset)
            enddo
          elseif(j.eq.2) then
            mforder(j) = j2mflg-1
            mgorder(j) = j2mglg-1
            m = j2mflg*j2mglg
            ioffset = j2mf*j2mg
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j2coeff(i+ioffset)
            enddo
          elseif(j.eq.3) then
            mforder(j) = j2mflg2-1
            mgorder(j) = j2mglg2-1
            m = j2mflg2*j2mglg2
            ioffset = j2mf*j2mg+j2mflg*j2mglg
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j2coeff(i+ioffset)
            enddo
          elseif(j.eq.4) then
            mforder(j) = j2mflf-1
            mgorder(j) = j2mglf-1
            m = j2mflf*j2mglf
            ioffset = j2mf*j2mg+j2mflg*j2mglg+j2mflg2*j2mglg2
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j2coeff(i+ioffset)
            enddo
          elseif(j.eq.5) then
            mforder(j) = k2mf-1
            mgorder(j) = k2mg-1
            m = k2mf*k2mg
            ioffset = 0
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = k2coeff(i+ioffset)
            enddo
          elseif(j.eq.6) then
            mforder(j) = k2mflg-1
            mgorder(j) = k2mglg-1
            m = k2mflg*k2mglg
            ioffset = k2mf*k2mg
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = k2coeff(i+ioffset)
            enddo
          endif
        elseif(index_exchange.eq.3) then
          if(j.eq.1) then
            mforder(j) = j3mf-1
            mgorder(j) = j3mg-1
            m = j3mf*j3mg
            ioffset = 0
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j3coeff(i+ioffset)
            enddo
          elseif(j.eq.2) then
            mforder(j) = j3mflg-1
            mgorder(j) = j3mglg-1
            m = j3mflg*j3mglg
            ioffset = j3mf*j3mg
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j3coeff(i+ioffset)
            enddo
          elseif(j.eq.3) then
            mforder(j) = j3mflg2-1
            mgorder(j) = j3mglg2-1
            m = j3mflg2*j3mglg2
            ioffset = j3mf*j3mg+j3mflg*j3mglg
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j3coeff(i+ioffset)
            enddo
          elseif(j.eq.4) then
            mforder(j) = j3mflf-1
            mgorder(j) = j3mglf-1
            m = j3mflf*j3mglf
            ioffset = j3mf*j3mg+j3mflg*j3mglg+j3mflg2*j3mglg2
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = j3coeff(i+ioffset)
            enddo
          elseif(j.eq.5) then
            mforder(j) = k3mf-1
            mgorder(j) = k3mg-1
            m = k3mf*k3mg
            ioffset = 0
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = k3coeff(i+ioffset)
            enddo
          elseif(j.eq.6) then
            mforder(j) = k3mflg-1
            mgorder(j) = k3mglg-1
            m = k3mflg*k3mglg
            ioffset = k3mf*k3mg
            if(m.gt.maxc) stop 'exchange_coeff: input maxc too small'
            do i = 1, m
              ccoeff(i,j) = k3coeff(i+ioffset)
            enddo
          endif
        else
          stop 'exchange_coeff: invalid ifexchange'
        endif
      enddo
      end
