C
C                          CS2TST
C                         11/20/98
C
C   This program computes interpolation errors using the
C scattered data package CSHEP2D for each of ten test
C functions and a 33 by 33 uniform grid of interpolation
C points in the unit square.
C
C   This program uses Subroutines TESTDT and TSTFN1 from
C ACM Algorithm SURVEY to generate a node set and and the
C test function values.
C



      SUBROUTINE TEST_RENKA790_SG
      USE interp_2d_lib_sg

      INTEGER   NMAX, NRMAX, NI
      PARAMETER (NMAX=100, NRMAX=10, NI=33)
C
C Array storage:
C
      REAL X(NMAX), Y(NMAX), W(NMAX), RW(NMAX),
     .                 A(9,NMAX), P(NI), FT(NI,NI)
      INTEGER          LCELL(NRMAX,NRMAX), LNEXT(NMAX)
C
      REAL DEL, DUM, DX, DY, ERMAX, ERMEAN, PW,
     .                 RMAX, SSA, SSE, SSM, SUM, XMIN, YMIN
C      REAL CS2VAL_sg
      INTEGER          I, IER, J, K, KF, KFF, KFL, KS, LOUT,
     .                 N, NC, NFUN, NP, NR, NSET, NW, NWMAX
C
C Data:
C
C LOUT = Logical unit number for output.
C NSET = Number of node sets.
C NFUN = Number of test functions.
C
      DATA LOUT/6/, NSET/5/, NFUN/10/
      
C Specify test parameters
C      
      KS = 4
      KFF = 1
      KFL = 10
      NC = 17
      NW = 30
      NR =  6
      CALL TESTDT_sg (KS, N,X,Y)

C
C Set up uniform grid points.
C
      DEL = 1./DBLE(NI-1)
      DO 5 I = 1,NI
        P(I) = DBLE(I-1)*DEL
    5   CONTINUE
C
C Initialize the average SSE/SSM value to zero.
C
      SSA = 0.
C
C Print a heading and loop on test functions.
C
      WRITE (LOUT,200) KS, N, NI, NC, NW, NR
      DO 11 KF = KFF,KFL
C
C   Compute true function values at the nodes.
C
        DO 6 K = 1,N
          CALL TSTFN1_sg (KF,X(K),Y(K),0, W(K),DUM,DUM)
    6     CONTINUE
C
C   Compute true function values FT on the uniform grid, and
C     accumulate the sum of values SUM and sum of squared
C     values SSM.
C
        SUM = 0.
        SSM = 0.
        DO 8 I = 1,NI
          DO 7 J = 1,NI
            CALL TSTFN1_sg (KF,P(I),P(J),0, FT(I,J),DUM,DUM)
            SUM = SUM + FT(I,J)
            SSM = SSM + FT(I,J)**2
    7       CONTINUE
    8     CONTINUE
C
C   Compute the sum of squared deviations from the mean SSM.
C
        SSM = SSM - SUM*SUM/DBLE(NI*NI)
C
C   Compute parameters A and RW defining the interpolant.
C
        CALL interp_CSHEP2_sg (N,X,Y,W,NC,NW,NR, LCELL,LNEXT,XMIN,
     .               YMIN,DX,DY,RMAX,RW,A,IER)
        IF (IER .NE. 0) GO TO 21
C
C   Compute interpolation errors.
C
        ERMEAN = 0.
        ERMAX = 0.
        SSE = 0.
        DO 10 I = 1,NI
          DO 9 J = 1,NI
            IER = 0
            PW = interp_CS2VAL_sg (P(I),P(J),N,X,Y,W,NR,LCELL,LNEXT,
     .                   XMIN,YMIN,DX,DY,RMAX,RW,A,IER) -
     .           FT(I,J)
            IF (IER .NE. 0) THEN
               WRITE(LOUT,*) 'IER nonzero from CS2VAL_sg'
               STOP 1
            END IF
            ERMEAN = ERMEAN + ABS(PW)
            ERMAX = MAX(ERMAX,ABS(PW))
            SSE = SSE + PW*PW
    9       CONTINUE
   10     CONTINUE
        NP = NI*NI
        ERMEAN = ERMEAN/DBLE(NP)
        SSE = SSE/SSM
        SSA = SSA + SSE
        WRITE (LOUT,210) KF, ERMAX, ERMEAN, SSE
   11   CONTINUE
C
C Print the average SSE/SSM value (averaged over the test
C   functions).
C
      RETURN
C
C N is outside its valid range.
C
   20 WRITE (LOUT,500) N, NMAX
      STOP 1
C
C Error in CSHEP2.
C
   21 IF (IER .EQ. 2) WRITE (LOUT,510)
      IF (IER .EQ. 3) WRITE (LOUT,520)
      STOP 1
C
C Print formats:
C
  200 FORMAT ('RENKA790_sg: Node set ',I2,4X,'N =',I4,4X,'NI = ',I2,
     .        4X,'NC = ',I2,4X,'NW = ',I2,4X,'NR = ',I2/
     .        1X,16X,'Function',4X,'Max Error',4X,
     .        'Mean Error',4X,'SSE/SSM')
  210 FORMAT (1X,19X,I2,9X,F7.4,6X,F8.5,2X,F9.6)
C
C Error message formats:
C
  500 FORMAT (///1X,10X,'*** Error in data -- N = ',I4,
     .        ', Maximum value =',I4,' ***')
  510 FORMAT (///1X,14X,'*** Error in CSHEP2 -- duplicate ',
     .        'nodes encountered ***')
  520 FORMAT (///1X,14X,'*** Error in CSHEP2 -- all nodes ',
     .        'are collinear ***')
      END
      
      
      
      
      SUBROUTINE TESTDT_sg (K, N,X,Y)
      real X(100), Y(100)
      INTEGER K, N
C
C***********************************************************
C
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/28/97
C
C   This subroutine returns one of five sets of nodes used
C for testing scattered data fitting methods.  All five sets
C approximately cover the unit square [0,1] X [0,1]:  the
C convex hulls of sets 1 and 3 extend slightly outside the
C square but do not completely cover it, those of sets 2 and
C 5 coincide with the unit square, and the convex hull of
C set 4 is a large subset of the unit square.
C
C On input:
C
C       K = Integer in the range 1 to 5 which determines the
C           choice of data set as follows:
C
C               K = 1 - Franke's 100-node set
C               K = 2 - Franke's 33-node set
C               K = 3 - Lawson's 25-node set
C               K = 4 - Random 100-node set
C               K = 5 - Gridded 81-node set
C
C       X,Y = Arrays of length at least N(K), where
C             N(1) = 100, N(2) = 33, N(3) = 25,
C             N(4) = 100, and N(5) = 81.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       N = Number of nodes in set K, or 0 if K is outside
C           its valid range.
C
C       X,Y = Nodal coordinates of node set K.
C
C Subprograms required by TESTDT:  None
C
C***********************************************************
C
      real X1(100), Y1(100),  X2(33), Y2(33),
     .                 X3(25), Y3(25),  X4(100), Y4(100),
     .                 X5(81), Y5(81)
      INTEGER I
C
C Node set 1:  Franke's 100-node set.
C
      DATA (X1(I),Y1(I), I = 1,20)/
     .      0.0227035, -0.0310206,  0.0539888,  0.1586742,
     .      0.0217008,  0.2576924,  0.0175129,  0.3414014,
     .      0.0019029,  0.4943596, -0.0509685,  0.5782854,
     .      0.0395408,  0.6993418, -0.0487061,  0.7470194,
     .      0.0315828,  0.9107649, -0.0418785,  0.9962890,
     .      0.1324189,  0.0501330,  0.1090271,  0.0918555,
     .      0.1254439,  0.2592973,  0.0934540,  0.3381592,
     .      0.0767578,  0.4171125,  0.1451874,  0.5615563,
     .      0.0626494,  0.6552235,  0.1452734,  0.7524066,
     .      0.0958668,  0.9146523,  0.0695559,  0.9632421/
      DATA (X1(I),Y1(I), I = 21,40)/
     .      0.2645602,  0.0292939,  0.2391645,  0.0602303,
     .      0.2088990,  0.2668783,  0.2767329,  0.3696044,
     .      0.1714726,  0.4801738,  0.2266781,  0.5940595,
     .      0.1909212,  0.6878797,  0.1867647,  0.8185576,
     .      0.2304634,  0.9046507,  0.2426219,  0.9805412,
     .      0.3663168,  0.0396955,  0.3857662,  0.0684484,
     .      0.3832392,  0.2389548,  0.3179087,  0.3124129,
     .      0.3466321,  0.4902989,  0.3776591,  0.5199303,
     .      0.3873159,  0.6445227,  0.3812917,  0.8203789,
     .      0.3795364,  0.8938079,  0.2803515,  0.9711719/
      DATA (X1(I),Y1(I), I = 41,60)/
     .      0.4149771, -0.0284618,  0.4277679,  0.1560965,
     .      0.4200010,  0.2262471,  0.4663631,  0.3175094,
     .      0.4855658,  0.3891417,  0.4092026,  0.5084949,
     .      0.4792578,  0.6324247,  0.4812279,  0.7511007,
     .      0.3977761,  0.8489712,  0.4027321,  0.9978728,
     .      0.5848691, -0.0271948,  0.5730076,  0.1272430,
     .      0.6063893,  0.2709269,  0.5013894,  0.3477728,
     .      0.5741311,  0.4259422,  0.6106955,  0.6084711,
     .      0.5990105,  0.6733781,  0.5380621,  0.7235242,
     .      0.6096967,  0.9242411,  0.5026188,  1.0308762/
      DATA (X1(I),Y1(I), I = 61,80)/
     .      0.6616928,  0.0255959,  0.6427836,  0.0707835,
     .      0.6396475,  0.2008336,  0.6703963,  0.3259843,
     .      0.7001181,  0.4890704,  0.6333590,  0.5096324,
     .      0.6908947,  0.6697880,  0.6895638,  0.7759569,
     .      0.6718889,  0.9366096,  0.6837675,  1.0064516,
     .      0.7736939,  0.0285374,  0.7635332,  0.1021403,
     .      0.7410424,  0.1936581,  0.8258981,  0.3235775,
     .      0.7306034,  0.4714228,  0.8086609,  0.6091595,
     .      0.8214531,  0.6685053,  0.7290640,  0.8022808,
     .      0.8076643,  0.8476790,  0.8170951,  1.0512371/
      DATA (X1(I),Y1(I), I = 81,100)/
     .      0.8424572,  0.0380499,  0.8684053,  0.0902048,
     .      0.8366923,  0.2083092,  0.9418461,  0.3318491,
     .      0.8478122,  0.4335632,  0.8599583,  0.5910139,
     .      0.9175700,  0.6307383,  0.8596328,  0.8144841,
     .      0.9279871,  0.9042310,  0.8512805,  0.9696030,
     .      1.0449820, -0.0120900,  0.9670631,  0.1334114,
     .      0.9857884,  0.2695844,  0.9676313,  0.3795281,
     .      1.0129299,  0.4396054,  0.9657040,  0.5044425,
     .      1.0019855,  0.6941519,  1.0359297,  0.7459923,
     .      1.0414677,  0.8682081,  0.9471506,  0.9801409/
C
C Node set 2:  Franke's 33-node set.
C
      DATA (X2(I),Y2(I), I = 1,33)/
     .      0.05,  0.45,  0.00,  0.50,
     .      0.00,  1.00,  0.00,  0.00,
     .      0.10,  0.15,  0.10,  0.75,
     .      0.15,  0.30,  0.20,  0.10,
     .      0.25,  0.20,  0.30,  0.35,
     .      0.35,  0.85,  0.50,  0.00,
     .      0.50,  1.00,  0.55,  0.95,
     .      0.60,  0.25,  0.60,  0.65,
     .      0.60,  0.85,  0.65,  0.70,
     .      0.70,  0.20,  0.70,  0.65,
     .      0.70,  0.90,  0.75,  0.10,
     .      0.75,  0.35,  0.75,  0.85,
     .      0.80,  0.40,  0.80,  0.65,
     .      0.85,  0.25,  0.90,  0.35,
     .      0.90,  0.80,  0.95,  0.90,
     .      1.00,  0.00,  1.00,  0.50,
     .      1.00,  1.00/
C
C Node set 3:  Lawson's 25-node set.
C
      DATA (X3(I),Y3(I), I = 1,25)/
     .      0.13750,  0.97500,   0.91250,  0.98750,
     .      0.71250,  0.76250,   0.22500,  0.83750,
     .     -0.05000,  0.41250,   0.47500,  0.63750,
     .      0.05000, -0.05000,   0.45000,  1.03750,
     .      1.08750,  0.55000,   0.53750,  0.80000,
     .     -0.03750,  0.75000,   0.18750,  0.57500,
     .      0.71250,  0.55000,   0.85000,  0.43750,
     .      0.70000,  0.31250,   0.27500,  0.42500,
     .      0.45000,  0.28750,   0.81250,  0.18750,
     .      0.45000, -0.03750,   1.00000,  0.26250,
     .      0.50000,  0.46250,   0.18750,  0.26250,
     .      0.58750,  0.12500,   1.05000, -0.06125,
     .      0.10000,  0.11250/
C
C Node set 4:  Random 100-node set.
C
      DATA (X4(I),Y4(I), I = 1,20)/
     .      0.0096326,  0.3083158,  0.0216348,  0.2450434,
     .      0.0298360,  0.8613847,  0.0417447,  0.0977864,
     .      0.0470462,  0.3648355,  0.0562965,  0.7156339,
     .      0.0646857,  0.5311312,  0.0740377,  0.9755672,
     .      0.0873907,  0.1781117,  0.0934832,  0.5452797,
     .      0.1032216,  0.1603881,  0.1110176,  0.7837139,
     .      0.1181193,  0.9982015,  0.1251704,  0.6910589,
     .      0.1327330,  0.1049580,  0.1439536,  0.8184662,
     .      0.1564861,  0.7086405,  0.1651043,  0.4456593,
     .      0.1786039,  0.1178342,  0.1886405,  0.3189021/
      DATA (X4(I),Y4(I), I = 21,40)/
     .      0.2016706,  0.9668446,  0.2099886,  0.7571834,
     .      0.2147003,  0.2016598,  0.2204141,  0.3232444,
     .      0.2343715,  0.4368583,  0.2409660,  0.8907869,
     .      0.2527740,  0.0647260,  0.2570839,  0.5692618,
     .      0.2733365,  0.2947027,  0.2853833,  0.4332426,
     .      0.2901755,  0.3347464,  0.2964854,  0.7436284,
     .      0.3019725,  0.1066265,  0.3125695,  0.8845357,
     .      0.3307163,  0.5158730,  0.3378504,  0.9425637,
     .      0.3439061,  0.4799701,  0.3529922,  0.1783069,
     .      0.3635507,  0.1146760,  0.3766172,  0.8225797/
      DATA (X4(I),Y4(I), I = 41,60)/
     .      0.3822429,  0.2270688,  0.3869838,  0.4073598,
     .      0.3973137,  0.8875080,  0.4170708,  0.7631616,
     .      0.4255588,  0.9972804,  0.4299218,  0.4959884,
     .      0.4372839,  0.3410421,  0.4705033,  0.2498120,
     .      0.4736655,  0.6409007,  0.4879299,  0.1058690,
     .      0.4940260,  0.5411969,  0.5055324,  0.0089792,
     .      0.5162593,  0.8784268,  0.5219219,  0.5515874,
     .      0.5348529,  0.4038952,  0.5483213,  0.1654023,
     .      0.5569571,  0.2965158,  0.5638611,  0.3660356,
     .      0.5784908,  0.0366554,  0.5863950,  0.9502420/
      DATA (X4(I),Y4(I), I = 61,80)/
     .      0.5929148,  0.2638101,  0.5987839,  0.9277386,
     .      0.6117561,  0.5377694,  0.6252296,  0.7374676,
     .      0.6331381,  0.4674627,  0.6399048,  0.9186109,
     .      0.6488972,  0.0416884,  0.6558537,  0.1291029,
     .      0.6677405,  0.6763676,  0.6814074,  0.8444238,
     .      0.6887812,  0.3273328,  0.6940896,  0.1893879,
     .      0.7061687,  0.0645923,  0.7160957,  0.0180147,
     .      0.7317445,  0.8904992,  0.7370798,  0.4160648,
     .      0.7462030,  0.4688995,  0.7566957,  0.2174508,
     .      0.7699998,  0.5734231,  0.7879347,  0.8853319/
      DATA (X4(I),Y4(I), I = 81,100)/
     .      0.7944014,  0.8018436,  0.8164468,  0.6388941,
     .      0.8192794,  0.8931002,  0.8368405,  0.1000558,
     .      0.8500993,  0.2789506,  0.8588255,  0.9082948,
     .      0.8646496,  0.3259159,  0.8792329,  0.8318747,
     .      0.8837536,  0.0508513,  0.8900077,  0.9708450,
     .      0.8969894,  0.5120548,  0.9044917,  0.2859716,
     .      0.9083947,  0.9581641,  0.9203972,  0.6183429,
     .      0.9347906,  0.3779934,  0.9434519,  0.4010423,
     .      0.9490328,  0.9478657,  0.9569571,  0.7425486,
     .      0.9772067,  0.8883287,  0.9983493,  0.5496750/
C
C Node set 5:  9 by 9 uniform grid.
C
      DATA (X5(I),Y5(I), I = 1,20)/
     .      0.125,  0.000,  0.000,  0.125,
     .      0.000,  0.250,  0.000,  0.375,
     .      0.000,  0.500,  0.000,  0.625,
     .      0.000,  0.750,  0.000,  0.875,
     .      0.000,  1.000,  0.000,  0.000,
     .      0.125,  0.125,  0.125,  0.250,
     .      0.125,  0.375,  0.125,  0.500,
     .      0.125,  0.625,  0.125,  0.750,
     .      0.125,  0.875,  0.125,  1.000,
     .      0.250,  0.000,  0.250,  0.125/
      DATA (X5(I),Y5(I), I = 21,40)/
     .      0.250,  0.250,  0.250,  0.375,
     .      0.250,  0.500,  0.250,  0.625,
     .      0.250,  0.750,  0.250,  0.875,
     .      0.250,  1.000,  0.375,  0.000,
     .      0.375,  0.125,  0.375,  0.250,
     .      0.375,  0.375,  0.375,  0.500,
     .      0.375,  0.625,  0.375,  0.750,
     .      0.375,  0.875,  0.375,  1.000,
     .      0.500,  0.000,  0.500,  0.125,
     .      0.500,  0.250,  0.500,  0.375/
      DATA (X5(I),Y5(I), I = 41,60)/
     .      0.500,  0.500,  0.500,  0.625,
     .      0.500,  0.750,  0.500,  0.875,
     .      0.500,  1.000,  0.625,  0.000,
     .      0.625,  0.125,  0.625,  0.250,
     .      0.625,  0.375,  0.625,  0.500,
     .      0.625,  0.625,  0.625,  0.750,
     .      0.625,  0.875,  0.625,  1.000,
     .      0.750,  0.000,  0.750,  0.125,
     .      0.750,  0.250,  0.750,  0.375,
     .      0.750,  0.500,  0.750,  0.625/
      DATA (X5(I),Y5(I), I = 61,81)/
     .      0.750,  0.750,  0.750,  0.875,
     .      0.750,  1.000,  0.875,  0.000,
     .      0.875,  0.125,  0.875,  0.250,
     .      0.875,  0.375,  0.875,  0.500,
     .      0.875,  0.625,  0.875,  0.750,
     .      0.875,  0.875,  0.875,  1.000,
     .      1.000,  0.000,  1.000,  0.125,
     .      1.000,  0.250,  1.000,  0.375,
     .      1.000,  0.500,  1.000,  0.625,
     .      1.000,  0.750,  1.000,  0.875,
     .      1.000,  1.000/
C
C Store node set K in (X,Y).
C
      IF (K .EQ. 1) THEN
        DO 1 I = 1,100
          X(I) = X1(I)
          Y(I) = Y1(I)
    1     CONTINUE
        N = 100
      ELSEIF (K .EQ. 2) THEN
        DO 2 I = 1,33
          X(I) = X2(I)
          Y(I) = Y2(I)
    2     CONTINUE
        N = 33
      ELSEIF (K .EQ. 3) THEN
        DO 3 I = 1,25
          X(I) = X3(I)
          Y(I) = Y3(I)
    3     CONTINUE
        N = 25
      ELSEIF (K .EQ. 4) THEN
        DO 4 I = 1,100
          X(I) = X4(I)
          Y(I) = Y4(I)
    4     CONTINUE
        N = 100
      ELSEIF (K .EQ. 5) THEN
        DO 5 I = 1,81
          X(I) = X5(I)
          Y(I) = Y5(I)
    5     CONTINUE
        N = 81
      ELSE
        N = 0
      ENDIF
      RETURN
      END
      
      SUBROUTINE TSTFN1_sg (K,X,Y,IFLAG, F,FX,FY)
      INTEGER K, IFLAG
      real X, Y, F, FX, FY
C
C***********************************************************
C
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/14/98
C
C   This subroutine computes the value and, optionally, the
C first partial derivatives of one of ten bivariate test
C functions.  The first six functions were chosen by Richard
C Franke to test interpolation software (See the reference
C below).  The last four functions represent more chal-
C lenging surface fitting problems.
C
C On input:
C
C       K = Integer in the range 1 to 10 which determines
C           the choice of function as follows:
C
C               K = 1 - Exponential
C               K = 2 - Cliff
C               K = 3 - Saddle
C               K = 4 - Gentle
C               K = 5 - Steep
C               K = 6 - Sphere
C               K = 7 - Trig
C               K = 8 - Synergistic Gaussian
C               K = 9 - Cloverleaf Asymmetric Peak/Valley
C               K = 10 - Cosine Peak
C
C   Note that function 6 is only defined inside a circle of
C radius 8/9 centered at (.5,.5).  Thus, if (X-.5)**2 +
C (Y-.5)**2 .GE. 64/81, the value (and partials if IFLAG=1)
C are set to 0 for this function.  Also, the first partial
C derivatives of function 10 are not defined at (.5,.5) --
C again, zeros are returned.
C
C       X,Y = Coordinates of the point at which the selected
C             function is to be evaluated.
C
C       IFLAG = Derivative option indicator:
C               IFLAG = 0 if only a function value is
C                         required.
C               IFLAG = 1 if both the function and its first
C                         partial derivatives are to be
C                         evaluated.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       F = Value of function K at (X,Y).
C
C       FX,FY = First partial derivatives of function K at
C               (X,Y) if IFLAG = 1, unaltered otherwise.
C
C Intrinsic functions called by TSTFN1:  COS, EXP, SIN,
C                                          SQRT, TANH
C
C Reference:  R. Franke, A Critical Comparison of Some
C               Methods for Interpolation of Scattered Data,
C               Naval Postgraduate School Technical Report,
C               NPS-53-79-003, 1979.
C
C***********************************************************
C
      real T1, T2, T3, T4
      IF (K .LT. 1  .OR.  K .GT. 10) RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10), K
C
C Exponential:
C
    1 F = .75*EXP(-((9.*X-2.)**2 + (9.*Y-2.)**2)/4.) +
     .    .75*EXP(-((9.*X+1.)**2)/49. - (9.*Y+1.)/10.) +
     .     .5*EXP(-((9.*X-7.)**2 + (9.*Y-3.)**2)/4.) -
     .     .2*EXP(-(9.*X-4.)**2 - (9.*Y-7.)**2)
      IF (IFLAG .NE. 1) RETURN
      T1 = EXP(-((9.*X-2.)**2 + (9.*Y-2.)**2)/4.)
      T2 = EXP(-((9.*X+1.)**2)/49. - (9.*Y+1.)/10.)
      T3 = EXP(-((9.*X-7.)**2 + (9.*Y-3.)**2)/4.)
      T4 = EXP(-(9.*X-4.)**2 - (9.*Y-7.)**2)
      FX = -3.375*(9.*X-2.)*T1 - (27./98.)*(9.*X+1.)*T2
     .     -2.25*(9.*X-7.)*T3 + 3.6*(9.*X-4.)*T4
      FY = -3.375*(9.*Y-2.)*T1 - .675*T2
     .     -2.25*(9.*Y-3.)*T3 + 3.6*(9.*Y-7.)*T4
      RETURN
C
C Cliff:
C
    2 F = (TANH(9.0*(Y-X)) + 1.0)/9.0
      IF (IFLAG .NE. 1) RETURN
      T1 = 18.0*(Y-X)
      FX = -4.0/(EXP(T1) + 2.0 + EXP(-T1))
      FY = -FX
      RETURN
C
C Saddle:
C
    3 F = (1.25 + COS(5.4*Y))/(6.0 + 6.0*(3.0*X-1.0)**2)
      IF (IFLAG .NE. 1) RETURN
      T1 = 5.4*Y
      T2 = 1.0 + (3.0*X-1.)**2
      FX = -(3.0*X-1.0)*(1.25 + COS(T1))/(T2**2)
      FY = -.9*SIN(T1)/T2
      RETURN
C
C Gentle:
C
    4 F = EXP(-5.0625*((X-.5)**2 + (Y-.5)**2))/3.0
      IF (IFLAG .NE. 1) RETURN
      T1 = X - .5
      T2 = Y - .5
      T3 = -3.375*EXP(-5.0625*(T1**2 + T2**2))
      FX = T1*T3
      FY = T2*T3
      RETURN
C
C Steep:
C
    5 F = EXP(-20.25*((X-.5)**2 + (Y-.5)**2))/3.0
      IF (IFLAG .NE. 1) RETURN
      T1 = X - .5
      T2 = Y - .5
      T3 = -13.5*EXP(-20.25*(T1**2 + T2**2))
      FX = T1*T3
      FY = T2*T3
      RETURN
C
C Sphere:
C
    6 T4 = 64.0 - 81.0*((X-.5)**2 + (Y-.5)**2)
      F = 0.
      IF (T4 .GE. 0.) F = SQRT(T4)/9.0 - .5
      IF (IFLAG .NE. 1) RETURN
      T1 = X - .5
      T2 = Y - .5
      T3 = 0.
      IF (T4 .GT. 0.) T3 = -9.0/SQRT(T4)
      FX = T1*T3
      FY = T2*T3
      RETURN
C
C Trig:
C
    7 F = 2.0*COS(10.0*X)*SIN(10.0*Y) + SIN(10.0*X*Y)
      IF (IFLAG .NE. 1) RETURN
      T1 = 10.0*X
      T2 = 10.0*Y
      T3 = 10.0*COS(10.0*X*Y)
      FX = -20.0*SIN(T1)*SIN(T2) + T3*Y
      FY = 20.0*COS(T1)*COS(T2) + T3*X
      RETURN
C
C Gaussx(1,.5,.1) + Gaussy(.75,.5,.1) + Gaussx(1,.5,.1)*
C   Gaussy(.75,.5,.1), where Gaussx(a,b,c) is the Gaussian
C   function of x with amplitude a, center (mean) b, and
C   width (standard deviation) c.
C
    8 T1 = 5.0 - 10.0*X
      T2 = 5.0 - 10.0*Y
      T3 = EXP(-.5*T1*T1)
      T4 = EXP(-.5*T2*T2)
      F = T3 + .75*T4*(1.0+T3)
      IF (IFLAG .NE. 1) RETURN
      FX = T1*T3*(10.0 + 7.5*T4)
      FY = T2*T4*(7.5 + 7.5*T3)
      RETURN
C
C Cloverleaf Asymmetric Hill/Valley:
C
    9 T1 = EXP((10.0 - 20.0*X)/3.0)
      T2 = EXP((10.0 - 20.0*Y)/3.0)
      T3 = 1.0/(1.0 + T1)
      T4 = 1.0/(1.0 + T2)
      F = ((20.0/3.0)**3 * T1*T2)**2 * (T3*T4)**5 *
     .    (T1-2.0*T3)*(T2-2.0*T4)
      IF (IFLAG .NE. 1) RETURN
      FX = ((20.0/3.0)*T1)**2 * ((20.0/3.0)*T3)**5 *
     .     (2.0*T1-3.0*T3-5.0+12.0*T3*T3)*T2*T2*T4**5 *
     .     (T2-2.0*T4)
      FY = ((20.0/3.0)*T1)**2 * ((20.0/3.0)*T3)**5 *
     .     (2.0*T2-3.0*T4-5.0+12.0*T4*T4)*T2*T2*T4**5 *
     .     (T1-2.0*T3)
      RETURN
C
C Cosine Peak:
C
   10 T1 = SQRT( (80.0*X - 40.0)**2 + (90.0*Y - 45.0)**2 )
      T2 = EXP(-.04*T1)
      T3 = COS(.15*T1)
      F = T2*T3
      IF (IFLAG .NE. 1) RETURN
      T4 = SIN(.15*T1)
      FX = 0.
      FY = 0.
      IF (T1 .EQ. 0.) RETURN
      T4 = SIN(.15*T1)
      FX = -T2*(12.0*T4 + 3.2*T3)*(80.0*X - 40.0)/T1
      FY = -T2*(13.5*T4 + 3.6*T3)*(90.0*Y - 45.0)/T1
      RETURN
      END
      
      
      
      SUBROUTINE TSTFN2_sg (K,X,Y,IFLAG, F,FX,FY,FXX,FXY,FYY)
      INTEGER K, IFLAG
      real X, Y, F, FX, FY, FXX, FXY, FYY
C
C***********************************************************
C
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/14/98
C
C   This subroutine computes the value and, optionally, the
C first and/or second partial derivatives of one of ten 
C bivariate test functions.  The first six functions were
C chosen by Richard Franke to test interpolation software.
C (See the reference below.)  The last four functions repre-
C sent more challenging surface fitting problems.
C
C On input:
C
C       K = Integer in the range 1 to 10 which determines
C           the choice of function as follows:
C
C               K = 1 - Exponential
C               K = 2 - Cliff
C               K = 3 - Saddle
C               K = 4 - Gentle
C               K = 5 - Steep
C               K = 6 - Sphere
C               K = 7 - Trig
C               K = 8 - Synergistic Gaussian
C               K = 9 - Cloverleaf Asymmetric Peak/Valley
C               K = 10 - Cosine Peak
C
C   Note that function 6 is only defined inside a circle of
C radius 8/9 centered at (.5,.5).  Thus, if (X-.5)**2 +
C (Y-.5)**2 .GE. 64/81, the value (and partials if IFLAG=1)
C are set to 0 for this function.  Also, the first partial
C derivatives of function 10 are not defined at (.5,.5) --
C again, zeros are returned.
C
C       X,Y = Coordinates of the point at which the selected
C             function is to be evaluated.
C
C       IFLAG = Derivative option indicator:
C               IFLAG = 0 if only a function value is
C                         required.
C               IFLAG = 1 if both the function and its first
C                         partial derivatives are to be
C                         evaluated.
C               IFLAG = 2 if the function, its first partial
C                         derivatives, and its second partial
C                         derivatives are to be evaluated.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       F = Value of function K at (X,Y).
C
C       FX,FY = First partial derivatives of function K at
C               (X,Y) if IFLAG >= 1, unaltered otherwise.
C
C       FXX,FXY,FYY = Second partial derivatives of function
C                     K at (X,Y) if IFLAG >= 2, unaltered
C                     otherwise.
C
C Intrinsic functions called by TSTFN2:  COS, EXP, SIN,
C                                          SQRT, TANH
C
C Reference:  R. Franke, A Critical Comparison of Some
C               Methods for Interpolation of Scattered Data,
C               Naval Postgraduate School Technical Report,
C               NPS-53-79-003, 1979.
C
C***********************************************************
C
      real T1, T2, T3, T4, T5, T6
      IF (K .LT. 1  .OR.  K .GT. 10) RETURN
      GO TO (1,2,3,4,5,6,7,8,9,10), K
C
C Exponential:
C
    1 F = .75*EXP(-((9.*X-2.)**2 + (9.*Y-2.)**2)/4.) +
     .    .75*EXP(-((9.*X+1.)**2)/49. - (9.*Y+1.)/10.) +
     .     .5*EXP(-((9.*X-7.)**2 + (9.*Y-3.)**2)/4.) -
     .     .2*EXP(-(9.*X-4.)**2 - (9.*Y-7.)**2)
      IF (IFLAG .LT. 1) RETURN
      T1 = EXP(-((9.*X-2.)**2 + (9.*Y-2.)**2)/4.)
      T2 = EXP(-((9.*X+1.)**2)/49. - (9.*Y+1.)/10.)
      T3 = EXP(-((9.*X-7.)**2 + (9.*Y-3.)**2)/4.)
      T4 = EXP(-(9.*X-4.)**2 - (9.*Y-7.)**2)
      FX = -3.375*(9.*X-2.)*T1 - (27./98.)*(9.*X+1.)*T2
     .     -2.25*(9.*X-7.)*T3 + 3.6*(9.*X-4.)*T4
      FY = -3.375*(9.*Y-2.)*T1 - .675*T2
     .     -2.25*(9.*Y-3.)*T3 + 3.6*(9.*Y-7.)*T4
      IF (IFLAG .LT. 2) RETURN
      FXX = 15.1875*((9.*X-2.)**2 - 2.)*T1 +
     .      60.75*((9.*X+1.)**2 - 24.5)*T2 +
     .      10.125*((9.*X-7.)**2 - 2.)*T3 -
     .      64.8*((9.*X-4.)**2 - .5)*T4
      FXY = 15.1875*(9.*X-2.)*(9.*Y-2.)*T1 +
     .      (243./980.)*(9.*X+1.)*T2 +
     .      10.125*(9.*X-7.)*(9.*Y-3.)*T3 -
     .      64.8*(9.*X-4.)*(9.*Y-7.)*T4
      FYY = 15.1875*((9.*Y-2.)**2 - 2.)*T1 +
     .      .6075*T2 +
     .      10.125*((9.*Y-3.)**2 - 2.)*T3 -
     .      64.8*((9.*Y-7.)**2 - .5)*T4
      RETURN
C
C Cliff:
C
    2 F = (TANH(9.0*(Y-X)) + 1.0)/9.0
      IF (IFLAG .LT. 1) RETURN
      T1 = 18.0*(Y-X)
      FX = -4.0/(EXP(T1) + 2.0 + EXP(-T1))
      FY = -FX
      IF (IFLAG .LT. 2) RETURN
      FXX = 18.0*TANH(0.5*T1)*FX
      FXY = -FXX
      FYY = FXX
      RETURN
C
C Saddle:
C
    3 F = (1.25 + COS(5.4*Y))/(6.0 + 6.0*(3.0*X-1.0)**2)
      IF (IFLAG .LT. 1) RETURN
      T1 = 5.4*Y
      T2 = 1.0 + (3.0*X-1.)**2
      FX = -(3.0*X-1.0)*(1.25 + COS(T1))/(T2**2)
      FY = -.9*SIN(T1)/T2
      IF (IFLAG .LT. 2) RETURN
      FXX = 3.0*(1.25 + COS(T1))*(3.0*T2-4.0)/(T2**3)
      FXY = 5.4*(3.0*X-1.0)*SIN(T1)/(T2**2)
      FYY = -4.86*COS(T1)/T2
      RETURN
C
C Gentle:
C
    4 F = EXP(-5.0625*((X-.5)**2 + (Y-.5)**2))/3.0
      IF (IFLAG .LT. 1) RETURN
      T1 = X - .5
      T2 = Y - .5
      T3 = -3.375*EXP(-5.0625*(T1**2 + T2**2))
      FX = T1*T3
      FY = T2*T3
      IF (IFLAG .LT. 2) RETURN
      FXX = (1.0 - 10.125*T1*T1)*T3
      FXY = -10.125*T1*T2*T3
      FYY = (1.0 - 10.125*T2*T2)*T3
      RETURN
C
C Steep:
C
    5 F = EXP(-20.25*((X-.5)**2 + (Y-.5)**2))/3.0
      IF (IFLAG .LT. 1) RETURN
      T1 = X - .5
      T2 = Y - .5
      T3 = -13.5*EXP(-20.25*(T1**2 + T2**2))
      FX = T1*T3
      FY = T2*T3
      IF (IFLAG .LT. 2) RETURN
      FXX = (1.0 - 40.5*T1*T1)*T3
      FXY = -40.5*T1*T2*T3
      FYY = (1.0 - 40.5*T2*T2)*T3
      RETURN
C
C Sphere:
C
    6 T4 = 64.0 - 81.0*((X-.5)**2 + (Y-.5)**2)
      F = 0.
      IF (T4 .GE. 0.) F = SQRT(T4)/9.0 - .5
      IF (IFLAG .LT. 1) RETURN
      T1 = X - .5
      T2 = Y - .5
      T3 = 0.
      IF (T4 .GT. 0.) T3 = -9.0/SQRT(T4)
      FX = T1*T3
      FY = T2*T3
      IF (IFLAG .LT. 2) RETURN
      FXX = (1.0 + FX*FX)*T3
      FXY = FX*FY
      FYY = (1.0 + FY*FY)*T3
      RETURN
C
C Trig:
C
    7 F = 2.0*COS(10.0*X)*SIN(10.0*Y) + SIN(10.0*X*Y)
      IF (IFLAG .LT. 1) RETURN
      T1 = 10.0*X
      T2 = 10.0*Y
      T3 = 10.0*COS(10.0*X*Y)
      FX = -20.0*SIN(T1)*SIN(T2) + T3*Y
      FY = 20.0*COS(T1)*COS(T2) + T3*X
      IF (IFLAG .LT. 2) RETURN
      T4 = 100.0*SIN(10.0*X*Y)
      FXX = -200.0*COS(T1)*SIN(T2) - T4*Y*Y
      FXY = -200.0*SIN(T1)*COS(T2) + T3 - T4*X*Y
      FYY = -200.0*COS(T1)*SIN(T2) - T4*X*X
      RETURN
C
C Gaussx(1,.5,.1) + Gaussy(.75,.5,.1) + Gaussx(1,.5,.1)*
C   Gaussy(.75,.5,.1), where Gaussx(a,b,c) is the Gaussian
C   function of x with amplitude a, center (mean) b, and
C   width (standard deviation) c.
C
    8 T1 = 5.0 - 10.0*X
      T2 = 5.0 - 10.0*Y
      T3 = EXP(-.5*T1*T1)
      T4 = EXP(-.5*T2*T2)
      F = T3 + .75*T4*(1.0+T3)
      IF (IFLAG .LT. 1) RETURN
      FX = T1*T3*(10.0 + 7.5*T4)
      FY = T2*T4*(7.5 + 7.5*T3)
      IF (IFLAG .LT. 2) RETURN
      FXX = T3*(T1*T1-1.0)*(100.0 + 75.0*T4)
      FXY = 75.0*T1*T2*T3*T4
      FYY = T4*(T2*T2-1.0)*(75.0 + 75.0*T3)
      RETURN
C
C Cloverleaf Asymmetric Hill/Valley:
C
    9 T1 = EXP((10.0 - 20.0*X)/3.0)
      T2 = EXP((10.0 - 20.0*Y)/3.0)
      T3 = 1.0/(1.0 + T1)
      T4 = 1.0/(1.0 + T2)
      T5 = 20.0/3.0
      F = (T5**3 * T1*T2)**2 * (T3*T4)**5 *
     .    (T1-2.0*T3)*(T2-2.0*T4)
      IF (IFLAG .LT. 1) RETURN
      T6 = (T5*T1*T2)**2 * (T5*T3*T4)**5
      FX = T6 * (T2-2.0*T4) *
     .     ((12.0*T3 - 3.0)*T3 + 2.0*T1 - 5.0)
      FY = T6 * (T1-2.0*T3) *
     .     ((12.0*T4 - 3.0)*T4 + 2.0*T2 - 5.0)
      IF (IFLAG .LT. 2) RETURN
      FXX = T5*T6 * (T2-2.0*T4) *
     .      (((-84.0*T3 + 78.0)*T3 + 23.0)*T3 + 4.0*T1-25.0)
      FXY = T5*T6 *
     .      ((12.0*T4 - 3.0)*T4 + 2.0*T2 - 5.0) *
     .      ((12.0*T3 - 3.0)*T3 + 2.0*T1 - 5.0)
      FYY = T5*T6 * (T1-2.0*T3) *
     .      (((-84.0*T4 + 78.0)*T4 + 23.0)*T4 + 4.0*T2-25.0)
      RETURN
C
C Cosine Peak:
C
   10 T1 = SQRT( (80.0*X - 40.0)**2 + (90.0*Y - 45.0)**2 )
      T2 = EXP(-.04*T1)
      T3 = COS(.15*T1)
      F = T2*T3
      IF (IFLAG .LT. 1) RETURN
      T4 = SIN(.15*T1)
      FX = 0.
      FY = 0.
      IF (T1 .EQ. 0.) RETURN
      T4 = SIN(.15*T1)
      FX = -T2*(12.0*T4 + 3.2*T3)*(80.0*X - 40.0)/T1
      FY = -T2*(13.5*T4 + 3.6*T3)*(90.0*Y - 45.0)/T1
      IF (IFLAG .LT. 2) RETURN
      FXX = 0.
      FXY = 0.
      FYY = 0.
      IF (T1 .EQ. 0.) RETURN
      T5 = T2/(T1**3)
      FXX = T5*(T1*(76.8*T4 - 133.76*T3)*(80.0*X-40.0)**2 -
     .             (960.0*T4 + 256.0*T3)*(90.0*Y-45.0)**2 )
      FXY = T5*(T1*(86.4*T4 - 150.48*T3) + 1080.0*T4 +
     .          288.0*T3)*(80.0*X-40.0)*(90.0*Y-45.0)
      FYY = T5*(T1*(97.2*T4 - 169.29*T3)*(90.0*Y-45.0)**2 -
     .             (1215.0*T4 + 324.0*T3)*(80.0*X-40.0)**2)
      RETURN
      END

