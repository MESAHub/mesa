      SUBROUTINE F01BRT(NN,ICN,A,LICN,LENR,LENRL,IDISP,IP,IQ,IRN,LIRN,
     *                  LENC,IFIRST,LASTR,NEXTR,LASTC,NEXTC,IPTR,IPC,U,
     *                  ABORT,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DERIVED FROM HARWELL LIBRARY ROUTINE MA30A
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  U
      INTEGER           IFAIL, LICN, LIRN, NN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LICN)
      INTEGER           ICN(LICN), IDISP(7), IFIRST(NN), IP(NN),
     *                  IPC(NN), IPTR(NN), IQ(NN), IRN(LIRN), LASTC(NN),
     *                  LASTR(NN), LENC(NN), LENR(NN), LENRL(NN),
     *                  NEXTC(NN), NEXTR(NN)
      LOGICAL           ABORT(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMAX, AU, UMAX, ZERO
      INTEGER           DISPC, I, I1, I2, IACTIV, IBEG, ICNCP, IDISPC,
     *                  IDUMMY, IEND, IFILL, IFIR, II, III, IJFIR, IJP1,
     *                  IJPOS, ILAST, INDROW, IOP, IPIV, IPOS, IRANK,
     *                  IRNCP, IROWS, ISAVE, ISING, ISTART, ISW, ISW1,
     *                  ITOP, J, J1, J2, JBEG, JCOST, JCOUNT, JDIFF,
     *                  JEND, JJ, JMORE, JNPOS, JOLD, JPIV, JPOS, JROOM,
     *                  JVAL, JZER, JZERO, K, KCOST, L, LC, LENPIV, LL,
     *                  LR, MINICN, MINIRN, MOREI, N, NBLOCK, NC, NERR,
     *                  NNM1, NR, NUM, NZ, NZ2, NZCOL, NZMIN, NZPC,
     *                  NZROW, OLDEND, OLDPIV, PIVEND, PIVOT, PIVROW,
     *                  ROWI
C     .. Local Arrays ..
      CHARACTER*65      REC(3)
C     .. External Subroutines ..
      EXTERNAL          F01BRS, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, MOD
C     .. Data statements ..
      DATA              UMAX/.9999D0/
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 0
C     NERR IS THE UNIT NUMBER FOR ERROR MESSAGES
      CALL X04AAF(0,NERR)
      MINIRN = 0
      MINICN = IDISP(1) - 1
      MOREI = 0
      IRANK = NN
      IRNCP = 0
      ICNCP = 0
C     RESET U IF NECESSARY.
      U = MIN(U,UMAX)
C     IBEG IS THE POSITION OF THE NEXT PIVOT ROW AFTER ELIMINATION
C     STEP USING IT.
      U = MAX(U,ZERO)
      IBEG = IDISP(1)
C     IACTIV IS THE POSITION OF THE FIRST ENTRY IN THE ACTIVE PART
C     OF A/ICN.
      IACTIV = IDISP(2)
C     NZROW IS CURRENT NUMBER OF NON-ZEROS IN ACTIVE AND UNPROCESSED
C     PART OF ROW FILE ICN.
      NZROW = LICN - IACTIV + 1
      MINICN = NZROW + MINICN
C
C     COUNT THE NUMBER OF DIAGONAL BLOCKS AND SET UP POINTERS TO
C     THE BEGINNINGS OF THE ROWS.
C     NUM IS THE NUMBER OF DIAGONAL BLOCKS.
      NUM = 1
      IPTR(1) = IACTIV
      IF (NN.EQ.1) GO TO 40
      NNM1 = NN - 1
      DO 20 I = 1, NNM1
         IF (IP(I).LT.0) NUM = NUM + 1
         IPTR(I+1) = IPTR(I) + LENR(I)
   20 CONTINUE
C     ILAST IS THE LAST ROW IN THE PREVIOUS BLOCK.
   40 ILAST = 0
C
C     ***********************************************
C     ****    LU DECOMPOSITION OF BLOCK NBLOCK   ****
C     ***********************************************
C
C     EACH PASS THROUGH THIS LOOP PERFORMS LU DECOMPOSITION ON ONE
C     OF THE DIAGONAL BLOCKS.
      DO 1820 NBLOCK = 1, NUM
         ISTART = ILAST + 1
         DO 60 IROWS = ISTART, NN
            IF (IP(IROWS).LT.0) GO TO 80
   60    CONTINUE
         IROWS = NN
   80    ILAST = IROWS
C        N IS THE NUMBER OF ROWS IN THE CURRENT BLOCK.
C        ISTART IS THE INDEX OF THE FIRST ROW IN THE CURRENT BLOCK.
C        ILAST IS THE INDEX OF THE LAST ROW IN THE CURRENT BLOCK.
C        IACTIV IS THE POSITION OF THE FIRST ELEMENT IN THE BLOCK.
C        ITOP IS THE POSITION OF THE LAST ELEMENT IN THE BLOCK.
         N = ILAST - ISTART + 1
         IF (N.NE.1) GO TO 160
C
C        CODE FOR DEALING WITH 1X1 BLOCK.
         LENRL(ILAST) = 0
         ISING = ISTART
         IF (LENR(ILAST).NE.0) GO TO 100
C        BLOCK IS STRUCTURALLY SINGULAR.
         IRANK = IRANK - 1
         ISING = -ISING
         IF (IFAIL.NE.-2 .AND. IFAIL.NE.5) IFAIL = -1
         IF ( .NOT. ABORT(1)) GO TO 140
         IFAIL = 1
         GO TO 2020
  100    IF (A(IACTIV).NE.ZERO) GO TO 120
C        BLOCK IS NUMERICALLY SINGULAR
         ISING = -ISING
         IRANK = IRANK - 1
         IPTR(ILAST) = 0
         IF (IFAIL.NE.5) IFAIL = -2
         IF ( .NOT. ABORT(2)) GO TO 120
         IFAIL = 2
         GO TO 2020
  120    A(IBEG) = A(IACTIV)
         ICN(IBEG) = ICN(IACTIV)
         IACTIV = IACTIV + 1
         IPTR(ISTART) = 0
         IBEG = IBEG + 1
         NZROW = NZROW - 1
  140    LASTR(ISTART) = ISTART
         LASTC(ISTART) = ISING
         GO TO 1820
C
C        NON-TRIVIAL BLOCK.
  160    ITOP = LICN
         IF (ILAST.NE.NN) ITOP = IPTR(ILAST+1) - 1
C
C        SET UP COLUMN ORIENTED STORAGE.
         DO 180 I = ISTART, ILAST
            LENRL(I) = 0
            LENC(I) = 0
  180    CONTINUE
         IF (ITOP-IACTIV.LT.LIRN) GO TO 200
         MINIRN = ITOP - IACTIV + 1
         PIVOT = ISTART - 1
         GO TO 1980
C
C        CALCULATE COLUMN COUNTS.
  200    DO 220 II = IACTIV, ITOP
            I = ICN(II)
            LENC(I) = LENC(I) + 1
  220    CONTINUE
C        SET UP COLUMN POINTERS SO THAT IPC(J) POINTS TO POSITION
C        AFTER END OF COLUMN J IN COLUMN FILE.
         IPC(ILAST) = LIRN + 1
         J1 = ISTART + 1
         DO 240 JJ = J1, ILAST
            J = ILAST - JJ + J1 - 1
            IPC(J) = IPC(J+1) - LENC(J+1)
  240    CONTINUE
         DO 280 INDROW = ISTART, ILAST
            J1 = IPTR(INDROW)
            J2 = J1 + LENR(INDROW) - 1
            IF (J1.GT.J2) GO TO 280
            DO 260 JJ = J1, J2
               J = ICN(JJ)
               IPOS = IPC(J) - 1
               IRN(IPOS) = INDROW
               IPC(J) = IPOS
  260       CONTINUE
  280    CONTINUE
C        DISPC IS THE LOWEST INDEXED ACTIVE LOCATION IN THE COLUMN
C        FILE.
         DISPC = IPC(ISTART)
         NZCOL = LIRN - DISPC + 1
         MINIRN = MAX(NZCOL,MINIRN)
         NZMIN = 1
C
C        INITIALIZE ARRAY IFIRST.  IFIRST(I) = +/- K INDICATES THAT
C        ROW/COL K HAS I NON-ZEROS.  IF IFIRST(I) = 0, THERE IS NO ROW
C        OR COLUMN WITH I NON ZEROS.
         DO 300 I = 1, N
            IFIRST(I) = 0
  300    CONTINUE
C
C        COMPUTE ORDERING OF ROW AND COLUMN COUNTS.
C        FIRST RUN THROUGH COLUMNS (FROM COLUMN N TO COLUMN 1).
         DO 340 JJ = ISTART, ILAST
            J = ILAST - JJ + ISTART
            NZ = LENC(J)
            IF (NZ.NE.0) GO TO 320
            IPC(J) = 0
            LASTC(J) = 0
            GO TO 340
  320       ISW = IFIRST(NZ)
            IFIRST(NZ) = -J
            LASTC(J) = 0
            NEXTC(J) = -ISW
            ISW1 = ABS(ISW)
            IF (ISW.NE.0) LASTC(ISW1) = J
  340    CONTINUE
C        NOW RUN THROUGH ROWS (AGAIN FROM N TO 1).
         DO 400 II = ISTART, ILAST
            I = ILAST - II + ISTART
            NZ = LENR(I)
            IF (NZ.NE.0) GO TO 360
            IPTR(I) = 0
            LASTR(I) = 0
            GO TO 400
  360       ISW = IFIRST(NZ)
            IFIRST(NZ) = I
            IF (ISW.GT.0) GO TO 380
            NEXTR(I) = 0
            LASTR(I) = ISW
            GO TO 400
  380       NEXTR(I) = ISW
            LASTR(I) = LASTR(ISW)
            LASTR(ISW) = I
  400    CONTINUE
C
C        **********************************************
C        ****    START OF MAIN ELIMINATION LOOP    ****
C        **********************************************
C
         DO 1780 PIVOT = ISTART, ILAST
C
C           FIRST FIND THE PIVOT USING MARKOWITZ CRITERION WITH
C           STABILITY CONTROL.
C           JCOST IS THE MARKOWITZ COST OF THE BEST PIVOT SO FAR,.. THIS
C           PIVOT IS IN ROW IPIV AND COLUMN JPIV.
            NZ2 = NZMIN
            JCOST = N*N
C
C           EXAMINE ROWS/COLUMNS IN ORDER OF ASCENDING COUNT.
            DO 660 L = 1, 2
               LL = L
C              A PASS WITH L EQUAL TO 2 IS ONLY PERFORMED IN THE CASE OF
C              SINGULARITY.
               DO 640 NZ = NZ2, N
                  IF (JCOST.LE.(NZ-1)**2) GO TO 820
                  IJFIR = IFIRST(NZ)
                  IF (IJFIR) 440, 420, 460
  420             IF (LL.EQ.1) NZMIN = NZ + 1
                  GO TO 640
  440             LL = 2
                  IJFIR = -IJFIR
                  GO TO 560
  460             LL = 2
C                 SCAN ROWS WITH NZ NON-ZEROS.
                  DO 520 IDUMMY = 1, N
                     IF (IJFIR.EQ.0) GO TO 540
C                    ROW IJFIR IS NOW EXAMINED.
                     I = IJFIR
                     IJFIR = NEXTR(I)
C                    FIRST CALCULATE MULTIPLIER THRESHOLD LEVEL.
                     AMAX = ZERO
                     J1 = IPTR(I) + LENRL(I)
                     J2 = IPTR(I) + LENR(I) - 1
                     DO 480 JJ = J1, J2
                        AMAX = MAX(AMAX,ABS(A(JJ)))
  480                CONTINUE
                     AU = AMAX*U
C                    SCAN ROW FOR POSSIBLE PIVOTS
                     DO 500 JJ = J1, J2
                        IF (ABS(A(JJ)).LE.AU .AND. L.EQ.1) GO TO 500
                        J = ICN(JJ)
                        KCOST = (NZ-1)*(LENC(J)-1)
                        IF (KCOST.GE.JCOST) GO TO 500
C                       BEST PIVOT SO FAR IS FOUND.
                        JCOST = KCOST
                        IJPOS = JJ
                        IPIV = I
                        JPIV = J
                        IF (JCOST.LE.(NZ-1)**2) GO TO 820
  500                CONTINUE
  520             CONTINUE
C
C                 COLUMNS WITH NZ NON-ZEROS NOW EXAMINED.
  540             IJFIR = IFIRST(NZ)
                  IJFIR = -LASTR(IJFIR)
  560             IF (JCOST.LE.NZ*(NZ-1)) GO TO 820
                  DO 620 IDUMMY = 1, N
                     IF (IJFIR.EQ.0) GO TO 640
                     J = IJFIR
                     IJFIR = NEXTC(IJFIR)
                     I1 = IPC(J)
                     I2 = I1 + NZ - 1
C                    SCAN COLUMN J.
                     DO 600 II = I1, I2
                        I = IRN(II)
                        KCOST = (NZ-1)*(LENR(I)-LENRL(I)-1)
                        IF (KCOST.GE.JCOST) GO TO 600
C                       PIVOT HAS BEST MARKOWITZ COUNT SO FAR ...
C                       NOW CHECK ITS SUITABILITY ON NUMERIC GROUNDS
C                       BY EXAMINING THE OTHER NON-ZEROS IN ITS ROW.
                        J1 = IPTR(I) + LENRL(I)
                        J2 = IPTR(I) + LENR(I) - 1
C                       WE NEED A STABILITY CHECK ON SINGLETON COLUMNS
C                       BECAUSE OF POSSIBLE PROBLEMS WITH
C                       UNDERDETERMINED SYSTEMS.
                        AMAX = ZERO
                        DO 580 JJ = J1, J2
                           AMAX = MAX(AMAX,ABS(A(JJ)))
                           IF (ICN(JJ).EQ.J) JPOS = JJ
  580                   CONTINUE
                        IF (ABS(A(JPOS)).LE.AMAX*U .AND. L.EQ.1)
     *                      GO TO 600
                        JCOST = KCOST
                        IPIV = I
                        JPIV = J
                        IJPOS = JPOS
                        IF (JCOST.LE.NZ*(NZ-1)) GO TO 820
  600                CONTINUE
  620             CONTINUE
  640          CONTINUE
C
C              MATRIX IS NUMERICALLY OR STRUCTURALLY SINGULAR  ...
C              WHICH IT IS WILL BE DIAGNOSED LATER.
               IRANK = IRANK - 1
  660       CONTINUE
C           ASSIGN REST OF ROWS AND COLUMNS TO ORDERING ARRAY.
C           MATRIX IS STRUCTURALLY SINGULAR.
            IF (IFAIL.NE.-2 .AND. IFAIL.NE.5) IFAIL = -1
            IRANK = IRANK - ILAST + PIVOT + 1
            IF ( .NOT. ABORT(1)) GO TO 680
            IFAIL = 1
            GO TO 2020
  680       K = PIVOT - 1
            DO 760 I = ISTART, ILAST
               IF (LASTR(I).NE.0) GO TO 760
               K = K + 1
               LASTR(I) = K
               IF (LENRL(I).EQ.0) GO TO 740
               MINICN = MAX(MINICN,NZROW+IBEG-1+MOREI+LENRL(I))
               IF (IACTIV-IBEG.GE.LENRL(I)) GO TO 700
               CALL F01BRS(A,ICN,IPTR(ISTART)
     *                     ,N,IACTIV,ITOP,.TRUE.,ICNCP)
C              CHECK NOW TO SEE IF F01BRS HAS CREATED ENOUGH AVAILABLE
C              SPACE.
               IF (IACTIV-IBEG.GE.LENRL(I)) GO TO 700
C              CREATE MORE SPACE BY DESTROYING PREVIOUSLY CREATED LU
C              FACTORS.
               MOREI = MOREI + IBEG - IDISP(1)
               IBEG = IDISP(1)
               IFAIL = 5
               IF (ABORT(3)) GO TO 2020
C              ** CODE FOR OUTPUT OF ERROR MESSAGE **
               IF (MOD(ISAVE/10,10).NE.0) THEN
                  WRITE (REC,FMT=99997)
                  CALL X04BAF(NERR,REC(1))
                  CALL X04BAF(NERR,REC(2))
               END IF
C              ** END OF CODE FOR OUTPUT OF ERROR MESSAGE **
  700          J1 = IPTR(I)
               J2 = J1 + LENRL(I) - 1
               IPTR(I) = 0
               DO 720 JJ = J1, J2
                  A(IBEG) = A(JJ)
                  ICN(IBEG) = ICN(JJ)
                  ICN(JJ) = 0
                  IBEG = IBEG + 1
  720          CONTINUE
               NZROW = NZROW - LENRL(I)
  740          IF (K.EQ.ILAST) GO TO 780
  760       CONTINUE
  780       K = PIVOT - 1
            DO 800 I = ISTART, ILAST
               IF (LASTC(I).NE.0) GO TO 800
               K = K + 1
               LASTC(I) = -K
               IF (K.EQ.ILAST) GO TO 1800
  800       CONTINUE
C
C           THE PIVOT HAS NOW BEEN FOUND IN POSITION (IPIV,JPIV) IN
C           LOCATION IJPOS IN ROW FILE.
C           UPDATE COLUMN AND ROW ORDERING ARRAYS TO CORRESPOND WITH
C           REMOVAL OF THE ACTIVE PART OF THE MATRIX.
  820       ISING = PIVOT
            IF (A(IJPOS).NE.ZERO) GO TO 840
C           NUMERICAL SINGULARITY IS RECORDED HERE.
            ISING = -ISING
            IF (IFAIL.NE.5) IFAIL = -2
            IF ( .NOT. ABORT(2)) GO TO 840
            IFAIL = 2
            GO TO 2020
  840       OLDPIV = IPTR(IPIV) + LENRL(IPIV)
            OLDEND = IPTR(IPIV) + LENR(IPIV) - 1
C           CHANGES TO COLUMN ORDERING.
            DO 880 JJ = OLDPIV, OLDEND
               J = ICN(JJ)
               LC = LASTC(J)
               NC = NEXTC(J)
               IF (NC.NE.0) LASTC(NC) = LC
               IF (LC.EQ.0) GO TO 860
               NEXTC(LC) = NC
               GO TO 880
  860          NZ = LENC(J)
               ISW = IFIRST(NZ)
               IF (ISW.GT.0) LASTR(ISW) = -NC
               IF (ISW.LT.0) IFIRST(NZ) = -NC
  880       CONTINUE
C           CHANGES TO ROW ORDERING.
            I1 = IPC(JPIV)
            I2 = I1 + LENC(JPIV) - 1
            DO 920 II = I1, I2
               I = IRN(II)
               LR = LASTR(I)
               NR = NEXTR(I)
               IF (NR.NE.0) LASTR(NR) = LR
               IF (LR.LE.0) GO TO 900
               NEXTR(LR) = NR
               GO TO 920
  900          NZ = LENR(I) - LENRL(I)
               IF (NR.NE.0) IFIRST(NZ) = NR
               IF (NR.EQ.0) IFIRST(NZ) = LR
  920       CONTINUE
C           RECORD THE COLUMN PERMUTATION IN LASTC(JPIV) AND THE ROW
C           PERMUTATION IN LASTR(IPIV).
            LASTC(JPIV) = ISING
            LASTR(IPIV) = PIVOT
C
C           MOVE PIVOT TO POSITION LENRL+1 IN PIVOT ROW AND MOVE PIVOT
C           ROW TO THE BEGINNING OF THE AVAILABLE STORAGE.
C           THE L PART AND THE PIVOT IN THE OLD COPY OF THE PIVOT ROW IS
C           NULLIFIED WHILE, IN THE STRICTLY UPPER TRIANGULAR PART, THE
C           COLUMN INDICES, J SAY, ARE OVERWRITTEN BY THE CORRESPONDING
C           ELEMENT OF IQ (IQ(J)) AND IQ(J) IS SET TO THE NEGATIVE OF
C           THE DISPLACEMENT OF THE COLUMN INDEX FROM THE PIVOT ELEMENT.
            IF (OLDPIV.EQ.IJPOS) GO TO 940
            AU = A(OLDPIV)
            A(OLDPIV) = A(IJPOS)
            A(IJPOS) = AU
            ICN(IJPOS) = ICN(OLDPIV)
            ICN(OLDPIV) = JPIV
C           CHECK TO SEE IF THERE IS SPACE IMMEDIATELY AVAILABLE IN
C           A/ICN TO HOLD NEW COPY OF PIVOT ROW.
  940       MINICN = MAX(MINICN,NZROW+IBEG-1+MOREI+LENR(IPIV))
            IF (IACTIV-IBEG.GE.LENR(IPIV)) GO TO 960
            CALL F01BRS(A,ICN,IPTR(ISTART),N,IACTIV,ITOP,.TRUE.,ICNCP)
            OLDPIV = IPTR(IPIV) + LENRL(IPIV)
            OLDEND = IPTR(IPIV) + LENR(IPIV) - 1
C           CHECK NOW TO SEE IF F01BRS HAS CREATED ENOUGH AVAILABLE
C           SPACE.
            IF (IACTIV-IBEG.GE.LENR(IPIV)) GO TO 960
C           CREATE MORE SPACE BY DESTROYING PREVIOUSLY CREATED LU
C           FACTORS.
            MOREI = MOREI + IBEG - IDISP(1)
            IBEG = IDISP(1)
            IFAIL = 5
            IF (ABORT(3)) GO TO 2020
C           ** CODE FOR OUTPUT OF ERROR MESSAGE **
            IF (MOD(ISAVE/10,10).NE.0) THEN
               WRITE (REC,FMT=99997)
               CALL X04BAF(NERR,REC(1))
               CALL X04BAF(NERR,REC(2))
            END IF
C           ** END OF CODE FOR OUTPUT OF ERROR MESSAGE **
            IF (IACTIV-IBEG.GE.LENR(IPIV)) GO TO 960
C           THERE IS STILL NOT ENOUGH ROOM IN A/ICN.
            IFAIL = 4
            GO TO 2020
C           COPY PIVOT ROW AND SET UP IQ ARRAY.
  960       IJPOS = 0
            J1 = IPTR(IPIV)
C
            DO 1020 JJ = J1, OLDEND
               A(IBEG) = A(JJ)
               ICN(IBEG) = ICN(JJ)
               IF (IJPOS.NE.0) GO TO 980
               IF (ICN(JJ).EQ.JPIV) IJPOS = IBEG
               ICN(JJ) = 0
               GO TO 1000
  980          K = IBEG - IJPOS
               J = ICN(JJ)
               ICN(JJ) = IQ(J)
               IQ(J) = -K
 1000          IBEG = IBEG + 1
 1020       CONTINUE
C
            IJP1 = IJPOS + 1
            PIVEND = IBEG - 1
            LENPIV = PIVEND - IJPOS
            NZROW = NZROW - LENRL(IPIV) - 1
            IPTR(IPIV) = OLDPIV + 1
            IF (LENPIV.EQ.0) IPTR(IPIV) = 0
C
C           REMOVE PIVOT ROW (INCLUDING PIVOT) FROM COLUMN
C           ORIENTED FILE.
            DO 1080 JJ = IJPOS, PIVEND
               J = ICN(JJ)
               I1 = IPC(J)
               LENC(J) = LENC(J) - 1
C              I2 IS LAST POSITION IN NEW COLUMN.
               I2 = IPC(J) + LENC(J) - 1
               IF (I2.LT.I1) GO TO 1060
               DO 1040 II = I1, I2
                  IF (IRN(II).NE.IPIV) GO TO 1040
                  IRN(II) = IRN(I2+1)
                  GO TO 1060
 1040          CONTINUE
 1060          IRN(I2+1) = 0
 1080       CONTINUE
            NZCOL = NZCOL - LENPIV - 1
C
C           GO DOWN THE PIVOT COLUMN AND FOR EACH ROW WITH A NON-ZERO
C           ADD THE APPROPRIATE MULTIPLE OF THE PIVOT ROW TO IT.
C           WE LOOP ON THE NUMBER OF NON-ZEROS IN THE PIVOT COLUMN SINCE
C           F01BRS MAY CHANGE ITS ACTUAL POSITION.
C
            NZPC = LENC(JPIV)
            IF (NZPC.EQ.0) GO TO 1640
            DO 1520 III = 1, NZPC
               II = IPC(JPIV) + III - 1
               I = IRN(II)
C              SEARCH ROW I FOR NON-ZERO TO BE ELIMINATED, CALCULATE
C              MULTIPLIER, AND PLACE IT IN POSITION LENRL+1 IN ITS ROW.
               J1 = IPTR(I) + LENRL(I)
               IEND = IPTR(I) + LENR(I) - 1
               DO 1100 JJ = J1, IEND
                  IF (ICN(JJ).NE.JPIV) GO TO 1100
C                 IF PIVOT IS ZERO, REST OF COLUMN IS AND SO MULTIPLIER
C                 IS ZERO.
                  AU = ZERO
                  IF (A(IJPOS).NE.ZERO) AU = -A(JJ)/A(IJPOS)
                  A(JJ) = A(J1)
                  A(J1) = AU
                  ICN(JJ) = ICN(J1)
                  ICN(J1) = JPIV
                  LENRL(I) = LENRL(I) + 1
                  GO TO 1120
 1100          CONTINUE
C              IF PIVOT ROW IS A SINGLETON, GO TO ...
 1120          IF (LENPIV.EQ.0) GO TO 1520
C              NOW PERFORM NECESSARY OPERATIONS ON REST OF NON-PIVOT
C              ROW I.
               ROWI = J1 + 1
               IOP = 0
C              IF ALL THE PIVOT ROW CAUSES FILL-IN GO TO 640
               IF (ROWI.GT.IEND) GO TO 1160
C              PERFORM OPERATIONS ON CURRENT NON-ZEROS IN ROW I.
C              INNERMOST LOOP.
               DO 1140 JJ = ROWI, IEND
                  J = ICN(JJ)
                  IF (IQ(J).GT.0) GO TO 1140
                  IOP = IOP + 1
                  PIVROW = IJPOS - IQ(J)
                  A(JJ) = A(JJ) + AU*A(PIVROW)
                  ICN(PIVROW) = -ICN(PIVROW)
 1140          CONTINUE
 1160          IFILL = LENPIV - IOP
C              IF THERE IS NO FILL-IN GO TO 740.
               IF (IFILL.EQ.0) GO TO 1360
C              NOW FOR THE FILL-IN.
               MINICN = MAX(MINICN,MOREI+IBEG-1+NZROW+IFILL+LENR(I))
C              SEE IF THERE IS ROOM FOR FILL-IN.
C              GET MAXIMUM SPACE FOR ROW I IN SITU.
               DO 1180 JDIFF = 1, IFILL
                  JNPOS = IEND + JDIFF
                  IF (JNPOS.GT.LICN) GO TO 1200
                  IF (ICN(JNPOS).NE.0) GO TO 1200
 1180          CONTINUE
C              THERE IS ROOM FOR ALL THE FILL-IN AFTER THE END OF
C              THE ROW SO IT CAN BE LEFT IN SITU.
C              NEXT AVAILABLE SPACE FOR FILL-IN.
               IEND = IEND + 1
               GO TO 1360
C              JMORE SPACES FOR FILL-IN ARE REQUIRED IN FRONT OF ROW.
 1200          JMORE = IFILL - JDIFF + 1
               I1 = IPTR(I)
C              WE NOW LOOK IN FRONT OF THE ROW TO SEE IF THERE IS
C              SPACE FOR THE REST OF THE FILL-IN.
               DO 1220 JDIFF = 1, JMORE
                  JNPOS = I1 - JDIFF
                  IF (JNPOS.LT.IACTIV) GO TO 1240
                  IF (ICN(JNPOS).NE.0) GO TO 1260
 1220          CONTINUE
 1240          JNPOS = I1 - JMORE
               GO TO 1280
C              WHOLE ROW MUST BE MOVED TO THE BEGINNING OF AVAILABLE
C              STORAGE.
 1260          JNPOS = IACTIV - LENR(I) - IFILL
C              IF THERE IS SPACE IMMEDIATELY AVAILABLE FOR THE SHIFTED
C              ROW GO TO ...
 1280          IF (JNPOS.GE.IBEG) GO TO 1320
               CALL F01BRS(A,ICN,IPTR(ISTART)
     *                     ,N,IACTIV,ITOP,.TRUE.,ICNCP)
               I1 = IPTR(I)
               IEND = I1 + LENR(I) - 1
               JNPOS = IACTIV - LENR(I) - IFILL
               IF (JNPOS.GE.IBEG) GO TO 1320
C              NO SPACE AVAILABLE SO TRY TO CREATE SOME BY THROWING AWAY
C              PREVIOUS LU DECOMPOSITION.
               MOREI = MOREI + IBEG - IDISP(1) - LENPIV - 1
               IFAIL = 5
               IF (ABORT(3)) GO TO 2020
C              ** CODE FOR OUTPUT OF ERROR MESSAGE **
               IF (MOD(ISAVE/10,10).NE.0) THEN
                  WRITE (REC,FMT=99997)
                  CALL X04BAF(NERR,REC(1))
                  CALL X04BAF(NERR,REC(2))
               END IF
C              ** END OF CODE FOR OUTPUT OF ERROR MESSAGE **
C              KEEP RECORD OF CURRENT PIVOT ROW.
               IBEG = IDISP(1)
               ICN(IBEG) = JPIV
               A(IBEG) = A(IJPOS)
               IJPOS = IBEG
               DO 1300 JJ = IJP1, PIVEND
                  IBEG = IBEG + 1
                  A(IBEG) = A(JJ)
                  ICN(IBEG) = ICN(JJ)
 1300          CONTINUE
               IJP1 = IJPOS + 1
               PIVEND = IBEG
               IBEG = IBEG + 1
               IF (JNPOS.GE.IBEG) GO TO 1320
C              THIS STILL DOES NOT GIVE ENOUGH ROOM.
               IFAIL = 4
               GO TO 2020
 1320          IACTIV = MIN(IACTIV,JNPOS)
C              MOVE NON-PIVOT ROW I.
               IPTR(I) = JNPOS
               DO 1340 JJ = I1, IEND
                  A(JNPOS) = A(JJ)
                  ICN(JNPOS) = ICN(JJ)
                  JNPOS = JNPOS + 1
                  ICN(JJ) = 0
 1340          CONTINUE
C              FIRST NEW AVAILABLE SPACE.
               IEND = JNPOS
 1360          NZROW = NZROW + IFILL
C              INNERMOST FILL-IN LOOP WHICH ALSO RESETS ICN.
               DO 1500 JJ = IJP1, PIVEND
                  J = ICN(JJ)
                  IF (J.LT.0) GO TO 1480
                  A(IEND) = AU*A(JJ)
                  ICN(IEND) = J
                  IEND = IEND + 1
C
C                 PUT NEW ENTRY IN COLUMN FILE.
                  MINIRN = MAX(MINIRN,NZCOL+LENC(J)+1)
                  JEND = IPC(J) + LENC(J)
                  JROOM = NZPC - III + 1 + LENC(J)
                  IF (JEND.GT.LIRN) GO TO 1380
                  IF (IRN(JEND).EQ.0) GO TO 1460
 1380             IF (JROOM.LT.DISPC) GO TO 1400
C                 COMPRESS COLUMN FILE TO OBTAIN SPACE FOR NEW
C                 COPY OF COLUMN.
                  CALL F01BRS(A,IRN,IPC(ISTART)
     *                        ,N,DISPC,LIRN,.FALSE.,IRNCP)
                  IF (JROOM.LT.DISPC) GO TO 1400
                  JROOM = DISPC - 1
                  IF (JROOM.GE.LENC(J)+1) GO TO 1400
C                 COLUMN FILE IS NOT LARGE ENOUGH.
                  GO TO 1980
C                 COPY COLUMN TO BEGINNING OF FILE.
 1400             JBEG = IPC(J)
                  JEND = IPC(J) + LENC(J) - 1
                  JZERO = DISPC - 1
                  DISPC = DISPC - JROOM
                  IDISPC = DISPC
                  DO 1420 II = JBEG, JEND
                     IRN(IDISPC) = IRN(II)
                     IRN(II) = 0
                     IDISPC = IDISPC + 1
 1420             CONTINUE
                  IPC(J) = DISPC
                  JEND = IDISPC
                  DO 1440 II = JEND, JZERO
                     IRN(II) = 0
 1440             CONTINUE
 1460             IRN(JEND) = I
                  NZCOL = NZCOL + 1
                  LENC(J) = LENC(J) + 1
C                 END OF ADJUSTMENT TO COLUMN FILE.
                  GO TO 1500
C
 1480             ICN(JJ) = -J
 1500          CONTINUE
               LENR(I) = LENR(I) + IFILL
C              END OF SCAN OF PIVOT COLUMN.
 1520       CONTINUE
C
C           REMOVE PIVOT COLUMN FROM COLUMN ORIENTED STORAGE AND UPDATE
C           ROW ORDERING ARRAYS.
            I1 = IPC(JPIV)
            I2 = IPC(JPIV) + LENC(JPIV) - 1
            NZCOL = NZCOL - LENC(JPIV)
            DO 1620 II = I1, I2
               I = IRN(II)
               IRN(II) = 0
               NZ = LENR(I) - LENRL(I)
               IF (NZ.NE.0) GO TO 1540
               LASTR(I) = 0
               GO TO 1620
 1540          IFIR = IFIRST(NZ)
               IFIRST(NZ) = I
               IF (IFIR) 1560, 1600, 1580
 1560          LASTR(I) = IFIR
               NEXTR(I) = 0
               GO TO 1620
 1580          LASTR(I) = LASTR(IFIR)
               NEXTR(I) = IFIR
               LASTR(IFIR) = I
               GO TO 1620
 1600          LASTR(I) = 0
               NEXTR(I) = 0
               NZMIN = MIN(NZMIN,NZ)
 1620       CONTINUE
C           RESTORE IQ AND NULLIFY U PART OF OLD PIVOT ROW.
 1640       IPC(JPIV) = 0
            IF (LENPIV.EQ.0) GO TO 1780
            NZROW = NZROW - LENPIV
            JVAL = IJP1
            JZER = IPTR(IPIV)
            IPTR(IPIV) = 0
            DO 1660 JCOUNT = 1, LENPIV
               J = ICN(JVAL)
               IQ(J) = ICN(JZER)
               ICN(JZER) = 0
               JVAL = JVAL + 1
               JZER = JZER + 1
 1660       CONTINUE
C           ADJUST COLUMN ORDERING ARRAYS.
            DO 1760 JJ = IJP1, PIVEND
               J = ICN(JJ)
               NZ = LENC(J)
               IF (NZ.NE.0) GO TO 1680
               LASTC(J) = 0
               GO TO 1760
 1680          IFIR = IFIRST(NZ)
               LASTC(J) = 0
               IF (IFIR) 1700, 1720, 1740
 1700          IFIRST(NZ) = -J
               IFIR = -IFIR
               LASTC(IFIR) = J
               NEXTC(J) = IFIR
               GO TO 1760
 1720          IFIRST(NZ) = -J
               NEXTC(J) = 0
               NZMIN = MIN(NZMIN,NZ)
               GO TO 1760
 1740          LC = -LASTR(IFIR)
               LASTR(IFIR) = -J
               NEXTC(J) = LC
               IF (LC.NE.0) LASTC(LC) = J
 1760       CONTINUE
 1780    CONTINUE
C
C        ********************************************
C        ****    END OF MAIN ELIMINATION LOOP    ****
C        ********************************************
C
C        RESET IACTIV TO POINT TO THE BEGINNING OF THE NEXT BLOCK.
 1800    IF (ILAST.NE.NN) IACTIV = IPTR(ILAST+1)
 1820 CONTINUE
C
C     ********************************************
C     ****    END OF DECOMPOSITION OF BLOCK    ***
C     ********************************************
C
C     RECORD SINGULARITY (IF ANY) IN IQ ARRAY.
      IF (IRANK.EQ.NN) GO TO 1860
      DO 1840 I = 1, NN
         IF (LASTC(I).GT.0) GO TO 1840
         ISING = -LASTC(I)
         IQ(ISING) = -IQ(ISING)
         LASTC(I) = ISING
 1840 CONTINUE
C
C     RUN THROUGH LU DECOMPOSITION CHANGING COLUMN INDICES TO THAT
C     OF NEW ORDER AND PERMUTING LENR AND LENRL ARRAYS ACCORDING TO
C     PIVOT PERMUTATIONS.
 1860 ISTART = IDISP(1)
      IEND = IBEG - 1
      DO 1880 JJ = ISTART, IEND
         JOLD = ICN(JJ)
         ICN(JJ) = LASTC(JOLD)
 1880 CONTINUE
      DO 1900 II = 1, NN
         I = LASTR(II)
         NEXTR(I) = LENR(II)
         NEXTC(I) = LENRL(II)
 1900 CONTINUE
      DO 1920 I = 1, NN
         LENRL(I) = NEXTC(I)
         LENR(I) = NEXTR(I)
 1920 CONTINUE
C
C     UPDATE PERMUTATION ARRAYS IP AND IQ.
      DO 1940 II = 1, NN
         I = LASTR(II)
         J = LASTC(II)
         NEXTR(I) = ABS(IP(II))
         NEXTC(J) = ABS(IQ(II))
 1940 CONTINUE
      DO 1960 I = 1, NN
         IF (IP(I).LT.0) NEXTR(I) = -NEXTR(I)
         IP(I) = NEXTR(I)
         IF (IQ(I).LT.0) NEXTC(I) = -NEXTC(I)
         IQ(I) = NEXTC(I)
 1960 CONTINUE
      IP(NN) = ABS(IP(NN))
      IDISP(2) = IEND
      IF (IFAIL.NE.5) GO TO 2060
C     ** CODE FOR OUTPUT OF ERROR MESSAGE **************************
      IF (MOD(ISAVE/10,10).NE.0) THEN
         WRITE (REC,FMT=99991) MINICN
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
         CALL X04BAF(NERR,REC(3))
      END IF
C     ** END OF CODE FOR OUTPUT OF ERROR MESSAGE *******************
      GO TO 2060
C
C     ***    ERROR RETURNS    ***
C     LIRN TOO SMALL - IFAIL = 3 OR 6
 1980 IF (IFAIL.EQ.5) GO TO 2000
      IFAIL = 3
      GO TO 2020
 2000 IFAIL = 6
C
 2020 CONTINUE
C     ** CODE FOR OUTPUT OF ERROR MESSAGES *************************
      IF (MOD(ISAVE/10,10).EQ.0) GO TO 2040
      IF (IFAIL.EQ.1) WRITE (REC,FMT=99999)
      IF (IFAIL.EQ.2) WRITE (REC,FMT=99998)
      IF (IFAIL.EQ.3) WRITE (REC,FMT=99996)
      IF (IFAIL.EQ.4) WRITE (REC,FMT=99995)
      IF (IFAIL.EQ.5) WRITE (REC,FMT=99994)
      IF (IFAIL.EQ.6) WRITE (REC,FMT=99993)
      IF (IFAIL.GE.1 .AND. IFAIL.LE.6) THEN
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
      END IF
      IF (IFAIL.LE.2) GO TO 2040
      PIVOT = PIVOT - ISTART + 1
      WRITE (REC,FMT=99992) PIVOT, NBLOCK, ISTART, ILAST
      CALL X04BAF(NERR,REC(1))
      CALL X04BAF(NERR,REC(2))
      IF (PIVOT.EQ.0) THEN
         WRITE (REC,FMT=99990) MINIRN
         CALL X04BAF(NERR,REC(1))
      END IF
C     ** END OF CODE FOR OUTPUT OF ERROR MESSAGES ******************
 2040 IDISP(2) = IACTIV
C
 2060 IDISP(3) = IRNCP
      IDISP(4) = ICNCP
      IDISP(5) = IRANK
      IDISP(6) = MINIRN
      IDISP(7) = MINICN
      RETURN
C
99999 FORMAT (/' MATRIX IS STRUCTURALLY SINGULAR - DECOMPOSITION ABORT',
     *  'ED')
99998 FORMAT (/' MATRIX IS NUMERICALLY SINGULAR - DECOMPOSITION ABORTED'
     *  )
99997 FORMAT (/' LU DECOMPOSITION DESTROYED TO CREATE MORE SPACE')
99996 FORMAT (/' LIRN TOO SMALL -')
99995 FORMAT (/' LICN MUCH TOO SMALL -')
99994 FORMAT (/' LICN TOO SMALL -')
99993 FORMAT (/' LICN AND LIRN TOO SMALL -')
99992 FORMAT ('  DECOMPOSITION ABORTED AT STAGE ',I5,' IN BLOCK ',I5,
     *  /'  WITH FIRST ROW ',I5,' AND LAST ROW ',I5)
99991 FORMAT (/' LICN TOO SMALL -',/'  FOR SUCCESSFUL DECOMPOSITION SE',
     *  'T LICN TO AT LEAST ',I8)
99990 FORMAT ('  TO CONTINUE SET LIRN TO AT LEAST ',I8)
      END
