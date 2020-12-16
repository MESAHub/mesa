      SUBROUTINE F01BRQ(N,ICN,A,LICN,LENR,LENRL,W)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DERIVED FROM HARWELL LIBRARY ROUTINE MC24A
C
C     OBTAINS AN ESTIMATE OF THE LARGEST ELEMENT ENCOUNTERED DURING
C     GAUSSIAN ELIMINATION ON A SPARSE MATRIX FROM THE LU FACTORS
C     OBTAINED.
C
C     .. Scalar Arguments ..
      INTEGER           LICN, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LICN), W(N)
      INTEGER           ICN(LICN), LENR(N), LENRL(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMAXL, AMAXU, WROWL, ZERO
      INTEGER           I, J, J0, J1, J2, JJ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
      AMAXL = ZERO
      DO 20 I = 1, N
         W(I) = ZERO
   20 CONTINUE
      J0 = 1
      DO 120 I = 1, N
         IF (LENR(I).EQ.0) GO TO 120
         J2 = J0 + LENR(I) - 1
         IF (LENRL(I).EQ.0) GO TO 60
C        CALCULATION OF 1-NORM OF L.
         J1 = J0 + LENRL(I) - 1
         WROWL = ZERO
         DO 40 JJ = J0, J1
C           AMAXL IS THE MAXIMUM NORM OF COLUMNS OF L SO FAR FOUND.
            WROWL = WROWL + ABS(A(JJ))
   40    CONTINUE
         AMAXL = MAX(AMAXL,WROWL)
         J0 = J1 + 1
C        CALCULATION OF NORMS OF COLUMNS OF U (MAX-NORMS).
   60    J0 = J0 + 1
         IF (J0.GT.J2) GO TO 100
         DO 80 JJ = J0, J2
            J = ICN(JJ)
            W(J) = MAX(ABS(A(JJ)),W(J))
   80    CONTINUE
  100    J0 = J2 + 1
  120 CONTINUE
C     AMAXU IS SET TO MAXIMUM MAX-NORM OF COLUMNS OF U.
      AMAXU = ZERO
      DO 140 I = 1, N
         AMAXU = MAX(AMAXU,W(I))
  140 CONTINUE
C     GROFAC IS MAX U MAX-NORM TIMES MAX L 1-NORM.
      W(1) = AMAXL*AMAXU
      RETURN
      END
