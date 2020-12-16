      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)
C
C     OVERWRITE COMPLEX CY WITH COMPLEX  CA*CX + CY.
C     FOR I = 0 TO N-1, REPLACE  CY(LY+I*INCY) WITH CA*CX(LX+I*INCX) +
C       CY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
C       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      COMPLEX CX(N),CY(N),CA
C
      CANORM = ABS(REAL(CA)) + ABS(AIMAG(CA))
      IF(N.LE.0.OR.CANORM.EQ.0.E0) RETURN
      IF(INCX.EQ.INCY.AND.INCX.GT.0) GO TO 20
      KX = 1
      KY = 1
      IF(INCX.LT.0) KX = 1+(1-N)*INCX
      IF(INCY.LT.0) KY = 1+(1-N)*INCY
          DO 10 I = 1,N
          CY(KY) = CY(KY) + CA*CX(KX)
          KX = KX + INCX
          KY = KY + INCY
   10 CONTINUE
      RETURN
   20 CONTINUE
      NS = N*INCX
          DO 30 I=1,NS,INCX
          CY(I) = CA*CX(I) + CY(I)
   30     CONTINUE
      RETURN
      END
      SUBROUTINE CCOPY(N,CX,INCX,CY,INCY)
C
C     COPY COMPLEX CX TO COMPLEX CY.
C     FOR I = 0 TO N-1, COPY CX(LX+I*INCX) TO CY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      COMPLEX CX(N),CY(N)
C
      IF(N .LE. 0)RETURN
      IF(INCX.EQ.INCY.AND.INCX.GT.0) GO TO 20
      KX = 1
      KY = 1
      IF(INCX.LT.0) KX = 1+(1-N)*INCX
      IF(INCY.LT.0) KY = 1+(1-N)*INCY
          DO 10 I = 1,N
          CY(KY) = CX(KX)
          KX = KX + INCX
          KY = KY + INCY
   10 CONTINUE
      RETURN
   20 CONTINUE
      NS = N*INCX
          DO 30 I=1,NS,INCX
          CY(I) = CX(I)
   30     CONTINUE
      RETURN
      END
      COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)
C
C     RETURNS THE DOT PRODUCT FOR COMPLEX CX AND CY, USES CONJUGATE(CX)
C     CDOTC = SUM FOR I = 0 TO N-1 OF CONJ(CX(LX+I*INCX))*CY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      COMPLEX CX(N),CY(N)
C
      CDOTC = (0.,0.)
      IF(N .LE. 0)RETURN
      IF(INCX.EQ.INCY.AND.INCX.GT.0) GO TO 20
      KX = 1
      KY = 1
      IF(INCX.LT.0) KX = 1+(1-N)*INCX
      IF(INCY.LT.0) KY = 1+(1-N)*INCY
          DO 10 I = 1,N
          CDOTC = CDOTC + CONJG(CX(KX))*CY(KY)
          KX = KX + INCX
          KY = KY + INCY
   10     CONTINUE
      RETURN
   20 CONTINUE
      NS = N*INCX
          DO 30 I=1,NS,INCX
          CDOTC = CONJG(CX(I))*CY(I) + CDOTC
   30     CONTINUE
      RETURN
      END
      COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY)
C
C     RETURNS THE DOT PRODUCT FOR COMPLEX CX AND CY, NO CONJUGATION
C     CDOTU = SUM FOR I = 0 TO N-1 OF  CX(LX+I*INCX) * CY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      COMPLEX CX(N),CY(N)
C
      CDOTU = (0.,0.)
      IF(N .LE. 0)RETURN
      IF(INCX.EQ.INCY.AND.INCX.GT.0) GO TO 20
      KX = 1
      KY = 1
      IF(INCX.LT.0) KX = 1+(1-N)*INCX
      IF(INCY.LT.0) KY = 1+(1-N)*INCY
          DO 10 I = 1,N
          CDOTU = CDOTU + CX(KX)*CY(KY)
          KX = KX + INCX
          KY = KY + INCY
   10     CONTINUE
      RETURN
   20 CONTINUE
      NS = N*INCX
          DO 30 I=1,NS,INCX
          CDOTU = CDOTU + CX(I)*CY(I)
   30     CONTINUE
      RETURN
      END
      SUBROUTINE CSCAL(N,CA,CX,INCX)
C
C     REPLACE COMPLEX CX BY COMPLEX CA*CX.
C     FOR I = 0 TO N-1, REPLACE CX(1+I*INCX) WITH  CA * CX(1+I*INCX)
C
      COMPLEX CA,CX(N)
C
      IF(N .LE. 0) RETURN
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          CX(I) = CA*CX(I)
   10     CONTINUE
      RETURN
      END
      SUBROUTINE CSSCAL(N,SA,CX,INCX)
C
C     REPLACE COMPLEX CX BY (SINGLE PRECISION SA) * (COMPLEX CX)
C     FOR I = 0 TO N-1, REPLACE CX(1+I*INCX) WITH  SA * CX(1+I*INCX)
C
      COMPLEX CX(N)
      REAL    SA
C
      IF(N .LE. 0) RETURN
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          CX(I) = SA*CX(I)
   10     CONTINUE
      RETURN
      END
      SUBROUTINE CSWAP(N,CX,INCX,CY,INCY)
C
C     INTERCHANGE COMPLEX CX AND COMPLEX CY
C     FOR I = 0 TO N-1, INTERCHANGE  CX(LX+I*INCX) AND CY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GT. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      COMPLEX CX(N),CY(N),CTEMP
C
      IF(N .LE. 0)RETURN
      IF(INCX.EQ.INCY.AND.INCX.GT.0) GO TO 20
      KX = 1
      KY = 1
      IF(INCX.LT.0) KX = 1+(1-N)*INCX
      IF(INCY.LT.0) KY = 1+(1-N)*INCY
          DO 10 I = 1,N
          CTEMP = CX(KX)
          CX(KX) = CY(KY)
          CY(KY) = CTEMP
          KX = KX + INCX
          KY = KY + INCY
   10 CONTINUE
      RETURN
   20 CONTINUE
      NS = N*INCX
          DO 30 I=1,NS,INCX
          CTEMP = CX(I)
          CX(I) = CY(I)
          CY(I) = CTEMP
   30     CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
C
C     RETURNS SUM OF MAGNITUDES OF DOUBLE PRECISION DX.
C     DASUM = SUM FROM 0 TO N-1 OF DABS(DX(1+I*INCX))
C
      DOUBLE PRECISION DX(N)
      DASUM = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I=1,NS,INCX
          DASUM = DASUM + DABS(DX(I))
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DASUM = DASUM + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
         DASUM = DASUM + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2))
     $   + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C     OVERWRITE DOUBLE PRECISION DY WITH DOUBLE PRECISION DA*DX + DY.
C     FOR I = 0 TO N-1, REPLACE  DY(LY+I*INCY) WITH DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
C       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      DOUBLE PRECISION DX(N),DY(N),DA
      IF(N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
C
C     COPY DOUBLE PRECISION DX TO DOUBLE PRECISION DY.
C     FOR I = 0 TO N-1, COPY DX(LX+I*INCX) TO DY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      DOUBLE PRECISION DX(N),DY(N)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS=N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DX(I)
   70     CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDIFF(DA,DB)
C1    ********************************* DDIFF **************************
C     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
C2
      DOUBLE PRECISION DA,DB
      DDIFF=DA-DB
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY.
C     DDOT = SUM FOR I = 0 TO N-1 OF  DX(LX+I*INCX) * DY(LY+I*INCY)
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      DOUBLE PRECISION DX(N),DY(N)
      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     $   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
     ~
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER          NEXT
      DOUBLE PRECISION   DX(N), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
     ~
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,DC,DS)
C
C     MULTIPLY THE 2 X 2 MATRIX  ( DC DS) TIMES THE 2 X N MATRIX (DX**T)
C                                (-DS DC)                        (DY**T)
C     WHERE **T INDICATES TRANSPOSE.    THE ELEMENTS OF DX ARE IN
C     DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
C     LX = (-INCX)*N, AND SIMILARLY FOR DY USING LY AND INCY.
      DOUBLE PRECISION DX,DY,DC,DS,ZERO,ONE,W,Z
      DIMENSION DX(N),DY(N)
C
      DATA ZERO,ONE/0.D0,1.D0/
      IF(N .LE. 0 .OR. (DS .EQ. ZERO .AND. DC .EQ. ONE)) GO TO 40
      IF(.NOT. (INCX .EQ. INCY .AND. INCX .GT. 0)) GO TO 20
C
           NSTEPS=INCX*N
           DO 10 I=1,NSTEPS,INCX
                W=DX(I)
                Z=DY(I)
                DX(I)=DC*W+DS*Z
                DY(I)=-DS*W+DC*Z
   10           CONTINUE
           GO TO 40
C
   20 CONTINUE
           KX=1
           KY=1
C
           IF(INCX .LT. 0) KX=1-(N-1)*INCX
           IF(INCY .LT. 0) KY=1-(N-1)*INCY
C
           DO 30 I=1,N
                W=DX(KX)
                Z=DY(KY)
                DX(KX)=DC*W+DS*Z
                DY(KY)=-DS*W+DC*Z
                KX=KX+INCX
                KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
C
      RETURN
      END
      SUBROUTINE DROTG(DA,DB,DC,DS)
C
C     DESIGNED BY C.L.LAWSON, JPL, 1977 SEPT 08
C
C
C     CONSTRUCT THE GIVENS TRANSFORMATION
C
C         ( DC  DS )
C     G = (        ) ,    DC**2 + DS**2 = 1 ,
C         (-DS  DC )
C
C     WHICH ZEROS THE SECOND ENTRY OF THE 2-VECTOR  (DA,DB)**T .
C
C     THE QUANTITY R = (+/-)DSQRT(DA**2 + DB**2) OVERWRITES DA IN
C     STORAGE.  THE VALUE OF DB IS OVERWRITTEN BY A VALUE Z WHICH
C     ALLOWS DC AND DS TO BE RECOVERED BY THE FOLLOWING ALGORITHM:
C           IF Z=1  SET  DC=0.D0  AND  DS=1.D0
C           IF DABS(Z) .LT. 1  SET  DC=DSQRT(1-Z**2)  AND  DS=Z
C           IF DABS(Z) .GT. 1  SET  DC=1/Z  AND  DS=DSQRT(1-DC**2)
C
C     NORMALLY, THE SUBPROGRAM DROT(N,DX,INCX,DY,INCY,DC,DS) WILL
C     NEXT BE CALLED TO APPLY THE TRANSFORMATION TO A 2 BY N MATRIX.
C
C ------------------------------------------------------------------
C
      DOUBLE PRECISION  DA, DB, DC, DS, U, V, R
      IF (DABS(DA) .LE. DABS(DB)) GO TO 10
C
C *** HERE DABS(DA) .GT. DABS(DB) ***
C
      U = DA + DA
      V = DB / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF DA
C
      R = DSQRT(.25D0 + V**2) * U
C
C     NOTE THAT DC IS POSITIVE
C
      DC = DA / R
      DS = V * (DC + DC)
      DB = DS
      DA = R
      RETURN
C
C *** HERE DABS(DA) .LE. DABS(DB) ***
C
   10 IF (DB .EQ. 0.D0) GO TO 20
      U = DB + DB
      V = DA / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF DB
C     (R IS IMMEDIATELY STORED IN DA)
C
      DA = DSQRT(.25D0 + V**2) * U
C
C     NOTE THAT DS IS POSITIVE
C
      DS = DB / DA
      DC = V * (DS + DS)
      IF (DC .EQ. 0.D0) GO TO 15
      DB = 1.D0 / DC
      RETURN
   15 DB = 1.D0
      RETURN
C
C *** HERE DA = DB = 0.D0 ***
C
   20 DC = 1.D0
      DS = 0.D0
      RETURN
C
      END
      SUBROUTINE DROTM (N,DX,INCX,DY,INCY,DPARAM)
C
C     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
C
C     (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
C     (DY**T)
C
C     DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
C     LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
C     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
C
C     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
C
C       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
C     H=(          )    (          )    (          )    (          )
C       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
C     SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
C
      DOUBLE PRECISION DFLAG,DH12,DH22,DX,TWO,Z,DH11,DH21,
     1 DPARAM,DY,W,ZERO
      DIMENSION DX(N),DY(N),DPARAM(5)
      DATA ZERO,TWO/0.D0,2.D0/
C
      DFLAG=DPARAM(1)
      IF(N .LE. 0 .OR.(DFLAG+TWO.EQ.ZERO)) GO TO 140
          IF(.NOT.(INCX.EQ.INCY.AND. INCX .GT.0)) GO TO 70
C
               NSTEPS=N*INCX
               IF(DFLAG) 50,10,30
   10          CONTINUE
               DH12=DPARAM(4)
               DH21=DPARAM(3)
                    DO 20 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W+Z*DH12
                    DY(I)=W*DH21+Z
   20               CONTINUE
               GO TO 140
   30          CONTINUE
               DH11=DPARAM(2)
               DH22=DPARAM(5)
                    DO 40 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z
                    DY(I)=-W+DH22*Z
   40               CONTINUE
               GO TO 140
   50          CONTINUE
               DH11=DPARAM(2)
               DH12=DPARAM(4)
               DH21=DPARAM(3)
               DH22=DPARAM(5)
                    DO 60 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z*DH12
                    DY(I)=W*DH21+Z*DH22
   60               CONTINUE
               GO TO 140
   70     CONTINUE
          KX=1
          KY=1
          IF(INCX .LT. 0) KX=1+(1-N)*INCX
          IF(INCY .LT. 0) KY=1+(1-N)*INCY
C
          IF(DFLAG)120,80,100
   80     CONTINUE
          DH12=DPARAM(4)
          DH21=DPARAM(3)
               DO 90 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W+Z*DH12
               DY(KY)=W*DH21+Z
               KX=KX+INCX
               KY=KY+INCY
   90          CONTINUE
          GO TO 140
  100     CONTINUE
          DH11=DPARAM(2)
          DH22=DPARAM(5)
               DO 110 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z
               DY(KY)=-W+DH22*Z
               KX=KX+INCX
               KY=KY+INCY
  110          CONTINUE
          GO TO 140
  120     CONTINUE
          DH11=DPARAM(2)
          DH12=DPARAM(4)
          DH21=DPARAM(3)
          DH22=DPARAM(5)
               DO 130 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z*DH12
               DY(KY)=W*DH21+Z*DH22
               KX=KX+INCX
               KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
          RETURN
          END
      SUBROUTINE DROTMG (DD1,DD2,DX1,DY1,DPARAM)
C
C     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
C     THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*
C     DY2)**T.
C     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
C
C     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
C
C       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
C     H=(          )    (          )    (          )    (          )
C       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
C     LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
C     RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
C     VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
C
C     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
C     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
C     OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
C
      DOUBLE PRECISION GAM,ONE,RGAMSQ,DD2,DH11,DH21,DPARAM,DP2,
     1 DQ2,DU,DY1,ZERO,GAMSQ,DD1,DFLAG,DH12,DH22,DP1,DQ1,
     2 DTEMP,DX1,TWO
      DIMENSION DPARAM(5)
C
      DATA ZERO,ONE,TWO /0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
      IF(.NOT. DD1 .LT. ZERO) GO TO 10
C       GO ZERO-H-D-AND-DX1..
          GO TO 60
   10 CONTINUE
C     CASE-DD1-NONNEGATIVE
      DP2=DD2*DY1
      IF(.NOT. DP2 .EQ. ZERO) GO TO 20
          DFLAG=-TWO
          GO TO 260
C     REGULAR-CASE..
   20 CONTINUE
      DP1=DD1*DX1
      DQ2=DP2*DY1
      DQ1=DP1*DX1
C
      IF(.NOT. DABS(DQ1) .GT. DABS(DQ2)) GO TO 40
          DH21=-DY1/DX1
          DH12=DP2/DP1
C
          DU=ONE-DH12*DH21
C
          IF(.NOT. DU .LE. ZERO) GO TO 30
C         GO ZERO-H-D-AND-DX1..
               GO TO 60
   30     CONTINUE
               DFLAG=ZERO
               DD1=DD1/DU
               DD2=DD2/DU
               DX1=DX1*DU
C         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF(.NOT. DQ2 .LT. ZERO) GO TO 50
C         GO ZERO-H-D-AND-DX1..
               GO TO 60
   50     CONTINUE
               DFLAG=ONE
               DH11=DP1/DP2
               DH22=DX1/DY1
               DU=ONE+DH11*DH22
               DTEMP=DD2/DU
               DD2=DD1/DU
               DD1=DTEMP
               DX1=DY1*DU
C         GO SCALE-CHECK
               GO TO 100
C     PROCEDURE..ZERO-H-D-AND-DX1..
   60 CONTINUE
          DFLAG=-ONE
          DH11=ZERO
          DH12=ZERO
          DH21=ZERO
          DH22=ZERO
C
          DD1=ZERO
          DD2=ZERO
          DX1=ZERO
C         RETURN..
          GO TO 220
C     PROCEDURE..FIX-H..
   70 CONTINUE
      IF(.NOT. DFLAG .GE. ZERO) GO TO 90
C
          IF(.NOT. DFLAG .EQ. ZERO) GO TO 80
          DH11=ONE
          DH22=ONE
          DFLAG=-ONE
          GO TO 90
   80     CONTINUE
          DH21=-ONE
          DH12=ONE
          DFLAG=-ONE
   90 CONTINUE
      GO TO IGO,(120,150,180,210)
C     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF(.NOT. DD1 .LE. RGAMSQ) GO TO 130
               IF(DD1 .EQ. ZERO) GO TO 160
               ASSIGN 120 TO IGO
C              FIX-H..
               GO TO 70
  120          CONTINUE
               DD1=DD1*GAM**2
               DX1=DX1/GAM
               DH11=DH11/GAM
               DH12=DH12/GAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF(.NOT. DD1 .GE. GAMSQ) GO TO 160
               ASSIGN 150 TO IGO
C              FIX-H..
               GO TO 70
  150          CONTINUE
               DD1=DD1/GAM**2
               DX1=DX1*GAM
               DH11=DH11*GAM
               DH12=DH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF(.NOT. DABS(DD2) .LE. RGAMSQ) GO TO 190
               IF(DD2 .EQ. ZERO) GO TO 220
               ASSIGN 180 TO IGO
C              FIX-H..
               GO TO 70
  180          CONTINUE
               DD2=DD2*GAM**2
               DH21=DH21/GAM
               DH22=DH22/GAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF(.NOT. DABS(DD2) .GE. GAMSQ) GO TO 220
               ASSIGN 210 TO IGO
C              FIX-H..
               GO TO 70
  210          CONTINUE
               DD2=DD2/GAM**2
               DH21=DH21*GAM
               DH22=DH22*GAM
          GO TO 200
  220 CONTINUE
          IF(DFLAG)250,230,240
  230     CONTINUE
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               GO TO 260
  240     CONTINUE
               DPARAM(2)=DH11
               DPARAM(5)=DH22
               GO TO 260
  250     CONTINUE
               DPARAM(2)=DH11
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               DPARAM(5)=DH22
  260 CONTINUE
          DPARAM(1)=DFLAG
          RETURN
      END
      SUBROUTINE DSCAL(N,DA,DX,INCX)
C
C     REPLACE DOUBLE PRECISION DX BY DOUBLE PRECISION DA*DX.
C     FOR I = 0 TO N-1, REPLACE DX(1+I*INCX) WITH  DA * DX(1+I*INCX)
C
      DOUBLE PRECISION DA,DX(N)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DSDOT(N,SX,INCX,SY,INCY)
C
C     RETURNS D.P. DOT PRODUCT ACCUMULATED IN D.P., FOR S.P. SX AND SY
C     DSDOT = SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL SX(N),SY(N)
C
      DSDOT = 0.D0
      IF(N .LE. 0)RETURN
      IF(INCX.EQ.INCY.AND.INCX.GT.0) GO TO 20
      KX = 1
      KY = 1
      IF(INCX.LT.0) KX = 1+(1-N)*INCX
      IF(INCY.LT.0) KY = 1+(1-N)*INCY
          DO 10 I = 1,N
          DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
          KX = KX + INCX
          KY = KY + INCY
   10 CONTINUE
      RETURN
   20 CONTINUE
      NS = N*INCX
          DO 30 I=1,NS,INCX
          DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
   30     CONTINUE
      RETURN
      END
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
C
C     INTERCHANGE DOUBLE PRECISION DX AND DOUBLE PRECISION DY.
C     FOR I = 0 TO N-1, INTERCHANGE  DX(LX+I*INCX) AND DY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      DOUBLE PRECISION DX(N),DY(N),DTEMP1,DTEMP2,DTEMP3
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
C
C     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
      NS = N*INCX
        DO 70 I=1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70   CONTINUE
      RETURN
      END
      INTEGER FUNCTION ICAMAX(N,CX,INCX)
C
C      RETURNS THE INDEX OF THE COMPONENT OF CX HAVING THE
C      LARGEST SUM OF MAGNITUDES OF REAL AND IMAGINARY PARTS.
C     ICAMAX = FIRST I, I = 1 TO N, TO MINIMIZE
C        ABS(REAL(CX(1-INCX+I*INCX))) + ABS(IMAG(CX(1-INCX+I*INCX)))
C
      COMPLEX CX(N)
C
      ICAMAX = 0
      IF(N.LE.0) RETURN
      ICAMAX = 1
      IF(N .LE. 1) RETURN
      NS = N*INCX
      II = 1
      SUMMAX = ABS(REAL(CX(1))) + ABS(AIMAG(CX(1)))
          DO 20 I=1,NS,INCX
          SUMRI = ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
          IF(SUMMAX.GE.SUMRI) GO TO 10
          SUMMAX = SUMRI
          ICAMAX = II
   10     II = II + 1
   20     CONTINUE
      RETURN
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF DOUBLE PRECISION DX.
C     IDAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(DX(1-INCX+I*INCX))
C
      DOUBLE PRECISION DX(N),DMAX,XMAG
      IDAMAX = 0
      IF(N.LE.0) RETURN
      IDAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      DMAX = DABS(DX(1))
      NS = N*INCX
      II = 1
          DO 10 I = 1,NS,INCX
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 5
          IDAMAX = II
          DMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = XMAG
   30 CONTINUE
      RETURN
      END
      INTEGER FUNCTION ISAMAX(N,SX,INCX)
C
C     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF SINGLE PRECISION SX.
C     ISAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(SX(1-INCX+I*INCX))
C
      REAL SX(N),SMAX,XMAG
      ISAMAX = 0
      IF(N.LE.0) RETURN
      ISAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      SMAX = ABS(SX(1))
      NS = N*INCX
      II = 1
          DO 10 I=1,NS,INCX
          XMAG = ABS(SX(I))
          IF(XMAG.LE.SMAX) GO TO 5
          ISAMAX = II
          SMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
         XMAG = ABS(SX(I))
         IF(XMAG.LE.SMAX) GO TO 30
         ISAMAX = I
         SMAX = XMAG
   30 CONTINUE
      RETURN
      END
      REAL FUNCTION SASUM(N,SX,INCX)
C
C     RETURNS SUM OF MAGNITUDES OF SINGLE PRECISION SX.
C     SASUM = SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
C
      REAL SX(N)
      SASUM = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I=1,NS,INCX
          SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SASUM = SASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        SASUM = SASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))
     $  + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
C
C     OVERWRITE SINGLE PRECISION SY WITH SINGLE PRECISION SA*SX +SY.
C     FOR I = 0 TO N-1, REPLACE  SY(LY+I*INCY) WITH SA*SX(LX+I*INCX) +
C       SY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
C       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL SX(N),SY(N),SA
      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END
      FUNCTION SCASUM(N,CX,INCX)
C     RETURNS SUMS OF MAGNITUDES OF REAL AND IMAGINARY PARTS OF
C     COMPONENTS OF CX.  NOTE THAT THIS IS NOT THE L1 NORM OF CX.
C     CASUM = SUM FROM 0 TO N-1 OF ABS(REAL(CX(1+I*INCX))) +
C             ABS(IMAG(CX(1+I*INCX)))
C
      COMPLEX CX(N)
C
      SCASUM=0.
      IF(N .LE. 0) RETURN
      NS = N*INCX
          DO 10 I=1,NS,INCX
          SCASUM = SCASUM + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
   10     CONTINUE
      RETURN
      END
      REAL FUNCTION SCNRM2( N, CX, INCX)
      LOGICAL IMAG, SCALE
      INTEGER          NEXT
      REAL         CUTLO, CUTHI, HITEST, SUM, XMAX, ABSX, ZERO, ONE
      COMPLEX      CX(N)
      DATA         ZERO, ONE /0.0E0, 1.0E0/
C
C     UNITARY NORM OF THE COMPLEX N-VECTOR STORED IN CX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON , 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C
      IF(N .GT. 0) GO TO 10
         SCNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      DO 210 I=1,NN,INCX
         ABSX = ABS(REAL(CX(I)))
         IMAG = .FALSE.
         GO TO NEXT,(30, 50, 70, 90, 110)
   30 IF( ABSX .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      SCALE = .FALSE.
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( ABSX .EQ. ZERO) GO TO 200
      IF( ABSX .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 ASSIGN 110 TO NEXT
      SUM = (SUM / ABSX) / ABSX
  105 SCALE = .TRUE.
      XMAX = ABSX
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( ABSX .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( ABSX .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / ABSX)**2
         XMAX = ABSX
         GO TO 200
C
  115 SUM = SUM + (ABSX/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
   85 ASSIGN 90 TO NEXT
      SCALE = .FALSE.
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
      HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
   90 IF(ABSX .GE. HITEST) GO TO 100
         SUM = SUM + ABSX**2
  200 CONTINUE
C                  CONTROL SELECTION OF REAL AND IMAGINARY PARTS.
C
      IF(IMAG) GO TO 210
         ABSX = ABS(AIMAG(CX(I)))
         IMAG = .TRUE.
      GO TO NEXT,(  50, 70, 90, 110 )
C
  210 CONTINUE
C
C              END OF MAIN LOOP.
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      SCNRM2 = SQRT(SUM)
      IF(SCALE) SCNRM2 = SCNRM2 * XMAX
  300 CONTINUE
      RETURN
      END
      SUBROUTINE  SCOPY(N,SX,INCX,SY,INCY)
C
C     COPY SINGLE PRECISION SX TO SINGLE PRECISION SY.
C     FOR I = 0 TO N-1, COPY  SX(LX+I*INCX) TO SY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL SX(N),SY(N)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I + 1) = SX(I + 1)
        SY(I + 2) = SX(I + 2)
        SY(I + 3) = SX(I + 3)
        SY(I + 4) = SX(I + 4)
        SY(I + 5) = SX(I + 5)
        SY(I + 6) = SX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SX(I)
   70     CONTINUE
      RETURN
      END
      FUNCTION SDIFF(SA,SB)
C1    ********************************* SDIFF **************************
C     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
C2
      SDIFF=SA-SB
      RETURN
      END
      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
C
C     RETURNS THE DOT PRODUCT OF SINGLE PRECISION SX AND SY.
C     SDOT = SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL SX(N),SY(N)
      SDOT = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SDOT = SDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SDOT = SDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SDOT = SDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     $   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
     ~
   50 CONTINUE
      RETURN
C
C        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        SDOT = SDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END
      REAL FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)
C
C     RETURNS S.P. RESULT WITH DOT PRODUCT ACCUMULATED IN D.P.
C     SDSDOT = SB + SUM FOR I = 0 TO N-1 OF SX(LX+I*INCX)*SY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL              SX(N),SY(N),SB
      DOUBLE PRECISION DSDOT
C
      DSDOT = DBLE(SB)
      IF(N .LE. 0) GO TO 30
      IF(INCX.EQ.INCY.AND.INCX.GT.0) GO TO 40
      KX = 1
      KY = 1
      IF(INCX.LT.0) KX = 1+(1-N)*INCX
      IF(INCY.LT.0) KY = 1+(1-N)*INCY
          DO 10 I = 1,N
          DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
          KX = KX + INCX
          KY = KY + INCY
   10     CONTINUE
   30 SDSDOT = SNGL(DSDOT)
      RETURN
   40 CONTINUE
      NS = N*INCX
          DO 50 I=1,NS,INCX
          DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
   50     CONTINUE
      SDSDOT = SNGL(DSDOT)
      RETURN
      END
      REAL FUNCTION SNRM2 ( N, SX, INCX)
      INTEGER          NEXT
      REAL   SX(N),  CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO, ONE
      DATA   ZERO, ONE /0.0E0, 1.0E0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN SX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C
      IF(N .GT. 0) GO TO 10
         SNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( ABS(SX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( SX(I) .EQ. ZERO) GO TO 200
      IF( ABS(SX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / SX(I)) / SX(I)
  105 XMAX = ABS(SX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( ABS(SX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( ABS(SX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / SX(I))**2
         XMAX = ABS(SX(I))
         GO TO 200
C
  115 SUM = SUM + (SX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(ABS(SX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + SX(J)**2
      SNRM2 = SQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      SNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END
      SUBROUTINE SROT(N,SX,INCX,SY,INCY,SC,SS)
C
C     MULTIPLY THE 2 X 2 MATRIX  ( SC SS) TIMES THE 2 X N MATRIX (SX**T)
C                                (-SS SC)                        (SY**T)
C     WHERE **T INDICATES TRANSPOSE.    THE ELEMENTS OF SX ARE IN
C     SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
C     LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
      REAL             SX,SY,SC,SS,ZERO,ONE,W,Z
      DIMENSION SX(N),SY(N)
C
      DATA ZERO,ONE/0.E0,1.E0/
      IF(N .LE. 0 .OR. (SS .EQ. ZERO .AND. SC .EQ. ONE)) GO TO 40
      IF(.NOT. (INCX .EQ. INCY .AND. INCX .GT. 0)) GO TO 20
C
           NSTEPS=INCX*N
           DO 10 I=1,NSTEPS,INCX
                W=SX(I)
                Z=SY(I)
                SX(I)=SC*W+SS*Z
                SY(I)=-SS*W+SC*Z
   10           CONTINUE
           GO TO 40
C
   20 CONTINUE
           KX=1
           KY=1
C
           IF(INCX .LT. 0) KX=1-(N-1)*INCX
           IF(INCY .LT. 0) KY=1-(N-1)*INCY
C
           DO 30 I=1,N
                W=SX(KX)
                Z=SY(KY)
                SX(KX)=SC*W+SS*Z
                SY(KY)=-SS*W+SC*Z
                KX=KX+INCX
                KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
C
      RETURN
      END
      SUBROUTINE SROTG(SA,SB,SC,SS)
C
C     DESIGNED BY C.L.LAWSON, JPL, 1977 SEPT 08
C
C
C     CONSTRUCT THE GIVENS TRANSFORMATION
C
C         ( SC  SS )
C     G = (        ) ,    SC**2 + SS**2 = 1 ,
C         (-SS  SC )
C
C     WHICH ZEROS THE SECOND ENTRY OF THE 2-VECTOR  (SA,SB)**T .
C
C     THE QUANTITY R = (+/-)SQRT(SA**2 + SB**2) OVERWRITES SA IN
C     STORAGE.  THE VALUE OF SB IS OVERWRITTEN BY A VALUE Z WHICH
C     ALLOWS SC AND SS TO BE RECOVERED BY THE FOLLOWING ALGORITHM:
C           IF Z=1  SET  SC=0.  AND  SS=1.
C           IF ABS(Z) .LT. 1  SET  SC=SQRT(1-Z**2)  AND  SS=Z
C           IF ABS(Z) .GT. 1  SET  SC=1/Z  AND  SS=SQRT(1-SC**2)
C
C     NORMALLY, THE SUBPROGRAM SROT(N,SX,INCX,SY,INCY,SC,SS) WILL
C     NEXT BE CALLED TO APPLY THE TRANSFORMATION TO A 2 BY N MATRIX.
C
C ------------------------------------------------------------------
C
      IF (ABS(SA) .LE. ABS(SB)) GO TO 10
C
C *** HERE ABS(SA) .GT. ABS(SB) ***
C
      U = SA + SA
      V = SB / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF SA
C
      R = SQRT(.25 + V**2) * U
C
C     NOTE THAT SC IS POSITIVE
C
      SC = SA / R
      SS = V * (SC + SC)
      SB = SS
      SA = R
      RETURN
C
C *** HERE ABS(SA) .LE. ABS(SB) ***
C
   10 IF (SB .EQ. 0.) GO TO 20
      U = SB + SB
      V = SA / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF SB
C     (R IS IMMEDIATELY STORED IN SA)
C
      SA = SQRT(.25 + V**2) * U
C
C     NOTE THAT SS IS POSITIVE
C
      SS = SB / SA
      SC = V * (SS + SS)
      IF (SC .EQ. 0.) GO TO 15
      SB = 1. / SC
      RETURN
   15 SB = 1.
      RETURN
C
C *** HERE SA = SB = 0. ***
C
   20 SC = 1.
      SS = 0.
      RETURN
C
      END
      SUBROUTINE SROTM (N,SX,INCX,SY,INCY,SPARAM)
C
C     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
C
C     (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
C     (DX**T)
C
C     SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
C     LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.
C     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
C
C     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
C
C       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
C     H=(          )    (          )    (          )    (          )
C       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
C     SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.
C
      DIMENSION SX(N),SY(N),SPARAM(5)
      DATA ZERO,TWO/0.E0,2.E0/
C
      SFLAG=SPARAM(1)
      IF(N .LE. 0 .OR.(SFLAG+TWO.EQ.ZERO)) GO TO 140
          IF(.NOT.(INCX.EQ.INCY.AND. INCX .GT.0)) GO TO 70
C
               NSTEPS=N*INCX
               IF(SFLAG) 50,10,30
   10          CONTINUE
               SH12=SPARAM(4)
               SH21=SPARAM(3)
                    DO 20 I=1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W+Z*SH12
                    SY(I)=W*SH21+Z
   20               CONTINUE
               GO TO 140
   30          CONTINUE
               SH11=SPARAM(2)
               SH22=SPARAM(5)
                    DO 40 I=1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W*SH11+Z
                    SY(I)=-W+SH22*Z
   40               CONTINUE
               GO TO 140
   50          CONTINUE
               SH11=SPARAM(2)
               SH12=SPARAM(4)
               SH21=SPARAM(3)
               SH22=SPARAM(5)
                    DO 60 I=1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W*SH11+Z*SH12
                    SY(I)=W*SH21+Z*SH22
   60               CONTINUE
               GO TO 140
   70     CONTINUE
          KX=1
          KY=1
          IF(INCX .LT. 0) KX=1+(1-N)*INCX
          IF(INCY .LT. 0) KY=1+(1-N)*INCY
C
          IF(SFLAG)120,80,100
   80     CONTINUE
          SH12=SPARAM(4)
          SH21=SPARAM(3)
               DO 90 I=1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W+Z*SH12
               SY(KY)=W*SH21+Z
               KX=KX+INCX
               KY=KY+INCY
   90          CONTINUE
          GO TO 140
  100     CONTINUE
          SH11=SPARAM(2)
          SH22=SPARAM(5)
               DO 110 I=1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W*SH11+Z
               SY(KY)=-W+SH22*Z
               KX=KX+INCX
               KY=KY+INCY
  110          CONTINUE
          GO TO 140
  120     CONTINUE
          SH11=SPARAM(2)
          SH12=SPARAM(4)
          SH21=SPARAM(3)
          SH22=SPARAM(5)
               DO 130 I=1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W*SH11+Z*SH12
               SY(KY)=W*SH21+Z*SH22
               KX=KX+INCX
               KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
          RETURN
          END
      SUBROUTINE SROTMG (SD1,SD2,SX1,SY1,SPARAM)
C
C     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
C     THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*
C     SY2)**T.
C     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
C
C     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
C
C       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
C     H=(          )    (          )    (          )    (          )
C       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
C     LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22
C     RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
C     VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)
C
C     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
C     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
C     OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
C
      DIMENSION SPARAM(5)
C
      DATA ZERO,ONE,TWO /0.E0,1.E0,2.E0/
      DATA GAM,GAMSQ,RGAMSQ/4096.E0,1.67772E7,5.96046E-8/
      IF(.NOT. SD1 .LT. ZERO) GO TO 10
C       GO ZERO-H-D-AND-SX1..
          GO TO 60
   10 CONTINUE
C     CASE-SD1-NONNEGATIVE
      SP2=SD2*SY1
      IF(.NOT. SP2 .EQ. ZERO) GO TO 20
          SFLAG=-TWO
          GO TO 260
C     REGULAR-CASE..
   20 CONTINUE
      SP1=SD1*SX1
      SQ2=SP2*SY1
      SQ1=SP1*SX1
C
      IF(.NOT. ABS(SQ1) .GT. ABS(SQ2)) GO TO 40
          SH21=-SY1/SX1
          SH12=SP2/SP1
C
          SU=ONE-SH12*SH21
C
          IF(.NOT. SU .LE. ZERO) GO TO 30
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   30     CONTINUE
               SFLAG=ZERO
               SD1=SD1/SU
               SD2=SD2/SU
               SX1=SX1*SU
C         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF(.NOT. SQ2 .LT. ZERO) GO TO 50
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   50     CONTINUE
               SFLAG=ONE
               SH11=SP1/SP2
               SH22=SX1/SY1
               SU=ONE+SH11*SH22
               STEMP=SD2/SU
               SD2=SD1/SU
               SD1=STEMP
               SX1=SY1*SU
C         GO SCALE-CHECK
               GO TO 100
C     PROCEDURE..ZERO-H-D-AND-SX1..
   60 CONTINUE
          SFLAG=-ONE
          SH11=ZERO
          SH12=ZERO
          SH21=ZERO
          SH22=ZERO
C
          SD1=ZERO
          SD2=ZERO
          SX1=ZERO
C         RETURN..
          GO TO 220
C     PROCEDURE..FIX-H..
   70 CONTINUE
      IF(.NOT. SFLAG .GE. ZERO) GO TO 90
C
          IF(.NOT. SFLAG .EQ. ZERO) GO TO 80
          SH11=ONE
          SH22=ONE
          SFLAG=-ONE
          GO TO 90
   80     CONTINUE
          SH21=-ONE
          SH12=ONE
          SFLAG=-ONE
   90 CONTINUE
      GO TO IGO,(120,150,180,210)
C     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF(.NOT. SD1 .LE. RGAMSQ) GO TO 130
               IF(SD1 .EQ. ZERO) GO TO 160
               ASSIGN 120 TO IGO
C              FIX-H..
               GO TO 70
  120          CONTINUE
               SD1=SD1*GAM**2
               SX1=SX1/GAM
               SH11=SH11/GAM
               SH12=SH12/GAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF(.NOT. SD1 .GE. GAMSQ) GO TO 160
               ASSIGN 150 TO IGO
C              FIX-H..
               GO TO 70
  150          CONTINUE
               SD1=SD1/GAM**2
               SX1=SX1*GAM
               SH11=SH11*GAM
               SH12=SH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF(.NOT. ABS(SD2) .LE. RGAMSQ) GO TO 190
               IF(SD2 .EQ. ZERO) GO TO 220
               ASSIGN 180 TO IGO
C              FIX-H..
               GO TO 70
  180          CONTINUE
               SD2=SD2*GAM**2
               SH21=SH21/GAM
               SH22=SH22/GAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF(.NOT. ABS(SD2) .GE. GAMSQ) GO TO 220
               ASSIGN 210 TO IGO
C              FIX-H..
               GO TO 70
  210          CONTINUE
               SD2=SD2/GAM**2
               SH21=SH21*GAM
               SH22=SH22*GAM
          GO TO 200
  220 CONTINUE
          IF(SFLAG)250,230,240
  230     CONTINUE
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               GO TO 260
  240     CONTINUE
               SPARAM(2)=SH11
               SPARAM(5)=SH22
               GO TO 260
  250     CONTINUE
               SPARAM(2)=SH11
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               SPARAM(5)=SH22
  260 CONTINUE
          SPARAM(1)=SFLAG
          RETURN
      END
      SUBROUTINE  SSCAL(N,SA,SX,INCX)
C
C     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX.
C     FOR I = 0 TO N-1, REPLACE SX(1+I*INCX) WITH  SA * SX(1+I*INCX)
C
      REAL SA,SX(N)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          SX(I) = SA*SX(I)
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE  SSWAP (N,SX,INCX,SY,INCY)
C
C     INTERCHANGE SINGLE PRECISION SX AND SINGLE PRECISION SY.
C     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL SX(N),SY(N),STEMP1,STEMP2,STEMP3
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        STEMP1 = SX(IX)
        SX(IX) = SY(IY)
        SY(IY) = STEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        STEMP1 = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        STEMP1 = SX(I)
        STEMP2 = SX(I+1)
        STEMP3 = SX(I+2)
        SX(I) = SY(I)
        SX(I+1) = SY(I+1)
        SX(I+2) = SY(I+2)
        SY(I) = STEMP1
        SY(I+1) = STEMP2
        SY(I+2) = STEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
C
C     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
      NS = N*INCX
        DO 70 I=1,NS,INCX
        STEMP1 = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP1
   70   CONTINUE
      RETURN
      END
