      INTEGER FUNCTION IQAMAX(N,DX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      REAL(16) DX(*)
*     ..
*
*  Purpose
*  =======
*
*     finds the index of element having max. absolute value.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      REAL(16) QMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS
*     ..
      IQAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IQAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      QMAX = ABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (ABS(DX(IX)).LE.QMAX) GO TO 5
          IQAMAX = I
          QMAX = ABS(DX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 QMAX = ABS(DX(1))
      DO 30 I = 2,N
          IF (ABS(DX(I)).LE.QMAX) GO TO 30
          IQAMAX = I
          QMAX = ABS(DX(I))
   30 CONTINUE
      RETURN
      END
