      SUBROUTINE F01BRV(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DERIVED FROM HARWELL LIBRARY ROUTINE MC13D
C
C     INTERFACE FOR F01BRV
C
C     .. Scalar Arguments ..
      INTEGER           LICN, N, NUM
C     .. Array Arguments ..
      INTEGER           IB(N), ICN(LICN), IOR(N), IP(N), IW(N,3),
     *                  LENR(N)
C     .. External Subroutines ..
      EXTERNAL          F01BRU
C     .. Executable Statements ..
      CALL F01BRU(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW(1,1),IW(1,2),IW(1,3))
      RETURN
      END
