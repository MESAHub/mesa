      SUBROUTINE F01BRX(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DERIVED FROM HARWELL LIBRARY ROUTINE MC21A
C
C     INTERFACE FOR F01BRW
C
C     .. Scalar Arguments ..
      INTEGER           LICN, N, NUMNZ
C     .. Array Arguments ..
      INTEGER           ICN(LICN), IP(N), IPERM(N), IW(N,5), LENR(N)
C     .. Local Scalars ..
      INTEGER           NP1
C     .. External Subroutines ..
      EXTERNAL          F01BRW
C     .. Executable Statements ..
      NP1 = N + 1
      CALL F01BRW(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),IW(1,3)
     *            ,IW(1,4),IW(1,5),NP1)
      RETURN
      END
