      SUBROUTINE X04ABF(I,NADV)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C      IF I = 0, SETS NADV TO CURRENT ADVISORY MESSAGE UNIT NUMBER
C     (STORED IN NADV1).
C     IF I = 1, CHANGES CURRENT ADVISORY MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NADV.
C
C     .. Scalar Arguments ..
      INTEGER           I, NADV
C     .. Local Scalars ..
      INTEGER           NADV1
C     .. Save statement ..
      SAVE              NADV1
C     .. Data statements ..
      DATA              NADV1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NADV = NADV1
      IF (I.EQ.1) NADV1 = NADV
      RETURN
      END
