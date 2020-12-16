      SUBROUTINE F01BRZ(NC,MAXA,A,INUM,JPTR,JNUM)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DERIVED FROM HARWELL LIBRARY ROUTINE MC20A
C
C     SORTS THE NON-ZEROS OF A SPARSE MATRIX FROM ARBITRARY ORDER
C     TO COLUMN ORDER, UNORDERED WITHIN EACH COLUMN.
C
C     .. Scalar Arguments ..
      INTEGER           MAXA, NC
C     .. Array Arguments ..
      DOUBLE PRECISION  A(MAXA)
      INTEGER           INUM(MAXA), JNUM(MAXA), JPTR(NC)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACE, ACEP
      INTEGER           I, ICE, ICEP, J, JA, JB, JCE, JCEP, JDISP, K,
     *                  KR, LOC, NULL
C     .. Executable Statements ..
      JDISP = 0
      NULL = -JDISP
C     CLEAR JPTR
      DO 20 J = 1, NC
         JPTR(J) = 0
   20 CONTINUE
C     COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN.
      DO 40 K = 1, MAXA
         J = JNUM(K) + JDISP
         JPTR(J) = JPTR(J) + 1
   40 CONTINUE
C     SET THE JPTR ARRAY
      K = 1
      DO 60 J = 1, NC
         KR = K + JPTR(J)
         JPTR(J) = K
         K = KR
   60 CONTINUE
C
C     REORDER THE ELEMENTS INTO COLUMN ORDER.  THE ALGORITHM
C     IS AN IN-PLACE SORT AND IS OF ORDER MAXA.
      DO 100 I = 1, MAXA
C        ESTABLISH THE CURRENT ENTRY.
         JCE = JNUM(I) + JDISP
         IF (JCE.EQ.0) GO TO 100
         ACE = A(I)
         ICE = INUM(I)
C        CLEAR THE LOCATION VACATED.
         JNUM(I) = NULL
C        CHAIN FROM CURRENT ENTRY TO STORE ITEMS.
         DO 80 J = 1, MAXA
C           CURRENT ENTRY NOT IN CORRECT POSITION.  DETERMINE CORRECT
C           POSITION TO STORE ENTRY.
            LOC = JPTR(JCE)
            JPTR(JCE) = JPTR(JCE) + 1
C           SAVE CONTENTS OF THAT LOCATION.
            ACEP = A(LOC)
            ICEP = INUM(LOC)
            JCEP = JNUM(LOC)
C           STORE CURRENT ENTRY.
            A(LOC) = ACE
            INUM(LOC) = ICE
            JNUM(LOC) = NULL
C           CHECK IF NEXT CURRENT ENTRY NEEDS TO BE PROCESSED.
            IF (JCEP.EQ.NULL) GO TO 100
C           IT DOES.  COPY INTO CURRENT ENTRY.
            ACE = ACEP
            ICE = ICEP
            JCE = JCEP + JDISP
   80    CONTINUE
  100 CONTINUE
C
C     **      RESET JPTR VECTOR.
      JA = 1
      DO 120 J = 1, NC
         JB = JPTR(J)
         JPTR(J) = JA
         JA = JB
  120 CONTINUE
      RETURN
      END
