      SUBROUTINE F01BRW(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,PREORD,ARP,CV,
     *                  OUT,NP1)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DERIVED FROM HARWELL LIBRARY ROUTINE MC21B
C
C     FINDS A ROW PERMUTATION OF A SPARSE MATRIX THAT MAKES ALL THE
C     DIAGONAL ELEMENTS NON-ZERO (OR AS MANY OF THEM AS POSSIBLE),
C     GIVEN THE PATTERN OF NON-ZEROS.
C
C     PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C     .     IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C     ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF
C     .     THE ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.
C     .     IN WHICH CASE (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ
C     .     ENTRIES.
C     CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C     .     WAS VISITED.
C     PREORD(I) IS THE ROW IN THE ORIGINAL MATRIX WHICH HAS THE
C     .     ITH SMALLEST NUMBER OF NON-ZEROS.
C     ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C     .     WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP
C     .     ASSIGNMENT.
C     OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C     .     WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE
C     .     MAIN LOOP.
C
C
C     INITIALIZATION OF ARRAYS.
C     .. Scalar Arguments ..
      INTEGER           LICN, N, NP1, NUMNZ
C     .. Array Arguments ..
      INTEGER           ARP(N), CV(NP1), ICN(LICN), IP(N), IPERM(N),
     *                  LENR(N), OUT(N), PR(N), PREORD(N)
C     .. Local Scalars ..
      INTEGER           I, II, IN1, IN2, IOUTK, IPTR, IROW, J, J1, JORD,
     *                  K, KK, NUM
C     .. Executable Statements ..
      DO 20 I = 1, N
         ARP(I) = LENR(I) - 1
         CV(I) = 0
         IPERM(I) = 0
   20 CONTINUE
      CV(NP1) = 0
      NUMNZ = 0
C
C     ROWS ORDERED IN ORDER OF INCREASING NUMBER OF NON-ZEROS,
C     USING  LIST SORT INVOLVING O(N) OPERATIONS.
      DO 40 I = 1, N
         IROW = LENR(I) + 1
         PR(I) = CV(IROW)
         CV(IROW) = I
   40 CONTINUE
      NUM = 1
      DO 100 I = 1, NP1
         IPTR = CV(I)
         DO 60 J = 1, N
            IF (IPTR.EQ.0) GO TO 80
            PREORD(NUM) = IPTR
            NUM = NUM + 1
            IPTR = PR(IPTR)
   60    CONTINUE
   80    CV(I) = 0
  100 CONTINUE
C
C     MAIN LOOP.
C     EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C     OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 280 JORD = 1, N
         J = PREORD(JORD)
         PR(J) = -1
         DO 220 K = 1, JORD
C           LOOK FOR A CHEAP ASSIGNMENT
            IN1 = ARP(J)
            IF (IN1.LT.0) GO TO 140
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
            DO 120 II = IN1, IN2
               I = ICN(II)
               IF (IPERM(I).EQ.0) GO TO 240
  120       CONTINUE
C           NO CHEAP ASSIGNMENT IN ROW.
            ARP(J) = -1
C           BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
  140       OUT(J) = LENR(J) - 1
C           INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
            DO 200 KK = 1, JORD
               IN1 = OUT(J)
               IF (IN1.LT.0) GO TO 180
               IN2 = IP(J) + LENR(J) - 1
               IN1 = IN2 - IN1
C              FORWARD SCAN.
               DO 160 II = IN1, IN2
                  I = ICN(II)
                  IF (CV(I).EQ.JORD) GO TO 160
C                 COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
                  J1 = J
                  J = IPERM(I)
                  CV(I) = JORD
                  PR(J) = J1
                  OUT(J1) = IN2 - II - 1
                  GO TO 220
  160          CONTINUE
C
C              BACKTRACKING STEP.
  180          J = PR(J)
               IF (J.EQ.-1) GO TO 280
  200       CONTINUE
  220    CONTINUE
C
C        NEW ASSIGNMENT IS MADE.
  240    IPERM(I) = J
         ARP(J) = IN2 - II - 1
         NUMNZ = NUMNZ + 1
         DO 260 K = 1, JORD
            J = PR(J)
            IF (J.EQ.-1) GO TO 280
            II = IP(J) + LENR(J) - OUT(J) - 2
            I = ICN(II)
            IPERM(I) = J
  260    CONTINUE
  280 CONTINUE
C
C     IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C     PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) RETURN
      DO 300 I = 1, N
         ARP(I) = 0
  300 CONTINUE
      K = 0
      DO 340 I = 1, N
         IF (IPERM(I).NE.0) GO TO 320
         K = K + 1
         OUT(K) = I
         GO TO 340
  320    J = IPERM(I)
         ARP(J) = I
  340 CONTINUE
      K = 0
      DO 360 I = 1, N
         IF (ARP(I).NE.0) GO TO 360
         K = K + 1
         IOUTK = OUT(K)
         IPERM(IOUTK) = I
  360 CONTINUE
      RETURN
      END
