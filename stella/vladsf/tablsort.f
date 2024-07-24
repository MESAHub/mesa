C --------------------------------------------------------------------------
C                  S U B R O U T I N E   T A B L S O R T
C --------------------------------------------------------------------------
      SUBROUTINE TABLSORT(N,ARRIN,INDX,IORDER)
      ENTRY DTBLSORT(N,ARRIN,INDX,IORDER)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER N, IORDER
C This subroutine is a essentially unmodified version of a routine published
C in 'Numerical Recipes' by Press et. al., which they called INDEXX.F.
C TABLSORT sorts the indecies of ARRIN into ascending order (for IORDER > 0)
C or decending order (for IORDER < 0). When IORDER > 0 (IORDER < 0), INDX(1)
C contains on output the index of the smallest (largest) element in ARRIN,
C INDX(2) contains the index of the second smallest (largest) element in
C ARRIN, etc. N is the length of ARRIN and INDX. The contents of ARRIN are
C left untouched.
C 
      DIMENSION ARRIN(N),INDX(N)
C
      IF (N .LE. 1) THEN
        INDX(1) = 1
        RETURN
      END IF
C
      DO 5 J = 1, N
        INDX(J) = J
   5  CONTINUE
C
C ---------------------------------------------------------------------------
C
C ------------------------------>
      IF (IORDER .GE. 0) THEN
C ------------------------------>
C
        L = N/2 + 1
        IR = N
C----------------
  10    CONTINUE
C----------------
C
        IF ( L .GT. 1) THEN
          L = L-1
          INDXT = INDX(L)
          Q = ARRIN(INDXT)
        ELSE
          INDXT = INDX(IR)
          Q = ARRIN(INDXT)
          INDX(IR) = INDX(1)
          IR = IR-1
C
          IF (IR .EQ. 1) THEN
            INDX(1) = INDXT
C
C           ******
            RETURN
C           ******
C
          END IF
C
        END IF
C
        I = L
        J = L + L
C
C---------------
  20    CONTINUE
C---------------
C
        IF (J .LE. IR) THEN
C
          IF (J .LT. IR) THEN
C
            IF (ARRIN(INDX(J)) .LT. ARRIN(INDX(J + 1))) THEN
              J = J + 1
            END IF
C
          END IF
C
          IF (Q .LT. ARRIN(INDX(J))) THEN
            INDX(I) = INDX(J)
            I = J
            J = J + J
C
          ELSE
            J = IR + 1
          END IF
C
          GO TO 20
        END IF
C
        INDX(I) = INDXT
C
C       >>>>>>>>
        GO TO 10
C       >>>>>>>>
C
C ------------------------------>
      ELSE
C ------------------------------>
        L = N/2 + 1
        IR = N
C----------------
 110    CONTINUE
C----------------
C
        IF ( L .GT. 1) THEN
          L = L-1
          INDXT = INDX(L)
          Q = ARRIN(INDXT)
        ELSE
          INDXT = INDX(IR)
          Q = ARRIN(INDXT)
          INDX(IR) = INDX(1)
          IR = IR-1
C
          IF (IR .EQ. 1) THEN
            INDX(1) = INDXT
C
C           ******
            RETURN
C           ******
C
          END IF
C
        END IF
C
        I = L
        J = L + L
C
C---------------
 120    CONTINUE
C---------------
C
        IF (J .LE. IR) THEN
C
          IF (J .LT. IR) THEN
C
            IF (ARRIN(INDX(J)) .GT. ARRIN(INDX(J + 1))) THEN
              J = J + 1
            END IF
C
          END IF
C
          IF (Q .GT. ARRIN(INDX(J))) THEN
            INDX(I) = INDX(J)
            I = J
            J = J + J
C
          ELSE
            J = IR + 1
          END IF
C
          GO TO 120
        END IF
C
        INDX(I) = INDXT
C
C       >>>>>>>>>
        GO TO 110
C       >>>>>>>>>
C
C ------------------------------>
      END IF
C ------------------------------>
C
      END
C --------------------------------------------------------------------------
C                  S U B R O U T I N E   T A B L S O R T
C --------------------------------------------------------------------------
      SUBROUTINE STBLSORT(N,ARRIN,INDX,IORDER)
      IMPLICIT REAL*4 (A-H, O-Z)
      INTEGER N, IORDER
C This subroutine is a essentially unmodified version of a routine published
C in 'Numerical Recipes' by Press et. al., which they called INDEXX.F.
C TABLSORT sorts the indecies of ARRIN into ascending order (for IORDER > 0)
C or decending order (for IORDER < 0). When IORDER > 0 (IORDER < 0), INDX(1)
C contains on output the index of the smallest (largest) element in ARRIN,
C INDX(2) contains the index of the second smallest (largest) element in
C ARRIN, etc. N is the length of ARRIN and INDX. The contents of ARRIN are
C left untouched.
C 
      DIMENSION ARRIN(N),INDX(N)
C
      IF (N .LE. 1) THEN
        INDX(1) = 1
        RETURN
      END IF
C
      DO 5 J = 1, N
        INDX(J) = J
   5  CONTINUE
C
C ---------------------------------------------------------------------------
C
C ------------------------------>
      IF (IORDER .GE. 0) THEN
C ------------------------------>
C
        L = N/2 + 1
        IR = N
C----------------
  10    CONTINUE
C----------------
C
        IF ( L .GT. 1) THEN
          L = L-1
          INDXT = INDX(L)
          Q = ARRIN(INDXT)
        ELSE
          INDXT = INDX(IR)
          Q = ARRIN(INDXT)
          INDX(IR) = INDX(1)
          IR = IR-1
C
          IF (IR .EQ. 1) THEN
            INDX(1) = INDXT
C
C           ******
            RETURN
C           ******
C
          END IF
C
        END IF
C
        I = L
        J = L + L
C
C---------------
  20    CONTINUE
C---------------
C
        IF (J .LE. IR) THEN
C
          IF (J .LT. IR) THEN
C
            IF (ARRIN(INDX(J)) .LT. ARRIN(INDX(J + 1))) THEN
              J = J + 1
            END IF
C
          END IF
C
          IF (Q .LT. ARRIN(INDX(J))) THEN
            INDX(I) = INDX(J)
            I = J
            J = J + J
C
          ELSE
            J = IR + 1
          END IF
C
          GO TO 20
        END IF
C
        INDX(I) = INDXT
C
C       >>>>>>>>
        GO TO 10
C       >>>>>>>>
C
C ------------------------------>
      ELSE
C ------------------------------>
        L = N/2 + 1
        IR = N
C----------------
 110    CONTINUE
C----------------
C
        IF ( L .GT. 1) THEN
          L = L-1
          INDXT = INDX(L)
          Q = ARRIN(INDXT)
        ELSE
          INDXT = INDX(IR)
          Q = ARRIN(INDXT)
          INDX(IR) = INDX(1)
          IR = IR-1
C
          IF (IR .EQ. 1) THEN
            INDX(1) = INDXT
C
C           ******
            RETURN
C           ******
C
          END IF
C
        END IF
C
        I = L
        J = L + L
C
C---------------
 120    CONTINUE
C---------------
C
        IF (J .LE. IR) THEN
C
          IF (J .LT. IR) THEN
C
            IF (ARRIN(INDX(J)) .GT. ARRIN(INDX(J + 1))) THEN
              J = J + 1
            END IF
C
          END IF
C
          IF (Q .GT. ARRIN(INDX(J))) THEN
            INDX(I) = INDX(J)
            I = J
            J = J + J
C
          ELSE
            J = IR + 1
          END IF
C
          GO TO 120
        END IF
C
        INDX(I) = INDXT
C
C       >>>>>>>>>
        GO TO 110
C       >>>>>>>>>
C
C ------------------------------>
      END IF
C ------------------------------>
C
      END
C --------------------------------------------------------------------------
C                  S U B R O U T I N E   I T A B L S R T
C --------------------------------------------------------------------------
      SUBROUTINE ITABLSRT(N,ARRIN,INDX,IORDER)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER N, IORDER
C This subroutine is a slightly modified version of a routine published
C in 'Numerical Recipes' by Press et. al., which they called INDEXX.F.
C ITABLSRT (similar to TABLSORT) sorts the indecies of ARRIN into 
C ascending order (for IORDER > 0) or decending order (for IORDER < 0).
C When IORDER > 0 (IORDER < 0), INDX(1) contains on output the index of 
C the smallest (largest) element in ARRIN, INDX(2) contains the index 
C of the second smallest (largest) element in ARRIN, etc. N is the 
C length of ARRIN and INDX. The contents of ARRIN are left untouched.
C In this version, ARRIN is assumed to be an integer.
C 
      INTEGER ARRIN(N),INDX(N)
C
      IF (N .LE. 1) THEN
        INDX(1) = 1
        RETURN
      END IF
C
      DO 5 J = 1, N
        INDX(J) = J
   5  CONTINUE
C
C ---------------------------------------------------------------------------
C
C ------------------------------>
      IF (IORDER .GE. 0) THEN
C ------------------------------>
C
        L = N/2 + 1
        IR = N
C----------------
  10    CONTINUE
C----------------
C
        IF ( L .GT. 1) THEN
          L = L-1
          INDXT = INDX(L)
          IQ = ARRIN(INDXT)
        ELSE
          INDXT = INDX(IR)
          IQ = ARRIN(INDXT)
          INDX(IR) = INDX(1)
          IR = IR-1
C
          IF (IR .EQ. 1) THEN
            INDX(1) = INDXT
C
C           ******
            RETURN
C           ******
C
          END IF
C
        END IF
C
        I = L
        J = L + L
C
C---------------
  20    CONTINUE
C---------------
C
        IF (J .LE. IR) THEN
C
          IF (J .LT. IR) THEN
C
            IF (ARRIN(INDX(J)) .LT. ARRIN(INDX(J + 1))) THEN
              J = J + 1
            END IF
C
          END IF
C
          IF (IQ .LT. ARRIN(INDX(J))) THEN
            INDX(I) = INDX(J)
            I = J
            J = J + J
C
          ELSE
            J = IR + 1
          END IF
C
          GO TO 20
        END IF
C
        INDX(I) = INDXT
C
C       >>>>>>>>
        GO TO 10
C       >>>>>>>>
C
C ------------------------------>
      ELSE
C ------------------------------>
        L = N/2 + 1
        IR = N
C----------------
 110    CONTINUE
C----------------
C
        IF ( L .GT. 1) THEN
          L = L-1
          INDXT = INDX(L)
          IQ = ARRIN(INDXT)
        ELSE
          INDXT = INDX(IR)
          IQ = ARRIN(INDXT)
          INDX(IR) = INDX(1)
          IR = IR-1
C
          IF (IR .EQ. 1) THEN
            INDX(1) = INDXT
C
C           ******
            RETURN
C           ******
C
          END IF
C
        END IF
C
        I = L
        J = L + L
C
C---------------
 120    CONTINUE
C---------------
C
        IF (J .LE. IR) THEN
C
          IF (J .LT. IR) THEN
C
            IF (ARRIN(INDX(J)) .GT. ARRIN(INDX(J + 1))) THEN
              J = J + 1
            END IF
C
          END IF
C
          IF (IQ .GT. ARRIN(INDX(J))) THEN
            INDX(I) = INDX(J)
            I = J
            J = J + J
C
          ELSE
            J = IR + 1
          END IF
C
          GO TO 120
        END IF
C
        INDX(I) = INDXT
C
C       >>>>>>>>>
        GO TO 110
C       >>>>>>>>>
C
C ------------------------------>
      END IF
C ------------------------------>
C
      END
