      SUBROUTINE   DAXPYI   ( NZ, A, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  DAXPYI -- INDEXED DOUBLE PRECISION ELEMENTARY           ====
C     ====            VECTOR OPERATION                              ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DAXPYI ADDS A DOUBLE PRECISION SCALAR MULTIPLE OF 
C             A DOUBLE PRECISION SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         TO  
C             A DOUBLE PRECISION VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST BE 
C         DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         A       DOUBLE      SCALAR MULTIPLIER OF  X.
C         X       DOUBLE      ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C
C     UPDATED ...
C
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ON OUTPUT
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN MODIFIED.
C
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C         
      DOUBLE PRECISION    Y (*), X (*), A
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      IF  ( A .EQ. 0.0D0 )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I))  = Y(INDX(I)) + A * X(I)
   10 CONTINUE
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION   DDOTI   ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  DDOTI -- DOUBLE PRECISION INDEXED DOT PRODUCT           ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DDOTI COMPUTES THE VECTOR INNER PRODUCT OF 
C             A DOUBLE PRECISION SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         WITH 
C             A DOUBLE PRECISION VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF Y WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       DOUBLE      ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         DDOTI   DOUBLE      DOUBLE PRECISION FUNCTION VALUE EQUAL TO
C                             THE VECTOR INNER PRODUCT.  
C                             IF  NZ .LE. 0  DDOTI IS SET TO ZERO.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      DOUBLE PRECISION    X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      DDOTI = 0.0D0
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          DDOTI = DDOTI  +  X(I) * Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE DGTHR ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  DGTHR -- DOUBLE PRECISION GATHER                        ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DGTHR GATHERS THE SPECIFIED ELEMENTS FROM 
C             A DOUBLE PRECISION VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A DOUBLE PRECISION VECTOR  X  IN COMPRESSED FORM (X,INDX).
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO 
C                             COMPRESSED FORM.
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS CORRESPONDING TO THE INDICES
C                             IN  INDX  WILL BE ACCESSED.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     OUTPUT ...
C
C         X       DOUBLE      ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
C
      INTEGER             NZ, INDX (*)
C
      DOUBLE PRECISION    Y (*), X (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          X(I) = Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE DGTHRZ ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  DGTHRZ -- DOUBLE PRECISION GATHER AND ZERO              ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DGTHRZ GATHERS THE SPECIFIED ELEMENTS FROM 
C             A DOUBLE PRECISION VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A DOUBLE PRECISION VECTOR  X  IN COMPRESSED FORM  (X,INDX).  
C         FURTHERMORE THE GATHERED ELEMENTS OF  Y  ARE SET TO ZERO.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     UPDATED ...
C
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  THE 
C                             GATHERED COMPONENTS IN  Y  ARE SET TO ZERO.  
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN ACCESSED.
C
C     OUTPUT ...
C
C         X       DOUBLE      ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      DOUBLE PRECISION    Y (*), X (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          X(I)       = Y(INDX(I))
          Y(INDX(I)) = 0.0D0
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE   DROTI   ( NZ, X, INDX, Y, C, S )
C
C     ==================================================================
C     ==================================================================
C     ====  DROTI  --  APPLY INDEXED DOUBLE PRECISION GIVENS        ====
C     ====             ROTATION                                     ====
C     ==================================================================
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
C     PURPOSE
C     -------
C
C         DROTI APPLIES A GIVENS ROTATION TO 
C             A SPARSE VECTOR   X  STORED IN COMPRESSED FORM  (X,INDX) 
C         AND 
C             ANOTHER VECTOR  Y  IN FULL STORAGE FORM.
C
C         DROTI DOES NOT HANDLE FILL-IN IN  X  AND THEREFORE, IT IS
C         ASSUMED THAT ALL NONZERO COMPONENTS OF  Y  ARE LISTED IN 
C         INDX.  ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN
C         INDX  ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST
C         BE DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICIES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C         C,S     DOUBLE      THE TWO SCALARS DEFINING THE GIVENS
C                             ROTATION.
C
C     UPDATED ...
C
C         X       DOUBLE      ARRAY CONTAINING THE VALUES OF THE 
C                             SPARSE VECTOR IN COMPRESSED FORM.
C         Y       DOUBLE      ARRAY WHICH CONTAINS THE VECTOR  Y
C                             IN FULL STORAGE FORM.  ONLY THE
C                             ELEMENTS WHOSE INDICIES ARE LISTED IN
C                             INDX  HAVE BEEN REFERENCED OR MODIFIED.
C
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      DOUBLE PRECISION    X (*), Y (*), C, S
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
      DOUBLE PRECISION    TEMP
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      IF  ( ( C .EQ. 1.0D0 ) .AND. ( S .EQ. 0.0D0 ) )  RETURN
C
      DO 10 I = 1, NZ
          TEMP        = - S * X (I)  +  C * Y (INDX(I))
          X (I)       =   C * X (I)  +  S * Y (INDX(I))
          Y (INDX(I)) = TEMP
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE DSCTR ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  DSCTR -- DOUBLE PRECISION SCATTER                       ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DSCTR SCATTERS THE COMPONENTS OF 
C             A SPARSE VECTOR  X  STORED IN COMPRESSED FORM  (X,INDX) 
C         INTO 
C             SPECIFIED COMPONENTS OF A DOUBLE PRECISION VECTOR  Y  
C             IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE MODIFIED.  THE VALUES IN  INDX  MUST BE DISTINCT TO
C         ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE SCATTERED FROM 
C                             COMPRESSED FORM.
C         X       DOUBLE      ARRAY CONTAINING THE VALUES TO BE 
C                             SCATTERED FROM COMPRESSED FORM INTO FULL 
C                             STORAGE FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE VALUES
C                             TO BE SCATTERED FROM COMPRESSED FORM.  
C                             IT IS ASSUMED THAT THE ELEMENTS IN  INDX 
C                             ARE DISTINCT.
C
C     OUTPUT ...
C
C         Y       DOUBLE      ARRAY WHOSE ELEMENTS SPECIFIED BY  INDX
C                             HAVE BEEN SET TO THE CORRESPONDING 
C                             ENTRIES OF  X.  ONLY THE ELEMENTS  
C                             CORRESPONDING TO THE INDICES IN  INDX
C                             HAVE BEEN MODIFIED.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      DOUBLE PRECISION    X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I)) = X(I)
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE   SAXPYI   ( NZ, A, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  SAXPYI -- INDEXED REAL ELEMENTARY VECTOR OPERATION      ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SAXPYI ADDS A REAL SCALAR MULTIPLE OF 
C             A REAL SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         TO  
C             A REAL VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST BE 
C         DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         A       REAL        SCALAR MULTIPLIER OF  X.
C         X       REAL        ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C
C     UPDATED ...
C
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ON OUTPUT
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN MODIFIED.
C
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C         
      REAL                Y (*), X (*), A
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      IF  ( A .EQ. 0.0E0 )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I))  = Y(INDX(I)) + A * X(I)
   10 CONTINUE
C
      RETURN
      END
      REAL FUNCTION   SDOTI   ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  SDOTI -- REAL INDEXED DOT PRODUCT                       ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SDOTI COMPUTES THE VECTOR INNER PRODUCT OF 
C             A REAL SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         WITH 
C             A REAL VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF Y WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       REAL        ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         SDOTI   REAL        REAL FUNCTION VALUE EQUAL TO THE 
C                             VECTOR INNER PRODUCT.  
C                             IF  NZ .LE. 0  SDOTI IS SET TO ZERO.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      REAL                X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      SDOTI = 0.0E0
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          SDOTI = SDOTI  +  X(I) * Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE SGTHR ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  SGTHR -- REAL GATHER                                    ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SGTHR GATHERS THE SPECIFIED ELEMENTS FROM 
C             A REAL VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A REAL VECTOR  X  IN COMPRESSED FORM (X,INDX).
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO 
C                             COMPRESSED FORM.
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS CORRESPONDING TO THE INDICES
C                             IN  INDX  WILL BE ACCESSED.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     OUTPUT ...
C
C         X       REAL        ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
C
      INTEGER             NZ, INDX (*)
C
      REAL                Y (*), X (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          X(I) = Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE SGTHRZ ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  SGTHRZ -- REAL GATHER AND ZERO                          ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SGTHRZ GATHERS THE SPECIFIED ELEMENTS FROM 
C             A REAL VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A REAL VECTOR  X  IN COMPRESSED FORM  (X,INDX).  
C         FURTHERMORE THE GATHERED ELEMENTS OF  Y  ARE SET TO ZERO.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     UPDATED ...
C
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  THE 
C                             GATHERED COMPONENTS IN  Y  ARE SET TO ZERO.  
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN ACCESSED.
C
C     OUTPUT ...
C
C         X       REAL        ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      REAL                Y (*), X (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          X(I)       = Y(INDX(I))
          Y(INDX(I)) = 0.0E0
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE   SROTI   ( NZ, X, INDX, Y, C, S )
C
C     ==================================================================
C     ==================================================================
C     ====  SROTI  --  APPLY INDEXED REAL GIVENS ROTATION           ====
C     ==================================================================
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
C     PURPOSE
C     -------
C
C         SROTI APPLIES A GIVENS ROTATION TO 
C             A SPARSE VECTOR   X  STORED IN COMPRESSED FORM  (X,INDX) 
C         AND 
C             ANOTHER VECTOR  Y  IN FULL STORAGE FORM.
C
C         SROTI DOES NOT HANDLE FILL-IN IN  X  AND THEREFORE, IT IS
C         ASSUMED THAT ALL NONZERO COMPONENTS OF  Y  ARE LISTED IN 
C         INDX.  ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN
C         INDX  ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST
C         BE DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICIES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C         C,S     REAL        THE TWO SCALARS DEFINING THE GIVENS
C                             ROTATION.
C
C     UPDATED ...
C
C         X       REAL        ARRAY CONTAINING THE VALUES OF THE 
C                             SPARSE VECTOR IN COMPRESSED FORM.
C         Y       REAL        ARRAY WHICH CONTAINS THE VECTOR  Y
C                             IN FULL STORAGE FORM.  ONLY THE
C                             ELEMENTS WHOSE INDICIES ARE LISTED IN
C                             INDX  HAVE BEEN REFERENCED OR MODIFIED.
C
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      REAL                X (*), Y (*), C, S
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
      REAL                TEMP
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      IF  ( ( C .EQ. 1.0E0 ) .AND. ( S .EQ. 0.0E0 ) )  RETURN
C
      DO 10 I = 1, NZ
          TEMP        = - S * X (I)  +  C * Y (INDX(I))
          X (I)       =   C * X (I)  +  S * Y (INDX(I))
          Y (INDX(I)) = TEMP
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE SSCTR ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  SSCTR -- REAL SCATTER                                   ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SSCTR SCATTERS THE COMPONENTS OF 
C             A SPARSE VECTOR  X  STORED IN COMPRESSED FORM  (X,INDX) 
C         INTO 
C             SPECIFIED COMPONENTS OF A REAL VECTOR  Y  
C             IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE MODIFIED.  THE VALUES IN  INDX  MUST BE DISTINCT TO
C         ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE SCATTERED FROM 
C                             COMPRESSED FORM.
C         X       REAL        ARRAY CONTAINING THE VALUES TO BE 
C                             SCATTERED FROM COMPRESSED FORM INTO FULL 
C                             STORAGE FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE VALUES
C                             TO BE SCATTERED FROM COMPRESSED FORM.  
C                             IT IS ASSUMED THAT THE ELEMENTS IN  INDX 
C                             ARE DISTINCT.
C
C     OUTPUT ...
C
C         Y       REAL        ARRAY WHOSE ELEMENTS SPECIFIED BY  INDX
C                             HAVE BEEN SET TO THE CORRESPONDING 
C                             ENTRIES OF  X.  ONLY THE ELEMENTS  
C                             CORRESPONDING TO THE INDICES IN  INDX
C                             HAVE BEEN MODIFIED.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      REAL                X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I)) = X(I)
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE   ZAXPYI   ( NZ, A, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  ZAXPYI -- INDEXED COMPLEX*16 ELEMENTARY VECTOR OPERATION ===
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         ZAXPYI ADDS A COMPLEX*16 SCALAR MULTIPLE OF 
C             A COMPLEX*16 SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         TO  
C             A COMPLEX*16 VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST BE 
C         DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         A       COMPLEX*16  SCALAR MULTIPLIER OF  X.
C         X       COMPLEX*16  ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C
C     UPDATED ...
C
C         Y       COMPLEX*16  ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ON OUTPUT
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN MODIFIED.
C
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C         
      COMPLEX*16          Y (*), X (*), A
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      IF  ( A .EQ. ( 0.0D0, 0.0D0 ) )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I))  = Y(INDX(I)) + A * X(I)
   10 CONTINUE
C
      RETURN
      END
      FUNCTION   ZDOTCI   ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  ZDOTCI -- COMPLEX*16 CONJUGATED INDEXED DOT PRODUCT     ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         ZDOTCI COMPUTES THE CONJUGATED VECTOR INNER PRODUCT OF 
C             A COMPLEX*16 SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         WITH 
C             A COMPLEX*16 VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF Y WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       COMPLEX*16  ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  
C         Y       COMPLEX*16  ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         ZDOTCI   COMPLEX*16 COMPLEX*16 FUNCTION VALUE EQUAL TO THE 
C                             CONJUGATED VECTOR INNER PRODUCT.  
C                             IF  NZ .LE. 0  ZDOTCI IS SET TO ZERO.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      COMPLEX*16          X (*), Y (*), ZDOTCI
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      ZDOTCI = ( 0.0D0, 0.0D0 )
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          ZDOTCI = ZDOTCI  +  CONJG ( X(I) ) * Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
      FUNCTION ZDOTUI ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  ZDOTUI -- COMPLEX*16 UNCONJUGATED INDEXED DOT PRODUCT   ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         ZDOTUI COMPUTES THE UNCONJUGATED VECTOR INNER PRODUCT OF 
C             A COMPLEX*16 SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         WITH 
C             A COMPLEX*16 VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       COMPLEX*16  ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  
C         Y       COMPLEX*16  ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         ZDOTUI   COMPLEX*16 COMPLEX*16 FUNCTION VALUE EQUAL TO THE 
C                             UNCONJUGATED VECTOR INNER PRODUCT.  
C                             IF  NZ .LE. 0  ZDOTCI IS SET TO ZERO.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      COMPLEX*16          X (*), Y (*), ZDOTUI
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      ZDOTUI = ( 0.0D0, 0.0D0 )
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          ZDOTUI = ZDOTUI  +  X(I) * Y(INDX(I))
   10 CONTINUE
C
          RETURN
      END
      SUBROUTINE ZGTHR ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  ZGTHR -- COMPLEX*16 GATHER                              ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         ZGTHR GATHERS THE SPECIFIED ELEMENTS FROM 
C             A COMPLEX*16 VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A COMPLEX*16 VECTOR  X  IN COMPRESSED FORM (X,INDX).
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO 
C                             COMPRESSED FORM.
C         Y       COMPLEX*16  ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS CORRESPONDING TO THE INDICES
C                             IN  INDX  WILL BE ACCESSED.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     OUTPUT ...
C
C         X       COMPLEX*16  ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
C
      INTEGER             NZ, INDX (*)
C
      COMPLEX*16          Y (*), X (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          X(I) = Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE ZGTHRZ ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  ZGTHRZ -- COMPLEX*16 GATHER AND ZERO                    ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         ZGTHRZ GATHERS THE SPECIFIED ELEMENTS FROM 
C             A COMPLEX*16 VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A COMPLEX*16 VECTOR  X  IN COMPRESSED FORM  (X,INDX).  
C         FURTHERMORE THE GATHERED ELEMENTS OF  Y  ARE SET TO ZERO.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     UPDATED ...
C
C         Y       COMPLEX*16  ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  THE 
C                             GATHERED COMPONENTS IN  Y  ARE SET TO ZERO.  
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN ACCESSED.
C
C     OUTPUT ...
C
C         X       COMPLEX*16  ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      COMPLEX*16          Y (*), X (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          X(I)       = Y(INDX(I))
          Y(INDX(I)) = ( 0.0D0, 0.0D0 )
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE ZSCTR ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  ZSCTR -- COMPLEX*16 SCATTER                             ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         ZSCTR SCATTERS THE COMPONENTS OF 
C             A SPARSE VECTOR  X  STORED IN COMPRESSED FORM  (X,INDX) 
C         INTO 
C             SPECIFIED COMPONENTS OF A COMPLEX*16 VECTOR  Y  
C             IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE MODIFIED.  THE VALUES IN  INDX  MUST BE DISTINCT TO
C         ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE SCATTERED FROM 
C                             COMPRESSED FORM.
C         X       COMPLEX*16  ARRAY CONTAINING THE VALUES TO BE 
C                             SCATTERED FROM COMPRESSED FORM INTO FULL 
C                             STORAGE FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE VALUES
C                             TO BE SCATTERED FROM COMPRESSED FORM.  
C                             IT IS ASSUMED THAT THE ELEMENTS IN  INDX 
C                             ARE DISTINCT.
C
C     OUTPUT ...
C
C         Y       COMPLEX*16  ARRAY WHOSE ELEMENTS SPECIFIED BY  INDX
C                             HAVE BEEN SET TO THE CORRESPONDING 
C                             ENTRIES OF  X.  ONLY THE ELEMENTS  
C                             CORRESPONDING TO THE INDICES IN  INDX
C                             HAVE BEEN MODIFIED.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      COMPLEX*16          X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I)) = X(I)
   10 CONTINUE
C
      RETURN
      END
C *****************************************
C Indexed assignment routines for integers.
C *****************************************
      SUBROUTINE IGATHER(N, A, B, INDEX)
C
      INTEGER A(N), B(N), INDEX(N)
C
      DO 10 I = 1, N
        A(I) = B(INDEX(I))
  10  CONTINUE
C
      RETURN
C
      END
C
      SUBROUTINE ISCATTER(N, A, INDEX, B)
C
      INTEGER A(N), B(N), INDEX(N)
C
      DO 10 I = 1, N
        A(INDEX(I)) = B(I)
  10  CONTINUE
C
      RETURN
C
      END
      SUBROUTINE ICOPY(N, DX, INCX, DY, INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C THIS VERSION MODIFIED FROM ONE WRITTEN BY  JACK DONGARRA, LINPACK, 6/17/77.
C RON EASTMAN, 1990.
      INTEGER I, INCX, INCY, IX, IY, M, MP1, N, DX(N), DY(N)
C
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N + 1) * INCX + 1
      IF (INCY .LT. 0) IY = (-N + 1) * INCY + 1
C
      DO 10 I = 1, N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
C
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
C
      DO 30 I = 1, M
        DY(I) = DX(I)
   30 CONTINUE
C
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
C
      DO 50 I = MP1, N, 7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
C
      RETURN
      END
