      subroutine do_CSHEP2_sg (N,X,Y,F,NC,NW,NR,LCELL,LNEXT,XMIN,
     &                   YMIN,DX,DY,RMAX,RW,A,IER)
      integer :: N, NC, NW, NR, LCELL(NR,NR), LNEXT(N), IER
      real :: X(N), Y(N), F(N), XMIN, YMIN, DX, DY, RMAX, RW(N), A(9,N)

! **********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/13/97

!   This subroutine computes a set of parameters defining a
! C2 (twice continuously differentiable) bivariate function
! C(X,Y) which interpolates data values F at a set of N
! arbitrarily distributed points (X,Y) in the plane (nodes).
! The interpolant C may be evaluated at an arbitrary point
! by function CS2VAL, and its first partial derivatives are
! computed by Subroutine CS2GRD.

!   The interpolation scheme is a modified Cubic Shepard
! method:

! C = [W(1)*C(1)+W(2)*C(2)+..+W(N)*C(N)]/[W(1)+W(2)+..+W(N)]

! for bivariate functions W(k) and C(k).  The nodal func-
! tions are given by

!  C(k)(x,y) = A(1,k)*(x-X(k))**3 +
!              A(2,k)*(x-X(k))**2*(y-Y(k)) +
!              A(3,k)*(x-X(k))*(y-Y(k))**2 +
!              A(4,k)*(y-Y(k))**3 + A(5,k)*(x-X(k))**2 +
!              A(6,k)*(x-X(k))*(y-Y(k)) + A(7,k)*(y-Y(k))**2
!              + A(8,k)*(x-X(k)) + A(9,k)*(y-Y(k)) + F(k) .

! Thus, C(k) is a cubic function which interpolates the data
! value at node k.  Its coefficients A(,k) are obtained by a
! weighted least squares fit to the closest NC data points
! with weights similar to W(k).  Note that the radius of
! influence for the least squares fit is fixed for each k,
! but varies with k.

! The weights are taken to be

!   W(k)(x,y) = ( (R(k)-D(k))+ / R(k)*D(k) )**3 ,

! where (R(k)-D(k))+ = 0 if R(k) < D(k), and D(k)(x,y) is
! the Euclidean distance between (x,y) and (X(k),Y(k)).  The
! radius of influence R(k) varies with k and is chosen so
! that NW nodes are within the radius.  Note that W(k) is
! not defined at node (X(k),Y(k)), but C(x,y) has limit F(k)
! as (x,y) approaches (X(k),Y(k)).

! On input:

!       N = Number of nodes and data values.  N >= 10.

!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.

!       F = Array of length N containing the data values
!           in one-to-one correspondence with the nodes.

!       NC = Number of data points to be used in the least
!            squares fit for coefficients defining the nodal
!            functions C(k).  Values found to be optimal for
!            test data sets ranged from 11 to 25.  A recom-
!            mended value for general data sets is NC = 17.
!            For nodes lying on (or close to) a rectangular
!            grid, the recommended value is NC = 11.  In any
!            case, NC must be in the range 9 to Min(40,N-1).

!       NW = Number of nodes within (and defining) the radii
!            of influence R(k) which enter into the weights
!            W(k).  For N sufficiently large, a recommended
!            value is NW = 30.  In general, NW should be
!            about 1.5*NC.  1 <= NW <= Min(40,N-1).

!       NR = Number of rows and columns in the cell grid de-
!            fined in Subroutine STORE2_sg.  A rectangle con-
!            taining the nodes is partitioned into cells in
!            order to increase search efficiency.  NR =
!            Sqrt(N/3) is recommended.  NR >= 1.

! The above parameters are not altered by this routine.

!       LCELL = Array of length >= NR**2.

!       LNEXT = Array of length >= N.

!       RW = Array of length >= N.

!       A = Array of length >= 9N.

! On output:

!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2_sg.

!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2_sg.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  Refer to Subroutine
!                         STORE2_sg.

!       RMAX = Largest element in RW -- maximum radius R(k).

!       RW = Array containing the the radii R(k) which enter
!            into the weights W(k).

!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.

!   Note that the output parameters described above are not
! defined unless IER = 0.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NC, NW, or NR is outside its
!                     valid range.
!             IER = 2 if duplicate nodes were encountered.
!             IER = 3 if all nodes are collinear.

! Modules required by CSHEP2:  GETNP2_sg, GIVENS_sg, ROTATE_sg,
!                                SETUP2_sg, STORE2_sg

! Intrinsic functions called by CSHEP2:  ABS, DBLE, MAX,
!                                          MIN, SQRT

! **********************************************************

      integer :: LMX
      parameter (LMX=40)
      integer :: I, IERR, IP1, IRM1, IROW, J, JP1, K, LMAX,
     &        LNP, NEQ, NN, NNC, NNR, NNW, NP, NPTS(LMX),
     &        NCWMAX
      real :: B(10,10), C, DDX, DDY, DMIN, DTOL,
     &                 FK, RC, RS, RSMX, RSOLD, RTOL, RWS,
     &                 S, SF, SFC, SFS, STF, SUM, T, XK,
     &                 XMN, YK, YMN

      DATA    RTOL/1.D-5/, DTOL/.01/

! Local parameters:

! B =          Transpose of the augmented regression matrix
! C =          First component of the plane rotation used to
!                zero the lower triangle of B**T -- computed
!                by Subroutine GIVENS_sg
! DDX,DDY =    Local variables for DX and DY
! DMIN =       Minimum of the magnitudes of the diagonal
!                elements of the regression matrix after
!                zeros are introduced below the diagonal
! DTOL =       Tolerance for detecting an ill-conditioned
!                system.  The system is accepted when
!                DMIN*RC >= DTOL.
! FK =         Data value at mode K -- F(K)
! I =          Index for A, B, and NPTS
! IERR =       Error flag for the call to Subroutine STORE2_sg
! IP1 =        I+1
! IRM1 =       IROW-1
! IROW =       Row index for B
! J =          Index for A and B
! JP1 =        J+1
! K =          Nodal function index and column index for A
! LMAX =       Maximum number of NPTS elements
! LMX =        Maximum value of LMAX
! LNP =        Current length of NPTS
! NEQ =        Number of equations in the least squares fit
! NN,NNC,NNR = Local copies of N, NC, and NR
! NNW =        Local copy of NW
! NP =         NPTS element
! NPTS =       Array containing the indexes of a sequence of
!                nodes to be used in the least squares fit
!                or to compute RW.  The nodes are ordered
!                by distance from K, and the last element
!                (usually indexed by LNP) is used only to
!                determine RC, or RW(K) if NW > NC.
! NCWMAX =     Max(NC,NW)
! RC =         Radius of influence which enters into the
!                weights for C(K) (see Subroutine SETUP2_sg)
! RS =         Squared distance between K and NPTS(LNP) --
!                used to compute RC and RW(K)
! RSMX =       Maximum squared RW element encountered
! RSOLD =      Squared distance between K and NPTS(LNP-1) --
!                used to compute a relative change in RS
!                between succeeding NPTS elements
! RTOL =       Tolerance for detecting a sufficiently large
!                relative change in RS.  If the change is
!                not greater than RTOL, the nodes are
!                treated as being the same distance from K
! RWS =        Current squared value of RW(K)
! S =          Second component of the plane rotation deter-
!                mined by subroutine GIVENS_sg
! SF =        Scale factor for the linear terms (columns 8
!               and 9) in the least squares fit -- inverse
!               of the root-mean-square distance between K
!               and the nodes (other than K) in the least
!               squares fit
! SFS =       Scale factor for the quadratic terms (columns
!               5, 6, and 7) in the least squares fit --
!               SF*SF
! SFC =       Scale factor for the cubic terms (first 4
!               columns) in the least squares fit -- SF**3
! STF =        Marquardt stabilization factor used to damp
!                out the first 4 solution components (third
!                partials of the cubic) when the system is
!                ill-conditioned.  As STF increases, the
!                fitting function approaches a quadratic
!                polynomial.
! SUM =        Sum of squared Euclidean distances between
!                node K and the nodes used in the least
!                squares fit (unless additional nodes are
!                added for stability)
! T =          Temporary variable for accumulating a scalar
!                product in the back solve
! XK,YK =      Coordinates of node K -- X(K), Y(K)
! XMN,YMN =    Local variables for XMIN and YMIN

      NN = N
      NNC = NC
      NNW = NW
      NNR = NR
      NCWMAX = MAX(NNC,NNW)
      LMAX = MIN(LMX,NN-1)
      if (NNC < 9 .or. NNW < 1 .or. NCWMAX >
     &   LMAX .or. NNR < 1) GOTO 21

! Create the cell data structure, and initialize RSMX.

      call STORE2_sg (NN,X,Y,NNR, LCELL,LNEXT,XMN,YMN,DDX,DDY,IERR)
      if (IERR /= 0) GOTO 23
      RSMX = 0.

! Outer loop on node K:

      do K = 1,NN
        XK = X(K)
        YK = Y(K)
        FK = F(K)

! Mark node K to exclude it from the search for nearest
!   neighbors.

        LNEXT(K) = -LNEXT(K)

! Initialize for loop on NPTS.

        RS = 0.
        SUM = 0.
        RWS = 0.
        RC = 0.
        LNP = 0

! Compute NPTS, LNP, RWS, NEQ, RC, and SFS.

    1   SUM = SUM + RS
          if (LNP == LMAX) GOTO 2
          LNP = LNP + 1
          RSOLD = RS
          call GETNP2_sg (XK,YK,X,Y,NNR,LCELL,LNEXT,XMN,YMN,
     &                 DDX,DDY, NP,RS)
          if (RS == 0.) GOTO 22
          NPTS(LNP) = NP
          if ( (RS-RSOLD)/RS < RTOL ) GOTO 1
          if (RWS == 0. .and. LNP > NNW) RWS = RS
          if (RC == 0. .and. LNP > NNC) then

!   RC = 0 (not yet computed) and LNP > NC.  RC = Sqrt(RS)
!     is sufficiently large to (strictly) include NC nodes.
!     The least squares fit will include NEQ = LNP - 1
!     equations for 9 <= NC <= NEQ < LMAX <= N-1.

            NEQ = LNP - 1
            RC = SQRT(RS)
            SFS = DBLE(NEQ)/SUM
          end if

!   Bottom of loop -- test for termination.

          if (LNP > NCWMAX) GOTO 3
          GOTO 1

! All LMAX nodes are included in NPTS.  RWS and/or RC**2 is
!   (arbitrarily) taken to be 10 percent larger than the
!   distance RS to the last node included.

    2   if (RWS == 0.) RWS = 1.1*RS
        if (RC == 0.) then
          NEQ = LMAX
          RC = SQRT(1.1*RS)
          SFS = DBLE(NEQ)/SUM
        end if

! Store RW(K), update RSMX if necessary, and compute SF
!   and SFC.

    3   RW(K) = SQRT(RWS)
        if (RWS > RSMX) RSMX = RWS
        SF = SQRT(SFS)
        SFC = SF*SFS

! A Q-R decomposition is used to solve the least squares
!   system.  The transpose of the augmented regression
!   matrix is stored in B with columns (rows of B) defined
!   as follows:  1-4 are the cubic terms, 5-7 are the quad-
!   ratic terms, 8 and 9 are the linear terms, and the last
!   column is the right hand side.

! Set up the equations and zero out the lower triangle with
!   Givens rotations.

        I = 0
    4     I = I + 1
          NP = NPTS(I)
          IROW = MIN(I,10)
          call SETUP2_sg (XK,YK,FK,X(NP),Y(NP),F(NP),SF,SFS,
     &                 SFC,RC, B(1,IROW))
          if (I == 1) GOTO 4
          IRM1 = IROW-1
          do 5 J = 1,IRM1
            JP1 = J + 1
            call GIVENS_sg (B(J,J),B(J,IROW),C,S)
            call ROTATE_sg (10-J,C,S,B(JP1,J),B(JP1,IROW))
    5       continue
          if (I < NEQ) GOTO 4

! Test the system for ill-conditioning.

        DMIN = MIN( ABS(B(1,1)),ABS(B(2,2)),ABS(B(3,3)),
     &              ABS(B(4,4)),ABS(B(5,5)),ABS(B(6,6)),
     &              ABS(B(7,7)),ABS(B(8,8)),ABS(B(9,9)) )
        if (DMIN*RC >= DTOL) GOTO 11
        if (NEQ == LMAX) GOTO 7

! Increase RC and add another equation to the system to
!   improve the conditioning.  The number of NPTS elements
!   is also increased if necessary.

    6   RSOLD = RS
        NEQ = NEQ + 1
        if (NEQ == LMAX) then
          RC = SQRT(1.1*RS)
          GOTO 4
        end if
        if (NEQ < LNP) then

!   NEQ < LNP.

          NP = NPTS(NEQ+1)
          RS = (X(NP)-XK)**2 + (Y(NP)-YK)**2
          if ( (RS-RSOLD)/RS < RTOL ) GOTO 6
          RC = SQRT(RS)
          GOTO 4
        end if

!   NEQ = LNP.  Add an element to NPTS.

        LNP = LNP + 1
        call GETNP2_sg (XK,YK,X,Y,NNR,LCELL,LNEXT,XMN,YMN,
     &               DDX,DDY, NP,RS)
        if (NP == 0) GOTO 22
        NPTS(LNP) = NP
        if ( (RS-RSOLD)/RS < RTOL ) GOTO 6
        RC = SQRT(RS)
        GOTO 4

! Stabilize the system by damping third partials -- add
!   multiples of the first four unit vectors to the first
!   four equations.

    7   STF = 1.0/RC
        do I = 1,4
          B(I,10) = STF
          IP1 = I + 1
          do J = IP1,10
            B(J,10) = 0.
          end do
          do J = I,9
            JP1 = J + 1
            call GIVENS_sg (B(J,J),B(J,10),C,S)
            call ROTATE_sg (10-J,C,S,B(JP1,J),B(JP1,10))
          end do
        end do

! Test the damped system for ill-conditioning.

        DMIN = MIN( ABS(B(5,5)),ABS(B(6,6)),ABS(B(7,7)),
     &              ABS(B(8,8)),ABS(B(9,9)) )
        if (DMIN*RC < DTOL) GOTO 23

! Solve the 9 by 9 triangular system for the coefficients.

   11   do I = 9,1,-1
          T = 0.
          if (I /= 9) then
            IP1 = I + 1
            do J = IP1,9
              T = T + B(J,I)*A(J,K)
            end do
          end if
          A(I,K) = (B(10,I)-T)/B(I,I)
        end do

! Scale the coefficients to adjust for the column scaling.

        do I = 1,4
          A(I,K) = A(I,K)*SFC
        end do
        A(5,K) = A(5,K)*SFS
        A(6,K) = A(6,K)*SFS
        A(7,K) = A(7,K)*SFS
        A(8,K) = A(8,K)*SF
        A(9,K) = A(9,K)*SF

! Unmark K and the elements of NPTS.

        LNEXT(K) = -LNEXT(K)
        do I = 1,LNP
          NP = NPTS(I)
          LNEXT(NP) = -LNEXT(NP)
        end do
      end do

! No errors encountered.

      XMIN = XMN
      YMIN = YMN
      DX = DDX
      DY = DDY
      RMAX = SQRT(RSMX)
      IER = 0
      return

! N, NC, NW, or NR is outside its valid range.

   21 IER = 1
      return

! Duplicate nodes were encountered by GETNP2_sg.

   22 IER = 2
      return

! No unique solution due to collinear nodes.

   23 XMIN = XMN
      YMIN = YMN
      DX = DDX
      DY = DDY
      IER = 3
      return
      end subroutine do_CSHEP2_sg

      real FUNCTION do_CS2VAL_sg (PX,PY,N,X,Y,F,NR,
     &                LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,IER)
      integer :: N, NR, LCELL(NR,NR), LNEXT(N), IER
      real :: PX, PY, X(N), Y(N), F(N), XMIN, YMIN,
     &                 DX, DY, RMAX, RW(N), A(9,N)

! **********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   This function returns the value C(PX,PY), where C is the
! weighted sum of cubic nodal functions defined in Subrou-
! tine CSHEP2.  CS2GRD may be called to compute a gradient
! of C along with the value, and/or to test for errors.
! CS2HES_sg may be called to compute a value, first partial
! derivatives, and second partial derivatives at a point.

! On input:

!       PX,PY = Cartesian coordinates of the point P at
!               which C is to be evaluated.

!       N = Number of nodes and data values defining C.
!           N >= 10.

!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.

!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2_sg.  NR >= 1.

!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2_sg.

!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2_sg.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2_sg.

!       RMAX = Largest element in RW -- maximum radius R(k).

!       RW = Array containing the the radii R(k) which enter
!            into the weights W(k) defining C.

!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.

!   Input parameters are not altered by this function.  The
! parameters other than PX and PY should be input unaltered
! from their values on output from CSHEP2.  This function
! should not be called if a nonzero error flag was returned
! by CSHEP2.

! On output:

!       CS2VAL = Function value C(PX,PY) unless N, NR, DX,
!                DY, or RMAX is invalid, in which case no
!                value is returned.

! Modules required by CS2VAL:  NONE

! Intrinsic functions called by CS2VAL:  INT, SQRT

! **********************************************************

      integer :: I, IMAX, IMIN, J, JMAX, JMIN, K, KP
      real :: D, DELX, DELY, R, SW, SWC, W, XP, YP

! Local parameters:

! D =         Distance between P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes
!               in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at P
! W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! XP,YP =     Local copies of PX and PY -- coordinates of P

      XP = PX
      YP = PY
      IER = -1
      if (N < 10 .or. NR < 1 .or. DX <= 0.  .or.
     &    DY <= 0. .or. RMAX < 0.) return
      IER = 0

! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
!   the range of the search for nodes whose radii include
!   P.  The cells which must be searched are those inter-
!   sected by (or contained in) a circle of radius RMAX
!   centered at P.

      IMIN = INT((XP-XMIN-RMAX)/DX) + 1
      IMAX = INT((XP-XMIN+RMAX)/DX) + 1
      if (IMIN < 1) IMIN = 1
      if (IMAX > NR) IMAX = NR
      JMIN = INT((YP-YMIN-RMAX)/DY) + 1
      JMAX = INT((YP-YMIN+RMAX)/DY) + 1
      if (JMIN < 1) JMIN = 1
      if (JMAX > NR) JMAX = NR

! The following is a test for no cells within the circle
!   of radius RMAX.

      if (IMIN > IMAX .or. JMIN > JMAX) GOTO 6

! Accumulate weight values in SW and weighted nodal function
!   values in SWC.  The weights are W(K) = ((R-D)+/(R*D))**3
!   for R = RW(K) and D = distance between P and node K.

      SW = 0.
      SWC = 0.

! Outer loop on cells (I,J).

      do 4 J = JMIN,JMAX
        do 3 I = IMIN,IMAX
          K = LCELL(I,J)
          if (K == 0) GOTO 3

! Inner loop on nodes K.

    1     DELX = XP - X(K)
          DELY = YP - Y(K)
          D = SQRT(DELX*DELX + DELY*DELY)
          R = RW(K)
          if (D >= R) GOTO 2
          if (D == 0.) GOTO 5
          W = (1.0/D - 1.0/R)*(1.0/D - 1.0/R)*(1.0/D - 1.0/R)
          SW = SW + W
          SWC = SWC + W*( ( (A(1,K)*DELX+A(2,K)*DELY+
     &                       A(5,K))*DELX + (A(3,K)*DELY+
     &                       A(6,K))*DELY + A(8,K) )*DELX +
     &                    ( (A(4,K)*DELY+A(7,K))*DELY +
     &                      A(9,K) )*DELY + F(K) )

! Bottom of loop on nodes in cell (I,J).

    2     KP = K
          K = LNEXT(KP)
          if (K /= KP) GOTO 1
    3     continue
    4   continue

! SW = 0 iff P is not within the radius R(K) for any node K.

      if (SW == 0.) GOTO 6
      do_CS2VAL_sg = SWC/SW
      return

! (PX,PY) = (X(K),Y(K)).

    5 do_CS2VAL_sg = F(K)
      return

! All weights are 0 at P.

    6 do_CS2VAL_sg = 0.
      return
      end function do_CS2VAL_sg

      subroutine do_CS2GRD_sg (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,
     &                   YMIN,DX,DY,RMAX,RW,A, C,CX,CY,IER)
      integer :: N, NR, LCELL(NR,NR), LNEXT(N), IER
      real :: PX, PY, X(N), Y(N), F(N), XMIN, YMIN,
     &                 DX, DY, RMAX, RW(N), A(9,N), C, CX, CY

! **********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   This subroutine computes the value and gradient at P =
! (PX,PY) of the interpolatory function C defined in Sub-
! routine CSHEP2.  C is a weighted sum of cubic nodal
! functions.

! On input:

!       PX,PY = Cartesian coordinates of the point P at
!               which C and its partial derivatives are
!               to be evaluated.

!       N = Number of nodes and data values defining C.
!           N >= 10.

!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.

!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2_sg.  NR >= 1.

!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2_sg.

!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2_sg.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2_sg.

!       RMAX = Largest element in RW -- maximum radius R(k).

!       RW = Array of length N containing the the radii R(k)
!            which enter into the weights W(k) defining C.

!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.

!   Input parameters are not altered by this subroutine.
! The parameters other than PX and PY should be input
! unaltered from their values on output from CSHEP2.  This
! subroutine should not be called if a nonzero error flag
! was returned by CSHEP2.

! On output:

!       C = Value of C at (PX,PY) unless IER == 1, in
!           which case no values are returned.

!       CX,CY = First partial derivatives of C at (PX,PY)
!               unless IER == 1.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NR, DX, DY or RMAX is invalid.
!             IER = 2 if no errors were encountered but
!                     (PX,PY) is not within the radius R(k)
!                     for any node k (and thus C=CX=CY=0).

! Modules required by CS2GRD:  None

! Intrinsic functions called by CS2GRD:  INT, SQRT

! **********************************************************

      integer :: I, IMAX, IMIN, J, JMAX, JMIN, K, KP
      real :: CK, CKX, CKY, D, DELX, DELY, R, SW,
     &                 SWC, SWCX, SWCY, SWS, SWX, SWY, T, W,
     &                 WX, WY, XP, YP

! Local parameters:

! CK =        Value of cubic nodal function C(K) at P
! CKX,CKY =   Partial derivatives of C(K) with respect to X
!               and Y, respectively
! D =         Distance between P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes
!               in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at P
! SWCX,SWCY = Partial derivatives of SWC with respect to X
!               and Y, respectively
! SWS =       SW**2
! SWX,SWY =   Partial derivatives of SW with respect to X
!               and Y, respectively
! T =         Temporary variable
! W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! WX,WY =     Partial derivatives of W with respect to X
!               and Y, respectively
! XP,YP =     Local copies of PX and PY -- coordinates of P

      XP = PX
      YP = PY
      if (N < 10 .or. NR < 1 .or. DX <= 0.  .or.
     &    DY <= 0. .or. RMAX < 0.) GOTO 6

! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
!   the range of the search for nodes whose radii include
!   P.  The cells which must be searched are those inter-
!   sected by (or contained in) a circle of radius RMAX
!   centered at P.

      IMIN = INT((XP-XMIN-RMAX)/DX) + 1
      IMAX = INT((XP-XMIN+RMAX)/DX) + 1
      if (IMIN < 1) IMIN = 1
      if (IMAX > NR) IMAX = NR
      JMIN = INT((YP-YMIN-RMAX)/DY) + 1
      JMAX = INT((YP-YMIN+RMAX)/DY) + 1
      if (JMIN < 1) JMIN = 1
      if (JMAX > NR) JMAX = NR

! The following is a test for no cells within the circle
!   of radius RMAX.

      if (IMIN > IMAX .or. JMIN > JMAX) GOTO 7

! C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is
!   from K = 1 to N, C(K) is the cubic nodal function value,
!   and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and dist-
!   ance D(K).  Thus

!        CX = (SWCX*SW - SWC*SWX)/SW**2  and
!        CY = (SWCY*SW - SWC*SWY)/SW**2

!   where SWCX and SWX are partial derivatives with respect
!   to X of SWC and SW, respectively.  SWCY and SWY are de-
!   fined similarly.

      SW = 0.
      SWX = 0.
      SWY = 0.
      SWC = 0.
      SWCX = 0.
      SWCY = 0.

! Outer loop on cells (I,J).

      do 4 J = JMIN,JMAX
        do 3 I = IMIN,IMAX
          K = LCELL(I,J)
          if (K == 0) GOTO 3

! Inner loop on nodes K.

    1     DELX = XP - X(K)
          DELY = YP - Y(K)
          D = SQRT(DELX*DELX + DELY*DELY)
          R = RW(K)
          if (D >= R) GOTO 2
          if (D == 0.) GOTO 5
          T = (1.0/D - 1.0/R)
          W = T*T*T
          T = -3.0*T*T/(D*D*D)
          WX = DELX*T
          WY = DELY*T
          T = A(2,K)*DELX + A(3,K)*DELY + A(6,K)
          CKY = ( 3.0*A(4,K)*DELY + A(3,K)*DELX +
     &            2.0*A(7,K) )*DELY + T*DELX + A(9,K)
          T = T*DELY + A(8,K)
          CKX = ( 3.0*A(1,K)*DELX + A(2,K)*DELY +
     &            2.0*A(5,K) )*DELX + T
          CK = ( (A(1,K)*DELX+A(5,K))*DELX + T )*DELX +
     &         ( (A(4,K)*DELY+A(7,K))*DELY + A(9,K) )*DELY +
     &         F(K)
          SW = SW + W
          SWX = SWX + WX
          SWY = SWY + WY
          SWC = SWC + W*CK
          SWCX = SWCX + WX*CK + W*CKX
          SWCY = SWCY + WY*CK + W*CKY

! Bottom of loop on nodes in cell (I,J).

    2     KP = K
          K = LNEXT(KP)
          if (K /= KP) GOTO 1
    3     continue
    4   continue

! SW = 0 iff P is not within the radius R(K) for any node K.

      if (SW == 0.) GOTO 7
      C = SWC/SW
      SWS = SW*SW
      CX = (SWCX*SW - SWC*SWX)/SWS
      CY = (SWCY*SW - SWC*SWY)/SWS
      IER = 0
      return

! (PX,PY) = (X(K),Y(K)).

    5 C = F(K)
      CX = A(8,K)
      CY = A(9,K)
      IER = 0
      return

! Invalid input parameter.

    6 IER = 1
      return

! No cells contain a point within RMAX of P, or
!   SW = 0 and thus D >= RW(K) for all K.

    7 C = 0.
      CX = 0.
      CY = 0.
      IER = 2
      return
      end subroutine do_CS2GRD_sg
      subroutine CS2HES_sg (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,
     &                   YMIN,DX,DY,RMAX,RW,A, C,CX,CY,CXX,
     &                   CXY,CYY,IER)
      integer :: N, NR, LCELL(NR,NR), LNEXT(N), IER
      real :: PX, PY, X(N), Y(N), F(N), XMIN, YMIN,
     &                 DX, DY, RMAX, RW(N), A(9,N), C, CX,
     &                 CY, CXX, CXY, CYY

! **********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   This subroutine computes the value, gradient, and
! Hessian at P = (PX,PY) of the interpolatory function C
! defined in Subroutine CSHEP2.  C is a weighted sum of
! cubic nodal functions.

! On input:

!       PX,PY = Cartesian coordinates of the point P at
!               which C and its partial derivatives are
!               to be evaluated.

!       N = Number of nodes and data values defining C.
!           N >= 10.

!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.

!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2_sg.  NR >= 1.

!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2_sg.

!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2_sg.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2_sg.

!       RMAX = Largest element in RW -- maximum radius R(k).

!       RW = Array of length N containing the the radii R(k)
!            which enter into the weights W(k) defining C.

!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.

!   Input parameters are not altered by this subroutine.
! The parameters other than PX and PY should be input
! unaltered from their values on output from CSHEP2.  This
! subroutine should not be called if a nonzero error flag
! was returned by CSHEP2.

! On output:

!       C = Value of C at (PX,PY) unless IER == 1, in
!           which case no values are returned.

!       CX,CY = First partial derivatives of C at (PX,PY)
!               unless IER == 1.

!       CXX,CXY,CYY = Second partial derivatives of C at
!                     (PX,PY) unless IER == 1.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NR, DX, DY or RMAX is invalid.
!             IER = 2 if no errors were encountered but
!                     (PX,PY) is not within the radius R(k)
!                     for any node k (and thus C = 0).

! Modules required by CS2HES_sg:  None

! Intrinsic functions called by CS2HES_sg:  INT, SQRT

! **********************************************************

      integer :: I, IMAX, IMIN, J, JMAX, JMIN, K, KP
      real :: CK, CKX, CKXX, CKXY, CKY, CKYY, D,
     &                 DELX, DELY, DXSQ, DYSQ, R, SW, SWC,
     &                 SWCX, SWCXX, SWCXY, SWCY, SWCYY, SWS,
     &                 SWX, SWXX, SWXY, SWY, SWYY, T1, T2,
     &                 T3, T4, W, WX, WXX, WXY, WY, WYY, XP,
     &                 YP, D6

! Local parameters:

! CK =        Value of cubic nodal function C(K) at P
! CKX,CKY =   Partial derivatives of C(K) with respect to X
!               and Y, respectively
! CKXX,CKXY,CKYY = Second partial derivatives of CK
! D =         Distance between P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! DXSQ,DYSQ = DELX**2, DELY**2
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes
!               in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at P
! SWCX,SWCY = Partial derivatives of SWC with respect to X
!               and Y, respectively
! SWCXX,SWCXY,SWCYY = Second partial derivatives of SWC
! SWS =       SW**2
! SWX,SWY =   Partial derivatives of SW with respect to X
!               and Y, respectively
! SWXX,SWXY,SWYY = Second partial derivatives of SW
! T1,T2,T3,T4 = Temporary variables
! W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! WX,WY =     Partial derivatives of W with respect to X
!               and Y, respectively
! WXX,WXY,WYY = Second partial derivatives of W
! XP,YP =     Local copies of PX and PY -- coordinates of P

      XP = PX
      YP = PY
      if (N < 10 .or. NR < 1 .or. DX <= 0.  .or.
     &    DY <= 0. .or. RMAX < 0.) GOTO 6

! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
!   the range of the search for nodes whose radii include
!   P.  The cells which must be searched are those inter-
!   sected by (or contained in) a circle of radius RMAX
!   centered at P.

      IMIN = INT((XP-XMIN-RMAX)/DX) + 1
      IMAX = INT((XP-XMIN+RMAX)/DX) + 1
      if (IMIN < 1) IMIN = 1
      if (IMAX > NR) IMAX = NR
      JMIN = INT((YP-YMIN-RMAX)/DY) + 1
      JMAX = INT((YP-YMIN+RMAX)/DY) + 1
      if (JMIN < 1) JMIN = 1
      if (JMAX > NR) JMAX = NR

! The following is a test for no cells within the circle
!   of radius RMAX.

      if (IMIN > IMAX .or. JMIN > JMAX) GOTO 7

! C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is
!   from K = 1 to N, C(K) is the cubic nodal function value,
!   and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and dist-
!   ance D(K).  Thus

!        CX = (SWCX*SW - SWC*SWX)/SW**2  and
!        CY = (SWCY*SW - SWC*SWY)/SW**2

!   where SWCX and SWX are partial derivatives with respect
!   to x of SWC and SW, respectively.  SWCY and SWY are de-
!   fined similarly.  The second partials are

!        CXX = ( SW*(SWCXX -    2*SWX*CX) - SWC*SWXX )/SW**2
!        CXY = ( SW*(SWCXY-SWX*CY-SWY*CX) - SWC*SWXY )/SW**2
!        CYY = ( SW*(SWCYY -    2*SWY*CY) - SWC*SWYY )/SW**2

!   where SWCXX and SWXX are second partials with respect
!   to x, SWCXY and SWXY are mixed partials, and SWCYY and
!   SWYY are second partials with respect to y.

      SW = 0.
      SWX = 0.
      SWY = 0.
      SWXX = 0.
      SWXY = 0.
      SWYY = 0.
      SWC = 0.
      SWCX = 0.
      SWCY = 0.
      SWCXX = 0.
      SWCXY = 0.
      SWCYY = 0.

! Outer loop on cells (I,J).

      do 4 J = JMIN,JMAX
        do 3 I = IMIN,IMAX
          K = LCELL(I,J)
          if (K == 0) GOTO 3

! Inner loop on nodes K.

    1     DELX = XP - X(K)
          DELY = YP - Y(K)
          DXSQ = DELX*DELX
          DYSQ = DELY*DELY
          D = SQRT(DXSQ + DYSQ)
          R = RW(K)
          if (D >= R) GOTO 2
          if (D == 0.) GOTO 5
          T1 = (1.0/D - 1.0/R)
          W = T1*T1*T1
          T2 = -3.0*T1*T1/(D*D*D)
          WX = DELX*T2
          WY = DELY*T2
          D6 = D*D*D*D*D*D
          T1 = 3.0*T1*(2.0+3.0*D*T1)/D6
          WXX = T1*DXSQ + T2
          WXY = T1*DELX*DELY
          WYY = T1*DYSQ + T2
          T1 = A(1,K)*DELX + A(2,K)*DELY + A(5,K)
          T2 = T1 + T1 + A(1,K)*DELX
          T3 = A(4,K)*DELY + A(3,K)*DELX + A(7,K)
          T4 = T3 + T3 + A(4,K)*DELY
          CK = (T1*DELX + A(6,K)*DELY + A(8,K))*DELX +
     &         (T3*DELY + A(9,K))*DELY + F(K)
          CKX = T2*DELX + (A(3,K)*DELY+A(6,K))*DELY + A(8,K)
          CKY = T4*DELY + (A(2,K)*DELX+A(6,K))*DELX + A(9,K)
          CKXX = T2 + 3.0*A(1,K)*DELX
          CKXY = 2.0*(A(2,K)*DELX + A(3,K)*DELY) + A(6,K)
          CKYY = T4 + 3.0*A(4,K)*DELY
          SW = SW + W
          SWX = SWX + WX
          SWY = SWY + WY
          SWXX = SWXX + WXX
          SWXY = SWXY + WXY
          SWYY = SWYY + WYY
          SWC = SWC + W*CK
          SWCX = SWCX + WX*CK + W*CKX
          SWCY = SWCY + WY*CK + W*CKY
          SWCXX = SWCXX + W*CKXX + 2.0*WX*CKX + CK*WXX
          SWCXY = SWCXY + W*CKXY + WX*CKY + WY*CKX + CK*WXY
          SWCYY = SWCYY + W*CKYY + 2.0*WY*CKY + CK*WYY

! Bottom of loop on nodes in cell (I,J).

    2     KP = K
          K = LNEXT(KP)
          if (K /= KP) GOTO 1
    3     continue
    4   continue

! SW = 0 iff P is not within the radius R(K) for any node K.

      if (SW == 0.) GOTO 7
      C = SWC/SW
      SWS = SW*SW
      CX = (SWCX*SW - SWC*SWX)/SWS
      CY = (SWCY*SW - SWC*SWY)/SWS
      CXX = (SW*(SWCXX-2.0*SWX*CX) - SWC*SWXX)/SWS
      CXY = (SW*(SWCXY-SWY*CX-SWX*CY) - SWC*SWXY)/SWS
      CYY = (SW*(SWCYY-2.0*SWY*CY) - SWC*SWYY)/SWS
      IER = 0
      return

! (PX,PY) = (X(K),Y(K)).

    5 C = F(K)
      CX = A(8,K)
      CY = A(9,K)
      CXX = 2.0*A(5,K)
      CXY = A(6,K)
      CYY = 2.0*A(7,K)
      IER = 0
      return

! Invalid input parameter.

    6 IER = 1
      return

! No cells contain a point within RMAX of P, or
!   SW = 0 and thus D >= RW(K) for all K.

    7 C = 0.
      CX = 0.
      CY = 0.
      CXX = 0.
      CXY = 0.
      CYY = 0.
      IER = 2
      return
      end subroutine CS2HES_sg
      subroutine GETNP2_sg (PX,PY,X,Y,NR,LCELL,LNEXT,XMIN,YMIN,
     &                   DX,DY, NP,DSQ)
      integer :: NR, LCELL(NR,NR), LNEXT(*), NP
      real :: PX, PY, X(*), Y(*), XMIN, YMIN, DX, DY, DSQ

! **********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   Given a set of N nodes and the data structure defined in
! Subroutine STORE2_sg, this subroutine uses the cell method to
! find the closest unmarked node NP to a specified point P.
! NP is then marked by setting LNEXT(NP) to -LNEXT(NP).  (A
! node is marked if and only if the corresponding LNEXT element
! is negative.  The absolute values of LNEXT elements,
! however, must be preserved.)  Thus, the closest M nodes to
! P may be determined by a sequence of M calls to this rou-
! tine.  Note that if the nearest neighbor to node K is to
! be determined (PX = X(K) and PY = Y(K)), then K should be
! marked before the call to this routine.

!   The search is begun in the cell containing (or closest
! to) P and proceeds outward in rectangular layers until all
! cells which contain points within distance R of P have
! been searched, where R is the distance from P to the first
! unmarked node encountered (infinite if no unmarked nodes
! are present).

!   This code is essentially unaltered from the subroutine
! of the same name in QSHEP2D.

! On input:

!       PX,PY = Cartesian coordinates of the point P whose
!               nearest unmarked neighbor is to be found.

!       X,Y = Arrays of length N, for N >= 2, containing
!             the Cartesian coordinates of the nodes.

!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2_sg.  NR >= 1.

!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2_sg.

!       LNEXT = Array of length N containing next-node
!               indexes (or their negatives).  Refer to
!               Subroutine STORE2_sg.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2_sg.

!   Input parameters other than LNEXT are not altered by
! this routine.  With the exception of (PX,PY) and the signs
! of LNEXT elements, these parameters should be unaltered
! from their values on output from Subroutine STORE2_sg.

! On output:

!       NP = Index (for X and Y) of the nearest unmarked
!            node to P, or 0 if all nodes are marked or NR
!            < 1 or DX <= 0 or DY <= 0.  LNEXT(NP)
!            < 0 if NP /= 0.

!       DSQ = Squared Euclidean distance between P and node
!             NP, or 0 if NP = 0.

! Modules required by GETNP2_sg:  None

! Intrinsic functions called by GETNP2_sg:  ABS, INT, SQRT

! **********************************************************

      integer :: I, I0, I1, I2, IMAX, IMIN, J, J0, J1, J2,
     &        JMAX, JMIN, L, LMIN, LN
      LOGICAL :: FIRST
      real :: DELX, DELY, R, RSMIN, RSQ, XP, YP

! Local parameters:

! DELX,DELY =   PX-XMIN, PY-YMIN
! FIRST =       Logical variable with value TRUE iff the
!                 first unmarked node has yet to be
!                 encountered
! I,J =         Cell indexes in the range [I1,I2] X [J1,J2]
! I0,J0 =       Indexes of the cell containing or closest
!                 to P
! I1,I2,J1,J2 = Range of cell indexes defining the layer
!                 whose intersection with the range
!                 [IMIN,IMAX] X [JMIN,JMAX] is currently
!                 being searched
! IMIN,IMAX =   Cell row indexes defining the range of the
!                 search
! JMIN,JMAX =   Cell column indexes defining the range of
!                 the search
! L,LN =        Indexes of nodes in cell (I,J)
! LMIN =        Current candidate for NP
! R =           Distance from P to node LMIN
! RSMIN =       Squared distance from P to node LMIN
! RSQ =         Squared distance from P to node L
! XP,YP =       Local copy of PX,PY -- coordinates of P

      XP = PX
      YP = PY

! Test for invalid input parameters.

      if (NR < 1 .or. DX <= 0. .or. DY <= 0.) GOTO 9

! Initialize parameters.

      FIRST = .true.
      IMIN = 1
      IMAX = NR
      JMIN = 1
      JMAX = NR
      DELX = XP - XMIN
      DELY = YP - YMIN
      I0 = INT(DELX/DX) + 1
      if (I0 < 1) I0 = 1
      if (I0 > NR) I0 = NR
      J0 = INT(DELY/DY) + 1
      if (J0 < 1) J0 = 1
      if (J0 > NR) J0 = NR
      I1 = I0
      I2 = I0
      J1 = J0
      J2 = J0

! Outer loop on layers, inner loop on layer cells, excluding
!   those outside the range [IMIN,IMAX] X [JMIN,JMAX].

    1 do 6 J = J1,J2
        if (J > JMAX) GOTO 7
        if (J < JMIN) GOTO 6
        do 5 I = I1,I2
          if (I > IMAX) GOTO 6
          if (I < IMIN) GOTO 5
          if (J /= J1 .and. J /= J2 .and. I /= I1
     &       .and. I /= I2) GOTO 5

! Search cell (I,J) for unmarked nodes L.

          L = LCELL(I,J)
          if (L == 0) GOTO 5

!   Loop on nodes in cell (I,J).

    2     LN = LNEXT(L)
          if (LN < 0) GOTO 4

!   Node L is not marked.

          RSQ = (X(L)-XP)**2 + (Y(L)-YP)**2
          if (.not. FIRST) GOTO 3

!   Node L is the first unmarked neighbor of P encountered.
!     Initialize LMIN to the current candidate for NP, and
!     RSMIN to the squared distance from P to LMIN.  IMIN,
!     IMAX, JMIN, and JMAX are updated to define the smal-
!     lest rectangle containing a circle of radius R =
!     Sqrt(RSMIN) centered at P, and contained in [1,NR] X
!     [1,NR] (except that, if P is outside the rectangle
!     defined by the nodes, it is possible that IMIN > NR,
!     IMAX < 1, JMIN > NR, or JMAX < 1).  FIRST is reset to
!     FALSE.

          LMIN = L
          RSMIN = RSQ
          R = SQRT(RSMIN)
          IMIN = INT((DELX-R)/DX) + 1
          if (IMIN < 1) IMIN = 1
          IMAX = INT((DELX+R)/DX) + 1
          if (IMAX > NR) IMAX = NR
          JMIN = INT((DELY-R)/DY) + 1
          if (JMIN < 1) JMIN = 1
          JMAX = INT((DELY+R)/DY) + 1
          if (JMAX > NR) JMAX = NR
          FIRST = .false.
          GOTO 4

!   Test for node L closer than LMIN to P.

    3     if (RSQ >= RSMIN) GOTO 4

!   Update LMIN and RSMIN.

          LMIN = L
          RSMIN = RSQ

!   Test for termination of loop on nodes in cell (I,J).

    4     if (ABS(LN) == L) GOTO 5
          L = ABS(LN)
          GOTO 2
    5     continue
    6   continue

! Test for termination of loop on cell layers.

    7 if (I1 <= IMIN .and. I2 >= IMAX  .and.
     &    J1 <= JMIN .and. J2 >= JMAX) GOTO 8
      I1 = I1 - 1
      I2 = I2 + 1
      J1 = J1 - 1
      J2 = J2 + 1
      GOTO 1

! Unless no unmarked nodes were encountered, LMIN is the
!   closest unmarked node to P.

    8 if (FIRST) GOTO 9
      NP = LMIN
      DSQ = RSMIN
      LNEXT(LMIN) = -LNEXT(LMIN)
      return

! Error:  NR, DX, or DY is invalid or all nodes are marked.

    9 NP = 0
      DSQ = 0.
      return
      end subroutine GETNP2_sg

      subroutine GIVENS_sg (A,B, C,S)
      real :: A, B, C, S

! **********************************************************

!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!   This subroutine constructs the Givens plane rotation,

!           ( C  S)
!       G = (     ) , where C*C + S*S = 1,
!           (-S  C)

! which zeros the second component of the vector (A,B)**T
! (transposed).  Subroutine ROTATE_sg may be called to apply
! the transformation to a 2 by N matrix.

!   This routine is identical to subroutine SROTG from the
! LINPACK BLAS (Basic Linear Algebra Subroutines).

! On input:

!       A,B = Components of the vector defining the rota-
!             tion.  These are overwritten by values R
!             and Z (described below) which define C and S.

! On output:

!       A = Signed Euclidean norm R of the input vector:
!           R = +/-SQRT(A*A + B*B)

!       B = Value Z such that:
!             C = SQRT(1-Z*Z) and S=Z if ABS(Z) <= 1, and
!             C = 1/Z and S = SQRT(1-C*C) if ABS(Z) > 1.

!       C = +/-(A/R) or 1 if R = 0.

!       S = +/-(B/R) or 0 if R = 0.

! Modules required by GIVENS_sg:  None

! Intrinsic functions called by GIVENS_sg:  ABS, SQRT

! **********************************************************

      real :: AA, BB, R, U, V

! Local parameters:

! AA,BB = Local copies of A and B
! R =     C*A + S*B = +/-SQRT(A*A+B*B)
! U,V =   Variables used to scale A and B for computing R

      AA = A
      BB = B
      if (ABS(AA) <= ABS(BB)) GOTO 1

! ABS(A) > ABS(B).

      U = AA + AA
      V = BB/U
      R = SQRT(.25 + V*V) * U
      C = AA/R
      S = V * (C + C)

! Note that R has the sign of A, C > 0, and S has
!   SIGN(A)*SIGN(B).

      B = S
      A = R
      return

! ABS(A) <= ABS(B).

    1 if (BB == 0.) GOTO 2
      U = BB + BB
      V = AA/U

! Store R in A.

      A = SQRT(.25 + V*V) * U
      S = BB/A
      C = V * (S + S)

! Note that R has the sign of B, S > 0, and C has
!   SIGN(A)*SIGN(B).

      B = 1.
      if (C /= 0.) B = 1./C
      return

! A = B = 0.

    2 C = 1.
      S = 0.
      return
      end subroutine GIVENS_sg

      subroutine ROTATE_sg (N,C,S, X,Y )
      integer :: N
      real :: C, S, X(N), Y(N)

! **********************************************************

!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!                                                ( C  S)
!   This subroutine applies the Givens rotation  (     )  to
!                                                (-S  C)
!                    (X(1) ... X(N))
! the 2 by N matrix  (             ) .
!                    (Y(1) ... Y(N))

!   This routine is identical to subroutine SROT from the
! LINPACK BLAS (Basic Linear Algebra Subroutines).

! On input:

!       N = Number of columns to be rotated.

!       C,S = Elements of the Givens rotation.  Refer to
!             subroutine GIVENS_sg.

! The above parameters are not altered by this routine.

!       X,Y = Arrays of length >= N containing the compo-
!             nents of the vectors to be rotated.

! On output:

!       X,Y = Arrays containing the rotated vectors (not
!             altered if N < 1).

! Modules required by ROTATE_sg:  None

! **********************************************************

      integer :: I
      real :: XI, YI

      do I = 1,N
        XI = X(I)
        YI = Y(I)
        X(I) = C*XI + S*YI
        Y(I) = -S*XI + C*YI
      end do
      return
      end subroutine ROTATE_sg
      subroutine SETUP2_sg (XK,YK,ZK,XI,YI,ZI,S1,S2,S3,R, ROW)
      real :: XK, YK, ZK, XI, YI, ZI, S1, S2, S3,
     &                 R, ROW(10)

! **********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   This subroutine sets up the I-th row of an augmented re-
! gression matrix for a weighted least squares fit of a
! cubic function f(x,y) to a set of data values z, where
! f(XK,YK) = ZK.  The first four columns (cubic terms) are
! scaled by S3, the next three columns (quadratic terms)
! are scaled by S2, and the eighth and ninth columns (lin-
! ear terms) are scaled by S1.

! On input:

!       XK,YK = Coordinates of node K.

!       ZK = Data value at node K to be interpolated by f.

!       XI,YI,ZI = Coordinates and data value at node I.

!       S1,S2,S3 = Scale factors.

!       R = Radius of influence about node K defining the
!           weight.

! The above parameters are not altered by this routine.

!       ROW = Array of length 10.

! On output:

!       ROW = Array containing a row of the augmented re-
!             gression matrix.

! Modules required by SETUP2_sg:  None

! Intrinsic function called by SETUP2_sg:  SQRT

! **********************************************************

      integer :: I
      real :: D, DX, DXSQ, DY, DYSQ, W, W1, W2, W3

! Local parameters:

! D =    Distance between nodes K and I
! DX =   XI - XK
! DXSQ = DX*DX
! DY =   YI - YK
! DYSQ = DY*DY
! I =    DO-loop index
! W =    Weight associated with the row:  (R-D)/(R*D)
!          (0 if D = 0 or D > R)
! W1 =   S1*W
! W2 =   S2*W
! W3 =   W3*W

      DX = XI - XK
      DY = YI - YK
      DXSQ = DX*DX
      DYSQ = DY*DY
      D = SQRT(DXSQ + DYSQ)
      if (D <= 0. .or. D >= R) GOTO 1
      W = (R-D)/R/D
      W1 = S1*W
      W2 = S2*W
      W3 = S3*W
      ROW(1) = DXSQ*DX*W3
      ROW(2) = DXSQ*DY*W3
      ROW(3) = DX*DYSQ*W3
      ROW(4) = DYSQ*DY*W3
      ROW(5) = DXSQ*W2
      ROW(6) = DX*DY*W2
      ROW(7) = DYSQ*W2
      ROW(8) = DX*W1
      ROW(9) = DY*W1
      ROW(10) = (ZI - ZK)*W
      return

! Nodes K and I coincide or node I is outside of the radius
!   of influence.  Set ROW to the zero vector.

    1 do I = 1,10
        ROW(I) = 0.
      end do
      return
      end subroutine SETUP2_sg

      subroutine STORE2_sg (N,X,Y,NR, LCELL,LNEXT,XMIN,YMIN,DX,
     &                   DY,IER)
      integer :: N, NR, LCELL(NR,NR), LNEXT(N), IER
      real :: X(N), Y(N), XMIN, YMIN, DX, DY

! **********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   03/28/97

!   Given a set of N arbitrarily distributed nodes in the
! plane, this subroutine creates a data structure for a
! cell-based method of solving closest-point problems.  The
! smallest rectangle containing the nodes is partitioned
! into an NR by NR uniform grid of cells, and nodes are as-
! sociated with cells.  In particular, the data structure
! stores the indexes of the nodes contained in each cell.
! For a uniform random distribution of nodes, the nearest
! node to an arbitrary point can be determined in constant
! expected time.

!   This code is essentially unaltered from the subroutine
! of the same name in QSHEP2D.

! On input:

!       N = Number of nodes.  N >= 2.

!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.

!       NR = Number of rows and columns in the grid.  The
!            cell density (average number of nodes per cell)
!            is D = N/(NR**2).  A recommended value, based
!            on empirical evidence, is D = 3 -- NR =
!            Sqrt(N/3).  NR >= 1.

! The above parameters are not altered by this routine.

!       LCELL = Array of length >= NR**2.

!       LNEXT = Array of length >= N.

! On output:

!       LCELL = NR by NR cell array such that LCELL(I,J)
!               contains the index (for X and Y) of the
!               first node (node with smallest index) in
!               cell (I,J), or LCELL(I,J) = 0 if no nodes
!               are contained in the cell.  The upper right
!               corner of cell (I,J) has coordinates (XMIN+
!               I*DX,YMIN+J*DY).  LCELL is not defined if
!               IER /= 0.

!       LNEXT = Array of next-node indexes such that
!               LNEXT(K) contains the index of the next node
!               in the cell which contains node K, or
!               LNEXT(K) = K if K is the last node in the
!               cell for K = 1,...,N.  (The nodes contained
!               in a cell are ordered by their indexes.)
!               If, for example, cell (I,J) contains nodes
!               2, 3, and 5 (and no others), then LCELL(I,J)
!               = 2, LNEXT(2) = 3, LNEXT(3) = 5, and
!               LNEXT(5) = 5.  LNEXT is not defined if
!               IER /= 0.

!       XMIN,YMIN = Cartesian coordinates of the lower left
!                   corner of the rectangle defined by the
!                   nodes (smallest nodal coordinates) un-
!                   less IER = 1.  The upper right corner is
!                   (XMAX,YMAX) for XMAX = XMIN + NR*DX and
!                   YMAX = YMIN + NR*DY.

!       DX,DY = Dimensions of the cells unless IER = 1.  DX
!               = (XMAX-XMIN)/NR and DY = (YMAX-YMIN)/NR,
!               where XMIN, XMAX, YMIN, and YMAX are the
!               extrema of X and Y.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N < 2 or NR < 1.
!             IER = 2 if DX = 0 or DY = 0.

! Modules required by STORE2_sg:  None

! Intrinsic functions called by STORE2_sg:  DBLE, INT

! **********************************************************

      integer :: I, J, K, L, NN, NNR
      real :: DELX, DELY, XMN, XMX, YMN, YMX

! Local parameters:

! DELX,DELY = Components of the cell dimensions -- local
!               copies of DX,DY
! I,J =       Cell indexes
! K =         Nodal index
! L =         Index of a node in cell (I,J)
! NN =        Local copy of N
! NNR =       Local copy of NR
! XMN,XMX =   Range of nodal X coordinates
! YMN,YMX =   Range of nodal Y coordinates

      NN = N
      NNR = NR
      if (NN < 2 .or. NNR < 1) GOTO 5

! Compute the dimensions of the rectangle containing the
!   nodes.

      XMN = X(1)
      XMX = XMN
      YMN = Y(1)
      YMX = YMN
      do K = 2,NN
        if (X(K) < XMN) XMN = X(K)
        if (X(K) > XMX) XMX = X(K)
        if (Y(K) < YMN) YMN = Y(K)
        if (Y(K) > YMX) YMX = Y(K)
      end do
      XMIN = XMN
      YMIN = YMN

! Compute cell dimensions and test for zero area.

      DELX = (XMX-XMN)/DBLE(NNR)
      DELY = (YMX-YMN)/DBLE(NNR)
      DX = DELX
      DY = DELY
      if (DELX == 0. .or. DELY == 0.) GOTO 6

! Initialize LCELL.

      do J = 1,NNR
        do I = 1,NNR
          LCELL(I,J) = 0
        end do
      end do

! Loop on nodes, storing indexes in LCELL and LNEXT.

      do K = NN,1,-1
        I = INT((X(K)-XMN)/DELX) + 1
        if (I > NNR) I = NNR
        J = INT((Y(K)-YMN)/DELY) + 1
        if (J > NNR) J = NNR
        L = LCELL(I,J)
        LNEXT(K) = L
        if (L == 0) LNEXT(K) = K
        LCELL(I,J) = K
      end do

! No errors encountered.

      IER = 0
      return

! Invalid input parameter.

    5 IER = 1
      return

! DX = 0 or DY = 0.

    6 IER = 2
      return
      end subroutine STORE2_sg
