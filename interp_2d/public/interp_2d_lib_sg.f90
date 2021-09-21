! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module interp_2d_lib_sg
      ! single precision library for 2D interpolation

      implicit none

      contains ! the procedure interface for the library
      ! client programs should only call these routines.

! Contents

! Rectangular-grid of data points

   ! interp_RGBI3P_sg -- point interpolation (Akima)
   ! interp_RGSF3P_sg -- surface interpolation (Akima)
   
   ! interp_mkbicub_sg -- bicubic splines
      
! Scattered set of data points

   ! interp_CS2VAL_sg -- point interpolation (Renka)
   ! interp_CS2GRD_sg -- point interpolation with gradients (Renka)



! see the documents in refs directory for more info.

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************
! ***********************************************************************



* Rectangular-grid bivariate interpolation and surface fitting
* from ACM Algorithm 760., ACM Trans. Math. Software (22) 1996, 357-361.
* Hiroshi Akima
* U.S. Department of Commerce, NTIA/ITS
* Version of 1995/08

* NOTE: versions for scattered data follow.

* see interp_2d_lib_db for double precision version.

      subroutine interp_RGBI3P_sg(MD,NXD,NYD,XD,YD,ZD,NIP,XI,YI,ZI,IER,WK)
         integer, intent(in) :: MD, NXD, NYD, NIP
         real, intent(in) :: XD(NXD), YD(NYD), ZD(NXD,NYD), XI(NIP), YI(NIP)
         real, intent(inout) :: ZI(NIP), WK(3,NXD,NYD)
         integer, intent(out) :: IER
*
* This subroutine performs interpolation of a bivariate function,
* z(x,y), on a rectangular grid in the x-y plane.  It is based on
* the revised Akima method.
*
* In this subroutine, the interpolating function is a piecewise
* function composed of a set of bicubic (bivariate third-degree)
* polynomials, each applicable to a rectangle of the input grid
* in the x-y plane.  Each polynomial is determined locally.
*
* This subroutine has the accuracy of a bicubic polynomial, i.e.,
* it interpolates accurately when all data points lie on a
* surface of a bicubic polynomial.
*
* The grid lines can be unevenly spaced.
*
* The input arguments are
*   MD  = mode of computation
*       = 1 for new XD, YD, or ZD data (default)
*       = 2 for old XD, YD, and ZD data,
*   NXD = number of the input-grid data points in the x
*         coordinate (must be 2 or greater),
*   NYD = number of the input-grid data points in the y
*         coordinate (must be 2 or greater),
*   XD  = array of dimension NXD containing the x coordinates
*         of the input-grid data points (must be in a
*         monotonic increasing order),
*   YD  = array of dimension NYD containing the y coordinates
*         of the input-grid data points (must be in a
*         monotonic increasing order),
*   ZD  = two-dimensional array of dimension NXD*NYD
*         containing the z(x,y) values at the input-grid data
*         points,
*   NIP = number of the output points at which interpolation
*         of the z value is desired (must be 1 or greater),
*   XI  = array of dimension NIP containing the x coordinates
*         of the output points,
*   YI  = array of dimension NIP containing the y coordinates
*         of the output points.
*
* The output arguments are
*   ZI  = array of dimension NIP where the interpolated z
*         values at the output points are to be stored,
*   IER = error flag
*       = 0 for no errors
*       = 1 for NXD = 1 or less
*       = 2 for NYD = 1 or less
*       = 3 for identical XD values or
*               XD values out of sequence
*       = 4 for identical YD values or
*               YD values out of sequence
*       = 5 for NIP = 0 or less.
*
* The other argument is
*   WK  = three dimensional array of dimension 3*NXD*NYD used
*         internally as a work area.
*
* The very first call to this subroutine and the call with a new
* XD, YD, and ZD array must be made with MD=1.  The call with MD=2
* must be preceded by another call with the same XD, YD, and ZD
* arrays.  Between the call with MD=2 and its preceding call, the
* WK array must not be disturbed.
*
         call do_RGBI3P_sg(MD,NXD,NYD,XD,YD,ZD,NIP,XI,YI,ZI,IER,WK)
      end subroutine interp_RGBI3P_sg
      
      
      subroutine interp_RGSF3P_sg(MD,NXD,NYD,XD,YD,ZD,NXI,XI,NYI,YI,ZI,IER,WK)
         integer, intent(in) :: MD, NXD, NYD, NXI, NYI
         real, intent(in) :: XD(NXD), YD(NYD), ZD(NXD,NYD), XI(NXI), YI(NYI)
         real, intent(inout) :: ZI(NXI,NYI), WK(3,NXD,NYD)
         integer, intent(out) :: IER
*
* Rectangular-grid surface fitting
* (a master subroutine of the RGBI3P/RGSF3P_sg subroutine package)
*
* Hiroshi Akima
* U.S. Department of Commerce, NTIA/ITS
* Version of 1995/08
*
* This subroutine performs surface fitting by interpolating
* values of a bivariate function, z(x,y), on a rectangular grid
* in the x-y plane.  It is based on the revised Akima method.
*
* In this subroutine, the interpolating function is a piecewise
* function composed of a set of bicubic (bivariate third-degree)
* polynomials, each applicable to a rectangle of the input grid
* in the x-y plane.  Each polynomial is determined locally.
*
* This subroutine has the accuracy of a bicubic polynomial, i.e.,
* it fits the surface accurately when all data points lie on a
* surface of a bicubic polynomial.
*
* The grid lines of the input and output data can be unevenly
* spaced.
*
* The input arguments are
*   MD  = mode of computation
*       = 1 for new XD, YD, or ZD data (default)
*       = 2 for old XD, YD, and ZD data,
*   NXD = number of the input-grid data points in the x
*         coordinate (must be 2 or greater),
*   NYD = number of the input-grid data points in the y
*         coordinate (must be 2 or greater),
*   XD  = array of dimension NXD containing the x coordinates
*         of the input-grid data points (must be in a
*         monotonic increasing order),
*   YD  = array of dimension NYD containing the y coordinates
*         of the input-grid data points (must be in a
*         monotonic increasing order),
*   ZD  = two-dimensional array of dimension NXD*NYD
*         containing the z(x,y) values at the input-grid data
*         points,
*   NXI = number of output grid points in the x coordinate
*         (must be 1 or greater),
*   XI  = array of dimension NXI containing the x coordinates
*         of the output grid points,
*   NYI = number of output grid points in the y coordinate
*         (must be 1 or greater),
*   YI  = array of dimension NYI containing the y coordinates
*         of the output grid points.
*
* The output arguments are
*   ZI  = two-dimensional array of dimension NXI*NYI, where
*         the interpolated z values at the output grid
*         points are to be stored,
*   IER = error flag
*       = 0 for no error
*       = 1 for NXD = 1 or less
*       = 2 for NYD = 1 or less
*       = 3 for identical XD values or
*               XD values out of sequence
*       = 4 for identical YD values or
*               YD values out of sequence
*       = 5 for NXI = 0 or less
*       = 6 for NYI = 0 or less.
*
* The other argument is
*   WK  = three-dimensional array of dimension 3*NXD*NYD used
*         internally as a work area.
*
* The very first call to this subroutine and the call with a new
* XD, YD, or ZD array must be made with MD=1.  The call with MD=2
* must be preceded by another call with the same XD, YD, and ZD
* arrays.  Between the call with MD=2 and its preceding call, the
* WK array must not be disturbed.

         call do_RGSF3P_sg(MD,NXD,NYD,XD,YD,ZD,NXI,XI,NYI,YI,ZI,IER,WK)
      end subroutine interp_RGSF3P_sg
      

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************
! ***********************************************************************



! cubic Shepard method for bivariate interpolation of scattered data.
! from ACM Algorithm 790., ACM Trans. Math. Software (25) 1999, 70-73.
! Robert J. Renka
! Dept. of Computer Science
! Univ. of North Texas
! renka@cs.unt.edu

! use CSHEP2_sg to set up the interpolation information.
! use CS2VAL_sg to evaluate it.
! use CS2GRD_sg to get value and derivatives.
! see interp_2d_lib_db for double precision versions.
      
      subroutine interp_CSHEP2_sg(N,X,Y,F,NC,NW,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,IER)
         integer, intent(in) :: N, NC, NW, NR
         integer, intent(out) :: LCELL(NR,NR), LNEXT(N), IER
         real, intent(in) :: X(N), Y(N), F(N)
         real, intent(inout) :: XMIN, YMIN, DX, DY, RMAX, RW(N), A(9,N)
!
! Here's the [slightly edited] documentation.
!
!   This subroutine computes a set of parameters defining a
! C2 (twice continuously differentiable) bivariate function
! C(X,Y) which interpolates data values F at a set of N
! arbitrarily distributed points (X,Y) in the plane (nodes).
! The interpolant C may be evaluated at an arbitrary point
! by function CS2VAL_sg, and its first partial derivatives are
! computed by Subroutine CS2GRD_sg.
!
!   The interpolation scheme is a modified Cubic Shepard
! method:
!
! C = [W(1)*C(1)+W(2)*C(2)+..+W(N)*C(N)]/[W(1)+W(2)+..+W(N)]
!
! for bivariate functions W(k) and C(k).  The nodal func-
! tions are given by
!
!  C(k)(x,y) = A(1,k)*(x-X(k))**3 +
!              A(2,k)*(x-X(k))**2*(y-Y(k)) +
!              A(3,k)*(x-X(k))*(y-Y(k))**2 +
!              A(4,k)*(y-Y(k))**3 + A(5,k)*(x-X(k))**2 +
!              A(6,k)*(x-X(k))*(y-Y(k)) + A(7,k)*(y-Y(k))**2
!              + A(8,k)*(x-X(k)) + A(9,k)*(y-Y(k)) + F(k) .
!
! Thus, C(k) is a cubic function which interpolates the data
! value at node k.  Its coefficients A(,k) are obtained by a
! weighted least squares fit to the closest NC data points
! with weights similar to W(k).  Note that the radius of
! influence for the least squares fit is fixed for each k,
! but varies with k.
!
! The weights are taken to be
!
!   W(k)(x,y) = ( (R(k)-D(k))+ / R(k)*D(k) )**3 ,
!
! where (R(k)-D(k))+ = 0 if R(k) < D(k), and D(k)(x,y) is
! the Euclidean distance between (x,y) and (X(k),Y(k)).  The
! radius of influence R(k) varies with k and is chosen so
! that NW nodes are within the radius.  Note that W(k) is
! not defined at node (X(k),Y(k)), but C(x,y) has limit F(k)
! as (x,y) approaches (X(k),Y(k)).
! On input:
!
!       N = Number of nodes and data values.  N .GE. 10.
!
!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.
!
!       F = Array of length N containing the data values
!           in one-to-one correspondence with the nodes.
!
!       NC = Number of data points to be used in the least
!            squares fit for coefficients defining the nodal
!            functions C(k).  Values found to be optimal for
!            test data sets ranged from 11 to 25.  A recom-
!            mended value for general data sets is NC = 17.
!            For nodes lying on (or close to) a rectangular
!            grid, the recommended value is NC = 11.  In any
!            case, NC must be in the range 9 to Min(40,N-1).
!
!       NW = Number of nodes within (and defining) the radii
!            of influence R(k) which enter into the weights
!            W(k).  For N sufficiently large, a recommended
!            value is NW = 30.  In general, NW should be
!            about 1.5*NC.  1 .LE. NW .LE. Min(40,N-1).
!
!       NR = Number of rows and columns in the cell grid de-
!            fined in Subroutine STORE2.  A rectangle con-
!            taining the nodes is partitioned into cells in
!            order to increase search efficiency.  NR =
!            Sqrt(N/3) is recommended.  NR .GE. 1.
!
! The above parameters are not altered by this routine.
!
!       LCELL = Array of length .GE. NR**2.
!
!       LNEXT = Array of length .GE. N.
!
!       RW = Array of length .GE. N.
!
!       A = Array of length .GE. 9N.
!
! On output:
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  Refer to Subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW -- maximum radius R(k).
!
!       RW = Array containing the the radii R(k) which enter
!            into the weights W(k).
!
!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.
!
!   Note that the output parameters described above are not
! defined unless IER = 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NC, NW, or NR is outside its
!                     valid range.
!             IER = 2 if duplicate nodes were encountered.
!             IER = 3 if all nodes are collinear.
   
         call do_CSHEP2_sg(N,X,Y,F,NC,NW,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,IER)

      end subroutine interp_CSHEP2_sg


      real function interp_CS2VAL_sg(PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,IER)
         integer, intent(in) :: N, NR, LCELL(NR,NR), LNEXT(N)
         real, intent(in) :: PX, PY, X(N), Y(N), F(N), XMIN, YMIN, DX, DY, RMAX, RW(N), A(9,N)
         integer, intent(out) :: IER
!
!   This function returns the value C(PX,PY), where C is the
! weighted sum of cubic nodal functions defined in Subrou-
! tine CSHEP2_sg.  CS2GRD_sg may be called to compute a gradient
! of C along with the value, and/or to test for errors.
! CS2HES may be called to compute a value, first partial
! derivatives, and second partial derivatives at a point.
!
! On input:
!
!       PX,PY = Cartesian coordinates of the point P at
!               which C is to be evaluated.
!
!       N = Number of nodes and data values defining C.
!           N .GE. 10.
!
!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.
!
!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR .GE. 1.
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW -- maximum radius R(k).
!
!       RW = Array containing the the radii R(k) which enter
!            into the weights W(k) defining C.
!
!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.
!
!   Input parameters are not altered by this function.  The
! parameters other than PX and PY should be input unaltered
! from their values on output from CSHEP2.  This function
! should not be called if a nonzero error flag was returned
! by CSHEP2_sg.
!
! On output:
!
!       CS2VAL_sg = Function value C(PX,PY) unless N, NR, DX,
!                DY, or RMAX is invalid, in which case IER /= 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!
         real :: do_CS2VAL_sg         
         interp_CS2VAL_sg = do_CS2VAL_sg(PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,IER)   
      end function interp_CS2VAL_sg
      

      subroutine interp_CS2GRD_sg(PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,C,CX,CY,IER)
         integer, intent(in) :: N, NR, LCELL(NR,NR), LNEXT(N)
         real, intent(in) :: PX, PY, X(N), Y(N), F(N), XMIN, YMIN, DX, DY, RMAX, RW(N), A(9,N)
         real, intent(out) :: C, CX, CY
         integer, intent(out) :: IER
!
!   This subroutine computes the value and gradient at P =
! (PX,PY) of the interpolatory function C defined in Sub-
! routine CSHEP2_sg.  C is a weighted sum of cubic nodal
! functions.
!
! On input:
!
!       PX,PY = Cartesian coordinates of the point P at
!               which C and its partial derivatives are
!               to be evaluated.
!
!       N = Number of nodes and data values defining C.
!           N .GE. 10.
!
!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.
!
!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR .GE. 1.
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW -- maximum radius R(k).
!
!       RW = Array of length N containing the the radii R(k)
!            which enter into the weights W(k) defining C.
!
!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.
!
!   Input parameters are not altered by this subroutine.
! The parameters other than PX and PY should be input
! unaltered from their values on output from CSHEP2_sg.  This
! subroutine should not be called if a nonzero error flag
! was returned by CSHEP2_sg.
!
! On output:
!
!       C = Value of C at (PX,PY) unless IER .EQ. 1, in
!           which case no values are returned.
!
!       CX,CY = First partial derivatives of C at (PX,PY)
!               unless IER .EQ. 1.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NR, DX, DY or RMAX is invalid.
!             IER = 2 if no errors were encountered but
!                     (PX,PY) is not within the radius R(k)
!                     for any node k (and thus C=CX=CY=0).

         call do_CS2GRD_sg(PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,C,CX,CY,IER)
      end subroutine interp_CS2GRD_sg



! ***********************************************************************
! ***********************************************************************
! ***********************************************************************
! ***********************************************************************


! bicubic splines and hermite

! from PSPLINE by Doug McCune (version as of February, 2004)
! PSPLINE Home Page:
! http://w3.pppl.gov/NTCC/PSPLINE
! Doug McCune, Princeton University
!     dmccune@pppl.gov
    
! use interp_mkbicub_sg to set up the interpolation information
! use interp_evbicub_sg to evaluate it
! see interp_2d_lib_db for double precision versions.

      subroutine interp_mkbicub_sg(x,nx,y,ny,f1,nf2,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ibcymin,bcymin,ibcymax,bcymax,
     >   ilinx,iliny,ier)

         use bicub_sg

         integer, intent(in) :: nx                        ! length of x vector
         integer, intent(in) :: ny                        ! length of y vector
         real, intent(in) :: x(:) ! (nx)                        ! x vector, strict ascending
         real, intent(in) :: y(:) ! (ny)                        ! y vector, strict ascending
         integer, intent(in) :: nf2                       ! 2nd dimension of f, nf2.ge.nx
         real, intent(inout), pointer :: f1(:) ! =(4,nf2,ny)               ! data & spline coefficients
!
!  on input:  f(1,i,j) = f(x(i),y(j))
!  on output:  f(1,i,j) unchanged
!              f(2,i,j) = d2f/dx2(x(i),y(j))
!              f(3,i,j) = d2f/dy2(x(i),y(j))
!              f(4,i,j) = d4f/dx2dy2(x(i),y(j))
!
!  and the interpolation formula for (x,y) in (x(i),x(i+1))^(y(j),y(j+1))
!  is:
!        hx = x(i+1)-x(i)   hy = y(j+1)-y(j)
!        dxp= (x-x(i))/hx   dxm= 1-dxp     dxp,dxm in (0,1)
!        dyp= (x-x(i))/hx   dym= 1-dyp     dyp,dym in (0,1)
!        dx3p = dxp**3-dxp  dx3m = dxm**3-dxm     dxp3,dxm3 in (0,1)
!
!   finterp = dxm*(dym*f(1,i,j)+dyp*f(1,i,j+1))
!            +dxp*(dym*f(1,i+1,j)+dyp*f(1,i+1,j+1))
!     +1/6*hx**2*
!            dx3m*(dym*f(2,i,j)+dyp*f(2,i,j+1))
!           +dx3p*(dym*f(2,i+1,j)+dyp*f(2,i+1,j+1))
!     +1/6*hy**2*
!            dxm*(dy3m*f(3,i,j)+dy3p*f(3,i,j+1))
!           +dxp*(dy3m*f(3,i+1,j)+dy3p*f(3,i+1,j+1))
!     +1/36*hx**2*hy**2*
!            dx3m*(dym*f(4,i,j)+dyp*f(4,i,j+1))
!           +dx3p*(dym*f(4,i+1,j)+dyp*f(4,i+1,j+1))
!
!  where the f(2:4,*,*) are cleverly evaluated to assure
!  (a) finterp is continuous and twice differentiable across all
!      grid cell boundaries, and
!  (b) all boundary conditions are satisfied.
!
!  input bdy condition data:
      integer, intent(in) :: ibcxmin                   ! bc flag for x=xmin
      real, intent(in) :: bcxmin(:) ! (ny)                   ! bc data vs. y at x=xmin
      integer, intent(in) :: ibcxmax                   ! bc flag for x=xmax
      real, intent(in) :: bcxmax(:) ! (ny)                   ! bc data vs. y at x=xmax
      integer, intent(in) :: ibcymin                   ! bc flag for y=ymin
      real, intent(in) :: bcymin(:) ! (nx)                   ! bc data vs. x at y=ymin
      integer, intent(in) :: ibcymax                   ! bc flag for y=ymax
      real, intent(in) :: bcymax(:) ! (nx)                   ! bc data vs. x at y=ymax
!
!  with interpretation:
!   ibcxmin -- indicator for boundary condition at x(1):
!    bcxmin(...) -- boundary condition data
!     =-1 -- periodic boundary condition
!     =0 -- use "not a knot"
!     =1 -- match slope, specified at x(1),th(ith) by bcxmin(ith)
!     =2 -- match 2nd derivative, specified at x(1),th(ith) by bcxmin(ith)
!     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all th(j)
!     =4 -- boundary condition is d2f/dx2=0 at x(1), all th(j)
!     =5 -- match 1st derivative to 1st divided difference
!     =6 -- match 2nd derivative to 2nd divided difference
!     =7 -- match 3rd derivative to 3rd divided difference
!           (for more detailed definition of B!s 5-7, see the
!           comments of subroutine mkspline)
!   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
!
!   ibcxmax -- indicator for boundary condition at x(nx):
!    bcxmax(...) -- boundary condition data
!     (interpretation as with ibcxmin, bcxmin)
!   NOTE:  if ibcxmin=-1, ibcxmax is ignored! ...and the B! is periodic.
!
!  and analogous interpretation for ibcymin,bcymin,ibcymax,bcymax
!  (df/dy or d2f/dy2 boundary conditions at y=ymin and y=ymax).
!
!  output linear grid flags and completion code (ier=0 is normal):
!
      integer, intent(out) :: ilinx                     ! =1: x grid is "nearly" equally spaced
      integer, intent(out) :: iliny                     ! =1: y grid is "nearly" equally spaced
!  ilinx and iliny are set to zero if corresponding grids are not equally
!  spaced
!
      integer, intent(out) :: ier                       ! =0 on exit if there is no error.
!
!  if there is an error, ier is set and a message is output on unit 6.
!  these are considered programming errors in the calling routine.
!
!  possible errors:
!    x(...) not strict ascending
!    y(...) not strict ascending
!    nx .lt. 4
!    ny .lt. 4
!    invalid boundary condition flag
!
!-----------------------

         call do_mkbicub(x,nx,y,ny,f1,nf2,
     >         ibcxmin,bcxmin,ibcxmax,bcxmax,
     >         ibcymin,bcymin,ibcymax,bcymax,
     >         ilinx,iliny,ier)

      end subroutine interp_mkbicub_sg
      
      
      
      subroutine interp_evbicub_sg(xget,yget,x,nx,y,ny,ilinx,iliny,f1,inf2,ict,fval,ier)

      use bicub_sg

      ! evaluate bicubic spline interpolant (single precision version)
!
!  use mkbicub to set up spline coefficients
!
!  evaluate a 2d cubic Spline interpolant on a rectilinear
!  grid -- this is C2 in both directions.
!
!  this subroutine calls two subroutines:
!     herm2xy  -- find cell containing (xget,yget)
!     fvbicub  -- evaluate interpolant function and (optionally) derivatives
      
!  input arguments:
!  ================
!
      integer, intent(in) :: nx,ny                     ! grid sizes
      real, intent(in) :: xget,yget                    ! target of this interpolation
      real, intent(in) :: x(:) ! (nx)                        ! ordered x grid
      real, intent(in) :: y(:) ! (ny)                        ! ordered y grid
      integer, intent(in) :: ilinx                     ! ilinx=1 => assume x evenly spaced
      integer, intent(in) :: iliny                     ! iliny=1 => assume y evenly spaced
!
      integer, intent(in) :: inf2
      real, intent(inout), pointer :: f1(:) ! function data
!
!       f 2nd dimension inf2 must be .ge. nx
!       contents of f:
!
!  f(0,i,j) = f @ x(i),y(j)
!  f(1,i,j) = d2f/dx2 @ x(i),y(j)
!  f(2,i,j) = d2f/dy2 @ x(i),y(j)
!  f(3,i,j) = d4f/dx2dy2 @ x(i),y(j)
!
!      (these are spline coefficients selected for continuous 2-
!      diffentiability, see mkbicub[w].for)
!
      integer, intent(in) :: ict(6)                    ! code specifying output desired
!
!  ict(1)=1 -- return f  (0, don't)
!  ict(2)=1 -- return df/dx  (0, don't)
!  ict(3)=1 -- return df/dy  (0, don't)
!  ict(4)=1 -- return d2f/dx2  (0, don't)
!  ict(5)=1 -- return d2f/dy2  (0, don't)
!  ict(6)=1 -- return d2f/dxdy (0, don't)
!                   the number of non zero values ict(1:6)
!                   determines the number of outputs...
!
!  new dmc December 2005 -- access to higher derivatives (even if not
!  continuous-- but can only go up to 3rd derivatives on any one coordinate.
!     if ict(1)=3 -- want 3rd derivatives
!          ict(2)=1 for d3f/dx3
!          ict(3)=1 for d3f/dx2dy
!          ict(4)=1 for d3f/dxdy2
!          ict(5)=1 for d3f/dy3
!               number of non-zero values ict(2:5) gives no. of outputs
!     if ict(1)=4 -- want 4th derivatives
!          ict(2)=1 for d4f/dx3dy
!          ict(3)=1 for d4f/dx2dy2
!          ict(4)=1 for d4f/dxdy3
!               number of non-zero values ict(2:4) gives no. of outputs
!     if ict(1)=5 -- want 5th derivatives
!          ict(2)=1 for d5f/dx3dy2
!          ict(3)=1 for d5f/dx2dy3
!               number of non-zero values ict(2:3) gives no. of outputs
!     if ict(1)=6 -- want 6th derivatives
!          d6f/dx3dy3 -- one value is returned.
!
! output arguments:
! =================
!
      real, intent(inout) :: fval(6)                      ! output data
      integer, intent(out) :: ier                       ! error code =0 ==> no error
!
!  fval(1) receives the first output (depends on ict(...) spec)
!  fval(2) receives the second output (depends on ict(...) spec)
!  fval(3) receives the third output (depends on ict(...) spec)
!  fval(4) receives the fourth output (depends on ict(...) spec)
!  fval(5) receives the fourth output (depends on ict(...) spec)
!  fval(6) receives the fourth output (depends on ict(...) spec)
!
!  examples:
!    on input ict = [1,1,1,0,0,1]
!   on output fval= [f,df/dx,df/dy,d2f/dxdy], elements 5 & 6 not referenced.
!
!    on input ict = [1,0,0,0,0,0]
!   on output fval= [f] ... elements 2 -- 6 never referenced.
!
!    on input ict = [0,0,0,1,1,0]
!   on output fval= [d2f/dx2,d2f/dy2] ... elements 3 -- 6 never referenced.
!
!    on input ict = [0,0,1,0,0,0]
!   on output fval= [df/dy] ... elements 2 -- 6 never referenced.
!
!  ier -- completion code:  0 means OK
!-------------------

      integer i,j
      integer ii(1), jj(1) 
      real xparam(1),yparam(1),hx(1),hxi(1),hy(1),hyi(1)
      call herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny,i,j,
     >         xparam(1),yparam(1),hx(1),hxi(1),hy(1),hyi(1),ier)
      if(ier.ne.0) return
      ii(1) = i
      jj(1) = j
      call fvbicub(ict,1,1,fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,f1,inf2,ny)

      end subroutine interp_evbicub_sg
      


! ***********************************************************************
! ***********************************************************************

! 2d piecewise monotonic interpolation -- values only, no slopes.
! does 4 1d interpolations in x followed by 1 1d interpolation in y.
! use interp_mkbipm_db to set up the interpolation information
! use interp_evbipm_db to evaluate it
      
      subroutine interp_mkbipm_sg(x,nx,y,ny,f1,nf2,ier)
         use bipm_sg, only: do_mkbipm_sg
         integer, intent(in) :: nx                        ! length of x vector
         integer, intent(in) :: ny                        ! length of y vector
         real, intent(in), pointer :: x(:) ! (nx)            ! x vector, strict ascending
         real, intent(in), pointer :: y(:) ! (ny)            ! y vector, strict ascending
         integer, intent(in) :: nf2                       ! 2nd dimension of f, nf2.ge.nx
         real, intent(inout), pointer :: f1(:) ! =(4,nf2,ny)   ! data & interpolant coefficients
         integer, intent(out) :: ier                      ! =0 on exit if there is no error.   
         call do_mkbipm_sg(x,nx,y,ny,f1,nf2,ier)
      end subroutine interp_mkbipm_sg


      subroutine interp_evbipm_sg(xget,yget,x,nx,y,ny,f1,nf2,z,ier)
         use bipm_sg, only: do_evbipm_sg
         integer, intent(in) :: nx,ny
         real, intent(in) :: xget,yget        ! target of this interpolation
         real, intent(in), pointer :: x(:) ! (nx)            ! ordered x grid
         real, intent(in), pointer :: y(:) ! (ny)            ! ordered y grid
         integer, intent(in) :: nf2
         real, intent(in), pointer :: f1(:) ! =(4,nf2,ny)      ! function data
         real, intent(out) :: z
         integer, intent(out) :: ier                      ! error code =0 ==> no error
         call do_evbipm_sg(xget,yget,x,nx,y,ny,f1,nf2,z,ier)   
      end subroutine interp_evbipm_sg



      end module interp_2d_lib_sg
