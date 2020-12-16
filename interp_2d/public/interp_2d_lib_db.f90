! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module interp_2d_lib_db
      ! real(dp) library for 2d interpolation
      
      use const_def, only: dp

      implicit none

      contains ! the procedure interface for the library
      ! client programs should only call these routines.


! contents

! rectangular-grid of data points
   ! interp_rgbi3p_db -- point interpolation (akima)
   ! interp_rgsf3p_db -- surface interpolation (akima)
   ! interp_mkbicub_db -- bicubic splines
   ! interp_mkbipm_db -- 2d piecewise monotonic
      
! scattered set of data points
   ! interp_cs2val_db -- point interpolation (renka)
   ! interp_cs2grd_db -- point interpolation with gradients (renka)



! ***********************************************************************
! ***********************************************************************
! ***********************************************************************
! ***********************************************************************



* rectangular-grid bivariate interpolation and surface fitting
* from acm algorithm 760., acm trans. math. software (22) 1996, 357-361.
* hiroshi akima
* u.s. department of commerce, ntia/its
* version of 1995/08

      ! interpolate at given points
      subroutine interp_rgbi3p_db(md,nxd,nyd,xd,yd,zd,nip,xi,yi,zi,ier,wk)
         integer, intent(in) :: md, nxd, nyd, nip
         real(dp), intent(in) :: xd(nxd), yd(nyd), zd(nxd,nyd), xi(nip), yi(nip)
         real(dp), intent(inout) :: zi(nip), wk(3,nxd,nyd)
         integer, intent(out) :: ier
         call do_rgbi3p_db(md,nxd,nyd,xd,yd,zd,nip,xi,yi,zi,ier,wk)
      end subroutine interp_rgbi3p_db

      ! interpolate surface for given grid in x and y
      subroutine interp_rgsf3p_db(md,nxd,nyd,xd,yd,zd,nxi,xi,nyi,yi,zi,ier,wk)
         integer, intent(in) :: md, nxd, nyd, nxi, nyi
         real(dp), intent(in) :: xd(nxd), yd(nyd), zd(nxd,nyd), xi(nxi), yi(nyi)
         real(dp), intent(inout) :: zi(nxi,nyi), wk(3,nxd,nyd)
         integer, intent(out) :: ier
         call do_rgsf3p_db(md,nxd,nyd,xd,yd,zd,nxi,xi,nyi,yi,zi,ier,wk)
      end subroutine interp_rgsf3p_db
      
      

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************
! ***********************************************************************



! cubic shepard method for bivariate interpolation of scattered data.
! from acm algorithm 790., acm trans. math. software (25) 1999, 70-73.
! robert j. renka
! dept. of computer science
! univ. of north texas
! renka@cs.unt.edu

! use cshep2_db to set up the interpolation information.
! use cs2val_db to evaluate it.
! use cs2grd_db to get value and derivatives.

! detailed documentation can be found in interp_2d_lib_sg.f
! these routines are exactly the same except they are real(dp).

      ! see cshep2_sg for documentation details.
      subroutine interp_cshep2_db(n,x,y,f,nc,nw,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rw,a,ier)
         integer, intent(in) :: n, nc, nw, nr
         integer, intent(out) :: lcell(nr,nr), lnext(n), ier
         real(dp), intent(in) :: x(n), y(n), f(n)
         real(dp), intent(inout) :: xmin, ymin, dx, dy, rmax, rw(n), a(9,n)
         call do_cshep2_db(n,x,y,f,nc,nw,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rw,a,ier)
      end subroutine interp_cshep2_db

      ! see cs2val_sg for documentation details.
      real(dp) function interp_cs2val_db(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rw,a,ier)
         integer, intent(in) :: n, nr, lcell(nr,nr), lnext(n)
         real(dp), intent(in) :: px, py, x(n), y(n), f(n), xmin, ymin, dx, dy, rmax, rw(n), a(9,n)
         integer, intent(out) :: ier
         real(dp) :: do_cs2val_db        
         interp_cs2val_db = do_cs2val_db(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rw,a,ier)   
      end function interp_cs2val_db

      ! see cs2grd_sg for documentation details.
      subroutine interp_cs2grd_db(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rw,a,c,cx,cy,ier)
         integer, intent(in) :: n, nr, lcell(nr,nr), lnext(n)
         real(dp), intent(in) :: px, py, x(n), y(n), f(n), xmin, ymin, dx, dy, rmax, rw(n), a(9,n)
         real(dp), intent(out) :: c, cx, cy
         integer, intent(out) :: ier
         call do_cs2grd_db(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rw,a,c,cx,cy,ier)
      end subroutine interp_cs2grd_db



! ***********************************************************************
! ***********************************************************************


! bicubic splines
! use interp_mkbicub_db to set up the interpolation information
! use interp_evbicub_db to evaluate it

! detailed documentation can be found in interp_2d_lib_sg.f
      
      ! see interp_mkbicub_sg for documentation details.
      subroutine interp_mkbicub_db(x,nx,y,ny,f1,nf2,
     >         ibcxmin,bcxmin,ibcxmax,bcxmax,
     >         ibcymin,bcymin,ibcymax,bcymax,
     >         ilinx,iliny,ier)
         use bicub_db, only: do_mkbicub_db
         integer, intent(in) :: nx                        ! length of x vector
         integer, intent(in) :: ny                        ! length of y vector
         real(dp), intent(in) :: x(:) ! (nx)            ! x vector, strict ascending
         real(dp), intent(in) :: y(:) ! (ny)            ! y vector, strict ascending
         integer, intent(in) :: nf2                       ! 2nd dimension of f, nf2.ge.nx
         real(dp), intent(inout), pointer :: f1(:) ! =(4,nf2,ny)   ! data & spline coefficients
         integer, intent(in) :: ibcxmin                   ! bc flag for x=xmin
         real(dp), intent(in) :: bcxmin(:) ! (ny)       ! bc data vs. y at x=xmin
         integer, intent(in) :: ibcxmax                   ! bc flag for x=xmax
         real(dp), intent(in) :: bcxmax(:) ! (ny)       ! bc data vs. y at x=xmax
         integer, intent(in) :: ibcymin                   ! bc flag for y=ymin
         real(dp), intent(in) :: bcymin(:) ! (nx)       ! bc data vs. x at y=ymin
         integer, intent(in) :: ibcymax                   ! bc flag for y=ymax
         real(dp), intent(in) :: bcymax(:) ! (nx)       ! bc data vs. x at y=ymax
         integer, intent(out) :: ilinx                    ! =1: x grid is "nearly" equally spaced
         integer, intent(out) :: iliny                    ! =1: y grid is "nearly" equally spaced
         integer, intent(out) :: ier                      ! =0 on exit if there is no error.   
         call do_mkbicub_db(x,nx,y,ny,f1,nf2,
     >         ibcxmin,bcxmin,ibcxmax,bcxmax,
     >         ibcymin,bcymin,ibcymax,bcymax,
     >         ilinx,iliny,ier)
         end subroutine interp_mkbicub_db


      ! see interp_evbicub_sg for documentation details.
      subroutine interp_evbicub_db(xget,yget,x,nx,y,ny,ilinx,iliny,f1,inf2,ict,fval,ier)
         use bicub_db, only: fvbicub_db, herm2xy_db
         integer, intent(in) :: nx,ny                     ! grid sizes
         real(dp), intent(in) :: xget,yget        ! target of this interpolation
         real(dp), intent(in) :: x(:) ! (nx)            ! ordered x grid
         real(dp), intent(in) :: y(:) ! (ny)            ! ordered y grid
         integer, intent(in) :: ilinx                     ! ilinx=1 => assume x evenly spaced
         integer, intent(in) :: iliny                     ! iliny=1 => assume y evenly spaced
         integer, intent(in) :: inf2
         real(dp), intent(inout), pointer :: f1(:) ! function data
         integer, intent(in) :: ict(6)                    ! code specifying output desired
         real(dp), intent(inout) :: fval(6)         ! output data
         integer, intent(out) :: ier                      ! error code =0 ==> no error
         integer i,j
         integer ii(1), jj(1) 
         real(dp) xparam(1),yparam(1),hx(1),hxi(1),hy(1),hyi(1)
 1       format(a40,1pe26.16)
         call herm2xy_db(xget,yget,x,nx,y,ny,ilinx,iliny,i,j,
     >         xparam(1),yparam(1),hx(1),hxi(1),hy(1),hyi(1),ier)
         if (ier .ne. 0) return
         ii(1) = i
         jj(1) = j
         call fvbicub_db(ict,1,1,fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,f1,inf2,ny)   
      end subroutine interp_evbicub_db
      
      
      ! this is used by do_mkbicub_db to get 2nd derivatives d_dx2 and d_dy2
      subroutine interp_mkspline_db(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,ilinx,ier)
         use bicub_db, only: mkspline_db
         integer, intent(in) :: nx ! no. of data points
         real(dp), intent(in) :: x(nx) ! x axis data, strict ascending order
         real(dp), intent(inout) :: fspl(2,nx) ! f(1,*): data in; f(2,*): coeffs out
         integer, intent(in) :: ibcxmin                   ! b.c. flag @ x=xmin=x(1)
         real(dp), intent(in) :: bcxmin                   ! b.c. data @xmin
         integer, intent(in) :: ibcxmax                   ! b.c. flag @ x=xmax=x(nx)
         real(dp), intent(in) :: bcxmax                   ! b.c. data @xmax
         integer, intent(in) :: ilinx                     ! ilinx=1 => assume x evenly spaced
         integer, intent(out) :: ier                      ! error code =0 ==> no error
         call mkspline_db(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,ilinx,ier)
      end subroutine interp_mkspline_db


! ***********************************************************************
! ***********************************************************************

! 2d piecewise monotonic interpolation -- values only, no slopes.
! does 4 1d interpolations in x followed by 1 1d interpolation in y.
! use interp_mkbipm_db to set up the interpolation information
! use interp_evbipm_db to evaluate it
      
      subroutine interp_mkbipm_db(x,nx,y,ny,f1,nf2,ier)
         use bipm_db, only: do_mkbipm_db
         integer, intent(in) :: nx                        ! length of x vector
         integer, intent(in) :: ny                        ! length of y vector
         real(dp), intent(in), pointer :: x(:) ! (nx)            ! x vector, strict ascending
         real(dp), intent(in), pointer :: y(:) ! (ny)            ! y vector, strict ascending
         integer, intent(in) :: nf2                       ! 2nd dimension of f, nf2.ge.nx
         real(dp), intent(inout), pointer :: f1(:) ! =(4,nf2,ny)   ! data & interpolant coefficients
         integer, intent(out) :: ier                      ! =0 on exit if there is no error.   
         call do_mkbipm_db(x,nx,y,ny,f1,nf2,ier)
      end subroutine interp_mkbipm_db


      subroutine interp_evbipm_db(xget,yget,x,nx,y,ny,f1,nf2,z,ier)
         use bipm_db, only: do_evbipm_db
         integer, intent(in) :: nx,ny
         real(dp), intent(in) :: xget,yget        ! target of this interpolation
         real(dp), intent(in), pointer :: x(:) ! (nx)            ! ordered x grid
         real(dp), intent(in), pointer :: y(:) ! (ny)            ! ordered y grid
         integer, intent(in) :: nf2
         real(dp), intent(in), pointer :: f1(:) ! =(4,nf2,ny)      ! function data
         real(dp), intent(out) :: z
         integer, intent(out) :: ier                      ! error code =0 ==> no error
         call do_evbipm_db(xget,yget,x,nx,y,ny,f1,nf2,z,ier)   
      end subroutine interp_evbipm_db


! ***********************************************************************
! ***********************************************************************


      end module interp_2d_lib_db
