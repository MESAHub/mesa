! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module bicub_db
      
      use const_def, only: dp

      !implicit none


      contains

!        from PSPLINE by Doug McCune (version as of February, 2004)


!        PSPLINE Home Page:
!        http://w3.pppl.gov/NTCC/PSPLINE

!        Doug McCune, Princeton University
!                dmccune@pppl.gov

!
!  bcspeval -- eval bicubic spline function and/or derivatives
!
      subroutine bcspeval_db(xget,yget,iselect,fval,
     >                    x,nx,y,ny,ilinx,iliny,f,inf3,ier)
!
      implicit none
      integer iselect(6)
      integer ilinx,iliny,nx,ny,inf3,ier
!
      real(dp) xget,yget
      real(dp) fval(6)
      real(dp) x(nx),y(ny),f(4,4,inf3,ny)
!
!  modification -- dmc 11 Jan 1999 -- remove SAVE stmts;
C    break routine into these parts:
C
C    bcspevxy -- find grid cell of target pt.
C    bcspevfn -- evaluate function using output of bcpsevxy
C
C    in cases where multiple functions are defined on the same grid,
C    time can be saved by using bcspevxy once and then bcspevfn
C    multiple times.
!
!  input:
!     (xget,yget)   location where interpolated value is desired
!                   x(1).le.xget.le.x(nx) expected
!                   y(1).le.yget.le.y(ny) expected
!
!     iselect       select desired output
!
!                     iselect(1)=1 -- want function value (f) itself
!                     iselect(2)=1 -- want  df/dx
!                     iselect(3)=1 -- want  df/dy
!                     iselect(4)=1 -- want  d2f/dx2
!                     iselect(5)=1 -- want  d2f/dy2
!                     iselect(6)=1 -- want  d2f/dxdy
!
!              example:  iselect(1)=iselect(2)=iselect(3)=1
!                            f, df/dx, and df/dy all evaluated
!                        iselect(4)=iselect(5)=iselect(6)=0
!                            2nd derivatives not evaluated.
!
!                   the number of non zero values iselect(1:6)
!                   determines the number of outputs...
!                   see fval (output) description.
!
!  new dmc December 2005 -- access to higher derivatives (even if not
!  continuous-- but can only go up to 3rd derivatives on any one coordinate.
!     if iselect(1)=3 -- want 3rd derivatives
!          iselect(2)=1 for d3f/dx3
!          iselect(3)=1 for d3f/dx2dy
!          iselect(4)=1 for d3f/dxdy2
!          iselect(5)=1 for d3f/dy3
!               number of non-zero values iselect(2:5) gives no. of outputs
!     if iselect(1)=4 -- want 4th derivatives
!          iselect(2)=1 for d4f/dx3dy
!          iselect(3)=1 for d4f/dx2dy2
!          iselect(4)=1 for d4f/dxdy3
!               number of non-zero values iselect(2:4) gives no. of outputs
!     if iselect(1)=5 -- want 5th derivatives
!          iselect(2)=1 for d5f/dx3dy2
!          iselect(3)=1 for d5f/dx2dy3
!               number of non-zero values iselect(2:3) gives no. of outputs
!     if iselect(1)=6 -- want 6th derivatives
!          d6f/dx3dy3 -- one value is returned.
!
!     x(1...nx)     independent coordinate x, strict ascending
!     y(1...ny)     independent coordinate y, strict ascending
!
!     ilinx  --  =1: flag that x is linearly spaced (avoid search for speed)
!     iliny  --  =1: flag that y is linearly spaced (avoid search for speed)
!
!  **CAUTION** actual even spacing of x, y is NOT CHECKED HERE!
!
!
!     f             the function values (at grid points) and spline coefs
!
!  evaluation formula:  for point x btw x(i) and x(i+1), dx=x-x(i)
!                             and y btw y(j) and y(j+1), dy=y-y(j),
!
!      spline value =
!        f(1,1,i,j) + dx*f(2,1,i,j) + dx**2*f(3,1,i,j) + dx**3*f(4,1,i,j)
!   +dy*(f(1,2,i,j) + dx*f(2,2,i,j) + dx**2*f(3,2,i,j) + dx**3*f(4,2,i,j))
!   +d2*(f(1,3,i,j) + dx*f(2,3,i,j) + dx**2*f(3,3,i,j) + dx**3*f(4,3,i,j))
!   +d3*(f(1,4,i,j) + dx*f(2,4,i,j) + dx**2*f(3,4,i,j) + dx**3*f(4,4,i,j))
!
!      where d2=dy**2 and d3=dy**3.
!
!  output:
!      up to 6 elements of fval, ordered as follows:
!        fval(1)=function value or lowest order derivative requested
!        fval(2)=next order derivative
!             etc
!        the ordering is a subset of the sequence given under the "iselect"
!        description.
!
!      ier = 0 -- successful completion; = 1 -- an error occurred.
!
!-------------------------------------------------------------------
!  local
!
      integer :: i=0
      integer :: j=0
!
      real(dp) dx,dy
!
!--------------------------
!
      call bcspevxy_db(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i,j,dx,dy,ier)
      if(ier.ne.0) return
!
      call bcspevfn_db(iselect,1,1,fval,(/i/),(/j/),
     <   (/dx/),(/dy/),f,inf3,ny)
!
      return
      end subroutine bcspeval_db

!
!-------------------------------------------------------------------------
!  bcspevxy -- look up x-y zone
!
!  this is the "first part" of bcspeval, see comments, above.
!
      subroutine bcspevxy_db(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i,j,dx,dy,ier)
!
      integer nx,ny                     ! array dimensions
!
      real(dp) xget,yget                    ! target point
      real(dp) x(nx),y(ny)                  ! indep. coords.
!
      integer ilinx                     ! =1:  assume x evenly spaced
      integer iliny                     ! =1:  assume y evenly spaced
!
!  output of bcspevxy
!
      integer i,j                       ! index to cell containing target pt
      real(dp) dx,dy                        ! displacement of target pt w/in cell
                                        ! dx=x-x(i)  dy=y-y(j)
C
      integer ier                       ! return ier.ne.0 on error
!
!------------------------------------
!
      real(dp) zxget, zyget
      ier=0
!
!  range check
!
      zxget=xget
      zyget=yget
 
      if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
         zxtol=4.0d-7*max(abs(x(1)),abs(x(nx)))
         if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
            ier=1
!            write(6,1001) xget,x(1),x(nx)
! 1001       format(' ?bcspeval:  xget=',1pe11.4,' out of range ',
!     >         1pe11.4,' to ',1pe11.4)
         else
!            if((xget.lt.x(1)-0.5*zxtol).or.
!     >         (xget.gt.x(nx)+0.5*zxtol))
!     >      write(6,1011) xget,x(1),x(nx)
! 1011       format(' %bcspeval:  xget=',1pe15.8,' beyond range ',
!     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(xget.lt.x(1)) then
               zxget=x(1)
            else
               zxget=x(nx)
            endif
         endif
      endif
      if((yget.lt.y(1)).or.(yget.gt.y(ny))) then
         zytol=4.0d-7*max(abs(y(1)),abs(y(ny)))
         if((yget.lt.y(1)-zytol).or.(yget.gt.y(ny)+zytol)) then
            ier=1
!            write(6,1002) yget,y(1),y(ny)
! 1002       format(' ?bcspeval:  yget=',1pe11.4,' out of range ',
!     >         1pe11.4,' to ',1pe11.4)
         else
!         if((yget.lt.y(1)-0.5*zytol).or.(yget.gt.y(ny)+0.5*zytol))
!     >      write(6,1012) yget,y(1),y(ny)
! 1012       format(' %bcspeval:  yget=',1pe15.8,' beyond range ',
!     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(yget.lt.y(1)) then
               zyget=y(1)
            else
               zyget=y(ny)
            endif
         endif
      endif
      if(ier.ne.0) return
!
!  now find interval in which target point lies..
!
      nxm=nx-1
      nym=ny-1
!
      if(ilinx.eq.1) then
         ii=1+nxm*(zxget-x(1))/(x(nx)-x(1))
         i=min(nxm, ii)
         if(zxget.lt.x(i)) then
            i=i-1
         else if(zxget.gt.x(i+1)) then
            i=i+1
         endif
      else
         if((1.le.i).and.(i.lt.nxm)) then
            if((x(i).le.zxget).and.(zxget.le.x(i+1))) then
               continue  ! already have the zone
            else
               call zonfind_db(x,nx,zxget,i)
            endif
         else
            call zonfind_db(x,nx,zxget,i)
         endif
      endif
!
      if(iliny.eq.1) then
         jj=1+nym*(zyget-y(1))/(y(ny)-y(1))
         j=min(nym, jj)
         if(zyget.lt.y(j)) then
            j=j-1
         else if(zyget.gt.y(j+1)) then
            j=j+1
         endif
      else
         if((1.le.j).and.(j.lt.nym)) then
            if((y(j).le.zyget).and.(zyget.le.y(j+1))) then
               continue  ! already have the zone
            else
               call zonfind_db(y,ny,zyget,j)
            endif
         else
            call zonfind_db(y,ny,zyget,j)
         endif
      endif
!
      dx=zxget-x(i)
      dy=zyget-y(j)
!
      return
      end subroutine bcspevxy_db

!------------------------------------------------------------------------
!  bcspevfn -- OK now evaluate the bicubic spline
!
      subroutine bcspevfn_db(ict,ivec,ivd,fval,iv,jv,dxv,dyv,f,inf3,ny)
!
!  input:
      integer ny
      integer ict(6)                    ! selector:
!        ict(1)=1 for f      (don't evaluate f if ict(1)=0)
!        ict(2)=1 for df/dx   ""
!        ict(3)=1 for df/dy   ""
!        ict(4)=1 for d2f/dx2
!        ict(5)=1 for d2f/dy2
!        ict(6)=1 for d2f/dxdy
!
!    note:  if ict(1)=-1, evaluate f,d2f/dx2,d2f/dy2,d4f/dx2dy2
!
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
      integer ivec,ivd                  ! vector dimensioning
!
!    ivec-- number of vector pts (spline values to look up)
!    ivd -- 1st dimension of fval, .ge.ivec
!
! output:
      real(dp) fval(ivd,6)                 ! output array
!
!    v = index to element in vector;
!  fval(v,1) = first item requested by ict(...),
!  fval(v,2) = 2nd item requested,  ...etc...
!
!  input:
      integer iv(ivec),jv(ivec)         ! grid cell indices -- vectors
      real(dp) dxv(ivec),dyv(ivec)          ! displacements w/in cell -- vectors
!
      integer inf3                      ! 3rd dimension of f -- .ge. nx
      real(dp) f(4,4,inf3,ny)               ! bicubic fcn spline coeffs array
!
!  usage example:
!
!  1.  for each element (xx(v),yy(v)) in a vector of (x,y) pairs,
!    find the x and y zone indices and displacements with respect
!    to the "lower left corner" of the zone; store these in vectors
!    iv,jv and dxv,dyv.
!
!  2.  set ict(1)=0, ict(2)=1, ict(3)=1, the rest zero -- get only
!      the 1st derivatives.
!
!  3.  ivec is the length of the vector; ivd is the 1st dimension
!      of the array fval to receive the output
!
!      real(dp) fval(ivd,6)
!      real(dp) xv(ivd),yv(ivd)
!      integer iv(ivd),jv(ivd)
!      real(dp) dxv(ivd),dyv(ivd)
!      integer ict(6)
!
!      real(dp) fspline(4,4,nx,ny)  ! spline coeffs
!      data ict/0,1,1,0,0,0/    ! this call:  want 1st derivatives
!                               ! only ... these will be output to
!                               ! fval(*,1) fval(*,2)
!      ...
!      do iv=1,ivec
!        ...                    ! find indices and displacements
!      enddo
!      call bcspevfn(ict,ivec,ivd,fval,iv,jv,dxv,dyv,fspline,nx,ny)
!
!-------------------------------------------------------------------
!  local:
!
      integer v                         ! vector element index
!
!  OK can now do evaluations
!
      iaval=0  ! fval addressing
!
      if(ict(1).le.2) then
         if((ict(1).gt.0).or.(ict(1).eq.-1)) then
!  evaluate f
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >       f(1,1,i,j)+dy*(f(1,2,i,j)+dy*(f(1,3,i,j)+dy*f(1,4,i,j)))
     >  +dx*(f(2,1,i,j)+dy*(f(2,2,i,j)+dy*(f(2,3,i,j)+dy*f(2,4,i,j)))
     >  +dx*(f(3,1,i,j)+dy*(f(3,2,i,j)+dy*(f(3,3,i,j)+dy*f(3,4,i,j)))
     >  +dx*(f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j))))))
            enddo
         endif
!
         if((ict(2).gt.0).and.(ict(1).ne.-1)) then
!  evaluate df/dx
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >         f(2,1,i,j)+dy*(f(2,2,i,j)+dy*(f(2,3,i,j)+dy*f(2,4,i,j)))
     >       +2.d0*dx*(
     >         f(3,1,i,j)+dy*(f(3,2,i,j)+dy*(f(3,3,i,j)+dy*f(3,4,i,j)))
     >       +1.5d0*dx*(
     >         f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j)))
     >              ))
            enddo
         endif
!
         if((ict(3).gt.0).and.(ict(1).ne.-1)) then
!  evaluate df/dy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >         f(1,2,i,j)+dy*(2.d0*f(1,3,i,j)+dy*3.d0*f(1,4,i,j))
     >      +dx*(f(2,2,i,j)+dy*(2.d0*f(2,3,i,j)+dy*3.d0*f(2,4,i,j))
     >      +dx*(f(3,2,i,j)+dy*(2.d0*f(3,3,i,j)+dy*3.d0*f(3,4,i,j))
     >      +dx*(f(4,2,i,j)+dy*(2.d0*f(4,3,i,j)+dy*3.d0*f(4,4,i,j))
     >              )))
            enddo
         endif
!
         if((ict(4).gt.0).or.(ict(1).eq.-1)) then
!  evaluate d2f/dx2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >              2.d0*(
     >              f(3,1,i,j)+dy*(f(3,2,i,j)+dy*(f(3,3,i,j)+dy*f(3,4,i,j))))
     >              +6.d0*dx*(
     >              f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j))))
            enddo
         endif
!
         if((ict(5).gt.0).or.(ict(1).eq.-1)) then
!  evaluate d2f/dy2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >              2.d0*f(1,3,i,j)+6.d0*dy*f(1,4,i,j)
     >              +dx*(2.d0*f(2,3,i,j)+6.d0*dy*f(2,4,i,j)
     >              +dx*(2.d0*f(3,3,i,j)+6.d0*dy*f(3,4,i,j)
     >              +dx*(2.d0*f(4,3,i,j)+6.d0*dy*f(4,4,i,j))))
            enddo
         endif
!
         if((ict(6).gt.0).and.(ict(1).ne.-1)) then
!  evaluate d2f/dxdy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >            f(2,2,i,j)+dy*(2.d0*f(2,3,i,j)+dy*3.d0*f(2,4,i,j))
     > +2.d0*dx*(f(3,2,i,j)+dy*(2.d0*f(3,3,i,j)+dy*3.d0*f(3,4,i,j))
     >+1.5d0*dx*(f(4,2,i,j)+dy*(2.d0*f(4,3,i,j)+dy*3.d0*f(4,4,i,j))
     >              ))
            enddo
         endif
!
         if(ict(1).eq.-1) then
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >              4.d0*f(3,3,i,j)+12.d0*dy*f(3,4,i,j)
     >              +dx*(12.d0*f(4,3,i,j)+36.d0*dy*f(4,4,i,j))
            enddo
         endif
!
!-----------------------------------
!  access to 3rd derivatives
!
      else if(ict(1).eq.3) then
         if(ict(2).eq.1) then
!  evaluate d3f/dx3 (not continuous)
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >              +6.d0*(
     >         f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j))))
            enddo
         endif
!
         if(ict(3).eq.1) then
!  evaluate d3f/dx2dy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >              2.d0*(
     >           f(3,2,i,j)+dy*(2.d0*f(3,3,i,j)+dy*3.d0*f(3,4,i,j)))
     >              +6.d0*dx*(
     >           f(4,2,i,j)+dy*(2.d0*f(4,3,i,j)+dy*3.d0*f(4,4,i,j)))
            enddo
         endif
!
         if(ict(4).eq.1) then
!  evaluate d3f/dxdy2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >              (2.d0*f(2,3,i,j)+6.d0*dy*f(2,4,i,j)
     >              +2.d0*dx*(2.d0*f(3,3,i,j)+6.d0*dy*f(3,4,i,j)
     >              +1.5d0*dx*(2.d0*f(4,3,i,j)+6.d0*dy*f(4,4,i,j))
     >              ))
            enddo
         endif

         if(ict(5).eq.1) then
!  evaluate d3f/dy3 (not continuous)
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               fval(v,iaval)=6.d0*(f(1,4,i,j)+
     >              dx*(f(2,4,i,j)+dx*(f(3,4,i,j)+dx*f(4,4,i,j))))
            enddo
         endif
!
!-----------------------------------
!  access to 4th derivatives
!
      else if(ict(1).eq.4) then
         if(ict(2).eq.1) then
!  evaluate d4f/dx3dy (not continuous)
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >              +6.d0*(
     >         f(4,2,i,j)+dy*2.d0*(f(4,3,i,j)+dy*1.5d0*f(4,4,i,j)))
            enddo
         endif
!
         if(ict(3).eq.1) then
!  evaluate d4f/dx2dy2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >              4.d0*f(3,3,i,j)+12.d0*dy*f(3,4,i,j)
     >              +dx*(12.d0*f(4,3,i,j)+36.d0*dy*f(4,4,i,j))
            enddo
         endif
!
         if(ict(4).eq.1) then
!  evaluate d4f/dxdy3 (not continuous)
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               fval(v,iaval)=
     >              6.d0*(f(2,4,i,j)
     >              +2.d0*dx*(f(3,4,i,j)+1.5d0*dx*f(4,4,i,j)))
            enddo
         endif
!
!-----------------------------------
!  access to 5th derivatives
!
      else if(ict(1).eq.5) then
         if(ict(2).eq.1) then
!  evaluate d5f/dx3dy2 (not continuous)
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dy=dyv(v)
               fval(v,iaval)=
     >              +12.d0*(f(4,3,i,j)+dy*3.d0*f(4,4,i,j))
            enddo
         endif
!
         if(ict(3).eq.1) then
!  evaluate d5f/dx3dy2 (not continuous)
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               dx=dxv(v)
               fval(v,iaval)=
     >              12.d0*(f(3,4,i,j)+dx*3.d0*f(4,4,i,j))
            enddo
         endif
!
!-----------------------------------
!  access to 6th derivatives
!
      else if(ict(1).eq.6) then
!  evaluate d6f/dx3dy3 (not continuous)
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            fval(v,iaval)=
     >              36.d0*f(4,4,i,j)
         enddo
      endif
!
      return
      end subroutine bcspevfn_db

!----------------------



!  bcspline -- dmc 30 May 1996
!
!  set up coefficients for bicubic spline with following BC's:
!  FULL BC CONTROL at all bdys
!
!  inhomogeneous explicit BCs -- this means setting of 1st or 2nd 
!  derivative at boundary to a non-zero value.
!
!  periodic, not-a-knot, zero derivative, and divided-difference based
!  BCs are "homogeneous"-- i.e. if splines s & t satisfy the BC then
!  the spline (c*s + t) formed as a linear combination of these two
!  splines, also satisfies the BC.
!
!  algorithm note -- handling of inhomogeneous explicit BC's while 
!  maintaining full C2 differentiability is delicate.  Basic method:  use 
!  a fully C2 method based on the "not-a-knot" BC, and then, correct to 
!  meet each user BC by calculating a C2 spline that is zero at all grid
!  points but satisfies a BC which is the difference btw the user spec
!  and the not-a-knot result; add the coeffs of this into the original.
!
!  for this more workspace is needed: nwk .ge. 4*inx*inth +5*max(inx,inth)
!
      subroutine bcspline_db(x,inx,th,inth,fspl,inf3,
     >                    ibcxmin,bcxmin,ibcxmax,bcxmax,
     >                    ibcthmin,bcthmin,ibcthmax,bcthmax,
     >                    wk,nwk,ilinx,ilinth,ier)
!
      implicit none
      integer inx, inth, inf3, nwk, ibcxmin, ibcxmax, ibcthmin, ibcthmax, ilinx,ilinth,ier
      real(dp) x(inx),th(inth),fspl(4,4,inf3,inth),wk(nwk)
      real(dp) bcxmin(inth),bcxmax(inth)
      real(dp) bcthmin(inx),bcthmax(inx)
!
!  input:
!    x(1...inx) -- abscissae, first dimension of data
!   th(1...inth) -- abscissae, second dimension of data  f(x,th)
!   fspl(1,1,1..inx,1..inth) -- function values
!   inf3 -- fspl dimensioning, inf3.ge.inx required.
!
!  boundary conditions input:
!   ibcxmin -- indicator for boundary condition at x(1):
!    bcxmin(...) -- boundary condition data
!     =-1 -- periodic boundary condition
!     =0 -- use "not a knot", bcxmin(...) ignored
!     =1 -- match slope, specified at x(1),th(ith) by bcxmin(ith)
!     =2 -- match 2nd derivative, specified at x(1),th(ith) by bcxmin(ith)
!     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all th(j)
!     =4 -- boundary condition is d2f/dx2=0 at x(1), all th(j)
!     =5 -- match 1st derivative df/dx to 1st divided difference
!     =6 -- match 2nd derivative d2f/dx2 to 2nd divided difference
!     =7 -- match 3rd derivative d3f/dx3 3rd divided difference
!           (for more detailed definition of BCs 5-7, see the
!           comments of subroutine mkspline)
!   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
!
!   ibcxmax -- indicator for boundary condition at x(nx):
!    bcxmax(...) -- boundary condition data
!     (interpretation as with ibcxmin, bcxmin)
!   NOTE:  if ibcxmin=-1, ibcxmax is ignored! ...and the BC is periodic.
!
!   ibcthmin -- indicator for boundary condition at th(1):
!    bcthmin(...) -- boundary condition data
!     (interpretation as with ibcxmin, bcxmin)
!   ibcthmax -- indicator for boundary condition at th(inth):
!    bcthmax(...) -- boundary condition data
!     (interpretation as with ibcxmin, bcxmin)
!   NOTE:  if ibcthmin=-1, ibcthmax is ignored! ...and the BC is periodic.
!
!   NOTE the bcxmin,bcxmax,bcthmin,bcthmax arrays are only used if the
!     corresponding boundary condition flags are set to 1 or 2.
!     Carefully note the dimensioning of these arrays!
!
!  output:
!   fspl(*,*,1..inx,1..inth) -- bicubic spline coeffs (4x4)
!   ...fspl(1,1,*,*) is not replaced.
!
!   ilinx -- =1 on output if x(inx) pts are nearly evenly spaced (tol=1e-3)
!   ilinth-- =1 on output if th(inth) evenly spaced (tol=1e-3)
!
!   ier -- completion code, 0 for normal
!
!  workspace:
!   wk -- must be at least 5*max(inx,inth) large
!                          5*max(inx,inth) + 4*inx*inth large
!                          if explicit non-zero d/dth or d2/dth2 BC's
!                          are supplied.
!  nwk -- size of workspace of workspace provided
!
!---------------------------------
!  in what follows, "f" is an abbreviation for "fspl"...
!
!  compute bicubic spline of 2d function, given values at the grid
!  grid crossing points, f(1,1,i,j)=f(x(i),th(j)).
!
!  on evaluation:  for point x btw x(i) and x(i+1), dx=x-x(i)
!                       and th btw th(j) and th(j+1), dt=th-th(j),
!
!      spline =
!        f(1,1,i,j) + dx*f(2,1,i,j) + dx**2*f(3,1,i,j) + dx**3*f(4,1,i,j)
!   +dt*(f(1,2,i,j) + dx*f(2,2,i,j) + dx**2*f(3,2,i,j) + dx**3*f(4,2,i,j))
!   +d2*(f(1,3,i,j) + dx*f(2,3,i,j) + dx**2*f(3,3,i,j) + dx**3*f(4,3,i,j))
!   +d3*(f(1,4,i,j) + dx*f(2,4,i,j) + dx**2*f(3,4,i,j) + dx**3*f(4,4,i,j))
!
!      where d2=dt**2 and d3=dt**3.
!
      integer iselect1(10)
      integer iselect2(10)
      
      
      integer iflg2, ix, itest, ierx, inxo, ith, jth, ii, iadr, ia5w, iaspl, ierth, intho, ic
      integer ibcthmina, ibcthmaxa, iasc, iinc, iawk, jx
      real(dp) xo2, xo6, zcur, zdiff1, zhxn, zhth, zdiff2
      real(dp) fval(6)
!
!---------------------------------
!
!  see if 2nd pass is needed due to "non-linear" d/dth bdy cond.
!
      iflg2=0
      if(ibcthmin.ne.-1) then
         if((ibcthmin.eq.1).or.(ibcthmin.eq.2)) then
            do ix=1,inx
               if (bcthmin(ix).ne.0.d0) iflg2=1
            enddo
         endif
         if((ibcthmax.eq.1).or.(ibcthmax.eq.2)) then
            do ix=1,inx
               if (bcthmax(ix).ne.0.d0) iflg2=1
            enddo
         endif
      endif
!
      ier=0
      itest=5*max(inx,inth)
      if(iflg2.eq.1) then
         itest=itest +4*inx*inth
      endif
!
      if(nwk.lt.itest) then
         write(6,9901) nwk,itest
 9901    format(' ?bcspline:  workspace too small:'/
     >          '  user supplied:  nwk=',i6,'; need at least:  ',i6/
     >          '  nwk=4*nx*ny +5*max(nx,ny) will work for any user'/
     >          '  choice of bdy conditions.')
         ier=1
      endif
      if(inx.lt.4) then
         write(6,'('' ?bcspline:  at least 4 x points required.'')')
         ier=1
      endif
      if(inth.lt.4) then
         write(6,'('' ?bcspline:  need at least 4 theta points.'')')
         ier=1
      endif
!
      call ibc_ck_db(ibcxmin,'bcspline','xmin',-1,7,ier)
      if(ibcxmin.ge.0) call ibc_ck_db(ibcxmax,'bcspline','xmax',0,7,ier)
      call ibc_ck_db(ibcthmin,'bcspline','thmin',-1,7,ier)
      if(ibcthmin.ge.0) call ibc_ck_db(ibcthmax,'bcspline','thmax',0,7,ier)
!
!  check ilinx & x vector
!
      call splinck_db(x,inx,ilinx,1.0d-3,ierx)
      if(ierx.ne.0) ier=2
!
      if(ier.eq.2) then
         write(6,'('' ?bcspline:  x axis not strict ascending'')')
      endif
!
!  check ilinth & th vector
!
      call splinck_db(th,inth,ilinth,1.0d-3,ierth)
      if(ierth.ne.0) ier=3
!
!      if(ier.eq.3) then
!         write(6,'('' ?bcspline:  th axis not strict ascending'')')
!      endif
!
      if(ier.ne.0) return
!
!------------------------------------
!
      xo2=0.5d0
      xo6=1.d0/6.d0
!
!  spline in x direction
!
      inxo=4*(inx-1)
      do ith=1,inth
!
!  copy the function in
!
         do ix=1,inx
            wk(4*(ix-1)+1)=fspl(1,1,ix,ith)
         enddo
!
         if(ibcxmin.eq.1) then
            wk(2)=bcxmin(ith)
         else if(ibcxmin.eq.2) then
            wk(3)=bcxmin(ith)
         endif
!
         if(ibcxmax.eq.1) then
            wk(inxo+2)=bcxmax(ith)
         else if(ibcxmax.eq.2) then
            wk(inxo+3)=bcxmax(ith)
         endif
!
!  use Wayne's routine
!
         call v_spline_db(ibcxmin,ibcxmax,inx,x,wk,wk(4*inx+1))
!
!  copy the coefficients out
!
         do ix=1,inx
            fspl(2,1,ix,ith)=wk(4*(ix-1)+2)
            fspl(3,1,ix,ith)=wk(4*(ix-1)+3)*xo2
            fspl(4,1,ix,ith)=wk(4*(ix-1)+4)*xo6
         enddo
!
      enddo
!
!-----------------------------------
!
!  spline in theta direction
!
      intho=4*(inth-1)
      do ix=1,inx
!
!  spline each x coeff
!
         do ic=1,4
!
!  copy ordinates in
!
            do ith=1,inth
               wk(4*(ith-1)+1)=fspl(ic,1,ix,ith)
            enddo
!
!  first pass:  use a linear BC -- if flag indicates BC correction
!  will be needed, it will be done later
!
            wk(2)=0.d0
            wk(3)=0.d0
            wk(intho+2)=0.d0
            wk(intho+3)=0.d0
!
            ibcthmina=ibcthmin
            ibcthmaxa=ibcthmax
            if(iflg2.eq.1) then
               if((ibcthmin.eq.1).or.(ibcthmin.eq.2)) ibcthmina=0
               if((ibcthmax.eq.1).or.(ibcthmax.eq.2)) ibcthmaxa=0
            endif
!
            call v_spline_db(ibcthmina,ibcthmaxa,inth,th,wk,wk(4*inth+1))
!
!  copy coeffs out
!
            do ith=1,inth
               fspl(ic,2,ix,ith)=wk(4*(ith-1)+2)
               fspl(ic,3,ix,ith)=wk(4*(ith-1)+3)*xo2
               fspl(ic,4,ix,ith)=wk(4*(ith-1)+4)*xo6
            enddo
!
         enddo
!
      enddo
!
!  now make correction for user BC's if needed
!
      if(iflg2.eq.1) then
!
         iasc=1                         ! wk addr for correction splines
         iinc=4*inth                    ! spacing btw correction splines
         iawk=iasc+4*inth*inx
!
!  last grid zone widths
!
         zhxn=x(inx)-x(inx-1)
         jx=inx-1
         zhth=th(inth)-th(inth-1)
         jth=inth-1
!
         do ii=1,10
            iselect1(ii)=0
            iselect2(ii)=0
         enddo
         if(ibcthmin.eq.1) iselect1(3)=1
         if(ibcthmin.eq.2) iselect1(5)=1
         if(ibcthmax.eq.1) iselect2(3)=1
         if(ibcthmax.eq.2) iselect2(5)=1
!
!  loop over BC's
!
         do ix=1,inx
!
!  (a) d/dth @ th(1) difference btw current BC and user request
!
            if(ibcthmin.eq.1) then
               if(ix.lt.inx) then
                  zcur=fspl(1,2,ix,1)   ! 1st deriv.
               else
                  zcur=fspl(1,2,jx,1)+zhxn*(fspl(2,2,jx,1)+zhxn*
     >               (fspl(3,2,jx,1)+zhxn*fspl(4,2,jx,1)))
               endif
               zdiff1=bcthmin(ix)-zcur
            else if(ibcthmin.eq.2) then
               if(ix.lt.inx) then
                  zcur=2.d0*fspl(1,3,ix,1) ! 2nd deriv.
               else
                  zcur=2.d0*(fspl(1,3,jx,1)+zhxn*(fspl(2,3,jx,1)+zhxn*
     >               (fspl(3,3,jx,1)+zhxn*fspl(4,3,jx,1))))
               endif
               zdiff1=bcthmin(ix)-zcur
            else
               zdiff1=0.d0
            endif
!
!  (b) d/dth @ th(inth) difference btw current BC and user request
!
            if(ibcthmax.eq.1) then
               if(ix.lt.inx) then
!  1st deriv.
                  zcur=fspl(1,2,ix,jth)+zhth*(2.d0*fspl(1,3,ix,jth)+
     >               zhth*3.d0*fspl(1,4,ix,jth))
               else
                  call bcspeval_db(x(inx),th(inth),iselect2,fval,
     >               x,inx,th,inth,ilinx,ilinth,fspl,inf3,ier)
                  zcur=fval(1)
                  if(ier.ne.0) return
               endif
               zdiff2=bcthmax(ix)-zcur
            else if(ibcthmax.eq.2) then
               if(ix.lt.inx) then
!  2nd deriv.
                  zcur=2.d0*fspl(1,3,ix,jth)+
     >               6.d0*zhth*fspl(1,4,ix,jth)
               else
                  call bcspeval_db(x(inx),th(inth),iselect2,fval,
     >               x,inx,th,inth,ilinx,ilinth,fspl,inf3,ier)
                  zcur=fval(1)
                  if(ier.ne.0) return
               endif
               zdiff2=bcthmax(ix)-zcur
            else
               zdiff2=0.d0
            endif
!
!  ok compute the theta spline with BC's to span the difference(s)
!  these theta "correction splines" are zero at all the grid points
!  but have at least one non-zero 1st or 2nd derivative BC
!
            iadr=iasc+(ix-1)*iinc
            do ith=1,inth
               wk(iadr+4*(ith-1))=0.d0
            enddo
!
            wk(iadr+1)=0.d0
            wk(iadr+2)=0.d0
            wk(iadr+intho+1)=0.d0
            wk(iadr+intho+2)=0.d0
!
            if(ibcthmin.eq.1) then
               wk(iadr+1)=zdiff1
            else if(ibcthmin.eq.2) then
               wk(iadr+2)=zdiff1
            endif
!
            if(ibcthmax.eq.1) then
               wk(iadr+intho+1)=zdiff2
            else if(ibcthmax.eq.2) then
               wk(iadr+intho+2)=zdiff2
            endif
!
            call v_spline_db(ibcthmin,ibcthmax,inth,th,wk(iadr),wk(iawk))
         enddo
!
!  add in results to main array -- th spline coef corrections
!
         do ix=1,inx
            iadr=iasc+(ix-1)*iinc
            do ith=1,inth-1
               wk(iadr+4*(ith-1)+2)=wk(iadr+4*(ith-1)+2)*xo2
               wk(iadr+4*(ith-1)+3)=wk(iadr+4*(ith-1)+3)*xo6
               if(ix.lt.inx) then
                  fspl(1,2,ix,ith)=fspl(1,2,ix,ith)+wk(iadr+4*(ith-1)+1)
                  fspl(1,3,ix,ith)=fspl(1,3,ix,ith)+wk(iadr+4*(ith-1)+2)
                  fspl(1,4,ix,ith)=fspl(1,4,ix,ith)+wk(iadr+4*(ith-1)+3)
               endif
            enddo
         enddo
!
!  compute the x splines of the th spline correction coeffs
!
         ia5w=iawk+4*inx
!
         do ith=1,inth-1
            do ic=2,4
               do ix=1,inx
                  iaspl=iasc+iinc*(ix-1)
                  wk(iawk+4*(ix-1))=wk(iaspl+4*(ith-1)+(ic-1))
               enddo
!
!  use zero BCs for this correction spline
!
               wk(iawk+1)=0.d0
               wk(iawk+2)=0.d0
               wk(iawk+inxo+1)=0.d0
               wk(iawk+inxo+2)=0.d0
!
!  periodic spline of correction spline higher coeffs (1st coeffs are
!  all zero by defn of the correction spline
!
               call v_spline_db(ibcxmin,ibcxmax,inx,x,wk(iawk),wk(ia5w))
!
               do ix=1,inx-1
                  fspl(2,ic,ix,ith)=fspl(2,ic,ix,ith)+
     >               wk(iawk+4*(ix-1)+1)
                  fspl(3,ic,ix,ith)=fspl(3,ic,ix,ith)+
     >               wk(iawk+4*(ix-1)+2)*xo2
                  fspl(4,ic,ix,ith)=fspl(4,ic,ix,ith)+
     >               wk(iawk+4*(ix-1)+3)*xo6
               enddo
!
            enddo
         enddo                          ! ith
!
      endif                             ! BC correction needs test
!
      return
      end subroutine bcspline_db



!  cspline -- dmc 15 Feb 1999
!
!  a standard interface to the 1d spline setup routine
!    modified dmc 3 Mar 2000 -- to use Wayne Houlberg's v_spline code.
!    new BC options added.
!
      subroutine cspline_db(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   wk,iwk,ilinx,ier)
!
      implicit none
      integer nx, iwk
      real(dp) x(nx)                        ! x axis (in)
      real(dp) fspl(4,nx)                   ! spline data (in/out)
      integer ibcxmin                   ! x(1) BC flag (in, see comments)
      real(dp) bcxmin                       ! x(1) BC data (in, see comments)
      integer ibcxmax                   ! x(nx) BC flag (in, see comments)
      real(dp) bcxmax                       ! x(nx) BC data (in, see comments)
      real(dp) wk(iwk)                      ! workspace of size at least nx
      integer ilinx                     ! even spacing flag (out)
      integer ier                       ! output, =0 means OK
!
!  ** note wk(...) array is not used unless ibcxmin=-1 (periodic spline
!  evaluation)
!
!  this routine computes spline coefficients for a 1d spline --
!  evaluation of the spline can be done by cspeval.for subroutines
!  or directly by inline code.
!
!  the input x axis x(1...nx) must be strictly ascending, i.e.
!  x(i+1).gt.x(i) is required for i=1 to nx-1.  This is checked and
!  ier=1 is set and the routine exits if the test is not satisfied.
!
!  on output, ilinx=1 is set if, to a reasonably close tolerance,
!  all grid spacings x(i+1)-x(i) are equal.  This allows a speedier
!  grid lookup algorithm on evaluation of the spline.  If on output
!  ilinx=2, this means the spline x axis is not evenly spaced.
!
!  the input data for the spline are given in f[j] = fspl(1,j).  The
!  output data are the spline coefficients fspl(2,j),fspl(3,j), and
!  fspl(4,j), j=1 to nx.  The result is a spline s(x) satisfying the
!  boundary conditions and with the properties
!
!     s(x(j)) = fspl(1,j)
!     s'(x) is continuous even at the grid points x(j)
!     s''(x) is continuous even at the grid points x(j)
!
!  the formula for evaluation of s(x) is:
!
!     let dx = x-x(i), where x(i).le.x.le.x(i+1).  Then,
!     s(x)=fspl(1,i) + dx*(fspl(2,i) +dx*(fspl(3,i) + dx*fspl(4,i)))
!
!  ==>boundary conditions.  Complete specification of a 1d spline
!  requires specification of boundary conditions at x(1) and x(nx).
!
!  this routine provides 4 options:
!
! -1 ***** PERIODIC BC
!  ibcxmin=-1  --  periodic boundary condition.  This means the
!    boundary conditions s'(x(1))=s'(x(nx)) and s''(x(1))=s''(x(nx))
!    are imposed.  Note that s(x(1))=s(x(nx)) (i.e. fspl(1,1)=fspl(1,nx))
!    is not required -- that is determined by the fspl array input data.
!    The periodic boundary condition is to be preferred for periodic
!    data.  When splining periodic data f(x) with period P, the relation
!    x(nx)=x(1)+n*P, n = the number of periods (usually 1), should hold.
!    (ibcxmax, bcxmin, bcxmax are ignored).
!
!  if a periodic boundary condition is set, this covers both boundaries.
!  for the other types of boundary conditions, the type of condition
!  chosen for the x(1) boundary need not be the same as the type chosen
!  for the x(nx) boundary.
!
!  0 ***** NOT A KNOT BC
!  ibcxmin=0 | ibcxmax=0 -- this specifies a "not a knot" boundary
!    condition -- see cubsplb.for.  This is a common way for inferring
!    a "good" spline boundary condition automatically from data in the
!    vicinity of the boundary.  (bcxmin | bcxmax are ignored).
!
!  1 ***** BC:  SPECIFIED SLOPE
!  ibcxmin=1 | ibcxmax=1 -- boundary condition is to have s'(x(1)) |
!    s'(x(nx)) match the passed value (bcxmin | bcxmax).
!
!  2 ***** BC:  SPECIFIED 2nd DERIVATIVE
!  ibcxmin=2 | ibcxmax=2 -- boundary condition is to have s''(x(1)) |
!    s''(x(nx)) match the passed value (bcxmin | bcxmax).
!
!  3 ***** BC:  SPECIFIED SLOPE = 0.0
!  ibcxmin=3 | ibcxmax=3 -- boundary condition is to have s'(x(1)) |
!    s'(x(nx)) equal to ZERO.
!
!  4 ***** BC:  SPECIFIED 2nd DERIVATIVE = 0.0
!  ibcxmin=4 | ibcxmax=4 -- boundary condition is to have s''(x(1)) |
!    s''(x(nx)) equal to ZERO.
!
!  5 ***** BC:  1st DIVIDED DIFFERENCE
!  ibcxmin=5 | ibcxmax=5 -- boundary condition is to have s'(x(1)) |
!    s'(x(nx)) equal to the slope from the 1st|last 2 points
!
!  6 ***** BC:  2nd DIVIDED DIFFERENCE
!  ibcxmin=6 | ibcxmax=6 -- boundary condition is to have s''(x(1)) |
!    s''(x(nx)) equal to the 2nd derivative from the 1st|last 3 points
!
!  7 ***** BC:  3rd DIVIDED DIFFERENCE
!  ibcxmin=7 | ibcxmax=7 -- boundary condition is to have s'''(x(1)) |
!    s'''(x(nx)) equal to the 3rd derivative from the 1st|last 4 points
!
!---------------------------------------------------------------------
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: sixth = 1d0/6d0
      integer inum,i,ierx
!
!  error checks
!
      ier = 0
      if(nx.lt.4) then
!         write(6,'('' ?cspline:  at least 4 x points required.'')')
         ier=1
      endif
      call ibc_ck_db(ibcxmin,'cspline','xmin',-1,7,ier)
      if(ibcxmin.ge.0) call ibc_ck_db(ibcxmax,'cspline','xmax',0,7,ier)
!
!  x axis check
!
      call splinck_db(x,nx,ilinx,1.0d-3,ierx)
      if(ierx.ne.0) ier=2
!
!      if(ier.eq.2) then
!         write(6,'('' ?cspline:  x axis not strict ascending'')')
!      endif
!
      if(ibcxmin.eq.-1) then
         inum=nx
         if(iwk.lt.inum) then
!            write(6,1009) inum,iwk,nx
! 1009       format(
!     >      ' ?cspline:  workspace too small.  need:  ',i6,' got:  ',i6/
!     >      '  (need = nx, nx=',i6)
            ier=3
         endif
      endif
!
      if(ier.ne.0) return
!
!  OK -- evaluate spline
!
      if(ibcxmin.eq.1) then
         fspl(2,1)=bcxmin
      else if(ibcxmin.eq.2) then
         fspl(3,1)=bcxmin
      endif
!
      if(ibcxmax.eq.1) then
         fspl(2,nx)=bcxmax
      else if(ibcxmax.eq.2) then
         fspl(3,nx)=bcxmax
      endif
!
      call v_spline_db(ibcxmin,ibcxmax,nx,x,fspl,wk)
!
      do i=1,nx
         fspl(3,i)=half*fspl(3,i)
         fspl(4,i)=sixth*fspl(4,i)
      enddo
!
      return
      end subroutine cspline_db


      subroutine evbicub_db(xget,yget,x,nx,y,ny,ilinx,iliny,
     >                   f1,inf2,ict,fval,ier)
C
C  use mkbicub to set up spline coefficients!
C
C  evaluate a 2d cubic Spline interpolant on a rectilinear
C  grid -- this is C2 in both directions.
C
C  this subroutine calls two subroutines:
C     herm2xy  -- find cell containing (xget,yget)
C     fvbicub  -- evaluate interpolant function and (optionally) derivatives
C
C  input arguments:
C  ================
C
      implicit none
      integer nx,ny                     ! grid sizes
      real(dp) xget,yget                    ! target of this interpolation
      real(dp) x(nx)                        ! ordered x grid
      real(dp) y(ny)                        ! ordered y grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
      integer iliny                     ! iliny=1 => assume y evenly spaced
C
      integer inf2
      real(dp), pointer :: f1(:)
C
C       f 2nd dimension inf2 must be .ge. nx
C       contents of f:
C
C  f(0,i,j) = f @ x(i),y(j)
C  f(1,i,j) = d2f/dx2 @ x(i),y(j)
C  f(2,i,j) = d2f/dy2 @ x(i),y(j)
C  f(3,i,j) = d4f/dx2dy2 @ x(i),y(j)
C
C      (these are spline coefficients selected for continuous 2-
C      diffentiability, see mkbicub[w].for)
C
      integer ict(6)                    ! code specifying output desired
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return df/dy  (0, don't)
C  ict(4)=1 -- return d2f/dx2  (0, don't)
C  ict(5)=1 -- return d2f/dy2  (0, don't)
C  ict(6)=1 -- return d2f/dxdy (0, don't)
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
C
C output arguments:
C =================
C
      real(dp) fval(6)                      ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C  fval(3) receives the third output (depends on ict(...) spec)
C  fval(4) receives the fourth output (depends on ict(...) spec)
C  fval(5) receives the fourth output (depends on ict(...) spec)
C  fval(6) receives the fourth output (depends on ict(...) spec)
C
C  examples:
C    on input ict = [1,1,1,0,0,1]
C   on output fval= [f,df/dx,df/dy,d2f/dxdy], elements 5 & 6 not referenced.
C
C    on input ict = [1,0,0,0,0,0]
C   on output fval= [f] ... elements 2 -- 6 never referenced.
C
C    on input ict = [0,0,0,1,1,0]
C   on output fval= [d2f/dx2,d2f/dy2] ... elements 3 -- 6 never referenced.
C
C    on input ict = [0,0,1,0,0,0]
C   on output fval= [df/dy] ... elements 2 -- 6 never referenced.
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i,j                       ! cell indices
C
C  normalized displacement from (x(i),y(j)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C    yparam=0 @y(j)  yparam=1 @y(j+1)
C
      real(dp) xparam,yparam
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
C
      real(dp) hx,hy
      real(dp) hxi,hyi
C
C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1
C
C  ** the interface is very similar to herm2ev.for; can use herm2xy **
C---------------------------------------------------------------------
C
      call herm2xy_db(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i,j,xparam,yparam,hx,hxi,hy,hyi,ier)
      if(ier.ne.0) return
!
      call fvbicub_db(ict,1,1,
     >   fval,(/i/),(/j/),(/xparam/),(/yparam/),
     <   (/hx/),(/hxi/),(/hy/),(/hyi/),f1,inf2,ny)
C
      return
      end subroutine evbicub_db

C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 2d fcn
C   --vectorized-- dmc 10 Feb 1999
C
C  use mkbicub to set up spline coefficients!
C
      subroutine fvbicub_db(ict,ivec,ivecd,
     >   fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >   f1,inf2,ny)
C
      integer ict(6)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj(ivec)         ! target cells (i,j)
      real(dp) xparam(ivec),yparam(ivec)
                          ! normalized displacements from (i,j) corners
C
      real(dp) hx(ivec),hy(ivec)            ! grid spacing, and
      real(dp) hxi(ivec),hyi(ivec)          ! inverse grid spacing 1/(x(i+1)-x(i))
                                        ! & 1/(y(j+1)-y(j))
C
      real(dp), pointer :: f1(:)
C
      real(dp) fval(ivecd,6)                ! output returned
C
C  for detailed description of fin, ict and fval see subroutine
C  evbicub comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii,jj and parameter inputs
C     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C  Spline evaluation consists of a "mixing" of the interpolant
C  data using the linear functionals xparam, xpi = 1-xparam,
C  yparam, ypi = 1-yparam, and the cubic functionals
C  xparam**3-xparam, xpi**3-xpi, yparam**3-yparam, ypi**3-ypi ...
C  and their derivatives as needed.
C
      integer v
      real(dp) sum
C
      real(dp), parameter :: sixth = 1.d0/6.d0

      real(dp), pointer :: fin(:,:,:)
      fin(0:3,1:inf2,1:ny) => f1(1:4*inf2*ny)


C
C---------------
C   ...in x direction
C
      z36th=sixth*sixth
      iadr=0
C
      if(ict(1).le.2) then
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
C  in x direction...
C
               xp=xparam(v)
               xpi=1.d0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.d0)
               cxi=xpi*(xpi2-1.d0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.d0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.d0)
               cyi=ypi*(ypi2-1.d0)
               hy2=hy(v)*hy(v)
C
               sum=xpi*(ypi*fin(0,i,j)  +yp*fin(0,i,j+1))+
     >              xp*(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1))
C
               sum=sum+sixth*hx2*(
     >              cxi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+
     >              cx*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
C
               sum=sum+sixth*hy2*(
     >              xpi*(cyi*fin(2,i,j)  +cy*fin(2,i,j+1))+
     >              xp*(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))
C
               sum=sum+z36th*hx2*hy2*(
     >              cxi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+
     >              cx*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
C  in x direction...
C
               xp=xparam(v)
               xpi=1.d0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.d0*xp2-1.d0
               cxdi=-3.d0*xpi2+1.d0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.d0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.d0)
               cyi=ypi*(ypi2-1.d0)
               hy2=hy(v)*hy(v)
C
               sum=hxi(v)*(
     >              -(ypi*fin(0,i,j)  +yp*fin(0,i,j+1))
     >              +(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1)))
C
               sum=sum+sixth*hx(v)*(
     >              cxdi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+
     >              cxd*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
C
               sum=sum+sixth*hxi(v)*hy2*(
     >              -(cyi*fin(2,i,j)  +cy*fin(2,i,j+1))
     >              +(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))
C
               sum=sum+z36th*hx(v)*hy2*(
     >              cxdi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+
     >              cxd*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
C  in x direction...
C
               xp=xparam(v)
               xpi=1.d0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.d0)
               cxi=xpi*(xpi2-1.d0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.d0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.d0*yp2-1.d0
               cydi=-3.d0*ypi2+1.d0
C
               sum=hyi(v)*(
     >              xpi*(-fin(0,i,j)  +fin(0,i,j+1))+
     >              xp*(-fin(0,i+1,j)+fin(0,i+1,j+1)))
C
               sum=sum+sixth*hx2*hyi(v)*(
     >              cxi*(-fin(1,i,j)  +fin(1,i,j+1))+
     >              cx*(-fin(1,i+1,j)+fin(1,i+1,j+1)))
C
               sum=sum+sixth*hy(v)*(
     >              xpi*(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1))+
     >              xp*(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))
C
               sum=sum+z36th*hx2*hy(v)*(
     >              cxi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+
     >              cx*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C
C  d2f/dx2:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
C  in x direction...
C
               xp=xparam(v)
               xpi=1.d0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.d0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.d0)
               cyi=ypi*(ypi2-1.d0)
               hy2=hy(v)*hy(v)
C
               sum=(
     >              xpi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+
     >              xp*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
C
               sum=sum+sixth*hy2*(
     >              xpi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+
     >              xp*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C
C  d2f/dy2:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
C  in x direction...
C
               xp=xparam(v)
               xpi=1.d0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.d0)
               cxi=xpi*(xpi2-1.d0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.d0-yp
C
               sum=(
     >              xpi*(ypi*fin(2,i,j)  +yp*fin(2,i,j+1))+
     >              xp*(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))
C
               sum=sum+sixth*hx2*(
     >              cxi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+
     >              cx*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C
C  d2f/dxdy:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
C  in x direction...
C
               xp=xparam(v)
               xpi=1.d0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.d0*xp2-1.d0
               cxdi=-3.d0*xpi2+1.d0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.d0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.d0*yp2-1.d0
               cydi=-3.d0*ypi2+1.d0
C
               sum=hxi(v)*hyi(v)*(
     >              fin(0,i,j)  -fin(0,i,j+1)
     >              -fin(0,i+1,j)+fin(0,i+1,j+1))
C
               sum=sum+sixth*hx(v)*hyi(v)*(
     >              cxdi*(-fin(1,i,j)  +fin(1,i,j+1))+
     >              cxd*(-fin(1,i+1,j)+fin(1,i+1,j+1)))
C
               sum=sum+sixth*hxi(v)*hy(v)*(
     >              -(cydi*fin(2,i,j)  +cyd*fin(2,i,j+1))
     >              +(cydi*fin(2,i+1,j)+cyd*fin(2,i+1,j+1)))
C
               sum=sum+z36th*hx(v)*hy(v)*(
     >              cxdi*(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))+
     >              cxd*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C-------------------------------------------------
C
      else if(ict(1).eq.3) then
         if(ict(2).eq.1) then
!  evaluate d3f/dx3 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               yp=yparam(v)
               ypi=1.d0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
               cy=yp*(yp2-1.d0)
               cyi=ypi*(ypi2-1.d0)
               hy2=hy(v)*hy(v)
               sum=hxi(v)*(
     >              -(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))
     >              +(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
C
               sum=sum+sixth*hy2*hxi(v)*(
     >              -(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))
     >              +(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
!
         if(ict(3).eq.1) then
!  evaluate d3f/dx2dy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               xp=xparam(v)
               xpi=1.d0-xp
               yp=yparam(v)
               ypi=1.d0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
               cyd=3.d0*yp2-1.d0
               cydi=-3.d0*ypi2+1.d0
C
               sum=hyi(v)*(
     >              xpi*(-fin(1,i,j)  +fin(1,i,j+1))+
     >              xp*(-fin(1,i+1,j) +fin(1,i+1,j+1)))
C
               sum=sum+sixth*hy(v)*(
     >              xpi*(cydi*fin(3,i,j) +cyd*fin(3,i,j+1))+
     >              xp*(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
!
         if(ict(4).eq.1) then
!  evaluate d3f/dxdy2
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               xp=xparam(v)
               xpi=1.d0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
               cxd=3.d0*xp2-1.d0
               cxdi=-3.d0*xpi2+1.d0
               yp=yparam(v)
               ypi=1.d0-yp
C
               sum=hxi(v)*(
     >              -(ypi*fin(2,i,j)  +yp*fin(2,i,j+1))
     >              +(ypi*fin(2,i+1,j)+yp*fin(2,i+1,j+1)))
C
               sum=sum+sixth*hx(v)*(
     >              cxdi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+
     >              cxd*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif

         if(ict(5).eq.1) then
!  evaluate d3f/dy3 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
               xp=xparam(v)
               xpi=1.d0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.d0)
               cxi=xpi*(xpi2-1.d0)
               hx2=hx(v)*hx(v)
C
               sum=hyi(v)*(
     >              xpi*(-fin(2,i,j)  +fin(2,i,j+1))+
     >              xp*(-fin(2,i+1,j) +fin(2,i+1,j+1)))
C
               sum=sum+sixth*hx2*hyi(v)*(
     >              cxi*(-fin(3,i,j)  +fin(3,i,j+1))+
     >              cx*(-fin(3,i+1,j) +fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
!
!-----------------------------------
!  access to 4th derivatives
!
      else if(ict(1).eq.4) then
         if(ict(2).eq.1) then
!  evaluate d4f/dx3dy (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               yp=yparam(v)
               ypi=1.d0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
               cyd=3.d0*yp2-1.d0
               cydi=-3.d0*ypi2+1.d0
C
               sum=hxi(v)*hyi(v)*(
     >              +( fin(1,i,j)  -fin(1,i,j+1))
     >              +(-fin(1,i+1,j)+fin(1,i+1,j+1)))
C
               sum=sum+sixth*hy(v)*hxi(v)*(
     >              -(cydi*fin(3,i,j)  +cyd*fin(3,i,j+1))
     >              +(cydi*fin(3,i+1,j)+cyd*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
!
         if(ict(3).eq.1) then
!  evaluate d4f/dx2dy2
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
               xp=xparam(v)
               xpi=1.d0-xp
               yp=yparam(v)
               ypi=1.d0-yp
C
               sum=xpi*(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))+
     >              xp*(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1))
C
               fval(v,iadr)=sum
            enddo
         endif
!
         if(ict(4).eq.1) then
!  evaluate d4f/dxdy3 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
               xp=xparam(v)
               xpi=1.d0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.d0*xp2-1.d0
               cxdi=-3.d0*xpi2+1.d0
C
               sum=hyi(v)*hxi(v)*(
     >              +( fin(2,i,j)  -fin(2,i,j+1))
     >              +(-fin(2,i+1,j)+fin(2,i+1,j+1)))
C
               sum=sum+sixth*hx(v)*hyi(v)*(
     >              cxdi*(-fin(3,i,j)  +fin(3,i,j+1))+
     >              cxd*(-fin(3,i+1,j) +fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
!
!-----------------------------------
!  access to 5th derivatives
!
      else if(ict(1).eq.5) then
         if(ict(2).eq.1) then
!  evaluate d5f/dx3dy2 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
               yp=yparam(v)
               ypi=1.d0-yp
C
               sum=hxi(v)*(
     >              -(ypi*fin(3,i,j)  +yp*fin(3,i,j+1))
     >              +(ypi*fin(3,i+1,j)+yp*fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
!
         if(ict(3).eq.1) then
!  evaluate d5f/dx2dy3 (not continuous)
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
C
               xp=xparam(v)
               xpi=1.d0-xp
C
               sum=hyi(v)*(
     >              xpi*(-fin(3,i,j)  +fin(3,i,j+1))+
     >              xp*(-fin(3,i+1,j)+fin(3,i+1,j+1)))
C
               fval(v,iadr)=sum
            enddo
         endif
!
!-----------------------------------
!  access to 6th derivatives
!
      else if(ict(1).eq.6) then
!  evaluate d6f/dx3dy3 (not continuous)
         iadr=iadr+1
         do v=1,ivec
            i=ii(v)
            j=jj(v)
            sum=hxi(v)*hyi(v)*(
     >              +( fin(3,i,j)  -fin(3,i,j+1))
     >              +(-fin(3,i+1,j)+fin(3,i+1,j+1)))
            fval(v,iadr)=sum
         enddo
      endif
!
      return
      end subroutine fvbicub_db


      subroutine herm2ev_db(xget,yget,x,nx,y,ny,ilinx,iliny,
     >                   f,inf2,ict,fval,ier)
!
!  evaluate a 2d cubic Hermite interpolant on a rectilinear
!  grid -- this is C1 in both directions.
!
!  this subroutine calls two subroutines:
!     herm2xy  -- find cell containing (xget,yget)
!     herm2fcn -- evaluate interpolant function and (optionally) derivatives
!
!  input arguments:
!  ================
!
      real(dp) xget,yget                    ! target of this interpolation
      real(dp) x(nx)                        ! ordered x grid
      real(dp) y(ny)                        ! ordered y grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
      integer iliny                     ! iliny=1 => assume y evenly spaced
!
      real(dp) f(0:3,inf2,ny)               ! function data
!
!       f 2nd dimension inf2 must be .ge. nx
!       contents of f:
!
!  f(0,i,j) = f @ x(i),y(j)
!  f(1,i,j) = df/dx @ x(i),y(j)
!  f(2,i,j) = df/dy @ x(i),y(j)
!  f(3,i,j) = d2f/dxdy @ x(i),y(j)
!
      integer ict(4)                    ! code specifying output desired
!
!  ict(1)=1 -- return f  (0, don't)
!  ict(2)=1 -- return df/dx  (0, don't)
!  ict(3)=1 -- return df/dy  (0, don't)
!  ict(4)=1 -- return d2f/dxdy (0, don't)
!
! output arguments:
! =================
!
      real(dp) fval(4)                      ! output data
      integer ier                       ! error code =0 ==> no error
!
!  fval(1) receives the first output (depends on ict(...) spec)
!  fval(2) receives the second output (depends on ict(...) spec)
!  fval(3) receives the third output (depends on ict(...) spec)
!  fval(4) receives the fourth output (depends on ict(...) spec)
!
!  examples:
!    on input ict = [1,1,1,1]
!   on output fval= [f,df/dx,df/dy,d2f/dxdy]
!
!    on input ict = [1,0,0,0]
!   on output fval= [f] ... elements 2 & 3 & 4 never referenced
!
!    on input ict = [0,1,1,0]
!   on output fval= [df/dx,df/dy] ... element 3 & 4 never referenced
!
!    on input ict = [0,0,1,0]
!   on output fval= [df/dy] ... elements 2 & 3 & 4 never referenced.
!
!  ier -- completion code:  0 means OK
!-------------------
!  local:
!
      integer i,j                       ! cell indices
!
!  normalized displacement from (x(i),y(j)) corner of cell.
!    xparam=0 @x(i)  xparam=1 @x(i+1)
!    yparam=0 @y(j)  yparam=1 @y(j+1)
!
      real(dp) xparam,yparam
!
!  cell dimensions and
!  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
!
      real(dp) hx,hy
      real(dp) hxi,hyi
!
!  0 .le. xparam .le. 1
!  0 .le. yparam .le. 1
!
!---------------------------------------------------------------------
!
      call herm2xy_db(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i,j,xparam,yparam,hx,hxi,hy,hyi,ier)
      if(ier.ne.0) return
!
      call herm2fcn_db(ict,1,1,
     >   fval,(/i/),(/j/),(/xparam/),(/yparam/),
     >   (/hx/),(/hxi/),(/hy/),(/hyi/),f,inf2,ny)
!
      return
      end subroutine herm2ev_db

!---------------------------------------------------------------------
!  herm2xy -- look up x-y zone
!
!  this is the "first part" of herm2ev, see comments, above.
!
      subroutine herm2xy_db(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i,j,xparam,yparam,hx,hxi,hy,hyi,ier)
!
!  input of herm2xy
!  ================
      implicit none
!
      integer nx,ny                     ! array dimensions
!
      real(dp) xget,yget                    ! target point
      real(dp) x(:) ! (nx)                  ! indep. coords., strict ascending
      real(dp) y(:) ! (ny)                  ! indep. coords., strict ascending
!
      integer ilinx                     ! =1:  x evenly spaced
      integer iliny                     ! =1:  y evenly spaced
!
!  output of herm2xy
!  =================
      integer i,j                       ! index to cell containing target pt
!          on exit:  1.le.i.le.nx-1   1.le.j.le.ny-1
!
!  normalized position w/in (i,j) cell (btw 0 and 1):
!
      real(dp) xparam                       ! (xget-x(i))/(x(i+1)-x(i))
      real(dp) yparam                       ! (yget-y(j))/(y(j+1)-y(j))
!
!  grid spacing
!
      real(dp) hx                           ! hx = x(i+1)-x(i)
      real(dp) hy                           ! hy = y(j+1)-y(j)
!
!  inverse grid spacing:
!
      real(dp) hxi                          ! 1/hx = 1/(x(i+1)-x(i))
      real(dp) hyi                          ! 1/hy = 1/(y(j+1)-y(j))
!
      integer ier                       ! return ier.ne.0 on error
!
!------------------------------------
      real(dp) zxget,zyget,zxtol,zytol
      integer nxm,nym,ii,jj                  

!
      ier=0
!
!  range check
!
      zxget=xget
      zyget=yget
      if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
         zxtol=4.0d-7*max(abs(x(1)),abs(x(nx)))
         if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
            ier=1
            !write(6,1001) xget,x(1),x(nx)
! 1001       format(' ?herm2ev:  xget=',1pe11.4,' out of range ',
!     >         1pe11.4,' to ',1pe11.4)
         else
            !if((xget.lt.x(1)-0.5*zxtol).or.(xget.gt.x(nx)+0.5*zxtol)) write(6,1011) xget,x(1),x(nx)
! 1011       format(' %herm2ev:  xget=',1pe15.8,' beyond range ',
!     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(xget.lt.x(1)) then
               zxget=x(1)
            else
               zxget=x(nx)
            endif
         endif
      endif
      if((yget.lt.y(1)).or.(yget.gt.y(ny))) then
         zytol=4.0d-7*max(abs(y(1)),abs(y(ny)))
         if((yget.lt.y(1)-zytol).or.(yget.gt.y(ny)+zytol)) then
            ier=1
!            write(6,1002) yget,y(1),y(ny)
! 1002       format(' ?herm2ev:  yget=',1pe11.4,' out of range ',
!     >         1pe11.4,' to ',1pe11.4)
         else
!            if((yget.lt.y(1)-0.5*zytol).or.
!     >         (yget.gt.y(ny)+0.5*zytol))
!     >      write(6,1012) yget,y(1),y(ny)
! 1012       format(' %herm2ev:  yget=',1pe15.8,' beyond range ',
!     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(yget.lt.y(1)) then
               zyget=y(1)
            else
               zyget=y(ny)
            endif
         endif
      endif
      if(ier.ne.0) return
!
!  now find interval in which target point lies..
!
      nxm=nx-1
      nym=ny-1
!
      if(ilinx.eq.1) then
         ii=1+nxm*(zxget-x(1))/(x(nx)-x(1))
         i=min(nxm, ii)
         if(zxget.lt.x(i)) then
            i=i-1
         else if(zxget.gt.x(i+1)) then
            i=i+1
         endif
      else
         if((1.le.i).and.(i.lt.nxm)) then
            if((x(i).le.zxget).and.(zxget.le.x(i+1))) then
               continue  ! already have the zone
            else
               call zonfind_db(x,nx,zxget,i)
            endif
         else
            call zonfind_db(x,nx,zxget,i)
         endif
      endif
!
      if(iliny.eq.1) then
         jj=1+nym*(zyget-y(1))/(y(ny)-y(1))
         j=min(nym, jj)
         if(zyget.lt.y(j)) then
            j=j-1
         else if(zyget.gt.y(j+1)) then
            j=j+1
         endif
      else
         if((1.le.j).and.(j.lt.nym)) then
            if((y(j).le.zyget).and.(zyget.le.y(j+1))) then
               continue  ! already have the zone
            else
               call zonfind_db(y,ny,zyget,j)
            endif
         else
            call zonfind_db(y,ny,zyget,j)
         endif
      endif
!
      hx=(x(i+1)-x(i))
      hy=(y(j+1)-y(j))
!
      hxi=1.d0/hx
      hyi=1.d0/hy
!
      xparam=(zxget-x(i))*hxi
      yparam=(zyget-y(j))*hyi
!
      return
      end subroutine herm2xy_db


!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 2d fcn
!   --vectorized-- dmc 10 Feb 1999
!
      subroutine herm2fcn_db(ict,ivec,ivecd,
     >   fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >   fin,inf2,ny)
!
      integer ict(4)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
!
      integer ii(ivec),jj(ivec)         ! target cells (i,j)
      real(dp) xparam(ivec),yparam(ivec)
                          ! normalized displacements from (i,j) corners
!
      real(dp) hx(ivec),hy(ivec)            ! grid spacing, and
      real(dp) hxi(ivec),hyi(ivec)          ! inverse grid spacing 1/(x(i+1)-x(i))
                                        ! & 1/(y(j+1)-y(j))
!
      real(dp) fin(0:3,inf2,ny)             ! interpolant data (cf "herm2ev")
!
      real(dp) fval(ivecd,4)                ! output returned
!
!  for detailed description of fin, ict and fval see subroutine
!  herm2ev comments.  Note ict is not vectorized; the same output
!  is expected to be returned for all input vector data points.
!
!  note that the index inputs ii,jj and parameter inputs
!     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
!     output array fval has a vector ** 1st dimension ** whose
!     size must be given as a separate argument
!
!  to use this routine in scalar mode, pass in ivec=ivecd=1
!
!---------------
!  Hermite cubic basis functions
!  -->for function value matching
!     a(0)=0 a(1)=1        a'(0)=0 a'(1)=0
!   abar(0)=1 abar(1)=0  abar'(0)=0 abar'(1)=0
!
!   a(x)=-2*x**3 + 3*x**2    = x*x*(-2.d0*x+3.0)
!   abar(x)=1-a(x)
!   a'(x)=-abar'(x)          = 6.d0*x*(1.0-x)
!
!  -->for derivative matching
!     b(0)=0 b(1)=0          b'(0)=0 b'(1)=1
!   bbar(0)=0 bbar(1)=0  bbar'(0)=1 bbar'(1)=0
!
!     b(x)=x**3-x**2         b'(x)=3*x**2-2*x
!     bbar(x)=x**3-2*x**2+x  bbar'(x)=3*x**2-4*x+1
!
      real(dp) sum
      integer v
!
      do v=1,ivec
         i=ii(v)
         j=jj(v)
!
!   ...in x direction
!
         xp=xparam(v)
         xpi=1.d0-xp
         xp2=xp*xp
         xpi2=xpi*xpi
         ax=xp2*(3.d0-2.d0*xp)
         axbar=1.d0-ax
         bx=-xp2*xpi
         bxbar=xpi2*xp
!
!   ...in y direction
!
         yp=yparam(v)
         ypi=1.d0-yp
         yp2=yp*yp
         ypi2=ypi*ypi
         ay=yp2*(3.d0-2.d0*yp)
         aybar=1.d0-ay
         by=-yp2*ypi
         bybar=ypi2*yp
!
!   ...derivatives...
!
         axp=6.d0*xp*xpi
         axbarp=-axp
         bxp=xp*(3.d0*xp-2.d0)
         bxbarp=xpi*(3.d0*xpi-2.d0)
!
         ayp=6.d0*yp*ypi
         aybarp=-ayp
         byp=yp*(3.d0*yp-2.d0)
         bybarp=ypi*(3.d0*ypi-2.d0)
!
         iadr=0
!
!  get desired values:
!
         if(ict(1).eq.1) then
!
!  function value:
!
            iadr=iadr+1
            sum=axbar*(aybar*fin(0,i,j)  +ay*fin(0,i,j+1))+
     >             ax*(aybar*fin(0,i+1,j)+ay*fin(0,i+1,j+1))
!
            sum=sum+hx(v)*(
     >         bxbar*(aybar*fin(1,i,j)  +ay*fin(1,i,j+1))+
     >            bx*(aybar*fin(1,i+1,j)+ay*fin(1,i+1,j+1)))
!
            sum=sum+hy(v)*(
     >         axbar*(bybar*fin(2,i,j)  +by*fin(2,i,j+1))+
     >            ax*(bybar*fin(2,i+1,j)+by*fin(2,i+1,j+1)))
!
            sum=sum+hx(v)*hy(v)*(
     >         bxbar*(bybar*fin(3,i,j)  +by*fin(3,i,j+1))+
     >            bx*(bybar*fin(3,i+1,j)+by*fin(3,i+1,j+1)))
!
            fval(v,iadr)=sum
         endif
!
         if(ict(2).eq.1) then
!
!  df/dx:
!
            iadr=iadr+1
!
            sum=hxi(v)*(
     >         axbarp*(aybar*fin(0,i,j)  +ay*fin(0,i,j+1))+
     >            axp*(aybar*fin(0,i+1,j)+ay*fin(0,i+1,j+1)))
!
            sum=sum+
     >         bxbarp*(aybar*fin(1,i,j)  +ay*fin(1,i,j+1))+
     >            bxp*(aybar*fin(1,i+1,j)+ay*fin(1,i+1,j+1))
!
            sum=sum+hxi(v)*hy(v)*(
     >         axbarp*(bybar*fin(2,i,j)  +by*fin(2,i,j+1))+
     >            axp*(bybar*fin(2,i+1,j)+by*fin(2,i+1,j+1)))
!
            sum=sum+hy(v)*(
     >         bxbarp*(bybar*fin(3,i,j)  +by*fin(3,i,j+1))+
     >            bxp*(bybar*fin(3,i+1,j)+by*fin(3,i+1,j+1)))
!
            fval(v,iadr)=sum
         endif
!
         if(ict(3).eq.1) then
!
!  df/dy:
!
            iadr=iadr+1
!
            sum=hyi(v)*(
     >         axbar*(aybarp*fin(0,i,j)  +ayp*fin(0,i,j+1))+
     >            ax*(aybarp*fin(0,i+1,j)+ayp*fin(0,i+1,j+1)))
!
            sum=sum+hx(v)*hyi(v)*(
     >         bxbar*(aybarp*fin(1,i,j)  +ayp*fin(1,i,j+1))+
     >            bx*(aybarp*fin(1,i+1,j)+ayp*fin(1,i+1,j+1)))
!
            sum=sum+
     >         axbar*(bybarp*fin(2,i,j)  +byp*fin(2,i,j+1))+
     >            ax*(bybarp*fin(2,i+1,j)+byp*fin(2,i+1,j+1))
!
            sum=sum+hx(v)*(
     >         bxbar*(bybarp*fin(3,i,j)  +byp*fin(3,i,j+1))+
     >            bx*(bybarp*fin(3,i+1,j)+byp*fin(3,i+1,j+1)))
!
            fval(v,iadr)=sum
         endif
!
         if(ict(4).eq.1) then
!
!  d2f/dxdy:
!
            iadr=iadr+1
!
            sum=hxi(v)*hyi(v)*(
     >         axbarp*(aybarp*fin(0,i,j)  +ayp*fin(0,i,j+1))+
     >            axp*(aybarp*fin(0,i+1,j)+ayp*fin(0,i+1,j+1)))
!
            sum=sum+hyi(v)*(
     >         bxbarp*(aybarp*fin(1,i,j)  +ayp*fin(1,i,j+1))+
     >            bxp*(aybarp*fin(1,i+1,j)+ayp*fin(1,i+1,j+1)))
!
            sum=sum+hxi(v)*(
     >         axbarp*(bybarp*fin(2,i,j)  +byp*fin(2,i,j+1))+
     >            axp*(bybarp*fin(2,i+1,j)+byp*fin(2,i+1,j+1)))
!
            sum=sum+
     >         bxbarp*(bybarp*fin(3,i,j)  +byp*fin(3,i,j+1))+
     >            bxp*(bybarp*fin(3,i+1,j)+byp*fin(3,i+1,j+1))
!
            fval(v,iadr)=sum
         endif
!
      enddo                             ! vector loop
!
      return
      end subroutine herm2fcn_db


!---------------------------------------------------------------------
!  evaluate C1 cubic Hermite function interpolation -- 2d fcn
!   --vectorized-- dmc 10 Feb 1999
!
      subroutine herm2fcn_mesa_db(ict,ivec,ivecd,
     >   fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >   f1,inf2,ny)
!
      integer ny, inf2
      integer ict(4)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
!
      integer ii(ivec),jj(ivec)         ! target cells (i,j)
      real(dp) :: xparam(ivec),yparam(ivec)
                          ! normalized displacements from (i,j) corners
!
      real(dp) :: hx(ivec),hy(ivec)            ! grid spacing, and
      real(dp) :: hxi(ivec),hyi(ivec)          ! inverse grid spacing 1/(x(i+1)-x(i))
                                        ! & 1/(y(j+1)-y(j))
!
      real(dp), pointer :: f1(:) ! =(0:3,inf2,ny)  interpolant data (cf "evbicub")
!
      real(dp) :: fval(ivecd,4)                ! output returned
!
!  for detailed description of fin, ict and fval see subroutine
!  herm2ev comments.  Note ict is not vectorized; the same output
!  is expected to be returned for all input vector data points.
!
!  note that the index inputs ii,jj and parameter inputs
!     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
!     output array fval has a vector ** 1st dimension ** whose
!     size must be given as a separate argument
!
!  to use this routine in scalar mode, pass in ivec=ivecd=1
!
!---------------
!  Hermite cubic basis functions
!  -->for function value matching
!     a(0)=0 a(1)=1        a'(0)=0 a'(1)=0
!   abar(0)=1 abar(1)=0  abar'(0)=0 abar'(1)=0
!
!   a(x)=-2*x**3 + 3*x**2    = x*x*(-2.d0*x+3.0)
!   abar(x)=1-a(x)
!   a'(x)=-abar'(x)          = 6.d0*x*(1.0-x)
!
!  -->for derivative matching
!     b(0)=0 b(1)=0          b'(0)=0 b'(1)=1
!   bbar(0)=0 bbar(1)=0  bbar'(0)=1 bbar'(1)=0
!
!     b(x)=x**3-x**2         b'(x)=3*x**2-2*x
!     bbar(x)=x**3-2*x**2+x  bbar'(x)=3*x**2-4*x+1
!
      real(dp) :: sum
      integer v,z36th,iadr,i,j
      real(dp) :: xp,yp,xpi,ypi,xp2,yp2,xpi2,ypi2
      real(dp) :: cx,cy,cyd,cxi,cyi,cydi,hx2,hy2,cxd,cxdi
      real(dp) :: ax,bx,axbar,bxbar,ay,by,aybar,bybar
      real(dp) :: axp,axbarp,bxp,bxbarp,ayp,aybarp,bybarp,byp
      real(dp), pointer :: fin(:,:,:)
!
      fin(0:3,1:inf2,1:ny) => f1(1:4*inf2*ny)
!
      do v=1,ivec
         i=ii(v)
         j=jj(v)
!
!   ...in x direction
!
         xp=xparam(v)
         xpi=1.d0-xp
         xp2=xp*xp
         xpi2=xpi*xpi
         ax=xp2*(3.d0-2.d0*xp)
         axbar=1.d0-ax
         bx=-xp2*xpi
         bxbar=xpi2*xp
!
!   ...in y direction
!
         yp=yparam(v)
         ypi=1.d0-yp
         yp2=yp*yp
         ypi2=ypi*ypi
         ay=yp2*(3.d0-2.d0*yp)
         aybar=1.d0-ay
         by=-yp2*ypi
         bybar=ypi2*yp
!
!   ...derivatives...
!
         axp=6.d0*xp*xpi
         axbarp=-axp
         bxp=xp*(3.d0*xp-2.d0)
         bxbarp=xpi*(3.d0*xpi-2.d0)
!
         ayp=6.d0*yp*ypi
         aybarp=-ayp
         byp=yp*(3.d0*yp-2.d0)
         bybarp=ypi*(3.d0*ypi-2.d0)
!
         iadr=0
!
!  get desired values:
!
         if(ict(1).eq.1) then
!
!  function value:
!
            iadr=iadr+1
            sum=axbar*(aybar*fin(0,i,j)  +ay*fin(0,i,j+1))+
     >             ax*(aybar*fin(0,i+1,j)+ay*fin(0,i+1,j+1))
!
            sum=sum+hx(v)*(
     >         bxbar*(aybar*fin(1,i,j)  +ay*fin(1,i,j+1))+
     >            bx*(aybar*fin(1,i+1,j)+ay*fin(1,i+1,j+1)))
!
            sum=sum+hy(v)*(
     >         axbar*(bybar*fin(2,i,j)  +by*fin(2,i,j+1))+
     >            ax*(bybar*fin(2,i+1,j)+by*fin(2,i+1,j+1)))
!
            sum=sum+hx(v)*hy(v)*(
     >         bxbar*(bybar*fin(3,i,j)  +by*fin(3,i,j+1))+
     >            bx*(bybar*fin(3,i+1,j)+by*fin(3,i+1,j+1)))
!
            fval(v,iadr)=sum
         endif
!
         if(ict(2).eq.1) then
!
!  df/dx:
!
            iadr=iadr+1
!
            sum=hxi(v)*(
     >         axbarp*(aybar*fin(0,i,j)  +ay*fin(0,i,j+1))+
     >            axp*(aybar*fin(0,i+1,j)+ay*fin(0,i+1,j+1)))
!
            sum=sum+
     >         bxbarp*(aybar*fin(1,i,j)  +ay*fin(1,i,j+1))+
     >            bxp*(aybar*fin(1,i+1,j)+ay*fin(1,i+1,j+1))
!
            sum=sum+hxi(v)*hy(v)*(
     >         axbarp*(bybar*fin(2,i,j)  +by*fin(2,i,j+1))+
     >            axp*(bybar*fin(2,i+1,j)+by*fin(2,i+1,j+1)))
!
            sum=sum+hy(v)*(
     >         bxbarp*(bybar*fin(3,i,j)  +by*fin(3,i,j+1))+
     >            bxp*(bybar*fin(3,i+1,j)+by*fin(3,i+1,j+1)))
!
            fval(v,iadr)=sum
         endif
!
         if(ict(3).eq.1) then
!
!  df/dy:
!
            iadr=iadr+1
!
            sum=hyi(v)*(
     >         axbar*(aybarp*fin(0,i,j)  +ayp*fin(0,i,j+1))+
     >            ax*(aybarp*fin(0,i+1,j)+ayp*fin(0,i+1,j+1)))
!
            sum=sum+hx(v)*hyi(v)*(
     >         bxbar*(aybarp*fin(1,i,j)  +ayp*fin(1,i,j+1))+
     >            bx*(aybarp*fin(1,i+1,j)+ayp*fin(1,i+1,j+1)))
!
            sum=sum+
     >         axbar*(bybarp*fin(2,i,j)  +byp*fin(2,i,j+1))+
     >            ax*(bybarp*fin(2,i+1,j)+byp*fin(2,i+1,j+1))
!
            sum=sum+hx(v)*(
     >         bxbar*(bybarp*fin(3,i,j)  +byp*fin(3,i,j+1))+
     >            bx*(bybarp*fin(3,i+1,j)+byp*fin(3,i+1,j+1)))
!
            fval(v,iadr)=sum
         endif
!
         if(ict(4).eq.1) then
!
!  d2f/dxdy:
!
            iadr=iadr+1
!
            sum=hxi(v)*hyi(v)*(
     >         axbarp*(aybarp*fin(0,i,j)  +ayp*fin(0,i,j+1))+
     >            axp*(aybarp*fin(0,i+1,j)+ayp*fin(0,i+1,j+1)))
!
            sum=sum+hyi(v)*(
     >         bxbarp*(aybarp*fin(1,i,j)  +ayp*fin(1,i,j+1))+
     >            bxp*(aybarp*fin(1,i+1,j)+ayp*fin(1,i+1,j+1)))
!
            sum=sum+hxi(v)*(
     >         axbarp*(bybarp*fin(2,i,j)  +byp*fin(2,i,j+1))+
     >            axp*(bybarp*fin(2,i+1,j)+byp*fin(2,i+1,j+1)))
!
            sum=sum+
     >         bxbarp*(bybarp*fin(3,i,j)  +byp*fin(3,i,j+1))+
     >            bxp*(bybarp*fin(3,i+1,j)+byp*fin(3,i+1,j+1))
!
            fval(v,iadr)=sum
         endif
!
      enddo                             ! vector loop
!
      return
      end subroutine herm2fcn_mesa_db


      subroutine ibc_ck_db(ibc,slbl,xlbl,imin,imax,ier)
!
!  check that spline routine ibc flag is in range
!
      implicit none
      integer ibc                       ! input -- flag value
      character*(*) slbl                ! input -- subroutine name
      character*(*) xlbl                ! input -- axis label
!
      integer imin                      ! input -- min allowed value
      integer imax                      ! input -- max allowed value
!
      integer ier                       ! output -- set =1 if error detected
!
!----------------------
!
      if((ibc.lt.imin).or.(ibc.gt.imax)) then
         ier=1
!         write(6,1001) slbl,xlbl,ibc,imin,imax
! 1001    format(' ?',a,' -- ibc',a,' = ',i9,' out of range ',
!     >      i2,' to ',i2)
      endif
!
      return
      end subroutine ibc_ck_db


      subroutine do_mkbicub_db(x,nx,y,ny,f1,nf2,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ibcymin,bcymin,ibcymax,bcymax,
     >   ilinx,iliny,ier)
!
!  setup bicubic spline, dynamic allocation of workspace
!  fortran-90 fixed form source
!
!  --NOTE-- dmc 22 Feb 2004 -- rewrite for direct calculation of
!  coefficients, to avoid large transient use of memory.
!
!
      implicit NONE
!
!  input:
      integer nx                        ! length of x vector
      integer ny                        ! length of y vector
      real(dp) x(:) ! (nx)                        ! x vector, strict ascending
      real(dp) y(:) ! (ny)                        ! y vector, strict ascending
!
      integer nf2                       ! 2nd dimension of f, nf2.ge.nx
!  input/output:
      real(dp), pointer :: f1(:) ! =(4,nf2,ny)                  ! data & spline coefficients
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
      integer ibcxmin                   ! bc flag for x=xmin
      real(dp) bcxmin(:) ! (ny)                   ! bc data vs. y at x=xmin
      integer ibcxmax                   ! bc flag for x=xmax
      real(dp) bcxmax(:) ! (ny)                   ! bc data vs. y at x=xmax
!
      integer ibcymin                   ! bc flag for y=ymin
      real(dp) bcymin(:) ! (nx)                   ! bc data vs. x at y=ymin
      integer ibcymax                   ! bc flag for y=ymax
      real(dp) bcymax(:) ! (nx)                   ! bc data vs. x at y=ymax
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
!           (for more detailed definition of BCs 5-7, see the
!           comments of subroutine mkspline)
!   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
!
!   ibcxmax -- indicator for boundary condition at x(nx):
!    bcxmax(...) -- boundary condition data
!     (interpretation as with ibcxmin, bcxmin)
!   NOTE:  if ibcxmin=-1, ibcxmax is ignored! ...and the BC is periodic.
!
!  and analogous interpretation for ibcymin,bcymin,ibcymax,bcymax
!  (df/dy or d2f/dy2 boundary conditions at y=ymin and y=ymax).
!
!  output linear grid flags and completion code (ier=0 is normal):
!
      integer ilinx                     ! =1: x grid is "nearly" equally spaced
      integer iliny                     ! =1: y grid is "nearly" equally spaced
!  ilinx and iliny are set to zero if corresponding grids are not equally
!  spaced
!
      integer ier                       ! =0 on exit if there is no error.
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
      integer ierx,iery
!
      real(dp), dimension(:,:), allocatable :: fwk
      real(dp) :: zbcmin,zbcmax
      integer ix,iy,ibcmin,ibcmax
!
      real(dp), dimension(:,:,:), allocatable :: fcorr
      integer iflg2
      real(dp) zdiff(2),hy
      
      real(dp), pointer :: f(:,:,:) ! =(4,nf2,ny)                  ! data & spline coefficients
      f(1:4,1:nf2,1:ny) => f1(1:4*nf2*ny)

!
!-----------------------
!
!  see if 2nd pass is needed due to inhomogeneous d/dy bdy cond.
!
      iflg2=0
      if(ibcymin.ne.-1) then
         if((ibcymin.eq.1).or.(ibcymin.eq.2)) then
            do ix=1,nx
               if (bcymin(ix).ne.0.d0) iflg2=1
            enddo
         endif
         if((ibcymax.eq.1).or.(ibcymax.eq.2)) then
            do ix=1,nx
               if (bcymax(ix).ne.0.d0) iflg2=1
            enddo
         endif
      endif
!
!  check boundary condition specifications
!
      ier=0
!
      call ibc_ck_db(ibcxmin,'bcspline','xmin',-1,7,ier)
      if(ibcxmin.ge.0) call ibc_ck_db(ibcxmax,'bcspline','xmax',0,7,ier)
      call ibc_ck_db(ibcymin,'bcspline','ymin',-1,7,ier)
      if(ibcymin.ge.0) call ibc_ck_db(ibcymax,'bcspline','ymax',0,7,ier)
!
!  check ilinx & x vector
!
      call splinck_db(x,nx,ilinx,1.0d-3,ierx)
      if(ierx.ne.0) ier=2
!
!      if(ier.eq.2) then
!         write(6,'('' ?bcspline:  x axis not strict ascending'')')
!      endif
!
!  check iliny & y vector
!
      call splinck_db(y,ny,iliny,1.0d-3,iery)
      if(iery.ne.0) ier=3
!
!      if(ier.eq.3) then
!         write(6,'('' ?bcspline:  y axis not strict ascending'')')
!      endif
!
      if(ier.ne.0) return
!
!------------------------------------
      allocate(fwk(2,max(nx,ny)))
!
!  evaluate fxx (spline in x direction)
!
      zbcmin=0
      zbcmax=0
      do iy=1,ny
         fwk(1,1:nx) = f(1,1:nx,iy)
         if((ibcxmin.eq.1).or.(ibcxmin.eq.2)) zbcmin=bcxmin(iy)
         if((ibcxmax.eq.1).or.(ibcxmax.eq.2)) zbcmax=bcxmax(iy)
         call mkspline_db(x,nx,fwk,
     >      ibcxmin,zbcmin,ibcxmax,zbcmax,ilinx,ier)
         if(ier.ne.0) then
            deallocate(fwk)
            return
         end if
         f(2,1:nx,iy)=fwk(2,1:nx)
      enddo
!
!  evaluate fyy (spline in y direction)
!  use homogeneous boundary condition; correction done later if necessary
!
      zbcmin=0
      zbcmax=0
      ibcmin=ibcymin
      ibcmax=ibcymax
      do ix=1,nx
         fwk(1,1:ny) = f(1,ix,1:ny)
         if(iflg2.eq.1) then
            if((ibcymin.eq.1).or.(ibcymin.eq.2)) ibcmin=0
            if((ibcymax.eq.1).or.(ibcymax.eq.2)) ibcmax=0
         endif
         call mkspline_db(y,ny,fwk,
     >      ibcmin,zbcmin,ibcmax,zbcmax,iliny,ier)
         if(ier.ne.0) then
            deallocate(fwk)
            return
         end if
         f(3,ix,1:ny)=fwk(2,1:ny)
      enddo
!
!  evaluate fxxyy (spline fxx in y direction; BC simplified; avg
!  d2(d2f/dx2)/dy2 and d2(df2/dy2)/dx2
!
      zbcmin=0
      zbcmax=0
      ibcmin=ibcymin
      ibcmax=ibcymax
      do ix=1,nx
         fwk(1,1:ny) = f(2,ix,1:ny)
         if(iflg2.eq.1) then
            if((ibcymin.eq.1).or.(ibcymin.eq.2)) ibcmin=0
            if((ibcymax.eq.1).or.(ibcymax.eq.2)) ibcmax=0
         endif
         call mkspline_db(y,ny,fwk,
     >      ibcmin,zbcmin,ibcmax,zbcmax,iliny,ier)
         if(ier.ne.0) then
            deallocate(fwk)
            return
         end if
         f(4,ix,1:ny)= fwk(2,1:ny)
      enddo
!
      if(iflg2.eq.1) then
         allocate(fcorr(2,nx,ny))
!
!  correct for inhomogeneous y boundary condition
!
         do ix=1,nx
            !  the desired inhomogenous BC is the difference btw the 
            !  requested derivative (1st or 2nd) and the current value

            zdiff(1)=0.d0
            if(ibcymin.eq.1) then
               hy=y(2)-y(1)
               zdiff(1)=(f(1,ix,2)-f(1,ix,1))/hy +
     >            hy*(-2*f(3,ix,1)-f(3,ix,2))/6
               zdiff(1)=bcymin(ix)-zdiff(1)
            else if(ibcymin.eq.2) then
               zdiff(1)=bcymin(ix)-f(3,ix,1)
            endif

            zdiff(2)=0.d0
            if(ibcymax.eq.1) then
               hy=y(ny)-y(ny-1)
               zdiff(2)=(f(1,ix,ny)-f(1,ix,ny-1))/hy + 
     >            hy*(2*f(3,ix,ny)+f(3,ix,ny-1))/6
               zdiff(2)=bcymax(ix)-zdiff(2)
            else if(ibcymax.eq.2) then
               zdiff(2)=bcymax(ix)-f(3,ix,ny)
            endif
!
            fwk(1,1:ny)=0.d0  ! values are zero; only BC is not
            call mkspline_db(y,ny,fwk,ibcymin,zdiff(1),ibcymax,zdiff(2),
     >         iliny,ier)
            if(ier.ne.0) then
               deallocate(fwk,fcorr)
               return
            end if
            fcorr(1,ix,1:ny)=fwk(2,1:ny)  ! fyy-correction
         enddo
!
         zbcmin=0
         zbcmax=0
         do iy=1,ny
            fwk(1,1:nx)=fcorr(1,1:nx,iy)
            call mkspline_db(x,nx,fwk,ibcxmin,zbcmin,ibcxmax,zbcmax,
     >         ilinx,ier)
            if(ier.ne.0) then
               deallocate(fwk,fcorr)
               return
            end if
            fcorr(2,1:nx,iy)=fwk(2,1:nx)  ! fxxyy-correction
         enddo
!
         f(3:4,1:nx,1:ny)=f(3:4,1:nx,1:ny)+fcorr(1:2,1:nx,1:ny)
!
         deallocate(fcorr)        
      endif
!
!  correction spline -- f=fxx=zero; fyy & fxxyy are affected
!
      deallocate(fwk)
!------------------------------------
!
!  thats all
!
      return
      end subroutine do_mkbicub_db



      subroutine mkspline_db(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ilinx,ier)
      implicit none
!
!  make a 2-coefficient 1d spline
!
!  only 2 coefficients, the data and its 2nd derivative, are needed to
!  fully specify a spline.  See e.g. Numerical Recipies in Fortran-77
!  (2nd edition) chapter 3, section on cubic splines.
!
!  input:
      integer nx                        ! no. of data points
      real(dp) x(nx)                        ! x axis data, strict ascending order
!
!  input/output:
      real(dp) fspl(2,nx)                   ! f(1,*):  data in; f(2,*):  coeffs out
!
!     f(1,j) = f(x(j))  on input (unchanged on output)
!     f(2,j) = f''(x(j)) (of interpolating spline) (on output).
!
!  ...boundary conditions...
!
!  input:
!
      integer ibcxmin                   ! b.c. flag @ x=xmin=x(1)
      real(dp) bcxmin                       ! b.c. data @xmin
!
      integer ibcxmax                   ! b.c. flag @ x=xmax=x(nx)
      real(dp) bcxmax                       ! b.c. data @xmax
!
!  ibcxmin=-1 -- periodic boundary condition
!                (bcxmin,ibcxmax,bcxmax are ignored)
!
!                the output spline s satisfies
!                s'(x(1))=s'(x(nx)) ..and.. s''(x(1))=s''(x(nx))
!
!  if non-periodic boundary conditions are used, then the xmin and xmax
!  boundary conditions can be specified independently:
!
!  ibcxmin (ibcxmax) = 0 -- this specifies a "not a knot" boundary
!                condition, see "cubsplb.for".  This is a common way
!                for inferring a "good" spline boundary condition
!                automatically from data in the vicinity of the
!                boundary.  ... bcxmin (bcxmax) are ignored.
!
!  ibcxmin (ibcxmax) = 1 -- boundary condition is to have s'(x(1))
!                ( s'(x(nx)) ) match the passed value bcxmin (bcxmax).
!
!  ibcxmin (ibcxmax) = 2 -- boundary condition is to have s''(x(1))
!                ( s''(x(nx)) ) match the passed value bcxmin (bcxmax).
!
!  ibcxmin (ibcxmax) = 3 -- boundary condition is to have s'(x(1))=0.d0
!                ( s'(x(nx))=0.d0 )
!
!  ibcxmin (ibcxmax) = 4 -- boundary condition is to have s''(x(1))=0.d0
!                ( s''(x(nx))=0.d0 )
!
!  ibcxmin (ibcxmax) = 5 -- boundary condition is to have s'(x(1))
!                ( s'(x(nx)) ) match the 1st divided difference
!                e.g. at x(1):  d(1)/h(1), where
!                           d(j)=f(1,j+1)-f(1,j)
!                           h(j)=x(j+1)-x(j)
!
!  ibcxmin (ibcxmax) = 6 -- BC is to have s''(x(1)) ( s''(x(nx)) )
!                match the 2nd divided difference
!                e.g. at x(1):
!                     e(1) = [d(2)/h(2) - d(1)/h(1)]/(0.5*(h(1)+h(2)))
!
!  ibcxmin (ibcxmax) = 7 -- BC is to have s'''(x(1)) ( s'''(x(nx)) )
!                match the 3rd divided difference
!                e.g. at x(1): [e(2)-e(1)]/(0.33333*(h(1)+h(2)+h(3)))
!
!  output:
!
      integer ilinx                     ! =1: hint, x axis is ~evenly spaced
!
!  let dx[avg] = (x(nx)-x(1))/(nx-1)
!  let dx[j] = x(j+1)-x(j), for all j satisfying 1.le.j.lt.nx
!
!  if for all such j, abs(dx[j]-dx[avg]).le.(1.0d-3*dx[avg]) then
!  ilinx=1 is returned, indicating the data is (at least nearly)
!  evenly spaced.  Even spacing is useful, for speed of zone lookup,
!  when evaluating a spline.
!
!  if the even spacing condition is not satisfied, ilinx=2 is returned.
!
      integer ier                       ! exit code, 0=OK
!
!  an error code is returned if the x axis is not strict ascending,
!  or if nx.lt.4, or if an invalid boundary condition specification was
!  input.
!
!------------------------------------
!
!  this routine calls traditional 4 coefficient spline software, and
!  translates the result to 2 coefficient form.
!
!  this could be done more efficiently but we decided out of conservatism
!  to use the traditional software.
!
!------------------------------------
      integer i, inwk
!  workspaces -- f90 dynamically allocated
!
      real(dp), dimension(:,:), allocatable :: fspl4 ! traditional 4-spline
      real(dp), dimension(:), allocatable :: wk ! cspline workspace
!
!------------------------------------
!
      allocate(fspl4(4,nx),wk(nx))
!
!  make the traditional call
!
      do i=1,nx
         fspl4(1,i)=fspl(1,i)
         fspl(2,i)=0.d0                  ! for now
      enddo
!
      inwk=nx
!
!  boundary conditions imposed by cspline...
!
      call cspline_db(x,nx,fspl4,ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   wk,inwk,ilinx,ier)
!
      if(ier.eq.0) then
!
!  copy the output -- careful of end point.
!
         do i=1,nx-1
            fspl(2,i)=2.d0*fspl4(3,i)
         enddo
         fspl(2,nx)=2.d0*fspl4(3,nx-1) +
     >        (x(nx)-x(nx-1))*6.d0*fspl4(4,nx-1)
      endif
!
      deallocate(fspl4,wk)
!
      return
      end subroutine mkspline_db


      subroutine splinck_db(x,inx,ilinx,ztol,ier)
      implicit none
!
!  check if a grid is strictly ascending and if it is evenly spaced
!  to w/in ztol
!
      integer inx
      real(dp) x(inx)                       ! input -- grid to check
!
      integer ilinx                     ! output -- =1 if evenly spaced =2 O.W.
!
      real(dp) ztol                         ! input -- spacing check tolerance
!
      integer ier                       ! output -- =0 if OK
!
!  ier=1 is returned if x(1...inx) is NOT STRICTLY ASCENDING...
!
!-------------------------------
!
      real(dp) zeps, zdiffx, zdiff, dxavg
      integer ix
      
      ier=0
      ilinx=1
      if(inx.le.1) return
!
      dxavg=(x(inx)-x(1))/(inx-1)
      zeps=abs(ztol*dxavg)
!
      do ix=2,inx
         zdiffx=(x(ix)-x(ix-1))
         if(zdiffx.le.0.0) ier=2
         zdiff=zdiffx-dxavg
         if(abs(zdiff).gt.zeps) then
            ilinx=2
         endif
      enddo
 10   continue
!
      return
      end subroutine splinck_db

      SUBROUTINE V_SPLINE_db(k_bc1,k_bcn,n,x,f,wk)
!***********************************************************************
!V_SPLINE evaluates the coefficients for a 1d cubic interpolating spline
!References:
!  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!    Computations, Prentice-Hall, 1977, p.76
!  Engeln-Muellges, Uhlig, Numerical Algorithms with Fortran, Springer,
!    1996, p.251
!  W.A.Houlberg, D.McCune 3/2000
!Input:
!  k_bc1-option for BC at x(1)
!       =-1 periodic, ignore k_bcn
!       =0 not-a-knot
!       =1 s'(x1) = input value of f(2,1)
!       =2 s''(x1) = input value of f(3,1)
!       =3 s'(x1) = 0.0
!       =4 s''(x1) = 0.0
!       =5 match first derivative to first 2 points
!       =6 match second derivative to first 3 points
!       =7 match third derivative to first 4 points
!       =else use not-a-knot
!  k_bcn-option for boundary condition at x(n)
!       =0 not-a-knot
!       =1 s'(x1) = input value of f(2,1)
!       =2 s''(x1) = input value of f(3,1)
!       =3 s'(x1) = 0.0
!       =4 s''(x1) = 0.0
!       =5 match first derivative to first 2 points
!       =6 match second derivative to first 3 points
!       =7 match third derivative to first 4 points
!       =else use knot-a-knot
!  n-number of data points or knots-(n.ge.2)
!  x(n)-abscissas of the knots in strictly increasing order
!  f(1,i)-ordinates of the knots
!  f(2,1)-input value of s'(x1) for k_bc1=1
!  f(2,n)-input value of s'(xn) for k_bcn=1
!  f(3,1)-input value of s''(x1) for k_bc1=2
!  f(3,n)-input value of s''(xn) for k_bcn=2
!  wk(n)-scratch work area for periodic BC
!Output:
!  f(2,i)=s'(x(i))
!  f(3,i)=s''(x(i))
!  f(4,i)=s'''(x(i))
!Comments:
!  s(x)=f(1,i)+f(2,i)*(x-x(i))+f(3,i)*(x-x(i))**2/2!
!       +f(4,i)*(x-x(i))**3/3! for x(i).le.x.le.x(i+1)
!  W_SPLINE can be used to evaluate the spline and its derivatives
!  The cubic spline is twice differentiable (C2)
!
!  bugfixes -- dmc 24 Feb 2004:
!    (a) fixed logic for not-a-knot:
!          !    Set f(3,1) for not-a-knot
!                    IF(k_bc1.le.0.or.k_bc1.gt.7) THEN ...
!        instead of
!          !    Set f(3,1) for not-a-knot
!                    IF(k_bc1.le.0.or.k_bc1.gt.5) THEN ...
!        and similarly for logic after cmt
!          !    Set f(3,n) for not-a-knot
!        as required since k_bc*=6 and k_bc*=7 are NOT not-a-knot BCs.
!
!    (b) the BCs to fix 2nd derivative at end points did not work if that
!        2nd derivative were non-zero.  The reason is that in those cases
!        the off-diagonal matrix elements nearest the corners are not
!        symmetric; i.e. elem(1,2).ne.elem(2,1) and 
!        elem(n-1,n).ne.elem(n,n-1) where I use "elem" to refer to
!        the tridiagonal matrix elements.  The correct values for the
!        elements is:   elem(1,2)=0, elem(2,1)=x(2)-x(1)
!                       elem(n,n-1)=0, elem(n-1,n)=x(n)-x(n-1)
!        the old code in effect had these as all zeroes.  Since this
!        meant the wrong set of linear equations was solved, the
!        resulting spline had a discontinuity in its 1st derivative
!        at x(2) and x(n-1).  Fixed by introducing elem21 and elemnn1
!        to represent the non-symmetric lower-diagonal values.  Since
!        elem21 & elemnn1 are both on the lower diagonals, logic to 
!        use them occurs in the non-periodic forward elimination loop
!        only.  DMC 24 Feb 2004.
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        k_bc1,                   k_bcn,
     &               n
      real(dp)           x(*),                    wk(*),
     &               f(4,*)
!Declaration in local variables
      INTEGER        i,                       ib,
     &               imax,                    imin
      real(dp)           a1,                      an,
     &               b1,                      bn,
     &               q,                       t,
     &               hn
      real(dp)           elem21,                  elemnn1    ! (dmc)

!Set default range
      imin=1
      imax=n
!Set first and second BC values
      a1=0.d0
      b1=0.d0
      an=0.d0
      bn=0.d0
      IF(k_bc1.eq.1) THEN
        a1=f(2,1)
      ELSEIF(k_bc1.eq.2) THEN
        b1=f(3,1)
      ELSEIF(k_bc1.eq.5) THEN
        a1=(f(1,2)-f(1,1))/(x(2)-x(1))
      ELSEIF(k_bc1.eq.6) THEN
        b1=2.d0*((f(1,3)-f(1,2))/(x(3)-x(2))
     &         -(f(1,2)-f(1,1))/(x(2)-x(1)))/(x(3)-x(1))
      ENDIF
      IF(k_bcn.eq.1) THEN
        an=f(2,n)
      ELSEIF(k_bcn.eq.2) THEN
        bn=f(3,n)
      ELSEIF(k_bcn.eq.5) THEN
        an=(f(1,n)-f(1,n-1))/(x(n)-x(n-1))
      ELSEIF(k_bcn.eq.6) THEN
        bn=2.d0*((f(1,n)-f(1,n-1))/(x(n)-x(n-1))
     &         -(f(1,n-1)-f(1,n-2))/(x(n-1)-x(n-2)))/(x(n)-x(n-2))
      ENDIF
!Clear f(2:4,n)
      f(2,n)=0.d0
      f(3,n)=0.d0
      f(4,n)=0.d0
      IF(n.eq.2) THEN
!Coefficients for n=2
        f(2,1)=(f(1,2)-f(1,1))/(x(2)-x(1))
        f(3,1)=0.d0
        f(4,1)=0.d0
        f(2,2)=f(2,1)
        f(3,2)=0.d0
        f(4,2)=0.d0
      ELSEIF(n.gt.2) THEN
!Set up tridiagonal system for A*y=B where y(i) are the second
!  derivatives at the knots
!  f(2,i) are the diagonal elements of A
!  f(4,i) are the off-diagonal elements of A
!  f(3,i) are the B elements/3, and will become c/3 upon solution
        f(4,1)=x(2)-x(1)
        f(3,2)=(f(1,2)-f(1,1))/f(4,1)
        DO i=2,n-1
          f(4,i)=x(i+1)-x(i)
          f(2,i)=2.d0*(f(4,i-1)+f(4,i))
          f(3,i+1)=(f(1,i+1)-f(1,i))/f(4,i)
          f(3,i)=f(3,i+1)-f(3,i)
        ENDDO
!
!  (dmc): save now:
!
        elem21=f(4,1)
        elemnn1=f(4,n-1)
!
!  BC's
!    Left
        IF(k_bc1.eq.-1) THEN
          f(2,1)=2.d0*(f(4,1)+f(4,n-1))
          f(3,1)=(f(1,2)-f(1,1))/f(4,1)-(f(1,n)-f(1,n-1))/f(4,n-1)
          wk(1)=f(4,n-1)
          DO i=2,n-3
            wk(i)=0.d0
          ENDDO
          wk(n-2)=f(4,n-2)
          wk(n-1)=f(4,n-1)
        ELSEIF(k_bc1.eq.1.or.k_bc1.eq.3.or.k_bc1.eq.5) THEN
          f(2,1)=2.d0*f(4,1)
          f(3,1)=(f(1,2)-f(1,1))/f(4,1)-a1
        ELSEIF(k_bc1.eq.2.or.k_bc1.eq.4.or.k_bc1.eq.6) THEN
          f(2,1)=2.d0*f(4,1)
          f(3,1)=f(4,1)*b1/3.d0
          f(4,1)=0.d0  ! upper diagonal only (dmc: cf elem21)
        ELSEIF(k_bc1.eq.7) THEN
          f(2,1)=-f(4,1)
          f(3,1)=f(3,3)/(x(4)-x(2))-f(3,2)/(x(3)-x(1))
          f(3,1)=f(3,1)*f(4,1)*f(4,1)/(x(4)-x(1))
        ELSE                             ! not a knot:
          imin=2
          f(2,2)=f(4,1)+2.d0*f(4,2)
          f(3,2)=f(3,2)*f(4,2)/(f(4,1)+f(4,2))
        ENDIF
!    Right
        IF(k_bcn.eq.1.or.k_bcn.eq.3.or.k_bcn.eq.5) THEN
          f(2,n)=2.d0*f(4,n-1)
          f(3,n)=-(f(1,n)-f(1,n-1))/f(4,n-1)+an
        ELSEIF(k_bcn.eq.2.or.k_bcn.eq.4.or.k_bcn.eq.6) THEN
          f(2,n)=2.d0*f(4,n-1)
          f(3,n)=f(4,n-1)*bn/3.d0
!xxx          f(4,n-1)=0.d0  ! dmc: preserve f(4,n-1) for back subst.
          elemnn1=0.d0  !  lower diaganol only (dmc)
        ELSEIF(k_bcn.eq.7) THEN
          f(2,n)=-f(4,n-1)
          f(3,n)=f(3,n-1)/(x(n)-x(n-2))-f(3,n-2)/(x(n-1)-x(n-3))
          f(3,n)=-f(3,n)*f(4,n-1)*f(4,n-1)/(x(n)-x(n-3))
        ELSEIF(k_bc1.ne.-1) THEN         ! not a knot:
          imax=n-1
          f(2,n-1)=2.d0*f(4,n-2)+f(4,n-1)
          f(3,n-1)=f(3,n-1)*f(4,n-2)/(f(4,n-1)+f(4,n-2))
        ENDIF
!  Limit solution for only three points in domain
        IF(n.eq.3) THEN
          f(3,1)=0.d0
          f(3,n)=0.d0
        ENDIF
        IF(k_bc1.eq.-1) THEN
!Solve system of equations for second derivatives at the knots
!  Periodic BC
!    Forward elimination
          DO i=2,n-2
            t=f(4,i-1)/f(2,i-1)
            f(2,i)=f(2,i)-t*f(4,i-1)
            f(3,i)=f(3,i)-t*f(3,i-1)
            wk(i)=wk(i)-t*wk(i-1)
            q=wk(n-1)/f(2,i-1)
            wk(n-1)=-q*f(4,i-1)
            f(2,n-1)=f(2,n-1)-q*wk(i-1)
            f(3,n-1)=f(3,n-1)-q*f(3,i-1)
          ENDDO
!    Correct the n-1 element
          wk(n-1)=wk(n-1)+f(4,n-2)
!    Complete the forward elimination
!    wk(n-1) and wk(n-2) are the off-diag elements of the lower corner
          t=wk(n-1)/f(2,n-2)
          f(2,n-1)=f(2,n-1)-t*wk(n-2)
          f(3,n-1)=f(3,n-1)-t*f(3,n-2)
!    Back substitution
          f(3,n-1)=f(3,n-1)/f(2,n-1)
          f(3,n-2)=(f(3,n-2)-wk(n-2)*f(3,n-1))/f(2,n-2)
          DO ib=3,n-1
            i=n-ib
            f(3,i)=(f(3,i)-f(4,i)*f(3,i+1)-wk(i)*f(3,n-1))/f(2,i)
          ENDDO
          f(3,n)=f(3,1)
        ELSE
!  Non-periodic BC
!    Forward elimination
!    For Not-A-Knot BC the off-diagonal end elements are not equal
          DO i=imin+1,imax
            IF((i.eq.n-1).and.(imax.eq.n-1)) THEN
              t=(f(4,i-1)-f(4,i))/f(2,i-1)
            ELSE
              if(i.eq.2) then
                 t=elem21/f(2,i-1)
              else if(i.eq.n) then
                 t=elemnn1/f(2,i-1)
              else
                 t=f(4,i-1)/f(2,i-1)
              endif
            ENDIF
            IF((i.eq.imin+1).and.(imin.eq.2)) THEN
              f(2,i)=f(2,i)-t*(f(4,i-1)-f(4,i-2))
            ELSE
              f(2,i)=f(2,i)-t*f(4,i-1)
            ENDIF
            f(3,i)=f(3,i)-t*f(3,i-1)
          ENDDO
!    Back substitution
          f(3,imax)=f(3,imax)/f(2,imax)
          DO ib=1,imax-imin
            i=imax-ib
            IF((i.eq.2).and.(imin.eq.2)) THEN
              f(3,i)=(f(3,i)-(f(4,i)-f(4,i-1))*f(3,i+1))/f(2,i)
            ELSE
              f(3,i)=(f(3,i)-f(4,i)*f(3,i+1))/f(2,i)
            ENDIF
          ENDDO
!    Reset d array to step size
          f(4,1)=x(2)-x(1)
          f(4,n-1)=x(n)-x(n-1)
!    Set f(3,1) for not-a-knot
          IF(k_bc1.le.0.or.k_bc1.gt.7) THEN
            f(3,1)=(f(3,2)*(f(4,1)+f(4,2))-f(3,3)*f(4,1))/f(4,2)
          ENDIF
!    Set f(3,n) for not-a-knot
          IF(k_bcn.le.0.or.k_bcn.gt.7) THEN
            f(3,n)=f(3,n-1)+(f(3,n-1)-f(3,n-2))*f(4,n-1)/f(4,n-2)
          ENDIF
        ENDIF
!f(3,i) is now the sigma(i) of the text and f(4,i) is the step size
!Compute polynomial coefficients
        DO i=1,n-1
          f(2,i)=
     >        (f(1,i+1)-f(1,i))/f(4,i)-f(4,i)*(f(3,i+1)+2.d0*f(3,i))
          f(4,i)=(f(3,i+1)-f(3,i))/f(4,i)
          f(3,i)=6.d0*f(3,i)
          f(4,i)=6.d0*f(4,i)
        ENDDO
        IF(k_bc1.eq.-1) THEN
          f(2,n)=f(2,1)
          f(3,n)=f(3,1)
          f(4,n)=f(4,1)
        ELSE
           hn=x(n)-x(n-1)
           f(2,n)=f(2,n-1)+hn*(f(3,n-1)+0.5d0*hn*f(4,n-1))
           f(3,n)=f(3,n-1)+hn*f(4,n-1)
           f(4,n)=f(4,n-1)
           IF(k_bcn.eq.1.or.k_bcn.eq.3.or.k_bcn.eq.5) THEN
              f(2,n)=an
           ELSE IF(k_bcn.eq.2.or.k_bcn.eq.4.or.k_bcn.eq.6) THEN
              f(3,n)=bn
           ENDIF
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE V_SPLINE_db


      subroutine zonfind_db(x,nx,zxget,i)
      implicit none
!
      integer nx
      real(dp) x(nx),zxget
      integer i
      
      integer nxm, i1, i2, ii, ij
      real(dp) dx
!
!  find index i such that x(i).le.zxget.le.x(i+1)
!
!  x(1...nx) is strict increasing and x(1).le.zxget.le.x(nx)
!  (this is assumed to already have been checked -- no check here!)
!
      nxm=nx-1
      if((i.lt.1).or.(i.gt.nxm)) then
         i1=1
         i2=nx-1
         go to 10
      endif
!
      if(x(i).gt.zxget) then
!  look down
         dx=x(i+1)-x(i)
         if((x(i)-zxget).gt.4*dx) then
            i1=1
            i2=i-1
            go to 10
         else
            i2=i-1
            do ij=i2,1,-1
               if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
                  i=ij
                  return
               endif
            enddo
            i=1
            return
         endif
      else if(x(i+1).lt.zxget) then
!  look up
         dx=x(i+1)-x(i)
         if((zxget-x(i+1)).gt.4*dx) then
            i1=i+1
            i2=nxm
            go to 10
         else
            i2=i+1
            do ij=i2,nxm
               if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
                  i=ij
                  return
               endif
            enddo
            ij=nxm
            return
         endif
      else
!  already there...
         return
      endif
!
!---------------------------
!  binary search
!
 10   continue
!
      if(i1.eq.i2) then
! found by proc. of elimination
         i=i1
         return
      endif
!
      ii=(i1+i2)/2
!
      if(zxget.lt.x(ii)) then
         i2=ii-1
      else if(zxget.gt.x(ii+1)) then
         i1=ii+1
      else
! found
         i=ii
         return
      endif
!
      go to 10
!
      end subroutine zonfind_db


      end module bicub_db
      