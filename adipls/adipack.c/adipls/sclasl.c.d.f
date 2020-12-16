      subroutine sclasl(x,y,iy,nn,ii,nev1,nw1,el,iplneq)
c
c  multiply solution by (x/x(nev1))**(el-1) in interior
c  or, in plane-parallel case, by exp(el*(x-x(nev1))
c  ****************************************************
c
c  original version: 14/2/89
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
c  Modified 9/7/05, changing only name to avoid conflict with
c  evolution code.
c
      implicit double precision (a-h,o-z)
      dimension x(1),y(iy,1)
      xev1=x(nev1)
c
      if(iplneq.ne.1) then
c
c  spherical case
c
        icxm=0
        if(el.gt.15) icxm=1
        do 10 n=nw1,nev1
        xx=x(n)/xev1
        do 10 i=1,ii
        yy=xlmult(y(i,n),xx,el-1,icxm,ierr)
   10   y(i,n)=yy
c
      else
c
c  plane case
c
        do 20 n=nw1,nev1
        yy=exp( el*(x(n)-xev1))
        do 20 i=1,ii
   20   y(i,n)=yy*y(i,n)
c
      end if
c
      return
      end
