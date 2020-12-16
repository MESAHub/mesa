      subroutine derive(x,y,dydx,nn,iy,idy,ipy,ipdy)
c     subroutine deriv2(x,y,dydx,nn,iy,idy,ipy,ipdy)
c
c
c     derive sets dydx(1,n) = first  derivative of y(1,n) w.r.t. x(n)
c     deriv2 sets dydx(1,n) = second derivative of y(1,n) w.r.t. x(n)
c                             n=1,nn
c
c     iy  is the first dimension of y    in the calling programme
c     idy is the first dimension of dydx in the calling programme
c
c     second order accuracy differences are used at interior points
c     at end points third order differences are used for first
c     derivative.  second derivative is obtained by quadratic
c     interpolation from the interior
c
c
c  revised on 6/9 1984 to include scaling by interval length,
c  to avoid underflow and consequent divide errors.
c
c                      ****************************************
c
c     notes on precision parameters:
c
c   The arguments ipy and ipdy are kept for consistency with previous
c   versions of routine. However they have no effect.
c
c  .............................................................................
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision(a-h,o-z)
      logical second
      dimension x(*), y(*),dydx(*)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data epsil/1.d-8/
c
c
      second=.false.
      go to 10
c
      entry deriv2(x,y,dydx,nn,iy,idy,ipy,ipdy)
      second=.true.
c
c
   10 ky=iy*ipy
      kdy=idy*ipdy
      if(ipdy.eq.1) go to 30
c
c     if dydx is real*8, set insignificant half to zero
      j=1-kdy
   20 do 21 n=1,nn
      j=j+kdy
   21 dydx(j)=0.d0
c
   30 n1=nn-1
      n=0
      i=1
      k=1+ky
      hn =x(2)-x(1)
      xn= abs(x(1))
      hna=abs(hn)
      e=y(k)-y(1)
      if(second) go to 37
c
c     first derivative at end points
      hn1=x(3)-x(2)
      h3=x(4)-x(3)
      nx1=0
      xn1= abs(x(4))
c
      hna1=abs(hn1)
      ha3=abs(h3)
c
c  rescale differences
c
      if(hn.eq.0) go to 361
      hnin=1.d0/hn
      hn1=hnin*hn1
      h3=hnin*h3
c
   31 hn2=1.d0+hn1
      hn3=hn2+h3
c     test size of intervals
      xxa=(xn+xn1)*epsil
      if(hna.lt.xxa.or.hna1.lt.xxa.or.ha3.lt.xxa) go to 361
      c=hn2*hn3*(hn2*((hn1+h3)*(hn3+1.d0)-2.d0*h3-hn1*hn3)-
     .  (h3*h3+hn1*hn3))
      if(second) go to 34
   32 d=hn2*hn2
      dysc=hnin
      f=1.d0
      b=f*h3
      f=f*hn1
      hn1=hn3*hn3
      a=d*hn1*h3
      b=-hn1*(b+f)
      d=d*f
      if(n) 35,35,33
   33 c=-c
      go to 36
   34 a=-hn2*hn3*h3*(hn2+hn3)
      b=hn3*(hn1+h3)*(hn3+1.d0)
      d=-hn1*hn2*(1.d0+hn2)
      c=0.5d0*c
      dysc=hnin*hnin
      if(n) 35,35,36
   35 j=1
      l=k+ky*2
   36 dydx(i)=dysc*(a*e+b*(y(k+ky)-y(j))+d*(y(l)-y(j)))/c
      go to 362
c     zero interval. write diagnostics
  361 nj1=nx1+1
      nj2=nx1+4
      if(istdpr.gt.0) 
     *  write(istdpr,1000) nj1,nj2,(x(nx1+j),j=1,4)
      dydx(i)=0.d0
c
  362 if(n.ne.0) return
c
c     derivative at interior points
   37 do 42 n=2,n1
      i=i+kdy
      j=k
      k=k+ky
      xn1=xn
      xn= abs(x(n))
      xxa=(xn1+xn)*epsil
      d=e
      e=y(k)-y(j)
      hn1=hn
      hn=x(n+1)-x(n)
      hna1=hna
      hna=abs(hn)
c     test size of intervals
      if(hna.ge.xxa.and.hna1.ge.xxa) go to 371
      dydx(i)=dydx(i-kdy)
      nj1=n-1
      nj2=n+1
      if(istdpr.gt.0) 
     *  write(istdpr,1000) nj1,nj2,x(nj1),x(n),x(nj2)
      go to 42
c
c  rescale differences
c
  371 hnin=1.d0/hn
      hn1=hnin*hn1
c
      c=hn1*(1.d0+hn1)
      if(second) go to 39
   38 a=hn1*hn1
      b=1.d0
      dysc=hnin
      go to 40
   39 a=hn1
      b=-1.d0
      c=0.5d0*c
      dysc=hnin*hnin
   40 dydx(i)=dysc*(a*e+b*d)/c
   42 continue
c
      h3=x(nn-2)-x(nn-3)
c
      ha3=abs(h3)
      h3=hnin*h3
c
      xn= abs(x(nn))
      xn1= abs(x(nn-3))
      nx1=nn-4
      if(second) go to 50
c
c     storage indices for first derivative at last point
      i=i+kdy
      j=k
      k=k-ky*3
      l=k
      e=-e
      go to 31
c
c     second derivative at end points
   50 j=i
      k=i-kdy
      l=k-kdy
      i=i+kdy
   51 a=1.d0+hn1
      b=a+h3
      c=hn1+h3
      xxa=(xn1+xn)*epsil
c     test size of intervals
      if(hna.ge.xxa.and.hna1.ge.xxa.and.ha3.ge.xxa) go to 52
51100 nj1=nx1+1
      nj2=nx1+4
      if(istdpr.gt.0) 
     *  write(istdpr,1000) nj1,nj2,(x(nx1+j),j=1,4)
      dydx(i)=0.d0
      go to 53
   52 dydx(i)=(a*b/(hn1*c))*dydx(j)-(b/(hn1*h3))*dydx(k)
     .       +(a/(c*h3))*dydx(l)
   53 if(i.eq.1) return
      i=1
      j=i+kdy
      k=j+kdy
      l=k+kdy
      hn=x(2)-x(1)
      hn1=x(3)-x(2)
      h3=x(4)-x(3)
      xn= abs(x(1))
      xn1= abs(x(4))
      nx1=0
      hna=abs(hn)
      hna1=abs(hn1)
      ha3=abs(h3)
c
c  rescale differences
c
      if(hn.eq.0) go to 51100
      hnin=1.d0/hn
      hn1=hnin*hn1
      h3=hnin*h3
c
      go to 51
c
 1000 format(' **** from derive: degeneracy among x(',i5,') - x(',
     .  i5,') = ',1p4e16.8)
      end
