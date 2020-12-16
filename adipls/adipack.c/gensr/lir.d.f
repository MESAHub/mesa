      subroutine   lir(z,zi,y,yi,ii,id,nt,l,inter)
c     subroutine  lir1(z,zi,y,yi,ii,id,nt,l,inter)
c
c
c
c                interpolation/extrapolation routine
c
c
c     for a such that z=zi(a),  sets y(i)=yi(i,a), i=1,ii
c
c     zi(n),yi(i,n) must be supplied for n=1,nt and i=1,ii
c     id is first dimension of yi
c
c     inter is set to 1 for interpolation and 0 for extrapolation
c     inter is returned as -1 in case of errors
c
c     if l.le.1, scan to find the zi(n) which immediately bound z
c                starts at n=1
c     if l.gt.1, scan starts from value of n from previous call of lir
c
c
c     lir use cubic interpolation/extrapolation unless nt.lt.4
c     lir1 use linear interpolation/extrapolation
c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
c
      implicit double precision(a-h,o-z)
      dimension zi(1),y(1),yi(1),a(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data n/-1/
c
      il=0
      go to 1
      entry lir1(z,zi,y,yi,ii,id,nt,l,inter)
      il=1
    1 continue
      ir=1
c
c     check nt and reset il if necessary
      if(nt.lt.2) go to 101
      if(nt.lt.4) il=1
c
c     addressing constants
      inter=1
      ir1=ir-1
      ird=ir*id
      iir=(ii-1)*ir+1
      j=(nt-1)*ir+1
      diff=zi(j)-zi(1)
c
c     set index for start of search
      n=(n-2)*ir+1
      if(l.le.1.or.n.lt.1) n=1
c
c     determine position of z within zi
    2 if(n.gt.j) go to 8
      if(diff) 4,102,3
    3 if(zi(n)-z) 5,6,9
    4 if(zi(n)-z) 9,6,5
    5 n=n+ir
      go to 2
c
c     set y when z lies on a mesh point
    6 j=(n-1)*id
      do 7 i=1,iir
      y(i)=yi(i+j)
    7 if(y(i).eq.0.d0) y(i+ir1)=0.d0
      go to 30
c
c     control when z does not lie on a mesh point
    8 inter=0
    9 if(n.le.1) inter=0
      if(il.eq.1) go to 20
c
c     cubic interpolation/extrapolation
c
c     pivotal point (m) and point (k) closest to z
   10 m=n
      k=3
      if(n.gt.1+ir) go to 11
      m=1+ir+ir
      k=n
   11 if(n.lt.j) go to 12
      m=j-ir
      k=4
c
c
c     weighting factors
   12 y1=zi(m-ir*2)
      y2=zi(m-ir)
      y3=zi(m)
      y4=zi(m+ir)
c
      z1=z-y1
      z2=z-y2
      z3=z-y3
      z4=z-y4
c
   13 z12=z1*z2
      z34=z3*z4
c
   14 a(1)=z2*z34/((y1-y2)*(y1-y3)*(y1-y4))
      a(2)=z1*z34/((y2-y1)*(y2-y3)*(y2-y4))
      a(3)=z12*z4/((y3-y1)*(y3-y2)*(y3-y4))
      a(4)=z12*z3/((y4-y1)*(y4-y2)*(y4-y3))
c
c     correct a(k)
   15 diff=a(1)+a(2)+a(3)+a(4)
      a(k)=(1.d0+a(k))-diff
c
c     compute y
   16 m=(m-1)/ir-3
      m=m*ird
      do 18 i=1,iir
      k=i+m
      yy=0.d0
      do 17 j=1,4
      k=k+ird
      diff=yi(k)
   17 yy=yy+a(j)*diff
      y(i)=yy
   18 if(y(i).eq.0.d0) y(i+ir1)=0.d0
      go to 30
c
c     linear interpolation/extrapolation
   20 if(n.eq.1) n=1+ir
      if(n.gt.j) n=j
      z1=zi(n)
      y1=(z1-z)/(z1-zi(n-ir))
      y2=1.0-y1
      j=(n-1)*id
      m=j-ird
      do 21 i=1,iir,ir
      y(i)=y1*yi(i+m)+y2*yi(i+j)
   21 if(y(i).eq.0.d0) y(i+ir1)=0.d0
c
c     reset n
   30 n=(n+ir-1)/ir
      return
c
c
c     diagnostics
  101 if(istdpr.gt.0) write(istdpr,1001) nt
      inter=-1
      return
  102 if(istdpr.gt.0) write(istdpr,1002) zi(1),nt,zi(j)
      inter=-1
      return
c
 1001 format(/1x,10('*'),5x,'there are fewer than two data points in',
     *      ' lir     nt =',i4,5x,10('*')/)
 1002 format(/1x,10('*'),5x,'extreme values of independent variable',
     *      ' equal in lir',5x,10('*')/16x,'zi(1) =',1pe13.5,',   ',
     *       'zi(',i4,') =',1pe13.5/)
c
      end
