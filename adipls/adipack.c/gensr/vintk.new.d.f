      subroutine vintk(x,a,y,ii,nn,idima,idimy,karg)
c
c      2k-1 -th order integration routine
c      ***********************************
c
c  sets integral of a(i,n) with respect to x into y(i,n),i=1,ii,n=1, nn.
c  the values of y(i,1), i = 1, ii, must be supplied.
c  idima and idimy are first dimensions of a and y respectively.
c  if x is found to be non-monotonic ii is set to 0.
c
c  if vintk is called with karg .le. 0, it is assumed that a table
c  of integration coefficients has been set up.
c
c  vintk uses as work space
c      wrk(2k*2k)
c  the size of wrk set in vintk is 100, which is sufficient, say, for
c  9-th order integration in 20 dependent variables. if more
c  work space is needed it must be set in common/wrklir/
c
c  note that this version of vintk also uses storage of size
c
c         dc(2*k*(nn-1))
c
c  in common /cintst/ to store integration coefficients. the
c  default size of 500 probably in general has to be increased in
c  the calling programme.
c
c  original version: 19/3/1986
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90.
c
c  Modified 26/4/99, fixing redefinition of k when nn is too small
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision(a-h,o-z)
      dimension x(1),a(idima,1),y(idimy,1)
      common/wrklir/ wrk(100)
      common/cintst/ dc(500)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save /cintst/,k,k2
c
c  test for setting up coefficients
c
      if(karg.le.0) go to 60
c
      k=karg
c
c  check order
c
      k2=2*k
      if(nn.lt.k2) then
        if(istdpr.gt.0) write(istdpr,100) k2,nn
        k=nn/2
        k2=2*k
      end if
c
c  first dimension of wrk, storage information
c
      iw=k2
      mst=iw*(iw-1)+2
c
c  total range
c
      diff=x(nn)-x(1)
      if(diff.eq.0) go to 95
c
c  start loop for setting coefficients
c
   20 jdc=0
c
      do 50 n=2,nn
      n1=n-1
      xtr=x(n1)
c
c  determine i
c
      i=min0(max0(0,n1-k),nn-k2)
c  length of interval
      dlx=x(i+k2)-x(i+1)
      if(diff*dlx.le.0) go to 97
c
c  set equations
c
   32 dlxi=1.d0/dlx
c
      m=0
      do 35 l=1,k2
      l1=i+l
      xl=(x(l1)-xtr)*dlxi
      aa=1.d0
      m=m+1
      wrk(m)=aa
      do 35 ir=2,k2
      m=m+1
      aa=xl*aa
   35 wrk(m)=aa
c
c  set right hand sides
c
      j1=jdc
      xl=(x(n)-xtr)*dlxi
      aa=1
c
      do 40 j=1,k2
      aa=aa*xl
      j1=j1+1
   40 dc(j1)=dlx*aa/j
c
c  solve equations
c
      call leq(wrk,dc(jdc+1),k2,1,iw,1,det)
c
      if(det.eq.0) go to 98
c
   50 jdc=jdc+k2
c
c  end setting coefficients
c
      if(istdpr.gt.0) write(istdpr,110) jdc
c
c            *****************************************
c
c  set integrals
c
   60 jdc=0
c
      do 70 n=2,nn
      n1=n-1
c
c  determine i
c
      i=min0(max0(0,n1-k),nn-k2)
c
      do 65 is=1,ii
      j1=jdc
      sum=0
      do 62 l=1,k2
      j1=j1+1
   62 sum=sum+dc(j1)*a(is,i+l)
c
   65 y(is,n)=y(is,n1)+sum
c
   70 jdc=jdc+k2
c
      return
c
c  diagnostics
c
   95 if(istdpr.gt.0) write(istdpr,120) x(1)
      ii=0
      return
   97 i1=i+1
      ik2=i+k2
      if(istdpr.gt.0) write(istdpr,130) x(1),x(nn),i1,x(i1),ik2,x(ik2)
      ii=0
      return
   98 i1=i+1
      ik2=i+k2
      if(istdpr.gt.0) write(istdpr,140) (l,x(l),l=i1,ik2)
      ii=0
      return
  100 format(1x,10('*'),' 2k+1 =',i4,' is greater than nn =',i4,
     *  ' in vintk.'/11x,'k has been reset to nn/2')
  110 format(//' storage needed in common/cintst/ in vintk:',
     *  i7)
  120 format(//1x,10('*'),' range is zero in vintk. x(1)=',
     *  1pe15.5,' = x(nn)'//)
  130 format(//1x,10('*'),' independent variable is not monotonic',
     *  ' in vintk'//' x(1)=',1pe15.5,' x(nn)=',e15.5,
     *  2(' x(',i4,')=',e15.5)//)
  140 format(//1x,10('*'),' points coincide in vintk'/
     *  (5(' x(',i4,')=',1pe13.4)))
      end
