      subroutine lsqpol(x,y,nn,ix,iy,kk,a,rmsres,idiag)
c
c  makes least-squares fit to polynomial of degree kk
c
c       y = a(1) + a(2)*x + ... + a(kk+1)*x**kk
c
c  to data points (x(1,i), y(1,i)), i = 1, ..., nn.
c  ix and iy are first dimensions of x and y.
c  returns rms residual around fit in rmsres.
c  if idiag .ge. 1 prints x, y, fitted curve and residual
c
c  ......................................................................
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 17/5/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision(a-h,o-z)
      dimension x(ix,1),y(iy,1),a(1),w(10000),aa(20,20),w1(10000)
c
      do 10 n=1,nn
   10 w(n)=1
c
c  to avoid over and underflow problems rescale x to
c  be between -1 and 1.
c
      ixmax=isamax(nn,x,ix)
      xscl=x(1,ixmax)
      xscli=1.d0/xscl
c
      do 12 n=1,nn
   12 w1(n)=xscli*x(1,n)
c
      kk1=kk+1
      kk2=kk+kk1
      kkb=kk1+1
c
      do 30 k=1,kk2
      sk=ssum(nn,w,1)
      if(k.le.kk1) aa(k,kkb)=sdot(nn,w,1,y,iy)
      i1=max0(1,k-kk1+1)
      i2=min0(k,kk1)
      do 15 i=i1,i2
   15 aa(i,k+1-i)=sk
      if(k.eq.kk2) go to 30
      do 20 n=1,nn
   20 w(n)=w1(n)*w(n)
   30 continue
c
c  test for output of equations
c
      if(idiag.lt.2) go to 36
c
      write(6,30090)
      do 35 i=1,kk1
   35 write(6,30091) (aa(i,j),j=1,kkb)
30090 format(///' aa:'/)
30091 format(1p10e13.5)
c
   36 call leq(aa,aa(1,kkb),kk1,1,20,20,det)
c
c
c  rescale coefficients
c
      scl=1
      do 38 k=1,kk1
      a(k)=scl*aa(k,kkb)
   38 scl=scl*xscli
c
      write(6,100) (a(k),k=1,kk1)
c
c  set residuals and their rms
c
      do 40 n=1,nn
      w(n)=1
   40 w1(n)=0
c
      do 50 k=1,kk1
      do 50 n=1,nn
      w1(n)=w1(n)+a(k)*w(n)
   50 w(n)=w(n)*x(1,n)
c
      do 55 n=1,nn
   55 w(n)=y(1,n)-w1(n)
c
      rmsres=sqrt(sdot(nn,w,1,w,1)/nn)
      write(6,120) rmsres
c
c  test for output of fit
c
      if(idiag.lt.1) return
c
      write(6,130) (n,x(1,n),y(1,n),w1(n),w(n),n=1,nn)
      return
  100 format(//' coefficients found in s/r lsqpol:'/
     *  (1p10e13.5))
  120 format(//' rms residual =',1pe13.5)
  130 format(//' n, x, y, fit, residual:'//(i4,1p4e13.5))
      end
