      subroutine gm1ker(idsgkr,x,y,aa,elarg,sig,iy,ia,nn,npgmkr)
c
c  calculates and outputs kernel gmk(1,n) for relative frequency
c  change caused by change in gamma1. thus
c
c    delta omega/omega = integral(gmk*delta gamma1 *dx)
c
c  mode of calculation depends on ivarf (for p modes use ivarf = 1,
c  for g modes use ivarf = 2), in the same way as the variational
c  integral calculation in s/r varfrq
c
c  Modified 21/7/94, to include unit number as argument
c
c  Modified 2/7/95 to move ivarf from argument list to common /varcon/
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension x(nn),y(iy,nn),aa(ia,nn)
      common/rhsdat/ el,ell,alb,els,el1
      common/worksp/ wrk(8,1)
      common/csumma/ cs(50)
      common/varcon/ ivarf,npvarf,k
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      iw=8
      if(k.lt.2.or.k.gt.3) k=2
c
c  find derivative of z1
c
      ii=1
      call derivk(x,y(1,1),wrk(1,1),ii,nn,iy,iw,k)
c
c  set elli to 1/(l(l+1)) in non-radial case, or zero for radial
c  modes
c
      elli=0
      if(el.gt.1.e-6) elli=1./ell
c
c  limit on x, to avoid overflow in x**(-l-1)
c
      xlim=epsufl**(1./(1.+el))
c
c  set first integrand
c
      wrk(1,1)=0
      wrk(2,1)=0
      do 10 n=2,nn
      xx=x(n)
c
c  test for x below limit
c
      if((el.le.1.or.xx.ge.xlim).and.abs(y(1,n)).gt.eprufl
     *  .and.aa(2,n).ne.0) go to 8
c
c  zero wrk
c
      wrk(1,n)=0
      wrk(2,n)=0
      go to 10
c
c  test for case
c
    8 if(ivarf.eq.2) go to 12
c
      d1=wrk(1,n)+(2*y(1,n)-y(2,n))/xx
      z2=y(2,n)
c
      go to 15
   12 d1=aa(2,n)*(y(1,n)-y(2,n)/(els*aa(1,n))-y(3,n))/xx
      z2=2*y(1,n)+xx*(wrk(1,n)-d1)
c
   15 qxu=aa(1,n)*aa(5,n)*xx*xx
      z1=y(1,n)
      wrk(2,n)=qxu*(z1*z1+elli*z2*z2)
      xd1=xx*d1
      wrk(1,n)=xd1*xd1*qxu*aa(1,n)/(aa(2,n)*aa(3,n))
c
   10 continue
c
c  zero integral at centre
c
      wrk(3,1)=0
      wrk(7,1)=(el+1)*y(1,1)
c
c   integration of denominator
c
      ii=1
      call vintk(x,wrk(2,1),wrk(3,1),ii,nn,iw,iw,k)
c
c  set scale factor
c
      gmkfct=1./(2*sig*wrk(3,nn))
c
c  set kernel
c
      do 20 n=1,nn
   20 wrk(1,n)=gmkfct*wrk(1,n)
c
c  output to file
c
      write(idsgkr) cs,nn,(x(n),wrk(1,n),n=1,nn)
c
c  test for printed output
c
      if(npgmkr.le.1) return
c
      nd=max0(1,(nn-1)/(npgmkr-1))
      if(istdpr.gt.0) write(istdpr,110) (n,x(n),wrk(1,n),n=1,nn,nd)
      return
c
  110 format(//' n, x, gamma1 kernel:'//(i5,0pf10.6,1pe13.5))
      end
