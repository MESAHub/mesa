      subroutine varfrq(x,y,data,aa,sig,iy,ia,nn,istsbc,sigv)
c
c  calculates squared frequency from variational integral
c
c  distinguishes between radial and non-radial case, depending on
c  value of el.
c  note however that by setting ivarf = 3, the non-radial formulation
c  may be used for radial modes. this option is  n o t  recommended.
c  apparently it causes rather severe cancellation between gravitational
c  and dynamical term, in particular for low-order modes.
c
c  modified 19/3/1986 to allow use of tabulated coefficients in s/r derivk
c  and vintk.
c
c  12/8/87: version from RECKU and version from Boulder (1985) 
c  merged.
c
c  modified 13/8/87 to standardize output
c
c  modified 21/2/89 to systematize diagnostic output
c
c  modified 10/5/89 to test for same mesh etc. before using tabulated
c  coefficients in s/r derivk and vintk.
c
c  modified 14/12/93, to take into account the (1-x)**(mu-1)
c  behaviour of some of the integrands at a polytropic surface
c  of index mu.
c  Note: this has not been done yet with sufficient care, particularly
c  for mu .lt. 1 where we have an integrable singularity.
c
c  Modified 20/1/95 for resetting of treatment of turbulent pressure
c
c  Modified 2/7/95 to move ivarf from argument list to common /varcon/
c
c  Modified 18/2/98, to correct radial variational integral when
c  g/(g tilde) is .ne. 1.
c
c  ...................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      logical usegam
      dimension x(nn),y(iy,nn),data(*),aa(ia,nn)
      common/rhsdat/ el,ell,alb,els,el1
      common/worksp/ wrk(8,1)
      common/varcon/ ivarf,npvarf,kvarfc
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data xpfrst, xphalf, xplast,nnp, kvarfp /-1., -1., -1., 0, 0/
      data iopfdg /0/
c
      save
c
      idiag=0
c
      if(idiag.gt.0.and.iopfdg.eq.0) then
        open(52,file='ttt.varfrq.diag',status='unknown')
        open(53,file='ttt.varfrq.diag.u',status='unknown',
     *    form='unformatted')
        iopfdg=1
      end if
c
      iw=8
c
      if(kvarfc.eq.1.or.kvarfc.gt.3) kvarfc=2
c
c  test for using stored coefficients
c
      if(xpfrst.ne.x(1).or.xphalf.ne.x(nn/2).or.xplast.ne.x(nn)
     *  .or.nn.ne.nnp.or.kvarfp.ne.kvarfc) then
        kvarfi=kvarfc
        kvarfp=kvarfc
        xpfrst=x(1)
        xphalf=x(nn/2)
        xplast=x(nn)
        nnp=nn
      end if
c
c  find derivative of z1
c
      ii=1
      call derivk(x,y(1,1),wrk(1,1),ii,nn,iy,iw,kvarfi)
c
c  in non-radial case, reset ivarf from 3 to 1
c
      if(el.gt.1.e-6.and.ivarf.eq.3) ivarf=1
      if(el.le.1.e-6.and.ivarf.ne.3) go to 30
c
c  nonradial case
c  **************
c
c  limit on x, to avoid overflow in x**(-l-1)
c
      xlim=epsufl**(1.d0/(1.d0+el))
c
c  set first integrand
c
      do 10 n=2,nn
      xx=x(n)
c
c  test for x below limit
c
      if((el.gt.1.and.xx.lt.xlim).or.abs(y(1,n)).le.eprufl) then
c
c  zero wrk
c
        do 5 i=3,7
    5   wrk(i,n)=0
c
      else
c
c  test for case
c
        if(ivarf.eq.1) then
          d1=wrk(1,n)+(2*y(1,n)-y(2,n))/xx
          z2=y(2,n)
        else if(ivarf.eq.2) then
          d1=aa(2,n)*(y(1,n)-y(2,n)/(els*aa(1,n))-y(3,n))/xx
          z2=2*y(1,n)+xx*(wrk(1,n)-d1)
        else if(ivarf.eq.3) then
          d1=wrk(1,n)+2*y(1,n)/xx
          z2=y(2,n)
        else
	  write(istdou,97) ivarf
	  if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,97) ivarf
	  stop
        end if
c
        d2=d1-(aa(2,n)+aa(4,n))*y(1,n)/xx
c
        qxu=aa(1,n)*aa(5,n)*xx*xx
        wrk(3,n)=d2*xx**el*qxu
        wrk(4,n)=d1
        wrk(5,n)=d2
        wrk(6,n)=qxu
        wrk(7,n)=z2
c
      end if
c
   10 continue
c
c  zero integrals at centre
c
      do 18 i=2,6
   18 wrk(i,1)=0
      wrk(7,1)=(el+1)*y(1,1)
c
c  reset integrand at surface, for polytropic surface of index below 1
c
      an=data(7)
      if(an.gt.0.and.an.le.1) wrk(3,nn)=wrk(3,nn-1)
c
c  first integration
c
      ii=1
      call vintk(x,wrk(3,1),wrk(2,1),ii,nn,iw,iw,kvarfi)
c
c  test for resetting the contribution from the last mesh interval,
c  in the case of a singular surface
c
      if(an.gt.0) then
	wrk(2,nn)=wrk(2,nn-1)+(1.d0-x(nn-1))*wrk(3,nn-1)/an
      end if
c
c  reset kvarfi to use tabulated coefficients in subsequent integrations
c  or differentiations
c
      kvarfi=0
c
c  set other integrands
c
      elli=0.d0
      if(el.gt.1.e-6) elli=1.d0/ell
      wrk(8,1)=0
      xl1=0.d0
c
      do 20 n=1,nn
      z1=y(1,n)
      z2=wrk(7,n)
      xx=x(n)
c
c  test for x or y below limit
c
      if(el.gt.1.and.(xx.lt.xlim.or.abs(z1).lt.eprufl)) then
c
        do 18050 i=6,8
18050   wrk(i,n)=0
c
      else
c
        xd1=wrk(4,n)*xx
        qxu=wrk(6,n)
        vg=aa(2,n)
        if(vg.eq.0) vg=epsofl
        wrk(6,n)=qxu*(z1*z1+z2*z2*elli)
        if(n.gt.1) xl1=xx**(-1.d0-el)
        wrk(7,n)=wrk(2,n)*wrk(5,n)*xl1*qxu
        if(n.gt.1) then
          wrk(8,n)=(xd1*xd1/vg-2*xd1*z1+(aa(2,n)+aa(4,n))*z1*z1)*qxu*
     *      aa(1,n)
	end if
      end if
   20 continue
c
c  reset integrand at surface, for polytropic surface of index below 1
c
      if(an.gt.0.and.an.le.1) wrk(8,nn)=wrk(8,nn-1)
c
c  the remaining integrations
c
      ii=3
      call vintk(x,wrk(6,1),wrk(3,1),ii,nn,iw,iw,kvarfi)
c
c  test for resetting the contribution from the last mesh interval,
c  in the case of a singular surface
c
      if(an.gt.0) then
	wrk(5,nn)=wrk(5,nn-1)+(1.d0-x(nn-1))*wrk(8,nn-1)/an
      end if
c
c  set variational sigma
c
      el21=2*el+1
      akv=wrk(5,nn)-2*wrk(4,nn)/el21
      uz=aa(5,nn)*z1
c
c  set p prime on surface, depending on whether mode is
c  non-radial or radial
c
      if(el.lt.1.e-6) go to 22
      pprs=z2/els+y(3,nn)
      go to 25
c
   22 pprs=sig*z2
c
   25 aks=uz*(pprs-(uz-2*wrk(2,nn))/el21)
      sigv=(akv+aks)/wrk(3,nn)
      sv=aks/akv
      if(istdpr.gt.0) then
        if(ivarf.ne.3) then
          write(istdpr,90) ivarf
        else
          write(istdpr,95)
        end if
c
        write(istdpr,100) sigv,sv
      end if
      sigv=abs(sigv)
c
c  diagnostic output
c
      if(npvarf.gt.0.and.istdpr.gt.0) then
c
c  print integrands and integrals
c
        nd=max0(1,nn/npvarf)
        write(istdpr,120) 
        do 27 n=1,nn,nd
   27   write(istdpr,125) n,x(n),wrk(2,n),(wrk(i,n),i=6,8),
     *    (wrk(i,n),i=3,5)
c
      end if
      go to 70
c
c  -----------------------------------------------------------
c
c  radial case
c  ***********
c
c  set integrands. for test purposes positive and negative parts of j
c  are calculated separately.
c  set test for use of gamma1
c
   30 usegam=aa(3,1).gt.0.1
c
c  when using gamma1, differentiate numerically (using low order method)
c
      if(usegam) call derive(x,aa(3,1),wrk(3,1),nn,ia,iw,1,1)
      do 40 n=2,nn
      xx=x(n)
      qtx=aa(1,n)*xx*xx
      z1=y(1,n)
      qxu=aa(10,n)*qtx*aa(5,n)
c
      wrk(2,n)=z1*z1*qxu
c
      qxu=qxu*qtx
      vg=aa(2,n)
      if(vg.eq.0) vg=epsofl
c  test for use of gamma1
      if(usegam) then
c  use gamma1
        z1=z1/xx
        d1=wrk(1,n)-z1
        gm1=aa(3,n)
        dgm1=wrk(3,n)
c
c  to avoid underflow problems, set dgm1 to zero when very small
c
        if(abs(dgm1).lt.1.e-10) dgm1=0
        wrk(3,n)=(d1*d1/vg+(3*gm1-4)*z1*z1)*qxu
        wrk(4,n)=-3*qxu*xx*dgm1*z1*z1/(vg*gm1)
c
c  for diagnostics, store dz
c
        wrk(8,n)=wrk(1,n)
      else
        z1=2*z1/xx
        d1=wrk(1,n)+z1
        wrk(3,n)=d1*d1*qxu/vg
        wrk(4,n)=z1*z1*qxu
      end if
   40 continue
c  central values
      do 45 i=2,7
   45 wrk(i,1)=0
c
c  integrate
c
      ii=3
      call vintk(x,wrk(2,1),wrk(5,1),ii,nn,iw,iw,kvarfi)
c
c  reset kvarfi to use tabulated coefficients in subsequent integrations
c  or differentiations
c
      kvarfi=0
c
c  set variational sigma**2
c  test for use of gamma1
c
      if(usegam) then
c  use gamma1
        z1=y(1,nn)
        ajs=z1*((3/vg-1)*z1+sig*y(2,nn))*aa(5,nn)
        dj=ajs+wrk(6,nn)+wrk(7,nn)
        ajs=ajs/dj
        sv=wrk(7,nn)/dj
        sigv=dj/wrk(5,nn)
        if(istdpr.gt.0) write(istdpr,130) sigv,sv,ajs
c
      else
c
        dj=wrk(6,nn)-wrk(7,nn)
c  for isothermal atmosphere include surface term
        if(istsbc.ne.0) then
          z1=y(1,nn)
          ajs=z1*(sig*y(2,nn)-z1)*aa(5,nn)
          dj=dj+ajs
          ajs=ajs/dj
        end if
c
        sigv=dj/wrk(5,nn)
c  output. relative contribution from negative part
        sv=wrk(7,nn)/dj
        if(istdpr.gt.0) write(istdpr,110) sigv,sv
c  isothermal atmosphere?
        if(istsbc.eq.1.and.istdpr.gt.0) write(istdpr,115) ajs
c
      end if
c
c  diagnostic output to file
c
      if(idiag.gt.0) then
        write(52,150) sig,sigv,sigv/sig-1,
     *    (n,x(n),y(1,n),wrk(8,n),wrk(2,n),wrk(3,n),wrk(5,n),wrk(6,n),
     *    n=1,nn)
      end if
c
c  diagnostic output
c
      if(npvarf.gt.0.and.istdpr.gt.0) then
c
c  print integrands and integrals
c
        nd=max0(1,nn/npvarf)
        write(istdpr,122) 
        do 60 n=1,nn,nd
   60   write(istdpr,125) n,x(n),(wrk(i,n),i=2,4),
     *    (wrk(i,n),i=5,7)
c
      end if
c
c  --------------------------------------------------
c
c  end calculation
c
   70 sigv=abs(sigv)
c
c  diagnostic output to file
c
      if(idiag.gt.0) then
	write(53) 9,nn,(x(n),(wrk(i,n),i=1,8),n=1,nn)
      end if
c
      return
   90 format(//' variational integral, formulation',i3)
   95 format(//' variational integral. non-radial formulation',
     *  ' for radial mode.')
   97 format(/' ***** Error in s/r varfrq. ivarf =',i5,' not allowed.')
  100 format(/' sigv =',1pe15.7,'.  as/av =',e13.5)
  110 format(//' radial variational integral.'//'  sigv =',1pe15.7,
     *  '.   (neg. part)/total =',e13.5)
  115 format(/' isothermal atmosphere. (surface term)/total =',
     *  1pe13.5)
  120 format(///' non-radial case in varfrq.'/' n, x,',
     *  ' integral from first integration, integrands(1-3)',
     *  ' and integrals(1-3) from second integration:'/)
  122 format(///' radial case in varfrq.'/' n, x,',
     *  ' integrands(1-3)',
     *  ' and integrals(1-3) from second integration:'/)
  125 format(i5,0pf10.5,1p8e13.5)
  130 format(//' radial variational integral.'//'  sigv =',1pe15.7,
     *  /' (der.of gamma1 part)/total =',e13.5/
     *  ' (surface term)/total =',e13.5)
  150 format('# eigenvalue  sigma**2 =',1pe15.7/
     *       '# variational sigma**2 =',e15.7/
     *       '# relative difference  =',e13.5/'#'/
     *       '# n, x, y, dy, integrand(1 - 2), integral(1 - 2):'/'#'/
     *        (i5,0pf12.7,1p6e13.5))
      end
