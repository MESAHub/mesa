      subroutine gravpo(x,y,data,aa,sig,iy,ia,nn,iasn,dsig,yp,iyp,icry)
c
c  calculates perturbation in gravitational potential and correction
c  to frequency
c
c  modified 13/8/87 to standardize output
c
c  modified 14/2/89 to allow step iasn .gt. 1 in model and solution
c
c  modified 12/7/90, to test for singular surface and take appropriate
c  action when iasm .gt. 1.
c
c  Note: 11/5/89 adipls was modified such that routine may be called 
c  also for truncated model. In this case it is assumed that 
c  contribution from missing part of model is zero.
c
c  Modified 20/1/95 for resetting of treatment of turbulent pressure
c  ....................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension x(nn),y(iy,nn),yp(iyp,nn),data(8),aa(ia,nn)
      common/worksp/ wrk(8,1)
      common/xarr1/ x1(1)
      common/xarr2/ x2(1)
      common/rhsdat/ el,ell,alb,els,el1
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data nwrt/0/
      data ncount /-1/
c
c  test and set actual number of meshpoints, depending on whether model
c  extends to centre and whether surface is singular
c
      if(ncount.ge.0) then
        ncount=ncount+1
        if(istdpr.gt.0) write(istdpr,*) 
     *    'Enter gravpo with ncount = ',ncount,' iasn =',iasn
      end if
      if(x(1).eq.0) then
        nc=2
      else
        nc=1
      end if
c
      if(data(7).lt.0) then
	nsr=nn
      else
	nsr=nn-1
      end if
c
      if(mod(nsr-nc,iasn).ne.0) then
        write(istdou,100) nn, nsr, nc, iasn
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,100) 
     *    nn, nsr, nc, iasn
        icry = -1
        return
      end if
c
c  set numbers of meshpoints for a possible submesh, possibly excluding
c  surface if singular
c
      if(iasn.eq.1) then
	nn1=nn
	nn2=nn
      else if(data(7).lt.0) then
        nn1=(nn-nc)/iasn+nc
	nn2=nn1
      else 
        nn1=(nsr-nc)/iasn+nc
	nn2=nn1+1
      end if
c
      iw=8
      call zero(wrk,8*nn2)
c
c  set first integrand
c
      icxm=0
      if(el.gt.10) icxm=1
      do 10 n1=2,nn1
      n=iasn*(n1-nc)+nc
      xx=x(n)
      d=(aa(4,n)*y(1,n)+aa(2,n)/(aa(1,n)*els)*y(2,n))*aa(1,n)*aa(5,n)*
     *  aa(10,n)/xx
      x1(n1)=xx
      wrk(1,n1)=d
    5 wrk(3,n1)=xlmult(d,xx,el+2,icxm,ierr)
   10 continue
      x1(1)=x(1)
      wrk(1,1)=0
      wrk(3,1)=0
c
c  test for setting x1(nn2) in case of submesh and singular 
c  surface
c
      if(nn1.lt.nn2) then
        x1(nn2)=1
      end if
c
c  first integration
c
      call vinta(x1,wrk(3,1),wrk(2,1),nn2,iw,iw)
c
c  set other integrands
c
      elli=1.d0/ell
      nnp1=nn+1
      do 20 n1=2,nn1
      n=iasn*(n1-nc)+nc
      xx=x(n)
      y1=y(1,n)*xx
      y2=y(2,n)*xx
      if(abs(y1).gt.eprufl) go to 15
      y1=0
      y2=0
   15 wrk(6,n1)=aa(1,n)*aa(5,n)*aa(10,n)*(y1*y1+y2*y2*elli)
      ww1=xlmult(wrk(1,n1),xx,1-el,icxm,ierr)
      wrk(7,n1)=wrk(2,n1)*ww1
      ns1=nn2+1-n1
      wrk(8,ns1)=-ww1
   20 x2(ns1)=xx
c
c  test for setting x2(1) in case of submesh and singular 
c  surface
c
      if(nn1.lt.nn2) then
        x2(1)=1
      end if
c
      wrk(6,1)=0
      wrk(7,1)=0
      wrk(8,nn2)=0
      x2(nn2)=0
c
      if(nwrt.eq.1.and.istdpr.gt.0) then
        write(istdpr,110)
        do 22 n1=1,nn1
        n=iasn*(n1-nc)+nc
   22   write(istdpr,115) n,x(n),y(1,n),y(2,n),(wrk(i,n1),
     *    i=1,3),(wrk(i,n1),i=6,7),wrk(8,nn1+1-n1)
      end if
c
c  the remaining integrations
c
      do 25 k=3,4
   25 call vinta(x1,wrk(k+3,1),wrk(k,1),nn2,iw,iw)
      call vinta(x2,wrk(8,1),wrk(5,1),nn2,iw,iw)
c  set correction to sigma
      el21=2*el+1
      akv=wrk(4,nn2)
      uz=aa(5,nn)*y(1,nn)
      aks=uz*uz-2*uz*wrk(2,nn1)
      dsig=-(2*akv+aks)/(wrk(3,nn1)*el21)
c  set y3
      yp(1,1)=0
      do 30 n1=2,nn2
      n=min0(nn,iasn*(n1-nc)+nc)
      xx=x1(n1)
      ww1=xlmult(uz+wrk(5,nn2+1-n1),xx,el-1,icxm,ierr)
      ww2=xlmult(wrk(2,n1),xx,-el-2,icxm,ierr)
      yp(1,n)=(ww1+ww2)/(el21*aa(1,n))
      if(nwrt.eq.1.and.istdpr.gt.0) then
        write(istdpr,120) n,xx,uz,wrk(5,nn1+1-n1),
     *    wrk(2,n1),ww1,ww2,yp(1,n)
      end if
   30 continue
      if(ncount.ge.0) then
        do 40 n1=2,nn1
        n=iasn*(n1-nc)+nc
        write(81,'(4i5,f10.5,1p16e13.5)') 
     *     ncount,iasn,n1,n,x(n),y(1,n),y(2,n),
     *    (aa(i,n),i=1,5),aa(10,n),(wrk(i,n1),i=1,8)
   40   continue
      end if
      nwrt=0
      return
  100 format(/' **** error in s/r gravpo.',
     *  ' number of steps is not integral multiplum of iasn'/
     *  '            nn, nsr, nc, iasn =',4i6)
  110 format(//' output from gravpo:'/)
  115 format(i5,0pf10.6,1p8e13.5)
  120 format(i5,f10.6,1p6e13.5)
      end
