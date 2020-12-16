      program main
c
c  redistributes mesh in physical model
c
c  modified 9/1/1985 to allow use of central expansion to reset
c  points in new model interior to first non-zero point in old
c  model. old expansion may still be used by setting ioldex = 1.
c
c  modified 17/5/1985. now interpolates in x**(-2)*aa(i,n)
c  for i = 2 and 4 (these tend to a constant at x = 0),
c  and sets aa(2,n) interior to first mesh point in original
c  model from expansion of p and other aa-s.
c
c  modified 21/5/1985 in treatment of convection zone boundary
c  (see s/r stcvzb)
c
c  modified 21/5/1985 to include contribution in stretching
c  from superadiabatic gradient, to ensure adequate resolution very
c  near surface.
c
c  modified 29/3 - 1/4 1989 to include various fixes for stretching
c  polytropic model
c
c  modified 30/5/92, to include term in change in buoyancy frequency,
c  and smoothing with running mean.
c
c  modified 25/7/92, to smooth A4 near centre, particularly
c  near hydrogen exhaustion
c
c  modified 11/11/93, setting proper values for expansion at
c  singular surface
c
c  modified 12/11/93, introducing csurf to increase distance
c  of first point from surface for singular surface
c
c  modified 3/12/93, to allow alternative interpolation
c  (with Numerical Recipies spline routines)
c
c  modified 9/12/93, to use explicitly the logarithmic form
c  of the stretch function near a singular surfac, when using
c  alternative interpolation
c
c  modified 10/10/95, allowing option of icvznb .gt. 10, to insert
c  a double point at base of convection zone, for Richardson
c  extrapolation
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      parameter (nnmax = 10001)
      logical singsf
      dimension x(nnmax),aa(6,nnmax),xsi(nnmax),xn(nnmax),an(6,nnmax),
     *  data(8),yspl(nnmax),alptil(nnmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c..      namelist/exec/ nn,icnmsh,
c..     *  icase,icvzbn,nsmth,
c..     *  cg,cx,ca,cdgr,cddgr,cdg,
c..     *  nout,cn,irsu,unew,iresa4,csurf,
c..     *  nmodel,kmodel,itsaml,ioldex
c  defaults in exec
c  nn: number of points in new mesh
      nn=601
      icnmsh=0
c  icase: controls redistribution.
c  icase is of the form icase = 100*iredst + icaser
c  Here iredst controls interpolation:
c    iredst = 0: old formulation, using lir throughout
c    iredst = 1: Use Num. Rec. routine to set x mesh,
c    lir to interpolate A
c  icaser controls redistribution parameters:
c  for icaser = 1, 2, 3, 11 or 12 set standard parameters.
c     icaser = 1, 2, 3 give old parameters, as used up until June 1988.
c     icaser = 11 gives new, optimized p mode parameters
c     icaser = 1  corresponds to p modes, old parameters
c     icaser = 2  corresponds to g modes, old parameters
c     icaser = 3  corresponds to intermediate modes, old parameters.
c     icaser = 11 corresponds to p modes, new parameters.
c     icaser = 12 corresponds to g modes, new parameters.
      icase=0
c  icvzbn: for icvzbn .gt. 0 ensure that there is a mesh point at
c     the lower boundary of the convective envelope.
c     for icvzbn = 1 interpolate normally in aa(i,.) for i .ne. 4,
c     and reset aa(4,.) by linear interpolation near boundary.
c     for icvzbn = 2 interpolate separately below and above
c     boundary.
c     icvznb = 11 and 12 are similar to icvznb = 1 and 2, but
c     two points are inserted at base of convection zone
      icvzbn=0
c  nsmth: for nsmth .gt. 1, smooth integrand with Gaussian-weighted
c     running mean over nsmth points
      nsmth=0
      cg=1
      cx=40.
      ca=0.01
      cdgr=0
      cddgr=0
c  cdg: weight for gradient in A in interior of model
c     if cdg .gt.0, limit contribution (before scaling by cdg)
c     to partial sum of other terms. If cdg .lt. 0, do not limit
c     (and use abs(cdg) as weight).
      cdg=0
c  alphsf: cuts off the singularity at the surface of a polytropic
c     model. In setting mesh, A2 is replaced by A2/(1+alphsf*A2).
      alphsf=0
      nout=100
      cn=2
      irsu=0
      unew=3
c  iresa4: if iresa4 = 1, reset A4 before redistribution,
c  to correct possible errors in A4 resulting from near exhaustion
c  of hydrogen
      iresa4=0
      csurf=0
c  nmodel: for nmodel .gt. 1 read and redistribute nmodel models.
c          otherwise only 1.
      nmodel=0
c  kmodel: for kmodel .gt. 0, start at model no. kmodel
      kmodel=0
c  itsaml: if itsaml = 1 test input model for conversion errors
      itsaml=0
c  ioldex: if ioldex = 1, old expansion (not using second derivatives
c     in data) is used.
      ioldex=0
c
c                   ......................................
c
      nrd=0
c
c  open files needed
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Files needed: input  model on unit 2'
      if(istdpr.gt.0) 
     *  write(istdpr,*) '              output model on unit 3'
      call ofiles
      call openf(2,'o','u')
      call openf(3,'u','u')
      open(81,file='ttt.redistrb.sum',status='unknown')
c
c  read exec
    5 continue
c..      read(5,exec,end=90)
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nn,icnmsh?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nn,icnmsh
      read(5,*,end=90,err=90) nn,icnmsh
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'icase,icvzbn,nsmth?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) icase,icvzbn,nsmth
      read(5,*) icase,icvzbn,nsmth
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'cg,cx,ca,cdgr,cddgr,cdg,alphsf?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) cg,cx,ca,cdgr,cddgr,cdg,alphsf
      read(5,*) cg,cx,ca,cdgr,cddgr,cdg,alphsf
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nout,cn,irsu,unew,iresa4,csurf   ?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nout,cn,irsu,unew,iresa4,csurf
      read(5,*) nout,cn,irsu,unew,iresa4,csurf
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nmodel,kmodel,itsaml,ioldex?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nmodel,kmodel,itsaml,ioldex
      read(5,*) nmodel,kmodel,itsaml,ioldex
c
c  set parameters for interpolation and redistribution
c
      icaser=mod(icase,100)
      iredst=mod(icase/100,10)
c
c  test for setting standard parameters
c
      if(icaser.eq.1) then
c  old p modes
        cg=0.001
        cx=0.1
        ca=0.01
        cn=2
      else if (icaser.eq.2) then
c  old g modes
        cg=1
        cx=40
        ca=0.01
        cn=2
      else if(icaser.eq.3) then
c  old intermediate modes
        cg=0.05
        cx=1
        ca=0.01
        cn=2
      else if(icaser.eq.11) then
c  new p modes
        cg=0.001
        cx=0.1
        ca=0.003
        cn=2
      else if (icaser.eq.12) then
c  new g modes
        cg=4
        cx=40
        ca=0.01
        cn=2
      end if
c
c..      write(6,exec)
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nn,icnmsh'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nn,icnmsh
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'icase,icvzbn,nsmth'
      if(istdpr.gt.0) 
     *  write(istdpr,*) icase,icvzbn,nsmth
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'cg,cx,ca,cdgr,cddgr,cdg,alphsf'
      if(istdpr.gt.0) 
     *  write(istdpr,*) cg,cx,ca,cdgr,cddgr,cdg,alphsf
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nout,cn,irsu,unew,iresa4,csurf   '
      if(istdpr.gt.0) 
     *  write(istdpr,*) nout,cn,irsu,unew,iresa4,csurf
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nmodel,kmodel,itsaml,ioldex'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nmodel,kmodel,itsaml,ioldex
c
c  output summary data
c
      write(81,105) iredst,cg, cx, ca, alphsf
c
c  zero model count
c
      if(nrd.gt.kmodel) then
        nrd=0
        rewind 2
      end if
c
      nwrt=0
c
c  read model
c
 9000 read(2,end=80,err=90) nmod,nh,data,(x(n),(aa(i,n),i=1,5),n=1,nh)
c
      nrd=nrd+1
c
c  test for correct number
c
      if(kmodel.gt.0.and.nrd.lt.kmodel) go to 9000
c
c  test for no model
c
      if(nh.le.0) then
        if(istdpr.gt.0) 
     *  write(istdpr,107) nh
        go to 5
      end if
c
c  if itsaml = 1 test for conversion errors
c
      if(itsaml.eq.1) call tstaml(x,aa,nh,6)
c
c  test for resetting A4 near centre
c
      if(iresa4.eq.1) call rseta4(x,aa,nh,data,6)
c
c  for icvzbn .ge. 1 add point at lower boundary of convective envelope
c
      icvzb1=0
c
      if(icvzbn.ge.1) then
c
        call stcvzb(x,aa,data,nh,6,ncbh,0,0.d0)
          if(ncbh.gt.0) then
c
          icvzb1=mod(icvzbn,10)
	  icvzb2=icvzbn/10
c
          ncbh1=ncbh-4
          ncbh2=ncbh+4
          xcb1=x(ncbh1)
          xcb2=x(ncbh2)
c
	  if(icvzb2.eq.1) then
	    nptcvz=2
          else
	    nptcvz=1
          end if
c
c  reduce nn (to ensure that the final model has nn points)
c
          nn=nn-nptcvz
c
	end if
      end if
c
c  test for congruent mesh
c
      if(icnmsh.eq.1) go to 251
c  reset u?
      if(irsu.gt.0) aa(5,2)=unew
c
      if(nout.gt.0) then
        nd=max0(1,nh/nout)
        if(istdpr.gt.0) 
     *  write(istdpr,100)
      else
	nd=0
      end if
c
c  set integrand.
c
      fnr=1./(1+cg)
      cgr=fnr*cg
      car=fnr*ca
c
      npha=0
c
      do 20 n=1,nh
      qx=aa(1,n)
      if(x(n).eq.0) then
        php=fnr*data(5)/qx
        pha=0
        phg=cgr*(data(6)-data(5))*qx
        pdgr=0
        pddgr=0
	pdg=0
c
      else
c
        x2=x(n)*x(n)
        if(aa(2,n).ne.0) then
          aar2=aa(2,n)
          aar2=aar2/(alphsf*aar2+1)
          aar4=abs(aa(4,n))
          aar4=aar4/(alphsf*aar4+1)
        else
          aar2=1./alphsf
          aar4=1./alphsf
        end if
c
        php=fnr*aar2/(qx*x2)
        pha=car*aar2*aar2/qx
        phg=cgr*aar4*qx/x2
c
	if(npha.eq.0.and.pha.ge.4.*php) npha=n
c
        if(aa(4,n).lt.0.and.cdgr.gt.0) then
          pdgr=cdgr*aa(4,n)*aa(4,n)
        else
          pdgr=0
        end if
c
        if(n.eq.1.or.cddgr.le.0) then
          pddgr=0
        else
          n1=n-1
          pddgr=1.e-5*(aa(4,n)-aa(4,n1))/(x(n)-x(n1))
          pddgr=cddgr*pddgr*pddgr
        end if
c
c  term from gradient in A4. To avoid unreasonable behaviour
c  limit basic term (before multiplication by cdg) to be
c  at most equal to remaining terms relevant to interior
c
c  Note: as a hack, do not limit, when cdg .lt. 0.
c
        if(n.eq.1.or.cdg.eq.0) then
          pdg=0
        else
	  partsum=php+pha+abs(phg)+cx
          n1=n-1
          daa4=(aa(4,n)-aa(4,n1))/(x(n)-x(n1))
	  if(cdg.gt.0) then
            pdg=(1-x(n))*(1-x(n))*daa4
	    pdg=min(pdg*pdg,partsum)
          else
            pdg=(1-x(n))*(1-x(n))*daa4
	    pdg=pdg*pdg
	  end if
          pdg=abs(cdg)*pdg
        end if
c
      end if
c
   15 dxsi=sqrt(php+pha+abs(phg)+pdgr+pddgr+pdg+cx)
c
      if(nout.gt.0) then
        if(mod(n-1,nd).eq.0.and.istdpr.gt.0) 
     *  write(istdpr,110) n,x(n),php,pha,phg,pdgr,pddgr,pdg,cx,dxsi
      end if
c
   20 xn(n)=dxsi
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) ' npha =',npha
c
c  test for smoothing
c
      if(nsmth.gt.1) then
	call rnmean(xn,xsi,nh,1,1,nsmth)
	do 22 n=1,nh
   22   xn(n)=xsi(n)
      end if
c
c  integrate and renormalize.
c  Normalisation depends on whether surface is singular or not
c
      singsf=data(7).ge.0
c
      call vinta(x,xn,xsi,nh,1,1)
c
c  test for resetting integral near surface, in logarithmic form
c
      if(npha.gt.0.and.iredst.ge.1) then
	do 22090 n=npha,nh
        xnrat=xn(n-1)/xn(n)
	alptil(n)=((1-x(n-1))*xnrat-1+x(n))/(1-xnrat)
	ampl=xn(n)*(1-x(n)+alptil(n))
	if(istdpr.gt.0) 
     *  write(istdpr,*) n, x(n), alptil(n), ampl
22090   xsi(n)=xsi(n-1)+
     *         ampl*log((1.-x(n-1)+alptil(n))/(1-x(n)+alptil(n)))
      end if
c
      xsi(1)=0
      if(singsf) then
        fnr=(nn+cn+csurf-1)/xsi(nh)
      else
        fnr=(nn+cn-1)/xsi(nh)
      end if
      do 23 n=1,nh
      an(1,n)=xsi(n)
      xsi(n)=1+fnr*xsi(n)
   23 aa(6,n)=x(n)
c
      if(nd.eq.1) then
	if(istdpr.gt.0) 
     *  write(istdpr,112) (xsi(n),n=1,nh)
      end if
c
c..      write(6,25091) fnr,(n,x(n),xn(n),an(1,n),xsi(n),n=1,nh)
c..25091 format(//' fnr =',1pe15.7//
c..     *  ' n, x, dxsi, integrand, xsi:'/(i5,0pf12.7,1p3e15.7))
c
      go to 27
c
c  set xsi for congruent mesh
c
  251 cn=x(2)/(x(3)-x(2))+1
      xsi(1)=0
      dxsi=(nn-1)/float(nh-1)
      xsi(2)=cn+1+dxsi
      do 253 n=1,nh
      if(n.le.2) go to 253
      xsi(n)=xsi(n-1)+dxsi
  253 aa(6,n)=x(n)
c
c  interpolate to new mesh
c  ***********************
c
   27 do 28 i=1,5
   28 an(i,1)=aa(i,1)
      xn(1)=x(1)
      nst=1
      nni=nh
      ni=1
c
c  reset aa(2,.) and aa(4,.) to interpolate only in
c  quantities going as constants at x = 0.
c
      if(x(1).eq.0) then
c
        n0=2
        xs=cn+1
        aa(2,1)=data(5)
        aa(4,1)=data(6)-data(5)
c
      else
        n0=1
        xs=0
      end if
c
      do 28100 n=n0,nh
      x2=x(n)*x(n)
      aa(2,n)=aa(2,n)/x2
28100 aa(4,n)=aa(4,n)/x2
c
c  test for resetting aa(2), aa(4) and aa(5) at singular surface
c
      if(singsf) then
        nh1=nh-1
        do 28200 n=1,nh1
        x1=1.-x(n)
        aa(2,n)=aa(2,n)*x1
        aa(4,n)=aa(4,n)*x1
        if(data(7).gt.0) then
          aa(5,n)=aa(5,n)/x1**data(7)
        end if
28200   continue
c
c  set surface values from proper expansion
c
        aa(2,nh)=(data(7)+1)/aa(3,nh)
        aa(4,nh)=data(7)-aa(2,nh)
c
c  set expansion coefficient for U
c
	dt=1.d0-x(nh1)
        u01 = aa(5,nh1)
        u02=u01/(1-(3-data(7))*dt)
        corr=(8.d0/((data(7)+1)*(data(7)+2)))*u02*
     *    (dt**(data(7)+1))/(1-(3-data(7))*dt)
        u03=2*u02/(1+sqrt(1+corr))
        u0=u03
        if(istdpr.gt.0) 
     *  write(istdpr,*) ' u0 =',u0
        aa(5,nh)=u0
      end if
c
c  prepare for separate interpolation for icvzbn = 2
c
      if(icvzb1.eq.2) then
        nncbh=nh+1-ncbh
        xscbh=xsi(ncbh)
      end if
c
      n1=0
c
c  test for initializing spline interpolation
c
      if(iredst.gt.0) then
	call spline(xsi(nst),x(nst),nni,1.d31,1.d31,yspl)
      end if
c
      do 30 nstep=n0,nn
      n=nstep
      xs=xs+1
c
c  test for separate interpolation for icvzbn = 2
c
      if(icvzb1.eq.2) then
        if(xs.le.xscbh) then
          nst=1
          nni=ncbh
        else
          if(nst.eq.1) ni=1
          nst=ncbh
          nni=nncbh
        end if
      end if
c
c  interpolation, depending on iredst
c
      if(iredst.eq.0) then
c
c  standard interpolation
c
        call lir(xs,xsi(nst),an(1,n),aa(1,nst),6,6,nni,ni,inter)
        ni=2
        xn(n)=an(6,n)
c
c  test that resulting xn is in fact monotonically increasing
c
        if(n.gt.n0.and.xn(n).le.xn(n-1)) then
	  if(istdpr.gt.0) 
     *  write(istdpr,115) n, xn(n), xn(n-1)
c
c  try linear interpolation in this neighbourhood instead
c
	  ni=1
          call lir1(xs-1,xsi(nst),an(1,n-1),aa(1,nst),6,6,nni,ni,inter)
          call lir1(xs,xsi(nst),an(1,n),aa(1,nst),6,6,nni,ni,inter)
	  xn(n-1)=an(6,n-1)
          xn(n)=an(6,n)
	  if(istdpr.gt.0) 
     *  write(istdpr,117) n-1,xn(n-1),n,xn(n)
	  if(xn(n).le.xn(n-1)) stop
        end if
c
      else
c
c  interpolation using spline package
c
        call splint(xsi(nst),x(nst),yspl,nni,xs,xn(n))
c
c  test for resetting using logarithmic form
c
	if(npha.gt.0.and.xs.ge.xsi(npha)) then
c
c  locate point in original mesh
c
	  do 29010 nx=npha,nh
	  if(xs.gt.xsi(nx-1).and.xs.le.xsi(nx)) then
	    atil=(xsi(nx) - xsi(nx-1))/
     *           log((1.-x(nx-1)+alptil(nx))/(1.-x(nx)+alptil(nx)))
	    xnnew=1+alptil(nx)-
     *            (1-x(nx-1)+alptil(nx))*exp(-(xs-xsi(nx-1))/atil)
	    if(istdpr.gt.0) 
     *  write(istdpr,*) 'n, xn, xnnew', n, xn(n), xnnew
	    xn(n)=xnnew
	    go to 29020
          end if
29010     continue
29020     continue
	end if
c
	call lir(xn(n),x(nst),an(1,n),aa(1,nst),5,6,nni,ni,inter)
      end if
c
c  if icvzbn = 1 reset aa(4,.) from linear interpolation around
c  lower edge of convective envelope
c
      if(icvzb1.eq.1.and.xn(n).ge.xcb1.and.xn(n).le.xcb2) then
c
c  locate interval in original model
c
 2910   if(xsi(n1).le.xs.and.xs.lt.xsi(n1+1)) go to 2920
        n1=n1+1
        go to 2910
c
c  interpolate
c
 2920   fct2=(xs-xsi(n1))/(xsi(n1+1)-xsi(n1))
        fct1=1-fct2
        an41=fct1*aa(4,n1)+fct2*aa(4,n1+1)
c
c..        write(6,29209) n,xn(n),n1,fct1,fct2,an41,an(4,n)
c..29209   format(i5,f12.7,i4,2f12.7,1p2e13.5)
c
        an(4,n)=an41
c
      end if
c
   30 continue
c
c  reset aa(2,.) and aa(4,.)
c
      do 35 n=n0,nn
      x2=xn(n)*xn(n)
      an(2,n)=x2*an(2,n)
   35 an(4,n)=x2*an(4,n)
c
      do 40 n=n0,nh
      x2=x(n)*x(n)
      aa(2,n)=x2*aa(2,n)
   40 aa(4,n)=x2*aa(4,n)
c
c  test for resetting aa(2) and aa(4) at singular surface
c
      if(singsf) then
        nn1=nn-1
        do 40100 n=1,nn1
        xn1=1.-xn(n)
        an(2,n)=an(2,n)/xn1
        an(4,n)=an(4,n)/xn1
        if(data(7).gt.0) then
          an(5,n)=an(5,n)*xn1**data(7)
        end if
40100   continue
	xn(nn)=1.d0
        an(2,nn)=0
        an(4,nn)=0
        if(data(7).ne.0) an(5,nn)=0
c
        do 40200 n=1,nh1
        x1=1.-x(n)
        aa(2,n)=aa(2,n)/x1
        aa(4,n)=aa(4,n)/x1
40200   aa(5,n)=aa(5,n)*x1**data(7)
        aa(2,nh)=0
        aa(4,nh)=0
        if(data(7).ne.0) aa(5,nh)=0
c
      end if
c
c  reset innermost points from expansion
c
      if(x(1).eq.0) call resexp(x,aa,nh,6,xn,an,nn,6,data,ioldex)
c
c  test for setting point at convection zone boundary
c
      if(icvzb1.ge.1) call stcvzb(xn,an,data,nn,6,ncb,nptcvz,1.d-8)
c
      if(nout.le.0) go to 50
      ndn=max0(1,nn/nout)
      if(istdpr.gt.0) 
     *  write(istdpr,120) data,(n,xn(n),(an(i,n),i=1,5),n=1,nn,ndn)
   50 write(3) nmod,nn,data,(xn(n),(an(i,n),i=1,5),n=1,nn)
c
c
c  test for reading next model
c
      nwrt=nwrt+1
      if(nwrt.lt.nmodel) go to 9000
c
      go to 5
c
c  diagnostics for end of file on model file
c
   80 if(istdpr.gt.0) 
     *  write(istdpr,180)
c
   90 continue
      stop
  100 format(///' n,x,php,pha,phg,pdgr,pddgr,pdg,cx,dxsi:'/)
  105 format('iredst, cg, cx, ca, alphsf:',i3,1p4e12.4)
  107 format(//' **** error in redistrb. nh =',i10)
  110 format(i5,f10.7,1p8e13.5)
  112 format(/' xsi:'/(10f8.2))
  115 format(//' ***** error in redistrb at n =',i6/
     *         '       xn(n) = ',f12.7,'  .le. xn(n-1) =',
     *         f12.7)
  117 format(/' Try resetting with linear interpolation. n, x:'/
     *  (i5,f10.5))
  120 format(///' model on new mesh.'//' data:',1p8e13.5//
     *  ' n,x,aa(1-5):'//(i5,0pf10.7,1p5e13.5))

  180 format(//' ********** end of file on d/s 2')
      end
      subroutine resexp(x,aa,nh,iaa,xn,an,nn,ian,data,ioldex)
c
c  reset innermost points in interpolated model from expansion
c  if ioldex .ne. 1, second derivatives of p and rho in data
c  are used.
c
c  for consistency with previous version of redistribution
c  programme, also tests that an(5,n) .le. 3. this should
c  not have any effect on recent, consistent models.
c
c  original version: 9/1/1985
c
c  modified 17/5/1985, to reset an(2,.) from expansion
c  of p and of other aa-s
c
      implicit double precision (a-h, o-z)
      dimension x(nh),aa(iaa,nh),xn(nn),an(ian,nn),data(8)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      x2=x(2)*x(2)
      x4=x2*x2
c
c  test for new expansion
c
      if(ioldex.ne.1) go to 20
c
      ca1=(aa(1,2)-aa(1,1))/x2
      cc1=0
      ca4=aa(4,2)/x2
      cc4=0
      ca5=(aa(5,2)-aa(5,1))/x2
      cc5=0
c
      go to 30
c
c  new expansion
c
   20 ca1=-0.3*aa(1,1)*data(6)
      cc1=(aa(1,2)-aa(1,1)-ca1*x2)/x4
      ca4=data(6)-data(5)
      cc4=(aa(4,2)-ca4*x2)/x4
      ca5=-0.6*data(6)
      cc5=(aa(5,2)-aa(5,1)-ca5*x2)/x4
c
c  expansion for aa(3,.) (same in the two cases)
c
   30 ca3=(aa(3,2)-aa(3,1))/x2
      cc3=0
c
      if(istdpr.gt.0) 
     *  write(istdpr,100) ca1,ca3,ca4,ca5,cc1,cc3,cc4,cc5
c
      xo2=x(2)
c
      amr2=data(1)/(data(2)*data(2))
      pfct=6.6732e-8*amr2**2/(16.*atan(1.))
      cc2=0.5*aa(3,1)*data(5)
c
      do 40 n=1,nn
      xx=xn(n)
      if(xx.ge.xo2) go to 40
      x2=xx*xx
      x4=x2*x2
c
      an(1,n)=aa(1,1)+ca1*x2+cc1*x4
      an(3,n)=aa(3,1)+ca3*x2
      an(4,n)=        ca4*x2+cc4*x4
      an(5,n)=aa(5,1)+ca5*x2+cc5*x4
c
      an(2,n)=pfct*x2*an(1,n)**2*an(5,n)/(an(3,n)*data(3)*(1-cc2*x2))
c
   40 continue
c
c  restrict u to be .le. 3
c
      do 50 n=1,nn
      if(an(5,n)-3.d0.gt.1.e-10) then
        if(istdpr.gt.0) write(istdpr,110) n,xn(n),an(5,n)
        an(5,n)=3
      end if
   50 continue
c
      return
  100 format(//' coefficients in s/r resexp:'//
     *  ' ca1, ca3, ca4, ca5 =',1p4e14.6/
     *  ' cc1, cc3, cc4, cc5 =',4e14.6)
  110 format(' **** warning. at n =',i5,'  x =',f12.7,
     *  '   aa(5,n) =',1pe15.7,' .gt. 3'/
     *  '      has been reset to 3.')
      end
      subroutine stcvzb(x,aa,data,nn,iaa,ncb1,npt,eps)
c
c  sets mesh point at the botton edge of the convective envelope
c  for adiabatic oscillation model.
c  position of edge is found by extrapolation in aa(4,.) from
c  radiative interior. values of aa(i,.) at this point is
c  found from linear interpolation between neighbouring points
c  for i = 1, 2, 3 and 5, and aa(4,.) is set to zero.
c
c  this was introduced 21/5/1985. before that date extrapolation
c  from radiative interior was used for all variables.
c
c  if extrapolated edge is outside or too close to next mesh
c  point, this is taken to be bottom of convection zone,
c  and no point is added
c
c  if npt .gt. 1, insert npt points, separated by eps,
c  but assuming the same value of the variables
c
c  ncb1 returns position of bottom edge, or 0 if no edge is found
c
      implicit double precision (a-h, o-z)
      dimension x(1),aa(iaa,1),data(1)
      dimension aacvz(5)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  locate bottom of convective envelope (to allow for convective overshoot,
c  define by small positive number)
c
      ncb=0
      epscb=1.e-4
      do 10 n=2,nn
      ncb=n
      if(aa(4,n).ge.epscb.and.aa(4,n+1).lt.epscb) go to 20
   10 continue
c
c  no boundary found. write diagnostics
c
      if(istdpr.gt.0) write(istdpr,100)
      ncb1=0
      return
c
c  extrapolate linearly to zero in aa(4,.)
c
   20 ncp=ncb-1
      fcp=aa(4,ncb)/(aa(4,ncb)-aa(4,ncp))
      fcb=1.-fcp
c
c  set x at boundary, check for location
c
      ncb1=ncb+1
      xcb=fcp*x(ncp)+fcb*x(ncb)
      if(xcb.ge.0.99*x(ncb1)+0.01*x(ncb)) then
c
        if(istdpr.gt.0) write(istdpr,105) xcb,x(ncb1)
        nn1=nn
        go to 40
      end if
c
      if(npt.le.0) then
        npt1=1
      else
        npt1=npt
      end if
c
c  set interpolation coefficients
c
      fcb1=(xcb-x(ncb))/(x(ncb1)-x(ncb))
      fcb=1-fcb1
c
c  set values at boundary
c
      do 23 i=1,5
   23 aacvz(i)=fcb*aa(i,ncb)+fcb1*aa(i,ncb1)
c
c  move x and aa
c
      nns=nn-ncb
      nn1=nn+1
      do 25 n=1,nns
      n1=nn1-n
      n2=n1+npt1
      x(n2)=x(n1)
      do 25 i=1,5
   25 aa(i,n2)=aa(i,n1)
c
c  set point(s) at boundary
c
      xcb=xcb+eps
      n=ncb1+npt1
      do 35 k=1,npt1
      xcb=xcb-eps
      n=n-1
      x(n)=xcb
      do 30 i=1,5
   30 aa(i,n)=aacvz(i)
c
   35 aa(4,n)=0
c
c  diagnostic output
c
   40 n1=ncb1-10
      n2=ncb1+10
      if(istdpr.gt.0) 
     *  write(istdpr,110) (n,x(n),(aa(i,n),i=1,5),n=n1,n2)
c
      nn=nn+npt1
c
      return
  100 format(///' ********* no lower edge of convective envelope ',
     *  'found in s/r stcvzb')
  105 format(///' extrapolated x at convection zone boundary =',
     *  f10.5,'  is outside or too close to next mesh point =',
     *  f10.5/' no point added.')
  110 format(///' points near lower boundary of convective envelope.'
     *  //' n, x, aa(1-5):'//(i5,0pf12.7,1p5e15.7))
      end
      subroutine tstaml(x,aa,nn,iaa)
c
c  test for conversion error leading to zero instead
c  of 1
c
      implicit double precision (a-h, o-z)
      dimension x(1),aa(iaa,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      do 10 n=2,nn
      if(x(n).ne.0) go to 10
      if(istdpr.gt.0) write(istdpr,100) n
      x(n)=1
   10 continue
      return
  100 format(/' ********** conversion error for adiabatic',
     *  ' model. x(',i4,') = 0 has been reset to 1')
      end
      subroutine rseta4(x,aa,nn,data,iaa)
c
c  Reset A4 near centre, to correct for problems near
c  end of hydrogen burning
c  Resetting is only applied if A4/x**2 is non-monotonic
c  near centre
c
      implicit double precision (a-h, o-z)
      dimension x(1),aa(iaa,1),data(8)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  reference values
c
      axc=data(6)-data(5)
      nref=4
      axref=aa(4,nref)/x(nref)**2
c
      ireset=0
      do 10 n=2,nref-1
      ax=aa(4,n)/x(n)**2
      if((axc-ax)*(ax-axref).lt.0) ireset=1
   10 continue
c
      if(ireset.eq.1) then
	if(istdpr.gt.0) write(istdpr,110)
	axcoef=(axref-axc)/x(nref)**2
	do 20 n=2,nref-1
	x2=x(n)**2
	ax=aa(4,n)/x2
	axnew=axc+axcoef*x2
	if(istdpr.gt.0) write(istdpr,115) n, x(n), ax, axnew
   20   aa(4,n)=axnew*x2
c
      end if
c
      return
  110 format(//' Reset A4 near centre. ',
     *   ' n, x, old, new values of A4/x**2:'/)
  115 format(i5,1p3e13.5)
      end
