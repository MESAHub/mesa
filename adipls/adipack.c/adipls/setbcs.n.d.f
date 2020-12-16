      subroutine setbcs(x,y,iy,nn,nibc,
     *  icncbc,istsb1,ii,ig,nd,ny,isig,istnrk)
c
c  when istnrk .ne. 1 sets values of solution at centre and surface
c
c  when istnrk = 1 sets up coefficients (in common/cnrkbc/) for
c  b.c. routine for nrk, and the coefficients relating the solution
c  at the singular points to those at the neighbouring point
c
c  Diagnostics are flagged by setting kdgbcs (in common/cdgsbr/):
c
c  kdgbcs = -1: ibotbc reset to 1 for radial modes
c  kdgbcs = -2: truncation point is in oscillatory region
c  kdgbcs = -10: beta+ solution not implemented for istsbs = 1
c  kdgbcs = -20: isothermal surface boundary condition applied 
c  	      above acoustical cut-off frequency.
c  Contributions can be added (may need a little more care)
c
c  when icncbc = 1 sets scaled dependent variables x**(1-l)*y(i,.)
c  in non-radial case for a full model. these variables tend to
c  non-zero constants at the centre.
c  when icncbc = 2 sets the original dependent variables y(i,.)
c
c  this modification was added 28/1/1985
c
c  test for singular surface modified 1/9/1985.
c  old test was that aa(2,nn) = 0 which was somewhat artificial.
c  new test is that data(7) .ge. 0.
c
c  modified 25/3/1986 to allow transition from delta p = 0 to p prime
c  = 0. this is controlled by parameter fctsbc going from 0 to 1.
c  note that previously fctsbc was used in an approximation to
c  mixing between beta + and beta - solutions.
c
c  12/8/87: version from RECKU and version from Boulder (1985) 
c  merged.
c
c  modified 13/8/87 to standardize output
c
c  modified 19/2/88 to include cowling approximation, and variable
c  lambda, for the radial case.
c
c  modified 30/11/88, to include option ibotbc = 2, corresponding
c  to bottom condition d(xir)/dx = 0, matching Stein & Nordlund.
c
c  modified 13/2/90 to fix surface boundary condition for isopycnic
c  model in radial case
c
c  modified 16/12/93, to make proper expansion, including terms
c  O(t**npol) and O(t**(npol+1)) in expansion at singular surface
c
c  Modified 20/1/95 for resetting of treatment of turbulent pressure
c
c  Modified 28/6/95 moving data, aa into common /rhsdat/ and
c  fctsbc, fctcbc, ibotbc into common /bcsdat/
c
c  Modified 2/7/95, adding (original) istsbc to common /bcsdat/,
c  replacing istsbc by istsb1 in bulk of routine.
c
c  Modified 14/1/04, adding option istsbc = 9, to set zero displacement
c  at surface (probably needs further checking!)
c
c  ....................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      logical radial,fulmod,singsf
      parameter(iaa=10)
      dimension x(nn),y(iy,nn)
      common/rhsdat/ el,ell,alb,els,el1,sig,anres,perfac,data(8),
     *  aa(iaa,1)
      common/bcsdat/ fctsbc, fcttbc, istsbc, ibotbc
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
      common/cnrkbc/ cbc(4,4),cavc(4,4),cavs(4,4),isngbt,isngsf,ibcnrk,
     *  nnq1,nnq2,idgnrk,ibotb1
      common/csexpc/ y11, y13, y21, y23, y31, y33, y41, y43,
     *  yt11, yt13, yt21, yt23, yt31, yt33, yt41, yt43, 
     *  yt241, yt243, ytt41, ytt43 
c
c  common containing diagnostic flags possibly set in various parts of the
c  code
c
      common/cdgsbr/ kdgbcs
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  as a start, initialize diagnostic flag
c
      kdgbcs=0
c
c  coefficients for expansion
      data ac11,ac31,ac12,ac32 /1.d0,0.d0,1.d0,1.d0/,
     *     as11,as31,as12,as32 /1.d0,0.d0,0.d0,1.d0/
c
      radial=el.le.1.e-6
      fulmod=nibc.eq.2
      singsf=data(7).ge.0
c
c  set upper limit for loop over solutions, depending on ig
c
      if(ig.eq.1) then
	kex=1
      else if(ig.eq.2) then
	kex=5
      else
	write(istdou,101) ig
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,101) ig
	stop
      end if
c
c  if aa(3,n) is not gamma 1, assume this value at centre and
c  surface
c
      gm1=1.6666666666667d0
c
c  test for initialization of nrk coefficients
c
      if(istnrk.eq.1) then
c
        isngbt=0
        isngsf=0
        do 5 i=1,4
        do 5 j=1,4
        cbc(i,j)=0
        cavc(i,j)=0
    5   cavs(i,j)=0
      end if
c
c                        ------------------------------------
c
c  conditions at innermost point
c  *****************************
c
c  test for truncated model
c
      if(fulmod) then
	ibotb1 = -1
	go to 20
      end if
c
      ibotb1 = ibotbc
c
c  conditions at base of truncated model
c
      if(radial.and.ibotbc.eq.0) then
        ibotbc=1
	ibotb1 = 1
        write(istdou,105)
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,105)
	kdgbcs=kdgbcs-1
      end if
c
      if(ibotbc.eq.1) then
c
c  set condition with fcttbc
c
        c21=fcttbc
        c11=1-fcttbc
c
      else if(ibotbc.eq.0) then
c
c  condition at truncation point to isolate solution growing towards
c  surface
c
        qxs=aa(1,nibc)/sig
        ak1=1-aa(2,nibc)/(qxs*ell)
        ak2=1-qxs*aa(4,nibc)
        akr=ak1*ak2
        if(akr.gt.0) then
c
c  set coefficients
c
          c11=1
          c21=sqrt(ell*akr)/ak1
c
        else
c
c  point in oscillatory region. use fcttbc
c
          if(isig.eq.0.and.istdpr.gt.0) write(istdpr,107) 
     *      x(nibc),sig,fcttbc
          c21=fcttbc
          c11=1-fcttbc
	  ibotb1 = 1
	  kdgbcs=kdgbcs-2
        end if
c
      else if(ibotbc.eq.2) then
c
c  set d(xir)/dx = 0, to match Stein & Nordlund
c
	if(radial) then
	  c11=-sig*aa(2,nibc)/(x(nibc)*aa(1,nibc))
	else
	  c11=1-sig*aa(2,nibc)/(ell*aa(1,nibc))
	end if
	c21=2-aa(2,nibc)
c
      else if(ibotbc.eq.3) then
c
c  set d(xir/x)/dx = 0, as suggested by JOP
c
	if(radial) then
	  c11=-sig*aa(2,nibc)/(x(nibc)*aa(1,nibc))
	else
	  c11=1-sig*aa(2,nibc)/(ell*aa(1,nibc))
	end if
	c21=3-aa(2,nibc)
c
      end if
c
c  test for setting nrk coefficients
c
      if(istnrk.eq.1) then
        cbc(1,1)=c21
        cbc(1,2)=-c11
      else
c
c  set solution at bottom
c
        y(1,nibc)=c11
        y(2,nibc)=c21
      end if
c
      go to 40
c
c                         --------------------------------
c
c  expansion at centre
c  *******************
c
   20 p2=data(5)
      rh2=data(6)
      a2=rh2-p2
      bc=aa(1,1)
      x1=x(2)
c
      if(radial) go to 30
c
c  nonradial case
c  **************
c
      ec=els*bc
      x2=0.5d0*x1*x1/(2.d0*el+3.d0)
      f2=1.d0-sig/(bc*el)
      c1=(p2*(el+2)-ec*a2)*x2
      c2=ell*(p2-(el+3.d0)*bc*a2/sig)*x2
c
c  set expansion coeffients
c
      c11=1.d0+c1*f2
      c13=-c1
      c21=el+1.d0+c2*f2
      c23=-c2
      c31=-3.d0*alb*(a2+(el+1)*p2/ec)*x2
      c3=3.d0*(0.4d0*(el-1.d0)*rh2+a2+(1.d0-alb)*p2)*x2
      c33=1.d0+c3
      c41=el*c31
      c43=el-2.d0+el*c3
c
c  test for setting nrk coefficients
c
      if(istnrk.ne.1) go to 24
c
c  set flag for singular bottom
c
      isngbt=1
c
c  test for ii = 2
c
      if(ii.eq.2) go to 22
c
c  set boundary condition coefficients in full case
c
      delta=c11*c33-c31*c13
      cbc(1,2)=-delta
      cbc(1,1)=c21*c33-c23*c31
      cbc(1,3)=c23*c11-c21*c13
      cbc(2,4)=-delta
      cbc(2,1)=c41*c33-c43*c31
      cbc(2,3)=c43*c11-c41*c13
c
c  set coefficients for central solution
c
      cavc(1,1)=c33/delta
      cavc(1,3)=-c13/delta
      cavc(2,1)=(el+1.d0)*c33/delta
      cavc(2,3)=-(el+1)*c13/delta
      cavc(3,1)=-c31/delta
      cavc(3,3)=c11/delta
      cavc(4,1)=-(el-2.d0)*c31/delta
      cavc(4,3)=(el-2.d0)*c11/delta
c
      go to 40
c
c  restricted case
c
   22 cbc(1,1)=c21
      cbc(1,2)=-c11
      cavc(1,1)=1.d0/c11
      cavc(2,1)=(el+1.d0)/c11
      go to 40
c
c  set solution at centre
c  **********************
c
c  first set of coefficients. to allow for cowling approximation
c  a30 is set explicitly to zero.
c
   24 a10=ac11
      a30=0.d0
c
      do 28 k=1,kex,4
      y(k,2)=c11*a10+c13*a30
      y(k+1,2)=c21*a10+c23*a30
c
      if(ii.eq.2) go to 25
      y(k+2,2)=c31*a10+c33*a30
      y(k+3,2)=c41*a10+c43*a30
c
c  central values
c
   25 y(k,1)=a10
      y(k+1,1)=(el+1.d0)*a10
c
      if(ii.eq.2) go to 26
      y(k+2,1)=a30
      y(k+3,1)=(el-2.d0)*a30
   26 a10=ac12
   28 a30=ac32
c
c  test for rescaling, to include x**(el-1)
c
      if(icncbc.ne.2.or.el.eq.1) go to 40
c
      xl1=x1**(el-1.d0)
      do 29 k=1,kex,4
      i1=k-1
      do 29 i=1,ii
      i1=i1+1
      y(i1,2)=xl1*y(i1,2)
      if(el.gt.1) y(i1,1)=0.d0
   29 continue
      go to 40
c
c  radial case
c  ***********
c
   30 x2=x1*x1
      if(aa(3,1).gt.0.1) gm1=aa(3,1)
      f2=p2*sig/bc
c
c  test for cowling approximation
c
      if(icow.eq.0) then
        c1=p2*(sig/bc+3*alb)
      else
        c1=f2
      end if
c
      c12=x1*f2*(-1.d0+0.1d0*(p2*(1.d0-3.d0*gm1)+c1)*x2)/3.d0
      c22=1.d0+x2*(0.5d0*a2-c1/6.d0)
c
c  test for setting nrk coefficients
c
      if(istnrk.eq.1) then
c
c  set nrk coefficients
c
        cbc(1,1)=c22
        cbc(1,2)=-c12
        cavc(2,2)=1.d0/c22
        isngbt=1
c
      else
c
c  set solution
c
c  central values
c
        y(1,1)=0
        y(2,1)=1.d0
c
        y(1,2)=c12
        y(2,2)=c22
c
      end if
c
c                        ---------------------------------
c
c  conditions at surface
c  *********************
c
   40 an=data(7)
      if(singsf)  go to 60
c
c  nonsingular surface
c  *******************
c
      nd=nn
      ny=nd+1
c
c  set parameter for pressure condition (delta p or p prime)
c
      fctsb1=1-fctsbc
c
c  test for inclusion of isothermal atmosphere
c
      if(istsb1.ne.1) go to 50
c
c  set coefficients for isothermal atmosphere condition
c
      vg=aa(2,nn)
      if(aa(3,1).gt.0.1) gm1=aa(3,nn)
c
      sigx=sig/aa(1,nn)
      elsx=ell/sigx
c  value of a for isothermal atmosphere
      a=vg*(gm1-1.d0)
      gamma=a-vg+4.d0
      gamma=gamma*gamma+4.d0*(elsx-vg)*(sigx-a)
c  test that gamma is positive
      if(gamma.ge.0) go to 45
c  write diagnostic, use standard condition
      write(istdou,115) sig,gamma
      if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,115) sig,gamma
      istsb1=0
      kdgbcs=kdgbcs-20
      go to 50
c
c  set coefficient to y1
c
   45 c21=((vg-a+sqrt(gamma))/2.d0-2.d0)/(vg-elsx)
c..      write(6,*) 'sig,sigx,vg,a,gamma,c21:'
c..      write(6,*) sig,sigx,vg,a,gamma,c21
c
c  set coefficients and possibly solution
c  **************************************
c
c
   50 continue
c
c  test for zero displacement
c
      if(istsb1.eq.9) then
	c11=0
	c13=0
	c21=1
	c23=0
	c31=0
	c33=1
	c41=0
	c43=-el-aa(5,nn)
	go to 70
      end if
c
      if(radial) go to 56
c
c  nonradial case
c  **************
c
      c11=1.d0
      c13=0.d0
      c31=0.d0
      c33=1.d0
      c41=aa(5,nn)
      c43=-el-c41
      if(istsb1.eq.1) go to 54
c
c  standard case
c
      etat=els*aa(1,nn)
      c21=fctsb1*etat
      c23=-aa(10,nn)*etat
      go to 70
c
c  isothermal atmosphere
c
   54 c21=elsx*c21
      c23=-elsx*(1.d0-(el+1.d0-elsx)/(vg+a))
      go to 70
c
c  radial case
c  ***********
c
   56 c11=1.d0
      c13=0.d0
      c23=0.d0
      if(istsb1.eq.1) go to 57
      c21=fctsb1*aa(1,nn)*x(nn)/sig
      go to 70
   57 c21=x(nn)*c21/sigx
      go to 70
c
c  singular surface
c  ****************
c
   60 nd=nn-1
      ny=nd+1
      ny1=nd+2
      ang=(an+1.d0)/gm1
      as=an-ang
      el2=el-2.d0
      dt=1.d0-x(nd)
c
c  set expansion coefficient for U
c
      u01 = aa(5,nd)/dt**an
      u02=u01/(1.d0-(3.d0-an)*dt)
      corr=(8.d0/((an+1.d0)*(an+2.d0)))*u02*(dt**(an+1.d0))/
     *  (1.d0-(3.d0-an)*dt)
      u03=2.d0*u02/(1.d0+sqrt(1.d0+corr))
      u0=u03
c..      write(6,*) ' u0 =',u0
c
      if(radial) go to 65
c
c  nonradial case
c  **************
c
c  test for an = 0
c
      if(an.gt.1.e-10) then
        dtn=dt/(an+1.d0)
c
c..        c11=1-dt-((els-3)*(1+as)+ang*(1+sig))*dtn
        y11=-1.d0-((els-3.d0)*(1.d0+as)+ang*(1.d0+sig))/(an+1.d0)
	yt11=(1.d0-alb)*u0/(gm1*(2.d0*an+1.d0))
	c11=1.d0+y11*dt+yt11*dt**(an+1.d0)
c..        c13=((1+as)*els+(el+1)*ang)*dtn
        y13=((1.d0+as)*els+(el+1.d0)*ang)/(an+1.d0)
	yt13=0.d0
	c13=y13*dt
c..        c21=els*(1-dt)+(els*as*(6-els)+(1+ang)*(2*els-ell))*dtn
        y21=-els+(els*as*(6.d0-els)+(1.d0+ang)*(2.d0*els-ell))/(an+1.d0)
	yt21=(1.d0-alb)*els*u0*as/((an+1.d0)*(2.d0*an+1.d0))
	c21=els+y21*dt+yt21*dt**(an+1.d0)
c..        c23=els*(-1+(as*(els-3-el)-2*(1+ang))*dtn+dt)
        y23=els*(1.d0+(as*(els-3.d0-el)-2.d0*(1.d0+ang))/(an+1.d0))
	ct23=0.d0
	c23=-els+y23*dt
c..        c31=0
	y31=0.d0
	yt31=-alb*u0/(an+1.d0)
        c31=yt31*dt**(an+1.d0)
c..        c33=1+(el-1)*dt
        y33=el-1.d0
	yt33=u0/(an+1.d0)
        c33=1.d0+y33*dt+yt33*dt**(an+1.d0)
c..        c41=0
        y41=0.d0
	yt41=u0*alb
c
	yt241=u0*(alb*as*y11+(alb*ang/els)*y21-(as+(1.d0-alb)*ang)*y31
     *      +alb*(an*(an-2.d0)-2.d0-3.d0*ang))/(1.d0+an)
c..	ytt41 = u0*u0*(alb/(an+1))*((1-alb)*(2/(2*an+1) + 1) +
c..     *          3*an/((an+1)*(an+2)))
	ytt41 = u0*u0*alb*(((1.d0-alb)/(2.d0*an+1.d0)**2.d0)*(4*an+1.d0-
     *          2.d0*(an+1.d0)/gm1)/gm1+
     *          (an+4.d0)/((an+1.d0)*(an+2.d0)))
	if(anres.eq.0) then
          c41=yt41*dt**an+yt241*dt**(an+1.d0)+ytt41*dt**(2.d0*an+1.d0)
        else
          c41=yt241*dt**(an+1.d0)+ytt41*dt**(2.d0*an+1.d0)
        end if
c..        c43=-el*c33
	y43=-el*(el-1.d0)
	yt43=-u0
	yt243=u0*(alb*as*y13+(alb*ang/els)*y23-(as+(1.d0-alb)*ang)*y33
     *      +(3.d0*alb*ang-an*(an-2.d0)+4.d0-2.d0*el))/(1.d0+an)
	ytt43=-u0*u0*(an+4.d0)/((an+1.d0)*(an+2.d0))
	if(anres.eq.0) then
	  c43=-el+y43*dt+yt43*dt**an+yt243*dt**(an+1.d0)
     *        +ytt43*dt**(2*an+1.d0)
        else
	  c43=-el+y43*dt+yt243*dt**(an+1.d0)
     *        +ytt43*dt**(2*an+1.d0)
	end if
        us=0
c
      else
c
c  iso-pycnic surface (corresponding to polytrope of index 0)
c
c..        c11=(els*(ang-1)-ang*(4+sig)+2)*dt
        y11=els*(ang-1.d0)-ang*(4.d0+sig)+2.d0
c..        c13=((1+el)*ang-els*(ang-1))*dt
        y13=(1.d0+el)*ang-els*(ang-1.d0)
c..        c21=(els*(els*ang+1-4*ang)-ell*(1+ang))*dt
        y21=els*(els*ang+1.d0-4.d0*ang)-ell*(1.d0+ang)
c..        c23=els*((1+el-els)*ang-1)*dt
        y23=els*((1.d0+el-els)*ang-1.d0)
c..        c31=-3*dt
        y31=-3.d0
c..        c33=(2+el)*dt
        y33=2.d0+el
c..        c41=3+3*ang*(-c11+c21/els+c31)+12*dt
        y41=3.d0*ang*(-y11+y21/els+y31)+12.d0
c..        c43=-(3+el)-(6+el*(el+5))*dt+3*ang*(-c13+c23/els+c33)
        y43=-(6.d0+el*(el+5.d0))+3.d0*ang*(-y13+y23/els+y33)
c..        c11=c11+1
        c11=1.d0 + y11*dt
	c13 = y13*dt
c..        c21=els+c21
        c21=els+y21*dt
c..        c23=-els+c23
        c23=-els+y23*dt
c..        c33=1+y33*dt
	c31 = y31*dt
        c33=1+y33*dt
	c41=3.d0+y41*dt
	c43=-(3.d0+el)+y43*dt
        us=3.d0
      end if
c
      go to 70
c
c  radial case
c  ***********
c
   65 continue
c..      c11=1+dt*(2*gm1-4-sig)/gm1
      y11=(2.d0*gm1-4.d0-sig)/gm1
      yt11 = 0.d0
      c11=1.d0+y11*dt
c..      c13=0
      y13=0.d0
      yt13=0.d0
      c13=0.d0
c
c  test for iso-pycnic surface
c
      if(an.le.1.e-10) then
c..        c21=(1-dt*(3+4/gm1+sig*(1+1/gm1)))/sig
        y21=-(3.d0+4.d0/gm1+sig*(1.d0+1.d0/gm1))/sig
	c21=1.d0/sig+y21*dt
	yt21=0.d0
      else
c..        c21=(1+dt*(4*as-sig*(ang+1))/(an+1))/sig
        y21=(4.d0*as-sig*(ang+1.d0))/((an+1.d0)*sig)
	yt21 = u0/(sig*(an+1.d0))
	c21=1.d0/sig+y21*dt+yt21*dt**(an+1.d0)
      end if
c
      y23=0.d0
      yt23=0.d0
      c23=0.d0
c
c  test for setting nrk coefficients
c
   70 if(istnrk.ne.1) go to 80
c
c  set nrk coefficients
c
c  test for setting flag for singular surface
c
      if(an.ge.0) isngsf=1
c
c  test for ii = 2
c
      if(ii.eq.2) go to 72
c
c  set boundary condition coefficients in full case
c
      delta=c11*c33-c31*c13
      cbc(3,2)=-delta
      cbc(3,1)=c21*c33-c23*c31
      cbc(3,3)=c23*c11-c21*c13
      cbc(4,4)=-delta
      cbc(4,1)=c41*c33-c43*c31
      cbc(4,3)=c43*c11-c41*c13
c
      if(isngsf.eq.0) return
c
c  set coefficients for surface solution
c
      cavs(1,1)=c33/delta
      cavs(1,3)=-c13/delta
      cavs(2,1)=els*(c31+c33)/delta
      cavs(2,3)=-els*(c11+c13)/delta
      cavs(3,1)=-c31/delta
      cavs(3,3)=c11/delta
      cavs(4,1)=(us*c33+(el+us)*c31)/delta
      cavs(4,3)=-(us*c13+(el+us)*c11)/delta
c
      return
c
c  restricted case
c
   72 cbc(3,1)=c21
      cbc(3,2)=-c11
c..      write(6,*) 'cbc(3,1),cbc(3,2):',cbc(3,1),cbc(3,2)
c
      if(isngsf.eq.0) return
c
      cavs(1,1)=1.d0/c11
c
c  test for radial case
c
      if(radial) go to 75
c
      cavs(2,1)=els/c11
      return
c
   75 cavs(2,1)=1.d0/(sig*c11)
      return
c
c                           ----------------------------------
c
c  set solution
c  ************
c
   80 a10=as11
c
c  to allow for cowling approximation force a30 to be zero in first
c  solution
c
      a30=0.d0
c
      do 90 k=1,kex,4
      y(k,ny)=a10*c11+a30*c13
      y(k+1,ny)=a10*c21+a30*c23
c
      if(ii.eq.2) go to 83
c
      y(k+2,ny)=a10*c31+a30*c33
      y(k+3,ny)=a10*c41+a30*c43
c
c  test for setting surface value in singular case
c
   83 if(an.lt.0) go to 88
      y(k,ny1)=a10
c
c  test for radial
c
      if(radial) go to 84
c
      y(k+1,ny1)=els*(a10-a30)
      go to 85
c
   84 y(k+1,ny1)=a10/sig
c
   85 if(ii.eq.2) go to 88
c
      y(k+2,ny1)=a30
      y(k+3,ny1)=us*a10-(el+us)*a30
c
c  new coefficients
c
   88 a30=as32
   90 a10=as12
c
      return
  101 format(//' ***** ig =',i4,' not allowed in s/r setbcs.'/
     *         '       Execution terminated.')
  105 format(//1x,10(1h*),' ibotbc = 0 not allowed for radial modes.'/
     *  11x,' ibotbc has been reset to 1')
  107 format(///1x,10(1h*),'  xtrnct =',f16.6,' is in oscillatory',
     *  ' interval for sig =',1pe13.5,'.  fcttbc =',e13.5,' used'//)
  110 format(//1x,10(1h*),' beta+ solution not implemented for',
     *  ' istsbc = 1. ')
  115 format(//1x,10(1h*),' for sig =',1pe13.5,'  gamma =',e13.5,
     *  ' is negative'/12x,'standard condition used, istsbc reset to 0')
      end
