      subroutine nrkint(x,y,iy,nw1,nibc,nn,nev1,iasn,ii,initnr,nfit,
     *  ncnorm,det,ddsol,icry)
c
c  solves adiabatic oscillation equations using nrkm. for n .le. nev1
c  it is assumed that the factor x**(l-1) is taken out of the solution.
c
c  det returns the value of the boundary  or matching condition that
c  is not forced to be satisfied. thus to get the eigenfrequency sig
c  one must iterate on det, as a function of sig, to get det = 0.
c
c  for nibc .lt. nfit .lt. nn nrkm is used with two regions, matched at
c  n = nfit. in this case all inner and outer boundary conditions are
c  satisfied, as well as the matching conditions at nfit for
c  y(1), y(3) and y(4). det is the discontinuity in
c  y(2) at the fitting point nfit.
c
c  for nfit .le. nibc or nfit .ge. nn nrkm is called with one region.
c  here either one of the inner or one of the outer boundary conditions
c  are left unsatisfied, and are used to iterate for sig.
c  for nfit .ge. nn det is set to the surface pressure boundary condition.
c  for nfit .le. nibc  det is set to the inner displacement boundary
c  condition.
c  in both cases the boundary condition is normalized by the rms of
c  the solution at point ncnorm.
c
c  icry is returned as -1 if an error is encountered in s/r nrk.
c
c  Modified 21/2/89 to allow use with Richardson extrapolation.
c  This required a streamlining of the mesh setting etc.
c
c  Modified 23/2/89, redefining interior match quantity to be
c  normalized by central and surface values
c
c  Modified 25/2/89, removing the setting of ddsol (previously
c  set to mean change in eigenfunction)
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      include 'adipls.c.d.incl'
      dimension x(*),y(iy,*)
      dimension ap(1),aq(1),zk(1),iv(8),ea(100)
      common/worksp/ x1(nnmax),y1(4,1)
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
      common/cnrkbc/ cbc(4,4),cavc(4,4),cavs(4,4),isngbt,isngsf,ibcnrk,
     *  nnq1,nnq2,idgnrk,ibotb1
c
c  parameters controlling mesh transformation. case flagged by imcase
c  imcase = 0: no fit at nev1, no fit at nfit
c  imcase = 1:    fit at nev1, no fit at nfit
c  imcase = 2: no fit at nev1,    fit at nfit
c  imcase = 3:    fit at nev1,    fit at nfit
c
      common/cnrkrh/ imcase,nibcc,nev1c,nfitc,nns,nev11,nev12,
     *  nfit11,nfit12,nns1,iasnc
      common/bcsdat/ fctsbc, fcttbc, istsbc, ibotbc
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external bcnrk,rhsnrm,bcnrm
      data ifrsol /1/
c
      save
c
c  return ddsol as 0
c
      ddsol = 0
c
      if(idgnrk.ge.1.and.istdpr.gt.0) write(istdpr,*) 
     *  'iy,nw1,nibc,nn,nev1,ii,initnr,nfit,ncnorm',
     *   iy,nw1,nibc,nn,nev1,ii,initnr,nfit,ncnorm
c
c  test for sufficient storage for x1
c
      if(nn+2.gt.nnmax) then
        write(istdou,105) nnmax, nn
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,105) nnmax, nn
        go to 80
      end if
c
c  initialize icry to 0
c
      icry=0
c
      iy1=4
c
c  set ibcnrk in common/cnrkbc/, depending on integration case
c  note that ibcnrk = 0 is used to flag interior fitting, for
c  nibc .lt. nfit .lt. nn. ibcnrk = 1 flags condition at surface,
c  and ibcnrk = 2 flags condition at centre.
c
      if(nfit.ge.nn) then
        ibcnrk=1
      else if(nfit.le.nibc) then
        ibcnrk=2
      else
        ibcnrk=0
      end if
c
c  output of boundary conditions
c
      if(idgnrk.gt.1.and.istdpr.gt.0) then
        write(istdpr,10091) ((cbc(i,j),j=1,4),i=1,4)
        write(istdpr,10092) ((cavc(i,j),j=1,4),i=1,4)
        write(istdpr,10093) ((cavs(i,j),j=1,4),i=1,4)
        write(istdpr,10094) isngbt,isngsf
10091   format(/' cbc:'/(1p4e16.8))
10092   format(/' cavc:'/(1p4e16.8))
10093   format(/' cavs:'/(1p4e16.8))
10094   format(//' isngbt, isngsf =',2i4)
      end if
c
c  make sure that ncnorm is in range.
c
      ncnorm=max0(nibc,min0(ncnorm,nn-5))
c
c  for initnr =0 skip initializations
c
      if(initnr.eq.0) go to 25
c
c  set storage indices and flag
c
      nibcc=nibc
      iasnc=iasn
c
      if(nev1.lt.nibc+5) then
        ifnev=0
        nev1c=nibc
        nev11=0
        nev12=1
      else
        ifnev=1
        nev1c=nev1
        nev11=(nev1-nibc)/iasn+1
        nev12=nev11+1
      end if
c
      if(ibcnrk.ne.0.or.nfit.eq.nev1) then
        ifnfit=0
        nfitc=nev1c
        nfit11=nev11
        nfit12=nev12
      else 
        ifnfit=1
        nfitc=nfit
        nfit11=(nfit-nev1c)/iasn+nev12
        nfit12=nfit11+1
      end if
c
c  set case number
c
      if(ifnev.eq.0.and.ifnfit.eq.0) then
        imcase = 0
      else if(ifnev.ne.0.and.ifnfit.eq.0) then
        imcase = 1
      else if(ifnev.eq.0.and.ifnfit.ne.0) then
        imcase = 2
      else if(ifnev.ne.0.and.ifnfit.ne.0) then
        imcase = 3
      end if
c
c  set outermost point, excluding possible singularity
c
      if(isngsf.eq.0) then
        nns=nn
      else
        nns=nn-1
      end if
c
c  set normalization point on new mesh
c
      if(ncnorm.le.nev1c) then
        ncnor1=(ncnorm-nibc)/iasn+1
      else if(ncnorm.le.nfitc) then
        ncnor1=(ncnorm-nev1)/iasn+nev12
      else
        ncnor1=(ncnorm-nfit)/iasn+nfit12
      end if
c  
      nns1=(nns-nfitc)/iasn+nfit12
c
c  set total number of points in one or two regions
c
      if(ibcnrk.eq.0) then
        nnq1=nfit11
        nnq2=nns1-nfit12+1
      else
        nnq1=nns1
      end if
      if(idgnrk.ge.1.and.istdpr.gt.0) write(istdpr,107)
     *  imcase,nibcc,nev1c,nfitc,nns,nev11,nev12,
     *  nfit11,nfit12,nns1,iasnc
c
c  set x1-array
c
      do 15 n1=1,nns1
      if(n1.le.nev11) then
        n=nibc+(n1-1)*iasn
        x1(n1)=x(n)
      else if(n1.eq.nev12.and.nev11.gt.0) then
        x1(n1)=x1(nev11)+epsprc
      else if(n1.le.nfit11) then
        n=nev1c+(n1-nev12)*iasn
        x1(n1)=x(n)
      else if(n1.eq.nfit12.and.nfit11.gt.0) then
        x1(n1)=x1(nfit11)
      else
        n=nfitc+(n1-nfit12)*iasn
        x1(n1)=x(n)
      end if
c
      do 15 i=1,ii
   15 y1(i,n1)=0
c
c  prepare for call of nrk
c
   21 kk=0
      ki=0
      ucy=1
c
c  test for using alternative setting of iv
c
      if(ibotb1.ne.1.or.fcttbc.ge.0.5) then
c
        do 22 i=1,ii
   22   iv(i)=i
c
      else
	if(ii.eq.4) then
	  write(istdou,108) 
	  if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,108) 
	  icry=-1
	  return
        else
	  iv(1)=2
	  iv(2)=1
	  if(istdpr.gt.0) write(istdpr,109) fcttbc
        end if
      end if
c
      if(ibcnrk.eq.0) then
	do 23 i=1,ii
   23   iv(i+ii)=i
      end if
c
c  set normalization variable depending on possible truncation
c
      if(ibotb1.ne.1.or.fcttbc.le.0.5) then
	icnorm = 1
      else 
	icnorm = 2
      end if
c
c  find solution
c
c  test for interior fitting (2 regions in nrkm)
c
   25 if(ibcnrk.ne.0) then
c
c  no fitting. one region
c
        if(idgnrk.ge.1.and.istdpr.gt.0) write(istdpr,*) 
     *    'one region. ii, kk, ki, nns1, iy1',
     *    ii, kk, ki, nns1, iy1
        call nrkm(x1,y1,zk,ap,aq,rhsnrm,bcnrk,
     *    ii,kk,1,ki,nns1,iy1,ucy,ea,detnrk,iv)
c
      else
c
c  interior fitting. two regions
c
        if(idgnrk.ge.1.and.istdpr.gt.0) write(istdpr,*) 
     *    'two regions. ii, kk, ki, nns1, iy1',
     *    ii, kk, ki, nns1, iy1
        call nrkm(x1,y1,zk,ap,aq,rhsnrm,bcnrm,
     *    ii,kk,2,ki,nns1,iy1,ucy,ea,detnrk,iv)
c
c  combine change in first part of ea
c
        iea2=2*iy1
        iea3=iea2+iea2
        do 26 i=1,ii
        ea(i)=0.5*(ea(i)+ea(i+iy1))
        eam1=ea(i+iea2)
        eam2=ea(i+iy1+iea2)
        iam1=ea(i+iea3)
        iam2=ea(i+iy1+iea3)
        if(eam1.ge.eam2) then
          ea(i+iy1)=eam1
          ea(i+iy1+iy1)=iam1
        else
c
          ea(i+iy1)=eam2
          ea(i+iy1+iy1)=iam2
        end if
   26   continue
c
      end if
c
c  test for error in nrk or nrkm
c
      if(iv(1).eq.0) go to 80
c
c  set and print changes
c
      eam1=0
      eam12=0
      do 27 i=1,ii
      eam1=eam1+ea(i)
   27 eam2=eam2+ea(i+iy1)
c
      eam1=eam1/ii
      eam2=eam2/ii
c
      if(idgnrk.gt.0.and.istdpr.gt.0) then
c
        write(istdpr,110) eam1,(ea(i),i=1,ii)
        write(istdpr,115) eam2,(ea(i+iy1),i=1,ii)
        write(istdpr,120)      (ea(i+iy1+iy1),i=1,ii)
c
      end if
c
c  test for printing complete solution in first pass
c
      if(((idgnrk.ge.3.and.ifrsol.ne.0).or.idgnrk.ge.4)
     *  .and.istdpr.gt.0) then
        ifrsol=0
        write(istdpr,122)
        do 30 n=1,nns1
   30   write(istdpr,123) n,x1(n),(y1(i,n),i=1,ii)
      end if
c
c  test for interior fitting
c
      if(ibcnrk.ne.0) then
c
c  eigenfrequency condition at one of the boundaries
c  set normalized boundary condition into det for sig iteration
c  (coefficients are stored in cbc by s/r setbcs)
c
        det=0
        snrm=0
c
c  test which condition is to be used
c
        if(ibcnrk.eq.1) then
c
c  pressure surface condition
c
          nbc=nns1
          icbc=3
        else
c
c  inner pressure-displacement condition
c
          nbc=1
          icbc=1
c
        end if
c
c  set det
c
        do 42 i=1,ii
        snrm=snrm+y1(i,ncnor1)*y1(i,ncnor1)
   42   det=det+cbc(icbc,i)*y1(i,nbc)
c
c  renormalize det with solution at point ncnorm
c
        det=det/sqrt(snrm)
c
      else
c
c  interior fitting. eigenfrequency condition at interior point
c
        if(y1(icnorm,1).eq.0) then
          write(istdou,127) icnorm
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,127) icnorm
          icry=-1
          return
        else
          det=(y1(2,nfit11)-y1(2,nfit12))/y1(icnorm,1)
        end if
c
        if(idgnrk.gt.1) then
c
          nwr1=nfit11-5
          nwr2=nfit12+5
          n1=1
          if(istdpr.gt.0) write(istdpr,125) n1,x1(1),(y1(i,1),i=1,4),
     *                      (n,x1(n),(y1(i,n),i=1,4),n=nwr1,nwr2),
     *                      nns1,x1(nns1),(y1(i,nns1),i=1,4)
c
        end if
c
      end if
c
c  shift solution
c
      if(idgnrk.gt.1.and.istdpr.gt.0) then
        nwr1=nev11-5
        nwr2=nev12+5
        write(istdpr,44090) (n,x1(n),y1(1,n),n=nwr1,nwr2)
44090   format(//' n, x1, y1(1,.):'//(i5,0pf11.6,1pe13.5))
      end if
c
      do 45 n1=1,nns1
      if(n1.lt.nev11) then
        n=nibc+(n1-1)*iasn
      else if(n1.ge.nev12.and.n1.lt.nfit11) then
        n=nev1c+(n1-nev12)*iasn
      else if(n1.ge.nfit12) then
        n=nfitc+(n1-nfit12)*iasn
      else
        n=0
      end if
c
      if(n.gt.0) then
        do 44 i=1,ii
   44   y(i,n)=y1(i,n1)
      end if
c
   45 continue
c
      if(idgnrk.gt.3.and.istdpr.gt.0) then
c
c  extra output
c
        write(istdpr,130)
        do 50 n=nw1,nn,20
   50   write(istdpr,135) n,x(n),(y(i,n),i=1,ii)
c
      end if
c
c  possibly set solution at singular points
c
c  inner point
c
      if(isngbt.ne.0) then
        do 55 i=1,ii
        sum=0
        do 54 j=1,ii
   54   sum=sum+cavc(i,j)*y(j,nibc)
   55   y(i,nw1)=sum
c
      end if
c
c  surface
c
      if(isngsf.ne.0) then
c
        do 60 i=1,ii
        sum=0
        do 58 j=1,ii
   58   sum=sum+cavs(i,j)*y(j,nns)
   60   y(i,nn)=sum
c
      end if
c
      return
c
c  set icry to -1 for error in nrk
c
   80 icry=-1
      write(istdou,180)
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,180)
      return
c
  105 format(/' **** error in s/r nrkint. assigned storage for x1, ',
     *  i6,' is insufficient for nn =',i6)
  107 format(/' imcase,nibcc,nev1c,nfitc,nns,nev11,nev12,',
     *  'nfit11,nfit12,nns1,iasnc ='/12i6)
  108 format(//' ***** Error in s/r nrkint:',
     *  ' ii = 4 not allowed in truncated model'/
     *         '       Execution terminated.')
  109 format(//' In s/r nrkint fcttcb =',f10.5,
     *  ' Set iv = 2,1')
  110 format(/' ea(.,1). mean and values:',1pe13.5,5x,4e13.5)
  115 format( ' ea(.,2). mean and values:',1pe13.5,5x,4e13.5)
  120 format( ' ea(.,3) (loc. of max. change):',12x,4f13.1)
  122 format(//' first nrkm solution:'//' n, x1, y1(1 - ii,.):'/)
  123 format(i6,0pf11.7,1p4e14.6)
  125 format(//' n, x1, y1 at centre, fitting point, surface:'//
     *      (i4,0pf15.10,1p4e15.5))
  127 format(//' **** Error in s/r nrkint. y1(icnorm,1) = 0 for',
     *  ' icnorm =',i2)
  130 format(//' n, x, y:'//)
  135 format(i5,0pf12.7,1p4e13.5)
  180 format(//' ********** error in call of s/r nrk from s/r',
     *  ' nrkint')
      end
      subroutine bcnrk(x1,x2,y1,y2,zk1,zk2,ap,aq,g,gd,ia,n,ka,iis,
     *  kk,nn,iq)
c
c  boundary condition routine for nrk, using coefficients set into
c  common/cnrkbc/ by previous call of s/r setbcs
c
c  this condition assumes that one of the boundary conditions is
c  used to iterate for the eigenfrequency, corresponding to
c  frfit = 0 or 1.
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension y1(*),y2(*),zk1(*),zk2(*),ap(*),aq(*),g(*),gd(ia,*)
c
      common/sysord/ ii
      common/cnrkbc/ cbc(4,4),cavc(4,4),cavs(4,4),isngbt,isngsf,ibcnrk,
     *  nnq1,nnq2,idgnrk,ibotb1
c
c  test for case
c
      if(ibcnrk.eq.2) go to 30
c
c                             ----------------------------------
c
c  eigenfrequency condition at outer boundary
c  ******************************************
c
c  test for interior or outer conditions
c
      if(iq.eq.2) go to 22
c
c  set interior conditions
c
      ka=ii/2
      do 20 i=1,ka
      sum=0
      do 15 j=1,ii
      gd(i,j)=cbc(i,j)
   15 sum=sum+gd(i,j)*y1(j)
   20 g(i)=sum
c
      kk=0
      iis=ii
      nn=nnq1
      return
c
c  set outer normalization condition
c
   22 ka1=ka+1
      g(1)=y2(1)-1
      gd(1,1)=1
c
c  test for second outer condition
c
      if(ii.eq.2) return
c
      sum=0
      do 25 j=1,ii
      gd(2,j)=cbc(4,j)
   25 sum=sum+gd(2,j)*y2(j)
      g(2)=sum
c
      return
c
c                          ------------------------------------
c
c  eigenfrequency condition at inner boundary
c
c  test for interior or outer conditions
c
   30 if(iq.eq.2) go to 40
c
c  set interior conditions
c
      ka=ii/2-1
      iis=ii
      kk=0
      nn=nnq1
c
      if(ka.eq.0) return
c
      sum=0
      do 35 j=1,ii
      gd(1,j)=cbc(2,j)
   35 sum=sum+gd(1,j)*y1(j)
      g(1)=sum
c
      return
c
c  set outer normalization condition
c
   40 g(1)=y2(1)-1
      gd(1,1)=1
c
c  set remaining outer conditions
c
      kb=ii/2
c
      i=1
      do 50 i1=1,kb
      i=i+1
      icbc=2+i1
      sum=0
      do 45 j=1,ii
      gd(i,j)=cbc(icbc,j)
   45 sum=sum+gd(i,j)*y2(j)
   50 g(i)=sum
c
      return
c
      end
      subroutine rhsnrm(x,y,zk,ap,aq,f,fd,h,hd,ia,n1,iq)
c
c  right hand side routine for nrm, calling previous right hand
c  side routine for adiabatic oscillations
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      parameter(iaa=10)
      dimension y(*),zk(*),ap(*),aq(*),f(*),fd(ia,*),h(*),hd(*),
     *  finh(100)
      common/rhsdat/ el,ell,alb,els,el1,sigcom,anres,perfac,data(8),
     *  aa(iaa,1)
      common/sysord/ ii
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
      common/cnrkbc/ cbc(4,4),cavc(4,4),cavs(4,4),isngbt,isngsf,ibcnrk,
     *  nnq1,nnq2,idgnrk,ibotb1
      common/cnrkrh/ imcase,nibcc,nev1c,nfitc,nns,nev11,nev12,
     *  nfit11,nfit12,nns1,iasnc
      common/xarra/ xin(1)
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
c
c  common controlling output
c
      common/coutpt/ nout,nprcen,iper,irotkr,nprtkr,igm1kr,npgmkr,
     *  nfmode,nfmesh,ispcpr,npout,nobs_stmx,nfmscn
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data nxdiag /0/
c
c  prepare for call of rhs routine for s/r linint
c
c  set point n in original model corresponding to point n1 in 
c  selected mesh
c
      if(idgnrk.gt.0.and.istdpr.gt.0) then
        write(istdpr,*) 'Enter rhsnrm at n, x =',n1, x
      end if
      if(n1.le.nev11) then
        n=nibcc+(n1-1)*iasnc
      else if(n1.le.nfit11) then
        n=nev1c+(n1-nev12)*iasnc
      else
        n=nfitc+(n1-nfit12)*iasnc
      end if
c
c  set el1 depending on whether or not in evanescent region
c
      if(n1.le.nev11) then
        if(iplneq.ne.1) then
          el1=el-1
        else
          el1=el
        end if
c
      else
c
        el1=0
      end if
c
c  test that xin(n)=x
c
      if(abs(xin(n)-x).gt.2*epsprc.and.nxdiag.le.100) then
        nxdiag=nxdiag+1
        write(istdou,120) iq,n1,n,x,xin(n)
        if(nxdiag.le.2) write(istdou,125)
     *    imcase,nibcc,nev1c,nfitc,nns,nev11,nev12,
     *    nfit11,nfit12,nns1,iasnc
        if(istdpr.ne.istdou.and.istdpr.gt.0) then
          write(istdpr,120) iq,n1,n,x,xin(n)
          if(nxdiag.le.2.and.istdpr.gt.0) write(istdpr,125)
     *      imcase,nibcc,nev1c,nfitc,nns,nev11,nev12,
     *      nfit11,nfit12,nns1,iasnc
        end if
      end if
c
   22 call adirhs(x,fd,finh,ia,n)
c..c
c..c  when outputting eigenfunctions in scan, reset Poisson's equation
c..c  to vacuum form
c..c
c..      if(ii.eq.4.and.nfmscn.eq.1.and.iq.eq.2) then
c..        fd(4,1)=0.d0
c..        fd(4,2)=0.d0
c..        fd(4,3)=ell
c..        fd(4,4)=2.d0
c..      end if
c
c  set f from fd and y
c
      do 30 i=1,ii
      sum=0
      do 25 j=1,ii
   25 sum=sum+fd(i,j)*y(j)
   30 f(i)=sum
c
      if(idgnrk.gt.1.and.istdpr.gt.0) write(istdpr,*) 
     *  'return from rhs. f =',(f(i),i=1,ii)
      if(idgnrk.gt.2.and.n.ge.nn-3.and.istdpr.gt.0) then
        write(istdpr,130) n
	do 40 i=1,ii
   40   write(istdpr,135) (fd(i,j),j=1,ii)
      end if
      return
  120 format(' *** error in s/r rhsnrm.',
     *  ' iq, n1, n, x, xin(n):',i2,2i6,2f15.10)
  125 format(/' imcase,nibcc,nev1c,nfitc,nns,nev11,nev12,',
     *  'nfit11,nfit12,nns1,iasnc ='/12i6)
  130 format(' fd at n =',i5)
  135 format(1p4e16.8)
      end
      subroutine bcnrm(x1,x2,y1,y2,zk1,zk2,ap,aq,g,gd,ia,n,ka,iis,
     *  kk,nn,iq)
c
c  boundary condition routine for nrkm, using coefficients set into
c  common/cnrkbc/ by previous call of s/r setbcs
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      parameter(iaa=10)
      dimension y1(*),y2(*),zk1(*),zk2(*),ap(*),aq(*),g(*),gd(ia,*)
c
      common/sysord/ ii
      common/rhsdat/ el,ell,alb,els,el1,sigcom,anres,perfac,data(8),
     *  aa(iaa,1)
      common/cnrkbc/ cbc(4,4),cavc(4,4),cavs(4,4),isngbt,isngsf,ibcnrk,
     *  nnq1,nnq2,idgnrk,ibotb1
c
c  common controlling output
c
      common/coutpt/ nout,nprcen,iper,irotkr,nprtkr,igm1kr,npgmkr,
     *  nfmode,nfmesh,ispcpr,npout,nobs_stmx,nfmscn
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(idgnrk.gt.0.and.istdpr.gt.0) then 
	write(istdpr,*) 'Enter bcnrm with x1, x2, iq =', x1, x2, iq
	write(istdpr,*) 'nnq1, nnq2 =',nnq1,nnq2
      end if
c
c  test for region
c
      go to (10,50,20), iq
c
c  set interior conditions
c
   10 ka=ii/2
      do 17 i=1,ka
      sum=0
      do 15 j=1,ii
      gd(i,j)=cbc(i,j)
   15 sum=sum+gd(i,j)*y1(j)
   17 g(i)=sum
c
      iis=ii
      kk=0
      nn=nnq1
c
      return
c
c                                  -------------------------
c
c
c  set outer conditions
c
c  set outer normalization condition
c
   20 g(1)=y2(1)-1
      gd(1,1)=1
c
c  remaining outer conditions
c
      kb=ii/2
      do 30 i=1,kb
      i1=i+1
      i2=i+2
c
      sum=0
      do 25 j=1,ii
      gd(i1,j)=cbc(i2,j)
   25 sum=sum+gd(i1,j)*y2(j)
   30 g(i1)=sum
c
      return
c
c  set matching conditions. When computing eigenfunctions for an 
c  unconverged solution (e.g. during a scan) replace matching on
c  y_4 with boundary condition on potential
c
   50 iis=ii
      kk=0
      ka=ii-1
      nn=nnq2
c
c  test for boundary condition on y_4 (switch off, 24/7/09)
c
c  simple matching
c
c..      if(nfmscn.ne.1.or.ii.eq.2) then
      do i=1,ka
        i1=i
        if(i.gt.1) i1=i+1
        g(i)=y1(i1)-y2(i1)
        gd(i,i1)=1
        gd(i,i1+ii)=-1
      end do
c
c..      else
c..c
c..c  apply boundary condition
c..c
c..        do i=1,ka-1
c..          i1=i
c..          if(i.gt.1) i1=i+1
c..          g(i)=y1(i1)-y2(i1)
c..          gd(i,i1)=1
c..          gd(i,i1+ii)=-1
c..        end do
c..c
c..c  As a check, try zero-density condition here.
c..c
c..c..        uu=aa(5,n)
c..        uu=0.d0
c..c
c..        g(ka)=-uu*y1(1)+(el+uu)*y1(3)+y1(4)
c..        gd(ka,1)=-uu
c..        gd(ka,3)=el+uu
c..        gd(ka,4)=1.d0
c..      end if
      return
c
      end
