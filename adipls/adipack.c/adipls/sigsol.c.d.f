      subroutine sigsol(sig,x,y,iy,nw1,nibc,nn,nnw,mdintg,nev1,nfit,
     *  istsb1,iscan,iord,icry,isolcv,isigcv,sigtst)
c
c  Finds eigenfrequency and eigenfunction of adiabatic oscillation,
c  starting from trial frequency sigma**2 = sig, and with degree given 
c  by el in common/rhsdat.
c  Possibly makes initial integration with Cowling approximation
c  and/or uses Richardson extrapolation.
c
c  ddsig, xfit1 and ddsol are returned from s/r sigitc and are 
c  therefore set in common/csumma/ for inclusion in grand summary.
c
c  Original version 1/7/95
c
c  Modified 30/3/00, dropping test for same order with Richarson
c  extrapolation if itsord = -1
c
c  Modified 30/4/02, generalizing changes of nfit
c
c  Modified 17/11/05, adding sigtst to argument list. This returns the
c  eigenvalue set on the thinned mesh, with Richardson extrapolation,
c  for test of proper interval with scan. Without Richardson extrapolation
c  the eigenvalue is returned.
c
c  Modified 8/3/06, adding reset of fsig to original value for iteration
c  on full mesh with Richardson extrapolation.
c  Also include option of dsigre .le. -1, to reset trial frequency
c  on full mesh with Richardson extrapolation when different orders
c  are obtained. The number of attempts is controlled by i_dsigremax
c  hardcoded below.
c
      implicit double precision (a-h, o-z)
      logical radial
      parameter(iaa=10,iyri=4)
      dimension x(*), y(iy,*)
      dimension csummm(50)
c
      common/sysord/ ii
      common/yyyyri/ yri(iyri,1)
      common/rhsdat/ el,ell,alb,els,el1,sigcom,anres,perfac,data(8),
     *  aa(iaa,1)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
      common/cincnt_new/ fsig0
      common/csumma/ xmod1,datsum(8),datmd1(2),xtrnct,dum13,xfit1,
     *  fsbcsm,fcbcsm,albsum,
     *  elsum,ordsum,sigsum,sigc,dmx,xmx,ekin,per,perv,frqv,
     *  ddsig,ddsol,ysum(4),dosdum(5),
     *  in,nnwcom,mdints,ivarfs,icase,iorign,idum7,idum8,mlname(4)
      common/cmdtst/ iordtr,icaswn, sigwn1, sigwn2, frqwn1, frqwn2, 
     *  iorwn1, iorwn2, frlwn1, frlwn2
c
c  common defining file input and output
c
      common/cdadsg/ imlds,idslog,idsgsm,idsssm,idsefn,
     *   idsrkr,idsgkr
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      equivalence(xmod1,csummm(1))
c
      save
c
c  store original dsigre
c
      dsigre0=dsigre
c
c  initialize istatus
c
      istatus=0
c
c  initialize parameters for resetting nfit
c
      nftrcs=nftmax/10
      nftmx1=mod(nftmax,10)
c
      radial = el.le.1.e-6
c
c  initialize counter for trying different fitting points
c  during Richardson extrapolation
c
      nftcnt=0
c
c  store nfit and trial frequency for possible retries 
c  with modified nfit
c
      nfitst=nfit
      sigst=sig
      xfit1=x(nfitst)
      xfitst=xfit1
c..      write(6,*) ' nfitst set to', nfitst
c
c  temporarily store icow in icowc, in case icow is modified
c  during solution (for icow = 1)
c
      icowc=icow
c
   10 dsgrcp = 0
c
c  Test for Richardson extrapolation, using thinned mesh initially
c
      if(iriche.eq.1) then
        ifrich = 1
        iasn = 2
      else
        ifrich = 0
        iasn = 1
      end if
c
c  test for radial case
c
      if(radial) then
        ii=2
        ig=1
        if(icow.eq.1) then
          write(istdou,110) 
          if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,110) 
          icow=0
        end if
c
c  nonradial case, test for full system
c  
      else if(icow.eq.0) then
        ii=4
        ig=2
c
c test for Cowling approximation for initial solution
c
      else if(icow.eq.1) then
        ii=2
        ig=1
c
        call sigitc(sig,x,y,iy,nw1,nibc,nn,
     *    mdintg,nev1,iasn,ii,ig,ifrich,nfit,istsb1,iscan,
     *    xfit1,ddsig,ddsol,iord,icry,isolcv,isigcv,istatus)
c
c  find correction to sig, and y3, when using cowling approximation
c
        call gravpo(x(nw1),y(1,nw1),data,aa,sig,iy,iaa,nnw,iasn,dsig,
     *    y(3,nw1),iy,icry)
        sig=sig+dsig
	if(istdpr.gt.0) then
          write(istdpr,120) dsig,sig
          write(istdpr,125)
	end if
        ii=4
        ig=2
        icow=0
c
      else
c
c  simple Cowling approximation
c
        ii=2
        ig=1
c
      end if
c
c  entry point for start of proper sigma iteration
c  entry point for retrying iteration with modified nfit
c
   20 call sigitc(sig,x,y,iy,nw1,nibc,nn,
     *  mdintg,nev1,iasn,ii,ig,ifrich,nfit,istsb1,iscan,
     *  xfit1,ddsig,ddsol,iord,icry,isolcv,isigcv,istatus)
c
c  now eigenfrequency and eigenfunction have been determined
c
c  set possibly corrected sigma**2, including effect of gravitational
c  potential, and y(3,.)
c
      sigc=sig
      dsig=0
c
      if(ig.eq.1.and..not.radial) then
c
c  find correction to sig, and y3, when using cowling approximation
c
        call gravpo(x(nw1),y(1,nw1),data,aa,sig,iy,iaa,nnw,iasn,dsig,
     *    y(3,nw1),iy,icry)
        sigc=sig+dsig
        if(istdpr.gt.0) write(istdpr,120) dsig,sigc
c
        do 25 n=nw1,nn
   25   y(4,n)=0
c
      end if
c
c  find order
c
c  for modes of high degree, exclude most of evanescent region in
c  determination of order.
c
      if(el.lt.50) then
        nwor1=nw1
      else
        nwor1=nev1
      end if
      nnwor=nn-nwor1+1
c
      call order(x(nwor1),y(1,nwor1),data,aa(1,nwor1),el,sig,icow,
     *  irsord,iy,iaa,iord,nnwor,iasn,mdintg,0)
c
c  test for continuing Richardson extrapolation case
c  note with Cowling approximation, uncorrected frequency
c  is used when icow = 3
c
      if(iriche.eq.1) then
c
c  test for convergence. otherwise skip to end
c
        if(isigcv.lt.0) go to 90
c
        if(icow.eq.3) then
          sigrc1 = sig
	  dsig1 = 0
        else
          sigrc1 = sigc
	  dsig1 = dsig
        end if
c
c  set sig for testing of proper interval in s/r sigscn
c
        sigtst = sig
c
        iord1=iord
c
c  store  for later comparison
c
        do 30 i=1,4
        do 30 n=nw1,nn
   30   yri(i,n)=y(i,n)
c
        ifrich = 2
        iasn = 1
c
c  reset fsig to original value
c 
	fsig = fsig0
c
c  reset trial frequency
c
        if(dsigre.gt.-1) sig=sig*(1.d0+dsigre)
c
c  Find solution on full mesh
c
        call sigitc(sig,x,y,iy,nw1,nibc,nn,
     *    mdintg,nev1,iasn,ii,ig,ifrich,nfit,istsb1,iscan,
     *    xfit1,ddsig,ddsol,iord,icry,isolcv,isigcv,istatus)
c
c  store order as possibly set by eigint
c
        iord2=iord
c
c  set possibly corrected sigma**2, including effect of gravitational
c  potential, and y(3,.)
c
        sigc=sig
	dsig=0
c
        if(ig.eq.1.and..not.radial) then
c
c  find correction to sig, and y3, when using cowling approximation
c
          call gravpo(x(nw1),y(1,nw1),data,aa,sig,iy,iaa,nnw,iasn,dsig,
     *      y(3,nw1),iy,icry)
          sigc=sig+dsig
          if(istdpr.gt.0) write(istdpr,120) dsig,sigc
c
          do 35 n=nw1,nn
   35     y(4,n)=0
c
        end if
c
c  test for convergence. otherwise skip to end
c
        if(isigcv.lt.0) go to 90
c
        if(icow.eq.3) then
          sigrc2 = sig
	  dsig2 = 0
        else
          sigrc2 = sigc
	  dsig2 = dsig
        end if
c
c  test that order is the same
c
        call order(x(nwor1),y(1,nwor1),data,aa(1,nwor1),el,sig,icow,
     *    irsord,iy,iaa,iord2,nnwor,iasn,mdintg,0)
c
        if(iord1.ne.iord2.and.itsord.ne.-1) then
c
c  test on frequency and eigenfunction change
c
          iasnp = 2
          idgtri = 0
          call testri(x(nw1),y(1,nw1),yri(1,nw1),iy,iyri,nnw,iasnp,
     *      aa,iaa,el,iord1,iord2,sigrc1,sigrc2,dsgrcp,dsolrc,
     *      isigcv,istatus,idgtri)
        end if
c
c  For Cowling approximation, test that corrections to sigma**2 are
c  similar
c
	if(dsig2.ne.0.and.isigcv.eq.1) then
	  if(abs(dsig1/dsig2-1).ge.0.2) then
	    write(istdou,127) dsig1, dsig2
	    if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,127) 
     *        dsig1, dsig2
	    isigcv=-3
	    istatus=9
          end if
        end if
c
c  test for resetting fitting point to attempt to get a solution
c  with order and/or eigenfunctions that agree
c
        if(isigcv.eq.1) then
c
c  Carry through Richardson extrapolation
c
          sigrce= (4*sigrc2 - sigrc1)/3
          perrce=perfac/sqrt(sigrce)
          frqrce=16.666666666667d0/perrce
          csummm(37)=frqrce
          write(istdou,130) sigrc1, sigrc2, sigrce, frqrce
          if(istdpr.ne.istdou.and.istdpr.gt.0) 
     *      write(istdpr,130) sigrc1, sigrc2, sigrce, frqrce
c
          dsgrcp=abs(sigrc1-sigrc2)
c
c  reset order to be for full mesh
c
          iord = iord2
c
        else if(i_dsigre.le.i_dsigremax.and.dsigre0.le.-1) then
	  i_dsigre = i_dsigre+1
	  if(iord2.gt.iord1) then
	    dsigre = -i_dsigre*0.5*abs(sigrc2/sigrc1-1)
          else
	    dsigre =  i_dsigre*0.5*abs(sigrc2/sigrc1-1)
          end if
	  write(istdou,'(/'' Reset dsigre to'',1pe13.5/)') dsigre
	  if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *      write(istdpr,'(/'' Reset dsigre to'',1pe13.5/)') dsigre
          sig=sigst
          istatus=6
          ifrich=1
          iasn=2
          go to 20
c
c  try to reset trial frequency for full mesh (by setting dsigre)
c
        else if(nftcnt.le.nftmx1) then
c
c  retry iteration with new nfit
c
          call rsnfit(x,nibc,nn,nev1,nfit,xfit1,xfitst,nftcnt,nftrcs)
c
          sig=sigst
          istatus=6
          ifrich=1
          iasn=2
          go to 20
        else
          write(istdou,145) nftmx1
          if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,145) nftmx1
	  if(istatus.ne.9) istatus=4
	  go to 90
        end if
c
c  End of Richardson extrapolation block
c
      else
c
c  set sig for testing of proper interval in s/r sigscn
c
        sigtst = sig
c
      end if
c
c  test for correct order, if itsord = 1
c
      if(itsord.eq.1.and.iord.ne.iordtr) then
        write(istdou,150) iord, iordtr
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,150) 
     *    iord, iordtr
        if(nftcnt.le.nftmx1) then
c
c  retry solution with new nfit
c
          call rsnfit(x,nibc,nn,nev1,nfit,xfit1,xfitst,nftcnt,nftrcs)
c
          sig=sigst
          istatus=8
          go to 10
        else
          write(istdou,155) nftmx1
          if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,155) nftmx1
	  istatus=7
        end if
      end if
c
      if(radial) then
c
c  in radial case set y(3,n) and y(4,n) to zero
c
        do 50 i=3,4
        do 50 n=nw1,nn
   50   y(i,n)=0
      end if
c
   90 continue
c
c  reset nfit, lest it has been modified to attempt getting same 
c  mode in Richardson extrapolation or correct order.
c
      if(nfit.ne.nfitst) then
        write(istdou,180) nfitst, nfit
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,180) 
     *    nfitst, nfit
        nfit=nfitst
      end if
c
c  reset icow, lest icow = 1 has been used
c
      if(icow.ne.icowc) then
        write(istdou,185) icowc, icow
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,185) 
     *    icowc, icow
        icow=icowc
      end if
c
c  test for adding to status log
c
      if(istatus.gt.0) call setlog(idslog,istatus,el,iord,iordtr,sig,
     *  eps,epssol,ddsig,ddsol,iord1,iord2,sigrc1,sigrc2,dsolrc,
     *  dsig1,dsig2,xfit1,xfitst)
c
      dsigre=dsigre0
      i_dsigre=0
      i_dsigremax=3
c
      return
  110 format(/' Warning. icow = 1 for radial mode.',
     *  ' icow temporarily reset to 0')
  120 format(//' correction to sig =',1pe13.5,'.  new sig =',
     *  e13.5)
  125 format(//' repeat iteration with full equations')
  127 format(
     *  ' ***** Error in s/r sigsol. With Richardson extrapolation,'/
     *  '       gravitational potential corrections =',1p2e13.5/
     *  '       differ excessively')
  130 format(//' results of Richardson extrapolation'/
     *  ' sigrc1 =',1pe13.5,' sigrc2 =',e13.5,
     *  ' extrapolated sig =',e13.5/' cyclic frequency =',e13.5,' mHz')
  145 format(//' ***** Richardson extrapolation failed after ',i3,
     *         ' attempts to change nfit'/
     *         '       Give up.')
  150 format(//' ***** Warning in s/r sigsol. Computed order =',i5,
     *  ' .ne. input order =',i5)
  155 format(//' ***** Getting correct order failed after ',i3,
     *         ' attempts to change nfit'/
     *         '       Give up.')
  180 format(/' Reset nfit to ',i5,' from temporary value ',i5)
  185 format(/' Reset icow to ',i2,' from temporary value ',i2)
      end
