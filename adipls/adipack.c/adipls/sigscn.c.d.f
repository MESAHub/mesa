      subroutine sigscn(sig1,sig2,iscan,nsig,itrsig,x,y,iy,
     *  nw1,nibc,nn,nnw,mdintg,nev1,nfit,istsb1,icry)
c
c  Scans in sig = sigma**2 between limits sig1 and sig2, with iscan steps,
c  for the value of l = el given in common /rhsdat/.
c  The type of scan is determined by nsig.
c  If itrsig = 1 in addition iterates for modes at zero crossings
c  of matching determinant and outputs results.
c
c  Scan is controlled by nsig:
c   nsig = 1: uniform in squared frequency (i.e., in sig)
c   nsig = 2: uniform in frequency (i.e., in sqrt(sig))
c   nsig = 3: uniform in 1/(frequency) (i.e., in 1/sqrt(sig))
c   nsig = 4: uniform in 1/(frequency) at low frequency and 
c             in frequency at high frequency.
c
c  Original version 2/7/95
c
c  Modified 18/10/00, setting sig before call of stevft.
c
c  Modified 4/5/02, resetting fsig before iteration
c
c  Modified 14/5/02, introducing possibility of fine scan near local
c  minimum or single-point sign change in fitting determinant.
c
c  Modified 7/3/06, to include option nsig = 4.
c
c  Modified 24/7/09, allowing output of unconverged eigenfunctions
c  during simple scan, for setting core phase, flagged by setting
c  nfmode.
c
      implicit double precision (a-h, o-z)
      logical radial, fulmod
      parameter(iaa=10)
      dimension x(1), y(iy,1)
      dimension ds1(4,5)
c
      common/sysord/ ii
      common/rhsdat/ el,ell,alb,els,el1,sig,anres,perfac,data(8),
     *  aa(iaa,1)
      common/bcsdat/ fctsbc, fcttbc, istsbc, ibotbc
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl,inomd1
c
c  common controlling output
c
      common/coutpt/ nout,nprcen,iper,irotkr,nprtkr,igm1kr,npgmkr,
     *  nfmode,nfmesh,ispcpr,npout,nobs_stmx,nfmscn
c
c  common containing diagnostic flags possibly set in various parts of the 
c  code
c
      common/cdgsbr/ kdgbcs
c
c  common controlling diagnostics
c
      common/cdiagn/ idgrhs, idgrh1, iprdet, itssol, idgtss
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
c  store original iscan and nsig
c
      iscan_orig=iscan
      nsig_orig=nsig
c
c  set ifine = -1 to disable using finer scan
c
      ifine=0
c
c  initialize integration parameters
c
      istsb1=istsbc
      radial = el.le.1.e-6
      fulmod = nibc.eq.2
c
c  Test for Richardson extrapolation. If so, use thinned mesh
c
      if(iriche.eq.1) then
	iasn = 2
      else
	iasn = 1
      end if
c
c  test for Cowling approximation (note that for icow = 1 scan is
c  done in Cowling approximation)
c
      if(radial.or.icow.gt.0) then
	ii = 2
	ig = 1
      else
	ii = 4
	ig = 2
      end if
c
c  set initial boundary of evanescent region
c
      els=ell/sig1
      sig=sig1
      call stevft(x,nn,nibc,nfit,iasn,radial,fulmod,nev1,xfit1)
c
      isig=0
      isigfn=1
c
c  Step for sigma*2 stepping.
c  case determined by nsig (as in function signew)
c
      iscan1=iscan-1
      if(nsig.eq.1) then
c
c  linear in squared frequency
c
        dsigs=(sig2-sig1)/iscan1
      else if(nsig.eq.2) then
c
c  linear in frequency
c
        dsigs=(sqrt(sig2)-sqrt(sig1))/iscan1
      else if(nsig.eq.3) then
c
c  linear in period
c
        dsigs=(1.d0/sqrt(sig2)-1.d0/sqrt(sig1))/iscan1
      else if(nsig.ne.4) then
c
	write(istdou,110) nsig
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,110) nsig
        stop
c
      end if
c
      det=0
      detp=0
      sigp=0
      istsbp=istsbc
      istspp=istsbc
c
c  initialize flag for stopping scan
c
      istop=0
c
c  set flag for initialization in nrkint (only used for mdintg=3)
c
      initnr=1
c
      sig=sig1
      if(istdpr.gt.0) write(istdpr,115)
c
   10 continue
c
c  test for stopping scan (e.g., due to diagnostics from boundary condition)
c  restoring original iscan and nsig
c
      if(istop.gt.0) then
	write(istdou,'(//'' Premature stop of scan''/)')
	if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *    write(istdpr,'(//'' Premature stop of scan''/)')
	iscan=iscan_orig
        nsig=nsig_orig
	return
      end if
c
c  entry for continuing scan
c
      isig=isig+1
      if(isig.ge.2) then
	sigpp=sigp
	sigp=sig
	if(nsig.ne.4) then
          sig=signew(1,nsig,1,0,sigp,dsigs,itrds,el,icry)
        else
          sig=signw_pg(nsig,sig1,sig2,iscan,isig,icry)
        end if
      end if
c
c  add diagnostics
c
c..      if(isig.ge.1300) iprdet=1
c
c  stop before end
c
c..      if(isig.ge.1400) stop
c
      els=ell/sig
c
      detpp=detp
      detp=det
      istspp=istsbp
      istsbp=istsb1
c
c  integration
c
      call sigint(x,y,iy,nw1,nibc,nn,mdintg,nev1,iasn,ii,ig,
     *  initnr,nfit,istsb1,isig,det,ddsol,ds1,iord,icry)
c
c  test for error in integration
c
      if(icry.lt.0) return
c
c  test for setting flag to stop scan
c
      if(inomd1.ge.2.and.kdgbcs.lt.0) istop=1
c
c  output of matching results
c
      if(mdintg.ne.2) then
        if(istdpr.gt.0) write(istdpr,120) isig,sig,det,ddsol
      else
        if(istdpr.gt.0) write(istdpr,121) isig,sig,det,ddsol,iord
      end if
c
c  test for iteration
c
      if(itrsig.eq.1.and.isig.ge.2) then
c
c  test for possible finer search (dummy write statement required
c  12/1/07 to ensure correct operation??!!)
c
        if(isig.gt.isigfn+1.and.detpp*det.gt.0.and.ifine.eq.0
     *    .and.2*detp/(detpp+det).le.0.6) then
c
          isigst=isig
	  sigfst=sig
	  sigpst=sigp
	  detst=det
	  detpst=detp
	  iscans=iscan
	  nsigst=nsig
	  dsigst=dsigs
	  nsig=1
	  iscan=50
	  dsigs=(sig-sigpp)/(iscan-1)
	  isig=1
	  isigfn=1
	  ifine=1
	  write(istdou,130)
	  if(istdou.ne.istdpr.and.istdpr.gt.0)
     *      write(istdpr,132) sigpp, sig
	  istsb1=istspp
	  sig=sigpp
	  go to 10
c
c  test for renewed fine scan after completing one set
c
        else if(isig.eq.isigfn+1.and.det*detp.lt.0.and.ifine.eq.0)
     *    then
c
          isigst=isig
	  sigfst=sig
	  sigpst=sigp
	  detst=det
	  detpst=detp
	  iscans=iscan
	  nsigst=nsig
	  dsigst=dsigs
	  nsig=1
	  iscan=25
	  dsigs=(sig-sigp)/(iscan-1)
	  isig=1
	  isigfn=1
	  ifine=1
	  write(istdou,130)
	  if(istdou.ne.istdpr.and.istdpr.gt.0)
     *      write(istdpr,132) sigp, sig
	  istsb1=istsbp
	  sig=sigp
	  det=detp
	  go to 10
c
c  set flag for possible mode
c  test for change of sign in det
c
        else if(isig.gt.isigfn+1.and.(detpp*detp.lt.0.or.detp.eq.0)) 
     *    then
	  ismode = 1
	  sigi1=sigpp
	  sigi2=sigp
	  sigi3=sig
	  deti1=detpp
	  deti2=detp
	  deti3=det
c
c  if latest scan exceeded the acoustical cut-off frequency, 
c  attempt to reset istsb1 to original value
c
	  if(istsb1.ne.istsbp) then
	    if(istdpr.gt.0) write(istdpr,134) istsbp
	    istsb1 = istsbp
          end if
        else if(isig.eq.iscan.and.(detp*det.lt.0.or.det.eq.0)) then
	  ismode = 1
	  sigi1=sig
	  sigi2=sigp
	  sigi3=sigpp
	  deti1=det
	  deti2=detp
	  deti3=detpp
	else
	  ismode = 0
        end if
c
c  for mdintg = 3 test for possible singularity instead of zero,
c  resetting ismode to 0
c
        if(mdintg.eq.3.and.isig.ge.3.and.
     *    abs(deti2).gt.abs(deti3)) ismode = 0
c
	if(ismode.eq.1) then
c
c  begin iteration
c
          if(istdpr.gt.0) write(istdpr,145)
          sigst = sig
	  fsigst=fsig
	  iretry=0
c
c  entry point after possible resetting of interval for iteration
c
   20     fsig=0.1d0*abs(sigi2-sigi1)/sigi1
	  if(deti2.ne.deti1) then
	    dsig=(sigi1-sigi2)*deti1/(deti2-deti1)
	    sig=sigi1+dsig
          end if
          sig0=sig
	  els=ell/sig
c
          call sigsol(sig,x,y,iy,nw1,nibc,nn,nnw,mdintg,nev1,nfit,
     *      istsb1,iscan,iord,icry,isolcv,isigcv,sigtst)
c
c  test that the solution is in the correct interval, and converged.
c  otherwise, try again.
c
	  if((sigtst.lt.min(sigi1,sigi2).or.sigtst.gt.max(sigi1,sigi2)
     *      .or.isolcv.lt.0).and.iretry.le.3) then
c
	    if(isolcv.ge.0) then
	      if(istdpr.gt.0) write(istdpr,146) sig,sigi1,sigi2
            else
	      if(istdpr.gt.0) write(istdpr,147)
            end if
c
	    iretry=iretry+1
	    sig=sig0
            call sigint(x,y,iy,nw1,nibc,nn,mdintg,nev1,iasn,ii,ig,
     *        initnr,nfit,istsb1,isig,det0,ddsol,ds1,iord,icry)
c
c  test for interval with change of sign
c
	    if(det0*deti1.lt.0) then
	      deti2=det0
	      sigi2=sig0
            else
	      deti1=det0
	      sigi1=sig0
            end if
c
	    go to 20
c
	  end if
c
	  fsig=fsigst
          if(icry.lt.0) return
c
c  call output routine
c
          call sigout(sig,x,y,iy,nw1,nibc,nn,nnw,mdintg,nev1,nfit,
     *      istsb1,iord,isolcv,isigcv)
c
c  continue scan after iteration
c
          sig = sigst
c
          istsb1=istsbc
          initnr=1
c
          if(istdpr.gt.0) write(istdpr,150)
c
        end if
c
c  end of section for iteration for solution
c  Test for output of eigenfunction during scan with no iteration
c
      else if(nfmscn.gt.0) then
c
c  call output routine to output eigenfunction
c  (might need some cleaning up later)
c
        isigcv=1
        isolcv=1
        call sigout(sig,x,y,iy,nw1,nibc,nn,nnw,mdintg,nev1,nfit,
     *    istsb1,iord,isolcv,isigcv)
        initnr=1
      end if
c
c  continue scan
c
      if(isig.lt.iscan) then
	go to 10
c
c  test for continuing after fine scan
c
      else if(ifine.eq.1) then
	ifine=0
	iscan=iscans
	isig=isigst
	isigfn=isig
	sig=sigfst
	sigp=sigpst
	det=detst
	detp=detpst
	dsigs=dsigst
	nsig=nsigst
	if(istdpr.gt.0) write(istdpr,155)
	go to 10
c
      end if
c
      return
  110 format(//' ***** Error in s/r sigscn. nsig = ',i4,' not allowed')
  115 format(///'  scan in sigma**2'//)
  120 format(' isig,sig,det,det/detnrm:',i6,f17.10,1p2e11.3)
  121 format(' isig,sig,det,det/detnrm,order:',i3,f16.9,1p2e11.3,i8)
  130 format(//' Start fine scan'/)
  132 format(//' Fine scan from sig =',1pe13.6,' to',e14.6)
  134 format(//' Try to reset istsb1 to',i2,
     *  ' to use isothermal b.c. for iteration')
  145 format(//' begin iteration for sig'/)
  146 format(//' ***** Warning. Converged sig =',1pe13.5/
     *         '       not in correct interval: (',2e13.5,')'//
     *         '       Try new iteration.'/)
  147 format(//' ***** Warning. Iteration not converged.'/
     *         '       Try new iteration.'/)
  150 format(//1x,30(1h*)//' continue scan in sig'/)
  155 format(//' End fine scan'/)
      end
      double precision function signw_pg(nsig,sig1,sig2,iscan,isig,icry)
c
c  Sets squared frequency corresponding to step isig in scan with
c  iscan steps between squared frequencies sig1 and sig2.
c  Stepping is approximately uniform in 1/frequency for squared frequencies
c  well below sigc, hardcoded below, and approximately uniform in frequency
c  for squared frequencies well above sigc.
c
c  Original version: 7/3/06
c
      implicit double precision (a-h, o-z)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  transitional squared dimensionless frequency
c
      sigc=4.d0
c
c  test for initialization
c
      if(isig.le.2) then
	iscan1=iscan-1
	beta=(sqrt(sig2)-sqrt(sig1)+
     *    sigc*(1.d0/sqrt(sig1)-1.d0/sqrt(sig2)))/iscan1
	gamma=(sigc/sqrt(sig1)-sqrt(sig1))/beta
	bsigc=4.d0*sigc/beta**2
      end if
c
      disig=isig-1-gamma
      sign=0.5d0*beta*(disig+sqrt(disig*disig+bsigc))
      signw_pg=sign*sign
      icry=0
      return
      end
