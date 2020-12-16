      subroutine sigitr(sig,x,y,iy,nw1,nibc,nn,
     *  mdintg,nev1,iasn,ii,ig,ifrich,nfit,istsb1,iscan,
     *  xfit1,det,ddsig,ddsol,ds1,iord,icry,isolcv,isigcv,istatus)
c
c  Iterates to determine sig = sigma**2 by integrating
c  adiabatic oscillation equations for given l = el.
c  (transmitted in common/rhsdat/). Integration method is determined
c  by mdintg.
c
c  Returns matching condition in det, absolute value of last 
c  relative change in sig in ddsig, possibly test of eigenfunction
c  continuity in ddsol and matching determinant, to be used in
c  s/r mchsol, in ds1 (ddsol and ds1 are not set for mdintg = 3).
c
c  Original version 30/6/95
c
c  Modified 22/10/96, to set flag for no convergence if det eq detp.
c
      implicit double precision (a-h, o-z)
      logical fulmod, radial
      parameter(iaa=10)
      dimension x(1), y(iy,1), ds1(4,1)
c
      common/rhsdat/ el,ell,alb,els,el1,sigcom,anres,perfac,data(8),
     *  aa(iaa,1)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  set flag for radial oscillations
c
      radial=el.le.1.e-6
      fulmod=nibc.eq.2
c
c  flags for convergence
c
      isolcv=1
      isigcv=1
c
c  initialize det and detp to zero
c
      det=0
      detp=0
c
c  take out factor x**l in solution up to point nev1 where solution
c  begins to vary more slowly than x**(l/2), or to nfit if
c  nfit .le. nev1.
c
c  For consistency with old version of programme, do not reset
c  boundary for iteration during scan, if irsevn = 0
c
      if(ifrich.ne.2.and.(iscan.le.1.or.irsevn.ge.1)) then
        els=ell/sig
        call stevft(x,nn,nibc,nfit,iasn,radial,fulmod,nev1,xfit1)
      end if
c
c  set flag for initialization in nrkint (only used for mdintg=3)
c
      initnr=1
c
      if(istdpr.gt.0) write(istdpr,137)
c
      isig=-1
c
c  begin integration with new sig
c  ******************************
c
   14 isig=isig+1
      detpp=detp
      detp=det
      els=ell/sig
      sigcom=sig
c
c  carry out integration for given sig and l
c
      call sigint(x,y,iy,nw1,nibc,nn,mdintg,nev1,iasn,ii,ig,
     *  initnr,nfit,istsb1,isig,det,ddsol,ds1,iord,icry)
c
c  test for error in integration
c
      if(icry.lt.0) return
c
c  output of matching results
c  **************************
c
      if(mdintg.ne.2) then
c
        if(istdpr.gt.0) write(istdpr,120) isig,sig,det,ddsol
      else
c
        if(istdpr.gt.0) write(istdpr,121) isig,sig,det,ddsol,iord
c
      end if
c
c  for itmax = 0, test continuity of eigenfunction
c
      if(itmax.eq.0) then
        if(ddsol.ge.epssol) then
c
c  solution discontinuous. set flags for no convergence
c
          write(istdou,159)
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,159)
          isolcv=-1
          isigcv=-1
        end if
c 
        go to 26
c
      end if
c
c  first step
c
      if(isig.eq.0) then
c
c  set second sigma**2
c
        sigp=sig
        dsig=fsig*sigp
        sig=sigp+dsig
        go to 14
c
      end if
c
c  subsequent steps, set change in sig
c
      if(abs(detp-det).le.epsprc*abs(detp)) then
c
c  no significant change in det. set dsig to 0 to force end of iteration
c  and set flags for no convergence.
c
        write(istdou,119) det,detp,sig,sigp
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,119) 
     *    det,detp,sig,sigp
        dsig=0
        isolcv=-1
        isigcv=-1
      else
c
        dsig=dsig*det/(detp-det)
c
      end if
c
      ddsig=abs(dsig/sig)
c
c  limit ddsig to less than dsigmx
c
      if(ddsig.ge.dsigmx) then
        dsig=sign(dsigmx,dsig)*abs(sig)
        ddsig=dsigmx
      end if
c
      sig=sig+dsig
c
c  test for convergence
c  ********************
c
      if(ddsig.gt.eps) then
c
c  frequency not converged. try again if maximum number of iterations
c  not exceeded
c
        if(isig.lt.itmax) then
          go to 14
        else
c
c  unconverged. do not file mode and/or write summary
c
          isigcv=-1
          isolcv=-1
	  istatus=1
          write(istdou,122) ddsig,eps
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,122) 
     *      ddsig,eps
          go to 26
c
        end if
c
      else
c
c  frequency converged. test for convergence of eigenfunction
c
        if(ddsol.lt.epssol) then
c
c  everything converged. finish for this mode
c
          go to 26
c
        else
c
c  eigenfunction not converged.
c
c  test for convergence of sig to full accuracy. otherwise try again.
c
          if(ddsig.gt.epsprc.and.isig.lt.itmax) then
            go to 14
c
c  give up on eigenfunction
c
          else
c
            write(istdou,124) ddsig,ddsol,epssol
            if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,124) 
     *        ddsig,ddsol,epssol
            isolcv=-1
	    istatus=2
c
c  end of tests for convergence
c
          end if
c
        end if
c
      end if
c
c  end of routine
c
   26 continue
c
      if(istdpr.ne.istdou.and.isigcv.eq.1) write(istdou,185) isig, sig
c
      return
  137 format(/)
  119 format(//' ********** no change in det in s/r adipls.'/
     *  12x,'det, detp =',1p2e23.15/
     *  12x,'sig, sigp =',2e23.15)
  120 format(' isig,sig,det,det/detnrm:',i4,f17.10,1p2e11.3)
  121 format(' isig,sig,det,det/detnrm,order:',i3,f16.9,1p2e11.3,i8)
  122 format(//1x,10(1h*),'  iteration unconverged. dsig =',
     *  1pe13.5,' .gt.  eps =',e13.5)
  159 format(//' **** solution discontinuous for',
     *  ' itmax = 0')
  124 format(//1x,10(1h*),' ddsig =',1pe13.5,' has converged,',
     *  ' but ddsol =',e13.5,' .gt. epssol =',e13.5/
     *  11x,' eigenfunction may be inaccurate')
  185 format(' Iteration converged after ',i3,
     *  ' steps.    sigma**2 =',1pe12.4)
      end
