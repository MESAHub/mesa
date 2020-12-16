      subroutine sigitc(sig,x,y,iy,nw1,nibc,nn,
     *  mdintg,nev1,iasn,ii,ig,ifrich,nfit,istsb1,iscan,
     *  xfit1,ddsig,ddsol,iord,icry,isolcv,isigcv,istatus)
c
c  Iterates to determine sig = sigma**2 by integrating
c  adiabatic oscillation equations.
c  This is a driver routine for s/r sigitr where the
c  actual iteration takes place. 
c  In this routine, possibly modify location nfit of
c  fitting point if solution has not converged.
c  Note: on exit from the routine nfit and xfit1 may therefore have 
c  been changed.
c
c  Also carry out matching of solution when mdintg .ne. 3. and
c  scaling of solution in evanescent region.
c
c  Returns absolute value of last 
c  relative change in sig in ddsig, possibly test of eigenfunction
c  continuity in ddsol (ddsol is not set for mdintg = 3).
c
c  Original version 1/7/95
c
c  Modified 30/4/02, generalizing changes of nfit
c
      implicit double precision (a-h, o-z)
      logical radial
      parameter(iaa=10)
      dimension x(1), y(iy,1)
      dimension ds1(4,5)
c
      common/rhsdat/ el,ell,alb,els,el1,sigcom,anres,perfac,data(8),
     *  aa(iaa,1)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
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
      radial = el.le.1.e-6
c
c  initialize parameters for resetting nfit
c
      nftrcs=nftmax/10
      nftmx1=mod(nftmax,10)
c
c  initialize counter for trying different fitting points
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
c  carry out iteration
c
   10 call sigitr(sig,x,y,iy,nw1,nibc,nn,
     *  mdintg,nev1,iasn,ii,ig,ifrich,nfit,istsb1,iscan,
     *  xfit1,det,ddsig,ddsol,ds1,iord,icry,isolcv,isigcv,istatus)
c
c  test for convergence or retry with modified nfit
c
      if(isigcv.eq.1.and.isolcv.eq.1) then
	go to 50
      else if(nftcnt.le.nftmx1) then
c
c  retry iteration with new nfit
c
        call rsnfit(x,nibc,nn,nev1,nfit,xfit1,xfitst,nftcnt,nftrcs)
c
	sig=sigst
	istatus=3
	go to 10
      else
	write(istdou,215) nftmx1
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,215) nftmx1
	if(isigcv.eq.-1) then
	  istatus=1
        else
	  istatus=2
        end if
      end if
c
c  skip matching solution when using nrkint
c
   50 if(mdintg.ne.3) then
        call mchsol(x,y,iy,aa,iaa,nw1,nibc,nn,iasn,ii,ig,
     *    nfit,npout,imissl,imjssl,imstsl,det,ds1,icry)
      end if
c
c  multiply solution by (x/x(nev1))**(el-1) in interior
c
      if(nev1.gt.nibc.and..not.radial) then
        call sclasl(x,y,iy,nn,ii,nev1,nw1,el,iplneq)
      end if
c
      return
c
  215 format(//' ***** Convergence failed after ',i3,
     *         ' attempts to change nfit'/
     *         '       Give up.')
  220 format(/' Reset nfit to ',i5,' from temporary value ',i5)
      end
