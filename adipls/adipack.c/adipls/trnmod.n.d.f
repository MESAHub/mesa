      subroutine trnmod(x,nn0,xtrnct,ntrnsf,ntrnct,xtrnc1,ntrns1,
     *  nw1,nibc,nn,fulmod)
c
c  sets truncation point to truncate model at x = xtrnct, taking into
c  account the possible needs of Richardson extrapolation.
c  May also set parameters to truncate the model at ntrnsf points
c  from the surface.
c
c  The actual mesh point and location of the interior truncation
c  point are returned in ntrnct and xtrnc1, and the actual 
c  surface truncation point is returned in ntrns1.
c
c  Original version: 5/7/95
c
      implicit double precision(a-h,o-z)
      logical fulmod, singsf
      dimension x(1)
      common/rhsdat/ el,ell,alb,els,el1,sig,anres,perfac,data(8)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      ntrnct=1
      xtrnc1=0
      ntrns1=0
c
      singsf = data(7).ge.0
      if(singsf) then
	nns=nn-1
      else
	nns=nn
      end if
c
      if(ntrnsf.ge.1) then
	nns=nns-ntrnsf
	ntrns1=ntrnsf
      else
	ntrns1=0
      end if
c
c  locate inner truncation point
c
      if(xtrnct.gt.0) then
	do 20 n=1,nn
	ntrnct=n
	if(x(n).ge.xtrnct) go to 25
   20   continue
c
      end if
c
c  test for resetting truncation points with Richardson extrapolation
c
   25 if(iriche.eq.1.and.mod(nns-ntrnct,2).eq.1) then
c
c  when truncating at the surface, take up change there
c
        if(ntrnsf.ge.1) then
          ntrns1=ntrns1+1
c
c  otherwise change inner truncation point
c
        else
          ntrnct=ntrnct+1
        end if
      end if
c
      if(ntrnct.le.1.and.x(1).le.1.e-9) then
        nw1=1
        nibc=2
        fulmod=.true.
      else
        fulmod=.false.
c
c  truncation point
c
        xtrnc1=x(ntrnct)
c
c  reset first output point and point of inner boundary condition
c
        nw1=ntrnct
        nibc=ntrnct
      end if
c
c  truncation at surface?
c
      if(ntrns1.gt.0) nn=nn-ntrns1
c
      return
      end
