      subroutine rsnfit(x,nibc,nn,nev1,nfit,xfit1,xfitst,nftcnt,nftrcs)
c
c  resets fitting point, depending on nftcrs:
c  nftcrs = 0: reset based on meshpoint number
c  nftrcs = 1: reset based on mesh in x
c  reset alternates between decreasing and increasing and decreasing,
c  depending on nftcnt
c
c  original version: 30/4/02
c
      implicit double precision(a-h, o-z)
      dimension x(1)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      nftcnt=nftcnt+1
c
c  test for mode of change
c
      if(nftrcs.eq.0) then 
	if(istdpr.gt.0) 
     *  write(istdpr,*) 'nftcnt, orig. nfit:',nftcnt, nfit
        if(mod(nftcnt,2).eq.1) then
          nfit=max(10,nint(nfit-0.05*nftcnt*nn))
        else
          nfit=min(nint(nfit+0.05*nftcnt*nn),nn-10)
        endif
c
      else
c  
c  base change on mesh in x
c
	dx=0.05d0*min(xfit1,1.d0-xfit1)
	if(mod(nftcnt,2).eq.1) then
	  xfit1=xfit1-nftcnt*dx
        else
	  xfit1=xfit1+nftcnt*dx
        end if
c
	do 20 n=1,nn
	nfit=n
	if(x(n).gt.xfit1) go to 25
   20   continue
c
c  limit range of nfit
c
   25   nfit=min(nn-10,max(10,nfit))
c
      end if
c
c  test for resetting nfit 
c
      if(mod(nfit-nibc,2).ne.0) then
        nfit=nfit+1
      end if
      xfit1=x(nfit)
      write(istdou,140) nfit, xfitst, xfit1
      if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,140) 
     *  nfit, xfitst, xfit1
      return
  140 format(//' Retry iteration, changing nfit to', i6/
     *  ' Changing xfit from',f10.5, ' to',f10.5)
      end
