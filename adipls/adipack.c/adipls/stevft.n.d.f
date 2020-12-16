      subroutine stevft(x,nn,nibc,nfit,iasn,radial,fulmod,nev1,xfit1)
c
c  take out factor x**l in solution up to point nev1 where solution
c  begins to vary more slowly than x**(l/2).
c  If irsevn .lt. 2 point may be restricted by nfit if nfit .le. nev1.
c  If irsevn .ge. 2, reset nfit to be just outside nev1, in that case
c
c  For truncated model, radial oscillations or if irsevn = -1
c  set point to central mesh point.
c
c  note. this was implemented on 24/4/83.
c
c  before this date nev1 was taken to be at the first boundary of
c  the evanescent region.
c
c  Note that routine assumes that nfit has been set
c  before call, but may reset it, depending on xfit
c  Also sets nev1, and possibly xfit1
c
c  original version: 13/2/89
c
c  Modified 15/5/91, introducing option of not using evanescent
c  region for irsevn = -1.
c
c  Modified 23/7/94, to allow resetting of nfit to be outside
c  evanescent region.
c
c  Modified 28/6/95, removing nfit1. Remove frfit1, replace test on frfit1
c  by test on xfit, now in common/cincnt/. Move aa, data to common /rhsdat/
c
c  Modified 19/7/95, introducing parameter xmnevn (in common/cincnt/)
c  Search for evanescent region is restricted to x .ge. xmnevn
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      logical radial, fulmod
      parameter(iaa = 10)
      dimension x(1)
      common/rhsdat/ el,ell,alb,els,el1,sig,anres,perfac,data(8),
     *  aa(iaa,1)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(radial.or.irsevn.eq.-1) then
        nev1=nibc
      else
        nev1=nfit
c
c  if fitting point is very close to centre, do not
c  locate evanescent point.
c
        if(nfit.gt.5.or.xfit.ne.-1) then
c
          if(fulmod) then
            nsev=5
          else
            nsev=nibc
          end if
          if(xfit.eq.-1.or.irsevn.ge.2) then
            nlev=nn
          else
            nlev=nfit
          end if
c
c  initialize flag for finding boundary of evanescent region
c
          ifnev = 0
          do 20 n=nsev,nlev
          nev1=n
          akr=(1-sig*aa(2,n)/(aa(1,n)*ell))*(1-aa(1,n)*aa(4,n)/sig)
          if(akr.lt.0.25.and.x(n).ge.xmnevn) then
            ifnev = 1
            go to 25
          end if
   20     continue
c
        end if
c
      end if
c
   25 if(nev1-nibc.le.5) nev1=nibc
c
c  test for nev1 being too close to surface (mainly relevant for
c  surface-truncated model)
c
      if(nev1.ge.nn-30) then
        write(istdou,110) nev1, nn-30
        if(istdpr.gt.0.and.istdpr.ne.istdou)
     *    write(istdpr,110) nev1, nn-30
        nev1=nn-30
      end if
c
c  for iasn.gt. 1 (i.e. when doing Richardson extrapolation)
c  reset nev1 to match step iasn
c
      if(mod(nev1-nibc,iasn).ne.0) then
        nstp=(nev1-nibc)/iasn+1
        nev1=nibc+nstp*iasn
      end if
c
c  for xfit = -1 set nfit to nev1
c
      if(xfit.eq.-1) then
        nfit=nev1
        xfit1=x(nfit)
        if(istdpr.gt.0) write(istdpr,120) nfit
      end if
c
      if(nev1.gt.nibc) then
        if(istdpr.gt.0) write(istdpr,130) x(nev1)
        if(nev1.eq.nfit.and.xfit.ne.-1.and.ifnev.eq.0) then
          if(istdpr.ne.istdou) then
            write(istdou,130) x(nev1)
            write(istdou,135)
          end if
          if(istdpr.gt.0) write(istdpr,135)
	else if(irsevn.ge.2.and.nev1.gt.nfit) then
	  nfit=nev1
	  xfit1=x(nfit)
	  write(istdou,150) nfit, xfit1
	  if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,150) 
     *      nfit, xfit1
        end if
c
      else
c
        if(istdpr.gt.0) write(istdpr,140)
c
      end if
c
      return
  110 format(//' ***** Warning in s/r stevft. nev1 =',i5,
     *  ' too close to surface'/
     *         '       nev1 reset to',i5/)
  120 format(///' xfit = -1. nfit reset to',i5,' in s/r stevft')
  130 format(//' factor x**(l-1) taken out up to x =',f10.6/)
  135 format(' *** Warning: this was restricted by the fitting point.')
  140 format(//' integration of ordinary equations on whole',
     *  ' interval'/)
  150 format(/
     *  ' *** Warning. fitting point reset to be at transition point'/
     *  '     New nfit, xfit =',i5,f10.5)
      end
