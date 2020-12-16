      subroutine sigint(x,y,iy,nw1,nibc,nn,mdintg,nev1,iasn,ii,ig,
     *  initnr,nfit,istsb1,isig,det,ddsol,ds1,iord,icry)
c
c  Integrates adiabatic oscillation equations for given sig and l
c  (transmitted in common/rhsdat/). Integration method is determined
c  by mdintg.
c
c  Returns matching condition in det, possibly test of eigenfunction
c  continuity in ddsol and matching determinant, to be used in
c  s/r mchsol, in ds1 (ddsol and ds1 are not set for mdintg = 3).
c
c  In case of errors, icry is returned as negative.
c
c  Original version: 28/6/95
c
      implicit double precision(a-h, o-z)
      parameter(iaa=10)
      dimension x(1), y(iy,1), ds1(4,1)
      common/rhsdat/ el,ell,alb,els,el1,sig,anres,perfac,data(8),
     *  aa(iaa,1)
      common/nrmchk/ nnww,ncfac,irsdif
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  set boundary conditions at centre and surface
c  *********************************************
c
      if(mdintg.eq.3) then
        istnrk=1
      else
        istnrk=0
      end if
c
c  set flag for taking out x**(l-1) behaviour at centre.
c
      if(nibc.eq.nev1.and.irsevn.lt.3) then
        icncbc=2
      else
        icncbc=1
      end if
c
c  set flag for taking out t**n behaviour near polytropic surface
c
      if(data(7).gt.0.and.data(7).lt.1.and.ii.eq.4
     *  .and.mdintg.eq.1) then
	anres=data(7)
      else
	anres=0
      end if
c
      call setbcs(x,y,iy,nn,nibc,
     *  icncbc,istsb1,ii,ig,nd,ny,isig,istnrk)
c
      anres=0
c
c                              --------------------------------------
c
c  the integration
c  ***************
c
c  test for using nrkm
c
      if(mdintg.eq.3) then
c
c  set normalization point
c
        ncnorm=nibc+fcnorm*(nn-nibc)
c
c  test for resetting nnww to zero to force nrkm to print actual storage
c  needed
c
        if(initnr.eq.1.and.nnwwin.eq.0) nnww=0
c
        call nrkint(x,y,iy,nw1,nibc,nn,nev1,iasn,ii,initnr,nfit,ncnorm,
     *    det,ddsol,icry)
c
        initnr=0
c
      else
c
c  use linint or eigint integration from centre and surface
c
        call shtint(x,y,iy,nw1,nibc,nn,nev1,nd,ny,iasn,ii,ig,
     *    nfit,mdintg,det,ddsol,ds1,iord,icry)
c
      end if
      return
      end
