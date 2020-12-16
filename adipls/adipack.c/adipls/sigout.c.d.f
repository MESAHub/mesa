      subroutine sigout(sig,x,y,iy,nw1,nibc,nn,nnw,mdintg,nev1,nfit,
     *  istsb1,iord,isolcv,isigcv)
c
c  Makes output to printer and files of results of adiabatic
c  oscillation calculation.
c  Includes setting of icase, for use in grand and short summaries.
c
c  Original version 1/7/95
c
c  Modified 26/4/99, replacing test on ig (which was undefined)
c
c  Modified 8/7/05, adding setting of mode quantities in
c  common /cobs_param/
c
      implicit double precision (a-h, o-z)
      logical radial, notwin, notwni, nscfil
      parameter(iaa=10, iaa1=10)
      dimension x(1), y(iy,1)
      dimension csummm(50),icsumm(8),ssummm(7),issumm(2),ssmmod(7)
c
      common/sysord/ ii
      common/worksp/ aa1(iaa1,1)
      common/rhsdat/ el,ell,alb,els,el1,sigcom,anres,perfac,data(8),
     *  aa(iaa,1)
      common/rotdat/ em, irotsl
      common/bcsdat/ fctsbc, fcttbc, istsbc, ibotbc
      common/varcon/ ivarf,npvarf,kvarfc
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
      common/csumma/ xmod1,datsum(8),datmd1(2),xtrnct,dum13,xfit1,
     *  fsbcsm,fcbcsm,albsum,
     *  elsum,ordsum,sigsum,sigc,dmx,xmx,ekin,per,perv,frqv,
     *  ddsig,ddsol,ysum(4),dosdum(5),
     *  in,nnwcom,mdints,ivarfs,icase,iorign,iekins,idum8,mlname(4)
      common/cmdtst/ iordtr,icaswn, sigwn1, sigwn2, frqwn1, frqwn2, 
     *  iorwn1, iorwn2, frlwn1, frlwn2
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
c
c  common for storage of modal parameters (degree, order, cyclic frequency,
c  inertia)
c
      common/cobs_param/ icobs_st, nobs_st, obs_st(10,1)
c
c  common controlling output
c
      common/coutpt/ nout,nprcen,iper,irotkr,nprtkr,igm1kr,npgmkr,
     *  nfmode,nfmesh,ispcpr,npout,nobs_stmx,nfmscn
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
      external rhsnrm,bcnrk
c
      equivalence (xmod1,csummm(1))
      equivalence (csummm(39),icsumm(1))
      equivalence (ssummm(6),issumm(1))
c
      radial = el.le.1.e-6
c
c  write squared frequency
c
      if(istdpr.gt.0) write(istdpr,110) sig
c
c  period
c
      per=perfac/sqrt(max(abs(sig),epsufl))
      if(istdpr.gt.0) write(istdpr,115) per
c
c  normalize to displacement = 1 at surface
c  ****************************************
c
      fcy=y(1,nn)
      if(.not.radial.and.icow.eq.0.and.abs(fcy).lt.1.e-5*abs(y(3,nn))) 
     *  fcy=y(3,nn)
      fcy1=fcy
      fcy=1
      if(fcy1.ne.0) fcy=1./fcy1
      if(abs(fcy-1).ge.1.e-8) then
c
        do 15 n=nw1,nn
        do 15 i=1,ii
   15   y(i,n)=fcy*y(i,n)
c
      end if
c
c  test for testing solution using nrkm right hand subroutine
c
      if(itssol.eq.1) then
c
c  test equations with standard equations throughout
c
        el1=0
	anres=0
c
        nnnrk=nn-nibc
        if(istdpr.gt.0) write(istdpr,120)
c
c  make test. use aa1 as dummy arrays for zk, ap, aq.
c
        call nrtssl(x(nibc),y(1,nibc),aa1,aa1,aa1,rhsnrm,bcnrk,ii,
     *    0,0,nnnrk,iy,nout,idgtss)
c
      end if
c
c  set energy normalized eigenfunction
c  (this badly needs clearing up)
c
      if(iplneq.eq.1) then
        do 20 n=nw1,nn
   20   y(6,n)=sqrt(aa(5,n))
      else
c
        do 25 n=nw1,nn
        y(6,n)=x(n)*x(n)*x(n)*aa(1,n)*aa(5,n)*aa(10,n)
        if(y(6,n).lt.0) then
c
c  write diagnostics for negative energy factor
c
          write(istdou,125) n,x(n),(aa(i,n),i=1,5),aa(10,n)
          if(istdpr.ne.istdou.and.istdpr.gt.0)
     *      write(istdpr,125) n,x(n),(aa(i,n),i=1,5),aa(10,n)
          y(6,n)=0
        else
c
          y(6,n)=sqrt(y(6,n))
        end if
c
   25   continue
c
      end if
c
      iwn=2
      if(radial) iwn=1
      fct=1
      do 30 i=1,iwn
      kwn=6+i
      if(i.eq.2) fct=1.d0/sqrt(ell)
      do 30 n=nw1,nn
   30 y(kwn,n)=fct*y(i,n)*y(6,n)
c
c  find maximum value
c  ******************
c
      call aramax(x(nw1),y(7,nw1),iy,nnw,xrmx,drmx,d2drmx)
      if(istdpr.gt.0) write(istdpr,130) drmx,xrmx
c
c  set values to grand summary
c
      csummm(34)=drmx
      csummm(35)=xrmx
c  normalize to maximum = 1
      if(drmx.eq.0) drmx=1
      fct=1.d0/drmx
      kwn=6+iwn
      do 35 k=7,kwn
      do 35 n=nw1,nn
   35 y(k,n)=fct*y(k,n)
c
      if(nout.gt.0.and.istdpr.gt.0) then
c
c  output
c  ******
c
	ndout=max0(1,nnw/nout)
c
        if(radial) then
          write(istdpr,135)
          iw1=2
          do 40 n=nw1,nn
          if(mod(n-nw1,ndout).eq.0.or.n.eq.nn.or.n.le.nprcen) 
     *      write(istdpr,140) n,x(n),(y(i,n),i=1,2),y(7,n)
   40     continue
c
        else
          write(istdpr,145)
          iw1=4
          do 45 n=nw1,nn
          if(mod(n-nw1,ndout).eq.0.or.n.eq.nn.or.n.le.nprcen) 
     *      write(istdpr,150) n,x(n),(y(i,n),i=1,4),y(7,n),y(8,n)
   45     continue
c
        end if
c
      end if
c
c  print order
c
      call order(x(nw1),y(1,nw1),data,aa(1,nw1),el,sig,icow,
     *  irsord,iy,iaa,iord,nnw,1,mdintg,-1)
      if(istdou.ne.istdpr) write(istdou,155) iord
c
c  find max/surface displacement
c
      call aramax(x(nw1),y(1,nw1),iy,nnw,xmx,dmx,d2dmx)
      if(istdpr.gt.0) write(istdpr,160)  dmx,xmx
c
c  find kinetic energy
c
      call kiner(x(nw1),y(1,nw1),aa(1,nw1),el,iy,iaa,nnw,ekin)
      if(istdpr.gt.0) write(istdpr,165) ekin
c
c  period and variational period
c
      if(iper.eq.1) then
        call varfrq(x(nw1),y(1,nw1),data,aa(1,nw1),sig,iy,iaa,nnw,
     *    istsb1,sigv)
        perv=perfac/sqrt(max(abs(sigv),epsufl))
        frqv=16.666666666667d0/perv
        if(istdpr.gt.0) write(istdpr,170) per,perv,frqv
        if(istdou.ne.istdpr) write(istdou,175) frqv
c
      end if
c
c                              --------------------------------------
c
c  test for file output, only in case of convergence
c
      if(isigcv.le.0) go to 80
c
c  set case number
c
      icase=10*iper+100*irotsl+10 000*istsb1+100 000*iplneq +
     *      1 000 000*iturpr
      if(icow.eq.2) then
        icase=icase+1
      else if(icow.eq.3) then
        icase=icase+2
      end if
      if(alb.ne.1.or.fctsbc.ne.0) icase=icase+1000
c
c  test for adding em (azimuthal order) to csummm
c
      if(irotsl.eq.1) then
	csummm(38)=em
      else
	csummm(38)=0.d0
      end if
c
c  test for windowing of modes for file output
c
      frmuhz=16666.666667d0/per
      flmuhz=frmuhz/(el+0.5)
      if(icaswn.ge.0.and.icase.ne.icaswn) then
	write(istdou,180) icase, icaswn
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,180) 
     *    icase, icaswn
	go to 80
      else if(notwin(sigwn1,sigwn2,sig)) then
	write(istdou,182) sig, sigwn1, sigwn2
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,182) 
     *    sig, sigwn1, sigwn2
	go to 80
      else if(notwin(frqwn1,frqwn2,frmuhz)) then
	write(istdou,184) frmuhz, frqwn1, frqwn2
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,184) 
     *    frmuhz, frqwn1, frqwn2
	go to 80
      else if(notwni(iorwn1,iorwn2,iord)) then
	write(istdou,186) iord, iorwn1, iorwn2
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,186) 
     *    iord, iorwn1, iorwn2
	go to 80
      else if(notwin(frlwn1,frlwn2,flmuhz)) then
	write(istdou,188) flmuhz, frlwn1, frlwn2
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,188) 
     *    flmuhz, frlwn1, frlwn2
	go to 80
      end if
c
c  set additional quantities in summary
c
      call store(data,datsum,8)
      fsbcsm=fctsbc
      fcbcsm=fcttbc
      albsum=alb
      elsum=el
      sigsum=sig
      ordsum=iord
      ivarfs=ivarf
      iekins=iekinr
c
c  test for storing dimensionless squared surface angular velocity
c
      if(irotsl.eq.1) then
        csummm(13)=1.5d0*(aa(10,nn)-1.d0)*aa(1,nn)
      end if
c
      do 50 i=1,4
   50 ysum(i)=y(i,nn)
c
      iorign=3
c
c  set mdints for summary, depending on whether Richardson extrapolation
c  was used
c
      if(iriche.eq.1) then
        mdints=mdintg+10
      else
        mdints=mdintg
      end if
c
c  test for rotational kernel
c
      if(irotkr.gt.0) call rotker(idsrkr,x(nw1),y(1,nw1),aa(1,nw1),
     *  el,sig,iy,iaa,nnw,nprtkr)
c
c  test for gamma1 kernel
c
      if(igm1kr.eq.1) call gm1ker(idsgkr,x(nw1),y(1,nw1),aa(1,nw1),
     *  el,sig,iy,iaa,nnw,npgmkr)
c
c  write mode on disc?
c
      if(((nfmode.ge.1.and.isolcv.ne.-1).or.nfmscn.eq.1)
     *    .and.nscfil(idsefn)) then
c
c  test for type of write
c
	if(nfmode.eq.1) then
          write(idsefn) csummm,nnw,(x(n),(y(i,n),i=1,4),(y(i,n),i=7,8),
     *      n=nw1,nn)
	else
	  if(nfmesh.ne.0) then
	    if(istdpr.gt.0) write(istdpr,*) 'nfmesh =',nfmesh
	    if(istdpr.gt.0) write(istdpr,*) 
     *        'write mode mesh, nw1, nn =',nw1, nn
	    write(idsefn) nnw,(x(n),n=nw1,nn)
	    nfmesh=0
	  end if
	  if(nfmode.eq.2) then
            write(idsefn) csummm,((y(i,n),i=1,2),n=nw1,nn)
          else
            write(idsefn) csummm,((y(i,n),i=7,8),n=nw1,nn)
	    if(istdpr.gt.0) write(istdpr,*) 
     *        'write amde, nw1, nn =',nw1, nn
	  end if
        end if
        call flush(idsefn)
        idsn=4
        if(istdpr.gt.0) write(istdpr,195) idsn, nfmode
      end if
c
c  output grand summary
c
      if(idsgsm.gt.0.and.nscfil(idsgsm)) then
        write(idsgsm) csummm
        call flush(idsgsm)
      end if
c
c  set and output short summary
c
      if(idsssm.gt.0.and.nscfil(idsssm)) then
        call setssm(csummm,icsumm,ssummm,issumm,ssmmod,irsmod)
        if(irsmod.eq.1) write(idsssm) ssmmod
        write(idsssm) ssummm
        call flush(idsssm)
      end if
c
c  set mode quantities in common/cobs_param/
c
      call setobs_st(csummm,icsumm,iriche)
c
c  test for special output
c
      if(ispcpr.ne.0) call spcout_adi(x,y,aa,data,nn,iy,iaa,ispcpr)
c
   80 continue
c
      return
  110 format(//' squared frequency =',f17.10)
  115 format(/' period =',1pe13.5,'  minutes')
  120 format(//' test solution with s/r rhsnrk:'/)
  125 format(' ***** negative energy factor. n, x, aa(1-5,10) =',
     *  i5,f10.5,1p6e11.3)
  130 format(//' maximum of energy normalized vertical displacement is',
     *  1pe13.5,' and occurs at x =',0pf12.7//
     *  ' renormalized to have maximum value 1')
  135 format(///' complete solution. n,x,y(1-2),energy normalized',
     *  ' vertical displacement:'/)
  140 format(i5,f11.7,1p2e13.5,5x, e13.5)
  145 format(///' complete solution. n,x,y(1-4),energy normalized',
     *  ' vertical and horizontal displacements:'/)
  150 format(i5,f11.7,1p4e13.5,5x,2e13.5)
  155 format(' mode order =',i4)
  160 format(/' max/surface displacement =',1pe13.5,
     *  ' and is found at x =',0pf10.7)
  165 format(/' normalized dimensionless kinetic energy =',
     *  1pe13.5)
  170 format(/' period from eigenvalue    =',f12.6,' minutes'//
     *        ' period from var. integral =',f12.6,' minutes'//
     *        ' freqv. from var. integral =',f12.6,' mhz')
  175 format(' Variational frequency =',1pe13.5,' mHz')
  180 format(/' ***** Warning: icase =',i10,' .ne. target case =',i10/
     *        '       No output to file.')
  182 format(/' ***** Warning: sigma**2 =',1pe13.5,' not in window',
     *                2e13.5/
     *        '       No output to file.')
  184 format(/' ***** Warning: nu (microHz) =',1pe13.5,' not in window',
     *                2e13.5/
     *        '       No output to file.')
  186 format(/' ***** Warning: order =',i5,' not in window',2i5/
     *        '       No output to file.')
  188 format(/' ***** Warning: nu/L (microHz) =',1pe13.5,
     *                ' not in window',2e13.5/
     *        '       No output to file.')
  195 format(//'  mode written on disc, dsn =',i4,' for case ',i2)
      end
