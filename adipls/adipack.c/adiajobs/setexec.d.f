      program main
c  set trial frequencies etc for adiabatic programme
c
c  based on grand summary or short summary of modes read from unit 2.
c
c  outputs short summary, containing trial frequencies,
c  on unit 3. this has only been implemented for ifill = 1,
c  and only for a single model.
c
c  namelist /exec/ is output on unit 4, in format ready
c  for adiabatic pulsation programme, apart from
c  global parameters controlling integration etc.
c  note: this is only relevant for versions of the pulsation
c  programme that use namelist.
c
c  when ifill .ne. 1 sets exec with nsel, dels, for modes of order
c  n between n1 and n2.
c  for l2 .le. 0 extrapolates from largest two l-values for each
c  n, and nsel is taken from input, possibly reduced such that
c  sig does not exceed sigmax.
c  for l2 .gt. 0 interpolates between l1 and l2, and nsel is
c  set by programme, as determined by l1, l2 and dels.
c
c                    -----------------------------
c
c  when ifill = 1 sets exec to fill gaps in summary, assuming that
c  this is ordered increasingly with l first.
c
c  when l2. ge. 0 only l-values between l1 and l2 are included.
c  when nfill1 .lt. nfill2 only modes with order between nfill1
c  and nfill2 are included.
c  when sigmin .gt. 0 only modes with sig .ge. sigmin are included.
c  when sigmax .gt. 0 only modes with sig .le. sigmax are included.
c
c  for sigmin .gt. 0 also extrapolates to get all sig down to sigmin.
c  for sigmax .gt. 0 also extrapolates to get all sig up to sigmax.
c
c  step is uniform in sig, sqrt(sig) or 1./sqrt(sig), depending on
c  whether nsig = 1, 2 or 3.
c
c  when ielfl1 = 1 eliminates f(l=1) mode
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      logical interp
      dimension cs(50),el(2,100),sig(2,100),ns(3),els(3),sigs(3),
     *  csmp(8),ss(7),iss(2)
      equivalence(ss(6),iss(1))
c
c  this common is probably not needed.
c..      common /cmexec/ n1,n2,l1,l2,nfill1,nfill2,
c..     *  nsel,dels,nsig,sigmax,ifill,
c..     *  sigmin,ielfl1
c..      namelist/exec/ icase,ifill,l1,l2,nfill1,nfill2,ielfl1,
c..     *  sigmin,sigmax,
c..     *  nsig,n1,n2,nsel,dels
c
      data csmp /8*-1.d0/
      data lp /-1/
c
c  icase: input type.
c  icase = 1: grand summary
c  icase = 2: short summary
      icase=1
c
      ifill=1
      l1=0
      l2=0
      nfill1=0
      nfill2=-1
      ielfl1=1
      sigmin=1
      sigmax=2500
      nsig=1
      n1=1
      n2=2
      nsel=20
      dels=1
c
c  set files
c
      write(6,*) 'Files needed:'
      write(6,*) 'Input grand summary on unit 2'
      write(6,*) 'Output short summary on unit 3'
      write(6,*) 'Output namelist on unit 4'
c
      call ofiles
      call openf(2,'o','u')
      call openf(3,'u','u')
      call openf(4,'u','f')
      open(12,file='ttt.setex.scan',status='unknown')
      open(21,file='setexec.log',status='unknown')
c
c  set iss (note: convention for iorign has yet to be established)
c
      iss(1)=0
      iss(2)=-5
c
   10 continue
c..      read(5,exec,end=90,err=90)
      write(6,*) 'icase,ifill,l1,l2,nfill1,nfill2,ielfl1,?'
      write(6,*) icase,ifill,l1,l2,nfill1,nfill2,ielfl1
      read(5,*,end=90,err=90) icase,ifill,l1,l2,nfill1,nfill2,ielfl1
      write(6,*) 'sigmin,sigmax,?'
      write(6,*) sigmin,sigmax
      read(5,*,end=90,err=90) sigmin,sigmax
      write(6,*) 'nsig,n1,n2,nsel,dels?'
      write(6,*) nsig,n1,n2,nsel,dels
      read(5,*,end=90,err=90) nsig,n1,n2,nsel,dels
c
c..      write(6,100)
c
      interp=l2.gt.0
      if(interp) nsel=max(0.5d0+dfloat(l2-l1-1)/dels,1.d0)
c
c..      write(6,exec)
c
      write(6,*) 'icase,ifill,l1,l2,nfill1,nfill2,ielfl1'
      write(6,*) icase,ifill,l1,l2,nfill1,nfill2,ielfl1
      write(6,*) 'sigmin,sigmax'
      write(6,*) sigmin,sigmax
      write(6,*) 'nsig,n1,n2,nsel,dels'
      write(6,*) nsig,n1,n2,nsel,dels
c
      if(sigmax.le.0) go to 14
c
c  set sigmx1
c
      go to (11,12,13), nsig
   11 sgmx1=sigmax
      go to 14
   12 sgmx1=sqrt(sigmax)
      go to 14
   13 sgmx1=1.d0/sqrt(sigmax)
c
c  if sigmin .gt. 0 set sgmn1
c
   14 if(sigmin.le.0) go to 15
c
 1410 go to (1411,1412,1413), nsig
c
 1411 sgmn1=sigmin
      go to 15
 1412 sgmn1=sqrt(sigmin)
      go to 15
 1413 sgmn1=1.d0/sqrt(sigmin)
c
   15 nmodel=0
      ifind=0
      rewind 2
c
c  entry point for new model
c
15500 do 16 k=1,100
      do 16 i=1,2
   16 el(i,k)=-1
c
      if(ifill.ne.1) go to 20
c
c  initialize counting indices and flags for ifill = 1
c
      is=0
      isp=0
      iexmin=0
      nmode=0
      write(6,150)
c
c
c  read mode data
c  **************
c
   20 inmod=0
c
      call rdfreq(icase,2,cs,l,n,s,frq,ekin,ierr)
c..   write(6,*) 'read mode, n, l, sigma, frq  =',n,l,s,frq
      if(ierr.ne.0) go to 30
c
c  test for new model
c
      if(icase.eq.1) then
        do 21 k=2,5
        if(abs(cs(k)/csmp(k)-1).gt.1.d-6) inmod=1
   21   csmp(k)=cs(k)
c
      else if(l.lt.0) then
        do 22 k=3,6
        if(abs(cs(k)/csmp(k)-1).gt.1.d-6) inmod=1
   22   csmp(k)=cs(k)
c
c  in this case, read next frequency record
c
        call rdfreq(icase,2,cs,l,n,s,frq,ekin,ierr)
        if(ierr.ne.0) go to 30
c
      end if
c
      if(inmod.eq.1) then
c
c  new model, test for setting ifind
c
        nmodel=nmodel+1
c
        ifindp=ifind
        ifind=2
c
        if(nmodel.gt.1) then
          backspace 2
          go to 30
        end if
c
      end if
c
      iendmd=0
c
      if(ifill.eq.1) go to 27500
c
c  ----------------------------------------------------------------
c
c  set arrays for ifill .ne. 1
c
      if(n.lt.n1.or.n.gt.n2) go to 20
c
      k=n-n1+1
c  test for interpolation
      if(interp) go to 25
c  test for largest l at given n
      if(l.le.el(2,k)) go to 20
      el(1,k)=el(2,k)
      el(2,k)=l
      sig(1,k)=sig(2,k)
      sig(2,k)=cs(20)
      go to 20
c  interpolation between l1 and l2
   25 if(l.ne.l1) go to 27
      el(1,k)=l
      sig(1,k)=cs(20)
      go to 20
   27 if(l.ne.l2) go to 20
      el(2,k)=l
      sig(2,k)=cs(20)
      go to 20
c
c  ----------------------------------------------------------------
c
c  ifill = 1: fill gaps in summary
c
c  test for l in range
c
27500 continue
c..   write(6,*) 
c..  *  'Enter testing section, l1, l2, l, nfill1, nfill2, n =',
c..  *  l1, l2, l, nfill1, nfill2, n
      if(l2.ge.0) then
        if(l.lt.l1) go to 20
      end if
c
c  test for n in range
c
      if(nfill1.lt.nfill2) then
        if(n.lt.nfill1.or.n.gt.nfill2) go to 20
      end if
c
c  test for initial point
c
c..   write(6,*) 'test for initial point, l, n, is =',l, n, is
      if(is.eq.0) go to 2950
c
c  test for s in range
c
c..   write(6,*) 'test for sigma. sigmin, sigmax, s =',sigmin, sigmax,
c..  *  s
      if(sigmin.gt.0.and.s.lt.sigmin) then
        iexmns=1
        go to 20
      end if
c
      if(sigmax.gt.sigmin.and.s.gt.sigmax) then
        iexmxs=1
        go to 20
      end if
c
c  test for non-monotonic order, skip identical points
c
      if(l.eq.lp) then
        if(n.eq.np) then
          go to 20
        else if(n.lt.np) then
          write(6,105) l,n
          go to 20
        end if
      end if
c  test for new l
c..   write(6,*) 'test on l, lp, iexmin =',l, lp, iexmin
      if(l.ne.lp) go to 29
c  test for extrapolating to sigmin (previous l-value was new)
      if(iexmin.eq.1) go to 2930
c  entry after extrapolating to sigmin. reset iexmin
 2805 iexmin=0
c
c  test for gap in order
c
c..   write(6,*) 'testing for gap, l, np, n =',l, np, n
      if(n-np.le.1) go to 2950
      istsig=n-np-1
c  test for skipping f(l=1) mode
      if(l.eq.1.and.np.lt.0.and.n.gt.0.and.ielfl1.eq.1) then
        istsig=istsig-1
        if(istsig.eq.0) go to 2950
      end if
c
      istord=1
      nord1=np+istord
c
c  set step and sig1, depending on nsig
c
      go to (2811,2812,2813), nsig
 2811 dfsig=(s-sp)/(n-np)
      sig1=sp+dfsig
      go to 2820
 2812 ssp=sqrt(sp)
      dfsig=(sqrt(s)-ssp)/(n-np)
      sig1=ssp+dfsig
      sig1=sig1*sig1
      go to 2820
 2813 ssp=1.d0/sqrt(sp)
      dfsig=(1.d0/sqrt(s)-ssp)/(n-np)
      sig1=ssp+dfsig
      sig1=1.d0/(sig1*sig1)
c  output
 2820 write(6,155) l,np,n,sig1,dfsig,istsig
      sigs1=sp
      sigs2=s
      nns=iabs(n-np)
      go to 2940
c
c  new l, test for extrapolating to sigmax
c
   29 iexmax=0
      if(sigmax.le.sp.or.iexmxs.eq.1) go to 2925
c  test for same l
      if(els(is).ne.els(isp)) then
c
c  diagnostics for different l-values in extrapolation
c
        write(6,180) lp,ns(is)
        iexmin=0
        go to 2950
      end if
c
c  set step, sig1 and istsig
c
 2910 se=sigs(is)
      sep=sigs(isp)
      ne=ns(is)
      nep=ns(isp)
c
      istord=isign(1,ne-nep)
      nord1=ne+istord
c
      go to (2911,2912,2913), nsig
 2911 dfsig=(se-sep)/(ne-nep)
      sig1=se+dfsig
      istsig=1+(sgmx1-se)/dfsig
      go to 2920
 2912 sse=sqrt(se)
      dfsig=(sse-sqrt(sep))/(ne-nep)
      sig1=sse+dfsig
      sig1=sig1*sig1
      istsig=1+(sgmx1-sse)/dfsig
      go to 2920
 2913 sse=1.d0/sqrt(se)
      dfsig=(sse-1.d0/sqrt(sep))/(ne-nep)
      sig1=sse+dfsig
      sig1=1.d0/(sig1*sig1)
      istsig=1+(sse-sgmx1)/dfsig
c
 2920 istsig=max0(1,istsig)
c
c  test for restriction when nfill1 .lt. nfill2
c
      if(nfill1.lt.nfill2) istsig=min0(istsig,nfill2-ne)
c
c  test for frequency exceeding sigmax
c
      if(sig1.gt.sigmax.or.istsig.eq.0) then
        if(iendmd.eq.1) then
          go to 30500
        else
          go to 2925
        end if
      end if
c
c
      write(6,160) lp,nep,ne,sig1,dfsig,istsig
      iexmax=1
c
c  set exec
c
      elw=lp
c
c  test for final extrapolation
c
      if(iendmd.eq.1) then
c
c  test for setting ifind, and outputting model record
c
        if(ifindp.eq.0) then
          ifindw=-1
          write(4,170) ifindw,elw,sig1,dfsig,istsig,nsig
        else
          ifindw=2
          xmod=nmodel
          if(inmod.eq.1) xmod=xmod-1
          write(4,175) ifindw,xmod,elw,sig1,dfsig,istsig,nsig
          ss(1)=-1
          ss(2)=0
          ss(3)=cs(2)
          ss(4)=cs(3)
          ss(5)=cs(4)
          write(3) ss
        end if
c
        call wrtssm(3,ss,elw,sig1,dfsig,istsig,nsig,nord1,istord,
     *    sigmin,sigmax,nfill1,nfill2)
	nmode=nmode+1
c
        go to 30500
c
      else
c
c  test for setting ifind
c
        if(ifind.eq.0) then
          ifindw=-1
          write(4,170) ifindw,elw,sig1,dfsig,istsig,nsig
        else
          ifindw=2
          xmod=nmodel
          ifind=0
          write(4,175) ifindw,xmod,elw,sig1,dfsig,istsig,nsig
          ss(1)=-1
          ss(2)=0
          ss(3)=cs(2)
          ss(4)=cs(3)
          ss(5)=cs(4)
          write(3) ss
        end if
c
        call wrtssm(3,ss,elw,sig1,dfsig,istsig,nsig,nord1,istord,
     *    sigmin,sigmax,nfill1,nfill2)
	nmode=nmode+1
      end if
c
c  new l. test for extrapolation to sigmin
c  (extrapolation is only done after next mode has been read in)
c
 2925 iexmin=0
      if(sigmin.gt.0.and.s.gt.sigmin.and.iexmns.ne.1.and.
     *  (l.gt.0.or.n.gt.1)) iexmin=1
c
      go to 2950
c
c  entry point for extrapolation to sigmin
c
c  test for same l-value
c
 2930 if(l.ne.lp) then
c
c  diagnostics for different l-values in extrapolation
c
        write(6,180) lp,ns(is)
        iexmin=0
        go to 2950
      end if
c
c  set step, sig1, istsig
c
      se=s
      sep=sp
      ne=n
      nep=np
c
      istord=isign(1,nep-ne)
      nord1=nep+istord
c
      go to (2931,2932,2933), nsig
c
 2931 dfsig=(se-sep)/(nep-ne)
      sig1=sep+dfsig
      istsig=1+(sgmn1-sep)/dfsig
      go to 2935
 2932 ssep=sqrt(sep)
      dfsig=(sqrt(se)-ssep)/(nep-ne)
      sig1=ssep+dfsig
      sig1=sig1*sig1
      istsig=1+(sgmn1-ssep)/dfsig
      go to 2935
 2933 ssep=1.d0/sqrt(sep)
      dfsig=(ssep-1.d0/sqrt(se))/(ne-nep)
      sig1=ssep+dfsig
      sig1=1.d0/(sig1*sig1)
      istsig=1+(sgmn1-ssep)/dfsig
c
 2935 istsig=max0(1,istsig)
c
c  test for restriction when nfill1 .lt. nfill2
c
      if(nfill1.lt.nfill2) istsig=min0(istsig,np-nfill1)
c
      if(sig1.lt.sigmin.or.istsig.eq.0) go to 2945
c
      write(6,165) l,np,n,sig1,dfsig,istsig
c
c  set exec
c
 2940 elw=lp
c
c  test for setting ifind
c
      if(ifind.ne.2) then
        ifindw=-1
        write(4,170) ifindw,elw,sig1,dfsig,istsig,nsig
      else
        ifindw=2
        xmod=nmodel
        ifind=0
        write(4,175) ifindw,xmod,elw,sig1,dfsig,istsig,nsig
        ss(1)=-1
        ss(2)=0
        ss(3)=cs(2)
        ss(4)=cs(3)
        ss(5)=cs(4)
        write(3) ss
      end if
c
      call wrtssm(3,ss,elw,sig1,dfsig,istsig,nsig,nord1,istord,
     *  sigmin,sigmax,nfill1,nfill2)
      nmode=nmode+1
c
c  output partial input for scan
c
      write(12,172) elw,sigs1,nsig,20*nns,sigs2
c
c  go back after extrapolating to sigmin
c
c  note that logic has not been fully thought through here
c
 2945 if(iexmin.eq.1) go to 2805
c
c  end for this mode
c  *****************
c
c  test for l outside range
c
 2950 if(l2.gt.0.and.l.gt.l2) go to 10
c
c  for initial mode, test for extrapolation to sigmin
c
      if(is.eq.0.and.sigmin.gt.0) then
        iexmin=0
        if(s.gt.sigmin.and.(l.gt.0.or.n.gt.1)) iexmin=1
      end if
c
c  store data
c
      lp=l
      np=n
      sp=s
      isp=is
      is=1+mod(is,3)
      ns(is)=n
      els(is)=l
      sigs(is)=s
c..   write(6,*) 'lp, np, isp, is set to ',lp, np, isp, is
c
c  initialize flags for suppressing extrapolation when appropriate mode
c  is already there, but is removed by the restriction of the range of s
c
      iexmns=0
      iexmxs=0
c
      go to 20
c
c  new model, or end scan through summary
c
   30 if(ifill.ne.1) go to 31
c
c  test for final extrapolation to sigmax
c
      if(s.lt.sigmax.and.iexmxs.eq.0.and.els(is).eq.els(isp)) then
        if(inmod.eq.0) ifindp=ifind
        lp=l
        iendmd=1
        go to 2910
      end if
c
c  entry after final extrapolation
c
30500 if(inmod.eq.0) then
        go to 10
      else
        go to 15500
      end if
c
c  ---------------------------------------------------------------
c
c  output for ifill .ne. 1
c  ***********************
c
   31 nntot=n2-n1+1
      write(6,110) (el(1,k),el(2,k),sig(1,k),sig(2,k),k=1,nntot)
c  set exec for output
      write(6,120) nsig
      i1=2
      if(interp) i1=1
      do 40 k=1,nntot
      if(el(1,k).lt.0.or.el(2,k).lt.0) go to 40
c  set dfsig and sig1, depending on nsig
      ddel=dels/(el(2,k)-el(1,k))
      go to (32,34,36), nsig
   32 ssig=sig(2,k)
      dfsig=(ssig-sig(1,k))*ddel
      if(interp) ssig=sig(1,k)
      sig1=ssig+dfsig
      go to 38
   34 ssig=sqrt(sig(2,k))
      ssig1=sqrt(sig(1,k))
      dfsig=(ssig-ssig1)*ddel
      if(interp) ssig=ssig1
      sig1=ssig+dfsig
      sig1=sig1*sig1
      go to 38
   36 ssig=1.d0/sqrt(sig(2,k))
      ssig1=1.d0/sqrt(sig(1,k))
      dfsig=(ssig-ssig1)*ddel
      if(interp) ssig=ssig1
      sig1=1.d0/(ssig+dfsig)
      sig1=sig1*sig1
c
   38 els1=el(i1,k)+dels
c
c  reset nsel to limit sig to less than sigmax
c
      nsel1=nsel
      if(sigmax.le.0) go to 39
      nsel1=1+(sgmx1-ssig)/dfsig
      nsel1=min0(nsel1,nsel)
c
   39 n=k+n1-1
      write(6,130) n,els1,sig1,dfsig,nsel1
      if(nsel1.eq.0) go to 40
c
c  test for setting ifind
c
      if(ifindp.ne.2) then
        ifindw=-1
        write(4,140) ifindw,els1,dels,nsel1,nsig,sig1,dfsig
      else
        ifindw=2
        xmod=nmodel-1
        ifindp=0
        write(4,145) ifindw,xmod,els1,dels,nsel1,nsig,sig1,dfsig
      end if
   40 continue
c
c  test for new model
c
      if(inmod.eq.0) then
        go to 10
      else
        go to 15500
      end if
c
   90 continue
      write(21,'(i5)') nmode
      write(6,200)
      stop
  100 format(////)
  105 format(/' ***** non-monotonic mode order at degree, order =',
     *  2i5)
  110 format(//' el(1-2),sig(1-2):'//(2f10.2,2f12.5))
  120 format(///' set trials with nsig =',i3//' n,els1,sig1,dfsig,',
     *  'nsel:'/)
  130 format(i4,f10.2,2f12.6,i4)
  140 format(' $exec ifind=',i2,',',/
     *  ' els1=',1pe13.5,',dels=',e13.5,',nsel=',i3,','/
     *  ' nsig=',i1,',sig1=',e15.7,',dfsig=',e15.7,',$end')
  145 format(' $exec ifind=',i2,',xmod=',f6.1,',',/
     *  ' els1=',1pe13.5,',dels=',e13.5,',nsel=',i3,','/
     *  ' nsig=',i1,',sig1=',e15.7,',dfsig=',e15.7,',$end')
  150 format(///' fill gaps in order'//' l, np, n, sig1, dfsig,',
     *  ' istsig:'/)
  155 format(/3i6,1p2e13.5,i5,'  fill gap')
  160 format(/3i6,1p2e13.5,i5,'  extrapolate to sigmax')
  165 format(/3i6,1p2e13.5,i5,'  extrapolate to sigmin')
  170 format(' $exec ifind=',i2,',',
     *  'el=',1pe13.5,',sig1=',e13.5,
     *  ',dfsig=',e13.5,','/' istsig=',i3,',nsig=',i1,
     *  ',itrsig=1,$end')
  172 format('cntrd:'/
     *  'osc     @'/
     *  'osc:'/
     *  'el,nsel,els1,dels,'/
     *  f10.3,',0,0,1,,,,,,,,,     @'/
     *  'itrsig,sig1,istsig,inomde,itrds,'/
     *  '1, ',1pe13.5,' ,    ,    1,10,,,,,,,,     @'/
     *  'dfsig,nsig,iscan,sig2,'/
     *  ',',2i5,e13.5,',,,,,,,,,,,,     @')
  175 format(' $exec ifind=',i2,',xmod=',f6.1,','/
     *  ' el=',1pe13.5,',sig1=',e13.5,
     *  ',dfsig=',e13.5,','/' istsig=',i3,',nsig=',i1,
     *  ',itrsig=1,$end')
  180 format(/2i6,37x,'  only one n-value for this l *****')
  200 format(1h1)
      end
      subroutine wrtssm(ids,ss,elw,sig1,dfsig,istsig,nsig,nord1,istord,
     *  sigmin,sigmax,nfill1,nfill2)
c
c  write short summary for modes corresponding to
c  sig1, dfsig, istsig, nsig
c
c  when sigmin .lt. sigmax only output modes with sig between
c  sigmin and sigmax
c  when nfill1 .lt. nfill2 only output modes with order between
c  nfill1 and nfill2
c
      implicit double precision (a-h, o-z)
      parameter (istmax =100)
      dimension ss(*),sigst(istmax),nordst(istmax)
c
      write(6,*) 'enter wrtssm with sig1 =',sig1,'  dfsig =',dfsig
c
      if(istsig.eq.0) return
c
      if(istsig.gt.istmax) then
        write(6,100) istsig, istmax
        istsig=istmax
      end if
c
      sigst(1)=sig1
      nordst(1)=nord1
c
      sigp=sig1
      nordp=nord1
c
      do 20 ist=2,istsig
c
      nordst(ist)=nordp+istord
c
c  find sig from sigp and dfsig
c
      if(nsig.eq.1) then
c  step in sig
        signew=sigp+dfsig
      else
c
        if(sigp.ge.0) then
          isig=1
        else
          isig=-1
        end if
        sig=sqrt(isig*sigp)
	write(6,*) 'sig, sigp =',sig,sigp
c
        if(nsig.eq.2) then
c  step in sqrt(sig)
          sig=sig+dfsig
        else
c  step in 1/sqrt(sig)
          sig=sig/(1+dfsig*sig)
        end if
c
        signew=isig*sig*sig
        write(6,*) 'signew =',signew
c
      end if
      sigp=signew
      nordp=nordst(ist)
   20 sigst(ist)=signew
c
c  output in order of increasing frequencies
c
      ss(1)=elw
      ss(4)=0
      ss(5)=0
c
      if((dfsig.gt.0.and.nsig.le.2).or.(dfsig.lt.0.and.nsig.eq.3)) then
        do 30 ist=1,istsig
        sig=sigst(ist)
        nord=nordst(ist)
        if((sigmin.ge.sigmax.or.(sig.ge.sigmin.and.sig.le.sigmax)) 
     *  .and.(nfill1.ge.nfill2.or.(nord.ge.nfill1.and.nord.le.nfill2))) 
     *    then
          ss(2)=nordst(ist)
          ss(3)=sigst(ist)
          write(ids) (ss(i),i=1,7)
        end if
   30   continue
      else
        istsg1=istsig+1
        do 40 ist=1,istsig
        sig=sigst(istsg1-ist)
        nord=nordst(istsg1-ist)
        if((sigmin.ge.sigmax.or.(sig.ge.sigmin.and.sig.le.sigmax)) 
     *  .and.(nfill1.ge.nfill2.or.(nord.ge.nfill1.and.nord.le.nfill2)))
     *    then
          ss(2)=nord
          ss(3)=sig
          write(ids) (ss(i),i=1,7)
        end if
   40   continue
c
      end if
c
      return
  100 format(//' **** warning. istsig = ',i5,' greater than istmax =',
     *  i5/
     *  16x,'istsig has been reset to istmax')
      end
