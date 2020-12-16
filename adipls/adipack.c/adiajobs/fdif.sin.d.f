      program main
c
c  differences between adiabatic frequencies in single grand summary.
c
c  modified 18/7/89, to allow approximate resetting of inertia scaling
c  to a normalization point other than the surface.
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
c  Modified 26/4/99, fixing various minor inconsistencies.
c
      implicit double precision (a-h, o-z)
      parameter(nmdmax=5000)
      character*280 file
      character tail*5, cnutyp*40
      dimension cs(50),itype(2),
     *  l(nmdmax),nord(nmdmax),frq(2,nmdmax),
     *  data(8,20),nm(20),datap(8),
     *  ekin(nmdmax),cnutyp(4)
      common/cofile/ nfiles, idsfil(20), file(20), iopen(20)
c
      data ifsrun /1/
      data idsinp /-1/
      data kmaxpl /0/
      data tail /'    @'/
      data cnutyp /
     *  ' (possibly corrected) eigenfrequency',
     *  ' (uncorrected) eigenfrequency',
     *  ' variational frequency',
     *  ' Richardson extrapolation frequency'/
c
c
c  iread: when iread = 1 (or for first run) and icmobs .lt. 2
c     read new theoretical sets of data from units idsin
c     (the latter only when icmobs = 0)
      iread=1
      idsin=11
c
c  irew: if irew = 1 rewind theoretical mode data sets before reading
      irew=1
c  lwin1, lwin2, nwin1, nwin2: windows in l and order when reading
c  theoretical modes. window in l (order) only used if lwin2 .ge. lwin1
c  (nwin2 .ge. nwin1).
      lwin1=0
      lwin2=-1
      nwin1=0
      nwin2=-1
c  frqwn1, frqwn2: window in frequency, in micorhz, defined as
c     for windows in degree and order.
      frqwn1=0
      frqwn2=-1
c  icdata, epsdt: when icdata .gt. 0 separate modes for different
c     models. 
c     models are assumed to be different 
c     if any of the first icdata elements in the array data 
c     (which is part of the summary) has a relative difference 
c     greater than epsdt.
      icdata=0
      epsdt=1.d-4
c  nmod: when nmod .gt. 0  take differences in specified model.
      nmod=0
c  itype(1 - 2): itype(i) determines type of frequency from data set
c  no i.
c  itype = 1: eigenfrequency, possibly corrected for perturbation
c     in gravitation potential by perturbation analysis
c  itype = 2: uncorrected eigenfrequency (i.e. would give proper
c     eigenfrequency in the cowling approximation)
c  itype = 3: variational frequency
c  itype = 4: frequency from Richardson extrapolation
      itype(1)=3
      itype(2)=3
c  ipfrq = 1: print results in terms of cyclic frequencies in muhz
c  ipfrq = 2: print results in terms of dimensionless frequencies
c  Note: ipfrq can currently only be used with itype = 1 and 2.
      ipfrq=1
c  idiff = 1: relative differences
c  idiff = 2: absolute differences
      idiff=2
c  iscdif: if iscdif = 1 scale differences by energy ratio at
c     fixed frequency, normalized with surface vertical displacement.
c     Energy ratio is between energy of mode and energy at same
c     frequency for a mode of degree lnorm.
c     if iscdif = 2, use normalization with total surface displacement.
c     (note that this only works for icsin = 1, i.e. grand summaries).
c     if iscdif = 3 or 4 normalize with energy, divided by energy
c     at degree lnorm and frequency frnorm. For iscdif = 3 normalize
c     with vertical surface displacement, for iscdif = 4 with total
c     surface displacement.
      iscdif=0
c  initsc: if initsc = 1 reinitialize energy as function of
c    frequency for first degree in input set of modes.
      initsc=0
c  xnorm: if xnorm .gt. 0. and xnorm .lt. 2, rescale inertia ratio
c     to xnorm using isothermal atmosphere approximate eigenfunctions.
      xnorm=0
c  xftc: assumed distance from surface to x = 1, in units of pressure
c     scale height.
      xftc=4.4833d0
c  frftc: assumed acoustical cut-off frequency in atmosphere.
c     if frftc .le. 0, reset xftc and frftc from value of V,g at
c     surface, assuming Gamma,1 = 5/3
      frftc=5255
c  lnorm, frnorm: normalization degree and (for iscdif = 3 or 4)
c     frequency for mode energy.
      lnorm=0
      frnorm=3000
c  iprdif: when iprdif = 1 print differences
      iprdif=1
c
c  ......................................................................
c
c  zero data for first model
c
      call zero(data,16)
c
c  set dataset names and numbers
c
      write(6,*) ' Files needed:'
      write(6,*) ' Mode input on d/s idsin (default 11)'
      write(6,*) ' Output of results on d/s 20'
      call ofiles
      call openf(20,'u','f')
c
10000 continue
c
c  reading input data
c
      write(6,*) 'iread,idsin,irew?'
      write(6,*) iread,idsin,irew
      read(5,*,end=85,err=85) iread,idsin,irew
      write(6,*) 'lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2?'
      write(6,*) lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2
      read(5,*) lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2
      write(6,*) 'icdata,epsdt,nmod?'
      write(6,*) icdata,epsdt,nmod
      read(5,*) icdata,epsdt,nmod
      write(6,*) 'itype(1-2),ipfrq,idiff,iscdif,initsc?'
      write(6,*) itype,ipfrq,idiff,iscdif,initsc
      read(5,*) itype,ipfrq,idiff,iscdif,initsc
      write(6,*) 'xnorm, xftc, frftc, lnorm, frnorm?'
      write(6,*) xnorm, xftc, frftc, lnorm, frnorm
      read(5,*) xnorm, xftc, frftc, lnorm, frnorm
      write(6,*) 'iprdif?'
      write(6,*) iprdif
      read(5,*) iprdif
c
      do 10010 i=1,2
      if(itype(i).ge.3) ipfrq=1
c
10010 continue
c
      icdata=icdata
      icdata=min0(icdata,6)
c
      write(6,100)
c
      write(6,*) 'lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2'
      write(6,*) lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2
      write(6,*) 'icdata,epsdt,nmod'
      write(6,*) icdata,epsdt,nmod
      write(6,*) 'itype(1-2),ipfrq,idiff,iscdif,initsc'
      write(6,*) itype,ipfrq,idiff,iscdif,initsc
      write(6,*) 'xnorm, xftc, frftc, lnorm, frnorm'
      write(6,*) xnorm, xftc, frftc, lnorm, frnorm,tail
      write(6,*) 'iprdif'
      write(6,*) iprdif
c
      nnnn=0
c
      mm=0
c
c  test for reading new data
c
      if(iread.ne.1.and.ifsrun.eq.0) go to 15500
c
      ifsrun=0
c
c  start reading data set
c
      ifrd=0
c
c  test for opening dataset
c
      if(idsin.ne.idsinp) then
        call openf(idsin,'o','u')
        idsinp=idsin
        ifrd=1
      end if
c
      if(irew.eq.1) then
        rewind idsin
        ifrd=1
      end if
c
c  initialize mode and model counters for this data set
c
      n=nnnn
      m=mm
c
c  diagnostic output
c
      nfrst=1
c
c  test for mode array full
c
    5 if(n.ge.nmdmax) then
        write(6,102) idsin,nmdmax
        go to 10
      end if
c
c  read mode
c
      read(idsin,end=10,err=10) cs
c
      lrd=nint(cs(18))
      nrd=nint(cs(19))
      frrq=16666.666667d0/cs(25)
c
c  test for windowing in l and order
c
      if(lwin2.lt.lwin1) go to 11010
c
c  use economical windowing in l when only one
c  model is considered
c
      if(icdata.eq.0) then
c
        if(lrd.lt.lwin1) go to 5
        if(lrd.gt.lwin2) go to 10
        go to 11010
      end if
c
      if(lrd.lt.lwin1.or.lrd.gt.lwin2) go to 5
c
11010 if(nwin2.ge.nwin1.and.(nrd.lt.nwin1.or.nrd.gt.nwin2)) go to 5
c
c  window in frequency
c
      if(frqwn2.ge.frqwn1.and.(ffrq.lt.frqwn1.or.ffrq.gt.frqwn2))
     *  go to 5
c
c  store mode
c
      n=n+1
      l(n)=lrd
      nord(n)=nrd
c
c  set frequencies, depending on itype
c
      do 11020 i=1,2
      call stfgms(itype(i),cs,ffrq,sigrd)
c
c  test for case
c
      if(ipfrq.eq.1) then
c
c  set frequency in muhz
c
	if(ffrq.eq.0.) ffrq=16666.667d0/cs(25)
c
      else
c
c  dimensionless frequency
c
         ffrq=sqrt(sigrd)
c
      end if
c
11020 frq(i,n)=ffrq
c
c  test for setting energy. normalization depends on value of iscdif
c
      if(iscdif.eq.1.or.lrd.eq.0) then
        ekin(n)=cs(24)
      else if(iscdif.eq.2) then
        yrat=cs(31)/cs(30)
        ekin(n)=cs(24)/(1+yrat**2/(lrd*(lrd+1)))
      end if
c
c  test for setting model parameters
c
      if(iscdif.ge.1.and.xnorm.gt.0.and.xnorm.lt.2) then
        vgs=cs(10)
        if(xnorm.ne.1.or.frftc.le.0.or.xftc.le.0) then
          if(n.eq.1) then
            frqfct=1.d6*sqrt(6.6732d-8*cs(2)/cs(3)**3)/6.283185307d0
            xs=cs(23)
          else
            xs=max(xs,cs(23))
          end if
        end if
      end if
c
c  setting model index
c
      if(icdata.gt.0) then
c
c  test for new model
c
        do 6 i=1,icdata
        i1=i+1
        if(nfrst.eq.1) go to 6
        if(abs(cs(i1)/datap(i)-1).gt.epsdt) go to 7
    6   datap(i)=cs(i1)
        nfrst=0
c
      end if
c
      if(m.gt.0) go to 5
c
c  for icdata = 0, set model data into data(.,1), unless this is
c  continuation of reading
c
    7 m=m+1
      nm(m)=n
      do 8 i=1,8
      datap(i)=cs(i+1)
    8 data(i,m)=cs(i+1)
      go to 5
c
c  end of file found on this dataset. set total number of modes
c
   10 nnnn=n
c
c  to simplify calculation of number of modes for last model,
c  set 1+nnnn into nmd
c
      nm(m+1)=n+1
c
      write(6,110) idsin,m,nnnn
c
c
c  total number of models in this dataset
c
   15 mm=m
c
c  test for setting scaled energy
c
      if(iscdif.ge.1) then
c
c  test for renormalizing energy
c
        if(xnorm.gt.0.and.xnorm.lt.2) then
          call rnekin(l,nord,frq,ekin,nnnn,xnorm,xftc,frftc,xs,vgs,
     *      frqfct)
        end if
        call scekin(l,nord,frq,ekin,nnnn,initsc,iscdif,lnorm,frnorm)
      end if
c
c  output file name to result file
c
15500 write(20,120)
      call stfile(idsin,nfin)
      write(20,122) file(nfin)
c
c  test for comparing specific models
c
      if(icdata.ne.0) then
c
c  print model indices
c
        mmm=mm
        write(20,132) (m,(data(i,m),i=1,8),nm(m),
     *    m=1,mmm)
        write(6,130) (m,(data(i,m),i=1,8),nm(m),
     *    m=1,mmm)
c
      end if
c
c  ----------------------------------------------------------
c
c  analyze model no. m1 on dataset
c  when model indices have been set up, various options are possible
c  for the determination of m1 and m2
c
c  if icdata = 0 use m1 = m2 = 1
c
 1620 if(icdata.gt.0) go to 1630
c
      m1=1
      go to 1750
c
c  test for models specified in nmod1(i) and nmod2(i)
c
 1630 imd=0
      m1=nmod
      if(m1.le.0.or.m1.gt.mm) go to 90
      imd=1
c
 1750 write(6,150) m1
      write(6,152) (data(i,m1),i=1,4)
      write(20,153) (data(i,m1),i=1,4)
c
      n1=nm(m1)
      nn1=nm(m1+1)-n1
c
c  now calculate and output frequency differences
c
   19 write(20,156) (i, cnutyp(itype(i)),i=1,2)
c
      call frqdf1(l,nord,frq,n1,nn1,ipfrq,idiff,iscdif,
     *  lnorm,frnorm,ekin,xnorm,xftc,frftc,xs,vgs,iprdif)
c
   80 continue
      go to 10000
c
   85 stop
c
c  specified model number not available
c
   90 write(6,160) m1,mm1
c
      stop 1
  100 format(///)
  102 format(///' ***** mode array full for unit',i4,
     *  '  number of modes is',i5)
  105 format(a)
  110 format(//' end reading on data set',i4/
     *  ' total number of models:',i3/
     *  ' total number of modes:',i6)
  120 format('# Data files:'/'#')
  122 format('#   ',a60)
  130 format(///' models on file'//' m, data(1-8), n start:'//
     *  (i4,1p8e13.5,i5))
  132 format(/'#'/'# models on file'/'#'/
     *  '# m, data(1-8), n start:'/'#'/
     *  ('#',i4,1p8e13.5,i5))
  150 format(///' analyze frequencies for model',i3,
     *  ' on file '/)
  152 format(' data:',1p4e13.5)
  153 format('# data:',1p4e13.5)
  156 format('#'/'#'/('# data type ',i3,' : ',a))
  160 format(///1x,10(1h*),' model m1 =',i4,
     *  '  not available. mm1 =',i4)
      end
      subroutine frqdf1(l,nord,frq,ns1,nn1,ipfrq,idiff,iscdif,
     *  lnorm,frnorm,ekin,xnorm,xftc,frftc,xs,vgs,iprdif)
c
c  prints and plots differences between frequencies in frq(1,n)
c  and frq(2,n).
c
c  the differencing, plotting and printing are controlled by the
c  remaining parameters, which are described in the calling
c  programme.
c
c  kmax returns the number of curves stored, but not plotted.
c
      implicit double precision (a-h, o-z)
      dimension itype(2),l(1),nord(1),frq(2,1),ekin(1)
c
c  set initial mode index
c
      if(iscdif.eq.1) then
        write(6,105) lnorm
        write(20,106) lnorm
      else if(iscdif.eq.2) then
        write(6,107) lnorm
        write(20,108) lnorm
      else if(iscdif.eq.3) then
        write(6,109) lnorm,frnorm
        write(20,110) lnorm,frnorm
      else if(iscdif.eq.4) then
        write(6,111) lnorm,frnorm
        write(20,112) lnorm,frnorm
      end if
c
      if(iscdif.ge.1.and.xnorm.gt.0.and.xnorm.lt.2) then
        write(6,115) xnorm,xftc,frftc,xs,vgs
        write(20,116) xnorm,xftc,frftc,xs,vgs
      end if
c
c  test for heading
c
      ihead=ipfrq+2*(idiff-1)
      if(iprdif.ne.1) go to 20
      go to (1911,1912,1913,1914), ihead
 1911 write(6,121)
      write(20,131)
      go to 20
 1912 write(6,122)
      write(20,132)
      go to 20
 1913 write(6,123)
      write(20,133)
      go to 20
 1914 write(6,124)
      write(20,134)
c
   20 continue
c
c  step in modes
c
      n1=ns1-1
c
      do 70 n=1,nn1
      n1=n1+1
   60 go to (61,62), idiff
c  relative difference
   61 dfrq=frq(2,n1)/frq(1,n1)-1
      go to 65
c  absolute difference
   62 dfrq=frq(2,n1)-frq(1,n1)
c
c  test for energy scaling of difference
c
   65 if(iscdif.ge.1) dfrq=ekin(n1)*dfrq
c
      if(iprdif.eq.1)
     *  write(6,140) n1,l(n1),nord(n1),frq(1,n1),frq(2,n1),dfrq
      write(20,142) l(n1),nord(n1),
     *  frq(1,n1),dfrq,frq(1,n1)/(l(n1)+0.5d0)
c
   70 continue
c
c  end comparison
c
   80 continue
      return
c
  102 format(///1x,a)
  105 format(//' frequency differences scaled by kinetic energy',
     *  ' ratio'/
     *  ' to l =',i5,' at fixed frequency'/
     *  ' normalized with vertical displacement')
  106 format('#'/'# frequency differences scaled by kinetic energy',
     *  ' ratio'/
     *  '# to l =',i5,' at fixed frequency'/
     *  '# normalized with vertical displacement')
  107 format(//' frequency differences scaled by kinetic energy',
     *  ' ratio'/
     *  ' to l =',i5,' at fixed frequency'/
     *  ' normalized with total displacement')
  108 format('#'/'# frequency differences scaled by kinetic energy',
     *  ' ratio'/
     *  '# to l =',i5,' at fixed frequency'/
     *  '# normalized with total displacement')
  109 format(//' frequency differences scaled by kinetic energy',
     *  ' ratio'/
     *  ' to l =',i5,' at frequency =',f12.3,' microHz'/
     *  ' normalized with vertical displacement')
  110 format('#'/'# frequency differences scaled by kinetic energy',
     *  ' ratio'/
     *  '# to l =',i5,' at frequency =',f12.3,' microHz'/
     *  '# normalized with vertical displacement')
  111 format(//' frequency differences scaled by kinetic energy',
     *  ' ratio'/
     *  ' to l =',i5,' at frequency =',f12.3,' microHz'/
     *  ' normalized with total displacement')
  112 format('#'/'# frequency differences scaled by kinetic energy',
     *  ' ratio'/
     *  '# to l =',i5,' at frequency =',f12.3,' microHz'/
     *  '# normalized with total displacement')
  115 format(/' Scaling renormalized to x =',f10.6,
     *  ' using isothermal atmosphere expression'/
     *  ' parameters: xftc =',f10.5,'  frftc =',f10.2,'  microHz',
     *  '   xs =',f10.6/' V,gs =',f10.2)
  116 format('#'/'# Scaling renormalized to x =',f10.6,
     *  ' using isothermal atmosphere expression'/
     *  '# parameters: xftc =',f10.5,'  frftc =',f10.2,'  microHz',
     *  '   xs =',f10.6/'# V,gs =',f10.2)
  121 format(///' n1,l1,ord1,cyclic frq1,frq2 (muhz),',
     *  'relative difference (frq2 - frq1)'/)
  122 format(///' n1,l1,ord1,dimensionless frq1,frq2 ,',
     *  'relative difference (frq2 - frq1)'/)
  123 format(///' n1,l1,ord1,cyclic frq1,frq2 (muhz),',
     *  'absolute difference (frq2 - frq1) (muhz)'/)
  124 format(///' n1,l1,ord1,dimensionless frq1,frq2 ,',
     *  'absolute difference (frq2 - frq1)'/)
  131 format('#'/'#'/'# l1,ord1,cyclic frq1 (muhz),',
     *  'rel. diff. (frq2 - frq1), ',
     *  '(cyclic frq1)/(l+1/2)'/'#')
  132 format('#'/'#'/'# l1,ord1,dimensionless frq1 ,',
     *  'rel. diff. (frq2 - frq1), ',
     *  '(dimensionless frq1)/(l+1/2)'/'#')
  133 format('#'/'#'/'# l1,ord1,cyclic frq1 (muhz),',
     *  'abs. diff. (frq2 - frq1) (muhz), ',
     *  '(cyclic frq1)/(l+1/2)'/'#')
  134 format('#'/'#'/'# l1,ord1,dimensionless frq1 ,',
     *  'abs. diff. (frq2 - frq1)',
     *  '(dimensionless frq1)/(l+1/2)'/'#')
  140 format(3i5,0p2f12.3,1pe13.5)
  142 format(2i5,0pf12.3,1p2e13.5)
      end
      subroutine scekin(l,nord,frq,ekin,nn,init,iscdif,lnorm,frnorm)
c
c  resets ekin to ratio of energies at fixed frequency
c
c  if init = 1, or in first call, ratio is taken to modes
c  with degree lnorm or, if lnorm .le. 0, 
c  with the first degree in the array l(.).
c  if init .ne. 1, ratio is taken to energy assumed to be
c  already in e1(.), as a function of frequency in fr1(.).
c
      implicit double precision (a-h, o-z)
      dimension l(1),nord(1),frq(2,1),ekin(1),fr1(100),
     *  e1(100),eint(1)
c
      data nne /0/
c
c  test for initializing
c
      if(init.ne.1.and.nne.gt.0) go to 20
c
c  set e1 to energy at lnorm, or first l-value
c
      if(lnorm.gt.0) then
        n1=0
        do 5 n=1,nn
        if(l(n).eq.lnorm) then
          n1=n
          go to 8
        end if
    5   continue
c
    8   if(n1.eq.0) then
          write(6,100) lnorm,l(1)
          lnorm=l(1)
          n1=1
        end if
c
      else
        n1=1
        if(l(1).ne.0) then
          write(6,100) lnorm,l(1)
          lnorm=l(1)
        end if
c
      end if
c
      lp=l(n1)
c
      nne=0
      ns=0
      do 10 n=n1,nn
      if(l(n).ne.lp) go to 15
      ns=ns+1
      fr1(ns)=frq(1,n)
   10 e1(ns)=log(ekin(n))
c
   15 nne=ns
c
   20 n1=1
c
c  interpolate in fr1 and set ratio
c
c  test for interpolating to fixed frequency
c
   30 if(iscdif.ge.3) then
        call lir1(frnorm,fr1,eint,e1,1,1,nne,1,inter)
      else
        frp=1.d10
      end if
c
      write(6,30091) (n, fr1(n), e1(n), n=1,nne)
30091 format(//' n, frnorm, enorm:'/(i5,f10.3,f10.5))
c
      do 40 n=1,nn
      if(iscdif.le.2) then
        int=2
        if(frq(1,n).lt.frp) int=1
c..        write(6,*) 'frp, frq(1,n), int:', frp, frq(1,n), int
        frp=frq(1,n)
        call lir(frp,fr1,eint,e1,1,1,nne,int,inter)
c..        write(6,*) eint(1)
      end if
   40 ekin(n)=ekin(n)/exp(eint(1))
c
c  diagnostic output
c
      nnw=min0(nn,300)
      if(iscdif.eq.1) then
        write(6,105)
      else 
        write(6,110)
      end if
c
      write(6,120) (n,l(n),nord(n),frq(1,n),ekin(n),n=1,nnw)
      return
  100 format(//'  ***** error in freqdif. lnorm = ',i5,
     * ' is not in mode set'/
     *         '        lnorm reset to',i5)
  105 format(//' scaling of energy.',
     *  ' normalize with surface vertical displacement.'/)
  110 format(//' scaling of energy.',
     *  ' normalize with surface total displacement.'/)
  120 format(' n, l, order, frequency, scaled energy:'/(3i5,0pf10.2,
     *  1pe13.5))
      end
      subroutine stfgms(itype,cs,frq,sig)
c
c  sets frequency data from grand summary in cs. case depends on itype.
c
c  itype = 1: eigenfrequency, possibly corrected for perturbation
c     in gravitation potential by perturbation analysis
c  itype = 2: uncorrected eigenfrequency (i.e. would give proper
c     eigenfrequency in the cowling approximation)
c  itype = 3: variational frequency
c  itype = 4: frequency from Richardson extrapolation
c
c  returns proper dimensionless frequency in sig, and
c  cyclic frequency (in microHz) frq 
c
c  original version: 15/2/89
c
c                     ................................
c
      implicit double precision (a-h, o-z)
      dimension cs(*)
      data icerr /0/
c
c  test for type
c
      if(itype.eq.1) then
c
        sig=cs(21)
        if(sig.le.0) sig=cs(20)
        frq=sqrt(sig/cs(20))*16666.66667d0/cs(25)
        frq1=1.5915494d5*sqrt(6.6732d-8*cs(2)*sig/cs(3)**3)
        if(abs(frq/frq1-1).gt.1.d-5.and.icerr.le.10) then
          icerr=icerr+1
          write(6,110) frq, frq1
        end if
c
      else if(itype.eq.2) then
c
c  uncorrected frequency
c  include a test of frequency factor, for the time being,
c  for cases where cgrav .ne. 6.6732e-8
c
        sig=cs(20)
        frq=16666.666667d0/cs(25)
c
      else if(itype.eq.3) then
c
        sig=cs(21)
        if(sig.le.0) sig=cs(20)
        if(cs(27).gt.0) then
          frq=1000*cs(27)
        else
          frq=16666.666667d0/cs(25)
        end if
c
      else if(itype.eq.4) then
c
c  frequency from cs(37) (Richardson extrapolation)
c  if not available, use variational frequency or eigenfrequency
c
        sig=cs(21)
        if(sig.le.0) sig=cs(20)
        if(cs(37).ne.0) then
          frq=1000*cs(37)
        else
          icerr=icerr+1
          if(icerr.le.20) write(6,120)
          if(cs(27).gt.0) then
            frq=1000*cs(27)
          else
            frq=16666.666667d0/cs(25)
          end if
        end if
      end if
c
      return
  110 format(' *** warning. in s/r rdfreq, with icase = 4, frq =',
     *  1pe13.6,' .ne. frq1 =',e13.6)
  120 format(' *** error in s/r rdfreq with icase = 5. cs(37) not set')
      end
      subroutine rnekin(l,nord,frq,ekin,nn,xnorm,xftc,frftc,xs,vgs,
     *  frqfct)
c
c  rescales mode inertia using approximate expression for isothermal
c  atmosphere.
c
c  original version: 18/7/89
c
      implicit double precision (a-h, o-z)
      dimension l(1),nord(1),frq(2,1),ekin(1)
c
      gm1=1.66666666667d0
c
c  test for setting xftc and frftc
c
	write(6,*) 'xnorm,xftc,frftc,xs,vgs',xnorm,xftc,frftc,xs,vgs
      if(frftc.le.0.or.xftc.le.0) then
        xftc=(xs-1)*vgs*gm1
        frftc=frqfct*gm1*sqrt(vgs)/2
        write(6,110) xftc, frftc
      end if
c
      v=vgs*gm1
c
c  set normalized distance from surface
c
      if(xnorm.ne.1) then
        xftcr=xftc*(xs-xnorm)/(xs-1)
      else
        xftcr=xftc
      end if
c
c  step through modes
c
      do 30 n=1,nn
      frqrat=frq(1,n)/frftc
c
c  expression including term in degree
c
      xx=frqrat*frqrat
      xx=1-xx+4*l(n)*(l(n)+1)*(1-4*(1-1/gm1)/(gm1*xx))/(v*v)
      if(xx.lt.0) then
        arat=1
        write(6,120) l(n), nord(n), frq(1,n)
      else
        xx=1-sqrt(xx)
        arat=exp(-xftcr*xx)
      end if
c
   30 ekin(n)=ekin(n)/arat
      return
  110 format(/' xftc, frftc reset from model quantities, to',
     *  f10.5,f10.2)
  120 format(' ***** error. isothermal ratio undefined for l =',i5,
     *  '  order =',i5,' frequency =',f10.2,' microHz')
      end
