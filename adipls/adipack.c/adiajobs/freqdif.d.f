      program main
c
c  differences between adiabatic frequencies. assumed sorted,
c  increasingly, in the same way determined by isort
c
c  frequencies may be given in the form of either grand summaries,
c  short summaries or observed data, as selected by icsin(1,2)
c  being 1, 2 or 3 respectively.
c
c  notes on dimensioning: the arrays l, nord and freq are
c  dimensioned as (2,nmdmax), where nmdmax is the maximum
c  number of modes that may be read in. in addition
c
c  modified 19/7/1985 to allow comparison of two sets of observational
c  frequencies.
c
c  modified 30/9/86, to take differences between two arbitrary sets
c  of data.
c
c  modified 18/7/89, to allow approximate resetting of inertia scaling
c  to a normalization point other than the surface.
c
c  modified 15/5/90, to allow scaling with interia normalized to a
c  single (l, nu) point.
c
c  modified 11/6/90, to output also relative differences of mode 
c  inertias on file ttt.dekin.out
c
c  modified 18/2/92, to allow calculation of relative difference
c  of dimensionless frequency for any itype
c
c  Modified 13/12/95, to allow input and output of errors on
c  observational frequencies
c
c  Modified 27/3/96, adding first line indicating type
c
c  Modified 3/5/97, adding option of scaling frequencies by
c  radius factor, to correct for difference in radii.
c
c  Modified 12/5/10, taking out reference to option icmobs = 1
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      parameter (nmdmax=10000)
      character*280 form,file
      character tail*5, cnutyp*40
      dimension form(7),
     *  cs(50),idsin(2),icsin(2),itype(2),idsinp(2),
     *  l(2,nmdmax),nord(2,nmdmax),frq(2,nmdmax),sig(2,nmdmax),nnnn(2),
     *  data(2,8,20),nm(2,20),datap(8),mm(2),nmod1(5),nmod2(5),
     *  errobs(10),lval(400),
     *  ekin(2,nmdmax),ekins(nmdmax),error(2,nmdmax),cnutyp(7)
      common/cofile/ nfiles, idsfil(20), file(20), iopen(20)
c
c..      namelist/exec/ iread,idsin,icsin,icntrd,irew,
c..     *  icmobs,idsobs,idiagr,
c..     *  lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2,
c..     *  isort,icdata,epsdt,imod1,nmod1,nmod2,
c..     *  itype,ipfrq,idiff,invdif,iscdif,initsc,
c..     *  xnorm,xftc,frftc, lnorm, frnorm, radfct
c..     *  iprdif
      data ifsrun /1/
      data kmaxpl /0/
      data idsinp /-1, -1/
      data form /'u','u','f','u','u','f','f'/
      data tail /'    @'/
      data cnutyp /
     *  ' (possibly corrected) eigenfrequency',
     *  ' (uncorrected) eigenfrequency',
     *  ' variational frequency',
     *  ' Richardson extrapolation frequency',
     *  ' observed frequency',
     *  ' observed frequency, with errors',
     *  ' observed frequency, with inertia'/
c
c  Gravitational constant
c
      data cgrav /6.67232d-8/
c
      twopi=8.d0*atan(1.d0)
c
c
c  iread: when iread = 1 (or for first run) and icmobs .lt. 2
c     read new theoretical sets of data from units idsin(1) and idsin(2)
c     (the latter only when icmobs = 0)
      iread=1
      idsin(1)=11
      idsin(2)=12
c  icsin(1-2): determine type of dataset read on idsin(1-2).
c  icsin(.) = 1: grand summary
c  icsin(.) = 2: short summary
c  icsin(.) = 3: observed data
c  icsin(.) = 4: grand summary, set frequency from eigenfrequency
c     in cs(20) (to allow using Cowling approximation)
c  icsin(.) = 5: grand summary, set frequency from Richardson 
c     extrapolation in cs(37) 
c  icsin(.) = 6: observed data, with errors as 4th column
c     In this case, the output differences are given with errors
c  icsin(.) = 7: `observed data', with mode inertia as 4th column
c  for icsin(.) .lt. 0, read data from single-precision dataset,
c     according to value of abs(icsin(.)), as above.
      icsin(1)=1
      icsin(2)=1
c  icntrd: if icntrd = 1 (2) go back to read new /exec/ and then
c     read additional frequencies to first (both first and second)
c     set of modes before taking differences.
c     note: this has not been consistently implemented to take several
c     models into account.
      icntrd=0
c
c  irew: if irew = 1 rewind theoretical mode data sets before reading
      irew=1
c  icmobs: redundant, but left for consistency with older input files
      icmobs=0
      idsobs=5
c  idiagr: if idiagr = 1 print diagnostics from rdofrq
      idiagr=0
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
c  isort = 1: sorted after l first
c  isort = 2: sorted after n first
      isort=1
c  icdata, epsdt: when icdata .gt. 0 separate modes for different
c     models. 
c     when reading grand summary, models are assumed to be different 
c     if any of the first icdata elements in the array data 
c     (which is part of the summary) has a relative difference 
c     greater than epsdt.
c     when reading short summary, a new model is flagged by a model
c     record in the input data.
c     model separation cannot be made with observed data.
      icdata=0
      epsdt=1.e-4
c  nmod1, nmod2: when nmod1(1) and nmod2(1) .gt. 0 
c  take differences between specified models.
c  for imod1 .ne. 1 take differences between frequencies for
c  models nmod1(i) and nmod2(i), i = 1,..., for as long as
c  nmod1(i) and nmod2(i) .gt. 0.
c  for imod1   =  1 take differences between frequencies for
c  model nmod1(1) and models nmod2(i), i = 1,..., for as long as
c  nmod2(i) .gt. 0.
      imod1=0
      call izero(nmod1,5)
      call izero(nmod2,5)
c  itype(1 - 2): itype(i) determines type of frequency from data set
c  no i.
c  itype = 1: eigenfrequency, possibly corrected for perturbation
c     in gravitation potential by perturbation analysis
c  itype = 2: uncorrected eigenfrequency (i.e. would give proper
c     eigenfrequency in the cowling approximation)
c  itype = 3: variational frequency
c  itype = 4: frequency from Richardson extrapolation
c  Note: itype = 2 and 4 can only be used with icase = 1
c     itype =1 - 4 is meaningless for icase = 3. However itype = 5 and 6
c     are used internally to flag for observed frequencies.
      itype(1)=3
      itype(2)=3
c  ipfrq = 1: print results in terms of cyclic frequencies in muhz
c  ipfrq = 2: print results in terms of dimensionless frequencies
c  Note: ipfrq can currently only be used with itype = 1 and 2.
      ipfrq=1
c  idiff = 1: relative differences
c  idiff = 2: absolute differences
c  idiff = 3: relative difference between dimensionless frequency,
c     as determined by itype. Note that this is useful for
c     comparisons of models of different radii.
      idiff=2
c  invdif: if invdif = 1 the difference is (first set) - (second set),
c     otherwise it is (second set) - (first set).
      invdif=0
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
      xftc=4.4833
c  frftc: assumed acoustical cut-off frequency in atmosphere.
c     if frftc .le. 0, reset xftc and frftc from value of V,g at
c     surface, assuming Gamma,1 = 5/3
      frftc=5255
c  lnorm, frnorm: normalization degree and (for iscdif = 3 or 4)
c     frequency for mode energy.
      lnorm=0
      frnorm=3000
c  radfct: If radfct .ne. 1, scale frequencies of set 1 by factor
c  radfct**(-1.5), to correct (homologously) for difference in
c  radii
      radfct=1.d0
c  iprdif: when iprdif = 1 print differences
      iprdif=1
c
c  ......................................................................
c
      kcntrd=0
c
c  zero data for first model
c
      call zero(data,16)
c
c  set dataset names and numbers
c
      write(6,*) ' Files needed:'
      write(6,*) ' Mode input on d/s idsin(1) and idsin(2)',
     *           ' (default 11 and 12)'
      write(6,*) ' Output of results on d/s 20'
      call ofiles
      call openf(20,'u','f')
c
c  open file for output of inertia differences
c
      open(30,file='ttt.dekin.out',status='unknown')
c
10000 continue
c..      read(5,exec,end=85,err=85)
c
c  reading input data
c
      write(6,*) 'iread,idsin(1-2),icsin(1-2),icntrd,irew?'
      write(6,*) iread,idsin,icsin,icntrd,irew
      read(5,*,end=85,err=85) iread,idsin,icsin,icntrd,irew
      write(6,*) 'icmobs,idsobs,idiagr?'
      write(6,*) icmobs,idsobs,idiagr
      read(5,*) icmobs,idsobs,idiagr
      write(6,*) 'lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2?'
      write(6,*) lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2
      read(5,*) lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2
      write(6,*) 'isort,icdata,epsdt,imod1,nmod1,nmod2?'
      write(6,*) isort,icdata,epsdt,imod1,nmod1,nmod2
      read(5,*) isort,icdata,epsdt,imod1,nmod1,nmod2
      write(6,*) 'itype(1-2),ipfrq,idiff,invdif,iscdif,initsc?'
      write(6,*) itype,ipfrq,idiff,invdif,iscdif,initsc
      read(5,*) itype,ipfrq,idiff,invdif,iscdif,initsc
      write(6,*) 'xnorm, xftc, frftc, lnorm, frnorm, radfct?'
      write(6,*) xnorm, xftc, frftc, lnorm, frnorm, radfct
      read(5,*) xnorm, xftc, frftc, lnorm, frnorm, radfct
      write(6,*) 'iprdif?'
      write(6,*) iprdif
      read(5,*) iprdif
c
c  test for consistent itype and icase
c
      ireset=0
c
      do 10010 i=1,2
      icsina=iabs(icsin(i))
      if(icsina.eq.1) then
        if(itype(i).eq.2) then
          icsina=4
          ireset=1
        else if(itype(i).eq.4) then
          icsina=5
          ireset=1
        end if
      else if(icsina.eq.2) then
        if(itype(i).eq.2) then
          itype(i)=1
          ireset=2
        else if(itype(i).eq.4) then
          itype(i)=3
          ireset=2
        end if
      else if(icsina.eq.3) then
        itype(i)=5
      else if(icsina.eq.6) then
        itype(i)=6
      else if(icsina.eq.7) then
        itype(i)=7
      end if
c
      if(ireset.eq.1) then
        icsin(i)=isign(icsina,icsin(i))
        write(6,102) i,i, icsin(i), itype(i)
      else if(ireset.eq.2) then
        write(6,103) i,i, icsin(i), itype(i)
      end if
c
      if(itype(i).ge.3) ipfrq=1
c
10010 continue
c
      icdata=icdata
      icdata=min0(icdata,6)
c
c  only allow iscdif = 2 or 4, if first dataset contains grand summaries
c
      icsina=iabs(icsin(1))
      if(icsina.ne.1.and.icsina.lt.4) then
        if(iscdif.eq.2) then
          iscdif=1
          write(6,104) 
        else if(iscdif.eq.4) then
          iscdif=3
          write(6,105) 
        end if
c
c  only allow xnorm .ne. 1, or setting renormalization parameters from
c  data, if first dataset contains grand summary
c
        if(xnorm.ne.1.and.xnorm.gt.0.and.xnorm.lt.2) then
          xnorm=1
          write(6,106)
        end if
c
        if(frftc.le.0) then
          frftc=5255
          xftc=4.4833
          write(6,107) xftc, frftc
        end if
      end if
c
      write(6,100)
c..      write(6,exec)
c
      write(6,*) 'iread,idsin(1-2),icsin(1-2),icntrd,irew'
      write(6,*) iread,idsin,icsin,icntrd,irew,tail
      write(6,*) 'icmobs,idsobs,idiagr'
      write(6,*) icmobs,idsobs,idiagr,tail
      write(6,*) 'lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2'
      write(6,*) lwin1,lwin2,nwin1,nwin2,frqwn1,frqwn2,tail
      write(6,*) 'isort,icdata,epsdt,imod1,nmod1,nmod2'
      write(6,*) isort,icdata,epsdt,imod1,nmod1,nmod2,tail
      write(6,*) 'itype(1-2),ipfrq,idiff,invdif,iscdif,initsc'
      write(6,*) itype,ipfrq,idiff,invdif,iscdif,initsc,tail
      write(6,*) 'xnorm, xftc, frftc, lnorm, frnorm, radfct'
      write(6,*) xnorm, xftc, frftc, lnorm, frnorm, radfct,tail
      write(6,*) 'iprdif'
      write(6,*) iprdif,tail
c
c  initialize parameters for renormalization
c
      xs=0
      vgs=0
c
c  test for continuing read
c
      if(kcntrd.gt.0) go to 10050
c
      nnnn(1)=0
      nnnn(2)=0
c
      mm(1)=0
      mm(2)=0
c
      irderr=0
c
c  test for reading new data
c
      if(iread.ne.1.and.ifsrun.eq.0) go to 15500
c
      ifsrun=0
c
c  input theoretical data
c
      idsmax=2
c
c  set possible scale factor for frequencies of first set
c
      frqscl=radfct**(-1.5d0)
c
c  start loop over data sets
c
10050 do 15 ids=1,idsmax
      idsrd=idsin(ids)
      icsrd=icsin(ids)
      itprd=itype(ids)
      icsrda=abs(icsrd)
c
c  test for resetting flag for errors
c
      if(icsrda.eq.6.or.icsrda.eq.7) irderr=1
c
      ifrd=0
c
c  test for opening dataset
c
      if(idsrd.ne.idsinp(ids)) then
        call openf(idsrd,'o',form(icsrda))
        idsinp(ids)=idsrd
        ifrd=1
      end if
c
      if(irew.eq.1) then
        rewind idsrd
        ifrd=1
      end if
c
c  test for first record on observational dataset
c
c  ***  surely this is taken care of in rdfreq  ****
c
c..      if(ifrd.eq.1.and.icsrda.eq.3) then
c..c
c..c  skip possible initialization records
c..c
c..        read(idsrd,*) lobsin
c..        if(lobsin.lt.0) then
c..          read(idsrd,*) nobsin
c..        else
c..          backspace idsrd
c..        end if
c..      end if
c
c  initialize mode and model counters for this data set
c
      n=nnnn(ids)
      m=mm(ids)
c
      if((icsrda.eq.3.or.icsrda.eq.6.or.icsrda.eq.7).and.m.eq.0) then
        m=1
        nm(ids,m)=1
      end if
c
c  diagnostic output
c
      write(6,10099) ids,idsrd,n,m
10099 format(//' ids, idsrd, n, m =',4i5)
      nfrst=1
c
c  test for mode array full
c
    5 if(n.ge.nmdmax) then
        write(6,112) idsrd,nmdmax
        go to 10
      end if
c
c  read mode. For icsrd = 6 read observed frequencies with errors
c  For icsrd = 7 read observed frequencies with inertia
c
      call rdfrqe(icsrd,idsrd,cs,lrd,nrd,sigrd,ffrq,ekinrd,errrd,
     *  ierr)
      if(ierr.gt.0) go to 10
c
c  test for model record in short summary
c
      if(lrd.lt.0) then
        if(icsrda.eq.2.and.(icdata.gt.0.or.m.eq.0)) then
c
c  add new or first (for m = 0) model to model index
c
          m=m+1
          nm(ids,m)=n+1
c
c  note that the first mode for this model is the following to be read
c
          do 5050 i=1,4
 5050     data(ids,i,m)=cs(i+2)
c
          go to 5
        end if
      end if
c
c  test for windowing in l and order
c
      if(lwin2.lt.lwin1) go to 11010
c
c  use economical windowing in l when isort = 1 and only one
c  model is considered
c
      if(isort.eq.1.and.icdata.eq.0) then
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
      l(ids,n)=lrd
      nord(ids,n)=nrd
c
c  test for case
c
      if(ipfrq.eq.1) then
c
c  set frequency in muhz
c  for itprd = 1 note that cs(25) gives period corresponding
c  to uncorrected eigenfrequency.
c
        if(icsrda.eq.1.and.(ffrq.eq.0.or.itprd.eq.1)) 
     *    ffrq=16666.666667d0*sqrt(sigrd/cs(20))/cs(25)
c
      else
c
c  dimensionless frequency
c
         ffrq=sqrt(sigrd)
c
      end if
c
      if(ids.eq.1) then
        frq(ids,n)=ffrq*frqscl
      else
        frq(ids,n)=ffrq
      end if
c
c  store dimensionless frequency, based on frq
c
      if(icsrda.ne.3.and.icsrda.ne.6) then
	sig(ids,n)=twopi*1.d-6*sqrt(cs(3)**3/(cgrav*cs(2)))*ffrq
      else
	sig(ids,n)=-1.e36
      end if
c
c  test for setting errors
c
      if(icsrda.eq.6) then
        error(ids,n)=errrd
        go to 5
      end if
c
c  store basic energy
c  
      ekin(ids,n)=ekinrd
c
c  test for setting energy. normalization depends on value of iscdif
c
      if(ids.eq.1.and.iscdif.ge.1) then
      iscdfm=mod(iscdif,2)
        if(iscdfm.eq.1.or.lrd.eq.0) then
          ekins(n)=ekinrd
        else if(iscdfm.eq.0) then
          yrat=cs(31)/cs(30)
          ekins(n)=ekinrd/(1+yrat**2/(lrd*(lrd+1)))
        end if
c
c  test for setting model parameters
c
        if(xnorm.gt.0.and.xnorm.lt.2) then
          vgs=cs(10)
          if(xnorm.ne.1.or.frftc.le.0.or.xftc.le.0) then
            if(n.eq.1) then
              frqfct=1.e6*sqrt(cgrav*cs(2)/cs(3)**3)/twopi
              xs=cs(23)
            else
              xs=max(xs,cs(23))
            end if
          end if
        end if
c
      end if
c
c  setting model index for icsrd = 1
c
      if(icsrda.ne.1.and.icsrda.lt.4) go to 5
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
c  for icdata = 0, set model data into data(ids,.,1), unless this is
c  continuation of reading
c
      if(m.gt.0.or.kcntrd.ne.0) go to 5
c
    7 m=m+1
      nm(ids,m)=n
      do 8 i=1,8
      datap(i)=cs(i+1)
    8 data(ids,i,m)=cs(i+1)
      go to 5
c
c  end of file found on this dataset. set total number of modes
c
   10 nnnn(ids)=n
c
c  to simplify calculation of number of modes for last model,
c  set 1+nnnn into nmd
c
      nm(ids,m+1)=n+1
c
      write(6,110) idsrd,m,nnnn(ids)
c
c  total number of models in this dataset
c
   15 mm(ids)=m
c
c  test for continued read
c
      if(icntrd.le.0) go to 15100
c
c  set controls for continued read
c
      kcntrd=icntrd
      idsmax=min0(idsmax,icntrd)
      go to 10000
c
c  end of continued read
c
15100 kcntrd=0
c
c  test for setting scaled energy
c
      if(iscdif.ge.1) then
c
c  test for renormalizing energy
c
        if(xnorm.gt.0.and.xnorm.lt.2) then
          call rnekin(l,nord,frq,ekins,nnnn(1),xnorm,xftc,frftc,xs,vgs,
     *      frqfct)
        end if
        call scekin(l,nord,frq,ekins,nnnn(1),initsc,iscdif,lnorm,frnorm)
      end if
c
c  end section for setting up data
c
c  ******************************************************
c
c  output data type and file names to result file
c
15500 idttyp=idiff
      if(iscdif.ge.1) idttyp=idttyp+10
      if(irderr.eq.1) idttyp=idttyp+100
c
      write(20,115) idttyp
      write(20,120)
      write(30,120)
      do 15510 ids=1,2
      call stfile(idsin(ids),nfin)
      write(20,122) file(nfin)
15510 write(30,122) file(nfin)
c
c  test for comparing specific models
c
      if(icdata.ne.0.and.idsmax.gt.0) then
c
c  print model indices
c
        do 16 ids=1,idsmax
        mmm=mm(ids)
        write(20,132) ids,(m,(data(ids,i,m),i=1,8),nm(ids,m),
     *    m=1,mmm)
        write(30,132) ids,(m,(data(ids,i,m),i=1,8),nm(ids,m),
     *    m=1,mmm)
   16   write(6,130) ids,(m,(data(ids,i,m),i=1,8),nm(ids,m),
     *    m=1,mmm)
c
      end if
c
      if(radfct.ne.1) then
	write(6,133) radfct
	write(20,134) radfct
	write(30,134) radfct
      end if
c
c  number of models on datasets
c
      mm1=mm(1)
      mm2=mm(2)
c
c  ----------------------------------------------------------
c
c  compare theoretical frequencies.
c
c  compare model no. m1 on dataset 1 with model no. m2 on dataset 2.
c  when model indices have been set up, various options are possible
c  for the determination of m1 and m2
c
c  if icdata = 0 use m1 = m2 = 1
c
 1620 if(icdata.gt.0) go to 1630
c
      m1=1
      m2=1
      go to 1750
c
c  test for models specified in nmod1(i) and nmod2(i)
c
 1630 imd=0
      m1=nmod1(1)
      m2=nmod2(1)
      if(m1.le.0.or.m2.le.0) go to 1680
      if(m1.gt.mm1.or.m2.gt.mm2) go to 90
      imd=1
      go to 1750
c
c  entry point for continuing pass through nmod1 and nmod2
c
 1650 imd=imd+1
      m2=nmod2(imd)
c  test for fixed model 1
      if(imod1.ne.1) m1=nmod1(imd)
      if(m1.le.0.or.m2.le.0) go to 80
      if(m1.gt.mm1.or.m2.gt.mm2) go to 90
c  compare models m1 and m2
      go to 1750
c
c   search for matching models on dataset 1 and dataset 2
c
 1680 m1=0
c
c  entry point for continuing search for matching models
c
   17 m1=m1+1
      if(m1.gt.mm1) go to 80
c  look for matching model on file 2
      do 1710 ms2=1,mm2
      m2=ms2
      do 1705 i=1,icdata
      if(abs(data(2,i,m2)/data(1,i,m1)-1).gt.epsdt) go to 1710
 1705 continue
c  match found
      go to 1750
 1710 continue
c  no match. write diagnostics
      write(6,140) m1
      go to 17
c
c   --------------------------------------------------------------
c
c  now model pair has been determined. set ranges for mode comparison
c
 1750 write(6,150) m1,m2
      do 1755 im=1,2
      if(data(im,1,m1).ne.0) then
        write(6,152) im,(data(im,i,m1),i=1,4)
        write(20,153) im,(data(im,i,m1),i=1,4)
        write(30,153) im,(data(im,i,m1),i=1,4)
      else
        write(6,154) im
        write(20,155) im
      end if
 1755 continue
c
      n1=nm(1,m1)
      n2=nm(2,m2)
      nn1=nm(1,m1+1)-n1
      nn2=nm(2,m2+1)-n2
c
c  now calculate and output frequency differences
c
   19 write(20,156) (i, cnutyp(itype(i)),i=1,2)
c
      call frqdif(l,nord,frq,sig,ekin,error,icsin,n1,n2,nn1,nn2,ipfrq,
     *  idiff,iscdif,lnorm,frnorm,ekins,invdif,isort,xnorm,xftc,frftc,
     *  xs,vgs,iprdif)
c
c  test for continuation
c
   75 if(icdata.eq.0) go to 80
c
c  test for continuing stepping through models specified by
c  nmod1 and nmod2
c
      if(imd.gt.0) go to 1650
c
c  continue searching for matching models
c
      go to 17
c
   80 continue
      go to 10000
c
   85 continue
c
c  test for plotting remaining plot variables
c
      if(kmaxpl.le.0) stop
c
      call frqdif(l,nord,frq,sig,ekin,error,icsin,n1,n2,nn1,nn2,ipfrq,
     *  idiff,iscdif,lnorm,frnorm,ekins,invdif,isort,xnorm,xftc,frftc,
     *  xs,vgs,iprdif)
c
      stop
c
c  specified model numbers not available
c
   90 write(6,160) m1,m2,mm1,mm2
c
      stop 1
  100 format(///)
  102 format(//' **** warning. icsin(',i1,'), itype(',i1,') reset to',
     *  2i3)
  103 format(//' **** error. icsin(',i1,'), itype(',i1,') reset to',
     *  2i3)
  104 format(//' ********* iscdif = 2 only allowed for grand summary.',
     *  ' iscdif reset to 1')
  105 format(//' ********* iscdif = 4 only allowed for grand summary.',
     *  ' iscdif reset to 3')
  106 format(//
     *  ' ********* xnorm .ne. 1 only allowed for grand summary.',
     *  ' xnorm reset to 1')
  107 format(//
     *  ' ********* setting frftc and xftc from model variables',
     *  ' only allowed for grand summary.'/
     *  10x,' xftc and frftc reset to default values =',f10.5,f10.2)
  110 format(//' end reading on data set',i4/
     *  ' total number of models:',i3/
     *  ' total number of modes:',i6)
  112 format(///' ***** mode array full for unit',i4,
     *  '  number of modes is',i5)
  115 format('# Data type :',i10)
  120 format('# Data files:'/'#')
  122 format('#   ',a60)
  130 format(///' models on file',i3//' m, data(1-8), n start:'//
     *  (i4,1p8e13.5,i5))
  132 format(/'#'/'# models on file',i3/'#'/
     *  '# m, data(1-8), n start:'/'#'/
     *  ('#',i4,1p8e13.5,i5))
  133 format(/'  Frequencies in first set have been rescaled'/
     *  ' corresponding to radius ratio =',f12.7/)
  134 format('#'/'#  Frequencies in first set have been rescaled'/
     *  '# corresponding to radius ratio =',f12.7/'#')
  140 format(///1x,10(1h*),' model no ',i2,' on file 1 has no',
     *  ' counterpart on file 2')
  150 format(///' compare frequencies for model',i3,
     *  ' on file 1, model',i3,' on file 2'/)
  152 format(' data',i1,':',1p4e13.5)
  153 format('# data',i1,':',1p4e13.5)
  154 format(' file ',i1,' is in the form of observed data')
  155 format('# file ',i1,' is in the form of observed data')
  156 format('#'/'#'/('# data type on dataset no',i3,' : ',a))
  158 format(///' compare frequencies for model',i3,
     *  ' on file 1 with observed frequencies.'//
     *  ' data1:',1p8e13.5)
  159 format('#'/'# compare frequencies for model',i3,
     *  ' on file 1 with observed frequencies.'/'#'/
     *  '# data1:',1p8e13.5)
  160 format(///1x,10(1h*),' model pair m1, m2 =',2i4,
     *  '  not available. mm1, mm2 =',2i4)
      end
      subroutine frqdif(l,nord,frq,sig,ekin,error,icsin,ns1,ns2,nn1,nn2,
     *  ipfrq,idiff,iscdif,lnorm,frnorm,ekins,invdif,isort,xnorm,xftc,
     *  frftc,xs,vgs,iprdif)
c
c  locates matching modes between the sets
c  (l(1,n),nord(1,n),frq(1,n), n=ns1, ..., ns1-1+nn1) and
c  (l(2,n),nord(2,n),frq(2,n), n=ns2, ..., ns2-1+nn2)
c  and prints and plots differences.
c
c  if isubmn =1 also determines mean difference curve from modes
c  with lsubmn .le. l. lsubmx and subtracts it from the differences.
c
c  the differencing, plotting and printing are controlled by the
c  remaining parameters, which are described in the calling
c  programme.
c
c  kmax returns the number of curves stored, but not plotted.
c
c  If icsin(1) and/or icsin(2) are 6, assumes errors set in error,
c  and produces error on frequency difference
c
c  Modified 14/12/95, to include error input and output
c
      implicit double precision (a-h, o-z)
      character*7 cerror
      dimension itype(2),l(2,1),nord(2,1),frq(2,1),sig(2,1),
     *  ekin(2,1),error(2,1), ekins(1), kl(1000),lk(200),lval(400),
     *  icsin(2)
c
c  set initial mode indices
c
      n1=ns1
      n2=ns2
      nf1=ns1+nn1-1
      nf2=ns2+nn2-1
c
c  set control for assumed ordering of modes
c
      isort1=isort-1
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
      if(iscdif.ge.1.and.xnorm.gt.0.and.xnorm.lt.2) then
        write(6,115) xnorm,xftc,frftc,xs,vgs
        write(20,116) xnorm,xftc,frftc,xs,vgs
      end if
      if(invdif.eq.1) then
        write(6,117)
        write(20,118)
      end if
c
c  initial value of k
c
      kinit=kmax+1
c
c  test for heading
c
      ihead=ipfrq+2*(idiff-1)
      if(icsin(1).eq.6.or.icsin(2).eq.6) then
        cerror=', error'
        lcerrr=7
        isterr=1
      else
        cerror=' '
        lcerrr=1
        isterr=0
      end if
      if(ihead.eq.1) then
        write(6,121) cerror(1:lcerrr)
        write(20,131) cerror(1:lcerrr)
      else if(ihead.eq.2) then
        write(6,122) cerror(1:lcerrr)
        write(20,132) cerror(1:lcerrr)
      else if(ihead.eq.3) then
        write(6,123) cerror(1:lcerrr)
        write(20,133) cerror(1:lcerrr)
      else if(ihead.eq.4) then
        write(6,124) cerror(1:lcerrr)
        write(20,134) cerror(1:lcerrr)
      else if(ihead.eq.5) then
        write(6,125) cerror(1:lcerrr)
        write(20,135) cerror(1:lcerrr)
      else if(ihead.eq.6) then
        write(6,126) cerror(1:lcerrr)
        write(20,136) cerror(1:lcerrr)
      end if
c
      write(30,140)
c
   20 continue
      if(isort1) 30,30,40
c
   30 if(n1.gt.nf1.or.n2.gt.nf2) go to 80
      if(l(1,n1)-l(2,n2)) 33,50,35
   33 n1=n1+1
      go to 20
   35 n2=n2+1
      go to 20
c
   40 if(n1.gt.nf1.or.n2.gt.nf2) go to 80
      if(nord(1,n1)-nord(2,n2)) 43,55,45
   43 n1=n1+1
      go to 20
   45 n2=n2+1
      go to 20
c
   50 if(isort1) 40,40,60
c
   55 if(isort1) 60,60,30
c
c  now l and nord should agree. set difference and possibly error
c
   60 if(idiff.eq.1) then
c
c  relative difference
c
        dfrq=frq(2,n2)/frq(1,n1)-1
c
c  test for error
c
        errsq=0
        if(icsin(1).eq.6) errsq=errsq+(error(1,n1)/frq(1,n1))**2
        if(icsin(2).eq.6) errsq=errsq+(error(2,n2)/frq(2,n2))**2
        err=sqrt(errsq)
      else if(idiff.eq.2) then
c  absolute difference
        dfrq=frq(2,n2)-frq(1,n1)
c
c  test for error
c
        errsq=0
        if(icsin(1).eq.6) errsq=errsq+error(1,n1)**2
        if(icsin(2).eq.6) errsq=errsq+error(2,n2)**2
        err=sqrt(errsq)
      else
c  relative difference, dimensionless frequency
	dfrq=sig(2,n2)/sig(1,n1)-1
        err=0
      end if
c
c  test for energy scaling of difference and possibly errors
c
   65 if(iscdif.ge.1) then
        dfrq=ekins(n1)*dfrq
        err =ekins(n1)*err
      end if
c
c  test for inverting difference
c
      if(invdif.eq.1) dfrq=-dfrq
c
c  output, depending on whether or not errors are included
c
      if(isterr.eq.0) then
        if(iprdif.eq.1)
     *    write(6,150) n1,n2,l(1,n1),l(2,n2),nord(1,n1),nord(2,n2),
     *    frq(1,n1),frq(2,n2),dfrq
        write(20,152) l(1,n1),nord(1,n1),
     *    frq(1,n1),dfrq,frq(1,n1)/(l(1,n1)+0.5)
c
      else
c
        if(iprdif.eq.1)
     *    write(6,150) n1,n2,l(1,n1),l(2,n2),nord(1,n1),nord(2,n2),
     *    frq(1,n1),frq(2,n2),dfrq,err
        write(20,152) l(1,n1),nord(1,n1),
     *    frq(1,n1),dfrq,frq(1,n1)/(l(1,n1)+0.5),err
      end if
c
c  output energy differences
c
      if(ekin(1,n1).ne.0.and.ekin(2,n2).ne.0) then
        dekin=(ekin(2,n2)-ekin(1,n1))/ekin(1,n1)
        write(30,152) l(1,n1), nord(1,n1),frq(1,n1),dekin
      end if 
c
c  increment n1 and try again
c
   68 n1=n1+1
      go to 20
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
  117 format(//' inverted difference. (first set) - (second set)')
  118 format('#'/'# inverted difference. (first set) - (second set)')
  121 format(///' n1,n2,l1,l2,ord1,ord2,cyclic frq1,frq2 (muhz),',
     *  'relative difference (frq2 - frq1)',a/)
  122 format(///' n1,n2,l1,l2,ord1,ord2,dimensionless frq1,frq2 ,',
     *  'relative difference (frq2 - frq1)',a/)
  123 format(///' n1,n2,l1,l2,ord1,ord2,cyclic frq1,frq2 (muhz),',
     *  'absolute difference (frq2 - frq1) (muhz)',a/)
  124 format(///' n1,n2,l1,l2,ord1,ord2,dimensionless frq1,frq2 ,',
     *  'absolute difference (frq2 - frq1)',a/)
  125 format(///' n1,n2,l1,l2,ord1,ord2,cyclic frq1,frq2 (muhz),',
     *  'relative diff., dimensionless (frq2 - frq1) (muhz)',a/)
  126 format(///' n1,n2,l1,l2,ord1,ord2,dimensionless frq1,frq2 ,',
     *  'relative diff., dimensionless (frq2 - frq1) (muhz)',a/)
  131 format('#'/'#'/'# l1,ord1,cyclic frq1 (muhz),',
     *  'rel. diff. (frq2 - frq1), ',
     *  '(cyclic frq1)/(l+1/2)',a/'#')
  132 format('#'/'#'/'# l1,ord1,dimensionless frq1 ,',
     *  'rel. diff. (frq2 - frq1), ',
     *  '(dimensionless frq1)/(l+1/2)',a/'#')
  133 format('#'/'#'/'# l1,ord1,cyclic frq1 (muhz),',
     *  'abs. diff. (frq2 - frq1) (muhz), ',
     *  '(cyclic frq1)/(l+1/2)',a/'#')
  134 format('#'/'#'/'# l1,ord1,dimensionless frq1 ,',
     *  'abs. diff. (frq2 - frq1)',
     *  '(dimensionless frq1)/(l+1/2)',a/'#')
  135 format('#'/'#'/'# l1,ord1,cyclic frq1 (muhz),',
     *  'rel. diff., dimensionless (frq2 - frq1), ',
     *  '(cyclic frq1)/(l+1/2)',a/'#')
  136 format('#'/'#'/'# l1,ord1,dimensionless frq1 ,',
     *  'rel. diff., dimensionless (frq2 - frq1), ',
     *  '(dimensionless frq1)/(l+1/2)',a/'#')
  140 format('#'/'# l1, ord1, cyclic frq1 (muhz), ',
     *  'rel. diff. (ekin2 - ekin1)',a/'#')
  150 format(6i5,0p2f12.3,1p2e13.5)
  152 format(i5,i6,0pf12.3,1p3e13.5)
  155 format(//' **** np =',i5,' exceeds nptmax =',i5,
     *  ' k, l =',2i5)
  170 format(4(1h$),1h',i4,1h',38(1h$))
  172 format(48(1h$))
      end
      subroutine scekin(l,nord,frq,ekin,nn,init,iscdif,lnorm,frnorm)
c
c  resets ekin to ratio of energies at fixed frequency
c
c  if init = 1, or in first call, ratio is taken to modes
c  with degree lnorm or, if lnorm .le. 0, 
c  the first degree in the array l(1,.).
c  if init .ne. 1, ratio is taken to energy assumed to be
c  already in e1(.), as a function of frequency in fr1(.).
c
      implicit double precision (a-h, o-z)
      dimension l(2,1),nord(2,1),frq(2,1),ekin(1),fr1(100),
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
        if(l(1,n).eq.lnorm) then
          n1=n
          go to 8
        end if
    5   continue
c
    8   if(n1.eq.0) then
          write(6,100) lnorm,l(1,1)
          lnorm=l(1,1)
          n1=1
        end if
c
      else
        n1=1
        if(l(1,1).ne.0) then
          write(6,100) lnorm,l(1,1)
          lnorm=l(1,1)
        end if
c
      end if
c
      lp=l(1,n1)
      nne=0
      ns=0
      do 10 n=n1,nn
      if(l(1,n).ne.lp) go to 15
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
        frp=1.e10
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
      nnw=min0(nn,1000)
      if(iscdif.eq.1) then
        write(6,105)
      else 
        write(6,110)
      end if
c
      write(6,120) (n,l(1,n),nord(1,n),frq(1,n),ekin(n),n=1,nnw)
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
      subroutine rdofrq
      return
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
      dimension l(2,1),nord(2,1),frq(2,1),ekin(1)
c
      gm1=1.6666667d0
c
c  test for setting xftc and frftc
c
c..	write(6,*) 'xnorm,xftc,frftc,xs,vgs',xnorm,xftc,frftc,xs,vgs
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
      xx=1-xx+4*l(1,n)*(l(1,n)+1)*(1-4*(1-1/gm1)/(gm1*xx))/(v*v)
      if(xx.lt.0) then
        arat=1
        write(6,120) l(1,n), nord(1,n), frq(1,n)
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
