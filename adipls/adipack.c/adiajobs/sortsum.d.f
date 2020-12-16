      program ordsum
c
c  merge and order summary (either grand or short)
c
c  file of summaries can contain data for several models,
c  but the data for each model must appear together.
c
c  data read from d/s iw (defined in /exec/, default 10)
c  grand summary output on d/s 2
c
c  Modified 10/12/93, to allow sorting after degree first, and
c  then eigenfrequency
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 11/4/90
c
      implicit double precision (a-h, o-z)
      parameter (ndmax=10000)
c
      dimension icase(ndmax),s1(ndmax),s2(ndmax),
     *  s3(ndmax),datg(50,ndmax),gs(50),igs(12),iss(2),
     *  gsmp(8),ssmod(7),
     .  sn(ndmax),ihh(ndmax),na(ndmax),mtn(ndmax),mode(50)
c  space saving equivalence
      equivalence (s3(1),ihh(1))
      equivalence (gs(39),igs(1))
      equivalence (gs(6),iss(1))
c
c..      namelist /exec/ icasrd,iw,elw1,elw2,iordw1,iordw2,
c..     *  mode,
c..     *  isort,igncse,incl,eps,epssol,
c..     *  icorcl
c
      data gsmp /8*-1.d0/
      data ssmod /7*0.d0/
c  defaults
c  icasrd: determines type of read
c  icasrd = 1: grand summaries
c  icasrd = 2: short summaries
      icasrd=1
c
      iw=10
c  window in l and order
      elw1=0.d0
      elw2=1.d5
      iordw1=-10000
      iordw2=10000
c  mode: (initialized to zero before each read of /exec/)
c        determines modes to be excluded. single modes may
c        be specified by their number. a range of modes, n1 to n2, say,
c        is indicated by setting mode(i) = n1, mode(i+1) = -n2.
c        thus mode = 2,4,5,-8,10,... excludes modes no 2,4,5,6,7,8,10,...
c
c  abs(isort) = 1: order after l first and then order
c  abs(isort) = 2: order after order first and then l
c  abs(isort) = 3: order after sig
c  abs(isort) = 4: order after l first and then eigenfrequency
c  isort .gt. 0: order increasingly
c  isort .lt. 0: order decreasingly
      isort=1
c  igncse: if igncse = 1 ignore case number when testing for same mode
      igncse=0
c  when the same mode occurs twice, but with different sig or variational
c  frequency (to accuracy eps) and/or ekin (to accuracy epssol)
c  the action taken depends on incl:
c   incl = 0  -  include both results
c   incl = 1  -  include first result
c   incl = 2  -  include second result
c  if the modes are deemed to be identical, the second one is
c  taken unless incl = 1.
      eps=1.d-8
      epssol=5.d-3
      incl=2
c  icorcl: if icorcl = 1, correct labelling of low-order l = 1 modes
c  assumed to be set incorrectly by adiabatic pulsation programme
c  for present sun.
      icorcl=0
c
c                  ---------------------------------------
c
      n1=1
      nr=0
      ifin=0
      id=11
c
c  set files
c
      write(6,*) 'Files needed:'
      write(6,*) 'Input from unit iw (default 10)'
      write(6,*) 'Output to unit 2'
c
      call ofiles
c
      call openf(2,'u','u')
c
c  zero mode before next read of exec
    5 do 7 k=1,50
    7 mode(k)=0
c
c
c..      read(5,exec,end=30,err=30)
      write(6,*) 'icasrd,iw,elw1,elw2,iordw1,iordw2?'
      write(6,*) icasrd,iw,elw1,elw2,iordw1,iordw2
      read(5,*,end=30,err=30) icasrd,iw,elw1,elw2,iordw1,iordw2
      write(6,*) 'mode?'
      call iprsgn(6,mode,50,' ')
      call rdilst(5,mode,50,ierr)
      if(ierr.gt.0) go to 30
      write(6,*) 'isort,igncse,incl,eps,epssol?'
      write(6,*) isort,igncse,incl,eps,epssol
      read(5,*,end=30,err=30) isort,igncse,incl,eps,epssol
      write(6,*) 'icorcl?'
      write(6,*) icorcl
      read(5,*,end=30,err=30) icorcl
c..      write(6,100)
c..      write(6,exec)
      write(6,*) 'icasrd,iw,elw1,elw2,iordw1,iordw2'
      write(6,*) icasrd,iw,elw1,elw2,iordw1,iordw2
      write(6,*) 'mode'
      call iprsgn(6,mode,50,' ')
      write(6,*) 'isort,igncse,incl,eps,epssol'
      write(6,*) isort,igncse,incl,eps,epssol
      write(6,*) 'icorcl'
      write(6,*) icorcl
c
      km=0
      nm1=0
      if(mode(1).le.0) nm1=1000000
      nm2=nm1
c
c  set parameters depending on icasrd
c
      if(icasrd.eq.1) then
        invar=50
        ivarl=18
        ivarn=19
        ivars=20
	ivare=24
	ivarf=27
      else
        invar=7
        ivarl=1
        ivarn=2
        ivars=3
	ivare=4
	ivarf=5
      end if
c
c  open input file
c
      call openf(iw,'o','u')
      rewind iw
c  read data
   10 ndmx=ndmax-nr
      if(ndmx.le.0) go to 20
c
c  read grand summaries
c
   15 nd=0
      nm=0
c
      ssmod(1)=0
c
c  read next mode
c
   16 inmod=0
      call rdfreq(icasrd,iw,gs,l,nord,sig,frq,ekin,ierr)
      if(ierr.ne.0) go to 18
c
c  test for new model, depending on icasrd
c
      if(icasrd.eq.1) then
        do 16050 k=2,5
        if(abs(gs(k)/gsmp(k)-1).gt.1.d-6) inmod=1
16050   gsmp(k)=gs(k)
c
      else if(l.lt.0) then
        do 16055 k=3,6
        if(abs(gs(k)/gsmp(k)-1).gt.1.d-6) inmod=1
16055   gsmp(k)=gs(k)
c
c  store model record for later output
c
        do 16057 k=1,7
16057   ssmod(k)=gs(k)
c
      end if
c
      if(inmod.eq.1.and.nm.gt.0) then
        backspace iw
        nr=nr+nd
        go to 35
      end if
c
      nm=nm+1
c  test for excluding modes
      if(nm.le.nm2) go to 1615
c  next value in mode
      km=km+1
      nm1=mode(km)
      nmm=mode(km+1)
      if(nmm.ge.0) go to 1612
c  range of modes
      km=km+1
      nm2=-nmm
      go to 1614
c  single mode
 1612 if(nm1.eq.0) nm1=1000000
      nm2=nm1
c
 1614 if(nm1.le.999999) write(6,105) nm1,nm2
c  test on nm
 1615 if(nm1.le.nm.and.nm.le.nm2) go to 16
c
      do 1691 i=1,invar
 1691 datg(i,n1)=gs(i)
c
c  set case number, depending on icasrd
c
      if(icasrd.eq.1) then
        icase(n1)=igs(5)
      else
        icase(n1)=iss(1)
      end if
c
      n1=n1+1
      nd=nd+1
      go to 16
c
   18 nr=nr+nd
      go to 5
c
c  data array is full. order and output results read so far.
c
   20 write(6,110) ndmax
      go to 35
c
c               -------------------------------------
c
c  no more input files. order and output
c
   30 ifin=1
c
   35 isor1=iabs(isort)
      m=0
      do 48 n=1,nr
c  grand summary
      el=datg(ivarl,n)
      ord=datg(ivarn,n)
      sig=datg(ivars,n)
c  correct classification for l=1, order=0,1
   38 if(icorcl.eq.1.and.el.eq.1.and.ord.ge.0.and.ord.lt.2) then
        ord=ord+1
        datg(ivarn,n)=ord
        write(6,102) ord
      end if
c  test that mode is in window
      if(el.ge.elw1.and.el.le.elw2.and.ord.ge.iordw1.and.ord.le.iordw2)
     .  then
        m=m+1
	if(isor1.eq.1) then
c  order after l first and then order
          s1(m)=el
          s2(m)=ord
	else if(isor1.eq.2) then
c  order after order first and then l
          s1(m)=ord
          s2(m)=el
	else if(isor1.eq.3) then
c  order after sig
          s1(m)=sig
          s2(m)=0
          s3(m)=0
	else if(isor1.eq.4) then
c  order after l first and then eigenfrequency
          s1(m)=el
          s2(m)=sig
        else
	  write(6,107) isor1
	  stop
        end if
c
        sn(m)=m
        mtn(m)=n
      end if
c
   48 continue
c  reset nr
      nr=m
c  sort after s1
      call sort(s1,sn,nr)
c  when sorting after sig, this is all
      if(isor1.eq.3) go to 60
c  rearrange s2 into s3
      do 50 n=1,nr
      n1=sn(n)
   50 s3(n)=s2(n1)
c  now sort after s3 within each group of fixed s1
      i1=1
   52 ss1=s1(i1)
      ssa= abs(ss1)+1.d-10
      i=i1
   54 i=i+1
      if(i.gt.nr) go to 55
      if(s1(i)-ss1.lt.1.d-15*ssa) go to 54
c  now group has been found. length of group
   55 i2=i-1
      il=i-i1
      if(il.gt.1) then
c  sort after s3
        call sort(s3(i1),sn(i1),il)
      end if
      i1=i
      if(i.le.nr) go to 52
c
c  sorting is now complete. commence output section
c
   60 nr1=nr+1
      m=0
      do 80 n1=1,nr
      n2=n1
      if(isort.lt.0) n2=nr1-n1
      n=sn(n2)
      ni=mtn(n)
      ic=icase(ni)
c
      if(isor1-3) 61,63,61
c
   61 ss1=s1(n2)
      ss2=s2(n)
      go to 65
c
   63 ss1=datg(ivarl,ni)
      ss2=datg(ivarn,ni)
c
c  test for same mode as previous
c
   65 if(n1.eq.1) go to 77
      ih=0
c
c  test for case number if igncse .ne. 1
c
      if(ic.ne.icp.and.igncse.ne.1) go to 75
      ds1= abs(ss1-ss1p)/( abs(ss1p)+1.d-10)
      if(ds1.gt.1.d-15) go to 74
c  same s1. test for same s2
      ds2= abs(ss2-ss2p)/( abs(ss2p)+1.d-10)
      if(ds2.gt.1.d-15) go to 77
c
c  same frequency?
c
  653 sig=datg(ivars,ni)
      frqv=datg(ivarf,ni)
      ekin=datg(ivare,ni)
      sigp=datg(ivars,np)
      frqvp=datg(ivarf,np)
      ekinp=datg(ivare,np)
c
  655 dsig= abs(sig-sigp)/( abs(sigp)+1.d-10)
      dfrqv= abs(frqv-frqvp)/( abs(frqvp)+1.d-5)
      dekin= abs(ekin-ekinp)/( abs(ekinp)+1.d-15)
c
c  if frequency and energy is same as for previous mode, skip mode
c
      if(dsig.ge.eps.or.dfrqv.ge.eps.or.dekin.ge.epssol)
     *  go to 660
c
c  select mode depending on incl
c
      if(incl-1) 69,67,69
c
c  otherwise print diagnostic
c
  660 write(6,120) np,ni
      call prtsum(6,datg(1,np),icasrd)
      call prtsum(6,datg(1,ni),icasrd)
c  include this and/or previous mode, depending on incl
      if(incl-1) 77,67,69
   67 if(ni-np) 73,73,80
   69 if(np-ni) 73,73,80
   73 m=m-1
c
c  test for change in case number (in case igncse was 1)
c
73500 if(ic.ne.icp) go to 75
c
      go to 77
c
c  change in icase or s1. print heading
c
   74 if(isor1.eq.3) go to 77
   75 ih=1
c  set controlling arrays
   77 m=m+1
      na(m)=ni
      ihh(m)=ih
      ss1p=ss1
      ss2p=ss2
      icp=ic
      np=ni
   80 continue
c  reset nr
      nr=m
c
c  output summary on d/s 2
c
c  test for output of model record
c
   86 if(icasrd.ne.1.and.ssmod(1).ne.0) write(2) ssmod
      do 87 m=1,nr
   87 write(2) (datg(i,na(m)),i=1,invar)
c
c            --------------------------------------
c
c  test for new model, or whether all data has been read
c
   88 if(ifin.eq.1.and.inmod.eq.0) go to 90
      n1=1
      nr=0
      go to 10
c
   90 continue
      stop
  100 format(////)
  102 format(///' for l = 1, order has been incremented to',f10.2)
  105 format(/' exclude modes',i5,'  to',i5)
  107 format(/' ***** Error in sortsum. isort = ',i5,' not defined')
  110 format(//1x,10(1h*),' number of modes on file greater than',
     .  i5,'.  order modes read so far')
  120 format(//1x,130(1h*)//' modes no.',i5,'  and',i5,
     .  '  have same classification but different frequencies ',
     .  'or energies. the modes are')
      end
      subroutine prtsum(iw,gs,icasrd)
c
c  outputs summary in array gs on iw
c  for icasrd = 1: grand summary
c  for icasrd = 2: short summary
c
      implicit double precision (a-h, o-z)
      dimension gs(*)
c
      if(icasrd.eq.1) then
	write(iw,100) gs(18),gs(19),gs(20),gs(24),gs(27)
      else
	write(iw,100) (gs(i),i=1,5)
      end if
      return
  100 format(/' l, order, sigma**2, ekin, nu:'/
     *  f10.2,f6.1,1p3e13.5)
      end
