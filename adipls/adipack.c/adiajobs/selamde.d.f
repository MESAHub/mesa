      program selsum
c
c  program to filter file of adiabatic eigenfunctions or kernels.
c
c  option for selecting modes corresponding to those in dataset
c  of grand summaries, short summaries or observed frequencies. 
c  all mode datasets must be sorted with order varying fastest.
c
c  outputs eigenfunction or kernel on d/s 10.
c  outputs grand summary on d/s 11.
c  outputs short summary on d/2 12.
c  output l, order and frequency on d/s 13, in same form as observed
c  results
c
c  if iselct .gt. 0, also outputs the corresponding data for
c  the selection set of modes on d/s 21, 22 and 23.
c
c  original version (based on programme selsum): 02/6/88.
c
      implicit double precision(a-h, o-z)
      character*280 form,ranges
      dimension cs(50),cssl(50),ics(24),icssl(24),nrange(20),
     *  ss(7),iss(2),sssl(7),isssl(2),ssmod(7),
     *  ranges(3),form(3),x(5001),y(6,5001)
      equivalence (cs(39),ics(1))
      equivalence (cssl(39),icssl(1))
      equivalence (ss(6),iss(1))
      equivalence (sssl(6),isssl(1))
c
      data idsinp, idsslp /-1, -1/
      data initrd /1/
      data ranges /'record number','mode order','mode degree'/
      data form /'unformatted','unformatted','formatted'/
c
c..      namelist /exec/ idsin,idssl,nw1,nw2,lw1,lw2,frqw1,frqw2,
c..     *  sigw1,sigw2,iselct,irewin,irewsl,irange,nrnrec,nrnord,idiag
c
c  defaults in exec
c
c  idsin: dataset containing input modes
      idsin=2
c  idssl: dataset containing modes for selection.
      idssl=3
c  icase: type of eigenfunction.
c     icase = 1: Full eigenfunction, as output from adipls
c     icase = 2: Reduced eigenfunction. First record gives nn and x,
c                subsequenct records give the eigenfunctions, in terms
c                of two variables.
c     icase = 3: Kernels for spherically symmetric rotation.
      icase = 1
c  icassl: type of modes used for selection.
c  icassl = 1: grand summary.
c  icassl = 2: short summary.
c  icassl = 3: observed frequencies.
      icassl=1
c  nw1, nw2, lw1, lw2, frqw1, frqw2, sigw1, sigw2: windows in
c  order, l, frequency (in microhz), and dimensionless squared
c  frequency.
c  window in, e.g., order is only applied if nw1 .le. nw2
c  note, however, that if lw2 = -2, modes with l = lw1 are selected.
      nw1=0
      nw2=-1
      lw1=0
      lw2=-1
      frqw1=0
      frqw2=-1
      sigw1=0
      sigw2=-1
c  iselct: if iselct = 1 select only modes matching modes from
c          file on d/s idssl of grand summaries, short summaries or
c          observed frequencies, depending on value of icassl.
      iselct=0
c  irewin: if irewin = 1 rewind input data set
      irewin=0
c  irewsl: if irewsl = 1 rewind selection data set
      irewsl=0
c  irange, nrange: selection of ranges in record number,
c     mode order or degree. ranges are
c     in the array nrange(1-20), which is initialized to zero
c     before /exec/ is read. see function inrnge for definition
c     of ranges.
c  irange = 1: impose ranges in input frequency record number.
c  irange = 2: impose ranges on mode orders. ranges are in the
c     array nrange. to allow modes with order .le. 0 (i.e. g and f modes)
c     the range is given as 1000 + order.
c     thus nrange = 990,-995,1000,1005,-1008 specifies the modes
c     with orders = -10 to -5, 0, and 5 to 8.
c  irange = 3: impose ranges on mode degree. 
      irange=0
c  idiag: if idiag .gt. 0, details are printed about modes
c     selected.
      idiag=1
c
c  .................................................................
c
c  initialize flags and counters for input and selection data sets
c
      nout=0
      idsinp=-1
      idsslp=-1
      initsl=0
      iselop=0
      isobin=0
      isobsl=0
c
      write(6,*) 'Files needed:'
      write(6,*) 'Input eigenfunctions on idsin (default 2)'
      write(6,*) 'Input of selection modes on idssl (default 3)'
      write(6,*) ' outputs eigenfunctions on d/s 10'
      write(6,*) ' outputs grand summary on d/s 11.'
      write(6,*) ' outputs short summary on d/s 12.'
      write(6,*) ' output l, order and frequency on d/s 13, '
      write(6,*) ' in same form as observed results'
      write(6,*) 
     *   'If selecting from other mode file, output grand summary,'
      write(6,*) 
     *   'short summary and observed results from selection set'
      write(6,*) 'on d/s 21, 22, 23.'
c
c  set up files
c
      call ofiles
c
c  open fixed datasets
c
      call openf(10,'u','u')
c
      do 3 i=1,3
      id=10+i
    3 call openf(id,'n',form(i))
c
c  zero range
c
    5 do 6 k=1,20
    6 nrange(k)=0
c
      nrange(1)=-1
c
      initr1=1
      initr2=1
      initr3=1
c
c..      read(5,exec,end=90)
c
      write(6,*) 'idsin,idssl?'
      write(6,*) idsin,idssl
      read(5,*,end=90) idsin,idssl
      write(6,*) 'icase,icassl?'
      write(6,*) icase,icassl
      read(5,*,end=90) icase,icassl
      write(6,*) 'nw1,nw2,lw1,lw2,frqw1,frqw2, sigw1, sigw2?'
      write(6,*) nw1,nw2,lw1,lw2,frqw1,frqw2, sigw1, sigw2
      read(5,*) nw1,nw2,lw1,lw2,frqw1,frqw2, sigw1, sigw2
      write(6,*) 'iselct,irewin,irewsl,irange?'
      write(6,*) iselct,irewin,irewsl,irange
      read(5,*) iselct,irewin,irewsl,irange
c
c  test for reading ranges
c
      if(irange.gt.0) then
        write(6,*) 'ranges in ',ranges(irange)
        write(6,*) (nrange(i),i=1,20)
        read(5,*) (nrange(i),i=1,20)
      end if
c
      write(6,*) 'idiag?'
      write(6,*) idiag
      read(5,*) idiag
c
c  output /exec/
c
      write(6,*) 'idsin,idssl'
      write(6,*) idsin,idssl
      write(6,*) 'icase,icassl'
      write(6,*) icase,icassl
      write(6,*) 'nw1,nw2,lw1,lw2,frqw1,frqw2, sigw1, sigw2'
      write(6,*) nw1,nw2,lw1,lw2,frqw1,frqw2,sigw1,sigw2
      write(6,*) 'iselct,irewin,irewsl,irange'
      write(6,*) iselct,irewin,irewsl,irange
c
c  test for ranges
c
      if(irange.ge.1) then
        write(6,*) 'ranges in ',ranges(irange)
        write(6,*) (nrange(i),i=1,20)
      end if
c
      write(6,*) 'idiag'
      write(6,*) idiag
      write(6,*)
      write(6,*) '******************'
      write(6,*)
c
c..      write(6,100)
c..      write(6,exec)
c
c  set number of input variables
c
      if(icase.eq.1) then
        ivar = 6
      else if(icase.eq.2) then
        ivar = 2
      else if(icase.eq.3) then
        ivar = 1
      end if
c
c  test for opening input files
c
      if(idsin.ne.idsinp) then
	write(6,*) 'Open input on d/s',idsin
        if(idsinp.gt.0) close(idsinp)
        call openf(idsin,'o','u')
        idsinp=idsin
c
c  for icase = 2, input, and possibly output, record with nn and x
c
        if(icase.eq.2) then
          read(idsin) nn,(x(n),n=1,nn)
	  write(6,*) 'End reading x. nn =',nn
          if(initrd.eq.1) then
            write(10) nn,(x(n),n=1,nn)
            initrd=0
          end if
        end if
c
      end if
c
      if(idssl.ne.idsslp.and.iselct.gt.0) then
        if(idsslp.gt.0) close(idsslp)
        call openf(idssl,'o',form(icassl))
        idsslp=idssl
      end if
c
c  test for opening output to selection dataset
c
      if(iselct.eq.1.and.iselop.eq.0) then
        do 6010 i=icassl,3
        id=20+i
 6010   call openf(id,'n',form(i))
c
        iselop=1
      end if
c
      if(idiag.gt.0) write(6,110)
c
      if(idsin.ne.idsinp) nin=0
      if(idssl.ne.idsslp) nsl=0
      idsinp=idsin
      idsslp=idssl
c
      lw11=lw1
      lw21=lw2
      if(lw2.eq.-2) lw21=lw1
c
c  test for rewind of input data set
c
      if(irewin.eq.1) then
        rewind idsin
        nin=0
      end if
c
c  test for rewind of selection data set
c
    7 if(iselct.eq.1.and.(initsl.ne.1.or.irewsl.eq.1)) then
        rewind idssl
        nsl=0
        initsl=1
        nords=-1
        ls=-1
c
c..        if(icassl.eq.3) read(idssl,115,end=55,err=55) lobsin,nnno
      end if
c
    8 lobsin=-2
      nnno=0
c
c  test for output of header for observational datasets
c
      if(isobin.eq.0) then
c..        write(13,115) lobsin,nnno
        isobin=1
      end if
c
      if(iselct.eq.1.and.isobsl.eq.0) then
c..        write(23,115) lobsin,nnno
        isobsl=1
      end if
c
c  read input record
c
   10 if(icase.ne.2) then
        read(idsin,end=50) cs,nn,(x(n),(y(i,n),i=1,ivar),n=1,nn)
      else
        read(idsin,end=50) cs,((y(i,n),i=1,2),n=1,nn)
      end if
c
      nin=nin+1
c
      l = nint(cs(18))
      nord = nint(cs(19))
      sig = cs(20)
      if(cs(27).gt.0) then
        frq = 1000.d0*cs(27)
      else
        frq = 16666.66667d0/cs(25)
      end if
c
c  test for range in record number
c
      if(irange.ne.1) go to 10500
      inr=inrnge(nin,nrange,20,initr1)
      initr1=0
      if(inr) 60,10,10500
c
c  windowing in l
c
10500 if(lw11.gt.lw21) go to 11
      if(l.lt.lw11) go to 10
      if(l.gt.lw21) go to 60
c
c  windowing in order, frequency and dimensionless squared
c  squared frquency 
c
   11 if(nw1.le.nw2.and.(nord.lt.nw1.or.nord.gt.nw2)) go to 10
      if(frqw1.le.frqw2.and.(frq.lt.frqw1.or.frq.gt.frqw2)) go to 10
      if(sigw1.le.sigw2.and.(sig.lt.sigw1.or.sig.gt.sigw2)) go to 10
c
c  test for range in mode order
c
      if(irange.ne.2) go to 11300
      inr=inrnge(nord+1000,nrange,20,initr2)
      initr2=0
      if(inr) 60,10,11300
c
c  test for range in mode degree
c
11300 if(irange.ne.3) go to 11500
      inr=inrnge(l,nrange,20,initr3)
      initr3=0
      if(inr) 60,10,11500
c
c  test for selection based on mode file
c
11500 if(iselct.ne.1) go to 30
c
c  test on l
c
   12 continue
c
      if(l-ls) 10,15,20
c
c  test on order
c
   15 if(nord-nords) 10,30,20
c
c  read next mode on selection data set
c
   20 call rdfreq(icassl,idssl,cssl,ls,nords,sigs,frqs,ekins,ierr)
c
      if(ierr.gt.0) go to 55
c
      if(l.ge.0) then
        nsl=nsl+1
      else if(icassl.eq.2) then
c
c  model record in short summary
c
        write(22) (cssl(i),i=1,7)
        go to 20
      end if
c
      go to 12
c
c  now mode has been located. output.
c
   30 nout=nout+1
c
c  test for output of selection mode
c
      if(iselct.eq.0) go to 40
c
c  output selection data sets
c
      if(icassl.eq.1) then
        write(21) cssl
c
c  set selection short summary
c
        call setssm(cssl,icssl,sssl,isssl,ssmod,irmod)
        if(irmod.eq.1) write(22) ssmod
c
        write(22) sssl
c
      else if (icassl.eq.2) then
        write(22) (cssl(i),i=1,7)
c
      end if
c
      write(23,120) l,nord,frqs
c
c  output eigenfunction
c
   40 if(icase.ne.2) then
        write(10) cs,nn,(x(n),(y(i,n),i=1,ivar),n=1,nn)
      else
        write(10) cs,((y(i,n),i=1,2),n=1,nn)
      end if
c
c  output summaries
c
      write(11) cs
c
c  set original short summary
c
      call setssm(cs,ics,ss,iss,ssmod,irmod)
      if(irmod.eq.1) write(12) ssmod
c
      write(12) ss
c
      write(13,120) l,nord,frq
c
      if(idiag.gt.0) write(6,130) nin,nsl,nout,l,nord,sig,frq
c
      go to 10
c
c
c  diagnostics for end of file on input or selection data sets
c
   50 write(6,140) idsin,nin
      go to 5
c
   55 write(6,145) idssl,nsl
      go to 5
c
c  exceeding window in l
c
   60 backspace idsin
      nin=nin-1
      go to 5
c
   90 continue
      stop
  100 format(1h1)
  110 format(///' modes selected. nin, nsl, nout, l, order, sigma**2,',
     *  ' frequency:'/)
  115 format(i5)
  120 format(2i5,f8.2)
  130 format(5i5,2f12.4)
  140 format(//' ***** eof reached on input data set',i3,' after',
     *  ' record no',i5)
  145 format(//' ***** eof reached on selection data set',i3,' after',
     *  ' record no',i5)
      end
