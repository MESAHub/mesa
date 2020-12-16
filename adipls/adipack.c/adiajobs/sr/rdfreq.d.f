      subroutine rdfreq(icase,ids,cs,l,nord,sig,frq,ekin,ierr)
c
c  read frequency data from d/s ids. case depends on icase.
c
c  icase = 1: grand summary. returned in cs(1-50). frequency from
c             variational frequency or, if not available from
c             eigenfrequency, possibly corrected for perturbation
c             in gravitational potential
c  icase = 2: short summary. returned in cs(1-7)
c  icase = 3: observed frequencies.
c  icase = 4: grand summary. returned in cs(1-50). sig and frq 
c             obtained from eigenfrequency in cs(20).
c             Note that this allows setting Cowling 
c             approximation frequency.
c  icase = 5: grand summary. returned in cs(1-50). frq obtained from
c             Richardson extrapolation frequency in cs(37), 
c             if this is set. Otherwise variational frequency is used.
c  icase = 6: grand summary. returned in cs(1-50). sig and frq 
c             obtained from eigenfrequency in cs(21),
c             hence possibly corrected for effect of 
c             gravitational potential, for Cowling calculation.
c
c  if icase .lt. 0, read from single-precision dataset, according
c  to value of abs(icase)
c
c  returns mode degree, order, dimensionless sigma**2, frequency
c  (in microHz) and dimensionless energy in l, nord, sig, frq and ekin.
c
c  Note: for icase = 1 returns corrected sigma**2, if frequency
c  calculation was in Cowling approximation. This is normally
c  also what is set in the short summary.
c
c  ierr is returned as 1 or 2, on end of file or error in read.
c
c  If data relevant for icase are not found, other data are used in the
c  order of preference: variational frequency, eigenfrequency. In that
c  case, ierr is returned as -icase_actual, where icase_actual is the
c  value of icase corresponding to the frequency used.
c
c  original version: 26/9/86.
c
c  Modified 1/1/08, to return negative ierr if proper data not 
c  available. So far implemented only for icase = 5.
c
c                     ................................
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      real css, sss
      dimension cs(*),csd(50),icsd(8), ssd(7),issd(2),
     *                css(50),icss(8), sss(7),isss(2)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      equivalence(csd(39), icsd(1)),(ssd(6), issd(1))
      equivalence(css(39), icss(1)),(sss(6), isss(1))
c
      data initrd, idsp /0, -1/
      data icerr /0/
c
      cgrav=6.67232d-8
      twopi=8.d0*atan(1.d0)
c
      icasea=iabs(icase)
c
      ierr = 0
c
      if(icasea.eq.2) then
c
c  short summary
c
        if(icase.eq.2) then
          read(ids,end=92,err=95) (cs(i),i=1,7)
        else
c
c  single precision read
c
          read(ids,end=92,err=95) (css(i),i=1,7)
          do 12 i=1,5
   12     csd(i)=css(i)
          do 14 i=1,2
   14     icsd(i)=icss(i)
          do 16 i=1,7
   16     cs(i)=csd(i)
        end if
c
        l=nint(cs(1))
        nord=nint(cs(2))
        sig=cs(3)
        frq=1000*cs(5)
        ekin=cs(4)
c
        go to 90
c
      else if(icasea.eq.3) then
c
c  observed frequencies
c
c  test for first read, file in the format for rdofrq
c
        if(initrd.eq.0.or.ids.ne.idsp) then
          call skpcom(ids)
          read(ids,*,end=92,err=95) xlobsn
          if(xlobsn.lt.0) then
            read(ids,*,end=92,err=95) xnnobs
          else
            backspace ids
          end if
c
          idsp=ids
          initrd=1
        end if
c
        read(ids,*,end=92,err=95) l,nord,frq
        sig=0
        ekin=0
c
        go to 90
c
      end if
c
c  other options assume grand summary input
c
c  read grand summary
c
   20 if(icase.gt.0) then
   22   read(ids,end=92,err=95) (cs(i),i=1,50)
c
c  test for proper summary (rather than x-record in mode file)
c
        if(cs(2).le.1) go to 22
      else
c
c  single precision read
c
   24   read(ids,end=92,err=95) (css(i),i=1,50)
c
c  test for proper summary (rather than x-record in mode file)
c
        if(css(2).le.1) go to 24
c
        do 32 i=1,38
   32   csd(i)=css(i)
        do 34 i=1,8
   34   icsd(i)=icss(i)
        do 36 i=1,50
   36   cs(i)=csd(i)
c
      end if
c
      l=nint(cs(18))
      nord=nint(cs(19))
      ekin=cs(24)
c
c  setting of frequencies depend on value of icasea
c
      if(icasea.eq.1) then
c
c  use variational frequency, if available
c
        sig=cs(21)
        if(sig.le.0) sig=cs(20)
        frq=1000*cs(27)
        if(frq.le.0) then
          frq=sqrt(sig/cs(20))*16666.6666667d0/cs(25)
          frq1=1.d6/twopi*sqrt(cgrav*cs(2)*sig/cs(3)**3)
          if(abs(frq/frq1-1).gt.1.e-5.and.icerr.le.20) then
            icerr=icerr+1
            if(istdpr.gt.0) write(istdpr,110) frq, frq1
            if(icerr.le.2.and.istdpr.gt.0) 
     *        write(istdpr,115) cs(2), cs(3), cs(20), cs(21),
     *        cs(25)
          end if
        end if
c
      else if(icasea.eq.4) then
c
c  grand summary, frequency from uncorrected eigenfrequency
c  (note that this is what is given in cs(25)
c  include a test of frequency factor, for the time being,
c  for cases where cgrav .ne. 6.67232e-8
c
        sig=cs(20)
        frq=16666.6666667d0/cs(25)
c
      else if(icasea.eq.5) then
c
c  grand summary, frequency from cs(37) (Richardson extrapolation)
c  if not available, use variational frequency or eigenfrequency
c
        sig=cs(21)
        if(sig.le.0) sig=cs(20)
        if(cs(37).ne.0) then
          frq=1000*cs(37)
        else
          icerr=icerr+1
          if(icerr.le.20.and.istdpr.gt.0) write(istdpr,120)
          if(icerr.le.1) write(istder,120)
          if(cs(27).gt.0) then
            frq=1000*cs(27)
            ierr=-1
          else
            frq=16666.666667d0/cs(25)
            ierr=-1
          end if
        end if
c
      else if(icasea.eq.6) then
c
c  use eigenfrequency, possibly corrected for potential perturbation
c
        sig=cs(21)
        if(sig.le.0) sig=cs(20)
        frq=sqrt(sig/cs(20))*16666.6666667d0/cs(25)
        frq1=1.d6/twopi*sqrt(cgrav*cs(2)*sig/cs(3)**3)
        if(abs(frq/frq1-1).gt.1.e-5.and.icerr.le.20) then
          icerr=icerr+1
          if(istdpr.gt.0) write(istdpr,110) icasea, frq, frq1
          if(icerr.le.2.and.istdpr.gt.0) 
     *      write(istdpr,115) cs(2), cs(3), cs(20), cs(21),
     *      cs(25)
        end if
c
      else
c
        if(istdpr.gt.0) write(istdpr,130) icase
        stop
      end if
c
   90 continue
c
c
      return
c
c  EOF or error
c
   92 ierr=1
      return
c
   95 ierr=2
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Error on reading d/s ',ids,'  in s/r rdfreq'
      return
  100 format(2i5,f8.2)
  110 format(' *** warning. in s/r rdfreq, with icase = ',i1,
     *  ', frq =',1pe13.6,' .ne. frq1 =',e13.6)
  115 format(' cs(2), cs(3) =',1p2e15.6,' cs(20), cs(21), cs(25) =',
     *  3e15.6)
  120 format(' *** error in s/r rdfreq with icase = 5. cs(37) not set')
  130 format(//' ***** error. icase =',i5,' is illegal in rdfreq')
      end
      subroutine rdobsf(icase,ids,l,nord,frq,error,ierr)
c
c  read observed frequency data from d/s ids. case depends on icase.
c  unused.
c
c  icase = 3: no error data. each record consists of
c             l  order  frequency
c
c  icase = 4: file contains error data. each record consists of
c             l  order  frequency  error
c
c  returns mode degree, order, frequency
c  (in microHz) and possibly error in l, nord, frq and error.
c
c  ierr is returned as 1 or 2, on end of file or error in read.
c
c  original version: 18/8/89.
c
c                     ................................
c
c  for the moment skip test for rdofrq format. think about later
c
      implicit double precision (a-h, o-z)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data initrd /0/
c
      ierr = 0
c
c  test for first read, file in the format for rdofrq
c
      if(initrd.eq.0) then
        call skpcom(ids)
        read(ids,*,end=92,err=95) xlobsn
        if(xlobsn.lt.0) then
          read(ids,*,end=92,err=95) xnnobs
        else
          backspace ids
        end if
c
        initrd=1
      end if
c
c  test for case
c
      if(icase.eq.3) then
c
        read(ids,*,end=92,err=95) l,nord,frq
        error=0
      else
        read(ids,*,end=92,err=95) l,nord,frq,error
      end if
c
   90 continue
c
c
      return
c
c  EOF or error
c
   92 ierr=1
      return
c
   95 ierr=2
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Error on reading d/s ',ids,'  in s/r rdobsf'
      return
      end
