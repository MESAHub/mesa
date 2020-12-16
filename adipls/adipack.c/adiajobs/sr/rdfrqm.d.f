      subroutine rdfrqm(icase,ids,nmodel,cs,l,nord,sig,frq,ekin,ierr)
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
c
c  for icase = 1, 2, 4 or 5 go through file, to read results for model
c  no nmodel, if this is .gt. 1
c  if nmodel .ge. 1, read only data for one model
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
c  original version: 26/9/86.
c
c                     ................................
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 17/5/90
c
      implicit double precision (a-h, o-z)
      logical newmod
      dimension cs(*), csmod(4), csmodp(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      save
      data csmod /4*-1./
      data csmodp /4*-1./
      data nmod /0/
c
c  for the moment skip test for rdofrq format. think about later
c
      data initrd /0/
      data icerr /0/
c
      ierr = 0
c
c  test for case
c
    5 go to (10, 20, 30, 40, 50), icase
c
c  grand summary
c
   10 read(ids,end=92,err=95) (cs(i),i=1,50)
c
c  set csmod
c
      do 12 i=1,4
   12 csmod(i)=cs(i+1)
c
      l=nint(cs(18))
      nord=nint(cs(19))
      sig=cs(21)
      if(sig.le.0) sig=cs(20)
      frq=1000*cs(27)
      if(frq.le.0) then
        frq=sqrt(sig/cs(20))*16666.66667/cs(25)
        frq1=1.5915494e5*sqrt(6.6732e-8*cs(2)*sig/cs(3)**3)
        if(abs(frq/frq1-1).gt.1.e-5.and.icerr.le.20) then
          icerr=icerr+1
          if(istdpr.gt.0) write(istdpr,110) frq, frq1
          if(icerr.le.2.and.istdpr.gt.0) 
     *      write(istdpr,115) cs(2), cs(3), cs(20), cs(21), cs(25)
        end if
      end if
      ekin=cs(24)
c
      go to 80
c
c  short summary
c
   20 read(ids,end=92,err=95) (cs(i),i=1,7)
c
      if(cs(1).lt.0) then
c
c  set csmod
c
        do 22 i=1,4
   22   csmod(i)=cs(i+1)
c
      end if
c
      l=nint(cs(1))
      nord=nint(cs(2))
      sig=cs(3)
      frq=1000*cs(5)
      ekin=cs(4)
c
      go to 80
c
c  observed frequencies
c
c  test for first read, file in the format for rdofrq
c
   30 continue
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
      read(ids,*,end=92,err=95) l,nord,frq
      sig=0
      ekin=0
c
      go to 80
c
c  grand summary, frequency from uncorrected eigenfrequency
c  (note that this is what is given in cs(25)
c  include a test of frequency factor, for the time being,
c  for cases where cgrav .ne. 6.6732e-8
c
   40 read(ids,end=92,err=95) (cs(i),i=1,50)
c
c  set csmod
c
      do 42 i=1,4
   42 csmod(i)=cs(i+1)
c
      l=nint(cs(18))
      nord=nint(cs(19))
      sig=cs(20)
      frq=16666.66667/cs(25)
      ekin=cs(24)
c
      go to 80
c
c  grand summary, frequency from cs(37) (Richardson extrapolation)
c  if not available, use variational frequency or eigenfrequency
c
   50 read(ids,end=92,err=95) (cs(i),i=1,50)
c
c  set csmod
c
      do 52 i=1,4
   52 csmod(i)=cs(i+1)
c
      l=nint(cs(18))
      nord=nint(cs(19))
      sig=cs(21)
      if(sig.le.0) sig=cs(20)
      if(cs(37).ne.0) then
        frq=1000*cs(37)
      else
        icerr=icerr+1
        if(icerr.le.20.and.istdpr.gt.0) write(istdpr,120)
        if(cs(27).gt.0) then
          frq=1000*cs(27)
        else
          frq=16666.667/cs(25)
        end if
      end if
      ekin=cs(24)
c
   80 continue
c
c  test for new model
c
      if(icase.ne.3.and.nmodel.ge.1) then
        do 82 i=1,4
        newmod=abs(csmodp(i)/csmod(i)-1).gt.1.e-4
   82   csmodp(i)=csmod(i)
        if(newmod) then
          nmod=nmod+1
	  if(istdpr.gt.0) 
     *      write(istdpr,*) 'model no',nmod,' csmod:',csmod
	  if(istdpr.gt.0) write(istdpr,*) 'nmodel =',nmodel
          if(nmod.gt.nmodel) then
            go to 92
          end if
        end if
c  test on model number
        if(nmod.lt.nmodel) then
          go to 5
        end if
      end if
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
  110 format(' *** warning. in s/r rdfreq, with icase = 4, frq =',
     *  1pe13.6,' .ne. frq1 =',e13.6)
  115 format(' cs(2), cs(3) =',1p2e15.6,' cs(20), cs(21), cs(25) =',
     *  3e15.6)
  120 format(' *** error in s/r rdfreq with icase = 5. cs(37) not set')
      end
