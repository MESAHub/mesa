      subroutine wrfreq(icase,ids,cs,l,nord,sig,frq,ekin,ierr)
c
c  writes frequency data to d/s ids. case depends on icase.
c
c  icase = 1: grand summary. from cs(1-50)
c  icase = 2: short summary. from cs(1-7)
c  icase = 3: observed frequencies, from l, nord, frq
c
c  when icase = 1 or 2, returns mode degree, order, 
c  dimensionless sigma**2, frequency (in microHz) and 
c  dimensionless energy in l, nord, sig, frq and ekin.
c
c  Note: for icase = 1 returns corrected sigma**2, if frequency
c  calculation was in Cowling approximation. This is normally
c  also what is set in the short summary.
c
c  ierr is currently not used, but is retained for compatibility
c  with s/r rdfreq.
c
c  original version: 24/8/87.
c
c                     ................................
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      dimension cs(*)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      ierr = 0
c
c  test for case
c
      go to (10, 20, 30), icase
c
c  grand summary
c
   10 write(ids) (cs(i),i=1,50)
c
      l=nint(cs(18))
      nord=nint(cs(19))
      sig=cs(21)
      if(sig.le.0) sig=cs(20)
      frq=1000*cs(27)
      if(frq.le.0) frq=16666.667/cs(25)
      ekin=cs(24)
c
      go to 40
c
c  short summary
c
   20 write(ids) (cs(i),i=1,7)
      l=nint(cs(1))
      nord=nint(cs(2))
      sig=cs(3)
      frq=1000*cs(5)
      ekin=cs(4)
c
      go to 40
c
c  observed frequencies
c
   30 write(ids,100) l,nord,frq
      sig=0
      ekin=0
c
   40 continue
c
c
      return
c
c  EOF or error (currently not used)
c
   50 ierr=1
      return
c
   60 ierr=2
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Error on reading d/s ',ids,'  in s/r rdfreq'
      return
  100 format(2i5,f8.2)
      end
      subroutine wrobsf(icase,ids,l,nord,frq,error,ierr)
c
c  writes frequency data to d/s ids. case depends on icase.
c
c  icase = 3: observed data, without errors.
c  icase = 4: observed data, with errors.
c
c  ierr is currently not used, but is retained for compatibility
c  with s/r rdfreq.
c
c  original version: 18/8/89.
c
c                     ................................
c
      implicit double precision (a-h, o-z)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      ierr = 0
c
c  test for case
c
c  observed frequencies, without errors
c
      if(icase.eq.3) then
        write(ids,100) l,nord,frq
      else
        write(ids,110) l,nord,frq,error
      end if
c
      return
c
c  EOF or error (currently not used)
c
   50 ierr=1
      return
c
   60 ierr=2
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Error on reading d/s ',ids,'  in s/r rdfreq'
      return
  100 format(2i5,f8.2)
  110 format(2i5,2f10.4)
      end
