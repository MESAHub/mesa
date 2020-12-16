      program main
c
c  brief summary of oscillation results
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
c  Modified 6/11/99, to allow separation of non-integer degrees
c  (flagged by ibrief = 2)
c
c  Modified 1/11/02, providing output of azimuthal order, for
c  results with rotation
c
c  Modified 23/5/08, stopping execution if istdpr .le. 0
c  (note that this would make no sense!)
c
      implicit double precision (a-h, o-z)
      character*280 fin
      logical newl
      dimension cs(50),csp(8),ics(8),icsp(8)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      equivalence (cs(39),ics(1))
c
      data csp /8*-1./
      data icsp /8*-1/
      data lp /-1/
      data elp /-1.d0/
c
c  set input file name
c
      if(istdpr.le.0) then
	write(istder,'(/'' ***** Calling scan-agsm.d with istdpr ='',
     *    i2)') istdpr
	stop
      end if
c
      write(istdpr,*) 'Enter input file name'
      read(5,'(a)',end=90,err=90) fin
      open(2,file=fin,status='old',form='unformatted')
      write(istdpr,105) ' Scan of file ',fin
c
c  set flag for only first mode at each l
c
      ibrief=0
      write(istdpr,*) 'Enter 1 or 2 for only first mode at each l'
      read(5,*,end=5,err=90) ibrief
c
    5 n=0
c
   10 read(2,end=20,err=20) cs
c
      n=n+1
c
c  test for new model
c
      inew=0
      do 12 k=2,5
      if(cs(k).eq.0) then
        if(k.le.3) inew=-1
      else if(inew.eq.0.and.abs(csp(k)/cs(k)-1).gt.1.e-5) then
        inew=1
      end if
c
      if(ics(2).ne.icsp(2).or.ics(3).ne.icsp(3)) inew=1
   12 continue
c
c  test for error in read
c
      if(inew.eq.-1) go to 20
c
      do 13 k=1,8
      icsp(k)=ics(k)
   13 csp(k)=cs(k)
c
      if(inew.eq.1) then
        write(istdpr,110) (cs(k),k=2,5),ics(2),ics(3)
	if(cs(12).gt.0.and.istdpr.gt.0) write(istdpr,111) cs(12)
      end if
c
c  printing depends on whether l is integer or not, and whether
c  or not cs contains Richardson extrapolated frequency
c
      el=cs(18)
      l=nint(cs(18))
      icase=ics(5)
      if(ibrief.eq.1) then
	newl = l.ne.lp
      else if(ibrief.eq.2) then
	newl = abs(el-elp).gt.1.e-6
      end if
c
c  set flag for inclusion of rotation
c
      irotsl=mod(icase/100,10)
      if(irotsl.eq.1) then
	em=cs(38)
	m=nint(em)
      end if
c
c  test for full output (for ibrief = -1)
c
      if(ibrief.eq.-1) then 
        write(istdpr,112) n
	do 15 i0=1,36,5
	imax=min0(37,i0+4)
   15   write(istdpr,113) i0,imax,(cs(i),i=i0,imax)
        write(istdpr,114) (ics(i),i=1,7)
      else if(ibrief.eq.0.or.newl) then
        nord=nint(cs(19))
        sig=cs(21)
        if(sig.le.0) sig=cs(20)
c
c  set frequency. variational, or, if not available, (possibly
c  Cowling-corrected) eigenfrequency
c
        frq=cs(27)
        if(frq.le.0) frq=sqrt(sig/cs(20))*16.66667/cs(25)
        if(cs(37).eq.0) then
          if(inew.eq.1) then
            if(irotsl.ne.1) then
              write(istdpr,116)
            else
              write(istdpr,117)
	    end if
	  end if
          if(abs(cs(18)-l).lt.1.e-5) then
	    if(irotsl.ne.1) then
              write(istdpr,120) n,l,nord,sig,cs(24),frq,icase
            else
              write(istdpr,121) n,l,m,nord,sig,cs(24),frq,icase
	    end if
          else
            write(istdpr,122) n,cs(18),nord,sig,cs(24),frq,icase
          end if
        else
          if(inew.eq.1) then
            if(irotsl.ne.1) then
              write(istdpr,118)
            else
              write(istdpr,119)
	    end if
	  end if
          if(abs(cs(18)-l).lt.1.e-5) then
	    if(irotsl.ne.1) then
              write(istdpr,130) n,l,nord,sig,cs(24),frq,cs(37),icase
            else
              write(istdpr,131) n,l,m,nord,sig,cs(24),frq,cs(37),icase
	    end if
          else
            write(istdpr,132) n,cs(18),nord,sig,cs(24),frq,cs(37),icase
          end if
        end if
c
      end if
c
      lp = l
      elp = el
c
      go to 10
c
   20 continue
c
   90 continue
      stop
  105 format(//a,a)
  110 format(///' data (1-4):',1p4e13.5/'       nn  :',i5,
     *  '       mdintg  :',i3)
  111 format(/' Model truncated at x =',f10.5)
  112 format(/' Mode no.',i5)
  113 format(' cs(',i2,'-',i2,'):',1p5e13.5)
  114 format(' ics:',8i8)
  116 format(//' n, l, order, sigma**2, E, frequency (mHz), case:'/)
  117 format(//' n, l, m, order, sigma**2, E, frequency (mHz), case:'/)
  118 format(//' n, l, order, sigma**2, E, v. frequency (mHz),',
     *         ' Ri. frequency (mHz), case:'/)
  119 format(//' n, l, m, order, sigma**2, E, v. frequency (mHz),',
     *         ' Ri. frequency (mHz), case:'/)
  120 format(i5,i6,3x,i6,f14.5,1pe13.5,0pf10.6,i14)
  121 format(i5,2i3,3x,i6,f14.5,1pe13.5,0pf10.6,i14)
  122 format(i5,f10.3,i5,f14.5,1pe13.5,0pf10.6,i14)
  130 format(i5,i6,3x,i6,f14.5,1pe13.5,0p2f9.6,i14)
  131 format(i5,2i3,3x,i5,f14.5,1pe13.5,0p2f9.6,i14)
  132 format(i5,f10.3,i5,f14.5,1pe13.5,0p2f9.6,i14)
      end
