      program main
c
c
c  change between binary and ASCII versions of adiabatic grand summary file
c
c  Original version: 8/10/97
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      character*280 fin, fout
      dimension cs(50),csp(8),ics(8),icsp(8)
      equivalence (cs(39),ics(1))
c
      data csp /8*-1.d0/
      data icsp /8*-1/
      data lp /-1/
c
      write(6,*) 'Enter 1 for transformation from binary to ASCII'
      write(6,*) '      2 for transformation from ASCII to binary'
      read(5,*) icase
c
c  set file names
c
      write(6,*) 'Enter input file name'
      read(5,'(a)',end=90,err=90) fin
      write(6,*) 'Enter output file name'
      read(5,'(a)',end=90,err=90) fout
c
      if(icase.eq.1) then
        open(2,file=fin,status='old',form='unformatted')
        open(3,file=fout,status='unknown')
      else
        open(2,file=fin,status='old')
        open(3,file=fout,status='unknown',form='unformatted')
      end if
      write(6,105) ' Transform file ',fin
c
c  set flag for only first mode at each l
c
      ibrief=0
      write(6,*) 'Enter 1 for only first mode at each l'
      read(5,*,end=5,err=90) ibrief
c
    5 n=0
c
   10 if(icase.eq.1) then
	read(2,end=20,err=20) cs
      else
	read(2,107,end=20,err=20) (cs(i),i=1,38)
	read(2,108,end=20,err=20) (ics(i),i=1,8)
      end if
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
        write(6,110) (cs(k),k=2,5),ics(2),ics(3)
	if(cs(12).gt.0) write(6,112) cs(12)
      end if
c
c  printing depends on whether l is integer or not, and whether
c  or not cs contains Richardson extrapolated frequency
c
      l=nint(cs(18))
c
      if(ibrief.ne.1.or.l.ne.lp) then
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
          if(inew.eq.1) write(6,115)
          if(abs(cs(18)-l).lt.1.e-5) then
            write(6,120) n,l,nord,sig,cs(24),frq,ics(5)
          else
            write(6,122) n,cs(18),nord,sig,cs(24),frq,ics(5)
          end if
        else
          if(inew.eq.1) write(6,117)
          if(abs(cs(18)-l).lt.1.e-5) then
            write(6,130) n,l,nord,sig,cs(24),frq,cs(37),ics(5)
          else
            write(6,132) n,cs(18),nord,sig,cs(24),frq,cs(37),ics(5)
          end if
        end if
c
      end if
c
      lp = l
c
c  output record
c
      if(icase.eq.1) then
	write(3,107) (cs(i),i=1,38)
	write(3,108) (ics(i),i=1,8)
      else
	write(3) cs
      end if
c
      go to 10
c
   20 continue
c
   90 continue
      stop
  105 format(//a,a)
  107 format(1p4e20.13)
  108 format(8i10)
  110 format(///' data (1-4):',1p4e13.5/'       nn  :',i5,
     *  '       mdintg  :',i3)
  112 format(/' Model truncated at x =',f10.5)
  115 format(//' n, l, order, sigma**2, E, frequency (mHz), case:'/)
  117 format(//' n, l, order, sigma**2, E, v. frequency (mHz),',
     *         ' Ri. frequency (mHz), case:'/)
  120 format(i5,i6,4x,i5,f14.5,1pe13.5,0pf10.6,i14)
  122 format(i5,f10.3,i5,f14.5,1pe13.5,0pf10.6,i14)
  130 format(i5,i6,4x,i5,f14.5,1pe13.5,0p2f9.6,i14)
  132 format(i5,f10.3,i5,f14.5,1pe13.5,0p2f9.6,i14)
      end
