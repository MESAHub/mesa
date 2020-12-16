      program main
c
c  set file corresponding to observed frequencies from grand summary
c  or short summary.
c  the type of dataset is determined by the parameter icasin.
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      character*280 fin, fout
      character head*40, ccase*5
      logical noprm
      dimension cs(50), head(4), ccase(6), ics(8)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      equivalence(cs(39),ics(1))
      data ccase /'grand', 'short', 'dummy', 'grand', 'grand', 'grand'/
      data head /'Data set from','      summary in file',' ',' '/
c
c  icasin: defines type of modes used for input, and controls output.
c  icasin = icase0 + 10*icase1 + 100*icase2
c  icase0 = 1: grand summary, variational frequency.
c  icase0 = 2: short summary.
c  icase0 = 4: grand summary, from eigenfrequency in cs(20).
c              Note that this allows setting Cowling 
c              approximation frequency.
c  icase0 = 5: grand summary, from Richardson extrapolation frequency 
c              in cs(37), if this is set. 
c              Otherwise variational frequency is used.
c  icase0 = 6: grand summary, from (possibly corrected) eigenfrequency 
c              in cs(21).
c 
c  If icase1 gt 0 include mode energy
c  If icase2 gt 0, and rotational solution is included, print 
c  m values separately, otherwise print only mode with m = 0.
c
      icasin=1
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter case number (1, 2, 4, 5 or 6)'
      read(5,*) icasin
      icase0=mod(icasin,10)
      icase1=mod(icasin/10,10)
      icase2=mod(icasin/100,10)
c..      write(6,*) 'icase0, icase1, icase2:', icase0, icase1, icase2
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter input file name'
      read(5,'(a)') fin
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter output file name'
      read(5,'(a)') fout
c
      open(2, file=fin, status='old',
     *  form='unformatted')
      open(10, file=fout, status='unknown')
      if(istdpr.gt.0) 
     *  write(istdpr,*) 
     *  'Enter 1 for output format f8.2, 2 for format f10.4,',
     *  ' 3 for format f12.6'
      iformt = 1
      read(5,*) iformt
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Input file:'
      if(istdpr.gt.0) 
     *  write(istdpr,*) fin
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Output file:'
      if(istdpr.gt.0) 
     *  write(istdpr,*) fout
c
c  set header
c
      head(2)(1:5) = ccase(icase0)
      write(head(3),110) fin
      if(icase0.eq.1) then
        head(4) = 'using variational frequencies'
      else if(icase0.eq.4) then
        head(4) = 'using (uncorrected) eigenfrequencies'
      else if(icase0.eq.5) then
        head(4) = 'using Richardson frequencies'
      else if(icase0.eq.6) then
        head(4) = 'using (corrected) eigenfrequencies'
      end if
c
c  start stepping through dataset
c
      nrd=0
      ierrri=0
c
   30 call rdfreq(icase0,2,cs,l,nord,sig,frq,ekin,ierr)
c
      if(ierr.gt.0) go to 60
c
c skip model quantity
c
      if(icase0.eq.2.and.l.lt.0) then
        go to 30
      end if
c
c  test for unavailability of Richardson frequencies
c
      if(ierrri.eq.0.and.ierr.eq.-1) then
	ierrri=1
	if(cs(27).gt.0) then
          head(4) = 'using variational frequencies'
        else
          head(4) = 'using (uncorrected) eigenfrequencies'
        end if
      end if
c
c  test for inclusion of individual m values
c
      noprm=.true.
c
      if(icase0.ne.2) then
        icase=ics(5)
        irotsl=mod(icase/100,10)
c..	write(6,*) 'icase, irotsl',  icase, irotsl
c
c  with rotational splitting, and icase2 = 0, only print m = 0 modes
c
        if(irotsl.eq.1) then
	  if(icase2.eq.1) then
	    noprm=.false.
	    m=nint(cs(38))
	  else if(cs(38).ne.0) then
	    go to 30
	  end if
        end if
      end if
c
      nrd=nrd+1
c
      if(iformt.eq.1) then
        if(nrd.le.4) then
	  if(noprm) then
	    if(icase1.eq.0) then
              write(10,120) l,nord,frq,head(nrd)
            else
              write(10,125) l,nord,frq,ekin,head(nrd)
            end if
          else
	    if(icase1.eq.0) then
              write(10,126) l,m,nord,frq,head(nrd)
            else
              write(10,127) l,m,nord,frq,ekin,head(nrd)
            end if
	  end if
        else
	  if(noprm) then
	    if(icase1.eq.0) then
              write(10,120) l,nord,frq
            else
              write(10,125) l,nord,frq,ekin
            end if
	  else
	    if(icase1.eq.0) then
              write(10,126) l,m,nord,frq
            else
              write(10,127) l,m,nord,frq,ekin
            end if
          end if
        end if
      else if(iformt.eq.2) then
        if(nrd.le.4) then
	  if(noprm) then
	    if(icase1.eq.0) then
              write(10,130) l,nord,frq,head(nrd)
            else
              write(10,135) l,nord,frq,ekin,head(nrd)
            end if
	  else
	    if(icase1.eq.0) then
              write(10,136) l,m,nord,frq,head(nrd)
            else
              write(10,137) l,m,nord,frq,ekin,head(nrd)
            end if
          end if
        else
	  if(noprm) then
	    if(icase1.eq.0) then
              write(10,130) l,nord,frq
            else
              write(10,135) l,nord,frq,ekin
            end if
	  else
	    if(icase1.eq.0) then
              write(10,136) l,m,nord,frq
            else
              write(10,137) l,m,nord,frq,ekin
            end if
          end if
        end if
      else 
        if(nrd.le.4) then
	  if(icase1.eq.0) then
            write(10,140) l,nord,frq,head(nrd)
          else
            write(10,145) l,nord,frq,ekin,head(nrd)
          end if
        else
	  if(icase1.eq.0) then
            write(10,140) l,nord,frq
          else
            write(10,145) l,nord,frq,ekin
          end if
        end if
      end if
      go to 30
c
   60 continue
      stop
  110 format(a40)
  120 format(i5,i7, f8.2,2x,a40)
  125 format(i5,i7, f8.2,1pe13.5,2x,a40)
  126 format(3i4,   f8.2,2x,a40)
  127 format(3i4,   f8.2,1pe13.5,2x,a40)
  130 format(i5,i7,f10.4,2x,a40)
  135 format(i5,i7,f10.4,1pe13.5,2x,a40)
  136 format(3i4,  f10.4,2x,a40)
  137 format(3i4,  f10.4,1pe13.5,2x,a40)
  140 format(i5,i7,f12.6,2x,a40)
  145 format(i5,i7,f12.6,1pe13.5,2x,a40)
      end
