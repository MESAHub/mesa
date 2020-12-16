      program main
c
c  set p-mode frequency spacings nu(n,l) - nu(n-delta n, l + delta l),
c  possibly scaled by 3/(2*l+3)
c
c  the type of dataset is determined by the parameter icasin.
c  input file is on d/s 2.
c  output is on file d/s 10.
c
c  original version: 30/6/89
c
c  modified 2/5/91, to allow possibility of reading observations
c  with standard deviations.
c
c  Modified 16/8/93, to allow non-integral step in order, to indicate 
c  interpolation.
c
c  Modified 20/7/97, to filter out f and g modes
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      parameter (nnmax=3000)
      character*280 fin, fout, formcs
      dimension cs(50), formcs(5), ns(nnmax),ls(nnmax), freqs(nnmax),
     *  errfrs(nnmax)
      data formcs /'unformatted','unformatted','formatted',
     *             'unformatted','unformatted'/
c
c  icasin: type of modes used for input.
c  icasin = 1: grand summary.
c  icasin = 2: short summary.
c  icasin = 3: observed frequencies.
c  icasin = 4: grand summary, use (uncorrected) eigenfrequency
c  icasin = 5: grand summary, use Richardson extrapolated
c              eigenfrequency
      icasin=1
c
      write(6,*) 'Enter case number (1 - 5)'
      read(5,*) icasin
c
      write(6,*) 'Enter input file name'
      read(5,'(a)') fin
      write(6,*) 'Enter output file name'
      read(5,'(a)') fout
c
      open(2, file=fin, status='old',
     *  form=formcs(icasin))
      open(10, file=fout, status='unknown')
c
      write(6,*) 'Enter maximum bottom degree'
      read(5,*) lmax
c
      xnstep=1
      lstep=2
      write(6,*) 'enter delta n, delta l (defaults: 1, 2)'
      read(5,*) xnstep, lstep
      iscale=0
      write(6,*) 
     *  'Enter 1 to scale with 3/(2*l+3), 2 to scale with 2/(l+2)'
      read(5,*) iscale
c
      irderr=0
      write(6,*) 'Enter 1 for reading frequency errors'
      write(6,*) '   (only for icasein = 3)'
      read(5,*,end=25) irderr
      write(6,*) 'irderr, icasin',irderr,icasin
      if(icasin.ne.3) irderr=0
c
c  start reading data
c
   25 n=0
c
   30 if(irderr.ne.1) then
        call rdfreq(icasin,2,cs,l,nord,sig,freq,ekin,ierr)
      else
        call rdobsf(4,2,l,nord,freq,errfr,ierr)
      end if
c
      if(ierr.gt.0) go to 35
c
      if(icasin.eq.2.and.cs(1).lt.0) go to 30
c
      if(nord.le.0) go to 30
c
      if(l.gt.lmax+lstep) go to 35
c
      n=n+1
      ls(n)=l
      ns(n)=nord
      freqs(n)=freq
      if(irderr.eq.1) errfrs(n)=errfr
      go to 30
c
c  end of file reached 
c
   35 nn=n
c
c  set parameters for possible interpolation
c
      nstep=xnstep
      fnstep=xnstep-nstep
      gnstep=1-fnstep
c
c  start output
c
      write(10,120) fin, icasin, lstep, xnstep
      if(fnstep.eq.0) then
        write(10,122) nstep, lstep
      else
        write(10,124) fnstep,nstep+1, lstep,gnstep,nstep,lstep
      end if
      if(iscale.eq.1) then
        write(10,128)
      else if(iscale.eq.2) then
        write(10,130)
      end if
      if(irderr.ne.1) then
        write(10,136)
      else
        write(10,138)
      end if
c
c  start stepping through dataset
c
      n1=1
      n2=1
c
   40 continue
c
c..      write(6,*) 'n1, n2, ls(n1), ns(n1), ls(n2), ns(n2):'
c..      write(6,*) n1, n2, ls(n1), ns(n1), ls(n2), ns(n2)
      if(ls(n2).lt.ls(n1)+lstep) then
        n2=n2+1
      else if(ls(n2).eq.ls(n1)+lstep) then
        if(ns(n2).lt.ns(n1)-nstep) then
          n2=n2+1
        else if(ns(n2).eq.ns(n1)-nstep) then
c
c  finally found a match. test for interpolation
c
	  if(fnstep.eq.0) then
            dfreq=freqs(n1)-freqs(n2)
	    isetdf=1
          else
	    if(ns(n2-1).eq.ns(n2)-1) then
c..	      write(6,*) ls(n1),ns(n1),ls(n2-1),ns(n2-1),ls(n2),ns(n2)
	      dfreq=freqs(n1)-fnstep*freqs(n2-1)-gnstep*freqs(n2)
	      isetdf=1
            else
	      isetdf=0
            end if
          end if
c
c  make output
c
	  if(isetdf.eq.1) then
c
            if(iscale.eq.1) then
              dfreq=3.d0*dfreq/dfloat(2*ls(n1)+3)
            else if(iscale.eq.2) then
              dfreq=2.d0*dfreq/dfloat(ls(n1)+2)
            end if
c
c  test for setting error
c
            if(irderr.ne.1) then
              write(10,140) ls(n1),ns(n1), freqs(n1), dfreq
              write(6,140) ls(n1),ns(n1), freqs(n1), dfreq
            else
              errdnl=sqrt(errfrs(n1)*errfrs(n1)+errfrs(n2)*errfrs(n2))
              if(iscale.eq.1) errdnl=errdnl*3.d0/dfloat(2*ls(n1)+3)
              write(10,140) ls(n1),ns(n1), freqs(n1), dfreq, errdnl
              write(6,140) ls(n1),ns(n1), freqs(n1), dfreq, errdnl
            end if
          end if
          n1=n1+1
c
        else if(ns(n2).gt.ns(n1)-nstep) then
          n1=n1+1
c
        end if
c
      else if(ls(n2).gt.ls(n1)+lstep) then
        n1=n1+1
      end if
c
c  test for continuing
c
      if(n1.le.nn.and.n2.le.nn) go to 40
c
      stop
  120 format('# input file: ',a/'#'/
     *  '# data type: ',i3/'#'/
     *  '# step in l =',i3,'   step in order =',f5.2/'#')
  122 format('# show separation nu(n,l) - nu(n -',i2,
     *  ', l +',i2,')')
  124 format('# show separation nu(n,l) - ',f5.2,'*nu(n -',i2,
     *  ', l +',i2,') - ',f5.2,'*nu(n -',i2,', l +',i2,')')
  128 format('# scaled by 3/(2l+3)'/'#')
  130 format('# scaled by 2/(l+2)'/'#')
  136 format('# l, order, nu(n,l), separation:'/'#')
  138 format('# l, order, nu(n,l), separation, error:'/'#')
  140 format(2i4,f11.3,2f11.5)
      end 
