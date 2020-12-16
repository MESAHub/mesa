      program main
c
c  merge ordered files of adiabatic grand summaries, short summaries
c  or observed frequencies.
c  the type of dataset is determined by the parameter icasin.
c  input files are on d/s 2 and 3.
c  output file is on d/s 10.
c
c  note: if the same mode is present on both files, the mode
c  from the  f i r s t  file is used.
c
c  if the same mode appears several times in one of the files,
c  only the  f i r s t  occurence is included.
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      character*280 fin1, fin2, fout, formcs
      dimension cs1(50),cs2(50), formcs(4)
c$nl      namelist /exec/ icasin,idiag,icmod
      data formcs /'unformatted','unformatted','formatted','formatted'/
c
c  icasin: type of modes used for input or selection.
c  icasin = 1: grand summary.
c  icasin = 2: short summary.
c  icasin = 3: observed frequencies, without errors.
c  icasin = 4: observed frequencies, with errors.
      icasin=1
c  idiag: if idiag = 1 print details of modes output
      idiag=1
c  icmod: if several models on files, order after equilibrium model
c     quantity determined by icmod, as
c     icmod = 1: mass
c     icmod = 2: radius
c     icmod = 3: central pressure
c     icmod = 4: central density
      icmod=3
c
c$nl      read(5,exec,end=10,err=10)
c
      write(6,*) 'Enter case number (1 - 4)'
      read(5,*) icasin
c
      write(6,*) 'Enter first input file name'
      read(5,'(a)') fin1
      write(6,*) 'Enter second input file name'
      read(5,'(a)') fin2
      write(6,*) 'Enter output file name'
      read(5,'(a)') fout
c
      open(2, file=fin1, status='old',
     *  form=formcs(icasin))
      open(3, file=fin2, status='old',
     *  form=formcs(icasin))
      open(10, file=fout, status='unknown',
     *  form=formcs(icasin))
c
      write(6,*) 'idiag,icmod?'
      write(6,*) idiag, icmod
      read(5,*,end=10,err=10) idiag, icmod
c
      write(6,*) 'Input files:'
      write(6,*) fin1
      write(6,*) fin2
      write(6,*) 'Output file:'
      write(6,*) fout
c
   10 if(idiag.gt.0) write(6,105)
c
c  start stepping through dataset
c
      nmd1=0
      nmd2=0
      nmd=1
      iend1=0
      iend2=0
c
      xlp1=-1
      nordp1=-1
      csmp1=-1
c
      xlp2=-1
      nordp2=-1
      csmp2=-1
c
      irdmd2=1
c
      idgmod=1
c
   30 if(icasin.le.3) then
        call rdfreq(icasin,2,cs1,l1,nord1,sig1,frq1,ekin1,ierr1)
      else
        call rdobsf(icasin,2,l1,nord1,frq1,error1,ierr1)
      end if
c..      write(6,*) 'read l1, nord1 =',l1, nord1
c
      if(ierr1.gt.0) go to 60
c
c set model quantity
c
      if(icasin.eq.1) then
        csm1=cs1(icmod+1)
      else if(icasin.eq.2) then
	if(l1.lt.0) csm1=cs1(icmod+2)
      else
	csm1=1
      end if
c
      if(csm1.eq.0) go to 60
c
      nmd1=nmd1+1
c
c  set degree, depending on icasin
c
      if(icasin.eq.1) then
        xl1=cs1(18)
      else if(icasin.eq.2) then
        xl1=cs1(1)
      else
        xl1=dfloat(l1)
      end if
c
c  test for the same mode
c
      if(abs(xl1-xlp1).lt.1.d-5.and.csm1.eq.csmp1) then
        if(nord1.eq.nordp1) then
          go to 30
        else if (nord1.lt.nordp1) then
c
c  diagnostics for non-monotonic order
c
          iwrd=1
          write(6,102) iwrd,nmd1,xl1,nord1
          go to 30
        end if
      end if
c
      xlp1=xl1
      nordp1=nord1
      csmp1=csm1
c
c
      if(iend2.eq.1) go to 50 
      if(irdmd2.eq.0) go to 40
c
   35 if(icasin.le.3) then
        call rdfreq(icasin,3,cs2,l2,nord2,sig2,frq2,ekin2,ierr2)
      else
        call rdobsf(icasin,3,l2,nord2,frq2,error2,ierr2)
      end if
c..      write(6,*) 'read l2, nord2 =',l2, nord2, 'iend1 =',iend1
c
      if(ierr2.gt.0) go to 65
c
c set model quantity
c
      if(icasin.eq.1) then
        csm2=cs2(icmod+1)
      else if(icasin.eq.2) then
	if(l2.lt.0) csm2=cs2(icmod+2)
      else
	csm2=1
      end if
c
      if(csm2.eq.0) go to 65
c
      nmd2=nmd2+1
c
c  set degree, depending on icasin
c
      if(icasin.eq.1) then
        xl2=cs2(18)
      else if(icasin.eq.2) then
        xl2=cs2(1)
      else
        xl2=dfloat(l2)
      end if
c
c  test for the same mode
c
      if(abs(xl2-xlp2).lt.1.d-5.and.csm2.eq.csmp2) then
        if(nord2.eq.nordp2) then
          go to 35
        else if (nord2.lt.nordp2) then
c
c  diagnostics for non-monotonic order
c
          iwrd=2
          write(6,102) iwrd,nmd2,xl2,nord2
          go to 35
        end if
      end if
c
      xlp2=xl2
      nordp2=nord2
      csmp2=csm2
c
      irdmd2=0
c
c..      write(6,*) 'test on iend1 =',iend1
c
      if(iend1.eq.1) go to 52 
c
c  test on model parameter
c
   40 dmdl=csm1/csm2-1
c..      write(6,*) 'test on dmdl2 =',dmdl
      if(abs(dmdl).lt.1.d-6) then
        idgmod=1
        go to 41
      else 
        if(idgmod.eq.1) then
          write(6,115) csm1, csm2
          write(6,105)
          idgmod=0
        end if
        if(dmdl.lt.0) then
          go to 50
        else
          go to 52
        end if
      end if
c
      if(dmdl) 50,41,52
c
c  same model, test on l
c
   41 del=xl1-xl2
      if(abs(del).lt.1.d-6) del=0
c..      write(6,*) 'test on xl1 - xl2 =',del
c
      if(del) 50,42,52
c
c  same l, test on order
c
   42 continue
c..      write(6,*) 'at 42, l1, n1, n2, l2',xl1,nord1,xl2,nord2
      if(nord1-nord2) 50,45,52
c
c  same mode. output mode 1 and set flag for reading new mode 2
c
   45 irdmd2=1
c
c  output mode 1
c
   50 if(icasin.le.3) then
        call wrfreq(icasin,10,cs1,l1,nord1,sig1,frq1,ekin1,ierr1)
      else
        call wrobsf(icasin,10,l1,nord1,frq1,error1,ierr1)
      end if
      iwr=1
      if(idiag.gt.0) write(6,110) nmd,iwr,nmd1,xl1,nord1,frq1
      go to 55
c
c  output mode 2
c
   52 if(icasin.le.3) then
        call wrfreq(icasin,10,cs2,l2,nord2,sig2,frq2,ekin2,ierr2)
      else
        call wrobsf(icasin,10,l2,nord2,frq2,error2,ierr2)
      end if
c..      write(6,*) 'write 2. l1, n1, l2, n2 =',l1, nord1, l2, nord2
      iwr=2
      if(idiag.gt.0) write(6,110) nmd,iwr,nmd2,xl2,nord2,frq2
c
   55 nmd=nmd+1
c
c  determine which mode to read next
c
      if(iend1.eq.1) go to 35 
      if(iend2.eq.1) go to 30 
      go to (30,35), iwr
c
c  end of file reached on one of the datasets
c
   60 iend1=1
      go to 70
c
   65 iend2=1
c
c  test for stopping
c
   70 if(iend1.eq.0) go to 50 
c
      if(iend2.eq.0) then
        if(irdmd2) 52,52,35
      end if
c
   90 continue
c
      write(6,120) nmd,nmd1,nmd2
c
      stop
  102 format(//' ***** non-monotonic order on d/s',i3/
     *  7x,'mode no.',i5,' degree, order =',f6.1,i5)
  105 format(///' n, data set, ni, l, order, frequency',
     *  ' (millihz)'/)
  110 format(3i5,f8.2,i5,f12.5)
  115 format(//' *****  Warning. different models.'/
     *         '        csm1 =',1pe15.7/
     *         '        csm2 =',1pe15.7)
  120 format(//' number of modes output =',i5/
     *  ' number of modes from file 1 =',i5/
     *  ' number of modes from file 2 =',i5)
      end 
