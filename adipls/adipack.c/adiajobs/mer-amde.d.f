      program main
c
c  merge ordered files of adiabatic eigenfunctions or rotationally
c  symmetric kernels.
c  the type of dataset is determined by the parameter icase.
c  input files are on d/s 2 and 3.
c  output file is on d/s 10.
c
c  note: if the same mode is present on both files, the mode
c  from the  f i r s t  file is used.
c
c  if the same mode appears several times in one of the files,
c  only the  f i r s t  occurence is included.
c
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 29/4/91
c
      implicit double precision (a-h, o-z)
      parameter (nnmax=2501)
      character*280 fin1, fin2, fout
      dimension cs1(50),cs2(50),x1(nnmax),x2(nnmax),y1(6,nnmax),
     *  y2(6,nnmax)
c$nl      namelist /exec/ icase,idiag,icmod
c
c  icase: type of eigenfunction.
c     icase = 1: Full eigenfunction, as output from adipls
c     icase = 2: Reduced eigenfunction. First record gives nn and x,
c                subsequenct records give the eigenfunctions, in terms
c                of two variables.
c     icase = 3: Kernels for spherically symmetric rotation.
      icase = 1
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
      write(6,*) 'Enter case number (1 - 3)'
      read(5,*) icase
c
      write(6,*) 'Enter first input file name'
      read(5,'(a)') fin1
      write(6,*) 'Enter second input file name'
      read(5,'(a)') fin2
      write(6,*) 'Enter output file name'
      read(5,'(a)') fout
c
      open(2, file=fin1, status='old',
     *  form='unformatted')
      open(3, file=fin2, status='old',
     *  form='unformatted')
      open(10, file=fout, status='unknown',
     *  form='unformatted')
c
      write(6,*) 'idiag,icmod?'
      write(6,*) idiag, icmod
      read(5,*,end=10) idiag, icmod
c
      write(6,*) 'Input files:'
      write(6,*) fin1
      write(6,*) fin2
      write(6,*) 'Output file:'
      write(6,*) fout
c
   10 if(idiag.gt.0) write(6,105)
c
c  for icase = 2, input, and possibly output, record with nn and x
c
      if(icase.eq.2) then
        read(2) nn1,(x1(n),n=1,nn1)
        read(3) nn2,(x2(n),n=1,nn2)
        write(10) nn1,(x1(n),n=1,nn1)
        ivar=2
      else
        ivar=6
      end if
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
   30 if(icase.ne.2) then
        read(2,end=60) cs1,nn1,(x1(n),(y1(i,n),i=1,ivar),n=1,nn1)
      else
        read(2,end=60) cs1,((y1(i,n),i=1,2),n=1,nn1)
      end if
c
      l1=nint(cs1(18))
      nord1=nint(cs1(19))
c
      if(cs1(27).gt.0) then
        frq1 = 1000*cs1(27)
      else
        frq1 = 16666.666667d0/cs1(25)
      end if
c
c set model quantity
c
      csm1=cs1(icmod+1)
c
      if(csm1.eq.0) go to 60
c
      nmd1=nmd1+1
c
c  set degree
c
      xl1=cs1(18)
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
   35 if(icase.ne.2) then
        read(3,end=65) cs2,nn2,(x2(n),(y2(i,n),i=1,ivar),n=1,nn2)
      else
        read(3,end=65) cs2,((y2(i,n),i=1,2),n=1,nn2)
      end if
c
      l2=nint(cs2(18))
      nord2=nint(cs2(19))
c
      if(cs2(27).gt.0) then
        frq2 = 1000*cs2(27)
      else
        frq2 = 16666.666667d0/cs2(25)
      end if
c
c set model quantity
c
      csm2=cs2(icmod+1)
c
      if(csm2.eq.0) go to 65
c
      nmd2=nmd2+1
c
c  set degree
c
      xl2=cs2(18)
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
      if(iend1.eq.1) go to 52 
c
c  test on model parameter
c
   40 dmdl=csm1/csm2-1
      if(abs(dmdl).lt.1.d-6) dmdl=0
c
      if(dmdl) 50,41,52
c
c  same model, test on l
c
   41 del=xl1-xl2
      if(abs(del).lt.1.d-6) del=0
c
      if(del) 50,42,52
c
c  same l, test on order
c
   42 if(nord1-nord2) 50,45,52
c
c  same mode. output mode 1 and set flag for reading new mode 2
c
   45 irdmd2=1
c
c  output mode 1
c
   50 if(icase.ne.2) then
        write(10) cs1,nn1,(x1(n),(y1(i,n),i=1,ivar),n=1,nn1)
      else
        write(10) cs1,((y1(i,n),i=1,2),n=1,nn1)
      end if
c
      iwr=1
      if(idiag.gt.0) write(6,110) nmd,iwr,nmd1,xl1,nord1,frq1
      go to 55
c
c  output mode 2
c
   52 if(icase.ne.2) then
        write(10) cs2,nn2,(x2(n),(y2(i,n),i=1,ivar),n=1,nn2)
      else
        write(10) cs2,((y2(i,n),i=1,2),n=1,nn2)
      end if
c
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
  120 format(//' number of modes output =',i5/
     *  ' number of modes from file 1 =',i5/
     *  ' number of modes from file 2 =',i5)
      end 
