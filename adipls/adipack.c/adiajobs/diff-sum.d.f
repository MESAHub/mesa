      program main
c
c  finds modes present in one but not the other of the two files of
c  adiabatic grand or short summaries, or observed frequencies,
c  on unit 2 and 3. 
c  non-overlapping modes are printed. 
c  furthermore a corresponding summary of the modes on unit 2
c  and not on unit 3 are output on unit 12, 
c  and a corresponding summary of the modes present on unit 3 
c  but not on unit 2 are output on unit 13.
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 29/3/90
c
      implicit double precision (a-h, o-z)
      character*280 fin1, fin2, fout1, fout2, form
      dimension cs1(50),cs2(50), form(3)
c
      data form /'unformatted','unformatted','formatted'/
c
c  set up data sets
c
      write(6,*) 'Enter the two case numbers (1-3)'
      read(5,*) icase1, icase2
c
      write(6,*) 'Enter first file of summaries'
      read(5,'(a)') fin1
      open(2,file=fin1,status='old',form=form(icase1))
c
      write(6,*) 'Enter second file of summaries'
      read(5,'(a)') fin2
      open(3,file=fin2,status='old',form=form(icase2))
c
      write(6,*) 'Enter file of excessive modes on first file'
      read(5,'(a)') fout1
      open(12,file=fout1,status='unknown',form=form(icase1))
c
      write(6,*) 'Enter file of excessive modes on second file'
      read(5,'(a)') fout2
      open(13,file=fout2,status='unknown',form=form(icase2))
c
c  start stepping through dataset
c
      nmd1=0
      nmd2=0
      nmd=1
      iend1=0
      iend2=0
c
      irdmd2=1
c
   30 call rdfreq(icase1,2,cs1,l1,nord1,sig1,frq1,ekin1,ierr)
      if(ierr.ne.0) go to 60
      el1=l1
      nmd1=nmd1+1
c
      if(iend2.eq.1) go to 50
      if(irdmd2.eq.0) go to 40
c
   35 call rdfreq(icase2,3,cs2,l2,nord2,sig2,frq2,ekin2,ierr)
      if(ierr.ne.0) go to 65
      el2=l2
      nmd2=nmd2+1
      irdmd2=0
c
      if(iend1.eq.1) go to 52
c
c  test on l
c
   40 del=el1-el2
      if(abs(del).lt.1.d-6) del=0
c
      if(del) 50,42,52
c
c  same l, test on order
c
   42 if(nord1-nord2) 50,45,52
c
c  same mode. read two new modes
c
   45 irdmd2=1
      go to 30
c
c  output mode 1
c
   50 call wrfreq(icase1,12,cs1,l1,nord1,sig1,frq1,ekin1,ierr)
      iwr=1
      write(6,110) nmd1,el1,nord1,frq1
      go to 55
c
c  output mode 2
c
   52 call wrfreq(icase2,13,cs2,l2,nord2,sig2,frq2,ekin2,ierr)
      iwr=2
      write(6,120) nmd2,el2,nord2,frq2
c
c  determine which mode to read next
c
   55 if(iend1.eq.1) go to 35
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
      stop
  110 format(' mode no',i5,' on d/s 2, with l =',f10.1,'  order =',i5,
     *  ' frequency =',f10.5,' microHz is not on d/s 3')
  120 format(' mode no',i5,' on d/s 3, with l =',f10.1,'  order =',i5,
     *  ' frequency =',f10.5,' microHz is not on d/s 2')
      end
