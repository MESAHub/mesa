       program main
c   
c  print data from short summary
c   
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      character*280 fin
      dimension ss(7),iss(2)
      equivalence (ss(6),iss(1))
c
c  read input file name
c
      write(6,*) ' input file'
      read(5,'(a)') fin
      write(6,*) ' input file is  ',fin
c
      open(2,file=fin,status='old',
     *  form='unformatted')
c
      nr=0
c   
   10 read(2,end=20,err=20) ss
c
      l=nint(ss(1))
c
c  test for model line
c
      if(l.lt.0) then
        write(6,110) (ss(i),i=3,6)
        go to 10
      end if
c
      nr=nr+1
c
      nord=nint(ss(2))
c   
c  output
c
      write(6,120) l, nord, ss(3),ss(4),ss(5),iss(1),iss(2),nr
c
      go to 10  
c   
   20 continue  
      stop  
  110 format(//' mass =',1pe13.5,'  radius =',e13.5/
     *         ' pc   =',  e13.5,'  rhoc   =',e13.5//
     *         '    l    n   sigma**2        E',
     *         '         nu (mHz)       icase     iorig        nr'/)
  120 format(2i5,1p3e13.5,3i10) 
      end   
