      program main
c   
c  change between binary and ASCII versions of adiabatic model file
c   
c  Original version: 8/10/97
c
      parameter(nnmax = 10000,iaa=10)
      implicit double precision (a-h, o-z)
      character*280 fin, fout
      dimension x(nnmax),aa(iaa,nnmax),data(8)
c
      write(6,*) 'Enter 1 for transformation from binary to ASCII'
      write(6,*) '      2 for transformation from ASCII to binary'
      read(5,*) icase
      write(6,*) 'Enter input file name'
      read(5,'(a)') fin
      write(6,*) 'Enter output file name'
      read(5,'(a)') fout
c   
      if(icase.eq.1) then
        open(2,file=fin,status='old',form='unformatted')
        open(3,file=fout,status='unknown')
      else
        open(2,file=fin,status='old')
        open(3,file=fout,status='unknown',form='unformatted')
      end if
c
      write(6,100) fin
      nrec=0
c
c  for binary file, make initial read to test for format
c
      if(icase.eq.1) then
        read(2,end=90,err=90) nmod,nn,data
        close(2)
        open(2,file=fin,status='old', form='unformatted')
        idata8 = idint(data(8)+0.1)
        if(idata8.ge.100) then
          ivar = 8
        else if(idata8.ge.10) then
          ivar = 6
        else
          ivar = 5
        end if
c
c  otherwise, read initial record
c
      else
	read(2,102) nmod,nn,ivar
      end if
c
      if(ivar.eq.8) then
	write(6,103) 
      else if(ivar.eq.6) then
	write(6,105) 
      end if
c   
   10 if(icase.eq.1) then
        read(2,end=90,err=90) nmod,nn,data,
     *    (x(n),(aa(i,n),i=1,ivar),n=1,nn)
      else 
	if(nrec.gt.0) read(2,102,end=90,err=90) nmod,nn,ivar
        read(2,107,end=90,err=90) data,
     *    (x(n),(aa(i,n),i=1,ivar),n=1,nn)
      end if
      nrec=nrec+1
c
      write(6,110) nrec,nn,(data(i),i=1,7)
c
c  output surface point and innermost points
c
      if(data(7).ge.0) then 
	nsurf=nn-1
      else
	nsurf=nn
      end if
      if(ivar.eq.5) then
        write(6,115) 1-x(nsurf),(aa(i,nsurf),i=1,5)
        write(6,116) x(1),(aa(i,1),i=1,5)
      else
        write(6,121) 1-x(nsurf),(aa(i,nsurf),i=1,6)
        write(6,122) x(1),(aa(i,1),i=1,6)
      end if
c
c  output model
c
      if(icase.eq.1) then
	write(3,102) nmod,nn,ivar
        write(3,107) data,(x(n),(aa(i,n),i=1,ivar),n=1,nn)
      else
        write(3) nmod,nn,data,(x(n),(aa(i,n),i=1,ivar),n=1,nn)
      end if
c   
      go to 10
c
   90 continue
      stop
c   
  100 format(' Input file: ',a)
  102 format(3i10)
  103 format(/' Model includes asymptotic scaling.')
  105 format(/' Model includes turbulent-pressure correction.')
  107 format(1p4e20.13)
  110 format(/' record no',i4,'    number of mesh points =',i6,
     *  '   data:'/1p7e11.3)
  115 format(' Last (nonsingular) point: 1 - r/R =',1pe13.5/
     *  ' A1 - A5:  ',5e11.3)
  116 format(' Innermost point: r/R =',1pe13.5/
     *  ' A1 - A5:  ',5e11.3)
  121 format(' Last (nonsingular) point: 1 - r/R =',1pe13.5/
     *  ' A1 - A6:  ',6e11.3)
  122 format(' Innermost point: r/R =',1pe13.5/
     *  ' A1 - A6:  ',6e11.3)
      end
