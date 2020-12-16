      program main
c   
c  select adiabatic model from file
c   
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      parameter(nnmax = 5000)
      implicit double precision (a-h, o-z)
      character fin*280, fout*280
      dimension x(nnmax),aa(6,nnmax),data(8)
c
      write(6,*) 'Enter input file name'
      read(5,'(a)') fin
c   
      open(2,file=fin,status='old',form='unformatted')
c
      write(6,*) 'Enter output file name'
      read(5,'(a)') fout
c   
      open(3,file=fout,status='unknown',form='unformatted')
c
    5 write(6,*) 'Enter model number'
      read(5,*,end=90) nout
c
      write(6,100) fin
      nrec=0
c   
   10 read(2,end=90,err=90) nmod,nn,data,(x(n),(aa(i,n),i=1,5),n=1,nn)
      nrec=nrec+1
      if(nrec.lt.nout) go to 10
c
      write(6,110) nrec,nn,(data(i),i=1,7)
      write(3) nmod,nn,data,(x(n),(aa(i,n),i=1,5),n=1,nn)
c
      go to 5
c   
   90 continue
      stop
c   
  100 format(' file: ',a)
  110 format(/' record no',i4,' output to file'/
     *  '    number of mesh points =',i6,
     *  '   data:'/1p7e11.3)
      end   
