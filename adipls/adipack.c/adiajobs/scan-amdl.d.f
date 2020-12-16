      program main
c   
c  scan adiabatic model file
c   
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      include 'adipr.incl'
      parameter(iaa=10)
      implicit double precision (a-h, o-z)
      character fin*280
      dimension x(nnmax),aa(iaa,nnmax),data(8)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter file name'
      read(5,'(a)') fin
c   
      open(2,file=fin,status='old',
     *  form='unformatted')
c
      if(istdpr.gt.0) 
     *  write(istdpr,100) fin
      nrec=0
c
c  make initial read to test for format
c
      read(2,end=90,err=90) nmod,nn,data
      close(2)
      open(2,file=fin,status='old', form='unformatted')
      idata8 = idint(data(8)+0.1)
      if(idata8.ge.100) then
        ivar = 8
	icase=2
	if(istdpr.gt.0) 
     *  write(istdpr,103) 
      else if(idata8.ge.10) then
        ivar = 6
	icase=1
      else
        ivar = 5
	icase=0
      end if
      icase1=mod(idata8/10,10)
      if(icase1.eq.1) then
	if(istdpr.gt.0) 
     *  write(istdpr,105) 
      else if(icase1.eq.2) then
	if(istdpr.gt.0) 
     *  write(istdpr,106) 
      end if
c   
   10 read(2,end=90,err=90) nmod,nn,data,
     *  (x(n),(aa(i,n),i=1,ivar),n=1,nn)
      nrec=nrec+1
c
      if(istdpr.gt.0) 
     *  write(istdpr,110) nrec,nn,(data(i),i=1,7)
c
c  output surface point and innermost points
c
      if(data(7).ge.0) then 
	nsurf=nn-1
      else
	nsurf=nn
      end if
      if(ivar.eq.5) then
        if(istdpr.gt.0) 
     *  write(istdpr,115) 1-x(nsurf),(aa(i,nsurf),i=1,5)
        if(istdpr.gt.0) 
     *  write(istdpr,116) x(1),(aa(i,1),i=1,5)
      else
        if(istdpr.gt.0) 
     *  write(istdpr,121) 1-x(nsurf),(aa(i,nsurf),i=1,6)
        if(istdpr.gt.0) 
     *  write(istdpr,122) x(1),(aa(i,1),i=1,6)
      end if
c   
      go to 10
c
   90 continue
      stop
c   
  100 format(' file: ',a)
  103 format(/' Model includes asymptotic scaling.')
  105 format(/' Model includes turbulent-pressure correction.')
  106 format(/' Model includes rotation-modified gravity')
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
