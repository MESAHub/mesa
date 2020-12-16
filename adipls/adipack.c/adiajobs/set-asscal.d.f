      program main
c   
c  set dimensionless and dimensional sound speed from adiabatic 
c  model variables
c  also sets dimensionless integral dr/c and asymptotic scaling
c  into amdl file (aa(7,.) and aa(8,.)) and output result
c  with proper resetting of data(8)
c   
      implicit double precision(a-h,o-z)
      parameter(nnmax=10010, ia=15,iw=4)
      character fin*280, fout*280
      dimension x(nnmax),aa(ia,nnmax),data(8),x1(nnmax),w(iw,nnmax),
     *  aai(ia)
c
      write(6,*) 'Enter unit   file:'
      write(6,*) '        2    input amdl'
      write(6,*) '        3    output amdl'
      write(6,*) 'end with -1 '''''
      call ofiles
c
      call openf(3,'u','u')
c
      write(6,*) 'enter outer radius for asymptotic scaling integrals'
      write(6,*) '(Default: 1)'
      xtrscl = 1.d0
      read(5,*) xtrscl
c   
      ids = 2
      call rdamdl(ids,x,aa,data,nn,nmod,ivar,icry,ia)
c
      cfct=sqrt(6.67232e-8*data(1)/data(2))
c   
      if(x(1).eq.0) then
        n0=2
        aa(10,1)=sqrt(data(3)*aa(3,1)/data(4))/cfct
        aa(11,1)=aa(10,1)*aa(10,1)
        aa(12,1)=1d0/aa(10,1)
      else
        n0=1
      end if
c
      do 10 n=n0,nn
c
      if(x(n).le.xtrscl) ntrscl = n
c   
      aa(11,n)=x(n)*x(n)*aa(1,n)/aa(2,n)
      aa(10,n)=sqrt(aa(11,n))
   10 aa(12,n)=1d0/aa(10,n)
c
      write(6,*) 'xtrscl, ntrscl', xtrscl, ntrscl
c
      call vinta(x,aa(12,1),aa(7,1),nn,ia,ia)   
      tau0tl=aa(7,nn)
      tau0=tau0tl*data(2)/cfct
      write(6,130) tau0, tau0tl
c
c  set asymptotic scale factor
c
c  restrict to xtrscl
c
      do 20 n=1,ntrscl
      x1(n) = x(n)
      w(1,n)=aa(11,n)
   20 w(2,n)=aa(12,n)
c
      if(x(ntrscl).eq.xtrscl) then
	tauscl=aa(7,ntrscl)
      else
	call lir(xtrscl,x,aai,aa,13,ia,nn,1,inter)
	ntrscl=ntrscl+1
	x1(ntrscl)=xtrscl
	w(1,ntrscl)=aai(11)
	w(2,ntrscl)=aai(12)
	tauscl=aai(7)
      end if
c
      do 30 nt=2,ntrscl-1
      xtsq=x(nt)*x(nt)/aa(11,nt)
c
      do 25 n=nt,ntrscl
   25 w(3,n)=1-xtsq*w(1,n)/(x(n)*x(n))
      call squint(x(nt),w(2,nt),w(3,nt),w(4,nt),ntrscl-nt+1,iw,iw,iw)
      aa(8,nt)=w(4,ntrscl)
c
   30 continue
c
c  set central value
c
      aa(8,1)=tauscl
c
c  possibly zero scaling outside ntrscl
c
      if(ntrscl.lt.nn) then
	do 35 n=ntrscl,nn
   35   aa(8,n)=0
      end if
c
c  output new amdl file
c
      idata8=data(8)
      data(8)=100+mod(idata8,100)
c
      write(3) nmod,nn,data,(x(n),(aa(i,n),i=1,8),n=1,nn)
c
   90 continue
      stop
  130 format(//' tau0 =',1pe13.5,'  sec'/
     *         ' tau0 tilde =',e13.5)
      end
