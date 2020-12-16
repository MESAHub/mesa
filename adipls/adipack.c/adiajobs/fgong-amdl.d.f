      program main
c
c  convert model given in standard ASCII GONG format into 
c  input model for adiabatic oscillation package.
c
c  Note: it is assumed that model corresponds to a complete 
c  model. Hence it is extended to zero radius, if that point is
c  not included.
c
c  If second derivatives are not set in glob(11) and glob(12),
c  they are estimated from the central values of the other
c  model variables
c
c  Original version: 25/7/97
c
      implicit double precision(a-h, o-z)
      parameter (nnmax=5000, iconmx=30, ivarmx=50, iaa=6)
      character*280 fin, fout,head
      dimension glob(iconmx), var(ivarmx,nnmax), var1(ivarmx,nnmax),
     *  data(8), aa(iaa,nnmax), ireset(16),q(nnmax),x(nnmax)
      data ireset /3,4,5,6,8,9,10,11,12,13,14,16,17,18,19,20/
c
      cgrav=6.67232d-8
      pi4=16.d0*atan(1.d0)
c
      write(6,*) 'Enter input GONG file'
      read(5,'(a)') fin
      write(6,*) 'Enter output amdl file'
      read(5,'(a)') fout
c
      open(2,file=fin,status='old')
      open(3,file=fout,status='unknown',form='unformatted')
c
      write(6,110)
      do 15 i=1,4
      read(2,'(a)') head
   15 write(6,'(a)') head
c
      read(2,120) nn, iconst, ivar, ivers
      read(2,130) (glob(i),i=1,iconst)
c
      do 20 n=1,nn
   20 read(2,130) (var(i,n),i=1,ivar)
c
      if(var(1,1).gt.var(1,nn)) then 
	nn1=nn+1
	do 30 i=1,ivar
	do 25 n=1,nn
   25   var1(i,n)=var(i,nn1-n)
	do 30 n=1,nn
   30   var(i,n)=var1(i,n)
      end if
c
      if(var(1,1).gt.1.d6) then 
	do 32 i=1,ivar
	do 32 n=1,nn
   32   var1(i,n+1)=var(i,n)
c
	do 34 i=1,ivar
   34   var1(i,1)=0
c
	do 36 ir=1,16
	i=ireset(ir)
   36   var1(i,1)=var1(i,2)
c
        nn=nn+1 
	do 38 i=1,ivar
	do 38 n=1,nn
   38   var(i,n)=var1(i,n)
      end if
c
      do 40 n=1,nn
      q(n)=exp(var(2,n))
   40 x(n)=var(1,n)/glob(2)
c
      x(1)=0
      q(1)=0
c
      do 45 n=2,nn
      aa(1,n)=x(n)
      aa(2,n)=q(n)/x(n)**3
      aa(3,n)=cgrav*glob(1)*q(n)*var(5,n)/ 
     *        (var(10,n)*var(4,n)*var(1,n))
      aa(4,n)=var(10,n)
      aa(5,n)=var(15,n)
   45 aa(6,n)=pi4*var(5,n)*var(1,n)**3/(glob(1)*q(n))
c
      aa(1,1)=0
      aa(2,1)=pi4/3.d0*var(5,1)*glob(2)**3/glob(1)
      aa(3,1)=0
      aa(4,1)=var(10,1)
      aa(5,1)=0
      aa(6,1)=3.d0
c
      if(aa(5,nn).le.10) then 
        nn=nn-1 
        write(6,*) 'Chop off outermost point' 
      end if
c
      data(1)=glob(1)
      data(2)=glob(2)
      data(3)=var(4,1)
      data(4)=var(5,1)
c
      if(glob(11).lt.0.and.glob(11).gt.-10000) then 
        data(5)=-glob(11)/var(10,1)
	data(6)=-glob(12) 
      else 
        data(5)=pi4/3.d0*cgrav*(var(5,1)*glob(2))**2/
     *          (var(4,1)*var(10,1))
	d2amax=0.d0
	do 50 n=2,nn
	d2amax=max(d2amax,aa(5,n)/x(n)**2)
	if(x(n).ge.0.05d0) go to 55
   50   continue
   55   data(6)=d2amax+data(5)
	write(6,140) data(5), data(6)
      end if
c
      data(7)=-1.d0
      data(8)=0.d0
c
      nmodel=1
      write(3) nmodel,nn,(data(i),i=1,8),((aa(i,n),i=1,6),n=1,nn)
      stop
  110 format(/' GONG model header:'/)
  120 format(4i10)
  130 format(5e16.9)
  140 format(' data(5), data(6) reset to ',1p2e13.5)
      end
