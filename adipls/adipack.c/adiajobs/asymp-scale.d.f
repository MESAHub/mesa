      program main
c
c  Do asymptotic scaling of frequency differences, as output from
c  freqdif, scaling with scale factor obtained from extended amdl
c  file. Allows both scaling with S and S/tau0 (the former
c  required for asymptotic inversion, the latter for direct
c  presentation).
c
c  Note: g and f modes are skipped
c
c  Original version: 27/3/96
c
      parameter(nnmax=3000, ia=8)
      implicit double precision(a-h, o-z)
      character*280 fin, fout, famdl, string
      integer getcas
      dimension data(8),x(nnmax), aa(ia,nnmax),ws(nnmax),ss(1)
c
      cgrav=6.67232d-8
c
c  set  files
c
      write(6,*) 'Enter frequency difference input file'
      read(5,'(a)') fin
      open(2,file=fin,status='old')
c
      write(6,*) 'Enter model input file'
      read(5,'(a)') famdl
      call openfs(3,famdl,'o','u')
c..      open(3,file=famdl,status='old',form='unformatted')
c
      write(6,*) 'Enter frequency difference output file'
      read(5,'(a)') fout
      open(10,file=fout,status='unknown')
c
c  set model data
c
      call rdamdl(3,x,aa,data,nn,nmod,ivar,icry,ia)
c
c  test that model includes scaling
c
      if(ivar.lt.8) then
	write(6,110) famdl
	stop
      end if
c
c  read and decode start of header
c
      read(2,'(a)') string
      itype=getcas(string)
      idiff=mod(itype,10)
      iscale=mod(itype/10,10)
      ierr=mod(itype/100,10)
c
c  test whether differences are already scaled
c
      if(iscale.ne.0) then
	write(6,120) fin
	stop
      end if
c
c  read control of scaling type
c
      write(6,*) 'Enter iscale = 2 for scaling with S'
      write(6,*) '      iscale = 3 for scaling with S/tau0'
      read(5,*) iscale
c
c  output new type, and remainder of header
c
      itype=idiff+10*iscale+100*ierr
      write(10,130) itype,fin
      if(iscale.eq.2) then
	write(10,135)
      else
	write(10,137)
      end if
      write(10,140) famdl
c
c  copy rest of input header
c
      write(10,142)
      call cpycom(2,10)
      write(10,143)
c
      if(ierr.eq.0) then
	write(10,145)
      else
	write(10,147)
      end if
c
c  set scale factor for phase speed and scaling
c
      frqfct=sqrt(cgrav*data(1)/data(2)**3)
      wfct=1.d6*frqfct/(8.d0*atan(1.d0))
      if(iscale.eq.2) then
	sfct=1.d0/frqfct
      else
	sfct=1.d0/aa(8,1)
      end if
c
c  set scaled w (to correspond to nu/(l+1/2) in microHz)
c
      do 10 n=2,nn
   10 ws(n)=wfct*sqrt(aa(1,n)/aa(2,n))
c
c  start going through the input
c  
      call skpcom(2)
c
      write(6,150)
c
   20 if(ierr.eq.0) then
	read(2,*,end=50) l, nord,freq,df,w
      else
	read(2,*,end=50) l, nord,freq,df,w,err
      end if
c
c  skip g and f modes
c
      if(nord.le.0) go to 20
c
c  reset w and interpolate
c
      w=freq/(l+0.5d0)
      call lir(w,ws(2),ss,aa(8,2),1,ia,nn-1,1,inter)
c
      ss(1)=sfct*ss(1)
c
      if(ierr.eq.0) then
	write(10,155) l, nord,freq,ss(1)*df,w
      else
	write(10,155) l, nord,freq,ss(1)*df,w,ss(1)*err
      end if
      write(6,155) l, nord, freq,ss(1)
      go to 20
c
   50 continue
      stop
  110 format(//' ***** Error in asymp-scale. Model contains no scaling'/
     *         '       Model file: ',a60)
  120 format(//' ***** Error in asymp-scale. ',
     *                 'Differences already scaled on file'/
     *         '       ',a60)
  130 format('# Data type :',i10/'#'/
     *       '# Asymptotic scaling'/'#'/
     *       '# Input file: ',a60)
  135 format('# Scale with S (in seconds)'/'#')
  137 format('# Scale with S/tau0'/'#')
  140 format('# Scaling model: ',a60/'#')
  142 format('# ----------------------------------------------------'/
     *       '# Original header:'/'#')
  143 format('# ----------------------------------------------------')
  145 format('# l, n, frequency (microHz), scaled delta freq.,',
     *       ' nu/(l+1/2)'/'#')
  147 format('# l, n, frequency (microHz), scaled delta freq.,',
     *       ' nu/(l+1/2), scaled error'/'#')
  150 format(//' l, n, frequency, scaling'/)
  155 format(2i5,f10.3,1p3e13.5)
      end
