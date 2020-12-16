      program main
c
c  scan adiabatic oscillation modes.
c
c  if icmode = 1 modes are read from full set of eigenfunctions, and
c  if icmode = 2 modes are read from restricted set.
c
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 5/10/90
c
      implicit double precision (a-h, o-z)
      character*280 fin
      dimension xin(10000),yin(6,10000),cs(50), ics(8)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      equivalence(cs(39), ics(1))
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter input file'
      read(5,'(a)') fin
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter 1 for full set, 2 for restricted set'
      read(5,*) icmode
c
      open(2,file=fin,status='old',form='unformatted')
c
      nm=0
c
c  test for case of reading
c
      if(icmode.eq.2) then
c
c  read x
c
        read(2,end=50) nnin,(xin(n),n=1,nnin)
c
      end if
c
c
   10 continue
c
c  test for case
c
      if(icmode.ne.2) then
        read(2,end=50) cs,nnin,(xin(n),(yin(i,n),i=1,6),n=1,nnin)
      else
        read(2,end=50) cs,((yin(i,n),i=1,2),n=1,nnin)
      end if
c
      l=cs(18)+0.5
      nord=cs(19)
c
c  set flag for inclusion of rotation
c
      icase=ics(5)
      irotsl=mod(icase/100,10)
      if(irotsl.eq.1) then
	em=cs(38)
	m=nint(em)
        if(nm.eq.0.and.istdpr.gt.0) write(istdpr,130) fin
      else
        if(nm.eq.0.and.istdpr.gt.0) write(istdpr,135) fin
      end if
c
      nm=nm+1
c
      if(istdpr.gt.0) then
        frq=1000*cs(27)
        if(irotsl.eq.1) then
          write(istdpr,140) nm,l,m,nord,frq
        else
          write(istdpr,145) nm,l,nord,frq
        end if
      end if
c
      go to 10
c
   50 continue
      stop
c
  130 format(' File name: ',a//' n, l, m, order, frequency (microHz)'/)
  135 format(' File name: ',a//' n, l, order, frequency (microHz)'/)
  140 format(4i6,f10.3)
  145 format(3i6,f10.3)
      end
