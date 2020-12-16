      program main
c
c  compare adiabatic oscillation variables at fixed r/rphot
c
c  modified 16/7/1985 to take into account changes in autograph package
c
c  modified 10/8/87 to take out all plot instructions. Now outputs
c  data for plot on file.
c
c  Modified 3/8/92, to allow comparing models from same file
c
c  Modified 28/4/97, allowing comparisons at fixed pressure 
c  or mass.
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 17/4/90
c
      implicit double precision (a-h, o-z)
      parameter(nnmax=5000,iaa=10)
      character title*80, absci*40, sound*5, tail*5,
     *  trim*40, file*80, ccase*4
      dimension x(nnmax,2),aa(iaa,nnmax,2),data(8,2),nn(2),
     *  xpl(nnmax),ypl(nnmax,8),aai(iaa),xdif(nnmax,2),
     *  modl(2),modp(2),w(5,nnmax),absci(3),sound(2),ccase(7)
      common/cofile/ nfiles, idsfil(20), file(20), iopen(20)
      common/ccgrav/ cgrav
c
      data ids1p, ids2p /-1, -1/
c
      data absci /'x','p','acoustical radius (min)'/
      data sound /'c','c**2'/
      data ccase /'r/R', 'q', 'r/Rs', 'r', 'p', 'm/Ms', 'm/Mr'/
      data tail /'    @'/
c
c  hardcoded gravitational constant
c
      cgrav=6.67232d-8
c
c  defaults in /exec/:
c  ******************
c
c  iread: if iread = 1 (or in first pass) read new models
      iread=0
c  ids1, ids2, mod1, mod2: take difference
c   (mod2 on d/s ids2) - (mod1 on d/s ids1)
      ids1=2
      ids2=3
      mod1=1
      mod2=1
c  nprt: if nprt .gt. 1 print differences at nprt points
      nprt=100
c  icase: determines what is held fixed in difference
c  icase = 1: differences at fixed r/R
c  icase = 2: differences at fixed q
c  icase = 3: differences at fixed r/r(last point)
c  icase = 4: differences at fixed r
c  icase = 5: differences at fixed p
c  icase = 6: differences at fixed m/M(last point)
c  icase = 7: differences at fixed m/M(fixed r/R)
c  (Note: not all of these are as yet implemented)
      icase=1
c  iabsci = 1: x as abscissa
c  iabsci = 2: p as abscissa, logarithmic scale
c  iabsci = 3: acoustical radius as abscissa
      iabsci=1
c  isound: when isound = 1 sound speed difference is output, otherwise
c          twice its value
      isound=2
c
c  ......................................................................
c
      iw=5
c
c  initialize indices for previous models
c
      modp(1)=0
      modp(2)=0
c
      ifirst=1
c
      write(6,*) 'Files needed:'
      write(6,*) 'Model files on units ids1 and ids2'
      write(6,*)  '(defaults ',ids1,' and',ids2,' )'
      write(6,*) 'Plot output on unit 11.'
      call ofiles
c
      call openf(11,'u','f')
c
10000 continue
c
      write(6,*) 'iread,ids1,ids2?'
      write(6,*) iread,ids1,ids2
      read(5,*,end=90,err=90) iread,ids1,ids2
      write(6,*) 'mod1,mod2,nprt?'
      write(6,*) mod1,mod2,nprt
      read(5,*,end=90,err=90) mod1,mod2,nprt
      write(6,*) 'icase,iabsci,isound?'
      write(6,*) icase,iabsci,isound
      read(5,*,end=90,err=90) icase,iabsci,isound
c
c      write(6,100)
c
      write(6,*) 'iread,ids1,ids2'
      write(6,*) iread,ids1,ids2,tail
      write(6,*) 'mod1,mod2,nprt'
      write(6,*) mod1,mod2,nprt,tail
      write(6,*) 'icase,iabsci,isound'
      write(6,*) icase,iabsci,isound,tail
c
c  test for skipping model read
c
      if(iread.eq.0.and.ifirst.ne.1) go to 41000
c
c
c  open model files
c
      if(ids1.ne.ids1p) then
        if(ids1p.gt.0) close(ids1p)
        call openf(ids1,'o','u')
        ids1p=ids1
      end if
c
      if(ids2.ne.ids2p.and.ids2.ne.ids1) then
        if(ids2p.gt.0) close(ids2p)
        call openf(ids2,'o','u')
        ids2p=ids2
      end if
c
      ifirst=0
c
c
c  store mod1 and mod2
c
      modl(1)=mod1
      modl(2)=mod2
c
c  read models
c
      ids=ids1
      do 10 i=1,2
c
c  test for rewinding
c
      if(modl(i).gt.modp(i)) go to 5
      rewind ids
      modp(i)=0
c
    5 nrd=modp(i)
c
    7 read(ids,end=90,err=90) nmod,nnr,(data(j,i),j=1,8),(x(n,i),
     *   (aa(j,n,i),j=1,5),n=1,nnr)
c
      nrd=nrd+1
      if(nrd.lt.modl(i)) go to 7
c
c  correct model found
c
      modp(i)=nrd
      nn(i)=nnr
c
      if(ids2.eq.ids1.and.i.eq.1) then
	modp(2)=nrd
      end if
   10 ids=ids2
c
c  output file names to result file
c
      write(11,103)
      ids=ids1
      do 11 i=1,2
      call stfile(ids,nfin)
      write(11,104) file(nfin)
   11 ids=ids2
c
      write(6,105) (ids,nn(ids),(data(i,ids),i=1,8),ids=1,2)
c
      write(11,106) (ids,nn(ids),(data(i,ids),i=1,8),ids=1,2)
c
      if(mod1.gt.1.or.mod2.gt.1) then
        write(6,140) mod2,ids2,mod1,ids1
        write(11,142) mod2,ids2,mod1,ids1
      end if
c
c  set acoustical radius and mass for models
c
      idgac=1
      do 15 i=1,2
      call acoust(x(1,i),aa(1,1,i),data(1,i),nn(i),w,iaa,iw,
     *  idgac,5)
c  store in aa
      nni=nn(i)
      do 12 n=1,nni
   12 aa(6,n,i)=w(1,n)
c
c  surface dimensionless value
c
      racs1=w(2,nni)
c
      call massfr(x(1,i),aa(1,1,i),data(1,i),nn(i),w,xpl,iaa,iw,
     *  idgac,5)
c  store in aa
      nni=nn(i)
      do 13 n=1,nni
   13 aa(7,n,i)=w(1,n)
c
c  set p
c
      pfct=data(1,i)/(data(2,i)*data(2,i))
      pfct=5.31036d-9*pfct*pfct
      if(x(1,i).eq.0) then
	n1=2
	aa(8,1,i)=data(3,i)
      else
	n1=1
      end if
      do 15 n=n1,nni
      xx2=x(n,i)*x(n,i)
   15 aa(8,n,i)=pfct*aa(5,n,i)*aa(1,n,i)*aa(1,n,i)*xx2/(aa(3,n,i)*
     *  aa(2,n,i))
c
      nn1=nn(1)
      nn2=nn(2)
c
c  test for zero surface, and, if so, reset nn1
c
      if(aa(5,nn1,1).le.1.d-30) nn1=nn1-1
c
c  if one of the models is an envelope model, shift acoustical
c  radius of model 2 to have same surface value as model 1
c
      if(x(1,1).gt.0.or.x(1,2).gt.0) then
c
        dacst=aa(6,nn1,1)-aa(6,nn2,2)
        write(6,107) dacst
        write(11,108) dacst
        do 23 n=1,nn2
   23   aa(6,n,2)=aa(6,n,2)+dacst
c
      end if
c
c  compare models at fixed x. interpolate to grid in model 1
c  when interpolating in full model, rescale aa(2,.) and aa(4,.)
c
      if(x(1,2).eq.0) then
	do 25 n=2,nn2
	xx2=x(n,2)*x(n,2)
	aa(2,n,2)=aa(2,n,2)/xx2
   25   aa(4,n,2)=aa(4,n,2)/xx2
	aa(2,1,2)=aa(2,2,2)
	aa(4,1,2)=aa(4,2,2)
      end if
c
c  set variable to hold fixed in differencing
c
      do 40 i=1,2
      nns=nn(i)
      if(icase.eq.1) then
	do 31 n=1,nns
   31   xdif(n,i)=x(n,i)
      else if(icase.eq.2) then
	do 32 n=1,nns
   32   xdif(n,i)=aa(7,n,i)
      else if(icase.eq.3) then
	do 33 n=1,nns
   33   xdif(n,i)=x(n,i)/x(nns,i)
      else if(icase.eq.4) then
	do 34 n=1,nns
   34   xdif(n,i)=x(n,i)*data(2,i)
      else if(icase.eq.5) then
	do 35 n=1,nns
   35   xdif(n,i)=aa(8,n,i)
      else if(icase.eq.6) then
	do 36 n=1,nns
   36   xdif(n,i)=aa(7,n,i)-aa(7,nns,i)
      else
	write(6,109) icase
	stop
      end if
c
   40 continue
c
      np=0
c
      nfirst=1
      do 50 n=1,nn1
c
c  do not extrapolate
c
      if((xdif(n,1)-xdif(nn2,2))*(xdif(n,1)-xdif(1,2)).gt.0) go to 50
c
      call lir(xdif(n,1),xdif(1,2),aai,aa(1,1,2),6,iaa,nn2,nfirst,inter)
      nfirst=2
c
c  if interpolating in full model, scale back
c
      xx2=x(n,1)*x(n,1)
      if(x(1,2).eq.0) then
	aai(2)=xx2*aai(2)
	aai(4)=xx2*aai(4)
      end if
c
      write(41,'(f10.5,1p10e13.5)') x(n,1),(aa(i,n,1),i=1,5),
     *  (aai(i),i=1,5)
c
c  set and print differences
c
      np=np+1
      do 48 i=1,5
      i1=i
      if(i.eq.4) i1=2
   48 ypl(np,i)=(aai(i)-aa(i,n,1))/max(aa(i1,n,1),1.d-37)
c
c  acoustical radius difference relative to surface value
c
      ypl(np,6)=(aai(6)-aa(6,n,1))/aa(6,nn1,1)
c
c  set sound speed or its square 
      if(xx2.gt.0) then
        c1sq=aa(1,n,1)*xx2/aa(2,n,1)
        c2sq=aai(1)*xx2/aai(2)
      else
        c1sq=data(2,1)*data(3,1)*aa(3,1,1)/(cgrav*data(1,1)*data(4,1))
        c2sq=data(2,2)*data(3,2)*aa(3,1,2)/(cgrav*data(1,2)*data(4,2))
      end if
      write(40,*) xx2,c1sq,c2sq
      c1=sqrt(c1sq)
      c2=sqrt(c2sq)
      if(isound.eq.1) then
        ypl(np,7)=(c2-c1)/max(c1,1.d-37)
      else
        ypl(np,7)=(c2sq-c1sq)/max(c1sq,1.d-37)
      end if
c
c  acoustical radius integrand
c
      ypl(np,8)=0.5d0*ypl(np,7)/(c1*racs1)
c
c  set possible abscissas in w
c
      w(3,np)=x(n,1)
      w(2,np)=aa(8,np,1)
      w(1,np)=aa(6,np,1)
   50 continue
c
      nnp=np
c
41000 do 52 n=1,nnp
c
c  set abscissa
c
      if(iabsci.eq.2) then
        xpl(n)=w(2,n)
      else if(iabsci.eq.3) then
        xpl(n)=w(1,n)
      else
        xpl(n)=w(3,n)
      end if
c
   52 continue
c
c  test for printing
c
      if(nprt.gt.1) then
c
        ndpr=max0(1,(nnp-1)/(nprt-1))
	write(6,110) ccase(icase)
        write(6,113)
        write(6,120) (n,w(3,n),xpl(n),(ypl(n,i),i=1,8),n=1,nnp,ndpr)
c
      end if
c
c  output plot variables
c
      write(11,125) ccase(icase)
      write(11,130) sound(isound),absci(iabsci)
      do 62 n=1,nnp
   62 write(11,135) xpl(n),(ypl(n,i),i=1,8)
c
      go to 10000
c
   90 continue
      stop
  100 format(1h1)
  101 format(a)
  102 format('# ',a)
  103 format('# Data files:'/'#')
  104 format('#   ',a60)
  105 format(///' adiabatic models compared.'//(' model no',i3,
     *  '   nn =',i4,'  data:',1p8e12.4))
  106 format('#'/'#'/'# adiabatic models compared.'/'#'/
     *   ('# model no',i3,'  nn =',i4,'  data:',1p8e12.4))
  107 format(//' one of the models is an envelope model.'/
     *  ' acoustical radius of model 2 shifted by',f10.5,' mins')
  108 format('#'/'#'/'# one of the models is an envelope model.'/
     *  '# acoustical radius of model 2 shifted by',f10.5,' mins')
  109 format(//' ***** Error in compamod. icase =',i5,
     *  ' not yet implemented')
  110 format(//'Differences at fixed ',a)
  113 format(///' n, x,abscissa,relative difference. a(1-5),',
     *  ' acoustical radius, c, acoustical integrand:'/)
  120 format(i4,0pf12.7,1p9e12.4)
  125 format('#'/'#'/'# Differences at fixed ',a)
  130 format('#'/'#'/'# relative differences for'/
     *  '# q/x3, Vg, Gamma1, A, U, t acous, ',a5,', i acoust'/'#'/
     *  '# ',a10,'  relative differences:'/'#')
  132 format('#'/'#',a,'  ',a/'#')
  135 format(1pe14.6,8e12.4)
  140 format(/'model',i3,' on d/s',i3,'  -  model',i3,
     *  ' on d/s',i3)
  142 format('#'/'# model',i3,' on d/s',i3,'  -  model',i3,
     *  ' on d/s',i3)
      end
      subroutine acoust(x,aa,data,nn,w,iaa,iw,idiag,intr)
c  sets acoustical distance from centre, in mins, into w(1,.), on the
c  basis of adiabatic oscillation variables in aa and data.
c  if idiag .ge. 1 acoustical radius is printed, and if in addition
c  idiag .ge. 2 acoustical distance is printed, with interval intr.
c
c  w(2,.) is used as work array.
c
      implicit double precision (a-h, o-z)
      dimension x(*),aa(iaa,*),w(iw,*),data(*)
      common/ccgrav/ cgrav
c
      taufct=sqrt(data(2)**3/(cgrav*data(1)))/60.d0
c  first and final point
      nn1=nn
      if(aa(5,nn).lt.1.d-30) nn1=nn-1
      n1=1
      if(x(1).eq.0) n1=2
c
      do 25 n=n1,nn1
   25 w(1,n)=aa(1,n)*x(n)*x(n)/aa(2,n)
c
      gmm=aa(3,1)
      if(n1.eq.1) go to 27
      w(1,1)=data(2)*gmm*data(3)/(cgrav*data(1)*data(4))
c
c  integration
c
c  set acoustical radius
   27 do 30 n=1,nn1
   30 w(1,n)=1.d0/sqrt(w(1,n))
      call vinta(x,w(1,1),w(2,1),nn1,iw,iw)
c  contribution from singular surface
      if(nn.eq.nn1) go to 40
      ts=2*(x(nn)-x(nn1))*w(1,nn1)
      w(2,nn)=w(2,nn1)+ts
c
   40 do 50 n=1,nn
   50 w(1,n)=w(2,n)*taufct
      if(idiag.lt.1) return
      write(6,100) w(2,nn),w(1,nn)
      if(idiag.ge.2.and.intr.gt.0) write(6,110) (n,x(n),w(2,n),
     *  w(1,n),n=1,nn,intr)
      return
  100 format(//' dimensionless acoustical radius =',1pe15.7,
     *  ' dimensional value =',1pe15.7,' mins')
  110 format(//' n, x, dimensionless, dimensional acoustical radius ',
     *  '(mins):'//(i4,0pf12.6,1p2e15.7))
      end
      subroutine massfr(x,aa,data,nn,w,xw,iaa,iw,idiag,intr)
c  sets log10(m/M) into w(1,.), on the
c  basis of adiabatic oscillation variables in aa and data.
c  m = M is assumed where r = R, as given in data(2)
c
c  w(2,.) is used as work array.
c
      implicit double precision (a-h, o-z)
      dimension x(*),aa(iaa,*),w(iw,*),xw(*),data(*)
      dimension wint(1)
c
      write(6,*) 'Entering massfr,idiag,nn =',idiag,nn
c
      if(x(1).eq.0) then 
	ns=2
      else
	ns=1
      end if
c
      nn1=nn+1
      do 10 n=ns,nn
      n1=nn1-n
      xw(n1)=x(n)
   10 w(2,n1)=x(n)*x(n)*aa(1,n)*aa(5,n)
c
      if(x(1).eq.0) then
	xw(nn)=0
	w(2,nn)=0
      end if
c
      call vinta(xw,w(2,1),w(1,1),nn,iw,iw)
c
c  shift to set zero-point at photosphere
c
      call lir(1.d0,xw,wint,w,1,iw,nn,1,inter)
      do 15 n=1,nn
   15 w(3,n)=-w(1,n)+wint(1)
c
c  rescale to logarithmic variables
c
      do 20 n=1,nn
      n1=nn1-n
      ammi=1.d0/log(10.d0)
      if(abs(w(3,n)).le.0.01d0) then
	w(1,n1)=-ammi*w(3,n)*(1+w(3,n)*(0.5d0+w(3,n)*(0.33333333333d0+
     *          0.25d0*w(3,n))))
      else if(w(3,n).lt.1.d0) then
	w(1,n1)=log10(1.d0-w(3,n))
      else
	w(1,n1)=-40
      end if
   20 continue
      if(idiag.lt.2) return
      write(6,110) (n,x(n),w(2,nn1-n),w(3,nn1-n),w(1,n),
     *  n=1,nn,intr)
  110 format(//' n, x, integrand, delta q, log10(q):'//
     *  (i4,0pf12.6,1p3e15.7))
      end
