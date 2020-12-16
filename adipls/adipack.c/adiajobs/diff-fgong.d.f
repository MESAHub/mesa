      program main
c
c  finds differences between formatted GONG models as supplied
c  by the participants.
c
c  original version: 2/9/88.
c
c  modified 20/1/89 to allow differences at fixed q.
c
c  modified 23/1/89 to interpolate in q rather than shifted ln(q)
c
c  modified 1/3/89 to do special interpolation at innermost meshpoint
c
c  modified 7/12/90, to allow differences at fixed pressure.
c
c  modified 23/9/05, adding q = m/M to output as column 2 (and hence
c  shifting everything else).
c
c  Note: before 15/3/89 difference no. 17 was erroneously labelled
c     delta(ln r/R). It is in fact delta(ln r)
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 13/3/90
c
      implicit double precision (a-h, o-z)
      include 'evolpr.incl'
      parameter(ndiff=17,ivarmx=40)
      character*280 cdata,fin1,fin2, fout, head, line
      character dyname*20, ccase*7
      dimension nmodel(2),varout(ivarmx,nnmax,2),cdata(4,2),
     *  varrd(ivarmx,nnmax),head(4,2),datout(15,2),nvarrd(2),nnrd(2),
     *  xab1(nnmax),xab2(nnmax),xr1(nnmax),xr2(nnmax),
     *  xq1(nnmax),xq2(nnmax),
     *  xabc1(nnmax),xabc2(nnmax),yint(ivarmx),varc1(ivarmx,20),
     *  dy(ndiff),xmax(ndiff),dymin(ndiff),dymax(ndiff),daymax(ndiff),
     *  dyname(ndiff),ccase(5)
c
      data dyname /
     *  'delta(ln q)',
     *  'delta(ln T)',
     *  'delta(ln p)',
     *  'delta(ln rho)',
     *  'delta(X)',
     *  'delta(ln L)',
     *  'delta(ln kappa)',
     *  'delta(ln epsilon)',
     *  'delta(ln GAMMA1)',
     *  'delta(ln ad. grad)',
     *  'delta(ln delta)',
     *  'delta(ln cp)',
     *  'delta(ln mue**(-1))',
     *  'delta(ln A)',
     *  'delta(ln rX)',
     *  'delta(ln c)',
     *  'delta(ln r)'/
      data ccase /'r/R', 'q', 'r/Rs', 'p','(R-r)/R'/
c
      fxp(x)=exp(min(max(x,-80.d0),80.d0))
c
      write(6,*) 'Enter first file name'
      read(5,'(a)') fin1
      open(2,file=fin1,status='old')
      write(6,*) 'Enter second file name'
      read(5,'(a)') fin2
      open(3,file=fin2,status='old')
c
      nmodel(1)=1
      nmodel(2)=1
      write(6,*) 'Enter model numbers on first and second file'
      read(5,*) nmodel
c
      write(6,100) nmodel(1),fin1,nmodel(2),fin2
c
      write(6,*) 'Enter output file'
      read(5,'(a)') fout
      open(10,file=fout,status='unknown')
      write(10,102) nmodel(1),fin1,nmodel(2),fin2
c
c  file for summary
c
      open(12,file='dgong-sum',status='unknown')
      write(6,108)
      write(12,100) nmodel(1),fin1,nmodel(2),fin2
c
      write(6,*) 'Enter 1 for differences at fixed r/R, 2 for fixed q'
      write(6,*) '3 for differences at fixed r/r(last point)'
      write(6,*) '4 for differences at fixed pressure'
      write(6,*) '5 for differences at fixed (R - r)/R'
      icase=1
      read(5,*) icase
c
      write(6,104) ccase(icase)
      write(10,105) ccase(icase)
      write(12,104) ccase(icase)
c
      write(6,*) 'Enter 1 for diagnostics from reading models'
      write(6,*) '      2 for diagnostics from interpolation'
      itest=0
      read(5,*) itest
c
      write(6,*) 'Enter 1 to allow extrapolation, 0 to block it'
      iextrp=0
      read(5,*) iextrp
c
c
      do 20 ids=2,3
      nrd=0
      ids1=ids-1
   15 read(ids,115,end=20) (head(i,ids1),i=1,4)
      read(ids,*) nn,iconst,ivar,ivers
c..      read(ids,117) (datout(i,ids1),i=1,iconst)
      read(ids,*) (datout(i,ids1),i=1,iconst)
      do 16 n=1,nn
c..      read(ids,117) (varrd(i,n),i=1,ivar)
      read(ids,*) (varrd(i,n),i=1,ivar)
      if(itest.eq.1) write(6,118) ids1,n
   16 continue
c
c  store such that r increases with n
c
      if(varrd(1,1).lt.varrd(1,nn)) then
        do 17 n=1,nn
        do 17 i=1,ivar
   17   varout(i,n,ids1)=varrd(i,n)
      else
        do 18 n=1,nn
        n1=nn+1-n
        do 18 i=1,ivar
   18   varout(i,n1,ids1)=varrd(i,n)
      end if
c
      nvarrd(ids1)=ivar
      nnrd(ids1)=nn
      nrd=nrd+1
      if(nrd.lt.nmodel(ids1)) go to 15
c
   20 continue
c
      write(6,110) nnrd
c
      nvar=min0(nvarrd(1),nvarrd(2))
c
      ndiffa=17
c
c  reset log q and log L, and set abscissa for the two models
c
      cq=2
      dq=1
c..      qshft  =fxp(varout(2,2,1))
c..      write(6,*) 'qshft =',qshft
c
      do 25 n=1,nnrd(1)
      if(icase.ne.3) then
        xr1(n)=varout(1,n,1)/datout(2,1)
      else
        xr1(n)=varout(1,n,1)/varout(1,nnrd(1),1)
      end if
      xq1(n)=exp(varout(2,n,1))
      if(icase.eq.1.or.icase.eq.3) then
        xab1(n)=xr1(n)
        xabc1(n)=xr1(n)*xr1(n)
      else if(icase.eq.2) then
        if(varout(1,n,1).gt.0) then
          xab1(n)=
     *      (fxp(varout(2,n,1)/3)+cq)*varout(2,n,1)/(dq-varout(2,n,1))
          xabc1(n)=fxp(0.6666667*varout(2,n,1))
        else
          xab1(n)=-cq
          xabc1(n)=0
        end if
      else if(icase.eq.4) then
        xab1(n)=-log10(varout(4,n,1))
      else
	xab1(n)=varout(17,n,1)/datout(2,1)
      end if
      varout(2,n,1)=fxp(varout(2,n,1))
   25 continue
c
c  for normalization, check that varout(1,*,2) is physical radius
c
      r2surf=varout(1,nnrd(2),2)
      do 27 n=1,nnrd(2)
      if(icase.ne.3.and.r2surf.ge.1.d5) then
        xr2(n)=varout(1,n,2)/datout(2,2)
      else
        xr2(n)=varout(1,n,2)/r2surf
      end if
      xq2(n)=exp(varout(2,n,2))
      if(icase.eq.1.or.icase.eq.3) then
        xab2(n)=xr2(n)
        xabc2(n)=xr2(n)*xr2(n)
      else if(icase.eq.2) then
        if(varout(1,n,2).gt.0) then
          xab2(n)=
     *      (fxp(varout(2,n,2)/3)+cq)*varout(2,n,2)/(dq-varout(2,n,2))
          xabc2(n)=fxp(0.6666667*varout(2,n,2))
        else
          xab2(n)=-cq
          xabc2(n)=0
        end if
c
c  the following line is most peculiar and almost certainly wrong!
c
c..        xab1(n)=-log10(varout(4,n,1))
      else if(icase.eq.4) then
        xab2(n)=-log10(varout(4,n,2))
      else
	xab2(n)=varout(17,n,2)/datout(2,2)
      end if
      varout(2,n,2)=fxp(varout(2,n,2))
   27 continue
c
c..      write(6,'(//'' xab1:''/(1p5e13.5))') (xab1(n),n=1,nnrd(1))
c..      write(6,'(//'' xab2:''/(1p5e13.5))') (xab2(n),n=1,nnrd(2))
c
c  set special variables in innermost part of model 1, for
c  special interpolation
c
      do 32 n=1,20
      do 32 i=1,ivar
   32 varc1(i,n)=varout(i,n,1)
c
      if(icase.eq.1.or.icase.eq.3) then
        varc1(2,1)=(4.1887902e33/datout(1,1))*varout(5,1,1)
        varc1(7,1)=4.1887902*varout(5,1,1)*varout(9,1,1)
        varc1(15,1)=1.e22*(datout(11,1)/varout(10,1,1)-datout(12,1))/
     *    datout(2,1)**2
        do 34 n=2,20
        varc1( 1,n)=varout(1,n,1)**2
        varc1( 2,n)=1.e33*varout( 2,n,1)/varout(1,n,1)**3
        varc1( 7,n)=      varout( 7,n,1)/varout(1,n,1)**3
   34   varc1(15,n)=1.e22*varout(15,n,1)/varout(1,n,1)**2
c
        write(6,112) (n,xabc1(n),varc1(2,n),varc1(7,n),varc1(15,n),
     *    n=1,20)
      else if(icase.eq.2) then
        varc1(1,1)=
     *    1.e-11*(datout(1,1)/(4.1887902*varout(5,1,1)))**0.3333333
        varc1(7,1)=1.e-33*datout(1,1)*varout(9,1,1)
        varc1(15,1)=(datout(11,1)/varout(10,1,1)-datout(12,1))*
     *    (datout(1,1)/(4.1887902*varout(5,1,1)))**0.6666667/
     *    datout(2,1)**2
        do 36 n=2,20
        varc1(1,n)=1.e-11*varout(1,n,1)/varout(2,n,1)**0.3333333
        varc1(2,n)=varout(2,n,1)**0.6666667
        varc1(7,n)=1.e-33*varout(7,n,1)/varout(2,n,1)
   36   varc1(15,n)=      varout(15,n,1)/xabc1(n)
c
        write(6,113) (n,xabc1(n),varc1(1,n),varc1(7,n),varc1(15,n),
     *    n=1,20)
      end if
c
      do 45 ids1=1,2
      write(10,120) ids1
      write(10,125) (head(i,ids1),i=1,4)
      write(12,126) ids1,(head(i,ids1),i=1,4)
      write(12,127) nnrd(ids1),(datout(i,ids1),i=1,10)
   45 write(10,130) (datout(i,ids1),i=1,iconst)
c
c  step through model 2, interpolating in model 1
c
      amm=log(10.d0)
      write(10,150) ndiffa,(i,dyname(i),i=1,ndiffa)
      write(10,155)
c
      nint=0
c
      do 60 n=1,nnrd(2)
c
c  interpolation in r/R or q
c
      if(icase.gt.3.or.xab2(n).gt.xab1(4)) then
        nint=nint+1
        call lir(xab2(n),xab1,yint,varout(1,1,1),ivar,ivarmx,nnrd(1),
     *    nint,inter)
c
      else
c
c  for innermost points, do interpolation in (r/R)**2 or q**(2/3).
c  to get more accurate values for q and L
c  inside innermost non-zero meshpoint use linear interpolation
c
	write(6,*) 'special interpolation at n =',n
        nint=0
        if(xab2(n).gt.xab1(2)) then
          call lir(xabc2(n),xabc1,yint,varc1(1,1),ivar,ivarmx,20,
     *      nint,inter)
        else
          call lir1(xabc2(n),xabc1,yint,varc1(1,1),ivar,ivarmx,20,
     *      nint,inter)
        end if
c
c  reset variables, depending on icase
c
        if(n.eq.1) then
          yint( 1) = 0
          yint( 2) = 0
          yint( 7) = 0
          yint(15) = 0
        else if(icase.eq.1.or.icase.eq.3) then
          yint( 1)=sqrt(yint(1))
          yint( 2)=1.e-33*yint( 2)*yint(1)**3
          yint( 7)=       yint( 7)*yint(1)**3
          yint(15)=1.e-22*yint(15)*yint(1)**2
        else
          yint(2)=yint(2)**1.5
          yint(1)=1.e11*yint(1)*yint(2)**0.3333333
          yint(7)=1.e33*yint(7)*yint(2)
          yint(15)=     yint(15)*xabc2(n)
        end if
      end if
c
      if(itest.eq.2) then
        write(6,140) (varout(i,n,2),i=1,ivar)
        write(6,142) (yint(i),i=1,ivar)
      end if
c
c  test for skipping extrapolation
c
      if(iextrp.eq.0.and.inter.eq.0) then
	write(6,145) xab2(n)
	go to 60
      end if
c
      if(varout(2,n,2).gt.0.and.yint(2).gt.0) then
        dy(1)=-log(yint(2)/varout(2,n,2))
      else
        dy(1)=0
      end if
      dy(2)=-log(yint(3)/varout(3,n,2))
      dy(3)=-log(yint(4)/varout(4,n,2))
      dy(4)=-log(yint(5)/varout(5,n,2))
      dy(5)=     -yint(6)+varout(6,n,2)
      if(varout(7,n,2).gt.0.and.yint(7).gt.0) then
        dy(6)=-log(yint(7)/varout(7,n,2))
      else
        dy(6)=0
      end if
      dy(7)=-log(yint(8)/varout(8,n,2))
      if(abs(varout(9,n,2)).gt.1.e-5) then
        dy(8)=-log(abs(yint(9)/varout(9,n,2))+1.e-20)
      else
        dy(8)=0
      end if
      dy(9)=-log(yint(10)/varout(10,n,2))
      dy(10)=-log(yint(11)/varout(11,n,2))
      if(varout(12,n,2).gt.0.and.yint(12).gt.0) then
        dy(11)=-log(yint(12)/varout(12,n,2))
      else
	dy(11)=0
      end if
      if(varout(13,n,2).gt.0.and.yint(13).gt.0) then
        dy(12)=-log(yint(13)/varout(13,n,2))
      else
	dy(12)=0
      end if
      if(varout(14,n,2).gt.1.e-20.and.yint(14).gt.0) then
        dy(13)=-log(yint(14)/varout(14,n,2))
      else
        dy(13)=0
      end if
      if(abs(varout(15,n,2)).gt.1.e-20.and.abs(yint(15)).gt.1.e-20) then
        dy(14)=-log(abs(yint(15)/varout(15,n,2))+1.e-20)
      else
        dy(14)=0
      end if
      if(abs(varout(16,n,2)).gt.1.e-20.and.abs(yint(16)).gt.1.e-20) then
        dy(15)=-log(abs(yint(16)/varout(16,n,2))+1.e-20)
      else
        dy(15)=0
      end if
      dy(16)=0.5*(dy(3)+dy(9)-dy(4))
      if(varout(1,n,2).gt.1.e-20.and.yint(1).gt.0) then
        dy(17)=-log(yint(1)/varout(1,n,2))
      else
        dy(17)=0
      end if
c
c  set extreme values
c
      if(n.eq.1) then
        do 55 i=1,ndiffa
        xmax(i)=xr2(1)
        dymax(i)=dy(i)
        dymin(i)=dy(i)
   55   daymax(i)=abs(dy(i))
c
      else
c
        do 57 i=1,ndiffa
        if(abs(dy(i)).gt.daymax(i)) then
          xmax(i)=xr2(n)
          daymax(i)=abs(dy(i))
        end if
        if(dy(i).gt.dymax(i)) then
          dymax(i)=dy(i)
        else if(dy(i).lt.dymin(i)) then
          dymin(i)=dy(i)
        end if
   57   continue
c
      end if
c
      write(10,160) xr2(n),xq2(n),(dy(i),i=1,ndiffa)
   60 continue
c
c  output maximum values
c
      write(6,170) (i,dyname(i),xmax(i),dymin(i),dymax(i),i=1,ndiffa)
      write(12,170) (i,dyname(i),xmax(i),dymin(i),dymax(i),i=1,ndiffa)
      write(0,
     *  '(/'' ***** Note: file output includes m/M as column 2'')')
c
   90 continue
c
      stop
  100 format(//' differences between GONG models.'/
     *  ' First  model: no.',i3,' on file ',a40/
     *  ' Second model: no.',i3,' on file ',a40)
  102 format('#'/'#'/'# differences between GONG models.'/
     *  '# First  model: no.',i3,' on file ',a40/
     *  '# Second model: no.',i3,' on file ',a40)
  104 format(//' differences at fixed ',a)
  105 format('#'/'#'/'# differences at fixed ',a)
  108 format(//' summary written to file dgong-sum')
  110 format(//' end reading models. nn1, nn2=',2i5)
  112 format(//' n, xabc1(n), varc1(2,n), varc1(7,n), varc1(15,n):'/
     *  (i4,1p4e14.6))
  113 format(//' n, xabc1(n), varc1(1,n), varc1(7,n), varc1(15,n):'/
     *  (i4,1p4e14.6))
  115 format(a)
  116 format(4i10)
  117 format(5e16.9)
  118 format(' data set ',i2,'  record ',i6)
  120 format('#'/'#  model no. ',i2)
  125 format('#'/'#'/'# head:'/('#',a79))
  126 format(//' head, model no.',i2,':'/(a79))
  127 format(' number of points =',i6,'   datout(1-10):'/(1p5e15.6))
  130 format('#'/'#'/'# datout:'/('#',1p5e15.6))
  140 format(' varout(i,n,2):'/(1p5e13.5))
  142 format(' yint(i):'/(1p5e13.5))
  145 format(/' ***** Warning. Extrapolation disallowed from x =',
     *  1pe13.5)
  150 format('#'/'#'/'# output format:'/
     *  '# r/R(1), m/M, dy(1-',i2,'),       where'/'#'/
     *  ('# dy(',i2,'): ',a20))
  155 format('#')
  160 format(0pf12.6,1pe13.6,20e11.3)
  170 format(//' extreme differences.'/
     *  '   variable                  xextreme          extremes:'/
     *  '                                            min        max'//
     *  (i4,1x,a20,0pf12.6,1p2e13.5))
      end
