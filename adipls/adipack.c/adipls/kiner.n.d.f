      subroutine kiner(x,y,aa,el,iy,ia,nn,ekin)
c  calculates kinetic energy integral
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
c  Modified 20/1/95 for resetting of treatment of turbulent pressure
c
c  Modified 5/8/95, to allow normalizing energy by squared norm of
c  total displacement, evaluated at phosotsphere, for iekinr = 1
c  For iekinr = 0 old normalization (with surface vertical displacement)
c  is used. iekinr is transmitted in common/cincnt/.
c
      implicit double precision (a-h,o-z)
      dimension x(nn),y(iy,nn),aa(ia,nn)
      dimension yi(2)
      common/worksp/ wrk(2,1)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
c
      save
c
      pi = 4.d0*atan(1.d0)
c
c  set integrand
c  test for plane-parallel case
      if(iplneq.eq.1) then
        elli=1.d0/(el*el)
        do 5 n=1,nn
        y1=y(1,n)
        y2=y(2,n)
c  avoid underflows
        if(abs(y1).lt.eprufl) y1=0
        if(abs(y2).lt.eprufl) y2=0
    5   wrk(1,n)=(y1*y1+elli*y2*y2)*aa(5,n)
c
      else if(el.gt.1.e-6) then
c
c  nonradial case
c
        elli=1.d0/(el*(el+1))
        do 10 n=1,nn
        y1=y(1,n)
        y2=y(2,n)
c  avoid underflows
        if(abs(y1).lt.eprufl) y1=0
        if(abs(y2).lt.eprufl) y2=0
        xx=x(n)
   10   wrk(1,n)=(y1*y1+y2*y2*elli)*aa(1,n)*aa(5,n)*aa(10,n)*xx*xx
c
      else
c
c  radial case
c
        do 25 n=1,nn
        xx=x(n)
        y1=y(1,n)
   25   wrk(1,n)=y1*y1*aa(1,n)*aa(5,n)*aa(10,n)*xx*xx
c
      end if
c
c  integrate
c
   30 call vinta(x,wrk,wrk(2,1),nn,2,2)
c
c  normalize result, depending on iekinr
c
      if(iekinr.ne.1) then
        if(y1.eq.0) y1=1
	ynorm=y1*y1
      else
	call lir(1.d0,x,yi,y,2,iy,nn,1,inter)
	if(el.le.1.e-6) then
	  ynorm=yi(1)*yi(1)
        else
	  ynorm=yi(1)*yi(1)+elli*yi(2)*yi(2)
        end if
      end if
      ekin=wrk(2,nn)/(4.d0*pi*ynorm)
      return
      end
c
c
