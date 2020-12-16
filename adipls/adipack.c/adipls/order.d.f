      subroutine order(x,y,data,aa,el,sig,icow,irsord,
     *  iy,iaa,iord,nn,iasn,mdintg,iprint)
c
c  calculates order (iord) of the oscillation according to scuflaire's
c  classification.
c
c  when irsord ge 1 and le 10, and icow .le. 1 order for l = 1 
c  and scuflaire order between  0 and irsord is incremented by 1,
c  to correct for problem with order. 
c  For traditional models of the present sun irsord = 1 is appropriate
c
c  if iabs(irsord) = 11, use instead classification 
c  proposed by Lee (PASJ, vol. 37, p. 279, 1985), based on 
c  -(delta Phi) and p'
c
c  if irsord = 20 use the Takata scheme (Takata,M. 2006, ESASP, 624, 26)
c
c  for irsord .lt. 0 evaluate both classifications and scream if they 
c  differ
c
c  if iprint = 1, print result in standard form
c  if iprint = -1, only print result in standard form
c
c  modified 11/7/1985 to use same calling sequence as in cambridge,
c  with test on icow.
c
c  modified 13/8/87 to standardize output
c
c  modified 16/7/90 to include option of not printing result, and for
c  allowing calling order on submesh, during Richardson extrapolation
c
c  modified 5/7/91, to include Lee classification
c
c  Modified 1/7/95, to do resetting also if icow = 1
c
c  Modified 5/11/02, including option for azimuthal order
c
c  Modified 18/3/05, changing Lee scheme to use y(2,.) rather than p'.
c  Also, interpolate y(2.) to zero of y(1,.) to determine appropriate
c  sign. Introduce irsord = 12, to use Lee scheme only for l = 1.
c
c  Modified 5/8/05, correcting setting of initial phase when
c  Lee scheme is used (still deserves a little checking!)
c
c  Modified 27/2/08 adding Gulnur's version of the Takata scheme.
c
c  ...................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      character*1 mtype
      dimension x(*),y(iy,*),data(*),aa(iaa,*),iordst(2)
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
      common/rotdat/ em, irotsl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      idiag=0
c
c  for iprint = -1, only print order
c
      if(iprint.eq.-1) go to 30
c
c  for mdintg = 2, iord is already set
c
      if(mdintg.eq.2) go to 25
c
c  test for using Takata scheme
c
      if(irsord.eq.20.and.abs(el-1.d0).le.1.d-10) then
        call takata(x,y,nn,iasn,data,aa,iy,iaa,iord,sig)
        if(istdpr.gt.0) write(istdpr,'(/'' Use Takata/Dogan order'')')
        go to 30
      end if
c
c  set control parameter for using Lee method (only applicable for
c  full nonradial case)
c  For irsord = 12 (or -12) only enforce Lee scheme for l = 1
c
      if((abs(el-1.d0).gt.1.d-10.and.
     *  (irsord.le.-12.or.irsord.ge.12)).or.icow.ne.0) then
	irsor1=0
      else
	irsor1=irsord
      end if
      if(irsor1.gt.10.and.istdpr.gt.0) 
     *  write(istdpr,*) 'Using Lee scheme'
c
      els=el*(el+1)/sig
c
c  set limit for underflow
c
      phxmin=eprufl*max(abs(y(1,nn)),abs(y(3,nn)))
c
c  test and set actual number of meshpoints, depending on whether model
c  extends to centre and whether surface is singular
c
      if(x(1).eq.0) then
        ns=2
      else
        ns=1
      end if
c
      if(data(7).lt.0) then
	nsr=nn
      else
	nsr=nn-1
      end if
c
      if(mod(nsr-ns,iasn).ne.0) then
        write(istdou,100) nn, nsr, ns, iasn
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,100) 
     *    nn, nsr, ns, iasn
        return
      end if
c
c  step up in n to find first point that avoids underflow
c
    2 if(abs(y(1,ns)).le.phxmin) then
        ns=ns+iasn
        go to 2
      end if
c
c  test for number of cases to be considered
c
      if(irsor1.lt.-10) then
	ncase=2
	iphx=2
	iphy=2
      else if(irsor1.gt.10) then
	ncase=1
	iphx=2
	iphy=2
      else
	ncase=1
	iphx=1
	iphy=1
      end if
c
c  find contribution to order from innermost point, depending on
c  quadrant
c
      do 22 icase=1,ncase
c
      if(iphx.eq.1) then
        phx2=y(1,ns)
      else
        phx2=y(1,ns) - y(3,ns)
      end if
c
      if(iphy.eq.1) then
        phy2=y(2,ns)
      else
c..        phy2=y(2,ns) + els*aa(1,ns)*y(3,ns)
        phy2=y(2,ns) 
      end if
      if(idiag.gt.0.and.istdpr.gt.0) then
        write(istdpr,'('' centre, x, phx, phy ='',i5,0pf12.7,1p2e11.3)')
     *      ns, x(ns), phx2, phy2
      end if
c
      if(phx2*phy2.lt.0.and.irsor1.le.10) then
	iordp=1
      else
        iordp=0
      end if
      iordg=0
      if(idiag.gt.0.and.istdpr.gt.0) 
     *  write(istdpr,*) 'Initial iordp =',iordp
c
c  step through solution
c
      ns=ns+iasn
c
      do 20 n=ns,nsr,iasn
      phx1=phx2
      phy1=phy2
c
      if(iphx.eq.1) then
        phx2=y(1,n)
      else
        phx2=y(1,n) - y(3,n)
      end if
c
      if(iphy.eq.1) then
        phy2=y(2,n)
      else
c..        phy2=y(2,n) + els*aa(1,n)*y(3,n)
        phy2=y(2,n) 
      end if
c
      if(abs(phx2).ge.phxmin.and.phx1*phx2.le.0) then
	phy2m=(phx1*phy2-phx2*phy1)/(phx1-phx2)
        dfyy=phy2m*(phx2-phx1)
        if(dfyy.lt.0) then
          iordp=iordp+1
        else if(dfyy.gt.0) then
          iordg=iordg+1
	end if
	if(idiag.gt.0.and.istdpr.gt.0) then
          write(istdpr,'('' n, x, phx, phy ='',i5,0pf12.7,1p2e11.3,
     *      '' new iordg, iordp ='',2i5)')
     *      n, x(n), phx2, phy2m, iordg, iordp
	end if
      end if
   20 continue
c
      if(iphx.eq.2.and.iordp.ge.iordg) then
	iord=iordp-iordg+1
      else
	iord=iordp-iordg
      end if
      iordst(icase)=iord
      iphx=1
      iphy=1
   22 continue
c
c  test for comparisons of methods
c
      if(ncase.eq.2.and.iordst(1).ne.iordst(2).and.istdpr.gt.0) then
	write(istdpr,102) iordst(2), iordst(1)
      end if
c
      iord=iordst(1)
c
c  test for resetting order
c
   25 if(irsor1.gt.0.and.irsor1.le.10.and.el.eq.1.and.icow.le.1.and.
     *  (iord.ge.0.and.iord.le.irsor1)) then
        iord=iord+1
        if(istdpr.gt.0) write(istdpr,130) iord,el
      end if
c
c  print order
   30 ial=el
      if(irotsl.eq.1) iam=nint(em)
      if((iprint.ne.1.and.iprint.ne.-1).or.istdpr.le.0) then
	return
      else if(iord.eq.0) then
c  f mode
	if(irotsl.ne.1) then
          write(istdpr,110) ial
	else
          write(istdpr,112) ial,iam
	end if
        return
      else if(iord.lt.0) then
c  g mode
        iord1=-iord
        mtype='g'
      else
c  p mode
        mtype='p'
	iord1=iord
      end if
c
      if(iord1.lt.100) then
	if(irotsl.ne.1) then
          write(istdpr,120) mtype,iord1,ial
	else
          write(istdpr,122) mtype,iord1,ial,iam
	end if
        return
      else
	if(irotsl.ne.1) then
          write(istdpr,124) mtype,iord1,ial
	else
          write(istdpr,126) mtype,iord1,ial,iam
	end if
        return
      end if
  100 format(/' **** error in s/r order.',
     *  ' number of steps is not integral multiplum of iasn'/
     *  '            nn, nsr, nc, iasn =',4i6)
  102 format(/' **** warning in s/r order. original order =',
     *  i4,' Lee order =',i4)
  110 format(///' this is an f(l =',i4,') mode')
  112 format(///' this is an f(l =',i4,', m =',i4,') mode')
  120 format(///' this is a ',a1,i2,'(l =',i4,') mode')
  122 format(///' this is a ',a1,i2,'(l =',i4,', m =',i4,') mode')
  124 format(///' this is a ',a1,i4,'(l =',i4,') mode')
  126 format(///' this is a ',a1,i4,'(l =',i4,', m =',i4,') mode')
  130 format(//' ***** in s/r order, order has been incremented',
     *  ' by 1 to',i4,' for l =',f6.1)
      end
