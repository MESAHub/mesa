      double precision function xlmult(y,x,el,icxm,ierr)
c
c  finds xlmult = y*x**el. when icxm = 0 the simple expression is used.
c  otherwise ln(xlmult) is calculated and checked for under- or over
c  flow. in case of under- or overflow ierr is set to 1 and xlmult to 0.
c  if x is negative and el is not an integer, xlmult is set to
c  y*abs(x)**el and ierr is set to -1.
c  if x = 0 and el .lt. 0, ierr is set to -2, and xlmult is
c  set to 0.
c  if no problems are encountered ierr is returned as 0.
c
c  modified on 25/1/1985.
c
c  limit of ln(real number)
c
c  ....................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      logical samex
      save
c
      data alimln /85.d0/
c     ********************
      data xp,elp,axl,icxp /0.d0,0.d0,0.d0,-1/
c
c  test for zero y
c
      if(y.eq.0) then
        xlm=0
        ierr=0
        go to 80
      end if
c
c  test for same x and el as in previous call
c
      samex=x.eq.xp.and.el.eq.elp.and.icxm.eq.icxp
      xp=x
      elp=el
      icxp=icxm
c
      ierr=0
c  sign of x
      if(x) 10,14,16
c
c  negative x. integer or non-integer el?
c
   10 ax=-x
      iel=el
      if(iel.eq.el) go to 12
      isx=1
      ierr=-1
      go to 20
c
   12 isx=1
      if(mod(iel,2).ne.0) isx=-1
      go to 20
c
c  x = 0.
c
   14 if(el) 14100,14200,14300
c
c  el .lt. 0, set error flag
c
14100 ierr=-2
      xlm=0
      go to 80
c
c  el = 0
c
14200 xlm=y
      go to 80
c
c  el .gt. 0
c
14300 xlm=0
      go to 80
c
c  positive x
c
   16 isx=1
      ax=x
c
   20 if(icxm.eq.0) go to 35
c  sign of result
      if(y) 22,24,26
   22 isx=-isx
      ay=-y
      go to 30
c  y = 0. set xlmult to 0.
   24 xlm=0
      go to 80
c
   26 ay=y
c  set ln(xlmult)
   30 if(samex) go to 32
      axl=el*log(ax)
   32 xml=log(ay)+axl
c  test for under- or overflow
      if(abs(xml).gt.alimln) go to 40
      xlm=isx*exp(xml)
      go to 80
c  moderately sized el
   35 if(samex) go to 37
      axl=ax**el
   37 xlm=isx*y*axl
      go to 80
c
c  under- or overflow
c
   40 xlm=0
      ierr=1
c
   80 continue
      xlmult=xlm
c
      return
      end
c
c
