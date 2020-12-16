C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: molecular_hydrogen.f 352 2006-04-13 02:12:47Z airwin $
C
C       For the latest version of this source code, please contact
C       Alan W. Irwin
C       Department of Physics and Astronomy
C       University of Victoria,
C       Box 3055
C       Victoria, B.C., Canada
C       V8W 3P6
C       e-mail: irwin@beluga.phys.uvic.ca.
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C
C       End of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C*******************************************************************************
      subroutine molecular_hydrogen(ifh2, ifh2plus, tl,
     &  qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt)
C       calculate partition function eos quantities for molecular hydogen and 
C       its positive ion.
C       ifh2 = 0  no h2
C       ifh2 = 1  vdb h2
C       ifh2 = 2  st h2
C       ifh2 = 3  irwin h2 (recommended)
C       ifh2 = 4  pteh h2
C       ifh2plus = 0  no h2plus
C       ifh2plus = 1  st h2plus
C       ifh2plus = 2  irwin h2plus (recommended)
C       qh2 is the ln partition function
C       qh2t is dlnq/dlnt
C       qh2tt is d2lnq/(dlnt)^2
C       qh2plus is the ln partition function
C       qh2plust is dlnq/dlnt
C       qh2plustt is d2lnq/(dlnt)^2
      implicit none
      include 'constants.h'
      integer ifh2, ifh2plus
      double precision tl,
     &  qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt
C       partition function coefficients
      double precision qh2vdb(2)
      data qh2vdb/4.874d-3, 1.371d0/
!  parameter(nqh2i87 = 8)
!  dimension qh2i87(nqh2i87), dqh2i87(nqh2i87-1), d2qh2i87(nqh2i87-2)
!  data qh2i87/
!  1 1.69179D+00, -1.7227D+00, 7.98033D-01, -1.57089D-01, -5.35313D-01,
!  2 1.75818D+00, -2.63895D+00, 1.35708D+00/
!  parameter(nqh2plusi89 = 10)
!  dimension qh2plusi89(nqh2plusi89), dqh2plusi89(nqh2plusi89-1),
!  1 d2qh2plusi89(nqh2plusi89-2)
!  data qh2plusi89/
!  1 2.564408220D+00, -2.15158152D+00, 4.046884D-01, 1.882055D+00,
!  2 -2.439623D+00, -3.47104D+00, 9.45846D+00, 6.971D-02,
!  3 -1.34950D+01, 8.59387D+00/
      integer nqh2i90
      parameter(nqh2i90 = 9)
      double precision qh2i90(nqh2i90), dqh2i90(nqh2i90-1), d2qh2i90(nqh2i90-2)
      data qh2i90/
     &  1.6918292822D+00,-1.72246845D+00,7.9631758D-01,
     &  -1.706903D-01,-4.574240D-01,1.825633D+00,
     &  -3.52468D+00,2.91102D+00,-8.40929D-01/
      integer nqh2plusi90
      parameter(nqh2plusi90 = 11)
      double precision qh2plusi90(nqh2plusi90), dqh2plusi90(nqh2plusi90-1),
     &  d2qh2plusi90(nqh2plusi90-2)
      data qh2plusi90/
     &  2.5644276281D+00,-2.152026894D+00,4.0264504D-01,
     &  1.9297089D+00,-2.432017D+00,-4.784831D+00,
     &  1.208765D+01,8.31227D+00,-4.94136D+01,5.46883D+01,-2.03073D+01/
      integer nqh2i90_hi
      parameter(nqh2i90_hi = 9)
      double precision qh2i90_hi(nqh2i90_hi), dqh2i90_hi(nqh2i90_hi-1),
     &  d2qh2i90_hi(nqh2i90_hi-2)
      data qh2i90_hi/
     &  1.6936496106D+00,-1.785296204D+00,-3.3412934D-01,
     &  -7.18046192D+00,-2.13716663D+01,-2.73213417D+01,
     &  -1.86763490D+01,-6.74105728D+00,-1.01516624D+00/
      integer nqh2plusi90_hi
      parameter(nqh2plusi90_hi = 8)
      double precision qh2plusi90_hi(nqh2plusi90_hi),
     &  dqh2plusi90_hi(nqh2plusi90_hi-1),
     &  d2qh2plusi90_hi(nqh2plusi90_hi-2)
      data qh2plusi90_hi/
     &  2.56934142D+00,-2.0517792D+00,1.239529D+00,
     &  5.653001D+00,6.798597D+00,4.355492D+00,
     &  1.508621D+00,2.226748D-01/
      integer nqh2st
      parameter(nqh2st = 4)
      double precision qh2st(nqh2st), dqh2st(nqh2st-1), d2qh2st(nqh2st-2)
      data qh2st/
     &  1.6498d0, -1.6265d0, 0.7472d0, -0.2751d0/
      integer nqh2plusst
      parameter(nqh2plusst = 5)
      double precision qh2plusst(nqh2plusst), dqh2plusst(nqh2plusst-1),
     &  d2qh2plusst(nqh2plusst-2)
      data qh2plusst/
     &  2.5410d0, -2.4336d0, 1.4979d0, 0.0192d0, -0.7483d0/
      integer nqh2pteh
      parameter (nqh2pteh = 5)
      double precision qh2pteh(nqh2pteh)
C       from PTEH paper with DH2 adjusted from 4.48 to 4.477 ev (Pols, private
C       communication, 1996).
      data qh2pteh /
     &  6608.8d0, 0.448d0, 0.1562d0, 0.0851d0, 4.477d0/
      integer *4 iffirst, ifprint
      data iffirst/1/  !flag for first entry into routine.
C       flag for one printout for h2 or h2+ extrapolation beyond one dex
C       above range.
      data ifprint/1/
      integer iorder, ifextrapolate, ifextrapolate_dex
      double precision rmod, t, dtl, ti, tc2, thetal, eh2, poly_sum,
     &  argzeta, zeta, dzeta, dzeta2, tmin, tmax
      parameter (tmin=0.9999999999999d3,
     &  tmax=1.000000000000001d5)
      save
      if(iffirst.eq.1) then
        iffirst = 0
        rmod = log(10.d0)
        do iorder = 1,nqh2i90-1
          dqh2i90(iorder) = dble(iorder)*qh2i90(iorder+1)
        enddo
        do iorder = 1,nqh2i90-2
          d2qh2i90(iorder) = dble(iorder)*dqh2i90(iorder+1)
        enddo
        do iorder = 1,nqh2plusi90-1
          dqh2plusi90(iorder) = dble(iorder)*qh2plusi90(iorder+1)
        enddo
        do iorder = 1,nqh2plusi90-2
          d2qh2plusi90(iorder) = dble(iorder)*dqh2plusi90(iorder+1)
        enddo
        do iorder = 1,nqh2i90_hi-1
          dqh2i90_hi(iorder) = dble(iorder)*qh2i90_hi(iorder+1)
        enddo
        do iorder = 1,nqh2i90_hi-2
          d2qh2i90_hi(iorder) = dble(iorder)*dqh2i90_hi(iorder+1)
        enddo
        do iorder = 1,nqh2plusi90_hi-1
          dqh2plusi90_hi(iorder) = dble(iorder)*qh2plusi90_hi(iorder+1)
        enddo
        do iorder = 1,nqh2plusi90_hi-2
          d2qh2plusi90_hi(iorder) =
     &      dble(iorder)*dqh2plusi90_hi(iorder+1)
        enddo
        do iorder = 1,nqh2st-1
          dqh2st(iorder) = dble(iorder)*qh2st(iorder+1)
        enddo
        do iorder = 1,nqh2st-2
          d2qh2st(iorder) = dble(iorder)*dqh2st(iorder+1)
        enddo
        do iorder = 1,nqh2plusst-1
          dqh2plusst(iorder) = dble(iorder)*qh2plusst(iorder+1)
        enddo
        do iorder = 1,nqh2plusst-2
          d2qh2plusst(iorder) = dble(iorder)*dqh2plusst(iorder+1)
        enddo
        qh2pteh(1) = log(qh2pteh(1))
        do iorder = 2, nqh2pteh
          qh2pteh(iorder) = (ergsperev/boltzmann)*qh2pteh(iorder)
        enddo
      endif
C       n.b. this is a local t and should not affect anything outside this
C       routine.
      t=exp(tl)
      if(ifh2.lt.0.or.ifh2.gt.4)
     &  stop 'molecular hydrogen: ifh2 must be in range from 0 to 4'
      if(.not.(ifh2plus.eq.ifh2-1.or.ifh2plus.eq.0.or.
     &  (ifh2plus.eq.2.and.ifh2.eq.4)))
     &  stop 'molecular_hydrogen: bad ifh2plus for input ifh2'
C       warning: the ifh2=1 and 2 options are simply compatibility modes
C       with old programmes.  Some effort has been made to keep these
C       options, current, but they are no longer tested, and they might
C       not work.
      if(ifh2.eq.1) then
        if(t.lt.tmin.or.t.gt.8.d3) then
          ifextrapolate = 1
          t=min(8.d3,max(tmin,t))
        else
          ifextrapolate = 0
        endif
        ifextrapolate_dex = ifextrapolate
      elseif(2.le.ifh2.and.ifh2.le.4) then
        if(ifh2.gt.2) then
C           H2 and H2+ are negligible at 1.d6, but second-order
C           extrapolation is ill-behaved beyond this temperature 
C           (underflows) so warn.
C           n.b. get this test done *before* t is modified below.
          if(t.lt.tmin.or.t.gt.10.d0*tmax) then
            ifextrapolate_dex = 1
          else
            ifextrapolate_dex = 0
          endif
        endif
        if(t.lt.tmin.or.t.gt.tmax) then
          ifextrapolate = 1
          t=min(tmax,max(tmin,t))
        else
          ifextrapolate = 0
        endif
        if(ifh2.eq.2) then
C           this mode has tremendous errors at T = 1.d5 = tmax
C           (original polynomial fit only to 9000 K) so warn.
C           n.b. this test must be made *after* ifextrapolate is
C           calculated.
          ifextrapolate_dex = ifextrapolate
        endif
      endif
      thetal = dlog10(5040.d0/t)
      ti=11605.5d0/t
      tc2 = c2/t
      if(ifh2.le.0.and.ifh2plus.gt.0)
     &  stop 'molecular_hydrogen: invalid ifh2, ifh2plus switches'
      if(ifh2.eq.1) then
C         vdb h2 partition function.
        eh2 = qh2vdb(2)/ti
        qh2 = log(qh2vdb(1)*t) + eh2
        qh2t = 1.d0+eh2
        qh2tt = eh2
      elseif(ifh2.eq.2) then
C         st h2 partition function.
        qh2 = poly_sum(thetal,qh2st,nqh2st)*rmod
        qh2t = -poly_sum(thetal,dqh2st,nqh2st-1)
        qh2tt = poly_sum(thetal,d2qh2st,nqh2st-2)/rmod
      elseif(ifh2.eq.3) then
C         irwin (1990) h2 partition function (published in 1996)
        if(t.lt.9.d3) then
          qh2 = poly_sum(thetal,qh2i90,nqh2i90)*rmod
          qh2t = -poly_sum(thetal,dqh2i90,nqh2i90-1)
          qh2tt = poly_sum(thetal,d2qh2i90,nqh2i90-2)/rmod
        else
          qh2 = poly_sum(thetal,qh2i90_hi,nqh2i90_hi)*rmod
          qh2t = -poly_sum(thetal,dqh2i90_hi,nqh2i90_hi-1)
          qh2tt = poly_sum(thetal,d2qh2i90_hi,nqh2i90_hi-2)/rmod
        endif
      elseif(ifh2.eq.4) then
        argzeta = qh2pteh(5)/t
        zeta = 1.d0 - exp(-argzeta)*(1.d0 + argzeta)
        qh2 = qh2pteh(1) + (qh2pteh(2) - (qh2pteh(3)*qh2pteh(3) -
     &    qh2pteh(4)*qh2pteh(4)*qh2pteh(4)/t)/t)/t -
     &    2.5d0*log(qh2pteh(5)/t) + log(zeta)
        dzeta = -exp(-argzeta)*argzeta*argzeta
        qh2t = -(qh2pteh(2) - (2.d0*qh2pteh(3)*qh2pteh(3) -
     &    3.d0*qh2pteh(4)*qh2pteh(4)*qh2pteh(4)/t)/t)/t + 2.5d0 +
     &    dzeta/zeta
        dzeta2 = -exp(-argzeta)*argzeta*argzeta*(argzeta-2.d0)
        qh2tt = (qh2pteh(2) - (4.d0*qh2pteh(3)*qh2pteh(3) -
     &    9.d0*qh2pteh(4)*qh2pteh(4)*qh2pteh(4)/t)/t)/t +
     &    (zeta*dzeta2 - dzeta*dzeta)/(zeta*zeta)
      endif
      if(ifh2.gt.0.and.ifextrapolate.eq.1) then
        if(ifh2.eq.3.or.ifh2.eq.4) then
C           Taylor series approach for ln q
C           n.b. t has been adjusted to maximum or minimum
          dtl = tl - log(t)
          qh2 = qh2 + dtl*qh2t + 0.5d0*dtl*dtl*qh2tt
          qh2t = qh2t + dtl*qh2tt
        else
C           don't bother with old crummy partition functions
          qh2t = 0.d0
          qh2tt = 0.d0
        endif
      endif
      if(ifh2.gt.0.and.ifextrapolate_dex.eq.1) then
        if(ifprint.eq.1) then
          ifprint = 0
!          write(0,'(a)') 
!     &      ' h2 p.f.extrapolated 1 dex beyond valid T range'
!          if(ifh2plus.gt.0) write(0,'(a)') 
!     &      ' h2+ p.f. extrapolated 1 dex beyond valid T range'
        endif
      endif
      if(ifh2plus.eq.1) then
C         st h2plus partition function.
        qh2plus = poly_sum(thetal,qh2plusst,nqh2plusst)*rmod
        qh2plust = -poly_sum(thetal,dqh2plusst,nqh2plusst-1)
        qh2plustt = poly_sum(thetal,d2qh2plusst,nqh2plusst-2)/rmod
      elseif(ifh2plus.eq.2) then
C         irwin (1990) h2plus partition function. (published in 1996)
        if(t.lt.9.d3) then
          qh2plus = poly_sum(thetal,qh2plusi90,nqh2plusi90)*rmod
          qh2plust = -poly_sum(thetal,dqh2plusi90,nqh2plusi90-1)
          qh2plustt = poly_sum(thetal,d2qh2plusi90,nqh2plusi90-2)/rmod
        else
          qh2plus = poly_sum(thetal,qh2plusi90_hi,nqh2plusi90_hi)*rmod
          qh2plust = -poly_sum(thetal,dqh2plusi90_hi,nqh2plusi90_hi-1)
          qh2plustt =
     &      poly_sum(thetal,d2qh2plusi90_hi,nqh2plusi90_hi-2)/rmod
        endif
      endif
      if(ifh2plus.gt.0.and.ifextrapolate.eq.1) then
        if(ifh2plus.eq.2) then
C           Taylor series approach for ln q
C           n.b. t has been adjusted to maximum or minimum
          dtl = tl - log(t)
          qh2plus = qh2plus + dtl*qh2plust +
     &      0.5d0*dtl*dtl*qh2plustt
          qh2plust = qh2plust + dtl*qh2plustt
        else
C           don't bother with old crummy partition functions
          qh2plust = 0.d0
          qh2plustt = 0.d0
        endif
      endif
      end
