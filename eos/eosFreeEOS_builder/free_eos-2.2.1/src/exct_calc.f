C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: exct_calc.f 815 2008-06-24 18:10:38Z airwin $
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
C      Note:
C      The original version of this routine is in the public domain.  It was
C      written in the late 60's or early 70's by Darwin Harwood for the US
C      government, and it was brought to the attention of AWI by Forrest J.
C      Rogers and Fritz J. Swenson. AWI relicensed the routine under the GPL
C      and put in many changes to convert to structured programming and
C      double precision, evaluate partial derivatives, and insert commentary.
C
      subroutine exct_calc(fp12, exct, dexct, d2exct, d3exct)
C       the subroutine exct_calc which returns exct and its first through
C       third derivatives wrt fp12.
C
C       exct is fp12^-2* integral from -inf to eta of square of 
C       derivative of fp12(eta') d eta'.  
C       N.B. For numerical
C       convenience, the integral is transformed so that the independent 
C       variable (upper limit) is fp12(eta) rather than eta, and this
C       is accounted for in the analytic verification below.

C       relative error less than 4.e-7 according to original notes.
C
C       we have numerically verified derivatives, and we have 
C       verified above definition by checking it for large negative eta
C       and by numerically verifying analytic derivative of above relationship
C       to 1.d-4 when eta > -20.
      implicit none
      double precision fp12, x, y, dy, d2y, d3y, z, exct, dexct, d2exct, d3exct
      x=fp12
      if(x.lt.0.d0) then
        write(*,*) 'x = ', x
        stop 'exct: bad x value'
      elseif(x.le.5.d0) then
        exct=.5d0+x*(-.117847745d0+x*(.033695787d0+x*(-.0102270934d0+x*
     &    (.0031103379d0+x*(-8.743415d-4+x*(2.0565698d-4+
     &    x*(-3.6687731d-5+x*(4.4980746d-6+x*(-3.3176255d-7+
     &    x*1.0992982d-8)))))))))
        dexct= -.117847745d0+x*(2.d0*.033695787d0+
     &    x*(3.d0*(-.0102270934d0)+x*(4.d0*.0031103379d0+
     &    x*(5.d0*(-8.743415d-4)+x*(6.d0*2.0565698d-4+
     &    x*(7.d0*(-3.6687731d-5)+x*(8.d0*4.4980746d-6+
     &    x*(9.d0*(-3.3176255d-7)+x*10.d0*1.0992982d-8))))))))
        d2exct= 2.d0*.033695787d0+
     &    x*(6.d0*(-.0102270934d0)+x*(12.d0*.0031103379d0+
     &    x*(20.d0*(-8.743415d-4)+x*(30.d0*2.0565698d-4+
     &    x*(42.d0*(-3.6687731d-5)+x*(56.d0*4.4980746d-6+
     &    x*(72.d0*(-3.3176255d-7)+x*90.d0*1.0992982d-8)))))))
        d3exct= 6.d0*(-.0102270934d0)+x*(24.d0*.0031103379d0+
     &    x*(60.d0*(-8.743415d-4)+x*(120.d0*2.0565698d-4+
     &    x*(210.d0*(-3.6687731d-5)+x*(336.d0*4.4980746d-6+
     &    x*(504.d0*(-3.3176255d-7)+x*720.d0*1.0992982d-8))))))
      else
        y = x**(-2.d0/3.d0)
        dy = (-2.d0/3.d0)*y/x
        d2y = (-5.d0/3.d0)*dy/x
        d3y = (-8.d0/3.d0)*d2y/x
        z=y*y
        exct=.93052568d0+z*(-.75731605d0+z*(1.2881026d0+
     &    z*(2.9520325d0+z*(660.40937d0+
     &    z*(-26107.983d0+z*(476965.95d0+
     &    z*(-5132815.7d0+z*(33305973.8d0+
     &    z*(-120795398.d0+z*188387439.d0)))))))))
        dexct=-.75731605d0+z*(2.d0*1.2881026d0+
     &    z*(3.d0*2.9520325d0+z*(4.d0*660.40937d0+
     &    z*(5.d0*(-26107.983d0)+z*(6.d0*476965.95d0+
     &    z*(7.d0*(-5132815.7d0)+z*(8.d0*33305973.8d0+
     &    z*(9.d0*(-120795398.d0)+z*10.d0*188387439.d0))))))))
        d2exct=2.d0*1.2881026d0+
     &    z*(6.d0*2.9520325d0+z*(12.d0*660.40937d0+
     &    z*(20.d0*(-26107.983d0)+z*(30.d0*476965.95d0+
     &    z*(42.d0*(-5132815.7d0)+z*(56.d0*33305973.8d0+
     &    z*(72.d0*(-120795398.d0)+z*90.d0*188387439.d0)))))))
        d3exct= 6.d0*2.9520325d0+z*(24.d0*660.40937d0+
     &    z*(60.d0*(-26107.983d0)+z*(120.d0*476965.95d0+
     &    z*(210.d0*(-5132815.7d0)+z*(336.d0*33305973.8d0+
     &    z*(504.d0*(-120795398.d0)+z*720.d0*188387439.d0))))))
C         transform from function of z to y
        d3exct = 2.d0*(6.d0*y*d2exct + 4.d0*z*y*d3exct)
        d2exct = 2.d0*(dexct + 2.d0*z*d2exct)
        dexct = 2.d0*y*dexct
C         transform to function of x
        d3exct = d3y*dexct + 3.d0*dy*d2y*d2exct + dy*dy*dy*d3exct
        d2exct = d2y*dexct + dy*dy*d2exct
        dexct = dy*dexct
C         transform to new function of x
        d3exct = d3y*exct + 3.d0*d2y*dexct +
     &    3.d0*dy*d2exct + y*d3exct -
     &    0.6981317008d0*(26.d0 - 24.d0*log(x))/(x*x*x*x*x)
        d2exct = d2y*exct + 2.d0*dy*dexct + y*d2exct -
     &    0.6981317008d0*(-5.d0 + 6.d0*log(x))/(x*x*x*x)
        dexct = dy*exct + y*dexct -
     &    0.6981317008d0*(1.d0 - 2.d0*log(x))/(x*x*x)
        exct = y*exct - 0.6981317008d0*log(x)/(x*x)
      endif
      end
