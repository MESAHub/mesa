C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fp12_calc.f 815 2008-06-24 18:10:38Z airwin $
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
      subroutine fp12_calc(eta, fp12, dfp12, d2fp12, d3fp12)
C************
C called by main
C  no calls
C************
C       the subroutine fp12_calc which returns fp12 and its first,
C       second, and third derivatives wrt eta.
C
C       fp12 = Cox and Giuli F(1/2,eta)/Gamma(3/2)
C
C       good to 5.e-8 according to original notes.
C
C       also checked above definition against analytic derivative 
C       of Cody-Thacher routine for F(3/2,eta) (used in stellar 
C       atmospheres EOS), and obtained agreement to 2.d-7.  
C       Not bad for a derivative of a fit.
      implicit none
      double precision eta, q, y, dy, d2y, d3y, x, r, fp12,
     &  dfp12, d2fp12, d3fp12
      q = eta
      if (q.lt.3.d0) then
        y = exp(q)
        if(q.le.-1.9375d0) then
          fp12 = y*(1.d0-y*(.35355283d0-y*(.19242767d0-y*(.12456909d0-
     &      y*(.085114507d0-y*.04551794d0)))))
          dfp12 = y*(1.d0-y*(2.d0*.35355283d0-y*(3.d0*.19242767d0-
     &      y*(4.d0*.12456909d0-
     &      y*(5.d0*.085114507d0-y*6.d0*.04551794d0)))))
          d2fp12 = y*(1.d0-y*(4.d0*.35355283d0-y*(9.d0*.19242767d0-
     &      y*(16.d0*.12456909d0-
     &      y*(25.d0*.085114507d0-y*36.d0*.04551794d0)))))
          d3fp12 = y*(1.d0-y*(8.d0*.35355283d0-y*(27.d0*.19242767d0-
     &      y*(64.d0*.12456909d0-
     &      y*(125.d0*.085114507d0-y*216.d0*.04551794d0)))))
        else
          x = q-.5d0
          fp12 = y*(.677695804d0-x*(.187773135d0+x*(2.16197521d-2-
     &      x*(9.23703807d-3+x*(1.71735167d-3-x*(6.07913775d-4+
     &      x*(1.1448629d-4-x*(4.5444320d-5+x*(6.4719368d-6-
     &      x*(3.7949830d-6+x*(1.7338029d-7-x*(3.5546516d-7-
     &      x*(3.7329191d-8+x*(3.3097822d-8-x*(8.3190193d-9+
     &      x*(2.2752769d-9-x*(7.836005d-10+x*(7.519551d-11-
     &      x*2.960006d-11))))))))))))))))))
          dfp12 = fp12 +
     &      y*(-(.187773135d0+x*(2.d0*2.16197521d-2-
     &      x*(3.d0*9.23703807d-3+x*(4.d0*1.71735167d-3-
     &      x*(5.d0*6.07913775d-4+
     &      x*(6.d0*1.1448629d-4-x*(7.d0*4.5444320d-5+
     &      x*(8.d0*6.4719368d-6-
     &      x*(9.d0*3.7949830d-6+x*(10.d0*1.7338029d-7-
     &      x*(11.d0*3.5546516d-7-
     &      x*(12.d0*3.7329191d-8+x*(13.d0*3.3097822d-8-
     &      x*(14.d0*8.3190193d-9+
     &      x*(15.d0*2.2752769d-9-x*(16.d0*7.836005d-10+
     &      x*(17.d0*7.519551d-11-
     &      x*18.d0*2.960006d-11))))))))))))))))))
          d2fp12 = 2.d0*dfp12 - fp12 +
     &      y*(-(2.d0*2.16197521d-2-
     &      x*(6.d0*9.23703807d-3+x*(12.d0*1.71735167d-3-
     &      x*(20.d0*6.07913775d-4+
     &      x*(30.d0*1.1448629d-4-x*(42.d0*4.5444320d-5+
     &      x*(56.d0*6.4719368d-6-
     &      x*(72.d0*3.7949830d-6+x*(90.d0*1.7338029d-7-
     &      x*(110.d0*3.5546516d-7-
     &      x*(132.d0*3.7329191d-8+x*(156.d0*3.3097822d-8-
     &      x*(182.d0*8.3190193d-9+
     &      x*(210.d0*2.2752769d-9-x*(240.d0*7.836005d-10+
     &      x*(272.d0*7.519551d-11-
     &      x*306.d0*2.960006d-11)))))))))))))))))
          d3fp12 = 3.d0*(d2fp12 - dfp12) + fp12 +
     &      y*( (6.d0*9.23703807d-3+x*(24.d0*1.71735167d-3-
     &      x*(60.d0*6.07913775d-4+
     &      x*(120.d0*1.1448629d-4-x*(210.d0*4.5444320d-5+
     &      x*(336.d0*6.4719368d-6-
     &      x*(504.d0*3.7949830d-6+x*(720.d0*1.7338029d-7-
     &      x*(990.d0*3.5546516d-7-
     &      x*(1320.d0*3.7329191d-8+x*(1716.d0*3.3097822d-8-
     &      x*(2184.d0*8.3190193d-9+
     &      x*(2730.d0*2.2752769d-9-x*(3360.d0*7.836005d-10+
     &      x*(4080.d0*7.519551d-11-
     &      x*4896.d0*2.960006d-11))))))))))))))))
        endif
      elseif(q.lt.10.d0) then
        x = q-6.5d0
        fp12 = 12.839811d0+x*(2.844774d0+x*(.114920926d0-
     &    x*(3.43733039d-3-x*(2.3980356d-4-x*(2.0201888d-5-
     &    x*(1.5219883d-6-x*(6.2770524d-8+x*(4.8830336d-9-
     &    x*(2.1031164d-9-x*(5.785753d-10-x*(7.233066d-11 -
     &    x*1.230727d-12)))))))))))
        dfp12 = 2.844774d0+x*(2.d0*.114920926d0-x*(3.d0*3.43733039d-3 -
     &    x*(4.d0*2.3980356d-4-x*(5.d0*2.0201888d-5-
     &    x*(6.d0*1.5219883d-6-x*(7.d0*6.2770524d-8 +
     &    x*(8.d0*4.8830336d-9-x*(9.d0*2.1031164d-9-
     &    x*(10.d0*5.785753d-10-x*(11.d0*7.233066d-11 -
     &    x*12.d0*1.230727d-12))))))))))
        d2fp12 = 2.d0*.114920926d0-x*(6.d0*3.43733039d-3 -
     &    x*(12.d0*2.3980356d-4-x*(20.d0*2.0201888d-5-
     &    x*(30.d0*1.5219883d-6-x*(42.d0*6.2770524d-8 +
     &    x*(56.d0*4.8830336d-9-x*(72.d0*2.1031164d-9-
     &    x*(90.d0*5.785753d-10-x*(110.d0*7.233066d-11 -
     &    x*132.d0*1.230727d-12)))))))))
        d3fp12 = -(6.d0*3.43733039d-3 -
     &    x*(24.d0*2.3980356d-4-x*(60.d0*2.0201888d-5-
     &    x*(120.d0*1.5219883d-6-x*(210.d0*6.2770524d-8 +
     &    x*(336.d0*4.8830336d-9-x*(504.d0*2.1031164d-9-
     &    x*(720.d0*5.785753d-10-x*(990.d0*7.233066d-11 -
     &    x*1320.d0*1.230727d-12)))))))))
      elseif(q.lt.20.d0) then
        x = q-14.5d0
        fp12 = 41.7799227d0+x*
     &    (4.2881461d0+x*(7.45407825d-2-x*(8.79243296d-4 -
     &    x*(2.38288861d-5-x*(8.82474867d-7-x*(3.82865217d-8 -
     &    x*(1.9274292d-9-x*(1.42248669d-10-x*(8.17019813d-12)))))))))
        dfp12 = 4.2881461d0+x*(2.d0*7.45407825d-2-
     &    x*(3.d0*8.79243296d-4-x*(4.d0*2.38288861d-5-
     &    x*(5.d0*8.82474867d-7-x*(6.d0*3.82865217d-8 -
     &    x*(7.d0*1.9274292d-9-x*(8.d0*1.42248669d-10-x*
     &    (9.d0*8.17019813d-12))))))))
        d2fp12 = 2.d0*7.45407825d-2-x*(6.d0*8.79243296d-4 -
     &    x*(12.d0*2.38288861d-5-x*(20.d0*8.82474867d-7-x*
     &    (30.d0*3.82865217d-8 -
     &    x*(42.d0*1.9274292d-9-x*(56.d0*1.42248669d-10-x*
     &    (72.d0*8.17019813d-12)))))))
        d3fp12 = -(6.d0*8.79243296d-4 -
     &    x*(24.d0*2.38288861d-5-x*(60.d0*8.82474867d-7-x*
     &    (120.d0*3.82865217d-8 -
     &    x*(210.d0*1.9274292d-9-x*(336.d0*1.42248669d-10-x*
     &    (504.d0*8.17019813d-12)))))))
      else
        r = 1.d0/q
        fp12 = 1.d0-r*(9.354d-7-r*(1.2338391d0-r*(6.77931d-3-
     &    r*1.17871643d0)))
        dfp12 = 
     &    r*r*(9.354d-7-r*(2.d0*1.2338391d0-r*(3.d0*6.77931d-3-
     &    r*4.d0*1.17871643d0)))
        d2fp12 = 
     &    -r*r*r*(2.d0*9.354d-7-r*(6.d0*1.2338391d0-r*(12.d0*6.77931d-3-
     &    r*20.d0*1.17871643d0)))
        d3fp12 = 
     &    r*r*r*r*(6.d0*9.354d-7-r*(24.d0*1.2338391d0-
     &    r*(60.d0*6.77931d-3-
     &    r*120.d0*1.17871643d0)))
        y = q*sqrt(q)/1.329340388d0
        dy = 1.5d0*y/q
        d2y = 0.5d0*dy/q
        d3y = -0.5d0*d2y/q
        d3fp12 = d3y*fp12 + 3.d0*d2y*dfp12 + 3.d0*dy*d2fp12 + y*d3fp12
        d2fp12 = d2y*fp12 + 2.d0*dy*dfp12 + y*d2fp12
        dfp12 = dy*fp12 + y*dfp12
        fp12 = y*fp12
      endif
      end
