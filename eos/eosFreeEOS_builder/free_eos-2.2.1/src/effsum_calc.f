C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: effsum_calc.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine effsum_calc(f, g, morder, norder, coeff,
     &  nderiv, result)
C      Calculate EFF-style double polynomial in f and g.
C      Return result and its ln f and ln g derivatives (in order f, g,
C      ff, fg, gg, fff, ffg, fgg, ggg), in result(0:nderiv)
      implicit none
      integer morder, norder, l, n, m, nderiv
      double precision f, g, coeff((morder+1)*(norder+1)), result(0:nderiv),
     &  sumf0, sumff, sumff2, sumff3
      if(nderiv.ne.2.and.nderiv.ne.9)
     &  stop 'effsum_calc: bad nderiv value'
      result(0) = 0.d0  !double sum over f and g index.
      result(1) = 0.d0  !logarithmic derivative wrt f.
      result(2) = 0.d0  !logarithmic derivative wrt g.
      if(nderiv.eq.9) then
        result(3) = 0.d0    !ff
        result(4) = 0.d0    !fg
        result(5) = 0.d0    !gg
        result(6) = 0.d0    !fff
        result(7) = 0.d0    !ffg
        result(8) = 0.d0    !fgg
        result(9) = 0.d0    !ggg
      endif
C       do double sum in reverse order using nested multiplication.
      l = (morder+1)*(norder+1)+1
      do n = norder,0,-1
        sumf0 = 0.d0  !sum over f index.
        sumff = 0.d0  !ln f derivative of same
        if(nderiv.eq.9) then
          sumff2 = 0.d0  !second ln f derivative of same
          sumff3 = 0.d0  !third ln f derivative of same
        endif
        do m = morder,0,-1
          l = l-1
          sumf0 = coeff(l) + sumf0*f
          sumff = dble(m)*coeff(l) + sumff*f
          if(nderiv.eq.9) then
            sumff2 = dble(m*m)*coeff(l) + sumff2*f
            sumff3 = dble(m*m*m)*coeff(l) + sumff3*f
          endif
        enddo
C           double sum
        result(0) = sumf0 + result(0)*g
C           ln f derivative of same
        result(1) = sumff + result(1)*g
C           ln g derivative of same
        result(2) = dble(n)*sumf0 + result(2)*g  
        if(nderiv.eq.9) then
C             lnf lnf
          result(3) = sumff2 + result(3)*g
C             ln f ln g
          result(4) = dble(n)*sumff + result(4)*g
C             ln g ln g
          result(5) = dble(n*n)*sumf0 + result(5)*g
C             ln f ln f ln f
          result(6) = sumff3 + result(6)*g
C             ln f ln f ln g
          result(7) = dble(n)*sumff2 + result(7)*g
C             ln f ln g ln g
          result(8) = dble(n*n)*sumff + result(8)*g
C             ln g ln g ln g
          result(9) = dble(n*n*n)*sumf0 + result(9)*g      
        endif
      enddo
      end
      
