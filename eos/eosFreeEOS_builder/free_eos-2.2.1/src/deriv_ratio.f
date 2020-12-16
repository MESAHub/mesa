C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: deriv_ratio.f 352 2006-04-13 02:12:47Z airwin $
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
      function deriv_ratio(x,y,n)
C       calculate nth derivative of ratio x/y of two functions x(t) and y(t).
C       all derivatives of x and y up to the nth
C       should be supplied in the arrays x(0:n), y(0:n)
      implicit none
      double precision inverse(0:100), square(0:100),
     &   deriv_ratio,x,y,deriv_product
      integer k,n
      dimension x(0:n),y(0:n)
      if(n.gt.100) stop 'n too large for dimensions in deriv_ratio'
C       form inverse derivatives.
      inverse(0) = 1.d0/y(0)
      do k = 1,n
        square(k-1) = deriv_product(inverse,inverse,k-1)
        inverse(k) = -deriv_product(square,y(1),k-1)
      enddo
      deriv_ratio = deriv_product(x,inverse,n)
      end
