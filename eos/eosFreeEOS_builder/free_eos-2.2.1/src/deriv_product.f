C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: deriv_product.f 352 2006-04-13 02:12:47Z airwin $
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
      function deriv_product(x,y,n)
C       calculate nth derivative of product of two functions x(t) and y(t)
C       using leibnitz's rule.  all derivatives of x and y up to the nth
C       should be supplied in the arrays x(0:n), y(0:n)
      implicit none
      double precision deriv_product,x,y,binom,sum
      integer n,k
      dimension x(0:n),y(0:n)
      binom = 1.d0
      sum = 0.d0
      do k = 0,n
        sum = sum + binom*x(n-k)*y(k)
C         binomial coefficient for next loop:
C         n!/((k+1)!(n-(k+1))!) is an integer so should be exact
C           until binom*n exceeds 10**17.  by stirling's formula
C           the maximum value of binom (at k=n/2) is
C           2**(n+1)/(2pi*n)**(1/2) = 10**17 at n > 50 or so.
        binom = (binom*(n-k))/(k+1)
      enddo
      deriv_product = sum
      end
