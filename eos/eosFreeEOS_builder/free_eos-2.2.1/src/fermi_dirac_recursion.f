C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fermi_dirac_recursion.f 352 2006-04-13 02:12:47Z airwin $
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
      subroutine fermi_dirac_recursion(pcoeff, ccoeff, morder)
C      calculate coefficients for thermodynamically
C      consistent set of fermi-dirac integrals following recursion
C      relationships in paper (see also recursion.mem).
      implicit none
      integer morder, morder_local
      parameter(morder_local=8)
C       n.b. ccoeff in order of rho, p, q.
      double precision ccoeff(0:morder,0:morder,3),
     &  plowcoeff(-2:morder_local,-1:morder_local),
     &  pcoeff(0:morder-1,0:morder-1)
      integer n, m
      save
      if(morder.gt.morder_local)
     &  stop 'fermi_dirac_recursion: invalid morder'
C       zero plowcoeff.
      do n=-1,morder
        do m=-2,morder
          plowcoeff(m,n) = 0.d0
        enddo
      enddo
C       fill in non-zero part of plowcoeff.
      do n = 0,morder-1
        do m = 0,morder-1
          plowcoeff(m,n) = pcoeff(m,n)
        enddo
      enddo
C      Calculate rhostar, pstar, and qstar coefficients from recursion
C      relations derived from thermodynamic consistency arguments.
      do n = 0,morder
        do m = 0,morder
C           rho coefficients.
          ccoeff(m,n,1) = 
     &      (1.d0+dble(m))*plowcoeff(m,n) + 
     &      (1.25d0+dble(m-morder)+0.5d0*dble(n))*plowcoeff(m-1,n) +
     &      (1.d0+dble(m))*plowcoeff(m,n-1) +
     &      (0.5d0*dble(4+n-morder)+dble(m-morder))*plowcoeff(m-1,n-1)
C           p coefficients.
          ccoeff(m,n,2) = 
     &      plowcoeff(m,n) + plowcoeff(m-1,n) +
     &      plowcoeff(m,n-1) + plowcoeff(m-1,n-1)
C           q coefficients.
          ccoeff(m,n,3) = 
     &      (2.5d0+dble(-2*(m+1)+n))*plowcoeff(m,n) + 
     &      (2.5d0+dble(-2*(2*m-morder)+n))*plowcoeff(m-1,n) +
     &      (dble(2-2*m-morder+n))*plowcoeff(m,n-1) + 
     &      (dble(4-2*(2*m-morder)-morder+n))*plowcoeff(m-1,n-1) +
     &      (dble(-2*(m-1-morder)))*plowcoeff(m-2,n) +
     &      (dble(-2*(m-1-morder)))*plowcoeff(m-2,n-1)
        enddo
      enddo
      end
