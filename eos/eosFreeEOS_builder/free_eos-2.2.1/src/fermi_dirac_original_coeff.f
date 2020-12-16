C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: fermi_dirac_original_coeff.f 370 2006-11-29 23:57:39Z airwin $
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
      subroutine fermi_dirac_original_coeff(ccoeff,morder)
C      subroutine to calculate coefficients for thermodynamically
C      consistent set of fermi-dirac integrals as in Eggleton et al.
C      1973, A&A 23,325.
C      have choice of original 3rd-order thermodynamically
C      consistent approximations based on P_{m,n} in Table 2 of EFF
C      results are closely equivalent (within rounding errors) to
C      EFF Table 3.

      implicit none
      integer morder, k, j, i
C       n.b. ccoeff in order of rho, p, q.
      double precision ccoeff(0:morder,0:morder,3),
     &  p2coeff(0:2,0:2), round_ln

C     Coefficients for mf, mg, m =     3    3    9
      data p2coeff/
     &  2.315474d0, 4.432633d0, 2.132280d0,
     &  5.522279d0, 9.169358d0, 3.345819d0,
     &  3.693282d0, 5.166589d0, 1.334126d0/
      logical ifround
      data ifround/.true./
!      data ifround/.false./
      save
      if(morder.eq.3) then
        call fermi_dirac_recursion(p2coeff, ccoeff, morder)
        if(ifround) then
C          For almost complete verisimilitude, round identically to original
C          eff and PTEH 6 figures after decimal.  This of course introduces
C          some thermodynamic inconsistency so set ifround
C          above to .false. in the unlikely event that you ever want to
C         _use_this for a normal EOS calculation.
          do k=1,3
            do j=0,morder
              do i=0,morder
                ccoeff(i,j,k) = round_ln(ccoeff(i,j,k),6)
              enddo
            enddo
          enddo
          write(0,'(a,/,(4f15.6))') 'Rounded EFF coefficients =',
     &      ccoeff
        endif
      else
        stop 'invalid morder argument to fermi_dirac_original_coeff'
      endif
      end
