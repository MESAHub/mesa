C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Irwin stellar interior equation of state code
C       Copyright (C) 1996 by Alan W. Irwin
C
C       $Id: abund_process.f 14 2004-09-09 23:07:54Z airwin $
C
C       For the latest version of this source code, please contact
C       Alan W. Irwin
C       Department of Physics and Astronomy
C       University of Victoria,
C       Box 3055
C       Victoria, B.C., Canada
C       V8W 3P6
C       e-mail: irwin@uvastro.phys.uvic.ca.
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
C       of the Irwin stellar interior equation of state code
C*******************************************************************************
      subroutine abund_process(eps, neps)
C       for standalone statef.
C       call abund_read to read in abundance and weight data, then
C       calculate eps(i) = (abund(i)/weight(i)),
      implicit none
      integer nzmix, i, neps
      parameter (nzmix=18)
      double precision eps(neps), xabund, yabund, zabund, zmix(nzmix)
      double precision xa(nzmix+2), xaglobal(neps)
      data xa/
     &  1.008d0, 4.0026d0, 12.0111d0, 14.0067d0, 15.9994d0, 
     &  20.179d0, 22.9898d0, 24.305d0, 26.9815d0, 28.086d0, 
     &  30.9738d0, 32.06d0, 35.453d0, 39.948d0, 40.08d0,
     &  47.9d0, 51.996d0, 54.938d0, 55.847d0, 58.71d0/
      integer ifcalled
      data ifcalled/0/
      save
      if(neps.ne.nzmix+2) stop 'abund_process: bad neps value'
      call abund_read(xabund, yabund, zabund, zmix, xa, nzmix)
C       abund_read returns properly normalized xabund, yabund, zabund, and
C       zmix.  Also may update xa.
      eps(1) = xabund/xa(1)
      eps(2) = yabund/xa(2)
      do i = 3, neps
        eps(i) = zabund*zmix(i-2)/xa(i)
      enddo
      write(*,'(a)') 'eps ='
      write(*,'(1p5d25.15)') eps
      ifcalled = 1
      return
      entry abund_weight(xaglobal, neps)
      if(neps.ne.nzmix+2) stop 'abund_weight: bad neps value'
      do i = 1, neps
        xaglobal(i) = xa(i)
      enddo
      end
