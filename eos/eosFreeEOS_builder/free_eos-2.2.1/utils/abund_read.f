C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Irwin stellar interior equation of state code
C       Copyright (C) 1996 by Alan W. Irwin
C
C       $Id: abund_read.f 372 2006-11-30 02:56:07Z airwin $
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
      subroutine abund_read(x, y, z, zmix, weight, nzmix)
C       subroutine for standalone statef.
C       read in from abund.dat, the abundances by weight, y,z, the atomic
C       weight scale (filled with defaults on call, but optionally can
C       be updated from abund.dat), and the
C       individual abundances (as fractions of z) zmix(i), i=1,18 in the
C       usual vdb and Los Alamos opacity order of 
C         6,  7,  8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24, 25, 26, 28
C         C,  N,  O, Ne, Na, Mg, Al, Si,  P,  S, Cl, Ar, Ca, Ti, Cr, Mn, Fe, Ni.
C       Calculate x = 1 - y - z.
C       (This basic scenario now changed a bit with the flag if_weight_abundance).
      implicit none
      integer nzmix
      double precision x, y, z, zmix(nzmix), weight(nzmix+2), sum
      integer idecide(28), i, number, if_weight_abundance,
     &  if_weight_default
      data idecide /1,2,0,0,0,3,4,5,0,6,7,8,9,10,11,12,13,14,0,15,0,16,
     &  0,17,18,19,0,20/
      character answer
      if(nzmix.ne.18) stop 'bad nzmix to abund_read'
      write(*,'(1x,''use pure element (y or n) '',$)')
      read(*,'(a)') answer
      if(answer.eq.'y') then
C        Pure element case
C         zero everything to start and (implicitly)_use_default 
C         atomic weight
        x=0.d0
        y=0.d0
        z=0.d0
        do i = 1,nzmix
          zmix(i) = 0.d0
        enddo
        write(*,'(1x,''atomic number? '',$)')
        read(*,*) number
C         adjust abundances for pure mixture with specified atomic number.
        if(number.lt.1.or.number.gt.28) then
          stop 'bad atomic number input'
        elseif(idecide(number).eq.0) then
          stop 'bad atomic number input'
        elseif(number.eq.1) then
          x=1.d0
        elseif(number.eq.2) then
          y=1.d0
        else
          z=1.d0
          zmix(idecide(number)-2) = 1.d0
        endif
      else
C        Mixture of elements case.
        open(unit=99,file='abund.dat',form='formatted',status='old')
        read(99,*) y, z, if_weight_abundance,
     &    if_weight_default
        read(99,*) zmix
        if(if_weight_default.ne.1) read(99,*) weight
        close (99)
        if(if_weight_abundance.gt.0) then
C          At least x, y and z are abundance by weight
          x = 1.d0 - y - z
          if(x.lt.1.d-10) then
            x = 0.d0
            y = 1.d0 - z
          endif
          if(if_weight_abundance.gt.1) then
C             this option has relative number abundances input for metals.
            do i = 1, nzmix
              zmix(i) = zmix(i)*weight(i+2)
            enddo
          endif
          sum = 0.d0
          do i = 1, nzmix
            sum = sum + zmix(i)
          enddo
          if(sum.gt.0.d0) then
            do i = 1, nzmix
              zmix(i) = zmix(i)/sum
            enddo
          elseif(z.eq.0.d0) then
C             okay to have zmix = 0 in this case.
          else
            stop 'abund_read: zero zmix and non-zero z is an error'
          endif
        else
C          for this zero or negative if_weight_abundance,
C          zmix is used to input relative
C          number abundance of the metals and y and z are misnomers
C          since they input the relative number abundance of H and He.
          if(if_weight_abundance.eq.-2) then
C            convert log input number abundance that must be renormalized
C            to number abundance that must be renormalized.
            y = 10.d0**y
            z = 10.d0**z
            do i = 1,nzmix
              zmix(i) = 10.d0**zmix(i)
            enddo
          endif
          sum = 0.d0
          do i = 1, nzmix
            zmix(i) = zmix(i)*weight(i+2)
            sum = sum + zmix(i)
          enddo
C           for this option, input y and z are misnomers.
C           they are actually used to store the
C           relative number abundances of hydrogen and helium.
          x = y*weight(1)
          y = z*weight(2)
          z = sum
C           normalize zmix to unity unless z = 0
          if(z.gt.0.d0) then
          do i = 1, nzmix
            zmix(i) = zmix(i)/z
          enddo
          endif
          if(if_weight_abundance.eq.-1) then
C             for this special case, input is eps array so don't renormalize.
          else
            sum = x + y + z
            x = x/sum
            y = y/sum
            z = z/sum
          endif
        endif
      endif
      write(*,'(a,f20.16)') ' x = ', x
      write(*,'(a,f20.16)') ' y = ', y
      write(*,'(a,f20.16)') ' z = ', z
      write(*,'(a,/,(5f20.16))') ' zmix = ', zmix
      write(*,'(a,/,(5f20.16))') ' atomic weights = ', weight
      end
