C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: round_ln.f 352 2006-04-13 02:12:47Z airwin $
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
      double precision function round_ln(value_ln, isig_in)
C      round value to isig_in digits beyond the decimal and return as
C      round_ln.  The absolute rounding error is typically 10**(-isig_in).
C      if input variable is the natural log of some quantity, then
C      then 10**(-isig_in) is the approximate relative precision of that
C      rounded quantity.
      implicit none
      double precision value_ln
      integer isig_in, isig
      character*30 string
C       execution time format not available on sun.....sob.
      character*8 format_string
      data format_string /'(f30.00)'/
      isig = max(1,min(17,isig_in))
      write(format_string(6:7),'(i2.2)') isig
      write(string,format_string) value_ln
      read(string,format_string) round_ln
      end
