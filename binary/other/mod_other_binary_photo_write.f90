! ***********************************************************************
!
!   Copyright (C) 2021 The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module mod_other_binary_photo_write

      use binary_def

      implicit none

      ! note: there is no flag "use_other_binary_photo_write".
      ! the other routine is always called when a photo is written.
      ! see private/binary_photos.f90

      contains


      subroutine default_other_binary_photo_write(binary_id, iounit)
         integer, intent(in) :: binary_id, iounit
         !write(iounit) stuff
      end subroutine default_other_binary_photo_write


      end module mod_other_binary_photo_write
