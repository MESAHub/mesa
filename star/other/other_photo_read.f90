! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   https://mesastar.org/
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

module other_photo_read

   implicit none

   ! note: there is no flag "use_other_photo_read".
   ! the other routine is always called when a photo is read.
   ! see private/model_in.f

contains

   subroutine default_other_photo_read(id, iounit, ierr)
      use star_def
      integer, intent(in) :: id, iounit
      integer, intent(out) :: ierr
      ierr = 0
      !read(iounit,iostat=ierr) stuff
   end subroutine default_other_photo_read

end module other_photo_read

