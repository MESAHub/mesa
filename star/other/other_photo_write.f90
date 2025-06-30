! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module other_photo_write

   use star_def

   implicit none

   ! note: there is no flag "use_other_photo_write".
   ! the other routine is always called when a photo is written.
   ! see private/model_out.f

contains

   subroutine default_other_photo_write(id, iounit)
      integer, intent(in) :: id, iounit
      !write(iounit) stuff
   end subroutine default_other_photo_write

end module other_photo_write

