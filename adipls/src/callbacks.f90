! ***********************************************************************
!
!   Copyright (C) 2025  The MESA Team
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

module adipls_callbacks
   implicit none

   procedure(), pointer :: spcout_adi_ptr, modmod_ptr, resdif_ptr => NULL()
end module adipls_callbacks

subroutine spcout_adi(x, y, aa, data, nn, iy, iaa, ispcpr)
   use const_def, only: dp
   use adipls_callbacks, only: spcout_adi_ptr
   integer :: nn, iy, iaa, ispcpr
   real(dp) :: x(1:nn), y(1:iy,1:nn), aa(1:iaa,1:nn), data(8)

   if (associated(spcout_adi_ptr)) then
      call spcout_adi_ptr(x, y, aa, data, nn, iy, iaa, ispcpr)
   end if
end subroutine spcout_adi

subroutine modmod(x, aa, data, nn, ivarmd, iaa, imdmod)
   use const_def, only: dp
   use adipls_callbacks, only: modmod_ptr
   integer :: nn, ivarmd, iaa, imdmod
   real(dp) :: x(nn), aa(iaa,nn), data(8)

   if (associated(modmod_ptr)) then
      call modmod_ptr(x, aa, data, nn, ivarmd, iaa, imdmod)
   end if
end subroutine modmod

subroutine resdif
   ! Complex signature, unused
   return
end subroutine resdif
