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

module other_diffusion

   ! consult star/other/README for general usage instructions
   ! control name: use_other_diffusion = .true.
   ! procedure pointer: s% other_diffusion => my_routine

   implicit none

contains

   subroutine null_other_diffusion(id, dt, ierr)
      use star_def
      integer, intent(in) :: id
      real(dp), intent(in) :: dt
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      !write(*,*) 'null_other_diffusion'
      !write(*,*) 'associated(s% edv)', associated(s% edv)
      s%edv(:, 1:s%nz) = 0
      !write(*,*) 'done null_other_diffusion'
   end subroutine null_other_diffusion

end module other_diffusion

