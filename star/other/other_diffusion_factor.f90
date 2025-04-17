! ***********************************************************************
!
!   Copyright (C)  2019 Evan Bauer & The MESA Team
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

module other_diffusion_factor

   ! consult star/other/README for general usage instructions
   ! control name: use_other_diffusion_factor = .true.
   ! procedure pointer: s% other_diffusion_factor => my_routine

   implicit none

contains

   subroutine default_other_diffusion_factor(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      integer :: num_classes
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (s%diffusion_use_full_net) then
         num_classes = s%species
      else
         num_classes = s%diffusion_num_classes
      end if

      s%extra_diffusion_factor(1:num_classes, 1:s%nz) = 1d0
   end subroutine default_other_diffusion_factor

end module other_diffusion_factor

