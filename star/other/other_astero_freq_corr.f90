! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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

module other_astero_freq_corr

   ! consult star/other/README for general usage instructions
   ! control name: use_other_astero_freq_corr = .true.
   ! procedure pointer: s% other_astero_freq_corr => my_routine

   implicit none

contains

   subroutine default_other_astero_freq_corr(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      ! see star/astero/src/astero_support.f, subroutine get_freq_corr.
      ! note that your routine can use both astero_data and astero_support.
   end subroutine default_other_astero_freq_corr

end module other_astero_freq_corr

