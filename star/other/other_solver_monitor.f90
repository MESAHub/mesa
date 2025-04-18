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

module other_solver_monitor

   ! consult star/other/README for general usage instructions
   ! control name: use_other_solver_monitor = .true.
   ! procedure pointer: s% other_solver_monitor => my_routine

   implicit none

contains

   subroutine default_other_solver_monitor( &
      id, iter, passed_tol_tests, &
      correction_norm, max_correction, &
      residual_norm, max_residual, ierr)
      use star_def
      integer, intent(in) :: id, iter
      ! iter is the number of the iteration we have just finished
      logical, intent(in) :: passed_tol_tests
      real(dp), intent(in) :: correction_norm, max_correction, &
                              residual_norm, max_residual
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine default_other_solver_monitor

end module other_solver_monitor

