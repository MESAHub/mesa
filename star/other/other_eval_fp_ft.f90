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

module other_eval_fp_ft

   ! consult star/other/README for general usage instructions
   ! control name: use_other_eval_fp_ft = .true.
   ! procedure pointer: s% other_eval_fp_ft => my_routine

   implicit none

contains

   subroutine null_other_eval_fp_ft( &
      id, nz, xm, r, rho, aw, ft, fp, r_polar, r_equatorial, report_ierr, ierr)
      use num_lib
      use star_utils
      use auto_diff_support
      use star_def
      integer, intent(in) :: id
      integer, intent(in) :: nz
      real(dp), intent(in) :: aw(:), r(:), rho(:), xm(:)  ! (nz)
      type(auto_diff_real_star_order1), intent(out) :: ft(:), fp(:)  ! (nz)
      real(dp), intent(inout) :: r_polar(:), r_equatorial(:)  ! (nz)
      logical, intent(in) :: report_ierr
      integer, intent(out) :: ierr

      write (*, *) 'no implementation for other_eval_fp_ft'
      ! must set fp, ft, r_polar and r_equatorial
      ierr = -1

   end subroutine null_other_eval_fp_ft

end module other_eval_fp_ft

