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

module other_eval_i_rot

   ! consult star/other/README for general usage instructions
   ! control name: use_other_eval_i_rot = .true.
   ! procedure pointer: s% other_eval_i_rot => my_routine

   implicit none

contains

   subroutine null_other_eval_i_rot(id, k, r00, w_div_w_crit_roche, i_rot)
      use star_def
      use auto_diff_support
      integer, intent(in) :: id, k
      real(dp), intent(in) :: r00, w_div_w_crit_roche
      type(auto_diff_real_star_order1), intent(out) :: i_rot

      i_rot = 0

      write (*, *) 'no implementation for other_eval_i_rot'
      stop

   end subroutine null_other_eval_i_rot

end module other_eval_i_rot

