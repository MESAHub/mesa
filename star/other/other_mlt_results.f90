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

module other_mlt_results

   ! consult star/other/README for general usage instructions
   ! control name: use_other_mlt_results = .true.
   ! procedure pointer: s% other_mlt_results => my_routine

   implicit none

contains

   subroutine null_other_mlt_results(id, k, MLT_option, &  ! NOTE: k=0 is a valid arg
                                     r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
                                     iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
                                     alpha_semiconvection, thermohaline_coeff, &
                                     mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
      use const_def, only: dp
      use auto_diff
      use star_def
      integer, intent(in) :: id
      integer, intent(in) :: k
      character(len=*), intent(in) :: MLT_option
      type(auto_diff_real_star_order1), intent(in) :: &
         r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height
      integer, intent(in) :: iso
      real(dp), intent(in) :: &
         XH1, cgrav, m, gradL_composition_term, &
         mixing_length_alpha, alpha_semiconvection, thermohaline_coeff
      integer, intent(out) :: mixing_type
      type(auto_diff_real_star_order1), intent(out) :: &
         gradT, Y_face, conv_vel, D, Gamma
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_other_mlt_results

end module other_mlt_results

