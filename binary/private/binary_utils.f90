! ***********************************************************************
!
!   Copyright (C) 2010-2019  Pablo Marchant & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
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


module binary_utils
   
   use const_def
   use math_lib
   use star_lib
   use star_def
   use binary_def
   
   implicit none

contains
   
   subroutine set_ignore_rlof_flag(binary_id, ignore_rlof_flag, ierr)
      integer, intent(in) :: binary_id
      logical, intent(in) :: ignore_rlof_flag
      integer, intent(out) :: ierr
      
      type (binary_info), pointer :: b
      
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return
      
      b% ignore_rlof_flag = ignore_rlof_flag
   end subroutine set_ignore_rlof_flag
   
   subroutine set_point_mass_i(binary_id, point_mass_i, ierr)
      integer, intent(in) :: binary_id
      integer, intent(in) :: point_mass_i
      integer, intent(out) :: ierr
      
      type (binary_info), pointer :: b
      
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return
      
      b% point_mass_i = point_mass_i
      if (point_mass_i == 1 .and. .not. b% have_star_2) then
         ierr = -1
         write(*, *) 'ERROR: setting point_mass_i=1 with have_star_2=.false.'
      else if (point_mass_i == 2 .and. .not. b% have_star_1) then
         ierr = -1
         write(*, *) 'ERROR: setting point_mass_i=2 with have_star_1=.false.'
      end if
      
      if (point_mass_i == 1 .and. b% d_i == 1) then
         b% d_i = 2
         b% a_i = 1
         b% s_donor => b% s2
         b% s_accretor => b% s1
         b% mtransfer_rate = 0d0
         b% change_factor = b% max_change_factor
      else if (point_mass_i == 2 .and. b% d_i == 2) then
         b% d_i = 1
         b% a_i = 2
         b% s_donor => b% s1
         b% s_accretor => b% s2
         b% mtransfer_rate = 0d0
         b% change_factor = b% max_change_factor
      end if
   end subroutine set_point_mass_i
   
   subroutine set_model_twins_flag(binary_id, model_twins_flag, ierr)
      integer, intent(in) :: binary_id
      logical, intent(in) :: model_twins_flag
      integer, intent(out) :: ierr
      
      type (binary_info), pointer :: b
      
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return
      
      b% model_twins_flag = model_twins_flag
      
      ! also need to set ignore_rlog_flag to true
      if (model_twins_flag) then
         call set_ignore_rlof_flag(binary_id, .true., ierr)
         if (ierr /= 0) return
         call set_point_mass_i(binary_id, 2, ierr)
      end if
   end subroutine set_model_twins_flag
   
   subroutine set_m1(binary_id, m1, ierr)
      integer, intent(in) :: binary_id
      real(dp), intent(in) :: m1
      integer, intent(out) :: ierr
      
      type (binary_info), pointer :: b
      
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return
      
      if (b% point_mass_i /= 1) then
         ierr = -1
         write(*, *) "ERROR: adjusting m1 with point_mass_i/=1 has no effect"
         return
      end if
      b% m(1) = m1 * Msun
      if (b% model_twins_flag) b% m(2) = m1
      call set_separation_eccentricity(binary_id, b% separation, b% eccentricity, ierr)
   end subroutine set_m1
   
   subroutine set_m2(binary_id, m2, ierr)
      integer, intent(in) :: binary_id
      real(dp), intent(in) :: m2
      integer, intent(out) :: ierr
      
      type (binary_info), pointer :: b
      
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return
      
      if (b% point_mass_i /= 2) then
         ierr = -1
         write(*, *) "ERROR: adjusting m2 with point_mass_i/=2 has no effect"
         return
      end if
      b% m(2) = m2 * Msun
      call set_separation_eccentricity(binary_id, b% separation, b% eccentricity, ierr)
   end subroutine set_m2
   
   subroutine set_period_eccentricity(binary_id, period, eccentricity, ierr)
      integer, intent(in) :: binary_id
      real(dp) :: period ! in seconds
      real(dp) :: eccentricity
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return
      
      b% eccentricity = eccentricity
      b% period = period
      b% separation = &
         pow((standard_cgrav * (b% m(1) + b% m(2))) * pow2(b% period / (2 * pi)), one_third)
      call set_angular_momentum_j(binary_id)
   
   end subroutine set_period_eccentricity
   
   subroutine set_separation_eccentricity(binary_id, separation, eccentricity, ierr)
      integer, intent(in) :: binary_id
      real(dp) :: separation ! in cm
      real(dp) :: eccentricity
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return
      
      b% eccentricity = eccentricity
      b% separation = separation
      b% period = &
         (2 * pi) * sqrt(b% separation * b% separation * b% separation / &
            (standard_cgrav * (b% m(1) + b% m(2))))
      call set_angular_momentum_j(binary_id)
   
   end subroutine set_separation_eccentricity
   
   subroutine set_angular_momentum_j(binary_id)
      ! Sets b% angular_momentum_j in terms of the masses, separation and eccentricity
      ! also sets the Roche lobe sizes and relative overflows
      integer, intent(in) :: binary_id
      type (binary_info), pointer :: b
      integer :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return
      
      b% angular_momentum_j = b% m(1) * b% m(2) * sqrt(standard_cgrav * &
         b% separation * (1 - pow2(b% eccentricity)) / (b% m(1) + b% m(2)))
      
      b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
      b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
      b% rl_relative_gap(1) = (b% r(1) - b% rl(1) * (1 - b% eccentricity)) / &
         b% rl(1) / (1 - b% eccentricity) ! gap < 0 means out of contact
      b% rl_relative_gap(2) = (b% r(2) - b% rl(2) * (1 - b% eccentricity)) / &
         b% rl(2) / (1 - b% eccentricity) ! gap < 0 means out of contact
      
      b% ignore_hard_limits_this_step = .true.
   
   end subroutine set_angular_momentum_j
   
   real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
      real(dp), intent(in) :: m1, m2, a
      real(dp) :: q
      q = pow(m1 / m2, one_third)
      ! Roche lobe size for star of mass m1 with a
      ! companion of mass m2 at separation a, according to
      ! the approximation of Eggleton 1983, apj 268:368-369
      rlobe = a * 0.49d0 * q * q / (0.6d0 * q * q + log1p(q))
   end function eval_rlobe

end module binary_utils

