! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton, Pablo Marchant & The MESA Team
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


      module binary_timestep

      use const_def
      use math_lib
      use star_lib
      use star_def
      use binary_def

      implicit none

      contains

      subroutine set_star_timesteps(b) ! sets the smallest next timestep for all stars
         type (binary_info), pointer :: b
         integer :: i, l
         real(dp) :: dt_min, rel_overlap
         type (star_info), pointer :: s
         integer :: ierr, num_stars
         ierr = 0
         dt_min = 1d99
         num_stars = 2
         if (b% point_mass_i > 0) num_stars = 1
         do l = 1, num_stars
            if (l == 1 .and. b% point_mass_i == 1) then
               i = 2
            else
               i = l
            end if
            call star_ptr(b% star_ids(i), s, ierr)
            if (ierr /= 0) then
                write(*, *) trim('star_ptr') // ' ierr', ierr
                return
            end if
            if (s% dt_next < dt_min) then
               dt_min = s% dt_next
            end if
         end do
         if (b% max_timestep <= dt_min) then
            dt_min = b% max_timestep
         else
            b% dt_why_reason = b_Tlim_comp
         end if
         ! just to be sure we dont cause a segfault
         if (b% dt_why_reason < 1 .or. b% dt_why_reason > b_numTlim) then
            dt_why_str(Tlim_binary) = " "
         else
            dt_why_str(Tlim_binary) = binary_dt_why_str(b% dt_why_reason)
         end if
         do l = 1, num_stars
            if (l == 1 .and. b% point_mass_i == 1) then
               i = 2
            else
               i = l
            end if
            call star_ptr(b% star_ids(i), s, ierr)
            if (ierr /= 0) then
                write(*, *) trim('star_ptr') // ' ierr', ierr
                return
            end if
            if (s% dt_next > dt_min) then
               s% dt_next = dt_min
               s% why_Tlim = Tlim_binary
            end if
         end do

         if (b% have_to_reduce_timestep_due_to_j) then
            ! lower timesteps after retries due to large changes in angular momentum
            if (b% point_mass_i /= 1) then
               b% s1% dt = b% s1% dt*b% dt_reduction_factor_for_j
            end if
            if (b% point_mass_i /= 2) then
               b% s2% dt = b% s2% dt*b% dt_reduction_factor_for_j
            end if
            b% have_to_reduce_timestep_due_to_j = .false.
         end if
         
      end subroutine set_star_timesteps

      integer function binary_pick_next_timestep(b)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         
         real(dp) :: &
            env_change, dtm, dtj, dta, dtr, dte, dtdm, &
            j_change, sep_change, rel_gap_change, e_change, set_dt, &
            rel_change

         include 'formats'

         dtm = 1d99
         dtj = 1d99
         dta = 1d99
         dtr = 1d99
         dte = 1d99
         dtdm = 1d99

         binary_pick_next_timestep = keep_going

         s => b% s_donor

         if (b% max_timestep < 0) b% max_timestep = b% s_donor% dt

         b% env = s% mstar/msun - s% he_core_mass 
         if (b% env_old /= 0) then
            env_change = b% env - b% env_old
         else
            env_change = 0
         end if
         
         if (b% rl_relative_gap_old(b% d_i) /= 0) then
            rel_gap_change = b% rl_relative_gap_old(b% d_i) - b% rl_relative_gap(b% d_i)
         else
            rel_gap_change = 0
         end if
         
         if (b% angular_momentum_j_old /= 0) then
            j_change = b% angular_momentum_j - b% angular_momentum_j_old
         else
            j_change = 0
         end if
         
         if (b% separation_old /= 0) then
            sep_change = b% separation - b% separation_old
         else
            sep_change = 0
         end if
         if (b% eccentricity_old /= 0) then
             e_change = b% eccentricity - b% eccentricity_old
         else
             e_change = 0
         end if
   
         ! get limits for dt based on relative changes
         if (b% fj > 0) then
            rel_change = abs(j_change/b% angular_momentum_j)
            if (.not. b% ignore_hard_limits_this_step .and. &
               b% fj_hard > 0d0 .and. rel_change > b% fj_hard) then
               write(*,*) "retry because of fj_hard limit,", &
                  "fj_hard:", b% fj_hard, "rel_change:", rel_change
               binary_pick_next_timestep = retry
               b% have_to_reduce_timestep_due_to_j = .true.
               return
            end if
            dtj = s% dt/secyer/(rel_change/b% fj+1d-99)
         end if

         if (b% fm > 0) then
            rel_change = abs(env_change/max(b% env, b% fm_limit))
            if (.not. b% ignore_hard_limits_this_step .and. &
               b% fm_hard > 0d0 .and. rel_change > b% fm_hard) then
               write(*,*) "retry because of fm_hard limit,", &
                  "fm_hard:", b% fm_hard, "rel_change:", rel_change
               binary_pick_next_timestep = retry
               return
            end if
            dtm = s% dt/secyer/(rel_change/b% fm+1d-99)
         end if
         
         if (b% fr > 0) then
            rel_change = abs(rel_gap_change/max(-b% rl_relative_gap(b% d_i), b% fr_limit))
            if (.not. b% ignore_hard_limits_this_step .and. &
               b% fr_hard > 0d0 .and. rel_change > b% fr_hard) then
               write(*,*) "retry because of fr_hard limit for donor,", &
                  "fr_hard:", b% fr_hard, "rel_change:", rel_change
               binary_pick_next_timestep = retry
               return
            end if
            dtr = s% dt/secyer/ &
                (rel_change/b% fr+1d-99)

            ! Check for accretor as well
            if (b% rl_relative_gap_old(b% a_i) /= 0) then
               rel_gap_change = b% rl_relative_gap_old(b% a_i) - b% rl_relative_gap(b% a_i)
            else
               rel_gap_change = 0
            end if
            rel_change = abs(rel_gap_change/max(-b% rl_relative_gap(b% a_i), b% fr_limit))
            if (.not. b% ignore_hard_limits_this_step .and. &
               b% fr_hard > 0d0 .and. rel_change > b% fr_hard) then
               write(*,*) "retry because of fr_hard limit for accretor,", &
                  "fr_hard:", b% fr_hard, "rel_change:", rel_change
               binary_pick_next_timestep = retry
               return
            end if
            dtr = min(dtr, s% dt/secyer/ &
                (rel_change/b% fr+1d-99))
         end if
         if (dtr < b% fr_dt_limit) dtr = b% fr_dt_limit

         if (b% fa > 0) then
            rel_change = abs(sep_change/b% separation)
            if (.not. b% ignore_hard_limits_this_step .and. &
               b% fa_hard > 0d0 .and. rel_change > b% fa_hard) then
               write(*,*) "retry because of fa_hard limit,", &
                  "fa_hard:", b% fa_hard, "rel_change:", rel_change
               binary_pick_next_timestep = retry
               return
            end if
            dta = s% dt/secyer/(rel_change/b% fa+1d-99)
         end if

         if (b% fe > 0) then
            rel_change = abs(e_change/ max( b% eccentricity, b% fe_limit ))
            if (.not. b% ignore_hard_limits_this_step .and. &
               b% fe_hard > 0d0 .and. rel_change > b% fe_hard) then
               write(*,*) "retry because of fe_hard limit,", &
                  "fe_hard:", b% fe_hard, "rel_change:", rel_change
               binary_pick_next_timestep = retry
               return
            end if
            dte = s% dt/secyer/(rel_change/b% fe+1d-99)
         end if

         if (b% fdm > 0d0) then
            rel_change = abs(b% m(b% d_i) - b% m_old(b% d_i))/b% m_old(b% d_i)
            if (.not. b% ignore_hard_limits_this_step .and. &
               b% fdm_hard > 0d0 .and. rel_change > b% fdm_hard) then
               write(*,*) "retry because of fdm_hard limit for donor,", &
                  "fdm_hard:", b% fdm_hard, "rel_change:", rel_change
               binary_pick_next_timestep = retry
               return
            end if
            dtdm = s% dt/secyer/(rel_change/b% fdm+1d-99)

            rel_change = abs(b% m(b% a_i) - b% m_old(b% a_i))/b% m_old(b% a_i)
            if (.not. b% ignore_hard_limits_this_step .and. &
               b% fdm_hard > 0d0 .and. rel_change > b% fdm_hard) then
               write(*,*) "retry because of fdm_hard limit for accretor,", &
                  "fdm_hard:", b% fdm_hard, "rel_change:", rel_change
               binary_pick_next_timestep = retry
               return
            end if
            dtdm = min(dtdm, s% dt/secyer/(rel_change/b% fdm+1d-99))
         end if

         set_dt = min(dtm, dtr, dtj, dta, dte, dtdm)
         
         if (set_dt == dtm) then
            b% dt_why_reason = b_Tlim_env
         else if (set_dt == dtr) then
            b% dt_why_reason = b_Tlim_roche
         else if (set_dt == dtj) then
            b% dt_why_reason = b_Tlim_jorb
         else if (set_dt == dta) then
            b% dt_why_reason = b_Tlim_sep
         else if (set_dt == dte) then
            b% dt_why_reason = b_Tlim_ecc
         else if (set_dt == dtdm) then
            b% dt_why_reason = b_Tlim_dm
         else
            stop 'Something wrong in binary timestep'
         end if

         if (set_dt < 1d-7) set_dt = 1d-7 ! there's a limit to everything

         b% max_timestep = exp10(b% dt_softening_factor*log10(set_dt*secyer) + &
             (1-b% dt_softening_factor)*log10(b% max_timestep))

         ! use variable varcontrols for different phases of evolution
         if (abs(b% mtransfer_rate)/Msun*secyer > 1d-20) then
            if (b% s_donor% center_h1 > 1d-12 .and. b% varcontrol_case_a > 0d0) then
               b% s_donor% varcontrol_target = b% varcontrol_case_a
               if (b% point_mass_i == 0) &
                   b% s_accretor% varcontrol_target = b% varcontrol_case_a
            else if (b% s_donor% center_h1 < 1d-12 .and. b% varcontrol_case_b > 0d0) then
               b% s_donor% varcontrol_target = b% varcontrol_case_b
               if (b% point_mass_i == 0) &
                   b% s_accretor% varcontrol_target = b% varcontrol_case_b
            end if
         else
            if (b% s_donor% center_h1 > 1d-12) then
               if (b% varcontrol_ms > 0d0) &
                   b% s_donor% varcontrol_target = b% varcontrol_ms
            else
               if (b% varcontrol_post_ms > 0d0) &
                   b% s_donor% varcontrol_target = b% varcontrol_post_ms
            end if

            if (b% point_mass_i == 0) then
               if (b% s_accretor% center_h1 > 1d-12) then
                  if (b% varcontrol_ms > 0d0) &
                      b% s_accretor% varcontrol_target = b% varcontrol_ms
               else
                  if (b% varcontrol_post_ms > 0d0) &
                      b% s_accretor% varcontrol_target = b% varcontrol_post_ms
               end if
            end if
         end if

         b% ignore_hard_limits_this_step = .false.
         
      end function binary_pick_next_timestep
      

      end module binary_timestep
