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


      module binary_evolve

      use const_def
      use math_lib
      use star_lib
      use star_def
      use binary_def
      use binary_utils, only: eval_rlobe, set_angular_momentum_j, &
         set_separation_eccentricity, set_period_eccentricity

      implicit none

      contains

         subroutine binarydata_init(b, doing_restart)
         use utils_lib, only: is_bad
         type (binary_info), pointer :: b
         logical, intent(in) :: doing_restart
         integer :: finish_step_result
         integer :: i, ierr
         real(dp) :: r_isco, Z1, Z2
         include 'formats'

         b% evolve_both_stars = b% job% evolve_both_stars
         b% warn_binary_extra = b% job% warn_binary_extra

         ! Initialize arrays for phase dependent calculations
         allocate(b% theta_co(b% anomaly_steps), b% time_co(b% anomaly_steps), &
            b% mdot_donor_theta(b% anomaly_steps))
         allocate(b% edot_theta(b% anomaly_steps), b% e1(b% anomaly_steps), &
            b% e2(b% anomaly_steps), b% e3(b% anomaly_steps))

         if (.not. doing_restart) then
            if (.not. b% evolve_both_stars) then
               b% point_mass_i = 2
            else if ((b% job% change_initial_model_twins_flag .or. b% job% change_model_twins_flag) &
                  .and. b% job% new_model_twins_flag) then
               b% point_mass_i = 2
            else
               b% point_mass_i = 0
            end if

            b% doing_first_model_of_run = .true.
            b% open_new_history_file = .true.

            b% s_donor => b% s1
            b% s_accretor => b% s2
            b% d_i = 1
            b% a_i = 2

            b% max_timestep = 1d99
            b% change_factor = b% initial_change_factor

            if (b% point_mass_i /= 1) then
               initial_mass(1) = b% s1% mstar / Msun
            else
               initial_mass(1) = b% m1
            end if
            if (b% point_mass_i /= 2) then
               initial_mass(2) = b% s2% mstar / Msun
            else
               if (.not. b% model_twins_flag) then
                  initial_mass(2) = b% m2
               else
                  initial_mass(2) = initial_mass(1)
               end if
            end if

            b% m(1) = initial_mass(1)*Msun
            b% m(2) = initial_mass(2)*Msun
            if (b% point_mass_i /= 1) then
               b% r(1) = Rsun*b% s1% photosphere_r
            else
               b% r(1) = 0
            end if
            if (b% point_mass_i /= 2) then
               b% r(2) = Rsun*b% s2% photosphere_r
            else
               if (.not. b% model_twins_flag) then
                  b% r(2) = 0
               else
                  b% r(2) = b% r(1)
               end if
            end if
            
            if (b% initial_period_in_days <= 0) then ! calculate from initial_separation_in_Rsuns
               call set_separation_eccentricity(b% binary_id, &
                  b% initial_separation_in_Rsuns*Rsun, b% initial_eccentricity, ierr)
                  if (ierr /= 0) then
                     return
                  end if
            else
               call set_period_eccentricity(b% binary_id, &
                  b% initial_period_in_days*secday, b% initial_eccentricity, ierr)
                  if (ierr /= 0) then
                     return
                  end if
            end if

         ! Set all parameters nessessary for integration over the binary orbit
         ! 1) true anomaly = polar angle from periastron 0 -> 2pi
            do i = 1,b% anomaly_steps 
               b% theta_co(i) = (i-1) * (2 * pi) / b% anomaly_steps
            end do
            ! 2) time between periastron and polar angle theta 0 -> 1 (fraction of the
            !    orbital period)
            do i = 1,b% anomaly_steps ! time between periastron and polar angle theta
               b% time_co(i) = ( 2 * atan( sqrt( (1-b% eccentricity)/(1 + b% eccentricity) ) * &
                               tan(b% theta_co(i)/2d0) ) - b% eccentricity * &
                               sqrt(1 - pow2(b% eccentricity)) * sin(b% theta_co(i)) / &
                               (1 + b% eccentricity * cos(b% theta_co(i)) ) ) /2.d0 /pi
               if (i > b% anomaly_steps/2+1) then
                  b% time_co(i) = b% time_co(i) + b% time_co(b% anomaly_steps/2+1) * 2
               end if
            end do
   
            if (is_bad(b% rl_relative_gap(1))) call mesa_error(__FILE__,__LINE__,'binarydata_init')
            if (is_bad(b% rl_relative_gap(2))) call mesa_error(__FILE__,__LINE__,'binarydata_init')
            b% using_jdot_mb(1) = .false.
            b% using_jdot_mb(2) = .false.

            if (b% point_mass_i /= 0) then
               ! this part is only relevant for BH accretors
               if (b% initial_bh_spin < 0d0) then
                  b% initial_bh_spin = 0d0
                  write(*,*) "initial_bh_spin is smaller than zero. It has been set to zero."
               else if (b% initial_bh_spin > 1d0) then
                  b% initial_bh_spin = 1d0
                  write(*,*) "initial_bh_spin is larger than one. It has been set to one."
               end if
               ! compute isco radius from eq. 2.21 of Bardeen et al. (1972), ApJ, 178, 347
               Z1 = 1d0 + pow(1d0 - pow2(b% initial_bh_spin),one_third) &
                  * (pow(1d0 + b% initial_bh_spin,one_third) + pow(1d0 - b% initial_bh_spin,one_third))
               Z2 = sqrt(3d0*pow2(b% initial_bh_spin) + pow2(Z1))
               r_isco = 3d0 + Z2 - sqrt((3d0 - Z1)*(3d0 + Z1 + 2d0*Z2))
               ! compute equivalent mass at zero spin from eq. (3+1/2) (ie. the equation between (3) and (4))
               ! of Bardeen (1970), Nature, 226, 65, taking values with subscript zero to correspond to
               ! zero spin (r_isco = 6).
               b% eq_initial_bh_mass = b% m(b% point_mass_i) * sqrt(r_isco/6d0)
            end if
            
            write(*,'(A)')
            write(*,1) 'm2', b% m2
            write(*,1) 'm1', b% m1
            write(*,1) 'initial_period_in_days', b% initial_period_in_days
            write(*,1) 'initial_separation_in_Rsun', b% separation/Rsun
            write(*,1) 'jdot_multiplier', b% jdot_multiplier
            write(*,1) 'fr', b% fr
            write(*,'(A)')

            min_binary_period = b% period
            b% min_binary_separation = b% separation
            initial_binary_period = b% period

            b% ixtra(:) = 0
            b% xtra(:) = 0d0
            b% lxtra(:) = .false.

            b% CE_num1 = 0
            b% CE_num2 = 0
            b% CE_lambda1 = 0d0
            b% CE_lambda2 = 0d0
            b% CE_Ebind1 = 0d0
            b% CE_Ebind2 = 0d0

            b% num_tries = 0

            finish_step_result = binary_finish_step(b)
         else
            if (b% d_i == b% point_mass_i) then
               write(*,*) "WARNING: restart has donor star set as point mass"
               write(*,*) "switching donor"
               b% d_i = b% a_i
               b% a_i = b% point_mass_i
            end if
            if (b% d_i == 1) then
               b% s_donor => b% s1
               b% s_accretor => b% s2
            else
               b% s_donor => b% s2
               b% s_accretor => b% s1
            end if
         end if
          
      end subroutine

      subroutine set_donor_star(b)
         type (binary_info), pointer :: b
         logical :: switch_donor
         real(dp) :: mdot_hi_temp
         include 'formats'

         switch_donor = .false.

         if (b% keep_donor_fixed .and. b% mdot_scheme /= "contact") return

         if (b% mdot_scheme == "roche_lobe" .and. &
            abs(b% mtransfer_rate/(Msun/secyer)) < b% mdot_limit_donor_switch .and. &
            b% rl_relative_gap_old(b% a_i) > b% rl_relative_gap_old(b% d_i)) then
            switch_donor = .true.
         else if (b% mtransfer_rate > 0d0) then
            switch_donor = .true.
            b% mtransfer_rate = - b% mtransfer_rate
            mdot_hi_temp = b% mdot_hi
            b% mdot_hi = - b% mdot_lo
            b% mdot_lo = - mdot_hi_temp
            if (.not. b% have_mdot_lo) then
               b% have_mdot_hi = .false.  
            end if
            b% have_mdot_lo = .true.
            b% fixed_delta_mdot = b% fixed_delta_mdot / 2.0d0
         else if (b% mdot_scheme == "contact" .and. &
            b% rl_relative_gap_old(b% a_i) > b% rl_relative_gap_old(b% d_i) .and. &
            b% rl_relative_gap_old(b% a_i) < - b% implicit_scheme_tolerance .and. &
            b% rl_relative_gap_old(b% d_i) < - b% implicit_scheme_tolerance .and. &
            abs(b% mtransfer_rate/(Msun/secyer)) < b% mdot_limit_donor_switch) then
            switch_donor = .true.
         end if

         if (switch_donor) then
            if (b% report_rlo_solver_progress) write(*,*) "switching donor"
            if (b% d_i == 2) then
               b% d_i = 1
               b% a_i = 2
               b% s_donor => b% s1
               b% s_accretor => b% s2
            else
               b% d_i = 2
               b% a_i = 1
               b% s_donor => b% s2
               b% s_accretor => b% s1
            end if
         end if
      end subroutine

      integer function binary_evolve_step(b)
         use utils_lib, only: is_bad
         use binary_jdot, only: get_jdot
         use binary_edot, only: get_edot
         type(binary_info), pointer :: b
         integer :: i
         
         include 'formats'

         ! store the final mdots used for each star
         ! for a point mass accretor this is already set in the subroutine adjust_mdots of binary_mdot
         b% component_mdot(b% d_i) = b% s_donor% mstar_dot
         if (b% point_mass_i == 0) then
            b% component_mdot(b% a_i) = b% s_accretor% mstar_dot
         else if (b% model_twins_flag) then
            b% component_mdot(b% a_i) = b% component_mdot(b% d_i)
         end if

         b% m(b% d_i) = b% s_donor% mstar
         b% time_step = b% s_donor% time_step
         if (b% point_mass_i == 0) then
            b% m(b% a_i) = b% s_accretor% mstar
         else if (.not. b% model_twins_flag) then
            b% m(b% a_i) = b% m(b% a_i) + b% component_mdot(b% a_i)*b% s_donor% dt
         else
            b% m(b% a_i) = b% m(b% d_i)
         end if
         
         if (b% point_mass_i /= 1) then
            b% r(1) = Rsun*b% s1% photosphere_r
         else
            b% r(1) = 0
         end if
         if (b% point_mass_i /= 2) then
            b% r(2) = Rsun*b% s2% photosphere_r
         else if (.not. b% model_twins_flag) then
            b% r(2) = 0
         else
            b% r(2) = b% r(1)
         end if

         ! solve the winds in the system for jdot calculation,
         ! these don't include mass lost due to mass_transfer_efficiency < 1.0
         ! Since s% mstar_dot is not just mass loss, but includes the contribution from
         ! RLO mass transfer and wind mass transfer from the other component, these
         ! need to be removed to get the actual wind. Also, the fraction of the wind
         ! that is fransferred to the other component does not leave the system, and
         ! needs to be removed as well.
         b% mdot_system_wind(b% d_i) = b% s_donor% mstar_dot - b% mtransfer_rate &
            + b% mdot_wind_transfer(b% a_i) - b% mdot_wind_transfer(b% d_i)
         if (b% point_mass_i == 0) then
            b% mdot_system_wind(b% a_i) = b% s_accretor% mstar_dot &
                + b% mtransfer_rate * b% fixed_xfer_fraction + b% mdot_wind_transfer(b% d_i) &
                - b% mdot_wind_transfer(b% a_i)
         else
            if (.not. b% model_twins_flag) then
               b% mdot_system_wind(b% a_i) = 0.0d0
            else
               b% mdot_system_wind(b% a_i) = b% mdot_system_wind(b% d_i)
            end if
         end if

         ! get jdot and update orbital J
         b% jdot = get_jdot(b)
         b% angular_momentum_j = b% angular_momentum_j + b% jdot*b% time_step*secyer

         if (b% angular_momentum_j <= 0) then
            write(*,*) 'Retry due to negative angular_momentum_j', b% angular_momentum_j
            b% have_to_reduce_timestep_due_to_j = .true.
            binary_evolve_step = retry
            return
         end if
         
         ! update the eccentricity (ignore in first step)
         if (.not. b% doing_first_model_of_run) then
            b% eccentricity = b% eccentricity + get_edot(b) *b% time_step*secyer
            if (b% eccentricity < b% min_eccentricity) b% eccentricity = b% min_eccentricity
            if (b% eccentricity > b% max_eccentricity) b% eccentricity = b% max_eccentricity
         end if
         
         !use new eccentricity to calculate new time coordinate
         do i = 1,b% anomaly_steps ! time between periastron and polar angle theta
            b% time_co(i) = ( 2 * atan( sqrt( (1-b% eccentricity)/(1 + b% eccentricity) ) * &
                            tan(b% theta_co(i)/2d0) ) - b% eccentricity * &
                            sqrt(1 - pow2(b% eccentricity)) * sin(b% theta_co(i)) / &
                            (1 + b% eccentricity * cos(b% theta_co(i)) ) ) /2.d0 /pi
            if (i > b% anomaly_steps/2+1) then
               b% time_co(i) = b% time_co(i) + b% time_co(b% anomaly_steps/2+1) * 2
            end if
         end do
         
         ! use the new j to calculate new separation
         b% separation = (pow2(b% angular_momentum_j/(b% m(1)*b% m(2)))) *&
             (b% m(1)+b% m(2)) / standard_cgrav * 1 / (1 - pow2(b% eccentricity))
         if (b% separation < b% min_binary_separation) &
            b% min_binary_separation = b% separation
         
         b% period = 2*pi*sqrt(pow3(b% separation)/&
               (standard_cgrav*(b% m(1)+b% m(2)))) 
         if (b% period < min_binary_period) min_binary_period = b% period
         
         ! use the new separation to calculate the new roche lobe radius
         
         b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
         b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
         b% rl_relative_gap(1) = (b% r(1) - b% rl(1) * (1 - b% eccentricity) ) / &
             b% rl(1) / (1 - b% eccentricity) ! gap < 0 means out of contact 
         b% rl_relative_gap(2) = (b% r(2) - b% rl(2) * (1 - b% eccentricity) ) / &
             b% rl(2) / (1 - b% eccentricity) ! gap < 0 means out of contact

         if (is_bad(b% rl_relative_gap(1)) .or. is_bad(b% rl_relative_gap(2))) then
            write(*,*) "rl_relative_gap for each component", b% rl_relative_gap(1), b% rl_relative_gap(2)
            write(*,*) "rl for each component", b% rl(1), b% rl(2)
            write(*,*) "r for each component", b% r(1), b% r(2)
            write(*,*) "m for each component", b% m(1), b% m(2)
            write(*,*) "separation, angular momentum", b% separation, b% angular_momentum_j
            write(*,*) "jdot, jdot_mb, jdot_gr, jdot_ml:", b% jdot, b% jdot_mb, b% jdot_gr, b% jdot_ml
            write(*,*) "jdot_ls, jdot_missing_wind, extra_jdot:", b% jdot_ls, b% jdot_missing_wind, b% extra_jdot
            write(*,*) 'error solving rl_rel_gap'
            binary_evolve_step = retry
            return
         end if

         b% model_number = b% model_number + 1
         b% binary_age = b% binary_age + b% time_step

         binary_evolve_step = keep_going

      end function binary_evolve_step

      integer function binary_check_model(b)
         use binary_mdot, only: rlo_mdot, check_implicit_rlo
         use binary_irradiation
         type (binary_info), pointer :: b

         integer :: i, j, ierr, id
         logical :: implicit_rlo
         real(dp) :: new_mdot, q


         include 'formats'

         binary_check_model = retry
         ierr = 0
         
         implicit_rlo = (b% max_tries_to_achieve > 0 .and. b% implicit_scheme_tolerance > 0d0)
         
         binary_check_model = keep_going
                  
         if (.not. b% ignore_rlof_flag) then
            if (implicit_rlo) then ! check agreement between new r and new rl
               if (.not. b% use_other_check_implicit_rlo) then
                  binary_check_model = check_implicit_rlo(b% binary_id, new_mdot)
               else
                  binary_check_model = b% other_check_implicit_rlo(b% binary_id, new_mdot)
               end if
               if (binary_check_model == keep_going) then
                  b% donor_started_implicit_wind = .false.
               end if
               b% donor_started_implicit_wind = b% donor_started_implicit_wind .or. &
                  b% s_donor% was_in_implicit_wind_limit
            else
               if (.not. b% use_other_rlo_mdot) then
                  call rlo_mdot(b% binary_id, new_mdot, ierr) ! grams per second
                  if (ierr /= 0) then
                     write(*,*) 'failed in rlo_mdot'
                     binary_check_model = retry
                     return
                  end if
               else
                  call b% other_rlo_mdot(b% binary_id, new_mdot, ierr)
                  if (ierr /= 0) then
                     write(*,*) 'failed in other rlo_mdot'
                     binary_check_model = retry
                     return
                  end if
               end if
               if (new_mdot > 0) then
                  new_mdot = 0.0d0
                  write(*,*) "WARNING: explicit computation of mass transfer results in accreting donor"
                  write(*,*) "Not transfering mass"
               end if
               ! smooth out the changes in mdot
               new_mdot = b% cur_mdot_frac*b% mtransfer_rate + (1-b% cur_mdot_frac)*new_mdot
               if (-new_mdot/(Msun/secyer) > b% max_explicit_abs_mdot) new_mdot = -b% max_explicit_abs_mdot*Msun/secyer 
            end if
            b% mtransfer_rate = new_mdot
         else
            b% mtransfer_rate = 0
         end if
         call adjust_irradiation(b)

         if (.not. b% CE_flag) then
            if ((b% point_mass_i == 0 .or. b% model_twins_flag) &
                  .and. b% rl_relative_gap(b% a_i) >= 0.0d0) then
               if (b% rl_relative_gap(b% a_i) >= b% accretor_overflow_terminate) then
                  binary_check_model = terminate
                  b% s_donor% termination_code = t_xtra1
                  termination_code_str(t_xtra1) = &
                      "Terminate because accretor (r-rl)/rl > accretor_overflow_terminate"
               end if
            end if
            if (b% doing_first_model_of_run .and. b% terminate_if_initial_overflow &
                  .and. (.not. b% ignore_rlof_flag .or. b% model_twins_flag)) then
               if (b% rl_relative_gap(b% d_i) >= 0.0d0 &
                     .or. (b% point_mass_i == 0 .and. b% rl_relative_gap(b% a_i) >= 0.0d0)) then
                  binary_check_model = terminate
                  b% s_donor% termination_code = t_xtra1
                  termination_code_str(t_xtra1) = &
                      "Terminate because of overflowing initial model"
               end if
            end if
            if ((b% point_mass_i == 0 .or. b% model_twins_flag) &
                  .and. b% terminate_if_L2_overflow) then
               if (b% m(1) > b% m(2)) then
                  q = b% m(2) / b% m(1)
                  id = 2
               else
                  q = b% m(1) / b% m(2)
                  id = 1
               end if
               if (b% rl_relative_gap(id) > 0.29858997d0*atan(1.83530121d0*pow(q,0.39661426d0))) then
                  binary_check_model = terminate
                  b% s_donor% termination_code = t_xtra1
                  termination_code_str(t_xtra1) = &
                      "Terminate because of L2 overflow"
               end if
            end if

            if (b% mtransfer_rate == -b% max_implicit_abs_mdot*Msun/secyer .and. &
               b% CE_begin_at_max_implicit_abs_mdot) then
               b% CE_flag = .true.
               b% CE_init = .false.
               binary_check_model = keep_going
            end if
         end if

      end function binary_check_model

      integer function binary_finish_step(b)
         type (binary_info), pointer :: b
         real(dp) :: spin_period

         binary_finish_step = keep_going
         ! update change factor in case mtransfer_rate has changed
         if(b% mtransfer_rate_old /= b% mtransfer_rate .and. &
             b% mtransfer_rate /= 0 .and. b% mtransfer_rate_old /= 0) then
            if(b% mtransfer_rate < b% mtransfer_rate_old) then
               b% change_factor = b% change_factor*(1d0-b% implicit_lambda) + b% implicit_lambda* &
                  (1+b% change_factor_fraction*(b% mtransfer_rate/b% mtransfer_rate_old-1))
            else
               b% change_factor = b% change_factor*(1d0-b% implicit_lambda) + b% implicit_lambda* &
                   (1+b% change_factor_fraction*(b% mtransfer_rate_old/b% mtransfer_rate-1))
            end if
            if(b% change_factor > b% max_change_factor) b% change_factor = b% max_change_factor
            if(b% change_factor < b% min_change_factor) b% change_factor = b% min_change_factor
         end if

         ! store all variables into "old"

         b% model_number_old = b% model_number
         b% binary_age_old = b% binary_age
         b% mtransfer_rate_old = b% mtransfer_rate
         b% angular_momentum_j_old = b% angular_momentum_j
         b% separation_old = b% separation
         b% eccentricity_old = b% eccentricity
         b% dt_old = b% dt
         b% env_old(1) = b% env(1)
         b% env_old(2) = b% env(2)
         b% period_old = b% period
         b% rl_relative_gap_old(1) = b% rl_relative_gap(1)
         b% rl_relative_gap_old(2) = b% rl_relative_gap(2)
         b% r_old(1) = b% r(1)
         b% r_old(2) = b% r(2)
         b% rl_old(1) = b% rl(1)
         b% rl_old(2) = b% rl(2)
         b% m_old(1) = b% m(1)
         b% m_old(2) = b% m(2)
         b% using_jdot_mb_old = b% using_jdot_mb
         b% max_timestep_old = b% max_timestep
         b% change_factor_old = b% change_factor

         b% d_i_old = b% d_i
         b% a_i_old = b% a_i
         b% point_mass_i_old = b% point_mass_i

         b% ignore_rlof_flag_old = b% ignore_rlof_flag
         b% model_twins_flag_old = b% model_twins_flag

         b% CE_flag_old = b% CE_flag
         b% CE_init_old = b% CE_init

         b% CE_num1_old = b% CE_num1
         b% CE_num2_old = b% CE_num2
         b% CE_lambda1_old = b% CE_lambda1
         b% CE_lambda2_old = b% CE_lambda2
         b% CE_Ebind1_old = b% CE_Ebind1
         b% CE_Ebind2_old = b% CE_Ebind2
         b% CE_years_detached_old = b% CE_years_detached

         b% dt_why_reason_old = b% dt_why_reason

         !set all xtra variables
         b% ixtra_old(:) = b% ixtra(:)
         b% xtra_old(:) = b% xtra(:)
         b% lxtra_old(:) = b% lxtra(:)

      end function binary_finish_step

      integer function binary_prepare_to_redo(b)
         type (binary_info), pointer :: b

         binary_prepare_to_redo = redo
         call binary_set_current_to_old(b)

      end function binary_prepare_to_redo

      integer function binary_prepare_to_retry(b)
         type (binary_info), pointer :: b

         b% num_tries = 0
         ! this call takes care of restoring variables
         call binary_set_current_to_old(b)
         binary_prepare_to_retry = retry

      end function binary_prepare_to_retry

      subroutine binary_set_current_to_old(b)
         type (binary_info), pointer :: b
         ! restore variables
         b% model_number = b% model_number_old
         b% binary_age = b% binary_age_old
         ! do not restore mtransfer_rate during implicit rlo
         if (b% num_tries == 0) b% mtransfer_rate = b% mtransfer_rate_old
         b% angular_momentum_j = b% angular_momentum_j_old
         b% separation = b% separation_old
         b% eccentricity = b% eccentricity_old
         b% dt = b% dt_old
         b% env(1) = b% env_old(1)
         b% env(2) = b% env_old(2)
         b% period = b% period_old
         b% rl_relative_gap(1) = b% rl_relative_gap_old(1)
         b% rl_relative_gap(2) = b% rl_relative_gap_old(2)
         b% r(1) = b% r_old(1)
         b% r(2) = b% r_old(2)
         b% rl(1) = b% rl_old(1)
         b% rl(2) = b% rl_old(2)
         b% m(1) = b% m_old(1)
         b% m(2) = b% m_old(2)
         b% using_jdot_mb = b% using_jdot_mb_old
         b% max_timestep = b% max_timestep_old

         ! the following need to be kept constant during implicit mdot
         if (b% num_tries == 0) then
            b% change_factor = b% change_factor_old
            b% d_i = b% d_i_old
            b% a_i = b% a_i_old
            b% point_mass_i = b% point_mass_i_old
            b% ignore_rlof_flag = b% ignore_rlof_flag_old
            b% model_twins_flag = b% model_twins_flag_old
            b% CE_flag = b% CE_flag_old
            b% CE_init = b% CE_init_old
         end if

         b% CE_num1 = b% CE_num1_old
         b% CE_num2 = b% CE_num2_old
         b% CE_lambda1 = b% CE_lambda1_old
         b% CE_lambda2 = b% CE_lambda2_old
         b% CE_Ebind1 = b% CE_Ebind1_old
         b% CE_Ebind2 = b% CE_Ebind2_old
         b% CE_years_detached = b% CE_years_detached_old

         b% dt_why_reason = b% dt_why_reason_old

         !set all xtra variables
         b% ixtra(:) = b% ixtra_old(:)
         b% xtra(:) = b% xtra_old(:)
         b% lxtra(:) = b% lxtra_old(:)
      end subroutine binary_set_current_to_old

      integer function binary_after_evolve(b)
         type (binary_info), pointer :: b
         binary_after_evolve = keep_going
         
         !take care of deallocating binary arrays here
         if (associated(b% theta_co)) then
            deallocate(b% theta_co)
            nullify(b% theta_co)
         end if
         if (associated(b% time_co)) then
            deallocate(b% time_co)
            nullify(b% time_co)
         end if
         if (associated(b% mdot_donor_theta)) then
            deallocate(b% mdot_donor_theta)
            nullify(b% mdot_donor_theta)
         end if
         if (associated(b% edot_theta)) then
            deallocate(b% edot_theta)
            nullify(b% edot_theta)
         end if
         if (associated(b% e1)) then
            deallocate(b% e1)
            nullify(b% e1)
         end if
         if (associated(b% e2)) then
            deallocate(b% e2)
            nullify(b% e2)
         end if
         if (associated(b% CE_m)) then
            deallocate(b% CE_m)
            nullify(b% CE_m)
         end if
         if (associated(b% CE_entropy)) then
            deallocate(b% CE_entropy)
            nullify(b% CE_entropy)
         end if
         if (associated(b% CE_U_in)) then
            deallocate(b% CE_U_in)
            nullify(b% CE_U_in)
         end if
         if (associated(b% CE_U_out)) then
            deallocate(b% CE_U_out)
            nullify(b% CE_U_out)
         end if
         if (associated(b% CE_Omega_in)) then
            deallocate(b% CE_Omega_in)
            nullify(b% CE_Omega_in)
         end if
         if (associated(b% CE_Omega_out)) then
            deallocate(b% CE_Omega_out)
            nullify(b% CE_Omega_out)
         end if
      end function binary_after_evolve

      end module binary_evolve
