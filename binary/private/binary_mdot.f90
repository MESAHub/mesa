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


      module binary_mdot

      use const_def
      use math_lib
      use star_lib
      use star_def
      use binary_def
      use binary_wind
      use binary_ce
      use utils_lib, only: mesa_error

      implicit none

      contains

      integer function check_implicit_rlo(binary_id, new_mdot)
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: new_mdot
         
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: function_to_solve, explicit_mdot, q, slope_contact
         integer :: ierr
         logical :: use_sum
         character (len=90) :: rlo_result
         
         include 'formats'
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         s => b% s_donor
         use_sum = .false.
         
         ! NOTE: keep in mind that for mass loss, mdot is negative.
         ! b% mtransfer_rate will be considered valid if function_to_solve = 0
         ! within the tolerance given by b% implicit_scheme_tolerance, i.e.
         ! |function_to_solve| < b% implicit_scheme_tolerance.
         !
         ! For the roche_lobe scheme, function_to_solve is set such that solutions
         ! will satisfy 0 > (r-rl)/rl > - b% implicit_scheme_tolerance. If
         ! (r-rl)/rl < 0
         ! and the new mass transfer rate that was attempted satisfies
         ! |new_mdot| < b% roche_min_mdot*Msun/secyer
         ! then the system is considered to be detached and mass loss is turned off.
         !
         ! For other schemes, function_to_solve is chosen as the difference between
         ! b% mtransfer_rate and the explicit transfer rate, divided by the
         ! explicit transfer rate.
         
         check_implicit_rlo = keep_going
         new_mdot = b% mtransfer_rate
         b% num_tries = b% num_tries + 1
         rlo_result = "missing message"

         if(b% report_rlo_solver_progress .and. b% num_tries == 1) then
            write(*,'(a)') &
               '  binary_step rlo_iter di ch_fact  curr_mdot' // &
               '  next_mdot    mdot_up   mdot_low    rl_gap1    rl_gap2          f              result'
            write(*,'(a)') &
               '  __________________________________________' // &
               '______________________________________________________________________________________'
         end if

         if (b% use_other_implicit_function_to_solve) then
            call b% other_implicit_function_to_solve(b% binary_id, function_to_solve, use_sum, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in other_implicit_function_to_solve'
               check_implicit_rlo = retry
               return
            end if
         else if (b% mdot_scheme == "roche_lobe" .and. .not. b% use_other_rlo_mdot) then
            function_to_solve = (b% rl_relative_gap(b% d_i) &
                + b% implicit_scheme_tolerance/2.0d0) * 2.0d0

            if (function_to_solve < 0 .and. abs(b% mtransfer_rate) == 0) then
               if (b% report_rlo_solver_progress) then
                  rlo_result = 'OK (detached)'
                  call report_rlo_iter
               end if
               return
            end if
         else if (b% mdot_scheme == "contact" .and. .not. b% use_other_rlo_mdot .and. .not. b% CE_flag) then
            if (b% point_mass_i /= 0) then
               new_mdot = 0d0
               write(*,*) "WARNING: contact scheme requires evolve_both_stars=.true."
               write(*,*) "Not transfering mass"
               return
            end if
            q = b% m(b% a_i) / b% m(b% d_i)
            slope_contact = pow(q, -0.52d0)
            ! If accretor is overflowing its Roche lobe, then the contact scheme needs to be used.
            ! Otherwise, if accretor radius is (within tolerance) below the equipotential
            ! of the donor, or donor is below tolerance for detachment, then use regular roche_lobe scheme.
            if (b% rl_relative_gap(b% a_i) < 0 .and. & 
                (b% rl_relative_gap(b% d_i)*slope_contact - b% rl_relative_gap(b% a_i) &
                 > b% implicit_scheme_tolerance .or. &
                 b% rl_relative_gap(b% d_i) < - b% implicit_scheme_tolerance)) then

               function_to_solve = (b% rl_relative_gap(b% d_i) &
                   + b% implicit_scheme_tolerance/2.0d0) * 2.0d0

               if (function_to_solve < 0 .and. abs(b% mtransfer_rate) == 0) then
                  if (b% report_rlo_solver_progress) then
                     rlo_result = 'OK (detached)'
                     call report_rlo_iter
                  end if
                  return
               end if
            else
               if (q < 1d0) then
                  function_to_solve = (b% rl_relative_gap(b% d_i)*slope_contact - b% rl_relative_gap(b% a_i))
               else
                  function_to_solve = (b% rl_relative_gap(b% d_i) - b% rl_relative_gap(b% a_i)/slope_contact)
               end if
               use_sum = .true.
            end if
         else
            if (b% CE_flag) then
               call CE_rlo_mdot(b% binary_id, explicit_mdot, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in CE_rlo_mdot'
                  check_implicit_rlo = retry
                  return
               end if
            else if (.not. b% use_other_rlo_mdot) then
               call rlo_mdot(b% binary_id, explicit_mdot, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in rlo_mdot'
                  check_implicit_rlo = retry
                  return
               end if
            else
               call b% other_rlo_mdot(b% binary_id, explicit_mdot, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in other_rlo_mdot'
                  check_implicit_rlo = retry
                  return
               end if
            end if
            if (explicit_mdot > 0d0) then
               new_mdot = 0d0
               write(*,*) "WARNING: explicit computation of mass transfer results in accreting donor"
               write(*,*) "Not transferring mass"
               return
            end if
            function_to_solve = 0d0
            if(explicit_mdot /= 0d0) then
               function_to_solve = (explicit_mdot - new_mdot) / explicit_mdot
            else if (new_mdot /= 0d0) then
               function_to_solve = (explicit_mdot - new_mdot) / new_mdot
            end if
            if (abs(explicit_mdot) <= b% min_mdot_for_implicit*Msun/secyer .and. &
               abs(new_mdot) <= b% min_mdot_for_implicit*Msun/secyer) then
               new_mdot = explicit_mdot
               if (b% report_rlo_solver_progress) then
                  rlo_result = 'OK (explicit)'
                  call report_rlo_iter
               end if
               return
            end if
         end if
         
         if (abs(function_to_solve) <= b% implicit_scheme_tolerance) then
            if (b% report_rlo_solver_progress) then
               rlo_result = 'OK (in tol)'
               call report_rlo_iter
               write(*,'(A)')
            end if
            return
         end if
         
         if (b% num_tries > b% max_tries_to_achieve) then
            check_implicit_rlo = retry
            if (b% report_rlo_solver_progress) then
               rlo_result = 'retry (>max tries)'
               call report_rlo_iter
            end if
            return
         end if
            
         if (b% num_tries == 1) then
            b% have_mdot_lo = .false.
            b% have_mdot_hi = .false.
            b% mdot_lo = 0
            b% mdot_hi = 0
            b% fixed_delta_mdot = b% mtransfer_rate * (1 - b% change_factor)
         end if
         
         new_mdot = pick_mdot_for_implicit_rlo(b, function_to_solve, b% mtransfer_rate, use_sum, ierr)

         !if this iteration is done using the maximum mass transfer rate,
         !and the mass transfer rate increases in the following iteration,
         !then stop at this point and accept the step
         !NOTE: both b% mtransfer_rate and new_mdot are negative
         if (b% mtransfer_rate == -b% max_implicit_abs_mdot*Msun/secyer &
            .and. b% mtransfer_rate > new_mdot) then
            if (b% report_rlo_solver_progress) then
               rlo_result = 'OK (max mdot)'
               call report_rlo_iter
            end if
            new_mdot = - b% max_implicit_abs_mdot*Msun/secyer
            return
         end if
         if (-new_mdot > b% max_implicit_abs_mdot*Msun/secyer) then
            !limit to maximum mdot, following iteration will be made
            !using this value
            new_mdot = - b% max_implicit_abs_mdot*Msun/secyer
         end if

         if (ierr /= 0) then
            check_implicit_rlo = retry
            if (b% report_rlo_solver_progress) then
               rlo_result = 'retry (ierr)'
               call report_rlo_iter
            end if
            return
         end if
         
         if (-new_mdot < b% roche_min_mdot*Msun/secyer .and. function_to_solve < 0 .and. &
             (b% mdot_scheme == "roche_lobe" .or. (b% mdot_scheme == "contact" .and. &
             .not. use_sum))) then
            check_implicit_rlo = keep_going
            new_mdot = 0d0
            if (b% report_rlo_solver_progress) then
               rlo_result = 'OK (detachment)'
               call report_rlo_iter
            end if
            return
         end if
         
         if (b% have_mdot_hi .and. b% have_mdot_lo) then
            if (abs(b% mdot_hi - b% mdot_lo) < &
                  b% implicit_scheme_tiny_factor*min(abs(b% mdot_hi),abs(b% mdot_lo))) then
               if (b% report_rlo_solver_progress) then
                  rlo_result = 'OK (tiny change)'
                  call report_rlo_iter
               end if
               new_mdot = b% mtransfer_rate
               check_implicit_rlo = keep_going
               return
            end if
         end if

         ! boost change_factor every num_tries_for_increase_change_factor tries, to avoid
         ! implicit scheme from becoming stuck with sudden changes.
         if (b% num_tries_for_increase_change_factor > 0 .and. b% change_factor_increase > 1d0 &
            .and. .not. (b% have_mdot_lo .and. b% have_mdot_hi)) then
            if(mod(b% num_tries, b% num_tries_for_increase_change_factor) == 0) then
               b% change_factor = min(b% max_change_factor, &
                  b% change_factor_increase * b% change_factor)
            end if
         end if
         if (b% report_rlo_solver_progress) then
            rlo_result = 'redo'
            call report_rlo_iter
         end if
         
         check_implicit_rlo = redo

         contains

         subroutine report_rlo_iter
            real(dp) :: mdot_hi, mdot_lo
            if (.not. b% have_mdot_hi) then
               mdot_hi = huge(mdot_hi)
            else
               mdot_hi = -b% mdot_hi/Msun*secyer
            end if
            if (.not. b% have_mdot_lo) then
               mdot_lo = -huge(mdot_lo)
            else
               mdot_lo = -b% mdot_lo/Msun*secyer
            end if
            write(*,'(i13,i9,i3,f8.4,7(1pe11.3,0p),a20)') &
               b% model_number, &
               b% num_tries, &
               b% d_i, &
               b% change_factor, &
               -b% mtransfer_rate/Msun*secyer, &
               -new_mdot/Msun*secyer, &
               mdot_hi, &
               mdot_lo, &
               b% rl_relative_gap(1), &
               b% rl_relative_gap(2), &
               function_to_solve, &
               trim(rlo_result)
         end subroutine report_rlo_iter

      end function check_implicit_rlo

      real(dp) function pick_mdot_for_implicit_rlo( &
            b, new_function_to_solve, mdot_current, use_sum, ierr) result(mdot_next)
         use num_lib, only: find0_quadratic, find0
         type(binary_info), pointer :: b
         real(dp), intent(in) :: new_function_to_solve, mdot_current
         logical, intent(in) :: use_sum
         integer, intent(out) :: ierr
         
         real(dp) :: starting_mdot, current_change_factor
         logical :: do_cubic
         include 'formats'
         
         ! NOTE: keep in mind that for mass loss, mdot is negative
         
         
         starting_mdot = -b% starting_mdot*Msun/secyer
         current_change_factor = pow(b% change_factor, b% num_tries+1)

         if (.not.(b% solver_type == "cubic" .or. &
                   b% solver_type == "bisect" .or. &
                   b% solver_type == "both")) then
            write(*,*) "ERROR: unrecognized b% solver_type", trim(b% solver_type)
            write(*,*) "setting b% solver_type to 'both'"
            b% solver_type = "both"
         end if

         if (b% have_mdot_lo .and. b% have_mdot_hi) then
            ierr = 0
            do_cubic = (mod(b% num_tries,2)==0 .and. b% solver_type == "both") .or. &
               b% solver_type == "cubic"
            if (do_cubic) then
               mdot_next = find0_quadratic( &
                   b% mdot_lo, b% implicit_function_lo, &
                   mdot_current, new_function_to_solve, &
                   b% mdot_hi, b% implicit_function_hi, ierr)
               if (ierr /= 0) then
                  mdot_next = find0(b% mdot_lo, b% implicit_function_lo, &
                      b% mdot_hi, b% implicit_function_hi)
                  ierr = 0
               end if
            end if
            if (new_function_to_solve >= 0) then
               b% mdot_lo = mdot_current
               b% implicit_function_lo = new_function_to_solve
            else
               b% mdot_hi = mdot_current
               b% implicit_function_hi = new_function_to_solve
            end if
            if (.not. do_cubic) then
               mdot_next = (b% mdot_hi+b% mdot_lo)/2.0d0
            end if
         else if (b% have_mdot_lo) then ! don't have mdot_hi
            if (new_function_to_solve < 0) then
               b% mdot_hi = mdot_current
               b% implicit_function_hi = new_function_to_solve
               b% have_mdot_hi = .true.
               mdot_next = (b% mdot_hi+b% mdot_lo)/2.0d0
               !mdot_next = find0(b% mdot_lo, b% implicit_function_lo, &
               !    b% mdot_hi, b% implicit_function_hi)
            else ! still too low
               b% mdot_lo = mdot_current
               b% implicit_function_lo = new_function_to_solve
               mdot_next = b% mdot_lo*b% change_factor
            end if
         else if (b% have_mdot_hi) then ! don't have mdot_lo
            if (new_function_to_solve >= 0) then
               b% mdot_lo = mdot_current
               b% implicit_function_lo = new_function_to_solve
               b% have_mdot_lo = .true.
               mdot_next = (b% mdot_hi+b% mdot_lo)/2.0d0
               !mdot_next = find0(b% mdot_lo, b% implicit_function_lo, &
               !    b% mdot_hi, b% implicit_function_hi)
            else ! mdot still too high
               b% mdot_hi = mdot_current
               b% implicit_function_hi = new_function_to_solve
               if (use_sum .and. mod(b% num_tries, 2) == 0) then
                  mdot_next = b% mdot_hi + b% fixed_delta_mdot
               else
                  mdot_next = b% mdot_hi/b% change_factor
               end if
            end if
         else ! don't have either
            if (mdot_current > starting_mdot .and. (.not. abs(mdot_current) > 0)) then ! recall that both are negative
               mdot_next = starting_mdot
            else if (new_function_to_solve >= 0) then
               b% mdot_lo = mdot_current
               b% implicit_function_lo = new_function_to_solve
               b% have_mdot_lo = .true.
               mdot_next = b% mdot_lo*b% change_factor
            else
               b% mdot_hi = mdot_current
               b% implicit_function_hi = new_function_to_solve
               b% have_mdot_hi = .true.
               if (use_sum .and. mod(b% num_tries, 2) == 0) then
                  mdot_next = b% mdot_hi + b% fixed_delta_mdot
               else
                  mdot_next = b% mdot_hi/b% change_factor
               end if
            end if
         end if
         
      end function pick_mdot_for_implicit_rlo


      subroutine eval_mdot_edd(binary_id, mdot_edd, mdot_edd_eta, ierr)
         use utils_lib, only: is_bad

         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot_edd ! eddington accretion rate
         real(dp), intent(out) :: mdot_edd_eta ! fraction of rest mass energy released as radiation
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         include 'formats'

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         if (b% point_mass_i == 0) then
            if (b% limit_retention_by_mdot_edd) then
               write(*,*) "Default mdot_edd calculation cannot be used when evolving both stars"
               write(*,*) "Maybe you want to set limit_retention_by_mdot_edd=.false. in binary_controls?"
               write(*,*) "Setting mdot_edd to zero"
            end if
            mdot_edd = 0
            mdot_edd_eta = 0
            return
         end if

         if (b% use_this_for_mdot_edd_eta > 0) then
            mdot_edd_eta = b% use_this_for_mdot_edd_eta
         else
            ! eg., eq. (6) of Podsiadlowski, Rappaport & Han 2003, MNRAS, 341, 385
            mdot_edd_eta = 1d0 &
                 - sqrt(1d0 - pow2(min(b% m(b% a_i),sqrt(6d0)*b% eq_initial_bh_mass)/(3d0*b% eq_initial_bh_mass)))
         end if

         if (b% use_this_for_mdot_edd > 0) then
            mdot_edd = b% use_this_for_mdot_edd*(Msun/secyer)
         else
            ! eg., eq. (9) of Podsiadlowski, Rappaport & Han 2003, MNRAS, 341, 385
            if (.not. b% use_es_opacity_for_mdot_edd) then
               mdot_edd = pi4*standard_cgrav*b% m(b% a_i) &
                  /(clight*b% s_donor% opacity(1)*mdot_edd_eta)
            else
               mdot_edd = pi4*standard_cgrav*b% m(b% a_i)&
                  /(clight*0.2d0*(1d0+b% s_donor% surface_h1)*mdot_edd_eta)
            end if
         end if

         if (is_bad(mdot_edd_eta) .or. mdot_edd_eta<0) then
            write(*,*) "ERROR while computing mdot_edd_eta"
            ierr = -1
            write(*,*) "mdot_edd_eta, b% m(b% a_i), b% eq_initial_bh_mass", &
               mdot_edd_eta, b% m(b% a_i), b% eq_initial_bh_mass
            stop
         end if

         if (is_bad(mdot_edd) .or. mdot_edd<0) then
            write(*,*) "ERROR while computing mdot_edd"
            ierr = -1
            write(*,*) "mdot_edd, b% m(b% a_i), b% s_donor% opacity(1), b% s_donor% surface_h1", &
               mdot_edd, b% m(b% a_i), b% s_donor% opacity(1), b% s_donor% surface_h1
            stop
         end if

      end subroutine eval_mdot_edd

      subroutine adjust_mdots(b)
         use binary_wind, only: eval_wind_xfer_fractions
         type (binary_info), pointer :: b

         real(dp) :: fixed_xfer_fraction, actual_mtransfer_rate
         integer :: ierr

         actual_mtransfer_rate = 0d0
         
         if (b% use_other_adjust_mdots) then
            call b% other_adjust_mdots(b% binary_id, ierr)
            if (ierr /= 0) then
               write(*,*) "Error in other_adjust_mdots"
               stop
            end if
            return
         end if 

         b% fixed_xfer_fraction = 1 - b% mass_transfer_alpha - b% mass_transfer_beta - &
            b% mass_transfer_delta

         if (.not. b% use_other_mdot_edd) then
            call eval_mdot_edd(b% binary_id, b% mdot_edd, b% mdot_edd_eta, ierr)
         else
            call b% other_mdot_edd(b% binary_id, b% mdot_edd, b% mdot_edd_eta, ierr)
         end if

         ! Add tidal enhancement of wind
         call Tout_enhance_wind(b, b% s_donor)
         if (b% point_mass_i == 0) then
            ! do not repeat if using the implicit wind
            if (.not. (b% num_tries >0 .and. b% s_accretor% was_in_implicit_wind_limit)) &
               call Tout_enhance_wind(b, b% s_accretor)
         end if

         ! solve wind mass transfer
         ! b% mdot_wind_transfer(b% d_i) is a negative number that gives the
         ! amount of mass transferred by unit time from the donor to the
         ! accretor.
         call eval_wind_xfer_fractions(b% binary_id, ierr)
         if (ierr/=0) then
            write(*,*) "Error in eval_wind_xfer_fractions"
            return
         end if
         b% mdot_wind_transfer(b% d_i) = b% s_donor% mstar_dot * &
            b% wind_xfer_fraction(b% d_i)
         if (b% point_mass_i == 0) then
            b% mdot_wind_transfer(b% a_i) = b% s_accretor% mstar_dot * &
               b% wind_xfer_fraction(b% a_i)
         else
            b% mdot_wind_transfer(b% a_i) = 0d0
         end if

         ! Set mdot for the donor
         b% s_donor% mstar_dot = b% s_donor% mstar_dot + b% mtransfer_rate - &
            b% mdot_wind_transfer(b% a_i)

         ! Set mdot for the accretor
         if (b% point_mass_i == 0 .and. .not. b% CE_flag) then
            ! do not repeat if using the implicit wind
            if (.not. (b% num_tries >0 .and. b% s_accretor% was_in_implicit_wind_limit)) then
               b% accretion_mode = 0
               b% acc_am_div_kep_am = 0.0d0
               b% s_accretor% mstar_dot = b% s_accretor% mstar_dot - &
                  b% mtransfer_rate*b% fixed_xfer_fraction - b% mdot_wind_transfer(b% d_i)

               !set angular momentum accretion as described in A.3.3 of de Mink et al. 2013
               if (b% do_j_accretion) then
                  if (.not. b% use_other_accreted_material_j) then
                     call eval_accreted_material_j(b% binary_id, ierr)
                  else
                     call b% other_accreted_material_j(b% binary_id, ierr)
                  end if
                  if (ierr /= 0) then
                     write(*,*) 'error in accreted_material_j'
                     return
                  end if
               end if
            end if

            b% accretion_luminosity = 0d0 !only set for point mass

         else if (.not. b% CE_flag) then
            ! accretor is a point mass
            if (.not. b% model_twins_flag) then
               !combine wind and RLOF mass transfer
               actual_mtransfer_rate = b% mtransfer_rate*b% fixed_xfer_fraction+b% mdot_wind_transfer(b% d_i) !defined negative
               b% component_mdot(b% a_i) = -actual_mtransfer_rate
               ! restrict accretion to the Eddington limit
               if (b% limit_retention_by_mdot_edd .and. b% component_mdot(b% a_i) > b% mdot_edd) then
                  b% component_mdot(b% a_i) = b% mdot_edd ! remove all accretion above the edd limit
               end if
               b% accretion_luminosity = &
                  b% mdot_edd_eta*b% component_mdot(b% a_i)*clight*clight
               ! remove rest mass radiated away
               if (b% use_radiation_corrected_transfer_rate) then
                  b% component_mdot(b% a_i) = (1 - b% mdot_edd_eta) * b% component_mdot(b% a_i)
               end if
            end if
         else
            !doing CE, just be sure to set mdot for a point mass to zero
            b % accretion_luminosity = 0d0
            if (b% point_mass_i /= 0) then
               b% component_mdot(b% a_i) = 0d0
            end if
         end if

         ! mdot_system_transfer is mass lost from the vicinity of each star
         ! due to inefficient rlof mass transfer, mdot_system_cct is mass lost
         ! from a circumbinary coplanar toroid.
         if (b% mtransfer_rate+b% mdot_wind_transfer(b% d_i) >= 0 .or. b% CE_flag) then
            b% mdot_system_transfer(b% d_i) = 0d0
            b% mdot_system_transfer(b% a_i) = 0d0
            b% mdot_system_cct = 0d0
         else 
            b% mdot_system_transfer(b% d_i) = b% mtransfer_rate * b% mass_transfer_alpha
            b% mdot_system_cct = b% mtransfer_rate * b% mass_transfer_delta
            if (b% point_mass_i == 0 .or. b% model_twins_flag) then
               b% mdot_system_transfer(b% a_i) = b% mtransfer_rate * b% mass_transfer_beta
            else
               ! do not compute mass lost from the vicinity using just mass_transfer_beta, as
               ! mass transfer can be stopped also by going past the Eddington limit
               b% mdot_system_transfer(b% a_i) = (actual_mtransfer_rate + b% component_mdot(b% a_i)) &
                  + b% mtransfer_rate * b% mass_transfer_beta
            end if
         end if

      end subroutine adjust_mdots

      subroutine rlo_mdot(binary_id, mdot, ierr) ! Adapted from a routine kindly provided by Anastasios Fragkos
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b

         include 'formats'

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         mdot = 0d0

         if (b% mdot_scheme == "roche_lobe") then
            write(*,*) "mdot_scheme = roche_lobe not applicable for explicit scheme"
            write(*,*) "Not transfering mass"
            mdot = 0
            return
         else if (b% mdot_scheme /= "Ritter" .and. b% mdot_scheme /= "Kolb" .and. b% mdot_scheme /= "Arras") then
            write(*,*) "mdot_scheme = " , b% mdot_scheme , " not recognized"
            write(*,*) "Not transfering mass"
            mdot = 0
            return
         end if

         if (b% mdot_scheme == "Kolb" .and. b% eccentricity <= 0.0d0) then
            call get_info_for_ritter(b)
            mdot = b% mdot_thin
            call get_info_for_kolb(b)
            mdot = mdot + b% mdot_thick
            
         else if (b% mdot_scheme == "Kolb" .and. b% eccentricity > 0.0d0) then
            call get_info_for_ritter_eccentric(b)
            mdot = b% mdot_thin
            call get_info_for_kolb_eccentric(b)
            mdot = mdot + b% mdot_thick
            
         else if (b% mdot_scheme == "Ritter" .and. b% eccentricity <= 0.0d0) then
            call get_info_for_ritter(b)
            mdot = b% mdot_thin
            
         else if (b% mdot_scheme == "Ritter" .and. b% eccentricity > 0.0d0) then
            call get_info_for_ritter_eccentric(b)
            mdot = b% mdot_thin

         end if
            
         if (b% mdot_scheme == "Arras") then
            if (b% eccentricity > 0d0) &
               write(*,*) "mdot_scheme = Arras is not properly implemented for e>0"
            call get_info_for_arras(b)
            mdot = b% mdot_thin

         end if

      end subroutine rlo_mdot

      subroutine get_info_for_ritter(b)
         type(binary_info), pointer :: b
         real(dp) :: rho_exponent, F1, q, rho, p, grav, hp, v_th, rl3, q_temp
         include 'formats'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% Peos(1) ! pressure at surface in dynes/cm^2
         grav = standard_cgrav*b% m(b% d_i)/pow2(b% r(b% d_i)) ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))

         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10(q_temp))
         rl3 = (b% rl(b% d_i))*(b% rl(b% d_i))*(b% rl(b% d_i))
         b% mdot_thin0 = (2.0D0*pi/exp(0.5d0)) * v_th*v_th*v_th * &
             rl3/(standard_cgrav*b% m(b% d_i)) * rho * F1
         !Once again, do not extrapolate! Eq. (7) of Ritter 1988
         q_temp = min(max(q,0.04d0),20d0)
         if (q_temp < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10(q_temp) - 0.038D0*pow2(log10(q_temp)) )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10(q_temp) + 0.114D0*pow2(log10(q_temp)) )
         end if

         b% ritter_exponent = (b% r(b% d_i)-b% rl(b% d_i))/b% ritter_h

         if (b% mdot_scheme == "Kolb") then
            if (b% ritter_exponent > 0) then
               b% mdot_thin = -b% mdot_thin0
            else
               b% mdot_thin = -b% mdot_thin0 * exp(b% ritter_exponent)
            end if
         else
            b% mdot_thin = -b% mdot_thin0 * exp(b% ritter_exponent)
         end if

      end subroutine get_info_for_ritter
      
      real(dp) function calculate_kolb_mdot_thick(b, indexR, rl_d) result(mdot_thick)
         real(dp), intent(in) :: rl_d
         integer, intent(in) :: indexR
         real(dp) :: F1, F3, G1, dP, q, rho, p, grav, hp, v_th, rl3, q_temp
         integer :: i
         type(binary_info), pointer :: b
         include 'formats'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! As described in Kolb and H. Ritter 1990, A&A 236,385-392
         
         ! compute integral in Eq. (A17 of Kolb & Ritter 1990)
         mdot_thick = 0d0
         do i=1,indexR-1
            G1 = b% s_donor% gamma1(i)
            F3 = sqrt(G1) * pow(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
            mdot_thick = mdot_thick + F3*sqrt(kerg * b% s_donor% T(i) / &
               (mp * b% s_donor% mu(i)))*(b% s_donor% Peos(i+1)-b% s_donor% Peos(i))
         end do
         ! only take a fraction of dP for last cell 
         G1 = b% s_donor% gamma1(i)
         F3 = sqrt(G1) * pow(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
         dP = (b% s_donor% r(indexR) - rl_d) / &
            (b% s_donor% r(indexR) - b% s_donor% r(indexR+1)) * (b% s_donor% Peos(i+1)-b% s_donor% Peos(i))
         mdot_thick = mdot_thick + F3*sqrt(kerg * b% s_donor% T(i) / (mp*b% s_donor% mu(i)))*dP

         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10(q_temp))
         mdot_thick = -2.0D0*pi*F1*rl_d*rl_d*rl_d/(standard_cgrav*b% m(b% d_i))*mdot_thick
      
      end function calculate_kolb_mdot_thick
      
      subroutine get_info_for_kolb(b)
         type(binary_info), pointer :: b
         real(dp) :: F3, FF, G1, x_L1, q, g
         real(dp) :: mdot_thick0,  R_gas, dP, rl, s_div_rl
         integer :: i, indexR
         include 'formats'

         !--------------------- Optically thick MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         ! First we need to find how deep inside the star the Roche lobe reaches. In other words the mesh point of the star at which R=R_RL
         b% mdot_thick = 0d0
         indexR=-1
         if(b% r(b% d_i)-b% rl(b% d_i) > 0.0d0) then
            i=1
            do while (b% s_donor% r(i) > b% rl(b% d_i))
               i=i+1
            end do
            
            if (i .eq. 1) then
               b% mdot_thick = 0d0
            else
               b% mdot_thick = calculate_kolb_mdot_thick(b, i-1, b% rl(b% d_i))
            end if
         end if

      end subroutine get_info_for_kolb

      subroutine get_info_for_arras(b)
         type(binary_info), pointer :: b
         real(dp) :: q, rho, p, grav, hp, v_th
         real(dp) :: area,Asl,G,ma,md,mfac1,mfac2,my_mdot_thin,my_ritter_exponent,Omega,&
            phi,phiL1,q13,rfac,rhoL1,rv,rvL1,sep
         include 'formats'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! Ritter 1988 but with better fits for the various formulas that work at extreme q

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% Peos(1) ! pressure at surface in dynes/cm^2
         grav = standard_cgrav*b% m(b% d_i)/pow2(b% r(b% d_i)) ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))

         q = b% m(b% a_i) / b% m(b% d_i)
         G = standard_cgrav
         md = b% m(b% d_i)
         ma = b% m(b% a_i)
         sep = b% separation
         Omega = 2.d0*pi / b% period
         rvL1 = b% rl(b% d_i)
         rv = b% r(b% d_i)
         mfac1 = 1d0 + ma/md ! (md+ma)/md
         ! mfac2 = ( (md+ma)**2 + 3d0*ma*(md+ma) + 9d0*ma**2 ) / md**2
         mfac2 = 1d0 + 5d0*ma/md + 13d0*ma*ma/(md*md)
         rfac=rvL1/sep
         phiL1 = - G*md/rvL1 &
              * ( 1.d0 +  mfac1*pow3(rfac)/3.d0 + 4.d0*mfac2*pow6(rfac)/45.d0  )
              rfac=rv/sep
         phi = - G*md/rv &
              * ( 1.d0 +  mfac1*pow3(rfac)/3.d0 + 4.d0*mfac2*pow6(rfac)/45.d0  )
         my_ritter_exponent = - (phiL1-phi)/(v_th*v_th)
         rhoL1 = rho/sqrt(exp(1.d0)) * exp( my_ritter_exponent )
              q13=pow(q,one_third)
              Asl = 4.d0 + 4.16d0/(-0.96d0 + q13 + 1.d0/q13)
         area = 2.d0 * pi * pow2(v_th/Omega) / sqrt( Asl*(Asl-1d0) )
         my_mdot_thin = - rhoL1 * v_th * area
         b% mdot_thin = my_mdot_thin

      end subroutine get_info_for_arras

      subroutine get_info_for_ritter_eccentric(b)
         type(binary_info), pointer :: b
         integer :: i
         real(dp) :: rho_exponent, F1, q, q_temp, rho, p, grav, hp, v_th, dm
         real(dp), DIMENSION(b% anomaly_steps):: mdot0, mdot, Erit, rl_d
         include 'formats'
         
         ! Optically thin MT rate adapted for eccentric orbits 
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% Peos(1) ! pressure at surface in dynes/cm^2
         grav = standard_cgrav*b% m(b% d_i)/pow2(b% r(b% d_i)) ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1))) ! kerg = Boltzmann's constant
         
         ! phase dependant RL radius
         do i = 1, b% anomaly_steps
            rl_d(i) = b% rl(b% d_i) * (1d0 - pow2(b% eccentricity)) / &
                 (1 + b% eccentricity * cos(b% theta_co(i)) )
         end do

         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10(q_temp))

         mdot0 = (2.0D0*pi/exp(0.5d0)) * pow3(v_th) * rl_d*rl_d*rl_d / &
             (standard_cgrav*b% m(b% d_i)) * rho * F1   
             
         q_temp = min(max(q,0.04d0),20d0)
         if (q_temp < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10(q_temp) - 0.038D0*pow2(log10(q_temp)) )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10(q_temp) + 0.114D0*pow2(log10(q_temp)) )
         end if

         Erit = (b% r(b% d_i)- rl_d) / b% ritter_h

         if (b% mdot_scheme == "Kolb") then
            do i = 1,b% anomaly_steps
               if (Erit(i) > 0) then
                  mdot(i) = -1 * mdot0(i)
               else
                  mdot(i) = -1 * mdot0(i) * exp(Erit(i))
               end if
            end do
         else
            do i = 1,b% anomaly_steps
               mdot(i) = -1 * mdot0(i) * exp(Erit(i))
            end do
         end if
         
         b% mdot_donor_theta = mdot
         
         !integrate to get total massloss
         dm = 0d0
         do i = 2,b% anomaly_steps ! trapezoidal integration
            dm = dm + 0.5d0 * (mdot(i-1) + mdot(i)) * (b% time_co(i) - b% time_co(i-1)) 
         end do
         
         b% mdot_thin = dm

      end subroutine get_info_for_ritter_eccentric
      
      subroutine get_info_for_kolb_eccentric(b)
         type(binary_info), pointer :: b
         real(dp) :: e, dm
         integer :: i, j
         real(dp), DIMENSION(b% anomaly_steps):: rl_d_i, mdot_thick_i
         include 'formats'
         
         ! Optically thick MT rate adapted for eccentric orbits
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         b% mdot_thick = 0d0
         e = b% eccentricity
         
         ! If the radius of the donor is smaller as the smallest RL radius,
         ! there is only atmospheric RLOF, thus return.
         if ( b% r(b% d_i) < b% rl(b% d_i) * (1-e*e)/(1+e) ) then
            return
         end if
         
         ! For each point in the orbit calculate mdot_thick 
         do i = 1,b% anomaly_steps
            ! phase dependent RL radius
            rl_d_i(i) = b% rl(b% d_i) * (1d0 - e*e) / &
                 (1 + e*cos(b% theta_co(i)) )
         
            ! find how deep in the star we are
            j=1
            do while (b% s_donor% r(j) > rl_d_i(i))
               j=j+1
            end do
            
            ! calculate mdot_thick
            if (j .eq. 1) then
               mdot_thick_i(i) = 0d0
            else
               mdot_thick_i(i) = calculate_kolb_mdot_thick(b, j-1, rl_d_i(i))
            end if
         end do
         
         b% mdot_donor_theta = b% mdot_donor_theta + mdot_thick_i
         
         ! Integrate mdot_thick over the orbit
         dm = 0d0
         do i = 2,b% anomaly_steps ! trapezoidal integration
            dm = dm + 0.5d0 * (mdot_thick_i(i-1) + mdot_thick_i(i)) * &
                              (b% time_co(i) - b% time_co(i-1)) 
         end do
         
         b% mdot_thick = dm
         
      end subroutine get_info_for_kolb_eccentric

      subroutine eval_accreted_material_j(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         real(dp) :: qratio, min_r
         logical, parameter :: dbg = .false.
         include 'formats'

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         qratio = b% m(b% a_i) / b% m(b% d_i)
         qratio = min(max(qratio,0.0667d0),15d0)
         min_r = 0.0425d0*b% separation*pow(qratio+qratio*qratio, 0.25d0)

         !TODO: MUST USE EQUATORIAL RADIUS
         if (dbg) write(*,*) "radius, impact_radius, separation: ", &
             b% r(b% a_i), min_r/rsun, b% separation/rsun
         if (b% r(b% a_i) < min_r) then
            b% accretion_mode = 2
            b% s_accretor% accreted_material_j = &
               sqrt(standard_cgrav * b% m(b% a_i) * b% r(b% a_i)) 
         else
            b% accretion_mode = 1
            b% s_accretor% accreted_material_j = &
               sqrt(standard_cgrav * b% m(b% a_i) * 1.7d0*min_r)
         end if
         b% acc_am_div_kep_am = b% s_accretor% accreted_material_j / &
             sqrt(standard_cgrav * b% m(b% a_i) * b% r(b% a_i))

          !TODO: when using wind mass transfer donor star can end up
          ! with positive mdot, need to properly set jdot in that case

      end subroutine eval_accreted_material_j

      subroutine set_accretion_composition(b, acc_index)
         use chem_def, only: chem_isos
         type (binary_info), pointer :: b
         integer, intent(in) :: acc_index ! index of star that gains mass

         integer j

         if (acc_index == b% a_i) then
            !set accreted material composition
            b% s_accretor% num_accretion_species = b% s_donor% species
            
            if(b% s_donor% species > size(b% s_accretor% accretion_species_id,dim=1)) then
               call mesa_error(__FILE__,__LINE__,'Nuclear network is too large for accretor, increase max_num_accretion_species')
            end if

            do j = 1, b% s_donor% species
               b% s_accretor% accretion_species_id(j) = chem_isos% name(b% s_donor% chem_id(j))
               b% s_accretor% accretion_species_xa(j) = b% s_donor% xa_removed(j)
            end do
         else
            ! also for the donor to account for wind mass transfer
            b% s_donor% num_accretion_species = b% s_accretor% species

            if(b% s_accretor% species > size(b% s_donor% accretion_species_id,dim=1)) then
               call mesa_error(__FILE__,__LINE__,'Nuclear network is too large for donor, increase max_num_accretion_species')
            end if

            do j = 1, b% s_accretor% species
               b% s_donor% accretion_species_id(j) = chem_isos% name(b% s_accretor% chem_id(j))
               b% s_donor% accretion_species_xa(j) = b% s_accretor% xa_removed(j)
            end do
         end if

      end subroutine set_accretion_composition

      end module binary_mdot
