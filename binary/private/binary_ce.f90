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


      module binary_ce

      use const_def
      use math_lib
      use star_lib
      use star_def
      use binary_def
      use interp_1d_def, only: pm_work_size
      use interp_1d_lib, only: interp_pm, interp_values, interp_value

      implicit none

      contains

      subroutine CE_init(b, restart, ierr)
         use chem_def, only: chem_isos
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interp_pm
         type (binary_info), pointer :: b
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp), pointer :: interp_work(:), adjusted_energy(:)
         integer :: i, k, op_err
         real(dp) :: rec_energy_HII_to_HI, &
                     rec_energy_HeII_to_HeI, &
                     rec_energy_HeIII_to_HeII, &
                     diss_energy_H2, &
                     frac_HI, frac_HII, &
                     frac_HeI, frac_HeII, frac_HeIII, &
                     avg_charge_He, energy_comp
         include 'formats'

         ! TODO: no care is taken in here when model_twins_flag is true
         ! not a priority, but needs to be sorted out whenever double core
         ! evolution is implemented

         ierr = 0

         if (b% use_other_CE_init) then
            call b% other_CE_init(b% binary_id, restart, ierr)
            return
         end if

         write(*,*) "TURNING ON CE"

         b% s_donor% mix_factor = 0d0
         b% s_donor% dxdt_nuc_factor = 0d0
         if (b% point_mass_i == 0) then
            b% s_accretor% mix_factor = 0d0
            b% s_accretor% dxdt_nuc_factor = 0d0
         end if

         if (b% d_i == 1) then
            b% CE_num1 = b% CE_num1 + 1
         else
            b% CE_num2 = b% CE_num2 + 1
         end if

         b% keep_donor_fixed = .true.

         b% CE_init = .true.

         if (.not. restart) then
            b% CE_years_detached = 0d0
            s => b% s_donor

            write(*,*) "Initiating common envelope phase"

            allocate(b% CE_m(s% nz), b% CE_entropy(4*s% nz), stat=ierr)
            if(ierr /= 0) then
               write(*,*) "CE_init: Error during allocation"
               return
            end if
            allocate(b% CE_U_in(4*s% nz), b% CE_U_out(4*s% nz), b% CE_Omega_in(4*s% nz), b% CE_Omega_out(4*s% nz), stat=ierr)
            if(ierr /= 0) then
               write(*,*) "CE_init: Error during allocation"
               return
            end if

            ! get energy from the EOS and adjust the different contributions from recombination/dissociation
            allocate(adjusted_energy(s% nz))
            do k=1, s% nz
               ! the following lines compute the fractions of HI, HII, HeI, HeII and HeIII
               ! things like ion_ifneut_H are defined in $MESA_DIR/ionization/public/ionization.def
               ! this file can be checked for additional ionization output available
               frac_HI = 0d0!get_ion_info(s,ion_ifneut_H,k)
               frac_HII = 1 - frac_HI

               ! ionization module provides neutral fraction and average charge of He.
               ! use these two to compute the mass fractions of HeI and HeII
               frac_HeI = 0d0!get_ion_info(s,ion_ifneut_He,k)
               avg_charge_He = 2d0!get_ion_info(s,ion_iZ_He,k)
               ! the following is the solution to the equations
               !   avg_charge_He = 2*fracHeIII + 1*fracHeII
               !               1 = fracHeI + fracHeII + fracHeIII
               frac_HeII = 2 - 2*frac_HeI - avg_charge_He
               frac_HeIII = 1 - frac_HeII - frac_HeI

               ! recombination energies from https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
               rec_energy_HII_to_HI = avo*13.59843449d0*frac_HII*ev2erg*s% X(k)
               diss_energy_H2 = avo*4.52d0/2d0*ev2erg*s% X(k)
               rec_energy_HeII_to_HeI = avo*24.58738880d0*(frac_HeII+frac_HeIII)*ev2erg*s% Y(k)/4d0
               rec_energy_HeIII_to_HeII = avo*54.4177650d0*frac_HeIII*ev2erg*s% Y(k)/4d0

               adjusted_energy(k) = s% energy(k) &
                                    - (1d0-b% CE_energy_factor_HII_toHI)*rec_energy_HII_to_HI &
                                    - (1d0-b% CE_energy_factor_HeII_toHeI)*rec_energy_HeII_to_HeI &
                                    - (1d0-b% CE_energy_factor_HeIII_toHeII)*rec_energy_HeIII_to_HeII &
                                    - (1d0-b% CE_energy_factor_H2)*diss_energy_H2

               if (adjusted_energy(k) < 0d0 .or. adjusted_energy(k) > s% energy(k)) then
                  write(*,*) "Error when computing adjusted energy in CE, ", &
                     "s% energy(k):", s% energy(k), " adjusted_energy, ", adjusted_energy(k)
                  stop
               end if

               if(.false.) then
                  ! for debug, check the mismatch between the EOS energy and that of a gas+radiation
                  energy_comp = 3*avo*boltzm*s% T(k)/(2*s% mu(k)) + crad*pow4(s% T(k))/s% rho(k) &
                                + rec_energy_HII_to_HI &
                                + rec_energy_HeII_to_HeI &
                                + rec_energy_HeIII_to_HeII &
                                + diss_energy_H2

                  write(*,*) "compare energies", k, s%m(k)/Msun, s% energy(k), energy_comp, &
                     (s% energy(k)-energy_comp)/s% energy(k)
               end if

            end do

            do k=1, s% nz
               b% CE_m(:) = s% m(:s% nz)
            end do

            ! setup values of starting model for the interpolant
            do k=1, s% nz
               b% CE_entropy(4*k-3) = exp(s% lnS(k))
            end do

            ! Compute internal and potential energies from the inside out, and in the opposite direction.
            b% CE_U_out(1) = adjusted_energy(1)*s% dm(1)
            b% CE_Omega_out(1) = - standard_cgrav*s% m(1)*s% dm_bar(1)/s% r(1)
            do k=2, s% nz
               b% CE_U_out(4*k-3) = b% CE_U_out(4*(k-1)-3) + adjusted_energy(k)*s% dm(k)
               b% CE_Omega_out(4*k-3) = b% CE_Omega_out(4*(k-1)-3) - standard_cgrav*s% m(k)*s% dm_bar(k)/s% r(k)
            end do
            b% CE_U_in(4*s% nz-3) = adjusted_energy(s% nz)*s% dm(s% nz)
            b% CE_Omega_in(4*s% nz-3) = - standard_cgrav*s% m(s% nz)*s% dm_bar(s% nz)/s% r(s% nz)
            do k=s% nz-1, 1, -1
               b% CE_U_in(4*k-3) = b% CE_U_in(4*(k+1)-3) + adjusted_energy(k)*s% dm(k)
               b% CE_Omega_in(4*k-3) = b% CE_Omega_in(4*(k+1)-3) - standard_cgrav*s% m(k)*s% dm_bar(k)/s% r(k)
            end do

            b% CE_initial_radius = b% r(b% d_i)
            b% CE_initial_separation = b% separation
            b% CE_initial_Mdonor = b% m(b% d_i)
            b% CE_initial_Maccretor = b% m(b% a_i)
            b% CE_initial_age = s% star_age
            b% CE_initial_model_number = s% model_number
            b% CE_b_initial_age = b% binary_age
            b% CE_b_initial_model_number = b% model_number
            b% CE_nz = s% nz

            allocate(interp_work(s% nz*pm_work_size), stat=ierr)
            call interp_pm(b% CE_m, s% nz, b% CE_entropy, pm_work_size, interp_work, 'entropy interpolant', op_err)
            call interp_pm(b% CE_m, s% nz, b% CE_U_in, pm_work_size, interp_work, 'U_in interpolant', op_err)
            call interp_pm(b% CE_m, s% nz, b% CE_U_out, pm_work_size, interp_work, 'U_out interpolant', op_err)
            call interp_pm(b% CE_m, s% nz, b% CE_Omega_in, pm_work_size, interp_work, 'Omega_in interpolant', op_err)
            call interp_pm(b% CE_m, s% nz, b% CE_Omega_out, pm_work_size, interp_work, 'Omega_out interpolant', op_err)
            if(op_err /= 0) then
               ierr = -1
               write(*,*) "CE_init: Error while creating interpolants"
               return
            end if
            deallocate(adjusted_energy,interp_work)
         end if
          
      end subroutine

      subroutine CE_rlo_mdot(binary_id, rlo_mdot, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: rlo_mdot
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: exp_factor

         ierr = 0

         call binary_ptr(binary_id, b, ierr)

         if (b% use_other_CE_rlo_mdot) then
            call b% other_CE_rlo_mdot(b% binary_id, rlo_mdot, ierr)
            return
         end if

         exp_factor = -log(b% CE_mass_loss_rate_low/b% CE_mass_loss_rate_high)

         if (b% r(b% d_i)-b% rl(b% d_i) > 0d0) then
            rlo_mdot = -Msun/secyer*b% CE_mass_loss_rate_high
         else if (b% r(b% d_i)-b% rl(b% d_i) < -b% CE_rel_rlo_for_detachment*b% rl(b% d_i)) then
            rlo_mdot = -Msun/secyer*b% CE_mass_loss_rate_low
         else
            rlo_mdot = -Msun/secyer*b% CE_mass_loss_rate_high * &
               exp(exp_factor*(b% r(b% d_i)-b% rl(b% d_i))/(b% rl(b% d_i)*b% CE_rel_rlo_for_detachment))
         end if

      end subroutine CE_rlo_mdot

      integer function CE_binary_evolve_step(b)
         use binary_utils, only:set_separation_eccentricity
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: Ebind, Ecore, Ecore_i, lambda, &
            U_removed ,Omega_removed, U_inold, Omega_inold
         real(dp) :: separation, initial_Eorb
         integer :: ierr, op_err, k

         if (b% use_other_CE_binary_evolve_step) then
            CE_binary_evolve_step = b% other_CE_binary_evolve_step(b% binary_id)
            return
         end if

         ! setup mdot_system_wind to get output right, it is not actually used.
         ! setup jdots to zero as well
         b% mdot_system_wind(b% d_i) = b% s_donor% mstar_dot - b% mtransfer_rate
         if (b% point_mass_i == 0) then
            b% mdot_system_wind(b% a_i) = b% s_accretor% mstar_dot
         else
            b% mdot_system_wind(b% a_i) = 0.0d0
         end if
         b% jdot = 0d0
         b% jdot_mb = 0d0
         b% jdot_gr = 0d0
         b% jdot_ml = 0d0
         b% jdot_missing_wind = 0d0
         b% extra_jdot = 0d0
         b% jdot_ls = 0d0

         s => b% s_donor

         if (b% CE_fixed_lambda < 0d0) then
            Ecore = 0
            do k=1, s% nz
               Ecore = Ecore + s% energy(k)*s% dm(k) - standard_cgrav*s% m(k)*s% dm_bar(k)/s% r(k)
            end do

            call interp_value(b% CE_m, b% CE_nz, b% CE_U_out, s% m(1), U_removed, op_err)
            call interp_value(b% CE_m, b% CE_nz, b% CE_Omega_out, s% m(1), Omega_removed, op_err)
            call interp_value(b% CE_m, b% CE_nz, b% CE_U_in, s% m(1), U_inold, op_err)
            call interp_value(b% CE_m, b% CE_nz, b% CE_Omega_in, s% m(1), Omega_inold, op_err)

            Ecore_i = U_inold + Omega_inold
            Ebind = b% CE_alpha_th*U_removed + Omega_removed
            lambda = -(standard_cgrav*b% CE_initial_Mdonor*(b% CE_initial_Mdonor - s% m(1))) &
               /(Ebind*b% CE_initial_radius)
         else
            lambda = b% CE_fixed_lambda
            Ebind = -(standard_cgrav*b% CE_initial_Mdonor*(b% CE_initial_Mdonor - s% m(1))) &
               /(lambda*b% CE_initial_radius)
         end if

         if (b% d_i == 1) then
            b% CE_Ebind1 = Ebind
            b% CE_lambda1 = lambda
         else
            b% CE_Ebind2 = Ebind
            b% CE_lambda2 = lambda
         end if
         
         initial_Eorb = -standard_cgrav*b% CE_initial_Mdonor*b% CE_initial_Maccretor/(2*b% CE_initial_separation)

         separation = -b% CE_alpha*standard_cgrav*s% m(1)*b% CE_initial_Maccretor &
            /(2*(Ebind+b% CE_alpha*initial_Eorb))

         b% m(b% d_i) = b% s_donor% mstar
         b% time_step = b% s_donor% time_step
         if (b% point_mass_i == 0) then
            b% m(b% a_i) = b% s_accretor% mstar
         end if
         
         if (b% point_mass_i /= 1) then
            b% r(1) = Rsun*b% s1% photosphere_r
         else
            b% r(1) = 0
         end if
         if (b% point_mass_i /= 2) then
            b% r(2) = Rsun*b% s2% photosphere_r
         else
            b% r(2) = 0
         end if

         !only change separation if its reduced from the initial value
         call set_separation_eccentricity(b% binary_id, &
            min(b% CE_initial_separation, separation), b% eccentricity, ierr)

         b% model_number = b% model_number + 1
         b% time_step = b% s_donor% time_step
         b% binary_age = b% binary_age + b% time_step

         if (b% r(b% d_i)-b% rl(b% d_i) < 0d0) then
            b% CE_years_detached = b% CE_years_detached + b% time_step
         else
            b% CE_years_detached = 0d0
         end if

         CE_binary_evolve_step = keep_going

      end function CE_binary_evolve_step

      integer function CE_binary_finish_step(b)
         use binary_utils, only: eval_rlobe
         type (binary_info), pointer :: b
         real(dp) :: h_diff, he_diff, rlobe
         logical :: terminate_CE
         integer :: k
         CE_binary_finish_step = keep_going

         if (b% use_other_CE_binary_finish_step) then
            CE_binary_finish_step = b% other_CE_binary_finish_step(b% binary_id)
            return
         end if

         terminate_CE = .false.
         if ((b% r(b% d_i)-b% rl(b% d_i))/(b% rl(b% d_i)*b% CE_rel_rlo_for_detachment) < -1d0) then
            terminate_CE = .true.
            write(*,*) "Have reached CE_rel_rlo_for_detachment"
         else if (b% CE_years_detached > b% CE_years_detached_to_terminate) then
            terminate_CE = .true.
            write(*,*) "Have reached CE_years_detached_to_terminate"
         end if

         if (terminate_CE) then
            b% CE_flag = .false.
            b% mtransfer_rate = 0d0
            b% s_donor% mix_factor = 1d0
            b% s_donor% dxdt_nuc_factor = 1d0
            b% s_donor% timestep_hold = b% s_donor% model_number
            if (b% point_mass_i == 0) then
               b% s_accretor% mix_factor = 1d0
               b% s_accretor% dxdt_nuc_factor = 1d0
               b% s_accretor% timestep_hold = b% s_accretor% model_number
            end if

            b% keep_donor_fixed = .false.
            write(*,*) "TURNING OFF CE"
         end if

         ! termination conditions

         h_diff = abs(b% s_donor% center_h1 - b% s_donor% surface_h1)
         he_diff = abs(b% s_donor% center_he4 - b% s_donor% surface_he4)
         if (h_diff < b% CE_xa_diff_to_terminate &
            .and. he_diff < b% CE_xa_diff_to_terminate) then
            write(*,*) "Central and surface abundances below CE_xa_diff_to_terminate"
            write(*,*) "Terminating evolution"
            CE_binary_finish_step = terminate
            return
         end if

         !terminate if, for the current orbital separation, stripping the star down to the point
         !where CE_xa_diff_to_terminate would apply, the remaining layers would be Roche lobe
         !overflowing.
         if (b% CE_terminate_when_core_overflows) then
            do k = 1, b% s_donor% nz
               h_diff = abs(b% s_donor% center_h1 - b% s_donor% X(k))
               he_diff = abs(b% s_donor% center_he4 - b% s_donor% Y(k))
               if (h_diff < b% CE_xa_diff_to_terminate &
                  .and. he_diff < b% CE_xa_diff_to_terminate) then
                  rlobe = eval_rlobe(b% s_donor% m(k), b% m(b% a_i), b% separation)
                  if (b% s_donor% r(k) > rlobe) then
                     write(*,*) "Terminate due to CE_terminate_when_core_overflows"
                     write(*,*) "Terminating evolution"
                     CE_binary_finish_step = terminate
                     return
                  end if
                  exit
               end if
            end do
         end if

         if (b% period < b% CE_min_period_in_minutes*60d0) then
            write(*,*) "Orbital period is below CE_min_period_in_minutes"
            write(*,*) "Terminating evolution"
            CE_binary_finish_step = terminate
            return
         end if

      end function CE_binary_finish_step

      end module binary_ce

