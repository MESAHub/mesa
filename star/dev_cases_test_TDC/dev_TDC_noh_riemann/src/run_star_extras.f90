! ***********************************************************************
!
!   Copyright (C) 2010-2026  The MESA Team
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

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib

      implicit none

      real(dp), parameter :: gamma_law = 1.66667d0
      integer, parameter :: max_energy_balance_model = 1000
      integer :: energy_balance_io = -1

      include 'test_suite_extras_def.inc'

      contains

      include 'test_suite_extras.inc'


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% eos_rq% use_other_eos_component = .true.
         s% eos_rq% other_eos_frac => gamma_eos_frac
         s% eos_rq% other_eos_component => gamma_eos_component
         s% use_other_cgrav = .true.
         s% other_cgrav => zero_cgrav
         s% use_other_kap = .true.
         s% other_kap_get => constant_kap

         if (.not. s% x_logical_ctrl(1)) return

         call create_noh_test_model(id, s, ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'failed in create_noh_test_model')
      end subroutine extras_controls


      subroutine create_noh_test_model(id, s, ierr)
         use chem_def, only: ih1, ihe4
         use chem_lib, only: basic_composition_info
         use eos_lib, only: eos_gamma_DP_get_ET
         integer, intent(in) :: id
         type(star_info), pointer :: s
         integer, intent(out) :: ierr

         logical :: log_mesh
         integer :: k, nz, h1, he4
         real(dp) :: xh, xhe, z, abar, zbar, z53bar, z2bar, ye, &
            mass_correction, sumx, R_min, R_max, dr, rho_0, v_0, &
            P, r00, rp1, lnR00, lnR_min, lnR_max, dlnR

         include 'formats'

         ierr = 0
         s% u_flag = .true.
         s% RTI_flag = .false.

         log_mesh = s% split_merge_amr_log_zoning
         nz = s% split_merge_amr_nz_baseline
         s% nz = nz
         v_0 = s% x_ctrl(1)
         R_min = s% x_ctrl(3)
         R_max = s% x_ctrl(4)
         rho_0 = s% x_ctrl(5)
         P = s% x_ctrl(6)

         lnR_min = log(R_min)
         lnR_max = log(R_max)
         dr = (R_max - R_min)/(nz - 1)
         dlnR = (lnR_max - lnR_min)/(nz - 1)

         call star_set_net(id, 'basic_plus_fe56.net', ierr)
         if (ierr /= 0) return
         call star_set_var_info(id, ierr)
         if (ierr /= 0) return
         call star_set_chem_names(id, ierr)
         if (ierr /= 0) return
         call star_allocate_arrays(id, ierr)
         if (ierr /= 0) return

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         s% M_center = 0d0
         s% R_center = 0d0

         do k = nz, 1, -1
            if (k == nz) then
               r00 = R_min
               rp1 = s% R_center
               lnR00 = lnR_min
            else
               rp1 = r00
               if (log_mesh) then
                  lnR00 = lnR00 + dlnR
                  r00 = exp(lnR00)
               else
                  r00 = rp1 + dr
                  lnR00 = log(r00)
               end if
            end if

            s% r(k) = r00
            s% lnR(k) = lnR00
            s% u(k) = v_0
            s% alpha_RTI(k) = 0d0
            s% rho(k) = rho_0
            s% lnd(k) = log(rho_0)
            s% dm(k) = rho_0*(4d0*pi/3d0)*(pow3(r00) - pow3(rp1))
            if (k == nz) then
               s% m(k) = s% dm(k)
            else
               s% m(k) = s% m(k + 1) + s% dm(k)
            end if
            s% L(k) = 0d0
            s% Peos(k) = P
            s% lnPeos(k) = log(P)
            s% xa(1:s% species,k) = 0d0
            if (k == nz) then
               s% xa(he4,k) = 1d0
            else
               s% xa(h1,k) = 1d0
            end if

            call basic_composition_info( &
               s% species, s% chem_id, s% xa(:,k), xh, xhe, z, &
               abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
            call eos_gamma_DP_get_ET( &
               abar, rho_0, P, gamma_law, s% energy(k), s% T(k), ierr)
            if (ierr /= 0) return
            s% abar(k) = abar
            s% zbar(k) = zbar
            s% z53bar(k) = z53bar
            s% lnE(k) = log(s% energy(k))
            s% lnT(k) = log(s% T(k))
         end do

         s% L_center = 0d0
         s% star_mass = s% m(1)/Msun
         s% mstar = s% m(1)
         s% xmstar = s% m(1)
         s% q(1) = 1d0
         do k = 1, nz - 1
            s% dq(k) = s% dm(k)/s% xmstar
            s% q(k + 1) = s% q(k) - s% dq(k)
         end do
         s% dq(nz) = s% q(nz)

         s% model_number = 0
         s% star_age = 0d0
         s% initial_z = 0d0

         call star_write_model(id, 'mods/noh_start.mod', ierr)
         if (ierr /= 0) return
         write(*,2) 'mods/noh_start.mod', nz
         stop
      end subroutine create_noh_test_model


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type(star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         if (ierr /= 0) return

         if (restart) then
            open(newunit=energy_balance_io, file='noh_energy_balance.txt', &
               status='unknown', position='append', action='write', iostat=ierr)
         else
            open(newunit=energy_balance_io, file='noh_energy_balance.txt', &
               status='replace', action='write', iostat=ierr)
         end if
         if (ierr /= 0) return

         write(energy_balance_io,'(a)') &
            '# Cell energy balance for the accepted hydrodynamic solve.'
         write(energy_balance_io,'(a)') &
            '# dU = dm*(energy-energy_start); dK = dt*dm*dkedt; dW = dt*dm*dwork_dm.'
         write(energy_balance_io,'(a)') &
            '# balance = dU+dK+dW; equation_error is the unscaled energy residual in ergs.'
         write(energy_balance_io,'(a)') &
            '# model k age_s dt_s r_cm dm_g rho u P energy dU_erg dK_erg dW_erg balance_erg ' // &
            'equation_error_erg balance_minus_error_erg relative_balance'
      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         integer :: k, min_k
         logical :: okay
         real(dp) :: min_dr

         include 'formats'

         ierr = 0
         if (energy_balance_io /= -1) then
            close(energy_balance_io)
            energy_balance_io = -1
         end if
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         okay = .true.
         min_k = minloc(abs(s% r(1:s% nz) - 0.1d0), dim=1)
         min_dr = abs(s% r(min_k) - 0.1d0)
         write(*,2) 'nearest cell to r = 0.1', min_k, min_dr

         do k = min(s% nz, min_k + 10), s% nz
            if (s% r(k) < 0.09d0) then
               if (abs(s% u(k)) > 1d-3) then
                  okay = .false.
                  write(*,2) 'u behind shock', k, s% u(k), s% r(k)
               end if
               exit
            end if
         end do
         do k = max(1, min_k - 10), 1, -1
            if (s% r(k) > 0.11d0) then
               if (abs(s% u(k) + 1d0) > 1d-3) then
                  okay = .false.
                  write(*,2) 'u ahead of shock', k, s% u(k), s% r(k)
               end if
               exit
            end if
         end do

         if (okay) then
            write(*,*) 'Noh velocities are within tolerance'
         else
            write(*,*) 'Noh velocities are outside tolerance'
         end if

         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type(star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         extras_check_model = keep_going
         if (ierr /= 0) return
      end function extras_check_model


      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type(star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         extras_finish_step = keep_going
         if (ierr /= 0) return

         call write_energy_balance(s, ierr)
         if (ierr /= 0) then
            extras_finish_step = terminate
            s% termination_code = t_extras_finish_step
         end if
      end function extras_finish_step


      subroutine write_energy_balance(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: dU, dK, dW, balance, equation_error, scale, &
            total_dU, total_dK, total_dW, total_balance

         ierr = 0
         if (s% model_number > max_energy_balance_model) return
         if (energy_balance_io == -1) return

         total_dU = 0d0
         total_dK = 0d0
         total_dW = 0d0
         do k = 1, s% nz
            dU = s% dm(k)*(s% energy(k) - s% energy_start(k))
            dK = s% dt*s% dm(k)*s% dkedt(k)
            dW = s% dt*s% dm(k)*s% dwork_dm(k)
            total_dU = total_dU + dU
            total_dK = total_dK + dK
            total_dW = total_dW + dW
         end do
         total_balance = total_dU + total_dK + total_dW

         write(energy_balance_io,'(a,i0,a,i0,5(a,es24.16e3))',iostat=ierr) &
            '# model ', s% model_number, ' nz ', s% nz, &
            ' age_s ', s% star_age*secyer, ' dU_erg ', total_dU, &
            ' dK_erg ', total_dK, ' dW_erg ', total_dW, &
            ' balance_erg ', total_balance
         if (ierr /= 0) return

         do k = s% nz, 1, -1
            dU = s% dm(k)*(s% energy(k) - s% energy_start(k))
            dK = s% dt*s% dm(k)*s% dkedt(k)
            dW = s% dt*s% dm(k)*s% dwork_dm(k)
            balance = dU + dK + dW
            equation_error = s% ergs_error(k)
            scale = max(abs(dU) + abs(dK) + abs(dW), tiny(1d0))
            write(energy_balance_io,'(2(i8,1x),15(es24.16e3,1x))',iostat=ierr) &
               s% model_number, k, s% star_age*secyer, s% dt, s% r(k), s% dm(k), &
               s% rho(k), s% u(k), s% Peos(k), s% energy(k), dU, dK, dW, balance, &
               equation_error, balance - equation_error, balance/scale
            if (ierr /= 0) return
         end do
         flush(unit=energy_balance_io, iostat=ierr)
      end subroutine write_energy_balance


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character(len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character(len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine data_for_extra_profile_columns


      subroutine gamma_eos_frac( &
            handle, species, chem_id, net_iso, xa, Rho, logRho, T, logT, &
            frac, dfrac_dlogRho, dfrac_dlogT, ierr)
         integer, intent(in) :: handle, species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:), Rho, logRho, T, logT
         real(dp), intent(out) :: frac, dfrac_dlogRho, dfrac_dlogT
         integer, intent(out) :: ierr

         ierr = 0
         frac = 1d0
         dfrac_dlogRho = 0d0
         dfrac_dlogT = 0d0
      end subroutine gamma_eos_frac


      subroutine gamma_eos_component( &
            handle, species, chem_id, net_iso, xa, Rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)
         use chem_lib, only: basic_composition_info
         use eos_lib, only: eos_gamma_DT_get
         integer, intent(in) :: handle, species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:), Rho, logRho, T, logT
         real(dp), intent(inout) :: res(:), d_dlnd(:), d_dlnT(:), d_dxa(:,:)
         integer, intent(out) :: ierr
         real(dp) :: xh, xhe, z, abar, zbar, z2bar, z53bar, ye, &
            mass_correction, sumx, Pgas, Prad, energy, entropy

         call basic_composition_info( &
            species, chem_id, xa, xh, xhe, z, abar, zbar, z2bar, &
            z53bar, ye, mass_correction, sumx)
         call eos_gamma_DT_get( &
            handle, abar, Rho, logRho, T, logT, gamma_law, res, &
            d_dlnd, d_dlnT, Pgas, Prad, energy, entropy, ierr)
         d_dxa = 0d0
      end subroutine gamma_eos_component


      subroutine zero_cgrav(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% cgrav(1:s% nz) = 0d0
      end subroutine zero_cgrav


      subroutine constant_kap( &
            id, k, handle, species, chem_id, net_iso, xa, &
            log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, &
            d_lnfree_e_dlnT, eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, &
            dln_kap_dxa, ierr)
         use kap_def, only: num_kap_fracs
         integer, intent(in) :: id, k, handle, species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:), log10_rho, log10_T, lnfree_e, &
            d_lnfree_e_dlnRho, d_lnfree_e_dlnT, eta, d_eta_dlnRho, &
            d_eta_dlnT
         real(dp), intent(out) :: kap_fracs(num_kap_fracs), kap, &
            dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa(:)
         integer, intent(out) :: ierr

         ierr = 0
         kap_fracs = 0d0
         kap = 0.2d0
         dln_kap_dlnRho = 0d0
         dln_kap_dlnT = 0d0
         dln_kap_dxa = 0d0
      end subroutine constant_kap


      end module run_star_extras
