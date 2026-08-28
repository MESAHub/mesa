! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
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

      module star_lna_support

      use star_private_def
      use star_lna_turbulence_closures, only: &
         chi_coefficient_for_star_LNA, d_mlt_Pturb_face_for_star_LNA, &
         frozen_flux_luminosity_resid_for_star_LNA, Ptot_for_star_LNA, &
         pressure_components_for_star_LNA, &
         rsp2_forces_non_turbulent_cell, &
         rsp2_luminosity_resid_for_star_LNA, rsp2_luminosity_terms_for_star_LNA, &
         rsp2_terms_for_star_LNA_audit, &
         rsp2_turbulent_energy_inertia_for_star_LNA, rsp2_turbulent_energy_rhs_for_star_LNA, &
         star_LNA_eval_dlnPdm_qhse, star_LNA_HSE_grav_term, &
         star_LNA_L_conv_closure_ad => star_LNA_L_conv_ad, &
         tdc_face_state_for_star_LNA, &
         tdc_lna_active, tdc_luminosity_resid_for_star_LNA, &
         tdc_luminosity_terms_for_star_LNA, tdc_relation_for_star_LNA, &
         tdc_zero_w_for_star_LNA, static_face_pressure_for_star_LNA, &
         turbulent_energy_inertia_for_star_LNA, turbulent_viscous_heating_for_star_LNA, &
         Uq_cell_for_star_LNA, Uq_face_for_star_LNA
      use const_def, only: &
         clight, convective_mixing, crad, dp, ln10, Msun, pi, Rsun, secyer, sqrt_2_div_3
      use eos_lib, only: Radiation_Pressure
      use math_lib, only: pow2, pow4
      use auto_diff
      use auto_diff_support, only: &
         shift_m1, wrap, wrap_d_00, wrap_Hp_00, wrap_L_00, wrap_L_p1, &
         wrap_lnd_00, wrap_lnPeos_00, wrap_lnPeos_m1, wrap_Peos_00, &
         wrap_lnT_00, wrap_lnT_m1, wrap_r_00, wrap_r_p1, wrap_r_m1, &
         wrap_T_00, wrap_T_m1, wrap_v_00, wrap_v_p1, wrap_v_m1
      use hydro_vars, only: set_Teff_info_for_eqns
      use hydro_riemann, only: eval_Riemann_dudt_rhs
      use hydro_rsp2, only: Hp_face_for_rsp2_eqn
      use reconstructed_face_support, only: &
         get_effective_gradr_factor_ad, get_Lrad_per_gradT_face_ad, &
         get_reconstructed_face_eos_kap_ad
      use star_utils, only: get_Cp_face, get_dke_dt_dpe_dt, get_grada_face, &
         get_gradr_face, get_kap_face, get_Peos_face, get_rho_face, &
         get_scale_height_face, get_T_face
      use utils_lib, only: folder_exists, is_bad, mkdir

      implicit none

      private
      public :: star_LNA_problem
      public :: star_LNA_var_map
      public :: star_LNA_matrix
      public :: check_star_LNA_model
      public :: setup_star_LNA_problem
      public :: report_star_LNA_setup
      public :: assemble_density_rows
      public :: assemble_radius_rows
      public :: assemble_momentum_rows
      public :: assemble_energy_rows
      public :: assemble_luminosity_rows
      public :: assemble_rsp2_turbulent_rows
      public :: assemble_tdc_turbulent_rows
      public :: write_star_LNA_matrix_summary
      public :: solve_dense_star_LNA
      public :: free_star_LNA_problem
      public :: star_LNA_L_conv_ad

      integer, parameter :: lna_var_lnd = 1
      integer, parameter :: lna_var_lnR = 2
      integer, parameter :: lna_var_v = 3
      integer, parameter :: lna_var_u = 4
      integer, parameter :: lna_var_lnT = 5
      integer, parameter :: lna_var_L = 6
      integer, parameter :: lna_var_w = 7
      integer, parameter :: lna_var_Hp = 8

      integer, parameter :: lna_eq_density = 1
      integer, parameter :: lna_eq_radius = 2
      integer, parameter :: lna_eq_zero_velocity = 3
      integer, parameter :: lna_eq_momentum = 4
      integer, parameter :: lna_eq_surface_momentum = 5
      integer, parameter :: lna_eq_surface_pressure = 6
      integer, parameter :: lna_eq_surface_fixed_velocity = 7
      integer, parameter :: lna_eq_energy = 8
      integer, parameter :: lna_eq_luminosity_rsp_surface = 9
      integer, parameter :: lna_eq_luminosity_rsp2 = 10
      integer, parameter :: lna_eq_luminosity_tdc = 11
      integer, parameter :: lna_eq_luminosity_frozen_flux = 12
      integer, parameter :: lna_eq_surface_temperature = 13
      integer, parameter :: lna_eq_temperature_gradient = 14
      integer, parameter :: lna_eq_rsp2_turbulent_energy = 15
      integer, parameter :: lna_eq_rsp2_zero_w = 16
      integer, parameter :: lna_eq_rsp2_Hp = 17
      integer, parameter :: lna_eq_tdc_velocity = 18
      integer, parameter :: lna_eq_tdc_zero_w = 19
      integer, parameter :: lna_eq_mlt_static_temperature_gradient = 20
      integer, parameter :: num_lna_equations = lna_eq_mlt_static_temperature_gradient

      character(len=1), parameter :: star_LNA_backslash = achar(92)

      type star_LNA_var_map
         integer :: nz = 0
         integer :: nvar_per_zone = 0
         integer :: nvar_total = 0
         integer, allocatable :: var_id(:)
      end type star_LNA_var_map

      type star_LNA_matrix
         integer :: n = 0
         real(dp), allocatable :: A(:,:)
         real(dp), allocatable :: B(:,:)
      end type star_LNA_matrix

      type star_LNA_problem
         type(star_LNA_var_map) :: map
         type(star_LNA_matrix) :: mtx
         integer, allocatable :: eq_id(:)
      end type star_LNA_problem

      contains

      ! Model setup checks.
      subroutine check_star_LNA_model(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: largest_active_kick_mode

         ierr = 0

         call check_static_star_LNA_background(s, ierr)
         if (ierr /= 0) return

         ! The Riemann face state includes every active turbulent pressure.
         if (s% u_flag .and. .not. s% star_LNA_perturb_turbulent_pressure .and. &
               ((s% RSP2_flag .and. s% RSP2_alfap /= 0d0 .and. &
                  s% mixing_length_alpha /= 0d0) .or. &
               s% mlt_Pturb_factor > 0d0)) then
            write(*,'(a)') &
               'u_flag star_LNA requires turbulent pressure perturbations for this background.'
            ierr = -1
            return
         end if

         ! RSP2 places Uq in u_face rather than the cell momentum source.
         if (s% u_flag .and. s% RSP2_flag .and. &
               s% star_LNA_perturb_eddy_viscosity .and. &
               s% RSP2_alfam*s% mixing_length_alpha /= 0d0) then
            write(*,'(a)') &
               'star_LNA does not support RSP2 eddy viscosity with u_flag.'
            ierr = -1
            return
         end if

         if (s% rotation_flag) then
            write(*,'(a)') 'star_LNA does not support rotation perturbations.'
            ierr = -1
            return
         end if

         if (s% RTI_flag) then
            write(*,'(a)') 'star_LNA does not support RTI perturbations.'
            ierr = -1
            return
         end if

         if (s% use_mass_corrections) then
            write(*,'(a)') 'star_LNA does not support mass corrections.'
            ierr = -1
            return
         end if

         if (s% use_other_momentum .or. s% use_other_momentum_implicit) then
            write(*,'(a)') 'star_LNA does not support other_momentum hooks.'
            ierr = -1
            return
         end if

         if (s% use_other_pressure) then
            write(*,'(a)') 'star_LNA does not support other_pressure hooks.'
            ierr = -1
            return
         end if

         if (s% use_compression_outer_BC) then
            write(*,'(a)') 'star_LNA does not support the compression outer BC.'
            ierr = -1
            return
         end if

         if (s% use_other_surface_PT) then
            write(*,'(a)') 'star_LNA does not support other_surface_PT hooks.'
            ierr = -1
            return
         end if

         if (s% drag_coefficient > 0d0 .or. s% v_drag_factor > 0d0) then
            write(*,'(a)') 'star_LNA does not support velocity drag.'
            ierr = -1
            return
         end if

         if (s% eps_grav_form_for_energy_eqn) then
            write(*,'(a)') 'star_LNA requires the dE/dt energy equation form.'
            ierr = -1
            return
         end if

         if (s% constant_L) then
            write(*,'(a)') 'star_LNA does not support constant_L.'
            ierr = -1
            return
         end if

         if (s% star_LNA_num_modes < 1) then
            write(*,'(a)') 'star_LNA_num_modes must be >= 1.'
            ierr = -1
            return
         end if

         if (s% star_LNA_min_mode_frequency_uHz < 0d0) then
            write(*,'(a)') 'star_LNA_min_mode_frequency_uHz must be >= 0.'
            ierr = -1
            return
         end if

         if (s% star_LNA_max_abs_mode_eta < 0d0) then
            write(*,'(a)') 'star_LNA_max_abs_mode_eta must be >= 0.'
            ierr = -1
            return
         end if

         if (s% star_LNA_set_initial_velocity) then
            if (.not. s% v_flag) then
               write(*,'(a)') 'star_LNA_set_initial_velocity requires v_flag = .true.'
               write(*,'(a)') &
                  'The LNA solve can run without hydro v_flag, but the velocity kick must be a hydro variable.'
               ierr = -1
               return
            end if

            if (s% use_fixed_vsurf_outer_BC) then
               write(*,'(a)') &
                  'star_LNA_set_initial_velocity is incompatible with use_fixed_vsurf_outer_BC.'
               write(*,'(a)') &
                  'The fixed surface velocity BC constrains the LNA surface velocity, so the requested scaling is undefined.'
               ierr = -1
               return
            end if

            if (s% star_LNA_mode_for_period < 0 .or. &
                  s% star_LNA_mode_for_period >= s% star_LNA_num_modes) then
               write(*,'(a,i0,a,i0)') 'star_LNA_mode_for_period = ', &
                  s% star_LNA_mode_for_period, &
                  ' but star_LNA_num_modes = ', s% star_LNA_num_modes
               ierr = -1
               return
            end if

            if (s% star_LNA_kick_mode_1 < 1) then
               write(*,'(a,i0)') 'star_LNA_kick_mode_1 must be >= 1; got ', &
                  s% star_LNA_kick_mode_1
               ierr = -1
               return
            end if

            if (s% star_LNA_kick_mode_2 < 0) then
               write(*,'(a,i0)') 'star_LNA_kick_mode_2 must be >= 0; got ', &
                  s% star_LNA_kick_mode_2
               ierr = -1
               return
            end if

            if (s% star_LNA_kick_mode_3 < 0) then
               write(*,'(a,i0)') 'star_LNA_kick_mode_3 must be >= 0; got ', &
                  s% star_LNA_kick_mode_3
               ierr = -1
               return
            end if

            if (s% star_LNA_kick_fraction_1 < 0d0 .or. &
                  s% star_LNA_kick_fraction_2 < 0d0 .or. &
                  s% star_LNA_kick_fraction_3 < 0d0) then
               write(*,'(a)') 'star_LNA kick fractions must be nonnegative.'
               ierr = -1
               return
            end if

            largest_active_kick_mode = s% star_LNA_kick_mode_1
            if (s% star_LNA_kick_mode_2 > 0 .and. s% star_LNA_kick_fraction_2 > 0d0) &
               largest_active_kick_mode = max(largest_active_kick_mode, s% star_LNA_kick_mode_2)
            if (s% star_LNA_kick_mode_3 > 0 .and. s% star_LNA_kick_fraction_3 > 0d0) &
               largest_active_kick_mode = max(largest_active_kick_mode, s% star_LNA_kick_mode_3)
            if (largest_active_kick_mode > s% star_LNA_num_modes) then
               write(*,'(a,i0,a,i0)') 'largest star_LNA kick mode is ', &
                  largest_active_kick_mode, &
                  ' but star_LNA_num_modes = ', s% star_LNA_num_modes
               ierr = -1
               return
            end if

            if (s% star_LNA_kick_fraction_1 + &
                  merge(s% star_LNA_kick_fraction_2, 0d0, &
                     s% star_LNA_kick_mode_2 > 0) + &
                  merge(s% star_LNA_kick_fraction_3, 0d0, &
                     s% star_LNA_kick_mode_3 > 0) <= 0d0) then
               write(*,'(a)') 'at least one active star_LNA kick fraction must be positive.'
               ierr = -1
               return
            end if
         end if

         if (trim(s% star_LNA_solver) /= 'dense') then
            write(*,'(a)') 'star_LNA_solver must be ''dense''.'
            ierr = -1
            return
         end if

         select case (trim(s% star_LNA_convection_treatment))
         case ('perturbed', 'frozen', 'frozen_flux', 'mlt_static')
         case default
            write(*,'(a)') &
               'star_LNA_convection_treatment must be ''perturbed'', ' // &
               '''frozen'', ''frozen_flux'', or ''mlt_static''.'
            ierr = -1
            return
         end select

         if (trim(s% star_LNA_convection_treatment) == 'mlt_static' .and. &
               (s% RSP2_flag .or. s% MLT_option == 'TDC')) then
            write(*,'(a)') &
               'star_LNA mlt_static requires an MLT background without TDC or RSP2.'
            ierr = -1
            return
         end if

         if (s% RSP2_flag .and. frozen_flux_lna_selected(s)) then
            write(*,'(a)') 'star_LNA frozen_flux supports only models without RSP2.'
            ierr = -1
            return
         end if

         if (s% RSP2_flag .and. .not. s% star_LNA_include_rsp2) then
            write(*,'(a)') 'RSP2 models require star_LNA_include_rsp2 = .true.'
            ierr = -1
            return
         end if

         if (s% RSP2_flag .and. .not. s% star_LNA_perturb_turbulent_energy) then
            write(*,'(a)') &
               'RSP2 models require star_LNA_perturb_turbulent_energy = .true.'
            ierr = -1
            return
         end if

      end subroutine check_star_LNA_model


      subroutine check_static_star_LNA_background(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr

         ierr = 0

         if (s% mstar_dot /= 0d0) then
            write(*,'(a)') 'star_LNA requires mstar_dot = 0.'
            ierr = -1
            return
         end if
      end subroutine check_static_star_LNA_background

      ! Problem setup and dense matrix allocation.
      subroutine setup_star_LNA_problem(s, problem, ierr)
         type(star_info), pointer :: s
         type(star_LNA_problem), intent(out) :: problem
         integer, intent(out) :: ierr

         ierr = 0
         call setup_star_LNA_var_map(s, problem% map, ierr)
         if (ierr /= 0) return

         call setup_star_LNA_equation_map(s, problem, ierr)
         if (ierr /= 0) then
            call free_star_LNA_problem(problem)
            return
         end if

         call allocate_star_LNA_matrix(problem% map, problem% mtx, ierr)
         if (ierr /= 0) call free_star_LNA_problem(problem)
      end subroutine setup_star_LNA_problem


      subroutine setup_star_LNA_var_map(s, map, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(out) :: map
         integer, intent(out) :: ierr
         integer :: iv

         ierr = 0
         map%nz = s% nz
         map%nvar_per_zone = 5
         if (s% RSP2_flag) then
            map%nvar_per_zone = map%nvar_per_zone + 2
         else if (tdc_lna_active(s)) then
            map%nvar_per_zone = map%nvar_per_zone + 1
         end if
         map%nvar_total = map%nz*map%nvar_per_zone

         allocate(map%var_id(map%nvar_per_zone), stat=ierr)
         if (ierr /= 0) return

         iv = 0
         call add_var(lna_var_lnd)
         call add_var(lna_var_lnR)
         if (s% u_flag) then
            call add_var(lna_var_u)
         else
            call add_var(lna_var_v)
         end if
         call add_var(lna_var_lnT)
         call add_var(lna_var_L)
         if (s% RSP2_flag) then
            call add_var(lna_var_w)
            call add_var(lna_var_Hp)
         else if (tdc_lna_active(s)) then
            call add_var(lna_var_w)
         end if

         contains

         subroutine add_var(var_id)
            integer, intent(in) :: var_id
            iv = iv + 1
            map%var_id(iv) = var_id
         end subroutine add_var
      end subroutine setup_star_LNA_var_map


      subroutine setup_star_LNA_equation_map(s, problem, ierr)
         type(star_info), pointer :: s
         type(star_LNA_problem), intent(inout) :: problem
         integer, intent(out) :: ierr
         integer :: eq_id, iv, k, row, var_id

         ierr = 0
         allocate(problem% eq_id(problem% map% nvar_total), stat=ierr)
         if (ierr /= 0) return

         do k = 1, problem% map% nz
            do iv = 1, problem% map% nvar_per_zone
               var_id = problem% map% var_id(iv)
               row = matrix_index(problem% map, k, var_id)
               eq_id = equation_id_for_star_LNA(s, k, var_id)
               if (row <= 0 .or. eq_id <= 0) then
                  write(*,'(a,2(i0,1x))') &
                     'star_LNA failed to identify equation for k and variable: ', k, var_id
                  ierr = -1
                  return
               end if
               problem% eq_id(row) = eq_id
            end do
         end do
      end subroutine setup_star_LNA_equation_map


      integer function equation_id_for_star_LNA(s, k, var_id) result(eq_id)
         type(star_info), pointer :: s
         integer, intent(in) :: k, var_id

         eq_id = 0
         select case (var_id)
         case (lna_var_lnd)
            eq_id = lna_eq_density
         case (lna_var_lnR)
            if (force_zero_velocity_for_star_LNA(s, k)) then
               eq_id = lna_eq_zero_velocity
            else
               eq_id = lna_eq_radius
            end if
         case (lna_var_v, lna_var_u)
            if (k > 1) then
               eq_id = lna_eq_momentum
            else if (s% use_fixed_vsurf_outer_BC) then
               eq_id = lna_eq_surface_fixed_velocity
            else if (use_surface_momentum_row_for_star_LNA(s)) then
               eq_id = lna_eq_surface_momentum
            else
               eq_id = lna_eq_surface_pressure
            end if
         case (lna_var_lnT)
            eq_id = lna_eq_energy
         case (lna_var_L)
            if (use_rsp_lsurf_row_for_star_LNA(s, k)) then
               eq_id = lna_eq_luminosity_rsp_surface
            else if (s% RSP2_flag) then
               eq_id = lna_eq_luminosity_rsp2
            else if (tdc_lna_active(s) .and. k > 1) then
               eq_id = lna_eq_luminosity_tdc
            else if (frozen_flux_lna_active(s) .and. k > 1) then
               eq_id = lna_eq_luminosity_frozen_flux
            else if (k == 1) then
               eq_id = lna_eq_surface_temperature
            else if (trim(s% star_LNA_convection_treatment) == 'mlt_static') then
               eq_id = lna_eq_mlt_static_temperature_gradient
            else
               eq_id = lna_eq_temperature_gradient
            end if
         case (lna_var_w)
            if (s% RSP2_flag) then
               if (rsp2_forces_non_turbulent_cell(s, k)) then
                  eq_id = lna_eq_rsp2_zero_w
               else
                  eq_id = lna_eq_rsp2_turbulent_energy
               end if
            else if (tdc_zero_w_for_star_LNA(s, k)) then
               eq_id = lna_eq_tdc_zero_w
            else
               eq_id = lna_eq_tdc_velocity
            end if
         case (lna_var_Hp)
            eq_id = lna_eq_rsp2_Hp
         end select
      end function equation_id_for_star_LNA


      subroutine allocate_star_LNA_matrix(map, mtx, ierr)
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(out) :: mtx
         integer, intent(out) :: ierr

         ierr = 0
         mtx%n = map%nvar_total
         allocate(mtx%A(mtx%n, mtx%n), mtx%B(mtx%n, mtx%n), stat=ierr)
         if (ierr /= 0) then
            write(*,'(a,i0,a)') &
               'star_LNA failed to allocate dense A/B matrices with dimension ', &
               mtx%n, '.'
            return
         end if
         mtx%A = 0d0
         mtx%B = 0d0
      end subroutine allocate_star_LNA_matrix

      ! Matrix row assembly.
      !
      ! The generalized problem is
      !
      !   A*x = sigma*B*x
      !
      ! Algebraic rows write only A. Dynamic rows put the static RHS
      ! perturbation in A and the time derivative inertia in B.

      ! Equation:
      !   rho_k*DeltaV_k = dm_k
      !
      ! Linearized form:
      !   delta ln(rho_k) + delta ln(DeltaV_k) = 0
      !
      ! Matrix:
      !   A(row,:) receives the AD partials of ln(rho_k) + ln(DeltaV_k).
      !   B(row,:) is unchanged because this is an algebraic closure.
      subroutine assemble_density_rows(s, map, mtx, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(out) :: ierr
         integer :: k, row
         type(auto_diff_real_star_order1) :: cell_volume_ad, density_resid_ad

         ierr = 0
         do k = 1, map%nz
            row = matrix_index(map, k, lna_var_lnd)
            if (row <= 0) then
               write(*,'(a,i0)') 'star_LNA density row is missing a variable at k = ', k
               ierr = -1
               return
            end if

            call cell_volume_for_star_LNA(s, k, cell_volume_ad, ierr)
            if (ierr /= 0) return

            density_resid_ad = wrap_lnd_00(s, k) + log(cell_volume_ad)
            call add_ad_partials_to_A(map, mtx, row, k, 1d0, density_resid_ad, ierr)
            if (ierr /= 0) return
         end do
      end subroutine assemble_density_rows


      subroutine cell_volume_for_star_LNA(s, k, cell_volume_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: cell_volume_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: r_outer_ad, r_inner_ad

         ierr = 0
         r_outer_ad = wrap_r_00(s, k)
         r_inner_ad = wrap_r_p1(s, k)
         cell_volume_ad = (4d0*pi/3d0)* &
            (r_outer_ad*r_outer_ad*r_outer_ad - r_inner_ad*r_inner_ad*r_inner_ad)

         if (cell_volume_ad%val <= 0d0) then
            write(*,'(a,i0,1pe14.6)') &
               'star_LNA found nonpositive cell volume at k = ', k, cell_volume_ad%val
            ierr = -1
            return
         end if
      end subroutine cell_volume_for_star_LNA


      ! Equation:
      !   d lnR_k/dt = v_face,k/r_k
      !
      ! Linearized form:
      !   delta(v_face,k/r_k) = sigma*delta lnR_k
      !
      ! Matrix:
      !   A(row,:) receives the AD partials of v_face,k/r_k.
      !   B(row,lnR_k) = 1.
      !   Cells with constrained velocity use the algebraic row delta velocity_k = 0.
      subroutine assemble_radius_rows(s, map, mtx, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(out) :: ierr
         integer :: k, row, col_lnR, col_velocity, velocity_var
         type(auto_diff_real_star_order1) :: dlnR_dt

         ierr = 0
         velocity_var = velocity_var_for_star_LNA(map)
         do k = 1, map%nz
            row = matrix_index(map, k, lna_var_lnR)
            col_lnR = matrix_index(map, k, lna_var_lnR)
            col_velocity = matrix_index(map, k, velocity_var)
            if (row <= 0 .or. col_lnR <= 0 .or. col_velocity <= 0) then
               write(*,'(a,i0)') 'star_LNA radius row is missing a variable at k = ', k
               ierr = -1
               return
            end if

            if (force_zero_velocity_for_star_LNA(s, k)) then
               mtx%A(row, col_velocity) = 1d0
            else
               if (s% r(k) <= 0d0) then
                  write(*,'(a,i0)') 'star_LNA found nonpositive radius at k = ', k
                  ierr = -1
                  return
               end if

               dlnR_dt = wrap_v_00(s, k)/s% r(k)
               call add_ad_partials_to_A(map, mtx, row, k, 1d0, dlnR_dt, ierr)
               if (ierr /= 0) return
               mtx%B(row, col_lnR) = 1d0
            end if
         end do
      end subroutine assemble_radius_rows


      ! Face velocity equation:
      !   dv_k/dt = grav_k + Uq_k
      !             - (DeltaPtot_k + DeltaPmlt_turb_k)/(dm_face_k/A_k)
      ! With u_flag, use the cell-centered Riemann momentum equation instead.
      !
      ! Linearized form:
      !   delta(RHS_k) = sigma*delta velocity_k
      !
      ! Matrix:
      !   A(row,:) receives the AD partials of RHS_k.
      !   B(row,velocity_k) = 1.
      !   The surface row may instead impose a pressure or fixed velocity BC.
      subroutine assemble_momentum_rows(s, map, mtx, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(out) :: ierr
         integer :: k, row, col_velocity, velocity_var
         type(auto_diff_real_star_order1) :: velocity_rhs_ad

         ierr = 0
         velocity_var = velocity_var_for_star_LNA(map)
         do k = 1, map%nz
            row = matrix_index(map, k, velocity_var)
            col_velocity = matrix_index(map, k, velocity_var)
            if (row <= 0 .or. col_velocity <= 0) then
               write(*,'(a,i0)') 'star_LNA momentum row is missing a variable at k = ', k
               ierr = -1
               return
            end if

            if (k == 1) then
               call assemble_surface_velocity_row(s, map, mtx, row, col_velocity, ierr)
               if (ierr /= 0) return
               cycle
            end if

            call momentum_rhs_for_star_LNA(s, k, velocity_rhs_ad, ierr)
            if (ierr /= 0) return
            call add_ad_partials_to_A(map, mtx, row, k, 1d0, velocity_rhs_ad, ierr)
            if (ierr /= 0) return
            mtx%B(row, col_velocity) = 1d0
         end do
      end subroutine assemble_momentum_rows


      ! Surface momentum form mirrors the interior equation with
      !   DeltaPtot_1 = P_bc - Ptot_1.
      ! If the active outer BC is not a momentum row, star_LNA imposes the
      ! algebraic pressure boundary delta(lnP_bc - lnPeos_1) = 0.
      subroutine assemble_surface_velocity_row(s, map, mtx, row, col_velocity, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(in) :: row, col_velocity
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            P_bc_ad, lnP_bc_ad, resid_ad, velocity_rhs_ad

         ierr = 0

         if (s% use_fixed_vsurf_outer_BC) then
            mtx%A(row, col_velocity) = 1d0
            return
         end if

         call surface_P_bc_for_star_LNA(s, P_bc_ad, lnP_bc_ad, ierr)
         if (ierr /= 0) return

         if (use_surface_momentum_row_for_star_LNA(s)) then
            call surface_momentum_rhs_for_star_LNA(s, P_bc_ad, velocity_rhs_ad, ierr)
            if (ierr /= 0) return
            call add_ad_partials_to_A(map, mtx, row, 1, 1d0, velocity_rhs_ad, ierr)
            if (ierr /= 0) return
            mtx%B(row, col_velocity) = 1d0
         else
            resid_ad = lnP_bc_ad - wrap_lnPeos_00(s, 1)
            call add_ad_partials_to_A(map, mtx, row, 1, 1d0, resid_ad, ierr)
            if (ierr /= 0) return
         end if
      end subroutine assemble_surface_velocity_row


      logical function use_surface_momentum_row_for_star_LNA(s) result(use_momentum_row)
         type(star_info), pointer :: s

         use_momentum_row = s% use_momentum_outer_BC .or. &
            s% use_zero_Pgas_outer_BC .or. s% use_fixed_Psurf_outer_BC
      end function use_surface_momentum_row_for_star_LNA


      ! Interior momentum RHS in acceleration units. The face-velocity path
      ! follows hydro_momentum.get1_momentum_eqn; u_flag reuses
      ! hydro_riemann.eval_Riemann_dudt_rhs.
      subroutine momentum_rhs_for_star_LNA(s, k, velocity_rhs_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: velocity_rhs_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: grav_ad, area_ad, dm_div_A_ad, &
            dPtot_ad, d_mlt_Pturb_ad, Uq_ad, P_surf_ad

         ierr = 0

         if (s% u_flag) then
            P_surf_ad = 0d0
            ! Omit finite-step centering and add the static LNA Uq below.
            call eval_Riemann_dudt_rhs( &
               s, k, P_surf_ad, use_time_centering=.false., &
               include_tdc_Uq=.false., dudt_expected_ad=velocity_rhs_ad, ierr=ierr)
            if (ierr /= 0) return
            call Uq_cell_for_star_LNA(s, k, Uq_ad, ierr)
            if (ierr /= 0) return
            velocity_rhs_ad = velocity_rhs_ad + Uq_ad
            return
         end if

         call star_LNA_HSE_grav_term(s, k, grav_ad, area_ad, ierr)
         if (ierr /= 0) return
         dm_div_A_ad = star_LNA_dm_face(s, k)/area_ad

         call dPtot_face_for_star_LNA(s, k, dPtot_ad, ierr)
         if (ierr /= 0) return
         call d_mlt_Pturb_face_for_star_LNA(s, k, d_mlt_Pturb_ad)

         call Uq_face_for_star_LNA(s, k, Uq_ad, ierr)
         if (ierr /= 0) return

         velocity_rhs_ad = grav_ad + Uq_ad - &
            (dPtot_ad + d_mlt_Pturb_ad)/dm_div_A_ad
      end subroutine momentum_rhs_for_star_LNA


      subroutine surface_momentum_rhs_for_star_LNA(s, P_bc_ad, velocity_rhs_ad, ierr)
         type(star_info), pointer :: s
         type(auto_diff_real_star_order1), intent(in) :: P_bc_ad
         type(auto_diff_real_star_order1), intent(out) :: velocity_rhs_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: grav_ad, area_ad, dm_div_A_ad, &
            Ptot00_ad, dPtot_ad, d_mlt_Pturb_ad, Uq_ad

         ierr = 0
         if (s% u_flag) then
            call eval_Riemann_dudt_rhs( &
               s, 1, P_bc_ad, use_time_centering=.false., &
               include_tdc_Uq=.false., dudt_expected_ad=velocity_rhs_ad, ierr=ierr)
            if (ierr /= 0) return
            call Uq_cell_for_star_LNA(s, 1, Uq_ad, ierr)
            if (ierr /= 0) return
            velocity_rhs_ad = velocity_rhs_ad + Uq_ad
            return
         end if

         call star_LNA_HSE_grav_term(s, 1, grav_ad, area_ad, ierr)
         if (ierr /= 0) return
         dm_div_A_ad = star_LNA_dm_face(s, 1)/area_ad

         call Ptot_for_star_LNA(s, 1, Ptot00_ad, ierr)
         if (ierr /= 0) return
         dPtot_ad = P_bc_ad - Ptot00_ad
         call d_mlt_Pturb_face_for_star_LNA(s, 1, d_mlt_Pturb_ad)

         call Uq_face_for_star_LNA(s, 1, Uq_ad, ierr)
         if (ierr /= 0) return

         velocity_rhs_ad = grav_ad + Uq_ad - &
            (dPtot_ad + d_mlt_Pturb_ad)/dm_div_A_ad
      end subroutine surface_momentum_rhs_for_star_LNA


      subroutine dPtot_face_for_star_LNA(s, k, dPtot_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dPtot_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Ptot00_ad, Ptotm1_ad

         ierr = 0
         call Ptot_for_star_LNA(s, k, Ptot00_ad, ierr)
         if (ierr /= 0) return
         call Ptot_for_star_LNA(s, k - 1, Ptotm1_ad, ierr)
         if (ierr /= 0) return
         Ptotm1_ad = shift_m1(Ptotm1_ad)
         dPtot_ad = Ptotm1_ad - Ptot00_ad
      end subroutine dPtot_face_for_star_LNA


      real(dp) function star_LNA_dm_face(s, k) result(dm_face)
         type(star_info), pointer :: s
         integer, intent(in) :: k

         if (s% use_mass_corrections) then
            if (k > 1) then
               dm_face = 0.5d0*( &
                  s% dm(k)*s% mass_correction(k) + &
                  s% dm(k - 1)*s% mass_correction(k - 1))
            else
               dm_face = 0.5d0*s% dm(k)*s% mass_correction(k)
            end if
         else if (k > 1) then
            dm_face = 0.5d0*(s% dm(k) + s% dm(k - 1))
         else
            dm_face = 0.5d0*s% dm(k)
         end if
      end function star_LNA_dm_face


      subroutine surface_P_bc_for_star_LNA(s, P_bc_ad, lnP_bc_ad, ierr)
         type(star_info), pointer :: s
         type(auto_diff_real_star_order1), intent(out) :: P_bc_ad, lnP_bc_ad
         integer, intent(out) :: ierr

         if (s% use_zero_Pgas_outer_BC) then
            call constant_surface_P_bc_for_star_LNA(Radiation_Pressure(s% T_start(1)), P_bc_ad, lnP_bc_ad)
            ierr = 0
            return
         end if

         if (s% use_fixed_Psurf_outer_BC) then
            call constant_surface_P_bc_for_star_LNA(s% fixed_Psurf, P_bc_ad, lnP_bc_ad)
            ierr = 0
            return
         end if

         call atmosphere_surface_P_bc_for_star_LNA(s, P_bc_ad, lnP_bc_ad, ierr)
      end subroutine surface_P_bc_for_star_LNA


      subroutine constant_surface_P_bc_for_star_LNA(P_bc, P_bc_ad, lnP_bc_ad)
         real(dp), intent(in) :: P_bc
         type(auto_diff_real_star_order1), intent(out) :: P_bc_ad, lnP_bc_ad

         P_bc_ad = P_bc
         lnP_bc_ad = log(P_bc)
      end subroutine constant_surface_P_bc_for_star_LNA


      subroutine atmosphere_surface_P_bc_for_star_LNA(s, P_bc_ad, lnP_bc_ad, ierr)
         type(star_info), pointer :: s
         type(auto_diff_real_star_order1), intent(out) :: P_bc_ad, lnP_bc_ad
         integer, intent(out) :: ierr
         real(dp) :: r, L, Teff, lnT_surf, lnP_surf, P_surf, P_bc, &
            dlnTsurf_dL, dlnTsurf_dlnR, dlnTsurf_dlnM, dlnTsurf_dlnkap, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnM, dlnPsurf_dlnkap, &
            dlnP_bc_dlnd, dlnP_bc_dlnT, dlnP_bc_dL, dlnP_bc_dlnR, &
            dlnkap_dlnd, dlnkap_dlnT, dP0, dP0_dlnR, &
            dlnP_bc_dlnPsurf, dlnP_bc_dP0
         logical :: offset_P_to_cell_center
         logical, parameter :: skip_partials = .false.
         logical, parameter :: need_P_surf = .true.
         logical, parameter :: need_T_surf = .false.

         ierr = 0
         call set_Teff_info_for_eqns(s, skip_partials, &
            need_P_surf, need_T_surf, r, L, Teff, &
            lnT_surf, dlnTsurf_dL, dlnTsurf_dlnR, dlnTsurf_dlnM, dlnTsurf_dlnkap, &
            lnP_surf, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnM, dlnPsurf_dlnkap, ierr)
         if (ierr /= 0) return

         P_surf = exp(lnP_surf)
         offset_P_to_cell_center = .not. s% use_momentum_outer_BC
         if (offset_P_to_cell_center) then
            dP0 = s% cgrav(1)*s% m_grav(1)*s% dm(1)/(8d0*pi*pow4(r))
            dP0_dlnR = -4d0*dP0
         else
            dP0 = 0d0
            dP0_dlnR = 0d0
         end if
         P_bc = P_surf + dP0

         dlnP_bc_dlnPsurf = P_surf/P_bc
         dlnP_bc_dP0 = 1d0/P_bc

         dlnkap_dlnd = s% d_opacity_dlnd(1)/s% opacity(1)
         dlnkap_dlnT = s% d_opacity_dlnT(1)/s% opacity(1)

         dlnP_bc_dlnd = dlnP_bc_dlnPsurf*dlnPsurf_dlnkap*dlnkap_dlnd
         dlnP_bc_dlnT = dlnP_bc_dlnPsurf*dlnPsurf_dlnkap*dlnkap_dlnT
         dlnP_bc_dL = dlnP_bc_dlnPsurf*dlnPsurf_dL
         dlnP_bc_dlnR = dlnP_bc_dlnPsurf*dlnPsurf_dlnR + dlnP_bc_dP0*dP0_dlnR

         call wrap(P_bc_ad, P_bc, &
            0d0, P_bc*dlnP_bc_dlnd, 0d0, &
            0d0, P_bc*dlnP_bc_dlnT, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, P_bc*dlnP_bc_dlnR, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, P_bc*dlnP_bc_dL, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0)

         call wrap(lnP_bc_ad, log(P_bc), &
            0d0, dlnP_bc_dlnd, 0d0, &
            0d0, dlnP_bc_dlnT, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, dlnP_bc_dlnR, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, dlnP_bc_dL, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0)
      end subroutine atmosphere_surface_P_bc_for_star_LNA


      ! Equation:
      !   d e_eff,k/dt = -dL_k/dm + sources_k - dwork_k/dm
      !
      ! Linearized form:
      !   delta[-dL/dm + sources - dwork/dm] = sigma*delta e_eff
      !
      ! For the active MESA dedt total energy form before the P*d(1/rho)
      ! time centering branch, e_eff includes internal, turbulent, kinetic, and
      ! potential specific energy.  Around a static background the kinetic
      ! derivative vanishes, but keeping the MESA helper here preserves the
      ! correct linearization if the background has small nonzero velocities.
      !
      ! Matrix:
      !   A(row,:) receives the RHS AD partials.
      !   B(row,:) receives the AD partials of e_eff.
      subroutine assemble_energy_rows(s, map, mtx, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(out) :: ierr
         integer :: k, row
         type(auto_diff_real_star_order1) :: energy_rhs_ad, energy_inertia_ad

         ierr = 0
         do k = 1, map%nz
            row = matrix_index(map, k, lna_var_lnT)
            if (row <= 0) then
               write(*,'(a,i0)') 'star_LNA energy row is missing a variable at k = ', k
               ierr = -1
               return
            end if

            call energy_rhs_for_star_LNA(s, k, energy_rhs_ad, ierr)
            if (ierr /= 0) return
            call add_ad_partials_to_A(map, mtx, row, k, 1d0, energy_rhs_ad, ierr)
            if (ierr /= 0) return

            call energy_inertia_for_star_LNA(s, k, energy_inertia_ad, ierr)
            if (ierr /= 0) return
            call add_ad_partials_to_B(map, mtx, row, k, 1d0, energy_inertia_ad, ierr)
            if (ierr /= 0) return
         end do
      end subroutine assemble_energy_rows


      subroutine energy_rhs_for_star_LNA(s, k, energy_rhs_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: energy_rhs_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: dL_dm_ad, sources_ad, dwork_dm_ad

         ierr = 0
         call dL_dm_for_star_LNA(s, k, dL_dm_ad)
         call energy_sources_for_star_LNA(s, k, sources_ad, ierr)
         if (ierr /= 0) return
         call dwork_dm_for_star_LNA(s, k, dwork_dm_ad, ierr)
         if (ierr /= 0) return

         energy_rhs_ad = -dL_dm_ad + sources_ad - dwork_dm_ad
      end subroutine energy_rhs_for_star_LNA


      subroutine dL_dm_for_star_LNA(s, k, dL_dm_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dL_dm_ad

         dL_dm_ad = (wrap_L_00(s, k) - wrap_L_p1(s, k))/s% dm(k)
      end subroutine dL_dm_for_star_LNA


      subroutine energy_sources_for_star_LNA(s, k, sources_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: sources_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: eps_nuc_ad, non_nuc_neu_ad, &
            extra_heat_ad, Eq_ad

         ierr = 0

         if (s% eps_nuc_factor == 0d0 .or. s% nonlocal_NiCo_decay_heat) then
            eps_nuc_ad = 0d0
         else if (s% op_split_burn .and. s% T_start(k) >= s% op_split_burn_min_T) then
            eps_nuc_ad = s% burn_avg_epsnuc(k)
         else
            eps_nuc_ad = 0d0
            eps_nuc_ad%val = s% eps_nuc(k)
            eps_nuc_ad%d1Array(i_lnd_00) = s% d_epsnuc_dlnd(k)
            eps_nuc_ad%d1Array(i_lnT_00) = s% d_epsnuc_dlnT(k)
         end if

         non_nuc_neu_ad = 0d0
         non_nuc_neu_ad%val = 0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))
         non_nuc_neu_ad%d1Array(i_lnd_00) = 0.5d0*s% d_nonnucneu_dlnd(k)
         non_nuc_neu_ad%d1Array(i_lnT_00) = 0.5d0*s% d_nonnucneu_dlnT(k)

         extra_heat_ad = s% extra_heat(k) + s% irradiation_heat(k)

         call turbulent_viscous_heating_for_star_LNA(s, k, Eq_ad, ierr)
         if (ierr /= 0) return

         sources_ad = eps_nuc_ad - non_nuc_neu_ad + extra_heat_ad + Eq_ad
      end subroutine energy_sources_for_star_LNA


      ! Static pressure work row:
      !   dwork/dm = [P*A*v]_k - [P*A*v]_{k+1}
      !
      ! For a static LNA background, v0 = 0, so the first order row is
      !   delta(dwork/dm) = [P0*delta(A*v)]_k - [P0*delta(A*v)]_{k+1}.
      ! This intentionally does not include delta(P)*A*v0 terms.
      subroutine dwork_dm_for_star_LNA(s, k, dwork_dm_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dwork_dm_ad
         integer, intent(out) :: ierr
         real(dp) :: P_out, P_in

         ierr = 0
         if (s% u_flag) then
            P_out = s% P_face_ad(k)%val
         else
            call static_face_pressure_for_star_LNA(s, k, P_out, ierr)
            if (ierr /= 0) return
         end if
         if (k < s% nz) then
            if (s% u_flag) then
               P_in = s% P_face_ad(k + 1)%val
            else
               call static_face_pressure_for_star_LNA(s, k + 1, P_in, ierr)
               if (ierr /= 0) return
            end if
            dwork_dm_ad = 4d0*pi*( &
               P_out*s% r(k)*s% r(k)*wrap_v_00(s, k) - &
               P_in*s% r(k + 1)*s% r(k + 1)*wrap_v_p1(s, k))/s% dm(k)
         else
            dwork_dm_ad = &
               4d0*pi*P_out*s% r(k)*s% r(k)*wrap_v_00(s, k)/s% dm(k)
         end if
      end subroutine dwork_dm_for_star_LNA


      subroutine energy_inertia_for_star_LNA(s, k, energy_inertia_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: energy_inertia_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: turbulent_inertia_ad, &
            mechanical_inertia_ad

         ierr = 0
         energy_inertia_ad = 0d0
         energy_inertia_ad%val = s% energy(k)
         energy_inertia_ad%d1Array(i_lnd_00) = &
            s% dE_dRho_for_partials(k)*s% rho(k)
         energy_inertia_ad%d1Array(i_lnT_00) = s% Cv_for_partials(k)*s% T(k)

         call turbulent_energy_inertia_for_star_LNA(s, k, turbulent_inertia_ad, ierr)
         if (ierr /= 0) return
         energy_inertia_ad = energy_inertia_ad + turbulent_inertia_ad

         call mechanical_energy_inertia_for_star_LNA(s, k, mechanical_inertia_ad, ierr)
         if (ierr /= 0) return
         energy_inertia_ad = energy_inertia_ad + mechanical_inertia_ad
      end subroutine energy_inertia_for_star_LNA


      subroutine mechanical_energy_inertia_for_star_LNA(s, k, inertia_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: inertia_ad
         integer, intent(out) :: ierr
         real(dp) :: dke_dt, d_dkedt_dv00, d_dkedt_dvp1, dpe_dt, &
            d_dpedt_dlnR00, d_dpedt_dlnRp1

         ierr = 0
         inertia_ad = 0d0
         if (trim(s% energy_eqn_option) /= 'dedt') return
         if (s% using_velocity_time_centering .and. &
               s% use_P_d_1_div_rho_form_of_work_when_time_centering_velocity) return

         call get_dke_dt_dpe_dt(s, k, 1d0, &
            dke_dt, d_dkedt_dv00, d_dkedt_dvp1, &
            dpe_dt, d_dpedt_dlnR00, d_dpedt_dlnRp1, ierr)
         if (ierr /= 0) return

         call wrap(inertia_ad, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, d_dpedt_dlnR00, d_dpedt_dlnRp1, &
            0d0, d_dkedt_dv00, d_dkedt_dvp1, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0)
      end subroutine mechanical_energy_inertia_for_star_LNA


      ! Luminosity/temperature slot. The variable in this slot is L, but the
      ! equation depends on the active physics:
      !   RSP surface:       L_1 - Lsurf(T_1,r_1) = 0
      !   RSP2:              Lr + Lc + Lt - L = 0
      !   TDC:               Lrad + Lconv - L = 0
      !   Frozen flux:       Lrad + Lconv0 - L = 0
      !   Surface without RSP: surface temperature boundary
      !   Static without TDC:  temperature gradient relation
      subroutine assemble_luminosity_rows(s, map, mtx, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(out) :: ierr
         integer :: k, row
         type(auto_diff_real_star_order1) :: luminosity_resid_ad

         ierr = 0
         do k = 1, map%nz
            row = matrix_index(map, k, lna_var_L)
            if (row <= 0) then
               write(*,'(a,i0)') 'star_LNA luminosity row is missing a variable at k = ', k
               ierr = -1
               return
            end if

            if (use_rsp_lsurf_row_for_star_LNA(s, k)) then
               call rsp_lsurf_resid_for_star_LNA(s, luminosity_resid_ad)
            else if (s% RSP2_flag) then
               call rsp2_luminosity_resid_for_star_LNA(s, k, luminosity_resid_ad)
            else if (tdc_lna_active(s) .and. k > 1) then
               call tdc_luminosity_resid_for_star_LNA(s, k, luminosity_resid_ad, ierr)
               if (ierr /= 0) return
            else if (frozen_flux_lna_active(s) .and. k > 1) then
               call frozen_flux_luminosity_resid_for_star_LNA(s, k, luminosity_resid_ad, ierr)
               if (ierr /= 0) return
            else if (k == 1) then
               call surface_temperature_resid_for_star_LNA(s, luminosity_resid_ad, ierr)
               if (ierr /= 0) return
            else
               call temperature_gradient_resid_for_star_LNA(s, k, luminosity_resid_ad, ierr)
               if (ierr /= 0) return
            end if

            call add_ad_partials_to_A(map, mtx, row, k, 1d0, luminosity_resid_ad, ierr)
            if (ierr /= 0) return
         end do
      end subroutine assemble_luminosity_rows


      logical function use_rsp_lsurf_row_for_star_LNA(s, k) result(use_rsp_lsurf_row)
         type(star_info), pointer :: s
         integer, intent(in) :: k

         use_rsp_lsurf_row = k == 1 .and. s% use_RSP_L_eqn_outer_BC
         if (use_rsp_lsurf_row .and. s% RSP2_flag .and. s% RSP2_use_L_eqn_at_surface) &
            use_rsp_lsurf_row = .false.
      end function use_rsp_lsurf_row_for_star_LNA


      subroutine rsp_lsurf_resid_for_star_LNA(s, luminosity_resid_ad)
         type(star_info), pointer :: s
         type(auto_diff_real_star_order1), intent(out) :: luminosity_resid_ad
         type(auto_diff_real_star_order1) :: L1_ad, r1_ad, area_ad, rhs_ad, T_surf_ad
         real(dp) :: scale

         L1_ad = wrap_L_00(s, 1)
         r1_ad = wrap_r_00(s, 1)
         T_surf_ad = wrap_T_00(s, 1)
         area_ad = 4d0*pi*r1_ad*r1_ad
         rhs_ad = s% RSP2_Lsurf_factor*area_ad*clight*(crad*pow4(T_surf_ad))

         scale = maxval(s% L_start(1:s% nz))
         luminosity_resid_ad = (L1_ad - rhs_ad)/scale
      end subroutine rsp_lsurf_resid_for_star_LNA


      subroutine surface_temperature_resid_for_star_LNA(s, luminosity_resid_ad, ierr)
         type(star_info), pointer :: s
         type(auto_diff_real_star_order1), intent(out) :: luminosity_resid_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: lnT_bc_ad

         call surface_lnT_bc_for_star_LNA(s, lnT_bc_ad, ierr)
         if (ierr /= 0) return
         luminosity_resid_ad = lnT_bc_ad - wrap_lnT_00(s, 1)
      end subroutine surface_temperature_resid_for_star_LNA


      subroutine temperature_gradient_resid_for_star_LNA(s, k, luminosity_resid_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: luminosity_resid_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: dlnPdm_ad, Ppoint_ad, gradT_ad, &
            dlnTdm_ad, Tm1_ad, T00_ad, dT_ad, Tpoint_ad, lnTdiff_ad
         real(dp) :: delm, alfa

         ierr = 0
         if (s% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn) then
            luminosity_resid_ad = &
               s% gradT_ad(k)*(wrap_lnPeos_m1(s, k) - wrap_lnPeos_00(s, k)) - &
               (wrap_lnT_m1(s, k) - wrap_lnT_00(s, k))
            return
         end if

         if (s% use_dPrad_dm_form_of_T_gradient_eqn) then
            call dPrad_dm_resid_for_star_LNA(s, k, luminosity_resid_ad, ierr)
            return
         end if

         call star_LNA_eval_dlnPdm_qhse(s, k, dlnPdm_ad, Ppoint_ad, ierr)
         if (ierr /= 0) return

         gradT_ad = s% gradT_ad(k)
         dlnTdm_ad = dlnPdm_ad*gradT_ad

         Tm1_ad = wrap_T_m1(s, k)
         T00_ad = wrap_T_00(s, k)
         dT_ad = Tm1_ad - T00_ad
         alfa = s% dm(k - 1)/(s% dm(k - 1) + s% dm(k))
         Tpoint_ad = alfa*T00_ad + (1d0 - alfa)*Tm1_ad
         lnTdiff_ad = dT_ad/Tpoint_ad
         delm = 0.5d0*(s% dm(k) + s% dm(k - 1))

         luminosity_resid_ad = delm*dlnTdm_ad - lnTdiff_ad
      end subroutine temperature_gradient_resid_for_star_LNA


      subroutine dPrad_dm_resid_for_star_LNA(s, k, luminosity_resid_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: luminosity_resid_ad
         integer, intent(out) :: ierr
         real(dp) :: dm_bar, scale
         type(auto_diff_real_star_order1) :: L_ad, r_ad, area_ad, Lrad_ad, &
            opacity_face_ad, kap_face_ad, L0_ad, gradr_factor_ad, &
            T_m1_ad, T_00_ad, T4_m1_ad, T4_00_ad, &
            d_P_rad_expected_ad, d_P_rad_actual_ad, &
            T_face_ad, rho_face_ad, P_face_ad, Cp_face_ad, &
            ChiRho_face_ad, ChiT_face_ad, grada_face_ad, &
            flxR_ad, flxLambda_ad

         ierr = 0
         dm_bar = s% dm_bar(k)
         scale = s% energy_start(k)*s% rho_start(k)
         L_ad = wrap_L_00(s, k)
         r_ad = wrap_r_00(s, k)
         area_ad = 4d0*pi*pow2(r_ad)

         if (s% use_face_reconstruction) then
            call get_reconstructed_face_eos_kap_ad( &
               s, k, T_face_ad, rho_face_ad, P_face_ad, Cp_face_ad, &
               ChiRho_face_ad, ChiT_face_ad, grada_face_ad, opacity_face_ad, ierr)
            if (ierr /= 0) return
         else
            T_face_ad = get_T_face(s, k)
            P_face_ad = get_Peos_face(s, k)
            opacity_face_ad = get_kap_face(s, k)
         end if

         ! Match hydro_temperature.do1_alt_dlnT_dm_eqn in convective cells.
         gradr_factor_ad = get_effective_gradr_factor_ad(s, k)
         if (s% lnT(k)/ln10 <= s% max_logT_for_mlt .and. &
               s% mlt_mixing_type(k) == convective_mixing .and. &
               abs(gradr_factor_ad%val) > 1d-20) then
            L0_ad = get_Lrad_per_gradT_face_ad( &
               s, k, T_face_ad, P_face_ad, opacity_face_ad, gradr_factor_ad)
            Lrad_ad = L0_ad*s% gradT_ad(k)
         else
            Lrad_ad = L_ad
         end if

         kap_face_ad = opacity_face_ad
         if (kap_face_ad%val < s% min_kap_for_dPrad_dm_eqn) &
            kap_face_ad = s% min_kap_for_dPrad_dm_eqn

         d_P_rad_expected_ad = &
            -dm_bar*kap_face_ad*Lrad_ad/(clight*pow2(area_ad))

         T_m1_ad = wrap_T_m1(s, k)
         T_00_ad = wrap_T_00(s, k)
         T4_m1_ad = pow4(T_m1_ad)
         T4_00_ad = pow4(T_00_ad)
         d_P_rad_actual_ad = (crad/3d0)*(T4_m1_ad - T4_00_ad)

         if (s% use_flux_limiting_with_dPrad_dm_form) then
            flxR_ad = area_ad*abs(T4_m1_ad - T4_00_ad)/dm_bar/ &
               (kap_face_ad*0.5d0*(T4_m1_ad + T4_00_ad))
            flxLambda_ad = (6d0 + 3d0*flxR_ad)/ &
               (6d0 + (3d0 + flxR_ad)*flxR_ad)
            d_P_rad_expected_ad = d_P_rad_expected_ad/flxLambda_ad
         end if

         luminosity_resid_ad = &
            (d_P_rad_expected_ad - d_P_rad_actual_ad)/scale
      end subroutine dPrad_dm_resid_for_star_LNA


      subroutine surface_lnT_bc_for_star_LNA(s, lnT_bc_ad, ierr)
         type(star_info), pointer :: s
         type(auto_diff_real_star_order1), intent(out) :: lnT_bc_ad
         integer, intent(out) :: ierr
         real(dp) :: r, L, Teff, lnT_surf, lnP_surf, T_surf, &
            dlnTsurf_dL, dlnTsurf_dlnR, dlnTsurf_dlnM, dlnTsurf_dlnkap, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnM, dlnPsurf_dlnkap, &
            dlnT_bc_dlnd, dlnT_bc_dlnT, dlnT_bc_dL, dlnT_bc_dlnR, &
            dlnkap_dlnd, dlnkap_dlnT, dPinv_dlnd, dPinv_dlnT, &
            dP0, dT0, T_bc, dT0_dlnR, dT0_dlnT, dT0_dlnd, dT0_dL, &
            dlnT_bc_dlnTsurf, dlnT_bc_dT0, d_gradT_dlnR, &
            d_gradT_dlnT00, d_gradT_dlnd00, d_gradT_dL
         logical, parameter :: skip_partials = .false.
         logical, parameter :: need_P_surf = .true.
         logical, parameter :: need_T_surf = .true.

         ierr = 0
         call set_Teff_info_for_eqns(s, skip_partials, &
            need_P_surf, need_T_surf, r, L, Teff, &
            lnT_surf, dlnTsurf_dL, dlnTsurf_dlnR, dlnTsurf_dlnM, dlnTsurf_dlnkap, &
            lnP_surf, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnM, dlnPsurf_dlnkap, ierr)
         if (ierr /= 0) return

         T_surf = exp(lnT_surf)

         dP0 = 0d0
         if (.not. s% use_momentum_outer_BC) &
            dP0 = s% cgrav(1)*s% m_grav(1)*s% dm(1)/(8d0*pi*pow4(r))
         dT0 = dP0*s% gradT(1)*s% T(1)/s% Peos(1)
         T_bc = T_surf + dT0

         dT0_dlnR = 0d0
         dT0_dlnT = 0d0
         dT0_dlnd = 0d0
         dT0_dL = 0d0
         if (dP0 /= 0d0) then
            d_gradT_dlnR = s% gradT_ad(1)%d1Array(i_lnR_00)
            d_gradT_dlnT00 = s% gradT_ad(1)%d1Array(i_lnT_00)
            d_gradT_dlnd00 = s% gradT_ad(1)%d1Array(i_lnd_00)
            d_gradT_dL = s% gradT_ad(1)%d1Array(i_L_00)
            dT0_dlnR = -4d0*dT0 + dP0*d_gradT_dlnR*s% T(1)/s% Peos(1)
            dPinv_dlnT = -s% chiT_for_partials(1)/s% Peos(1)
            dT0_dlnT = dT0 + dP0*d_gradT_dlnT00*s% T(1)/s% Peos(1) + &
               dP0*s% gradT(1)*s% T(1)*dPinv_dlnT
            dPinv_dlnd = -s% chiRho_for_partials(1)/s% Peos(1)
            dT0_dlnd = dP0*d_gradT_dlnd00*s% T(1)/s% Peos(1) + &
               dP0*s% gradT(1)*s% T(1)*dPinv_dlnd
            dT0_dL = dP0*d_gradT_dL*s% T(1)/s% Peos(1)
         end if

         dlnT_bc_dT0 = 1d0/T_bc
         dlnT_bc_dlnTsurf = T_surf/T_bc
         dlnkap_dlnd = s% d_opacity_dlnd(1)/s% opacity(1)
         dlnkap_dlnT = s% d_opacity_dlnT(1)/s% opacity(1)

         dlnT_bc_dlnT = dlnT_bc_dlnTsurf*dlnTsurf_dlnkap*dlnkap_dlnT + &
            dlnT_bc_dT0*dT0_dlnT
         dlnT_bc_dlnd = dlnT_bc_dlnTsurf*dlnTsurf_dlnkap*dlnkap_dlnd + &
            dlnT_bc_dT0*dT0_dlnd
         dlnT_bc_dL = dlnT_bc_dlnTsurf*dlnTsurf_dL + dlnT_bc_dT0*dT0_dL
         dlnT_bc_dlnR = dlnT_bc_dlnTsurf*dlnTsurf_dlnR + dlnT_bc_dT0*dT0_dlnR

         call wrap(lnT_bc_ad, log(T_bc), &
            0d0, dlnT_bc_dlnd, 0d0, &
            0d0, dlnT_bc_dlnT, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, dlnT_bc_dlnR, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, dlnT_bc_dL, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0)
      end subroutine surface_lnT_bc_for_star_LNA


      subroutine assemble_rsp2_turbulent_rows(s, map, mtx, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(out) :: ierr
         integer :: k

         ierr = 0
         if (.not. s% RSP2_flag) return

         do k = 1, map%nz
            call assemble_rsp2_w_row(s, map, mtx, k, ierr)
            if (ierr /= 0) return
            call assemble_rsp2_Hp_row(s, map, mtx, k, ierr)
            if (ierr /= 0) return
         end do
      end subroutine assemble_rsp2_turbulent_rows


      ! RSP2 turbulent energy row:
      !   d etrb_k/dt = COUPL_k - dLt_k/dm - Ptrb_k*dVdt_k/dm
      !
      ! Static background linearized form:
      !   delta[COUPL - dLt/dm - Ptrb0*dVdt/dm] = sigma*delta etrb.
      ! The nonlinear hydro residual has an Eq term, but Eq is quadratic in the
      ! velocity gradient perturbation and is omitted from the static LNA.
      ! Forced nonturbulent cells use delta w = 0.
      subroutine assemble_rsp2_w_row(s, map, mtx, k, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         integer :: row, col_w
         type(auto_diff_real_star_order1) :: rhs_ad, inertia_ad

         ierr = 0
         row = matrix_index(map, k, lna_var_w)
         col_w = matrix_index(map, k, lna_var_w)
         if (row <= 0 .or. col_w <= 0) then
            write(*,'(a,i0)') 'star_LNA missing RSP2 w row variable at k = ', k
            ierr = -1
            return
         end if

         if (rsp2_forces_non_turbulent_cell(s, k)) then
            mtx%A(row, col_w) = 1d0
            return
         end if

         call rsp2_turbulent_energy_rhs_for_star_LNA(s, k, rhs_ad, ierr)
         if (ierr /= 0) return
         call add_ad_partials_to_A(map, mtx, row, k, 1d0, rhs_ad, ierr)
         if (ierr /= 0) return

         call rsp2_turbulent_energy_inertia_for_star_LNA(s, k, inertia_ad, ierr)
         if (ierr /= 0) return
         call add_ad_partials_to_B(map, mtx, row, k, 1d0, inertia_ad, ierr)
         if (ierr /= 0) return
      end subroutine assemble_rsp2_w_row


      ! RSP2 algebraic pressure scale height row:
      !   Hp_expected_k - Hp_k = 0
      !
      ! Matrix:
      !   A(row,:) receives the scaled AD partials of this residual.
      !   B(row,:) is unchanged.
      subroutine assemble_rsp2_Hp_row(s, map, mtx, k, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         integer :: row
         real(dp) :: scale
         type(auto_diff_real_star_order1) :: Hp_expected_ad, Hp_resid_ad

         ierr = 0
         row = matrix_index(map, k, lna_var_Hp)
         if (row <= 0) then
            write(*,'(a,i0)') 'star_LNA RSP2 Hp row is missing a variable at k = ', k
            ierr = -1
            return
         end if

         Hp_expected_ad = Hp_face_for_rsp2_eqn(s, k, ierr)
         if (ierr /= 0) return
         scale = 1d0/max(1d0, abs(Hp_expected_ad%val), abs(s% Hp_face(k)))
         Hp_resid_ad = scale*(Hp_expected_ad - wrap_Hp_00(s, k))
         call add_ad_partials_to_A(map, mtx, row, k, 1d0, Hp_resid_ad, ierr)
         if (ierr /= 0) return
      end subroutine assemble_rsp2_Hp_row


      subroutine assemble_tdc_turbulent_rows(s, map, mtx, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(out) :: ierr
         integer :: k

         ierr = 0
         if (.not. tdc_lna_active(s)) return

         do k = 1, map%nz
            call assemble_tdc_w_row(s, map, mtx, k, ierr)
            if (ierr /= 0) return
         end do
      end subroutine assemble_tdc_turbulent_rows


      ! TDC internal w row:
      !   velocity_rhs(A,rho,T,...) = sigma*velocity_inertia(A,rho,...)
      !
      ! with A = conv_vel/sqrt(2/3). The TDC terms are returned by
      ! turb:set_TDC_LNA and are internal to the LNA; this does not add a normal
      ! MESA/star TDC state variable. Forced nonturbulent faces use delta w = 0.
      subroutine assemble_tdc_w_row(s, map, mtx, k, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         integer :: row, col_w
         type(auto_diff_real_star_order1) :: luminosity_resid_ad, rhs_ad, &
            inertia_ad, L_rad_ad, L_conv_ad

         ierr = 0
         row = matrix_index(map, k, lna_var_w)
         col_w = matrix_index(map, k, lna_var_w)
         if (row <= 0 .or. col_w <= 0) then
            write(*,'(a,i0)') 'star_LNA missing TDC w row variable at k = ', k
            ierr = -1
            return
         end if

         if (tdc_zero_w_for_star_LNA(s, k)) then
            mtx%A(row, col_w) = 1d0
            return
         end if

         call tdc_relation_for_star_LNA(s, k, luminosity_resid_ad, rhs_ad, &
            inertia_ad, L_rad_ad, L_conv_ad, ierr)
         if (ierr /= 0) return
         call add_ad_partials_to_A(map, mtx, row, k, 1d0, rhs_ad, ierr)
         if (ierr /= 0) return
         call add_ad_partials_to_B(map, mtx, row, k, 1d0, inertia_ad, ierr)
         if (ierr /= 0) return
      end subroutine assemble_tdc_w_row

      ! Automatic differentiation to matrix plumbing.
      subroutine add_ad_partials_to_A(map, mtx, row, k, scale, ad, ierr)
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(in) :: row, k
         real(dp), intent(in) :: scale
         type(auto_diff_real_star_order1), intent(in) :: ad
         integer, intent(out) :: ierr

         call add_ad_partials_to_matrix(map, mtx%A, row, k, scale, ad, ierr)
      end subroutine add_ad_partials_to_A


      subroutine add_ad_partials_to_B(map, mtx, row, k, scale, ad, ierr)
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(inout) :: mtx
         integer, intent(in) :: row, k
         real(dp), intent(in) :: scale
         type(auto_diff_real_star_order1), intent(in) :: ad
         integer, intent(out) :: ierr

         call add_ad_partials_to_matrix(map, mtx%B, row, k, scale, ad, ierr)
      end subroutine add_ad_partials_to_B


      subroutine add_ad_partials_to_matrix(map, mat, row, k, scale, ad, ierr)
         type(star_LNA_var_map), intent(in) :: map
         real(dp), intent(inout) :: mat(:,:)
         integer, intent(in) :: row, k
         real(dp), intent(in) :: scale
         type(auto_diff_real_star_order1), intent(in) :: ad
         integer, intent(out) :: ierr
         integer :: iad, kk, var_id, col
         real(dp) :: coeff

         ierr = 0
         do iad = 1, auto_diff_star_num_vars
            coeff = scale*ad%d1Array(iad)
            if (coeff == 0d0) cycle

            call ad_index_to_star_LNA_var(map, k, iad, kk, var_id)
            if (var_id == 0) then
               write(*,'(3a)') 'star_LNA does not map automatic differentiation variable ', &
                  trim(auto_diff_star_d1_names(iad)), ' into the LNA matrix.'
               ierr = -1
               return
            end if
            if (kk < 1 .or. kk > map%nz) cycle

            col = matrix_index(map, kk, var_id)
            if (col <= 0) then
               write(*,'(3a)') 'star_LNA active equation has unsupported variable ', &
                  trim(var_name(var_id)), ' in its automatic differentiation partials.'
               ierr = -1
               return
            end if
            mat(row, col) = mat(row, col) + coeff
         end do
      end subroutine add_ad_partials_to_matrix


      subroutine ad_index_to_star_LNA_var(map, k, iad, kk, var_id)
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: k, iad
         integer, intent(out) :: kk, var_id

         kk = k
         var_id = 0
         select case (iad)
         case (i_lnd_m1)
            kk = k - 1; var_id = lna_var_lnd
         case (i_lnd_00)
            kk = k; var_id = lna_var_lnd
         case (i_lnd_p1)
            kk = k + 1; var_id = lna_var_lnd
         case (i_lnT_m1)
            kk = k - 1; var_id = lna_var_lnT
         case (i_lnT_00)
            kk = k; var_id = lna_var_lnT
         case (i_lnT_p1)
            kk = k + 1; var_id = lna_var_lnT
         case (i_w_m1)
            kk = k - 1; var_id = lna_var_w
         case (i_w_00)
            kk = k; var_id = lna_var_w
         case (i_w_p1)
            kk = k + 1; var_id = lna_var_w
         case (i_lnR_m1)
            kk = k - 1; var_id = lna_var_lnR
         case (i_lnR_00)
            kk = k; var_id = lna_var_lnR
         case (i_lnR_p1)
            kk = k + 1; var_id = lna_var_lnR
         case (i_v_m1)
            kk = k - 1; var_id = velocity_var_for_star_LNA(map)
         case (i_v_00)
            kk = k; var_id = velocity_var_for_star_LNA(map)
         case (i_v_p1)
            kk = k + 1; var_id = velocity_var_for_star_LNA(map)
         case (i_L_m1)
            kk = k - 1; var_id = lna_var_L
         case (i_L_00)
            kk = k; var_id = lna_var_L
         case (i_L_p1)
            kk = k + 1; var_id = lna_var_L
         case (i_Hp_m1)
            kk = k - 1; var_id = lna_var_Hp
         case (i_Hp_00)
            kk = k; var_id = lna_var_Hp
         case (i_Hp_p1)
            kk = k + 1; var_id = lna_var_Hp
         end select
      end subroutine ad_index_to_star_LNA_var

      ! Dense generalized eigensolver and algebraic elimination.
      subroutine solve_dense_star_LNA(s, map, mtx, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(in) :: mtx
         integer, intent(out) :: ierr
         type(star_LNA_matrix) :: scaled_mtx
         integer :: info, lwork, ndyn, nalg
         integer, allocatable :: dyn_idx(:), alg_idx(:)
         real(dp), allocatable :: Ared(:,:), Bred(:,:), alg_from_dyn(:,:), &
            Awork(:,:), Bwork(:,:), alphar(:), alphai(:), beta(:), &
            vl(:,:), vr_red(:,:), work(:), full_col_scale(:)

         ierr = 0
         call scale_star_LNA_matrix_for_solver(mtx, scaled_mtx, full_col_scale, ierr)
         if (ierr /= 0) return
         call build_reduced_star_LNA_problem( &
            map, scaled_mtx, Ared, Bred, alg_from_dyn, dyn_idx, alg_idx, ndyn, nalg, ierr)
         call free_star_LNA_matrix(scaled_mtx)
         if (ierr /= 0) then
            deallocate(full_col_scale)
            return
         end if
         write(*,'(a,i0,a,i0)') 'star_LNA: reduced dynamic variables = ', &
            ndyn, ', algebraic variables eliminated = ', nalg

         allocate( &
            Awork(ndyn, ndyn), Bwork(ndyn, ndyn), &
            alphar(ndyn), alphai(ndyn), beta(ndyn), &
            vl(ndyn, ndyn), vr_red(ndyn, ndyn), stat=ierr)
         if (ierr /= 0) then
            deallocate(full_col_scale)
            deallocate(Ared, Bred, alg_from_dyn, dyn_idx, alg_idx)
            return
         end if

         lwork = max(1, 8*ndyn)
         allocate(work(lwork), stat=ierr)
         if (ierr /= 0) then
            deallocate(Awork, Bwork, alphar, alphai, beta, vl, vr_red)
            deallocate(full_col_scale)
            deallocate(Ared, Bred, alg_from_dyn, dyn_idx, alg_idx)
            return
         end if

         Awork = Ared
         Bwork = Bred
         call DGGEV( &
            'N', 'V', ndyn, Awork, ndyn, Bwork, ndyn, &
            alphar, alphai, beta, vl, ndyn, vr_red, ndyn, work, lwork, info)
         if (info /= 0) then
            write(*,'(a,i0)') 'star_LNA: LAPACK/DGGEV failed with info = ', info
            ierr = -1
         else
            call finish_dense_star_LNA_solution( &
               s, map, alphar, alphai, beta, vr_red, dyn_idx, alg_idx, &
               alg_from_dyn, full_col_scale, ierr)
         end if

         deallocate(Awork, Bwork, alphar, alphai, beta, vl, vr_red, work)
         deallocate(full_col_scale)
         deallocate(Ared, Bred, alg_from_dyn, dyn_idx, alg_idx)
      end subroutine solve_dense_star_LNA


      subroutine finish_dense_star_LNA_solution( &
            s, map, alphar, alphai, beta, vr_red, dyn_idx, alg_idx, &
            alg_from_dyn, full_col_scale, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         real(dp), intent(in) :: alphar(:), alphai(:), beta(:), vr_red(:,:)
         real(dp), intent(in) :: alg_from_dyn(:,:), full_col_scale(:)
         integer, intent(in) :: dyn_idx(:), alg_idx(:)
         integer, intent(out) :: ierr
         integer :: num_modes
         integer, allocatable :: mode_indices(:)
         real(dp), allocatable :: vr_full(:,:)

         ierr = 0
         call reconstruct_star_LNA_eigenvectors( &
            map, dyn_idx, alg_idx, alg_from_dyn, vr_red, vr_full, ierr)
         if (ierr == 0) call unscale_star_LNA_full_eigenvectors(vr_full, full_col_scale)
         if (ierr == 0) &
            call select_star_LNA_modes(s, alphar, alphai, beta, mode_indices, num_modes, ierr)
         if (ierr == 0 .and. s% star_LNA_write_period_growth) &
            call write_star_LNA_raw_modes(s, map, alphar, alphai, beta, ierr)
         if (ierr == 0) call print_star_LNA_period_growth( &
            s, alphar, alphai, beta, mode_indices, num_modes)
         if (ierr == 0 .and. s% star_LNA_write_period_growth) &
            call write_star_LNA_period_growth( &
               s, map, alphar, alphai, beta, mode_indices, num_modes, ierr)
         if (ierr == 0 .and. &
               (s% star_LNA_write_eigenfunctions .or. s% star_LNA_write_work_integrals)) &
            call write_star_LNA_eigenfunctions( &
               s, map, alphar, alphai, beta, vr_full, mode_indices, num_modes, ierr)
         if (ierr == 0 .and. s% star_LNA_set_initial_velocity) &
            call set_star_LNA_initial_velocity( &
               s, map, alphar, alphai, beta, vr_full, mode_indices, num_modes, ierr)

         if (allocated(mode_indices)) deallocate(mode_indices)
         if (allocated(vr_full)) deallocate(vr_full)
      end subroutine finish_dense_star_LNA_solution


      subroutine scale_star_LNA_matrix_for_solver(mtx, scaled_mtx, full_col_scale, ierr)
         type(star_LNA_matrix), intent(in) :: mtx
         type(star_LNA_matrix), intent(out) :: scaled_mtx
         real(dp), allocatable, intent(out) :: full_col_scale(:)
         integer, intent(out) :: ierr
         integer :: i, j, n
         real(dp) :: row_norm, col_norm, scale

         ierr = 0
         n = mtx%n
         scaled_mtx%n = n
         allocate(scaled_mtx%A(n, n), scaled_mtx%B(n, n), full_col_scale(n), stat=ierr)
         if (ierr /= 0) return

         scaled_mtx%A = mtx%A
         scaled_mtx%B = mtx%B
         full_col_scale = 1d0

         do i = 1, n
            row_norm = max(maxval(abs(scaled_mtx%A(i, :))), maxval(abs(scaled_mtx%B(i, :))))
            if (row_norm > tiny(1d0)) then
               scale = 1d0/row_norm
               scaled_mtx%A(i, :) = scale*scaled_mtx%A(i, :)
               scaled_mtx%B(i, :) = scale*scaled_mtx%B(i, :)
            end if
         end do

         do j = 1, n
            col_norm = max(maxval(abs(scaled_mtx%A(:, j))), maxval(abs(scaled_mtx%B(:, j))))
            if (col_norm > tiny(1d0)) then
               full_col_scale(j) = 1d0/col_norm
               scaled_mtx%A(:, j) = full_col_scale(j)*scaled_mtx%A(:, j)
               scaled_mtx%B(:, j) = full_col_scale(j)*scaled_mtx%B(:, j)
            end if
         end do
      end subroutine scale_star_LNA_matrix_for_solver


      ! Algebraic elimination:
      !   Aaa*x_alg + Aad*x_dyn = 0
      !   x_alg = -inv(Aaa)*Aad*x_dyn
      !
      ! Reduced problem:
      !   (Add - Ada*inv(Aaa)*Aad)*x_dyn =
      !      sigma*(Bdd - Bda*inv(Aaa)*Aad)*x_dyn
      !
      ! The minus sign is applied later when reconstructing the full eigenvector.
      subroutine build_reduced_star_LNA_problem( &
            map, mtx, Ared, Bred, alg_from_dyn, dyn_idx, alg_idx, ndyn, nalg, ierr)
         type(star_LNA_var_map), intent(in) :: map
         type(star_LNA_matrix), intent(in) :: mtx
         real(dp), allocatable, intent(out) :: Ared(:,:), Bred(:,:), alg_from_dyn(:,:)
         integer, allocatable, intent(out) :: dyn_idx(:), alg_idx(:)
         integer, intent(out) :: ndyn, nalg
         integer, intent(out) :: ierr
         integer :: info
         integer, allocatable :: ipiv(:)
         real(dp), allocatable :: Aaa(:,:), A_ad(:,:), A_da(:,:), B_da(:,:)

         ierr = 0
         call partition_star_LNA_indices(map, dyn_idx, alg_idx, ndyn, nalg, ierr)
         if (ierr /= 0) return

         if (ndyn < 1) then
            write(*,'(a)') 'star_LNA: no dynamic variables remain after algebraic partition.'
            ierr = -1
            return
         end if

         allocate(Ared(ndyn, ndyn), Bred(ndyn, ndyn), alg_from_dyn(nalg, ndyn), stat=ierr)
         if (ierr /= 0) return

         Ared = mtx%A(dyn_idx, dyn_idx)
         Bred = mtx%B(dyn_idx, dyn_idx)
         if (nalg < 1) then
            alg_from_dyn = 0d0
            return
         end if
         if (any(mtx%B(alg_idx, :) /= 0d0)) then
            write(*,'(a)') 'star_LNA: algebraic rows contain time derivative terms.'
            ierr = -1
            deallocate(Ared, Bred, alg_from_dyn)
            return
         end if

         allocate( &
            Aaa(nalg, nalg), A_ad(nalg, ndyn), A_da(ndyn, nalg), B_da(ndyn, nalg), &
            ipiv(nalg), stat=ierr)
         if (ierr /= 0) return

         Aaa = mtx%A(alg_idx, alg_idx)
         A_ad = mtx%A(alg_idx, dyn_idx)
         A_da = mtx%A(dyn_idx, alg_idx)
         B_da = mtx%B(dyn_idx, alg_idx)
         alg_from_dyn = A_ad

         call DGESV(nalg, ndyn, Aaa, nalg, ipiv, alg_from_dyn, nalg, info)
         if (info /= 0) then
            write(*,'(a,i0)') &
               'star_LNA: LAPACK/DGESV failed while eliminating algebraic rows, info = ', info
            ierr = -1
            deallocate(Aaa, A_ad, A_da, B_da, ipiv)
            return
         end if

         Ared = Ared - matmul(A_da, alg_from_dyn)
         if (any(B_da /= 0d0)) Bred = Bred - matmul(B_da, alg_from_dyn)

         deallocate(Aaa, A_ad, A_da, B_da, ipiv)
      end subroutine build_reduced_star_LNA_problem


      subroutine partition_star_LNA_indices(map, dyn_idx, alg_idx, ndyn, nalg, ierr)
         type(star_LNA_var_map), intent(in) :: map
         integer, allocatable, intent(out) :: dyn_idx(:), alg_idx(:)
         integer, intent(out) :: ndyn, nalg
         integer, intent(out) :: ierr
         integer :: k, slot, idx, var_id

         ierr = 0
         ndyn = 0
         nalg = 0
         do k = 1, map%nz
            do slot = 1, map%nvar_per_zone
               var_id = map%var_id(slot)
               if (star_LNA_var_is_dynamic(var_id)) then
                  ndyn = ndyn + 1
               else
                  nalg = nalg + 1
               end if
            end do
         end do

         allocate(dyn_idx(ndyn), alg_idx(nalg), stat=ierr)
         if (ierr /= 0) return

         ndyn = 0
         nalg = 0
         do k = 1, map%nz
            do slot = 1, map%nvar_per_zone
               var_id = map%var_id(slot)
               idx = (k - 1)*map%nvar_per_zone + slot
               if (star_LNA_var_is_dynamic(var_id)) then
                  ndyn = ndyn + 1
                  dyn_idx(ndyn) = idx
               else
                  nalg = nalg + 1
                  alg_idx(nalg) = idx
               end if
            end do
         end do
      end subroutine partition_star_LNA_indices


      logical function star_LNA_var_is_dynamic(var_id) result(is_dynamic)
         integer, intent(in) :: var_id

         select case (var_id)
         case (lna_var_lnR, lna_var_v, lna_var_u, lna_var_lnT, lna_var_w)
            is_dynamic = .true.
         case default
            is_dynamic = .false.
         end select
      end function star_LNA_var_is_dynamic


      subroutine reconstruct_star_LNA_eigenvectors( &
            map, dyn_idx, alg_idx, alg_from_dyn, vr_red, vr_full, ierr)
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: dyn_idx(:), alg_idx(:)
         real(dp), intent(in) :: alg_from_dyn(:,:), vr_red(:,:)
         real(dp), allocatable, intent(out) :: vr_full(:,:)
         integer, intent(out) :: ierr
         integer :: j, ndyn, nalg

         ierr = 0
         ndyn = size(dyn_idx)
         nalg = size(alg_idx)
         allocate(vr_full(map%nvar_total, size(vr_red, 2)), stat=ierr)
         if (ierr /= 0) return

         vr_full = 0d0
         do j = 1, size(vr_red, 2)
            vr_full(dyn_idx, j) = vr_red(1:ndyn, j)
            if (nalg > 0) vr_full(alg_idx, j) = -matmul(alg_from_dyn, vr_red(1:ndyn, j))
         end do
      end subroutine reconstruct_star_LNA_eigenvectors


      subroutine unscale_star_LNA_full_eigenvectors(vr_full, full_col_scale)
         real(dp), intent(inout) :: vr_full(:,:)
         real(dp), intent(in) :: full_col_scale(:)
         integer :: i

         do i = 1, size(vr_full, 1)
            vr_full(i, :) = full_col_scale(i)*vr_full(i, :)
         end do
      end subroutine unscale_star_LNA_full_eigenvectors

      ! Mode selection, period/growth output, eigenfunctions, and kicks.
      subroutine select_star_LNA_modes(s, alphar, alphai, beta, mode_indices, num_modes, ierr)
         type(star_info), pointer :: s
         real(dp), intent(in) :: alphar(:), alphai(:), beta(:)
         integer, allocatable, intent(out) :: mode_indices(:)
         integer, intent(out) :: num_modes
         integer, intent(out) :: ierr
         integer :: best, j, op_err
         real(dp) :: sigma_re, sigma_im, omega, best_omega, frequency_uHz, logKE_per_cycle
         logical :: found_first_mode
         logical, allocatable :: used(:)

         ierr = 0
         num_modes = 0
         allocate(mode_indices(min(s% star_LNA_num_modes, size(alphar))), stat=ierr)
         if (ierr /= 0) return

         allocate(used(size(alphar)), stat=ierr)
         if (ierr /= 0) then
            deallocate(mode_indices)
            return
         end if
         used = .false.
         found_first_mode = .false.
         do
            if (num_modes >= size(mode_indices)) exit
            best = 0
            best_omega = huge(1d0)
            do j = 1, size(alphar)
               if (used(j)) cycle
               call finite_sigma_from_star_LNA_eigenvalue(alphar(j), alphai(j), beta(j), &
                  sigma_re, sigma_im, op_err)
               if (op_err /= 0) cycle
               omega = sigma_im
               if (omega <= 1d-99) cycle
               frequency_uHz = frequency_uHz_for_star_LNA(sigma_im)
               if (frequency_uHz < s% star_LNA_min_mode_frequency_uHz) cycle
               if (omega < best_omega) then
                  best = j
                  best_omega = omega
               end if
            end do
            if (best == 0) exit

            used(best) = .true.
            call finite_sigma_from_star_LNA_eigenvalue(alphar(best), alphai(best), beta(best), &
               sigma_re, sigma_im, op_err)
            if (op_err /= 0) cycle
            logKE_per_cycle = rsp_logKE_per_cycle_for_star_LNA(sigma_re, sigma_im)
            if (s% star_LNA_max_abs_mode_eta > 0d0 .and. &
                  abs(logKE_per_cycle) > s% star_LNA_max_abs_mode_eta) cycle
            if (.not. found_first_mode) then
               if (logKE_per_cycle <= s% star_LNA_min_first_mode_eta) cycle
               found_first_mode = .true.
            end if
            if (logKE_per_cycle < s% star_LNA_min_mode_eta) cycle

            num_modes = num_modes + 1
            mode_indices(num_modes) = best
         end do

         deallocate(used)
      end subroutine select_star_LNA_modes


      subroutine print_star_LNA_period_growth(s, alphar, alphai, beta, mode_indices, num_modes)
         type(star_info), pointer :: s
         real(dp), intent(in) :: alphar(:), alphai(:), beta(:)
         integer, intent(in) :: mode_indices(:), num_modes
         integer :: mode, idx, op_err
         real(dp) :: sigma_re, sigma_im, period_days, logKE_per_cycle
         real(dp) :: ke_growth_per_period, grekm_growth_per_period, amp_growth_per_period

         write(*,'(a)') 'star_LNA modes:'
         write(*,'(a)') &
            ' mode     P(days)      logKE/cyc      KE frac        GREKM       amp frac'
         if (num_modes < 1) then
            write(*,'(a)') 'star_LNA: no modes passed the frequency/logKE_per_cycle selection.'
            return
         end if
         do mode = 1, num_modes
            idx = mode_indices(mode)
            call finite_sigma_from_star_LNA_eigenvalue(alphar(idx), alphai(idx), beta(idx), &
               sigma_re, sigma_im, op_err)
            if (op_err /= 0) cycle
            period_days = period_days_for_star_LNA(sigma_im)
            logKE_per_cycle = rsp_logKE_per_cycle_for_star_LNA(sigma_re, sigma_im)
            ke_growth_per_period = ke_fractional_growth_per_period_for_star_LNA(sigma_re, sigma_im)
            grekm_growth_per_period = rsp_grekm_growth_per_period_for_star_LNA(sigma_re, sigma_im)
            amp_growth_per_period = amp_fractional_growth_per_period_for_star_LNA(sigma_re, sigma_im)
            write(*,'(1x,i5,5(1x,1pe13.5))') mode - 1, period_days, logKE_per_cycle, &
               ke_growth_per_period, grekm_growth_per_period, amp_growth_per_period
         end do
      end subroutine print_star_LNA_period_growth


      subroutine write_star_LNA_raw_modes(s, map, alphar, alphai, beta, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         real(dp), intent(in) :: alphar(:), alphai(:), beta(:)
         integer, intent(out) :: ierr
         integer :: io, raw_rank, best, j, op_err
         real(dp) :: sigma_re, sigma_im, omega, best_omega
         real(dp) :: period_days, pulsation_Q_days, logKE_per_cycle
         real(dp) :: ke_growth_per_period, grekm_growth_per_period, amp_growth_per_period
         logical, allocatable :: used(:)
         character(len=512) :: filename

         ierr = 0
         call star_LNA_output_filename(s, '_raw_modes.data', filename)
         open(newunit=io, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) return

         write(io,'(a)') '# star_LNA raw finite positive frequency eigenvalues'
         write(io,'(a)') '# sorted by increasing sigma_imag before acoustic branch selection'
         write(io,'(a,1x,i0)') '# nvar_total', map%nvar_total
         write(io,'(a,1x,i0)') '# nvar_dynamic', size(alphar)
         write(io,'(a)') '# eigen indices refer to the reduced dynamic variable eigenproblem'
         write(io,'(a)') &
            '# raw_rank reduced_eigen_index sigma_real sigma_imag ' // &
            'period_days pulsation_constant_Q_days W_rad_per_sec logKE_per_cycle beta ' // &
            'ke_fractional_growth_per_period grekm_growth_per_period ' // &
            'amplitude_fractional_growth_per_period'

         allocate(used(size(alphar)), stat=ierr)
         if (ierr /= 0) then
            close(io)
            return
         end if
         used = .false.
         raw_rank = 0

         do
            best = 0
            best_omega = huge(1d0)
            do j = 1, size(alphar)
               if (used(j)) cycle
               call finite_sigma_from_star_LNA_eigenvalue(alphar(j), alphai(j), beta(j), &
                  sigma_re, sigma_im, op_err)
               if (op_err /= 0) cycle
               omega = sigma_im
               if (omega <= 1d-99) cycle
               if (omega < best_omega) then
                  best = j
                  best_omega = omega
               end if
            end do
            if (best == 0) exit

            used(best) = .true.
            call finite_sigma_from_star_LNA_eigenvalue(alphar(best), alphai(best), beta(best), &
               sigma_re, sigma_im, op_err)
            if (op_err /= 0) cycle
            raw_rank = raw_rank + 1
            omega = abs(sigma_im)
            period_days = period_days_for_star_LNA(sigma_im)
            pulsation_Q_days = pulsation_constant_Q_days_for_star_LNA(s, period_days)
            logKE_per_cycle = rsp_logKE_per_cycle_for_star_LNA(sigma_re, sigma_im)
            ke_growth_per_period = ke_fractional_growth_per_period_for_star_LNA(sigma_re, sigma_im)
            grekm_growth_per_period = rsp_grekm_growth_per_period_for_star_LNA(sigma_re, sigma_im)
            amp_growth_per_period = amp_fractional_growth_per_period_for_star_LNA(sigma_re, sigma_im)
            write(io,'(2(i0,1x),10(1pe24.16,1x))') &
               raw_rank, best, &
               sigma_re, sigma_im, period_days, pulsation_Q_days, omega, logKE_per_cycle, beta(best), &
               ke_growth_per_period, grekm_growth_per_period, amp_growth_per_period
         end do

         deallocate(used)
         close(io)
         ierr = 0
      end subroutine write_star_LNA_raw_modes


      subroutine write_star_LNA_period_growth(s, map, alphar, alphai, beta, &
            mode_indices, num_modes, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         real(dp), intent(in) :: alphar(:), alphai(:), beta(:)
         integer, intent(in) :: mode_indices(:), num_modes
         integer, intent(out) :: ierr
         integer :: io, mode, idx, op_err, selected
         real(dp) :: sigma_re, sigma_im, omega, period_sec, period_days, &
            pulsation_Q_days, logKE_per_cycle
         real(dp) :: ke_growth_per_period, grekm_growth_per_period, amp_growth_per_period
         character(len=512) :: filename

         ierr = 0
         call star_LNA_output_filename(s, '_period_growth.data', filename)
         open(newunit=io, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) return

         write(io,'(a)') '# star_LNA mode summary'
         write(io,'(a)') '# reduced generalized eigenproblem: A*x = sigma*B*x'
         write(io,'(a,1x,i0)') '# nvar_total', map%nvar_total
         write(io,'(a,1x,i0)') '# nvar_dynamic', size(alphar)
         write(io,'(a,1x,1pe24.16)') '# star_LNA_min_mode_frequency_uHz', &
            s% star_LNA_min_mode_frequency_uHz
         write(io,'(a,1x,1pe24.16)') '# star_LNA_min_first_mode_eta', &
            s% star_LNA_min_first_mode_eta
         write(io,'(a,1x,1pe24.16)') '# star_LNA_min_mode_eta', &
            s% star_LNA_min_mode_eta
         write(io,'(a,1x,1pe24.16)') '# star_LNA_max_abs_mode_eta', &
            s% star_LNA_max_abs_mode_eta
         write(io,'(a)') '# legacy *_eta controls filter logKE_per_cycle'
         write(io,'(a,1x,l1)') '# star_LNA_set_initial_velocity', &
            s% star_LNA_set_initial_velocity
         write(io,'(a,1x,i0)') '# star_LNA_mode_for_period', &
            s% star_LNA_mode_for_period
         write(io,'(a,1x,i0)') '# star_LNA_kick_mode_1', &
            s% star_LNA_kick_mode_1
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_fraction_1', &
            s% star_LNA_kick_fraction_1
         write(io,'(a,1x,i0)') '# star_LNA_kick_mode_2', &
            s% star_LNA_kick_mode_2
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_fraction_2', &
            s% star_LNA_kick_fraction_2
         write(io,'(a,1x,i0)') '# star_LNA_kick_mode_3', &
            s% star_LNA_kick_mode_3
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_fraction_3', &
            s% star_LNA_kick_fraction_3
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_vsurf_km_per_sec', &
            s% star_LNA_kick_vsurf_km_per_sec
         write(io,'(a)') '# mode and star_LNA_kick_mode_* start at one; control_mode_index starts at zero'
         write(io,'(a)') &
            '# logKE_per_cycle = 4*pi*sigma_real/sigma_imag, the RSP LINA logarithmic kinetic energy convention'
         write(io,'(a)') &
            '# ke_fractional_growth_per_period is exp(logKE_per_cycle)-1; ' // &
            'grekm_growth_per_period is 2*tanh(logKE_per_cycle/2)'
         write(io,'(a)') &
            '# amplitude_fractional_growth_per_period is exp(logKE_per_cycle/2)-1'
         write(io,'(a)') '# eigen indices refer to the reduced dynamic variable eigenproblem'
         write(io,'(a)') &
            '# pulsation_constant_Q_days follows RSP LINA: P_days*sqrt((M/Msun)*(Rsun/R)^3)'
         write(io,'(a)') &
            '# mode control_mode_index selected_for_initial_velocity reduced_eigen_index ' // &
            'sigma_real sigma_imag period_days pulsation_constant_Q_days W_rad_per_sec ' // &
            'logKE_per_cycle beta ' // &
            'ke_fractional_growth_per_period grekm_growth_per_period ' // &
            'amplitude_fractional_growth_per_period'

         do mode = 1, num_modes
            idx = mode_indices(mode)
            call finite_sigma_from_star_LNA_eigenvalue(alphar(idx), alphai(idx), beta(idx), &
               sigma_re, sigma_im, op_err)
            if (op_err /= 0) cycle
            period_sec = period_sec_for_star_LNA(sigma_im)
            omega = abs(sigma_im)
            period_days = period_sec/86400d0
            pulsation_Q_days = pulsation_constant_Q_days_for_star_LNA(s, period_days)
            logKE_per_cycle = rsp_logKE_per_cycle_for_star_LNA(sigma_re, sigma_im)
            ke_growth_per_period = ke_fractional_growth_per_period_for_star_LNA(sigma_re, sigma_im)
            grekm_growth_per_period = rsp_grekm_growth_per_period_for_star_LNA(sigma_re, sigma_im)
            amp_growth_per_period = amp_fractional_growth_per_period_for_star_LNA(sigma_re, sigma_im)
            selected = merge(1, 0, selected_for_star_LNA_initial_velocity(s, mode))
            write(io,'(4(i0,1x),10(1pe24.16,1x))') &
               mode, mode - 1, selected, idx, &
               sigma_re, sigma_im, period_days, pulsation_Q_days, omega, logKE_per_cycle, beta(idx), &
               ke_growth_per_period, grekm_growth_per_period, amp_growth_per_period
         end do

         close(io)
         ierr = 0
      end subroutine write_star_LNA_period_growth


      subroutine write_star_LNA_eigenfunctions(s, map, alphar, alphai, beta, vr, &
            mode_indices, num_modes, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         real(dp), intent(in) :: alphar(:), alphai(:), beta(:), vr(:,:)
         integer, intent(in) :: mode_indices(:), num_modes
         integer, intent(out) :: ierr
         integer :: mode, idx, op_err
         real(dp) :: sigma_re, sigma_im
         complex(dp), allocatable :: eigenvector(:)

         ierr = 0
         allocate(eigenvector(map%nvar_total), stat=ierr)
         if (ierr /= 0) return

         do mode = 1, num_modes
            idx = mode_indices(mode)
            call finite_sigma_from_star_LNA_eigenvalue(alphar(idx), alphai(idx), beta(idx), &
               sigma_re, sigma_im, op_err)
            if (op_err /= 0) cycle
            call load_star_LNA_eigenvector(vr, idx, alphai(idx), eigenvector, op_err)
            if (op_err /= 0) cycle
            call normalize_star_LNA_eigenvector(map, eigenvector)
            if (s% star_LNA_write_eigenfunctions) then
               call write_one_star_LNA_eigenfunction(s, map, mode, idx, &
                  sigma_re, sigma_im, eigenvector, ierr)
               if (ierr /= 0) exit
            end if
            if (s% star_LNA_write_work_integrals) then
               call write_one_star_LNA_work(s, map, mode, idx, &
                  sigma_re, sigma_im, eigenvector, ierr)
               if (ierr /= 0) exit
            end if
         end do

         deallocate(eigenvector)
      end subroutine write_star_LNA_eigenfunctions


      subroutine set_star_LNA_initial_velocity(s, map, alphar, alphai, beta, vr, &
            mode_indices, num_modes, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         real(dp), intent(in) :: alphar(:), alphai(:), beta(:), vr(:,:)
         integer, intent(in) :: mode_indices(:), num_modes
         integer, intent(out) :: ierr
         integer :: period_mode, idx, op_err, largest_active_kick_mode
         real(dp) :: sigma_re, sigma_im, period_sec, period_days
         real(dp) :: active_fraction_sum
         complex(dp), allocatable :: kick_eigenvector(:)

         ierr = 0
         period_mode = s% star_LNA_mode_for_period + 1
         if (period_mode < 1 .or. period_mode > num_modes) then
            write(*,'(a,i0,a,i0)') 'star_LNA_mode_for_period = ', &
               s% star_LNA_mode_for_period, ' but selected num_modes = ', num_modes
            ierr = -1
            return
         end if

         largest_active_kick_mode = s% star_LNA_kick_mode_1
         if (s% star_LNA_kick_mode_2 > 0 .and. s% star_LNA_kick_fraction_2 > 0d0) &
            largest_active_kick_mode = max(largest_active_kick_mode, s% star_LNA_kick_mode_2)
         if (s% star_LNA_kick_mode_3 > 0 .and. s% star_LNA_kick_fraction_3 > 0d0) &
            largest_active_kick_mode = max(largest_active_kick_mode, s% star_LNA_kick_mode_3)
         if (largest_active_kick_mode > num_modes) then
            write(*,'(a,i0,a,i0)') 'largest star_LNA kick mode is ', &
               largest_active_kick_mode, &
               ' but selected num_modes = ', num_modes
            ierr = -1
            return
         end if

         allocate(kick_eigenvector(map%nvar_total), stat=ierr)
         if (ierr /= 0) return
         kick_eigenvector = (0d0, 0d0)

         active_fraction_sum = s% star_LNA_kick_fraction_1
         call add_star_LNA_kick_component( &
            s, map, alphai, vr, mode_indices, s% star_LNA_kick_mode_1, &
            s% star_LNA_kick_fraction_1, kick_eigenvector, ierr)
         if (ierr /= 0) then
            deallocate(kick_eigenvector)
            return
         end if

         if (s% star_LNA_kick_mode_2 > 0 .and. s% star_LNA_kick_fraction_2 > 0d0) then
            active_fraction_sum = active_fraction_sum + s% star_LNA_kick_fraction_2
            call add_star_LNA_kick_component( &
               s, map, alphai, vr, mode_indices, s% star_LNA_kick_mode_2, &
               s% star_LNA_kick_fraction_2, kick_eigenvector, ierr)
            if (ierr /= 0) then
               deallocate(kick_eigenvector)
               return
            end if
         end if

         if (s% star_LNA_kick_mode_3 > 0 .and. s% star_LNA_kick_fraction_3 > 0d0) then
            active_fraction_sum = active_fraction_sum + s% star_LNA_kick_fraction_3
            call add_star_LNA_kick_component( &
               s, map, alphai, vr, mode_indices, s% star_LNA_kick_mode_3, &
               s% star_LNA_kick_fraction_3, kick_eigenvector, ierr)
            if (ierr /= 0) then
               deallocate(kick_eigenvector)
               return
            end if
         end if

         if (active_fraction_sum <= 0d0) then
            write(*,'(a)') 'star_LNA active kick fractions sum to zero.'
            ierr = -1
            deallocate(kick_eigenvector)
            return
         end if
         kick_eigenvector = kick_eigenvector/active_fraction_sum

         idx = mode_indices(period_mode)
         call finite_sigma_from_star_LNA_eigenvalue(alphar(idx), alphai(idx), beta(idx), &
            sigma_re, sigma_im, op_err)
         if (op_err /= 0) then
            ierr = -1
         else
            call set_star_LNA_velocity_from_eigenvector(s, map, kick_eigenvector, ierr)
         end if

         if (ierr == 0) then
            period_sec = period_sec_for_star_LNA(sigma_im)
            period_days = period_sec/86400d0
            if (s% RSP_flag) then
               s% rsp_period = period_sec
               s% RSP_have_set_velocities = .true.
            end if
            if (s% RSP2_flag) s% RSP2_period = period_sec
            write(*,'(a,3(i0,1x),a,3(1pe14.6,1x),a,1pe14.6,a,1pe14.6)') &
               'star_LNA: set initial velocity from kick modes ', &
               s% star_LNA_kick_mode_1, s% star_LNA_kick_mode_2, s% star_LNA_kick_mode_3, &
               ', fractions ', &
               s% star_LNA_kick_fraction_1, s% star_LNA_kick_fraction_2, &
               s% star_LNA_kick_fraction_3, &
               ', period_days = ', period_days, &
               ', surface_kick_km_per_sec = ', s% star_LNA_kick_vsurf_km_per_sec
         end if

         deallocate(kick_eigenvector)
      end subroutine set_star_LNA_initial_velocity


      subroutine add_star_LNA_kick_component( &
            s, map, alphai, vr, mode_indices, mode, fraction, kick_eigenvector, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         real(dp), intent(in) :: alphai(:), vr(:,:), fraction
         integer, intent(in) :: mode_indices(:), mode
         complex(dp), intent(inout) :: kick_eigenvector(:)
         integer, intent(out) :: ierr
         integer :: idx, k, col_v, op_err
         real(dp) :: dR_amp, phase_angle, sign
         complex(dp) :: surf_dlnR, surf_dR, rel_dlnR, dlnR
         complex(dp), allocatable :: eigenvector(:)

         ierr = 0
         if (fraction == 0d0) return

         allocate(eigenvector(map%nvar_total), stat=ierr)
         if (ierr /= 0) return

         idx = mode_indices(mode)
         call load_star_LNA_eigenvector(vr, idx, alphai(idx), eigenvector, op_err)
         if (op_err /= 0) then
            ierr = -1
            deallocate(eigenvector)
            return
         end if

         surf_dlnR = star_LNA_eigen_component(map, eigenvector, 1, lna_var_lnR)
         surf_dR = s% r(1)*surf_dlnR
         if (abs(surf_dR) <= 1d-99) then
            write(*,'(a,i0)') &
               'star_LNA selected kick mode has zero surface displacement amplitude: ', mode
            ierr = -1
            deallocate(eigenvector)
            return
         end if

         do k = 1, map%nz
            dlnR = star_LNA_eigen_component(map, eigenvector, k, lna_var_lnR)
            rel_dlnR = dlnR/surf_dlnR
            dR_amp = abs(s% r(k)*dlnR)/abs(surf_dR)
            phase_angle = atan2(aimag(rel_dlnR), dble(rel_dlnR))
            sign = 1d0
            if (abs(phase_angle) > 0.5d0*pi) sign = -1d0
            col_v = matrix_index(map, k, lna_var_v)
            kick_eigenvector(col_v) = kick_eigenvector(col_v) + fraction*sign*dR_amp
         end do
         deallocate(eigenvector)
      end subroutine add_star_LNA_kick_component


      subroutine set_star_LNA_velocity_from_eigenvector(s, map, eigenvector, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         complex(dp), intent(in) :: eigenvector(:)
         integer, intent(out) :: ierr
         integer :: k
         integer :: max_k, max_v_div_cs_k
         real(dp) :: kick_cms, surf_amp, max_abs_v, max_v_div_cs, v_div_cs
         complex(dp) :: surf_v, phase, dv

         ierr = 0
         kick_cms = 1d5*s% star_LNA_kick_vsurf_km_per_sec
         surf_v = star_LNA_eigen_component(map, eigenvector, 1, lna_var_v)
         surf_amp = abs(surf_v)
         if (surf_amp <= 1d-99) then
            write(*,'(a)') &
               'star_LNA combined kick has zero surface velocity amplitude.'
            ierr = -1
            return
         end if

         phase = surf_v/surf_amp
         s% v_center = 0d0
         max_k = 1
         max_v_div_cs_k = 1
         max_abs_v = 0d0
         max_v_div_cs = 0d0
         do k = 1, map%nz
            dv = star_LNA_eigen_component(map, eigenvector, k, lna_var_v)/phase
            s% v(k) = kick_cms*dble(dv)/surf_amp
            if (is_bad(s% v(k))) then
               write(*,'(a,i0)') 'star_LNA kick produced a bad velocity at k = ', k
               ierr = -1
               return
            end if
            if (abs(s% v(k)) > max_abs_v) then
               max_abs_v = abs(s% v(k))
               max_k = k
            end if
            if (s% csound(k) > 0d0) then
               v_div_cs = abs(s% v(k))/s% csound(k)
               if (v_div_cs > max_v_div_cs) then
                  max_v_div_cs = v_div_cs
                  max_v_div_cs_k = k
               end if
            end if
            if (s% i_v > 0) then
               s% xh(s% i_v, k) = s% v(k)
               s% xh_start(s% i_v, k) = s% v(k)
            end if
            s% v_start(k) = s% v(k)
            s% dxh_v(k) = 0d0
         end do
         write(*,'(a,1pe14.6,a,i0,a,1pe14.6,a,i0)') &
            'star_LNA: kick max |v| = ', max_abs_v/1d5, ' km/s at k = ', max_k, &
            ', max |v|/cs = ', max_v_div_cs, ' at k = ', max_v_div_cs_k
         if (abs(kick_cms) > 0d0 .and. max_abs_v > 10d0*abs(kick_cms)) then
            write(*,'(a,1pe14.6)') &
               'star_LNA warning: interior kick amplitude exceeds 10 times the requested surface kick; ratio = ', &
               max_abs_v/abs(kick_cms)
         end if
      end subroutine set_star_LNA_velocity_from_eigenvector


      logical function selected_for_star_LNA_initial_velocity(s, mode) result(selected)
         type(star_info), pointer :: s
         integer, intent(in) :: mode

         selected = .false.
         if (.not. s% star_LNA_set_initial_velocity) return
         if (s% star_LNA_kick_fraction_1 > 0d0 .and. &
               mode == s% star_LNA_kick_mode_1) selected = .true.
         if (s% star_LNA_kick_mode_2 > 0 .and. &
               s% star_LNA_kick_fraction_2 > 0d0 .and. &
               mode == s% star_LNA_kick_mode_2) selected = .true.
         if (s% star_LNA_kick_mode_3 > 0 .and. &
               s% star_LNA_kick_fraction_3 > 0d0 .and. &
               mode == s% star_LNA_kick_mode_3) selected = .true.
      end function selected_for_star_LNA_initial_velocity


      subroutine load_star_LNA_eigenvector(vr, idx, alphai_idx, eigenvector, ierr)
         real(dp), intent(in) :: vr(:,:), alphai_idx
         integer, intent(in) :: idx
         complex(dp), intent(out) :: eigenvector(:)
         integer, intent(out) :: ierr

         ierr = 0
         if (alphai_idx > 0d0) then
            if (idx >= size(vr, 2)) then
               ierr = -1
               return
            end if
            eigenvector = cmplx(vr(:, idx), vr(:, idx + 1), kind=dp)
         else if (alphai_idx < 0d0) then
            if (idx <= 1) then
               ierr = -1
               return
            end if
            eigenvector = cmplx(vr(:, idx - 1), -vr(:, idx), kind=dp)
         else
            eigenvector = cmplx(vr(:, idx), 0d0, kind=dp)
         end if
      end subroutine load_star_LNA_eigenvector


      subroutine normalize_star_LNA_eigenvector(map, eigenvector)
         type(star_LNA_var_map), intent(in) :: map
         complex(dp), intent(inout) :: eigenvector(:)
         complex(dp) :: norm
         real(dp) :: max_abs
         integer :: i, col_lnR

         col_lnR = matrix_index(map, 1, lna_var_lnR)
         if (col_lnR > 0) then
            norm = eigenvector(col_lnR)
         else
            norm = (0d0, 0d0)
         end if

         if (abs(norm) <= 1d-99) then
            max_abs = 0d0
            do i = 1, size(eigenvector)
               if (abs(eigenvector(i)) > max_abs) then
                  max_abs = abs(eigenvector(i))
                  norm = eigenvector(i)
               end if
            end do
         end if

         if (abs(norm) > 1d-99) eigenvector = eigenvector/norm
      end subroutine normalize_star_LNA_eigenvector


      subroutine write_one_star_LNA_eigenfunction(s, map, mode, eigen_index, &
            sigma_re, sigma_im, eigenvector, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: mode, eigen_index
         real(dp), intent(in) :: sigma_re, sigma_im
         complex(dp), intent(in) :: eigenvector(:)
         integer, intent(out) :: ierr
         integer :: io, k, velocity_var
         character(len=512) :: filename
         character(len=16) :: mode_string
         complex(dp) :: lnd, lnR, v, lnT, L, w, Hp, dLr, dLc, dLt, dL_div_L0

         ierr = 0
         velocity_var = velocity_var_for_star_LNA(map)
         write(mode_string,'(i0)') mode
         call star_LNA_output_filename( &
            s, '_eigenfunction_' // trim(mode_string) // '.data', filename)
         open(newunit=io, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) return

         write(io,'(a)') '# star_LNA eigenfunction'
         write(io,'(a,1x,i0)') '# mode', mode
         write(io,'(a,1x,i0)') '# control_mode_index', mode - 1
         write(io,'(a,1x,l1)') '# selected_for_initial_velocity', &
            selected_for_star_LNA_initial_velocity(s, mode)
         write(io,'(a)') '# eigen_index refers to the reduced dynamic variable eigenproblem'
         write(io,'(a,1x,i0)') '# eigen_index', eigen_index
         write(io,'(a,1x,a)') '# velocity_variable', &
            trim(var_name(velocity_var))
         write(io,'(a,1x,1pe24.16)') '# sigma_real', sigma_re
         write(io,'(a,1x,1pe24.16)') '# sigma_imag', sigma_im
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_vsurf_km_per_sec', &
            s% star_LNA_kick_vsurf_km_per_sec
         write(io,'(a,1x,i0)') '# star_LNA_kick_mode_1', s% star_LNA_kick_mode_1
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_fraction_1', &
            s% star_LNA_kick_fraction_1
         write(io,'(a,1x,i0)') '# star_LNA_kick_mode_2', s% star_LNA_kick_mode_2
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_fraction_2', &
            s% star_LNA_kick_fraction_2
         write(io,'(a,1x,i0)') '# star_LNA_kick_mode_3', s% star_LNA_kick_mode_3
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_fraction_3', &
            s% star_LNA_kick_fraction_3
         write(io,'(a)') &
            '# normalized by surface delta_lnR; uses the largest component if abs(surface delta_lnR) <= 1d-99'
         write(io,'(a)') &
            '# k q m r re_lnd im_lnd re_lnR im_lnR re_v im_v ' // &
            're_lnT im_lnT re_L im_L re_w im_w re_Hp im_Hp ' // &
            're_dLr im_dLr re_dLc im_dLc re_dLt im_dLt ' // &
            'abs_dlnR phase_dlnR abs_dlnT phase_dlnT ' // &
            'abs_dL_div_L0 phase_dL_div_L0 abs_w phase_w'

         do k = 1, map%nz
            lnd = star_LNA_eigen_component(map, eigenvector, k, lna_var_lnd)
            lnR = star_LNA_eigen_component(map, eigenvector, k, lna_var_lnR)
            v = star_LNA_eigen_component(map, eigenvector, k, velocity_var)
            lnT = star_LNA_eigen_component(map, eigenvector, k, lna_var_lnT)
            L = star_LNA_eigen_component(map, eigenvector, k, lna_var_L)
            w = star_LNA_eigen_component(map, eigenvector, k, lna_var_w)
            Hp = star_LNA_eigen_component(map, eigenvector, k, lna_var_Hp)
            call star_LNA_luminosity_perturbations(s, map, k, eigenvector, dLr, dLc, dLt, ierr)
            if (ierr /= 0) then
               close(io)
               return
            end if
            dL_div_L0 = (0d0, 0d0)
            if (abs(s% L(1)) > 0d0) dL_div_L0 = L/s% L(1)
            write(io,'(i8,1x,31(1pe24.16,1x))') &
               k, s% q(k), s% m(k), s% r(k), &
               dble(lnd), aimag(lnd), dble(lnR), aimag(lnR), &
               dble(v), aimag(v), dble(lnT), aimag(lnT), &
               dble(L), aimag(L), dble(w), aimag(w), dble(Hp), aimag(Hp), &
               dble(dLr), aimag(dLr), dble(dLc), aimag(dLc), dble(dLt), aimag(dLt), &
               abs(lnR), phase_for_star_LNA(lnR), abs(lnT), phase_for_star_LNA(lnT), &
               abs(dL_div_L0), phase_for_star_LNA(dL_div_L0), abs(w), phase_for_star_LNA(w)
         end do

         close(io)
      end subroutine write_one_star_LNA_eigenfunction


      ! Diagnostic work output convention:
      !   d(1/rho) = -dlnrho/rho
      !   W_P = -pi*dm*Im(conjg(delta P)*delta(1/rho))
      !   W_L = -pi*dm*Im(conjg(delta lnT)*delta(dL/dm))/omega
      !
      ! RSP LINA writes work divided by modal kinetic energy.  star_LNA does the
      ! same here; the unnormalized eigenfunction is arbitrary, but the normalized
      ! work profile is comparable to RSP LINA_work*.data.
      subroutine write_one_star_LNA_work(s, map, mode, eigen_index, &
            sigma_re, sigma_im, eigenvector, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: mode, eigen_index
         real(dp), intent(in) :: sigma_re, sigma_im
         complex(dp), intent(in) :: eigenvector(:)
         integer, intent(out) :: ierr
         integer :: io, k
         real(dp) :: kinetic_energy, work_norm, &
            pressure_work, turb_pressure_work, eddy_visc_work, &
            rad_lum_work, conv_lum_work, turb_lum_work, total_work, &
            c_pressure_work, c_turb_pressure_work, c_eddy_visc_work, &
            c_rad_lum_work, c_conv_lum_work, c_turb_lum_work, c_total_work
         character(len=512) :: filename
         character(len=16) :: mode_string
         complex(dp) :: dlnrho, dlnT, dspec_vol, dPeos, dPtrb, &
            dLr_00, dLc_00, dLt_00, dLr_p1, dLc_p1, dLt_p1, &
            dLr_dm, dLc_dm, dLt_dm

         ierr = 0
         write(mode_string,'(i0)') mode
         call star_LNA_output_filename(s, '_work_' // trim(mode_string) // '.data', filename)
         open(newunit=io, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) return

         kinetic_energy = kinetic_energy_for_star_LNA_work(s, map, eigenvector)
         work_norm = kinetic_energy
         if (work_norm <= 0d0) work_norm = 1d0

         write(io,'(a)') '# star_LNA work terms'
         write(io,'(a,1x,i0)') '# mode', mode
         write(io,'(a,1x,i0)') '# control_mode_index', mode - 1
         write(io,'(a,1x,l1)') '# selected_for_initial_velocity', &
            selected_for_star_LNA_initial_velocity(s, mode)
         write(io,'(a)') '# eigen_index refers to the reduced dynamic variable eigenproblem'
         write(io,'(a,1x,i0)') '# eigen_index', eigen_index
         write(io,'(a,1x,1pe24.16)') '# sigma_real', sigma_re
         write(io,'(a,1x,1pe24.16)') '# sigma_imag', sigma_im
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_vsurf_km_per_sec', &
            s% star_LNA_kick_vsurf_km_per_sec
         write(io,'(a,1x,i0)') '# star_LNA_kick_mode_1', s% star_LNA_kick_mode_1
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_fraction_1', &
            s% star_LNA_kick_fraction_1
         write(io,'(a,1x,i0)') '# star_LNA_kick_mode_2', s% star_LNA_kick_mode_2
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_fraction_2', &
            s% star_LNA_kick_fraction_2
         write(io,'(a,1x,i0)') '# star_LNA_kick_mode_3', s% star_LNA_kick_mode_3
         write(io,'(a,1x,1pe24.16)') '# star_LNA_kick_fraction_3', &
            s% star_LNA_kick_fraction_3
         write(io,'(a)') '# normalized with the same eigenvector as the eigenfunction file'
         write(io,'(a,1x,1pe24.16)') '# kinetic_energy', kinetic_energy
         write(io,'(a,1x,1pe24.16)') '# work_normalization', work_norm
         write(io,'(a)') '# work columns are divided by work_normalization'
         write(io,'(a)') &
            '# pressure_work and turb_pressure_work use the operator pressure closure'
         write(io,'(a)') &
            '# eddy_visc_work and luminosity work have not been compared with RSP/RSP2 output'
         write(io,'(a)') &
            '# total_work is the sum of the six work columns; it is not constrained by sigma_real'
         write(io,'(a)') &
            '# k logT r_div_R pressure_work turb_pressure_work eddy_visc_work ' // &
            'rad_lum_work conv_lum_work turb_lum_work total_work ' // &
            'c_pressure_work c_turb_pressure_work c_eddy_visc_work ' // &
            'c_rad_lum_work c_conv_lum_work c_turb_lum_work c_total_work'

         c_pressure_work = 0d0
         c_turb_pressure_work = 0d0
         c_eddy_visc_work = 0d0
         c_rad_lum_work = 0d0
         c_conv_lum_work = 0d0
         c_turb_lum_work = 0d0
         c_total_work = 0d0

         do k = 1, map%nz
            dlnrho = star_LNA_eigen_component(map, eigenvector, k, lna_var_lnd)
            dlnT = star_LNA_eigen_component(map, eigenvector, k, lna_var_lnT)
            dspec_vol = -dlnrho/s% rho(k)

            call star_LNA_pressure_perturbations(s, map, k, eigenvector, dPeos, dPtrb, ierr)
            if (ierr /= 0) then
               close(io)
               return
            end if

            call star_LNA_luminosity_perturbations(s, map, k, eigenvector, &
               dLr_00, dLc_00, dLt_00, ierr)
            if (ierr /= 0) then
               close(io)
               return
            end if
            if (k < map%nz) then
               call star_LNA_luminosity_perturbations(s, map, k + 1, eigenvector, &
                  dLr_p1, dLc_p1, dLt_p1, ierr)
               if (ierr /= 0) then
                  close(io)
                  return
               end if
            else
               dLr_p1 = (0d0, 0d0)
               dLc_p1 = (0d0, 0d0)
               dLt_p1 = (0d0, 0d0)
            end if

            dLr_dm = (dLr_00 - dLr_p1)/s% dm(k)
            dLc_dm = (dLc_00 - dLc_p1)/s% dm(k)
            dLt_dm = (dLt_00 - dLt_p1)/s% dm(k)

            pressure_work = -pi*s% dm(k)*aimag(conjg(dPeos)*dspec_vol)
            turb_pressure_work = -pi*s% dm(k)*aimag(conjg(dPtrb)*dspec_vol)
            eddy_visc_work = eddy_viscous_work_for_star_LNA( &
               s, map, k, sigma_im, eigenvector, ierr)
            if (ierr /= 0) then
               close(io)
               return
            end if
            rad_lum_work = luminosity_work_for_star_LNA(s, k, sigma_im, dlnT, dLr_dm)
            conv_lum_work = luminosity_work_for_star_LNA(s, k, sigma_im, dlnT, dLc_dm)
            turb_lum_work = luminosity_work_for_star_LNA(s, k, sigma_im, dlnT, dLt_dm)
            pressure_work = pressure_work/work_norm
            turb_pressure_work = turb_pressure_work/work_norm
            eddy_visc_work = eddy_visc_work/work_norm
            rad_lum_work = rad_lum_work/work_norm
            conv_lum_work = conv_lum_work/work_norm
            turb_lum_work = turb_lum_work/work_norm
            total_work = pressure_work + turb_pressure_work + eddy_visc_work + &
               rad_lum_work + conv_lum_work + turb_lum_work

            c_pressure_work = c_pressure_work + pressure_work
            c_turb_pressure_work = c_turb_pressure_work + turb_pressure_work
            c_eddy_visc_work = c_eddy_visc_work + eddy_visc_work
            c_rad_lum_work = c_rad_lum_work + rad_lum_work
            c_conv_lum_work = c_conv_lum_work + conv_lum_work
            c_turb_lum_work = c_turb_lum_work + turb_lum_work
            c_total_work = c_total_work + total_work

            write(io,'(i8,1x,16(1pe24.16,1x))') &
               k, s% lnT(k)/ln10, s% r(k)/s% r(1), &
               pressure_work, turb_pressure_work, eddy_visc_work, &
               rad_lum_work, conv_lum_work, turb_lum_work, total_work, &
               c_pressure_work, c_turb_pressure_work, c_eddy_visc_work, &
               c_rad_lum_work, c_conv_lum_work, c_turb_lum_work, c_total_work
         end do

         close(io)
      end subroutine write_one_star_LNA_work


      real(dp) function kinetic_energy_for_star_LNA_work(s, map, eigenvector) &
            result(kinetic_energy)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         complex(dp), intent(in) :: eigenvector(:)
         integer :: k, velocity_var
         real(dp) :: dm_velocity
         complex(dp) :: dv

         kinetic_energy = 0d0
         velocity_var = velocity_var_for_star_LNA(map)
         do k = 1, map%nz
            dv = star_LNA_eigen_component(map, eigenvector, k, velocity_var)
            if (velocity_var == lna_var_u) then
               dm_velocity = s% dm(k)
            else
               dm_velocity = s% dm_bar(k)
            end if
            kinetic_energy = kinetic_energy + 0.5d0*dm_velocity*pow2(abs(dv))
         end do
      end function kinetic_energy_for_star_LNA_work


      subroutine star_LNA_pressure_perturbations(s, map, k, eigenvector, dPeos, dPtrb, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: k
         complex(dp), intent(in) :: eigenvector(:)
         complex(dp), intent(out) :: dPeos, dPtrb
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Peos_ad, Ptrb_ad

         call pressure_components_for_star_LNA(s, k, Peos_ad, Ptrb_ad, ierr)
         if (ierr /= 0) return

         call star_LNA_perturbation_from_ad(map, k, Peos_ad, eigenvector, dPeos, ierr)
         if (ierr /= 0) return
         call star_LNA_perturbation_from_ad(map, k, Ptrb_ad, eigenvector, dPtrb, ierr)
      end subroutine star_LNA_pressure_perturbations


      real(dp) function eddy_viscous_work_for_star_LNA( &
            s, map, k, sigma_im, eigenvector, ierr) result(eddy_work)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: k
         real(dp), intent(in) :: sigma_im
         complex(dp), intent(in) :: eigenvector(:)
         integer, intent(out) :: ierr
         real(dp) :: chi_coeff
         complex(dp) :: d_v_div_r

         ierr = 0
         eddy_work = 0d0
         if (.not. s% star_LNA_perturb_eddy_viscosity) return
         if (abs(sigma_im) <= 1d-99) return

         chi_coeff = chi_coefficient_for_star_LNA(s, k, ierr)
         if (ierr /= 0) return
         if (chi_coeff <= 0d0) return

         call delta_d_v_div_r_for_star_LNA(s, map, k, eigenvector, d_v_div_r, ierr)
         if (ierr /= 0) return
         eddy_work = -4d0*pi*pi*chi_coeff*abs(d_v_div_r)*abs(d_v_div_r)/abs(sigma_im)
      end function eddy_viscous_work_for_star_LNA


      logical function frozen_flux_lna_active(s) result(active)
         type(star_info), pointer :: s

         active = .not. s% RSP2_flag .and. frozen_flux_lna_selected(s)
      end function frozen_flux_lna_active


      logical function frozen_flux_lna_selected(s) result(selected)
         type(star_info), pointer :: s

         selected = trim(s% star_LNA_convection_treatment) == 'frozen_flux'
      end function frozen_flux_lna_selected


      subroutine delta_d_v_div_r_for_star_LNA( &
            s, map, k, eigenvector, d_v_div_r, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: k
         complex(dp), intent(in) :: eigenvector(:)
         complex(dp), intent(out) :: d_v_div_r
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: v00_ad, vp1_ad
         complex(dp) :: v00, vp1
         real(dp) :: rp1

         ierr = 0
         v00_ad = wrap_v_00(s, k)
         vp1_ad = wrap_v_p1(s, k)
         call star_LNA_perturbation_from_ad(map, k, v00_ad, eigenvector, v00, ierr)
         if (ierr /= 0) return
         call star_LNA_perturbation_from_ad(map, k, vp1_ad, eigenvector, vp1, ierr)
         if (ierr /= 0) return
         if (k < s% nz) then
            rp1 = s% r(k + 1)
         else
            rp1 = s% r_center
         end if
         if (rp1 == 0d0) rp1 = 1d0

         d_v_div_r = v00/s% r(k) - vp1/rp1
      end subroutine delta_d_v_div_r_for_star_LNA


      real(dp) function luminosity_work_for_star_LNA(s, k, sigma_im, dlnT, dL_dm) &
            result(lum_work)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: sigma_im
         complex(dp), intent(in) :: dlnT, dL_dm

         if (abs(sigma_im) <= 1d-99) then
            lum_work = 0d0
         else
            lum_work = -pi*s% dm(k)*aimag(conjg(dlnT)*dL_dm)/sigma_im
         end if
      end function luminosity_work_for_star_LNA


      subroutine star_LNA_luminosity_perturbations(s, map, k, eigenvector, &
            dLr, dLc, dLt, ierr)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: k
         complex(dp), intent(in) :: eigenvector(:)
         complex(dp), intent(out) :: dLr, dLc, dLt
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lr_ad, Lc_ad, Lt_ad

         ierr = 0
         if (s% RSP2_flag) then
            call rsp2_luminosity_terms_for_star_LNA(s, k, Lr_ad, Lc_ad, Lt_ad)
         else
            Lc_ad = star_LNA_L_conv_ad(s, k)
            Lr_ad = wrap_L_00(s, k) - Lc_ad
            Lt_ad = 0d0
         end if

         call star_LNA_perturbation_from_ad(map, k, Lr_ad, eigenvector, dLr, ierr)
         if (ierr /= 0) return
         call star_LNA_perturbation_from_ad(map, k, Lc_ad, eigenvector, dLc, ierr)
         if (ierr /= 0) return
         call star_LNA_perturbation_from_ad(map, k, Lt_ad, eigenvector, dLt, ierr)
      end subroutine star_LNA_luminosity_perturbations


      subroutine star_LNA_perturbation_from_ad(map, k, ad, eigenvector, perturbation, ierr)
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: ad
         complex(dp), intent(in) :: eigenvector(:)
         complex(dp), intent(out) :: perturbation
         integer, intent(out) :: ierr
         integer :: iad, kk, var_id, col
         real(dp) :: coeff

         ierr = 0
         perturbation = (0d0, 0d0)
         do iad = 1, auto_diff_star_num_vars
            coeff = ad%d1Array(iad)
            if (coeff == 0d0) cycle
            call ad_index_to_star_LNA_var(map, k, iad, kk, var_id)
            if (var_id == 0) then
               ierr = -1
               return
            end if
            if (kk < 1 .or. kk > map%nz) cycle
            col = matrix_index(map, kk, var_id)
            if (col <= 0) then
               ierr = -1
               return
            end if
            perturbation = perturbation + coeff*eigenvector(col)
         end do
      end subroutine star_LNA_perturbation_from_ad


      complex(dp) function star_LNA_eigen_component(map, eigenvector, k, var_id) result(component)
         type(star_LNA_var_map), intent(in) :: map
         complex(dp), intent(in) :: eigenvector(:)
         integer, intent(in) :: k, var_id
         integer :: col

         component = (0d0, 0d0)
         col = matrix_index(map, k, var_id)
         if (col > 0) component = eigenvector(col)
      end function star_LNA_eigen_component


      subroutine finite_sigma_from_star_LNA_eigenvalue(alphar, alphai, beta, &
            sigma_re, sigma_im, ierr)
         real(dp), intent(in) :: alphar, alphai, beta
         real(dp), intent(out) :: sigma_re, sigma_im
         integer, intent(out) :: ierr

         ierr = 0
         sigma_re = 0d0
         sigma_im = 0d0
         if (abs(beta) <= 1d-99) then
            ierr = -1
            return
         end if
         sigma_re = alphar/beta
         sigma_im = alphai/beta
      end subroutine finite_sigma_from_star_LNA_eigenvalue


      real(dp) function rsp_logKE_per_cycle_for_star_LNA(sigma_re, sigma_im) &
            result(logKE_per_cycle)
         real(dp), intent(in) :: sigma_re, sigma_im

         if (abs(sigma_im) <= 1d-99) then
            logKE_per_cycle = 0d0
         else
            logKE_per_cycle = 4d0*pi*sigma_re/sigma_im
         end if
      end function rsp_logKE_per_cycle_for_star_LNA


      real(dp) function ke_fractional_growth_per_period_for_star_LNA(sigma_re, sigma_im) &
            result(growth_per_period)
         real(dp), intent(in) :: sigma_re, sigma_im
         real(dp) :: logKE_per_cycle

         logKE_per_cycle = rsp_logKE_per_cycle_for_star_LNA(sigma_re, sigma_im)
         growth_per_period = safe_expm1_for_star_LNA(logKE_per_cycle)
      end function ke_fractional_growth_per_period_for_star_LNA


      real(dp) function rsp_grekm_growth_per_period_for_star_LNA(sigma_re, sigma_im) &
            result(growth_per_period)
         real(dp), intent(in) :: sigma_re, sigma_im
         real(dp) :: logKE_per_cycle

         logKE_per_cycle = rsp_logKE_per_cycle_for_star_LNA(sigma_re, sigma_im)
         growth_per_period = 2d0*tanh(0.5d0*logKE_per_cycle)
      end function rsp_grekm_growth_per_period_for_star_LNA


      real(dp) function amp_fractional_growth_per_period_for_star_LNA(sigma_re, sigma_im) &
            result(growth_per_period)
         real(dp), intent(in) :: sigma_re, sigma_im
         real(dp) :: logKE_per_cycle

         logKE_per_cycle = rsp_logKE_per_cycle_for_star_LNA(sigma_re, sigma_im)
         growth_per_period = safe_expm1_for_star_LNA(0.5d0*logKE_per_cycle)
      end function amp_fractional_growth_per_period_for_star_LNA


      real(dp) function safe_expm1_for_star_LNA(x) result(expm1_value)
         real(dp), intent(in) :: x

         if (x > log(huge(1d0))) then
            expm1_value = huge(1d0)
         else if (x < log(tiny(1d0))) then
            expm1_value = -1d0
         else if (abs(x) < 1d-5) then
            expm1_value = x*(1d0 + x*(0.5d0 + x*(1d0/6d0 + x*(1d0/24d0))))
         else
            expm1_value = exp(x) - 1d0
         end if
      end function safe_expm1_for_star_LNA


      real(dp) function period_sec_for_star_LNA(sigma_im) result(period_sec)
         real(dp), intent(in) :: sigma_im

         if (abs(sigma_im) <= 1d-99) then
            period_sec = huge(1d0)
         else
            period_sec = 2d0*pi/abs(sigma_im)
         end if
      end function period_sec_for_star_LNA


      real(dp) function period_days_for_star_LNA(sigma_im) result(period_days)
         real(dp), intent(in) :: sigma_im

         period_days = period_sec_for_star_LNA(sigma_im)/86400d0
      end function period_days_for_star_LNA


      real(dp) function pulsation_constant_Q_days_for_star_LNA(s, period_days) result(Q_days)
         type(star_info), pointer :: s
         real(dp), intent(in) :: period_days
         real(dp) :: radius_ratio

         Q_days = 0d0
         if (s% m(1) <= 0d0 .or. s% r(1) <= 0d0) return

         radius_ratio = Rsun/s% r(1)
         Q_days = period_days*sqrt((s% m(1)/Msun)*pow2(radius_ratio)*radius_ratio)
      end function pulsation_constant_Q_days_for_star_LNA


      real(dp) function phase_for_star_LNA(z) result(phase)
         complex(dp), intent(in) :: z

         if (abs(z) <= 0d0) then
            phase = 0d0
         else
            phase = atan2(aimag(z), dble(z))
         end if
      end function phase_for_star_LNA


      real(dp) function frequency_uHz_for_star_LNA(sigma_im) result(frequency_uHz)
         real(dp), intent(in) :: sigma_im

         frequency_uHz = abs(sigma_im)/(2d0*pi)*1d6
      end function frequency_uHz_for_star_LNA

      ! Matrix and row audit diagnostics.
      subroutine write_star_LNA_matrix_summary(s, problem, ierr)
         type(star_info), pointer :: s
         type(star_LNA_problem), intent(in) :: problem
         integer, intent(out) :: ierr
         integer :: io
         real(dp) :: diag, min_row_norm, max_row_norm, min_col_norm, max_col_norm
         character(len=512) :: filename

         ierr = 0
         call star_LNA_output_filename(s, '_matrix_summary.data', filename)
         open(newunit=io, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) return

         write(io,'(a)') '# star_LNA matrix summary'
         write(io,'(a,1x,i0)') 'nz', problem% map% nz
         write(io,'(a,1x,i0)') 'nvar_per_zone', problem% map% nvar_per_zone
         write(io,'(a,1x,i0)') 'nvar_total', problem% map% nvar_total
         write(io,'(a,1x,i0)') 'nnz_A', count(problem% mtx% A /= 0d0)
         write(io,'(a,1x,i0)') 'nnz_B', count(problem% mtx% B /= 0d0)
         call write_star_LNA_equation_counts(io, problem)
         call matrix_norm_range_for_star_LNA( &
            problem% mtx, min_row_norm, max_row_norm, min_col_norm, max_col_norm)
         write(io,'(a,1x,1pe24.16)') 'min_nonzero_raw_row_norm_AB', min_row_norm
         write(io,'(a,1x,1pe24.16)') 'max_raw_row_norm_AB', max_row_norm
         write(io,'(a,1x,1pe24.16)') 'min_nonzero_raw_col_norm_AB', min_col_norm
         write(io,'(a,1x,1pe24.16)') 'max_raw_col_norm_AB', max_col_norm
         write(io,'(a,1x,1pe24.16)') 'max_abs_v_div_csound', max_star_LNA_v_div_csound(s)
         write(io,'(a,1x,1pe24.16)') &
            'max_abs_density_mass_resid', max_abs_density_mass_resid_for_star_LNA(s)
         call max_abs_momentum_row_resid_for_star_LNA(s, diag, ierr)
         if (ierr /= 0) then
            close(io)
            return
         end if
         write(io,'(a,1x,1pe24.16)') 'max_abs_momentum_row_resid', diag
         call max_abs_energy_row_resid_for_star_LNA(s, diag, ierr)
         if (ierr /= 0) then
            close(io)
            return
         end if
         write(io,'(a,1x,1pe24.16)') 'max_abs_energy_row_resid', diag
         call max_abs_luminosity_row_resid_for_star_LNA(s, diag, ierr)
         if (ierr /= 0) then
            close(io)
            return
         end if
         write(io,'(a,1x,1pe24.16)') 'max_abs_luminosity_row_resid', diag
         call max_abs_rsp2_w_row_resid_for_star_LNA(s, diag, ierr)
         if (ierr /= 0) then
            close(io)
            return
         end if
         write(io,'(a,1x,1pe24.16)') 'max_abs_rsp2_w_row_resid', diag
         call max_abs_tdc_w_row_resid_for_star_LNA(s, diag, ierr)
         if (ierr /= 0) then
            close(io)
            return
         end if
         write(io,'(a,1x,1pe24.16)') 'max_abs_tdc_w_row_resid', diag
         call max_abs_rsp2_Hp_row_resid_for_star_LNA(s, diag, ierr)
         if (ierr /= 0) then
            close(io)
            return
         end if
         write(io,'(a,1x,1pe24.16)') 'max_abs_rsp2_Hp_row_resid', diag
         write(io,'(a)') &
            '# dense solve row/column scales A and B before algebraic elimination'
         write(io,'(a)') '# eigenvalues still solve the original A*x = sigma*B*x problem'
         close(io)

         call write_star_LNA_row_structure(s, problem, ierr)
         if (ierr /= 0) return
         if (s% RSP2_flag) call write_star_LNA_rsp2_term_audit(s, ierr)
         if (ierr /= 0) return
         if (tdc_lna_active(s)) call write_star_LNA_tdc_face_audit(s, ierr)
      end subroutine write_star_LNA_matrix_summary


      subroutine write_star_LNA_equation_counts(io, problem)
         integer, intent(in) :: io
         type(star_LNA_problem), intent(in) :: problem
         integer :: eq_id, num_rows

         write(io,'(a)') '# equation row counts'
         do eq_id = 1, num_lna_equations
            num_rows = count(problem% eq_id == eq_id)
            if (num_rows == 0) cycle
            write(io,'(2a,1x,i0)') 'equation_count ', &
               trim(equation_name_for_star_LNA(eq_id)), num_rows
         end do
      end subroutine write_star_LNA_equation_counts


      subroutine write_star_LNA_row_structure(s, problem, ierr)
         type(star_info), pointer :: s
         type(star_LNA_problem), intent(in) :: problem
         integer, intent(out) :: ierr
         integer :: io, row, k, slot, nnz_A, nnz_B, dom_A_col, dom_B_col, &
            dom_A_k, dom_A_slot, dom_B_k, dom_B_slot
         real(dp) :: max_abs_A, max_abs_B
         character(len=512) :: filename

         ierr = 0
         call star_LNA_output_filename(s, '_row_structure.data', filename)
         open(newunit=io, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) return

         write(io,'(a)') '# star_LNA row structure'
         write(io,'(a)') '# row and column indices start at one and refer to the full matrix'
         write(io,'(a)') &
            '# row zone row_var equation nnz_A nnz_B max_abs_A max_abs_B ' // &
            'dom_A_col dom_A_zone dom_A_var dom_B_col dom_B_zone dom_B_var'

         do row = 1, problem% mtx% n
            call matrix_index_to_zone_slot(problem% map, row, k, slot)
            call dominant_row_entry(problem% mtx% A(row, :), nnz_A, max_abs_A, dom_A_col)
            call dominant_row_entry(problem% mtx% B(row, :), nnz_B, max_abs_B, dom_B_col)
            call matrix_index_to_zone_slot(problem% map, dom_A_col, dom_A_k, dom_A_slot)
            call matrix_index_to_zone_slot(problem% map, dom_B_col, dom_B_k, dom_B_slot)
            write(io,'(2(i0,1x),2(a,1x),2(i0,1x),2(1pe24.16,1x),2(i0,1x,i0,1x,a,1x))') &
               row, k, trim(var_name_from_slot(problem% map, slot)), &
               trim(equation_name_for_star_LNA(problem% eq_id(row))), &
               nnz_A, nnz_B, max_abs_A, max_abs_B, &
               dom_A_col, dom_A_k, trim(var_name_from_slot(problem% map, dom_A_slot)), &
               dom_B_col, dom_B_k, trim(var_name_from_slot(problem% map, dom_B_slot))
         end do

         close(io)
      end subroutine write_star_LNA_row_structure


      subroutine write_star_LNA_rsp2_term_audit(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: io, k
         character(len=512) :: filename
         type(auto_diff_real_star_order1) :: PII_ad, Lc_ad, Lt_ad, source_ad, &
            damping_ad, rad_damping_ad, Ptrb_div_etrb_ad, dwork_dm_ad, &
            dLt_dm_ad, rhs_ad, inertia_ad

         ierr = 0
         call star_LNA_output_filename(s, '_rsp2_term_audit.data', filename)
         open(newunit=io, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) return

         write(io,'(a)') '# star_LNA RSP2 background term audit'
         write(io,'(a)') '# k=1 is the surface; PII is face centered at face k'
         write(io,'(a)') &
            '# k forced_non_turbulent PII_face L_conv L_turb source dissipation ' // &
            'radiative_damping Ptrb_div_etrb Ptrb_dVdt_div_dm dLt_div_dm rhs inertia'

         do k = 1, s% nz
            call rsp2_terms_for_star_LNA_audit(s, k, PII_ad, Lc_ad, Lt_ad, &
               source_ad, damping_ad, rad_damping_ad, Ptrb_div_etrb_ad, &
               dwork_dm_ad, dLt_dm_ad, rhs_ad, inertia_ad, ierr)
            if (ierr /= 0) then
               close(io)
               return
            end if
            write(io,'(i8,1x,l1,1x,11(1pe24.16,1x))') k, &
               rsp2_forces_non_turbulent_cell(s, k), PII_ad% val, Lc_ad% val, &
               Lt_ad% val, source_ad% val, damping_ad% val, rad_damping_ad% val, &
               Ptrb_div_etrb_ad% val, dwork_dm_ad% val, dLt_dm_ad% val, &
               rhs_ad% val, inertia_ad% val
         end do

         close(io)
      end subroutine write_star_LNA_rsp2_term_audit


      subroutine write_star_LNA_tdc_face_audit(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: io, k
         character(len=512) :: filename
         type(auto_diff_real_star_order1) :: T_face_ad, rho_face_ad, &
            Peos_face_ad, energy_face_ad, Cp_face_ad, chiRho_face_ad, &
            chiT_face_ad, grada_face_ad, opacity_face_ad, scale_height_face_ad, &
            gradr_face_ad, luminosity_resid_ad, velocity_rhs_ad, &
            velocity_inertia_ad, L_rad_ad, L_conv_ad, stored_T_face_ad, &
            stored_rho_face_ad, stored_Peos_face_ad, stored_Cp_face_ad, &
            stored_opacity_face_ad, stored_scale_height_face_ad, &
            stored_grada_face_ad, stored_gradr_face_ad

         ierr = 0
         call star_LNA_output_filename(s, '_tdc_face_audit.data', filename)
         open(newunit=io, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) return

         write(io,'(a)') '# star_LNA TDC face state audit'
         write(io,'(a)') '# k=1 is the surface; all thermodynamic quantities are at face k'
         write(io,'(a)') &
            '# k use_face_reconstruction zero_w active_T active_rho active_P active_Cp ' // &
            'active_kap active_Hp active_grada active_gradr stored_T stored_rho stored_P ' // &
            'stored_Cp stored_kap stored_Hp stored_grada stored_gradr mlt_vc gradT ' // &
            'background_Lconv closure_Lrad closure_Lconv luminosity_resid velocity_rhs inertia'

         do k = 1, s% nz
            call tdc_face_state_for_star_LNA(s, k, &
               T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
               chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
               scale_height_face_ad, gradr_face_ad, ierr)
            if (ierr /= 0) then
               close(io)
               return
            end if
            stored_T_face_ad = get_T_face(s, k)
            stored_rho_face_ad = get_rho_face(s, k)
            stored_Peos_face_ad = get_Peos_face(s, k)
            stored_Cp_face_ad = get_Cp_face(s, k)
            stored_opacity_face_ad = get_kap_face(s, k)
            stored_scale_height_face_ad = get_scale_height_face(s, k)
            stored_grada_face_ad = get_grada_face(s, k)
            stored_gradr_face_ad = get_gradr_face(s, k)

            luminosity_resid_ad = 0d0
            velocity_rhs_ad = 0d0
            velocity_inertia_ad = 0d0
            L_rad_ad = 0d0
            L_conv_ad = 0d0
            if (k > 1) then
               call tdc_relation_for_star_LNA(s, k, luminosity_resid_ad, &
                  velocity_rhs_ad, velocity_inertia_ad, L_rad_ad, L_conv_ad, ierr)
               if (ierr /= 0) then
                  close(io)
                  return
               end if
            end if

            write(io,'(i8,1x,2(l1,1x),24(1pe24.16,1x))') k, &
               s% use_face_reconstruction, tdc_zero_w_for_star_LNA(s, k), &
               T_face_ad% val, rho_face_ad% val, Peos_face_ad% val, Cp_face_ad% val, &
               opacity_face_ad% val, scale_height_face_ad% val, grada_face_ad% val, &
               gradr_face_ad% val, stored_T_face_ad% val, stored_rho_face_ad% val, &
               stored_Peos_face_ad% val, stored_Cp_face_ad% val, stored_opacity_face_ad% val, &
               stored_scale_height_face_ad% val, stored_grada_face_ad% val, &
               stored_gradr_face_ad% val, s% mlt_vc(k), s% gradT(k), s% L_conv(k), &
               L_rad_ad% val, L_conv_ad% val, luminosity_resid_ad% val, &
               velocity_rhs_ad% val, velocity_inertia_ad% val
         end do

         close(io)
      end subroutine write_star_LNA_tdc_face_audit


      subroutine dominant_row_entry(row_values, nnz, max_abs_value, dominant_col)
         real(dp), intent(in) :: row_values(:)
         integer, intent(out) :: nnz, dominant_col
         real(dp), intent(out) :: max_abs_value
         integer :: j
         real(dp) :: abs_value

         nnz = 0
         dominant_col = 0
         max_abs_value = 0d0
         do j = 1, size(row_values)
            if (row_values(j) == 0d0) cycle
            nnz = nnz + 1
            abs_value = abs(row_values(j))
            if (abs_value > max_abs_value) then
               max_abs_value = abs_value
               dominant_col = j
            end if
         end do
      end subroutine dominant_row_entry


      subroutine matrix_index_to_zone_slot(map, idx, k, slot)
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: idx
         integer, intent(out) :: k, slot

         if (idx < 1 .or. idx > map%nvar_total) then
            k = 0
            slot = 0
            return
         end if
         k = (idx - 1)/map%nvar_per_zone + 1
         slot = mod(idx - 1, map%nvar_per_zone) + 1
      end subroutine matrix_index_to_zone_slot


      function var_name_from_slot(map, slot) result(name)
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: slot
         character(len=8) :: name

         if (slot < 1 .or. slot > map%nvar_per_zone) then
            name = 'none'
         else
            name = var_name(map%var_id(slot))
         end if
      end function var_name_from_slot


      subroutine matrix_norm_range_for_star_LNA( &
            mtx, min_row_norm, max_row_norm, min_col_norm, max_col_norm)
         type(star_LNA_matrix), intent(in) :: mtx
         real(dp), intent(out) :: min_row_norm, max_row_norm, min_col_norm, max_col_norm
         integer :: i, j
         real(dp) :: norm

         min_row_norm = huge(1d0)
         max_row_norm = 0d0
         do i = 1, mtx%n
            norm = max(maxval(abs(mtx%A(i, :))), maxval(abs(mtx%B(i, :))))
            if (norm > 0d0) min_row_norm = min(min_row_norm, norm)
            max_row_norm = max(max_row_norm, norm)
         end do
         if (min_row_norm == huge(1d0)) min_row_norm = 0d0

         min_col_norm = huge(1d0)
         max_col_norm = 0d0
         do j = 1, mtx%n
            norm = max(maxval(abs(mtx%A(:, j))), maxval(abs(mtx%B(:, j))))
            if (norm > 0d0) min_col_norm = min(min_col_norm, norm)
            max_col_norm = max(max_col_norm, norm)
         end do
         if (min_col_norm == huge(1d0)) min_col_norm = 0d0
      end subroutine matrix_norm_range_for_star_LNA


      real(dp) function max_abs_density_mass_resid_for_star_LNA(s) result(max_resid)
         type(star_info), pointer :: s
         integer :: k, ierr
         real(dp) :: ratio
         type(auto_diff_real_star_order1) :: cell_volume_ad

         max_resid = 0d0
         do k = 1, s% nz
            call cell_volume_for_star_LNA(s, k, cell_volume_ad, ierr)
            if (ierr /= 0 .or. cell_volume_ad%val <= 0d0 .or. s% dm(k) <= 0d0) then
               max_resid = huge(1d0)
               return
            end if
            ratio = s% rho(k)*cell_volume_ad%val/s% dm(k)
            if (ratio <= 0d0) then
               max_resid = huge(1d0)
               return
            end if
            max_resid = max(max_resid, abs(log(ratio)))
         end do
      end function max_abs_density_mass_resid_for_star_LNA


      subroutine max_abs_momentum_row_resid_for_star_LNA(s, max_resid, ierr)
         type(star_info), pointer :: s
         real(dp), intent(out) :: max_resid
         integer, intent(out) :: ierr
         integer :: k
         type(auto_diff_real_star_order1) :: resid_ad, P_bc_ad, lnP_bc_ad

         ierr = 0
         max_resid = 0d0
         do k = 1, s% nz
            if (k == 1) then
               if (s% use_fixed_vsurf_outer_BC) then
                  if (s% u_flag) then
                     resid_ad = s% u(1)
                  else
                     resid_ad = s% v(1)
                  end if
               else
                  call surface_P_bc_for_star_LNA(s, P_bc_ad, lnP_bc_ad, ierr)
                  if (ierr /= 0) return
                  if (use_surface_momentum_row_for_star_LNA(s)) then
                     call surface_momentum_rhs_for_star_LNA(s, P_bc_ad, resid_ad, ierr)
                     if (ierr /= 0) return
                  else
                     resid_ad = lnP_bc_ad - wrap_lnPeos_00(s, 1)
                  end if
               end if
            else
               call momentum_rhs_for_star_LNA(s, k, resid_ad, ierr)
               if (ierr /= 0) return
            end if
            max_resid = max(max_resid, abs(resid_ad%val))
         end do
      end subroutine max_abs_momentum_row_resid_for_star_LNA


      subroutine max_abs_energy_row_resid_for_star_LNA(s, max_resid, ierr)
         type(star_info), pointer :: s
         real(dp), intent(out) :: max_resid
         integer, intent(out) :: ierr
         integer :: k
         type(auto_diff_real_star_order1) :: resid_ad

         ierr = 0
         max_resid = 0d0
         do k = 1, s% nz
            call energy_rhs_for_star_LNA(s, k, resid_ad, ierr)
            if (ierr /= 0) return
            max_resid = max(max_resid, abs(resid_ad%val))
         end do
      end subroutine max_abs_energy_row_resid_for_star_LNA


      subroutine max_abs_luminosity_row_resid_for_star_LNA(s, max_resid, ierr)
         type(star_info), pointer :: s
         real(dp), intent(out) :: max_resid
         integer, intent(out) :: ierr
         integer :: k
         type(auto_diff_real_star_order1) :: resid_ad

         ierr = 0
         max_resid = 0d0
         do k = 1, s% nz
            if (use_rsp_lsurf_row_for_star_LNA(s, k)) then
               call rsp_lsurf_resid_for_star_LNA(s, resid_ad)
            else if (s% RSP2_flag) then
               call rsp2_luminosity_resid_for_star_LNA(s, k, resid_ad)
               resid_ad = resid_ad/max(1d0, maxval(abs(s% L(1:s% nz))))
            else if (tdc_lna_active(s) .and. k > 1) then
               call tdc_luminosity_resid_for_star_LNA(s, k, resid_ad, ierr)
               if (ierr /= 0) return
               resid_ad = resid_ad/max(1d0, maxval(abs(s% L(1:s% nz))))
            else if (frozen_flux_lna_active(s) .and. k > 1) then
               call frozen_flux_luminosity_resid_for_star_LNA(s, k, resid_ad, ierr)
               if (ierr /= 0) return
               resid_ad = resid_ad/max(1d0, maxval(abs(s% L(1:s% nz))))
            else if (k == 1) then
               call surface_temperature_resid_for_star_LNA(s, resid_ad, ierr)
               if (ierr /= 0) return
            else
               call temperature_gradient_resid_for_star_LNA(s, k, resid_ad, ierr)
               if (ierr /= 0) return
            end if
            max_resid = max(max_resid, abs(resid_ad%val))
         end do
      end subroutine max_abs_luminosity_row_resid_for_star_LNA


      subroutine max_abs_rsp2_w_row_resid_for_star_LNA(s, max_resid, ierr)
         type(star_info), pointer :: s
         real(dp), intent(out) :: max_resid
         integer, intent(out) :: ierr
         integer :: k
         type(auto_diff_real_star_order1) :: resid_ad

         ierr = 0
         max_resid = 0d0
         if (.not. s% RSP2_flag) return
         do k = 1, s% nz
            if (rsp2_forces_non_turbulent_cell(s, k)) then
               resid_ad = s% w(k)/max(1d0, s% csound(k))
            else
               call rsp2_turbulent_energy_rhs_for_star_LNA(s, k, resid_ad, ierr)
               if (ierr /= 0) return
            end if
            max_resid = max(max_resid, abs(resid_ad%val))
         end do
      end subroutine max_abs_rsp2_w_row_resid_for_star_LNA


      subroutine max_abs_tdc_w_row_resid_for_star_LNA(s, max_resid, ierr)
         type(star_info), pointer :: s
         real(dp), intent(out) :: max_resid
         integer, intent(out) :: ierr
         integer :: k
         type(auto_diff_real_star_order1) :: luminosity_resid_ad, resid_ad, &
            inertia_ad, L_rad_ad, L_conv_ad

         ierr = 0
         max_resid = 0d0
         if (.not. tdc_lna_active(s)) return
         do k = 1, s% nz
            if (tdc_zero_w_for_star_LNA(s, k)) then
               resid_ad = s% mlt_vc(k)/sqrt_2_div_3
            else
               call tdc_relation_for_star_LNA(s, k, luminosity_resid_ad, resid_ad, &
                  inertia_ad, L_rad_ad, L_conv_ad, ierr)
               if (ierr /= 0) return
            end if
            max_resid = max(max_resid, abs(resid_ad%val))
         end do
      end subroutine max_abs_tdc_w_row_resid_for_star_LNA


      subroutine max_abs_rsp2_Hp_row_resid_for_star_LNA(s, max_resid, ierr)
         type(star_info), pointer :: s
         real(dp), intent(out) :: max_resid
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: scale
         type(auto_diff_real_star_order1) :: Hp_expected_ad, resid_ad

         ierr = 0
         max_resid = 0d0
         if (.not. s% RSP2_flag) return
         do k = 1, s% nz
            Hp_expected_ad = Hp_face_for_rsp2_eqn(s, k, ierr)
            if (ierr /= 0) return
            scale = 1d0/max(1d0, abs(Hp_expected_ad%val), abs(s% Hp_face(k)))
            resid_ad = scale*(Hp_expected_ad - wrap_Hp_00(s, k))
            max_resid = max(max_resid, abs(resid_ad%val))
         end do
      end subroutine max_abs_rsp2_Hp_row_resid_for_star_LNA

      ! Small indexing, reporting, and public helper utilities.
      subroutine star_LNA_output_filename(s, suffix, filename)
         type(star_info), pointer :: s
         character(len=*), intent(in) :: suffix
         character(len=*), intent(out) :: filename
         character(len=512) :: output_directory
         character(len=512) :: file_prefix

         call ensure_star_LNA_output_directory(s, output_directory)
         call resolve_star_LNA_output_file_prefix(s, file_prefix)
         if (trim(output_directory) == '/') then
            filename = '/' // trim(file_prefix) // trim(suffix)
         else
            filename = trim(output_directory) // '/' // trim(file_prefix) // trim(suffix)
         end if
      end subroutine star_LNA_output_filename


      subroutine ensure_star_LNA_output_directory(s, output_directory)
         type(star_info), pointer :: s
         character(len=*), intent(out) :: output_directory

         call resolve_star_LNA_output_directory(s, output_directory)
         if (.not. folder_exists(trim(output_directory))) call mkdir(trim(output_directory))
      end subroutine ensure_star_LNA_output_directory


      subroutine resolve_star_LNA_output_directory(s, output_directory)
         type(star_info), pointer :: s
         character(len=*), intent(out) :: output_directory

         if (len_trim(s% star_LNA_output_directory) > 0) then
            output_directory = trim(s% star_LNA_output_directory)
         else if (len_trim(s% log_directory) > 0) then
            output_directory = trim(s% log_directory)
         else
            output_directory = 'LOGS'
         end if
         call normalize_star_LNA_output_directory(output_directory)
      end subroutine resolve_star_LNA_output_directory


      subroutine resolve_star_LNA_output_file_prefix(s, file_prefix)
         type(star_info), pointer :: s
         character(len=*), intent(out) :: file_prefix

         if (len_trim(s% star_LNA_output_file_prefix) > 0) then
            file_prefix = trim(s% star_LNA_output_file_prefix)
         else
            file_prefix = 'star_LNA'
         end if
         call normalize_star_LNA_output_file_prefix(file_prefix)
         if (len_trim(file_prefix) == 0) file_prefix = 'star_LNA'
      end subroutine resolve_star_LNA_output_file_prefix


      subroutine normalize_star_LNA_output_directory(output_directory)
         character(len=*), intent(inout) :: output_directory
         integer :: n

         output_directory = adjustl(output_directory)
         n = len_trim(output_directory)
         do while (n > 1)
            if (output_directory(n:n) /= '/' .and. &
                  output_directory(n:n) /= star_LNA_backslash) exit
            output_directory(n:n) = ' '
            n = len_trim(output_directory)
         end do
         if (len_trim(output_directory) == 0) output_directory = 'LOGS'
      end subroutine normalize_star_LNA_output_directory


      subroutine normalize_star_LNA_output_file_prefix(file_prefix)
         character(len=*), intent(inout) :: file_prefix
         integer :: i, n

         file_prefix = adjustl(file_prefix)
         n = len_trim(file_prefix)
         do while (n > 0)
            if (file_prefix(1:1) /= '/' .and. file_prefix(1:1) /= star_LNA_backslash) exit
            if (n == 1) then
               file_prefix = ''
            else
               file_prefix = file_prefix(2:n)
            end if
            n = len_trim(file_prefix)
         end do
         do while (n > 0)
            if (file_prefix(n:n) /= '/' .and. file_prefix(n:n) /= star_LNA_backslash) exit
            file_prefix(n:n) = ' '
            n = len_trim(file_prefix)
         end do
         do i = 1, n
            if (file_prefix(i:i) == '/' .or. file_prefix(i:i) == star_LNA_backslash) &
               file_prefix(i:i) = '_'
         end do
      end subroutine normalize_star_LNA_output_file_prefix


      integer function matrix_index(map, k, var_id) result(idx)
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: k, var_id
         integer :: slot

         idx = 0
         if (k < 1 .or. k > map%nz) return
         slot = var_slot(map, var_id)
         if (slot <= 0) return
         idx = (k - 1)*map%nvar_per_zone + slot
      end function matrix_index


      integer function var_slot(map, var_id) result(slot)
         type(star_LNA_var_map), intent(in) :: map
         integer, intent(in) :: var_id
         integer :: i

         slot = 0
         do i = 1, map%nvar_per_zone
            if (map%var_id(i) == var_id) then
               slot = i
               return
            end if
         end do
      end function var_slot


      integer function velocity_var_for_star_LNA(map) result(var_id)
         type(star_LNA_var_map), intent(in) :: map

         if (var_slot(map, lna_var_u) > 0) then
            var_id = lna_var_u
         else
            var_id = lna_var_v
         end if
      end function velocity_var_for_star_LNA


      logical function force_zero_velocity_for_star_LNA(s, k) result(force_zero_v)
         type(star_info), pointer :: s
         integer, intent(in) :: k

         force_zero_v = (s% q(k) > s% velocity_q_upper_bound) .or. &
            (s% tau(k) < s% velocity_tau_lower_bound) .or. &
            (s% lnT_start(k)/ln10 < s% velocity_logT_lower_bound .and. &
               s% dt < secyer*s% max_dt_yrs_for_velocity_logT_lower_bound)
      end function force_zero_velocity_for_star_LNA


      real(dp) function max_star_LNA_v_div_csound(s) result(max_v_div_csound)
         type(star_info), pointer :: s

         if (s% u_flag) then
            max_v_div_csound = maxval( &
               abs(s% u(1:s% nz))/max(1d-99, abs(s% csound(1:s% nz))))
         else
            max_v_div_csound = maxval( &
               abs(s% v(1:s% nz))/max(1d-99, abs(s% csound_face(1:s% nz))))
         end if
      end function max_star_LNA_v_div_csound


      subroutine report_star_LNA_setup(s, map)
         type(star_info), pointer :: s
         type(star_LNA_var_map), intent(in) :: map
         integer :: iv
         character(len=512) :: output_directory
         character(len=512) :: file_prefix

         write(*,'(a,i0)') 'star_LNA: nz = ', map%nz
         write(*,'(a,i0)') 'star_LNA: variables per zone = ', map%nvar_per_zone
         write(*,'(a,i0)') 'star_LNA: total perturbation variables = ', map%nvar_total
         write(*,'(a,i0)') 'star_LNA: requested modes = ', s% star_LNA_num_modes
         write(*,'(a,1pe12.4,a)') 'star_LNA: min selected frequency = ', &
            s% star_LNA_min_mode_frequency_uHz, ' microHz'
         write(*,'(a,1pe12.4)') 'star_LNA: first selected mode min logKE_per_cycle = ', &
            s% star_LNA_min_first_mode_eta
         write(*,'(a,1pe12.4)') 'star_LNA: later selected modes min logKE_per_cycle = ', &
            s% star_LNA_min_mode_eta
         if (s% star_LNA_max_abs_mode_eta > 0d0) then
            write(*,'(a,1pe12.4)') 'star_LNA: selected mode max abs logKE_per_cycle = ', &
               s% star_LNA_max_abs_mode_eta
         else
            write(*,'(a)') 'star_LNA: selected mode max abs logKE_per_cycle filter disabled'
         end if
         call resolve_star_LNA_output_directory(s, output_directory)
         call resolve_star_LNA_output_file_prefix(s, file_prefix)
         write(*,'(a,a)') 'star_LNA: output directory = ', trim(output_directory)
         write(*,'(a,a)') 'star_LNA: output file prefix = ', trim(file_prefix)
         write(*,'(a,a)') 'star_LNA: convection treatment = ', &
            trim(s% star_LNA_convection_treatment)
         write(*,'(a)', advance='no') 'star_LNA: variable set = '
         do iv = 1, map%nvar_per_zone
            if (iv > 1) write(*,'(a)', advance='no') ', '
            write(*,'(a)', advance='no') trim(var_name(map%var_id(iv)))
         end do
         write(*,*)

         if (s% RSP2_flag) then
            write(*,'(a)') 'star_LNA: RSP2 perturbations selected for Lr, Lc, Lt, w, and Hp.'
         else if (tdc_lna_active(s)) then
            write(*,'(a)') &
               'star_LNA: TDC perturbations selected for Lr, Lc, and internal w.'
         else if (frozen_flux_lna_active(s)) then
            write(*,'(a)') &
               'star_LNA: frozen_flux selected; perturbing Lrad with delta Lconv = 0.'
         else if (s% MLT_option == 'TDC' .and. s% star_LNA_include_tdc) then
            write(*,'(a)') 'star_LNA: TDC perturbations frozen by star_LNA_convection_treatment.'
         else if (s% MLT_option == 'TDC') then
            write(*,'(a)') 'star_LNA: TDC perturbations disabled by star_LNA_include_tdc.'
         else if (trim(s% star_LNA_convection_treatment) == 'mlt_static') then
            write(*,'(a)') &
               'star_LNA: static MLT temperature gradient rows; no convective velocity variable.'
         else
            write(*,'(a)') 'star_LNA: no dynamic convection variable selected.'
         end if
         if (s% u_flag) then
            write(*,'(a)') &
               'star_LNA: using cell velocity u and the Riemann face reconstruction.'
         else if (.not. s% v_flag) then
            write(*,'(a)') &
               'star_LNA: hydro v_flag is inactive; LNA uses its own face velocity perturbation.'
         end if
         if (s% use_TDC_enthalpy_flux_limiter .and. &
               (s% RSP2_flag .or. tdc_lna_active(s))) &
            write(*,'(a)') &
               'star_LNA: hydro enthalpy flux limiter is active; LNA ignores it.'
         if (s% using_velocity_time_centering) &
            write(*,'(a)') &
               'star_LNA: hydro velocity time centering is active; LNA uses continuous equations without time centering.'
         if (s% use_Pvsc_art_visc) &
            write(*,'(a)') &
               'star_LNA: hydro artificial viscosity pressure is active; LNA ignores it like RSP LINA.'
      end subroutine report_star_LNA_setup


      function var_name(var_id) result(name)
         integer, intent(in) :: var_id
         character(len=8) :: name

         select case (var_id)
         case (lna_var_lnd)
            name = 'lnd'
         case (lna_var_lnR)
            name = 'lnR'
         case (lna_var_v)
            name = 'v'
         case (lna_var_u)
            name = 'u'
         case (lna_var_lnT)
            name = 'lnT'
         case (lna_var_L)
            name = 'L'
         case (lna_var_w)
            name = 'w'
         case (lna_var_Hp)
            name = 'Hp'
         case default
            name = 'unknown'
         end select
      end function var_name


      function equation_name_for_star_LNA(eq_id) result(name)
         integer, intent(in) :: eq_id
         character(len=32) :: name

         select case (eq_id)
         case (lna_eq_density)
            name = 'density'
         case (lna_eq_radius)
            name = 'radius'
         case (lna_eq_zero_velocity)
            name = 'zero_velocity'
         case (lna_eq_momentum)
            name = 'momentum'
         case (lna_eq_surface_momentum)
            name = 'surface_momentum'
         case (lna_eq_surface_pressure)
            name = 'surface_pressure'
         case (lna_eq_surface_fixed_velocity)
            name = 'surface_fixed_velocity'
         case (lna_eq_energy)
            name = 'energy'
         case (lna_eq_luminosity_rsp_surface)
            name = 'luminosity_rsp_surface'
         case (lna_eq_luminosity_rsp2)
            name = 'luminosity_rsp2'
         case (lna_eq_luminosity_tdc)
            name = 'luminosity_tdc'
         case (lna_eq_luminosity_frozen_flux)
            name = 'luminosity_frozen_flux'
         case (lna_eq_surface_temperature)
            name = 'surface_temperature'
         case (lna_eq_temperature_gradient)
            name = 'temperature_gradient'
         case (lna_eq_rsp2_turbulent_energy)
            name = 'rsp2_turbulent_energy'
         case (lna_eq_rsp2_zero_w)
            name = 'rsp2_zero_w'
         case (lna_eq_rsp2_Hp)
            name = 'rsp2_Hp'
         case (lna_eq_tdc_velocity)
            name = 'tdc_velocity'
         case (lna_eq_tdc_zero_w)
            name = 'tdc_zero_w'
         case (lna_eq_mlt_static_temperature_gradient)
            name = 'mlt_static_temperature_gradient'
         case default
            name = 'unknown'
         end select
      end function equation_name_for_star_LNA


      subroutine free_star_LNA_problem(problem)
         type(star_LNA_problem), intent(inout) :: problem

         if (allocated(problem% eq_id)) deallocate(problem% eq_id)
         call free_star_LNA_matrix(problem% mtx)
         call free_star_LNA_var_map(problem% map)
      end subroutine free_star_LNA_problem


      subroutine free_star_LNA_var_map(map)
         type(star_LNA_var_map), intent(inout) :: map
         if (allocated(map%var_id)) deallocate(map%var_id)
         map%nz = 0
         map%nvar_per_zone = 0
         map%nvar_total = 0
      end subroutine free_star_LNA_var_map


      subroutine free_star_LNA_matrix(mtx)
         type(star_LNA_matrix), intent(inout) :: mtx
         if (allocated(mtx%A)) deallocate(mtx%A)
         if (allocated(mtx%B)) deallocate(mtx%B)
         mtx%n = 0
      end subroutine free_star_LNA_matrix


      function star_LNA_L_conv_ad(s, k) result(L_conv_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: L_conv_ad

         L_conv_ad = star_LNA_L_conv_closure_ad(s, k)
      end function star_LNA_L_conv_ad


      end module star_lna_support
