! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module init

      use star_private_def
      use const_def

      implicit none

      private
      public :: alloc_star_data, set_starting_star_data, do_star_init, &
         do_starlib_shutdown, set_kap_and_eos_handles, load_zams_model, &
         create_pre_ms_model, create_initial_model, create_RSP_model, &
         doing_restart, load_restart_photo, load_saved_model, &
         load_saved_RSP_model, do_garbage_collection

      integer, parameter :: do_create_pre_ms_model = 0
      integer, parameter :: do_load_zams_model = 1
      integer, parameter :: do_load_saved_model = 2
      integer, parameter :: do_create_initial_model = 3
      integer, parameter :: do_create_RSP_model = 4
      integer, parameter :: do_load_saved_RSP_model = 5
      

      logical :: have_done_starlib_init = .false.

      contains


      subroutine set_kap_and_eos_handles(id, ierr)
         use kap_lib, only: alloc_kap_handle_using_inlist, kap_ptr
         use eos_lib, only: alloc_eos_handle_using_inlist, eos_ptr
         integer, intent(in) :: id
         integer, intent(out) :: ierr ! 0 means AOK.
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_star_ptr failed in alloc_eos_handle'
            return
         end if
         if (s% eos_handle == 0) then
            s% eos_handle = alloc_eos_handle_using_inlist(s% inlist_fname, ierr)
            if (ierr /= 0) then
               write(*,*) 'set_kap_and_eos_handles failed in alloc_eos_handle_using_inlist'
               return
            end if
            call eos_ptr(s% eos_handle, s% eos_rq, ierr)
            if (ierr /= 0) then
               write(*,*) 'eos_ptr failed in alloc_eos_handle'
               return
            end if            
         end if
         if (s% kap_handle == 0) then
            s% kap_handle = alloc_kap_handle_using_inlist(s% inlist_fname, ierr)
            if (ierr /= 0) then
               write(*,*) 'set_kap_and_eos_handles failed in alloc_kap_handle_using_inlist'
               return
            end if
            call kap_ptr(s% kap_handle,s% kap_rq,ierr)
            if (ierr /= 0) then
               write(*,*) 'kap_ptr failed in alloc_kap_handle'
               return
            end if
         end if
      end subroutine set_kap_and_eos_handles


      subroutine do_star_init( &
            my_mesa_dir, chem_isotopes_filename, &
            net_reaction_filename, jina_reaclib_filename, &
            use_suzuki_weak_rates, &
            use_special_weak_rates, special_weak_states_file, special_weak_transitions_file, &
            reaclib_min_T9, &
            rate_tables_dir, rates_cache_suffix, &
            eosDT_cache_dir, &
            kap_cache_dir, rates_cache_dir, &
            color_num_files,color_file_names,color_num_colors,&
            ierr)
         use paquette_coeffs, only: initialise_collision_integrals
         use hydro_rotation, only: init_rotation
         use alloc, only: init_alloc
         character (len=*), intent(in) :: &
            my_mesa_dir, chem_isotopes_filename, net_reaction_filename, &
            jina_reaclib_filename, rate_tables_dir, &
            special_weak_states_file, special_weak_transitions_file, &
            rates_cache_suffix, &
            eosDT_cache_dir, &
            kap_cache_dir, rates_cache_dir
         logical, intent(in) :: use_suzuki_weak_rates, use_special_weak_rates
         real(dp), intent(in) :: reaclib_min_T9
         integer, intent(in) :: color_num_files
         character (len=*), intent(in) :: color_file_names(:)
         integer , intent(in):: color_num_colors(:)
         integer, intent(out) :: ierr
         integer :: iam, nprocs, nprow, npcol, i, n
         integer, dimension(:), allocatable :: seed
         include 'formats'
         ierr = 0
         if (have_done_starlib_init) return
         call initialise_collision_integrals
         call init_alloc
         call stardata_init( &
            my_mesa_dir, chem_isotopes_filename, &
            net_reaction_filename, jina_reaclib_filename, &
            use_suzuki_weak_rates, &
            use_special_weak_rates, special_weak_states_file, special_weak_transitions_file, &
            reaclib_min_T9, &
            rate_tables_dir, rates_cache_suffix, &
            eosDT_cache_dir, &
            kap_cache_dir, rates_cache_dir, &
            color_num_files,color_file_names,color_num_colors,&
            ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in stardata_init'
            return
         end if

         call init_rotation(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in init_rotation'
            return
         end if

         have_done_starlib_init = .true.

      end subroutine do_star_init


      subroutine do_starlib_shutdown

        use alloc, only: shutdown_alloc
        use micro, only: shutdown_microphys
        use atm_lib, only: atm_shutdown
        use colors_lib, only: colors_shutdown
        use chem_lib, only: chem_shutdown
        use rates_lib, only: rates_shutdown
        use star_profile_def, only: profile_column_names_shutdown
        use star_history_def, only: history_column_names_shutdown
        use paquette_coeffs, only: free_collision_integrals

        call rates_shutdown()
        call atm_shutdown()
        call shutdown_microphys()
        call colors_shutdown()
        call chem_shutdown()
        call profile_column_names_shutdown()
        call history_column_names_shutdown()

        call free_collision_integrals()

        call shutdown_alloc()

        have_done_starlib_init = .false.

      end subroutine do_starlib_shutdown


      subroutine alloc_star_data(id, ierr)
         use rates_def, only: rates_reaction_id_max, rates_NACRE_if_available
         use chem_def, only: num_categories
         use net, only: default_set_which_rates, default_set_rate_factors, &
            default_set_op_mono_factors


         integer, intent(out) :: id, ierr

         type (star_info), pointer :: s
         integer, parameter :: init_alloc_nvar = 20
         character (len=32) :: extra_name
         integer :: i

         ierr = 0

         call alloc_star(id, ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_star_data failed in alloc_star'
            return
         end if

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_star_data failed in get_star_ptr'
            return
         end if

         nullify(s% dq)
         nullify(s% xa)
         nullify(s% xh)
         
         nullify( &
            s% op_mono_umesh1, s% op_mono_semesh1, s% op_mono_ff1, &
            s% op_mono_rs1)

         nullify(s% atm_structure)
         s% atm_structure_num_pts = 0

         nullify(s% chem_id)
         nullify(s% xa_removed)

         nullify(s% which_rates)
         s% set_which_rates => default_set_which_rates

         nullify(s% rate_factors)
         s% set_rate_factors => default_set_rate_factors

         nullify(s% op_mono_factors)
         s% set_op_mono_factors => default_set_op_mono_factors

         allocate(s% nameofvar(init_alloc_nvar),stat=ierr)
         if (ierr /= 0) return

         allocate(s% nameofequ(init_alloc_nvar),stat=ierr)
         if (ierr /= 0) return

         nullify(s% history_values)
         nullify(s% history_value_is_integer)
         nullify(s% history_names)
         nullify(s% history_names_dict)

         nullify(s% prev_mesh_xh)
         nullify(s% prev_mesh_xa)
         nullify(s% prev_mesh_j_rot)
         nullify(s% prev_mesh_omega)
         nullify(s% prev_mesh_dq)

         nullify(s% other_star_info)

         nullify(s% bcyclic_odd_storage)
         nullify(s% bcyclic_odd_storage_qp)

         nullify(s% solver_iwork)
         nullify(s% solver_work)
         nullify(s% AF1)

         s% net_name = ''
         s% species = 0
         s% num_reactions = 0

         s% M_center = 0
         s% R_center = 0
         s% v_center = 0
         s% L_center = 0

         nullify(s% profile_column_spec)
         nullify(s% history_column_spec)

         s% num_conv_boundaries = 0
         nullify(s% conv_bdy_q)
         nullify(s% conv_bdy_loc)
         nullify(s% top_conv_bdy)

         s% num_mix_boundaries = 0
         nullify(s% mix_bdy_q)
         nullify(s% mix_bdy_loc)
         nullify(s% top_mix_bdy)

         nullify(s% burn_h_conv_region)
         nullify(s% burn_he_conv_region)
         nullify(s% burn_z_conv_region)

         s% have_burner_storage = .false.
         s% burner_storage_sz_per_thread = 0
         s% burner_num_threads = 0
         nullify(s% burner_storage)

         s% doing_timing = .false.
         s% time_evolve_step = 0
         s% time_remesh = 0
         s% time_adjust_mass = 0
         s% time_conv_premix = 0
         s% time_element_diffusion = 0
         s% time_struct_burn_mix = 0
         s% time_solver_matrix = 0
         s% time_solve_mix = 0
         s% time_solve_burn = 0
         s% time_solve_omega_mix = 0
         s% time_eos = 0
         s% time_neu_kap = 0
         s% time_nonburn_net = 0
         s% time_mlt = 0
         s% time_set_hydro_vars = 0
         s% time_set_mixing_info = 0
         s% time_total = 0

         s% timing_num_get_eos_calls = 0
         s% timing_num_solve_eos_calls = 0
         s% timing_num_get_kap_calls = 0

         s% model_profile_filename = ''
         s% most_recent_profile_filename = ''

         s% model_controls_filename = ''
         s% most_recent_controls_filename = ''

         s% most_recent_photo_name = ''

         s% doing_flash_wind = .false.
         s% doing_rlo_wind = .false.

         s% phase_of_evolution = phase_starting
         s% recent_log_header = -1000

         s% tau_base = 2d0/3d0
         s% tau_factor = 1

         s% solver_iter = 0

         s% using_gold_tolerances = .false.
         s% using_velocity_time_centering = .false.

         s% using_revised_net_name = .false.
         s% revised_net_name = ''

         s% using_revised_max_yr_dt = .false.
         s% revised_max_yr_dt = 0
 
         s% astero_using_revised_max_yr_dt = .false.
         s% astero_revised_max_yr_dt = 0

         s% cumulative_energy_error = 0
         s% cumulative_extra_heating = 0

         s% have_initial_energy_integrals = .false.

         s% num_solver_iterations = 0         
         s% bad_max_corr_cnt = 0

         s% mesh_call_number = 0
         s% solver_call_number = 0
         s% diffusion_call_number = 0
         s% model_number = 0
         s% RSP_have_set_velocities = .false.
         s% RSP_just_set_velocities = .false.
         
         s% dt = 0d0
         s% mstar_dot = 0d0
         s% boost_mlt_alfa = 0

         s% power_nuc_burn = -1
         s% power_h_burn = -1
         s% power_he_burn = -1
         s% power_z_burn = -1
         s% power_photo = -1

         s% k_const_mass = 1
         s% k_below_just_added = 1
         s% k_below_const_q = 1

         s% why_Tlim = Tlim_struc
         s% dt_why_count(:) = 0
         s% dt_why_retry_count(:) = 0

         s% len_extra_iwork = 0
         s% len_extra_work = 0

         s% eos_handle = 0
         s% kap_handle = 0

         do i = 1, max_num_profile_extras
            if (i < 10) then
               write(extra_name,'(a,i1)') 'profile_extra_', i
            else if (i < 100) then
               write(extra_name,'(a,i2)') 'profile_extra_', i
            else
               write(extra_name,'(a,i3)') 'profile_extra_', i
            end if
            s% profile_extra_name(i) = trim(extra_name)
         end do

         call set_starting_star_data(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_star_data failed in set_starting_star_data'
            return
         end if

      end subroutine alloc_star_data


      subroutine null_other_new_generation(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_other_new_generation


      subroutine null_other_set_current_to_old(id)
         integer, intent(in) :: id
      end subroutine null_other_set_current_to_old


      subroutine set_starting_star_data(s, ierr)
         use other_wind, only: null_other_wind
         use other_accreting_state, only: null_other_accreting_state
         use other_adjust_mdot, only: null_other_adjust_mdot
         use other_j_for_adjust_J_lost, only: null_other_j_for_adjust_J_lost
         use other_eval_fp_ft, only: null_other_eval_fp_ft
         use other_eval_i_rot, only: null_other_eval_i_rot
         use other_torque, only: default_other_torque
         use other_torque_implicit, only: default_other_torque_implicit
         use other_remove_surface, only: default_other_remove_surface
         use other_momentum, only: default_other_momentum
         use other_momentum_implicit, only: default_other_momentum_implicit
         use other_pressure, only: default_other_pressure
         use other_energy, only: default_other_energy
         use other_energy_implicit, only: default_other_energy_implicit
         use other_D_mix, only: null_other_D_mix
         use other_am_mixing, only: null_other_am_mixing
         use other_brunt, only: default_other_brunt
         use other_brunt_smoothing, only: null_other_brunt_smoothing
         use other_adjust_mlt_gradT_fraction, only: &
            default_other_adjust_mlt_gradT_fraction
         use other_after_set_mixing_info, only: &
            default_other_after_set_mixing_info
         use other_after_solver_setmatrix, only: &
            default_other_after_solver_setmatrix
         use other_diffusion, only: null_other_diffusion
         use other_diffusion_factor, only: default_other_diffusion_factor
         use other_neu, only: null_other_neu
         use other_net_get, only: null_other_net_get
         use other_cgrav, only: default_other_cgrav
         use other_mesh_delta_coeff_factor, only: default_other_mesh_delta_coeff_factor
         use other_alpha_mlt, only: default_other_alpha_mlt
         use other_mlt_results, only: null_other_mlt_results
         use other_opacity_factor, only: default_other_opacity_factor
         use other_pgstar_plots, only: null_other_pgstar_plots_info
         use other_mesh_functions
         use other_eps_grav, only: null_other_eps_grav
         use other_rsp_build_model, only: null_other_rsp_build_model
         use other_rsp_linear_analysis, only: null_other_rsp_linear_analysis
         use other_gradr_factor, only: null_other_gradr_factor
         use other_surface_PT, only: null_other_surface_PT
         use other_timestep_limit, only: null_other_timestep_limit
         use other_overshooting_scheme, only: null_other_overshooting_scheme
         use other_photo_write, only: default_other_photo_write
         use other_photo_read, only: default_other_photo_read
         use other_set_pgstar_controls, only: default_other_set_pgstar_controls
         use other_eos
         use other_kap
         use pgstar_decorator
         use star_utils, only: init_random

         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         ! note: keep the handles for eos, kap, and net

         ierr = 0

         s% model_number = 0
         s% time = 0
         s% dt = 0
         
         s% total_num_solver_iterations = 0
         s% total_num_solver_relax_iterations = 0
         s% total_num_solver_calls_made = 0
         s% total_num_solver_relax_calls_made = 0
         s% total_num_solver_calls_converged = 0
         s% total_num_solver_relax_calls_converged = 0
         
         s% num_solver_iterations = 0
         s% num_skipped_setvars = 0
         s% num_setvars = 0
         s% num_solver_setvars = 0
         s% num_retries = 0
         s% num_hydro_merges = 0
         s% num_hydro_splits = 0
         s% timestep_hold = -1
         
         s% mesh_call_number = 0
         s% solver_call_number = 0
         s% diffusion_call_number = 0
         s% model_number_for_last_retry = 0
         s% dt_limit_ratio = 0
         s% force_tau_factor = 0
         s% force_Tsurf_factor = 0
         s% force_opacity_factor = 0

         s% generations = 0

         s% nvar_hydro = 0
         s% nvar_chem = 0
         s% nvar_total = 0

         s% nz = 0
         s% prev_mesh_nz = 0
         s% net_name = ''
         s% species = 0
         s% num_reactions = 0

         s% v_flag = .false.
         s% u_flag = .false.
         s% rotation_flag = .false.
         s% RTI_flag = .false.
         s% conv_vel_flag = .false.
         s% w_div_wc_flag = .false.
         s% D_omega_flag = .false.
         s% am_nu_rot_flag = .false.
         s% RSP_flag = .false.
         s% RSP2_flag = .false.
         s% using_RSP2 = .false.
         s% using_TDC = .false.
         
         s% have_mixing_info = .false.
         s% doing_solver_iterations = .false.
         s% need_to_setvars = .true.
         s% okay_to_set_mixing_info = .true.
         s% okay_to_set_mlt_vc = .false. ! not until have set mlt_cv_old
         s% need_to_reset_w = .false.

         s% just_wrote_terminal_header = .false.
         s% doing_relax = .false.
         s% mstar_dot = 0

         s% surf_lnT = 0
         s% surf_lnd = 0
         s% surf_lnR = 0
         s% surf_v = 0
         s% surf_lnS = 0

         s% termination_code = -1

         s% prev_create_atm_R0_div_R = -1

         s% screening_mode_value = -1

         s% dt = -1
         s% dt_next = -1
         s% dt_next_unclipped = -1

         s% i_lnd = 0
         s% i_lnT = 0
         s% i_lnR = 0
         s% i_lum = 0
         s% i_v = 0
         s% i_u = 0
         s% i_alpha_RTI = 0

         s% i_dv_dt = 0
         s% i_du_dt = 0
         s% i_equL = 0
         s% i_dlnd_dt = 0
         s% i_dlnE_dt = 0
         s% i_dlnR_dt = 0
         s% i_dalpha_RTI_dt = 0

         s% op_mono_nptot = 0
         s% op_mono_ipe = 0
         s% op_mono_nrad = 0
         s% op_mono_n = 0

         s% number_of_history_columns = -1
         s% model_number_of_history_values = -1
         s% need_to_set_history_names_etc = .true.
         s% doing_finish_load_model = .false.

         nullify(s% finish_relax_step)
         nullify(s% finished_relax)

         s% how_many_extra_profile_header_items => &
            null_how_many_extra_header_items
         s% data_for_extra_profile_header_items => &
            null_data_for_extra_header_items

         s% how_many_extra_history_header_items => &
            null_how_many_extra_header_items
         s% data_for_extra_history_header_items => &
            null_data_for_extra_header_items

         s% other_wind => null_other_wind
         s% other_accreting_state => null_other_accreting_state
         s% other_adjust_mdot => null_other_adjust_mdot
         s% other_j_for_adjust_J_lost => null_other_j_for_adjust_J_lost
         s% other_D_mix => null_other_D_mix
         s% other_am_mixing => null_other_am_mixing
         s% other_eval_fp_ft => null_other_eval_fp_ft
         s% other_eval_i_rot => null_other_eval_i_rot
         s% other_torque => default_other_torque
         s% other_torque_implicit => default_other_torque_implicit
         s% other_remove_surface => default_other_remove_surface
         s% other_momentum => default_other_momentum
         s% other_momentum_implicit => default_other_momentum_implicit
         s% other_pressure => default_other_pressure
         s% other_energy => default_other_energy
         s% other_energy_implicit => default_other_energy_implicit
         s% other_cgrav => default_other_cgrav
         s% other_mesh_delta_coeff_factor => default_other_mesh_delta_coeff_factor
         s% other_alpha_mlt => default_other_alpha_mlt
         s% other_opacity_factor => default_other_opacity_factor
         s% other_brunt => default_other_brunt
         s% other_brunt_smoothing => null_other_brunt_smoothing
         s% other_adjust_mlt_gradT_fraction => default_other_adjust_mlt_gradT_fraction
         s% other_after_set_mixing_info => default_other_after_set_mixing_info
         s% other_after_solver_setmatrix => default_other_after_solver_setmatrix
         s% other_diffusion => null_other_diffusion
         s% other_diffusion_factor => default_other_diffusion_factor
         s% other_mlt_results => null_other_mlt_results
         s% other_neu => null_other_neu
         s% other_net_get => null_other_net_get
         s% other_eps_grav => null_other_eps_grav
         s% other_rsp_build_model => null_other_rsp_build_model
         s% other_rsp_linear_analysis => null_other_rsp_linear_analysis
         s% other_gradr_factor => null_other_gradr_factor

         s% other_eosDT_get => null_other_eosDT_get
         s% other_eosDT_get_T => null_other_eosDT_get_T
         s% other_eosDT_get_Rho => null_other_eosDT_get_Rho

         s% other_eosDE_get => null_other_eosDE_get

         s% other_kap_get => null_other_kap_get
         s% other_kap_get_op_mono => null_other_kap_get_op_mono

         s% other_surface_PT => null_other_surface_PT

         s% other_timestep_limit => null_other_timestep_limit
         s% other_overshooting_scheme => null_other_overshooting_scheme

         s% other_pgstar_plots_info => null_other_pgstar_plots_info
         s% how_many_other_mesh_fcns => null_how_many_other_mesh_fcns
         s% other_mesh_fcn_data => null_other_mesh_fcn_data

         s% other_photo_write => default_other_photo_write
         s% other_photo_read => default_other_photo_read

         s% other_set_pgstar_controls => default_other_set_pgstar_controls

         s% other_new_generation => null_other_new_generation
         s% other_set_current_to_old => null_other_set_current_to_old

         nullify(s% extra_profile_col_names)
         nullify(s% extra_profile_col_vals)
         s% num_extra_profile_cols = 0

         s% binary_id = 0
         s% include_binary_history_in_log_file = .false.
         s% how_many_binary_history_columns => null_how_many_binary_history_columns
         s% data_for_binary_history_columns => null_data_for_binary_history_columns
         s% how_many_extra_binary_history_columns => null_how_many_extra_binary_history_columns
         s% data_for_extra_binary_history_columns => null_data_for_extra_binary_history_columns
         
         s% Abundance_pgstar_decorator => null_pgstar_decorator

         s% generations = 0

         s% nz = 0

         s% nvar_hydro = 0
         s% nvar_chem = 0
         s% nvar_total = 0

         s% species = 0
         s% num_reactions = 0

         s% model_number = 0
         s% mstar = 0
         s% xmstar = 0
         s% M_center = 0
         s% v_center = 0
         s% R_center = 0
         s% L_center = 0
         s% time = 0
         s% total_angular_momentum = 0
         s% prev_create_atm_R0_div_R = 0

         s% dt = 0

         s% have_previous_conv_vel = .false.
         s% have_mlt_vc = .false.

         s% net_name = ''

         s% mstar_dot = 0
         s% v_surf = 0
         s% L_nuc_burn_total = 0
         s% L_by_category = 0
         s% dt_limit_ratio = 0
         s% L_phot = 0
         s% T_surf = 0
         s% P_surf = 0
         
         s% gradT_excess_alpha = 0

         s% h1_czb_mass = 0

         s% he_core_mass = 0
         s% co_core_mass = 0
         s% Teff = -1 ! need to calculate it
         s% center_eps_nuc = 0
         s% Lrad_div_Ledd_avg_surf = 0
         s% w_div_w_crit_avg_surf = 0         
         s% total_internal_energy = 0d0
         s% total_gravitational_energy = 0d0
         s% total_radial_kinetic_energy = 0d0
         s% total_turbulent_energy = 0d0
         s% total_rotational_kinetic_energy = 0d0
         s% total_energy = 0d0

         s% n_conv_regions = 0
         s% cz_bot_mass(:) = 0
         s% cz_top_mass(:) = 0
         s% dt_next = 0

         s% model_profile_filename = ''
         s% model_controls_filename = ''
         s% model_data_filename = ''

         s% most_recent_profile_filename = ''
         s% most_recent_controls_filename = ''

         s% most_recent_model_data_filename = ''

         s% recent_log_header = -1000
         s% phase_of_evolution = 0

         s% num_solver_iterations = 0
         s% num_skipped_setvars = 0
         s% num_setvars = 0
         s% num_retries = 0

         s% num_hydro_merges = 0
         s% num_hydro_splits = 0

         s% timestep_hold = 0
         s% model_number_for_last_retry = 0

         s% mesh_call_number = 0
         s% solver_call_number = 0
         s% diffusion_call_number = 0

         s% num_rotation_solver_steps = 0
         s% num_diffusion_solver_steps = 0
         s% initial_timestep = 0
         s% why_Tlim = 0
         s% result_reason = 0

         s% have_new_generation = .false.
         s% have_new_cz_bdy_info = .false.
         s% need_to_update_history_now = .false.
         s% need_to_save_profiles_now = .false.
         s% save_profiles_model_priority = 0

         s% doing_flash_wind = .false.
         s% doing_rlo_wind = .false.
         s% most_recent_photo_name = ''

         s% len_extra_iwork = 0
         s% len_extra_work = 0
         
         call init_random(s)

      end subroutine set_starting_star_data


      subroutine create_pre_ms_model(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         character (len=0) :: model_dir
         call model_builder( &
            id, model_dir, do_create_pre_ms_model, &
            .false., 'restart_photo', ierr)
      end subroutine create_pre_ms_model


      subroutine create_initial_model(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         character (len=0) :: model_dir
         call model_builder( &
            id, model_dir, do_create_initial_model, &
            .false., 'restart_photo', ierr)
      end subroutine create_initial_model


      subroutine create_RSP_model(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         character (len=0) :: model_dir
         call model_builder( &
            id, model_dir, do_create_RSP_model, &
            .false., 'restart_photo', ierr)
      end subroutine create_RSP_model


      subroutine load_zams_model(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         call model_builder( &
            id, '', do_load_zams_model, &
            .false., 'restart_photo', ierr)
      end subroutine load_zams_model


      subroutine load_saved_RSP_model(id, model_fname, ierr)
         integer, intent(in) :: id
         character (len=*), intent(in) :: model_fname
         integer, intent(out) :: ierr
         integer :: l
         l = len_trim(model_fname)
         call model_builder( &
            id, model_fname, do_load_saved_RSP_model, &
            .false., 'restart_photo', ierr)
      end subroutine load_saved_RSP_model


      subroutine load_saved_model(id, model_fname, ierr)
         integer, intent(in) :: id
         character (len=*), intent(in) :: model_fname
         integer, intent(out) :: ierr
         integer :: l
         l = len_trim(model_fname)
         call model_builder( &
            id, model_fname, do_load_saved_model, &
            .false., 'restart_photo', ierr)
      end subroutine load_saved_model


      subroutine load_restart_photo(id, restart_filename, ierr)
         integer, intent(in) :: id
         character (len=*), intent(in) :: restart_filename
         integer, intent(out) :: ierr
         call model_builder( &
            id, '', do_load_zams_model, .true., restart_filename, ierr)
      end subroutine load_restart_photo


      ! for both zams and pre-main-sequence
      subroutine model_builder( &
            id, model_info, do_which, restart, restart_filename, ierr)
         use net, only: set_net, do_micro_change_net
         use alloc, only: set_var_info
         use photo_in, only: read_star_photo
         use init_model, only: get_zams_model
         use star_utils, only: set_phase_of_evolution, yrs_for_init_timestep, &
            eval_deltaM_total_energy_integrals
         use adjust_xyz, only: set_z, set_y
         use pre_ms_model, only: build_pre_ms_model
         use create_initial_model, only: build_initial_model
         use rsp, only: build_rsp_model, rsp_total_energy_integrals
         use rsp_def, only: init_def
         use read_model, only: do_read_saved_model, &
            finish_load_model
         use relax, only: do_relax_to_limit, do_relax_mass, &
            do_relax_mass_scale, do_relax_num_steps, do_relax_to_radiative_core
         use auto_diff_support
         integer, intent(in) :: id, do_which
         character (len=*), intent(in) :: model_info, restart_filename
         logical, intent(in) :: restart
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         real(dp) :: initial_mass, initial_z, dlgm_per_step
         logical :: want_rsp_model, is_rsp_model
         real(dp), parameter :: lg_max_abs_mdot = -1000 ! use default
         real(dp), parameter :: change_mass_years_for_dt = 1
         real(dp), parameter :: min_mass_for_create_pre_ms = 0.03d0
         logical :: restore_at_end
         real(dp) :: xm, total_radiation, warning_limit_for_max_residual
         integer :: k, num_trace_history_values
         real(dp) :: save_Pextra_factor
         character (len=256):: save_atm_option, &
            save_atm_T_tau_relation, save_atm_T_tau_opacity
         real(dp), allocatable, dimension(:) :: total_energy_profile

         include 'formats'

         ierr = 0

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         initial_mass = s% initial_mass
         initial_z = s% initial_z
         s% dt = 0
         s% termination_code = -1
         want_rsp_model = .false.
         is_rsp_model = .false.

         if (restart) then
            s% doing_first_model_of_run = .false.
            s% doing_first_model_after_restart = .true.
            call read_star_photo(s, restart_filename, ierr)
            if (ierr /= 0) return
            call check_initials
            call set_net(s, s% net_name, ierr)
            if (ierr /= 0) return
            if (s% rotation_flag) s% have_j_rot = .true.
            call init_def(s) ! RSP
            call finish_load_model(s, restart, want_rsp_model, is_rsp_model, ierr)
            if (s% max_years_for_timestep > 0) &
               s% dt_next = min(s% dt_next, secyer*s% max_years_for_timestep)
            return
         end if

         num_trace_history_values = s% num_trace_history_values
         s% num_trace_history_values = 0
         
         warning_limit_for_max_residual = s% warning_limit_for_max_residual
         s% warning_limit_for_max_residual = 1d0
         
         s% doing_first_model_of_run = .true.
         s% doing_first_model_after_restart = .false.
         
         if (do_which == do_load_saved_model .or. do_which == do_load_saved_RSP_model) then
            s% dt_next = -1
            want_rsp_model = (do_which == do_load_saved_RSP_model)
            call do_read_saved_model(s, model_info, want_rsp_model, is_rsp_model, ierr)
            if (ierr /= 0) then
               write(*,*) 'load failed in do_read_saved_model'
               return
            end if
            call check_initials
            if (s% dt_next < 0) s% dt_next = yrs_for_init_timestep(s)*secyer
         else
            s% net_name = s% default_net_name
            s% species = 0
            s% v_flag = .false.
            s% u_flag = .false.
            s% rotation_flag = .false.
            s% D_omega_flag = .false.
            s% am_nu_rot_flag = .false.
            s% star_mass = s% initial_mass
            s% mstar = s% initial_mass*Msun
            s% M_center = s% mstar - s% xmstar
            call set_var_info(s, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in set_var_info'
               return
            end if
            select case (do_which)
               case (do_create_initial_model)
                  if (s% initial_model_change_net) then
                     call do_micro_change_net(s, s% initial_model_new_net_name, ierr)
                  else
                     call set_net(s, s% net_name, ierr)
                  end if
                  if (ierr /= 0) then
                     write(*,*) 'failed in set_net'
                     return
                  end if
                  call build_initial_model(s, ierr)
                  if (ierr /= 0) then
                     write(*,*) 'failed in build_initial_model'
                     return
                  end if
                  s% generations = 1
                  s% dt_next = 1d-5*secyer
               case (do_create_pre_ms_model)
                  if (s% initial_mass < min_mass_for_create_pre_ms) then
                     write(*,*)
                     write(*,*)
                     write(*,*)
                     write(*,'(a,1x,f5.2)') 'sorry: cannot create pre-ms smaller than', &
                        min_mass_for_create_pre_ms
                     write(*,'(a)') &
                        'please create pre-ms and then relax to lower mass as a separate operation'
                     write(*,*)
                     write(*,'(a)') 'here is an example:'
                     write(*,'(a)') 'in your inlist &controls section, set initial_mass = 0.03'
                     write(*,'(a)') 'in the &star_job section, add something like these lines'
                     write(*,'(a)') '  relax_mass_scale = .true.'
                     write(*,'(a)') '  dlgm_per_step = 1d-3 ! log10(delta M/Msun/step)'
                     write(*,'(a)') '  new_mass = 2.863362d-3 ! 3 Mjupiter in Msun units'
                     write(*,'(a)') '  change_mass_years_for_dt = 1'
                     write(*,*)
                     write(*,*)
                     write(*,*)
                     ierr = -1
                     return
                  end if
                  if (s% pre_ms_change_net) then
                     call do_micro_change_net(s, s% pre_ms_new_net_name, ierr)
                  else
                     call set_net(s, s% net_name, ierr)
                  end if
                  if (ierr /= 0) then
                     write(*,*) 'failed in set_net'
                     return
                  end if
                  write(*,2) 'species  mass', s% species, s% mstar/Msun

                  call build_pre_ms_model(id, s, s% nvar_hydro, s% species, ierr)
                  if (ierr /= 0) then
                     write(*,*) 'failed in build_pre_ms_model'
                     return
                  end if
                  s% generations = 1
                  s% dt_next = 1d-5*secyer
               case (do_load_zams_model)
                  s% generations = 1
                  call get_zams_model(s, s% zams_filename, ierr)
                  if (ierr /= 0) then
                     write(*,*) 'failed in get_zams_model'
                     return
                  end if
                  if (s% dt_next <= 0d0) then
                     s% dt_next = yrs_for_init_timestep(s)*secyer
                  end if
               case (do_create_RSP_model)
                  call build_rsp_model(s, ierr) ! like build_pre_ms_model
                  if (ierr /= 0) then
                     write(*,*) 'failed in build_rsp_model'
                     return
                  end if
                  s% generations = 1
                  s% dt_next = 1d-2*secyer
               case default
                  write(*,*) 'bad value for do_which in build_model'
                  ierr = -1
                  return
            end select
         end if
         
         do k=1,s% nz
            s% extra_heat(k) = 0
         end do

         call finish_load_model(s, restart, want_rsp_model, is_rsp_model, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in finish_load_model'
            return
         end if
         
         if (s% max_years_for_timestep > 0) &
            s% dt_next = min(s% dt_next, secyer*s% max_years_for_timestep)
         call set_phase_of_evolution(s)

         if (s% RSP_flag) then
            call rsp_total_energy_integrals(s, &
               s% total_internal_energy, &
               s% total_gravitational_energy, &
               s% total_radial_kinetic_energy, &
               s% total_rotational_kinetic_energy, &
               s% total_turbulent_energy, &
               s% total_energy, total_radiation)
         else
            allocate(total_energy_profile(1:s% nz))
            call eval_deltaM_total_energy_integrals( &
               s, 1, s% nz, s% mstar, .false., &
               total_energy_profile, &
               s% total_internal_energy, &
               s% total_gravitational_energy, &
               s% total_radial_kinetic_energy, &
               s% total_rotational_kinetic_energy, &
               s% total_turbulent_energy, &
               s% total_energy)
         end if

         if (do_which == do_create_pre_ms_model) then
            call setup_for_relax_after_create_pre_ms_model
            if (s% mstar > s% initial_mass*Msun) then ! need to reduce mass
               write(*,1) 'reduce mass to', s% initial_mass
               call do_relax_mass(s% id, s% initial_mass, lg_max_abs_mdot, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in do_relax_mass'
                  return
               end if
            else if (s% mstar < s% initial_mass*Msun) then ! need to increase mass
               write(*,1) 'increase mass to', s% initial_mass
               call do_relax_mass_scale( &
                  s% id, s% initial_mass, s% job% dlgm_per_step, s% job% change_mass_years_for_dt, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in do_relax_mass'
                  return
               end if
            end if
            call do_relax_num_steps( &
               s% id, s% pre_ms_relax_num_steps, s% dt_next, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in do_relax_num_steps'
               return
            end if
            
            if (s% job% pre_ms_relax_to_start_radiative_core) then               
               call do_relax_to_radiative_core(s% id, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in do_relax_to_radiative_core'
                  return
               end if
            end if
            call done_relax_after_create_pre_ms_model
         else if (do_which == do_create_initial_model .and. &
                  s% initial_model_relax_num_steps > 0) then
            call do_relax_num_steps( &
               s% id, s% initial_model_relax_num_steps, s% dt_next, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in do_relax_num_steps'
               return
            end if
         end if

         s% doing_first_model_of_run = .true.
         s% num_trace_history_values = num_trace_history_values
         s% warning_limit_for_max_residual = warning_limit_for_max_residual

         contains
         
         subroutine setup_for_relax_after_create_pre_ms_model
            save_atm_option = s% atm_option
            save_atm_T_tau_relation = s% atm_T_tau_relation
            save_atm_T_tau_opacity = s% atm_T_tau_opacity
            save_Pextra_factor = s% Pextra_factor
            s% atm_option = 'T_tau'
            s% atm_T_tau_relation = 'Eddington'
            s% atm_T_tau_opacity = 'fixed'
            s% Pextra_factor = 2
         end subroutine setup_for_relax_after_create_pre_ms_model

         subroutine done_relax_after_create_pre_ms_model
            s% atm_option = save_atm_option
            s% atm_T_tau_relation = save_atm_T_tau_relation
            s% atm_T_tau_opacity = save_atm_T_tau_opacity
            s% Pextra_factor = save_Pextra_factor
         end subroutine done_relax_after_create_pre_ms_model

         subroutine check_initials
            include 'formats'
            if (abs(initial_mass - s% initial_mass) > &
                  1d-3*initial_mass .and. initial_mass > 0) then
               write(*,1) "WARNING -- inlist initial_mass ignored", initial_mass
               write(*,1) "using saved initial_mass instead", s% initial_mass
               write(*,*)
            end if
            if (abs(initial_z - s% initial_z) > &
                  1d-3*initial_z .and. initial_z > 0) then
               write(*,1) "WARNING -- inlist initial_z ignored", initial_z
               write(*,1) "using saved initial_z instead", s% initial_z
               write(*,*)
            end if
         end subroutine check_initials

      end subroutine model_builder


      logical function doing_restart(restart_filename)
         character (len=*) :: restart_filename
         inquire(file=restart_filename, exist=doing_restart)
      end function doing_restart


      integer function null_how_many_extra_header_items(id)
         integer, intent(in) :: id
         null_how_many_extra_header_items = 0
      end function null_how_many_extra_header_items


      subroutine null_data_for_extra_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=80) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_data_for_extra_header_items


      integer function null_how_many_binary_history_columns(binary_id)
         integer, intent(in) :: binary_id
         null_how_many_binary_history_columns = 0
      end function null_how_many_binary_history_columns


      subroutine null_data_for_binary_history_columns( &
            binary_id, n, names, vals, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id, n
         character (len=80) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_data_for_binary_history_columns
      
      
      integer function null_how_many_extra_binary_history_columns(binary_id)
         integer, intent(in) :: binary_id
         null_how_many_extra_binary_history_columns = 0
      end function null_how_many_extra_binary_history_columns


      subroutine null_data_for_extra_binary_history_columns( &
            binary_id, n, names, vals, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id, n
         character (len=80) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_data_for_extra_binary_history_columns


      subroutine do_garbage_collection(eosDT_cache_dir, ierr)
         use eos_lib, only: eos_init, eos_shutdown
         use eos_def, only: use_cache_for_eos
         integer, intent(inout) :: ierr
         character (len=*), intent(in) :: eosDT_cache_dir
         ! Remove existing eos data 
         call eos_shutdown()
         ! Re-initliaze eos
         call eos_init(eosDT_cache_dir,&
               ! After first init this is set in eos_def
               use_cache_for_eos,&
               ierr)
      end subroutine do_garbage_collection
      
      
      end module init
