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

      module photo_in

      use star_private_def

      implicit none

      private
      public :: read_star_photo

      contains


      subroutine read_star_photo(s, fname, ierr)
         use utils_lib, only: integer_dict_define, integer_dict_create_hash
         use rates_def, only: num_rvs
         use alloc, only: set_var_info, allocate_star_info_arrays, &
            alloc_star_info_old_arrays
         use rsp_def, only: rsp_photo_in
         type (star_info), pointer :: s
         character (len=*), intent(in) :: fname
         integer, intent(out) :: ierr

         integer :: iounit, i, j, k, version, part_number, &
            len_history_col_spec, nz, nz_old
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         part_number = 0 ! part_numbers are just a consistency check

         write(*, *) 'read ', trim(fname)
         open(newunit=iounit, file=trim(fname), action='read', &
            status='old', iostat=ierr, form='unformatted')
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'Failed to open ', trim(fname)
            return
         end if

         read(iounit, iostat=ierr) version
         if (failed('version')) return
         if (version /= star_def_version) then
            write(*,'(/,a,/)') ' FAILURE: the restart data' // &
               ' is from a previous version of the code and is no longer usable.'
            ierr = -1
            return
         end if

         call read_part_number(iounit)
         if (failed('generations')) return

         read(iounit, iostat=ierr) &
            s% generations, s% total_num_solver_iterations, &
            s% nz, s% nz_old, &
            s% nvar_hydro, s% nvar_chem, s% nvar, s% init_model_number, &
            s% v_flag, s% u_flag, s% rotation_flag, s% Eturb_flag, &
            s% RTI_flag, s% conv_vel_flag, s% w_div_wc_flag, s% j_rot_flag, s% D_omega_flag, s% am_nu_rot_flag, &
            s% rsp_flag, s% prev_Lmax, s% species, s% num_reactions, &
            s% model_number, s% model_number_old, s% star_mass, &
            s% mstar, s% xmstar, s% M_center, s% v_center, s% R_center, s% L_center, &
            s% total_radiation, s% min_kap_floor, s% time, s% total_angular_momentum, &
            s% prev_create_atm_R0_div_R, s% dt, s% have_previous_conv_vel, &
            s% initial_mass, s% initial_z, s% initial_y, &
            s% was_in_implicit_wind_limit, &
            s% model_number_for_last_retry, &
            s% using_revised_net_name, &
            s% revised_net_name, &
            s% using_revised_max_yr_dt, &
            s% revised_max_yr_dt, &
            s% astero_using_revised_max_yr_dt, &
            s% astero_revised_max_yr_dt, &
            s% cumulative_energy_error, &
            s% force_tau_factor, s% force_Tsurf_factor, s% force_opacity_factor
         if (failed('initial_y')) return
         
         if (s% force_tau_factor > 0 .and. s% tau_factor /= s% force_tau_factor .and. &
               s% tau_factor /= s% job% set_to_this_tau_factor) then
            s% tau_factor = s% force_tau_factor
            write(*,1) 'set tau_factor to photo value', s% tau_factor
         end if

         if (s% force_Tsurf_factor > 0 .and. s% Tsurf_factor /= s% force_Tsurf_factor .and. &
               s% Tsurf_factor /= s% job% set_to_this_Tsurf_factor) then
            s% Tsurf_factor = s% force_Tsurf_factor
            write(*,1) 'set Tsurf_factor to photo value', s% Tsurf_factor
         end if

         if (s% force_opacity_factor > 0 .and. s% opacity_factor /= s% force_opacity_factor .and. &
               s% opacity_factor /= s% job% relax_to_this_opacity_factor) then
            s% opacity_factor = s% force_opacity_factor
            write(*,1) 'set opacity_factor to photo value', s% opacity_factor
         end if

         if (s% using_revised_net_name)  &
            s% net_name = s% revised_net_name

         if (s% using_revised_max_yr_dt) &
            s% max_years_for_timestep = s% revised_max_yr_dt

         if (s% astero_using_revised_max_yr_dt) &
            s% max_years_for_timestep = s% astero_revised_max_yr_dt

         nz = s% nz
         nz_old = s% nz_old

         read(iounit, iostat=ierr) s% net_name
         if (failed('nz_old')) return

         call set_var_info(s, ierr)
         if (failed('set_var_info')) return

         ! allocate everything
         call allocate_star_info_arrays(s, ierr)
         if (failed('allocate_star_info_arrays')) return

         call alloc_star_info_old_arrays(s, ierr)
         if (failed('alloc_star_info_old_arrays')) return

         call read_part_number(iounit)
         if (failed('dq')) return

         read(iounit, iostat=ierr) &
            s% dq(1:nz), s% xa(:,1:nz), s% xh(:,1:nz), &
            s% omega(1:nz), s% j_rot(1:nz)

!         read(iounit, iostat=ierr) &
!            s% q(1:nz), s% m(1:nz), s% dm(1:nz), s% dm_bar(1:nz), &
!            s% w_div_w_crit_roche(1:nz), &
!            s% D_omega(1:nz), s% am_nu_rot(1:nz), &
!            s% dlnd_dt(1:nz), s% dlnT_dt(1:nz), &
!            s% eps_grav(1:nz), s% conv_vel(1:nz), s% lnT(1:nz)
            
!         call read_part_number(iounit)
!         if (failed('*_old')) return
!         if (s% generations > 1 .and. .not. s% rsp_flag) then
!            read(iounit, iostat=ierr) &
!               s% omega_old(1:nz_old), &
!               s% j_rot_old(1:nz_old), &
!               s% dq_old(1:nz_old), &
!               s% xh_old(:,1:nz_old), &
!               s% xa_old(:,1:nz_old)
!            if (failed('xh_old')) return
!         end if

         call read_part_number(iounit)
         if (failed('mstar_dot')) return

         read(iounit, iostat=ierr) &
            s% rsp_num_periods, s% rsp_dt, s% rsp_period, s% RSP_have_set_velocities, &
            s% mstar_dot, s% v_surf, s% gradT_excess_alpha, s% dt_limit_ratio, &
            s% L_phot, s% T_surf, s% P_surf, s% L_surf, s% h1_czb_mass, &
            s% he_core_mass, s% c_core_mass, s% tau_base, s% Teff, &
            s% center_eps_nuc, s% Lrad_div_Ledd_avg_surf, s% w_div_w_crit_avg_surf, &
            s% n_conv_regions, s% cz_bot_mass(:), s% cz_top_mass(:)
         if (failed('cz_top_mass_old')) return

         read(iounit, iostat=ierr) &
            s% i_lnd, s% i_lnT, s% i_lnR, s% i_lum, s% i_eturb_RSP, s% i_erad_RSP, s% i_Fr_RSP, &
            s% i_v, s% i_u, s% i_alpha_RTI, s% i_ln_cvpv0, s% i_eturb, &
            s% i_w_div_wc, s% i_j_rot, &
            s% i_dv_dt, s% i_equL, s% i_dlnd_dt, s% i_dlnE_dt, &
            s% i_deturb_RSP_dt, s% i_derad_RSP_dt, s% i_dFr_RSP_dt, &
            s% i_du_dt, s% i_dlnR_dt, s% i_dln_cvpv0_dt, s% i_dalpha_RTI_dt, s% i_deturb_dt
         if (failed('i_dalpha_RTI_dt')) return

         read(iounit, iostat=ierr) &
            s% model_profile_filename, s% model_controls_filename, s% model_data_filename, &
            s% most_recent_profile_filename, s% most_recent_controls_filename, &
            s% most_recent_model_data_filename
         if (failed('most_recent_model_data_filename')) return

         call read_part_number(iounit)
         if (failed('recent_log_header')) return

         read(iounit, iostat=ierr) &
            s% recent_log_header, s% phase_of_evolution, &
            s% dt_next, s% dt_next_unclipped !, &
!            s% eps_nuc(1:nz), &
!            s% d_epsnuc_dlnd(1:nz), &
!            s% d_epsnuc_dlnT(1:nz), &
!            s% d_epsnuc_dx(:,1:nz), &
!            s% eps_nuc_categories(:,1:nz), &
!            s% dxdt_nuc(:,1:nz), &
!            s% d_dxdt_nuc_dRho(:,1:nz), &
!            s% d_dxdt_nuc_dT(:,1:nz), &
!            s% d_dxdt_nuc_dx(:,:,1:nz), &
!            s% eps_nuc_neu_total(1:nz)
         if (failed('eps_nuc_neu_total')) return
         
         if (s% dt_next <= 0d0) s% dt_next = s% dt_next_unclipped

         call read_part_number(iounit)
         if (failed('read_part_number')) return

         read(iounit, iostat=ierr) &
            s% num_solver_iterations, s% num_skipped_setvars, s% num_retries, s% num_setvars, &  
                   
            s% total_num_solver_iterations, &
            s% total_num_solver_relax_iterations, &
            s% total_num_solver_calls_made, &
            
            s% total_num_solver_relax_calls_made, &
            s% total_num_solver_calls_converged, &
            s% total_num_solver_relax_calls_converged, &

            s% total_step_attempts, s% total_relax_step_attempts, &
            s% total_step_retries, s% total_relax_step_retries, &
            s% total_step_redos, s% total_relax_step_redos, &
            s% total_steps_finished, s% total_relax_steps_finished, &

            s% num_hydro_merges, s% num_hydro_splits, &
            s% initial_L_center, s% initial_R_center, s% initial_v_center, s% tau_center, &
            s% timestep_hold, s% model_number_for_last_retry, s% bad_max_corr_cnt, &
            s% mesh_call_number, s% solver_call_number, s% diffusion_call_number, &
            s% Tlim_dXnuc_species, s% Tlim_dXnuc_cell, &
            s% Tlim_dXnuc_drop_species, s% Tlim_dXnuc_drop_cell, &
            s% Tlim_dX_species, s% Tlim_dX_cell, &
            s% Tlim_dt_div_min_dr_div_cs_cell, &
            s% Tlim_dX_div_X_species, s% Tlim_dX_div_X_cell, &
            s% Tlim_dlgL_nuc_category, s% Tlim_dlgL_nuc_cell, &
            s% why_Tlim, s% dt_why_retry_count, s% dt_why_count, &
            s% initial_timestep, s% result_reason, s% need_to_update_history_now, &
            s% dt_why_count(1:numTlim), s% dt_why_retry_count(1:numTlim), &
            s% need_to_save_profiles_now, s% save_profiles_model_priority, &
            s% doing_flash_wind, s% doing_rlo_wind, s% doing_nova_wind, s% most_recent_photo_name, &
            s% rand_i97, s% rand_j97, s% rand_u(1:rand_u_len), s% rand_c, s% rand_cd, s% rand_cm
         if (failed('most_recent_photo_name')) return

         call read_part_number(iounit)
         if (failed('len_extra_iwork')) return

         read(iounit, iostat=ierr) s% len_extra_iwork, s% len_extra_work
         if (failed('len_extra_work')) return

         if (s% len_extra_iwork > 0) then
            allocate( &
               s% extra_iwork(s% len_extra_iwork), &
               s% extra_iwork_old(s% len_extra_iwork), &
               stat=ierr)
            if (failed('allocate extra_iwork')) return
            read(iounit, iostat=ierr) s% extra_iwork(1:s% len_extra_iwork)
            if (failed('read extra_iwork')) return
            read(iounit, iostat=ierr) s% extra_iwork_old(1:s% len_extra_iwork)
            if (failed('allocate extra_iwork_old')) return
         else
            nullify(s% extra_iwork, s% extra_iwork_old)
         end if

         if (s% len_extra_work > 0) then
            allocate( &
               s% extra_work(s% len_extra_work), &
               s% extra_work_old(s% len_extra_work), &
               stat=ierr)
            if (failed('allocate extra_work')) return
            read(iounit, iostat=ierr) s% extra_work(1:s% len_extra_work)
            if (failed('read extra_work')) return
            read(iounit, iostat=ierr) s% extra_work_old(1:s% len_extra_work)
            if (failed('read extra_work_old')) return
         else
            nullify(s% extra_work, s% extra_work_old)
         end if

         read(iounit, iostat=ierr) s% ixtra
         if (failed('ixtra')) return
         read(iounit, iostat=ierr) s% xtra
         if (failed('xtra')) return
         read(iounit, iostat=ierr) s% lxtra
         if (failed('lxtra')) return
         
         read(iounit, iostat=ierr) s% ixtra_old
         if (failed('ixtra_old')) return
         read(iounit, iostat=ierr) s% xtra_old
         if (failed('xtra_old')) return
         read(iounit, iostat=ierr) s% lxtra_old
         if (failed('lxtra_old')) return

         read(iounit, iostat=ierr) len_history_col_spec
         if (failed('len_history_col_spec')) return
         if (len_history_col_spec > 0) then
            allocate(s% history_column_spec(len_history_col_spec), stat=ierr)
            if (failed('alloc history_column_spec')) return
            read(iounit, iostat=ierr) s% history_column_spec(1:len_history_col_spec)
            if (failed('read history_column_spec')) return
         end if

         read(iounit, iostat=ierr) &
            s% number_of_history_columns, s% model_number_of_history_values, &
            s% need_to_set_history_names_etc
         if (failed('number_of_history_columns')) return

         if (s% number_of_history_columns > 0) then

            allocate(s% history_value_is_integer(s% number_of_history_columns), stat=ierr)
            if (failed('alloc history_value_is_integer')) return
            read(iounit, iostat=ierr) s% history_value_is_integer(1:s% number_of_history_columns)
            if (failed('read history_value_is_integer')) return

            allocate(s% history_names(s% number_of_history_columns), stat=ierr)
            if (failed('alloc history_names')) return
            do k=1,s% number_of_history_columns
               read(iounit, iostat=ierr) s% history_names(k)
               if (failed('read history_names')) return
            end do

            ! rebuild the history_names_dict
            do j = 1, s% number_of_history_columns
               call integer_dict_define(s% history_names_dict, s% history_names(j), j, ierr)
               if (failed('integer_dict_define history_names_dict')) return
            end do
            call integer_dict_create_hash(s% history_names_dict, ierr)
            if (failed('integer_dict_create_hash history_names_dict')) return

         end if
         
         if (s% rsp_flag) then
            call rsp_photo_in(s, iounit, ierr)
            if (failed('after rsp_photo_in')) return
         end if

         call read_part_number(iounit)
         if (failed('before other_photo_read')) return

         call s% other_photo_read(s% id, iounit, ierr)
         if (failed('after other_photo_read')) return

         call read_part_number(iounit)
         if (failed('final read_part_number')) return
         
         s% need_to_setvars = .true. ! set this after photo out or photo in

         close(iounit)

         contains

         subroutine read_part_number(iounit)
            integer, intent(in) :: iounit
            integer :: i
            part_number = part_number + 1
            read(iounit, iostat=ierr) i
            if (ierr /= 0) return
            if (i /= part_number) ierr = -1
            !write(*,*) 'part_number', part_number
         end subroutine read_part_number

         logical function failed(str)
            character (len=*), intent(in) :: str
            i = i+1
            if (ierr /= 0) then
               write(*, *) 'read_star_photo failed for ' // trim(str)
               failed = .true.
               return
            end if
            failed = .false.
         end function failed


      end subroutine read_star_photo


      end module photo_in
