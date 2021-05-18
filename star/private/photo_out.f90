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

      module photo_out

      use star_private_def

      implicit none

      private
      public :: output_star_photo

      contains


      subroutine output_star_photo(s,iounit,ierr)
         use rates_def,only:num_rvs
         use rsp_def, only: rsp_photo_out
         use alloc
         type (star_info), pointer :: s
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr

         integer :: part_number, nz, len_history_col_spec
         integer :: k

         include 'formats'

         ierr = 0
         nz = s% nz

         part_number = 0 ! part_numbers are just a consistency check on the data file

         write(iounit) star_def_version
         
         call write_part_number(iounit)
         write(iounit) &
            s% initial_z, & ! need this since read_model can change what is in the inlist
            s% total_num_solver_iterations, &
            s% nz, s% nvar_hydro, s% nvar_chem, s% nvar_total, &
            s% v_flag, s% u_flag, s% rotation_flag, s% RSP2_flag, s% RSP_flag, &
            s% RTI_flag, s% conv_vel_flag, s% w_div_wc_flag, s% j_rot_flag, s% D_omega_flag, s% am_nu_rot_flag, &
            s% have_mlt_vc, s% species, s% num_reactions, &
            s% model_number, s% star_mass, &
            s% mstar, s% xmstar, s% M_center, s% v_center, s% R_center, s% L_center, &
            s% time, s% dt, s% have_previous_conv_vel, &
            s% was_in_implicit_wind_limit, &
            s% using_revised_net_name, &
            s% revised_net_name, &
            s% using_revised_max_yr_dt, &
            s% revised_max_yr_dt, &
            s% astero_using_revised_max_yr_dt, &
            s% astero_revised_max_yr_dt, &
            s% using_RSP2, s% previous_step_was_using_RSP2, &
            s% cumulative_energy_error, s% cumulative_extra_heating, &
            s% have_initial_energy_integrals, s% total_energy_initial, &
            s% force_tau_factor, s% force_Tsurf_factor, s% force_opacity_factor

         write(iounit) s% net_name

         call write_part_number(iounit)
         write(iounit) &
            s% dq(1:nz), s% xa(:,1:nz), s% xh(:,1:nz), &
            s% omega(1:nz), s% j_rot(1:nz), s% mlt_vc(1:nz), s% conv_vel(1:nz)

         call write_part_number(iounit)
         write(iounit) &
            s% rsp_num_periods, s% rsp_dt, s% rsp_period, s% RSP_have_set_velocities, &
            s% dt_limit_ratio, s% tau_base

         write(iounit) &
            s% i_lnd, s% i_lnT, s% i_lnR, s% i_lum, s% i_Et_RSP, s% i_erad_RSP, s% i_Fr_RSP, &
            s% i_v, s% i_u, s% i_alpha_RTI, s% i_ln_cvpv0, s% i_w, s% i_Hp, s% i_w_div_wc, s% i_j_rot, &
            s% i_dv_dt, s% i_equL, s% i_dlnd_dt, s% i_dlnE_dt, &
            s% i_dEt_RSP_dt, s% i_derad_RSP_dt, s% i_dFr_RSP_dt, s% i_du_dt, s% i_dlnR_dt, &
            s% i_dln_cvpv0_dt, s% i_dalpha_RTI_dt, s% i_detrb_dt, s% i_equ_Hp

         write(iounit) &
            s% model_controls_filename, s% model_data_filename, &
            s% most_recent_profile_filename, s% most_recent_controls_filename, &
            s% most_recent_model_data_filename

         call write_part_number(iounit)
         write(iounit) &
            s% recent_log_header, s% phase_of_evolution, s% dt_next, s% dt_next_unclipped

         call write_part_number(iounit) 
         write(iounit) &
            s% num_skipped_setvars, s% num_retries, s% num_setvars, &  
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
            s% num_hydro_merges, s% num_hydro_splits, s% num_solver_setvars, &
            s% mesh_call_number, s% solver_call_number, s% diffusion_call_number, &
            s% gradT_excess_alpha, s% Teff, s% power_nuc_burn, s% power_h_burn, s% power_he_burn, s% power_z_burn, s% power_photo, &
            s% why_Tlim, s% dt_why_count(1:numTlim), s% dt_why_retry_count(1:numTlim), &
            s% timestep_hold, s% model_number_for_last_retry, s% model_number_for_last_retry_old, &
            s% init_model_number, s% most_recent_photo_name, &
            s% rand_i97, s% rand_j97, s% rand_u(1:rand_u_len), s% rand_c, s% rand_cd, s% rand_cm

         call write_part_number(iounit)
         write(iounit) s% len_extra_iwork, s% len_extra_work
         if (s% len_extra_iwork > 0) then
            write(iounit) s% extra_iwork(1:s% len_extra_iwork)
            !write(iounit) s% extra_iwork_old(1:s% len_extra_iwork)
         end if
         if (s% len_extra_work > 0) then
            write(iounit) s% extra_work(1:s% len_extra_work)
            !write(iounit) s% extra_work_old(1:s% len_extra_work)
         end if

         write(iounit) s% ixtra
         write(iounit) s% xtra
         write(iounit) s% lxtra

         if (associated(s% history_column_spec)) then
            len_history_col_spec = size(s% history_column_spec)
            write(iounit) len_history_col_spec
            write(iounit) s% history_column_spec(1:len_history_col_spec)
         else
            write(iounit) 0 ! len_log_col_spec
         end if

         write(iounit)  &
            s% number_of_history_columns, s% model_number_of_history_values, &
            s% need_to_set_history_names_etc
         if (s% number_of_history_columns > 0) then
            write(iounit) s% history_value_is_integer(1:s% number_of_history_columns)
            do k=1,s% number_of_history_columns
               write(iounit) s% history_names(k)
            end do
         end if
         
         if (s% rsp_flag) call rsp_photo_out(s, iounit)

         call write_part_number(iounit)

         call s% other_photo_write(s% id, iounit)

         call write_part_number(iounit)
         
         s% need_to_setvars = .true.

         contains

         subroutine write_part_number(iounit)
            integer, intent(in) :: iounit
            part_number = part_number + 1
            write(iounit) part_number
         end subroutine write_part_number


      end subroutine output_star_photo


      end module photo_out
