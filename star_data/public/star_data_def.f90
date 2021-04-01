! ***********************************************************************
!
!   Copyright (C) 2019  Bill Paxton and MESA Team
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
!
! ***********************************************************************

      module star_data_def

      use rates_def, only: rates_reaction_id_max, other_screening_interface
      use utils_def, only: integer_dict
      use chem_def, only: num_categories, iso_name_length
      use const_def, only: sp, dp, qp, strlen
      use rates_def, only: maxlen_reaction_Name
      use eos_def, only: EoS_General_Info
      use kap_def, only: Kap_General_Info
      use net_def, only: Net_General_Info
      use colors_def, only:  max_num_color_files, max_num_bcs_per_file
      use auto_diff, only: auto_diff_real_star_order1
      
      implicit none      
      
      include "star_data_def.inc"
      
      integer, parameter :: max_extras_params = 20, max_extras_cpar_len = strlen
      integer, parameter :: max_num_special_rate_factors = 20

      integer, parameter :: star_num_xtra_vals = 30

      type star_job_controls
         include "star_job_controls.inc"
         include "star_job_controls_dev.inc"
         real(dp) :: &
             step_loop_timing, after_step_timing, before_step_timing, &
             check_time_start, check_time_end, elapsed_time, &
             check_step_loop_timing, check_after_step_timing, check_before_step_timing
         integer(8) :: time0, time1, clock_rate, time0_extra, time1_extra, time0_initial
      end type star_job_controls
      

      type star_info
         
         include "star_data.inc"
   
         ! handles
            integer :: eos_handle
            integer :: kap_handle
            integer :: net_handle
               
         ! star id
            integer :: id ! unique identifier for each star_info instance
            
         ! Name of the main inlist used
            character (len=strlen) :: inlist_fname
         
         ! private
            logical :: in_use
            logical :: do_burn, do_mix
            logical :: used_extra_iter_in_solver_for_accretion
            integer :: retry_cnt, redo_cnt
            type (EoS_General_Info), pointer :: eos_rq ! from call eos_ptr(s% eos_handle,s% eos_rq,ierr)
            type (Kap_General_Info), pointer :: kap_rq ! from call kap_ptr(s% kap_handle,s% kap_rq,ierr)
            type (Net_General_Info), pointer :: net_rq ! from call net_ptr(s% net_handle,s% net_rq, ierr)
            
            ! parameters for create pre ms -- set in run_star before calling star_create_pre_ms_model
            real(dp) :: pre_ms_T_c, pre_ms_guess_rho_c, &
               pre_ms_d_log10_P, pre_ms_logT_surf_limit, pre_ms_logP_surf_limit
            integer :: pre_ms_initial_zfracs, pre_ms_relax_num_steps
            logical :: pre_ms_change_net, pre_ms_dump_missing_heaviest
            character (len=net_name_len) :: pre_ms_new_net_name
            
            ! parameters for create initial model
            real(dp) :: & 
               radius_in_cm_for_create_initial_model, &
               mass_in_gm_for_create_initial_model, &
               center_logP_1st_try_for_create_initial_model, &
               entropy_1st_try_for_create_initial_model, &
               abs_e01_tolerance_for_create_initial_model, &
               abs_e02_tolerance_for_create_initial_model           
            integer :: initial_zfracs_for_create_initial_model, &
               max_tries_for_create_initial_model
            integer :: initial_model_relax_num_steps
            real(dp) :: initial_model_eps
            logical :: initial_model_change_net, initial_dump_missing_heaviest
            character (len=net_name_len) :: initial_model_new_net_name
         
            ! extra profile entries for developer debugging
            real(dp), dimension(:,:), pointer :: profile_extra ! (nz,max_num_profile_extras)
            character (len=64) :: profile_extra_name(max_num_profile_extras)
            
         ! controls
            type (star_job_controls) :: job ! separate type to avoid name clashes
            include "star_controls.inc"
            include "star_controls_dev.inc"
            include "pgstar_controls.inc"
         
      end type star_info

      logical :: have_initialized_star_handles = .false.
      integer, parameter :: max_star_handles = 10 ! this can be increased as necessary
      type (star_info), target, save :: star_handles(max_star_handles) 
         ! gfortran seems to require "save" here.  at least it did once upon a time.


      contains


      subroutine star_ptr(id, s, ierr)
         integer, intent(in) :: id
         type (star_info), pointer, intent(inout) :: s
         integer, intent(out) :: ierr
         call get_star_ptr(id, s, ierr)
      end subroutine star_ptr


      subroutine get_star_ptr(id,s,ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr         
         if (id < 1 .or. id > max_star_handles) then
            ierr = -1
            return
         end if
         s => star_handles(id)
         ierr = 0
      end subroutine get_star_ptr
      
      
      subroutine result_reason_init         
         result_reason_str(result_reason_normal) = 'normal'
         result_reason_str(dt_is_zero) = 'dt_is_zero'
         result_reason_str(nonzero_ierr) = 'nonzero_ierr'
         result_reason_str(hydro_failed_to_converge) = 'hydro_failed'
         result_reason_str(do_burn_failed) = 'do_burn_failed'
         result_reason_str(diffusion_failed) = 'element_diffusion_failed'
         result_reason_str(too_many_steps_for_burn) = 'too_many_steps_for_burn'
         result_reason_str(too_many_steps_for_diffusion) = 'too_many_steps_for_diffusion'
         result_reason_str(too_many_steps_for_hydro) = 'too_many_steps_for_hydro'
         result_reason_str(adjust_mesh_failed) = 'adjust_mesh_failed'
         result_reason_str(adjust_mass_failed) = 'adjust_mass_failed'
         result_reason_str(core_dump_model_number) = 'core_dump_model_number'
         result_reason_str(timestep_limits) = 'cannot find acceptable model'
         result_reason_str(variable_change_limits) = 'variable_change_limits'
         result_reason_str(explicit_hydro_failed) = 'explicit_hydro_failed'
         result_reason_str(abs_rel_run_E_err) = 'abs_rel_run_E_err'
         result_reason_str(forced_stop) = 'forced_stop'
      end subroutine result_reason_init

      
      subroutine do_star_def_init(mesa_dir_init, ierr)
         character (len=*), intent(in) :: mesa_dir_init
         integer, intent(out) :: ierr
         ierr = 0
         call result_reason_init
      end subroutine do_star_def_init


      end module star_data_def

