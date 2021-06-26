! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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
 
      module run_star_support

      use star_lib
      use star_def
      use chem_def
      use chem_lib
      use const_def
      use math_lib
      use eos_lib
      use kap_def
      use net_def
      use net_lib
      use other_extras
      use rates_lib, only: set_which_rate_1212

      implicit none
      
      integer :: id_from_read_star_job = 0

      ! Set MESA_INLIST_RESOLVED to true when you no longer want the routine
      ! resolve_inlist_fname to look at the MESA_INLIST environment variable
      logical :: MESA_INLIST_RESOLVED = .false.

      private
      public :: do_read_star_job, do_read_star_job_and_return_id
      public :: run1_star
      public :: start_run1_star
      public :: do_evolve_one_step
      public :: after_evolve_loop
      public :: failed
      public :: id_from_read_star_job
      public :: MESA_INLIST_RESOLVED
      
      ! deprecated, but kept around for use by binary
      public :: before_evolve_loop, after_step_loop, before_step_loop, do_saves, &
         resolve_inlist_fname, terminate_normal_evolve_loop
      
      contains 
            

      subroutine run1_star( &
            do_alloc_star, do_free_star, okay_to_restart, &
            id, restart, &
            extras_controls, &
            ierr, &
            inlist_fname_arg)
         
         logical, intent(in) :: do_alloc_star, do_free_star, okay_to_restart
         integer, intent(inout) :: id ! input if not do_alloc_star
         logical, intent(inout) :: restart ! input if not do_alloc_star
         character (len=*) :: inlist_fname_arg
         integer, intent(out) :: ierr
         optional inlist_fname_arg
         
         interface

            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      
     
         end interface
         
         logical :: continue_evolve_loop
         type (star_info), pointer :: s
         character (len=strlen) :: restart_filename
            
         logical, parameter :: pgstar_ok = .true.
         logical, parameter :: dbg = .false.
         
         1 format(a35, 99(1pe26.16))
         2 format(a55, i7, 1pe26.16)
         3 format(a15, 2x, f15.6)
         4 format(a15, 2x, e15.6)

         11 format(a35, f20.10)

         restart_filename = 'restart_photo'
         call start_run1_star( &
            do_alloc_star, do_free_star, okay_to_restart, &
            id, restart, restart_filename, pgstar_ok, dbg, &
            extras_controls, ierr, inlist_fname_arg)
         if (failed('do_before_evolve_loop',ierr)) return

         call star_ptr(id, s, ierr)
         if (failed('star_ptr',ierr)) return

         continue_evolve_loop = .true.

         if (dbg) write(*,*) 'start evolve_loop'
         evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop
         
            continue_evolve_loop = do_evolve_one_step(s, dbg, ierr)
            if (failed('do_evolve_one_step',ierr)) return

         end do evolve_loop

         call after_evolve_loop(s% id, do_free_star, ierr)
         if (failed('after_evolve_loop',ierr)) return

      end subroutine run1_star  
      
      
      subroutine start_run1_star( &
            do_alloc_star, do_free_star, okay_to_restart, &
            id, restart, restart_filename, pgstar_ok, dbg, &
            extras_controls, ierr, inlist_fname_arg)
            
         logical, intent(in) :: do_alloc_star, do_free_star, okay_to_restart
         integer, intent(inout) :: id ! input if not do_alloc_star
         logical, intent(inout) :: restart ! input if not do_alloc_star
         logical, intent(in) :: pgstar_ok, dbg
         character (len=*) :: restart_filename, inlist_fname_arg
         optional inlist_fname_arg
         integer, intent(out) :: ierr
         
         interface

            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      
     
         end interface

         type (star_info), pointer :: s
         character (len=strlen) :: inlist_fname
         
         include 'formats'

         ierr = 0

         call resolve_inlist_fname(inlist_fname,inlist_fname_arg)

         ! star is initialized here
         call do_before_evolve_loop( &
              do_alloc_star, okay_to_restart, restart, pgstar_ok, &
              null_binary_controls, extras_controls, &
              id_from_read_star_job, inlist_fname, restart_filename, &
              dbg, 0, id, ierr)
         if (failed('do_before_evolve_loop',ierr)) return

         call star_ptr(id, s, ierr)
         if (failed('star_ptr',ierr)) return

         s% doing_timing = .false.
         s% job% check_before_step_timing = 0
         s% job% check_step_loop_timing = 0
         s% job% check_after_step_timing = 0
         s% job% time0_initial = 0

      end subroutine start_run1_star
      
      
      logical function do_evolve_one_step(s, dbg, ierr) result(continue_evolve_loop)
         type (star_info), pointer :: s
         logical, intent(in) :: dbg
         integer, intent(out) :: ierr
         
         logical :: first_try
         integer :: id
         integer :: result, model_number, source, nsteps, j, ci, nz
         
         include 'formats'
         
         ierr = 0
         id = s% id
         continue_evolve_loop = .true.

         call before_step_loop(s% id, ierr)
         if (failed('before_step_loop',ierr)) return

         result = s% extras_start_step(id)  
         if (result /= keep_going) then
            continue_evolve_loop = .false.
            return 
         end if        

         first_try = .true.
      
         step_loop: do ! may need to repeat this loop
         
            if (stop_is_requested(s)) then
               continue_evolve_loop = .false.
               result = terminate
               exit
            end if
         
            result = star_evolve_step(id, first_try)
            if (result == keep_going) result = star_check_model(id)
            if (result == keep_going) result = s% extras_check_model(id)
            if (result == keep_going) result = star_pick_next_timestep(id)            
            if (result == keep_going) exit step_loop
            
            model_number = get_model_number(id, ierr)
            if (failed('get_model_number',ierr)) return
                           
            if (result == retry .and. s% job% report_retries) then
               write(*,'(i6,3x,a,/)') model_number, &
                  'retry reason ' // trim(result_reason_str(s% result_reason))
            end if
            
            if (result == redo) then
               result = star_prepare_to_redo(id)
            end if
            if (result == retry) then
               result = star_prepare_to_retry(id)
            end if
            if (result == terminate) then
               continue_evolve_loop = .false.
               exit step_loop
            end if
            first_try = .false.
            
         end do step_loop
         
         ! once we get here, the only options are keep_going or terminate.
         ! redo or retry must be done inside the step_loop
         
         call after_step_loop(s% id, s% inlist_fname, &
             dbg, result, ierr)
         if (failed('after_step_loop',ierr)) return
            
         if (result /= keep_going) then
            if (result /= terminate) then
               write(*,2) 'ERROR in result value in run_star_extras: model', &
                  s% model_number
               write(*,2) 'extras_finish_step must return keep_going or terminate'
               write(*,2) 'result', result
               continue_evolve_loop = .false.
               return 
            end if
            if (s% result_reason == result_reason_normal) then
               call terminate_normal_evolve_loop(s% id, &
                  dbg, result, ierr)
               if (failed('terminate_normal_evolve_loop',ierr)) return
            end if
            continue_evolve_loop = .false.
            return 
         end if
         
         call do_saves(id, ierr)
         if (failed('do_saves',ierr)) return

         if (s% doing_timing) then
            call system_clock(s% job% time1_extra,s% job% clock_rate)
            s% job% after_step_timing = s% job% after_step_timing + &
               dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate
            s% job% check_time_end = eval_total_times(s% id, ierr)
            s% job% check_after_step_timing = s% job% check_after_step_timing + &
               (s% job% check_time_end - s% job% check_time_start)
         end if
         
      end function do_evolve_one_step
            

      subroutine null_binary_controls(id, binary_id, ierr)
         integer, intent(in) :: id, binary_id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_binary_controls


      ! Binary requires to set some controls in here, which is why
      ! binary_controls and binary_id are arguments. These do nothing
      ! for the case of single star evolution.
      subroutine before_evolve_loop( &
              do_alloc_star, okay_to_restart, restart, &
              binary_controls, extras_controls, &
              id_from_read_star_job, inlist_fname, restart_filename, &
              dbg, binary_id, id, ierr)
         logical, intent(in) :: do_alloc_star, okay_to_restart
         logical :: restart
         interface
            subroutine binary_controls(id, binary_id, ierr)
               integer, intent(in) :: id, binary_id
               integer, intent(out) :: ierr
            end subroutine binary_controls     
            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      
         end interface
         integer :: id_from_read_star_job
         character (len=*) :: inlist_fname, restart_filename
         character (len=512) :: temp_fname
         logical, intent(in) :: dbg
         integer, intent(in) :: binary_id
         integer, intent(out) :: id, ierr
         call do_before_evolve_loop( &
              do_alloc_star, okay_to_restart, restart, .true., &
              binary_controls, extras_controls, &
              id_from_read_star_job, inlist_fname, restart_filename, &
              dbg, binary_id, id, ierr)
      end subroutine before_evolve_loop
      
      
      subroutine do_before_evolve_loop( &
              do_alloc_star, okay_to_restart, restart, pgstar_ok, &
              binary_controls, extras_controls, &
              id_from_read_star_job, inlist_fname, restart_filename, &
              dbg, binary_id, id, ierr)
         logical, intent(in) :: do_alloc_star, okay_to_restart, pgstar_ok
         logical :: restart
         interface
            subroutine binary_controls(id, binary_id, ierr)
               integer, intent(in) :: id, binary_id
               integer, intent(out) :: ierr
            end subroutine binary_controls     
            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      
         end interface
         integer :: id_from_read_star_job
         character (len=*) :: inlist_fname, restart_filename
         character (len=512) :: temp_fname
         logical, intent(in) :: dbg
         integer, intent(in) :: binary_id
         integer, intent(out) :: id, ierr

         type (star_info), pointer :: s         
         
         include 'formats'
         
         if (do_alloc_star) then           
            if (id_from_read_star_job /= 0) then 
               ! already allocated by read_star_job
               id = id_from_read_star_job
               id_from_read_star_job = 0
            else
               call alloc_star(id, ierr)
               if (failed('alloc_star',ierr)) return
            end if         
            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return
         else
            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return
            call init_starting_star_data(s, ierr)
            if (failed('init_starting_star_data',ierr)) return
         end if
         
         s% inlist_fname = inlist_fname
         
         if (dbg) write(*,*) 'call starlib_init'
         call starlib_init(s, ierr) ! okay to do extra calls on this
         if (failed('star_init',ierr)) return
         
         if (dbg) write(*,*) 'call star_set_kap_and_eos_handles'
         call star_set_kap_and_eos_handles(id, ierr)
         if (failed('set_star_kap_and_eos_handles',ierr)) return
         
         if (dbg) write(*,*) 'call star_setup'
         call star_setup(id, inlist_fname, ierr)
         if (failed('star_setup',ierr)) return
         
         if(dbg) write(*,*) 'call add_fpe_checks'
         call add_fpe_checks(id, s, ierr)
         if (failed('add_fpe_checks',ierr)) return
         
         if(dbg) write(*,*) 'call multiply_tolerances'
         call multiply_tolerances(id, s, ierr)
         if (failed('multiply_tolerances',ierr)) return

         if(dbg) write(*,*) 'call pgstar_env_check'
         call pgstar_env_check(id, s, ierr)
         if (failed('pgstar_env_check',ierr)) return        

         ! testing module-level (atm/eos/kap/net) partials requires single-threaded execution
         if (s% solver_test_atm_partials .or. s% solver_test_eos_partials .or. &
               s% solver_test_kap_partials .or. s% solver_test_net_partials) then
            if (s% solver_test_partials_k > 0 .and. s% solver_test_partials_dx_0 > 0) then
               write(*,*) 'Forcing single-thread mode for testing of module-level partials'
               call omp_set_num_threads(1)
            end if
         end if
         
         if (len_trim(s% op_mono_data_path) == 0) &
            call get_environment_variable( &
               "MESA_OP_MONO_DATA_PATH", s% op_mono_data_path)
         
         if (len_trim(s% op_mono_data_cache_filename) == 0) &
            call get_environment_variable( &
               "MESA_OP_MONO_DATA_CACHE_FILENAME", s% op_mono_data_cache_filename)         
         if (restart_filename /= "restart_photo") then
            temp_fname  = trim(s% photo_directory) // '/' // trim(restart_filename)
            restart_filename  = trim(temp_fname)
         end if

         if (okay_to_restart) then
            restart = doing_a_restart(restart_filename)
         else
            restart = .false.
         end if
         
         if (s% job% show_log_description_at_start .and. .not. restart) then
            write(*,*)
            call show_log_description(id, ierr)
            if (failed('show_log_description',ierr)) return
         end if

         s% extras_startup => null_extras_startup
         s% extras_check_model => null_extras_check_model
         s% extras_start_step => null_extras_start_step
         s% extras_finish_step => null_extras_finish_step
         s% extras_after_evolve => null_extras_after_evolve
         s% how_many_extra_history_columns => null_how_many_extra_history_columns
         s% data_for_extra_history_columns => null_data_for_extra_history_columns
         s% how_many_extra_profile_columns => null_how_many_extra_profile_columns
         s% data_for_extra_profile_columns => null_data_for_extra_profile_columns

         if (dbg) write(*,*) 'call extras_controls'
         call extras_controls(id, ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call binary_controls'
         call binary_controls(id, binary_id, ierr)
         if (ierr /= 0) return
         
         if (dbg) write(*,*) 'call do_star_job_controls_before'
         call do_star_job_controls_before(id, s, restart, ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call do_load1_star'
         call do_load1_star(id, s, restart, restart_filename, ierr)
         if (failed('do_load1_star',ierr)) return
         
         if (dbg) write(*,*) 'call do_star_job_controls_after'
         call do_star_job_controls_after(id, s, restart, pgstar_ok, ierr)
         if (failed('do_star_job_controls_after',ierr)) return

         write(*,*)
         write(*,*)
         
         if (.not. restart) then
            if (dbg) write(*,*) 'call before_evolve'
            call before_evolve(id, ierr)
            if (failed('before_evolve',ierr)) return
         else
            call show_terminal_header(id, ierr)
            if (failed('show_terminal_header',ierr)) return
         end if
         
         if (dbg) write(*,*) 'call extras_startup'
         call s% extras_startup(id, restart, ierr)
         if (failed('extras_startup',ierr)) return

         if (s% job% profile_starting_model .and. .not. restart) then
            call star_set_vars(id, 0d0, ierr)
            if (failed('star_set_vars',ierr)) return
            write(*, '(a, i12)') 'save profile for model number ', s% model_number            
            call save_profile(id,3,ierr)
            if (failed('save_profile',ierr)) return
         end if

         if (s% model_number == s% job% save_model_number) then
            call star_set_vars(id, 0d0, ierr)
            if (failed('star_set_vars',ierr)) return
            write(*, '(a, i12)') 'write initial model ', s% model_number            
            call star_write_model(id, 'initial.mod', ierr)
            if (failed('star_write_model',ierr)) return
            write(*, *) 'saved to ' // 'initial.mod' ! trim(s% job% save_model_filename)
         end if
         
         if (len_trim(s% job% echo_at_start) > 0) then
            write(*,*)
            write(*,'(a)') trim(s% job% echo_at_start)
            write(*,*)
         end if

      end subroutine do_before_evolve_loop


      subroutine before_step_loop(id, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s         
         integer, intent(out) :: ierr
         integer :: k, model_number, j
         real(dp) :: gamma1_integral, integral_norm
         integer :: num_DT, num_PT, num_FreeEOS
         
         1 format(a35, 99(1pe26.16))
         2 format(a35, i7, 1pe26.16)
         3 format(a15, 2x, f15.6)
         4 format(a15, 2x, e15.6)

         11 format(a35, f20.10)

         call star_ptr(id, s, ierr)
         if (ierr/=0) return

         s% result_reason = result_reason_normal
         
         if (s% job% first_model_for_timing >= 0 .and. &
               s% model_number >= s% job% first_model_for_timing .and. &
               .not. s% doing_timing) then
            s% doing_timing = .true.
            write(*,*) 'start timing', s% model_number
            write(*,*)
            call system_clock(s% job% time0, s% job% clock_rate)
            s% job% time0_initial = s% job% time0
            s% job% step_loop_timing = 0
            s% job% after_step_timing = 0
            s% job% before_step_timing = 0
         end if
         
         if (s% doing_timing) then
            call system_clock(s% job% time0_extra,s% job% clock_rate)
            s% job% check_time_start = eval_total_times(s% id, ierr)
         end if
         
         if(s% job% num_steps_for_garbage_collection > 0 .and. s% model_number > 1) then
            if(mod(s% model_number, s% job% num_steps_for_garbage_collection) == 0)then
               if (s% job% report_garbage_collection) then
                  call num_eos_files_loaded( &
                     num_DT, num_FreeEOS)
                  write(*,*) "Start garbage collection model_number", s%model_number,"num eosDT", num_DT, &
                              "num FreeEOS",num_FreeEOS
               end if
               call star_do_garbage_collection(s% id,ierr)
               if (failed('star_do_garbage_collection',ierr)) return
            end if
            
            ! If reporting, we want to look at the step and the next step (to see the difference)
            if(mod(s% model_number-1, s% job% num_steps_for_garbage_collection) == 0 &
                  .and. s% job% report_garbage_collection)then
                  call num_eos_files_loaded( &
                     num_DT, num_FreeEOS)
                  write(*,*) "End garbage collection model_number  ", s%model_number,"num eosDT", num_DT, &
                              "num FreeEOS",num_FreeEOS
            end if
         end if
         
         if (s% job% enable_adaptive_network) then
            call star_adjust_net(s% id, &
               s% job% min_x_for_keep, &
               s% job% min_x_for_n, &
               s% job% min_x_for_add, &
               s% job% max_Z_for_add, &
               s% job% max_N_for_add, &
               s% job% max_A_for_add, &
               ierr)
            if (failed('star_adjust_net',ierr)) return
         end if
         
         if (s% job% auto_extend_net) then
            call extend_net(s, ierr)
            if (failed('extend_net',ierr)) return
         end if
         
         if (s% use_other_remove_surface) then
            call s% other_remove_surface(id, ierr, j)
            if (failed('other_remove_surface',ierr)) return
            if (j > 0) then
               call star_remove_surface_at_cell_k(s% id, j, ierr)
            end if
            call do_remove_surface(id, s, ierr)
         else
            call do_remove_surface(id, s, ierr)
            if (failed('do_remove_surface',ierr)) return
         end if
         
         if (s% job% remove_fallback_at_each_step) then
            call star_remove_fallback(id,ierr)
            if (failed('star_remove_fallback',ierr)) return
         end if
         
         if (s% job% limit_center_logP_at_each_step > -1d90) then
            call star_limit_center_logP( &
               id, s% job% limit_center_logP_at_each_step, ierr)
            if (failed('star_limit_center_logP',ierr)) return
         end if
         
         if (s% job% remove_center_logRho_limit > -1d90) then
            call star_remove_center_by_logRho( &
               id, s% job% remove_center_logRho_limit, ierr)
            if (failed('star_remove_center_by_logRho',ierr)) return
         end if
         
         if (s% center_ye <= s% job% center_ye_limit_for_v_flag &
               .and. (.not. s% v_flag) .and. (.not. s% u_flag)) then
            write(*,1) 'have reached center ye limit', &
               s% center_ye, s% job% center_ye_limit_for_v_flag
            write(*,1) 'set v_flag true'
            call star_set_v_flag(id, .true., ierr)
            if (failed('star_set_v_flag',ierr)) return
            if (ierr /= 0) return
         end if
         
         if (s% log_max_temperature >= s% job% logT_for_conv_vel_flag &
               .and. (.not. s% conv_vel_flag)) then
            write(*,1) 'have reached logT_for_conv_vel_flag', &
               s% log_max_temperature, s% job% logT_for_conv_vel_flag
            write(*,1) 'set conv_vel_flag true'
            call star_set_conv_vel_flag(id, .true., ierr)
            if (failed('star_set_conv_vel_flag',ierr)) return
            if (ierr /= 0) return
         end if
         
         if (s% job% change_RSP2_flag_at_model_number == s% model_number) then
            write(*,*) 'have reached model number for new_RSP2_flag', &
               s% model_number, s% job% new_RSP2_flag
            call star_set_RSP2_flag(id, s% job% new_RSP2_flag, ierr)
            if (failed('star_set_RSP2_flag',ierr)) return
         end if
         
         if (s% job% report_mass_not_fe56) call do_report_mass_not_fe56(s)
         if (s% job% report_cell_for_xm > 0) call do_report_cell_for_xm(s)
         
         model_number = get_model_number(id, ierr)
         if (failed('get_model_number',ierr)) return
         
         if (s% star_age < s% job% set_cumulative_energy_error_each_step_if_age_less_than) then
            if (mod(model_number, s% terminal_interval) == 0) &
               write(*,1) 'cumulative_energy_error reset to', s% job% new_cumulative_energy_error
            s% cumulative_energy_error = s% job% new_cumulative_energy_error
         end if
         
         if (s% doing_timing) then
         
            call system_clock(s% job% time1_extra, s% job% clock_rate)
            s% job% before_step_timing = &
               s% job% before_step_timing + &
                  dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate
            
            s% job% check_time_end = eval_total_times(s% id, ierr)
            s% job% check_before_step_timing = &
               s% job% check_before_step_timing + &
                  (s% job% check_time_end - s% job% check_time_start)

            s% job% time0_extra = s% job% time1_extra
            s% job% check_time_start = s% job% check_time_end

         end if

      end subroutine before_step_loop


      subroutine after_step_loop(id, inlist_fname, &
             dbg, result, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s         
         character (len=*) :: inlist_fname
         logical, intent(in) :: dbg
         integer, intent(out) :: result, ierr
         logical :: will_read_pgstar_inlist

         real(dp) :: tmp
         
         include 'formats'

         call star_ptr(id, s, ierr)
         if (ierr/=0) return

         if (s% doing_timing) then            
            call system_clock(s% job% time1_extra,s% job% clock_rate)
            s% job% step_loop_timing = s% job% step_loop_timing + &
               dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate               
            s% job% check_time_end = eval_total_times(s% id, ierr)
            s% job% check_step_loop_timing = s% job% check_step_loop_timing + &
                (s% job% check_time_end - s% job% check_time_start)
            s% job% time0_extra = s% job% time1_extra
            s% job% check_time_start = s% job% check_time_end               
         end if
         
         if (s% model_number == s% job% set_cumulative_energy_error_at_step) then
            write(*,1) 'set_cumulative_energy_error', s% job% new_cumulative_energy_error
            s% cumulative_energy_error = s% job% new_cumulative_energy_error
         end if
         
         if(s% total_energy_end .ne. 0d0) then
            if (abs(s% cumulative_energy_error/s% total_energy_end) > &
                  s% warn_when_large_rel_run_E_err) then
               write(*,2) 'WARNING: rel_run_E_err', &
                  s% model_number, abs(s% cumulative_energy_error/s% total_energy_end)
            end if
         end if
         
         if (.not. (s% rotation_flag .or. s% u_flag .or. s% use_mass_corrections &
               .or. s% v_flag .or. s% m_center > 0 .or. s% star_mdot /= 0d0)) then
            tmp = abs(1d0 + s% total_gravitational_energy_end/s% virial_thm_P_avg)
            if (tmp > s% warn_when_large_virial_thm_rel_err) then
               write(*,2) 'WARNING: virial_thm_rel_err', &
                  s% model_number, tmp, s% warn_when_large_virial_thm_rel_err, &
                  abs(s% total_gravitational_energy_end), s% virial_thm_P_avg
            end if
         end if
                     
         if (result == keep_going) then 
            if (s% job% pgstar_flag) then
                will_read_pgstar_inlist = .false.
                if (s% pgstar_interval <= 0) then
                    will_read_pgstar_inlist = .true.
                else if(mod(s% model_number, s% pgstar_interval) == 0) then
                    will_read_pgstar_inlist  = .true.
                end if
                if(will_read_pgstar_inlist) then
                  call read_pgstar_inlist(s, inlist_fname, ierr) 
                  if (failed('read_pgstar_controls',ierr)) return
               end if
            end if
         end if
      
         if (result == keep_going) then
            result = s% extras_finish_step(id)      
         else if (result == terminate) then 
            ! call extras_finish_step one last time before terminate
            result = s% extras_finish_step(id)      
            result = terminate
         end if
         
         if (result == keep_going) then
            if (dbg) write(*,*) 'call star_finish_step'
            result = star_finish_step(id, ierr)
            if (failed('star_finish_step',ierr)) return
         end if
                     
         if (result == keep_going .and. s% job% pgstar_flag) then
            if (dbg) write(*,*) 'call update_pgstar_plots'
            call update_pgstar_plots(s, .false., ierr)
            if (failed('update_pgstar_plots',ierr)) return
         end if
         
         if (result == keep_going) then 
            call adjust_tau_factor(s)
            if (s% L_nuc_burn_total/s% L_phot >= s% Lnuc_div_L_zams_limit &
                  .and. .not. s% rotation_flag) then  
               call do_rotation_near_zams(s,ierr)
               if (ierr /= 0) return
            end if           
            if (s% rotation_flag) then      
               call do_rotation(s,ierr)
               if (ierr /= 0) return
            end if 
         end if

      end subroutine after_step_loop


      subroutine terminate_normal_evolve_loop(id, &
             dbg, result, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s         
         logical, intent(in) :: dbg
         integer, intent(out) :: result, ierr
         integer :: i
         include 'formats'

         call star_ptr(id, s, ierr)
         if (ierr/=0) return

         if (dbg) write(*,*) 'call star_pick_next_timestep'
         result = star_pick_next_timestep(id) ! for saved model if any  
         if (dbg) write(*,*) 'call save_profile'
         call save_profile(id, 3, ierr)
         s% need_to_save_profiles_now = .false.
         s% need_to_update_history_now = .true.
         if (dbg) write(*,*) 'call star_finish_step'
         result = star_finish_step(id, ierr)
         if (failed('star_finish_step',ierr)) return         
         if (s% job% save_photo_when_terminate .and. termination_code_string_okay()) &
            s% job% save_photo_number = s% model_number
         if (s% job% save_model_when_terminate .and. termination_code_string_okay()) &
            s% job% save_model_number = s% model_number 
         if (s% job% save_pulse_data_when_terminate) &
            s% job% save_pulse_data_for_model_number = s% model_number
         if (s% job% write_profile_when_terminate) then
            if (len_trim(s% job% filename_for_profile_when_terminate) > 0) then
               call star_write_profile_info( &
                  id, s% job% filename_for_profile_when_terminate, &
                  ierr)
               if (failed('star_write_profile_info',ierr)) return
            else
               write(*,*) "filename_for_profile_when_terminate must be non empty"
               ierr = -1 
               return
            end if
         end if
         if (s% job% show_retry_counts_when_terminate) then
            do i=1,numTlim
               if (s% dt_why_retry_count(i) > 0) then
                  write(*,2) trim(dt_why_str(i)) // ' retries', s% dt_why_retry_count(i)
               end if
            end do
            write(*,*)
         end if
         if (s% job% show_timestep_limit_counts_when_terminate) then
            do i=1,numTlim
               if (s% dt_why_count(i) > 0) then
                  write(*,2) trim(dt_why_str(i)) // ' dt limit', s% dt_why_count(i)
               end if
            end do
            write(*,*)
         end if
         call do_saves(id, ierr)
         if (failed('do_saves terminate_normal_evolve_loop',ierr)) return
         
         contains
         
         logical function termination_code_string_okay()
            integer :: j, n
            termination_code_string_okay = .true.
            if (s% termination_code == 0) return
            n = num_termination_code_strings
            j = maxval(len_trim(s% job% required_termination_code_string(1:n)))
            if (j == 0) return
            termination_code_string_okay = .false.
            do j=1,num_termination_code_strings
               if (s% job% required_termination_code_string(j) == &
                   termination_code_str(s% termination_code)) then
                  termination_code_string_okay = .true.
                  return
               end if
            end do
         end function termination_code_string_okay

      end subroutine terminate_normal_evolve_loop


      subroutine after_evolve_loop(id, &
             do_free_star, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s         
         logical, intent(in) :: do_free_star
         integer, intent(out) :: ierr

         call star_ptr(id, s, ierr)
         if (ierr/=0) return

         if (s% doing_timing) then
            call system_clock(s% job% time1,s% job% clock_rate)
            s% job% elapsed_time = &
                dble(s% job% time1 - s% job% time0_initial) / s% job% clock_rate
            call show_times(id,s)
         end if
         
         if (s% result_reason /= result_reason_normal) then
            write(*, '(a)') 'terminated evolution: ' // &
               trim(result_reason_str(s% result_reason))
         end if
         
         if (s% termination_code > 0 .and. s% termination_code <= num_termination_codes) then
            write(*, '(a)') 'termination code: ' // &
               trim(termination_code_str(s% termination_code))
         end if
         
         if (s% job% pause_before_terminate) then
            write(*,'(a)') 'pause_before_terminate: hit RETURN to continue'
            read(*,*)
         end if

         call s% extras_after_evolve(id, ierr)
         if (failed('after_evolve_extras',ierr)) return

         call star_after_evolve(id, ierr)
         if (failed('star_after_evolve',ierr)) return
         
         if (s% result_reason == result_reason_normal) then
         
            if (s% job% pgstar_flag) &
               call update_pgstar_plots( &
                  s, s% job% save_pgstar_files_when_terminate, &
                  ierr)
            if (failed('update_pgstar_plots',ierr)) return

            call show_terminal_header(id, ierr)
            if (failed('show_terminal_header',ierr)) return
            
            call write_terminal_summary(id, ierr)
            if (failed('write_terminal_summary',ierr)) return
         
         end if
         
         if (len_trim(s% job% echo_at_end) > 0) then
            write(*,*)
            write(*,'(a)') trim(s% job% echo_at_end)
            write(*,*)
         end if
         
         if (do_free_star) then
            call free_star(id, ierr)
            if (failed('free_star',ierr)) return
         end if

      end subroutine after_evolve_loop
         
      
      subroutine adjust_tau_factor(s)
         type (star_info), pointer :: s         
         include 'formats'
         
         if (s% job% adjust_tau_factor_to_surf_density .and. &
               s% job% base_for_adjust_tau_factor_to_surf_density > 0d0) then
            s% tau_factor = s% rho(1)/s% job% base_for_adjust_tau_factor_to_surf_density
            !write(*,1) 'adjust_tau_factor_to_surf_density', s% tau_factor
            s% need_to_setvars = .true.
         end if
         
         if (s% job% set_tau_factor_after_core_He_burn > 0 .and. &
               abs(s% tau_factor - s% job% set_to_this_tau_factor) > &
                  1d-6*max(s% tau_factor, s% job% set_to_this_tau_factor)) then
            if (check_for_after_He_burn(s, s% job% set_tau_factor_after_core_He_burn)) then
               s% tau_factor = s% job% set_to_this_tau_factor
               write(*,1) 'set_tau_factor_after_core_He_burn', s% tau_factor
               s% need_to_setvars = .true.
            end if
         end if
   
         if (s% job% set_tau_factor_after_core_C_burn > 0 .and. &
               abs(s% tau_factor - s% job% set_to_this_tau_factor) > &
                  1d-6*max(s% tau_factor, s% job% set_to_this_tau_factor)) then
            if (check_for_after_C_burn(s, s% job% set_tau_factor_after_core_C_burn)) then
               s% tau_factor = s% job% set_to_this_tau_factor
               write(*,1) 'set_tau_factor_after_core_C_burn', s% tau_factor
               s% need_to_setvars = .true.
            end if
         end if
   
         if (s% job% relax_tau_factor_after_core_He_burn > 0 .and. &
               abs(s% tau_factor - s% job% relax_to_this_tau_factor) > &
                  1d-6*max(s% tau_factor, s% job% relax_to_this_tau_factor)) then
            if (check_for_after_He_burn(s, s% job% relax_tau_factor_after_core_He_burn)) &
               call relax_tau_factor(s)
         end if
   
         if (s% job% relax_tau_factor_after_core_C_burn > 0 .and. &
               abs(s% tau_factor - s% job% relax_to_this_tau_factor) > &
                  1d-6*max(s% tau_factor, s% job% relax_to_this_tau_factor)) then
            if (check_for_after_C_burn(s, s% job% relax_tau_factor_after_core_C_burn)) &
               call relax_tau_factor(s)
         end if


      end subroutine adjust_tau_factor

            
      subroutine do_rotation(s,ierr)
         type (star_info), pointer :: s         
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         
         if (s% model_number <= s% job% set_surf_rotation_v_step_limit) then
            s% job% new_omega = s% job% new_surface_rotation_v*1d5/(s% photosphere_r*Rsun)
            write(*,2) 'surface_rotation_v', s% model_number, s% job% new_surface_rotation_v
            write(*,2) 'omega', s% model_number, s% job% new_omega
            call star_set_uniform_omega(s% id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         
         else if (s% model_number <= s% job% set_omega_step_limit) then
            write(*,2) 'omega', s% model_number, s% job% new_omega
            if (failed('star_surface_omega_crit',ierr)) return
            call star_set_uniform_omega(s% id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         
         else if (s% model_number <= s% job% set_omega_div_omega_crit_step_limit) then
            s% job% new_omega = &
               s% job% new_omega_div_omega_crit*star_surface_omega_crit(s% id, ierr)
            write(*,2) 'omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit
            write(*,2) 'omega', s% model_number, s% job% new_omega
            if (failed('star_surface_omega_crit',ierr)) return
            call star_set_uniform_omega(s% id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if    
      end subroutine do_rotation  
                     
      
      subroutine do_rotation_near_zams(s,ierr)
         type (star_info), pointer :: s         
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
                      
         if (s% job% set_near_zams_surface_rotation_v_steps > 0 .and. &
                  s% job% new_surface_rotation_v /= 0d0) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            s% job% set_surf_rotation_v_step_limit = &
               s% model_number + s% job% set_near_zams_surface_rotation_v_steps - 1
            write(*,2) 'near zams: set_surf_rotation_v_step_limit', &
               s% job% set_surf_rotation_v_step_limit

         else if (s% job% set_near_zams_omega_steps > 0 .and. &
                  s% job% new_omega /= 0d0) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            s% job% set_omega_step_limit = &
               s% model_number + s% job% set_near_zams_omega_steps - 1
            write(*,2) 'near zams: set_omega_step_limit', s% job% set_omega_step_limit

         else if (s% job% set_near_zams_omega_div_omega_crit_steps > 0 .and. &
                  s% job% new_omega_div_omega_crit /= 0d0) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            s% job% set_omega_div_omega_crit_step_limit = &
               s% model_number + s% job% set_near_zams_omega_div_omega_crit_steps - 1
            write(*,2) 'near zams: set_omega_div_omega_crit_step_limit', &
               s% job% set_omega_div_omega_crit_step_limit

         else if (s% job% near_zams_relax_omega .and. &
                  s% job% new_omega /= 0d0) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            write(*,2) 'new_omega', s% model_number, s% job% new_omega
            call star_relax_uniform_omega( &
               s% id, relax_to_new_omega, s% job% new_omega, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return

         else if (s% job% near_zams_relax_omega_div_omega_crit .and. &
                  s% job% new_omega_div_omega_crit /= 0d0) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit
            call star_relax_uniform_omega( &
               s% id, relax_to_new_omega_div_omega_crit, s% job% new_omega_div_omega_crit, &
               s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return

         else if (s% job% near_zams_relax_initial_surface_rotation_v .and. &
                  s% job% new_surface_rotation_v /= 0d0) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            write(*,2) 'new_surface_rotation_v', &
               s% model_number, s% job% new_surface_rotation_v
            call star_relax_uniform_omega( &
               s% id, relax_to_new_surface_rotation_v, s% job% new_surface_rotation_v, &
               s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return

         end if  
      end subroutine do_rotation_near_zams   

         
      subroutine relax_tau_factor(s)
         type (star_info), pointer :: s         
         real(dp) :: next
         include 'formats'
         write(*,*) 'relax_to_this_tau_factor < s% tau_factor', &
            s% job% relax_to_this_tau_factor < s% tau_factor
         write(*,1) 'relax_to_this_tau_factor', s% job% relax_to_this_tau_factor
         write(*,1) 's% tau_factor', s% tau_factor
         if (s% job% relax_to_this_tau_factor < s% tau_factor) then
           next = exp10(safe_log10(s% tau_factor) - s% job% dlogtau_factor)
           if (next < s% job% relax_to_this_tau_factor) &
              next = s% job% relax_to_this_tau_factor
         else
           next = exp10(safe_log10(s% tau_factor) + s% job% dlogtau_factor)
           if (next > s% job% relax_to_this_tau_factor) &
              next = s% job% relax_to_this_tau_factor
         end if
         if (next /= s% tau_factor) then
            s% tau_factor = next
            write(*,1) 'relax_tau_factor', next, s% job% relax_to_this_tau_factor
            s% need_to_setvars = .true.
         end if
      end subroutine relax_tau_factor

         
      subroutine relax_Tsurf_factor(s)
         type (star_info), pointer :: s         
         real(dp) :: next
         include 'formats'
         write(*,*) 'relax_to_this_Tsurf_factor < s% Tsurf_factor', &
            s% job% relax_to_this_Tsurf_factor < s% Tsurf_factor
         write(*,1) 'relax_to_this_Tsurf_factor', s% job% relax_to_this_Tsurf_factor
         write(*,1) 's% Tsurf_factor', s% Tsurf_factor
         if (s% job% relax_to_this_Tsurf_factor < s% Tsurf_factor) then
            next = exp10(safe_log10(s% Tsurf_factor) - s% job% dlogTsurf_factor)
            if (next < s% job% relax_to_this_Tsurf_factor) &
               next = s% job% relax_to_this_Tsurf_factor
         else
            next = exp10(safe_log10(s% Tsurf_factor) + s% job% dlogTsurf_factor)
            if (next > s% job% relax_to_this_Tsurf_factor) &
               next = s% job% relax_to_this_Tsurf_factor
         end if
         s% Tsurf_factor = next
         write(*,1) 'relax_Tsurf_factor', next, s% job% relax_to_this_Tsurf_factor
      end subroutine relax_Tsurf_factor

         
      subroutine check_if_want_to_stop_warnings(s)
         use utils_lib
         type (star_info), pointer :: s    
         character (len=200) :: fname
         integer :: iounit, ierr
         ierr = 0
         if (s% warn_when_large_rel_run_E_err < 1d2) then
            fname = trim(mesa_dir) // '/stop_warnings_for_rel_E_err'
            open(newunit=iounit, file=trim(fname), &
               status='old', action='read', iostat=ierr)
            if (ierr == 0) then
               close(iounit)
               s% warn_when_large_rel_run_E_err = 1d99
               write(*,*) 'turn off warnings for rel_run_E_err'
            end if
         end if
         ierr = 0
      end subroutine check_if_want_to_stop_warnings     
      
      
      logical function stop_is_requested(s)
         type (star_info), pointer :: s         
         integer :: ierr
         logical :: file_exists
         stop_is_requested = .false.
         if (mod(s% model_number,100) /= 0) return
         if (len_trim(s% job% stop_if_this_file_exists) == 0) return
         inquire(file=trim(s% job% stop_if_this_file_exists), exist=file_exists)
         if (.not. file_exists) return
         write(*,*) 'stopping because found file ' // &
            trim(s% job% stop_if_this_file_exists)
         stop_is_requested = .true.
      end function stop_is_requested
      
      
      logical function failed(str,ierr)
         character (len=*), intent(in) :: str
         integer, intent(in) :: ierr
         failed = (ierr /= 0)
         if (failed) write(*, *) trim(str) // ' ierr', ierr
      end function failed
      
      
      subroutine show_times(id, s)
         use utils_lib, only: utils_OMP_GET_MAX_THREADS
         use num_lib, only: qsort
         
         integer, intent(in) :: id
         type (star_info), pointer :: s

         integer, parameter :: max_num_items = 50
         character(len=60) :: item_names(max_num_items)
         real(dp) :: item_values(max_num_items)
         integer, target :: index_arry(max_num_items) 
         integer, pointer :: index(:) 
         integer :: item_order(max_num_items)
         integer :: ierr, omp_num_threads, item_num, num_items, i, j
         real(dp) :: total, misc, tmp
         include 'formats'
         ierr = 0
         omp_num_threads = utils_OMP_GET_MAX_THREADS()
         s% time_total = s% job% check_before_step_timing + &
             s% job% check_step_loop_timing + s% job% check_after_step_timing
         
         write(*,*)
         write(*,'(a50,i18)') 'nz', s% nz
         write(*,'(a50,i18)') 'nvar_total', s% nvar_total
         write(*,'(a50,i18)') trim(s% net_name) // ' species', s% species
         write(*,'(a50,i18)') 'total_num_solver_iterations', &
            s% total_num_solver_iterations
         write(*,'(a50,i18)') 'timing_num_get_eos_calls', &
            s% timing_num_get_eos_calls
         write(*,'(a50,i18)') 'timing_num_solve_eos_calls', &
            s% timing_num_solve_eos_calls
         write(*,'(a50,i18)') 'timing_num_get_kap_calls', &
            s% timing_num_get_kap_calls
         write(*,*)
         write(*,'(a50,i18)') 'threads', omp_num_threads
         total = 0
         item_num = 0
         call save1('remesh', s% time_remesh, total)
         call save1('adjust_mass', s% time_adjust_mass, total)
         call save1('conv_premix', s% time_conv_premix, total)
         call save1('element_diffusion', s% time_element_diffusion, total)
         call save1('burn', s% time_solve_burn, total)
         call save1('mix', s% time_solve_mix, total)
         call save1('solve', s% time_struct_burn_mix, total)
         call save1('matrix', s% time_solver_matrix, total)
         call save1('omega_mix', s% time_solve_omega_mix, total)
         call save1('eos', s% time_eos, total)
         call save1('neu_and_kap', s% time_neu_kap, total)
         call save1('net', s% time_nonburn_net, total)
         call save1('mlt', s% time_mlt, total)
         call save1('hydro_vars', s% time_set_hydro_vars, total)
         call save1('mixing_info', s% time_set_mixing_info, total)
         call save1('evolve_step', s% time_evolve_step, total)
         call save1('run1_star', s% job% elapsed_time - total, total)
         tmp = 0
         call save1('total', total, tmp)
         
         num_items = item_num
         index(1:num_items) => index_arry(1:num_items)
         call qsort(index, num_items, item_values)
         
         write(*,*)
         write(*,*)
         do i=1,num_items
            j = index(num_items+1-i)
            if (item_values(j) == 0d0) cycle
            write(*,'(a50,2f9.3)') trim(item_names(j)), &
               item_values(j), item_values(j)/total
            if (j == num_items) write(*,*)
         end do
         
         if (s% job% step_loop_timing/s% job% elapsed_time < 0.9d0) then
            write(*,*)
            write(*,*)
            write(*,1) 'before_step', s% job% before_step_timing/s% job% elapsed_time
            write(*,1) 'step_loop', s% job% step_loop_timing/s% job% elapsed_time
            write(*,1) 'after_step', s% job% after_step_timing/s% job% elapsed_time
            write(*,*)
         end if
         write(*,*)
         write(*,*)
         
         
         contains
         
         
         subroutine save1(name, value, total)
            use utils_lib, only: is_bad_num
            character (len=*), intent(in) :: name
            real(dp), intent(in) :: value
            real(dp), intent(inout) :: total
            include 'formats'
            item_num = item_num + 1
            item_names(item_num) = name
            item_values(item_num) = value
            total = total + value
         end subroutine save1
         

      end subroutine show_times
      
      
      subroutine do_saves( &
            id, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s

         integer :: ierr
         ierr = 0
      
         call star_ptr(id, s, ierr)
         if (ierr/=0) return        
      
         if (s% model_number == s% job% save_model_number) then
            call star_write_model(id, s% job% save_model_filename, ierr)
            if (failed('star_write_model',ierr)) return
            write(*, *) 'model saved to ' // trim(s% job% save_model_filename)
         end if
      
         if (s% model_number == s% job% save_photo_number) then
            call star_write_photo(id, s% job% save_photo_filename, ierr)
            if (failed('star_write_photo',ierr)) return
            if (len_trim(s% job% save_photo_filename) > 0) &
               write(*, *) 'photo saved to ' // trim(s% job% save_photo_filename)
         end if
         
         if (s% model_number == s% job% save_pulse_data_for_model_number) then
            call star_export_pulse_data(id, s%pulse_data_format, s%job%save_pulse_data_filename, &
                 s%add_center_point_to_pulse_data, s%keep_surface_point_for_pulse_data, &
                 s%add_atmosphere_to_pulse_data, ierr)
            if (failed('star_export_pulse_data',ierr)) return
            write(*, *) 'pulsation data saved to ' // &
               trim(s% job% save_pulse_data_filename)
         end if
         
         if (s% model_number == s% job% profile_model_number) then
            write(*, '(a, i7)') 'save profile for model number', s% model_number
            call save_profile(id, 3, ierr)
            if (failed('save_profile',ierr)) return
         end if
         
      end subroutine do_saves

                  
      subroutine write_colors_info(id, s, ierr)
         use colors_lib
         use colors_def
         use chem_def, only: zsol
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer :: io, i, j
         character (len=strlen) :: fname
         real(dp)  :: log_Teff ! log10 of surface temp
         real(dp)  :: log_L ! log10 of luminosity in solar units
         real(dp)  :: mass ! mass in solar units
         real(dp)  :: Fe_H ! [Fe/H]
         ! output
         real(dp),dimension(bc_total_num_colors) :: results
         real(dp) :: log_g
         
         character(len=strlen),dimension(bc_total_num_colors) :: names
         
         ierr = 0
         
         call get_all_bc_names(names,ierr)
         if (ierr /= 0) then 
            ierr=-1
            call cleanup
            return
         end if
         
         fname = 'colors.log'
         !if (s% doing_first_model_of_run) then
         if (.false.) then
            open(newunit=io, file=trim(fname), action='write', status='replace', iostat=ierr)
            ! write column numbers
            j = 1
            write(io,fmt='(i10)',advance='no') j
            j = j+1
            do i=1,4+bc_total_num_colors
               write(io,fmt='(i25)',advance='no') j
               j = j+1
            end do
            write(io,fmt='(i25)') j
            ! write column labels
            write(io,fmt='(a10)',advance='no') 'model'
            write(io,fmt='(a25)',advance='no') 'log_Teff'
            write(io,fmt='(a25)',advance='no') 'log_L'
            write(io,fmt='(a25)',advance='no') 'mass'
            write(io,fmt='(a25)',advance='no') 'Fe_H'
            do i=1,bc_total_num_colors
               write(io,fmt='(a25)',advance='no') trim(names(i))
            end do
            write(io,fmt='(a25)') 'log_g'
         else
            open(newunit=io, file=trim(fname), action='write', position='append', iostat=ierr)
         end if
         if (ierr /= 0) then
            write(*,*) 'failed to open colors.log'
            call cleanup
            return
         end if
         
         log_Teff = log10(s% Teff)
         log_L = s% log_surface_luminosity
         mass = s% star_mass
         Fe_H = safe_log10(get_current_z_at_point(id, 1, ierr) / zsol)
         log_g = safe_log10(s% grav(1))
         if (ierr /= 0) then
            write(*,*) 'failed in get_current_z_at_point'
            call cleanup
            return
         end if
         
         call get_bcs_all(log_Teff, log_g, Fe_H, results, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in colors_get'
            call cleanup
            return
         end if
         
         1 format(1x,f24.12)
         write(io,fmt='(i10)',advance='no') s% model_number
         write(io,fmt=1,advance='no') log_Teff
         write(io,fmt=1,advance='no') log_L
         write(io,fmt=1,advance='no') mass
         write(io,fmt=1,advance='no') Fe_H
         do i=1,bc_total_num_colors
            write(io,fmt=1,advance='no') results(i)
         end do
         write(io,1) log_g
         
         call cleanup
         
         contains
         
         subroutine cleanup
            close(io)
         end subroutine cleanup
      
      end subroutine write_colors_info
      
      
      subroutine read_masses(filename, masses, nmasses, ierr)
         character (len=*), intent(in) :: filename
         real(dp), pointer, intent(inout) :: masses(:)
         integer, intent(out) :: nmasses, ierr
         call read_items(filename, masses, nmasses, 'masses', ierr)
      end subroutine read_masses
      
      
      subroutine read_items(filename, items, nitems, name, ierr)
         use utils_lib
         use utils_def
         character (len=*), intent(in) :: filename, name
         real(dp), pointer, intent(inout) :: items(:)
         integer, intent(out) :: nitems, ierr
         
         integer :: iounit, n, i, t, capacity
         character (len=strlen) :: buffer, string
         
         nitems = 0
         if (.not. associated(items)) then
            capacity = 10
            allocate(items(capacity))
         else
            capacity = size(items,dim=1)
         end if
         
         ierr = 0

         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open file ' // trim(filename)
            return
         end if
         
         n = 0
         i = 0
         
         do
            t = token(iounit, n, i, buffer, string)
            select case(t)
               case(name_token)
                  if (string == name) then
                     call do_read_items(ierr)
                     if (ierr /= 0) then
                        return
                     end if
                     exit ! for now, nothing else to be read
                  end if
                  call error; return
               case(eof_token)
                  exit
               case default
                  call error; return
            end select
            
         end do
         
         close(iounit)
         
         contains
         
         
         subroutine error
            ierr = -1
            write(*,*) 'error in reading file' // trim(filename)
            close(iounit)
         end subroutine error
         
         
         subroutine do_read_items(ierr)
            integer, intent(out) :: ierr
            real(dp) :: mass
            ierr = 0
            t = token(iounit, n, i, buffer, string)
            if (t /= left_paren_token) then
               call error; return
            end if
         mass_loop: do
               t = token(iounit, n, i, buffer, string)
               if (t /= name_token) then
                  call error; return
               end if
               read(string,fmt=*,iostat=ierr) mass
               if (ierr /= 0) then
                  call error; return
               end if
               nitems = nitems+1
               if (nitems > capacity) then
                  capacity = capacity + 10
                  call realloc_double(items,capacity,ierr)
                  if (ierr /= 0) then
                     call error; return
                  end if
               end if
               items(nitems) = mass
               t = token(iounit, n, i, buffer, string)
               if (t == right_paren_token) exit mass_loop
               if (t /= comma_token) then
                  call error; return
               end if
            end do mass_loop
         end subroutine do_read_items
         
      
      end subroutine read_items
      
      
      subroutine do_report_mass_not_fe56(s)
         use const_def
         type (star_info), pointer :: s
         integer :: k, fe56
         real(dp) :: sumdq
         include 'formats'
         fe56 = s% net_iso(ife56)
         if (fe56 == 0) return
         sumdq = 0
         do k = 1, s% nz
            sumdq = sumdq + s% dq(k)*(1-s% xa(fe56,k))
         end do
         write(*,1) 'R', s% r(1)
         write(*,1) 'g', s% cgrav(1)*s% mstar/(s% r(1)*s% r(1))
         write(*,1) 'mass non fe56', s% xmstar*sumdq, sumdq
         write(*,1) 'M_center (Msun)', s% M_center/Msun
         write(*,1) 'xmstar (g)', s% xmstar
         do k=1,s% nz
            if (fe56 == maxloc(s% xa(:,k),dim=1)) then
               write(*,2) 'mass exterior to fe56 (g)', k, (1d0 - s% q(k))*s% xmstar
               write(*,2) 'mass coord top of fe56 (g)', k, s% q(k)*s% xmstar
               return
            end if
         end do
      end subroutine do_report_mass_not_fe56
      
      
      subroutine do_report_cell_for_xm(s)
         use const_def
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: sumdq, dq
         include 'formats'
         dq = s% job% report_cell_for_xm/s% xmstar
         if (dq > 1) then
            write(*,2) 'report_cell_for_xm > xmstar', s% nz
            return
         end if
         sumdq = 0
         do k = 1, s% nz
            sumdq = sumdq + s% dq(k)
            if (sumdq >= dq) then
               write(*,*)
               write(*,2) 'total mass in cells from 1 to k', k, sumdq*s% xmstar
               write(*,2) 'logT(k)', k, s% lnT(k)/ln10
               write(*,2) 'logRho(k)', k, s% lnd(k)/ln10
               write(*,2) 'entropy(k)', k, exp(s% lnS(k))*amu/kerg
               write(*,2) 'xmstar*q(k)', k, s% xmstar*s% q(k)
               write(*,2) 'q(k)', k, s% q(k)
               write(*,*)
               return
            end if
         end do
         write(*,2) 'total mass in cells from 1 to nz', s% nz, s% xmstar
      end subroutine do_report_cell_for_xm
      
      
      subroutine set_which_rates(id, ierr)
         use rates_def
         use rates_lib
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: which_rate
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% job% set_rates_preference) then
            write(*,*) 'change rates preference to', s% job% new_rates_preference
            s% which_rates(:) = s% job% new_rates_preference
         else
            s% which_rates(:) = rates_NACRE_if_available
         end if
         
         if (len_trim(s% job% set_rate_c12ag) > 0) then
            if (s% job% set_rate_c12ag == 'NACRE') then
               which_rate = use_rate_c12ag_NACRE
            else if (s% job% set_rate_c12ag == 'jina reaclib') then
               which_rate = use_rate_c12ag_JR
            else if (s% job% set_rate_c12ag == 'Kunz') then
               which_rate = use_rate_c12ag_Kunz
            else if (s% job% set_rate_c12ag == 'CF88') then
               which_rate = use_rate_c12ag_CF88
            else if (s% job% set_rate_c12ag == 'Buchmann') then
               write(*,*) 'Buchmann rate for c12ag is not in the current jina reaclib'
               write(*,*) 'to use it, switch to the old jina file '
               write(*,*) 'and use set_rate_c12ag == "jina reaclib"'
               write(*,*) '.'
               ierr = -1
               return
            else
               write(*,*) 'invalid string for set_rate_c12ag ' // trim(s% job% set_rate_c12ag)
               write(*,*) 'options are NACRE, jina reaclib, Kunz, CF88'
               ierr = -1
               return
            end if
            call set_which_rate_c12ag(s% which_rates, which_rate)
         end if
         
         if (len_trim(s% job% set_rate_n14pg) > 0) then
            if (s% job% set_rate_n14pg == 'NACRE') then
               which_rate = use_rate_n14pg_NACRE
            else if (s% job% set_rate_n14pg == 'jina reaclib') then
               which_rate = use_rate_n14pg_JR
            else if (s% job% set_rate_n14pg == 'CF88') then
               which_rate = use_rate_n14pg_CF88
            else if (s% job% set_rate_n14pg == 'Imbriani') then
               write(*,*) 'Imbriani rate for n14pg is not in the current jina reaclib'
               write(*,*) 'to use it, switch to the old jina file '
               write(*,*) 'and use set_rate_n14pg == "jina reaclib"'
               write(*,*) '.'
               ierr = -1
               return
            else
               write(*,*) 'invalid string for set_rate_n14pg ' // trim(s% job% set_rate_n14pg)
               write(*,*) 'options are NACRE, jina reaclib, CF88'
               ierr = -1
               return
            end if
            call set_which_rate_n14pg(s% which_rates, which_rate)
         end if
         
         if (len_trim(s% job% set_rate_3a) > 0) then
            if (s% job% set_rate_3a == 'NACRE') then
               which_rate = use_rate_3a_NACRE
            else if (s% job% set_rate_3a == 'jina reaclib') then
               which_rate = use_rate_3a_JR
            else if (s% job% set_rate_3a == 'CF88') then
               which_rate = use_rate_3a_CF88
            else if (s% job% set_rate_3a == 'FL87') then
               which_rate = use_rate_3a_FL87
            else
               write(*,*) 'invalid string for set_rate_3a ' // trim(s% job% set_rate_3a)
               write(*,*) 'options are NACRE, jina reaclib, CF88, FL87'
               ierr = -1
               return
            end if
            call set_which_rate_3a(s% which_rates, which_rate)
         end if
         
         if (len_trim(s% job% set_rate_1212) > 0) then
            if (s% job% set_rate_1212 == 'CF88_basic_1212') then
               which_rate = use_rate_1212_CF88_basic
            else if (s% job% set_rate_1212 == 'CF88_multi_1212') then
               which_rate = use_rate_1212_CF88_multi
            else
               write(*,*) 'invalid string for set_rate_1212 ' // trim(s% job% set_rate_1212)
               ierr = -1
               return
            end if
            call set_which_rate_1212(s% which_rates, which_rate)
         end if

      end subroutine set_which_rates
      
      
      subroutine set_rate_factors(id, ierr)
         use net_lib, only: get_net_reaction_table_ptr
         use rates_lib, only: rates_reaction_id
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: j, i, ir
         integer, pointer :: net_reaction_ptr(:) 
         logical :: error
         
         include 'formats'
         
         ierr = 0
         error = .false.
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% rate_factors(:) = 1
         if (s% job% num_special_rate_factors <= 0) return
         
         call get_net_reaction_table_ptr(s% net_handle, net_reaction_ptr, ierr)
         if (ierr /= 0) return
         
         do i=1,s% job% num_special_rate_factors
            if (len_trim(s% job% reaction_for_special_factor(i)) == 0) cycle
            ir = rates_reaction_id(s% job% reaction_for_special_factor(i))
            j = 0
            if (ir > 0) j = net_reaction_ptr(ir)
            if (j <= 0) then
               write(*,2) 'Failed to find reaction_for_special_factor ' // &
               trim(s% job% reaction_for_special_factor(i)), &
               j, s% job% special_rate_factor(i)
               error = .true.
               cycle
            end if
            s% rate_factors(j) = s% job% special_rate_factor(i)
            write(*,2) 'set special rate factor for ' // &
                  trim(s% job% reaction_for_special_factor(i)), &
                  j, s% job% special_rate_factor(i)
         end do

         if(error) call mesa_error(__FILE__,__LINE__)
         
      end subroutine set_rate_factors


      subroutine do_star_job_controls_before(id, s, restart, ierr)

         use rates_lib, only: rates_warning_init
         use atm_support, only: get_atm_tau_base

         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         logical, parameter :: kap_use_cache = .true.
         logical :: save_flag
         include 'formats'
      
         ierr = 0

         s% set_which_rates => set_which_rates ! will be called after net is defined
         s% set_rate_factors => set_rate_factors ! will be called after net is defined
         
         call get_atm_tau_base(s, s% tau_base, ierr)
         if (failed('atm_tau_base',ierr)) return

         call rates_warning_init( &
            s% warn_rates_for_high_temp, s% max_safe_logT_for_rates)

      end subroutine do_star_job_controls_before

      
      subroutine do_read_star_job_and_return_id(filename, id, ierr)
         character(*), intent(in) :: filename
         integer, intent(out) :: id  
         integer, intent(out) :: ierr  
         type (star_info), pointer :: s
         character(len=strlen) :: inlist_fname
         
         include 'formats'
         ierr = 0   

         if (id_from_read_star_job /= 0) then
            write(*,2) 'id_from_read_star_job', id_from_read_star_job
            ierr = -1
            return
         end if
         
         call alloc_star(id, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_read_star_job failed in alloc_star'
            return
         end if
         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_read_star_job failed in star_ptr'
            return
         end if
         
         call resolve_inlist_fname(inlist_fname,filename)
         call read_star_job(s, inlist_fname, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr from read_star_job ' // trim(inlist_fname)
            return
         end if
         
         id_from_read_star_job = id

         if (s% job% save_star_job_namelist) then
            call write_star_job(s, s% job% star_job_namelist_name, ierr)
            if (ierr /= 0) then
               write(*,*) 'ierr from write_star_job ' // &
                  trim(s% job% star_job_namelist_name)
               return
            end if
         end if
         
      end subroutine do_read_star_job_and_return_id
      
      ! in a perfect world, we'd pass s as an arg to this routine.
      ! but for backward compatibility for a large number of users
      ! we do it this strange way instead.
      subroutine do_read_star_job(filename, ierr)
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr  
         integer :: id
         call do_read_star_job_and_return_id(filename, id, ierr)
      end subroutine do_read_star_job
      
      
      subroutine do_load1_star(id, s, restart, restart_filename, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         character (len=*), intent(in) :: restart_filename
         integer, intent(out) :: ierr
      
         if (restart) then
            call star_load_restart_photo(id, restart_filename, ierr)
            if (failed('star_load_restart_photo',ierr)) return
         else if (s% job% load_saved_photo) then
            write(*,'(a)') 'load saved photo ' // trim(s% job% saved_photo_name)
            write(*,*)
            call star_load_restart_photo(id, s% job% saved_photo_name, ierr)
            if (failed('star_load_restart_photo',ierr)) return
         else if (s% job% load_saved_model) then
            if (s% job% create_merger_model) then
               write(*,*) 'you have both load_saved_model and create_merger_model set true'
               write(*,*) 'please pick one and try again'
               call mesa_error(__FILE__,__LINE__)
            end if
            if (s% job% create_pre_main_sequence_model) then
               write(*,*) 'you have both load_saved_model and ' // &
                  'create_pre_main_sequence_model set true'
               write(*,*) 'please pick one and try again'
               call mesa_error(__FILE__,__LINE__)
            end if
            if (s% job% create_initial_model) then
               write(*,*) 'you have both load_saved_model and create_initial_model set true'
               write(*,*) 'please pick one and try again'
               call mesa_error(__FILE__,__LINE__)
            end if
            write(*,'(a)') 'load saved model ' // trim(s% job% load_model_filename)
            write(*,*)
            call star_read_model(id, s% job% load_model_filename, ierr)
            if (failed('star_read_model',ierr)) return
         else if (s% job% create_merger_model) then
            call create_merger_model(s, ierr)
            if (failed('create_merger_model',ierr)) return
         else if (s% job% create_pre_main_sequence_model) then
            if (.not. restart) write(*, *) 'create pre-main-sequence model'
            if (s% job% create_initial_model) then
               write(*,*) 'you have both create_pre_main_sequence_model ' // &
                  'and create_initial_model set true'
               write(*,*) 'please pick one and try again'
               call mesa_error(__FILE__,__LINE__)
            end if
            call star_create_pre_ms_model( &
               id, s% job% pre_ms_T_c, s% job% pre_ms_guess_rho_c, &
               s% job% pre_ms_d_log10_P, s% job% pre_ms_logT_surf_limit, &
               s% job% pre_ms_logP_surf_limit, s% job% initial_zfracs, &
               s% job% dump_missing_metals_into_heaviest, &
               (s% job% change_net .or. (s% job% change_initial_net .and. .not. restart)), &
               s% job% new_net_name, s% job% pre_ms_relax_num_steps, ierr)
            if (failed('star_create_pre_ms_model',ierr)) return
         else if (s% job% create_RSP_model) then
            if (.not. restart) write(*, *) 'create initial RSP model'
            call star_create_RSP_model(id, ierr)
            if (failed('star_create_RSP_model',ierr)) return
         else if (s% job% create_RSP2_model) then
            if (.not. restart) write(*, *) 'create initial RSP2 model'
            call star_create_RSP2_model(id, ierr)
            if (failed('star_create_RSP_model',ierr)) return
         else if (s% job% create_initial_model) then
            if (.not. restart) write(*, *) 'create initial model'
            if (s% job% create_pre_main_sequence_model) then
               write(*,*) 'you have both create_initial_model and ' // &
                  'create_pre_main_sequence_model set true'
               write(*,*) 'please pick one and try again'
               call mesa_error(__FILE__,__LINE__)
            end if
            call star_create_initial_model(id, &
               s% job% radius_in_cm_for_create_initial_model, &
               s% job% mass_in_gm_for_create_initial_model, &
               s% job% center_logP_1st_try_for_create_initial_model, &
               s% job% entropy_1st_try_for_create_initial_model, &
               s% job% max_tries_for_create_initial_model, &
               s% job% abs_e01_tolerance_for_create_initial_model, &
               s% job% abs_e02_tolerance_for_create_initial_model, &
               s% job% initial_zfracs, &
               s% job% dump_missing_metals_into_heaviest, &
               (s% job% change_net .or. (s% job% change_initial_net .and. .not. restart)), &
               s% job% new_net_name, s% job% initial_model_relax_num_steps, &
               s% job% initial_model_eps, &
               ierr)
            if (failed('star_create_initial_model',ierr)) return
         else
            call star_load_zams(id, ierr)
            if (failed('star_load_zams',ierr)) return
         end if

      end subroutine do_load1_star


      subroutine create_merger_model(s, ierr)
         use ctrls_io, only : store_controls
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: id, id_aux, i, j, k
         type (star_info), pointer :: s_aux
         real(dp), pointer :: xq(:), xa(:,:)
         real(dp) :: total_mass, partial_mass

         include 'formats'
         ierr = 0
         id = s% id

         if (s% job% create_pre_main_sequence_model) then
            write(*,*) 'you have both load_saved_model and ' // &
               'create_pre_main_sequence_model set true'
            write(*,*) 'please pick one and try again'
            call mesa_error(__FILE__,__LINE__)
         end if
         if (s% job% create_initial_model) then
            write(*,*) 'you have both load_saved_model and create_initial_model set true'
            write(*,*) 'please pick one and try again'
            call mesa_error(__FILE__,__LINE__)
         end if
         !load first star
         call star_read_model(id, s% job% saved_model_for_merger_1, ierr)
         if (failed('star_read_model',ierr)) return

         !load second star
         call alloc_star(id_aux, ierr)
         if (failed('alloc_star',ierr)) return
         call star_ptr(id_aux, s_aux, ierr)
         if (failed('star_ptr',ierr)) return
         call init_starting_star_data(s_aux, ierr)
         if (failed('init_starting_star_data',ierr)) return
         call star_set_kap_and_eos_handles(id_aux, ierr)
         if (failed('set_star_kap_and_eos_handles',ierr)) return
         call store_controls(s_aux, ierr)
         if (failed('store_controls',ierr)) return
         call do_star_job_controls_before(id_aux, s_aux, .false., ierr)
         if (ierr /= 0) return
         s_aux% job% set_rate_c12ag = s% job% set_rate_c12ag
         s_aux% job% set_rate_n14pg = s% job% set_rate_n14pg
         s_aux% job% set_rate_3a = s% job% set_rate_3a
         s_aux% job% set_rate_1212 = s% job% set_rate_1212
         call star_read_model(id_aux, s% job% saved_model_for_merger_2, ierr)
         if (failed('star_read_model',ierr)) return

         ! create composition and q array through an entropy sorting
         total_mass = s% mstar + s_aux% mstar
         partial_mass = 0
         i = 1
         j = 1
         allocate(xq(s% nz + s_aux% nz), xa(s% species, s% nz + s_aux% nz))
         do while (i <= s% nz .or. j <= s_aux% nz)
            if (j > s_aux% nz .or. (i <= s% nz .and. &
               s% entropy(i) >= s_aux% entropy(j))) then
                  partial_mass = partial_mass + s% dm(i)
                  do k=1, s% species
                     xa(k, i+j-1) = s% xa(k, i)
                  end do
                  i = i + 1
            else if (i > s% nz .or. (j <= s_aux% nz .and. &
               s_aux% entropy(j) > s% entropy(i))) then
                  partial_mass = partial_mass + s_aux% dm(j)
                  do k=1, s% species
                     xa(k, i+j-1) = s_aux% xa(k, j)
                  end do
                  j = j + 1
            end if
            xq(i+j-2) = partial_mass / total_mass
            !write(*,*) "check", i+j-2, xq(i+j-2), xa(1, i+j-2), xa(2, i+j-2), xa(3, i+j-2)
         end do
         ! Relax composition first, then composition mass
         ! Turn off rotation for relaxation
         call star_set_rotation_flag(id, .false., ierr)
         if (failed('star_set_rotation_flag',ierr)) then
            deallocate(xq,xa)
            return
         end if
         write(*,*) "Relaxing composition to merger composition"
         call star_relax_composition( &
            id, s% job% num_steps_to_relax_composition, s% nz + s_aux% nz, s% species, xa, xq, ierr)
         if (failed('star_relax_composition',ierr)) then
            deallocate(xq,xa)
            return
         end if
         write(*,*) "Relaxing star mass to total merger mass"
         call star_relax_mass_scale( &
            id, total_mass/Msun, s% job% dlgm_per_step, &
            s% job% change_mass_years_for_dt, ierr)
         deallocate(xq,xa)
         if (failed('star_relax_mass_scale',ierr)) return
         
      end subroutine create_merger_model
      

      subroutine extend_net(s, ierr)
         use net_def
         use chem_def
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp), parameter :: tiny = 1d-10, small = 1d-2
         
         real(dp) :: cntr_h, cntr_he
         
         include 'formats'
         
         ierr = 0
         
         !write(*,2) 'extend_net: current net ' // trim(s% net_name), s% model_number
         
         if (s% net_name == s% job% adv_net) return

         if (s% net_name == s% job% co_net) then
            if (s% log_max_temperature > 9d0 .or. s% log_center_density > 9d0) then
               call change_net(s% job% adv_net)
               if (len_trim(s% job% profile_columns_file) > 0) &
                  write(*,*) 'read ' // trim(s% job% profile_columns_file)
               call star_set_profile_columns( &
                  s% id, s% job% profile_columns_file, .true., ierr)
            end if
            return
         end if
         
         if (s% net_name == s% job% h_he_net) then
            cntr_h = current_abundance_at_point(s% id, ih1, s% nz, ierr)
            !write(*,2) 'cntr_h', s% model_number, cntr_h, tiny
            if (ierr /= 0) return
            if (cntr_h > tiny) return
            cntr_he = current_abundance_at_point(s% id, ihe4, s% nz, ierr)
            !write(*,2) 'cntr_he', s% model_number, cntr_he, small
            if (ierr /= 0) return
            if (cntr_he > small) return
            if (s% log_max_temperature > 8.3d0 .or. s% log_center_density > 8.5d0) then
               call change_net(s% job% co_net)
               if (len_trim(s% job% profile_columns_file) > 0) &
                  write(*,*) 'read ' // trim(s% job% profile_columns_file)
               call star_set_profile_columns( &
                  s% id, s% job% profile_columns_file, .true., ierr)
            end if
         end if
         
         
         contains
                  
         
         subroutine change_net(net_name)
            use const_def
            character (len=*), intent(in) :: net_name
            integer :: j
            
            include 'formats'
            
            call star_change_to_new_net( &
               s% id, s% job% adjust_abundances_for_new_isos, net_name, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_change_to_new_net ' // trim(net_name)
               stop 'change_net'
               return
            end if
            
            if (net_name /= s% net_name) then
               write(*,*) '   new net_name ', trim(net_name)
               write(*,*) 'old s% net_name ', trim(s% net_name)
               write(*,*) 'failed to change'
               stop 'change_net'
            end if

            write(*,'(a)') ' new net = ' // trim(s% net_name)
            !do j=1,s% species
            !   write(*,fmt='(a,x)',advance='no') trim(chem_isos% name(s% chem_id(j)))
            !end do
            !write(*,*)
            s% dt_next = s% dt_next/5
            !write(*,1) 'reduce timestep', log10(s% dt_next/secyer)
            write(*,*)
         end subroutine change_net
         
         
      end subroutine extend_net         


      subroutine before_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine before_evolve         
       

      subroutine do_star_job_controls_after(id, s, restart, pgstar_ok, ierr)
         use const_def
         use rates_def
         use rates_lib
         use utils_lib, only: utils_OMP_GET_MAX_THREADS

         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart, pgstar_ok
         integer, intent(out) :: ierr
         
         real(dp) :: log_m, log_lifetime, max_dt, max_timestep, minq, maxq
         integer :: i, j, k, nzlo, nzhi, chem_id, chem_id1, chem_id2
         logical :: change_v, change_u
         include 'formats'
         
         if (len_trim(s% job% history_columns_file) > 0) &
            write(*,*) 'read ' // trim(s% job% history_columns_file)
         call star_set_history_columns(id, s% job% history_columns_file, .true., ierr)
         if (failed('star_set_history_columns',ierr)) return
         
         if (len_trim(s% job% profile_columns_file) > 0) &
            write(*,*) 'read ' // trim(s% job% profile_columns_file)
         call star_set_profile_columns(id, s% job% profile_columns_file, .true., ierr)
         if (failed('star_set_profile_columns',ierr)) return
         
         if (pgstar_ok) then
            if (s% job% clear_pgstar_history .or. &
                  (s% job% clear_initial_pgstar_history .and. .not. restart)) then
               call start_new_run_for_pgstar(s, ierr)
               if (failed('start_new_run_for_pgstar',ierr)) return
            else
               call restart_run_for_pgstar(s, ierr)
               if (failed('restart_run_for_pgstar',ierr)) return
            end if
         end if
         
         if (s% job% set_tau_factor .or. &
               (s% job% set_initial_tau_factor .and. .not. restart)) then
            write(*,1) 'set_tau_factor', s% job% set_to_this_tau_factor
            s% tau_factor = s% job% set_to_this_tau_factor
         end if
         
         if (s% job% set_Tsurf_factor .or. &
               (s% job% set_initial_Tsurf_factor .and. .not. restart)) then
            write(*,1) 'set_Tsurf_factor', s% job% set_to_this_Tsurf_factor
            s% Tsurf_factor = s% job% set_to_this_Tsurf_factor
         end if
         
         if (s% job% set_initial_age .and. .not. restart) then
            write(*,1) 'set_initial_age', s% job% initial_age ! in years
            call star_set_age(id, s% job% initial_age, ierr)
            if (failed('star_set_age',ierr)) return
         end if

         if (s% job% set_initial_dt .and. .not. restart) then
            if (s% job% years_for_initial_dt > 0d0) then
               write(*,1) 'set_initial_dt (years)', s% job% years_for_initial_dt
               s% dt_next = s% job% years_for_initial_dt*secyer
            else if (s% job% seconds_for_initial_dt > 0d0) then
               write(*,1) 'set_initial_dt (seconds)', s% job% seconds_for_initial_dt
               s% dt_next = s% job% seconds_for_initial_dt
            end if
         end if

         if (s% job% limit_initial_dt .and. .not. restart) then
            if (s% job% years_for_initial_dt > 0d0) then
               write(*,1) 'limit_initial_dt (years)', s% job% years_for_initial_dt
               s% dt_next = min(s% dt_next, s% job% years_for_initial_dt*secyer)
            else if (s% job% seconds_for_initial_dt > 0d0) then
               write(*,1) 'limit_initial_dt (seconds)', s% job% seconds_for_initial_dt
               s% dt_next = min(s% dt_next, s% job% seconds_for_initial_dt)
            end if

         end if

         ! enforce max_timestep on first step

         if (s% max_years_for_timestep > 0) then
            max_timestep = secyer*s% max_years_for_timestep
            if (s% max_timestep > 0 .and. s% max_timestep < max_timestep) &
                 max_timestep = s% max_timestep
         else
            max_timestep = s% max_timestep
         end if

         if (max_timestep > 0 .and. max_timestep < s% dt_next) then
            write(*,1) 'max_timestep (seconds)', max_timestep
            s% dt_next = max_timestep
         endif

         if (s% job% set_initial_model_number .and. .not. restart) then
            write(*,2) 'set_initial_model_number', s% job% initial_model_number
            s% model_number = s% job% initial_model_number
            s% init_model_number = s% model_number
         end if

         if (s% job% set_initial_number_retries .and. .not. restart) then
            write(*,2) 'set_initial_number_retries', s% job% initial_number_retries
            s% num_retries = s% job% initial_number_retries
         end if

         if (s% job% steps_to_take_before_terminate >= 0) then
            s% max_model_number = s% model_number + s% job% steps_to_take_before_terminate
            write(*,2) 'steps_to_take_before_terminate', &
               s% job% steps_to_take_before_terminate
            write(*,2) 'max_model_number', s% max_model_number
         end if

         if (s% job% steps_before_start_timing > 0) then
            s% job% first_model_for_timing = s% model_number + s% job% steps_before_start_timing
            write(*,2) 'steps_before_start_timing', &
               s% job% steps_before_start_timing
         end if
         
         if (s% job% change_net .or. (s% job% change_initial_net .and. .not. restart)) then         
            call star_change_to_new_net( &
               id, s% job% adjust_abundances_for_new_isos, s% job% new_net_name, ierr)
            if (failed('star_change_to_new_net',ierr)) return
         end if

         if (s% job% change_small_net .or. &
               (s% job% change_initial_small_net .and. .not. restart)) then         
            write(*,*) 'change small net to ' // trim(s% job% new_small_net_name)
            call star_change_to_new_small_net( &
               id, s% job% adjust_abundances_for_new_isos, s% job% new_small_net_name, ierr)
            if (failed('star_change_to_new_small_net',ierr)) return
            write(*,*) 'number of species', s% species
         end if
         
         if (abs(s% job% T9_weaklib_full_off - T9_weaklib_full_off) > 1d-6) then
            write(*,1) 'set T9_weaklib_full_off', s% job% T9_weaklib_full_off
            T9_weaklib_full_off = s% job% T9_weaklib_full_off
         end if
         
         if (abs(s% job% T9_weaklib_full_on - T9_weaklib_full_on) > 1d-6) then
            write(*,1) 'set T9_weaklib_full_on', s% job% T9_weaklib_full_on
            T9_weaklib_full_on = s% job% T9_weaklib_full_on
         end if
         
         if (s% job% weaklib_blend_hi_Z /= weaklib_blend_hi_Z) then
            write(*,1) 'set weaklib_blend_hi_Z', s% job% weaklib_blend_hi_Z
            weaklib_blend_hi_Z = s% job% weaklib_blend_hi_Z
         end if
         
         if (abs(s% job% T9_weaklib_full_off_hi_Z - T9_weaklib_full_off_hi_Z) > 1d-6) then
            write(*,1) 'set T9_weaklib_full_off_hi_Z', s% job% T9_weaklib_full_off_hi_Z
            T9_weaklib_full_off_hi_Z = s% job% T9_weaklib_full_off_hi_Z
         end if
         
         if (abs(s% job% T9_weaklib_full_on_hi_Z - T9_weaklib_full_on_hi_Z) > 1d-6) then
            write(*,1) 'set T9_weaklib_full_on_hi_Z', s% job% T9_weaklib_full_on_hi_Z
            T9_weaklib_full_on_hi_Z = s% job% T9_weaklib_full_on_hi_Z
         end if

         ! set up coulomb corrections for the special weak rates
         which_mui_coulomb = get_mui_value(s% job% ion_coulomb_corrections)
         which_vs_coulomb = get_vs_value(s% job% electron_coulomb_corrections)
         
         change_v = s% job% change_v_flag .or. &
               (s% job% change_initial_v_flag .and. .not. restart)
         change_u = s% job% change_u_flag .or. &
               (s% job% change_initial_u_flag .and. .not. restart)
         if (change_v .or. change_u) then
            ! do add new before remove old so can set initial values
            if (change_v .and. s% job% new_v_flag) then
               write(*,*) 'new_v_flag', s% job% new_v_flag
               call star_set_v_flag(id, s% job% new_v_flag, ierr)
               if (failed('star_set_v_flag',ierr)) return
            end if
            if (change_u .and. s% job% new_u_flag) then
               write(*,*) 'new_u_flag', s% job% new_u_flag
               call star_set_u_flag(id, s% job% new_u_flag, ierr)
               if (failed('star_set_u_flag',ierr)) return
            end if
            if (change_v .and. .not. s% job% new_v_flag) then
               write(*,*) 'new_v_flag', s% job% new_v_flag
               call star_set_v_flag(id, s% job% new_v_flag, ierr)
               if (failed('star_set_v_flag',ierr)) return
            end if
            if (change_u .and. .not. s% job% new_u_flag) then
               write(*,*) 'new_u_flag', s% job% new_u_flag
               call star_set_u_flag(id, s% job% new_u_flag, ierr)
               if (failed('star_set_u_flag',ierr)) return
            end if
         end if

         if (s% job% change_RTI_flag .or. &
               (s% job% change_initial_RTI_flag .and. .not. restart)) then
            write(*,*) 'new_RTI_flag', s% job% new_RTI_flag
            call star_set_RTI_flag(id, s% job% new_RTI_flag, ierr)
            if (failed('star_set_RTI_flag',ierr)) return
         end if
         
         if (s% job% change_RSP2_flag .or. &
               (s% job% change_initial_RSP2_flag .and. .not. restart)) then
            write(*,*) 'new_RSP2_flag', s% job% new_RSP2_flag
            call star_set_RSP2_flag(id, s% job% new_RSP2_flag, ierr)
            if (failed('star_set_RSP2_flag',ierr)) return
         end if

         if (s% job% change_RSP_flag .or. &
               (s% job% change_initial_RSP_flag .and. .not. restart)) then
            write(*,*) 'new_RSP_flag', s% job% new_RSP_flag
            call star_set_RSP_flag(id, s% job% new_RSP_flag, ierr)
            if (failed('star_set_RSP_flag',ierr)) return
         end if

         if (s% job% change_conv_vel_flag .or. &
               (s% job% change_initial_conv_vel_flag .and. .not. restart)) then
            write(*,*) 'new_conv_vel_flag', s% job% new_conv_vel_flag
            call star_set_conv_vel_flag(id, s% job% new_conv_vel_flag, ierr)
            if (failed('star_set_conv_vel_flag',ierr)) return
         end if

         if (s% job% change_w_div_wc_flag .or. &
               (s% job% change_initial_w_div_wc_flag .and. .not. restart)) then
            write(*,*) 'new_w_div_wc_flag', s% job% new_w_div_wc_flag
            call star_set_w_div_wc_flag(id, s% job% new_w_div_wc_flag, ierr)
            if (failed('star_set_w_div_wc_flag',ierr)) return
         end if

         if (s% job% change_j_rot_flag .or. &
               (s% job% change_initial_j_rot_flag .and. .not. restart)) then
            write(*,*) 'new_j_rot_flag', s% job% new_j_rot_flag
            call star_set_j_rot_flag(id, s% job% new_j_rot_flag, ierr)
            if (failed('star_set_j_rot_flag',ierr)) return
         end if

         if (s% job% change_D_omega_flag .or. &
               (s% job% change_initial_D_omega_flag .and. .not. restart)) then
            call star_set_D_omega_flag(id, s% job% new_D_omega_flag, ierr)
            if (failed('star_set_D_omega_flag',ierr)) return
         end if

         if (s% job% change_am_nu_rot_flag .or. &
               (s% job% change_initial_am_nu_rot_flag .and. .not. restart)) then
            call star_set_am_nu_rot_flag(id, s% job% new_am_nu_rot_flag, ierr)
            if (failed('star_set_am_nu_rot_flag',ierr)) return
         end if

         if (s% job% change_rotation_flag .or. &
               (s% job% change_initial_rotation_flag .and. .not. restart)) then
            write(*,*) 'new_rotation_flag', s% job% new_rotation_flag
            call star_set_rotation_flag(id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% set_omega) then
            write(*,1) 'new_omega', s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% set_initial_omega .and. .not. restart) then
            write(*,1) 'new_omega', s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% set_surface_rotation_v) then
            s% job% new_omega = s% job% new_surface_rotation_v*1d5/s% r(1)
            write(*,1) 'new_surface_rotation_v', &
               s% job% new_surface_rotation_v, s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. &
             s% job% set_initial_surface_rotation_v .and. .not. restart) then
            s% job% new_omega = s% job% new_surface_rotation_v*1d5/s% r(1)
            write(*,2) 'new_surface_rotation_v', &
               s% model_number, s% job% new_surface_rotation_v, s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% set_omega_div_omega_crit) then
            s% job% new_omega = &
               s% job% new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
            if (failed('star_surface_omega_crit',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit, s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. &
             s% job% set_initial_omega_div_omega_crit .and. .not. restart) then
            s% job% new_omega = &
               s% job% new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
            if (failed('star_surface_omega_crit',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit, s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if
         
         if (s% job% set_to_xa_for_accretion .or. &
               (s% job% set_initial_to_xa_for_accretion .and. .not. restart)) then
            write(*,*) 'set_to_xa_for_accretion'
            call change_to_xa_for_accretion(id, s% job% set_nzlo, s% job% set_nzhi, ierr)
            if (failed('set_to_xa_for_accretion',ierr)) return
         end if
         
         if (s% job% first_model_for_timing > 0) &
            write(*,2) 'first_model_for_timing', s% job% first_model_for_timing
         
         if (s% job% set_uniform_initial_composition .and. .not. restart) then
            write(*,*)
            write(*,1) 'set_uniform_initial_composition'
            write(*,1) 'initial_h1', s% job% initial_h1
            write(*,1) 'initial_h2', s% job% initial_h2
            write(*,1) 'initial_he3', s% job% initial_he3
            write(*,1) 'initial_he4', s% job% initial_he4
            select case(s% job% initial_zfracs)
               case (AG89_zfracs)
                  write(*,1) 'metals AG89'
               case (GN93_zfracs)
                  write(*,1) 'metals GN93'
               case (GS98_zfracs)
                  write(*,1) 'metals GS98'
               case (L03_zfracs)
                  write(*,1) 'metals L03'
               case (AGS05_zfracs)
                  write(*,1) 'metals AGS05'
               case (AGSS09_zfracs)
                  write(*,1) 'metals AGSS09'
               case (L09_zfracs)
                  write(*,1) 'metals L09'
               case (A09_Prz_zfracs)
                  write(*,1) 'metals A09_Prz'
               case default
                  write(*,2) 'unknown value for initial_zfracs', s% job% initial_zfracs
            end select
            call star_set_standard_composition( &
               id, s% job% initial_h1, s% job% initial_h2, &
               s% job% initial_he3, s% job% initial_he4, s% job% initial_zfracs, &
               s% job% dump_missing_metals_into_heaviest, ierr)
            if (failed('set_uniform_initial_composition',ierr)) return
         end if
         
         if (s% job% relax_initial_composition .and. .not. restart) then
            call do_relax_initial_composition(ierr)
            if (failed('do_relax_initial_composition',ierr)) return
         end if
         
         if (s% job% relax_initial_to_xaccrete .and. .not. restart) then
            call star_relax_to_xaccrete(id, s% job% num_steps_to_relax_composition, ierr)
            if (failed('star_relax_to_xaccrete',ierr)) return
         end if

         if (s% job% set_uniform_xa_from_file) then
            call star_uniform_xa_from_file(id, s% job% file_for_uniform_xa, ierr)
            if (failed('star_uniform_xa_from_file',ierr)) return
         end if
         
         if (s% job% relax_initial_angular_momentum .and. .not. restart) then
            call do_relax_initial_angular_momentum(ierr)
            if (failed('do_relax_initial_angular_momentum',ierr)) return
         end if
         
         if (s% job% relax_initial_entropy .and. .not. restart) then
            call do_relax_initial_entropy(ierr)
            if (failed('do_relax_initial_entropy',ierr)) return
         end if
         
         if (s% job% mix_section .or. &
               (s% job% mix_initial_section .and. .not. restart)) then
            write(*,*) 'mix_section'
            call uniform_mix_section( &
               id, s% job% mix_section_nzlo, s% job% mix_section_nzhi, ierr)
            if (failed('uniform_mix_section',ierr)) return
         end if

         if (s% job% mix_initial_envelope_down_to_T > 0d0 .and. .not. restart) then
            call uniform_mix_envelope_down_to_T(id, s% job% mix_initial_envelope_down_to_T, ierr)
            if (failed('uniform_mix_envelope_down_to_T',ierr)) return
         end if

         if (s% job% mix_envelope_down_to_T > 0d0) then
            call uniform_mix_envelope_down_to_T(id, s% job% mix_envelope_down_to_T, ierr)
            if (failed('uniform_mix_envelope_down_to_T',ierr)) return
         end if

         if (s% job% mix_initial_envelope_down_to_T > 0d0) then
            call uniform_mix_envelope_down_to_T(id, s% job% mix_initial_envelope_down_to_T, ierr)
            if (failed('uniform_mix_envelope_down_to_T',ierr)) return
         end if

         if (s% job% set_uniform_initial_xa_from_file .and. .not. restart) then
            call star_uniform_xa_from_file(id, s% job% file_for_uniform_xa, ierr)
            if (failed('star_uniform_xa_from_file',ierr)) return
         end if
         
         ! do change Z before change Y since changing Z can change Y
         if (s% job% change_Z) then
            call star_set_z(id, s% job% new_Z, ierr)
            if (failed('star_set_z',ierr)) return
         end if

         if (s% job% change_initial_Z .and. .not. restart) then
            call star_set_z(id, s% job% new_Z, ierr)
            if (failed('star_set_z',ierr)) return
         end if

         if (s% job% change_Y) then
            call star_set_y(id, s% job% new_Y, ierr)
            if (failed('change_Y',ierr)) return
         end if

         if (s% job% change_initial_Y .and. .not. restart) then
            call star_set_y(id, s% job% new_Y, ierr)
            if (failed('change_initial_Y',ierr)) return
         end if

         if (s% job% zero_alpha_RTI .or. &
               (s% job% zero_initial_alpha_RTI .and. .not. restart)) then
            call star_zero_alpha_RTI(id, ierr)
            if (failed('star_zero_alpha_RTI',ierr)) return
         end if

         if (s% job% set_abundance .or. &
               (s% job% set_initial_abundance .and. .not. restart)) then
            nzlo = s% job% set_abundance_nzlo
            nzhi = s% job% set_abundance_nzhi
            if (nzhi <= 0) nzhi = s% nz
            if (nzlo <= 0) nzlo = 1
            write(*, *) 'set_abundance of ', &
               trim(s% job% chem_name), s% job% new_frac, nzlo, nzhi 
            chem_id = get_nuclide_index(s% job% chem_name)
            if (chem_id <= 0) then
               write(*,*) 'failed to find ' // trim(s% job% chem_name)
               write(*,*) 'check valid chem_isos% names in chem/public/chem_def.f'
            end if
            call set_abundance_in_section(id, chem_id, s% job% new_frac, nzlo, nzhi, ierr)
            if (failed('set_abundance_in_section',ierr)) return
         end if
         
         if (s% job% replace_element .or. &
               (s% job% replace_initial_element .and. .not. restart)) then
            write(*, *) 'replace_element ', &
               trim(s% job% chem_name1), ' by ', trim(s% job% chem_name2)
            chem_id1 = get_nuclide_index(s% job% chem_name1)
            chem_id2 = get_nuclide_index(s% job% chem_name2)
            if (chem_id1 <= 0) then
               write(*,*) 'failed to find ' // trim(s% job% chem_name1)
               write(*,*) 'check valid chem_isos% names in chem/public/chem_def.f'
            end if
            if (chem_id2 <= 0) then
               write(*,*) 'failed to find ' // trim(s% job% chem_name2)
               write(*,*) 'check valid chem_isos% names in chem/public/chem_def.f'
            end if
            nzhi = s% job% replace_element_nzhi
            nzlo = s% job% replace_element_nzlo
            if (nzhi <= 0) nzhi = s% nz
            if (nzlo <= 0) nzlo = 1
            write(*, *) 'in section', nzlo, nzhi
            call replace_element_in_section( &
               id, chem_id1, chem_id2, nzlo, nzhi, ierr)
            if (failed('replace_element_in_section',ierr)) return
         end if

         if (s% job% set_irradiation .or. &
               (s% job% set_initial_irradiation .and. .not. restart)) then
            write(*,2) 'set_irradiation'
            s% irradiation_flux = s% job% set_to_this_irrad_flux
            s% column_depth_for_irradiation = s% job% irrad_col_depth
         end if
         
         if (s% job% do_special_test) then
            write(*, *) 'do_special_test'
            call star_special_test(id, ierr)
            if (failed('star_special_test',ierr)) return
         end if
         
         if (s% job% set_v_center .or. &
               (s% job% set_initial_v_center .and. .not. restart)) then
            write(*, 1) 'set_v_center', s% job% new_v_center
            s% v_center = s% job% new_v_center
         end if
         
         if (s% job% set_L_center .or. &
               (s% job% set_initial_L_center .and. .not. restart)) then
            write(*, 1) 'set_L_center', s% job% new_L_center
            s% L_center = s% job% new_L_center*Lsun
         end if

         ! do "set" before "relax"
         
         ! must do relax Z before relax Y since relax Z can change Y
         ! (Warrick Ball pointed out this requirement)
         if (s% job% relax_initial_Z .and. .not. restart) then
            write(*,1) 'relax_initial_Z', s% job% new_Z
            call star_relax_Z(id, s% job% new_Z, s% relax_dlnZ, &
               s% job% relax_Z_minq, s% job% relax_Z_maxq, ierr)
            if (failed('star_relax_Z',ierr)) return
            write(*, 1) 'new z', get_current_z(id, ierr)
            if (failed('get_current_z',ierr)) return
         end if

         if (s% job% relax_Z) then
            write(*,1) 'relax_Z', s% job% new_Z
            call star_relax_Z(id, s% job% new_Z, s% relax_dlnZ, &
               s% job% relax_Z_minq, s% job% relax_Z_maxq, ierr)
            if (failed('star_relax_Z',ierr)) return
            write(*, 1) 'new z', get_current_z(id, ierr)
            if (failed('get_current_z',ierr)) return
         end if

         if (s% job% relax_initial_Y .and. .not. restart) then
            write(*,1) 'relax_initial_Y', s% job% new_Y
            call star_relax_Y(id, s% job% new_Y, s% relax_dY, &
               s% job% relax_Y_minq, s% job% relax_Y_maxq, ierr)
            if (failed('star_relax_Y',ierr)) return
            write(*, 1) 'new y', get_current_y(id, ierr)
            if (failed('get_current_y',ierr)) return
         end if

         if (s% job% relax_Y) then
            write(*,1) 'relax_Y', s% job% new_Y
            call star_relax_Y(id, s% job% new_Y, s% relax_dY, &
               s% job% relax_Y_minq, s% job% relax_Y_maxq, ierr)
            if (failed('star_relax_Y',ierr)) return
            write(*, 1) 'new y', get_current_y(id, ierr)
            if (failed('get_current_y',ierr)) return
         end if

         if (s% job% relax_mass) then
            write(*, 1) 'relax_mass', s% job% new_mass
            call star_relax_mass(id, s% job% new_mass, s% job% lg_max_abs_mdot, ierr)
            if (failed('star_relax_mass',ierr)) return
         end if

         if (s% job% relax_mass_to_remove_H_env) then
            write(*, 1) 'relax_mass_to_remove_H_env_mass'
            call star_relax_mass_to_remove_H_env( &
               id, s% job% extra_mass_retained_by_remove_H_env, s% job% lg_max_abs_mdot, ierr)
            if (failed('star_relax_mass_to_remove_H_env',ierr)) return
         end if

         if (s% job% relax_dxdt_nuc_factor .or. &
               (s% job% relax_initial_dxdt_nuc_factor .and. .not. restart)) then
            write(*, 1) 'relax_dxdt_nuc_factor', s% job% new_dxdt_nuc_factor
            call star_relax_dxdt_nuc_factor( &
               id, s% job% new_dxdt_nuc_factor, s% job% dxdt_nuc_factor_multiplier, ierr)
            if (failed('star_relax_dxdt_nuc_factor',ierr)) return
         end if

         if (s% job% relax_eps_nuc_factor .or. &
               (s% job% relax_initial_eps_nuc_factor .and. .not. restart)) then
            write(*, 1) 'relax_eps_nuc_factor', s% job% new_eps_nuc_factor
            call star_relax_eps_nuc_factor( &
               id, s% job% new_eps_nuc_factor, s% job% eps_nuc_factor_multiplier, ierr)
            if (failed('star_relax_eps_nuc_factor',ierr)) return
         end if

         if (s% job% relax_opacity_max .or. &
               (s% job% relax_initial_opacity_max .and. .not. restart)) then
            write(*, 1) 'relax_opacity_max', s% job% new_opacity_max
            call star_relax_opacity_max( &
               id, s% job% new_opacity_max, s% job% opacity_max_multiplier, ierr)
            if (failed('star_relax_opacity_max',ierr)) return
         end if

         if (s% job% relax_max_surf_dq .or. &
               (s% job% relax_initial_max_surf_dq .and. .not. restart)) then
            write(*, 1) 'relax_max_surf_dq', s% job% new_max_surf_dq
            call star_relax_max_surf_dq( &
               id, s% job% new_max_surf_dq, s% job% max_surf_dq_multiplier, ierr)
            if (failed('star_relax_max_surf_dq',ierr)) return
         end if

         if (s% job% relax_initial_mass .and. .not. restart) then
            write(*, 1) 'relax_initial_mass to new_mass', s% job% new_mass
            call star_relax_mass(id, s% job% new_mass, s% job% lg_max_abs_mdot, ierr)
            if (failed('relax_initial_mass',ierr)) return
         end if

         if (s% job% relax_initial_mass_to_remove_H_env .and. .not. restart) then
            write(*, 1) 'relax_initial_mass_to_remove_H_env'
            call star_relax_mass_to_remove_H_env( &
               id, s% job% extra_mass_retained_by_remove_H_env, s% job% lg_max_abs_mdot, ierr)
            if (failed('relax_initial_mass_to_remove_H_env',ierr)) return
         end if

         if (s% job% relax_mass_scale .or. &
               (s% job% relax_initial_mass_scale .and. .not. restart)) then
            write(*, 1) 'relax_mass_scale', s% job% new_mass
            call star_relax_mass_scale( &
               id, s% job% new_mass, s% job% dlgm_per_step, &
               s% job% change_mass_years_for_dt, ierr)
            if (failed('star_relax_mass_scale',ierr)) return
         end if

         if (s% job% relax_core .or. &
               (s% job% relax_initial_core .and. .not. restart)) then
            write(*, 1) 'relax_core', s% job% new_core_mass
            call star_relax_core( &
               id, s% job% new_core_mass, s% job% dlg_core_mass_per_step, &
               s% job% relax_core_years_for_dt, &
               s% job% core_avg_rho, s% job% core_avg_eps, ierr)
            if (failed('star_relax_core',ierr)) return
         end if

         call do_remove_center(id, s, restart, ierr)
         if (ierr /= 0) return
         
         if (s% job% relax_M_center .or. &
               (s% job% relax_initial_M_center .and. .not. restart)) then
            write(*, 1) 'relax_M_center', s% job% new_mass
            call star_relax_M_center( &
               id, s% job% new_mass, s% job% dlgm_per_step, s% job% relax_M_center_dt, ierr)
            if (failed('star_relax_M_center',ierr)) return
         end if
         
         if (s% job% relax_R_center .or. &
               (s% job% relax_initial_R_center .and. .not. restart)) then
            write(*, 1) 'relax_R_center', s% job% new_R_center
            call star_relax_R_center( &
               id, s% job% new_R_center, s% job% dlgR_per_step, s% job% relax_R_center_dt, ierr)
            if (failed('star_relax_R_center',ierr)) return
         end if
         
         if (s% job% relax_v_center .or. &
               (s% job% relax_initial_v_center .and. .not. restart)) then
            write(*, 1) 'relax_v_center', s% job% new_v_center
            call star_relax_v_center( &
               id, s% job% new_v_center, s% job% dv_per_step, s% job% relax_v_center_dt, ierr)
            if (failed('star_relax_v_center',ierr)) return
         end if
         
         if (s% job% relax_L_center .or. &
               (s% job% relax_initial_L_center .and. .not. restart)) then
            write(*, 1) 'relax_L_center', s% job% new_L_center
            call star_relax_L_center( &
               id, s% job% new_L_center, s% job% dlgL_per_step, s% job% relax_L_center_dt, ierr)
            if (failed('star_relax_L_center',ierr)) return
         end if
         
         if (s% job% relax_Tsurf_factor .or. &
               (s% job% relax_initial_Tsurf_factor .and. .not. restart)) then
            write(*,1) 'relax_Tsurf_factor', s% job% relax_to_this_Tsurf_factor
            call star_relax_Tsurf_factor( &
               id, s% job% relax_to_this_Tsurf_factor, s% job% dlogTsurf_factor, ierr)
            if (failed('star_relax_Tsurf_factor',ierr)) return
         end if
         
         if (s% job% relax_tau_factor .or. &
               (s% job% relax_initial_tau_factor .and. .not. restart)) then
            write(*,1) 'relax_tau_factor', s% job% relax_to_this_tau_factor
            call star_relax_tau_factor( &
               id, s% job% relax_to_this_tau_factor, s% job% dlogtau_factor, ierr)
            if (failed('star_relax_tau_factor',ierr)) return
         end if
         
         if (s% job% relax_opacity_factor .or. &
               (s% job% relax_initial_opacity_factor .and. .not. restart)) then
            write(*,1) 'relax_opacity_factor', s% job% relax_to_this_opacity_factor
            call star_relax_opacity_factor( &
               id, s% job% relax_to_this_opacity_factor, s% job% d_opacity_factor, ierr)
            if (failed('star_relax_opacity_factor',ierr)) return
         end if

         if (s% job% relax_irradiation .or. &
               (s% job% relax_initial_irradiation .and. .not. restart)) then
            write(*,2) 'relax_irradiation -- min steps', s% job% relax_irradiation_min_steps
            write(*,1) 'relax_irradiation -- max yrs dt', s% job% relax_irradiation_max_yrs_dt
            call star_relax_irradiation(id, &
               s% job% relax_irradiation_min_steps, &
               s% job% relax_to_this_irrad_flux, s% job% irrad_col_depth, &
               s% job% relax_irradiation_max_yrs_dt, ierr)
            if (failed('star_relax_irradiation',ierr)) return
         end if

         if (s% job% relax_mass_change .or. &
               (s% job% relax_initial_mass_change .and. .not. restart)) then
            write(*,2) 'relax_mass_change -- min steps', &
               s% job% relax_mass_change_min_steps
            write(*,1) 'relax_mass_change -- max yrs dt', &
               s% job% relax_mass_change_max_yrs_dt
            write(*,1) 'relax_mass_change -- initial_mass_change', &
               s% job% relax_mass_change_init_mdot
            write(*,1) 'relax_mass_change -- final_mass_change', &
               s% job% relax_mass_change_final_mdot
            call star_relax_mass_change(id, &
               s% job% relax_mass_change_min_steps, &
               s% job% relax_mass_change_init_mdot, &
               s% job% relax_mass_change_final_mdot, &
               s% job% relax_mass_change_max_yrs_dt, ierr)
            if (failed('star_relax_mass_change',ierr)) return
         end if
         
         call do_remove_initial_surface(id, s, restart, ierr)
         if (ierr /= 0) return
         
         call do_remove_surface(id, s, ierr)
         if (ierr /= 0) return

         if (s% rotation_flag .and. s% job% relax_omega) then
            write(*,1) 'new_omega', s% job% new_omega
            call star_relax_uniform_omega( &
               id, relax_to_new_omega, &
               s% job% new_omega, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% relax_initial_omega .and. .not. restart) then
            call star_relax_uniform_omega( &
               id, relax_to_new_omega, &
               s% job% new_omega, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            write(*,1) 'new_omega', s% job% new_omega
         end if

         if (s% rotation_flag .and. s% job% relax_omega_div_omega_crit) then
            if (failed('star_surface_omega_crit',ierr)) return
            call star_relax_uniform_omega( &
               id, relax_to_new_omega_div_omega_crit, &
               s% job% new_omega_div_omega_crit, &
               s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit
         end if

         if (s% rotation_flag .and. &
               s% job% relax_initial_omega_div_omega_crit .and. .not. restart) then
            if (failed('star_surface_omega_crit',ierr)) return
            call star_relax_uniform_omega( &
               id, relax_to_new_omega_div_omega_crit, &
               s% job% new_omega_div_omega_crit, &
               s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit
         end if

         if (s% rotation_flag .and. s% job% relax_surface_rotation_v) then
            call star_relax_uniform_omega( &
               id, relax_to_new_surface_rotation_v, &
               s% job% new_surface_rotation_v, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            s% job% new_omega = s% job% new_surface_rotation_v*1d5/s% r(1)
            write(*,1) 'new_surface_rotation_v', &
               s% job% new_surface_rotation_v, s% job% new_omega
         end if

         if (s% rotation_flag .and. &
               s% job% relax_initial_surface_rotation_v .and. .not. restart) then
            write(*,1) 'new_omega', s% job% new_omega
            write(*,*) 'call star_relax_uniform_omega'
            call star_relax_uniform_omega( &
               id, relax_to_new_surface_rotation_v, &
               s% job% new_surface_rotation_v, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            write(*,2) 'new_surface_rotation_v', &
               s% model_number, s% job% new_surface_rotation_v
         end if

        if (s% job% set_max_dt_to_frac_lifetime) then
           log_m = log10(s% star_mass) ! in Msun units
           log_lifetime = 9.921d0 - (3.6648d0 + (1.9697d0 - 0.9369d0*log_m)*log_m)*log_m
           ! Iben & Laughlin (1989) as quoted in H&K (eqn 2.3)
           max_dt = s% job% max_frac_of_lifetime_per_step*secyer*exp10(log_lifetime)
           if (max_dt < s% max_timestep) then
              s% max_timestep = max_dt
              write(*, *) 'set_max_dt_to_frac_lifetime: lg(maxdt/secyer)', &
                 log10(s% max_timestep/secyer)
           end if
        end if
         
         ! print out info about selected non-standard parameter settings
         
         write(*,*) 'net name ' // trim(s% net_name)

         if (s% do_element_diffusion) &
            write(*,*) 'do_element_diffusion', s% do_element_diffusion
         
         if (s% RSP_flag) &
            write(*,*) 'RSP_flag', s% RSP_flag
         
         if (s% v_flag) &
            write(*,*) 'v_flag', s% v_flag
         
         if (s% u_flag) &
            write(*,*) 'u_flag', s% u_flag
         
         if (s% rotation_flag) &
            write(*,*) 'rotation_flag', s% rotation_flag
         
         if (s% conv_vel_flag) &
            write(*,*) 'conv_vel_flag', s% conv_vel_flag
         
         if (s% w_div_wc_flag) &
            write(*,*) 'w_div_wc_flag', s% w_div_wc_flag
         
         if (s% j_rot_flag) &
            write(*,*) 'j_rot_flag', s% j_rot_flag
         
         if (s% mix_factor /= 1d0) &
            write(*,1) 'mix_factor', s% mix_factor
            
         if (abs(s% tau_base - 2d0/3d0) > 1d-4) &
            write(*,1) 'tau_base', s% tau_base
            
         if (abs(s% tau_factor - 1) > 1d-4) &
            write(*,1) 'tau_factor', s% tau_factor
            
         if (s% eps_grav_factor /= 1) &
            write(*,1) 'eps_grav_factor', s% eps_grav_factor
            
         if (s% eps_mdot_factor /= 1) &
            write(*,1) 'eps_mdot_factor', s% eps_mdot_factor

         if (s% dxdt_nuc_factor /= 1) &
            write(*,1) 'dxdt_nuc_factor', s% dxdt_nuc_factor
            
         if (.NOT. ( &
              s% atm_option == 'T_tau' .AND. &
              s% atm_T_tau_relation == 'Eddington' .AND. &
              s% atm_T_tau_opacity == 'fixed')) &
            write(*,1) 'atm_option: ' // trim(s% atm_option)
           
         if (s% M_center /= 0) then
            write(*,1) 'xmstar/mstar', s% xmstar/s% mstar
            write(*,1) 'xmstar (g)', s% xmstar
            write(*,1) 'M_center (g)', s% M_center
            write(*,1) 'xmstar/Msun', s% xmstar/Msun
            write(*,1) 'M_center/Msun', s% M_center/Msun
         end if
         
         if (s% v_flag .or. s% u_flag) then
            if (s% v_center /= 0) &
               write(*,1) 'v_center (cm/s)', s% v_center
         end if
            
         if (s% R_center /= 0) then
            write(*,1) 'R_center (cm)', s% R_center
            write(*,1) 'R_center/Rsun', s% R_center/Rsun
            write(*,1) 'core density', &
               s% M_center/(4*pi/3*s% R_center*s% R_center*s% R_center)
         end if
         
         if (s% L_center /= 0) &
            write(*,1) 'L_center/Lsun', s% L_center/Lsun
                     
         if (s% opacity_max > 0) &
            write(*,1) 'opacity_max', s% opacity_max
         
         if (s% job% show_net_reactions_info) then
            write(*,'(a)') ' net reactions '
            call show_net_reactions_and_info(s% net_handle, 6, ierr)
            if (failed('show_net_reactions_and_info',ierr)) return
         end if
         
         if (s% job% list_net_reactions) then
            write(*,'(a)') ' net reactions '
            call show_net_reactions(s% net_handle, 6, ierr)
            if (failed('show_net_reactions',ierr)) return
         end if
         
         if (s% job% set_cumulative_energy_error .or. &
               (s% job% set_initial_cumulative_energy_error .and. .not. restart) .or. &
               (s% model_number == s% job% set_cumulative_energy_error_at_step)) then
            write(*,1) 'set_cumulative_energy_error', s% job% new_cumulative_energy_error
            s% cumulative_energy_error = s% job% new_cumulative_energy_error
         end if
         
         if (s% job% show_net_species_info) then
            write(*,'(a)') ' species'
            do j=1,s% species
               write(*,'(i6,3x,a)') j, chem_isos% name(s% chem_id(j))
            end do
            write(*,*)
         end if
         
         if (s% job% show_eqns_and_vars_names) then
            do i=1,s% nvar_total
               write(*,*) i, s% nameofvar(i), s% nameofequ(i)
            end do
            write(*,*)
         end if         
         
         write(*,*) 'kap_option ' // trim(kap_option_str(s% kap_rq% kap_option))
         write(*,*) 'kap_CO_option ' // trim(kap_CO_option_str(s% kap_rq% kap_CO_option))
         write(*,*) 'kap_lowT_option ' // trim(kap_lowT_option_str(s% kap_rq% kap_lowT_option))
         write(*,2) 'OMP_NUM_THREADS', utils_omp_get_max_threads()
               
         call check_if_want_to_stop_warnings(s)
         
         contains

         subroutine do_relax_initial_composition(ierr)
            use utils_lib
            integer, intent(out) :: ierr
            real(dp), pointer :: xq(:), xa(:,:)
            integer :: num_pts, num_species, i, iounit
            include 'formats'
            
            write(*,*)
            write(*,1) 'relax_initial_composition'

            open(newunit=iounit, file=trim(s% job% relax_composition_filename), &
                  status='old', action='read', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'open failed', ierr, iounit
               write(*, '(a)') 'failed to open ' // trim(s% job% relax_composition_filename)
               return
            end if
            read(iounit, *, iostat=ierr) num_pts, num_species
            if (ierr /= 0) then
               close(iounit)
               write(*, '(a)') 'failed while trying to read 1st line of ' // &
                  trim(s% job% relax_composition_filename)
               return
            end if
            if(num_species .ne. s% species) then
               write(*,*) 'Error in ',trim(s% job% relax_composition_filename)
               write(*,'(a,I4,a)') 'got ',num_species,' species'
               write(*,'(a,I4,a)') 'expected ', s% species,' species'
               write(*,*)
               ierr=-1
               return
            end if
            allocate(xq(num_pts), xa(num_species,num_pts))
            do i = 1, num_pts
               read(iounit,*,iostat=ierr) xq(i), xa(1:num_species,i)
               if (ierr /= 0) then
                  close(iounit)
                  write(*, '(a)') &
                     'failed while trying to read ' // trim(s% job% relax_composition_filename)
                  write(*,*) 'line', i+1
                  write(*,*) 'perhaps wrong info in 1st line?'
                  write(*,*) '1st line must have num_pts and num_species in that order'
                  deallocate(xq,xa)
                  return
               end if
            end do
            close(iounit)
            
            call star_relax_composition( &
               id, s% job% num_steps_to_relax_composition, num_pts, num_species, xa, xq, ierr)
            deallocate(xq,xa)
            
         end subroutine do_relax_initial_composition
         
         subroutine do_relax_initial_angular_momentum(ierr)
            use utils_lib
            integer, intent(out) :: ierr
            real(dp), pointer :: xq(:), angular_momentum(:)
            integer :: num_pts, i, iounit
            include 'formats'
            
            write(*,*)
            write(*,1) 'relax_initial_angular_momentum'

            open(newunit=iounit, file=trim(s% job% relax_angular_momentum_filename), &
                  status='old', action='read', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'open failed', ierr, iounit
               write(*, '(a)') 'failed to open "' // trim(s% job% relax_angular_momentum_filename)//'"'
               return
            end if
            read(iounit, *, iostat=ierr) num_pts
            if (ierr /= 0) then
               close(iounit)
               write(*, '(a)') 'failed while trying to read 1st line of ' // &
                  trim(s% job% relax_angular_momentum_filename)
               return
            end if
            allocate(xq(num_pts), angular_momentum(num_pts))
            do i = 1, num_pts
               read(iounit,*,iostat=ierr) xq(i), angular_momentum(i)
               if (ierr /= 0) then
                  close(iounit)
                  write(*, '(a)') &
                     'failed while trying to read ' // trim(s% job% relax_angular_momentum_filename)
                  write(*,*) 'line', i+1
                  write(*,*) 'perhaps wrong info in 1st line?'
                  write(*,*) '1st line must have num_pts'
                  deallocate(xq,angular_momentum)
                  return
               end if
            end do
            close(iounit)
            call star_relax_angular_momentum(id, s% job% max_steps_to_relax_angular_momentum, &
               num_pts, angular_momentum, xq, ierr)
            deallocate(xq,angular_momentum)
         end subroutine do_relax_initial_angular_momentum
         
         subroutine do_relax_initial_entropy(ierr)
            use utils_lib
            use eos_def
            integer, intent(out) :: ierr
            ! arrays into which data from the input file is read.
            ! in case any of the eos* options is used, input from the
            ! file is read into var1 and var2, and the chosen eos function
            ! is used to extract the entropy from that pair.
            real(dp), pointer :: xq(:), entropy(:)
            real(dp) :: var1, var2
            integer :: num_pts, i, k, iounit
            ! these are needed to call eosPT_get
            real(dp) :: Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
            real(dp) :: T, log10T
            ! these are needed to call eosDT_get_T
            real(dp) :: T_guess_gas, T_guess_rad, logT_guess
            integer :: eos_calls
            ! these are used for all eos calls
            real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
            real(dp), dimension(num_eos_d_dxa_results, s% species) :: d_dxa
            real(dp), parameter :: logT_tol = 1d-8, logE_tol = 1d-8
            integer, parameter :: MAX_ITERS = 20
            include 'formats'
            
            write(*,*)
            write(*,1) 'relax_initial_entropy'

            open(newunit=iounit, file=trim(s% job% relax_entropy_filename), &
                  status='old', action='read', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'open failed', ierr, iounit
               write(*, '(a)') 'failed to open "' // trim(s% job% relax_entropy_filename)//'"'
               return
            end if
            read(iounit, *, iostat=ierr) num_pts
            if (ierr /= 0) then
               close(iounit)
               write(*, '(a)') 'failed while trying to read 1st line of ' // &
                  trim(s% job% relax_entropy_filename)
               return
            end if
            if (.not. (s% job% get_entropy_for_relax_from_eos == '' .or. &
                  s% job% get_entropy_for_relax_from_eos == 'eosDT' .or. &
                  s% job% get_entropy_for_relax_from_eos == 'eosPT' .or. &
                  s% job% get_entropy_for_relax_from_eos == 'eosDE')) then
               ierr = 1
               write(*,*) 'invalid value for get_entropy_for_relax_from_eos =', &
                  s% job% get_entropy_for_relax_from_eos
            end if
            allocate(xq(num_pts), entropy(num_pts))
            do i = 1, num_pts
               if (s% job% get_entropy_for_relax_from_eos == '') then
                  read(iounit,*,iostat=ierr) xq(i), entropy(i)
               else
                  read(iounit,*,iostat=ierr) xq(i), var1, var2
                  ! get nearest value matching xq for the composition TODO: interpolate
                  do k=1, s% nz-1
                     if(1-s% q(k) <= xq(i) .and. 1-s% q(k+1) >= xq(i)) then
                        exit
                     end if
                  end do
                  ! get entropy
                  if (s% job% get_entropy_for_relax_from_eos == 'eosDT') then
                     call eosDT_get( &
                        s% eos_handle, &
                        s% species, s% chem_id, s% net_iso, s% xa(:,k), &
                        var1, log10(var1), var2, log10(var2), &
                        res, d_dlnd, d_dlnT, d_dxa, ierr)
                     if (ierr /= 0) then
                        write(*,*) "failed in eosDT_get"
                        return
                     end if
                     entropy(i) = exp(res(i_lnS))
                  else if (s% job% get_entropy_for_relax_from_eos == 'eosPT') then
                     call eosPT_get( &
                        s% eos_handle, &
                        s% species, s% chem_id, s% net_iso, s% xa(:,k), &
                        var1, log10(var1), var2, log10(var2), &
                        Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
                        res, d_dlnd, d_dlnT, d_dxa, ierr)
                     if (ierr /= 0) then
                        write(*,*) "failed in eosPT_get"
                        return
                     end if
                     entropy(i) = exp(res(i_lnS))
                  else
                     T_guess_gas = 2*var2*s% abar(k)*mp/(3*kerg*(1+s% zbar(k))) ! ideal gas (var2=energy)
                     T_guess_rad = pow(var2/crad,0.25d0)
                     logT_guess = log10(max(T_guess_gas,T_guess_rad))
                     call eosDT_get_T( &
                        s% eos_handle, &
                        s% species, s% chem_id, s% net_iso, s% xa(:,k), &
                        log10(var1), i_lnE, log10(var2)*ln10, &
                        logT_tol, logE_tol*ln10, MAX_ITERS, logT_guess, &
                        arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
                        log10T, res, d_dlnd, d_dlnT, d_dxa, &
                        eos_calls, ierr)
                     if (ierr /= 0) then
                        write(*,*) "failed in eosDT_get_T (as eosDE)"
                        return
                     end if
                     entropy(i) = exp(res(i_lnS))
                  end if
               end if
               if (ierr /= 0) then
                  close(iounit)
                  write(*, '(a)') &
                     'failed while trying to read ' // trim(s% job% relax_entropy_filename)
                  write(*,*) 'line', i+1
                  write(*,*) 'perhaps wrong info in 1st line?'
                  write(*,*) '1st line must have num_pts'
                  deallocate(xq,entropy)
                  return
               end if
            end do
            close(iounit)
            call star_relax_entropy(id, s% job% max_steps_to_relax_entropy, num_pts, entropy, xq, ierr)
            deallocate(xq,entropy)
         end subroutine do_relax_initial_entropy

      end subroutine do_star_job_controls_after
       

      subroutine do_remove_center(id, s, restart, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         include 'formats'
         
         if (s% job% remove_center_by_temperature > 0) then
            write(*, 1) 'remove_center_by_temperature', s% job% remove_center_by_temperature
            call star_remove_center_by_temperature( &
               id, s% job% remove_center_by_temperature, ierr)
            if (failed('star_remove_center_by_temperature',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_temperature > 0 .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_temperature', &
               s% job% remove_initial_center_by_temperature
            call star_remove_center_by_temperature( &
               id, s% job% remove_initial_center_by_temperature, ierr)
            if (failed('star_remove_center_by_temperature',ierr)) return
         end if
         
         if (s% job% remove_center_by_radius_cm > s% R_center .and. &
               s% job% remove_center_by_radius_cm < s% r(1)) then
            write(*, 1) 'remove_center_by_radius_cm', &
               s% job% remove_center_by_radius_cm
            call star_remove_center_by_radius_cm( &
               id, s% job% remove_center_by_radius_cm, ierr)
            if (failed('star_remove_center_by_radius_cm',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_radius_cm > s% R_center .and. &
               s% job% remove_initial_center_by_radius_cm < s% r(1) .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_radius_cm', &
               s% job% remove_initial_center_by_radius_cm
            call star_remove_center_by_radius_cm( &
               id, s% job% remove_initial_center_by_radius_cm, ierr)
            if (failed('star_remove_center_by_radius_cm',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_he4 > 0d0 .and. &
               s% job% remove_initial_center_by_he4 < 1d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_he4', &
               s% job% remove_initial_center_by_he4
            call star_remove_center_by_he4( &
               id, s% job% remove_initial_center_by_he4, ierr)
            if (failed('star_remove_initial_center_by_he4',ierr)) return
         end if
         
         if (s% job% remove_center_by_he4 > 0d0 .and. &
               s% job% remove_center_by_he4 < 1d0) then
            write(*, 1) 'remove_center_by_he4', &
               s% job% remove_center_by_he4
            call star_remove_center_by_he4( &
               id, s% job% remove_center_by_he4, ierr)
            if (failed('star_remove_center_by_he4',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_c12_o16 > 0d0 .and. &
               s% job% remove_initial_center_by_c12_o16 < 1d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_c12_o16', &
               s% job% remove_initial_center_by_c12_o16
            call star_remove_center_by_c12_o16( &
               id, s% job% remove_initial_center_by_c12_o16, ierr)
            if (failed('star_remove_initial_center_by_c12_o16',ierr)) return
         end if
         
         if (s% job% remove_center_by_c12_o16 > 0d0 .and. &
               s% job% remove_center_by_c12_o16 < 1d0) then
            write(*, 1) 'remove_center_by_c12_o16', &
               s% job% remove_center_by_c12_o16
            call star_remove_center_by_c12_o16( &
               id, s% job% remove_center_by_c12_o16, ierr)
            if (failed('star_remove_center_by_c12_o16',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_si28 > 0d0 .and. &
               s% job% remove_initial_center_by_si28 < 1d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_si28', &
               s% job% remove_initial_center_by_si28
            call star_remove_center_by_si28( &
               id, s% job% remove_initial_center_by_si28, ierr)
            if (failed('star_remove_initial_center_by_si28',ierr)) return
         end if
         
         if (s% job% remove_center_by_si28 > 0d0 .and. &
               s% job% remove_center_by_si28 < 1d0) then
            write(*, 1) 'remove_center_by_si28', &
               s% job% remove_center_by_si28
            call star_remove_center_by_si28( &
               id, s% job% remove_center_by_si28, ierr)
            if (failed('star_remove_center_by_si28',ierr)) return
         end if
         
         if (s% job% remove_initial_center_to_reduce_co56_ni56 > 0d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_to_reduce_co56_ni56', &
               s% job% remove_initial_center_to_reduce_co56_ni56
            call star_remove_center_to_reduce_co56_ni56( &
               id, s% job% remove_initial_center_to_reduce_co56_ni56, ierr)
            if (failed('star_remove_initial_center_to_reduce_co56_ni56',ierr)) return
         end if
         
         if (s% job% remove_center_to_reduce_co56_ni56 > 0d0) then
            write(*, 1) 'remove_center_to_reduce_co56_ni56', &
               s% job% remove_center_to_reduce_co56_ni56
            call star_remove_center_to_reduce_co56_ni56( &
               id, s% job% remove_center_to_reduce_co56_ni56, ierr)
            if (failed('star_remove_center_to_reduce_co56_ni56',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_ye > 0d0 .and. &
               s% job% remove_initial_center_by_ye < 1d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_ye', &
               s% job% remove_initial_center_by_ye
            call star_remove_center_by_ye( &
               id, s% job% remove_initial_center_by_ye, ierr)
            if (failed('star_remove_initial_center_by_ye',ierr)) return
         end if
         
         if (s% job% remove_center_by_ye > 0d0 .and. &
               s% job% remove_center_by_ye < 1d0) then
            write(*, 1) 'remove_center_by_ye', &
               s% job% remove_center_by_ye
            call star_remove_center_by_ye( &
               id, s% job% remove_center_by_ye, ierr)
            if (failed('star_remove_center_by_ye',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_entropy > 0d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_entropy', &
               s% job% remove_initial_center_by_entropy
            call star_remove_center_by_entropy( &
               id, s% job% remove_initial_center_by_entropy, ierr)
            if (failed('star_remove_initial_center_by_entropy',ierr)) return
         end if
         
         if (s% job% remove_center_by_entropy > 0d0) then
            write(*, 1) 'remove_center_by_entropy', &
               s% job% remove_center_by_entropy
            call star_remove_center_by_entropy( &
               id, s% job% remove_center_by_entropy, ierr)
            if (failed('star_remove_center_by_entropy',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_infall_kms /= 0d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_infall_kms', &
               s% job% remove_initial_center_by_infall_kms
            call star_remove_center_by_infall_kms( &
               id, s% job% remove_initial_center_by_infall_kms, ierr)
            if (failed('star_remove_initial_center_by_infall_kms',ierr)) return
         end if
         
         if (s% job% remove_center_by_infall_kms /= 0d0) then
            write(*, 1) 'remove_center_by_infall_kms', &
               s% job% remove_center_by_infall_kms
            call star_remove_center_by_infall_kms( &
               id, s% job% remove_center_by_infall_kms, ierr)
            if (failed('star_remove_center_by_infall_kms',ierr)) return
         end if
         
         if (s% job% remove_initial_center_at_inner_max_abs_v &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_at_inner_max_abs_v'
            call star_remove_center_at_inner_max_abs_v(id, ierr)
            if (failed('remove_center_at_inner_max_abs_v',ierr)) return
         end if
         
         if (s% job% remove_center_at_inner_max_abs_v) then
            write(*, 1) 'remove_initial_center_at_inner_max_abs_v'
            call star_remove_center_at_inner_max_abs_v(id, ierr)
            if (failed('remove_center_at_inner_max_abs_v',ierr)) return
         end if
         
         if (s% job% remove_initial_fe_core .and. .not. restart) then
            write(*, 1) 'remove_initial_fe_core'
            call star_remove_fe_core(id, ierr)
            if (failed('remove_fe_core',ierr)) return
         end if
         
         if (s% job% remove_fe_core) then
            write(*, 1) 'remove_initial_fe_core'
            call star_remove_fe_core(id, ierr)
            if (failed('remove_fe_core',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_mass_fraction_q > 0d0 .and. &
               s% job% remove_initial_center_by_mass_fraction_q < 1d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_mass_fraction_q', &
               s% job% remove_initial_center_by_mass_fraction_q
            call star_remove_center_by_mass_fraction_q( &
               id, s% job% remove_initial_center_by_mass_fraction_q, ierr)
            if (failed('star_remove_initial_center_by_mass_fraction_q',ierr)) return
         end if
         
         if (s% job% remove_center_by_mass_fraction_q > 0d0 .and. &
               s% job% remove_center_by_mass_fraction_q < 1d0) then
            write(*, 1) 'remove_center_by_mass_fraction_q', &
               s% job% remove_center_by_mass_fraction_q
            call star_remove_center_by_mass_fraction_q( &
               id, s% job% remove_center_by_mass_fraction_q, ierr)
            if (failed('star_remove_center_by_mass_fraction_q',ierr)) return
         end if
         
         if (s% job% remove_center_by_delta_mass_gm > 0) then
            write(*, 1) 'remove_center_by_delta_mass_gm', &
               s% job% remove_center_by_delta_mass_gm
            call star_remove_center_by_mass_gm(id, &
               s% M_center + s% job% remove_center_by_delta_mass_gm, ierr)
            if (failed('star_remove_center_by_mass',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_delta_mass_gm > 0 .and. &
               .not. restart) then
            write(*, 1) 'remove_initial_center_by_delta_mass_gm', &
               s% job% remove_initial_center_by_delta_mass_gm
            call star_remove_center_by_mass_gm(id, &
               s% M_center + s% job% remove_initial_center_by_delta_mass_gm, ierr)
            if (failed('star_remove_center_by_mass',ierr)) return
         end if
         
         if (s% job% remove_center_by_delta_mass_Msun > 0) then
            write(*, 1) 'remove_center_by_delta_mass_Msun', &
               s% job% remove_center_by_delta_mass_Msun
            call star_remove_center_by_mass_gm(id, &
               s% M_center + s% job% remove_center_by_delta_mass_Msun*Msun, ierr)
            if (failed('star_remove_center_by_mass',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_delta_mass_Msun > 0 .and. &
               .not. restart) then
            write(*, 1) 'remove_initial_center_by_delta_mass_Msun', &
               s% job% remove_initial_center_by_delta_mass_Msun
            call star_remove_center_by_mass_gm(id, &
               s% M_center + s% job% remove_initial_center_by_delta_mass_Msun*Msun, ierr)
            if (failed('star_remove_center_by_mass',ierr)) return
         end if
         
         if (s% job% remove_center_by_mass_gm > s% M_center .and. &
               s% job% remove_center_by_mass_gm < s% m(1)) then
            write(*, 1) 'remove_center_by_mass_gm', &
               s% job% remove_center_by_mass_gm
            call star_remove_center_by_mass_gm( &
               id, s% job% remove_center_by_mass_gm, ierr)
            if (failed('star_remove_center_by_mass_gm',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_mass_gm > s% M_center .and. &
               s% job% remove_initial_center_by_mass_gm < s% m(1) .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_mass_gm', &
               s% job% remove_initial_center_by_mass_gm
            call star_remove_center_by_mass_gm( &
               id, s% job% remove_initial_center_by_mass_gm, ierr)
            if (failed('star_remove_center_by_mass_gm',ierr)) return
         end if
         
         if (s% job% remove_center_by_mass_Msun > s% M_center/Msun .and. &
               s% job% remove_center_by_mass_Msun < s% m(1)/Msun) then
            write(*, 1) 'remove_center_by_mass_Msun', &
               s% job% remove_center_by_mass_Msun
            call star_remove_center_by_mass_gm( &
               id, s% job% remove_center_by_mass_Msun*Msun, ierr)
            if (failed('star_remove_center_by_mass_Msun',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_mass_Msun > s% M_center/Msun .and. &
               s% job% remove_initial_center_by_mass_Msun < s% m(1)/Msun .and. &
               .not. restart) then
            write(*, 1) 'remove_initial_center_by_mass_Msun', &
               s% job% remove_initial_center_by_mass_Msun
            call star_remove_center_by_mass_gm( &
               id, s% job% remove_initial_center_by_mass_Msun*Msun, ierr)
            if (failed('star_remove_center_by_mass_Msun',ierr)) return
         end if
         
         if (s% job% remove_center_by_radius_Rsun > s% R_center/Rsun .and. &
               s% job% remove_center_by_radius_Rsun < s% r(1)/Rsun) then
            write(*, 1) 'remove_center_by_radius_Rsun', &
               s% job% remove_center_by_radius_Rsun
            call star_remove_center_by_radius_cm( &
               id, s% job% remove_center_by_radius_Rsun*Rsun, ierr)
            if (failed('star_remove_center_by_radius_Rsun',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_radius_Rsun > s% R_center/Rsun .and. &
               s% job% remove_initial_center_by_radius_Rsun < s% r(1)/Rsun .and. &
               .not. restart) then
            write(*, 1) 'remove_initial_center_by_radius_Rsun', &
               s% job% remove_initial_center_by_radius_Rsun
            call star_remove_center_by_radius_cm( &
               id, s% job% remove_initial_center_by_radius_Rsun*Rsun, ierr)
            if (failed('star_remove_center_by_radius_Rsun',ierr)) return
         end if

         if (s% job% remove_initial_center_at_cell_k > 0 .and. .not. restart .and. &
               s% job% remove_initial_center_at_cell_k <= s% nz) then
            write(*, 2) 'remove_initial_center_at_cell_k', s% job% remove_initial_center_at_cell_k
            call star_remove_center_at_cell_k( &
               id, s% job% remove_initial_center_at_cell_k, ierr)
            if (failed('star_remove_center_at_cell_k',ierr)) return
         end if

         if (s% job% remove_center_at_cell_k > 0 .and. &
               s% job% remove_center_at_cell_k <= s% nz) then
            write(*, 2) 'remove_center_at_cell_k', s% job% remove_center_at_cell_k
            call star_remove_center_at_cell_k(id, s% job% remove_center_at_cell_k, ierr)
            if (failed('star_remove_center_at_cell_k',ierr)) return
         end if
         
      end subroutine do_remove_center
       

      subroutine do_remove_initial_surface(id,s,restart,ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr

         include 'formats'
         
         ierr = 0

         if (s% job% remove_initial_surface_at_he_core_boundary > 0 .and. .not. restart) then
            write(*, 1) 'remove_initial_surface_at_he_core_boundary', &
               s% job% remove_initial_surface_at_he_core_boundary
            call star_remove_surface_at_he_core_boundary( &
               id, s% job% remove_initial_surface_at_he_core_boundary, ierr)
            if (failed('star_remove_surface_at_he_core_boundary',ierr)) return
         end if

         if (s% job% remove_initial_surface_by_optical_depth > 0 .and. .not. restart) then
            write(*, 1) 'remove_initial_surface_by_optical_depth', &
               s% job% remove_initial_surface_by_optical_depth
            call star_remove_surface_by_optical_depth( &
               id, s% job% remove_initial_surface_by_optical_depth, ierr)
            if (failed('star_remove_surface_by_optical_depth',ierr)) return
         end if
         
         if (s% job% remove_initial_surface_by_density > 0 .and. .not. restart) then
            write(*, 1) 'remove_initial_surface_by_density', &
               s% job% remove_initial_surface_by_density
            call star_remove_surface_by_density( &
               id, s% job% remove_initial_surface_by_density, ierr)
            if (failed('star_remove_surface_by_density',ierr)) return
         end if
         
         if (s% job% remove_initial_surface_by_pressure > 0 .and. .not. restart) then
            write(*, 1) 'remove_initial_surface_by_pressure', &
               s% job% remove_initial_surface_by_pressure
            call star_remove_surface_by_pressure( &
               id, s% job% remove_initial_surface_by_pressure, ierr)
            if (failed('star_remove_surface_by_pressure',ierr)) return
         end if
         
         if (s% job% remove_initial_surface_by_radius_cm > s% R_center .and. &
               s% job% remove_initial_surface_by_radius_cm < s% r(1) .and. .not. restart) then
            write(*, 1) 'remove_initial_surface_by_radius_cm', &
               s% job% remove_initial_surface_by_radius_cm
            call star_remove_surface_by_radius_cm( &
               id, s% job% remove_initial_surface_by_radius_cm, ierr)
            if (failed('star_remove_surface_by_radius_cm',ierr)) return
         end if
         
         if (s% job% remove_initial_surface_by_mass_fraction_q > 0d0 .and. &
               s% job% remove_initial_surface_by_mass_fraction_q < 1d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_surface_by_mass_fraction_q', &
               s% job% remove_initial_surface_by_mass_fraction_q
            call star_remove_surface_by_mass_fraction_q( &
               id, s% job% remove_initial_surface_by_mass_fraction_q, ierr)
            if (failed('star_remove_initial_surface_by_mass_fraction_q',ierr)) return
         end if
         
         if (s% job% remove_initial_surface_by_mass_gm > s% M_center .and. &
               s% job% remove_initial_surface_by_mass_gm < s% m(1) .and. .not. restart) then
            write(*, 1) 'remove_initial_surface_by_mass_gm', &
               s% job% remove_initial_surface_by_mass_gm
            call star_remove_surface_by_mass_gm( &
               id, s% job% remove_initial_surface_by_mass_gm, ierr)
            if (failed('star_remove_surface_by_mass_gm',ierr)) return
         end if
         
         if (s% job% remove_initial_surface_by_radius_Rsun > s% R_center/Rsun .and. &
               s% job% remove_initial_surface_by_radius_Rsun < s% r(1)/Rsun .and. &
               .not. restart) then
            write(*, 1) 'remove_initial_surface_by_radius_Rsun', &
               s% job% remove_initial_surface_by_radius_Rsun
            call star_remove_surface_by_radius_cm( &
               id, s% job% remove_initial_surface_by_radius_Rsun*Rsun, ierr)
            if (failed('star_remove_surface_by_radius_Rsun',ierr)) return
         end if
         
         if (s% job% remove_initial_surface_by_mass_Msun > s% M_center/Msun .and. &
               s% job% remove_initial_surface_by_mass_Msun < s% m(1)/Msun .and. &
               .not. restart) then
            write(*, 1) 'remove_initial_surface_by_mass_Msun', &
               s% job% remove_initial_surface_by_mass_Msun
            call star_remove_surface_by_mass_gm( &
               id, s% job% remove_initial_surface_by_mass_Msun*Msun, ierr)
            if (failed('star_remove_surface_by_mass_Msun',ierr)) return
         end if

         if (s% job% remove_initial_surface_by_v_surf_km_s > 0 .and. .not. restart) then
            write(*, 2) 'remove_initial_surface_by_v_surf_km_s', &
               s% job% remove_initial_surface_by_v_surf_km_s
            call star_remove_surface_by_v_surf_km_s( &
               id, s% job% remove_initial_surface_by_v_surf_km_s, ierr)
            if (failed('star_remove_surface_by_v_surf_km_s',ierr)) return
         end if

         if (s% job% remove_initial_surface_by_v_surf_div_cs > 0 .and. .not. restart) then
            write(*, 2) 'remove_initial_surface_by_v_surf_div_cs', &
               s% job% remove_initial_surface_by_v_surf_div_cs
            call star_remove_surface_by_v_surf_div_cs( &
               id, s% job% remove_initial_surface_by_v_surf_div_cs, ierr)
            if (failed('star_remove_surface_by_v_surf_div_cs',ierr)) return
         end if

         if (s% job% remove_initial_surface_by_v_surf_div_v_escape > 0 .and. .not. restart) then
            write(*, 2) 'remove_initial_surface_by_v_surf_div_v_escape', &
               s% job% remove_initial_surface_by_v_surf_div_v_escape
            call star_remove_surface_by_v_surf_div_v_escape( &
               id, s% job% remove_initial_surface_by_v_surf_div_v_escape, ierr)
            if (failed('star_remove_surface_by_v_surf_div_v_escape',ierr)) return
         end if

         if (s% job% remove_initial_surface_at_cell_k > 0 .and. .not. restart .and. &
               s% job% remove_initial_surface_at_cell_k <= s% nz) then
            write(*, 2) 'remove_initial_surface_at_cell_k', s% job% remove_initial_surface_at_cell_k
            call star_remove_surface_at_cell_k( &
               id, s% job% remove_initial_surface_at_cell_k, ierr)
            if (failed('star_remove_surface_at_cell_k',ierr)) return
         end if

      end subroutine do_remove_initial_surface
             

      subroutine do_remove_surface(id,s,ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         include 'formats'
         
         ierr = 0

         if (s% job% remove_surface_at_he_core_boundary > 0) then
            !write(*, 1) 'remove_surface_at_he_core_boundary', s% job% remove_surface_at_he_core_boundary
            call star_remove_surface_at_he_core_boundary( &
               id, s% job% remove_surface_at_he_core_boundary, ierr)
            if (failed('star_remove_surface_at_he_core_boundary',ierr)) return
         end if

         if (s% job% remove_surface_by_optical_depth > 0) then
            !write(*, 1) 'remove_surface_by_optical_depth', s% job% remove_surface_by_optical_depth
            call star_remove_surface_by_optical_depth( &
               id, s% job% remove_surface_by_optical_depth, ierr)
            if (failed('star_remove_surface_by_optical_depth',ierr)) return
         end if
         
         if (s% job% remove_surface_by_density > 0) then
            !write(*, 1) 'remove_surface_by_density', s% job% remove_surface_by_density
            call star_remove_surface_by_density( &
               id, s% job% remove_surface_by_density, ierr)
            if (failed('star_remove_surface_by_density',ierr)) return
         end if
         
         if (s% job% remove_surface_by_pressure > 0) then
            !write(*, 1) 'remove_surface_by_pressure', s% job% remove_surface_by_pressure
            call star_remove_surface_by_pressure( &
               id, s% job% remove_surface_by_pressure, ierr)
            if (failed('star_remove_surface_by_pressure',ierr)) return
         end if
         
         if (s% job% remove_surface_by_radius_cm > s% R_center .and. &
               s% job% remove_surface_by_radius_cm < s% r(1)) then
            !write(*, 1) 'remove_surface_by_radius_cm', s% job% remove_surface_by_radius_cm
            call star_remove_surface_by_radius_cm( &
               id, s% job% remove_surface_by_radius_cm, ierr)
            if (failed('star_remove_surface_by_radius_cm',ierr)) return
         end if
         
         if (s% job% remove_surface_by_mass_fraction_q > 0d0 .and. &
               s% job% remove_surface_by_mass_fraction_q < 1d0) then
            !write(*, 1) 'remove_surface_by_mass_fraction_q', &
            !   s% job% remove_surface_by_mass_fraction_q
            call star_remove_surface_by_mass_fraction_q( &
               id, s% job% remove_surface_by_mass_fraction_q, ierr)
            if (failed('star_remove_surface_by_mass_fraction_q',ierr)) return
         end if
         
         if (s% job% remove_surface_by_mass_gm > s% M_center .and. &
               s% job% remove_surface_by_mass_gm < s% m(1)) then
            !write(*, 1) 'remove_surface_by_mass_gm', &
            !   s% job% remove_surface_by_mass_gm
            call star_remove_surface_by_mass_gm( &
               id, s% job% remove_surface_by_mass_gm, ierr)
            if (failed('star_remove_surface_by_mass_gm',ierr)) return
         end if
         
         if (s% job% remove_surface_by_radius_Rsun > s% R_center/Rsun .and. &
               s% job% remove_surface_by_radius_Rsun < s% r(1)/Rsun) then
            !write(*, 1) 'remove_surface_by_radius_Rsun', &
            !   s% job% remove_surface_by_radius_Rsun
            call star_remove_surface_by_radius_cm( &
               id, s% job% remove_surface_by_radius_Rsun*Rsun, ierr)
            if (failed('star_remove_surface_by_radius_Rsun',ierr)) return
         end if
         
         if (s% job% remove_surface_by_mass_Msun > s% M_center/Msun .and. &
               s% job% remove_surface_by_mass_Msun < s% m(1)/Msun) then
            !write(*, 1) 'remove_surface_by_mass_Msun', &
            !   s% job% remove_surface_by_mass_Msun
            call star_remove_surface_by_mass_gm( &
               id, s% job% remove_surface_by_mass_Msun*Msun, ierr)
            if (failed('star_remove_surface_by_mass_Msun',ierr)) return
         end if

         if (s% job% remove_surface_by_v_surf_km_s > 0) then
            !write(*, 2) 'remove_surface_by_v_surf_km_s', s% job% remove_surface_by_v_surf_km_s
            call star_remove_surface_by_v_surf_km_s(id, s% job% remove_surface_by_v_surf_km_s, ierr)
            if (failed('star_remove_surface_by_v_surf_km_s',ierr)) return
         end if

         if (s% job% remove_surface_by_v_surf_div_cs > 0) then
            !write(*, 1) 'remove_surface_by_v_surf_div_cs', s% job% remove_surface_by_v_surf_div_cs
            call star_remove_surface_by_v_surf_div_cs(id, s% job% remove_surface_by_v_surf_div_cs, ierr)
            if (failed('star_remove_surface_by_v_surf_div_cs',ierr)) return
         end if

         if (s% job% remove_surface_by_v_surf_div_v_escape > 0) then
            !write(*, 2) 'remove_surface_by_v_surf_div_v_escape', s% job% remove_surface_by_v_surf_div_v_escape
            call star_remove_surface_by_v_surf_div_v_escape(id, s% job% remove_surface_by_v_surf_div_v_escape, ierr)
            if (failed('star_remove_surface_by_v_surf_div_v_escape',ierr)) return
         end if

         if (s% job% remove_surface_at_cell_k > 0 .and. &
               s% job% remove_surface_at_cell_k <= s% nz) then
            !write(*, 2) 'remove_surface_at_cell_k', s% job% remove_surface_at_cell_k
            call star_remove_surface_at_cell_k(id, s% job% remove_surface_at_cell_k, ierr)
            if (failed('star_remove_surface_at_cell_k',ierr)) return
         end if

      end subroutine do_remove_surface


      subroutine resolve_inlist_fname(inlist_out,inlist_opt)

        use ISO_FORTRAN_ENV

        character(len=*),intent(out) :: inlist_out
        character(len=*),optional   :: inlist_opt

        integer :: status

        ! initialize inlist_out as empty
        inlist_out = ''

         if (.not. MESA_INLIST_RESOLVED) then
            if (COMMAND_ARGUMENT_COUNT() >= 1) then

              ! Get filename from the first command-line argument

              call GET_COMMAND_ARGUMENT(1, inlist_out, STATUS=status)
              if (status /= 0) inlist_out = ''

            else

              ! Get filename from the MESA_INLIST environment variable

              call GET_ENVIRONMENT_VARIABLE('MESA_INLIST', inlist_out, STATUS=status)
              if (status /= 0) inlist_out = ''

            endif
         end if

        if (inlist_out == '') then

           if (PRESENT(inlist_opt)) then
              inlist_out = inlist_opt
           else
              inlist_out = 'inlist'
           endif

        endif

        ! Finish

        return

      end subroutine resolve_inlist_fname
      
      
      subroutine add_fpe_checks(id, s, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         character(len=1) fpe_check
         integer :: status

         include 'formats'

         ierr = 0

         call GET_ENVIRONMENT_VARIABLE('MESA_FPE_CHECKS_ON', fpe_check, STATUS=status)
         if (status /= 0) return

         if (fpe_check(1:1)=="1") then
            write(*,*) "FPE checking is on"
            s% fill_arrays_with_nans = .true.
         end if

      end subroutine add_fpe_checks
      
      
      
      subroutine multiply_tolerances(id, s, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: status
         
         real(dp) :: test_suite_res_factor = 1
         character(len=20) :: test_suite_resolution_factor_str
         
         include 'formats'
         
         ierr = 0
         call GET_ENVIRONMENT_VARIABLE('MESA_TEST_SUITE_RESOLUTION_FACTOR', &
            test_suite_resolution_factor_str, STATUS=status)
         if (status /= 0) return
         
         if (test_suite_resolution_factor_str .ne. "") then 
            read(test_suite_resolution_factor_str, *) test_suite_res_factor
            write(*,*) ""
            write(*,*) "***"
            write(*,*) "MESA_TEST_SUITE_RESOLUTION_FACTOR set to", test_suite_res_factor
            write(*,*) "***"
            write(*,*) "Warning: This environment variable is for testing purposes"
            write(*,*) "          and should be set to 1 during normal MESA use."
            write(*,*) "***"
            write(*,*) "Multiplying mesh_delta_coeff and time_delta_coeff by this factor,"
            write(*,*) "and max_model_number by its inverse:"
            write(*,*) ""
            write(*,*)    "   old mesh_delta_coeff = ",   s% mesh_delta_coeff
            s% mesh_delta_coeff = test_suite_res_factor * s% mesh_delta_coeff
            write(*,*)    "   new mesh_delta_coeff = ",   s% mesh_delta_coeff
            write(*,*)    ""
            write(*,*)    "   old time_delta_coeff = ",   s% time_delta_coeff
            s% time_delta_coeff = test_suite_res_factor * s% time_delta_coeff
            write(*,*)    "   new time_delta_coeff = ",   s% time_delta_coeff
            write(*,*)    ""
            write(*,*)    "   old max_model_number = ",   s% max_model_number
            s% max_model_number = s% max_model_number / test_suite_res_factor
            write(*,*)    "   new max_model_number = ",   s% max_model_number
            write(*,*)    ""
         end if
      
      end subroutine multiply_tolerances
      
      
      subroutine pgstar_env_check(id, s, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         character(len=5) :: flag
         integer :: status

         include 'formats'

         ierr = 0

         call get_environment_variable('MESA_FORCE_PGSTAR_FLAG', flag, STATUS=status)
         if (status /= 0) return

         select case (trim(flag))
         case ("TRUE", "true")
            write(*,*) "PGSTAR forced on"
            s% job% pgstar_flag = .true.
         case ("FALSE", "false")
            write(*,*) "PGSTAR forced off"
            s% job% pgstar_flag = .false.     
         end select

      end subroutine pgstar_env_check

      end module run_star_support
      
      
      
      
