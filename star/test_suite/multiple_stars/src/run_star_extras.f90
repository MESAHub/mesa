! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use run_star_support
      
      implicit none
      
      integer, parameter :: max_num_stars=10 ! make as large as needed
      integer :: num_stars
      logical :: restart_flag
      character (len=256), dimension(max_num_stars) :: &
         inlist_names, restart_names
      integer :: which_for_pgstar
      real(dp) :: stopping_age
      
      namelist /multi_stars_job/ &
         num_stars, inlist_names, restart_flag, &
         restart_names, which_for_pgstar, stopping_age
         
      integer :: star_ids(max_num_stars)
      
      
      include "test_suite_extras_def.inc"

         
      contains

      include "test_suite_extras.inc"


      subroutine do_run
         
         integer :: id, ierr, i, i_prev, result, result_reason, model_number
         type (star_info), pointer :: s
         character (len=64) :: inlist_fname
         logical :: first_try, continue_evolve_loop
         real(dp) :: sum_times
         real(dp) :: dt
         
         include 'formats'
         
         write(*,*) 'do multi stars'

         ierr = 0
         call test_suite_startup(s, .false., ierr)

         inlist_fname = 'inlist_multi_stars_job'    
         call read_controls(inlist_fname,ierr)
         if (ierr /= 0) return
         
         if (num_stars < 1) then
            write(*,*) 'need to set num_stars >= 1'
            return
         end if
         
         if (num_stars > max_num_stars) then
            write(*,*) 'need to set num_stars <= max_num_stars or rebuild with larger max'
            return
         end if
         
         write(*,*)
         write(*,*)
         
         do i = 1, num_stars
         
            call do_read_star_job(inlist_names(i), ierr)
            if (failed('do_read_star_job')) return
            
            id = id_from_read_star_job ! star allocated by do_read_star_job
            id_from_read_star_job = 0
            
            call star_ptr(id, s, ierr)
            if (failed('star_ptr')) return
         
            call starlib_init(s, ierr) ! okay to do extra calls on this
            if (failed('star_init')) return

            s% inlist_fname = inlist_names(i)

            call star_set_kap_and_eos_handles(id, ierr)
            if (failed('set_star_kap_and_eos_handles')) return
            
            call star_setup(id, inlist_names(i), ierr)
            if (failed('star_setup')) return
            
            star_ids(i) = id
            call extras_controls(s% id, ierr)
            if (failed('extras_controls')) return
            
            call do_star_job_controls_before(id, s, restart_flag, ierr)
            if (failed('do_star_job_controls_before')) return    
                 
            call do_load1_star(id, s, restart_flag, restart_names(i), ierr)
            if (failed('do_load1_star')) return         
            
            call do_star_job_controls_after(id, s, restart_flag, ierr)
            if (failed('do_star_job_controls_after')) return
            
            if (.not. restart_flag) then
               call before_evolve(id, ierr)
               if (failed('before_evolve')) return
            end if
            
            if (i == which_for_pgstar .or. which_for_pgstar < 0) then
               if (.not. restart_flag) then
                  call start_new_run_for_pgstar(s, ierr)
                  if (failed('start_new_run_for_pgstar')) return
               else
                  call show_terminal_header(id, ierr)
                  if (failed('show_terminal_header')) return
                  call restart_run_for_pgstar(s, ierr)
                  if (failed('restart_run_for_pgstar')) return
               end if
            end if
            
            call s% extras_startup(id, restart_flag, ierr)
            if (failed('extras_startup')) return
         
            if (s% job% profile_starting_model) then
               write(*, '(a, i12)') 'save profile for model number ', s% model_number
               call save_profile(id, &
                  3, ierr)
            end if
            
            s% doing_timing = .false.
            
            write(*,*)
            write(*,*)

         end do
         
         continue_evolve_loop = .true.
         i_prev = 0

         evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop
            
            i = select_youngest_star()
            id = star_ids(i)
            call star_ptr(id, s, ierr)
            if (failed('star_ptr')) return
            
            if (s% model_number == s% job% first_model_for_timing) then
               s% doing_timing = .true.
            end if
            
            if (s% job% auto_extend_net) then
               call extend_net(s, ierr)
               if (failed('extend_net')) return
            end if
            
            if (s% center_ye <= s% job% center_ye_limit_for_v_flag .and. .not. s% v_flag) then
               write(*,1) 'have reached center ye limit', &
                  s% center_ye, s% job% center_ye_limit_for_v_flag
               write(*,1) 'set v_flag true'
               call star_set_v_flag(id, .true., ierr)
               if (failed('star_set_v_flag')) return
               if (ierr /= 0) return
            end if
         
            if (s% job% report_mass_not_fe56) call do_report_mass_not_fe56(s)
            if (s% job% report_cell_for_xm > 0) call do_report_cell_for_xm(s)
         
            first_try = .true.
            
            model_number = get_model_number(id, ierr)
            if (failed('get_model_number')) return

            step_loop: do ! may need to repeat this loop
            
               if (stop_now(s, i, id)) then
                  result = terminate
                  result_reason = 0
                  exit step_loop
               end if

               result = star_evolve_step(id, first_try)
               if (result == keep_going) result = check_model(s, id)
               if (result == keep_going) result = star_pick_next_timestep(id)            
               if (result == keep_going) exit step_loop
               
               model_number = get_model_number(id, ierr)
               if (failed('get_model_number')) return
               
               result_reason = get_result_reason(id, ierr)
               if (result == retry) then
                  if (failed('get_result_reason')) return
                  if (s% job% report_retries) &
                     write(*,'(i6,3x,a,/)') model_number, &
                        'retry reason ' // trim(result_reason_str(result_reason))
               end if
               
               if (result == redo) result = star_prepare_to_redo(id)
               if (result == retry) result = star_prepare_to_retry(id)
               if (result == terminate) then
                  if (result_reason == result_reason_normal) then
                     write(*, '(a, i12)') 'save profile for model number ', s% model_number
                     call save_profile(id, &
                        3, ierr)
                  end if
                  continue_evolve_loop = .false.
                  exit step_loop
               end if
               first_try = .false.
               
            end do step_loop
                        
            if (result == keep_going) then
               ! if you have data that needs to be saved and restored for restarts, 
               ! save it in s% extra_iwork and s% extra_work
               ! before calling star_finish_step
               if (s% job% pgstar_flag .and. &
                     (i == which_for_pgstar) .or. (which_for_pgstar < 0)) &
                  call read_pgstar_inlist(s, inlist_names(i),ierr)
               if (failed('read_pgstar_inlist')) return
               result = s% extras_finish_step(id)
               if (result /= keep_going) exit evolve_loop
               result = star_finish_step(id, ierr)
               if (failed('star_finish_step')) return
               if (result /= keep_going) exit evolve_loop
               if (s% job% pgstar_flag .and. (i == which_for_pgstar) .or. (which_for_pgstar < 0)) &
                  call update_pgstar_plots(s, .false., ierr)
               if (failed('update_pgstar_plots')) return
            else if (result == terminate) then
               if (result_reason == result_reason_normal) then
                  result = star_finish_step(id, ierr)
                  if (failed('star_finish_step')) return
                  call do_saves( &
                     id, ierr)
                     if (failed('do_saves evolve_loop')) return
               end if
               exit evolve_loop
            end if
            
            call do_saves( &
               id, ierr)
            if (failed('do_saves')) return
            
         end do evolve_loop

         do i = 1, num_stars
         
            id = star_ids(i)
            
            call star_ptr(id, s, ierr)
            if (failed('star_ptr')) return

            if (s% doing_timing) call show_times(id,s)
         
            result_reason = get_result_reason(id, ierr)
            if (result_reason /= result_reason_normal) then
               write(*, *) 
               write(*, *) 'terminated evolution because ' // trim(result_reason_str(result_reason))
               write(*, *)
            end if

            call s% extras_after_evolve(id, ierr)
            if (failed('after_evolve_extras')) return

            if (s% job% pgstar_flag .and. &
                  (i == which_for_pgstar .or. which_for_pgstar < 0)) &
               call update_pgstar_plots( &
                  s, s% job% save_pgstar_files_when_terminate, &
                  ierr)
            if (failed('update_pgstar_plots')) return

            call star_after_evolve(id, ierr)
            if (failed('star_after_evolve')) return

            call write_terminal_summary(id, ierr)
            if (failed('write_terminal_summary')) return
         
            call free_star(id, ierr)
            if (failed('free_star')) return
            
         end do
         
         call starlib_shutdown
         
         
         call test_suite_after_evolve(s, ierr)


         contains
         

         integer function select_youngest_star()
            integer :: i, i_min
            real(dp) :: age_min
            type (star_info), pointer :: s
            i_min = 0
            age_min = 1d99
            do i = 1, num_stars
               call star_ptr(star_ids(i), s, ierr)
               if (failed('star_ptr')) return
               if (s% star_age < age_min) then
                  age_min = s% star_age
                  i_min = i
               end if
            end do
            select_youngest_star = i_min
            write(*,'(99a20)') 'star', 'model', 'age', 'mass', 'last photo'
            do i = 1, num_stars
               call star_ptr(star_ids(i), s, ierr)
               if (failed('star_ptr')) return
               if (i == i_min) then
                  write(*,'(a18,i2,i20,2(4x,1pe16.9),8x,a)') 'next >', i, &
                     s% model_number, s% star_age, s% star_mass, &
                     trim(s% most_recent_photo_name)
               else
                  write(*,'(2i20,2(4x,1pe16.9),8x,a)') i, &
                     s% model_number, s% star_age, s% star_mass, &
                     trim(s% most_recent_photo_name)
               end if
            end do
            select_youngest_star = i_min
         end function select_youngest_star
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed



      end subroutine do_run    

      
      subroutine read_controls(filename,ierr)
         use utils_lib
         character (len=*) :: filename
         integer, intent(out) :: ierr

         
         character (len=256) :: message
         integer :: unit
         
         ! set defaults
         num_stars = 0
         restart_flag = .false.
         restart_names(:) = 'undefined'
         inlist_names(:) = 'undefined'
         which_for_pgstar = 1
         stopping_age = 1d99
         
         open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=multi_stars_job, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=multi_stars_job)
               close(unit)
            end if  
         end if

      end subroutine read_controls
      
      
      integer function check_model(s, id)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         
         check_model = s% extras_check_model(id)
         if (check_model /= keep_going) return

         check_model = star_check_model(id)
         if (check_model /= keep_going) then
            return
         end if
                              
      end function check_model
      
      
      logical function stop_now(s, i, id)
         type (star_info), pointer :: s
         integer, intent(in) :: i, id
         
         stop_now = (s% star_age > stopping_age)
         if (stop_now) write(*,'(a)') 'all stars have reached stopping age'
      
      end function stop_now
      
      include 'standard_run_star_extras.inc'
      

      end module run_star_extras
      
