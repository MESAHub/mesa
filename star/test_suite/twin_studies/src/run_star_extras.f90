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
      
      
      include 'test_suite_extras_def.inc'
      
      contains

      include 'test_suite_extras.inc'
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         real(dp) :: target_period, rel_run_E_err
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         if (s% x_integer_ctrl(1) <= 0) return
         if (s% rsp_num_periods < s% x_integer_ctrl(1)) return
         write(*,*)
         write(*,*)
         write(*,*)
         target_period = s% x_ctrl(1)
         rel_run_E_err = s% cumulative_energy_error/s% total_energy
         write(*,*) 'rel_run_E_err', rel_run_E_err
         if (s% total_energy /= 0d0 .and. abs(rel_run_E_err) > 1d-5) then
            write(*,*) '*** BAD rel_run_E_error ***', &
            s% cumulative_energy_error/s% total_energy
         else if (abs(s% rsp_period/(24*3600) - target_period) > 1d-2) then
            write(*,*) '*** BAD ***', s% rsp_period/(24*3600) - target_period, &
               s% rsp_period/(24*3600), target_period
         else
            write(*,*) 'good match for period', &
               s% rsp_period/(24*3600), target_period
         end if
         write(*,*)
         write(*,*)
         write(*,*)
         extras_finish_step = terminate
      end function extras_finish_step
      
      
      subroutine do_run
         use utils_lib, only: cp
         integer :: id, ierr, i, i_prev, result, result_reason, model_number
         type (star_info), pointer :: s
         character (len=64) :: inlist_fname
         logical :: first_try, continue_evolve_loop
         real(dp) :: sum_times
         real(dp) :: dt
         logical :: restart, pgstar_ok
         logical, parameter :: &
            do_alloc_star = .true., &
            do_free_star = .true., &
            okay_to_restart = .true., &
            dbg = .false.
         
         include 'formats'

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
         
         do i = 1, num_stars
            
            id_from_read_star_job = 0
            call do_read_star_job_and_return_id(inlist_names(i), id, ierr)
            if (failed('do_read_star_job',ierr)) return
            star_ids(i) = id
            
            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return
            
            pgstar_ok = (i == which_for_pgstar .or. which_for_pgstar < 0)

            call start_run1_star( &
               do_alloc_star, do_free_star, okay_to_restart, &
               id, restart, pgstar_ok, dbg, &
               extras_controls, &
               ierr, inlist_names(i))
            if (failed('before_evolve_loop',ierr)) return

            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return
            
            s% doing_timing = .false.

         end do
         
         write(*,*) 'finished startup for stars'
         write(*,*)
         
         continue_evolve_loop = .true.
         i_prev = 0

         evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop
            
            i = select_youngest_star()
            if (i == 0) then
               stop 'failed to find youngest'
            end if
            call star_ptr(star_ids(i), s, ierr)
            if (failed('star_ptr',ierr)) return
         
            continue_evolve_loop = do_evolve_one_step(s, dbg, ierr)
            if (failed('do_evolve_one_step',ierr)) return
            
         end do evolve_loop

         do i = 1, num_stars
            
            call star_ptr(star_ids(i), s, ierr)
            if (failed('star_ptr',ierr)) return

            call after_evolve_loop(s% id, do_free_star, ierr)
            if (failed('after_evolve_loop',ierr)) return
            
         end do
         
         call starlib_shutdown
         
         
         call test_suite_after_evolve(s, ierr)


         contains
         

         integer function select_youngest_star()
            integer :: i
            real(dp) :: age_min
            type (star_info), pointer :: s
            select_youngest_star = 0
            age_min = 1d99
            do i = 1, num_stars
               call star_ptr(star_ids(i), s, ierr)
               if (failed('star_ptr',ierr)) return
               if (s% star_age < age_min) then
                  age_min = s% star_age
                  select_youngest_star = i
               end if
            end do
         end function select_youngest_star

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


      end subroutine do_run    


      end module run_star_extras
      
