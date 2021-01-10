! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton, Aaron Dotter
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

      
      include "test_suite_extras_def.inc"
      
      logical :: start_from_photo, save_last_model
      character (len=256) :: masses_filename, pre_zahb_model, pre_zahb_photo
      
      namelist /create_zahb/ &
         masses_filename, pre_zahb_model, pre_zahb_photo, start_from_photo, save_last_model
         
      real(dp), pointer :: masses(:)
      integer :: nmasses

      contains

      include "test_suite_extras.inc"
      
      
      subroutine zahb_extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
         
         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         s% extras_check_model => isochrone_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => no_extra_history_columns
         s% data_for_extra_history_columns => none_for_extra_history_columns
         s% how_many_extra_profile_columns => no_extra_profile_columns
         s% data_for_extra_profile_columns => none_for_extra_profile_columns  
         
      end subroutine zahb_extras_controls


      subroutine do_run
         
         integer :: id, ierr
         logical :: restart
         type (star_info), pointer :: s
         real(dp) :: dt
         
         write(*,*) 'do create_zahb'

         ierr = 0
         call test_suite_startup(s, restart, ierr)
         
         call do_read_star_job('inlist', ierr)
         if (failed('do_read_star_job',ierr)) return

         id = id_from_read_star_job
         id_from_read_star_job = 0
         
         call star_ptr(id, s, ierr)
         if (failed('star_ptr',ierr)) return
      
         call starlib_init(s, ierr)
         if (failed('star_init',ierr)) return

         s% inlist_fname = 'inlist'
         
         call star_set_kap_and_eos_handles(id, ierr)
         if (failed('set_star_kap_and_eos_handles',ierr)) return
                  
         call star_setup(id, 'inlist', ierr)
         if (failed('star_setup',ierr)) return

         call zahb_extras_controls(s% id, ierr)
         if (failed('extras_controls',ierr)) return

         call do_star_job_controls_before(id, s, restart, ierr)
         if (failed('do_star_job_controls_before',ierr)) return

         call do_isochrone( s, id, ierr)
         
         call test_suite_after_evolve(s, ierr)


      end subroutine do_run    

      
      subroutine read_controls(filename,ierr)
         use utils_lib
         character (len=*) :: filename
         integer, intent(out) :: ierr

         
         character (len=256) :: message
         integer :: unit
         
         ! set defaults
         start_from_photo = .false.
         pre_zahb_photo = 'pre_zahb.photo'
         pre_zahb_model = 'pre_zahb.mod'
         masses_filename = 'zahb_masses.list'
         save_last_model = .false.
         
         open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=create_zahb, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=create_zahb)
               close(unit)
            end if  
         end if

      end subroutine read_controls
      

      subroutine do_isochrone( s, id, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         call do_multi_track( &
            s, id, & !isochrone_inlist, history_columns_file, profile_columns_file, &
!             isochrone_check_model, & !create_pre_main_sequence_model, &
            !pre_ms_T_c, pre_ms_guess_rho_c, pre_ms_d_log10_P, &
            ierr)
      end subroutine do_isochrone
      

      subroutine do_multi_track( s, id, ierr)
         use run_star_support
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         
         integer :: i, j, result
         logical :: first_try
         character (len=256) :: history_name, profiles_index_name, &
            profile_data_prefix, zahb_inlist, last_model

         real(dp) :: next_mass

         logical, parameter :: dbg = .true.
         
         include 'formats'
         
         ierr = 0
         
         zahb_inlist = 'inlist_create_zahb'
         call read_controls(zahb_inlist,ierr)
         if (ierr /= 0) return

         nullify(masses)
         call read_masses(trim(masses_filename), masses, nmasses, ierr)
         if (ierr /= 0) return
      
         if (dbg) then
            write(*,*)
            write(*,2) 'done read_masses', nmasses
            write(*,*)
            do i=1, nmasses
               write(*,2) 'mass', i, masses(i)
            end do
            write(*,*)
            !stop 'do_isochrone'
         end if
         
         do i=1, nmasses
            
            next_mass = masses(i)
            if (next_mass <= 0) cycle ! skip this one
            
            if (i < 10) then
               write(history_name, '(a,i1,a)') 'i00', i, '_history.data'
               write(profiles_index_name, '(a,i1,a)') 'i00', i, '_profiles.index'
               write(profile_data_prefix, '(a,i1,a)') 'i00', i, '_log'
            else if (i < 100) then
               write(history_name, '(a,i2,a)') 'i0', i, '_history.data'
               write(profiles_index_name, '(a,i2,a)') 'i0', i, '_profiles.index'
               write(profile_data_prefix, '(a,i2,a)') 'i0', i, '_log'
            else
               write(history_name, '(a,i3,a)') 'i', i, '_history.data'
               write(profiles_index_name, '(a,i3,a)') 'i', i, '_profiles.index'
               write(profile_data_prefix, '(a,i3,a)') 'i', i, '_log'
            end if 
            s% star_history_name = trim(history_name)
            s% profiles_index_name = trim(profiles_index_name)
            s% profile_data_prefix = trim(profile_data_prefix)
           
            if(start_from_photo)then
               write(*,2) 'reading photo'
               call star_load_restart_photo(id, trim(pre_zahb_photo), ierr)
               if(failed('star_load_restart_photo',ierr)) return
            else 
               write(*,2) 'reading saved model'
               call star_read_model(id, trim(pre_zahb_model), ierr)
               if(failed('star_read_model',ierr)) return
            endif

            call do_star_job_controls_after(id, s, .false., ierr)
            if (failed('do_star_job_controls_after',ierr)) return

            write(*,1) 'relax_mass', next_mass
            s% job% lg_max_abs_mdot=-5d0
            call star_relax_mass(id, next_mass, s% job% lg_max_abs_mdot, ierr)
            if(failed('star_relax_mass',ierr)) return

         
            if (len_trim(s% job% history_columns_file) > 0) &
               write(*,*) 'read ' // trim(s% job% history_columns_file)
            call star_set_history_columns(id, s% job% history_columns_file, .true., ierr)
            if (failed('star_set_history_columns',ierr)) return
         
            if (len_trim(s% job% profile_columns_file) > 0) &
               write(*,*) 'read ' // trim(s% job% profile_columns_file)
            call star_set_profile_columns(id, s% job% profile_columns_file, .true., ierr)
            if (failed('star_set_profile_columns',ierr)) return

            evolve_loop: do ! evolve one step per loop
               first_try = .true.
               step_loop: do ! repeat for retry
                  result = star_evolve_step(id, first_try)
                  if (result == keep_going) result = s% extras_check_model(id)
                  if (result == keep_going) result = star_pick_next_timestep(id)            
                  if (result == keep_going) exit step_loop
                  if (result == redo) result = star_prepare_to_redo(id)
                  if (result == retry) result = star_prepare_to_retry(id)
                  if (result == terminate) exit evolve_loop
                  first_try = .false.
               end do step_loop
               result = star_finish_step(id, ierr)
               if (result /= keep_going) exit evolve_loop         
            end do evolve_loop

            call save_profile(id, &
                  5, ierr)
            if (failed('star_after_evolve',ierr)) return

            if(save_last_model)then
            ! set filename for isochrone mass number
               if (i < 10) then
                  write(last_model, '(a,i1,a)') 'i00', i, '_star.mod'
               else if (i < 100) then
                  write(last_model, '(a,i2,a)') 'i0', i, '_star.mod'
               else
                  write(last_model, '(a,i3,a)') 'i', i, '_star.mod'
               end if 
               call star_write_model(id, last_model, ierr)
               write(*,'(a25,a13)') '      save last model to ', trim(last_model)
               if(failed('star_write_model',ierr)) return
            endif

            call show_terminal_header(id, ierr)
            if (failed('show_terminal_header',ierr)) return

            call write_terminal_summary(id, ierr)
            if (failed('write_terminal_summary',ierr)) return
            
            call star_after_evolve(id, ierr)
            if (failed('star_after_evolve',ierr)) return
            
            do j = 1, 10
               write(*,*)
            end do
                        
         end do
         
         write(*,*)

         write (*, *)
         write (*, *)
         write(*, *)
         write(*, '(a)') 'finished'
         write(*, *)
         

      end subroutine do_multi_track
      

      integer function isochrone_check_model(id)
         integer, intent(in) :: id
         include 'formats'
         isochrone_check_model = star_check_model(id)
      end function isochrone_check_model

      
      integer function no_extra_history_columns(id)
         integer, intent(in) :: id
         no_extra_history_columns = 0
      end function no_extra_history_columns
      
      
      subroutine none_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine none_for_extra_history_columns


      integer function no_extra_profile_columns(id)
         integer, intent(in) :: id
         no_extra_profile_columns = 0
      end function no_extra_profile_columns
      
      
      subroutine none_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine none_for_extra_profile_columns


      include 'standard_run_star_extras.inc'

      end module run_star_extras
