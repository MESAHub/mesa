! ***********************************************************************
!
!   Copyright (C) 2020  The MESA Team
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
 
      module astero_run_support

      use star_lib
      use star_def
      use const_def
      use astero_support
      
      implicit none
      
      
      logical, parameter :: scale_simplex_params = .false.


      contains
      
      
      subroutine do_run_star_astero( &
            extras_controls, inlist_astero_search_controls_fname)
         use run_star_support
         use extras_support
         use adipls_support
         use gyre_support, only: gyre_is_enabled, init_gyre
         interface
            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      
         end interface
         character (len=256) :: inlist_astero_search_controls_fname
         optional inlist_astero_search_controls_fname

         type (star_info), pointer :: s
         integer :: id, ierr
         character (len=256) :: inlist_fname
         
         include 'formats'

         ierr = 0
         call do_read_star_job('inlist', ierr) ! this does alloc_star
         ! and saves the id in id_from_read_star_job
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         id = id_from_read_star_job
         id_from_read_star_job = 0
         star_id = id        

         call star_setup(id, 'inlist', ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_setup'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         okay_to_restart = .true.
         
         star_astero_procs% extras_controls => extras_controls

         if (present(inlist_astero_search_controls_fname)) then
            inlist_astero_fname = inlist_astero_search_controls_fname
         else
            inlist_astero_fname = 'inlist_astero_search_controls'
         end if
         write(*,*) 'read ' // trim(inlist_astero_fname)
         call read_astero_search_controls(inlist_astero_fname, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in read_astero_search_controls'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (Y_depends_on_Z .and. vary_Y) then
            vary_Y = .false.
            write(*,*) &
               'WARNING: vary_Y has been changed to false since Y_depends_on_Z is true'
         end if

         nullify(el, order, cyclic_freq, inertia)

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% job% pgstar_flag) then
            call read_astero_pgstar_controls(inlist_astero_fname, ierr)
            if (failed('read_astero_pgstar_controls',ierr)) return
         end if

         
         if (oscillation_code == 'gyre') then
         
            if (gyre_is_enabled) then            
               call init_gyre(gyre_input_file, ierr)
               if (ierr /= 0) return
            ! else give caller a chance to respond before quitting.
            end if
         
         else if (oscillation_code == 'adipls') then 

            call run_adipls(s, .true., .false., &
               add_center_point, keep_surface_point, add_atmosphere, &
               do_redistribute_mesh, ierr)
            if (ierr /= 0) return
            
         else
         
            write(*,'(a)') 'invalid oscillation_code: ' // trim(oscillation_code)
            ierr = -1
            return
         
         end if

         if (save_controls) then
            call write_astero_search_controls(save_controls_filename, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in write_astero_search_controls'
               call mesa_error(__FILE__,__LINE__)
            end if
         end if
         
         call check_search_controls(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in check_search_controls'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         nu_max_sun = s% nu_max_sun
         delta_nu_sun = s% delta_nu_sun
         call init_obs_data(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in init_obs_data'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         next_Y_to_try = -1
         next_FeH_to_try = -1
         next_mass_to_try = -1
         next_alpha_to_try = -1
         next_f_ov_to_try = -1
         next_my_param1_to_try = -1
         next_my_param2_to_try = -1
         next_my_param3_to_try = -1
         sample_number = 0
         max_num_samples = 0
         num_chi2_too_big = 0
         avg_age_top_samples = 1d99
         avg_age_sigma = 1d99
         avg_model_number_top_samples = 1d99
         avg_model_number_sigma = 1d99
         nvar = 0
         total_time_in_oscillation_code = 0d0
         my_var1 = 0d0
         my_var2 = 0d0
         my_var3 = 0d0
         
         call init_sample_ptrs
         
         write(*,*) 'search_type == ' // trim(search_type)
         
         if (search_type == 'use_first_values' .or. &
               s% job% astero_just_call_my_extras_check_model) then
            vary_Y = .false.
            vary_FeH = .false.
            vary_mass = .false.
            vary_alpha = .false.
            vary_f_ov = .false.
            vary_my_param1 = .false.
            vary_my_param2 = .false.
            vary_my_param3 = .false.
            chi2 = eval1(id,ierr)
         else if (search_type == 'simplex') then
            call do_simplex(ierr)
         else if (search_type == 'newuoa') then
            call do_bobyqa_or_newuoa(.true.,ierr)
         else if (search_type == 'bobyqa') then
            call do_bobyqa_or_newuoa(.false.,ierr)
         else if (search_type == 'scan_grid') then
            call do_scan_grid(s, ierr)
         else if (search_type == 'from_file') then
            call do_get_parameters_from_file(s, ierr)
         else 
            write(*,*) 'bad value for search_type ' // trim(search_type)
            ierr = -1
         end if


      end subroutine do_run_star_astero

      
      real(dp) function eval1(id_in,ierr)         
         use run_star_support, only: run1_star
         use extras_support
         
         integer, intent(in) :: id_in
         integer, intent(out) :: ierr
         
         logical, parameter :: &
            do_alloc_star = .false., &
            do_free_star = .false.
            
         type (star_info), pointer :: s
         logical :: restart
         integer :: id, i      

         include 'formats'
         
         ierr = 0
         id = id_in
         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         eval1 = -1
         
         ! init for start of run
         best_chi2 = -1    
         num_chi2_too_big = 0     
         astero_max_dt_next = 1d99
         
         call run1_star( &
            do_alloc_star, do_free_star, okay_to_restart, &
            id, restart, &
            astero_extras_controls, &
            ierr)
         if (ierr /= 0) return
         
         s% max_years_for_timestep = initial_max_years_for_timestep
         s% astero_using_revised_max_yr_dt = .false.
         s% astero_revised_max_yr_dt = s% max_years_for_timestep         
         
         okay_to_restart = .false. ! only allow restart on 1st call to run1_star
         
         eval1 = best_chi2
         
         if (s% job% astero_just_call_my_extras_check_model) return
         
         if (best_chi2 < 0) then
            write(*,*) 'failed to find chi^2 for this run'
            call zero_best_info
            best_chi2 = 999999d0
            return
         end if
         
         sample_number = sample_number + 1
         write(*,*)         
         call show_best(6)
         
         if (write_best_model_data_for_each_sample) &
            call write_best(sample_number)
         
      end function eval1
      
      
      subroutine do_get_parameters_from_file(s, ierr)
         use utils_lib
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer, parameter :: max_col_num = 500
         real(dp) :: filedata(max_col_num)
         integer :: iounit, num_to_read
         
         include 'formats'
         
         
         write(*,*) 'do_get_parameters_from_file'

         sample_number = 0
         ierr = 0

         num_to_read = 0
         
         if (vary_FeH) then
            if (file_column_for_FeH < 1 .or. file_column_for_FeH > max_col_num) then
               write(*,1) 'need to set file_column_for_FeH'
               ierr = -1
               return
            end if
            if (file_column_for_FeH > num_to_read) num_to_read = file_column_for_FeH
         end if
         
         if (vary_Y) then
            if (file_column_for_Y < 1 .or. file_column_for_Y > max_col_num) then
               write(*,1) 'need to set file_column_for_Y'
               ierr = -1
               return
            end if
            if (file_column_for_Y > num_to_read) num_to_read = file_column_for_Y
         end if
         
         if (vary_f_ov) then
            if (file_column_for_f_ov < 1 .or. file_column_for_f_ov > max_col_num) then
               write(*,1) 'need to set file_column_for_f_ov'
               ierr = -1
               return
            end if
            if (file_column_for_f_ov > num_to_read) num_to_read = file_column_for_f_ov
         end if

         if (vary_alpha) then
            if (file_column_for_alpha < 1 .or. file_column_for_alpha > max_col_num) then
               write(*,1) 'need to set file_column_for_alpha'
               ierr = -1
               return
            end if
            if (file_column_for_alpha > num_to_read) num_to_read = file_column_for_alpha
         end if

         if (vary_mass) then
            if (file_column_for_mass < 1 .or. file_column_for_mass > max_col_num) then
               write(*,1) 'need to set file_column_for_mass'
               ierr = -1
               return
            end if
            if (file_column_for_mass > num_to_read) num_to_read = file_column_for_mass
         end if

         if (vary_my_param1) then
            if (file_column_for_my_param1 < 1 .or. file_column_for_my_param1 > max_col_num) then
               write(*,1) 'need to set file_column_for_my_param1'
               ierr = -1
               return
            end if
            if (file_column_for_my_param1 > num_to_read) num_to_read = file_column_for_my_param1
         end if

         if (vary_my_param2) then
            if (file_column_for_my_param2 < 1 .or. file_column_for_my_param2 > max_col_num) then
               write(*,1) 'need to set file_column_for_my_param2'
               ierr = -1
               return
            end if
            if (file_column_for_my_param2 > num_to_read) num_to_read = file_column_for_my_param2
         end if

         if (vary_my_param3) then
            if (file_column_for_my_param3 < 1 .or. file_column_for_my_param3 > max_col_num) then
               write(*,1) 'need to set file_column_for_my_param3'
               ierr = -1
               return
            end if
            if (file_column_for_my_param3 > num_to_read) num_to_read = file_column_for_my_param3
         end if

         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return
         
         open(iounit, file=trim(filename_for_parameters), &
            action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'failed to open filename_for_parameters: ' // &
               trim(filename_for_parameters)
            call free_iounit(iounit)
            return
         end if
         
         write(*,*) 'reading ' // trim(filename_for_parameters)
         write(*,2) 'max_num_from_file', max_num_from_file
         
         read(iounit,*) ! skip 1st line
         
         do while (sample_number < max_num_from_file .or. max_num_from_file < 0)
         
            read(iounit,*,iostat=ierr) filedata(1:num_to_read)
            if (ierr /= 0) then
               write(*,2) 'read failed: sample_number', sample_number
               exit
            end if
            
            if (vary_FeH) then
               next_FeH_to_try = filedata(file_column_for_FeH)
               write(*,1) 'next_FeH_to_try', next_FeH_to_try
            end if
            if (vary_Y) then
               next_Y_to_try = filedata(file_column_for_Y) 
               write(*,1) 'next_Y_to_try', next_Y_to_try
            end if
            if (vary_f_ov) then
               next_f_ov_to_try = filedata(file_column_for_f_ov) 
               write(*,1) 'next_f_ov_to_try', next_f_ov_to_try
            end if
            if (vary_alpha) then
               next_alpha_to_try = filedata(file_column_for_alpha) 
               write(*,1) 'next_alpha_to_try', next_alpha_to_try
            end if
            if (vary_mass) then
               next_mass_to_try = filedata(file_column_for_mass) 
               write(*,1) 'next_mass_to_try', next_mass_to_try
            end if
            if (vary_my_param1) then
               next_my_param1_to_try = filedata(file_column_for_my_param1) 
               write(*,1) 'next_my_param1_to_try', next_my_param1_to_try
            end if
            if (vary_my_param2) then
               next_my_param2_to_try = filedata(file_column_for_my_param2) 
               write(*,1) 'next_my_param2_to_try', next_my_param2_to_try
            end if
            if (vary_my_param3) then
               next_my_param3_to_try = filedata(file_column_for_my_param3) 
               write(*,1) 'next_my_param3_to_try', next_my_param3_to_try
            end if
            write(*,*)
            
            call do1_grid(ierr)
            if (ierr /= 0) then
               write(*,2) 'do1_grid failed: sample_number', sample_number
               exit
            end if
         
         end do
                  
         close(iounit)
         call free_iounit(iounit)


         contains
         
         
         subroutine do1_grid(ierr)
            integer, intent(out) :: ierr
            real(dp) :: chi2
            
            include 'formats'
            ierr = 0
                     
            chi2 = eval1(s% id,ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in eval1'
               return
            end if
            
            call save_best_for_sample(sample_number, 0)

            call save_sample_results_to_file(-1,from_file_output_filename,ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in save_sample_results_to_file'
               return
            end if
            
         end subroutine do1_grid

         
      end subroutine do_get_parameters_from_file
      
      
      subroutine do_scan_grid(s, ierr)
         use utils_lib
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer :: num_FeH, num_Y, num_alpha, num_mass, num_f_ov, &
              num_my_param1, num_my_param2, num_my_param3
         integer :: i_total
         real(dp) :: FeH, Y, alpha, mass, f_ov, chi2, &
              my_param1, my_param2, my_param3
         real(dp), parameter :: eps = 1d-6
         logical :: just_counting
         
         include 'formats'
         
         ierr = 0
         call set_starting_values
         
         num_FeH = 0
         num_Y = 0
         num_alpha = 0
         num_mass = 0
         num_f_ov = 0
         num_my_param1 = 0
         num_my_param2 = 0
         num_my_param3 = 0
         
         just_counting = .true.
         call do_my_param3(ierr)
         i_total = sample_number
         
         write(*,2) 'grid total', i_total
         write(*,2) 'num_FeH', num_FeH
         write(*,2) 'num_Y', num_Y
         write(*,2) 'num_alpha', num_alpha
         write(*,2) 'num_mass', num_mass
         write(*,2) 'num_f_ov', num_f_ov
         write(*,2) 'num_my_param1', num_my_param1
         write(*,2) 'num_my_param2', num_my_param2
         write(*,2) 'num_my_param3', num_my_param3
         write(*,*)
         
         sample_number = 0
         just_counting = .false.
         
         if (restart_scan_grid_from_file) then
            call read_samples_from_file(scan_grid_output_filename, ierr)
            if (ierr /= 0) return
            scan_grid_skip_number = sample_number
            sample_number = 0
            write(*,2) 'scan_grid_skip_number', scan_grid_skip_number
         else
            scan_grid_skip_number = 0
         end if

         call do_my_param3(ierr)
                 
                 
         contains
         
         
         subroutine set_starting_values
            FeH = min_FeH
            Y = min_Y
            alpha = min_alpha
            mass = min_mass
            f_ov = min_f_ov
            my_param1 = min_my_param1
            my_param2 = min_my_param2
            my_param3 = min_my_param3
         end subroutine set_starting_values
         
         
         subroutine do_my_param3(ierr)
            integer, intent(out) :: ierr
            integer :: cnt
            cnt = 0
            ierr = 0
            do while (my_param3 <= max_my_param3 + eps .or. .not. vary_my_param3)          
               if (vary_my_param3) next_my_param3_to_try = my_param3
               call do_my_param2(ierr)
               if (ierr /= 0) return
               cnt = cnt+1
               if (delta_my_param3 <= 0 .or. .not. vary_my_param3) exit
               my_param3 = my_param3 + delta_my_param3         
            end do
            if (num_my_param3 == 0) num_my_param3 = cnt
            my_param3 = min_my_param3
         end subroutine do_my_param3
         
         
         subroutine do_my_param2(ierr)
            integer, intent(out) :: ierr
            integer :: cnt
            cnt = 0
            ierr = 0
            do while (my_param2 <= max_my_param2 + eps .or. .not. vary_my_param2)          
               if (vary_my_param2) next_my_param2_to_try = my_param2
               call do_my_param1(ierr)
               if (ierr /= 0) return
               cnt = cnt+1
               if (delta_my_param2 <= 0 .or. .not. vary_my_param2) exit
               my_param2 = my_param2 + delta_my_param2         
            end do
            if (num_my_param2 == 0) num_my_param2 = cnt
            my_param2 = min_my_param2
         end subroutine do_my_param2


         subroutine do_my_param1(ierr)
            integer, intent(out) :: ierr
            integer :: cnt
            cnt = 0
            ierr = 0
            do while (my_param1 <= max_my_param1 + eps .or. .not. vary_my_param1)          
               if (vary_my_param1) next_my_param1_to_try = my_param1
               call do_f_ov(ierr)
               if (ierr /= 0) return
               cnt = cnt+1
               if (delta_my_param1 <= 0 .or. .not. vary_my_param1) exit
               my_param1 = my_param1 + delta_my_param1         
            end do
            if (num_my_param1 == 0) num_my_param1 = cnt
            my_param1 = min_my_param1
         end subroutine do_my_param1

         
         subroutine do_f_ov(ierr)
            integer, intent(out) :: ierr
            integer :: cnt
            cnt = 0
            ierr = 0
            do while (f_ov <= max_f_ov + eps .or. .not. vary_f_ov)          
               if (vary_f_ov) next_f_ov_to_try = f_ov
               call do_alpha(ierr)
               if (ierr /= 0) return
               cnt = cnt+1
               if (delta_f_ov <= 0 .or. .not. vary_f_ov) exit
               f_ov = f_ov + delta_f_ov         
            end do
            if (num_f_ov == 0) num_f_ov = cnt
            f_ov = min_f_ov
         end subroutine do_f_ov
         
         
         subroutine do_alpha(ierr)
            integer, intent(out) :: ierr
            integer :: cnt
            cnt = 0
            ierr = 0
            do while (alpha <= max_alpha + eps .or. .not. vary_alpha)                 
               if (vary_alpha) next_alpha_to_try = alpha
               call do_FeH(ierr)
               if (ierr /= 0) return
               cnt = cnt+1
               if (delta_alpha <= 0 .or. .not. vary_alpha) exit
               alpha = alpha + delta_alpha                  
            end do
            if (num_alpha == 0) num_alpha = cnt
            alpha = min_alpha
         end subroutine do_alpha
         
         
         subroutine do_FeH(ierr)
            integer, intent(out) :: ierr
            integer :: cnt
            cnt = 0
            ierr = 0
            do while (FeH <= max_FeH + eps .or. .not. vary_FeH)          
               if (vary_FeH) next_FeH_to_try = FeH
               call do_Y(ierr)
               if (ierr /= 0) return
               cnt = cnt+1
               if (delta_FeH <= 0 .or. .not. vary_FeH) exit
               FeH = FeH + delta_FeH            
            end do
            if (num_FeH == 0) num_FeH = cnt
            FeH = min_FeH
         end subroutine do_FeH
         
         
         subroutine do_Y(ierr)
            integer, intent(out) :: ierr
            integer :: cnt
            cnt = 0
            ierr = 0
            do while (Y <= max_Y + eps .or. .not. vary_Y)           
               if (vary_Y) next_Y_to_try = Y  
               call do_mass(ierr)
               if (ierr /= 0) return
               cnt = cnt+1
               if (delta_Y <= 0 .or. .not. vary_Y) exit
               Y = Y + delta_Y               
            end do
            if (num_Y == 0) num_Y = cnt
            Y = min_Y
         end subroutine do_Y
         
         
         subroutine do_mass(ierr)
            integer, intent(out) :: ierr
            integer :: cnt
            cnt = 0
            ierr = 0
            do while (mass <= max_mass + eps .or. .not. vary_mass)             
               if (vary_mass) next_mass_to_try = mass
               call do1_grid(ierr)
               if (ierr /= 0) return
               if (delta_mass <= 0 .or. .not. vary_mass) exit
               cnt = cnt+1
               mass = mass + delta_mass
            end do
            if (num_mass == 0) num_mass = cnt
            mass = min_mass
         end subroutine do_mass
         
         
         subroutine do1_grid(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            
            if (just_counting) then
               sample_number = sample_number + 1
               return
            end if
            
            if (sample_number < scan_grid_skip_number) then
               sample_number = sample_number + 1
               if (mod(sample_number,20) == 1) &
                  write(*,'(67x,99a26)') 'mass', 'Y', 'FeH', 'alpha', 'f_ov', &
                    'my_param1', 'my_param2', 'my_param3'
               write(*,2) 'restore sample from file', sample_number, mass, Y, FeH, alpha, f_ov, &
                    my_param1, my_param2, my_param3
               return
            end if
            
            write(*,2) 'eval1 sample_number', sample_number+1, mass, Y, FeH, alpha, f_ov, &
                 my_param1, my_param2, my_param3
            
            chi2 = eval1(s% id,ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in eval1'
               return
            end if
            
            if (best_chi2 > 9d5) then
               sample_number = sample_number + 1
               write(*,2) 'failed to get chi2 for grid point', sample_number
            else
               write(*,2) 'save best sample for grid point', sample_number
            end if
            
            call save_best_for_sample(sample_number, 0)

            call save_sample_results_to_file(i_total,scan_grid_output_filename,ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in save_sample_results_to_file'
               return
            end if
            
         end subroutine do1_grid

         
      end subroutine do_scan_grid


      subroutine bobyqa_fun(n,x,f)
         integer, intent(in) :: n
         double precision, intent(in) :: x(*)
         double precision, intent(out) :: f
         
         integer :: ierr

         call bobyqa_or_newuoa_fun(n,x,f)

         write(*,*)
         ierr = 0
         call save_sample_results_to_file(-1,bobyqa_output_filename,ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in save_sample_results_to_file'
            call mesa_error(__FILE__,__LINE__,'bobyqa_fun')
         end if
         
      end subroutine bobyqa_fun
      

      subroutine newuoa_fun(n,x,f)
         integer, intent(in) :: n
         double precision, intent(in) :: x(*)
         double precision, intent(out) :: f
         
         integer :: ierr
         
         call bobyqa_or_newuoa_fun(n,x,f)

         write(*,*)
         ierr = 0
         call save_sample_results_to_file(-1,newuoa_output_filename,ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in save_sample_results_to_file'
            call mesa_error(__FILE__,__LINE__,'newuoa_fun')
         end if
         
      end subroutine newuoa_fun


      subroutine bobyqa_or_newuoa_fun(n,x,f)
         integer, intent(in) :: n
         double precision, intent(in) :: x(*)
         double precision, intent(out) :: f
         integer :: ierr, prev_sample_number
         include 'formats'
         
         ierr = 0
         
         if (vary_Y) then
            next_Y_to_try = bobyqa_param( &
               x(i_Y), first_Y, min_Y, max_Y)
            write(*,1) 'next_Y_to_try', next_Y_to_try, x(i_Y)
         end if

         if (vary_FeH) then
            next_FeH_to_try = bobyqa_param( &
               x(i_FeH), first_FeH, min_FeH, max_FeH)
            write(*,1) 'next_FeH_to_try', next_FeH_to_try, x(i_FeH)
         end if
         
         if (vary_mass) then
            next_mass_to_try = bobyqa_param( &
               x(i_mass), first_mass, min_mass, max_mass)
            write(*,1) 'next_mass_to_try', next_mass_to_try
         end if
         
         if (vary_alpha) then
            next_alpha_to_try = bobyqa_param( &
               x(i_alpha), first_alpha, min_alpha, max_alpha)
            write(*,1) 'next_alpha_to_try', next_alpha_to_try, x(i_alpha)
         end if
         
         if (vary_f_ov) then
            next_f_ov_to_try = bobyqa_param( &
               x(i_f_ov), first_f_ov, min_f_ov, max_f_ov)
            write(*,1) 'next_f_ov_to_try', next_f_ov_to_try, x(i_f_ov)
         end if

         if (vary_my_param1) then
            next_my_param1_to_try = bobyqa_param( &
               x(i_my_param1), first_my_param1, min_my_param1, max_my_param1)
            write(*,1) 'next_my_param1_to_try', next_my_param1_to_try, x(i_my_param1)
         end if
         
         if (vary_my_param2) then
            next_my_param2_to_try = bobyqa_param( &
               x(i_my_param2), first_my_param2, min_my_param2, max_my_param2)
            write(*,1) 'next_my_param2_to_try', next_my_param2_to_try, x(i_my_param2)
         end if
         
         if (vary_my_param3) then
            next_my_param3_to_try = bobyqa_param( &
               x(i_my_param3), first_my_param3, min_my_param3, max_my_param3)
            write(*,1) 'next_my_param3_to_try', next_my_param3_to_try, x(i_my_param3)
         end if
         
         prev_sample_number = sample_number
         f = eval1(star_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'got ierr from eval1'
            call mesa_error(__FILE__,__LINE__,'bobyqa_fun')
         end if
         if (sample_number == prev_sample_number) then
            if (sample_number <= 0) then ! failed on 1st try
               write(*,*) 'failed to find chi^2 on 1st try'
               write(*,*) 'must give "first" values that yield a chi^2 result'
               call mesa_error(__FILE__,__LINE__)
            end if
            return ! failed to get new chi^2
         end if
         
         call save_best_for_sample(sample_number, 0)

         write(*,*)
         write(*,*) 'current set of sample results'
         call show_all_sample_results(6,-1,ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in show_all_sample_results'
            call mesa_error(__FILE__,__LINE__,'bobyqa_fun')
         end if
         
         min_sample_chi2_so_far = minval(sample_chi2(1:sample_number))
         
      end subroutine bobyqa_or_newuoa_fun

      
      real(dp) function bobyqa_param(x, first, min, max)
         real(dp), intent(in) :: x, first, min, max
         if (x > 0) then
            bobyqa_param = first + x*(max-first)
         else
            bobyqa_param = first + x*(first-min)
         end if
      end function bobyqa_param
      
      
      subroutine do_bobyqa_or_newuoa(newuoa_flag, ierr)
         use num_lib
         logical, intent(in) :: newuoa_flag
         integer, intent(out) :: ierr
         integer, parameter :: maxfun = 1000, iprint = 0
         real(dp), pointer, dimension(:) :: xl, xu, x, w
         real(dp) :: min_chi2, rhobeg, max_value
         integer :: i, npt
         include 'formats'
         ierr = 0
         
         write(*,*)
         write(*,*)
         
         if (vary_Y) then
            nvar = nvar+1; i_Y = nvar
            if (min_Y >= max_Y) then
               write(*,1) 'min_Y >= max_Y', min_Y, max_Y
               ierr = -1
            end if
         end if
         
         if (vary_FeH) then
            nvar = nvar+1; i_FeH = nvar
            if (min_FeH >= max_FeH) then
               write(*,1) 'min_FeH >= max_FeH', min_FeH, max_FeH
               ierr = -1
            end if
         end if
                  
         if (vary_mass) then
            nvar = nvar+1; i_mass = nvar
            if (min_mass >= max_mass) then
               write(*,1) 'min_mass >= max_mass', &
                  min_mass, max_mass
               ierr = -1
            end if
         end if
         
         if (vary_alpha) then
            nvar = nvar+1; i_alpha = nvar
            if (min_alpha >= max_alpha) then
               write(*,1) 'min_alpha >= max_alpha', &
                  min_alpha, max_alpha
               ierr = -1
            end if
         end if
         
         if (vary_f_ov) then
            nvar = nvar+1; i_f_ov = nvar
            if (min_f_ov >= max_f_ov) then
               write(*,1) 'min_f_ov >= max_f_ov', &
                  min_f_ov, max_f_ov
               ierr = -1
            end if
         end if

         if (vary_my_param1) then
            nvar = nvar+1; i_my_param1 = nvar
            if (min_my_param1 >= max_my_param1) then
               write(*,1) 'min_my_param1 >= max_my_param1', &
                  min_my_param1, max_my_param1
               ierr = -1
            end if
         end if

         if (vary_my_param2) then
            nvar = nvar+1; i_my_param2 = nvar
            if (min_my_param2 >= max_my_param2) then
               write(*,1) 'min_my_param2 >= max_my_param2', &
                  min_my_param2, max_my_param2
               ierr = -1
            end if
         end if

         if (vary_my_param3) then
            nvar = nvar+1; i_my_param3 = nvar
            if (min_my_param3 >= max_my_param3) then
               write(*,1) 'min_my_param3 >= max_my_param3', &
                  min_my_param3, max_my_param3
               ierr = -1
            end if
         end if

         if (ierr /= 0) return
         
         npt = 2*nvar + 1
         
         allocate( &
            xl(nvar), xu(nvar), x(nvar), w((npt+5)*(npt+nvar)+3*nvar*(nvar+5)/2))
         
         XL(1:nvar) = 0
         X(1:nvar) = 0
         XU(1:nvar) = 1
         
!       RHOBEG and bobyqa_rhoend must be set to the initial and final values of a trust
!       region radius, so both must be positive with bobyqa_rhoend no greater than
!       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
!       expected change to a variable, while bobyqa_rhoend should indicate the
!       accuracy that is required in the final values of the variables. An
!       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
!       is less than 2*RHOBEG.
         rhobeg = 0.45d0
         max_value = 1d6
         
         if (newuoa_flag) then
            call newuoa( &
               nvar,npt,x,rhobeg,newuoa_rhoend,iprint,&
               maxfun,w,newuoa_fun,max_value)
         else
            call bobyqa( &
               nvar,npt,x,xl,xu,rhobeg,bobyqa_rhoend,iprint,&
               maxfun,w,bobyqa_fun,max_value)
         end if
         
         if (vary_Y) &
            final_Y = bobyqa_param( &
               x(i_Y), first_Y, min_Y, max_Y)

         if (vary_FeH) &
            final_FeH = bobyqa_param( &
               x(i_FeH), first_FeH, min_FeH, max_FeH)
      
         if (vary_mass) &
            final_mass = bobyqa_param( &
               x(i_mass), first_mass, min_mass, max_mass)
         
         if (vary_alpha) &
            final_alpha = bobyqa_param( &
               x(i_alpha), first_alpha, min_alpha, max_alpha)
         
         if (vary_f_ov) &
            final_f_ov = bobyqa_param( &
               x(i_f_ov), first_f_ov, min_f_ov, max_f_ov)
         
         if (vary_my_param1) &
            final_my_param1 = bobyqa_param( &
               x(i_my_param1), first_my_param1, min_my_param1, max_my_param1)
         
         if (vary_my_param2) &
            final_my_param2 = bobyqa_param( &
               x(i_my_param2), first_my_param2, min_my_param2, max_my_param2)
         
         if (vary_my_param3) &
            final_my_param3 = bobyqa_param( &
               x(i_my_param3), first_my_param3, min_my_param3, max_my_param3)
         
         deallocate(xl, xu, x, w)
      
      
      end subroutine do_bobyqa_or_newuoa


      ! for nelder-mead minimization
      real(dp) function simplex_f( &
            n, x, lrpar, rpar, lipar, ipar, op_code, ierr)
         use const_def, only: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(:) ! (n)
         integer, intent(in) :: lrpar, lipar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: op_code
         integer, intent(out) :: ierr
         
         integer :: prev_sample_number
         include 'formats'
         
         ierr = 0
         
         write(*,*)
         write(*,*)
         
         if (vary_Y) then
            next_Y_to_try = simplex_param( &
               x(i_Y), first_Y, min_Y, max_Y)
            write(*,1) 'next_Y_to_try', next_Y_to_try, x(i_Y)
            if (next_Y_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if

         if (vary_FeH) then
            next_FeH_to_try = simplex_param( &
               x(i_FeH), first_FeH, min_FeH, max_FeH)
            write(*,1) 'next_FeH_to_try', &
               next_FeH_to_try, x(i_FeH)
         end if
         
         if (vary_mass) then
            next_mass_to_try = simplex_param(&
               x(i_mass), first_mass, min_mass, max_mass)
            write(*,1) 'next_mass_to_try', &
               next_mass_to_try, x(i_mass)
            if (next_mass_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         if (vary_alpha) then
            next_alpha_to_try = simplex_param( &
               x(i_alpha), first_alpha, min_alpha, max_alpha)
            write(*,1) 'next_alpha_to_try', &
               next_alpha_to_try, x(i_alpha)
            if (next_alpha_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         if (vary_f_ov) then
            next_f_ov_to_try = simplex_param( &
               x(i_f_ov), first_f_ov, min_f_ov, max_f_ov)
            write(*,1) 'next_f_ov_to_try', &
               next_f_ov_to_try, x(i_f_ov)
            if (next_f_ov_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         if (vary_my_param1) then
            next_my_param1_to_try = simplex_param( &
               x(i_my_param1), first_my_param1, min_my_param1, max_my_param1)
            write(*,1) 'next_my_param1_to_try', &
               next_my_param1_to_try, x(i_my_param1)
            if (next_my_param1_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if

         if (vary_my_param2) then
            next_my_param2_to_try = simplex_param( &
               x(i_my_param2), first_my_param2, min_my_param2, max_my_param2)
            write(*,1) 'next_my_param2_to_try', &
               next_my_param2_to_try, x(i_my_param2)
            if (next_my_param2_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if

         if (vary_my_param3) then
            next_my_param3_to_try = simplex_param( &
               x(i_my_param3), first_my_param3, min_my_param3, max_my_param3)
            write(*,1) 'next_my_param3_to_try', &
               next_my_param3_to_try, x(i_my_param3)
            if (next_my_param3_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         prev_sample_number = sample_number
         simplex_f = eval1(star_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'got ierr from eval1 in astero -- try again'
            ierr = 0
            simplex_f = 1d99
            return
         end if
         
         if (sample_number == prev_sample_number) then
            write(*,*) 'failed to get new chi^2 -- try again'
            simplex_f = 1d99
            return
         end if
         
         call save_best_for_sample(sample_number, op_code)

         call save_sample_results_to_file(-1, simplex_output_filename, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in save_sample_results_to_file'
            call mesa_error(__FILE__,__LINE__)
         end if

         if (simplex_f < simplex_chi2_tol) then
            write(*,*) 'chi2 < simplex_chi2_tol; stopping further iteration'
            ierr = -1
            return
         endif
         
      end function simplex_f

      
      real(dp) function simplex_param(x, first, min, max) result(param)
         real(dp), intent(in) :: x, first, min, max
         if (.not. scale_simplex_params) then
            param = x
         else
            if (x > 0) then
               param = first + x*(max-first)
            else
               param = first + x*(first-min)
            end if
         end if
      end function simplex_param

      
      real(dp) function simplex_inverse(param, first, min, max) result(x)
         real(dp), intent(in) :: param, first, min, max
         if (.not. scale_simplex_params) then
            x = param
         else
            if (param > first) then
               if (max == first) then
                  x = 1d0
               else
                  x = (param-first)/(max-first)
               end if
            else
               if (first == min) then
                  x = -1d0
               else
                  x = (param-first)/(first-min)
               end if
            end if
         end if
      end function simplex_inverse
      
      
      subroutine do_simplex(ierr)
         use num_lib
         integer, intent(out) :: ierr
         
         real(dp) :: final_mass, final_alpha, final_Y, final_FeH, &
              final_my_param1, final_my_param2, final_my_param3
         real(dp), dimension(:), pointer :: x_first, x_lower, x_upper, x_final
         real(dp), pointer :: simplex(:,:), f(:)
         real(dp) :: f_final
         integer :: lrpar, lipar
         integer, pointer :: ipar(:) ! (lipar)
         real(dp), pointer :: rpar(:) ! (lrpar)
         integer :: num_iters, num_fcn_calls, &
            num_fcn_calls_for_ars, num_accepted_for_ars
         integer :: seed, i, j, k, num_samples
         logical :: start_from_given_simplex_and_f
         
         include 'formats'
         
         ierr = 0
         
         if (vary_mass) then
            nvar = nvar+1; i_mass = nvar
            if (min_mass >= max_mass) then
               write(*,1) 'min_mass >= max_mass', &
                  min_mass, max_mass
               ierr = -1
            end if
         end if

         if (vary_Y) then
            nvar = nvar+1; i_Y = nvar
            if (min_Y >= max_Y) then
               write(*,1) 'min_Y >= max_Y', min_Y, max_Y
               ierr = -1
            end if
         end if
         
         if (vary_FeH) then
            nvar = nvar+1; i_FeH = nvar
            if (min_FeH >= max_FeH) then
               write(*,1) 'min_FeH >= max_FeH', &
                  min_FeH, max_FeH
               ierr = -1
            end if
         end if
         
         if (vary_alpha) then
            nvar = nvar+1; i_alpha = nvar
            if (min_alpha >= max_alpha) then
               write(*,1) 'min_alpha >= max_alpha', &
                  min_alpha, max_alpha
               ierr = -1
            end if
         end if
         
         if (vary_f_ov) then
            nvar = nvar+1; i_f_ov = nvar
            if (min_f_ov >= max_f_ov) then
               write(*,1) 'min_f_ov >= max_f_ov', &
                  min_f_ov, max_f_ov
               ierr = -1
            end if
         end if

         if (vary_my_param1) then
            nvar = nvar+1; i_my_param1 = nvar
            if (min_my_param1 >= max_my_param1) then
               write(*,1) 'min_my_param1 >= max_my_param1', &
                  min_my_param1, max_my_param1
               ierr = -1
            end if
         end if

         if (vary_my_param2) then
            nvar = nvar+1; i_my_param2 = nvar
            if (min_my_param2 >= max_my_param2) then
               write(*,1) 'min_my_param2 >= max_my_param2', &
                  min_my_param2, max_my_param2
               ierr = -1
            end if
         end if

         if (vary_my_param3) then
            nvar = nvar+1; i_my_param3 = nvar
            if (min_my_param3 >= max_my_param3) then
               write(*,1) 'min_my_param3 >= max_my_param3', &
                  min_my_param3, max_my_param3
               ierr = -1
            end if
         end if

         if (ierr /= 0) return
         
         lrpar = 0; lipar = 0

         allocate( &
            rpar(lrpar), ipar(lipar), simplex(nvar,nvar+1), f(nvar+1), &
            x_lower(nvar), x_upper(nvar), x_first(nvar), x_final(nvar))
         
         if (.not. scale_simplex_params) then
            call set_xs
         else ! values are scaled to -1..1 with first at 0
            x_lower(1:nvar) = -1
            x_upper(1:nvar) = 1
            x_first(1:nvar) = 0
         end if
                  
         if (restart_simplex_from_file) then
            call read_samples_from_file(simplex_output_filename, ierr)
            if (ierr /= 0) return
            if (sample_number < nvar+1) then
               write(*,2) 'sorry: too few points. for simplex restart need at least', nvar+1
               ierr = -1
               return
            end if
            num_samples = sample_number
            call setup_simplex_and_f(ierr)
            if (ierr /= 0) return            
            start_from_given_simplex_and_f = .true.
            call set_sample_averages
         else
            start_from_given_simplex_and_f = .false.
         end if
         
         call NM_simplex( &
            nvar, x_lower, x_upper, x_first, x_final, f_final, &
            simplex, f, start_from_given_simplex_and_f, simplex_f, &
            simplex_x_atol, simplex_x_rtol, &
            simplex_itermax, simplex_fcn_calls_max, &
            simplex_centroid_weight_power, simplex_enforce_bounds, &
            simplex_adaptive_random_search, simplex_seed, &
            simplex_alpha, simplex_beta, &
            simplex_gamma, simplex_delta, &
            lrpar, rpar, lipar, ipar, &
            num_iters, num_fcn_calls, &
            num_fcn_calls_for_ars, num_accepted_for_ars, ierr)
         
         if (vary_Y) &
            final_Y = simplex_param( &
               x_final(i_Y), first_Y, min_Y, max_Y)

         if (vary_FeH) &
            final_FeH = simplex_param( &
               x_final(i_FeH), first_FeH, &
               min_FeH, max_FeH)
         
         if (vary_mass) &
            final_mass = simplex_param( &
               x_final(i_mass), first_mass, &
               min_mass, max_mass)
         
         if (vary_alpha) &
            final_alpha = simplex_param( &
               x_final(i_alpha), first_alpha, &
               min_alpha, max_alpha)
         
         if (vary_f_ov) &
            final_f_ov = simplex_param( &
               x_final(i_f_ov), first_f_ov, &
               min_f_ov, max_f_ov)
         
         if (vary_my_param1) &
            final_my_param1 = simplex_param( &
               x_final(i_my_param1), first_my_param1, &
               min_my_param1, max_my_param1)
         
         if (vary_my_param2) &
            final_my_param2 = simplex_param( &
               x_final(i_my_param2), first_my_param2, &
               min_my_param2, max_my_param2)
         
         if (vary_my_param3) &
            final_my_param3 = simplex_param( &
               x_final(i_my_param3), first_my_param3, &
               min_my_param3, max_my_param3)

         deallocate( &
            rpar, ipar, simplex, f, x_lower, x_upper, x_first, x_final)
            
            
         contains
         
         
         subroutine set_xs ! x_first, x_lower, x_upper
            if (vary_Y) then
               x_first(i_Y) = first_Y
               x_lower(i_Y) = min_Y
               x_upper(i_Y) = max_Y
            end if
            if (vary_FeH) then
               x_first(i_FeH) = first_FeH
               x_lower(i_FeH) = min_FeH
               x_upper(i_FeH) = max_FeH
            end if
            if (vary_mass) then
               x_first(i_mass) = first_mass
               x_lower(i_mass) = min_mass
               x_upper(i_mass) = max_mass
            end if
            if (vary_alpha) then
               x_first(i_alpha) = first_alpha
               x_lower(i_alpha) = min_alpha
               x_upper(i_alpha) = max_alpha
            end if
            if (vary_f_ov) then
               x_first(i_f_ov) = first_f_ov
               x_lower(i_f_ov) = min_f_ov
               x_upper(i_f_ov) = max_f_ov
            end if         
            if (vary_my_param1) then
               x_first(i_my_param1) = first_my_param1
               x_lower(i_my_param1) = min_my_param1
               x_upper(i_my_param1) = max_my_param1
            end if         
            if (vary_my_param2) then
               x_first(i_my_param2) = first_my_param2
               x_lower(i_my_param2) = min_my_param2
               x_upper(i_my_param2) = max_my_param2
            end if         
            if (vary_my_param3) then
               x_first(i_my_param3) = first_my_param3
               x_lower(i_my_param3) = min_my_param3
               x_upper(i_my_param3) = max_my_param3
            end if         
         end subroutine set_xs
         
         
         subroutine setup_simplex_and_f(ierr)
            use num_lib, only: qsort
            integer, intent(out) :: ierr
            
            integer :: j, i, k, max_i, jj
            integer, pointer :: index(:)
            ! sort results by increasing sample_chi2
            
            include 'formats'
            
            ierr = 0
            allocate(index(num_samples), stat=ierr)
            if (ierr /= 0) then
               call mesa_error(__FILE__,__LINE__,'failed in allocate before calling qsort from show_all_sample_results')
            end if
            call qsort(index, num_samples, sample_chi2)
            max_i = 0
            do j=1,nvar+1
               i = index(j)
               if (i > max_i) max_i = i ! max sample restored
               write(*,3) 'restore simplex', j, i
               f(j) = sample_chi2(i)
               write(*,3) 'chi2', j, i, f(j)
               if (vary_Y) then
                  simplex(i_Y,j) = &
                     simplex_inverse(sample_init_Y(i), first_Y, min_Y, max_Y)
                  write(*,3) 'Y', j, i, sample_init_Y(i)
               end if
               if (vary_FeH) then
                  simplex(i_FeH,j) = &
                     simplex_inverse(sample_init_FeH(i), first_FeH, min_FeH, max_FeH)
                  write(*,3) 'FeH', j, i, sample_init_FeH(i)
               end if         
               if (vary_mass) then
                  simplex(i_mass,j) = &
                     simplex_inverse(sample_mass(i), first_mass, min_mass, max_mass)
                  write(*,3) 'mass', j, i, sample_mass(i)
               end if        
               if (vary_alpha) then
                  simplex(i_alpha,j) = &
                     simplex_inverse(sample_alpha(i), first_alpha, min_alpha, max_alpha)
                  write(*,3) 'alpha', j, i, sample_alpha(i)
               end if         
               if (vary_f_ov) then
                  simplex(i_f_ov,j) = &
                     simplex_inverse(sample_f_ov(i), first_f_ov, min_f_ov, max_f_ov)
                  write(*,3) 'f_ov', j, i, sample_f_ov(i)
               end if
               if (vary_my_param1) then
                  simplex(i_my_param1,j) = &
                     simplex_inverse(sample_my_param1(i), first_my_param1, min_my_param1, max_my_param1)
                  write(*,3) 'my_param1', j, i, sample_my_param1(i)
               end if
               if (vary_my_param2) then
                  simplex(i_my_param2,j) = &
                     simplex_inverse(sample_my_param2(i), first_my_param2, min_my_param2, max_my_param2)
                  write(*,3) 'my_param2', j, i, sample_my_param2(i)
               end if
               if (vary_my_param3) then
                  simplex(i_my_param3,j) = &
                     simplex_inverse(sample_my_param3(i), first_my_param3, min_my_param3, max_my_param3)
                  write(*,3) 'my_param3', j, i, sample_my_param3(i)
               end if
               write(*,*)
            end do
            
            deallocate(index)

            write(*,2) 'num_samples', max_i
            j = 0
            do i=max_i+1,num_samples
               j = j+1
               do jj=1,j
                  write(*,'(a)',advance='no') ' damn, '
               end do
               write(*,'(a)',advance='no') 'will have to rerun sample '
               if (i < 10) then
                  write(*,'(i1)') i
               else if (i < 100) then
                  write(*,'(i2)') i
               else if (i < 1000) then
                  write(*,'(i3)') i
               else if (i < 10000) then
                  write(*,'(i4)') i
               else if (i < 100000) then
                  write(*,'(i5)') i
               else
                  write(*,'(i6)') i
               end if
            end do
            write(*,*)
            write(*,*)
            num_samples = max_i
            
         end subroutine setup_simplex_and_f
         
      
      end subroutine do_simplex
      
      
      subroutine save_best_for_sample(i, op_code)
         integer, intent(in) :: i, op_code
         integer :: ierr
                
         if (i <= 0) return
         if (i > max_num_samples) then
            call alloc_sample_ptrs(ierr)
            if (ierr /= 0) then
               write(*,*) 'ERROR -- failed to allocate for samples'
               call mesa_error(__FILE__,__LINE__,'save_best_for_sample')
               return
            end if
         end if
                  
         sample_op_code(i) = op_code
         sample_chi2(i) = best_chi2
         sample_chi2_seismo(i) = best_chi2_seismo
         sample_chi2_spectro(i) = best_chi2_spectro
         
         sample_age(i) = best_age
         sample_init_Y(i) = current_Y
         sample_init_FeH(i) = current_FeH
         sample_mass(i) = current_mass
         sample_alpha(i) = current_alpha
         sample_f_ov(i) = current_f_ov
         sample_my_param1(i) = current_my_param1
         sample_my_param2(i) = current_my_param2
         sample_my_param3(i) = current_my_param3

         sample_init_h1(i) = current_h1
         sample_init_he3(i) = current_he3
         sample_init_he4(i) = current_he4
         sample_init_Z(i) = current_Z

         sample_radius(i) = best_radius
         sample_logL(i) = best_logL
         sample_Teff(i) = best_Teff
         sample_logg(i) = best_logg
         sample_FeH(i) = best_FeH
         
         sample_logR(i) = best_logR
         sample_surface_Z_div_X(i) = best_surface_Z_div_X
         sample_surface_He(i) = best_surface_He
         sample_Rcz(i) = best_Rcz
         sample_my_var1(i) = best_my_var1
         sample_my_var2(i) = best_my_var2
         sample_my_var3(i) = best_my_var3
         
         sample_surf_coef1(i) = best_surf_coef1
         sample_surf_coef2(i) = best_surf_coef2
         sample_delta_nu(i) = best_delta_nu
         sample_nu_max(i) = best_nu_max
         sample_model_number(i) = best_model_number

         sample_order(0,:,i) = best_order(0,:)
         sample_freq(0,:,i) = best_freq(0,:)
         sample_freq_corr(0,:,i) = best_freq_corr(0,:)
         sample_inertia(0,:,i) = best_inertia(0,:)
      
         sample_order(1,:,i) = best_order(1,:)
         sample_freq(1,:,i) = best_freq(1,:)
         sample_freq_corr(1,:,i) = best_freq_corr(1,:)
         sample_inertia(1,:,i) = best_inertia(1,:)
      
         sample_order(2,:,i) = best_order(2,:)
         sample_freq(2,:,i) = best_freq(2,:)
         sample_freq_corr(2,:,i) = best_freq_corr(2,:)
         sample_inertia(2,:,i) = best_inertia(2,:)
      
         sample_order(3,:,i) = best_order(3,:)
         sample_freq(3,:,i) = best_freq(3,:)
         sample_freq_corr(3,:,i) = best_freq_corr(3,:)
         sample_inertia(3,:,i) = best_inertia(3,:)
         
         sample_ratios_r01(:,i) = best_ratios_r01(:)
         sample_ratios_r10(:,i) = best_ratios_r10(:)
         sample_ratios_r02(:,i) = best_ratios_r02(:)

         call set_sample_averages

      end subroutine save_best_for_sample
      
      
      subroutine set_sample_averages
         integer :: ierr, jj, j, n
         real(dp) :: avg_age_top_samples2, avg_model_number_top_samples2
         
         include 'formats'
      
         call set_sample_index_by_chi2
         n = min(sample_number, max_num_samples_for_avg)
         if (n < max(2,min_num_samples_for_avg)) then
            avg_age_top_samples = 1d99
            avg_age_sigma = 1d99
            avg_model_number_top_samples = 1d99
            avg_model_number_sigma = 1d99
            return
         end if
         avg_age_top_samples = 0
         avg_model_number_top_samples = 0
         avg_age_top_samples2 = 0
         avg_model_number_top_samples2 = 0
         do jj=1,n
            j = sample_index_by_chi2(jj)
            !write(*,3) 'j, jj', j, jj
            avg_age_top_samples = &
               avg_age_top_samples + sample_age(j)
            avg_age_top_samples2 = &
               avg_age_top_samples2 + sample_age(j)*sample_age(j)
            avg_model_number_top_samples = &
               avg_model_number_top_samples + sample_model_number(j)
            avg_model_number_top_samples2 = &
               avg_model_number_top_samples2 + &
               sample_model_number(j)*sample_model_number(j)
            !write(*,2) 'avg_age_top_samples', j, avg_age_top_samples/j
            !write(*,2) 'avg_age_top_samples2', j, avg_age_top_samples2/j
            !write(*,2) 'avg_model_number_top_samples', j, avg_model_number_top_samples/j
            !write(*,2) 'avg_model_number_top_samples2', j, avg_model_number_top_samples2/j
         end do
         avg_age_sigma = &
            sqrt(max(0d0,(avg_age_top_samples2 - avg_age_top_samples*avg_age_top_samples/n)/(n-1)))
         avg_age_top_samples = avg_age_top_samples/n
         avg_model_number_sigma = &
            sqrt(max(0d0,(avg_model_number_top_samples2 - &
                  avg_model_number_top_samples*avg_model_number_top_samples/n)/(n-1)))
         avg_model_number_top_samples = avg_model_number_top_samples/n
         
         write(*,*)
         write(*,2) 'n for averages', n
         write(*,1) 'avg_age_top_samples', avg_age_top_samples
         write(*,1) 'avg_age_sigma', avg_age_sigma
         write(*,1) 'age limit', avg_age_top_samples + avg_age_sigma_limit*avg_age_sigma
         write(*,1) 'avg_model_number_top_samples', avg_model_number_top_samples
         write(*,1) 'avg_model_number_sigma', avg_model_number_sigma
         write(*,1) 'model number limit', &
            avg_model_number_top_samples + &
               avg_model_number_sigma_limit*avg_model_number_sigma
         write(*,*)
         !call mesa_error(__FILE__,__LINE__,'set_sample_averages')
         
      end subroutine set_sample_averages
      
      
      subroutine zero_best_info
         best_chi2 = 0
         best_chi2_seismo = 0
         best_chi2_spectro = 0
         best_init_h1 = 0
         best_init_he3 = 0
         best_init_he4 = 0
         best_init_Z = 0
         best_age = 0
         best_radius = 0
         best_logL = 0
         best_Teff = 0
         best_logg = 0
         best_FeH = 0
         best_logR = 0
         best_surface_Z_div_X = 0
         best_surface_He = 0
         best_Rcz = 0
         best_my_var1 = 0
         best_my_var2 = 0
         best_my_var3 = 0
         best_my_param1 = 0
         best_my_param2 = 0
         best_my_param3 = 0
         best_delta_nu = 0
         best_nu_max = 0
         best_surf_coef1 = 0
         best_surf_coef2 = 0
         best_model_number = 0
         best_order = 0
         best_freq = 0
         best_freq_corr = 0
         best_inertia = 0
         best_ratios_r01 = 0
         best_ratios_r10 = 0
         best_ratios_r02 = 0
      end subroutine zero_best_info


      end module astero_run_support
