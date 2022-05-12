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
         integer :: id, i, ierr
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

         num_constraints = 0
         do i = 1, max_constraints
            if (constraint_name(i) /= '') num_constraints = num_constraints + 1
         end do
         
         num_parameters = 0
         do i = 1, max_parameters
            if (param_name(i) /= '') num_parameters = num_parameters + 1
         end do

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
         
         next_param_to_try(1:max_parameters) = -1
         sample_number = 0
         max_num_samples = 0
         num_chi2_too_big = 0
         avg_age_top_samples = 1d99
         avg_age_sigma = 1d99
         avg_model_number_top_samples = 1d99
         avg_model_number_sigma = 1d99
         nvar = 0
         total_time_in_oscillation_code = 0d0
         constraint_value = 0d0
         
         call init_sample_ptrs
         
         write(*,*) 'search_type == ' // trim(search_type)
         
         if (search_type == 'use_first_values' .or. &
               s% job% astero_just_call_my_extras_check_model) then
            vary_param(1:max_parameters) = .false.
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
         integer :: iounit, num_to_read, i
         
         include 'formats'
         
         
         write(*,*) 'do_get_parameters_from_file'

         sample_number = 0
         ierr = 0

         num_to_read = 0

         do i = 1, max_parameters
            if (param_name(i) /= '' .and. vary_param(i)) then
               if (file_column_for_param(i) < 1 .or. file_column_for_param(i) > max_col_num) then
                  write(*,1) 'need to set file_column_for_param '//trim(param_name(i)), i
                  ierr = -1
                  return
               end if
               if (file_column_for_param(i) > num_to_read) num_to_read = file_column_for_param(i)
            end if
         end do

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
            
            do i = 1, max_parameters
               if (param_name(i) /= '' .and. vary_param(i)) then
                  next_param_to_try(i) = filedata(file_column_for_param(i))
                  write(*,2) 'next_param_to_try '//trim(param_name(i)), i, next_param_to_try(i)
               end if
            end do
            
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
            character(len=256) :: filename
            
            include 'formats'
            ierr = 0
                     
            chi2 = eval1(s% id,ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in eval1'
               return
            end if
            
            call save_best_for_sample(sample_number, 0)

            if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))
            filename = trim(astero_results_directory) // '/' // trim(from_file_output_filename)
         
            call save_sample_results_to_file(-1, filename, ierr)
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
         
         integer :: num_param(1:max_parameters)
         integer :: i_total, i
         real(dp) :: chi2, param(1:max_parameters)
         real(dp), parameter :: eps = 1d-6
         logical :: just_counting
         
         include 'formats'
         
         ierr = 0
         call set_starting_values
         
         num_param(1:max_parameters) = 0
         
         just_counting = .true.
         call do_param(max_parameters, ierr)
         i_total = sample_number
         
         write(*,2) 'grid total', i_total

         do i = 1, max_parameters
            if (param_name(i) /= '') write(*,3) 'num_param(i) '//trim(param_name(i)), i, num_param(i)
         end do

         write(*,'(A)')
         
         sample_number = 0
         just_counting = .false.
         
         if (restart_scan_grid_from_file) then
            call read_samples_from_file( &
               trim(astero_results_directory) // '/' // trim(scan_grid_output_filename), ierr)
            if (ierr /= 0) return
            scan_grid_skip_number = sample_number
            sample_number = 0
            write(*,2) 'scan_grid_skip_number', scan_grid_skip_number
         else
            scan_grid_skip_number = 0
         end if

         call do_param(max_parameters, ierr)
                 
                 
         contains
         
         
         subroutine set_starting_values
            ! this looks like a bug: only set varied parameters to mins; rest to first?
            param(1:max_parameters) = min_param(1:max_parameters)
         end subroutine set_starting_values


         recursive subroutine do_param(k, ierr)
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            integer :: cnt
            cnt = 0
            ierr = 0
            do while (param(k) <= max_param(k) + eps .or. .not. vary_param(k))
               if (vary_param(k)) next_param_to_try(k) = param(k)
               if (k == 1) then
                  call do1_grid(ierr)
               else
                  call do_param(k-1, ierr)
               end if
               if (ierr /= 0) return
               cnt = cnt+1
               if (delta_param(k) <= 0 .or. .not. vary_param(k)) exit
               param(k) = param(k) + delta_param(k)
            end do
            if (num_param(k) == 0) num_param(k) = cnt
            param(k) = min_param(k)
         end subroutine do_param
         
         
         subroutine do1_grid(ierr)
            integer, intent(out) :: ierr
            character(len=256) :: filename
            include 'formats'
            ierr = 0
            
            if (just_counting) then
               sample_number = sample_number + 1
               return
            end if
            
            if (sample_number < scan_grid_skip_number) then
               sample_number = sample_number + 1
               if (mod(sample_number,20) == 1) then
                  write(*,'(66x,a1)',advance='no') ' '
                  do i = 1, max_parameters
                     if (param_name(i) /= '') write(*,'(a26)',advance='no') trim(param_name(i))
                  end do
                  write(*,*) ''
               end if

               write(*,2,advance='no') 'restore sample from file', sample_number
               do i = 1, max_parameters
                  if (param_name(i) /= '') write(*,'(1pd26.16)',advance='no') param(i)
               end do
               write(*,*) ''

               return
            end if

            write(*,2) 'eval1 sample_number', sample_number+1
            write(*,'(a)') ''

            do i = 1, max_parameters
               if (param_name(i) /= '') write(*,1) trim(param_name(i)), param(i)
            end do

            write(*,'(a)') ''

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

            if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))
            filename = trim(astero_results_directory) // '/' // trim(scan_grid_output_filename)

            call save_sample_results_to_file(i_total, filename, ierr)
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
         
         character(len=256) :: filename
         integer :: ierr

         call bobyqa_or_newuoa_fun(n,x,f)

         write(*,'(A)')
         ierr = 0

         if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))
         filename = trim(astero_results_directory) // '/' // trim(bobyqa_output_filename)

         call save_sample_results_to_file(-1, filename, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in save_sample_results_to_file'
            call mesa_error(__FILE__,__LINE__,'bobyqa_fun')
         end if
         
      end subroutine bobyqa_fun
      

      subroutine newuoa_fun(n,x,f)
         integer, intent(in) :: n
         double precision, intent(in) :: x(*)
         double precision, intent(out) :: f
         
         character(len=256) :: filename
         integer :: ierr
         
         call bobyqa_or_newuoa_fun(n,x,f)

         write(*,'(A)')
         ierr = 0

         if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))
         filename = trim(astero_results_directory) // '/' // trim(newuoa_output_filename)

         call save_sample_results_to_file(-1, filename, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in save_sample_results_to_file'
            call mesa_error(__FILE__,__LINE__,'newuoa_fun')
         end if
         
      end subroutine newuoa_fun


      subroutine bobyqa_or_newuoa_fun(n,x,f)
         integer, intent(in) :: n
         double precision, intent(in) :: x(*)
         double precision, intent(out) :: f
         integer :: ierr, prev_sample_number, i
         include 'formats'
         
         ierr = 0

         do i = 1, max_parameters
            if (param_name(i) /= '' .and. vary_param(i)) then
               next_param_to_try(i) = bobyqa_param( &
                  x(i_param(i)), first_param(i), min_param(i), max_param(i))
               write(*,2) 'next_param_to_try '//trim(param_name(i)), i, next_param_to_try(i), x(i_param(i))
            end if
         end do
         
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

         write(*,'(A)')
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
         
         write(*,'(A)')
         write(*,'(A)')

         do i = 1, max_parameters
            if (param_name(i) /= '' .and. vary_param(i)) then
               nvar = nvar+1; i_param(i) = nvar
               if (min_param(i) >= max_param(i)) then
                  write(*,2) 'min_param >= max_param', &
                     i, min_param(i), max_param(i)
                  ierr = -1
               end if
            end if
         end do

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

         do i = 1, max_parameters
            if (param_name(i) /= '' .and. vary_param(i)) then
               final_param(i) = bobyqa_param( &
                  x(i_param(i)), first_param(i), min_param(i), max_param(i))
            end if
         end do
         
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
         
         character(len=256) :: filename
         integer :: prev_sample_number, i
         include 'formats'
         
         ierr = 0
         
         write(*,'(A)')
         write(*,'(A)')

         do i = 1, max_parameters
            if (vary_param(i)) then
               next_param_to_try(i) = simplex_param( &
                  x(i_param(i)), first_param(i), min_param(i), max_param(i))
               write(*,2) 'next_param_to_try '//trim(param_name(i)), &
                  i, next_param_to_try(i), x(i_param(i))
            end if
         end do
         
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

         if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))
         filename = trim(astero_results_directory) // '/' // trim(simplex_output_filename)

         call save_sample_results_to_file(-1, filename, ierr)
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
              final_param(1:max_parameters)
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
         
         do i = 1, max_parameters
            if (param_name(i) /= '' .and. vary_param(i)) then
               nvar = nvar+1; i_param(i) = nvar
               if (min_param(i) >= max_param(i)) then
                  write(*,2) 'min_param >= max_param '//trim(param_name(i)), &
                     i, min_param(i), max_param(i)
                  ierr = -1
               end if
            end if
         end do

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
            call read_samples_from_file( &
               trim(astero_results_directory) // '/' // trim(simplex_output_filename), ierr)
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

         do i = 1, max_parameters
            if (param_name(i) /= '' .and. vary_param(i)) then
               final_param(i) = simplex_param( &
                  x_final(i_param(i)), first_param(i), &
                  min_param(i), max_param(i))
            end if
         end do

         deallocate( &
            rpar, ipar, simplex, f, x_lower, x_upper, x_first, x_final)
            
            
         contains
         
         
         subroutine set_xs ! x_first, x_lower, x_upper

            do i = 1, max_parameters
               if (vary_param(i)) then
                  x_first(i_param(i)) = first_param(i)
                  x_lower(i_param(i)) = min_param(i)
                  x_upper(i_param(i)) = max_param(i)
               end if
            end do

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

               do k = 1, max_parameters
                  if (param_name(k) /= '' .and. vary_param(k)) then
                     simplex(i_param(k),j) = &
                        simplex_inverse(sample_param(k,i), first_param(k), min_param(k), max_param(k))
                     write(*,4) 'param '//trim(param_name(i)), k, j, i, sample_param(k,i)
                  end if
               end do

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
            write(*,'(A)')
            write(*,'(A)')
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
                  
         sample_constraint_value(1:max_constraints,i) = best_constraint_value(1:max_constraints)

         sample_op_code(i) = op_code
         sample_chi2(i) = best_chi2
         sample_chi2_seismo(i) = best_chi2_seismo
         sample_chi2_spectro(i) = best_chi2_spectro
         
         sample_age(i) = best_age

         sample_param(1:max_parameters,i) = current_param(1:max_parameters)

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
         
         write(*,'(A)')
         write(*,2) 'n for averages', n
         write(*,1) 'avg_age_top_samples', avg_age_top_samples
         write(*,1) 'avg_age_sigma', avg_age_sigma
         write(*,1) 'age limit', avg_age_top_samples + avg_age_sigma_limit*avg_age_sigma
         write(*,1) 'avg_model_number_top_samples', avg_model_number_top_samples
         write(*,1) 'avg_model_number_sigma', avg_model_number_sigma
         write(*,1) 'model number limit', &
            avg_model_number_top_samples + &
               avg_model_number_sigma_limit*avg_model_number_sigma
         write(*,'(A)')
         !call mesa_error(__FILE__,__LINE__,'set_sample_averages')
         
      end subroutine set_sample_averages
      
      
      subroutine zero_best_info
         best_chi2 = 0
         best_chi2_seismo = 0
         best_chi2_spectro = 0
         best_age = 0
         best_constraint_value(1:max_constraints) = 0
         best_param(1:max_parameters) = 0
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
