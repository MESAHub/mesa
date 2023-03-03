! ***********************************************************************
!
!   Copyright (C) 2013-2022  The MESA Team, Pablo Marchant & Matthias Fabry
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

module run_binary_support

   use star_lib
   use star_def
   use const_def
   use utils_lib
   use binary_def
   use binary_private_def
   use binary_ctrls_io, only : do_one_binary_setup
   use binary_do_one_utils
   use binary_ce
   use binary_photos

   implicit none

contains

   subroutine do_run1_binary(tst, &
      ! star extras
      extras_controls, &
      ! binary extras
      extras_binary_controls, &
      ierr, &
      inlist_fname_arg)

      use binary_job_ctrls_io
      use binary_mdot, only : adjust_mdots, set_accretion_composition
      use binary_tides, only : sync_spin_orbit_torque
      use binary_evolve
      use mod_other_rlo_mdot
      use mod_other_implicit_rlo
      use mod_other_tsync
      use mod_other_sync_spin_to_orbit
      use mod_other_mdot_edd
      use mod_other_adjust_mdots
      use mod_other_accreted_material_j
      use mod_other_binary_jdot
      use mod_other_binary_wind_transfer
      use mod_other_binary_edot
      use mod_other_binary_ce
      use mod_other_binary_extras
      use mod_other_binary_photo_read
      use mod_other_binary_photo_write
      use mod_other_e2
      use mod_other_pgbinary_plots
      use binary_timestep
      use binary_history
      use binary_history_specs
      use run_star_support
      use pgbinary, only : read_pgbinary_inlist, update_pgbinary_plots, &
         start_new_run_for_pgbinary, restart_run_for_pgbinary

      logical, intent(in) :: tst

      interface

         subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine extras_controls

         subroutine extras_binary_controls(binary_id, ierr)
            integer :: binary_id
            integer, intent(out) :: ierr
         end subroutine extras_binary_controls

      end interface

      integer, intent(out) :: ierr
      character (len = *) :: inlist_fname_arg
      optional inlist_fname_arg

      integer :: id, binary_id, i, j, k, l, i_prev, result, partial_result, &
         result_reason, model_number, iounit, binary_startup, model, num_stars
      type (star_info), pointer :: s
      character (len = 256) :: restart_filename, photo_filename
      integer(8) :: total, time0, time1, clock_rate
      logical :: doing_restart, first_try, continue_evolve_loop, &
          get_history_info, write_history, write_terminal, will_read_pgbinary_inlist
      real(dp) :: sum_times, dt, timestep_factor
      type (binary_info), pointer :: b
      character (len = strlen) :: inlist_fname

      include 'formats'

      ierr = 0
      call system_clock(time0, clock_rate)

      call resolve_inlist_fname(inlist_fname, inlist_fname_arg)
      MESA_INLIST_RESOLVED = .true. ! Now any call to resolve_inlist_fname will only return inlist_fname_arg or 'inlist'

      ! Find out if this is a restart
      open(newunit = iounit, file = '.restart', status = 'old', action = 'read', iostat = ierr)
      doing_restart = (ierr == 0)
      if (doing_restart) then
         ! photo_filename will be like 'x100'
         ! the photo for star 1 is '1_x100', for star 2 '2_x100' and for binary 'b_x100'
         read(iounit, '(a)', iostat = ierr) photo_filename
         if (ierr /= 0) then
            stop "Problem while reading restart info"
         end if
      else
         ierr = 0
      end if

      binary_id = alloc_binary(ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in alloc_binary'
         return
      end if

      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      if (.not. doing_restart) then
         b% model_number = 0
         b% model_number_old = 0
         b% binary_age = 0
         b% binary_age_old = 0
      end if

      call do_read_binary_job(b, inlist_fname, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in read_binary_job'
         return
      end if

      if (.not. b% job% evolve_both_stars .or. &
         ((b% job% change_initial_model_twins_flag .or. b% job% change_model_twins_flag) &
            .and. b% job% new_model_twins_flag)) then
         b% have_star_1 = .true.
         b% have_star_2 = .false.
      else
         b% have_star_1 = .true.
         b% have_star_2 = .true.
      end if

      write(*, '(A)')
      write(*, '(A)')

      result_reason = 0

      ! Setup null hooks
      b% other_rlo_mdot => null_other_rlo_mdot
      b% other_implicit_function_to_solve => null_other_implicit_function_to_solve
      b% other_check_implicit_rlo => null_other_check_implicit_rlo
      b% other_tsync => null_other_tsync
      b% other_sync_spin_to_orbit => null_other_sync_spin_to_orbit
      b% other_mdot_edd => null_other_mdot_edd
      b% other_adjust_mdots => null_other_adjust_mdots
      b% other_accreted_material_j => null_other_accreted_material_j
      b% other_jdot_gr => null_other_jdot_gr
      b% other_jdot_ml => null_other_jdot_ml
      b% other_jdot_ls => null_other_jdot_ls
      b% other_jdot_missing_wind => null_other_jdot_missing_wind
      b% other_jdot_mb => null_other_jdot_mb
      b% other_extra_jdot => null_other_extra_jdot
      b% other_binary_wind_transfer => null_other_binary_wind_transfer
      b% other_edot_tidal => null_other_edot_tidal
      b% other_edot_enhance => null_other_edot_enhance
      b% other_extra_edot => null_other_extra_edot
      b% other_CE_init => null_other_CE_init
      b% other_CE_rlo_mdot => null_other_CE_rlo_mdot
      b% other_CE_binary_evolve_step => null_other_CE_binary_evolve_step
      b% other_CE_binary_finish_step => null_other_CE_binary_finish_step
      b% other_e2 => null_other_e2
      b% other_pgbinary_plots_info => null_other_pgbinary_plots_info

      b% extras_binary_startup => null_extras_binary_startup
      b% extras_binary_start_step => null_extras_binary_start_step
      b% extras_binary_check_model => null_extras_binary_check_model
      b% extras_binary_finish_step => null_extras_binary_finish_step
      b% extras_binary_after_evolve => null_extras_binary_after_evolve
      b% how_many_extra_binary_history_columns => null_how_many_extra_binary_history_columns
      b% data_for_extra_binary_history_columns => null_data_for_extra_binary_history_columns
      b% how_many_extra_binary_history_header_items => null_how_many_extra_binary_history_header_items
      b% data_for_extra_binary_history_header_items => null_data_for_extra_binary_history_header_items

      b% other_binary_photo_read => default_other_binary_photo_read
      b% other_binary_photo_write => default_other_binary_photo_write

      b% ignore_hard_limits_this_step = .false.

      call do_one_binary_setup(b, inlist_fname, ierr)
      ! extras_binary_controls is defined in run_binary_extras.f and hooks can
      ! be specified there
      call extras_binary_controls(b% binary_id, ierr)

      ! load binary photo if this is a restart
      if (doing_restart) then
         restart_filename = trim(trim(b% photo_directory) // '/b_' // photo_filename)
         call binary_load_photo(b, restart_filename, ierr)
         if (failed('binary_load_photo', ierr)) return
      end if

      b% donor_id = -1
      b% accretor_id = -1
      do i = 1, 2

         if (i==1 .and. .not. b% have_star_1) then
            cycle
         else if (i==2 .and. .not. b% have_star_2) then
            cycle
         end if

         call do_read_star_job(b% job% inlist_names(i), ierr)
         if (failed('do_read_star_job', ierr)) return

         if (i==1) then
            restart_filename = trim('1_' // photo_filename)
         else
            restart_filename = trim('2_' // photo_filename)
         end if

         ! the star is initialized in this call
         call before_evolve_loop(.true., doing_restart, doing_restart, &
            binary_controls, extras_controls, &
            id_from_read_star_job, b% job% inlist_names(i), restart_filename, &
            .false., binary_id, id, ierr)

         if (ierr /= 0) return

         call star_ptr(id, s, ierr)
         if (failed('star_ptr', ierr)) return
         b% star_ids(i) = id

         s% include_binary_history_in_log_file = b% append_to_star_history

         s% how_many_binary_history_columns => how_many_binary_history_columns
         s% data_for_binary_history_columns => data_for_binary_history_columns

         s% how_many_extra_binary_history_columns => b% how_many_extra_binary_history_columns
         s% data_for_extra_binary_history_columns => b% data_for_extra_binary_history_columns

         ! additional settings for mass transfer and tides
         if (b% do_j_accretion) then
            s% use_accreted_material_j = .true.
         end if
         s% accrete_given_mass_fractions = .true.
         s% accrete_same_as_surface = .false.
         s% binary_other_torque => sync_spin_orbit_torque

         s% doing_timing = .false.

         write(*, '(A)')
         write(*, '(A)')

      end do

      ! Error if users set the saved model to the same filename
      if(b% have_star_1 .and. b% have_star_2) then
         if(b% s1% job% save_model_when_terminate .and. b% s2% job% save_model_when_terminate) then
            if(len_trim(b% s1% job% save_model_filename) > 0 .and. len_trim(b% s2% job% save_model_filename) >0) then
               if(trim(b% s1% job% save_model_filename) == trim(b% s2% job% save_model_filename)) then
                  write(*, *) "ERROR: Both stars are set to write save_model_filename to the same file"
                  call mesa_error(__FILE__, __LINE__)
               end if
            end if
         end if
      end if



      ! binary data must be initiated after stars, such that masses are available
      ! if using saved models
      call binarydata_init(b, doing_restart)
      call binary_private_def_init
      call binary_history_column_names_init(ierr)
      call set_binary_history_columns(b, b% job% binary_history_columns_file, ierr)

      ! setup pgbinary
      if (.not. doing_restart) then
         call start_new_run_for_pgbinary(b, ierr)
         if (failed('start_new_run_for_pgbinary', ierr)) return
      else
         call restart_run_for_pgbinary(b, ierr)
         if (failed('restart_run_for_pgbinary', ierr)) return
      end if

      if (b% job% show_binary_log_description_at_start .and. .not. doing_restart) then
         write(*, '(A)')
         call do_show_binary_log_description(id, ierr)
         if (failed('show_log_description', ierr)) return
      end if

      if (b% point_mass_i /= 1 .and. b% do_initial_orbit_sync_1 .and. &
         .not. doing_restart) then
         call star_relax_uniform_omega(&
            b% s1% id, 0, (2 * pi) / b% period, b% s1% job% num_steps_to_relax_rotation, &
            b% s1% job% relax_omega_max_yrs_dt, ierr)
         if (ierr /= 0) then
            write(*, *) 'failed in initial orbital sync'
            return
         end if
      end if

      if (b% point_mass_i /= 2 .and. b% do_initial_orbit_sync_2 .and. &
         .not. doing_restart) then
         call star_relax_uniform_omega(&
            b% s2% id, 0, (2 * pi) / b% period, b% s2% job% num_steps_to_relax_rotation, &
            b% s2% job% relax_omega_max_yrs_dt, ierr)
         if (ierr /= 0) then
            write(*, *) 'failed in initial orbital sync'
            return
         end if
      end if

      continue_evolve_loop = .true.
      s% doing_timing = .false.
      i_prev = 0

      binary_startup = b% extras_binary_startup(b% binary_id, doing_restart, ierr)
      if (ierr /= 0) return

      ! perform changes specified by binary_job variables
      call do_binary_job_controls_after(b% binary_id, b, doing_restart, ierr)
      if (ierr /= 0) return

      ! setup things for common envelope in case of a restart
      if (doing_restart .and. b% CE_flag .and. b% CE_init) then
         call ce_init(b, doing_restart, ierr)
      end if

      evolve_loop : do while(continue_evolve_loop) ! evolve one step per loop

         if (b% point_mass_i /= 0) then
            num_stars = 1
         else
            num_stars = 2
         end if

         do l = 1, num_stars
            if (l == 1 .and. b% point_mass_i == 1) then
               i = 2
            else
               i = l
            end if

            id = b% star_ids(i)
            call star_ptr(id, s, ierr)
            if (failed('star_ptr', ierr)) return
            call before_step_loop(id, ierr)
            if (ierr /= 0) return

            result = s% extras_start_step(id)
            if (result /= keep_going) then
               continue_evolve_loop = .false.
               exit evolve_loop
            end if
         end do

         first_try = .true.
         b% donor_started_implicit_wind = .false.
         b% num_tries = 0

         if (b% CE_flag .and. .not. b% CE_init) then
            call ce_init(b, .false., ierr)
            write(*, *) "CE flag is on!!"
         end if

         step_loop : do ! may need to repeat this loop

            result = b% extras_binary_start_step(b% binary_id, ierr)
            if (ierr /= 0) then
               write(*, *) "error in extras_binary_start_step"
               return
            end if
            if (result /= keep_going) exit evolve_loop
            !update num_stars in case point_mass_i has changed
            if (b% point_mass_i /= 0) then
               num_stars = 1
            else
               num_stars = 2
            end if

            call set_donor_star(b)
            call set_star_timesteps(b)
            result = keep_going

            ! Store mtransfer_rate used in a step, as it is rewritten by binary_check_model and
            ! that can produce inconsistent output.
            b% step_mtransfer_rate = b% mtransfer_rate

            if (.not. b% CE_flag) then
               ! if donor reaches implicit wind limit during mass transfer, set was_in_implicit_wind_limit = .true.
               ! so that further iterations of the implicit wind do not start far off from the required mdot.
               ! system will likely detach suddenly, so set change_factor to double its maximum
               if (b% donor_started_implicit_wind) then
                  b% s_donor% was_in_implicit_wind_limit = .true.
                  b% change_factor = 2 * b% max_change_factor
               end if
            end if

            do l = 1, num_stars

               if (b% point_mass_i /= 0) then  ! if any point mass, evolve the other one
                  j = 3 - b% point_mass_i
               else  ! both stars
                  if (l == 1) then
                     j = b% d_i  ! evolve the donor first
                  else
                     j = b% a_i
                  end if
               end if

               id = b% star_ids(j)

               ! Avoid repeting the accretor when using the implicit scheme plus
               ! rotation and implicit winds. When this happens the accretor won't
               ! usually care about the result of the evolution of the donor.
               if (j == b% a_i .and. b% num_tries >0 .and. s% was_in_implicit_wind_limit) &
                  cycle

               ! fix sync timescales to zero. If synching stars these will be
               ! updated at each star_evolve_step
               if (j == 1) then
                  b% t_sync_1 = 0
               else if (j == 2) then
                  b% t_sync_2 = 0
               end if

               result = worst_result(result, &
                  star_evolve_step_part1(id, first_try))

            end do

            ! modify mdots to account for mass transfer
            call adjust_mdots(b)
            ! if both stars are accreting mass at this point, then there is
            ! something very wrong! If one star loses and the other gains mass,
            ! then the mass losing star must be evolved first
            k = b% d_i
            if (b% point_mass_i == 0 .and. b% s_donor% mstar_dot > 0) then
               if (b% s_accretor% mstar_dot > 0) then
                  write(*, *) "ERROR: both stars accreting, terminating evolution"
                  result = terminate
                  exit step_loop
               end if
               k = b% a_i  ! donor is gaining while accretor is losing, accretor plays donor in this step
            else if (b% point_mass_i /= 0 .and. b% s_donor% mstar_dot > 0) then
               write(*, *) "ERROR: donor accreting, terminating evolution"
               result = terminate
               exit step_loop
            end if

            if (result == keep_going) then
               do l = 1, num_stars
                  ! if there's any point mass, evolve the non-point mass, num_stars should be 1 here
                  ! so no risk of evolving things twice
                  if (b% point_mass_i /= 0) then
                     k = 3 - b% point_mass_i
                  else  ! two stars present
                     if (l == 1) then  ! evolve donor in first loop pass
                        j = k
                     else  ! evolve acc in second loop pass
                        j = 3 - k
                     end if
                  end if
                  id = b% star_ids(j)

                  if (j == b% a_i .and. b% num_tries >0 .and. s% was_in_implicit_wind_limit) &
                     cycle

                  ! set accretion composition
                  if (.not. b% CE_flag) then
                     if (i == 2) call set_accretion_composition(b, j)
                  end if

                  result = worst_result(result, &
                     star_evolve_step_part2(id, first_try))

               end do
            end if

            ! do not evolve binary step on failure, its unnecesary and some variables are not properly
            ! set when the newton solver fails.
            if (result == keep_going) then
               if (.not. b% CE_flag) then
                  result = worst_result(result, binary_evolve_step(b))
               else
                  result = worst_result(result, CE_binary_evolve_step(b))
               end if
            end if

            if (result == keep_going) then
               result = worst_result(result, binary_check_model(b))
            end if

            do l = 1, num_stars
               if (l == 1 .and. b% point_mass_i == 1) then
                  i = 2
               else
                  i = l
               end if
               if (i == 1 .and. b% point_mass_i == 1) i = 2
               if (result == keep_going) then
                  id = b% star_ids(i)
                  call star_ptr(id, s, ierr)
                  result = worst_result(result, s% extras_check_model(id))
                  result = worst_result(result, star_check_model(id))
               end if
            end do

            partial_result = b% extras_binary_check_model(b% binary_id)
            result = worst_result(result, partial_result)

            ! solve first binary timestep limit because star_pick_next_timestep needs it
            ! check redos as well to stop the implicit solver if a hard limit is hit
            if (result == keep_going .or. result == redo) then
               result = worst_result(result, binary_pick_next_timestep(b))
            end if

            if (result == keep_going) then
               do l = 1, num_stars
                  if (l == 1 .and. b% point_mass_i == 1) then
                     i = 2
                  else
                     i = l
                  end if
                  id = b% star_ids(i)
                  call star_ptr(id, s, ierr)
                  if (failed('star_ptr', ierr)) return
                  result = worst_result(result, star_pick_next_timestep(id))
               end do
            end if
            if (result == keep_going) then
               exit step_loop
            end if

            do l = 1, num_stars
               if (l == 1 .and. b% point_mass_i == 1) then
                  i = 2
               else
                  i = l
               end if

               id = b% star_ids(i)
               model_number = get_model_number(id, ierr)
               if (failed('get_model_number', ierr)) return

               result_reason = get_result_reason(id, ierr)
               if (result == retry) then
                  if (failed('get_result_reason', ierr)) return
                  if (s% job% report_retries) &
                     write(*, '(i6,3x,a,/)') model_number, &
                        'retry reason ' // trim(result_reason_str(result_reason))
               end if

            end do

            b% generations = b% s_donor% generations
            if (result == redo) then
               do l = 1, num_stars
                  if (l == 1 .and. b% point_mass_i == 1) then
                     i = 2
                  else
                     i = l
                  end if
                  id = b% star_ids(i)

                  ! Avoid repeting the accretor when using the implicit scheme plus
                  ! rotation and implicit winds. When this happens the accretor won't
                  ! usually care about the result of the evolution of the donor.

                  if (i == b% a_i .and. b% num_tries >0 .and. s% was_in_implicit_wind_limit) &
                     cycle
                  result = worst_result(result, star_prepare_to_redo(id))
               end do
               result = worst_result(result, binary_prepare_to_redo(b))
            end if
            if (result == retry) then
               do l = 1, num_stars
                  if (l == 1 .and. b% point_mass_i == 1) then
                     i = 2
                  else
                     i = l
                  end if
                  id = b% star_ids(i)
                  result = worst_result(result, star_prepare_to_retry(id))
               end do
               result = worst_result(result, binary_prepare_to_retry(b))
            end if

            !correct pointers to stars after retry
            if (b% d_i == 1) then
               if (b% have_star_1) then
                  b% s_donor => b% s1
                  if (b% point_mass_i == 0) then
                     if (b% have_star_2) then
                        b% s_accretor => b% s2
                     else
                        call mesa_error(__FILE__, __LINE__, 'ERROR: missing star pointer for accretor')
                     end if
                  end if
               else
                  call mesa_error(__FILE__, __LINE__, 'ERROR: missing star pointer for donor')
               end if
            else
               if (b% have_star_2) then
                  b% s_donor => b% s2
                  if (b% point_mass_i == 0) then
                     if (b% have_star_1) then
                        b% s_accretor => b% s1
                     else
                        call mesa_error(__FILE__, __LINE__, 'ERROR: missing star pointer for accretor')
                     end if
                  end if
               else
                  call mesa_error(__FILE__, __LINE__, 'ERROR: missing star pointer for donor')
               end if
            end if

            if (result == terminate) then
               continue_evolve_loop = .false.
               exit step_loop
            end if
            first_try = .false.

         end do step_loop

         if(result == keep_going) result = binary_finish_step(b)
         if (b% CE_flag .and. b% CE_init .and. result == keep_going) then
            result = worst_result(result, CE_binary_finish_step(b))
         end if

         partial_result = b% extras_binary_finish_step(b% binary_id)
         result = worst_result(result, partial_result)

         if (result == keep_going) then
            ! write terminal info
            model = b% model_number
            if (b% history_interval > 0) then
               write_history = (mod(model, b% history_interval) == 0)
            else
               write_history = .false.
            end if
            if (s% terminal_interval > 0) then
               write_terminal = (mod(model, b% terminal_interval) == 0)
            else
               write_terminal = .false.
            end if
            if (write_history) b% need_to_update_binary_history_now = .true.
            get_history_info = b% need_to_update_binary_history_now .or. write_terminal
            if (get_history_info) then
               if (b% write_header_frequency * b% terminal_interval > 0) then
                  if (mod(model, b% write_header_frequency * b% terminal_interval) == 0 &
                     .and. .not. b% doing_first_model_of_run) then
                     write(*, '(A)')
                     call write_binary_terminal_header(b)
                  end if
               end if
               if (write_terminal) call do_binary_terminal_summary(b)
               if (b% need_to_update_binary_history_now) call write_binary_history_info(b, ierr)
               b% need_to_update_binary_history_now = .false.
            end if
         end if

         do l = 1, num_stars
            if (l == 1 .and. b% point_mass_i == 1) then
               i = 2
            else
               i = l
            end if
            id = b% star_ids(i)
            call star_ptr(id, s, ierr)
            partial_result = result
            call after_step_loop(id, b% job% inlist_names(i), &
               .false., partial_result, ierr)
            if (ierr /= 0) return
            result = worst_result(result, partial_result)
         end do

         ! do pgbinary stuff
         if (result == keep_going .and. b% job% pgbinary_flag) then
            will_read_pgbinary_inlist = .false.
            if (b% pg% pgbinary_interval <= 0) then
               will_read_pgbinary_inlist = .true.
            else if(mod(b% model_number, b% pg% pgbinary_interval) == 0) then
               will_read_pgbinary_inlist = .true.
            end if
            if (will_read_pgbinary_inlist) then
               call read_pgbinary_inlist(b, inlist_fname, ierr)
               if (failed('read_pgbinary_controls', ierr)) return
               if (b% point_mass_i /= 1) &
                  call read_pgstar_inlist(b% s1, b% job% inlist_names(1), ierr)
               if (failed('read_pgstar_controls 1', ierr)) return
               if (b% point_mass_i /= 2) &
                  call read_pgstar_inlist(b% s2, b% job% inlist_names(2), ierr)
               if (failed('read_pgstar_controls 2', ierr)) return
            end if
         end if
         if (result == keep_going .and. b% job% pgbinary_flag) then
            call update_pgbinary_plots(b, .false., ierr)
            if (failed('update_pgbinary_plots', ierr)) return
         end if

         if (result /= keep_going) then
            if (result /= terminate) then
               write(*, 2) 'ERROR in result value in run_star_extras: model', &
                  s% model_number
               write(*, 2) 'result', result
               exit evolve_loop
            end if
            do l = 1, num_stars
               if (l == 1 .and. b% point_mass_i == 1) then
                  i = 2
               else
                  i = l
               end if
               id = b% star_ids(i)
               call star_ptr(id, s, ierr)
               if (s% result_reason == result_reason_normal) then

                  partial_result = result
                  call terminate_normal_evolve_loop(id, &
                     .false., partial_result, ierr)
                  if (ierr /= 0) return
                  result = worst_result(result, partial_result)

               end if
            end do
            call write_binary_history_info(b, ierr)
            call do_binary_terminal_summary(b)
            exit evolve_loop
         end if

         do l = 1, num_stars
            if (l == 1 .and. b% point_mass_i == 1) then
               i = 2
            else
               i = l
            end if
            id = b% star_ids(i)
            call star_ptr(id, s, ierr)

            call do_saves(&
               id, ierr)
            if (ierr /= 0) return

            if (s% doing_timing) then
               call system_clock(s% job% time1_extra, s% job% clock_rate)
               s% job% after_step_timing = s% job% after_step_timing + &
                  dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate
               s% job% check_time_end = eval_total_times(s% id, ierr)
               s% job% check_after_step_timing = s% job% check_after_step_timing + &
                  (s% job% check_time_end - s% job% check_time_start)
            end if
         end do
         if (b% photo_interval > 0 .and. mod(b% model_number, b% photo_interval) == 0) then
            call do_saves_for_binary(b, ierr)
            if (ierr /= 0) return
         end if

         if (b% doing_first_model_of_run) b% doing_first_model_of_run = .false.

      end do evolve_loop

      ! save photos at end of the run
      call do_saves_for_binary(b, ierr)
      if (ierr /= 0) return

      ! deallocate arrays used for calculation of phase dependent variables
      deallocate(b% theta_co, b% time_co, b% mdot_donor_theta)
      deallocate(b% edot_theta, b% e1, b% e2, b% e3)

      ierr = binary_after_evolve(b)
      if (ierr /= 0) then
         write(*, *) "Error in binary_after_evolve"
         return
      end if
      call b% extras_binary_after_evolve(b% binary_id, ierr)
      if (ierr /= 0) then
         write(*, *) "Error in extras_binary_after_evolve"
         return
      end if

      do i = 1, num_stars
         if (i == 1 .and. b% point_mass_i == 1) then
            l = 2
         else
            l = i
         end if

         id = b% star_ids(l)

         call star_ptr(id, s, ierr)
         if (failed('star_ptr', ierr)) then
            ierr = 0
            cycle
         end if

         call after_evolve_loop(id, .true., ierr)

         if (s% doing_timing) then
            call system_clock(s% job% time1_extra, s% job% clock_rate)
            s% job% after_step_timing = s% job% after_step_timing + &
               dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate
            s% job% check_time_end = eval_total_times(s% id, ierr)
            s% job% check_after_step_timing = s% job% check_after_step_timing + &
               (s% job% check_time_end - s% job% check_time_start)
         end if

      end do

      call starlib_shutdown

   end subroutine do_run1_binary

   integer function worst_result(result1, result2)
      integer, intent(in) :: result1, result2

      if(result1 == terminate .or. result2 == terminate) then
         worst_result = terminate
         return
      end if

      if(result1 == retry .or. result2 == retry) then
         worst_result = retry
         return
      end if

      if(result1 == redo .or. result2 == redo) then
         worst_result = redo
         return
      end if

      worst_result = keep_going
      return

   end function worst_result

   subroutine binary_controls(id, binary_id, ierr)
      integer, intent(in) :: id, binary_id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      ierr = 0

      call star_ptr(id, s, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in star_ptr'
         return
      end if
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      s% binary_id = binary_id
      ! binary takes care of writing photos
      s% photo_interval = -1
      s% job% save_photo_when_terminate = .false.
      s% photo_digits = b% photo_digits

      if (b% donor_id == -1) then
         b% donor_id = id
         b% s1 => s
         s% initial_mass = b% m1
      else
         b% accretor_id = id
         b% s2 => s
         s% initial_mass = b% m2
      end if
   end subroutine binary_controls

   subroutine do_binary_job_controls_after(binary_id, b, restart, ierr)
      use binary_utils
      use binary_evolve, only : binary_finish_step
      integer, intent(in) :: binary_id
      type (binary_info), pointer :: b
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      logical :: need_to_fix_generations

      ! If something is changed in here, fix generations to zero so a retry wont
      ! erase change
      need_to_fix_generations = .false.

      if (b% job% change_ignore_rlof_flag .or. &
         (b% job% change_initial_ignore_rlof_flag .and. .not. restart)) then
         if (b% job% new_ignore_rlof_flag .neqv. b% ignore_rlof_flag) then
            call set_ignore_rlof_flag(binary_id, b% job% new_ignore_rlof_flag, ierr)
            if (ierr /= 0) then
               write(*, *) "error in set_ignore_rlof_flag"
               return
            end if
            need_to_fix_generations = .true.
         end if
      end if

      if (b% job% change_model_twins_flag .or. &
         (b% job% change_initial_model_twins_flag .and. .not. restart)) then
         if (b% job% new_model_twins_flag .neqv. b% model_twins_flag) then
            call set_model_twins_flag(binary_id, b% job% new_model_twins_flag, ierr)
            if (ierr /= 0) then
               write(*, *) "error in set_model_twins_flag"
               return
            end if
            need_to_fix_generations = .true.
         end if
      end if

      if (b% job% change_point_mass_i .or. &
         (b% job% change_initial_point_mass_i .and. .not. restart)) then
         if (b% job% new_point_mass_i /= b% point_mass_i) then
            call set_point_mass_i(binary_id, b% job% new_point_mass_i, ierr)
            if (ierr /= 0) then
               write(*, *) "error in set_point_mass_i"
               return
            end if
            need_to_fix_generations = .true.
         end if
      end if

      if (b% job% change_m1 .or. &
         (b% job% change_initial_m1 .and. .not. restart)) then
         if (b% m(1) /= b% job% new_m1) then
            call set_m1(binary_id, b% job% new_m1, ierr)
            if (ierr /= 0) then
               write(*, *) "error in set_m1"
               return
            end if
            need_to_fix_generations = .true.
         end if
      end if

      if (b% job% change_m2 .or. &
         (b% job% change_initial_m2 .and. .not. restart)) then
         if (b% m(2) /= b% job% new_m2) then
            call set_m2(binary_id, b% job% new_m2, ierr)
            if (ierr /= 0) then
               write(*, *) "error in set_m2"
               return
            end if
            need_to_fix_generations = .true.
         end if
      end if

      if (b% job% change_period_eccentricity .or. &
         (b% job% change_initial_period_eccentricity .and. .not. restart)) then
         if (b% period /= b% job% new_period * secday .or. &
            b% eccentricity /= b% job% new_eccentricity) then
            call set_period_eccentricity(binary_id, &
               b% job% new_period * secday, b% job% new_eccentricity, ierr)
            if (ierr /= 0) then
               write(*, *) "error in set_period_eccentricity"
               return
            end if
            need_to_fix_generations = .true.
         end if
      end if

      if (b% job% change_separation_eccentricity .or. &
         (b% job% change_initial_separation_eccentricity .and. .not. restart)) then
         if (b% separation /= b% job% new_separation * Rsun .or. &
            b% eccentricity /= b% job% new_eccentricity) then
            call set_separation_eccentricity(binary_id, &
               b% job% new_separation * Rsun, b% job% new_eccentricity, ierr)
            if (ierr /= 0) then
               write(*, *) "error in set_separation_eccentricity"
               return
            end if
            need_to_fix_generations = .true.
         end if
      end if

      if (need_to_fix_generations) then
         if (b% point_mass_i /= 1) then
            b% s1% generations = 0
         end if
         if (b% point_mass_i /= 2) then
            b% s2% generations = 0
         end if
         b% generations = 0

         ! this updates 'old' values in case there is a retry on the first step
         ierr = binary_finish_step(b)

         if (ierr == keep_going) then
            ierr = 0
         else
            ierr = -1
         end if
      end if

   end subroutine

end module run_binary_support
