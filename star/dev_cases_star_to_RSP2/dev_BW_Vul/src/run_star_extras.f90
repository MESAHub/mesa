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
      use gyre_lib
      
      implicit none

      include 'test_suite_extras_def.inc'

      integer :: num_periods, prev_cycle_run_num_steps, &
         run_num_iters_prev_period, run_num_retries_prev_period, best_model_number
      real(dp) :: period, time_started, period_r_min, period_max_vsurf_div_cs, best_G_div_P, &
         prev_period, prev_growth, prev_period_delta_r, prev_period_max_vsurf_div_cs, &
         best_growth, best_period
         

!alpha_mlt_routine
         !alpha_H = s% x_ctrl(21)
         !alpha_other = s% x_ctrl(22)
         !H_limit = s% x_ctrl(23)

!gyre
      !x_logical_ctrl(37) = .false. ! if true, then run GYRE
      !x_integer_ctrl(1) = 2 ! output GYRE info at this step interval
      !x_logical_ctrl(1) = .false. ! save GYRE info whenever save profile
      !x_integer_ctrl(2) = 2 ! max number of modes to output per call
      !x_logical_ctrl(2) = .false. ! output eigenfunction files
      !x_integer_ctrl(3) = 0 ! mode l (e.g. 0 for p modes, 1 for g modes)
      !x_integer_ctrl(4) = 1 ! order
      !x_ctrl(1) = 0.158d-05 ! freq ~ this (Hz)
      !x_ctrl(2) = 0.33d+03 ! growth < this (days)
            
      
      contains


      include 'test_suite_extras.inc'
            
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (ierr /= 0) return                  
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% other_alpha_mlt => alpha_mlt_routine       
         s% other_photo_write => photo_write
         s% other_photo_read => photo_read
      end subroutine extras_controls


      subroutine alpha_mlt_routine(id, ierr)
         use chem_def, only: ih1
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, h1
         real(dp) :: alpha_H, alpha_other, H_limit
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         alpha_H = s% x_ctrl(21)
         alpha_other = s% x_ctrl(22)
         H_limit = s% x_ctrl(23)
         h1 = s% net_iso(ih1)
         !write(*,1) 'alpha_H', alpha_H
         !write(*,1) 'alpha_other', alpha_other
         !write(*,1) 'H_limit', H_limit
         !write(*,2) 'h1', h1
         !write(*,2) 's% nz', s% nz
         if (alpha_H <= 0 .or. alpha_other <= 0 .or. h1 <= 0) return
         do k=1,s% nz
            if (s% xa(h1,k) >= H_limit) then
               s% alpha_mlt(k) = alpha_H
            else
               s% alpha_mlt(k) = alpha_other
            end if
         end do
         !stop
      end subroutine alpha_mlt_routine


      subroutine photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) num_periods, prev_cycle_run_num_steps, &
            run_num_iters_prev_period, run_num_retries_prev_period, &
            period, time_started, period_r_min, period_max_vsurf_div_cs, &
            prev_period, prev_growth, prev_period_delta_r, prev_period_max_vsurf_div_cs, &
            best_model_number, best_G_div_P, best_growth, best_period
      end subroutine photo_write


      subroutine photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         read(iounit, iostat=ierr) num_periods, prev_cycle_run_num_steps, &
            run_num_iters_prev_period, run_num_retries_prev_period, &
            period, time_started, period_r_min, period_max_vsurf_div_cs, &
            prev_period, prev_growth, prev_period_delta_r, prev_period_max_vsurf_div_cs, &
            best_model_number, best_G_div_P, best_growth, best_period
      end subroutine photo_read

      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         if (ierr /= 0) return
         if (.not. restart) then     
            num_periods = 0
            period = 0
            time_started = 0
            prev_cycle_run_num_steps = 0
            run_num_iters_prev_period = 0
            run_num_retries_prev_period = 0
            period_r_min = 0
            period_max_vsurf_div_cs = 0
            prev_growth = 0
            prev_period_delta_r = 0
            prev_period = 0
            prev_period_max_vsurf_div_cs = 0
            best_model_number = 0
            best_G_div_P = 0
            best_growth = 0
            best_period = 0
         end if
         if (.not. s% x_logical_ctrl(5)) then
            call gyre_init('gyre.in')
            call gyre_set_constant('G_GRAVITY', standard_cgrav)
            call gyre_set_constant('C_LIGHT', clight)
            call gyre_set_constant('A_RADIATION', crad)
            call gyre_set_constant('M_SUN', Msun)
            call gyre_set_constant('R_SUN', Rsun)
            call gyre_set_constant('L_SUN', Lsun)
            call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')         
         else
            call gyre_linear_analysis_and_set_velocities(s,restart,ierr)
         end if    
         
      end subroutine extras_startup


      subroutine gyre_linear_analysis_and_set_velocities(s,restart,ierr)
         use const_def
         use math_lib
         use gyre_lib
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         
         real(dp), allocatable     :: global_data(:)
         real(dp), allocatable     :: point_data(:,:)
         integer                   :: ipar(5), mode_l
         real(dp)                  :: rpar(1)
         
         integer, parameter :: modes = 3
         integer :: npts(modes), nz, i, k, i_v
         real(dp), pointer :: vel(:)
         real(dp), allocatable, dimension(:,:) :: r, v
         real(dp) :: v_surf, v1, amix1, amix2, amixF, &
            period(modes)
         
         include 'formats'
         
         if (restart) return
         
         write(*,*) 'set gyre starting velocities'
         
         nz = s% nz
         allocate(r(modes,nz+10), v(modes,nz+10))
         npts = 0

         call gyre_init('gyre.in')

         call gyre_set_constant('G_GRAVITY', standard_cgrav)
         call gyre_set_constant('C_LIGHT', clight)
         call gyre_set_constant('A_RADIATION', crad)

         call gyre_set_constant('M_SUN', Msun)
         call gyre_set_constant('R_SUN', Rsun)
         call gyre_set_constant('L_SUN', Lsun)

         call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')
         
         mode_l = 0 ! mode l (e.g. 0 for p modes, 1 for g modes)
                        ! should match gyre.in mode l
         
         !write(*,*) 'call star_get_pulse_data'
         call star_get_pulse_data(s%id, 'GYRE', &
            .FALSE., .FALSE., .FALSE., global_data, point_data, ierr)
         if (ierr /= 0) then
            print *,'Failed when calling get_pulse_data'
            return
         end if
         
         !write(*,*) 'call star_write_pulse_data'
         call star_write_pulse_data(s%id, &
            'GYRE', 'gyre.data', global_data, point_data, ierr)
         if (ierr /= 0) return

         !write(*,*) 'call gyre_set_model'
         call gyre_set_model(global_data, point_data, s%gyre_data_schema)

         write(*, 100) 'order', 'freq (Hz)', 'P (sec)', &
           'P (min)', 'P (day)', 'growth (day)', 'growth/P'
100      format(A8,A16,A16,A14,A12,A16,A14)

         rpar(1) = 0.5d-6 ! freq < this (Hz)
         ipar(1) = s% model_number
         ipar(2) = 1 ! order_target
         ipar(3) = 1 ! 1 means output eigenfunction files
         ipar(4) = 3 ! max number of modes to output per call
         ipar(5) = 0 ! num_written

         !write(*,*) 'call gyre_get_modes'
         call gyre_get_modes(mode_l, process_mode_, ipar, rpar)
         
         amix1 = s% x_ctrl(4) ! s% RSP_fraction_1st_overtone
         amix2 = s% x_ctrl(5) ! s% RSP_fraction_2nd_overtone
         if((amix1+amix2) > 1d0) then
            write(*,*) 'AMIX DO NOT ADD UP RIGHT' 
            call mesa_error(__FILE__,__LINE__,'set_gyre_linear_analysis')
         end if
         amixF = 1d0 - (amix1 + amix2)
         
         if (amixF > 0d0 .and. npts(1) /= nz-1) then
            write(*,3) 'amixF > 0d0 .and. npts(1) /= nz-1', npts(1)
            write(*,*) 'cannot use fundamental for setting starting velocities'
            write(*,*) 'need to add code to interpolate from gyre grid to model'
            call mesa_error(__FILE__,__LINE__,'set_gyre_linear_analysis')
            ierr = -1
            return
         end if
         
         if (AMIX1 > 0d0 .and. npts(2) /= nz-1) then
            write(*,3) 'AMIX1 > 0d0 .and. npts(2) /= nz-1', npts(2)
            write(*,*) 'cannot use 1st overtone for setting starting velocities'
            write(*,*) 'need to add code to interpolate from gyre grid to model'
            call mesa_error(__FILE__,__LINE__,'set_gyre_linear_analysis')
            ierr = -1
            return
         end if
         
         if (AMIX2 > 0d0 .and. npts(2) /= nz-1) then
            write(*,3) 'AMIX2 > 0d0 .and. npts(3) /= nz-1', npts(3)
            write(*,*) 'cannot use 2nd overtone for setting starting velocities'
            write(*,*) 'need to add code to interpolate from gyre grid to model'
            call mesa_error(__FILE__,__LINE__,'set_gyre_linear_analysis')
            ierr = -1
            return
         end if
         
         v_surf = amixF*v(1,nz-1) + AMIX1*v(2,nz-1) + AMIX2*v(3,nz-1)   
         v1 = 1d5/v_surf
         if (s% x_ctrl(6) > 0d0) v1 = v1*s% x_ctrl(6)
         
         if (s% v_flag) then
            vel => s% v
            i_v = s% i_v
         else if (s% u_flag) then
            vel => s% u
            i_v = s% i_u
         else
            call mesa_error(__FILE__,__LINE__,'set_gyre_linear_analysis vel')
         end if
         
         do i=nz-1,1,-1
            k = nz+1-i ! v(1) from gyre => vel(nz) in star
            vel(k) = v1*(amixF*v(1,i) + AMIX1*v(2,i) + AMIX2*v(3,i))
            write(*,2) 'vel', k, vel(k)
         end do
         vel(1) = vel(2)
         s% v_center = 0d0
         
         do k=1,nz
            s% xh(i_v,k) = vel(k)
         end do
         
         write(*,'(A)')
         write(*,1) 'v_surf F 1 2', v_surf, v(1,nz-1), v(2,nz-1), v(3,nz-1)
         write(*,1) 'amixF amix1 amix2', amixF, amix1, amix2
         write(*,'(A)')
         write(*,2) 'nz', nz
         write(*,1) 'v(1)/1d5', vel(1)/1d5    
         write(*,1) 'T(nz)', s% T(s%nz)             
         write(*,1) 'L_center/Lsun', s% L_center/Lsun           
         write(*,1) 'R_center/Rsun', s% R_center/Rsun           
         write(*,1) 'M_center/Msun', s% M_center/Msun           
         write(*,1) 'L(1)/Lsun', s% L(1)/Lsun           
         write(*,1) 'R(1)/Rsun', s% r(1)/Rsun           
         write(*,1) 'M(1)/Msun', s% m(1)/Msun           
         write(*,1) 'X(1)', s% X(1)      
         write(*,1) 'Y(1)', s% Y(1)      
         write(*,1) 'Z(1)', s% Z(1)      
         write(*,1) 'tau_factor', s% tau_factor   
         write(*,1) 'tau_base', s% tau_base   
         write(*,1) 'Teff', s% Teff
         write(*,*) 

         contains

         subroutine process_mode_ (md, ipar, rpar, retcode)

            type(mode_t), intent(in) :: md
            integer, intent(inout)   :: ipar(:)
            real(dp), intent(inout)  :: rpar(:)
            integer, intent(out)     :: retcode

            character(LEN=128) :: filename
            integer               :: unit, k, model_number, order_target, num_written, max_to_write
            complex(dp)           :: cfreq
            real(dp)              :: freq, growth, freq_lim
            logical               :: write_flag
            type(grid_t)          :: gr

            max_to_write = ipar(4)
            num_written = ipar(5)
            if (num_written >= max_to_write) return
            ipar(5) = num_written + 1

            model_number = ipar(1)
            order_target = ipar(2)
            freq_lim = rpar(1)
            write_flag = (ipar(3) == 1)

            cfreq = md% freq('HZ')
            freq = REAL(cfreq)

            if (AIMAG(cfreq) > 0._dp) then ! unstable
               growth = 1d0/(2*pi*24*3600*AIMAG(cfreq))
               write(*, 100)  md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600), &
                  growth, growth*freq*24*3600
100               format(I8,E16.4,F16.4,F14.4,F12.4,E16.4,E14.4)
            else ! stable
               write(*, 110) md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600), 'stable'
110         format(I8,E16.4,F16.4,F14.4,F12.4,A16)
            end if
            
            if (md%n_pg > modes) return

            gr = md%grid()

            period(md%n_pg) = 1d0/freq
            npts(md%n_pg) = md%n_k
            do k = 1, md%n_k
               r(md%n_pg,k) = gr%pt(k)%x
               v(md%n_pg,k) = md%xi_r(k)            
            end do

            if (write_flag) then
               ! Write the mode radial & horizontal eigenfunctions, together with the differential work
               write(filename, 120) 'eigfunc.', md%n_pg, '.dat'
120            format(A,I0,A)
               !print *,'Writing eigenfunction to ', TRIM(filename)
               !write(*,*)
               open(NEWUNIT=unit, FILE=filename, STATUS='REPLACE')
               write(unit, 130) 'x=r/R', 'Real(xi_r/R)', 'Imag(xi_r/R)', 'Real(xi_h/R)', 'Imag(xi_h/R)', 'dW/dx'
130               format(6(1X,A24))
               do k = 1, md%n_k
                  write(unit, 140) gr%pt(k)%x, md%xi_r(k), md%xi_h(k), md%dW_dx(k)
140               format(6(1X,E24.16))
               end do
               close(unit)
            end if

            retcode = 0

         end subroutine process_mode_
         
      end subroutine gyre_linear_analysis_and_set_velocities

      
      include 'gyre_in_mesa_extras_finish_step.inc'
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_finish_step(id)
         use chem_def
         integer, intent(in) :: id
         integer :: ierr, i
         type (star_info), pointer :: s
         real(dp) :: rel_run_E_err, target_period, time_ended, period_r_max, &
            v_surf, v_surf_start, growth, period_delta_r
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         if (s% x_integer_ctrl(1) > 0) then
             if (MOD(s% model_number, s% x_integer_ctrl(1)) == 0) then
               extras_finish_step = gyre_in_mesa_extras_finish_step(id)
               if (s% ixtra3_array(1) > 0) then ! get the GYRE results
                  do i=1,s% ixtra3_array(1)
                     if (s% xtra1_array(i) == 0d0 .or. s% ixtra1_array(i) /= s% x_integer_ctrl(4)) cycle
                     if (best_G_div_P == 0d0 .or. s% xtra1_array(i)/s% xtra2_array(i) < best_G_div_P) then
                        best_G_div_P = s% xtra1_array(i)/s% xtra2_array(i)
                        best_growth = s% xtra1_array(i)
                        best_period = s% xtra2_array(i)
                        best_model_number = s% model_number
                     end if
                  end do
                  if (best_period > 0) &
                     write(*,*) 'best_model_number best_G_div_P period(d) growth(d)', &
                        best_model_number, best_G_div_P, best_period, best_growth
               end if
            end if
         end if
         if (extras_finish_step == terminate) then
            s% termination_code = t_extras_finish_step
            return
         end if
         if (extras_finish_step /= keep_going) then
            return
         end if
         if (.not. s% x_logical_ctrl(7)) return
         if (s% x_ctrl(7) <= 0d0) return ! must give expected period > 0
         if (s% v_flag) then
            v_surf = s% v(1)
            v_surf_start = s% v_start(1)
         else if (s% u_flag) then
            v_surf = s% u_face_val(1)
            v_surf_start = s% u_face_start(1)
         else
            call mesa_error(__FILE__,__LINE__,'extras_finish_step: both v_flag and u_flag are false')
         end if
         if (v_surf/s% csound(1) > period_max_vsurf_div_cs) &
            period_max_vsurf_div_cs = v_surf/s% csound(1)
         ! check_cycle_completed when v_surf goes from positive to negative
         if (v_surf*v_surf_start < 0d0 .and. v_surf > 0d0) then
            period_r_min = s% r(1)
            return
         end if
         ! at max radius when v_surf goes from positive to negative
         if (v_surf*v_surf_start > 0d0 .or. v_surf > 0d0) return
         period_r_max = s% r(1)
         ! either start of 1st cycle, or end of current
         if (time_started == 0) then
            time_started = s% time
            prev_cycle_run_num_steps = s% model_number
            run_num_iters_prev_period = s% total_num_solver_iterations
            run_num_retries_prev_period = s% num_retries
            write(*,*) 'first maximum radius, period calculations start at model, day', &
               s% model_number, s% time/(24*3600)
            return
         end if
         time_ended = s% time
         if (abs(v_surf - v_surf_start) > 1d-10) & ! tweak the end time to match when v_surf == 0
            time_ended = s% time - v_surf*s% dt/(v_surf - v_surf_start)
         period = time_ended - time_started
         if (period/(24*3600) < 0.1d0*s% x_ctrl(7)) return ! reject as bogus if < 10% expected
         num_periods = num_periods + 1
         
         period_delta_r = period_r_max - period_r_min
         if (period_delta_r > 0d0 .and. prev_period_delta_r > 0d0) then
            growth = period/log(period_delta_r/prev_period_delta_r) ! seconds
         else
            growth = 0d0
         end if
         prev_growth = 0.1d0*growth + 0.9d0*prev_growth ! smoothing
         write(*,'(i4,a12,f9.4,a12,e13.4,a14,f9.4,a14,f9.4,a9,i3,a7,i6,a16,f8.3,a6,i7,a9,f10.3)')  &
            num_periods,'period (d)',  period/(24*3600), &
            'growth (d)', prev_growth/(24*3600), &
            'delta R/Rsun', period_delta_r/Rsun, &
            'max vsurf/cs', period_max_vsurf_div_cs, &
            'retries', s% num_retries - run_num_retries_prev_period,     &
            'steps', s% model_number - prev_cycle_run_num_steps, &
            'avg iters/step',  &
               dble(s% total_num_solver_iterations - run_num_iters_prev_period)/ &
               dble(s% model_number - prev_cycle_run_num_steps), &
            'step', s% model_number, 'age (d)', s% time/(24*3600)
         time_started = time_ended
         prev_cycle_run_num_steps = s% model_number
         prev_period = period
         prev_period_delta_r = period_delta_r
         prev_period_max_vsurf_div_cs = period_max_vsurf_div_cs
         run_num_iters_prev_period = s% total_num_solver_iterations
         run_num_retries_prev_period = s% num_retries
         period_max_vsurf_div_cs = 0d0
         if (num_periods < s% x_integer_ctrl(7) .or. s% x_integer_ctrl(7) <= 0) return
         write(*,'(A)')
         write(*,'(A)')
         write(*,'(A)')
         target_period = s% x_ctrl(7)
         rel_run_E_err = s% cumulative_energy_error/s% total_energy
         write(*,*) 'rel_run_E_err', rel_run_E_err
         if (s% total_energy /= 0d0 .and. abs(rel_run_E_err) > 1d-5) then
            write(*,*) '*** BAD rel_run_E_error ***', &
            s% cumulative_energy_error/s% total_energy
         else if (abs(period/(24*3600) - target_period) > 1d-2) then
            write(*,*) '*** BAD period ***', period/(24*3600) - target_period, &
               period/(24*3600), target_period
         else
            write(*,*) 'good match for period', &
               period/(24*3600), target_period
         end if
         write(*,'(A)')
         write(*,'(A)')
         write(*,'(A)')
         extras_finish_step = terminate
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
         call gyre_final()
      end subroutine extras_after_evolve
      

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
         how_many_extra_history_columns = 5
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s, s_other
         integer :: id_other
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = 'num_periods'
         vals(1) = num_periods
         names(2) = 'period'
         vals(2) = prev_period/(24*3600)
         names(3) = 'growth'
         vals(3) = prev_growth/(24*3600)
         names(4) = 'delta_r'
         vals(4) = prev_period_delta_r/Rsun
         names(5) = 'max_vsurf_div_cs'
         vals(5) = prev_period_max_vsurf_div_cs
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0 ! 6
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
         return
         
         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'xtra1'
         names(2) = 'xtra2'
         names(3) = 'xtra3'
         names(4) = 'xtra4'
         names(5) = 'xtra5'
         names(6) = 'xtra6'

         do k=1,nz            
            vals(k,1) = s% xtra1_array(k)
            vals(k,2) = s% xtra2_array(k)
            vals(k,3) = s% xtra3_array(k)
            vals(k,4) = s% xtra4_array(k)
            vals(k,5) = s% xtra5_array(k)
            vals(k,6) = s% xtra6_array(k)            
         end do
            
      end subroutine data_for_extra_profile_columns
      

      end module run_star_extras
      
