

      subroutine rsp2_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) num_periods, run_num_steps_end_prev, &
            run_num_iters_end_prev, run_num_retries_end_prev, &
            period, KE_growth, KE_growth_avg, prev_KE_max, &
            period_max_vsurf_div_cs, period_delta_R, period_delta_Teff, &
            period_delta_logL, period_delta_Mag, time_started, v_div_cs_max, &
            KE_min, KE_max, R_min, R_max, L_min, L_max, T_min, T_max, &
            best_period, best_model_number, best_order, best_cycles_to_double
      end subroutine rsp2_photo_write


      subroutine rsp2_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         read(iounit, iostat=ierr) num_periods, run_num_steps_end_prev, &
            run_num_iters_end_prev, run_num_retries_end_prev, &
            period, KE_growth, KE_growth_avg, prev_KE_max, &
            period_max_vsurf_div_cs, period_delta_R, period_delta_Teff, &
            period_delta_logL, period_delta_Mag, time_started, v_div_cs_max, &
            KE_min, KE_max, R_min, R_max, L_min, L_max, T_min, T_max, &
            best_period, best_model_number, best_order, best_cycles_to_double
      end subroutine rsp2_photo_read


      subroutine rsp2_extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. restart) then     
            num_periods = 0
            run_num_steps_end_prev = 0
            run_num_iters_end_prev = 0
            run_num_retries_end_prev = 0
            period = 0
            KE_growth = 0
            KE_growth_avg = 0
            prev_KE_max = 0
            period_max_vsurf_div_cs = 0
            period_delta_R = 0
            period_delta_Teff = 0
            period_delta_logL = 0
            period_delta_Mag = 0
            time_started = 0
            v_div_cs_max = 0
            KE_min = 0
            KE_max = 0
            R_min = 0
            R_max = 0
            L_min = 0
            L_max = 0
            T_min = 0
            T_max = 0
            best_period = 0
            best_model_number = 0    
            best_order = 0
            best_cycles_to_double = 0                    
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
      end subroutine rsp2_extras_startup
      
      
      integer function rsp2_extras_finish_step(id)
         use chem_def
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr, gyre_interval, test_period
         real(dp) :: target_period
         logical :: doing_pulses
         include 'formats'
         
         rsp2_extras_finish_step = terminate
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         rsp2_extras_finish_step = keep_going
         
         gyre_interval = s% x_integer_ctrl(1)
         if (gyre_interval > 0) then
            if (MOD(s% model_number, gyre_interval) == 0) &
               call get_gyre_info_for_this_step
            if (rsp2_extras_finish_step == terminate) &
               s% termination_code = t_extras_finish_step
            if (rsp2_extras_finish_step /= keep_going) return
         end if
         
         doing_pulses = s% x_logical_ctrl(7)
         if (.not. doing_pulses) return
         target_period = s% x_ctrl(7)
         if (target_period <= 0d0) return
         if (.not. get_period_info()) return
         
         test_period = s% x_integer_ctrl(7)
         if (num_periods < test_period .or. test_period <= 0) return
         
         ! have finished test run
         call report_test_results
         rsp2_extras_finish_step = terminate
         
         contains
      
         subroutine get_gyre_info_for_this_step
            integer :: i
            rsp2_extras_finish_step = gyre_in_mesa_extras_finish_step(id)
            if (s% ixtra3_array(1) > 0) then ! unpack the GYRE results
               do i=1,s% ixtra3_array(1)
                  if (s% xtra1_array(i) == 0d0 .or. &
                     (s% ixtra1_array(i) /= s% x_integer_ctrl(4) .and. s% x_integer_ctrl(4) > 0)) cycle
                  if (s% xtra3_array(i) > 0d0 .and. &
                      (best_cycles_to_double == 0d0 .or. s% xtra3_array(i) < best_cycles_to_double)) then
                     !best_growth = s% xtra1_array(i)
                     best_period = 1d0/s% xtra2_array(i) ! xtra2_array = freq (s^-1)
                     best_period = best_period/(24*3600) ! change to days
                     best_cycles_to_double = s% xtra3_array(i)
                     best_order = s% ixtra1_array(i)
                     best_model_number = s% model_number
                  end if
               end do
               if (.false. .and. best_cycles_to_double > 0) &
                  write(*,'(a,i6,a,i2,a,f12.0,a,f7.4,a,f9.4)') &
                     'best cycles to double: model_number', best_model_number, &
                     'order', best_order, 'P(s)', 24*3600*best_period, &
                     'P(d)', best_period, 'cycles to double', best_cycles_to_double
            end if
         end subroutine get_gyre_info_for_this_step
         
         logical function get_period_info()
            real(dp) :: v_surf, v_surf_start, KE, KE_avg, min_period, time_ended, &
               delta_R, min_deltaR_for_periods, KE_growth_avg_frac_new, &
               min_period_div_target, cs
            include 'formats'
            get_period_info = .false.
         
            if (s% r(1) < R_min) R_min = s% r(1)
            if (s% r(1) > R_max) R_max = s% r(1)
            if (s% L(1) < L_min) L_min = s% L(1)
            if (s% L(1) > L_max) L_max = s% L(1)
            if (s% Teff < T_min) T_min = s% Teff
            if (s% Teff > T_max) T_max = s% Teff
            KE = s% total_radial_kinetic_energy_end
            if (KE > KE_max) KE_max = KE
            if (KE < KE_min) KE_min = KE
            
            if (s% v_flag) then
               v_surf = s% v(1)
               v_surf_start = s% v_start(1)
            else if (s% u_flag) then
               v_surf = s% u_face_val(1)
               v_surf_start = s% u_face_start(1)
            else
               stop 'rsp2_extras_finish_step: both v_flag and u_flag are false'
            end if
            cs = s% csound(1)
            if (v_surf > v_div_cs_max*cs) v_div_cs_max = v_surf/cs
               
            ! period is completed when v_surf goes from positive to negative during step
            if (v_surf > 0d0 .or. v_surf_start < 0d0) return
            
            if (time_started == 0) then ! start of 1st cycle
               time_started = s% time
               run_num_steps_end_prev = s% model_number
               run_num_iters_end_prev = s% total_num_solver_iterations
               run_num_retries_end_prev = s% num_retries
               prev_KE_max = 0d0
               call init_min_max_info
               write(*,*) 'first maximum radius, period calculations starting at model, day', &
                  s% model_number, s% time/(24*3600)
               return
            end if
            
            delta_R = R_max - R_min
            min_deltaR_for_periods = s% x_ctrl(8)*Rsun
            if (min_deltaR_for_periods > 0d0) then
               if (delta_R < min_deltaR_for_periods) return ! filter out glitches
            end if
         
            time_ended = s% time
            if (abs(v_surf - v_surf_start) > 1d-10) & ! tweak the end time to match when v_surf == 0
               time_ended = s% time - v_surf*s% dt/(v_surf - v_surf_start)
            min_period_div_target = s% x_ctrl(10)
            min_period = target_period*(24*3600)*min_period_div_target
            if (min_period > 0d0 .and. &
                time_ended - time_started < min_period) return ! filter out glitches

            period = time_ended - time_started
            num_periods = num_periods + 1
         
            if (num_periods > 1) then
               KE_avg = 0.5d0*(KE_max + prev_KE_max)
               KE_growth = (KE_max - prev_KE_max)/KE_avg
               KE_growth_avg_frac_new = s% x_ctrl(9)
               KE_growth_avg = KE_growth_avg_frac_new*KE_growth + &
                  (1d0 - KE_growth_avg_frac_new)*KE_growth_avg
            end if
         
            period_delta_Teff = T_max - T_min
            period_delta_R = R_max - R_min
            period_delta_logL = log10(L_max/L_min)
            period_delta_Mag = 2.5d0*period_delta_logL
            period_max_vsurf_div_cs = v_div_cs_max
            prev_KE_max = KE_max
            write(*,'(i4,a14,i6,a14,f8.3,a14,f10.3,99(a14,f10.5))')  &
               num_periods, &
               'steps/cycle', s% model_number - run_num_steps_end_prev, &
               'iters/step',  &
                  dble(s% total_num_solver_iterations - run_num_iters_end_prev)/ &
                  dble(s% model_number - run_num_steps_end_prev), &
               'period (d)', period/(24*3600), &
               'KE growth', KE_growth_avg, &
               'delta R/Rsun', period_delta_R/Rsun, &
               'delta logL', period_delta_logL, &
               'delta Teff', period_delta_Teff, &
               'max vsurf/cs', period_max_vsurf_div_cs

            time_started = time_ended
            run_num_steps_end_prev = s% model_number
            run_num_iters_end_prev = s% total_num_solver_iterations
            run_num_retries_end_prev = s% num_retries
            call init_min_max_info
            get_period_info = .true.

         end function get_period_info
         
         subroutine init_min_max_info
            v_div_cs_max = 0d0
            KE_min = 1d99
            KE_max  = -1d99
            R_min = 1d99
            R_max = -1d99
            L_min = 1d99
            L_max = -1d99
            T_min = 1d99
            T_max = -1d99
         end subroutine init_min_max_info
         
         subroutine report_test_results
            real(dp) :: rel_run_E_err
            write(*,*)
            write(*,*)
            write(*,*)
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
            write(*,*)
            write(*,*)
            write(*,*)
         end subroutine report_test_results
           
      end function rsp2_extras_finish_step
           
      
      include 'gyre_in_mesa_extras_finish_step.inc'


      subroutine rsp2_extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
         call gyre_final()
      end subroutine rsp2_extras_after_evolve


      integer function rsp2_how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         rsp2_how_many_extra_history_columns = 8
      end function rsp2_how_many_extra_history_columns
      
      
      subroutine rsp2_data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: i
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         i = 1
         names(i) = 'num_periods'; vals(i) = num_periods; i=i+1
         names(i) = 'period'; vals(i) = period/(24*3600); i=i+1
         names(i) = 'growth'; vals(i) = KE_growth_avg; i=i+1
         names(i) = 'max_v_div_cs'; vals(i) = period_max_vsurf_div_cs; i=i+1
         names(i) = 'delta_R'; vals(i) = period_delta_R/Rsun; i=i+1
         names(i) = 'delta_Teff'; vals(i) = period_delta_Teff; i=i+1
         names(i) = 'delta_logL'; vals(i) = period_delta_logL; i=i+1
         names(i) = 'delta_Mag'; vals(i) = period_delta_Mag; i=i+1
      end subroutine rsp2_data_for_extra_history_columns


      integer function rsp2_how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         rsp2_how_many_extra_profile_columns = 0 ! 6
      end function rsp2_how_many_extra_profile_columns
      
      
      subroutine rsp2_data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
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
            
      end subroutine rsp2_data_for_extra_profile_columns


      subroutine gyre_linear_analysis_and_set_velocities(s,restart,ierr)
         use const_def
         use math_lib
         use gyre_lib
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interpolate_vector_pm
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         
         real(dp), allocatable     :: global_data(:)
         real(dp), allocatable     :: point_data(:,:)
         integer                   :: ipar(5), mode_l
         real(dp)                  :: rpar(1)
         
         integer, parameter :: modes = 3
         integer :: npts(modes), nz, nz_max, i, k, i_v, j, nz_gyre, mode
         real(dp), pointer :: vel(:)
         real(dp), allocatable, dimension(:,:) :: r, v
         real(dp) :: v_surf, v1, period(modes)
         real(dp), allocatable, dimension(:) :: &
            r_old, r_new, v_old, v_new
         real(dp), pointer :: work1(:) ! =(nz_max, pm_work_size)
         
         include 'formats'
         
         if (restart) return
         
         write(*,*) 'set gyre starting velocities'
         
         nz = s% nz
         nz_max = nz+100
         allocate( &
            r(modes,nz_max), v(modes,nz_max), &
            r_old(nz_max), r_new(nz_max), v_old(nz_max), v_new(nz_max))
         npts(:) = 0
         
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
            write(*,*) 'Failed when calling star_get_pulse_data'
            return
         end if
         
         !write(*,*) 'call star_write_pulse_data'
         call star_write_pulse_data(s%id, &
            'GYRE', 'gyre.data', global_data, point_data, ierr)
         if (ierr /= 0) then
            write(*,*) 'Failed when calling star_write_pulse_data'
            return
         end if

         !write(*,*) 'call gyre_set_model'
         call gyre_set_model(global_data, point_data, 101)

          write(*, 100) 'order', 'freq (Hz)', &
             'P (sec)', 'P (day)', 'cycles to double'
100       format(A8,8A20)

         rpar(1) = 0.5d-6 ! freq < this (Hz)
         ipar(1) = s% model_number
         ipar(2) = 1 ! order_target
         ipar(3) = 1 ! 1 means output eigenfunction files
         ipar(4) = 3 ! max number of modes to output per call
         ipar(5) = 0 ! num_written

         call gyre_get_modes(mode_l, process_mode_, ipar, rpar)
         
         if (s% x_ctrl(4) == 1d0) then
            mode = 2
         else if (s% x_ctrl(5) == 1d0) then
            mode = 3
         else
            mode = 1
         end if
         nz_gyre = npts(mode)
         
         v_surf = v(mode,nz_gyre)
         v1 = 1d5/v_surf ! rescale v's
         if (s% x_ctrl(6) > 0d0) v1 = v1*s% x_ctrl(6)

         do i=1,nz_gyre ! change to star index order
            k = nz_gyre-i+1
            v_old(k) = v1*v(mode,i)
            r_old(k) = r(mode,i)*s% r(1)
         end do
         
         do k=1,nz
            r_new(k) = s% r(k)
         end do
         
         allocate(work1(nz_max*pm_work_size))
         call interpolate_vector_pm( &
            nz_gyre, r_old, nz, r_new, v_old, v_new, work1, 'get_pulse_data', ierr)
         deallocate(work1)
         if (ierr /= 0) then
            write(*,*) 'Failed when calling interpolate_vector_pm for pulse data'
            return
         end if
                  
         if (s% v_flag) then
            do k=1,nz
               s% xh(s% i_v,k) = v_new(k)
            end do
         else if (s% u_flag) then
            do k=1,nz
               s% xh(s% i_u,k) = v_new(k)
            end do
         else
            stop 'set_gyre_linear_analysis failed to set velocities'
         end if
         
         return
         
         do k=1,nz_gyre
            write(*,2) 'gyre r v', k, r_old(k)/Rsun, v_old(k)
            if (k <= nz) write(*,2) 'star r v', k, r_new(k)/Rsun, v_new(k)
            write(*,*)
         end do
         
         do k=nz_gyre+1,nz
            write(*,2) 'star r v', k, r_new(k), v_new(k)
         end do
         
         stop 'gyre_linear_analysis_and_set_velocities'

         contains

         subroutine process_mode_ (md, ipar, rpar, retcode)

            type(mode_t), intent(in) :: md
            integer, intent(inout)   :: ipar(:)
            real(dp), intent(inout)  :: rpar(:)
            integer, intent(out)     :: retcode

            character(LEN=128) :: filename
            integer               :: unit, k, model_number, order_target, num_written, max_to_write
            complex(dp)           :: cfreq
            real(dp)              :: freq, growth, freq_lim, per
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
            growth = AIMAG(cfreq)
            if (freq /= 0d0) then
               per = 1d0/freq
            else
               per = 0d0
            end if

            if (growth > 0._dp) then ! unstable
               write(*, 100) md%n_pg, &
                  freq, per, per/(24*3600), freq/(2d0*pi*growth)
100            format(I8,E20.4,3F20.4)
            else ! stable
               write(*, 110) md%n_pg, &
                  freq, per, per/(24*3600), 'stable'
110            format(I8,E20.4,2F20.4,A20)
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
