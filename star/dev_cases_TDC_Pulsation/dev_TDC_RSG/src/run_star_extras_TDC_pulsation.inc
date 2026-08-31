

      subroutine TDC_pulsation_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) num_periods, run_num_steps_end_prev, &
            run_num_iters_end_prev, run_num_retries_end_prev, &
            period, KE_growth, KE_growth_avg, prev_KE_max, &
            delta_R_growth, delta_R_growth_avg, prev_delta_R, &
            period_max_v_div_vesc, period_max_v_div_cs, period_delta_R, &
            period_delta_Rphot, period_delta_V, period_delta_Vphot, &
            period_delta_Teff, period_delta_logTeff, &
            period_delta_logL, period_delta_Mag, &
            time_started, v_div_cs_max, v_div_vesc_max, &
            KE_min, KE_max, R_min, R_max, Rphot_min, Rphot_max, &
            V_min, V_max, Vphot_min, Vphot_max, L_min, L_max, T_min, T_max
      end subroutine TDC_pulsation_photo_write


      subroutine TDC_pulsation_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         read(iounit, iostat=ierr) num_periods, run_num_steps_end_prev, &
            run_num_iters_end_prev, run_num_retries_end_prev, &
            period, KE_growth, KE_growth_avg, prev_KE_max, &
            delta_R_growth, delta_R_growth_avg, prev_delta_R, &
            period_max_v_div_vesc, period_max_v_div_cs, period_delta_R, &
            period_delta_Rphot, period_delta_V, period_delta_Vphot, &
            period_delta_Teff, period_delta_logTeff, &
            period_delta_logL, period_delta_Mag, &
            time_started, v_div_cs_max, v_div_vesc_max, &
            KE_min, KE_max, R_min, R_max, Rphot_min, Rphot_max, &
            V_min, V_max, Vphot_min, Vphot_max, L_min, L_max, T_min, T_max
      end subroutine TDC_pulsation_photo_read


      subroutine TDC_pulsation_extras_startup(id, restart, ierr)
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
            delta_R_growth = 0
            delta_R_growth_avg = 0
            prev_delta_R = 0
            period_max_v_div_cs = 0
            period_max_v_div_vesc = 0
            period_delta_R = 0
            period_delta_Rphot = 0
            period_delta_V = 0
            period_delta_Vphot = 0
            period_delta_Teff = 0
            period_delta_logTeff = 0
            period_delta_logL = 0
            period_delta_Mag = 0
            time_started = 0
            v_div_cs_max = 0
            v_div_vesc_max = 0
            KE_min = 0
            KE_max = 0
            R_min = 0
            R_max = 0
            Rphot_min = 0
            Rphot_max = 0
            V_min = 0
            V_max = 0
            Vphot_min = 0
            Vphot_max = 0
            L_min = 0
            L_max = 0
            T_min = 0
            T_max = 0
         end if
      end subroutine TDC_pulsation_extras_startup


      integer function TDC_pulsation_extras_finish_step(id)
         use utils_lib, only: is_bad
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr, test_period
         real(dp) :: target_period
         logical :: doing_pulses
         include 'formats'

         TDC_pulsation_extras_finish_step = terminate
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         TDC_pulsation_extras_finish_step = keep_going

         doing_pulses = s% x_logical_ctrl(7)
         if (.not. doing_pulses) return
         if (period > 0d0 .and. num_periods > 1) then
            target_period = period/(24d0*3600d0)
         else if (s% x_ctrl(7) > 0d0) then
            target_period = s% x_ctrl(7)
         else
            target_period = 0d0
            if (.not. is_bad(s% dynamic_timescale) .and. s% dynamic_timescale > 0d0) &
               target_period = s% dynamic_timescale/(24d0*3600d0)
         end if
         if (target_period <= 0d0 .or. &
               .not. (s% v_flag .or. s% u_flag)) return
         if (s% star_LNA_model_number >= 0 .and. &
               s% model_number <= s% star_LNA_model_number) return
         if (.not. get_period_info()) return

         test_period = s% x_integer_ctrl(7)
         if (num_periods < test_period .or. test_period <= 0) return

         ! have finished test run
         call report_test_results
         TDC_pulsation_extras_finish_step = terminate

         contains

         logical function get_period_info()
            real(dp) :: v_surf, v_surf_start, KE, min_period, time_ended, &
               Rphot, v_phot, delta_R, min_deltaR_for_periods, growth_avg_frac_new, &
               min_period_div_target, cs, vesc
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
               v_surf = s% u(1)
               v_surf_start = s% u_start(1)
            else
               v_surf = 0d0
               v_surf_start = 0d0
            end if
            Rphot = s% photosphere_r*Rsun
            v_phot = s% photosphere_v
            if (Rphot < Rphot_min) Rphot_min = Rphot
            if (Rphot > Rphot_max) Rphot_max = Rphot
            if (v_surf < V_min) V_min = v_surf
            if (v_surf > V_max) V_max = v_surf
            if (v_phot < Vphot_min) Vphot_min = v_phot
            if (v_phot > Vphot_max) Vphot_max = v_phot
            cs = s% csound(1)
            if (v_surf > v_div_cs_max*cs) v_div_cs_max = v_surf/cs
            vesc = sqrt(2*s% cgrav(1)*s% m(1)/(s% r(1)))
            if (v_surf > v_div_vesc_max*vesc) v_div_vesc_max = v_surf/vesc

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
               growth_avg_frac_new = s% x_ctrl(9)
               KE_growth = (KE_max - prev_KE_max)/prev_KE_max
               KE_growth_avg = growth_avg_frac_new*KE_growth + &
                  (1d0 - growth_avg_frac_new)*KE_growth_avg
               delta_R_growth = (delta_R - prev_delta_R)/prev_delta_R
               delta_R_growth_avg = growth_avg_frac_new*delta_R_growth + &
                  (1d0 - growth_avg_frac_new)*delta_R_growth_avg
            end if

            period_delta_Teff = T_max - T_min
            period_delta_logTeff = log10(T_max/T_min)
            period_delta_R = R_max - R_min
            period_delta_Rphot = Rphot_max - Rphot_min
            period_delta_V = V_max - V_min
            period_delta_Vphot = Vphot_max - Vphot_min
            period_delta_logL = log10(L_max/L_min)
            period_delta_Mag = 2.5d0*period_delta_logL
            period_max_v_div_cs = v_div_cs_max
            period_max_v_div_vesc = v_div_vesc_max
            prev_KE_max = KE_max
            prev_delta_R = period_delta_R
            !                 1       2       3         4       5        6       7          8        9
            write(*,'(i4,a14,i6,a13,f8.3,a13,f9.3,a9,f9.4,a15,f10.4,a13,f9.4,a13,f10.4,a11,f9.4,a13,f9.4)')  &
               num_periods, &
               'steps/cycle', s% model_number - run_num_steps_end_prev, &              ! 1   a14,i6
               'iters/step',  &
                  dble(s% total_num_solver_iterations - run_num_iters_end_prev)/ &
                  dble(s% model_number - run_num_steps_end_prev), &                    ! 2   a13,f8.3
               'period (d)', period/(24*3600), &                                       ! 3   a13,f9.3
               'growth', delta_R_growth_avg, &                                         ! 4   a9,f9.4
               'delta R/Rsun', period_delta_R/Rsun, &                                  ! 5   a15,f10.4
               'delta logL', period_delta_logL, &                                      ! 6   a13,f9.4
               'delta Teff', period_delta_Teff, &                                      ! 7   a13,f10.4
               'max v/cs', period_max_v_div_cs, &                                      ! 8   a11,f9.4
               'max v/vesc', period_max_v_div_vesc                                     ! 9   a13,f9.4

            time_started = time_ended
            run_num_steps_end_prev = s% model_number
            run_num_iters_end_prev = s% total_num_solver_iterations
            run_num_retries_end_prev = s% num_retries
            call init_min_max_info
            get_period_info = .true.

         end function get_period_info

         subroutine init_min_max_info
            v_div_cs_max = 0d0
            v_div_vesc_max = 0d0
            KE_min = 1d99
            KE_max  = -1d99
            R_min = 1d99
            R_max = -1d99
            Rphot_min = 1d99
            Rphot_max = -1d99
            V_min = 1d99
            V_max = -1d99
            Vphot_min = 1d99
            Vphot_max = -1d99
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

      end function TDC_pulsation_extras_finish_step


      subroutine TDC_pulsation_extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         ierr = 0
      end subroutine TDC_pulsation_extras_after_evolve


      integer function TDC_pulsation_how_many_extra_history_columns(id)
         integer, intent(in) :: id

         TDC_pulsation_how_many_extra_history_columns = 15
      end function TDC_pulsation_how_many_extra_history_columns


      subroutine TDC_pulsation_data_for_extra_history_columns(id, n, names, vals, ierr)
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
         names(i) = 'growth'; vals(i) = delta_R_growth_avg; i=i+1
         names(i) = 'max_v_div_cs'; vals(i) = period_max_v_div_cs; i=i+1
         names(i) = 'max_v_div_vesc'; vals(i) = period_max_v_div_vesc; i=i+1
         names(i) = 'delta_R'; vals(i) = period_delta_R/Rsun; i=i+1
         names(i) = 'delta_rphot'; vals(i) = period_delta_Rphot/Rsun; i=i+1
         names(i) = 'delta_V'; vals(i) = period_delta_V/1d5; i=i+1
         names(i) = 'delta_Vphot'; vals(i) = period_delta_Vphot/1d5; i=i+1
         names(i) = 'delta_Teff'; vals(i) = period_delta_Teff; i=i+1
         names(i) = 'delta_logTeff'; vals(i) = period_delta_logTeff; i=i+1
         names(i) = 'delta_logL'; vals(i) = period_delta_logL; i=i+1
         names(i) = 'delta_Mag'; vals(i) = period_delta_Mag; i=i+1
         names(i) = 'KE_growth'; vals(i) = KE_growth; i = i+1
         names(i) = 'KE_growth_avg'; vals(i) = KE_growth_avg; i = i+1

      end subroutine TDC_pulsation_data_for_extra_history_columns


      integer function TDC_pulsation_how_many_extra_profile_columns(id)
         integer, intent(in) :: id

         TDC_pulsation_how_many_extra_profile_columns = 0
      end function TDC_pulsation_how_many_extra_profile_columns


      subroutine TDC_pulsation_data_for_extra_profile_columns( &
            id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character(len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz, n)
         integer, intent(out) :: ierr

         ierr = 0
      end subroutine TDC_pulsation_data_for_extra_profile_columns
