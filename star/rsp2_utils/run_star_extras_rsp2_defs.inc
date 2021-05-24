      
      ! GYRE "best" info
      real(dp) :: best_period, best_cycles_to_double
      integer :: best_model_number, best_order

      ! summary info at time of recently completely period
      integer :: num_periods, run_num_steps_end_prev, &
         run_num_iters_end_prev, run_num_retries_end_prev
      real(dp) :: period, KE_growth, KE_growth_avg, prev_KE_max, &
         delta_R_growth, delta_R_growth_avg, prev_delta_R, &
         period_max_vsurf_div_cs, period_delta_R, period_delta_Teff, &
         period_delta_logL, period_delta_Mag
      ! info for period in progress
      real(dp) :: time_started, v_div_cs_max, &
         KE_min, KE_max, R_min, R_max, L_min, L_max, T_min, T_max
