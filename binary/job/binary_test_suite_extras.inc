subroutine test_suite_testhub_start(b, ierr)
   use binary_def, only : binary_info
   type (binary_info), pointer :: b
   integer, intent(out) :: ierr
   
   integer :: iounit
   
   ierr = 0
   
   open(newunit = iounit, file = 'testhub.yml', position = 'append')
   close(iounit)

end subroutine test_suite_testhub_start


subroutine test_suite_testhub_end(b, ierr)
   use binary_def, only : binary_info
   type (binary_info), pointer :: b
   integer, intent(out) :: ierr
   
   integer :: iounit
   integer :: steps, retries, redos, solver_calls_made, solver_calls_failed, solver_iterations
   real(dp) :: val1, val2
   
   ierr = 0
   
   ! write quantities that should exist for every inlist
   open(newunit = iounit, file = 'testhub.yml', position = 'append')
   
   ! timing
   write(iounit, '(6X, A, F6.2)') 'runtime_minutes: ', test_suite_runtime
   
   ! calculate log_rel_run_E_err
   val1 = -99
   if (b% s1% total_energy_end /= 0d0) &
         val1 = safe_log10(abs(b% s1% cumulative_energy_error / b% s1% total_energy_end))
   
   val2 = -99
   if (b% evolve_both_stars) then
      if (b% s2% total_energy_end /= 0d0) &
            val2 = safe_log10(abs(b% s2% cumulative_energy_error / b% s2% total_energy_end))
   end if
   
   ! model information
   ! these quantities can be cumulative as they are stored in saved models
   write(iounit, '(6X, A, I12)') 'model_number: ', b% s1% model_number
   write(iounit, '(6X, A, ES26.16)') 'star_age: ', b% s1% star_age
   write(iounit, '(6X, A, I12)') 'num_retries: ', b% s1% num_retries
   write(iounit, '(6X, A, F26.16)') 'log_rel_run_E_err: ', max(val1, val2)
   
   ! per-part performance counters
   if (b% evolve_both_stars) then
      steps = max(b% s1% total_steps_finished, b% s2% total_steps_finished) ! should be equal
      retries = b% s1% total_step_retries + b% s2% total_step_retries ! total from both stars; does not count mdot
      redos = b% s1% total_step_redos + b% s2% total_step_redos ! total from both stars; does not count mdot
      solver_calls_made = b% s1% total_num_solver_calls_made + b% s2% total_num_solver_calls_made  ! total from both stars
      solver_calls_failed = (b% s1% total_num_solver_calls_made - b% s1% total_num_solver_calls_converged) + &
            (b% s2% total_num_solver_calls_made - b% s2% total_num_solver_calls_converged) ! total from both stars
      solver_iterations = b% s1% total_num_solver_iterations + b% s2% total_num_solver_iterations ! total from both stars
   else
      steps = b% s1% total_steps_finished
      retries = b% s1% total_step_retries
      redos = b% s1% total_step_redos
      solver_calls_made = b% s1% total_num_solver_calls_made
      solver_calls_failed = (b% s1% total_num_solver_calls_made - b% s1% total_num_solver_calls_converged)
      solver_iterations = b% s1% total_num_solver_iterations
   end if
   
   write(iounit, '(6X, A, I12)') 'steps: ', steps
   write(iounit, '(6X, A, I12)') 'retries: ', retries
   write(iounit, '(6X, A, I12)') 'redos: ', redos
   write(iounit, '(6X, A, I12)') 'solver_calls_made: ', solver_calls_made
   write(iounit, '(6X, A, I12)') 'solver_calls_failed: ', solver_calls_failed
   write(iounit, '(6X, A, I12)') 'solver_iterations: ', solver_iterations
   
   ! user-defined quantities
   call write_testhub_extras(iounit)
   
   close(iounit)

end subroutine test_suite_testhub_end


subroutine write_testhub_extras(iounit)
   
   integer, intent(in) :: iounit
   integer :: i
   
   do i = 1, NUM_TESTHUB_EXTRAS
      if (testhub_extras_names(i) == '') exit
      if (i == 1) write(iounit, '(6X, A)') 'extra_testhub_names:'
      write(iounit, '(12X, A)') "- '" // trim(testhub_extras_names(i)) // "'"
   end do
   
   do i = 1, NUM_TESTHUB_EXTRAS
      if (testhub_extras_names(i) == '') exit
      if (i == 1) write(iounit, '(6X, A)') 'extra_testhub_vals:'
      write(iounit, '(12X, A, ES26.16)') "- ", testhub_extras_vals(i)
   end do

end subroutine write_testhub_extras


subroutine test_suite_startup(b, restart, ierr)
   use binary_def, only : binary_info
   type (binary_info), pointer :: b
   logical, intent(in) :: restart
   integer, intent(out) :: ierr
   ierr = 0
   call system_clock(test_suite_time0, test_suite_clock_rate)
   test_suite_restart_flag = restart
   if (.not. test_suite_restart_flag) call test_suite_testhub_start(b, ierr)
   testhub_extras_names(1:NUM_TESTHUB_EXTRAS) = ''
   testhub_extras_vals(1:NUM_TESTHUB_EXTRAS) = 0
end subroutine test_suite_startup


subroutine test_suite_after_evolve(b, ierr)
   use binary_def, only : binary_info
   type (binary_info), pointer :: b
   integer, intent(out) :: ierr
   ierr = 0
   if (.not. test_suite_restart_flag) then
      call system_clock(test_suite_time1, test_suite_clock_rate)
      test_suite_runtime = real(test_suite_time1 - test_suite_time0, dp) / test_suite_clock_rate / 60
      write(*, '(/,a50,f12.2,99i10/)') 'runtime (minutes), retries, steps', &
            test_suite_runtime, b% s1% num_retries, b% s1% model_number
   end if
   if (.not. test_suite_restart_flag) call test_suite_testhub_end(b, ierr)
   ierr = 0
end subroutine test_suite_after_evolve
