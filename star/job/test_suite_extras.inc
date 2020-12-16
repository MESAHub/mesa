subroutine test_suite_testhub_start(s, ierr)
   use star_def, only: star_info
   type (star_info), pointer :: s
   integer, intent(out) :: ierr

   integer :: iounit

   ierr = 0

   open(newunit=iounit, file='testhub.yml', position='append')
   close(iounit)

end subroutine test_suite_testhub_start


subroutine test_suite_testhub_end(s, ierr)
   use star_def, only: star_info
   type (star_info), pointer :: s
   integer, intent(out) :: ierr

   integer :: iounit
   real(dp) :: val

   ierr = 0

   ! write quantities that should exist for every inlist
   open(newunit=iounit, file='testhub.yml', position='append')

   ! timing
   write(iounit,'(6X, A, F6.2)') 'runtime_minutes: ', test_suite_runtime


   ! performance counters
   write(iounit,'(6X, A, I6)') 'steps: ', s% total_steps_finished
   write(iounit,'(6X, A, I6)') 'retries: ', s% total_step_retries
   write(iounit,'(6X, A, I6)') 'redos: ', s% total_step_redos
   write(iounit,'(6X, A, I6)') 'solver_calls_made: ', s% total_num_solver_calls_made
   write(iounit,'(6X, A, I6)') 'solver_calls_failed: ', (s% total_num_solver_calls_made - s% total_num_solver_calls_converged)
   write(iounit,'(6X, A, I6)') 'solver_iterations: ', s% total_num_solver_iterations


   ! log_rel_run_E_err
   val = -99
   if (s% total_energy_end /= 0d0) &
        val = safe_log10(abs(s% cumulative_energy_error/s% total_energy_end))
   write(iounit,'(6X, A, F26.16)') 'log_rel_run_E_err: ', val

   ! user-defined quantities
   call write_testhub_extras(iounit)

   close(iounit)

end subroutine test_suite_testhub_end


subroutine write_testhub_extras(iounit)

   integer, intent(in) :: iounit
   integer :: i

   do i = 1, NUM_TESTHUB_EXTRAS
      if (testhub_extras_names(i) == '') exit
      if (i == 1) write(iounit,'(6X, A)') 'extra_testhub_names:'
      write(iounit,'(12X, A)') "- '" // trim(testhub_extras_names(i)) // "'"
   end do

   do i = 1, NUM_TESTHUB_EXTRAS
      if (testhub_extras_names(i) == '') exit
      if (i == 1) write(iounit,'(6X, A)') 'extra_testhub_vals:'
      write(iounit,'(12X, A, ES26.16)') "- ", testhub_extras_vals(i)
   end do

end subroutine write_testhub_extras


subroutine test_suite_startup(s, restart, ierr)
   use star_def, only: star_info
   type (star_info), pointer :: s
   logical, intent(in) :: restart
   integer, intent(out) :: ierr
   ierr = 0
   call system_clock(test_suite_time0,test_suite_clock_rate)
   test_suite_restart_flag = restart
   if (.not. test_suite_restart_flag) call test_suite_testhub_start(s, ierr)
   testhub_extras_names(1:NUM_TESTHUB_EXTRAS) = ''
   testhub_extras_vals(1:NUM_TESTHUB_EXTRAS) = 0
end subroutine test_suite_startup


subroutine test_suite_after_evolve(s, ierr)
   use star_def, only: star_info
   type (star_info), pointer :: s
   integer, intent(out) :: ierr
   ierr = 0
   if (.not. test_suite_restart_flag) then
      call system_clock(test_suite_time1,test_suite_clock_rate)
      test_suite_runtime = real(test_suite_time1 - test_suite_time0, dp) / test_suite_clock_rate / 60
      write(*,'(/,a50,f12.2,99i10/)') 'runtime (minutes), retries, steps', &
           test_suite_runtime, s% num_retries, s% model_number
   end if
   if (.not. test_suite_restart_flag) call test_suite_testhub_end(s, ierr)
   ierr = 0
end subroutine test_suite_after_evolve
