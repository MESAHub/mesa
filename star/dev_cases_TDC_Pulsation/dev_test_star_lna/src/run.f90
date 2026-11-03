program run
   use run_star_support, only: do_read_star_job, start_run1_star
   use run_star_extras, only: extras_controls
   use star_lib, only: star_ptr, star_export_pulse_data, free_star, starlib_shutdown
   use star_def, only: star_info
   use const_def, only: dp
   use utils_lib, only: mesa_error

   implicit none

   integer :: id, ierr
   logical :: restart, analyzed
   real(dp) :: saved_dt
   type(star_info), pointer :: s

   id = 0
   restart = .false.
   analyzed = .false.
   call do_read_star_job('inlist', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'failed to read inlist')

   ! Load the model and run LNA without entering the evolution loop.
   call start_run1_star( &
      do_alloc_star=.true., do_free_star=.true., okay_to_restart=.false., &
      id=id, restart=restart, restart_filename='restart_photo', &
      pgstar_ok=.false., dbg=.false., extras_controls=extras_controls, &
      ierr=ierr, inlist_fname_arg='inlist', stop_after_star_LNA=analyzed)
   if (ierr /= 0 .or. .not. analyzed) &
      call mesa_error(__FILE__, __LINE__, 'startup LNA failed')

   call star_ptr(id, s, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'failed to get star pointer')
   ! The .mod file supplies no previous thermal state for eps_grav.
   saved_dt = s% dt
   s% dt = 0d0
   call star_export_pulse_data(id, 'GYRE', 'gyre.data', &
      s% add_center_point_to_pulse_data, s% keep_surface_point_for_pulse_data, &
      s% add_atmosphere_to_pulse_data, ierr)
   s% dt = saved_dt
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'GYRE export failed')

   call free_star(id, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'failed to free star')
   call starlib_shutdown
   write(*,'(a)') 'star_LNA analysis and GYRE export complete; no evolution steps taken.'
end program run
