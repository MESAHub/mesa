module run_star_extras
   use star_lib, only: star_ptr
   use star_def, only: star_info

   implicit none

contains

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s

      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (.not. s% job% load_saved_model) then
         write(*,'(a)') 'dev_test_star_lna requires load_saved_model.'
         ierr = -1
         return
      end if
      s% star_LNA_flag = .true.
      s% star_LNA_model_number = -1
      s% star_LNA_stop_after_run = .true.
      s% star_LNA_set_initial_velocity = .false.
   end subroutine extras_controls

end module run_star_extras
