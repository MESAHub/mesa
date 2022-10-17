subroutine extras_controls(id, ierr)
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   
   ! this is the place to set any procedure pointers you want to change
   ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
   
   
   ! the extras functions in this file will not be called
   ! unless you set their function pointers as done below.
   ! otherwise we use a null_ version which does nothing (except warn).
   
   s% extras_startup => extras_startup
   s% extras_start_step => extras_start_step
   s% extras_check_model => extras_check_model
   s% extras_finish_step => extras_finish_step
   s% extras_after_evolve => extras_after_evolve
   s% how_many_extra_history_columns => how_many_extra_history_columns
   s% data_for_extra_history_columns => data_for_extra_history_columns
   s% how_many_extra_profile_columns => how_many_extra_profile_columns
   s% data_for_extra_profile_columns => data_for_extra_profile_columns
   
   s% how_many_extra_history_header_items => how_many_extra_history_header_items
   s% data_for_extra_history_header_items => data_for_extra_history_header_items
   s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
   s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

end subroutine extras_controls


subroutine extras_startup(id, restart, ierr)
   integer, intent(in) :: id
   logical, intent(in) :: restart
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
end subroutine extras_startup


integer function extras_start_step(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   extras_start_step = 0
end function extras_start_step


! returns either keep_going, retry, or terminate.
integer function extras_check_model(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   extras_check_model = keep_going
   if (.false. .and. s% star_mass_h1 < 0.35d0) then
      ! stop when star hydrogen mass drops to specified level
      extras_check_model = terminate
      write(*, *) 'have reached desired hydrogen mass'
      return
   end if
   
   
   ! if you want to check multiple conditions, it can be useful
   ! to set a different termination code depending on which
   ! condition was triggered.  MESA provides 9 customizeable
   ! termination codes, named t_xtra1 .. t_xtra9.  You can
   ! customize the messages that will be printed upon exit by
   ! setting the corresponding termination_code_str value.
   ! termination_code_str(t_xtra1) = 'my termination condition'
   
   ! by default, indicate where (in the code) MESA terminated
   if (extras_check_model == terminate) s% termination_code = t_extras_check_model
end function extras_check_model


integer function how_many_extra_history_columns(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   how_many_extra_history_columns = 0
end function how_many_extra_history_columns


subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
   integer, intent(in) :: id, n
   character (len = maxlen_history_column_name) :: names(n)
   real(dp) :: vals(n)
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   
   ! note: do NOT add the extras names to history_columns.list
   ! the history_columns.list is only for the built-in history column options.
   ! it must not include the new column names you are adding here.

end subroutine data_for_extra_history_columns


integer function how_many_extra_profile_columns(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   how_many_extra_profile_columns = 0
end function how_many_extra_profile_columns


subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
   integer, intent(in) :: id, n, nz
   character (len = maxlen_profile_column_name) :: names(n)
   real(dp) :: vals(nz, n)
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   integer :: k
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   
   ! note: do NOT add the extra names to profile_columns.list
   ! the profile_columns.list is only for the built-in profile column options.
   ! it must not include the new column names you are adding here.
   
   ! here is an example for adding a profile column
   !if (n /= 1) stop 'data_for_extra_profile_columns'
   !names(1) = 'beta'
   !do k = 1, nz
   !   vals(k,1) = s% Pgas(k)/s% P(k)
   !end do

end subroutine data_for_extra_profile_columns


integer function how_many_extra_history_header_items(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   how_many_extra_history_header_items = 0
end function how_many_extra_history_header_items


subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
   integer, intent(in) :: id, n
   character (len = maxlen_history_column_name) :: names(n)
   real(dp) :: vals(n)
   type(star_info), pointer :: s
   integer, intent(out) :: ierr
   ierr = 0
   call star_ptr(id, s, ierr)
   if(ierr/=0) return
   
   ! here is an example for adding an extra history header item
   ! also set how_many_extra_history_header_items
   ! names(1) = 'mixing_length_alpha'
   ! vals(1) = s% mixing_length_alpha

end subroutine data_for_extra_history_header_items


integer function how_many_extra_profile_header_items(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   how_many_extra_profile_header_items = 0
end function how_many_extra_profile_header_items


subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
   integer, intent(in) :: id, n
   character (len = maxlen_profile_column_name) :: names(n)
   real(dp) :: vals(n)
   type(star_info), pointer :: s
   integer, intent(out) :: ierr
   ierr = 0
   call star_ptr(id, s, ierr)
   if(ierr/=0) return
   
   ! here is an example for adding an extra profile header item
   ! also set how_many_extra_profile_header_items
   ! names(1) = 'mixing_length_alpha'
   ! vals(1) = s% mixing_length_alpha

end subroutine data_for_extra_profile_header_items


! returns either keep_going or terminate.
! note: cannot request retry; extras_check_model can do that.
integer function extras_finish_step(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   extras_finish_step = keep_going
   
   ! to save a profile,
   ! s% need_to_save_profiles_now = .true.
   ! to update the star log,
   ! s% need_to_update_history_now = .true.
   
   ! see extras_check_model for information about custom termination codes
   ! by default, indicate where (in the code) MESA terminated
   if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
end function extras_finish_step


subroutine extras_after_evolve(id, ierr)
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
end subroutine extras_after_evolve
