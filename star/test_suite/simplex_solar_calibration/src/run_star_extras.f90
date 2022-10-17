! ***********************************************************************
!
!   Copyright (C) 2020 The MESA Team
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
   use auto_diff
   use utils_lib
   
   implicit none
   
   include "test_suite_extras_def.inc"
   
   ! you can add your own data declarations here.

contains

include "test_suite_extras.inc"
   
   
   subroutine extras_controls(id, ierr)
      use simplex_search_data
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      s% extras_startup => extras_startup
      s% how_many_extra_history_columns => how_many_extra_history_columns
      s% data_for_extra_history_columns => data_for_extra_history_columns
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns
      star_simplex_procs% set_my_vars => set_my_vars
      star_simplex_procs% will_set_my_param => will_set_my_param
      star_simplex_procs% extras_check_model => extras_check_model
      star_simplex_procs% extras_finish_step => extras_finish_step
      star_simplex_procs% extras_after_evolve => extras_after_evolve
      
      s% how_many_other_mesh_fcns => how_many_mesh_fcns
      s% other_mesh_fcn_data => gradr_grada_mesh_fcn_data
   end subroutine extras_controls
   
   
   subroutine set_my_vars(id, ierr) ! called from star_simplex code
      use simplex_search_data, only : include_my_var1_in_chi2, my_var1
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ! my_var's are predefined in the simplex_search_data.
      ! this routine's job is to assign those variables to current value in the model.
      ! it is called whenever a new value of chi2 is calculated.
      ! only necessary to set the my_var's you are actually using.
      ierr = 0
      if (include_my_var1_in_chi2) then
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         my_var1 = s% Teff
      end if
   end subroutine set_my_vars
   
   
   subroutine will_set_my_param(id, i, new_value, ierr) ! called from star_simplex code
      use simplex_search_data, only : vary_my_param1
      integer, intent(in) :: id
      integer, intent(in) :: i ! which of my_param's will be set
      real(dp), intent(in) :: new_value
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ! old value has not yet been changed.
      ! do whatever is necessary for this new value.
      ! i.e. change whatever mesa params you need to adjust.
      ! as example, my_param1 is alpha_mlt
      ierr = 0
      if (i == 1 .and. vary_my_param1) then
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% mixing_length_alpha = new_value
      end if
   end subroutine will_set_my_param
   
   
   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call test_suite_startup(s, restart, ierr)
   end subroutine extras_startup
   
   
   ! returns either keep_going, retry, or terminate.
   integer function extras_check_model(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going
      ! by default, indicate where (in the code) MESA terminated
      if (extras_check_model == terminate) s% termination_code = t_extras_check_model
   end function extras_check_model
   
   
   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      how_many_extra_history_columns = 0
   end function how_many_extra_history_columns
   
   
   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len = maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      
      !note: do NOT add these names to history_columns.list
      ! the history_columns.list is only for the built-in log column options.
      ! it must not include the new column names you are adding here.
      
      ierr = 0
   end subroutine data_for_extra_history_columns
   
   
   integer function how_many_extra_profile_columns(id)
      integer, intent(in) :: id
      how_many_extra_profile_columns = 0
   end function how_many_extra_profile_columns
   
   
   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr
      integer :: k
      ierr = 0
      
      !note: do NOT add these names to profile_columns.list
      ! the profile_columns.list is only for the built-in profile column options.
      ! it must not include the new column names you are adding here.
      
      ! here is an example for adding a profile column
      !if (n /= 1) call mesa_error(__FILE__,__LINE__,'data_for_extra_profile_columns')
      !names(1) = 'beta'
      !do k = 1, nz
      !   vals(k,1) = s% Pgas(k)/s% P(k)
      !end do
   end subroutine data_for_extra_profile_columns
   
   
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
      real(dp) :: dt
      type (star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call test_suite_after_evolve(s, ierr)
   end subroutine extras_after_evolve
   
   
   subroutine how_many_mesh_fcns(id, n)
      integer, intent(in) :: id
      integer, intent(out) :: n
      n = 1
   end subroutine how_many_mesh_fcns
   
   
   subroutine gradr_grada_mesh_fcn_data(&
      id, nfcns, names, gval_is_xa_function, vals1, ierr)
      integer, intent(in) :: id
      integer, intent(in) :: nfcns
      character (len = *) :: names(:)
      logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
      real(dp), pointer :: vals1(:) ! =(nz, nfcns)
      integer, intent(out) :: ierr
      
      real(dp), pointer :: vals(:, :)
      type (star_info), pointer :: s
      integer :: nz, k
      real(dp) :: weight, width, center
      real(dp), parameter :: maxval = 700d0 ! max value for tanh from crlibm
      
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      weight = s% x_ctrl(1)
      width = s% x_ctrl(2)
      center = s% x_ctrl(3)
      
      names(1) = 'gradr_grada_function'
      gval_is_xa_function(1) = .false.
      
      nz = s% nz
      vals(1:nz, 1:nfcns) => vals1(1:nz * nfcns)
      
      do k = 1, nz
         vals(k, 1) = weight * tanh(min(maxval, (s% gradr(k) - s% grada(k) - center) / width)) * width
      end do
   
   end subroutine gradr_grada_mesh_fcn_data

end module run_star_extras
