! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

! Photo restart test: at the start of every step, fill every user variable of
! each star with values that depend on the star and on their previous values,
! so a variable that a restart does not restore makes the runs diverge:
!   - s% xtra, s% ixtra, s% lxtra                  (saved in photos by MESA)
!   - s% extra_work, s% extra_iwork                (saved in photos by MESA)
!   - pt_hook_int, pt_hook_real (module variables) (saved via other_photo_write/read)
!   - s% x_ctrl, x_integer_ctrl, x_logical_ctrl, x_character_ctrl
!     (controls: re-read from the inlist on restart, never in photos; so they are
!      set from saved state only, and entries set in the inlist are left alone)

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use binary_def
      use utils_lib, only: mesa_error

      implicit none

      integer, parameter :: len_iwork = 10, len_work = 10
      integer, parameter :: max_ids = 10

      ! x_*ctrl entries not set in the inlist (true = this file may fill them)
      logical :: free_x_ctrl(num_x_ctrls, max_ids), free_x_integer_ctrl(num_x_ctrls, max_ids), &
         free_x_logical_ctrl(num_x_ctrls, max_ids), free_x_character_ctrl(num_x_ctrls, max_ids)

      ! user data saved only through other_photo_write / other_photo_read
      integer :: pt_hook_int(max_ids) = 0
      real(dp) :: pt_hook_real(max_ids) = 0d0

      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (id > max_ids) call mesa_error(__FILE__, __LINE__, 'star id > max_ids')

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% other_photo_write => photo_write
         s% other_photo_read => photo_read

         ! the inlist has been read: remember which x_*ctrl it left at default
         free_x_ctrl(:, id) = (s% x_ctrl == 0d0)
         free_x_integer_ctrl(:, id) = (s% x_integer_ctrl == 0)
         free_x_logical_ctrl(:, id) = (.not. s% x_logical_ctrl)
         free_x_character_ctrl(:, id) = (len_trim(s% x_character_ctrl) == 0)
      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: i
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (restart) return  ! everything below was restored from the photo
         call star_alloc_extras(id, len_iwork, len_work, ierr)
         if (ierr /= 0) return
         do i = 1, len_iwork
            s% extra_iwork(i) = 10*id + i
         end do
         do i = 1, len_work
            s% extra_work(i) = 0.1d0*id + 0.01d0*i
         end do
         s% xtra(:) = 0d0
         s% ixtra(:) = 0
         s% lxtra(:) = .false.
         pt_hook_int(id) = 0
         pt_hook_real(id) = 1d0/3d0
      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr, i
         type (star_info), pointer :: s
         real(dp) :: m
         character(len=32) :: str
         ierr = 0
         extras_start_step = terminate
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         m = s% star_mass

         ! saved in photos: depend on their previous values
         do i = 1, size(s% xtra)
            s% xtra(i) = 0.9d0*s% xtra(i) + i + 1d-3*m
            s% ixtra(i) = s% ixtra(i) + i
            s% lxtra(i) = .not. s% lxtra(i)
         end do
         do i = 1, s% len_extra_work
            s% extra_work(i) = 0.5d0*s% extra_work(i) + sqrt(dble(i)) + 1d-3*m
         end do
         do i = 1, s% len_extra_iwork
            s% extra_iwork(i) = s% extra_iwork(i) + 2*i
         end do
         pt_hook_int(id) = pt_hook_int(id) + 1
         pt_hook_real(id) = 0.7d0*pt_hook_real(id) + 1d-2*m

         ! controls: from saved state only
         do i = 1, num_x_ctrls
            if (free_x_ctrl(i, id)) s% x_ctrl(i) = s% model_number + 1d-3*i + 1d-6*m
            if (free_x_integer_ctrl(i, id)) s% x_integer_ctrl(i) = s% model_number + i
            if (free_x_logical_ctrl(i, id)) s% x_logical_ctrl(i) = (mod(s% model_number + i, 2) == 0)
            if (free_x_character_ctrl(i, id)) then
               write(str, '(a, i0, a, i0)') 'x_char_', i, '_model_', s% model_number
               s% x_character_ctrl(i) = trim(str)
            end if
         end do

         extras_start_step = keep_going
      end function extras_start_step


      subroutine photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) pt_hook_int(id), pt_hook_real(id)
      end subroutine photo_write


      subroutine photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         read(iounit, iostat=ierr) pt_hook_int(id), pt_hook_real(id)
      end subroutine photo_read


      integer function extras_check_model(id)
         integer, intent(in) :: id
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine data_for_extra_profile_columns


      integer function extras_finish_step(id)
         integer, intent(in) :: id
         extras_finish_step = keep_going
      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine extras_after_evolve

      end module run_star_extras
