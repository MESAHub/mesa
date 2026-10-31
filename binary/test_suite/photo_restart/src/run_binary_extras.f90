! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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

! Photo restart test: at the start of every binary step, fill every user
! variable of the binary with values that depend on the binary and on their
! previous values, so a variable that a restart does not restore makes the runs
! diverge:
!   - b% xtra, b% ixtra, b% lxtra                  (saved in photos by MESA)
!   - pt_hook_int, pt_hook_real (module variables) (saved via other_binary_photo_write/read)
!   - b% x_ctrl, x_integer_ctrl, x_logical_ctrl, x_character_ctrl
!     (controls: re-read from the inlist on restart, never in photos; so they are
!      set from saved state only, and entries set in the inlist are left alone)

      module run_binary_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use binary_def

      implicit none

      include "binary_test_suite_extras_def.inc"

      ! x_*ctrl entries not set in the inlist (true = this file may fill them)
      logical :: free_x_ctrl(binary_num_x_ctrls), free_x_integer_ctrl(binary_num_x_ctrls), &
         free_x_logical_ctrl(binary_num_x_ctrls), free_x_character_ctrl(binary_num_x_ctrls)

      ! user data saved only through other_binary_photo_write / other_binary_photo_read
      integer :: pt_hook_int = 0
      real(dp) :: pt_hook_real = 0d0

      contains

      include "binary_test_suite_extras.inc"

      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve
         ! set here: the binary photo is read right after this routine
         b% other_binary_photo_write => photo_write
         b% other_binary_photo_read => photo_read

         b% warn_binary_extra =.false.

         ! the inlist has been read: remember which x_*ctrl it left at default
         free_x_ctrl = (b% x_ctrl == 0d0)
         free_x_integer_ctrl = (b% x_integer_ctrl == 0)
         free_x_logical_ctrl = (.not. b% x_logical_ctrl)
         free_x_character_ctrl = (len_trim(b% x_character_ctrl) == 0)
      end subroutine extras_binary_controls


      subroutine photo_write(binary_id, iounit)
         integer, intent(in) :: binary_id, iounit
         write(iounit) pt_hook_int, pt_hook_real
      end subroutine photo_write


      subroutine photo_read(binary_id, iounit, ierr)
         integer, intent(in) :: binary_id, iounit
         integer, intent(out) :: ierr
         read(iounit, iostat=ierr) pt_hook_int, pt_hook_real
      end subroutine photo_read


      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items


      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_columns


      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         extras_binary_startup = terminate
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if

         call test_suite_startup(b, restart, ierr)

         if (.not. restart) then  ! on restart these come from the photo
            b% xtra(:) = 0d0
            b% ixtra(:) = 0
            b% lxtra(:) = .false.
            pt_hook_int = 0
            pt_hook_real = 1d0/7d0
         end if

         extras_binary_startup = keep_going
      end function extras_binary_startup


      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer :: i
         character(len=32) :: str

         extras_binary_start_step = terminate
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if

         ! saved in photos: depend on their previous values
         do i = 1, size(b% xtra)
            b% xtra(i) = 0.9d0*b% xtra(i) + i + 1d-3*b% period/secday
            b% ixtra(i) = b% ixtra(i) + i
            b% lxtra(i) = .not. b% lxtra(i)
         end do
         pt_hook_int = pt_hook_int + 1
         pt_hook_real = 0.7d0*pt_hook_real + 1d-2*b% separation/rsun

         ! controls: from saved state only
         do i = 1, binary_num_x_ctrls
            if (free_x_ctrl(i)) b% x_ctrl(i) = b% model_number + 1d-3*i + 1d-6*b% period/secday
            if (free_x_integer_ctrl(i)) b% x_integer_ctrl(i) = b% model_number + i
            if (free_x_logical_ctrl(i)) b% x_logical_ctrl(i) = (mod(b% model_number + i, 2) == 0)
            if (free_x_character_ctrl(i)) then
               write(str, '(a, i0, a, i0)') 'x_char_', i, '_model_', b% model_number
               b% x_character_ctrl(i) = trim(str)
            end if
         end do

         extras_binary_start_step = keep_going
      end function  extras_binary_start_step


      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if
         extras_binary_check_model = keep_going
      end function extras_binary_check_model


      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if
         extras_binary_finish_step = keep_going
      end function extras_binary_finish_step


      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then  ! failure in  binary_ptr
            return
         end if
         call test_suite_after_evolve(b, ierr)
      end subroutine extras_binary_after_evolve

      end module run_binary_extras
