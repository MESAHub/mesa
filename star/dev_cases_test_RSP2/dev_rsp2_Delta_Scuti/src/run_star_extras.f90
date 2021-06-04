! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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
      use gyre_lib
      
      implicit none

      include 'test_suite_extras_def.inc'
      include '../../../rsp2_utils/run_star_extras_rsp2_defs.inc'


!alpha_mlt_routine
         !alpha_H = s% x_ctrl(21)
         !alpha_other = s% x_ctrl(22)
         !H_limit = s% x_ctrl(23)

!gyre
      !x_logical_ctrl(37) = .false. ! if true, then run GYRE
      !x_integer_ctrl(1) = 2 ! output GYRE info at this step interval
      !x_logical_ctrl(1) = .false. ! save GYRE info whenever save profile
      !x_integer_ctrl(2) = 2 ! max number of modes to output per call
      !x_logical_ctrl(2) = .false. ! output eigenfunction files
      !x_integer_ctrl(3) = 0 ! mode l (e.g. 0 for p modes, 1 for g modes)
      !x_integer_ctrl(4) = 1 ! order
      !x_ctrl(1) = 0.158d-05 ! freq ~ this (Hz)
      !x_ctrl(2) = 0.33d+03 ! growth < this (days)
            
      
      contains


      include 'test_suite_extras.inc'
      include '../../../rsp2_utils/run_star_extras_rsp2.inc'
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (ierr /= 0) return                  
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% other_alpha_mlt => alpha_mlt_routine       
         s% other_photo_write => photo_write
         s% other_photo_read => photo_read
      end subroutine extras_controls


      subroutine alpha_mlt_routine(id, ierr)
         use chem_def, only: ih1
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, h1
         real(dp) :: alpha_H, alpha_other, H_limit
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         alpha_H = s% x_ctrl(21)
         alpha_other = s% x_ctrl(22)
         H_limit = s% x_ctrl(23)
         h1 = s% net_iso(ih1)
         !write(*,1) 'alpha_H', alpha_H
         !write(*,1) 'alpha_other', alpha_other
         !write(*,1) 'H_limit', H_limit
         !write(*,2) 'h1', h1
         !write(*,2) 's% nz', s% nz
         if (alpha_H <= 0 .or. alpha_other <= 0 .or. h1 <= 0) return
         do k=1,s% nz
            if (s% xa(h1,k) >= H_limit) then
               s% alpha_mlt(k) = alpha_H
            else
               s% alpha_mlt(k) = alpha_other
            end if
         end do
         !stop
      end subroutine alpha_mlt_routine


      subroutine photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         call rsp2_photo_write(id, iounit)
      end subroutine photo_write


      subroutine photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         call rsp2_photo_read(id, iounit, ierr)
      end subroutine photo_read

      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         if (ierr /= 0) return
         call rsp2_extras_startup(id, restart, ierr)
      end subroutine extras_startup


      ! returns either keep_going, retry, or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         extras_finish_step = rsp2_extras_finish_step(id)
      end function extras_finish_step

      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
         call rsp2_extras_after_evolve(id, ierr)
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model
      
      
      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         how_many_extra_history_columns = rsp2_how_many_extra_history_columns(id)
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         call rsp2_data_for_extra_history_columns(id, n, names, vals, ierr)
      end subroutine data_for_extra_history_columns
      
      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         how_many_extra_profile_columns = rsp2_how_many_extra_profile_columns(id)
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         call rsp2_data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      end subroutine data_for_extra_profile_columns


      end module run_star_extras
      
