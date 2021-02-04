! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
      
      implicit none
      
      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% other_set_pgstar_controls => set_pgstar_controls       
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return         
         call test_suite_startup(s, restart, ierr)         
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if 
         if (s% x_logical_ctrl(7) .and. & ! inlist_finish
             len_trim(s% x_character_ctrl(1)) > 0) then
            call star_read_controls(id, s% x_character_ctrl(1), ierr)
            if (ierr /= 0) then
               write(*,*) 'failed reading ' // trim(s% x_character_ctrl(1))
               return
            end if         
         end if
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         use num_lib, only: find0
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt, m
         integer :: k, nz
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
         if (ierr /= 0) return
         
         if (s% x_logical_ctrl(6)) & ! inlist_prepare
            call check_termination_code( &
               t_log_max_temp_upper_limit, t_non_fe_core_infall_limit, -1)
         
         if (s% x_logical_ctrl(7)) & ! inlist_finish
            call check_termination_code( &
               t_fe_core_infall_limit, t_log_Rsurf_upper_limit, t_bound_mass_min_limit)
         
         contains
         
         subroutine check_termination_code(tc1, tc2, tc3)
            integer, intent(in) :: tc1, tc2, tc3
            if (s% termination_code /= tc1 .and. &
                s% termination_code /= tc2 .and. &
                s% termination_code /= tc3) then
               if (s% termination_code > 0 .and. s% termination_code <= num_termination_codes) then
                  write(*,*) 'Failed to get a valid termination reason. got this instead: ' // &
                     trim(termination_code_str(s% termination_code))
               else
                  write(*,*) 'Failed to get a valid termination reason'
               end if
               ierr = -1
            else
               write(*,*) 'Has a valid termination reason: ' // &
                  trim(termination_code_str(s% termination_code))
            end if
         end subroutine check_termination_code
         
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
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: lgTmax, lgRsurf
         include 'formats'
         extras_start_step = keep_going    
         ierr = 0
         call star_ptr(id, s, ierr)
         if (failed('star_ptr',ierr)) return
         if (.not. s% x_logical_ctrl(7)) return
         ! check logR
         lgRsurf = log10(exp(s% xh(s% i_lnR,1))/Rsun)
         !write(*,2) 'lgRsurf', s% model_number, lgRsurf
         if (lgRsurf > s% x_ctrl(19)) then
            s% mass_change = -2d0
            write(*,*) 'surface logR > limit: turn on wind'
         
         
            !write(*,*) 'surface logR > limit ', lgRsurf, s% x_ctrl(19)
            !write(*,*) 'prune surface to v/v_esc = ', s% x_ctrl(20)
            !call star_remove_surface_by_v_surf_div_v_escape(id, s% x_ctrl(20), ierr)
            !if (failed('star_remove_surface_by_v_surf_div_v_escape',ierr)) return
         else
            s% mass_change = 0d0
         end if
         ! check u_flag vs v_flag choice
         lgTmax = maxval(s% xh(s% i_lnT,1:s% nz))/ln10
         !write(*,2) 'lgTmax', s% model_number, lgTmax
         if (s% u_flag) then
            s% dt_div_min_dr_div_cs_limit = s% x_ctrl(23)
            s% dt_div_min_dr_div_cs_hard_limit = s% x_ctrl(24)
            !write(*,2) 'lower limit', s% model_number, s% x_ctrl(21)
            if (lgTmax < s% x_ctrl(21)) then ! switch to v_flag
               ! do add new before remove old so can set initial values
               write(*,*) 'new_v_flag', .true.
               call star_set_v_flag(s% id, .true., ierr)
               if (failed('star_set_v_flag',ierr)) return
               write(*,*) 'new_u_flag', .false.
               call star_set_u_flag(s% id, .false., ierr)
               if (failed('star_set_u_flag',ierr)) return
               s% use_avQ_art_visc = .true.
            end if
         else if (s% v_flag) then
            s% dt_div_min_dr_div_cs_limit = 1d99
            s% dt_div_min_dr_div_cs_hard_limit = 1d99
            !write(*,2) 'upper limit', s% model_number, s% x_ctrl(22)
            if (lgTmax > s% x_ctrl(22)) then ! switch to u_flag
               ! do add new before remove old so can set initial values
               write(*,*) 'new_u_flag', .true.
               call star_set_u_flag(id, .true., ierr)
               if (failed('star_set_u_flag',ierr)) return
               write(*,*) 'new_v_flag', .false.
               call star_set_v_flag(s% id, .false., ierr)
               if (failed('star_set_v_flag',ierr)) return
               s% use_avQ_art_visc = .false.
            end if
         end if
         
         contains
      
         logical function failed(str,ierr)
            character (len=*), intent(in) :: str
            integer, intent(in) :: ierr
            failed = (ierr /= 0)
            if (.not. failed) return
            write(*,*) 'extras_start_step failed in ' // trim(str)
            extras_start_step = terminate
         end function failed
         
      end function extras_start_step


      subroutine set_pgstar_controls(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% History_Panels2_xaxis_name = 'model_number'
         if (.not. s% x_logical_ctrl(7)) return
         ! check age and set History_Panels2_xaxis_name
         if (s% star_age > 1d0) then
            s% History_Panels2_xaxis_name = 'star_age'
         else if (s% star_age*secyer > 24*60*60) then
            s% History_Panels2_xaxis_name = 'star_age_day'
         else if (s% star_age*secyer > 60*60) then
            s% History_Panels2_xaxis_name = 'star_age_hr'
         else if (s% star_age*secyer > 60) then
            s% History_Panels2_xaxis_name = 'star_age_min'
         else
            s% History_Panels2_xaxis_name = 'star_age_sec'
         end if
      end subroutine set_pgstar_controls
   

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: v_esc
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)         
      end function extras_finish_step

      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg
         !call move_int(vsurf_gt_cs_count)   
         !call move_flg(using_fixed_outer_BCs)
         num_ints = i
         
         i = 0
         !call move_dbl(vsurf)   
         !call move_dbl(Lsurf)   
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            double precision :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info


      end module run_star_extras
      
