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
      
      real(dp) :: Psurf, Tsurf, Lsurf

! here are the x controls used below

!alpha_mlt_routine
         !alpha_H = s% x_ctrl(21)
         !alpha_other = s% x_ctrl(22)
         !H_limit = s% x_ctrl(23)

      
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
         s% other_alpha_mlt => alpha_mlt_routine       
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
            !write(*,2) 'alpha_mlt', k, s% alpha_mlt(k), 
         end do
         !stop
      end subroutine alpha_mlt_routine
      
      
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
            Psurf = s% P(1)
            Tsurf = s% T(1)
            Lsurf = s% L(1)
            call alloc_extra_info(s)
            if (s% x_logical_ctrl(1)) then
               if (s% fixed_L_for_BB_outer_BC < 0d0) then
                  s% fixed_L_for_BB_outer_BC = Lsurf
               else
                  Lsurf = s% fixed_L_for_BB_outer_BC
               end if
               call switch_BCs(s)
            end if
         else ! it is a restart
            call unpack_extra_info(s)
            if (s% x_logical_ctrl(1)) call switch_BCs(s)
         end if
         if (s% x_ctrl(16) > 0d0) &
            s% x_ctrl(2) = s% star_mass - s% x_ctrl(16)
         
      end subroutine extras_startup
      
      
      subroutine switch_BCs(s)
         type (star_info), pointer :: s
         !s% use_compression_outer_BC = .true.
         !s% use_momentum_outer_BC = .true.
         s% use_T_black_body_outer_BC = .true.
         s% use_fixed_L_for_BB_outer_BC = .true.
         s% fixed_L_for_BB_outer_BC = Lsurf
      end subroutine switch_BCs
      
      
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
         nz = s% nz
         write(*,*)
         do k=nz-1,1,-1
            if (s% ye(k) >= 0.495d0) then
               m = find0( &
                  s% m(k+1), s% ye(k+1) - 0.495d0, &
                  s% m(k), s% ye(k) - 0.495d0)
               if (m > 0d0) write(*,2) 'ye = 0.495', k, m/Msun
               exit
            end if
         end do
         call test_suite_after_evolve(s, ierr)
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
         integer :: k, k0, k1
         real(dp) :: v_esc
         real(dp),pointer, dimension(:) :: vel
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = keep_going        
         if (.not. s% x_logical_ctrl(2)) return 
         if (s% u_flag) then
            vel => s% u
         else if (s% v_flag) then
            vel => s% v
         else
            write(*,*) 'extras_start_step: Must have either v_flag or u_flag enabled'
            stop
         end if
         k0 = s% nz+1
         do k = 1, s% nz
            if (s% q(k) < s% x_ctrl(19)) then
               !write(*,2) 'nothing in outer layer with v >= v_esc', k, s% q(k)
               exit ! only check outer layer
            end if
            v_esc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
            if (vel(k) > v_esc) then
               k0 = k
               exit
            end if
         end do
         if (k0 >= s% nz) return
         ! k0 is outermost location below surface with u > escape velocity.
         ! remove inward from k0 where u large enough compared to escape velocity.
         k1 = k0
         do k = k0+1, s% nz
            v_esc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
            if (vel(k) < s% x_ctrl(17)*v_esc) then ! stop removing here
               k1 = k-1
               exit
            end if
         end do
         if (s% q(k1) > s% x_ctrl(18)) then
            write(*,2) 'v > vesc, but too little to bother with', k1, s% q(k1)
            return ! too little to bother with
         end if
         k1 = max(k1, s% x_integer_ctrl(20))
         do while (s% L(k1) <= 0)
            k1 = k1 + 1
         end do
         call star_remove_surface_at_cell_k(s% id, k1, ierr)
         if (ierr /= 0) then
            write(*,*) 'extras_start_step failed in star_remove_surface_at_cell_k'
            write(*,2) 'at q', k1, s% q(k1)
            extras_start_step = terminate
            return
         end if
         if (s% x_logical_ctrl(1)) then
            Psurf = s% P(1)
            Tsurf = s% T(1)
            Lsurf = s% L(1)
            call switch_BCs(s)
         end if
         extras_start_step = keep_going
      end function extras_start_step


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
         if (s% x_logical_ctrl(1)) call switch_BCs(s)
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
         !call move_flg(in_period_with_vsurf_gt_cs)
         num_ints = i
         
         i = 0
         call move_dbl(Psurf)   
         call move_dbl(Tsurf)   
         call move_dbl(Lsurf)   
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
      
