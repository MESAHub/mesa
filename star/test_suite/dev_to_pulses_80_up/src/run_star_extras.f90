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
      use gyre_lib
      
      implicit none
      
      include "test_suite_extras_def.inc"

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
         s% other_remove_surface => remove_ejecta_one_cell_per_step
         s% use_other_remove_surface = .true.
      end subroutine extras_controls
      
      
      subroutine remove_ejecta_one_cell_per_step(id, ierr, j)
         integer, intent(in) :: id
         integer, intent(out) :: ierr, j
         ierr = 0
         if (star_ejecta_mass(id) > 0.1d0*Msun) then
            call star_remove_surface_at_cell_k(id, 2, ierr)
            write(*,2) 'remove_ejecta_one_cell_per_step', s% model_number
         end if
      end  subroutine remove_ejecta_one_cell_per_step
      
      
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
         
         if (.not. s% x_logical_ctrl(37)) return
         
         ! Initialize GYRE

         call gyre_init('gyre.in')

         ! Set constants

         call gyre_set_constant('G_GRAVITY', standard_cgrav)
         call gyre_set_constant('C_LIGHT', clight)
         call gyre_set_constant('A_RADIATION', crad)

         call gyre_set_constant('M_SUN', Msun)
         call gyre_set_constant('R_SUN', Rsun)
         call gyre_set_constant('L_SUN', Lsun)

         call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')

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
         if (.not. s% x_logical_ctrl(37)) return
         call gyre_final()
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
         extras_start_step = keep_going    
      end function extras_start_step

  
      include 'gyre_in_mesa_extras_finish_step.inc'
   

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
         if (.not. s% x_logical_ctrl(37)) return
         extras_finish_step = gyre_in_mesa_extras_finish_step(id)
         if (extras_finish_step == terminate) &
             s% termination_code = t_extras_finish_step
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
      
