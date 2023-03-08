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
      use auto_diff
      use gyre_lib
      
      implicit none

      ! (Gamma1 - 4/3) at the center when integral_gamma1-4/3 first drops below 0
      real(dp) :: gamma1_cntr_pulse_start 

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
         s% other_remove_surface => remove_ejecta_one_cell_per_step
         !s% use_other_remove_surface = .true.

         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write

      end subroutine extras_controls
      
      
      subroutine remove_ejecta_one_cell_per_step(id, ierr, j)
         integer, intent(in) :: id
         integer, intent(out) :: ierr, j
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return         
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
         
         if(.not.restart) then
            gamma1_cntr_pulse_start = HUGE(gamma1_cntr_pulse_start)
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

         select case (s% x_integer_ctrl(1))
         case(7)
            testhub_extras_names(1) = 'gamma1_cntr_pulse_start'
            testhub_extras_vals(1) = gamma1_cntr_pulse_start
         end select

         call test_suite_after_evolve(s, ierr)
         if (ierr /= 0) return         
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

         how_many_extra_history_columns = 1

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

         names(1) = 'integral_gamma1'
         vals(1) = gamma1_integral(s)

      end subroutine data_for_extra_history_columns


      real(dp) function gamma1_integral(s)
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: integral_norm
         ! Pressure weighted average of Gamma1-4/3
         ! See https://ui.adsabs.harvard.edu/abs/2020A%26A...640A..56R/abstract

         integral_norm = 0d0
         gamma1_integral = 0d0
         do k=1,s% nz
            integral_norm = integral_norm + s% Peos(k)*s% dm(k)/s% rho(k)
            gamma1_integral = gamma1_integral + &
               (s% gamma1(k)-4.d0/3.d0)*s% Peos(k)*s% dm(k)/s% rho(k)
         end do
         gamma1_integral = gamma1_integral/max(1d-99,integral_norm)

      end function gamma1_integral

      
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

         if (s% total_energy > 0) then
            extras_finish_step = terminate
            termination_code_str(t_xtra1) = 'Star is unbound'
            s% termination_code = t_xtra1
         end if

         select case (s% x_integer_ctrl(1))
         case(7)
            if(gamma1_cntr_pulse_start > 1d50 .and. gamma1_integral(s) < 0.0) then
               gamma1_cntr_pulse_start = s% gamma1(s% nz)-4.d0/3.d0
            end if
         end select

      end function extras_finish_step


      subroutine extras_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
 
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
 
         select case (s% x_integer_ctrl(1))
         case(7)
            read(iounit,iostat=ierr) gamma1_cntr_pulse_start
         end select
 
       end subroutine extras_photo_read
 
       subroutine extras_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
 
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
 
         select case (s% x_integer_ctrl(1))
         case(7)
            write(iounit) gamma1_cntr_pulse_start
         end select
 
       end subroutine extras_photo_write



   end module run_star_extras
      
