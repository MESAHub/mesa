! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
      
      ! these routines are called by the standard run_star check_model
      contains

      include "test_suite_extras.inc"
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

         s% other_adjust_mdot => accretor_adjust_mdot
         s% lxtra(1) = .false.
         s% lxtra(2) = .false.
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         
      end subroutine extras_controls

      subroutine accretor_adjust_mdot(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         call star_ptr(id, s, ierr)
         s% accreted_material_j = &
              s% x_ctrl(1)*sqrt(s% cgrav(1) * s% mstar * s% photosphere_r*Rsun)

         !write(*,*) "debug", s% mstar_dot/Msun*secyer, 10**(s% x_ctrl(2))
         s% mstar_dot = s% mstar_dot + pow(10d0, s% x_ctrl(2))*Msun/secyer

      end subroutine accretor_adjust_mdot
      
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
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: am_error
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         

         ! Check accretion of angular momentum
         if (s% model_number > 1 .and. s% mstar_dot > 0) then
            if (.false.) then ! off for now
               write(*,*) "Total accreted J should be:", s% accreted_material_j*s% mstar_dot*s% dt
               write(*,*) "Current J, old J, delta J:", s% total_angular_momentum, &
                   s% total_angular_momentum_old, &
                   s% total_angular_momentum - s% total_angular_momentum_old
            end if
            am_error = (s% accreted_material_j*s% mstar_dot*s% dt &
                - (s% total_angular_momentum - s% total_angular_momentum_old)) &
                / (s% accreted_material_j*s% mstar_dot*s% dt)
            if (.true.) write(*,*) "Relative diff in accreted J vs expected:", am_error
            ! Ignore large error if per step change is only small part of total angular momentum
            ! otherwise can run into floating point precision
            if (abs(am_error) > 1d-5 .and. s% accreted_material_j*s% mstar_dot*s% dt > 1e-5*s% total_angular_momentum) then
                extras_check_model = terminate
                write(*,*) "Error in accreted J is too high!"
            end if
         end if

         ! Check switches from mass accretion to mass loss
         ! (and vice-versa) in implicit mdot calculation
         if (s% model_number > 1 .and. s% generations > 1) then
             if (s% mstar_dot < 0 .and. s% mstar_dot_old > 0) s% lxtra(1) = .true.
             if (s% mstar_dot > 0 .and. s% mstar_dot_old < 0) s% lxtra(2) = .true.
         end if

         if (s% center_h1 < 1d-6) then
             ! Check if implicit wind has switched from mass accretion to mass loss
             if (s% lxtra(1)) then
                 write(*,*) "Test has finished with no problems."
                 s% termination_code = t_xa_central_lower_limit
             else
                 write(*,*) "Apparent problem with implicit wind.", s% lxtra(1), s% lxtra(2)
                 s% termination_code = t_extras_check_model
             end if
             extras_check_model = terminate
         end if

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depenending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination conditon'

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
         
         !note: do NOT add these names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
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

         include 'formats'

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated

         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      

      end module run_star_extras
      
