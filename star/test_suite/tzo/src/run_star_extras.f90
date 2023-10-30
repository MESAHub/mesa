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
      
      implicit none
      
      include 'test_suite_extras_def.inc'
      
      logical :: use_hydro
      integer :: inlist_part
      real(dp) :: eta_ledd, eta_medd,target_mass
      
      contains

      include 'test_suite_extras.inc'
      
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


         s% other_cgrav => my_other_cgrav

         ! Store user controls

         use_hydro = s% x_logical_ctrl(1)
         inlist_part = s%x_integer_ctrl(1)
         eta_ledd = s% x_ctrl(1)
         eta_medd = s% x_ctrl(2)
         target_mass = s% x_ctrl(3)


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


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: l_center, m_center
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0

         select case(inlist_part)

            case(4)
               ! Turn on hydro after 100 years
               if(s% star_age>100d0 .and. .not. s% v_flag .and. use_hydro) then
                  call star_set_v_flag(s% id, .true., ierr)
                  if(ierr/=0) return
               end if
      
            ! Eddington L
            l_center =  (4d0 * pi * standard_cgrav * s% m_center * clight)/s% opacity(s% nz)
      
            ! eddington limited accretion
            m_center = eta_medd * l_center / (clight*clight)
      
            s% m_center = s% m_center + m_center * s% dt
            s% xmstar = s% mstar - s% M_center
      
            s% L_center = eta_ledd * l_center
   
         end select


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
         
         if(s% m_center/msun > target_mass .and. inlist_part == 4) then
            termination_code_str(t_xtra1) = 'PASS: Have reached requested NS mass'
            s% termination_code = t_xtra1
            extras_check_model = terminate
         end if

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
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step
      

   ! use the Tolman–Oppenheimer–Volkoff (TOV) equation. 
   ! See first equation in https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation.  
   ! want to replace -G*m/r^2 by -G*m/r^2*(1 + P/(rho c^2))(1 + 4 pi r^3 P /(m c^2))/(1 - 2 G m/(r c^2))
      subroutine my_other_cgrav(id, ierr)
         use star_def
         use utils_lib, only: is_bad
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         integer :: k
         real(dp) :: r, G, m, rho, P, f1, f2, f3
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         do k=1,s% nz
            r = s% r(k)
            G = standard_cgrav
            m = s% m(k) ! m_grav hasn't been set when other_cgrav is called
            if (k > 1) then ! get approximate P and rho at face
               rho = 0.5d0*(s% rho(k-1) + s% rho(k))
               P = 0.5d0*(s% Peos(k-1) + s% Peos(k))
               f1 = 1d0 + P/(rho*clight**2)
               f2 = 1d0 + 4d0*pi*r**3*P/(m*clight**2)
            else ! k == 1
               f1 = 1d0
               f2 = 1d0
            end if
            f3 = 1d0 - 2d0*G*m/(r*clight**2)
            s% cgrav(k) = G*f1*f2/f3
            if (is_bad(s% cgrav(k))) then
               write(*,2) 's% cgrav(k)', k, s% cgrav(k)
               write(*,2) 'f1', k, f1
               write(*,2) 'f2', k, f2
               write(*,2) 'f3', k, f3
               write(*,2) 'r', k, r
               write(*,2) 'G', k, G
               write(*,2) 'm', k, m
               if (k > 1) then
                  write(*,2) 'P', k, P
                  write(*,2) 'rho', k, rho
               end if
               stop 'my_other_cgrav'
            end if
         end do
         !write(*,*) 'done my_other_cgrav', s% model_number
      end subroutine my_other_cgrav
      

      end module run_star_extras
      
