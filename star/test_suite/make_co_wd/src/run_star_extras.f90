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
      use auto_diff
      
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
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         s% lxtra(1) = .false.
         s% lxtra(2) = restart

         ! Only allow diffusion after burning from thermal pulse has settled down
         if(s% ctrl% x_integer_ctrl(1) == 3) then
            if(s% power_he_burn > 1d3) then
               s% ctrl% do_element_diffusion = .false.
            else
               s% ctrl% do_element_diffusion = .true.
            end if
         end if

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
         use chem_def, only: icc
         integer, intent(in) :: id
         integer :: ierr
         
         real(dp), parameter :: Blocker_scaling_factor_after_TP = 5d0
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         include 'formats'
         
         extras_check_model = keep_going

         ! check for c12_c12 burning that will make an ONe WD.
         if(s% L_by_category(icc) > 1d4) then
            write(*,'(A)')
            write(*,*) "This model is too massive." 
            write(*,*) "Carbon has ignited in the interior and will produce an ONe WD."
            write(*,'(A)')
            extras_check_model = terminate
            s% termination_code = t_extras_check_model
         end if

         ! Only allow diffusion after burning from thermal pulse has settled down
         if(s% ctrl% x_integer_ctrl(1) == 3) then
            if(s% power_he_burn > 1d3) then
               s% ctrl% do_element_diffusion = .false.
            else
               s% ctrl% do_element_diffusion = .true.
            end if
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
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         real(dp) :: bottom_mass, H_env_limit
         integer :: ierr, k, kbot
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         if (s% ctrl% x_integer_ctrl(1) == 2) then ! part2
            H_env_limit = s% ctrl% x_ctrl(1)
            if (.not. (s% lxtra(1) .or. s% lxtra(2))) then
               ! find mass coordinate for cut
               kbot = 1
               do k = s% nz, 1, -1
                  if (s% X(k) > s% ctrl% he_core_boundary_h1_fraction) then
                     bottom_mass = s% m(k)
                     kbot = k
                     exit
                  end if
               end do
               do k = kbot, 1, -1
                  if (s% m(k) > bottom_mass + H_env_limit*Msun) then
                  
                     write(*,2) 'call star_remove_surface_at_cell_k', k, s% m(k)/Msun
                     call star_remove_surface_at_cell_k(s% id, k, ierr)
                     if (ierr /= 0) then
                        write(*,*) 'failed in star_remove_surface_at_cell_k', k
                        extras_finish_step = terminate
                        return
                     end if
                     
                     !call star_relax_to_star_cut(s% id, k, .false., .true., .true., ierr)
                     
                     
                     s% lxtra(1) = .true.
                     exit
                  end if
               end do
            end if
         end if
         
      end function extras_finish_step
      
      

      end module run_star_extras
      
