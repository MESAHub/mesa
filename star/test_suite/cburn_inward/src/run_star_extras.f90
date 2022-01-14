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

      ! Tracks quanties when the flame ignited
      real(dp) :: ign_mass, ign_density, ign_co_core_mass, flame_mass

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

         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write

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
            ign_mass = -1
            ign_density = -1
            ign_co_core_mass = -1
            flame_mass = -1
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


         select case (s% x_integer_ctrl(1)) 

         case(2) ! inlist_cburn_inward

            ! Information for testhub
            testhub_extras_names(1) = 'ign_mass'
            testhub_extras_names(2) = 'ign_log_density'
            testhub_extras_names(3) = 'ign_co_core_mass'

            testhub_extras_vals(1) = ign_mass/msun
            testhub_extras_vals(2) = safe_log10(ign_density)
            testhub_extras_vals(3) = ign_co_core_mass

         end select



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
         how_many_extra_history_columns = 4
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
         names(1) = 'ign_mass'
         names(2) = 'log(ign_den)'
         names(3) = 'ign_co_core_mass'
         names(4) = 'flame_mass'
         vals(1) = ign_mass/msun
         vals(2) = safe_log10(ign_density)
         vals(3) = ign_co_core_mass
         vals(4) = flame_mass/msun
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
         integer :: ierr, k
         type (star_info), pointer :: s
         integer :: flame_cell

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      
         
         ! Store initial flame location 
         select case (s% x_integer_ctrl(1)) 

         case(2) ! inlist_cburn_inward
            flame_cell = -1
            ! Check to see if carbon has ignited
             do k=s%nz, 1, -1
                if (has_ignited(s, k)) then 
                   flame_cell = k 
                   exit
                end if
            end do

            if(flame_cell>0) then
               ! Initial flame location
               if(ign_mass < 0) then
                  ign_mass = s% m(flame_cell)
                  ign_density = s% rho(flame_cell)
                  ign_co_core_mass = s% co_core_mass
               end if

               flame_mass = s%m(flame_cell)

               ! Final flame location
               if(ign_mass > 0d0 .and. s% m(flame_cell) < 0.5d0*ign_mass) then
                  extras_finish_step = terminate
                  write(*,'(a)') "Terminate as flame reached half way" 
                  s% termination_code = t_extras_finish_step
               end if
            end if

         end select

      end function extras_finish_step
      
      
      logical function has_ignited(s, k)
         use net_def
         use chem_def
         use chem_lib
         implicit none
         type (star_info), pointer,intent(in) :: s
         integer,intent(in) :: k
         real(dp) :: neAbun,naAbun,mgAbun,heAbun
         real(dp) :: netEng,ne_burn,o_burn
         
         has_ignited = .false.
         if(s% co_core_mass > 0d0) then
            if(s%m(k)/Msun < s%co_core_mass)THEN
               netEng = star_get_profile_output(s,'net_nuclear_energy',k)
         
               if(netEng >= 0.0)THEN
                  neAbun = s%xa(s%net_iso(chem_get_iso_id('ne20')),k)
                  naAbun = s%xa(s%net_iso(chem_get_iso_id('na23')),k)
                  mgAbun = s%xa(s%net_iso(chem_get_iso_id('mg24')),k)
                  heAbun = s%xa(s%net_iso(chem_get_iso_id('he4')),k)
                  if(neAbun > naAbun .and. naAbun > mgAbun .and. heAbun < 1d-8)THEN
                     has_ignited=.true.
                     return
                  end if
               end if
            end if
         end if
      end function has_ignited

      subroutine extras_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
 
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
 
         select case (s% x_integer_ctrl(1))
         case(2)
            read(iounit,iostat=ierr) ign_mass, ign_density, ign_co_core_mass,flame_mass
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
         case(2)
            write(iounit) ign_mass, ign_density, ign_co_core_mass,flame_mass
         end select
 
       end subroutine extras_photo_write




      end module run_star_extras
      
