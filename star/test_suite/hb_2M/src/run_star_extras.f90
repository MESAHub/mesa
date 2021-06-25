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
      
      implicit none

      real(dp) :: mass_conv_core_y050
      
      include "test_suite_extras_def.inc"
      
      contains

      include "test_suite_extras.inc"
      
      
      subroutine extras_photo_read(id, iounit, ierr)
        integer, intent(in) :: id, iounit
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        ierr = 0

        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        select case (s% x_integer_ctrl(1))
        case(3)
           read(iounit,iostat=ierr) mass_conv_core_y050
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
        case(3)
           write(iounit) mass_conv_core_y050
        end select

      end subroutine extras_photo_write


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write
         
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
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! initialize mass_conv_core_y050 to "unset" value
         if (.not. restart) then
            mass_conv_core_y050 = -1
         end if

         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: min_mass_conv_core, max_mass_conv_core
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         include 'formats'

         select case (s% x_integer_ctrl(1))
         case (3) ! inlist_hb_2M

            ! put target info in TestHub output
            testhub_extras_names(1) = 'mass_conv_core_y050'; testhub_extras_vals(1) = mass_conv_core_y050

            ! display value to user
            write(*,*)
            write(*,'(A70, F8.3)') '[TestHub] Convective core mass at Yc = 0.5 (Msun): ', mass_conv_core_y050
            write(*,*)

            ! get target range from inlist
            min_mass_conv_core = s% x_ctrl(1)
            max_mass_conv_core = s% x_ctrl(2)

            ! check if value is outside of the target range
            if ((mass_conv_core_y050 .lt. min_mass_conv_core) .or. (mass_conv_core_y050 .gt. max_mass_conv_core)) then
               write(*,*) 'bad value for mass_conv_core_y050'
               write(*,1) 'min allowed value', min_mass_conv_core
               write(*,1) 'mass_conv_core_y050', mass_conv_core_y050
               write(*,1) 'max allowed value', max_mass_conv_core
            else
               write(*,'(a)') 'all values are within tolerance'
            end if
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
         vals(1) = mass_conv_core_y050
         names(1) = 'mass_conv_core_y050'
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
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
         names(1) = 'zbar_div_abar'
         do k=1,s% nz
            vals(k,1) = s% zbar(k)/s% abar(k)
         end do
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

         ! during part 3: CHeB
         if (s% x_integer_ctrl(1) == 3) then
            ! save core mass the first time Yc < 0.5
            if (s% center_he4 .lt. 0.5d0 .and. mass_conv_core_y050 .lt. 0) then
               mass_conv_core_y050 = s% mass_conv_core
            end if
         end if

      end function extras_finish_step
      
      

      end module run_star_extras
      
