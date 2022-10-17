! ***********************************************************************
!
!   Copyright (C) 2017-2019  Rob Farmer & The MESA Team
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
   use colors_lib
   
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
      
      s% extras_startup => extras_startup
      s% extras_check_model => extras_check_model
      s% extras_finish_step => extras_finish_step
      s% extras_after_evolve => extras_after_evolve
      s% how_many_extra_history_columns => how_many_extra_history_columns
      s% data_for_extra_history_columns => data_for_extra_history_columns
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns

#ifdef USE_PGPLOT
         !Add custom decorator to pgplots
         s% color_magnitude1_pgstar_decorator => col_mag1_decorator
#endif
   
   end subroutine extras_controls

#ifdef USE_PGPLOT
      subroutine col_mag1_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
         integer, intent(in) :: id
         !Not dp
         real,intent(in) :: xmin, xmax, ymin, ymax 
         real :: xcenter,ycenter,dx,dy,a
         integer, intent(in) :: plot_num
         integer, intent(out) :: ierr
         integer :: i
         type (star_info), pointer :: s
         character(len=20) :: temp
        
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
        
         dx=(xmax-xmin)
         dy=(ymax-ymin)
       
         xcenter = xmin + dx*0.5d0
         ycenter = ymin + dy*0.5d0
     
         call pgsci(clr_Coral)
        
        !Add stuff to the top panel in color_magnitude1
        if(plot_num==1) Then        
            a = 0.1d0
            call pgline(5, (/xcenter-a*dx,xcenter-a*dx,xcenter+a*dx,xcenter+a*dx,xcenter-a*dx/),&
                            (/ycenter-a*dy,ycenter+a*dy,ycenter+a*dy,ycenter-a*dy,ycenter-a*dy/))
            
        else
        ! Second or higher panel, this function gets called once per panel for the color_magnitude1 plot, so
        ! num_panel distinguishes between each panel
            xcenter=xmin+dx*0.75d0
            ycenter=ymin+dy*0.25d0
            write(temp, '(f10.2)') log10(s%T(s%nz))
            
            !xcenter, ycenter is the position, The 45 rotates the text, 0.0 is a padding number then the string follows that
            call PGPTXT(xcenter,ycenter,45.0d0,0.0d0, 'log T\dc\u='//trim(temp))
        end if
        
        
      end subroutine col_mag1_decorator    
#endif      
   
   
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
      
      write(*, '(a)') 'finished custom colors'
      
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
      character (len = maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      
      !Here we get the fake av data
      names(1) = 'av_v' ! Same name as used in the fake_av_v.txt file for the column
      
      vals(1) = get_bc_by_name(names(1), safe_log10(s% T(1)), &
         safe_log10(s% grav(1)), &
         ! Normally we have the metalicity as the third parameter here,
         ! but that is not required. We do not need the Teff or logg either,
         ! we could do the interpolation over three other parameters or inputs.
         s%job%extras_rpar(1), ierr)
   
   end subroutine data_for_extra_history_columns
   
   
   integer function how_many_extra_profile_columns(id)
      use star_def, only : star_info
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 0
   end function how_many_extra_profile_columns
   
   
   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      use star_def, only : star_info, maxlen_profile_column_name
      use const_def, only : dp
      integer, intent(in) :: id, n, nz
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine data_for_extra_profile_columns
   
   
   ! returns either keep_going, retry, or terminate.
   integer function extras_finish_step(id)
      use chem_def
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_finish_step = keep_going
   end function extras_finish_step


end module run_star_extras
      
