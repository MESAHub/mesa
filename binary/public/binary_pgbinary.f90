! ***********************************************************************
!
!   Copyright (C) 2022  Matthias Fabry
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

module binary_pgbinary

   use const_def
   use star_pgstar

   abstract interface

      subroutine pgbinary_plot_interface(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
      end subroutine pgbinary_plot_interface

      subroutine other_do_plot_in_binary_grid_interface(&
         id, device_id, xleft, xright, ybot, ytop, txt_scale, ierr)
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         integer, intent(out) :: ierr
      end subroutine other_do_plot_in_binary_grid_interface

      subroutine pgbinary_decorator_interface(id, xmin, xmax, ymin, ymax, plot_num, ierr)
         integer, intent(in) :: id
         !Not doubles
         real, intent(in) :: xmin, xmax, ymin, ymax
         integer, intent(in) :: plot_num
         integer, intent(out) :: ierr
      end subroutine pgbinary_decorator_interface

   end interface

   type pgbinary_win_file_data
      integer :: id
      character (len = 64) :: name
      logical :: win_flag, file_flag, do_win, do_file
      integer :: id_win, id_file, file_interval
      real :: win_width, prev_win_width
      real :: win_aspect_ratio, prev_win_aspect_ratio
      real :: file_width, file_aspect_ratio
      character (len = strlen) :: file_dir, file_prefix, most_recent_filename
      character (len = strlen) :: file_dir_for_previous_mkdir
      logical :: have_called_mkdir
      procedure (pgbinary_plot_interface), pointer, nopass :: plot => null()
      ! the following make it possible to use "other" plots in grids
      logical :: okay_to_call_do_plot_in_binary_grid
      procedure (other_do_plot_in_binary_grid_interface), pointer, nopass :: &
         do_plot_in_binary_grid => null()
   end type pgbinary_win_file_data

   integer, parameter :: max_num_Binary_History_Panels = 10
   integer, parameter :: max_num_Summary_Binary_History_Lines = 16
   integer, parameter :: max_num_Binary_Other_plots = 10

   integer, parameter :: i_Summary_Binary_History = 1

   integer, parameter :: i_Binary_Text_Summary1 = i_Summary_Binary_History + 1
   integer, parameter :: i_Binary_Text_Summary2 = i_Binary_Text_Summary1 + 1
   integer, parameter :: i_Binary_Text_Summary3 = i_Binary_Text_Summary2 + 1
   integer, parameter :: i_Binary_Text_Summary4 = i_Binary_Text_Summary3 + 1
   integer, parameter :: i_Binary_Text_Summary5 = i_Binary_Text_Summary4 + 1
   integer, parameter :: i_Binary_Text_Summary6 = i_Binary_Text_Summary5 + 1
   integer, parameter :: i_Binary_Text_Summary7 = i_Binary_Text_Summary6 + 1
   integer, parameter :: i_Binary_Text_Summary8 = i_Binary_Text_Summary7 + 1
   integer, parameter :: i_Binary_Text_Summary9 = i_Binary_Text_Summary8 + 1

   integer, parameter :: i_Binary_Hist_Track1 = i_Binary_Text_Summary9 + 1
   integer, parameter :: i_Binary_Hist_Track2 = i_Binary_Hist_Track1 + 1
   integer, parameter :: i_Binary_Hist_Track3 = i_Binary_Hist_Track2 + 1
   integer, parameter :: i_Binary_Hist_Track4 = i_Binary_Hist_Track3 + 1
   integer, parameter :: i_Binary_Hist_Track5 = i_Binary_Hist_Track4 + 1
   integer, parameter :: i_Binary_Hist_Track6 = i_Binary_Hist_Track5 + 1
   integer, parameter :: i_Binary_Hist_Track7 = i_Binary_Hist_Track6 + 1
   integer, parameter :: i_Binary_Hist_Track8 = i_Binary_Hist_Track7 + 1
   integer, parameter :: i_Binary_Hist_Track9 = i_Binary_Hist_Track8 + 1

   integer, parameter :: i_Binary_Hist_Panels1 = i_Binary_Hist_Track9 + 1
   integer, parameter :: i_Binary_Hist_Panels2 = i_Binary_Hist_Panels1 + 1
   integer, parameter :: i_Binary_Hist_Panels3 = i_Binary_Hist_Panels2 + 1
   integer, parameter :: i_Binary_Hist_Panels4 = i_Binary_Hist_Panels3 + 1
   integer, parameter :: i_Binary_Hist_Panels5 = i_Binary_Hist_Panels4 + 1
   integer, parameter :: i_Binary_Hist_Panels6 = i_Binary_Hist_Panels5 + 1
   integer, parameter :: i_Binary_Hist_Panels7 = i_Binary_Hist_Panels6 + 1
   integer, parameter :: i_Binary_Hist_Panels8 = i_Binary_Hist_Panels7 + 1
   integer, parameter :: i_Binary_Hist_Panels9 = i_Binary_Hist_Panels8 + 1

   integer, parameter :: i_Binary_Grid1 = i_Binary_Hist_Panels9 + 1
   integer, parameter :: i_Binary_Grid2 = i_Binary_Grid1 + 1
   integer, parameter :: i_Binary_Grid3 = i_Binary_Grid2 + 1
   integer, parameter :: i_Binary_Grid4 = i_Binary_Grid3 + 1
   integer, parameter :: i_Binary_Grid5 = i_Binary_Grid4 + 1
   integer, parameter :: i_Binary_Grid6 = i_Binary_Grid5 + 1
   integer, parameter :: i_Binary_Grid7 = i_Binary_Grid6 + 1
   integer, parameter :: i_Binary_Grid8 = i_Binary_Grid7 + 1
   integer, parameter :: i_Binary_Grid9 = i_Binary_Grid8 + 1

   integer, parameter :: i_Binary_Star1 = i_Binary_Grid9 + 1
   integer, parameter :: i_Binary_Star2 = i_Binary_Star1 + 1

   integer, parameter :: i_Binary_Orbit = i_Binary_Star2 + 1
   integer, parameter :: i_Binary_Other = i_Binary_Orbit + 1

   integer, parameter :: num_pgbinary_plots = i_Binary_Other + &
      max_num_Binary_Other_plots

   integer, parameter :: max_num_pgbinary_Grid_plots = 10

   integer, parameter :: max_num_pgbinary_trace_Binary_History_values = 20

   type pgbinary_hist_node
      real(dp) :: age
      integer :: step
      real(dp), pointer :: vals(:) => null() ! values of items in history_columns list
      type (pgbinary_hist_node), pointer :: next => null()
      ! list kept in strictly decreasing order by age & step
   end type pgbinary_hist_node

   type pgbinary_controls
      include "pgbinary_controls.inc"

      type (pgbinary_hist_node), pointer :: pgbinary_hist => null()

      procedure(pgbinary_decorator_interface), pointer, nopass :: &
            History_Panels1_pgbinary_decorator => null(), &
            History_Panels2_pgbinary_decorator => null(), &
            History_Panels3_pgbinary_decorator => null(), &
            History_Panels4_pgbinary_decorator => null(), &
            History_Panels5_pgbinary_decorator => null(), &
            History_Panels6_pgbinary_decorator => null(), &
            History_Panels7_pgbinary_decorator => null(), &
            History_Panels8_pgbinary_decorator => null(), &
            History_Panels9_pgbinary_decorator => null(), &
            History_Track1_pgbinary_decorator => null(), &
            History_Track2_pgbinary_decorator => null(), &
            History_Track3_pgbinary_decorator => null(), &
            History_Track4_pgbinary_decorator => null(), &
            History_Track5_pgbinary_decorator => null(), &
            History_Track6_pgbinary_decorator => null(), &
            History_Track7_pgbinary_decorator => null(), &
            History_Track8_pgbinary_decorator => null(), &
            History_Track9_pgbinary_decorator => null(), &
            summary_history_pgbinary_decorator => null()

   end type pgbinary_controls

end module binary_pgbinary
