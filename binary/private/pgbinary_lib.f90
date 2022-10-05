! ***********************************************************************
!
!   Copyright (C) 2013-2022  Bill Paxton, Pablo Marchant, Matthias Fabry
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


module pgbinary_lib

   implicit none

contains

   ! pgbinary interface
   subroutine start_new_run_for_pgbinary(b, ierr) ! reset logs
      use pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      call do_start_new_run_for_pgbinary(b, ierr)
   end subroutine start_new_run_for_pgbinary


   subroutine restart_run_for_pgbinary(b, ierr)
      use pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      call do_restart_run_for_pgbinary(b, ierr)
   end subroutine restart_run_for_pgbinary


   subroutine read_pgbinary_controls(b, ierr)
      use pgbinary, only : do_read_pgbinary_controls
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      call do_read_pgbinary_controls(b, 'inlist', ierr)
   end subroutine read_pgbinary_controls


   subroutine read_pgbinary_inlist(b, inlist_fname, ierr)
      use pgbinary, only : do_read_pgbinary_controls
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character(*), intent(in) :: inlist_fname
      integer, intent(out) :: ierr
      call do_read_pgbinary_controls(b, inlist_fname, ierr)
   end subroutine read_pgbinary_inlist


   subroutine update_pgbinary_plots(b, must_write_files, ierr)
      use pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      logical, intent(in) :: must_write_files
      integer, intent(out) :: ierr
      call do_pgbinary_plots(&
         b, must_write_files, &
         ierr)
   end subroutine update_pgbinary_plots


   subroutine create_pgbinary_file_name(b, dir, prefix, name)
      use pgbinary, only : do_create_file_name
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: dir, prefix
      character (len = *), intent(out) :: name
      call do_create_file_name(b, dir, prefix, name)
   end subroutine create_pgbinary_file_name


   subroutine pgbinary_write_plot_to_file(b, p, filename, ierr)
      use binary_def, only : pgbinary_win_file_data, binary_info
      use pgbinary, only : do_write_plot_to_file
      type (binary_info), pointer :: b
      type (pgbinary_win_file_data), pointer :: p
      character (len = *), intent(in) :: filename
      integer, intent(out) :: ierr
      call do_write_plot_to_file(b, p, filename, ierr)
   end subroutine pgbinary_write_plot_to_file


   subroutine show_pgbinary_annotations(&
      b, show_annotation1, show_annotation2, show_annotation3)
      use pgbinary, only : do_show_pgbinary_annotations
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      logical, intent(in) :: &
         show_annotation1, show_annotation2, show_annotation3
      call do_show_pgbinary_annotations(&
         b, show_annotation1, show_annotation2, show_annotation3)
   end subroutine show_pgbinary_annotations


   subroutine pgbinary_show_box(b, str1, str2)
      use pgbinary, only : show_box_pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: str1, str2
      call show_box_pgbinary(b, str1, str2)
   end subroutine pgbinary_show_box


   subroutine pgbinary_show_title(b, title, pad)
      use pgbinary, only : show_title_pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: title
      real, intent(in) :: pad
      optional pad
      real :: pad_arg
      pad_arg = 0
      if (present(pad)) pad_arg = pad
      call show_title_pgbinary(b, title, pad_arg)
   end subroutine pgbinary_show_title


   subroutine pgbinary_show_xaxis_label(b, label, pad)
      use pgbinary, only : show_xaxis_label_pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
      real :: pad_arg
      pad_arg = 0
      if (present(pad)) pad_arg = pad
      call show_xaxis_label_pgbinary(b, label, pad_arg)
   end subroutine pgbinary_show_xaxis_label


   subroutine pgbinary_show_left_yaxis_label(b, label, pad)
      use pgbinary, only : show_left_yaxis_label_pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
      real :: pad_arg
      pad_arg = 0
      if (present(pad)) pad_arg = pad
      call show_left_yaxis_label_pgbinary(b, label, pad_arg)
   end subroutine pgbinary_show_left_yaxis_label


   subroutine pgbinary_show_right_yaxis_label(b, label, pad)
      use pgbinary, only : show_right_yaxis_label_pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
      real :: pad_arg
      pad_arg = 0
      if (present(pad)) pad_arg = pad
      call show_right_yaxis_label_pgbinary(b, label, pad_arg)
   end subroutine pgbinary_show_right_yaxis_label


   subroutine pgbinary_show_left_axis_label_pgmtxt(&
      b, coord, fjust, label, pad)
      use pgbinary, only : show_left_yaxis_label_pgmtxt_pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: pad_arg
      pad_arg = 0
      if (present(pad)) pad_arg = pad
      call show_left_yaxis_label_pgmtxt_pgbinary(&
         b, coord, fjust, label, pad)
   end subroutine pgbinary_show_left_axis_label_pgmtxt


   subroutine pgbinary_show_right_axis_label_pgmtxt(&
      b, coord, fjust, label, pad)
      use pgbinary, only : show_right_yaxis_label_pgmtxt_pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: pad_arg
      pad_arg = 0
      if (present(pad)) pad_arg = pad
      call show_right_yaxis_label_pgmtxt_pgbinary(&
         b, coord, fjust, label, pad)
   end subroutine pgbinary_show_right_axis_label_pgmtxt


   subroutine pgbinary_show_model_number(b)
      use pgbinary, only : show_model_number_pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      call show_model_number_pgbinary(b)
   end subroutine pgbinary_show_model_number


   subroutine pgbinary_show_age(b)
      use pgbinary, only : show_age_pgbinary
      use binary_def, only : binary_info
      type (binary_info), pointer :: b
      call show_age_pgbinary(b)
   end subroutine pgbinary_show_age


   subroutine binary_shutdown_pgbinary(id, ierr)
      use pgbinary, only : shutdown_pgbinary
      use binary_def, only : binary_info, binary_ptr
      integer, intent(in) :: id ! id for star
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call shutdown_pgbinary(b)
   end subroutine binary_shutdown_pgbinary

end module pgbinary_lib
