! ***********************************************************************
!
!   Copyright (C) 2013-2022  The MESA Team, Bill Paxton & Matthias Fabry
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module pgbinary_hist_track

   use binary_private_def
   use const_def
   use pgbinary_support

   implicit none


contains


   subroutine History_Track1_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track1_plot(b, id, device_id, &
         b% pg% History_Track_xleft(1), b% pg% History_Track_xright(1), &
         b% pg% History_Track_ybot(1), b% pg% History_Track_ytop(1), .false., &
         b% pg% History_Track_title(1), b% pg% History_Track_txt_scale(1), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track1_plot


   subroutine do_History_Track1_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Track_xname(1), b% pg% History_Track_yname(1), &
         b% pg% History_Track_xaxis_label(1), b% pg% History_Track_yaxis_label(1), &
         b% pg% History_Track_xmin(1), b% pg% History_Track_xmax(1), &
         b% pg% History_Track_xmargin(1), b% pg% History_Track_dxmin(1), &
         b% pg% History_Track_ymin(1), b% pg% History_Track_ymax(1), &
         b% pg% History_Track_ymargin(1), b% pg% History_Track_dymin(1), &
         b% pg% History_Track_step_min(1), b% pg% History_Track_step_max(1), &
         b% pg% History_Track_reverse_xaxis(1), b% pg% History_Track_reverse_yaxis(1), &
         b% pg% History_Track_log_xaxis(1), b% pg% History_Track_log_yaxis(1), &
         b% pg% show_History_Track_target_box(1), b% pg% History_Track_n_sigma(1), &
         b% pg% History_Track_xtarget(1), b% pg% History_Track_ytarget(1), &
         b% pg% History_Track_xsigma(1), b% pg% History_Track_ysigma(1), &
         b% pg% show_History_Track_annotation1(1), &
         b% pg% show_History_Track_annotation2(1), &
         b% pg% show_History_Track_annotation3(1), &
         b% pg% History_Track_fname(1), &
         b% pg% History_Track_use_decorator(1), b% pg% History_Track1_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_History_Track1_plot


      subroutine History_Track2_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track2_plot(b, id, device_id, &
         b% pg% History_Track_xleft(2), b% pg% History_Track_xright(2), &
         b% pg% History_Track_ybot(2), b% pg% History_Track_ytop(2), .false., &
         b% pg% History_Track_title(2), b% pg% History_Track_txt_scale(2), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track2_plot


   subroutine do_History_Track2_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Track_xname(2), b% pg% History_Track_yname(2), &
         b% pg% History_Track_xaxis_label(2), b% pg% History_Track_yaxis_label(2), &
         b% pg% History_Track_xmin(2), b% pg% History_Track_xmax(2), &
         b% pg% History_Track_xmargin(2), b% pg% History_Track_dxmin(2), &
         b% pg% History_Track_ymin(2), b% pg% History_Track_ymax(2), &
         b% pg% History_Track_ymargin(2), b% pg% History_Track_dymin(2), &
         b% pg% History_Track_step_min(2), b% pg% History_Track_step_max(2), &
         b% pg% History_Track_reverse_xaxis(2), b% pg% History_Track_reverse_yaxis(2), &
         b% pg% History_Track_log_xaxis(2), b% pg% History_Track_log_yaxis(2), &
         b% pg% show_History_Track_target_box(2), b% pg% History_Track_n_sigma(2), &
         b% pg% History_Track_xtarget(2), b% pg% History_Track_ytarget(2), &
         b% pg% History_Track_xsigma(2), b% pg% History_Track_ysigma(2), &
         b% pg% show_History_Track_annotation1(2), &
         b% pg% show_History_Track_annotation2(2), &
         b% pg% show_History_Track_annotation3(2), &
         b% pg% History_Track_fname(2), &
         b% pg% History_Track_use_decorator(2), b% pg% History_Track2_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_History_Track2_plot

      subroutine History_Track3_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track3_plot(b, id, device_id, &
         b% pg% History_Track_xleft(3), b% pg% History_Track_xright(3), &
         b% pg% History_Track_ybot(3), b% pg% History_Track_ytop(3), .false., &
         b% pg% History_Track_title(3), b% pg% History_Track_txt_scale(3), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track3_plot


   subroutine do_History_Track3_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Track_xname(3), b% pg% History_Track_yname(3), &
         b% pg% History_Track_xaxis_label(3), b% pg% History_Track_yaxis_label(3), &
         b% pg% History_Track_xmin(3), b% pg% History_Track_xmax(3), &
         b% pg% History_Track_xmargin(3), b% pg% History_Track_dxmin(3), &
         b% pg% History_Track_ymin(3), b% pg% History_Track_ymax(3), &
         b% pg% History_Track_ymargin(3), b% pg% History_Track_dymin(3), &
         b% pg% History_Track_step_min(3), b% pg% History_Track_step_max(3), &
         b% pg% History_Track_reverse_xaxis(3), b% pg% History_Track_reverse_yaxis(3), &
         b% pg% History_Track_log_xaxis(3), b% pg% History_Track_log_yaxis(3), &
         b% pg% show_History_Track_target_box(3), b% pg% History_Track_n_sigma(3), &
         b% pg% History_Track_xtarget(3), b% pg% History_Track_ytarget(3), &
         b% pg% History_Track_xsigma(3), b% pg% History_Track_ysigma(3), &
         b% pg% show_History_Track_annotation1(3), &
         b% pg% show_History_Track_annotation2(3), &
         b% pg% show_History_Track_annotation3(3), &
         b% pg% History_Track_fname(3), &
         b% pg% History_Track_use_decorator(3), b% pg% History_Track3_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_History_Track3_plot

      subroutine History_Track4_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track4_plot(b, id, device_id, &
         b% pg% History_Track_xleft(4), b% pg% History_Track_xright(4), &
         b% pg% History_Track_ybot(4), b% pg% History_Track_ytop(4), .false., &
         b% pg% History_Track_title(4), b% pg% History_Track_txt_scale(4), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track4_plot


   subroutine do_History_Track4_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Track_xname(4), b% pg% History_Track_yname(4), &
         b% pg% History_Track_xaxis_label(4), b% pg% History_Track_yaxis_label(4), &
         b% pg% History_Track_xmin(4), b% pg% History_Track_xmax(4), &
         b% pg% History_Track_xmargin(4), b% pg% History_Track_dxmin(4), &
         b% pg% History_Track_ymin(4), b% pg% History_Track_ymax(4), &
         b% pg% History_Track_ymargin(4), b% pg% History_Track_dymin(4), &
         b% pg% History_Track_step_min(4), b% pg% History_Track_step_max(4), &
         b% pg% History_Track_reverse_xaxis(4), b% pg% History_Track_reverse_yaxis(4), &
         b% pg% History_Track_log_xaxis(4), b% pg% History_Track_log_yaxis(4), &
         b% pg% show_History_Track_target_box(4), b% pg% History_Track_n_sigma(4), &
         b% pg% History_Track_xtarget(4), b% pg% History_Track_ytarget(4), &
         b% pg% History_Track_xsigma(4), b% pg% History_Track_ysigma(4), &
         b% pg% show_History_Track_annotation1(4), &
         b% pg% show_History_Track_annotation2(4), &
         b% pg% show_History_Track_annotation3(4), &
         b% pg% History_Track_fname(4), &
         b% pg% History_Track_use_decorator(4), b% pg% History_Track4_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_History_Track4_plot

      subroutine History_Track5_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track5_plot(b, id, device_id, &
         b% pg% History_Track_xleft(5), b% pg% History_Track_xright(5), &
         b% pg% History_Track_ybot(5), b% pg% History_Track_ytop(5), .false., &
         b% pg% History_Track_title(5), b% pg% History_Track_txt_scale(5), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track5_plot


   subroutine do_History_Track5_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Track_xname(5), b% pg% History_Track_yname(5), &
         b% pg% History_Track_xaxis_label(5), b% pg% History_Track_yaxis_label(5), &
         b% pg% History_Track_xmin(5), b% pg% History_Track_xmax(5), &
         b% pg% History_Track_xmargin(5), b% pg% History_Track_dxmin(5), &
         b% pg% History_Track_ymin(5), b% pg% History_Track_ymax(5), &
         b% pg% History_Track_ymargin(5), b% pg% History_Track_dymin(5), &
         b% pg% History_Track_step_min(5), b% pg% History_Track_step_max(5), &
         b% pg% History_Track_reverse_xaxis(5), b% pg% History_Track_reverse_yaxis(5), &
         b% pg% History_Track_log_xaxis(5), b% pg% History_Track_log_yaxis(5), &
         b% pg% show_History_Track_target_box(5), b% pg% History_Track_n_sigma(5), &
         b% pg% History_Track_xtarget(5), b% pg% History_Track_ytarget(5), &
         b% pg% History_Track_xsigma(5), b% pg% History_Track_ysigma(5), &
         b% pg% show_History_Track_annotation1(5), &
         b% pg% show_History_Track_annotation2(5), &
         b% pg% show_History_Track_annotation3(5), &
         b% pg% History_Track_fname(5), &
         b% pg% History_Track_use_decorator(5), b% pg% History_Track5_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_History_Track5_plot

      subroutine History_Track6_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track6_plot(b, id, device_id, &
         b% pg% History_Track_xleft(6), b% pg% History_Track_xright(6), &
         b% pg% History_Track_ybot(6), b% pg% History_Track_ytop(6), .false., &
         b% pg% History_Track_title(6), b% pg% History_Track_txt_scale(6), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track6_plot


   subroutine do_History_Track6_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Track_xname(6), b% pg% History_Track_yname(6), &
         b% pg% History_Track_xaxis_label(6), b% pg% History_Track_yaxis_label(6), &
         b% pg% History_Track_xmin(6), b% pg% History_Track_xmax(6), &
         b% pg% History_Track_xmargin(6), b% pg% History_Track_dxmin(6), &
         b% pg% History_Track_ymin(6), b% pg% History_Track_ymax(6), &
         b% pg% History_Track_ymargin(6), b% pg% History_Track_dymin(6), &
         b% pg% History_Track_step_min(6), b% pg% History_Track_step_max(6), &
         b% pg% History_Track_reverse_xaxis(6), b% pg% History_Track_reverse_yaxis(6), &
         b% pg% History_Track_log_xaxis(6), b% pg% History_Track_log_yaxis(6), &
         b% pg% show_History_Track_target_box(6), b% pg% History_Track_n_sigma(6), &
         b% pg% History_Track_xtarget(6), b% pg% History_Track_ytarget(6), &
         b% pg% History_Track_xsigma(6), b% pg% History_Track_ysigma(6), &
         b% pg% show_History_Track_annotation1(6), &
         b% pg% show_History_Track_annotation2(6), &
         b% pg% show_History_Track_annotation3(6), &
         b% pg% History_Track_fname(6), &
         b% pg% History_Track_use_decorator(6), b% pg% History_Track6_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_History_Track6_plot

      subroutine History_Track7_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track7_plot(b, id, device_id, &
         b% pg% History_Track_xleft(7), b% pg% History_Track_xright(7), &
         b% pg% History_Track_ybot(7), b% pg% History_Track_ytop(7), .false., &
         b% pg% History_Track_title(7), b% pg% History_Track_txt_scale(7), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track7_plot


   subroutine do_History_Track7_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Track_xname(7), b% pg% History_Track_yname(7), &
         b% pg% History_Track_xaxis_label(7), b% pg% History_Track_yaxis_label(7), &
         b% pg% History_Track_xmin(7), b% pg% History_Track_xmax(7), &
         b% pg% History_Track_xmargin(7), b% pg% History_Track_dxmin(7), &
         b% pg% History_Track_ymin(7), b% pg% History_Track_ymax(7), &
         b% pg% History_Track_ymargin(7), b% pg% History_Track_dymin(7), &
         b% pg% History_Track_step_min(7), b% pg% History_Track_step_max(7), &
         b% pg% History_Track_reverse_xaxis(7), b% pg% History_Track_reverse_yaxis(7), &
         b% pg% History_Track_log_xaxis(7), b% pg% History_Track_log_yaxis(7), &
         b% pg% show_History_Track_target_box(7), b% pg% History_Track_n_sigma(7), &
         b% pg% History_Track_xtarget(7), b% pg% History_Track_ytarget(7), &
         b% pg% History_Track_xsigma(7), b% pg% History_Track_ysigma(7), &
         b% pg% show_History_Track_annotation1(7), &
         b% pg% show_History_Track_annotation2(7), &
         b% pg% show_History_Track_annotation3(7), &
         b% pg% History_Track_fname(7), &
         b% pg% History_Track_use_decorator(7), b% pg% History_Track7_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_History_Track7_plot

      subroutine History_Track8_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track8_plot(b, id, device_id, &
         b% pg% History_Track_xleft(8), b% pg% History_Track_xright(8), &
         b% pg% History_Track_ybot(8), b% pg% History_Track_ytop(8), .false., &
         b% pg% History_Track_title(8), b% pg% History_Track_txt_scale(8), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track8_plot


   subroutine do_History_Track8_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Track_xname(8), b% pg% History_Track_yname(8), &
         b% pg% History_Track_xaxis_label(8), b% pg% History_Track_yaxis_label(8), &
         b% pg% History_Track_xmin(8), b% pg% History_Track_xmax(8), &
         b% pg% History_Track_xmargin(8), b% pg% History_Track_dxmin(8), &
         b% pg% History_Track_ymin(8), b% pg% History_Track_ymax(8), &
         b% pg% History_Track_ymargin(8), b% pg% History_Track_dymin(8), &
         b% pg% History_Track_step_min(8), b% pg% History_Track_step_max(8), &
         b% pg% History_Track_reverse_xaxis(8), b% pg% History_Track_reverse_yaxis(8), &
         b% pg% History_Track_log_xaxis(8), b% pg% History_Track_log_yaxis(8), &
         b% pg% show_History_Track_target_box(8), b% pg% History_Track_n_sigma(8), &
         b% pg% History_Track_xtarget(8), b% pg% History_Track_ytarget(8), &
         b% pg% History_Track_xsigma(8), b% pg% History_Track_ysigma(8), &
         b% pg% show_History_Track_annotation1(8), &
         b% pg% show_History_Track_annotation2(8), &
         b% pg% show_History_Track_annotation3(8), &
         b% pg% History_Track_fname(8), &
         b% pg% History_Track_use_decorator(8), b% pg% History_Track8_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_History_Track8_plot

      subroutine History_Track9_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track9_plot(b, id, device_id, &
         b% pg% History_Track_xleft(9), b% pg% History_Track_xright(9), &
         b% pg% History_Track_ybot(9), b% pg% History_Track_ytop(9), .false., &
         b% pg% History_Track_title(9), b% pg% History_Track_txt_scale(9), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track9_plot


   subroutine do_History_Track9_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Track_xname(9), b% pg% History_Track_yname(9), &
         b% pg% History_Track_xaxis_label(9), b% pg% History_Track_yaxis_label(9), &
         b% pg% History_Track_xmin(9), b% pg% History_Track_xmax(9), &
         b% pg% History_Track_xmargin(9), b% pg% History_Track_dxmin(9), &
         b% pg% History_Track_ymin(9), b% pg% History_Track_ymax(9), &
         b% pg% History_Track_ymargin(9), b% pg% History_Track_dymin(9), &
         b% pg% History_Track_step_min(9), b% pg% History_Track_step_max(9), &
         b% pg% History_Track_reverse_xaxis(9), b% pg% History_Track_reverse_yaxis(9), &
         b% pg% History_Track_log_xaxis(9), b% pg% History_Track_log_yaxis(9), &
         b% pg% show_History_Track_target_box(9), b% pg% History_Track_n_sigma(9), &
         b% pg% History_Track_xtarget(9), b% pg% History_Track_ytarget(9), &
         b% pg% History_Track_xsigma(9), b% pg% History_Track_ysigma(9), &
         b% pg% show_History_Track_annotation1(9), &
         b% pg% show_History_Track_annotation2(9), &
         b% pg% show_History_Track_annotation3(9), &
         b% pg% History_Track_fname(9), &
         b% pg% History_Track_use_decorator(9), b% pg% History_Track9_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_History_Track9_plot


   subroutine null_decorate(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_decorate


   subroutine do_Hist_Track(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
      xname, yname, xaxis_label, yaxis_label, &
      given_xmin, given_xmax, xmargin, dxmin, &
      given_ymin, given_ymax, ymargin, dymin, &
      step_min, step_max, &
      reverse_xaxis, reverse_yaxis, &
      log_xaxis, log_yaxis, &
      show_target_box, n_sigma, &
      xtarget, ytarget, xsigma, ysigma, &
      show_annotation1, show_annotation2, show_annotation3, &
      fname, &
      use_decorator, pgbinary_decorator, &
      decorate, ierr)

      use utils_lib
      use pgstar_support, only : set_xleft_xright, set_ytop_ybot

      type (binary_info), pointer :: b
      integer, intent(in) :: &
         id, device_id, step_min, step_max, n_sigma
      real, intent(in) :: &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale, &
         xtarget, ytarget, xsigma, ysigma, &
         given_xmin, given_xmax, xmargin, dxmin, &
         given_ymin, given_ymax, ymargin, dymin
      character (len = *), intent(in) :: &
         title, xname, yname, xaxis_label, yaxis_label, fname
      logical, intent(in) :: subplot, &
         reverse_xaxis, reverse_yaxis, log_xaxis, log_yaxis, show_target_box, &
         show_annotation1, show_annotation2, show_annotation3, use_decorator
      interface
         subroutine decorate(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine decorate
      end interface
      integer, intent(out) :: ierr
      procedure(pgbinary_decorator_interface), pointer :: pgbinary_decorator

      real :: xmin, xmax, ymin, ymax, xleft, xright, ybot, ytop
      integer :: i, j, j_min, j_max
      real :: dx, dy, xplus, xminus, yplus, yminus
      real, dimension(:), pointer :: xvec, yvec
      character (len = strlen) :: str
      integer :: k, n
      integer :: ix, iy
      integer :: file_data_len
      real, pointer, dimension(:) :: file_data_xvec, file_data_yvec

      logical, parameter :: dbg = .false.

      include 'formats'

      ierr = 0

      call integer_dict_lookup(b% binary_history_names_dict, xname, ix, ierr)
      if (ierr /= 0) ix = -1
      if (ix <= 0) then
         write(*, *)
         write(*, *) 'ERROR: failed to find ' // &
            trim(xname) // ' in history data'
         write(*, *)
         ierr = -1
      end if

      call integer_dict_lookup(b% binary_history_names_dict, yname, iy, ierr)
      if (ierr /= 0) iy = -1
      if (iy <= 0) then
         write(*, *)
         write(*, *) 'ERROR: failed to find ' // &
            trim(yname) // ' in history data'
         write(*, *)
         ierr = -1
      end if
      if (ierr /= 0) return

      n = count_binary_hist_points(b, step_min, step_max)
      allocate(xvec(n), yvec(n), stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'allocate failed for PGSTAR'
         return
      end if

      call get_binary_hist_points(b, step_min, step_max, n, ix, xvec)
      call get_binayr_hist_points(b, step_min, step_max, n, iy, yvec)

      if (log_xaxis) then
         do k = 1, n
            xvec(k) = log10(max(tiny(xvec(k)), abs(xvec(k))))
         end do
      end if

      if (log_yaxis) then
         do k = 1, n
            yvec(k) = log10(max(tiny(yvec(k)), abs(yvec(k))))
         end do
      end if

      call set_xleft_xright(&
         n, xvec, given_xmin, given_xmax, xmargin, &
         reverse_xaxis, dxmin, xleft, xright)

      call set_ytop_ybot(&
         n, yvec, given_ymin, given_ymax, -101.0, ymargin, &
         reverse_yaxis, dymin, ybot, ytop)

      call pgsave
      call pgsch(txt_scale)
      call pgsvp(vp_xleft, vp_xright, vp_ybot, vp_ytop)
      call pgswin(xleft, xright, ybot, ytop)
      call pgscf(1)
      call pgsci(1)
      call show_box_pgbinary(b, 'BCNST1', 'BCNSTV1')

      if (log_xaxis) then
         call show_xaxis_label_pgbinary(b, 'log ' // xaxis_label)
      else
         call show_xaxis_label_pgbinary(b, xaxis_label)
      end if

      if (log_yaxis) then
         call show_left_yaxis_label_pgbinary(b, 'log ' // yaxis_label)
      else
         call show_left_yaxis_label_pgbinary(b, yaxis_label)
      end if

      if (.not. subplot) then
         call show_model_number_pgbinary(b)
         call show_age_pgbinary(b)
      end if
      call show_title_pgbinary(b, title)

      call pgslw(b% pg% pgbinary_lw)

      call show_file_track

      if (show_target_box) then
         call pgsci(clr_Silver)
         if (n_sigma >= 0) then
            j_min = n_sigma
            j_max = n_sigma
         else
            j_min = 1
            j_max = -n_sigma
         end if
         do j = j_min, j_max
            dx = xsigma * j
            xplus = xtarget + dx
            xminus = xtarget - dx
            if (log_xaxis) then
               xplus = log10(max(tiny(xplus), xplus))
               xminus = log10(max(tiny(xminus), xminus))
            end if
            dy = ysigma * j
            yplus = ytarget + dy
            yminus = ytarget - dy
            if (log_yaxis) then
               yplus = log10(max(tiny(yplus), yplus))
               yminus = log10(max(tiny(yminus), yminus))
            end if
            call pgmove(xminus, yminus)
            call pgdraw(xplus, yminus)
            call pgdraw(xplus, yplus)
            call pgdraw(xminus, yplus)
            call pgdraw(xminus, yminus)
         end do
      end if

      call pgsci(clr_Teal)
      call pgline(n, xvec, yvec)
      call pgsci(clr_Crimson)
      call pgsch(2.8 * txt_scale)
      call pgpt1(xvec(n), yvec(n), 0902)

      call show_annotations(b, &
         show_annotation1, &
         show_annotation2, &
         show_annotation3)

      call decorate(b% binary_id, ierr)

      call pgunsa

      call show_pgbinary_decorator(b% binary_id, use_decorator, pgbinary_decorator, 0, ierr)

      deallocate(xvec, yvec)

   contains


      subroutine show_file_track
         integer :: k
         if (len_trim(fname) == 0) return
         if (.not. read_values_from_file(fname, &
            file_data_xvec, file_data_yvec, file_data_len)) then
            write(*, *) &
               'bad filename for History tracks plot ' // trim(fname)
            return
         end if
         if (log_xaxis) then
            do k = 1, file_data_len
               file_data_xvec(k) = log10(max(tiny(file_data_xvec(k)), abs(file_data_xvec(k))))
            end do
         end if
         if (log_yaxis) then
            do k = 1, file_data_len
               file_data_yvec(k) = log10(max(tiny(file_data_yvec(k)), abs(file_data_yvec(k))))
            end do
         end if
         call pgsci(clr_Goldenrod)
         call pgline(file_data_len, file_data_xvec, file_data_yvec)
         deallocate(file_data_xvec, file_data_yvec)
      end subroutine show_file_track


   end subroutine do_Hist_Track


end module pgbinary_hist_track

