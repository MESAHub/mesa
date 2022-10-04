! ***********************************************************************
!
!   Copyright (C) 2013-2022  Bill Paxton, Matthias Fabry
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
         b% History_Track1_xleft, b% History_Track1_xright, &
         b% History_Track1_ybot, b% History_Track1_ytop, .false., &
         b% History_Track1_title, b% History_Track1_txt_scale, ierr)
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
         b% History_Track1_xname, &
         b% History_Track1_yname, &
         b% History_Track1_xaxis_label, &
         b% History_Track1_yaxis_label, &
         b% History_Track1_xmin, &
         b% History_Track1_xmax, &
         b% History_Track1_xmargin, &
         b% History_Track1_dxmin, &
         b% History_Track1_ymin, &
         b% History_Track1_ymax, &
         b% History_Track1_ymargin, &
         b% History_Track1_dymin, &
         b% History_Track1_step_min, &
         b% History_Track1_step_max, &
         b% History_Track1_reverse_xaxis, &
         b% History_Track1_reverse_yaxis, &
         b% History_Track1_log_xaxis, &
         b% History_Track1_log_yaxis, &
         b% show_History_Track1_target_box, &
         b% History_Track1_n_sigma, &
         b% History_Track1_xtarget, &
         b% History_Track1_ytarget, &
         b% History_Track1_xsigma, &
         b% History_Track1_ysigma, &
         b% show_History_Track1_annotation1, &
         b% show_History_Track1_annotation2, &
         b% show_History_Track1_annotation3, &
         b% History_Track1_fname, &
         b% History_Track1_use_decorator, &
         b% History_Track1_pgbinary_decorator, &
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
         b% History_Track2_xleft, b% History_Track2_xright, &
         b% History_Track2_ybot, b% History_Track2_ytop, .false., &
         b% History_Track2_title, b% History_Track2_txt_scale, ierr)
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
         b% History_Track2_xname, &
         b% History_Track2_yname, &
         b% History_Track2_xaxis_label, &
         b% History_Track2_yaxis_label, &
         b% History_Track2_xmin, &
         b% History_Track2_xmax, &
         b% History_Track2_xmargin, &
         b% History_Track2_dxmin, &
         b% History_Track2_ymin, &
         b% History_Track2_ymax, &
         b% History_Track2_ymargin, &
         b% History_Track2_dymin, &
         b% History_Track2_step_min, &
         b% History_Track2_step_max, &
         b% History_Track2_reverse_xaxis, &
         b% History_Track2_reverse_yaxis, &
         b% History_Track2_log_xaxis, &
         b% History_Track2_log_yaxis, &
         b% show_History_Track2_target_box, &
         b% History_Track2_n_sigma, &
         b% History_Track2_xtarget, &
         b% History_Track2_ytarget, &
         b% History_Track2_xsigma, &
         b% History_Track2_ysigma, &
         b% show_History_Track2_annotation1, &
         b% show_History_Track2_annotation2, &
         b% show_History_Track2_annotation3, &
         b% History_Track2_fname, &
         b% History_Track2_use_decorator, &
         b% History_Track2_pgbinary_decorator, &
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
         b% History_Track3_xleft, b% History_Track3_xright, &
         b% History_Track3_ybot, b% History_Track3_ytop, .false., &
         b% History_Track3_title, b% History_Track3_txt_scale, ierr)
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
         b% History_Track3_xname, &
         b% History_Track3_yname, &
         b% History_Track3_xaxis_label, &
         b% History_Track3_yaxis_label, &
         b% History_Track3_xmin, &
         b% History_Track3_xmax, &
         b% History_Track3_xmargin, &
         b% History_Track3_dxmin, &
         b% History_Track3_ymin, &
         b% History_Track3_ymax, &
         b% History_Track3_ymargin, &
         b% History_Track3_dymin, &
         b% History_Track3_step_min, &
         b% History_Track3_step_max, &
         b% History_Track3_reverse_xaxis, &
         b% History_Track3_reverse_yaxis, &
         b% History_Track3_log_xaxis, &
         b% History_Track3_log_yaxis, &
         b% show_History_Track3_target_box, &
         b% History_Track3_n_sigma, &
         b% History_Track3_xtarget, &
         b% History_Track3_ytarget, &
         b% History_Track3_xsigma, &
         b% History_Track3_ysigma, &
         b% show_History_Track3_annotation1, &
         b% show_History_Track3_annotation2, &
         b% show_History_Track3_annotation3, &
         b% History_Track3_fname, &
         b% History_Track3_use_decorator, &
         b% History_Track3_pgbinary_decorator, &
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
         b% History_Track4_xleft, b% History_Track4_xright, &
         b% History_Track4_ybot, b% History_Track4_ytop, .false., &
         b% History_Track4_title, b% History_Track4_txt_scale, ierr)
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
         b% History_Track4_xname, &
         b% History_Track4_yname, &
         b% History_Track4_xaxis_label, &
         b% History_Track4_yaxis_label, &
         b% History_Track4_xmin, &
         b% History_Track4_xmax, &
         b% History_Track4_xmargin, &
         b% History_Track4_dxmin, &
         b% History_Track4_ymin, &
         b% History_Track4_ymax, &
         b% History_Track4_ymargin, &
         b% History_Track4_dymin, &
         b% History_Track4_step_min, &
         b% History_Track4_step_max, &
         b% History_Track4_reverse_xaxis, &
         b% History_Track4_reverse_yaxis, &
         b% History_Track4_log_xaxis, &
         b% History_Track4_log_yaxis, &
         b% show_History_Track4_target_box, &
         b% History_Track4_n_sigma, &
         b% History_Track4_xtarget, &
         b% History_Track4_ytarget, &
         b% History_Track4_xsigma, &
         b% History_Track4_ysigma, &
         b% show_History_Track4_annotation1, &
         b% show_History_Track4_annotation2, &
         b% show_History_Track4_annotation3, &
         b% History_Track4_fname, &
         b% History_Track4_use_decorator, &
         b% History_Track4_pgbinary_decorator, &
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
         b% History_Track5_xleft, b% History_Track5_xright, &
         b% History_Track5_ybot, b% History_Track5_ytop, .false., &
         b% History_Track5_title, b% History_Track5_txt_scale, ierr)
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
         b% History_Track5_xname, &
         b% History_Track5_yname, &
         b% History_Track5_xaxis_label, &
         b% History_Track5_yaxis_label, &
         b% History_Track5_xmin, &
         b% History_Track5_xmax, &
         b% History_Track5_xmargin, &
         b% History_Track5_dxmin, &
         b% History_Track5_ymin, &
         b% History_Track5_ymax, &
         b% History_Track5_ymargin, &
         b% History_Track5_dymin, &
         b% History_Track5_step_min, &
         b% History_Track5_step_max, &
         b% History_Track5_reverse_xaxis, &
         b% History_Track5_reverse_yaxis, &
         b% History_Track5_log_xaxis, &
         b% History_Track5_log_yaxis, &
         b% show_History_Track5_target_box, &
         b% History_Track5_n_sigma, &
         b% History_Track5_xtarget, &
         b% History_Track5_ytarget, &
         b% History_Track5_xsigma, &
         b% History_Track5_ysigma, &
         b% show_History_Track5_annotation1, &
         b% show_History_Track5_annotation2, &
         b% show_History_Track5_annotation3, &
         b% History_Track5_fname, &
         b% History_Track5_use_decorator, &
         b% History_Track5_pgbinary_decorator, &
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
         b% History_Track6_xleft, b% History_Track6_xright, &
         b% History_Track6_ybot, b% History_Track6_ytop, .false., &
         b% History_Track6_title, b% History_Track6_txt_scale, ierr)
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
         b% History_Track6_xname, &
         b% History_Track6_yname, &
         b% History_Track6_xaxis_label, &
         b% History_Track6_yaxis_label, &
         b% History_Track6_xmin, &
         b% History_Track6_xmax, &
         b% History_Track6_xmargin, &
         b% History_Track6_dxmin, &
         b% History_Track6_ymin, &
         b% History_Track6_ymax, &
         b% History_Track6_ymargin, &
         b% History_Track6_dymin, &
         b% History_Track6_step_min, &
         b% History_Track6_step_max, &
         b% History_Track6_reverse_xaxis, &
         b% History_Track6_reverse_yaxis, &
         b% History_Track6_log_xaxis, &
         b% History_Track6_log_yaxis, &
         b% show_History_Track6_target_box, &
         b% History_Track6_n_sigma, &
         b% History_Track6_xtarget, &
         b% History_Track6_ytarget, &
         b% History_Track6_xsigma, &
         b% History_Track6_ysigma, &
         b% show_History_Track6_annotation1, &
         b% show_History_Track6_annotation2, &
         b% show_History_Track6_annotation3, &
         b% History_Track6_fname, &
         b% History_Track6_use_decorator, &
         b% History_Track6_pgbinary_decorator, &
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
         b% History_Track7_xleft, b% History_Track7_xright, &
         b% History_Track7_ybot, b% History_Track7_ytop, .false., &
         b% History_Track7_title, b% History_Track7_txt_scale, ierr)
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
         b% History_Track7_xname, &
         b% History_Track7_yname, &
         b% History_Track7_xaxis_label, &
         b% History_Track7_yaxis_label, &
         b% History_Track7_xmin, &
         b% History_Track7_xmax, &
         b% History_Track7_xmargin, &
         b% History_Track7_dxmin, &
         b% History_Track7_ymin, &
         b% History_Track7_ymax, &
         b% History_Track7_ymargin, &
         b% History_Track7_dymin, &
         b% History_Track7_step_min, &
         b% History_Track7_step_max, &
         b% History_Track7_reverse_xaxis, &
         b% History_Track7_reverse_yaxis, &
         b% History_Track7_log_xaxis, &
         b% History_Track7_log_yaxis, &
         b% show_History_Track7_target_box, &
         b% History_Track7_n_sigma, &
         b% History_Track7_xtarget, &
         b% History_Track7_ytarget, &
         b% History_Track7_xsigma, &
         b% History_Track7_ysigma, &
         b% show_History_Track7_annotation1, &
         b% show_History_Track7_annotation2, &
         b% show_History_Track7_annotation3, &
         b% History_Track7_fname, &
         b% History_Track7_use_decorator, &
         b% History_Track7_pgbinary_decorator, &
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
         b% History_Track8_xleft, b% History_Track8_xright, &
         b% History_Track8_ybot, b% History_Track8_ytop, .false., &
         b% History_Track8_title, b% History_Track8_txt_scale, ierr)
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
         b% History_Track8_xname, &
         b% History_Track8_yname, &
         b% History_Track8_xaxis_label, &
         b% History_Track8_yaxis_label, &
         b% History_Track8_xmin, &
         b% History_Track8_xmax, &
         b% History_Track8_xmargin, &
         b% History_Track8_dxmin, &
         b% History_Track8_ymin, &
         b% History_Track8_ymax, &
         b% History_Track8_ymargin, &
         b% History_Track8_dymin, &
         b% History_Track8_step_min, &
         b% History_Track8_step_max, &
         b% History_Track8_reverse_xaxis, &
         b% History_Track8_reverse_yaxis, &
         b% History_Track8_log_xaxis, &
         b% History_Track8_log_yaxis, &
         b% show_History_Track8_target_box, &
         b% History_Track8_n_sigma, &
         b% History_Track8_xtarget, &
         b% History_Track8_ytarget, &
         b% History_Track8_xsigma, &
         b% History_Track8_ysigma, &
         b% show_History_Track8_annotation1, &
         b% show_History_Track8_annotation2, &
         b% show_History_Track8_annotation3, &
         b% History_Track8_fname, &
         b% History_Track8_use_decorator, &
         b% History_Track8_pgbinary_decorator, &
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
         b% History_Track9_xleft, b% History_Track9_xright, &
         b% History_Track9_ybot, b% History_Track9_ytop, .false., &
         b% History_Track9_title, b% History_Track9_txt_scale, ierr)
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
         b% History_Track9_xname, &
         b% History_Track9_yname, &
         b% History_Track9_xaxis_label, &
         b% History_Track9_yaxis_label, &
         b% History_Track9_xmin, &
         b% History_Track9_xmax, &
         b% History_Track9_xmargin, &
         b% History_Track9_dxmin, &
         b% History_Track9_ymin, &
         b% History_Track9_ymax, &
         b% History_Track9_ymargin, &
         b% History_Track9_dymin, &
         b% History_Track9_step_min, &
         b% History_Track9_step_max, &
         b% History_Track9_reverse_xaxis, &
         b% History_Track9_reverse_yaxis, &
         b% History_Track9_log_xaxis, &
         b% History_Track9_log_yaxis, &
         b% show_History_Track9_target_box, &
         b% History_Track9_n_sigma, &
         b% History_Track9_xtarget, &
         b% History_Track9_ytarget, &
         b% History_Track9_xsigma, &
         b% History_Track9_ysigma, &
         b% show_History_Track9_annotation1, &
         b% show_History_Track9_annotation2, &
         b% show_History_Track9_annotation3, &
         b% History_Track9_fname, &
         b% History_Track9_use_decorator, &
         b% History_Track9_pgbinary_decorator, &
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

      n = count_hist_points(b, step_min, step_max)
      allocate(xvec(n), yvec(n), stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'allocate failed for PGSTAR'
         return
      end if

      call get_hist_points(b, step_min, step_max, n, ix, xvec)
      call get_hist_points(b, step_min, step_max, n, iy, yvec)

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

      call pgslw(b% pgbinary_lw)

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

