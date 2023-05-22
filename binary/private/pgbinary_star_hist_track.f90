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

module pgbinary_star_hist_track

   use binary_private_def
   use const_def
   use pgbinary_support

   implicit none


contains


   subroutine Star_history_track1_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star_history_track1_plot(b, id, device_id, &
         b% pg% Star_history_track1_xleft, b% pg% Star_history_track1_xright, &
         b% pg% Star_history_track1_ybot, b% pg% Star_history_track1_ytop, .false., &
         b% pg% Star_history_track1_title, b% pg% Star_history_track1_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Star_history_track1_plot


   subroutine do_Star_history_track1_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Star_hist_track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% Star_history_track1_xname, &
         b% pg% Star_history_track1_yname, &
         b% pg% Star_history_track1_xaxis_label, &
         b% pg% Star_history_track1_yaxis_label, &
         b% pg% Star_history_track1_xmin, &
         b% pg% Star_history_track1_xmax, &
         b% pg% Star_history_track1_xmargin, &
         b% pg% Star_history_track1_dxmin, &
         b% pg% Star_history_track1_ymin, &
         b% pg% Star_history_track1_ymax, &
         b% pg% Star_history_track1_ymargin, &
         b% pg% Star_history_track1_dymin, &
         b% pg% Star_history_track1_step_min, &
         b% pg% Star_history_track1_step_max, &
         b% pg% Star_history_track1_reverse_xaxis, &
         b% pg% Star_history_track1_reverse_yaxis, &
         b% pg% Star_history_track1_log_xaxis, &
         b% pg% Star_history_track1_log_yaxis, &
         b% pg% show_Star_history_track1_annotation1, &
         b% pg% show_Star_history_track1_annotation2, &
         b% pg% show_Star_history_track1_annotation3, &
         b% pg% Star_history_track1_fname, &
         b% pg% Star_history_track1_use_decorator, &
         b% pg% Star_history_track1_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_Star_history_track1_plot


   subroutine Star_history_track2_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star_history_track2_plot(b, id, device_id, &
         b% pg% Star_history_track2_xleft, b% pg% Star_history_track2_xright, &
         b% pg% Star_history_track2_ybot, b% pg% Star_history_track2_ytop, .false., &
         b% pg% Star_history_track2_title, b% pg% Star_history_track2_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Star_history_track2_plot


   subroutine do_Star_history_track2_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Star_hist_track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% Star_history_track2_xname, &
         b% pg% Star_history_track2_yname, &
         b% pg% Star_history_track2_xaxis_label, &
         b% pg% Star_history_track2_yaxis_label, &
         b% pg% Star_history_track2_xmin, &
         b% pg% Star_history_track2_xmax, &
         b% pg% Star_history_track2_xmargin, &
         b% pg% Star_history_track2_dxmin, &
         b% pg% Star_history_track2_ymin, &
         b% pg% Star_history_track2_ymax, &
         b% pg% Star_history_track2_ymargin, &
         b% pg% Star_history_track2_dymin, &
         b% pg% Star_history_track2_step_min, &
         b% pg% Star_history_track2_step_max, &
         b% pg% Star_history_track2_reverse_xaxis, &
         b% pg% Star_history_track2_reverse_yaxis, &
         b% pg% Star_history_track2_log_xaxis, &
         b% pg% Star_history_track2_log_yaxis, &
         b% pg% show_Star_history_track2_annotation1, &
         b% pg% show_Star_history_track2_annotation2, &
         b% pg% show_Star_history_track2_annotation3, &
         b% pg% Star_history_track2_fname, &
         b% pg% Star_history_track2_use_decorator, &
         b% pg% Star_history_track2_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_Star_history_track2_plot


   subroutine Star_history_track3_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star_history_track3_plot(b, id, device_id, &
         b% pg% Star_history_track3_xleft, b% pg% Star_history_track3_xright, &
         b% pg% Star_history_track3_ybot, b% pg% Star_history_track3_ytop, .false., &
         b% pg% Star_history_track3_title, b% pg% Star_history_track3_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Star_history_track3_plot


   subroutine do_Star_history_track3_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Star_hist_track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% Star_history_track3_xname, &
         b% pg% Star_history_track3_yname, &
         b% pg% Star_history_track3_xaxis_label, &
         b% pg% Star_history_track3_yaxis_label, &
         b% pg% Star_history_track3_xmin, &
         b% pg% Star_history_track3_xmax, &
         b% pg% Star_history_track3_xmargin, &
         b% pg% Star_history_track3_dxmin, &
         b% pg% Star_history_track3_ymin, &
         b% pg% Star_history_track3_ymax, &
         b% pg% Star_history_track3_ymargin, &
         b% pg% Star_history_track3_dymin, &
         b% pg% Star_history_track3_step_min, &
         b% pg% Star_history_track3_step_max, &
         b% pg% Star_history_track3_reverse_xaxis, &
         b% pg% Star_history_track3_reverse_yaxis, &
         b% pg% Star_history_track3_log_xaxis, &
         b% pg% Star_history_track3_log_yaxis, &
         b% pg% show_Star_history_track3_annotation1, &
         b% pg% show_Star_history_track3_annotation2, &
         b% pg% show_Star_history_track3_annotation3, &
         b% pg% Star_history_track3_fname, &
         b% pg% Star_history_track3_use_decorator, &
         b% pg% Star_history_track3_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_Star_history_track3_plot


   subroutine Star_history_track4_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star_history_track4_plot(b, id, device_id, &
         b% pg% Star_history_track4_xleft, b% pg% Star_history_track4_xright, &
         b% pg% Star_history_track4_ybot, b% pg% Star_history_track4_ytop, .false., &
         b% pg% Star_history_track4_title, b% pg% Star_history_track4_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Star_history_track4_plot


   subroutine do_Star_history_track4_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Star_hist_track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% Star_history_track4_xname, &
         b% pg% Star_history_track4_yname, &
         b% pg% Star_history_track4_xaxis_label, &
         b% pg% Star_history_track4_yaxis_label, &
         b% pg% Star_history_track4_xmin, &
         b% pg% Star_history_track4_xmax, &
         b% pg% Star_history_track4_xmargin, &
         b% pg% Star_history_track4_dxmin, &
         b% pg% Star_history_track4_ymin, &
         b% pg% Star_history_track4_ymax, &
         b% pg% Star_history_track4_ymargin, &
         b% pg% Star_history_track4_dymin, &
         b% pg% Star_history_track4_step_min, &
         b% pg% Star_history_track4_step_max, &
         b% pg% Star_history_track4_reverse_xaxis, &
         b% pg% Star_history_track4_reverse_yaxis, &
         b% pg% Star_history_track4_log_xaxis, &
         b% pg% Star_history_track4_log_yaxis, &
         b% pg% show_Star_history_track4_annotation1, &
         b% pg% show_Star_history_track4_annotation2, &
         b% pg% show_Star_history_track4_annotation3, &
         b% pg% Star_history_track4_fname, &
         b% pg% Star_history_track4_use_decorator, &
         b% pg% Star_history_track4_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_Star_history_track4_plot


   subroutine Star_history_track5_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star_history_track5_plot(b, id, device_id, &
         b% pg% Star_history_track5_xleft, b% pg% Star_history_track5_xright, &
         b% pg% Star_history_track5_ybot, b% pg% Star_history_track5_ytop, .false., &
         b% pg% Star_history_track5_title, b% pg% Star_history_track5_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Star_history_track5_plot


   subroutine do_Star_history_track5_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Star_hist_track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% Star_history_track5_xname, &
         b% pg% Star_history_track5_yname, &
         b% pg% Star_history_track5_xaxis_label, &
         b% pg% Star_history_track5_yaxis_label, &
         b% pg% Star_history_track5_xmin, &
         b% pg% Star_history_track5_xmax, &
         b% pg% Star_history_track5_xmargin, &
         b% pg% Star_history_track5_dxmin, &
         b% pg% Star_history_track5_ymin, &
         b% pg% Star_history_track5_ymax, &
         b% pg% Star_history_track5_ymargin, &
         b% pg% Star_history_track5_dymin, &
         b% pg% Star_history_track5_step_min, &
         b% pg% Star_history_track5_step_max, &
         b% pg% Star_history_track5_reverse_xaxis, &
         b% pg% Star_history_track5_reverse_yaxis, &
         b% pg% Star_history_track5_log_xaxis, &
         b% pg% Star_history_track5_log_yaxis, &
         b% pg% show_Star_history_track5_annotation1, &
         b% pg% show_Star_history_track5_annotation2, &
         b% pg% show_Star_history_track5_annotation3, &
         b% pg% Star_history_track5_fname, &
         b% pg% Star_history_track5_use_decorator, &
         b% pg% Star_history_track5_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_Star_history_track5_plot


   subroutine Star_history_track6_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star_history_track6_plot(b, id, device_id, &
         b% pg% Star_history_track6_xleft, b% pg% Star_history_track6_xright, &
         b% pg% Star_history_track6_ybot, b% pg% Star_history_track6_ytop, .false., &
         b% pg% Star_history_track6_title, b% pg% Star_history_track6_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Star_history_track6_plot


   subroutine do_Star_history_track6_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Star_hist_track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% Star_history_track6_xname, &
         b% pg% Star_history_track6_yname, &
         b% pg% Star_history_track6_xaxis_label, &
         b% pg% Star_history_track6_yaxis_label, &
         b% pg% Star_history_track6_xmin, &
         b% pg% Star_history_track6_xmax, &
         b% pg% Star_history_track6_xmargin, &
         b% pg% Star_history_track6_dxmin, &
         b% pg% Star_history_track6_ymin, &
         b% pg% Star_history_track6_ymax, &
         b% pg% Star_history_track6_ymargin, &
         b% pg% Star_history_track6_dymin, &
         b% pg% Star_history_track6_step_min, &
         b% pg% Star_history_track6_step_max, &
         b% pg% Star_history_track6_reverse_xaxis, &
         b% pg% Star_history_track6_reverse_yaxis, &
         b% pg% Star_history_track6_log_xaxis, &
         b% pg% Star_history_track6_log_yaxis, &
         b% pg% show_Star_history_track6_annotation1, &
         b% pg% show_Star_history_track6_annotation2, &
         b% pg% show_Star_history_track6_annotation3, &
         b% pg% Star_history_track6_fname, &
         b% pg% Star_history_track6_use_decorator, &
         b% pg% Star_history_track6_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_Star_history_track6_plot


   subroutine Star_history_track7_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star_history_track7_plot(b, id, device_id, &
         b% pg% Star_history_track7_xleft, b% pg% Star_history_track7_xright, &
         b% pg% Star_history_track7_ybot, b% pg% Star_history_track7_ytop, .false., &
         b% pg% Star_history_track7_title, b% pg% Star_history_track7_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Star_history_track7_plot


   subroutine do_Star_history_track7_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Star_hist_track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% Star_history_track7_xname, &
         b% pg% Star_history_track7_yname, &
         b% pg% Star_history_track7_xaxis_label, &
         b% pg% Star_history_track7_yaxis_label, &
         b% pg% Star_history_track7_xmin, &
         b% pg% Star_history_track7_xmax, &
         b% pg% Star_history_track7_xmargin, &
         b% pg% Star_history_track7_dxmin, &
         b% pg% Star_history_track7_ymin, &
         b% pg% Star_history_track7_ymax, &
         b% pg% Star_history_track7_ymargin, &
         b% pg% Star_history_track7_dymin, &
         b% pg% Star_history_track7_step_min, &
         b% pg% Star_history_track7_step_max, &
         b% pg% Star_history_track7_reverse_xaxis, &
         b% pg% Star_history_track7_reverse_yaxis, &
         b% pg% Star_history_track7_log_xaxis, &
         b% pg% Star_history_track7_log_yaxis, &
         b% pg% show_Star_history_track7_annotation1, &
         b% pg% show_Star_history_track7_annotation2, &
         b% pg% show_Star_history_track7_annotation3, &
         b% pg% Star_history_track7_fname, &
         b% pg% Star_history_track7_use_decorator, &
         b% pg% Star_history_track7_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_Star_history_track7_plot


   subroutine Star_history_track8_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star_history_track8_plot(b, id, device_id, &
         b% pg% Star_history_track8_xleft, b% pg% Star_history_track8_xright, &
         b% pg% Star_history_track8_ybot, b% pg% Star_history_track8_ytop, .false., &
         b% pg% Star_history_track8_title, b% pg% Star_history_track8_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Star_history_track8_plot


   subroutine do_Star_history_track8_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Star_hist_track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% Star_history_track8_xname, &
         b% pg% Star_history_track8_yname, &
         b% pg% Star_history_track8_xaxis_label, &
         b% pg% Star_history_track8_yaxis_label, &
         b% pg% Star_history_track8_xmin, &
         b% pg% Star_history_track8_xmax, &
         b% pg% Star_history_track8_xmargin, &
         b% pg% Star_history_track8_dxmin, &
         b% pg% Star_history_track8_ymin, &
         b% pg% Star_history_track8_ymax, &
         b% pg% Star_history_track8_ymargin, &
         b% pg% Star_history_track8_dymin, &
         b% pg% Star_history_track8_step_min, &
         b% pg% Star_history_track8_step_max, &
         b% pg% Star_history_track8_reverse_xaxis, &
         b% pg% Star_history_track8_reverse_yaxis, &
         b% pg% Star_history_track8_log_xaxis, &
         b% pg% Star_history_track8_log_yaxis, &
         b% pg% show_Star_history_track8_annotation1, &
         b% pg% show_Star_history_track8_annotation2, &
         b% pg% show_Star_history_track8_annotation3, &
         b% pg% Star_history_track8_fname, &
         b% pg% Star_history_track8_use_decorator, &
         b% pg% Star_history_track8_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_Star_history_track8_plot


   subroutine Star_history_track9_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Star_history_track9_plot(b, id, device_id, &
         b% pg% Star_history_track9_xleft, b% pg% Star_history_track9_xright, &
         b% pg% Star_history_track9_ybot, b% pg% Star_history_track9_ytop, .false., &
         b% pg% Star_history_track9_title, b% pg% Star_history_track9_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Star_history_track9_plot


   subroutine do_Star_history_track9_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Star_hist_track(b, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% Star_history_track9_xname, &
         b% pg% Star_history_track9_yname, &
         b% pg% Star_history_track9_xaxis_label, &
         b% pg% Star_history_track9_yaxis_label, &
         b% pg% Star_history_track9_xmin, &
         b% pg% Star_history_track9_xmax, &
         b% pg% Star_history_track9_xmargin, &
         b% pg% Star_history_track9_dxmin, &
         b% pg% Star_history_track9_ymin, &
         b% pg% Star_history_track9_ymax, &
         b% pg% Star_history_track9_ymargin, &
         b% pg% Star_history_track9_dymin, &
         b% pg% Star_history_track9_step_min, &
         b% pg% Star_history_track9_step_max, &
         b% pg% Star_history_track9_reverse_xaxis, &
         b% pg% Star_history_track9_reverse_yaxis, &
         b% pg% Star_history_track9_log_xaxis, &
         b% pg% Star_history_track9_log_yaxis, &
         b% pg% show_Star_history_track9_annotation1, &
         b% pg% show_Star_history_track9_annotation2, &
         b% pg% show_Star_history_track9_annotation3, &
         b% pg% Star_history_track9_fname, &
         b% pg% Star_history_track9_use_decorator, &
         b% pg% Star_history_track9_pgbinary_decorator, &
         null_decorate, ierr)
   end subroutine do_Star_history_track9_plot


   subroutine null_decorate(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_decorate


   subroutine do_Star_hist_track(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
      xname, yname, xaxis_label, yaxis_label, &
      given_xmin, given_xmax, xmargin, dxmin, &
      given_ymin, given_ymax, ymargin, dymin, &
      step_min, step_max, &
      reverse_xaxis, reverse_yaxis, &
      log_xaxis, log_yaxis, &
      show_annotation1, show_annotation2, show_annotation3, &
      fname, &
      use_decorator, pgbinary_decorator, &
      decorate, ierr)

      use utils_lib
      use pgstar_support, only: set_xleft_xright, set_ytop_ybot

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
      integer :: i, j, j_min, j_max, k, n, star, ix, iy, file_data_len
      real :: dx, dy, xplus, xminus, yplus, yminus
      real, dimension(:), pointer :: xvec, yvec
      character (len = strlen) :: str
      real, pointer, dimension(:) :: file_data_xvec, file_data_yvec
      type (star_info), pointer :: s

      logical, parameter :: dbg = .false.

      include 'formats'

      ierr = 0

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

      do star = 1, 2
         if (star == 1) then
            if (.not. b% have_star_1) then
               cycle
            else
               s => b% s1
            end if
         else if (star == 2) then
            if (.not. b% have_star_2) then
               cycle
            else
               s => b% s2
            end if
         end if

         call integer_dict_lookup(s% history_names_dict, xname, ix, ierr)
         if (ierr /= 0) ix = -1
         if (ix <= 0) then
            write(*, *)
            write(*, *) 'ERROR: failed to find ' // &
               trim(xname) // ' in history data'
            write(*, *)
            ierr = -1
         end if

         call integer_dict_lookup(s% history_names_dict, yname, iy, ierr)
         if (ierr /= 0) iy = -1
         if (iy <= 0) then
            write(*, *)
            write(*, *) 'ERROR: failed to find ' // &
               trim(yname) // ' in history data'
            write(*, *)
            ierr = -1
         end if
         if (ierr /= 0) return

         n = count_hist_points(s, step_min, step_max)
         allocate(xvec(n), yvec(n), stat = ierr)
         if (ierr /= 0) then
            write(*, *) 'allocate failed for PGSTAR'
            return
         end if

         call get_hist_points(s, step_min, step_max, n, ix, xvec)
         call get_hist_points(s, step_min, step_max, n, iy, yvec)

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

         if (star == 1)
            call pgsci(b% pg% star_1_color)
         else
            call pgsci(b% pg% star_2_color)
         end if
         call pgline(n, xvec, yvec)  ! history track
         call pgsci(clr_Crimson)
         call pgsch(2.8 * txt_scale)
         call pgpt1(xvec(n), yvec(n), 0902)  ! current point highlighted

         deallocate(xvec, yvec)  ! size could potentially be different, so to be safe alloc again for star 2

      end do

      call pgunsa

      call show_annotations(b, &
         show_annotation1, &
         show_annotation2, &
         show_annotation3)

      call decorate(b% binary_id, ierr)

      call show_pgbinary_decorator(b% binary_id, use_decorator, pgbinary_decorator, 0, ierr)

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

   end subroutine do_Star_hist_track


end module pgbinary_star_hist_track

