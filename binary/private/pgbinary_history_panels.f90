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

module pgbinary_history_panels

   use binary_private_def
   use const_def
   use pgbinary_support

   implicit none


contains


   subroutine History_Panels1_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels1_plot(b, id, device_id, &
         b% pg% History_Panels1_xleft, b% pg% History_Panels1_xright, &
         b% pg% History_Panels1_ybot, b% pg% History_Panels1_ytop, .false., &
         b% pg% History_Panels1_title, b% pg% History_Panels1_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels1_plot


   subroutine do_History_Panels1_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_history_panels_plot(&
         id, b, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Panels1_xaxis_name, &
         b% pg% History_Panels1_xmin, &
         b% pg% History_Panels1_xmax, &
         b% pg% History_Panels1_dxmin, &
         b% pg% History_Panels1_xmargin, &
         b% pg% History_Panels1_max_width, &
         b% pg% History_Panels1_num_panels, &
         b% pg% History_Panels1_other_ymin, &
         b% pg% History_Panels1_other_ymax, &
         b% pg% History_Panels1_yaxis_reversed, &
         b% pg% History_Panels1_other_yaxis_log, &
         b% pg% History_Panels1_other_dymin, &
         b% pg% History_Panels1_points_name, &
         b% pg% History_Panels1_other_ymargin, &
         b% pg% History_Panels1_other_yaxis_name, &
         b% pg% History_Panels1_ymin, &
         b% pg% History_Panels1_ymax, &
         b% pg% History_Panels1_xaxis_reversed, &
         b% pg% History_Panels1_yaxis_reversed, &
         b% pg% History_Panels1_xaxis_log, &
         b% pg% History_Panels1_yaxis_log, &
         b% pg% History_Panels1_dymin, &
         b% pg% History_Panels1_ymargin, &
         b% pg% History_Panels1_yaxis_name, &
         b% pg% History_Panels1_use_decorator, &
         b% pg% History_Panels1_pgbinary_decorator, &
         ierr)
   end subroutine do_History_Panels1_plot


   subroutine History_Panels2_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels2_plot(b, id, device_id, &
         b% pg% History_Panels2_xleft, b% pg% History_Panels2_xright, &
         b% pg% History_Panels2_ybot, b% pg% History_Panels2_ytop, .false., &
         b% pg% History_Panels2_title, b% pg% History_Panels2_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels2_plot


   subroutine do_History_Panels2_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_history_panels_plot(&
         id, b, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Panels2_xaxis_name, &
         b% pg% History_Panels2_xmin, &
         b% pg% History_Panels2_xmax, &
         b% pg% History_Panels2_dxmin, &
         b% pg% History_Panels2_xmargin, &
         b% pg% History_Panels2_max_width, &
         b% pg% History_Panels2_num_panels, &
         b% pg% History_Panels2_other_ymin, &
         b% pg% History_Panels2_other_ymax, &
         b% pg% History_Panels2_yaxis_reversed, &
         b% pg% History_Panels2_other_yaxis_log, &
         b% pg% History_Panels2_other_dymin, &
         b% pg% History_Panels2_points_name, &
         b% pg% History_Panels2_other_ymargin, &
         b% pg% History_Panels2_other_yaxis_name, &
         b% pg% History_Panels2_ymin, &
         b% pg% History_Panels2_ymax, &
         b% pg% History_Panels2_xaxis_reversed, &
         b% pg% History_Panels2_yaxis_reversed, &
         b% pg% History_Panels2_xaxis_log, &
         b% pg% History_Panels2_yaxis_log, &
         b% pg% History_Panels2_dymin, &
         b% pg% History_Panels2_ymargin, &
         b% pg% History_Panels2_yaxis_name, &
         b% pg% History_Panels2_use_decorator, &
         b% pg% History_Panels2_pgbinary_decorator, &
         ierr)
   end subroutine do_History_Panels2_plot


   subroutine History_Panels3_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels3_plot(b, id, device_id, &
         b% pg% History_Panels3_xleft, b% pg% History_Panels3_xright, &
         b% pg% History_Panels3_ybot, b% pg% History_Panels3_ytop, .false., &
         b% pg% History_Panels3_title, b% pg% History_Panels3_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels3_plot


   subroutine do_History_Panels3_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_history_panels_plot(&
         id, b, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Panels3_xaxis_name, &
         b% pg% History_Panels3_xmin, &
         b% pg% History_Panels3_xmax, &
         b% pg% History_Panels3_dxmin, &
         b% pg% History_Panels3_xmargin, &
         b% pg% History_Panels3_max_width, &
         b% pg% History_Panels3_num_panels, &
         b% pg% History_Panels3_other_ymin, &
         b% pg% History_Panels3_other_ymax, &
         b% pg% History_Panels3_yaxis_reversed, &
         b% pg% History_Panels3_other_yaxis_log, &
         b% pg% History_Panels3_other_dymin, &
         b% pg% History_Panels3_points_name, &
         b% pg% History_Panels3_other_ymargin, &
         b% pg% History_Panels3_other_yaxis_name, &
         b% pg% History_Panels3_ymin, &
         b% pg% History_Panels3_ymax, &
         b% pg% History_Panels3_xaxis_reversed, &
         b% pg% History_Panels3_yaxis_reversed, &
         b% pg% History_Panels3_xaxis_log, &
         b% pg% History_Panels3_yaxis_log, &
         b% pg% History_Panels3_dymin, &
         b% pg% History_Panels3_ymargin, &
         b% pg% History_Panels3_yaxis_name, &
         b% pg% History_Panels3_use_decorator, &
         b% pg% History_Panels3_pgbinary_decorator, &
         ierr)
   end subroutine do_History_Panels3_plot


   subroutine History_Panels4_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels4_plot(b, id, device_id, &
         b% pg% History_Panels4_xleft, b% pg% History_Panels4_xright, &
         b% pg% History_Panels4_ybot, b% pg% History_Panels4_ytop, .false., &
         b% pg% History_Panels4_title, b% pg% History_Panels4_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels4_plot


   subroutine do_History_Panels4_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_history_panels_plot(&
         id, b, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Panels4_xaxis_name, &
         b% pg% History_Panels4_xmin, &
         b% pg% History_Panels4_xmax, &
         b% pg% History_Panels4_dxmin, &
         b% pg% History_Panels4_xmargin, &
         b% pg% History_Panels4_max_width, &
         b% pg% History_Panels4_num_panels, &
         b% pg% History_Panels4_other_ymin, &
         b% pg% History_Panels4_other_ymax, &
         b% pg% History_Panels4_yaxis_reversed, &
         b% pg% History_Panels4_other_yaxis_log, &
         b% pg% History_Panels4_other_dymin, &
         b% pg% History_Panels4_points_name, &
         b% pg% History_Panels4_other_ymargin, &
         b% pg% History_Panels4_other_yaxis_name, &
         b% pg% History_Panels4_ymin, &
         b% pg% History_Panels4_ymax, &
         b% pg% History_Panels4_xaxis_reversed, &
         b% pg% History_Panels4_yaxis_reversed, &
         b% pg% History_Panels4_xaxis_log, &
         b% pg% History_Panels4_yaxis_log, &
         b% pg% History_Panels4_dymin, &
         b% pg% History_Panels4_ymargin, &
         b% pg% History_Panels4_yaxis_name, &
         b% pg% History_Panels4_use_decorator, &
         b% pg% History_Panels4_pgbinary_decorator, &
         ierr)
   end subroutine do_History_Panels4_plot


   subroutine History_Panels5_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels5_plot(b, id, device_id, &
         b% pg% History_Panels5_xleft, b% pg% History_Panels5_xright, &
         b% pg% History_Panels5_ybot, b% pg% History_Panels5_ytop, .false., &
         b% pg% History_Panels5_title, b% pg% History_Panels5_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels5_plot


   subroutine do_History_Panels5_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_history_panels_plot(&
         id, b, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Panels5_xaxis_name, &
         b% pg% History_Panels5_xmin, &
         b% pg% History_Panels5_xmax, &
         b% pg% History_Panels5_dxmin, &
         b% pg% History_Panels5_xmargin, &
         b% pg% History_Panels5_max_width, &
         b% pg% History_Panels5_num_panels, &
         b% pg% History_Panels5_other_ymin, &
         b% pg% History_Panels5_other_ymax, &
         b% pg% History_Panels5_yaxis_reversed, &
         b% pg% History_Panels5_other_yaxis_log, &
         b% pg% History_Panels5_other_dymin, &
         b% pg% History_Panels5_points_name, &
         b% pg% History_Panels5_other_ymargin, &
         b% pg% History_Panels5_other_yaxis_name, &
         b% pg% History_Panels5_ymin, &
         b% pg% History_Panels5_ymax, &
         b% pg% History_Panels5_xaxis_reversed, &
         b% pg% History_Panels5_yaxis_reversed, &
         b% pg% History_Panels5_xaxis_log, &
         b% pg% History_Panels5_yaxis_log, &
         b% pg% History_Panels5_dymin, &
         b% pg% History_Panels5_ymargin, &
         b% pg% History_Panels5_yaxis_name, &
         b% pg% History_Panels5_use_decorator, &
         b% pg% History_Panels5_pgbinary_decorator, &
         ierr)
   end subroutine do_History_Panels5_plot


   subroutine History_Panels6_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels6_plot(b, id, device_id, &
         b% pg% History_Panels6_xleft, b% pg% History_Panels6_xright, &
         b% pg% History_Panels6_ybot, b% pg% History_Panels6_ytop, .false., &
         b% pg% History_Panels6_title, b% pg% History_Panels6_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels6_plot


   subroutine do_History_Panels6_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_history_panels_plot(&
         id, b, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Panels6_xaxis_name, &
         b% pg% History_Panels6_xmin, &
         b% pg% History_Panels6_xmax, &
         b% pg% History_Panels6_dxmin, &
         b% pg% History_Panels6_xmargin, &
         b% pg% History_Panels6_max_width, &
         b% pg% History_Panels6_num_panels, &
         b% pg% History_Panels6_other_ymin, &
         b% pg% History_Panels6_other_ymax, &
         b% pg% History_Panels6_yaxis_reversed, &
         b% pg% History_Panels6_other_yaxis_log, &
         b% pg% History_Panels6_other_dymin, &
         b% pg% History_Panels6_points_name, &
         b% pg% History_Panels6_other_ymargin, &
         b% pg% History_Panels6_other_yaxis_name, &
         b% pg% History_Panels6_ymin, &
         b% pg% History_Panels6_ymax, &
         b% pg% History_Panels6_xaxis_reversed, &
         b% pg% History_Panels6_yaxis_reversed, &
         b% pg% History_Panels6_xaxis_log, &
         b% pg% History_Panels6_yaxis_log, &
         b% pg% History_Panels6_dymin, &
         b% pg% History_Panels6_ymargin, &
         b% pg% History_Panels6_yaxis_name, &
         b% pg% History_Panels6_use_decorator, &
         b% pg% History_Panels6_pgbinary_decorator, &
         ierr)
   end subroutine do_History_Panels6_plot


   subroutine History_Panels7_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels7_plot(b, id, device_id, &
         b% pg% History_Panels7_xleft, b% pg% History_Panels7_xright, &
         b% pg% History_Panels7_ybot, b% pg% History_Panels7_ytop, .false., &
         b% pg% History_Panels7_title, b% pg% History_Panels7_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels7_plot


   subroutine do_History_Panels7_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_history_panels_plot(&
         id, b, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Panels7_xaxis_name, &
         b% pg% History_Panels7_xmin, &
         b% pg% History_Panels7_xmax, &
         b% pg% History_Panels7_dxmin, &
         b% pg% History_Panels7_xmargin, &
         b% pg% History_Panels7_max_width, &
         b% pg% History_Panels7_num_panels, &
         b% pg% History_Panels7_other_ymin, &
         b% pg% History_Panels7_other_ymax, &
         b% pg% History_Panels7_yaxis_reversed, &
         b% pg% History_Panels7_other_yaxis_log, &
         b% pg% History_Panels7_other_dymin, &
         b% pg% History_Panels7_points_name, &
         b% pg% History_Panels7_other_ymargin, &
         b% pg% History_Panels7_other_yaxis_name, &
         b% pg% History_Panels7_ymin, &
         b% pg% History_Panels7_ymax, &
         b% pg% History_Panels7_xaxis_reversed, &
         b% pg% History_Panels7_yaxis_reversed, &
         b% pg% History_Panels7_xaxis_log, &
         b% pg% History_Panels7_yaxis_log, &
         b% pg% History_Panels7_dymin, &
         b% pg% History_Panels7_ymargin, &
         b% pg% History_Panels7_yaxis_name, &
         b% pg% History_Panels7_use_decorator, &
         b% pg% History_Panels7_pgbinary_decorator, &
         ierr)
   end subroutine do_History_Panels7_plot


   subroutine History_Panels8_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels8_plot(b, id, device_id, &
         b% pg% History_Panels8_xleft, b% pg% History_Panels8_xright, &
         b% pg% History_Panels8_ybot, b% pg% History_Panels8_ytop, .false., &
         b% pg% History_Panels8_title, b% pg% History_Panels8_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels8_plot


   subroutine do_History_Panels8_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_history_panels_plot(&
         id, b, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Panels8_xaxis_name, &
         b% pg% History_Panels8_xmin, &
         b% pg% History_Panels8_xmax, &
         b% pg% History_Panels8_dxmin, &
         b% pg% History_Panels8_xmargin, &
         b% pg% History_Panels8_max_width, &
         b% pg% History_Panels8_num_panels, &
         b% pg% History_Panels8_other_ymin, &
         b% pg% History_Panels8_other_ymax, &
         b% pg% History_Panels8_yaxis_reversed, &
         b% pg% History_Panels8_other_yaxis_log, &
         b% pg% History_Panels8_other_dymin, &
         b% pg% History_Panels8_points_name, &
         b% pg% History_Panels8_other_ymargin, &
         b% pg% History_Panels8_other_yaxis_name, &
         b% pg% History_Panels8_ymin, &
         b% pg% History_Panels8_ymax, &
         b% pg% History_Panels8_xaxis_reversed, &
         b% pg% History_Panels8_yaxis_reversed, &
         b% pg% History_Panels8_xaxis_log, &
         b% pg% History_Panels8_yaxis_log, &
         b% pg% History_Panels8_dymin, &
         b% pg% History_Panels8_ymargin, &
         b% pg% History_Panels8_yaxis_name, &
         b% pg% History_Panels8_use_decorator, &
         b% pg% History_Panels8_pgbinary_decorator, &
         ierr)
   end subroutine do_History_Panels8_plot


   subroutine History_Panels9_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels9_plot(b, id, device_id, &
         b% pg% History_Panels9_xleft, b% pg% History_Panels9_xright, &
         b% pg% History_Panels9_ybot, b% pg% History_Panels9_ytop, .false., &
         b% pg% History_Panels9_title, b% pg% History_Panels9_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels9_plot


   subroutine do_History_Panels9_plot(b, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_history_panels_plot(&
         id, b, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         b% pg% History_Panels9_xaxis_name, &
         b% pg% History_Panels9_xmin, &
         b% pg% History_Panels9_xmax, &
         b% pg% History_Panels9_dxmin, &
         b% pg% History_Panels9_xmargin, &
         b% pg% History_Panels9_max_width, &
         b% pg% History_Panels9_num_panels, &
         b% pg% History_Panels9_other_ymin, &
         b% pg% History_Panels9_other_ymax, &
         b% pg% History_Panels9_other_yaxis_reversed, &
         b% pg% History_Panels9_other_yaxis_log, &
         b% pg% History_Panels9_other_dymin, &
         b% pg% History_Panels9_points_name, &
         b% pg% History_Panels9_other_ymargin, &
         b% pg% History_Panels9_other_yaxis_name, &
         b% pg% History_Panels9_ymin, &
         b% pg% History_Panels9_ymax, &
         b% pg% History_Panels9_xaxis_reversed, &
         b% pg% History_Panels9_yaxis_reversed, &
         b% pg% History_Panels9_xaxis_log, &
         b% pg% History_Panels9_yaxis_log, &
         b% pg% History_Panels9_dymin, &
         b% pg% History_Panels9_ymargin, &
         b% pg% History_Panels9_yaxis_name, &
         b% pg% History_Panels9_use_decorator, &
         b% pg% History_Panels9_pgbinary_decorator, &
         ierr)
   end subroutine do_History_Panels9_plot


   subroutine do_history_panels_plot(&
      id, b, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
      hist_xaxis_name, hist_xmin_in, hist_xmax, dxmin, hist_xmargin, &
      hist_max_width, hist_num_panels, &
      hist_other_ymin, hist_other_ymax, &
      hist_other_yaxis_reversed, hist_other_yaxis_log, &
      hist_other_dymin, hist_points_name, &
      hist_other_ymargin, hist_other_yaxis_name, &
      hist_ymin, hist_ymax, &
      hist_xaxis_reversed, hist_yaxis_reversed, &
      hist_xaxis_log, hist_yaxis_log, &
      hist_dymin, hist_ymargin, hist_yaxis_name, &
      use_decorator, pgbinary_decorator, &
      ierr)

      use utils_lib
      use pgstar_support, only : set_xleft_xright, set_ytop_ybot

      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id, hist_num_panels
      logical, intent(in) :: subplot, hist_xaxis_reversed, hist_xaxis_log
      character (len = *), intent(in) :: title, hist_xaxis_name
      real, intent(in) :: &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale, &
         hist_xmin_in, hist_xmax, hist_max_width, hist_xmargin, dxmin
      real, intent(in), dimension(:) :: &
         hist_other_ymin, hist_other_ymax, &
         hist_other_dymin, hist_other_ymargin, &
         hist_ymin, hist_ymax, hist_dymin, hist_ymargin
      logical, intent(in), dimension(:) :: &
         hist_other_yaxis_reversed, hist_other_yaxis_log, &
         hist_yaxis_reversed, hist_yaxis_log
      logical, intent(in) :: use_decorator
      character (len = *), intent(in), dimension(:) :: &
         hist_points_name, hist_other_yaxis_name, hist_yaxis_name
      integer, intent(out) :: ierr
      procedure(pgbinary_decorator_interface), pointer :: pgbinary_decorator

      character (len = strlen) :: yname, other_yname
      real, pointer, dimension(:) :: xvec, yvec, other_yvec
      real, pointer, dimension(:) :: yfile_xdata, other_yfile_xdata
      real, pointer, dimension(:) :: yfile_ydata, other_yfile_ydata
      integer :: i, ii, n, j, k, max_width, step_min, step_max, &
         y_color, other_y_color, yaxis_id, other_yaxis_id, &
         clr_sav, npts, yfile_data_len, other_yfile_data_len
      real :: hist_xmin, xmin, xmax, dx, xleft, xright, &
         ymargin, panel_dy, panel_ytop, panel_ybot, &
         ymin, ymax, dy, ybot, ytop, xpt, ypt, errpt, &
         other_ymin, other_ymax, other_ybot, other_ytop
      logical :: have_yaxis, have_other_yaxis

      integer :: grid_min, grid_max
      integer :: ix, iounit, ishape, num_pts

      include 'formats'

      ierr = 0

      call integer_dict_lookup(b% binary_history_names_dict, hist_xaxis_name, ix, ierr)
      if (ierr /= 0) ix = -1
      if (ix <= 0) then
         write(*, *)
         write(*, *) 'ERROR: failed to find ' // &
            trim(hist_xaxis_name) // ' in history data'
         write(*, *)
         ierr = -1
      end if

      hist_xmin = hist_xmin_in

      if (hist_xaxis_name == 'model_number') then
         max_width = int(hist_max_width)
         step_min = int(hist_xmin)
         if (step_min <= 0) step_min = 1
         step_max = int(hist_xmax)
         if (step_max <= 0) step_max = b% model_number
         if (step_min >= b% model_number) step_min = 1
         if (max_width > 0) step_min = max(step_min, step_max - max_width)
      else
         step_min = 1
         step_max = b% model_number
      end if

      n = count_hist_points(b, step_min, step_max)
      allocate(xvec(n), yvec(n), other_yvec(n), stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'allocate failed for PGBINARY'
         return
      end if

      call get_hist_points(b, step_min, step_max, n, ix, xvec)

      if (hist_xaxis_log) then
         do k = 1, n
            xvec(k) = log10(max(tiny(xvec(k)), abs(xvec(k))))
         end do
      end if

      if (hist_max_width > 0d0 .and. hist_xaxis_name /= 'model_number') then
         ii = n
         do i = ii, 1, -1
            if (xvec(ii) - xvec(i) > hist_max_width) then
               do j = i, ii
                  xvec(j - i + 1) = xvec(j)
               end do
               n = n - i + 1
               hist_xmin = xvec(n) - hist_max_width
               exit
            end if
         end do
      end if

      call set_xleft_xright(&
         n, xvec, hist_xmin, hist_xmax, hist_xmargin, &
         hist_xaxis_reversed, dxmin, xleft, xright)

      call pgsave
      call pgsch(txt_scale)

      ymargin = 0.05
      y_color = clr_Goldenrod
      other_y_color = clr_LightSkyBlue

      panel_dy = (vp_ytop - vp_ybot) / real(hist_num_panels)

      do j = 1, hist_num_panels

         yfile_data_len = 0
         other_yfile_data_len = 0

         yname = hist_yaxis_name(j)
         if (len_trim(yname) == 0) then
            have_yaxis = .false.
         else
            have_yaxis = get1_yvec(yname, yvec)
            if (.not. have_yaxis) then
               if (.not. read_values_from_file(yname, &
                  yfile_xdata, yfile_ydata, yfile_data_len)) then
                  write(*, *) &
                     'bad yaxis for History panels plot ' // trim(yname)
                  cycle
               end if
               have_yaxis = .true.
            end if
         end if

         other_yname = hist_other_yaxis_name(j)
         if (len_trim(hist_other_yaxis_name(j)) == 0) then
            have_other_yaxis = .false.
         else
            have_other_yaxis = get1_yvec(other_yname, other_yvec)
            if (.not. have_other_yaxis) then
               if (.not. read_values_from_file(other_yname, &
                  other_yfile_xdata, other_yfile_ydata, other_yfile_data_len)) then
                  write(*, *) &
                     'bad other yaxis for History panels plot ' // trim(other_yname)
                  cycle
               end if
               have_other_yaxis = .true.
            end if
         end if

         if ((.not. have_yaxis) .and. (.not. have_other_yaxis)) cycle

         panel_ytop = vp_ytop - real(j - 1) * panel_dy
         panel_ybot = panel_ytop - panel_dy

         call pgsvp(vp_xleft, vp_xright, panel_ybot, panel_ytop)

         if (j == 1) then
            if (.not. subplot) then
               call show_model_number_pgbinary(b)
               call show_age_pgbinary(b)
            end if
            call show_title_pgbinary(b, title)
         end if

         if (have_other_yaxis) then
            if (other_yfile_data_len > 0) then
               if (hist_other_yaxis_log(j)) then
                  do k = 1, other_yfile_data_len
                     other_yfile_ydata(k) = &
                        log10(max(tiny(other_yfile_ydata(k)), abs(other_yfile_ydata(k))))
                  end do
               end if
               call set_ytop_ybot(&
                  other_yfile_data_len, other_yfile_ydata, &
                  hist_other_ymin(j), hist_other_ymax(j), -101.0, &
                  hist_other_ymargin(j), hist_other_yaxis_reversed(j), &
                  hist_other_dymin(j), other_ybot, other_ytop)
            else
               if (hist_other_yaxis_log(j)) then
                  do k = 1, n
                     other_yvec(k) = log10(max(tiny(other_yvec(k)), abs(other_yvec(k))))
                  end do
               end if
               call set_ytop_ybot(&
                  n, other_yvec, hist_other_ymin(j), hist_other_ymax(j), -101.0, &
                  hist_other_ymargin(j), hist_other_yaxis_reversed(j), &
                  hist_other_dymin(j), other_ybot, other_ytop)
            end if
            !write(*,1) trim(other_yname), other_ybot, other_ytop
            call pgswin(xleft, xright, other_ybot, other_ytop)
            call pgscf(1)
            call pgsci(1)
            call show_box_pgbinary(b, '', 'CMSTV')
            call pgsci(other_y_color)
            if (hist_other_yaxis_log(j)) then
               call show_right_yaxis_label_pgbinary(b, 'log ' // other_yname)
            else
               call show_right_yaxis_label_pgbinary(b, other_yname)
            end if
            call pgslw(b% pg% pgbinary_lw)
            if (other_yfile_data_len > 0) then
               call pgline(&
                  other_yfile_data_len, other_yfile_xdata, other_yfile_ydata)
               deallocate(other_yfile_xdata, other_yfile_ydata)
               nullify(other_yfile_xdata, other_yfile_ydata)
            else
               call pgline(n, xvec, other_yvec)
            end if
            call pgslw(1)
         end if

         if (have_yaxis) then
            if (yfile_data_len > 0) then
               if (hist_yaxis_log(j)) then
                  do k = 1, yfile_data_len
                     yfile_ydata(k) = log10(max(tiny(yfile_ydata(k)), abs(yfile_ydata(k))))
                  end do
               end if
               call set_ytop_ybot(&
                  yfile_data_len, yfile_ydata, hist_ymin(j), hist_ymax(j), -101.0, &
                  hist_ymargin(j), hist_yaxis_reversed(j), &
                  hist_dymin(j), ybot, ytop)
            else
               if (hist_yaxis_log(j)) then
                  do k = 1, n
                     yvec(k) = log10(max(tiny(yvec(k)), abs(yvec(k))))
                  end do
               end if
               call set_ytop_ybot(&
                  n, yvec, hist_ymin(j), hist_ymax(j), -101.0, &
                  hist_ymargin(j), hist_yaxis_reversed(j), &
                  hist_dymin(j), ybot, ytop)
            end if
            !write(*,1) trim(yname), ybot, ytop
            call pgswin(xleft, xright, ybot, ytop)
            call pgscf(1)
            call pgsci(1)
            if (j < hist_num_panels) then
               if (.not. have_other_yaxis) then
                  call show_box_pgbinary(b, 'BCST1', 'BCMNSTV1')
               else
                  call show_box_pgbinary(b, 'BCST', 'BNSTV')
               end if
            else
               if (.not. have_other_yaxis) then
                  call show_box_pgbinary(b, 'BCNST1', 'BCMNSTV1')
               else
                  call show_box_pgbinary(b, 'BCNST', 'BNSTV')
               end if
            end if

            if (len_trim(hist_points_name(j)) > 0) then
               iounit = 33
               open(unit = iounit, file = trim(hist_points_name(j)), &
                  status = 'old', action = 'read', iostat = ierr)
               if (ierr /= 0) then
                  write(*, '(a)') 'failed to open ' // trim(hist_points_name(j))
                  return
               end if
               read(iounit, *) num_pts
               ishape = b% pg% History_Panel_points_marker ! 5
               call pgsave
               call pgsci(b% pg% History_Panel_points_ci) !1)
               call pgslw(b% pg% History_Panel_points_lw) !2)
               call pgsch(b% pg% History_Panel_points_ch) !1.0)
               do k = 1, num_pts
                  if (b% pg% History_Panel_points_error_bars) then
                     read(iounit, *) xpt, ypt, errpt
                  else
                     read(iounit, *) xpt, ypt
                  end if
                  if (mod(k - 1, b% pg% History_Panel_points_interval) == 0) then
                     if (b% pg% History_Panel_points_error_bars) then
                        call pgmove(xpt, ypt + errpt)
                        call pgdraw(xpt, ypt - errpt)
                     else
                        call pgpt1(xpt, ypt, ishape)
                     end if
                  end if
               end do
               close(iounit)
               call pgunsa
            end if

            call pgsci(y_color)
            if (hist_yaxis_log(j)) then
               call show_left_yaxis_label_pgbinary(b, 'log ' // yname)
            else
               call show_left_yaxis_label_pgbinary(b, yname)
            end if
            call pgslw(b% pg% pgbinary_lw)
            if (yfile_data_len > 0) then
               call pgline(yfile_data_len, yfile_xdata, yfile_ydata)
               deallocate(yfile_xdata, yfile_ydata)
               nullify(yfile_xdata, yfile_ydata)
            else
               call pgline(n, xvec, yvec)
            end if
            call pgslw(1)
         end if

         call pgsci(1)
         call show_pgbinary_decorator(b% binary_id, use_decorator, pgbinary_decorator, j, ierr)
      end do

      if (hist_xaxis_log) then
         call show_xaxis_label_pgbinary(b, 'log ' // hist_xaxis_name)
      else
         call show_xaxis_label_pgbinary(b, hist_xaxis_name)
      end if

      deallocate(xvec, yvec, other_yvec)

      call pgunsa

   contains


      logical function get1_yvec(name, vec)
         character (len = *) :: name
         real, dimension(:), pointer :: vec
         get1_yvec = get1_hist_yvec(b, step_min, step_max, n, name, vec)
      end function get1_yvec


   end subroutine do_history_panels_plot


end module pgbinary_history_panels

