! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
!
! ***********************************************************************

module pgstar_Color_Magnitude

   use star_private_def
   use const_def, only: dp, strlen
   use pgstar_support
   use star_pgstar

   implicit none

contains

   ! ==========================================================================
   ! WINDOW PANEL WRAPPERS 1 THROUGH 9
   ! ==========================================================================

   subroutine Color_Magnitude1_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude1_plot(s, id, device_id, &
         s% pg% Color_Magnitude1_xleft, s% pg% Color_Magnitude1_xright, &
         s% pg% Color_Magnitude1_ybot, s% pg% Color_Magnitude1_ytop, .false., &
         s% pg% Color_Magnitude1_title, s% pg% Color_Magnitude1_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude1_plot

   subroutine do_Color_Magnitude1_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude1_xaxis1_name, s% pg% Color_Magnitude1_xaxis2_name, &
         s% pg% Color_Magnitude1_xmin, s% pg% Color_Magnitude1_xmax, s% pg% Color_Magnitude1_dxmin, &
         s% pg% Color_Magnitude1_xmargin, s% pg% Color_Magnitude1_max_width, s% pg% Color_Magnitude1_num_panels, &
         s% pg% Color_Magnitude1_other_ymin, s% pg% Color_Magnitude1_other_ymax, s% pg% Color_Magnitude1_other_yaxis_reversed, &
         s% pg% Color_Magnitude1_other_yaxis_log, s% pg% Color_Magnitude1_other_dymin, s% pg% Color_Magnitude1_other_ymargin, &
         s% pg% Color_Magnitude1_other_yaxis1_name, s% pg% Color_Magnitude1_other_yaxis2_name, &
         s% pg% Color_Magnitude1_ymin, s% pg% Color_Magnitude1_ymax, s% pg% Color_Magnitude1_xaxis_reversed, &
         s% pg% Color_Magnitude1_yaxis_reversed, s% pg% Color_Magnitude1_xaxis_log, s% pg% Color_Magnitude1_yaxis_log, &
         s% pg% Color_Magnitude1_dymin, s% pg% Color_Magnitude1_ymargin, s% pg% Color_Magnitude1_yaxis1_name, &
         s% pg% Color_Magnitude1_yaxis2_name, s% pg% Color_Magnitude1_use_decorator, s% pg% Color_Magnitude1_pgstar_decorator, &
         ierr)
   end subroutine do_Color_Magnitude1_plot

   subroutine Color_Magnitude2_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude2_plot(s, id, device_id, &
         s% pg% Color_Magnitude2_xleft, s% pg% Color_Magnitude2_xright, &
         s% pg% Color_Magnitude2_ybot, s% pg% Color_Magnitude2_ytop, .false., &
         s% pg% Color_Magnitude2_title, s% pg% Color_Magnitude2_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude2_plot

   subroutine do_Color_Magnitude2_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude2_xaxis1_name, s% pg% Color_Magnitude2_xaxis2_name, &
         s% pg% Color_Magnitude2_xmin, s% pg% Color_Magnitude2_xmax, s% pg% Color_Magnitude2_dxmin, &
         s% pg% Color_Magnitude2_xmargin, s% pg% Color_Magnitude2_max_width, s% pg% Color_Magnitude2_num_panels, &
         s% pg% Color_Magnitude2_other_ymin, s% pg% Color_Magnitude2_other_ymax, s% pg% Color_Magnitude2_other_yaxis_reversed, &
         s% pg% Color_Magnitude2_other_yaxis_log, s% pg% Color_Magnitude2_other_dymin, s% pg% Color_Magnitude2_other_ymargin, &
         s% pg% Color_Magnitude2_other_yaxis1_name, s% pg% Color_Magnitude2_other_yaxis2_name, &
         s% pg% Color_Magnitude2_ymin, s% pg% Color_Magnitude2_ymax, s% pg% Color_Magnitude2_xaxis_reversed, &
         s% pg% Color_Magnitude2_yaxis_reversed, s% pg% Color_Magnitude2_xaxis_log, s% pg% Color_Magnitude2_yaxis_log, &
         s% pg% Color_Magnitude2_dymin, s% pg% Color_Magnitude2_ymargin, s% pg% Color_Magnitude2_yaxis1_name, &
         s% pg% Color_Magnitude2_yaxis2_name, s% pg% Color_Magnitude2_use_decorator, s% pg% Color_Magnitude2_pgstar_decorator, &
         ierr)
   end subroutine do_Color_Magnitude2_plot

   subroutine Color_Magnitude3_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude3_plot(s, id, device_id, &
         s% pg% Color_Magnitude3_xleft, s% pg% Color_Magnitude3_xright, &
         s% pg% Color_Magnitude3_ybot, s% pg% Color_Magnitude3_ytop, .false., &
         s% pg% Color_Magnitude3_title, s% pg% Color_Magnitude3_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude3_plot

   subroutine do_Color_Magnitude3_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude3_xaxis1_name, s% pg% Color_Magnitude3_xaxis2_name, &
         s% pg% Color_Magnitude3_xmin, s% pg% Color_Magnitude3_xmax, s% pg% Color_Magnitude3_dxmin, &
         s% pg% Color_Magnitude3_xmargin, s% pg% Color_Magnitude3_max_width, s% pg% Color_Magnitude3_num_panels, &
         s% pg% Color_Magnitude3_other_ymin, s% pg% Color_Magnitude3_other_ymax, s% pg% Color_Magnitude3_other_yaxis_reversed, &
         s% pg% Color_Magnitude3_other_yaxis_log, s% pg% Color_Magnitude3_other_dymin, s% pg% Color_Magnitude3_other_ymargin, &
         s% pg% Color_Magnitude3_other_yaxis1_name, s% pg% Color_Magnitude3_other_yaxis2_name, &
         s% pg% Color_Magnitude3_ymin, s% pg% Color_Magnitude3_ymax, s% pg% Color_Magnitude3_xaxis_reversed, &
         s% pg% Color_Magnitude3_yaxis_reversed, s% pg% Color_Magnitude3_xaxis_log, s% pg% Color_Magnitude3_yaxis_log, &
         s% pg% Color_Magnitude3_dymin, s% pg% Color_Magnitude3_ymargin, s% pg% Color_Magnitude3_yaxis1_name, &
         s% pg% Color_Magnitude3_yaxis2_name, s% pg% Color_Magnitude3_use_decorator, s% pg% Color_Magnitude3_pgstar_decorator, &
         ierr)
   end subroutine do_Color_Magnitude3_plot

   subroutine Color_Magnitude4_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude4_plot(s, id, device_id, &
         s% pg% Color_Magnitude4_xleft, s% pg% Color_Magnitude4_xright, &
         s% pg% Color_Magnitude4_ybot, s% pg% Color_Magnitude4_ytop, .false., &
         s% pg% Color_Magnitude4_title, s% pg% Color_Magnitude4_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude4_plot

   subroutine do_Color_Magnitude4_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude4_xaxis1_name, s% pg% Color_Magnitude4_xaxis2_name, &
         s% pg% Color_Magnitude4_xmin, s% pg% Color_Magnitude4_xmax, s% pg% Color_Magnitude4_dxmin, &
         s% pg% Color_Magnitude4_xmargin, s% pg% Color_Magnitude4_max_width, s% pg% Color_Magnitude4_num_panels, &
         s% pg% Color_Magnitude4_other_ymin, s% pg% Color_Magnitude4_other_ymax, s% pg% Color_Magnitude4_other_yaxis_reversed, &
         s% pg% Color_Magnitude4_other_yaxis_log, s% pg% Color_Magnitude4_other_dymin, s% pg% Color_Magnitude4_other_ymargin, &
         s% pg% Color_Magnitude4_other_yaxis1_name, s% pg% Color_Magnitude4_other_yaxis2_name, &
         s% pg% Color_Magnitude4_ymin, s% pg% Color_Magnitude4_ymax, s% pg% Color_Magnitude4_xaxis_reversed, &
         s% pg% Color_Magnitude4_yaxis_reversed, s% pg% Color_Magnitude4_xaxis_log, s% pg% Color_Magnitude4_yaxis_log, &
         s% pg% Color_Magnitude4_dymin, s% pg% Color_Magnitude4_ymargin, s% pg% Color_Magnitude4_yaxis1_name, &
         s% pg% Color_Magnitude4_yaxis2_name, s% pg% Color_Magnitude4_use_decorator, s% pg% Color_Magnitude4_pgstar_decorator, &
         ierr)
   end subroutine do_Color_Magnitude4_plot

   subroutine Color_Magnitude5_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude5_plot(s, id, device_id, &
         s% pg% Color_Magnitude5_xleft, s% pg% Color_Magnitude5_xright, &
         s% pg% Color_Magnitude5_ybot, s% pg% Color_Magnitude5_ytop, .false., &
         s% pg% Color_Magnitude5_title, s% pg% Color_Magnitude5_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude5_plot

   subroutine do_Color_Magnitude5_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude5_xaxis1_name, s% pg% Color_Magnitude5_xaxis2_name, &
         s% pg% Color_Magnitude5_xmin, s% pg% Color_Magnitude5_xmax, s% pg% Color_Magnitude5_dxmin, &
         s% pg% Color_Magnitude5_xmargin, s% pg% Color_Magnitude5_max_width, s% pg% Color_Magnitude5_num_panels, &
         s% pg% Color_Magnitude5_other_ymin, s% pg% Color_Magnitude5_other_ymax, s% pg% Color_Magnitude5_other_yaxis_reversed, &
         s% pg% Color_Magnitude5_other_yaxis_log, s% pg% Color_Magnitude5_other_dymin, s% pg% Color_Magnitude5_other_ymargin, &
         s% pg% Color_Magnitude5_other_yaxis1_name, s% pg% Color_Magnitude5_other_yaxis2_name, &
         s% pg% Color_Magnitude5_ymin, s% pg% Color_Magnitude5_ymax, s% pg% Color_Magnitude5_xaxis_reversed, &
         s% pg% Color_Magnitude5_yaxis_reversed, s% pg% Color_Magnitude5_xaxis_log, s% pg% Color_Magnitude5_yaxis_log, &
         s% pg% Color_Magnitude5_dymin, s% pg% Color_Magnitude5_ymargin, s% pg% Color_Magnitude5_yaxis1_name, &
         s% pg% Color_Magnitude5_yaxis2_name, s% pg% Color_Magnitude5_use_decorator, s% pg% Color_Magnitude5_pgstar_decorator, &
         ierr)
   end subroutine do_Color_Magnitude5_plot

   subroutine Color_Magnitude6_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude6_plot(s, id, device_id, &
         s% pg% Color_Magnitude6_xleft, s% pg% Color_Magnitude6_xright, &
         s% pg% Color_Magnitude6_ybot, s% pg% Color_Magnitude6_ytop, .false., &
         s% pg% Color_Magnitude6_title, s% pg% Color_Magnitude6_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude6_plot

   subroutine do_Color_Magnitude6_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude6_xaxis1_name, s% pg% Color_Magnitude6_xaxis2_name, &
         s% pg% Color_Magnitude6_xmin, s% pg% Color_Magnitude6_xmax, s% pg% Color_Magnitude6_dxmin, &
         s% pg% Color_Magnitude6_xmargin, s% pg% Color_Magnitude6_max_width, s% pg% Color_Magnitude6_num_panels, &
         s% pg% Color_Magnitude6_other_ymin, s% pg% Color_Magnitude6_other_ymax, s% pg% Color_Magnitude6_other_yaxis_reversed, &
         s% pg% Color_Magnitude6_other_yaxis_log, s% pg% Color_Magnitude6_other_dymin, s% pg% Color_Magnitude6_other_ymargin, &
         s% pg% Color_Magnitude6_other_yaxis1_name, s% pg% Color_Magnitude6_other_yaxis2_name, &
         s% pg% Color_Magnitude6_ymin, s% pg% Color_Magnitude6_ymax, s% pg% Color_Magnitude6_xaxis_reversed, &
         s% pg% Color_Magnitude6_yaxis_reversed, s% pg% Color_Magnitude6_xaxis_log, s% pg% Color_Magnitude6_yaxis_log, &
         s% pg% Color_Magnitude6_dymin, s% pg% Color_Magnitude6_ymargin, s% pg% Color_Magnitude6_yaxis1_name, &
         s% pg% Color_Magnitude6_yaxis2_name, s% pg% Color_Magnitude6_use_decorator, s% pg% Color_Magnitude6_pgstar_decorator, &
         ierr)
   end subroutine do_Color_Magnitude6_plot

   subroutine Color_Magnitude7_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude7_plot(s, id, device_id, &
         s% pg% Color_Magnitude7_xleft, s% pg% Color_Magnitude7_xright, &
         s% pg% Color_Magnitude7_ybot, s% pg% Color_Magnitude7_ytop, .false., &
         s% pg% Color_Magnitude7_title, s% pg% Color_Magnitude7_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude7_plot

   subroutine do_Color_Magnitude7_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude7_xaxis1_name, s% pg% Color_Magnitude7_xaxis2_name, &
         s% pg% Color_Magnitude7_xmin, s% pg% Color_Magnitude7_xmax, s% pg% Color_Magnitude7_dxmin, &
         s% pg% Color_Magnitude7_xmargin, s% pg% Color_Magnitude7_max_width, s% pg% Color_Magnitude7_num_panels, &
         s% pg% Color_Magnitude7_other_ymin, s% pg% Color_Magnitude7_other_ymax, s% pg% Color_Magnitude7_other_yaxis_reversed, &
         s% pg% Color_Magnitude7_other_yaxis_log, s% pg% Color_Magnitude7_other_dymin, s% pg% Color_Magnitude7_other_ymargin, &
         s% pg% Color_Magnitude7_other_yaxis1_name, s% pg% Color_Magnitude7_other_yaxis2_name, &
         s% pg% Color_Magnitude7_ymin, s% pg% Color_Magnitude7_ymax, s% pg% Color_Magnitude7_xaxis_reversed, &
         s% pg% Color_Magnitude7_yaxis_reversed, s% pg% Color_Magnitude7_xaxis_log, s% pg% Color_Magnitude7_yaxis_log, &
         s% pg% Color_Magnitude7_dymin, s% pg% Color_Magnitude7_ymargin, s% pg% Color_Magnitude7_yaxis1_name, &
         s% pg% Color_Magnitude7_yaxis2_name, s% pg% Color_Magnitude7_use_decorator, s% pg% Color_Magnitude7_pgstar_decorator, &
         ierr)
   end subroutine do_Color_Magnitude7_plot

   subroutine Color_Magnitude8_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude8_plot(s, id, device_id, &
         s% pg% Color_Magnitude8_xleft, s% pg% Color_Magnitude8_xright, &
         s% pg% Color_Magnitude8_ybot, s% pg% Color_Magnitude8_ytop, .false., &
         s% pg% Color_Magnitude8_title, s% pg% Color_Magnitude8_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude8_plot

   subroutine do_Color_Magnitude8_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude8_xaxis1_name, s% pg% Color_Magnitude8_xaxis2_name, &
         s% pg% Color_Magnitude8_xmin, s% pg% Color_Magnitude8_xmax, s% pg% Color_Magnitude8_dxmin, &
         s% pg% Color_Magnitude8_xmargin, s% pg% Color_Magnitude8_max_width, s% pg% Color_Magnitude8_num_panels, &
         s% pg% Color_Magnitude8_other_ymin, s% pg% Color_Magnitude8_other_ymax, s% pg% Color_Magnitude8_other_yaxis_reversed, &
         s% pg% Color_Magnitude8_other_yaxis_log, s% pg% Color_Magnitude8_other_dymin, s% pg% Color_Magnitude8_other_ymargin, &
         s% pg% Color_Magnitude8_other_yaxis1_name, s% pg% Color_Magnitude8_other_yaxis2_name, &
         s% pg% Color_Magnitude8_ymin, s% pg% Color_Magnitude8_ymax, s% pg% Color_Magnitude8_xaxis_reversed, &
         s% pg% Color_Magnitude8_yaxis_reversed, s% pg% Color_Magnitude8_xaxis_log, s% pg% Color_Magnitude8_yaxis_log, &
         s% pg% Color_Magnitude8_dymin, s% pg% Color_Magnitude8_ymargin, s% pg% Color_Magnitude8_yaxis1_name, &
         s% pg% Color_Magnitude8_yaxis2_name, s% pg% Color_Magnitude8_use_decorator, s% pg% Color_Magnitude8_pgstar_decorator, &
         ierr)
   end subroutine do_Color_Magnitude8_plot

   subroutine Color_Magnitude9_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude9_plot(s, id, device_id, &
         s% pg% Color_Magnitude9_xleft, s% pg% Color_Magnitude9_xright, &
         s% pg% Color_Magnitude9_ybot, s% pg% Color_Magnitude9_ytop, .false., &
         s% pg% Color_Magnitude9_title, s% pg% Color_Magnitude9_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude9_plot

   subroutine do_Color_Magnitude9_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude9_xaxis1_name, s% pg% Color_Magnitude9_xaxis2_name, &
         s% pg% Color_Magnitude9_xmin, s% pg% Color_Magnitude9_xmax, s% pg% Color_Magnitude9_dxmin, &
         s% pg% Color_Magnitude9_xmargin, s% pg% Color_Magnitude9_max_width, s% pg% Color_Magnitude9_num_panels, &
         s% pg% Color_Magnitude9_other_ymin, s% pg% Color_Magnitude9_other_ymax, s% pg% Color_Magnitude9_other_yaxis_reversed, &
         s% pg% Color_Magnitude9_other_yaxis_log, s% pg% Color_Magnitude9_other_dymin, s% pg% Color_Magnitude9_other_ymargin, &
         s% pg% Color_Magnitude9_other_yaxis1_name, s% pg% Color_Magnitude9_other_yaxis2_name, &
         s% pg% Color_Magnitude9_ymin, s% pg% Color_Magnitude9_ymax, s% pg% Color_Magnitude9_xaxis_reversed, &
         s% pg% Color_Magnitude9_yaxis_reversed, s% pg% Color_Magnitude9_xaxis_log, s% pg% Color_Magnitude9_yaxis_log, &
         s% pg% Color_Magnitude9_dymin, s% pg% Color_Magnitude9_ymargin, s% pg% Color_Magnitude9_yaxis1_name, &
         s% pg% Color_Magnitude9_yaxis2_name, s% pg% Color_Magnitude9_use_decorator, s% pg% Color_Magnitude9_pgstar_decorator, &
         ierr)
   end subroutine do_Color_Magnitude9_plot

   ! ==========================================================================
   ! CORE LAYOUT PROCESSING UTILITY
   ! ==========================================================================

   subroutine do_Color_Magnitude_plot( &
         id, s, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         color_xaxis1_name, color_xaxis2_name, &
         color_xmin_in, color_xmax, dxmin, color_xmargin, &
         color_max_width, color_num_panels, &
         color_other_ymin, color_other_ymax, &
         color_other_yaxis_reversed, color_other_yaxis_log, &
         color_other_dymin, color_other_ymargin, &
         color_other_yaxis1_name, color_other_yaxis2_name, &
         color_ymin, color_ymax, &
         color_xaxis_reversed, color_yaxis_reversed, &
         color_xaxis_log, color_yaxis_log, &
         color_dymin, color_ymargin, &
         color_yaxis1_name, color_yaxis2_name, &
         color_use_decorator, color_pgstar_decorator, &
         ierr)
      use utils_lib
      use star_def
      use pgstar_colors

      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id, color_num_panels
      logical, intent(in) :: subplot, color_xaxis_reversed, color_xaxis_log
      character(len=*), intent(in) :: title, color_xaxis1_name, color_xaxis2_name
      real, intent(in) :: &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale, &
         color_xmin_in, color_xmax, color_max_width, color_xmargin, dxmin
      real, intent(in), dimension(:) :: &
         color_other_ymin, color_other_ymax, &
         color_other_dymin, color_other_ymargin, &
         color_ymin, color_ymax, color_dymin, color_ymargin
      logical, intent(in), dimension(:) :: &
         color_other_yaxis_reversed, color_other_yaxis_log, &
         color_yaxis_reversed, color_yaxis_log
      logical, intent(in) :: color_use_decorator
      character(len=*), intent(in), dimension(:) :: &
         color_other_yaxis1_name, color_other_yaxis2_name, &
         color_yaxis1_name, color_yaxis2_name
      integer, intent(out) :: ierr
      procedure(pgstar_decorator_interface), pointer :: color_pgstar_decorator

      real, allocatable, dimension(:) :: xvec, yvec, other_yvec, tmp_vec
      integer :: n, j, step_min, step_max, y_color, other_y_color
      real :: color_xmin, xleft, xright, panel_dy, panel_ytop, panel_ybot, ybot, ytop, other_ybot, other_ytop
      logical :: have_yaxis, have_other_yaxis
      character(len=6) :: bottom_xopt

      include 'formats'
      ierr = 0
      step_min = 1
      step_max = s% model_number
      n = count_hist_points(s, step_min, step_max)

      allocate(xvec(n), yvec(n), other_yvec(n), tmp_vec(n), stat=ierr)
      if (ierr /= 0) then
         write(*,*) 'allocate failed for PGSTAR'
         return
      end if

      if (.not. get_axis_vec(color_xaxis1_name, color_xaxis2_name, xvec)) then
         write(*,'(A)')
         write(*,*) 'ERROR: failed to find ' // trim(axis_label(color_xaxis1_name, color_xaxis2_name)) // &
            ' in history data'
         write(*,'(A)')
         call dealloc
         return
      end if

      color_xmin = color_xmin_in
      call set_xleft_xright(n, xvec, color_xmin, color_xmax, color_xmargin, &
                            color_xaxis_reversed, dxmin, xleft, xright)

      ! force decimal x labels only for model numbers; let pgplot format the rest (e.g. star_age)
      if (adjustl(color_xaxis1_name) == 'model_number' .and. len_trim(color_xaxis2_name) == 0) then
         bottom_xopt = 'BCNST1'
      else
         bottom_xopt = 'BCNST'
      end if

      call pgsave
      call pgsch(txt_scale)
      y_color = clr_Goldenrod
      other_y_color = clr_LightSkyBlue
      panel_dy = (vp_ytop - vp_ybot)/real(color_num_panels)

      do j = 1, color_num_panels

         have_yaxis = .false.
         if (len_trim(color_yaxis1_name(j)) > 0) &
            have_yaxis = get_axis_vec(color_yaxis1_name(j), color_yaxis2_name(j), yvec)

         have_other_yaxis = .false.
         if (len_trim(color_other_yaxis1_name(j)) > 0) &
            have_other_yaxis = get_axis_vec(color_other_yaxis1_name(j), color_other_yaxis2_name(j), other_yvec)

         if ((.not. have_yaxis) .and. (.not. have_other_yaxis)) cycle

         ! make sure limits are sensible for plotting
         if (have_yaxis) yvec = max(-100.0, min(100.0, yvec))
         if (have_other_yaxis) other_yvec = max(-100.0, min(100.0, other_yvec))

         panel_ytop = vp_ytop - real(j - 1)*panel_dy
         panel_ybot = panel_ytop - panel_dy
         call pgsvp(vp_xleft, vp_xright, panel_ybot, panel_ytop)

         if (j == 1) then
            if (.not. subplot) then
               call show_model_number_pgstar(s)
               call show_age_pgstar(s)
            end if
            call show_title_pgstar(s, title)
         end if

         if (have_other_yaxis) then
            call set_ytop_ybot(n, other_yvec, color_other_ymin(j), color_other_ymax(j), -101.0, &
                               color_other_ymargin(j), color_other_yaxis_reversed(j), &
                               color_other_dymin(j), other_ybot, other_ytop)
            call pgswin(xleft, xright, other_ybot, other_ytop)
            call pgscf(1)
            call pgsci(clr_Foreground)
            call show_box_pgstar(s, '', 'CMSTV')
            call pgsci(other_y_color)
            call show_right_yaxis_label_pgstar(s, &
               trim(axis_label(color_other_yaxis1_name(j), color_other_yaxis2_name(j))))
            call pgslw(s% pg% pgstar_lw)
            call pgline(n, xvec, other_yvec)
            call pgslw(1)
         end if

         if (have_yaxis) then
            call set_ytop_ybot(n, yvec, color_ymin(j), color_ymax(j), -101.0, &
                               color_ymargin(j), color_yaxis_reversed(j), &
                               color_dymin(j), ybot, ytop)
            call pgswin(xleft, xright, ybot, ytop)
            call pgscf(1)
            call pgsci(clr_Foreground)
            if (j < color_num_panels) then
               if (.not. have_other_yaxis) then
                  call show_box_pgstar(s, 'BCST1', 'BCMNSTV1')
               else
                  call show_box_pgstar(s, 'BCST', 'BNSTV')
               end if
            else
               if (.not. have_other_yaxis) then
                  call show_box_pgstar(s, trim(bottom_xopt), 'BCMNSTV1')
               else
                  call show_box_pgstar(s, trim(bottom_xopt), 'BNSTV')
               end if
            end if
            call pgsci(y_color)
            call show_left_yaxis_label_pgstar(s, &
               trim(axis_label(color_yaxis1_name(j), color_yaxis2_name(j))))
            call pgslw(s% pg% pgstar_lw)
            call pgline(n, xvec, yvec)
            call pgslw(1)
         end if

         call pgsci(clr_Foreground)
         call show_pgstar_decorator(s% id, color_use_decorator, color_pgstar_decorator, j, ierr)
      end do

      call show_xaxis_label_pgstar(s, trim(axis_label(color_xaxis1_name, color_xaxis2_name)))
      call pgunsa
      call dealloc

   contains

      subroutine dealloc
         deallocate(xvec, yvec, other_yvec, tmp_vec)
      end subroutine dealloc

      ! name1 - name2, or name1 if name2 is blank.
      ! each name is a history column or the difference of two columns written as 'B - V'.
      logical function get_axis_vec(name1, name2, vec) result(ok)
         character(len=*), intent(in) :: name1, name2
         real, intent(out) :: vec(:)
         ok = get_expr_vec(name1, vec)
         if (.not. ok .or. len_trim(name2) == 0) return
         ok = get_expr_vec(name2, tmp_vec)
         if (ok) vec = vec - tmp_vec
      end function get_axis_vec

      logical function get_expr_vec(expr, vec) result(ok)
         character(len=*), intent(in) :: expr
         real, intent(out) :: vec(:)
         real, allocatable :: vec2(:)
         integer :: i
         ! a column name may itself contain '-', so try the whole name first
         ok = get_column_vec(expr, vec)
         if (ok) return
         allocate(vec2(size(vec)))
         do i = 2, len_trim(expr) - 1
            if (expr(i:i) /= '-') cycle
            if (.not. get_column_vec(expr(1:i - 1), vec)) cycle
            if (.not. get_column_vec(expr(i + 1:), vec2)) cycle
            vec = vec - vec2
            ok = .true.
            exit
         end do
         deallocate(vec2)
      end function get_expr_vec

      logical function get_column_vec(name, vec) result(ok)
         character(len=*), intent(in) :: name
         real, intent(out) :: vec(:)
         integer :: ix, ierr
         ok = .false.
         if (len_trim(name) == 0) return
         call integer_dict_lookup(s% history_names_dict, trim(adjustl(name)), ix, ierr)
         if (ierr /= 0 .or. ix <= 0) return
         ierr = 0
         call get_hist_points(s, step_min, step_max, n, ix, vec, ierr)
         ok = (ierr == 0)
      end function get_column_vec

      function axis_label(name1, name2) result(label)
         character(len=*), intent(in) :: name1, name2
         character(len=strlen) :: label
         if (len_trim(name2) > 0) then
            label = trim(adjustl(name1)) // ' - ' // trim(adjustl(name2))
         else
            label = adjustl(name1)
         end if
      end function axis_label

   end subroutine do_Color_Magnitude_plot

end module pgstar_Color_Magnitude
