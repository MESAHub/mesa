! ***********************************************************************
!
!   Copyright (C) 2013  The MESA Team
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

module pgstar_Color_Magnitude

   use star_private_def
   use const_def
   use pgstar_support
   use star_pgstar

   implicit none


contains


   subroutine Color_Magnitude_plot(id, device_id, array_ix, ierr)
      integer, intent(in) :: id, device_id, array_ix
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0

      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Color_Magnitude_plot(s, id, device_id, array_ix, &
         s% pg% Color_Magnitude_xleft(array_ix), s% pg% Color_Magnitude_xright(array_ix), &
         s% pg% Color_Magnitude_ybot(array_ix), s% pg% Color_Magnitude_ytop(array_ix), .false., &
         s% pg% Color_Magnitude_title(array_ix), s% pg% Color_Magnitude_txt_scale(array_ix), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Color_Magnitude_plot


   subroutine do_Color_Magnitude_plot(s, id, device_id, array_ix, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: id, device_id, array_ix
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      procedure(pgstar_decorator_interface), pointer :: decorator

      select case (array_ix)  ! decorator interfaces cannot be arrayified
      case(1)
         decorator => s% pg% Color_Magnitude1_pgstar_decorator
      case(2)
         decorator => s% pg% Color_Magnitude2_pgstar_decorator
      case(3)
         decorator => s% pg% Color_Magnitude3_pgstar_decorator
      case(4)
         decorator => s% pg% Color_Magnitude4_pgstar_decorator
      case(5)
         decorator => s% pg% Color_Magnitude5_pgstar_decorator
      case(6)
         decorator => s% pg% Color_Magnitude6_pgstar_decorator
      case(7)
         decorator => s% pg% Color_Magnitude7_pgstar_decorator
      case(8)
         decorator => s% pg% Color_Magnitude8_pgstar_decorator
      case(9)
         decorator => s% pg% Color_Magnitude9_pgstar_decorator
      case default
         write(*, *) "select appropriate color magnitude window (1..9)"
         return
      end select

      call do_one_cm_plot(&
         id, s, device_id, array_ix, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Color_Magnitude_max_width(array_ix), s% pg% Color_Magnitude_num_panels(array_ix), &

         s% pg% Color_Magnitude_xaxis1_name(array_ix), s% pg% Color_Magnitude_xaxis2_name(array_ix), &
         s% pg% Color_Magnitude_xmin(array_ix), s% pg% Color_Magnitude_xmax(array_ix), &
         s% pg% Color_Magnitude_dxmin(array_ix),  s% pg% Color_Magnitude_xmargin(array_ix), &
         s% pg% Color_Magnitude_xaxis_reversed(array_ix), s% pg% Color_Magnitude_xaxis_log(array_ix), &

         s% pg% Color_Magnitude_yaxis1_name(array_ix, :), s% pg% Color_Magnitude_yaxis2_name(array_ix, :), &
         s% pg% Color_Magnitude_ymin(array_ix, :), s% pg% Color_Magnitude_ymax(array_ix, :), &
         s% pg% Color_Magnitude_dymin(array_ix, :), s% pg% Color_Magnitude_ymargin(array_ix, :), &
         s% pg% Color_Magnitude_yaxis_reversed(array_ix, :), s% pg% Color_Magnitude_yaxis_log(array_ix, :), &

         s% pg% Color_Magnitude_other_yaxis1_name(array_ix, :), s% pg% Color_Magnitude_other_yaxis2_name(array_ix, :), &
         s% pg% Color_Magnitude_other_ymin(array_ix, :), s% pg% Color_Magnitude_other_ymax(array_ix, :), &
         s% pg% Color_Magnitude_other_dymin(array_ix, :), s% pg% Color_Magnitude_other_ymargin(array_ix, :), &
         s% pg% Color_Magnitude_other_yaxis_reversed(array_ix, :), s% pg% Color_Magnitude_other_yaxis_log(array_ix, :), &

         s% pg% Color_Magnitude_use_decorator(array_ix), decorator, &
         ierr)
   end subroutine do_Color_Magnitude_plot


   subroutine do_one_cm_plot(&
      id, s, device_id, array_ix, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
      color_max_width, color_num_panels, &

      color_xaxis1_name, color_xaxis2_name, &
      color_xmin_in, color_xmax, &
      color_dxmin, color_xmargin, &
      color_xaxis_reversed, color_xaxis_log, &

      color_yaxis1_name, color_yaxis2_name, &
      color_ymin, color_ymax, &
      color_dymin, color_ymargin, &
      color_yaxis_reversed, color_yaxis_log, &

      color_other_yaxis1_name, color_other_yaxis2_name, &
      color_other_ymin, color_other_ymax, &
      color_other_dymin, color_other_ymargin, &
      color_other_yaxis_reversed, color_other_yaxis_log, &

      color_use_decorator, color_pgstar_decorator, &
      ierr)
      use utils_lib
      use star_def

      type (star_info), pointer :: s
      integer, intent(in) :: id, device_id, array_ix, color_num_panels
      logical, intent(in) :: subplot, color_xaxis_reversed, color_xaxis_log
      character (len = *), intent(in) :: title, color_xaxis1_name, color_xaxis2_name
      real, intent(in) :: &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale, &
         color_xmin_in, color_xmax, color_max_width, color_xmargin, color_dxmin
      real, intent(in), dimension(:) :: &
         color_other_ymin, color_other_ymax, &
         color_other_dymin, color_other_ymargin, &
         color_ymin, color_ymax, color_dymin, color_ymargin
      logical, intent(in), dimension(:) :: &
         color_other_yaxis_reversed, color_other_yaxis_log, &
         color_yaxis_reversed, color_yaxis_log
      logical, intent(in) :: color_use_decorator
      character (len = *), intent(in), dimension(:) :: &
         color_other_yaxis1_name, color_other_yaxis2_name, &
         color_yaxis1_name, color_yaxis2_name
      integer, intent(out) :: ierr
      procedure(pgstar_decorator_interface), pointer :: color_pgstar_decorator

      character (len = strlen) :: yname, other_yname
      character (len = strlen) :: yname1, yname2, &
         other_yname1, other_yname2
      real, allocatable, dimension(:) :: xvec1, xvec2, yvec1, yvec2, &
         other_yvec1, other_yvec2
      real, allocatable, dimension(:) :: xvec, yvec, other_yvec

      integer :: i, ii, n, j, k, max_width, step_min, step_max, &
         y_color, other_y_color, yaxis_id, other_yaxis_id, &
         clr_sav, npts, yfile_data_len, other_yfile_data_len
      real :: color_xmin, xmin, xmax, dx, xleft, xright, &
         ymargin, panel_dy, panel_ytop, panel_ybot, &
         ymin, ymax, dy, ybot, ytop, &
         other_ymin, other_ymax, other_ybot, other_ytop
      logical :: have_yaxis1, have_other_yaxis1, have_yaxis2, have_other_yaxis2, have_xaxis2
      logical :: have_yaxis, have_other_yaxis

      integer :: grid_min, grid_max
      integer :: ix1, ix2

      include 'formats'
      ierr = 0
      have_xaxis2 = .False.

      step_min = 1
      step_max = s% model_number

      n = count_hist_points(s, step_min, step_max)

      allocate(xvec1(n), xvec2(n), yvec1(n), yvec2(n), other_yvec1(n), other_yvec2(n), stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'allocate failed for PGSTAR'
         return
      end if

      allocate(xvec(n), yvec(n), other_yvec(n), stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'allocate failed for PGSTAR'
         return
      end if

      call integer_dict_lookup(s% history_names_dict, color_xaxis1_name, ix1, ierr)
      if (ierr /= 0) ix1 = -1
      if (ix1 <= 0) then
         write(*, '(A)')
         write(*, *) 'ERROR: failed to find ' // &
            trim(color_xaxis1_name) // ' in history data'
         write(*, '(A)')
         ierr = -1
      end if

      if(len_trim(color_xaxis2_name)>0) THEN
         have_xaxis2 = .True.
         call integer_dict_lookup(s% history_names_dict, color_xaxis2_name, ix2, ierr)
         if (ierr /= 0) ix2 = -1
         if (ix2 <= 0) then
            write(*, '(A)')
            write(*, *) 'ERROR: failed to find ' // &
               trim(color_xaxis2_name) // ' in history data'
            write(*, '(A)')
            ierr = -1
         end if
      end if

      color_xmin = color_xmin_in

      call get_hist_points(s, step_min, step_max, n, ix1, xvec1, ierr)
      if (ierr /= 0) then
         write(*, *) 'pgstar get_hist_points failed ' // trim(color_xaxis1_name)
         deallocate(xvec1, xvec2)
         call dealloc
         ierr = 0
         return
      end if
      if (have_xaxis2) call get_hist_points(s, step_min, step_max, n, ix2, xvec2, ierr)
      if (ierr /= 0) then
         write(*, *) 'pgstar get_hist_points failed ' // trim(color_xaxis2_name)
         deallocate(xvec1, xvec2)
         call dealloc
         ierr = 0
         return
      end if

      xvec = xvec1
      if(have_xaxis2) xvec = xvec - xvec2

      deallocate(xvec1, xvec2)

      call set_xleft_xright(&
         n, xvec, color_xmin, color_xmax, color_xmargin, &
         color_xaxis_reversed, color_dxmin, xleft, xright)

      call pgsave
      call pgsch(txt_scale)

      ymargin = 0.05
      y_color = clr_Goldenrod
      other_y_color = clr_LightSkyBlue

      panel_dy = (vp_ytop - vp_ybot) / real(color_num_panels)

      do j = 1, color_num_panels

         yname1 = color_yaxis1_name(j)
         if (len_trim(yname1) == 0) then
            have_yaxis1 = .false.
         else
            have_yaxis1 = get1_yvec(yname1, yvec1)
         end if

         yname2 = color_yaxis2_name(j)
         if (len_trim(yname2) == 0) then
            have_yaxis2 = .false.
         else
            have_yaxis2 = get1_yvec(yname2, yvec2)
         end if

         other_yname1 = color_other_yaxis1_name(j)
         if (len_trim(color_other_yaxis1_name(j)) == 0) then
            have_other_yaxis1 = .false.
         else
            have_other_yaxis1 = get1_yvec(other_yname1, other_yvec1)
         end if

         other_yname2 = color_other_yaxis2_name(j)
         if (len_trim(color_other_yaxis2_name(j)) == 0) then
            have_other_yaxis2 = .false.
         else
            have_other_yaxis2 = get1_yvec(other_yname2, other_yvec2)
         end if

         if (have_yaxis1) then
            have_yaxis = .True.
         else
            have_yaxis = .false.
         end if

         if (have_other_yaxis1) then
            have_other_yaxis = .True.
         else
            have_other_yaxis = .false.
         end if

         if ((.not. have_yaxis) .and. (.not. have_other_yaxis)) cycle

         yvec = yvec1
         if(have_yaxis2) yvec = yvec - yvec2

         ! Make sure limits are sensible for plotting
         do i = lbound(yvec, dim = 1), ubound(yvec, dim = 1)
            if (yvec(i)>100) yvec(i) = 100
            if (yvec(i)<-100) yvec(i) = -100
         end do

         other_yvec = other_yvec1
         if (have_other_yaxis2) other_yvec = other_yvec - other_yvec2

         ! Make sure limits are sensible for plotting
         do i = lbound(other_yvec, dim = 1), ubound(other_yvec, dim = 1)
            if (other_yvec(i)>100) other_yvec(i) = 100
            if (other_yvec(i)<-100) other_yvec(i) = -100
         end do

         panel_ytop = vp_ytop - real(j - 1) * panel_dy
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
            call set_ytop_ybot(&
               n, other_yvec, color_other_ymin(j), color_other_ymax(j), -101.0, &
               color_other_ymargin(j), color_other_yaxis_reversed(j), &
               color_other_dymin(j), other_ybot, other_ytop)
            call pgswin(xleft, xright, other_ybot, other_ytop)
            call pgscf(1)
            call pgsci(1)
            call show_box_pgstar(s, '', 'CMSTV')
            call pgsci(other_y_color)

            if (have_other_yaxis2) then
               call show_right_yaxis_label_pgstar(s, trim(create_label(other_yname1, other_yname2)))
            else
               call show_right_yaxis_label_pgstar(s, trim(create_label(other_yname1, '')))
            end if

            call pgslw(s% pg% pgstar_lw)
            call pgline(n, xvec, other_yvec)
            call pgslw(1)
         end if

         if (have_yaxis) then
            call set_ytop_ybot(&
               n, yvec, color_ymin(j), color_ymax(j), -101.0, &
               color_ymargin(j), color_yaxis_reversed(j), &
               color_dymin(j), ybot, ytop)
            call pgswin(xleft, xright, ybot, ytop)
            call pgscf(1)
            call pgsci(1)
            if (j < color_num_panels) then
               if (.not. have_other_yaxis) then
                  call show_box_pgstar(s, 'BCST1', 'BCMNSTV1')
               else
                  call show_box_pgstar(s, 'BCST', 'BNSTV')
               end if
            else
               if (.not. have_other_yaxis) then
                  call show_box_pgstar(s, 'BCNST1', 'BCMNSTV1')
               else
                  call show_box_pgstar(s, 'BCNST', 'BNSTV')
               end if
            end if
            call pgsci(y_color)
            if (have_yaxis2) then
               call show_left_yaxis_label_pgstar(s, trim(create_label(yname1, yname2)))
            else
               call show_left_yaxis_label_pgstar(s, trim(create_label(yname1, '')))
            end if
            call pgslw(s% pg% pgstar_lw)
            call pgline(n, xvec, yvec)
            call pgslw(1)
         end if

         call pgsci(1)

         call show_pgstar_decorator(s% id, color_use_decorator, color_pgstar_decorator, j, ierr)

      end do

      if (have_xaxis2) then
         call show_xaxis_label_pgstar(s, trim(create_label(color_xaxis1_name, color_xaxis2_name)))
      else
         call show_xaxis_label_pgstar(s, trim(create_label(color_xaxis1_name, '')))
      end if

      call pgunsa

      call dealloc

   contains

      subroutine dealloc
         deallocate(xvec, yvec, other_yvec, yvec1, yvec2, other_yvec1, other_yvec2)
      end subroutine dealloc

      logical function get1_yvec(name, vec)
         character (len = *) :: name
         real, dimension(:), allocatable :: vec
         get1_yvec = get1_hist_yvec(s, step_min, step_max, n, name, vec)
      end function get1_yvec

      function create_label(str1, str2) result(new_str)
         character(len = *) :: str1, str2
         integer :: len1, len2, endStr1
         character(len = strlen) :: new_str

         new_str = ''

         len1 = len_trim(str1)
         len2 = len_trim(str2)

         !Strings start with bc_
         if(str1(1:3)=='bc_') then
            new_str = str1(4:len1)
            !String start with abs_mag_
         else if(str1(1:8)=='abs_mag_') then
            new_str = 'M\d' // str1(9:len1) // '\u'
         else
            new_str = str1(1:len1)
         end if
         endStr1 = len_trim(new_str)

         if(len2>0)then
            !Strings start with bc_
            if(str1(1:1)=='b') then
               new_str(endStr1 + 1:) = ' - ' // str2(4:len2)
               !String start with abs_mag_
            else if(str1(1:1)=='a') then
               new_str(endStr1 + 1:) = ' - M\d' // str2(9:len2) // '\u'
            else
               new_str(endStr1 + 1:) = str2(1:len2)
            end if
         end if

      end function create_label

   end subroutine do_one_cm_plot

end module pgstar_Color_Magnitude

