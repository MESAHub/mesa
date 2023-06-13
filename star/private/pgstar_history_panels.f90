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

module pgstar_history_panels

   use star_private_def
   use const_def
   use pgstar_support
   use star_pgstar

   implicit none


contains


   subroutine History_Panels_plot(id, device_id, array_ix, ierr)
      integer, intent(in) :: id, device_id, array_ix
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Panels_plot(s, id, device_id, array_ix, &
         s% pg% History_Panels_xleft(array_ix), s% pg% History_Panels_xright(array_ix), &
         s% pg% History_Panels_ybot(array_ix), s% pg% History_Panels_ytop(array_ix), .false., &
         s% pg% History_Panels_title(array_ix), s% pg% History_Panels_txt_scale(array_ix), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Panels_plot


   subroutine do_History_Panels_plot(s, id, device_id, array_ix, &
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
         decorator => s% pg% History_Panels1_pgstar_decorator
      case(2)
         decorator => s% pg% History_Panels2_pgstar_decorator
      case(3)
         decorator => s% pg% History_Panels3_pgstar_decorator
      case(4)
         decorator => s% pg% History_Panels4_pgstar_decorator
      case(5)
         decorator => s% pg% History_Panels5_pgstar_decorator
      case(6)
         decorator => s% pg% History_Panels6_pgstar_decorator
      case(7)
         decorator => s% pg% History_Panels7_pgstar_decorator
      case(8)
         decorator => s% pg% History_Panels8_pgstar_decorator
      case(9)
         decorator => s% pg% History_Panels9_pgstar_decorator
      case default
         write(*, *) "select appropriate history panels window (1..9)"
         return
      end select

      call do_one_panels_plot(&
         id, s, device_id, array_ix, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% History_Panels_max_width(array_ix), s% pg% History_Panels_num_panels(array_ix), &

         s% pg% History_Panels_xaxis_name(array_ix), &
         s% pg% History_Panels_xmin(array_ix), s% pg% History_Panels_xmax(array_ix), &
         s% pg% History_Panels_dxmin(array_ix), s% pg% History_Panels_xmargin(array_ix), &
         s% pg% History_Panels_xaxis_reversed(array_ix), s% pg% History_Panels_xaxis_log(array_ix), &

         s% pg% History_Panels_yaxis_name(array_ix, :), &
         s% pg% History_Panels_ymin(array_ix, :), s% pg% History_Panels_ymax(array_ix, :), &
         s% pg% History_Panels_dymin(array_ix, :), s% pg% History_Panels_ymargin(array_ix, :), &
         s% pg% History_Panels_yaxis_reversed(array_ix, :), s% pg% History_Panels_yaxis_log(array_ix, :), &

         s% pg% History_Panels_other_yaxis_name(array_ix, :), &
         s% pg% History_Panels_other_ymin(array_ix, :), s% pg% History_Panels_other_ymax(array_ix, :), &
         s% pg% History_Panels_other_dymin(array_ix, :), s% pg% History_Panels_other_ymargin(array_ix, :), &
         s% pg% History_Panels_other_yaxis_reversed(array_ix, :), s% pg% History_Panels_other_yaxis_log(array_ix, :), &

         s% pg% History_Panels_same_yaxis_range(array_ix, :), &
         s% pg% History_Panels_points_name(array_ix, :), &
         s% pg% History_Panels_automatic_star_age_units(array_ix), &

         s% pg% History_Panels_use_decorator(array_ix), decorator, &
         ierr)
   end subroutine do_History_Panels_plot


   subroutine do_one_panels_plot( &
      id, s, device_id, array_ix, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
      hist_max_width, hist_num_panels, &
      hist_xaxis_name_in, &
      hist_xmin_in, hist_xmax,&
      dxmin, hist_xmargin, &
      hist_xaxis_reversed, hist_xaxis_log, &
      hist_yaxis_name, &
      hist_ymin, hist_ymax, &
      hist_dymin, hist_ymargin,  &
      hist_yaxis_reversed, hist_yaxis_log, &
      hist_other_yaxis_name, &
      hist_other_ymin, hist_other_ymax, &
      hist_other_dymin, hist_other_ymargin, &
      hist_other_yaxis_reversed, hist_other_yaxis_log, &

      hist_same_yaxis_range, &
      hist_points_name, &
      automatic_star_age_units, &

      use_decorator, pgstar_decorator, &
      ierr)
      use utils_lib
      use chem_def
      use net_def
      use net_lib, only: get_net_reaction_table
      use rates_def, only: rates_reaction_id_max
      use const_def, only: Msun, Rsun

      type (star_info), pointer :: s
      integer, intent(in) :: id, device_id, array_ix, hist_num_panels
      logical, intent(in) :: subplot, automatic_star_age_units, hist_xaxis_reversed, hist_xaxis_log
      character (len = *), intent(in) :: title, hist_xaxis_name_in
      real, intent(in) :: &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale, &
         hist_xmin_in, hist_xmax, hist_max_width, hist_xmargin, dxmin
      real, intent(in), dimension(:) :: &
         hist_other_ymin, hist_other_ymax, &
         hist_other_dymin, hist_other_ymargin, &
         hist_ymin, hist_ymax, hist_dymin, hist_ymargin
      logical, intent(in), dimension(:) :: &
         hist_other_yaxis_reversed, hist_other_yaxis_log, &
         hist_yaxis_reversed, hist_yaxis_log, hist_same_yaxis_range
      logical, intent(in) :: use_decorator
      character (len = *), intent(in), dimension(:) :: &
         hist_points_name, hist_other_yaxis_name, hist_yaxis_name
      integer, intent(out) :: ierr
      procedure(pgstar_decorator_interface), pointer :: pgstar_decorator

      character (len = strlen) :: yname, other_yname, hist_xaxis_name
      real, allocatable, dimension(:) :: xvec, yvec, other_yvec
      real, allocatable, dimension(:) :: yfile_xdata, other_yfile_xdata
      real, allocatable, dimension(:) :: yfile_ydata, other_yfile_ydata
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

      hist_xaxis_name = hist_xaxis_name_in
      hist_xmin = hist_xmin_in

      step_min = 1
      step_max = s% model_number

      if (hist_xaxis_name == 'model_number') then
         max_width = int(hist_max_width)
         step_min = int(hist_xmin)
         if (step_min <= 0) step_min = 1
         step_max = int(hist_xmax)
         if (step_max <= 0) step_max = s% model_number
         if (step_min >= s% model_number) step_min = 1
         if (max_width > 0) step_min = max(step_min, step_max - max_width)
      else if (hist_xaxis_name == 'star_age' .and. automatic_star_age_units) then
         if (s% star_age > 1d0) then
            hist_xaxis_name = 'star_age_yr'
         else if (s% star_age * secyer > secday) then
            hist_xaxis_name = 'star_age_day'
         else if (s% star_age * secyer > 60 * 60) then
            hist_xaxis_name = 'star_age_hr'
         else if (s% star_age * secyer > 60) then
            hist_xaxis_name = 'star_age_min'
         else
            hist_xaxis_name = 'star_age_sec'
         end if
      end if

      call integer_dict_lookup(s% history_names_dict, hist_xaxis_name, ix, ierr)
      if (ierr /= 0) ix = -1
      if (ix <= 0) then
         write(*, '(A)')
         write(*, *) 'ERROR: failed to find ' // &
            trim(hist_xaxis_name) // ' in history data'
         write(*, '(A)')
         ierr = -1
      end if

      n = count_hist_points(s, step_min, step_max)
      allocate(xvec(n), yvec(n), other_yvec(n), stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'allocate failed for PGSTAR'
         return
      end if

      call get_hist_points(s, step_min, step_max, n, ix, xvec, ierr)
      if (ierr /= 0) then
         write(*, *) 'pgstar do_history_panels_plot get_hist_points failed ' // trim(hist_xaxis_name)
         ierr = 0
         !return
      end if

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
               call show_model_number_pgstar(s)
               call show_age_pgstar(s)
            end if
            call show_title_pgstar(s, title)
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
         end if

         if (have_yaxis .and. have_other_yaxis .and. hist_same_yaxis_range(j)) then
            if (other_ybot < ybot) ybot = other_ybot
            if (ybot < other_ybot) other_ybot = ybot
            if (other_ytop > ytop) ytop = other_ytop
            if (ytop > other_ytop) other_ytop = ytop
         end if

         if (have_other_yaxis) then
            !write(*,1) trim(other_yname), other_ybot, other_ytop
            call pgswin(xleft, xright, other_ybot, other_ytop)
            call pgscf(1)
            call pgsci(1)
            call show_box_pgstar(s, '', 'CMSTV')
            call pgsci(other_y_color)
            if (hist_other_yaxis_log(j)) then
               call show_right_yaxis_label_pgstar(s, 'log ' // other_yname)
            else
               call show_right_yaxis_label_pgstar(s, other_yname)
            end if
            call pgslw(s% pg% pgstar_lw)
            if (other_yfile_data_len > 0) then
               call pgline(&
                  other_yfile_data_len, other_yfile_xdata, other_yfile_ydata)
               deallocate(other_yfile_xdata, other_yfile_ydata)
            else
               call pgline(n, xvec, other_yvec)
            end if
            call pgslw(1)
         end if

         if (have_yaxis) then
            !write(*,1) trim(yname), ybot, ytop
            call pgswin(xleft, xright, ybot, ytop)
            call pgscf(1)
            call pgsci(1)
            if (j < hist_num_panels) then
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

            if (len_trim(hist_points_name(j)) > 0) then
               iounit = 33
               open(unit = iounit, file = trim(hist_points_name(j)), &
                  status = 'old', action = 'read', iostat = ierr)
               if (ierr /= 0) then
                  write(*, '(a)') 'failed to open ' // trim(hist_points_name(j))
                  return
               end if
               read(iounit, *) num_pts
               ishape = s% pg% History_Panel_points_marker ! 5
               call pgsave
               call pgsci(s% pg% History_Panel_points_ci) !1)
               call pgslw(s% pg% History_Panel_points_lw) !2)
               call pgsch(s% pg% History_Panel_points_ch) !1.0)
               do k = 1, num_pts
                  if (s% pg% History_Panel_points_error_bars) then
                     read(iounit, *) xpt, ypt, errpt
                  else
                     read(iounit, *) xpt, ypt
                  end if
                  if (mod(k - 1, s% pg% History_Panel_points_interval) == 0) then
                     if (s% pg% History_Panel_points_error_bars) then
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
               call show_left_yaxis_label_pgstar(s, 'log ' // yname)
            else
               call show_left_yaxis_label_pgstar(s, yname)
            end if
            call pgslw(s% pg% pgstar_lw)
            if (yfile_data_len > 0) then
               call pgsls(s% pg% pgstar_history_line_style)
               call pgline(yfile_data_len, yfile_xdata, yfile_ydata)
               call pgsls(1)
               deallocate(yfile_xdata, yfile_ydata)
            else
               call pgsls(s% pg% pgstar_history_line_style)
               call pgline(n, xvec, yvec)
               call pgsls(1)
            end if
            call pgslw(1)
         end if

         call pgsci(1)
         call show_pgstar_decorator(s%id, use_decorator, pgstar_decorator, j, ierr)
      end do

      if (hist_xaxis_log) then
         call show_xaxis_label_pgstar(s, 'log ' // hist_xaxis_name)
      else
         call show_xaxis_label_pgstar(s, hist_xaxis_name)
      end if

      deallocate(xvec, yvec, other_yvec)

      call pgunsa

   contains


      logical function get1_yvec(name, vec)
         character (len = *) :: name
         real, dimension(:), allocatable :: vec
         get1_yvec = get1_hist_yvec(s, step_min, step_max, n, name, vec)
      end function get1_yvec


   end subroutine do_one_panels_plot


end module pgstar_history_panels

