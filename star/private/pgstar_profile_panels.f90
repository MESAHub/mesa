! ***********************************************************************
!
!   Copyright (C) 2013-2019  The MESA Team
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

module pgstar_profile_panels

   use star_private_def
   use const_def
   use pgstar_support
   use pgstar_trho_profile
   use star_pgstar

   implicit none


contains


   subroutine Profile_Panels_plot(id, device_id, array_ix, ierr)
      integer, intent(in) :: id, device_id, array_ix
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_profile_panels_plot(s, id, device_id, array_ix, &
         s% pg% Profile_Panels_xleft(array_ix), s% pg% Profile_Panels_xright(array_ix), &
         s% pg% Profile_Panels_ybot(array_ix), s% pg% Profile_Panels_ytop(array_ix), .false., &
         s% pg% Profile_Panels_title(array_ix), s% pg% Profile_Panels_txt_scale(array_ix), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Profile_Panels_plot


   subroutine do_profile_panels_plot(s, id, device_id, array_ix, &
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
         decorator => s% pg% Profile_Panels1_pgstar_decorator
      case(2)
         decorator => s% pg% Profile_Panels2_pgstar_decorator
      case(3)
         decorator => s% pg% Profile_Panels3_pgstar_decorator
      case(4)
         decorator => s% pg% Profile_Panels4_pgstar_decorator
      case(5)
         decorator => s% pg% Profile_Panels5_pgstar_decorator
      case(6)
         decorator => s% pg% Profile_Panels6_pgstar_decorator
      case(7)
         decorator => s% pg% Profile_Panels7_pgstar_decorator
      case(8)
         decorator => s% pg% Profile_Panels8_pgstar_decorator
      case(9)
         decorator => s% pg% Profile_Panels9_pgstar_decorator
      case default
         write(*, *) "select appropriate profile panels window (1..9)"
         return
      end select

      call Pro_panels_plot(s, device_id, array_ix, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% Profile_Panels_num_panels(array_ix), &

         s% pg% Profile_Panels_xaxis_name(array_ix), &
         s% pg% Profile_Panels_xmin(array_ix), s% pg% Profile_Panels_xmax(array_ix), &
         s% pg% Profile_Panels_xmargin(array_ix), &
         s% pg% Profile_Panels_xaxis_reversed(array_ix), &

         s% pg% Profile_Panels_yaxis_name(array_ix, :), s% pg% Profile_Panels_ycenter(array_ix, :), &
         s% pg% Profile_Panels_ymin(array_ix, :), s% pg% Profile_Panels_ymax(array_ix, :), &
         s% pg% Profile_Panels_dymin(array_ix, :), s% pg% Profile_Panels_ymargin(array_ix, :), &
         s% pg% Profile_Panels_yaxis_reversed(array_ix, :), s% pg% Profile_Panels_yaxis_log(array_ix, :), &

         s% pg% Profile_Panels_other_yaxis_name(array_ix, :), s% pg% Profile_Panels_other_ycenter(array_ix, :), &
         s% pg% Profile_Panels_other_ymin(array_ix, :), s% pg% Profile_Panels_other_ymax(array_ix, :), &
         s% pg% Profile_Panels_other_dymin(array_ix, :), s% pg% Profile_Panels_other_ymargin(array_ix, :), &
         s% pg% Profile_Panels_other_yaxis_reversed(array_ix, :), s% pg% Profile_Panels_other_yaxis_log(array_ix, :), &

         s% pg% Profile_Panels_same_yaxis_range(array_ix, :), &
         s% pg% Profile_Panels_show_mix_regions_on_xaxis(array_ix), &
         s% pg% Profile_Panels_show_grid(array_ix), &

         s% pg% Profile_Panels_use_decorator(array_ix), decorator, &
         ierr)
   end subroutine do_profile_panels_plot


   subroutine Pro_panels_plot(s, device_id, array_ix, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
      panels_num_panels, &
      panels_xaxis_name, &
      panels_xmin_in, panels_xmax_in, &
      panels_xmargin_in, &
      panels_xaxis_reversed, &

      panels_yaxis_name, panels_ycenter, &
      panels_ymin, panels_ymax, &
      panels_dymin, panels_ymargin, &
      panels_yaxis_reversed, panels_yaxis_log, &

      panels_other_yaxis_name, panels_other_ycenter, &
      panels_other_ymin, panels_other_ymax, &
      panels_other_dymin, panels_other_ymargin, &
      panels_other_yaxis_reversed, panels_other_yaxis_log, &

      panels_same_yaxis_range, &
      show_mix_regions, &
      panels_show_grid, &

      use_decorator, pgstar_decorator, &
      ierr)

      use pgstar_abundance, only: do_abundance_panel
      use pgstar_power, only: do_power_panel
      use pgstar_dynamo, only: do_Dynamo_panel
      use pgstar_mixing_Ds, only: do_Mixing_panel
      use pgstar_mode_prop, only: do_mode_propagation_panel
      use pgstar_summary_profile, only: do_summary_profile_panel
      use utils_lib
      use profile_getval, only: get_profile_val, get_profile_id

      type (star_info), pointer :: s
      integer, intent(in) :: &
         device_id, panels_num_panels, array_ix
      real, intent(in) :: &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale, &
         panels_xmin_in, panels_xmax_in, panels_xmargin_in
      logical, intent(in) :: subplot, show_mix_regions
      logical, intent(in), dimension(:) :: &
         panels_yaxis_log, &
         panels_other_yaxis_log, &
         panels_same_yaxis_range
      real, intent(in), dimension(:) :: &
         panels_ymin, panels_other_ymin, &
         panels_ymax, panels_other_ymax, &
         panels_ycenter, panels_other_ycenter, &
         panels_ymargin, panels_other_ymargin, &
         panels_dymin, panels_other_dymin
      character (len = *), intent(in) :: &
         title, panels_xaxis_name
      character (len = *), intent(in), dimension(:) :: &
         panels_yaxis_name, &
         panels_other_yaxis_name
      logical, intent(in) :: &
         panels_xaxis_reversed, use_decorator
      logical, intent(in), dimension(:) :: &
         panels_yaxis_reversed, &
         panels_other_yaxis_reversed
      logical, intent(in) :: panels_show_grid
      integer, intent(out) :: ierr
      procedure(pgstar_decorator_interface), pointer :: pgstar_decorator

      integer :: &
         j, k, k0, k_max, i, nz, kmin, kmax, cnt, y_color, clr_sav, id, &
         other_y_color, grid_min, grid_max, npts, yaxis_id, other_yaxis_id, ishape
      real :: del, xpos, ypos, panel_dy, panel_ybot, panel_ytop, &
         dx, dy, xpos0, dxpos, dypos, dxval, other_ytop, other_ybot, &
         ybot, ytop, xmin, xmax, xleft, xright, xwidth, panels_xmargin, &
         panels_xmin, panels_xmax, xwidth_left_frac, xwidth_right_frac, &
         xwidth_left_of_shock, xwidth_right_of_shock, shock_xmin, shock_xmax, &
         xshock_sp
      real(dp) :: cs, x00, xp1, ms, photosphere_logxm, xshock
      character (len = strlen) :: str, xname, yname, other_yname
      logical :: found_shock
      real, allocatable, dimension(:) :: xvec, yvec, other_yvec, unshifted_xvec
      real, allocatable, dimension(:) :: yfile_xdata, other_yfile_xdata
      real, allocatable, dimension(:) :: yfile_ydata, other_yfile_ydata
      integer :: yfile_data_len, other_yfile_data_len, xaxis_id

      include 'formats'

      ierr = 0

      id = s% id

      call pgsave
      call pgsch(txt_scale)

      y_color = clr_Goldenrod
      other_y_color = clr_LightSkyBlue

      panel_dy = (vp_ytop - vp_ybot) / real(panels_num_panels)

      nz = s% nz
      allocate (xvec(nz), yvec(nz), other_yvec(nz), unshifted_xvec(nz))

      panels_xmin = panels_xmin_in
      panels_xmax = panels_xmax_in
      panels_xmargin = panels_xmargin_in

      xwidth_left_frac = s% pg% Profile_Panels_xwidth_left_div_shock_value
      xwidth_right_frac = s% pg% Profile_Panels_xwidth_right_div_shock_value
      xwidth_left_of_shock = s% pg% Profile_Panels_xwidth_left_of_shock
      xwidth_right_of_shock = s% pg% Profile_Panels_xwidth_right_of_shock

      xshock = 0
      found_shock = .false.

      xaxis_id = get_profile_id(s, panels_xaxis_name)
      if (xaxis_id > 0 .and. s% v_center >= 0 .and. (&
         s% pg% Profile_Panels_show_Mach_1_location .or. &
            xwidth_left_frac > 0 .or. xwidth_right_frac > 0 .or. &
            xwidth_left_of_shock > 0 .or. xwidth_right_of_shock > 0)) then
         found_shock = find_shock(s, xaxis_id, xshock)
         if (found_shock .and. xshock <= 0) then
            write(*, *) 'Panel:', array_ix, ' shock location on xaxis must be positive for tracking location in plot.'
            if (panels_xaxis_name == 'logR') write(*, *) 'perhaps use logR_cm instead of logR?'
            write(*, '(A)')
            found_shock = .false.
         end if
      end if
      xshock_sp = real(xshock, kind = kind(xshock_sp))

      if (found_shock .and. (&
         xwidth_left_frac > 0 .or. xwidth_right_frac > 0 .or. &
            xwidth_left_of_shock > 0 .or. xwidth_right_of_shock > 0)) then

         xmin = get_profile_val(s, xaxis_id, nz)
         xmax = get_profile_val(s, xaxis_id, 1)
         if (xmin > xmax) then ! switch
            dx = xmin; xmin = xmax; xmax = dx
         end if

         shock_xmin = -HUGE(shock_xmin)
         if (xwidth_left_frac > 0) &
            shock_xmin = xshock * (1.0 - xwidth_left_frac)
         if (xwidth_left_of_shock > 0) &
            shock_xmin = max(shock_xmin, xshock_sp - xwidth_left_of_shock)
         if (shock_xmin > xmin) then
            panels_xmin = shock_xmin
            panels_xmargin = 0
         end if

         shock_xmax = HUGE(shock_xmax)
         if (xwidth_right_frac > 0) &
            shock_xmax = xshock * (1.0 + xwidth_right_frac)
         if (xwidth_right_of_shock > 0) &
            shock_xmax = min(shock_xmax, xshock_sp + xwidth_right_of_shock)
         if (shock_xmax < xmax) then
            panels_xmax = shock_xmax
            panels_xmargin = 0
         end if

      end if

      call set_xaxis_bounds(&
         s, panels_xaxis_name, panels_xmin, &
         panels_xmax, panels_xaxis_reversed, panels_xmargin, &
         xvec, xmin, xmax, xleft, xright, dx, &
         grid_min, grid_max, npts, ierr)
      if (ierr /= 0) then
         write(*, *) 'Panel:', array_ix, ' set_xaxis_bounds error in Profile panels -- please check ' // &
            trim(panels_xaxis_name)
         write(*, 1) 'xleft', xleft
         write(*, 1) 'xright', xright
         stop
         return
      end if
      if (xleft == xright) return

      if (found_shock .and. s% pg% Profile_Panels_show_Mach_1_location) then
         ! mark shock location
         call pgsave
         call pgsvp(vp_xleft, vp_xright, vp_ybot, vp_ytop)
         call pgswin(0.0, 1.0, 0.0, 1.0)
         call pgsci(clr_Gray)
         call pgsls(Line_Type_Dash)
         call pgslw(1)
         dx = (xshock - xmin) / (xmax - xmin)
         call pgmove(dx, 0.0)
         call pgdraw(dx, 1.0)
         call pgunsa
      end if

      if (s% pg% Profile_Panels_show_photosphere_location .and. &
         (panels_xaxis_name == 'mass' .or. &
            panels_xaxis_name == 'logxm' .or. &
            panels_xaxis_name == 'radius' .or. &
            panels_xaxis_name == 'radius_cm')) then
         call pgsave
         call pgsvp(vp_xleft, vp_xright, vp_ybot, vp_ytop)
         call pgswin(0.0, 1.0, 0.0, 1.0)
         call pgsci(clr_Gray)
         call pgsls(Line_Type_Dash)
         call pgslw(5)
         if (panels_xaxis_name == 'radius') then
            dx = (s% photosphere_r - xmin) / (xmax - xmin)
         else if (panels_xaxis_name == 'radius_cm') then
            dx = (s% photosphere_r * Rsun - xmin) / (xmax - xmin)
         else if (panels_xaxis_name == 'mass') then
            dx = (s% photosphere_m - xmin) / (xmax - xmin)
         else
            if (s% star_mass > s% photosphere_m) then
               photosphere_logxm = log10(s% star_mass - s% photosphere_m)
            else
               photosphere_logxm = -99
            end if
            dx = 1 - (photosphere_logxm - xmin) / (xmax - xmin)
            write(*, 2) 'Panel:', array_ix, ' photosphere_xm xmin photo_x xmax dx', s% model_number, &
               s% star_mass - s% photosphere_m, &
               xmin, photosphere_logxm, xmax, dx
         end if
         call pgmove(dx, 0.0)
         call pgdraw(dx, 1.0)
         call pgunsa
      end if

      do k = 1, nz
         unshifted_xvec(k) = xvec(k)
      end do
      if (grid_min > 1) then
         do k = 1, npts
            xvec(k) = xvec(k + grid_min - 1)
         end do
      end if

      do j = 1, panels_num_panels

         panel_ytop = vp_ytop - real(j - 1) * panel_dy
         panel_ybot = panel_ytop - panel_dy

         call pgsvp(vp_xleft, vp_xright, panel_ybot, panel_ytop)

         if (j == 1) then
            if (.not. subplot) then
               call show_model_number_pgstar(s)
               call show_age_pgstar(s)
            end if
            call show_title_pgstar(s, title, 0.4)
         end if

         yname = panels_yaxis_name(j)

         if (trim(yname) == 'Abundance') then

            call do_abundance_panel(s, id, device_id, &
               vp_xleft, vp_xright, panel_ybot, panel_ytop, .true., '', txt_scale, &
               panels_xaxis_name, xmin, xmax, panels_xaxis_reversed, &
               panels_ymin(j), panels_ymax(j), &
               .true., (j == panels_num_panels), ierr)
            if (ierr /= 0) then
               write(*, *) 'Panel:', array_ix, ' panels failed in do_abundance_panel'
               stop
            end if
            cycle

         else if (trim(yname) == 'Power') then

            call do_power_panel(s, id, device_id, &
               vp_xleft, vp_xright, panel_ybot, panel_ytop, .true., '', txt_scale, &
               panels_xaxis_name, xmin, xmax, panels_xaxis_reversed, &
               panels_ymin(j), panels_ymax(j), &
               .true., (j == panels_num_panels), ierr)
            if (ierr /= 0) then
               write(*, *) 'Panel:', array_ix, ' panels failed in do_power_panel'
               stop
            end if
            cycle

         else if (trim(yname) == 'Dynamo') then

            call do_Dynamo_panel(s, id, device_id, &
               vp_xleft, vp_xright, panel_ybot, panel_ytop, .true., '', txt_scale, &
               panels_xaxis_name, xmin, xmax, panels_xaxis_reversed, &
               .true., (j == panels_num_panels), ierr)
            if (ierr /= 0) then
               write(*, *) 'Panel:', array_ix, ' panels failed in do_dynamo_panel'
               stop
            end if
            cycle

         else if (trim(yname) == 'Mixing') then

            call do_Mixing_panel(s, id, device_id, &
               vp_xleft, vp_xright, panel_ybot, panel_ytop, .true., '', txt_scale, &
               panels_xaxis_name, xmin, xmax, panels_xaxis_reversed, &
               panels_ymin(j), panels_ymax(j), &
               .true., (j == panels_num_panels), ierr)
            if (ierr /= 0) then
               write(*, *) 'Panel:', array_ix, ' panels failed in do_mixing_panel'
               stop
            end if
            cycle

         else if (trim(yname) == 'Mode_Prop') then

            call do_mode_propagation_panel(s, id, device_id, &
               vp_xleft, vp_xright, panel_ybot, panel_ytop, .true., '', txt_scale, &
               panels_xaxis_name, xmin, xmax, panels_xaxis_reversed, &
               .true., (j == panels_num_panels), ierr)
            if (ierr /= 0) then
               write(*, *) 'Panel:', array_ix, ' panels failed in do_mode_propagation_panel'
               stop
            end if
            cycle

         else if (trim(yname) == 'Summary_Profile') then

            call do_summary_profile_panel(s, id, device_id, &
               vp_xleft, vp_xright, panel_ybot, panel_ytop, .true., '', txt_scale, &
               panels_xaxis_name, xmin, xmax, panels_xaxis_reversed, &
               .true., (j == panels_num_panels), ierr)
            if (ierr /= 0) then
               write(*, *) 'Panel:', array_ix, ' panels failed in do_summary_profile_panel'
               stop
            end if
            cycle

         end if

         other_yaxis_id = 0
         yfile_data_len = 0
         other_yfile_data_len = 0

         if (len_trim(panels_other_yaxis_name(j)) > 0) then
            other_yname = panels_other_yaxis_name(j)
            other_yaxis_id = get_profile_id(s, other_yname)
            if (other_yaxis_id <= 0) then
               if (.not. read_values_from_file(other_yname, &
                  other_yfile_xdata, other_yfile_ydata, other_yfile_data_len)) then
                  write(*, '(A,1X,I1,A,1X,I1,A)') &
                     'Panel:', array_ix, ' bad other yaxis(', j, ') for Profile panels plot name=', trim(other_yname)
                  return
               end if
               if (panels_other_yaxis_log(j)) then
                  do k = 1, other_yfile_data_len
                     other_yfile_ydata(k) = log10(max(tiny(other_yfile_ydata(k)), abs(other_yfile_ydata(k))))
                  end do
               end if
               call set_ytop_ybot(&
                  other_yfile_data_len, other_yfile_ydata, panels_other_ymin(j), &
                  panels_other_ymax(j), panels_other_ycenter(j), panels_other_ymargin(j), &
                  panels_other_yaxis_reversed(j), panels_other_dymin(j), other_ybot, other_ytop)
            else
               do k = 1, npts
                  other_yvec(k) = get_profile_val(s, other_yaxis_id, k + grid_min - 1)
               end do
               if (panels_other_yaxis_log(j)) then
                  do k = 1, npts
                     other_yvec(k) = log10(max(tiny(other_yvec(k)), abs(other_yvec(k))))
                  end do
               end if
               call set_ytop_ybot(&
                  npts, other_yvec, panels_other_ymin(j), panels_other_ymax(j), panels_other_ycenter(j), &
                  panels_other_ymargin(j), panels_other_yaxis_reversed(j), &
                  panels_other_dymin(j), other_ybot, other_ytop)
            end if
         end if

         yaxis_id = get_profile_id(s, yname)
         if (yaxis_id <= 0) then
            if (.not. read_values_from_file(&
               yname, yfile_xdata, yfile_ydata, yfile_data_len)) then
               write(*, '(A,1X,I1,A,1X,I1,A)') &
                  'Panel:', array_ix, ' bad yaxis(', j, ') for Profile panels plot name=', trim(yname)
               return
            end if
            if (panels_yaxis_log(j)) then
               do k = 1, yfile_data_len
                  yfile_ydata(k) = log10(max(tiny(yfile_ydata(k)), abs(yfile_ydata(k))))
               end do
            end if
            call set_ytop_ybot(&
               yfile_data_len, yfile_xdata, panels_ymin(j), panels_ycenter(j), panels_ymax(j), &
               panels_ymargin(j), panels_yaxis_reversed(j), panels_dymin(j), ybot, ytop)
         else
            do k = 1, npts
               yvec(k) = get_profile_val(s, yaxis_id, k + grid_min - 1)
            end do
            if (panels_yaxis_log(j)) then
               do k = 1, npts
                  yvec(k) = log10(max(tiny(yvec(k)), abs(yvec(k))))
               end do
            end if
            call set_ytop_ybot(&
               npts, yvec, panels_ymin(j), panels_ymax(j), panels_ycenter(j), panels_ymargin(j), &
               panels_yaxis_reversed(j), panels_dymin(j), ybot, ytop)
         end if

         if (panels_same_yaxis_range(j) .and. len_trim(panels_other_yaxis_name(j)) > 0) then
            if (other_ybot < ybot) ybot = other_ybot
            if (ybot < other_ybot) other_ybot = ybot
            if (other_ytop > ytop) ytop = other_ytop
            if (ytop > other_ytop) other_ytop = ytop
         end if

         if (len_trim(panels_other_yaxis_name(j)) > 0) then
            call pgswin(xleft, xright, other_ybot, other_ytop)
            call pgscf(1)
            call pgsci(1)
            call show_box_pgstar(s, '', 'CMSTV')
            call pgsci(other_y_color)
            if (panels_other_yaxis_log(j)) then
               call show_right_yaxis_label_pgstar(s, 'log ' // other_yname)
            else
               call show_right_yaxis_label_pgstar(s, other_yname)
            end if
            call pgslw(s% pg% pgstar_lw)
            call pgsci(other_y_color)
            if (other_yfile_data_len > 0) then
               call pgline(&
                  other_yfile_data_len, other_yfile_xdata, other_yfile_ydata)
               deallocate(other_yfile_xdata, other_yfile_ydata)
            else
               call pgline(npts, xvec, other_yvec)
               if (panels_show_grid) then
                  ishape = 21
                  do k = 1, npts
                     call pgpt1(xvec(k), other_yvec(k), ishape)
                  end do
               end if
            end if
            call pgslw(1)
         end if

         call pgswin(xleft, xright, ybot, ytop)
         call pgscf(1)
         call pgsci(1)
         if (j < panels_num_panels) then
            if (other_yaxis_id <= 0 .and. other_yfile_data_len <= 0) then
               call show_box_pgstar(s, 'BCST', 'BCMNSTV')
            else
               call show_box_pgstar(s, 'BCST', 'BNSTV')
            end if
         else
            if (other_yaxis_id <= 0 .and. other_yfile_data_len <= 0) then
               call show_box_pgstar(s, 'BCNST', 'BCMNSTV')
            else
               call show_box_pgstar(s, 'BCNST', 'BNSTV')
            end if
         end if

         call pgqci(clr_sav)
         call pgsci(y_color)
         if (panels_yaxis_log(j)) then
            call show_left_yaxis_label_pgstar(s, 'log ' // yname)
         else
            call show_left_yaxis_label_pgstar(s, yname)
         end if
         call pgslw(s% pg% pgstar_lw)
         if (yfile_data_len > 0) then
            call pgsls(s% pg% pgstar_profile_line_style)
            call pgline(yfile_data_len, yfile_xdata, yfile_ydata)
            call pgsls(1)
            deallocate(yfile_xdata, yfile_ydata)
         else
            call pgsls(s% pg% pgstar_profile_line_style)
            call pgline(npts, xvec, yvec)
            call pgsls(1)
            if (panels_show_grid) then
               ishape = 21
               do k = 1, npts
                  call pgpt1(xvec(k), yvec(k), ishape)
               end do
            end if
         end if
         call pgslw(1)
         call pgsci(1)

         call show_pgstar_decorator(s% id, use_decorator, pgstar_decorator, j, ierr)
      end do

      xname = trim(panels_xaxis_name)
      call show_xaxis_label_pgstar(s, xname)

      if (show_mix_regions) then ! show mix regions at bottom of plot
         call pgslw(10)
         call show_mix_regions_on_xaxis(&
            s, ybot, ytop, grid_min, grid_max, unshifted_xvec)
      end if

      deallocate(xvec, yvec, other_yvec, unshifted_xvec)

      call pgunsa

   end subroutine Pro_panels_plot


end module pgstar_profile_panels

