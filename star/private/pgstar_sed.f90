! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

! plots the SED that the colors module computed for the current model,
! with the filter transmission curves drawn over it

module pgstar_sed

   use star_private_def
   use const_def, only: dp
   use pgstar_support
   use star_pgstar

   implicit none

contains

   subroutine SED1_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_SED1_plot(s, id, device_id, &
         s% pg% SED1_xleft, s% pg% SED1_xright, &
         s% pg% SED1_ybot, s% pg% SED1_ytop, .false., &
         s% pg% SED1_title, s% pg% SED1_txt_scale, ierr)
      call pgebuf()
   end subroutine SED1_plot

   subroutine do_SED1_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% SED1_lambda_min, s% pg% SED1_lambda_max, &
         s% pg% SED1_flux_min, s% pg% SED1_flux_max, &
         s% pg% SED1_filter_scale, ierr)
   end subroutine do_SED1_plot

   subroutine SED2_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_SED2_plot(s, id, device_id, &
         s% pg% SED2_xleft, s% pg% SED2_xright, &
         s% pg% SED2_ybot, s% pg% SED2_ytop, .false., &
         s% pg% SED2_title, s% pg% SED2_txt_scale, ierr)
      call pgebuf()
   end subroutine SED2_plot

   subroutine do_SED2_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% SED2_lambda_min, s% pg% SED2_lambda_max, &
         s% pg% SED2_flux_min, s% pg% SED2_flux_max, &
         s% pg% SED2_filter_scale, ierr)
   end subroutine do_SED2_plot

   subroutine SED3_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_SED3_plot(s, id, device_id, &
         s% pg% SED3_xleft, s% pg% SED3_xright, &
         s% pg% SED3_ybot, s% pg% SED3_ytop, .false., &
         s% pg% SED3_title, s% pg% SED3_txt_scale, ierr)
      call pgebuf()
   end subroutine SED3_plot

   subroutine do_SED3_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% SED3_lambda_min, s% pg% SED3_lambda_max, &
         s% pg% SED3_flux_min, s% pg% SED3_flux_max, &
         s% pg% SED3_filter_scale, ierr)
   end subroutine do_SED3_plot

   subroutine SED4_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_SED4_plot(s, id, device_id, &
         s% pg% SED4_xleft, s% pg% SED4_xright, &
         s% pg% SED4_ybot, s% pg% SED4_ytop, .false., &
         s% pg% SED4_title, s% pg% SED4_txt_scale, ierr)
      call pgebuf()
   end subroutine SED4_plot

   subroutine do_SED4_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% SED4_lambda_min, s% pg% SED4_lambda_max, &
         s% pg% SED4_flux_min, s% pg% SED4_flux_max, &
         s% pg% SED4_filter_scale, ierr)
   end subroutine do_SED4_plot

   subroutine SED5_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_SED5_plot(s, id, device_id, &
         s% pg% SED5_xleft, s% pg% SED5_xright, &
         s% pg% SED5_ybot, s% pg% SED5_ytop, .false., &
         s% pg% SED5_title, s% pg% SED5_txt_scale, ierr)
      call pgebuf()
   end subroutine SED5_plot

   subroutine do_SED5_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% SED5_lambda_min, s% pg% SED5_lambda_max, &
         s% pg% SED5_flux_min, s% pg% SED5_flux_max, &
         s% pg% SED5_filter_scale, ierr)
   end subroutine do_SED5_plot

   subroutine SED6_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_SED6_plot(s, id, device_id, &
         s% pg% SED6_xleft, s% pg% SED6_xright, &
         s% pg% SED6_ybot, s% pg% SED6_ytop, .false., &
         s% pg% SED6_title, s% pg% SED6_txt_scale, ierr)
      call pgebuf()
   end subroutine SED6_plot

   subroutine do_SED6_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% SED6_lambda_min, s% pg% SED6_lambda_max, &
         s% pg% SED6_flux_min, s% pg% SED6_flux_max, &
         s% pg% SED6_filter_scale, ierr)
   end subroutine do_SED6_plot

   subroutine SED7_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_SED7_plot(s, id, device_id, &
         s% pg% SED7_xleft, s% pg% SED7_xright, &
         s% pg% SED7_ybot, s% pg% SED7_ytop, .false., &
         s% pg% SED7_title, s% pg% SED7_txt_scale, ierr)
      call pgebuf()
   end subroutine SED7_plot

   subroutine do_SED7_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% SED7_lambda_min, s% pg% SED7_lambda_max, &
         s% pg% SED7_flux_min, s% pg% SED7_flux_max, &
         s% pg% SED7_filter_scale, ierr)
   end subroutine do_SED7_plot

   subroutine SED8_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_SED8_plot(s, id, device_id, &
         s% pg% SED8_xleft, s% pg% SED8_xright, &
         s% pg% SED8_ybot, s% pg% SED8_ytop, .false., &
         s% pg% SED8_title, s% pg% SED8_txt_scale, ierr)
      call pgebuf()
   end subroutine SED8_plot

   subroutine do_SED8_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% SED8_lambda_min, s% pg% SED8_lambda_max, &
         s% pg% SED8_flux_min, s% pg% SED8_flux_max, &
         s% pg% SED8_filter_scale, ierr)
   end subroutine do_SED8_plot

   subroutine SED9_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_SED9_plot(s, id, device_id, &
         s% pg% SED9_xleft, s% pg% SED9_xright, &
         s% pg% SED9_ybot, s% pg% SED9_ytop, .false., &
         s% pg% SED9_title, s% pg% SED9_txt_scale, ierr)
      call pgebuf()
   end subroutine SED9_plot

   subroutine do_SED9_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      integer, intent(out) :: ierr
      call do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% SED9_lambda_min, s% pg% SED9_lambda_max, &
         s% pg% SED9_flux_min, s% pg% SED9_flux_max, &
         s% pg% SED9_filter_scale, ierr)
   end subroutine do_SED9_plot

   subroutine do_SED_plot(s, id, device_id, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         lambda_min, lambda_max, flux_min, flux_max, filter_scale, ierr)
      use colors_def, only: Colors_General_Info
      use colors_lib, only: colors_ptr
      use pgstar_colors

      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      real, intent(in) :: lambda_min, lambda_max, flux_min, flux_max, filter_scale
      integer, intent(out) :: ierr

      ! colors for the filter curves, cycled through in filter order
      integer, parameter :: num_filter_clrs = 8
      integer, parameter :: filter_clrs(num_filter_clrs) = [ &
         clr_DodgerBlue, clr_LimeGreen, clr_Crimson, clr_Goldenrod, &
         clr_MediumSlateBlue, clr_Teal, clr_Coral, clr_Plum]

      type(Colors_General_Info), pointer :: cs
      real, allocatable :: xvec(:), yvec(:), filter_x(:), filter_y(:)
      logical, allocatable :: in_range(:)
      real :: xleft, xright, ybot, ytop, dy, peak
      integer :: i, k, n, i_peak
      character(len=100) :: filter_name

      ierr = 0

      call pgsave
      call pgsch(txt_scale)
      call pgsvp(vp_xleft, vp_xright, vp_ybot, vp_ytop)
      if (.not. subplot) then
         call show_model_number_pgstar(s)
         call show_age_pgstar(s)
      end if
      call show_title_pgstar(s, title)

      ! a missing SED is not an error, e.g. before the first step or with use_colors = .false.
      nullify(cs)
      if (s% colors_handle > 0) then
         call colors_ptr(s% colors_handle, cs, ierr)
         if (ierr /= 0) nullify(cs)
         ierr = 0
      end if
      if (associated(cs)) then
         if (.not. (cs% pgstar_sed_valid .and. allocated(cs% current_sed_wavelengths) &
               .and. allocated(cs% current_sed_fluxes))) nullify(cs)
      end if
      if (.not. associated(cs)) then
         call show_message('no SED available (needs use_colors = .true.)')
         return
      end if

      n = size(cs% current_sed_wavelengths)
      allocate(xvec(n), yvec(n), in_range(n))
      xvec = real(cs% current_sed_wavelengths)

      if (lambda_min == -101) then
         xleft = minval(xvec)
      else
         xleft = lambda_min
      end if
      if (lambda_max == -101) then
         xright = maxval(xvec)
      else
         xright = lambda_max
      end if

      in_range = (cs% current_sed_fluxes > 0) .and. (xvec >= xleft) .and. (xvec <= xright)
      if (xright <= xleft .or. .not. any(in_range)) then
         call show_message('no SED points in the wavelength range')
         deallocate(xvec, yvec, in_range)
         return
      end if

      where (cs% current_sed_fluxes > 0)
         yvec = real(log10(cs% current_sed_fluxes))
      elsewhere
         yvec = -huge(yvec)
      end where

      if (flux_max == -101) then
         ytop = maxval(yvec, mask=in_range)
      else
         ytop = flux_max
      end if
      if (flux_min == -101) then
         ybot = minval(yvec, mask=in_range)
      else
         ybot = flux_min
      end if
      dy = max(0.05*(ytop - ybot), 0.1)
      if (flux_max == -101) ytop = ytop + dy
      if (flux_min == -101) ybot = ybot - dy
      if (ytop <= ybot) ytop = ybot + 1

      ! points with no flux go off the bottom of the plot
      where (yvec < ybot) yvec = ybot - (ytop - ybot)

      call pgswin(xleft, xright, ybot, ytop)
      call pgscf(1)
      call pgsci(clr_Foreground)
      call show_box_pgstar(s, 'BCNST', 'BCNSTV')
      call show_xaxis_label_pgstar(s, 'wavelength (\A)')
      call show_left_yaxis_label_pgstar(s, 'log F\d\gl\u (erg s\u-1\d cm\u-2\d \A\u-1\d)')

      call pgslw(s% pg% pgstar_lw)
      call pgline(n, xvec, yvec)
      call pgslw(1)

      ! filter curves hang from the top of the plot, scaled to their peak
      if (cs% filters_loaded .and. filter_scale > 0) then
         call pgsls(Line_Type_Dash)
         do i = 1, size(cs% filters)
            if (.not. (allocated(cs% filters(i)% wavelengths) .and. &
                  allocated(cs% filters(i)% transmission))) cycle
            k = size(cs% filters(i)% wavelengths)
            if (k < 2) cycle
            peak = real(maxval(cs% filters(i)% transmission))
            if (peak <= 0) cycle

            allocate(filter_x(k), filter_y(k))
            filter_x = real(cs% filters(i)% wavelengths)
            filter_y = ytop - (1 - real(cs% filters(i)% transmission)/peak)*(ytop - ybot)*filter_scale

            call pgsci(filter_clrs(mod(i - 1, num_filter_clrs) + 1))
            call pgline(k, filter_x, filter_y)

            i_peak = maxloc(cs% filters(i)% transmission, dim=1)
            if (filter_x(i_peak) > xleft .and. filter_x(i_peak) < xright) then
               filter_name = cs% filters(i)% name
               k = len_trim(filter_name)
               if (k > 4) then
                  if (filter_name(k - 3:k) == '.dat') filter_name(k - 3:k) = ''
               end if
               call pgsch(0.7*txt_scale)
               call pgptxt(filter_x(i_peak), ytop - 0.08*(ytop - ybot), 0.0, 0.5, trim(filter_name))
               call pgsch(txt_scale)
            end if

            deallocate(filter_x, filter_y)
         end do
         call pgsls(Line_Type_Solid)
      end if

      call pgsci(clr_Foreground)
      call pgunsa

      deallocate(xvec, yvec, in_range)

   contains

      subroutine show_message(msg)
         character(len=*), intent(in) :: msg
         call pgswin(0.0, 1.0, 0.0, 1.0)
         call pgsci(clr_Foreground)
         call show_box_pgstar(s, 'BC', 'BC')
         call pgptxt(0.5, 0.5, 0.0, 0.5, msg)
         call pgunsa
      end subroutine show_message

   end subroutine do_SED_plot

end module pgstar_sed
