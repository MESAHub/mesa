! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module pgstar_HR

      use star_private_def
      use const_def, only: dp
      use pgstar_support
      use star_pgstar

      implicit none


      contains


      subroutine HR_Plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_HR_Plot(s, id, device_id, &
            s% pg% HR_xleft, s% pg% HR_xright, &
            s% pg% HR_ybot, s% pg% HR_ytop, .false., &
            s% pg% HR_title, s% pg% HR_txt_scale, ierr)
         if (ierr /= 0) return

         call pgebuf()

      end subroutine HR_Plot


      subroutine HR_decorate(id, ierr)
         use pgstar_colors
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real :: logT1, logT2, logL1, logL2
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% pg% show_HR_radius_lines) call HR_radius_lines()
         if (s% pg% show_HR_Mira_instability_region) then
            logL1 = 2.8
            logL2 = 4.0
            logT2 = 3.7
            logT1 = 3.45

            call pgsls(Line_Type_Solid)
            call pgslw(s% pg% pgstar_lw)
            call pgsci(clr_Blue)
            call pgmove(logT1, logL1)
            call pgdraw(logT1, logL2)
            call pgdraw(logT2, logL2)
            call pgdraw(logT2, logL1)
            call pgdraw(logT1, logL1)
         end if
         if (s% pg% show_HR_classical_instability_strip) then

            call pgsls(Line_Type_Solid)
            call pgslw(s% pg% pgstar_lw)

            ! approximate edges

            ! blue edge
            logT1 = 3.70
            logL1 = 5.5

            logT2 = 3.93
            logL2 = 1.0

            call pgsci(clr_Blue)
            call pgmove(logT1, logL1)
            call pgdraw(logT2, logL2)

            ! red edge
            logT1 = 3.60
            logL1 = 5.5

            logT2 = 3.83
            logL2 = 1.0

            call pgsci(clr_FireBrick)
            call pgmove(logT1, logL1)
            call pgdraw(logT2, logL2)
         end if
         if (s% pg% show_HR_WD_instabilities) then
            ! from Winget & Kepler, Annu. Rev. Astron. Astrophys., 2008, Fig 3.
            call pgsave
            call pgsls(Line_Type_Solid)
            call pgslw(s% pg% pgstar_lw)
            call pgsci(clr_Silver)
            call pgmove(5.1, 4.3)  ! DOV
            call pgdraw(4.8, 0.3)
            call pgdraw(4.95, 0.9)
            call pgdraw(5.38, 4.2)
            call pgdraw(5.1, 4.3)
            call pgmove(4.42, -0.7)  ! DBV
            call pgdraw(4.34, -1.2)
            call pgdraw(4.38, -1.7)
            call pgdraw(4.45, -1.4)
            call pgdraw(4.42, -0.7)
            call pgmove(4.03, -2.4)  ! DAV
            call pgdraw(4.03, -3.4)
            call pgdraw(4.1, -3.2)
            call pgdraw(4.1, -2.2)
            call pgdraw(4.03, -2.4)
            call pgsch(s% pg% HR_txt_scale*0.7)
            call pgslw(1)
            call pgptxt(5.1 - 0.05, 4.3, 0.0, 0.0, 'DOV')
            call pgptxt(4.42 - 0.05, -0.7, 0.0, 0.0, 'DBV')
            call pgptxt(4.03 - 0.05, -2.4, 0.0, 0.0, 'DAV')
            call pgunsa
         end if
      end subroutine HR_decorate


      subroutine HR_radius_lines()
         use const_def, only: pi, boltz_sigma, Rsun, Lsun
         use math_lib, only: log10
         use pgstar_colors, only: clr_Gray
         character(len=48) :: label
         character(len=16) :: radius
         real :: xmin, xmax, ymin, ymax
         real :: x(2), y(2), xlo, xhi, ylo, yhi, dx, dy
         real :: intercept, label_lo, label_hi, label_x, char_height
         real :: solar_intercept
         integer :: i, imin, imax

         call pgqwin(xmin, xmax, ymin, ymax)
         ! HR temperature axes run from hot to cool; handle either ordering.
         xlo = min(xmin, xmax)
         xhi = max(xmin, xmax)
         ylo = min(ymin, ymax)
         yhi = max(ymin, ymax)
         dx = xhi - xlo
         dy = yhi - ylo
         if (dx <= 0.0 .or. dy <= 0.0) return

         ! save the current drawing settings -- we'll restore them later
         ! so other decorations don't inherit this style
         call pgsave

         call pgsci(clr_Gray)
         call pgsls(Line_Type_Dash)
         call pgslw(1)
         call pgqch(char_height)
         ! do_Hist_Track leaves the multiplied the character size by 2.8, so undo
         ! that here.
         call pgsch(0.7*char_height/2.8)

         ! log10(L/Lsun) = 4*log10(Teff) + 2*log10(R/Rsun) + solar_intercept.
         solar_intercept = real(log10(4.0_dp*pi*boltz_sigma*pow2(Rsun)/Lsun)) ! cast to real b/c pgstar uses single precision

         ! figure out which radii (powers of 10) are visible
         ! smallest visible line is at highest Teff and lowest L (xhi, ylo)
         ! largest visible line is at lowest Teff and highest L (xlo, yhi)
         imin = ceiling(0.5*(ylo - 4.0*xhi - solar_intercept))
         imax = floor(0.5*(yhi - 4.0*xlo - solar_intercept))

         do i = imin, imax
            ! get the intercept for the current radius
            intercept = solar_intercept + 2.0*i

            ! Clip the line to the current axes before drawing it.
            x(1) = max(xlo, (ylo-intercept)/4.0)
            x(2) = min(xhi, (yhi-intercept)/4.0)
            if (x(1) >= x(2)) cycle
            y = 4.0*x + intercept
            call pgline(2, x, y)

            ! Label an interior segment so labels follow changing axis limits.
            label_lo = max(xlo + 0.12*dx, (ylo + 0.08*dy-intercept)/4.0)
            label_hi = min(xhi - 0.12*dx, (yhi - 0.12*dy-intercept)/4.0)

            if (label_lo >= label_hi) cycle

            label_x = 0.5*(label_lo + label_hi)

            if (i >= 0 .and. i <= 3) then
               write(radius, '(i0)') 10**i                    ! 1, 10, 100, 1000
               label = trim(radius) // ' R\d\(2281)\u'
            else if (i < 0 .and. i >= -3) then
               radius = '0.' // repeat('0', -i-1) // '1'      ! 0.1, 0.01, 0.001
               label = trim(radius) // ' R\d\(2281)\u'
            else
               write(radius, '(i0)') i                        ! i0 includes the minus sign
               label = '10\u' // trim(radius) // '\d R\d\(2281)\u'
            end if

            call pgptxt(label_x, 4.0*label_x + intercept + 0.02*dy, &
               0.0, 0.5, trim(label))
         end do

         ! restore the drawing settings
         call pgunsa

      end subroutine HR_radius_lines


      subroutine do_HR_Plot(s, id, device_id, &
            xleft, xright, ybot, ytop, subplot, &
            title, txt_scale, ierr)
         use pgstar_hist_track, only: do_Hist_Track
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         integer, intent(out) :: ierr
         logical, parameter :: &
            reverse_xaxis = .true., reverse_yaxis = .false.
         ierr = 0
         call do_Hist_Track(s, id, device_id, &
            xleft, xright, ybot, ytop, subplot, title, txt_scale, &
            'log_Teff', 'log_L', &
            'log Teff', 'log L/L\d\(2281)', &
            s% pg% HR_logT_min, s% pg% HR_logT_max, &
            s% pg% HR_logT_margin, s% pg% HR_dlogT_min, &
            s% pg% HR_logL_min, s% pg% HR_logL_max, &
            s% pg% HR_logL_margin, s% pg% HR_dlogL_min, &
            s% pg% HR_step_min, s% pg% HR_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            s% pg% show_HR_target_box, s% pg% HR_target_n_sigma, &
            s% pg% HR_target_logT, s% pg% HR_target_logL, &
            s% pg% HR_target_logT_sigma, s% pg% HR_target_logL_sigma, &
            s% pg% show_HR_annotation1, &
            s% pg% show_HR_annotation2, &
            s% pg% show_HR_annotation3, &
            s% pg% HR_fname, &
            s% pg% HR_use_decorator, &
            s% pg% HR_pgstar_decorator, &
            HR_decorate, ierr)
      end subroutine do_HR_Plot


      end module pgstar_HR

