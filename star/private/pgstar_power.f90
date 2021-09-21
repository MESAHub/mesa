! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

      module pgstar_power

      use star_private_def
      use const_def
      use pgstar_support

      implicit none


      contains


      subroutine power_plot(id, device_id, ierr)
         implicit none
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr

         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_power_plot(s, id, device_id, &
            s% Power_xleft, s% Power_xright, &
            s% Power_ybot, s% Power_ytop, .false., &
            s% Power_title, s% Power_txt_scale, ierr)

         call pgebuf()

      end subroutine power_plot


      subroutine do_power_plot(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         integer, intent(out) :: ierr
         call do_power_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, s% Power_xaxis_name, &
            s% Power_xmin, s% Power_xmax, s% Power_xaxis_reversed, &
            s% Power_ymin, s% Power_ymax, .false., .true., ierr)
      end subroutine do_power_plot


      subroutine do_power_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            xaxis_name, xaxis_min, xaxis_max, xaxis_reversed, ymin_in, ymax_in, &
            panel_flag, xaxis_numeric_labels_flag, ierr)
         use utils_lib
         use chem_def
         use net_def
         use const_def, only: Msun, Rsun
         implicit none

         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: &
            winxmin, winxmax, winymin, winymax, xaxis_min, xaxis_max, &
            ymin_in, ymax_in
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title, xaxis_name
         real, intent(in) :: txt_scale
         logical, intent(in) :: &
            xaxis_reversed, panel_flag, xaxis_numeric_labels_flag
         integer, intent(out) :: ierr

         character (len=strlen) :: str
         real, allocatable, dimension(:) :: xvec, yvec
         real :: xmin, xmax, xleft, xright, dx, chScale, windy, &
            ymin, ymax, exp10_ymin, xmargin, &
            legend_xmin, legend_xmax, legend_ymin, legend_ymax
         integer :: lw, lw_sav, grid_min, grid_max, npts, nz
         integer, parameter :: num_colors = 20
         integer :: colors(num_colors)

         include 'formats'
         ierr = 0
         nz = s% nz

         colors(:) = (/ &
               clr_MediumSlateBlue, clr_LightSkyBlue, clr_Goldenrod, clr_Lilac, &
               clr_Coral, clr_Crimson, clr_LightSkyGreen, clr_DarkGray, &
               clr_Tan, clr_IndianRed, clr_Gold, &
               clr_Teal, clr_Silver, clr_BrightBlue, clr_FireBrick, &
               clr_RoyalPurple, clr_SlateGray, clr_LightSteelBlue, &
               clr_Gray, clr_RoyalBlue /)

         chScale = txt_scale

         windy = winymax - winymin

         legend_xmin = winxmax
         legend_xmax = min(1.0, winxmax + (winxmax - winxmin)/5)
         legend_ymin = winymin
         legend_ymax = winymax

         allocate (xvec(nz), yvec(nz))

         xmargin = 0
         call set_xaxis_bounds(s, xaxis_name, xaxis_min, xaxis_max, &
            xaxis_reversed, xmargin, xvec, xmin, xmax, xleft, xright, dx, &
            grid_min, grid_max, npts, ierr)

         call pgsave
         call pgsch(txt_scale)
         if (ierr == 0) call plot(ierr)
         call pgunsa

         deallocate(xvec, yvec)

         contains


         subroutine plot(ierr)
            use rates_def
            integer, intent(out) :: ierr

            integer :: ii, jj, i, j, cnt
            logical, parameter :: dbg = .false.
            real(dp) :: max_power(num_categories), max_power_copy(num_categories)
            real :: ybot

            include 'formats'

            ymax = -1e9
            do i = 1, num_categories
               max_power(i) = maxval(s% eps_nuc_categories(i,grid_min:grid_max))
               if (max_power(i) > ymax) ymax = max_power(i)
               max_power_copy(i) = max_power(i)
            end do

            if (ymax < 1e-29) ymax = 1e-29
            ymax = log10(dble(ymax))
            if (ymax <= 4) then
               ymax = 4.3
               ymin = -4.1
            else
               ymax = ymax*1.1
               ymin = -0.1
            end if

            if (ymax_in /= -101) ymax = ymax_in
            if (ymin_in /= -101) ymin = ymin_in

            exp10_ymin = exp10(dble(ymin))

            lw = s% pgstar_lw
            call pgqlw(lw_sav)

            call pgsave
            call pgsvp(legend_xmin, legend_xmax, legend_ymin, legend_ymax)
            call pgswin(0.0, 1.0, ymin, ymax)
            cnt = 0
            do j=1,num_categories
               i = maxloc(max_power(:),dim=1)
               cnt = power_line_legend(cnt,i)
               max_power(i) = -1d99
            end do
            call pgunsa

            call pgsvp(winxmin, winxmax, winymin, winymax)
            if (.not. panel_flag) then
               if (.not. subplot) then
                  call show_model_number_pgstar(s)
                  call show_age_pgstar(s)
               end if
               call show_title_pgstar(s, title)
               call pgsci(1)
               call show_xaxis_name(s,xaxis_name,ierr)
               if (ierr /= 0) return
            end if

            ybot = -0.05
            call pgswin(xleft, xright, ymin+ybot, ymax)
            call pgscf(1)
            call pgsci(1)
            if (xaxis_numeric_labels_flag) then
               call show_box_pgstar(s,'BCNST','BCNSTV')
            else
               call show_box_pgstar(s,'BCST','BCNSTV')
            end if
            call show_left_yaxis_label_pgstar(s,'log ergs/g/s')

            call pgslw(lw)
            cnt = 0
            do j=1,num_categories
               i = maxloc(max_power_copy(:),dim=1)
               cnt = power_line(cnt,i)
               max_power_copy(i) = -1d99
            end do
            call pgslw(lw_sav)

            if (.not. panel_flag) then ! show mix regions at bottom of plot
               call pgslw(10)
               call show_mix_regions_on_xaxis( &
                  s,ymin+ybot,ymax,grid_min,grid_max,xvec)
            end if

         call show_pgstar_decorator(s%id,s% power_use_decorator,s% power_pgstar_decorator, 0, ierr)


         end subroutine plot


         integer function power_line(cnt, icat)
            integer, intent(in) :: cnt, icat
            real :: ymx, xpos, dx, ypos, xpts(2), ypts(2)
            integer :: iclr, k
            power_line = cnt
            ymx = maxval(s% eps_nuc_categories(icat,grid_min:grid_max))
            if (ymx < exp10_ymin) return
            iclr = modulo(icat,num_colors) +1
            power_line = cnt + 1
            call pgsci(colors(iclr))
            do k=1,nz
               yvec(k) = safe_log10(s% eps_nuc_categories(icat,k))
            end do
            call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
         end function power_line


         integer function power_line_legend(cnt, icat)
            integer, intent(in) :: cnt, icat
            real :: ymx, dx, dyline, ypos, xpts(2), ypts(2)
            integer :: iclr, num_max
            num_max = min(num_categories, s% Power_legend_max_cnt)
            power_line_legend = cnt
            if (cnt >= num_max) return
            ymx = maxval(s% eps_nuc_categories(icat,grid_min:grid_max))
            if (ymx < exp10_ymin) return
            iclr = modulo(icat,num_colors) +1
            power_line_legend = cnt + 1
            call pgsci(colors(iclr))
            dx = 0.1
            dyline = (ymax-ymin)/num_max
            ypos = ymax - (cnt+0.5)*dyline
            xpts(1) = 1.3*dx
            xpts(2) = xpts(1) + 2.3*dx
            ypts = ypos + dyline*0.1
            call pgslw(lw)
            call pgline(2, xpts, ypts)
            call pgslw(lw_sav)
            call pgsci(1)
            call pgsch(txt_scale*s% Power_legend_txt_scale_factor)
            call pgptxt(xpts(2) + dx, ypos, 0.0, 0.0, &
               trim(adjustl(category_name(icat))))
            power_line_legend = cnt + 1
         end function power_line_legend


      end subroutine do_power_panel


      end module pgstar_power

