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

      module pgstar_abundance

      use star_private_def
      use const_def
      use pgstar_support
      use star_pgstar

      implicit none


      contains


      subroutine abundance_plot(id, device_id, ierr)
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

         call do_abundance_plot(s, id, device_id, &
            s% pg% Abundance_xleft, s% pg% Abundance_xright, &
            s% pg% Abundance_ybot, s% pg% Abundance_ytop, .false., &
            s% pg% Abundance_title, s% pg% Abundance_txt_scale, ierr)

         call pgebuf()

      end subroutine abundance_plot


      subroutine do_abundance_plot(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         integer, intent(out) :: ierr
         call do_abundance_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, s% pg% Abundance_xaxis_name, &
            s% pg% Abundance_xmin, s% pg% Abundance_xmax, s% pg% Abundance_xaxis_reversed, &
            s% pg% abundance_log_mass_frac_min, s% pg% abundance_log_mass_frac_max, &
            .false., .true., ierr)
      end subroutine do_abundance_plot


      subroutine do_abundance_panel(s, id, device_id, &
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
            winxmin, winxmax, winymin, winymax, xaxis_min, xaxis_max, ymin_in, ymax_in
         character (len=*), intent(in) :: title, xaxis_name
         real, intent(in) :: txt_scale
         logical, intent(in) :: subplot, &
            xaxis_reversed, panel_flag, xaxis_numeric_labels_flag
         integer, intent(out) :: ierr

         character (len=strlen) :: str
         real, allocatable, dimension(:) :: xvec, yvec
         real :: xmin, xmax, xleft, xright, dx, dylbl, chScale, windy, xmargin, &
            ymin, ymax, legend_xmin, legend_xmax, legend_ymin, legend_ymax
         integer :: lw, lw_sav, grid_min, grid_max, npts, i, nz
         integer, parameter :: num_colors = 14
         integer :: colors(num_colors)
         integer, parameter :: max_num_labels = 30
         integer :: num_labels, iloc_abundance_label(max_num_labels)
         real :: xloc_abundance_label(max_num_labels)

         include 'formats'
         ierr = 0
         nz = s% nz

         colors(:) = (/ &
               clr_Gold, clr_LightSkyBlue, clr_Crimson, clr_Goldenrod, clr_MediumSlateBlue, &
               clr_Coral, clr_LightSkyGreen, clr_DarkGray, clr_Lilac, &
               clr_Tan, clr_IndianRed, clr_Teal, clr_Silver, clr_BrightBlue /)

         chScale = txt_scale

         ymin = ymin_in
         ymax = ymax_in
         if (abs(ymin+101.0) < 1e-6) ymin = s% pg% abundance_log_mass_frac_min
         if (abs(ymax+101.0) < 1e-6) ymax = s% pg% abundance_log_mass_frac_max

         windy = winymax - winymin

         legend_xmin = winxmax
         legend_xmax = min(1.0, winxmax + (winxmax - winxmin)/5)

         legend_ymin = winymin
         legend_ymax = winymax

         allocate(xvec(nz), yvec(nz))

         xmargin = 0
         call set_xaxis_bounds( &
            s, xaxis_name, xaxis_min, xaxis_max, xaxis_reversed, &
            xmargin, xvec, xmin, xmax, xleft, xright, dx, &
            grid_min, grid_max, npts, ierr)
         if (ierr == 0) call plot(ierr)

         deallocate(xvec, yvec)

         contains


         subroutine plot(ierr)
            use rates_def
            integer, intent(out) :: ierr

            integer :: ii, jj, i, k, xaxis_id
            logical, parameter :: dbg = .false.
            logical :: found_shock
            real(dp) :: xshock, photosphere_logxm
            real :: lgz, x, y, ybot

            include 'formats'
            ierr = 0

            if (ymin > 0) then
               ymin = -5.1
               lgz = log10(s% ctrl% initial_z + 1e-9) - 1
               if (lgz-1 < ymin) ymin = lgz
            end if
            if (ymax >= 100) then
               if (ymin < -8) then
                  ymax = 0.5
               else
                  ymax = 0.25
               end if
            end if

            num_labels = max(0,min(max_num_labels, s% pg% num_abundance_line_labels))
            
            iloc_abundance_label = -HUGE(grid_min)
            xloc_abundance_label = -HUGE(grid_min)
            do i=1,num_labels
               x = xmin + (i-0.5)*dx/num_labels
               do k=2,nz
                  if ((xvec(k-1)-x)*(x-xvec(k)) >= 0) then
                     iloc_abundance_label(i) = k
                     xloc_abundance_label(i) = x
                     exit
                  end if
               end do
            end do

            dylbl = (ymax - ymin)*0.015

            lw = s% pg% pgstar_lw
            call pgqlw(lw_sav)

            call pgsave
            call pgsch(txt_scale)
            call pgsvp(legend_xmin, legend_xmax, legend_ymin, legend_ymax)
            call pgswin(0.0, 1.0, ymin, ymax)
            call do_all(.true.)
            call pgunsa

            call pgsave
            call pgsch(txt_scale)

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
            call show_left_yaxis_label_pgstar(s, 'log mass fraction')

            call pgsave
            call pgsch(txt_scale*1.05)
            call do_all(.false.)
            call pgunsa

            if (.not. panel_flag) then ! show mix regions at bottom of plot
               call pgslw(10)
               call show_mix_regions_on_xaxis( &
                  s,ymin+ybot,ymax,grid_min,grid_max,xvec)
            end if

            call pgunsa
            
            if (s% pg% Abundance_show_photosphere_location .and. &
                  (xaxis_name == 'mass' .or. &
                   xaxis_name == 'logxm' .or. &
                   xaxis_name == 'radius')) then
               call pgsave
               call pgsvp(winxmin, winxmax, winymin, winymax)
               call pgswin(0.0, 1.0, 0.0, 1.0)
               call pgsci(clr_Gray)
               call pgsls(Line_Type_Dash)
               call pgslw(5)
               if (xaxis_name == 'radius') then
                  dx = (s% photosphere_r - xmin)/(xmax - xmin)
               else if (xaxis_name == 'radius_cm') then
                  dx = (s% photosphere_r*Rsun - xmin)/(xmax - xmin)
               else if (xaxis_name == 'mass') then
                  dx = (s% photosphere_m - xmin)/(xmax - xmin)
               else
                  if (s% star_mass > s% photosphere_m) then
                     photosphere_logxm = log10(s% star_mass - s% photosphere_m)
                  else
                     photosphere_logxm = -99
                  end if
                  !dx = (xmax - photosphere_logxm)/(xmin - xmax)
                  dx = 1 - (photosphere_logxm - xmin)/(xmax - xmin)
                  write(*,2) 'photosphere_xm xmin photo_x xmax dx', s% model_number, &
                     s% star_mass - s% photosphere_m, &
                     xmin, photosphere_logxm, xmax, dx
               end if
               call pgmove(dx, 0.0)
               call pgdraw(dx, 1.0)
               call pgunsa
            end if
            
         call show_pgstar_decorator(s%id,s% pg% Abundance_use_decorator,s% pg% Abundance_pgstar_decorator,0, ierr)

         end subroutine plot


         subroutine do_all(legend_flag)
            logical, intent(in) :: legend_flag
            integer :: cnt, num_to_show, i, j, k, jmax
            real(dp) :: max_abund(s% species)
            include 'formats'
            cnt = 0
            num_to_show = s% pg% Abundance_num_isos_to_show
            do j=1, s% species
               max_abund(j) = maxval(s% xa(j,grid_min:grid_max))
            end do
            if (num_to_show < 0) then ! show as many as fit
               if (legend_flag) then
                  jmax = min(s% species, max_num_labels)
               else
                  jmax = s% species
               end if
               do j = 1, jmax
                  i = maxloc(max_abund(:),dim=1)
                  cnt = do1(cnt, chem_isos% name(s% chem_id(i)), legend_flag)
                  max_abund(i) = -1
               end do
            else
               do i = 1, num_to_show
                  cnt = do1(cnt, s% pg% Abundance_which_isos_to_show(i), legend_flag)
               end do
            end if
         end subroutine do_all


         integer function do1(cnt, str, legend_flag)
            use chem_lib
            integer, intent(in) :: cnt
            character (len=*), intent(in) :: str
            logical, intent(in) :: legend_flag
            integer :: i, cid, k
            include 'formats'
            do1 = cnt
            if (len_trim(str) == 0) return
            cid = chem_get_iso_id(str)
            if (cid <= 0) return
            i = s% net_iso(cid)
            if (i == 0) return
            do k=1,nz
               yvec(k) = safe_log10(s% xa(i,k))
            end do
            if (s% pg% abundance_log_mass_frac_min < 0) then
               if (maxval(yvec(grid_min:grid_max)) < &
                     s% pg% abundance_log_mass_frac_min) return
            end if
            if (legend_flag) then
               do1 = abundance_line_legend(cnt, str)
            else
               do1 = abundance_line(cnt, i, str)
            end if
         end function do1


         subroutine set_line_style(cnt)
            integer, intent(in) :: cnt
            integer :: iclr, itype
            iclr = cnt - num_colors*(cnt/num_colors) + 1
            call pgsci(colors(iclr))
            if (cnt >= num_colors) then
               itype = Line_Type_Dot
            else
               itype = Line_Type_Solid
            end if
            call pgsls(itype)
         end subroutine set_line_style


         integer function abundance_line(cnt, j, str)
            integer, intent(in) :: cnt, j
            character (len=*), intent(in) :: str
            integer :: i, ii
            real :: x, frac, y, dx
            include 'formats'
            call set_line_style(cnt)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
            call pgslw(lw_sav)
            call pgsch(txt_scale*s% pg% Abundance_line_txt_scale_factor)
            abundance_line = cnt + 1
            do i=1,num_labels
               ii = iloc_abundance_label(i)
               if (ii > grid_max .or. ii < grid_min) cycle
               x = xloc_abundance_label(i)
               if (ii > grid_min) then
                  dx = xvec(ii)-xvec(ii-1)
                  if (abs(dx) > 1e-20) then
                     frac = (x-xvec(ii-1))/dx
                  else
                     frac = 1.0
                  end if
                  y = frac*yvec(ii) + (1.0-frac)*yvec(ii-1) + dylbl
               else
                  y = yvec(ii) + dylbl
               end if
               if (y < ymin .or. y > ymax) cycle
               call pgptxt(x, y, 0.0, 0.5, str)
            end do
         end function abundance_line


         integer function abundance_line_legend(cnt, str)
            integer, intent(in) :: cnt
            character (len=*), intent(in) :: str
            real :: ymx, dx, dyline, ypos, xpts(2), ypts(2)
            integer :: iclr, max_cnt
            max_cnt = min(max_num_labels, s% pg% Abundance_legend_max_cnt)
            if (cnt >= max_cnt) then
               abundance_line_legend = cnt
               return
            end if
            call set_line_style(cnt)
            dx = 0.1
            dyline = (ymax-ymin)/max_cnt
            ypos = ymax - (cnt+0.5)*dyline
            xpts(1) = 1.3*dx
            xpts(2) = xpts(1) + 2.3*dx
            ypts = ypos + dyline*0.1
            call pgslw(lw)
            call pgline(2, xpts, ypts)
            call pgslw(lw_sav)
            call pgsci(1)
            call pgsch(txt_scale*s% pg% Abundance_legend_txt_scale_factor)
            call pgptxt(xpts(2) + dx, ypos, 0.0, 0.0, trim(str))
            abundance_line_legend = cnt + 1
         end function abundance_line_legend


      end subroutine do_abundance_panel


      end module pgstar_abundance

