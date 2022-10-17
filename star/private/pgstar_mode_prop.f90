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

module pgstar_mode_prop
   
   use star_private_def
   use const_def
   use pgstar_support
   use star_pgstar
   
   implicit none


contains
   
   
   subroutine mode_propagation_plot(id, device_id, ierr)
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
      
      call do_mode_propagation_plot(s, id, device_id, &
         s% pg% Mode_Prop_xleft, s% pg% Mode_Prop_xright, &
         s% pg% Mode_Prop_ybot, s% pg% Mode_Prop_ytop, .false., &
         s% pg% Mode_Prop_title, s% pg% Mode_Prop_txt_scale, ierr)
      
      call pgebuf()
   
   end subroutine mode_propagation_plot
   
   
   subroutine do_mode_propagation_plot(s, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      real, intent(in) :: txt_scale
      integer, intent(out) :: ierr
      call do_mode_propagation_panel(s, id, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, &
         title, txt_scale, s% pg% Mode_Prop_xaxis_name, &
         s% pg% Mode_Prop_xmin, s% pg% Mode_Prop_xmax, &
         s% pg% Mode_Prop_xaxis_reversed, .false., .true., ierr)
   end subroutine do_mode_propagation_plot
   
   
   subroutine do_mode_propagation_panel(s, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
      xaxis_name, xaxis_min, xaxis_max, xaxis_reversed, &
      panel_flag, xaxis_numeric_labels_flag, ierr)
      use utils_lib
      use chem_def
      use net_def
      use const_def, only : Msun, Rsun
      
      type (star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: &
         winxmin, winxmax, winymin, winymax, xaxis_min, xaxis_max
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title, xaxis_name
      real, intent(in) :: txt_scale
      logical, intent(in) :: &
         xaxis_reversed, panel_flag, xaxis_numeric_labels_flag
      integer, intent(out) :: ierr
      
      character (len = strlen) :: str
      real, allocatable, dimension(:) :: xvec, log_brunt_nu, &
         log_lamb_Sl1, log_lamb_Sl2, log_lamb_Sl3, temp_vec
      real :: xmin, xmax, xleft, xright, dx, chScale, windy, &
         ymin, ymax, exp10_ymin, xmargin, &
         legend_xmin, legend_xmax, legend_ymin, legend_ymax
      integer :: lw, lw_sav, grid_min, grid_max, npts, nz
      integer, parameter :: num_colors = 20
      integer :: colors(num_colors)
      
      include 'formats'
      
      ierr = 0
      if (.not. s% calculate_Brunt_N2) &
         call mesa_error(__FILE__, __LINE__, 'pgstar mode_propagation: must have calculate_Brunt_N2 = .true.')
      
      nz = s% nz
      
      colors(:) = (/ &
         clr_MediumSlateBlue, clr_Goldenrod, clr_LightSkyBlue, clr_Lilac, &
            clr_Coral, clr_Crimson, clr_LightSkyGreen, clr_DarkGray, &
            clr_Tan, clr_IndianRed, clr_Gold, &
            clr_Teal, clr_Silver, clr_BrightBlue, clr_FireBrick, &
            clr_RoyalPurple, clr_SlateGray, clr_LightSteelBlue, &
            clr_Gray, clr_RoyalBlue /)
      
      chScale = txt_scale
      
      windy = winymax - winymin
      
      legend_xmin = winxmax - 0.01
      legend_xmax = 0.99
      legend_ymin = winymin
      legend_ymax = winymax
      
      allocate (xvec(nz), log_brunt_nu(nz), &
         log_lamb_Sl1(nz), log_lamb_Sl2(nz), log_lamb_Sl3(nz), temp_vec(nz))
      
      xmargin = 0
      call set_xaxis_bounds(&
         s, xaxis_name, xaxis_min, xaxis_max, xaxis_reversed, xmargin, &
         xvec, xmin, xmax, xleft, xright, dx, &
         grid_min, grid_max, npts, ierr)
      if (ierr == 0) then
         call pgsave
         call pgsch(txt_scale)
         call plot(ierr)
         call pgunsa
      end if
      
      deallocate(xvec, log_brunt_nu, &
         log_lamb_Sl1, log_lamb_Sl2, log_lamb_Sl3, temp_vec)
   
   contains
      
      
      subroutine plot(ierr)
         use rates_def
         integer, intent(out) :: ierr
         
         integer :: ii, jj, i, cnt, k
         logical, parameter :: dbg = .false.
         real :: ybot, nu_max, lg_nu_max, lg_2pt0_nu_max, lg_0pt5_nu_max, lg_nu_max_obs
         real, parameter :: teff_sun = 5777.0, nu_max_sun = 3100.0
         
         include 'formats'
         
         do k = grid_min, grid_max
            log_brunt_nu(k) = safe_log10((1d6 / (2 * pi)) * sqrt(max(0d0, s% brunt_N2(k))))
            log_lamb_Sl1(k) = safe_log10((1d6 / (2 * pi)) * sqrt(2d0) * s% csound_face(k) / s% r(k))
            log_lamb_Sl2(k) = safe_log10((1d6 / (2 * pi)) * sqrt(6d0) * s% csound_face(k) / s% r(k))
            log_lamb_Sl3(k) = safe_log10((1d6 / (2 * pi)) * sqrt(12d0) * s% csound_face(k) / s% r(k))
         end do
         
         nu_max = nu_max_sun * s% star_mass / (pow2(s% photosphere_r) * sqrt(s% Teff / teff_sun))
         lg_nu_max = log10(dble(nu_max))
         lg_2pt0_nu_max = log10(dble(2.0 * nu_max))
         lg_0pt5_nu_max = log10(dble(0.5 * nu_max))
         lg_nu_max_obs = safe_log10(dble(s% pg% Mode_Prop_nu_max_obs))
         
         ymax = max(1.33 * lg_2pt0_nu_max, maxval(log_brunt_nu(grid_min:grid_max)))
         ymin = 0.5 * lg_0pt5_nu_max
         ymax = ymax + 0.1 * (ymax - ymin)
         
         if (s% pg% Mode_Prop_ymax /= -101) ymax = s% pg% Mode_Prop_ymax
         if (s% pg% Mode_Prop_ymin /= -101) ymin = s% pg% Mode_Prop_ymin
         
         lw = s% pg% pgstar_lw
         call pgqlw(lw_sav)
         
         call pgsave
         call pgsvp(legend_xmin, legend_xmax, legend_ymin, legend_ymax)
         call pgswin(0.0, 1.0, ymin, ymax)
         cnt = 0
         cnt = mode_propagation_line_legend(cnt, 'N\dBV\u')
         cnt = mode_propagation_line_legend(cnt, 'S\dl=1\u')
         cnt = mode_propagation_line_legend(cnt, 'S\dl=2\u')
         cnt = mode_propagation_line_legend(cnt, 'S\dl=3\u')
         cnt = mode_propagation_line_legend(cnt, '2\(2723)\(2139)\dmax\u')
         cnt = mode_propagation_line_legend(cnt, '\(2139)\dmax\u')
         call pgsls(4) ! dotted
         cnt = mode_propagation_line_legend(cnt, '\(2139)\dmax\uobs')
         call pgsls(1) ! solid
         cnt = mode_propagation_line_legend(cnt, '0.5\(2723)\(2139)\dmax\u')
         call pgunsa
         
         call pgsave
         call pgsvp(winxmin, winxmax, winymin, winymax)
         if (.not. panel_flag) then
            if (.not. subplot) then
               call show_model_number_pgstar(s)
               call show_age_pgstar(s)
            end if
            call show_title_pgstar(s, title)
         end if
         
         ybot = -0.05
         call pgswin(xleft, xright, ymin + ybot, ymax)
         call pgscf(1)
         call pgsci(1)
         if (xaxis_numeric_labels_flag) then
            call show_box_pgstar(s, 'BCNST', 'BCNSTV')
         else
            call show_box_pgstar(s, 'BCST', 'BCNSTV')
         end if
         call show_left_yaxis_label_pgstar(s, 'log \(2139) (\(2138)Hz)')
         
         call pgslw(lw)
         cnt = 0
         cnt = mode_propagation_line(cnt, log_brunt_nu)
         cnt = mode_propagation_line(cnt, log_lamb_Sl1)
         cnt = mode_propagation_line(cnt, log_lamb_Sl2)
         cnt = mode_propagation_line(cnt, log_lamb_Sl3)
         temp_vec(1:nz) = lg_2pt0_nu_max
         cnt = mode_propagation_line(cnt, temp_vec)
         temp_vec(1:nz) = lg_nu_max
         cnt = mode_propagation_line(cnt, temp_vec)
         call pgsls(4) ! dotted
         temp_vec(1:nz) = lg_nu_max_obs
         cnt = mode_propagation_line(cnt, temp_vec)
         call pgsls(1) ! solid
         temp_vec(1:nz) = lg_0pt5_nu_max
         cnt = mode_propagation_line(cnt, temp_vec)
         call pgslw(lw_sav)
         
         if (.not. panel_flag) then
            call pgsci(1)
            call show_xaxis_name(s, xaxis_name, ierr)
            if (ierr == 0) then ! show mix regions at bottom of plot
               call pgslw(10)
               call show_mix_regions_on_xaxis(&
                  s, ymin + ybot, ymax, grid_min, grid_max, xvec)
            end if
         end if
         
         call pgunsa
         
         call show_pgstar_decorator(s%id, s% pg% mode_prop_use_decorator, &
            s% pg% mode_prop_pgstar_decorator, 0, ierr)
      
      end subroutine plot
      
      
      integer function mode_propagation_line(cnt, yvec)
         integer, intent(in) :: cnt
         real, intent(in) :: yvec(:)
         integer :: iclr
         iclr = cnt - num_colors * (cnt / num_colors) + 1
         mode_propagation_line = cnt + 1
         call pgsci(colors(iclr))
         call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
      end function mode_propagation_line
      
      
      integer function mode_propagation_line_legend(cnt, name)
         integer, intent(in) :: cnt
         character (len = *), intent(in) :: name
         real :: dx, dyline, ypos, xpts(2), ypts(2)
         integer :: iclr, num_max
         num_max = 10
         mode_propagation_line_legend = cnt
         iclr = cnt - num_colors * (cnt / num_colors) + 1
         call pgsci(colors(iclr))
         dx = 0.1
         dyline = (ymax - ymin) / num_max
         ypos = ymax - (cnt + 1.5) * dyline
         xpts(1) = 1.3 * dx
         xpts(2) = xpts(1) + 2.3 * dx
         ypts = ypos + dyline * 0.1
         call pgslw(lw)
         call pgline(2, xpts, ypts)
         call pgslw(lw_sav)
         call pgsci(1)
         call pgsch(txt_scale * 0.70)
         call pgptxt(xpts(2) + dx, ypos, 0.0, 0.0, name)
         mode_propagation_line_legend = cnt + 1
      end function mode_propagation_line_legend
   
   
   end subroutine do_mode_propagation_panel


end module pgstar_mode_prop

