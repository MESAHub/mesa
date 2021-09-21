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

      module pgstar_mixing_Ds

      use star_private_def
      use const_def
      use pgstar_support
      use pgstar_trho_profile

      implicit none


      contains


      subroutine Mixing_plot(id, device_id, ierr)
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

         call do_Mixing_plot(s, id, device_id, &
            s% Mixing_xleft, s% Mixing_xright, &
            s% Mixing_ybot, s% Mixing_ytop, .false., &
            s% Mixing_title, s% Mixing_txt_scale, ierr)

         call pgebuf()

      end subroutine Mixing_plot


      subroutine do_Mixing_plot(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         integer, intent(out) :: ierr
         call do_Mixing_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            s% Mixing_xaxis_name, s% Mixing_xmin, s% Mixing_xmax, &
            s% Mixing_xaxis_reversed, s% Mixing_ymin, s% Mixing_ymax, &
            .false., .true., ierr)
      end subroutine do_Mixing_plot


      subroutine do_Mixing_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            xaxis_name, xmin, xmax, xaxis_reversed, ymin, ymax, &
            panel_flag, xaxis_numeric_labels_flag, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: &
            winxmin, winxmax, winymin, winymax, xmin, xmax, ymin, ymax
         character (len=*), intent(in) :: title, xaxis_name
         real, intent(in) :: txt_scale
         logical, intent(in) :: subplot, &
            xaxis_reversed, panel_flag, xaxis_numeric_labels_flag
         integer, intent(out) :: ierr
         call MixDs_plot(s, device_id, &
            s% Mixing_win_flag, s% Mixing_file_flag, &
            s% do_Mixing_win, s% do_Mixing_file, &
            s% id_Mixing_win, s% id_Mixing_file, s% Mixing_file_interval, &
            s% Mixing_file_dir, s% Mixing_file_prefix, &
            s% show_Mixing_annotation1, s% show_Mixing_annotation2, &
            s% show_Mixing_annotation3, &
            xaxis_name, xmin, xmax, xaxis_reversed, &
            ymin, ymax, s% Mixing_dymin, &
            s% Mixing_win_width, s% Mixing_win_aspect_ratio, &
            s% prev_Mixing_win_width, s% prev_Mixing_win_ratio, &
            s% Mixing_file_width, s% Mixing_file_aspect_ratio, &
            s% prev_Mixing_file_width, s% prev_Mixing_file_ratio, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            panel_flag, xaxis_numeric_labels_flag, ierr)
      end subroutine do_Mixing_panel



      subroutine MixDs_plot(s, device_id, &
            MixDs_win_flag, MixDs_file_flag, &
            do_MixDs_win, do_MixDs_file, &
            id_MixDs_win, id_MixDs_file, MixDs_file_interval, &
            MixDs_file_dir, MixDs_file_prefix, &
            show_MixDs_annotation1, show_MixDs_annotation2, show_MixDs_annotation3, &
            MixDs_xaxis_name, MixDs_xmin, MixDs_xmax, MixDs_xaxis_reversed, &
            MixDs_ymin, MixDs_ymax, MixDs_dymin, &
            MixDs_win_width, MixDs_win_aspect_ratio, &
            prev_MixDs_win_width, prev_MixDs_win_ratio, &
            MixDs_file_width, MixDs_file_aspect_ratio, &
            prev_MixDs_file_width, prev_MixDs_file_ratio, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            panel_flag, xaxis_numeric_labels_flag, ierr)

         use utils_lib
         implicit none

         type (star_info), pointer :: s
         integer, intent(in) :: device_id
         logical, intent(in) :: MixDs_win_flag, MixDs_file_flag
         logical, intent(in) :: do_MixDs_win, do_MixDs_file
         integer, intent(in) :: id_MixDs_win, id_MixDs_file, MixDs_file_interval
         character (len=strlen), intent(in) :: MixDs_file_dir, MixDs_file_prefix
         logical, intent(in) :: show_MixDs_annotation1, show_MixDs_annotation2, show_MixDs_annotation3
         character (len=*), intent(in) :: MixDs_xaxis_name, title
         real, intent(in) :: &
            MixDs_xmin, MixDs_xmax, &
            MixDs_ymin, MixDs_ymax, MixDs_dymin, &
            MixDs_win_width, MixDs_win_aspect_ratio, &
            prev_MixDs_win_width, prev_MixDs_win_ratio, &
            MixDs_file_width, MixDs_file_aspect_ratio, &
            prev_MixDs_file_width, prev_MixDs_file_ratio
         real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
         logical, intent(in) :: subplot, &
            MixDs_xaxis_reversed, panel_flag, xaxis_numeric_labels_flag
         integer, intent(out) :: ierr

         real :: chScale, windy, xmargin
         real :: xmin, xmax, xleft, xright, dx, tmp, ymin, ymax, ymin2, ymax2, dy, &
            legend_xmin, legend_xmax, legend_ymin, legend_ymax
         integer :: grid_min, grid_max, npts, nz, number_of_legend_lines
         real, allocatable, dimension(:) :: &
            xvec, yvec, y_conv, y_left, y_sc, y_ovr, y_th, y_min_mix, y_RTI_mix

         include 'formats'
         ierr = 0

         xmargin = 0

         chScale = txt_scale

         nz = s% nz
         allocate (xvec(nz), yvec(nz), y_conv(nz), y_left(nz), &
            y_sc(nz), y_ovr(nz), y_th(nz), y_min_mix(nz), y_RTI_mix(nz))

         call set_xaxis_bounds( &
            s, MixDs_xaxis_name, MixDs_xmin, MixDs_xmax, &
            MixDs_xaxis_reversed, xmargin, &
            xvec, xmin, xmax, xleft, xright, dx, &
            grid_min, grid_max, npts, ierr)
         if (ierr /= 0) return

         legend_xmin = winxmax
         legend_xmax = min(1.0, winxmax + (winxmax - winxmin)/5)
         legend_ymin = winymin
         legend_ymax = winymax

         call pgsave
         call pgsch(txt_scale)
         call plot(ierr)
         if (ierr == 0) call show_annotations(s, &
            show_MixDs_annotation1, show_MixDs_annotation2, show_MixDs_annotation3)
         call pgunsa

         deallocate(xvec, yvec, y_conv, y_left, y_sc, y_ovr, y_th, y_min_mix, y_RTI_mix)

         contains


         subroutine plot(ierr)
            integer, intent(out) :: ierr

            integer :: lw, lw_sav
            real :: val
            character (len=128) :: str
            integer :: i, ii, k, cnt
            logical :: rotation
            real(dp) :: &
               D_visc_factor, &
               D_DSI_factor, &
               D_SH_factor, &
               D_SSI_factor, &
               D_ES_factor, &
               D_GSF_factor, &
               D_ST_factor

            include 'formats'
            ierr = 0

            rotation = s% rotation_flag

            lw = 8
            call pgqlw(lw_sav)

            if (.not. panel_flag) then
               call pgsvp(winxmin, winxmax, winymin, winymax)
               if (.not. subplot) then
                  call show_model_number_pgstar(s)
                  call show_age_pgstar(s)
               end if
               call show_title_pgstar(s, title)
               call pgsci(1)
               call show_xaxis_name(s,MixDs_xaxis_name,ierr)
               if (ierr /= 0) return
            end if

            y_conv(:) = -100
            y_left(:) = -100
            y_sc(:) = -100
            y_ovr(:) = -100
            y_th(:) = -100
            y_RTI_mix(:) = -100
            y_min_mix(:) = -100
            do k=grid_min, grid_max
               val = safe_log10(s% D_mix_non_rotation(k))
               select case (s% mixing_type(k))
                  case (convective_mixing)
                     y_conv(k) = val
                  case (leftover_convective_mixing)
                     y_left(k) = val
                  case (semiconvective_mixing)
                     y_sc(k) = val
                  case (overshoot_mixing)
                     y_ovr(k) = val
                  case (thermohaline_mixing)
                     y_th(k) = val
                  case (rayleigh_taylor_mixing)
                     y_RTI_mix(k) = val
                  case (minimum_mixing)
                     y_min_mix(k) = val
               end select
            end do

            if (MixDs_ymax /= -101) then
               ymax = MixDs_ymax
            else
               ymax = max(18.0,maxval(y_conv(grid_min:grid_max)))
               ymax2 = max(18.0,maxval(y_left(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
               ymax2 = max(18.0,maxval(y_sc(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
               ymax2 = max(18.0,maxval(y_ovr(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
               ymax2 = max(18.0,maxval(y_th(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
               ymax2 = max(18.0,maxval(y_RTI_mix(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
               ymax2 = max(18.0,maxval(y_min_mix(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
               if (rotation) then
                  ymax2 = max(18.0,real(safe_log10(maxval(s% D_visc(grid_min:grid_max))),kind=kind(ymax2)))
                  if (ymax2 > ymax) ymax = ymax2
                  ymax2 = max(18.0,real(safe_log10(maxval(s% D_DSI(grid_min:grid_max))),kind=kind(ymax2)))
                  if (ymax2 > ymax) ymax = ymax2
                  ymax2 = max(18.0,real(safe_log10(maxval(s% D_SH(grid_min:grid_max))),kind=kind(ymax2)))
                  if (ymax2 > ymax) ymax = ymax2
                  ymax2 = max(18.0,real(safe_log10(maxval(s% D_SSI(grid_min:grid_max))),kind=kind(ymax2)))
                  if (ymax2 > ymax) ymax = ymax2
                  ymax2 = max(18.0,real(safe_log10(maxval(s% D_ES(grid_min:grid_max))),kind=kind(ymax2)))
                  if (ymax2 > ymax) ymax = ymax2
                  ymax2 = max(18.0,real(safe_log10(maxval(s% D_GSF(grid_min:grid_max))),kind=kind(ymax2)))
                  if (ymax2 > ymax) ymax = ymax2
                  ymax2 = max(18.0,real(safe_log10(maxval(s% D_ST(grid_min:grid_max))),kind=kind(ymax2)))
                  if (ymax2 > ymax) ymax = ymax2
                  if (s% D_omega_flag) then
                     ymax2 = max(18.0,real(safe_log10(maxval(s% D_omega(grid_min:grid_max))),kind=kind(ymax2)))
                     if (ymax2 > ymax) ymax = ymax2
                  end if
               end if
            end if

            if (MixDs_ymin /= -101) then
               ymin = MixDs_ymin
            else
               ymin = max(0.0,minval(y_conv(grid_min:grid_max)))
               ymin2 = max(0.0,minval(y_left(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
               ymin2 = max(0.0,minval(y_sc(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
               ymin2 = max(0.0,minval(y_ovr(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
               ymin2 = max(0.0,minval(y_th(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
               ymin2 = max(0.0,minval(y_RTI_mix(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
               ymin2 = max(0.0,minval(y_min_mix(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
               if (rotation) then
                  ymin2 = max(0.0,real(safe_log10(minval(s% D_visc(grid_min:grid_max))),kind=kind(ymin2)))
                  if (ymin2 < ymin) ymin = ymin2
                  ymin2 = max(0.0,real(safe_log10(minval(s% D_DSI(grid_min:grid_max))),kind=kind(ymin2)))
                  if (ymin2 < ymin) ymin = ymin2
                  ymin2 = max(0.0,real(safe_log10(minval(s% D_SH(grid_min:grid_max))),kind=kind(ymin2)))
                  if (ymin2 < ymin) ymin = ymin2
                  ymin2 = max(0.0,real(safe_log10(minval(s% D_SSI(grid_min:grid_max))),kind=kind(ymin2)))
                  if (ymin2 < ymin) ymin = ymin2
                  ymin2 = max(0.0,real(safe_log10(minval(s% D_ES(grid_min:grid_max))),kind=kind(ymin2)))
                  if (ymin2 < ymin) ymin = ymin2
                  ymin2 = max(0.0,real(safe_log10(minval(s% D_GSF(grid_min:grid_max))),kind=kind(ymin2)))
                  if (ymin2 < ymin) ymin = ymin2
                  ymin2 = max(0.0,real(safe_log10(minval(s% D_ST(grid_min:grid_max))),kind=kind(ymin2)))
                  if (ymin2 < ymin) ymin = ymin2
                  if (s% D_omega_flag) then
                     ymin2 = max(0.0,real(safe_log10(minval(s% D_omega(grid_min:grid_max))),kind=kind(ymin2)))
                     if (ymin2 < ymin) ymin = ymin2
                  end if
               end if
            end if

            dy = ymax-ymin
            if (MixDs_dymin /= -101) dy = MixDs_dymin
            ymax = ymax + 0.05*dy
            ymin = ymin - 0.05*dy

            call pgswin(xleft, xright, ymin, ymax)

            call pgsci(1)
            if (xaxis_numeric_labels_flag) then
               call show_box_pgstar(s,'BCNST','BCNSTV')
            else
               call show_box_pgstar(s,'BCST','BCNSTV')
            end if
            call show_left_yaxis_label_pgstar(s,'log D (cm\u2\d s\u-1\d)')

            call pgsch(txt_scale*s% Mixing_legend_txt_scale_factor)
            
            if (rotation) then
            
               if (s% D_omega_flag) then         
                  
                  do k=grid_min, grid_max
                     yvec(k) = safe_log10(s% D_mix_rotation(k))
                  end do
                  call pgsci(clr_rotation)
                  call pgslw(lw)
                  call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
                  call pgslw(lw_sav)               
                  
               else if (.not. s% Mixing_show_rotation_details) then
               
                  do k=grid_min, grid_max
                     yvec(k) = safe_log10( &
                        s% D_DSI_factor  * s% D_DSI(k)  + &
                        s% D_SH_factor   * s% D_SH(k)   + &
                        s% D_SSI_factor  * s% D_SSI(k)  + &
                        s% D_ES_factor   * s% D_ES(k)   + &
                        s% D_GSF_factor  * s% D_GSF(k)  + &
                        s% D_ST_factor   * s% D_ST(k))
                  end do
                  call pgsci(clr_rotation)
                  call pgslw(lw)
                  call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
                  call pgslw(lw_sav)
                  
               end if
            
               if (s% Mixing_show_rotation_details) then

                  D_DSI_factor = s% D_DSI_factor
                  D_SH_factor = s% D_SH_factor
                  D_SSI_factor = s% D_SSI_factor
                  D_ES_factor = s% D_ES_factor
                  D_GSF_factor = s% D_GSF_factor
                  D_ST_factor = s% D_ST_factor
                  D_visc_factor = s% D_visc_factor

                  do k=grid_min, grid_max
                     yvec(k) = safe_log10(D_ST_factor*s% D_ST(k))
                  end do
                  call pgsci(clr_IndianRed)
                  call pgslw(lw)
                  call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
                  call pgslw(lw_sav)

                  do k=grid_min, grid_max
                     yvec(k) = safe_log10(D_visc_factor*s% D_visc(k))
                  end do
                  call pgsci(clr_Silver)
                  call pgslw(lw)
                  call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
                  call pgslw(lw_sav)

                  do k=grid_min, grid_max
                     yvec(k) = safe_log10(D_DSI_factor*s% D_DSI(k))
                  end do
                  call pgsci(clr_Goldenrod)
                  call pgslw(lw)
                  call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
                  call pgslw(lw_sav)

                  do k=grid_min, grid_max
                     yvec(k) = safe_log10(D_SH_factor*s% D_SH(k))
                  end do
                  call pgsci(clr_Lilac)
                  call pgslw(lw)
                  call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
                  call pgslw(lw_sav)

                  do k=grid_min, grid_max
                     yvec(k) = safe_log10(D_SSI_factor*s% D_SSI(k))
                  end do
                  call pgsci(clr_Coral)
                  call pgslw(lw)
                  call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
                  call pgslw(lw_sav)

                  do k=grid_min, grid_max
                     yvec(k) = safe_log10(D_ES_factor*s% D_ES(k))
                  end do
                  call pgsci(clr_FireBrick)
                  !call show_right_yaxis_label_pgmtxt_pgstar(s,0.6,0.5,'ES',1.5)
                  call pgslw(lw)
                  call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
                  call pgslw(lw_sav)

                  do k=grid_min, grid_max
                     yvec(k) = safe_log10(D_GSF_factor*s% D_GSF(k))
                  end do
                  call pgsci(clr_MediumSlateBlue)
                  call pgslw(lw)
                  call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
                  call pgslw(lw_sav)
               
               end if
            
            end if
            
            call pgsci(clr_convection)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), y_conv(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgsci(clr_leftover_convection)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), y_left(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgsci(clr_semiconvection)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), y_sc(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgsci(clr_thermohaline)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), y_th(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgsci(clr_overshoot)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), y_ovr(grid_min:grid_max))
            call pgslw(lw_sav)

            if (s% RTI_flag) then
               call pgsci(clr_rayleigh_taylor)
               call pgslw(lw)
               call pgline(npts, xvec(grid_min:grid_max), y_RTI_mix(grid_min:grid_max))
               call pgslw(lw_sav)
            end if
            
            call pgsci(clr_minimum)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), y_min_mix(grid_min:grid_max))
            call pgslw(lw_sav)
            
            ! now do legend lines

            call pgsave
            call pgsvp(legend_xmin, legend_xmax, legend_ymin, legend_ymax)
            call pgswin(0.0, 1.0, ymin, ymax)
            number_of_legend_lines = 6
            if (rotation .and. s% Mixing_show_rotation_details) &
                  number_of_legend_lines = number_of_legend_lines + 7
            cnt = 1
            cnt = mixing_line_legend(cnt, clr_convection, &
               lw, lw_sav, txt_scale, 'convection')
            cnt = mixing_line_legend(cnt, clr_overshoot, &
               lw, lw_sav, txt_scale, 'overshooting')
            if (s% alpha_semiconvection > 0) then
               cnt = mixing_line_legend(cnt, clr_semiconvection, &
                  lw, lw_sav, txt_scale, 'semiconvection')
            end if
            if (s% thermohaline_coeff > 0) then
               cnt = mixing_line_legend(cnt, clr_thermohaline, &
                  lw, lw_sav, txt_scale, 'thermohaline')
            end if
            if (rotation) then
               if (s% D_omega_flag) then         
                  cnt = mixing_line_legend(cnt, clr_rotation, &
                     lw, lw_sav, txt_scale, 'D_mix_rotation')
               else
                  cnt = mixing_line_legend(cnt, clr_rotation, &
                     lw, lw_sav, txt_scale, 'rotation')
               end if
            end if
            if (s% RTI_flag) then
               cnt = mixing_line_legend(cnt, clr_rayleigh_taylor, &
                  lw, lw_sav, txt_scale, 'RTI')
            end if
            if (rotation .and. s% Mixing_show_rotation_details) then
               cnt = mixing_line_legend(cnt, clr_IndianRed, &
                  lw, lw_sav, txt_scale, 'ST')
               cnt = mixing_line_legend(cnt, clr_Silver, &
                  lw, lw_sav, txt_scale, 'visc')
               cnt = mixing_line_legend(cnt, clr_Goldenrod, &
                  lw, lw_sav, txt_scale, 'DSI')
               cnt = mixing_line_legend(cnt, clr_Lilac, &
                  lw, lw_sav, txt_scale, 'SH')
               cnt = mixing_line_legend(cnt, clr_Coral, &
                  lw, lw_sav, txt_scale, 'SSI')
               cnt = mixing_line_legend(cnt, clr_FireBrick, &
                  lw, lw_sav, txt_scale, 'ES')
               cnt = mixing_line_legend(cnt, clr_MediumSlateBlue, &
                  lw, lw_sav, txt_scale, 'GSF')
            end if
            
            call pgunsa

            call show_pgstar_decorator(s%id,s% mixing_use_decorator,s% mixing_pgstar_decorator, 0, ierr)

         end subroutine plot


         integer function mixing_line_legend( &
               cnt, clr, lw, lw_sav, txt_scale, str)
            integer, intent(in) :: cnt, clr, lw, lw_sav
            real, intent(in) :: txt_scale
            character (len=*), intent(in) :: str
            real :: dx, dyline, ypos, xpts(2), ypts(2)
            call pgsci(clr)
            dx = 0.1
            dyline = (ymax-ymin)/number_of_legend_lines
            ypos = ymax - cnt*dyline
            xpts(1) = 1.3*dx
            xpts(2) = xpts(1) + 2.3*dx
            ypts = ypos + dyline*0.1
            call pgslw(lw)
            call pgline(2, xpts, ypts)
            call pgslw(lw_sav)
            call pgsci(1)
            call pgsch(txt_scale*s% Mixing_legend_txt_scale_factor)
            call pgptxt(xpts(2) + dx, ypos, 0.0, 0.0, trim(str))
            mixing_line_legend = cnt + 1
         end function mixing_line_legend


      end subroutine MixDs_plot


      end module pgstar_mixing_Ds

