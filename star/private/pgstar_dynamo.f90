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

      module pgstar_dynamo

      use star_private_def
      use const_def
      use pgstar_support
      use star_pgstar

      implicit none


      contains


      subroutine Dynamo_plot(id, device_id, ierr)
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

         call do_Dynamo_plot(s, id, device_id, &
            s% pg% Dynamo_xleft, s% pg% Dynamo_xright, &
            s% pg% Dynamo_ybot, s% pg% Dynamo_ytop, .false., &
            s% pg% Dynamo_title, s% pg% Dynamo_txt_scale, ierr)

         call pgebuf()

      end subroutine Dynamo_plot


      subroutine do_Dynamo_plot(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         integer, intent(out) :: ierr
         call do_Dynamo_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, s% pg% Dynamo_xaxis_name, &
            s% pg% Dynamo_xmin, s% pg% Dynamo_xmax, &
            s% pg% Dynamo_xaxis_reversed, .false., .true., ierr)
      end subroutine do_Dynamo_plot


      subroutine do_Dynamo_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            xaxis_name, xmin, xmax, reverse_xaxis, &
            panel_flag, xaxis_numeric_labels_flag, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: &
            winxmin, winxmax, winymin, winymax, xmin, xmax
         character (len=*), intent(in) :: title, xaxis_name
         real, intent(in) :: txt_scale
         logical, intent(in) :: subplot, &
            reverse_xaxis, panel_flag, xaxis_numeric_labels_flag
         integer, intent(out) :: ierr
         call Dyn_plot(s, device_id, &
            s% pg% show_Dynamo_annotation1, s% pg% show_Dynamo_annotation2, &
            s% pg% show_Dynamo_annotation3, &
            xaxis_name, xmin, xmax, reverse_xaxis, &
            s% pg% Dynamo_ymin_left, s% pg% Dynamo_ymax_left, s% pg% Dynamo_dymin_left, &
            s% pg% Dynamo_ymin_right, s% pg% Dynamo_ymax_right, s% pg% Dynamo_dymin_right, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            panel_flag, xaxis_numeric_labels_flag, ierr)
      end subroutine do_Dynamo_panel


      subroutine Dyn_plot(s, device_id, &
            show_Dyn_annotation1, show_Dyn_annotation2, show_Dyn_annotation3, &
            Dyn_xaxis_name, Dyn_xmin, Dyn_xmax, Dyn_reverse_xaxis, &
            Dyn_ymin_left, Dyn_ymax_left, Dyn_dymin_left, &
            Dyn_ymin_right, Dyn_ymax_right, Dyn_dymin_right, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            panel_flag, xaxis_numeric_labels_flag, ierr)

         use utils_lib
         implicit none

         type (star_info), pointer :: s
         integer, intent(in) :: device_id
         logical, intent(in) :: subplot, &
            show_Dyn_annotation1, show_Dyn_annotation2, show_Dyn_annotation3
         character (len=*), intent(in) :: Dyn_xaxis_name, title
         real, intent(in) :: &
            Dyn_xmin, Dyn_xmax, &
            Dyn_ymin_left, Dyn_ymax_left, Dyn_dymin_left, &
            Dyn_ymin_right, Dyn_ymax_right, Dyn_dymin_right
         real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
         logical, intent(in) :: &
            Dyn_reverse_xaxis, panel_flag, xaxis_numeric_labels_flag
         integer, intent(out) :: ierr

         real :: windy, xmargin
         real :: xmin, xmax, xleft, xright, dx, tmp, ymin, ymax, ymin2, ymax2, dy
         integer :: grid_min, grid_max, npts, nz
         real, allocatable, dimension(:) :: xvec, yvec, yvec2, yvec3

         include 'formats'
         ierr = 0

         if (.not. s% rotation_flag) return

         xmargin = 0

         nz = s% nz
         allocate (xvec(nz), yvec(nz), yvec2(nz), yvec3(nz))

         call set_xaxis_bounds( &
            s, Dyn_xaxis_name, Dyn_xmin, Dyn_xmax, Dyn_reverse_xaxis, xmargin, &
            xvec, xmin, xmax, xleft, xright, dx, &
            grid_min, grid_max, npts, ierr)
         if (ierr /= 0) return

         call pgsave
         call pgsch(txt_scale)

         call plot(ierr)
         if (ierr == 0) call show_annotations(s, &
            show_Dyn_annotation1, show_Dyn_annotation2, show_Dyn_annotation3)

         call pgunsa

         deallocate(xvec, yvec, yvec2, yvec3)
         
         contains


         subroutine plot(ierr)
            use pgstar_support, only: show_convective_section, show_semiconvective_section, &
               show_thermohaline_section, show_overshoot_section
            integer, intent(out) :: ierr

            integer :: lw, lw_sav, k
            real :: ybot, eps, &
               default_ymax_left, default_ymin_left, &
               default_ymax_right, default_ymin_right
            character (len=128) :: str

            include 'formats'
            ierr = 0

            lw = s% pg% pgstar_lw
            call pgqlw(lw_sav)

            if (.not. panel_flag) then
               call pgsvp(winxmin, winxmax, winymin, winymax)
               call show_title_pgstar(s, title)
               if (.not. subplot) then
                  call show_model_number_pgstar(s)
                  call show_age_pgstar(s)
               end if
               call pgsci(1)
               call show_xaxis_name(s,Dyn_xaxis_name,ierr)
               if (ierr /= 0) return
            end if

            default_ymax_left = 10
            default_ymin_left = -2

            do k=1,nz
               yvec(k) = safe_log10(s% dynamo_B_phi(k))
               yvec2(k) = safe_log10(s% dynamo_B_r(k))
            end do

            if (Dyn_ymax_left /= -101) then
               ymax = Dyn_ymax_left
            else
               ymax = max(default_ymax_left,maxval(yvec(grid_min:grid_max)))
               ymax2 = max(default_ymax_left,maxval(yvec2(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
            end if

            if (Dyn_ymin_left /= -101) then
               ymin = Dyn_ymin_left
            else
               ymin = max(default_ymin_left,minval(yvec(grid_min:grid_max)))
               ymin2 = max(default_ymin_left,minval(yvec2(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
            end if

            dy = ymax-ymin
            if (Dyn_dymin_left /= -101) dy = Dyn_dymin_left

            ymax = ymax + 0.1*dy
            ymin = ymin - 0.1*dy

            call pgswin(xleft, xright, ymin, ymax)
            call pgscf(1)
            call pgsci(1)
            call show_box_pgstar(s,'','BNSTV')

            call pgsci(clr_Teal)
            call pgsch(txt_scale*s% pg% Dynamo_legend_txt_scale_factor)
            call show_left_yaxis_label_pgstar(s,'log B\dphi\u (Gauss)',-0.5)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgsci(clr_Coral)
            call pgsch(txt_scale*s% pg% Dynamo_legend_txt_scale_factor)
            call show_left_yaxis_label_pgstar(s,'log B\dr\u (Gauss)',1.3)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec2(grid_min:grid_max))
            call pgslw(lw_sav)


            ! right axis

            lw = 8

            default_ymax_right = 0
            default_ymin_right = -10

            ! log omega
            do k=1,nz
               yvec(k) = safe_log10(s% omega(k))
               yvec2(k) = safe_log10(s% j_rot(k)) - 20d0
            end do

            if (Dyn_ymax_right /= -101) then
               ymax = Dyn_ymax_right
            else
               ymax = max(default_ymax_right,maxval(yvec(grid_min:grid_max)))
               ymax2 = max(default_ymax_right,maxval(yvec2(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
            end if
            if (Dyn_ymin_right /= -101) then
               ymin = Dyn_ymin_right
            else
               ymin = max(default_ymin_right,minval(yvec(grid_min:grid_max)))
               ymin2 = max(default_ymin_right,minval(yvec2(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
            end if

            dy = ymax-ymin
            if (Dyn_dymin_right /= -101) dy = Dyn_dymin_right

            ymax = ymax + 0.1*dy
            ymin = ymin - 0.1*dy

            call pgswin(xleft, xright, ymin, ymax)

            call pgscf(1)
            call pgsci(1)
            if (xaxis_numeric_labels_flag) then
               call show_box_pgstar(s,'BCNST','CMSTV')
            else
               call show_box_pgstar(s,'BCST','CMSTV')
            end if

            call pgsci(clr_FireBrick)
            call pgsch(txt_scale*s% pg% Dynamo_legend_txt_scale_factor)
            call show_right_yaxis_label_pgstar( &
                  s,'log j (10\u20\d cm\u2\d/s)',1.3)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec2(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgsci(clr_RoyalBlue)
            call pgsch(txt_scale*s% pg% Dynamo_legend_txt_scale_factor)
            call show_right_yaxis_label_pgstar(s,'log \(0650) (rad/s)',-0.5)
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
            call pgslw(lw_sav)


            if (.not. panel_flag) then ! show mix regions at bottom of plot
               ybot = -0.05
               call pgswin(xleft, xright, ybot, 0.85)
               call pgslw(10)
               call show_mix_regions_on_xaxis( &
                  s,ybot,0.85,grid_min,grid_max,xvec)
            end if

         call show_pgstar_decorator(s%id,s% pg% dynamo_use_decorator, &
            s% pg% dynamo_pgstar_decorator, 0, ierr)


         end subroutine plot


      end subroutine Dyn_plot


      end module pgstar_dynamo

