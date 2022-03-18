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

      module pgstar_summary_profile

      use star_private_def
      use const_def
      use pgstar_support
      use star_pgstar

      implicit none


      contains


      subroutine summary_profile_plot(id, device_id, ierr)
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

         call do_summary_profile_plot(s, id, device_id, &
            s% pg% Summary_Profile_xleft, s% pg% Summary_Profile_xright, &
            s% pg% Summary_Profile_ybot, s% pg% Summary_Profile_ytop, .false., &
            s% pg% Summary_Profile_title, s% pg% Summary_Profile_txt_scale, ierr)

         call pgebuf()

      end subroutine summary_profile_plot


      subroutine do_summary_profile_plot(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         integer, intent(out) :: ierr
         call do_summary_profile_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, s% pg% Summary_Profile_xaxis_name, &
            s% pg% Summary_Profile_xmin, s% pg% Summary_Profile_xmax, &
            s% pg% Summary_Profile_xaxis_reversed, &
            .false., .true., ierr)
      end subroutine do_summary_profile_plot


      subroutine do_summary_profile_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            xaxis_name, xaxis_min, xaxis_max, xaxis_reversed, &
            panel_flag, xaxis_numeric_labels_flag, ierr)
         use utils_lib
         use chem_def
         use net_def
         use const_def, only: Msun, Rsun

         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: &
            winxmin, winxmax, winymin, winymax, xaxis_min, xaxis_max
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title, xaxis_name
         real, intent(in) :: txt_scale
         logical, intent(in) :: &
            xaxis_reversed, panel_flag, xaxis_numeric_labels_flag
         integer, intent(out) :: ierr

         character (len=strlen) :: yname
         real, allocatable, dimension(:) :: xvec, yvec, unshifted_xvec
         real :: xmin, xmax, xleft, xright, dx, windy, &
            ymin, ymax, xmargin, &
            legend_xmin, legend_xmax, legend_ymin, legend_ymax
         integer :: lw, lw_sav, grid_min, grid_max, npts, nz, num_lines
         integer, parameter :: num_colors = 20
         integer :: colors(num_colors)

         include 'formats'

         ierr = 0

         nz = s% nz

         num_lines = s% pg% Summary_Profile_num_lines

         colors(:) = (/ &
               clr_MediumSlateBlue, clr_Goldenrod, clr_LightSkyBlue, clr_Lilac, &
               clr_Coral, clr_Crimson, clr_LightSkyGreen, clr_DarkGray, &
               clr_Tan, clr_IndianRed, clr_Gold, &
               clr_Teal, clr_Silver, clr_BrightBlue, clr_FireBrick, &
               clr_RoyalPurple, clr_SlateGray, clr_LightSteelBlue, &
               clr_Gray, clr_RoyalBlue /)

         windy = winymax - winymin

         legend_xmin = winxmax - 0.01
         legend_xmax = 0.99
         legend_ymin = winymin
         legend_ymax = winymax

         allocate(xvec(nz), yvec(nz),unshifted_xvec(nz))

         xmargin = 0
         call set_xaxis_bounds( &
            s, xaxis_name, xaxis_min, xaxis_max, xaxis_reversed, xmargin, &
            xvec, xmin, xmax, xleft, xright, dx, &
            grid_min, grid_max, npts, ierr)

         if (ierr == 0) then
            call pgsave
            call pgsch(txt_scale)
            call plot(ierr)
            call pgunsa
         end if

         deallocate(xvec, yvec,unshifted_xvec)


         contains


         subroutine plot(ierr)
            use rates_def
            use profile_getval, only : get_profile_val,get_profile_id
            integer, intent(out) :: ierr

            integer :: j, ii, jj, i, cnt, k, yaxis_id
            logical :: show(num_lines)
            logical, parameter :: dbg = .false.
            real :: ybot, yvec_min, yvec_max

            include 'formats'

            ymax = 1.02
            ymin = 0.0

            lw = s% pg% pgstar_lw
            call pgqlw(lw_sav)

            call pgsvp(winxmin, winxmax, winymin, winymax)
            if (.not. panel_flag) then
               if (.not. subplot) then
                  call show_model_number_pgstar(s)
                  call show_age_pgstar(s)
               end if
               call show_title_pgstar(s, title)
            end if

            ybot = -0.02
            call pgswin(xleft, xright, ymin+ybot, ymax)
            call pgscf(1)
            call pgsci(1)
            if (xaxis_numeric_labels_flag) then
               call show_box_pgstar(s,'BCNST','BCNSTV')
            else
               call show_box_pgstar(s,'BCST','BCNSTV')
            end if

            do k=1,nz
               unshifted_xvec(k) = xvec(k)
            end do
            if (grid_min > 1) then
               do k=1,npts
                  xvec(k) = xvec(k+grid_min-1)
               end do
            end if

            cnt = 0
            do j = 1, num_lines

               yname = s% pg% Summary_Profile_name(j)
               if (len_trim(yname) == 0 .or. trim(yname) == trim(xaxis_name)) then
                  show(j) = .false.
                  cycle
               end if

               yaxis_id = get_profile_id(s, yname)
               if (yaxis_id <= 0) then
                  write(*,*) &
                     'bad yaxis for Profile panels plot ' // trim(yname)
                  return
               end if

               do k=1,npts
                  yvec(k) = get_profile_val(s, yaxis_id, k+grid_min-1)
               end do

               if (s% pg% Summary_Profile_scaled_value(j)) then ! scale yvec

                  yvec_max = maxval(yvec(1:npts))
                  yvec_min = minval(yvec(1:npts))
                  show(j) = (yvec_max > yvec_min)
                  if (.not. show(j)) then
                     cycle
                  end if
                  do k=1,npts
                     yvec(k) = (yvec(k) - yvec_min)/(yvec_max - yvec_min)
                  end do

               else

                  show(j) = .true.

               end if

               call pgslw(lw)
               cnt = summary_profile_line(cnt, yvec)
               call pgslw(lw_sav)

            end do

            if (.not. panel_flag) then ! show xaxis info
               call pgsci(1)
               call show_xaxis_name(s,xaxis_name,ierr)
               if (ierr == 0) then ! show mix regions at bottom of plot
                  call pgslw(10)
                  call show_mix_regions_on_xaxis( &
                     s,ymin+ybot,ymax,grid_min,grid_max,unshifted_xvec)
               end if
            end if

            ! show the legend
            call pgsave
            call pgsvp(legend_xmin, legend_xmax, legend_ymin, legend_ymax)
            call pgswin(0.0, 1.0, ymin, ymax)
            cnt = 0
            do j=1,num_lines
               if (.not. show(j)) cycle
               if (len_trim(s% pg% Summary_Profile_legend(j)) == 0) then
                  cnt = summary_profile_line_legend( &
                           cnt,s% pg% Summary_Profile_name(j))
               else
                  cnt = summary_profile_line_legend( &
                           cnt,s% pg% Summary_Profile_legend(j))
               end if
            end do
            call pgunsa

         call show_pgstar_decorator(s%id, s% pg% summary_profile_use_decorator, &
               s% pg% summary_profile_pgstar_decorator, 0, ierr)


         end subroutine plot


         integer function summary_profile_line(cnt, yvec)
            integer, intent(in) :: cnt
            real, intent(in) :: yvec(:)
            integer :: iclr
            iclr = cnt - num_colors*(cnt/num_colors) + 1
            summary_profile_line = cnt + 1
            call pgsci(colors(iclr))
            call pgline(npts, xvec, yvec)
         end function summary_profile_line


         integer function summary_profile_line_legend(cnt, name)
            integer, intent(in) :: cnt
            character (len=*), intent(in) :: name
            real :: dx, dyline, ypos, xpts(2), ypts(2)
            integer :: iclr, num_max
            num_max = max_num_Summary_Profile_Lines
            summary_profile_line_legend = cnt
            iclr = cnt - num_colors*(cnt/num_colors) + 1
            call pgsci(colors(iclr))
            dx = 0.1
            dyline = (ymax-ymin)/num_max
            ypos = ymax - (cnt+1.5)*dyline
            xpts(1) = 1.3*dx
            xpts(2) = xpts(1) + 2.3*dx
            ypts = ypos + dyline*0.1
            call pgslw(lw)
            call pgline(2, xpts, ypts)
            call pgslw(lw_sav)
            call pgsci(1)
            call pgsch(txt_scale*0.70)
            call pgptxt(xpts(2) + dx, ypos, 0.0, 0.0, name)
            summary_profile_line_legend = cnt + 1
         end function summary_profile_line_legend


      end subroutine do_summary_profile_panel


      end module pgstar_summary_profile

