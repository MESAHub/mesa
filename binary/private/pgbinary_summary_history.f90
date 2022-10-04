! ***********************************************************************
!
!   Copyright (C) 2013-2022  Bill Paxton, Matthias Fabry
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

module pgbinary_summary_history

   use binary_private_def
   use const_def
   use pgbinary_support

   implicit none


contains


   subroutine summary_history_plot(id, device_id, ierr)
      implicit none
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr

      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return

      call pgslct(device_id)
      call pgbbuf()
      call pgeras()

      call do_summary_history_plot(b, id, device_id, &
         b% pg% Summary_History_xleft, b% pg% Summary_History_xright, &
         b% pg% Summary_History_ybot, b% pg% Summary_History_ytop, .false., &
         b% pg% Summary_History_title, b% pg% Summary_History_txt_scale, ierr)

      call pgebuf()

   end subroutine summary_history_plot


   subroutine do_summary_history_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)

      use utils_lib
      use chem_def
      use net_def
      use const_def, only : Msun, Rsun

      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      real, intent(in) :: txt_scale
      integer, intent(out) :: ierr

      character (len = strlen) :: yname
      real, pointer, dimension(:) :: xvec, yvec
      real :: xmin, xmax, windy, ymin, ymax, xmargin, &
         legend_xmin, legend_xmax, legend_ymin, legend_ymax
      integer :: lw, lw_sav, num_lines, &
         npts, step_min, step_max
      integer, parameter :: num_colors = 20
      integer :: colors(num_colors)

      include 'formats'

      ierr = 0

      step_min = b% pg% Summary_History_xmin
      if (step_min <= 0) step_min = 1
      step_max = b% pg% Summary_History_xmax
      if (step_max <= 0) step_max = b% model_number

      if (step_min >= b% model_number) step_min = 1

      if (b% pg% Summary_History_max_width > 0) &
         step_min = max(step_min, step_max - b% pg% Summary_History_max_width)

      npts = count_hist_points(b, step_min, step_max)
      if (npts <= 1) return

      xmin = real(max(1, step_min))
      xmax = real(min(b% model_number, step_max))

      num_lines = b% pg% Summary_History_num_lines

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

      allocate(xvec(npts), yvec(npts))

      call set_hist_points_steps(&
         b, step_min, step_max, npts, xvec, ierr)
      if (ierr /= 0) then
         write(*, *) 'set_hist_points_steps failed for PGSTAR Summary History'
         return
      end if

      if (ierr == 0) then
         call pgsave
         call pgsch(txt_scale)
         call plot(ierr)
         call pgunsa
      end if

      deallocate(xvec, yvec)


   contains


      subroutine plot(ierr)
         use rates_def
         integer, intent(out) :: ierr

         integer :: j, ii, jj, i, cnt, k, yaxis_id
         logical :: show(num_lines)
         logical, parameter :: dbg = .false.
         real :: ybot, yvec_min, yvec_max

         include 'formats'

         ymax = 1.02
         ymin = 0.0

         lw = b% pg% pgbinary_lw
         call pgqlw(lw_sav)

         call pgsvp(winxmin, winxmax, winymin, winymax)
         if (.not. subplot) then
            call show_model_number_pgbinary(b)
            call show_age_pgbinary(b)
         end if
         call show_title_pgbinary(b, title)

         ybot = 0
         call pgswin(xmin, xmax, ymin + ybot, ymax)
         call pgscf(1)
         call pgsci(1)
         call show_box_pgbinary(b, 'BCNST', 'BCNSTV')
         call show_left_yaxis_label_pgbinary(b, 'rel=(val-min)/(max-min)')

         cnt = 0
         do j = 1, num_lines

            yname = b% pg% Summary_History_name(j)
            if (len_trim(yname) == 0) then
               show(j) = .false.
               cycle
            end if

            show(j) = get1_yvec(yname, yvec)
            if (.not. show(j)) then
               write(*, *) 'failed to find history information for ' // trim(yname)
               cycle
            end if

            if (b% pg% Summary_History_scaled_value(j)) then ! scale yvec

               yvec_max = maxval(yvec(1:npts))
               yvec_min = minval(yvec(1:npts))
               show(j) = (yvec_max > yvec_min)
               if (.not. show(j)) then
                  write(*, 1) trim(yname) // ' same min max', yvec_max
                  cycle
               end if
               !write(*,1) 'relative ' // trim(yname), yvec_min, yvec_max
               do k = 1, npts
                  yvec(k) = (yvec(k) - yvec_min) / (yvec_max - yvec_min)
               end do

            else

               show(j) = .true.
               !write(*,1) 'absolute ' // trim(yname), yvec_min, yvec_max

            end if

            call pgslw(lw)
            cnt = summary_history_line(cnt, yvec)
            call pgslw(lw_sav)

         end do

         call pgsci(1)
         call show_xaxis_label_pgbinary(b, 'model number')

         ! show the legend
         call pgsave
         call pgsvp(legend_xmin, legend_xmax, legend_ymin, legend_ymax)
         call pgswin(0.0, 1.0, ymin, ymax)
         cnt = 0
         do j = 1, num_lines
            if (.not. show(j)) cycle
            if (len_trim(b% pg% Summary_History_legend(j)) == 0) then
               cnt = summary_history_line_legend(&
                  cnt, b% pg% Summary_History_name(j))
            else
               cnt = summary_history_line_legend(&
                  cnt, b% pg% Summary_History_legend(j))
            end if
         end do
         call pgunsa

         call show_pgbinary_decorator(b% binary_id, b% pg% Summary_history_use_decorator, &
            b% pg% Summary_history_pgbinary_decorator, 0, ierr)

      end subroutine plot


      logical function get1_yvec(name, vec)
         character (len = *) :: name
         real, dimension(:), pointer :: vec
         get1_yvec = get1_hist_yvec(b, step_min, step_max, npts, name, vec)
      end function get1_yvec


      integer function summary_history_line(cnt, yvec)
         integer, intent(in) :: cnt
         real, intent(in) :: yvec(:)
         integer :: iclr
         iclr = cnt - num_colors * (cnt / num_colors) + 1
         summary_history_line = cnt + 1
         call pgsci(colors(iclr))
         call pgline(npts, xvec(1:npts), yvec(1:npts))
      end function summary_history_line


      integer function summary_history_line_legend(cnt, name)
         integer, intent(in) :: cnt
         character (len = *), intent(in) :: name
         real :: dx, dyline, ypos, xpts(2), ypts(2)
         integer :: iclr, num_max
         num_max = max_num_Summary_History_Lines
         summary_history_line_legend = cnt
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
         summary_history_line_legend = cnt + 1
      end function summary_history_line_legend


   end subroutine do_summary_history_plot


end module pgbinary_summary_history

