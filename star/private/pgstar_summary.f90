! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module pgstar_summary

      use star_private_def
      use const_def
      use pgstar_support
      use star_pgstar

      implicit none


      contains


      subroutine Text_Summary_plot(id, device_id, array_ix, ierr)
         integer, intent(in) :: id, device_id, array_ix
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         if (.not. present(array_ix)) then
            ierr = -1
            return
         end if

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call pgslct(device_id)
         call pgbbuf()
         call pgeras()
         call do_Text_Summary_plot(s, id, device_id, array_ix, &
            s% pg% Text_Summary_xleft(array_ix), s% pg% Text_Summary_xright(array_ix), &
            s% pg% Text_Summary_ybot(array_ix), s% pg% Text_Summary_ytop(array_ix), .false., &
            s% pg% Text_Summary_title(array_ix), s% pg% Text_Summary_txt_scale(array_ix), s% pg% Text_Summary_dxval(array_ix), ierr)
         if (ierr /= 0) return
         call pgebuf()
      end subroutine Text_Summary_plot


      subroutine do_Text_Summary_plot(s, id, device_id, array_ix, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, dxval, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id, array_ix
         real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale, dxval
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         integer, intent(out) :: ierr

         call Summary_plot(s, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, dxval, &
            s% pg% Text_Summary_num_rows(array_ix), s% pg% Text_Summary_num_cols(array_ix), &
            s% pg% Text_Summary_name(array_ix, :, :), ierr)
      end subroutine do_Text_Summary_plot


      subroutine Summary_plot(s, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, dxval, &
            Text_Summary_num_rows, Text_Summary_num_cols, &
            Text_Summary_name, ierr)

         use utils_lib
         use chem_def
         use net_def

         type (star_info), pointer :: s
         integer, intent(in) :: device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale, dxval
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         integer, intent(in) :: Text_Summary_num_rows, Text_Summary_num_cols
         character (len=*), intent(in) :: Text_Summary_name(:, :)
         integer, intent(out) :: ierr

         integer :: col, num_cols, num_rows

         include 'formats'

         ierr = 0

         num_rows = Text_Summary_num_rows
         num_cols = Text_Summary_num_cols
         if (num_rows <= 0 .or. num_cols <= 0) return

         call pgsave
         call pgsch(txt_scale)

         call pgsvp(winxmin, winxmax, winymin, winymax)
         call pgsci(1)
         call pgscf(1)
         call pgswin(0.0,1.0,0.0,1.0)
         call show_title_pgstar(s, title)
         call pgsch(txt_scale*0.8)

         do col = 1, num_cols
            call show_column(col, num_rows)
         end do

         call pgunsa


         contains


         subroutine show_column(col, num_rows)
            use history, only: get_history_specs, get_history_values, get1_hist_value
            integer, intent(in) :: col, num_rows

            real(dp) :: values(num_rows)
            integer :: int_values(num_rows), specs(num_rows), int_val
            logical :: is_int_value(num_rows)
            logical :: failed_to_find_value(num_rows)

            integer :: i, cnt
            real :: xpos0, dxpos, ypos, dypos
            real(dp) :: val

            call get_history_specs(s, num_rows, Text_Summary_name(:,col), specs, .false.)
            call get_history_values( &
               s, num_rows, specs, &
               is_int_value, int_values, values, failed_to_find_value)

            xpos0 = (real(col) - 0.5)/real(num_cols)

            dxpos = 0.00

            ypos = 0.90
            dypos = -0.95/num_rows

            do i=1,num_rows
               if (i > 1) ypos = ypos + dypos
               if (failed_to_find_value(i)) then
                  if (.not. get1_hist_value(s, Text_Summary_name(i,col), val)) then
                     if (len_trim(Text_Summary_name(i,col)) > 0) &
                        write(*,'(a)') 'failed_to_find_value ' // trim(Text_Summary_name(i,col)) &
                           // '. check that it is in your history_columns.list'
                     cycle
                  end if
                  int_val = int(val)
                  if (abs(val - dble(int_val)) < 1d-10*max(1d-10,abs(val))) then
                     cnt = write_info_line_int(0, ypos, xpos0, dxpos, dxval, &
                           Text_Summary_name(i,col), int_val)
                     cycle
                  else
                     values(i) = val
                  end if
               else if (is_int_value(i)) then
                  cnt = write_info_line_int(0, ypos, xpos0, dxpos, dxval, &
                        Text_Summary_name(i,col), int_values(i))
                  cycle
               end if
               if (values(i) == 0d0) then
                  cnt = write_info_line_int(0, ypos, xpos0, dxpos, dxval, &
                        Text_Summary_name(i,col), 0)
               else if (abs(values(i)) > 1d-3 .and. abs(values(i)) < 1d3) then
                  cnt = write_info_line_flt(0, ypos, xpos0, dxpos, dxval, &
                        Text_Summary_name(i,col), values(i))
               else
                  cnt = write_info_line_exp(0, ypos, xpos0, dxpos, dxval, &
                        Text_Summary_name(i,col), values(i))
               end if
            end do

         end subroutine show_column


      end subroutine Summary_plot


      end module pgstar_summary

