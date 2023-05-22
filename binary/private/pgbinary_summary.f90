! ***********************************************************************
!
!   Copyright (C) 2010-2022  The MESA Team, Bill Paxton & Matthias Fabry
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

module pgbinary_summary

   use binary_private_def
   use const_def
   use pgbinary_support

   implicit none


contains


   subroutine Text_Summary1_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Text_Summary_plot(b, id, device_id, &
         b% pg% Text_Summary_xleft(1), b% pg% Text_Summary_xright(1), &
         b% pg% Text_Summary_ybot(1), b% pg% Text_Summary_ytop(1), .false., &
         b% pg% Text_Summary_title(1), b% pg% Text_Summary_txt_scale(1), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Text_Summary1_plot


   subroutine do_Text_Summary1_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call Summary_plot(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
         b% pg% Text_Summary_num_rows(1), b% pg% Text_Summary_num_cols(1), &
         b% pg% Text_Summary_name(1, :, :), ierr)
   end subroutine do_Text_Summary1_plot


   subroutine Text_Summary2_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Text_Summary2_plot(b, id, device_id, &
         b% pg% Text_Summary_xleft(2), b% pg% Text_Summary_xright(2), &
         b% pg% Text_Summary_ybot(2), b% pg% Text_Summary_ytop(2), .false., &
         b% pg% Text_Summary_title(2), b% pg% Text_Summary_txt_scale(2), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Text_Summary2_plot


   subroutine do_Text_Summary2_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call Summary_plot(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
         b% pg% Text_Summary_num_rows(2), b% pg% Text_Summary_num_cols(2), &
         b% pg% Text_Summary_name(2, :, :), ierr)
   end subroutine do_Text_Summary2_plot


   subroutine Text_Summary3_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Text_Summary3_plot(b, id, device_id, &
         b% pg% Text_Summary_xleft(3), b% pg% Text_Summary_xright(3), &
         b% pg% Text_Summary_ybot(3), b% pg% Text_Summary_ytop(3), .false., &
         b% pg% Text_Summary_title(3), b% pg% Text_Summary_txt_scale(3), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Text_Summary3_plot


   subroutine do_Text_Summary3_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call Summary_plot(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
         b% pg% Text_Summary_num_rows(3), b% pg% Text_Summary_num_cols(3), &
         b% pg% Text_Summary_name(3, :, :), ierr)
   end subroutine do_Text_Summary3_plot


   subroutine Text_Summary4_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Text_Summary4_plot(b, id, device_id, &
         b% pg% Text_Summary_xleft(4), b% pg% Text_Summary_xright(4), &
         b% pg% Text_Summary_ybot(4), b% pg% Text_Summary_ytop(4), .false., &
         b% pg% Text_Summary_title(4), b% pg% Text_Summary_txt_scale(4), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Text_Summary4_plot


   subroutine do_Text_Summary4_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call Summary_plot(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
         b% pg% Text_Summary_num_rows(4), b% pg% Text_Summary_num_cols(4), &
         b% pg% Text_Summary_name(4, :, :), ierr)
   end subroutine do_Text_Summary4_plot


   subroutine Text_Summary5_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Text_Summary5_plot(b, id, device_id, &
         b% pg% Text_Summary_xleft(5), b% pg% Text_Summary_xright(5), &
         b% pg% Text_Summary_ybot(5), b% pg% Text_Summary_ytop(5), .false., &
         b% pg% Text_Summary_title(5), b% pg% Text_Summary_txt_scale(5), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Text_Summary5_plot


   subroutine do_Text_Summary5_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call Summary_plot(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
         b% pg% Text_Summary_num_rows(5), b% pg% Text_Summary_num_cols(5), &
         b% pg% Text_Summary_name(5, :, :), ierr)
   end subroutine do_Text_Summary5_plot


   subroutine Text_Summary6_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Text_Summary6_plot(b, id, device_id, &
         b% pg% Text_Summary_xleft(6), b% pg% Text_Summary_xright(6), &
         b% pg% Text_Summary_ybot(6), b% pg% Text_Summary_ytop(6), .false., &
         b% pg% Text_Summary_title(6), b% pg% Text_Summary_txt_scale(6), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Text_Summary6_plot


   subroutine do_Text_Summary6_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call Summary_plot(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
         b% pg% Text_Summary_num_rows(6), b% pg% Text_Summary_num_cols(6), &
         b% pg% Text_Summary_name(6, :, :), ierr)
   end subroutine do_Text_Summary6_plot


   subroutine Text_Summary7_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Text_Summary7_plot(b, id, device_id, &
         b% pg% Text_Summary_xleft(7), b% pg% Text_Summary_xright(7), &
         b% pg% Text_Summary_ybot(7), b% pg% Text_Summary_ytop(7), .false., &
         b% pg% Text_Summary_title(7), b% pg% Text_Summary_txt_scale(7), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Text_Summary7_plot


   subroutine do_Text_Summary7_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call Summary_plot(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
         b% pg% Text_Summary_num_rows(7), b% pg% Text_Summary_num_cols(7), &
         b% pg% Text_Summary_name(7, :, :), ierr)
   end subroutine do_Text_Summary7_plot


   subroutine Text_Summary8_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Text_Summary8_plot(b, id, device_id, &
         b% pg% Text_Summary_xleft(8), b% pg% Text_Summary_xright(8), &
         b% pg% Text_Summary_ybot(8), b% pg% Text_Summary_ytop(8), .false., &
         b% pg% Text_Summary_title(8), b% pg% Text_Summary_txt_scale(8), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Text_Summary8_plot


   subroutine do_Text_Summary8_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call Summary_plot(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
         b% pg% Text_Summary_num_rows(8), b% pg% Text_Summary_num_cols(8), &
         b% pg% Text_Summary_name(8, :, :), ierr)
   end subroutine do_Text_Summary8_plot


   subroutine Text_Summary9_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Text_Summary9_plot(b, id, device_id, &
         b% pg% Text_Summary_xleft(9), b% pg% Text_Summary_xright(9), &
         b% pg% Text_Summary_ybot(9), b% pg% Text_Summary_ytop(9), .false., &
         b% pg% Text_Summary_title(9), b% pg% Text_Summary_txt_scale(9), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Text_Summary9_plot


   subroutine do_Text_Summary9_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call Summary_plot(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
         b% pg% Text_Summary_num_rows(9), b% pg% Text_Summary_num_cols(9), &
         b% pg% Text_Summary_name(9, :, :), ierr)
   end subroutine do_Text_Summary9_plot


   subroutine Summary_plot(b, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
      Text_Summary_num_rows, Text_Summary_num_cols, &
      Text_Summary_name, ierr)

      use utils_lib
      use pgstar_support, only : write_info_line_exp, write_info_line_flt, &
         write_info_line_int

      type (binary_info), pointer :: b
      integer, intent(in) :: device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(in) :: Text_Summary_num_rows, Text_Summary_num_cols
      character (len = *), intent(in) :: Text_Summary_name(:, :)
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
      call pgswin(0.0, 1.0, 0.0, 1.0)
      call show_title_pgbinary(b, title)
      call pgsch(txt_scale * 0.8)

      do col = 1, num_cols
         call show_column(col, num_rows)
      end do

      call pgunsa


   contains


      subroutine show_column(col, num_rows)
         use binary_history, only : get_binary_history_specs, get_binary_history_values, get1_binary_hist_value
         integer, intent(in) :: col, num_rows

         real(dp) :: values(num_rows)
         integer :: int_values(num_rows), specs(num_rows)
         logical :: is_int_value(num_rows)
         logical :: failed_to_find_value(num_rows)

         integer :: i, cnt
         real :: xpos0, dxpos, dxval, ypos, dypos
         real(dp) :: val

         call get_binary_history_specs(&
            b, num_rows, Text_Summary_name(:, col), specs)
         call get_binary_history_values(&
            b, num_rows, specs, &
            is_int_value, int_values, values, failed_to_find_value)

         xpos0 = (real(col) - 0.5) / real(num_cols)

         dxpos = 0.00
         dxval = 0.02

         ypos = 0.90
         dypos = -0.95 / num_rows

         do i = 1, num_rows
            if (i > 1) ypos = ypos + dypos
            if (failed_to_find_value(i)) then
               if (.not. get1_binary_hist_value(b, Text_Summary_name(i, col), val)) then
                  if (len_trim(Text_Summary_name(i, col)) > 0) &
                     write(*, '(a)') 'failed_to_find_value ' // trim(Text_Summary_name(i, col)) &
                        // '. check that it is in your history_columns.list'
                  cycle
               end if
               values(i) = val
            else if (is_int_value(i)) then
               cnt = write_info_line_int(0, ypos, xpos0, dxpos, dxval, &
                  Text_Summary_name(i, col), int_values(i))
               cycle
            end if
            if (values(i) == 0d0) then
               cnt = write_info_line_int(0, ypos, xpos0, dxpos, dxval, &
                  Text_Summary_name(i, col), 0)
            else if (abs(values(i)) > 1d-3 .and. abs(values(i)) < 1d3) then
               cnt = write_info_line_flt(0, ypos, xpos0, dxpos, dxval, &
                  Text_Summary_name(i, col), values(i))
            else
               cnt = write_info_line_exp(0, ypos, xpos0, dxpos, dxval, &
                  Text_Summary_name(i, col), values(i))
            end if
         end do

      end subroutine show_column


   end subroutine Summary_plot


end module pgbinary_summary

