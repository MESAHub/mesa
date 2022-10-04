! ***********************************************************************
!
!   Copyright (C) 2010-2022  Bill Paxton, Matthias Fabry
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
      call do_Text_Summary1_plot(b, id, device_id, &
         b% Text_Summary1_xleft, b% Text_Summary1_xright, &
         b% Text_Summary1_ybot, b% Text_Summary1_ytop, .false., &
         b% Text_Summary1_title, b% Text_Summary1_txt_scale, ierr)
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
         b% Text_Summary1_num_rows, b% Text_Summary1_num_cols, &
         b% Text_Summary1_name, ierr)
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
         b% Text_Summary2_xleft, b% Text_Summary2_xright, &
         b% Text_Summary2_ybot, b% Text_Summary2_ytop, .false., &
         b% Text_Summary2_title, b% Text_Summary2_txt_scale, ierr)
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
         b% Text_Summary2_num_rows, b% Text_Summary2_num_cols, &
         b% Text_Summary2_name, ierr)
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
         b% Text_Summary3_xleft, b% Text_Summary3_xright, &
         b% Text_Summary3_ybot, b% Text_Summary3_ytop, .false., &
         b% Text_Summary3_title, b% Text_Summary3_txt_scale, ierr)
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
         b% Text_Summary3_num_rows, b% Text_Summary3_num_cols, &
         b% Text_Summary3_name, ierr)
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
         b% Text_Summary4_xleft, b% Text_Summary4_xright, &
         b% Text_Summary4_ybot, b% Text_Summary4_ytop, .false., &
         b% Text_Summary4_title, b% Text_Summary4_txt_scale, ierr)
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
         b% Text_Summary4_num_rows, b% Text_Summary4_num_cols, &
         b% Text_Summary4_name, ierr)
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
         b% Text_Summary5_xleft, b% Text_Summary5_xright, &
         b% Text_Summary5_ybot, b% Text_Summary5_ytop, .false., &
         b% Text_Summary5_title, b% Text_Summary5_txt_scale, ierr)
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
         b% Text_Summary5_num_rows, b% Text_Summary5_num_cols, &
         b% Text_Summary5_name, ierr)
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
         b% Text_Summary6_xleft, b% Text_Summary6_xright, &
         b% Text_Summary6_ybot, b% Text_Summary6_ytop, .false., &
         b% Text_Summary6_title, b% Text_Summary6_txt_scale, ierr)
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
         b% Text_Summary6_num_rows, b% Text_Summary6_num_cols, &
         b% Text_Summary6_name, ierr)
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
         b% Text_Summary7_xleft, b% Text_Summary7_xright, &
         b% Text_Summary7_ybot, b% Text_Summary7_ytop, .false., &
         b% Text_Summary7_title, b% Text_Summary7_txt_scale, ierr)
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
         b% Text_Summary7_num_rows, b% Text_Summary7_num_cols, &
         b% Text_Summary7_name, ierr)
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
         b% Text_Summary8_xleft, b% Text_Summary8_xright, &
         b% Text_Summary8_ybot, b% Text_Summary8_ytop, .false., &
         b% Text_Summary8_title, b% Text_Summary8_txt_scale, ierr)
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
         b% Text_Summary8_num_rows, b% Text_Summary8_num_cols, &
         b% Text_Summary8_name, ierr)
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
         b% Text_Summary9_xleft, b% Text_Summary9_xright, &
         b% Text_Summary9_ybot, b% Text_Summary9_ytop, .false., &
         b% Text_Summary9_title, b% Text_Summary9_txt_scale, ierr)
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
         b% Text_Summary9_num_rows, b% Text_Summary9_num_cols, &
         b% Text_Summary9_name, ierr)
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

