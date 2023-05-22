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

module pgbinary_support

   use binary_private_def
   use const_def
   use rates_def, only : i_rate
   use utils_lib
   use pgstar_support, only : Set_Colours, do1_pgmtxt, &
      clr_no_mixing, clr_convection, clr_leftover_convection, clr_semiconvection, &
      clr_thermohaline, clr_overshoot, clr_rotation, clr_minimum, clr_rayleigh_taylor, &
      clr_anonymous, colormap_offset, colormap_last, colormap_size, &
      colormap, Line_Type_Solid, Line_Type_Dash, Line_Type_Dash_Dot, Line_Type_Dot_Dash, &
      Line_Type_Dot ! inherit colors/linetypes + some routines from pgstar
   use star_pgstar

   implicit none

   logical :: have_initialized_pgbinary = .false.


contains


   subroutine add_to_pgbinary_hist(b, pg_hist_new)
      type (binary_info), pointer :: b
      type (pgbinary_hist_node), pointer :: pg_hist_new
      type (pgbinary_hist_node), pointer :: pg_hist => null(), next => null()
      integer :: step
      step = pg_hist_new% step
      do
         if (.not. associated(b% pg% pgbinary_hist)) then
            b% pg% pgbinary_hist => pg_hist_new
            nullify(pg_hist_new% next)
            return
         end if
         if (step > b% pg% pgbinary_hist% step) then
            pg_hist_new% next => b% pg% pgbinary_hist
            b% pg% pgbinary_hist => pg_hist_new
            return
         end if
         ! discard item
         next => b% pg% pgbinary_hist% next
         deallocate(b% pg% pgbinary_hist% vals)
         deallocate(b% pg% pgbinary_hist)
         b% pg% pgbinary_hist => next
      end do
   end subroutine add_to_pgbinary_hist


   subroutine pgbinary_clear(b)
      type (binary_info), pointer :: b
      integer :: i, num
      type (pgbinary_win_file_data), pointer :: p
      type (pgbinary_hist_node), pointer :: pg_hist => null(), next => null()
      pg_hist => b% pg% pgbinary_hist
      do while(associated(pg_hist))
         if (associated(pg_hist% vals)) deallocate(pg_hist% vals)
         next => pg_hist% next
         deallocate(pg_hist)
         pg_hist => next
      end do
      nullify(b% pg% pgbinary_hist)
      if (have_initialized_pgbinary) return
      do i = 1, num_pgbinary_plots
         p => b% pg% pgbinary_win_file_ptr(i)
         p% id_win = 0
         p% have_called_mkdir = .false.
         p% file_dir_for_previous_mkdir = ''
      end do
   end subroutine pgbinary_clear


   subroutine init_pgbinary(ierr)
      use pgstar_support, only : init_pgstar
      integer, intent(out) :: ierr

      call init_pgstar(ierr)

      if (ierr /= 0) then
         write(*, *) 'failed to init pgstar, required for pgbinary'
         return
      end if

      b% pg% star_1_color = clr_Goldenrod
      b% pg% star_2_color = clr_LightSkyBlue

      have_initialized_pgbinary = .true.
   end subroutine init_pgbinary

   subroutine check_window(b, p, ierr)
      type (binary_info), pointer :: b
      type (pgbinary_win_file_data), pointer :: p
      integer, intent(out) :: ierr
      ierr = 0
      if (p% do_win .and. (.not. p% win_flag)) then
         p% do_win = .false.
         if (p% id_win > 0) then
            call pgslct(p% id_win)
            call pgclos
            p% id_win = 0
         endif
      else if (p% win_flag .and. (.not. p% do_win)) then
         if (p% id_win == 0) &
            call open_device(b, p, .false., '/xwin', p% id_win, ierr)
         if (ierr == 0 .and. p% id_win > 0) p% do_win = .true.
      end if
      if (p% do_win .and. p% id_win > 0 .and. &
         (p% win_width /= p% prev_win_width .or. &
            p% win_aspect_ratio /= p% prev_win_aspect_ratio)) then
         call pgslct(p% id_win)
         call pgpap(p% win_width, p% win_aspect_ratio)
         p% prev_win_width = p% win_width
         p% prev_win_aspect_ratio = p% win_aspect_ratio
      end if
   end subroutine check_window


   subroutine check_file(b, p, ierr)
      use utils_lib, only : mkdir
      type (binary_info), pointer :: b
      type (pgbinary_win_file_data), pointer :: p
      integer, intent(out) :: ierr
      character (len = strlen) :: name
      ierr = 0
      if (p% do_file .and. (.not. p% file_flag)) then
         p% do_file = .false.
      else if (p% file_flag .and. (.not. p% do_file)) then
         if (p% id_file == 0) then
            if (.not. p% have_called_mkdir .or. &
               p% file_dir /= p% file_dir_for_previous_mkdir) then
               call mkdir(p% file_dir)
               p% have_called_mkdir = .true.
               p% file_dir_for_previous_mkdir = p% file_dir
            end if
            call create_file_name(b, p% file_dir, p% file_prefix, name)
            name = trim(name) // '/' // trim(b% pg% file_device)
            call open_device(b, p, .true., name, p% id_file, ierr)
            if (ierr /= 0) return
            p% most_recent_filename = name
         end if
         p% do_file = .true.
      end if
   end subroutine check_file


   subroutine create_file_name(b, dir, prefix, name)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: dir, prefix
      character (len = *), intent(out) :: name
      character (len = strlen) :: num_str, fstring
      write(fstring, '( "(i",i2.2,".",i2.2,")" )') &
         b% pg% file_digits, b% pg% file_digits
      write(num_str, fstring) b% model_number
      if (len_trim(dir) > 0) then
         name = trim(dir) // '/' // trim(prefix)
      else
         name = prefix
      end if
      name = trim(name) // trim(num_str) // '.' // trim(b% pg% file_extension)
   end subroutine create_file_name


   subroutine write_plot_to_file(b, p, filename, ierr)
      type (binary_info), pointer :: b
      type (pgbinary_win_file_data), pointer :: p
      character (len = *), intent(in) :: filename
      integer, intent(out) :: ierr
      character (len = strlen) :: name
      ierr = 0
      !name = trim(filename) // '/' // trim(b% file_device)
      name = trim(filename) // '/png'
      write(*, '(a)') 'write_plot_to_file device: ' // trim(name)
      call open_device(b, p, .true., trim(name), p% id_file, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in open_device'
         return
      end if
      call p% plot(b% binary_id, p% id_file, ierr)
      call pgclos
      p% id_file = 0
      p% do_file = .false.
   end subroutine write_plot_to_file


   subroutine open_device(b, p, is_file, dev, id, ierr)
      type (binary_info), pointer :: b
      type (pgbinary_win_file_data), pointer :: p
      logical, intent(in) :: is_file
      character (len = *), intent(in) :: dev
      integer, intent(out) :: id
      integer, intent(out) :: ierr

      integer :: pgopen, system
      character (len = strlen) :: dir, cmd
      logical :: white_on_black_flag
      real :: width, ratio

      if (is_file) then
         dir = p% file_dir
         white_on_black_flag = b% pg% file_white_on_black_flag
      else
         dir = ''
         white_on_black_flag = b% pg% win_white_on_black_flag
      end if

      ierr = 0
      id = -1
      id = pgopen(trim(dev))
      if (id <= 0) return

      !      write(*,*) 'open device <' // trim(dev) // '> ' // trim(p% name), id
      if (is_file) then
         width = p% file_width; if (width < 0) width = p% win_width
         ratio = p% file_aspect_ratio; if (ratio < 0) ratio = p% win_aspect_ratio
         call pgpap(width, ratio)
      else
         call pgpap(p% win_width, p% win_aspect_ratio)
         p% prev_win_width = p% win_width
         p% prev_win_aspect_ratio = p% win_aspect_ratio
      end if
      call Set_Colours(white_on_black_flag, ierr)
   end subroutine open_device


   integer function count_hist_points(b, step_min, step_max) result(numpts)
      type (binary_info), pointer :: b
      integer, intent(in) :: step_min, step_max
      type (pgbinary_hist_node), pointer :: pg
      include 'formats'
      numpts = 0
      pg => b% pg% pgbinary_hist
      do ! recall that hist list is decreasing by age (and step)
         if (.not. associated(pg)) return
         if (pg% step < step_min) return
         if (pg% step <= step_max .or. step_max <= 0) numpts = numpts + 1
         pg => pg% next
      end do
   end function count_hist_points


   logical function get1_hist_yvec(b, step_min, step_max, n, name, vec)
      use utils_lib, only : integer_dict_lookup
      type (binary_info), pointer :: b
      integer, intent(in) :: step_min, step_max, n ! n = count_hist_points
      character (len = *) :: name
      real, dimension(:), pointer :: vec
      integer :: i, ierr, cnt
      character (len = 64) :: key_name
      include 'formats'
      cnt = 0
      do i = 1, len(key_name)
         key_name(i:i) = ' '
      end do
      do i = 1, len_trim(name)
         if (name(i:i) == ' ') then
            cnt = cnt + 1
            key_name(i:i) = '_'
         else
            key_name(i:i) = name(i:i)
         end if
      end do
      call integer_dict_lookup(b% binary_history_names_dict, key_name, i, ierr)
      if (ierr /= 0 .or. i <= 0) then ! didn't find it
         get1_hist_yvec = .false.
         return
      end if
      call get_hist_points(b, step_min, step_max, n, i, vec)
      get1_hist_yvec = .true.
   end function get1_hist_yvec


   subroutine set_hist_points_steps(&
      b, step_min, step_max, numpts, vec, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: step_min, step_max, numpts
      real, intent(out) :: vec(:)
      integer, intent(out) :: ierr
      integer :: i, n
      type (pgbinary_hist_node), pointer :: pg
      ierr = 0
      if (numpts == 0) return
      pg => b% pg% pgbinary_hist
      i = numpts
      do ! recall that hist list is decreasing by age (and step)
         if (.not. associated(pg)) then
            ierr = -1
            return
         end if
         if (pg% step < step_min) then
            ierr = -1
            return
         end if
         if (pg% step <= step_max) then
            vec(i) = real(pg% step)
            i = i - 1
            if (i == 0) return
         end if
         pg => pg% next
      end do
   end subroutine set_hist_points_steps


   integer function get_hist_index(b, spec) result(index)
      type (binary_info), pointer :: b
      integer, intent(in) :: spec
      integer :: i, num
      ! note: this doesn't include "extra" columns
      num = size(b% binary_history_column_spec, dim = 1)
      do i = 1, num
         if (b% binary_history_column_spec(i) == spec) then
            index = i
            return
         end if
      end do
      index = -1
   end function get_hist_index


   subroutine get_hist_points(&
      b, step_min, step_max, numpts, index, vec)
      type (binary_info), pointer :: b
      integer, intent(in) :: step_min, step_max, numpts, index
      real, intent(out) :: vec(:)
      integer :: i, n
      type (pgbinary_hist_node), pointer :: pg => null()
      include 'formats'
      if (numpts == 0) return
      pg => b% pg% pgbinary_hist
      i = numpts
      vec = 0
      do ! recall that hist list is decreasing by age (and step)
         if (.not. associated(pg)) return
         if (pg% step < step_min) then
            ! this will not happen if have correct numpts
            return
         end if
         if (pg% step <= step_max .or. step_max <= 0) then
            if (.not. associated(pg% vals)) return
            if (size(pg% vals, dim = 1) < index) return
            vec(i) = pg% vals(index)
            i = i - 1
            if (i == 0) return
         end if
         pg => pg% next
      end do
   end subroutine get_hist_points


   subroutine show_annotations(b, show_annotation1, show_annotation2, show_annotation3)
      type (binary_info), pointer :: b
      logical, intent(in) :: show_annotation1, show_annotation2, show_annotation3
      if (show_annotation1 .and. len_trim(b% pg% annotation1_text) > 0) then
         call pgsci(b% pg% annotation1_ci)
         call pgscf(b% pg% annotation1_cf)
         call do1_pgmtxt(b% pg% annotation1_side, b% pg% annotation1_disp, &
            b% pg% annotation1_coord, b% pg% annotation1_fjust, b% pg% annotation1_text, &
            b% pg% annotation1_ch, b% pg% annotation1_lw)
      end if
      if (show_annotation2 .and. len_trim(b% pg% annotation2_text) > 0) then
         call pgsci(b% pg% annotation2_ci)
         call pgscf(b% pg% annotation2_cf)
         call do1_pgmtxt(b% pg% annotation2_side, b% pg% annotation2_disp, &
            b% pg% annotation2_coord, b% pg% annotation2_fjust, b% pg% annotation2_text, &
            b% pg% annotation2_ch, b% pg% annotation2_lw)
      end if
      if (show_annotation3 .and. len_trim(b% pg% annotation3_text) > 0) then
         call pgsci(b% pg% annotation3_ci)
         call pgscf(b% pg% annotation3_cf)
         call do1_pgmtxt(b% pg% annotation3_side, b% pg% annotation3_disp, &
            b% pg% annotation3_coord, b% pg% annotation3_fjust, b% pg% annotation3_text, &
            b% pg% annotation3_ch, b% pg% annotation3_lw)
      end if
   end subroutine show_annotations


   subroutine show_box_pgbinary(b, str1, str2)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: str1, str2
      real :: ch
      integer :: lw
      call pgqch(ch)
      call pgqlw(lw)
      call pgsch(b% pg% pgbinary_num_scale * ch)
      call pgslw(b% pg% pgbinary_box_lw)
      call pgbox(str1, 0.0, 0, str2, 0.0, 0)
      call pgsch(ch)
      call pgslw(lw)
   end subroutine show_box_pgbinary


   subroutine draw_rect()
      real, dimension(5) :: xs, ys
      call pgsave
      call pgsci(1)
      xs = (/0.0, 0.0, 1.0, 1.0, 0.0/)
      ys = (/0.0, 1.0, 1.0, 0.0, 0.0/)
      call pgswin(0.0, 1.0, 0.0, 1.0)
      call pgmove(0.0, 0.0)
      call pgline(5, xs, ys)
      call pgunsa
   end subroutine draw_rect


   subroutine show_grid_title_pgbinary(b, title, pad)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: title
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      if (.not. b% pg% pgbinary_grid_show_title) return
      if (len_trim(title) == 0) return
      call pgqch(ch)
      disp = b% pg% pgbinary_grid_title_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('T', disp, &
         b% pg% pgbinary_grid_title_coord, b% pg% pgbinary_grid_title_fjust, title, &
         b% pg% pgbinary_grid_title_scale * ch, b% pg% pgbinary_grid_title_lw)
   end subroutine show_grid_title_pgbinary


   subroutine show_title_pgbinary(b, title, pad)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: title
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      if (.not. b% pg% pgbinary_show_title) return
      if (len_trim(title) == 0) return
      call pgqch(ch)
      disp = b% pg% pgbinary_title_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('T', disp, &
         b% pg% pgbinary_title_coord, b% pg% pgbinary_title_fjust, title, &
         b% pg% pgbinary_title_scale * ch, b% pg% pgbinary_title_lw)
   end subroutine show_title_pgbinary


   subroutine show_title_label_pgmtxt_pgbinary(&
      b, coord, fjust, label, pad)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: disp
      disp = b% pg% pgbinary_title_disp
      if (present(pad)) disp = disp + pad
      call pgmtxt('T', disp, coord, fjust, label)
   end subroutine show_title_label_pgmtxt_pgbinary


   subroutine show_xaxis_label_pgbinary(b, label, pad)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      call pgqch(ch)
      disp = b% pg% pgbinary_xaxis_label_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('B', disp, 0.5, 0.5, label, &
         b% pg% pgbinary_xaxis_label_scale * ch, b% pg% pgbinary_xaxis_label_lw)
   end subroutine show_xaxis_label_pgbinary


   subroutine show_xaxis_label_pgmtxt_pgbinary(&
      b, coord, fjust, label, pad)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: disp
      disp = b% pg% pgbinary_xaxis_label_disp
      if (present(pad)) disp = disp + pad
      call pgmtxt('B', disp, coord, fjust, label)
   end subroutine show_xaxis_label_pgmtxt_pgbinary


   subroutine show_left_yaxis_label_pgbinary(b, label, pad)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      call pgqch(ch)
      disp = b% pg% pgbinary_left_yaxis_label_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('L', disp, 0.5, 0.5, label, &
         b% pg% pgbinary_left_yaxis_label_scale * ch, b% pg% pgbinary_left_yaxis_label_lw)
   end subroutine show_left_yaxis_label_pgbinary


   subroutine show_right_yaxis_label_pgbinary(b, label, pad)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      call pgqch(ch)
      disp = b% pg% pgbinary_right_yaxis_label_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('R', disp, 0.5, 0.5, label, &
         b% pg% pgbinary_right_yaxis_label_scale * ch, b% pg% pgbinary_right_yaxis_label_lw)
   end subroutine show_right_yaxis_label_pgbinary


   subroutine show_left_yaxis_label_pgmtxt_pgbinary(&
      b, coord, fjust, label, pad)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: ch, disp
      call pgqch(ch)
      call pgsch(1.1 * ch)
      disp = b% pg% pgbinary_left_yaxis_label_disp
      if (present(pad)) disp = disp + pad
      call pgmtxt('L', disp, coord, fjust, label)
      call pgsch(ch)
   end subroutine show_left_yaxis_label_pgmtxt_pgbinary


   subroutine show_right_yaxis_label_pgmtxt_pgbinary(&
      b, coord, fjust, label, pad)
      type (binary_info), pointer :: b
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: ch, disp
      call pgqch(ch)
      call pgsch(1.1 * ch)
      disp = b% pg% pgbinary_right_yaxis_label_disp
      if (present(pad)) disp = disp + pad
      call pgmtxt('R', disp, coord, fjust, label)
      call pgsch(ch)
   end subroutine show_right_yaxis_label_pgmtxt_pgbinary


   subroutine show_model_number_pgbinary(b)
      type (binary_info), pointer :: b
      character (len = 32) :: str
      real :: ch
      if (.not. b% pg% pgbinary_show_model_number) return
      write(str, '(i9)') b% model_number
      str = 'model ' // trim(adjustl(str))
      call pgqch(ch)
      call do1_pgmtxt('T', &
         b% pg% pgbinary_model_disp, b% pg% pgbinary_model_coord, &
         b% pg% pgbinary_model_fjust, str, &
         b% pg% pgbinary_model_scale * ch, b% pg% pgbinary_model_lw)
   end subroutine show_model_number_pgbinary


   subroutine show_age_pgbinary(b)
      type (binary_info), pointer :: b
      character (len = 32) :: age_str, units_str
      real(dp) :: age
      real :: ch
      integer :: len, i, j, iE, n
      if (.not. b% pg% pgbinary_show_age) return
      age = b% binary_age
      if (b% pg% pgbinary_show_age_in_seconds) then
         age = age * secyer
         units_str = 'secs'
      else if (b% pg% pgbinary_show_age_in_minutes) then
         age = age * secyer / 60
         units_str = 'mins'
      else if (b% pg% pgbinary_show_age_in_hours) then
         age = age * secyer / (60 * 60)
         units_str = 'hrs'
      else if (b% pg% pgbinary_show_age_in_days) then
         age = age * secyer / (60 * 60 * 24)
         units_str = 'days'
      else if (b% pg% pgbinary_show_age_in_years) then
         !age = age
         units_str = 'yrs'
      else if (b% pg% pgbinary_show_log_age_in_years) then
         age = log10(max(1d-99, age))
         units_str = 'log yrs'
      else if (age * secyer < 60) then
         age = age * secyer
         units_str = 'secs'
      else if (age * secyer < 60 * 60) then
         age = age * secyer / 60
         units_str = 'mins'
      else if (age * secyer < 60 * 60 * 24) then
         age = age * secyer / (60 * 60)
         units_str = 'hrs'
      else if (age * secyer < 60 * 60 * 24 * 500) then
         age = age * secyer / (60 * 60 * 24)
         units_str = 'days'
      else
         !age = age
         units_str = 'yrs'
      end if
      if (abs(age) > 1e-3 .and. abs(age) < 1e3) then
         write(age_str, '(f14.6)') age
      else
         write(age_str, '(1pe14.6)') age
         len = len_trim(age_str)
         iE = 0
         do i = 1, len
            if (age_str(i:i) == 'E') then
               iE = i
               age_str(i:i) = 'e'
               exit
            end if
         end do
         if (iE > 0) then
            i = iE + 1
            if (age_str(i:i) == '+') then
               do j = i, len - 1
                  age_str(j:j) = age_str(j + 1:j + 1)
               end do
               age_str(len:len) = ' '
               len = len - 1
            else
               i = i + 1
            end if
            if (age_str(i:i) == '0') then
               do j = i, len - 1
                  age_str(j:j) = age_str(j + 1:j + 1)
               end do
               age_str(len:len) = ' '
               len = len - 1
            end if
         end if
      end if
      age_str = adjustl(age_str)
      age_str = 'age ' // trim(age_str) // ' ' // trim(units_str)
      call pgqch(ch)
      call do1_pgmtxt('T', &
         b% pg% pgbinary_age_disp, b% pg% pgbinary_age_coord, &
         b% pg% pgbinary_age_fjust, age_str, &
         b% pg% pgbinary_age_scale * ch, b% pg% pgbinary_age_lw)
   end subroutine show_age_pgbinary


   logical function read_values_from_file(fname, x_data, y_data, data_len)
      character(len = *), intent(in) :: fname
      real, pointer, dimension(:) :: x_data, y_data
      integer, intent(out) :: data_len
      integer :: iounit, ierr, i
      include 'formats'
      read_values_from_file = .false.
      ierr = 0
      open(newunit = iounit, file = trim(fname), action = 'read', status = 'old', iostat = ierr)
      if (ierr /= 0) then
         !write(*, *) 'failed to open ' // trim(fname)
         return
      end if
      read(iounit, *, iostat = ierr) data_len
      if (ierr /= 0) then
         write(*, *) 'failed to read num points on 1st line ' // trim(fname)
         return
      end if
      !write(*,2) trim(fname) // ' data_len', data_len
      allocate(x_data(data_len), y_data(data_len))
      do i = 1, data_len
         read(iounit, *, iostat = ierr) x_data(i), y_data(i)
         if (ierr /= 0) then
            write(*, *) 'failed to read data ' // trim(fname)
            deallocate(x_data, y_data)
            return
         end if
      end do
      close(iounit)
      read_values_from_file = .true.
   end function read_values_from_file


   subroutine show_pgbinary_decorator(binary_id, use_flag, pgbinary_decorator, plot_num, ierr)
      logical, intent(in) :: use_flag
      real :: xmin, xmax, ymin, ymax
      integer, intent(in) :: binary_id, plot_num
      integer, intent(inout) :: ierr
      procedure(pgbinary_decorator_interface), pointer :: pgbinary_decorator

      if(use_flag)then
         if(associated(pgbinary_decorator))then
            call pgsave
            call PGQWIN(xmin, xmax, ymin, ymax)
            call pgbinary_decorator(binary_id, xmin, xmax, ymin, ymax, plot_num, ierr)
            call pgunsa
            if(ierr/=0)then
               write(*, *) "Error in pgbinary_decorator"
            end if
         end if
      end if

   end subroutine show_pgbinary_decorator

end module pgbinary_support

