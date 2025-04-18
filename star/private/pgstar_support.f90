! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module pgstar_support

   use star_private_def
   use const_def, only: dp, secday, secyer, mesa_data_dir, &
      overshoot_mixing, rotation_mixing, thermohaline_mixing, semiconvective_mixing, &
      leftover_convective_mixing, convective_mixing, no_mixing
   use rates_def, only : i_rate
   use utils_lib
   use star_pgstar

   implicit none

   integer, parameter :: category_offset = 1000
   integer, parameter :: abundance_offset = 2000
   integer, parameter :: extras_offset = 3000

   logical :: have_initialized_pgstar = .false.

   real :: sum_dHR_since_last_file_write


   ! lines for TRho profile
   real, dimension(:), allocatable :: hydrogen_burn_logT, hydrogen_burn_logRho
   real, dimension(:), allocatable :: helium_burn_logT, helium_burn_logRho
   real, dimension(:), allocatable :: carbon_burn_logT, carbon_burn_logRho
   real, dimension(:), allocatable :: oxygen_burn_logT, oxygen_burn_logRho
   real, dimension(:), allocatable :: psi4_logT, psi4_logRho
   real, dimension(:), allocatable :: elect_data_logT, elect_data_logRho
   real, dimension(:), allocatable :: gamma_4_thirds_logT, gamma_4_thirds_logRho
   real, dimension(:), allocatable :: kap_rad_cond_eq_logT, kap_rad_cond_eq_logRho
   real, dimension(:), allocatable :: opal_clip_logT, opal_clip_logRho
   real, dimension(:), allocatable :: scvh_clip_logT, scvh_clip_logRho

   ! Tioga line types
   integer, parameter :: Line_Type_Solid = 1
   integer, parameter :: Line_Type_Dash = 2
   integer, parameter :: Line_Type_Dot_Dash = 3
   integer, parameter :: Line_Type_Dot = 4

contains


   subroutine add_to_pgstar_hist(s, pg_hist_new)
      type (star_info), pointer :: s
      type (pgstar_hist_node), pointer :: pg_hist_new
      type (pgstar_hist_node), pointer :: next => null()
      integer :: step
      step = pg_hist_new% step
      do
         if (.not. associated(s% pg% pgstar_hist)) then
            s% pg% pgstar_hist => pg_hist_new
            nullify(pg_hist_new% next)
            return
         end if
         if (step > s% pg% pgstar_hist% step) then
            pg_hist_new% next => s% pg% pgstar_hist
            s% pg% pgstar_hist => pg_hist_new
            return
         end if
         ! discard item
         next => s% pg% pgstar_hist% next
         deallocate(s% pg% pgstar_hist% vals)
         deallocate(s% pg% pgstar_hist)
         s% pg% pgstar_hist => next
      end do
   end subroutine add_to_pgstar_hist


   subroutine pgstar_clear(s)
      type (star_info), pointer :: s
      integer :: i
      type (pgstar_win_file_data), pointer :: p
      type (pgstar_hist_node), pointer :: pg_hist => null(), next => null()
      pg_hist => s% pg% pgstar_hist
      do while(associated(pg_hist))
         if (associated(pg_hist% vals)) deallocate(pg_hist% vals)
         next => pg_hist% next
         deallocate(pg_hist)
         pg_hist => next
      end do
      nullify(s% pg% pgstar_hist)
      if (have_initialized_pgstar) return
      do i = 1, num_pgstar_plots
         p => s% pg% pgstar_win_file_ptr(i)
         p% id_win = 0
         p% have_called_mkdir = .false.
         p% file_dir_for_previous_mkdir = ''
      end do
   end subroutine pgstar_clear


   subroutine init_pgstar(ierr)
      integer, intent(out) :: ierr

      if (have_initialized_pgstar) return

      call read_support_info(ierr)
      if (failed('read_support_info')) return

      have_initialized_pgstar = .true.
      sum_dHR_since_last_file_write = 0

   contains

      logical function failed(str)
         character (len = *), intent(in) :: str
         failed = (ierr /= 0)
         if (failed) then
            write(*, *) trim(str) // ' ierr', ierr
         end if
      end function failed

   end subroutine init_pgstar

   subroutine check_window(s, p, ierr)
      type (star_info), pointer :: s
      type (pgstar_win_file_data), pointer :: p
      integer, intent(out) :: ierr
      ierr = 0
      if (p% do_win .and. (.not. p% win_flag)) then
         p% do_win = .false.
         if (p% id_win > 0) then
            call pgslct(p% id_win)
            call pgclos
            p% id_win = 0
         end if
      else if (p% win_flag .and. (.not. p% do_win)) then
         if (p% id_win == 0) &
            call open_device(s, p, .false., '/xwin', p% id_win, ierr)
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


   subroutine check_file(s, p, ierr)
      use utils_lib, only : mkdir
      type (star_info), pointer :: s
      type (pgstar_win_file_data), pointer :: p
      integer, intent(out) :: ierr
      character (len = strlen) :: name
      ierr = 0
      if (p% do_file .and. (.not. p% file_flag)) then
         p% do_file = .false.
      else if (p% file_flag .and. (.not. p% do_file)) then
         if (p% id_file == 0) then
            if (.not. p% have_called_mkdir .or. &
               p% file_dir /= p% file_dir_for_previous_mkdir) then
               if(.not. folder_exists(trim(p% file_dir))) call mkdir(trim(p% file_dir))
               p% have_called_mkdir = .true.
               p% file_dir_for_previous_mkdir = p% file_dir
            end if
            call create_file_name(s, p% file_dir, p% file_prefix, name)
            name = trim(name) // '/' // trim(s% pg% file_device)
            call open_device(s, p, .true., name, p% id_file, ierr)
            if (ierr /= 0) return
            p% most_recent_filename = name
         end if
         p% do_file = .true.
      end if
   end subroutine check_file


   subroutine create_file_name(s, dir, prefix, name)
      use star_utils, only : get_string_for_model_number
      type (star_info), pointer :: s
      character (len = *), intent(in) :: dir, prefix
      character (len = *), intent(out) :: name
      character (len = strlen) :: num_str, fstring
      character (len = 32) :: file_extension
      write(fstring, '( "(i",i2.2,".",i2.2,")" )') s% pg% file_digits, s% pg% file_digits
      write(num_str, fstring) s% model_number
      if (len_trim(dir) > 0) then
         name = trim(dir) // '/' // trim(prefix)
      else
         name = prefix
      end if
      if (s%pg%file_device=='vcps') then
         file_extension = 'ps'
      else
         file_extension = s%pg%file_device  ! e.g.: png, ps
      end if
      name = trim(name) // trim(num_str) // '.' // trim(file_extension)
   end subroutine create_file_name


   subroutine write_plot_to_file(s, p, filename, ierr)
      type (star_info), pointer :: s
      type (pgstar_win_file_data), pointer :: p
      character (len = *), intent(in) :: filename
      integer, intent(out) :: ierr
      character (len = strlen) :: name
      ierr = 0
      !name = trim(filename) // '/' // trim(s% pg% file_device)
      name = trim(filename) // '/png'
      write(*, '(a)') 'write_plot_to_file device: ' // trim(name)
      call open_device(s, p, .true., trim(name), p% id_file, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in open_device'
         return
      end if
      call p% plot(s% id, p% id_file, ierr)
      call pgclos
      p% id_file = 0
      p% do_file = .false.
   end subroutine write_plot_to_file


   subroutine open_device(s, p, is_file, dev, id, ierr)
      use pgstar_colors, only: set_device_colors
      type (star_info), pointer :: s
      type (pgstar_win_file_data), pointer :: p
      logical, intent(in) :: is_file
      character (len = *), intent(in) :: dev
      integer, intent(out) :: id
      integer, intent(out) :: ierr

      integer :: pgopen
      character (len = strlen) :: dir
      logical :: white_on_black_flag
      real :: width, ratio

      if (is_file) then
         dir = p% file_dir
         white_on_black_flag = s% pg% file_white_on_black_flag
      else
         dir = ''
         white_on_black_flag = s% pg% win_white_on_black_flag
      end if

      ierr = 0
      id = -1
      id = pgopen(trim(dev))
      if (id <= 0) return

      !write(*,*) 'open device <' // trim(dev) // '> ' // trim(p% name), id
      if (is_file) then
         width = p% file_width; if (width < 0) width = p% win_width
         ratio = p% file_aspect_ratio; if (ratio < 0) ratio = p% win_aspect_ratio
         call pgpap(width, ratio)
      else
         call pgpap(p% win_width, p% win_aspect_ratio)
         p% prev_win_width = p% win_width
         p% prev_win_aspect_ratio = p% win_aspect_ratio
      end if
      call set_device_colors(white_on_black_flag)
   end subroutine open_device


   subroutine read_support_info(ierr)
      integer, intent(out) :: ierr

      ierr = 0

      call read_TRho_data(&
         'hydrogen_burn.data', hydrogen_burn_logT, hydrogen_burn_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading hydrogen burn data'
         return
      end if

      call read_TRho_data(&
         'helium_burn.data', helium_burn_logT, helium_burn_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading helium burn data'
         return
      end if

      call read_TRho_data(&
         'carbon_burn.data', carbon_burn_logT, carbon_burn_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading carbon burn data'
         return
      end if

      call read_TRho_data(&
         'oxygen_burn.data', oxygen_burn_logT, oxygen_burn_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading oxygen burn data'
         return
      end if

      call read_TRho_data(&
         'psi4.data', psi4_logT, psi4_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading psi4 data'
         return
      end if

      call read_TRho_data(&
         'elect.data', elect_data_logT, elect_data_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading elect data'
         return
      end if

      call read_TRho_data(&
         'gamma_4_thirds.data', gamma_4_thirds_logT, gamma_4_thirds_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading gamma_4_thirds data'
         return
      end if

      call read_TRho_data(&
         'kap_rad_cond_eq.data', kap_rad_cond_eq_logT, kap_rad_cond_eq_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading kap_rad_cond_eq data'
         return
      end if

      call read_TRho_data(&
         'opal_clip.data', opal_clip_logT, opal_clip_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading opal_clip data'
         return
      end if

      call read_TRho_data(&
         'scvh_clip.data', scvh_clip_logT, scvh_clip_logRho, ierr)
      if (ierr /= 0) then
         write(*, *) 'PGSTAR failed in reading scvh_clip data'
         return
      end if

   end subroutine read_support_info

   subroutine read_TRho_data(fname, logTs, logRhos, ierr)
      use utils_lib
      character (len = *), intent(in) :: fname
      real, dimension(:), allocatable :: logTs, logRhos  ! will allocate
      integer, intent(out) :: ierr

      character (len = strlen) :: filename
      real :: logT, logRho
      integer :: iounit, i, sz, cnt

      filename = trim(mesa_data_dir) // '/star_data/plot_info/' // trim(fname)

      open(newunit = iounit, file = trim(filename), status = 'old', action = 'read', iostat = ierr)
      if (ierr/=0) then
         write(*, *) 'failed to open ' // trim(filename)
         call done
         return
      end if

      sz = 0
      do
         read(iounit, *, iostat = ierr)
         if(ierr/=0) exit
         sz = sz + 1
      end do

      rewind(iounit)

      allocate(logTs(sz), logRhos(sz))

      cnt = 0
      do i = 1, sz
         read(iounit, *, iostat = ierr) logRho, logT
         if (ierr /= 0) then
            ierr = 0; exit
         end if
         logRhos(i) = logRho
         logTs(i) = logT
         cnt = i
      end do

      call done


   contains


      subroutine done
         close(iounit)
      end subroutine done

   end subroutine read_TRho_data


   integer function write_info_line_str(cnt, ypos, xpos0, dxpos, str)
      integer, intent(in) :: cnt
      real, intent(in) :: ypos, xpos0, dxpos
      character (len = *), intent(in) :: str
      real :: xpos
      xpos = cnt * dxpos + xpos0
      call pgptxt(xpos, ypos, 0.0, 0.5, trim(adjustl(str)))
      write_info_line_str = cnt + 1
   end function write_info_line_str


   integer function write_info_line_int(cnt, ypos, xpos0, dxpos, dxval, label, val)
      integer, intent(in) :: cnt, val
      real, intent(in) :: ypos, xpos0, dxpos, dxval
      character (len = *), intent(in) :: label

      character (len = 128) :: str
      real :: xpos

      write(str, '(a)') trim(label)
      xpos = cnt * dxpos + xpos0
      call pgptxt(xpos, ypos, 0.0, 1.0, trim(adjustl(str)))
      write(str, '(i9)') val
      xpos = xpos + dxval
      call pgptxt(xpos, ypos, 0.0, 0.0, trim(adjustl(str)))

      write_info_line_int = cnt + 1
   end function write_info_line_int


   integer function write_info_line_flt(&
      cnt, ypos, xpos0, dxpos, dxval, label, val)
      integer, intent(in) :: cnt
      real, intent(in) :: ypos, xpos0, dxpos, dxval
      real(dp), intent(in) :: val
      character (len = *), intent(in) :: label

      character (len = 128) :: str
      real :: xpos
      integer :: ierr

      write_info_line_flt = cnt + 1
      write(str, '(a)')   trim(label)
      xpos = cnt * dxpos + xpos0
      call pgptxt(xpos, ypos, 0.0, 1.0, trim(adjustl(str)))
      ierr = 0
      write(str, '(f12.7)', iostat = ierr) val
      if (ierr /= 0) then
         ierr = 0
         write(str, '(e10.3)', iostat = ierr) val
         if (ierr /= 0) then
            write(*, *) trim(label), val
            write(*, *) 'problem in write_info_line_flt'
            return
         end if
      end if
      xpos = xpos + dxval
      call pgptxt(xpos, ypos, 0.0, 0.0, trim(adjustl(str)))

   end function write_info_line_flt


   integer function write_info_line_flt2(cnt, ypos, xpos0, dxpos, dxval, label, val)
      integer, intent(in) :: cnt
      real, intent(in) :: ypos, xpos0, dxpos, dxval
      real(dp), intent(in) :: val
      character (len = *), intent(in) :: label

      character (len = 128) :: str
      real :: xpos
      integer :: ierr

      write_info_line_flt2 = cnt + 1

      write(str, '(a)')   trim(label)
      xpos = cnt * dxpos + xpos0
      call pgptxt(xpos, ypos, 0.0, 1.0, trim(adjustl(str)))
      ierr = 0
      write(str, '(f12.3)', iostat = ierr) val
      if (ierr /= 0) then
         ierr = 0
         write(str, '(e10.3)', iostat = ierr) val
         if (ierr /= 0) then
            write(*, *) trim(label), val
            write(*, *) 'problem in write_info_line_flt2'
            return
         end if
      end if
      xpos = xpos + dxval
      call pgptxt(xpos, ypos, 0.0, 0.0, trim(adjustl(str)))

   end function write_info_line_flt2


   integer function write_info_line_exp(cnt, ypos, xpos0, dxpos, dxval, label, val)
      integer, intent(in) :: cnt
      real, intent(in) :: ypos, xpos0, dxpos, dxval
      real(dp), intent(in) :: val
      character (len = *), intent(in) :: label

      character (len = 128) :: str
      real :: xpos

      write(str, '(a)')   trim(label)
      xpos = cnt * dxpos + xpos0
      call pgptxt(xpos, ypos, 0.0, 1.0, trim(adjustl(str)))
      write(str, '(1pe10.3)') val
      xpos = xpos + dxval
      call pgptxt(xpos, ypos, 0.0, 0.0, trim(adjustl(str)))

      write_info_line_exp = cnt + 1
   end function write_info_line_exp


   integer function count_hist_points(&
      s, step_min, step_max) result(numpts)
      type (star_info), pointer :: s
      integer, intent(in) :: step_min, step_max
      type (pgstar_hist_node), pointer :: pg
      include 'formats'
      numpts = 0
      pg => s% pg% pgstar_hist
      do  ! recall that hist list is decreasing by age (and step)
         if (.not. associated(pg)) return
         if (pg% step < step_min) return
         if (pg% step <= step_max .or. step_max <= 0) numpts = numpts + 1
         pg => pg% next
      end do
   end function count_hist_points


   logical function get1_hist_yvec(s, step_min, step_max, n, name, vec)
      use utils_lib, only : integer_dict_lookup
      type (star_info), pointer :: s
      integer, intent(in) :: step_min, step_max, n  ! n = count_hist_points
      character (len = *) :: name
      real, dimension(:), allocatable :: vec
      integer :: i, cnt, ierr
      character (len = 64) :: key_name
      include 'formats'
      cnt = 0
      ierr = 0
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
      call integer_dict_lookup(s% history_names_dict, key_name, i, ierr)
      if (ierr /= 0 .or. i <= 0) then  ! didn't find it
         get1_hist_yvec = .false.
         return
      end if
      call get_hist_points(s, step_min, step_max, n, i, vec, ierr)
      if (ierr /= 0) then  ! didn't get them
         get1_hist_yvec = .false.
         return
      end if
      get1_hist_yvec = .true.
   end function get1_hist_yvec


   subroutine set_hist_points_steps(&
      s, step_min, step_max, numpts, vec, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: step_min, step_max, numpts
      real, intent(out) :: vec(:)
      integer, intent(out) :: ierr
      integer :: i
      type (pgstar_hist_node), pointer :: pg
      ierr = 0
      if (numpts == 0) return
      pg => s% pg% pgstar_hist
      i = numpts
      do  ! recall that hist list is decreasing by age (and step)
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


   integer function get_hist_index(s, spec) result(index)
      type (star_info), pointer :: s
      integer, intent(in) :: spec
      integer :: i, num
      ! note: this doesn't include "extra" columns
      num = size(s% history_column_spec, dim = 1)
      do i = 1, num
         if (s% history_column_spec(i) == spec) then
            index = i
            return
         end if
      end do
      index = -1
   end function get_hist_index


   subroutine get_hist_points(&
      s, step_min, step_max, numpts, index, vec, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: step_min, step_max, numpts, index
      real, intent(out) :: vec(:)
      integer, intent(out) :: ierr
      integer :: i
      type (pgstar_hist_node), pointer :: pg => null()
      include 'formats'
      if (numpts == 0) return
      pg => s% pg% pgstar_hist
      i = numpts
      vec = 0
      ierr = 0
      do  ! recall that hist list is decreasing by age (and step)
         if (.not. associated(pg)) return
         if (pg% step < step_min) then
            ! this will not happen if have correct numpts
            return
         end if
         if (pg% step <= step_max .or. step_max <= 0) then
            if (.not. associated(pg% vals)) then
               !ierr = -1
               !write(*,6) 'failed in get_hist_points: not associated', &
               !   s% model_number, index, numpts, step_min, step_max
               !call mesa_error(__FILE__,__LINE__,'get_hist_points')
               return
            end if
            if (size(pg% vals, dim = 1) < index) then
               !ierr = -1
               !write(*,7) 'failed in get_hist_points: size < index', &
               !   s% model_number, size(pg% vals,dim=1), index, numpts, step_min, step_max
               !call mesa_error(__FILE__,__LINE__,'get_hist_points')
               return
            end if
            vec(i) = pg% vals(index)
            i = i - 1
            if (i == 0) return
         end if
         pg => pg% next
      end do
   end subroutine get_hist_points


   logical function find_shock(s, xaxis_id, xshock) result(found_shock)
      use profile, only : get_profile_val
      use num_lib, only : find0
      type (star_info), pointer :: s
      integer, intent(in) :: xaxis_id
      real(dp), intent(out) :: xshock
      integer :: k, nz
      real(dp) :: cs, x00, xp1, ms
      real(dp), pointer :: v(:) => null()
      include 'formats'
      nz = s% nz
      if (s% u_flag) then
         v => s% u
      else
         v => s% v
      end if
      ! search in from surface
      do k = 1, nz - 1
         cs = s% csound(k)  ! cell center
         if (v(k + 1) >= cs .and. v(k) < cs) then
            found_shock = .true.
            exit
         end if
      end do
      if (.not. found_shock) then
         ! search out from center
         do k = nz - 1, 1, -1
            cs = s% csound(k)  ! cell center
            if (v(k + 1) >= -cs .and. v(k) < -cs) then
               found_shock = .true.
               exit
            end if
         end do
      end if
      if (found_shock) then
         x00 = get_profile_val(s, xaxis_id, k)
         xp1 = get_profile_val(s, xaxis_id, k + 1)
         cs = s% csound(k)
         ms = find0(0d0, s% dm(k), v(k + 1) - cs, v(k) - cs)
         xshock = xp1 + (x00 - xp1) * ms / s% dm(k)
         if (is_bad(xshock)) then
            write(*, 2) 'xshock nz', nz, xshock
            write(*, 2) 's% dm(k)', k, s% dm(k)
            write(*, 2) 's% csound(k)', k, s% csound(k)
            write(*, 2) 'v(k)', k, v(k)
            write(*, 2) 'v(k+1)', k + 1, v(k + 1)
            write(*, 1) 'ms', ms
            write(*, 1) 'x00', x00
            write(*, 1) 'xp1', xp1
            nullify(v)
            call mesa_error(__FILE__, __LINE__, 'find_shock')
         end if
      end if
      nullify(v)
   end function find_shock


   subroutine set_grid_minmax(&
      nz, xvec, xmin, xmax, xleft, xright, xaxis_by, &
      given_xmin, given_xmax, margin_in, reversed, &
      grid_min, grid_max, dxmin)
      integer, intent(in) :: nz
      real, intent(in), dimension(:) :: xvec
      real, intent(out) :: xmin, xmax, xleft, xright
      character (len = *), intent(in) :: xaxis_by
      real, intent(in) :: given_xmin, given_xmax, margin_in
      real, intent(in) :: dxmin
      logical, intent(in) :: reversed
      integer, intent(out) :: grid_min, grid_max
      integer :: k
      real :: dx, margin
      logical :: use_given_xmin, use_given_xmax

      include 'formats'

      margin = max(0.0, margin_in)

      ! use given if it isn't = 101
      use_given_xmin = abs(given_xmin + 101.0) > 1e-6
      if (xaxis_by == 'mass' .and. given_xmin < 0 .and. use_given_xmin) then
         xmin = maxval(xvec(1:nz)) + given_xmin
      else if (use_given_xmin) then
         xmin = given_xmin
      else if (xaxis_by == 'logxm' .or. xaxis_by == 'logxq') then
         xmin = minval(xvec(2:nz))
      else
         xmin = minval(xvec(1:nz))
      end if

      use_given_xmax = abs(given_xmax + 101.0) > 1e-6
      if (xaxis_by == 'mass' .and. given_xmax < 0 .and. use_given_xmax) then
         xmax = maxval(xvec(1:nz)) + given_xmax
      else if (use_given_xmax) then
         xmax = given_xmax
      else
         xmax = maxval(xvec(1:nz))
      end if
      dx = xmax - xmin

      if (.not. use_given_xmin) xmin = xmin - margin * dx
      if (.not. use_given_xmax) xmax = xmax + margin * dx

      dx = xmax - xmin
      if (dx < dxmin) then
         dx = dxmin
         xmax = (xmax + xmin) / 2 + dx / 2
         xmin = xmax - dx
      end if

      if (xmin == xmax) then
         xmin = xmin - margin / 2
         xmax = xmax + margin / 2
      end if

      if (reversed) then
         xright = xmin; xleft = xmax
      else
         xright = xmax; xleft = xmin
      end if

      if (xvec(1) < xvec(nz)) then  ! increasing xs
         grid_max = nz
         do k = nz - 1, 1, -1  ! in decreasing order
            if (xvec(k) < xmax) then  ! this is the first one < xmax
               grid_max = k + 1
               exit
            end if
         end do
         grid_min = 1
         do k = grid_max, 1, -1
            if (xvec(k) <= xmin) then  ! this is the first one <= xmin
               grid_min = k
               exit
            end if
         end do
      else  ! decreasing
         grid_min = 1
         do k = 2, nz  ! in decreasing order
            if (xvec(k) < xmax) then  ! this is the first one < xmax
               grid_min = k - 1
               exit
            end if
         end do
         grid_max = nz
         do k = grid_min, nz
            if (xvec(k) <= xmin) then  ! this is the first one <= xmin
               grid_max = k
               exit
            end if
         end do
      end if

   end subroutine set_grid_minmax


   subroutine set_xleft_xright(&
      npts, xvec, given_xmin, given_xmax, xmargin_in, &
      reversed, dxmin, xleft, xright)
      integer, intent(in) :: npts
      real, intent(in), dimension(:) :: xvec
      real, intent(in) :: given_xmin, given_xmax, xmargin_in
      logical, intent(in) :: reversed
      real, intent(in) :: dxmin
      real, intent(out) :: xleft, xright
      call set_ytop_ybot(&
         npts, xvec, given_xmin, given_xmax, -101.0, xmargin_in, &
         reversed, dxmin, xleft, xright)
   end subroutine set_xleft_xright


   subroutine set_ytop_ybot(&
      npts, yvec, given_ymin, given_ymax, given_ycenter, &
      ymargin_in, reversed, dymin, ybot, ytop)
      integer, intent(in) :: npts
      real, intent(in), dimension(:) :: yvec
      real, intent(in) :: given_ymin, given_ymax, given_ycenter, ymargin_in
      logical, intent(in) :: reversed
      real, intent(in) :: dymin
      real, intent(out) :: ybot, ytop

      real :: dy, ymax, ymin, ymargin
      logical :: use_given_ymin, use_given_ymax
      real, parameter :: dymin_min = 1d-34

      include 'formats'

      ymargin = max(0.0, ymargin_in)

      use_given_ymin = abs(given_ymin + 101.0) > 1e-6
      if (use_given_ymin) then
         ymin = given_ymin
      else
         ymin = minval(yvec(1:npts))
      end if

      use_given_ymax = abs(given_ymax + 101.0) > 1e-6
      if (use_given_ymax) then
         ymax = given_ymax
      else
         ymax = maxval(yvec(1:npts))
      end if
      dy = ymax - ymin

      if (.not. use_given_ymin) ymin = ymin - ymargin * dy
      if (.not. use_given_ymax) ymax = ymax + ymargin * dy

      if (abs(given_ycenter + 101.0) > 1e-6 .and. &
         ymax > given_ycenter .and. given_ycenter > ymin) then
         dy = 2d0 * max(given_ycenter - ymin, ymax - given_ycenter)
         ymax = given_ycenter + 0.5d0 * dy
         ymin = given_ycenter - 0.5d0 * dy
      end if

      dy = ymax - ymin
      if (dy == 0.0) then
         ymax = ymax + dy
         ymin = ymin - dy
      else if (dy < max(dymin, dymin_min)) then
         !dymin_min prevents graphical glitches and segmentation faults when dy is too small
         dy = max(dymin, dymin_min)
         ymax = (ymax + ymin) / 2 + dy / 2
         ymin = ymax - dy
      end if

      if (ymin == ymax) then
         ymin = ymin - max(1.0, ymargin) / 2
         ymax = ymax + max(1.0, ymargin) / 2
      end if

      if (ymin == ymax) then  ! round off problems
         dy = 1e-6 * abs(ymax)
         ymin = ymin - dy
         ymax = ymax + dy
      end if

      if (reversed) then
         ytop = ymin; ybot = ymax
      else
         ytop = ymax; ybot = ymin
      end if

   end subroutine set_ytop_ybot


   subroutine do1_pgmtxt(side, disp, coord, fjust, label, ch, lw)
      character (len = *), intent(in) :: side, label
      real, intent(in) :: disp, coord, fjust, ch
      integer, intent(in) :: lw
      real :: sav_ch
      integer :: sav_lw
      call pgqch(sav_ch)
      call pgqlw(sav_lw)
      call pgslw(lw)
      call pgsch(ch)
      call pgmtxt(side, disp, coord, fjust, label)
      call pgslw(sav_lw)
      call pgsch(sav_ch)
   end subroutine do1_pgmtxt


   subroutine show_annotations(s, show_annotation1, show_annotation2, show_annotation3)
      type (star_info), pointer :: s
      logical, intent(in) :: show_annotation1, show_annotation2, show_annotation3
      if (show_annotation1 .and. len_trim(s% pg% annotation1_text) > 0) then
         call pgsci(s% pg% annotation1_ci)
         call pgscf(s% pg% annotation1_cf)
         call do1_pgmtxt(s% pg% annotation1_side, s% pg% annotation1_disp, &
            s% pg% annotation1_coord, s% pg% annotation1_fjust, s% pg% annotation1_text, &
            s% pg% annotation1_ch, s% pg% annotation1_lw)
      end if
      if (show_annotation2 .and. len_trim(s% pg% annotation2_text) > 0) then
         call pgsci(s% pg% annotation2_ci)
         call pgscf(s% pg% annotation2_cf)
         call do1_pgmtxt(s% pg% annotation2_side, s% pg% annotation2_disp, &
            s% pg% annotation2_coord, s% pg% annotation2_fjust, s% pg% annotation2_text, &
            s% pg% annotation2_ch, s% pg% annotation2_lw)
      end if
      if (show_annotation3 .and. len_trim(s% pg% annotation3_text) > 0) then
         call pgsci(s% pg% annotation3_ci)
         call pgscf(s% pg% annotation3_cf)
         call do1_pgmtxt(s% pg% annotation3_side, s% pg% annotation3_disp, &
            s% pg% annotation3_coord, s% pg% annotation3_fjust, s% pg% annotation3_text, &
            s% pg% annotation3_ch, s% pg% annotation3_lw)
      end if
   end subroutine show_annotations

   subroutine set_xaxis_bounds(&
      s, xaxis_by, win_xmin_in, win_xmax_in, xaxis_reversed, xmargin, &
      xvec, xmin, xmax, xleft, xright, dx, &
      grid_min, grid_max, npts, ierr)
      use profile_getval, only : get_profile_id, get_profile_val
      type (star_info), pointer :: s
      character (len = *), intent(in) :: xaxis_by
      real, intent(in) :: win_xmin_in, win_xmax_in, xmargin
      logical, intent(in) :: xaxis_reversed
      real, allocatable, dimension(:) :: xvec
      real, intent(out) :: xmin, xmax, xleft, xright, dx
      integer, intent(out) :: grid_min, grid_max, npts
      integer, intent(out) :: ierr

      integer :: k, nz, xaxis_id
      real :: win_xmin, win_xmax

      include 'formats'

      ierr = 0
      win_xmin = win_xmin_in
      win_xmax = win_xmax_in
      nz = s% nz

      xaxis_id = get_profile_id(s, xaxis_by)
      if (xaxis_id <= 0) then
         write(*, '(a)') &
            'pgstar inlist problem: bad value for xaxis_by: <' // trim(xaxis_by) // '>'
         ierr = -1
         return
      end if
      do k = 1, nz
         xvec(k) = get_profile_val(s, xaxis_id, k)
      end do

      call set_grid_minmax(&
         nz, xvec, xmin, xmax, xleft, xright, xaxis_by, &
         win_xmin, win_xmax, xmargin, xaxis_reversed, grid_min, grid_max, 0.0)
      dx = xmax - xmin
      npts = grid_max - grid_min + 1
      if (npts <= 0) then
         write(*, *) 'invalid x axis bounds for xaxis_by = ' // trim(xaxis_by)
         write(*, 1) 'xmax', xmax
         write(*, 1) 'xmin', xmin
         write(*, 1) 'dx', dx
         write(*, 1) 'xleft', xleft
         write(*, 1) 'xright', xright
         write(*, 1) 'win_xmin', win_xmin
         write(*, 1) 'win_xmax', win_xmax
         write(*, 1) 'xmargin', xmargin
         write(*, 2) 'grid_min', grid_min
         write(*, 2) 'grid_max', grid_max
         write(*, 2) 'npts', npts
         write(*, 2) 'nz', nz
         write(*, 1) 'maxval(xvec(1:nz))', maxval(xvec(1:nz))
         write(*, 1) 'minval(xvec(1:nz))', minval(xvec(1:nz))
         ierr = -1
      end if

   end subroutine set_xaxis_bounds


   subroutine show_xaxis_name(s, name, ierr)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: name
      integer, intent(out) :: ierr
      ierr = 0
      if (name == 'mass') then
         call show_xaxis_label_pgstar(s, 'm/M\d\(2281)')
      else if (name == 'grid') then
         call show_xaxis_label_pgstar(s, 'grid')
      else if (name == 'radius') then
         call show_xaxis_label_pgstar(s, 'r/R\d\(2281)')
      else if (name == 'logR') then
         call show_xaxis_label_pgstar(s, 'log r/R\d\(2281)')
      else if (name == 'logT') then
         call show_xaxis_label_pgstar(s, 'log T')
      else if (name == 'logP') then
         call show_xaxis_label_pgstar(s, 'log P')
      else if (name == 'logxq') then
         if (s% M_center == 0) then
            call show_xaxis_label_pgstar(s, 'log(1-q) q=fraction of total mass')
         else
            call show_xaxis_label_pgstar(s, 'log(1-q) q=fraction of envelope mass')
         end if
      else if (name == 'logxm') then
         call show_xaxis_label_pgstar(s, 'log((Mstar-m)/M\d\(2281)\u)')
      else if (name == 'r_div_R') then
         call show_xaxis_label_pgstar(s, 'r/R')
      else if (name == 'log_column_depth') then
         call show_xaxis_label_pgstar(s, 'log column depth (g cm\u-2\d)')
      else
         call show_xaxis_label_pgstar(s, name)
      end if
   end subroutine show_xaxis_name


   subroutine show_mix_regions_on_xaxis(s, ybot_in, ytop, grid_min, grid_max, xvec)
      type (star_info), pointer :: s
      real, intent(in) :: ybot_in, ytop
      integer, intent(in) :: grid_min, grid_max
      real, allocatable, dimension(:) :: xvec
      real :: ybot
      ybot = ybot_in + 0.001 * (ytop - ybot_in)
      call show_no_mixing_section(s, ybot, grid_min, grid_max, xvec)
      call show_convective_section(s, ybot, grid_min, grid_max, xvec)
      call show_leftover_convective_section(s, ybot, grid_min, grid_max, xvec)
      call show_semiconvective_section(s, ybot, grid_min, grid_max, xvec)
      call show_thermohaline_section(s, ybot, grid_min, grid_max, xvec)
      call show_rotation_section(s, ybot, grid_min, grid_max, xvec)
      call show_overshoot_section(s, ybot, grid_min, grid_max, xvec)
   end subroutine show_mix_regions_on_xaxis


   subroutine show_no_mixing_section(s, ybot, grid_min, grid_max, xvec)
      use pgstar_colors
      type (star_info), pointer :: s
      real, intent(in) :: ybot
      integer, intent(in) :: grid_min, grid_max
      real, allocatable, dimension(:) :: xvec
      call show_mixing_section(s, ybot, grid_min, grid_max, xvec, no_mixing, clr_no_mixing)
   end subroutine show_no_mixing_section


   subroutine show_convective_section(s, ybot, grid_min, grid_max, xvec)
      use pgstar_colors
      type (star_info), pointer :: s
      real, intent(in) :: ybot
      integer, intent(in) :: grid_min, grid_max
      real, allocatable, dimension(:) :: xvec
      call show_mixing_section(&
         s, ybot, grid_min, grid_max, xvec, convective_mixing, clr_convection)
   end subroutine show_convective_section


   subroutine show_leftover_convective_section(s, ybot, grid_min, grid_max, xvec)
      use pgstar_colors
      type (star_info), pointer :: s
      real, intent(in) :: ybot
      integer, intent(in) :: grid_min, grid_max
      real, allocatable, dimension(:) :: xvec
      call show_mixing_section(&
         s, ybot, grid_min, grid_max, xvec, leftover_convective_mixing, clr_leftover_convection)
   end subroutine show_leftover_convective_section


   subroutine show_semiconvective_section(s, ybot, grid_min, grid_max, xvec)
      use pgstar_colors
      type (star_info), pointer :: s
      real, intent(in) :: ybot
      integer, intent(in) :: grid_min, grid_max
      real, allocatable, dimension(:) :: xvec
      call show_mixing_section(&
         s, ybot, grid_min, grid_max, xvec, semiconvective_mixing, clr_semiconvection)
   end subroutine show_semiconvective_section


   subroutine show_thermohaline_section(s, ybot, grid_min, grid_max, xvec)
      use pgstar_colors
      type (star_info), pointer :: s
      real, intent(in) :: ybot
      integer, intent(in) :: grid_min, grid_max
      real, allocatable, dimension(:) :: xvec
      call show_mixing_section(&
         s, ybot, grid_min, grid_max, xvec, thermohaline_mixing, clr_thermohaline)
   end subroutine show_thermohaline_section


   subroutine show_rotation_section(s, ybot, grid_min, grid_max, xvec)
      use pgstar_colors
      type (star_info), pointer :: s
      real, intent(in) :: ybot
      integer, intent(in) :: grid_min, grid_max
      real, allocatable, dimension(:) :: xvec
      call show_mixing_section(&
         s, ybot, grid_min, grid_max, xvec, rotation_mixing, clr_rotation)
   end subroutine show_rotation_section


   subroutine show_overshoot_section(s, ybot, grid_min, grid_max, xvec)
      use pgstar_colors
      type (star_info), pointer :: s
      real, intent(in) :: ybot
      integer, intent(in) :: grid_min, grid_max
      real, allocatable, dimension(:) :: xvec
      call show_mixing_section(&
         s, ybot, grid_min, grid_max, xvec, overshoot_mixing, clr_overshoot)
   end subroutine show_overshoot_section


   subroutine show_mixing_section(s, ybot, grid_min, grid_max, xvec, mixing_type, clr)
      type (star_info), pointer :: s
      real, intent(in) :: ybot
      real, allocatable, dimension(:) :: xvec
      integer, intent(in) :: mixing_type, clr, grid_min, grid_max

      integer :: k, first, last
      logical :: inside
      include 'formats'
      inside = (s% mixing_type(grid_min) == mixing_type)
      first = grid_min
      call pgsci(clr)
      do k = grid_min, grid_max  ! 2,s% nz
         if (.not. inside) then
            if (s% mixing_type(k) == mixing_type) then  ! starting
               inside = .true.
               first = k
            end if
         else  ! inside
            if (s% mixing_type(k) /= mixing_type) then  ! ending
               last = k - 1
               call pgmove(xvec(first), ybot)
               call pgdraw(xvec(last), ybot)
               inside = .false.
            end if
         end if
      end do
      if (inside) then
         last = grid_max
         call pgmove(xvec(first), ybot)
         call pgdraw(xvec(last), ybot)
      end if
   end subroutine show_mixing_section


   subroutine show_profile_line(&
      s, xvec, yvec, txt_scale, xmin, xmax, ymin, ymax, &
      show_legend, legend_coord, legend_disp1, legend_del_disp, legend_fjust, &
      show_mass_pts)
      use pgstar_colors
      type (star_info), pointer :: s
      real, intent(in) :: xvec(:), yvec(:), txt_scale, xmin, xmax, ymin, ymax, &
         legend_coord, legend_disp1, legend_del_disp, legend_fjust
      logical, intent(in) :: show_legend, show_mass_pts

      real :: disp
      integer :: nz
      logical :: has_convection, has_leftover_convection, has_overshoot, &
         has_semiconvection, has_thermohaline, has_rotation

      include 'formats'

      call pgsave

      nz = s% nz
      call pgsch(s% pg% TRho_Profile_legend_txt_scale * txt_scale)

      call pgsci(clr_Gold)
      call pgslw(14)
      call do_show_eps_nuc_section(1d0)
      call pgslw(1)
      disp = legend_disp1 + 2 * legend_del_disp
      if (show_legend) &
         call pgmtxt('T', disp, legend_coord, legend_fjust, '> 1 erg g\u-1\d s\u-1\d')

      call pgsci(clr_Coral)
      call pgslw(18)
      call do_show_eps_nuc_section(1d3)
      call pgslw(1)
      disp = legend_disp1 + legend_del_disp
      if (show_legend) &
         call pgmtxt('T', disp, legend_coord, legend_fjust, '> 1000 erg g\u-1\d s\u-1\d')

      call pgsci(clr_Crimson)
      call pgslw(20)
      call do_show_eps_nuc_section(1d7)
      call pgslw(1)
      disp = legend_disp1
      if (show_legend) &
         call pgmtxt('T', disp, legend_coord, legend_fjust, '> 10\u7\d erg g\u-1\d s\u-1\d')

      disp = legend_disp1 + 2 * legend_del_disp

      call pgsci(clr_no_mixing)
      call pgslw(10)
      call pgline(nz, xvec, yvec)
      has_convection = do_show_convective_section()
      has_leftover_convection = do_show_leftover_convective_section()
      has_overshoot = do_show_overshoot_section()
      has_semiconvection = do_show_semiconvective_section()
      has_thermohaline = do_show_thermohaline_section()
      if (s% rotation_flag) then
         has_rotation = do_show_rotation_section()
      else
         has_rotation = .false.
      end if
      call pgslw(1)
      if (show_legend) then
         call pgslw(1)
         disp = disp + legend_del_disp
         call show_legend_text(clr_no_mixing, 'no mixing')
         if (s% rotation_flag) then
            disp = disp + legend_del_disp
            call show_legend_text(clr_rotation, 'rotation')
         end if
         disp = disp + legend_del_disp
         call show_legend_text(clr_convection, 'convection')
         disp = disp + legend_del_disp
         call show_legend_text(clr_overshoot, 'overshoot')
         disp = disp + legend_del_disp
         call show_legend_text(clr_semiconvection, 'semiconvection')
         disp = disp + legend_del_disp
         call show_legend_text(clr_thermohaline, 'thermohaline')
         if(s% pg% show_TRho_accretion_mesh_borders) then
            disp = disp + legend_del_disp
            call show_legend_text(clr_RoyalPurple, 'Lagrangian Outer Border')
            disp = disp + legend_del_disp
            call show_legend_text(clr_RoyalBlue, 'Homologous Inner Boundary')
            disp = disp + legend_del_disp
            call show_legend_text(clr_Tan, 'Mass Added This Step')
         end if
         call pgslw(10)
      end if

      if (show_mass_pts) &
         call show_mass_points(s, xvec, yvec, txt_scale, xmin, xmax, ymin, ymax)

      call pgunsa


   contains

      subroutine show_legend_text(clr, txt)
         integer, intent(in) :: clr
         character (len = *), intent(in) :: txt
         call pgsci(clr)
         call pgmtxt('T', disp, legend_coord, legend_fjust, txt)
      end subroutine show_legend_text


      subroutine do_show_eps_nuc_section(eps)
         real(dp), intent(in) :: eps
         integer :: k, first, last
         logical :: inside
         inside = (s% eps_nuc(1) > eps)
         if (inside) first = 1
         do k = 2, s% nz
            if (.not. inside) then
               if (s% eps_nuc(k) > eps) then  ! starting
                  inside = .true.
                  first = k
               end if
            else  ! inside
               if (s% eps_nuc(k) <= eps) then  ! ending
                  last = k - 1
                  call pgline(k - first, xvec(first:last), yvec(first:last))
                  inside = .false.
               end if
            end if
         end do
         if (inside) then
            last = nz
            call pgline(k - first, xvec(first:last), yvec(first:last))
         end if
      end subroutine do_show_eps_nuc_section


      logical function do_show_convective_section()
         do_show_convective_section = do_show_mixing_section(convective_mixing, clr_convection)
      end function do_show_convective_section


      logical function do_show_leftover_convective_section()
         do_show_leftover_convective_section = do_show_mixing_section(leftover_convective_mixing, clr_leftover_convection)
      end function do_show_leftover_convective_section


      logical function do_show_semiconvective_section()
         do_show_semiconvective_section = do_show_mixing_section(semiconvective_mixing, clr_semiconvection)
      end function do_show_semiconvective_section


      logical function do_show_thermohaline_section()
         do_show_thermohaline_section = do_show_mixing_section(thermohaline_mixing, clr_thermohaline)
      end function do_show_thermohaline_section


      logical function do_show_rotation_section()
         do_show_rotation_section = do_show_mixing_section(rotation_mixing, clr_rotation)
      end function do_show_rotation_section


      logical function do_show_overshoot_section()
         do_show_overshoot_section = do_show_mixing_section(overshoot_mixing, clr_overshoot)
      end function do_show_overshoot_section


      logical function do_show_mixing_section(mixing_type, clr)
         integer, intent(in) :: mixing_type, clr
         integer :: k, first, last
         logical :: inside
         include 'formats'
         call pgsave
         call pgsci(clr)
         inside = (s% mixing_type(1) == mixing_type)
         if (inside) first = 1
         do_show_mixing_section = .false.
         do k = 2, s% nz
            if (.not. inside) then
               if (s% mixing_type(k) == mixing_type) then  ! starting
                  inside = .true.
                  first = k
               end if
            else  ! inside
               if (s% mixing_type(k) /= mixing_type) then  ! ending
                  last = k - 1
                  call pgline(k - first, xvec(first:last), yvec(first:last))
                  do_show_mixing_section = .true.
                  inside = .false.
               end if
            end if
         end do
         if (inside) then
            last = nz
            call pgline(k - first, xvec(first:last), yvec(first:last))
            do_show_mixing_section = .true.
         end if
         call pgunsa
      end function do_show_mixing_section


   end subroutine show_profile_line


   subroutine show_mass_points(s, xvec, yvec, txt_scale, xmin, xmax, ymin, ymax)
      type (star_info), pointer :: s
      real, intent(in) :: xvec(:), yvec(:), txt_scale, xmin, xmax, ymin, ymax
      integer :: i
      do i = 1, s% pg% num_profile_mass_points
         call show_mass_point(&
            s, xvec, yvec, txt_scale, xmin, xmax, ymin, ymax, &
            s% pg% profile_mass_point_q(i), &
            s% pg% profile_mass_point_color_index(i), &
            s% pg% profile_mass_point_symbol(i), &
            s% pg% profile_mass_point_symbol_scale(i), &
            s% pg% profile_mass_point_str(i), &
            s% pg% profile_mass_point_str_clr(i), &
            s% pg% profile_mass_point_str_scale(i))
      end do
   end subroutine show_mass_points


   subroutine show_mass_point(&
      s, xvec, yvec, txt_scale, xmin, xmax, ymin, ymax, &
      q_in, clr_index, symbol, symbol_scale, str, str_clr, str_scale)
      type (star_info), pointer :: s
      real, intent(in) :: xvec(:), yvec(:), txt_scale, q_in, &
         xmin, xmax, ymin, ymax, symbol_scale, str_scale
      integer, intent(in) :: clr_index, symbol, str_clr
      character (len = *), intent(in) :: str
      real :: q, q0, q1, x, y, dy
      integer :: nz, i, j, k
      include 'formats'
      q = max(0.0, min(1.0, q_in))
      nz = s% nz
      i = nz
      dy = ymax - ymin
      do k = 1, s% nz - 1
         if (s% q(k) >= q .and. q > s% q(k + 1)) then
            i = k; exit
         end if
      end do
      j = i + 1
      if (j >= nz) j = i
      q0 = s% q(i)
      q1 = s% q(j)
      if ((q0 - q) * (q - q1) < 0) then
         j = i - 1
         q1 = s% q(j)
      end if
      x = find0(xvec(i), q0 - q, xvec(j), q1 - q)
      if (x > xmax .or. x < xmin) return
      y = find0(yvec(i), q0 - q, yvec(j), q1 - q)
      if (y > ymax .or. y < ymin) return
      call pgsave
      call pgscf(1)
      call pgslw(1)
      call pgsci(clr_index)
      call pgsch(symbol_scale * txt_scale)
      call pgpt(1, x, y, symbol)
      call pgsci(str_clr)
      call pgsch(str_scale * txt_scale)
      call pgptxt(x, y - 0.015 * dy, 0.0, 0.0, trim(str))
      call pgunsa
   end subroutine show_mass_point


   real function find0(xx1, yy1, xx2, yy2)
      real :: xx1, yy1, xx2, yy2
      real :: a, b, xz
      ! returns x where y is 0 on line connecting the points (xx1,yy1) and (xx2,yy2)
      a = (xx1 * yy2) - (xx2 * yy1)
      b = yy2 - yy1
      if ((abs(a) >= abs(b) * 1e30) .and. ((yy1 >= 0 .and. yy2 <= 0) &
         .or. (yy1 <= 0 .and. yy2 > 0))) then
         xz = 0.5 * (xx1 + xx2)
      else
         xz = a / b
      end if
      find0 = xz
   end function find0


   subroutine show_box_pgstar(s, str1, str2)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: str1, str2
      real :: ch
      integer :: lw
      call pgqch(ch)
      call pgqlw(lw)
      call pgsch(s% pg% pgstar_num_scale * ch)
      call pgslw(s% pg% pgstar_box_lw)
      call pgbox(str1, 0.0, 0, str2, 0.0, 0)
      call pgsch(ch)
      call pgslw(lw)
   end subroutine show_box_pgstar


   subroutine show_grid_title_pgstar(s, title, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: title
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      if (.not. s% pg% pgstar_grid_show_title) return
      if (len_trim(title) == 0) return
      call pgqch(ch)
      disp = s% pg% pgstar_grid_title_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('T', disp, &
         s% pg% pgstar_grid_title_coord, s% pg% pgstar_grid_title_fjust, title, &
         s% pg% pgstar_grid_title_scale * ch, s% pg% pgstar_grid_title_lw)
   end subroutine show_grid_title_pgstar


   subroutine show_title_pgstar(s, title, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: title
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      if (.not. s% pg% pgstar_show_title) return
      if (len_trim(title) == 0) return
      call pgqch(ch)
      disp = s% pg% pgstar_title_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('T', disp, &
         s% pg% pgstar_title_coord, s% pg% pgstar_title_fjust, title, &
         s% pg% pgstar_title_scale * ch, s% pg% pgstar_title_lw)
   end subroutine show_title_pgstar


   subroutine show_title_label_pgmtxt_pgstar(&
      s, coord, fjust, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: disp
      disp = s% pg% pgstar_title_disp
      if (present(pad)) disp = disp + pad
      call pgmtxt('T', disp, coord, fjust, label)
   end subroutine show_title_label_pgmtxt_pgstar


   subroutine show_xaxis_label_pgstar(s, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      call pgqch(ch)
      disp = s% pg% pgstar_xaxis_label_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('B', disp, 0.5, 0.5, label, &
         s% pg% pgstar_xaxis_label_scale * ch, s% pg% pgstar_xaxis_label_lw)
   end subroutine show_xaxis_label_pgstar


   subroutine show_xaxis_label_pgmtxt_pgstar(&
      s, coord, fjust, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: disp
      disp = s% pg% pgstar_xaxis_label_disp
      if (present(pad)) disp = disp + pad
      call pgmtxt('B', disp, coord, fjust, label)
   end subroutine show_xaxis_label_pgmtxt_pgstar


   subroutine show_left_yaxis_label_pgstar(s, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      call pgqch(ch)
      disp = s% pg% pgstar_left_yaxis_label_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('L', disp, 0.5, 0.5, label, &
         s% pg% pgstar_left_yaxis_label_scale * ch, s% pg% pgstar_left_yaxis_label_lw)
   end subroutine show_left_yaxis_label_pgstar


   subroutine show_right_yaxis_label_pgstar(s, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad
      optional pad
      real :: ch, disp
      call pgqch(ch)
      disp = s% pg% pgstar_right_yaxis_label_disp
      if (present(pad)) disp = disp + pad
      call do1_pgmtxt('R', disp, 0.5, 0.5, label, &
         s% pg% pgstar_right_yaxis_label_scale * ch, s% pg% pgstar_right_yaxis_label_lw)
   end subroutine show_right_yaxis_label_pgstar


   subroutine show_left_yaxis_label_pgmtxt_pgstar(&
      s, coord, fjust, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: ch, disp
      call pgqch(ch)
      call pgsch(1.1 * ch)
      disp = s% pg% pgstar_left_yaxis_label_disp
      if (present(pad)) disp = disp + pad
      call pgmtxt('L', disp, coord, fjust, label)
      call pgsch(ch)
   end subroutine show_left_yaxis_label_pgmtxt_pgstar


   subroutine show_right_yaxis_label_pgmtxt_pgstar(&
      s, coord, fjust, label, pad)
      type (star_info), pointer :: s
      character (len = *), intent(in) :: label
      real, intent(in) :: pad, coord, fjust
      optional pad
      real :: ch, disp
      call pgqch(ch)
      call pgsch(1.1 * ch)
      disp = s% pg% pgstar_right_yaxis_label_disp
      if (present(pad)) disp = disp + pad
      call pgmtxt('R', disp, coord, fjust, label)
      call pgsch(ch)
   end subroutine show_right_yaxis_label_pgmtxt_pgstar


   subroutine show_model_number_pgstar(s)
      type (star_info), pointer :: s
      character (len = 32) :: str
      real :: ch
      if (.not. s% pg% pgstar_show_model_number) return
      write(str, '(i9)') s% model_number
      str = 'model ' // trim(adjustl(str))
      call pgqch(ch)
      call do1_pgmtxt('T', &
         s% pg% pgstar_model_disp, s% pg% pgstar_model_coord, &
         s% pg% pgstar_model_fjust, str, &
         s% pg% pgstar_model_scale * ch, s% pg% pgstar_model_lw)
   end subroutine show_model_number_pgstar


   subroutine show_age_pgstar(s)
      type (star_info), pointer :: s
      character (len = 32) :: age_str, units_str
      real(dp) :: age
      real :: ch
      integer :: len, i, j, iE
      if (.not. s% pg% pgstar_show_age) return
      age = s% star_age
      if (s% pg% pgstar_show_age_in_seconds) then
         age = age * secyer
         units_str = 'secs'
      else if (s% pg% pgstar_show_age_in_minutes) then
         age = age * secyer / 60
         units_str = 'mins'
      else if (s% pg% pgstar_show_age_in_hours) then
         age = age * secyer / (60 * 60)
         units_str = 'hrs'
      else if (s% pg% pgstar_show_age_in_days) then
         age = age * secyer / secday
         units_str = 'days'
      else if (s% pg% pgstar_show_age_in_years) then
         !age = age
         units_str = 'yrs'
      else if (s% pg% pgstar_show_log_age_in_years) then
         age = log10(max(1d-99, age))
         units_str = 'log yrs'
      else if (age * secyer < 60) then
         age = age * secyer
         units_str = 'secs'
      else if (age * secyer < 60 * 60) then
         age = age * secyer / 60
         units_str = 'mins'
      else if (age * secyer < secday) then
         age = age * secyer / (60 * 60)
         units_str = 'hrs'
      else if (age * secyer < secday * 500) then
         age = age * secyer / secday
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
         s% pg% pgstar_age_disp, s% pg% pgstar_age_coord, &
         s% pg% pgstar_age_fjust, age_str, &
         s% pg% pgstar_age_scale * ch, s% pg% pgstar_age_lw)
   end subroutine show_age_pgstar


   logical function read_values_from_file(fname, x_data, y_data, data_len)
      character(len = *), intent(in) :: fname
      real, allocatable, dimension(:) :: x_data, y_data
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


   subroutine show_pgstar_decorator(id, use_flag, pgstar_decorator, plot_num, ierr)
      logical, intent(in) :: use_flag
      real :: xmin, xmax, ymin, ymax
      integer, intent(in) :: id, plot_num
      integer, intent(inout) :: ierr
      procedure(pgstar_decorator_interface), pointer :: pgstar_decorator

      if(use_flag)then
         if(associated(pgstar_decorator))then
            call pgsave
            call PGQWIN(xmin, xmax, ymin, ymax)
            call pgstar_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
            call pgunsa
            if(ierr/=0)then
               write(*, *) "Error in pgstar_decorator"
            end if
         end if
      end if

   end subroutine show_pgstar_decorator

end module pgstar_support

