! ***********************************************************************
!
!   Copyright (C) 2010-2022  Bill Paxton, Pablo Marchant, Matthias Fabry
!                                                          & The MESA Team
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

module binary_history

   use const_def, only: dp, pi, lsun, msun, rsun, secyer, secday, two_thirds
   use math_lib
   use binary_def
   use binary_private_def
   use binary_history_specs

   implicit none

contains

   integer function how_many_binary_history_columns(binary_id)
      integer, intent(in) :: binary_id
      integer :: numcols, ierr
      type (binary_info), pointer :: b

      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         numcols = 0
         return
      end if

      if (.not. associated(b% binary_history_column_spec)) then
         numcols = 0
      else
         numcols = size(b% binary_history_column_spec, dim = 1)
      end if

      how_many_binary_history_columns = numcols
   end function how_many_binary_history_columns


   subroutine data_for_binary_history_columns(&
      binary_id, n, names, vals, ierr)
      integer, intent(in) :: binary_id, n
      character (len = 80) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr

      type (binary_info), pointer :: b
      integer :: c, int_val, j
      logical :: is_int_val
      real(dp) :: val

      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      do j = 1, n
         c = b% binary_history_column_spec(j)
         names(j) = trim(binary_history_column_name(c))
         call binary_history_getval(&
            b, c, val, int_val, is_int_val, ierr)
         if (ierr /= 0) then
            write(*, *) "Unknown binary_history_columns.list column"
            return
         end if
         if (is_int_val) then
            vals(j) = int_val
         else
            vals(j) = val
         end if
      end do
   end subroutine data_for_binary_history_columns

   subroutine write_binary_history_info(b, ierr)
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      logical, parameter :: write_flag = .true.

      call do_binary_history_info(b, write_flag, ierr)
   end subroutine write_binary_history_info

   subroutine do_get_data_for_binary_history_columns(b, ierr)
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr
      logical, parameter :: write_flag = .false.

      call do_binary_history_info(b, write_flag, ierr)
   end subroutine do_get_data_for_binary_history_columns

   subroutine do_binary_history_info(b, write_flag, ierr)
      use utils_lib, only : integer_dict_create_hash, integer_dict_free
      type (binary_info), pointer :: b
      logical, intent(in) :: write_flag
      integer, intent(out) :: ierr

      character (len = strlen) :: fname, dbl_fmt, int_fmt, txt_fmt
      integer :: numcols, io, i, col, j, i0, n

      integer :: num_extra_header_items, num_extra_cols

      character (len = maxlen_history_column_name), pointer, dimension(:) :: &
         extra_header_item_names, extra_col_names
      real(dp), pointer, dimension(:) :: &
         extra_header_item_vals, extra_col_vals

      logical :: binary_history_file_exists
      character (len = maxlen_history_column_name), pointer :: &
         names(:)  ! (num_history_columns)
      real(dp), pointer :: vals(:)  ! (num_history_columns)
      logical, pointer :: is_int(:)  ! (num_history_columns)

      include 'formats'

      extra_header_item_names => null()
      extra_header_item_vals => null()

      extra_col_names => null()
      extra_col_vals => null()

      dbl_fmt = b% history_dbl_format
      int_fmt = b% history_int_format
      txt_fmt = b% history_txt_format

      ierr = 0

      if (.not. associated(b% binary_history_column_spec)) then
         numcols = 0
      else
         numcols = size(b% binary_history_column_spec, dim = 1)
      end if

      num_extra_cols = b% how_many_extra_binary_history_columns(b% binary_id)
      n = numcols + num_extra_cols

      if (n == 0) then
         write(*, *) 'WARNING: do not have any output specified for binary logs.'
         return
      end if
!      write(*, *) " got num of cols"
      if (b% number_of_binary_history_columns < 0) then
         b% number_of_binary_history_columns = n
      else if (b% number_of_binary_history_columns /= n) then
         if (associated(b% binary_history_values)) then
            deallocate(b% binary_history_values)
            nullify(b% binary_history_values)
         end if
         if (associated(b% binary_history_names)) then
            deallocate(b% binary_history_names)
            nullify(b% binary_history_names)
         end if
         if (associated(b% binary_history_value_is_integer)) then
            deallocate(b% binary_history_value_is_integer)
            nullify(b% binary_history_value_is_integer)
         end if
         if (associated(b% binary_history_names_dict)) then
            call integer_dict_free(b% binary_history_names_dict)
            nullify(b% binary_history_names_dict)
         end if
         b% need_to_set_binary_history_names_etc = .true.
         b% number_of_binary_history_columns = n
      end if

      if (.not. associated(b% binary_history_values)) then
         allocate(b% binary_history_values(n))
      else if (size(b% binary_history_values, dim = 1) /= n) then
         ierr = -1
         write(*, 3) 'bad size b% binary_history_values', &
            size(b% binary_history_values, dim = 1), n
      end if
      vals => b% binary_history_values

      if (.not. associated(b% binary_history_names)) then
         allocate(b% binary_history_names(n))
      else if (size(b% binary_history_names, dim = 1) /= n) then
         ierr = -1
         write(*, 3) 'bad size b% binary_history_names', &
            size(b% binary_history_names, dim = 1), n
      end if

      if (b% need_to_set_binary_history_names_etc) then
         names => b% binary_history_names
      else
         nullify(names)
      end if

      if (.not. associated(b% binary_history_value_is_integer)) then
         allocate(b% binary_history_value_is_integer(n))
      else if (size(b% binary_history_value_is_integer, dim = 1) /= n) then
         ierr = -1
         write(*, 2) 'bad size b% binary_history_value_is_integer', &
            size(b% binary_history_value_is_integer, dim = 1), n
      end if
      if (b% need_to_set_binary_history_names_etc) then
         is_int => b% binary_history_value_is_integer
      else
         nullify(is_int)
      end if
!      write(*, *) " associated arrays"

      nullify(extra_col_names)
      nullify(extra_col_vals)

      if (num_extra_cols > 0) then
         allocate(&
            extra_col_names(num_extra_cols), extra_col_vals(num_extra_cols), stat = ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if
         extra_col_names(1:num_extra_cols) = 'unknown'
         extra_col_vals(1:num_extra_cols) = -1d99
         call b% data_for_extra_binary_history_columns(&
            b% binary_id, num_extra_cols, extra_col_names, extra_col_vals, ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if
         do i = 1, num_extra_cols
            if(trim(extra_col_names(i))=='unknown') then
               write(*, *) "Warning empty history name for extra_binary_history_column ", i
            end if
         end do
      end if

!      write(*, *) " starting write loop ", write_flag
      i0 = 1
      if (write_flag .and. (open_close_log .or. b% s_donor% model_number == -100)) then
         fname = trim(b% log_directory) // '/' // trim(b% history_name)
         inquire(file = trim(fname), exist = binary_history_file_exists)
         if ((.not. binary_history_file_exists) .or. b% open_new_history_file) then
            ierr = 0
            open(newunit = io, file = trim(fname), action = 'write', iostat = ierr)
            b% open_new_history_file = .false.
         else
            i0 = 3
            open(newunit = io, file = trim(fname), action = 'write', position = 'append', iostat = ierr)
         end if
         if (ierr /= 0) then
            write(*, *) 'failed to open ' // trim(fname)
            call dealloc
            return
         end if
      end if

      if (write_flag .and. i0 == 1) then  ! write parameters at start of log

         num_extra_header_items = b% how_many_extra_binary_history_header_items(b% binary_id)

         if (num_extra_header_items > 0) then
            allocate(&
               extra_header_item_names(num_extra_header_items), &
               extra_header_item_vals(num_extra_header_items), stat = ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            extra_header_item_names(1:num_extra_header_items) = 'unknown'
            extra_header_item_vals(1:num_extra_header_items) = -1d99
            call b% data_for_extra_binary_history_header_items(&
               b% binary_id, num_extra_header_items, &
               extra_header_item_names, extra_header_item_vals, ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            do i = 1, num_extra_header_items
               if(trim(extra_header_item_names(i))=='unknown') then
                  write(*, *) "Warning empty history name for extra_binary_history_header ", i
               end if
            end do
         end if

         do i = 1, 3
            col = 0
            call write_string(io, col, i, 'version_number', version_number)
            call write_val(io, col, i, 'initial_don_mass', initial_mass(1))
            call write_val(io, col, i, 'initial_acc_mass', initial_mass(2))
            call write_val(io, col, i, 'initial_period_days', &
               initial_binary_period / (3600 * 24))

            call write_string(io, col, i, 'compiler', compiler_name)
            call write_string(io, col, i, 'build', compiler_version_name)
            call write_string(io, col, i, 'MESA_SDK_version', mesasdk_version_name)
            call write_string(io, col, i, 'date', date)

            do j = 1, num_extra_header_items
               call write_val(io, col, i, &
                  extra_header_item_names(j), extra_header_item_vals(j))
            end do

            write(io, *)
         end do
         write(io, *)
      end if

      do i = i0, 3  ! add a row to the log
         col = 0
!         write(*, *) "doing cols pass", i
         do j = 1, numcols
            call do_col(i, j)
         end do
         do j = 1, num_extra_cols
            call do_extra_col(i, j)
         end do
         if (write_flag) write(io, *)
      end do
!      write (*, *) "cols handled"
      if (open_close_log) close(io)

      call dealloc
!      write(*, *) "history written"

      b% model_number_of_binary_history_values = b% model_number

      if (b% need_to_set_binary_history_names_etc) then
!         write(*, *) " creating hash of the history dict"
         call integer_dict_create_hash(b% binary_history_names_dict, ierr)
         if (ierr /= 0) then
            write(*, *) "hash failed"
            return
         end if
      end if

      b% need_to_set_binary_history_names_etc = .false.

   contains


      subroutine dealloc
         if (associated(extra_header_item_names)) deallocate(extra_header_item_names)
         if (associated(extra_header_item_vals)) deallocate(extra_header_item_vals)
         if (associated(extra_col_names)) deallocate(extra_col_names)
         if (associated(extra_col_vals)) deallocate(extra_col_vals)
      end subroutine dealloc


      subroutine do_extra_col(pass, j)
         integer, intent(in) :: pass, j
         if (pass == 1) then
            if (write_flag) write(io, fmt = int_fmt, advance = 'no') j + numcols
         else if (pass == 2) then
            call do_name(j + numcols, extra_col_names(j))
         else if (pass == 3) then
            call do_val(j + numcols, extra_col_vals(j))
         end if
      end subroutine do_extra_col


      subroutine do_name(j, col_name)
         use utils_lib, only: integer_dict_define
         integer, intent(in) :: j
         character (len = *), intent(in) :: col_name
         if (write_flag) write(io, fmt = txt_fmt, advance = 'no') trim(col_name)
         if (associated(names)) names(j) = trim(col_name)
         if (b% need_to_set_binary_history_names_etc) then
            call integer_dict_define(b% binary_history_names_dict, col_name, j, ierr)
            if (ierr /= 0) write(*, *) 'failed in dict define ' // trim(col_name)
         end if

      end subroutine do_name


      subroutine do_col(pass, j)
         integer, intent(in) :: pass, j
         if (pass == 1) then
            call do_col_pass1
         else if (pass == 2) then
            call do_col_pass2(j)
         else if (pass == 3) then
            call do_col_pass3(b% binary_history_column_spec(j))
         end if
      end subroutine do_col


      subroutine do_col_pass1  ! write the column number
         col = col + 1
         if (write_flag) write(io, fmt = int_fmt, advance = 'no') col
      end subroutine do_col_pass1


      subroutine do_col_pass2(j)  ! get the column name
         integer, intent(in) :: j
         character (len = 100) :: col_name
         integer :: c
         c = b% binary_history_column_spec(j)
         col_name = trim(binary_history_column_name(c))
         call do_name(j, col_name)
      end subroutine do_col_pass2


      subroutine do_col_pass3(c)  ! get the column value
         integer, intent(in) :: c
         integer :: k, int_val
         logical :: is_int_val
         real(dp) :: val
         int_val = 0; val = 0; is_int_val = .false.
         call binary_history_getval(&
            b, c, val, int_val, is_int_val, ierr)
         if (ierr /= 0) then
            write(*, *) 'missing log info for ' // trim(binary_history_column_name(c)), j, k
            return
         end if
         if (is_int_val) then
            call do_int_val(j, int_val)
         else
            call do_val(j, val)
         end if
      end subroutine do_col_pass3


      subroutine do_val(j, val)
         use utils_lib, only : is_bad
         integer, intent(in) :: j
         real(dp), intent(in) :: val
         if (write_flag) then
            if (is_bad(val)) then
               write(io, fmt = dbl_fmt, advance = 'no') -1d99
            else
               write(io, fmt = dbl_fmt, advance = 'no') val
            end if
         end if
         if (associated(vals)) vals(j) = val
         if (associated(is_int)) is_int(j) = .false.
      end subroutine do_val


      subroutine do_int_val(j, val)
         integer, intent(in) :: j
         integer, intent(in) :: val
         if (write_flag) write(io, fmt = int_fmt, advance = 'no') val
         if (associated(vals)) vals(j) = dble(val)
         if (associated(is_int)) is_int(j) = .true.
      end subroutine do_int_val


      subroutine write_integer(io, col, pass, name, val)
         integer, intent(in) :: io, pass
         integer, intent(inout) :: col
         character (len = *), intent(in) :: name
         integer, intent(in) :: val
         if (pass == 1) then
            col = col + 1
            write(io, fmt = int_fmt, advance = 'no') col
         else if (pass == 2) then
            write(io, fmt = txt_fmt, advance = 'no') trim(name)
         else if (pass == 3) then
            write(io, fmt = int_fmt, advance = 'no') val
         end if
      end subroutine write_integer


      subroutine write_val(io, col, pass, name, val)  ! for header items only
         integer, intent(in) :: io, pass
         integer, intent(inout) :: col
         character (len = *), intent(in) :: name
         real(dp), intent(in) :: val
         if (pass == 1) then
            col = col + 1
            write(io, fmt = int_fmt, advance = 'no') col
         else if (pass == 2) then
            write(io, fmt = txt_fmt, advance = 'no') trim(name)
         else if (pass == 3) then
            write(io, fmt = dbl_fmt, advance = 'no') val
         end if
      end subroutine write_val


      subroutine write_string(io, col, pass, name, val)  !for header items only
         integer, intent(in) :: io, pass
         integer, intent(inout) :: col
         character(len = *), intent(in) :: name, val
         character(len = strlen) :: my_val

         my_val = '"' // trim(val) // '"'
         if (pass == 1) then
            col = col + 1
            write(io, fmt = int_fmt, advance = 'no') col
         else if (pass == 2) then
            write(io, fmt = txt_fmt, advance = 'no') trim(name)
         else if (pass == 3) then
            write(io, fmt = txt_fmt, advance = 'no') trim(my_val)
         end if
      end subroutine write_string


   end subroutine do_binary_history_info


   subroutine binary_history_getval(b, c, val, int_val, is_int_val, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: c
      real(dp), intent(out) :: val
      integer, intent(out) :: int_val
      logical, intent(out) :: is_int_val
      integer, intent(out) :: ierr

      include 'formats'

      ierr = 0
      is_int_val = .false.
      int_val = 0
      val = 0
      select case(c)

      case(bh_model_number)
         int_val = b% model_number
         is_int_val = .true.
      case(bh_age)
         val = b% binary_age
      case(bh_donor_index)
         int_val = b% d_i
         is_int_val = .true.
      case(bh_period_days)
         val = b% period / secday
      case(bh_period_hr)
         val = b% period / (60d0 * 60d0)
      case(bh_period_minutes)
         val = b% period / 60d0
      case(bh_lg_separation)
         val = safe_log10(b% separation)
      case(bh_binary_separation)
         val = b% separation / Rsun
      case(bh_eccentricity)
         val = b% eccentricity
      case(bh_star_1_radius)
         val = b% r(1) / Rsun
      case(bh_star_2_radius)
         val = b% r(2) / Rsun
      case(bh_rl_1)
         val = b% rl(1) / Rsun
      case(bh_rl_2)
         val = b% rl(2) / Rsun
      case(bh_rl_overflow_1)
         val = (b% r(1) - b% rl(1)) / Rsun
      case(bh_rl_overflow_2)
         val = (b% r(2) - b% rl(2)) / Rsun
      case(bh_rl_relative_overflow_1)
         val = b% rl_relative_gap(1)
      case(bh_rl_relative_overflow_2)
         val = b% rl_relative_gap(2)
      case(bh_P_rot_div_P_orb_1)
         if (b% point_mass_i /= 1) then
            val = 2 * pi / b% s1% omega_avg_surf / b% period
         else
            val = 0.0d0
         end if
      case(bh_P_rot_div_P_orb_2)
         if (b% point_mass_i /= 2) then
            val = 2 * pi / b% s2% omega_avg_surf / b% period
         else
            if (.not. b% model_twins_flag) then
               val = 0.0d0
            else
               val = 2 * pi / b% s1% omega_avg_surf / b% period
            end if
         end if
      case(bh_lg_t_sync_1)
         val = safe_log10(abs(b% t_sync_1) / secyer)
      case(bh_lg_t_sync_2)
         val = safe_log10(abs(b% t_sync_2) / secyer)
      case(bh_star_1_mass)
         val = b% m(1) / Msun
      case(bh_lg_star_1_mass)
         val = safe_log10(b% m(1) / Msun)
      case(bh_star_2_mass)
         val = b% m(2) / Msun
      case(bh_lg_star_2_mass)
         val = safe_log10(b% m(2) / Msun)
      case(bh_sum_of_masses)
         val = (b% m(1) + b% m(2)) / Msun
      case(bh_mass_ratio)
         val = b% m(2) / b% m(1)
      case(bh_obs_mass_ratio)
         val = min(b% m(2) / b% m(1), b% m(1) / b% m(2))
      case(bh_lg_mtransfer_rate)
         val = safe_log10(abs(b% step_mtransfer_rate) / Msun * secyer)
      case(bh_lg_mstar_dot_1)
         val = safe_log10(abs(b% component_mdot(1)) / Msun * secyer)
      case(bh_lg_mstar_dot_2)
         val = safe_log10(abs(b% component_mdot(2)) / Msun * secyer)
      case(bh_lg_system_mdot_1)
         val = safe_log10(abs(b% mdot_system_transfer(1)) / Msun * secyer)
      case(bh_lg_system_mdot_2)
         val = safe_log10(abs(b% mdot_system_transfer(2)) / Msun * secyer)
      case(bh_lg_wind_mdot_1)
         val = safe_log10(abs(b% mdot_system_wind(1)) / Msun * secyer)
      case(bh_lg_wind_mdot_2)
         val = safe_log10(abs(b% mdot_system_wind(2)) / Msun * secyer)
      case(bh_star_1_div_star_2_mass)
         val = b% m(1) / b% m(2)
      case(bh_delta_star_1_mass)
         val = b% m(1) - initial_mass(1)
      case(bh_delta_star_2_mass)
         val = b% m(2) - initial_mass(2)
      case(bh_lg_F_irr)
         val = safe_log10(b% s_donor% irradiation_flux)
      case(bh_fixed_xfer_fraction)
         val = b% fixed_xfer_fraction
      case(bh_eff_xfer_fraction)
         if (b% component_mdot(b% d_i) == 0d0) then
            val = 1d0
         else
            val = (-b% component_mdot(b% a_i)) / (b% component_mdot(b% d_i))
         end if
      case(bh_lg_mdot_edd)
         if (b% limit_retention_by_mdot_edd) then
            val = safe_log10(b% mdot_edd / Msun * secyer)
         else
            val = safe_log10(0d0)
         end if
      case(bh_mdot_edd_eta)
         if (b% limit_retention_by_mdot_edd) then
            val = b% mdot_edd_eta
         else
            val = 0d0
         end if
      case(bh_lg_accretion_luminosity)
         val = safe_log10(b% accretion_luminosity / Lsun)
      case(bh_bh_spin)
         if (b% point_mass_i /= 0) then
            val = sqrt(two_thirds) &
               * (b% eq_initial_bh_mass / min(b% m(b% point_mass_i), sqrt(6d0) * b% eq_initial_bh_mass)) &
               * (4d0 - sqrt(18d0 * pow2(b% eq_initial_bh_mass / &
                  min(b% m(b% point_mass_i), sqrt(6d0) * b% eq_initial_bh_mass)) - 2d0))
         else
            val = 0
         end if
      case(bh_v_orb_1)
         val = 2.0d0 * pi * b% m(2) / (b% m(1) + b% m(2)) * b% separation / b% period / 1.0d5
      case(bh_v_orb_2)
         val = 2.0d0 * pi * b% m(1) / (b% m(1) + b% m(2)) * b% separation / b% period / 1.0d5
      case(bh_J_orb)
         val = b% angular_momentum_j
      case(bh_J_spin_1)
         if (b% point_mass_i /= 1) then
            val = b% s1% total_angular_momentum
         else
            val = 0d0
         end if
      case(bh_J_spin_2)
         if (b% point_mass_i /= 2) then
            val = b% s2% total_angular_momentum
         else
            if (.not. b% model_twins_flag) then
               val = 0d0
            else
               val = b% s1% total_angular_momentum
            end if
         end if
      case(bh_J_total)
         val = b% angular_momentum_j
         if (b% point_mass_i /= 1) &
            val = val + b% s1% total_angular_momentum
         if (b% point_mass_i /= 2) then
            val = val + b% s2% total_angular_momentum
         else if (b% model_twins_flag) then
            val = val + b% s1% total_angular_momentum
         end if
         val = val
      case(bh_Jdot)
         val = b% jdot
      case(bh_jdot_mb)
         val = b% jdot_mb
      case(bh_jdot_gr)
         val = b% jdot_gr
      case(bh_jdot_ml)
         val = b% jdot_ml
      case(bh_jdot_ls)
         val = b% jdot_ls
      case(bh_jdot_missing_wind)
         val = b% jdot_missing_wind
      case(bh_extra_jdot)
         val = b% extra_jdot
      case(bh_accretion_mode)
         int_val = b% accretion_mode
         is_int_val = .true.
      case(bh_acc_am_div_kep_am)
         val = b% acc_am_div_kep_am
      case(bh_edot)
         val = b% edot
      case(bh_edot_tidal)
         val = b% edot_tidal
      case(bh_edot_enhance)
         val = b% edot_enhance
      case(bh_extra_edot)
         val = b% extra_edot
      case(bh_point_mass_index)
         is_int_val = .true.
         int_val = b% point_mass_i
      case(bh_ignore_rlof_flag)
         is_int_val = .true.
         if (b% ignore_rlof_flag) then
            int_val = 1d0
         else
            int_val = 0d0
         end if
      case(bh_model_twins_flag)
         is_int_val = .true.
         if (b% model_twins_flag) then
            int_val = 1d0
         else
            int_val = 0d0
         end if
      case(bh_CE_flag)
         is_int_val = .true.
         if (b% CE_flag) then
            int_val = 1d0
         else
            int_val = 0d0
         end if
      case(bh_CE_lambda1)
         val = b% CE_lambda1
      case(bh_CE_lambda2)
         val = b% CE_lambda2
      case(bh_CE_Ebind1)
         val = b% CE_Ebind1
      case(bh_CE_Ebind2)
         val = b% CE_Ebind2
      case(bh_CE_num1)
         is_int_val = .true.
         int_val = b% CE_num1
      case(bh_CE_num2)
         is_int_val = .true.
         int_val = b% CE_num2

      case default
         ierr = -1

      end select

   end subroutine binary_history_getval

   subroutine get_binary_history_specs(b, num, names, specs, report)

      use utils_lib
      use utils_def

      type (binary_info), pointer :: b
      integer, intent(in) :: num
      character (len = *), intent(in) :: names(:)
      integer, intent(out) :: specs(:)
      logical, intent(in) :: report

      integer :: i, ierr, n, j, iounit, t
      character (len = strlen) :: buffer, string

      include 'formats'
      ierr = 0
      if (num <= 0) return
      iounit = -1
      specs(1:num) = 0
      do i = 1, num
         buffer = names(i)
         n = len_trim(buffer) + 1
         buffer(n:n) = ' '
         j = 0
         t = token(iounit, n, j, buffer, string)
         if (t /= name_token) then
            if (len_trim(names(i)) > 0 .and. report) &
               write(*, *) 'bad value for name of history item ' // trim(names(i))
            specs(i) = -1
            ierr = 0
            cycle
         end if
         specs(i) = do1_binary_history_spec(&
            iounit, t, n, j, string, buffer, report, ierr)
         if (ierr /= 0) then
            if (report) write(*, *) 'get_binary_history_specs failed for ' // trim(names(i))
            specs(i) = -1
            ierr = 0
         end if
      end do

   end subroutine get_binary_history_specs


   subroutine get_binary_history_values(b, num, specs, &
      is_int_value, int_values, values, failed_to_find_value)
      ! note: this doesn't handle user-defined extra columns

      use utils_lib
      use utils_def

      type (binary_info), pointer :: b
      integer, intent(in) :: num
      integer, intent(in) :: specs(:)
      logical, intent(out) :: is_int_value(:)
      integer, intent(out) :: int_values(:)
      real(dp), intent(inout) :: values(:)
      logical, intent(out) :: failed_to_find_value(:)

      integer :: i, c, ierr

      include 'formats'
      ierr = 0

      do i = 1, num
         failed_to_find_value(i) = .false.
         c = specs(i)
         if (c <= 0) then
            failed_to_find_value(i) = .true.
         else
            call binary_history_getval(&
               b, c, values(i), int_values(i), is_int_value(i), ierr)
            if (ierr /= 0) then
               failed_to_find_value(i) = .true.
               ierr = 0
            end if
         end if
      end do

   end subroutine get_binary_history_values

   logical function get1_binary_hist_value(b, name, val)
      ! includes other_history_columns from run_binary_extras
      use utils_lib, only : integer_dict_lookup
      type (binary_info), pointer :: b
      character (len = *) :: name
      real(dp), intent(out) :: val
      integer :: i, ierr, num_extra_cols
      character (len = 80), pointer, dimension(:) :: extra_col_names
      real(dp), pointer, dimension(:) :: extra_col_vals
      include 'formats'

      get1_binary_hist_value = .false.
      call integer_dict_lookup(b% binary_history_names_dict, name, i, ierr)
      if (ierr /= 0 .or. i <= 0) return  ! didn't find it
      if (associated(b% pg% pgbinary_hist)) then
         if (associated(b% pg% pgbinary_hist% vals)) then
            if (size(b% pg% pgbinary_hist% vals, dim = 1) >= i) then
               val = b% pg% pgbinary_hist% vals(i)
               get1_binary_hist_value = .true.
               return
            end if
         end if
      end if

      ! try extras
      if (associated(b% how_many_extra_binary_history_columns) .and. &
         associated(b% data_for_extra_binary_history_columns)) then
         num_extra_cols = b% how_many_extra_binary_history_columns(b% binary_id)
         if (num_extra_cols > 0) then
            allocate(&
               extra_col_names(num_extra_cols), &
               extra_col_vals(num_extra_cols), stat = ierr)
            call b% data_for_extra_binary_history_columns(&
               b% binary_id, num_extra_cols, extra_col_names, extra_col_vals, ierr)
            do i = 1, num_extra_cols
               if (extra_col_names(i) == name) then
                  val = extra_col_vals(i)
                  get1_binary_hist_value = .true.
                  exit
               end if
            end do
            deallocate(extra_col_names, extra_col_vals)
            if (get1_binary_hist_value) return
         end if
      end if

   end function get1_binary_hist_value

end module binary_history
