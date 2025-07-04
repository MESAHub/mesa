! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton & The MESA Team
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

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib

      implicit none

      include "test_suite_extras_def.inc"


      ! these routines are called by the standard run_star check_model
      contains

      include "test_suite_extras.inc"


      subroutine extras_controls(id, ierr)
         use astero_def, only: star_astero_procs
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve

         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         include 'set_star_astero_procs.inc'
      end subroutine extras_controls


      subroutine set_constraint_value(id, name, val, ierr)  ! called from star_astero code
         use astero_def, only: Z_div_X_solar

         integer, intent(in) :: id
         character(len=strlen), intent(in) :: name
         real(dp), intent(out) :: val
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: i
         real(dp) :: X, Y, Z
         ! constraints are predefined in the simplex_search_data.
         ! this routine's job is to assign those variables to current value in the model.
         ! it is called whenever a new value of chi2 is calculated.
         ! only necessary to set the constraints you are actually using.
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         X = max(s% surface_h1, 1d-10)
         Y = s% surface_he3 + s% surface_he4
         Z = max(1d-99, min(1d0, 1-X-Y))

         select case (name)
            ! for custom constraints, create a case with the name of your constraint e.g.
            ! case ('delta_Pg')
            !    val = s% delta_Pg
            ! fall back to history column if user doesn't define name
            case ('M_H')
               val = log10(Z/X/Z_div_X_solar)

            case ('surface_Z_div_X')
               val = Z/X

            case ('surface_He')
               val = Y

            case ('Rcz')
               do i = 1, s% nz-1  ! locate bottom of solar convective zone
                  if (s% mixing_type(i+1) /= convective_mixing &
                        .and. s% mixing_type(i) == convective_mixing) then
                     if (s% r(i+1) > 0.25d0*Rsun .and. s% r(i) < 0.9d0*Rsun) then
                        val = s% r(i)/Rsun
                        exit
                     end if
                  end if
               end do

            case default
               val = star_get_history_output(s, name, ierr)
               if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
         end select

      end subroutine set_constraint_value


      subroutine set_param(id, name, val, ierr)  ! called from star_astero code
         use astero_def, only: f0_ov_div_f_ov, Y_frac_he3, Z_div_X_solar, &
            Y_depends_on_Z, dYdZ, Y0

         integer, intent(in) :: id
         character(len=strlen), intent(in) :: name  ! which of param's will be set
         real(dp), intent(in) :: val
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         real(dp) :: X, Y, Z_div_X, c

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         select case (name)
            case ('initial_mass')
               s% job% new_mass = val
               s% job% relax_initial_mass = .true.
               s% initial_mass = val
            case ('initial_Y')
               s% job% initial_he3 = Y_frac_he3*val
               s% job% initial_he4 = val - s% job% initial_he3
               s% job% set_uniform_initial_composition = .true.
            case ('initial_FeH')
               ! this only works because initial_Y is an earlier parameter than FeH
               ! so set first and we can use that value to infer Z and X
               Z_div_X = Z_div_X_solar*exp10(val)  ! (Z/X) = (Z/X)sun * 10^[Fe/H]

               if (Y_depends_on_Z) then
                  c = 1d0 + Z_div_X*(1d0 + dYdZ)
                  X = (1d0 - Y0)/c
                  Y = (Y0 + Z_div_X*(dYdZ + Y0))/c

                  s% job% initial_he3 = Y_frac_he3*Y
                  s% job% initial_he4 = Y - s% job% initial_he3
               else
                  Y = s% job% initial_he3/Y_frac_he3  ! get Y, which we know has been set
                  X = (1d0 - Y)/(1d0 + Z_div_X)  ! X = (1-Y)/(1 + (Z/X))
               end if

               s% job% initial_h1 = X
               s% job% initial_h2 = 0
               s% job% set_uniform_initial_composition = .true.
            case ('alpha')
               s% mixing_length_alpha = val
            case ('f_ov')
               if (val > 0._dp) then
                  s% overshoot_scheme(1) = 'exponential'
                  s% overshoot_zone_type(1) = 'any'
                  s% overshoot_zone_loc(1) = 'any'
                  s% overshoot_bdy_loc(1) = 'any'
                  s% overshoot_f(1) = val
                  s% overshoot_f0(1) = f0_ov_div_f_ov*val
               else
                  s% overshoot_scheme(1) = ''
                  s% overshoot_zone_type(1) = ''
                  s% overshoot_zone_loc(1) = ''
                  s% overshoot_bdy_loc(1) = ''
                  s% overshoot_f(1) = 0d0
                  s% overshoot_f0(1) = 0d0
               end if
            case ('semiconv')
               s% alpha_semiconvection = val
               write(*,*) 'set semiconvection to', val
            case ('thermohal')
               s% thermohaline_coeff = val
               write(*,*) 'set thermohaline to', val
            case ('diffusion')
               s% diffusion_class_factor(:) = val
               write(*,*) 'set diffusion to', val
            case default
               ierr = -1
               write(*,*) 'invalid name in set_param', name
         end select

      end subroutine set_param


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         use astero_def
         use utils_lib, only: mv

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         character (len=256) :: format_string, num_string, basename
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)

         ! demonstrate how to move some files generated for each
         ! sample
         write(format_string,'( "(i",i2.2,".",i2.2,")" )') num_digits, num_digits
         write(num_string,format_string) sample_number+1  ! sample number hasn't been incremented yet
         basename = trim(sample_results_prefix) // trim(num_string)
         call move(best_model_fgong_filename, '.fgong')
         call move(best_model_gyre_filename, '.gyre')
         call move(best_model_profile_filename, '.profile')
         call move(best_model_save_model_filename, '.mod')

         contains

         subroutine move(filename, extension)
            character (len=*), intent(in) :: filename, extension
            character (len=256) :: src, tgt

            src = trim(astero_results_directory) // '/' // trim(filename)
            tgt = trim(astero_results_directory) // '/' // trim(basename) // extension

            ! file won't exist if sample quit early because of bad chi-squared,
            ! so skip errors
            call mv(trim(src), trim(tgt), skip_errors=.true.)

            write(*,*) 'moved ' // trim(src) // ' to ' // trim(tgt)
         end subroutine move

      end subroutine extras_after_evolve


      integer function extras_check_model(id)
         use astero_lib, only: astero_gyre_is_enabled

         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         extras_check_model = keep_going
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         include 'formats'

         if (.not. astero_gyre_is_enabled) then
            write(*,*) 'not using gyre: pretending iteration worked and exiting'
            write(*,*) 'save_sample_results_to_file outputs/from_file_results.data'
            extras_check_model = terminate
            return
         end if

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step


      end module run_star_extras
