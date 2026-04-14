! ***********************************************************************
!
!   Copyright (C) 2025  The MESA Team
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

module colors_history

   use const_def, only: dp, mbolsun, Teffsun, rsun
   use utils_lib, only: mesa_error
   use colors_def, only: Colors_General_Info, get_colors_ptr, num_color_filters, color_filter_names
   use colors_utils, only: remove_dat, resolve_path
   use bolometric, only: calculate_bolometric
   use synthetic, only: calculate_synthetic
   use mod_colors, only: Eval_Colors, thead_all, num_thead, lgt_list, max_num_bcs_per_file

   implicit none

contains

   integer function how_many_colors_history_columns(colors_handle)
      integer, intent(in) :: colors_handle
      integer :: num_cols, ierr
      type(Colors_General_Info), pointer :: colors_settings

      ierr = 0
      call get_colors_ptr(colors_handle, colors_settings, ierr)
      if (ierr /= 0) then
         write (*, *) 'failed in colors_ptr'
         how_many_colors_history_columns = 0
         return
      end if

      if (.not. colors_settings%use_colors) then
         num_cols = 0
      else
         num_cols = 3 + num_color_filters
      end if

      how_many_colors_history_columns = num_cols
   end function how_many_colors_history_columns

   subroutine data_for_colors_history_columns( &
      t_eff, log_g, R, metallicity, model_number, &
      colors_handle, n, names, vals, ierr)
      use const_def, only: mbolsun, Teffsun, rsun
      use mod_colors, only: Eval_Colors, thead_all, num_thead, lgt_list, max_num_bcs_per_file
      real(dp), intent(in) :: t_eff, log_g, R, metallicity
      integer, intent(in) :: colors_handle, n
      character(len=80) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      integer, intent(in) :: model_number

      type(Colors_General_Info), pointer :: cs
      integer :: i, j, k, filter_offset
      real(dp) :: d, bolometric_magnitude, bolometric_flux, interpolation_radius
      real(dp) :: zero_point, lum
      character(len=256) :: sed_filepath
      character(len=80) :: filter_name
      logical :: make_sed
      real(dp), dimension(max_num_bcs_per_file) :: bcs
      type(lgt_list), pointer :: thead

      real(dp), dimension(:), allocatable :: wavelengths, fluxes

      ierr = 0
      call get_colors_ptr(colors_handle, cs, ierr)
      if (ierr /= 0) then
         write (*, *) 'failed in get_colors_ptr'
         return
      end if

      !write(0,*) 'DEBUG colors_history: n=', n, 'num_color_filters=', num_color_filters, &
      !           'photo_precompute=', cs%photo_precompute
      !flush(0)

      if (cs%photo_precompute) then

         lum = (R / rsun)**2 * (t_eff / Teffsun)**4
         bolometric_magnitude = mbolsun - 2.5_dp * log10(lum)

         names(1) = 'Mag_bol'
         vals(1) = bolometric_magnitude
         names(2) = 'Flux_bol'
         vals(2) = -1.0_dp
         names(3) = 'Interp_rad'
         vals(3) = -1.0_dp
         filter_offset = 3

         k = 0
         do j = 1, num_thead
            thead => thead_all(j)%thead
            call Eval_Colors(log10(t_eff), log_g, metallicity, bcs, thead, thead_all(j)%n_colors, ierr)
            do i = 1, thead_all(j)%n_colors
               k = k + 1
               filter_name = trim(remove_dat(color_filter_names(k)))
               names(k + filter_offset) = filter_name
               if (ierr /= 0) then
                  vals(k + filter_offset) = -1.0_dp
               else
                  vals(k + filter_offset) = bolometric_magnitude - bcs(i)
               end if
            end do
         end do

      else

         ! verify data was loaded at initialization
         if (.not. cs%lookup_loaded) then
            write (*, *) 'colors error: lookup table not loaded'
            ierr = -1
            return
         end if
         if (.not. cs%filters_loaded) then
            write (*, *) 'colors error: filter data not loaded'
            ierr = -1
            return
         end if

         d = cs%distance
         sed_filepath = trim(resolve_path(cs%stellar_atm))
         make_sed = cs%make_csv

         call calculate_bolometric(cs, t_eff, log_g, metallicity, R, d, &
                                   bolometric_magnitude, bolometric_flux, wavelengths, fluxes, &
                                   sed_filepath, interpolation_radius)

         names(1) = "Mag_bol"
         vals(1) = bolometric_magnitude
         names(2) = "Flux_bol"
         vals(2) = bolometric_flux
         names(3) = "Interp_rad"
         vals(3) = interpolation_radius
         filter_offset = 3

         if (n == num_color_filters + filter_offset) then
            do i = 1, num_color_filters
               filter_name = trim(remove_dat(color_filter_names(i)))
               names(i + filter_offset) = filter_name

               ! Negative [M/H] values are valid for metal-poor atmosphere grids.
               ! so we do not apply a limit on the "metallicity" parameter.
               if (t_eff >= 0) then
                  ! Select precomputed zero-point based on magnitude system
                  select case (trim(cs%mag_system))
                  case ('VEGA', 'Vega', 'vega')
                     zero_point = cs%filters(i)%vega_zero_point
                  case ('AB', 'ab')
                     zero_point = cs%filters(i)%ab_zero_point
                  case ('ST', 'st')
                     zero_point = cs%filters(i)%st_zero_point
                  case default
                     write (*, *) 'colors error: unknown magnitude system: ', trim(cs%mag_system)
                     zero_point = -1.0_dp
                  end select

                  vals(i + filter_offset) = calculate_synthetic(t_eff, log_g, metallicity, ierr, &
                                                                wavelengths, fluxes, &
                                                                cs%filters(i)%wavelengths, &
                                                                cs%filters(i)%transmission, &
                                                                zero_point, &
                                                                color_filter_names(i), &
                                                                make_sed, cs%sed_per_model, &
                                                                cs%colors_results_directory, model_number)

                  if (ierr /= 0) vals(i + filter_offset) = -1.0_dp
               else
                  vals(i + filter_offset) = -1.0_dp
                  ierr = 1
               end if
            end do
         else
            ierr = 1
            call mesa_error(__FILE__, __LINE__, 'colors: data_for_colors_history_columns array size mismatch')
         end if

         ! clean up
         if (allocated(wavelengths)) deallocate (wavelengths)
         if (allocated(fluxes)) deallocate (fluxes)

      end if

   end subroutine data_for_colors_history_columns

end module colors_history
