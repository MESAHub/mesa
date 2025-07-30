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

   use const_def, only: dp, mesa_dir
   use utils_lib, only: mesa_error
   use colors_def, only: Colors_General_Info, get_colors_ptr, num_color_filters, color_filter_names
   use colors_utils, only: remove_dat
   use bolometric, only: calculate_bolometric
   use synthetic, only: calculate_synthetic

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
         num_cols = 0
         return
      end if

      if (.not. colors_settings%use_colors) then
         num_cols = 0
      else
         num_cols = 2 + num_color_filters
      end if

      how_many_colors_history_columns = num_cols
   end function how_many_colors_history_columns

   subroutine data_for_colors_history_columns( &
      t_eff, log_g, R, metallicity, &
      colors_handle, n, names, vals, ierr)
      real(dp), intent(in) :: t_eff, log_g, R, metallicity
      integer, intent(in) :: colors_handle, n
      character(len=80) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr

      type(Colors_General_Info), pointer :: colors_settings
      integer :: i, filter_offset
      real(dp) :: d, bolometric_magnitude, bolometric_flux
      character(len=256) :: sed_filepath, filter_filepath, filter_name, filter_dir, vega_filepath
      logical :: make_sed

      real(dp), dimension(:), allocatable :: wavelengths, fluxes, filter_wavelengths, filter_trans

      ierr = 0
      call get_colors_ptr(colors_handle, colors_settings, ierr)
      if (ierr /= 0) then
         write (*, *) 'failed in get_colors_ptr'
         return
      end if

      !metallicity = colors_settings% metallicity
      d = colors_settings%distance
      sed_filepath = trim(mesa_dir)//colors_settings%stellar_atm
      filter_dir = trim(mesa_dir)//colors_settings%instrument
      vega_filepath = trim(mesa_dir)//colors_settings%vega_sed
      make_sed = colors_settings%make_csv

      call calculate_bolometric(t_eff, log_g, metallicity, R, d, &
                                bolometric_magnitude, bolometric_flux, wavelengths, fluxes, sed_filepath)
      names(1) = "Mag_bol"
      vals(1) = bolometric_magnitude
      names(2) = "Flux_bol"
      vals(2) = bolometric_flux
      filter_offset = 2

      if (n == num_color_filters + filter_offset) then
         do i = 1, num_color_filters
            filter_name = trim(remove_dat(color_filter_names(i)))
            names(i + filter_offset) = filter_name

            filter_filepath = trim(filter_dir)//"/"//color_filter_names(i)

            if (t_eff >= 0 .and. metallicity >= 0) then
               vals(i + filter_offset) = calculate_synthetic(t_eff, log_g, metallicity, ierr, &
                                                             wavelengths, fluxes, filter_wavelengths, filter_trans, &
                                                             filter_filepath, vega_filepath, color_filter_names(i), &
                                                             make_sed, colors_settings%colors_results_directory)
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

   end subroutine data_for_colors_history_columns

end module colors_history
