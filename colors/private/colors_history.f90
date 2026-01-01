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
         num_cols = 3 + num_color_filters
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

      type(Colors_General_Info), pointer :: cs  ! colors_settings
      integer :: i, filter_offset
      real(dp) :: d, bolometric_magnitude, bolometric_flux, interpolation_radius
      real(dp) :: zero_point
      character(len=256) :: sed_filepath, filter_name
      logical :: make_sed

      real(dp), dimension(:), allocatable :: wavelengths, fluxes

      ierr = 0
      call get_colors_ptr(colors_handle, cs, ierr)
      if (ierr /= 0) then
         write (*, *) 'failed in get_colors_ptr'
         return
      end if

      ! Safety check: verify data was loaded at initialization
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
      sed_filepath = trim(mesa_dir)//cs%stellar_atm
      make_sed = cs%make_csv

      ! Calculate bolometric magnitude using cached lookup table
      call calculate_bolometric(t_eff, log_g, metallicity, R, d, &
                                bolometric_magnitude, bolometric_flux, wavelengths, fluxes, &
                                sed_filepath, interpolation_radius, &
                                cs%lu_file_names, cs%lu_teff, cs%lu_logg, cs%lu_meta)

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

            if (t_eff >= 0 .and. metallicity >= 0) then
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

               ! Calculate synthetic magnitude using cached filter data and precomputed zero-point
               vals(i + filter_offset) = calculate_synthetic(t_eff, log_g, metallicity, ierr, &
                                                             wavelengths, fluxes, &
                                                             cs%filters(i)%wavelengths, &
                                                             cs%filters(i)%transmission, &
                                                             zero_point, &
                                                             color_filter_names(i), &
                                                             make_sed, cs%colors_results_directory)
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

      ! Clean up allocated arrays from calculate_bolometric
      if (allocated(wavelengths)) deallocate(wavelengths)
      if (allocated(fluxes)) deallocate(fluxes)

   end subroutine data_for_colors_history_columns

end module colors_history