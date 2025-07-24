! ***********************************************************************
!
!   Copyright (C) 2025  Niall Miller & The MESA Team
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

module synthetic

   use const_def, only: dp
   use utils_lib, only: mkdir, folder_exists
   use shared_funcs, only: remove_dat, romberg_integration, load_filter, load_vega_sed, load_lookup_table
   use knn_interp, only: interpolate_array

   implicit none

   private
   public :: calculate_synthetic

contains

   !****************************
   ! Calculate Synthetic Photometry Using SED and Filter
   !****************************
   real(dp) function calculate_synthetic(temperature, gravity, metallicity, ierr, &
                                         wavelengths, fluxes, filter_wavelengths, &
                                         filter_trans, &
                                         filter_filepath, vega_filepath, &
                                         filter_name, make_sed, colors_results_directory)
      ! Input arguments
      real(dp), intent(in) :: temperature, gravity, metallicity
      character(len=*), intent(in) :: filter_filepath, filter_name, vega_filepath, colors_results_directory
      integer, intent(out) :: ierr
      character(len=1000) :: line

      real(dp), dimension(:), INTENT(INOUT) :: wavelengths, fluxes
      real(dp), dimension(:), allocatable, INTENT(INOUT) :: filter_wavelengths, filter_trans
      logical, intent(in) :: make_sed

      ! Local variables
      real(dp), dimension(:), allocatable :: convolved_flux
      character(len=100) :: csv_file
      real(dp) :: synthetic_flux, vega_flux
      integer :: max_size, i
      real(dp) :: wv, fl, cf, fwv, ftr

      if (.not. folder_exists(trim(colors_results_directory))) call mkdir(trim(colors_results_directory))
      csv_file = trim(colors_results_directory)//'/'//trim(remove_dat(filter_name))//'_SED.csv'
      ierr = 0

      ! Load filter data
      call load_filter(filter_filepath, filter_wavelengths, filter_trans)

      ! Check for invalid gravity input
      if (gravity <= 0.0_dp) then
         ierr = 1
         calculate_synthetic = -1.0_dp
         return
      end if

      ! Perform SED convolution
      allocate (convolved_flux(size(wavelengths)))
      call convolve_sed(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)

      ! Write SED to CSV if requested
      if (make_sed) then
         ! Determine the maximum size among all arrays
         max_size = max(size(wavelengths), size(filter_wavelengths), &
                        size(fluxes), size(convolved_flux), size(filter_trans))

         ! Open the CSV file for writing
         open (unit=10, file=csv_file, status='REPLACE', action='write', iostat=ierr)
         if (ierr /= 0) then
            print *, "Error opening file for writing"
            return
         end if

         ! Write headers to the CSV file
         write (10, '(A)') "wavelengths,fluxes,convolved_flux,filter_wavelengths,filter_trans"

         ! Loop through data and safely write values, ensuring no out-of-bounds errors
         do i = 1, max_size
            ! Initialize values to zero in case they are out of bounds
            wv = 0.0_dp
            fl = 0.0_dp
            cf = 0.0_dp
            fwv = 0.0_dp
            ftr = 0.0_dp

            ! Assign actual values only if within valid indices
            if (i <= size(wavelengths)) wv = wavelengths(i)
            if (i <= size(fluxes)) fl = fluxes(i)
            if (i <= size(convolved_flux)) cf = convolved_flux(i)
            if (i <= size(filter_wavelengths)) fwv = filter_wavelengths(i)
            if (i <= size(filter_trans)) ftr = filter_trans(i)

            ! Write the formatted output
            write (line, '(ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6)') &
               wv, fl, cf, fwv, ftr
            write (10, '(A)') trim(line)
         end do

         ! Close the file
         close (10)
      end if

      ! Calculate Vega flux for zero point calibration
      vega_flux = calculate_vega_flux(vega_filepath, filter_wavelengths, filter_trans, &
                                      filter_name, make_sed, colors_results_directory)

      ! Calculate synthetic flux
      call calculate_synthetic_flux(wavelengths, convolved_flux, synthetic_flux, &
                                    filter_wavelengths, filter_trans)

      ! Calculate magnitude using Vega zero point
      if (vega_flux > 0.0_dp) then
         calculate_synthetic = -2.5*LOG10(synthetic_flux/vega_flux)
      else
         print *, "Error: Vega flux is zero, magnitude calculation is invalid."
         calculate_synthetic = HUGE(1.0_dp)
      end if

      ! Clean up
      deallocate (convolved_flux)
   end function calculate_synthetic

   !-----------------------------------------------------------------------
   ! Internal functions for synthetic photometry
   !-----------------------------------------------------------------------

   !****************************
   ! Convolve SED With Filter
   !****************************
   subroutine convolve_sed(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)
      real(dp), dimension(:), INTENT(INOUT) :: wavelengths, fluxes
      real(dp), dimension(:), INTENT(INOUT) :: filter_wavelengths, filter_trans
      real(dp), dimension(:), allocatable, intent(out) :: convolved_flux
      real(dp), dimension(:), allocatable :: interpolated_filter
      integer :: n

      n = size(wavelengths)

      ! Allocate arrays
      allocate (interpolated_filter(n))

      ! Interpolate the filter transmission onto the wavelengths array
      call interpolate_array(filter_wavelengths, filter_trans, wavelengths, interpolated_filter)

      ! Perform convolution (element-wise multiplication)
      convolved_flux = fluxes*interpolated_filter

      ! Deallocate temporary arrays
      deallocate (interpolated_filter)
   end subroutine convolve_sed

   !****************************
   ! Calculate Synthetic Flux
   !****************************
   subroutine calculate_synthetic_flux(wavelengths, fluxes, synthetic_flux, &
                                       filter_wavelengths, filter_trans)

      real(dp), dimension(:), intent(in) :: wavelengths, fluxes
      real(dp), dimension(:), INTENT(INOUT) :: filter_wavelengths, filter_trans
      real(dp), intent(out) :: synthetic_flux
      integer :: i
      real(dp) :: integrated_flux, integrated_filter

      ! Validate inputs
      do i = 1, size(wavelengths) - 1
         if (wavelengths(i) <= 0.0 .OR. fluxes(i) < 0.0) then
            print *, "synthetic Invalid input at index", i, ":", wavelengths(i), fluxes(i)
            stop
         end if
      end do

      call romberg_integration(wavelengths, fluxes*wavelengths, integrated_flux)
      call romberg_integration(filter_wavelengths, &
                               filter_trans*filter_wavelengths, integrated_filter)
      ! Store the total flux
      if (integrated_filter > 0.0) then
         synthetic_flux = integrated_flux/integrated_filter
      else
         print *, "Error: Integrated filter transmission is zero."
         synthetic_flux = -1.0_dp
         return
      end if
   end subroutine calculate_synthetic_flux

   !****************************
   ! Calculate Vega Flux for Zero Point
   !****************************
   function calculate_vega_flux(vega_filepath, filt_wave, filt_trans, &
                                filter_name, make_sed, colors_results_directory) result(vega_flux)
      character(len=*), intent(in) :: vega_filepath, filter_name, colors_results_directory
      character(len=100) :: output_csv
      real(dp), dimension(:), INTENT(INOUT) :: filt_wave, filt_trans
      real(dp) :: vega_flux
      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: vega_wave(:), vega_flux_arr(:), conv_flux(:)
      logical, intent(in) :: make_sed
      integer :: i, unit, max_size
      real(dp) :: wv, fl, cf, fwv, ftr
      integer:: ierr
      character(len=1000) :: line

      ! Load the Vega SED
      call load_vega_sed(vega_filepath, vega_wave, vega_flux_arr)

      ! Convolve the Vega SED with the filter transmission
      allocate (conv_flux(size(vega_wave)))
      call convolve_sed(vega_wave, vega_flux_arr, filt_wave, filt_trans, conv_flux)

      ! Integrate the convolved Vega SED and the filter transmission
      call romberg_integration(vega_wave, vega_wave*conv_flux, int_flux)
      call romberg_integration(filt_wave, filt_wave*filt_trans, int_filter)

      if (int_filter > 0.0_dp) then
         vega_flux = int_flux/int_filter
      else
         vega_flux = -1.0_dp
      end if

      ! Write Vega SED to CSV if requested
      if (make_sed) then
         ! Determine the maximum size among all arrays
         max_size = max(size(vega_wave), size(vega_flux_arr), size(conv_flux), &
                        size(filt_wave), size(filt_trans))

         if (.not. folder_exists(trim(colors_results_directory))) call mkdir(trim(colors_results_directory))
         output_csv = trim(colors_results_directory)//'/VEGA_'//trim(remove_dat(filter_name))//'_SED.csv'

         ! Open the CSV file for writing
         open (unit=10, file=output_csv, status='REPLACE', action='write', iostat=ierr)
         if (ierr /= 0) then
            print *, "Error opening file for writing"
            return
         end if

         write (10, '(A)') "wavelengths,fluxes,convolved_flux,filter_wavelengths,filter_trans"

         ! Loop through data and safely write values, ensuring no out-of-bounds errors
         do i = 1, max_size
            ! Initialize values to zero in case they are out of bounds
            wv = 0.0_dp
            fl = 0.0_dp
            cf = 0.0_dp
            fwv = 0.0_dp
            ftr = 0.0_dp

            ! Assign actual values only if within valid indices
            if (i <= size(vega_wave)) wv = vega_wave(i)
            if (i <= size(vega_flux_arr)) fl = vega_flux_arr(i)
            if (i <= size(conv_flux)) cf = conv_flux(i)
            if (i <= size(filt_wave)) fwv = filt_wave(i)
            if (i <= size(filt_trans)) ftr = filt_trans(i)

            ! Write the formatted output
            write (line, '(ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6)') &
               wv, fl, cf, fwv, ftr
            write (10, '(A)') trim(line)
         end do

         ! Close the file
         close (10)
      end if

      ! Clean up
      deallocate (conv_flux, vega_wave, vega_flux_arr)
   end function calculate_vega_flux

end module synthetic
