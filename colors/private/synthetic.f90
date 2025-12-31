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

module synthetic

   ! Added clight for AB conversions
   use const_def, only: dp, clight
   use utils_lib, only: mkdir, folder_exists
   use colors_utils, only: remove_dat, romberg_integration
   use knn_interp, only: interpolate_array

   implicit none

   private
   public :: calculate_synthetic

contains

   !****************************
   ! Calculate Synthetic Photometry Using SED and Filter
   ! Now accepts cached filter and Vega data instead of file paths
   !****************************
   real(dp) function calculate_synthetic(temperature, gravity, metallicity, ierr, &
                                         wavelengths, fluxes, &
                                         filter_wavelengths, filter_trans, &
                                         vega_wavelengths, vega_fluxes, &
                                         filter_name, make_sed, colors_results_directory, &
                                         mag_system)
      ! Input arguments
      real(dp), intent(in) :: temperature, gravity, metallicity
      character(len=*), intent(in) :: filter_name, colors_results_directory
      character(len=*), intent(in) :: mag_system  ! 'Vega', 'AB', or 'ST'
      integer, intent(out) :: ierr
      character(len=1000) :: line

      real(dp), dimension(:), intent(in) :: wavelengths, fluxes
      ! Cached filter data (passed in from colors_settings)
      real(dp), dimension(:), intent(in) :: filter_wavelengths, filter_trans
      ! Cached Vega SED (passed in from colors_settings)
      real(dp), dimension(:), intent(in) :: vega_wavelengths, vega_fluxes
      logical, intent(in) :: make_sed

      ! Local variables
      real(dp), dimension(:), allocatable :: convolved_flux
      ! Local copies for routines that need intent(inout)
      real(dp), dimension(:), allocatable :: local_filt_wave, local_filt_trans
      real(dp), dimension(:), allocatable :: local_wavelengths, local_fluxes
      character(len=100) :: csv_file
      real(dp) :: synthetic_flux, zero_point_flux
      integer :: max_size, i
      real(dp) :: wv, fl, cf, fwv, ftr

      if (.not. folder_exists(trim(colors_results_directory))) call mkdir(trim(colors_results_directory))
      csv_file = trim(colors_results_directory)//'/'//trim(remove_dat(filter_name))//'_SED.csv'
      ierr = 0

      ! Make local copies since some routines need intent(inout)
      allocate(local_wavelengths(size(wavelengths)))
      allocate(local_fluxes(size(fluxes)))
      allocate(local_filt_wave(size(filter_wavelengths)))
      allocate(local_filt_trans(size(filter_trans)))
      local_wavelengths = wavelengths
      local_fluxes = fluxes
      local_filt_wave = filter_wavelengths
      local_filt_trans = filter_trans

      ! Perform SED convolution (Source Object)
      allocate(convolved_flux(size(wavelengths)))
      call convolve_sed(local_wavelengths, local_fluxes, local_filt_wave, local_filt_trans, convolved_flux)

      ! Write SED to CSV if requested
      if (make_sed) then
         max_size = max(size(wavelengths), size(filter_wavelengths), &
                        size(fluxes), size(convolved_flux), size(filter_trans))

         open (unit=10, file=csv_file, status='REPLACE', action='write', iostat=ierr)
         if (ierr /= 0) then
            print *, "Error opening file for writing"
            return
         end if

         write (10, '(A)') "wavelengths,fluxes,convolved_flux,filter_wavelengths,filter_trans"

         do i = 1, max_size
            wv = 0.0_dp; fl = 0.0_dp; cf = 0.0_dp; fwv = 0.0_dp; ftr = 0.0_dp
            if (i <= size(wavelengths)) wv = wavelengths(i)
            if (i <= size(fluxes)) fl = fluxes(i)
            if (i <= size(convolved_flux)) cf = convolved_flux(i)
            if (i <= size(filter_wavelengths)) fwv = filter_wavelengths(i)
            if (i <= size(filter_trans)) ftr = filter_trans(i)

            write (line, '(ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6)') &
               wv, fl, cf, fwv, ftr
            write (10, '(A)') trim(line)
         end do
         close (10)
      end if

      select case (trim(mag_system))
      case ('VEGA', 'Vega', 'vega')
         zero_point_flux = calculate_vega_flux(vega_wavelengths, vega_fluxes, &
                                               local_filt_wave, local_filt_trans, &
                                               filter_name, make_sed, colors_results_directory)
      case ('AB', 'ab')
         zero_point_flux = calculate_ab_zero_point(local_filt_wave, local_filt_trans)
      case ('ST', 'st')
         zero_point_flux = calculate_st_zero_point(local_filt_wave, local_filt_trans)
      case default
         print *, "Error: Unknown magnitude system: ", mag_system
         calculate_synthetic = huge(1.0_dp)
         return
      end select

      call calculate_synthetic_flux(local_wavelengths, convolved_flux, synthetic_flux, &
                                    local_filt_wave, local_filt_trans)

      if (zero_point_flux > 0.0_dp) then
         calculate_synthetic = -2.5d0*log10(synthetic_flux/zero_point_flux)
      else
         print *, "Error: Zero point flux is zero, magnitude calculation is invalid."
         calculate_synthetic = huge(1.0_dp)
      end if

      ! Clean up
      deallocate(convolved_flux)
      deallocate(local_wavelengths, local_fluxes)
      deallocate(local_filt_wave, local_filt_trans)
   end function calculate_synthetic

   !****************************
   ! Convolve SED With Filter
   !****************************
   subroutine convolve_sed(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)
      real(dp), dimension(:), intent(inout) :: wavelengths, fluxes
      real(dp), dimension(:), intent(inout) :: filter_wavelengths, filter_trans
      real(dp), dimension(:), allocatable, intent(out) :: convolved_flux
      real(dp), dimension(:), allocatable :: interpolated_filter
      integer :: n

      n = size(wavelengths)
      allocate (interpolated_filter(n))
      call interpolate_array(filter_wavelengths, filter_trans, wavelengths, interpolated_filter)
      convolved_flux = fluxes*interpolated_filter
      deallocate (interpolated_filter)
   end subroutine convolve_sed

   !****************************
   ! Calculate Synthetic Flux (Integration)
   !****************************
   subroutine calculate_synthetic_flux(wavelengths, fluxes, synthetic_flux, &
                                       filter_wavelengths, filter_trans)

      real(dp), dimension(:), intent(in) :: wavelengths, fluxes
      real(dp), dimension(:), intent(inout) :: filter_wavelengths, filter_trans
      real(dp), intent(out) :: synthetic_flux
      integer :: i
      real(dp) :: integrated_flux, integrated_filter

      real(dp), dimension(:), allocatable :: filter_on_sed_grid

      allocate (filter_on_sed_grid(size(wavelengths)))

      ! Validate inputs
      do i = 1, size(wavelengths) - 1
         if (wavelengths(i) <= 0.0_dp .or. fluxes(i) < 0.0_dp) then
            print *, "synthetic Invalid input at index", i, ":", wavelengths(i), fluxes(i)
            stop
         end if
      end do

      call interpolate_array(filter_wavelengths, filter_trans, wavelengths, filter_on_sed_grid)

      call romberg_integration(wavelengths, fluxes*wavelengths, integrated_flux)
      call romberg_integration(wavelengths, filter_on_sed_grid*wavelengths, integrated_filter)

      if (integrated_filter > 0.0_dp) then
         synthetic_flux = integrated_flux/integrated_filter
      else
         print *, "Error: Integrated filter transmission is zero."
         synthetic_flux = -1.0_dp
         return
      end if
   end subroutine calculate_synthetic_flux

   !****************************
   ! Calculate Vega Flux for Zero Point
   ! Now accepts cached Vega SED instead of file path
   !****************************
   function calculate_vega_flux(vega_wave, vega_flux_arr, filt_wave, filt_trans, &
                                filter_name, make_sed, colors_results_directory) result(vega_flux)
      real(dp), dimension(:), intent(in) :: vega_wave, vega_flux_arr
      real(dp), dimension(:), intent(inout) :: filt_wave, filt_trans
      character(len=*), intent(in) :: filter_name, colors_results_directory
      logical, intent(in) :: make_sed

      real(dp) :: vega_flux
      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: conv_flux(:)
      real(dp), allocatable :: filt_trans_on_vega_grid(:)
      real(dp), allocatable :: local_vega_wave(:), local_vega_flux(:)
      character(len=100) :: output_csv
      integer :: i, max_size, ierr
      real(dp) :: wv, fl, cf, fwv, ftr
      character(len=1000) :: line

      ! Make local copies for routines that need intent(inout)
      allocate(local_vega_wave(size(vega_wave)))
      allocate(local_vega_flux(size(vega_flux_arr)))
      local_vega_wave = vega_wave
      local_vega_flux = vega_flux_arr

      allocate(conv_flux(size(vega_wave)))
      call convolve_sed(local_vega_wave, local_vega_flux, filt_wave, filt_trans, conv_flux)

      allocate(filt_trans_on_vega_grid(size(vega_wave)))
      call interpolate_array(filt_wave, filt_trans, local_vega_wave, filt_trans_on_vega_grid)

      call romberg_integration(local_vega_wave, local_vega_wave*conv_flux, int_flux)
      call romberg_integration(local_vega_wave, local_vega_wave*filt_trans_on_vega_grid, int_filter)

      if (int_filter > 0.0_dp) then
         vega_flux = int_flux/int_filter
      else
         vega_flux = -1.0_dp
      end if

      ! Write Vega SED to CSV if requested
      if (make_sed) then
         max_size = max(size(vega_wave), size(vega_flux_arr), size(conv_flux), &
                        size(filt_wave), size(filt_trans))

         if (.not. folder_exists(trim(colors_results_directory))) call mkdir(trim(colors_results_directory))
         output_csv = trim(colors_results_directory)//'/VEGA_'//trim(remove_dat(filter_name))//'_SED.csv'

         open (unit=10, file=output_csv, status='REPLACE', action='write', iostat=ierr)
         if (ierr /= 0) then
            print *, "Error opening file for writing"
            return
         end if

         write (10, '(A)') "wavelengths,fluxes,convolved_flux,filter_wavelengths,filter_trans"

         do i = 1, max_size
            wv = 0.0_dp; fl = 0.0_dp; cf = 0.0_dp; fwv = 0.0_dp; ftr = 0.0_dp
            if (i <= size(vega_wave)) wv = vega_wave(i)
            if (i <= size(vega_flux_arr)) fl = vega_flux_arr(i)
            if (i <= size(conv_flux)) cf = conv_flux(i)
            if (i <= size(filt_wave)) fwv = filt_wave(i)
            if (i <= size(filt_trans)) ftr = filt_trans(i)

            write (line, '(ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6)') &
               wv, fl, cf, fwv, ftr
            write (10, '(A)') trim(line)
         end do
         close (10)
      end if

      deallocate(conv_flux, local_vega_wave, local_vega_flux, filt_trans_on_vega_grid)
   end function calculate_vega_flux

   !****************************
   ! Calculate AB Zero Point Flux
   ! f_nu = 3631 Jy = 3.631e-20 erg/s/cm^2/Hz
   ! f_lambda = f_nu * c / lambda^2
   !****************************
   function calculate_ab_zero_point(filt_wave, filt_trans) result(ab_flux)
      real(dp), dimension(:), intent(inout) :: filt_wave, filt_trans
      real(dp) :: ab_flux
      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: ab_sed_flux(:)
      integer :: i

      allocate (ab_sed_flux(size(filt_wave)))

      ! Construct AB Spectrum (f_lambda) on the filter wavelength grid
      ! Assumes wavelengths are in Angstroms and clight is in Angstroms/sec
      ! 3631 Jy = 3.631E-20 erg/s/cm^2/Hz
      do i = 1, size(filt_wave)
         if (filt_wave(i) > 0.0_dp) then
            ab_sed_flux(i) = 3.631d-20*((clight*1.0e8)/(filt_wave(i)**2))
         else
            ab_sed_flux(i) = 0.0_dp
         end if
      end do

      call romberg_integration(filt_wave, ab_sed_flux*filt_trans*filt_wave, int_flux)
      call romberg_integration(filt_wave, filt_wave*filt_trans, int_filter)

      if (int_filter > 0.0_dp) then
         ab_flux = int_flux/int_filter
      else
         ab_flux = -1.0_dp
      end if

      deallocate (ab_sed_flux)
   end function calculate_ab_zero_point

   !****************************
   ! Calculate ST Zero Point Flux
   ! f_lambda = 3.63e-9 erg/s/cm^2/A (Constant)
   !****************************
   function calculate_st_zero_point(filt_wave, filt_trans) result(st_flux)
      real(dp), dimension(:), intent(inout) :: filt_wave, filt_trans
      real(dp) :: st_flux
      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: st_sed_flux(:)

      allocate (st_sed_flux(size(filt_wave)))
      st_sed_flux = 3.63d-9

      call romberg_integration(filt_wave, st_sed_flux*filt_trans*filt_wave, int_flux)
      call romberg_integration(filt_wave, filt_wave*filt_trans, int_filter)

      if (int_filter > 0.0_dp) then
         st_flux = int_flux/int_filter
      else
         st_flux = -1.0_dp
      end if

      deallocate (st_sed_flux)
   end function calculate_st_zero_point

end module synthetic