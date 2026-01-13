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

   use const_def, only: dp, clight
   use utils_lib, only: mkdir, folder_exists
   use colors_utils, only: remove_dat, romberg_integration
   use knn_interp, only: interpolate_array

   implicit none

   private
   public :: calculate_synthetic
   ! Export zero-point computation functions for precomputation at initialization
   public :: compute_vega_zero_point, compute_ab_zero_point, compute_st_zero_point

   ! ============================================================
   ! INTERNAL SWITCH: Set to .true. to write unique SED files for 
   ! each history row. Each SED file will have _NNNNNN suffix 
   ! matching model_number, allowing linking to history file rows.
   ! ============================================================
character(len=256), save :: first_filter_name = ''

contains

   !****************************
   ! Calculate Synthetic Photometry Using SED and Filter
   ! Uses precomputed zero-point from filter data
   !****************************
real(dp) function calculate_synthetic(temperature, gravity, metallicity, ierr, &
                                      wavelengths, fluxes, &
                                      filter_wavelengths, filter_trans, &
                                      zero_point_flux, &
                                      filter_name, make_sed, sed_per_model, colors_results_directory, model_number)
      ! Input arguments
      real(dp), intent(in) :: temperature, gravity, metallicity
      character(len=*), intent(in) :: filter_name, colors_results_directory
      integer, intent(out) :: ierr

      real(dp), dimension(:), intent(in) :: wavelengths, fluxes
      real(dp), dimension(:), intent(in) :: filter_wavelengths, filter_trans
      real(dp), intent(in) :: zero_point_flux  ! precomputed at initialization
      logical, intent(in) :: make_sed, sed_per_model
      integer, intent(in) :: model_number

      ! Local variables
      real(dp), dimension(:), allocatable :: convolved_flux, filter_on_sed_grid
      character(len=256) :: csv_file
      character(len=20) :: model_str
      character(len=1000) :: line
      real(dp) :: synthetic_flux
      integer :: max_size, i
      real(dp) :: wv, fl, cf, fwv, ftr

      ierr = 0

      ! Allocate working arrays
      allocate(convolved_flux(size(wavelengths)))
      allocate(filter_on_sed_grid(size(wavelengths)))

      ! Interpolate filter onto SED wavelength grid
      call interpolate_array(filter_wavelengths, filter_trans, wavelengths, filter_on_sed_grid)

      ! Convolve SED with filter
      convolved_flux = fluxes * filter_on_sed_grid

      ! Write SED to CSV if requested
      if (make_sed) then
         if (.not. folder_exists(trim(colors_results_directory))) call mkdir(trim(colors_results_directory))
         
         ! Track model number internally when write_sed_per_model is enabled
         if (sed_per_model) then
            
               write(model_str, '(I8.8)') model_number
               csv_file = trim(colors_results_directory)//'/'//trim(remove_dat(filter_name))//'_SED_'//trim(model_str)//'.csv'

         else
            csv_file = trim(colors_results_directory)//'/'//trim(remove_dat(filter_name))//'_SED.csv'
         end if

         max_size = max(size(wavelengths), size(filter_wavelengths))

         open (unit=10, file=csv_file, status='REPLACE', action='write', iostat=ierr)
         if (ierr /= 0) then
            print *, "Error opening file for writing"
            deallocate(convolved_flux, filter_on_sed_grid)
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

      ! Calculate synthetic flux using photon-counting integration
      call calculate_synthetic_flux(wavelengths, convolved_flux, filter_on_sed_grid, synthetic_flux)

      ! Calculate magnitude
      if (zero_point_flux > 0.0_dp .and. synthetic_flux > 0.0_dp) then
         calculate_synthetic = -2.5d0 * log10(synthetic_flux / zero_point_flux)
      else
         if (zero_point_flux <= 0.0_dp) then
            print *, "Error: Zero point flux is zero or negative for filter ", trim(filter_name)
         end if
         if (synthetic_flux <= 0.0_dp) then
            print *, "Error: Synthetic flux is zero or negative for filter ", trim(filter_name)
         end if
         calculate_synthetic = huge(1.0_dp)
      end if

      deallocate(convolved_flux, filter_on_sed_grid)
   end function calculate_synthetic

   !****************************
   ! Calculate Synthetic Flux (photon-counting integration)
   !****************************
   subroutine calculate_synthetic_flux(wavelengths, convolved_flux, filter_on_sed_grid, synthetic_flux)
      real(dp), dimension(:), intent(in) :: wavelengths, convolved_flux, filter_on_sed_grid
      real(dp), intent(out) :: synthetic_flux

      real(dp) :: integrated_flux, integrated_filter

      ! Photon-counting: weight by wavelength
      call romberg_integration(wavelengths, convolved_flux * wavelengths, integrated_flux)
      call romberg_integration(wavelengths, filter_on_sed_grid * wavelengths, integrated_filter)

      if (integrated_filter > 0.0_dp) then
         synthetic_flux = integrated_flux / integrated_filter
      else
         print *, "Error: Integrated filter transmission is zero."
         synthetic_flux = -1.0_dp
      end if
   end subroutine calculate_synthetic_flux

   !****************************
   ! Compute Vega Zero Point Flux
   ! Called once at initialization, result cached in filter_data
   !****************************
   real(dp) function compute_vega_zero_point(vega_wave, vega_flux, filt_wave, filt_trans)
      real(dp), dimension(:), intent(in) :: vega_wave, vega_flux
      real(dp), dimension(:), intent(in) :: filt_wave, filt_trans

      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: filt_on_vega_grid(:), conv_flux(:)

      allocate(filt_on_vega_grid(size(vega_wave)))
      allocate(conv_flux(size(vega_wave)))

      ! Interpolate filter onto Vega wavelength grid
      call interpolate_array(filt_wave, filt_trans, vega_wave, filt_on_vega_grid)

      ! Convolve Vega with filter
      conv_flux = vega_flux * filt_on_vega_grid

      ! Photon-counting integration
      call romberg_integration(vega_wave, vega_wave * conv_flux, int_flux)
      call romberg_integration(vega_wave, vega_wave * filt_on_vega_grid, int_filter)

      if (int_filter > 0.0_dp) then
         compute_vega_zero_point = int_flux / int_filter
      else
         compute_vega_zero_point = -1.0_dp
      end if

      deallocate(filt_on_vega_grid, conv_flux)
   end function compute_vega_zero_point

   !****************************
   ! Compute AB Zero Point Flux
   ! f_nu = 3631 Jy = 3.631e-20 erg/s/cm^2/Hz
   ! f_lambda = f_nu * c / lambda^2
   ! Called once at initialization, result cached in filter_data
   !****************************
   real(dp) function compute_ab_zero_point(filt_wave, filt_trans)
      real(dp), dimension(:), intent(in) :: filt_wave, filt_trans

      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: ab_sed_flux(:)
      integer :: i

      allocate(ab_sed_flux(size(filt_wave)))

      ! Construct AB spectrum (f_lambda) on the filter wavelength grid
      ! 3631 Jy = 3.631E-20 erg/s/cm^2/Hz
      ! clight in cm/s, wavelength in Angstroms, need to convert
      do i = 1, size(filt_wave)
         if (filt_wave(i) > 0.0_dp) then
            ab_sed_flux(i) = 3.631d-20 * ((clight * 1.0d8) / (filt_wave(i)**2))
         else
            ab_sed_flux(i) = 0.0_dp
         end if
      end do

      ! Photon-counting integration
      call romberg_integration(filt_wave, ab_sed_flux * filt_trans * filt_wave, int_flux)
      call romberg_integration(filt_wave, filt_wave * filt_trans, int_filter)

      if (int_filter > 0.0_dp) then
         compute_ab_zero_point = int_flux / int_filter
      else
         compute_ab_zero_point = -1.0_dp
      end if

      deallocate(ab_sed_flux)
   end function compute_ab_zero_point

   !****************************
   ! Compute ST Zero Point Flux
   ! f_lambda = 3.63e-9 erg/s/cm^2/A (Constant)
   ! Called once at initialization, result cached in filter_data
   !****************************
   real(dp) function compute_st_zero_point(filt_wave, filt_trans)
      real(dp), dimension(:), intent(in) :: filt_wave, filt_trans

      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: st_sed_flux(:)

      allocate(st_sed_flux(size(filt_wave)))
      st_sed_flux = 3.63d-9

      ! Photon-counting integration
      call romberg_integration(filt_wave, st_sed_flux * filt_trans * filt_wave, int_flux)
      call romberg_integration(filt_wave, filt_wave * filt_trans, int_filter)

      if (int_filter > 0.0_dp) then
         compute_st_zero_point = int_flux / int_filter
      else
         compute_st_zero_point = -1.0_dp
      end if

      deallocate(st_sed_flux)
   end function compute_st_zero_point

end module synthetic