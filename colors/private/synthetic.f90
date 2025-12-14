! ***********************************************************************
!
!   Copyright (C) 2025  Niall Miller & The MESA Team
!   Modified to include AB and ST magnitude systems.
!
! ***********************************************************************

module synthetic

   ! Added clight for AB conversions
   use const_def, only: dp, clight
   use utils_lib, only: mkdir, folder_exists
   use colors_utils, only: remove_dat, romberg_integration, load_filter, load_vega_sed, load_lookup_table
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
                                         filter_name, make_sed, colors_results_directory, &
                                         mag_system)
      ! Input arguments
      real(dp), intent(in) :: temperature, gravity, metallicity
      character(len=*), intent(in) :: filter_filepath, filter_name, vega_filepath, colors_results_directory
      character(len=*), intent(in) :: mag_system  ! NEW: Arguments for 'Vega', 'AB', or 'ST'
      integer, intent(out) :: ierr
      character(len=1000) :: line

      real(dp), dimension(:), intent(inout) :: wavelengths, fluxes
      real(dp), dimension(:), allocatable, intent(inout) :: filter_wavelengths, filter_trans
      logical, intent(in) :: make_sed

      ! Local variables
      real(dp), dimension(:), allocatable :: convolved_flux
      character(len=100) :: csv_file
      real(dp) :: synthetic_flux, zero_point_flux
      integer :: max_size, i
      real(dp) :: wv, fl, cf, fwv, ftr

      if (.not. folder_exists(trim(colors_results_directory))) call mkdir(trim(colors_results_directory))
      csv_file = trim(colors_results_directory)//'/'//trim(remove_dat(filter_name))//'_SED.csv'
      ierr = 0

      ! Load filter data
      call load_filter(filter_filepath, filter_wavelengths, filter_trans)

      ! Perform SED convolution (Source Object)
      allocate (convolved_flux(size(wavelengths)))
      call convolve_sed(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)

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

      ! ------------------------------------------------------------------
      ! Calculate Zero Point Flux based on System Selection
      ! ------------------------------------------------------------------
      select case (trim(mag_system))
      case ('VEGA', 'Vega', 'vega')
         zero_point_flux = calculate_vega_flux(vega_filepath, filter_wavelengths, filter_trans, &
                                               filter_name, make_sed, colors_results_directory)
      case ('AB', 'ab')
         zero_point_flux = calculate_ab_zero_point(filter_wavelengths, filter_trans)
      case ('ST', 'st')
         zero_point_flux = calculate_st_zero_point(filter_wavelengths, filter_trans)
      case default
         print *, "Error: Unknown magnitude system: ", mag_system
         calculate_synthetic = huge(1.0_dp)
         return
      end select

      ! Calculate synthetic flux (Source Object)
      call calculate_synthetic_flux(wavelengths, convolved_flux, synthetic_flux, &
                                    filter_wavelengths, filter_trans)

      ! Calculate magnitude using the selected zero point
      if (zero_point_flux > 0.0_dp) then
         calculate_synthetic = -2.5d0*log10(synthetic_flux/zero_point_flux)
      else
         print *, "Error: Zero point flux is zero, magnitude calculation is invalid."
         calculate_synthetic = huge(1.0_dp)
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

      ! Validate inputs
      do i = 1, size(wavelengths) - 1
         if (wavelengths(i) <= 0.0_dp .or. fluxes(i) < 0.0_dp) then
            print *, "synthetic Invalid input at index", i, ":", wavelengths(i), fluxes(i)
            stop
         end if
      end do

      call romberg_integration(wavelengths, fluxes*wavelengths, integrated_flux)
      call romberg_integration(filter_wavelengths, &
                               filter_trans*filter_wavelengths, integrated_filter)

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
   !****************************
   function calculate_vega_flux(vega_filepath, filt_wave, filt_trans, &
                                filter_name, make_sed, colors_results_directory) result(vega_flux)
      character(len=*), intent(in) :: vega_filepath, filter_name, colors_results_directory
      character(len=100) :: output_csv
      real(dp), dimension(:), intent(inout) :: filt_wave, filt_trans
      real(dp) :: vega_flux
      real(dp) :: int_flux, int_filter
      real(dp), allocatable :: vega_wave(:), vega_flux_arr(:), conv_flux(:)
      logical, intent(in) :: make_sed
      integer :: i, max_size, ierr
      real(dp) :: wv, fl, cf, fwv, ftr
      character(len=1000) :: line

      ! Load the Vega SED
      call load_vega_sed(vega_filepath, vega_wave, vega_flux_arr)

      ! Convolve the Vega SED with the filter transmission
      allocate (conv_flux(size(vega_wave)))
      call convolve_sed(vega_wave, vega_flux_arr, filt_wave, filt_trans, conv_flux)

      ! Integrate
      call romberg_integration(vega_wave, vega_wave*conv_flux, int_flux)
      call romberg_integration(filt_wave, filt_wave*filt_trans, int_filter)

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
      deallocate (conv_flux, vega_wave, vega_flux_arr)
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

      allocate(ab_sed_flux(size(filt_wave)))

      ! Construct AB Spectrum (f_lambda) on the filter wavelength grid
      ! Assumes wavelengths are in Angstroms and clight is in Angstroms/sec
      ! 3631 Jy = 3.631E-20 erg/s/cm^2/Hz
      do i = 1, size(filt_wave)
         if (filt_wave(i) > 0.0_dp) then
            ab_sed_flux(i) = 3.631d-20 * (clight / (filt_wave(i)**2))
         else
            ab_sed_flux(i) = 0.0_dp
         endif
      end do

      ! Integrate using same method as source (f_lambda * T * lambda)
      ! Note: We multiply by filt_wave inside the integration because the
      ! romberg helper expects (flux * lambda)
      call romberg_integration(filt_wave, ab_sed_flux * filt_trans * filt_wave, int_flux)
      call romberg_integration(filt_wave, filt_wave * filt_trans, int_filter)

      if (int_filter > 0.0_dp) then
         ab_flux = int_flux / int_filter
      else
         ab_flux = -1.0_dp
      end if

      deallocate(ab_sed_flux)
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

      ! For ST system, flux is constant in wavelength
      ! However, to maintain exact consistency with how the source is integrated
      ! (numerical integration over the filter grid), we integrate the constant array.

      allocate(st_sed_flux(size(filt_wave)))
      st_sed_flux = 3.63d-9

      call romberg_integration(filt_wave, st_sed_flux * filt_trans * filt_wave, int_flux)
      call romberg_integration(filt_wave, filt_wave * filt_trans, int_filter)

      if (int_filter > 0.0_dp) then
         st_flux = int_flux / int_filter
      else
         st_flux = -1.0_dp
      end if

      deallocate(st_sed_flux)
   end function calculate_st_zero_point

end module synthetic