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

module colors_lib

  use const_def, only: dp, strlen
  use colors_def
  use math_lib
  use hermite_interp
  use knn_interp
  use linear_interp
  use shared_funcs

  implicit none
  
  ! Make public interface explicit
  public :: colors_init, colors_shutdown
  public :: calculate_bolometric, calculate_synthetic
  public :: load_lookuptable, remove_dat
  
  ! Keep internals private
  private :: calculate_bolometric_phot
  
  contains

      ! call this routine to initialize the colors module.
      ! only needs to be done once at start of run.
      ! Reads data from the 'colors' directory in the data_dir.
      ! If use_cache is true and there is a 'colors/cache' directory, it will try that first.
      ! If it doesn't find what it needs in the cache,
      ! it reads the data and writes the cache for next time.
  subroutine colors_init(use_cache, colors_cache_dir, ierr)
    use colors_def, only : colors_def_init, colors_use_cache, colors_is_initialized
    logical, intent(in) :: use_cache
    character (len=*), intent(in) :: colors_cache_dir  ! blank means use default
    integer, intent(out) :: ierr  ! 0 means AOK.
    ierr = 0
    if (colors_is_initialized) return
    call colors_def_init(colors_cache_dir)
    colors_use_cache = use_cache
    colors_is_initialized = .true.
 end subroutine colors_init


 subroutine colors_shutdown
    use colors_def, only: do_free_colors_tables, colors_is_initialized
    call do_free_colors_tables()
    colors_is_initialized = .false.
 end subroutine colors_shutdown


      ! after colors_init has finished, you can allocate a "handle".

 integer function alloc_colors_handle(ierr) result(handle)
 integer, intent(out) :: ierr  ! 0 means AOK.
 character (len=0) :: inlist
 handle = alloc_colors_handle_using_inlist(inlist, ierr)
end function alloc_colors_handle

integer function alloc_colors_handle_using_inlist(inlist,ierr) result(handle)
 use colors_def, only: do_alloc_colors, colors_is_initialized
 use colors_ctrls_io, only: read_namelist
 character (len=*), intent(in) :: inlist  ! empty means just use defaults.
 integer, intent(out) :: ierr  ! 0 means AOK.
 ierr = 0
 if (.not. colors_is_initialized) then
    ierr=-1
    return
 endif
 handle = do_alloc_colors(ierr)
 if (ierr /= 0) return
 call read_namelist(handle, inlist, ierr)
 if (ierr /= 0) return
 call colors_setup_tables(handle, ierr)
 call colors_setup_hooks(handle, ierr)
end function alloc_colors_handle_using_inlist

subroutine free_colors_handle(handle)
 ! frees the handle and all associated data
 use colors_def,only: colors_General_Info, do_free_colors
 integer, intent(in) :: handle
 call do_free_colors(handle)
end subroutine free_colors_handle


subroutine colors_ptr(handle,rq,ierr)
 use colors_def,only:Colors_General_Info,get_colors_ptr,colors_is_initialized
 integer, intent(in) :: handle  ! from alloc_colors_handle
 type (colors_General_Info), pointer :: rq
 integer, intent(out):: ierr
 if (.not. colors_is_initialized) then
    ierr=-1
    return
 endif
 call get_colors_ptr(handle,rq,ierr)
end subroutine colors_ptr


subroutine colors_setup_tables(handle, ierr)
 use colors_def, only : colors_General_Info, get_colors_ptr
 ! TODO: use load_colors, only : Setup_colors_Tables
 integer, intent(in) :: handle
 integer, intent(out):: ierr

 type (colors_General_Info), pointer :: rq
 logical, parameter :: use_cache = .true.
 logical, parameter :: load_on_demand = .true.

 ierr = 0
 call get_colors_ptr(handle,rq,ierr)
 ! TODO: call Setup_colors_Tables(rq, use_cache, load_on_demand, ierr)

end subroutine colors_setup_tables


subroutine colors_setup_hooks(handle, ierr)
 use colors_def, only : colors_General_Info, get_colors_ptr
 integer, intent(in) :: handle
 integer, intent(out):: ierr

 type (colors_General_Info), pointer :: rq

 ierr = 0
 call get_colors_ptr(handle,rq,ierr)

 ! TODO: currently does nothing. See kap if this feature is needed

end subroutine colors_setup_hooks




  !-----------------------------------------------------------------------
  ! Main public interface functions
  !-----------------------------------------------------------------------
  
  !****************************
  ! Calculate Bolometric Photometry Using Multiple SEDs
  !****************************
  SUBROUTINE calculate_bolometric(teff, log_g, metallicity, R, d, bolometric_magnitude, &
                                 bolometric_flux, wavelengths, fluxes, sed_filepath)
    REAL(dp), INTENT(IN) :: teff, log_g, metallicity, R, d
    CHARACTER(LEN=*), INTENT(IN) :: sed_filepath
    REAL(dp), INTENT(OUT) :: bolometric_magnitude, bolometric_flux
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, fluxes
    
    REAL(dp), ALLOCATABLE :: lu_logg(:), lu_meta(:), lu_teff(:)
    CHARACTER(LEN=100), ALLOCATABLE :: file_names(:)
    REAL, DIMENSION(:,:), ALLOCATABLE :: lookup_table
    CHARACTER(LEN=256) :: lookup_file

    lookup_file = TRIM(sed_filepath) // '/lookup_table.csv'

    ! Load the lookup table
    CALL load_lookuptable(lookup_file, lookup_table, file_names, lu_logg, lu_meta, lu_teff)
    
    ! Use KNN interpolation for constructing the SED
    CALL constructsed_linear(teff, log_g, metallicity, R, d, file_names, &
                         lu_teff, lu_logg, lu_meta, sed_filepath, wavelengths, fluxes)
    
    ! Calculate bolometric flux and magnitude
    CALL calculate_bolometric_phot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
  END SUBROUTINE calculate_bolometric

  !****************************
  ! Calculate Synthetic Photometry Using SED and Filter
  !****************************
  REAL(dp) FUNCTION calculate_synthetic(temperature, gravity, metallicity, ierr, &
                                      wavelengths, fluxes, filter_wavelengths, &
                                      filter_trans, filter_filepath, vega_filepath, &
                                      filter_name, make_sed)
    ! Input arguments
    REAL(dp), INTENT(IN) :: temperature, gravity, metallicity
    CHARACTER(LEN=*), INTENT(IN) :: filter_filepath, filter_name, vega_filepath
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=1000) :: line

    REAL(dp), DIMENSION(:), INTENT(INOUT) :: wavelengths, fluxes
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: filter_wavelengths, filter_trans
    LOGICAL, INTENT(IN) :: make_sed
    
    ! Local variables
    REAL(dp), DIMENSION(:), ALLOCATABLE :: convolved_flux
    CHARACTER(LEN=100) :: csv_file
    REAL(dp) :: synthetic_flux, vega_flux
    INTEGER :: max_size, i
    REAL(dp) :: wv, fl, cf, fwv, ftr

    csv_file = 'LOGS/SED/' // TRIM(remove_dat(filter_name)) // '_SED.csv'
    
    ! Initialize error flag
    ierr = 0

    ! Load filter data
    CALL loadfilter(filter_filepath, filter_wavelengths, filter_trans)

    ! Check for invalid gravity input
    IF (gravity <= 0.0_dp) THEN
        ierr = 1
        calculate_synthetic = -1.0_dp
        RETURN
    END IF

    ! Perform SED convolution
    ALLOCATE(convolved_flux(SIZE(wavelengths)))
    CALL convolvesed(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)

    ! Write SED to CSV if requested
    IF (make_sed) THEN
      ! Determine the maximum size among all arrays
      max_size = MAX(SIZE(wavelengths), SIZE(filter_wavelengths), &
                     SIZE(fluxes), SIZE(convolved_flux), SIZE(filter_trans))

      ! Open the CSV file for writing
      OPEN(UNIT=10, FILE=csv_file, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
      IF (ierr /= 0) THEN
          PRINT *, "Error opening file for writing"
          RETURN
      END IF

      ! Write headers to the CSV file
      WRITE(10, '(A)') "wavelengths,fluxes,convolved_flux,filter_wavelengths,filter_trans"

      ! Loop through data and safely write values, ensuring no out-of-bounds errors
      DO i = 1, max_size
          ! Initialize values to zero in case they are out of bounds
          wv = 0.0_dp
          fl = 0.0_dp
          cf = 0.0_dp
          fwv = 0.0_dp
          ftr = 0.0_dp

          ! Assign actual values only if within valid indices
          IF (i <= SIZE(wavelengths)) wv = wavelengths(i)
          IF (i <= SIZE(fluxes)) fl = fluxes(i)
          IF (i <= SIZE(convolved_flux)) cf = convolved_flux(i)
          IF (i <= SIZE(filter_wavelengths)) fwv = filter_wavelengths(i)
          IF (i <= SIZE(filter_trans)) ftr = filter_trans(i)

          ! Write the formatted output
          WRITE(line, '(ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6)') &
              wv, fl, cf, fwv, ftr
          WRITE(10, '(A)') TRIM(line)
      END DO

      ! Close the file
      CLOSE(10)
    END IF

    ! Calculate Vega flux for zero point calibration
    vega_flux = calculatevegaflux(vega_filepath, filter_wavelengths, filter_trans, &
                                 filter_name, make_sed)

    ! Calculate synthetic flux
    CALL calculate_syntheticflux(wavelengths, convolved_flux, synthetic_flux, &
                               filter_wavelengths, filter_trans)

    ! Calculate magnitude using Vega zero point
    IF (vega_flux > 0.0_dp) THEN
      calculate_synthetic = -2.5 * LOG10(synthetic_flux / vega_flux)
    ELSE
      PRINT *, "Error: Vega flux is zero, magnitude calculation is invalid."
      calculate_synthetic = HUGE(1.0_dp)
    END IF
    
    ! Clean up
    DEALLOCATE(convolved_flux)
  END FUNCTION calculate_synthetic

  !-----------------------------------------------------------------------
  ! Internal functions for synthetic photometry
  !-----------------------------------------------------------------------

  !****************************
  ! Convolve SED With Filter
  !****************************
  SUBROUTINE convolvesed(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: wavelengths, fluxes
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: filter_wavelengths, filter_trans
    REAL(dp), DIMENSION(:), ALLOCATABLE :: convolved_flux
    REAL(dp), DIMENSION(:), ALLOCATABLE :: interpolated_filter
    INTEGER :: n

    n = SIZE(wavelengths)

    ! Allocate arrays
    ALLOCATE(interpolated_filter(n))

    ! Interpolate the filter transmission onto the wavelengths array
    CALL interpolatearray(filter_wavelengths, filter_trans, wavelengths, interpolated_filter)

    ! Perform convolution (element-wise multiplication)
    convolved_flux = fluxes * interpolated_filter

    ! Deallocate temporary arrays
    DEALLOCATE(interpolated_filter)
  END SUBROUTINE convolvesed
  
  !****************************
  ! Calculate Synthetic Flux
  !****************************
  SUBROUTINE calculate_syntheticflux(wavelengths, fluxes, synthetic_flux, filter_wavelengths, filter_trans)
    REAL(dp), DIMENSION(:), INTENT(IN) :: wavelengths, fluxes
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: filter_wavelengths, filter_trans
    REAL(dp), INTENT(OUT) :: synthetic_flux
    INTEGER :: i
    REAL(dp) :: integrated_flux, integrated_filter

    ! Validate inputs
    DO i = 1, SIZE(wavelengths) - 1
      IF (wavelengths(i) <= 0.0 .OR. fluxes(i) < 0.0) THEN
        PRINT *, "synthetic Invalid input at index", i, ":", wavelengths(i), fluxes(i)
        STOP
      END IF
    END DO

    CALL rombergintegration(wavelengths, fluxes * wavelengths, integrated_flux)
    CALL rombergintegration(filter_wavelengths, filter_trans * filter_wavelengths, integrated_filter)

    ! Store the total flux
    IF (integrated_filter > 0.0) THEN
        synthetic_flux = integrated_flux / integrated_filter
    ELSE
        PRINT *, "Error: Integrated filter transmission is zero."
        synthetic_flux = -1.0_dp
        RETURN
    END IF
  END SUBROUTINE calculate_syntheticflux

  !****************************
  ! Calculate Bolometric Magnitude and Flux
  !****************************
  SUBROUTINE calculate_bolometric_phot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: wavelengths, fluxes
    REAL(dp), INTENT(OUT) :: bolometric_magnitude, bolometric_flux
    INTEGER :: i

    ! Validate inputs and replace invalid wavelengths with 0
    DO i = 1, SIZE(wavelengths) - 1
      IF (wavelengths(i) <= 0.0 .OR. fluxes(i) < 0.0) THEN
        PRINT *, "bolometric Invalid input at index", i, ":", wavelengths(i), fluxes(i)
        fluxes(i) = 0.0  ! Replace invalid wavelength with 0
      END IF
    END DO

    ! Call Romberg integration
    CALL rombergintegration(wavelengths, fluxes, bolometric_flux)

    ! Validate integration result
    IF (bolometric_flux <= 0.0) THEN
      PRINT *, "Error: Flux integration resulted in non-positive value."
      bolometric_magnitude = 99.0
      RETURN
    END IF

    ! Calculate bolometric magnitude
    IF (bolometric_flux <= 0.0) THEN
      PRINT *, "Error: Flux integration resulted in non-positive value."
      bolometric_magnitude = 99.0
      RETURN
    ELSE IF (bolometric_flux < 1.0E-10) THEN
      PRINT *, "Warning: Flux value is very small, precision might be affected."
    END IF

    bolometric_magnitude = fluxtomagnitude(bolometric_flux)
  END SUBROUTINE calculate_bolometric_phot
  
  !****************************
  ! Convert Flux to Magnitude
  !****************************
  REAL(dp) FUNCTION fluxtomagnitude(flux)
    REAL(dp), INTENT(IN) :: flux
    IF (flux <= 0.0) THEN
      PRINT *, "Error: Flux must be positive to calculate magnitude."
      fluxtomagnitude = 99.0  ! Return an error value
    ELSE
      fluxtomagnitude = -2.5 * LOG10(flux)
    END IF
  END FUNCTION fluxtomagnitude

  !****************************
  ! Calculate Vega Flux for Zero Point
  !****************************
  FUNCTION calculatevegaflux(vega_filepath, filt_wave, filt_trans, filter_name, make_sed) RESULT(vega_flux)
    CHARACTER(LEN=*), INTENT(IN) :: vega_filepath, filter_name
    CHARACTER(len = 100) :: output_csv
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: filt_wave, filt_trans
    REAL(dp) :: vega_flux
    REAL(dp) :: int_flux, int_filter
    REAL(dp), ALLOCATABLE :: vega_wave(:), vega_flux_arr(:), conv_flux(:)
    LOGICAL, INTENT(IN) :: make_sed
    INTEGER :: i, unit, max_size
    REAL(dp) :: wv, fl, cf, fwv, ftr
    INTEGER:: ierr
    CHARACTER(LEN=1000) :: line

    ! Load the Vega SED
    CALL loadvegased(vega_filepath, vega_wave, vega_flux_arr)

    ! Convolve the Vega SED with the filter transmission
    ALLOCATE(conv_flux(SIZE(vega_wave)))
    CALL convolvesed(vega_wave, vega_flux_arr, filt_wave, filt_trans, conv_flux)

    ! Integrate the convolved Vega SED and the filter transmission
    CALL rombergintegration(vega_wave, vega_wave*conv_flux, int_flux)
    CALL rombergintegration(filt_wave, filt_wave*filt_trans, int_filter)

    IF (int_filter > 0.0_dp) THEN
      vega_flux = int_flux / int_filter
    ELSE
      vega_flux = -1.0_dp
    END IF

    ! Write Vega SED to CSV if requested
    IF (make_sed) THEN
      ! Determine the maximum size among all arrays
      max_size = MAX(SIZE(vega_wave), SIZE(vega_flux_arr), SIZE(conv_flux), &
                     SIZE(filt_wave), SIZE(filt_trans))

      output_csv = 'LOGS/SED/VEGA_' //TRIM(remove_dat(filter_name)) // '_SED.csv'

      ! Open the CSV file for writing
      OPEN(UNIT=10, FILE=output_csv, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
      IF (ierr /= 0) THEN
        PRINT *, "Error opening file for writing"
        RETURN
      END IF

      WRITE(10, '(A)') "wavelengths,fluxes,convolved_flux,filter_wavelengths,filter_trans"

      ! Loop through data and safely write values, ensuring no out-of-bounds errors
      DO i = 1, max_size
        ! Initialize values to zero in case they are out of bounds
        wv = 0.0_dp
        fl = 0.0_dp
        cf = 0.0_dp
        fwv = 0.0_dp
        ftr = 0.0_dp

        ! Assign actual values only if within valid indices
        IF (i <= SIZE(vega_wave)) wv = vega_wave(i)
        IF (i <= SIZE(vega_flux_arr)) fl = vega_flux_arr(i)
        IF (i <= SIZE(conv_flux)) cf = conv_flux(i)
        IF (i <= SIZE(filt_wave)) fwv = filt_wave(i)
        IF (i <= SIZE(filt_trans)) ftr = filt_trans(i)

        ! Write the formatted output
        WRITE(line, '(ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6, ",", ES14.6)') &
            wv, fl, cf, fwv, ftr
        WRITE(10, '(A)') TRIM(line)
      END DO

      ! Close the file
      CLOSE(10)
    END IF

    ! Clean up
    DEALLOCATE(conv_flux, vega_wave, vega_flux_arr)
  END FUNCTION calculatevegaflux




  !-----------------------------------------------------------------------
  ! Bolometric correction interface (stub implementations)
  !-----------------------------------------------------------------------
  
  ! These functions are defined as stubs since they appear to be
  ! placeholder interfaces for future implementation
  
  real(dp) function get_bc_by_name(name, log_Teff, log_g, M_div_h, ierr)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: log_Teff  ! log10 of surface temp
    real(dp), intent(in) :: log_g  ! log_10 of surface gravity
    real(dp), intent(in) :: M_div_h  ! [M/H]
    integer, intent(inout) :: ierr
    
    get_bc_by_name = -99.9d0
    ierr = 0
  end function get_bc_by_name

  real(dp) function get_bc_by_id(id, log_Teff, log_g, M_div_h, ierr)
    integer, intent(in) :: id
    real(dp), intent(in) :: log_Teff  ! log10 of surface temp
    real(dp), intent(in) :: log_g  ! log_10 of surface gravity
    real(dp), intent(in) :: M_div_h  ! [M/H]
    integer, intent(inout) :: ierr
    
    get_bc_by_id = -99.9d0
    ierr = 0
  end function get_bc_by_id

  integer function get_bc_id_by_name(name, ierr)
    character(len=*), intent(in) :: name
    integer, intent(inout) :: ierr
    
    get_bc_id_by_name = -1
    ierr = 0
  end function get_bc_id_by_name

  character(len=strlen) function get_bc_name_by_id(id, ierr)
    integer, intent(in) :: id
    integer, intent(inout) :: ierr
    
    get_bc_name_by_id = ''
    ierr = 0
  end function get_bc_name_by_id

  real(dp) function get_abs_bolometric_mag(lum)
    use const_def
    real(dp), intent(in) :: lum  ! Luminosity in lsun units
    
    get_abs_bolometric_mag = -99.9d0
  end function get_abs_bolometric_mag

  real(dp) function get_abs_mag_by_name(name, log_Teff, log_g, M_div_h, lum, ierr)
    character(len=*) :: name
    real(dp), intent(in) :: log_Teff  ! log10 of surface temp
    real(dp), intent(in) :: M_div_h  ! [M/H]
    real(dp), intent(in) :: log_g  ! log_10 of surface gravity
    real(dp), intent(in) :: lum  ! Luminosity in lsun units
    integer, intent(inout) :: ierr
    
    ierr = 0
    get_abs_mag_by_name = -99.9d0
  end function get_abs_mag_by_name

  real(dp) function get_abs_mag_by_id(id, log_Teff, log_g, M_div_h, lum, ierr)
    integer, intent(in) :: id
    real(dp), intent(in) :: log_Teff  ! log10 of surface temp
    real(dp), intent(in) :: log_g  ! log_10 of surface gravity
    real(dp), intent(in) :: M_div_h  ! [M/H]
    real(dp), intent(in) :: lum  ! Luminosity in lsun units
    integer, intent(inout) :: ierr
    
    ierr = 0
    get_abs_mag_by_id = -99.9d0
  end function get_abs_mag_by_id

  subroutine get_bcs_all(log_Teff, log_g, M_div_h, results, ierr)
    real(dp), intent(in) :: log_Teff  ! log10 of surface temp
    real(dp), intent(in) :: M_div_h  ! [M/H]
    real(dp), dimension(:), intent(out) :: results
    real(dp), intent(in) :: log_g
    integer, intent(inout) :: ierr
    
    ierr = 0
    results(:) = -99.d0
  end subroutine get_bcs_all

  real(dp) function get_lum_band_by_name(name, log_Teff, log_g, M_div_h, lum, ierr)
    character(len=*) :: name
    real(dp), intent(in) :: log_Teff  ! log10 of surface temp
    real(dp), intent(in) :: M_div_h  ! [M/H]
    real(dp), intent(in) :: log_g  ! log_10 of surface gravity
    real(dp), intent(in) :: lum  ! Total luminosity in lsun units
    integer, intent(inout) :: ierr
    
    ierr = 0
    get_lum_band_by_name = -99.d0
  end function get_lum_band_by_name

  real(dp) function get_lum_band_by_id(id, log_Teff, log_g, M_div_h, lum, ierr)
    integer, intent(in) :: id
    real(dp), intent(in) :: log_Teff  ! log10 of surface temp
    real(dp), intent(in) :: log_g  ! log_10 of surface gravity
    real(dp), intent(in) :: M_div_h  ! [M/H]
    real(dp), intent(in) :: lum  ! Total luminosity in lsun units
    integer, intent(inout) :: ierr
    
    ierr = 0
    get_lum_band_by_id = -99.d0
  end function get_lum_band_by_id

  !-----------------------------------------------------------------------
  ! Test suite interface (stub implementations)
  !-----------------------------------------------------------------------
  
  subroutine test_suite_startup(restart, ierr)
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    ierr = 0
  end subroutine test_suite_startup

  subroutine test_suite_after_evolve(ierr)
    integer, intent(out) :: ierr
    ierr = 0
  end subroutine test_suite_after_evolve

end module colors_lib