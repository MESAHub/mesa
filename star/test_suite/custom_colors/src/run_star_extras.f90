
! ***********************************************************************
!
!   Copyright (C) 2017-2019  Rob Farmer & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
 

 
 
module run_star_extras

  use star_lib
  use star_def
  use const_def
  use math_lib
  use auto_diff
  use colors_lib

  implicit none


!DEFINE ALL GLOCBAL VARIABLE HERE



  include "test_suite_extras_def.inc"

  ! these routines are called by the standard run_star check_model
  contains

  include "test_suite_extras.inc"


  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
       print *, "Extras startup routine"
           
    call process_color_files(id, ierr)
    s% extras_startup => extras_startup
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns

    print *, "Sellar atmosphere:", s% x_character_ctrl(1)
    print *, "Instrument:", s% x_character_ctrl(2)         

  end subroutine extras_controls

                
      
  

!###########################################################
!## THINGS I HAVE NOT TOUCHED
!###########################################################
  
  subroutine process_color_files(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type(star_info), pointer :: s
    integer :: i

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

  end subroutine process_color_files


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
     integer, intent(in) :: id
     integer, intent(out) :: ierr
     type (star_info), pointer :: s
     real(dp) :: dt
     ierr = 0
     call star_ptr(id, s, ierr)
     if (ierr /= 0) return

     write(*,'(a)') 'finished custom colors'
     
     call test_suite_after_evolve(s, ierr)
     
  end subroutine extras_after_evolve


  ! returns either keep_going, retry, or terminate.
  integer function extras_check_model(id)
     integer, intent(in) :: id
     integer :: ierr
     type (star_info), pointer :: s
     ierr = 0
     call star_ptr(id, s, ierr)
     if (ierr /= 0) return
     extras_check_model = keep_going         
  end function extras_check_model


  INTEGER FUNCTION how_many_extra_profile_columns(id)
     USE star_def, ONLY: star_info
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: id

     INTEGER :: ierr
     TYPE(star_info), POINTER :: s

     ierr = 0
     CALL star_ptr(id, s, ierr)
     IF (ierr /= 0) RETURN

     how_many_extra_profile_columns = 0
  END FUNCTION how_many_extra_profile_columns


  SUBROUTINE data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
     USE star_def, ONLY: star_info, maxlen_profile_column_name
     USE const_def, ONLY: DP
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: id, n, nz
     CHARACTER(LEN=maxlen_profile_column_name) :: names(n)
     REAL(DP) :: vals(nz, n)
     INTEGER, INTENT(OUT) :: ierr

     TYPE(star_info), POINTER :: s

     ierr = 0
     CALL star_ptr(id, s, ierr)
     IF (ierr /= 0) RETURN

  END SUBROUTINE data_for_extra_profile_columns


  ! Returns either keep_going, retry, or terminate
  INTEGER FUNCTION extras_finish_step(id)
     USE chem_def
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: id

     INTEGER :: ierr
     TYPE(star_info), POINTER :: s

     ierr = 0
     CALL star_ptr(id, s, ierr)
     IF (ierr /= 0) RETURN

     extras_finish_step = keep_going
  END FUNCTION extras_finish_step






!###########################################################
!## MESA STUFF
!###########################################################
  
  !FUNCTIONS FOR OPENING LOOKUP FILE AND FINDING THE NUMBER OF FILES AND THIER FILE PATHS
  integer function how_many_extra_history_columns(id)
      ! Determines how many extra history columns are added based on a file
      integer, intent(in) :: id
      integer :: ierr, n
      character(len=100), allocatable :: strings(:)
      type(star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) then
          how_many_extra_history_columns = 0
          return
      end if

      ! Read strings from the file
      call read_strings_from_file(strings, n, id)

      ! Number of columns is the size of the strings array
      how_many_extra_history_columns = n + 2

      !print *, "This many columns added to history file:", n

      if (allocated(strings)) deallocate(strings)
  end function how_many_extra_history_columns


  function basename(path) result(base)
      ! Extracts the base name from a given file path
      character(len=*), intent(in) :: path
      character(len=512) :: base
      integer :: last_slash

      ! Find the position of the last slash
      last_slash = len_trim(path)
      do while (last_slash > 0 .and. path(last_slash:last_slash) /= '/')
          last_slash = last_slash - 1
      end do

      ! Extract the base name
      base = path(last_slash+1:)
  end function basename

  function remove_dat(path) result(base)
      ! Extracts the portion of the string after the first dot
      character(len=*), intent(in) :: path
      character(len=512) :: base
      integer :: first_dot

      ! Find the position of the first dot
      first_dot = 0
      do while (first_dot < len_trim(path) .and. path(first_dot+1:first_dot+1) /= '.')
          first_dot = first_dot + 1
      end do

      ! Check if an dot was found
      if (first_dot < len_trim(path)) then
          ! Extract the part after the dot
          base = path(:first_dot)
      else
          ! No dot found, return the input string
          base = path
      end if
  end function remove_dat


  subroutine read_strings_from_file(strings, n, id)
      ! Reads strings from a file into an allocatable array
      implicit none
      integer, intent(in) :: id
      character(len=512) :: filename
      character(len=100), allocatable :: strings(:)
      integer, intent(out) :: n
      integer :: unit, i, status
      character(len=100) :: line
      integer :: ierr
      type(star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! Construct the filename
      filename = trim(s%x_character_ctrl(2)) // "/" // trim(basename(s%x_character_ctrl(2)))

      ! Initialize
      n = 0

      ! Open the file
      unit = 10
      open(unit, file=filename, status='old', action='read', iostat=status)
      if (status /= 0) then
          print *, "Error: Could not open file", filename
          stop
      end if

      ! Count lines in the file to determine the size of the array
      do
          read(unit, '(A)', iostat=status) line
          if (status /= 0) exit
          n = n + 1 ! for bolometric correctionms
      end do
      rewind(unit)

      ! Allocate the array and read the strings
      if (allocated(strings)) deallocate(strings)
      allocate(strings(n))
      do i = 1, n
          read(unit, '(A)') strings(i)
      end do

      close(unit)
  end subroutine read_strings_from_file



  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      ! Populates data for the extra history columns
      integer, intent(in) :: id, n
      integer, intent(out) :: ierr
      character(len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer :: i, num_strings
      character(len=100), allocatable :: array_of_strings(:)
      real(dp) :: teff, log_g, metallicity, R, d,  bolometric_magnitude, bolometric_flux
      character(len=256) :: sed_filepath, filter_filepath, filter_name, filter_dir, vega_filepath
      real(dp), dimension(:), allocatable :: wavelengths, fluxes, filter_wavelengths, filter_trans
      logical :: make_sed
      
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! Extract input parameters
      teff = s%T(1)
      log_g = LOG10(s%grav(1))
      R = s%R(1)! * 1d3
      metallicity = s%job%extras_rpar(1)
      d = s%job%extras_rpar(2)      

      sed_filepath = s%x_character_ctrl(1)
      filter_dir = s%x_character_ctrl(2)
      vega_filepath = s%x_character_ctrl(3)
      make_sed = trim(adjustl(s%x_character_ctrl(4))) == 'true'
      
      ! Read filters from file
      if (allocated(array_of_strings)) deallocate(array_of_strings)
      allocate(array_of_strings(n))
      call read_strings_from_file(array_of_strings, num_strings, id)

      !PRINT *, "################################################"
      
      ! Compute bolometric values
      CALL CalculateBolometric(teff, log_g, metallicity, R, d,  bolometric_magnitude, bolometric_flux, wavelengths, fluxes, sed_filepath)
      names(1) = "Mag_bol"
      vals(1) = bolometric_magnitude
      names(2) = "Flux_bol"
      vals(2) = bolometric_flux
      
      ! Populate history columns
      if (allocated(array_of_strings)) then
          do i = 3, how_many_extra_history_columns(id)
              filter_name = "Unknown"
              if (i <= num_strings + 2) filter_name = trim(remove_dat(array_of_strings(i - 2)))
              names(i) = filter_name
              filter_filepath = trim(filter_dir) // "/" // array_of_strings(i - 2)
              
              if (teff >= 0 .and. log_g >= 0 .and. metallicity >= 0) then
                  vals(i) = CalculateSynthetic(teff, log_g, metallicity, ierr, wavelengths, fluxes, filter_wavelengths, filter_trans, filter_filepath, vega_filepath, array_of_strings(i - 2), make_sed)
                  if (ierr /= 0) vals(i) = -1.0_dp
              else
                  vals(i) = -1.0_dp
                  ierr = 1
              end if
              !PRINT *, names(i), vals(i)
          end do
      else
          ierr = 1 ! Indicate an error if array_of_strings is not allocated
      end if
      
      if (allocated(array_of_strings)) deallocate(array_of_strings)
  end subroutine data_for_extra_history_columns




!###########################################################
!## CUSTOM COLOURS
!###########################################################

!****************************
!Calculate Bolometric Photometry Using Multiple SEDs
!****************************

  SUBROUTINE CalculateBolometric(teff, log_g, metallicity, R, d, bolometric_magnitude, bolometric_flux, wavelengths, fluxes, sed_filepath)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: teff, log_g, metallicity, R, d
    CHARACTER(LEN=*), INTENT(IN) :: sed_filepath
    REAL(DP), INTENT(OUT) :: bolometric_magnitude, bolometric_flux

    REAL (8), ALLOCATABLE :: lu_logg(:), lu_meta(:), lu_teff(:)
    CHARACTER(LEN=100), ALLOCATABLE :: file_names(:)
    REAL, DIMENSION(:,:), ALLOCATABLE :: lookup_table
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, fluxes
    CHARACTER(LEN=256) :: lookup_file

    lookup_file = TRIM(sed_filepath) // '/lookup_table.csv'

    ! Call to load the lookup table
    CALL LoadLookupTable(lookup_file, lookup_table, file_names, lu_logg, lu_meta, lu_teff)
    !print *, 'logg', lu_logg
    !print *,  'meta', lu_meta
    !print *, 'teff', lu_teff
    ! Interpolate Spectral Energy Distribution
    !CALL ConstructSED_Interpolated(teff, log_g, metallicity, R, d, file_names, lu_teff, lu_logg, lu_meta, sed_filepath, wavelengths, fluxes)
    CALL ConstructSED(teff, log_g, metallicity, R, d, file_names, lu_teff, lu_logg, lu_meta, sed_filepath, wavelengths, fluxes)    

    ! Calculate bolometric flux and magnitude
    CALL CalculateBolometricPhot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
  END SUBROUTINE CalculateBolometric



!****************************
!Construct SED With Combination of SEDs
!****************************

SUBROUTINE ConstructSED(teff, log_g, metallicity, R, d, file_names, lu_teff, lu_logg, lu_meta, stellar_model_dir, wavelengths, fluxes)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: teff, log_g, metallicity, R, d
  REAL(8), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
  CHARACTER(LEN=*), INTENT(IN) :: stellar_model_dir
  CHARACTER(LEN=100), INTENT(IN) :: file_names(:)
  REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, fluxes

  INTEGER, DIMENSION(4) :: closest_indices
  REAL(DP), DIMENSION(:), ALLOCATABLE :: temp_wavelengths, temp_flux, common_wavelengths
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: model_fluxes
  REAL(DP), DIMENSION(4) :: weights, distances
  INTEGER :: i, n_points
  REAL(DP) :: sum_weights
  REAL(DP), DIMENSION(:), ALLOCATABLE :: diluted_flux

  ! Get the four closest stellar models
  CALL GetClosestStellarModels(teff, log_g, metallicity, lu_teff, lu_logg, lu_meta, closest_indices)

  ! Load the first SED to define the wavelength grid
  CALL LoadSED(TRIM(stellar_model_dir) // TRIM(file_names(closest_indices(1))), closest_indices(1), temp_wavelengths, temp_flux)
  n_points = SIZE(temp_wavelengths)
  ALLOCATE(common_wavelengths(n_points))
  common_wavelengths = temp_wavelengths

  ! Allocate flux array for the models (4 models, n_points each)
  ALLOCATE(model_fluxes(4, n_points))
  CALL InterpolateArray(temp_wavelengths, temp_flux, common_wavelengths, model_fluxes(1, :))

  ! Load and interpolate remaining SEDs
  DO i = 2, 4
    CALL LoadSED(TRIM(stellar_model_dir) // TRIM(file_names(closest_indices(i))), closest_indices(i), temp_wavelengths, temp_flux)
    CALL InterpolateArray(temp_wavelengths, temp_flux, common_wavelengths, model_fluxes(i, :))
  END DO

  ! Compute distances and weights for the four models
  DO i = 1, 4
    distances(i) = SQRT((lu_teff(closest_indices(i)) - teff)**2 + &
                        (lu_logg(closest_indices(i)) - log_g)**2 + &
                        (lu_meta(closest_indices(i)) - metallicity)**2)
    IF (distances(i) == 0.0) distances(i) = 1.0E-10  ! Prevent division by zero
    weights(i) = 1.0 / distances(i)
  END DO

  ! Normalize weights
  sum_weights = SUM(weights)
  weights = weights / sum_weights

  ! Allocate output arrays
  ALLOCATE(wavelengths(n_points), fluxes(n_points))
  wavelengths = common_wavelengths
  fluxes = 0.0

  ! Perform weighted combination of the model fluxes (still at the stellar surface)
  DO i = 1, 4
    fluxes = fluxes + weights(i) * model_fluxes(i, :)
  END DO

  ! Now, apply the dilution factor (R/d)^2 to convert the surface flux density
  ! into the observed flux density at Earth.
  ALLOCATE(diluted_flux(n_points))
  CALL dilute_flux(fluxes, R, d, diluted_flux)
  fluxes = diluted_flux

  ! Deallocate temporary arrays
  DEALLOCATE(temp_wavelengths, temp_flux, common_wavelengths, diluted_flux)

END SUBROUTINE ConstructSED





SUBROUTINE dilute_flux(surface_flux, R, d, calibrated_flux)
  IMPLICIT NONE
  ! Define the double precision kind if not already defined
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  
  ! Input: surface_flux is an array of flux values at the stellar surface
  REAL(DP), INTENT(IN)  :: surface_flux(:)
  REAL(DP), INTENT(IN)  :: R, d  ! R = stellar radius, d = distance (both in the same units, e.g., cm)
  
  ! Output: calibrated_flux will be the flux observed at Earth
  REAL(DP), INTENT(OUT) :: calibrated_flux(:)
  
  ! Check that the output array has the same size as the input
  IF (SIZE(calibrated_flux) /= SIZE(surface_flux)) THEN
    PRINT *, "Error in dilute_flux: Output array must have the same size as input array."
    STOP 1
  END IF
  
  ! Apply the dilution factor (R/d)^2 to each element
  calibrated_flux = surface_flux * ( (R / d)**2 )
  
END SUBROUTINE dilute_flux






!****************************
!Identify The Four Closest Stellar Models
!****************************

SUBROUTINE GetClosestStellarModels(teff, log_g, metallicity, lu_teff, lu_logg, lu_meta, closest_indices)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: teff, log_g, metallicity
  REAL(8), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
  INTEGER, DIMENSION(4), INTENT(OUT) :: closest_indices

  INTEGER :: i, n, j
  REAL(DP) :: distance, norm_teff, norm_logg, norm_meta
  REAL(DP), DIMENSION(:), ALLOCATABLE :: scaled_lu_teff, scaled_lu_logg, scaled_lu_meta
  REAL(DP), DIMENSION(4) :: min_distances
  INTEGER, DIMENSION(4) :: indices
  REAL(DP) :: teff_min, teff_max, logg_min, logg_max, meta_min, meta_max, teff_dist, logg_dist, meta_dist

  n = SIZE(lu_teff)
  min_distances = HUGE(1.0)
  indices = -1

  ! Find min and max for normalization
  teff_min = MINVAL(lu_teff)
  teff_max = MAXVAL(lu_teff)
  logg_min = MINVAL(lu_logg)
  logg_max = MAXVAL(lu_logg)
  meta_min = MINVAL(lu_meta)
  meta_max = MAXVAL(lu_meta)

  ! Allocate and scale lookup table values
  ALLOCATE(scaled_lu_teff(n), scaled_lu_logg(n), scaled_lu_meta(n))

  IF (teff_max - teff_min > 0.00) THEN
    scaled_lu_teff = (lu_teff - teff_min) / (teff_max - teff_min)
  END IF

  IF (logg_max - logg_min > 0.00) THEN
    scaled_lu_logg = (lu_logg - logg_min) / (logg_max - logg_min)
  END IF

  IF (meta_max - meta_min > 0.00) THEN    
    scaled_lu_meta = (lu_meta - meta_min) / (meta_max - meta_min)
  END IF

  ! Normalize input parameters
  norm_teff = (teff - teff_min) / (teff_max - teff_min)
  norm_logg = (log_g - logg_min) / (logg_max - logg_min)
  norm_meta = (metallicity - meta_min) / (meta_max - meta_min)

  ! Debug: !PRINT normalized input parameters
  !PRINT *, "Normalized parameters for target:"
  !PRINT *, "  teff = ", teff, "  logg = ", log_g, "  meta = ", metallicity, n

  ! Find closest models
  DO i = 1, n

    teff_dist = 0.0
    logg_dist = 0.0
    meta_dist = 0.0

    IF (teff_max - teff_min > 0.00) THEN
      teff_dist = scaled_lu_teff(i) - norm_teff
    END IF

    IF (logg_max - logg_min > 0.00) THEN
      logg_dist = scaled_lu_logg(i) - norm_logg
    END IF

    IF (meta_max - meta_min > 0.00) THEN    
      meta_dist = scaled_lu_meta(i) - norm_meta
    END IF


    distance = SQRT(teff_dist**2 + logg_dist**2 + meta_dist**2)

    ! Check if this distance is smaller than any in the current top 4
    !PRINT *, distance
    !PRINT *, scaled_lu_teff(i)
    !PRINT *, norm_teff
    !PRINT *, scaled_lu_logg(i)
    !PRINT *, norm_logg
    !PRINT *, scaled_lu_meta(i)
    !PRINT *, norm_meta

    DO j = 1, 4
      IF (distance < min_distances(j)) THEN
        ! Shift larger distances down
        IF (j < 4) THEN
          min_distances(j+1:4) = min_distances(j:3)
          indices(j+1:4) = indices(j:3)
        END IF
        min_distances(j) = distance
        indices(j) = i
        EXIT
      END IF
    END DO
  END DO

  closest_indices = indices
  ! Deallocate arrays
  DEALLOCATE(scaled_lu_teff, scaled_lu_logg, scaled_lu_meta)
END SUBROUTINE GetClosestStellarModels




!****************************
!Calculate Bolometric Magnitude and Flux
!****************************

  SUBROUTINE CalculateBolometricPhot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: wavelengths, fluxes
    REAL(DP), INTENT(OUT) :: bolometric_magnitude, bolometric_flux
    INTEGER :: i

    ! Validate inputs and replace invalid wavelengths with 0
    DO i = 1, SIZE(wavelengths) - 1
      IF (wavelengths(i) <= 0.0 .OR. fluxes(i) < 0.0) THEN
        PRINT *, "bolometric Invalid input at index", i, ":", wavelengths(i), fluxes(i)
        fluxes(i) = 0.0  ! Replace invalid wavelength with 0
      END IF
    END DO


    ! Perform trapezoidal integration
    ! Debug: Print the first few wavelengths and flux values
    !PRINT *, "Wavelengths (first 5):", wavelengths(1:MIN(5, SIZE(wavelengths)))
    !PRINT *, "Fluxes (first 5):", fluxes(1:MIN(5, SIZE(fluxes)))

    ! Call trapezoidal integration
    CALL RombergIntegration(wavelengths, fluxes, bolometric_flux)

    ! Debug: Check the integration result
    !PRINT *, "Integrated Flux:", bolometric_flux

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

  bolometric_magnitude = FluxToMagnitude(bolometric_flux)

  END SUBROUTINE CalculateBolometricPhot












!###########################################################
!## Synthetic Photometry
!###########################################################

!****************************
!Calculate Synthetic Photometry Using SED and Filter
!****************************



REAL(DP) FUNCTION CalculateSynthetic(temperature, gravity, metallicity, ierr, wavelengths, fluxes, filter_wavelengths, filter_trans, filter_filepath, vega_filepath, filter_name, make_sed)
    IMPLICIT NONE

    ! Input arguments
    REAL(DP), INTENT(IN) :: temperature, gravity, metallicity
    CHARACTER(LEN=*), INTENT(IN) :: filter_filepath, filter_name, vega_filepath
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=1000) :: line

    REAL(DP), DIMENSION(:), INTENT(INOUT) :: wavelengths, fluxes
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: filter_wavelengths, filter_trans
    LOGICAL, INTENT(IN) :: make_sed
    LOGICAL :: dir_exists
    ! Local variables
    REAL(DP), DIMENSION(:), ALLOCATABLE :: convolved_flux, interpolated_filter
    CHARACTER(LEN=100) :: csv_file
    REAL(DP) :: synthetic_magnitude, synthetic_flux, vega_flux
    INTEGER :: max_size, i
    REAL(DP) :: magnitude
    REAL(DP) :: wv, fl, cf, fwv, ftr

    csv_file =  'LOGS/SED/' //TRIM(remove_dat(filter_name)) // '_SED.csv'
    ! Initialize error flag
    ierr = 0

    ! Load filter data
    CALL LoadFilter(filter_filepath, filter_wavelengths, filter_trans)

    ! Check for invalid gravity input
    IF (gravity <= 0.0_DP) THEN
        ierr = 1
        CalculateSynthetic = -1.0_DP
        RETURN
    END IF

    ! Allocate interpolated_filter if not already allocated
    IF (.NOT. ALLOCATED(interpolated_filter)) THEN
        ALLOCATE(interpolated_filter(SIZE(wavelengths)))
        interpolated_filter = 0.0_DP
    END IF

      ! Perform SED convolution
  ALLOCATE(convolved_flux(SIZE(wavelengths)))
  CALL ConvolveSED(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)

  IF (make_sed) THEN
! Determine the maximum size among all arrays
max_size = MAX(SIZE(wavelengths), SIZE(filter_wavelengths), SIZE(fluxes), SIZE(convolved_flux), SIZE(filter_trans))

! Open the CSV file for writing
OPEN(UNIT=10, FILE=csv_file, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
IF (ierr /= 0) THEN
    PRINT *, "Error opening file for writing"
    STOP
END IF

! Write headers to the CSV file
WRITE(10, '(A)') "wavelengths,fluxes,convolved_flux,filter_wavelengths,filter_trans"

! Loop through data and safely write values, ensuring no out-of-bounds errors
DO i = 1, max_size
    ! Initialize values to zero in case they are out of bounds
    wv = 0.0_DP
    fl = 0.0_DP
    cf = 0.0_DP
    fwv = 0.0_DP
    ftr = 0.0_DP

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

    ! Inform the user of successful writing
    !PRINT *, "Data written to ", csv_file
    vega_flux = CalculateVegaFlux(vega_filepath, filter_wavelengths, filter_trans, filter_name, make_sed)

    ! Calculate synthetic flux and magnitude
    CALL CalculateSyntheticFlux(wavelengths, convolved_flux, synthetic_flux, filter_wavelengths, filter_trans)

    !PRINT *, "VEGA zero point:", vega_flux

    IF (vega_flux > 0.0_DP) THEN
      CalculateSynthetic = -2.5 * LOG10(synthetic_flux / vega_flux)
    ELSE
      PRINT *, "Error: Vega flux is zero, magnitude calculation is invalid."
      CalculateSynthetic = HUGE(1.0_DP)
    END IF

END FUNCTION CalculateSynthetic




!****************************
!Convolve SED With Filter
!****************************

  SUBROUTINE ConvolveSED(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: wavelengths, fluxes
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: filter_wavelengths, filter_trans
    REAL(DP), DIMENSION(:), ALLOCATABLE :: convolved_flux
    REAL(DP), DIMENSION(:), ALLOCATABLE :: interpolated_filter
    INTEGER :: n

    n = SIZE(wavelengths)

    ! Allocate arrays
    ALLOCATE(interpolated_filter(n))
    !ALLOCATE(convolved_flux(n))

    ! Interpolate the filter transmission onto the wavelengths array
    CALL InterpolateArray(filter_wavelengths, filter_trans, wavelengths, interpolated_filter)

    ! Perform convolution (element-wise multiplication)
    convolved_flux = fluxes * interpolated_filter

    ! Deallocate arrays (optional, depending on context)
    DEALLOCATE(interpolated_filter)
  END SUBROUTINE ConvolveSED



!****************************
!Calculate Synthetic Flux and Magnitude
!****************************
  SUBROUTINE CalculateSyntheticFlux(wavelengths, fluxes, synthetic_flux, filter_wavelengths, filter_trans)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: wavelengths, fluxes
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: filter_wavelengths, filter_trans
    REAL(DP), INTENT(OUT) :: synthetic_flux
    INTEGER :: i
    REAL(DP) :: integrated_flux, integrated_filter
    CHARACTER(LEN=256) :: vega_filepath



    ! Validate inputs
    DO i = 1, SIZE(wavelengths) - 1
      IF (wavelengths(i) <= 0.0 .OR. fluxes(i) < 0.0) THEN
        PRINT *, "synthetic Invalid input at index", i, ":", wavelengths(i), fluxes(i)
        STOP
      END IF
    END DO

    CALL RombergIntegration(wavelengths, fluxes* wavelengths, integrated_flux)
    CALL RombergIntegration(filter_wavelengths, filter_trans * filter_wavelengths, integrated_filter)

    ! Store the total flux
    IF (integrated_filter > 0.0) THEN
        synthetic_flux = integrated_flux / integrated_filter
    ELSE
        PRINT *, "Error: Integrated filter transmission is zero."
        synthetic_flux = -1.0_DP
        RETURN
    END IF

  END SUBROUTINE CalculateSyntheticFlux



  REAL(DP) FUNCTION FluxToMagnitude(flux)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: flux
    !print *, 'flux:', flux
    IF (flux <= 0.0) THEN
      PRINT *, "Error: Flux must be positive to calculate magnitude."
      FluxToMagnitude = 99.0  ! Return an error value
    ELSE
      FluxToMagnitude = -2.5 * LOG10(flux) 
    END IF
  END FUNCTION FluxToMagnitude






FUNCTION CalculateVegaFlux(vega_filepath, filt_wave, filt_trans, filter_name, make_sed) RESULT(vega_flux)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: vega_filepath, filter_name
  CHARACTER(len = 100) :: output_csv
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: filt_wave, filt_trans
  REAL(DP) :: vega_flux
  REAL(DP) :: int_flux, int_filter
  REAL(DP), ALLOCATABLE :: vega_wave(:), vega_flux_arr(:), conv_flux(:)
  LOGICAL, INTENT(IN) :: make_sed
  INTEGER :: i, unit, max_size
  REAL(DP) :: wv, fl, cf, fwv, ftr
  INTEGER:: ierr
  CHARACTER(LEN=1000) :: line

  ! Load the Vega SED using the custom routine.
  CALL LoadVegaSED(vega_filepath, vega_wave, vega_flux_arr)

  ! Convolve the Vega SED with the filter transmission.
  CALL ConvolveSED(vega_wave, vega_flux_arr, filt_wave, filt_trans, conv_flux)

  ! Integrate the convolved Vega SED and the filter transmission.
  CALL RombergIntegration(vega_wave, vega_wave*conv_flux, int_flux)
  CALL RombergIntegration(filt_wave, filt_wave*filt_trans, int_filter)

  IF (int_filter > 0.0_DP) THEN
    vega_flux = int_flux / int_filter
  ELSE
    vega_flux = -1.0_DP
  END IF




  IF (make_sed) THEN
    ! Determine the maximum size among all arrays
    max_size = MAX(SIZE(vega_wave), SIZE(vega_flux_arr), SIZE(conv_flux), SIZE(filt_wave), SIZE(filt_trans))


    output_csv = 'LOGS/SED/VEGA_' //TRIM(remove_dat(filter_name)) // '_SED.csv'

    ! Open the CSV file for writing
    OPEN(UNIT=10, FILE=output_csv, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
    IF (ierr /= 0) THEN
        PRINT *, "Error opening file for writing"
        STOP
    END IF

    WRITE(10, '(A)') "wavelengths,fluxes,convolved_flux,filter_wavelengths,filter_trans"


    ! Loop through data and safely write values, ensuring no out-of-bounds errors
    DO i = 1, max_size
        ! Initialize values to zero in case they are out of bounds
        wv = 0.0_DP
        fl = 0.0_DP
        cf = 0.0_DP
        fwv = 0.0_DP
        ftr = 0.0_DP

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





  DEALLOCATE(conv_flux, vega_wave, vega_flux_arr)
END FUNCTION CalculateVegaFlux










SUBROUTINE LoadVegaSED(filepath, wavelengths, flux)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: filepath
  REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, flux
  CHARACTER(LEN=512) :: line
  INTEGER :: unit, n_rows, status, i
  REAL(DP) :: temp_wave, temp_flux

  unit = 20
  OPEN(unit, FILE=TRIM(filepath), STATUS='OLD', ACTION='READ', IOSTAT=status)
  IF (status /= 0) THEN
    PRINT *, "Error: Could not open Vega SED file ", TRIM(filepath)
    STOP
  END IF

  ! Skip header line.
  READ(unit, '(A)', IOSTAT=status) line
  IF (status /= 0) THEN
    PRINT *, "Error: Could not read header from Vega SED file ", TRIM(filepath)
    STOP
  END IF

  ! Count the number of data lines.
  n_rows = 0
  DO
    READ(unit, '(A)', IOSTAT=status) line
    IF (status /= 0) EXIT
    n_rows = n_rows + 1
  END DO

  REWIND(unit)
  READ(unit, '(A)', IOSTAT=status) line  ! Skip header again

  ALLOCATE(wavelengths(n_rows))
  ALLOCATE(flux(n_rows))
  
  i = 0
  DO
    READ(unit, *, IOSTAT=status) temp_wave, temp_flux  ! Ignore any extra columns.
    IF (status /= 0) EXIT
    i = i + 1
    wavelengths(i) = temp_wave
    flux(i) = temp_flux
  END DO

  CLOSE(unit)
END SUBROUTINE LoadVegaSED

                          


!###########################################################
!## FILE IO
!###########################################################

!****************************
!Load SED File
!****************************

  SUBROUTINE LoadSED(directory, index, wavelengths, flux)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: directory
    INTEGER, INTENT(IN) :: index
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, flux

    CHARACTER(LEN=512) :: line
    INTEGER :: unit, n_rows, status, i
    REAL :: temp_wavelength, temp_flux

    ! Open the file
    unit = 20
    OPEN(unit, FILE=TRIM(directory), STATUS='OLD', ACTION='READ', IOSTAT=status)
    IF (status /= 0) THEN
      PRINT *, "Error: Could not open file ", TRIM(directory)
      STOP
    END IF

    ! Skip header lines
    DO
      READ(unit, '(A)', IOSTAT=status) line
      IF (status /= 0) THEN
        PRINT *, "Error: Could not read the file", TRIM(directory)
        STOP
      END IF
      IF (line(1:1) /= "#") EXIT
    END DO

    ! Count rows in the file
    n_rows = 0
    DO
      READ(unit, '(A)', IOSTAT=status) line
      IF (status /= 0) EXIT
      n_rows = n_rows + 1
    END DO

    ! Allocate arrays
    ALLOCATE(wavelengths(n_rows))
    ALLOCATE(flux(n_rows))

    ! Rewind to the first non-comment line
    REWIND(unit)
    DO
      READ(unit, '(A)', IOSTAT=status) line
      IF (status /= 0) THEN
        PRINT *, "Error: Could not rewind file", TRIM(directory)
        STOP
      END IF
      IF (line(1:1) /= "#") EXIT
    END DO

    ! Read and parse data
    i = 0
    DO
      READ(unit, *, IOSTAT=status) temp_wavelength, temp_flux
      IF (status /= 0) EXIT
      i = i + 1
      ! Convert f_lambda to f_nu
      wavelengths(i) = temp_wavelength
      flux(i) = temp_flux
    END DO

    CLOSE(unit)
  END SUBROUTINE LoadSED



!****************************
!Load Filter File
!****************************

  SUBROUTINE LoadFilter(directory, filter_wavelengths, filter_trans)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: directory
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: filter_wavelengths, filter_trans

    CHARACTER(LEN=512) :: line
    INTEGER :: unit, n_rows, status, i
    REAL :: temp_wavelength, temp_trans

    ! Open the file
    unit = 20
    OPEN(unit, FILE=TRIM(directory), STATUS='OLD', ACTION='READ', IOSTAT=status)
    IF (status /= 0) THEN
      PRINT *, "Error: Could not open file ", TRIM(directory)
      STOP
    END IF

    ! Skip header line
    READ(unit, '(A)', IOSTAT=status) line
    IF (status /= 0) THEN
      PRINT *, "Error: Could not read the file", TRIM(directory)
      STOP
    END IF

    ! Count rows in the file
    n_rows = 0
    DO
      READ(unit, '(A)', IOSTAT=status) line
      IF (status /= 0) EXIT
      n_rows = n_rows + 1
    END DO

    ! Allocate arrays
    ALLOCATE(filter_wavelengths(n_rows))
    ALLOCATE(filter_trans(n_rows))

    ! Rewind to the first non-comment line
    REWIND(unit)
    DO
      READ(unit, '(A)', IOSTAT=status) line
      IF (status /= 0) THEN
        PRINT *, "Error: Could not rewind file", TRIM(directory)
        STOP
      END IF
      IF (line(1:1) /= "#") EXIT
    END DO

    ! Read and parse data
    i = 0
    DO
      READ(unit, *, IOSTAT=status) temp_wavelength, temp_trans
      IF (status /= 0) EXIT
      i = i + 1
      
      filter_wavelengths(i) = temp_wavelength
      filter_trans(i) = temp_trans
    END DO

    CLOSE(unit)
  END SUBROUTINE LoadFilter


!****************************
!Load Lookup Table For Identifying Stellar Atmosphere Models
!****************************


  SUBROUTINE LoadLookupTable(lookup_file, lookup_table, out_file_names, out_logg, out_meta, out_teff)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: lookup_file
    REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: lookup_table
    CHARACTER(LEN=100), ALLOCATABLE, INTENT(INOUT) :: out_file_names(:)
    REAL(8), ALLOCATABLE, INTENT(INOUT) :: out_logg(:), out_meta(:), out_teff(:)

    INTEGER :: i, n_rows, status, unit
    CHARACTER(LEN=512) :: line
    CHARACTER(LEN=*), PARAMETER :: delimiter = ","
    CHARACTER(LEN=100), ALLOCATABLE :: columns(:), headers(:)
    INTEGER :: logg_col, meta_col, teff_col

    ! Open the file
    unit = 10
    OPEN(unit, FILE=lookup_file, STATUS='old', ACTION='read', IOSTAT=status)
    IF (status /= 0) THEN
      PRINT *, "Error: Could not open file", lookup_file
      STOP
    END IF

    ! Read header line
    READ(unit, '(A)', IOSTAT=status) line
    IF (status /= 0) THEN
      PRINT *, "Error: Could not read header line"
      STOP
    END IF

    CALL SplitLine(line, delimiter, headers)

    ! Determine column indices for logg, meta, and teff
    logg_col = GetColumnIndex(headers, "logg")
    teff_col = GetColumnIndex(headers, "teff")

    meta_col = GetColumnIndex(headers, "meta")
    IF (meta_col < 0) THEN
      meta_col = GetColumnIndex(headers, "feh")
    END IF 

    n_rows = 0
    DO
      READ(unit, '(A)', IOSTAT=status) line
      IF (status /= 0) EXIT
      n_rows = n_rows + 1
    END DO
    REWIND(unit)

    ! Skip header
    READ(unit, '(A)', IOSTAT=status) line

    ! Allocate output arrays
    ALLOCATE(out_file_names(n_rows))
    ALLOCATE(out_logg(n_rows), out_meta(n_rows), out_teff(n_rows))

    ! Read and parse the file
    i = 0
    DO
      READ(unit, '(A)', IOSTAT=status) line
      IF (status /= 0) EXIT
      i = i + 1

      CALL SplitLine(line, delimiter, columns)

      ! Populate arrays
      out_file_names(i) = columns(1)
      !PRINT *, columns

      IF (logg_col > 0) THEN
        IF (columns(logg_col) /= "") THEN
          READ(columns(logg_col), *) out_logg(i)
        ELSE
          out_logg(i) = 0.0
        END IF
      ELSE
        out_logg(i) = 0.0
      END IF

      IF (meta_col > 0) THEN
        IF (columns(meta_col) /= "") THEN
          READ(columns(meta_col), *) out_meta(i)
        ELSE
          out_meta(i) = 0.0
        END IF
      ELSE
        out_meta(i) = 0.0
      END IF

      IF (teff_col > 0) THEN
        IF (columns(teff_col) /= "") THEN
          READ(columns(teff_col), *) out_teff(i)
        ELSE
          out_teff(i) = 0.0
        END IF
      ELSE
        out_teff(i) = 0.0
      END IF

    END DO

    CLOSE(unit)

  CONTAINS

    FUNCTION GetColumnIndex(headers, target) RESULT(index)
      CHARACTER(LEN=100), INTENT(IN) :: headers(:)
      CHARACTER(LEN=*), INTENT(IN) :: target
      INTEGER :: index, i
      CHARACTER(LEN=100) :: clean_header, clean_target

      index = -1
      clean_target = TRIM(ADJUSTL(target))  ! Clean the target string

      DO i = 1, SIZE(headers)
        clean_header = TRIM(ADJUSTL(headers(i)))  ! Clean each header
        IF (clean_header == clean_target) THEN
          index = i
          EXIT
        END IF
      END DO
    END FUNCTION GetColumnIndex

    SUBROUTINE SplitLine(line, delimiter, tokens)
      CHARACTER(LEN=*), INTENT(IN) :: line, delimiter
      CHARACTER(LEN=100), ALLOCATABLE, INTENT(OUT) :: tokens(:)
      INTEGER :: num_tokens, pos, start, len_delim

      len_delim = LEN_TRIM(delimiter)
      start = 1
      num_tokens = 0
      IF (ALLOCATED(tokens)) DEALLOCATE(tokens)

      DO
        pos = INDEX(line(start:), delimiter)

        IF (pos == 0) EXIT
        num_tokens = num_tokens + 1
        CALL AppendToken(tokens, line(start:start + pos - 2))
        start = start + pos + len_delim - 1
      END DO

      num_tokens = num_tokens + 1
      CALL AppendToken(tokens, line(start:))
    END SUBROUTINE SplitLine

    SUBROUTINE AppendToken(tokens, token)
      CHARACTER(LEN=*), INTENT(IN) :: token
      CHARACTER(LEN=100), ALLOCATABLE, INTENT(INOUT) :: tokens(:)
      CHARACTER(LEN=100), ALLOCATABLE :: temp(:)
      INTEGER :: n

      IF (.NOT. ALLOCATED(tokens)) THEN
        ALLOCATE(tokens(1))
        tokens(1) = token
      ELSE
        n = SIZE(tokens)
        ALLOCATE(temp(n))
        temp = tokens  ! Backup the current tokens
        DEALLOCATE(tokens)  ! Deallocate the old array
        ALLOCATE(tokens(n + 1))  ! Allocate with one extra space
        tokens(1:n) = temp  ! Restore old tokens
        tokens(n + 1) = token  ! Add the new token
        DEALLOCATE(temp)  ! Clean up temporary array
      END IF
    END SUBROUTINE AppendToken

  END SUBROUTINE LoadLookupTable
  





  !###########################################################
  !## MATHS
  !###########################################################

!****************************
!Trapezoidal and Simpson Integration For Flux Calculation
!****************************

  SUBROUTINE TrapezoidalIntegration(x, y, result)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
    REAL(DP), INTENT(OUT) :: result

    INTEGER :: i, n
    REAL :: sum

    n = SIZE(x)
    sum = 0.0

    ! Validate input sizes
    IF (SIZE(x) /= SIZE(y)) THEN
      PRINT *, "Error: x and y arrays must have the same size."
      STOP
    END IF

    IF (SIZE(x) < 2) THEN
      PRINT *, "Error: x and y arrays must have at least 2 elements."
      STOP
    END IF

    ! Perform trapezoidal integration
    DO i = 1, n - 1
      sum = sum + 0.5 * (x(i + 1) - x(i)) * (y(i + 1) + y(i))
    END DO

    result = sum
  END SUBROUTINE TrapezoidalIntegration


SUBROUTINE SimpsonIntegration(x, y, result)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
  REAL(DP), INTENT(OUT) :: result

  INTEGER :: i, n
  REAL(DP) :: sum, h1, h2, f1, f2, f0

  n = SIZE(x)
  sum = 0.0_DP

  ! Validate input sizes
  IF (SIZE(x) /= SIZE(y)) THEN
    PRINT *, "Error: x and y arrays must have the same size."
    STOP
  END IF

  IF (SIZE(x) < 2) THEN
    PRINT *, "Error: x and y arrays must have at least 2 elements."
    STOP
  END IF

  ! Perform adaptive Simpson’s rule
  DO i = 1, n - 2, 2
    h1 = x(i+1) - x(i)       ! Step size for first interval
    h2 = x(i+2) - x(i+1)     ! Step size for second interval

    f0 = y(i)
    f1 = y(i+1)
    f2 = y(i+2)

    ! Simpson's rule: (h/3) * (f0 + 4f1 + f2)
    sum = sum + (h1 + h2) / 6.0_DP * (f0 + 4.0_DP * f1 + f2)
  END DO

  ! Handle the case where n is odd (last interval)
  IF (MOD(n,2) == 0) THEN
    sum = sum + 0.5_DP * (x(n) - x(n-1)) * (y(n) + y(n-1))
  END IF

  result = sum
END SUBROUTINE SimpsonIntegration

SUBROUTINE BooleIntegration(x, y, result)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
  REAL(DP), INTENT(OUT) :: result

  INTEGER :: i, n
  REAL(DP) :: sum, h, f0, f1, f2, f3, f4

  n = SIZE(x)
  sum = 0.0_DP

  ! Validate input sizes
  IF (SIZE(x) /= SIZE(y)) THEN
    PRINT *, "Error: x and y arrays must have the same size."
    STOP
  END IF

  IF (SIZE(x) < 5) THEN
    PRINT *, "Error: x and y arrays must have at least 5 elements."
    STOP
  END IF

  ! Apply Boole's rule
  DO i = 1, n - 4, 4
    h = (x(i+4) - x(i)) / 4.0_DP   ! Step size

    f0 = y(i)
    f1 = y(i+1)
    f2 = y(i+2)
    f3 = y(i+3)
    f4 = y(i+4)

    ! Boole's Rule: (2h/45) * (7f0 + 32f1 + 12f2 + 32f3 + 7f4)
    sum = sum + (2.0_DP * h / 45.0_DP) * (7.0_DP * f0 + 32.0_DP * f1 + 12.0_DP * f2 + 32.0_DP * f3 + 7.0_DP * f4)
  END DO

  ! Handle leftover intervals
  IF (MOD(n, 4) /= 1) THEN
    PRINT *, "Warning: Remaining points not fitting Boole's rule. Applying Simpson's rule for last interval."
    sum = sum + (x(n) - x(n-1)) / 6.0_DP * (y(n-1) + 4.0_DP * y(n-1) + y(n))
  END IF

  result = sum
END SUBROUTINE BooleIntegration

SUBROUTINE RombergIntegration(x, y, result)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
  REAL(DP), INTENT(OUT) :: result

  INTEGER :: i, j, k, n, m
  REAL(DP), DIMENSION(:), ALLOCATABLE :: R
  REAL(DP) :: h, sum, factor

  n = SIZE(x)
  m = INT(LOG(REAL(n, DP)) / LOG(2.0_DP)) + 1  ! Number of refinement levels

  ! Validate input sizes
  IF (SIZE(x) /= SIZE(y)) THEN
    PRINT *, "Error: x and y arrays must have the same size."
    STOP
  END IF

  IF (n < 2) THEN
    PRINT *, "Error: x and y arrays must have at least 2 elements."
    STOP
  END IF

  ALLOCATE(R(m))

  ! Compute initial trapezoidal rule estimate
  h = x(n) - x(1)
  R(1) = 0.5_DP * h * (y(1) + y(n))

  ! Refinement using Romberg's method
  DO j = 2, m
    sum = 0.0_DP
    DO i = 1, 2**(j-2)
      sum = sum + y(1 + (2*i - 1) * (n-1) / (2**(j-1)))
    END DO

    h = h / 2.0_DP
    R(j) = 0.5_DP * R(j-1) + h * sum

    ! Richardson extrapolation
    factor = 4.0_DP
    DO k = j, 2, -1
      R(k-1) = (factor * R(k) - R(k-1)) / (factor - 1.0_DP)
      factor = factor * 4.0_DP
    END DO
  END DO

  result = R(1)
  DEALLOCATE(R)
END SUBROUTINE RombergIntegration

!****************************
!Linear Interpolation For SED Construction
!****************************

  SUBROUTINE LinearInterpolate(x, y, x_val, y_val)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x(:), y(:), x_val
    REAL(DP), INTENT(OUT) :: y_val
    INTEGER :: i
    REAL(DP) :: slope

    ! Validate input sizes
    IF (SIZE(x) < 2) THEN
      PRINT *, "Error: x array has fewer than 2 points."
      y_val = 0.0_DP
      RETURN
    END IF

    IF (SIZE(x) /= SIZE(y)) THEN
      PRINT *, "Error: x and y arrays have different sizes."
      y_val = 0.0_DP
      RETURN
    END IF

    ! Handle out-of-bounds cases
    IF (x_val < MINVAL(x)) THEN
      y_val = y(1)
      RETURN
    ELSE IF (x_val > MAXVAL(x)) THEN
      y_val = y(SIZE(y))
      RETURN
    END IF

    ! Perform interpolation
    DO i = 1, SIZE(x) - 1
      IF (x_val >= x(i) .AND. x_val <= x(i + 1)) THEN
        slope = (y(i + 1) - y(i)) / (x(i + 1) - x(i))
        y_val = y(i) + slope * (x_val - x(i))
        RETURN
      END IF
    END DO

    y_val = 0.0_DP
  END SUBROUTINE LinearInterpolate


!****************************
!Array Interpolation For SED Construction
!****************************

  SUBROUTINE InterpolateArray(x_in, y_in, x_out, y_out)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x_in(:), y_in(:), x_out(:)
    REAL(DP), INTENT(OUT) :: y_out(:)
    INTEGER :: i

    ! Validate input sizes
    IF (SIZE(x_in) < 2 .OR. SIZE(y_in) < 2) THEN
      PRINT *, "Error: x_in or y_in arrays have fewer than 2 points."
      STOP
    END IF

    IF (SIZE(x_in) /= SIZE(y_in)) THEN
      PRINT *, "Error: x_in and y_in arrays have different sizes."
      STOP
    END IF

    IF (SIZE(x_out) <= 0) THEN
      PRINT *, "Error: x_out array is empty."
      STOP
    END IF

    DO i = 1, SIZE(x_out)
      CALL LinearInterpolate(x_in, y_in, x_out(i), y_out(i))
    END DO
  END SUBROUTINE InterpolateArray





!****************************
!SED Interpolation attemps
!****************************


SUBROUTINE ConstructSED_Interpolated(teff, log_g, metallicity, R, d,   &
     &         file_names, lu_teff, lu_logg, lu_meta, stellar_model_dir,  &
     &         wavelengths, fluxes)
  IMPLICIT NONE
  ! Inputs:
  REAL(8), INTENT(IN) :: teff, log_g, metallicity, R, d
  REAL(8), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
  CHARACTER(LEN=*), INTENT(IN) :: stellar_model_dir
  CHARACTER(LEN=100), INTENT(IN) :: file_names(:)
  ! Outputs:
  REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, fluxes

  ! Local variables
  INTEGER :: i0, i1, j0, j1, k0, k1
  INTEGER :: n_points, i
  REAL(DP) :: f_teff, f_logg, f_meta
  REAL(DP), DIMENSION(:), ALLOCATABLE :: common_wavelengths
  REAL(DP), DIMENSION(:), ALLOCATABLE :: interp_surface_flux
  REAL(DP), DIMENSION(:), ALLOCATABLE :: diluted_flux
  ! sed_grid(2,2,2, :) will hold the SED flux arrays (for each wavelength) for the 8 grid corners
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: sed_grid
  REAL(DP), DIMENSION(:), ALLOCATABLE :: temp_wavelengths, temp_flux

  !--------------------------------------------------------------------
  ! Find the bounding indices in each parameter dimension.
  ! These routines should find i0 and i1 such that:
  !    lu_teff(i0) <= teff <= lu_teff(i1)
  ! Similarly for log_g and metallicity.
  ! (You will need to implement these or adapt your existing routines.)
  CALL FindBoundingIndices(teff, lu_teff, i0, i1)
  CALL FindBoundingIndices(log_g, lu_logg, j0, j1)
  CALL FindBoundingIndices(metallicity, lu_meta, k0, k1)
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Load the eight SEDs corresponding to the grid points:
  ! (i0,j0,k0), (i1,j0,k0), (i0,j1,k0), (i1,j1,k0),
  ! (i0,j0,k1), (i1,j0,k1), (i0,j1,k1), (i1,j1,k1)
  ! We assume that file_names is ordered such that the index in lu_teff, etc., matches the SED filename.
  ! (You may need to adapt this if your file naming is more complicated.)
  ! First, load one SED (say, for (i0,j0,k0)) to define the common wavelength grid.
  CALL LoadSED(TRIM(stellar_model_dir)//TRIM(file_names(i0)), i0, temp_wavelengths, temp_flux)
  n_points = SIZE(temp_wavelengths)
  ALLOCATE(common_wavelengths(n_points))
  common_wavelengths = temp_wavelengths
  ! Allocate sed_grid: dimensions 2 x 2 x 2 x n_points
  ALLOCATE(sed_grid(2,2,2, n_points))

  ! Now load each of the eight SEDs and interpolate to the common grid.
  ! For brevity, we assume a helper subroutine "LoadAndInterpolateSED" that loads an SED
  ! from a given file (using LoadSED) and then interpolates it onto common_wavelengths (using InterpolateArray).
  CALL LoadAndInterpolateSED(TRIM(stellar_model_dir)//TRIM(file_names(i0)), i0, common_wavelengths, sed_grid(1,1,1,:))
  CALL LoadAndInterpolateSED(TRIM(stellar_model_dir)//TRIM(file_names(i1)), i1, common_wavelengths, sed_grid(2,1,1,:))
  CALL LoadAndInterpolateSED(TRIM(stellar_model_dir)//TRIM(file_names(j0)), j0, common_wavelengths, sed_grid(1,2,1,:))
  CALL LoadAndInterpolateSED(TRIM(stellar_model_dir)//TRIM(file_names(j1)), j1, common_wavelengths, sed_grid(2,2,1,:))
  CALL LoadAndInterpolateSED(TRIM(stellar_model_dir)//TRIM(file_names(k0)), k0, common_wavelengths, sed_grid(1,1,2,:))
  CALL LoadAndInterpolateSED(TRIM(stellar_model_dir)//TRIM(file_names(k1)), k1, common_wavelengths, sed_grid(2,1,2,:))
  CALL LoadAndInterpolateSED(TRIM(stellar_model_dir)//TRIM(file_names(k0+1)), k0+1, common_wavelengths, sed_grid(1,2,2,:))
  CALL LoadAndInterpolateSED(TRIM(stellar_model_dir)//TRIM(file_names(k1+1)), k1+1, common_wavelengths, sed_grid(2,2,2,:))
  !--------------------------------------------------------------------
  ! Note: The mapping of indices in file_names to (teff, logg, meta) grid coordinates
  ! must be arranged consistently in your system.
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Compute fractional distances along each parameter axis:
  f_teff = (teff - lu_teff(i0)) / (lu_teff(i1) - lu_teff(i0))
  f_logg = (log_g - lu_logg(j0)) / (lu_logg(j1) - lu_logg(j0))
  f_meta = (metallicity - lu_meta(k0)) / (lu_meta(k1) - lu_meta(k0))
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Perform trilinear interpolation at each wavelength point.
  ALLOCATE(interp_surface_flux(n_points))
  DO i = 1, n_points
    interp_surface_flux(i) = &
         (1.0D0 - f_teff)*(1.0D0 - f_logg)*(1.0D0 - f_meta)*sed_grid(1,1,1,i) + &
         f_teff*(1.0D0 - f_logg)*(1.0D0 - f_meta)*sed_grid(2,1,1,i) + &
         (1.0D0 - f_teff)*f_logg*(1.0D0 - f_meta)*sed_grid(1,2,1,i) + &
         f_teff*f_logg*(1.0D0 - f_meta)*sed_grid(2,2,1,i) + &
         (1.0D0 - f_teff)*(1.0D0 - f_logg)*f_meta*sed_grid(1,1,2,i) + &
         f_teff*(1.0D0 - f_logg)*f_meta*sed_grid(2,1,2,i) + &
         (1.0D0 - f_teff)*f_logg*f_meta*sed_grid(1,2,2,i) + &
         f_teff*f_logg*f_meta*sed_grid(2,2,2,i)
  END DO
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Apply the dilution factor to convert the interpolated surface flux density
  ! into an observed flux density at Earth.
  ALLOCATE(diluted_flux(n_points))
  CALL dilute_flux(interp_surface_flux, R, d, diluted_flux)
  !--------------------------------------------------------------------

  ! Set the output arrays.
  ALLOCATE(wavelengths(n_points), fluxes(n_points))
  wavelengths = common_wavelengths
  fluxes = diluted_flux

  ! Deallocate temporary arrays.
  DEALLOCATE(common_wavelengths, interp_surface_flux, diluted_flux)
  ! (Also deallocate sed_grid when done, if appropriate.)
  
END SUBROUTINE ConstructSED_Interpolated


SUBROUTINE FindBoundingIndices(target, grid, i0, i1)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  REAL(8), INTENT(IN) :: target
  REAL(8), INTENT(IN) :: grid(:)
  INTEGER, INTENT(OUT) :: i0, i1
  INTEGER :: i, n

  n = SIZE(grid)

  ! If target is below the grid, return first two indices.
  IF (target <= grid(1)) THEN
    i0 = 1
    i1 = 2
    RETURN
  END IF

  ! If target is above the grid, return the last two indices.
  IF (target >= grid(n)) THEN
    i0 = n - 1
    i1 = n
    RETURN
  END IF

  ! Otherwise, find the indices such that grid(i0) <= target <= grid(i1)
  DO i = 1, n - 1
    IF (grid(i) <= target .AND. target <= grid(i+1)) THEN
      i0 = i
      i1 = i + 1
      RETURN
    END IF
  END DO
END SUBROUTINE FindBoundingIndices



SUBROUTINE LoadAndInterpolateSED(filename, index, common_wavelengths, flux_out)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: index
  REAL(DP), INTENT(IN) :: common_wavelengths(:)
  REAL(DP), INTENT(OUT) :: flux_out(:)
  
  REAL(DP), DIMENSION(:), ALLOCATABLE :: temp_wavelengths, temp_flux

  ! Use your existing LoadSED subroutine to load the SED from the file.
  CALL LoadSED(TRIM(filename), index, temp_wavelengths, temp_flux)
  
  ! Use your existing InterpolateArray to re-grid the loaded SED onto common_wavelengths.
  CALL InterpolateArray(temp_wavelengths, temp_flux, common_wavelengths, flux_out)
  
  ! Clean up temporary arrays.
  DEALLOCATE(temp_wavelengths, temp_flux)
END SUBROUTINE LoadAndInterpolateSED

















end module run_star_extras
