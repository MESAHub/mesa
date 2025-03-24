! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
!
! ***********************************************************************

   module colors_lib
      use math_lib
      use colors_def
      use const_def
      ! library for calculating theoretical estimates of magnitudes and colors
      ! from Teff, L, M, and [M/H].

      ! Color-magnitude data shipped with MESA is from:
      ! Lejeune, Cuisinier, Buser (1998) A&AS 130, 65-75.
      ! However, you add your own bolometric corrections files for mesa to use

      ! The data interface for the library is defined in colors_def
      ! Th easiest way to get output is to add the columns to your history_columns.list file

      ! The preferred way for users (in a run_star_extras routine) for accessing the colors data is to
      ! call either get_by_by_name, get_abs_mag_by_name or get_abs_bolometric_mag. Other routines are there
      ! to hook into the rest of MESA.

      ! Routines get_bc will return the coefficients from interpolating over log Teff, log g, [M/H]
      ! even though the tables are defined as Teff, log g, [M/H]. get_abs_mag routines return
      ! data thats been turned into an absolute magnitude. A color can be computed by taking the difference between
      ! two get_bc or two get_abs_mag calls.

      ! Names for the filters should be unique across all data files (left to the user to enforce this).
      ! Name matching is performed in a case sensitive manner.
      ! The names themselves are not important as far as MESA is concerned, you can name each filter (including the
      ! ones MESA ships by defaults) by what ever name you want by editing the data file(s) and changing the names in the header.
      ! MESA does not rely on any particlaur band existing.

      implicit none


      contains  ! the procedure interface for the library
      ! client programs should only call these routines.




   !###########################################################
   !## CUSTOM COLOURS
   !###########################################################

   !****************************
   !Calculate Bolometric Photometry Using Multiple SEDs
   !****************************

     SUBROUTINE calculatebolometric(teff, log_g, metallicity, R, d, bolometric_magnitude, bolometric_flux, wavelengths, fluxes, sed_filepath)
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
       CALL loadlookuptable(lookup_file, lookup_table, file_names, lu_logg, lu_meta, lu_teff)
       !print *, 'logg', lu_logg
       !print *,  'meta', lu_meta
       !print *, 'teff', lu_teff
       ! Interpolate Spectral Energy Distribution
       !CALL constructsed_Robust(teff, log_g, metallicity, R, d, file_names, lu_teff, lu_logg, lu_meta, sed_filepath, wavelengths, fluxes)
       CALL constructsed(teff, log_g, metallicity, R, d, file_names, lu_teff, lu_logg, lu_meta, sed_filepath, wavelengths, fluxes)

       ! Calculate bolometric flux and magnitude
       CALL calculatebolometricphot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
     END SUBROUTINE calculatebolometric



   !****************************
   !Construct SED With Combination of SEDs
   !****************************

   SUBROUTINE constructsed(teff, log_g, metallicity, R, d, file_names, lu_teff, lu_logg, lu_meta, stellar_model_dir, wavelengths, fluxes)
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
     CALL getcloseststellarmodels(teff, log_g, metallicity, lu_teff, lu_logg, lu_meta, closest_indices)

     ! Load the first SED to define the wavelength grid
     CALL loadsed(TRIM(stellar_model_dir) // TRIM(file_names(closest_indices(1))), closest_indices(1), temp_wavelengths, temp_flux)
     n_points = SIZE(temp_wavelengths)
     ALLOCATE(common_wavelengths(n_points))
     common_wavelengths = temp_wavelengths

     ! Allocate flux array for the models (4 models, n_points each)
     ALLOCATE(model_fluxes(4, n_points))
     CALL interpolatearray(temp_wavelengths, temp_flux, common_wavelengths, model_fluxes(1, :))

     ! Load and interpolate remaining SEDs
     DO i = 2, 4
       CALL loadsed(TRIM(stellar_model_dir) // TRIM(file_names(closest_indices(i))), closest_indices(i), temp_wavelengths, temp_flux)
       CALL interpolatearray(temp_wavelengths, temp_flux, common_wavelengths, model_fluxes(i, :))
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

   END SUBROUTINE constructsed


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



   SUBROUTINE dilute_flux(surface_flux, R, d, calibrated_flux)
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

   SUBROUTINE getcloseststellarmodels(teff, log_g, metallicity, lu_teff, lu_logg, lu_meta, closest_indices)
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


       !distance = SQRT(teff_dist**2 + logg_dist**2 + meta_dist**2)
       distance = teff_dist**2 + logg_dist**2 + meta_dist**2   !SQRT is a monotonic transform so pointless?

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
   END SUBROUTINE getcloseststellarmodels




   !****************************
   !Calculate Bolometric Magnitude and Flux
   !****************************

     SUBROUTINE calculatebolometricphot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
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
       CALL rombergintegration(wavelengths, fluxes, bolometric_flux)

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

     bolometric_magnitude = fluxtomagnitude(bolometric_flux)

     END SUBROUTINE calculatebolometricphot












   !###########################################################
   !## Synthetic Photometry
   !###########################################################

   !****************************
   !Calculate Synthetic Photometry Using SED and Filter
   !****************************



   REAL(DP) FUNCTION calculatesynthetic(temperature, gravity, metallicity, ierr, wavelengths, fluxes, filter_wavelengths, filter_trans, filter_filepath, vega_filepath, filter_name, make_sed)

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
       CALL loadfilter(filter_filepath, filter_wavelengths, filter_trans)

       ! Check for invalid gravity input
       IF (gravity <= 0.0_DP) THEN
           ierr = 1
           calculatesynthetic = -1.0_DP
           RETURN
       END IF

       ! Allocate interpolated_filter if not already allocated
       IF (.NOT. ALLOCATED(interpolated_filter)) THEN
           ALLOCATE(interpolated_filter(SIZE(wavelengths)))
           interpolated_filter = 0.0_DP
       END IF

         ! Perform SED convolution
     ALLOCATE(convolved_flux(SIZE(wavelengths)))
     CALL convolvesed(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)

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
       vega_flux = calculatevegaflux(vega_filepath, filter_wavelengths, filter_trans, filter_name, make_sed)

       ! Calculate synthetic flux and magnitude
       CALL calculatesyntheticflux(wavelengths, convolved_flux, synthetic_flux, filter_wavelengths, filter_trans)

       !PRINT *, "VEGA zero point:", vega_flux

       IF (vega_flux > 0.0_DP) THEN
         calculatesynthetic = -2.5 * LOG10(synthetic_flux / vega_flux)
       ELSE
         PRINT *, "Error: Vega flux is zero, magnitude calculation is invalid."
         calculatesynthetic = HUGE(1.0_DP)
       END IF

   END FUNCTION calculatesynthetic




   !****************************
   !Convolve SED With Filter
   !****************************

     SUBROUTINE convolvesed(wavelengths, fluxes, filter_wavelengths, filter_trans, convolved_flux)
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
       CALL interpolatearray(filter_wavelengths, filter_trans, wavelengths, interpolated_filter)

       ! Perform convolution (element-wise multiplication)
       convolved_flux = fluxes * interpolated_filter

       ! Deallocate arrays (optional, depending on context)
       DEALLOCATE(interpolated_filter)
     END SUBROUTINE convolvesed



   !****************************
   !Calculate Synthetic Flux and Magnitude
   !****************************
     SUBROUTINE calculatesyntheticflux(wavelengths, fluxes, synthetic_flux, filter_wavelengths, filter_trans)
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

       CALL rombergintegration(wavelengths, fluxes* wavelengths, integrated_flux)
       CALL rombergintegration(filter_wavelengths, filter_trans * filter_wavelengths, integrated_filter)

       ! Store the total flux
       IF (integrated_filter > 0.0) THEN
           synthetic_flux = integrated_flux / integrated_filter
       ELSE
           PRINT *, "Error: Integrated filter transmission is zero."
           synthetic_flux = -1.0_DP
           RETURN
       END IF

     END SUBROUTINE calculatesyntheticflux



     REAL(DP) FUNCTION fluxtomagnitude(flux)
       REAL(DP), INTENT(IN) :: flux
       !print *, 'flux:', flux
       IF (flux <= 0.0) THEN
         PRINT *, "Error: Flux must be positive to calculate magnitude."
         fluxtomagnitude = 99.0  ! Return an error value
       ELSE
         fluxtomagnitude = -2.5 * LOG10(flux)
       END IF
     END FUNCTION fluxtomagnitude






   FUNCTION calculatevegaflux(vega_filepath, filt_wave, filt_trans, filter_name, make_sed) RESULT(vega_flux)
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
     CALL loadvegased(vega_filepath, vega_wave, vega_flux_arr)

     ! Convolve the Vega SED with the filter transmission.
     CALL convolvesed(vega_wave, vega_flux_arr, filt_wave, filt_trans, conv_flux)

     ! Integrate the convolved Vega SED and the filter transmission.
     CALL rombergintegration(vega_wave, vega_wave*conv_flux, int_flux)
     CALL rombergintegration(filt_wave, filt_wave*filt_trans, int_filter)

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
   END FUNCTION calculatevegaflux










   SUBROUTINE loadvegased(filepath, wavelengths, flux)
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
   END SUBROUTINE loadvegased




   !###########################################################
   !## FILE IO
   !###########################################################


   !****************************
   !Load Filter File
   !****************************

     SUBROUTINE loadfilter(directory, filter_wavelengths, filter_trans)
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
     END SUBROUTINE loadfilter


   !****************************
   !Load Lookup Table For Identifying Stellar Atmosphere Models
   !****************************


     SUBROUTINE loadlookuptable(lookup_file, lookup_table, out_file_names, out_logg, out_meta, out_teff)
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

       CALL splitline(line, delimiter, headers)

       ! Determine column indices for logg, meta, and teff
       logg_col = getcolumnindex(headers, "logg")
       teff_col = getcolumnindex(headers, "teff")

       meta_col = getcolumnindex(headers, "meta")
       IF (meta_col < 0) THEN
         meta_col = getcolumnindex(headers, "feh")
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

         CALL splitline(line, delimiter, columns)

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

       FUNCTION getcolumnindex(headers, target) RESULT(index)
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
       END FUNCTION getcolumnindex

       SUBROUTINE splitline(line, delimiter, tokens)
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
       END SUBROUTINE splitline

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

     END SUBROUTINE loadlookuptable






     !###########################################################
     !## MATHS
     !###########################################################

   !****************************
   !Trapezoidal and Simpson Integration For Flux Calculation
   !****************************

     SUBROUTINE trapezoidalintegration(x, y, result)
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
     END SUBROUTINE trapezoidalintegration


   SUBROUTINE SimpsonIntegration(x, y, result)
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

   SUBROUTINE rombergintegration(x, y, result)
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
   END SUBROUTINE rombergintegration

   !****************************
   !Linear Interpolation For SED Construction
   !****************************

   !will be removed in the future so long as binary search proves consistently faster
     SUBROUTINE linearinterpolate_linearsearch(x, y, x_val, y_val)
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
     END SUBROUTINE linearinterpolate_linearsearch


   SUBROUTINE linearinterpolate(x, y, x_val, y_val)
     REAL(DP), INTENT(IN) :: x(:), y(:), x_val
     REAL(DP), INTENT(OUT) :: y_val
     INTEGER :: low, high, mid

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
     IF (x_val <= x(1)) THEN
       y_val = y(1)
       RETURN
     ELSE IF (x_val >= x(SIZE(x))) THEN
       y_val = y(SIZE(y))
       RETURN
     END IF

     ! Binary search to find the proper interval [x(low), x(low+1)]
     low = 1
     high = SIZE(x)
     DO WHILE (high - low > 1)
       mid = (low + high) / 2
       IF (x(mid) <= x_val) THEN
         low = mid
       ELSE
         high = mid
       END IF
     END DO

     ! Linear interpolation between x(low) and x(low+1)
     y_val = y(low) + (y(low+1) - y(low)) / (x(low+1) - x(low)) * (x_val - x(low))
   END SUBROUTINE linearinterpolate




   !****************************
   !Array Interpolation For SED Construction
   !****************************

     SUBROUTINE interpolatearray(x_in, y_in, x_out, y_out)
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
         CALL linearinterpolate(x_in, y_in, x_out(i), y_out(i))
       END DO
     END SUBROUTINE interpolatearray





   FUNCTION det3(M) RESULT(d)
     IMPLICIT NONE
     REAL(8), INTENT(IN) :: M(3,3)
     REAL(8) :: d
     d = M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) - &
         M(1,2)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) + &
         M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
   END FUNCTION det3

   SUBROUTINE ComputeBarycentrics(P, P0, P1, P2, P3, bary)
     IMPLICIT NONE
     REAL(8), INTENT(IN) :: P(3), P0(3), P1(3), P2(3), P3(3)
     REAL(8), INTENT(OUT) :: bary(4)
     REAL(8) :: M(3,3), d, d0, d1, d2, d3
     REAL(8) :: rhs(3)

     ! Build matrix M with columns = P1-P0, P2-P0, P3-P0
     M(:,1) = P1 - P0
     M(:,2) = P2 - P0
     M(:,3) = P3 - P0

     d = det3(M)
     IF (ABS(d) < 1.0D-12) THEN
       bary = -1.0D0  ! signal degenerate
       RETURN
     END IF

     ! Solve M * [u, v, w]^T = P - P0 using Cramer's rule
     rhs = P - P0
     d0 = det3(reshape([rhs(1), M(1,2), M(1,3), &
                         rhs(2), M(2,2), M(2,3), &
                         rhs(3), M(3,2), M(3,3)], [3,3]))
     d1 = det3(reshape([M(1,1), rhs(1), M(1,3), &
                         M(2,1), rhs(2), M(2,3), &
                         M(3,1), rhs(3), M(3,3)], [3,3]))
     d2 = det3(reshape([M(1,1), M(1,2), rhs(1), &
                         M(2,1), M(2,2), rhs(2), &
                         M(3,1), M(3,2), rhs(3)], [3,3]))
     ! The barycentrics: w0 = 1 - u - v - w, w1 = u, w2 = v, w3 = w.
     bary(2) = d0/d
     bary(3) = d1/d
     bary(4) = d2/d
     bary(1) = 1.0D0 - bary(2) - bary(3) - bary(4)
   END SUBROUTINE ComputeBarycentrics



   SUBROUTINE findenclosingsimplex(teff, log_g, metallicity, lu_teff, lu_logg, lu_meta, &
                                     simplex_indices, bary_weights)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = KIND(1.0D0)
     REAL(8), INTENT(IN) :: teff, log_g, metallicity
     REAL(8), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
     INTEGER, ALLOCATABLE, INTENT(OUT) :: simplex_indices(:)
     REAL(DP), ALLOCATABLE, INTENT(OUT) :: bary_weights(:)
     
     INTEGER :: i, num_points, j, temp_index, k
     REAL(8), ALLOCATABLE :: dists(:)
     REAL(8), DIMENSION(3) :: P, P0, P1, P2, P3
     REAL(8), DIMENSION(4) :: bary
     REAL(8) :: tol, sumw
     REAL(8) :: temp_w(4)
     ! Normalization factors
     REAL(8) :: teff_min, teff_max, logg_min, logg_max, meta_min, meta_max
     REAL(8) :: t_norm, g_norm, m_norm
     REAL(8), DIMENSION(3) :: pt, p0t, p1t, p2t, p3t

     ! Set a tolerance appropriate for normalized values (e.g., 1e-3)
     tol = 1.0D-3
     
     num_points = SIZE(lu_teff)
     ALLOCATE(dists(num_points))
     
     ! Compute min and max from the lookup arrays
     teff_min = MINVAL(lu_teff)
     teff_max = MAXVAL(lu_teff)
     logg_min = MINVAL(lu_logg)
     logg_max = MAXVAL(lu_logg)
     meta_min = MINVAL(lu_meta)
     meta_max = MAXVAL(lu_meta)
     
     ! Normalize the query point
     t_norm = (teff - teff_min) / (teff_max - teff_min)
     g_norm = (log_g - logg_min) / (logg_max - logg_min)
     m_norm = (metallicity - meta_min) / (meta_max - meta_min)
     P = [ t_norm, g_norm, m_norm ]
     
     ! Compute distances for each lookup point in normalized space.
     DO i = 1, num_points
       pt(1) = (lu_teff(i) - teff_min) / (teff_max - teff_min)
       pt(2) = (lu_logg(i) - logg_min) / (logg_max - logg_min)
       pt(3) = (lu_meta(i) - meta_min) / (meta_max - meta_min)
       dists(i) = SQRT( (pt(1) - t_norm)**2 + (pt(2) - g_norm)**2 + (pt(3) - m_norm)**2 )
     END DO

     ! For debugging: you might print the smallest 10 distances here if needed.
     ! [Your debug code would go here.]

     ! Find indices of the 4 smallest distances (simple selection sort for 4 elements)
     ALLOCATE(simplex_indices(4))
     DO i = 1, 4
       simplex_indices(i) = i
     END DO
     DO i = 5, num_points
       IF (dists(i) < dists(simplex_indices(4))) THEN
         simplex_indices(4) = i
         ! Re-sort the 4 indices by distance (simple bubble sort)
         DO j = 1, 3
            IF (dists(simplex_indices(j)) > dists(simplex_indices(j+1))) THEN
               temp_index = simplex_indices(j)
               simplex_indices(j) = simplex_indices(j+1)
               simplex_indices(j+1) = temp_index
            END IF
         END DO
       END IF
     END DO

     ! Now form the normalized coordinates for the 4 candidate vertices:
     p0t(1) = (lu_teff(simplex_indices(1)) - teff_min) / (teff_max - teff_min)
     p0t(2) = (lu_logg(simplex_indices(1)) - logg_min) / (logg_max - logg_min)
     p0t(3) = (lu_meta(simplex_indices(1)) - meta_min) / (meta_max - meta_min)

     p1t(1) = (lu_teff(simplex_indices(2)) - teff_min) / (teff_max - teff_min)
     p1t(2) = (lu_logg(simplex_indices(2)) - logg_min) / (logg_max - logg_min)
     p1t(3) = (lu_meta(simplex_indices(2)) - meta_min) / (meta_max - meta_min)

     p2t(1) = (lu_teff(simplex_indices(3)) - teff_min) / (teff_max - teff_min)
     p2t(2) = (lu_logg(simplex_indices(3)) - logg_min) / (logg_max - logg_min)
     p2t(3) = (lu_meta(simplex_indices(3)) - meta_min) / (meta_max - meta_min)

     p3t(1) = (lu_teff(simplex_indices(4)) - teff_min) / (teff_max - teff_min)
     p3t(2) = (lu_logg(simplex_indices(4)) - logg_min) / (logg_max - logg_min)
     p3t(3) = (lu_meta(simplex_indices(4)) - meta_min) / (meta_max - meta_min)
     
     P0 = p0t
     P1 = p1t
     P2 = p2t
     P3 = p3t

     ! Compute barycentrics in normalized space
     CALL ComputeBarycentrics(P, P0, P1, P2, P3, bary)
     
     ! If any barycentric is less than -tol, consider the tetrahedron degenerate
     IF ( ANY(bary < -tol) ) THEN
       PRINT *, "Warning: Degenerate tetrahedron. Using inverse-distance weighting fallback."
       ALLOCATE(bary_weights(4))
       sumw = 0.0D0
       DO k = 1, 4
         temp_w(k) = 1.0D0 / (dists(simplex_indices(k)) + 1.0D-12)
         sumw = sumw + temp_w(k)
       END DO
       bary_weights = temp_w / sumw
       RETURN
     ELSE
       ALLOCATE(bary_weights(4))
       bary_weights = bary
       RETURN
     END IF

   END SUBROUTINE findenclosingsimplex




   !--------------------------------------------------------------------
   ! A robust SED interpolation using scattered data interpolation.
   ! Ideally, this uses Delaunay triangulation with barycentric interpolation.
   ! The subroutine findenclosingsimplex is the heart of the method.
   !--------------------------------------------------------------------
   SUBROUTINE constructsed_Robust(teff, log_g, metallicity, R, d,  &
            file_names, lu_teff, lu_logg, lu_meta, stellar_model_dir,  &
            wavelengths, fluxes)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = KIND(1.0D0)
     ! Inputs
     REAL(8), INTENT(IN) :: teff, log_g, metallicity, R, d
     REAL(8), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
     CHARACTER(LEN=*), INTENT(IN) :: stellar_model_dir
     CHARACTER(LEN=100), INTENT(IN) :: file_names(:)
     ! Outputs
     REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, fluxes

     ! Local variables
     INTEGER, ALLOCATABLE, DIMENSION(:) :: simplex_indices
     REAL(DP), ALLOCATABLE, DIMENSION(:) :: bary_weights
     INTEGER :: num_vertices, n_points, i
     REAL(DP), DIMENSION(:), ALLOCATABLE :: common_wavelengths
     REAL(DP), DIMENSION(:), ALLOCATABLE :: interp_flux, temp_wavelength, temp_flux

     !--------------------------------------------------------------------
     ! Step 1: Find the simplex that encloses (teff, log_g, metallicity)
     CALL findenclosingsimplex(teff, log_g, metallicity, lu_teff, lu_logg, lu_meta, &
                               simplex_indices, bary_weights)
     num_vertices = SIZE(simplex_indices)
     PRINT *, "Enclosing simplex indices: ", simplex_indices
     PRINT *, "Barycentric weights: ", bary_weights

     !--------------------------------------------------------------------
     ! Step 2: Define a common wavelength grid.
     ! If we have only one vertex (nearest neighbor fallback), just use that SED.
     CALL loadsed(TRIM(stellar_model_dir)//TRIM(file_names(simplex_indices(1))), &
                  simplex_indices(1), temp_wavelength, temp_flux)
     n_points = SIZE(temp_wavelength)
     IF (n_points <= 0) THEN
       PRINT *, "Error: Loaded SED from ", TRIM(file_names(simplex_indices(1))), " has no wavelengths."
       STOP
     END IF
     ALLOCATE(common_wavelengths(n_points))
     common_wavelengths = temp_wavelength
     DEALLOCATE(temp_wavelength, temp_flux)

     !--------------------------------------------------------------------
     ! Step 3: Compute the SED.
     ALLOCATE(interp_flux(n_points))
     interp_flux = 0.0D0
     IF (num_vertices == 1) THEN
       ! Nearest neighbor: just load the SED without interpolation.
       CALL loadandinterpolatesed(TRIM(stellar_model_dir)//TRIM(file_names(simplex_indices(1))), simplex_indices(1), common_wavelengths, interp_flux)
     ELSE
       DO i = 1, num_vertices
         !print *, TRIM(stellar_model_dir)//TRIM(file_names(simplex_indices(i))), simplex_indices(i), common_wavelengths, temp_flux
         CALL loadandinterpolatesed(TRIM(stellar_model_dir)//TRIM(file_names(simplex_indices(i))), simplex_indices(i), common_wavelengths, temp_flux)
         ! Check that temp_flux has the expected size:
         IF (SIZE(temp_flux) /= n_points) THEN
            PRINT *, "Error: SED from ", TRIM(file_names(simplex_indices(i))), " has mismatched wavelength grid."
            STOP
         END IF
         interp_flux = interp_flux + bary_weights(i) * temp_flux
         DEALLOCATE(temp_flux)
       END DO
     END IF

     !--------------------------------------------------------------------
     ! Step 4: Apply the dilution factor (R/d)^2.
     ALLOCATE(fluxes(n_points))
     CALL dilute_flux(interp_flux, R, d, fluxes)
     ALLOCATE(wavelengths(n_points))
     wavelengths = common_wavelengths

     ! Clean up
     DEALLOCATE(common_wavelengths, interp_flux)
     
   END SUBROUTINE constructsed_Robust




   SUBROUTINE loadandinterpolatesed(filename, index, common_wavelengths, flux_out)
     INTEGER, PARAMETER :: DP = KIND(1.0D0)
     CHARACTER(LEN=*), INTENT(IN) :: filename
     INTEGER, INTENT(IN) :: index
     REAL(DP), INTENT(IN) :: common_wavelengths(:)
     REAL(DP), ALLOCATABLE, INTENT(OUT) :: flux_out(:)

     REAL(DP), DIMENSION(:), ALLOCATABLE :: temp_wavelengths, temp_flux
     INTEGER :: n, i, n_print

     ! Load the SED from the file.
     CALL loadsed(TRIM(filename), index, temp_wavelengths, temp_flux)
     
     ! Check that the loaded data arrays have at least 2 points.
     IF (SIZE(temp_wavelengths) < 2 .OR. SIZE(temp_flux) < 2) THEN
       PRINT *, "Error: Loaded SED arrays are too small."
       STOP
     END IF

     ! Allocate flux_out to match the size of the common wavelength grid.
     n = SIZE(common_wavelengths)
     ALLOCATE(flux_out(n))
     
     ! Interpolate the loaded SED onto the common wavelength grid.
     CALL interpolatearray(temp_wavelengths, temp_flux, common_wavelengths, flux_out)

     ! Clean up temporary arrays.
     DEALLOCATE(temp_wavelengths, temp_flux)
   END SUBROUTINE loadandinterpolatesed

   !****************************
   !Load SED File
   !****************************

     SUBROUTINE loadsed(directory, index, wavelengths, flux)
       CHARACTER(LEN=*), INTENT(IN) :: directory
       INTEGER, INTENT(IN) :: index
       REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, flux

       CHARACTER(LEN=512) :: line
       INTEGER :: unit, n_rows, status, i
       REAL(DP) :: temp_wavelength, temp_flux

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

     END SUBROUTINE loadsed




!things needed to make it compile


subroutine colors_init(num_files,fnames,num_colors,ierr)
   integer, intent(in) :: num_files
   integer, dimension(:), intent(in) :: num_colors
   character(len=*), dimension(:), intent(in) :: fnames
   integer, intent(out) :: ierr
   ierr=0
end subroutine colors_init



subroutine colors_shutdown ()
 
end subroutine colors_shutdown



real(dp) function get_bc_by_name(name,log_Teff,log_g, M_div_h, ierr)
! input
character(len=*),intent(in) :: name
real(dp), intent(in) :: log_Teff ! log10 of surface temp
real(dp), intent(in) :: log_g ! log_10 of surface gravity
real(dp), intent(in) :: M_div_h ! [M/H]
integer, intent(inout) :: ierr
integer :: i,j,n_colors

get_bc_by_name=-99.9d0
ierr=0

end function get_bc_by_name





real(dp) function get_bc_by_id(id,log_Teff,log_g, M_div_h, ierr)
! input
integer, intent(in) :: id
real(dp), intent(in) :: log_Teff ! log10 of surface temp
real(dp), intent(in) :: log_g ! log_10 of surface gravity
real(dp), intent(in) :: M_div_h ! [M/H]
integer, intent(inout) :: ierr
character(len=strlen) :: name

get_bc_by_id=-99.9d0
ierr=0

end function get_bc_by_id




integer function get_bc_id_by_name(name,ierr)
! input
character(len=*), intent(in) :: name
integer, intent(inout) :: ierr
integer :: i,j,k

get_bc_id_by_name=-1
ierr=0

end function get_bc_id_by_name




character(len=strlen) function get_bc_name_by_id(id,ierr)
! input
integer, intent(in) :: id
integer, intent(inout) :: ierr
integer :: i,j,k

get_bc_name_by_id=''
ierr=0

end function get_bc_name_by_id



real(dp) function get_abs_bolometric_mag(lum)
use const_def
real(dp), intent(in) :: lum ! Luminsoity in lsun units

get_abs_bolometric_mag = -99.9d0

end function get_abs_bolometric_mag





real(dp) function get_abs_mag_by_name(name,log_Teff,log_g, M_div_h,lum, ierr)
! input
character(len=*) :: name
real(dp), intent(in) :: log_Teff ! log10 of surface temp
real(dp), intent(in) :: M_div_h ! [M/H]
real(dp), intent(in) :: log_g ! log_10 of surface gravity
real(dp), intent(in) :: lum ! Luminsoity in lsun units
integer, intent(inout) :: ierr

ierr=0
get_abs_mag_by_name=-99.9d0

end function get_abs_mag_by_name





real(dp) function get_abs_mag_by_id(id,log_Teff,log_g, M_div_h,lum, ierr)
! input
integer, intent(in) :: id
real(dp), intent(in) :: log_Teff ! log10 of surface temp
real(dp), intent(in) :: log_g ! log_10 of surface gravity
real(dp), intent(in) :: M_div_h ! [M/H]
real(dp), intent(in) :: lum ! Luminsoity in lsun units
integer, intent(inout) :: ierr
character(len=strlen) :: name

ierr=0
get_abs_mag_by_id=-99.9d0

end function get_abs_mag_by_id




subroutine get_bcs_all(log_Teff, log_g, M_div_h, results, ierr)
! input
real(dp), intent(in) :: log_Teff ! log10 of surface temp
real(dp), intent(in) :: M_div_h ! [M/H]
! output
real(dp),dimension(:), intent(out) :: results
real(dp), intent(in) :: log_g
integer, intent(inout) :: ierr
integer :: i,iStart,iEnd

ierr=0
results(:)=-99.d0

end subroutine get_bcs_all



!Returns in lsun units
real(dp) function get_lum_band_by_name(name,log_Teff,log_g, M_div_h, lum, ierr)
! input
character(len=*) :: name
real(dp), intent(in) :: log_Teff ! log10 of surface temp
real(dp), intent(in) :: M_div_h ! [M/H]
real(dp), intent(in) :: log_g ! log_10 of surface gravity
real(dp), intent(in) :: lum ! Total luminsoity in lsun units
real(dp) :: solar_abs_mag, star_abs_mag
integer, intent(inout) :: ierr

ierr=0
get_lum_band_by_name=-99.d0

end function get_lum_band_by_name

!Returns in lsun units
real(dp) function get_lum_band_by_id(id,log_Teff,log_g, M_div_h, lum, ierr)
! input
integer, intent(in) :: id
real(dp), intent(in) :: log_Teff ! log10 of surface temp
real(dp), intent(in) :: log_g ! log_10 of surface gravity
real(dp), intent(in) :: M_div_h ! [M/H]
real(dp), intent(in) :: lum ! Total luminsoity in lsun units
real(dp) :: solar_abs_mag, star_abs_mag
integer, intent(inout) :: ierr

ierr=0
get_lum_band_by_id=-99.d0

end function get_lum_band_by_id

   end module colors_lib

