! ***********************************************************************
! K-Nearest Neighbors interpolation module for spectral energy distributions (SEDs)
! ***********************************************************************

MODULE knn_interp
  USE const_def, ONLY: dp
  USE shared_funcs, only: dilute_flux, loadsed
  implicit none

  PRIVATE
  PUBLIC :: constructsed_knn, loadsed, interpolatearray, dilute_flux

CONTAINS

  !---------------------------------------------------------------------------
  ! Main entry point: Construct a SED using KNN interpolation
  !---------------------------------------------------------------------------
  SUBROUTINE constructsed_knn(teff, log_g, metallicity, R, d, file_names, &
                         lu_teff, lu_logg, lu_meta, stellar_model_dir, &
                         wavelengths, fluxes)
    REAL(dp), INTENT(IN) :: teff, log_g, metallicity, R, d
    REAL(dp), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
    CHARACTER(LEN=*), INTENT(IN) :: stellar_model_dir
    CHARACTER(LEN=100), INTENT(IN) :: file_names(:)
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, fluxes

    INTEGER, DIMENSION(4) :: closest_indices
    REAL(dp), DIMENSION(:), ALLOCATABLE :: temp_wavelengths, temp_flux, common_wavelengths
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: model_fluxes
    REAL(dp), DIMENSION(4) :: weights, distances
    INTEGER :: i, n_points
    REAL(dp) :: sum_weights
    REAL(dp), DIMENSION(:), ALLOCATABLE :: diluted_flux

    ! Get the four closest stellar models
    CALL getcloseststellarmodels(teff, log_g, metallicity, lu_teff, &
                                lu_logg, lu_meta, closest_indices)

    ! Load the first SED to define the wavelength grid
    CALL loadsed(TRIM(stellar_model_dir) // TRIM(file_names(closest_indices(1))), &
                closest_indices(1), temp_wavelengths, temp_flux)

    n_points = SIZE(temp_wavelengths)
    ALLOCATE(common_wavelengths(n_points))
    common_wavelengths = temp_wavelengths

    ! Allocate flux array for the models (4 models, n_points each)
    ALLOCATE(model_fluxes(4, n_points))
    CALL interpolatearray(temp_wavelengths, temp_flux, common_wavelengths, model_fluxes(1, :))

    ! Load and interpolate remaining SEDs
    DO i = 2, 4
      CALL loadsed(TRIM(stellar_model_dir) // TRIM(file_names(closest_indices(i))), &
                  closest_indices(i), temp_wavelengths, temp_flux)

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

  END SUBROUTINE constructsed_knn

  !---------------------------------------------------------------------------
  ! Identify the four closest stellar models
  !---------------------------------------------------------------------------
  SUBROUTINE getcloseststellarmodels(teff, log_g, metallicity, lu_teff, &
                                  lu_logg, lu_meta, closest_indices)
    REAL(dp), INTENT(IN) :: teff, log_g, metallicity
    REAL(dp), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
    INTEGER, DIMENSION(4), INTENT(OUT) :: closest_indices

    INTEGER :: i, n, j
    REAL(dp) :: distance, norm_teff, norm_logg, norm_meta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: scaled_lu_teff, scaled_lu_logg, scaled_lu_meta
    REAL(dp), DIMENSION(4) :: min_distances
    INTEGER, DIMENSION(4) :: indices
    REAL(dp) :: teff_min, teff_max, logg_min, logg_max, meta_min, meta_max
    REAL(dp) :: teff_dist, logg_dist, meta_dist

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

      ! Using squared distance without sqrt (monotonic transform)
      distance = teff_dist**2 + logg_dist**2 + meta_dist**2

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
  END SUBROUTINE getcloseststellarmodels

  !---------------------------------------------------------------------------
  ! Linear interpolation (binary search version for efficiency)
  !---------------------------------------------------------------------------
  SUBROUTINE linearinterpolate(x, y, x_val, y_val)
    REAL(dp), INTENT(IN) :: x(:), y(:), x_val
    REAL(dp), INTENT(OUT) :: y_val
    INTEGER :: low, high, mid

    ! Validate input sizes
    IF (SIZE(x) < 2) THEN
      PRINT *, "Error: x array has fewer than 2 points."
      y_val = 0.0_dp
      RETURN
    END IF

    IF (SIZE(x) /= SIZE(y)) THEN
      PRINT *, "Error: x and y arrays have different sizes."
      y_val = 0.0_dp
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

  !---------------------------------------------------------------------------
  ! Array interpolation for SED construction
  !---------------------------------------------------------------------------
  SUBROUTINE interpolatearray(x_in, y_in, x_out, y_out)
    REAL(dp), INTENT(IN) :: x_in(:), y_in(:), x_out(:)
    REAL(dp), INTENT(OUT) :: y_out(:)
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

END MODULE knn_interp