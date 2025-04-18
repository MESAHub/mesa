
MODULE tetra_interp
  USE const_def, ONLY: dp
  USE shared_funcs
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: constructsed_tetra

CONTAINS

  !---------------------------------------------------------------------------
  ! Main entry point: Construct a SED using Hermite tensor interpolation
  !---------------------------------------------------------------------------
  SUBROUTINE constructsed_tetra(teff, log_g, metallicity, R, d, file_names, &
                                  lu_teff, lu_logg, lu_meta, stellar_model_dir, &
                                  wavelengths, fluxes)
    REAL(dp), INTENT(IN) :: teff, log_g, metallicity, R, d
    REAL(dp), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
    CHARACTER(LEN=*), INTENT(IN) :: stellar_model_dir
    CHARACTER(LEN=100), INTENT(IN) :: file_names(:)
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, fluxes

    INTEGER :: i, n_lambda
    REAL(dp), DIMENSION(:), ALLOCATABLE :: interp_flux, diluted_flux
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: flux_cube
    
    ! Parameter grids
    INTEGER :: n_teff, n_logg, n_meta
    REAL(dp), ALLOCATABLE :: teff_grid(:), logg_grid(:), meta_grid(:)

    ! Create parameter grids from available model parameters
    PRINT *, 'Sorting unique Teff...'
    CALL get_unique_sorted(lu_teff, teff_grid)
    PRINT *, 'Unique Teff count:', SIZE(teff_grid)
    
    PRINT *, 'Sorting unique logg...'
    CALL get_unique_sorted(lu_logg, logg_grid)
    PRINT *, 'Unique logg count:', SIZE(logg_grid)
    
    PRINT *, 'Sorting unique metallicity...'
    CALL get_unique_sorted(lu_meta, meta_grid)
    PRINT *, 'Unique metallicity count:', SIZE(meta_grid)

    n_teff = SIZE(teff_grid)
    n_logg = SIZE(logg_grid)
    n_meta = SIZE(meta_grid)
    
    ! Allocate 3D grid for flux values at each wavelength point
    ALLOCATE(flux_cube(n_teff, n_logg, n_meta))

    ! Load first SED to get wavelength grid
    PRINT *, 'Loading first SED for wavelength grid:'
    PRINT *, TRIM(stellar_model_dir) // TRIM(file_names(1))
    CALL loadsed(TRIM(stellar_model_dir) // TRIM(file_names(1)), 1, wavelengths, fluxes)
    PRINT *, 'Wavelength count:', SIZE(wavelengths)
    n_lambda = SIZE(wavelengths)
    
  END SUBROUTINE constructsed_tetra

  !****************************
  ! Find Enclosing Simplex for Tetrahedral Interpolation
  !****************************
  SUBROUTINE findenclosingsimplex(teff, log_g, metallicity, lu_teff, lu_logg, lu_meta, &
                                  simplex_indices, bary_weights)
    REAL(dp), INTENT(IN) :: teff, log_g, metallicity
    REAL(dp), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: simplex_indices(:)
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: bary_weights(:)

    INTEGER :: i, num_points, j, temp_index, k
    REAL(dp), ALLOCATABLE :: dists(:)
    REAL(dp), DIMENSION(3) :: P, P0, P1, P2, P3
    REAL(dp), DIMENSION(4) :: bary
    REAL(dp) :: tol, sumw
    REAL(dp) :: temp_w(4)
    ! Normalization factors
    REAL(dp) :: teff_min, teff_max, logg_min, logg_max, meta_min, meta_max
    REAL(dp) :: t_norm, g_norm, m_norm
    REAL(dp), DIMENSION(3) :: pt, p0t, p1t, p2t, p3t

    ! Set a tolerance appropriate for normalized values
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

    ! Compute distances for each lookup point in normalized space
    DO i = 1, num_points
      pt(1) = (lu_teff(i) - teff_min) / (teff_max - teff_min)
      pt(2) = (lu_logg(i) - logg_min) / (logg_max - logg_min)
      pt(3) = (lu_meta(i) - meta_min) / (meta_max - meta_min)
      dists(i) = SQRT( (pt(1) - t_norm)**2 + (pt(2) - g_norm)**2 + (pt(3) - m_norm)**2 )
    END DO

    ! Find indices of the 4 smallest distances
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

    ! Now form the normalized coordinates for the 4 candidate vertices
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

  !****************************
  ! Calculate Barycentric Coordinates for Tetrahedral Interpolation
  !****************************
  SUBROUTINE ComputeBarycentrics(P, P0, P1, P2, P3, bary)
    REAL(dp), INTENT(IN) :: P(3), P0(3), P1(3), P2(3), P3(3)
    REAL(dp), INTENT(OUT) :: bary(4)
    REAL(dp) :: M(3,3), d, d0, d1, d2, d3
    REAL(dp) :: rhs(3)

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

  !****************************
  ! Determinant of 3x3 Matrix
  !****************************
  FUNCTION det3(M) RESULT(d)
    REAL(dp), INTENT(IN) :: M(3,3)
    REAL(dp) :: d
    d = M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) - &
        M(1,2)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) + &
        M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
  END FUNCTION det3

  end module tetra_interp