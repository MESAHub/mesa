! ***********************************************************************
!
! Hermite interpolation module for spectral energy distributions (SEDs)
! ***********************************************************************

MODULE hermite_interp
  USE const_def, ONLY: dp
  USE shared_func
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: constructsed_hermite, hermite_tensor_interp3d

CONTAINS

  !---------------------------------------------------------------------------
  ! Main entry point: Construct a SED using Hermite tensor interpolation
  !---------------------------------------------------------------------------
  SUBROUTINE constructsed_hermite(teff, log_g, metallicity, R, d, file_names, &
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
    
    ! Allocate space for interpolated flux
    ALLOCATE(interp_flux(n_lambda))

    ! Process each wavelength point
    PRINT *, 'Beginning interpolation over all wavelengths...'
    DO i = 1, n_lambda
      IF (MOD(i, 100) == 0) PRINT *, '  Interpolating wavelength index:', i
      
      ! Build the 3D grid for this wavelength
      CALL build_flux_grid_at_lambda(i, file_names, stellar_model_dir, &
                                    teff_grid, logg_grid, meta_grid, flux_cube)
      
      ! Interpolate at the target parameters
      interp_flux(i) = hermite_tensor_interp3d(teff, log_g, metallicity, &
                                              teff_grid, logg_grid, meta_grid, flux_cube)
    END DO

    ! Apply distance dilution to get observed flux
    ALLOCATE(diluted_flux(n_lambda))
    CALL dilute_flux(interp_flux, R, d, diluted_flux)
    fluxes = diluted_flux
    
    ! Clean up
    DEALLOCATE(diluted_flux, interp_flux)
  END SUBROUTINE constructsed_hermite

  !---------------------------------------------------------------------------
  ! Build the 3D flux grid for a single wavelength point
  !---------------------------------------------------------------------------
  SUBROUTINE build_flux_grid_at_lambda(lambda_index, file_names, stellar_model_dir, &
                                      teff_grid, logg_grid, meta_grid, flux_cube)
    INTEGER, INTENT(IN) :: lambda_index
    CHARACTER(LEN=100), INTENT(IN) :: file_names(:)
    CHARACTER(LEN=*), INTENT(IN) :: stellar_model_dir
    REAL(dp), INTENT(IN) :: teff_grid(:), logg_grid(:), meta_grid(:)
    REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: flux_cube

    INTEGER :: i, j, k, model_index
    CHARACTER(LEN=256) :: model_path
    REAL(dp), DIMENSION(:), ALLOCATABLE :: temp_wavelengths, temp_flux

    model_index = 1

    DO i = 1, SIZE(teff_grid)
      DO j = 1, SIZE(logg_grid)
        DO k = 1, SIZE(meta_grid)
          IF (MOD(model_index, 500) == 0) THEN
            PRINT *, '  Loading model ', model_index, ' for lambda=', lambda_index
          END IF
          
          IF (model_index <= SIZE(file_names)) THEN
            model_path = TRIM(stellar_model_dir) // TRIM(file_names(model_index))
            CALL loadsed(model_path, model_index, temp_wavelengths, temp_flux)
            
            IF (ALLOCATED(temp_flux) .AND. lambda_index <= SIZE(temp_flux)) THEN
              flux_cube(i,j,k) = temp_flux(lambda_index)
            ELSE
              flux_cube(i,j,k) = 0.0_dp
            END IF
          ELSE
            flux_cube(i,j,k) = 0.0_dp
          END IF
          
          model_index = model_index + 1
        END DO
      END DO
    END DO

    IF (ALLOCATED(temp_wavelengths)) DEALLOCATE(temp_wavelengths)
    IF (ALLOCATED(temp_flux)) DEALLOCATE(temp_flux)
    
    PRINT *, '  Finished lambda ', lambda_index
  END SUBROUTINE build_flux_grid_at_lambda

  !---------------------------------------------------------------------------
  ! Get sorted unique values from an array
  !---------------------------------------------------------------------------
  SUBROUTINE get_unique_sorted(input_array, output_array)
    REAL(dp), INTENT(IN) :: input_array(:)
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: output_array(:)

    INTEGER :: i, n, unique_count
    REAL(dp), ALLOCATABLE :: temp(:)
    REAL(dp), PARAMETER :: TOLERANCE = 1.0D-8

    n = SIZE(input_array)
    ALLOCATE(temp(n))
    temp = input_array

    ! Sort the array
    CALL sort_array(temp)

    ! Count unique values (with tolerance)
    unique_count = 1
    DO i = 2, n
      IF (ABS(temp(i) - temp(unique_count)) > TOLERANCE) THEN
        unique_count = unique_count + 1
        temp(unique_count) = temp(i)
      END IF
    END DO

    ! Create output array with unique values
    ALLOCATE(output_array(unique_count))
    output_array = temp(1:unique_count)
    DEALLOCATE(temp)
  END SUBROUTINE get_unique_sorted

  !---------------------------------------------------------------------------
  ! Sort an array using selection sort (simple but effective for small arrays)
  !---------------------------------------------------------------------------
  SUBROUTINE sort_array(array)
    REAL(dp), INTENT(INOUT) :: array(:)
    INTEGER :: i, j, min_idx
    REAL(dp) :: temp

    DO i = 1, SIZE(array)-1
      min_idx = i
      DO j = i+1, SIZE(array)
        IF (array(j) < array(min_idx)) min_idx = j
      END DO
      IF (min_idx /= i) THEN
        temp = array(i)
        array(i) = array(min_idx)
        array(min_idx) = temp
      END IF
    END DO
  END SUBROUTINE sort_array

  !---------------------------------------------------------------------------
  ! Main 3D Hermite interpolation function
  !---------------------------------------------------------------------------
  FUNCTION hermite_tensor_interp3d(x_val, y_val, z_val, x_grid, y_grid, z_grid, f_values) RESULT(f_interp)
    REAL(dp), INTENT(IN) :: x_val, y_val, z_val
    REAL(dp), INTENT(IN) :: x_grid(:), y_grid(:), z_grid(:)
    REAL(dp), INTENT(IN) :: f_values(:,:,:)
    REAL(dp) :: f_interp
    
    INTEGER :: i_x, i_y, i_z
    REAL(dp) :: t_x, t_y, t_z
    REAL(dp) :: dx, dy, dz
    REAL(dp) :: dx_values(2,2,2), dy_values(2,2,2), dz_values(2,2,2)
    REAL(dp) :: h00_x, h10_x, h01_x, h11_x
    REAL(dp) :: h00_y, h10_y, h01_y, h11_y
    REAL(dp) :: h00_z, h10_z, h01_z, h11_z
    REAL(dp) :: values(2,2,2)
    REAL(dp) :: sum
    INTEGER :: ix, iy, iz
    
    ! Find containing cell and parameter values
    CALL find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                              i_x, i_y, i_z, t_x, t_y, t_z)
    
    ! If outside grid, use nearest point
    IF (i_x < 1 .OR. i_x >= SIZE(x_grid) .OR. &
        i_y < 1 .OR. i_y >= SIZE(y_grid) .OR. &
        i_z < 1 .OR. i_z >= SIZE(z_grid)) THEN
      
      ! Find nearest grid point
      CALL find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                             i_x, i_y, i_z)
      
      ! Use the value at the nearest point
      f_interp = f_values(i_x, i_y, i_z)
      RETURN
    END IF
    
    ! Grid cell spacing
    dx = x_grid(i_x+1) - x_grid(i_x)
    dy = y_grid(i_y+1) - y_grid(i_y)
    dz = z_grid(i_z+1) - z_grid(i_z)
    
    ! Extract the local 2x2x2 grid cell
    DO iz = 0, 1
      DO iy = 0, 1
        DO ix = 0, 1
          values(ix+1, iy+1, iz+1) = f_values(i_x+ix, i_y+iy, i_z+iz)
          
          ! Compute derivatives at each corner using finite differences
          CALL compute_derivatives_at_point(f_values, i_x+ix, i_y+iy, i_z+iz, &
                                           SIZE(x_grid), SIZE(y_grid), SIZE(z_grid), &
                                           dx, dy, dz, &
                                           dx_values(ix+1,iy+1,iz+1), &
                                           dy_values(ix+1,iy+1,iz+1), &
                                           dz_values(ix+1,iy+1,iz+1))
        END DO
      END DO
    END DO
    
    ! Evaluate Hermite basis functions
    h00_x = h00(t_x)
    h10_x = h10(t_x)
    h01_x = h01(t_x)
    h11_x = h11(t_x)
    
    h00_y = h00(t_y)
    h10_y = h10(t_y)
    h01_y = h01(t_y)
    h11_y = h11(t_y)
    
    h00_z = h00(t_z)
    h10_z = h10(t_z)
    h01_z = h01(t_z)
    h11_z = h11(t_z)
    
    ! Perform tensor product interpolation
    sum = 0.0_dp
    
    ! Function values at corners
    DO iz = 1, 2
      DO iy = 1, 2
        DO ix = 1, 2
          ! Function values
          IF (ix == 1) THEN
            IF (iy == 1) THEN
              IF (iz == 1) THEN
                sum = sum + h00_x * h00_y * h00_z * values(ix,iy,iz)
              ELSE
                sum = sum + h00_x * h00_y * h01_z * values(ix,iy,iz)
              END IF
            ELSE
              IF (iz == 1) THEN
                sum = sum + h00_x * h01_y * h00_z * values(ix,iy,iz)
              ELSE
                sum = sum + h00_x * h01_y * h01_z * values(ix,iy,iz)
              END IF
            END IF
          ELSE
            IF (iy == 1) THEN
              IF (iz == 1) THEN
                sum = sum + h01_x * h00_y * h00_z * values(ix,iy,iz)
              ELSE
                sum = sum + h01_x * h00_y * h01_z * values(ix,iy,iz)
              END IF
            ELSE
              IF (iz == 1) THEN
                sum = sum + h01_x * h01_y * h00_z * values(ix,iy,iz)
              ELSE
                sum = sum + h01_x * h01_y * h01_z * values(ix,iy,iz)
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    
    ! x-derivatives
    DO iz = 1, 2
      DO iy = 1, 2
        DO ix = 1, 2
          IF (ix == 1) THEN
            IF (iy == 1) THEN
              IF (iz == 1) THEN
                sum = sum + h10_x * h00_y * h00_z * dx * dx_values(ix,iy,iz)
              ELSE
                sum = sum + h10_x * h00_y * h01_z * dx * dx_values(ix,iy,iz)
              END IF
            ELSE
              IF (iz == 1) THEN
                sum = sum + h10_x * h01_y * h00_z * dx * dx_values(ix,iy,iz)
              ELSE
                sum = sum + h10_x * h01_y * h01_z * dx * dx_values(ix,iy,iz)
              END IF
            END IF
          ELSE
            IF (iy == 1) THEN
              IF (iz == 1) THEN
                sum = sum + h11_x * h00_y * h00_z * dx * dx_values(ix,iy,iz)
              ELSE
                sum = sum + h11_x * h00_y * h01_z * dx * dx_values(ix,iy,iz)
              END IF
            ELSE
              IF (iz == 1) THEN
                sum = sum + h11_x * h01_y * h00_z * dx * dx_values(ix,iy,iz)
              ELSE
                sum = sum + h11_x * h01_y * h01_z * dx * dx_values(ix,iy,iz)
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    
    ! y-derivatives
    DO iz = 1, 2
      DO iy = 1, 2
        DO ix = 1, 2
          IF (ix == 1) THEN
            IF (iy == 1) THEN
              IF (iz == 1) THEN
                sum = sum + h00_x * h10_y * h00_z * dy * dy_values(ix,iy,iz)
              ELSE
                sum = sum + h00_x * h10_y * h01_z * dy * dy_values(ix,iy,iz)
              END IF
            ELSE
              IF (iz == 1) THEN
                sum = sum + h00_x * h11_y * h00_z * dy * dy_values(ix,iy,iz)
              ELSE
                sum = sum + h00_x * h11_y * h01_z * dy * dy_values(ix,iy,iz)
              END IF
            END IF
          ELSE
            IF (iy == 1) THEN
              IF (iz == 1) THEN
                sum = sum + h01_x * h10_y * h00_z * dy * dy_values(ix,iy,iz)
              ELSE
                sum = sum + h01_x * h10_y * h01_z * dy * dy_values(ix,iy,iz)
              END IF
            ELSE
              IF (iz == 1) THEN
                sum = sum + h01_x * h11_y * h00_z * dy * dy_values(ix,iy,iz)
              ELSE
                sum = sum + h01_x * h11_y * h01_z * dy * dy_values(ix,iy,iz)
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    
    ! z-derivatives
    DO iz = 1, 2
      DO iy = 1, 2
        DO ix = 1, 2
          IF (ix == 1) THEN
            IF (iy == 1) THEN
              IF (iz == 1) THEN
                sum = sum + h00_x * h00_y * h10_z * dz * dz_values(ix,iy,iz)
              ELSE
                sum = sum + h00_x * h00_y * h11_z * dz * dz_values(ix,iy,iz)
              END IF
            ELSE
              IF (iz == 1) THEN
                sum = sum + h00_x * h01_y * h10_z * dz * dz_values(ix,iy,iz)
              ELSE
                sum = sum + h00_x * h01_y * h11_z * dz * dz_values(ix,iy,iz)
              END IF
            END IF
          ELSE
            IF (iy == 1) THEN
              IF (iz == 1) THEN
                sum = sum + h01_x * h00_y * h10_z * dz * dz_values(ix,iy,iz)
              ELSE
                sum = sum + h01_x * h00_y * h11_z * dz * dz_values(ix,iy,iz)
              END IF
            ELSE
              IF (iz == 1) THEN
                sum = sum + h01_x * h01_y * h10_z * dz * dz_values(ix,iy,iz)
              ELSE
                sum = sum + h01_x * h01_y * h11_z * dz * dz_values(ix,iy,iz)
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    
    f_interp = sum
  END FUNCTION hermite_tensor_interp3d

  !---------------------------------------------------------------------------
  ! Find the cell containing the interpolation point 
  !---------------------------------------------------------------------------
  SUBROUTINE find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                  i_x, i_y, i_z, t_x, t_y, t_z)
    REAL(dp), INTENT(IN) :: x_val, y_val, z_val
    REAL(dp), INTENT(IN) :: x_grid(:), y_grid(:), z_grid(:)
    INTEGER, INTENT(OUT) :: i_x, i_y, i_z
    REAL(dp), INTENT(OUT) :: t_x, t_y, t_z
    
    ! Find x interval
    CALL find_interval(x_grid, x_val, i_x, t_x)
    
    ! Find y interval
    CALL find_interval(y_grid, y_val, i_y, t_y)
    
    ! Find z interval
    CALL find_interval(z_grid, z_val, i_z, t_z)
  END SUBROUTINE find_containing_cell

  !---------------------------------------------------------------------------
  ! Find the interval in a sorted array containing a value
  !---------------------------------------------------------------------------
  SUBROUTINE find_interval(x, val, i, t)
    REAL(dp), INTENT(IN) :: x(:), val
    INTEGER, INTENT(OUT) :: i
    REAL(dp), INTENT(OUT) :: t
    
    INTEGER :: n, lo, hi, mid
    
    n = SIZE(x)
    
    ! Handle out-of-bounds cases
    IF (val <= x(1)) THEN
      i = 1
      t = 0.0_dp
      RETURN
    ELSE IF (val >= x(n)) THEN
      i = n-1
      t = 1.0_dp
      RETURN
    END IF
    
    ! Binary search to find interval
    lo = 1
    hi = n
    DO WHILE (hi - lo > 1)
      mid = (lo + hi) / 2
      IF (val >= x(mid)) THEN
        lo = mid
      ELSE
        hi = mid
      END IF
    END DO
    
    i = lo
    t = (val - x(i)) / (x(i+1) - x(i))
  END SUBROUTINE find_interval

  !---------------------------------------------------------------------------
  ! Find the nearest grid point
  !---------------------------------------------------------------------------
  SUBROUTINE find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z)
    REAL(dp), INTENT(IN) :: x_val, y_val, z_val
    REAL(dp), INTENT(IN) :: x_grid(:), y_grid(:), z_grid(:)
    INTEGER, INTENT(OUT) :: i_x, i_y, i_z
    
    INTEGER :: i
    REAL(dp) :: min_dist, dist
    
    ! Find nearest x grid point
    min_dist = ABS(x_val - x_grid(1))
    i_x = 1
    DO i = 2, SIZE(x_grid)
      dist = ABS(x_val - x_grid(i))
      IF (dist < min_dist) THEN
        min_dist = dist
        i_x = i
      END IF
    END DO
    
    ! Find nearest y grid point
    min_dist = ABS(y_val - y_grid(1))
    i_y = 1
    DO i = 2, SIZE(y_grid)
      dist = ABS(y_val - y_grid(i))
      IF (dist < min_dist) THEN
        min_dist = dist
        i_y = i
      END IF
    END DO
    
    ! Find nearest z grid point
    min_dist = ABS(z_val - z_grid(1))
    i_z = 1
    DO i = 2, SIZE(z_grid)
      dist = ABS(z_val - z_grid(i))
      IF (dist < min_dist) THEN
        min_dist = dist
        i_z = i
      END IF
    END DO
  END SUBROUTINE find_nearest_point

  !---------------------------------------------------------------------------
  ! Compute derivatives at a grid point
  !---------------------------------------------------------------------------
  SUBROUTINE compute_derivatives_at_point(f, i, j, k, nx, ny, nz, dx, dy, dz, &
                                          df_dx, df_dy, df_dz)
    REAL(dp), INTENT(IN) :: f(:,:,:)
    INTEGER, INTENT(IN) :: i, j, k, nx, ny, nz
    REAL(dp), INTENT(IN) :: dx, dy, dz
    REAL(dp), INTENT(OUT) :: df_dx, df_dy, df_dz
    
    ! Compute x derivative using centered differences where possible
    IF (i > 1 .AND. i < nx) THEN
      df_dx = (f(i+1,j,k) - f(i-1,j,k)) / (2.0_dp * dx)
    ELSE IF (i == 1) THEN
      df_dx = (f(i+1,j,k) - f(i,j,k)) / dx
    ELSE ! i == nx
      df_dx = (f(i,j,k) - f(i-1,j,k)) / dx
    END IF
    
    ! Compute y derivative using centered differences where possible
    IF (j > 1 .AND. j < ny) THEN
      df_dy = (f(i,j+1,k) - f(i,j-1,k)) / (2.0_dp * dy)
    ELSE IF (j == 1) THEN
      df_dy = (f(i,j+1,k) - f(i,j,k)) / dy
    ELSE ! j == ny
      df_dy = (f(i,j,k) - f(i,j-1,k)) / dy
    END IF
    
    ! Compute z derivative using centered differences where possible
    IF (k > 1 .AND. k < nz) THEN
      df_dz = (f(i,j,k+1) - f(i,j,k-1)) / (2.0_dp * dz)
    ELSE IF (k == 1) THEN
      df_dz = (f(i,j,k+1) - f(i,j,k)) / dz
    ELSE ! k == nz
      df_dz = (f(i,j,k) - f(i,j,k-1)) / dz
    END IF
  END SUBROUTINE compute_derivatives_at_point


  !---------------------------------------------------------------------------
  ! Hermite basis functions
  !---------------------------------------------------------------------------
  FUNCTION h00(t) RESULT(h)
    REAL(dp), INTENT(IN) :: t
    REAL(dp) :: h
    h = 2.0_dp*t**3 - 3.0_dp*t**2 + 1.0_dp
  END FUNCTION h00
  
  FUNCTION h10(t) RESULT(h)
    REAL(dp), INTENT(IN) :: t
    REAL(dp) :: h
    h = t**3 - 2.0_dp*t**2 + t
  END FUNCTION h10
  
  FUNCTION h01(t) RESULT(h)
    REAL(dp), INTENT(IN) :: t
    REAL(dp) :: h
    h = -2.0_dp*t**3 + 3.0_dp*t**2
  END FUNCTION h01
  
  FUNCTION h11(t) RESULT(h)
    REAL(dp), INTENT(IN) :: t
    REAL(dp) :: h
    h = t**3 - t**2
  END FUNCTION h11












END MODULE hermite_interp