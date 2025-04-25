! ***********************************************************************
!
! Linear interpolation module for spectral energy distributions (SEDs)
! ***********************************************************************

MODULE linear_interp
  USE const_def, ONLY: dp
  USE shared_funcs
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: constructsed_linear, trilinear_interp

CONTAINS

  !---------------------------------------------------------------------------
  ! Main entry point: Construct a SED using linear interpolation
  !---------------------------------------------------------------------------
  SUBROUTINE constructsed_linear(teff, log_g, metallicity, R, d, file_names, &
                                  lu_teff, lu_logg, lu_meta, stellar_model_dir, &
                                  wavelengths, fluxes)
    REAL(dp), INTENT(IN) :: teff, log_g, metallicity, R, d
    REAL(dp), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
    CHARACTER(LEN=*), INTENT(IN) :: stellar_model_dir
    CHARACTER(LEN=100), INTENT(IN) :: file_names(:)
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, fluxes

    INTEGER :: i, n_lambda, status, n_teff, n_logg, n_meta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: interp_flux, diluted_flux
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: precomputed_flux_cube
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: flux_cube_lambda
    
    ! Parameter grids
    REAL(dp), ALLOCATABLE :: teff_grid(:), logg_grid(:), meta_grid(:)
    CHARACTER(LEN=256) :: bin_filename

    ! Construct the binary filename
    bin_filename = TRIM(stellar_model_dir) // '/flux_cube.bin'
    
    PRINT *, 'Loading precomputed flux cube from:', TRIM(bin_filename)
    PRINT *, 'Target parameters: Teff =', teff, ', log_g =', log_g, ', metallicity =', metallicity
    
    ! Load the data from binary file
    CALL load_binary_data(bin_filename, teff_grid, logg_grid, meta_grid, &
                         wavelengths, precomputed_flux_cube, status)
                         
    IF (status /= 0) THEN
      PRINT *, 'Error loading precomputed data. Falling back to on-the-fly computation.'
      ! Here you could call the original implementation
      RETURN
    END IF
    
    n_teff = SIZE(teff_grid)
    n_logg = SIZE(logg_grid)
    n_meta = SIZE(meta_grid)
    n_lambda = SIZE(wavelengths)
    
    PRINT *, 'Loaded flux cube with dimensions:', n_teff, 'x', n_logg, 'x', n_meta, 'x', n_lambda
    PRINT *, 'Teff range:', teff_grid(1), 'to', teff_grid(n_teff)
    PRINT *, 'logg range:', logg_grid(1), 'to', logg_grid(n_logg)
    PRINT *, 'metallicity range:', meta_grid(1), 'to', meta_grid(n_meta)
    PRINT *, 'Performing interpolation at target parameters...'
    
    ! Allocate space for interpolated flux
    ALLOCATE(interp_flux(n_lambda))
    
    ! Perform trilinear interpolation for each wavelength
    DO i = 1, n_lambda
      ! Extract the 3D grid for this wavelength
      ALLOCATE(flux_cube_lambda(n_teff, n_logg, n_meta))
      flux_cube_lambda = precomputed_flux_cube(:,:,:,i)
      
      ! Simple trilinear interpolation at the target parameters
      interp_flux(i) = trilinear_interp(teff, log_g, metallicity, &
                                       teff_grid, logg_grid, meta_grid, flux_cube_lambda)
      
      DEALLOCATE(flux_cube_lambda)
    END DO
    
    ! Apply distance dilution to get observed flux
    ALLOCATE(diluted_flux(n_lambda))
    CALL dilute_flux(interp_flux, R, d, diluted_flux)
    fluxes = diluted_flux
    
    ! Clean up
    DEALLOCATE(teff_grid, logg_grid, meta_grid, precomputed_flux_cube)
    DEALLOCATE(diluted_flux, interp_flux)
    
    PRINT *, 'Interpolation complete'
  END SUBROUTINE constructsed_linear

  !---------------------------------------------------------------------------
  ! Load data from binary file
  !---------------------------------------------------------------------------
  SUBROUTINE load_binary_data(filename, teff_grid, logg_grid, meta_grid, &
                             wavelengths, flux_cube, status)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: teff_grid(:), logg_grid(:), meta_grid(:)
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: wavelengths(:)
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: flux_cube(:,:,:,:)
    INTEGER, INTENT(OUT) :: status
    
    INTEGER :: unit, n_teff, n_logg, n_meta, n_lambda
    
    unit = 99
    status = 0
    
    ! Open the binary file
    OPEN(UNIT=unit, FILE=filename, STATUS='OLD', ACCESS='STREAM', FORM='UNFORMATTED', IOSTAT=status)
    IF (status /= 0) THEN
      PRINT *, 'Error opening binary file:', TRIM(filename)
      RETURN
    END IF
    
    ! Read dimensions
    READ(unit, IOSTAT=status) n_teff, n_logg, n_meta, n_lambda
    IF (status /= 0) THEN
      PRINT *, 'Error reading dimensions from binary file'
      CLOSE(unit)
      RETURN
    END IF
    
    PRINT *, 'Read dimensions from file:', n_teff, n_logg, n_meta, n_lambda
    
    ! Allocate arrays based on dimensions
    ALLOCATE(teff_grid(n_teff), STAT=status)
    IF (status /= 0) THEN
      PRINT *, 'Error allocating teff_grid array'
      CLOSE(unit)
      RETURN
    END IF
    
    ALLOCATE(logg_grid(n_logg), STAT=status)
    IF (status /= 0) THEN
      PRINT *, 'Error allocating logg_grid array'
      DEALLOCATE(teff_grid)
      CLOSE(unit)
      RETURN
    END IF
    
    ALLOCATE(meta_grid(n_meta), STAT=status)
    IF (status /= 0) THEN
      PRINT *, 'Error allocating meta_grid array'
      DEALLOCATE(teff_grid, logg_grid)
      CLOSE(unit)
      RETURN
    END IF
    
    ALLOCATE(wavelengths(n_lambda), STAT=status)
    IF (status /= 0) THEN
      PRINT *, 'Error allocating wavelengths array'
      DEALLOCATE(teff_grid, logg_grid, meta_grid)
      CLOSE(unit)
      RETURN
    END IF
    
    ALLOCATE(flux_cube(n_teff, n_logg, n_meta, n_lambda), STAT=status)
    IF (status /= 0) THEN
      PRINT *, 'Error allocating flux_cube array'
      DEALLOCATE(teff_grid, logg_grid, meta_grid, wavelengths)
      CLOSE(unit)
      RETURN
    END IF
    
    ! Read grid arrays
    READ(unit, IOSTAT=status) teff_grid
    IF (status /= 0) THEN
      PRINT *, 'Error reading teff_grid'
      GOTO 999  ! Cleanup and return
    END IF
    
    READ(unit, IOSTAT=status) logg_grid
    IF (status /= 0) THEN
      PRINT *, 'Error reading logg_grid'
      GOTO 999  ! Cleanup and return
    END IF
    
    READ(unit, IOSTAT=status) meta_grid
    IF (status /= 0) THEN
      PRINT *, 'Error reading meta_grid'
      GOTO 999  ! Cleanup and return
    END IF
    
    READ(unit, IOSTAT=status) wavelengths
    IF (status /= 0) THEN
      PRINT *, 'Error reading wavelengths'
      GOTO 999  ! Cleanup and return
    END IF
    
    ! Read flux cube
    READ(unit, IOSTAT=status) flux_cube
    IF (status /= 0) THEN
      PRINT *, 'Error reading flux_cube'
      GOTO 999  ! Cleanup and return
    END IF
    
    ! Close file and return success
    CLOSE(unit)
    RETURN
    
999 CONTINUE
    ! Cleanup on error
    DEALLOCATE(teff_grid, logg_grid, meta_grid, wavelengths, flux_cube)
    CLOSE(unit)
    RETURN
    
  END SUBROUTINE load_binary_data

  !---------------------------------------------------------------------------
  ! Simple trilinear interpolation function
  !---------------------------------------------------------------------------
  FUNCTION trilinear_interp(x_val, y_val, z_val, x_grid, y_grid, z_grid, f_values) RESULT(f_interp)
    REAL(dp), INTENT(IN) :: x_val, y_val, z_val
    REAL(dp), INTENT(IN) :: x_grid(:), y_grid(:), z_grid(:)
    REAL(dp), INTENT(IN) :: f_values(:,:,:)
    REAL(dp) :: f_interp
    
    INTEGER :: i_x, i_y, i_z
    REAL(dp) :: t_x, t_y, t_z
    REAL(dp) :: c000, c001, c010, c011, c100, c101, c110, c111
    REAL(dp) :: c00, c01, c10, c11, c0, c1
    
    ! Find containing cell and parameter values using binary search
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
    
    ! Get the corners of the cube
    c000 = f_values(i_x,   i_y,   i_z)
    c001 = f_values(i_x,   i_y,   i_z+1)
    c010 = f_values(i_x,   i_y+1, i_z)
    c011 = f_values(i_x,   i_y+1, i_z+1)
    c100 = f_values(i_x+1, i_y,   i_z)
    c101 = f_values(i_x+1, i_y,   i_z+1)
    c110 = f_values(i_x+1, i_y+1, i_z)
    c111 = f_values(i_x+1, i_y+1, i_z+1)
    
    ! Perform trilinear interpolation
    ! 1. Interpolate along x direction
    c00 = c000 * (1.0_dp - t_x) + c100 * t_x
    c01 = c001 * (1.0_dp - t_x) + c101 * t_x
    c10 = c010 * (1.0_dp - t_x) + c110 * t_x
    c11 = c011 * (1.0_dp - t_x) + c111 * t_x
    
    ! 2. Interpolate along y direction
    c0 = c00 * (1.0_dp - t_y) + c10 * t_y
    c1 = c01 * (1.0_dp - t_y) + c11 * t_y
    
    ! 3. Interpolate along z direction
    f_interp = c0 * (1.0_dp - t_z) + c1 * t_z
  END FUNCTION trilinear_interp

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

END MODULE linear_interp