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

  SUBROUTINE constructsed_linear(teff, log_g, metallicity, R, d, file_names, lu_teff, lu_logg, lu_meta, stellar_model_dir, wavelengths, fluxes)
    REAL(dp), INTENT(IN) :: teff, log_g, metallicity, R, d
    REAL(dp), INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
    CHARACTER(LEN=*), INTENT(IN) :: stellar_model_dir
    CHARACTER(LEN=100), INTENT(IN) :: file_names(:)
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, fluxes

    INTEGER :: i, n_lambda, status, n_teff, n_logg, n_meta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: interp_flux, diluted_flux
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: precomputed_flux_cube
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: flux_cube_lambda
    REAL(dp) :: min_flux, max_flux, mean_flux, progress_pct

REAL(dp), ALLOCATABLE :: log_flux_cube(:,:,:), flux_slice(:,:,:)
REAL(dp) :: interp_log_flux
    ! Parameter grids
    REAL(dp), ALLOCATABLE :: teff_grid(:), logg_grid(:), meta_grid(:)
    CHARACTER(LEN=256) :: bin_filename, clean_path
    LOGICAL :: file_exists

    ! Clean up any double slashes in the path
    clean_path = TRIM(stellar_model_dir)
    IF (clean_path(LEN_TRIM(clean_path):LEN_TRIM(clean_path)) == '/') THEN
      bin_filename = TRIM(clean_path) // 'flux_cube.bin'
    ELSE
      bin_filename = TRIM(clean_path) // '/flux_cube.bin'
    END IF
    
    !PRINT *, '============================================================'
    !PRINT *, 'Loading precomputed flux cube from:', TRIM(bin_filename)
    !PRINT *, 'Target parameters:'
    !PRINT *, '  Teff    =', teff
    !PRINT *, '  log_g   =', log_g
    !PRINT *, '  [M/H]   =', metallicity
    !PRINT *, '  Radius  =', R, ' (solar radii)'
    !PRINT *, '  Distance=', d, ' (cm)'
    
    ! Check if file exists first
    INQUIRE(FILE=bin_filename, EXIST=file_exists)
    
    IF (.NOT. file_exists) THEN
      !PRINT *, 'ERROR: Required binary file not found:', TRIM(bin_filename)
      !PRINT *, 'Please run the precompute_flux_cube.py script to generate this file.'
      !PRINT *, 'Sample command: python precompute_flux_cube.py --model_dir=', TRIM(stellar_model_dir)
      STOP 'Missing required binary file for interpolation'
    END IF
    
    ! Load the data from binary file
    CALL load_binary_data(bin_filename, teff_grid, logg_grid, meta_grid, &
                         wavelengths, precomputed_flux_cube, status)
                         
    IF (status /= 0) THEN
      !PRINT *, 'ERROR: Failed to load binary data from:', TRIM(bin_filename)
      !PRINT *, 'The file exists but may be corrupted or in the wrong format.'
      STOP 'Binary data loading error'
    END IF
    
    n_teff = SIZE(teff_grid)
    n_logg = SIZE(logg_grid)
    n_meta = SIZE(meta_grid)
    n_lambda = SIZE(wavelengths)
    
    !PRINT *, 'Flux cube dimensions:', n_teff, 'x', n_logg, 'x', n_meta, 'x', n_lambda
    !PRINT *, 'Parameter ranges:'
    !PRINT *, '  Teff:     ', teff, ' in range [', teff_grid(1), ',', teff_grid(n_teff), ']'
    !PRINT *, '  log_g:    ', log_g, ' in range [', logg_grid(1), ',', logg_grid(n_logg), ']'
    !PRINT *, '  [M/H]:    ', metallicity, ' in range [', meta_grid(1), ',', meta_grid(n_meta), ']'
    !PRINT *, 'Wavelength range: [', wavelengths(1), ',', wavelengths(n_lambda), '] Å'
    
    ! Allocate space for interpolated flux
    ALLOCATE(interp_flux(n_lambda))
    


    ALLOCATE(log_flux_cube(SIZE(teff_grid), SIZE(logg_grid), SIZE(meta_grid)))
    ALLOCATE(flux_slice(SIZE(teff_grid), SIZE(logg_grid), SIZE(meta_grid)))

    DO i = 1, n_lambda
      flux_slice = f_cube(:, :, :, i)
      log_flux_cube = LOG10(MAX(flux_slice, 1d-99))  ! safe log10
      interp_log_flux = trilinear_interp(teff, log_g, metallicity, &
                                          teff_grid, logg_grid, meta_grid, log_flux_cube)
      interp_flux(i) = 10.0_dp ** interp_log_flux
    END DO

    
    ! Calculate statistics for validation
    min_flux = MINVAL(interp_flux)
    max_flux = MAXVAL(interp_flux)
    mean_flux = SUM(interp_flux) / n_lambda
    
    !PRINT *, 'Interpolation completed successfully.'
    !PRINT *, 'Interpolated flux statistics:'
    !PRINT *, '  Min flux:', min_flux
    !PRINT *, '  Max flux:', max_flux
    !PRINT *, '  Mean flux:', mean_flux
    
    ! Apply distance dilution to get observed flux
    ALLOCATE(diluted_flux(n_lambda))
    CALL dilute_flux(interp_flux, R, d, diluted_flux)
    fluxes = diluted_flux
    
    ! Calculate statistics after dilution
    min_flux = MINVAL(diluted_flux)
    max_flux = MAXVAL(diluted_flux)
    mean_flux = SUM(diluted_flux) / n_lambda
    
    !PRINT *, 'Diluted flux statistics (observed from Earth):'
    !PRINT *, '  Min flux:', min_flux
    !PRINT *, '  Max flux:', max_flux
    !PRINT *, '  Mean flux:', mean_flux
    
    ! Clean up
    DEALLOCATE(teff_grid, logg_grid, meta_grid, precomputed_flux_cube)
    DEALLOCATE(diluted_flux, interp_flux)
    
    !PRINT *, 'SED construction complete'
    !PRINT *, '============================================================'
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
      !PRINT *, 'Error opening binary file:', TRIM(filename)
      RETURN
    END IF
    
    ! Read dimensions
    READ(unit, IOSTAT=status) n_teff, n_logg, n_meta, n_lambda
    IF (status /= 0) THEN
      !PRINT *, 'Error reading dimensions from binary file'
      CLOSE(unit)
      RETURN
    END IF
    
    !PRINT *, 'Read dimensions from file:', n_teff, n_logg, n_meta, n_lambda
    
    ! Allocate arrays based on dimensions
    ALLOCATE(teff_grid(n_teff), STAT=status)
    IF (status /= 0) THEN
      !PRINT *, 'Error allocating teff_grid array'
      CLOSE(unit)
      RETURN
    END IF
    
    ALLOCATE(logg_grid(n_logg), STAT=status)
    IF (status /= 0) THEN
      !PRINT *, 'Error allocating logg_grid array'
      DEALLOCATE(teff_grid)
      CLOSE(unit)
      RETURN
    END IF
    
    ALLOCATE(meta_grid(n_meta), STAT=status)
    IF (status /= 0) THEN
      !PRINT *, 'Error allocating meta_grid array'
      DEALLOCATE(teff_grid, logg_grid)
      CLOSE(unit)
      RETURN
    END IF
    
    ALLOCATE(wavelengths(n_lambda), STAT=status)
    IF (status /= 0) THEN
      !PRINT *, 'Error allocating wavelengths array'
      DEALLOCATE(teff_grid, logg_grid, meta_grid)
      CLOSE(unit)
      RETURN
    END IF
    
    ALLOCATE(flux_cube(n_teff, n_logg, n_meta, n_lambda), STAT=status)
    IF (status /= 0) THEN
      !PRINT *, 'Error allocating flux_cube array'
      DEALLOCATE(teff_grid, logg_grid, meta_grid, wavelengths)
      CLOSE(unit)
      RETURN
    END IF
    
    ! Read grid arrays
    READ(unit, IOSTAT=status) teff_grid
    IF (status /= 0) THEN
      !PRINT *, 'Error reading teff_grid'
      GOTO 999  ! Cleanup and return
    END IF
    
    READ(unit, IOSTAT=status) logg_grid
    IF (status /= 0) THEN
      !PRINT *, 'Error reading logg_grid'
      GOTO 999  ! Cleanup and return
    END IF
    
    READ(unit, IOSTAT=status) meta_grid
    IF (status /= 0) THEN
      !PRINT *, 'Error reading meta_grid'
      GOTO 999  ! Cleanup and return
    END IF
    
    READ(unit, IOSTAT=status) wavelengths
    IF (status /= 0) THEN
      !PRINT *, 'Error reading wavelengths'
      GOTO 999  ! Cleanup and return
    END IF
    
    ! Read flux cube
    READ(unit, IOSTAT=status) flux_cube
    IF (status /= 0) THEN
      !PRINT *, 'Error reading flux_cube'
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

! After reading the grid arrays
!PRINT *, 'Teff grid min/max:', MINVAL(teff_grid), MAXVAL(teff_grid)
!PRINT *, 'logg grid min/max:', MINVAL(logg_grid), MAXVAL(logg_grid)
!PRINT *, 'meta grid min/max:', MINVAL(meta_grid), MAXVAL(meta_grid)    

  END SUBROUTINE load_binary_data

  !---------------------------------------------------------------------------
  ! Simple trilinear interpolation function
  !---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Log-space trilinear interpolation function with normalization
!---------------------------------------------------------------------------
FUNCTION trilinear_interp(x_val, y_val, z_val, x_grid, y_grid, z_grid, f_values) RESULT(f_interp)
  REAL(dp), INTENT(IN) :: x_val, y_val, z_val
  REAL(dp), INTENT(IN) :: x_grid(:), y_grid(:), z_grid(:)
  REAL(dp), INTENT(IN) :: f_values(:,:,:)
  REAL(dp) :: f_interp
      ! Compute log-space result
  REAL(dp) :: log_result
  INTEGER :: i_x, i_y, i_z
  REAL(dp) :: t_x, t_y, t_z
  REAL(dp) :: c000, c001, c010, c011, c100, c101, c110, c111
  REAL(dp) :: c00, c01, c10, c11, c0, c1
  REAL(dp) :: tiny_value = 1.0e-10_dp  ! Increased from 1e-30 to avoid underflow
  
  ! Find containing cell and parameter values using binary search
  CALL find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                            i_x, i_y, i_z, t_x, t_y, t_z)
  
  ! Boundary safety check
  IF (i_x < 1) i_x = 1
  IF (i_y < 1) i_y = 1
  IF (i_z < 1) i_z = 1
  IF (i_x >= SIZE(x_grid)) i_x = SIZE(x_grid)-1
  IF (i_y >= SIZE(y_grid)) i_y = SIZE(y_grid)-1
  IF (i_z >= SIZE(z_grid)) i_z = SIZE(z_grid)-1
  
  ! Force interpolation parameters to be in [0,1]
  t_x = MAX(0.0_dp, MIN(1.0_dp, t_x))
  t_y = MAX(0.0_dp, MIN(1.0_dp, t_y))
  t_z = MAX(0.0_dp, MIN(1.0_dp, t_z))
  
  ! Get the corners of the cube with safety checks
  c000 = MAX(tiny_value, f_values(i_x,   i_y,   i_z))
  c001 = MAX(tiny_value, f_values(i_x,   i_y,   i_z+1))
  c010 = MAX(tiny_value, f_values(i_x,   i_y+1, i_z))
  c011 = MAX(tiny_value, f_values(i_x,   i_y+1, i_z+1))
  c100 = MAX(tiny_value, f_values(i_x+1, i_y,   i_z))
  c101 = MAX(tiny_value, f_values(i_x+1, i_y,   i_z+1))
  c110 = MAX(tiny_value, f_values(i_x+1, i_y+1, i_z))
  c111 = MAX(tiny_value, f_values(i_x+1, i_y+1, i_z+1))
  
  ! Try standard linear interpolation first (safer)
  c00 = c000 * (1.0_dp - t_x) + c100 * t_x
  c01 = c001 * (1.0_dp - t_x) + c101 * t_x
  c10 = c010 * (1.0_dp - t_x) + c110 * t_x
  c11 = c011 * (1.0_dp - t_x) + c111 * t_x
  
  c0 = c00 * (1.0_dp - t_y) + c10 * t_y
  c1 = c01 * (1.0_dp - t_y) + c11 * t_y
  
  f_interp = c0 * (1.0_dp - t_z) + c1 * t_z
  
  ! If the linear result is valid and non-zero, try log space
  IF (f_interp > tiny_value) THEN
    ! Perform log-space interpolation
    c00 = LOG(c000) * (1.0_dp - t_x) + LOG(c100) * t_x
    c01 = LOG(c001) * (1.0_dp - t_x) + LOG(c101) * t_x
    c10 = LOG(c010) * (1.0_dp - t_x) + LOG(c110) * t_x
    c11 = LOG(c011) * (1.0_dp - t_x) + LOG(c111) * t_x
    
    c0 = c00 * (1.0_dp - t_y) + c10 * t_y
    c1 = c01 * (1.0_dp - t_y) + c11 * t_y
    

    log_result = c0 * (1.0_dp - t_z) + c1 * t_z
    
    ! Only use the log-space result if it's valid
    IF (log_result == log_result) THEN  ! NaN check
      f_interp = EXP(log_result)
    END IF
  END IF
  
  ! Final sanity check
  IF (f_interp /= f_interp .OR. f_interp <= 0.0_dp) THEN
    ! If we somehow still got an invalid result, use nearest neighbor
    CALL find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, i_x, i_y, i_z)
    f_interp = MAX(tiny_value, f_values(i_x, i_y, i_z))
  END IF
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