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

! ***********************************************************************
! Hermite interpolation module for spectral energy distributions (SEDs)
! ***********************************************************************

MODULE hermite_interp
  USE const_def, ONLY: dp
  USE shared_funcs, only: dilute_flux
  implicit none

  PRIVATE
  PUBLIC :: construct_sed_hermite, hermite_tensor_interp3d

CONTAINS

  !---------------------------------------------------------------------------
  ! Main entry point: Construct a SED using Hermite tensor interpolation
  !---------------------------------------------------------------------------
SUBROUTINE construct_sed_hermite(teff, log_g, metallicity, R, d, file_names, &
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

    ! Load the data from binary file
    CALL load_binary_data(bin_filename, teff_grid, logg_grid, meta_grid, &
                         wavelengths, precomputed_flux_cube, status)

    n_teff = SIZE(teff_grid)
    n_logg = SIZE(logg_grid)
    n_meta = SIZE(meta_grid)
    n_lambda = SIZE(wavelengths)

    ! Allocate space for interpolated flux
    ALLOCATE(interp_flux(n_lambda))

    ! Process each wavelength point
    DO i = 1, n_lambda
      ALLOCATE(flux_cube_lambda(n_teff, n_logg, n_meta))
      flux_cube_lambda = precomputed_flux_cube(:,:,:,i)

      interp_flux(i) = hermite_tensor_interp3d(teff, log_g, metallicity, &
                                              teff_grid, logg_grid, meta_grid, flux_cube_lambda)

      DEALLOCATE(flux_cube_lambda)
    END DO

    ! Apply distance dilution to get observed flux
    ALLOCATE(diluted_flux(n_lambda))
    CALL dilute_flux(interp_flux, R, d, diluted_flux)
    fluxes = diluted_flux

  END SUBROUTINE construct_sed_hermite


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
    CLOSE(unit)
    RETURN
  END IF

  ALLOCATE(meta_grid(n_meta), STAT=status)
  IF (status /= 0) THEN
    PRINT *, 'Error allocating meta_grid array'
    CLOSE(unit)
    RETURN
  END IF

  ALLOCATE(wavelengths(n_lambda), STAT=status)
  IF (status /= 0) THEN
    PRINT *, 'Error allocating wavelengths array'
    CLOSE(unit)
    RETURN
  END IF

  ALLOCATE(flux_cube(n_teff, n_logg, n_meta, n_lambda), STAT=status)
  IF (status /= 0) THEN
    PRINT *, 'Error allocating flux_cube array'
    CLOSE(unit)
    RETURN
  END IF

  ! Read grid arrays
  READ(unit, IOSTAT=status) teff_grid
  IF (status /= 0) THEN
    PRINT *, 'Error reading teff_grid'
    CLOSE(unit)
    RETURN
  END IF

  READ(unit, IOSTAT=status) logg_grid
  IF (status /= 0) THEN
    PRINT *, 'Error reading logg_grid'
    CLOSE(unit)
    RETURN
  END IF

  READ(unit, IOSTAT=status) meta_grid
  IF (status /= 0) THEN
    PRINT *, 'Error reading meta_grid'
    CLOSE(unit)
    RETURN
  END IF

  READ(unit, IOSTAT=status) wavelengths
  IF (status /= 0) THEN
    PRINT *, 'Error reading wavelengths'
    CLOSE(unit)
    RETURN
  END IF

  ! Read flux cube
  READ(unit, IOSTAT=status) flux_cube
  IF (status /= 0) THEN
    PRINT *, 'Error reading flux_cube'
    CLOSE(unit)
    RETURN
  END IF

  ! Close file and return success
  CLOSE(unit)
END SUBROUTINE load_binary_data


!---------------------------------------------------------------------------
  ! Main 3D Hermite interpolation function
  !---------------------------------------------------------------------------
FUNCTION hermite_tensor_interp3d(x_val, y_val, z_val, x_grid, y_grid, &
                                z_grid, f_values) RESULT(f_interp)
    REAL(dp), INTENT(IN) :: x_val, y_val, z_val
    REAL(dp), INTENT(IN) :: x_grid(:), y_grid(:), z_grid(:)
    REAL(dp), INTENT(IN) :: f_values(:,:,:)
    REAL(dp) :: f_interp

    INTEGER :: i_x, i_y, i_z
    REAL(dp) :: t_x, t_y, t_z
    REAL(dp) :: dx, dy, dz
    REAL(dp) :: dx_values(2,2,2), dy_values(2,2,2), dz_values(2,2,2)
    REAL(dp) :: values(2,2,2)
    REAL(dp) :: sum
    INTEGER :: ix, iy, iz
    REAL(dp) :: h_x(2), h_y(2), h_z(2)

    ! Find containing cell and parameter values
    CALL find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                              i_x, i_y, i_z, t_x, t_y, t_z)

    ! If outside grid, use nearest point
    IF (i_x < 1 .OR. i_x >= SIZE(x_grid) .OR. &
        i_y < 1 .OR. i_y >= SIZE(y_grid) .OR. &
        i_z < 1 .OR. i_z >= SIZE(z_grid)) THEN

      CALL find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                             i_x, i_y, i_z)
      f_interp = f_values(i_x, i_y, i_z)
      RETURN
    END IF

    ! Grid cell spacing
    dx = x_grid(i_x+1) - x_grid(i_x)
    dy = y_grid(i_y+1) - y_grid(i_y)
    dz = z_grid(i_z+1) - z_grid(i_z)

    ! Extract the local 2x2x2 grid cell and compute derivatives
    DO iz = 0, 1
      DO iy = 0, 1
        DO ix = 0, 1
          values(ix+1, iy+1, iz+1) = f_values(i_x+ix, i_y+iy, i_z+iz)
          CALL compute_derivatives_at_point(f_values, i_x+ix, i_y+iy, i_z+iz, &
                                           SIZE(x_grid), SIZE(y_grid), SIZE(z_grid), &
                                           dx, dy, dz, &
                                           dx_values(ix+1,iy+1,iz+1), &
                                           dy_values(ix+1,iy+1,iz+1), &
                                           dz_values(ix+1,iy+1,iz+1))
        END DO
      END DO
    END DO

    sum = 0.0_dp

    ! Function values
    h_x = [h00(t_x), h01(t_x)]
    h_y = [h00(t_y), h01(t_y)]
    h_z = [h00(t_z), h01(t_z)]
    DO iz = 1, 2
      DO iy = 1, 2
        DO ix = 1, 2
          sum = sum + h_x(ix) * h_y(iy) * h_z(iz) * values(ix,iy,iz)
        END DO
      END DO
    END DO

    ! x-derivatives
    h_x = [h10(t_x), h11(t_x)]
    DO iz = 1, 2
      DO iy = 1, 2
        DO ix = 1, 2
          sum = sum + h_x(ix) * h_y(iy) * h_z(iz) * dx * dx_values(ix,iy,iz)
        END DO
      END DO
    END DO

    ! y-derivatives
    h_x = [h00(t_x), h01(t_x)]
    h_y = [h10(t_y), h11(t_y)]
    DO iz = 1, 2
      DO iy = 1, 2
        DO ix = 1, 2
          sum = sum + h_x(ix) * h_y(iy) * h_z(iz) * dy * dy_values(ix,iy,iz)
        END DO
      END DO
    END DO

    ! z-derivatives
    h_y = [h00(t_y), h01(t_y)]
    h_z = [h10(t_z), h11(t_z)]
    DO iz = 1, 2
      DO iy = 1, 2
        DO ix = 1, 2
          sum = sum + h_x(ix) * h_y(iy) * h_z(iz) * dz * dz_values(ix,iy,iz)
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
