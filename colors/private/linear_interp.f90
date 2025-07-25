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
! Linear interpolation module for spectral energy distributions (SEDs)
! ***********************************************************************

module linear_interp
   use const_def, only: dp
   use colors_utils, only: dilute_flux
   implicit none

   private
   public :: construct_sed_linear, trilinear_interp

contains

   !---------------------------------------------------------------------------
   ! Main entry point: Construct a SED using linear interpolation
   !---------------------------------------------------------------------------

   subroutine construct_sed_linear(teff, log_g, metallicity, R, d, file_names, &
                                   lu_teff, lu_logg, lu_meta, stellar_model_dir, &
                                   wavelengths, fluxes)

      real(dp), intent(in) :: teff, log_g, metallicity, R, d
      real(dp), intent(in) :: lu_teff(:), lu_logg(:), lu_meta(:)
      character(len=*), intent(in) :: stellar_model_dir
      character(len=100), intent(in) :: file_names(:)
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, fluxes

      integer :: i, n_lambda, status, n_teff, n_logg, n_meta
      real(dp), dimension(:), allocatable :: interp_flux, diluted_flux
      real(dp), dimension(:, :, :, :), allocatable :: precomputed_flux_cube
      real(dp), dimension(:, :, :), allocatable :: flux_cube_lambda
      real(dp) :: min_flux, max_flux, mean_flux, progress_pct

      ! Parameter grids
      real(dp), allocatable :: teff_grid(:), logg_grid(:), meta_grid(:)
      character(len=256) :: bin_filename, clean_path
      logical :: file_exists

      ! Clean up any double slashes in the path
      clean_path = trim(stellar_model_dir)
      if (clean_path(len_trim(clean_path):len_trim(clean_path)) == '/') then
         bin_filename = trim(clean_path)//'flux_cube.bin'
      else
         bin_filename = trim(clean_path)//'/flux_cube.bin'
      end if

      ! Check if file exists first
      INQUIRE (file=bin_filename, EXIST=file_exists)

      if (.not. file_exists) then
         stop 'Missing required binary file for interpolation'
      end if

      ! Load the data from binary file
      call load_binary_data(bin_filename, teff_grid, logg_grid, meta_grid, &
                            wavelengths, precomputed_flux_cube, status)

      if (status /= 0) then
         stop 'Binary data loading error'
      end if

      n_teff = size(teff_grid)
      n_logg = size(logg_grid)
      n_meta = size(meta_grid)
      n_lambda = size(wavelengths)

      ! Allocate space for interpolated flux
      allocate (interp_flux(n_lambda))

      ! Perform trilinear interpolation for each wavelength
      do i = 1, n_lambda
         ! !print progress updates at regular intervals

         ! Extract the 3D grid for this wavelength
         allocate (flux_cube_lambda(n_teff, n_logg, n_meta))
         flux_cube_lambda = precomputed_flux_cube(:, :, :, i)

         ! Simple trilinear interpolation at the target parameters
         interp_flux(i) = trilinear_interp(teff, log_g, metallicity, &
                                           teff_grid, logg_grid, meta_grid, flux_cube_lambda)

      end do

      ! Calculate statistics for validation
      min_flux = minval(interp_flux)
      max_flux = maxval(interp_flux)
      mean_flux = sum(interp_flux)/n_lambda

      ! Apply distance dilution to get observed flux
      allocate (diluted_flux(n_lambda))
      call dilute_flux(interp_flux, R, d, diluted_flux)
      fluxes = diluted_flux

      ! Calculate statistics after dilution
      min_flux = minval(diluted_flux)
      max_flux = maxval(diluted_flux)
      mean_flux = sum(diluted_flux)/n_lambda

   end subroutine construct_sed_linear

   !---------------------------------------------------------------------------
   ! Load data from binary file
   !---------------------------------------------------------------------------
   subroutine load_binary_data(filename, teff_grid, logg_grid, meta_grid, &
                               wavelengths, flux_cube, status)
      character(len=*), intent(in) :: filename
      real(dp), allocatable, intent(out) :: teff_grid(:), logg_grid(:), meta_grid(:)
      real(dp), allocatable, intent(out) :: wavelengths(:)
      real(dp), allocatable, intent(out) :: flux_cube(:, :, :, :)
      integer, intent(out) :: status

      integer :: unit, n_teff, n_logg, n_meta, n_lambda

      unit = 99
      status = 0

      ! Open the binary file
      open (unit=unit, file=filename, status='OLD', ACCESS='STREAM', FORM='UNFORMATTED', iostat=status)
      if (status /= 0) then
         !print *, 'Error opening binary file:', trim(filename)
         return
      end if

      ! Read dimensions
      read (unit, iostat=status) n_teff, n_logg, n_meta, n_lambda
      if (status /= 0) then
         !print *, 'Error reading dimensions from binary file'
         close (unit)
         return
      end if

      ! Allocate arrays based on dimensions
      allocate (teff_grid(n_teff), STAT=status)
      if (status /= 0) then
         !print *, 'Error allocating teff_grid array'
         close (unit)
         return
      end if

      allocate (logg_grid(n_logg), STAT=status)
      if (status /= 0) then
         !print *, 'Error allocating logg_grid array'
         close (unit)
         return
      end if

      allocate (meta_grid(n_meta), STAT=status)
      if (status /= 0) then
         !print *, 'Error allocating meta_grid array'
         close (unit)
         return
      end if

      allocate (wavelengths(n_lambda), STAT=status)
      if (status /= 0) then
         !print *, 'Error allocating wavelengths array'
         close (unit)
         return
      end if

      allocate (flux_cube(n_teff, n_logg, n_meta, n_lambda), STAT=status)
      if (status /= 0) then
         !print *, 'Error allocating flux_cube array'
         close (unit)
         return
      end if

      ! Read grid arrays
      read (unit, iostat=status) teff_grid
      if (status /= 0) then
         !print *, 'Error reading teff_grid'
         GOTO 999  ! Cleanup and return
      end if

      read (unit, iostat=status) logg_grid
      if (status /= 0) then
         !print *, 'Error reading logg_grid'
         GOTO 999  ! Cleanup and return
      end if

      read (unit, iostat=status) meta_grid
      if (status /= 0) then
         !print *, 'Error reading meta_grid'
         GOTO 999  ! Cleanup and return
      end if

      read (unit, iostat=status) wavelengths
      if (status /= 0) then
         !print *, 'Error reading wavelengths'
         GOTO 999  ! Cleanup and return
      end if

      ! Read flux cube
      read (unit, iostat=status) flux_cube
      if (status /= 0) then
         !print *, 'Error reading flux_cube'
         GOTO 999  ! Cleanup and return
      end if

      ! Close file and return success
      close (unit)
      return

999   CONTINUE
      ! Cleanup on error
      close (unit)
      return

! After reading the grid arrays
!print *, 'Teff grid min/max:', minval(teff_grid), maxval(teff_grid)
!print *, 'logg grid min/max:', minval(logg_grid), maxval(logg_grid)
!print *, 'meta grid min/max:', minval(meta_grid), maxval(meta_grid)

   end subroutine load_binary_data

   !---------------------------------------------------------------------------
   ! Simple trilinear interpolation function
   !---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Log-space trilinear interpolation function with normalization
!---------------------------------------------------------------------------
   function trilinear_interp(x_val, y_val, z_val, x_grid, y_grid, z_grid, f_values) result(f_interp)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in) :: f_values(:, :, :)
      real(dp) :: f_interp
      ! Compute log-space result
      real(dp) :: log_result
      integer :: i_x, i_y, i_z
      real(dp) :: t_x, t_y, t_z
      real(dp) :: c000, c001, c010, c011, c100, c101, c110, c111
      real(dp) :: c00, c01, c10, c11, c0, c1
      real(dp), parameter :: tiny_value = 1.0e-10_dp

      ! Find containing cell and parameter values using binary search
      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      ! Boundary safety check
      if (i_x < 1) i_x = 1
      if (i_y < 1) i_y = 1
      if (i_z < 1) i_z = 1
      if (i_x >= size(x_grid)) i_x = size(x_grid) - 1
      if (i_y >= size(y_grid)) i_y = size(y_grid) - 1
      if (i_z >= size(z_grid)) i_z = size(z_grid) - 1

      ! Force interpolation parameters to be in [0,1]
      t_x = max(0.0_dp, MIN(1.0_dp, t_x))
      t_y = max(0.0_dp, MIN(1.0_dp, t_y))
      t_z = max(0.0_dp, MIN(1.0_dp, t_z))

      ! Get the corners of the cube with safety checks
      c000 = max(tiny_value, f_values(i_x, i_y, i_z))
      c001 = max(tiny_value, f_values(i_x, i_y, i_z + 1))
      c010 = max(tiny_value, f_values(i_x, i_y + 1, i_z))
      c011 = max(tiny_value, f_values(i_x, i_y + 1, i_z + 1))
      c100 = max(tiny_value, f_values(i_x + 1, i_y, i_z))
      c101 = max(tiny_value, f_values(i_x + 1, i_y, i_z + 1))
      c110 = max(tiny_value, f_values(i_x + 1, i_y + 1, i_z))
      c111 = max(tiny_value, f_values(i_x + 1, i_y + 1, i_z + 1))

      ! Try standard linear interpolation first (safer)
      c00 = c000*(1.0_dp - t_x) + c100*t_x
      c01 = c001*(1.0_dp - t_x) + c101*t_x
      c10 = c010*(1.0_dp - t_x) + c110*t_x
      c11 = c011*(1.0_dp - t_x) + c111*t_x

      c0 = c00*(1.0_dp - t_y) + c10*t_y
      c1 = c01*(1.0_dp - t_y) + c11*t_y

      f_interp = c0*(1.0_dp - t_z) + c1*t_z

      ! If the linear result is valid and non-zero, try log space
      if (f_interp > tiny_value) then
         ! Perform log-space interpolation
         c00 = log(c000)*(1.0_dp - t_x) + log(c100)*t_x
         c01 = log(c001)*(1.0_dp - t_x) + log(c101)*t_x
         c10 = log(c010)*(1.0_dp - t_x) + log(c110)*t_x
         c11 = log(c011)*(1.0_dp - t_x) + log(c111)*t_x

         c0 = c00*(1.0_dp - t_y) + c10*t_y
         c1 = c01*(1.0_dp - t_y) + c11*t_y

         log_result = c0*(1.0_dp - t_z) + c1*t_z

         ! Only use the log-space result if it's valid
         if (log_result == log_result) then  ! NaN check
            f_interp = EXP(log_result)
         end if
      end if

      ! Final sanity check
      if (f_interp /= f_interp .or. f_interp <= 0.0_dp) then
         ! If we somehow still got an invalid result, use nearest neighbor
         call find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, i_x, i_y, i_z)
         f_interp = max(tiny_value, f_values(i_x, i_y, i_z))
      end if
   end function trilinear_interp

   !---------------------------------------------------------------------------
   ! Find the cell containing the interpolation point
   !---------------------------------------------------------------------------
   subroutine find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                   i_x, i_y, i_z, t_x, t_y, t_z)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      integer, intent(out) :: i_x, i_y, i_z
      real(dp), intent(out) :: t_x, t_y, t_z

      ! Find x interval
      call find_interval(x_grid, x_val, i_x, t_x)

      ! Find y interval
      call find_interval(y_grid, y_val, i_y, t_y)

      ! Find z interval
      call find_interval(z_grid, z_val, i_z, t_z)
   end subroutine find_containing_cell

   !---------------------------------------------------------------------------
   ! Find the interval in a sorted array containing a value
   !---------------------------------------------------------------------------
   subroutine find_interval(x, val, i, t)
      real(dp), intent(in) :: x(:), val
      integer, intent(out) :: i
      real(dp), intent(out) :: t

      integer :: n, lo, hi, mid

      n = size(x)

      ! Handle out-of-bounds cases
      if (val <= x(1)) then
         i = 1
         t = 0.0_dp
         return
      else if (val >= x(n)) then
         i = n - 1
         t = 1.0_dp
         return
      end if

      ! Binary search to find interval
      lo = 1
      hi = n
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (val >= x(mid)) then
            lo = mid
         else
            hi = mid
         end if
      end do

      i = lo
      t = (val - x(i))/(x(i + 1) - x(i))
   end subroutine find_interval

   !---------------------------------------------------------------------------
   ! Find the nearest grid point
   !---------------------------------------------------------------------------
   subroutine find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      integer, intent(out) :: i_x, i_y, i_z

      integer :: i
      real(dp) :: min_dist, dist

      ! Find nearest x grid point
      min_dist = abs(x_val - x_grid(1))
      i_x = 1
      do i = 2, size(x_grid)
         dist = abs(x_val - x_grid(i))
         if (dist < min_dist) then
            min_dist = dist
            i_x = i
         end if
      end do

      ! Find nearest y grid point
      min_dist = abs(y_val - y_grid(1))
      i_y = 1
      do i = 2, size(y_grid)
         dist = abs(y_val - y_grid(i))
         if (dist < min_dist) then
            min_dist = dist
            i_y = i
         end if
      end do

      ! Find nearest z grid point
      min_dist = abs(z_val - z_grid(1))
      i_z = 1
      do i = 2, size(z_grid)
         dist = abs(z_val - z_grid(i))
         if (dist < min_dist) then
            min_dist = dist
            i_z = i
         end if
      end do
   end subroutine find_nearest_point

end module linear_interp
