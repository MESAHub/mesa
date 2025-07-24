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

module hermite_interp
   use const_def, only: dp
   use colors_utils, only: dilute_flux
   implicit none

   private
   public :: construct_sed_hermite, hermite_tensor_interp3d

contains

   !---------------------------------------------------------------------------
   ! Main entry point: Construct a SED using Hermite tensor interpolation
   !---------------------------------------------------------------------------
   subroutine construct_sed_hermite(teff, log_g, metallicity, R, d, file_names, &
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

      ! Parameter grids
      real(dp), allocatable :: teff_grid(:), logg_grid(:), meta_grid(:)
      character(len=256) :: bin_filename

      ! Construct the binary filename
      bin_filename = trim(stellar_model_dir)//'/flux_cube.bin'

      ! Load the data from binary file
      call load_binary_data(bin_filename, teff_grid, logg_grid, meta_grid, &
                            wavelengths, precomputed_flux_cube, status)

      n_teff = size(teff_grid)
      n_logg = size(logg_grid)
      n_meta = size(meta_grid)
      n_lambda = size(wavelengths)

      ! Allocate space for interpolated flux
      allocate (interp_flux(n_lambda))

      ! Process each wavelength point
      do i = 1, n_lambda
         allocate (flux_cube_lambda(n_teff, n_logg, n_meta))
         flux_cube_lambda = precomputed_flux_cube(:, :, :, i)

         interp_flux(i) = hermite_tensor_interp3d(teff, log_g, metallicity, &
                                                  teff_grid, logg_grid, meta_grid, flux_cube_lambda)

         deallocate (flux_cube_lambda)
      end do

      ! Apply distance dilution to get observed flux
      allocate (diluted_flux(n_lambda))
      call dilute_flux(interp_flux, R, d, diluted_flux)
      fluxes = diluted_flux

   end subroutine construct_sed_hermite

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
         print *, 'Error opening binary file:', trim(filename)
         return
      end if

      ! Read dimensions
      read (unit, iostat=status) n_teff, n_logg, n_meta, n_lambda
      if (status /= 0) then
         print *, 'Error reading dimensions from binary file'
         close (unit)
         return
      end if

      ! Allocate arrays based on dimensions
      allocate (teff_grid(n_teff), STAT=status)
      if (status /= 0) then
         print *, 'Error allocating teff_grid array'
         close (unit)
         return
      end if

      allocate (logg_grid(n_logg), STAT=status)
      if (status /= 0) then
         print *, 'Error allocating logg_grid array'
         close (unit)
         return
      end if

      allocate (meta_grid(n_meta), STAT=status)
      if (status /= 0) then
         print *, 'Error allocating meta_grid array'
         close (unit)
         return
      end if

      allocate (wavelengths(n_lambda), STAT=status)
      if (status /= 0) then
         print *, 'Error allocating wavelengths array'
         close (unit)
         return
      end if

      allocate (flux_cube(n_teff, n_logg, n_meta, n_lambda), STAT=status)
      if (status /= 0) then
         print *, 'Error allocating flux_cube array'
         close (unit)
         return
      end if

      ! Read grid arrays
      read (unit, iostat=status) teff_grid
      if (status /= 0) then
         print *, 'Error reading teff_grid'
         close (unit)
         return
      end if

      read (unit, iostat=status) logg_grid
      if (status /= 0) then
         print *, 'Error reading logg_grid'
         close (unit)
         return
      end if

      read (unit, iostat=status) meta_grid
      if (status /= 0) then
         print *, 'Error reading meta_grid'
         close (unit)
         return
      end if

      read (unit, iostat=status) wavelengths
      if (status /= 0) then
         print *, 'Error reading wavelengths'
         close (unit)
         return
      end if

      ! Read flux cube
      read (unit, iostat=status) flux_cube
      if (status /= 0) then
         print *, 'Error reading flux_cube'
         close (unit)
         return
      end if

      ! Close file and return success
      close (unit)
   end subroutine load_binary_data

!---------------------------------------------------------------------------
   ! Main 3D Hermite interpolation function
   !---------------------------------------------------------------------------
   function hermite_tensor_interp3d(x_val, y_val, z_val, x_grid, y_grid, &
                                    z_grid, f_values) result(f_interp)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in) :: f_values(:, :, :)
      real(dp) :: f_interp

      integer :: i_x, i_y, i_z
      real(dp) :: t_x, t_y, t_z
      real(dp) :: dx, dy, dz
      real(dp) :: dx_values(2, 2, 2), dy_values(2, 2, 2), dz_values(2, 2, 2)
      real(dp) :: values(2, 2, 2)
      real(dp) :: sum
      integer :: ix, iy, iz
      real(dp) :: h_x(2), h_y(2), h_z(2)

      ! Find containing cell and parameter values
      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      ! If outside grid, use nearest point
      if (i_x < 1 .OR. i_x >= size(x_grid) .OR. &
          i_y < 1 .OR. i_y >= size(y_grid) .OR. &
          i_z < 1 .OR. i_z >= size(z_grid)) then

         call find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
         f_interp = f_values(i_x, i_y, i_z)
         return
      end if

      ! Grid cell spacing
      dx = x_grid(i_x + 1) - x_grid(i_x)
      dy = y_grid(i_y + 1) - y_grid(i_y)
      dz = z_grid(i_z + 1) - z_grid(i_z)

      ! Extract the local 2x2x2 grid cell and compute derivatives
      do iz = 0, 1
         do iy = 0, 1
            do ix = 0, 1
               values(ix + 1, iy + 1, iz + 1) = f_values(i_x + ix, i_y + iy, i_z + iz)
               call compute_derivatives_at_point(f_values, i_x + ix, i_y + iy, i_z + iz, &
                                                 size(x_grid), size(y_grid), size(z_grid), &
                                                 dx, dy, dz, &
                                                 dx_values(ix + 1, iy + 1, iz + 1), &
                                                 dy_values(ix + 1, iy + 1, iz + 1), &
                                                 dz_values(ix + 1, iy + 1, iz + 1))
            end do
         end do
      end do

      sum = 0.0_dp

      ! Function values
      h_x = [h00(t_x), h01(t_x)]
      h_y = [h00(t_y), h01(t_y)]
      h_z = [h00(t_z), h01(t_z)]
      do iz = 1, 2
         do iy = 1, 2
            do ix = 1, 2
               sum = sum + h_x(ix)*h_y(iy)*h_z(iz)*values(ix, iy, iz)
            end do
         end do
      end do

      ! x-derivatives
      h_x = [h10(t_x), h11(t_x)]
      do iz = 1, 2
         do iy = 1, 2
            do ix = 1, 2
               sum = sum + h_x(ix)*h_y(iy)*h_z(iz)*dx*dx_values(ix, iy, iz)
            end do
         end do
      end do

      ! y-derivatives
      h_x = [h00(t_x), h01(t_x)]
      h_y = [h10(t_y), h11(t_y)]
      do iz = 1, 2
         do iy = 1, 2
            do ix = 1, 2
               sum = sum + h_x(ix)*h_y(iy)*h_z(iz)*dy*dy_values(ix, iy, iz)
            end do
         end do
      end do

      ! z-derivatives
      h_y = [h00(t_y), h01(t_y)]
      h_z = [h10(t_z), h11(t_z)]
      do iz = 1, 2
         do iy = 1, 2
            do ix = 1, 2
               sum = sum + h_x(ix)*h_y(iy)*h_z(iz)*dz*dz_values(ix, iy, iz)
            end do
         end do
      end do

      f_interp = sum
   end function hermite_tensor_interp3d

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
   SUBROUTINE find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
     REAL(dp), INTENT(IN) :: x_val, y_val, z_val
     REAL(dp), INTENT(IN) :: x_grid(:), y_grid(:), z_grid(:)
     INTEGER, INTENT(OUT) :: i_x, i_y, i_z

     ! Find nearest grid points using intrinsic MINLOC
     i_x = MINLOC(ABS(x_val - x_grid), 1)
     i_y = MINLOC(ABS(y_val - y_grid), 1)
     i_z = MINLOC(ABS(z_val - z_grid), 1)
   END SUBROUTINE find_nearest_point

   !---------------------------------------------------------------------------
   ! Compute derivatives at a grid point
   !---------------------------------------------------------------------------
   subroutine compute_derivatives_at_point(f, i, j, k, nx, ny, nz, dx, dy, dz, &
                                           df_dx, df_dy, df_dz)
      real(dp), intent(in) :: f(:, :, :)
      integer, intent(in) :: i, j, k, nx, ny, nz
      real(dp), intent(in) :: dx, dy, dz
      real(dp), intent(out) :: df_dx, df_dy, df_dz

      ! Compute x derivative using centered differences where possible
      if (i > 1 .AND. i < nx) then
         df_dx = (f(i + 1, j, k) - f(i - 1, j, k))/(2.0_dp*dx)
      else if (i == 1) then
         df_dx = (f(i + 1, j, k) - f(i, j, k))/dx
      else ! i == nx
         df_dx = (f(i, j, k) - f(i - 1, j, k))/dx
      end if

      ! Compute y derivative using centered differences where possible
      if (j > 1 .AND. j < ny) then
         df_dy = (f(i, j + 1, k) - f(i, j - 1, k))/(2.0_dp*dy)
      else if (j == 1) then
         df_dy = (f(i, j + 1, k) - f(i, j, k))/dy
      else ! j == ny
         df_dy = (f(i, j, k) - f(i, j - 1, k))/dy
      end if

      ! Compute z derivative using centered differences where possible
      if (k > 1 .AND. k < nz) then
         df_dz = (f(i, j, k + 1) - f(i, j, k - 1))/(2.0_dp*dz)
      else if (k == 1) then
         df_dz = (f(i, j, k + 1) - f(i, j, k))/dz
      else ! k == nz
         df_dz = (f(i, j, k) - f(i, j, k - 1))/dz
      end if
   end subroutine compute_derivatives_at_point

   !---------------------------------------------------------------------------
   ! Hermite basis functions
   !---------------------------------------------------------------------------
   function h00(t) result(h)
      real(dp), intent(in) :: t
      real(dp) :: h
      h = 2.0_dp*t**3 - 3.0_dp*t**2 + 1.0_dp
   end function h00

   function h10(t) result(h)
      real(dp), intent(in) :: t
      real(dp) :: h
      h = t**3 - 2.0_dp*t**2 + t
   end function h10

   function h01(t) result(h)
      real(dp), intent(in) :: t
      real(dp) :: h
      h = -2.0_dp*t**3 + 3.0_dp*t**2
   end function h01

   function h11(t) result(h)
      real(dp), intent(in) :: t
      real(dp) :: h
      h = t**3 - t**2
   end function h11

end module hermite_interp
