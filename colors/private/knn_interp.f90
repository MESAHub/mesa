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
! K-Nearest Neighbors interpolation module for spectral energy distributions (SEDs)
! ***********************************************************************

module knn_interp
   use const_def, only: dp
   use colors_utils, only: dilute_flux, load_sed
   implicit none

   private
   public :: construct_sed_knn, load_sed, interpolate_array, dilute_flux

contains

   !---------------------------------------------------------------------------
   ! Main entry point: Construct a SED using KNN interpolation
   !---------------------------------------------------------------------------
   subroutine construct_sed_knn(teff, log_g, metallicity, R, d, file_names, &
                                lu_teff, lu_logg, lu_meta, stellar_model_dir, &
                                wavelengths, fluxes)
      real(dp), intent(in) :: teff, log_g, metallicity, R, d
      real(dp), intent(in) :: lu_teff(:), lu_logg(:), lu_meta(:)
      character(len=*), intent(in) :: stellar_model_dir
      character(len=100), intent(in) :: file_names(:)
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, fluxes

      integer, dimension(4) :: closest_indices
      real(dp), dimension(:), allocatable :: temp_wavelengths, temp_flux, common_wavelengths
      real(dp), dimension(:, :), allocatable :: model_fluxes
      real(dp), dimension(4) :: weights, distances
      integer :: i, n_points
      real(dp) :: sum_weights
      real(dp), dimension(:), allocatable :: diluted_flux

      ! Get the four closest stellar models
      call get_closest_stellar_models(teff, log_g, metallicity, lu_teff, &
                                      lu_logg, lu_meta, closest_indices)

      ! Load the first SED to define the wavelength grid
      call load_sed(trim(stellar_model_dir)//trim(file_names(closest_indices(1))), &
                    closest_indices(1), temp_wavelengths, temp_flux)

      n_points = size(temp_wavelengths)
      allocate (common_wavelengths(n_points))
      common_wavelengths = temp_wavelengths

      ! Allocate flux array for the models (4 models, n_points each)
      allocate (model_fluxes(4, n_points))
      call interpolate_array(temp_wavelengths, temp_flux, common_wavelengths, model_fluxes(1, :))

      ! Load and interpolate remaining SEDs
      do i = 2, 4
         call load_sed(trim(stellar_model_dir)//trim(file_names(closest_indices(i))), &
                       closest_indices(i), temp_wavelengths, temp_flux)

         call interpolate_array(temp_wavelengths, temp_flux, common_wavelengths, model_fluxes(i, :))
      end do

      ! Compute distances and weights for the four models
      do i = 1, 4
         distances(i) = sqrt((lu_teff(closest_indices(i)) - teff)**2 + &
                             (lu_logg(closest_indices(i)) - log_g)**2 + &
                             (lu_meta(closest_indices(i)) - metallicity)**2)
         if (distances(i) == 0.0) distances(i) = 1.0d-10  ! Prevent division by zero
         weights(i) = 1.0/distances(i)
      end do

      ! Normalize weights
      sum_weights = sum(weights)
      weights = weights/sum_weights

      ! Allocate output arrays
      allocate (wavelengths(n_points), fluxes(n_points))
      wavelengths = common_wavelengths
      fluxes = 0.0

      ! Perform weighted combination of the model fluxes (still at the stellar surface)
      do i = 1, 4
         fluxes = fluxes + weights(i)*model_fluxes(i, :)
      end do

      ! Now, apply the dilution factor (R/d)^2 to convert the surface flux density
      ! into the observed flux density at Earth.
      allocate (diluted_flux(n_points))
      call dilute_flux(fluxes, R, d, diluted_flux)
      fluxes = diluted_flux

   end subroutine construct_sed_knn

   !---------------------------------------------------------------------------
   ! Identify the four closest stellar models
   !---------------------------------------------------------------------------
   subroutine get_closest_stellar_models(teff, log_g, metallicity, lu_teff, &
                                         lu_logg, lu_meta, closest_indices)
      real(dp), intent(in) :: teff, log_g, metallicity
      real(dp), intent(in) :: lu_teff(:), lu_logg(:), lu_meta(:)
      integer, dimension(4), intent(out) :: closest_indices
      logical :: use_teff_dim, use_logg_dim, use_meta_dim

      integer :: i, n, j
      real(dp) :: distance, norm_teff, norm_logg, norm_meta
      real(dp), dimension(:), allocatable :: scaled_lu_teff, scaled_lu_logg, scaled_lu_meta
      real(dp), dimension(4) :: min_distances
      integer, dimension(4) :: indices
      real(dp) :: teff_min, teff_max, logg_min, logg_max, meta_min, meta_max
      real(dp) :: teff_dist, logg_dist, meta_dist

      n = size(lu_teff)
      min_distances = huge(1.0)
      indices = -1

      ! Find min and max for normalization
      teff_min = minval(lu_teff)
      teff_max = maxval(lu_teff)
      logg_min = minval(lu_logg)
      logg_max = maxval(lu_logg)
      meta_min = minval(lu_meta)
      meta_max = maxval(lu_meta)

      ! Allocate and scale lookup table values
      allocate (scaled_lu_teff(n), scaled_lu_logg(n), scaled_lu_meta(n))

      if (teff_max - teff_min > 0.00) then
         scaled_lu_teff = (lu_teff - teff_min)/(teff_max - teff_min)
      end if

      if (logg_max - logg_min > 0.00) then
         scaled_lu_logg = (lu_logg - logg_min)/(logg_max - logg_min)
      end if

      if (meta_max - meta_min > 0.00) then
         scaled_lu_meta = (lu_meta - meta_min)/(meta_max - meta_min)
      end if

      ! Normalize input parameters
      norm_teff = (teff - teff_min)/(teff_max - teff_min)
      norm_logg = (log_g - logg_min)/(logg_max - logg_min)
      norm_meta = (metallicity - meta_min)/(meta_max - meta_min)

      ! Find closest models
      do i = 1, n
         teff_dist = 0.0
         logg_dist = 0.0
         meta_dist = 0.0

         if (teff_max - teff_min > 0.00) then
            teff_dist = scaled_lu_teff(i) - norm_teff
         end if

         if (logg_max - logg_min > 0.00) then
            logg_dist = scaled_lu_logg(i) - norm_logg
         end if

         if (meta_max - meta_min > 0.00) then
            meta_dist = scaled_lu_meta(i) - norm_meta
         end if

         ! Detect dummy axes once

         use_teff_dim = .not. (all(lu_teff == 0.0_dp) .or. all(lu_teff == 999.0_dp) .or. all(lu_teff == -999.0_dp))
         use_logg_dim = .not. (all(lu_logg == 0.0_dp) .or. all(lu_logg == 999.0_dp) .or. all(lu_logg == -999.0_dp))
         use_meta_dim = .not. (all(lu_meta == 0.0_dp) .or. all(lu_meta == 999.0_dp) .or. all(lu_meta == -999.0_dp))

         ! Inside the loop:
         distance = 0.0_dp
         if (use_teff_dim) distance = distance + teff_dist**2
         if (use_logg_dim) distance = distance + logg_dist**2
         if (use_meta_dim) distance = distance + meta_dist**2

         distance = teff_dist**2 + logg_dist**2 + meta_dist**2

         do j = 1, 4
            if (distance < min_distances(j)) then
               ! Shift larger distances down
               if (j < 4) then
                  min_distances(j + 1:4) = min_distances(j:3)
                  indices(j + 1:4) = indices(j:3)
               end if
               min_distances(j) = distance
               indices(j) = i
               exit
            end if
         end do
      end do

      closest_indices = indices
   end subroutine get_closest_stellar_models

   !---------------------------------------------------------------------------
   ! Linear interpolation (binary search version for efficiency)
   !---------------------------------------------------------------------------
   subroutine linear_interpolate(x, y, x_val, y_val)
      real(dp), intent(in) :: x(:), y(:), x_val
      real(dp), intent(out) :: y_val
      integer :: low, high, mid

      ! Validate input sizes
      if (size(x) < 2) then
         print *, "Error: x array has fewer than 2 points."
         y_val = 0.0_dp
         return
      end if

      if (size(x) /= size(y)) then
         print *, "Error: x and y arrays have different sizes."
         y_val = 0.0_dp
         return
      end if

      ! Handle out-of-bounds cases
      if (x_val <= x(1)) then
         y_val = y(1)
         return
      else if (x_val >= x(size(x))) then
         y_val = y(size(y))
         return
      end if

      ! Binary search to find the proper interval [x(low), x(low+1)]
      low = 1
      high = size(x)
      do while (high - low > 1)
         mid = (low + high)/2
         if (x(mid) <= x_val) then
            low = mid
         else
            high = mid
         end if
      end do

      ! Linear interpolation between x(low) and x(low+1)
      y_val = y(low) + (y(low + 1) - y(low))/(x(low + 1) - x(low))*(x_val - x(low))
   end subroutine linear_interpolate

   !---------------------------------------------------------------------------
   ! Array interpolation for SED construction
   !---------------------------------------------------------------------------
   subroutine interpolate_array(x_in, y_in, x_out, y_out)
      real(dp), intent(in) :: x_in(:), y_in(:), x_out(:)
      real(dp), intent(out) :: y_out(:)
      integer :: i

      ! Validate input sizes
      if (size(x_in) < 2 .or. size(y_in) < 2) then
         print *, "Error: x_in or y_in arrays have fewer than 2 points."
         stop
      end if

      if (size(x_in) /= size(y_in)) then
         print *, "Error: x_in and y_in arrays have different sizes."
         stop
      end if

      if (size(x_out) <= 0) then
         print *, "Error: x_out array is empty."
         stop
      end if

      do i = 1, size(x_out)
         call linear_interpolate(x_in, y_in, x_out(i), y_out(i))
      end do
   end subroutine interpolate_array

end module knn_interp
