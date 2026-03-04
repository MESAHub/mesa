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
!
! Supports two data-loading strategies, selected by rq%cube_loaded:
!   .true.  -> extract neighbor SEDs directly from the preloaded 4-D cube
!   .false. -> load individual SED files via the lookup table (fallback)
! ***********************************************************************

module knn_interp
   use const_def, only: dp
   use colors_def, only: Colors_General_Info
   use colors_utils, only: dilute_flux, load_sed_cached
   use utils_lib, only: mesa_error
   implicit none

   private
   public :: construct_sed_knn, interpolate_array

contains

   !---------------------------------------------------------------------------
   ! Main entry point: Construct a SED using KNN interpolation.
   ! Data loading strategy is determined by rq%cube_loaded (set at init):
   !   cube_loaded = .true.  -> use the preloaded 4-D cube on the handle
   !   cube_loaded = .false. -> load individual SED files via the lookup table
   !---------------------------------------------------------------------------
   subroutine construct_sed_knn(rq, teff, log_g, metallicity, R, d, &
                                stellar_model_dir, wavelengths, fluxes)
      type(Colors_General_Info), intent(inout) :: rq
      real(dp), intent(in) :: teff, log_g, metallicity, R, d
      character(len=*), intent(in) :: stellar_model_dir
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, fluxes

      integer :: n_lambda
      real(dp), dimension(:), allocatable :: interp_flux, diluted_flux

      if (rq%cube_loaded) then
         ! ---- Fast path: extract neighbors from preloaded cube ----
         call construct_sed_from_cube(rq, teff, log_g, metallicity, &
                                      interp_flux, wavelengths)
         n_lambda = size(wavelengths)
      else
         ! ---- Fallback path: load individual SED files ----
         call construct_sed_from_files(rq, teff, log_g, metallicity, &
                                       stellar_model_dir, interp_flux, wavelengths)
         n_lambda = size(wavelengths)
      end if

      ! Apply distance dilution to get observed flux
      allocate (diluted_flux(n_lambda))
      call dilute_flux(interp_flux, R, d, diluted_flux)
      fluxes = diluted_flux

   end subroutine construct_sed_knn

   !---------------------------------------------------------------------------
   ! Cube path: Find 4 nearest grid points in the structured cube grid,
   ! extract their SEDs directly from cube_flux, and blend by inverse
   ! distance weighting.  No file I/O required.
   !---------------------------------------------------------------------------
   subroutine construct_sed_from_cube(rq, teff, log_g, metallicity, &
                                       interp_flux, wavelengths)
      type(Colors_General_Info), intent(inout) :: rq
      real(dp), intent(in) :: teff, log_g, metallicity
      real(dp), dimension(:), allocatable, intent(out) :: interp_flux, wavelengths

      integer :: n_lambda, k
      integer, dimension(4) :: nbr_it, nbr_ig, nbr_im
      real(dp), dimension(4) :: distances, weights
      real(dp) :: sum_weights

      n_lambda = size(rq%cube_wavelengths)
      allocate (wavelengths(n_lambda))
      wavelengths = rq%cube_wavelengths

      ! Find the 4 nearest grid points in the structured cube
      call get_closest_grid_points(teff, log_g, metallicity, &
                                    rq%cube_teff_grid, rq%cube_logg_grid, &
                                    rq%cube_meta_grid, &
                                    nbr_it, nbr_ig, nbr_im, distances)

      ! Compute inverse-distance weights
      do k = 1, 4
         if (distances(k) == 0.0_dp) distances(k) = 1.0e-10_dp
         weights(k) = 1.0_dp / distances(k)
      end do
      sum_weights = sum(weights)
      weights = weights / sum_weights

      ! Blend neighbor SEDs from cube
      allocate (interp_flux(n_lambda))
      interp_flux = 0.0_dp
      do k = 1, 4
         interp_flux = interp_flux + weights(k) * &
            rq%cube_flux(nbr_it(k), nbr_ig(k), nbr_im(k), :)
      end do

   end subroutine construct_sed_from_cube

   !---------------------------------------------------------------------------
   ! Fallback path: Find 4 nearest models in the flat lookup table,
   ! load their SEDs via the memory cache, and blend by inverse
   ! distance weighting.
   !---------------------------------------------------------------------------
   subroutine construct_sed_from_files(rq, teff, log_g, metallicity, &
                                        stellar_model_dir, interp_flux, wavelengths)
      use colors_utils, only: resolve_path
      type(Colors_General_Info), intent(inout) :: rq
      real(dp), intent(in) :: teff, log_g, metallicity
      character(len=*), intent(in) :: stellar_model_dir
      real(dp), dimension(:), allocatable, intent(out) :: interp_flux, wavelengths

      integer, dimension(4) :: closest_indices
      real(dp), dimension(:), allocatable :: temp_wavelengths, temp_flux, common_wavelengths
      real(dp), dimension(:, :), allocatable :: model_fluxes
      real(dp), dimension(4) :: weights, distances
      integer :: i, n_points
      real(dp) :: sum_weights
      character(len=512) :: resolved_dir

      resolved_dir = trim(resolve_path(stellar_model_dir))

      ! Get the four closest stellar models from the flat lookup table
      call get_closest_stellar_models(teff, log_g, metallicity, &
                                      rq%lu_teff, rq%lu_logg, rq%lu_meta, &
                                      closest_indices)

      ! Load the first SED to define the wavelength grid (using cache)
      call load_sed_cached(rq, resolved_dir, closest_indices(1), temp_flux)

      ! Get wavelengths from canonical copy on the handle
      if (rq%fallback_wavelengths_set) then
         n_points = size(rq%fallback_wavelengths)
         allocate (common_wavelengths(n_points))
         common_wavelengths = rq%fallback_wavelengths
      else
         ! Should not happen — load_sed_cached sets this on first call
         print *, 'KNN fallback: wavelengths not set after first SED load'
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Allocate flux array for the models (4 models, n_points each)
      allocate (model_fluxes(4, n_points))
      model_fluxes(1, :) = temp_flux(1:n_points)
      if (allocated(temp_flux)) deallocate(temp_flux)

      ! Load and store remaining SEDs
      do i = 2, 4
         call load_sed_cached(rq, resolved_dir, closest_indices(i), temp_flux)
         model_fluxes(i, :) = temp_flux(1:n_points)
         if (allocated(temp_flux)) deallocate(temp_flux)
      end do

      ! Compute distances and weights for the four models
      do i = 1, 4
         distances(i) = sqrt((rq%lu_teff(closest_indices(i)) - teff)**2 + &
                             (rq%lu_logg(closest_indices(i)) - log_g)**2 + &
                             (rq%lu_meta(closest_indices(i)) - metallicity)**2)
         if (distances(i) == 0.0_dp) distances(i) = 1.0e-10_dp
         weights(i) = 1.0_dp / distances(i)
      end do

      ! Normalize weights
      sum_weights = sum(weights)
      weights = weights / sum_weights

      ! Allocate output arrays
      allocate (wavelengths(n_points))
      wavelengths = common_wavelengths

      allocate (interp_flux(n_points))
      interp_flux = 0.0_dp

      ! Perform weighted combination of the model fluxes
      do i = 1, 4
         interp_flux = interp_flux + weights(i) * model_fluxes(i, :)
      end do

   end subroutine construct_sed_from_files

   !---------------------------------------------------------------------------
   ! Find the 4 closest grid points in the structured cube grid.
   ! Searches over all (i_t, i_g, i_m) combinations using normalised
   ! Euclidean distance (same scaling logic as get_closest_stellar_models).
   !---------------------------------------------------------------------------
   subroutine get_closest_grid_points(teff, log_g, metallicity, &
                                       teff_grid, logg_grid, meta_grid, &
                                       nbr_it, nbr_ig, nbr_im, distances)
      real(dp), intent(in) :: teff, log_g, metallicity
      real(dp), intent(in) :: teff_grid(:), logg_grid(:), meta_grid(:)
      integer, dimension(4), intent(out) :: nbr_it, nbr_ig, nbr_im
      real(dp), dimension(4), intent(out) :: distances

      integer :: it, ig, im, j
      real(dp) :: dist, norm_teff, norm_logg, norm_meta
      real(dp) :: teff_min, teff_max, logg_min, logg_max, meta_min, meta_max
      real(dp) :: scaled_t, scaled_g, scaled_m, dt, dg, dm
      logical :: use_teff_dim, use_logg_dim, use_meta_dim

      distances = huge(1.0_dp)
      nbr_it = 1; nbr_ig = 1; nbr_im = 1

      ! Normalisation ranges
      teff_min = minval(teff_grid); teff_max = maxval(teff_grid)
      logg_min = minval(logg_grid); logg_max = maxval(logg_grid)
      meta_min = minval(meta_grid); meta_max = maxval(meta_grid)

      ! Detect dummy axes
      use_teff_dim = .not. (all(teff_grid == 0.0_dp) .or. &
                            all(teff_grid == 999.0_dp) .or. all(teff_grid == -999.0_dp))
      use_logg_dim = .not. (all(logg_grid == 0.0_dp) .or. &
                            all(logg_grid == 999.0_dp) .or. all(logg_grid == -999.0_dp))
      use_meta_dim = .not. (all(meta_grid == 0.0_dp) .or. &
                            all(meta_grid == 999.0_dp) .or. all(meta_grid == -999.0_dp))

      ! Normalised target values
      norm_teff = 0.0_dp; norm_logg = 0.0_dp; norm_meta = 0.0_dp
      if (use_teff_dim .and. teff_max - teff_min > 0.0_dp) &
         norm_teff = (teff - teff_min) / (teff_max - teff_min)
      if (use_logg_dim .and. logg_max - logg_min > 0.0_dp) &
         norm_logg = (log_g - logg_min) / (logg_max - logg_min)
      if (use_meta_dim .and. meta_max - meta_min > 0.0_dp) &
         norm_meta = (metallicity - meta_min) / (meta_max - meta_min)

      do it = 1, size(teff_grid)
         if (use_teff_dim .and. teff_max - teff_min > 0.0_dp) then
            scaled_t = (teff_grid(it) - teff_min) / (teff_max - teff_min)
         else
            scaled_t = 0.0_dp
         end if
         dt = 0.0_dp
         if (use_teff_dim) dt = (scaled_t - norm_teff)**2

         do ig = 1, size(logg_grid)
            if (use_logg_dim .and. logg_max - logg_min > 0.0_dp) then
               scaled_g = (logg_grid(ig) - logg_min) / (logg_max - logg_min)
            else
               scaled_g = 0.0_dp
            end if
            dg = 0.0_dp
            if (use_logg_dim) dg = (scaled_g - norm_logg)**2

            do im = 1, size(meta_grid)
               if (use_meta_dim .and. meta_max - meta_min > 0.0_dp) then
                  scaled_m = (meta_grid(im) - meta_min) / (meta_max - meta_min)
               else
                  scaled_m = 0.0_dp
               end if
               dm = 0.0_dp
               if (use_meta_dim) dm = (scaled_m - norm_meta)**2

               dist = dt + dg + dm

               ! Insert into sorted top-4 if closer
               do j = 1, 4
                  if (dist < distances(j)) then
                     if (j < 4) then
                        distances(j + 1:4) = distances(j:3)
                        nbr_it(j + 1:4) = nbr_it(j:3)
                        nbr_ig(j + 1:4) = nbr_ig(j:3)
                        nbr_im(j + 1:4) = nbr_im(j:3)
                     end if
                     distances(j) = dist
                     nbr_it(j) = it
                     nbr_ig(j) = ig
                     nbr_im(j) = im
                     exit
                  end if
               end do
            end do
         end do
      end do

      ! Convert squared distances to actual distances for weighting
      do j = 1, 4
         distances(j) = sqrt(distances(j))
      end do

   end subroutine get_closest_grid_points

   !---------------------------------------------------------------------------
   ! Identify the four closest stellar models in the flat lookup table
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

      if (teff_max - teff_min > 0.0_dp) then
         scaled_lu_teff = (lu_teff - teff_min)/(teff_max - teff_min)
      end if

      if (logg_max - logg_min > 0.0_dp) then
         scaled_lu_logg = (lu_logg - logg_min)/(logg_max - logg_min)
      end if

      if (meta_max - meta_min > 0.0_dp) then
         scaled_lu_meta = (lu_meta - meta_min)/(meta_max - meta_min)
      end if

      ! Normalize input parameters
      norm_teff = (teff - teff_min)/(teff_max - teff_min)
      norm_logg = (log_g - logg_min)/(logg_max - logg_min)
      norm_meta = (metallicity - meta_min)/(meta_max - meta_min)

      ! Detect dummy axes once (outside the loop)
      use_teff_dim = .not. (all(lu_teff == 0.0_dp) .or. all(lu_teff == 999.0_dp) .or. all(lu_teff == -999.0_dp))
      use_logg_dim = .not. (all(lu_logg == 0.0_dp) .or. all(lu_logg == 999.0_dp) .or. all(lu_logg == -999.0_dp))
      use_meta_dim = .not. (all(lu_meta == 0.0_dp) .or. all(lu_meta == 999.0_dp) .or. all(lu_meta == -999.0_dp))

      ! Find closest models
      do i = 1, n
         teff_dist = 0.0_dp
         logg_dist = 0.0_dp
         meta_dist = 0.0_dp

         if (teff_max - teff_min > 0.0_dp) then
            teff_dist = scaled_lu_teff(i) - norm_teff
         end if

         if (logg_max - logg_min > 0.0_dp) then
            logg_dist = scaled_lu_logg(i) - norm_logg
         end if

         if (meta_max - meta_min > 0.0_dp) then
            meta_dist = scaled_lu_meta(i) - norm_meta
         end if

         ! Compute distance using only valid dimensions
         distance = 0.0_dp
         if (use_teff_dim) distance = distance + teff_dist**2
         if (use_logg_dim) distance = distance + logg_dist**2
         if (use_meta_dim) distance = distance + meta_dist**2

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
         call mesa_error(__FILE__, __LINE__)
      end if

      if (size(x_in) /= size(y_in)) then
         print *, "Error: x_in and y_in arrays have different sizes."
         call mesa_error(__FILE__, __LINE__)
      end if

      if (size(x_out) <= 0) then
         print *, "Error: x_out array is empty."
         call mesa_error(__FILE__, __LINE__)
      end if

      do i = 1, size(x_out)
         call linear_interpolate(x_in, y_in, x_out(i), y_out(i))
      end do
   end subroutine interpolate_array

end module knn_interp