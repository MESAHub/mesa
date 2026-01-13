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

module bolometric

   use const_def, only: dp
   use colors_utils, only: romberg_integration
   use hermite_interp, only: construct_sed_hermite
   use linear_interp, only: construct_sed_linear
   use knn_interp, only: construct_sed_knn

   implicit none

   private
   public :: calculate_bolometric

contains

   !****************************
   ! Calculate Bolometric Photometry Using Multiple SEDs
   ! Accepts cached lookup table data instead of loading from file
   !****************************
   subroutine calculate_bolometric(teff, log_g, metallicity, R, d, bolometric_magnitude, &
                                   bolometric_flux, wavelengths, fluxes, sed_filepath, interpolation_radius, &
                                   lu_file_names, lu_teff, lu_logg, lu_meta)
      real(dp), intent(in) :: teff, log_g, metallicity, R, d
      character(len=*), intent(in) :: sed_filepath
      real(dp), intent(out) :: bolometric_magnitude, bolometric_flux, interpolation_radius
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, fluxes

      ! Cached lookup table data (passed in from colors_settings)
      character(len=100), intent(in) :: lu_file_names(:)
      real(dp), intent(in) :: lu_teff(:), lu_logg(:), lu_meta(:)

      character(len=32) :: interpolation_method

      interpolation_method = 'Hermite'   ! or 'Linear' / 'KNN' later

      ! Quantify how far (teff, log_g, metallicity) is from the grid points
      interpolation_radius = compute_interp_radius(teff, log_g, metallicity, &
                                                   lu_teff, lu_logg, lu_meta)

      select case (interpolation_method)
      case ('Hermite', 'hermite', 'HERMITE')
         call construct_sed_hermite(teff, log_g, metallicity, R, d, lu_file_names, &
                                    lu_teff, lu_logg, lu_meta, sed_filepath, &
                                    wavelengths, fluxes)

      case ('Linear', 'linear', 'LINEAR')
         call construct_sed_linear(teff, log_g, metallicity, R, d, lu_file_names, &
                                   lu_teff, lu_logg, lu_meta, sed_filepath, &
                                   wavelengths, fluxes)

      case ('KNN', 'knn', 'Knn')
         call construct_sed_knn(teff, log_g, metallicity, R, d, lu_file_names, &
                                lu_teff, lu_logg, lu_meta, sed_filepath, &
                                wavelengths, fluxes)

      case default
         ! Fallback: Hermite
         call construct_sed_hermite(teff, log_g, metallicity, R, d, lu_file_names, &
                                    lu_teff, lu_logg, lu_meta, sed_filepath, &
                                    wavelengths, fluxes)
      end select

      ! Calculate bolometric flux and magnitude
      call calculate_bolometric_phot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
   end subroutine calculate_bolometric

   !****************************
   ! Calculate Bolometric Magnitude and Flux
   !****************************
   subroutine calculate_bolometric_phot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
      real(dp), dimension(:), intent(inout) :: wavelengths, fluxes
      real(dp), intent(out) :: bolometric_magnitude, bolometric_flux
      integer :: i

      ! Validate inputs and replace invalid values with 0
      do i = 1, size(wavelengths) - 1
         if (wavelengths(i) <= 0.0d0 .or. fluxes(i) < 0.0d0) then
            fluxes(i) = 0.0d0
         end if
      end do

      ! Integrate to get bolometric flux
      call romberg_integration(wavelengths, fluxes, bolometric_flux)

      ! Validate and calculate magnitude
      if (bolometric_flux <= 0.0d0) then
         print *, "Error: Flux integration resulted in non-positive value."
         bolometric_magnitude = 99.0d0
         return
      else if (bolometric_flux < 1.0d-10) then
         print *, "Warning: Flux value is very small, precision might be affected."
      end if

      bolometric_magnitude = flux_to_magnitude(bolometric_flux)
   end subroutine calculate_bolometric_phot

   !****************************
   ! Convert Flux to Magnitude
   !****************************
   real(dp) function flux_to_magnitude(flux)
      real(dp), intent(in) :: flux
      if (flux <= 0.0d0) then
         print *, "Error: Flux must be positive to calculate magnitude."
         flux_to_magnitude = 99.0d0
      else
         flux_to_magnitude = -2.5d0 * log10(flux)
      end if
   end function flux_to_magnitude

   !--------------------------------------------------------------------
   ! Scalar metric: distance to nearest grid point in normalized space
   !--------------------------------------------------------------------
   real(dp) function compute_interp_radius(teff, log_g, metallicity, &
                                           lu_teff, lu_logg, lu_meta)

      real(dp), intent(in) :: teff, log_g, metallicity
      real(dp), intent(in) :: lu_teff(:), lu_logg(:), lu_meta(:)

      real(dp) :: teff_min, teff_max, logg_min, logg_max, meta_min, meta_max
      real(dp) :: teff_range, logg_range, meta_range
      real(dp) :: norm_teff, norm_logg, norm_meta
      real(dp) :: grid_teff, grid_logg, grid_meta
      real(dp) :: d, d_min
      integer  :: i, n

      logical :: use_teff, use_logg, use_meta
      real(dp), parameter :: eps = 1.0d-12

      ! Detect dummy columns (entire axis is 0 or ±999)
      use_teff = .not. (all(lu_teff == 0.0d0) .or. &
                        all(lu_teff == 999.0d0) .or. &
                        all(lu_teff == -999.0d0))

      use_logg = .not. (all(lu_logg == 0.0d0) .or. &
                        all(lu_logg == 999.0d0) .or. &
                        all(lu_logg == -999.0d0))

      use_meta = .not. (all(lu_meta == 0.0d0) .or. &
                        all(lu_meta == 999.0d0) .or. &
                        all(lu_meta == -999.0d0))

      ! Compute min/max for valid axes
      if (use_teff) then
         teff_min = minval(lu_teff)
         teff_max = maxval(lu_teff)
         teff_range = max(teff_max - teff_min, eps)
         norm_teff = (teff - teff_min) / teff_range
      end if

      if (use_logg) then
         logg_min = minval(lu_logg)
         logg_max = maxval(lu_logg)
         logg_range = max(logg_max - logg_min, eps)
         norm_logg = (log_g - logg_min) / logg_range
      end if

      if (use_meta) then
         meta_min = minval(lu_meta)
         meta_max = maxval(lu_meta)
         meta_range = max(meta_max - meta_min, eps)
         norm_meta = (metallicity - meta_min) / meta_range
      end if

      ! Find minimum distance to any grid point
      d_min = huge(1.0d0)
      n = size(lu_teff)

      do i = 1, n
         d = 0.0d0

         if (use_teff) then
            grid_teff = (lu_teff(i) - teff_min) / teff_range
            d = d + (norm_teff - grid_teff)**2
         end if

         if (use_logg) then
            grid_logg = (lu_logg(i) - logg_min) / logg_range
            d = d + (norm_logg - grid_logg)**2
         end if

         if (use_meta) then
            grid_meta = (lu_meta(i) - meta_min) / meta_range
            d = d + (norm_meta - grid_meta)**2
         end if

         d = sqrt(d)
         if (d < d_min) d_min = d
      end do

      compute_interp_radius = d_min

   end function compute_interp_radius

end module bolometric