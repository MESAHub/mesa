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
!
! Supports two data-loading strategies, selected by rq%cube_loaded:
!   .true.  -> use the preloaded 4-D flux cube on the handle
!   .false. -> load individual SED files via the lookup table (fallback)
! ***********************************************************************

module linear_interp
   use const_def, only: dp
   use colors_def, only: Colors_General_Info
   use colors_utils, only: dilute_flux, find_containing_cell, find_interval, &
                          find_nearest_point, find_bracket_index, &
                          load_sed_cached, load_stencil
   use utils_lib, only: mesa_error
   implicit none

   private
   public :: construct_sed_linear, trilinear_interp

contains

   !---------------------------------------------------------------------------
   ! Main entry point: Construct a SED using trilinear interpolation.
   ! Data loading strategy is determined by rq%cube_loaded (set at init):
   !   cube_loaded = .true.  -> use the preloaded 4-D cube on the handle
   !   cube_loaded = .false. -> load individual SED files via the lookup table
   !---------------------------------------------------------------------------
   subroutine construct_sed_linear(rq, teff, log_g, metallicity, R, d, &
                                   stellar_model_dir, wavelengths, fluxes)
      type(Colors_General_Info), intent(inout) :: rq
      real(dp), intent(in) :: teff, log_g, metallicity, R, d
      character(len=*), intent(in) :: stellar_model_dir
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, fluxes

      integer :: n_lambda
      real(dp), dimension(:), allocatable :: interp_flux, diluted_flux

      if (rq%cube_loaded) then
         ! ---- Fast path: use preloaded cube from handle ----
         n_lambda = size(rq%cube_wavelengths)

         ! Copy wavelengths to output
         allocate (wavelengths(n_lambda))
         wavelengths = rq%cube_wavelengths

         ! Vectorised interpolation over all wavelengths in one pass —
         ! cell location is computed once and reused.
         allocate (interp_flux(n_lambda))
         call trilinear_interp_vector(teff, log_g, metallicity, &
                                       rq%cube_teff_grid, rq%cube_logg_grid, &
                                       rq%cube_meta_grid, &
                                       rq%cube_flux, n_lambda, interp_flux)
      else
         ! ---- Fallback path: load individual SED files from lookup table ----
         call construct_sed_from_files(rq, teff, log_g, metallicity, &
                                       stellar_model_dir, interp_flux, wavelengths)
         n_lambda = size(wavelengths)
      end if

      ! Apply distance dilution to get observed flux
      allocate (diluted_flux(n_lambda))
      call dilute_flux(interp_flux, R, d, diluted_flux)
      fluxes = diluted_flux

   end subroutine construct_sed_linear

   !---------------------------------------------------------------------------
   ! Fallback: Build a local 2x2x2 sub-cube from individual SED files,
   ! then trilinear-interpolate all wavelengths in a single pass.
   !
   ! Unlike the Hermite method, trilinear interpolation does not require
   ! derivative context, so the stencil is exactly the 2x2x2 cell corners
   ! (no extension needed).
   !---------------------------------------------------------------------------
   subroutine construct_sed_from_files(rq, teff, log_g, metallicity, &
                                        stellar_model_dir, interp_flux, wavelengths)
      use colors_utils, only: resolve_path, build_grid_to_lu_map
      type(Colors_General_Info), intent(inout) :: rq
      real(dp), intent(in) :: teff, log_g, metallicity
      character(len=*), intent(in) :: stellar_model_dir
      real(dp), dimension(:), allocatable, intent(out) :: interp_flux, wavelengths

      integer :: i_t, i_g, i_m   ! bracketing indices in unique grids
      integer :: lo_t, hi_t, lo_g, hi_g, lo_m, hi_m   ! stencil bounds
      integer :: nt, ng, nm, n_lambda
      character(len=512) :: resolved_dir
      logical :: need_reload

      resolved_dir = trim(resolve_path(stellar_model_dir))

      ! Ensure the grid-to-lu mapping exists (built once, then reused)
      if (.not. rq%grid_map_built) call build_grid_to_lu_map(rq)

      ! Find bracketing cell in the unique grids
      call find_bracket_index(rq%u_teff, teff, i_t)
      call find_bracket_index(rq%u_logg, log_g, i_g)
      call find_bracket_index(rq%u_meta, metallicity, i_m)

      ! Check if the stencil cache is still valid for this cell
      need_reload = .true.
      if (rq%stencil_valid .and. &
          i_t == rq%stencil_i_t .and. &
          i_g == rq%stencil_i_g .and. &
          i_m == rq%stencil_i_m) then
         need_reload = .false.
      end if

      if (need_reload) then
         ! Trilinear needs exactly the 2x2x2 cell corners — no extension.
         nt = size(rq%u_teff)
         ng = size(rq%u_logg)
         nm = size(rq%u_meta)

         if (nt < 2) then
            lo_t = 1; hi_t = 1
         else
            lo_t = i_t
            hi_t = min(nt, i_t + 1)
         end if

         if (ng < 2) then
            lo_g = 1; hi_g = 1
         else
            lo_g = i_g
            hi_g = min(ng, i_g + 1)
         end if

         if (nm < 2) then
            lo_m = 1; hi_m = 1
         else
            lo_m = i_m
            hi_m = min(nm, i_m + 1)
         end if

         ! Load SEDs for every stencil point (using memory cache)
         call load_stencil(rq, resolved_dir, lo_t, hi_t, lo_g, hi_g, lo_m, hi_m)

         ! Store subgrid arrays on the handle
         if (allocated(rq%stencil_teff)) deallocate(rq%stencil_teff)
         if (allocated(rq%stencil_logg)) deallocate(rq%stencil_logg)
         if (allocated(rq%stencil_meta)) deallocate(rq%stencil_meta)

         allocate (rq%stencil_teff(hi_t - lo_t + 1))
         allocate (rq%stencil_logg(hi_g - lo_g + 1))
         allocate (rq%stencil_meta(hi_m - lo_m + 1))
         rq%stencil_teff = rq%u_teff(lo_t:hi_t)
         rq%stencil_logg = rq%u_logg(lo_g:hi_g)
         rq%stencil_meta = rq%u_meta(lo_m:hi_m)

         rq%stencil_i_t = i_t
         rq%stencil_i_g = i_g
         rq%stencil_i_m = i_m
         rq%stencil_valid = .true.
      end if

      ! Copy wavelengths to output
      n_lambda = size(rq%stencil_wavelengths)
      allocate (wavelengths(n_lambda))
      wavelengths = rq%stencil_wavelengths

      ! Interpolate all wavelengths using precomputed stencil
      allocate (interp_flux(n_lambda))
      call trilinear_interp_vector(teff, log_g, metallicity, &
                                    rq%stencil_teff, rq%stencil_logg, rq%stencil_meta, &
                                    rq%stencil_fluxes, n_lambda, interp_flux)

   end subroutine construct_sed_from_files

   !---------------------------------------------------------------------------
   ! Vectorised trilinear interpolation over all wavelengths.
   !
   ! The cell location (i_x, i_y, i_z, t_x, t_y, t_z) depends only on
   ! (teff, logg, meta) and the grids — not on wavelength.  Computing
   ! it once and reusing across all n_lambda samples eliminates redundant
   ! binary searches.
   !---------------------------------------------------------------------------
   subroutine trilinear_interp_vector(x_val, y_val, z_val, &
                                       x_grid, y_grid, z_grid, &
                                       f_values_4d, n_lambda, result_flux)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in) :: f_values_4d(:,:,:,:)   ! (nx, ny, nz, n_lambda)
      integer, intent(in) :: n_lambda
      real(dp), intent(out) :: result_flux(n_lambda)

      integer :: i_x, i_y, i_z, lam
      real(dp) :: t_x, t_y, t_z
      integer :: nx, ny, nz
      real(dp) :: c000, c001, c010, c011, c100, c101, c110, c111
      real(dp) :: c00, c01, c10, c11, c0, c1
      real(dp) :: lin_result, log_result
      real(dp), parameter :: tiny_value = 1.0e-10_dp

      nx = size(x_grid)
      ny = size(y_grid)
      nz = size(z_grid)

      ! Locate the cell once
      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      ! Boundary safety check
      if (i_x < 1) i_x = 1
      if (i_y < 1) i_y = 1
      if (i_z < 1) i_z = 1
      if (i_x >= nx) i_x = max(1, nx - 1)
      if (i_y >= ny) i_y = max(1, ny - 1)
      if (i_z >= nz) i_z = max(1, nz - 1)

      ! Clamp interpolation parameters to [0,1]
      t_x = max(0.0_dp, min(1.0_dp, t_x))
      t_y = max(0.0_dp, min(1.0_dp, t_y))
      t_z = max(0.0_dp, min(1.0_dp, t_z))

      ! Loop over wavelengths with the same cell location
      do lam = 1, n_lambda
         ! Get the 8 corners of the cube with safety floor
         c000 = max(tiny_value, f_values_4d(i_x,     i_y,     i_z,     lam))
         c001 = max(tiny_value, f_values_4d(i_x,     i_y,     i_z + 1, lam))
         c010 = max(tiny_value, f_values_4d(i_x,     i_y + 1, i_z,     lam))
         c011 = max(tiny_value, f_values_4d(i_x,     i_y + 1, i_z + 1, lam))
         c100 = max(tiny_value, f_values_4d(i_x + 1, i_y,     i_z,     lam))
         c101 = max(tiny_value, f_values_4d(i_x + 1, i_y,     i_z + 1, lam))
         c110 = max(tiny_value, f_values_4d(i_x + 1, i_y + 1, i_z,     lam))
         c111 = max(tiny_value, f_values_4d(i_x + 1, i_y + 1, i_z + 1, lam))

         ! Standard linear interpolation first (safer)
         c00 = c000*(1.0_dp - t_x) + c100*t_x
         c01 = c001*(1.0_dp - t_x) + c101*t_x
         c10 = c010*(1.0_dp - t_x) + c110*t_x
         c11 = c011*(1.0_dp - t_x) + c111*t_x

         c0 = c00*(1.0_dp - t_y) + c10*t_y
         c1 = c01*(1.0_dp - t_y) + c11*t_y

         lin_result = c0*(1.0_dp - t_z) + c1*t_z

         ! If valid, try log-space interpolation (smoother for flux)
         if (lin_result > tiny_value) then
            c00 = log(c000)*(1.0_dp - t_x) + log(c100)*t_x
            c01 = log(c001)*(1.0_dp - t_x) + log(c101)*t_x
            c10 = log(c010)*(1.0_dp - t_x) + log(c110)*t_x
            c11 = log(c011)*(1.0_dp - t_x) + log(c111)*t_x

            c0 = c00*(1.0_dp - t_y) + c10*t_y
            c1 = c01*(1.0_dp - t_y) + c11*t_y

            log_result = c0*(1.0_dp - t_z) + c1*t_z

            ! Only use log-space result if valid
            if (log_result == log_result) then  ! NaN check
               lin_result = exp(log_result)
            end if
         end if

         ! Final sanity check — fall back to nearest neighbour
         if (lin_result /= lin_result .or. lin_result <= 0.0_dp) then
            call find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                    i_x, i_y, i_z)
            lin_result = max(tiny_value, f_values_4d(i_x, i_y, i_z, lam))
         end if

         result_flux(lam) = lin_result
      end do

   end subroutine trilinear_interp_vector

   !---------------------------------------------------------------------------
   ! Scalar trilinear interpolation function (for external callers)
   !
   ! Retained for backward compatibility and for use cases that need
   ! a single-wavelength interpolation on an arbitrary 3-D grid.
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
      if (i_x < lbound(x_grid,1)) i_x = lbound(x_grid,1)
      if (i_y < lbound(y_grid,1)) i_y = lbound(y_grid,1)
      if (i_z < lbound(z_grid,1)) i_z = lbound(z_grid,1)
      if (i_x >= ubound(x_grid,1)) i_x = ubound(x_grid,1) - 1
      if (i_y >= ubound(y_grid,1)) i_y = ubound(y_grid,1) - 1
      if (i_z >= ubound(z_grid,1)) i_z = ubound(z_grid,1) - 1

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

end module linear_interp