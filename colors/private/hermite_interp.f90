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
   use colors_def, only: Colors_General_Info, sed_mem_cache_cap
   use colors_utils, only: dilute_flux, load_sed
   implicit none

   private
   public :: construct_sed_hermite, hermite_tensor_interp3d

contains

   !---------------------------------------------------------------------------
   ! Main entry point: Construct a SED using Hermite tensor interpolation.
   ! Data loading strategy is determined by rq%cube_loaded (set at init):
   !   cube_loaded = .true.  -> use the preloaded 4-D cube on the handle
   !   cube_loaded = .false. -> load individual SED files via the lookup table
   !---------------------------------------------------------------------------
   subroutine construct_sed_hermite(rq, teff, log_g, metallicity, R, d, &
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
         ! cell location is computed once and reused, no per-wavelength
         ! allocation or 3-D slice extraction needed.
         allocate (interp_flux(n_lambda))
         call hermite_interp_vector(teff, log_g, metallicity, &
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

   end subroutine construct_sed_hermite

   !---------------------------------------------------------------------------
   ! Fallback: Build a local sub-cube from individual SED files with enough
   ! context for Hermite derivative computation, then interpolate all
   ! wavelengths in a single pass.
   !
   ! Correctness guarantee
   ! ---------------------
   ! The cube path passes the FULL grid arrays to hermite_tensor_interp3d,
   ! so compute_derivatives_at_point can use centred differences at interior
   ! nodes.  To replicate this exactly, we load not just the 2x2x2 cell
   ! corners but also one extra grid point on each side when available
   ! (the "derivative stencil").  The resulting sub-grid has 2-4 points per
   ! axis, and the interpolator sees the same derivative context it would
   ! get from the full cube.
   !
   ! Performance
   ! -----------
   ! * grid_to_lu map      -- O(1) lookup per stencil point (vs O(n_lu) scan)
   ! * SED memory cache     -- file I/O only on first visit to a grid point
   ! * Stencil cache        -- no work at all when the cell hasn't changed
   ! * Vectorised wavelength loop -- cell geometry computed once, reused
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
      integer :: it, ig, im, lu_idx, i
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
         ! Determine the extended stencil bounds:
         ! For each axis, include one point before and after the cell
         ! when available, so that centred differences match the cube.
         nt = size(rq%u_teff)
         ng = size(rq%u_logg)
         nm = size(rq%u_meta)

         if (nt < 2) then
            lo_t = 1; hi_t = 1
         else
            lo_t = max(1,  i_t - 1)
            hi_t = min(nt, i_t + 2)
         end if

         if (ng < 2) then
            lo_g = 1; hi_g = 1
         else
            lo_g = max(1,  i_g - 1)
            hi_g = min(ng, i_g + 2)
         end if

         if (nm < 2) then
            lo_m = 1; hi_m = 1
         else
            lo_m = max(1,  i_m - 1)
            hi_m = min(nm, i_m + 2)
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
      call hermite_interp_vector(teff, log_g, metallicity, &
                                  rq%stencil_teff, rq%stencil_logg, rq%stencil_meta, &
                                  rq%stencil_fluxes, n_lambda, interp_flux)

   end subroutine construct_sed_from_files

   !---------------------------------------------------------------------------
   ! Load SEDs for every point in the stencil, using the memory cache.
   ! Populates rq%stencil_fluxes and rq%stencil_wavelengths.
   !
   ! Wavelengths are stored once on the handle (rq%fallback_wavelengths)
   ! on the first disk read and reused thereafter — all SEDs in a given
   ! atmosphere grid share the same wavelength array.
   !---------------------------------------------------------------------------
   subroutine load_stencil(rq, resolved_dir, lo_t, hi_t, lo_g, hi_g, lo_m, hi_m)
      type(Colors_General_Info), intent(inout) :: rq
      character(len=*), intent(in) :: resolved_dir
      integer, intent(in) :: lo_t, hi_t, lo_g, hi_g, lo_m, hi_m

      integer :: st, sg, sm, n_lambda, lu_idx
      integer :: it, ig, im
      real(dp), dimension(:), allocatable :: sed_flux

      st = hi_t - lo_t + 1
      sg = hi_g - lo_g + 1
      sm = hi_m - lo_m + 1

      ! Free previous stencil flux data
      if (allocated(rq%stencil_fluxes)) deallocate(rq%stencil_fluxes)

      n_lambda = 0

      do it = lo_t, hi_t
         do ig = lo_g, hi_g
            do im = lo_m, hi_m
               lu_idx = rq%grid_to_lu(it, ig, im)

               call load_sed_cached(rq, resolved_dir, lu_idx, sed_flux)

               if (n_lambda == 0) then
                  n_lambda = size(sed_flux)
                  allocate (rq%stencil_fluxes(st, sg, sm, n_lambda))
               end if

               rq%stencil_fluxes(it - lo_t + 1, ig - lo_g + 1, im - lo_m + 1, :) = &
                  sed_flux(1:n_lambda)

               if (allocated(sed_flux)) deallocate(sed_flux)
            end do
         end do
      end do

      ! Set stencil wavelengths from the canonical copy on the handle
      if (allocated(rq%stencil_wavelengths)) deallocate(rq%stencil_wavelengths)
      allocate (rq%stencil_wavelengths(n_lambda))
      rq%stencil_wavelengths = rq%fallback_wavelengths(1:n_lambda)

   end subroutine load_stencil

   !---------------------------------------------------------------------------
   ! Retrieve an SED flux from the memory cache, or load from disk on miss.
   ! Uses a bounded circular buffer (sed_mem_cache_cap slots).
   !
   ! On the first disk read, the wavelength array is stored once on the
   ! handle as rq%fallback_wavelengths (all SEDs in a given atmosphere
   ! grid share the same wavelength array).  Only the flux is cached
   ! and returned.
   !---------------------------------------------------------------------------
   subroutine load_sed_cached(rq, resolved_dir, lu_idx, flux)
      type(Colors_General_Info), intent(inout) :: rq
      character(len=*), intent(in) :: resolved_dir
      integer, intent(in) :: lu_idx
      real(dp), dimension(:), allocatable, intent(out) :: flux

      integer :: slot, n_lam
      character(len=512) :: filepath
      real(dp), dimension(:), allocatable :: sed_wave

      ! Initialise the cache on first call
      if (.not. rq%sed_mcache_init) then
         allocate (rq%sed_mcache_keys(sed_mem_cache_cap))
         rq%sed_mcache_keys = 0   ! 0 means empty slot
         rq%sed_mcache_count = 0
         rq%sed_mcache_next = 1
         rq%sed_mcache_nlam = 0
         rq%sed_mcache_init = .true.
      end if

      ! Search for a cache hit (linear scan over a small array)
      do slot = 1, rq%sed_mcache_count
         if (rq%sed_mcache_keys(slot) == lu_idx) then
            ! Hit — return cached flux
            n_lam = rq%sed_mcache_nlam
            allocate (flux(n_lam))
            flux = rq%sed_mcache_data(:, slot)
            return
         end if
      end do

      ! Miss — load from disk
      filepath = trim(resolved_dir)//'/'//trim(rq%lu_file_names(lu_idx))
      call load_sed(filepath, lu_idx, sed_wave, flux)

      ! Store the canonical wavelength array on the handle (once only)
      if (.not. rq%fallback_wavelengths_set) then
         n_lam = size(sed_wave)
         allocate (rq%fallback_wavelengths(n_lam))
         rq%fallback_wavelengths = sed_wave
         rq%fallback_wavelengths_set = .true.
      end if
      if (allocated(sed_wave)) deallocate(sed_wave)

      ! Store flux in the cache
      n_lam = size(flux)
      if (rq%sed_mcache_nlam == 0) then
         ! First SED ever loaded — set the wavelength count and allocate data
         rq%sed_mcache_nlam = n_lam
         allocate (rq%sed_mcache_data(n_lam, sed_mem_cache_cap))
      end if

      ! Write to the next slot (circular)
      slot = rq%sed_mcache_next
      rq%sed_mcache_keys(slot) = lu_idx
      rq%sed_mcache_data(:, slot) = flux(1:rq%sed_mcache_nlam)

      if (rq%sed_mcache_count < sed_mem_cache_cap) then
         rq%sed_mcache_count = rq%sed_mcache_count + 1
      end if
      rq%sed_mcache_next = mod(slot, sed_mem_cache_cap) + 1

   end subroutine load_sed_cached

   !---------------------------------------------------------------------------
   ! Vectorised Hermite interpolation over all wavelengths.
   !
   ! The cell location (i_x, i_y, i_z, t_x, t_y, t_z) depends only on
   ! (teff, logg, meta) and the sub-grids — not on wavelength.  Computing
   ! it once and reusing across all n_lambda samples eliminates redundant
   ! binary searches and basis-function evaluations.
   !---------------------------------------------------------------------------
   subroutine hermite_interp_vector(x_val, y_val, z_val, &
                                     x_grid, y_grid, z_grid, &
                                     f_values_4d, n_lambda, result_flux)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in) :: f_values_4d(:,:,:,:)   ! (nx, ny, nz, n_lambda)
      integer, intent(in) :: n_lambda
      real(dp), intent(out) :: result_flux(n_lambda)

      integer :: i_x, i_y, i_z
      real(dp) :: t_x, t_y, t_z
      real(dp) :: dx, dy, dz
      integer :: nx, ny, nz
      integer :: ix, iy, iz, lam
      real(dp) :: h_x(2), h_y(2), h_z(2)
      real(dp) :: hx_d(2), hy_d(2), hz_d(2)
      real(dp) :: val, df_dx, df_dy, df_dz, s
      real(dp) :: wx, wy, wz, wxd, wyd, wzd

      nx = size(x_grid)
      ny = size(y_grid)
      nz = size(z_grid)

      ! Find containing cell (done once for all wavelengths)
      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      ! If outside grid, use nearest point for all wavelengths
      if (i_x < 1 .or. i_x >= nx .or. &
          i_y < 1 .or. i_y >= ny .or. &
          i_z < 1 .or. i_z >= nz) then

         call find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
         do lam = 1, n_lambda
            result_flux(lam) = f_values_4d(i_x, i_y, i_z, lam)
         end do
         return
      end if

      ! Grid cell spacing
      dx = x_grid(i_x + 1) - x_grid(i_x)
      dy = y_grid(i_y + 1) - y_grid(i_y)
      dz = z_grid(i_z + 1) - z_grid(i_z)

      ! Precompute Hermite basis functions (same for all wavelengths)
      h_x  = [h00(t_x), h01(t_x)]
      hx_d = [h10(t_x), h11(t_x)]
      h_y  = [h00(t_y), h01(t_y)]
      hy_d = [h10(t_y), h11(t_y)]
      h_z  = [h00(t_z), h01(t_z)]
      hz_d = [h10(t_z), h11(t_z)]

      ! Loop over wavelengths — the hot loop
      do lam = 1, n_lambda
         s = 0.0_dp
         do iz = 0, 1
            wz  = h_z(iz + 1)
            wzd = hz_d(iz + 1)
            do iy = 0, 1
               wy  = h_y(iy + 1)
               wyd = hy_d(iy + 1)
               do ix = 0, 1
                  wx  = h_x(ix + 1)
                  wxd = hx_d(ix + 1)

                  val = f_values_4d(i_x + ix, i_y + iy, i_z + iz, lam)

                  call compute_derivatives_at_point_4d( &
                     f_values_4d, i_x + ix, i_y + iy, i_z + iz, lam, &
                     nx, ny, nz, dx, dy, dz, df_dx, df_dy, df_dz)

                  s = s + wx*wy*wz     * val &
                        + wxd*wy*wz    * dx * df_dx &
                        + wx*wyd*wz    * dy * df_dy &
                        + wx*wy*wzd    * dz * df_dz
               end do
            end do
         end do
         result_flux(lam) = s
      end do

   end subroutine hermite_interp_vector

   !---------------------------------------------------------------------------
   ! Compute derivatives directly from the 4-D array at a given wavelength,
   ! avoiding the need to extract a 3-D slice first.
   !---------------------------------------------------------------------------
   subroutine compute_derivatives_at_point_4d(f4d, i, j, k, lam, nx, ny, nz, &
                                               dx, dy, dz, df_dx, df_dy, df_dz)
      real(dp), intent(in) :: f4d(:,:,:,:)
      integer, intent(in) :: i, j, k, lam, nx, ny, nz
      real(dp), intent(in) :: dx, dy, dz
      real(dp), intent(out) :: df_dx, df_dy, df_dz

      ! x derivative
      if (dx < 1.0e-30_dp) then
         df_dx = 0.0_dp
      else if (i > 1 .and. i < nx) then
         df_dx = (f4d(i + 1, j, k, lam) - f4d(i - 1, j, k, lam)) / (2.0_dp * dx)
      else if (i == 1) then
         df_dx = (f4d(i + 1, j, k, lam) - f4d(i, j, k, lam)) / dx
      else
         df_dx = (f4d(i, j, k, lam) - f4d(i - 1, j, k, lam)) / dx
      end if

      ! y derivative
      if (dy < 1.0e-30_dp) then
         df_dy = 0.0_dp
      else if (j > 1 .and. j < ny) then
         df_dy = (f4d(i, j + 1, k, lam) - f4d(i, j - 1, k, lam)) / (2.0_dp * dy)
      else if (j == 1) then
         df_dy = (f4d(i, j + 1, k, lam) - f4d(i, j, k, lam)) / dy
      else
         df_dy = (f4d(i, j, k, lam) - f4d(i, j - 1, k, lam)) / dy
      end if

      ! z derivative
      if (dz < 1.0e-30_dp) then
         df_dz = 0.0_dp
      else if (k > 1 .and. k < nz) then
         df_dz = (f4d(i, j, k + 1, lam) - f4d(i, j, k - 1, lam)) / (2.0_dp * dz)
      else if (k == 1) then
         df_dz = (f4d(i, j, k + 1, lam) - f4d(i, j, k, lam)) / dz
      else
         df_dz = (f4d(i, j, k, lam) - f4d(i, j, k - 1, lam)) / dz
      end if

   end subroutine compute_derivatives_at_point_4d

   !---------------------------------------------------------------------------
   ! Find the bracketing index i such that grid(i) <= val < grid(i+1)
   !---------------------------------------------------------------------------
   subroutine find_bracket_index(grid, val, idx)
      real(dp), intent(in) :: grid(:), val
      integer, intent(out) :: idx
      integer :: n

      n = size(grid)
      if (n < 2) then
         idx = 1
         return
      end if

      if (val <= grid(1)) then
         idx = 1
      else if (val >= grid(n)) then
         idx = n - 1
      else
         idx = 1
         do while (idx < n - 1 .and. grid(idx + 1) <= val)
            idx = idx + 1
         end do
      end if
   end subroutine find_bracket_index



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
      real(dp) :: hx_d(2), hy_d(2), hz_d(2)

      ! Find containing cell and parameter values
      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      ! If outside grid, use nearest point
      if (i_x < 1 .or. i_x >= size(x_grid) .or. &
          i_y < 1 .or. i_y >= size(y_grid) .or. &
          i_z < 1 .or. i_z >= size(z_grid)) then

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

      ! Precompute Hermite basis functions and derivatives
      h_x  = [h00(t_x), h01(t_x)]
      hx_d = [h10(t_x), h11(t_x)]
      h_y  = [h00(t_y), h01(t_y)]
      hy_d = [h10(t_y), h11(t_y)]
      h_z  = [h00(t_z), h01(t_z)]
      hz_d = [h10(t_z), h11(t_z)]

      ! Final interpolation sum
      sum = 0.0_dp
      do iz = 1, 2
         do iy = 1, 2
            do ix = 1, 2
               sum = sum + h_x(ix)*h_y(iy)*h_z(iz)     * values(ix, iy, iz)
               sum = sum + hx_d(ix)*h_y(iy)*h_z(iz)    * dx * dx_values(ix, iy, iz)
               sum = sum + h_x(ix)*hy_d(iy)*h_z(iz)    * dy * dy_values(ix, iy, iz)
               sum = sum + h_x(ix)*h_y(iy)*hz_d(iz)    * dz * dz_values(ix, iy, iz)
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
      logical :: dummy_axis

      n = size(x)

      ! Detect dummy axis: all values == 0, 999, or -999
      dummy_axis = all(x == 0.0_dp) .or. all(x == 999.0_dp) .or. all(x == -999.0_dp)

      if (dummy_axis) then
         ! Collapse axis: always use first point, no interpolation
         i = 1
         t = 0.0_dp
         return
      end if

      ! ---------- ORIGINAL CODE BELOW ----------------

      if (val <= x(1)) then
         i = 1
         t = 0.0_dp
         return
      else if (val >= x(n)) then
         i = n - 1
         t = 1.0_dp
         return
      end if

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
      if (abs(x(i + 1) - x(i)) < 1.0e-30_dp) then
         t = 0.0_dp  ! degenerate interval — no interpolation needed
      else
         t = (val - x(i))/(x(i + 1) - x(i))
      end if
   end subroutine find_interval


   !---------------------------------------------------------------------------
   ! Find the nearest grid point
   !---------------------------------------------------------------------------
   subroutine find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      integer, intent(out) :: i_x, i_y, i_z

      ! Find nearest grid points using intrinsic minloc
      i_x = minloc(abs(x_val - x_grid), 1)
      i_y = minloc(abs(y_val - y_grid), 1)
      i_z = minloc(abs(z_val - z_grid), 1)
   end subroutine find_nearest_point

   !---------------------------------------------------------------------------
   ! Compute derivatives at a grid point (3-D version, used by scalar path)
   !---------------------------------------------------------------------------
   subroutine compute_derivatives_at_point(f, i, j, k, nx, ny, nz, dx, dy, dz, &
                                           df_dx, df_dy, df_dz)
      real(dp), intent(in) :: f(:, :, :)
      integer, intent(in) :: i, j, k, nx, ny, nz
      real(dp), intent(in) :: dx, dy, dz
      real(dp), intent(out) :: df_dx, df_dy, df_dz

      ! Compute x derivative using centered differences where possible
      if (dx < 1.0e-30_dp) then
         df_dx = 0.0_dp  ! degenerate axis
      else if (i > 1 .and. i < nx) then
         df_dx = (f(i + 1, j, k) - f(i - 1, j, k))/(2.0_dp*dx)
      else if (i == 1) then
         df_dx = (f(i + 1, j, k) - f(i, j, k))/dx
      else ! i == nx
         df_dx = (f(i, j, k) - f(i - 1, j, k))/dx
      end if

      ! Compute y derivative using centered differences where possible
      if (dy < 1.0e-30_dp) then
         df_dy = 0.0_dp  ! degenerate axis
      else if (j > 1 .and. j < ny) then
         df_dy = (f(i, j + 1, k) - f(i, j - 1, k))/(2.0_dp*dy)
      else if (j == 1) then
         df_dy = (f(i, j + 1, k) - f(i, j, k))/dy
      else ! j == ny
         df_dy = (f(i, j, k) - f(i, j - 1, k))/dy
      end if

      ! Compute z derivative using centered differences where possible
      if (dz < 1.0e-30_dp) then
         df_dz = 0.0_dp  ! degenerate axis
      else if (k > 1 .and. k < nz) then
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