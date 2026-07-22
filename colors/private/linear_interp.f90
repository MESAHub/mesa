! ***********************************************************************
!
!   Copyright (C) 2026  Niall Miller & The MESA Team
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

! linear interpolation for SEDs
!
! data-loading strategy selected by rq%cube_loaded:
!   .true.  -> use the preloaded 4-D flux cube on the handle
!   .false. -> load individual SED files via the lookup table (fallback)

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

   ! main entry point -- construct a SED using trilinear interpolation
   ! strategy controlled by rq%cube_loaded (set at init)
   subroutine construct_sed_linear(rq, teff, log_g, metallicity, R, d, &
                                   stellar_model_dir, wavelengths, fluxes)
      type(Colors_General_Info), intent(inout) :: rq
      real(dp), intent(in) :: teff, log_g, metallicity, R, d
      character(len=*), intent(in) :: stellar_model_dir
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, fluxes

      integer :: n_lambda
      real(dp), dimension(:), allocatable :: interp_flux, diluted_flux

      if (rq%cube_loaded) then
         ! fast path: use preloaded cube from handle
         n_lambda = size(rq%cube_wavelengths)

         allocate (wavelengths(n_lambda))
         wavelengths = rq%cube_wavelengths

         ! vectorised interpolation -- cell location computed once, reused across n_lambda
         allocate (interp_flux(n_lambda))
         call trilinear_interp_vector(teff, log_g, metallicity, &
                                      rq%cube_teff_grid, rq%cube_logg_grid, &
                                      rq%cube_meta_grid, &
                                      rq%cube_flux, n_lambda, interp_flux)
      else
         ! fallback path: load individual SED files from lookup table
         call construct_sed_from_files(rq, teff, log_g, metallicity, &
                                       stellar_model_dir, interp_flux, wavelengths)
         n_lambda = size(wavelengths)
      end if

      allocate (diluted_flux(n_lambda))
      call dilute_flux(interp_flux, R, d, diluted_flux)
      fluxes = diluted_flux

   end subroutine construct_sed_linear

   ! fallback: build a 2x2x2 sub-cube from SED files, then trilinear-interpolate
   ! unlike hermite, no derivative context needed -- stencil is exactly the 2x2x2 cell corners
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

      ! ensure the grid-to-lu mapping exists (built once, then reused)
      if (.not. rq%grid_map_built) call build_grid_to_lu_map(rq)

      ! find bracketing cell in the unique grids
      call find_bracket_index(rq%u_teff, teff, i_t)
      call find_bracket_index(rq%u_logg, log_g, i_g)
      call find_bracket_index(rq%u_meta, metallicity, i_m)

      ! check if the stencil cache is still valid for this cell
      need_reload = .true.
      if (rq%stencil_valid .and. &
          i_t == rq%stencil_i_t .and. &
          i_g == rq%stencil_i_g .and. &
          i_m == rq%stencil_i_m) then
         need_reload = .false.
      end if

      if (need_reload) then
         ! trilinear needs exactly the 2x2x2 cell corners -- no extension
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

         ! load SEDs for every stencil point (using memory cache)
         call load_stencil(rq, resolved_dir, lo_t, hi_t, lo_g, hi_g, lo_m, hi_m)

         ! store subgrid arrays on the handle
         if (allocated(rq%stencil_teff)) deallocate (rq%stencil_teff)
         if (allocated(rq%stencil_logg)) deallocate (rq%stencil_logg)
         if (allocated(rq%stencil_meta)) deallocate (rq%stencil_meta)

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

      n_lambda = size(rq%stencil_wavelengths)
      allocate (wavelengths(n_lambda))
      wavelengths = rq%stencil_wavelengths

      allocate (interp_flux(n_lambda))
      call trilinear_interp_vector(teff, log_g, metallicity, &
                                   rq%stencil_teff, rq%stencil_logg, rq%stencil_meta, &
                                   rq%stencil_fluxes, n_lambda, interp_flux)

   end subroutine construct_sed_from_files

   ! vectorised trilinear interpolation over all wavelengths
   ! cell location depends only on (teff, logg, meta) -- computed once, reused across n_lambda
   subroutine trilinear_interp_vector(x_val, y_val, z_val, &
                                      x_grid, y_grid, z_grid, &
                                      f_values_4d, n_lambda, result_flux)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in) :: f_values_4d(:, :, :, :)   ! (nx, ny, nz, n_lambda)
      integer, intent(in) :: n_lambda
      real(dp), intent(out) :: result_flux(n_lambda)

      integer :: i_x, i_y, i_z, lam
      integer :: nx, ny, nz
      integer :: ix_max, iy_max, iz_max
      integer :: ox, oy, oz
      real(dp) :: t_x, t_y, t_z
      real(dp) :: w_x(0:1), w_y(0:1), w_z(0:1)
      real(dp) :: accum, weight
      real(dp), parameter :: tiny_value = 1.0d-300

      nx = size(x_grid)
      ny = size(y_grid)
      nz = size(z_grid)

      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      call setup_axis(nx, x_val, x_grid, i_x, t_x, ix_max, w_x)
      call setup_axis(ny, y_val, y_grid, i_y, t_y, iy_max, w_y)
      call setup_axis(nz, z_val, z_grid, i_z, t_z, iz_max, w_z)

      do lam = 1, n_lambda
         accum = 0.0_dp

         do oz = 0, iz_max
            do oy = 0, iy_max
               do ox = 0, ix_max
                  weight = w_x(ox)*w_y(oy)*w_z(oz)
                  accum = accum + weight*f_values_4d(i_x + ox, i_y + oy, i_z + oz, lam)
               end do
            end do
         end do

         result_flux(lam) = max(tiny_value, accum)
      end do

   end subroutine trilinear_interp_vector

   ! scalar trilinear interpolation (external callers / single-wavelength use)
   ! retained for backward compatibility
   function trilinear_interp(x_val, y_val, z_val, x_grid, y_grid, z_grid, f_values) result(f_interp)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in) :: f_values(:, :, :)
      real(dp) :: f_interp

      integer :: i_x, i_y, i_z
      integer :: nx, ny, nz
      integer :: ix_max, iy_max, iz_max
      integer :: ox, oy, oz
      real(dp) :: t_x, t_y, t_z
      real(dp) :: w_x(0:1), w_y(0:1), w_z(0:1)
      real(dp) :: accum, weight
      real(dp), parameter :: tiny_value = 1.0d-300

      nx = size(x_grid)
      ny = size(y_grid)
      nz = size(z_grid)

      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      call setup_axis(nx, x_val, x_grid, i_x, t_x, ix_max, w_x)
      call setup_axis(ny, y_val, y_grid, i_y, t_y, iy_max, w_y)
      call setup_axis(nz, z_val, z_grid, i_z, t_z, iz_max, w_z)

      accum = 0.0_dp

      do oz = 0, iz_max
         do oy = 0, iy_max
            do ox = 0, ix_max
               weight = w_x(ox)*w_y(oy)*w_z(oz)
               accum = accum + weight*f_values(i_x + ox, i_y + oy, i_z + oz)
            end do
         end do
      end do

      f_interp = max(tiny_value, accum)

   end function trilinear_interp

   subroutine setup_axis(n, x_val, grid, i, t, i_max, w)
      integer, intent(in) :: n
      real(dp), intent(in) :: x_val
      real(dp), intent(in) :: grid(:)
      integer, intent(inout) :: i
      real(dp), intent(inout) :: t
      integer, intent(out) :: i_max
      real(dp), intent(out) :: w(0:1)

      if (n <= 1) then
         i = 1
         t = 0.0_dp
         i_max = 0
         w(0) = 1.0_dp
         w(1) = 0.0_dp
         return
      end if

      if (x_val <= grid(1)) then
         i = 1
         t = 0.0_dp
      else if (x_val >= grid(n)) then
         i = n - 1
         t = 1.0_dp
      else
         if (i < 1) i = 1
         if (i >= n) i = n - 1
         t = max(0.0_dp, min(1.0_dp, t))
      end if

      i_max = 1
      w(0) = 1.0_dp - t
      w(1) = t
   end subroutine setup_axis

end module linear_interp
