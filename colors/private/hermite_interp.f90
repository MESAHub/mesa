! ***********************************************************************
!   Copyright (C) 2026  Niall Miller & The MESA Team
! ***********************************************************************

! ***********************************************************************
! Hermite interpolation module for spectral energy distributions (SEDs)
! ***********************************************************************

module hermite_interp
   use const_def, only: dp
   use colors_def, only: Colors_General_Info
   use colors_utils, only: dilute_flux, find_containing_cell, &
                          find_nearest_point, find_bracket_index, load_stencil
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
      real(dp), dimension(:), allocatable :: interp_flux

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
      allocate(fluxes(n_lambda))
      call dilute_flux(interp_flux, R, d, fluxes)

   end subroutine construct_sed_hermite

   !---------------------------------------------------------------------------
   ! Fallback: Build a local sub-cube from individual SED files with enough
   ! context for Hermite derivative computation, then interpolate all
   ! wavelengths in a single pass.
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
      real(dp) :: val, df_dx, df_dy, df_dz
      integer :: nx, ny, nz
      integer :: ix, iy, iz, lam
      real(dp) :: h_x(2), h_y(2), h_z(2)
      real(dp) :: hx_d(2), hy_d(2), hz_d(2)
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


      ! stencil loop — weights are invariant over lambda, so lambda is innermost
      result_flux = 0.0_dp
      do iz = 0, 1
         wz  = h_z(iz + 1)
         wzd = hz_d(iz + 1)
         do iy = 0, 1
            wy  = h_y(iy + 1)
            wyd = hy_d(iy + 1)
            do ix = 0, 1
               wx  = h_x(ix + 1)
               wxd = hx_d(ix + 1)
               do lam = 1, n_lambda
                  val = f_values_4d(i_x + ix, i_y + iy, i_z + iz, lam)

                  call compute_derivatives_at_point_4d( &
                     f_values_4d, i_x + ix, i_y + iy, i_z + iz, lam, &
                     nx, ny, nz, x_grid, y_grid, z_grid, df_dx, df_dy, df_dz)

                  result_flux(lam) = result_flux(lam) &
                        + wx*wy*wz     * val &
                        + wxd*wy*wz    * dx * df_dx &
                        + wx*wyd*wz    * dy * df_dy &
                        + wx*wy*wzd    * dz * df_dz
               end do
            end do
         end do
      end do

   end subroutine hermite_interp_vector

   !---------------------------------------------------------------------------
   ! Compute derivatives directly from the 4-D array at a given wavelength,
   ! avoiding the need to extract a 3-D slice first.
   !---------------------------------------------------------------------------
   subroutine compute_derivatives_at_point_4d(f4d, i, j, k, lam, nx, ny, nz, &
                                               x_grid, y_grid, z_grid, df_dx, df_dy, df_dz)
      real(dp), intent(in) :: f4d(:,:,:,:)
      integer, intent(in) :: i, j, k, lam, nx, ny, nz
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(out) :: df_dx, df_dy, df_dz

      ! x derivative
      if (nx == 1) then
         df_dx = 0.0_dp
      else if (i > 1 .and. i < nx) then
         df_dx = (f4d(i+1,j,k,lam) - f4d(i-1,j,k,lam)) / (x_grid(i+1) - x_grid(i-1))
      else if (i == 1) then
         df_dx = (f4d(i+1,j,k,lam) - f4d(i,j,k,lam)) / (x_grid(i+1) - x_grid(i))
      else
         df_dx = (f4d(i,j,k,lam) - f4d(i-1,j,k,lam)) / (x_grid(i) - x_grid(i-1))
      end if

      ! y derivative
      if (ny == 1) then
         df_dy = 0.0_dp
      else if (j > 1 .and. j < ny) then
         df_dy = (f4d(i,j+1,k,lam) - f4d(i,j-1,k,lam)) / (y_grid(j+1) - y_grid(j-1))
      else if (j == 1) then
         df_dy = (f4d(i,j+1,k,lam) - f4d(i,j,k,lam)) / (y_grid(j+1) - y_grid(j))
      else
         df_dy = (f4d(i,j,k,lam) - f4d(i,j-1,k,lam)) / (y_grid(j) - y_grid(j-1))
      end if

      ! z derivative
      if (nz == 1) then
         df_dz = 0.0_dp
      else if (k > 1 .and. k < nz) then
         df_dz = (f4d(i,j,k+1,lam) - f4d(i,j,k-1,lam)) / (z_grid(k+1) - z_grid(k-1))
      else if (k == 1) then
         df_dz = (f4d(i,j,k+1,lam) - f4d(i,j,k,lam)) / (z_grid(k+1) - z_grid(k))
      else
         df_dz = (f4d(i,j,k,lam) - f4d(i,j,k-1,lam)) / (z_grid(k) - z_grid(k-1))
      end if

   end subroutine compute_derivatives_at_point_4d


   function hermite_tensor_interp3d(x_val, y_val, z_val, x_grid, y_grid, &
                                    z_grid, f_values) result(f_interp)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(in) :: f_values(:, :, :)
      real(dp) :: f_interp

      integer :: i_x, i_y, i_z, ix, iy, iz
      integer :: nx, ny, nz
      real(dp) :: t_x, t_y, t_z, dx, dy, dz, f_sum
      real(dp) :: dx_values(2, 2, 2), dy_values(2, 2, 2), dz_values(2, 2, 2)
      real(dp) :: values(2, 2, 2)
      real(dp) :: h_x(2), h_y(2), h_z(2)
      real(dp) :: hx_d(2), hy_d(2), hz_d(2)

      nx = size(x_grid)
      ny = size(y_grid)
      nz = size(z_grid)

      call find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                i_x, i_y, i_z, t_x, t_y, t_z)

      if (i_x < 1 .or. i_x >= nx .or. &
          i_y < 1 .or. i_y >= ny .or. &
          i_z < 1 .or. i_z >= nz) then
         call find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
         f_interp = f_values(i_x, i_y, i_z)
         return
      end if

      dx = x_grid(i_x + 1) - x_grid(i_x)
      dy = y_grid(i_y + 1) - y_grid(i_y)
      dz = z_grid(i_z + 1) - z_grid(i_z)

      do iz = 0, 1
         do iy = 0, 1
            do ix = 0, 1
               values(ix+1, iy+1, iz+1) = f_values(i_x+ix, i_y+iy, i_z+iz)
               call compute_derivatives_at_point(f_values, i_x+ix, i_y+iy, i_z+iz, &
                                                 nx, ny, nz, x_grid, y_grid, z_grid, &
                                                 dx_values(ix+1, iy+1, iz+1), &
                                                 dy_values(ix+1, iy+1, iz+1), &
                                                 dz_values(ix+1, iy+1, iz+1))
            end do
         end do
      end do

      h_x  = [h00(t_x), h01(t_x)]
      hx_d = [h10(t_x), h11(t_x)]
      h_y  = [h00(t_y), h01(t_y)]
      hy_d = [h10(t_y), h11(t_y)]
      h_z  = [h00(t_z), h01(t_z)]
      hz_d = [h10(t_z), h11(t_z)]

      f_sum = 0.0_dp
      do iz = 1, 2
         do iy = 1, 2
            do ix = 1, 2
               f_sum = f_sum + h_x(ix)*h_y(iy)*h_z(iz)  * values(ix, iy, iz)
               f_sum = f_sum + hx_d(ix)*h_y(iy)*h_z(iz) * dx * dx_values(ix, iy, iz)
               f_sum = f_sum + h_x(ix)*hy_d(iy)*h_z(iz) * dy * dy_values(ix, iy, iz)
               f_sum = f_sum + h_x(ix)*h_y(iy)*hz_d(iz) * dz * dz_values(ix, iy, iz)
            end do
         end do
      end do

      f_interp = f_sum
   end function hermite_tensor_interp3d


   subroutine compute_derivatives_at_point(f, i, j, k, nx, ny, nz, &
                                           x_grid, y_grid, z_grid, &
                                           df_dx, df_dy, df_dz)
      real(dp), intent(in) :: f(:, :, :)
      integer, intent(in) :: i, j, k, nx, ny, nz
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      real(dp), intent(out) :: df_dx, df_dy, df_dz

      if (nx == 1) then
         df_dx = 0.0_dp
      else if (i > 1 .and. i < nx) then
         df_dx = (f(i+1,j,k) - f(i-1,j,k)) / (x_grid(i+1) - x_grid(i-1))
      else if (i == 1) then
         df_dx = (f(i+1,j,k) - f(i,j,k))   / (x_grid(i+1) - x_grid(i))
      else
         df_dx = (f(i,j,k)   - f(i-1,j,k)) / (x_grid(i)   - x_grid(i-1))
      end if

      if (ny == 1) then
         df_dy = 0.0_dp
      else if (j > 1 .and. j < ny) then
         df_dy = (f(i,j+1,k) - f(i,j-1,k)) / (y_grid(j+1) - y_grid(j-1))
      else if (j == 1) then
         df_dy = (f(i,j+1,k) - f(i,j,k))   / (y_grid(j+1) - y_grid(j))
      else
         df_dy = (f(i,j,k)   - f(i,j-1,k)) / (y_grid(j)   - y_grid(j-1))
      end if

      if (nz == 1) then
         df_dz = 0.0_dp
      else if (k > 1 .and. k < nz) then
         df_dz = (f(i,j,k+1) - f(i,j,k-1)) / (z_grid(k+1) - z_grid(k-1))
      else if (k == 1) then
         df_dz = (f(i,j,k+1) - f(i,j,k))   / (z_grid(k+1) - z_grid(k))
      else
         df_dz = (f(i,j,k)   - f(i,j,k-1)) / (z_grid(k)   - z_grid(k-1))
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