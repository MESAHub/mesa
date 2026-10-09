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

! Regression test for the MESA colors module.
!
! Like the other module unit tests (see e.g. neu/test, interp_2d/test),
! this program only prints computed values -- it does not assert against
! hardcoded expected values itself. Correctness is verified externally by
! `make check`, which diffs the printed output against test_output.

program test_colors

   use const_lib, only: const_init
   use math_lib, only: math_init, pow2
   use colors_lib, only: &
      colors_init, colors_shutdown, &
      alloc_colors_handle_using_inlist, free_colors_handle, colors_ptr, &
      colors_setup_tables, colors_setup_hooks, how_many_colors_history_columns, &
      data_for_colors_history_columns
   use colors_def, only: Colors_General_Info
   use const_def, only: dp, pi, pc, rsun, Lsun
   use utils_lib, only: mesa_error
   use synthetic, only: calculate_synthetic, compute_vega_zero_point
   use hermite_interp, only: hermite_tensor_interp3d
   use hermite_interp_bounded, only: construct_sed_hermite_bounded
   use bolometric, only: calculate_bolometric_phot

   implicit none

   integer, parameter :: n_cases = 4

   real(dp), parameter :: test_teff(n_cases) = &
      [5778d0, 15000d0, 4000d0, 5778d0]
   real(dp), parameter :: test_logg(n_cases) = &
      [4.44d0, 4.0d0, 2.0d0, 4.44d0]
   real(dp), parameter :: test_meta(n_cases) = &
      [0.0d0, 0.0d0, 0.0d0, -2.0d0]
   real(dp), parameter :: test_R(n_cases) = &
      [rsun, 5d0*rsun, 20d0*rsun, rsun]

   character(len=12), parameter :: labels(n_cases) = &
      ['solar       ', 'hot_ms      ', 'cool_giant  ', 'metal_poor  ']

   character(len=32) :: my_mesa_dir
   integer :: handle, ierr, n_cols, j, k
   integer :: model_num
   type(Colors_General_Info), pointer :: cs
   character(len=80), allocatable :: col_names(:)
   real(dp), allocatable :: col_vals(:)

   my_mesa_dir = '../..'
   call const_init(my_mesa_dir, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

   call math_init()

   call colors_init(.false., '', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

   handle = alloc_colors_handle_using_inlist('', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

   call colors_ptr(handle, cs, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

   cs%use_colors = .true.
   cs%mag_system = 'Vega'

   call colors_setup_tables(handle, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

   call colors_setup_hooks(handle, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

   n_cols = how_many_colors_history_columns(handle)
   if (n_cols <= 0) call mesa_error(__FILE__, __LINE__)

   allocate (col_names(n_cols), col_vals(n_cols))
   model_num = 0

   do j = 1, n_cases
      model_num = model_num + 1
      call data_for_colors_history_columns( &
         test_teff(j), test_logg(j), test_R(j), test_meta(j), model_num, &
         handle, n_cols, col_names, col_vals, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

      ! col_names is only populated by the call above, so the header can
      ! only be printed once the first case has run.
      if (j == 1) then
         write (*, '(a)', advance='no') 'columns: Teff logg MH R_Rsun'
         do k = 1, n_cols
            write (*, '(1x,a)', advance='no') trim(col_names(k))
         end do
         write (*, '(a)') ''
      end if

      write (*, '(a,":",1x,4(1pe16.8,1x))', advance='no') &
         trim(labels(j)), test_teff(j), test_logg(j), test_meta(j), &
         test_R(j)/rsun
      do k = 1, n_cols
         write (*, '(1pe16.8,1x)', advance='no') col_vals(k)
      end do
      write (*, '(a)') ''
   end do

   call check_filter_support()
   call check_nonuniform_hermite_derivatives()
   call check_bounded_hermite_small_negative()
   call check_bounded_hermite_fallback()
   call check_bolometric_non_destructive()
   call check_bolometric_zero_point()

   deallocate (col_names, col_vals)
   call free_colors_handle(handle)
   call colors_shutdown()

contains

   subroutine check_filter_support()
      real(dp) :: sed_wave(5), filter_wave(3), sed_flux(5), transmission(3)
      real(dp) :: magnitude, zero_point
      integer :: local_ierr

      ! Filter is tabulated only over the middle three wavelengths; the
      ! bright outer two samples must not leak into either result.
      sed_wave = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      filter_wave = [2.0_dp, 3.0_dp, 4.0_dp]
      sed_flux = [100.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 100.0_dp]
      transmission = 1.0_dp

      magnitude = calculate_synthetic( &
         0.0_dp, 0.0_dp, 0.0_dp, local_ierr, sed_wave, sed_flux, &
         filter_wave, transmission, 1.0_dp, 'test.dat', .false., &
         .false., '.', 0)
      if (local_ierr /= 0) call mesa_error(__FILE__, __LINE__)

      zero_point = compute_vega_zero_point( &
         sed_wave, sed_flux, filter_wave, transmission)

      write (*, '(a,":",1x,2(1pe23.13,1x))') &
         'filter_support', magnitude, zero_point
   end subroutine check_filter_support

   subroutine check_nonuniform_hermite_derivatives()
      type(Colors_General_Info) :: bounded_cs
      real(dp) :: x_grid(4), singleton(1), values(4, 1, 1), actual
      real(dp), allocatable :: actual_wave(:), actual_flux(:)
      integer :: lam

      x_grid = [0.0_dp, 1.0_dp, 3.0_dp, 6.0_dp]
      singleton = [0.0_dp]
      values(:, 1, 1) = 100.0_dp + x_grid*x_grid

      actual = hermite_tensor_interp3d( &
         1.5_dp, 0.0_dp, 0.0_dp, x_grid, singleton, singleton, values)
      write (*, '(a,":",1x,1pe23.13)') 'nonuniform_hermite_3d', actual

      bounded_cs%cube_loaded = .true.
      bounded_cs%cube_teff_grid = x_grid
      bounded_cs%cube_logg_grid = singleton
      bounded_cs%cube_meta_grid = singleton
      bounded_cs%cube_wavelengths = [1.0_dp, 2.0_dp, 4.0_dp]
      allocate (bounded_cs%cube_flux(4, 1, 1, 3))
      do lam = 1, 3
         bounded_cs%cube_flux(:, 1, 1, lam) = values(:, 1, 1)
      end do

      call construct_sed_hermite_bounded( &
         bounded_cs, 1.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, '', &
         actual_wave, actual_flux)
      write (*, '(a,":",1x,*(1pe23.13,1x))') &
         'nonuniform_hermite_4d', actual_flux
   end subroutine check_nonuniform_hermite_derivatives

   subroutine check_bounded_hermite_small_negative()
      type(Colors_General_Info) :: bounded_cs
      real(dp), allocatable :: actual_wave(:), actual_flux(:)

      bounded_cs%cube_loaded = .true.
      bounded_cs%cube_teff_grid = [0.0_dp, 1.0_dp, 2.0_dp]
      bounded_cs%cube_logg_grid = [0.0_dp]
      bounded_cs%cube_meta_grid = [0.0_dp]
      bounded_cs%cube_wavelengths = [1.0_dp, 2.0_dp, 4.0_dp]
      allocate (bounded_cs%cube_flux(3, 1, 1, 3))
      bounded_cs%cube_flux(:, 1, 1, 1) = 100.0_dp
      bounded_cs%cube_flux(:, 1, 1, 2) = [1.0_dp, 1.0_dp, 17.000016_dp]
      bounded_cs%cube_flux(:, 1, 1, 3) = 100.0_dp

      call construct_sed_hermite_bounded( &
         bounded_cs, 0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, '', &
         actual_wave, actual_flux)
      write (*, '(a,":",1x,*(1pe23.13,1x))') &
         'bounded_hermite_small_negative', actual_flux
   end subroutine check_bounded_hermite_small_negative

   subroutine check_bounded_hermite_fallback()
      type(Colors_General_Info) :: bounded_cs
      real(dp), allocatable :: actual_wave(:), actual_flux(:)
      integer :: lam

      bounded_cs%cube_loaded = .true.
      bounded_cs%cube_teff_grid = [0.0_dp, 1.0_dp, 2.0_dp]
      bounded_cs%cube_logg_grid = [0.0_dp]
      bounded_cs%cube_meta_grid = [0.0_dp]
      bounded_cs%cube_wavelengths = [1.0_dp, 2.0_dp, 4.0_dp]
      allocate (bounded_cs%cube_flux(3, 1, 1, 3))
      do lam = 1, 3
         bounded_cs%cube_flux(:, 1, 1, lam) = [1.0_dp, 1.0_dp, 100.0_dp]
      end do

      call construct_sed_hermite_bounded( &
         bounded_cs, 0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, '', &
         actual_wave, actual_flux)
      write (*, '(a,":",1x,*(1pe23.13,1x))') &
         'bounded_hermite_fallback', actual_flux
   end subroutine check_bounded_hermite_fallback

   subroutine check_bolometric_non_destructive()
      real(dp) :: local_wavelengths(3), local_fluxes(3), original_fluxes(3)
      real(dp) :: local_mag, local_bol_flux

      local_wavelengths = [1.0_dp, 2.0_dp, 3.0_dp]
      local_fluxes = [-1.0_dp, 2.0_dp, 3.0_dp]
      original_fluxes = local_fluxes

      call calculate_bolometric_phot( &
         local_wavelengths, local_fluxes, local_mag, local_bol_flux)

      ! local_fluxes must equal original_fluxes: the negative sample is
      ! sanitized for the integral only, never written back to the caller.
      write (*, '(a,":",1x,*(1pe23.13,1x))') &
         'bolometric_input_unchanged', local_fluxes - original_fluxes
      write (*, '(a,":",1x,1pe23.13)') &
         'bolometric_sanitized_flux', local_bol_flux
   end subroutine check_bolometric_non_destructive

   subroutine check_bolometric_zero_point()
      real(dp) :: local_wavelengths(3), local_fluxes(3)
      real(dp) :: solar_flux, local_mag, local_bol_flux

      solar_flux = Lsun/(4.0_dp*pi*pow2(10.0_dp*pc))
      local_wavelengths = [1.0_dp, 2.0_dp, 3.0_dp]
      local_fluxes = solar_flux/2.0_dp

      call calculate_bolometric_phot( &
         local_wavelengths, local_fluxes, local_mag, local_bol_flux)
      write (*, '(a,":",1x,1pe23.13)') &
         'bolometric_solar_zero_point', local_mag
   end subroutine check_bolometric_zero_point

end program test_colors
