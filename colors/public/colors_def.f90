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

module colors_def

   use const_def, only: dp

   implicit none

   public

   ! max number of SEDs to keep in the memory cache
   ! each slot holds one wavelength array (~1200 doubles ~ 10 KB),
   ! so 256 slots ~ 2.5 MB -- negligible even when the full cube
   ! cannot be allocated
   integer, parameter :: sed_mem_cache_cap = 256

   type :: filter_data
      character(len=100) :: name
      real(dp), allocatable :: wavelengths(:)
      real(dp), allocatable :: transmission(:)
      ! precomputed zero-point fluxes, computed once at init
      real(dp) :: vega_zero_point = -1.0_dp
      real(dp) :: ab_zero_point = -1.0_dp
      real(dp) :: st_zero_point = -1.0_dp
   end type filter_data

   type :: Colors_General_Info
      character(len=256) :: instrument
      character(len=256) :: vega_sed
      character(len=256) :: stellar_atm
      character(len=256) :: colors_results_directory
      character(len=256) :: mag_system
      real(dp) :: metallicity
      real(dp) :: distance
      logical :: make_csv
      logical :: sed_per_model
      logical :: use_colors
      integer :: handle
      logical :: in_use

      ! cached lookup table data
      logical :: lookup_loaded = .false.
      character(len=100), allocatable :: lu_file_names(:)
      real(dp), allocatable :: lu_logg(:)
      real(dp), allocatable :: lu_meta(:)
      real(dp), allocatable :: lu_teff(:)

      ! cached vega SED
      logical :: vega_loaded = .false.
      real(dp), allocatable :: vega_wavelengths(:)
      real(dp), allocatable :: vega_fluxes(:)

      ! cached filter data (includes precomputed zero-points)
      logical :: filters_loaded = .false.
      type(filter_data), allocatable :: filters(:)

      ! cached flux cube
      logical :: cube_loaded = .false.
      real(dp), allocatable :: cube_flux(:, :, :, :)    ! (n_teff, n_logg, n_meta, n_lambda)
      real(dp), allocatable :: cube_teff_grid(:)
      real(dp), allocatable :: cube_logg_grid(:)
      real(dp), allocatable :: cube_meta_grid(:)
      real(dp), allocatable :: cube_wavelengths(:)

      ! unique sorted grids (built once from lookup table at init)
      logical :: unique_grids_built = .false.
      real(dp), allocatable :: u_teff(:), u_logg(:), u_meta(:)

      ! grid_to_lu(i_t, i_g, i_m) gives the lookup-table row index for
      ! (u_teff(i_t), u_logg(i_g), u_meta(i_m)) -- avoids O(n_lu)
      ! nearest-neighbour searches at runtime
      logical :: grid_map_built = .false.
      integer, allocatable :: grid_to_lu(:, :, :)

      ! fallback-path caches (used only when cube_loaded == .false.)

      ! stencil cache: the extended neighbourhood around the current
      ! interpolation cell, includes derivative-context points
      ! (i-1 .. i+2 per axis, clamped to boundaries) so that
      ! hermite_tensor_interp3d gives the same result as the cube path
      logical :: stencil_valid = .false.
      integer :: stencil_i_t = -1, stencil_i_g = -1, stencil_i_m = -1
      real(dp), allocatable :: stencil_fluxes(:, :, :, :)   ! (st, sg, sm, n_lambda)
      real(dp), allocatable :: stencil_wavelengths(:)     ! (n_lambda)
      real(dp), allocatable :: stencil_teff(:)            ! subgrid values (st)
      real(dp), allocatable :: stencil_logg(:)            ! subgrid values (sg)
      real(dp), allocatable :: stencil_meta(:)            ! subgrid values (sm)

      ! canonical wavelength grid for fallback SEDs (set once on first disk
      ! read -- all SEDs in a given atmosphere grid share the same wavelengths)
      logical :: fallback_wavelengths_set = .false.
      real(dp), allocatable :: fallback_wavelengths(:)    ! (n_lambda)

      ! bounded SED memory cache (circular buffer, keyed by lu index)
      ! avoids re-reading text files for SEDs we've already parsed
      logical :: sed_mcache_init = .false.
      integer :: sed_mcache_count = 0
      integer :: sed_mcache_next = 1
      integer :: sed_mcache_nlam = 0
      integer, allocatable :: sed_mcache_keys(:)          ! (sed_mem_cache_cap)
      real(dp), allocatable :: sed_mcache_data(:, :)       ! (n_lambda, sed_mem_cache_cap)

   end type Colors_General_Info

   ! Global filter name list (shared across handles)
   integer :: num_color_filters
   character(len=100), allocatable :: color_filter_names(:)

   integer, parameter :: max_colors_handles = 10
   type(Colors_General_Info), target :: colors_handles(max_colors_handles)

   logical :: colors_is_initialized = .false.

   character(len=1000) :: colors_dir, colors_cache_dir, colors_temp_cache_dir
   logical :: colors_use_cache = .true.

contains

   subroutine colors_def_init(colors_cache_dir_in)
      use utils_lib, only: mkdir
      use const_def, only: mesa_data_dir, mesa_caches_dir, mesa_temp_caches_dir, use_mesa_temp_cache
      character(*), intent(in) :: colors_cache_dir_in
      integer :: i

      if (len_trim(colors_cache_dir_in) > 0) then
         colors_cache_dir = colors_cache_dir_in
      else if (len_trim(mesa_caches_dir) > 0) then
         colors_cache_dir = trim(mesa_caches_dir)//'/colors_cache'
      else
         colors_cache_dir = trim(mesa_data_dir)//'/colors_data/cache'
      end if
      call mkdir(colors_cache_dir)

      do i = 1, max_colors_handles
         colors_handles(i)%handle = i
         colors_handles(i)%in_use = .false.
         colors_handles(i)%lookup_loaded = .false.
         colors_handles(i)%vega_loaded = .false.
         colors_handles(i)%filters_loaded = .false.
         colors_handles(i)%cube_loaded = .false.
         colors_handles(i)%unique_grids_built = .false.
         colors_handles(i)%grid_map_built = .false.
         colors_handles(i)%stencil_valid = .false.
         colors_handles(i)%stencil_i_t = -1
         colors_handles(i)%stencil_i_g = -1
         colors_handles(i)%stencil_i_m = -1
         colors_handles(i)%sed_mcache_init = .false.
         colors_handles(i)%sed_mcache_count = 0
         colors_handles(i)%sed_mcache_next = 1
         colors_handles(i)%sed_mcache_nlam = 0
         colors_handles(i)%fallback_wavelengths_set = .false.
      end do

      colors_temp_cache_dir = trim(mesa_temp_caches_dir)//'/colors_cache'
      if (use_mesa_temp_cache) call mkdir(colors_temp_cache_dir)

   end subroutine colors_def_init

   integer function do_alloc_colors(ierr)
      integer, intent(out) :: ierr
      integer :: i
      ierr = 0
      do_alloc_colors = -1
      !$omp critical (colors_handle)
      do i = 1, max_colors_handles
         if (.not. colors_handles(i)%in_use) then
            colors_handles(i)%in_use = .true.
            do_alloc_colors = i
            exit
         end if
      end do
      !$omp end critical (colors_handle)
      if (do_alloc_colors == -1) then
         ierr = -1
         return
      end if
      if (colors_handles(do_alloc_colors)%handle /= do_alloc_colors) then
         ierr = -1
         return
      end if
   end function do_alloc_colors

   subroutine do_free_colors(handle)
      integer, intent(in) :: handle
      if (handle >= 1 .and. handle <= max_colors_handles) then
         colors_handles(handle)%in_use = .false.
         call free_colors_cache(handle)
      end if
   end subroutine do_free_colors

   subroutine free_colors_cache(handle)
      integer, intent(in) :: handle
      integer :: i

      if (handle < 1 .or. handle > max_colors_handles) return

      if (allocated(colors_handles(handle)%lu_file_names)) &
         deallocate (colors_handles(handle)%lu_file_names)
      if (allocated(colors_handles(handle)%lu_logg)) &
         deallocate (colors_handles(handle)%lu_logg)
      if (allocated(colors_handles(handle)%lu_meta)) &
         deallocate (colors_handles(handle)%lu_meta)
      if (allocated(colors_handles(handle)%lu_teff)) &
         deallocate (colors_handles(handle)%lu_teff)
      colors_handles(handle)%lookup_loaded = .false.

      if (allocated(colors_handles(handle)%vega_wavelengths)) &
         deallocate (colors_handles(handle)%vega_wavelengths)
      if (allocated(colors_handles(handle)%vega_fluxes)) &
         deallocate (colors_handles(handle)%vega_fluxes)
      colors_handles(handle)%vega_loaded = .false.

      if (allocated(colors_handles(handle)%filters)) then
         do i = 1, size(colors_handles(handle)%filters)
            if (allocated(colors_handles(handle)%filters(i)%wavelengths)) &
               deallocate (colors_handles(handle)%filters(i)%wavelengths)
            if (allocated(colors_handles(handle)%filters(i)%transmission)) &
               deallocate (colors_handles(handle)%filters(i)%transmission)
         end do
         deallocate (colors_handles(handle)%filters)
      end if
      colors_handles(handle)%filters_loaded = .false.

      if (allocated(colors_handles(handle)%cube_flux)) &
         deallocate (colors_handles(handle)%cube_flux)
      if (allocated(colors_handles(handle)%cube_teff_grid)) &
         deallocate (colors_handles(handle)%cube_teff_grid)
      if (allocated(colors_handles(handle)%cube_logg_grid)) &
         deallocate (colors_handles(handle)%cube_logg_grid)
      if (allocated(colors_handles(handle)%cube_meta_grid)) &
         deallocate (colors_handles(handle)%cube_meta_grid)
      if (allocated(colors_handles(handle)%cube_wavelengths)) &
         deallocate (colors_handles(handle)%cube_wavelengths)
      colors_handles(handle)%cube_loaded = .false.

      if (allocated(colors_handles(handle)%u_teff)) &
         deallocate (colors_handles(handle)%u_teff)
      if (allocated(colors_handles(handle)%u_logg)) &
         deallocate (colors_handles(handle)%u_logg)
      if (allocated(colors_handles(handle)%u_meta)) &
         deallocate (colors_handles(handle)%u_meta)
      colors_handles(handle)%unique_grids_built = .false.

      if (allocated(colors_handles(handle)%grid_to_lu)) &
         deallocate (colors_handles(handle)%grid_to_lu)
      colors_handles(handle)%grid_map_built = .false.

      if (allocated(colors_handles(handle)%stencil_fluxes)) &
         deallocate (colors_handles(handle)%stencil_fluxes)
      if (allocated(colors_handles(handle)%stencil_wavelengths)) &
         deallocate (colors_handles(handle)%stencil_wavelengths)
      if (allocated(colors_handles(handle)%stencil_teff)) &
         deallocate (colors_handles(handle)%stencil_teff)
      if (allocated(colors_handles(handle)%stencil_logg)) &
         deallocate (colors_handles(handle)%stencil_logg)
      if (allocated(colors_handles(handle)%stencil_meta)) &
         deallocate (colors_handles(handle)%stencil_meta)
      colors_handles(handle)%stencil_valid = .false.
      colors_handles(handle)%stencil_i_t = -1
      colors_handles(handle)%stencil_i_g = -1
      colors_handles(handle)%stencil_i_m = -1

      if (allocated(colors_handles(handle)%sed_mcache_keys)) &
         deallocate (colors_handles(handle)%sed_mcache_keys)
      if (allocated(colors_handles(handle)%sed_mcache_data)) &
         deallocate (colors_handles(handle)%sed_mcache_data)
      colors_handles(handle)%sed_mcache_init = .false.
      colors_handles(handle)%sed_mcache_count = 0
      colors_handles(handle)%sed_mcache_next = 1
      colors_handles(handle)%sed_mcache_nlam = 0

      if (allocated(colors_handles(handle)%fallback_wavelengths)) &
         deallocate (colors_handles(handle)%fallback_wavelengths)
      colors_handles(handle)%fallback_wavelengths_set = .false.

   end subroutine free_colors_cache

   subroutine get_colors_ptr(handle, rq, ierr)
      integer, intent(in) :: handle
      type(Colors_General_Info), pointer, intent(out) :: rq
      integer, intent(out):: ierr
      if (handle < 1 .or. handle > max_colors_handles) then
         ierr = -1
         return
      end if
      rq => colors_handles(handle)
      ierr = 0
   end subroutine get_colors_ptr

   subroutine do_free_colors_tables
      integer :: i

      if (allocated(color_filter_names)) deallocate (color_filter_names)

      do i = 1, max_colors_handles
         call free_colors_cache(i)
      end do

   end subroutine do_free_colors_tables

end module colors_def