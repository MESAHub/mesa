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

   ! Make everything in this module public by default
   public

   ! Type to hold individual filter data (cached)
   type :: filter_data
      character(len=100) :: name
      real(dp), allocatable :: wavelengths(:)
      real(dp), allocatable :: transmission(:)
   end type filter_data

   ! Colors Module control parameters
   type :: Colors_General_Info
      character(len=256) :: instrument
      character(len=256) :: vega_sed
      character(len=256) :: stellar_atm
      character(len=256) :: colors_results_directory
      character(len=256) :: mag_system
      real(dp) :: metallicity
      real(dp) :: distance
      logical :: make_csv
      logical :: use_colors
      ! bookkeeping
      integer :: handle
      logical :: in_use

      ! Cached lookup table data
      logical :: lookup_loaded = .false.
      character(len=100), allocatable :: lu_file_names(:)
      real(dp), allocatable :: lu_logg(:)
      real(dp), allocatable :: lu_meta(:)
      real(dp), allocatable :: lu_teff(:)

      ! Cached Vega SED
      logical :: vega_loaded = .false.
      real(dp), allocatable :: vega_wavelengths(:)
      real(dp), allocatable :: vega_fluxes(:)

      ! Cached filter data
      logical :: filters_loaded = .false.
      type(filter_data), allocatable :: filters(:)

   end type Colors_General_Info

   ! TODO: Use handles/caching in the future once we have more colors tables
   ! For now, we will just point to a single file
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

      ! Free lookup table arrays
      if (allocated(colors_handles(handle)%lu_file_names)) &
         deallocate(colors_handles(handle)%lu_file_names)
      if (allocated(colors_handles(handle)%lu_logg)) &
         deallocate(colors_handles(handle)%lu_logg)
      if (allocated(colors_handles(handle)%lu_meta)) &
         deallocate(colors_handles(handle)%lu_meta)
      if (allocated(colors_handles(handle)%lu_teff)) &
         deallocate(colors_handles(handle)%lu_teff)
      colors_handles(handle)%lookup_loaded = .false.

      ! Free Vega SED arrays
      if (allocated(colors_handles(handle)%vega_wavelengths)) &
         deallocate(colors_handles(handle)%vega_wavelengths)
      if (allocated(colors_handles(handle)%vega_fluxes)) &
         deallocate(colors_handles(handle)%vega_fluxes)
      colors_handles(handle)%vega_loaded = .false.

      ! Free filter data arrays
      if (allocated(colors_handles(handle)%filters)) then
         do i = 1, size(colors_handles(handle)%filters)
            if (allocated(colors_handles(handle)%filters(i)%wavelengths)) &
               deallocate(colors_handles(handle)%filters(i)%wavelengths)
            if (allocated(colors_handles(handle)%filters(i)%transmission)) &
               deallocate(colors_handles(handle)%filters(i)%transmission)
         end do
         deallocate(colors_handles(handle)%filters)
      end if
      colors_handles(handle)%filters_loaded = .false.

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

      ! Free the filter names array
      if (allocated(color_filter_names)) deallocate(color_filter_names)

      ! Free cached data for all handles
      do i = 1, max_colors_handles
         call free_colors_cache(i)
      end do

   end subroutine do_free_colors_tables

end module colors_def