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

module colors_lib

   use const_def, only: dp, strlen
   use bolometric, only: calculate_bolometric
   use synthetic, only: calculate_synthetic
   use colors_utils, only: read_strings_from_file
   use colors_history, only: how_many_colors_history_columns, data_for_colors_history_columns

   implicit none

   private

   public :: colors_init, colors_shutdown
   public :: alloc_colors_handle, alloc_colors_handle_using_inlist, free_colors_handle
   public :: colors_ptr
   public :: colors_setup_tables, colors_setup_hooks
   ! Main functions
   public :: calculate_bolometric, calculate_synthetic
   public :: how_many_colors_history_columns, data_for_colors_history_columns
   ! Old bolometric correction functions that MESA expects (stub implementations, remove later):
   public :: get_bc_id_by_name, get_lum_band_by_id, get_abs_mag_by_id
   public :: get_bc_by_id, get_bc_name_by_id, get_bc_by_name
   public :: get_abs_bolometric_mag, get_abs_mag_by_name, get_bcs_all
   public :: get_lum_band_by_name
contains

   ! call this routine to initialize the colors module.
   ! only needs to be done once at start of run.
   ! Reads data from the 'colors' directory in the data_dir.
   ! If use_cache is true and there is a 'colors/cache' directory, it will try that first.
   ! If it doesn't find what it needs in the cache,
   ! it reads the data and writes the cache for next time.
   subroutine colors_init(use_cache, colors_cache_dir, ierr)
      use colors_def, only: colors_def_init, colors_use_cache, colors_is_initialized
      logical, intent(in) :: use_cache
      character(len=*), intent(in) :: colors_cache_dir  ! blank means use default
      integer, intent(out) :: ierr  ! 0 means AOK.
      ierr = 0
      if (colors_is_initialized) return
      call colors_def_init(colors_cache_dir)
      colors_use_cache = use_cache
      colors_is_initialized = .true.
   end subroutine colors_init

   subroutine colors_shutdown
      use colors_def, only: do_free_colors_tables, colors_is_initialized
      call do_free_colors_tables()
      colors_is_initialized = .false.
   end subroutine colors_shutdown

   ! after colors_init has finished, you can allocate a "handle".
   integer function alloc_colors_handle(ierr) result(handle)
      integer, intent(out) :: ierr  ! 0 means AOK.
      character(len=0) :: inlist
      handle = alloc_colors_handle_using_inlist(inlist, ierr)
   end function alloc_colors_handle

   integer function alloc_colors_handle_using_inlist(inlist, ierr) result(handle)
      use colors_def, only: do_alloc_colors, colors_is_initialized
      use colors_ctrls_io, only: read_namelist
      character(len=*), intent(in) :: inlist  ! empty means just use defaults.
      integer, intent(out) :: ierr  ! 0 means AOK.
      ierr = 0
      if (.not. colors_is_initialized) then
         ierr = -1
         return
      end if
      handle = do_alloc_colors(ierr)
      if (ierr /= 0) return
      call read_namelist(handle, inlist, ierr)
      if (ierr /= 0) return
      call colors_setup_tables(handle, ierr)
      call colors_setup_hooks(handle, ierr)
   end function alloc_colors_handle_using_inlist

   subroutine free_colors_handle(handle)
      ! frees the handle and all associated data
      use colors_def, only: colors_General_Info, do_free_colors
      integer, intent(in) :: handle
      call do_free_colors(handle)
   end subroutine free_colors_handle

   subroutine colors_ptr(handle, rq, ierr)

      use colors_def, only: Colors_General_Info, get_colors_ptr, colors_is_initialized

      type(colors_General_Info), pointer, intent(out) :: rq
      integer, intent(in) :: handle
      integer, intent(out):: ierr

      if (.not. colors_is_initialized) then
         ierr = -1
         return
      end if

      call get_colors_ptr(handle, rq, ierr)

   end subroutine colors_ptr

   subroutine colors_setup_tables(handle, ierr)
      use colors_def, only: colors_General_Info, get_colors_ptr, color_filter_names, num_color_filters
      ! TODO: use load_colors, only: Setup_colors_Tables
      integer, intent(in) :: handle
      integer, intent(out):: ierr

      type(colors_General_Info), pointer :: rq
      logical, parameter :: use_cache = .true.
      logical, parameter :: load_on_demand = .true.

      ierr = 0
      call get_colors_ptr(handle, rq, ierr)
      ! TODO: call Setup_colors_Tables(rq, use_cache, load_on_demand, ierr)

      ! TODO: For now, don't use cache (future feature)
      ! but rely on user specifying a single filters directory, and read it here
      call read_strings_from_file(rq, color_filter_names, num_color_filters, ierr)

   end subroutine colors_setup_tables

   subroutine colors_setup_hooks(handle, ierr)
      use colors_def, only: colors_General_Info, get_colors_ptr
      integer, intent(in) :: handle
      integer, intent(out):: ierr

      type(colors_General_Info), pointer :: rq

      ierr = 0
      call get_colors_ptr(handle, rq, ierr)

      ! TODO: currently does nothing. See kap if this feature is needed

   end subroutine colors_setup_hooks

   !-----------------------------------------------------------------------
   ! Bolometric correction interface (stub implementations)
   !-----------------------------------------------------------------------

   real(dp) function get_bc_by_name(name, log_Teff, log_g, M_div_h, ierr)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: log_Teff  ! log10 of surface temp
      real(dp), intent(in) :: log_g  ! log_10 of surface gravity
      real(dp), intent(in) :: M_div_h  ! [M/H]
      integer, intent(inout) :: ierr

      get_bc_by_name = -99.9d0
      ierr = 0
   end function get_bc_by_name

   real(dp) function get_bc_by_id(id, log_Teff, log_g, M_div_h, ierr)
      integer, intent(in) :: id
      real(dp), intent(in) :: log_Teff  ! log10 of surface temp
      real(dp), intent(in) :: log_g  ! log_10 of surface gravity
      real(dp), intent(in) :: M_div_h  ! [M/H]
      integer, intent(inout) :: ierr

      get_bc_by_id = -99.9d0
      ierr = 0
   end function get_bc_by_id

   integer function get_bc_id_by_name(name, ierr)
      character(len=*), intent(in) :: name
      integer, intent(inout) :: ierr

      get_bc_id_by_name = -1
      ierr = 0
   end function get_bc_id_by_name

   character(len=strlen) function get_bc_name_by_id(id, ierr)
      integer, intent(in) :: id
      integer, intent(inout) :: ierr

      get_bc_name_by_id = ''
      ierr = 0
   end function get_bc_name_by_id

   real(dp) function get_abs_bolometric_mag(lum)
      use const_def, only: dp
      real(dp), intent(in) :: lum  ! Luminosity in lsun units

      get_abs_bolometric_mag = -99.9d0
   end function get_abs_bolometric_mag

   real(dp) function get_abs_mag_by_name(name, log_Teff, log_g, M_div_h, lum, ierr)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: log_Teff  ! log10 of surface temp
      real(dp), intent(in) :: M_div_h  ! [M/H]
      real(dp), intent(in) :: log_g  ! log_10 of surface gravity
      real(dp), intent(in) :: lum  ! Luminosity in lsun units
      integer, intent(inout) :: ierr

      ierr = 0
      get_abs_mag_by_name = -99.9d0
   end function get_abs_mag_by_name

   real(dp) function get_abs_mag_by_id(id, log_Teff, log_g, M_div_h, lum, ierr)
      integer, intent(in) :: id
      real(dp), intent(in) :: log_Teff  ! log10 of surface temp
      real(dp), intent(in) :: log_g  ! log_10 of surface gravity
      real(dp), intent(in) :: M_div_h  ! [M/H]
      real(dp), intent(in) :: lum  ! Luminosity in lsun units
      integer, intent(inout) :: ierr

      ierr = 0
      get_abs_mag_by_id = -99.9d0
   end function get_abs_mag_by_id

   subroutine get_bcs_all(log_Teff, log_g, M_div_h, results, ierr)
      real(dp), intent(in) :: log_Teff  ! log10 of surface temp
      real(dp), intent(in) :: M_div_h  ! [M/H]
      real(dp), dimension(:), intent(out) :: results
      real(dp), intent(in) :: log_g
      integer, intent(inout) :: ierr

      ierr = 0
      results(:) = -99.d0
   end subroutine get_bcs_all

   real(dp) function get_lum_band_by_name(name, log_Teff, log_g, M_div_h, lum, ierr)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: log_Teff  ! log10 of surface temp
      real(dp), intent(in) :: M_div_h  ! [M/H]
      real(dp), intent(in) :: log_g  ! log_10 of surface gravity
      real(dp), intent(in) :: lum  ! Total luminosity in lsun units
      integer, intent(inout) :: ierr

      ierr = 0
      get_lum_band_by_name = -99.d0
   end function get_lum_band_by_name

   real(dp) function get_lum_band_by_id(id, log_Teff, log_g, M_div_h, lum, ierr)
      integer, intent(in) :: id
      real(dp), intent(in) :: log_Teff  ! log10 of surface temp
      real(dp), intent(in) :: log_g  ! log_10 of surface gravity
      real(dp), intent(in) :: M_div_h  ! [M/H]
      real(dp), intent(in) :: lum  ! Total luminosity in lsun units
      integer, intent(inout) :: ierr

      ierr = 0
      get_lum_band_by_id = -99.d0
   end function get_lum_band_by_id

end module colors_lib
