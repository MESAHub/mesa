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

   ! Colors Module control parameters
   type :: Colors_General_Info
      character(len=256) :: instrument
      character(len=256) :: vega_sed
      character(len=256) :: stellar_atm
      real(dp) :: metallicity
      real(dp) :: distance
      logical :: make_csv
      logical :: use_colors
      ! bookkeeping
      integer :: handle
      logical :: in_use
   end type Colors_General_Info

   integer, parameter :: max_colors_handles = 10
   type (Colors_General_Info), target :: colors_handles(max_colors_handles)

   logical :: colors_is_initialized = .false.

   character (len=1000) :: colors_dir, colors_cache_dir, colors_temp_cache_dir
   logical :: colors_use_cache = .true.

contains

subroutine colors_def_init(colors_cache_dir_in)
  use utils_lib, only : mkdir
  use const_def, only: mesa_data_dir, mesa_caches_dir, mesa_temp_caches_dir, use_mesa_temp_cache
  character (*), intent(in) :: colors_cache_dir_in
  integer :: i

  if (len_trim(colors_cache_dir_in) > 0) then
     colors_cache_dir = colors_cache_dir_in
  else if (len_trim(mesa_caches_dir) > 0) then
     colors_cache_dir = trim(mesa_caches_dir) // '/colors_cache'
  else
     colors_cache_dir = trim(mesa_data_dir) // '/colors_data/cache'
  end if
  call mkdir(colors_cache_dir)

  do i=1,max_colors_handles
      colors_handles(i)% handle = i
      colors_handles(i)% in_use = .false.
  end do

  colors_temp_cache_dir=trim(mesa_temp_caches_dir)//'/colors_cache'
  if(use_mesa_temp_cache) call mkdir(colors_temp_cache_dir)

end subroutine colors_def_init


integer function do_alloc_colors(ierr)
  integer, intent(out) :: ierr
  integer :: i
  ierr = 0
  do_alloc_colors = -1
  !$omp critical (colors_handle)
  do i = 1, max_colors_handles
     if (.not. colors_handles(i)% in_use) then
      colors_handles(i)% in_use = .true.
        do_alloc_colors = i
        exit
     end if
  end do
  !$omp end critical (colors_handle)
  if (do_alloc_colors == -1) then
     ierr = -1
     return
  end if
  if (colors_handles(do_alloc_colors)% handle /= do_alloc_colors) then
     ierr = -1
     return
  end if
end function do_alloc_colors


subroutine do_free_colors(handle)
  integer, intent(in) :: handle
  if (handle >= 1 .and. handle <= max_colors_handles) &
    colors_handles(handle)% in_use = .false.
end subroutine do_free_colors


subroutine get_colors_ptr(handle,rq,ierr)
  integer, intent(in) :: handle
  type (Colors_General_Info), pointer, intent(out) :: rq
  integer, intent(out):: ierr
  if (handle < 1 .or. handle > max_colors_handles) then
     ierr = -1
     return
  end if
  rq => colors_handles(handle)
  ierr = 0
end subroutine get_colors_ptr


subroutine do_free_colors_tables

  ! TODO: implement me if needed, see kap

end subroutine do_free_colors_tables



end module colors_def