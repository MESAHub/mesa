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

module colors_ctrls_io

   use const_def, only: dp, strlen, max_extra_inlists
   use colors_def, only: Colors_General_Info, get_colors_ptr

   implicit none

   public :: read_colors_namelist, write_namelist, get_colors_controls, set_colors_controls

   private

   logical, dimension(max_extra_inlists) :: read_extra_colors_inlist
   character(len=strlen), dimension(max_extra_inlists) :: extra_colors_inlist_name

   character(len=256) :: instrument
   character(len=256) :: vega_sed
   character(len=256) :: stellar_atm
   character(len=256) :: colors_results_directory
   character(len=32) :: mag_system

   real(dp) :: distance
   logical :: make_csv
   logical :: use_colors

   namelist /colors/ &
      instrument, &
      vega_sed, &
      stellar_atm, &
      distance, &
      make_csv, &
      mag_system, &
      colors_results_directory, &
      use_colors, &
      read_extra_colors_inlist, &
      extra_colors_inlist_name

contains

! read a "namelist" file and set parameters
   subroutine read_colors_namelist(handle, inlist, ierr)
      use utils_namelist, only: read_namelist, missing_namelist_warning
      integer, intent(in) :: handle
      character(len=*), intent(in) :: inlist
      integer, intent(out) :: ierr  ! 0 means AOK.
      type(Colors_General_Info), pointer :: rq

      call get_colors_ptr(handle, rq, ierr)

      if (ierr /= 0) return

      call set_default_controls
      call read_namelist(inlist, read_colors_file, "colors", ierr, missing_namelist_warning)

      if (ierr /= 0) return

      call store_controls(rq)
   end subroutine read_colors_namelist

   subroutine read_colors_file(unit, iostat, iomsg, extra_inlists, extra_inlists_mask)
      use const_def, only: strlen, max_extra_inlists

      integer, intent(in) :: unit
      integer, intent(out) :: iostat
      character(len=strlen), intent(out) :: iomsg
      character(len=strlen), dimension(max_extra_inlists), intent(out) :: extra_inlists
      logical, dimension(max_extra_inlists), intent(out) :: extra_inlists_mask

      integer :: i

      read(unit, nml=colors, iostat=iostat, iomsg=iomsg)

      if (iostat /= 0) then
         return
      end if

      do i=1, max_extra_inlists
         extra_inlists(i) = extra_colors_inlist_name(i)
         extra_inlists_mask(i) = read_extra_colors_inlist(i)
      end do

   end subroutine read_colors_file

   subroutine set_default_controls
      include 'colors.defaults'
   end subroutine set_default_controls

   subroutine store_controls(rq)
      type(Colors_General_Info), pointer, intent(inout) :: rq

      integer :: i

      rq%instrument = instrument
      rq%vega_sed = vega_sed
      rq%stellar_atm = stellar_atm
      rq%distance = distance
      rq%make_csv = make_csv
      rq%colors_results_directory = colors_results_directory
      rq%use_colors = use_colors
      rq%mag_system = mag_system

   end subroutine store_controls

   subroutine write_namelist(handle, filename, ierr)
      integer, intent(in) :: handle
      character(*), intent(in) :: filename
      integer, intent(out) :: ierr
      type(Colors_General_Info), pointer :: rq
      integer :: iounit
      open (newunit=iounit, file=trim(filename), &
            action='write', status='replace', iostat=ierr)
      if (ierr /= 0) then
         write (*, *) 'failed to open '//trim(filename)
         return
      end if
      call get_colors_ptr(handle, rq, ierr)
      if (ierr /= 0) then
         close (iounit)
         return
      end if
      call set_controls_for_writing(rq)
      write (iounit, nml=colors, iostat=ierr)
      close (iounit)
   end subroutine write_namelist

   subroutine set_controls_for_writing(rq)
      type(Colors_General_Info), pointer, intent(inout) :: rq

      instrument = rq%instrument
      vega_sed = rq%vega_sed
      stellar_atm = rq%stellar_atm
      distance = rq%distance
      make_csv = rq%make_csv
      colors_results_directory = rq%colors_results_directory
      use_colors = rq%use_colors
      mag_system = rq%mag_system

   end subroutine set_controls_for_writing

   subroutine get_colors_controls(rq, name, val, ierr)
      use utils_lib, only: StrUpCase
      type(Colors_General_Info), pointer, intent(inout) :: rq
      character(len=*), intent(in) :: name
      character(len=512), intent(out) :: val
      integer, intent(out) :: ierr

      character(len(name) + 1) :: upper_name
      character(len=512) :: str
      integer :: iounit, iostat, ind, i

      ierr = 0

      ! First save current controls
      call set_controls_for_writing(rq)

      ! Write namelist to temporary file
      open (newunit=iounit, status='scratch')
      write (iounit, nml=colors)
      rewind (iounit)

      ! Namelists get written in capitals
      upper_name = trim(StrUpCase(name))//'='
      val = ''
      ! Search for name inside namelist
      do
         read (iounit, '(A)', iostat=iostat) str
         ind = index(trim(str), trim(upper_name))
         if (ind /= 0) then
            val = str(ind + len_trim(upper_name):len_trim(str) - 1)  ! Remove final comma and starting =
            do i = 1, len(val)
               if (val(i:i) == '"') val(i:i) = ' '
            end do
            exit
         end if
         if (is_iostat_end(iostat)) exit
      end do

      if (len_trim(val) == 0 .and. ind == 0) ierr = -1

      close (iounit)

   end subroutine get_colors_controls

   subroutine set_colors_controls(rq, name, val, ierr)
      type(Colors_General_Info), pointer, intent(inout) :: rq
      character(len=*), intent(in) :: name, val
      character(len=len(name) + len(val) + 8) :: tmp
      integer, intent(out) :: ierr

      ierr = 0

      ! First save current colors_controls
      call set_controls_for_writing(rq)

      tmp = ''
      tmp = '&colors '//trim(name)//'='//trim(val)//' /'

      ! Load into namelist
      read (tmp, nml=colors)

      ! Add to colors
      call store_controls(rq)
   end subroutine set_colors_controls

end module colors_ctrls_io
