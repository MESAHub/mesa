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

   public :: read_namelist, write_namelist, get_colors_controls, set_colors_controls

   private

   logical, dimension(max_extra_inlists) :: read_extra_colors_inlist
   character(len=strlen), dimension(max_extra_inlists) :: extra_colors_inlist_name

   character(len=256) :: instrument
   character(len=256) :: vega_sed
   character(len=256) :: stellar_atm
   character(len=256) :: colors_results_directory
   real(dp) :: metallicity
   real(dp) :: distance
   logical :: make_csv
   logical :: use_colors

   namelist /colors/ &
      instrument, &
      vega_sed, &
      stellar_atm, &
      metallicity, &
      distance, &
      make_csv, &
      colors_results_directory, &
      use_colors, &
      read_extra_colors_inlist, &
      extra_colors_inlist_name

contains

! read a "namelist" file and set parameters
   subroutine read_namelist(handle, inlist, ierr)
      integer, intent(in) :: handle
      character(len=*), intent(in) :: inlist
      integer, intent(out) :: ierr  ! 0 means AOK.
      type(Colors_General_Info), pointer :: rq
      include 'formats'
      call get_colors_ptr(handle, rq, ierr)
      if (ierr /= 0) return
      call set_default_controls
      call read_controls_file(rq, inlist, 1, ierr)
      if (ierr /= 0) return
   end subroutine read_namelist

   recursive subroutine read_controls_file(rq, filename, level, ierr)
      use iso_fortran_env, only: iostat_end
      character(*), intent(in) :: filename
      integer, intent(in) :: level
      type(Colors_General_Info), pointer, intent(inout) :: rq
      integer, intent(out) :: ierr
      logical, dimension(max_extra_inlists) :: read_extra
      character(len=strlen), dimension(max_extra_inlists) :: extra
      integer :: unit, i

      ierr = 0
      if (level >= 10) then
         write (*, *) 'ERROR: too many levels of nested extra controls inlist files'
         ierr = -1
         return
      end if

      if (len_trim(filename) > 0) then
         open (newunit=unit, file=trim(filename), &
               action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            if (level == 1) then
               ierr = 0  ! no inlist file so just use defaults
               call store_controls(rq, ierr)
            else
               write (*, *) 'Failed to open colors namelist file ', trim(filename)
            end if
            return
         end if
         read (unit, nml=colors, iostat=ierr)
         close (unit)
         if (ierr == IOSTAT_END) then  ! end-of-file means didn't find an &colors namelist
            ierr = 0
            write (*, *) 'WARNING: Failed to find colors namelist in file: ', trim(filename)
            call store_controls(rq, ierr)
            close (unit)
            return
         else if (ierr /= 0) then
            write (*, *)
            write (*, *)
            write (*, *)
            write (*, *)
            write (*, '(a)') 'Failed while trying to read colors namelist file: '//trim(filename)
            write (*, '(a)') 'Perhaps the following runtime error message will help you find the problem.'
            write (*, *)
            open (newunit=unit, file=trim(filename), action='read', &
                  delim='quote', status='old', iostat=ierr)
            read (unit, nml=colors)
            close (unit)
            return
         end if
      end if

      call store_controls(rq, ierr)

      if (len_trim(filename) == 0) return

      ! recursive calls to read other inlists
      do i = 1, max_extra_inlists
         read_extra(i) = read_extra_colors_inlist(i)
         read_extra_colors_inlist(i) = .false.
         extra(i) = extra_colors_inlist_name(i)
         extra_colors_inlist_name(i) = 'undefined'

         if (read_extra(i)) then
            call read_controls_file(rq, extra(i), level + 1, ierr)
            if (ierr /= 0) return
         end if
      end do

   end subroutine read_controls_file

   subroutine set_default_controls
      include 'colors.defaults'
   end subroutine set_default_controls

   subroutine store_controls(rq, ierr)
      type(Colors_General_Info), pointer, intent(inout) :: rq

      integer :: i
      integer, intent(out) :: ierr

      rq%instrument = instrument
      rq%vega_sed = vega_sed
      rq%stellar_atm = stellar_atm
      rq%metallicity = metallicity
      rq%distance = distance
      rq%make_csv = make_csv
      rq%colors_results_directory = colors_results_directory
      rq%use_colors = use_colors

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
      metallicity = rq%metallicity
      distance = rq%distance
      make_csv = rq%make_csv
      colors_results_directory = rq%colors_results_directory
      use_colors = rq%use_colors

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
      call store_controls(rq, ierr)
      if (ierr /= 0) return

   end subroutine set_colors_controls

end module colors_ctrls_io
