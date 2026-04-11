! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
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

!> Reading nested namelists
module utils_namelist
   implicit none
   private

   integer, parameter :: max_nested_inlists = 10

   abstract interface
      !> Read a single inlist
      !>
      !> Implementations of this interface should only read one namelist (with the unit, iostat, and iomsg passed to read) and
      !> optionally set the extra_inlists and extra_inlists_mask arguments. Each element of extra_inlists for which
      !> extra_inlists_mask is set to true will also be read in by read_namelist. If there is no need to read in extra inlists,
      !> just set all elements of extra_inlists_mask to false.
      subroutine reader(unit, iostat, iomsg, extra_inlists, extra_inlists_mask)
         use const_def, only: max_extra_inlists, strlen
         implicit none
         integer, intent(in) :: unit
         integer, intent(out) :: iostat
         character(len=strlen), intent(out) :: iomsg
         character(len=strlen), dimension(max_extra_inlists), intent(out) :: extra_inlists
         logical, dimension(max_extra_inlists), intent(out) :: extra_inlists_mask
      end subroutine reader
   end interface

   type missing_namelist
      integer, private :: action
   end type missing_namelist

   type(missing_namelist), public, parameter :: missing_namelist_error = missing_namelist(0)
   type(missing_namelist), public, parameter :: missing_namelist_warning = missing_namelist(1)
   type(missing_namelist), public, parameter :: missing_namelist_silent = missing_namelist(2)

   public :: read_namelist, missing_namelist, reader

   contains
      !> Read a nested set of namelists starting from a single file.
      !>
      !> This also handles error reporting to the user. Missing namelist
      !> entries are handled based on the value of the `missing` argument.
      subroutine read_namelist(file, r, namelist_name, ierr, missing)
         character(len=*), intent(in) :: file
         procedure(reader) :: r
         character(len=*), intent(in) :: namelist_name
         integer, intent(out) :: ierr
         type(missing_namelist), intent(in) :: missing

         call read_one_namelist(file, r, namelist_name, 1, ierr, missing)
      end subroutine read_namelist

      recursive subroutine read_one_namelist(file, r, namelist_name, level, ierr, missing)
         use const_def, only: strlen, max_extra_inlists

         character(len=*), intent(in) :: file
         procedure(reader) :: r
         character(len=*), intent(in) :: namelist_name
         integer, intent(in) :: level
         integer, intent(out) :: ierr
         type(missing_namelist), intent(in) :: missing

         integer :: iostat, unit, i
         character(len=strlen) :: iomsg
         character(len=strlen), dimension(max_extra_inlists) :: extra_inlists
         logical, dimension(max_extra_inlists) :: extra_inlists_mask

         if (level >= max_nested_inlists) then
            write(*, *) '[ERROR]: too many levels of nested ', namelist_name, ' inlist files'
            ierr = -1
            return
         end if

         open(newunit = unit, file = trim(file), action = 'read', &
            delim = 'quote', status = 'old', iostat = iostat, iomsg = iomsg)

         if (iostat /= 0) then
            write(*, *) '[ERROR]: Failed to open ', namelist_name, &
               ' namelist file "', trim(file), '". Error message: "', trim(iomsg), '"'
            ierr = -1
            return
         end if

         call r(unit, iostat, iomsg, extra_inlists, extra_inlists_mask)

         close(unit)

         if (iostat /= 0) then
            if (is_iostat_end(iostat)) then
               select case(missing%action)
               case(missing_namelist_error%action)
                  write(*, *) '[ERROR]: Failed to read ', namelist_name, &
                     ' namelist from "', trim(file), '". Namelist ', namelist_name, ' is not found'
                  ierr = -1
               case(missing_namelist_warning%action)
                  write(*, *) '[WARNING]: Failed to read ', namelist_name, &
                     ' namelist from "', trim(file), '". Namelist ', namelist_name, ' is not found'
               case(missing_namelist_silent%action)
                  ! Do nothing
               end select
               extra_inlists_mask(:) = .false.
            else
               write(*, *) '[ERROR]: Failed to read ', namelist_name, &
                  ' namelist from "', trim(file), '". Error message: "', trim(iomsg), '"'
               ierr = -1
            end if
            return
         end if

         do i=1, max_extra_inlists
            if (extra_inlists_mask(i) .and. len_trim(extra_inlists(i)) /= 0) then
               call read_one_namelist(extra_inlists(i), r, namelist_name, level + 1, ierr, missing)

               if (ierr /= 0) then
                  return
               end if
            end if
         end do

      end subroutine read_one_namelist
end module utils_namelist
