! ***********************************************************************
!
!   Copyright (C) 2010  Ed Brown & The MESA Team
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

module lodders_mod

   use const_def, only : dp, mesa_data_dir

   implicit none

   contains

   subroutine read_lodders03_data(datafile,ierr)
      use iso_fortran_env, only : iostat_end
      use chem_def
      use utils_lib, only : integer_dict_define

      character(len=*), intent(in) :: datafile
      integer, intent(out) :: ierr
      integer, parameter :: lodders_header_length = 5, max_number_isotopes = 500
      integer :: Z, A
      real(dp), dimension(max_number_isotopes) :: percent
      integer :: nentries
      real(dp) :: NSi
      integer :: iounit, ios, i
      character(len=2) :: el
      character(len=iso_name_length), dimension(max_number_isotopes) :: lodders03_isotopes
      character(len=256) :: filename

      ierr = 0
      filename = trim(mesa_data_dir)//'/chem_data/'//trim(datafile)
      open(newunit=iounit, file=trim(filename), iostat=ierr, status="old", action="read")
      if ( ierr /= 0 ) then
         write(*,*) 'chem_init: Error opening file containing Lodders (2003) table'
         write(*,*) 'filename ' // trim(filename)
         return
      end if

      ! skip the header
      do i = 1, lodders_header_length
         read(iounit,*)
      end do

      ! read in the file, setting bookmarks as we go.
      nentries = 0   ! accumulates number of spaces to hold the percentages
      do i = 1, max_number_isotopes
         read(iounit,*,iostat=ios) Z, el, A, percent(i), NSi
         if (ios == iostat_end) exit
         nentries = nentries + 1
         write(lodders03_isotopes(i), '(a,i0)') trim(el_name(Z)),A

      end do

      allocate(lodders03_tab6% isotopic_percent(nentries))
      lodders03_tab6% isotopic_percent(1:nentries) = percent(1:nentries)
      do i = 1, nentries
         call integer_dict_define(lodders03_tab6% name_dict, lodders03_isotopes(i), i, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: read_lodders03_data failed in integer_dict_define'
            call mesa_error(__FILE__,__LINE__)
         end if
      end do

      close(iounit)
   end subroutine read_lodders03_data

   function get_lodders03_isotopic_abundance(nuclei,ierr) result(percent)
      use chem_def
      use utils_lib, only : integer_dict_lookup
      character(len=*), intent(in) :: nuclei
      integer, intent(out) :: ierr
      real(dp) :: percent
      integer :: indx

      percent = 0.0d0
      if (.not.chem_has_been_initialized) then
         ierr = -9
         return
      end if

      ierr = 0
      call integer_dict_lookup(lodders03_tab6% name_dict, nuclei, indx, ierr)
      if (ierr /= 0) then
         ierr = 0
         return
      end if
      percent = lodders03_tab6% isotopic_percent(indx)

   end function get_lodders03_isotopic_abundance

end module lodders_mod
