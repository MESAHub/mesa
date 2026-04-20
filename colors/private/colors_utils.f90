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

module colors_utils
   use const_def, only: dp, strlen, mesa_dir
   use colors_def, only: Colors_General_Info
   use utils_lib, only: mesa_error

   implicit none

   public :: dilute_flux, trapezoidal_integration, romberg_integration, &
             simpson_integration, load_sed, load_filter, load_vega_sed, &
             load_lookup_table, remove_dat
contains

   !---------------------------------------------------------------------------
   ! Apply dilution factor to convert surface flux to observed flux
   !---------------------------------------------------------------------------
   subroutine dilute_flux(surface_flux, R, d, calibrated_flux)
      real(dp), intent(in)  :: surface_flux(:)
      real(dp), intent(in)  :: R, d  ! R = stellar radius, d = distance (both in the same units, e.g., cm)
      real(dp), intent(out) :: calibrated_flux(:)

      ! Check that the output array has the same size as the input
      if (size(calibrated_flux) /= size(surface_flux)) then
         print *, "Error in dilute_flux: Output array must have the same size as input array."
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Apply the dilution factor (R/d)^2 to each element
      calibrated_flux = surface_flux*((R/d)**2)
   end subroutine dilute_flux

   !###########################################################
   !## MATHS
   !###########################################################

   !****************************
   !Trapezoidal and Simpson Integration For Flux Calculation
   !****************************

   subroutine trapezoidal_integration(x, y, result)
      real(dp), dimension(:), intent(in) :: x, y
      real(dp), intent(out) :: result

      integer :: i, n
      real(dp) :: sum

      n = size(x)
      sum = 0.0_dp

      ! Validate input sizes
      if (size(x) /= size(y)) then
         print *, "Error: x and y arrays must have the same size."
         call mesa_error(__FILE__, __LINE__)
      end if

      if (size(x) < 2) then
         print *, "Error: x and y arrays must have at least 2 elements."
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Perform trapezoidal integration
      do i = 1, n - 1
         sum = sum + 0.5_dp*(x(i + 1) - x(i))*(y(i + 1) + y(i))
      end do

      result = sum
   end subroutine trapezoidal_integration

   subroutine simpson_integration(x, y, result)
      real(dp), dimension(:), intent(in) :: x, y
      real(dp), intent(out) :: result

      integer :: i, n
      real(dp) :: sum, h1, h2, f1, f2, f0

      n = size(x)
      sum = 0.0_dp

      ! Validate input sizes
      if (size(x) /= size(y)) then
         print *, "Error: x and y arrays must have the same size."
         call mesa_error(__FILE__, __LINE__)
      end if

      if (size(x) < 2) then
         print *, "Error: x and y arrays must have at least 2 elements."
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Perform adaptive Simpson's rule
      do i = 1, n - 2, 2
         h1 = x(i + 1) - x(i)       ! Step size for first interval
         h2 = x(i + 2) - x(i + 1)     ! Step size for second interval

         f0 = y(i)
         f1 = y(i + 1)
         f2 = y(i + 2)

         ! Simpson's rule: (h/3) * (f0 + 4f1 + f2)
         sum = sum + (h1 + h2)/6.0_dp*(f0 + 4.0_dp*f1 + f2)
      end do

      ! Handle the case where n is odd (last interval)
      if (MOD(n, 2) == 0) then
         sum = sum + 0.5_dp*(x(n) - x(n - 1))*(y(n) + y(n - 1))
      end if

      result = sum
   end subroutine simpson_integration

   subroutine romberg_integration(x, y, result)
      real(dp), dimension(:), intent(in) :: x, y
      real(dp), intent(out) :: result

      integer :: i, j, k, n, m
      real(dp), dimension(:), allocatable :: R
      real(dp) :: h, sum, factor

      n = size(x)
      m = int(log(real(n, DP))/log(2.0_dp)) + 1  ! Number of refinement levels

      ! Validate input sizes
      if (size(x) /= size(y)) then
         print *, "Error: x and y arrays must have the same size."
         call mesa_error(__FILE__, __LINE__)
      end if

      if (n < 2) then
         print *, "Error: x and y arrays must have at least 2 elements."
         call mesa_error(__FILE__, __LINE__)
      end if

      allocate (R(m))

      ! Compute initial trapezoidal rule estimate
      h = x(n) - x(1)
      R(1) = 0.5_dp*h*(y(1) + y(n))

      ! Refinement using Romberg's method
      do j = 2, m
         sum = 0.0_dp
         do i = 1, 2**(j - 2)
            sum = sum + y(1 + (2*i - 1)*(n - 1)/(2**(j - 1)))
         end do

         h = h/2.0_dp
         R(j) = 0.5_dp*R(j - 1) + h*sum

         ! Richardson extrapolation
         factor = 4.0_dp
         do k = j, 2, -1
            R(k - 1) = (factor*R(k) - R(k - 1))/(factor - 1.0_dp)
            factor = factor*4.0_dp
         end do
      end do

      result = R(1)
   end subroutine romberg_integration

   !-----------------------------------------------------------------------
   ! File I/O functions
   !-----------------------------------------------------------------------

   !****************************
   ! Load Vega SED for Zero Point Calculation
   !****************************
   subroutine load_vega_sed(filepath, wavelengths, flux)
      character(len=*), intent(in) :: filepath
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, flux
      character(len=512) :: line
      integer :: unit, n_rows, status, i
      real(dp) :: temp_wave, temp_flux

      unit = 20
      open (unit, file=trim(filepath), status='OLD', action='READ', iostat=status)
      if (status /= 0) then
         print *, "Error: Could not open Vega SED file ", trim(filepath)
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Skip header line
      read (unit, '(A)', iostat=status) line
      if (status /= 0) then
         print *, "Error: Could not read header from Vega SED file ", trim(filepath)
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Count the number of data lines
      n_rows = 0
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         n_rows = n_rows + 1
      end do

      rewind (unit)
      read (unit, '(A)', iostat=status) line  ! Skip header again

      allocate (wavelengths(n_rows))
      allocate (flux(n_rows))

      i = 0
      do
         read (unit, *, iostat=status) temp_wave, temp_flux  ! Ignore any extra columns
         if (status /= 0) exit
         i = i + 1
         wavelengths(i) = temp_wave
         flux(i) = temp_flux
      end do

      close (unit)
   end subroutine load_vega_sed

   !****************************
   ! Load Filter File
   !****************************
   subroutine load_filter(directory, filter_wavelengths, filter_trans)
      character(len=*), intent(in) :: directory
      real(dp), dimension(:), allocatable, intent(out) :: filter_wavelengths, filter_trans

      character(len=512) :: line
      integer :: unit, n_rows, status, i
      real(dp) :: temp_wavelength, temp_trans

      ! Open the file
      unit = 20
      open (unit, file=trim(directory), status='OLD', action='READ', iostat=status)
      if (status /= 0) then
         print *, "Error: Could not open file ", trim(directory)
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Skip header line
      read (unit, '(A)', iostat=status) line
      if (status /= 0) then
         print *, "Error: Could not read the file", trim(directory)
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Count rows in the file
      n_rows = 0
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         n_rows = n_rows + 1
      end do

      ! Allocate arrays
      allocate (filter_wavelengths(n_rows))
      allocate (filter_trans(n_rows))

      ! Rewind to the first non-comment line
      rewind (unit)
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) then
            print *, "Error: Could not rewind file", trim(directory)
            call mesa_error(__FILE__, __LINE__)
         end if
         if (line(1:1) /= "#") exit
      end do

      ! Read and parse data
      i = 0
      do
         read (unit, *, iostat=status) temp_wavelength, temp_trans
         if (status /= 0) exit
         i = i + 1

         filter_wavelengths(i) = temp_wavelength
         filter_trans(i) = temp_trans
      end do

      close (unit)
   end subroutine load_filter

   !****************************
   ! Load Lookup Table For Identifying Stellar Atmosphere Models
   !****************************
   subroutine load_lookup_table(lookup_file, lookup_table, out_file_names, &
                                out_logg, out_meta, out_teff)

      character(len=*), intent(in) :: lookup_file
      REAL, dimension(:, :), allocatable, intent(out) :: lookup_table
      character(len=100), allocatable, intent(inout) :: out_file_names(:)
      real(dp), allocatable, intent(inout) :: out_logg(:), out_meta(:), out_teff(:)

      integer :: i, n_rows, status, unit
      character(len=512) :: line
      character(len=*), parameter :: delimiter = ","
      character(len=100), allocatable :: columns(:), headers(:)
      integer :: logg_col, meta_col, teff_col

      ! Open the file
      unit = 10
      open (unit, file=lookup_file, status='old', action='read', iostat=status)
      if (status /= 0) then
         print *, "Error: Could not open file", lookup_file
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Read header line
      read (unit, '(A)', iostat=status) line
      if (status /= 0) then
         print *, "Error: Could not read header line"
         call mesa_error(__FILE__, __LINE__)
      end if

      call split_line(line, delimiter, headers)

      ! Determine column indices for logg, meta, and teff
      logg_col = get_column_index(headers, "logg")
      teff_col = get_column_index(headers, "teff")

      meta_col = get_column_index(headers, "meta")
      if (meta_col < 0) then
         meta_col = get_column_index(headers, "feh")
      end if

      n_rows = 0
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         n_rows = n_rows + 1
      end do
      rewind (unit)

      ! Skip header
      read (unit, '(A)', iostat=status) line

      ! Allocate output arrays
      allocate (out_file_names(n_rows))
      allocate (out_logg(n_rows), out_meta(n_rows), out_teff(n_rows))

      ! Read and parse the file
      i = 0
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         i = i + 1

         call split_line(line, delimiter, columns)

         ! Populate arrays
         out_file_names(i) = columns(1)

         if (logg_col > 0) then
            if (columns(logg_col) /= "") then
               read (columns(logg_col), *) out_logg(i)
            else
               out_logg(i) = -999.0
            end if
         else
            out_logg(i) = -999.0
         end if

         if (meta_col > 0) then
            if (columns(meta_col) /= "") then
               read (columns(meta_col), *) out_meta(i)
            else
               out_meta(i) = 0.0
            end if
         else
            out_meta(i) = 0.0
         end if

         if (teff_col > 0) then
            if (columns(teff_col) /= "") then
               read (columns(teff_col), *) out_teff(i)
            else
               out_teff(i) = 0.0
            end if
         else
            out_teff(i) = 0.0
         end if
      end do

      close (unit)

   contains

      function get_column_index(headers, target) result(index)
         character(len=100), intent(in) :: headers(:)
         character(len=*), intent(in) :: target
         integer :: index, i
         character(len=100) :: clean_header, clean_target

         index = -1
         clean_target = trim(adjustl(target))  ! Clean the target string

         do i = 1, size(headers)
            clean_header = trim(adjustl(headers(i)))  ! Clean each header
            if (clean_header == clean_target) then
               index = i
               exit
            end if
         end do
      end function get_column_index

      subroutine split_line(line, delimiter, tokens)
         character(len=*), intent(in) :: line, delimiter
         character(len=100), allocatable, intent(out) :: tokens(:)
         integer :: num_tokens, pos, start, len_delim

         len_delim = len_trim(delimiter)
         start = 1
         num_tokens = 0
         if (allocated(tokens)) deallocate (tokens)

         do
            pos = index(line(start:), delimiter)

            if (pos == 0) exit
            num_tokens = num_tokens + 1
            call append_token(tokens, line(start:start + pos - 2))
            start = start + pos + len_delim - 1
         end do

         num_tokens = num_tokens + 1
         call append_token(tokens, line(start:))
      end subroutine split_line

      subroutine append_token(tokens, token)
         character(len=*), intent(in) :: token
         character(len=100), allocatable, intent(inout) :: tokens(:)
         character(len=100), allocatable :: temp(:)
         integer :: n

         if (.not. allocated(tokens)) then
            allocate (tokens(1))
            tokens(1) = token
         else
            n = size(tokens)
            allocate (temp(n))
            temp = tokens  ! Backup the current tokens
            deallocate (tokens)  ! Deallocate the old array
            allocate (tokens(n + 1))  ! Allocate with one extra space
            tokens(1:n) = temp  ! Restore old tokens
            tokens(n + 1) = token  ! Add the new token
            deallocate (temp)  ! unsure if this is till needed.
         end if
      end subroutine append_token

   end subroutine load_lookup_table

   subroutine load_sed(directory, index, wavelengths, flux)
      character(len=*), intent(in) :: directory
      integer, intent(in) :: index
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, flux

      character(len=512) :: line
      integer :: unit, n_rows, status, i
      real(dp) :: temp_wavelength, temp_flux

      ! Open the file
      unit = 20
      open (unit, file=trim(directory), status='OLD', action='READ', iostat=status)
      if (status /= 0) then
         print *, "Error: Could not open file ", trim(directory)
         call mesa_error(__FILE__, __LINE__)
      end if

      ! Skip header lines
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) then
            print *, "Error: Could not read the file", trim(directory)
            call mesa_error(__FILE__, __LINE__)
         end if
         if (line(1:1) /= "#") exit
      end do

      ! Count rows in the file
      n_rows = 0
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         n_rows = n_rows + 1
      end do

      ! Allocate arrays
      allocate (wavelengths(n_rows))
      allocate (flux(n_rows))

      ! Rewind to the first non-comment line
      rewind (unit)
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) then
            print *, "Error: Could not rewind file", trim(directory)
            call mesa_error(__FILE__, __LINE__)
         end if
         if (line(1:1) /= "#") exit
      end do

      ! Read and parse data
      i = 0
      do
         read (unit, *, iostat=status) temp_wavelength, temp_flux
         if (status /= 0) exit
         i = i + 1
         ! Convert f_lambda to f_nu
         wavelengths(i) = temp_wavelength
         flux(i) = temp_flux
      end do

      close (unit)

   end subroutine load_sed

   !-----------------------------------------------------------------------
   ! Helper function for file names
   !-----------------------------------------------------------------------

   function remove_dat(path) result(base)
      ! Extracts the portion of the string before the first dot
      character(len=*), intent(in) :: path
      character(len=strlen) :: base
      integer :: first_dot

      ! Find the position of the first dot
      first_dot = 0
      do while (first_dot < len_trim(path) .and. path(first_dot + 1:first_dot + 1) /= '.')
         first_dot = first_dot + 1
      end do

      ! Check if a dot was found
      if (first_dot < len_trim(path)) then
         ! Extract the part before the dot
         base = path(:first_dot)
      else
         ! No dot found, return the input string
         base = path
      end if
   end function remove_dat

   function basename(path) result(name)
      character(len=*), intent(in) :: path
      character(len=strlen) :: name
      integer :: i
      if (len_trim(path) == 0) then
         name = ''
         return
      end if
      i = index(path, '/', back=.true.)
      name = path(i + 1:)
   end function basename

   subroutine read_strings_from_file(colors_settings, strings, n, ierr)
      character(len=512) :: filename
      character(len=100), allocatable, intent(out) :: strings(:)
      integer, intent(out) :: n, ierr
      integer :: unit, i, status
      character(len=100) :: line
      type(Colors_General_Info), pointer :: colors_settings

      ierr = 0

      filename = trim(mesa_dir)//trim(colors_settings%instrument)//"/"// &
                 trim(basename(colors_settings%instrument))

      !filename = trim(mesa_dir)//trim(colors_settings%instrument)//"/"

      n = 0
      unit = 10
      open (unit, file=filename, status='old', action='read', iostat=status)
      if (status /= 0) then
         ierr = -1
         print *, "Error: Could not open file", filename
         call mesa_error(__FILE__, __LINE__)
      end if

      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         n = n + 1
      end do
      rewind (unit)
      if (allocated(strings)) deallocate (strings)
      allocate (strings(n))
      do i = 1, n
         read (unit, '(A)') strings(i)
      end do
      close (unit)
   end subroutine read_strings_from_file

end module colors_utils
