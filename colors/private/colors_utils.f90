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
   use colors_def, only: Colors_General_Info, sed_mem_cache_cap
   use utils_lib, only: mesa_error

   implicit none

   public :: dilute_flux, trapezoidal_integration, &
             simpson_integration, load_sed, load_filter, load_vega_sed, &
             load_lookup_table, remove_dat, load_flux_cube, build_unique_grids, &
             build_grid_to_lu_map, &
             find_containing_cell, find_interval, find_nearest_point, &
             find_bracket_index, load_sed_cached, load_stencil
contains

   ! apply dilution factor (R/d)^2 to convert surface flux to observed flux
   subroutine dilute_flux(surface_flux, R, d, calibrated_flux)
      real(dp), intent(in)  :: surface_flux(:)
      real(dp), intent(in)  :: R, d  ! R = stellar radius, d = distance (both in the same units, e.g., cm)
      real(dp), intent(out) :: calibrated_flux(:)

      if (size(calibrated_flux) /= size(surface_flux)) then
         print *, "Error in dilute_flux: Output array must have the same size as input array."
         call mesa_error(__FILE__, __LINE__)
      end if

      calibrated_flux = surface_flux*((R/d)**2)
   end subroutine dilute_flux

   subroutine trapezoidal_integration(x, y, result)
      real(dp), dimension(:), intent(in) :: x, y
      real(dp), intent(out) :: result

      integer :: i, n
      real(dp) :: sum

      n = size(x)
      sum = 0.0_dp

      if (size(x) /= size(y)) then
         print *, "Error: x and y arrays must have the same size."
         call mesa_error(__FILE__, __LINE__)
      end if

      if (size(x) < 2) then
         print *, "Error: x and y arrays must have at least 2 elements."
         call mesa_error(__FILE__, __LINE__)
      end if

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

      if (size(x) /= size(y)) then
         print *, "Error: x and y arrays must have the same size."
         call mesa_error(__FILE__, __LINE__)
      end if

      if (size(x) < 2) then
         print *, "Error: x and y arrays must have at least 2 elements."
         call mesa_error(__FILE__, __LINE__)
      end if

      ! adaptive Simpson's rule
      do i = 1, n - 2, 2
         h1 = x(i + 1) - x(i)
         h2 = x(i + 2) - x(i + 1)

         f0 = y(i)
         f1 = y(i + 1)
         f2 = y(i + 2)

         ! Simpson's rule: (h/3) * (f0 + 4f1 + f2)
         sum = sum + (h1 + h2)/6.0_dp * ( &
            f0*(2.0_dp*h1 - h2)/h1 + &
            f1*(h1 + h2)**2/(h1*h2) + &
            f2*(2.0_dp*h2 - h1)/h2)

      end do

      ! handle the last interval if n is even (odd number of points)
      if (MOD(n, 2) == 0) then
         sum = sum + 0.5_dp*(x(n) - x(n - 1))*(y(n) + y(n - 1))
      end if

      result = sum
   end subroutine simpson_integration

   subroutine load_vega_sed(filepath, wavelengths, flux)
      character(len=*), intent(in) :: filepath
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, flux
      character(len=512) :: line
      integer :: unit, n_rows, status, i
      real(dp) :: temp_wave, temp_flux

      open (newunit=unit, file=trim(filepath), status='OLD', action='READ', iostat=status)
      if (status /= 0) then
         print *, "Error: Could not open Vega SED file ", trim(filepath)
         call mesa_error(__FILE__, __LINE__)
      end if

      read (unit, '(A)', iostat=status) line
      if (status /= 0) then
         print *, "Error: Could not read header from Vega SED file ", trim(filepath)
         call mesa_error(__FILE__, __LINE__)
      end if

      n_rows = 0
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         n_rows = n_rows + 1
      end do

      rewind (unit)
      read (unit, '(A)', iostat=status) line  ! skip header again

      allocate (wavelengths(n_rows))
      allocate (flux(n_rows))

      i = 0
      do
         read (unit, *, iostat=status) temp_wave, temp_flux  ! ignore any extra columns
         if (status /= 0) exit
         i = i + 1
         wavelengths(i) = temp_wave
         flux(i) = temp_flux
      end do

      close (unit)
   end subroutine load_vega_sed

   subroutine load_filter(directory, filter_wavelengths, filter_trans)
      character(len=*), intent(in) :: directory
      real(dp), dimension(:), allocatable, intent(out) :: filter_wavelengths, filter_trans

      character(len=512) :: line
      integer :: unit, n_rows, status, i
      real(dp) :: temp_wavelength, temp_trans

      open (newunit=unit, file=trim(directory), status='OLD', action='READ', iostat=status)
      if (status /= 0) then
         print *, "Error: Could not open file ", trim(directory)
         call mesa_error(__FILE__, __LINE__)
      end if

      read (unit, '(A)', iostat=status) line
      if (status /= 0) then
         print *, "Error: Could not read the file", trim(directory)
         call mesa_error(__FILE__, __LINE__)
      end if

      n_rows = 0
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         n_rows = n_rows + 1
      end do

      allocate (filter_wavelengths(n_rows))
      allocate (filter_trans(n_rows))

      ! rewind to the first non-comment line
      rewind (unit)
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) then
            print *, "Error: Could not rewind file", trim(directory)
            call mesa_error(__FILE__, __LINE__)
         end if
         if (line(1:1) /= "#") exit
      end do

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

   ! parses a csv lookup table mapping atmosphere grid parameters to SED filenames
   subroutine load_lookup_table(lookup_file, lookup_table, out_file_names, &
                                out_logg, out_meta, out_teff)

      character(len=*), intent(in) :: lookup_file
      REAL, dimension(:, :), allocatable, intent(out) :: lookup_table
      character(len=100), allocatable, intent(inout) :: out_file_names(:)
      real(dp), allocatable, intent(inout) :: out_logg(:), out_meta(:), out_teff(:)

      integer :: i, n_rows, status, unit, ios
      character(len=512) :: line
      character(len=*), parameter :: delimiter = ","
      character(len=100), allocatable :: columns(:), headers(:)
      character(len=256) :: token
      integer :: logg_col, meta_col, teff_col

      open (newunit=unit, file=lookup_file, status='old', action='read', iostat=status)
      if (status /= 0) then
         print *, "Error: Could not open file", lookup_file
         call mesa_error(__FILE__, __LINE__)
      end if

      read (unit, '(A)', iostat=status) line
      if (status /= 0) then
         print *, "Error: Could not read header line"
         call mesa_error(__FILE__, __LINE__)
      end if

      call split_line(line, delimiter, headers)

      ! determine column indices -- try all plausible header name variants
      logg_col = get_column_index(headers, "logg")
      if (logg_col < 0) logg_col = get_column_index(headers, "log_g")
      if (logg_col < 0) logg_col = get_column_index(headers, "log(g)")
      if (logg_col < 0) logg_col = get_column_index(headers, "log10g")
      if (logg_col < 0) logg_col = get_column_index(headers, "log10_g")

      teff_col = get_column_index(headers, "teff")
      if (teff_col < 0) teff_col = get_column_index(headers, "t_eff")
      if (teff_col < 0) teff_col = get_column_index(headers, "t(eff)")
      if (teff_col < 0) teff_col = get_column_index(headers, "temperature")
      if (teff_col < 0) teff_col = get_column_index(headers, "temp")

      meta_col = get_column_index(headers, "meta")
      if (meta_col < 0) meta_col = get_column_index(headers, "feh")
      if (meta_col < 0) meta_col = get_column_index(headers, "fe_h")
      if (meta_col < 0) meta_col = get_column_index(headers, "[fe/h]")
      if (meta_col < 0) meta_col = get_column_index(headers, "mh")
      if (meta_col < 0) meta_col = get_column_index(headers, "[m/h]")
      if (meta_col < 0) meta_col = get_column_index(headers, "m_h")
      if (meta_col < 0) meta_col = get_column_index(headers, "z")
      if (meta_col < 0) meta_col = get_column_index(headers, "logz")
      if (meta_col < 0) meta_col = get_column_index(headers, "metallicity")
      if (meta_col < 0) meta_col = get_column_index(headers, "metals")

      n_rows = 0
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         n_rows = n_rows + 1
      end do
      rewind (unit)

      read (unit, '(A)', iostat=status) line  ! skip header

      allocate (out_file_names(n_rows))
      allocate (out_logg(n_rows), out_meta(n_rows), out_teff(n_rows))

      i = 0
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         i = i + 1

         call split_line(line, delimiter, columns)

         out_file_names(i) = columns(1)

         ! robust numeric parsing: never crash on bad/missing values
         if (logg_col > 0) then
            token = trim(adjustl(columns(logg_col)))
            if (len_trim(token) == 0 .or. token == '""') then
               out_logg(i) = -999.0_dp
            else
               read(token, *, iostat=ios) out_logg(i)
               if (ios /= 0) out_logg(i) = -999.0_dp
            end if
         else
            out_logg(i) = -999.0
         end if

         if (meta_col > 0) then
            token = trim(adjustl(columns(meta_col)))
            if (len_trim(token) == 0 .or. token == '""') then
               out_meta(i) = 0.0_dp
            else
               read(token, *, iostat=ios) out_meta(i)
               if (ios /= 0) out_meta(i) = 0.0_dp
            end if
         else
            out_meta(i) = 0.0
         end if

         if (teff_col > 0) then
            token = trim(adjustl(columns(teff_col)))
            if (len_trim(token) == 0 .or. token == '""') then
               out_teff(i) = 0.0_dp
            else
               read(token, *, iostat=ios) out_teff(i)
               if (ios /= 0) out_teff(i) = 0.0_dp
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
         clean_target = trim(adjustl(target))

         do i = 1, size(headers)
            clean_header = trim(adjustl(headers(i)))
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
            temp = tokens
            deallocate (tokens)
            allocate (tokens(n + 1))
            tokens(1:n) = temp
            tokens(n + 1) = token
            deallocate (temp)  ! unsure if this is till needed.
         end if
      end subroutine append_token

   end subroutine load_lookup_table




   subroutine load_sed(directory, index, wavelengths, flux)
      character(len=*), intent(in) :: directory
      integer, intent(in) :: index
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, flux

      character(len=512) :: line
      integer :: unit, n_rows, status, i, header_lines
      real(dp) :: temp_wavelength, temp_flux


      header_lines = 0
      open (newunit=unit, file=trim(directory), status='OLD', action='READ', iostat=status)
      if (status /= 0) then
         print *, "Error: Could not open file ", trim(directory)
         call mesa_error(__FILE__, __LINE__)
      end if

      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         if (line(1:1) /= "#") exit
         header_lines = header_lines + 1
      end do

      ! count data rows (we've already read one; count it)
      n_rows = 1
      do
         read (unit, '(A)', iostat=status) line
         if (status /= 0) exit
         n_rows = n_rows + 1
      end do

      allocate (wavelengths(n_rows))
      allocate (flux(n_rows))

      rewind (unit)
      ! skip exactly header_lines lines
      do i = 1, header_lines
         read (unit, '(A)', iostat=status) line
      end do

      i = 0
      do
         read (unit, *, iostat=status) temp_wavelength, temp_flux
         if (status /= 0) exit
         i = i + 1
         wavelengths(i) = temp_wavelength
         flux(i) = temp_flux
      end do

      close (unit)

   end subroutine load_sed

   function remove_dat(path) result(base)
      ! returns the portion of the string before the first dot
      character(len=*), intent(in) :: path
      character(len=strlen) :: base
      integer :: first_dot

      first_dot = 0
      do while (first_dot < len_trim(path) .and. path(first_dot + 1:first_dot + 1) /= '.')
         first_dot = first_dot + 1
      end do

      if (first_dot < len_trim(path)) then
         base = path(:first_dot)
      else
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

   function resolve_path(path) result(full_path)
      use const_def, only: mesa_dir
      character(len=*), intent(in) :: path
      character(len=512) :: full_path
      character(len=:), allocatable :: p
      logical :: exists
      integer :: n

      exists = .false.
      p = trim(adjustl(path))
      n = len_trim(p)

      if (n >= 2 .and. p(1:2) == './') then
         full_path = p
      else if (n >= 3 .and. p(1:3) == '../') then
         full_path = p

      else if (n >= 1 .and. p(1:1) == '/') then
         inquire (file=p, exist=exists)
         if (.not. exists) inquire (file=p//'/.', exist=exists)

         if (exists) then
            full_path = p
         else
            !it might be nice to warn the user that their filepath implies absolute
            !silently correcting this WILL lead to confusion
            !warning every step seems a smidge over the top though
            !write (*, *) trim(p), " not found. Trying ", trim(mesa_dir)//trim(p)
            full_path = trim(mesa_dir)//trim(p)
         end if

      else
         full_path = trim(mesa_dir)//'/'//trim(p)
      end if
   end function resolve_path

   subroutine read_strings_from_file(colors_settings, strings, n, ierr)
      character(len=512) :: filename
      character(len=100), allocatable, intent(out) :: strings(:)
      integer, intent(out) :: n, ierr
      integer :: unit, i, status
      character(len=100) :: line
      type(Colors_General_Info), pointer :: colors_settings

      ierr = 0

      filename = trim(resolve_path(colors_settings%instrument))//"/"// &
                 trim(basename(colors_settings%instrument))

      n = 0
      open (newunit=unit, file=filename, status='old', action='read', iostat=status)
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

   ! load flux cube from binary file into handle at initialization.
   ! if the file cannot be opened or the large flux_cube array cannot be
   ! allocated, sets cube_loaded = .false. so the runtime will fall back to
   ! loading individual SED files via the lookup table.
   ! grids and wavelengths are always loaded (small); only the 4-D cube
   ! allocation is treated as the fallback trigger.
   subroutine load_flux_cube(rq, stellar_model_dir)
      type(Colors_General_Info), intent(inout) :: rq
      character(len=*), intent(in) :: stellar_model_dir

      character(len=512) :: bin_filename
      integer :: unit, status, n_teff, n_logg, n_meta, n_lambda
      real(dp) :: cube_mb

      rq%cube_loaded = .false.

      bin_filename = trim(resolve_path(stellar_model_dir))//'/flux_cube.bin'

      open (newunit=unit, file=trim(bin_filename), status='OLD', &
            access='STREAM', form='UNFORMATTED', iostat=status)
      if (status /= 0) then
         ! no binary cube available -- will use individual SED files
         write (*, '(a)') 'colors: no flux_cube.bin found; using per-file SED loading'
         return
      end if

      read (unit, iostat=status) n_teff, n_logg, n_meta, n_lambda
      if (status /= 0) then
         close (unit)
         return
      end if

      ! attempt the large allocation first -- this is the one that may fail.
      ! doing it before the small grid allocations avoids partial cleanup.
      allocate (rq%cube_flux(n_teff, n_logg, n_meta, n_lambda), stat=status)
      if (status /= 0) then
         if (allocated(rq%cube_flux)) deallocate (rq%cube_flux)
         cube_mb = real(n_teff, dp)*n_logg*n_meta*n_lambda*8.0_dp/(1024.0_dp**2)
         write (*, '(a,f0.1,a)') &
            'colors: flux cube allocation failed (', cube_mb, &
            ' MB); falling back to per-file SED loading'
         close (unit)
         return
      end if

      ! grid arrays are small -- always expected to succeed
      allocate (rq%cube_teff_grid(n_teff), stat=status)
      if (status /= 0) goto 900

      allocate (rq%cube_logg_grid(n_logg), stat=status)
      if (status /= 0) goto 900

      allocate (rq%cube_meta_grid(n_meta), stat=status)
      if (status /= 0) goto 900

      allocate (rq%cube_wavelengths(n_lambda), stat=status)
      if (status /= 0) goto 900

      read (unit, iostat=status) rq%cube_teff_grid
      if (status /= 0) goto 900

      read (unit, iostat=status) rq%cube_logg_grid
      if (status /= 0) goto 900

      read (unit, iostat=status) rq%cube_meta_grid
      if (status /= 0) goto 900

      read (unit, iostat=status) rq%cube_wavelengths
      if (status /= 0) goto 900

      read (unit, iostat=status) rq%cube_flux
      if (status /= 0) goto 900

      close (unit)
      rq%cube_loaded = .true.

      cube_mb = real(n_teff, dp)*n_logg*n_meta*n_lambda*8.0_dp/(1024.0_dp**2)
      write (*, '(a,i0,a,i0,a,i0,a,i0,a,f0.1,a)') &
         'colors: flux cube loaded (', &
         n_teff, ' x ', n_logg, ' x ', n_meta, ' x ', n_lambda, &
         ', ', cube_mb, ' MB)'
      return

      ! error cleanup -- deallocate everything that may have been allocated
900   continue
      write (*, '(a)') 'colors: error reading flux_cube.bin; falling back to per-file SED loading'
      if (allocated(rq%cube_flux)) deallocate (rq%cube_flux)
      if (allocated(rq%cube_teff_grid)) deallocate (rq%cube_teff_grid)
      if (allocated(rq%cube_logg_grid)) deallocate (rq%cube_logg_grid)
      if (allocated(rq%cube_meta_grid)) deallocate (rq%cube_meta_grid)
      if (allocated(rq%cube_wavelengths)) deallocate (rq%cube_wavelengths)
      close (unit)

   end subroutine load_flux_cube

   ! build unique sorted grids from lookup table and store on handle.
   ! called once at init so the fallback interpolation path never rebuilds these.
   subroutine build_unique_grids(rq)
      type(Colors_General_Info), intent(inout) :: rq
      logical :: found

      if (rq%unique_grids_built) return
      if (.not. rq%lookup_loaded) return

      call extract_unique_sorted(rq%lu_teff, rq%u_teff)
      call extract_unique_sorted(rq%lu_logg, rq%u_logg)
      call extract_unique_sorted(rq%lu_meta, rq%u_meta)

      rq%unique_grids_built = .true.

   contains

      subroutine extract_unique_sorted(arr, unique)
         real(dp), intent(in) :: arr(:)
         real(dp), allocatable, intent(out) :: unique(:)
         real(dp), allocatable :: buf(:)
         integer :: ii, jj, nn, nnu
         real(dp) :: sw

         nn = size(arr)
         allocate (buf(nn))
         nnu = 0

         do ii = 1, nn
            found = .false.
            do jj = 1, nnu
               if (abs(arr(ii) - buf(jj)) < 1.0e-10_dp) then
                  found = .true.
                  exit
               end if
            end do
            if (.not. found) then
               nnu = nnu + 1
               buf(nnu) = arr(ii)
            end if
         end do

         ! insertion sort (grids are small)
         do ii = 2, nnu
            sw = buf(ii)
            jj = ii - 1
            do while (jj >= 1 .and. buf(jj) > sw)
               buf(jj + 1) = buf(jj)
               jj = jj - 1
            end do
            buf(jj + 1) = sw
         end do

         allocate (unique(nnu))
         unique = buf(1:nnu)
         deallocate (buf)
      end subroutine extract_unique_sorted

   end subroutine build_unique_grids

   ! build a 3-D mapping from unique-grid indices to lookup-table rows.
   ! grid_to_lu(i_t, i_g, i_m) = row in lu_* whose parameters match
   ! (u_teff(i_t), u_logg(i_g), u_meta(i_m)).
   ! called once at init; replaces the O(n_lu) nearest-neighbour scan
   ! that previously ran per corner at runtime.
   subroutine build_grid_to_lu_map(rq)
      type(Colors_General_Info), intent(inout) :: rq

      integer :: nt, ng, nm, n_lu
      integer :: it, ig, im, idx
      integer :: best_idx
      real(dp) :: best_dist, dist
      real(dp), parameter :: tol = 1.0e-10_dp

      if (rq%grid_map_built) return
      if (.not. rq%unique_grids_built) return
      if (.not. rq%lookup_loaded) return

      nt = size(rq%u_teff)
      ng = size(rq%u_logg)
      nm = size(rq%u_meta)
      n_lu = size(rq%lu_teff)

      allocate (rq%grid_to_lu(nt, ng, nm))
      rq%grid_to_lu = 0

      do it = 1, nt
         do ig = 1, ng
            do im = 1, nm
               best_idx = 1
               best_dist = huge(1.0_dp)
               do idx = 1, n_lu
                  dist = abs(rq%lu_teff(idx) - rq%u_teff(it)) + &
                         abs(rq%lu_logg(idx) - rq%u_logg(ig)) + &
                         abs(rq%lu_meta(idx) - rq%u_meta(im))
                  if (dist < best_dist) then
                     best_dist = dist
                     best_idx = idx
                  end if
                  if (best_dist < tol) exit  ! exact match found
               end do
               rq%grid_to_lu(it, ig, im) = best_idx
            end do
         end do
      end do

      rq%grid_map_built = .true.

   end subroutine build_grid_to_lu_map

   subroutine find_containing_cell(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                   i_x, i_y, i_z, t_x, t_y, t_z)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      integer, intent(out) :: i_x, i_y, i_z
      real(dp), intent(out) :: t_x, t_y, t_z

      call find_interval(x_grid, x_val, i_x, t_x)
      call find_interval(y_grid, y_val, i_y, t_y)
      call find_interval(z_grid, z_val, i_z, t_z)
   end subroutine find_containing_cell

   ! find the interval in a sorted array containing a value.
   ! returns index i such that x(i) <= val <= x(i+1) and fractional position t in [0,1].
   ! detects dummy axes (all zeros / 999 / -999) and collapses them to i=1, t=0.
   subroutine find_interval(x, val, i, t)
      real(dp), intent(in) :: x(:), val
      integer, intent(out) :: i
      real(dp), intent(out) :: t

      integer :: n, lo, hi, mid
      logical :: dummy_axis

      n = size(x)

      ! detect dummy axis: all values == 0, 999, or -999
      dummy_axis = all(x == 0.0_dp) .or. all(x == 999.0_dp) .or. all(x == -999.0_dp)

      if (dummy_axis) then
         i = 1
         t = 0.0_dp
         return
      end if

      if (val <= x(1)) then
         i = 1
         t = 0.0_dp
         return
      else if (val >= x(n)) then
         i = n - 1
         t = 1.0_dp
         return
      end if

      lo = 1
      hi = n
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (val >= x(mid)) then
            lo = mid
         else
            hi = mid
         end if
      end do

      i = lo
      if (abs(x(i + 1) - x(i)) < 1.0e-30_dp) then
         t = 0.0_dp  ! degenerate interval -- no interpolation needed
      else
         t = (val - x(i))/(x(i + 1) - x(i))
      end if
   end subroutine find_interval

   subroutine find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      integer, intent(out) :: i_x, i_y, i_z

      i_x = minloc(abs(x_val - x_grid), 1)
      i_y = minloc(abs(y_val - y_grid), 1)
      i_z = minloc(abs(z_val - z_grid), 1)
   end subroutine find_nearest_point

   ! returns idx such that grid(idx) <= val < grid(idx+1), clamped to bounds
   subroutine find_bracket_index(grid, val, idx)
      real(dp), intent(in) :: grid(:), val
      integer, intent(out) :: idx

      integer :: n, lo, hi, mid

      n = size(grid)
      if (n < 2) then
         idx = 1
         return
      end if

      if (val <= grid(1)) then
         idx = 1
         return
      else if (val >= grid(n)) then
         idx = n - 1
         return
      end if

      lo = 1
      hi = n
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (val >= grid(mid)) then
            lo = mid
         else
            hi = mid
         end if
      end do
      idx = lo
   end subroutine find_bracket_index

   ! fallback SED cache and stencil loader

   ! load a stencil sub-cube of SED fluxes for the given index ranges.
   ! uses load_sed_cached so repeated visits to the same grid point are
   ! served from memory rather than disk.
   subroutine load_stencil(rq, resolved_dir, lo_t, hi_t, lo_g, hi_g, lo_m, hi_m)
      type(Colors_General_Info), intent(inout) :: rq
      character(len=*), intent(in) :: resolved_dir
      integer, intent(in) :: lo_t, hi_t, lo_g, hi_g, lo_m, hi_m

      integer :: st, sg, sm, n_lambda, lu_idx
      integer :: it, ig, im
      real(dp), dimension(:), allocatable :: sed_flux

      st = hi_t - lo_t + 1
      sg = hi_g - lo_g + 1
      sm = hi_m - lo_m + 1

      ! free previous stencil flux data
      if (allocated(rq%stencil_fluxes)) deallocate (rq%stencil_fluxes)

      n_lambda = 0

      do it = lo_t, hi_t
         do ig = lo_g, hi_g
            do im = lo_m, hi_m
               lu_idx = rq%grid_to_lu(it, ig, im)
               call load_sed_cached(rq, resolved_dir, lu_idx, sed_flux)

               if (n_lambda == 0) then
                  n_lambda = size(sed_flux)
                  allocate (rq%stencil_fluxes(st, sg, sm, n_lambda))
               else if (size(sed_flux) /= n_lambda) then
                  ! SED files have inconsistent wavelength counts — this is a
                  ! fatal grid inconsistency; crash with a clear message.
                  write (*, '(a,i0,a,i0,a,i0)') &
                     'colors ERROR: SED at lu_idx=', lu_idx, &
                     ' has ', size(sed_flux), ' wavelength points; expected ', n_lambda
                  write (*, '(a)') &
                     'colors: bt-settl (or other) grid files have non-uniform wavelength grids.'
                  call mesa_error(__FILE__, __LINE__)
               end if

               rq%stencil_fluxes(it - lo_t + 1, ig - lo_g + 1, im - lo_m + 1, :) = &
                  sed_flux(1:n_lambda)
               if (allocated(sed_flux)) deallocate (sed_flux)
            end do
         end do
      end do

      ! set stencil wavelengths from the canonical copy on the handle
      if (allocated(rq%stencil_wavelengths)) deallocate (rq%stencil_wavelengths)
      allocate (rq%stencil_wavelengths(n_lambda))
      rq%stencil_wavelengths = rq%fallback_wavelengths(1:n_lambda)

   end subroutine load_stencil

   ! retrieve an SED flux from the memory cache, or load from disk on miss.
   ! uses a bounded circular buffer (sed_mem_cache_cap slots).
   ! on the first disk read, the wavelength array is stored once on the handle
   ! as rq%fallback_wavelengths -- all SEDs in a given atmosphere grid share
   ! the same wavelengths, so only flux is cached and returned.
   subroutine load_sed_cached(rq, resolved_dir, lu_idx, flux)
      type(Colors_General_Info), intent(inout) :: rq
      character(len=*), intent(in) :: resolved_dir
      integer, intent(in) :: lu_idx
      real(dp), dimension(:), allocatable, intent(out) :: flux

      integer :: slot, n_lam, status
      character(len=512) :: filepath
      real(dp), dimension(:), allocatable :: sed_wave
      real(dp), dimension(:), allocatable :: flux_interp

      ! initialise the cache on first call
      if (.not. rq%sed_mcache_init) then
         allocate (rq%sed_mcache_keys(sed_mem_cache_cap))
         rq%sed_mcache_keys = 0   ! 0 means empty slot
         rq%sed_mcache_count = 0
         rq%sed_mcache_next = 1
         rq%sed_mcache_nlam = 0
         rq%sed_mcache_init = .true.
      end if

      ! search for a cache hit (linear scan over a small array)
      do slot = 1, rq%sed_mcache_count
         if (rq%sed_mcache_keys(slot) == lu_idx) then
            ! hit -- return cached flux
            n_lam = rq%sed_mcache_nlam
            allocate (flux(n_lam))
            flux = rq%sed_mcache_data(:, slot)
            return
         end if
      end do

      ! miss -- load from disk
      filepath = trim(resolved_dir)//'/'//trim(rq%lu_file_names(lu_idx))
      call load_sed(filepath, lu_idx, sed_wave, flux)

      ! store the canonical wavelength array on the handle (once only)
      if (.not. rq%fallback_wavelengths_set) then
         n_lam = size(sed_wave)
         allocate (rq%fallback_wavelengths(n_lam))
         rq%fallback_wavelengths = sed_wave
         rq%fallback_wavelengths_set = .true.
      end if



      ! store flux in the cache
      n_lam = size(flux)
      if (rq%sed_mcache_nlam == 0) then
         rq%sed_mcache_nlam = n_lam
         allocate (rq%sed_mcache_data(n_lam, sed_mem_cache_cap), stat=status)
         if (status /= 0) then
            write (*, '(a,f0.1,a)') 'colors: SED memory cache allocation failed (', &
               real(n_lam, dp)*sed_mem_cache_cap*8.0_dp/1024.0_dp**2, ' MB); cache disabled'
            rq%sed_mcache_init = .false.   ! fall through to disk-only path
            if (allocated(sed_wave)) deallocate (sed_wave)
            return
         end if
      else if (n_lam /= rq%sed_mcache_nlam) then
         ! non-uniform wavelength grids across SED files: interpolate onto the canonical grid
         if (.not. rq%fallback_wavelengths_set) then
            allocate (rq%fallback_wavelengths(n_lam))
            rq%fallback_wavelengths = sed_wave
            rq%fallback_wavelengths_set = .true.
         end if

         allocate (flux_interp(rq%sed_mcache_nlam), stat=status)
         if (status /= 0) then
            write (*, '(a)') 'colors: WARNING: could not allocate interpolation buffer; SED cache disabled'
            rq%sed_mcache_init = .false.
            if (allocated(sed_wave)) deallocate (sed_wave)
            return
         end if

         call interp_linear_internal(rq%fallback_wavelengths, sed_wave, flux, flux_interp)

         deallocate (flux)
         allocate (flux(rq%sed_mcache_nlam))
         flux = flux_interp
         deallocate (flux_interp)

         n_lam = rq%sed_mcache_nlam
      end if

      ! write to the next slot (circular)
      slot = rq%sed_mcache_next
      rq%sed_mcache_keys(slot) = lu_idx
      rq%sed_mcache_data(:, slot) = flux(1:n_lam)

      ! advance the circular pointer
      if (rq%sed_mcache_count < sed_mem_cache_cap) &
         rq%sed_mcache_count = rq%sed_mcache_count + 1
      rq%sed_mcache_next = mod(slot, sed_mem_cache_cap) + 1

      if (allocated(sed_wave)) deallocate (sed_wave)
   end subroutine load_sed_cached


   ! simple 1D linear interpolation (np.interp-like): clamps to endpoints
   subroutine interp_linear_internal(x_out, x_in, y_in, y_out)
      real(dp), intent(in) :: x_out(:), x_in(:), y_in(:)
      real(dp), intent(out) :: y_out(:)
      integer :: i, n_in, n_out, lo, hi, mid
      real(dp) :: t, denom

      n_in = size(x_in)
      n_out = size(x_out)

      if (n_out <= 0) return

      if (n_in <= 0) then
         y_out = 0.0_dp
         return
      end if

      if (n_in == 1) then
         y_out = y_in(1)
         return
      end if

      do i = 1, n_out
         if (x_out(i) <= x_in(1)) then
            y_out(i) = y_in(1)
         else if (x_out(i) >= x_in(n_in)) then
            y_out(i) = y_in(n_in)
         else
            lo = 1
            hi = n_in
            do while (hi - lo > 1)
               mid = (lo + hi)/2
               if (x_out(i) >= x_in(mid)) then
                  lo = mid
               else
                  hi = mid
               end if
            end do

            denom = x_in(lo + 1) - x_in(lo)
            if (abs(denom) <= 0.0_dp) then
               y_out(i) = y_in(lo)
            else
               t = (x_out(i) - x_in(lo))/denom
               y_out(i) = y_in(lo) + t*(y_in(lo + 1) - y_in(lo))
            end if
         end if
      end do

   end subroutine interp_linear_internal
end module colors_utils