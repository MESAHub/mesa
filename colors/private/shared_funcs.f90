MODULE shared_funcs
   USE const_def, ONLY: dp, strlen
   USE colors_def, ONLY: Colors_General_Info, sed_mem_cache_cap
   USE colors_utils, ONLY: load_sed
   USE utils_lib, ONLY: mesa_error
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: dilute_flux, trapezoidalintegration, rombergintegration, SimpsonIntegration, loadsed, &
             loadfilter, loadvegased, load_lookuptable, remove_dat, &
             find_containing_cell, find_interval, find_nearest_point, find_bracket_index, &
             load_sed_cached, load_stencil

CONTAINS

   ! ====================================================================
   !  Grid-search utilities (shared by hermite, linear, and knn modules)
   ! ====================================================================

   !---------------------------------------------------------------------------
   ! Find the cell containing the interpolation point
   !---------------------------------------------------------------------------
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

   !---------------------------------------------------------------------------
   ! Find the interval in a sorted array containing a value.
   ! Returns index i such that x(i) <= val < x(i+1), and the fractional
   ! position t in [0,1].  Dummy-axis detection collapses degenerate axes.
   !---------------------------------------------------------------------------
   subroutine find_interval(x, val, i, t)
      real(dp), intent(in) :: x(:), val
      integer, intent(out) :: i
      real(dp), intent(out) :: t

      integer :: n, lo, hi, mid
      logical :: dummy_axis

      n = size(x)

      ! Detect dummy axis: all values == 0, 999, or -999
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

   !---------------------------------------------------------------------------
   ! Find the nearest grid point (3-D)
   !---------------------------------------------------------------------------
   subroutine find_nearest_point(x_val, y_val, z_val, x_grid, y_grid, z_grid, &
                                 i_x, i_y, i_z)
      real(dp), intent(in) :: x_val, y_val, z_val
      real(dp), intent(in) :: x_grid(:), y_grid(:), z_grid(:)
      integer, intent(out) :: i_x, i_y, i_z

      i_x = minloc(abs(x_val - x_grid), 1)
      i_y = minloc(abs(y_val - y_grid), 1)
      i_z = minloc(abs(z_val - z_grid), 1)
   end subroutine find_nearest_point

   !---------------------------------------------------------------------------
   ! Find the bracketing index i such that grid(i) <= val < grid(i+1)
   !---------------------------------------------------------------------------
   subroutine find_bracket_index(grid, val, idx)
      real(dp), intent(in) :: grid(:), val
      integer, intent(out) :: idx
      integer :: n

      n = size(grid)
      if (n < 2) then
         idx = 1
         return
      end if

      if (val <= grid(1)) then
         idx = 1
      else if (val >= grid(n)) then
         idx = n - 1
      else
         idx = 1
         do while (idx < n - 1 .and. grid(idx + 1) <= val)
            idx = idx + 1
         end do
      end if
   end subroutine find_bracket_index

   ! ====================================================================
   !  Fallback-path shared infrastructure
   !  (stencil loading + SED memory cache, used by all interp modules)
   ! ====================================================================

   !---------------------------------------------------------------------------
   ! Load SEDs for every point in a stencil sub-cube, using the memory cache.
   ! Populates rq%stencil_fluxes and rq%stencil_wavelengths.
   !
   ! Wavelengths are stored once on the handle (rq%fallback_wavelengths)
   ! on the first disk read and reused thereafter -- all SEDs in a given
   ! atmosphere grid share the same wavelength array.
   !---------------------------------------------------------------------------
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

      ! Free previous stencil flux data
      if (allocated(rq%stencil_fluxes)) deallocate(rq%stencil_fluxes)

      n_lambda = 0

      do it = lo_t, hi_t
         do ig = lo_g, hi_g
            do im = lo_m, hi_m
               lu_idx = rq%grid_to_lu(it, ig, im)

               call load_sed_cached(rq, resolved_dir, lu_idx, sed_flux)

               if (n_lambda == 0) then
                  n_lambda = size(sed_flux)
                  allocate (rq%stencil_fluxes(st, sg, sm, n_lambda))
               end if

               rq%stencil_fluxes(it - lo_t + 1, ig - lo_g + 1, im - lo_m + 1, :) = &
                  sed_flux(1:n_lambda)

               if (allocated(sed_flux)) deallocate(sed_flux)
            end do
         end do
      end do

      ! Set stencil wavelengths from the canonical copy on the handle
      if (allocated(rq%stencil_wavelengths)) deallocate(rq%stencil_wavelengths)
      allocate (rq%stencil_wavelengths(n_lambda))
      rq%stencil_wavelengths = rq%fallback_wavelengths(1:n_lambda)

   end subroutine load_stencil

   !---------------------------------------------------------------------------
   ! Retrieve an SED flux from the memory cache, or load from disk on miss.
   ! Uses a bounded circular buffer (sed_mem_cache_cap slots).
   !
   ! On the first disk read, the wavelength array is stored once on the
   ! handle as rq%fallback_wavelengths (all SEDs in a given atmosphere
   ! grid share the same wavelength array).  Only the flux is cached
   ! and returned.
   !---------------------------------------------------------------------------
   subroutine load_sed_cached(rq, resolved_dir, lu_idx, flux)
      type(Colors_General_Info), intent(inout) :: rq
      character(len=*), intent(in) :: resolved_dir
      integer, intent(in) :: lu_idx
      real(dp), dimension(:), allocatable, intent(out) :: flux

      integer :: slot, n_lam
      character(len=512) :: filepath
      real(dp), dimension(:), allocatable :: sed_wave

      ! Initialise the cache on first call
      if (.not. rq%sed_mcache_init) then
         allocate (rq%sed_mcache_keys(sed_mem_cache_cap))
         rq%sed_mcache_keys = 0   ! 0 means empty slot
         rq%sed_mcache_count = 0
         rq%sed_mcache_next = 1
         rq%sed_mcache_nlam = 0
         rq%sed_mcache_init = .true.
      end if

      ! Search for a cache hit (linear scan over a small array)
      do slot = 1, rq%sed_mcache_count
         if (rq%sed_mcache_keys(slot) == lu_idx) then
            ! Hit -- return cached flux
            n_lam = rq%sed_mcache_nlam
            allocate (flux(n_lam))
            flux = rq%sed_mcache_data(:, slot)
            return
         end if
      end do

      ! Miss -- load from disk
      filepath = trim(resolved_dir)//'/'//trim(rq%lu_file_names(lu_idx))
      call load_sed(filepath, lu_idx, sed_wave, flux)

      ! Store the canonical wavelength array on the handle (once only)
      if (.not. rq%fallback_wavelengths_set) then
         n_lam = size(sed_wave)
         allocate (rq%fallback_wavelengths(n_lam))
         rq%fallback_wavelengths = sed_wave
         rq%fallback_wavelengths_set = .true.
      end if
      if (allocated(sed_wave)) deallocate(sed_wave)

      ! Store flux in the cache
      n_lam = size(flux)
      if (rq%sed_mcache_nlam == 0) then
         ! First SED ever loaded -- set the wavelength count and allocate data
         rq%sed_mcache_nlam = n_lam
         allocate (rq%sed_mcache_data(n_lam, sed_mem_cache_cap))
      end if

      ! Write to the next slot (circular)
      slot = rq%sed_mcache_next
      rq%sed_mcache_keys(slot) = lu_idx
      rq%sed_mcache_data(:, slot) = flux(1:rq%sed_mcache_nlam)

      if (rq%sed_mcache_count < sed_mem_cache_cap) then
         rq%sed_mcache_count = rq%sed_mcache_count + 1
      end if
      rq%sed_mcache_next = mod(slot, sed_mem_cache_cap) + 1

   end subroutine load_sed_cached

   ! ====================================================================
   !  Flux dilution
   ! ====================================================================

   SUBROUTINE dilute_flux(surface_flux, R, d, calibrated_flux)
      REAL(dp), INTENT(IN)  :: surface_flux(:)
      REAL(dp), INTENT(IN)  :: R, d
      REAL(dp), INTENT(OUT) :: calibrated_flux(:)

      IF (SIZE(calibrated_flux) /= SIZE(surface_flux)) THEN
         PRINT *, "Error in dilute_flux: Output array must have the same size as input array."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      calibrated_flux = surface_flux*((R/d)**2)
   END SUBROUTINE dilute_flux

   ! ====================================================================
   !  Integration routines
   ! ====================================================================

   SUBROUTINE trapezoidalintegration(x, y, result)
      REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
      REAL(DP), INTENT(OUT) :: result
      INTEGER :: i, n
      REAL(DP) :: sum

      n = SIZE(x)
      sum = 0.0_dp

      IF (SIZE(x) /= SIZE(y)) THEN
         PRINT *, "Error: x and y arrays must have the same size."
         CALL mesa_error(__FILE__, __LINE__)
      END IF
      IF (SIZE(x) < 2) THEN
         PRINT *, "Error: x and y arrays must have at least 2 elements."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      DO i = 1, n - 1
         sum = sum + 0.5_dp*(x(i + 1) - x(i))*(y(i + 1) + y(i))
      END DO
      result = sum
   END SUBROUTINE trapezoidalintegration

   SUBROUTINE SimpsonIntegration(x, y, result)
      REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
      REAL(DP), INTENT(OUT) :: result
      INTEGER :: i, n
      REAL(DP) :: sum, h1, h2, f0, f1, f2

      n = SIZE(x)
      sum = 0.0_DP

      IF (SIZE(x) /= SIZE(y)) THEN
         PRINT *, "Error: x and y arrays must have the same size."
         CALL mesa_error(__FILE__, __LINE__)
      END IF
      IF (SIZE(x) < 2) THEN
         PRINT *, "Error: x and y arrays must have at least 2 elements."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      DO i = 1, n - 2, 2
         h1 = x(i + 1) - x(i)
         h2 = x(i + 2) - x(i + 1)
         f0 = y(i); f1 = y(i + 1); f2 = y(i + 2)
         sum = sum + (h1 + h2)/6.0_DP*(f0 + 4.0_DP*f1 + f2)
      END DO
      IF (MOD(n, 2) == 0) THEN
         sum = sum + 0.5_DP*(x(n) - x(n - 1))*(y(n) + y(n - 1))
      END IF
      result = sum
   END SUBROUTINE SimpsonIntegration

   SUBROUTINE rombergintegration(x, y, result)
      REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
      REAL(DP), INTENT(OUT) :: result
      INTEGER :: i, j, k, n, m
      REAL(DP), DIMENSION(:), ALLOCATABLE :: R
      REAL(DP) :: h, sum, factor

      n = SIZE(x)
      m = INT(LOG(REAL(n, DP))/LOG(2.0_DP)) + 1

      IF (SIZE(x) /= SIZE(y)) THEN
         PRINT *, "Error: x and y arrays must have the same size."
         CALL mesa_error(__FILE__, __LINE__)
      END IF
      IF (n < 2) THEN
         PRINT *, "Error: x and y arrays must have at least 2 elements."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ALLOCATE (R(m))
      h = x(n) - x(1)
      R(1) = 0.5_DP*h*(y(1) + y(n))

      DO j = 2, m
         sum = 0.0_DP
         DO i = 1, 2**(j - 2)
            sum = sum + y(1 + (2*i - 1)*(n - 1)/(2**(j - 1)))
         END DO
         h = h/2.0_DP
         R(j) = 0.5_DP*R(j - 1) + h*sum
         factor = 4.0_DP
         DO k = j, 2, -1
            R(k - 1) = (factor*R(k) - R(k - 1))/(factor - 1.0_DP)
            factor = factor*4.0_DP
         END DO
      END DO
      result = R(1)
      DEALLOCATE (R)
   END SUBROUTINE rombergintegration

   ! ====================================================================
   !  Legacy file I/O (retained for backward compatibility)
   ! ====================================================================

   SUBROUTINE loadvegased(filepath, wavelengths, flux)
      CHARACTER(LEN=*), INTENT(IN) :: filepath
      REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, flux
      CHARACTER(LEN=512) :: line
      INTEGER :: unit, n_rows, status, i
      REAL(dp) :: temp_wave, temp_flux

      OPEN (newunit=unit, FILE=TRIM(filepath), STATUS='OLD', ACTION='READ', IOSTAT=status)
      IF (status /= 0) THEN
         PRINT *, "Error: Could not open Vega SED file ", TRIM(filepath)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      READ (unit, '(A)', IOSTAT=status) line
      IF (status /= 0) THEN
         PRINT *, "Error: Could not read header from Vega SED file ", TRIM(filepath)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      n_rows = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         n_rows = n_rows + 1
      END DO

      REWIND (unit)
      READ (unit, '(A)', IOSTAT=status) line

      ALLOCATE (wavelengths(n_rows))
      ALLOCATE (flux(n_rows))

      i = 0
      DO
         READ (unit, *, IOSTAT=status) temp_wave, temp_flux
         IF (status /= 0) EXIT
         i = i + 1
         wavelengths(i) = temp_wave
         flux(i) = temp_flux
      END DO
      CLOSE (unit)
   END SUBROUTINE loadvegased

   SUBROUTINE loadfilter(directory, filter_wavelengths, filter_trans)
      CHARACTER(LEN=*), INTENT(IN) :: directory
      REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: filter_wavelengths, filter_trans
      CHARACTER(LEN=512) :: line
      INTEGER :: unit, n_rows, status, i
      REAL(dp) :: temp_wavelength, temp_trans

      OPEN (newunit=unit, FILE=TRIM(directory), STATUS='OLD', ACTION='READ', IOSTAT=status)
      IF (status /= 0) THEN
         PRINT *, "Error: Could not open file ", TRIM(directory)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      READ (unit, '(A)', IOSTAT=status) line
      IF (status /= 0) THEN
         PRINT *, "Error: Could not read the file", TRIM(directory)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      n_rows = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         n_rows = n_rows + 1
      END DO

      ALLOCATE (filter_wavelengths(n_rows))
      ALLOCATE (filter_trans(n_rows))

      REWIND (unit)
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) THEN
            PRINT *, "Error: Could not rewind file", TRIM(directory)
            CALL mesa_error(__FILE__, __LINE__)
         END IF
         IF (line(1:1) /= "#") EXIT
      END DO

      i = 0
      DO
         READ (unit, *, IOSTAT=status) temp_wavelength, temp_trans
         IF (status /= 0) EXIT
         i = i + 1
         filter_wavelengths(i) = temp_wavelength
         filter_trans(i) = temp_trans
      END DO
      CLOSE (unit)
   END SUBROUTINE loadfilter

   SUBROUTINE load_lookuptable(lookup_file, lookup_table, out_file_names, out_logg, out_meta, out_teff)
      CHARACTER(LEN=*), INTENT(IN) :: lookup_file
      REAL, DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: lookup_table
      CHARACTER(LEN=100), ALLOCATABLE, INTENT(INOUT) :: out_file_names(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: out_logg(:), out_meta(:), out_teff(:)

      INTEGER :: i, n_rows, status, unit
      CHARACTER(LEN=512) :: line
      CHARACTER(LEN=*), PARAMETER :: delimiter = ","
      CHARACTER(LEN=100), ALLOCATABLE :: columns(:), headers(:)
      INTEGER :: logg_col, meta_col, teff_col

      OPEN (newunit=unit, FILE=lookup_file, STATUS='old', ACTION='read', IOSTAT=status)
      IF (status /= 0) THEN
         PRINT *, "Error: Could not open file", lookup_file
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      READ (unit, '(A)', IOSTAT=status) line
      IF (status /= 0) THEN
         PRINT *, "Error: Could not read header line"
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      CALL splitline(line, delimiter, headers)
      logg_col = getcolumnindex(headers, "logg")
      teff_col = getcolumnindex(headers, "teff")
      meta_col = getcolumnindex(headers, "meta")
      IF (meta_col < 0) meta_col = getcolumnindex(headers, "feh")

      n_rows = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         n_rows = n_rows + 1
      END DO
      REWIND (unit)
      READ (unit, '(A)', IOSTAT=status) line

      ALLOCATE (out_file_names(n_rows))
      ALLOCATE (out_logg(n_rows), out_meta(n_rows), out_teff(n_rows))

      i = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         i = i + 1
         CALL splitline(line, delimiter, columns)
         out_file_names(i) = columns(1)

         IF (logg_col > 0 .and. columns(logg_col) /= "") THEN
            READ (columns(logg_col), *) out_logg(i)
         ELSE
            out_logg(i) = 0.0
         END IF
         IF (meta_col > 0 .and. columns(meta_col) /= "") THEN
            READ (columns(meta_col), *) out_meta(i)
         ELSE
            out_meta(i) = 0.0
         END IF
         IF (teff_col > 0 .and. columns(teff_col) /= "") THEN
            READ (columns(teff_col), *) out_teff(i)
         ELSE
            out_teff(i) = 0.0
         END IF
      END DO
      CLOSE (unit)

   CONTAINS

      FUNCTION getcolumnindex(headers, target) RESULT(index)
         CHARACTER(LEN=100), INTENT(IN) :: headers(:)
         CHARACTER(LEN=*), INTENT(IN) :: target
         INTEGER :: index, i
         CHARACTER(LEN=100) :: clean_header, clean_target
         index = -1
         clean_target = TRIM(ADJUSTL(target))
         DO i = 1, SIZE(headers)
            clean_header = TRIM(ADJUSTL(headers(i)))
            IF (clean_header == clean_target) THEN
               index = i
               EXIT
            END IF
         END DO
      END FUNCTION getcolumnindex

      SUBROUTINE splitline(line, delimiter, tokens)
         CHARACTER(LEN=*), INTENT(IN) :: line, delimiter
         CHARACTER(LEN=100), ALLOCATABLE, INTENT(OUT) :: tokens(:)
         INTEGER :: num_tokens, pos, start, len_delim
         len_delim = LEN_TRIM(delimiter)
         start = 1; num_tokens = 0
         IF (ALLOCATED(tokens)) DEALLOCATE (tokens)
         DO
            pos = INDEX(line(start:), delimiter)
            IF (pos == 0) EXIT
            num_tokens = num_tokens + 1
            CALL AppendToken(tokens, line(start:start + pos - 2))
            start = start + pos + len_delim - 1
         END DO
         num_tokens = num_tokens + 1
         CALL AppendToken(tokens, line(start:))
      END SUBROUTINE splitline

      SUBROUTINE AppendToken(tokens, token)
         CHARACTER(LEN=*), INTENT(IN) :: token
         CHARACTER(LEN=100), ALLOCATABLE, INTENT(INOUT) :: tokens(:)
         CHARACTER(LEN=100), ALLOCATABLE :: temp(:)
         INTEGER :: n
         IF (.NOT. ALLOCATED(tokens)) THEN
            ALLOCATE (tokens(1)); tokens(1) = token
         ELSE
            n = SIZE(tokens)
            ALLOCATE (temp(n)); temp = tokens
            DEALLOCATE (tokens); ALLOCATE (tokens(n + 1))
            tokens(1:n) = temp; tokens(n + 1) = token
            DEALLOCATE (temp)
         END IF
      END SUBROUTINE AppendToken

   END SUBROUTINE load_lookuptable

   SUBROUTINE loadsed(directory, index, wavelengths, flux)
      CHARACTER(LEN=*), INTENT(IN) :: directory
      INTEGER, INTENT(IN) :: index
      REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, flux
      CHARACTER(LEN=512) :: line
      INTEGER :: unit, n_rows, status, i
      REAL(DP) :: temp_wavelength, temp_flux

      OPEN (newunit=unit, FILE=TRIM(directory), STATUS='OLD', ACTION='READ', IOSTAT=status)
      IF (status /= 0) THEN
         PRINT *, "Error: Could not open file ", TRIM(directory)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) THEN
            PRINT *, "Error: Could not read the file", TRIM(directory)
            CALL mesa_error(__FILE__, __LINE__)
         END IF
         IF (line(1:1) /= "#") EXIT
      END DO

      n_rows = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         n_rows = n_rows + 1
      END DO

      ALLOCATE (wavelengths(n_rows)); ALLOCATE (flux(n_rows))

      REWIND (unit)
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) THEN
            PRINT *, "Error: Could not rewind file", TRIM(directory)
            CALL mesa_error(__FILE__, __LINE__)
         END IF
         IF (line(1:1) /= "#") EXIT
      END DO

      i = 0
      DO
         READ (unit, *, IOSTAT=status) temp_wavelength, temp_flux
         IF (status /= 0) EXIT
         i = i + 1
         wavelengths(i) = temp_wavelength
         flux(i) = temp_flux
      END DO
      CLOSE (unit)
   END SUBROUTINE loadsed

   function remove_dat(path) result(base)
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

END MODULE shared_funcs