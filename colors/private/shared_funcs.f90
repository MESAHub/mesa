
MODULE shared_funcs
   USE const_def, ONLY: dp, strlen
   USE utils_lib, ONLY: mesa_error
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: dilute_flux, trapezoidalintegration, rombergintegration, SimpsonIntegration, loadsed, &
             loadfilter, loadvegased, load_lookuptable, remove_dat

CONTAINS

   !---------------------------------------------------------------------------
   ! Apply dilution factor to convert surface flux to observed flux
   !---------------------------------------------------------------------------
   SUBROUTINE dilute_flux(surface_flux, R, d, calibrated_flux)
      REAL(dp), INTENT(IN)  :: surface_flux(:)
      REAL(dp), INTENT(IN)  :: R, d  ! R = stellar radius, d = distance (both in the same units, e.g., cm)
      REAL(dp), INTENT(OUT) :: calibrated_flux(:)

      ! Check that the output array has the same size as the input
      IF (SIZE(calibrated_flux) /= SIZE(surface_flux)) THEN
         PRINT *, "Error in dilute_flux: Output array must have the same size as input array."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ! Apply the dilution factor (R/d)^2 to each element
      calibrated_flux = surface_flux*((R/d)**2)
   END SUBROUTINE dilute_flux

   !###########################################################
   !## MATHS
   !###########################################################

   !****************************
   !Trapezoidal and Simpson Integration For Flux Calculation
   !****************************

   SUBROUTINE trapezoidalintegration(x, y, result)
      REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
      REAL(DP), INTENT(OUT) :: result

      INTEGER :: i, n
      REAL(DP) :: sum

      n = SIZE(x)
      sum = 0.0_dp

      ! Validate input sizes
      IF (SIZE(x) /= SIZE(y)) THEN
         PRINT *, "Error: x and y arrays must have the same size."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      IF (SIZE(x) < 2) THEN
         PRINT *, "Error: x and y arrays must have at least 2 elements."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ! Perform trapezoidal integration
      DO i = 1, n - 1
         sum = sum + 0.5_dp*(x(i + 1) - x(i))*(y(i + 1) + y(i))
      END DO

      result = sum
   END SUBROUTINE trapezoidalintegration

   SUBROUTINE SimpsonIntegration(x, y, result)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
      REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
      REAL(DP), INTENT(OUT) :: result

      INTEGER :: i, n
      REAL(DP) :: sum, h1, h2, f1, f2, f0

      n = SIZE(x)
      sum = 0.0_DP

      ! Validate input sizes
      IF (SIZE(x) /= SIZE(y)) THEN
         PRINT *, "Error: x and y arrays must have the same size."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      IF (SIZE(x) < 2) THEN
         PRINT *, "Error: x and y arrays must have at least 2 elements."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ! Perform adaptive Simpson's rule
      DO i = 1, n - 2, 2
         h1 = x(i + 1) - x(i)       ! Step size for first interval
         h2 = x(i + 2) - x(i + 1)     ! Step size for second interval

         f0 = y(i)
         f1 = y(i + 1)
         f2 = y(i + 2)

         ! Simpson's rule: (h/3) * (f0 + 4f1 + f2)
         sum = sum + (h1 + h2)/6.0_DP*(f0 + 4.0_DP*f1 + f2)
      END DO

      ! Handle the case where n is odd (last interval)
      IF (MOD(n, 2) == 0) THEN
         sum = sum + 0.5_DP*(x(n) - x(n - 1))*(y(n) + y(n - 1))
      END IF

      result = sum
   END SUBROUTINE SimpsonIntegration

   SUBROUTINE rombergintegration(x, y, result)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
      REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
      REAL(DP), INTENT(OUT) :: result

      INTEGER :: i, j, k, n, m
      REAL(DP), DIMENSION(:), ALLOCATABLE :: R
      REAL(DP) :: h, sum, factor

      n = SIZE(x)
      m = INT(LOG(REAL(n, DP))/LOG(2.0_DP)) + 1  ! Number of refinement levels

      ! Validate input sizes
      IF (SIZE(x) /= SIZE(y)) THEN
         PRINT *, "Error: x and y arrays must have the same size."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      IF (n < 2) THEN
         PRINT *, "Error: x and y arrays must have at least 2 elements."
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ALLOCATE (R(m))

      ! Compute initial trapezoidal rule estimate
      h = x(n) - x(1)
      R(1) = 0.5_DP*h*(y(1) + y(n))

      ! Refinement using Romberg's method
      DO j = 2, m
         sum = 0.0_DP
         DO i = 1, 2**(j - 2)
            sum = sum + y(1 + (2*i - 1)*(n - 1)/(2**(j - 1)))
         END DO

         h = h/2.0_DP
         R(j) = 0.5_DP*R(j - 1) + h*sum

         ! Richardson extrapolation
         factor = 4.0_DP
         DO k = j, 2, -1
            R(k - 1) = (factor*R(k) - R(k - 1))/(factor - 1.0_DP)
            factor = factor*4.0_DP
         END DO
      END DO

      result = R(1)
      DEALLOCATE (R)
   END SUBROUTINE rombergintegration

   !-----------------------------------------------------------------------
   ! File I/O functions
   !-----------------------------------------------------------------------

   !****************************
   ! Load Vega SED for Zero Point Calculation
   !****************************
   SUBROUTINE loadvegased(filepath, wavelengths, flux)
      CHARACTER(LEN=*), INTENT(IN) :: filepath
      REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, flux
      CHARACTER(LEN=512) :: line
      INTEGER :: unit, n_rows, status, i
      REAL(dp) :: temp_wave, temp_flux

      unit = 20
      OPEN (unit, FILE=TRIM(filepath), STATUS='OLD', ACTION='READ', IOSTAT=status)
      IF (status /= 0) THEN
         PRINT *, "Error: Could not open Vega SED file ", TRIM(filepath)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ! Skip header line
      READ (unit, '(A)', IOSTAT=status) line
      IF (status /= 0) THEN
         PRINT *, "Error: Could not read header from Vega SED file ", TRIM(filepath)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ! Count the number of data lines
      n_rows = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         n_rows = n_rows + 1
      END DO

      REWIND (unit)
      READ (unit, '(A)', IOSTAT=status) line  ! Skip header again

      ALLOCATE (wavelengths(n_rows))
      ALLOCATE (flux(n_rows))

      i = 0
      DO
         READ (unit, *, IOSTAT=status) temp_wave, temp_flux  ! Ignore any extra columns
         IF (status /= 0) EXIT
         i = i + 1
         wavelengths(i) = temp_wave
         flux(i) = temp_flux
      END DO

      CLOSE (unit)
   END SUBROUTINE loadvegased

   !****************************
   ! Load Filter File
   !****************************
   SUBROUTINE loadfilter(directory, filter_wavelengths, filter_trans)
      CHARACTER(LEN=*), INTENT(IN) :: directory
      REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: filter_wavelengths, filter_trans

      CHARACTER(LEN=512) :: line
      INTEGER :: unit, n_rows, status, i
      REAL(dp) :: temp_wavelength, temp_trans

      ! Open the file
      unit = 20
      OPEN (unit, FILE=TRIM(directory), STATUS='OLD', ACTION='READ', IOSTAT=status)
      IF (status /= 0) THEN
         PRINT *, "Error: Could not open file ", TRIM(directory)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ! Skip header line
      READ (unit, '(A)', IOSTAT=status) line
      IF (status /= 0) THEN
         PRINT *, "Error: Could not read the file", TRIM(directory)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ! Count rows in the file
      n_rows = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         n_rows = n_rows + 1
      END DO

      ! Allocate arrays
      ALLOCATE (filter_wavelengths(n_rows))
      ALLOCATE (filter_trans(n_rows))

      ! Rewind to the first non-comment line
      REWIND (unit)
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) THEN
            PRINT *, "Error: Could not rewind file", TRIM(directory)
            CALL mesa_error(__FILE__, __LINE__)
         END IF
         IF (line(1:1) /= "#") EXIT
      END DO

      ! Read and parse data
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

   !****************************
   ! Load Lookup Table For Identifying Stellar Atmosphere Models
   !****************************
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

      ! Open the file
      unit = 10
      OPEN (unit, FILE=lookup_file, STATUS='old', ACTION='read', IOSTAT=status)
      IF (status /= 0) THEN
         PRINT *, "Error: Could not open file", lookup_file
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ! Read header line
      READ (unit, '(A)', IOSTAT=status) line
      IF (status /= 0) THEN
         PRINT *, "Error: Could not read header line"
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      CALL splitline(line, delimiter, headers)

      ! Determine column indices for logg, meta, and teff
      logg_col = getcolumnindex(headers, "logg")
      teff_col = getcolumnindex(headers, "teff")

      meta_col = getcolumnindex(headers, "meta")
      IF (meta_col < 0) THEN
         meta_col = getcolumnindex(headers, "feh")
      END IF

      n_rows = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         n_rows = n_rows + 1
      END DO
      REWIND (unit)

      ! Skip header
      READ (unit, '(A)', IOSTAT=status) line

      ! Allocate output arrays
      ALLOCATE (out_file_names(n_rows))
      ALLOCATE (out_logg(n_rows), out_meta(n_rows), out_teff(n_rows))

      ! Read and parse the file
      i = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         i = i + 1

         CALL splitline(line, delimiter, columns)

         ! Populate arrays
         out_file_names(i) = columns(1)

         IF (logg_col > 0) THEN
            IF (columns(logg_col) /= "") THEN
               READ (columns(logg_col), *) out_logg(i)
            ELSE
               out_logg(i) = 0.0
            END IF
         ELSE
            out_logg(i) = 0.0
         END IF

         IF (meta_col > 0) THEN
            IF (columns(meta_col) /= "") THEN
               READ (columns(meta_col), *) out_meta(i)
            ELSE
               out_meta(i) = 0.0
            END IF
         ELSE
            out_meta(i) = 0.0
         END IF

         IF (teff_col > 0) THEN
            IF (columns(teff_col) /= "") THEN
               READ (columns(teff_col), *) out_teff(i)
            ELSE
               out_teff(i) = 0.0
            END IF
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
         clean_target = TRIM(ADJUSTL(target))  ! Clean the target string

         DO i = 1, SIZE(headers)
            clean_header = TRIM(ADJUSTL(headers(i)))  ! Clean each header
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
         start = 1
         num_tokens = 0
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
            ALLOCATE (tokens(1))
            tokens(1) = token
         ELSE
            n = SIZE(tokens)
            ALLOCATE (temp(n))
            temp = tokens  ! Backup the current tokens
            DEALLOCATE (tokens)  ! Deallocate the old array
            ALLOCATE (tokens(n + 1))  ! Allocate with one extra space
            tokens(1:n) = temp  ! Restore old tokens
            tokens(n + 1) = token  ! Add the new token
            DEALLOCATE (temp)  ! Clean up temporary array
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

      ! Open the file
      unit = 20
      OPEN (unit, FILE=TRIM(directory), STATUS='OLD', ACTION='READ', IOSTAT=status)
      IF (status /= 0) THEN
         PRINT *, "Error: Could not open file ", TRIM(directory)
         CALL mesa_error(__FILE__, __LINE__)
      END IF

      ! Skip header lines
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) THEN
            PRINT *, "Error: Could not read the file", TRIM(directory)
            CALL mesa_error(__FILE__, __LINE__)
         END IF
         IF (line(1:1) /= "#") EXIT
      END DO

      ! Count rows in the file
      n_rows = 0
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         n_rows = n_rows + 1
      END DO

      ! Allocate arrays
      ALLOCATE (wavelengths(n_rows))
      ALLOCATE (flux(n_rows))

      ! Rewind to the first non-comment line
      REWIND (unit)
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) THEN
            PRINT *, "Error: Could not rewind file", TRIM(directory)
            CALL mesa_error(__FILE__, __LINE__)
         END IF
         IF (line(1:1) /= "#") EXIT
      END DO

      ! Read and parse data
      i = 0
      DO
         READ (unit, *, IOSTAT=status) temp_wavelength, temp_flux
         IF (status /= 0) EXIT
         i = i + 1
         ! Convert f_lambda to f_nu
         wavelengths(i) = temp_wavelength
         flux(i) = temp_flux
      END DO

      CLOSE (unit)

   END SUBROUTINE loadsed

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

END MODULE shared_funcs
