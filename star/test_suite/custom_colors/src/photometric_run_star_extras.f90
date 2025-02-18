! ***********************************************************************
!
!   Copyright (C) 2017-2019  Rob Farmer & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************





module run_star_extras

  use star_lib
  use star_def
  use const_def
  use math_lib
  use auto_diff
  use colors_lib

  implicit none


!DEFINE ALL GLOCBAL VARIABLE HERE



  include "test_suite_extras_def.inc"

  ! these routines are called by the standard run_star check_model
  contains

  include "test_suite_extras.inc"


  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
       print *, "Extras startup routine"

    call process_color_files(id, ierr)
    s% extras_startup => extras_startup
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns

    print *, "Sellar atmosphere:", s% x_character_ctrl(1)
    print *, "Instrument:", s% x_character_ctrl(2)

  end subroutine extras_controls





!###########################################################
!## THINGS I HAVE NOT TOUCHED
!###########################################################

  subroutine process_color_files(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type(star_info), pointer :: s
    integer :: i

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

  end subroutine process_color_files


  subroutine extras_startup(id, restart, ierr)
     integer, intent(in) :: id
     logical, intent(in) :: restart
     integer, intent(out) :: ierr
     type (star_info), pointer :: s
     ierr = 0
     call star_ptr(id, s, ierr)
     if (ierr /= 0) return
     call test_suite_startup(s, restart, ierr)
  end subroutine extras_startup


  subroutine extras_after_evolve(id, ierr)
     integer, intent(in) :: id
     integer, intent(out) :: ierr
     type (star_info), pointer :: s
     real(dp) :: dt
     ierr = 0
     call star_ptr(id, s, ierr)
     if (ierr /= 0) return

     write(*,'(a)') 'finished custom colors'

     call test_suite_after_evolve(s, ierr)

  end subroutine extras_after_evolve


  ! returns either keep_going, retry, or terminate.
  integer function extras_check_model(id)
     integer, intent(in) :: id
     integer :: ierr
     type (star_info), pointer :: s
     ierr = 0
     call star_ptr(id, s, ierr)
     if (ierr /= 0) return
     extras_check_model = keep_going
  end function extras_check_model


  INTEGER FUNCTION how_many_extra_profile_columns(id)
     USE star_def, ONLY: star_info
     INTEGER, INTENT(IN) :: id

     INTEGER :: ierr
     TYPE(star_info), POINTER :: s

     ierr = 0
     CALL star_ptr(id, s, ierr)
     IF (ierr /= 0) RETURN

     how_many_extra_profile_columns = 0
  END FUNCTION how_many_extra_profile_columns


  SUBROUTINE data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
     USE star_def, ONLY: star_info, maxlen_profile_column_name
     USE const_def, ONLY: DP
     INTEGER, INTENT(IN) :: id, n, nz
     CHARACTER(LEN=maxlen_profile_column_name) :: names(n)
     REAL(DP) :: vals(nz, n)
     INTEGER, INTENT(OUT) :: ierr

     TYPE(star_info), POINTER :: s

     ierr = 0
     CALL star_ptr(id, s, ierr)
     IF (ierr /= 0) RETURN

  END SUBROUTINE data_for_extra_profile_columns


  ! Returns either keep_going, retry, or terminate
  INTEGER FUNCTION extras_finish_step(id)
     USE chem_def
     INTEGER, INTENT(IN) :: id

     INTEGER :: ierr
     TYPE(star_info), POINTER :: s

     ierr = 0
     CALL star_ptr(id, s, ierr)
     IF (ierr /= 0) RETURN

     extras_finish_step = keep_going
  END FUNCTION extras_finish_step






!###########################################################
!## MESA STUFF
!###########################################################

SUBROUTINE read_strings_from_file(strings, values, n, id)
    INTEGER, INTENT(IN) :: id
    CHARACTER(LEN=512) :: filename
    CHARACTER(LEN=100), ALLOCATABLE :: strings(:)
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: values(:,:)
    INTEGER, INTENT(OUT) :: n
    INTEGER :: unit, status, num_cols, num_rows, i, j, ierr
    CHARACTER(LEN=100), ALLOCATABLE :: headers(:)
    CHARACTER(LEN=1024) :: line
    TYPE(star_info), POINTER :: s

    ! Get file name from star_info structure
    ierr = 0
    CALL star_ptr(id, s, ierr)
    IF (ierr /= 0) RETURN
    filename = TRIM(s%x_character_ctrl(1))

    ! Open the file
    unit = 10
    print *, filename
    OPEN(unit, FILE=filename, STATUS='old', ACTION='read', IOSTAT=status)
    IF (status /= 0) THEN
        PRINT *, "Error: Could not open file", filename
        STOP
    END IF

    ! Skip comment lines until we reach column headers
    DO
        READ(unit, '(A)', IOSTAT=status) line
        IF (status /= 0) EXIT
        IF (line(1:1) /= "#") EXIT  ! Found the header row
    END DO

    ! Parse the header line to get column names
    CALL SplitLine(line, " ", headers, num_cols)

    ! Ensure headers is allocated before using it
    IF (.NOT. ALLOCATED(headers)) THEN
        PRINT *, "Error: Headers array not allocated!"
        STOP
    END IF

    ! Determine how many filter names exist
    n = num_cols - 9
    IF (n < 1) THEN
        PRINT *, "Error: No filter names found in header!"
        STOP
    END IF

    ! Allocate array for filter names
    IF (ALLOCATED(strings)) DEALLOCATE(strings)
    ALLOCATE(strings(n))

    ! Store filter names
    DO i = 1, n
        strings(i) = TRIM(headers(i + 9))
    END DO

    ! Determine number of data rows
    num_rows = 0
    DO
        READ(unit, '(A)', IOSTAT=status) line
        IF (status /= 0) EXIT
        num_rows = num_rows + 1
    END DO

    ! Allocate memory for filter data
    IF (ALLOCATED(values)) DEALLOCATE(values)
    ALLOCATE(values(n, num_rows))

    ! Rewind file and skip headers
    REWIND(unit)
    DO
        READ(unit, '(A)', IOSTAT=status) line
        IF (status /= 0) EXIT
        IF (line(1:1) /= "#") EXIT
    END DO

    ! Read data
    DO j = 1, num_rows
        READ(unit, *, IOSTAT=status) (values(i, j), i = 1, n)
        IF (status /= 0) EXIT
    END DO

    CLOSE(unit)

END SUBROUTINE read_strings_from_file

SUBROUTINE SplitLine(line, delimiter, tokens, num_tokens)
    CHARACTER(LEN=*), INTENT(IN) :: line, delimiter
    CHARACTER(LEN=100), ALLOCATABLE, INTENT(OUT) :: tokens(:)
    INTEGER, INTENT(OUT) :: num_tokens
    INTEGER :: pos, start, end_pos, len_line, count
    CHARACTER(LEN=100), ALLOCATABLE :: temp_tokens(:)

    len_line = LEN_TRIM(line)
    start = 1
    count = 0

    ! Count number of valid tokens while handling consecutive delimiters
    DO WHILE (start <= len_line)
        ! Skip leading delimiters
        DO WHILE (start <= len_line .AND. line(start:start) == delimiter)
            start = start + 1
        END DO
        IF (start > len_line) EXIT

        ! Find end of the token
        pos = INDEX(line(start:), delimiter)
        IF (pos == 0) THEN
            count = count + 1
            EXIT
        ELSE
            count = count + 1
            start = start + pos
        END IF
    END DO

    ! Allocate the tokens array
    IF (ALLOCATED(tokens)) DEALLOCATE(tokens)
    ALLOCATE(tokens(count))
    num_tokens = count

    ! Extract tokens while skipping multiple delimiters
    count = 0
    start = 1

    DO WHILE (start <= len_line)
        ! Skip leading delimiters
        DO WHILE (start <= len_line .AND. line(start:start) == delimiter)
            start = start + 1
        END DO
        IF (start > len_line) EXIT

        ! Find end of the token
        pos = INDEX(line(start:), delimiter)
        IF (pos == 0) THEN
            count = count + 1
            IF (count <= num_tokens) tokens(count) = TRIM(line(start:))
            EXIT
        ELSE
            count = count + 1
            IF (count <= num_tokens) tokens(count) = TRIM(line(start:start + pos - 2))
            start = start + pos
        END IF
    END DO
END SUBROUTINE SplitLine




  !FUNCTIONS FOR OPENING LOOKUP FILE AND FINDING THE NUMBER OF FILES AND THIER FILE PATHS
  integer function how_many_extra_history_columns(id)
      ! Determines how many extra history columns are added based on a file
      integer, intent(in) :: id
      integer :: ierr, n
      character(len=100), allocatable :: strings(:)
      type(star_info), pointer :: s
      REAL(DP), ALLOCATABLE :: filter_values(:,:)

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) then
          how_many_extra_history_columns = 0
          return
      end if

      ! Read strings from the file
      CALL read_strings_from_file(strings, filter_values, n, id)

      ! Number of columns is the size of the strings array
      how_many_extra_history_columns = n

      !print *, "This many columns added to history file:", n

      if (allocated(strings)) deallocate(strings)
  end function how_many_extra_history_columns


SUBROUTINE data_for_extra_history_columns(id, n, names, vals, ierr)
    INTEGER, INTENT(IN) :: id, n
    CHARACTER(LEN=maxlen_history_column_name) :: names(n)
    REAL(dp) :: vals(n)
    INTEGER, INTENT(OUT) :: ierr
    TYPE(star_info), POINTER :: s
    INTEGER :: i
    CHARACTER(LEN=100), ALLOCATABLE :: array_of_strings(:)
    REAL(DP), ALLOCATABLE :: filter_values(:,:)
    INTEGER :: num_strings, num_rows
    REAL(DP) :: teff, log_g, metallicity, filter_data
    CHARACTER(LEN=256) :: filter_name

    ierr = 0
    CALL star_ptr(id, s, ierr)
    IF (ierr /= 0) RETURN

    ! Input parameters
    teff = s%T(1)
    log_g = LOG10(s%grav(1))
    metallicity = s%job%extras_rpar(1)

    print *, teff, log_g, metallicity



    ! Read filter names and data
    CALL read_strings_from_file(array_of_strings, filter_values, num_strings, id)

    ! Assign filter names and data
    IF (ALLOCATED(array_of_strings)) THEN
        DO i = 1, how_many_extra_history_columns(id)
            IF (i <= num_strings) THEN
                filter_name = TRIM(array_of_strings(i))
                filter_data = filter_values(i, 1)  ! Use first row of data
            ELSE
                filter_name = "Unknown"
                filter_data = -1.0_dp
            END IF

            names(i) = filter_name

            ! Use filter data in magnitude calculation
            IF (s%T(1) >= 0 .AND. LOG10(s%grav(1)) >= 0) THEN
                CALL CalculateSyntheticMagnitude(teff, log_g, metallicity, filter_data, ierr)
                IF (ierr /= 0) vals(i) = -1.0_dp
            ELSE
                vals(i) = -1.0_dp
                ierr = 1
            END IF
        END DO
    ELSE
        ierr = 1  ! Indicate an error if array_of_strings is not allocated
    END IF

    ! Cleanup
    IF (ALLOCATED(array_of_strings)) DEALLOCATE(array_of_strings)
    IF (ALLOCATED(filter_values)) DEALLOCATE(filter_values)
END SUBROUTINE data_for_extra_history_columns




  !###########################################################
  !## CUSTOM COLOURS
  !###########################################################

SUBROUTINE CalculateSyntheticMagnitude(teff, log_g, metallicity, filter_data, ierr)
    REAL(8), INTENT(IN) :: teff, log_g, metallicity
    REAL(DP), INTENT(IN) :: filter_data
    INTEGER, INTENT(OUT) :: ierr

    print *, filter_data
    ! Stub implementation - does nothing, just ensures compilation
    ierr = 0  ! No error
    RETURN
END SUBROUTINE CalculateSyntheticMagnitude




  !****************************
  !## LOAD LOOKUP TABLE FOR CLOSEST MATCH FOR MODEL
  !****************************


SUBROUTINE LoadMISTLookupTable(lookup_file, out_teff, out_logg, out_meta, n_rows)
  CHARACTER(LEN=*), INTENT(IN) :: lookup_file
  REAL, ALLOCATABLE, INTENT(OUT) :: out_teff(:), out_logg(:), out_meta(:)
  INTEGER, INTENT(OUT) :: n_rows

  INTEGER :: i, status, unit, num_cols
  CHARACTER(LEN=1024) :: line
  CHARACTER(LEN=100), ALLOCATABLE :: columns(:), headers(:)
  INTEGER :: teff_col, logg_col, meta_col

  ! Open the file
  unit = 10
  OPEN(unit, FILE=lookup_file, STATUS='old', ACTION='read', IOSTAT=status)
  IF (status /= 0) THEN
    PRINT *, "Error: Could not open file", lookup_file
    STOP
  END IF

  ! Skip all lines until we reach the header
  DO
    READ(unit, '(A)', IOSTAT=status) line
    IF (status /= 0) EXIT
    IF (line(1:1) /= "#") EXIT  ! Header line reached
  END DO

  ! Parse the header line
  CALL SplitLine(line, " ", headers, num_cols)

  ! Determine column indices
  teff_col = GetColumnIndex(headers, "log_Teff")
  logg_col = GetColumnIndex(headers, "log_g")
  meta_col = GetColumnIndex(headers, "[Fe/H]")

  IF (teff_col == -1 .OR. logg_col == -1 .OR. meta_col == -1) THEN
    PRINT *, "Error: Required columns not found in the file!"
    STOP
  END IF

  ! Count the number of data rows
  n_rows = 0
  DO
    READ(unit, '(A)', IOSTAT=status) line
    IF (status /= 0) EXIT
    n_rows = n_rows + 1
  END DO
  REWIND(unit)

  ! Skip header again
  DO
    READ(unit, '(A)', IOSTAT=status) line
    IF (line(1:1) /= "#") EXIT
  END DO

  ! Allocate arrays
  ALLOCATE(out_teff(n_rows), out_logg(n_rows), out_meta(n_rows))

  ! Read and parse the file
  i = 0
  DO
    READ(unit, '(A)', IOSTAT=status) line
    IF (status /= 0) EXIT
    i = i + 1

    CALL SplitLine(line, " ", columns, num_cols)

    ! Read values
    READ(columns(teff_col), *) out_teff(i)
    READ(columns(logg_col), *) out_logg(i)
    READ(columns(meta_col), *) out_meta(i)
  END DO

  CLOSE(unit)

CONTAINS

  FUNCTION GetColumnIndex(headers, target) RESULT(index)
    CHARACTER(LEN=100), INTENT(IN) :: headers(:)
    CHARACTER(LEN=*), INTENT(IN) :: target
    INTEGER :: index, j

    index = -1
    DO j = 1, SIZE(headers)
      IF (TRIM(ADJUSTL(headers(j))) == target) THEN
        index = j
        EXIT
      END IF
    END DO
  END FUNCTION GetColumnIndex

  SUBROUTINE SplitLine(line, delimiter, tokens, num_tokens)
    CHARACTER(LEN=*), INTENT(IN) :: line, delimiter
    CHARACTER(LEN=100), ALLOCATABLE, INTENT(OUT) :: tokens(:)
    INTEGER, INTENT(OUT) :: num_tokens
    INTEGER :: pos, start, len_delim

    len_delim = LEN_TRIM(delimiter)
    start = 1
    num_tokens = 0
    IF (ALLOCATED(tokens)) DEALLOCATE(tokens)

    DO
      pos = INDEX(line(start:), delimiter)

      IF (pos == 0) EXIT
      num_tokens = num_tokens + 1
      CALL AppendToken(tokens, line(start:start + pos - 2))
      start = start + pos + len_delim - 1
    END DO

    num_tokens = num_tokens + 1
    CALL AppendToken(tokens, line(start:))
  END SUBROUTINE SplitLine

  SUBROUTINE AppendToken(tokens, token)
    CHARACTER(LEN=*), INTENT(IN) :: token
    CHARACTER(LEN=100), ALLOCATABLE, INTENT(INOUT) :: tokens(:)
    CHARACTER(LEN=100), ALLOCATABLE :: temp(:)
    INTEGER :: n

    IF (.NOT. ALLOCATED(tokens)) THEN
      ALLOCATE(tokens(1))
      tokens(1) = token
    ELSE
      n = SIZE(tokens)
      ALLOCATE(temp(n))
      temp = tokens  ! Backup the current tokens
      DEALLOCATE(tokens)
      ALLOCATE(tokens(n + 1))
      tokens(1:n) = temp
      tokens(n + 1) = token
      DEALLOCATE(temp)
    END IF
  END SUBROUTINE AppendToken

END SUBROUTINE LoadMISTLookupTable





SUBROUTINE GetClosestStellarModels(teff, log_g, metallicity, lu_teff, lu_logg, lu_meta, closest_indices)
  REAL(8), INTENT(IN) :: teff, log_g, metallicity
  REAL, INTENT(IN) :: lu_teff(:), lu_logg(:), lu_meta(:)
  INTEGER, DIMENSION(4), INTENT(OUT) :: closest_indices

  INTEGER :: i, n, j
  REAL :: distance, norm_teff, norm_logg, norm_meta
  REAL, DIMENSION(:), ALLOCATABLE :: scaled_lu_teff, scaled_lu_logg, scaled_lu_meta
  REAL, DIMENSION(4) :: min_distances
  INTEGER, DIMENSION(4) :: indices
  REAL :: teff_min, teff_max, logg_min, logg_max, meta_min, meta_max

  n = SIZE(lu_teff)
  min_distances = HUGE(1.0)
  indices = -1

  ! Find min and max for normalization
  teff_min = MINVAL(lu_teff)
  teff_max = MAXVAL(lu_teff)
  logg_min = MINVAL(lu_logg)
  logg_max = MAXVAL(lu_logg)
  meta_min = MINVAL(lu_meta)
  meta_max = MAXVAL(lu_meta)

  ! Allocate and scale lookup table values
  ALLOCATE(scaled_lu_teff(n), scaled_lu_logg(n), scaled_lu_meta(n))

  scaled_lu_teff = (lu_teff - teff_min) / (teff_max - teff_min)
  scaled_lu_logg = (lu_logg - logg_min) / (logg_max - logg_min)
  scaled_lu_meta = (lu_meta - meta_min) / (meta_max - meta_min)

  ! Normalize input parameters
  norm_teff = (teff - teff_min) / (teff_max - teff_min)
  norm_logg = (log_g - logg_min) / (logg_max - logg_min)
  norm_meta = (metallicity - meta_min) / (meta_max - meta_min)

  ! Debug: !PRINT normalized input parameters
  !PRINT *, "Normalized parameters for target:"
  !PRINT *, "  norm_teff = ", norm_teff, "  norm_logg = ", norm_logg, "  norm_meta = ", norm_meta

  ! Find closest models
  DO i = 1, n
    distance = SQRT((scaled_lu_teff(i) - norm_teff)**2 + &
                    (scaled_lu_logg(i) - norm_logg)**2 + &
                    (scaled_lu_meta(i) - norm_meta)**2)

    ! Check if this distance is smaller than any in the current top 4
    DO j = 1, 4
      IF (distance < min_distances(j)) THEN
        ! Shift larger distances down
        IF (j < 4) THEN
          min_distances(j+1:4) = min_distances(j:3)
          indices(j+1:4) = indices(j:3)
        END IF
        min_distances(j) = distance
        indices(j) = i
        EXIT
      END IF
    END DO
  END DO

  closest_indices = indices

  DEALLOCATE(scaled_lu_teff, scaled_lu_logg, scaled_lu_meta)
END SUBROUTINE GetClosestStellarModels



  SUBROUTINE LinearInterpolate(x, y, x_val, y_val)
    REAL(DP), INTENT(IN) :: x(:), y(:), x_val
    REAL(DP), INTENT(OUT) :: y_val
    INTEGER :: i
    REAL(DP) :: slope

    ! Validate input sizes
    IF (SIZE(x) < 2) THEN
      PRINT *, "Error: x array has fewer than 2 points."
      y_val = 0.0_DP
      RETURN
    END IF

    IF (SIZE(x) /= SIZE(y)) THEN
      PRINT *, "Error: x and y arrays have different sizes."
      y_val = 0.0_DP
      RETURN
    END IF

    ! Handle out-of-bounds cases
    IF (x_val < MINVAL(x)) THEN
      y_val = y(1)
      RETURN
    ELSE IF (x_val > MAXVAL(x)) THEN
      y_val = y(SIZE(y))
      RETURN
    END IF

    ! Perform interpolation
    DO i = 1, SIZE(x) - 1
      IF (x_val >= x(i) .AND. x_val <= x(i + 1)) THEN
        slope = (y(i + 1) - y(i)) / (x(i + 1) - x(i))
        y_val = y(i) + slope * (x_val - x(i))
        RETURN
      END IF
    END DO

    y_val = 0.0_DP
  END SUBROUTINE LinearInterpolate


  SUBROUTINE InterpolateArray(x_in, y_in, x_out, y_out)
    REAL(DP), INTENT(IN) :: x_in(:), y_in(:), x_out(:)
    REAL(DP), INTENT(OUT) :: y_out(:)
    INTEGER :: i

    ! Validate input sizes
    IF (SIZE(x_in) < 2 .OR. SIZE(y_in) < 2) THEN
      PRINT *, "Error: x_in or y_in arrays have fewer than 2 points."
      STOP
    END IF

    IF (SIZE(x_in) /= SIZE(y_in)) THEN
      PRINT *, "Error: x_in and y_in arrays have different sizes."
      STOP
    END IF

    IF (SIZE(x_out) <= 0) THEN
      PRINT *, "Error: x_out array is empty."
      STOP
    END IF

    DO i = 1, SIZE(x_out)
      CALL LinearInterpolate(x_in, y_in, x_out(i), y_out(i))
    END DO
  END SUBROUTINE InterpolateArray
end module run_star_extras
