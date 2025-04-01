
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

  !FUNCTIONS FOR OPENING LOOKUP FILE AND FINDING THE NUMBER OF FILES AND THIER FILE PATHS
  integer function how_many_extra_history_columns(id)
      ! Determines how many extra history columns are added based on a file
      integer, intent(in) :: id
      integer :: ierr, n
      character(len=100), allocatable :: strings(:)
      type(star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) then
          how_many_extra_history_columns = 0
          return
      end if

      ! Read strings from the file
      call read_strings_from_file(strings, n, id)

      ! Number of columns is the size of the strings array
      how_many_extra_history_columns = n + 2

      !print *, "This many columns added to history file:", n

      if (allocated(strings)) deallocate(strings)
  end function how_many_extra_history_columns


  function basename(path) result(base)
      ! Extracts the base name from a given file path
      character(len=*), intent(in) :: path
      character(len=512) :: base
      integer :: last_slash

      ! Find the position of the last slash
      last_slash = len_trim(path)
      do while (last_slash > 0 .and. path(last_slash:last_slash) /= '/')
          last_slash = last_slash - 1
      end do

      ! Extract the base name
      base = path(last_slash+1:)
  end function basename


  subroutine read_strings_from_file(strings, n, id)
      ! Reads strings from a file into an allocatable array
      integer, intent(in) :: id
      character(len=512) :: filename
      character(len=100), allocatable :: strings(:)
      integer, intent(out) :: n
      integer :: unit, i, status
      character(len=100) :: line
      integer :: ierr
      type(star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! Construct the filename
      filename = trim(s%x_character_ctrl(2)) // "/" // trim(basename(s%x_character_ctrl(2)))

      ! Initialize
      n = 0

      ! Open the file
      unit = 10
      open(unit, file=filename, status='old', action='read', iostat=status)
      if (status /= 0) then
          print *, "Error: Could not open file", filename
          stop
      end if

      ! Count lines in the file to determine the size of the array
      do
          read(unit, '(A)', iostat=status) line
          if (status /= 0) exit
          n = n + 1  ! for bolometric correctionms
      end do
      rewind(unit)

      ! Allocate the array and read the strings
      if (allocated(strings)) deallocate(strings)
      allocate(strings(n))
      do i = 1, n
          read(unit, '(A)') strings(i)
      end do

      close(unit)
  end subroutine read_strings_from_file



  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      ! Populates data for the extra history columns
      integer, intent(in) :: id, n
      integer, intent(out) :: ierr
      character(len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer :: i, num_strings
      character(len=100), allocatable :: array_of_strings(:)
      real(dp) :: teff, log_g, metallicity, R, d,  bolometric_magnitude, bolometric_flux
      character(len=256) :: sed_filepath, filter_filepath, filter_name, filter_dir, vega_filepath
      real(dp), dimension(:), allocatable :: wavelengths, fluxes, filter_wavelengths, filter_trans
      logical :: make_sed

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! Extract input parameters
      teff = s%T(1)
      log_g = LOG10(s%grav(1))
      R = s%R(1)  ! * 1d3
      metallicity = s%job%extras_rpar(1)
      d = s%job%extras_rpar(2)

      sed_filepath = s%x_character_ctrl(1)
      filter_dir = s%x_character_ctrl(2)
      vega_filepath = s%x_character_ctrl(3)
      make_sed = trim(adjustl(s%x_character_ctrl(4))) == 'true'

      ! Read filters from file
      if (allocated(array_of_strings)) deallocate(array_of_strings)
      allocate(array_of_strings(n))
      call read_strings_from_file(array_of_strings, num_strings, id)

      !PRINT *, "################################################"

      ! Compute bolometric values
      CALL calculatebolometric(teff, log_g, metallicity, R, d,  bolometric_magnitude, bolometric_flux, wavelengths, fluxes, sed_filepath)
      names(1) = "Mag_bol"
      vals(1) = bolometric_magnitude
      names(2) = "Flux_bol"
      vals(2) = bolometric_flux

      ! Populate history columns
      if (allocated(array_of_strings)) then
          do i = 3, how_many_extra_history_columns(id)
              filter_name = "Unknown"
              if (i <= num_strings + 2) filter_name = trim(remove_dat(array_of_strings(i - 2)))
              names(i) = filter_name
              filter_filepath = trim(filter_dir) // "/" // array_of_strings(i - 2)

              if (teff >= 0 .and. log_g >= 0 .and. metallicity >= 0) then
                  vals(i) = calculatesynthetic(teff, log_g, metallicity, ierr, wavelengths, fluxes, filter_wavelengths, filter_trans, filter_filepath, vega_filepath, array_of_strings(i - 2), make_sed)
                  if (ierr /= 0) vals(i) = -1.0_dp
              else
                  vals(i) = -1.0_dp
                  ierr = 1
              end if
              !PRINT *, names(i), vals(i)
          end do
      else
          ierr = 1  ! Indicate an error if array_of_strings is not allocated
      end if

      if (allocated(array_of_strings)) deallocate(array_of_strings)
  end subroutine data_for_extra_history_columns



end module run_star_extras