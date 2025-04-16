! Modified run_star_extras.f90
module run_star_extras
  use star_lib
  use star_def
  use const_def
  use math_lib
  use auto_diff
  use colors_def
  use colors_lib

  implicit none

  ! Remove the col pointer declaration entirely

  include "test_suite_extras_def.inc"
  ! these routines are called by the standard run_star check_model
  contains

subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    print *, "Extras startup routine"

    ! Initialize the colors module
    call colors_init(ierr)
    if (ierr /= 0) then
        print *, "Error initializing colors module"
        return
    end if

    ! Set the colors parameters directly in star_info
    s%color_instrument = 'data/filters/GAIA/GAIA'
    s%color_vega_sed = 'data/stellar_models/vega_flam.csv'
    s%color_atm = 'data/stellar_models/Kurucz2003all/'
    s%color_z = 0.0d0
    s%color_d = 3.0857d17  ! 10pc for abs mag
    s%color_make_csv = .false.

    ! Register callbacks
    s%extras_startup         => extras_startup
    s%extras_check_model     => extras_check_model
    s%extras_finish_step     => extras_finish_step
    s%extras_after_evolve    => extras_after_evolve
    s%how_many_extra_history_columns => how_many_extra_history_columns
    s%data_for_extra_history_columns => data_for_extra_history_columns
    s%how_many_extra_profile_columns => how_many_extra_profile_columns
    s%data_for_extra_profile_columns => data_for_extra_profile_columns
end subroutine extras_controls

  subroutine process_color_files(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type(star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! Debug print statements to check that values are properly loaded from inlist
    write(*,*) 'DEBUG: color_instrument = ', trim(s%color_instrument)
    write(*,*) 'DEBUG: color_vega_sed = ', trim(s%color_vega_sed)
    write(*,*) 'DEBUG: color_atm = ', trim(s%color_atm)
    write(*,*) 'DEBUG: color_z = ', s%color_z
    write(*,*) 'DEBUG: color_d = ', s%color_d
    write(*,*) 'DEBUG: color_make_csv = ', s%color_make_csv
  end subroutine process_color_files

  subroutine extras_startup(id, restart, ierr)
     integer, intent(in) :: id
     logical, intent(in) :: restart
     integer, intent(out) :: ierr
     type (star_info), pointer :: s
     ierr = 0
     call star_ptr(id, s, ierr)
     if (ierr /= 0) return
     call test_suite_startup(restart, ierr)
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

    ! Remove cleanup for col pointer

    call test_suite_after_evolve(ierr)
  end subroutine extras_after_evolve

  
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

  integer function how_many_extra_history_columns(id)
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
   write(*,*) 'how_many: color_vega_sed = ', trim(s%color_vega_sed)
      call read_strings_from_file(strings, n, id)
      how_many_extra_history_columns = n + 2
      if (allocated(strings)) deallocate(strings)
  end function how_many_extra_history_columns


function basename(path) result(name)
   character(len=*), intent(in) :: path
   character(len=strlen) :: name
   integer :: i
   if (len_trim(path) == 0) then
      name = ''
      return
   end if
   i = index(path, '/', back=.true.)
   name = path(i+1:)
end function

  subroutine read_strings_from_file(strings, n, id)
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


      !s%color_instrument = s%x_character_ctrl(2)
      write(*,*) 's%color_vega_sed = ', trim(s%color_vega_sed)
      write(*,*) 's%color_instrument = ', trim(s%color_instrument)

      filename = trim(s%color_instrument) // "/" // &
      trim(basename(s%color_instrument))
      n = 0
      unit = 10
      open(unit, file=filename, status='old', action='read', iostat=status)
      if (status /= 0) then
          print *, "Error: Could not open file", filename
          stop
      end if

      do
          read(unit, '(A)', iostat=status) line
          if (status /= 0) exit
          n = n + 1
      end do
      rewind(unit)
      if (allocated(strings)) deallocate(strings)
      allocate(strings(n))
      do i = 1, n
          read(unit, '(A)') strings(i)
      end do
      close(unit)
  end subroutine read_strings_from_file

  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      integer, intent(out) :: ierr
      character(len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer :: i, num_strings
      character(len=100), allocatable :: array_of_strings(:)
      real(dp) :: teff, log_g, metallicity, R, d, bolometric_magnitude,&
       bolometric_flux
      character(len=256) :: sed_filepath, filter_filepath, filter_name,&
       filter_dir, vega_filepath
      real(dp), dimension(:), allocatable :: wavelengths, fluxes, &
      filter_wavelengths, filter_trans
      logical :: make_sed

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      teff   = s%T(1)
      log_g  = LOG10(s%grav(1))
      R      = s%R(1)


      metallicity = s%color_z
      d = s%color_d
      sed_filepath = s%color_atm
      filter_dir = s%color_instrument
      vega_filepath = s%color_vega_sed
      make_sed = s%color_make_csv

      ! Use s% just like s%
      print *, "Stellar atmosphere:", trim(s%color_atm)
      print *, "Instrument:", trim(s%color_instrument)

      if (allocated(array_of_strings)) deallocate(array_of_strings)
      allocate(array_of_strings(n))
      call read_strings_from_file(array_of_strings, num_strings, id)

      CALL calculatebolometric(teff, log_g, metallicity, R, d,&
       bolometric_magnitude, bolometric_flux, wavelengths, fluxes, sed_filepath)
      names(1) = "Mag_bol"
      vals(1) = bolometric_magnitude
      names(2) = "Flux_bol"
      vals(2) = bolometric_flux

      if (allocated(array_of_strings)) then
          do i = 3, how_many_extra_history_columns(id)
              filter_name = "Unknown"
              if (i <= num_strings + 2) filter_name =&
               trim(remove_dat(array_of_strings(i - 2)))
              names(i) = filter_name
              filter_filepath = trim(filter_dir) // "/" // array_of_strings(i - 2)

              if (teff >= 0 .and. log_g >= 0 .and. metallicity >= 0) then
                  vals(i) = calculatesynthetic(teff, log_g, metallicity, ierr,&
                   wavelengths, fluxes, filter_wavelengths, filter_trans, &
                   filter_filepath, vega_filepath, array_of_strings(i - 2),&
                    make_sed)
                  if (ierr /= 0) vals(i) = -1.0_dp
              else
                  vals(i) = -1.0_dp
                  ierr = 1
              end if
          end do
      else
          ierr = 1
      end if

      if (allocated(array_of_strings)) deallocate(array_of_strings)
  end subroutine data_for_extra_history_columns

end module run_star_extras
