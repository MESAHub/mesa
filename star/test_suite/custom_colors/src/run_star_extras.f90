module run_star_extras
  use star_lib, only: star_ptr
  use star_def, only: star_info, maxlen_history_column_name, maxlen_profile_column_name, keep_going
  use const_def, only: dp, strlen, mesa_dir
  ! TODO: below are things we will need to incorporate into the main code
  ! So that we do not need to create a custom run_stars_extras file
  use utils_lib, only: mkdir  ! Add this line for directory creation
  use colors_def, only: Colors_General_Info
  use colors_lib, only: colors_init, colors_ptr, alloc_colors_handle_using_inlist, &
                        calculate_bolometric, calculate_synthetic, remove_dat
  implicit none
  private
  public :: extras_controls  ! Add this line
  type (colors_General_Info), pointer :: colors_settings
  include "test_suite_extras_def.inc"

  contains

subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s

    integer :: colors_handle

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    print *, "Extras startup routine"

    ! Initialize the colors module
    ! TODO: we will have to do this somewhere in MESA, dig around how/where kap does this.
    ! this is a little hacky
    call colors_init(.false., '', ierr)
    colors_handle = alloc_colors_handle_using_inlist(s% inlist_fname, ierr)
    call colors_ptr(colors_handle,colors_settings,ierr)
    if (ierr /= 0) then
        print *, "Error initializing colors module"
        return
    end if


    call mkdir('LOGS')
    call mkdir('LOGS/SED')

    ! Debug print statements to check that values are properly loaded from inlist
    write(*,*) 'DEBUG: colors_handle = ', colors_handle
    write(*,*) 'DEBUG: colors instrument = ', trim(colors_settings% instrument)
    write(*,*) 'DEBUG: colors vega_sed = ', trim(colors_settings% vega_sed)
    write(*,*) 'DEBUG: colors stellar_atm = ', trim(colors_settings% stellar_atm)
    write(*,*) 'DEBUG: colors metallicity = ', colors_settings% metallicity
    write(*,*) 'DEBUG: colors distance = ', colors_settings% distance
    write(*,*) 'DEBUG: colors make_csv = ', colors_settings% make_csv
    write(*,*) 'DEBUG: colors use_colors = ', colors_settings% use_colors

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

  subroutine extras_startup(id, restart, ierr)
     integer, intent(in) :: id
     logical, intent(in) :: restart
     integer, intent(out) :: ierr
     type (star_info), pointer :: s
     ierr = 0
     call star_ptr(id, s, ierr)
     if (ierr /= 0) return
     !XXXcall test_suite_startup(restart, ierr)
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

    !XXXcall test_suite_after_evolve(ierr)
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
     INTEGER, INTENT(IN) :: id

     INTEGER :: ierr
     TYPE(star_info), POINTER :: s

     ierr = 0
     CALL star_ptr(id, s, ierr)
     IF (ierr /= 0) RETURN

     how_many_extra_profile_columns = 0
  END FUNCTION how_many_extra_profile_columns

  SUBROUTINE data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
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
end function basename

  subroutine read_strings_from_file(strings, n, id)
      integer, intent(in) :: id
      character(len=512) :: filename
      character(len=100), allocatable, intent(out) :: strings(:)
      integer, intent(out) :: n
      integer :: unit, i, status
      character(len=100) :: line
      integer :: ierr
      type(star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return


      !write(*,*) 'vega_sed = ', trim(colors_settings% vega_sed)
      !write(*,*) 'instrument = ', trim(colors_settings% instrument)

      filename = trim(mesa_dir) // trim(colors_settings% instrument) // "/" // &
      trim(basename(colors_settings% instrument))
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
      metallicity = s%kap_rq%Zbase

      !metallicity = colors_settings% metallicity
      d = colors_settings% distance
      sed_filepath = trim(mesa_dir) // colors_settings% stellar_atm
      filter_dir = trim(mesa_dir) // colors_settings% instrument
      vega_filepath = trim(mesa_dir) // colors_settings% vega_sed
      make_sed = colors_settings% make_csv

      if (allocated(array_of_strings)) deallocate(array_of_strings)
      allocate(array_of_strings(n))
      call read_strings_from_file(array_of_strings, num_strings, id)

      CALL calculate_bolometric(teff, log_g, metallicity, R, d,&
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
                  vals(i) = calculate_synthetic(teff, log_g, metallicity, ierr,&
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
