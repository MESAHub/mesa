! ***********************************************************************
!
!   run_star_extras file for custom colors
!
! ***********************************************************************

module run_star_extras
  use star_lib
  use star_def
  use const_def
  use math_lib
  ! Custom colors imports
  use colors_def, only: Colors_General_Info
  use colors_lib, only: colors_init, colors_ptr, alloc_colors_handle_using_inlist, calculate_bolometric, calculate_synthetic, remove_dat

  implicit none

  ! Colors module settings
  type (colors_General_Info), pointer :: colors_settings

  ! These routines are called by the standard run_star check_model
  contains

  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s

    integer :: colors_handle

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    
    ! Initialize the colors module
    call colors_init(.false., '', ierr)
    colors_handle = alloc_colors_handle_using_inlist(s% inlist_fname, ierr)
    call colors_ptr(colors_handle, colors_settings, ierr)
    if (ierr /= 0) then
        print *, "Error initializing colors module"
        return
    end if

    ! Debug info about colors configuration
    write(*,*) 'Custom colors initialized with:'
    write(*,*) '  instrument = ', trim(colors_settings% instrument)
    write(*,*) '  vega_sed = ', trim(colors_settings% vega_sed)
    write(*,*) '  stellar_atm = ', trim(colors_settings% stellar_atm)
    write(*,*) '  metallicity = ', colors_settings% metallicity
    write(*,*) '  distance = ', colors_settings% distance
    write(*,*) '  make_csv = ', colors_settings% make_csv
    write(*,*) '  use_colors = ', colors_settings% use_colors

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

    ! Any additional startup operations can go here
  end subroutine extras_startup

  integer function extras_check_model(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going
  end function extras_check_model

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

  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    integer, intent(out) :: ierr
    character(len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer :: i, num_strings
    character(len=100), allocatable :: array_of_strings(:)
    real(dp) :: teff, log_g, metallicity, R, d, bolometric_magnitude, bolometric_flux
    character(len=256) :: sed_filepath, filter_filepath, filter_name, filter_dir, vega_filepath
    real(dp), dimension(:), allocatable :: wavelengths, fluxes, filter_wavelengths, filter_trans
    logical :: make_sed

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! Get stellar parameters for colors calculation
    teff   = s%T(1)
    log_g  = LOG10(s%grav(1))
    R      = s%R(1)
    metallicity = s%kap_rq%Zbase

    ! Get colors settings
    d = colors_settings% distance
    sed_filepath = colors_settings% stellar_atm
    filter_dir = colors_settings% instrument
    vega_filepath = colors_settings% vega_sed
    make_sed = colors_settings% make_csv

    if (allocated(array_of_strings)) deallocate(array_of_strings)
    allocate(array_of_strings(n))
    call read_strings_from_file(array_of_strings, num_strings, id)

    ! Calculate bolometric values
    CALL calculate_bolometric(teff, log_g, metallicity, R, d, &
     bolometric_magnitude, bolometric_flux, wavelengths, fluxes, sed_filepath)
    names(1) = "Mag_bol"
    vals(1) = bolometric_magnitude
    names(2) = "Flux_bol"
    vals(2) = bolometric_flux

    ! Calculate filter-specific magnitudes
    if (allocated(array_of_strings)) then
      do i = 3, how_many_extra_history_columns(id)
        filter_name = "Unknown"
        if (i <= num_strings + 2) filter_name = trim(remove_dat(array_of_strings(i - 2)))
        names(i) = filter_name
        filter_filepath = trim(filter_dir) // "/" // array_of_strings(i - 2)

        if (teff >= 0 .and. log_g >= 0 .and. metallicity >= 0) then
          vals(i) = calculate_synthetic(teff, log_g, metallicity, ierr, &
           wavelengths, fluxes, filter_wavelengths, filter_trans, &
           filter_filepath, vega_filepath, array_of_strings(i - 2), make_sed)
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

  integer function how_many_extra_profile_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 0
  end function how_many_extra_profile_columns

  subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
    integer, intent(in) :: id, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    ! No profile columns needed for custom colors
  end subroutine data_for_extra_profile_columns

  integer function extras_finish_step(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going
  end function extras_finish_step

  subroutine extras_after_evolve(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    write(*,'(a)') 'Finished custom colors evolution'
  end subroutine extras_after_evolve

  ! Helper functions
  
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
    character(len=100), allocatable :: strings(:)
    integer, intent(out) :: n
    integer :: unit, i, status
    character(len=100) :: line
    integer :: ierr
    type(star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    filename = trim(colors_settings% instrument) // "/" // &
    trim(basename(colors_settings% instrument))
    n = 0
    unit = 10
    open(unit, file=filename, status='old', action='read', iostat=status)
    if (status /= 0) then
      print *, "Error: Could not open file", filename
      return
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

end module run_star_extras
