module run_star_extras
  use star_lib, only: star_ptr
  use star_def, only: star_info, maxlen_history_column_name, maxlen_profile_column_name, keep_going
  use const_def, only: dp, i8, mesa_dir
  ! TODO: below are things we will need to incorporate into the main code
  ! So that we do not need to create a custom run_stars_extras file
  use colors_def, only: Colors_General_Info
  use colors_lib, only: colors_ptr, calculate_bolometric, calculate_synthetic, remove_dat, read_strings_from_file
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

    call colors_ptr(s% colors_handle, colors_settings, ierr)
    if (ierr /= 0) then
        print *, "Error getting colors_settings"
        return
    end if
    ! Debug print statements to check that values are properly loaded from inlist
    write(*,*) 'DEBUG: colors_handle = ', s% colors_handle
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
  end subroutine extras_startup

   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
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
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
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

      call read_strings_from_file(colors_settings, strings, n)
      ! TODO: move this into mesa
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
      real(dp) :: teff, log_g, metallicity, R, d, bolometric_magnitude,&
       bolometric_flux
      character(len=256) :: sed_filepath, filter_filepath, filter_name,&
       filter_dir, vega_filepath
      real(dp), dimension(:), allocatable :: wavelengths, fluxes, &
      filter_wavelengths, filter_trans
      logical :: make_sed
      integer :: n_tmp

      ! TODO: move this into mesa

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
      call read_strings_from_file(colors_settings, array_of_strings, n_tmp)
      num_strings = n

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
                    make_sed, colors_settings% colors_results_directory)
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
