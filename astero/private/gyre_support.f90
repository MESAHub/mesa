! ***********************************************************************
!
!   Copyright (C) 2013-2019  Rich Townsend & The MESA TEAM
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
 
module gyre_support

  ! Uses
      
  use astero_def
  use star_lib
  use star_def
  use const_def
  use utils_lib
  use gyre_lib

  ! No implicit typing
      
  implicit none

  ! Module variables
      
  logical, parameter :: GYRE_IS_ENABLED = .true.

  ! Access specifiers

  private

  public :: GYRE_IS_ENABLED
  public :: init_gyre
  public :: do_gyre_get_modes
  public :: null_gyre_call_back
  public :: store_model_for_gyre
  public :: gyre_call_back
  public :: save_gyre_mode_info
  
  ! Procedures

contains
      
  subroutine init_gyre (gyre_file, ierr)

    use const_def

    character(*), intent(in) :: gyre_file
    integer, intent(out)     :: ierr

    ierr = 0

    ! Initialize GYRE

    call gyre_init(gyre_file)

    ! Set constants

    call gyre_set_constant('G_GRAVITY', standard_cgrav)
    call gyre_set_constant('C_LIGHT', clight)
    call gyre_set_constant('A_RADIATION', crad)

    call gyre_set_constant('M_SUN', Msun)
    call gyre_set_constant('R_SUN', Rsun)
    call gyre_set_constant('L_SUN', Lsun)

    call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')

    write(*,*) 'done init_gyre'

    ! Finish

    return

  end subroutine init_gyre

  !****
      
  subroutine do_gyre_get_modes (s, el, store_model, ierr)

    type (star_info), pointer :: s
    integer, intent(in)       :: el
    logical, intent(in)       :: store_model
    integer, intent(out)      :: ierr

    integer  :: ipar(1)
    real(dp) :: rpar(1)
    integer  :: i
    integer(8) :: time0, time1, clock_rate
    real(dp) :: time
         
    include 'formats'
         
    ierr = 0

    ! If necessary, store the model
         
    if (store_model) then
       call store_model_for_gyre( &
            s, add_center_point, keep_surface_point, add_atmosphere, ierr)
       if (ierr /= 0) then
          write(*,*) 'failed in store_model_for_gyre'
          return
       end if
    end if

    ! Get modes
         
    if (trace_time_in_oscillation_code) then
       call system_clock(time0, clock_rate)
    end if
         
    call gyre_get_modes(el, gyre_call_back, ipar, rpar)  
         
    if (trace_time_in_oscillation_code) then
       call system_clock(time1, clock_rate)
       time = dble(time1-time0)/clock_rate
       total_time_in_oscillation_code = total_time_in_oscillation_code + time
       write(*,1) 'time_in_oscillation_code and total', time, total_time_in_oscillation_code
    end if

    ! Finish

    return

  end subroutine do_gyre_get_modes

  !****
  
  subroutine null_gyre_call_back (md, ipar, rpar, ierr)         
    type(mode_t), intent(in) :: md
    integer, intent(inout)   :: ipar(:)
    real(dp), intent(inout)  :: rpar(:)
    integer, intent(out)     :: ierr
    ierr = 0
  end subroutine null_gyre_call_back

  !****
      
  subroutine store_model_for_gyre (s, add_center_point, keep_surface_point, add_atmosphere, ierr)
            
    type (star_info), intent(in) :: s
    logical, intent(in)          :: add_center_point
    logical, intent(in)          :: keep_surface_point
    logical, intent(in)          :: add_atmosphere
    integer, intent(out)         :: ierr

    real(dp), allocatable     :: global_data(:)
    real(dp), allocatable     :: point_data(:,:)
    !character(:), allocatable :: filename
    character (len=1000)       :: filename  ! temporary until gfortran stops giving bogus warning
         
    logical, parameter :: dbg = .false.
            
    include 'formats'

    ! If necessary, write an FGONG file
         
    if (write_fgong_for_each_model) then

       if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))

       filename = TRIM(astero_results_directory) // '/' // TRIM(fgong_prefix)//TRIM(num_string(s%model_number))//TRIM(fgong_postfix)

       call star_export_pulse_data(s%id, 'FGONG', filename, &
            add_center_point, keep_surface_point, add_atmosphere, ierr)

       if (ierr /= 0) then
          write(*,*) 'failed in export_pulse_data'
          call mesa_error(__FILE__,__LINE__)
       end if

    end if

    ! Get the model data

    if (dbg) write(*,2) 'call star_get_pulse_data (GYRE)', s%model_number

    call star_get_pulse_data(s%id, 'GYRE', &
         add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)

    if (ierr /= 0) then
       write(*,*) 'failed in star_get_pulse_data'
       call mesa_error(__FILE__,__LINE__)
    end if

    if (dbg) write(*,2) 'done star_pulse_data (GYRE)', s%model_number

    ! Pass the data to GYRE
         
    if (dbg) write(*,2) 'call gyre_set_model', s%model_number

    call gyre_set_model(global_data, point_data, GYRE_MODEL_VERSION)

    if (dbg) write(*,2) 'done gyre_set_model', s%model_number

    ! If necessary, write a GYRE file
         
    if (write_gyre_for_each_model) then

       if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))

       filename = TRIM(astero_results_directory) // '/' // TRIM(gyre_prefix)//TRIM(num_string(s%model_number))//TRIM(gyre_postfix)

       call star_write_pulse_data(s%id, 'GYRE', filename, &
            global_data, point_data, ierr)

       if (ierr /= 0) then
          write(*,*) 'failed in star_write_pulse_data (GYRE)'
          stop
       end if

    end if

    ! Finish

    return

  contains

    function num_string (n)

      integer, intent(in) :: n
      character(12)       :: num_string

      character(6) :: format_string

      ! Convert the number to a string

       write(format_string, 100) model_num_digits, model_num_digits
100    format('I',I2.2,'.',I2.2)
               
       write(num_string, format_string) s%model_number

       ! Finish

       return

     end function num_string

  end subroutine store_model_for_gyre

  !****

  subroutine gyre_call_back(md, ipar, rpar, ierr)

    use astero_def, only: store_new_oscillation_results
         
    type(mode_t), intent(in) :: md
    integer, intent(inout)   :: ipar(:)
    real(dp), intent(inout)  :: rpar(:)
    integer, intent(out)     :: ierr
    
    integer  :: new_el, new_order, new_em
    real(dp) :: new_inertia, new_cyclic_freq, new_growth_rate
    
    include 'formats'

    ierr = 0

    new_el = md% md_p% l
    new_order = md% n_pg
    new_inertia = md% E_norm()
    new_cyclic_freq = REAL(md% freq('UHZ'))
    new_growth_rate = AIMAG(md% freq('RAD_PER_SEC'))
    new_em = 0
         
    call store_new_oscillation_results( &
         new_el, new_order, new_em, new_inertia, new_cyclic_freq, new_growth_rate, ierr)
         
    call save_gyre_mode_info( &
         new_el, new_order, new_em, new_inertia, new_cyclic_freq, new_growth_rate, &
         md, ipar, rpar, ierr)

    ! Finish

    return

  end subroutine gyre_call_back

  !****
      
  subroutine save_gyre_mode_info( &
       new_el, new_order, new_em, new_inertia, new_cyclic_freq, new_growth_rate, &
       md, ipar, rpar, ierr)

    integer, intent(in)      :: new_el, new_order, new_em
    real(dp), intent(in)     :: new_inertia, new_cyclic_freq, new_growth_rate
    type(mode_t), intent(in) :: md
    integer, intent(inout)   :: ipar(:)
    real(dp), intent(inout)  :: rpar(:)
    integer, intent(out)     :: ierr

    integer :: iounit, i, j, skip, nn
    real(dp) :: y_r, y_h

    include 'formats'
         
    !if (use_other_gyre_mode_info) then
    !   call astero_other_procs% other_gyre_mode_info(md, ipar, rpar, ierr)
    !end if
         
    if (star_model_number /= save_mode_model_number) return
    if (new_el /= el_to_save .or. new_order /= order_to_save) return
         
    if (len_trim(save_mode_filename) <= 0) save_mode_filename = 'save_mode.data'
    write(*,*) 'save eigenfunction info to file ' // trim(save_mode_filename)
    write(*,'(3a8,99a20)') 'el', 'order', 'em', 'freq (microHz)', 'inertia', 'growth rate (s)'
    write(*,'(3i8,f20.10,2e20.10,i20)') &
         new_el, new_order, new_em, new_cyclic_freq, new_inertia, new_growth_rate
    ierr = 0
    iounit = alloc_iounit(ierr)
    if (ierr /= 0) return
    open(unit=iounit, file=trim(save_mode_filename), action='write', iostat=ierr)
    if (ierr /= 0) return
    write(iounit,'(3a8,99a20)') 'el', 'order', 'em', 'freq (microHz)', 'inertia', 'growth rate (s)'!, 'nn'
    write(iounit,'(3i8,f20.10,2e20.10,i20)') &
         new_el, new_order, new_em, new_cyclic_freq, new_inertia, new_growth_rate

    close(iounit)

    call free_iounit(iounit)         

  end subroutine save_gyre_mode_info

end module gyre_support
