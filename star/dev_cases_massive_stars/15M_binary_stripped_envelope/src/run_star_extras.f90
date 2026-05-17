! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
  use chem_def
  use num_lib
  use binary_def
  use ionization_def

  implicit none

  ! these routines are called by the standard run_star check_model

contains


  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! this is the place to set any procedure pointers you want to change
    ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


    ! the extras functions in this file will not be called
    ! unless you set their function pointers as done below.
    ! otherwise we use a null_ version which does nothing (except warn).

    s% extras_startup => extras_startup
    s% extras_start_step => extras_start_step
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns

    s% how_many_extra_history_header_items => how_many_extra_history_header_items
    s% data_for_extra_history_header_items => data_for_extra_history_header_items
    s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
    s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

  end subroutine extras_controls

  subroutine extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
 
    ! Initialize variables on startup.
    s% lxtra(11) = .true. ! do we still need to read inlist_to_CC? -> yes.

  end subroutine extras_startup


  integer function extras_start_step(id)
    integer, intent(in) :: id
    integer :: ierr
    integer :: binary_id
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_start_step = 0

   ! these two if statements below could be combined. -EbF
    if ((s% center_h1 < 1.0d-2) .and. &
        (s% center_he4 < 1.0d-4 .or. s% log_center_temperature >= 9.3d0) &
        .and. (s% center_c12 < 2.0d-2)) then

       if(s% x_logical_ctrl(1) .and. s% lxtra(11)) then !check for central carbon depletion, only in case we run single stars.
          print *,  "*** Single star depleted carbon ***"
          print *, "read inlist_to_CC"
          call read_star_job(s, "inlist_to_cc", ierr)
          if (ierr /= 0) then
             print *, "Failed reading star_job in inlist_to_CC"
             return
          end if
          print *, "read star_job from inlist_to_CC"

          call star_read_controls(id, "inlist_to_cc", ierr)
          if (ierr /= 0) then
             print *, "Failed reading controls in inlist_to_CC"
             return
          end if
          print *, "read controls from inlist_to_CC"
          s% lxtra(11) = .false. ! avoid re-entering here
       end if

    end if

    ! We can also soften the surface BC by switching to a hydrodynamic BC.
    if (s% log_center_temperature >= 9.5d0 .and. s% center_si28 <1d-3) then
      s% use_compression_outer_BC = .true.
    end if

  end function extras_start_step


  ! returns either keep_going, retry, backup, or terminate.
  integer function extras_check_model(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    real(dp) :: error, atol, rtol
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going


    ! by default, indicate where (in the code) MESA terminated
    if (extras_check_model == terminate) s% termination_code = t_extras_check_model
  end function extras_check_model


  integer function how_many_extra_history_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_columns = 1
  end function how_many_extra_history_columns


  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
    use chem_def, only: chem_isos
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer :: k
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    names(1) = 'fe_core_infall_speed'

    k = 0
    if (s% fe_core_mass > 0.0d0) then
       k = s% nz
       do while (s% m(k) <= s% fe_core_mass * Msun)
          k = k-1 ! loop outwards
       end do
        ! We multiply by -1, since v is negative. This yields a positive infall speed, e.g. 50 km/s.
        vals(1) = - min(0d0, minval(s%v(k:s%nz))/1d5)
     else
        vals(1) = 0d0
     end if

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
    integer :: k
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine data_for_extra_profile_columns


  integer function how_many_extra_history_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_header_items = 0
  end function how_many_extra_history_header_items


  subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    integer :: i
  end subroutine data_for_extra_history_header_items


  integer function how_many_extra_profile_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_header_items = 0
  end function how_many_extra_profile_header_items


  subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return
  end subroutine data_for_extra_profile_header_items


  ! returns either keep_going or terminate.
  ! note: cannot request retry or backup; extras_check_model can do that.
  integer function extras_finish_step(id)
    integer, intent(in) :: id
    integer :: ierr, k
    real(dp) :: fe_core_infall_limit, current_fe_core_infall
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going

    ! Initialize variables
    k = 0
    current_fe_core_infall = 0
    fe_core_infall_limit = max(500d0,s% x_ctrl(1)) ! in case s% x_ctrl(1) is not read on restart.
    ! exit condition for star #2
    if (.not. s% x_logical_ctrl(1)) then ! stop
       print *,  "*** Single star depleted carbon ***"
       extras_finish_step = terminate
    end if

!  Custom Fe core collapse condition from Mathieu Renzo.
    if (s% fe_core_mass > 0.0d0) then
       k = s% nz
       do while (s% m(k) <= s% fe_core_mass * Msun)
          k = k-1 ! loop outwards
       end do
       ! k is now the outer index of the fe core
       current_fe_core_infall = - min(0d0, minval(s%v(k:s%nz))/1d5)
       if (current_fe_core_infall > 2d0) then ! avoid print out until the infall speed exceeds 2 km/s
          write(*,*) 'fe_core_infall = ', current_fe_core_infall, 'km/s'! (-)
          write (*,*) 'fe_core_infall limit', fe_core_infall_limit, 'km/s'
          if (current_fe_core_infall >= fe_core_infall_limit) then
             s% termination_code = t_fe_core_infall_limit
             write(*, '(/,a,/, 99e20.10)') &
                  'stop because fe_core_infall > fe_core_infall_limit', &
                  current_fe_core_infall,  fe_core_infall_limit
             print *, "treshold v used", maxval(abs(s%v(k:s%nz)))
             extras_finish_step = terminate
          end if
       end if
    end if

    if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
  end function extras_finish_step


  subroutine extras_after_evolve(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine extras_after_evolve
end module run_star_extras
