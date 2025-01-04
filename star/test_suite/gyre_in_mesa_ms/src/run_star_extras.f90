! ***********************************************************************
!
!   Copyright (C) 2018-2019  Rich Townsend, The MESA Team
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

  ! Uses

  use star_lib
  use star_def
  use const_def
  use math_lib
  use auto_diff

  use gyre_mesa_m

  ! No implicit typing

  implicit none
  include "test_suite_extras_def.inc"

  ! Procedures

contains

  include "test_suite_extras.inc"

  subroutine extras_controls(id, ierr)

    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! Set up hooks

    s% extras_startup => extras_startup
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve

    s% job% warn_run_star_extras =.false.

  end subroutine extras_controls

  !****

  subroutine extras_startup(id, restart, ierr)

    integer, intent(in)  :: id
    logical, intent(in)  :: restart
    integer, intent(out) :: ierr

    type (star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    call test_suite_startup(s, restart, ierr)

    if (s% job% create_pre_main_sequence_model) return

    ! Initialize GYRE

    call init('gyre.in')

    ! Set constants

    call set_constant('G_GRAVITY', standard_cgrav)
    call set_constant('C_LIGHT', clight)
    call set_constant('A_RADIATION', crad)

    call set_constant('M_SUN', Msun)
    call set_constant('R_SUN', Rsun)
    call set_constant('L_SUN', Lsun)

    call set_constant('GYRE_DIR', TRIM(mesa_dir)//'/build/gyre/src')

  end subroutine extras_startup

  !****

  include 'gyre_in_mesa_extras_finish_step.inc'

  integer function extras_finish_step(id)
    integer, intent(in) :: id
    integer :: ierr

    type (star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    extras_finish_step = keep_going

    if (s% job% create_pre_main_sequence_model) return

    extras_finish_step = gyre_in_mesa_extras_finish_step(id)

    if (extras_finish_step == terminate) &
       s% termination_code = t_extras_finish_step

  end function extras_finish_step

  !****

  subroutine extras_after_evolve(id, ierr)

    integer, intent(in)  :: id
    integer, intent(out) :: ierr

    type (star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    call test_suite_after_evolve(s, ierr)

    if (s% job% create_pre_main_sequence_model) return

    ! Finalize GYRE

    call final()

  end subroutine extras_after_evolve

end module run_star_extras
