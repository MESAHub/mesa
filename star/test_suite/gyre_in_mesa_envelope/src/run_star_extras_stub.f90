! ***********************************************************************
!
!   Copyright (C) 2018-2019  Rich Townsend, Bill Paxton & The MESA Team
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

    write(*,*) 'GYRE is not installed: this test will run without '
    write(*,*) 'calling GYRE and pretend to pass'

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

  end subroutine extras_startup

  !****

  integer function extras_finish_step(id)
    integer, intent(in) :: id
    integer :: ierr

    type (star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    extras_finish_step = keep_going

  end function extras_finish_step

  !****

  subroutine extras_after_evolve(id, ierr)

    integer, intent(in)  :: id
    integer, intent(out) :: ierr

    type (star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    write(*,*) 'matched target'

    call test_suite_after_evolve(s, ierr)

  end subroutine extras_after_evolve

end module run_star_extras
