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
   
   use star_lib
   use star_def
   use const_def
   use gyre_lib
   use math_lib
   
   implicit none
   
   include "test_suite_extras_def.inc"

contains

include "test_suite_extras.inc"

include 'run_star_extras.inc'
      
      
      subroutine extras_controls(id, ierr)
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   
   write(*, *) 'cannot run rsp without gyre.'
   write(*, *) 'this test was intentionally skipped'
   write(*, *) 'good match for period', -1d0, -1d0
   
   open(unit = 30, file = 'final.mod', action ='write', status = 'replace')
   write(30, *) 'fake final.mod'
   close(30)
   
   ierr = -1
   return
   
   end subroutine extras_controls


end module run_star_extras
      
