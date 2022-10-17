! ***********************************************************************
!
!   Copyright (C) 2013  The MESA Team
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
   
   use star_lib
   use star_def
   
   implicit none
   
   logical, parameter :: GYRE_IS_ENABLED = .false.


contains
   
   
   subroutine init_gyre(gyre_file, ierr)
      character(*), intent(in) :: gyre_file
      integer, intent(out) :: ierr
      ierr = -1
   end subroutine init_gyre
   
   
   subroutine do_gyre_get_modes (s, el, store_model, ierr)
      
      type (star_info), pointer :: s
      integer, intent(in) :: el
      logical, intent(in) :: store_model
      integer, intent(out) :: ierr
      
      ierr = -1
   end subroutine do_gyre_get_modes


end module gyre_support
