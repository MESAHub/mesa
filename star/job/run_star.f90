! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module run_star
      implicit none

      contains

      subroutine do_run_star(inlist_fname_arg)
         use run_star_support, only: run1_star
         use run_star_extras, only: extras_controls
         use star_lib, only: starlib_shutdown
         character (len=*) :: inlist_fname_arg
         optional inlist_fname_arg
         logical, parameter :: &
            do_alloc_star = .true., &
            do_free_star = .true., &
            restart_okay = .true.
         integer :: id, ierr
         logical :: restart
         call run1_star( &
            do_alloc_star, &
            do_free_star, &
            restart_okay, &
            id, restart, &
            extras_controls, &
            ierr, &
            inlist_fname_arg)
         call starlib_shutdown
      end subroutine do_run_star

      end module run_star
