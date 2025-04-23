! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

program run

   use run_star_support, only: do_read_star_job
   use run_star, only: do_run_star
   use star_def

   implicit none

   integer :: ierr

   ierr = 0
   call do_read_star_job('inlist', ierr)
   if (ierr /= 0) stop 'run'

   call do_run_star

end program run
