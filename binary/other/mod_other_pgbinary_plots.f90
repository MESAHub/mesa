! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton & The MESA Team
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

module mod_other_pgbinary_plots

   ! consult star/other/README.md for general usage instructions
   ! control name: use_other_pgbinary_plots = .true.
   ! procedure pointer: s% other_pgbinary_plots_info => my_routine

   ! see star/other/sample_pgbinary_plot.f90

   implicit none

contains

   ! default does nothing
   subroutine null_other_pgbinary_plots_info(id, ierr)
      use binary_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_other_pgbinary_plots_info

end module mod_other_pgbinary_plots
