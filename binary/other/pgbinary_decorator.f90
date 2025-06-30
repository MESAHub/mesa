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

module pgbinary_decorator

   ! NOTE: remember to set X_use_decorator = .true. to enable this,
   ! where X is the name of the pgbinary plot
   ! and set s% X_pgbinary_decorator => your_function in your
   ! run_binary_extras.f

   ! List of pgplot routines: http://www.astro.caltech.edu/~tjp/pgplot/annlist.html

   implicit none

contains

   ! default does nothing
   ! xmin, xmax, ymin, ymax: current plot boundary
   ! plot_num: If a plot has multiple sub-panels, then this tells you which panel is being called
   subroutine null_pgbinary_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
      use binary_def
      use const_def, only: dp
      integer, intent(in) :: id
      !Not doubles
      real, intent(in) :: xmin, xmax, ymin, ymax
      integer, intent(in) :: plot_num
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_pgbinary_decorator

end module pgbinary_decorator
