! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

module pgstar_decorator

   ! NOTE: remember to set X_use_decorator = .true. to enable this,
   ! where X is the name of the pgstar plot
   ! and set s% pg% X_pgstar_decorator => your_function in your
   ! run_star_extras.f

   ! List of pgplot routines: http://www.astro.caltech.edu/~tjp/pgplot/annlist.html

   implicit none

contains

   ! default does nothing
   ! xmin, xmax, ymin, ymax: current plot boundary
   ! plot_num: If a plot has multiple sub-panels, then this tells you which panel is being called
   subroutine null_pgstar_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
      use star_def
      use const_def, only: dp
      integer, intent(in) :: id
      !Not doubles
      real, intent(in) :: xmin, xmax, ymin, ymax
      integer, intent(in) :: plot_num
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_pgstar_decorator

! Example function to add squares and some text to the abundance plot
!      subroutine Abundance_pgstar_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
!         use star_def
!         use const_def, only: dp
!         use pgstar_colors
!         integer, intent(in) :: id
!         !Not dp
!         real,intent(in) :: xmin, xmax, ymin, ymax
!         real :: xcenter,ycenter,dx,dy,a
!         integer, intent(out) :: ierr
!         integer :: i
!         type (star_info), pointer :: s

!         ierr = 0
!         call star_ptr(id, s, ierr)
!         if (ierr /= 0) return

!         dx=(xmax-xmin)
!         dy=(ymax-ymin)

!         xcenter=xmin+dx*0.5
!         ycenter=ymin+dy*0.5

!         call pgsci(clr_Coral)

!         do i=1,4
!            a=(i/10.0)
!            call pgline(5, (/xcenter-a*dx,xcenter-a*dx,xcenter+a*dx,xcenter+a*dx,xcenter-a*dx/),&
!                           (/ycenter-a*dy,ycenter+a*dy,ycenter+a*dy,ycenter-a*dy,ycenter-a*dy/))
!         end do

!         call pgptxt(xcenter,ycenter, 0.0, 1.0, 'Some added text on this plot')

!      end subroutine Abundance_pgstar_decorator

end module pgstar_decorator

