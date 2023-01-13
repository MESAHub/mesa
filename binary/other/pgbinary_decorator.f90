! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module pgbinary_decorator

   ! NOTE: remember to set X_use_decorator = .true. to enable this, 
   ! where X is the name of the pgbinary plot
   ! and set s% X_pgbinary_decorator => your_function in your
   ! run_binary_extras.f

   ! List of pgplot routines: http://www.astro.caltech.edu/~tjp/pgplot/annlist.html

   use binary_def
   use const_def

   implicit none


contains


   ! default does nothing
   ! xmin, xmax, ymin, ymax: current plot boundary
   ! plot_num: If a plot has multiple sub-panels, then this tells you which panel is being called
   subroutine null_pgbinary_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
      integer, intent(in) :: id
      !Not doubles
      real, intent(in) :: xmin, xmax, ymin, ymax
      integer, intent(in) :: plot_num
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_pgbinary_decorator

end module pgbinary_decorator
      
      
      
      
