! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
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

module mod_integrate
   use const_def
   use math_lib
   use num_def
   use utils_lib, only : mesa_error
   
   implicit none

contains
   
   recursive function integrator(func, minx, maxx, args, atol, rtol, max_steps, ierr) result(res)
      procedure(integrand) :: func
      real(dp), intent(in) :: minx, maxx ! Min and max values to integrate over
      real(dp), intent(in) :: args(:) ! Extra args passed to func
      real(dp), intent(in) :: atol, rtol ! Absolute and relative tolerances
      integer, intent(in) :: max_steps ! Max number of sub-steps
      integer, intent(inout) :: ierr ! Error code
      real(dp) :: res
      
      real(dp) :: val1, val2
      real(dp) :: xlow, xhigh, xmid
      
      if(max_steps < 1) then
         ierr = -1
         return
      end if
      
      ierr = 0
      
      xlow = minx
      xhigh = maxx
      xmid = (xhigh + xlow) / 2.d0
      
      if(xhigh < xlow) then
         ierr = -1
         return
      end if
      
      val1 = simp38(func, xlow, xhigh, args, ierr)
      if(ierr/=0) return
      
      val2 = simp38(func, xlow, xmid, args, ierr) + simp38(func, xmid, xhigh, args, ierr)
      if(ierr/=0) return
      
      if(val1==0d0 .or. val2 == 0d0) then
         res = val2
         return
      end if
      
      if(abs(val1 - val2) < atol .or. abs(val1 - val2) / val1 < rtol) then
         res = val2
      else
         val1 = integrator(func, xlow, xmid, args, atol, rtol, max_steps - 1, ierr)
         val2 = integrator(func, xmid, xhigh, args, atol, rtol, max_steps - 1, ierr)
         
         res = val1 + val2
         if(ierr/=0) return
      end if
   
   end function integrator
   
   
   real(dp) function simp38(func, minx, maxx, args, ierr)
      ! Simpsons 3/8 rule
      procedure(integrand) :: func
      real(dp), intent(in) :: minx, maxx ! Min and max values to integrate over
      real(dp), intent(in) :: args(:) ! Extra args passed to func
      integer, intent(inout) :: ierr
      
      real(dp) :: x
      
      ierr = 0
      
      x = minx
      simp38 = func(x, args, ierr)
      if(ierr/=0) return
      
      x = (2 * minx + maxx) / 3.d0
      simp38 = simp38 + 3 * func(x, args, ierr)
      if(ierr/=0) return
      
      x = (minx + 2 * maxx) / 3.d0
      simp38 = simp38 + 3 * func(x, args, ierr)
      if(ierr/=0) return
      
      x = maxx
      simp38 = simp38 + func(x, args, ierr)
      if(ierr/=0) return
      
      simp38 = (maxx - minx) / 8.d0 * simp38
   
   end function simp38


end module mod_integrate
