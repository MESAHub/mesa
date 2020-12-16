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


      module mod_solout_root 
      use const_def, only: dp
      
      
      implicit none


      contains


      real(dp) function do_solout_root(
     >         f,x_guess,x1_in,x3_in,y1_in,y3_in,imax,epsx,epsy,
     >         rwork_y,iwork_y,interp_y,lrpar,rpar,lipar,ipar,ierr)
         use const_def, only: arg_not_provided
         interface
#include "num_solout_root_fcn.dek"
         end interface
         real(dp), intent(in) :: x_guess, x1_in, x3_in 
         real(dp), intent(in) :: y1_in, y3_in ! f(x1) and f(x3)
         integer, intent(in) :: lipar, lrpar
         real(dp), intent(inout), target :: rwork_y(*)
         integer, intent(inout), target :: iwork_y(*)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         interface
#include "num_interp_y.dek"
         end interface
         integer, intent(in) :: imax
         real(dp), intent(in) :: epsx, epsy
         integer, intent(out) :: ierr
         
         real(dp) :: x1,x2,x3,y1,y2,y3,xm,ym,tmp,y21,y32,y31,b,c,dydx,dx,dxold,xfirst
         integer :: i
         logical :: converged, try_xm
         
         do_solout_root = 0
         
         x1 = x1_in
         x3 = x3_in
         y1 = y1_in
         y3 = y3_in

         dx = abs(x3-x1)
         dydx = 0
         converged = .false.
         ierr = 0
         
         if (y1 == 0) then
            do_solout_root = x1; return
         end if
         
         if (y3 == 0) then
            do_solout_root = x3; return
         end if
         
         if (y1*y3 > 0) then
            do_solout_root = 0
            ierr = -1
            return
         end if
         
         if (x_guess == arg_not_provided) then
            xfirst = (x1+x3)/2
         else
            xfirst = x_guess
         end if
         
         do i=1,imax
            
            if (i == 1) then
               x2 = xfirst
            else
               x2 = (x1+x3)/2
            end if
            dxold = dx
            dx = abs(x3-x1) 
               ! don't check dx until after evaluate y
               ! so that when return, will have just evaluated y for the accepted root.
            y2 = get_y(x2,dydx)
            if (ierr /= 0) return
            if (converged .or.dx <= epsx) then
               do_solout_root = x2; return
            end if
         
            if (y1*y2 > 0) then ! exchange points 1 and 3
               tmp = x1; x1 = x3; x3 = tmp
               tmp = y1; y1 = y3; y3 = tmp
            end if
            
            ! at this point, y1 and y2 are opposite sign
            ! so root is between x1 and x2
            
            try_xm = .false.
            xm = x2; ym = y2
            do while (dydx /= 0 .and. abs(2*ym) < abs(dxold*dydx))
               ! try newton
               xm = xm - ym/dydx 
               try_xm = ((xm-x1)*(xm-x2) < 0) ! xm in bounds
               if (.not. try_xm) exit
               call eval_xm
               if (converged .or. ierr /= 0) return
            end do
         
            if (.not. try_xm) then ! try parabolic interpolation
               y21 = y2-y1
               y32 = y3-y2
               if (y32 /= 0) then
                  y31 = y3-y1
                  b = (x2-x1)/y21
                  c = (y21-y32)/(y32*y31)
                  if (y3*y31 >= 2*y2*y21) then
                     xm = x1-b*y1*(1-c*y2)
                     try_xm = .true.
                  end if
               end if
            end if
            
            if (try_xm) then
               call eval_xm
               if (converged .or. ierr /= 0) return               
            end if
            
            x3 = x2
            y3 = y2
         
         end do
         
         ierr = -1
         do_solout_root = 0
         
         contains
         
         subroutine eval_xm
            ym = get_y(xm,dydx)
            if (ierr /= 0) return
            if (converged) then
               do_solout_root = xm; return
            end if
            if (ym*y1 < 0) then
               x2 = xm; y2 = ym
            else
               x1 = xm; y1 = ym
            end if
         end subroutine eval_xm
         
         real(dp) function get_y(x,dydx)
            real(dp), intent(in) :: x
            real(dp), intent(out) :: dydx
            get_y = f(x, dydx,rwork_y,iwork_y,interp_y,lrpar,rpar,lipar,ipar,ierr)
            converged = (abs(get_y) <= epsy)
         end function get_y

      end function do_solout_root
      
      
      subroutine do_solout_brackets(
     >         x,dx,x1,x3,f,y1,y3,imax,rwork_y,iwork_y,interp_y,lrpar,rpar,lipar,ipar,ierr)
         real(dp), intent(in) :: x, dx ! x is initial guess and dx is increment for searching
         real(dp), intent(out) :: x1, x3 ! bounds
         real(dp), intent(out) :: y1, y3 ! f(x1) and f(x3)
         interface
#include "num_solout_root_fcn.dek"
         end interface
         integer, intent(in) :: imax
         integer, intent(in) :: lipar, lrpar
         real(dp), intent(inout), target :: rwork_y(*)
         integer, intent(inout), target :: iwork_y(*)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         interface
#include "num_interp_y.dek"
         end interface
         integer, intent(out) :: ierr
         
         real(dp) :: jump, dfdx
         integer :: i
         logical :: move_x1, move_x3

         ierr = -1
         y1 = 0
         y3 = 0
         
         x1 = x; x3 = x
         ! after this, keep x1 < x3

         do i = 1, imax
            
            jump = (2**i)*dx
            
            if (i == 1) then
               x1 = x1 - jump
               x3 = x3 + jump
               move_x1 = .true.
               move_x3 = .true.
            else if (y1 > 0) then ! both positive. move x for smaller
               if (y1 < y3) then
                  x1 = x1 - jump
                  move_x1 = .true.
                  move_x3 = .false.
               else
                  x3 = x3 + jump
                  move_x3 = .true.
                  move_x1 = .false.
               end if
            else ! both negative. move x for larger
               if (y1 > y3) then
                  x1 = x1 - jump
                  move_x1 = .true.
                  move_x3 = .false.
               else
                  x3 = x3 + jump
                  move_x3 = .true.
                  move_x1 = .false.
               end if
            end if
            
            if (move_x1) then
               y1 = f(x1,dfdx,rwork_y,iwork_y,interp_y,lrpar,rpar,lipar,ipar,ierr)
               if (ierr /= 0) return
            end if
            
            if (move_x3) then
               y3 = f(x3,dfdx,rwork_y,iwork_y,interp_y,lrpar,rpar,lipar,ipar,ierr)
               if (ierr /= 0) return
            end if
            
            if (y1*y3 <= 0) then
               ierr = 0; return
            end if
         
         end do
      
      end subroutine do_solout_brackets
      

      end module mod_solout_root
      
      
      
      
      
