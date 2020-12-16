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


      module mod_root 
      use const_def, only: dp, arg_not_provided
      
      
      implicit none


      contains


      real(dp) function do_safe_root_with_brackets(f,x1_in,x3_in,y1_in,y3_in,
     >        imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr)
         use mod_brent, only: eval_brent_safe_zero
         interface
! f provides function values
#include "num_root_fcn.dek" 
         end interface
         real(dp), intent(in) :: x1_in, x3_in 
         real(dp), intent(in) :: y1_in, y3_in ! f(x1) and f(x3)
         integer, intent(in) :: lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: imax
         real(dp), intent(in) :: epsx, epsy
         integer, intent(out) :: ierr
         
         do_safe_root_with_brackets = eval_brent_safe_zero( 
     >         x1_in, x3_in, 1d-14, epsx, epsy, f, y1_in, y3_in, 
     >         lrpar, rpar, lipar, ipar, ierr )

      end function do_safe_root_with_brackets

      
      real(dp) function do_safe_root_with_guess(
     >      f, x_guess, dx, x1_in, x3_in, y1_in, y3_in, newt_imax, imax, 
     >      epsx, epsy, lrpar, rpar, lipar, ipar, ierr)
         interface
! f provides function values and optional derivatives
#include "num_root_fcn.dek" 
         end interface
         real(dp), intent(in) :: x_guess ! initial guess for the root (required)
         real(dp), intent(in) :: dx, x1_in, x3_in, y1_in, y3_in
            ! dx is increment for searching for brackets (not used if x1 and x3 are given)
            ! x1 and x3 bracket solution (can be arg_not_provided)
            ! y1 and y3 are f(x1) and f(x3) (can be arg_not_provided)
         integer, intent(in) :: lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: newt_imax, imax
         real(dp), intent(in) :: epsx, epsy 
         ! stop seaching when x is determined to within epsx
         ! or when abs(f(x)) is less than epsy
         integer, intent(out) :: ierr
         
         integer :: i
         logical :: have_x1, have_x3, have_y1, have_y3
         real(dp) :: x, y, x1, x3, y1, y3, dydx, absy, absy_prev
         
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         
         ierr = 0
         do_safe_root_with_guess = 0
         x1 = x1_in
         x3 = x3_in
         y1 = y1_in
         y3 = y3_in
         have_x1 = (x1 /= arg_not_provided)
         have_x3 = (x3 /= arg_not_provided)
         have_y1 = (y1 /= arg_not_provided)
         have_y3 = (y3 /= arg_not_provided)
         x = x_guess
         absy_prev = 0
         
         do i=1, newt_imax
            y = f(x,dydx,lrpar,rpar,lipar,ipar,ierr)
            if (ierr /= 0) then
               if (dbg) then
                  write(*,3) 'x y dydx', ierr, i, x, y, dydx
                  stop
               end if
               ierr = 0; exit ! try safe_root
            end if
            if (dbg) then
               write(*,2) 'x y dydx', i, x, y, dydx
            end if
            absy = abs(y)
            ! converged if abs(y) < epsy or abs(y/dydx) < epsx
            !write(*,2) 'y, epsy, y/dydx, epsx', i, y, epsy, y/dydx, epsx
            if (absy < max(epsy, abs(dydx)*epsx)) then
               if (dbg) write(*,1) 'converged', x
               do_safe_root_with_guess = x
               return
            end if
            if (dydx == 0) exit ! try safe_root
            if (i > 1 .and. absy > absy_prev) exit ! not converging
            absy_prev = absy
            if (y < 0) then
               x1 = x; y1 = y; have_x1 = .true.
            else
               x3 = x; y3 = y; have_x3 = .true.
            end if
            x = x - y/dydx
         end do
         
         if (.not. (have_x1 .and. have_x3)) then
            call do_look_for_brackets( 
     >         x_guess,dx,x1,x3,f,y1,y3,imax,lrpar,rpar,lipar,ipar,ierr)
            if (ierr /= 0) then
               if (dbg) then
                  write(*,*) 'failed in do_look_for_brackets'
                  stop 'do_safe_root_with_guess'
               end if
               return
            end if
         else
            if (.not. have_y1) then
               y1 = f(x1,dydx,lrpar,rpar,lipar,ipar,ierr)
               if (ierr /= 0) then
                  if (dbg) then
                     write(*,*) 'failed evaluating f(x1)'
                     stop 'do_safe_root_with_guess'
                  end if
                  return
               end if
            end if
            if (.not. have_y3) then
               y3 = f(x3,dydx,lrpar,rpar,lipar,ipar,ierr)
               if (ierr /= 0) then
                  if (dbg) then
                     write(*,*) 'failed evaluating f(x3)'
                     stop 'do_safe_root_with_guess'
                  end if
                  return
               end if
            end if
         end if
         
         do_safe_root_with_guess =  
     >      do_safe_root_with_brackets( 
     >         f,x1,x3,y1,y3,imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr) 
         if (ierr /= 0) then
            if (dbg) then
               write(*,*) 'failed in do_safe_root'
               stop 'do_safe_root_with_guess'
            end if
            return
         end if

         if (dbg) then
            write(*,1) 'do_safe_root result:', do_safe_root_with_guess
            write(*,*)
         end if
      
      end function do_safe_root_with_guess

      
      
      ! safe_root requires bracketing values for the root.
      ! if you don't have them, you can use this routine to do a (not too dumb) search.
      subroutine do_look_for_brackets(x,dx,x1,x3,f,y1,y3,imax,lrpar,rpar,lipar,ipar,ierr)
         real(dp), intent(in) :: x, dx ! x is initial guess and dx is increment for searching
         real(dp), intent(out) :: x1, x3 ! bounds
         real(dp), intent(out) :: y1, y3 ! f(x1) and f(x3)
         interface
! f provides function values
#include "num_root_fcn.dek" 
         end interface
         integer, intent(in) :: imax
         integer, intent(in) :: lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         
         real(dp) :: jump, dfdx
         integer :: i
         logical :: move_x1, move_x3
         
         include 'formats.dek'

         ierr = -1
         y1 = 0
         y3 = 0
         
         x1 = x; x3 = x
         ! after this, keep x1 < x3
         
         !write(*,2) 'do_look_for_brackets imax x dx', imax, x, dx

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
               y1 = f(x1,dfdx,lrpar,rpar,lipar,ipar,ierr)
               if (ierr /= 0) then
                  !write(*,*) 'ierr from f(x1)'
                  return
               end if
            end if
            
            if (move_x3) then
               y3 = f(x3,dfdx,lrpar,rpar,lipar,ipar,ierr)
               if (ierr /= 0) then
                  !write(*,*) 'ierr from f(x3)'
                  return
               end if
            end if
            
            !write(*,'(a,i4,4(3x,e18.8))') 'look_for_brackets', i, x1, y1, x3, y3
            
            if (y1*y3 <= 0) then
               ierr = 0
               !write(*,1) 'done do_look_for_brackets', y1, y3
               return
            end if
         
         end do
         
         !write(*,1) 'exit do_look_for_brackets'
      
      end subroutine do_look_for_brackets


      end module mod_root
      
      
      
      
      
