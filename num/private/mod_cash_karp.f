! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module mod_cash_karp
      use const_def, only: dp
      use math_lib
      
      
      contains

      subroutine do_cash_karp(
     &         n,fcn,x0,y0,xend,h0,max_step_size_in,max_steps,
     &         rtol,atol,itol,solout,iout,work,lwork,iwork,liwork,
     &         lrpar,rpar,lipar,ipar,lout,idid)
      integer, intent(in) :: n ! the dimension of the system
      interface
#include "num_fcn.dek"
#include "num_solout.dek"
      end interface
      real(dp), intent(inout) :: x0
      real(dp), intent(inout), pointer :: y0(:) ! (n)
      real(dp), intent(in) :: xend, h0, max_step_size_in
      real(dp), intent(in) :: rtol(*)
      real(dp), intent(in) :: atol(*)
      integer, intent(in) :: itol, iout, liwork, lwork, max_steps
      integer, intent(inout), pointer :: iwork(:) ! (liwork)
      real(dp), intent(inout), pointer :: work(:) ! (lwork)
      integer, intent(in) :: lrpar, lipar
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      integer, intent(in) :: lout
      integer, intent(out) :: idid

      real(dp), dimension(:), pointer :: k1, k2, k3, k4, k5, k6,
     >   y, y6, y5, yd, dydx, err, tol
      real(dp) :: x, h, maxerr, delta, rwork_y(0), max_step_size
      real(dp), parameter :: SAFETY=0.9D0 ! for Runge Kutta Fehlberg tolerance check
      integer :: nstep, i, nfcn, naccpt, nrejct, irtrn, ierr, iwork_y(0)
      logical :: increasing
      
      include 'formats'
      
      idid = 0
      call setup
      if (idid /= 0) return
      
      do i=1,n
         y(i) = y0(i)
      end do
      x = x0
      h = h0
      max_step_size = abs(max_step_size_in)
      ierr = 0
      idid = -2
      naccpt = 0
      nrejct = 0
      nfcn = 0
      increasing = (xend > x0)
      if ((.not. increasing) .and. h > 0) h = -h ! want negative h for decreasing x
      
      do nstep = 0, max_steps
         
         if (increasing) then
            if (x >= xend) then
               idid = 1
               exit
            end if
         else ! decreasing x
            if (x <= xend) then
               idid = 1
               exit
            end if
         end if

         call fcn(n, x, h, y, dydx, lrpar, rpar, lipar, ipar, ierr)
         nfcn = nfcn+1
         do i=1,n
            k1(i) = h*dydx(i)
         end do

         if (ierr == 0) then
            do i=1,n
               yd(i) = y(i) + k1(i)/5d0
            end do
            call fcn(n, x + h/5d0, h, yd, dydx, lrpar, rpar, lipar, ipar, ierr)
            nfcn = nfcn+1
            do i=1,n
               k2(i) = h*dydx(i)
            end do
         end if

         if (ierr == 0) then
            do i=1,n
               yd(i) = y(i) + (3d0*k1(i) + 9d0*k2(i))/40d0
            end do
            call fcn(n, x + 3d0*h/10d0, h, yd, dydx, lrpar, rpar, lipar, ipar, ierr)
            nfcn = nfcn+1
            do i=1,n
               k3(i) = h*dydx(i)
            end do
         end if

         if (ierr == 0) then
            do i=1,n
               yd(i) = y(i) + (3d0*k1(i) - 9d0*k2(i) + 12d0*k3(i))/10d0
            end do
            call fcn(n, x + 3d0*h/5d0, h, yd, dydx, lrpar, rpar, lipar, ipar, ierr)
            nfcn = nfcn+1
            do i=1,n
               k4(i) = h*dydx(i)
            end do
         end if

         if (ierr == 0) then
            do i=1,n
               yd(i) = y(i) -11d0*k1(i)/54d0 + 5d0*k2(i)/2d0 - 70d0*k3(i)/27d0 + 35d0*k4(i)/27d0
            end do
            call fcn(n, x + h, h, yd, dydx, lrpar, rpar, lipar, ipar, ierr)
            do i=1,n
               k5(i) = h*dydx(i)
            end do
         end if

         if (ierr == 0) then
            do i=1,n
               yd(i) = y(i) + 1631d0*k1(i)/55296d0 + 175d0*k2(i)/512d0 + 
     >            575d0*k3(i)/13824d0 + 44275d0*k4(i)/110592d0 + 253d0*k5(i)/4096d0
            end do
            call fcn(n, x + 7d0*h/8d0, h, yd, dydx, lrpar, rpar, lipar, ipar, ierr)
            nfcn = nfcn+1
            do i=1,n
               k6(i) = h*dydx(i)
            end do
         end if
         
         ! if ierr /= 0 then reduce stepsize and try again

         maxerr=0d0
         if (ierr == 0) then
            do i=1,n
               y5(i) = y(i) + 2825d0*k1(i)/27648d0 + 18575d0*k3(i)/48384d0 + 
     >            13525d0*k4(i)/55296d0 + 277d0*k5(i)/14336d0 + k6(i)/4d0
               y6(i) = y(i) + 37*k1(i)/378d0 + 250d0*k3(i)/621d0 + 
     >            125d0*k4(i)/594d0  + 512d0*k6(i)/1771d0
               err(i) = abs(y6(i) - y5(i))
            end do
            if (itol == 0) then
               do i=1,n
                  tol(i)= abs(y5(i))*rtol(1) +  atol(1)
               end do
            else
               do i=1,n 
                  tol(i)= abs(y5(i))*rtol(i) +  atol(i)
               end do
            end if
            maxerr = 0
            do i=1,n
               if (err(i)/tol(i) > maxerr) maxerr = err(i)/tol(i)
            end do
         end if

         if (ierr == 0 .and. maxerr <= 1) then ! okay.  accept step.
         
            naccpt = naccpt + 1
            do i=1,n
               y(i) = y5(i)
            end do
            x = x + h
            if (iout > 0) then
               call solout(naccpt, x-h, x, n, y, rwork_y, iwork_y, 
     >               null_interp_y, lrpar, rpar, lipar, ipar, irtrn)
               if (irtrn < 0) then
                  idid = 2
                  exit
               end if
            end if
            
            delta = SAFETY*pow(maxerr,-0.2d0)
            if (delta >= 4d0) then ! max increase is factor of 4
               h = h*4d0
            else if (delta > 1) then
               h = h*delta
            end if
         
            if (increasing) then
               if (x+h > xend) h = xend-x
               if (h > max_step_size .and. max_step_size > 0) h = max_step_size
            else ! decreasing x
               if (x+h < xend) h = xend-x
               if (-h > max_step_size .and. max_step_size > 0) h = -max_step_size
            end if
            
         else ! reject and retry with smaller step
         
            nrejct = nrejct + 1
            delta = SAFETY*pow(maxerr,-0.25d0)
            if (delta < 0.1d0) then
               h = h*0.1d0
            else 
               h = h*delta
            end if
            ierr = 0
            
         end if
         
      end do
      
      do i=1,n
         y0(i) = y(i)
      end do
      x0 = x

      ! statistics returned in iwork
      iwork(1) = nfcn
      iwork(2) = nstep
      iwork(3) = naccpt
      iwork(4) = nrejct
      

      contains
      
      subroutine setup
         integer :: i
         i = 0
         k1 => work(i+1:i+n); i = i+n
         k2 => work(i+1:i+n); i = i+n
         k3 => work(i+1:i+n); i = i+n
         k4 => work(i+1:i+n); i = i+n
         k5 => work(i+1:i+n); i = i+n
         k6 => work(i+1:i+n); i = i+n
         y => work(i+1:i+n); i = i+n
         y6 => work(i+1:i+n); i = i+n
         y5 => work(i+1:i+n); i = i+n
         yd => work(i+1:i+n); i = i+n
         dydx => work(i+1:i+n); i = i+n
         err => work(i+1:i+n); i = i+n
         tol => work(i+1:i+n); i = i+n
         if (i > lwork) then
            idid = -1
            if (lout > 0) write(lout,*) ' insufficient storage for work, min. lwork=',i
         end if
         if (liwork < 4) then
            idid = -1
            if (lout > 0) write(lout,*) ' insufficient storage for iwork, min. liwork=',4
         end if
      end subroutine setup
      
      end subroutine do_cash_karp
      

      real(dp) function null_interp_y(i, s, rwork_y, iwork_y, ierr)
         integer, intent(in) :: i
         real(dp), intent(in) :: s
         real(dp), intent(inout), target :: rwork_y(*)
         integer, intent(inout), target :: iwork_y(*)
         integer, intent(out) :: ierr
         ierr = -1
         null_interp_y = 0
         write(*,*) 'sorry: cash_karp integrator does not support interpolation'
      end function null_interp_y



      end module mod_cash_karp


