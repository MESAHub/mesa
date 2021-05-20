! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module interp_1d_misc_sg

      implicit none
      
      contains
      
      subroutine do_integrate_values_sg(init_x, nx, f1, nv, x, vals, ierr)
         real, intent(in) :: init_x(:) ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx ! length of init_x vector
         real, intent(in), pointer :: f1(:) ! =(4, nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv ! length of new x vector and vals vector
         real, intent(in) :: x(:) ! (nv)
            ! strictly monotonic in same way as init_x
         real, intent(inout) :: vals(:) ! (nv)
            ! for i > 1, vals(i) = integral of interpolating poly from x(i-1) to x(i)
            ! vals(1) = 0
         integer, intent(out) :: ierr ! 0 means aok
   
         integer :: k_old, k_new
         real :: xk_old, xkp1_old, xk_new, xk_prev, sum
         logical :: increasing              
         real, pointer :: f(:,:) ! (4, nx)  ! data & interpolation coefficients
         f(1:4,1:nx) => f1(1:4*nx)

         increasing = (init_x(1) < init_x(nx))
         
         if (increasing .and. (x(1) < init_x(1) .or. x(nv) > init_x(nx)) &
             .or. ((.not. increasing) .and. (x(1) > init_x(1) .or. x(nv) < init_x(nx)))) then
            ierr = -1
            return
         end if
         
         ierr = 0
         
         k_old = 1; xk_old = init_x(k_old); xkp1_old = init_x(k_old+1)
         sum = 0; xk_prev = x(1); vals(1) = 0
         
         do k_new = 2, nv
         
            xk_new = x(k_new)
            do while ((increasing .and. xk_new > xkp1_old) .or. ((.not. increasing) .and. xk_new < xkp1_old))
               k_old = k_old + 1
               if (k_old >= nx) then
                  k_old = k_old - 1
                  xk_new = xkp1_old
                  exit
               end if
               call add_to_integral(k_old - 1, xkp1_old)
               xk_old = xkp1_old
               xkp1_old = init_x(k_old+1)
            end do
            
            call add_to_integral(k_old, xk_new)
            vals(k_new) = sum
            sum = 0
            
         end do
         
         contains
         
         subroutine add_to_integral(k, x2)
            integer, intent(in) :: k
            real, intent(in) :: x2
            
            real :: x0, x1, a1, a2, d1, d2, area
            
            x0 = init_x(k)
            x1 = xk_prev
            if (x1 == x2) return
            d2 = x2 - x0
            a2 = d2*(f(1, k) + d2*(f(2, k)/2  &
                     + d2*(f(3, k)/3 + d2*f(4, k)/4)))
            if (x1 > x0) then
               d1 = x1 - x0
               a1 = d1*(f(1, k) + d1*(f(2, k)/2  &
                           + d1*(f(3, k)/3 + d1*f(4, k)/4)))
               area = a2 - a1
            else
               d1 = 0; a1 = 0; area = a2
            end if
            sum = sum + area
            xk_prev = x2
         
         end subroutine add_to_integral
                     
      
      end subroutine do_integrate_values_sg
      
      
      subroutine do_interp_values_sg(init_x, nx, f1, nv, x, vals, ierr)
         real, intent(in) :: init_x(:) ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx ! length of init_x vector
         real, intent(in), pointer :: f1(:) ! =(4, nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv ! length of new x vector and vals vector
         real, intent(in) :: x(:) ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         real, intent(inout) :: vals(:) ! (nv)
         integer, intent(out) :: ierr ! 0 means aok
   
         integer :: k_old, k_new
         real :: xk_old, xkp1_old, xk_new, delta
         logical :: increasing
         real, pointer :: f(:,:) ! (4, nx)  ! data & interpolation coefficients
         f(1:4,1:nx) => f1(1:4*nx)
         
         ierr = 0      
         
         if (nx == 1) then
            vals(1:nv) = f1(1)
            return
         end if

         f(1:4,1:nx) => f1(1:4*nx)

         increasing = (init_x(1) < init_x(nx))
         
         k_old = 1; xk_old = init_x(k_old); xkp1_old = init_x(k_old+1)
                  
         do k_new = 1, nv
      
            xk_new = x(k_new)
            if (increasing) then
               if (xk_new > init_x(nx)) then
                  xk_new = init_x(nx)
               else if (xk_new < init_x(1)) then
                  xk_new = init_x(1)
               end if
            else ! decreasing
               if (xk_new < init_x(nx)) then
                  xk_new = init_x(nx)
               else if (xk_new > init_x(1)) then
                  xk_new = init_x(1)
               end if
            end if
            do while ((increasing .and. xk_new > xkp1_old) .or. ((.not. increasing) .and. xk_new < xkp1_old))
               k_old = k_old + 1
               if (k_old >= nx) then
                  k_old = k_old - 1
                  xk_new = xkp1_old
                  exit
               end if
               xk_old = xkp1_old
               xkp1_old = init_x(k_old+1)
            end do
         
            delta = xk_new - xk_old
         
            vals(k_new) =  &
                  f(1, k_old) + delta*(f(2, k_old)  &
                     + delta*(f(3, k_old) + delta*f(4, k_old)))
            
         end do

   
      end subroutine do_interp_values_sg
      
      
      subroutine do_interp_values_and_slopes_sg(init_x, nx, f1, nv, x, vals, slopes, ierr)
         real, intent(in) :: init_x(:) ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx ! length of init_x vector
         real, intent(in), pointer :: f1(:) ! =(4, nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv ! length of new x vector and vals vector
         real, intent(in) :: x(:) ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         real, intent(inout) :: vals(:) ! (nv)
         real, intent(inout) :: slopes(:) ! (nv)
         integer, intent(out) :: ierr ! 0 means aok
   
         integer :: k_old, k_new
         real :: xk_old, xkp1_old, xk_new, delta                  
         logical :: increasing   
         real, pointer :: f(:,:) ! (4, nx)  ! data & interpolation coefficients
         f(1:4,1:nx) => f1(1:4*nx)
         
         ierr = 0           
         
         k_old = 1; xk_old = init_x(k_old); xkp1_old = init_x(k_old+1)

         increasing = (init_x(1) < init_x(nx))
         
         do k_new = 1, nv
         
            xk_new = x(k_new)
            if (increasing) then
               if (xk_new > init_x(nx)) then
                  xk_new = init_x(nx)
               else if (xk_new < init_x(1)) then
                  xk_new = init_x(1)
               end if
            else ! decreasing
               if (xk_new < init_x(nx)) then
                  xk_new = init_x(nx)
               else if (xk_new > init_x(1)) then
                  xk_new = init_x(1)
               end if
            end if
            do while ((increasing .and. xk_new > xkp1_old) .or. ((.not. increasing) .and. xk_new < xkp1_old))
               k_old = k_old + 1
               if (k_old >= nx) then
                  k_old = k_old - 1
                  xk_new = xkp1_old
                  exit
               end if
               xk_old = xkp1_old
               xkp1_old = init_x(k_old+1)
            end do
            
            delta = xk_new - xk_old
            
            vals(k_new) =  &
                  f(1, k_old) + delta*(f(2, k_old)  &
                     + delta*(f(3, k_old) + delta*f(4, k_old)))
            
            slopes(k_new) =  &
                  f(2, k_old) + 2*delta*(f(3, k_old) + 1.5*delta*f(4, k_old))
            
         end do

   
      end subroutine do_interp_values_and_slopes_sg
      
            
      real function minmod1_sg(f1, f2)
         real, intent(in) :: f1, f2       
         minmod1_sg = 0.5 * (sign(1.0, f1) + sign(1.0, f2)) * min(abs(f1), abs(f2))    
      end function minmod1_sg
      
      
      real function median1_sg(f1, f2, f3)
         real, intent(in) :: f1, f2, f3
         median1_sg = f1 + minmod1_sg(f2 - f1, f3 - f1)
      end function median1_sg

      
      subroutine minmod_sg(z, n, f1, f2)
         real, intent(inout) :: z(:) ! (n)     
         integer, intent(in) :: n       ! length of vectors
         real, intent(in) :: f1(:), f2(:) ! (n)       
         z(1:n) = 0.5 * (sign(1.0, f1(1:n)) + sign(1.0, f2(1:n))) * min(abs(f1(1:n)), abs(f2(1:n)))      
      end subroutine minmod_sg

      
      subroutine minmod4_sg(z, n, f1, f2, f3, f4)
         real, intent(inout) :: z(:) ! (n)     
         integer, intent(in) :: n       ! length of vectors
         real, intent(in) :: f1(:), f2(:), f3(:), f4(:)
         call minmod_sg(z, n, f1, f2)
         call minmod_sg(z, n, z, f3)
         call minmod_sg(z, n, z, f4)
      end subroutine minmod4_sg
      
      
      subroutine median_sg(z, n, f1, f2, f3)
         real, intent(out) :: z(:)     
         integer, intent(in) :: n       ! length of vectors
         real, intent(in) :: f1(:), f2(:), f3(:)
         real, target :: tmp1_ary(n), tmp2_ary(n)
         real, pointer :: tmp1(:), tmp2(:)
         tmp1 => tmp1_ary
         tmp2 => tmp2_ary
         tmp1(1:n) = f2(1:n) - f1(1:n)
         tmp2(1:n) = f3(1:n) - f1(1:n)
         call minmod_sg(z(1:n), n, tmp1(1:n), tmp2(1:n))
         z(1:n) = z(1:n) + f1(1:n)
      end subroutine median_sg


      end module interp_1d_misc_sg
