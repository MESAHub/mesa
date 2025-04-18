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

      module interp_1d_misc

      use const_def, only: dp
      use auto_diff

      implicit none

      contains

      subroutine do_integrate_values(init_x, nx, f1, nv, x, vals, ierr)
         real(dp), intent(in) :: init_x(:)  ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx  ! length of init_x vector
         real(dp), intent(in), pointer :: f1(:)  ! =(4, nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv  ! length of new x vector and vals vector
         real(dp), intent(in) :: x(:)  ! (nv)
            ! strictly monotonic in same way as init_x
         real(dp), intent(inout) :: vals(:)  ! (nv)
            ! for i > 1, vals(i) = integral of interpolating poly from x(i-1) to x(i)
            ! vals(1) = 0
         integer, intent(out) :: ierr  ! 0 means aok

         integer :: k_old, k_new
         real(dp) :: xk_old, xkp1_old, xk_new, xk_prev, sum
         logical :: increasing
         real(dp), pointer :: f(:,:)  ! (4, nx)  ! data & interpolation coefficients
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
            real(dp), intent(in) :: x2

            real(dp) :: x0, x1, a1, a2, d1, d2, area

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


      end subroutine do_integrate_values


      subroutine do_interp_values(init_x, nx, f1, nv, x, vals, ierr)
         real(dp), intent(in) :: init_x(:)  ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx  ! length of init_x vector
         real(dp), intent(in), pointer :: f1(:)  ! =(4, nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv  ! length of new x vector and vals vector
         real(dp), intent(in) :: x(:)  ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         real(dp), intent(inout) :: vals(:)  ! (nv)
         integer, intent(out) :: ierr  ! 0 means aok

         integer :: k_old, k_new
         real(dp) :: xk_old, xkp1_old, xk_new, delta
         logical :: increasing
         real(dp), pointer :: f(:,:)  ! (4, nx)  ! data & interpolation coefficients

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
            else  ! decreasing
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


      end subroutine do_interp_values

      subroutine do_interp_values_autodiff(init_x, nx, f1, nv, x, vals, ierr)
         type(auto_diff_real_2var_order1), intent(in) :: init_x(:)  ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx  ! length of init_x vector
         type(auto_diff_real_2var_order1), intent(in), pointer :: f1(:)  ! =(4, nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv  ! length of new x vector and vals vector
         type(auto_diff_real_2var_order1), intent(in) :: x(:)  ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         type(auto_diff_real_2var_order1), intent(inout) :: vals(:)  ! (nv)
         integer, intent(out) :: ierr  ! 0 means aok

         integer :: k_old, k_new
         type(auto_diff_real_2var_order1) :: xk_old, xkp1_old, xk_new, delta
         logical :: increasing
         type(auto_diff_real_2var_order1), pointer :: f(:,:)  ! (4, nx)  ! data & interpolation coefficients

         ierr = 0

         if (nx == 1) then
            vals(1:nv) = f1(1)
            return
         end if

         f(1:4,1:nx) => f1(1:4*nx)

         if(init_x(1) < init_x(nx)) then
            increasing = .true.
         end if

         k_old = 1; xk_old = init_x(k_old); xkp1_old = init_x(k_old+1)

         do k_new = 1, nv

            xk_new = x(k_new)
            if (increasing) then
               if (xk_new > init_x(nx)) then
                  xk_new = init_x(nx)
               else if (xk_new < init_x(1)) then
                  xk_new = init_x(1)
               end if
            else  ! decreasing
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


      end subroutine do_interp_values_autodiff

      subroutine do_interp_values_and_slopes(init_x, nx, f1, nv, x, vals, slopes, ierr)
         real(dp), intent(in) :: init_x(:)  ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx  ! length of init_x vector
         real(dp), intent(in), pointer :: f1(:)  ! =(4, nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv  ! length of new x vector and vals vector
         real(dp), intent(in) :: x(:)  ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         real(dp), intent(inout) :: vals(:)  ! (nv)
         real(dp), intent(inout) :: slopes(:)  ! (nv)
         integer, intent(out) :: ierr  ! 0 means aok

         integer :: k_old, k_new
         real(dp) :: xk_old, xkp1_old, xk_new, delta
         logical :: increasing
         real(dp), pointer :: f(:,:)  ! (4, nx)  ! data & interpolation coefficients
         f(1:4,1:nx) => f1(1:4*nx)

         ierr = 0

         if (nx == 1) then
            vals(1:nv) = f(1,1)
            slopes(1:nv) = 0
            return
         end if

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
            else  ! decreasing
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
                  f(2, k_old) + 2*delta*(f(3, k_old) + 1.5d0*delta*f(4, k_old))

         end do


      end subroutine do_interp_values_and_slopes


      subroutine do_interp2_values_and_slopes( &
            init_x, nx, f1_1, f1_2, nv, x, vals_1, slopes_1, vals_2, slopes_2, ierr)
         real(dp), intent(in) :: init_x(:)  ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx  ! length of init_x vector
         real(dp), intent(in), pointer :: f1_1(:), f1_2(:)  ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv  ! length of new x vector and vals vector
         real(dp), intent(in) :: x(:)  ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         real(dp), intent(inout) :: vals_1(:), vals_2(:)  ! (nv)
         real(dp), intent(inout) :: slopes_1(:), slopes_2(:)  ! (nv)
         integer, intent(out) :: ierr  ! 0 means aok

         integer :: k_old, k_new
         real(dp) :: xk_old, xkp1_old, xk_new, delta
         logical :: increasing
         real(dp), pointer :: f_1(:,:), f_2(:,:)  ! (4, nx)  ! data & interpolation coefficients
         f_1(1:4,1:nx) => f1_1(1:4*nx)
         f_2(1:4,1:nx) => f1_2(1:4*nx)

         ierr = 0

         if (nx == 1) then
            vals_1(1:nv) = f_1(1,1)
            slopes_1(1:nv) = 0
            vals_2(1:nv) = f_2(1,1)
            slopes_2(1:nv) = 0
            return
         end if

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
            else  ! decreasing
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

            vals_1(k_new) =  &
                  f_1(1, k_old) + delta*(f_1(2, k_old)  &
                     + delta*(f_1(3, k_old) + delta*f_1(4, k_old)))

            slopes_1(k_new) =  &
                  f_1(2, k_old) + 2*delta*(f_1(3, k_old) + 1.5d0*delta*f_1(4, k_old))

            vals_2(k_new) =  &
                  f_2(1, k_old) + delta*(f_2(2, k_old)  &
                     + delta*(f_2(3, k_old) + delta*f_2(4, k_old)))

            slopes_2(k_new) =  &
                  f_2(2, k_old) + 2*delta*(f_2(3, k_old) + 1.5d0*delta*f_2(4, k_old))

         end do


      end subroutine do_interp2_values_and_slopes


      subroutine do_interp3_values_and_slopes( &
            init_x, nx, f1_1, f1_2, f1_3, nv, x, &
            vals_1, slopes_1, vals_2, slopes_2, vals_3, slopes_3, ierr)
         real(dp), intent(in) :: init_x(:)  ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx  ! length of init_x vector
         real(dp), intent(in), pointer :: f1_1(:), f1_2(:), f1_3(:)  ! =(4,nx)  ! data & interpolation coefficients
         integer, intent(in) :: nv  ! length of new x vector and vals vector
         real(dp), intent(in) :: x(:)  ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         real(dp), intent(inout) :: vals_1(:), vals_2(:), vals_3(:)  ! (nv)
         real(dp), intent(inout) :: slopes_1(:), slopes_2(:), slopes_3(:)  ! (nv)
         integer, intent(out) :: ierr  ! 0 means aok

         integer :: k_old, k_new
         real(dp) :: xk_old, xkp1_old, xk_new, delta
         logical :: increasing
         real(dp), pointer :: f_1(:,:), f_2(:,:), f_3(:,:)  ! (4, nx)  ! data & interpolation coefficients
         f_1(1:4,1:nx) => f1_1(1:4*nx)
         f_2(1:4,1:nx) => f1_2(1:4*nx)
         f_3(1:4,1:nx) => f1_3(1:4*nx)

         ierr = 0

         if (nx == 1) then
            vals_1(1:nv) = f_1(1,1)
            slopes_1(1:nv) = 0
            vals_2(1:nv) = f_2(1,1)
            slopes_2(1:nv) = 0
            vals_3(1:nv) = f_3(1,1)
            slopes_3(1:nv) = 0
            return
         end if

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
            else  ! decreasing
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

            vals_1(k_new) =  &
                  f_1(1, k_old) + delta*(f_1(2, k_old)  &
                     + delta*(f_1(3, k_old) + delta*f_1(4, k_old)))
            slopes_1(k_new) =  &
                  f_1(2, k_old) + 2*delta*(f_1(3, k_old) + 1.5d0*delta*f_1(4, k_old))

            vals_2(k_new) =  &
                  f_2(1, k_old) + delta*(f_2(2, k_old)  &
                     + delta*(f_2(3, k_old) + delta*f_2(4, k_old)))
            slopes_2(k_new) =  &
                  f_2(2, k_old) + 2*delta*(f_2(3, k_old) + 1.5d0*delta*f_2(4, k_old))

            vals_3(k_new) =  &
                  f_3(1, k_old) + delta*(f_3(2, k_old)  &
                     + delta*(f_3(3, k_old) + delta*f_3(4, k_old)))
            slopes_3(k_new) =  &
                  f_3(2, k_old) + 2*delta*(f_3(3, k_old) + 1.5d0*delta*f_3(4, k_old))

         end do


      end subroutine do_interp3_values_and_slopes


      subroutine do_interp6_values_and_slopes( &
            init_x, nx, f1_1, f1_2, f1_3, f1_4, f1_5, f1_6, nv, x, &
            vals_1, slopes_1, vals_2, slopes_2, vals_3, slopes_3, &
            vals_4, slopes_4, vals_5, slopes_5, vals_6, slopes_6, &
            ierr)
         real(dp), intent(in) :: init_x(:)  ! (nx) ! junction points, strictly monotonic
         integer, intent(in) :: nx  ! length of init_x vector
         real(dp), intent(in), pointer, dimension(:) :: f1_1, f1_2, f1_3, f1_4, f1_5, f1_6
         integer, intent(in) :: nv  ! length of new x vector and vals vector
         real(dp), intent(in) :: x(:)  ! (nv)  ! locations where want interpolated values
            ! strictly monotonic in same way as init_x
            ! values out of range of init_x's are clipped to boundaries of init_x's
         real(dp), intent(inout), dimension(:) :: &
            vals_1, slopes_1, vals_2, slopes_2, vals_3, slopes_3, &
            vals_4, slopes_4, vals_5, slopes_5, vals_6, slopes_6
         integer, intent(out) :: ierr  ! 0 means aok

         integer :: k_old, k_new
         real(dp) :: xk_old, xkp1_old, xk_new, delta
         logical :: increasing
         real(dp), pointer, dimension(:,:) :: f_1, f_2, f_3, f_4, f_5, f_6
         f_1(1:4,1:nx) => f1_1(1:4*nx)
         f_2(1:4,1:nx) => f1_2(1:4*nx)
         f_3(1:4,1:nx) => f1_3(1:4*nx)
         f_4(1:4,1:nx) => f1_4(1:4*nx)
         f_5(1:4,1:nx) => f1_5(1:4*nx)
         f_6(1:4,1:nx) => f1_6(1:4*nx)

         ierr = 0

         if (nx == 1) then
            vals_1(1:nv) = f_1(1,1); slopes_1(1:nv) = 0
            vals_2(1:nv) = f_2(1,1); slopes_2(1:nv) = 0
            vals_3(1:nv) = f_3(1,1); slopes_3(1:nv) = 0
            vals_4(1:nv) = f_4(1,1); slopes_4(1:nv) = 0
            vals_5(1:nv) = f_5(1,1); slopes_5(1:nv) = 0
            vals_6(1:nv) = f_6(1,1); slopes_6(1:nv) = 0
            return
         end if

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
            else  ! decreasing
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

            vals_1(k_new) =  &
                  f_1(1, k_old) + delta*(f_1(2, k_old)  &
                     + delta*(f_1(3, k_old) + delta*f_1(4, k_old)))
            slopes_1(k_new) =  &
                  f_1(2, k_old) + 2*delta*(f_1(3, k_old) + 1.5d0*delta*f_1(4, k_old))

            vals_2(k_new) =  &
                  f_2(1, k_old) + delta*(f_2(2, k_old)  &
                     + delta*(f_2(3, k_old) + delta*f_2(4, k_old)))
            slopes_2(k_new) =  &
                  f_2(2, k_old) + 2*delta*(f_2(3, k_old) + 1.5d0*delta*f_2(4, k_old))

            vals_3(k_new) =  &
                  f_3(1, k_old) + delta*(f_3(2, k_old)  &
                     + delta*(f_3(3, k_old) + delta*f_3(4, k_old)))
            slopes_3(k_new) =  &
                  f_3(2, k_old) + 2*delta*(f_3(3, k_old) + 1.5d0*delta*f_3(4, k_old))

            vals_4(k_new) =  &
                  f_4(1, k_old) + delta*(f_4(2, k_old)  &
                     + delta*(f_4(3, k_old) + delta*f_4(4, k_old)))
            slopes_4(k_new) =  &
                  f_4(2, k_old) + 2*delta*(f_4(3, k_old) + 1.5d0*delta*f_4(4, k_old))

            vals_5(k_new) =  &
                  f_5(1, k_old) + delta*(f_5(2, k_old)  &
                     + delta*(f_5(3, k_old) + delta*f_5(4, k_old)))
            slopes_5(k_new) =  &
                  f_5(2, k_old) + 2*delta*(f_5(3, k_old) + 1.5d0*delta*f_5(4, k_old))

            vals_6(k_new) =  &
                  f_6(1, k_old) + delta*(f_6(2, k_old)  &
                     + delta*(f_6(3, k_old) + delta*f_6(4, k_old)))
            slopes_6(k_new) =  &
                  f_6(2, k_old) + 2*delta*(f_6(3, k_old) + 1.5d0*delta*f_6(4, k_old))

         end do

      end subroutine do_interp6_values_and_slopes


      real(dp) function minmod1(f1, f2)
         real(dp), intent(in) :: f1, f2
         minmod1 = 0.5d0 * (sign(1d0, f1) + sign(1d0, f2)) * min(abs(f1), abs(f2))
      end function minmod1


      real(dp) function median1(f1, f2, f3)
         real(dp), intent(in) :: f1, f2, f3
         median1 = f1 + minmod1(f2 - f1, f3 - f1)
      end function median1


      subroutine minmod(z, n, f1, f2)
         real(dp), intent(inout) :: z(:)
         integer, intent(in) :: n       ! length of vectors
         real(dp), intent(in) :: f1(:), f2(:)
         z(1:n) = 0.5d0 * (sign(1d0, f1(1:n)) + sign(1d0, f2(1:n))) * min(abs(f1(1:n)), abs(f2(1:n)))
      end subroutine minmod


      subroutine minmod4(z, n, f1, f2, f3, f4)
         real(dp), intent(inout) :: z(:)
         integer, intent(in) :: n       ! length of vectors
         real(dp), intent(in) :: f1(:), f2(:), f3(:), f4(:)
         call minmod(z, n, f1, f2)
         call minmod(z, n, z, f3)
         call minmod(z, n, z, f4)
      end subroutine minmod4


      subroutine median(z, n, f1, f2, f3)
         real(dp), intent(inout) :: z(:)
         integer, intent(in) :: n       ! length of vectors
         real(dp), intent(in) :: f1(:), f2(:), f3(:)
         real(dp), target :: tmp1_ary(n), tmp2_ary(n)
         real(dp), pointer :: tmp1(:), tmp2(:)
         tmp1 => tmp1_ary
         tmp2 => tmp2_ary
         tmp1(1:n) = f2(1:n) - f1(1:n)
         tmp2(1:n) = f3(1:n) - f1(1:n)
         call minmod(z(1:n), n, tmp1(1:n), tmp2(1:n))
         z(1:n) = z(1:n) + f1(1:n)
      end subroutine median


      type(auto_diff_real_2var_order1) function minmod1_autodiff(f1, f2)
         use auto_diff
         type(auto_diff_real_2var_order1), intent(in) :: f1, f2
         minmod1_autodiff = 0.5d0 * (sign(f1) + sign(f2)) * min(abs(f1), abs(f2))
      end function minmod1_autodiff


      type(auto_diff_real_2var_order1) function median1_autodiff(f1, f2, f3)
         use auto_diff
         type(auto_diff_real_2var_order1), intent(in) :: f1, f2, f3
         median1_autodiff = f1 + minmod1_autodiff(f2 - f1, f3 - f1)
      end function median1_autodiff


      subroutine minmod_autodiff(z, n, f1, f2)
         use auto_diff
         type(auto_diff_real_2var_order1), intent(inout) :: z(:)
         integer, intent(in) :: n       ! length of vectors
         integer :: k
         type(auto_diff_real_2var_order1), intent(in) :: f1(:), f2(:)
         do k =1, n
         z(k) = 0.5d0 * (sign(f1(k)) + sign(f2(k))) * min(abs(f1(k)), abs(f2(k)))
         end do
      end subroutine minmod_autodiff

      subroutine minmod4_autodiff(z, n, f1, f2, f3, f4)
         use auto_diff
         type(auto_diff_real_2var_order1), intent(inout) :: z(:)
         integer, intent(in) :: n       ! length of vectors
         type(auto_diff_real_2var_order1), intent(in) :: f1(:), f2(:), f3(:), f4(:)
         call minmod_autodiff(z, n, f1, f2)
         call minmod_autodiff(z, n, z, f3)
         call minmod_autodiff(z, n, z, f4)
      end subroutine minmod4_autodiff


      subroutine median_autodiff(z, n, f1, f2, f3)
         use auto_diff
         type(auto_diff_real_2var_order1), intent(inout) :: z(:)
         integer, intent(in) :: n       ! length of vectors
         integer :: k
         type(auto_diff_real_2var_order1), intent(in) :: f1(:), f2(:), f3(:)
         type(auto_diff_real_2var_order1), target :: tmp1_ary(n), tmp2_ary(n)
         type(auto_diff_real_2var_order1), pointer :: tmp1(:), tmp2(:)
         tmp1 => tmp1_ary
         tmp2 => tmp2_ary
         do k =1,n
         tmp1(k) = f2(k) - f1(k)
         tmp2(k) = f3(k) - f1(k)
         end do
         call minmod_autodiff(z, n, tmp1, tmp2)
         do k =1,n
         z(k) = z(k) + f1(k)
         end do
      end subroutine median_autodiff

      end module interp_1d_misc
