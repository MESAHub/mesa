! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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

      module mod_simplex
      
      use const_def, only: dp
      use math_lib
      use num_def
      

      contains

      
      subroutine do_simplex( &
            n, x_lower, x_upper, x_first, x_final, f_final, &
            simplex, f, start_from_given_simplex_and_f, &
            fcn, x_atol, x_rtol, iter_max, fcn_calls_max, &
            centroid_weight_power, enforce_bounds, &
            adaptive_random_search, seed, &
            alpha, beta, gamma, delta, &
            lrpar, rpar, lipar, ipar, &
            num_iters, num_fcn_calls, &
            num_fcn_calls_for_ars, num_accepted_for_ars, ierr)
         integer, intent(in) :: n ! number of dimensions
         real(dp), intent(in) :: x_lower(:), x_upper(:), x_first(:) ! (n)
         real(dp), intent(inout) :: x_final(:) ! (n)
         real(dp), intent(inout) :: simplex(:,:) ! (n,n+1)
         real(dp), intent(inout) :: f(:) ! (n+1)
         logical, intent(in) :: start_from_given_simplex_and_f
         interface
#include "num_simplex_fcn.dek"
         end interface
         real(dp), intent(in) :: x_atol, x_rtol, centroid_weight_power
         integer, intent(inout) :: seed
         real(dp), intent(in) :: alpha, beta, gamma, delta
         integer, intent(in) :: iter_max, fcn_calls_max
         logical, intent(in) :: enforce_bounds, adaptive_random_search
         integer, intent(in) :: lrpar, lipar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp), intent(out) :: f_final
         integer, intent(out) :: num_iters, num_fcn_calls, &
            num_fcn_calls_for_ars, num_accepted_for_ars, ierr
         
         real(dp), dimension(n) :: c, x_reflect, x_expand, x_contract, x_ars
         real(dp) :: f_reflect, f_expand, f_contract, f_ars, rand01, dist, min_dist, &
            dx, term1, xmid, weight, sum_weight, term_val_x
         integer :: h, s, l, i, j, k, iter, min_k
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         num_fcn_calls = 0
         num_fcn_calls_for_ars = 0
         num_accepted_for_ars = 0
         num_iters = 0
         
         if (.not. start_from_given_simplex_and_f) then
            call set_initial_simplex(ierr)
            if (ierr /= 0) then
               write(*,*) 'ierr while evaluating initial simplex'
               return
            end if
         end if
                  
         do 
            
            num_iters = num_iters + 1
            
            if (dbg) write(*,*)
            if (dbg) write(*,2) 'iter', num_iters
            
            ! h = index of max f
            ! s = index of 2nd max f
            ! l = index of min f
            h = 1; l = 1
            do j=2,n+1
               if (f(j) > f(h)) h = j
               if (f(j) < f(l)) l = j
            end do
            s = 0
            do j=1,n+1
               if (j == h) cycle
               if (s <= 0) then
                  s = j
               else if (f(j) > f(s)) then
                  s = j
               end if
            end do
            
            if (dbg) write(*,2) 'worst', h, f(h), simplex(1:n,h)
            if (dbg) write(*,2) '2nd worst', s, f(s), simplex(1:n,s)
            if (dbg) write(*,2) 'best', l, f(l), simplex(1:n,l)

            ! check for domain convergence
            term_val_x = 0
            do j=1,n+1
               if (j == l) cycle
               do i=1,n
                  term1 = abs(simplex(i,j)-simplex(i,l)) / &
                     (x_atol + x_rtol*max(abs(simplex(i,j)), abs(simplex(i,l))))
                  if (term1 > term_val_x) term_val_x = term1
               end do
            end do
            if (dbg) write(*,1) 'term_val_x', term_val_x
            if (term_val_x <= 1d0) exit
            
            ! check for failure to converge in allowed iterations or function calls
            if (num_iters > iter_max .or. num_fcn_calls > fcn_calls_max) exit
         
            ! c = centroid excluding worst point
            c(1:n) = 0
            sum_weight = 0d0
            do j=1,n+1
               if (j == h) cycle
               if (centroid_weight_power == 0d0) then
                  weight = 1d0
               else
                  weight = 1/f(j)
                  if (centroid_weight_power /= 1d0) &
                     weight = pow(weight,centroid_weight_power)
               end if
               do i=1,n
                  c(i) = c(i) + simplex(i,j)*weight
               end do
               sum_weight = sum_weight + weight
            end do
            do i=1,n
               c(i) = c(i)/sum_weight
            end do
            if (dbg) write(*,1) 'c', c(1:n)
         
            ! transform the simplex
            
            call reflect(ierr)
            if (ierr /= 0) return
            if (dbg) write(*,1) 'reflect', f_reflect, x_reflect(1:n)
         
            if (f_reflect < f(s)) then ! accept reflect
               if (dbg) write(*,2) 'accept reflect', num_iters
               do i=1,n
                  simplex(i,h) = x_reflect(i)
               end do
               f(h) = f_reflect
               if (f_reflect < f(l)) then ! try to expand
                  call expand(ierr)
                  if (ierr /= 0) return
                  if (dbg) write(*,1) 'expand', f_expand, x_expand(1:n)
                  if (f_expand < f(l)) then ! accept expand
                     ! note: to keep the simplex large in a good direction,
                     ! we take x_expand even if f_expand > f_reflect
                     do i=1,n
                        simplex(i,h) = x_expand(i)
                     end do
                     f(h) = f_expand
                     if (dbg) write(*,2) 'accept expand', num_iters
                  end if
               end if
            else ! try to contract
               call contract(ierr)
               if (ierr /= 0) return
               if (dbg) write(*,1) 'contract', f_contract, x_contract(1:n)
               if (f_contract < min(f(h),f_reflect)) then ! accept contraction
                  if (dbg) write(*,2) 'accept contraction', num_iters
                  do i=1,n
                     simplex(i,h) = x_contract(i)
                  end do
                  f(h) = f_contract
               else if (adaptive_random_search) then
                  call ARS(ierr)
                  if (ierr /= 0) return
               else
                  call shrink
               end if
            end if
         
         end do
         
         f_final = f(l)
         do i=1,n
            x_final(i) = simplex(i,l)
         end do

         if (dbg) write(*,*)
         if (dbg) write(*,1) 'final', f_final, x_final(1:n)
         if (dbg) write(*,*)
         
         
         contains
         
         
         subroutine set_initial_simplex(ierr)
            integer, intent(out) :: ierr
            integer :: j, i, k
            logical :: okay
            include 'formats'
            
            ierr = 0
            
            do i=1,n
               simplex(i,n+1) = x_first(i)
            end do
            f(n+1) = get_val(simplex(:,n+1), simplex_initial, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,2) 'failed to get value for first simplex point'
               return
            end if
            
            do j=1,n
               do i=1,n
                  simplex(i,j) = x_first(i)
               end do
               if (x_upper(j) - x_first(j) < x_first(j) - x_lower(j)) then
                  ! closer to upper, so displace toward lower
                  simplex(j,j) = x_first(j) - 0.25d0*(x_upper(j) - x_lower(j))
               else
                  simplex(j,j) = x_first(j) + 0.25d0*(x_upper(j) - x_lower(j))
               end if
               okay = .false.
               do k=1,20
                  f(j) = get_val(simplex(:,j), simplex_initial, ierr)
                  if (ierr == 0) then
                     okay = .true.
                     exit
                  end if
                  ! move closer to x_first and retry
                  ierr = 0
                  simplex(j,j) = 0.5d0*(simplex(j,j) + x_first(j))
                  if (abs(simplex(j,j) - x_first(j)) < &
                        1d-12*(1d0 + abs(x_first(j)))) exit
               end do
               if (.not. okay) then
                  ierr = -1
                  if (dbg) write(*,2) 'failed to get value for initial simplex point'
                  return
               end if
            end do

         end subroutine set_initial_simplex
         
         
         real(dp) function get_val(x, op_code, ierr) result(f)
            real(dp), intent(in) :: x(:)
            integer, intent(in) :: op_code ! what nelder-mead is doing for this call
            integer, intent(out) :: ierr
            integer :: i
            include 'formats'
            ierr = 0
            if (enforce_bounds) then
               do i=1,n
                  if (x(i) > x_upper(i) .or. x(i) < x_lower(i)) then
                     if (dbg) write(*,2) 'out of bounds', &
                        num_iters, x(i), x_lower(i), x_upper(i)
                     f = 1d99
                     return
                  end if
               end do
            end if
            num_fcn_calls = num_fcn_calls + 1
            f = fcn(n, x, lrpar, rpar, lipar, ipar, op_code, ierr)
         end function get_val
         
         
         subroutine reflect(ierr)
            integer, intent(out) :: ierr
            integer :: i
            include 'formats'
            ierr = 0
            do i=1,n
               x_reflect(i) = c(i) + alpha*(c(i) - simplex(i,h))
            end do
            f_reflect = get_val(x_reflect, simplex_reflect, ierr)
         end subroutine reflect
         
         
         subroutine expand(ierr)
            integer, intent(out) :: ierr
            integer :: i
            include 'formats'
            ierr = 0
            do i=1,n
               x_expand(i) = c(i) + beta*(x_reflect(i) - c(i))
            end do
            f_expand = get_val(x_expand, simplex_expand, ierr)
         end subroutine expand
         
         
         subroutine contract(ierr)
            integer, intent(out) :: ierr
            integer :: i, op_code
            include 'formats'
            ierr = 0
            if (f_reflect < f(h)) then ! outside contraction
               if (dbg) write(*,1) 'outside contraction'
               op_code = simplex_outside
               do i=1,n
                  x_contract(i) = c(i) + gamma*(x_reflect(i) - c(i))
               end do
            else ! inside contraction
               if (dbg) write(*,1) 'inside contraction'
               op_code = simplex_inside
               do i=1,n
                  x_contract(i) = c(i) + gamma*(simplex(i,h) - c(i))
               end do
            end if
            f_contract = get_val(x_contract, op_code, ierr)
         end subroutine contract
         
         
         subroutine ARS(ierr)
            integer, intent(out) :: ierr
            integer :: i, k, k_max
            include 'formats'
            ierr = 0
            k_max = 100
            do k=1,k_max ! keep trying until find a better random point
               if (num_fcn_calls > fcn_calls_max) exit           
               call get_point_for_ars(ierr)
               if (ierr /= 0) return
               if (dbg) write(*,2) 'adaptive_random_search', num_iters, f_ars, x_ars(1:n)
               if (f_ars <= f(h)) then ! accept adaptive random search
                  if (dbg) write(*,2) 'accept adaptive random search', num_iters, f_ars, x_ars(1:n)
                  do i=1,n
                     simplex(i,h) = x_ars(i)
                  end do
                  f(h) = f_ars
                  num_accepted_for_ars = num_accepted_for_ars + 1
                  return
               end if
               if (dbg) write(*,2) 'reject adaptive random search', num_iters, f_ars, x_ars(1:n)
            end do
         end subroutine ARS
         
         
         subroutine get_point_for_ars(ierr)
            use mod_random, only: r8_uniform_01
            integer, intent(out) :: ierr
            integer :: i
            real(dp) :: rand01
            include 'formats'
            ierr = 0
            do i=1,n
               rand01 = r8_uniform_01(seed)
               x_ars(i) = x_lower(i) + rand01*(x_upper(i) - x_lower(i))
            end do
            if (dbg) write(*,1) 'adaptive random search', x_ars(1:n)
            num_fcn_calls_for_ars = num_fcn_calls_for_ars + 1
            f_ars = get_val(x_ars, simplex_random, ierr)
         end subroutine get_point_for_ars
         
         
         subroutine shrink ! shrink the simplex towards the best point
            integer :: j, i
            include 'formats'
            do j=1,n+1
               if (j == l) cycle
               do i=1,n
                  simplex(i,j) = simplex(i,l) + delta*(simplex(i,j) - simplex(i,l))
               end do
               f(j) = get_val(simplex(:,j), simplex_shrink, ierr)
               if (ierr /= 0) return
               if (dbg) write(*,2) 'shrink', j, f(j), simplex(1:n,j)
            end do
         end subroutine shrink

      
      end subroutine do_simplex



      end module mod_simplex



























