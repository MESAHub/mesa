! ***********************************************************************
!
!   Copyright (C) 2018-2020  The MESA Team
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

module cubic_interpolator

   use const_def, only: dp

   implicit none

   private
   public :: interpolator

   !4 points on either dimension for cubic interpolation
   integer, parameter :: num_points = 16
   integer, parameter :: max_order = 3

   type interpolator
      real(dp) :: matrix(num_points, num_points)
      real(dp) :: vector(num_points)
      real(dp) :: solution(num_points)
      logical :: solved
      integer :: points_added
   contains
      private
      procedure, public :: initialize
      procedure, public :: add_point
      procedure, public :: solve_matrix
      procedure, public :: evaluate
      procedure, public :: evaluate_deriv
   end type interpolator

contains

   subroutine initialize(self)
      class(interpolator), intent(inout) :: self
      self%points_added = 0
      self%solved = .false.
   end subroutine initialize

   subroutine add_point(self, logT, logRho, val)
      class(interpolator), intent(inout) :: self
      real(dp), intent(in) :: logT
      real(dp), intent(in) :: logRho
      real(dp), intent(in) :: val
      integer :: column
      integer :: m, n

      column = 1
      do m = 0, max_order
         do n = 0, max_order
            self%matrix(self%points_added + 1, column) = logT**m*logRho**n
            column = column + 1
         end do
      end do

      self%vector(self%points_added + 1) = val
      self%points_added = self%points_added + 1
   end subroutine add_point

   subroutine solve_matrix(self)
      class(interpolator), intent(inout) :: self
      integer :: info, ipiv(num_points)
      real(dp) :: matrix_copy(num_points, num_points)
      real(dp) :: vector_copy(num_points)

      if (self%points_added /= num_points) then
         print *, "Need 16 points for cubic interpolation"
         stop
      end if

      matrix_copy = self%matrix
      vector_copy = self%vector

      !write(*,*) "Matrix", matrix_copy
      !write(*,*) "Vector", vector_copy
      call dgesv(num_points, 1, matrix_copy, num_points, ipiv, vector_copy, num_points, info)
      if (info /= 0) then
         write (*, *) "Matrix solve failed in interpolator", info
         write (*, *) "Matrix(:,2)", self%matrix(:, 2)
         write (*, *) "Matrix(:,5)", self%matrix(:, 5)
         stop
      end if

      self%solved = .true.
      self%solution = vector_copy
   end subroutine solve_matrix

   function evaluate(self, logT, logRho)
      class(interpolator), intent(inout) :: self
      real(dp), intent(in) :: logT
      real(dp), intent(in) :: logRho
      real(dp) :: evaluate
      real(dp) :: sum
      integer :: m, n, column

      if (.not. self%solved) then
         call self%solve_matrix()
      end if

      sum = 0
      column = 1

      do m = 0, max_order
         do n = 0, max_order
            sum = sum + self%solution(column)*logT**m*logRho**n
            column = column + 1
         end do
      end do

      evaluate = sum
   end function evaluate

   function evaluate_deriv(self, logT, logRho, deriv_logT, deriv_logRho)
      class(interpolator), intent(inout) :: self
      real(dp), intent(in) :: logT
      real(dp), intent(in) :: logRho
      real(dp) :: evaluate_deriv
      real(dp) :: result, product
      integer :: m, n, column
      logical :: deriv_logT, deriv_logRho

      if (.not. self%solved) then
         call self%solve_matrix()
      end if

      if (deriv_logT .and. deriv_logRho) then
         print *, "May choose only one derivative"
         stop
      end if

      if (.not. deriv_logT .and. .not. deriv_logRho) then
         print *, "Must choose one derivative"
         stop
      end if

      result = 0
      column = 1

      do m = 0, max_order
         do n = 0, max_order
            product = 1

            if (deriv_logT) then
               if (m == 0) then
                  product = 0
               else
                  product = product*m*logT**(m - 1)*logRho**n
               end if
            end if

            if (deriv_logRho) then
               if (n == 0) then
                  product = 0
               else
                  product = product*n*logT**m*logRho**(n - 1)
               end if
            end if

            product = product*self%solution(column)
            result = result + product
            column = column + 1
         end do
      end do

      evaluate_deriv = result

   end function evaluate_deriv

end module cubic_interpolator
