module cubic_interpolator
  use const_def, only : dp
  implicit none
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
    self% points_added = 0
    self% solved = .false.
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
          self% matrix(self% points_added + 1, column) = logT**m * logRho**n
          column = column + 1
       enddo
    enddo

    self% vector(self% points_added + 1) = val
    self% points_added = self% points_added + 1
  end subroutine add_point

  subroutine solve_matrix(self)
    class(interpolator), intent(inout) :: self
    integer :: info, ipiv(num_points)
    real(dp) :: matrix_copy(num_points, num_points)
    real(dp) :: vector_copy(num_points)

    if (self% points_added /= num_points) then
       print *, "Need 16 points for cubic interpolation"
       stop
    endif

    matrix_copy = self% matrix
    vector_copy = self% vector

    !write(*,*) "Matrix", matrix_copy
    !write(*,*) "Vector", vector_copy
    call dgesv(num_points, 1, matrix_copy, num_points, ipiv, vector_copy, num_points, info)
    if (info /= 0) then
       write(*,*) "Matrix solve failed in interpolator", info
       write(*,*) "Matrix(:,2)", self% matrix(:,2)
       write(*,*) "Matrix(:,5)", self% matrix(:,5)
       stop
    endif

    self% solved = .true.
    self% solution = vector_copy
  end subroutine solve_matrix

  function evaluate(self, logT, logRho)
    class(interpolator), intent(inout) :: self
    real(dp), intent(in) :: logT
    real(dp), intent(in) :: logRho
    real(dp) :: evaluate
    real(dp) :: sum
    integer :: m, n, column

    if (.not. self% solved) then
       call self% solve_matrix()
    endif

    sum = 0
    column = 1

    do m = 0, max_order
       do n = 0, max_order
          sum = sum + self% solution(column) * logT**m * logRho**n
          column = column + 1
       enddo
    enddo

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

    if (.not. self% solved) then
       call self% solve_matrix()
    endif

    if (deriv_logT .and. deriv_logRho) then
       print *, "May choose only one derivative"
       stop
    endif

    if (.not. deriv_logT .and. .not. deriv_logRho) then
       print *, "Must choose one derivative"
       stop
    endif

    result = 0
    column = 1

    do m = 0, max_order
       do n = 0, max_order
          product = 1

          if (deriv_logT) then
             if (m == 0) then
                product = 0
             else
                product = product * m * logT**(m-1) * logRho**n
             endif
          endif

          if (deriv_logRho) then
             if (n == 0) then
                product = 0
             else
                product = product * n * logT**m * logRho**(n-1)
             end if
          end if

          product = product * self% solution(column)
          result = result + product
          column = column + 1
       enddo
    enddo

    evaluate_deriv = result
  end function evaluate_deriv

end module cubic_interpolator
