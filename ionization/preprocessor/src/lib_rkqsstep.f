!
! rkqsstep
! francois hebert, 09/09/01
! 
! code taken from numerical recipes -- integrates by one step
!
! i have added a new argument, args, meant to be a struct for passing
! data to the derivative routines
!


      module lib_rkqsstep

      use def_args
      use def_type, only: dp
      use lib_alert

      implicit none

      real (dp), private, parameter :: PSHRNK = -0.25_dp
      real (dp), private, parameter :: PGROW  = -0.2_dp
      real (dp), private, parameter :: SAFETY = 0.9_dp
      real (dp), private, parameter :: ERRCON = 1.89e-4_dp

      ! the following are constants for use in rkck. will it help the code
      ! to declare them out here, so that they don't need to be reassigned
      ! at each call of rkck?

      real (dp), private, parameter :: a2  = 0.2_dp
      real (dp), private, parameter :: a3  = 0.3_dp
      real (dp), private, parameter :: a4  = 0.6_dp
      real (dp), private, parameter :: a5  = 1.0_dp
      real (dp), private, parameter :: a6  = 0.875_dp
      real (dp), private, parameter :: b21 = 0.2_dp
      real (dp), private, parameter :: b31 = 3.0_dp/40.0_dp
      real (dp), private, parameter :: b32 = 9.0_dp/40.0_dp
      real (dp), private, parameter :: b41 = 0.3_dp
      real (dp), private, parameter :: b42 = -0.9_dp
      real (dp), private, parameter :: b43 = 1.2_dp
      real (dp), private, parameter :: b51 = -11.0_dp/54.0_dp
      real (dp), private, parameter :: b52 = 2.5_dp
      real (dp), private, parameter :: b53 = -70.0_dp/27.0_dp
      real (dp), private, parameter :: b54 = 35.0_dp/27.0_dp
      real (dp), private, parameter :: b61 = 1631.0_dp/55296.0_dp
      real (dp), private, parameter :: b62 = 175.0_dp/512.0_dp
      real (dp), private, parameter :: b63 = 575.0_dp/13824.0_dp
      real (dp), private, parameter :: b64 = 44275.0_dp/110592.0_dp
      real (dp), private, parameter :: b65 = 253.0_dp/4096.0_dp
      real (dp), private, parameter :: c1  = 37.0_dp/378.0_dp
      real (dp), private, parameter :: c3  = 250.0_dp/621.0_dp
      real (dp), private, parameter :: c4  = 125.0_dp/594.0_dp
      real (dp), private, parameter :: c6  = 512.0_dp/1771.0_dp
      real (dp), private, parameter :: dc1 = c1 - 2825.0_dp/27648.0_dp
      real (dp), private, parameter :: dc3 = c3 - 18575.0_dp/48384.0_dp
      real (dp), private, parameter :: dc4 = c4 - 13525.0_dp/55296.0_dp
      real (dp), private, parameter :: dc5 = -277.0_dp/14336.0_dp
      real (dp), private, parameter :: dc6 = c6 - 0.25_dp

      private

      public :: rkqs

      contains



      subroutine rkqs(x, dxin, dxout, dxnext, y, dydx, yscal, eps, derivs, args)

         real (dp), intent(inout) :: x
         real (dp), intent(in) :: dxin, eps
         real (dp), intent(out) :: dxout, dxnext
         real (dp), dimension(:), intent(in) :: dydx, yscal
         real (dp), dimension(:), intent(inout) :: y
         type(datastruct), intent(inout) :: args

         interface
            subroutine derivs(xs, ys, dydxs, args)
               use def_args
               use def_type, only: dp
               implicit none
               real (dp), intent(in) :: xs
               real (dp), dimension(:), intent(in) :: ys
               real (dp), dimension(:), intent(out) :: dydxs
               type(datastruct), intent(inout) :: args
            end subroutine derivs
         end interface

         integer :: int
         real (dp) :: errmax, dx, dxtemp, xnew
         real (dp), dimension(size(y)) :: yerr, ytemp

         if ( size(y) /= size(dydx) .or. size(y) /= size(yscal) ) &
            call alert(1, '(rkqs) different sized vectors in input')

         dx = dxin

         do
            call rkck(x, dx, y, dydx, ytemp, yerr, derivs, args)  ! take a step.
            errmax = maxval( abs(yerr(:)/yscal(:)) )/eps          ! evaluate accuracy.
            if (errmax <= 1.0_dp) exit                            ! step succeeded.
            dxtemp = SAFETY * dx * (errmax**PSHRNK)               ! truncation error too large, reduce stepsize.
            dx = sign( max( abs(dxtemp), 0.1_dp * abs(dx) ), dx)  ! no more than a factor of 10.
            xnew = x + dx
            if (xnew == x) call alert(1, '(rkqs) step size underflow')
         end do

         if (errmax > ERRCON) then                 ! compute size of next step.
            dxnext = SAFETY * dx * (errmax**PGROW)
         else                                      ! no more than a factor of 5 increase.
            dxnext = 5.0_dp * dx
         end if

         dxout = dx
         x = x + dx
         y(:) = ytemp(:)



         contains


      !  
      ! rkck
      !
      ! takes a single dumb rk step, and returns error estimate
      !
         subroutine rkck(x, dx, y, dydx, yout, yerr, derivs, args)

            real (dp), intent(in) :: x, dx
            real (dp), dimension(:), intent(in) :: y, dydx
            real (dp), dimension(:), intent(out) :: yout, yerr
            type(datastruct), intent(inout) :: args

            interface
               subroutine derivs(xs, ys, dydxs, args)
                  use def_args
                  use def_type, only: dp
                  implicit none
                  real (dp), intent(in) :: xs
                  real (dp), dimension(:), intent(in) :: ys
                  real (dp), dimension(:), intent(out) :: dydxs
                  type(datastruct), intent(inout) :: args
               end subroutine derivs
            end interface

            real (dp), dimension(size(y)) :: ak2, ak3, ak4, ak5, ak6, ytemp

            ! while rkck remains a private subroutine called by rkqs, we know these arrays
            ! will have the same length. however, this check should be uncommented if rkck  
            ! becomes a public subroutine!
            !if ( size(y) /= size(dydx) .or. size(dydx) /= size(yout) &
            !  .or. size(yout) /= size(yerr) ) then
            !  call alert(1, '(rkck) different sized vectors in input')
            !end if  

            ! first step
            ytemp = y + dx*b21*dydx

            ! second step
            call derivs(x+a2*dx,ytemp,ak2,args)
            ytemp = y + dx*(b31*dydx + b32*ak2)

            ! third step
            call derivs(x+a3*dx,ytemp,ak3,args)
            ytemp = y + dx*(b41*dydx + b42*ak2 + b43*ak3)

            ! fourth step
            call derivs(x+a4*dx,ytemp,ak4,args)
            ytemp = y + dx*(b51*dydx + b52*ak2 + b53*ak3 + b54*ak4)

            ! fifth step
            call derivs(x+a5*dx,ytemp,ak5,args)
            ytemp = y + dx*(b61*dydx + b62*ak2 + b63*ak3 + b64*ak4 + b65*ak5)

            ! sixth step
            call derivs(x+a6*dx,ytemp,ak6,args)

            ! accumulate increments with proper weights
            yout = y + dx*(c1*dydx + c3*ak3 + c4*ak4 + c6*ak6)
            ! estimate error as difference between fourth and fifth order methods
            yerr = dx * (dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)
   
         end subroutine rkck
   
      end subroutine rkqs

      end module lib_rkqsstep

