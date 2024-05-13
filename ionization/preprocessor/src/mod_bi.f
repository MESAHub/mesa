!
! mod_bi
! francois hebert, jul 17 2010
!
! integrates the bound-electron contribution to the fermi-dirac integral
!

      module mod_bi

      use def_args
      use def_type, only: dp
      use lib_alert
      use lib_rkqsstep

      implicit none

      contains


      real (dp) function bound_integral(args)

         type(datastruct), intent(inout) :: args

         integer, parameter :: size = 1
         integer, parameter :: max_iter = 100000
         real (dp), parameter :: eps = 1.0e-12_dp
         real (dp), parameter :: tiny = 1.0e-30_dp

         integer :: istep
         real (dp) :: xi, cut, ylast 
         real (dp) :: x, xinit, xfinal, dx, dxdid, dxnext
         real (dp), dimension(size) :: y, dydx, yscal

         bound_integral = 0.0_dp

         ! assign initial values to all these variables

         xi  = args % xi
         cut = args % cut

         xinit  = 0.0_dp
         xfinal = xi + cut
         dx = (xfinal - xinit)/100
         x = xinit

         y(1) = 0.0_dp

         if (xfinal <= xinit) return

         do, istep=1, max_iter
   
            call d_fermidirac(x, y, dydx, args)
            yscal(:) = abs(y(:)) + abs(dx*dydx(:)) + tiny
            if ((x+dx-xfinal)*(x+dx-xinit) > 0) dx = xfinal - x

            call rkqs(x, dx, dxdid, dxnext, y, dydx, yscal, eps, d_fermidirac, args)

            ! done when integral underflows or reaches final bound
            if (y(1) == ylast .and. istep > 0 .or. (x-xfinal)*(xfinal-xinit) >= 0) then
               bound_integral = y(1)
               return
            end if

            dx = dxnext
            ylast = y(1)
         end do

         call alert(1, '(bound_integral) too many steps!')

      end function bound_integral


   
      subroutine d_fermidirac(x, y, dydx, args)

         real (dp), intent(in) :: x
         real (dp), intent(in), dimension(:) :: y
         real (dp), intent(out), dimension(:) :: dydx
         type (datastruct), intent(inout) :: args

         real (dp) :: chi, tau, xi

         chi = args % chi
         tau = args % tau
         xi  = args % xi

         dydx(1) = (1.0 + tau*x) * sqrt(x + tau*x*x/2.0) / (exp(x-chi-xi) + 1.0)

      end subroutine d_fermidirac



      end module mod_bi


   
   
   