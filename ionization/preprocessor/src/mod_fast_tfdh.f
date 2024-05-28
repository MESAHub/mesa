!
! tfdh_fast
! francois hebert, jul 17 2010
!
! a stripped-down tfdh integrator which does not attempt to calculate bound
! charge (the speed-limiting step). it is called by root_vex to find the
! value of potential at the center which matches the b.c. at infinity.
!

      module mod_fast_tfdh

      use def_args
      use def_const
      use def_type, only: dp
      use lib_rkqsstep
      use mod_density

      implicit none

      contains


      logical function fasttfdh(args)

         type(datastruct), intent(inout) :: args

         integer, parameter :: size = 2
         integer, parameter :: max_iter = 10000000
         real (dp), parameter :: eps = 1.0e-13_dp
         real (dp), parameter :: tiny = 1.0e-30_dp

         integer :: istep
         real (dp) :: ztr, dv0
         real (dp) :: xlast, vlast, dlast !,xstop
         real (dp) :: x, xinit, xfinal, dx, dxdid, dxnext
         real (dp), dimension(size) :: y, dydx, yscal

         ! assign initial values to all these variables

         ztr = args % ztr
         dv0 = args % dv0

         xinit  = eps
         xfinal = 1.e4_dp * rbohr
         dx = xinit/100
         x = xinit

         y(1) = qe*ztr + xinit*dv0
         y(2) = dv0

         xlast = 0.0_dp
         vlast = 0.0_dp
         dlast = 0.0_dp

         do, istep=1, max_iter
      
            call d_fasttfdh(x, y, dydx, args)
            yscal(:) = abs(y(:)) + abs(dx*dydx(:)) + tiny
            if ((x+dx-xfinal)*(x+dx-xinit) > 0) dx = xfinal - x

            call rkqs(x, dx, dxdid, dxnext, y, dydx, yscal, eps, d_fasttfdh, args)

            !stop if potential negative, or if potential starts increasing
            if (y(1) < 0 .and. istep > 1) then
               !xstop = (y(1)*xlast - vlast*x)/(y(1)-vlast)
               fasttfdh = .false.
               return
            else if (y(2) > 0 .and. istep > 1) then
               !xstop = (y(2)*xlast - dlast*x)/(y(2)-dlast)
               fasttfdh = .true.
               return
            end if

            if ((x-xfinal)*(xfinal-xinit) >= 0) then
               call alert(0, '(tfdh_fast) reached xfinal before meeting end condition')
               !x_stop = x_final
               fasttfdh = .true.
               return
            end if

            dx = dxnext
            xlast = x
            vlast = y(1)
            dlast = y(2)
         end do

         call alert(1, '(tfdh_fast) too many steps!')
   
      end function fasttfdh



      subroutine d_fasttfdh(x, y, dydx, args)
   
         real (dp), intent(in) :: x
         real (dp), intent(in), dimension(:) :: y
         real (dp), intent(out), dimension(:) :: dydx
         type (datastruct), intent(inout) :: args

         real (dp) :: xi, ne_loc, qi_sum

         xi = qe * y(1) / (x * args%kt)
         args % xi = merge(xi, 0.0_dp, xi>0.0_dp)

         ! charge (NOT number) densities
         qi_sum = qi_total(args)
         ne_loc = ne_total(args)

         ! dV/dr
         dydx(1) = y(2)
         dydx(2) = - 4.0*pi*qe * x * (qi_sum - ne_loc)

      end subroutine d_fasttfdh



      end module mod_fast_tfdh


   
   
   