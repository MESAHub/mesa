!
! tfdh
! francois hebert, jul 20 2010
!
! the full version of the tfdh integrator
!

      module mod_tfdh

      use def_args
      use def_const
      use def_type, only: dp
      use lib_alert
      use lib_rkqsstep
      use lib_util
      use mod_density

      implicit none

      contains



      subroutine tfdh(args, profile_flag)

         type(datastruct), intent(inout) :: args
         logical, intent(in) :: profile_flag

         integer, parameter :: firstlength = 64
         integer, parameter :: size = 3
         integer, parameter :: funit = 1
         integer, parameter :: max_iter = 10000000
         real (dp), parameter :: eps = 1.0e-13_dp
         real (dp), parameter :: tiny = 1.0e-30_dp

         logical :: rws_flag
         integer :: istep, i, array_length
         real (dp) :: ztr, dv0, rws
         real (dp) :: ne, nb, qi
         real (dp) :: blast, nlast, vlast, xlast
         real (dp) :: x, xinit, xfinal, dx, dxdid, dxnext
         real (dp), dimension(size) :: y, dydx, yscal

         real (dp), dimension(args%species) :: v_rex
         logical, dimension(args%species) :: rex_flag

         ! assign initial values to all these variables

         ztr = args % ztr
         dv0 = args % dv0
         rws = args % rws

         xinit  = eps
         xfinal = 1e4_dp * rbohr
         dx = xinit/100
         x = xinit

         y(1) = qe*ztr + xinit*dv0
         y(2) = dv0
         y(3) = ztr

         rws_flag = .false.
         rex_flag(:) = .false.

         ! for now define the exclusion radius for each species
         v_rex(:) = args%kt / (qe*args%zs(:))

         xlast = 0.0_dp
         vlast = 0.0_dp
         blast = 0.0_dp

         ! if we're doing file output, need arrays for passing variables up
         if (profile_flag) then
            array_length = firstlength
            call reallocate_dp(args%r,  firstlength)
            call reallocate_dp(args%v,  firstlength)
            call reallocate_dp(args%z,  firstlength)
            call reallocate_dp(args%zb, firstlength)
            call reallocate_dp(args%ne, firstlength)
            call reallocate_dp(args%nb, firstlength)
            call reallocate_dp(args%qitot, firstlength)
            call reallocate_dp(args%qi, args%species*firstlength)
         end if




         do, istep=1, max_iter

            call d_tfdh(x, y, dydx, args)
            yscal(:) = abs(y(:)) + abs(dx*dydx(:)) + tiny
            if ((x+dx-xfinal)*(x+dx-xinit) > 0) dx = xfinal - x

            ! write data to output arrays
            if (profile_flag) then
               args % r(istep)  = x
               args % v(istep)  = y(1)/x
               args % z(istep)  = x/qe * (y(1)/x - y(2))
               args % zb(istep) = y(3)
               args % ne(istep) = args % neloc
               args % nb(istep) = args % nbloc
               args % qitot(istep) = args % qitotloc
               do i=1, args%species
                  args % qi(args%species*(istep-1) + i) = args % qiloc(i)
               end do
            end if

            call rkqs(x, dx, dxdid, dxnext, y, dydx, yscal, eps, d_tfdh, args)

            ! look for exclusion radius for each species
            do, i=1, args%species
               if (.not. rex_flag(i)) then
                  if (y(1)/x < v_rex(i)) then
                     args%rex(i) = (x*(vlast/xlast - v_rex(i)) + xlast*(v_rex(i) - y(1)/x))/(vlast/xlast - y(1)/x)
                     rex_flag(i) = .true.
                  end if
               end if
            end do
               
            ! look for local properties at wigner-seitz radius 
            if (.not. rws_flag) then
               if (x > rws) then
                  !args%nienc = (nithisstep*(rws-xlast)+nilaststep*(x-rws))/(x-xlast)
                  args%nerws = (args%neloc*(rws-xlast)+nlast*(x-rws))/(x-xlast)
                  rws_flag = .true.
               end if
            end if

            ! stop if potential negative, if potential starts increasing, or if bound charge underflows
            if (istep > 1 .and. (y(1) < 0 .or. y(2) > 0 .or. y(3) == blast)) then
               args % output_length = istep
               args % znet = y(3)
               args % rfinal = x
               return
            end if

            if ((x-xfinal)*(xfinal-xinit) >= 0) then
               call alert(1, '(tfdh) reached xfinal before meeting end condition')
            end if

            dx = dxnext

            xlast = x
            vlast = y(1)
            nlast = args%neloc
            blast = y(3)
      
            ! check to see that we don't go past bounds of output arrays
            ! if so, we need bigger arrays to hold our data
            if (profile_flag) then
               if (istep == array_length) then
                  array_length = 2*array_length
                  call reallocate_dp(args%r,  array_length)
                  call reallocate_dp(args%v,  array_length)
                  call reallocate_dp(args%z,  array_length)
                  call reallocate_dp(args%zb, array_length)
                  call reallocate_dp(args%ne, array_length)
                  call reallocate_dp(args%nb, array_length)
                  call reallocate_dp(args%qitot, array_length)
                  call reallocate_dp(args%qi, args%species * array_length)
               end if
            end if
      
         end do
   
         if (istep == max_iter) call alert(1, '(tfdh) too many steps!')
   
      end subroutine tfdh



      subroutine d_tfdh(x, y, dydx, args)
   
         real (dp), intent(in) :: x
         real (dp), intent(in), dimension(:) :: y
         real (dp), intent(out), dimension(:) :: dydx
         type (datastruct), intent(inout) :: args

         real (dp) :: xi, ne_loc, nb_loc, qi_sum

         xi = qe * y(1) / (x * args%kt)
         args % xi = merge(xi, 0.0_dp, xi>0.0_dp)

         ! charge (NOT number) densities
         call qi_vector(args)
         qi_sum = sum(args%qiloc(:))
         ne_loc = ne_total(args)
         nb_loc = ne_bound(args)
         
         ! dV/dr
         dydx(1) = y(2)
         dydx(2) = - 4.0*pi*qe * x * (qi_sum - ne_loc)

         ! dZb/dr
         dydx(3) = - 4.0*pi * x*x * nb_loc

         args % neloc = ne_loc
         args % nbloc = nb_loc
         args % qitotloc = qi_sum
      
      end subroutine d_tfdh

      end module mod_tfdh



   
   