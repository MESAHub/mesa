!
! gfdi
! francois hebert, jul 15 2010
!
! this module implements the analytical approximation to
! the generalized fermi-dirac integrals (gfdi) presented in
! chabrier & potekhin (phys rev e 1998)
!

      module mod_gfdi

      use def_const
      use def_type, only: dp
      use lib_alert

      implicit none

      contains



      real (dp) function gfdi(order, chi, tau)

         real (dp), intent(in) :: order, chi, tau

         integer :: i, k, dk
         real (dp) :: temp
!<@cond DEV

         real (dp), dimension(5), parameter :: &
            x  = (/ 7.265351e-2_dp, 0.2694608_dp, 0.533122_dp, 0.7868801_dp, 0.9569313_dp /),         &
            xi = (/ 0.26356032_dp, 1.4134031_dp, 3.5964258_dp, 7.0858100_dp, 12.640801_dp /),         &
            h  = (/ 3.818735e-2_dp, 0.1256732_dp, 0.1986308_dp, 0.1976334_dp, 0.1065420_dp /),        &
            v  = (/ 0.29505869_dp, 0.32064856_dp, 7.3915570e-2_dp, 3.6087389e-3_dp, 2.3369894e-5_dp /)

         ! these arrays are 5x3 in chabrier & potekhin, but i implement them as 15x1
         ! the elements are in order: c(1,1)...c(5,1), c(1,2)...c(5,2), c(1,3)...c(5,3)
         real (dp), dimension(15), parameter :: &
            c   = (/ 0.37045057_dp, 0.41258437_dp, 9.777982e-2_dp, 5.3734153e-3_dp, 3.8746281e-5_dp,  &
                     0.39603109_dp, 0.69468795_dp, 0.22322760_dp, 1.5262934e-2_dp, 1.3081939e-4_dp,   &
                     0.76934619_dp, 1.7891437_dp, 0.70754974_dp, 5.6755672e-2_dp, 5.5571480e-4_dp /), &
            khi = (/ 0.43139881_dp, 1.7597537_dp, 4.1044654_dp, 7.7467038_dp, 13.457678_dp,           &
                     0.81763176_dp, 2.4723339_dp, 5.1160061_dp, 9.0441465_dp, 15.049882_dp,           &
                     1.2558461_dp, 3.2070406_dp, 6.1239082_dp, 10.316126_dp, 16.597079_dp /)
!<@endcond DEV


         if (tau > 100.0) call alert(1, "(gfdi) 'tau' outside of convergence region")
         if (order /= 1/2.0_dp .and. order /= 3/2.0_dp .and. order /= 5/2.0_dp ) then
            call alert(1, "(gfdi) invalid 'order'")
         end if


         temp = 0.0_dp
         gfdi = 0.0_dp
         k = floor(order)
         dk = 5*k

         if (chi <= 0.6_dp) then
            do, i=1, 5
               temp = sqrt(1 + khi(i+dk)*tau/2) / (exp(-khi(i+dk)) + exp(-chi))
               gfdi = gfdi + c(i+dk) * temp
            end do
         else if (chi < 14.0_dp) then
            do, i=1, 5
               temp = h(i) * x(i)**k * chi**(k+3/2.0_dp) * sqrt(1 + chi*x(i)*tau/2)
               temp = temp / (1 + exp(chi*(x(i)-1)))
               temp = temp + v(i) * (xi(i) + chi)**(k+1/2.0_dp) * sqrt(1 + (xi(i) + chi)*tau/2)
               gfdi = gfdi + temp
            end do
         else
            ! recycle temp as chabrier & potekhin's r variable
            temp = sqrt( chi*(1 + chi*tau/2) )
            gfdi = f_gfdi(k, chi, tau, temp) + pi*pi/6.0_dp * chi**k * (k + 1/2.0_dp + (k+1.0_dp)*chi*tau/2) / temp
         end if




         contains
            real (dp) function f_gfdi(k, chi, tau, r)

               integer, intent(in) :: k
               real (dp), intent(in) :: chi, tau, r

               f_gfdi = 0.0_dp

               if (chi*tau < 1.e-6 .and. chi > 0.) then
                  f_gfdi = chi**(k+3/2.0_dp)/(k+3/2.0_dp)
               else
                  f_gfdi = (chi + 1/tau)*r/2 - (2*tau)**(-3/2.0_dp) * log(1 + tau*chi + sqrt(2*tau)*r)
                  if (k > 0) then
                     f_gfdi = (2/3.0_dp * r**3 - f_gfdi) / tau
                     if (k > 1) then
                        f_gfdi = (2*chi * r**3 - 5*f_gfdi) / (4*tau)
                     end if
                  end if
               end if
            
            end function f_gfdi
   
      end function gfdi

      end module mod_gfdi

