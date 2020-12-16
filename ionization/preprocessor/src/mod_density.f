!
! mod_density
! francois hebert, jul 17 2010
!
! module containing functions for calculating electron+ion densities
!
! NOTE!!
! n's refer to number densities
! q's refer to charge denisties
!
! for electrons these are the same in magnitude,
! so use 'n' and avoid potential ambiguity over electron charge sign
!

      module mod_density

      use lib_alert
      use def_args
      use def_type, only: dp
      use mod_gfdi
      use mod_bi

      implicit none

      real (dp), private, parameter :: C_NE  = 3.359080052588e+39_dp !C_NE  = (2 me)^(3/2) / (2 pi^2 hbar^3)
      real (dp), private, parameter :: C_MCC = 8.187104786845e-7_dp  !C_MCC = me c^2
      real (dp), private, parameter :: XI_FACTOR = 1e4_dp


      contains

      ! total ion charge density
      real (dp) function qi_total(args)

         type (datastruct), intent(inout) :: args

         qi_total = sum( args%zs(:)*args%niinf(:) * exp( - args%zs(:)*args%xi ) )
      
      end function qi_total


      ! list of charge densities for the different ion species
      subroutine qi_vector(args)

         type (datastruct), intent(inout) :: args

         if (size(args%qiloc) /= args%species) call alert(1, '(qi_vector) size(args%qiloc) .ne. args%species')

         args%qiloc(:) = args%zs(:)*args%niinf(:) * exp( - args%zs(:)*args%xi )

      end subroutine qi_vector


      ! total electron density
      real (dp) function ne_total(args)

         type (datastruct), intent(in) :: args

         real (dp) :: chi, tau, xi

         chi = args % chi
         tau = args % tau
         xi  = args % xi

         if (xi < 0) then
            ne_total = 0.0_dp
         else
            ne_total = gfdi(1/2.0_dp, chi+xi, tau) + tau * gfdi(3/2.0_dp, chi+xi, tau)
            ne_total = ne_total * C_NE * (tau*C_MCC)**(3/2.0_dp)
         end if
   
      end function ne_total


      ! total plasma (= at inifinity, no potential) electron density
      real (dp) function ne_plasma(chi, tau)
      
         real (dp), intent(in) :: chi, tau
      
         ne_plasma = gfdi(1/2.0_dp, chi, tau) + tau * gfdi(3/2.0_dp, chi, tau)
         ne_plasma = ne_plasma * C_NE * (tau*C_MCC)**(3/2.0_dp)
      
      end function ne_plasma


      ! total density of electrons with E < 0, ie bound
      real (dp) function ne_bound(args)

         type (datastruct), intent(inout) :: args

         real (dp) :: chi, tau, xi

         chi = args % chi
         tau = args % tau
         xi  = args % xi

         if (xi < 0) then
            ne_bound = 0.0_dp
         else if (xi > (abs(chi)+1.0_dp)*XI_FACTOR) then
            ne_bound = gfdi(1/2.0_dp, xi, tau) + tau * gfdi(3/2.0_dp, xi, tau)
            ne_bound = ne_bound * C_NE * (tau*C_MCC)**(3/2.0_dp)
         else
            ne_bound = bound_integral(args)
            ne_bound = ne_bound * C_NE * (tau*C_MCC)**(3/2.0_dp)
         end if
   
      end function ne_bound


      end module mod_density
   
   
   