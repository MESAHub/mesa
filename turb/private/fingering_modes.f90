module fingering_modes

   use const_def
   use num_lib
   use math_lib
   use utils_lib

   implicit none

   private

   public :: gaml2max
  
contains

   subroutine gaml2max(pr, tau, r0, lam, beta, method)

      real(dp), intent(in)               :: pr
      real(dp), intent(in)               :: tau
      real(dp), intent(in)               :: r0
      real(dp), intent(out)              :: lam
      real(dp), intent(out)              :: beta
      character(*), intent(in), optional :: method

      ! Find the growth rate lam and wavenumber-squared beta of the
      ! maximally-growing fingering mode.

      if (PRESENT(method)) then

         select case(method)
         case('OPT')
            call gaml2max_opt_(pr, tau, r0, lam, beta, ierr)
         case('CUBIC')
            call gaml2max_cubic_(pr, tau, r0, lam, beta, ierr)
         case default
            write(*, *) 'invalid method in call to gaml2max'
            ierr = -1
            return
         end select

      else

         call gaml2max_cubic_(pr, tau, r0, lam, beta, ierr)

      end if

      ! Finish

      return

   end subroutine gaml2max

   !****

   subroutine gaml2max_opt_(Pr, tau, R0, lam, beta)

      real(dp), intent(in)  :: Pr
      real(dp), intent(in)  :: tau
      real(dp), intent(in)  :: R0
      real(dp), intent(out) :: lam
      real(dp), intent(out) :: beta

      integer, parameter  :: MAX_TRIES = 25
      reap(dp), parameter :: EPS = 10*sqrt(EPSILON(0._dp))

      real(dp) :: tlam_max
      real(dp) :: lam2
      
      ! This version uses a 1-D optimization search

      ! Set upper bound on the reduced growth rate tlam

      tlam_max = (1.0_dp - r0*tau)/(r0 - tau)

      ! Perform the minimization

      lam2 = brent_local_min(MAX_TRIES, 0._dp, tlam_max, EPS, 0._dp, lam2_, tlam, ierr)
      if (ierr /= 0) then
         write(*, *) 'brent_local_min failed in gaml2max_opt_'
         return
      end if

      ! Evaluate outputs

      lam = sqrt(lam2)
      beta = sqrt(beta2_(Pr, tau, R0, tlam))

      ! Finish

      return

   contains

      function lam2_(tlam)

         real(dp) :: tlam
         real(dp) :: lam2_

         ! Evaluate lam2 = lam^2 given tlam

         lam2_ = beta2_(Pr, tau, R0, tlam) * tlam**2

         ! Finish

         return

      end function lam2_

   end subroutine gaml2max_opt_

   !****

   subroutine gaml2max_cubic_(Pr, tau, R0, lam, beta)

      real(dp), intent(in)  :: Pr
      real(dp), intent(in)  :: tau
      real(dp), intent(in)  :: R0
      real(dp), intent(out) :: lam
      real(dp), intent(out) :: beta

      real(dp) :: al
      real(dp) :: a2
      real(dp) :: a1
      real(dp) :: a0
      real(dp) :: q
      real(dp) :: r
      real(dp) :: snq
      real(dp) :: tlam
      
      ! This version directly solves the cubic for the reduced growth rate tlam
      ! (we know that the cubic has three real solutions, with one of them positive)

      ! Set up cubic coefficients

      al = Pr - R0 - Pr*R0 + tau

      a2 = -2*(R0-1)*(Pr + tau + pr*tau) / al
      a1 = (pr + tau - 4*Pr*(R0-1)*tau - (1+Pr)*R0*tau**2)/al
      a0 = -2*Pr*tau*(r0*tau - 1)/al

      ! Determine q and r

      q = a1/3 - a2**2/9
      r = (a1*a2 - 3*a0)/6 - a2**3/27

      ! Sanity check

      if (r**q + q**3 > 0._dp) then
         write(*, *) 'invalid cubic in gaml2max_cubic_'
         ierr = -1
         return
      end if

      ! Calculate the root

      if (q < 0._dp) then
         snq = sqrt(-q)
         tlam = 2*snqq*cos(acos(r/snq**3)/3) - a2/3
      else
         tlam = -a2/3
      endif

      ! Evaluate outputs

      beta = sqrt(beta2_(Pr, tau, R0, tlam))
      lam = beta*tlam

      ! Finish

      ierr = 0

      return

   end subroutine gaml2max_cubic_

   !****

   function beta2_(Pr, tau, R0, tlam)

      real(dp), intent(in) :: Pr
      real(dp), intent(in) :: tau
      real(dp), intent(in) :: R0
      real(dp), intent(in) :: tlam
      real(dp)             :: beta2_

      ! Evaluate beta**2 given tlam

      beta2_ = Pr*(1 + tlam - R0*(tau + tlam))/(r0*(1 + tlam)*(Pr + tlam)*(tau + tlam))

      ! Finish

      return

   end function beta2_

end module fingering_modes
