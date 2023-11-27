! ***********************************************************************
!
!   Copyright (C) 2010-2023  The MESA Team
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

module fingering_modes

   use const_def
   use num_lib
   use math_lib
   use utils_lib

   implicit none

   private

   public :: gaml2max
  
contains

   subroutine gaml2max(pr, tau, r0, lam, beta, ierr, method)

      real(dp), intent(in)               :: pr
      real(dp), intent(in)               :: tau
      real(dp), intent(in)               :: r0
      real(dp), intent(out)              :: lam
      real(dp), intent(out)              :: beta
      integer, intent(out)               :: ierr
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

   subroutine gaml2max_opt_(Pr, tau, R0, lam, beta, ierr)

      real(dp), intent(in)  :: Pr
      real(dp), intent(in)  :: tau
      real(dp), intent(in)  :: R0
      real(dp), intent(out) :: lam
      real(dp), intent(out) :: beta
      integer, intent(out)  :: ierr

      integer, parameter  :: MAX_TRIES = 25
      real(dp), parameter :: EPS = 10*sqrt(EPSILON(0._dp))

      real(dp) :: tlam_max
      real(dp) :: tlam
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

         real(dp), intent(in) :: tlam
         real(dp)             :: lam2_

         ! Evaluate lam2 = lam^2 given tlam

         lam2_ = beta2_(Pr, tau, R0, tlam) * tlam**2

         ! Finish

         return

      end function lam2_

   end subroutine gaml2max_opt_

   !****

   subroutine gaml2max_cubic_(Pr, tau, R0, lam, beta, ierr)

      real(dp), intent(in)  :: Pr
      real(dp), intent(in)  :: tau
      real(dp), intent(in)  :: R0
      real(dp), intent(out) :: lam
      real(dp), intent(out) :: beta
      integer, intent(out)  :: ierr

      real(dp) :: a0, a1, a2, a3
      real(dp) :: q
      real(dp) :: r
      real(dp) :: snq
      real(dp) :: tlam
      
      ! This version directly solves the cubic for the reduced growth rate tlam
      ! (we know that the cubic has three real solutions, with one of them positive)

      ! Set up cubic coefficients

      a3 = Pr - R0 - Pr*R0 + tau
      a2 = -2*(R0-1)*(Pr + tau + pr*tau)
      a1 = pr + tau - 4*Pr*(R0-1)*tau - (1+Pr)*R0*tau**2
      a0 = -2*Pr*tau*(r0*tau - 1)

      ! Determine q and r

      q = a1/(3*a3) - (a2/(3*a3))**2
      r = a1*a2/(6*a3**2) - a0/(2*a3) - (a2/(3*a3))**3

      ! Sanity check (ensures that the cubic has three real roots)

      if (r**2 + q**3 > 0._dp) then
         write(*, *) 'invalid cubic in gaml2max_cubic_'
         ierr = -1
         return
      end if

      ! Calculate the root

      if (q < 0._dp) then
         snq = sqrt(-q)
         tlam = 2*snq*cos(acos(r/snq**3)/3) - a2/3
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
