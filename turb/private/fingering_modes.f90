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

   public :: gaml2max, calc_mode_properties
  
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
         tlam = 2*snq*cos(acos(r/snq**3)/3) - a2/(3*a3)
      else
         tlam = -a2/(3*a3)
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

   
   ! Routines for the 2D root numerical root find method
   subroutine calc_mode_properties(R0,prandtl,tau,maxl2,maxl,lambdamax)
      real(dp), intent(in) :: R0,prandtl,tau
      real(dp), intent(out) :: maxl2,maxl,lambdamax
      real(dp) :: myvars(2), r_th
      integer :: ierr, iter, max_iters

      r_th = (R0 - 1d0)/(1d0/tau - 1d0)

      ! Initialize guess using estimates from Brown et al. 2013
      call analytical_estimate_th(maxl,lambdamax,r_th,prandtl,tau)
            
      myvars(1) = maxl
      myvars(2) = lambdamax

      !Call Newton relaxation algorithm
      call NR(myvars,prandtl,tau,R0,ierr)
      
      !If the growth rate is negative, then try another set of parameters as first guess.  
      !Repeat as many times as necessary until convergence is obtained.
      iter = 1
      max_iters = 200
      do while(iter<=max_iters .and. ((myvars(2)<0).or.(ierr /= 0))) 
         !write(*,*) 'Alternative', r_th,prandtl,tau,iter
         !Reset guess values
         myvars(1) = maxl
         myvars(2) = lambdamax
         !Call relaxation for slightly different Pr, tau, R0.
         call NR(myvars,prandtl*(1d0+iter*1.d-2),tau,R0/(1d0+iter*1.d-2),ierr)
         !If it converged this time, call NR for the real parameters.
         if(ierr.eq.0) call NR(myvars,prandtl,tau,R0,ierr)
         !write(*,*) prandtl,tau,R0,myvars(1),myvars(2),ierr
         !Otherwise, increase counter and try again.
         iter = iter + 1            
      enddo
      
      if((myvars(2)<0).or.(ierr /= 0)) then
         write(*,*) "WARNING: thermohaline Newton relaxation failed to converge, falling back to estimate"
         maxl2 = maxl*maxl
      else ! NR succeeded, so use results in myvars
         !Plug solution into "l^2" and lambda.
         maxl2 = myvars(1)*myvars(1)
         lambdamax = myvars(2)
         !write(*,*) prandtl,tau,r_th,maxl2,lambdamax
      end if
      
   end subroutine calc_mode_properties

   subroutine thermohaline_rhs(myx,myf,myj,prandtl,diffratio,R0)
      ! This routine is needed for the NR solver.
      ! Inputs the two following equations for lambda and maxl2:
      ! lambda^3 + a_2 lambda^2 + a_1 lambda + a_0 = 0 (eq. 19 of Brown et al.)
      ! b_2 lambda^2 + b_1 lambda + b_0 = 0 (eq. 20 of Brown et al.)
      ! Inputs f, the equations, and j, their jacobian.
      ! Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 

      real(dp), intent(in) :: myx(2),  prandtl, diffratio, R0
      real(dp), intent(out) :: myf(2), myj(2,2)
      real(dp) :: a_2,a_1,a_0,b_2,b_1,b_0,myterm,myx1_2,myx1_3,myx1_4

      !This inputs the coefficients.
      b_2 = 1d0+prandtl+diffratio
      myx1_2 = myx(1)*myx(1)
      myx1_3 = myx1_2*myx(1)
      myx1_4 = myx1_3*myx(1)
      a_2 = myx1_2*b_2
      myterm = diffratio*prandtl+prandtl+diffratio
      b_1 = 2*myx1_2*myterm
      a_1 = myx1_4*myterm + prandtl*(1. - (1d0/R0))
      b_0 = 3.d0*myx1_4*diffratio*prandtl + prandtl*(diffratio - (1d0/R0))
      a_0 = myx1_4*myx1_2*diffratio*prandtl + myx1_2*prandtl*(diffratio - (1d0/R0))

   !         write(*,*) a_2,a_1,a_0,b_2,b_1,b_0

      !These are equations 19 and 20
      myf(1) = ((myx(2) + a_2)*myx(2) + a_1)*myx(2) + a_0
      myf(2) = b_2*myx(2)*myx(2) + b_1*myx(2) + b_0

      !These are their Jacobians for the NR relaxation.
      myj(1,1) = 2*myx(1)*b_2*myx(2)*myx(2) + &
         4*myx1_3*myterm*myx(2) + 6*myx1_4*myx(1)*diffratio*prandtl   &
           + 2*myx(1)*prandtl*(diffratio - (1d0/R0))
      myj(1,2) = 3*myx(2)*myx(2) + 2*a_2*myx(2) + a_1
      myj(2,1) = 4*myx(1)*myterm*myx(2) + 12.d0*myx1_3*diffratio*prandtl
      myj(2,2) = 2*b_2*myx(2) + b_1

      return
   end subroutine thermohaline_rhs               


   subroutine analytical_estimate_th(maxl,lambdamax,r_th,prandtl,diffratio)
      !Inputs analytical estimates for l and lambda from Brown et al. 2013.

      real(dp) :: prandtl, diffratio, maxl, lambdamax, r_th, phi, maxl4, maxl6

      phi = diffratio/prandtl

      if(r_th .lt. 0.5d0) then
         if(r_th .gt. prandtl) then
            maxl = pow((1.d0/(1.d0+phi)) - 2.d0*dsqrt(r_th*phi)/pow(1d0+phi,2.5d0),0.25d0)   
               ! Equation (B14)
            maxl4 = maxl*maxl*maxl*maxl
            maxl6 = maxl4*maxl*maxl
            lambdamax = 2*prandtl*phi*maxl6/(1d0-(1d0+phi)*maxl4)    ! Equation (B11)
         else
            maxl = dsqrt(dsqrt(1d0/(1d0+phi)) - dsqrt(prandtl)*(1d0+phi/((1d0+phi)*(1d0+phi))))  
               ! Equation (B5)
            lambdamax = dsqrt(prandtl) - prandtl*dsqrt(1d0+phi)   !Equation (B5)
         endif
      else
         maxl = pow((one_third)*(1d0-r_th) + (1d0-r_th)*(1d0-r_th)*(5d0-4d0*phi)/27d0,0.25d0)
            ! Equation (B19) carried to next order (doesn't work well otherwise)
         maxl4 = maxl*maxl*maxl*maxl
         maxl6 = maxl4*maxl*maxl
         lambdamax = 2d0*prandtl*phi*maxl6/(1d0-(1d0+phi)*maxl4) ! Equation (B11)
      endif
      if(lambdamax<0) then   ! shouldn't be needed, but just as precaution
         maxl = 0.5d0
         lambdamax = 0.5d0
      endif

      return
   end subroutine analytical_estimate_th


   subroutine NR(xrk,prandtl,diffratio,R0,ierr)
      ! Newton Relaxation routine used to solve cubic & quadratic in thermohaline case.
      ! Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 

      real(dp), parameter :: acy = 1.d-13 ! accuracy of NR solution.
      integer, parameter :: niter = 20  ! max number of iterations allowed before giving up.
      integer, parameter :: &  !array dimension input parameters for dgesvx
            n = 2, &
            nrhs = 1, &
            lda = n, &
            ldaf = n, &
            ldb = n, &
            ldx = n

      integer :: iter,ierr
      real(dp) :: xrk(2), f(2) ! Functions f 
      real(dp) :: j(2,2) ! Jacobian
      real(dp) :: err,errold ! Error at each iteration
      real(dp) :: x1_sav,x2_sav
      real(dp) :: prandtl, diffratio, R0
      real(dp) :: A(lda,n), AF(ldaf,n), R(n), C(n), B(ldb,nrhs), X(ldx,nrhs), &
            rcond, ferr(nrhs), berr(nrhs), work(4*n)
      character :: fact, trans, equed
      integer :: ipiv(n), iwork(n)

      include 'formats'

      !Initialize flags and other counters.
      ierr = 0
      iter = 0
      err = 1d99
      errold = 1d99
      !Save input guess
      x1_sav = xrk(1)
      x2_sav = xrk(2)

      !While error is too large .and. decreasing, iterate.
      do while ((err.gt.acy).and.(ierr.eq.0).and.(iter.lt.niter))
         call thermohaline_rhs(xrk,f,j,prandtl,diffratio,R0)    
      
         fact = 'E'
         trans = 'N'
         equed = ''
         
         A  = j
         B(1,1) = f(1)
         B(2,1) = f(2)

         call dgesvx( fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, &
            equed, r, c, B, ldb, x, ldx, rcond, ferr, berr, &
            work, iwork, ierr )

         if (ierr /= 0) then
            !write(*,*) 'dgesvx failed in thermohaline routine', iter
            !write(*,2) j(1,1),j(1,2)
            !write(*,2) j(2,1),j(2,2)
         else
            iter = iter + 1
            f(1) = X(1,1)
            f(2) = X(2,1)
            err = dsqrt(f(1)*f(1)+f(2)*f(2)) ! Calculate the new error
            ! If, after a while, the error is still not decreasing, give up and exit NR.
            ! Otherwise, continue.
            if((iter.gt.5).and.(err.gt.errold)) then              
               ! Write(*,2) 'Error not decreasing at iter', iter, err, errold
               ierr = 1
               ! Reset xs and exit loop.
               xrk(1) = x1_sav
               xrk(2) = x2_sav                   
            else
               xrk = xrk - f ! The solution is now in f, so update x 
               errold = err
            endif
         endif
      enddo
      
      if(err<=acy) then
         ierr = 0
      else
         ! write(*,2) 'NR failed to converge', err, iter
         ierr = 1
      end if

      return
   end subroutine NR
   
end module fingering_modes
