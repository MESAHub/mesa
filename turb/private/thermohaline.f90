! ***********************************************************************
!
!   Copyright (C) 2010-2021  The MESA Team
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


module thermohaline
   
   use const_def
   use num_lib
   use utils_lib
   use auto_diff
   
   implicit none
   
   private
   public :: get_D_thermohaline

contains
   
   !> Computes the diffusivity of thermohaline mixing when the
   !! thermal gradient is stable and the composition gradient is unstable.
   !!
   !! @param thermohaline_option A string specifying which thermohaline prescription to use.
   !! @param grada Adiabatic gradient dlnT/dlnP
   !! @param gradr Radiative temperature gradient dlnT/dlnP, equals the actual gradient because there's no convection
   !! @param T Temperature
   !! @param opacity opacity
   !! @param rho Density
   !! @param Cp Heat capacity at constant pressure
   !! @param gradL_composition_term dlnMu/dlnP where Mu is the mean molecular weight.
   !! @param iso The index of the species that drives thermohaline mixing.
   !! @param XH1 Mass fraction of H1.
   !! @param thermohaline_coeff Free parameter multiplying the thermohaline diffusivity.
   !! @param D_thrm Output, diffusivity.
   !! @param ierr Output, error index.
   subroutine get_D_thermohaline(thermohaline_option, &
      grada, gradr, T, opacity, rho, Cp, gradL_composition_term, &
      iso, XH1, thermohaline_coeff, D_thrm, ierr)
      character(len = *) :: thermohaline_option
      real(dp), intent(in) :: &
         grada, gradr, T, opacity, rho, Cp, gradL_composition_term, XH1, &
         thermohaline_coeff
      integer, intent(in) :: iso
      real(dp), intent(out) :: D_thrm
      integer, intent(out) :: ierr
      real(dp) :: dgrad, K_therm, K_T, K_mu, nu, R0, Pr, tau, r_th
      include 'formats'
      dgrad = max(1d-40, grada - gradr) ! positive since Schwarzschild stable               
      K_therm = 4d0 * crad * clight * pow3(T) / (3d0 * opacity * rho) ! thermal conductivity
      if (thermohaline_option == 'Kippenhahn') then
         ! Kippenhahn, R., Ruschenplatt, G., & Thomas, H.-C. 1980, A&A, 91, 175
         D_thrm = -3d0 * K_therm / (2 * rho * cp) * gradL_composition_term / dgrad
      else if (thermohaline_option == 'Traxler_Garaud_Stellmach_11' .or. &
         thermohaline_option == 'Brown_Garaud_Stellmach_13') then
         call get_diff_coeffs(K_therm, Cp, rho, T, opacity, iso, XH1, K_T, K_mu, nu)
         R0 = (gradr - grada) / gradL_composition_term
         Pr = nu / K_T
         tau = K_mu / K_T
         r_th = (R0 - 1d0) / (1d0 / tau - 1d0)
         if (r_th >= 1d0) then ! stable if R0 >= 1/tau
            D_thrm = 0d0
         else if (Pr < 0d0) then
            ! Bad results from get_diff_coeffs will just result in NaNs from thermohaline options, so skip
            D_thrm = 0d0
         else if (thermohaline_option == 'Traxler_Garaud_Stellmach_11') then
            ! Traxler, Garaud, & Stellmach, ApJ Letters, 728:L29 (2011).
            ! also see Denissenkov. ApJ 723:563–579, 2010.
            D_thrm = 101d0 * sqrt(K_mu * nu) * exp(-3.6d0 * r_th) * pow(1d0 - r_th, 1.1d0) ! eqn 24
         else ! if (s% thermohaline_option == 'Brown_Garaud_Stellmach_13') then
            D_thrm = K_mu * (Numu(R0, r_th, pr, tau) - 1d0)
         endif
      else
         D_thrm = 0
         ierr = -1
         write(*, *) 'unknown for MLT thermohaline_option' // trim(thermohaline_option)
      end if
      D_thrm = thermohaline_coeff * D_thrm
   end subroutine get_D_thermohaline
   
   
   subroutine get_diff_coeffs(K_therm, Cp, rho, T, opacity, iso, XH1, kt, kmu, vis)
      use chem_def, only : chem_isos
      real(dp), intent(in) :: K_therm, Cp, rho, T, opacity, XH1
      integer, intent(in) :: iso
      real(dp), intent(out) :: kt, kmu, vis
      real(dp) :: loglambdah, loglambdacx, loglambdacy, ccx, ccy, qe4
      real(dp) :: Bcoeff, chemA, chemZ, acx, acy, nu_mol, nu_rad
      real(dp), parameter :: sqrt5 = sqrt(5d0)
      kt = K_therm / (Cp * rho)       ! thermal diffusivity (assumes radiatively dominated)
      qe4 = pow4(qe)
      
      ! Log Lambda for pure H (equation 10 from Proffitt Michaud 93)
      loglambdah = -19.26d0 - 0.5d0 * log(rho) + 1.5d0 * log(T) - 0.5d0 * log(1d0 + 0.5d0 * (1 + XH1))
      nu_rad = 4d0 * crad * pow4(T) / (15d0 * clight * opacity * pow2(rho)) ! radiative viscosity
      nu_mol = 0.406d0 * sqrt(amu) * pow(boltzm * T, 2.5d0) / (qe4 * loglambdah * rho)
      ! From Spitzer "Physics of Fully Ionized Gases equation 5-54
      ! Assumes pure H. Still trying to work out what it would be for a mixture. 
      vis = nu_mol + nu_rad   ! total viscosity
      
      ! The following is from Proffitt & Michaud, 1993.
      ! Their constant B (equation 15)
      Bcoeff = (15.d0 / 16.d0) * sqrt(2.d0 * amu / (5 * pi)) * pow(boltzm, 2.5d0) / qe4
      ! Extract what species drives the thermohaline concvection
      chemA = chem_isos%Z_plus_N(iso)
      chemZ = chem_isos%Z(iso)
      
      if(chemZ.gt.2) then
         ! This is if the driving chemical is NOT He.
         ! Log Lambda for H-dominant chem mixture (equation 10)
         loglambdacx = loglambdah - log(chemz)
         ! Log Lambda for He-dominant chem mixture (equation 10)
         loglambdacy = loglambdah - log(2.d0 * chemz)
         ! Calculation of C_ij coeffs (equation 12)
         ccx = log(exp(1.2d0 * loglambdacx) + 1.) / 1.2d0
         ccy = log(exp(1.2d0 * loglambdacy) + 1.) / 1.2d0
         ! Reduced masses (I had to guess, from Bahcall & Loeb 1990), with H and He
         acx = (1.d0 * chemA) / (1.d0 + chemA)
         acy = 4 * chemA / (4.d0 + chemA)
         ! My formula (see notes) based on Proffitt and Michaud 1993
         kmu = 2 * Bcoeff * pow(T, 2.5d0) / (sqrt5 * rho * chemZ * chemZ) / &
            (XH1 * sqrt(acx) * ccx + (1 - XH1) * sqrt(acy) * ccy)
      
      else
         ! Log Lambda for H-He mixture (equation 10)
         loglambdah = -19.26d0 - log(2d0) - 0.5d0 * log(rho) + &
            1.5d0 * log(T) - 0.5d0 * log(1d0 + 0.5d0 * (1 + XH1))
         ! Calculation of C_ij coeffs (equation 12)
         ccy = log(exp(1.2d0 * loglambdah) + 1d0) / 1.2d0
         ! My formula (see notes) based on Proffitt and Michaud 1993
         kmu = (Bcoeff * pow(T, 2.5d0) / (rho * ccy)) * (3 + XH1) / ((1 + XH1) * (3 + 5 * XH1) * (0.7d0 + 0.3d0 * XH1))
      
      endif
      ! write(57,*) kt,kmu,vis,chemZ
   
   end subroutine get_diff_coeffs
   
   
   real(dp) function numu(R0, r_th, prandtl, diffratio)
      !Function calculates Nu_mu from input parameters, following Brown et al. 2013.
      !Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 
      
      real(dp), intent(in) :: R0, r_th, prandtl, diffratio
      real(dp) :: maxl2, maxl, lambdamax
      real(dp) :: myvars(2)
      integer :: ierr, iter, max_iters
      
      ! Initialize guess using estimates from Brown et al. 2013
      call analytical_estimate_th(maxl, lambdamax, r_th, prandtl, diffratio)
      
      myvars(1) = maxl
      myvars(2) = lambdamax
      
      !Call Newton relaxation algorithm
      call NR(myvars, prandtl, diffratio, R0, ierr)
      
      !If the growth rate is negative, then try another set of parameters as first guess.
      !Repeat as many times as necessary until convergence is obtained.
      iter = 1
      max_iters = 200
      do while(iter<=max_iters .and. ((myvars(2)<0).or.(ierr /= 0)))
         !write(*,*) 'Alternative', r_th,prandtl,diffratio,iter
         !Reset guess values
         myvars(1) = maxl
         myvars(2) = lambdamax
         !Call relaxation for slightly different Pr, tau, R0.
         call NR(myvars, prandtl * (1d0 + iter * 1.d-2), diffratio, R0 / (1d0 + iter * 1.d-2), ierr)
         !If it converged this time, call NR for the real parameters.
         if(ierr.eq.0) call NR(myvars, prandtl, diffratio, R0, ierr)
         !write(*,*) prandtl,diffratio,R0,myvars(1),myvars(2),ierr
         !Otherwise, increase counter and try again.
         iter = iter + 1
      enddo
      
      if((myvars(2)<0).or.(ierr /= 0)) then
         write(*, *) "WARNING: thermohaline Newton relaxation failed to converge, falling back to estimate"
         maxl2 = maxl * maxl
      else ! NR succeeded, so use results in myvars
         !Plug solution into "l^2" and lambda.
         maxl2 = myvars(1) * myvars(1)
         lambdamax = myvars(2)
         !write(*,*) prandtl,diffratio,r_th,maxl2,lambdamax
      end if
      
      !Calculate Nu_mu using Formula (33) from Brown et al, with C = 7.
      numu = 1.d0 + 49.d0 * lambdamax * lambdamax / (diffratio * maxl2 * (lambdamax + diffratio * maxl2))
      
      return
   end function numu
   
   
   subroutine thermohaline_rhs(myx, myf, myj, prandtl, diffratio, R0)
      ! This routine is needed for the NR solver.
      ! Inputs the two following equations for lambda and maxl2:
      ! lambda^3 + a_2 lambda^2 + a_1 lambda + a_0 = 0 (eq. 19 of Brown et al.)
      ! b_2 lambda^2 + b_1 lambda + b_0 = 0 (eq. 20 of Brown et al.)
      ! Inputs f, the equations, and j, their jacobian.
      ! Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 
      
      real(dp), intent(in) :: myx(2), prandtl, diffratio, R0
      real(dp), intent(out) :: myf(2), myj(2, 2)
      real(dp) :: a_2, a_1, a_0, b_2, b_1, b_0, myterm, myx1_2, myx1_3, myx1_4
      
      !This inputs the coefficients.
      b_2 = 1d0 + prandtl + diffratio
      myx1_2 = myx(1) * myx(1)
      myx1_3 = myx1_2 * myx(1)
      myx1_4 = myx1_3 * myx(1)
      a_2 = myx1_2 * b_2
      myterm = diffratio * prandtl + prandtl + diffratio
      b_1 = 2 * myx1_2 * myterm
      a_1 = myx1_4 * myterm + prandtl * (1. - (1d0 / R0))
      b_0 = 3.d0 * myx1_4 * diffratio * prandtl + prandtl * (diffratio - (1d0 / R0))
      a_0 = myx1_4 * myx1_2 * diffratio * prandtl + myx1_2 * prandtl * (diffratio - (1d0 / R0))
      
      !         write(*,*) a_2,a_1,a_0,b_2,b_1,b_0
      
      !These are equations 19 and 20
      myf(1) = ((myx(2) + a_2) * myx(2) + a_1) * myx(2) + a_0
      myf(2) = b_2 * myx(2) * myx(2) + b_1 * myx(2) + b_0
      
      !These are their Jacobians for the NR relaxation.
      myj(1, 1) = 2 * myx(1) * b_2 * myx(2) * myx(2) + &
         4 * myx1_3 * myterm * myx(2) + 6 * myx1_4 * myx(1) * diffratio * prandtl   &
         + 2 * myx(1) * prandtl * (diffratio - (1d0 / R0))
      myj(1, 2) = 3 * myx(2) * myx(2) + 2 * a_2 * myx(2) + a_1
      myj(2, 1) = 4 * myx(1) * myterm * myx(2) + 12.d0 * myx1_3 * diffratio * prandtl
      myj(2, 2) = 2 * b_2 * myx(2) + b_1
      
      return
   end subroutine thermohaline_rhs
   
   
   subroutine analytical_estimate_th(maxl, lambdamax, r_th, prandtl, diffratio)
      !Inputs analytical estimates for l and lambda from Brown et al. 2013.
      
      real(dp) :: prandtl, diffratio, maxl, lambdamax, r_th, phi, maxl4, maxl6
      
      phi = diffratio / prandtl
      
      if(r_th .lt. 0.5d0) then
         if(r_th .gt. prandtl) then
            maxl = pow((1.d0 / (1.d0 + phi)) - 2.d0 * dsqrt(r_th * phi) / pow(1d0 + phi, 2.5d0), 0.25d0)
            ! Equation (B14)
            maxl4 = maxl * maxl * maxl * maxl
            maxl6 = maxl4 * maxl * maxl
            lambdamax = 2 * prandtl * phi * maxl6 / (1d0 - (1d0 + phi) * maxl4)    ! Equation (B11)
         else
            maxl = dsqrt(dsqrt(1d0 / (1d0 + phi)) - dsqrt(prandtl) * (1d0 + phi / ((1d0 + phi) * (1d0 + phi))))
            ! Equation (B5)
            lambdamax = dsqrt(prandtl) - prandtl * dsqrt(1d0 + phi)   !Equation (B5)
         endif
      else
         maxl = pow((one_third) * (1d0 - r_th) + (1d0 - r_th) * (1d0 - r_th) * (5d0 - 4d0 * phi) / 27d0, 0.25d0)
         ! Equation (B19) carried to next order (doesn't work well otherwise)
         maxl4 = maxl * maxl * maxl * maxl
         maxl6 = maxl4 * maxl * maxl
         lambdamax = 2d0 * prandtl * phi * maxl6 / (1d0 - (1d0 + phi) * maxl4) ! Equation (B11)
      endif
      if(lambdamax<0) then   ! shouldn't be needed, but just as precaution
         maxl = 0.5d0
         lambdamax = 0.5d0
      endif
      
      return
   end subroutine analytical_estimate_th
   
   
   subroutine NR(xrk, prandtl, diffratio, R0, ierr)
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
      
      integer :: iter, ierr
      real(dp) :: xrk(2), f(2) ! Functions f 
      real(dp) :: j(2, 2) ! Jacobian
      real(dp) :: err, errold ! Error at each iteration
      real(dp) :: x1_sav, x2_sav
      real(dp) :: prandtl, diffratio, R0
      real(dp) :: A(lda, n), AF(ldaf, n), R(n), C(n), B(ldb, nrhs), X(ldx, nrhs), &
         rcond, ferr(nrhs), berr(nrhs), work(4 * n)
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
         call thermohaline_rhs(xrk, f, j, prandtl, diffratio, R0)
         
         fact = 'E'
         trans = 'N'
         equed = ''
         
         A = j
         B(1, 1) = f(1)
         B(2, 1) = f(2)
         
         call dgesvx(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, &
            equed, r, c, B, ldb, x, ldx, rcond, ferr, berr, &
            work, iwork, ierr)
         
         if (ierr /= 0) then
            !write(*,*) 'dgesvx failed in thermohaline routine', iter
            !write(*,2) j(1,1),j(1,2)
            !write(*,2) j(2,1),j(2,2)
         else
            iter = iter + 1
            f(1) = X(1, 1)
            f(2) = X(2, 1)
            err = dsqrt(f(1) * f(1) + f(2) * f(2)) ! Calculate the new error
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

end module thermohaline
