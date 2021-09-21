   ! ***********************************************************************
   !
   !   Copyright (C) 2018  Sam Jones, Robert Farmer & The MESA Team
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
   
   ! Implement screening a la Chugunov, DeWitt & Yakovlev 2007, PhRvD, 76, 025028
   
   module screening_chugunov
      use math_lib
      use const_def
      use rates_def, only: screen_info
   
      implicit none
      private
   
      ! there are various realms of applicability in the Chugunov paper, for
      ! the particle-in-cell Monte-Carlo calculations, for the MKB
      ! approximation, and for the Potekhin-Cabrier-based fit to the MKB
      ! results for the mean field potential H(r).
      ! hence, here are parameters for smoothly limiting the domain of application
      real(dp), parameter ::   gamfitlim = 600.0d0 !< upper limit for applicability of fit for h(gamma)
      real(dp), parameter ::   g0 = 590.0d0        !< lower limit to start blending from
      real(dp), parameter ::   deltagam = gamfitlim - g0 ! How much to blend over the gamma boundaries
      real(dp), parameter ::   tp2 = 0.2d0         !< minimum allowed ratio T/T_p before fading out
      real(dp), parameter ::   tp1 = 0.1d0         !< floor value of T/T_p (T_p = plasma temperature)
      real(dp), parameter ::   deltatp = tp2 - tp1         !< Blend over the tp boundary
   
      ! coefficients from chugunov
      real(dp), parameter ::   c_a1 = 2.7822d0
      real(dp), parameter ::   c_a2 = 98.34d0
      real(dp), parameter ::   c_a3 = sqrt(3.0d0) - c_a1/sqrt(c_a2)
      real(dp), parameter ::   c_b1 = -1.7476d0
      real(dp), parameter ::   c_b2 = 66.07d0
      real(dp), parameter ::   c_b3 = 1.12d0
      real(dp), parameter ::   c_b4 = 65.d0
      real(dp), parameter ::   alfa = 0.022d0
   
      real(dp), parameter :: x13   = 1.0d0/3.0d0 
      real(dp), parameter :: x12   = 1.0d0/2.0d0 
      real(dp), parameter :: x23   = 2.0d0/3.0d0 
      real(dp), parameter :: x32   = 3.0d0/2.0d0   
      real(dp), parameter :: x43   = 4.0d0/3.0d0 
      real(dp), parameter :: x53   = 5.0d0/3.0d0 
   
      real(dp), parameter :: h0fitlim = 300d0, h0fit0 = 295d0
      real(dp), parameter :: deltah0fit = h0fitlim - h0fit0
   
      real(dp), dimension(:), allocatable :: z13 ! Z^1/3
      logical :: have_initialization = .false.
   
   
      public  eval_screen_chugunov, screen_chugunov_init, free_chugunov
   
   contains
   

      subroutine screen_chugunov_init()
         integer :: i
      
         if(have_initialization) return
!$omp critical  (omp_critical_screen_chugunov_init)
         if(.not. have_initialization) then
            allocate(z13(0:150))
            do i=lbound(z13,dim=1), ubound(z13,dim=1)
               z13(i) = pow(i*1d0,x13)
            end do
            have_initialization = .true.
         end if
!$omp end critical  (omp_critical_screen_chugunov_init) 
      
      end subroutine  screen_chugunov_init
   

      subroutine free_chugunov()

!$omp critical  (omp_critical_screen_free_chugunov)
         if(allocated(z13)) deallocate(z13)
         have_initialization = .false.
!$omp end critical  (omp_critical_screen_free_chugunov)

      end subroutine free_chugunov

   
      subroutine eval_screen_chugunov(sc, z1, z2, a1, a2, screen, dscreendt, dscreendd, ierr)
         implicit none
   
         type (Screen_Info)  :: sc
         real(dp),intent(in) ::    z1, z2      !< charge numbers of reactants
         real(dp),intent(in) ::    a1, a2     !< mass numbers of reactants
         real(dp),intent(out) ::   screen     !< on return, screening factor for this reaction
         real(dp),intent(out) ::   dscreendt     !< on return, temperature derivative of the screening factor
         real(dp),intent(out) ::   dscreendd    !< on return, density derivative of the screening factor
         integer, intent(out) ::   ierr
   
         real(dp) ::   z1z2       !< z1*z2
         real(dp) ::   t1, dt1dt, dt1dd 
         real(dp) ::   t2, dt2dt, dt2dd
         real(dp) ::   t3, dt3dt, dt3dd !< the three terms in the fitting formula for h0fit
         real(dp) ::   h0fit, dh0fitdt, dh0fitdd !< mean field fit
         real(dp) ::   denom, ddenomdt, ddenomdd, denom2
         real(dp) ::   tp, dtpdd, wp, dwpdd     !< plasma temperature and frequency
         real(dp) ::   tn, dtndd, dtndt         !< normalised temperature (tk/tp)
         real(dp) ::   beta, dbetadt, dbetadd, gama, dgamadt, dgamadd !< coeffs of zeta formula
         real(dp) ::   gam, dgamdt, dgamdd
         real(dp) ::   gamtild, dgamtilddt, dgamtilddd    !< "effective" coulomb coupling parameter
         real(dp) ::   gamtild2   !< "effective" coulomb coupling parameter squared
         real(dp) ::   zeta, dzetadt, dzetadd       !< function of plasma temperature
         real(dp) ::   zeta2, zeta3
         real(dp) ::   a_1, da_1dd, a_2,da_2dd   !< ion sphere radii
         real(dp) ::   a_av, da_avdd       !< average ion sphere radii
         real(dp) ::   a_e        !< electron sphere radius
         real(dp) ::   a_b        !< bohr radius
         real(dp) ::   rs         !< ion sphere radius normalised to bohr radius
         real(dp) ::   m1, m2     !< ion masses
         real(dp) ::   n1, n2     !< number density of ions of types 1 and 2
         real(dp) ::   ntot, dntotdd       !< total number density of these ions
         real(dp) ::   s          !< evaluated sigmoid function
         real(dp) ::   mav           !< average ion mass
         real(dp) ::   A, B, C, U, dAdt, dAdd, dBdt, dBdd, dCdt, dCdd
         real(dp) ::   temp, rho, abar, zbar, rr
         real(dp) ::   alpha, dalphadgam, dbetadgam, dalphadh0, dbetadh0, dalphadtn,dbetadtn
         real(dp) ::   tk, dtkdt, dtkdd, dtk2_dtp2dd
   
         ! check whether both reactants are charged ions
   
         screen = 1d0
         dscreendt = 0d0
         dscreendd = 0d0
         ierr = 0
   
         ! Must be charged ions
         if (z1 <= 0d0 .or. z2  <= 0d0) then
            return
         end if
         
   
         rho = sc% den
         temp = sc% temp
         zbar = sc% zbar
         abar = sc% abar
         rr = sc% rr ! den/abar
         
         ! ion masses and number densities
         mav   = abar * amu
         ! ntot  = den / (amu*abar)
         ntot  = sc% ntot
         dntotdd = 1d0/(amu*abar)
   
         ! ion and electron sphere radii (itoh 1979 eq 1-3)
         !a_e = pow((3.d0 /(pi4 * zbar * ntot)),x13)
         a_e = sc% a_e

         a_1 = a_e * z13(int(z1)) !pow(z1,x13)
         a_2 = a_e * z13(int(z2)) !pow(z2,x13)
         
         da_1dd = -x13 * a_1/rho
         da_2dd = -x13 * a_2/rho
            
         a_av  = 0.5d0 * (a_1 + a_2)
         da_avdd = 0.5d0 * (da_1dd + da_2dd) 

         z1z2 = z1 * z2
            
         ! bohr radius and normalised ion sphere radius
   
         a_b = rbohr/z1z2/ntot
         rs  = a_av / a_b
   
         ! plasma frequency and temperature (chugunov 2007 eq 2)
   
         wp = sqrt((pi4 * z1z2 * qe * qe * ntot/mav))
         dwpdd = x12 * wp/rho
         
         tp = hbar*wp/kerg
         dtpdd = (hbar/kerg) * dwpdd
         
         tn = temp/tp
         dtndt = 1d0/tp
         dtndd = -tn/tp * dtpdd
         
         ! Revise temperature used if nearing low tp values
         tk = temp
         dtkdt = 1d0
         dtkdd = 0d0
         denom = tp * (tk * dtpdd - tp * dtkdd)/(tk * tk * tk)

         ! Blend out for low Tp values
         if ( tn .le. tp1)then
            tk = tp1 * tp
            dtkdt = 0d0
            dtkdd = tp1 * dtpdd
            denom = 0d0 ! Set as zero overwise floating point issues cause small error in the subtraction (tk * dtpdd - tp * dtkdd)
         
         else if (tn .gt.tp1 .and. tn .le. tp2) then
         
            alpha = 0.5d0 * (1d0 - cospi((tn-tp1)/deltatp))
            dalphadtn = 0.5d0 * (pi/deltatp) * sinpi((tn-tp1)/deltatp)
            
            beta = 1.d0 - alpha
            dbetadtn = -dalphadtn
            
            dtkdt = (dbetadtn * dtndt * tp1 * tp) + (dalphadtn * dtndt * tp2 * tp)
            dtkdd = (beta * tp1 * dtpdd) + (dbetadtn * dtndd * tp1 * tp) + &
                     (alpha * tp2 * dtpdd) + (dalphadtn * dtndd * tp2 * tp)
         
            tk = beta * (tp1 * tp) + alpha * (tp2 * tp)
            denom = tp * (tk * dtpdd - tp * dtkdd)/(tk * tk * tk)
         end if
         
         ! zeta (chugunov 2007 eq 3)
         U = four_thirds*(tp * tp)/(pi2 * tk * tk)
         zeta = pow(U,x13)
         dzetadt = -x23 * (zeta/tk) * dtkdt 
         dzetadd =  x13 * (zeta/U) * (four_thirds/pi2) * 2.d0 * denom

         ! coulomb coupling parameter gamma (itoh 1979 eq 4)
         gam = z1z2 * qe * qe/(a_av * tk * kerg)
         dgamdt = -(gam / tk) * dtkdt 
         dgamdd = -(gam/(a_av * tk)) * (tk * da_avdd + a_av * dtkdd)

         if (gam >=gamfitlim) then
         
            gam = gamfitlim
            dgamdt = 0d0
            dgamdd = 0d0
            
         else if (gam .le. gamfitlim .and. gam .ge. g0 ) then
            alpha = 0.5d0 * (1d0 - cospi((gam-g0)/deltagam))
            dalphadgam = 0.5d0 * pi * sinpi((gam-g0)/deltagam) * (1d0/deltagam)
            beta = 1.d0 - alpha
            dbetadgam = -dalphadgam
            
            dgamdt = beta * dgamdt + (dbetadgam * dgamdt) * gam + (dalphadgam * dgamdt) * gamfitlim 
            dgamdd = beta * dgamdd + (dbetadgam * dgamdd) * gam + (dalphadgam * dgamdd) * gamfitlim 
            
            gam = beta * gam + alpha * gamfitlim
            
         end if
         
         ! coefficients of zeta dependent on the ion coulomb coupling
         ! parameter (gam) (chugunov 2007, just after eq 21)
   
         beta = 0.41d0 - 0.6d0/gam
         dbetadt = 0.6d0/(gam * gam) * dgamdt
         dbetadd = 0.6d0/(gam * gam) * dgamdd
         
         gama = 0.06d0 + 2.2d0/gam
         dgamadt = -2.2d0/(gam * gam) * dgamdt
         dgamadd = -2.2d0/(gam * gam) * dgamdd
   
         ! gamma tilda (chugunov 2007 eq 21)
         zeta2 = zeta * zeta
         zeta3 = zeta2 * zeta
   
         U = (1d0 + alfa * zeta + beta * zeta2 + gama * zeta3)
         denom = pow(U,-x13)
         denom2 = denom/U ! pow(U,-x43)
         ddenomdt = -x13 * denom2 * &
                     ((alfa * dzetadt)  + (zeta2 * dbetadt + 2d0 * dzetadt * beta * zeta) + &
                     (zeta3 * dgamadt + 3d0 * dzetadt * gama * zeta2 ))
         ddenomdd = -x13 * denom2 * &
                     ((alfa * dzetadd)  + (zeta2 * dbetadd + 2d0 * dzetadd * beta * zeta) + &
                     (zeta3 * dgamadd + 3d0 * dzetadd * gama * zeta2 ))
                     
         gamtild  = gam * denom
         dgamtilddt =  ddenomdt * gam + dgamdt * denom
         dgamtilddd =  ddenomdd * gam + dgamdd * denom
   
         gamtild2 = gamtild * gamtild
         ! for for mean field potential H (chugunov 2007 eq 19)
   
         A = gamtild * sqrt(gamtild)
         dAdt = x32 * (A/gamtild) * dgamtilddt
         dAdd = x32 * (A/gamtild) * dgamtilddd
         
         B = c_a1/sqrt(c_a2 + gamtild)
         dBdt = -1d0/2d0 * B/(c_a2 + gamtild) * dgamtilddt
         dBdd = -1d0/2d0 * B/(c_a2 + gamtild) * dgamtilddd
         
         C = c_a3/(1d0 + gamtild)
         dCdt = -C/(1d0 + gamtild) * dgamtilddt
         dCdd = -C/(1d0 + gamtild) * dgamtilddd
         
         t1 = A * ( B + C )
         dt1dt = A * ( dBdt + dCdt ) + dAdt * ( B + C )
         dt1dd = A * ( dBdd + dCdd ) + dAdd * ( B + C )
         
         denom = c_b2 + gamtild
         ddenomdt = dgamtilddt
         ddenomdd = dgamtilddd
         t2 = c_b1 * gamtild2 / denom
         dt2dt = (( denom * c_b1 * 2d0 * gamtild * dgamtilddt ) - (c_b1 * gamtild2 * ddenomdt)) / (denom * denom)
         dt2dd = (( denom * c_b1 * 2d0 * gamtild * dgamtilddd ) - (c_b1 * gamtild2 * ddenomdd))/ (denom * denom)
         
         denom  = c_b4 + gamtild2
         ddenomdt = 2d0 * gamtild * dgamtilddt
         ddenomdd = 2d0 * gamtild * dgamtilddd
         t3 = c_b3 * gamtild2 / denom
         dt3dt = (( denom * c_b3 * 2d0 * gamtild * dgamtilddt ) - (c_b3 * gamtild2 * ddenomdt)) / (denom * denom)
         dt3dd = (( denom * c_b3 * 2d0 * gamtild * dgamtilddd ) - (c_b3 * gamtild2 * ddenomdd)) / (denom * denom)
         
         h0fit = t1 + t2 + t3
         dh0fitdt = dt1dt + dt2dt + dt3dt
         dh0fitdd = dt1dd + dt2dd + dt3dd
         
         dscreendt = dh0fitdt * screen 
         dscreendd = dh0fitdd * screen 
         
         ! limit screening factor to h0fitlim
         if (h0fit > h0fitlim) then
            h0fit = h0fitlim
            dh0fitdt = 0d0
            dh0fitdd = 0d0
         end if

         screen = exp(h0fit)
         dscreendt = dh0fitdt * screen
         dscreendd = dh0fitdd * screen
   
      end subroutine eval_screen_chugunov
   
   end module screening_chugunov
   
