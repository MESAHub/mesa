! ***********************************************************************
!
!   Copyright (C) 2015-2019  Evan Bauer, Bill Paxton & The MESA Team
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
!
! ***********************************************************************

      module diffusion_support

      use const_def
      use chem_def
      use utils_lib, only: is_bad_num
      use star_private_def

      implicit none

      
      real(dp), parameter :: Xlim = 1d-14
      real(dp), parameter :: tiny_mass = 1d3 ! a kilogram
      real(dp), parameter :: tinyX = 1d-50
      real(dp), parameter :: smallX = 1d-20


      contains

      subroutine get_matrix_coeffs( &
            s, nz, nc, m, nzlo, nzhi, ih1, ihe4, pure_Coulomb, &
            dt, v_advection_max, tiny_C, diffusion_factor, &
            A, X, Z, rho_face, T_face, four_pi_r2_rho_face, &
            xm_face, cell_dm, dm_bar, dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, &
            r_face, r_mid, gamma_T_limit_coeffs_face, alfa_face, &
            rad_accel_face, log10_g_rad, g_rad, &
            min_T_for_radaccel, max_T_for_radaccel, &
            X_init, X_face, C, C_div_X, C_div_X_face, &
            e_ap, e_at, e_ar, e_ax, E_field_face, &
            g_ap, g_at, g_ar, g_ax, g_field_face, &
            v_advection_face, v_total_face, vlnP_face, vlnT_face, v_rad_face, &
            GT_face, D_self_face, AD_face, SIG_face, sigma_lnC, ierr)
         
         type (star_info), pointer :: s
         integer, intent(in) :: &
            nz, nc, m, nzlo, nzhi, ih1, ihe4
         logical, intent(in) :: pure_Coulomb
         real(dp), intent(in) :: &
            dt, v_advection_max, tiny_C, min_T_for_radaccel, max_T_for_radaccel
         real(dp), dimension(:), intent(in) :: &
            diffusion_factor, A, rho_face, T_face, four_pi_r2_rho_face, &
            xm_face, cell_dm, dm_bar, dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, &
            r_face, r_mid, gamma_T_limit_coeffs_face, alfa_face
         real(dp), dimension(:,:), intent(in) :: &
            Z, X_init, rad_accel_face, log10_g_rad, g_rad
         real(dp), dimension(:,:), intent(inout) :: X

         real(dp), dimension(:), intent(out) :: AD_face, &
            e_ap, e_at, e_ar, E_field_face, &
            g_ap, g_at, g_ar, g_field_face
         real(dp), dimension(:,:), intent(out) :: &
            X_face, C, C_div_X, C_div_X_face, v_advection_face, v_total_face, &
            vlnP_face, vlnT_face, v_rad_face, GT_face, e_ax, g_ax, D_self_face
         real(dp), dimension(:,:,:), intent(out) :: SIG_face, sigma_lnC
         integer, intent(out) :: ierr
         
         integer :: i, j, jj, k, op_err, im
         real(dp) :: dv_im, alfa, beta, cc, tmp, tinyX, dlamch, sfmin, &
            AD_dm_full_on, AD_dm_full_off, AD_boost_factor, sum_dm, &
            Vlimit_dm_full_on, Vlimit_dm_full_off, Vlimit, sigmax, &
            SIG_factor, GT_factor
         real(dp), dimension(m) :: C_face, Z_face, dC_dr_face
         real(dp), dimension(nc) :: total_diffusion_factor
         real(dp) :: dlnne_dr_face
         
         include 'formats'
         
         ierr = 0
         sfmin = dlamch('S')  
            
         tinyX = 1d-50
         do k=nzlo,nzhi
            do j=1,nc
               X(j,k) = max(X(j,k),tinyX)
            end do
            tmp = sum(Z(1:nc,k)*X(1:nc,k)/A(1:nc))
            do j=1,nc
               C_div_X(j,k) = 1d0/(A(j)*tmp)
               C(j,k) = X(j,k)*C_div_X(j,k)
            end do
            C(m,k) = 1d0
            X(m,k) = A(m)/dot_product(A(1:nc),C(1:nc,k)) 
         end do
         
         Vlimit_dm_full_on = s% diffusion_Vlimit_dm_full_on*Msun
         Vlimit_dm_full_off = s% diffusion_Vlimit_dm_full_off*Msun
         Vlimit = s% diffusion_Vlimit

         
!$OMP PARALLEL DO PRIVATE(k, j, i, total_diffusion_factor, op_err, C_face, Z_face, dC_dr_face, dlnne_dr_face, tmp) SCHEDULE(dynamic,2)

         do k = nzlo+1, nzhi

            ! Total diffusion scaling factor is product of
            ! diffusion_class_factor (const throughout the star), and
            ! extra_diffusion_factor (profile set in other_diffusion_factor)
            if(s% use_other_diffusion_factor) then
               do j=1,nc
                  total_diffusion_factor(j) = diffusion_factor(j)*s% extra_diffusion_factor(j,k)
               end do
            else
               total_diffusion_factor(1:nc) = diffusion_factor(1:nc)
            end if
            
            call get1_CXZn_face( &
               k, nz, nc, m, nzlo, nzhi, C, X, Z, A, alfa_face, tiny_C, &
               four_pi_r2_rho_face(k)/dm_bar(k), dlnRho_dr_face(k), &
               C_face, X_face, Z_face, C_div_X_face, dC_dr_face, dlnne_dr_face)
            
            op_err = 0
            call get1_coeffs_face( &
               s, k, nz, nc, m, nzlo, nzhi, ih1, ihe4, pure_Coulomb, &
               rho_face(k), T_face(k), gamma_T_limit_coeffs_face(k), &
               four_pi_r2_rho_face(k), dm_bar(k), v_advection_max, tiny_C, sfmin, &
               dlnP_dr_face(k), dlnT_dr_face(k), dlnRho_dr_face(k), &
               s% grav(k), dlnne_dr_face, &
               Vlimit_dm_full_on, Vlimit_dm_full_off, Vlimit, xm_face(k), r_mid, dt, &
               A, X_face(:,k), Z_face, C_face, C_div_X_face(:,k), &
               total_diffusion_factor, rad_accel_face(:,k), &
               s% diffusion_use_cgs_solver, s% eta(k), &
               s% cgs_thermal_diffusion_eta_full_on, s% cgs_thermal_diffusion_eta_full_off, &
               (T_face(k) <= max_T_for_radaccel .and. T_face(k) >= min_T_for_radaccel), &
               v_advection_face(:,k), vlnP_face(:,k), vlnT_face(:,k), v_rad_face(:,k), &
               e_ap(k), e_at(k), e_ar(k), e_ax(:,k), &
               g_ap(k), g_at(k), g_ar(k), g_ax(:,k), &
               sigma_lnC(:,:,k), op_err)
            if (op_err /= 0) ierr = op_err

            if(s% diffusion_use_cgs_solver) then
               ! Electric Field from Iben & MacDonald solve.
               E_field_face(k) = &
                    boltzm*T_face(k)*dlnT_dr_face(k)*e_at(k) + &
                    amu*(s% grav(k))*e_ap(k) ! + e_ar(k) should be the radiation term, skipping for now
               do j=1,nc
                  E_field_face(k) = E_field_face(k) + &
                       boltzm*T_face(k)*e_ax(j,k)*(dC_dr_face(j)&
                       &/C_face(j) + dlnne_dr_face)
               end do
               ! Because we solved for eE as the unknown rather than E in the vector containing diffusion velocities.
               E_field_face(k) = E_field_face(k)/qe
               g_field_face(k) = s% grav(k)
            else
               ! Electric and gravitational field from Thoul solve.
               E_field_face(k) = &
                    e_ap(k)*dlnP_dr_face(k) + &
                    e_at(k)*dlnT_dr_face(k) ! + e_ar(k)*......  skipping the radiation term for now
               g_field_face(k) = &
                    g_ap(k)*dlnP_dr_face(k) + &
                    g_at(k)*dlnT_dr_face(k) ! + g_ar(k)*......  skipping the radiation term for now
               do j=1,nc
                  if (C_face(j) < 1d-20) cycle
                  E_field_face(k) = E_field_face(k) + &
                       e_ax(j,k)*dC_dr_face(j)/C_face(j)
                  g_field_face(k) = g_field_face(k) + &
                       g_ax(j,k)*dC_dr_face(j)/C_face(j)
               end do
               ! Convert to cgs units
               tmp = sum(Z_face(1:nc)*X_face(1:nc,k)/A(1:nc))
               g_field_face(k) = g_field_face(k)*(1.144d-40)*T_face(k)*tmp/(amu*amu)
               e_field_face(k) = e_field_face(k)*(1.144d-40)*T_face(k)*tmp/(qe*amu)
            end if

            do i = 1,nc
               v_total_face(i,k) = v_advection_face(i,k)
               do j = 1,nc
                  v_total_face(i,k) = v_total_face(i,k) - sigma_lnC(i,j,k)*dC_dr_face(j)/C_face(j)
               end do
            end do
            
         end do

!$OMP END PARALLEL DO

         if (ierr /= 0) return
         sum_dm = cell_dm(nzlo)
         
         AD_dm_full_on = s% diffusion_AD_dm_full_on*Msun
         AD_dm_full_off = s% diffusion_AD_dm_full_off*Msun
         AD_boost_factor = s% diffusion_AD_boost_factor
         
         SIG_factor = s% diffusion_SIG_factor
         GT_factor = s% diffusion_GT_factor
         
         !write(*,1) 'GT_factor SIG_factor', GT_factor, SIG_factor

         do k = nzlo+1, nzhi         
            call get1_flow_coeffs( &
               k, nc, m, v_advection_face(:,k), v_advection_max, &
               SIG_factor, GT_factor, sigma_lnC(:,:,k), &
               four_pi_r2_rho_face(k), dm_bar(k), &
               C_div_X_face(:,k), GT_face(:,k), D_self_face(:,k), SIG_face(:,:,k))            
            if (sum_dm >= AD_dm_full_off) then
               AD_face(k) = 0d0
            else
               sigmax = 0d0
               do j=1,nc
                  if (SIG_face(j,j,k) > sigmax) sigmax = SIG_face(j,j,k)
               end do
               AD_face(k) = AD_boost_factor*sigmax
               if (sum_dm > AD_dm_full_on) &
                  AD_face(k) = AD_face(k) * &
                     (sum_dm - AD_dm_full_off)/&
                        (AD_dm_full_on - AD_dm_full_off)
               !write(*,2) 'boost factor AD_face', k, AD_face(k)/sigmax, AD_face(k)
            end if
            sum_dm = sum_dm + cell_dm(k)              
         end do
         
         do j=1,nc ! not used, but copy just for sake of plotting
            D_self_face(j,nzlo) = D_self_face(j,nzlo+1)
            v_advection_face(j,nzlo) = v_advection_face(j,nzlo+1)
            v_total_face(j,nzlo) = v_total_face(j,nzlo+1)
            vlnP_face(j,nzlo) = vlnP_face(j,nzlo+1)
            vlnT_face(j,nzlo) = vlnT_face(j,nzlo+1)
            v_rad_face(j,nzlo) = v_rad_face(j,nzlo+1)
            GT_face(j,nzlo) = GT_face(j,nzlo+1)
            do i=1,nc
               SIG_face(i,j,nzlo) = SIG_face(i,j,nzlo+1)
            end do
         end do
                              
      end subroutine get_matrix_coeffs
      
      
      subroutine get1_coeffs_face( &
            s, k, nz, nc, m, nzlo, nzhi, ih1, ihe4, pure_Coulomb, &
            rho_face, T_face, gamma_T_limit_coeff_face, &
            four_pi_r2_rho_face, dm_bar, &
            v_advection_max, tiny_C, sfmin, &
            dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, grav, dlnne_dr_face, &
            Vlimit_dm_full_on, Vlimit_dm_full_off, Vlimit, xm_face, r_mid, dt, &
            A, X_face, Z_face, C_face, C_div_X_face, &
            diffusion_factor, rad_accel_face, &
            use_cgs_solver, eta, eta_on, eta_off, rad, &
            v_advection_face, vlnP_face, vlnT_face, v_rad_face, &
            e_ap, e_at, e_ar, e_ax, &
            g_ap, g_at, g_ar, g_ax, &
            sigma_lnC, ierr)
            
         type (star_info), pointer :: s
         integer, intent(in) :: k, nz, nc, m, nzlo, nzhi, ih1, ihe4
         logical, intent(in) :: pure_Coulomb
         real(dp), intent(in) :: rho_face, T_face, gamma_T_limit_coeff_face, &
            xm_face, r_mid(:), Vlimit_dm_full_on, Vlimit_dm_full_off, Vlimit, dt
         real(dp), intent(in) :: four_pi_r2_rho_face, dm_bar
         real(dp), intent(in) :: v_advection_max, tiny_C, sfmin
         real(dp), intent(in) :: dlnP_dr_face, dlnT_dr_face, &
              dlnRho_dr_face, grav, dlnne_dr_face
         real(dp), intent(in), dimension(:) :: &
            A, X_face, Z_face, C_face, C_div_X_face, &
            rad_accel_face, diffusion_factor 
         logical, intent(in) :: use_cgs_solver, rad
         real(dp), intent(in) :: eta, eta_on, eta_off
         real(dp), intent(inout), dimension(:) :: &
            v_advection_face, vlnP_face, vlnT_face, v_rad_face ! (nc)
         real(dp), intent(out) :: e_ap, e_at, e_ar, g_ap, g_at, g_ar
         real(dp), intent(inout) :: e_ax(:), g_ax(:) ! (m)
         real(dp), intent(inout) :: sigma_lnC(:,:) ! (nc,nc)
         integer, intent(out) :: ierr
         
         integer :: i, j
         real(dp), dimension(m) :: AP, AT, AR
         real(dp), dimension(m,m) :: kappa_st, Zdiff, Zdiff1, Zdiff2, AX
         
         include 'formats'
         
         ierr = 0
            
         call get1_burgers_coeffs( &
            s, k, nc, m, A, Z_face, X_face, C_face, &
            rho_face, T_face, pure_Coulomb, &
            kappa_st, Zdiff, Zdiff1, Zdiff2)
         
         call get1_gradient_coeffs( &
            k, m, sfmin, A, Z_face, X_face, C_face, rho_face, T_face, &
            use_cgs_solver, eta, eta_on, eta_off, &
            rad, rad_accel_face, kappa_st, Zdiff, Zdiff1, Zdiff2, &
            AP, AT, AR, AX, &
            e_ap, e_at, e_ar, e_ax, &
            g_ap, g_at, g_ar, g_ax, &
            ierr)
         if (ierr /= 0) return
         
         call get1_diffusion_velocities( &
            k, nc, m, nzlo, nzhi, AP, AT, AR, AX, rho_face, T_face, &
            dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, &
            grav, dlnne_dr_face, X_face, &
            Vlimit_dm_full_on, Vlimit_dm_full_off, Vlimit, xm_face, r_mid, s% dt, &
            gamma_T_limit_coeff_face, v_advection_max, diffusion_factor, &
            use_cgs_solver, &
            v_advection_face, vlnP_face, vlnT_face, v_rad_face, sigma_lnC)  
         
      end subroutine get1_coeffs_face

            
      subroutine get1_CXZn_face( &
            k, nz, nc, m, nzlo, nzhi, C, X, Z, A, alfa_face, tiny_C, &
            d_dr_factor, dlnRho_dr_face, C_face, X_face, Z_face, C_div_X_face, &
            dC_dr_face, dlnne_dr_face)
         integer, intent(in) :: k, nc, m, nz, nzlo, nzhi         
         real(dp), dimension(:,:), intent(in) :: C, X, Z ! (m,nz)
         real(dp), intent(in) :: A(:) ! (m) atomic number
         real(dp), intent(in) :: alfa_face(:), d_dr_factor, dlnRho_dr_face
         real(dp), intent(in) :: tiny_C
         real(dp), dimension(:), intent(out) :: C_face, Z_face, dC_dr_face ! (m)
         real(dp), intent(out) :: dlnne_dr_face
         real(dp), dimension(:,:), intent(out) :: X_face, C_div_X_face ! (m,nz)
         integer :: j
         real(dp) :: tmp, tmp1, tmp2, dlntmp_dr_face, alfa, beta

         alfa = alfa_face(k)
         beta = 1d0 - alfa
         do j = 1, m
            X_face(j,k) = alfa*X(j,k) + beta*X(j,k-1)
            Z_face(j) = alfa*Z(j,k) + beta*Z(j,k-1)
         end do

         ! "tmp" at the face of the zone
         tmp = sum(Z_face(1:nc)*X_face(1:nc,k)/A(1:nc))
         do j = 1, m
            C_div_X_face(j,k) = 1/(A(j)*tmp)
            C_face(j) = X_face(j,k)*C_div_X_face(j,k)

            ! Old way of calculating the derivative. Fixed below.
            ! dC_dr_face(j) = (X(j,k-1) - X(j,k))*C_div_X_face(j,k)*d_dr_factor

            dC_dr_face(j) = (C(j,k-1) - C(j,k))*d_dr_factor
         end do

         ! Calculate the electron number density ln gradient used in
         ! converting between the Thoul and Iben/MacDonald notations.
         ! This is accomplished by adding two logarithmic derivatives
         ! evaluated at the face of the zone, one of which we already
         ! know, and one of which must be calculated.

         ! "tmp" averaged over zone k
         tmp1 = sum(Z(1:nc,k)*X(1:nc,k)/A(1:nc))
         ! "tmp" averaged over zone k-1
         tmp2 = sum(Z(1:nc,k-1)*X(1:nc,k-1)/A(1:nc))
         dlntmp_dr_face = ((tmp2-tmp1)/tmp)*d_dr_factor
         ! ne = rho*tmp/amu, so...
         dlnne_dr_face = dlnRho_dr_face + dlntmp_dr_face

      end subroutine get1_CXZn_face
      
      
      subroutine get1_burgers_coeffs( &
            s, k, nc, m, A, Z, X, C, rho, T, pure_Coulomb, &
            kappa_st, Zdiff, Zdiff1, Zdiff2)
         
         use paquette_coeffs, only: paquette_coefficients
         
         type (star_info), pointer :: s
         integer, intent(in) :: k, nc, m
         real(dp), intent(in) :: rho, T
         logical, intent(in) :: pure_Coulomb
         real(dp), intent(in), dimension(:) :: A, X, Z, C ! (m)
         real(dp), intent(inout), dimension(:,:) :: &
            kappa_st, Zdiff, Zdiff1, Zdiff2 ! (m,m)

         integer :: i, j
         real(dp) :: ac, ni, cz, xij, ne, ao, lambdad, lambda
         real(dp), dimension(m) :: charge, na
         real(dp), dimension(m,m) :: cl, Ath, Ddiff, Kdiff
         real(dp) :: kappa_SM ! For diagnosis right now.
            
         do i = 1, nc
            charge(i) = max(1d0, Z(i)) ! assume some ionization
         end do
         charge(m) = Z(m)
         
         if (.not. pure_Coulomb) then ! use Paquette coeffs
            ! Get number densities (per cm^3)
            do i = 1, nc
               na(i) = rho*X(i)/(A(i)*amu)   
            end do         
            na(m) = 0.d0      
            do i = 1, nc
               na(m) = na(m) + charge(i)*na(i)
            end do
            ! Compute resistance coefficients from Paquette&al (1986)   
            call paquette_coefficients( &
               rho, T, m, A, charge, na, Ddiff, Kdiff, Zdiff, Zdiff1, Zdiff2, Ath)

            if(.not. s% diffusion_use_paquette .and. .not. s% use_other_diffusion_coefficients) then
               call get_SM_coeffs(nc,m,rho,T,A,charge,na,Kdiff,Zdiff,Zdiff1,Zdiff2,kappa_SM)
               ! This must get called after Paquette because it doesn't calculate
               ! the electron entries (m). It leaves them untouched while calculating
               ! and changing all the ion-ion terms (1:nc), so after calling this
               ! routine we have all ion-ion coefficients from Stanton&Murillo and
               ! all ion-electron coefficients from Paquette&al.
            end if

            if(s% use_other_diffusion_coefficients) then
               call s% other_diffusion_coefficients( &
                  s% id, k, nc, m, rho, T,  A, X, Z, C, charge, na, &
                  Ddiff, Kdiff, Zdiff, Zdiff1, Zdiff2, Ath)
            end if

            ! Unit conversion conveniently applies to both Paquette and Stanton&Murillo
            kappa_st(:,:) = Kdiff(:,:)/(1.41D-25*pow(T,-1.5D0)*na(m)*na(m))
               ! = kappa_st of eq 37, Thoul&al 1994 
            return
         end if
         
         ! calculate density of electrons (ne) from mass density (rho):
         ac=0.d0
         do i=1, m
            ac=ac+a(i)*c(i)
         end do   
         ne=rho/(mp*ac) 
         ! calculate interionic distance (ao): 
         ni=0.d0
         do i=1, nc
            ni=ni+c(i)*ne
         end do
         ao=pow(0.23873d0/ni,one_third) 
         ! calculate debye length (lambdad):
         cz=0.d0
         do i=1, m
            cz=cz+c(i)*charge(i)*charge(i)
         end do
         lambdad=6.9010d0*sqrt(t/(ne*cz))
         ! calculate lambda to use in coulomb logarithm:
         lambda=max(lambdad, ao)
         ! calculate coulomb logarithms:
         do i=1, m
            do j=1, m
               xij=2.3939d3*t*lambda/abs(charge(i)*charge(j))
               cl(i,j)=0.81245d0*log1p(0.18769d0*pow(xij,1.2d0))
            end do
         end do

         ! set coeffs for pure Coulomb potential
         do i=1, m
            do j=1, m
               Zdiff(i,j) = 0.6d0
               Zdiff1(i,j) = 1.3d0
               Zdiff2(i,j) = 2d0
               kappa_st(i,j) = &
                  cl(i,j)*sqrt(a(i)*a(j)/(a(i)+a(j)))* &
                     c(i)*c(j)*charge(i)*charge(i)*charge(j)*charge(j)
            end do
         end do
         
      end subroutine get1_burgers_coeffs


      subroutine get1_gradient_coeffs( &
            k, m, sfmin, A, Z, X, C, rho, T, use_cgs_solver, eta, eta_on, eta_off, &
            rad, rad_accel, kappa_st, Zdiff, Zdiff1, Zdiff2, &
            AP, AT, AR, AX, &
            e_ap, e_at, e_ar, e_ax, &
            g_ap, g_at, g_ar, g_ax, &
            ierr)
         integer, intent(in) :: k, m
         real(dp), intent(in) :: sfmin, rho, T
         real(dp), intent(in) :: eta, eta_on, eta_off
         real(dp), intent(in), dimension(:) :: A, X, Z, C, rad_accel ! (m)
         logical, intent(in) :: use_cgs_solver, rad
         real(dp), dimension(:,:), intent(in) :: &
            kappa_st, Zdiff, Zdiff1, Zdiff2 ! (m,m)
         real(dp), dimension(:), intent(out) :: AP, AT, AR ! (m)
         real(dp), intent(inout) :: AX(:,:) ! (m,m)
         real(dp), intent(inout) :: e_ap, e_at, e_ar, e_ax(:) ! (m)
         real(dp), intent(inout) :: g_ap, g_at, g_ar, g_ax(:) ! (m)
         integer, intent(out) :: ierr
           
         integer :: i, j
         real(dp) :: charge(m), nd(m), Kdiff(m,m), alfa, beta

         ! For blending solves with and without thermal diffusion.
         real(dp) :: AP1(m), AT1(m), AR1(m), AX1(m,m)
         real(dp) :: e_ap1, e_at1, e_ar1, e_ax1(m)
         real(dp) :: AP2(m), AT2(m), AR2(m), AX2(m,m)
         real(dp) :: e_ap2, e_at2, e_ar2, e_ax2(m)
         
         include 'formats'
         
         ierr = 0
            
         do i=1,m-1
            charge(i) = max(1d0, Z(i))
         end do
         charge(m) = Z(m)

         if(use_cgs_solver) then ! Use the cgs solver
            ! Get number densities using info contained in X,A,Z,rho
            nd(1:m) = 0d0
            do i=1,m-1
               nd(i) = rho*X(i)/(A(i)*amu)
               nd(m) = nd(m) + nd(i)*charge(i) ! Electron Number Density satisfies charge neutrality
            end do
            
            Kdiff(:,:) = kappa_st(:,:)*(1.41D-25*pow(T,-1.5D0)*nd(m)*nd(m))

            if(eta < eta_on) then
               call solve_burgers_cgs_with_thermal(2*m+1,m,A,charge,nd,rad_accel,rad, &
                    Kdiff, Zdiff, Zdiff1, Zdiff2, &
                    AP,AT,AR,AX, &
                    e_ap,e_at,e_ar,e_ax,ierr)
            else if(eta > eta_off) then
               call solve_burgers_cgs_no_thermal(m+1,m,A,charge,nd,rad_accel,rad, &
                    Kdiff,AP,AT,AR,AX, &
                    e_ap,e_at,e_ar,e_ax,ierr)
            else
               ! Call both and do a linear blend of all coefficients.
               alfa = (eta - eta_on)/(eta_off - eta_on) ! alfa = 1 means no thermal diffusion.
               beta = 1d0 - alfa ! beta = 1 means full thermal diffusion.
               
               call solve_burgers_cgs_no_thermal(m+1,m,A,charge,nd,rad_accel,rad, &
                    Kdiff,AP1,AT1,AR1,AX1, &
                    e_ap1,e_at1,e_ar1,e_ax1,ierr)

               if (ierr /= 0) then
                  !return      
                  write(*,2) 'solve_burgers_cgs_no_thermal failed', k
                  do i=1,m-1
                     write(*,2) 'A X Z C', i, A(i), X(i), Z(i), C(i)
                  end do
                  stop 'get1_gradient_coeffs'
               end if

               call solve_burgers_cgs_with_thermal(2*m+1,m,A,charge,nd,rad_accel,rad, &
                    Kdiff, Zdiff, Zdiff1, Zdiff2, &
                    AP2,AT2,AR2,AX2, &
                    e_ap2,e_at2,e_ar2,e_ax2,ierr)

               ! Check how much different the two solutions are.
               ! if( abs((AT1(3) - AT2(3))/AT1(3)) > 2d0 ) then ! 3 is index for Helium in basic.net
               !    print *, "Thermal diffusion changing temperature coefficient by more than factor of two."
               !    print *, "Relative difference: ", abs((AT1(3) - AT2(3))/AT1(3))
               ! end if
               
               ! Blending between the two solutions.
               do i = 1,m
                  AP(i) = alfa*AP1(i) + beta*AP2(i)
                  AT(i) = alfa*AT1(i) + beta*AT2(i)
                  AR(i) = alfa*AR1(i) + beta*AR2(i)
                  do j = 1,m
                     AX(i,j) = alfa*AX1(i,j) + beta*AX2(i,j)
                  end do
                  e_ax(i) = alfa*e_ax1(i) + beta*e_ax2(i)
               end do
               e_ap = alfa*e_ap1 + beta*e_ap2
               e_at = alfa*e_at1 + beta*e_at2
               e_ar = alfa*e_ar1 + beta*e_ar2

            end if
            ! Gravity isn't being calculated by this version of the
            ! diffusion routine, so set the coefficients to zero.
            g_ap = 0d0
            g_at = 0d0
            g_ar = 0d0
            g_ax(1:m) = 0d0
            
            if (ierr /= 0) then
               !return      
               write(*,2) 'solve_burgers_cgs failed', k
               do i=1,m-1
                  write(*,2) 'A X Z C', i, A(i), X(i), Z(i), C(i)
               end do
               stop 'get1_gradient_coeffs'
            end if
         else ! Use the Thoul solver
            call do1_solve_thoul_hu( &
               2*m+2, m, sfmin, A, charge, X, C, rad_accel, rad, &
               kappa_st, Zdiff, Zdiff1, Zdiff2, &
               AP, AT, AR, AX, &
               e_ap, e_at, e_ar, e_ax, &
               g_ap, g_at, g_ar, g_ax, &
               ierr)

            if (ierr /= 0) then
               !return      
               write(*,2) 'do1_solve_thoul_hu failed', k
               do i=1,m-1
                  write(*,2) 'A X Z C', i, A(i), X(i), Z(i), C(i)
               end do
               stop 'get1_gradient_coeffs'
            end if
         end if

         
      end subroutine get1_gradient_coeffs
      
      
      subroutine get1_diffusion_velocities( &
            k, nc, m, nzlo, nzhi, AP, AT, AR, AX, rho, T, &
            dlnP_dr, dlnT_dr, dlnRho_dr, grav, dlnne_dr, X_face, &
            Vlimit_dm_full_on, Vlimit_dm_full_off, Vlimit, xm_face, r_mid, dt, &
            limit_coeff, v_advection_max, diffusion_factor, use_cgs_solver, &
            vgt, vlnP, vlnT, vrad, sigma_lnC)
         integer, intent(in) :: k, nc, m, nzlo, nzhi
         real(dp), intent(in), dimension(:) :: AP, AT, AR, r_mid
         real(dp), intent(in) :: AX(:,:) ! (m,m)
         real(dp), intent(in) :: rho, T, limit_coeff, v_advection_max, &
            Vlimit_dm_full_on, Vlimit_dm_full_off, Vlimit, xm_face, dt
         real(dp), intent(in) :: dlnP_dr, dlnT_dr, dlnRho_dr, grav, dlnne_dr
         real(dp), intent(in) :: X_face(:)
         real(dp), intent(in) :: diffusion_factor(:)
         logical, intent(in) :: use_cgs_solver
         real(dp), intent(inout), dimension(:) :: vgt, vlnP, vlnT, vrad
         real(dp), intent(inout) :: sigma_lnC(:,:) ! (nc,nc)
         
         integer :: i, j, im
         real(dp) :: coef, coef_vrad, dv_im, dr, T2pt5, &
            vcross, vmax, vmax_limit, frac, alfa, beta
         real(dp) :: tau0  ! = 6d13*secyer, solar diffusion time (seconds)
         real(dp), parameter :: rho_unit = 1d2
         real(dp), parameter :: T_unit = 1d7
         
         include 'formats'

         if (limit_coeff <= 0) then
            vgt(:) = 0
            sigma_lnC(:,:) = 0
            return
         end if
         
         dr = r_mid(k-1) - r_mid(k)
         vcross = dr/dt
         if (xm_face >= Vlimit_dm_full_off .or. Vlimit <= 0d0) then
            vmax_limit = 1d99
            alfa = 0d0
            beta = 1d0
         else if (xm_face <= Vlimit_dm_full_on) then
            vmax_limit = vcross*Vlimit
            alfa = 1d0
            beta = 0d0
         else ! combine
            alfa = (xm_face - Vlimit_dm_full_off)/&
                        (Vlimit_dm_full_on - Vlimit_dm_full_off)
            beta = 1d0 - alfa ! fraction of normal v when it is > vmax
            vmax_limit = vcross*Vlimit/alfa ! Want to scale to no limit at alfa = 0
         end if
         
         if(use_cgs_solver) then ! Converts coefficients to velocities
            ! assuming cgs routine.
            do i=1,nc
               do j=1,nc
                  sigma_lnC(j,i) = -limit_coeff*diffusion_factor(i)*boltzm*T*AX(j,i)
               end do
               vlnP(i) = AP(i)*amu*grav*diffusion_factor(i)*limit_coeff
               vlnT(i) = AT(i)*boltzm*T*dlnT_dr*diffusion_factor(i)*limit_coeff
               vrad(i) = AR(i)*diffusion_factor(i)*limit_coeff ! AR already contains all constants.
               vgt(i) = vlnP(i) + vlnT(i) + vrad(i)
            end do
            
            do i = 1,nc
               ! Converting from Iben/MacDonald notation to Thoul
               ! notation using electron number density gradient.
               vgt(i) = vgt(i) - dlnne_dr*sum(sigma_lnC(i,1:nc))
               if (X_face(i) < 1d-15) vgt(i) = 0d0
            end do
            
         else ! converts coefficients to velocities assuming Thoul.
            tau0 = 6d13*secyer
            T2pt5 = pow(T/T_unit,2.5d0)
            coef = limit_coeff*Rsun*T2pt5/(rho/rho_unit)*(Rsun/tau0)
            coef_vrad = (limit_coeff/T)*Rsun*Rsun*T2pt5/(rho/rho_unit)/tau0
            ! converts to cgs units
            ! T^(5/2)/rho takes care of the last part from Thoul equations
            ! (21) and (28). Rsun/tau0 converts the velocity from units of
            ! Rsun/tau0 to cm/s. The extra Rsun accounts for the fact that
            ! MESA's gradients (e.g. dlnP_dr) have units of cm^-1, whereas
            ! Thoul's gradients in eqn (21) are Rsun^-1.

            do i=1,nc
               do j=1,nc
                  sigma_lnC(j,i) = -diffusion_factor(i)*coef*AX(j,i)
               end do
               vlnP(i) = AP(i)*dlnP_dr*diffusion_factor(i)*coef
               vlnT(i) = AT(i)*dlnT_dr*diffusion_factor(i)*coef
               vrad(i) = AR(i)*diffusion_factor(i)*coef_vrad
               vgt(i) = vlnP(i) + vlnT(i) + vrad(i)
               if (X_face(i) < 1d-15) vgt(i) = 0d0
            end do
         end if
         
         ! final fixup for vgt of most abundant so it gives baryon conservation.
         im = maxloc(X_face(1:nc),dim=1)
         dv_im = -dot_product(X_face(1:nc), vgt(1:nc))/X_face(im)
         vgt(im) = vgt(im) + dv_im
         
         vmax = maxval(abs(vgt(1:nc)))
         if (vmax > v_advection_max) then
            frac = v_advection_max/vmax
            do i=1,nc
               vgt(i) = vgt(i)*frac
               ! Need to rescale the sigma terms by the same factor
               ! or diffusive equilibrium can break in WDs.
               do j=1,nc
                  sigma_lnC(j,i) = sigma_lnC(j,i)*frac
               end do
            end do
            vmax = v_advection_max
            !write(*,3) 'vmax > v_advection_max', im, k, vmax, v_advection_max
            !stop 'get1_diffusion_velocities'
         end if

         if (alfa > 0d0 .and. vmax > vmax_limit) then
            frac = vmax_limit/vmax
            do i=1,nc
               vgt(i) = vgt(i)*frac
               ! Need to rescale the sigma terms by the same factor
               ! or diffusive equilibrium can break in WDs.
               do j=1,nc
                  sigma_lnC(j,i) = sigma_lnC(j,i)*frac
               end do
            end do
         end if
         
      end subroutine get1_diffusion_velocities
      
      
      subroutine get1_flow_coeffs( &
            k, nc, m, &
            v_advection_face, v_advection_max, SIG_factor, GT_factor, &
            sigma_lnC_face, four_pi_r2_rho_face, &
            dm_bar, C_div_X_face, GT_face, D_self_face, SIG_face)
         integer, intent(in) :: k, nc, m
         real(dp), intent(in) :: &
            v_advection_max, SIG_factor, GT_factor, v_advection_face(:) ! (nc)
         real(dp), intent(in) :: sigma_lnC_face(:,:) ! (nc,nc)
         real(dp), intent(in) :: four_pi_r2_rho_face, dm_bar
         real(dp), intent(in), dimension(:) :: C_div_X_face ! (m)
         real(dp), intent(inout) :: GT_face(:) ! (nc)
         real(dp), intent(inout) :: D_self_face(:) ! (nc)
         real(dp), intent(inout) :: SIG_face(:,:) ! (nc,nc)
         
         integer :: i, j
         real(dp) :: c, boost
         
         include 'formats'

         c = SIG_factor*four_pi_r2_rho_face*four_pi_r2_rho_face/dm_bar
         do i = 1, nc
            GT_face(i) = GT_factor*four_pi_r2_rho_face*v_advection_face(i)
            D_self_face(i) = sigma_lnC_face(i,i)  
            do j = 1, nc
               SIG_face(i,j) = c*sigma_lnC_face(i,j)/C_div_X_face(j)               
            end do
         end do
         
      end subroutine get1_flow_coeffs


!*************************************************************
! Original of this routine was written by Anne A. Thoul, at the Institute
! for Advanced Study, Princeton, NJ 08540.
! See Thoul et al., Ap.J. 421, p. 828 (1994)

! With modifications by Hali Hu for non Coulomb and rad levitation.
!*************************************************************
! This routine inverses the burgers equations.
!
! The system contains N equations with N unknowns. 
! The equations are: the M momentum equations, 
!                    the M energy equations, 
!                    two constraints: the current neutrality 
!                                     the zero fluid velocity.
! The unknowns are: the M diffusion velocities,
!                   the M heat fluxes,
!                   the electric field E
!                   the gravitational force g.
!
!**************************************************


! comments from Anne's original version
   ! Inverse the system for each possible right-hand-side, i.e.,
   ! if alpha is the r.h.s., we obtain the coefficient A_p
   ! if nu    ---------------------------------------- A_T
   ! if gamma(i,j) ----------------------------------- A_Cj
   ! 
   ! If I=1, we obtain the hydrogen diffusion velocity
   ! If I=2, ------------- helium   ------------------
   ! If I=3,M-1, --------- heavy element -------------
   ! If I=M, ------------- electrons -----------------
   ! For I=M,2M, we get the heat fluxes
   ! For I=N-1, we get the electric field
   ! For I=N, we get the gravitational force g



      subroutine do1_solve_thoul_hu( &
            n, m, sfmin, a, z, x, c, rad_accel, rad, &
            kappa_st, Zdiff, Zdiff1, Zdiff2, &
            ap, at, ar, ax, &
            e_ap, e_at, e_ar, e_ax, &
            g_ap, g_at, g_ar, g_ax, &
            ierr)

         ! the parameter m is the number of fluids considered (ions+electrons)
         ! the parameter n is the number of equations (2*m+2).
         !
         ! the vectors a,z and x contain the atomic mass numbers, 
         ! the charges (ionization), and the mass fractions, of the elements.
         ! note: since m is the electron fluid, its mass and charge must be
         !      a(m)=m_e/m_u
         !      z(m)=-1.
         !
         ! the array cl contains the values of the coulomb logarithms.
         ! the vector ap, at, and array ax contains the results for the diffusion 
         ! coefficients.

         integer, intent(in) :: m,n
         real(dp), intent(in) :: sfmin
         real(dp), intent(in), dimension(:) :: A, Z, X, C, rad_accel ! (m)
         logical, intent(in) :: rad
         real(dp), intent(in), dimension(:,:) :: &
            kappa_st, Zdiff, Zdiff1, Zdiff2 ! (m,m)
!           kappa_st from the resistance coefficient Kdiff with eq (37) Thoul&al.
!           Zdiff, Zdiff1, Zdiff2 = arrays of resistance coefficients,
         real(dp), intent(inout), dimension(:) :: ap, at, ar ! (m)
         real(dp), intent(inout) :: ax(:,:) ! (m,m)
         real(dp), intent(inout) :: e_ap, e_at, e_ar, e_ax(:) ! (m)
         real(dp), intent(inout) :: g_ap, g_at, g_ar, g_ax(:) ! (m)
         integer, intent(out) :: ierr

         integer :: i, j, l, indx(n), nmax
         real(dp) :: aamax, cc, ac, temp, ko, d, f
         real(dp), dimension(m,m) :: xx, y, yy, k
         real(dp), dimension(n) :: alpha, nu, ga, beta
         real(dp), dimension(n,n) :: delta, gamma

         ! the vector c contains the concentrations
         ! cc is the total concentration: cc=sum(c_s)
         ! ac is proportional to the mass density: ac=sum(a_s c_s)
         ! the arrays xx,y,yy and k are various parameters which appear in 
         ! burgers equations.
         ! the vectors and arrays alpha, nu, gamma, delta, and ga represent
         ! the "right- and left-hand-sides" of burgers equations, and later 
         ! the diffusion coefficients.
      
         ! initialize:

         ierr = 0
         ko = 2d0  
         indx(1:n) = 0    

         ! calculate cc and ac:
      
         cc=sum(c(1:m))
         ac=dot_product(a(1:m),c(1:m))

         ! calculate the coefficients of the burgers equations

         do i=1,m
            do j=1,m
               xx(i,j)=a(j)/(a(i)+a(j))
               y(i,j)=a(i)/(a(i)+a(j))
               yy(i,j) = 3D0*y(i,j) + Zdiff1(i,j)*xx(i,j)*a(j)/a(i)
               k(i,j) = kappa_st(i,j)
            end do
         end do

         ! write the burgers equations and the two constraints as
         ! alpha_s dp + nu_s dt + sum_t(not ihe or m) gamma_st dc_t 
         !                     = sum_t delta_st w_t

         do i=1,m
            alpha(i)=c(i)/cc
            nu(i)=0d0
            gamma(i,1:n)=0d0
            if (rad) then
               beta(i) = -(amu/boltzm)*alpha(i)*a(i)*rad_accel(i)  
            else
               beta(i) = 0d0
            end if
            do j=1,m
               if (j /= m) then ! HH: Include He gradient
                  gamma(i,j) = -c(j)/cc
                  if (j == i) gamma(i,j) = gamma(i,j) + 1d0
                  gamma(i,j) = gamma(i,j)*c(i)/cc
               end if
            end do
         end do
      
         do i=m+1,n-2
            alpha(i)=0d0
            nu(i)=2.5d0*c(i-m)/cc
            beta(i) = 0d0
            gamma(i,1:n)=0d0
         end do
      
         alpha(n-1)=0d0
         nu(n-1)=0d0
         beta(n-1)=0d0
         gamma(n-1,1:n)=0d0
      
         alpha(n)=0d0
         nu(n)=0d0
         beta(n)=0d0
         gamma(n,1:n)=0d0
      
         delta(1:n,1:n) = 0d0
      
         do i=1,m
         
            do j=1,m
               if (j == i) then
                  do l=1,m
                     if (l /= i) then
                        delta(i,j)=delta(i,j)-k(i,l)
                     end if
                  end do
               else
                  delta(i,j)=k(i,j)
               end if
            end do
         
            do j=m+1,n-2
               if (j-m == i) then
                  do l=1,m
                     if (l /= i) &
                        delta(i,j) = delta(i,j) + Zdiff(i,l)*xx(i,l)*k(i,l)
                  end do
               else
                  delta(i,j) = -Zdiff(i,j-m)*y(i,j-m)*k(i,j-m)
               end if
            end do
         
            delta(i,n-1)=c(i)*z(i)
         
            delta(i,n)=-c(i)*a(i)
            
         end do
      
         do i=m+1,n-2
         
            do j=1,m
               if (j == i-m) then
                  do l=1,m
                     if (l /= i-m) delta(i,j) = &
                        delta(i,j) + 2.5D0*Zdiff(i-m,l)*xx(i-m,l)*k(i-m,l)
                  end do
               else
                  delta(i,j) = -(2.5d0*Zdiff(i-m,j))*xx(i-m,j)*k(i-m,j)
               end if
            end do
         
            do j=m+1,n-2
               if (j-m == i-m) then
                  do l=1,m
                     if (l /= i-m) delta(i,j) = delta(i,j) - &
                           y(i-m,l)*k(i-m,l)*(0.8D0*Zdiff2(i-m,l)*xx(i-m,l)+yy(i-m,l))
                  end do
                  delta(i,j) = delta(i,j) - 0.4D0*Zdiff2(i-m,i-m)*k(i-m,i-m)
               else
                  delta(i,j) = k(i-m,j-m)*xx(i-m,j-m)*y(i-m,j-m) * &
                        (3D0 + Zdiff1(i-m,j-m) - 0.8D0*Zdiff2(i-m,j-m))
               end if
            end do
         
            delta(i,n-1:n)=0d0
            
         end do
      
         do j=1,m
            delta(n-1,j) = c(j)*z(j)
         end do
         delta(n-1,m+1:n) = 0d0
      
         do j=1,m
            delta(n,j) = c(j)*a(j)
         end do
         delta(n,m+1:n) = 0d0
         
         call dgetrf(n, n, delta, n, indx, ierr)
         if (ierr /= 0) return
      
         call dgetrs( 'n', n, 1, delta, n, indx, alpha, n, ierr )
         if (ierr /= 0) return
      
         call dgetrs( 'n', n, 1, delta, n, indx, nu, n, ierr )
         if (ierr /= 0) return
      
         if (rad) then
            call dgetrs( 'n', n, 1, delta, n, indx, beta, n, ierr )
            if (ierr /= 0) return
         end if
      
         do j=1,n
            do i=1,n
               ga(i)=gamma(i,j)
            end do
            call dgetrs( 'n', n, 1, delta, n, indx, ga, n, ierr )
            if (ierr /= 0) return
            do i=1,n
               gamma(i,j)=ga(i)
            end do
         end do
         
         f = ko*ac*cc
         do j=1,m
            ap(j)=alpha(j)*f
            at(j)=nu(j)*f
            ar(j)=beta(j)*f
            do i=1,m
               ax(i,j)=gamma(i,j)*f
            end do
         end do
         
         e_ap=alpha(n-1)*f
         g_ap=alpha(n)*f

         e_at=nu(n-1)*f
         g_at=nu(n)*f

         e_ar=beta(n-1)*f
         g_ar=beta(n)*f
         
         do i=1,m
            e_ax(i)=gamma(n-1,i)*f
            g_ax(i)=gamma(n,i)*f
         end do

      end subroutine do1_solve_thoul_hu




      subroutine solve_burgers_cgs_no_thermal( &
           n, m, A, Z, nd, rad_accel, rad, &
           Kdiff, ap, at, ar, ax, &
           e_ap, e_at, e_ar, e_ax, ierr)
        
        ! nd = array of number densities
        ! m = # of species including electrons
        ! n = m+1 without thermal diffusion
        !   is the number of equations and unkowns
        !   Thermal diffusion off: (m-1) diffusion equations,
        !                          2 conservation equations.
        !   Thermal diffusion  on: 2m diffusion equations (maybe 2*m-1)
        !                          2 conservation equations.
        
        integer, intent(in) :: m,n
        real(dp), intent(in), dimension(:) :: A, Z, nd, rad_accel ! (m)
        logical, intent(in) :: rad
        real(dp), intent(in), dimension(:,:) :: Kdiff ! (m,m)
        real(dp), intent(inout), dimension(:) :: ap, at, ar ! (m)
        real(dp), intent(inout) :: ax(:,:) ! (m,m)
        real(dp), intent(inout) :: e_ap, e_at, e_ar, e_ax(:)
        integer, intent(out) :: ierr

        integer :: i,j,l,indx(n),nmax
        real(dp), dimension(n) :: alpha, beta, nu, ga
        ! Add in beta later for rad lev. ga is the temp holder for the
        ! columns of gamm when doing matrix solve one column at a time.
        real(dp), dimension(n,n) :: delta, gamm

        ! initialize
        ierr = 0
        indx(1:n) = 0
        delta(1:n,1:n) = 0d0
        alpha(1:n) = 0d0
        beta(1:n) = 0d0
        nu(1:n) = 0d0
        ga(1:n) = 0d0
        gamm(1:n,1:n) = 0d0

        e_ap = 0d0
        e_at = 0d0
        e_ar = 0d0
        e_ax(1:m) = 0d0
        
        ! Assign the RHS Matrix multiplying the unkown quantities.
        ! Right now this is for thermal diffusion off, assuming gravity
        ! is a known, so there are m uknown diffusion velocities and
        ! one unkown electric field, corresponding to (m-1) ion equaitions
        ! (no electron equation) and 2 conservation equations.
        do i = 1,(m-1)

           do j = 1,m
              if (j == i) then
                 do l = 1,m
                    if (l /= j) then
                       delta(i,j) = delta(i,j) - Kdiff(i,l)
                    end if
                 end do
              else
                 delta(i,j) = Kdiff(i,j)
              end if
           end do

           delta(i,n) = nd(i)*Z(i)

        end do

        do j = 1,m
           delta(m,j) = A(j)*nd(j)
           delta(n,j) = Z(j)*nd(j)
        end do
        ! Delta assignment finished.

        ! Assign LHS quantities. (Everything not assigned here stays zero thanks to earlier initialization)
        do i = 1,(m-1)
           alpha(i) = nd(i)*A(i)
           nu(i) = nd(i)
           gamm(i,i) = nd(i)
           if (rad) then
              beta(i) = -1d0*alpha(i)*amu*rad_accel(i)
           end if
        end do

        ! alpha,nu,gamm assignment finished


        ! Invert the matrix and solve to get the contribution of each
        ! vector/column on the LHS. This takes alpha, nu, ga from being
        ! the knowns on the LHS to vectors containing information on
        ! contributions for the unknowns. These then get assigned to the
        ! output arrays ap, at, ax, which will then be multiplied by
        ! gravity, k_b*T, and number density gradients to provide the
        ! full solution for diffusion velocities.

        call dgetrf(n,n,delta,n,indx,ierr)

        if( ierr /= 0 ) then
           ! print *, "Factoring failed!"
           return
        end if
        
        call dgetrs('N',n,1,delta,n,indx,alpha,n,ierr)
        if( ierr /= 0 ) then
           ! print *, "solve failed on alpha"
           return
        end if

        call dgetrs('N',n,1,delta,n,indx,nu,n,ierr)
        if( ierr /= 0 ) then
           ! print *, "solve failed on nu"
           return
        end if

        if (rad) then
           call dgetrs('N',n,1,delta,n,indx,beta,n,ierr)
           if( ierr /= 0 ) then
              ! print *, "solve failed on beta"
              return
           end if
        else
           beta(1:n) = 0d0
        end if
        
        do j=1,n
           do i=1,n
              ga(i) = gamm(i,j)
           end do
           call dgetrs('N',n,1,delta,n,indx,ga,n,ierr)
           if( ierr /= 0 ) then
              ! print *, "solve failed on gamma"
              return
           end if
           do i=1,n
              gamm(i,j) = ga(i)
           end do
        end do
                
        ! Assign the results of the matrix solve to the output
        ! arrays/matrix.
        
        do j = 1,m
           ap(j) = alpha(j)
           at(j) = nu(j)
           ar(j) = beta(j)
           do i=1,m
              ax(i,j) = gamm(i,j)
           end do
        end do
        
        e_ap = alpha(n)
        e_at = nu(n)
        e_ar = beta(n)
        do i=1,m
           e_ax(i) = gamm(n,i)
        end do
        
      end subroutine solve_burgers_cgs_no_thermal


      subroutine solve_burgers_cgs_with_thermal( &
           n, m, A, Z, nd, rad_accel, rad, &
           Kdiff, zdiff, zdiff1, zdiff2, &
           ap, at, ar, ax, &
           e_ap, e_at, e_ar, e_ax, ierr)
        
        ! nd = array of number densities
        ! m = # of species including electrons
        ! n = 2*m+1 with thermal diffusion
        !   is the number of equations and unkowns
        !   Thermal diffusion on: 2m-1 diffusion equations (m-1 momentum + m energy)
        !                         2 conservation equations.
        ! There is one less momentum equation because we are neglecting the electron
        ! equation (which is wrong when electrons are degenerate),
        ! and making up for losing that equation by treating grav as known.
        ! We can't do the same with the energy equation because there are no more
        ! unknowns that can be dropped, so this solver shouldn't really be used when
        ! there is degeneracy. Having it in this form makes it easier to transition
        ! between ideal gas (where this solver is valid) and solve_burgers_cgs
        ! (which is much better when things are degenerate).
        
        integer, intent(in) :: m,n
        real(dp), intent(in), dimension(:) :: A, Z, nd, rad_accel ! (m)
        logical, intent(in) :: rad
        real(dp), intent(in), dimension(:,:) :: Kdiff ! (m,m)
        real(dp), intent(in), dimension(:,:) :: zdiff, zdiff1, zdiff2 ! (m,m)
        real(dp), intent(inout), dimension(:) :: ap, at, ar ! (m)
        real(dp), intent(inout) :: ax(:,:) ! (m,m)
        real(dp), intent(inout) :: e_ap, e_at, e_ar, e_ax(:)
        integer, intent(out) :: ierr

        integer :: i,j,l,indx(n),nmax, rightshift, downshift
        real(dp), dimension(n) :: alpha, beta, nu, ga
        ! Add in beta later for rad lev. ga is the temp holder for the
        ! columns of gamm when doing matrix solve one column at a time.
        real(dp), dimension(n,n) :: delta, gamm

        ! initialize
        ierr = 0
        indx(1:n) = 0
        delta(1:n,1:n) = 0d0
        alpha(1:n) = 0d0
        beta(1:n) = 0d0
        nu(1:n) = 0d0
        ga(1:n) = 0d0
        gamm(1:n,1:n) = 0d0

        e_ap = 0d0
        e_at = 0d0
        e_ar = 0d0
        e_ax(1:m) = 0d0
        ! Just a way to make it clear what's going on with the thermal terms in delta
        ! by making block matrices and then shifting them into the proper position.
        ! This makes the comparisons for checking the subdiagonals easier, as well
        ! as indexing of the coefficients.
        rightshift = m 
        downshift = m-1 ! Because electron momentum equation dropped.

        ! Assign the RHS Matrix multiplying the unkown quantities.
        ! Assuming gravity is a known, so there are m uknown
        ! diffusion velocities, m unknown heat flow vectors, and
        ! one unkown electric field, corresponding to (m-1) momentum equations
        ! (no electron equation), m energy equations, and 2 conservation equations.

        ! Terms corresponding to the momentum equations.
        do i = 1,(m-1)

           ! Terms that multiply the diffusion velocities.
           do j = 1,m
              if (j == i) then
                 do l = 1,m
                    if (l /= j) then
                       delta(i,j) = delta(i,j) - Kdiff(i,l)
                    end if
                 end do
              else
                 delta(i,j) = Kdiff(i,j)
              end if
           end do

           ! Terms that multiply the heat flow vectors. 
           do j = 1,m
              if (j == i) then
                 do l = 1,m
                    if (l /= j) then
                       delta(i,j+rightshift) = delta(i,j+rightshift) + Kdiff(i,l)*zdiff(i,l)*A(l)/(A(i)+A(l))
                    end if
                 end do
              else
                 delta(i,j+rightshift) = -1d0*Kdiff(i,j)*zdiff(i,j)*A(i)/(A(i) + A(j))
              end if
           end do

           ! Term multiplying the electric field.
           delta(i,n) = nd(i)*Z(i)

        end do

        ! Terms corresponding to the energy equations.
        do i = 1,m ! All these entries get shifted lower into the i+downshift position.
           
           ! Terms that multiply the diffusion velocities.
           do j = 1,m
              if (j == i) then
                 do l = 1,m
                    if (l /= j) then
                       delta(i+downshift,j) = delta(i+downshift,j) + 2.5d0*Kdiff(i,l)*zdiff(i,l)*A(l)/(A(i) + A(l))
                    end if
                 end do
              else
                 delta(i+downshift,j) = -2.5d0*Kdiff(i,j)*zdiff(i,j)*A(j)/(A(i) + A(j))
              end if
           end do

           ! Terms that multiply the heat flow vectors.
           do j = 1,m
              ! Shift these entries to be in the lower-right half of the matrix.
              if (j == i) then
                 delta(i+downshift,j+rightshift) = -0.4d0*Kdiff(i,j)*zdiff2(i,j) ! 1st term on RHS of eqn (90) in MESA3
                 do l=1,m ! second line of eqn (90) in MESA3
                    if (l /= j) then
                       delta(i+downshift,j+rightshift) = delta(i+downshift,j+rightshift) - &
                            Kdiff(i,l)* &
                            ( &
                            (3d0*pow2(A(i)) + pow2(A(l))*zdiff1(i,l))/pow2(A(i)+A(l)) + &
                            0.8d0*zdiff2(i,l)*A(i)*A(l)/pow2(A(i) + A(l)) &
                            )
                    end if
                 end do
              else ! third line of eqn (90) in MESA3
                 delta(i+downshift,j+rightshift) = Kdiff(i,j)* &
                      ( 3d0 + zdiff1(i,j) - 0.8d0*zdiff2(i,j) )* &
                      A(i)*A(j)/pow2(A(i)+A(j))
              end if              
           end do
           
           ! Term multiplying the electric field. (doesn't appear in energy equations)
           delta(i+downshift,n) = 0d0
        end do

        ! Terms from the mass and charge conservation laws, the last two rows of the matrix.
        ! These only involve the diffusion velocities, not heat flow or e-field.
        do j = 1,m
           delta(n-1,j) = A(j)*nd(j)
           delta(n,j) = Z(j)*nd(j)
        end do
        ! Delta assignment finished.

        ! Assign LHS quantities. (Everything not assigned here stays zero thanks to earlier initialization)
        ! LHS of momentum equations.
        do i = 1,(m-1)
           alpha(i) = nd(i)*A(i)
           nu(i) = nd(i)
           gamm(i,i) = nd(i)
           if (rad) then
              beta(i) = -1d0*alpha(i)*amu*rad_accel(i)
           end if
        end do
        ! LHS of energy equations.
        do i = 1,m
           nu(i+downshift) = 2.5d0*nd(i)
        end do

        ! alpha,nu,gamm assignment finished


        ! Invert the matrix and solve to get the contribution of each
        ! vector/column on the LHS. This takes alpha, nu, ga from being
        ! the knowns on the LHS to vectors containing information on
        ! contributions for the unknowns. These then get assigned to the
        ! output arrays ap, at, ax, which will then be multiplied by
        ! gravity, k_b*T, and number density gradients to provide the
        ! full solution for diffusion velocities.

        call dgetrf(n,n,delta,n,indx,ierr)

        if( ierr /= 0 ) then
           ! print *, "Factoring failed!"
           return
        end if
        
        call dgetrs('N',n,1,delta,n,indx,alpha,n,ierr)
        if( ierr /= 0 ) then
           ! print *, "solve failed on alpha"
           return
        end if

        call dgetrs('N',n,1,delta,n,indx,nu,n,ierr)
        if( ierr /= 0 ) then
           ! print *, "solve failed on nu"
           return
        end if

        if (rad) then
           call dgetrs('N',n,1,delta,n,indx,beta,n,ierr)
           if( ierr /= 0 ) then
              ! print *, "solve failed on beta"
              return
           end if
        else
           beta(1:n) = 0d0
        end if

        do j=1,n
           do i=1,n
              ga(i) = gamm(i,j)
           end do
           call dgetrs('N',n,1,delta,n,indx,ga,n,ierr)
           if( ierr /= 0 ) then
              ! print *, "solve failed on gamma"
              return
           end if
           do i=1,n
              gamm(i,j) = ga(i)
           end do
        end do
        

        ! Assign the results of the matrix solve from the diffusion
        ! velocity part of the arrays to the solution vectors.
        
        do j = 1,m
           ap(j) = alpha(j)
           at(j) = nu(j)
           ar(j) = beta(j)
           do i=1,m
              ax(i,j) = gamm(i,j)
           end do
        end do

        ! Entries m+1-2m are heat flow vectors, so those don't get output.
        ! Skip to final entry (n) for electric field.
        e_ap = alpha(n)
        e_at = nu(n)
        e_ar = beta(n)
        do i=1,m
           e_ax(i) = gamm(n,i)
        end do
        
      end subroutine solve_burgers_cgs_with_thermal


      
      ! Calculate coefficients given in Appendix C.3 of Stanton & Murillo, PR E 93, 043203 (2016)
      subroutine get_SM_coeffs(nc,m,rho,T,A,Z,nd,Kdiff,zdiff,zdiff1,zdiff2,kappa)
        integer, intent(in) :: nc, m
        real(dp), intent(in) :: rho, T
        real(dp), dimension(:), intent(in) :: A, Z, nd ! m
        real(dp), intent(out) :: kappa ! ion separation relative to electron screening length.
        real(dp), dimension(:,:), intent(inout) :: Kdiff,zdiff,zdiff1,zdiff2 ! (m,m) but electron entries shouldn't be used.
        ! Only the ion terms are modified by this routine right now.

        real(dp), dimension(nc,nc) :: g_plasma, mu
        real(dp), dimension(nc,nc,2,3) :: Keff, Omega ! Indexed by all species and different orders
        ! Fitting coefficients from Stanton & Murillo
        real(dp), dimension(2,3) :: a1,a2,a3,a4,a5,b0,b1,b2,b3,b4
        real(dp) :: lambda ! The screening length
        integer :: i,j,k,facmo,no,mo ! Last two are used for different orders of collisions.

        real(dp) :: lgp, gp1, gp2, gp3, gp4, gp5, kbT32, tmp

        ! Initialize all to 0, since some entries never get set or used.
        a1(:,:) = 0d0
        a2(:,:) = 0d0
        a3(:,:) = 0d0
        a4(:,:) = 0d0
        a5(:,:) = 0d0
        b0(:,:) = 0d0
        b1(:,:) = 0d0
        b2(:,:) = 0d0
        b3(:,:) = 0d0
        b4(:,:) = 0d0

        ! Table IV, coefficients for fits (C23)-(C24)
        a1(1,1) = 1.4660d0
        a1(1,2) = 0.52094d0
        a1(1,3) = 0.30346d0
        a1(2,2) = 0.85401d0

        a2(1,1) = -1.7836d0
        a2(1,2) = 0.25153d0
        a2(1,3) = 0.23739d0
        a2(2,2) = -0.22898d0

        a3(1,1) = 1.4313d0
        a3(1,2) = -1.1337d0
        a3(1,3) = -0.62167d0
        a3(2,2) = -0.60059d0

        a4(1,1) = -0.55833d0
        a4(1,2) = 1.2155d0
        a4(1,3) = 0.56110d0
        a4(2,2) = 0.80591d0

        a5(1,1) = 0.061162d0
        a5(1,2) = -0.43784d0
        a5(1,3) = -0.18046d0
        a5(2,2) = -0.30555d0

        b0(1,1) = 0.081033d0
        b0(1,2) = 0.20572d0
        b0(1,3) = 0.68375d0
        b0(2,2) = 0.43475d0

        b1(1,1) = -0.091336d0
        b1(1,2) = -0.16536d0
        b1(1,3) = -0.38459d0
        b1(2,2) = -0.21147d0

        b2(1,1) = 0.051760d0
        b2(1,2) = 0.061572d0
        b2(1,3) = 0.10711d0
        b2(2,2) = 0.11116d0

        b3(1,1) = -0.50026d0
        b3(1,2) = -0.12770d0
        b3(1,3) = 0.10649d0
        b3(2,2) = 0.19665d0

        b4(1,1) = 0.17044d0
        b4(1,2) = 0.066993d0
        b4(1,3) = 0.028760d0
        b4(2,2) = 0.15195d0

        ! Get the screening length
        call lam_SM(nc,m,rho,T,Z,nd,lambda,kappa)

        ! Calculate g_plasma
        do i=1,nc
           do j=1,nc
              g_plasma(i,j) = Z(i)*Z(j)*qe*qe/(lambda*boltzm*T)
           end do
        end do

        ! Calculate Keff using Eqns (C22-C24)
        do i=1,nc
           do j=1,nc
              ! cache g powers
              lgp = safe_log(g_plasma(i,j))
              gp1 = g_plasma(i,j)
              gp2 = pow2(gp1)
              gp3 = pow3(gp1)
              gp4 = pow4(gp1)
              gp5 = pow5(gp1)
              if( g_plasma(i,j) < 0d0) then ! Don't calculate for attractive potentials, set to 0
                 Keff(i,j,:,:) = 0d0
              else if( g_plasma(i,j) < 1d0) then ! Use eqn C23 for weakly coupled
                 do no=1,2
                    do mo=1,3 ! Implementing the (m-1)! term with a simple if statement.
                       if(mo .eq. 3) then
                          facmo = 2 ! (3-1)!
                       else
                          facmo = 1 ! (1-1)! and (2-1)!
                       end if
                       Keff(i,j,no,mo) = (-1d0*no/4d0)*facmo*safe_log( &
                            a1(no,mo)*gp1 &
                            + a2(no,mo)*gp2 &
                            + a3(no,mo)*gp3 &
                            + a4(no,mo)*gp4 &
                            + a5(no,mo)*gp5)
                    end do
                 end do
              else ! Use eqn C24 for strongly coupled
                 do no=1,2
                    do mo=1,3
                       Keff(i,j,no,mo) = &
                            (b0(no,mo) + b1(no,mo)*lgp + b2(no,mo)*pow2(lgp) ) / &
                            (1d0 + b3(no,mo)*gp1 + b4(no,mo)*gp2 )
                    end do
                 end do
              end if
           end do
        end do

        ! Calculate the collision integrals using (C19)
        kBT32 = pow(boltzm*T,1.5d0)
        do i=1,nc
           do j=1,nc
              mu(i,j) = amu*A(i)*A(j)/(A(i) + A(j)) ! Reduced mass for collision
              tmp = sqrt(2d0*pi/(mu(i,j))) * (pow2(Z(i)*Z(j)*qe*qe)/kBT32)
              do no=1,2
                 do mo=1,3
                    Omega(i,j,no,mo) = tmp * Keff(i,j,no,mo)
                 end do
              end do
           end do
        end do

        ! Resistance coefficient is (16/3) ni nj mu Omega11
        ! (compare Paquette eqn 22 with SM eq B8)
        ! The zdiffs are given in terms of collision integrals in Paquette eqns (14-16), (23-25)
        do i=1,nc
           do j=1,nc
              Kdiff(i,j) = (16d0/3d0)*nd(i)*nd(j)*mu(i,j)*Omega(i,j,1,1)
              zdiff(i,j) = 1d0 - 2d0*Omega(i,j,1,2)/(5d0*Omega(i,j,1,1))
              zdiff1(i,j) = 2.5d0 + (2d0*Omega(i,j,1,3) - 10d0*Omega(i,j,1,2))/(5d0*Omega(i,j,1,1))
              zdiff2(i,j) = Omega(i,j,2,2)/Omega(i,j,1,1)
           end do
        end do

        ! Note that SM don't give fits for attractive potentials, so we haven't
        ! touched the electron entries. They exit this routine unchanged, so they
        ! either need to be initialized before this routine is called or somehow
        ! calculated later.
        
      end subroutine get_SM_coeffs

      ! Screening Length according to Stanton & Murillo
      subroutine lam_SM(nc,m,rho,T,Z,nd,lam_eff,kappa)
        integer, intent(in) :: nc, m
        real(dp), intent(in) :: rho, T
        real(dp), dimension(:), intent(in) :: Z, nd ! m. charges, number densities of all species
        real(dp), intent(out) :: lam_eff
        real(dp), intent(out) :: kappa ! ion separation relative to electron screening length.

        real(dp) :: tiny_n, ne, EF, lam_e, lam_sum, rhotot, ni_tot
        real(dp), dimension(nc) :: ai, gam_is, lam_i ! ion sphere radius, coupling, screening for each type of ion
        integer :: i

        tiny_n = 1d-20 ! g/cc
        ne = nd(m)
        ! Electron Fermi energy
        EF = (hbar*hbar*pow(3d0*pi*pi*ne,2d0/3d0))/(2d0*me)
        ! Electron screening length accounting for degeneracy correction
        lam_e = pow(pi4*qe*qe*ne/sqrt(pow2(boltzm*T) + pow2((2d0/3d0)*EF)),-0.5d0)

        ! Compute kappa
        ni_tot = 0d0
        do i = 1,nc
           ni_tot = ni_tot + nd(i)
        end do
        kappa = pow(3d0/(pi4*ni_tot),one_third)/lam_e

        rhotot = 0d0
        do i = 1,nc
           ! rhotot is a CHARGE density. Just trying to follow the notation of SM eq (34)
           rhotot = rhotot + Z(i)*qe*nd(i) ! = qe*ne?
        end do
        do i = 1,nc
           ai(i) = pow(Z(i)*qe/(four_thirds_pi*rhotot), one_third)
           gam_is(i) = Z(i)*Z(i)*qe*qe/(ai(i)*boltzm*T)
           ! Number densities that are 0 or tiny cause screening length to diverge.
           ! This is physical; nothing is there to screen anything.
           ! But numerically, I don't want to rely on fortran to handle division by zero and infinity,
           ! so just set these screening lengths to be huge and they won't
           ! contribute anything to overall screening length.
           if(nd(i) < tiny_n) then
              lam_i(i) = 1d99 ! cm. This won't contribute to any screening.
           else
              lam_i(i) = sqrt(boltzm*T/( pi4*pow2(Z(i)*qe)*nd(i) ))
           end if
        end do

        lam_sum = 1d0/pow2(lam_e) ! The electron part of the screening length.
        do i = 1,nc ! Sum over all the ions.
           lam_sum = lam_sum + pow( pow2(lam_i(i))*(1d0+3d0*gam_is(i)), -1d0)
        end do
        lam_eff = pow(lam_sum,-0.5d0)
      end subroutine lam_SM


      end module diffusion_support

