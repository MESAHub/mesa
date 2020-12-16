! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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


      module mlt_get_results

      use star_private_def
      use const_def
      use num_lib
      use utils_lib

      implicit none

      private
      public :: Get_results

      logical, parameter :: dbg = .false.
      integer, parameter :: kdbg = -1
      
      integer, parameter :: nvbs = num_mlt_partials

      contains
      

      subroutine Get_results(ss, kz, &
            cgrav, m, mstar, r, L, xh, &            
            T, rho, P, chiRho, chiT, Cp, opacity, grada, &            
            a_00, a_m1, &
            T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
            chiRho_for_partials_00, chiT_for_partials_00, &
            chiRho_for_partials_m1, chiT_for_partials_m1, &
            chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
            chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
            chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
            chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
            Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
            Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
            opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
            opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
            grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
            grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &            
            gradr_factor, d_gradr_factor_dw, gradL_composition_term, &
            alpha_semiconvection, semiconvection_option, &
            thermohaline_coeff, thermohaline_option, &
            dominant_iso_for_thermohaline, &
            mixing_length_alpha, alt_scale_height, remove_small_D_limit, &
            MLT_option, Henyey_y_param, Henyey_nu_param, &
            normal_mlt_gradT_factor, &
            prev_conv_vel, max_conv_vel, g_theta, dt, tau, just_gradr, mixing_type, &
            gradT, d_gradT_dvb, &
            gradr, d_gradr_dvb, &
            gradL, d_gradL_dvb, &
            scale_height, d_scale_height_dvb, &
            Lambda, d_Lambda_dvb, &
            conv_vel, d_conv_vel_dvb, & ! convection velocity
            D, d_D_dvb, &
            D_semi, d_D_semi_dvb, &
            D_thrm, d_D_thrm_dvb, &
            Gamma, d_Gamma_dvb, &
            ierr)

         use utils_lib, only: is_bad
         use chem_def, only: chem_isos
         
         type (star_info), pointer :: ss
         integer, intent(in) :: kz
         real(dp), intent(in) :: &
            cgrav, m, mstar, r, L, xh, &            
            T, rho, P, chiRho, chiT, Cp, opacity, grada, &            
            a_00, a_m1, &
            T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
            chiRho_for_partials_00, chiT_for_partials_00, &
            chiRho_for_partials_m1, chiT_for_partials_m1, &
            chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
            chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
            chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
            chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
            Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
            Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
            opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
            opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
            grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
            grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &
            gradr_factor, d_gradr_factor_dw, gradL_composition_term, &
            alpha_semiconvection, thermohaline_coeff, mixing_length_alpha, &
            Henyey_y_param, Henyey_nu_param, &
            prev_conv_vel, max_conv_vel, g_theta, dt, tau, remove_small_D_limit, &
            normal_mlt_gradT_factor
            
         logical, intent(in) :: alt_scale_height
         character (len=*), intent(in) :: thermohaline_option, MLT_option, semiconvection_option
         integer, intent(in) :: dominant_iso_for_thermohaline
         logical, intent(in) :: just_gradr
         
         integer, intent(out) :: mixing_type
         real(dp), intent(inout) :: gradT, d_gradT_dvb(:)
         real(dp), intent(inout) :: gradr, d_gradr_dvb(:)
         real(dp), intent(inout) :: gradL, d_gradL_dvb(:)
         real(dp), intent(inout) :: scale_height, d_scale_height_dvb(:)
         real(dp), intent(inout) :: Lambda, d_Lambda_dvb(:)
         real(dp), intent(inout) :: conv_vel, d_conv_vel_dvb(:)
         real(dp), intent(inout) :: D, d_D_dvb(:), D_semi, d_D_semi_dvb(:), D_thrm, d_D_thrm_dvb(:)
         real(dp), intent(inout) :: Gamma, d_Gamma_dvb(:) ! convective efficiency
         
         integer, intent(out) :: ierr

         real(dp) :: scale_height1, scale_height2
         real(dp) :: Pg, Pr, dP_dvb(nvbs), dPg_dvb(nvbs), dPr_dvb(nvbs), dRho_dvb(nvbs)
         real(dp) :: dT_dvb(nvbs), dL_dvb(nvbs), alpha, phi, dgrad, denom, tmp
         
         real(dp) :: grav, d_grav_dvb(nvbs)
         real(dp) :: diff_grads, d_diff_grads_dvb(nvbs)
         real(dp) :: convective_conductivity, d_cc_dvb(nvbs)
         real(dp) :: radiative_conductivity, d_rc_dvb(nvbs)
         real(dp) :: surf, dsurf_dvb(nvbs)
         real(dp) :: beta, d_beta_dvb(nvbs)
         real(dp) :: chi, d_chi_dvb(nvbs)
         real(dp) :: D_div_B, d_D_div_B_dvb(nvbs)
         real(dp) :: Q, dQ_dvb(nvbs)
         real(dp) :: A, dA_dvb(nvbs)
         real(dp) :: Bcubed, d_Bcubed_dvb(nvbs)
         real(dp) :: Zeta, d_Zeta_dvb(nvbs)
         real(dp) :: d_Cp_dvb(nvbs)
         real(dp) :: dR_dvb(nvbs)
         real(dp) :: d_opacity_dvb(nvbs)
         real(dp) :: d_grada_dvb(nvbs)
         real(dp) :: Dconv, d_Dconv_dvb(nvbs)
         real(dp) :: delta, d_delta_dvb(nvbs)
         real(dp) :: f, f0, d_f0_dvb(nvbs)
         real(dp) :: f1, d_f1_dvb(nvbs)
         real(dp) :: f2, d_f2_dvb(nvbs)
         real(dp) :: x, d_x_dvb(nvbs)

         real(dp) :: d_chiT_dvb(nvbs), d_chiRho_dvb(nvbs)

         real(dp) :: a0, omega, theta, s_gradr, f_gradr, dilute_factor
         real(dp) :: d_omega_dvb(nvbs), d_a0_dvb(nvbs), d_theta_dvb(nvbs)
         
         integer :: i
         real(dp), parameter :: tiny = 1d-30, min_D_th = 1d-3
         character (len=256) :: message        
         logical ::  quit
         real(dp) :: diff_grad, K, gamma0, L_ratio, frac, s, &
            dilution_factor, conv_tau, init_conv_vel
         real(dp) :: K_T, K_mu, nu_rad, nu_mol, nu, grad_mu, R0, r_th, H_P

         real(dp) :: scale_factor, interp_factor, dinterp_factor
         real(dp) :: gradT_temp, d_gradT_temp_dvb(nvbs)
         
         logical :: test_partials, debug

         ! These variables are to scale gradr
         real(dp) :: Lrad_div_Ledd, alfa0, diff_grads_factor, Gamma_factor, grad_scale, &
            Gamma_inv_threshold, Gamma_term, Gamma_limit, scale_value1, scale_value2,diff_grads_limit, &
            reduction_limit, lambda_limit, exp_limit
         real(dp) :: d_Lrad_div_Ledd_dvb(nvbs), d_alfa0_dvb(nvbs), d_diff_grads_factor_dvb(nvbs), &
            d_Gamma_factor_dvb(nvbs), d_grad_scale_dvb(nvbs), d_Gamma_inv_threshold_dvb(nvbs), &
            d_Gamma_term_dvb(nvbs)

         include 'formats.dek'

         !test_partials = (kz == ss% solver_test_partials_k)
         test_partials = .false.
      
         ierr = 0
         gradT_temp = 0
         debug = .false.
         
         dL_dvb(:) = 0d0
         dL_dvb(mlt_dL) = 1d0
         
         dR_dvb(:) = 0d0
         dR_dvb(mlt_dlnR) = r
         
         call set1_dvb(dT_dvb, &
            0d0, a_00*T_00, 0d0, a_m1*T_m1)
            
         call set1_dvb(dRho_dvb, &
            a_00*rho_00, 0d0, a_m1*rho_m1, 0d0)
            
         call set1_dvb(dP_dvb, &
            a_00*P_00*chiRho_for_partials_00, a_00*P_00*chiT_for_partials_00, &
            a_m1*P_m1*chiRho_for_partials_m1, a_m1*P_m1*chiT_for_partials_m1)
            
         call set1_dvb(d_chiT_dvb, &
            a_00*d_chiT_00_dlnd, a_00*d_chiT_00_dlnT, &
            a_m1*d_chiT_m1_dlnd, a_m1*d_chiT_m1_dlnT)
            
         call set1_dvb(d_chiRho_dvb, &
            a_00*d_chiRho_00_dlnd, a_00*d_chiRho_00_dlnT, &
            a_m1*d_chiRho_m1_dlnd, a_m1*d_chiRho_m1_dlnT)
            
         call set1_dvb(d_Cp_dvb, &
            a_00*d_Cp_00_dlnd, a_00*d_Cp_00_dlnT, &
            a_m1*d_Cp_m1_dlnd, a_m1*d_Cp_m1_dlnT)
            
         call set1_dvb(d_opacity_dvb, &
            a_00*d_opacity_00_dlnd, a_00*d_opacity_00_dlnT, &
            a_m1*d_opacity_m1_dlnd, a_m1*d_opacity_m1_dlnT)
            
         call set1_dvb(d_grada_dvb, &
            a_00*d_grada_00_dlnd, a_00*d_grada_00_dlnT, &
            a_m1*d_grada_m1_dlnd, a_m1*d_grada_m1_dlnT)

         if (is_bad(d_grada_dvb(mlt_dlnd00))) then
            ierr = -1
            if (.not. ss% report_ierr) return
!$OMP critical (mlt_info_crit5)
            write(*,2) 'd_grada_dvb(mlt_dlnd00)', kz, d_grada_dvb(mlt_dlnd00)
            write(*,2) 'a_00', kz, a_00
            write(*,2) 'a_m1', kz, a_m1
            write(*,2) 'grada', kz, grada
            write(*,2) 'd_grada_00_dlnd', kz, d_grada_00_dlnd
            write(*,2) 'd_grada_00_dlnT', kz, d_grada_00_dlnT
            write(*,2) 'd_grada_m1_dlnd', kz, d_grada_m1_dlnd
            write(*,2) 'd_grada_m1_dlnT', kz, d_grada_m1_dlnT
            call mesa_error(__FILE__,__LINE__)
!$OMP end critical (mlt_info_crit5)
            return
         end if

         Pr = one_third*crad*T*T*T*T
         if (debug) write(*,1) 'Pr', Pr
         call set1_dvb(dPr_dvb, &
            0d0, 4d0*Pr*a_00*T_00/T, &
            0d0, 4d0*Pr*a_m1*T_m1/T)
         
         !gradr = eval_Paczynski_gradr(P,opacity,L,m,cgrav,Pr,tau,T,r,rho)
         gradr = P*opacity*L / (16*pi*clight*m*cgrav*Pr)
         if (tau < 2d0/3d0) then ! B. Paczynski, 1969, Acta Astr., vol. 19, 1., eqn 14.
            s_gradr = (2d0*crad*T*T*T*sqrt(r))/(3d0*cgrav*m*rho)*pow(L/(8d0*pi*boltz_sigma), 0.25d0) ! eqn 15
            f_gradr = 1d0 - 1.5d0*tau ! Paczynski, 1969, eqn 8
            dilute_factor = (1 + f_gradr*s_gradr*(4*pi*cgrav*clight*m)/(opacity*L))/(1 + f_gradr*s_gradr)
            gradr = gradr*dilute_factor
         end if
         gradr = gradr*gradr_factor
         d_gradr_dvb = gradr*(dP_dvb/P + d_opacity_dvb/opacity - dPr_dvb/Pr)
         d_gradr_dvb(mlt_dL) = gradr_factor*P*opacity / (16*pi*clight*m*cgrav*Pr)
         if (ss% w_div_wc_flag) then
            d_gradr_dvb(mlt_w_div_wc_var) = d_gradr_factor_dw*P*opacity*L / (16*pi*clight*m*cgrav*Pr)
         end if
         
         if (test_partials) then
            ss% solver_test_partials_val = opacity
            ss% solver_test_partials_var = ss% i_lnd
            ss% solver_test_partials_dval_dx = d_opacity_dvb(mlt_dlnd00)
            write(*,*) 'mlt get_results   opacity', ss% solver_test_partials_var
         end if
         
         
         
         if (is_bad(gradr)) then
            ierr = -1
            if (.not. ss% report_ierr) return
!$OMP critical (mlt_info_crit6)
            write(*,2) 'gradr', kz, gradr
            write(*,2) 'P', kz, P
            write(*,2) 'L', kz, L
            write(*,2) 'opacity', kz, opacity
            write(*,2) 'm', kz, m
            write(*,2) 'cgrav', kz, cgrav
            write(*,2) 'standard_cgrav', kz, standard_cgrav
            write(*,2) 'Pr', kz, Pr
            write(*,2) '16*pi*clight*m*cgrav*Pr', kz, 16*pi*clight*m*cgrav*Pr
            write(*,2) 'P*opacity*L', kz, P*opacity*L
            write(*,2) 'gradr_factor', kz, gradr_factor
            write(*,2) 'tau', kz, tau
            !write(*,2) '', kz, 
            call mesa_error(__FILE__,__LINE__)
!$OMP end critical (mlt_info_crit6)
         end if
         
         if (just_gradr) return

         grav = cgrav*m / (r*r)
         d_grav_dvb = 0d0
         d_grav_dvb(mlt_dlnR) = -2*grav
         
         if (grav < 0) then
            write(*,1) 'grav', grav
            write(*,1) 'cgrav', cgrav
            write(*,1) 'm', m
            write(*,1) 'r', r
            call mesa_error(__FILE__,__LINE__)
         end if

         scale_height = P / (grav*rho)
         d_scale_height_dvb = scale_height*(dP_dvb/P - d_grav_dvb/grav - dRho_dvb/Rho)
         if (alt_scale_height) then
            ! consider sound speed*hydro time scale as an alternative scale height
            ! (this comes from Eggleton's code.)
            scale_height2 = sqrt(P/cgrav)/rho
            if (scale_height2 < scale_height) then
               scale_height = scale_height2
               d_scale_height_dvb = scale_height*(0.5d0*dP_dvb/P - dRho_dvb/Rho)
            end if
         end if
         H_P = scale_height

         if (scale_height <= 0d0 .or. is_bad(scale_height)) then
            ierr = -1
            return
!$OMP critical (mlt_info_crit7)
            write(*,1) 'scale_height', scale_height
            stop 'set_convective_mixing'
!$OMP end critical (mlt_info_crit7)
         end if

         if (is_bad(d_scale_height_dvb(mlt_dlnd00))) then
            ierr = -1
            return
!$OMP critical (mlt_info_crit8)
            write(*,1) 'd_scale_height_dvb(mlt_dlnd00)', d_scale_height_dvb(mlt_dlnd00)
            stop 'set_convective_mixing'
!$OMP end critical (mlt_info_crit8)
            return
         end if
         
         surf = pi4*r*r
         if (debug) write(*,1) 'surf', surf
         dsurf_dvb = 8*pi*r*dR_dvb
         
         ! Ledoux temperature gradient (same as Schwarzschild if composition term = 0)
         gradL = grada + gradL_composition_term
         d_gradL_dvb = d_grada_dvb ! ignore partials of composition term
         
         diff_grads = gradr - gradL ! convective if this is > 0
         d_diff_grads_dvb = d_gradr_dvb - d_gradL_dvb
         if (is_bad(d_diff_grads_dvb(mlt_dlnT00))) then
            ierr = -1
            return
!$omp critical (mlt_info_crit9)
            write(*,1) 'd_grada_dvb(mlt_dlnT00)', d_grada_dvb(mlt_dlnT00)
            write(*,1) 'd_gradr_dvb(mlt_dlnT00)', d_gradr_dvb(mlt_dlnT00)
            write(*,1) 'd_gradL_dvb(mlt_dlnT00)', d_gradL_dvb(mlt_dlnT00)
            write(*,1) 'd_diff_grads_dvb(mlt_dlnT00)', d_diff_grads_dvb(mlt_dlnT00)
            call mesa_error(__FILE__,__LINE__)
!$omp end critical (mlt_info_crit9)
         end if

         Pg = P - Pr
         if (debug) write(*,1) 'Pg', Pg
         if (Pg < tiny) then
            call set_no_mixing
            return
         end if
         
         dPg_dvb = dP_dvb - dPr_dvb

         beta = Pg / P
         if (debug) write(*,1) 'beta', beta
         d_beta_dvb = beta*(dPg_dvb/Pg - dP_dvb/P)
         
         if (debug) write(*,1) 'scale_height', scale_height

         ! mixing length, Lambda
         Lambda = mixing_length_alpha*scale_height
         if (debug) write(*,1) 'Lambda', Lambda
         d_Lambda_dvb = mixing_length_alpha*d_scale_height_dvb
                  
         if (mixing_length_alpha <= 0) then
            call set_no_mixing
            return
         end if
         
         if (MLT_option == 'none') then
            call set_no_mixing
            return
         end if
         
         if (opacity < 1d-10 .or. P < 1d-20 .or. T < 1d-10 .or. Rho < 1d-20 &
               .or. m < 1d-10 .or. r < 1d-10 .or. cgrav < 1d-10 .or. &
               max_conv_vel == 0d0) then
            if (.false.) then
               write(*,2) 'special set no mixing', kz
               write(*,*) 'opacity < 1d-10', opacity < 1d-10
               write(*,*) 'P < 1d-20', P < 1d-20
               write(*,*) 'T < 1d-10', T < 1d-10
               write(*,*) 'Rho < 1d-20', Rho < 1d-20
               write(*,*) 'm < 1d-10', m < 1d-10
               write(*,*) 'r < 1d-10', r < 1d-10
               write(*,*) 'cgrav < 1d-10', cgrav < 1d-10
               write(*,*) 'max_conv_vel == 0d0', max_conv_vel == 0d0       
               write(*,*) "MLT_option == 'none' ", MLT_option == 'none'      
               call mesa_error(__FILE__,__LINE__)
            end if
            call set_no_mixing
            return
         end if

         ! 'Q' param  C&G 14.24
         Q = chiT/chiRho
         dQ_dvb = Q*( d_chiT_dvb/chiT - d_chiRho_dvb/chiRho )
         if (Q <= 0) then
            call set_no_mixing
            return
         end if
                     
         radiative_conductivity = (4*crad*clight/3)*T*T*T / (opacity*rho) ! erg / (K cm sec)
         if (debug) write(*,1) 'radiative_conductivity', radiative_conductivity
         d_rc_dvb = radiative_conductivity*(3d0*dT_dvb/T - dRho_dvb/rho - d_opacity_dvb/opacity)
         
         if (diff_grads <= 0d0) then ! not convective (Ledoux stable)    
            call set_no_mixing ! also sets gradT = gradr    
            if (gradL_composition_term < 0) then ! composition unstable
               call set_thermohaline
               D = D_thrm
               d_D_dvb = d_D_thrm_dvb
               if (debug) write(*,1) 'after set_thermohaline D_thrm', D_thrm
               if (ss% conv_vel_flag .and. ss% conv_vel_ignore_thermohaline) then
                  conv_vel = 0d0
                  d_conv_vel_dvb = 0d0
                  !mixing_type = no_mixing
               end if
            else if (gradr > grada) then ! Schw unstable
               call set_semiconvection
               D = D_semi
               d_D_dvb = d_D_semi_dvb
               if (debug) write(*,1) 'after set_semiconvection D_semi', D_semi
               if (ss% conv_vel_flag .and. ss% conv_vel_ignore_semiconvection) then
                  conv_vel = 0d0
                  d_conv_vel_dvb = 0d0
                  !mixing_type = no_mixing
               end if
            else
               call time_limit_conv_vel
            end if
            if (debug) write(*,1) 'remove_small_D_limit', remove_small_D_limit
            if (D < remove_small_D_limit .or. is_bad(D)) then
               call set_no_mixing
            end if
            if (debug) write(*,1) 'final D', D
            if (conv_vel > 0d0) then
               D = conv_vel*Lambda/3     ! diffusion coefficient [cm^2/sec]
               if (debug) write(*,1) 'D', D
               d_D_dvb = (d_conv_vel_dvb*Lambda + conv_vel*d_Lambda_dvb)/3
               if (mixing_type == no_mixing) then
                  mixing_type = leftover_convective_mixing
               end if
            end if
            if (ss% conv_vel_flag) then
               if (normal_mlt_gradT_factor > 0d0) then
                  gradT_temp = gradT
                  d_gradT_temp_dvb = d_gradT_dvb
               end if
               if (normal_mlt_gradT_factor < 1d0) then
                  call revise_using_cv_var_variable
                  if (ierr /= 0) return
               end if
               if (normal_mlt_gradT_factor > 0d0 .and. normal_mlt_gradT_factor < 1d0) then
                  scale_factor = normal_mlt_gradT_factor
                  interp_factor = 10d0*pow(scale_factor,3d0) &
                     -15*pow(scale_factor,4d0)+6*pow(scale_factor,5d0)
                  gradT = (1-interp_factor)*gradT + interp_factor*gradT_temp
                  d_gradT_dvb = (1-interp_factor)*d_gradT_dvb + interp_factor*d_gradT_temp_dvb
               end if
               if (kz > 0) then
                  if (ss% conv_vel(kz) > 0d0 .and. mixing_type == no_mixing) &
                     mixing_type = leftover_convective_mixing
               end if
            end if
            return            
         end if
         
         call set_convective_mixing
         if (quit) return

         ! Reduce gradr-grada
         if (ss% use_superad_reduction) then

            Gamma_limit = ss% superad_reduction_Gamma_limit
            scale_value1 = ss% superad_reduction_Gamma_limit_scale
            scale_value2 = ss% superad_reduction_Gamma_inv_scale
            diff_grads_limit = ss% superad_reduction_diff_grads_limit
            reduction_limit = ss% superad_reduction_limit
            Lrad_div_Ledd = 4d0*crad/3d0*pow4(T)/P*gradT
            d_Lrad_div_Ledd_dvb = Lrad_div_Ledd*(4*dT_dvb/T-dP_dvb/P+d_gradT_dvb/gradT)
            Gamma_inv_threshold = 4*(1-beta)/chiT
            d_Gamma_inv_threshold_dvb = Gamma_inv_threshold*(-d_beta_dvb/(1-beta)-d_chiT_dvb/chiT)

            Gamma_factor = 1d0
            d_Gamma_factor_dvb = 0d0
            if (gradT-grada > 0d0) then
               if (Lrad_div_Ledd > Gamma_limit .or. Lrad_div_Ledd > Gamma_inv_threshold) then
                  alfa0 = (gradT-grada)/diff_grads_limit
                  if (alfa0 < 1d0) then
                     diff_grads_factor = -alfa0*alfa0*alfa0*(-10d0 + alfa0*(15d0 - 6d0*alfa0))
                     d_diff_grads_factor_dvb = &
                        30d0*(alfa0 - 1d0)*(alfa0 - 1d0)*alfa0*alfa0*(d_gradT_dvb-d_grada_dvb)/diff_grads_limit
                  else
                     diff_grads_factor = 1d0
                     d_diff_grads_factor_dvb = 0d0
                  end if

                  Gamma_term = 0d0
                  d_Gamma_term_dvb = 0d0
                  if (Lrad_div_Ledd > Gamma_limit) then
                     Gamma_term = Gamma_term + scale_value1*pow2(Lrad_div_Ledd/Gamma_limit-1d0)
                     d_Gamma_term_dvb = d_Gamma_term_dvb + scale_value1*2d0*(Lrad_div_Ledd/Gamma_limit-1d0)*d_Lrad_div_Ledd_dvb/Gamma_limit
                  end if
                  if (Lrad_div_Ledd > Gamma_inv_threshold) then
                     Gamma_term = Gamma_term + scale_value2*pow2(Lrad_div_Ledd/Gamma_inv_threshold-1d0)
                     d_Gamma_term_dvb= d_Gamma_term_dvb + scale_value2*2d0*(Lrad_div_Ledd/Gamma_inv_threshold-1d0)*&
                        (d_Lrad_div_Ledd_dvb/Gamma_inv_threshold &
                           -Lrad_div_Ledd/pow2(Gamma_inv_threshold)*d_Gamma_inv_threshold_dvb)
                  end if
                  
                  if (Gamma_term > 0d0) then
                     Gamma_factor = Gamma_term/beta*diff_grads_factor
                     d_Gamma_factor_dvb = &
                        Gamma_factor*(d_Gamma_term_dvb/Gamma_term -d_beta_dvb/beta + d_diff_grads_factor_dvb/diff_grads_factor)
                     Gamma_factor = Gamma_factor + 1d0
                     if (reduction_limit > 1d0) then
                        lambda_limit = 2d0/(reduction_limit-1d0)
                        exp_limit = exp(-lambda_limit*(Gamma_factor-1d0))
                        d_Gamma_factor_dvb = 2*(reduction_limit-1)*lambda_limit*exp_limit/(1+exp_limit)**2*d_Gamma_factor_dvb
                        Gamma_factor = 2*(reduction_limit-1)*(1d0/(1+exp_limit)-0.5)+1d0
                     end if
                  end if
               end if
            end if 
            if (kz /= 0) ss% superad_reduction_factor(kz) = Gamma_factor
            if (Gamma_factor > 1d0) then
               grad_scale = (gradr-grada)/(Gamma_factor*gradr) + grada/gradr
               d_grad_scale_dvb = (gradr-grada)/(Gamma_factor*gradr)* &
                  ((d_gradr_dvb-d_grada_dvb)/(gradr-grada) - (d_Gamma_factor_dvb*gradr + Gamma_factor*d_gradr_dvb)/(Gamma_factor*gradr))
               d_grad_scale_dvb = d_grad_scale_dvb + d_grada_dvb/gradr - grada*d_gradr_dvb/pow2(gradr)
               d_gradr_dvb = grad_scale*d_gradr_dvb + d_grad_scale_dvb*gradr
               gradr = grad_scale*gradr
               diff_grads = gradr - grada
               d_diff_grads_dvb = d_gradr_dvb - d_grada_dvb
               call set_convective_mixing
               if (quit) return
            end if
         end if
         
         D = conv_vel*Lambda/3     ! diffusion coefficient [cm^2/sec]
         if (debug) write(*,1) 'D', D
         d_D_dvb = (d_conv_vel_dvb*Lambda + conv_vel*d_Lambda_dvb)/3

         mixing_type = convective_mixing
         
         if (debug .or. D < 0) then
            write(*,*) 'get_gradT: convective_mixing'
            write(*,1) 'D', D
            write(*,1) 'conv_vel', conv_vel
            write(*,1) 'Pg/P', Pg/P
            write(*,1) 'H_P', H_P
            write(*,1) 'scale_height', scale_height
            write(*,1) 'scale_height/H_P', scale_height/H_P
            write(*,1) 'm/Msun', m/Msun
            write(*,1) 'r/Rsun', r/Rsun
            write(*,1) 'T', T
            write(*,1) 'rho', rho
            write(*,1) 'grada', grada
            write(*,1) 'chiT', chiT
            write(*,1) 'chiRho', chiRho
            write(*,2) 'mixing_type', mixing_type
            write(*,*)
            if (.not. debug) stop 'MLT: get_gradT'
         end if

         if (D < remove_small_D_limit .or. is_bad(D)) then
            call set_no_mixing
         end if

         if (ss% conv_vel_flag) then
            if (normal_mlt_gradT_factor > 0d0) then
               gradT_temp = gradT
               d_gradT_temp_dvb = d_gradT_dvb
            end if
            if (normal_mlt_gradT_factor < 1d0) then
               call revise_using_cv_var_variable
               if (ierr /= 0) return
            end if
            if (normal_mlt_gradT_factor > 0d0 .and. normal_mlt_gradT_factor < 1d0) then
               scale_factor = normal_mlt_gradT_factor
               interp_factor = 10d0*pow(scale_factor,3d0) &
                  -15*pow(scale_factor,4d0)+6*pow(scale_factor,5d0)
               gradT = (1-interp_factor)*gradT + interp_factor*gradT_temp
               d_gradT_dvb = (1-interp_factor)*d_gradT_dvb + interp_factor*d_gradT_temp_dvb
            end if
            if (kz > 0) then ! e.g. pre_ms_model calls with kz = 0
               if (ss% conv_vel(kz) > 0d0 .and. mixing_type == no_mixing) then
                  mixing_type = leftover_convective_mixing
               end if
            end if
         end if
         
         contains

         
         subroutine set1_dvb(dvb, dlnd00, dlnT00, dlndm1, dlndTm1)
            real(dp), intent(inout) :: dvb(nvbs)
            real(dp), intent(in) :: dlnd00, dlnT00, dlndm1, dlndTm1
            dvb(mlt_dlnd00) = dlnd00
            dvb(mlt_dlnT00) = dlnT00
            dvb(mlt_dlndm1) = dlndm1
            dvb(mlt_dlnTm1) = dlndTm1
            dvb(mlt_dlnR) = 0d0
            dvb(mlt_dL) = 0d0
            dvb(mlt_cv_var) = 0d0
            dvb(mlt_w_div_wc_var) = 0d0
         end subroutine set1_dvb
         
         subroutine set_thermohaline
            real(dp) :: kipp_D, tau, Pr, BGS_C, Nu_mu, l2, lmda, phi, r_guess, &
            sqrt_1_plus_phi, sqrt_Pr
   
            logical, parameter :: dbg = .false.
            include 'formats'
                        
            if (dbg) write(*,*) 'set_thermohaline ' // trim(thermohaline_option)
            
            diff_grad = max(1d-40, grada - gradr) ! positive since Schwarzschild stable               
            K = 4*crad*clight*T*T*T/(3*opacity*rho) ! thermal conductivity
            
            if (thermohaline_option == 'Kippenhahn') then
            
               ! Kippenhahn, R., Ruschenplatt, G., & Thomas, H.-C. 1980, A&A, 91, 175
               D_thrm = -thermohaline_coeff*3*K/(2*rho*cp)*gradL_composition_term/diff_grad
            
            else if (thermohaline_option == 'Brown_Garaud_Stellmach_13' .or. &
                     thermohaline_option == 'Traxler_Garaud_Stellmach_11') then
                     
               call get_diff_coeffs(K_T,K_mu,nu)

               R0 = (gradr - grada)/gradL_composition_term
               Pr = nu/K_T
               tau = K_mu/K_T
               r_th = (R0 - 1d0)/(1d0/tau - 1d0)

               if (r_th >= 1d0) then ! stable if R0 >= 1/tau
                  D_thrm = 0d0
               else if (thermohaline_option == 'Traxler_Garaud_Stellmach_11') then 
                  ! Traxler, Garaud, & Stellmach, ApJ Letters, 728:L29 (2011).
                  ! also see Denissenkov. ApJ 723:563–579, 2010.
                  D_thrm = 101d0*sqrt(K_mu*nu)*exp(-3.6d0*r_th)*pow(1d0 - r_th,1.1d0) ! eqn 24
               else             
                  ! if (thermohaline_option == 'Brown_Garaud_Stellmach_13') then
                  D_thrm = K_mu*(Numu(R0,r_th,pr,tau) - 1d0)
                  ! evbauer 07/18: changed from K_mu*Numu(R0,r_th,pr,tau), Pascale signed off
               endif
               D_thrm = thermohaline_coeff*D_thrm
               
            else
                 
               D_thrm = 0
               ierr = -1
               write(*,*) 'unknown value for MLT thermohaline_option' // trim(thermohaline_option)
               return   
               
            end if
            
            conv_vel = 3*D_thrm/Lambda
            d_conv_vel_dvb = 0
            d_D_thrm_dvb = 0
            
            call time_limit_conv_vel
            if (conv_vel > max_conv_vel) conv_vel = max_conv_vel
            if (init_conv_vel /= conv_vel) D_thrm = conv_vel*Lambda/3
            
            if (D_thrm < min_D_th .or. D_thrm <= 0) then
               call set_no_mixing
               return
            end if
            
            mixing_type = thermohaline_mixing 

         end subroutine set_thermohaline
         
         subroutine time_limit_conv_vel
            real(dp), dimension(nvbs) :: d_vconv_accel_dvb, d_brunt_timescale_dvb
            real(dp), dimension(nvbs) :: d_diff_grads_s_dvb, d_tau_dvb
            real(dp) :: new_conv_vel, vconv_accel, l
            real(dp) :: brunt_timescale, diff_grads_s, tau
            include 'formats'

            init_conv_vel = conv_vel
            if (dt <= 0d0 .or. prev_conv_vel < 0d0) return

            if (conv_vel <= prev_conv_vel) return
            
            if (g_theta > 0) then ! max accel is grav*g_theta
               if (conv_vel > prev_conv_vel) then
                  !increase convective velocity as needed
                  new_conv_vel = prev_conv_vel + dt*grav*g_theta
                  if (new_conv_vel < conv_vel) then
                     conv_vel = new_conv_vel
                     d_conv_vel_dvb = dt*d_grav_dvb*g_theta
                  end if
               else
                  !reduce convective velocity as needed
                  new_conv_vel = prev_conv_vel - dt*grav*g_theta
                  if (new_conv_vel > conv_vel) then 
                     conv_vel = new_conv_vel
                     d_conv_vel_dvb = -dt*d_grav_dvb*g_theta
                  end if
               end if
            else
               ! Arnett, W.D., 1969, Ap. and Space Sci, 5, 180.
               l = Lambda
               if (conv_vel > 0d0) then
                  vconv_accel = 2d0*(conv_vel*conv_vel - prev_conv_vel*prev_conv_vel)/l
                  d_vconv_accel_dvb = 4d0*conv_vel*d_conv_vel_dvb/l &
                     -2d0*(conv_vel*conv_vel-prev_conv_vel*prev_conv_vel)*d_Lambda_dvb/l*l
                  if (conv_vel > prev_conv_vel) then
                     !increase convective velocity as needed
                     new_conv_vel = prev_conv_vel + dt*vconv_accel
                     if (new_conv_vel < conv_vel) then
                        conv_vel = new_conv_vel
                        d_conv_vel_dvb = dt*d_vconv_accel_dvb
                     end if
                  else
                     !reduce convective velocity as needed
                     new_conv_vel = prev_conv_vel + dt*vconv_accel
                     if (new_conv_vel > conv_vel) then
                        conv_vel = new_conv_vel
                        d_conv_vel_dvb = dt*d_vconv_accel_dvb
                     end if
                  end if
               else ! conv_vel == 0d0
                  if (prev_conv_vel == 0d0) return
                  diff_grads_s = diff_grads
                  d_diff_grads_s_dvb = d_diff_grads_dvb
                  if (diff_grads_s < 0d0) then
                     brunt_timescale = 1/sqrt(-chiT/chiRho*diff_grads_s*grav/scale_height)
                  else
                     brunt_timescale = 1d99
                  end if
                  if (l/prev_conv_vel > brunt_timescale) then
                     ! reduce conv_vel on the brunt_timescale
                     new_conv_vel = prev_conv_vel*exp(-dt/brunt_timescale)
                     d_brunt_timescale_dvb = 0.5d0*brunt_timescale*(&
                        -d_chiT_dvb/chiT + d_chiRho_dvb/chiRho - d_diff_grads_s_dvb/diff_grads_s &
                        -d_grav_dvb/grav + d_scale_height_dvb/scale_height)
                     d_conv_vel_dvb = d_brunt_timescale_dvb*&
                        (dt*exp(-dt/brunt_timescale)*prev_conv_vel)/pow(brunt_timescale,2)
                  else
                     new_conv_vel = prev_conv_vel*l/(2*prev_conv_vel*dt+l)
                     d_conv_vel_dvb = (2*dt*prev_conv_vel*prev_conv_vel)/&
                        (4*dt*dt*prev_conv_vel*prev_conv_vel+4*dt*l*prev_conv_vel+l*l) &
                        *d_Lambda_dvb
                  end if
                  conv_vel = new_conv_vel
               end if
            end if

            if (conv_vel > max_conv_vel) then
               conv_vel = max_conv_vel
               d_conv_vel_dvb = 0d0
            end if

         end subroutine time_limit_conv_vel

         subroutine get_diff_coeffs(kt,kmu,vis)

         use chem_def, only: chem_isos

         real(dp) :: kt,kmu,vis,qe4
         real(dp) :: loglambdah,loglambdacx,loglambdacy,ccx,ccy
         real(dp) :: Bcoeff
         real(dp) :: chemA,chemZ,acx,acy       
         real(dp), parameter :: sqrt5 = sqrt(5d0)
                 
         kt = K/(Cp*rho)       ! thermal diffusivity (assumes radiatively dominated)
         qe4=qe*qe*qe*qe
         
         ! Log Lambda for pure H (equation 10 from Proffitt Michaud 93)
         loglambdah = -19.26d0 - 0.5d0*log(rho) + 1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+xh)) 
         nu_rad = 4d0*crad*T*T*T*T/(15d0*clight*opacity*rho*rho) ! radiative viscosity
         nu_mol = 0.406d0*sqrt(amu)*pow(boltzm*T,2.5d0)/(qe4*loglambdah*rho) 
         ! From Spitzer "Physics of Fully Ionized Gases equation 5-54
         ! Assumes pure H. Still trying to work out what it would be for a mixture. 
         vis = nu_mol + nu_rad   ! total viscosity

         ! The following is from Proffitt & Michaud, 1993.
         ! Their constant B (equation 15)
         Bcoeff = (15.d0/16.d0)*sqrt(2.d0*amu/(5*pi))*pow(boltzm,2.5d0)/qe4
         ! Extract what species drives the thermohaline concvection
         chemA = chem_isos%Z_plus_N(dominant_iso_for_thermohaline)
         chemZ = chem_isos%Z(dominant_iso_for_thermohaline)

         if(chemZ.gt.2) then
         ! This is if the driving chemical is NOT He.
            ! Log Lambda for H-dominant chem mixture (equation 10)
            loglambdacx = loglambdah - log(chemz)  
            ! Log Lambda for He-dominant chem mixture (equation 10)
            loglambdacy = loglambdah - log(2.d0*chemz)
            ! Calculation of C_ij coeffs (equation 12)
            ccx = log(exp(1.2d0*loglambdacx)+1.)/1.2d0
            ccy = log(exp(1.2d0*loglambdacy)+1.)/1.2d0
            ! Reduced masses (I had to guess, from Bahcall & Loeb 1990), with H and He
            acx = (1.d0*chemA)/(1.d0+chemA)
            acy = 4*chemA/(4.d0+chemA)
            ! My formula (see notes) based on Proffitt and Michaud 1993
            kmu = 2*Bcoeff*pow(T,2.5d0)/(sqrt5*rho*chemZ*chemZ)/ &
               (xh*sqrt(acx)*ccx + (1-xh)*sqrt(acy)*ccy)

         else
            ! Log Lambda for H-He mixture (equation 10)
            loglambdah = -19.26d0 - log(2d0) - 0.5d0*log(rho) + &
               1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+xh)) 
            ! Calculation of C_ij coeffs (equation 12)
            ccy = log(exp(1.2d0*loglambdah)+1d0)/1.2d0
            ! My formula (see notes) based on Proffitt and Michaud 1993
            kmu = (Bcoeff*pow(T,2.5d0)/(rho*ccy))*(3+xh)/((1+xh)*(3+5*xh)*(0.7d0+0.3d0*xh))
            
         endif
         ! write(57,*) kt,kmu,vis,chemZ

         end subroutine get_diff_coeffs

         real(dp) function numu(R0,r_th,prandtl,diffratio)
         !Function calculates Nu_mu from input parameters, following Brown et al. 2013.
         !Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 

            real(dp), intent(in) :: R0,r_th,prandtl,diffratio
            real(dp) :: maxl2,maxl,lambdamax
            real(dp) :: myvars(2)
            integer :: ierr, iter

            ! Initialize guess using estimates from Brown et al. 2013
            call analytical_estimate_th(maxl,lambdamax,r_th,prandtl,diffratio)
                  
            myvars(1) = maxl
            myvars(2) = lambdamax

           !Call Newton relaxation algorithm
           call NR(myvars,prandtl,diffratio,R0,ierr)
          
           !If the growth rate is negative, then try another set of parameters as first guess.  
           !Repeat as many times as necessary until convergence is obtained.
           iter = 1
           do while((myvars(2)<0).or.(ierr /= 0)) 
              !write(*,*) 'Alternative', r_th,prandtl,diffratio,iter
           !Reset guess values
              myvars(1) = maxl
              myvars(2) = lambdamax
           !Call relaxation for slightly different Pr, tau, R0.
              call NR(myvars,prandtl*(1d0+iter*1.d-2),diffratio,R0/(1d0+iter*1.d-2),ierr)
           !If it converged this time, call NR for the real parameters.
              if(ierr.eq.0) call NR(myvars,prandtl,diffratio,R0,ierr)
              !write(*,*) prandtl,diffratio,R0,myvars(1),myvars(2),ierr
              !Otherwise, increase counter and try again.
              iter = iter + 1            
           enddo

           !Plug solution into "l^2" and lambda.
           maxl2 = myvars(1)*myvars(1)
           lambdamax = myvars(2) 
           !write(*,*) prandtl,diffratio,r_th,maxl2,lambdamax

           !Calculate Nu_mu using Formula (33) from Brown et al, with C = 7.
           numu = 1.d0 + 49.d0*lambdamax*lambdamax/(diffratio*maxl2*(lambdamax+diffratio*maxl2))

         return
         end function numu 

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
            maxl = pow((1d0/3d0)*(1d0-r_th) + (1d0-r_th)*(1d0-r_th)*(5d0-4d0*phi)/27d0,0.25d0)
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
         err = 0d0
         errold = 0d0
         !Save input guess (probably not necessary here but useful in other routines)
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
            
         !if(iter.eq.niter) write(*,2) 'Failed to converge'
         return
         end subroutine NR
         
         subroutine set_convective_mixing
            ! need to set gradT, d_gradT_dvb, conv_vel, d_conv_vel_dvb
            include 'formats.dek'
            real(dp) ff1, ff2, ff3, ff4, ff5, aa, bb, y0, xres, a1, a2, sqrt_x
            real(dp) :: A_0, A_1, A_2, A_numerator, A_denom, inv_sqrt_x
            real(dp), dimension(nvbs) :: &
               dA_0_dvb, dA_1_dvb, dA_2_dvb, dA_numerator_dvb, dA_denom_dvb, &
               d_inv_sqrt_x_dvb
            
            real(qp) :: q1, q2, q3
            real(qp), dimension(nvbs) :: qd_dvb1, qd_dvb2, qd_dvb3
            
            real(dp), parameter :: two_13 = 1.2599210498948730d0 ! = pow(2d0,1d0/3d0)
            real(dp), parameter :: four_13 = 1.5874010519681994d0 ! = pow(4d0,1d0/3d0)

! options for MLT_option are:
!    'ML1'        Bohm-Vitense 1958 MLT
!    'ML2'        Bohm and Cassinelli 1971 MLT
!    'Mihalas'    Mihalas 1978, Kurucz 1979 MLT
!    'Henyey'     Henyey, Rardya, and Bodenheimer 1965 MLT
! Values of the f1..f4 coefficients are taken from Table 1 of Ludwig et al. 1999, A&A, 346, 111
! with the following exception: their value of f3 for Henyey convection is f4/8 when it should be
! 8*f4, i.e., f3=32*pi**2/3 and f4=4*pi**2/3. f3 and f4 are related to the henyey y parameter, so
! for the 'Henyey' case they are set based on the value of Henyey_y_param. The f1..f4 parameters
! have been renamed with a double ff, i.e., ff1..ff4, to avoid a collision of variable names with
! f1 in the cubic root solver.
            
            quit = .false.

            x = Q*Rho / (2*P)
            d_x_dvb = x*(drho_dvb/rho + dQ_dvb/Q - dP_dvb/P)
         
            convective_conductivity = Cp*grav*Lambda*Lambda*Rho*(sqrt(x)) / 9 ! erg / (K cm sec)
            
            if (convective_conductivity < 0) then
               ierr = -1
               return
               write(*,1) 'MLT error: convective_conductivity', convective_conductivity
               write(*,1) 'Cp', Cp
               write(*,1) 'grav', grav
               write(*,1) 'Lambda', Lambda
               write(*,1) 'Rho', Rho
               write(*,1) 'x', x
               call mesa_error(__FILE__,__LINE__)
            end if
            
            if (debug) write(*,1) 'convective_conductivity', convective_conductivity
            d_cc_dvb = convective_conductivity* &
                 (d_Cp_dvb/Cp + d_grav_dvb/grav + &
                     2*d_Lambda_dvb/Lambda + dRho_dvb/rho + d_x_dvb/(2*x))

            if (MLT_option == 'Cox') then ! this assumes optically thick

               a0 = 9d0/4d0
               d_a0_dvb = 0d0

               ! 'A' param is ratio of convective to radiative conductivities   C&G 14.98
               A = convective_conductivity / radiative_conductivity !  unitless.

               if (debug) write(*,1) 'A', A
               dA_dvb = (d_cc_dvb - d_rc_dvb*A) / radiative_conductivity
               
               if (A < 0 .or. is_bad(A)) then
                  write(*,*) "MLT_option == 'Cox'", MLT_option == 'Cox'
                  write(*,1) 'A', A
                  write(*,1) 'convective_conductivity', convective_conductivity
                  write(*,1) 'radiative_conductivity', radiative_conductivity
                  call mesa_error(__FILE__,__LINE__)
               end if

            else
            
               select case(trim(MLT_option))
               case ('Henyey')
                  ff1=1.0d0/Henyey_nu_param
                  ff2=0.5d0 ! 1d0/2.
                  ! popular values for y are 1/3 or 3/(4*pi**2)
                  ff3=8.0d0/Henyey_y_param
                  ff4=1.0d0/Henyey_y_param
               case ('ML1')
                  ff1=0.125d0 ! 1/8
                  ff2=0.5d0 ! 1/2
                  ff3=24.0d0
                  ff4=0.0d0
               case ('ML2')
                  ff1=1.0d0
                  ff2=2.0d0
                  ff3=16.0d0
                  ff4=0.0d0
               case ('Mihalas')
                  ff1=0.125d0 ! 1/8
                  ff2=0.5d0 ! 1/2
                  ff3=16.0d0
                  ff4=2.0d0
               case default
                  write(*,'(3a)') 'Error: ',trim(MLT_option), &
                     ' is not an allowed MLT version for convection'
                  write(*,*)
                  return
               end select
            
               omega = Lambda*Rho*opacity !dimensionless
               d_omega_dvb = omega*( d_Lambda_dvb/Lambda + dRho_dvb/Rho + d_opacity_dvb/opacity)

               ! the variable theta in no longer needed
               ! theta = omega / ( 1d0 + Henyey_y_param*omega**2 )
               ! d_theta_dvb = d_omega_dvb*(1d0 - Henyey_y_param*omega**2 ) /
               ! ( ( 1d0 + Henyey_y_param*omega**2 )**2 )

               ! a0 = 0.75d0*omega*theta
               !d_a0_dvb = a0*( d_omega_dvb/omega + d_theta_dvb/theta )
               a0 = (3d0/16d0)*ff2*ff3/(1d0+ff4/(omega*omega))
               d_a0_dvb = a0*2*ff4*d_omega_dvb/(ff4 + omega*omega)/omega
               !ignore d_a0_dvb
               d_a0_dvb = 0d0

               ! A = sqrt(P*Q*rho/Henyey_nu_param)*(Cp*mixing_length_alpha)/
               !        (2*crad*clight*T**3*theta)               
               A_0 = sqrt(ff1*P*Q*rho)
               dA_0_dvb = ff1*(dP_dvb*Q*rho + P*dQ_dvb*rho + P*Q*drho_dvb)/(2*A_0)
               
               A_1 = 4*A_0*Cp
               dA_1_dvb = 4*(dA_0_dvb*Cp + A_0*d_Cp_dvb)
               
               A_2 = mixing_length_alpha*omega*(1.d0+ff4/(omega*omega))
               dA_2_dvb = mixing_length_alpha*(1-ff4/(omega*omega))*d_omega_dvb

               if (is_bad(dA_2_dvb(mlt_dlnT00))) then
                  ierr = -1
                  write(*,1) 'dA_2_dvb(mlt_dlnT00)', dA_2_dvb(mlt_dlnT00)
                  call mesa_error(__FILE__,__LINE__)
                  return
               end if
               
               A_numerator = A_1*A_2
               dA_numerator_dvb = dA_1_dvb*A_2 + A_1*dA_2_dvb

               if (is_bad(dA_numerator_dvb(mlt_dlnT00))) then
                  ierr = -1
                  return
                  write(*,1) 'dA_numerator_dvb(mlt_dlnT00)', dA_numerator_dvb(mlt_dlnT00)
                  call mesa_error(__FILE__,__LINE__)
               end if
                     
               A_denom = ff3*crad*clight*T*T*T
               dA_denom_dvb = A_denom*3*dT_dvb/T
               
               A = A_numerator/A_denom                     
               dA_dvb = dA_numerator_dvb/A_denom - A_numerator*dA_denom_dvb/(A_denom*A_denom)

               if (is_bad(dA_dvb(mlt_dlnT00))) then
                  ierr = -1
                  return
                  write(*,1) 'dA_dvb(mlt_dlnT00)', dA_dvb(mlt_dlnT00)
                  call mesa_error(__FILE__,__LINE__)
               end if
               
               if (A < 0) then
                  ierr = -1
                  return
                  write(*,1) 'A', A
                  write(*,1) 'A_numerator', A_numerator
                  write(*,1) 'A_denom', A_denom
                  write(*,1) 'A_1', A_1
                  write(*,1) 'ff3', ff3
                  write(*,1) 'A_0', A_0
                  call mesa_error(__FILE__,__LINE__)
               end if
            
            end if

            ! 'B' param  C&G 14.81
            Bcubed = (A*A / a0)*diff_grads         
            d_Bcubed_dvb = (A*A / a0)*d_diff_grads_dvb + (2*A*dA_dvb / a0)*diff_grads &
               - Bcubed*d_a0_dvb/a0

            if (is_bad(d_Bcubed_dvb(mlt_dlnT00))) then
               ierr = -1
               return
!$omp critical (mlt_info_crit10)
               write(*,1) 'd_diff_grads_dvb(mlt_dlnT00)', d_diff_grads_dvb(mlt_dlnT00)
               write(*,1) 'dA_dvb(mlt_dlnT00)', dA_dvb(mlt_dlnT00)
               write(*,1) 'd_Bcubed_dvb(mlt_dlnT00)', d_Bcubed_dvb(mlt_dlnT00)
               write(*,1) 'a0', a0
               write(*,1) 'diff_grads', diff_grads
               write(*,1) 'A', A
               write(*,1) 'Bcubed', Bcubed
               call mesa_error(__FILE__,__LINE__)
!$omp end critical (mlt_info_crit10)
            end if
         
            if (debug) write(*,1) 'Bcubed', Bcubed

            ! now solve cubic equation for convective efficiency, Gamma
            ! a0*Gamma^3 + Gamma^2 + Gamma - a0*Bcubed == 0   C&G 14.82, 
            ! rewritten in terms of Gamma
            ! leave it to Mathematica to find an expression for the root we want (with a0 = 9/4)
         
            delta = a0*Bcubed
            d_delta_dvb = a0*d_Bcubed_dvb + Bcubed*d_a0_dvb
         
            if (debug) write(*,1) 'a0', a0
            if (debug) write(*,1) 'delta', delta
      
            f = -2 + 9*a0 + 27*a0*a0*delta
            if (debug) write(*,1) 'f', f
            if (f > 1d100) then
               f0 = f
               d_f0_dvb = 27*a0*a0*d_delta_dvb + (9+54*a0*delta)*d_a0_dvb
            else
               f0 = f*f + 4*(-1 + 3*a0)*(-1 + 3*a0)*(-1 + 3*a0)
               if (f0 <= 0d0) then
                  f0 = f
                  d_f0_dvb = 27*a0*a0*d_delta_dvb + (9+54*a0*delta)*d_a0_dvb
               else
                  f0 = sqrt(f0)         
                  !d_f0_dvb = (f*(27*a0*a0*d_delta_dvb + (9+54*a0*delta)*d_a0_dvb) &
                  !   + 18*(-1 + 3*a0)*(-1 + 3*a0)*d_a0_dvb)/f0!27*a0*a0*f*d_delta_dvb / f0
                  d_f0_dvb = 27*a0*a0*f*d_delta_dvb / f0 &
                     + (f*(9+54*a0*delta)+ 18*(-1 + 3*a0)*(-1 + 3*a0))*d_a0_dvb/f0
               end if
            end if
         
            if (debug) write(*,1) 'f0', f0

            f1 = -2 + 9*a0 + 27*a0*a0*delta + f0  

            if (is_bad(f1)) then
               ierr = -1
               return
!$omp critical (mlt_info_crit11)
               write(*,1) 'f1', f1
               write(*,1) 'a0', a0
               write(*,1) 'delta', delta
               write(*,1) 'f0', f0
               stop 'MLT: bad f1'
!$omp end critical (mlt_info_crit11)
            end if   

            if (f1 < 0) then
               call set_no_mixing
               call time_limit_conv_vel
               return
            end if   
            f1 = pow(f1,one_third)     
            d_f1_dvb = (27*a0*a0*d_delta_dvb + (9+54*a0*delta)*d_a0_dvb + d_f0_dvb) / (3*f1*f1)
            f2 = 2*two_13*(1 - 3*a0) / f1       
            d_f2_dvb = -f2*d_f1_dvb / f1 - 6*two_13*d_a0_dvb / f1

            Gamma = (four_13*f1 + f2 - 2) / (6*a0)
            d_Gamma_dvb = (four_13*d_f1_dvb + d_f2_dvb) / (6*a0) - Gamma*d_a0_dvb/a0

            if (is_bad(Gamma)) then
               ierr = -1
               return
!$omp critical (mlt_info_crit12)
               write(*,1) 'Gamma', Gamma
               write(*,1) 'f1', f1
               write(*,1) 'f2', f2
               write(*,1) 'a0', a0
               write(*,1) 'd_f1_dvb', d_f1_dvb
               write(*,1) 'd_f2_dvb', d_f2_dvb
               stop 'MLT: bad f1'
!$omp end critical (mlt_info_crit12)
!               call set_no_mixing
!               quit = .true.
!               return
            end if

            if (Gamma < 0) then
               call set_no_mixing
               call time_limit_conv_vel
               return
            end if
            
            ! average convection velocity, vbar   C&G 14.86b
            ! vbar = vsound*Sqrt(Q)*alpha*Gamma / (2*Sqrt(2*Gamma1)*A)
            ! vsound = Sqrt(Gamma1*P / rho), so
            ! vbar = Sqrt(Q*P / (8*rho))*alpha*Gamma / A

            x = Q*P / (8*rho)
            sqrt_x = sqrt(x)
            conv_vel = mixing_length_alpha*sqrt_x*Gamma / A
            if (conv_vel > max_conv_vel) then
               conv_vel = max_conv_vel
               d_conv_vel_dvb = 0
            else
               d_conv_vel_dvb = 0.5d0*conv_vel* &
                 (-2*dA_dvb/A + 2*d_Gamma_dvb/Gamma + &
                 dP_dvb/P + dQ_dvb/Q - drho_dvb/rho)
            end if

            conv_tau = 1d99
            call time_limit_conv_vel
            
            if (init_conv_vel /= conv_vel .or. conv_vel == max_conv_vel) then
               ! need to recalculate Gamma to match modified conv_vel
               if (A <= 1d-99 .or. sqrt_x <= 1d-99) then
                  Gamma = 1d25
                  d_Gamma_dvb = 0
                  if (dbg) write(*,1) 'A or sqrt_x too small', A, sqrt_x
               else      
                  if (dbg) write(*,*) 'recalculate Gamma to match modified conv_vel'
                  inv_sqrt_x = 1d0/sqrt_x
                  d_inv_sqrt_x_dvb = &
                     (Q*P*drho_dvb/rho - Q*dP_dvb - P*dQ_dvb)/(16d0*x*sqrt_x*rho)
                  Gamma = conv_vel*A*inv_sqrt_x/mixing_length_alpha  
                  !d_Gamma_dvb = ( &
                  !   d_conv_vel_dvb*A*inv_sqrt_x + &
                  !   conv_vel*dA_dvb*inv_sqrt_x + &
                  !   conv_vel*A*d_inv_sqrt_x_dvb)/mixing_length_alpha  
                  ! the "correct" form breaks edep in example_ccsn_IIp.
                  d_Gamma_dvb = 0d0
               end if
            end if

            if (dbg) &
               write(*,1) 'prev/init init/final conv_vel', &
                  prev_conv_vel/init_conv_vel, &
                  prev_conv_vel, init_conv_vel, conv_vel
            
            if (debug) write(*,1) 'conv_vel', conv_vel
            if (conv_vel < 0) then
               ierr = -1
               return
!$omp critical (mlt_info_crit13)
               write(*,1) 'conv_vel', conv_vel
               write(*,1) 'mixing_length_alpha', mixing_length_alpha
               write(*,1) 'x', x
               write(*,1) 'A', A
               write(*,1) 'Gamma', Gamma
               stop 'MLT: set_convective_mixing'
!$omp end critical (mlt_info_crit13)
            end if
            
            Zeta = Gamma*Gamma*Gamma/Bcubed  ! C&G 14.80
            
            ! quad didn't help for the case i was trying to fix.      
            !q1 = Gamma; qd_dvb1 = d_Gamma_dvb
            !q2 = Bcubed; qd_dvb2 = d_Bcubed_dvb
            !q3 = Zeta
            !qd_dvb3 = (3d0*q1*q1*qd_dvb1 - qd_dvb2*q3)/q2
            !d_Zeta_dvb = qd_dvb3
            
            d_Zeta_dvb = (3d0*Gamma*Gamma*d_Gamma_dvb - d_Bcubed_dvb*Zeta)/Bcubed
            
            
            !Zeta = a0*Gamma*Gamma/(1+Gamma*(1+a0*Gamma)) ! C&G, eq. (14.78)
            !d_Zeta_dvb = d_Gamma_dvb*(((Gamma**2+2*Gamma)*a0)/(pow4(Gamma)*a0**2&
            !   +(2*pow3(Gamma)+2*Gamma**2)*a0+Gamma**2+2*Gamma+1)) + &
            !   d_a0_dvb*(pow2(Gamma)/(Gamma*(Gamma*a0+1)+1)&
            !             -(pow4(Gamma)*a0)/pow2(Gamma*(Gamma*a0+1)+1))

            !write(*,*) "Compare Zeta", kz, log10(T), Zeta, &
            !   a0*Gamma*Gamma/(1+Gamma*(1+a0*Gamma)) ! C&G, eq. (14.78)
            
            ! Zeta must be >= 0 and < 1
            if (Zeta < 0d0) then
               Zeta = 0
               d_Zeta_dvb = 0
            else if (Zeta >= 1d0) then
               Zeta = 1d0
               d_Zeta_dvb = 0
            end if
            
            !gradT = (1d0 - Zeta)*gradr + Zeta*gradL ! C&G 14.79 with gradL for grada
            !d_gradT_dvb = (1d0 - Zeta)*d_gradr_dvb + Zeta*d_gradL_dvb + &
            !            (gradL - gradr)*d_Zeta_dvb
            gradT = (1d0 - Zeta)*gradr + Zeta*grada ! C&G 14.79
            d_gradT_dvb = (1d0 - Zeta)*d_gradr_dvb + Zeta*d_grada_dvb + &
                        (grada - gradr)*d_Zeta_dvb

            

            if (test_partials) then
               ss% solver_test_partials_val = gradT
               ss% solver_test_partials_var = ss% i_lnd
               ss% solver_test_partials_dval_dx = d_gradT_dvb(mlt_dlnd00)
               write(*,*) 'mlt get_results', ss% solver_test_partials_var
            end if
            
            
            
            if (is_bad(gradT)) then
               call set_no_mixing
               quit = .true.
               return
            end if
         
         end subroutine set_convective_mixing   

         subroutine set_semiconvection ! Langer 1983 & 1985
            real(dp) :: alpha, bc, LG, &
               a0, a1, a2, a3, a4, a5, a6, a, &
               b1, b2, b3, b4, b5, b6, b7, b, div, bsq
            real(dp), dimension(nvbs) :: &
               d_bc_dvb, d_LG_dvb, d_a0_dvb, d_a1_dvb, d_a2_dvb, d_a3_dvb, d_a4_dvb, &
               d_a5_dvb, d_a6_dvb, d_a_dvb, d_b1_dvb, d_b2_dvb, d_b3_dvb, d_b4_dvb, &
               d_b5_dvb, d_b6_dvb, d_b7_dvb, d_b_dvb, d_div_dvb
            
            include 'formats.dek'
            if (dbg) write(*,*) 'check for semiconvection'
            call set_no_mixing ! sets gradT = gradr
            D_semi = alpha_semiconvection*radiative_conductivity/(6*Cp*rho) &
                  *(gradr - grada)/(gradL - gradr)
            if (D_semi <= 0) then
               if (dbg) then
                  write(*,1) 'set_no_mixing D_semi', D_semi
                  write(*,1) 'alpha_semiconvection', alpha_semiconvection
                  write(*,1) 'radiative_conductivity', radiative_conductivity
                  write(*,1) 'gradr - grada', gradr - grada
                  write(*,1) 'gradL - gradr', gradL - gradr
                  stop
               end if
               call set_no_mixing
               return
            end if
            d_D_semi_dvb = 0 ! not used, so skip for now.
            conv_vel = 3*D_semi/Lambda 
            d_conv_vel_dvb = 0
            call time_limit_conv_vel
            if (conv_vel > max_conv_vel) conv_vel = max_conv_vel
            if (init_conv_vel /= conv_vel) D_semi = conv_vel*Lambda/3
            if (D_semi <= 0) then
               call set_no_mixing
               return
            end if
            
            mixing_type = semiconvective_mixing
            if (dbg) write(*,2) 'mixing_type', mixing_type
            
            if (semiconvection_option == 'Langer_85 mixing; gradT = gradr') return
            if (semiconvection_option /= 'Langer_85') then
               write(*,*) 'MLT: unknown values for semiconvection_option ' // &
                  trim(semiconvection_option)
               ierr = -1
               return
            end if
            
            
!            Solve[{
!                  L/Lrad - Lsc/Lrad - 1 == 0, 
!                  Lrad == grad LG, 
!                  gradMu == (4 - 3*beta)/beta*gradL_composition_term,
!                  Lsc/Lrad == alpha (grad - gradA)/(2 grad (gradL - grad))
!                              (grad - gradA - (beta (8 - 3 beta))/bc gradMu)}, 
!                  grad, {Lsc, Lrad, gradMu}] // Simplify
                  
            alpha = min(1d0, alpha_semiconvection)

            bc = 32 - 24*beta - beta*beta
            d_bc_dvb = - 24*d_beta_dvb - 2*d_beta_dvb*beta
            
            LG = (16d0/3d0*pi*clight*m*cgrav*crad*T*T*T*T)/(P*opacity)
            d_LG_dvb = LG*(4d0*dT_dvb/T - dP_dvb/P - d_opacity_dvb/opacity)
            
            a0 = alpha*gradL_composition_term*LG
            d_a0_dvb = alpha*gradL_composition_term*d_LG_dvb
            
            a1 = -2*bc*L
            d_a1_dvb = -2*L*d_bc_dvb
            d_a1_dvb(mlt_dL) = d_a1_dvb(mlt_dL) - 2*bc
            
            a2 = 2*alpha*bc*grada*LG
            d_a2_dvb = 2*alpha*(d_bc_dvb*grada*LG + bc*d_grada_dvb*LG + bc*grada*d_LG_dvb)
            
            a3 = -2*bc*gradL*LG
            d_a3_dvb = -2*(d_bc_dvb*gradL*LG + bc*d_gradL_dvb*LG + bc*gradL*d_LG_dvb)
            
            a4 = 32*a0
            d_a4_dvb = 32*d_a0_dvb
            
            a5 = -36*beta*a0
            d_a5_dvb = -36*(d_beta_dvb*a0 + beta*d_a0_dvb)
            
            a6 = 9*beta*beta*a0
            d_a6_dvb = 9*(2*beta*d_beta_dvb*a0 + beta*beta*d_a0_dvb)
            
            a = a1 + a2 + a3 + a4 + a5 + a6
            d_a_dvb = d_a1_dvb + d_a2_dvb + d_a3_dvb + d_a4_dvb + d_a5_dvb + d_a6_dvb 
                           
            b1 = 32 - 36*beta + 9*beta*beta
            d_b1_dvb = - 36*d_beta_dvb + 18*beta*d_beta_dvb
            
            b2 = b1*a0
            d_b2_dvb = d_b1_dvb*a0 + b1*d_a0_dvb
            
            b3 = -2*gradL*L + alpha*grada*grada*LG
            d_b3_dvb = -2*d_gradL_dvb*L + alpha*(2*grada*d_grada_dvb*LG + grada*grada*d_LG_dvb)
            d_b3_dvb(mlt_dL) = d_b3_dvb(mlt_dL) - 2*gradL
            
            b4 = (-alpha*gradA + gradL)*LG
            d_b4_dvb = (-alpha*d_grada_dvb + d_gradL_dvb)*LG + (-alpha*gradA + gradL)*d_LG_dvb
            
            b5 = -b2 + 2*bc*(L + b4)
            d_b5_dvb = -d_b2_dvb + 2*d_bc_dvb*(L + b4) + 2*bc*d_b4_dvb
            d_b5_dvb(mlt_dL) = d_b5_dvb(mlt_dL) + 2*bc
            
            b6 = b2*grada + bc*b3
            d_b6_dvb = d_b2_dvb*grada + b2*d_grada_dvb + d_bc_dvb*b3 + bc*d_b3_dvb
            
            b7 = -4*(-2 + alpha)*bc*LG*b6
            d_b7_dvb = -4*(-2 + alpha)*(d_bc_dvb*LG*b6 + bc*d_LG_dvb*b6 + bc*LG*d_b6_dvb)
            
            b = b7 + b5*b5
            d_b_dvb = d_b7_dvb + 2*b5*d_b5_dvb
            
            div = 2*(-2 + alpha)*bc*LG
            d_div_dvb = 2*(-2 + alpha)*(d_bc_dvb*LG + bc*d_LG_dvb)

            bsq = sqrt(b)
            gradT = (a + bsq)/div
            d_gradT_dvb = -gradT*d_div_dvb/div + d_a_dvb/div + 0.5d0*d_b_dvb/(div*bsq)
            
         end subroutine set_semiconvection
                  
         subroutine set_no_mixing
            ! assumes have set gradr, scale_height, gradL, and Lambda.
            mixing_type = no_mixing
            gradT = gradr
            d_gradT_dvb = d_gradr_dvb
            conv_vel = 0
            d_conv_vel_dvb = 0
            D = 0
            d_D_dvb = 0
            D_semi = 0
            d_D_semi_dvb = 0
            D_thrm = 0
            d_D_thrm_dvb = 0
            Gamma = 0
            d_Gamma_dvb = 0
         end subroutine set_no_mixing         
         
         subroutine show_args
 1          format(a30,1pe26.16)
            
            write(*,1) 'cgrav = ', cgrav
            write(*,1) 'm = ', m
            write(*,1) 'r = ', r 
            write(*,1) 'T = ', T 
            write(*,1) 'Rho = ', Rho 
            write(*,1) 'L  = ', L 
            write(*,1) 'P = ', P
            write(*,1) 'chiRho = ', chiRho 
            write(*,1) 'chiT = ', chiT
            write(*,1) 'Cp = ', Cp 
            write(*,1) 'xh = ', xh
            write(*,1) 'opacity = ', opacity 
            write(*,1) 'grada = ', grada
            write(*,1) 'mixing_length_alpha = ', mixing_length_alpha
            
         end subroutine show_args


         subroutine revise_using_cv_var_variable
         
            ! only changes D, d_D_dvb, gradT, d_gradT_dvb
            ! does NOT change any others such as mlt_vc
            
            include 'formats.dek'
            real(dp) ff1, ff2, ff3, ff4, sqrt_x, tmp
            real(dp) :: cv_var, A_0, A_1, A_2, A_numerator, A_denom, &
               inv_sqrt_x, save_gradT
            real(dp), dimension(nvbs) :: &
               dA_0_dvb, dA_1_dvb, dA_2_dvb, dA_numerator_dvb, dA_denom_dvb, &
               d_inv_sqrt_x_dvb, d_cv_var_dvb, d_save_dvb
            real(qp) :: q1, q2, q3
            real(qp), dimension(nvbs) :: dq1_dvb, dq2_dvb, dq3_dvb
            integer :: j

            quit = .false.
            if (kz == 0) return
            
           !!Pablo: TODO, not sure if this helps, use velocity from middle of the step
           !!NOTE, needed to change solver_vars as well to make conv_vel_start available
           if (ss% conv_vel_flag) then
              cv_var = 0.5d0*(ss% conv_vel(kz)+ss% conv_vel_start(kz))
              d_cv_var_dvb = 0
              d_cv_var_dvb(mlt_cv_var) = 0.5d0
           else if (ss% cv_flag) then
              cv_var = ss% cv(kz)
              d_cv_var_dvb = 0
              d_cv_var_dvb(mlt_cv_var) = 1d0
           end if
                        
            if (cv_var < 0d0 .or. is_bad(cv_var)) then
               ierr = -1
               if (ss% report_ierr) then
                  write(*,2) 'cv_var', kz, cv_var
                  if (ss% stop_for_bad_nums) stop 'revise_using_cv_var_variable'
               else
                  return
               end if
            end if
            
            d_D_dvb(mlt_cv_var) = 0d0
            d_gradT_dvb(mlt_cv_var) = 0d0

            D = cv_var*Lambda/3     ! diffusion coefficient [cm^2/sec]
            if (debug) write(*,1) 'D', D
            d_D_dvb = (d_cv_var_dvb*Lambda + cv_var*d_Lambda_dvb)/3

            x = Q*Rho / (2d0*P) ! using x as a temporary variable here
            d_x_dvb = x*(drho_dvb/rho + dQ_dvb/Q - dP_dvb/P)
         
            convective_conductivity = Cp*grav*Lambda*Lambda*Rho*(sqrt(x)) / 9 ! erg / (K cm sec)
                        
            if (is_bad(convective_conductivity)) then
               ierr = -1
               if (ss% report_ierr) then
                  write(*,2) 'convective_conductivity', kz, convective_conductivity
                  if (ss% stop_for_bad_nums) stop 'revise_using_cv_var_variable'
               else
                  return
               end if
            end if
            
            if (convective_conductivity < 0) then
               ierr = -1
               if (ss% report_ierr) then
                  write(*,1) 'MLT error: convective_conductivity', convective_conductivity
                  write(*,1) 'Cp', Cp
                  write(*,1) 'grav', grav
                  write(*,1) 'Lambda', Lambda
                  write(*,1) 'Rho', Rho
                  write(*,1) 'x', x
                  call mesa_error(__FILE__,__LINE__)
               else
                  return
               end if
            end if
            
            if (debug) write(*,1) 'convective_conductivity', convective_conductivity
            d_cc_dvb = convective_conductivity* &
                 (d_Cp_dvb/Cp + d_grav_dvb/grav + &
                     2*d_Lambda_dvb/Lambda + dRho_dvb/rho + d_x_dvb/(2*x))

            if (MLT_option == 'Cox') then ! this assumes optically thick

               a0 = 9d0/4d0
               d_a0_dvb = 0d0

               ! 'A' param is ratio of convective to radiative conductivities   C&G 14.98
               A = convective_conductivity / radiative_conductivity !  unitless.

               if (debug) write(*,1) 'A', A
               dA_dvb = (d_cc_dvb - d_rc_dvb*A) / radiative_conductivity
               
               if (A < 0 .or. is_bad(A)) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,*) "MLT_option == 'Cox'", MLT_option == 'Cox'
                     write(*,1) 'A', A
                     write(*,1) 'convective_conductivity', convective_conductivity
                     write(*,1) 'radiative_conductivity', radiative_conductivity
                     call mesa_error(__FILE__,__LINE__)
                  else
                     return
                  end if
               end if

            else
            
               select case(trim(MLT_option))
               case ('Henyey')
                  ff1=1d0/Henyey_nu_param
                  ff2=0.5d0
                  ! popular values for y are 1/3 or 3/(4*pi**2)
                  ff3=8d0/Henyey_y_param
                  ff4=1d0/Henyey_y_param
               case ('ML1')
                  ff1=1d0/8d0
                  ff2=1d0/2d0
                  ff3=24d0
                  ff4=0d0
               case ('ML2')
                  ff1=1d0
                  ff2=2d0
                  ff3=16d0
                  ff4=0d0
               case ('Mihalas')
                  ff1=1d0/8d0
                  ff2=1d0/2d0
                  ff3=16d0
                  ff4=2d0
               case default
                  write(*,'(3a)') 'Error: ',trim(MLT_option), &
                     ' is not an allowed MLT version for convection'
                  write(*,*)
                  return
               end select
            
               omega = Lambda*Rho*opacity !dimensionless
               d_omega_dvb = omega*( d_Lambda_dvb/Lambda + dRho_dvb/Rho + d_opacity_dvb/opacity)

               ! the variable theta in no longer needed
               ! theta = omega / ( 1d0 + Henyey_y_param*omega**2 )
               ! d_theta_dvb = d_omega_dvb*(1d0 - Henyey_y_param*omega**2 ) /
               ! ( ( 1d0 + Henyey_y_param*omega**2 )**2 )

               ! a0 = 0.75d0*omega*theta
               !d_a0_dvb = a0*( d_omega_dvb/omega + d_theta_dvb/theta )
               a0 = (3d0/16d0)*ff2*ff3/(1d0+ff4/(omega*omega))
               d_a0_dvb = a0*2d0*ff4*d_omega_dvb/(ff4 + omega*omega)/omega

               ! A = sqrt(P*Q*rho/Henyey_nu_param)*(Cp*mixing_length_alpha)/
               !        (2*crad*clight*T**3*theta)               
               A_0 = sqrt(ff1*P*Q*rho)
               dA_0_dvb = ff1*(dP_dvb*Q*rho + P*dQ_dvb*rho + P*Q*drho_dvb)/(2*A_0)
               
               A_1 = 4d0*A_0*Cp
               dA_1_dvb = 4d0*(dA_0_dvb*Cp + A_0*d_Cp_dvb)
               
               A_2 = mixing_length_alpha*omega*(1d0+ff4/(omega*omega))
               dA_2_dvb = mixing_length_alpha*(1d0-ff4/(omega*omega))*d_omega_dvb

               if (is_bad(dA_2_dvb(mlt_dlnT00))) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,1) 'dA_2_dvb(mlt_dlnT00)', dA_2_dvb(mlt_dlnT00)
                     call mesa_error(__FILE__,__LINE__)
                  else
                     return
                  end if
               end if
               
               A_numerator = A_1*A_2
               dA_numerator_dvb = dA_1_dvb*A_2 + A_1*dA_2_dvb

               if (is_bad(dA_numerator_dvb(mlt_dlnT00))) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,1) 'dA_numerator_dvb(mlt_dlnT00)', dA_numerator_dvb(mlt_dlnT00)
                     call mesa_error(__FILE__,__LINE__)
                  else
                     return
                  end if
               end if
                     
               A_denom = ff3*crad*clight*T*T*T
               dA_denom_dvb = A_denom*3d0*dT_dvb/T
               
               A = A_numerator/A_denom                     
               dA_dvb = dA_numerator_dvb/A_denom - A_numerator*dA_denom_dvb/(A_denom*A_denom)

               if (is_bad(dA_dvb(mlt_dlnT00))) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,1) 'dA_dvb(mlt_dlnT00)', dA_dvb(mlt_dlnT00)
                     call mesa_error(__FILE__,__LINE__)
                  else
                     return
                  end if
               end if
               
               if (A < 0) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,1) 'A', A
                     write(*,1) 'A_numerator', A_numerator
                     write(*,1) 'A_denom', A_denom
                     write(*,1) 'A_1', A_1
                     write(*,1) 'ff3', ff3
                     write(*,1) 'A_0', A_0
                     call mesa_error(__FILE__,__LINE__)
                  else
                     return
                  end if
               end if
            
            end if
                        
            if (is_bad(A)) then
               ierr = -1
               if (ss% report_ierr) then
                  write(*,2) 'A', kz, A
                  if (ss% stop_for_bad_nums) stop 'revise_using_cv_var_variable'
               else
                  return
               end if
            end if

            ! average convection velocity, vbar   C&G 14.86b
            ! cv_var = vsound*Sqrt(Q)*alpha*Gamma / (2*Sqrt(2*Gamma1)*A)
            ! vsound = Sqrt(Gamma1*P / rho), so
            ! cv_var = Sqrt(Q*P / (8*rho))*alpha*Gamma / A
            ! cv_var = sqrt_x*alpha*Gamma / A
            ! Gamma = cv_var*A/(sqrt_x*alpha)

            x = Q*P / (8d0*rho) ! using x as a temporary variable here
            sqrt_x = sqrt(x)
            
            if (A <= 1d-99 .or. sqrt_x <= 1d-99) then
               Gamma = 1d25
               d_Gamma_dvb = 0
               !if (dbg) &
                  write(*,2) 'mlt: A or sqrt_x too small', kz, A, sqrt_x
            else      
               if (dbg) write(*,*) 'calculate Gamma using cv_var'
               inv_sqrt_x = 1d0/sqrt_x
               d_inv_sqrt_x_dvb = &
                  (Q*P*drho_dvb/rho - Q*dP_dvb - P*dQ_dvb)/(16d0*x*sqrt_x*rho)
               if (.false.) then
                  write(*,2) 'rel_diff new old Gamma', kz, &
                     (cv_var*A*inv_sqrt_x/mixing_length_alpha - Gamma)/Gamma, &
                     cv_var*A*inv_sqrt_x/mixing_length_alpha, Gamma
                  write(*,2) 'old d_Gamma_dvb(mlt_dL)', kz, d_Gamma_dvb(mlt_dL)
                  write(*,2) 'old d_Gamma_dvb(mlt_dlnR)', kz, d_Gamma_dvb(mlt_dlnR)
               end if
               Gamma = cv_var*A*inv_sqrt_x/mixing_length_alpha  
               d_Gamma_dvb = ( &
                  d_cv_var_dvb*A*inv_sqrt_x + &
                  cv_var*dA_dvb*inv_sqrt_x + &
                  cv_var*A*d_inv_sqrt_x_dvb)/mixing_length_alpha  
               if (is_bad(Gamma)) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,2) 'Gamma', kz, Gamma
                     write(*,2) 'ss% conv_vel', kz, ss% conv_vel(kz)
                     write(*,2) 'cv_var', kz, cv_var
                     write(*,2) 'A', kz, A
                     write(*,2) 'inv_sqrt_x', kz, inv_sqrt_x
                     write(*,2) 'x', kz, x
                     if (ss% stop_for_bad_nums) stop 'revise_using_cv_var_variable'
                  else
                     return
                  end if
               end if
            end if

                        
            
            if (.true.) then
            
               ! C&G, eq. (14.78), but rewritten in a way that prevents
               ! the multiplication of terms that go as ~Gamma. This is because
               ! Gamma can be very large, and Gamma^2 can actually overflow and
               ! produce a NaN.
               tmp = 1d0/(1d0 + a0*Gamma - a0*Gamma/(Gamma+1d0))
               Zeta = 1d0 - tmp
               d_Zeta_dvb =d_Gamma_dvb*(((a0*Gamma/(Gamma+1d0))*tmp)&
                  *((Gamma+2d0)/(Gamma+1d0)*tmp)) + &
                  d_a0_dvb*((Gamma/(Gamma+1))*(Gamma*tmp)*tmp)
            
            else if (.false.) then ! quad precision
            
               q1 = Gamma
               dq1_dvb = d_Gamma_dvb
               q2 = a0
               q3 = q2*q1*q1/(1+q1*(1+q2*q1)) ! C&G, eq. (14.78)
               dq3_dvb = dq1_dvb*(((q1*q1+2*q1)*q2)/(q1*q1*q1*q1*q2*q2&
                  +(2*q1*q1*q1+2*q1*q1)*q2+q1*q1+2*q1+1))
               Zeta = q3
               d_Zeta_dvb = dq3_dvb
               
            else

               Zeta = a0*Gamma*Gamma/(1+Gamma*(1+a0*Gamma)) ! C&G, eq. (14.78)
               d_Zeta_dvb = d_Gamma_dvb*(((Gamma*Gamma+2*Gamma)*a0)/(pow4(Gamma)*a0*a0&
                  +(2*pow3(Gamma)+2*Gamma*Gamma)*a0+Gamma*Gamma+2*Gamma+1))

            end if

                        
            if (is_bad(Zeta)) then
               ierr = -1
               if (ss% report_ierr) then
                  write(*,2) 'Zeta', kz, Zeta
                  if (ss% stop_for_bad_nums) stop 'revise_using_cv_var_variable'
               else
                  return
               end if
            end if
            
            save_gradT = gradT
            d_save_dvb = d_gradT_dvb        

            gradT = (1d0 - Zeta)*gradr + Zeta*grada ! C&G 14.79
            d_gradT_dvb = (1d0 - Zeta)*d_gradr_dvb + Zeta*d_grada_dvb + &
                        (grada - gradr)*d_Zeta_dvb
            
            if (.false. .and. Zeta > 0.45d0 .and. Zeta < 0.55d0) then
!$OMP critical (mlt_info_crit14)
               write(*,2) 'Zeta', kz, Zeta
               write(*,2) 'rel_diff new old gradT', kz, (gradT - save_gradT)/gradT, gradT, save_gradT
               do j=1,num_mlt_partials
                  tmp = d_gradT_dvb(mlt_cv_var)*d_conv_vel_dvb(j) + d_gradT_dvb(j)
                  write(*,3) 'rel_diff new old gradT partial ' // trim(mlt_partial_str(j)), j, kz, &
                     (tmp-d_save_dvb(j))/tmp, tmp, d_save_dvb(j)
               end do
               stop 'revise_using_cv_var_variable'
!$OMP end critical (mlt_info_crit14)
            end if

            if (is_bad(gradT)) then
               ierr = -1
               if (ss% report_ierr) then
                  write(*,2) 'gradT', kz, gradT
                  if (ss% stop_for_bad_nums) stop 'revise_using_cv_var_variable'
               else
                  return
               end if
            end if
         
         end subroutine revise_using_cv_var_variable


      end subroutine Get_results


      end module mlt_get_results
