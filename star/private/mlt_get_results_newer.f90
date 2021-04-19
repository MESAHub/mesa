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


      module mlt_get_results_newer

      use star_private_def
      use const_def
      use num_lib
      use utils_lib
      use auto_diff_support
      use star_utils

      implicit none
      
      private
      public :: do1_mlt_eval_newer, Get_results_newer, get_gradT_newer

      contains
      
      
      subroutine get_gradT_newer(s, MLT_option, & ! used to create models
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            gradT, mixing_type, ierr)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: MLT_option
         real(dp), intent(in) :: &
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
            XH1, cgrav, m, gradL_composition_term, mixing_length_alpha
         integer, intent(in) :: iso
         real(dp), intent(out) :: gradT
         integer, intent(out) :: mixing_type, ierr 
         type(auto_diff_real_star_order1) :: &
            gradr_ad, grada_ad, scale_height_ad, gradT_ad, Y_face_ad, mlt_vc_ad, D_ad, &
            Gamma_ad, r_ad, L_ad, T_ad, P_ad, opacity_ad, rho_ad, chiRho_ad, chiT_ad, Cp_ad
         ierr = 0
         r_ad = r
         L_ad = L
         T_ad = T
         P_ad = P
         opacity_ad = opacity
         rho_ad = rho
         chiRho_ad = chiRho
         chiT_ad = chiT
         Cp_ad = Cp
         gradr_ad = gradr
         grada_ad = grada
         scale_height_ad = scale_height
         if (s% use_other_mlt) then
            !call s% other_mlt()
         else         
            call Get_results_newer(s, 0, MLT_option, &
               r_ad, L_ad, T_ad, P_ad, opacity_ad, rho_ad, chiRho_ad, chiT_ad, Cp_ad, &
               gradr_ad, grada_ad, scale_height_ad, &
               iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
               s% alpha_semiconvection, s% thermohaline_coeff, s% using_TDC, &
               mixing_type, gradT_ad, Y_face_ad, mlt_vc_ad, D_ad, Gamma_ad, ierr)
         end if
         gradT = gradT_ad%val
      end subroutine get_gradT_newer
      
         
      subroutine do1_mlt_eval_newer(s, k, MLT_option, just_gradr, gradL_composition_term, &
            gradr, grada, scale_height, mixing_length_alpha, &
            mixing_type, gradT, Y_face, mlt_vc, D, Gamma, ierr)
         use chem_def, only: ih1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         character (len=*), intent(in) :: MLT_option
         logical, intent(in) :: just_gradr
         type(auto_diff_real_star_order1), intent(in) :: &
            gradr, grada, scale_height
         real(dp), intent(in) :: &
            gradL_composition_term, mixing_length_alpha
         integer, intent(out) :: mixing_type
         type(auto_diff_real_star_order1), intent(out) :: &
            gradT, Y_face, mlt_vc, D, Gamma
         integer, intent(out) :: ierr 
                 
         real(dp) :: cgrav, m, XH1, gradL_old, grada_face_old, alpha_semiconvection
         integer :: iso, old_mix_type
         type(auto_diff_real_star_order1) :: &
            r, L, T, P, opacity, rho, chiRho, chiT, Cp
         include 'formats'
         ierr = 0
         if (just_gradr) then
            mixing_type = no_mixing
            gradT = gradr
            Y_face = gradT - (grada+gradL_composition_term)
            mlt_vc = 0d0
            D = 0d0
            Gamma = 0d0
            return
         end if
         cgrav = s% cgrav(k)
         m = s% m_grav(k)
         L = wrap_L_00(s,k)
         T = get_T_face(s,k)
         P = get_Peos_face(s,k)
         r = wrap_r_00(s,k)
         opacity = get_kap_face(s,k)
         rho = get_Rho_face(s,k)
         chiRho = get_ChiRho_face(s,k)
         chiT = get_ChiT_face(s,k)
         Cp = get_Cp_face(s,k)
         iso = s% dominant_iso_for_thermohaline(k)
         XH1 = s% xa(s% net_iso(ih1),k)
         alpha_semiconvection = s% alpha_semiconvection
         if (s% center_h1 > s% semiconvection_upper_limit_center_h1) alpha_semiconvection = 0
         if (s% use_other_mlt) then
            !call s% other_mlt(s% id, k, &               
            !   gradr_factor, gradL_composition_term, &
            !   mixing_length_alpha, alpha_semiconvection, &
            !   mixing_type, gradT, gradr, mlt_vc, D, Gamma, scale_height, ierr)
         else         
            call Get_results_newer(s, k, MLT_option, &
               r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
               iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
               alpha_semiconvection, s% thermohaline_coeff, s% using_TDC, &
               mixing_type, gradT, Y_face, mlt_vc, D, Gamma, ierr)
         end if

      end subroutine do1_mlt_eval_newer


      subroutine Get_results_newer(s, k, MLT_option, &  ! NOTE: k=0 is a valid arg
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            alpha_semiconvection, thermohaline_coeff, using_TDC, &
            mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
         use star_utils
         type (star_info), pointer :: s
         integer, intent(in) :: k
         character (len=*), intent(in) :: MLT_option
         type(auto_diff_real_star_order1), intent(in) :: &
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height
         integer, intent(in) :: iso
         real(dp), intent(in) :: &
            XH1, cgrav, m, gradL_composition_term, &
            mixing_length_alpha, alpha_semiconvection, thermohaline_coeff
         logical, intent(in) :: using_TDC
         integer, intent(out) :: mixing_type
         type(auto_diff_real_star_order1), intent(out) :: &
            gradT, Y_face, conv_vel, D, Gamma
         integer, intent(out) :: ierr
         
         type(auto_diff_real_star_order1) :: &
            Pr, Pg, grav, scale_height2, Lambda, gradL, beta, Y_guess
         character (len=256) :: message        
         logical ::  okay_to_use_TDC, test_partials, compare_TDC_to_MLT
         include 'formats'
         
         compare_TDC_to_MLT = .false.

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         ierr = 0          
         
         mixing_type = no_mixing        
         Pr = crad*pow4(T)/3d0
         Pg = P - Pr
         beta = Pg / P
         ! Ledoux temperature gradient (same as Schwarzschild if composition term = 0)
         gradL = grada + gradL_composition_term
         grav = cgrav*m/pow2(r)
         
         if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
               s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
            write(*,2) 'enter Get_results_newer gradr grada scale_height ' // trim(MLT_option), k, &
               gradr%val, grada%val, scale_height%val
         end if
               
         
         call set_no_mixing('')
         if (MLT_option == 'none' .or. beta < 1d-10 .or. mixing_length_alpha <= 0d0) then
            return
         end if

         if (opacity%val < 1d-10 .or. P%val < 1d-20 .or. T%val < 1d-10 .or. Rho%val < 1d-20 &
               .or. m < 1d-10 .or. r%val < 1d-10 .or. cgrav < 1d-10) then
            call set_no_mixing('vals too smqll')
            return
         end if
           
         conv_vel = 0d0 ! will be set below if needed
         if (k > 0) then ! will be set below if TDC
            s% SOURCE(k) = 0d0
            s% DAMP(k) = 0d0
            s% DAMPR(k) = 0d0
            s% COUPL(k) = 0d0
         end if
         
         ! check if this particular k can be done with TDC
         if (using_TDC .and. k > 0 .and. s% dt > 0d0) then
            okay_to_use_TDC = (s% X(k) <= s% max_X_for_TDC)
         else
            okay_to_use_TDC = .false.
         end if
         
         if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
            write(*,2) 'gradr gradL grada comp_term', k, gradr%val, gradL%val, grada%val, gradL_composition_term
         end if
         Lambda = mixing_length_alpha*scale_height
         if (okay_to_use_TDC .and. .not. compare_TDC_to_MLT) then 
            ! this means TDC must do all types:
            ! convective, radiative, thermohaline, semiconvective
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'call set_TDC', k
            end if
            call set_TDC
         else if (gradr > gradL) then ! convective
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'call set_MLT', k
            end if
            call set_MLT
         else if (gradL_composition_term < 0) then
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'call set_thermohaline', k
            end if
            call set_thermohaline
         else if (gradr > grada) then
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'call set_semiconvection', k
            end if
            call set_semiconvection
         end if         
         
         if (D%val < s% remove_small_D_limit .or. is_bad(D%val)) then
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'D < s% remove_small_D_limit', k, D%val, s% remove_small_D_limit
            end if
            mixing_type = no_mixing
         end if
         
         if (mixing_type == no_mixing) call set_no_mixing('final mixing_type == no_mixing')
         
         if (okay_to_use_TDC .and. compare_TDC_to_MLT) call do_compare_TDC_to_MLT
         
         contains


         subroutine set_no_mixing(str)
            character (len=*) :: str
            include 'formats'            
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                  s% solver_iter == s% x_integer_ctrl(20) &
                  .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0) &
                  .and. len_trim(str) > 0) then
               write(*,2) 'Get_results_newer set_no_mixing ' // trim(str), k
            end if
            mixing_type = no_mixing
            gradT = gradr
            Y_face = gradT - gradL
            conv_vel = 0d0
            D = 0d0
            Gamma = 0d0
         end subroutine set_no_mixing
         

         subroutine do_compare_TDC_to_MLT
            ! for compare_TDC_to_MLT
            integer :: std_mixing_type
            type(auto_diff_real_star_order1) :: &
               std_gradT, std_gradr, std_conv_vel, std_D, std_Gamma, std_scale_height
            include 'formats'            
            std_mixing_type = mixing_type
            std_gradT = gradT
            std_gradr = gradr
            std_conv_vel = conv_vel
            std_D = D
            std_Gamma = Gamma
            std_scale_height = scale_height
            if (.false.) then
               write(*,2) 'std_mixing_type', std_mixing_type
               write(*,2) 'grada', k, grada%val
               write(*,2) 'std_gradT', k, std_gradT%val
               write(*,2) 'std_gradr', k, std_gradr%val
               write(*,2) 'std_conv_vel', k, std_conv_vel%val
               write(*,2) 'std_Gamma', k, std_Gamma%val
               write(*,2) 'mlt_vc_old', k, s% mlt_vc_old(k)
               write(*,2) 'std_Y_face', k, std_gradT%val - grada%val
            end if
            call set_TDC
            if (std_mixing_type /= mixing_type .or. &
                abs(std_gradT%val - gradT%val) > 1d-4) then
               write(*,4) 'mix type gradT', k, std_mixing_type, mixing_type, &
                  std_gradT%val - gradT%val
               stop 'do_compare_TDC_to_MLT'
            end if
         end subroutine do_compare_TDC_to_MLT


         subroutine set_TDC
            include 'formats'
            if (k == 1) then
               call set_no_mixing('set_TDC')
            else
               gradT = (wrap_lnT_m1(s,k) - wrap_lnT_00(s,k)) / &
                       (wrap_lnPeos_m1(s,k) - wrap_lnPeos_00(s,k))
               Y_guess = gradT - grada ! use actual gradT to guess Y_face
               call get_TDC_solution(s, k, &
                  mixing_length_alpha, cgrav, m, Y_guess, &
                  mixing_type, L, r, P, T, rho, Cp, opacity, &
                  scale_height, grada, conv_vel, Y_face, ierr)
               if (ierr /= 0) then
                  write(*,2) 'get_TDC_solution failed in set_TDC', k
                  stop 'get_TDC_solution failed in set_TDC'
               end if
            end if
            gradT = Y_face + grada
         end subroutine set_TDC
         
         
         subroutine set_MLT
            real(dp) :: ff1, ff2, ff3, ff4
            type(auto_diff_real_star_order1) :: &
               Q, omega, a0, ff4_omega2_plus_1, A_1, A_2, &
               A_numerator, A_denom, A, Bcubed, delta, Zeta, &
               f, f0, f1, f2, radiative_conductivity, convective_conductivity
            include 'formats' 
            Q = chiT/chiRho ! 'Q' param  C&G 14.24
            if (MLT_option == 'Cox') then ! this assumes optically thick
               a0 = 9d0/4d0
               convective_conductivity = &
                  Cp*grav*pow2(Lambda)*rho*(sqrt(Q*rho/(2d0*P)))/9d0 ! erg / (K cm sec)
               radiative_conductivity = &
                  (4d0/3d0*crad*clight)*pow3(T)/(opacity*rho) ! erg / (K cm sec)
               A = convective_conductivity / radiative_conductivity !  unitless.
            else
               select case(trim(MLT_option))
               case ('Henyey')
                  ff1=1.0d0/s% Henyey_MLT_nu_param
                  ff2=0.5d0 
                  ff3=8.0d0/s% Henyey_MLT_y_param
                  ff4=1.0d0/s% Henyey_MLT_y_param
               case ('ML1')
                  ff1=0.125d0 
                  ff2=0.5d0 
                  ff3=24.0d0
                  ff4=0.0d0
               case ('ML2')
                  ff1=1.0d0
                  ff2=2.0d0
                  ff3=16.0d0
                  ff4=0.0d0
               case ('Mihalas')
                  ff1=0.125d0 
                  ff2=0.5d0 
                  ff3=16.0d0
                  ff4=2.0d0
               case default
                  write(*,'(3a)') 'Error: ',trim(MLT_option), &
                     ' is not an allowed MLT version for convection'
                  write(*,*)
                  return
               end select
               omega = Lambda*rho*opacity
               ff4_omega2_plus_1 = ff4/pow2(omega) + 1d0
               a0 = (3d0/16d0)*ff2*ff3/ff4_omega2_plus_1
               A_1 = 4d0*Cp*sqrt(ff1*P*Q*rho)
               A_2 = mixing_length_alpha*omega*ff4_omega2_plus_1
               A_numerator = A_1*A_2
               A_denom = ff3*crad*clight*pow3(T)
               A = A_numerator/A_denom   
            end if  
            ! 'B' param  C&G 14.81
            Bcubed = (pow2(A)/a0)*(gradr - gradL)         
            ! now solve cubic equation for convective efficiency, Gamma
            ! a0*Gamma^3 + Gamma^2 + Gamma - a0*Bcubed == 0   C&G 14.82, 
            ! leave it to Mathematica to find an expression for the root we want      
            delta = a0*Bcubed               
            f = -2d0 + 9d0*a0 + 27d0*a0*a0*delta
            if (f > 1d100) then
               f0 = f
            else
               f0 = pow2(f) + 4d0*(-1d0 + 3d0*a0)*(-1d0 + 3d0*a0)*(-1d0 + 3d0*a0)
               if (f0 <= 0d0) then
                  f0 = f
               else
                  f0 = sqrt(f0)         
               end if
            end if
            f1 = -2d0 + 9d0*a0 + 27d0*a0*a0*delta + f0  
            if (f1 < 0) return
            f1 = pow(f1,one_third)     
            f2 = 2d0*two_13*(1d0 - 3d0*a0) / f1       
            Gamma = (four_13*f1 + f2 - 2d0) / (6d0*a0)
            if (Gamma < 0) return
            ! average convection velocity   C&G 14.86b
            conv_vel = mixing_length_alpha*sqrt(Q*P/(8d0*rho))*Gamma / A
            D = conv_vel*Lambda/3d0     ! diffusion coefficient [cm^2/sec]
            !Zeta = pow3(Gamma)/Bcubed  ! C&G 14.80     
            Zeta = exp(3d0*log(Gamma) - log(Bcubed)) ! write it this way to avoid overflow problems
            ! Zeta must be >= 0 and <= 1
            if (is_bad(Zeta%val) .or. Zeta < 0d0) then
               Zeta = 0d0
            else if (Zeta > 1d0) then
               Zeta = 1d0
            end if            
            gradT = (1d0 - Zeta)*gradr + Zeta*grada ! C&G 14.79      
            Y_face = gradT - grada
            mixing_type = convective_mixing
            if (k > 0) then
               s% xtra1_array(k) = conv_vel%val
               s% xtra2_array(k) = gradT%val
            end if

            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'set_MLT Zeta gradr grada gradT dgradT_dlnd', k, &
                  Zeta%val, gradr%val, grada%val, gradT%val
               write(*,2) 'set_MLT d_dlnd_00 Zeta gradr grada gradT', k, &
                  Zeta%d1Array(i_lnd_00), gradr%d1Array(i_lnd_00), grada%d1Array(i_lnd_00), gradT%d1Array(i_lnd_00)
            end if

         end subroutine set_MLT   


         subroutine set_semiconvection ! Langer 1983 & 1985
            type(auto_diff_real_star_order1) :: bc, LG, &
               radiative_conductivity, a0, a1, a2, a3, a4, a5, a6, a, &
               b1, b2, b3, b4, b5, b6, b7, b, div, bsq    
            real(dp) :: alpha
            include 'formats'
            radiative_conductivity = &
               (4d0/3d0*crad*clight)*pow3(T)/(opacity*rho) ! erg / (K cm sec)
            D = alpha_semiconvection*radiative_conductivity/(6d0*Cp*rho) &
                  *(gradr - grada)/(gradL - gradr)
            if (D%val <= 0) return
            
            if (s% semiconvection_option == 'Langer_85 mixing; gradT = gradr') then
               gradT = gradr
               Y_face = gradT - grada
               conv_vel = 3d0*D/Lambda             
               mixing_type = semiconvective_mixing
               return
            end if
            
            if (s% semiconvection_option /= 'Langer_85') then
               write(*,*) 'MLT: unknown values for semiconvection_option ' // &
                  trim(s% semiconvection_option)
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

            ! here's what Mathematica tells us to do
            alpha = min(1d0, alpha_semiconvection)
            bc = 32d0 - 24d0*beta - beta*beta            
            LG = (16d0*pi*clight*m*cgrav*Pr)/(P*opacity)            
            a0 = alpha*gradL_composition_term*LG            
            a1 = -2d0*bc*L            
            a2 = 2d0*alpha*bc*grada*LG            
            a3 = -2d0*bc*gradL*LG            
            a4 = 32d0*a0            
            a5 = -36d0*beta*a0            
            a6 = 9d0*beta*beta*a0            
            a = a1 + a2 + a3 + a4 + a5 + a6                           
            b1 = 32d0 - 36d0*beta + 9d0*beta*beta            
            b2 = b1*a0            
            b3 = -2d0*gradL*L + alpha*grada*grada*LG            
            b4 = (-alpha*gradA + gradL)*LG            
            b5 = -b2 + 2d0*bc*(L + b4)            
            b6 = b2*grada + bc*b3            
            b7 = -4d0*(alpha - 2d0)*bc*LG*b6            
            b = b7 + b5*b5            
            div = 2d0*(alpha - 2d0)*bc*LG
            bsq = sqrt(b)
            gradT = (a + bsq)/div
            Y_face = gradT - grada
            conv_vel = 3d0*D/Lambda             
            mixing_type = semiconvective_mixing
         end subroutine set_semiconvection
         
         
         subroutine set_thermohaline
            real(dp), parameter :: min_D_th = 1d-3
            real(dp) :: D_thrm
            call get_D_thermohaline(s, &
               grada%val, gradr%val, T%val, opacity%val, rho%val, &
               Cp%val, gradL_composition_term, &
               iso, XH1, thermohaline_coeff, D_thrm, ierr)
            if (D_thrm < min_D_th) return
            D = D_thrm
            gradT = gradr
            Y_face = gradT - grada
            conv_vel = 3d0*D/Lambda
            mixing_type = thermohaline_mixing 
         end subroutine set_thermohaline


      end subroutine Get_results_newer
      
      
      subroutine get_TDC_solution(s, k, &
            mixing_length_alpha, cgrav, m, Y_guess, &
            mixing_type, L, r, P, T, rho, Cp, kap, Hp, grada, cv, Y_face, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha, cgrav, m
         integer, intent(out) :: mixing_type
         type(auto_diff_real_star_order1), intent(in) :: &
            L, r, P, T, rho, Cp, kap, Hp, grada, Y_guess
         type(auto_diff_real_star_order1),intent(out) :: cv, Y_face
         integer, intent(out) :: ierr
         
         type(auto_diff_real_star_order1) :: A0, c0, L0
         type(auto_diff_real_tdc) :: Af, Y, Z, Q, Z_new, dQdZ
         real(dp) ::  tolerance, gradT, Lr, Lc, scale
         integer :: iter, max_iter
         logical :: converged, Y_is_positive
         include 'formats'
         ierr = 0
         if (mixing_length_alpha == 0d0 .or. k <= 1 .or. s% dt <= 0d0) then
            stop 'bad call to TDC get_TDC_solution'
         end if         
         c0 = mixing_length_alpha*rho*T*Cp*4d0*pi*pow2(r)
         L0 = (16d0*pi*crad*clight/3d0)*cgrav*m*pow4(T)/(P*kap) ! assumes QHSE for dP/dm
         if (s% okay_to_set_mlt_vc) then ! is also ok to use mlt_vc_old   
            A0 = s% mlt_vc_old(k)/sqrt_2_div_3
         else
            A0 = 0d0
         end if
         ! Newton's method to find solution Y
         Y = convert(Y_guess)
         if (Y == 0d0) Y = 1d-30
         Y%d1val1 = Y%val ! Fill in starting dY/dZ. Using Y = \pm exp(Z) we find dY/dZ = Y.
         Y_is_positive = (Y > 0d0)
         tolerance = 1d-8 ! ??
         max_iter = 20 ! ??
         converged = .false.
         scale = max(abs(s% L_start(k)), 1d-3*maxval(s% L_start(1:s% nz)))
         do iter = 1, max_iter
            call compute_Q(s, k, mixing_length_alpha, &
               Y, c0, L, L0, A0, T, rho, Cp, kap, Hp, grada, Q, Af)
            if (is_bad(Q%val)) exit
            if (abs(Q%val)/scale <= tolerance) then
               converged = .true.
               exit
            end if
            dQdZ = differentiate_1(Q)
            if (is_bad(dQdZ%val) .or. abs(dQdZ%val) < 1d-99) then
               exit
            end if
            Z_new = Z - Q / dQdZ
            Z_new%d1val1 = 1d0            
            Z = Z_new
            if (Y_is_positive) then
               Y = exp(Z)
            else
               Y = -exp(Z)
            end if
         end do
         if (.not. converged) then
            write(*,4) 'get_TDC_solution failed to converge', s% model_number, k, iter
            write(*,2) 'Q', k, Q%val
            write(*,2) 'scale', k, scale
            write(*,2) 'Q/scale', k, Q%val/scale
            write(*,2) 'tolerance', k, tolerance
            write(*,2) 'dQdZ', k, dQdZ%val
            write(*,2) 'Y', k, Y%val
            write(*,2) 'exp(Z)', k, exp(Z%val)
            write(*,2) 'Z', k, Z%val
            write(*,2) 'Af', k, Af%val
            write(*,2) 'A0', k, A0%val
            write(*,2) 'c0', k, c0%val
            write(*,2) 'L0', k, L0%val
            write(*,*)
            stop 'get_TDC_solution failed to converge'
         end if
         cv = sqrt_2_div_3*unconvert(Af)   
         Y_face = unconvert(Y)
         gradT = Y_face%val + grada%val
         Lr = L0%val*gradT
         Lc = L%val - Lr
         if (Lc > L%val*1d-4) then
            mixing_type = convective_mixing
         else
            mixing_type = no_mixing
         end if
      end subroutine get_TDC_solution
            
            
      subroutine compute_Q(s, k, mixing_length_alpha, &
            Y, c0_in, L_in, L0_in, A0, T, rho, Cp, kap, Hp, grada_in, Q, Af)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha
         type(auto_diff_real_star_order1), intent(in) :: &
            c0_in, L_in, L0_in, A0, T, rho, Cp, kap, Hp, grada_in
         type(auto_diff_real_tdc), intent(in) :: Y
         type(auto_diff_real_tdc), intent(out) :: Q, Af
         type(auto_diff_real_tdc) :: xi0, xi1, xi2, c0, L0, L, grada
         call eval_xis(s, k, mixing_length_alpha, &
            Y, T, rho, Cp, kap, Hp, grada_in, xi0, xi1, xi2) 
         Af = eval_Af(s, k, A0, xi0, xi1, xi2)
         L = convert(L_in)
         L0 = convert(L0_in)
         grada = convert(grada_in)
         c0 = convert(c0_in)
         Q = (L - L0*grada) - (L0 + c0*Af)*Y
      end subroutine compute_Q


      subroutine eval_xis(s, k, mixing_length_alpha, &
            Y, T, rho, Cp, kap, Hp, grada, xi0, xi1, xi2) 
         ! eval_xis sets up Y with partial wrt Z
         ! so results come back with partials wrt Z
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha
         type(auto_diff_real_star_order1), intent(in) :: T, rho, Cp, kap, Hp, grada
         type(auto_diff_real_tdc), intent(in) :: Y
         type(auto_diff_real_tdc), intent(out) :: xi0, xi1, xi2
         type(auto_diff_real_tdc) :: S0, D0, DR0
         type(auto_diff_real_star_order1) :: gammar_div_alfa
         real(dp), parameter :: x_ALFAS = (1.d0/2.d0)*sqrt_2_div_3
         real(dp), parameter :: x_CEDE  = (8.d0/3.d0)*sqrt_2_div_3
         real(dp), parameter :: x_GAMMAR = 2.d0*sqrt(3.d0)
         include 'formats'
         S0 = convert(x_ALFAS*mixing_length_alpha*Cp*T*grada/Hp)
         S0 = S0*Y
         D0 = convert(s% alpha_TDC_DAMP*x_CEDE/(mixing_length_alpha*Hp))
         if (s% alpha_TDC_DAMPR == 0d0) then
            DR0 = 0d0
         else
            gammar_div_alfa = s% alpha_TDC_DAMPR*x_GAMMAR/(mixing_length_alpha*Hp)
            DR0 = convert(4d0*boltz_sigma*pow2(gammar_div_alfa)*pow3(T)/(pow2(rho)*Cp*kap))
         end if
         xi0 = S0
         xi1 = -DR0
         xi2 = -D0         
      end subroutine eval_xis
      
      
      function eval_Af(s, k, A0_in, xi0, xi1, xi2) result(Af)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: A0_in
         type(auto_diff_real_tdc), intent(in) :: xi0, xi1, xi2
         type(auto_diff_real_tdc) :: Af ! output
         type(auto_diff_real_tdc) :: J2, J, Jt, Jt4, num, den, y_for_atan, root, A0        
         include 'formats'
         J2 = pow2(xi1) - 4d0 * xi0 * xi2
         A0 = convert(A0_in)
         if (J2 > 0d0) then ! Hyperbolic branch
            J = sqrt(J2)
            Jt = s%dt * J
            Jt4 = 0.25d0 * Jt
            num = tanh(Jt4) * (2d0 * xi0 + A0 * xi1) + A0 * J
            den = tanh(Jt4) * (xi1 + 2d0 * A0 * xi2) - J
            Af = num / den
            if (Af < 0d0) then
               Af = -Af
            end if
         else if (J2 < 0d0) then ! Trigonometric branch
            J = sqrt(-J2)
            Jt = s%dt * J
            ! This branch contains decaying solutions that reach A = 0, at which point
            ! they switch onto the 'zero' branch. So we have to calculate the position of
            ! the first root to check it against dt.
            y_for_atan = xi1 + 2d0 * A0 * xi2
            ! We had a choice above to pick which of +-I to use in switching branches.
            ! That choice has to be consistent with a decaying solution, which we check now.
            if (y_for_atan > 0d0) then
               J = -J
            end if
            root = two_var_pos_atan(J, y_for_atan) - two_var_pos_atan(J, xi1)
            if (0.25d0 * Jt < root) then
               num = -xi1 + J * tan(0.25d0 * Jt + atan(y_for_atan / J)) 
               den = 2d0 * xi2
               Af = num / den
            else
               Af = 0d0
            end if
         else ! if (J2 == 0d0) then         
            Af = A0            
         end if         
         if (k > 0) then ! save for plots
            s% SOURCE(k) = xi0%val*Af%val
            s% DAMP(k) = -xi1%val*Af%val
            s% DAMPR(k) = -xi2%val*pow2(Af%val)
            s% COUPL(k) = s% SOURCE(k) - s% DAMP(k) - s% DAMPR(k)
         end if
      end function eval_Af
      
      
      !> Returns the smallest positive z such that tan(z) = y/x
      type(auto_diff_real_tdc) function two_var_pos_atan(x,y) result(z)
         type(auto_diff_real_tdc), intent(in) :: x,y
         type(auto_diff_real_tdc) :: x1, y1
         x1 = abs(x) + 1d-50
         y1 = abs(y) + 1d-50
         z = atan(y1/x1)
         if (z < 0d0) then
            z = z + pi
         end if
      end function two_var_pos_atan
      
      
      function convert(K_in) result(K)
         type(auto_diff_real_star_order1), intent(in) :: K_in
         type(auto_diff_real_tdc) :: K
         K%val = K_in%val
         K%d1Array(1:auto_diff_star_num_vars) = K_in%d1Array(1:auto_diff_star_num_vars)
         K%d1val1 = 0d0
         K%d1val1_d1Array(1:auto_diff_star_num_vars) = 0d0
      end function convert
      
      
      function unconvert(K_in) result(K)
         type(auto_diff_real_tdc), intent(in) :: K_in
         type(auto_diff_real_star_order1) :: K
         K%val = K_in%val
         K%d1Array(1:auto_diff_star_num_vars) = K_in%d1Array(1:auto_diff_star_num_vars)
      end function unconvert
      

      function get_Eq_0_face(s,k) result(Eq0)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Eq0
         type(auto_diff_real_star_order1) :: &
            v_div_r_m1, v_div_r_00, v_div_r_p1, d_v_div_r_dm_face, &
            rho2_face, r6_face, Hp_face, Chi_div_w_face
         real(dp) :: f
         include 'formats'
         if (k == 1 .or. k == s% nz) then
            Eq0 = 0d0
            return
         end if
         v_div_r_m1 = wrap_v_m1(s,k)/wrap_r_m1(s,k)
         if (.true.) then ! one-sided to avoid partials wrt p1
            v_div_r_00 = wrap_v_00(s,k)/wrap_r_00(s,k)
            d_v_div_r_dm_face = (v_div_r_m1 - v_div_r_00)/s% dm(k-1)
         else
            v_div_r_p1 = wrap_v_p1(s,k)/wrap_r_p1(s,k)
            d_v_div_r_dm_face = (v_div_r_m1 - v_div_r_p1)/(s% dm(k-1) + s% dm(k))
         end if
         f = (16d0/3d0)*pi*s% alpha_TDC_eddy_viscosity*s% mixing_length_alpha
         rho2_face = pow2(get_rho_face(s,k))
         r6_face = pow6(wrap_r_00(s,k))
         Hp_face = get_Hp_face(s,k)
         Chi_div_w_face = f*rho2_face*r6_face*Hp_face*d_v_div_r_dm_face 
         Eq0 = 4d0*pi*Chi_div_w_face*d_v_div_r_dm_face
      end function get_Eq_0_face


      function get_Pturb_0_face(s,k) result(Pturb_0)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Pturb_0
         real(dp), parameter :: x_ALFAP = 2.d0/3.d0
         Pturb_0 = s% alpha_TDC_turbulent_pressure*x_ALFAP*get_rho_face(s,k)
      end function get_Pturb_0_face


      function get_dVdt_face(s,k) result(dVdt)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: dVdt
         real(dp) :: alfa, beta
         type(auto_diff_real_star_order1) :: dV_00, dV_m1
         if (k == 1) then
            dVdt = 0
            return
         end if
         call get_face_weights(s, k, alfa, beta)
         ! dV = 1d0/d_00 - 1d0/s% rho_start(k) ! = (rho_start - rho)/(rho*rho_start)
         ! dlnd = wrap_dxh_lnd = lnd(k) - lnd_start(k) from solver without subtraction
         ! expm1(dlnd) = (rho - rho_start)/rho_start
         dV_00 = -expm1(wrap_dxh_lnd(s,k))/wrap_d_00(s,k)
         dV_m1 = -expm1(wrap_dxh_lnd(s,k-1))/wrap_d_m1(s,k)
         dVdt = (alfa*dV_00 + beta*dV_m1)/s% dt
      end function get_dVdt_face
      
      
      function get_d_v_div_r_cell(s,k) result(d_v_div_r) ! s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: d_v_div_r
         type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
         include 'formats'
         v_00 = wrap_v_00(s,k)
         v_p1 = wrap_v_p1(s,k)
         r_00 = wrap_r_00(s,k)
         r_p1 = wrap_r_p1(s,k)
         if (r_p1%val == 0d0) r_p1 = 1d0
         d_v_div_r = v_00/r_00 - v_p1/r_p1 ! units s^-1
      end function get_d_v_div_r_cell
      
      
      function get_Chi_cell(s, k) result(Chi_cell) 
         ! eddy viscosity energy (Kuhfuss 1986) [erg]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Chi_cell
         type(auto_diff_real_star_order1) :: &
            rho2, r6_cell, d_v_div_r, Hp_cell, w_cell
         real(dp) :: f, ALFAM_ALFA
         integer :: j
         include 'formats'
         ALFAM_ALFA = s% alpha_TDC_eddy_viscosity*s% mixing_length_alpha
         if (ALFAM_ALFA == 0d0) then
            Chi_cell = 0d0
         else
            f = (16d0/3d0)*pi*ALFAM_ALFA/s% dm(k)  
            rho2 = pow2(wrap_d_00(s,k))
            r6_cell = 0.5d0*(pow6(wrap_r_00(s,k)) + pow6(wrap_r_p1(s,k)))
            d_v_div_r = get_d_v_div_r_cell(s, k)
            Hp_cell = get_Hp_cell(s,k)
            Chi_cell = wrap_w_00(s,k)*f*rho2*r6_cell*d_v_div_r*Hp_cell
            ! use solver variable w_cell for consistency of Eq and Uq and 3 point stencil
         end if
         s% Chi(k) = Chi_cell%val
      end function get_Chi_cell


      function get_Eq_cell(s,k) result(Eq) ! for energy equation
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Eq
         type(auto_diff_real_star_order1) :: Chi_cell, d_v_div_r, Eq_div_w
         integer :: ierr
         logical :: test_partials
         include 'formats'
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         if (k == 1) then
            Eq = 0d0
         else
            Chi_cell = get_Chi_cell(s,k)
            d_v_div_r = get_d_v_div_r_cell(s,k)
            Eq = 4d0*pi*Chi_cell*d_v_div_r/s% dm(k) ! erg s^-1 g^-1
         end if
         s% Eq(k) = Eq%val
         if (test_partials) then
            s% solver_test_partials_val = Eq%val
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = Eq%d1Array(i_lnR_00)
            write(*,4) 'get_Eq_cell', s% solver_test_partials_var
         end if      
      end function get_Eq_cell


      function get_Uq_face(s,k) result(Uq) ! for momentum equation
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Uq
         type(auto_diff_real_star_order1) :: dChi_dm_bar, Chi_00, Chi_out
         integer :: ierr
         logical :: test_partials
         include 'formats'
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         if (k == 1) then
            Uq = 0d0
         else
            Chi_00 = get_Chi_cell(s,k)
            if (k > 1) then
               Chi_out = shift_m1(get_Chi_cell(s,k-1)) ! , 'get_Uq_face')
            else
               Chi_out = 0d0
            end if
            dChi_dm_bar = (Chi_out - Chi_00)/s% dm_bar(k)
            Uq = 4d0*pi*dChi_dm_bar/wrap_r_00(s,k)
         end if
         ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2, acceleration
         s% Uq(k) = Uq%val      
         if (test_partials) then
            s% solver_test_partials_val = s% Uq(k) 
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = Uq%d1Array(i_lnd_00)
            write(*,4) 'get_Uq_face', s% solver_test_partials_var
         end if      
      end function get_Uq_face


      function get_Hp_face(s, k) result(Hp_face)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Hp_face
         type(auto_diff_real_star_order1) :: &
            r_00, r2, Peos_00, d_00, Peos_m1, d_m1, Peos_div_rho, &
            rho_face, area, dlnPeos, Peos_face, alt_Hp_face, A
         real(dp) :: alfa, beta, cgrav, m
         include 'formats'
         if (k > s% nz) then
            Hp_face = 0d0 ! not used
            Hp_face%val = 1d0
            s% Hp_face(k) = Hp_face%val
            return
         end if
         cgrav = s% cgrav(k)
         m = s% m(k)
         r_00 = wrap_r_00(s,k)
         r2 = pow2(r_00)
         d_00 = wrap_d_00(s,k)
         Peos_00 = wrap_Peos_00(s,k)
         if (k == 1) then
            Peos_div_rho = Peos_00/d_00
            Hp_face = r2*Peos_div_rho/(cgrav*m)
         else
            d_m1 = wrap_d_m1(s,k)
            Peos_m1 = wrap_Peos_m1(s,k)
            call get_face_weights(s, k, alfa, beta)
            Peos_div_rho = alfa*Peos_00/d_00 + beta*Peos_m1/d_m1
            Hp_face = r2*Peos_div_rho/(cgrav*m)
            if (s% alt_scale_height_flag) then
               ! consider sound speed*hydro time scale as an alternative scale height
               rho_face = alfa*d_00 + beta*d_m1
               Peos_face = alfa*Peos_00 + beta*Peos_m1
               alt_Hp_face = sqrt(Peos_face/cgrav)/rho_face
               if (alt_Hp_face%val < Hp_face%val) then ! blend
                  A = pow2(alt_Hp_face/Hp_face) ! 0 <= A%val < 1
                  Hp_face = A*Hp_face + (1d0 - A)*alt_Hp_face
               end if
            end if
         end if
         s% Hp_face(k) = Hp_face%val
      end function get_Hp_face

      
      function get_Hp_cell(s, k) result(Hp_cell) ! cm
         ! instead of 0.5d0*(Hp_face(k) + Hp_face(k+1)) to keep block tridiagonal
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Hp_cell
         type(auto_diff_real_star_order1) :: r_mid, r_00, r_p1, &
            Peos_00, d_00, alt_Hp_cell, alfa
         real(dp) :: cgrav_00, cgrav_p1, cgrav_mid, m_00, m_p1, m_mid
         include 'formats'
         r_00 = wrap_r_00(s,k)
         cgrav_00 = s% cgrav(k)
         m_00 = s% m(k)
         d_00 = wrap_d_00(s,k)
         Peos_00 = wrap_Peos_00(s,k)
         r_p1 = wrap_r_p1(s,k)
         if (k < s% nz) then
            cgrav_p1 = s% cgrav(k+1)
            m_p1 = s% m(k+1)
         else
            cgrav_p1 = s% cgrav(k)
            m_p1 = s% m_center
         end if
         cgrav_mid = 0.5d0*(cgrav_00 + cgrav_p1)
         m_mid = 0.5d0*(m_00 + m_p1)
         r_mid = 0.5d0*(r_00 + r_p1)
         Hp_cell = pow2(r_mid)*Peos_00 / (d_00*cgrav_mid*m_mid)
         if (s% alt_scale_height_flag) then
            ! consider sound speed*hydro time scale as an alternative scale height
            alt_Hp_cell = sqrt(Peos_00/cgrav_mid)/d_00
            if (alt_Hp_cell%val < Hp_cell%val) then ! blend
               alfa = pow2(alt_Hp_cell/Hp_cell) ! 0 <= alfa%val < 1
               Hp_cell = alfa*Hp_cell + (1d0 - alfa)*alt_Hp_cell
            end if
         end if
      end function get_Hp_cell


!------------------------------


      subroutine get_D_thermohaline(s, &
            grada, gradr, T, opacity, rho, Cp, gradL_composition_term, &
            iso, XH1, thermohaline_coeff, D_thrm, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: &
            grada, gradr, T, opacity, rho, Cp, gradL_composition_term, XH1, &
            thermohaline_coeff
         integer, intent(in) :: iso      
         real(dp), intent(out) :: D_thrm
         integer, intent(out) :: ierr
         real(dp) :: dgrad, K_therm, K_T, K_mu, nu, R0, Pr, tau, r_th            
         include 'formats'     
         dgrad = max(1d-40, grada - gradr) ! positive since Schwarzschild stable               
         K_therm = 4d0*crad*clight*pow3(T)/(3d0*opacity*rho) ! thermal conductivity
         if (s% thermohaline_option == 'Kippenhahn') then
            ! Kippenhahn, R., Ruschenplatt, G., & Thomas, H.-C. 1980, A&A, 91, 175
            D_thrm = -3d0*K_therm/(2*rho*cp)*gradL_composition_term/dgrad
         else if (s% thermohaline_option == 'Brown_Garaud_Stellmach_13' .or. &
                  s% thermohaline_option == 'Traxler_Garaud_Stellmach_11') then
            call get_diff_coeffs(s, &
               K_therm, Cp, rho, T, opacity, iso, XH1, K_T, K_mu, nu)
            R0 = (gradr - grada)/gradL_composition_term
            Pr = nu/K_T
            tau = K_mu/K_T
            r_th = (R0 - 1d0)/(1d0/tau - 1d0)
            if (r_th >= 1d0) then ! stable if R0 >= 1/tau
               D_thrm = 0d0
            else if (s% thermohaline_option == 'Traxler_Garaud_Stellmach_11') then 
               ! Traxler, Garaud, & Stellmach, ApJ Letters, 728:L29 (2011).
               ! also see Denissenkov. ApJ 723:563–579, 2010.
               D_thrm = 101d0*sqrt(K_mu*nu)*exp(-3.6d0*r_th)*pow(1d0 - r_th,1.1d0) ! eqn 24
            else ! if (s% thermohaline_option == 'Brown_Garaud_Stellmach_13') then
               D_thrm = K_mu*(Numu(R0,r_th,pr,tau) - 1d0)
               ! evbauer 07/18: changed from K_mu*Numu(R0,r_th,pr,tau), Pascale signed off
            endif
         else
            D_thrm = 0
            ierr = -1
            write(*,*) 'unknown for MLT thermohaline_option' // trim(s% thermohaline_option)
         end if
         D_thrm = thermohaline_coeff*D_thrm
      end subroutine get_D_thermohaline


      subroutine get_diff_coeffs(s, &
            K_therm, Cp, rho, T, opacity, iso, XH1, kt, kmu, vis)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         real(dp), intent(in) :: K_therm, Cp, rho, T, opacity, XH1
         integer, intent(in) :: iso      
         real(dp), intent(out) :: kt, kmu, vis
         real(dp) :: loglambdah, loglambdacx, loglambdacy, ccx, ccy, qe4
         real(dp) :: Bcoeff, chemA, chemZ, acx, acy, nu_mol, nu_rad      
         real(dp), parameter :: sqrt5 = sqrt(5d0)           
         kt = K_therm/(Cp*rho)       ! thermal diffusivity (assumes radiatively dominated)
         qe4=pow4(qe)
      
         ! Log Lambda for pure H (equation 10 from Proffitt Michaud 93)
         loglambdah = -19.26d0 - 0.5d0*log(rho) + 1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+XH1)) 
         nu_rad = 4d0*crad*pow4(T)/(15d0*clight*opacity*pow2(rho)) ! radiative viscosity
         nu_mol = 0.406d0*sqrt(amu)*pow(boltzm*T,2.5d0)/(qe4*loglambdah*rho) 
         ! From Spitzer "Physics of Fully Ionized Gases equation 5-54
         ! Assumes pure H. Still trying to work out what it would be for a mixture. 
         vis = nu_mol + nu_rad   ! total viscosity

         ! The following is from Proffitt & Michaud, 1993.
         ! Their constant B (equation 15)
         Bcoeff = (15.d0/16.d0)*sqrt(2.d0*amu/(5*pi))*pow(boltzm,2.5d0)/qe4
         ! Extract what species drives the thermohaline concvection
         chemA = chem_isos%Z_plus_N(iso)
         chemZ = chem_isos%Z(iso)

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
               (XH1*sqrt(acx)*ccx + (1-XH1)*sqrt(acy)*ccy)

         else
            ! Log Lambda for H-He mixture (equation 10)
            loglambdah = -19.26d0 - log(2d0) - 0.5d0*log(rho) + &
               1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+XH1)) 
            ! Calculation of C_ij coeffs (equation 12)
            ccy = log(exp(1.2d0*loglambdah)+1d0)/1.2d0
            ! My formula (see notes) based on Proffitt and Michaud 1993
            kmu = (Bcoeff*pow(T,2.5d0)/(rho*ccy))*(3+XH1)/((1+XH1)*(3+5*XH1)*(0.7d0+0.3d0*XH1))
         
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


      end module mlt_get_results_newer
