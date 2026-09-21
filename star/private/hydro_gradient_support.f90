! ***********************************************************************
!
!   Copyright (C) 2018-2019  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module hydro_gradient_support

      use star_private_def
      use const_def, only: dp, pi4, crad, clight
      use auto_diff
      use auto_diff_support
      use reconstructed_face_support, only: get_reconstructed_face_eos_kap_ad
      use utils_lib, only: is_bad, mesa_error

      implicit none
      private
      public :: expected_HSE_grav_term, eval_dlnPdm_qhse
      public :: get_RSP2_alfa_beta_face_weights, get_rsp2_face_eos
      public :: get_dPrad_dm_factors, get_rsp2_Lrad_coeff, get_rsp2_thermal_gradient

      contains

      ! Returns -G*m/r^2 with possible modifications for rotation. MESA 2, eqn 22.
      subroutine expected_HSE_grav_term(s, k, grav, area, ierr, use_time_centering)
         use star_utils, only: get_area_info_opt_time_center
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: area, grav
         integer, intent(out) :: ierr

         logical, intent(in), optional :: use_time_centering
         type(auto_diff_real_star_order1) :: inv_R2
         logical :: test_partials, time_centering

         include 'formats'
         ierr = 0

         time_centering = .true.
         if (present(use_time_centering)) time_centering = use_time_centering
         if (time_centering) then
            call get_area_info_opt_time_center(s, k, area, inv_R2, ierr)
            if (ierr /= 0) return
         else
            inv_R2 = 1d0/pow2(wrap_r_00(s,k))
            area = pi4*pow2(wrap_r_00(s,k))
         end if

         if (s% rotation_flag .and. s% use_gravity_rotation_correction) then
            grav = -s% cgrav(k)*s% m_grav(k)*inv_R2*s% fp_rot(k)
         else
            grav = -s% cgrav(k)*s% m_grav(k)*inv_R2
         end if

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = 0
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'expected_HSE_grav_term', s% solver_test_partials_var
         end if

      end subroutine expected_HSE_grav_term

      subroutine eval_dlnPdm_qhse(s, k, &  ! calculate the expected dlnPdm for HSE
            dlnPdm_qhse, Ppoint, ierr, use_time_centering)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dlnPdm_qhse, Ppoint
         integer, intent(out) :: ierr

         logical, intent(in), optional :: use_time_centering
         logical :: time_centering
         real(dp) :: alfa, P_theta
         type(auto_diff_real_star_order1) :: grav, area, P00, Pm1, mlt_Ptrb00, mlt_Ptrbm1, mlt_Ptrb_face
         type(auto_diff_real_star_order1) :: T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, opacity_face
         include 'formats'

         ierr = 0
         time_centering = .true.
         if (present(use_time_centering)) time_centering = use_time_centering

         ! basic eqn is dP/dm = -G m / (4 pi r^4)
         ! divide by Ppoint to make it unitless

         ! for rotation, multiply gravity by factor fp.  MESA 2, eqn 22.
         call expected_HSE_grav_term(s, k, grav, area, ierr, time_centering)
         if (ierr /= 0) return

         if (time_centering .and. s% using_velocity_time_centering .and. &
               s% include_P_in_velocity_time_centering) then
            P_theta = s% P_theta_for_velocity_time_centering
         else
            P_theta = 1d0
         end if

         if (s% use_face_reconstruction) then
            if (s% reconstructed_face_state_valid(k)) then
               rho_face = s% reconstructed_rho_face_ad(k)
               Ppoint = s% reconstructed_P_face_ad(k)
            else
               call get_reconstructed_face_eos_kap_ad( &
                  s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, opacity_face, ierr)
               if (ierr /= 0) return
               Ppoint = P_face
            end if
            if (P_theta /= 1d0) then
               Ppoint = P_theta*Ppoint + (1d0 - P_theta)*s% reconstructed_P_face_start(k)
            end if
            if (s% have_mlt_vc .and. s% okay_to_set_mlt_vc .and. s% include_mlt_Pturb_in_thermodynamic_gradients &
               .and. s% mlt_Pturb_factor > 0d0) then
               ! Keep the lagged convective velocity, but form the pressure term from the same
               ! face density used by the reconstructed face thermodynamic quantities.
               mlt_Ptrb_face = s% mlt_Pturb_factor*pow2(s% mlt_vc_old(k))*rho_face/3d0
               Ppoint = Ppoint + mlt_Ptrb_face
            end if
         else
            ! mlt_pturb in thermodynamic gradients does not currently support time centering because it is timelagged.
            ! replace mlt_vc check with s% mlt_vc_old(k) >0 check.
            if ((s% have_mlt_vc .and. s% okay_to_set_mlt_vc) .and. s% include_mlt_Pturb_in_thermodynamic_gradients &
               .and. s% mlt_Pturb_factor > 0d0) then
               if (k ==1) then
                  mlt_Ptrb00 = s% mlt_Pturb_factor*pow2(s% mlt_vc_old(k))*wrap_d_00(s,k)/3d0
                  mlt_Ptrbm1 = 0d0
               else
                  mlt_Ptrb00 = s% mlt_Pturb_factor*pow2(s% mlt_vc_old(k))*wrap_d_00(s,k)/3d0
                  mlt_Ptrbm1 = s% mlt_Pturb_factor*pow2(s% mlt_vc_old(k))*wrap_d_m1(s,k)/3d0
               end if
            else  ! no mlt_pturb
               mlt_Ptrb00 = 0d0
               mlt_Ptrbm1 = 0d0
            end if

            P00 = wrap_Peos_00(s,k)

            ! mlt Pturb doesn't support time centering yet.
            if (P_theta /= 1d0) P00 = P_theta*P00 + (1d0 - P_theta)*s% Peos_start(k)

            if (k == 1) then
               Pm1 = 0d0
               Ppoint = P00 + mlt_Ptrb00
            else
               Pm1 = wrap_Peos_m1(s,k)
               if (P_theta /= 1d0) Pm1 = P_theta*Pm1 + (1d0 - P_theta)*s% Peos_start(k-1)
               Pm1 = Pm1 + mlt_Ptrbm1  ! include mlt Ptrb in k-1
               P00 = P00 + mlt_Ptrb00  ! include mlt Ptrb in k
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
               Ppoint = alfa*P00 + (1d0-alfa)*Pm1
            end if
         end if

         dlnPdm_qhse = grav/(area*Ppoint)  ! note that expected_HSE_grav_term is negative

         if (is_bad(dlnPdm_qhse%val)) then
            ierr = -1
            s% retry_message = 'eval_dlnPdm_qhse: is_bad(dlnPdm_qhse)'
            if (s% report_ierr) then
!$OMP critical (hydro_vars_crit1)
               write(*,*) 'eval_dlnPdm_qhse: is_bad(dlnPdm_qhse)'
               stop
!$OMP end critical (hydro_vars_crit1)
            end if
            if (s% stop_for_bad_nums) then
               write(*,2) 'dlnPdm_qhse', k, dlnPdm_qhse
               call mesa_error(__FILE__,__LINE__,'eval_dlnPdm_qhse')
            end if
            return
         end if

      end subroutine eval_dlnPdm_qhse

      subroutine get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: alfa, beta
         ! face_value = alfa*cell_value(k) + beta*cell_value(k-1)
         if (k == 1) then
            alfa = 1d0
            beta = 0d0
            return
         end if
         if (s% RSP2_use_mass_interp_face_values) then
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1d0 - alfa
         else
            alfa = 0.5d0
            beta = 0.5d0
         end if
      end subroutine get_RSP2_alfa_beta_face_weights

      subroutine get_rsp2_face_eos( &
            s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face
         integer, intent(out) :: ierr
         real(dp) :: alfa, beta

         ierr = 0
         if (s% use_face_reconstruction .and. k > 1) then
            call get_reconstructed_face_eos_kap_ad( &
               s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
            return
         end if
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         T_face = alfa*wrap_T_00(s,k) + beta*wrap_T_m1(s,k)
         rho_face = alfa*wrap_d_00(s,k) + beta*wrap_d_m1(s,k)
         P_face = alfa*wrap_Peos_00(s,k) + beta*wrap_Peos_m1(s,k)
         Cp_face = alfa*wrap_Cp_00(s,k) + beta*wrap_Cp_m1(s,k)
         ChiRho_face = alfa*wrap_chiRho_00(s,k) + beta*wrap_chiRho_m1(s,k)
         ChiT_face = alfa*wrap_chiT_00(s,k) + beta*wrap_chiT_m1(s,k)
         grad_ad = alfa*wrap_grad_ad_00(s,k) + beta*wrap_grad_ad_m1(s,k)
         kap_face = alfa*wrap_kap_00(s,k) + beta*wrap_kap_m1(s,k)
      end subroutine get_rsp2_face_eos


      subroutine get_dPrad_dm_factors(s, k, opacity_face, kap_face, flxR, flxLambda, dm_bar)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: opacity_face
         type(auto_diff_real_star_order1), intent(out) :: kap_face, flxR, flxLambda
         real(dp), intent(out) :: dm_bar
         type(auto_diff_real_star_order1) :: area, T4_m1, T4_00

         dm_bar = s% dm_bar(k)
         ! Use ordinary cell-center spacing next to an excised inner boundary.
         if (s% R_center > 0d0 .and. k == s% nz) &
            dm_bar = 0.5d0*(s% dm(k - 1) + s% dm(k))
         kap_face = opacity_face
         if (kap_face%val < s% min_kap_for_dPrad_dm_eqn) kap_face = s% min_kap_for_dPrad_dm_eqn
         flxR = 0d0
         flxLambda = 1d0
         if (.not. s% use_flux_limiting_with_dPrad_dm_form) return
         area = pi4*pow2(wrap_r_00(s,k))
         T4_m1 = pow4(wrap_T_m1(s,k))
         T4_00 = pow4(wrap_T_00(s,k))
         flxR = area*abs(T4_m1 - T4_00)/dm_bar/ &
            (kap_face*0.5d0*(T4_m1 + T4_00))
         flxLambda = (6d0 + 3d0*flxR)/(6d0 + (3d0 + flxR)*flxR)
      end subroutine get_dPrad_dm_factors


      function get_rsp2_Lrad_coeff(s, k, ierr) result(Lrad_coeff)
         use tdc_hydro, only: get_TDC_Hp_face
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lrad_coeff, &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, Hp_face, area

         Lrad_coeff = 0d0
         call get_rsp2_face_eos( &
            s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
         if (ierr /= 0) return
         Hp_face = get_TDC_Hp_face(s, k, ierr)
         if (ierr /= 0) return
         if (Hp_face%val <= 0d0 .or. kap_face%val <= 0d0 .or. rho_face%val <= 0d0) then
            ierr = -1
            s% retry_message = 'invalid RSP2 radiative coefficient'
            return
         end if
         area = pi4*pow2(wrap_r_00(s,k))
         Lrad_coeff = area*(4d0*crad*clight/3d0)*pow4(T_face)/(kap_face*rho_face*Hp_face)
      end function get_rsp2_Lrad_coeff


      subroutine get_rsp2_thermal_gradient(s, k, grad_ad, gradL, entropy_gradient, ierr, use_time_centering)
         use star_utils, only: get_kap_face
         use turb_support, only: get_TDC_dynamical_gradL
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: grad_ad, gradL, entropy_gradient
         integer, intent(out) :: ierr
         logical, intent(in), optional :: use_time_centering
         type(auto_diff_real_star_order1) :: &
            T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, kap_face, &
            area, grav, area_qhse, dlnPdm, Ppoint, Tpoint, T_00, T_m1, P_00, P_m1, &
            pressure_jump, pressure_jump_reference, pressure_gradient_reference, &
            temperature_gradient_per_gradT, mass_to_radius, &
            opacity_row, kap_row, flxR, flxLambda, Lrad_coeff, Prad_face, Y_face, dynamical_gradL
         real(dp) :: alfa, dm_temperature, dm_momentum, composition_term
         logical :: dynamical, time_centering

         entropy_gradient = 0d0
         call get_rsp2_face_eos( &
            s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grad_ad, kap_face, ierr)
         if (ierr /= 0) return
         composition_term = 0d0
         if (s% use_Ledoux_criterion) composition_term = s% gradL_composition_term(k)
         gradL = grad_ad + composition_term
         if (k <= 1) return
         time_centering = .true.
         if (present(use_time_centering)) time_centering = use_time_centering
         dynamical = s% TDC_use_dynamical_gradL .and. (s% u_flag .or. s% v_flag)
         dm_temperature = 0.5d0*(s% dm(k - 1) + s% dm(k))
         area = pi4*pow2(wrap_r_00(s,k))
         T_00 = wrap_T_00(s,k)
         T_m1 = wrap_T_m1(s,k)
         P_00 = wrap_Peos_00(s,k)
         P_m1 = wrap_Peos_m1(s,k)
         pressure_jump = P_00 - P_m1
         Y_face = wrap_Y_00(s,k)
         if (T_face%val <= 0d0 .or. P_face%val <= 0d0 .or. rho_face%val <= 0d0 .or. &
               P_00%val <= 0d0 .or. P_m1%val <= 0d0 .or. dm_temperature <= 0d0) then
            ierr = -1
            s% retry_message = 'invalid RSP2 thermal-gradient state'
            return
         end if
         mass_to_radius = area*rho_face/dm_temperature

         if (s% constant_L) then
            ! There is no temperature row to eliminate in this path.
            entropy_gradient = Cp_face*mass_to_radius* &
               ((T_00 - T_m1)/T_face - grad_ad*pressure_jump/P_face)
            if (dynamical) then
               call get_TDC_dynamical_gradL(s, k, gradL, dynamical_gradL, ierr)
               if (ierr /= 0) return
               gradL = dynamical_gradL
            end if
         else if (s% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn) then
            ! gradT already uses the resolved pressure gradient in this row.
            entropy_gradient = Cp_face*mass_to_radius*log1p(pressure_jump/P_m1)*(Y_face + composition_term)
         else
            if (s% use_dPrad_dm_form_of_T_gradient_eqn) then
               opacity_row = kap_face
               if (.not. s% use_face_reconstruction) opacity_row = get_kap_face(s,k)
               call get_dPrad_dm_factors(s, k, opacity_row, kap_row, flxR, flxLambda, dm_temperature)
               mass_to_radius = area*rho_face/dm_temperature
               Lrad_coeff = get_rsp2_Lrad_coeff(s, k, ierr)
               if (ierr /= 0) return
               Prad_face = (crad/3d0)*pow4(T_face)
               temperature_gradient_per_gradT = &
                  rho_face*kap_row*Lrad_coeff/(4d0*clight*area*flxLambda*Prad_face)
            else
               call eval_dlnPdm_qhse(s, k, dlnPdm, Ppoint, ierr, time_centering)
               if (ierr /= 0) return
               alfa = s% dm(k - 1)/(s% dm(k - 1) + s% dm(k))
               Tpoint = alfa*T_00 + (1d0 - alfa)*T_m1
               temperature_gradient_per_gradT = -area*rho_face*dlnPdm*Tpoint/T_face
            end if
            if (temperature_gradient_per_gradT%val <= 0d0 .or. is_bad(temperature_gradient_per_gradT%val)) then
               ierr = -1
               s% retry_message = 'invalid RSP2 temperature-gradient coefficient'
               return
            end if
            if (dynamical) then
               pressure_jump_reference = pressure_jump
            else
               call expected_HSE_grav_term(s, k, grav, area_qhse, ierr, time_centering)
               if (ierr /= 0) return
               dm_momentum = 0.5d0*(s% dm(k - 1) + s% dm(k))
               if (s% use_mass_corrections) dm_momentum = 0.5d0* &
                  (s% dm(k - 1)*s% mass_correction(k - 1) + s% dm(k)*s% mass_correction(k))
               pressure_jump_reference = -grav*dm_momentum/area_qhse
            end if
            pressure_gradient_reference = mass_to_radius*pressure_jump_reference/P_face
            gradL = (grad_ad + composition_term)*pressure_gradient_reference/temperature_gradient_per_gradT
            ! Cancel the neutral contribution before adding the independent Y.
            entropy_gradient = Cp_face*(temperature_gradient_per_gradT*Y_face + &
               pressure_gradient_reference*composition_term)
            if (.not. dynamical) entropy_gradient = entropy_gradient + &
               Cp_face*grad_ad*mass_to_radius*(pressure_jump_reference - pressure_jump)/P_face
         end if
         if (is_bad(gradL%val) .or. is_bad(entropy_gradient%val)) then
            ierr = -1
            s% retry_message = 'bad RSP2 thermal gradient'
         end if
      end subroutine get_rsp2_thermal_gradient


      end module hydro_gradient_support
