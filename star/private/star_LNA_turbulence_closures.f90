! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
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

      module star_lna_turbulence_closures

      use star_private_def
      use const_def, only: boltz_sigma, clight, crad, dp, no_mixing, pi, &
         sqrt_2_div_3
      use math_lib, only: pow2, pow3, pow4
      use auto_diff
      use auto_diff_support, only: &
         shift_m1, shift_p1, wrap, wrap_Cp_00, wrap_Cp_m1, wrap_chiRho_00, &
         wrap_chiT_00, wrap_d_00, wrap_d_m1, wrap_etrb_00, wrap_Hp_00, &
         wrap_Hp_p1, wrap_kap_00, wrap_L_00, wrap_lnPeos_00, &
         wrap_lnPeos_m1, wrap_lnT_00, wrap_lnT_m1, wrap_Peos_00, &
         wrap_Peos_m1, wrap_r_00, wrap_T_00, wrap_T_m1, wrap_v_00, &
         wrap_v_p1, wrap_w_00, wrap_w_m1
      use hydro_rsp2, only: get_RSP2_alfa_beta_face_weights
      use reconstructed_face_support, only: get_reconstructed_face_state_ad
      use star_utils, only: get_rho_face_val
      use turb, only: set_TDC_LNA

      implicit none

      private
      public :: star_LNA_HSE_grav_term
      public :: star_LNA_eval_dlnPdm_qhse
      public :: star_LNA_perturb_convective_luminosity
      public :: rsp2_Ptrb_for_star_LNA
      public :: rsp2_luminosity_resid_for_star_LNA
      public :: rsp2_luminosity_terms_for_star_LNA
      public :: tdc_luminosity_resid_for_star_LNA
      public :: tdc_luminosity_terms_for_star_LNA
      public :: tdc_radiative_luminosity_for_star_LNA
      public :: frozen_flux_luminosity_resid_for_star_LNA
      public :: star_LNA_L_conv_ad
      public :: pressure_components_for_star_LNA
      public :: Ptot_for_star_LNA
      public :: d_mlt_Pturb_face_for_star_LNA
      public :: Uq_face_for_star_LNA
      public :: static_face_pressure_for_star_LNA
      public :: turbulent_viscous_heating_for_star_LNA
      public :: turbulent_energy_inertia_for_star_LNA
      public :: rsp2_turbulent_energy_rhs_for_star_LNA
      public :: rsp2_turbulent_energy_inertia_for_star_LNA
      public :: rsp2_forces_non_turbulent_cell
      public :: tdc_zero_w_for_star_LNA
      public :: tdc_relation_for_star_LNA
      public :: tdc_face_state_for_star_LNA
      public :: tdc_conv_vel_for_star_LNA
      public :: turbulent_tdc_energy_inertia_for_star_LNA
      public :: chi_coefficient_for_star_LNA
      public :: rsp2_chi_coefficient_for_star_LNA
      public :: tdc_chi_coefficient_for_star_LNA
      public :: tdc_forces_non_turbulent_cell
      public :: tdc_lna_active
      public :: rsp2_terms_for_star_LNA_audit

      contains

      ! -----------
      ! Shared LNA convection switches and hydrostatic-gradient helpers.
      ! -----------
      logical function star_LNA_perturb_convective_luminosity(s) result(perturb_luminosity)
         type(star_info), pointer :: s

         perturb_luminosity = s% star_LNA_perturb_convective_flux .and. &
            trim(s% star_LNA_convection_treatment) == 'perturbed'

         if (s% MLT_option == 'TDC') then
            perturb_luminosity = perturb_luminosity .and. s% star_LNA_include_tdc
         end if
      end function star_LNA_perturb_convective_luminosity


      logical function tdc_lna_active(s) result(active)
         type(star_info), pointer :: s

         active = .not. s% RSP2_flag .and. &
            s% MLT_option == 'TDC' .and. &
            s% star_LNA_include_tdc .and. &
            trim(s% star_LNA_convection_treatment) == 'perturbed'
      end function tdc_lna_active


      subroutine star_LNA_HSE_grav_term(s, k, grav_ad, area_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: grav_ad, area_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: r_ad, r2_ad

         ierr = 0
         r_ad = wrap_r_00(s, k)
         if (r_ad%val <= 0d0) then
            write(*,'(a,i0)') 'star_LNA found non-positive radius in momentum row at k = ', k
            ierr = -1
            return
         end if

         r2_ad = r_ad*r_ad
         area_ad = 4d0*pi*r2_ad
         grav_ad = -s% cgrav(k)*s% m_grav(k)/r2_ad
      end subroutine star_LNA_HSE_grav_term


      subroutine star_LNA_eval_dlnPdm_qhse(s, k, dlnPdm_qhse, Ppoint, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dlnPdm_qhse, Ppoint
         integer, intent(out) :: ierr
         real(dp) :: alfa
         type(auto_diff_real_star_order1) :: conv_vel_ad, grav, area, P00, Pm1, &
            mlt_Ptrb00, mlt_Ptrbm1

         ierr = 0
         call star_LNA_HSE_grav_term(s, k, grav, area, ierr)
         if (ierr /= 0) return

         if ((s% have_mlt_vc .and. s% okay_to_set_mlt_vc) .and. &
               s% include_mlt_Pturb_in_thermodynamic_gradients .and. &
               s% mlt_Pturb_factor > 0d0) then
            if (tdc_lna_active(s)) then
               call tdc_conv_vel_for_star_LNA(s, k, conv_vel_ad)
               mlt_Ptrb00 = s% mlt_Pturb_factor*pow2(conv_vel_ad)* &
                  wrap_d_00(s, k)/3d0
               if (k == 1) then
                  mlt_Ptrbm1 = 0d0
               else
                  mlt_Ptrbm1 = s% mlt_Pturb_factor*pow2(conv_vel_ad)* &
                     wrap_d_m1(s, k)/3d0
               end if
            else
               mlt_Ptrb00 = s% mlt_Pturb_factor*pow2(s% mlt_vc_old(k))* &
                  wrap_d_00(s, k)/3d0
               if (k == 1) then
                  mlt_Ptrbm1 = 0d0
               else
                  mlt_Ptrbm1 = s% mlt_Pturb_factor*pow2(s% mlt_vc_old(k))* &
                     wrap_d_m1(s, k)/3d0
               end if
            end if
         else
            mlt_Ptrb00 = 0d0
            mlt_Ptrbm1 = 0d0
         end if

         P00 = wrap_Peos_00(s, k)
         if (k == 1) then
            Ppoint = P00 + mlt_Ptrb00
         else
            Pm1 = wrap_Peos_m1(s, k)
            Pm1 = Pm1 + mlt_Ptrbm1
            P00 = P00 + mlt_Ptrb00
            alfa = s% dq(k - 1)/(s% dq(k - 1) + s% dq(k))
            Ppoint = alfa*P00 + (1d0 - alfa)*Pm1
         end if

         dlnPdm_qhse = grav/(area*Ppoint)
      end subroutine star_LNA_eval_dlnPdm_qhse


      subroutine actual_gradT_for_star_LNA(s, k, gradT_actual_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: gradT_actual_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: dlnPdm_ad, Ppoint_ad, &
            Tm1_ad, T00_ad, dT_ad, Tpoint_ad, lnTdiff_ad
         real(dp) :: delm, alfa

         ierr = 0
         if (k == 1) then
            gradT_actual_ad = s% gradT_ad(1)
            return
         end if

         if (s% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn) then
            gradT_actual_ad = (wrap_lnT_m1(s, k) - wrap_lnT_00(s, k))/ &
               (wrap_lnPeos_m1(s, k) - wrap_lnPeos_00(s, k))
            return
         end if

         call star_LNA_eval_dlnPdm_qhse(s, k, dlnPdm_ad, Ppoint_ad, ierr)
         if (ierr /= 0) return

         Tm1_ad = wrap_T_m1(s, k)
         T00_ad = wrap_T_00(s, k)
         dT_ad = Tm1_ad - T00_ad
         alfa = s% dm(k - 1)/(s% dm(k - 1) + s% dm(k))
         Tpoint_ad = alfa*T00_ad + (1d0 - alfa)*Tm1_ad
         lnTdiff_ad = dT_ad/Tpoint_ad
         delm = 0.5d0*(s% dm(k) + s% dm(k - 1))
         gradT_actual_ad = lnTdiff_ad/(delm*dlnPdm_ad)
      end subroutine actual_gradT_for_star_LNA

      ! -----------
      ! RSP2 pressure and luminosity closures.
      ! -----------
      subroutine rsp2_Ptrb_for_star_LNA(s, k, Ptrb_ad, Ptrb_div_etrb)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Ptrb_ad, Ptrb_div_etrb
         real(dp), parameter :: x_ALFAP = 2d0/3d0
         type(auto_diff_real_star_order1) :: rho_ad, etrb_ad

         Ptrb_ad = 0d0
         Ptrb_div_etrb = 0d0
         if (s% RSP2_alfap == 0d0 .or. s% mixing_length_alpha == 0d0 .or. &
               rsp2_forces_non_turbulent_cell(s, k)) return

         rho_ad = wrap_d_00(s, k)
         etrb_ad = wrap_etrb_00(s, k)
         Ptrb_div_etrb = s% RSP2_alfap*x_ALFAP*rho_ad
         Ptrb_ad = Ptrb_div_etrb*etrb_ad
      end subroutine rsp2_Ptrb_for_star_LNA


      subroutine pressure_components_for_star_LNA(s, k, Peos_ad, Ptrb_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Peos_ad, Ptrb_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Ptrb_div_etrb

         ierr = 0
         Peos_ad = wrap_Peos_00(s, k)
         Ptrb_ad = 0d0
         if (s% star_LNA_perturb_turbulent_pressure .and. s% RSP2_flag) &
            call rsp2_Ptrb_for_star_LNA(s, k, Ptrb_ad, Ptrb_div_etrb)
      end subroutine pressure_components_for_star_LNA


      subroutine Ptot_for_star_LNA(s, k, Ptot_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Ptot_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Peos_ad, Ptrb_ad

         call pressure_components_for_star_LNA(s, k, Peos_ad, Ptrb_ad, ierr)
         if (ierr /= 0) return
         Ptot_ad = Peos_ad + Ptrb_ad
      end subroutine Ptot_for_star_LNA


      subroutine d_mlt_Pturb_face_for_star_LNA(s, k, d_mlt_Pturb_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: d_mlt_Pturb_ad
         type(auto_diff_real_star_order1) :: conv_vel_ad, rho_00, rho_m1

         if (s% star_LNA_perturb_turbulent_pressure .and. &
               s% mlt_Pturb_factor > 0d0 .and. tdc_lna_active(s)) then
            call tdc_conv_vel_for_star_LNA(s, k, conv_vel_ad)
            rho_00 = wrap_d_00(s, k)
            rho_m1 = wrap_d_m1(s, k)
            d_mlt_Pturb_ad = &
               s% mlt_Pturb_factor*pow2(conv_vel_ad)*(rho_m1 - rho_00)/3d0
         else if (s% star_LNA_perturb_turbulent_pressure .and. &
               s% mlt_Pturb_factor > 0d0 .and. s% mlt_vc_old(k) > 0d0) then
            rho_00 = wrap_d_00(s, k)
            rho_m1 = wrap_d_m1(s, k)
            d_mlt_Pturb_ad = &
               s% mlt_Pturb_factor*s% mlt_vc_old(k)*s% mlt_vc_old(k)*(rho_m1 - rho_00)/3d0
         else
            d_mlt_Pturb_ad = 0d0
         end if
      end subroutine d_mlt_Pturb_face_for_star_LNA


      subroutine static_face_pressure_for_star_LNA(s, k, P_face, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: P_face
         integer, intent(out) :: ierr
         real(dp) :: alfa, beta
         type(auto_diff_real_star_order1) :: P_left_ad, P_right_ad, &
            T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
            chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
            scale_height_face_ad, gradr_face_ad

         ierr = 0
         if (k < 1 .or. k > s% nz) then
            write(*,'(a,i0)') 'star_LNA bad face index for work pressure: ', k
            ierr = -1
            return
         end if

         call Ptot_for_star_LNA(s, k, P_left_ad, ierr)
         if (ierr /= 0) return
         if (k == 1) then
            P_face = P_left_ad%val
         else
            call Ptot_for_star_LNA(s, k - 1, P_right_ad, ierr)
            if (ierr /= 0) return
            alfa = s% dq(k - 1)/(s% dq(k - 1) + s% dq(k))
            beta = 1d0 - alfa
            P_face = alfa*P_left_ad%val + beta*P_right_ad%val
         end if

         if (s% mlt_Pturb_factor > 0d0 .and. tdc_lna_active(s) .and. &
               s% mlt_vc(k) > 0d0) then
            call tdc_face_state_for_star_LNA(s, k, &
               T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
               chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
               scale_height_face_ad, gradr_face_ad, ierr)
            if (ierr /= 0) return
            P_face = P_face + &
               s% mlt_Pturb_factor*s% mlt_vc(k)*s% mlt_vc(k)* &
               rho_face_ad% val/3d0
         else if (s% mlt_Pturb_factor > 0d0 .and. s% mlt_vc_old(k) > 0d0) then
            P_face = P_face + &
               s% mlt_Pturb_factor*s% mlt_vc_old(k)*s% mlt_vc_old(k)* &
               get_rho_face_val(s, k)/3d0
         end if
      end subroutine static_face_pressure_for_star_LNA


      subroutine rsp2_luminosity_resid_for_star_LNA(s, k, luminosity_resid_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: luminosity_resid_ad
         type(auto_diff_real_star_order1) :: Lr_ad, Lc_ad, Lt_ad

         call rsp2_luminosity_terms_for_star_LNA(s, k, Lr_ad, Lc_ad, Lt_ad)
         luminosity_resid_ad = Lr_ad + Lc_ad + Lt_ad - wrap_L_00(s, k)
      end subroutine rsp2_luminosity_resid_for_star_LNA


      subroutine rsp2_luminosity_terms_for_star_LNA(s, k, Lr_ad, Lc_ad, Lt_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Lr_ad, Lc_ad, Lt_ad

         Lr_ad = s% Lr_ad(k)
         if (star_LNA_perturb_convective_luminosity(s)) then
            call rsp2_convective_luminosity_for_star_LNA(s, k, Lc_ad)
         else
            call rsp2_convective_luminosity_for_star_LNA(s, k, Lc_ad)
            Lc_ad = Lc_ad%val
         end if

         Lt_ad = s% Lt_ad(k)
      end subroutine rsp2_luminosity_terms_for_star_LNA


      subroutine rsp2_convective_luminosity_for_star_LNA(s, k, Lc_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Lc_ad
         real(dp), parameter :: x_ALFAS = 0.5d0*sqrt_2_div_3
         real(dp), parameter :: x_ALFAC = 0.5d0*sqrt_2_div_3
         real(dp) :: alfa, beta
         type(auto_diff_real_star_order1) :: r_ad, area_ad, T_rho_face_ad, &
            PII_face_ad, w_face_ad

         Lc_ad = 0d0
         if (k <= 1 .or. k > s% nz .or. rsp2_forces_non_turbulent_cell(s, k)) return

         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         r_ad = wrap_r_00(s, k)
         area_ad = 4d0*pi*pow2(r_ad)
         T_rho_face_ad = alfa*wrap_T_00(s, k)*wrap_d_00(s, k) + &
            beta*wrap_T_m1(s, k)*wrap_d_m1(s, k)
         call rsp2_PII_face_for_star_LNA(s, k, PII_face_ad)
         w_face_ad = alfa*wrap_w_00(s, k) + beta*wrap_w_m1(s, k)
         Lc_ad = w_face_ad*area_ad*(x_ALFAC/x_ALFAS)*T_rho_face_ad*PII_face_ad
      end subroutine rsp2_convective_luminosity_for_star_LNA


      subroutine rsp2_PII_face_for_star_LNA(s, k, PII_face_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: PII_face_ad
         real(dp), parameter :: x_ALFAS = 0.5d0*sqrt_2_div_3
         real(dp) :: alfa, beta
         type(auto_diff_real_star_order1) :: Cp_face_ad

         PII_face_ad = 0d0
         if (k <= 1 .or. k >= s% nz) return

         ! LNA uses the unsaturated convective response.  The nonlinear
         ! enthalpy-flux limiter is a hydro safeguard and is intentionally
         ! excluded from the linear operator.
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         Cp_face_ad = alfa*wrap_Cp_00(s, k) + beta*wrap_Cp_m1(s, k)
         PII_face_ad = x_ALFAS*s% mixing_length_alpha*Cp_face_ad*s% Y_face_ad(k)
      end subroutine rsp2_PII_face_for_star_LNA

      ! -----------
      ! TDC luminosity and velocity closures.
      ! -----------
      subroutine tdc_face_state_for_star_LNA(s, k, &
            T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
            chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
            scale_height_face_ad, gradr_face_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: &
            T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
            chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
            scale_height_face_ad, gradr_face_ad
         integer, intent(out) :: ierr

         call get_reconstructed_face_state_ad( &
            s, k, T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
            chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
            scale_height_face_ad, gradr_face_ad, ierr)
      end subroutine tdc_face_state_for_star_LNA


      subroutine tdc_luminosity_resid_for_star_LNA(s, k, luminosity_resid_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: luminosity_resid_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: L_rad_ad, L_conv_ad

         call tdc_luminosity_terms_for_star_LNA(s, k, L_rad_ad, L_conv_ad, ierr)
         if (ierr /= 0) return
         luminosity_resid_ad = L_rad_ad + L_conv_ad - wrap_L_00(s, k)
      end subroutine tdc_luminosity_resid_for_star_LNA


      subroutine tdc_luminosity_terms_for_star_LNA(s, k, L_rad_ad, L_conv_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: L_rad_ad, L_conv_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: luminosity_resid_ad, velocity_rhs_ad, &
            velocity_inertia_ad

         call tdc_relation_for_star_LNA(s, k, luminosity_resid_ad, velocity_rhs_ad, &
            velocity_inertia_ad, L_rad_ad, L_conv_ad, ierr)
      end subroutine tdc_luminosity_terms_for_star_LNA


      subroutine tdc_relation_for_star_LNA(s, k, luminosity_resid_ad, &
            velocity_rhs_ad, velocity_inertia_ad, L_rad_ad, L_conv_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: luminosity_resid_ad, &
            velocity_rhs_ad, velocity_inertia_ad, L_rad_ad, L_conv_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: A_ad, Eq_div_w_ad, gradT_actual_ad, &
            r_ad, grav_ad, T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, &
            Cp_face_ad, chiRho_face_ad, chiT_face_ad, grada_face_ad, &
            opacity_face_ad, scale_height_face_ad, gradr_face_ad

         ierr = 0
         if (s% mixing_length_alpha <= 0d0) then
            call tdc_radiative_luminosity_terms_for_star_LNA(s, k, &
               L_rad_ad, L_conv_ad, ierr)
            if (ierr /= 0) return
            luminosity_resid_ad = L_rad_ad - wrap_L_00(s, k)
            velocity_rhs_ad = 0d0
            velocity_inertia_ad = 0d0
            return
         end if

         call tdc_face_state_for_star_LNA(s, k, &
            T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
            chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
            scale_height_face_ad, gradr_face_ad, ierr)
         if (ierr /= 0) return
         call actual_gradT_for_star_LNA(s, k, gradT_actual_ad, ierr)
         if (ierr /= 0) return
         call tdc_A_for_star_LNA(s, k, A_ad)
         call tdc_Eq_div_w_for_star_LNA(s, k, Eq_div_w_ad, ierr)
         if (ierr /= 0) return

         r_ad = wrap_r_00(s, k)
         grav_ad = s% cgrav(k)*s% m_grav(k)/pow2(r_ad)
         ! LNA uses the unsaturated TDC relation; do not linearize through the
         ! nonlinear enthalpy-flux limiter.
         call set_TDC_LNA( &
            s% mixing_length_alpha, s% TDC_alpha_D, s% TDC_alpha_R, &
            s% TDC_alpha_Pt, s% cgrav(k), s% m_grav(k), &
            chiT_face_ad, chiRho_face_ad, wrap_L_00(s, k), &
            gradT_actual_ad, r_ad, Peos_face_ad, T_face_ad, &
            rho_face_ad, Cp_face_ad, opacity_face_ad, &
            scale_height_face_ad, s% gradL_ad(k), grada_face_ad, &
            A_ad, Eq_div_w_ad, grav_ad, s% include_mlt_corr_to_TDC, &
            s% TDC_alpha_C, s% TDC_alpha_S, .false., &
            energy_face_ad, luminosity_resid_ad, velocity_rhs_ad, &
            velocity_inertia_ad, L_rad_ad, L_conv_ad, ierr)
      end subroutine tdc_relation_for_star_LNA


      subroutine tdc_radiative_luminosity_terms_for_star_LNA(s, k, &
            L_rad_ad, L_conv_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: L_rad_ad, L_conv_ad
         integer, intent(out) :: ierr

         call tdc_radiative_luminosity_for_star_LNA(s, k, L_rad_ad, ierr)
         if (ierr /= 0) return
         L_conv_ad = 0d0
      end subroutine tdc_radiative_luminosity_terms_for_star_LNA


      subroutine tdc_radiative_luminosity_for_star_LNA(s, k, L_rad_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: L_rad_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: gradT_actual_ad, L0_ad, &
            T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
            chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
            scale_height_face_ad, gradr_face_ad

         ierr = 0
         call tdc_face_state_for_star_LNA(s, k, &
            T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
            chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
            scale_height_face_ad, gradr_face_ad, ierr)
         if (ierr /= 0) return
         call actual_gradT_for_star_LNA(s, k, gradT_actual_ad, ierr)
         if (ierr /= 0) return
         L0_ad = (16d0*pi*crad*clight/3d0)*s% cgrav(k)*s% m_grav(k)* &
            pow4(T_face_ad)/(Peos_face_ad*opacity_face_ad)
         L_rad_ad = L0_ad*gradT_actual_ad
      end subroutine tdc_radiative_luminosity_for_star_LNA


      subroutine frozen_flux_luminosity_resid_for_star_LNA(s, k, luminosity_resid_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: luminosity_resid_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: L_rad_ad, L_conv0_ad

         call tdc_radiative_luminosity_for_star_LNA(s, k, L_rad_ad, ierr)
         if (ierr /= 0) return

         ! Frozen-flux LNA keeps the equilibrium convective luminosity but sets
         ! its perturbation to zero, so only the radiative luminosity carries
         ! AD partials.
         L_conv0_ad = s% L_conv(k)
         luminosity_resid_ad = L_rad_ad + L_conv0_ad - wrap_L_00(s, k)
      end subroutine frozen_flux_luminosity_resid_for_star_LNA


      function star_LNA_L_conv_ad(s, k) result(L_conv_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: L_conv_ad
         type(auto_diff_real_star_order1) :: L_ad, L_rad_ad
         integer :: ierr

         if (k < 1 .or. k > s% nz) then
            L_conv_ad = 0d0
            return
         end if

         if (.not. star_LNA_perturb_convective_luminosity(s)) then
            L_conv_ad = s% L_conv(k)
            return
         end if

         if (tdc_lna_active(s)) then
            call tdc_luminosity_terms_for_star_LNA(s, k, L_rad_ad, L_conv_ad, ierr)
            if (ierr /= 0) L_conv_ad = 0d0
            return
         end if

         if (s% mlt_mixing_type(k) == no_mixing .or. abs(s% gradr(k)) < 1d-20) then
            L_conv_ad = 0d0
            return
         end if

         L_ad = wrap_L_00(s, k)
         L_conv_ad = L_ad*(1d0 - s% gradT_ad(k)/s% gradr_ad(k))
      end function star_LNA_L_conv_ad


      logical function tdc_zero_w_for_star_LNA(s, k) result(force_zero_w)
         type(star_info), pointer :: s
         integer, intent(in) :: k

         force_zero_w = s% mixing_length_alpha <= 0d0 .or. &
            tdc_forces_non_turbulent_cell(s, k) .or. &
            s% mlt_mixing_type(k) == no_mixing .or. &
            s% mlt_vc(k) <= 0d0
      end function tdc_zero_w_for_star_LNA


      subroutine tdc_A_for_star_LNA(s, k, A_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: A_ad

         if (tdc_zero_w_for_star_LNA(s, k)) then
            A_ad = 0d0
            return
         end if

         call wrap(A_ad, s% mlt_vc(k)/sqrt_2_div_3, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 1d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0, &
            0d0, 0d0, 0d0)
      end subroutine tdc_A_for_star_LNA


      subroutine tdc_conv_vel_for_star_LNA(s, k, conv_vel_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: conv_vel_ad
         type(auto_diff_real_star_order1) :: A_ad

         call tdc_A_for_star_LNA(s, k, A_ad)
         conv_vel_ad = sqrt_2_div_3*A_ad
      end subroutine tdc_conv_vel_for_star_LNA


      subroutine tdc_Eq_div_w_for_star_LNA(s, k, Eq_div_w_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Eq_div_w_ad
         integer, intent(out) :: ierr

         ierr = 0
         Eq_div_w_ad = 0d0
         if (.not. associated(s)) then
            ierr = -1
            return
         end if
         if (k < 1 .or. k > s% nz) return
         ! Eddy-viscous heating is quadratic in the velocity-gradient
         ! perturbation, so it is absent from the static first-order LNA.
         ! Keep this helper explicit so the TDC source equation stays aligned
         ! with the normal TDC notation without adding a nonlinear source.
      end subroutine tdc_Eq_div_w_for_star_LNA

      ! -----------
      ! RSP2 and TDC turbulent-energy closures.
      ! -----------
      subroutine rsp2_turbulent_energy_rhs_for_star_LNA(s, k, rhs_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: rhs_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: coupling_ad, dLt_dm_ad, dwork_dm_ad

         ierr = 0
         call rsp2_coupling_for_star_LNA(s, k, coupling_ad)
         call rsp2_dLt_dm_for_star_LNA(s, k, dLt_dm_ad)
         call rsp2_turbulent_pressure_work_dm_for_star_LNA(s, k, dwork_dm_ad, ierr)
         if (ierr /= 0) return
         rhs_ad = coupling_ad - dLt_dm_ad - dwork_dm_ad
      end subroutine rsp2_turbulent_energy_rhs_for_star_LNA


      subroutine turbulent_viscous_heating_for_star_LNA(s, k, Eq_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Eq_ad
         integer, intent(out) :: ierr

         ierr = 0
         Eq_ad = 0d0
      end subroutine turbulent_viscous_heating_for_star_LNA


      subroutine turbulent_energy_inertia_for_star_LNA(s, k, turbulent_inertia_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: turbulent_inertia_ad
         integer, intent(out) :: ierr

         ierr = 0
         turbulent_inertia_ad = 0d0
         if (.not. s% star_LNA_perturb_turbulent_energy) return

         if (s% RSP2_flag) then
            turbulent_inertia_ad = wrap_etrb_00(s, k)
         else if (s% star_LNA_include_tdc .and. &
               s% MLT_option == 'TDC' .and. &
               s% TDC_include_eturb_in_energy_equation) then
            if (tdc_lna_active(s)) then
               call turbulent_tdc_energy_inertia_for_star_LNA(s, k, &
                  turbulent_inertia_ad, ierr)
               return
            end if
            if (k < s% nz) then
               turbulent_inertia_ad = 0.75d0*(pow2(s% mlt_vc_ad(k)) + &
                  pow2(shift_p1(s% mlt_vc_ad(k + 1))))
            else
               turbulent_inertia_ad = 0.75d0*pow2(s% mlt_vc_ad(k))
            end if
         end if
      end subroutine turbulent_energy_inertia_for_star_LNA


      subroutine rsp2_coupling_for_star_LNA(s, k, coupling_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: coupling_ad
         type(auto_diff_real_star_order1) :: source_ad, damping_ad, rad_damping_ad

         coupling_ad = 0d0
         if (rsp2_forces_non_turbulent_cell(s, k)) return

         call rsp2_source_for_star_LNA(s, k, source_ad)
         call rsp2_damping_for_star_LNA(s, k, damping_ad)
         call rsp2_radiative_damping_for_star_LNA(s, k, rad_damping_ad)
         coupling_ad = source_ad - damping_ad - rad_damping_ad
      end subroutine rsp2_coupling_for_star_LNA


      subroutine rsp2_source_for_star_LNA(s, k, source_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: source_ad
         type(auto_diff_real_star_order1) :: Hp_face_00_ad, Hp_face_p1_ad, &
            PII_face_00_ad, PII_face_p1_ad, PII_div_Hp_cell_ad, QQ_00_ad, &
            P_QQ_div_Cp_ad

         source_ad = 0d0
         if (rsp2_forces_non_turbulent_cell(s, k)) return

         Hp_face_00_ad = wrap_Hp_00(s, k)
         call rsp2_PII_face_for_star_LNA(s, k, PII_face_00_ad)
         if (k == s% nz) then
            PII_div_Hp_cell_ad = PII_face_00_ad/Hp_face_00_ad
         else
            Hp_face_p1_ad = wrap_Hp_p1(s, k)
            call rsp2_PII_face_for_star_LNA(s, k + 1, PII_face_p1_ad)
            PII_face_p1_ad = shift_p1(PII_face_p1_ad)
            PII_div_Hp_cell_ad = 0.5d0*( &
               PII_face_00_ad/Hp_face_00_ad + PII_face_p1_ad/Hp_face_p1_ad)
         end if

         QQ_00_ad = wrap_chiT_00(s, k)/( &
            wrap_d_00(s, k)*wrap_T_00(s, k)*wrap_chiRho_00(s, k))
         P_QQ_div_Cp_ad = wrap_Peos_00(s, k)*QQ_00_ad/wrap_Cp_00(s, k)
         source_ad = (wrap_w_00(s, k) + s% RSP2_source_seed)* &
            PII_div_Hp_cell_ad*wrap_T_00(s, k)*P_QQ_div_Cp_ad
      end subroutine rsp2_source_for_star_LNA


      subroutine rsp2_damping_for_star_LNA(s, k, damping_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: damping_ad
         real(dp), parameter :: x_CEDE = (8d0/3d0)*sqrt_2_div_3
         type(auto_diff_real_star_order1) :: Hp_cell_ad, w_00_ad

         damping_ad = 0d0
         if (s% mixing_length_alpha == 0d0) return

         Hp_cell_ad = rsp2_Hp_cell_for_star_LNA(s, k)
         w_00_ad = wrap_w_00(s, k)
         damping_ad = (s% RSP2_alfad*x_CEDE/s% mixing_length_alpha)* &
            (pow3(w_00_ad) - pow3(s% RSP2_w_min_for_damping))/Hp_cell_ad
      end subroutine rsp2_damping_for_star_LNA


      subroutine rsp2_radiative_damping_for_star_LNA(s, k, rad_damping_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: rad_damping_ad
         real(dp), parameter :: x_GAMMAR = 2d0*sqrt(3d0)
         real(dp) :: alpha, gammar
         type(auto_diff_real_star_order1) :: Hp_cell_ad, POM2_ad, w_00_ad

         rad_damping_ad = 0d0
         alpha = s% mixing_length_alpha
         gammar = s% RSP2_alfar*x_GAMMAR
         if (alpha == 0d0 .or. gammar == 0d0) return

         Hp_cell_ad = rsp2_Hp_cell_for_star_LNA(s, k)
         w_00_ad = wrap_w_00(s, k)
         POM2_ad = pow3(wrap_T_00(s, k))/( &
            pow2(wrap_d_00(s, k))*wrap_Cp_00(s, k)*wrap_kap_00(s, k))
         rad_damping_ad = pow2(w_00_ad)* &
            4d0*boltz_sigma*pow2(gammar/alpha)*POM2_ad/pow2(Hp_cell_ad)
      end subroutine rsp2_radiative_damping_for_star_LNA


      function rsp2_Hp_cell_for_star_LNA(s, k) result(Hp_cell_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Hp_cell_ad, grav_cell_ad

         grav_cell_ad = 0.5d0*s% cgrav(k)*s% m_grav(k)/pow2(wrap_r_00(s, k))
         if (k < s% nz) &
            grav_cell_ad = grav_cell_ad + &
               0.5d0*s% cgrav(k+1)*s% m_grav(k+1)/pow2(wrap_r_p1(s, k))
         Hp_cell_ad = wrap_Peos_00(s, k)/(wrap_d_00(s, k)*grav_cell_ad)
      end function rsp2_Hp_cell_for_star_LNA


      subroutine rsp2_turbulent_pressure_work_dm_for_star_LNA(s, k, dwork_dm_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dwork_dm_ad
         integer, intent(out) :: ierr
         real(dp) :: Ptrb
         type(auto_diff_real_star_order1) :: Ptrb_ad, Ptrb_div_etrb, dVdt_ad

         ierr = 0
         dwork_dm_ad = 0d0
         if (.not. s% star_LNA_perturb_turbulent_pressure) return
         call rsp2_Ptrb_for_star_LNA(s, k, Ptrb_ad, Ptrb_div_etrb)
         Ptrb = Ptrb_ad%val
         if (Ptrb == 0d0) return

         call cell_dVdt_for_star_LNA(s, k, dVdt_ad)
         dwork_dm_ad = Ptrb*dVdt_ad/s% dm(k)
      end subroutine rsp2_turbulent_pressure_work_dm_for_star_LNA


      subroutine rsp2_dLt_dm_for_star_LNA(s, k, dLt_dm_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dLt_dm_ad
         type(auto_diff_real_star_order1) :: Lt_p1_ad

         if (k < s% nz) then
            Lt_p1_ad = shift_p1(s% Lt_ad(k + 1))
         else
            Lt_p1_ad = 0d0
         end if
         dLt_dm_ad = (s% Lt_ad(k) - Lt_p1_ad)/s% dm(k)
      end subroutine rsp2_dLt_dm_for_star_LNA


      subroutine rsp2_turbulent_energy_inertia_for_star_LNA(s, k, inertia_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: inertia_ad
         integer, intent(out) :: ierr

         ierr = 0
         inertia_ad = wrap_etrb_00(s, k)
      end subroutine rsp2_turbulent_energy_inertia_for_star_LNA


      subroutine rsp2_terms_for_star_LNA_audit(s, k, PII_ad, Lc_ad, Lt_ad, &
            source_ad, damping_ad, rad_damping_ad, Ptrb_div_etrb_ad, &
            dwork_dm_ad, dLt_dm_ad, rhs_ad, inertia_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: PII_ad, Lc_ad, Lt_ad, &
            source_ad, damping_ad, rad_damping_ad, Ptrb_div_etrb_ad, &
            dwork_dm_ad, dLt_dm_ad, rhs_ad, inertia_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lr_ad, Ptrb_ad

         ierr = 0
         call rsp2_PII_face_for_star_LNA(s, k, PII_ad)
         call rsp2_luminosity_terms_for_star_LNA(s, k, Lr_ad, Lc_ad, Lt_ad)
         call rsp2_source_for_star_LNA(s, k, source_ad)
         call rsp2_damping_for_star_LNA(s, k, damping_ad)
         call rsp2_radiative_damping_for_star_LNA(s, k, rad_damping_ad)
         call rsp2_Ptrb_for_star_LNA(s, k, Ptrb_ad, Ptrb_div_etrb_ad)
         call rsp2_turbulent_pressure_work_dm_for_star_LNA(s, k, dwork_dm_ad, ierr)
         if (ierr /= 0) return
         call rsp2_dLt_dm_for_star_LNA(s, k, dLt_dm_ad)
         call rsp2_turbulent_energy_rhs_for_star_LNA(s, k, rhs_ad, ierr)
         if (ierr /= 0) return
         call rsp2_turbulent_energy_inertia_for_star_LNA(s, k, inertia_ad, ierr)
      end subroutine rsp2_terms_for_star_LNA_audit


      subroutine turbulent_tdc_energy_inertia_for_star_LNA(s, k, inertia_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: inertia_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: vc_00_ad, vc_p1_ad

         ierr = 0
         call tdc_conv_vel_for_star_LNA(s, k, vc_00_ad)
         if (k < s% nz) then
            call tdc_conv_vel_for_star_LNA(s, k + 1, vc_p1_ad)
            vc_p1_ad = shift_p1(vc_p1_ad)
            inertia_ad = 0.75d0*(pow2(vc_00_ad) + pow2(vc_p1_ad))
         else
            inertia_ad = 0.75d0*pow2(vc_00_ad)
         end if
      end subroutine turbulent_tdc_energy_inertia_for_star_LNA


      subroutine cell_dVdt_for_star_LNA(s, k, dVdt_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dVdt_ad
         type(auto_diff_real_star_order1) :: v_inner_ad
         real(dp) :: r_inner

         if (k < s% nz) then
            r_inner = s% r(k + 1)
            v_inner_ad = wrap_v_p1(s, k)
         else
            r_inner = s% r_center
            v_inner_ad = 0d0
         end if

         dVdt_ad = 4d0*pi*(s% r(k)*s% r(k)*wrap_v_00(s, k) - &
            r_inner*r_inner*v_inner_ad)
      end subroutine cell_dVdt_for_star_LNA


      subroutine Uq_face_for_star_LNA(s, k, Uq_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Uq_ad
         integer, intent(out) :: ierr

         ierr = 0
         if (.not. s% star_LNA_perturb_eddy_viscosity) then
            Uq_ad = 0d0
         else
            call static_eddy_Uq_face_for_star_LNA(s, k, Uq_ad, ierr)
         end if
      end subroutine Uq_face_for_star_LNA


      subroutine static_eddy_Uq_face_for_star_LNA(s, k, Uq_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Uq_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_00, Chi_m1, r_ad

         ierr = 0
         call static_eddy_Chi_cell_for_star_LNA(s, k, Chi_00, ierr)
         if (ierr /= 0) return
         if (k > 1) then
            call static_eddy_Chi_cell_for_star_LNA(s, k - 1, Chi_m1, ierr)
            if (ierr /= 0) return
            Chi_m1 = shift_m1(Chi_m1)
         else
            Chi_m1 = 0d0
         end if

         r_ad = wrap_r_00(s, k)
         Uq_ad = 4d0*pi*(Chi_m1 - Chi_00)/(r_ad*s% dm_bar(k))
      end subroutine static_eddy_Uq_face_for_star_LNA


      subroutine static_eddy_Chi_cell_for_star_LNA(s, k, Chi_ad, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: Chi_ad
         integer, intent(out) :: ierr
         real(dp) :: chi_coeff
         type(auto_diff_real_star_order1) :: d_v_div_r_ad

         ierr = 0
         Chi_ad = 0d0
         chi_coeff = chi_coefficient_for_star_LNA(s, k, ierr)
         if (ierr /= 0) return
         if (chi_coeff <= 0d0) return

         call static_d_v_div_r_for_star_LNA(s, k, d_v_div_r_ad)
         Chi_ad = chi_coeff*d_v_div_r_ad
      end subroutine static_eddy_Chi_cell_for_star_LNA


      subroutine static_d_v_div_r_for_star_LNA(s, k, d_v_div_r_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: d_v_div_r_ad
         real(dp) :: rp1

         if (k < s% nz) then
            rp1 = s% r(k + 1)
         else
            rp1 = s% r_center
         end if
         if (rp1 == 0d0) rp1 = 1d0

         d_v_div_r_ad = wrap_v_00(s, k)/s% r(k) - wrap_v_p1(s, k)/rp1
      end subroutine static_d_v_div_r_for_star_LNA

      ! -----------
      ! Forced-zero zones and eddy-viscosity closures.
      ! -----------
      logical function rsp2_forces_non_turbulent_cell(s, k) result(force_zero_w)
         type(star_info), pointer :: s
         integer, intent(in) :: k

         force_zero_w = s% mixing_length_alpha == 0d0 .or. &
            k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
            k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)
      end function rsp2_forces_non_turbulent_cell


      logical function tdc_forces_non_turbulent_cell(s, k) result(force_zero_w)
         type(star_info), pointer :: s
         integer, intent(in) :: k

         force_zero_w = k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
            k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent
      end function tdc_forces_non_turbulent_cell


      real(dp) function chi_coefficient_for_star_LNA(s, k, ierr) result(chi_coeff)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr

         ierr = 0
         chi_coeff = 0d0
         if (s% RSP2_flag) then
            chi_coeff = rsp2_chi_coefficient_for_star_LNA(s, k)
         else if (s% star_LNA_include_tdc .and. s% MLT_option == 'TDC') then
            chi_coeff = tdc_chi_coefficient_for_star_LNA(s, k, ierr)
         end if
      end function chi_coefficient_for_star_LNA


      real(dp) function rsp2_chi_coefficient_for_star_LNA(s, k) result(chi_coeff)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: alfam_alpha

         chi_coeff = 0d0
         alfam_alpha = s% RSP2_alfam*s% mixing_length_alpha
         if (alfam_alpha == 0d0 .or. rsp2_forces_non_turbulent_cell(s, k)) return

         chi_coeff = chi_geometry_coefficient_for_star_LNA(s, k, alfam_alpha, &
            Hp_cell_for_rsp2_chi_for_star_LNA(s, k), s% w(k))
      end function rsp2_chi_coefficient_for_star_LNA


      real(dp) function tdc_chi_coefficient_for_star_LNA(s, k, ierr) result(chi_coeff)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: alfam_alpha, Hp_cell

         ierr = 0
         chi_coeff = 0d0
         alfam_alpha = s% TDC_alpha_M*s% mixing_length_alpha
         if (alfam_alpha == 0d0 .or. tdc_forces_non_turbulent_cell(s, k)) return

         call Hp_cell_for_tdc_chi_for_star_LNA(s, k, Hp_cell, ierr)
         if (ierr /= 0) return
         chi_coeff = chi_geometry_coefficient_for_star_LNA(s, k, alfam_alpha, &
            Hp_cell, tdc_w_for_star_LNA(s, k))
      end function tdc_chi_coefficient_for_star_LNA


      real(dp) function chi_geometry_coefficient_for_star_LNA( &
            s, k, alfam_alpha, Hp_cell, w_cell) result(chi_coeff)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: alfam_alpha, Hp_cell, w_cell
         real(dp) :: r00, rp1, r6_cell

         if (w_cell <= 0d0 .or. Hp_cell <= 0d0) then
            chi_coeff = 0d0
            return
         end if

         r00 = s% r(k)
         if (k < s% nz) then
            rp1 = s% r(k + 1)
         else
            rp1 = s% r_center
         end if
         r6_cell = 0.5d0*(r00**6 + rp1**6)
         chi_coeff = (16d0/3d0)*pi*alfam_alpha*s% rho(k)*s% rho(k)* &
            r6_cell*Hp_cell*w_cell/s% dm(k)
      end function chi_geometry_coefficient_for_star_LNA


      real(dp) function Hp_cell_for_rsp2_chi_for_star_LNA(s, k) result(Hp_cell)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: grav_cell

         grav_cell = 0.5d0*s% cgrav(k)*s% m_grav(k)/pow2(s% r(k))
         if (k < s% nz) &
            grav_cell = grav_cell + &
               0.5d0*s% cgrav(k+1)*s% m_grav(k+1)/pow2(s% r(k+1))
         Hp_cell = s% Peos(k)/(s% rho(k)*grav_cell)
      end function Hp_cell_for_rsp2_chi_for_star_LNA


      subroutine Hp_cell_for_tdc_chi_for_star_LNA(s, k, Hp_cell, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: Hp_cell
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
            chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
            scale_height_face_ad, gradr_face_ad

         ierr = 0
         call tdc_face_state_for_star_LNA(s, k, &
            T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
            chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
            scale_height_face_ad, gradr_face_ad, ierr)
         if (ierr /= 0) return
         Hp_cell = 0.5d0*scale_height_face_ad% val

         if (k + 1 < s% nz) then
            call tdc_face_state_for_star_LNA(s, k + 1, &
               T_face_ad, rho_face_ad, Peos_face_ad, energy_face_ad, Cp_face_ad, &
               chiRho_face_ad, chiT_face_ad, grada_face_ad, opacity_face_ad, &
               scale_height_face_ad, gradr_face_ad, ierr)
            if (ierr /= 0) return
            Hp_cell = Hp_cell + 0.5d0*scale_height_face_ad% val
         end if
      end subroutine Hp_cell_for_tdc_chi_for_star_LNA


      real(dp) function tdc_w_for_star_LNA(s, k) result(w_cell)
         type(star_info), pointer :: s
         integer, intent(in) :: k

         if (k < s% nz) then
            if (s% okay_to_set_mlt_vc .and. &
                  s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then
               w_cell = 0.5d0*(s% mlt_vc_old(k) + s% mlt_vc_old(k + 1))/sqrt_2_div_3
            else
               w_cell = 0.5d0*(s% mlt_vc(k) + s% mlt_vc(k + 1))/sqrt_2_div_3
            end if
         else
            if (s% okay_to_set_mlt_vc .and. &
                  s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then
               w_cell = 0.5d0*s% mlt_vc_old(k)/sqrt_2_div_3
            else
               w_cell = 0.5d0*s% mlt_vc(k)/sqrt_2_div_3
            end if
         end if
      end function tdc_w_for_star_LNA

      end module star_lna_turbulence_closures
