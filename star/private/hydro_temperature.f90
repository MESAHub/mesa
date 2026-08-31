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

      module hydro_temperature

      use star_private_def
      use const_def, only: dp, ln10, pi4, crad, clight, convective_mixing
      use reconstructed_face_support, only: get_reconstructed_face_eos_kap_ad, &
         get_effective_gradr_factor_ad, get_Lrad_per_gradT_face_ad
      use utils_lib, only: mesa_error, is_bad
      use auto_diff
      use auto_diff_support

      implicit none

      private
      public :: do1_alt_dlnT_dm_eqn
      public :: do1_gradT_eqn
      public :: do1_dlnT_dm_eqn
      public :: do1_constant_L_eqn

      contains

      ! just relate L_rad to T gradient.
      ! d_P_rad/dm = -<opacity_face>*L_rad/(clight*area^2) -- see, e.g., K&W (5.12)
      ! P_rad = (1/3)*crad*T^4
      ! d_P_rad/dm = (crad/3)*(T(k-1)^4 - T(k)^4)/dm_bar
      ! L_rad = L - L_non_rad, L_non_rad = L_start - L_rad_start
      ! L_rad_start = (-d_P_rad/dm_bar*clight*area^2/<opacity_face>)_start
      subroutine do1_alt_dlnT_dm_eqn(s, k, nvar, ierr)
         use eos_def
         use star_utils, only: save_eqn_residual_info, get_T_face, get_Peos_face, get_kap_face
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr

         real(dp) :: scale, dm_bar
         type(auto_diff_real_star_order1) :: L_ad, r_00, area, area2, Lrad_ad, &
            opacity_face, kap_face, L0_ad, gradr_factor, &
            d_P_rad_expected_ad, T_m1, T4_m1, T_00, T4_00, &
            P_rad_m1, P_rad_00, d_P_rad_actual_ad, resid
         type(auto_diff_real_star_order1) :: T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face
         type(auto_diff_real_star_order1) :: flxR, flxLambda

         integer :: i_equL
         logical :: dbg
         logical :: test_partials

         include 'formats'
         ierr = 0
         i_equL = s% i_equL
         if (i_equL == 0) return

         if (.not. s% use_dPrad_dm_form_of_T_gradient_eqn) then
            ierr = -1
            return
         end if

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         dbg = .false.

         scale = s% energy_start(k)*s% rho_start(k)
         dm_bar = s% dm_bar(k)
         L_ad = wrap_L_00(s,k)
         r_00 = wrap_r_00(s,k)
         area = pi4*pow2(r_00); area2 = pow2(area)

         if (s% use_face_reconstruction) then
            call get_reconstructed_face_eos_kap_ad( &
               s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, opacity_face, ierr)
            if (ierr /= 0) return
         else
            T_face = get_T_face(s, k)
            P_face = get_Peos_face(s, k)
            opacity_face = get_kap_face(s, k)
         end if

         gradr_factor = get_effective_gradr_factor_ad(s, k)
         ! RTI can replace the final chemical-mixing label while MLT remains active.
         if (s% RSP2_flag) then
            Lrad_ad = s% Lr_ad(k)
         else if (s% lnT(k)/ln10 <= s% max_logT_for_mlt &
               .and. s% mlt_mixing_type(k) == convective_mixing &
               .and. abs(gradr_factor%val) > 1d-20) then
            ! Evaluate L/gradr analytically so the split is finite at zero luminosity.
            L0_ad = get_Lrad_per_gradT_face_ad( &
               s, k, T_face, P_face, opacity_face, gradr_factor)
            Lrad_ad = L0_ad*s% gradT_ad(k)  ! C&G 14.109
         else
            Lrad_ad = L_ad
         end if

         kap_face = opacity_face
         if (kap_face%val < s% min_kap_for_dPrad_dm_eqn) &
            kap_face = s% min_kap_for_dPrad_dm_eqn

         ! calculate expected d_P_rad from current L_rad
         d_P_rad_expected_ad = -dm_bar*kap_face*Lrad_ad/(clight*area2)

         ! calculate actual d_P_rad in current model
         T_m1 = wrap_T_m1(s,k); T4_m1 = pow4(T_m1)
         T_00 = wrap_T_00(s,k); T4_00 = pow4(T_00)

         P_rad_m1 = (crad/3._dp)*T4_m1
         P_rad_00 = (crad/3._dp)*T4_00
         d_P_rad_actual_ad = P_rad_m1 - P_rad_00

         ! enable flux-limited radiation transport derived by Levermore & Pomraning 1981
         s% flux_limit_R(k) = 0._dp
         s% flux_limit_lambda(k) =0._dp
         if (s% use_flux_limiting_with_dPrad_dm_form) then
            ! calculate the flux ratio R
            flxR = area * abs(T4_m1 - T4_00) / dm_bar / &
                  (kap_face * 0.5_dp * (T4_m1 + T4_00))

            s% flux_limit_R(k) = flxR%val

            ! calculate the flux limiter lambda
            flxLambda = (6._dp + 3._dp*flxR) / (6._dp + (3._dp + flxR)*flxR)

            s% flux_limit_lambda(k) = flxLambda%val

            ! calculate d_P_rad given the flux limiter
            d_P_rad_expected_ad = d_P_rad_expected_ad / flxLambda
         end if

         ! residual
         resid = (d_P_rad_expected_ad - d_P_rad_actual_ad)/scale
         s% equ(i_equL, k) = resid%val

         if (is_bad(resid%val)) then
!$OMP critical (star_alt_dlntdm_bad_num)
            write(*,2) 'resid%val', k, resid%val
            if (s% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'do1_alt_dlnT_dm_eqn')
!$OMP end critical (star_alt_dlntdm_bad_num)
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% gradT(k)
         end if

         call save_eqn_residual_info( &
            s, k, nvar, i_equL, resid, 'do1_alt_dlnT_dm_eqn', ierr)

         if (test_partials) then
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do1_alt_dlnT_dm_eqn', s% solver_test_partials_var
         end if

         contains

      end subroutine do1_alt_dlnT_dm_eqn


      subroutine do1_gradT_eqn(s, k, nvar, ierr)
         use eos_def
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: &
            resid, gradT, dlnT, dlnP
         integer :: i_equL
         logical :: test_partials

         include 'formats'
         ierr = 0

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         i_equL = s% i_equL
         if (i_equL == 0) return

         gradT = s% gradT_ad(k)
         dlnT = wrap_lnT_m1(s,k) - wrap_lnT_00(s,k)
         dlnP = wrap_lnPeos_m1(s,k) - wrap_lnPeos_00(s,k)

         resid = gradT*dlnP - dlnT
         s% equ(i_equL, k) = resid%val

         if (is_bad(s% equ(i_equL, k))) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            if (s% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'do1_gradT_eqn')
            return
            write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            write(*,2) 'gradT', k, gradT
            write(*,2) 'dlnT', k, dlnT
            write(*,2) 'dlnP', k, dlnP
            call mesa_error(__FILE__,__LINE__,'do1_gradT_eqn')
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% equ(i_equL,k)
         end if

         call save_eqn_residual_info( &
            s, k, nvar, i_equL, resid, 'do1_gradT_eqn', ierr)

         !call set_xtras

         contains

         subroutine set_xtras
            use auto_diff_support
            use star_utils, only: get_Lrad
            type(auto_diff_real_star_order1) :: &
               T4m1, T400, kap_m1, kap_00, alfa, beta, kap_face, &
               diff_T4_div_kap
            T4m1 = pow4(wrap_T_m1(s,k))
            T400 = pow4(wrap_T_00(s,k))
            kap_m1 = wrap_kap_m1(s,k)
            kap_00 = wrap_kap_00(s,k)
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1d0 - alfa
            kap_face = alfa*kap_00 + beta*kap_m1
            diff_T4_div_kap = (T4m1 - T400)/kap_face
            s% xtra1_array(k) = s% T_start(k)
            s% xtra2_array(k) = T4m1%val - T400%val
            s% xtra3_array(k) = kap_face%val
            s% xtra4_array(k) = diff_T4_div_kap%val
            s% xtra5_array(k) = get_Lrad(s,k)
            s% xtra6_array(k) = 1
         end subroutine set_xtras

      end subroutine do1_gradT_eqn


      subroutine do1_dlnT_dm_eqn(s, k, nvar, ierr)
         use eos_def
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: resid, &
            dlnPdm, Ppoint, gradT, dlnTdm, T00, Tm1, dT, Tpoint, lnTdiff
         real(dp) :: delm, alfa
         integer :: i_equL
         logical :: test_partials

         include 'formats'
         ierr = 0

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         i_equL = s% i_equL
         if (i_equL == 0) return

         if (s% constant_L) then
            call do1_constant_L_eqn(s, k, nvar, ierr)
            return
         end if

         if (k ==1 .and. s% use_RSP_L_eqn_outer_BC) then
            call set_RSP_Lsurf_BC(s, nvar, ierr)
            return
         end if

         if (s% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn) then
            call do1_gradT_eqn(s, k, nvar, ierr)
            return
         end if

         if (s% use_dPrad_dm_form_of_T_gradient_eqn) then
            call do1_alt_dlnT_dm_eqn(s, k, nvar, ierr)
            return
         end if

         ! dT/dm = dP/dm * T/P * grad_T, grad_T = dlnT/dlnP from MLT.
         ! but use hydrostatic value for dP/dm in this.
         ! this is because of limitations of MLT for calculating grad_T.
         ! (MLT assumes hydrostatic equilibrium)
         ! see comment in K&W chpt 9.1.

         call eval_dlnPdm_qhse(s, k, dlnPdm, Ppoint, ierr)
         if (ierr /= 0) return

         gradT = s% gradT_ad(k)
         dlnTdm = dlnPdm*gradT

         Tm1 = wrap_T_m1(s,k)
         T00 = wrap_T_00(s,k)
         dT = Tm1 - T00
         alfa = s% dm(k-1)/(s% dm(k-1) + s% dm(k))
         Tpoint = alfa*T00 + (1d0 - alfa)*Tm1
         lnTdiff = dT/Tpoint  ! use this in place of lnT(k-1)-lnT(k)
         delm = (s% dm(k) + s% dm(k-1))/2

         resid = delm*dlnTdm - lnTdiff
         s% equ(i_equL, k) = resid%val

         if (is_bad(s% equ(i_equL, k))) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            if (s% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'hydro eqns')
            return
            write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            write(*,2) 'lnTdiff', k, lnTdiff
            write(*,2) 'delm', k, delm
            write(*,2) 'dlnPdm', k, dlnPdm
            write(*,2) 'gradT', k, gradT
            call mesa_error(__FILE__,__LINE__,'i_equL')
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% equ(i_equL,k)
         end if

         call save_eqn_residual_info( &
            s, k, nvar, i_equL, resid, 'do1_dlnT_dm_eqn', ierr)

      end subroutine do1_dlnT_dm_eqn


      subroutine do1_constant_L_eqn(s, k, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type(star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: resid
         real(dp) :: scale

         ierr = 0
         scale = max(1d0, abs(s% L_center), abs(s% L_start(k)))
         if (k < s% nz) scale = max(scale, abs(s% L_start(k + 1)))
         resid = (wrap_L_00(s, k) - wrap_L_p1(s, k))/scale
         s% equ(s% i_equL, k) = resid%val
         call save_eqn_residual_info( &
            s, k, nvar, s% i_equL, resid, 'do1_constant_L_eqn', ierr)
      end subroutine do1_constant_L_eqn



      subroutine set_RSP_Lsurf_BC(s, nvar, ierr)
         use const_def, only: crad, clight, pi4
         use eos_def
         use star_utils, only: save_eqn_residual_info, get_area_info_opt_time_center
         use auto_diff_support
         implicit none

         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         integer, intent(in) :: nvar

         type(auto_diff_real_star_order1) :: L1_ad, r1_ad, area_ad, rhs_ad, lhs_ad, resid_ad, inv_R2
         type(auto_diff_real_star_order1) :: T_surf, Erad_ad
         integer :: i_equL
         real(dp) :: factor, scale, L_theta
         logical :: debug

         ierr = 0
         debug = .false.

         i_equL = s% i_equL

         if (s%nz < 1) then
            write(*,*) 'ERROR: Insufficient zones (nz < 1)'
            ierr = -1
            return
         end if

         if (debug) write(*,*) 'RSP zone 1 surface BC being set'

         call get_area_info_opt_time_center(s, 1, area_ad, inv_R2, ierr)
         ! no time centering the surface equations.
         L1_ad = wrap_L_00(s, 1)
         T_surf = wrap_T_00(s,1)

         if (debug) then
            write(*,*) 'T_surf =', T_surf%val, ' r_surf =', r1_ad%val, ' area =', area_ad%val
         end if

         ! rsp equation, zone 1
         rhs_ad = s%RSP2_Lsurf_factor * area_ad * clight * (crad * pow4(T_surf)) ! missing Lc at the moment, so only radiative surface

         if (debug) then
            write(*,*) 'RSP_Lsurf_factor =', s%RSP2_Lsurf_factor
            write(*,*) 'rhs_ad (RSP BC) =', rhs_ad%val
         end if

         ! residual
         lhs_ad = L1_ad
         resid_ad = lhs_ad - rhs_ad

         scale =maxval(s% L_start(1:s% nz))
         resid_ad = resid_ad / scale

         if (debug) then
            write(*,*) 'lhs (L1) =', lhs_ad%val
            write(*,*) 'scaled residual =', resid_ad%val
         end if

         s%equ(i_equL,1) = resid_ad%val

         if (is_bad(resid_ad%val)) then
            write(*,*) 'ERROR: NaN or Inf residual:', resid_ad%val
            ierr = -1
         end if

      call save_eqn_residual_info( &
         s, 1, nvar, i_equL, resid_ad, 'do1_dlnT_dm_eqn', ierr)


      end subroutine set_RSP_Lsurf_BC

      ! only used for dlnT_dm equation
      subroutine eval_dlnPdm_qhse(s, k, &  ! calculate the expected dlnPdm for HSE
            dlnPdm_qhse, Ppoint, ierr)
         use hydro_momentum, only: expected_HSE_grav_term
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dlnPdm_qhse, Ppoint
         integer, intent(out) :: ierr

         real(dp) :: alfa, P_theta
         type(auto_diff_real_star_order1) :: grav, area, P00, Pm1, inv_R2, mlt_Ptrb00, mlt_Ptrbm1, mlt_Ptrb_face
         type(auto_diff_real_star_order1) :: T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, opacity_face
         include 'formats'

         ierr = 0

         ! basic eqn is dP/dm = -G m / (4 pi r^4)
         ! divide by Ppoint to make it unitless

         ! for rotation, multiply gravity by factor fp.  MESA 2, eqn 22.
         call expected_HSE_grav_term(s, k, grav, area, ierr) ! note that expected_HSE_grav_term is negative

         if (s% using_velocity_time_centering .and. &
               s% include_P_in_velocity_time_centering .and. &
               s% lnT(k)/ln10 <= s% max_logT_for_include_P_and_L_in_velocity_time_centering) then
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
            else ! no mlt_pturb
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
               Pm1 = Pm1 + mlt_Ptrbm1 ! include mlt Ptrb in k-1
               P00 = P00 + mlt_Ptrb00 ! include mlt Ptrb in k
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

      end module hydro_temperature
