! ***********************************************************************
!
!   Copyright (C) 2015-2019  The MESA Team
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

      module hydro_riemann

      use star_private_def
      use const_def, only: dp, pi
      use star_utils, only: em1, e00, ep1
      use utils_lib
      use auto_diff
      use auto_diff_support

      implicit none

      ! Cheng, J, Shu, C-W, and Zeng, Q.,
      ! "A Conservative Lagrangian Scheme for Solving
      !  Compressible Fluid Flows with Multiple Internal Energy Equations",
      ! Commun. Comput. Phys., 12, pp 1307-1328, 2012.

      ! Cheng, J. and Shu, C-W,
      ! "Positivity-preserving Lagrangian scheme for multi-material
      !  compressible flow", J. Comp. Phys., 257 (2014), 143-168.

      ! Kappeli, R. and Mishra, S.,
      ! "Well-balanced schemes for the Euler equations with gravitation",
      ! J. Comp. Phys., 259 (2014), 199-219.

      private
      public :: do_surf_Riemann_dudt_eqn, do1_Riemann_momentum_eqn, &
         do_uface_and_Pface, get_Riemann_shock_diagnostics, &
         get_RTI_momentum_diffusion, eval_Riemann_dudt_rhs
         ! Riemann energy eqn is now part of the standard energy equation
         ! Riemann dlnR_dt rqn is now part of the standard radius equation

      contains

      subroutine do_surf_Riemann_dudt_eqn(s, P_surf_ad, nvar, ierr)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1), intent(in) :: P_surf_ad
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         call do1_dudt_eqn(s, 1, P_surf_ad, nvar, ierr)
      end subroutine do_surf_Riemann_dudt_eqn


      subroutine do1_Riemann_momentum_eqn(s, k, nvar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: P_surf_ad
         P_surf_ad = 0
         call do1_dudt_eqn(s, k, P_surf_ad, nvar, ierr)
      end subroutine do1_Riemann_momentum_eqn


      subroutine do1_dudt_eqn( &
            s, k, P_surf_ad, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: P_surf_ad  ! only for k=1
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         integer :: i_du_dt
         type(auto_diff_real_star_order1) :: &
            dudt_expected_ad, dudt_actual_ad, resid_ad
         real(dp) :: dt, ie_plus_ke, scal, residual
         logical :: test_partials

         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (s% use_other_momentum) &
            call mesa_error(__FILE__,__LINE__,'Riemann dudt does not support use_other_momentum')
         if (s% use_other_momentum_implicit) &
            call mesa_error(__FILE__,__LINE__,'Riemann dudt does not support use_other_momentum_implicit')
         if (s% use_mass_corrections) &
            call mesa_error(__FILE__,__LINE__,'Riemann dudt does not support use_mass_corrections')

         ierr = 0
         i_du_dt = s% i_du_dt
         dt = s% dt
         call eval_Riemann_dudt_rhs( &
            s, k, P_surf_ad, .true., .true., dudt_expected_ad, ierr)
         if (ierr /= 0) return

         ! make residual units be relative difference in energy
         ie_plus_ke = s% energy_start(k) + 0.5d0*s% u_start(k)*s% u_start(k)
         scal = dt*max(abs(s% u_start(k)),s% csound_start(k))/ie_plus_ke
         if (k == 1) scal = scal*1d-2

         dudt_actual_ad = 0d0
         dudt_actual_ad%val = s% dxh_u(k)/dt
         dudt_actual_ad%d1Array(i_v_00) = 1d0/dt

         resid_ad = scal*(dudt_expected_ad - dudt_actual_ad)
         residual = resid_ad%val
         s% equ(i_du_dt, k) = residual

         if (is_bad(residual)) then
            ierr = -1
            return
!$omp critical (dudt_eqn)
            write(*,2) 'residual', k, residual
            call mesa_error(__FILE__,__LINE__,'do1_dudt_eqn')
!$omp end critical (dudt_eqn)
         end if

         call save_eqn_residual_info(s, k, nvar, i_du_dt, resid_ad, 'do1_dudt_eqn', ierr)

         if (test_partials) then
            s% solver_test_partials_val = resid_ad% val
         end if

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = resid_ad% d1Array(i_lnR_00)
            !write(*,*) 'do1_dudt_eqn', s% solver_test_partials_var
            end if

      end subroutine do1_dudt_eqn


      subroutine eval_Riemann_dudt_rhs( &
            s, k, P_surf_ad, use_time_centering, include_tdc_Uq, &
            dudt_expected_ad, ierr)
         use accurate_sum_auto_diff_star_order1
         use star_utils, only: get_area_info_opt_time_center
         use tdc_hydro, only: compute_tdc_Uq_dm_cell
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: P_surf_ad
         logical, intent(in) :: use_time_centering, include_tdc_Uq
         type(auto_diff_real_star_order1), intent(out) :: dudt_expected_ad
         integer, intent(out) :: ierr
         integer :: nz
         type(auto_diff_real_star_order1) :: &
            flux_in_ad, flux_out_ad, diffusion_source_ad, &
            geometry_source_ad, gravity_source_ad, &
            area_00, area_p1, inv_R2_00, inv_R2_p1, Uq_cell
         type(accurate_auto_diff_real_star_order1) :: sum_ad
         real(dp) :: dm, v_drag, drag_factor, drag_fraction

         ierr = 0
         nz = s% nz
         dm = s% dm(k)

         call get_area_info(k, area_00, inv_R2_00, ierr)
         if (ierr /= 0) return
         if (k < nz) then
            call get_area_info(k + 1, area_p1, inv_R2_p1, ierr)
            if (ierr /= 0) return
            area_p1 = shift_p1(area_p1)
            inv_R2_p1 = shift_p1(inv_R2_p1)
         end if

         call setup_momentum_flux
         call setup_geometry_source(ierr); if (ierr /= 0) return
         call setup_gravity_source
         call setup_diffusion_source

         ! compute_tdc_Uq_dm_cell returns Uq*dm for the force sum.
         Uq_cell = 0d0
         if (include_tdc_Uq .and. &
               s% MLT_option == 'TDC' .and. s% TDC_alpha_M > 0d0) then
            Uq_cell = compute_tdc_Uq_dm_cell(s, k, ierr)
            if (ierr /= 0) return
         end if

         sum_ad = flux_in_ad
         sum_ad = sum_ad - flux_out_ad
         sum_ad = sum_ad + geometry_source_ad
         sum_ad = sum_ad + gravity_source_ad
         sum_ad = sum_ad + diffusion_source_ad
         sum_ad = sum_ad + Uq_cell
         dudt_expected_ad = sum_ad
         dudt_expected_ad = dudt_expected_ad/dm

         drag_factor = s% v_drag_factor
         v_drag = s% v_drag
         if (s% q(k) < s% q_for_v_drag_full_off) then
            drag_fraction = 0d0
         else if (s% q(k) > s% q_for_v_drag_full_on) then
            drag_fraction = 1d0
         else
            drag_fraction = (s% q(k) - s% q_for_v_drag_full_off)&
                               /(s% q_for_v_drag_full_on - s% q_for_v_drag_full_off)
         end if
         drag_factor = drag_factor*drag_fraction

         if (drag_factor > 0d0) then
            if (s% u(k) > v_drag) then
               dudt_expected_ad = dudt_expected_ad - drag_factor*pow2(s% u(k) - v_drag)/s% r(k)
            else if (s% u(k) < -v_drag) then
               dudt_expected_ad = dudt_expected_ad + drag_factor*pow2(s% u(k) + v_drag)/s% r(k)
            end if
         end if

         contains

         subroutine get_area_info(kk, area_ad, inv_R2_ad, ierr)
            integer, intent(in) :: kk
            type(auto_diff_real_star_order1), intent(out) :: area_ad, inv_R2_ad
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: r_ad, r2_ad

            if (use_time_centering) then
               call get_area_info_opt_time_center(s, kk, area_ad, inv_R2_ad, ierr)
               return
            end if

            ierr = 0
            r_ad = wrap_r_00(s, kk)
            r2_ad = pow2(r_ad)
            area_ad = 4d0*pi*r2_ad
            inv_R2_ad = 1d0/r2_ad
         end subroutine get_area_info

         subroutine setup_momentum_flux
            if (k == 1) then
               flux_out_ad = P_surf_ad*area_00
            else
               flux_out_ad = s% P_face_ad(k)*area_00
            end if
            if (k < nz) then
               flux_in_ad = shift_p1(s% P_face_ad(k+1))*area_p1
            else
               flux_in_ad = 0d0
            end if
         end subroutine setup_momentum_flux

         subroutine setup_geometry_source(ierr)
            use star_utils, only: calc_Ptot_ad_tw
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: P
            real(dp), dimension(s% species) :: d_Ptot_dxa
            logical, parameter :: skip_Peos = .false., skip_mlt_Pturb = .false.
            ierr = 0
            ! use same P here as the cell pressure in P_face calculation
            call calc_Ptot_ad_tw( &
               s, k, skip_Peos, skip_mlt_Pturb, P, d_Ptot_dxa, ierr, &
               use_time_centering)
            if (ierr /= 0) return
            if (k == nz) then
               ! no flux in from left, so only have geometry source on right
               ! this matters for cases with R_center > 0.
               geometry_source_ad = P*area_00
            else
               geometry_source_ad = P*(area_00 - area_p1)
            end if
         end subroutine setup_geometry_source

         subroutine setup_gravity_source
            type(auto_diff_real_star_order1) :: G00, Gp1, gsL, gsR
            real(dp) :: mR, mL
            ! left 1/2 of dm gets gravity force at left face
            ! right 1/2 of dm gets gravity force at right face.
            ! this form is to match the gravity force equilibrium reconstruction.
            mR = s% m(k)
            if (k == nz) then
               mL = s% M_center
            else
               mL = s% m(k+1)
            end if
            call get_G(s, k, G00)
            gsR = -G00*mR*0.5d0*dm*inv_R2_00
            if (k == nz) then
               gsL = 0d0
            else
               call get_G(s, k+1, Gp1)
               Gp1 = shift_p1(Gp1)
               gsL = -Gp1*mL*0.5d0*dm*inv_R2_p1
            end if
            gravity_source_ad = gsL + gsR  ! total gravitational force on cell
         end subroutine setup_gravity_source

         subroutine setup_diffusion_source
            type(auto_diff_real_star_order1) :: dissipation_ad
            call get_RTI_momentum_diffusion(s, k, diffusion_source_ad, dissipation_ad)
            s% dudt_RTI(k) = diffusion_source_ad%val/dm
         end subroutine setup_diffusion_source

      end subroutine eval_Riemann_dudt_rhs


      subroutine get_RTI_momentum_diffusion(s, k, force_ad, dissipation_ad)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: force_ad, dissipation_ad

         real(dp) :: sig00, sigp1
         type(auto_diff_real_star_order1) :: &
            u_m1, u_00, u_p1, ubar_m1, ubar_00, ubar_p1, &
            du00, dup1, dubar00, dubarp1

         force_ad = 0d0
         dissipation_ad = 0d0
         if (.not. s% RTI_flag .or. s% dudt_RTI_diffusion_factor <= 0d0) return

         u_m1 = 0d0
         u_p1 = 0d0
         ubar_m1 = 0d0
         ubar_p1 = 0d0
         sig00 = 0d0
         sigp1 = 0d0

         u_00 = wrap_u_00(s,k)
         ubar_00 = 0.5d0*(u_00 + s% u_start(k))
         if (k > 1) then
            u_m1 = wrap_u_m1(s,k)
            ubar_m1 = 0.5d0*(u_m1 + s% u_start(k-1))
            sig00 = s% dudt_RTI_diffusion_factor*s% sig_RTI(k)
         end if
         if (k < s% nz) then
            u_p1 = wrap_u_p1(s,k)
            ubar_p1 = 0.5d0*(u_p1 + s% u_start(k+1))
            sigp1 = s% dudt_RTI_diffusion_factor*s% sig_RTI(k+1)
         end if

         du00 = u_m1 - u_00
         dup1 = u_00 - u_p1
         dubar00 = ubar_m1 - ubar_00
         dubarp1 = ubar_00 - ubar_p1
         force_ad = sig00*du00 - sigp1*dup1
         ! Share the kinetic energy dissipated at each interface between its cells.
         dissipation_ad = 0.5d0*(sig00*du00*dubar00 + sigp1*dup1*dubarp1)
      end subroutine get_RTI_momentum_diffusion


      subroutine do_uface_and_Pface( &
            s, ierr, include_rsp2_Uq, use_time_centering)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         logical, intent(in), optional :: include_rsp2_Uq
         logical, intent(in), optional :: use_time_centering
         integer :: k, op_err
         logical :: include_Uq, time_center
         include 'formats'
         ierr = 0
         include_Uq = .true.
         if (present(include_rsp2_Uq)) include_Uq = include_rsp2_Uq
         time_center = .true.
         if (present(use_time_centering)) time_center = use_time_centering
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k = 1, s% nz
            op_err = 0
            call do1_uface_and_Pface(s, k, include_Uq, time_center, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
      end subroutine do_uface_and_Pface


      subroutine get_G(s, k, G)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: G
         real(dp) :: cgrav
         cgrav = s% cgrav(k)
         G = cgrav
         if (s% rotation_flag .and. s% use_gravity_rotation_correction) &
            G = G*s% fp_rot(k)
      end subroutine get_G


      subroutine get_Riemann_shock_diagnostics( &
            s, k, compression, pressure_jump, shock_strength, D_mix_factor, ierr)
         use math_lib, only: pow3
         use star_utils, only: calc_Ptot_ad_tw
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: &
            compression, pressure_jump, shock_strength, D_mix_factor
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: &
            r_ad, area_ad, PL_ad, PR_ad, G_ad, dPdm_grav_ad, &
            csL_ad, csR_ad
         real(dp), dimension(s% species) :: d_Ptot_dxa
         real(dp) :: cs_face, P_face, P_min, onset, full_on, reduction, x
         real(dp) :: delta_m, uL, uR, rhoL, rhoR, Sl, Sr, Ss
         real(dp) :: numerator, denominator, P_face_L, P_face_R
         logical, parameter :: skip_Peos = .false., skip_mlt_Pturb = .false.

         ierr = 0
         compression = 0d0
         pressure_jump = 0d0
         shock_strength = 0d0
         D_mix_factor = 1d0
         if (.not. s% u_flag .or. k <= 1 .or. k > s% nz) return

         call calc_Ptot_ad_tw(s, k, skip_Peos, skip_mlt_Pturb, &
            PL_ad, d_Ptot_dxa, ierr)
         if (ierr /= 0) return
         call calc_Ptot_ad_tw(s, k-1, skip_Peos, skip_mlt_Pturb, &
            PR_ad, d_Ptot_dxa, ierr)
         if (ierr /= 0) return
         PR_ad = shift_m1(PR_ad)

         if (PL_ad%val <= 0d0 .or. PR_ad%val <= 0d0) return
         csL_ad = sqrt(wrap_gamma1_00(s,k)*PL_ad/wrap_d_00(s,k))
         csR_ad = sqrt(wrap_gamma1_m1(s,k)*PR_ad/wrap_d_m1(s,k))
         cs_face = 0.5d0*(csL_ad%val + csR_ad%val)
         if (cs_face <= 0d0 .or. is_bad(cs_face)) return
         uL = s% u(k)
         uR = s% u(k-1)
         rhoL = s% rho(k)
         rhoR = s% rho(k-1)
         Sl = min(uL - csL_ad%val, uR - csR_ad%val)
         Sr = max(uR + csR_ad%val, uL + csL_ad%val)

         r_ad = wrap_r_00(s,k)
         area_ad = 4d0*pi*pow2(r_ad)
         call get_G(s, k, G_ad)
         dPdm_grav_ad = -G_ad*s% m_grav(k)/(pow2(r_ad)*area_ad)

         delta_m = 0.5d0*s% dm(k)
         PL_ad = PL_ad + delta_m*dPdm_grav_ad
         delta_m = -0.5d0*s% dm(k-1)
         PR_ad = PR_ad + delta_m*dPdm_grav_ad

         if (PL_ad%val <= 0d0 .or. PR_ad%val <= 0d0) return

         numerator = uR*rhoR*(Sr-uR) + uL*rhoL*(uL-Sl) + &
            (PL_ad%val - PR_ad%val)
         denominator = rhoR*(Sr-uR) + rhoL*(uL-Sl)
         if (denominator == 0d0 .or. is_bad(denominator)) return
         Ss = numerator/denominator

         P_face_L = rhoL*(uL-Sl)*(uL-Ss) + PL_ad%val
         P_face_R = rhoR*(uR-Sr)*(uR-Ss) + PR_ad%val
         P_face = 0.5d0*(P_face_L + P_face_R)
         if (P_face <= 0d0 .or. is_bad(P_face)) return

         compression = max(0d0, (uL-uR)/cs_face)
         P_min = min(PL_ad%val, PR_ad%val)
         pressure_jump = max(0d0, P_face/P_min - 1d0)
         shock_strength = min(compression, pressure_jump)

         onset = s% Riemann_shock_D_mix_reduction_on
         full_on = s% Riemann_shock_D_mix_reduction_full_on
         if (full_on <= onset .or. shock_strength <= onset) return

         x = min(1d0, (shock_strength-onset)/(full_on-onset))
         reduction = pow3(x)*(10d0 + x*(-15d0 + 6d0*x))
         D_mix_factor = 1d0 - reduction

      end subroutine get_Riemann_shock_diagnostics


      subroutine do1_uface_and_Pface( &
            s, k, include_rsp2_Uq, use_time_centering, ierr)
         use eos_def, only: i_gamma1, i_lnfree_e, i_lnPgas
         use star_utils, only: calc_Ptot_ad_tw, get_face_weights
         use hydro_rsp2, only: compute_Uq_face
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in) :: include_rsp2_Uq, use_time_centering
         integer, intent(out) :: ierr
         logical :: test_partials

         type(auto_diff_real_star_order1) :: &
            r_ad, A_ad, PL_ad, PR_ad, uL_ad, uR_ad, rhoL_ad, rhoR_ad, &
            gamma1L_ad, gamma1R_ad, csL_ad, csR_ad, G_ad, dPdm_grav_ad, &
            Sl1_ad, Sl2_ad, Sr1_ad, Sr2_ad, numerator_ad, denominator_ad, &
            Sl_ad, Sr_ad, Ss_ad, P_face_L_ad, P_face_R_ad, du_ad, Uq_ad
         real(dp), dimension(s% species) :: d_Ptot_dxa  ! skip this
         logical, parameter :: skip_Peos = .false., skip_mlt_Pturb = .false.
         real(dp) :: delta_m, f

         include 'formats'

         ierr = 0
         test_partials = .false.
         !test_partials = (k == s% solver_test_partials_k)

         s% RTI_du_diffusion_kick(k) = 0d0
         s% d_uface_domega(k) = 0

         if (k == 1) then
            s% u_face_ad(k) = wrap_u_00(s,k)
            s% P_face_ad(k) = wrap_Peos_00(s,k)
            return
         end if

         r_ad = wrap_r_00(s,k)
         A_ad = 4d0*pi*pow2(r_ad)

         call calc_Ptot_ad_tw( &
            s, k, skip_Peos, skip_mlt_Pturb, PL_ad, d_Ptot_dxa, ierr, &
            use_time_centering)
         if (ierr /= 0) return
         call calc_Ptot_ad_tw( &
            s, k - 1, skip_Peos, skip_mlt_Pturb, PR_ad, d_Ptot_dxa, ierr, &
            use_time_centering)
         if (ierr /= 0) return
         PR_ad = shift_m1(PR_ad)

         uL_ad = wrap_u_00(s,k)
         uR_ad = wrap_u_m1(s,k)

         rhoL_ad = wrap_d_00(s,k)
         rhoR_ad = wrap_d_m1(s,k)

         gamma1L_ad = wrap_gamma1_00(s,k)
         gamma1R_ad = wrap_gamma1_m1(s,k)

         csL_ad = sqrt(gamma1L_ad*PL_ad/rhoL_ad)
         csR_ad = sqrt(gamma1R_ad*PR_ad/rhoR_ad)

         ! change PR and PL for gravity
         call get_G(s, k, G_ad)

         dPdm_grav_ad = -G_ad*s% m_grav(k)/(pow2(r_ad)*A_ad)  ! cm^-1 s^-2

         delta_m = 0.5d0*s% dm(k)  ! positive delta_m from left center to edge
         PL_ad = PL_ad + delta_m*dPdm_grav_ad

         delta_m = -0.5d0*s% dm(k-1)  ! negative delta_m from right center to edge
         PR_ad = PR_ad + delta_m*dPdm_grav_ad

         ! acoustic wavespeeds (eqn 2.38)
         Sl1_ad = uL_ad - csL_ad
         Sl2_ad = uR_ad - csR_ad

         ! take Sl = min(Sl1, Sl2)
         if (Sl1_ad%val < Sl2_ad%val) then
            Sl_ad = Sl1_ad
         else
            Sl_ad = Sl2_ad
         end if

         Sr1_ad = uR_ad + csR_ad
         Sr2_ad = uL_ad + csL_ad

         ! take Sr = max(Sr1, Sr2)
         if (Sr1_ad%val > Sr2_ad%val) then
            Sr_ad = Sr1_ad
         else
            Sr_ad = Sr2_ad
         end if

         ! contact velocity (eqn 2.20)
         numerator_ad = uR_ad*rhoR_ad*(Sr_ad - uR_ad) + uL_ad*rhoL_ad*(uL_ad - Sl_ad) + (PL_ad - PR_ad)
         denominator_ad = rhoR_ad*(Sr_ad - uR_ad) + rhoL_ad*(uL_ad - Sl_ad)

         if (denominator_ad%val == 0d0 .or. is_bad(denominator_ad%val)) then
            ierr = -1
            if (s% report_ierr) then
               write(*,2) 'u_face denominator bad', k, denominator_ad%val
            end if
            return
         end if

         Ss_ad = numerator_ad/denominator_ad

         s% u_face_ad(k) = Ss_ad
         s% d_uface_domega(k) = s% u_face_ad(k)%d1Array(i_L_00)

         ! contact pressure (eqn 2.19)
         P_face_L_ad = rhoL_ad*(uL_ad-Sl_ad)*(uL_ad-Ss_ad) + PL_ad
         P_face_R_ad = rhoR_ad*(uR_ad-Sr_ad)*(uR_ad-Ss_ad) + PR_ad

         s% P_face_ad(k) = 0.5d0*(P_face_L_ad + P_face_R_ad)  ! these are ideally equal

         if (k < s% nz .and. s% RTI_flag) then
             if (s% eta_RTI(k) > 0d0 .and. &
                   s% dlnddt_RTI_diffusion_factor > 0d0 .and. s% dt > 0d0) then
                f = s% dlnddt_RTI_diffusion_factor*s% eta_RTI(k)/s% dm_bar(k)
                du_ad = f*A_ad*(rhoL_ad - rhoR_ad)  ! bump uface in direction of lower density
                s% RTI_du_diffusion_kick(k) = du_ad%val
                s% u_face_ad(k) = s% u_face_ad(k) + du_ad
             end if
         end if


         ! RSP2 currently applies Uq to the reconstructed face velocity.
         if (s% RSP2_flag .and. include_rsp2_Uq) then
            Uq_ad = compute_Uq_face(s, k, ierr)
            if (ierr /= 0) return
            s% u_face_ad(k) = s% u_face_ad(k) + Uq_ad
         end if

         s% u_face_val(k) = s% u_face_ad(k)%val

         if (s% P_face_start(k) < 0d0) then
            s% u_face_start(k) = s% u_face_val(k)
            s% P_face_start(k) = s% P_face_ad(k)%val
         end if

         if (test_partials) then
            s% solver_test_partials_val = PL_ad% val
            s% solver_test_partials_var = s% i_w_div_wc
            s% solver_test_partials_dval_dx = PL_ad% d1Array(i_w_div_wc_00)
            write(*,*) 'do1_uface_and_Pface', s% solver_test_partials_var, PL_ad% val
         end if

      end subroutine do1_uface_and_Pface

      end module hydro_riemann
