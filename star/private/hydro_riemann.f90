! ***********************************************************************
!
!   Copyright (C) 2015-2019  The MESA Team
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

      module hydro_riemann
      
      use star_private_def
      use const_def
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
         do_uface_and_Pface
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
         use accurate_sum_auto_diff_star_order1
         use star_utils, only: get_area_info_opt_time_center, save_eqn_residual_info
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: P_surf_ad ! only for k=1
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
      
         integer :: j, nz, i_du_dt, i_u
         type(auto_diff_real_star_order1) :: &
            flux_in_ad, flux_out_ad, diffusion_source_ad, &
            geometry_source_ad, gravity_source_ad, &
            area_00, area_p1, inv_R2_00, inv_R2_p1, &
            dudt_expected_ad, dudt_actual_ad, resid_ad
         type(accurate_auto_diff_real_star_order1) :: sum_ad       
         real(dp) :: dt, dm, ie_plus_ke, scal, residual
         logical :: dbg, do_diffusion, test_partials

         include 'formats'
         dbg = .false.

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         if (s% use_other_momentum) &
            stop 'Riemann dudt does not support use_other_momentum'
         if (s% use_other_momentum_implicit) &
            stop 'Riemann dudt does not support use_other_momentum_implicit'
         if (s% use_mass_corrections) &
            stop 'Riemann dudt does not support use_mass_corrections'
            
         ierr = 0
         nz = s% nz
         i_du_dt = s% i_du_dt
         dt = s% dt
         dm = s% dm(k)     
             
         call get_area_info_opt_time_center(s, k, area_00, inv_R2_00, ierr)
         if (ierr /= 0) return
         if (k < nz) then
            call get_area_info_opt_time_center(s, k+1, area_p1, inv_R2_p1, ierr)
            if (ierr /= 0) return
            area_p1 = shift_p1(area_p1)
            inv_R2_p1 = shift_p1(inv_R2_p1)
         end if

         call setup_momentum_flux
         call setup_geometry_source(ierr); if (ierr /= 0) return
         call setup_gravity_source
         call setup_diffusion_source

         sum_ad = flux_in_ad - flux_out_ad + &
            geometry_source_ad + gravity_source_ad + diffusion_source_ad
         dudt_expected_ad = sum_ad
         dudt_expected_ad = dudt_expected_ad/dm
         
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
            stop 'do1_dudt_eqn'
!$omp end critical (dudt_eqn)
         end if

         if (test_partials) then
            s% solver_test_partials_val = residual
         end if
         
         call save_eqn_residual_info(s, k, nvar, i_du_dt, resid_ad, 'do1_dudt_eqn', ierr)

         if (test_partials) then
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do1_dudt_eqn', s% solver_test_partials_var
         end if
         
         contains
         
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
            call calc_Ptot_ad_tw(s, k, skip_Peos, skip_mlt_Pturb, P, d_Ptot_dxa, ierr)
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
            real(dp) :: mR, mL, dG00_dw, dGp1_dw
            ! left 1/2 of dm gets gravity force at left face
            ! right 1/2 of dm gets gravity force at right face.
            ! this form is to match the gravity force equilibrium reconstruction.
            mR = s% m(k)
            if (k == nz) then
               mL = s% M_center
            else
               mL = s% m(k+1)
            end if
            call get_G(s, k, G00, dG00_dw)
            G00%d1Array(i_L_00) = dG00_dw
            gsR = -G00*mR*0.5d0*dm*inv_R2_00
            if (k == nz) then
               gsL = 0d0
            else
               call get_G(s, k+1, Gp1, dGp1_dw)
               Gp1 = shift_p1(Gp1)
               Gp1%d1Array(i_L_p1) = dGp1_dw
               gsL = -Gp1*mL*0.5d0*dm*inv_R2_p1
            end if
            gravity_source_ad = gsL + gsR ! total gravitational force on cell

         end subroutine setup_gravity_source

         subroutine setup_diffusion_source
            type(auto_diff_real_star_order1) :: u_m1, u_00, u_p1
            real(dp) :: sig00, sigp1
            do_diffusion = s% RTI_flag .and. s% dudt_RTI_diffusion_factor > 0d0
            if (do_diffusion) then ! add diffusion source term to dudt
               u_p1 = 0d0 ! sets val and d1Array to 0
               if (k < nz) then
                  sigp1 = s% dudt_RTI_diffusion_factor*s% sig_RTI(k+1)
                  u_p1%val = s% u(k+1)
                  u_p1%d1Array(i_v_p1) = 1d0
               else
                  sigp1 = 0
               end if
               u_m1 = 0d0 ! sets val and d1Array to 0
               if (k > 1) then
                  sig00 = s% dudt_RTI_diffusion_factor*s% sig_RTI(k)
                  u_m1%val = s% u(k-1)
                  u_m1%d1Array(i_v_m1) = 1d0
               else
                  sig00 = 0
               end if
               u_00 = 0d0 ! sets val and d1Array to 0
               u_00%val = s% u(k)
               u_00%d1Array(i_v_00) = 1d0
               diffusion_source_ad = sig00*(u_m1 - u_00) - sigp1*(u_00 - u_p1)
            else
               diffusion_source_ad = 0d0
            end if
            s% dudt_RTI(k) = diffusion_source_ad%val/dm
         end subroutine setup_diffusion_source
         
      end subroutine do1_dudt_eqn
      

      subroutine do_uface_and_Pface(s, ierr)
         type (star_info), pointer :: s 
         integer, intent(out) :: ierr
         integer :: k, op_err
         include 'formats'
         ierr = 0
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k = 1, s% nz
            op_err = 0
            call do1_uface_and_Pface(s, k, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
      end subroutine do_uface_and_Pface


      subroutine get_G(s, k, G, dG_dw_div_wc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: G
         real(dp), intent(out) :: dG_dw_div_wc
         real(dp) :: cgrav
         cgrav = s% cgrav(k)
         if (s% rotation_flag .and. s% use_gravity_rotation_correction) &
            cgrav = cgrav*s% fp_rot(k)
         G = cgrav
         if (s% rotation_flag .and. s% use_gravity_rotation_correction &
               .and. s% w_div_wc_flag) then
            dG_dw_div_wc = G%val/s% fp_rot(k)*s% dfp_rot_dw_div_wc(k)
         else
            dG_dw_div_wc = 0d0
         end if
      end subroutine get_G

      
      subroutine do1_uface_and_Pface(s, k, ierr)
         use eos_def, only: i_gamma1, i_lnfree_e, i_lnPgas
         use star_utils, only: calc_Ptot_ad_tw, get_face_weights
         use hydro_rsp2, only: compute_Uq_face
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         logical :: test_partials

         type(auto_diff_real_star_order1) :: &
            r_ad, A_ad, PL_ad, PR_ad, uL_ad, uR_ad, rhoL_ad, rhoR_ad, &
            gamma1L_ad, gamma1R_ad, csL_ad, csR_ad, G_ad, dPdm_grav_ad, &
            Sl1_ad, Sl2_ad, Sr1_ad, Sr2_ad, numerator_ad, denominator_ad, &
            Sl_ad, Sr_ad, Ss_ad, P_face_L_ad, P_face_R_ad, du_ad, Uq_ad
         real(dp), dimension(s% species) :: d_Ptot_dxa ! skip this
         logical, parameter :: skip_Peos = .false., skip_mlt_Pturb = .false.
         real(dp) :: dG_dw_div_wc, delta_m, f
            
         include 'formats'
         
         ierr = 0
         test_partials = .false.
         !test_partials = (k == s% solver_test_partials_k)
         
         s% RTI_du_diffusion_kick(k) = 0d0
         s% d_uface_domega(k) = 0
         s% d_Pface_domega(k) = 0
                            
         if (k == 1) then
            s% u_face_ad(k) = wrap_u_00(s,k)
            s% P_face_ad(k) = wrap_Peos_00(s,k)
            return            
         end if
      
         r_ad = wrap_r_00(s,k)
         A_ad = 4d0*pi*pow2(r_ad)
         
         call calc_Ptot_ad_tw(s, k, skip_Peos, skip_mlt_Pturb, PL_ad, d_Ptot_dxa, ierr)
         if (ierr /= 0) return
         call calc_Ptot_ad_tw(s, k-1, skip_Peos, skip_mlt_Pturb, PR_ad, d_Ptot_dxa, ierr)
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
         call get_G(s, k, G_ad, dG_dw_div_wc)
         G_ad%d1Array(i_L_00) = dG_dw_div_wc
         
         dPdm_grav_ad = -G_ad*s% m_grav(k)/(pow2(r_ad)*A_ad)  ! cm^-1 s^-2
         
         delta_m = 0.5d0*s% dm(k) ! positive delta_m from left center to edge
         PL_ad = PL_ad + delta_m*dPdm_grav_ad

         delta_m = -0.5d0*s% dm(k-1) ! negative delta_m from right center to edge
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
         s% u_face_ad(k)%d1Array(i_L_00) = 0d0

         ! contact pressure (eqn 2.19)
         P_face_L_ad = rhoL_ad*(uL_ad-Sl_ad)*(uL_ad-Ss_ad) + PL_ad         
         P_face_R_ad = rhoR_ad*(uR_ad-Sr_ad)*(uR_ad-Ss_ad) + PR_ad
         
         s% P_face_ad(k) = 0.5d0*(P_face_L_ad + P_face_R_ad) ! these are ideally equal
         s% d_Pface_domega(k) = s% P_face_ad(k)%d1Array(i_L_00)
         s% P_face_ad(k)%d1Array(i_L_00) = 0d0

         if (k < s% nz .and. s% RTI_flag) then
             if (s% eta_RTI(k) > 0d0 .and. &
                   s% dlnddt_RTI_diffusion_factor > 0d0 .and. s% dt > 0d0) then
                f = s% dlnddt_RTI_diffusion_factor*s% eta_RTI(k)/s% dm_bar(k)
                du_ad = f*A_ad*(rhoL_ad - rhoR_ad) ! bump uface in direction of lower density
                s% RTI_du_diffusion_kick(k) = du_ad%val
                s% u_face_ad(k) = s% u_face_ad(k) + du_ad
             end if
         end if
         
         if (s% using_RSP2) then ! include Uq in u_face
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
            s% solver_test_partials_val = 0
            s% solver_test_partials_var = 0        
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do1_uface_and_Pface', s% solver_test_partials_var
         end if
     
      end subroutine do1_uface_and_Pface

         
      end module hydro_riemann

