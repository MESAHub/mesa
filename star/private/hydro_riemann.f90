! ***********************************************************************
!
!   Copyright (C) 2015-2019  Bill Paxton & The MESA Team
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
         do_uface_and_Pface, get_G
         ! Riemann energy eqn is now part of the standard energy equation
         ! Riemann dlnR_dt rqn is now part of the standard radius equation

      contains


      subroutine do_surf_Riemann_dudt_eqn( &
            s, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            skip_partials, nvar, ierr)
         type (star_info), pointer :: s         
         real(dp), intent(in) :: P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         include 'formats'
         call do1_dudt_eqn( &
            s, 1, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            skip_partials, nvar, ierr)
      end subroutine do_surf_Riemann_dudt_eqn
      

      subroutine do1_Riemann_momentum_eqn( &
            s, k, P_surf, skip_partials, nvar, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp), intent(in) :: P_surf ! only used if k==1
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         call do1_dudt_eqn( &
            s, k, P_surf, 0d0, 0d0, 0d0, 0d0, &
            skip_partials, nvar, ierr)
      end subroutine do1_Riemann_momentum_eqn
         

      subroutine do1_dudt_eqn( &
            s, k, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            skip_partials, nvar, ierr)
         use auto_diff_support
         use accurate_sum_auto_diff_18var_order1
         use star_utils, only: get_area_info, store_partials
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp), intent(in) :: P_surf ! only used if k==1
         real(dp), intent(in) :: &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
      
         integer :: j, nz, i_du_dt, i_u
         type(auto_diff_real_18var_order1) :: &
            flux_in_18, flux_out_18, diffusion_source_18, &
            geometry_source_18, gravity_source_18, &
            dudt_expected_18, dudt_actual_18, resid_18
         type(accurate_auto_diff_real_18var_order1) :: sum_18       
         real(dp) :: dt, dm, ie_plus_ke, scal, &
            d_momflux00_dw00, d_momfluxp1_dwp1, residual, &
            area_00, d_area_00_dlnR, inv_R2_00, d_inv_R2_00_dlnR, &
            area_p1, d_area_p1_dlnR, inv_R2_p1, d_inv_R2_p1_dlnR
         real(dp), dimension(nvar) :: d_dm1, d_d00, d_dp1
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
         d_dm1 = 0d0; d_d00 = 0d0; d_dp1 = 0d0
             
         call get_area_info(s, k, & ! using_velocity_time_centering
            area_00, d_area_00_dlnR, inv_R2_00, d_inv_R2_00_dlnR, ierr)
         if (ierr /= 0) return
         if (k < nz) then
            call get_area_info(s, k+1, &
               area_p1, d_area_p1_dlnR, inv_R2_p1, d_inv_R2_p1_dlnR, ierr)
            if (ierr /= 0) return
         end if

         call setup_momentum_flux
         call setup_geometry_source
         call setup_gravity_source
         call setup_diffusion_source

         sum_18 = flux_in_18 - flux_out_18 + &
            geometry_source_18 + gravity_source_18 + diffusion_source_18
         dudt_expected_18 = sum_18
         dudt_expected_18 = dudt_expected_18/dm
         
         ! make residual units be relative difference in energy
         ie_plus_ke = s% energy_start(k) + 0.5d0*s% u_start(k)*s% u_start(k)
         scal = dt*max(abs(s% u_start(k)),s% csound_start(k))/ie_plus_ke
         if (k == 1) scal = scal*1d-2
         
         dudt_actual_18 = 0d0
         dudt_actual_18%val = s% du_dt(k)
         dudt_actual_18%d1Array(i_v_00) = 1d0/dt
         
         resid_18 = scal*(dudt_expected_18 - dudt_actual_18)
         residual = resid_18%val
         s% equ(i_du_dt, k) = residual
         s% u_residual(k) = residual
         
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
         
         if (skip_partials) return
         call unpack_res18(resid_18)
         call store_partials(s, k, i_du_dt, nvar, d_dm1, d_d00, d_dp1)

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = d_d00(s% solver_test_partials_var)  
            write(*,*) 'do1_dudt_eqn', s% solver_test_partials_var
         end if
         
         contains
         
         subroutine setup_momentum_flux
            if (k == 1) then
               call eval_surf_momentum_flux_18(s, P_surf, &
                  dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
                  flux_out_18, d_momflux00_dw00,ierr)   
               if (ierr /= 0) return
            else
               call eval1_momentum_flux_18(s, k, &
                  flux_out_18, d_momflux00_dw00,ierr)   
               if (ierr /= 0) return
            end if
            flux_out_18%d1Array(i_L_00) = d_momflux00_dw00
            if (k < nz) then
               call eval1_momentum_flux_18(s, k+1, &
                  flux_in_18, d_momfluxp1_dwp1, ierr)   
               if (ierr /= 0) return
               flux_in_18 = shift_p1(flux_in_18)
               flux_in_18%d1Array(i_L_p1) = d_momfluxp1_dwp1
            else
               flux_in_18 = 0d0
            end if                  
         end subroutine setup_momentum_flux

         subroutine setup_geometry_source
            type(auto_diff_real_18var_order1) :: P, areaR, areaL
            P = wrap_p_00(s,k)
            areaR = 0d0
            areaR%val = area_00
            areaR%d1Array(i_lnR_00) = d_area_00_dlnR
            if (k == nz) then 
               ! no flux in from left, so only have geometry source on right
               ! this matters for cases with R_center > 0.
               geometry_source_18 = P*areaR
            else
               areaL = 0d0
               areaL%val = area_p1
               areaL%d1Array(i_lnR_p1) = d_area_p1_dlnR
               geometry_source_18 = P*(areaR - areaL)
            end if
         end subroutine setup_geometry_source
         
         subroutine setup_gravity_source
            type(auto_diff_real_18var_order1) :: G, inv_R2, gsR, gsL
            real(dp) :: G00, dG00_dlnR, dG00_dw, Gp1, dGp1_dlnR, dGp1_dw, mR, mL
            ! left 1/2 of dm gets gravity force at left face
            ! right 1/2 of dm gets gravity force at right face.
            ! this form is to match the gravity force equilibrium reconstruction.
            mR = s% m(k)
            if (k == nz) then
               mL = s% M_center
            else
               mL = s% m(k+1)
            end if
            call get_G(s, k, G00, dG00_dlnR, dG00_dw)
            if (dG00_dw /= 0d0) stop 'need to fix dG00_dw for riemann hydro, setup_gravity_source'
            G = 0d0
            G%val = G00
            G%d1Array(i_lnR_00) = dG00_dlnR
            G%d1Array(i_L_00) = dG00_dw
            inv_R2 = 0d0
            inv_R2%val = inv_R2_00
            inv_R2%d1Array(i_lnR_00) = d_inv_R2_00_dlnR
            gsR = -G*mR*0.5d0*dm*inv_R2
            if (k == nz) then
               gsL = 0d0
            else
               call get_G(s, k+1, Gp1, dGp1_dlnR, dGp1_dw)
               G = 0d0
               G%val = Gp1
               G%d1Array(i_lnR_p1) = dGp1_dlnR
               G%d1Array(i_L_p1) = dGp1_dw
               inv_R2 = 0d0
               inv_R2%val = inv_R2_p1
               inv_R2%d1Array(i_lnR_p1) = d_inv_R2_p1_dlnR
               gsL = -G*mL*0.5d0*dm*inv_R2
            end if
            gravity_source_18 = gsL + gsR ! total gravitational force on cell

         end subroutine setup_gravity_source

         subroutine setup_diffusion_source
            type(auto_diff_real_18var_order1) :: u_m1, u_00, u_p1
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
               diffusion_source_18 = sig00*(u_m1 - u_00) - sigp1*(u_00 - u_p1)
            else
               diffusion_source_18 = 0d0
            end if
            s% dudt_RTI(k) = diffusion_source_18%val/dm
         end subroutine setup_diffusion_source
         
         subroutine unpack_res18(res18)
            use star_utils, only: unpack_res18_partials
            type(auto_diff_real_18var_order1) :: res18
            include 'formats'
            call unpack_res18_partials(s, k, nvar, i_du_dt, &
               res18, d_dm1, d_d00, d_dp1)
            if (s% w_div_wc_flag) then
               call e00(s, i_du_dt, s% i_w_div_wc, k, nvar, d_d00(s% i_lum))
               d_d00(s% i_lum) = 0d0
               if (k < nz) then
                  call ep1(s, s% i_du_dt, s% i_w_div_wc, k, nvar, d_dp1(s% i_lum))
                  d_dp1(s% i_lum) = 0d0
               end if
            end if
         end subroutine unpack_res18
         
      end subroutine do1_dudt_eqn

      
      ! (g cm/s)/s from cell(k) to cell(k-1)
      subroutine eval1_momentum_flux_18(s, k, momflux_18, d_momflux_dw, ierr)
         use star_utils, only: get_area_info
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(out) :: momflux_18
         real(dp), intent(out) :: d_momflux_dw
         integer, intent(out) :: ierr
         real(dp) :: area, d_area_dlnR, inv_R2, d_inv_R2_dlnR
         integer :: nz
         logical :: test_partials
         include 'formats'
         
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         ierr = 0
         d_momflux_dw = 0 
         nz = s% nz
         if (k > nz) then
            momflux_18 = 0d0
            return    
         end if
         
         call get_area_info(s, k, &
            area, d_area_dlnR, inv_R2, d_inv_R2_dlnR, ierr)
         if (ierr /= 0) return
         
         momflux_18 = area*s% P_face_18(k)     
         momflux_18%d1Array(i_lnR_00) = &
            momflux_18%d1Array(i_lnR_00) + d_area_dlnR*s% P_face_18(k)%val  
         d_momflux_dw = area*s% d_Pface_dw(k)
              
         if (test_partials) then
            s% solver_test_partials_val = momflux_18%val
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = momflux_18%d1Array(i_lnR_00)
            write(*,*) 'eval1_momentum_flux_18', s% solver_test_partials_var
         end if
         
      end subroutine eval1_momentum_flux_18
      
      
      subroutine eval_surf_momentum_flux_18(s, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            momflux_18, d_momflux_dw, ierr)
         type (star_info), pointer :: s 
         real(dp), intent(in) :: P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
         type(auto_diff_real_18var_order1), intent(out) :: momflux_18
         real(dp), intent(out) :: d_momflux_dw
         integer, intent(out) :: ierr
         integer :: k
         real(qp) :: r, A, momflux
         include 'formats'
         ierr = 0
         k = 1
         r = s% r(k)
         A = 4d0*pi*r*r
         momflux = A*P_surf
         momflux_18%val = momflux
         momflux_18%d1Array(:) = 0d0
         momflux_18%d1Array(i_lnR_00) = momflux*(2 + dlnPsurf_dlnR)
         momflux_18%d1Array(i_L_00) = momflux*dlnPsurf_dL
         momflux_18%d1Array(i_lnd_00) = momflux*dlnPsurf_dlnd
         momflux_18%d1Array(i_lnT_00) = momflux*dlnPsurf_dlnT
         d_momflux_dw = 0
      end subroutine eval_surf_momentum_flux_18
      

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


      subroutine get_G(s, k, G, dG_dlnR, dG_dw)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: G, dG_dlnR, dG_dw
         real(dp) :: cgrav, gr_factor, d_gr_factor_dlnR
         cgrav = s% cgrav(k)
         if (s% rotation_flag .and. s% use_gravity_rotation_correction) &
            cgrav = cgrav*s% fp_rot(k)
         G = cgrav
         dG_dlnR = 0d0
         if (s% rotation_flag .and. s% use_gravity_rotation_correction &
               .and. s% w_div_wc_flag) then
            dG_dw = G/s% fp_rot(k)*s% dfp_rot_dw_div_wc(k)
         else
            dG_dw = 0d0
         end if
      end subroutine get_G

      
      subroutine do1_uface_and_Pface(s, k, ierr)
         use eos_def, only: i_gamma1, i_lnfree_e, i_lnPgas
         use auto_diff_support
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         logical :: test_partials

         type(auto_diff_real_18var_order1) :: &
            r_18, A_18, PL_18, PR_18, uL_18, uR_18, rhoL_18, rhoR_18, &
            gamma1L_18, gamma1R_18, csL_18, csR_18, G_18, dPdm_grav_18, &
            Sl1_18, Sl2_18, Sr1_18, Sr2_18, numerator_18, denominator_18, &
            Sl_18, Sr_18, Ss_18, P_face_L_18, P_face_R_18, du_18
         ! use d1Array(i_L_00) for w_00
         real(dp) :: G, dG_dlnR, dG_dw, delta_m, f
            
         include 'formats'
         
         ierr = 0
         test_partials = .false.
         !test_partials = (k == s% solver_test_partials_k)
         
         s% RTI_du_diffusion_kick(k) = 0d0
         s% d_uface_dw(k) = 0
         s% d_Pface_dw(k) = 0
                            
         if (k == 1) then
            s% u_face_18(k) = wrap_u_00(s,k)
            s% P_face_18(k) = wrap_p_00(s,k)
            return            
         end if
      
         r_18 = wrap_r_00(s,k)
         A_18 = 4d0*pi*r_18**2
         
         PL_18 = wrap_p_00(s,k)
         PR_18 = wrap_p_m1(s,k)

         uL_18 = wrap_u_00(s,k)
         uR_18 = wrap_u_m1(s,k)
      
         rhoL_18 = wrap_d_00(s,k)
         rhoR_18 = wrap_d_m1(s,k)
         
         gamma1L_18 = wrap_gamma1_00(s,k)
         gamma1R_18 = wrap_gamma1_m1(s,k)
      
         csL_18 = sqrt(gamma1L_18*PL_18/rhoL_18)
         csR_18 = sqrt(gamma1R_18*PR_18/rhoR_18)
         
         ! change PR and PL for gravity
         call get_G(s, k, G, dG_dlnR, dG_dw)
         G_18 = 0d0
         G_18%val = G
         G_18%d1Array(i_lnR_00) = dG_dlnR
         G_18%d1Array(i_L_00) = dG_dw
         
         dPdm_grav_18 = -G_18*s% m_grav(k)/(r_18**2*A_18)  ! cm^-1 s^-2
         
         delta_m = 0.5d0*s% dm(k) ! positive delta_m from left center to edge
         PL_18 = PL_18 + delta_m*dPdm_grav_18

         delta_m = -0.5d0*s% dm(k-1) ! negative delta_m from right center to edge
         PR_18 = PR_18 + delta_m*dPdm_grav_18
            
         ! acoustic wavespeeds (eqn 2.38)
         Sl1_18 = uL_18 - csL_18
         Sl2_18 = uR_18 - csR_18

         ! take Sl = min(Sl1, Sl2)
         if (Sl1_18%val < Sl2_18%val) then
            Sl_18 = Sl1_18
         else
            Sl_18 = Sl2_18
         end if

         Sr1_18 = uR_18 + csR_18         
         Sr2_18 = uL_18 + csL_18
         
         ! take Sr = max(Sr1, Sr2)
         if (Sr1_18%val > Sr2_18%val) then
            Sr_18 = Sr1_18
         else
            Sr_18 = Sr2_18
         end if
         
         ! contact velocity (eqn 2.20)
         numerator_18 = uR_18*rhoR_18*(Sr_18 - uR_18) + uL_18*rhoL_18*(uL_18 - Sl_18) + (PL_18 - PR_18)         
         denominator_18 = rhoR_18*(Sr_18 - uR_18) + rhoL_18*(uL_18 - Sl_18)         

         if (denominator_18%val == 0d0 .or. is_bad(denominator_18%val)) then
            ierr = -1
            if (s% report_ierr) then
               write(*,2) 'u_face denominator bad', k, denominator_18%val
            end if
            return
         end if
         
         Ss_18 = numerator_18/denominator_18
         
         s% u_face_18(k) = Ss_18
         s% d_uface_dw(k) = s% u_face_18(k)%d1Array(i_L_00)
         s% u_face_18(k)%d1Array(i_L_00) = 0d0

         ! contact pressure (eqn 2.19)
         P_face_L_18 = rhoL_18*(uL_18-Sl_18)*(uL_18-Ss_18) + PL_18         
         P_face_R_18 = rhoR_18*(uR_18-Sr_18)*(uR_18-Ss_18) + PR_18
         
         s% P_face_18(k) = 0.5d0*(P_face_L_18 + P_face_R_18) ! these are ideally equal
         s% d_Pface_dw(k) = s% P_face_18(k)%d1Array(i_L_00)
         s% P_face_18(k)%d1Array(i_L_00) = 0d0

         if (k < s% nz .and. s% RTI_flag) then
             if (s% eta_RTI(k) > 0d0 .and. &
                   s% dlnddt_RTI_diffusion_factor > 0d0 .and. s% dt > 0d0) then
                f = s% dlnddt_RTI_diffusion_factor*s% eta_RTI(k)/s% dm_bar(k)
                du_18 = f*A_18*(rhoL_18 - rhoR_18) ! bump uface in direction of lower density
                s% RTI_du_diffusion_kick(k) = du_18%val
                s% u_face_18(k) = s% u_face_18(k) + du_18
             end if
         end if

         if (s% P_face_start(k) < 0d0) then
            s% u_face_start(k) = s% u_face_18(k)%val
            s% P_face_start(k) = s% P_face_18(k)%val
         end if

         if (test_partials) then
            s% solver_test_partials_val = 0
            s% solver_test_partials_var = 0        
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do1_uface_and_Pface', s% solver_test_partials_var
         end if
     
      end subroutine do1_uface_and_Pface

         
      end module hydro_riemann

