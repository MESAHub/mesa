! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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


      module hydro_momentum

      use star_private_def
      use const_def
      use utils_lib, only: mesa_error, is_bad
      use auto_diff
      use star_utils, only: em1, e00, ep1

      implicit none

      private
      public :: do1_momentum_eqn, do_surf_momentum_eqn, do1_radius_eqn


      contains

      
      subroutine do_surf_momentum_eqn(s, P_surf_ad, skip_partials, nvar, ierr)
         use star_utils, only: store_partials
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1), intent(in) :: P_surf_ad
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         real(dp) :: d_dm1(nvar), d_d00(nvar), d_dp1(nvar)
         include 'formats'
         ierr = 0
         call get1_momentum_eqn( &
            s, 1, P_surf_ad, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'ierr /= 0 for do_surf_momentum_eqn'
            return
         end if         
         if (skip_partials) return
         call store_partials(s, 1, s% i_dv_dt, nvar, d_dm1, d_d00, d_dp1)
      end subroutine do_surf_momentum_eqn

      
      subroutine do1_momentum_eqn(s, k, skip_partials, nvar, ierr)
         use star_utils, only: store_partials
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         real(dp) :: d_dm1(nvar), d_d00(nvar), d_dp1(nvar)
         type(auto_diff_real_star_order1) :: P_surf_ad ! only used if k == 1
         include 'formats'
         P_surf_ad = 0d0
         call get1_momentum_eqn( &
            s, k, P_surf_ad, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'ierr /= 0 for get1_momentum_eqn', k
            return
         end if         
         if (skip_partials) return
         call store_partials(s, k, s% i_dv_dt, nvar, d_dm1, d_d00, d_dp1)
      end subroutine do1_momentum_eqn


      subroutine get1_momentum_eqn( &
            s, k, P_surf_ad, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         use chem_def, only: chem_isos
         use accurate_sum_auto_diff_star_order1
         use auto_diff_support

         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: P_surf_ad ! only used if k == 1
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         real(dp), intent(out) :: d_dm1(nvar), d_d00(nvar), d_dp1(nvar)
         integer, intent(out) :: ierr
         
         real(dp) :: residual, dm_face, d_grav_dw, dXP, iXPavg, dm_div_A
         real(dp), dimension(s% species) :: &
            d_dXP_dxam1, d_dXP_dxa00, d_iXPavg_dxam1, d_iXPavg_dxa00, &
            d_residual_dxam1, d_residual_dxa00
         integer :: nz, j, i_dv_dt, i_lnd, i_lnT, i_lnR, i_lum, i_v
         logical :: test_partials
         
         type(auto_diff_real_star_order1) :: resid1_ad, resid_ad, &
            other_ad, dm_div_A_ad, grav_ad, area_ad, dXP_ad, d_mlt_Pturb_ad, &
            iXPavg_ad, other_dm_div_A_ad, grav_dm_div_A_ad, &
            RTI_terms_ad, RTI_terms_dm_div_A_ad
         type(accurate_auto_diff_real_star_order1) :: residual_sum_ad

         include 'formats'
         
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         ierr = 0
         call init

!   dv/dt = - G*m/r^2 - (dXP_ad + d_mlt_Pturb_ad)*area/dm + extra_grav + Uq + RTI_diffusion + RTI_kick
! 
!   grav_ad = expected_HSE_grav_term = -G*m/r^2 with possible modifications for rotation
!   other_ad = expected_non_HSE_term = extra_grav - dv/dt + Uq
!   extra_grav is from the other_momentum hook
!   dXP_ad = pressure difference across face from center to center of adjacent cells (excluding mlt_Pturb effects)
!        XP = P_ad + avQ_ad + Pt_ad + extra_pressure, with time centering
!   iXPavg_ad = 1/(avg XP).  for normalizing equation
!   d_mlt_Pturb_ad = difference in MLT convective pressure across face
!   RTI_terms_ad = RTI_diffusion + RTI_kick
!   dm_div_A_ad = dm/area
! 
!   0  = extra_grav - dv/dt + Uq - G*m/r^2 - RTI_diffusion - RTI_kick - (dXP_ad + d_mlt_Pturb_ad)*area/dm
!   0  = other + grav - RTI_terms - (dXP_ad + d_mlt_Pturb_ad)*area/dm
!   0  = (other + grav - RTI_terms)*dm/area - dXP_ad - d_mlt_Pturb_ad
!   0  = other_dm_div_A_ad + grav_dm_div_A_ad - dXP_ad - d_mlt_Pturb_ad + RTI_terms_dm_div_A_ad

         call setup_HSE(d_grav_dw, dm_div_A, ierr); if (ierr /= 0) return ! grav_ad and dm_div_A_ad
         call setup_non_HSE(ierr); if (ierr /= 0) return ! other = s% extra_grav(k) - s% dv_dt(k)
         call setup_dXP(ierr); if (ierr /= 0) return ! dXP_ad, iXPavg_ad
         call setup_d_mlt_Pturb(ierr); if (ierr /= 0) return ! d_mlt_Pturb_ad
         call setup_RTI_terms(ierr); if (ierr /= 0) return ! RTI_terms_ad
         
         other_dm_div_A_ad = other_ad*dm_div_A_ad
         grav_dm_div_A_ad = grav_ad*dm_div_A_ad
         RTI_terms_dm_div_A_ad = RTI_terms_ad*dm_div_A_ad
         
         if (.false.) then
            if (is_bad(other_dm_div_A_ad%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 other_dm_div_A_ad', k, other_dm_div_A_ad%d1Array(i_lnd_m1)
               stop 'get1_momentum_eqn'
            end if
            if (is_bad(grav_dm_div_A_ad%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 grav_dm_div_A_ad', k, grav_dm_div_A_ad%d1Array(i_lnd_m1)
               stop 'get1_momentum_eqn'
            end if
            if (is_bad(dXP_ad%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 dXP_ad', k, dXP_ad%d1Array(i_lnd_m1)
               stop 'get1_momentum_eqn'
            end if
            if (is_bad(d_mlt_Pturb_ad%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 d_mlt_Pturb_ad', k, d_mlt_Pturb_ad%d1Array(i_lnd_m1)
               stop 'get1_momentum_eqn'
            end if
         end if
         
         ! sum terms in residual_sum_ad using accurate_auto_diff_real_star_order1
         residual_sum_ad = other_dm_div_A_ad + grav_dm_div_A_ad - dXP_ad - d_mlt_Pturb_ad + RTI_terms_dm_div_A_ad
         
         resid1_ad = residual_sum_ad ! convert back to auto_diff_real_star_order1
         resid_ad = resid1_ad*iXPavg_ad ! scaling
         residual = resid_ad%val
         s% equ(i_dv_dt, k) = residual      
         s% v_residual(k) = residual

         if (is_bad(residual)) then
!$omp critical (hydro_momentum_crit1)
            write(*,2) 'momentum eqn residual', k, residual
            stop 'get1_momentum_eqn'
!$omp end critical (hydro_momentum_crit1)
         end if
         if (test_partials) then
            s% solver_test_partials_val = residual
         end if
         if (skip_partials) return
         call unpack_res18(resid_ad)

         if (test_partials) then
            s% solver_test_partials_var = i_lnR
            s% solver_test_partials_dval_dx = d_d00(s% solver_test_partials_var)
            write(*,*) 'get1_momentum_eqn', s% solver_test_partials_var
         end if
         
         contains
         
         subroutine init
            i_dv_dt = s% i_dv_dt
            i_lnd = s% i_lnd
            i_lnT = s% i_lnT
            i_lnR = s% i_lnR
            i_lum = s% i_lum
            i_v = s% i_v
            nz = s% nz
            if (k > 1) then
               dm_face = (s% dm(k) + s% dm(k-1))/2d0
            else ! k == 1
               dm_face = s% dm(k)/2d0
            end if
            d_dm1 = 0d0; d_d00 = 0d0; d_dp1 = 0d0
         end subroutine init
         
         subroutine setup_HSE(d_grav_dw, dm_div_A, ierr)
            real(dp), intent(out) :: d_grav_dw, dm_div_A
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            call expected_HSE_grav_term(s, k, grav_ad, d_grav_dw, area_ad, ierr)
            if (ierr /= 0) return
            dm_div_A_ad = dm_face/area_ad
            dm_div_A = dm_div_A_ad%val
         end subroutine setup_HSE
         
         subroutine setup_non_HSE(ierr)
            integer, intent(out) :: ierr
            real(dp) :: other
            include 'formats'
            ierr = 0
            ! other = extra_grav - dv/dt
            call expected_non_HSE_term(s, k, other_ad, other, ierr)
         end subroutine setup_non_HSE

         subroutine setup_dXP(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            ! dXP = pressure difference across face from center to center of adjacent cells.
            ! iXPavg = average pressure at face for normalization of the equation to something like dlnP/dm
            call get_dXP_face_info(s, k, P_surf_ad, &
               dXP_ad, dXP, d_dXP_dxam1, d_dXP_dxa00, &
               iXPavg_ad, iXPavg, d_iXPavg_dxam1, d_iXPavg_dxa00, ierr)
            if (ierr /= 0) return
         end subroutine setup_dXP
                  
         subroutine setup_d_mlt_Pturb(ierr)
            integer, intent(out) :: ierr
            real(dp) :: d_mlt_Pturb, d_dmltPturb_dlndm1, d_dmltPturb_dlnd00
            ierr = 0

            ! d_mlt_Pturb = difference in MLT convective pressure across face
            if (s% mlt_Pturb_factor > 0d0 .and. s% mlt_vc_start(k) > 0d0 .and. k > 1) then
               d_mlt_Pturb = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*(s% rho(k-1) - s% rho(k))/3d0
               d_dmltPturb_dlndm1 = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*s% rho(k-1)/3d0
               d_dmltPturb_dlnd00 = -s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*s% rho(k)/3d0
            else
               d_mlt_Pturb = 0d0
               d_dmltPturb_dlndm1 = 0d0
               d_dmltPturb_dlnd00 = 0d0
            end if
            d_mlt_Pturb_ad = 0d0
            d_mlt_Pturb_ad%val = d_mlt_Pturb
            d_mlt_Pturb_ad%d1Array(i_lnd_m1) = d_dmltPturb_dlndm1
            d_mlt_Pturb_ad%d1Array(i_lnd_00) = d_dmltPturb_dlnd00
         end subroutine setup_d_mlt_Pturb         
                  
         subroutine setup_RTI_terms(ierr)
            use auto_diff_support
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: v_p1, v_00, v_m1, dvdt_diffusion, &
               f, rho_00, rho_m1, dvdt_kick
            real(dp) :: sigm1, sig00
            ierr = 0
            RTI_terms_ad = 0d0
            if (.not. s% RTI_flag) return
            if (k >= s% nz .or. k <= 1) return
            ! diffusion of specific momentum (i.e. v)
            if (s% dudt_RTI_diffusion_factor > 0d0) then ! add diffusion source term to dvdt
               ! sigmid_RTI(k) is mixing flow at center k in (gm sec^1)
               sigm1 = s% dudt_RTI_diffusion_factor*s% sigmid_RTI(k-1)
               sig00 = s% dudt_RTI_diffusion_factor*s% sigmid_RTI(k)
               v_p1 = wrap_v_p1(s, k)
               v_00 = wrap_v_00(s, k)
               v_m1 = wrap_v_m1(s, k)
               dvdt_diffusion = sig00*(v_p1 - v_00) - sigm1*(v_00 - v_m1) ! (g/s)*(cm/s)
               dvdt_diffusion = dvdt_diffusion/s% dm_bar(k) ! divide by g to get units of cm/s^2
            else
               dvdt_diffusion = 0d0
            end if
            ! kick to adjust densities
            if (s% eta_RTI(k) > 0d0 .and. &
               s% dlnddt_RTI_diffusion_factor > 0d0 .and. s% dt > 0d0) then
               f = s% dlnddt_RTI_diffusion_factor*s% eta_RTI(k)/dm_div_A_ad
               rho_00 = wrap_d_00(s, k)
               rho_m1 = wrap_d_m1(s, k)
               dvdt_kick = f*(rho_00 - rho_m1)/s% dt ! change v according to direction of lower density
            else
               dvdt_kick = 0d0
            end if            
            RTI_terms_ad = dvdt_diffusion + dvdt_kick            
         end subroutine setup_RTI_terms
         
         subroutine unpack_res18(res18)
            use star_utils, only: unpack_res18_partials
            type(auto_diff_real_star_order1) :: res18
            real(dp) :: resid1
            integer :: j
            include 'formats'
            call unpack_res18_partials(s, k, nvar, i_dv_dt, &
               res18, d_dm1, d_d00, d_dp1)
            if (s% rotation_flag .and. s% w_div_wc_flag .and. s% use_gravity_rotation_correction) then
               call e00(s, i_dv_dt, s% i_w_div_wc, k, nvar, iXPavg*d_grav_dw*dm_div_A)            
            end if            
            ! do partials wrt composition   
            resid1 = resid1_ad%val         
            do j=1,s% species
               d_residual_dxa00(j) = resid1*d_iXPavg_dxa00(j) - iXPavg*d_dXP_dxa00(j)
               call e00(s, i_dv_dt, j+s% nvar_hydro, k, nvar, d_residual_dxa00(j))
            end do
            if (k > 1) then 
               do j=1,s% species
                  d_residual_dxam1(j) = resid1*d_iXPavg_dxam1(j) - iXPavg*d_dXP_dxam1(j)
                  call em1(s, i_dv_dt, j+s% nvar_hydro, k, nvar, d_residual_dxam1(j))
               end do
            end if            
         end subroutine unpack_res18
            
      end subroutine get1_momentum_eqn
      
      
      ! returns -G*m/r^2 with possible modifications for rotation.  MESA 2, eqn 22.
      subroutine expected_HSE_grav_term(s, k, grav, d_grav_dw_div_wc, area, ierr)
         use star_utils, only: get_area_info
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: area, grav
         real(dp), intent(out) :: d_grav_dw_div_wc
         integer, intent(out) :: ierr
         
         type(auto_diff_real_star_order1) :: inv_R2
         logical :: test_partials

         include 'formats'
         ierr = 0
         
         call get_area_info(s, k, area, inv_R2, ierr)
         if (ierr /= 0) return

         grav = -s% cgrav(k)*s% m_grav(k)*inv_R2
         
         d_grav_dw_div_wc = 0d0
         if (s% rotation_flag .and. s% use_gravity_rotation_correction) then
            if (s% w_div_wc_flag) d_grav_dw_div_wc = grav%val*s% dfp_rot_dw_div_wc(k)
            grav = grav*s% fp_rot(k) 
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
      
      
      ! other = s% extra_grav(k) - s% dv_dt(k)
      subroutine expected_non_HSE_term(s, k, other_ad, other, ierr)
         use hydro_tdc, only: compute_Uq_face
         use accurate_sum_auto_diff_star_order1
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: other_ad
         real(dp), intent(out) :: other
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            extra_ad, accel_ad, v_00, Uq_ad
         real(dp) :: accel, d_accel_dv, fraction_on
         logical :: test_partials, local_v_flag

         include 'formats'

         ierr = 0
         
         extra_ad = 0d0
         if (s% use_other_momentum .or. s% use_other_momentum_implicit) then
            if (s% use_other_momentum_implicit) then
               call wrap(extra_ad, s% extra_grav(k), &
                  s% d_extra_grav_dlndm1(k), s% d_extra_grav_dlnd00(k), 0d0, &
                  s% d_extra_grav_dlnTm1(k), s% d_extra_grav_dlnT00(k), 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, s% d_extra_grav_dlnR(k), 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, s% d_extra_grav_dL(k), 0d0)
            else
               extra_ad%val = s% extra_grav(k)
            end if
         end if
         
         accel_ad = 0d0
         if (s% v_flag) then
            
            if (s% i_lnT == 0) then
               local_v_flag = .true.
            else
               local_v_flag = &
                  (s% xh_old(s% i_lnT,k)/ln10 >= s% velocity_logT_lower_bound)
            end if

            if (local_v_flag) then
               accel = s% dv_dt(k)
               d_accel_dv = s% dVARDOT_dVAR
            else ! assume vstart(k) = 0 and
               ! constant acceleration dv_dt so vfinal(k) = dv_dt*dt
               ! v(k) = dr/dt = average velocity =
               !      (vstart + vfinal)/2 = dv_dt*dt/2 when vstart = 0
               ! so (1/2)*dv_dt*dt = v(k)
               accel = 2d0*s% v(k)*s% dVARDOT_dVAR
               d_accel_dv = 2d0*s% dVARDOT_dVAR
            end if
            accel_ad%val = accel
            accel_ad%d1Array(i_v_00) = d_accel_dv
         
         end if ! v_flag

         Uq_ad = 0d0
         if (s% TDC_flag) then ! Uq(k) is turbulent viscosity drag at face k
            Uq_ad = compute_Uq_face(s, k, ierr)
            if (ierr /= 0) return
         end if
         
         other_ad = extra_ad - accel_ad + Uq_ad
         other = other_ad%val
         
         if (.false.) then
            if (is_bad(extra_ad%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 extra_ad', k, extra_ad%d1Array(i_lnd_m1)
               stop 'expected_non_HSE_term'
            end if
            if (is_bad(accel_ad%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 accel_ad', k, accel_ad%d1Array(i_lnd_m1)
               stop 'expected_non_HSE_term'
            end if
            if (is_bad(Uq_ad%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 Uq_ad', k, Uq_ad%d1Array(i_lnd_m1)
               stop 'expected_non_HSE_term'
            end if
         end if
         
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         if (test_partials) then
            s% solver_test_partials_val = other
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = 0d0
            write(*,*) 'expected_non_HSE_term', s% solver_test_partials_var
         end if
        
      end subroutine expected_non_HSE_term

      ! dXP = pressure difference across face from center to center of adjacent cells.
      ! excluding mlt_Pturb effects
      subroutine get_dXP_face_info(s, k, P_surf_ad, &
            dXP_ad, dXP, d_dXP_dxam1, d_dXP_dxa00, &
            iXPavg_ad, iXPavg, d_iXPavg_dxam1, d_iXPavg_dxa00, ierr)
         use star_utils, only: calc_XP_ad_tw
         use accurate_sum_auto_diff_star_order1
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: P_surf_ad ! only used if k == 1
         type(auto_diff_real_star_order1), intent(out) :: dXP_ad, iXPavg_ad
         real(dp), intent(out) :: dXP, iXPavg
         real(dp), intent(out), dimension(s% species) :: &
            d_dXP_dxam1, d_dXP_dxa00, d_iXPavg_dxam1, d_iXPavg_dxa00
         integer, intent(out) :: ierr
         
         real(dp) :: XPm1, XP00, XPavg, alfa, beta
         real(dp), dimension(s% species) :: &
            d_XPm1_dxam1, d_XP00_dxa00, d_XPavg_dxam1, d_XPavg_dxa00
         type(auto_diff_real_star_order1) :: &
            XP00_ad, XPm1_ad, XPavg_ad
         integer :: j
         logical, parameter :: skip_P = .false., skip_mlt_Pturb = .true.
         logical :: test_partials

         include 'formats'

         ierr = 0
         
         call calc_XP_ad_tw( &
            s, k, skip_P, skip_mlt_Pturb, XP00_ad, d_XP00_dxa00, ierr)
         if (ierr /= 0) return
         XP00 = XP00_ad%val
            
         if (k > 1) then
            call calc_XP_ad_tw( &
               s, k-1, skip_P, skip_mlt_Pturb, XPm1_ad, d_XPm1_dxam1, ierr)
            if (ierr /= 0) return
            XPm1_ad = shift_m1(XPm1_ad)
         else ! k == 1
            XPm1_ad = P_surf_ad
         end if
         XPm1 = XPm1_ad%val
            
         dXP_ad = XPm1_ad - XP00_ad
         dXP = XPm1 - XP00
         do j=1,s% species
            d_dXP_dxam1(j) = d_XPm1_dxam1(j)
            d_dXP_dxa00(j) = -d_XP00_dxa00(j)
         end do

         if (k == 1) then
            XPavg_ad = XP00_ad
            do j=1,s% species
               d_XPavg_dxam1(j) = 0d0  
               d_XPavg_dxa00(j) = d_XP00_dxa00(j)
            end do
         else
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1d0 - alfa
            XPavg_ad = alfa*XP00_ad + beta*XPm1_ad
            do j=1,s% species
               d_XPavg_dxam1(j) = beta*d_XPm1_dxam1(j)
               d_XPavg_dxa00(j) = alfa*d_XP00_dxa00(j)
            end do
         end if
         XPavg = XPavg_ad%val
         
         iXPavg_ad = 1d0/XPavg_ad
         iXPavg = 1d0/XPavg         
         do j=1,s% species
            d_iXPavg_dxam1(j) = -iXPavg*d_XPavg_dxam1(j)/XPavg   
            d_iXPavg_dxa00(j) = -iXPavg*d_XPavg_dxa00(j)/XPavg
         end do

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (test_partials) then
            s% solver_test_partials_val = XP00
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = 0d0
            write(*,*) 'get_dXP_face_info', s% solver_test_partials_var
         end if
        
      end subroutine get_dXP_face_info      


      subroutine do1_radius_eqn( &
            s, k, skip_partials, nvar, ierr)
         use auto_diff_support, only: unwrap
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: uc_ad
         real(dp) :: dt, r, r0, r_div_r0, cs, v_expected, v_factor, residual, &
            dr_div_r0_actual, dr_div_r0_expected, uc_factor, unused, &
            d_uface_dlnR, d_uface_du00, d_uface_dum1, d_dlnR00, d_dv00, &
            d_uface_dlnd00, d_uface_dlndm1,  d_uface_dlnT00, d_uface_dlnTm1, &
            residual_old, d_dlnR00_old, d_dv00_old
         integer :: nz, i_dlnR_dt, i_v, i_u, i_lnR, i_w_div_wc
         logical :: test_partials, force_zero_v

         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         ierr = 0
         dt = s% dt
         nz = s% nz
         i_dlnR_dt = s% i_dlnR_dt
         i_v = s% i_v
         i_u = s% i_u
         i_lnR = s% i_lnR
         i_w_div_wc = s% i_w_div_wc
         
         if (i_v == 0 .and. i_u == 0) stop 'must have either v or u for do1_radius_eqn'

         r = s% r(k)
         r0 = s% r_start(k)
         r_div_r0 = r/r0

         force_zero_v = (s% q(k) > s% velocity_q_upper_bound)
         
         if (s% i_lnT /= 0 .and. .not. force_zero_v) &
            force_zero_v = &
               (s% xh_old(s% i_lnT,k)/ln10 < s% velocity_logT_lower_bound .and. &
                  s% dt < secyer*s% max_dt_yrs_for_velocity_logT_lower_bound)
                  
         if (force_zero_v) then
            cs = s% csound_start(k)
            if (i_v /= 0) then
               s% equ(i_dlnR_dt, k) = s% v(k)/cs ! this makes v(k) => 0
               s% lnR_residual(k) = s% equ(i_dlnR_dt, k)
               if (skip_partials) return
               call e00(s, i_dlnR_dt, i_v, k, nvar, 1d0/cs)
               return
            else if (i_u /= 0) then
               s% equ(i_dlnR_dt, k) = s% u(k)/cs ! this makes u(k) => 0
               s% lnR_residual(k) = s% equ(i_dlnR_dt, k)
               if (skip_partials) return
               call e00(s, i_dlnR_dt, i_u, k, nvar, 1d0/cs)
               return
            end if
         end if

         if (i_u /= 0) then
            if (s% using_velocity_time_centering) then
               uc_ad = 0.5d0*(s% u_face_ad(k) + s% u_face_start(k))
            else
               uc_ad = s% u_face_ad(k)
            end if
            call unwrap(uc_ad, v_expected, &
               d_uface_dlndm1, d_uface_dlnd00, unused, &
               d_uface_dlnTm1, d_uface_dlnT00, unused, &
               unused, unused, unused, &
               unused, d_uface_dlnR, unused, &
               d_uface_dum1, d_uface_du00, unused, &
               unused, unused, unused)
         else
            v_expected = s% vc(k)
         end if
         v_factor = s% d_vc_dv

         ! dr = r - r0 = v_expected*dt
         ! eqn: dr/r0 = v_expected*dt/r0
         ! (r - r0)/r0 = r/r0 - 1 = exp(lnR)/exp(lnR0) - 1
         ! = exp(lnR - lnR0) - 1 = exp(dlnR) - 1 = exp(dlnR_dt*dt) - 1
         ! eqn becomes: v_expected*dt/r0 = expm1(dlnR)
         dr_div_r0_actual = expm1(s% dxh_lnR(k)) ! expm1(x) = E^x - 1
         dr_div_r0_expected = v_expected*dt/r0
         
         residual = dr_div_r0_expected - dr_div_r0_actual
         s% equ(i_dlnR_dt, k) = residual
         s% lnR_residual(k) = residual

         if (test_partials) then
            s% solver_test_partials_val = residual
         end if

         if (skip_partials) return

         ! partials of dr_div_r0_expected
         if (i_v /= 0) then            
            call e00(s, i_dlnR_dt, i_v, k, nvar, v_factor*dt/r0)            
         else if (i_u /= 0) then
            uc_factor = v_factor*dt/r0
            call e00(s, i_dlnR_dt, i_lnR, k, nvar, uc_factor*d_uface_dlnR)
            call e00(s, i_dlnR_dt, i_u, k, nvar, uc_factor*d_uface_du00)         
            call e00(s, i_dlnR_dt, s% i_lnd, k, nvar, uc_factor*d_uface_dlnd00)
            if (s% do_struct_thermo) &
               call e00(s, i_dlnR_dt, s% i_lnT, k, nvar, uc_factor*d_uface_dlnT00)         
            if (k > 1) then
               call em1(s, i_dlnR_dt, i_u, k, nvar, uc_factor*d_uface_dum1)            
               call em1(s, i_dlnR_dt, s% i_lnd, k, nvar, uc_factor*d_uface_dlndm1)
               if (s% do_struct_thermo) &
                  call em1(s, i_dlnR_dt, s% i_lnT, k, nvar, uc_factor*d_uface_dlnTm1)
            end if         
            if (s% w_div_wc_flag) then
               call e00(s, i_dlnR_dt, i_w_div_wc, k, nvar, uc_factor*s% d_uface_domega(k))
            end if
         end if

         ! partial of -dr_div_r0_actual wrt lnR
         call e00(s, i_dlnR_dt, i_lnR, k, nvar, -r_div_r0) 

         if (test_partials) then   
            s% solver_test_partials_var = i_lnR
            s% solver_test_partials_dval_dx = -r_div_r0
            write(*,*) 'do1_radius_eqn', s% solver_test_partials_var
         end if

      end subroutine do1_radius_eqn


      end module hydro_momentum

