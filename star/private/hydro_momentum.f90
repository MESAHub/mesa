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
      public :: do1_momentum_eqn, do_surf_momentum_eqn


      contains

      
      subroutine do_surf_momentum_eqn( &
            s, P_surf_18, xscale, equ, skip_partials, nvar, ierr)
         use star_utils, only: store_partials
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1), intent(in) :: P_surf_18
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         real(dp) :: d_dm1(nvar), d_d00(nvar), d_dp1(nvar)
         include 'formats'
         ierr = 0
         call get1_momentum_eqn( &
            s, 1, P_surf_18, xscale, equ, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'ierr /= 0 for do_surf_momentum_eqn'
            return
         end if         
         if (skip_partials) return
         call store_partials(s, 1, xscale, s% i_dv_dt, nvar, d_dm1, d_d00, d_dp1)
      end subroutine do_surf_momentum_eqn

      
      subroutine do1_momentum_eqn(s, k, xscale, equ, skip_partials, nvar, ierr)
         use star_utils, only: store_partials
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         real(dp) :: d_dm1(nvar), d_d00(nvar), d_dp1(nvar)
         type(auto_diff_real_18var_order1) :: P_surf_18 ! only used if k == 1
         include 'formats'
         P_surf_18 = 0d0
         call get1_momentum_eqn( &
            s, k, P_surf_18, xscale, equ, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'ierr /= 0 for get1_momentum_eqn', k
            return
         end if         
         if (skip_partials) return
         call store_partials(s, k, xscale, s% i_dv_dt, nvar, d_dm1, d_d00, d_dp1)
      end subroutine do1_momentum_eqn


      subroutine get1_momentum_eqn( &
            s, k, P_surf_18, xscale, equ, skip_partials, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         use chem_def, only: chem_isos
         use accurate_sum_auto_diff_18var_order1
         use auto_diff_support

         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(in) :: P_surf_18 ! only used if k == 1
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
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
         
         type(auto_diff_real_18var_order1) :: resid1_18, resid_18, &
            other_18, dm_div_A_18, grav_18, area_18, dXP_18, d_mlt_Pturb_18, &
            iXPavg_18, other_dm_div_A_18, grav_dm_div_A_18, &
            RTI_terms_18, RTI_terms_dm_div_A_18
         type(accurate_auto_diff_real_18var_order1) :: residual_sum_18

         include 'formats'
         
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         ierr = 0
         call init

!   dv/dt = - G*m/r^2 - (dXP_18 + d_mlt_Pturb_18)*area/dm + extra_grav + Uq + RTI_diffusion + RTI_kick
! 
!   grav_18 = expected_HSE_grav_term = -G*m/r^2 with possible modifications for rotation
!   other_18 = expected_non_HSE_term = extra_grav - dv/dt + Uq
!   extra_grav is from the other_momentum hook
!   dXP_18 = pressure difference across face from center to center of adjacent cells (excluding mlt_Pturb effects)
!        XP = P_18 + avQ_18 + Pt_18 + extra_pressure, with time centering
!   iXPavg_18 = 1/(avg XP).  for normalizing equation
!   d_mlt_Pturb_18 = difference in MLT convective pressure across face
!   RTI_terms_18 = RTI_diffusion + RTI_kick
!   dm_div_A_18 = dm/area
! 
!   0  = extra_grav - dv/dt + Uq - G*m/r^2 - RTI_diffusion - RTI_kick - (dXP_18 + d_mlt_Pturb_18)*area/dm
!   0  = other + grav - RTI_terms - (dXP_18 + d_mlt_Pturb_18)*area/dm
!   0  = (other + grav - RTI_terms)*dm/area - dXP_18 - d_mlt_Pturb_18
!   0  = other_dm_div_A_18 + grav_dm_div_A_18 - dXP_18 - d_mlt_Pturb_18 + RTI_terms_dm_div_A_18

         call setup_HSE(d_grav_dw, dm_div_A, ierr); if (ierr /= 0) return ! grav_18 and dm_div_A_18
         call setup_non_HSE(ierr); if (ierr /= 0) return ! other = s% extra_grav(k) - s% dv_dt(k)
         call setup_dXP(ierr); if (ierr /= 0) return ! dXP_18, iXPavg_18
         call setup_d_mlt_Pturb(ierr); if (ierr /= 0) return ! d_mlt_Pturb_18
         call setup_RTI_terms(ierr); if (ierr /= 0) return ! RTI_terms_18
         
         other_dm_div_A_18 = other_18*dm_div_A_18
         grav_dm_div_A_18 = grav_18*dm_div_A_18
         RTI_terms_dm_div_A_18 = RTI_terms_18*dm_div_A_18
         
         if (.false.) then
            if (is_bad(other_dm_div_A_18%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 other_dm_div_A_18', k, other_dm_div_A_18%d1Array(i_lnd_m1)
               stop 'get1_momentum_eqn'
            end if
            if (is_bad(grav_dm_div_A_18%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 grav_dm_div_A_18', k, grav_dm_div_A_18%d1Array(i_lnd_m1)
               stop 'get1_momentum_eqn'
            end if
            if (is_bad(dXP_18%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 dXP_18', k, dXP_18%d1Array(i_lnd_m1)
               stop 'get1_momentum_eqn'
            end if
            if (is_bad(d_mlt_Pturb_18%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 d_mlt_Pturb_18', k, d_mlt_Pturb_18%d1Array(i_lnd_m1)
               stop 'get1_momentum_eqn'
            end if
         end if
         
         ! sum terms in residual_sum_18 using accurate_auto_diff_real_18var_order1
         residual_sum_18 = other_dm_div_A_18 + grav_dm_div_A_18 - dXP_18 - d_mlt_Pturb_18 + RTI_terms_dm_div_A_18
         
         resid1_18 = residual_sum_18 ! convert back to auto_diff_real_18var_order1
         resid_18 = resid1_18*iXPavg_18 ! scaling
         residual = resid_18%val
         equ(i_dv_dt, k) = residual      
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
         call unpack_res18(resid_18)

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
         
         subroutine setup_HSE(d_grav_dw, dm_div_A,ierr)
            real(dp), intent(out) :: d_grav_dw, dm_div_A
            integer, intent(out) :: ierr
            real(dp) :: grav, d_grav_dlnR, area, d_area_dlnR
            include 'formats'
            ierr = 0
            call expected_HSE_grav_term(s, k, &
               grav, d_grav_dlnR, d_grav_dw, area, d_area_dlnR, ierr)
            if (ierr /= 0) return
            grav_18 = 0d0
            grav_18%val = grav
            grav_18%d1Array(i_lnR_00) = d_grav_dlnR
            area_18 = 0d0
            area_18%val = area
            area_18%d1Array(i_lnR_00) = d_area_dlnR
            dm_div_A_18 = dm_face/area_18
            dm_div_A = dm_div_A_18%val
         end subroutine setup_HSE
         
         subroutine setup_non_HSE(ierr)
            integer, intent(out) :: ierr
            real(dp) :: other
            include 'formats'
            ierr = 0
            ! other = extra_grav - dv/dt
            call expected_non_HSE_term(s, k, other_18, other, ierr)
         end subroutine setup_non_HSE

         subroutine setup_dXP(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            ! dXP = pressure difference across face from center to center of adjacent cells.
            ! iXPavg = average pressure at face for normalization of the equation to something like dlnP/dm
            call get_dXP_face_info(s, k, P_surf_18, &
               dXP_18, dXP, d_dXP_dxam1, d_dXP_dxa00, &
               iXPavg_18, iXPavg, d_iXPavg_dxam1, d_iXPavg_dxa00, ierr)
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
            d_mlt_Pturb_18 = 0d0
            d_mlt_Pturb_18%val = d_mlt_Pturb
            d_mlt_Pturb_18%d1Array(i_lnd_m1) = d_dmltPturb_dlndm1
            d_mlt_Pturb_18%d1Array(i_lnd_00) = d_dmltPturb_dlnd00
         end subroutine setup_d_mlt_Pturb         
                  
         subroutine setup_RTI_terms(ierr)
            use auto_diff_support
            integer, intent(out) :: ierr
            type(auto_diff_real_18var_order1) :: v_p1, v_00, v_m1, dvdt_diffusion, &
               f, rho_00, rho_m1, dvdt_kick
            real(dp) :: sigm1, sig00
            RTI_terms_18 = 0d0
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
               f = s% dlnddt_RTI_diffusion_factor*s% eta_RTI(k)/dm_div_A_18
               rho_00 = wrap_d_00(s, k)
               rho_m1 = wrap_d_m1(s, k)
               dvdt_kick = f*(rho_00 - rho_m1)/s% dt ! change v according to direction of lower density
            else
               dvdt_kick = 0d0
            end if            
            RTI_terms_18 = dvdt_diffusion + dvdt_kick            
         end subroutine setup_RTI_terms
         
         subroutine unpack_res18(res18)
            use star_utils, only: unpack_res18_partials
            type(auto_diff_real_18var_order1) :: res18
            real(dp) :: resid1
            integer :: j
            include 'formats'
            call unpack_res18_partials(s, k, nvar, xscale, i_dv_dt, &
               res18, d_dm1, d_d00, d_dp1)
            if (s% rotation_flag .and. s% w_div_wc_flag .and. s% use_gravity_rotation_correction) then
               call e00(s, xscale, i_dv_dt, s% i_w_div_wc, k, nvar, iXPavg*d_grav_dw*dm_div_A)            
            end if            
            ! do partials wrt composition   
            resid1 = resid1_18%val         
            do j=1,s% species
               d_residual_dxa00(j) = resid1*d_iXPavg_dxa00(j) - iXPavg*d_dXP_dxa00(j)
               call e00(s, xscale, i_dv_dt, j+s% nvar_hydro, k, nvar, d_residual_dxa00(j))
            end do
            if (k > 1) then 
               do j=1,s% species
                  d_residual_dxam1(j) = resid1*d_iXPavg_dxam1(j) - iXPavg*d_dXP_dxam1(j)
                  call em1(s, xscale, i_dv_dt, j+s% nvar_hydro, k, nvar, d_residual_dxam1(j))
               end do
            end if            
         end subroutine unpack_res18
            
      end subroutine get1_momentum_eqn
      
      
      ! returns -G*m/r^2 with possible modifications for rotation
      subroutine expected_HSE_grav_term(s, k, &
            grav, d_grav_dlnR, d_grav_dw, area, d_area_dlnR, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: grav, d_grav_dlnR, d_grav_dw, area, d_area_dlnR
         integer, intent(out) :: ierr
         
         real(dp) :: inv_R2, d_inv_R2_dlnR, m
         logical :: test_partials

         include 'formats'
      
         ! using_Fraley_time_centering
         ! use_gravity_rotation_correction

         ierr = 0

         m = s% m_grav(k)
         
         if (s% using_Fraley_time_centering) then
            area = 4d0*pi*(s% r(k)**2 + s% r(k)*s% r_start(k) + s% r_start(k)**2)/3d0
            d_area_dlnR = 4d0*pi*s% r(k)*(2d0*s% r(k) + s% r_start(k))/3d0
            inv_R2 = 1d0/(s% r(k)*s% r_start(k))
            d_inv_R2_dlnR = -1d0*inv_R2
         else
            area = 4d0*pi*s% r(k)**2
            d_area_dlnR = 2d0*area
            inv_R2 = 1d0/s% r(k)**2
            d_inv_R2_dlnR = -2d0*inv_R2
         end if

         grav = -s% cgrav(k)*m*inv_R2
         d_grav_dlnR = -s% cgrav(k)*m*d_inv_R2_dlnR

         if (s% rotation_flag .and. s% use_gravity_rotation_correction) then
            grav = grav*s% fp_rot(k) 
            d_grav_dlnR = d_grav_dlnR*s% fp_rot(k)
            if (s% w_div_wc_flag) then
               d_grav_dw = grav/s% fp_rot(k)*s% dfp_rot_dw_div_wc(k)
            else
               d_grav_dw = 0d0
            end if
         else
            d_grav_dw = 0d0
         end if

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         if (test_partials) then
            s% solver_test_partials_val = grav
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = d_grav_dlnR
            write(*,*) 'expected_HSE_grav_term', s% solver_test_partials_var
         end if
      
      end subroutine expected_HSE_grav_term
      
      
      ! other = s% extra_grav(k) - s% dv_dt(k)
      subroutine expected_non_HSE_term(s, k, other_18, other, ierr)
         use hydro_eturb, only: calc_Uq_18
         use accurate_sum_auto_diff_18var_order1
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(out) :: other_18
         real(dp), intent(out) :: other
         integer, intent(out) :: ierr
         type(auto_diff_real_18var_order1) :: &
            extra_18, accel_18, v_00, Uq_18
         real(dp) :: accel, d_accel_dv, fraction_on
         logical :: test_partials, local_v_flag

         include 'formats'

         ierr = 0
         
         extra_18 = 0d0
         if (s% use_other_momentum .or. s% use_other_momentum_implicit) then
            if (s% use_other_momentum_implicit) then
               call wrap(extra_18, s% extra_grav(k), &
                  s% d_extra_grav_dlndm1(k), s% d_extra_grav_dlnd00(k), 0d0, &
                  s% d_extra_grav_dlnTm1(k), s% d_extra_grav_dlnT00(k), 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, s% d_extra_grav_dlnR(k), 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, s% d_extra_grav_dL(k), 0d0)
            else
               extra_18%val = s% extra_grav(k)
            end if
         end if
         
         accel_18 = 0d0
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
            accel_18%val = accel
            accel_18%d1Array(i_v_00) = d_accel_dv
         
         end if ! v_flag

         Uq_18 = 0d0
         if (s% et_flag) then ! Uq(k) is turbulent viscosity drag at face k
            call calc_Uq_18(s, k, Uq_18, ierr)
            if (ierr /= 0) return
         end if
         
         other_18 = extra_18 - accel_18 + Uq_18
         other = other_18%val
         
         if (.false.) then
            if (is_bad(extra_18%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 extra_18', k, extra_18%d1Array(i_lnd_m1)
               stop 'expected_non_HSE_term'
            end if
            if (is_bad(accel_18%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 accel_18', k, accel_18%d1Array(i_lnd_m1)
               stop 'expected_non_HSE_term'
            end if
            if (is_bad(Uq_18%d1Array(i_lnd_m1))) then
               write(*,2) 'lnd_m1 Uq_18', k, Uq_18%d1Array(i_lnd_m1)
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
      subroutine get_dXP_face_info(s, k, P_surf_18, &
            dXP_18, dXP, d_dXP_dxam1, d_dXP_dxa00, &
            iXPavg_18, iXPavg, d_iXPavg_dxam1, d_iXPavg_dxa00, ierr)
         use star_utils, only: calc_XP_18_tw
         use accurate_sum_auto_diff_18var_order1
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(in) :: P_surf_18 ! only used if k == 1
         type(auto_diff_real_18var_order1), intent(out) :: dXP_18, iXPavg_18
         real(dp), intent(out) :: dXP, iXPavg
         real(dp), intent(out), dimension(s% species) :: &
            d_dXP_dxam1, d_dXP_dxa00, d_iXPavg_dxam1, d_iXPavg_dxa00
         integer, intent(out) :: ierr
         
         real(dp) :: XPm1, XP00, XPavg, alfa, beta
         real(dp), dimension(s% species) :: &
            d_XPm1_dxam1, d_XP00_dxa00, d_XPavg_dxam1, d_XPavg_dxa00
         type(auto_diff_real_18var_order1) :: &
            XP00_18, XPm1_18, XPavg_18
         integer :: j
         logical, parameter :: skip_P = .false., skip_mlt_Pturb = .true.
         logical :: test_partials

         include 'formats'

         ierr = 0
         
         call calc_XP_18_tw(s, k, skip_P, skip_mlt_Pturb, XP00_18, d_XP00_dxa00, ierr)
         if (ierr /= 0) return
         XP00 = XP00_18%val
            
         if (k > 1) then
            call calc_XP_18_tw(s, k-1, skip_P, skip_mlt_Pturb, XPm1_18, d_XPm1_dxam1, ierr)
            if (ierr /= 0) return
            XPm1_18 = shift_m1(XPm1_18)
         else ! k == 1
            XPm1_18 = P_surf_18
         end if
         XPm1 = XPm1_18%val
            
         dXP_18 = XPm1_18 - XP00_18
         dXP = XPm1 - XP00
         do j=1,s% species
            d_dXP_dxam1(j) = d_XPm1_dxam1(j)
            d_dXP_dxa00(j) = -d_XP00_dxa00(j)
         end do

         if (k == 1) then
            XPavg_18 = XP00_18
            do j=1,s% species
               d_XPavg_dxam1(j) = 0d0  
               d_XPavg_dxa00(j) = d_XP00_dxa00(j)
            end do
         else
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1d0 - alfa
            XPavg_18 = alfa*XP00_18 + beta*XPm1_18
            do j=1,s% species
               d_XPavg_dxam1(j) = beta*d_XPm1_dxam1(j)
               d_XPavg_dxa00(j) = alfa*d_XP00_dxa00(j)
            end do
         end if
         XPavg = XPavg_18%val
         
         iXPavg_18 = 1d0/XPavg_18
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


      end module hydro_momentum

