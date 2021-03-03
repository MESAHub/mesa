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
      public :: do1_momentum_eqn, do_surf_momentum_eqn, do1_radius_eqn, &
         expected_HSE_grav_term


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
         call store_partials( &
            s, 1, s% i_dv_dt, nvar, d_dm1, d_d00, d_dp1, 'do_surf_momentum_eqn', ierr)
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
         call store_partials( &
            s, k, s% i_dv_dt, nvar, d_dm1, d_d00, d_dp1, 'do1_momentum_eqn', ierr)
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
         
         real(dp) :: residual, dm_face, dXP, iXPavg, dm_div_A
         real(dp), dimension(s% species) :: &
            d_dXP_dxam1, d_dXP_dxa00, d_iXPavg_dxam1, d_iXPavg_dxa00, &
            d_residual_dxam1, d_residual_dxa00
         integer :: nz, j, i_dv_dt, i_lum, i_v
         logical :: test_partials
         
         type(auto_diff_real_star_order1) :: resid1_ad, resid_ad, &
            other_ad, dm_div_A_ad, grav_ad, area_ad, dXP_ad, d_mlt_Pturb_ad, &
            iXPavg_ad, other_dm_div_A_ad, grav_dm_div_A_ad, &
            RTI_terms_ad, RTI_terms_dm_div_A_ad, accel_ad, Uq_ad
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

         call setup_HSE(dm_div_A, ierr); if (ierr /= 0) return ! grav_ad and dm_div_A_ad
         call setup_non_HSE(ierr); if (ierr /= 0) return ! other = s% extra_grav(k) - dv_dt
         call setup_dXP(ierr); if (ierr /= 0) return ! dXP_ad, iXPavg_ad
         call setup_d_mlt_Pturb(ierr); if (ierr /= 0) return ! d_mlt_Pturb_ad
         call setup_RTI_terms(ierr); if (ierr /= 0) return ! RTI_terms_ad
         
         other_dm_div_A_ad = other_ad*dm_div_A_ad
         grav_dm_div_A_ad = grav_ad*dm_div_A_ad
         RTI_terms_dm_div_A_ad = RTI_terms_ad*dm_div_A_ad
         
         ! sum terms in residual_sum_ad using accurate_auto_diff_real_star_order1
         residual_sum_ad = &
            other_dm_div_A_ad + grav_dm_div_A_ad - dXP_ad - d_mlt_Pturb_ad + RTI_terms_dm_div_A_ad
         
         resid1_ad = residual_sum_ad ! convert back to auto_diff_real_star_order1
         resid_ad = resid1_ad*iXPavg_ad ! scaling
         residual = resid_ad%val
         s% equ(i_dv_dt, k) = residual      
         
         !s% xtra1_array(k) = accel_ad%val
         !s% xtra2_array(k) = safe_log10(-dXP_ad%val/dm_div_A_ad%val)
         !s% xtra3_array(k) = safe_log10(-grav_ad%val)
         !s% xtra4_array(k) = Uq_ad%val
         !s% xtra5_array(k) = s% Chi(k)

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
         call unpack_res18(s% species, resid_ad)

         if (test_partials) then
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = d_d00(s% solver_test_partials_var)
            write(*,*) 'get1_momentum_eqn', s% solver_test_partials_var
         end if
         
         contains
         
         subroutine init
            i_dv_dt = s% i_dv_dt
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
         
         subroutine setup_HSE(dm_div_A, ierr)
            real(dp), intent(out) :: dm_div_A
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            call expected_HSE_grav_term(s, k, grav_ad, area_ad, ierr)
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
            call expected_non_HSE_term(s, k, other_ad, other, accel_ad, Uq_ad, ierr)
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
               if (.not. s% solver_use_lnd) then
                  d_dmltPturb_dlndm1 = d_dmltPturb_dlndm1/s% rho(k-1)
                  d_dmltPturb_dlnd00 = d_dmltPturb_dlnd00/s% rho(k)
               end if
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
         
         subroutine unpack_res18(species, res18)
            use star_utils, only: save_eqn_dxa_partials, unpack_residual_partials
            integer, intent(in) :: species
            type(auto_diff_real_star_order1) :: res18
            real(dp) :: resid1, dxap1(species)
            logical, parameter :: checking = .true.
            integer :: j
            include 'formats'
            ! do partials wrt composition   
            resid1 = resid1_ad%val         
            do j=1,species
               d_residual_dxa00(j) = resid1*d_iXPavg_dxa00(j) - iXPavg*d_dXP_dxa00(j)
               if (checking) call check_dequ(d_dXP_dxa00(j),'d_dXP_dxa00(j)')
               if (checking) call check_dequ(d_iXPavg_dxa00(j),'d_iXPavg_dxa00(j)')
            end do
            if (k > 1) then 
               do j=1,species
                  d_residual_dxam1(j) = resid1*d_iXPavg_dxam1(j) - iXPavg*d_dXP_dxam1(j)
                  if (checking) call check_dequ(d_dXP_dxam1(j),'d_dXP_dxam1(j)')
                  if (checking) call check_dequ(d_iXPavg_dxam1(j),'d_iXPavg_dxam1(j)')
               end do
            else
               d_residual_dxam1 = 0d0
            end if            
            dxap1 = 0d0
            call save_eqn_dxa_partials(&
               s, k, nvar, i_dv_dt, species, &
               d_residual_dxam1, d_residual_dxa00, dxap1, 'get1_momentum_eqn', ierr)
            call unpack_residual_partials(s, k, nvar, i_dv_dt, &
               res18, d_dm1, d_d00, d_dp1)
         end subroutine unpack_res18

         subroutine check_dequ(dequ, str)
            real(dp), intent(in) :: dequ
            character (len=*), intent(in) :: str
            include 'formats'
            if (is_bad(dequ)) then
!$omp critical (hydro_momentum_crit2)
               ierr = -1
               if (s% report_ierr) then
                  write(*,2) 'get1_momentum_eqn: bad ' // trim(str), k, dequ
               end if
               if (s% stop_for_bad_nums) stop 'get1_momentum_eqn'
!$omp end critical (hydro_momentum_crit2)
               return
            end if
         end subroutine check_dequ
            
      end subroutine get1_momentum_eqn
      
      
      ! returns -G*m/r^2 with possible modifications for rotation.  MESA 2, eqn 22.
      subroutine expected_HSE_grav_term(s, k, grav, area, ierr)
         use star_utils, only: get_area_info
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: area, grav
         integer, intent(out) :: ierr
         
         type(auto_diff_real_star_order1) :: inv_R2
         logical :: test_partials

         include 'formats'
         ierr = 0
         
         call get_area_info(s, k, area, inv_R2, ierr)
         if (ierr /= 0) return

         grav = -s% cgrav(k)*s% m_grav(k)*inv_R2
         
         if (s% rotation_flag .and. s% use_gravity_rotation_correction) &
            grav = grav*s% fp_rot(k) 

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
      subroutine expected_non_HSE_term( &
            s, k, other_ad, other, accel_ad, Uq_ad, ierr)
         use hydro_tdc, only: compute_Uq_face
         use accurate_sum_auto_diff_star_order1
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: &
            other_ad, accel_ad,Uq_ad
         real(dp), intent(out) :: other
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: extra_ad, v_00
         real(dp) :: accel, d_accel_dv, fraction_on, dlnR00, &
            dlnTm1, dlnT00, dlndm1, dlnd00
         logical :: test_partials, local_v_flag

         include 'formats'

         ierr = 0
         
         extra_ad = 0d0
         if (s% use_other_momentum .or. s% use_other_momentum_implicit) then
            if (s% use_other_momentum_implicit) then
               dlnR00 = s% d_extra_grav_dlnR(k)
               if (.not. s% solver_use_lnR) dlnR00 = dlnR00/s% r(k)
               dlnTm1 = s% d_extra_grav_dlnTm1(k)
               dlnT00 = s% d_extra_grav_dlnT00(k)
               if (.not. s% solver_use_lnT) then
                  dlnTm1 = dlnTm1/s% T(k-1)
                  dlnT00 = dlnT00/s% T(k)
               end if
               dlndm1 = s% d_extra_grav_dlndm1(k)
               dlnd00 = s% d_extra_grav_dlnd00(k)
               if (.not. s% solver_use_lnd) then
                  dlndm1 = dlndm1/s% rho(k-1)
                  dlnd00 = dlnd00/s% rho(k)
               end if
               call wrap(extra_ad, s% extra_grav(k), &
                  dlndm1, dlnd00, 0d0, &
                  dlnTm1, dlnT00, 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, dlnR00, 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, s% d_extra_grav_dL(k), 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, 0d0, 0d0)
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
               accel = s% dxh_v(k)/s% dt
               d_accel_dv = 1d0/s% dt
            else ! assume vstart(k) = 0 and
               ! constant acceleration dv_dt so vfinal(k) = dv_dt*dt
               ! v(k) = dr/dt = average velocity =
               !      (vstart + vfinal)/2 = dv_dt*dt/2 when vstart = 0
               ! so (1/2)*dv_dt*dt = v(k)
               accel = 2d0*s% v(k)/s% dt
               d_accel_dv = 2d0/s% dt
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
         
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         if (test_partials) then
            s% solver_test_partials_val = 0
            s% solver_test_partials_var = 0
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


      subroutine do1_radius_eqn(s, k, skip_partials, nvar, ierr)
         use auto_diff_support
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            v00, r_actual, r_expected, dxh_lnR, resid_ad, &
            dr_div_r0_actual, dr_div_r0_expected, dr
         logical :: test_partials, force_zero_v
         include 'formats'
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         ierr = 0         
         if (.not. (s% u_flag .or. s% v_flag)) stop 'must have either v or u for do1_radius_eqn'
         
         force_zero_v = (s% q(k) > s% velocity_q_upper_bound) .or. &
            (s% lnT_start(k)/ln10 < s% velocity_logT_lower_bound .and. &
               s% dt < secyer*s% max_dt_yrs_for_velocity_logT_lower_bound)                  
         if (force_zero_v) then
            if (s% u_flag) then
               v00 = wrap_u_00(s,k)
            else
               v00 = wrap_v_00(s,k)
            end if
            resid_ad = v00/s% csound_start(k)
            if (skip_partials) return            
            call save_eqn_residual_info( &
               s, k, nvar, s% i_dlnR_dt, resid_ad, 'do1_radius_eqn', ierr)           
            return
         end if
         
         if (s% u_flag) then
            v00 = s% u_face_ad(k)
            if (s% using_velocity_time_centering) &
               v00 = 0.5d0*(v00 + s% u_face_start(k))
         else
            v00 = wrap_opt_time_center_v_00(s,k)
         end if         
         
         if (s% solver_use_lnR) then
            ! dr = r - r0 = v00*dt
            ! eqn: dr/r0 = v00*dt/r0
            ! (r - r0)/r0 = r/r0 - 1 = exp(lnR)/exp(lnR0) - 1
            ! = exp(lnR - lnR0) - 1 = exp(dlnR) - 1 = exp(dlnR_dt*dt) - 1
            ! eqn becomes: v00*dt/r0 = expm1(dlnR)
            dxh_lnR = wrap_dxh_lnR(s,k) ! lnR - lnR_start
            dr_div_r0_actual = expm1(dxh_lnR) ! expm1(x) = E^x - 1
            dr_div_r0_expected = v00*s% dt/s% r_start(k)
            resid_ad = dr_div_r0_expected - dr_div_r0_actual
         else ! not solver_use_lnR
            dr = wrap_dxh_lnR(s,k) ! r - r_start
            resid_ad = (dr - v00*s% dt)/s% r_start(k)
         end if
         
         s% equ(s% i_dlnR_dt, k) = resid_ad%val
         
         if (test_partials) then
            s% solver_test_partials_val = 0
         end if
         if (skip_partials) return            
         call save_eqn_residual_info( &
            s, k, nvar, s% i_dlnR_dt, resid_ad, 'do1_radius_eqn', ierr)           
         if (test_partials) then   
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do1_radius_eqn', s% solver_test_partials_var
         end if
      end subroutine do1_radius_eqn


      end module hydro_momentum

