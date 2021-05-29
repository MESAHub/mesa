! ***********************************************************************
!
!   Copyright (C) 2010-2020  The MESA Team
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

      module hydro_rsp2_support

      use star_private_def
      use const_def
      use utils_lib, only: is_bad
      use auto_diff
      use auto_diff_support
      use accurate_sum_auto_diff_star_order1
      use star_utils

      implicit none

      private
      public :: remesh_for_RSP2

      contains

      
      subroutine remesh_for_RSP2(s,ierr)
         ! uses these controls
         !  RSP2_nz = 150
         !  RSP2_nz_outer = 40
         !  RSP2_T_anchor = 11d3
         !  RSP2_dq_1_factor = 2d0
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interpolate_vector_pm
         type (star_info), pointer :: s      
         integer, intent(out) :: ierr    
         integer :: k, j, nz_old, nz
         real(dp) :: xm_anchor, P_surf, T_surf, &
            old_logT_cntr, old_logT_1, old_L1, old_r1
         real(dp), allocatable, dimension(:) :: &
            xm_old, xm, xm_mid_old, xm_mid, v_old, v_new
         real(dp), pointer :: work1(:) ! =(nz_old+1, pm_work_size)
         include 'formats'         
         ierr = 0
         nz_old = s% nz
         nz = s% RSP2_nz
         if (nz == nz_old) return ! assume have already done remesh for RSP2
         if (nz > nz_old) stop 'remesh_for_RSP2 cannot increase nz'
         call setvars(ierr)
         if (ierr /= 0) stop 'remesh_for_RSP2 failed in setvars'
         old_L1 = s% L(1)
         old_r1 = s% r(1)
         call set_phot_info(s) ! sets Teff
         call get_PT_surf(P_surf, T_surf, ierr)
         if (ierr /= 0) stop 'remesh_for_RSP2 failed in get_PT_surf'
         allocate(&
            xm_old(nz_old+1), xm_mid_old(nz_old), v_old(nz_old+1), &
            xm(nz+1), xm_mid(nz), v_new(nz+1), work1((nz_old+1)*pm_work_size))
         call set_xm_old
         call find_xm_anchor
         call set_xm_new
         call interpolate1_face_val(s% i_lnR, log(max(1d0,s% r_center)))
         call check_new_lnR
         call interpolate1_face_val(s% i_lum, s% L_center)
         if (s% i_v /= 0) call interpolate1_face_val(s% i_v, s% v_center)
         call set_new_lnd
         call interpolate1_cell_val(s% i_lnT)
         call interpolate1_cell_val(s% i_w)
         do j=1,s% species
            call interpolate1_xa(j)
         end do
         call rescale_xa
         call revise_lnT_for_QHSE(P_surf, ierr)
         if (ierr /= 0) stop 'remesh_for_RSP2 failed in revise_lnT_for_QHSE'
         do k=1,nz
            call set_Hp_face(k)
         end do
         deallocate(work1)  
         s% nz = nz
         write(*,1) 'new old L_surf/Lsun', s% xh(s% i_lum,1)/Lsun, old_L1/Lsun
         write(*,1) 'new old R_surf/Rsun', exp(s% xh(s% i_lnR,1))/Rsun, old_r1/Rsun
         write(*,*)
         !stop 'remesh_for_RSP2'
         
         contains

         subroutine setvars(ierr)
            use hydro_vars, only: unpack_xh, set_hydro_vars
            integer, intent(out) :: ierr
            logical, parameter :: &
               skip_basic_vars = .false., &
               skip_micro_vars = .false., &
               skip_m_grav_and_grav = .false., &
               skip_net = .true., &
               skip_neu = .true., &
               skip_kap = .false., &
               skip_grads = .true., &
               skip_rotation = .true., &
               skip_brunt = .true., &
               skip_other_cgrav = .true., &
               skip_mixing_info = .true., &
               skip_set_cz_bdy_mass = .true., &
               skip_mlt = .true., &
               skip_eos = .false.
            ierr = 0
            call unpack_xh(s,ierr)
            if (ierr /= 0) stop 'remesh_for_RSP2 failed in unpack_xh'
            call set_hydro_vars( &
               s, 1, nz_old, skip_basic_vars, &
               skip_micro_vars, skip_m_grav_and_grav, skip_eos, skip_net, skip_neu, &
               skip_kap, skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
               skip_mixing_info, skip_set_cz_bdy_mass, skip_mlt, ierr)
            if (ierr /= 0) stop 'remesh_for_RSP2 failed in set_hydro_vars'
         end subroutine setvars
         
         subroutine get_PT_surf(P_surf, T_surf, ierr)
            use atm_support, only: get_atm_PT
            real(dp), intent(out) :: P_surf, T_surf
            integer, intent(out) :: ierr
            real(dp) :: &
               Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
            logical, parameter :: skip_partials = .true.
            include 'formats'
            ierr = 0
            call set_phot_info(s) ! sets s% Teff
            Teff = s% Teff
            call get_atm_PT( & ! this uses s% opacity(1)
                 s, s% tau_factor*s% tau_base, s% L(1), s% r(1), s% m(1), s% cgrav(1), skip_partials, &
                 Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
                 lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
            if (ierr /= 0) stop 'get_P_surf failed in get_atm_PT'
            P_surf = exp(lnP_surf)
            T_surf = exp(lnT_surf)
            return
            
            write(*,1) 'get_PT_surf P_surf', P_surf
            write(*,1) 'get_PT_surf T_surf', T_surf
            write(*,1) 'get_PT_surf Teff', Teff
            write(*,1) 'get_PT_surf opacity(1)', s% opacity(1)
            write(*,1)
            !stop 'get_PT_surf'
         end subroutine get_PT_surf
         
         subroutine set_xm_old
            xm_old(1) = 0d0
            do k=2,nz_old
               xm_old(k) = xm_old(k-1) + s% dm(k-1)
            end do
            xm_old(nz_old+1) = s% xmstar
            do k=1,nz_old
               xm_mid_old(k) = xm_old(k) + 0.5d0*s% dm(k)
            end do
         end subroutine set_xm_old
         
         subroutine find_xm_anchor
            real(dp) :: lnT_anchor, xmm1, xm00, lnTm1, lnT00
            include 'formats'
            lnT_anchor = log(s% RSP2_T_anchor)
            if (lnT_anchor <= s% xh(s% i_lnT,1)) then
               write(*,1) 'T_anchor < T_surf', s% RSP2_T_anchor, exp(s% xh(s% i_lnT,1))
               stop 'find_xm_anchor'
            end if
            xm_anchor = xm_old(nz_old)
            do k=2,nz_old
               if (s% xh(s% i_lnT,k) >= lnT_anchor) then
                  xmm1 = xm_old(k-1)
                  xm00 = xm_old(k)
                  lnTm1 = s% xh(s% i_lnT,k-1)
                  lnT00 = s% xh(s% i_lnT,k)
                  xm_anchor = xmm1 + &
                     (xm00 - xmm1)*(lnT_anchor - lnTm1)/(lnT00 - lnTm1)
                  if (is_bad(xm_anchor) .or. xm_anchor <= 0d0) then
                     write(*,2) 'bad xm_anchor', k, xm_anchor, xmm1, xm00, lnTm1, lnT00, lnT_anchor, s% lnT(1)
                     stop 'find_xm_anchor'
                  end if
                  return
               end if
            end do
         end subroutine find_xm_anchor         
         
         subroutine set_xm_new ! sets xm, dm, m, dq, q
            integer :: nz_outer, n_inner, iter, k, j
            real(dp) :: dq_1_factor, dxm_outer, lnx, dlnx
            include 'formats'
            nz_outer = s% RSP2_nz_outer
            dq_1_factor = s% RSP2_dq_1_factor
            dxm_outer = xm_anchor/(nz_outer - 1d0 + dq_1_factor)
            !write(*,2) 'dxm_outer', nz_outer, dxm_outer, xm_anchor
            xm(1) = 0d0
            xm(2) = dxm_outer*dq_1_factor
            s% dm(1) = xm(2)
            do k=3,nz_outer+1
               xm(k) = xm(k-1) + dxm_outer
               s% dm(k-1) = dxm_outer
            end do
            lnx = log(xm(nz_outer+1))
            if (is_bad(lnx)) then
               write(*,2) 'bad lnx', nz_outer+1, lnx, xm(nz_outer+1)
               stop 'set_xm_new'
            end if
            dlnx = (log(s% xmstar) - lnx)/(nz - nz_outer)
            do k=nz_outer+2,nz
               lnx = lnx + dlnx
               xm(k) = exp(lnx)
               s% dm(k-1) = xm(k) - xm(k-1)
            end do
            s% dm(nz) = s% xmstar - xm(nz)
            do k=1,nz-1
               xm_mid(k) = 0.5d0*(xm(k) + xm(k+1))
            end do
            xm_mid(nz) = 0.5d0*(xm(nz) + s% xmstar)
            s% m(1) = s% mstar
            s% q(1) = 1d0
            s% dq(1) = s% dm(1)/s% xmstar
            do k=2,nz
               s% m(k) = s% m(k-1) - s% dm(k-1)
               s% dq(k) = s% dm(k)/s% xmstar
               s% q(k) = s% q(k-1) - s% dq(k-1)
            end do
            call set_dm_bar(s, s% nz, s% dm, s% dm_bar)
            return
            
            do k=2,nz
               write(*,2) 'dm(k)/dm(k-1) m(k)', k, s%dm(k)/s%dm(k-1), s%m(k)/Msun
            end do
            write(*,1) 'm_center', s% m_center/msun
            stop 'set_xm_new'
         end subroutine set_xm_new
         
         subroutine interpolate1_face_val(i, cntr_val)
            integer, intent(in) :: i
            real(dp), intent(in) :: cntr_val
            do k=1,nz_old
               v_old(k) = s% xh(i,k)
            end do
            v_old(nz_old+1) = cntr_val
            call interpolate_vector_pm( &
               nz_old+1, xm_old, nz+1, xm, v_old, v_new, work1, 'remesh_for_RSP2', ierr)
            do k=1,nz
               s% xh(i,k) = v_new(k)
            end do
         end subroutine interpolate1_face_val
         
         subroutine check_new_lnR
            include 'formats'
            do k=1,nz
               s% lnR(k) = s% xh(s% i_lnR,k)
               s% r(k) = exp(s% lnR(k))
            end do
            do k=1,nz-1
               if (s% r(k) <= s% r(k+1)) then
                  write(*,2) 'bad r', k, s% r(k), s% r(k+1)
                  stop 'check_new_lnR remesh rsp2'
               end if
            end do
            if (s% r(nz) <= s% r_center) then
               write(*,2) 'bad r center', nz, s% r(nz), s% r_center
               stop 'check_new_lnR remesh rsp2'
            end if
         end subroutine check_new_lnR
         
         subroutine set_new_lnd
            real(dp) :: vol, r300, r3p1
            include 'formats'
            do k=1,nz
               r300 = pow3(s% r(k))
               if (k < nz) then
                  r3p1 = pow3(s% r(k+1))
               else
                  r3p1 = pow3(s% r_center)
               end if
               vol = (4d0*pi/3d0)*(r300 - r3p1)
               s% rho(k) = s% dm(k)/vol
               s% lnd(k) = log(s% rho(k))
               s% xh(s% i_lnd,k) = s% lnd(k)
               if (is_bad(s% lnd(k))) then
                  write(*,2) 'bad lnd vol dm r300 r3p1', k, s% lnd(k), vol, s% dm(k), r300, r3p1
                  stop 'remesh for rsp2'
               end if
            end do
         end subroutine set_new_lnd
         
         subroutine interpolate1_cell_val(i)
            integer, intent(in) :: i
            do k=1,nz_old
               v_old(k) = s% xh(i,k)
            end do
            call interpolate_vector_pm( &
               nz_old, xm_mid_old, nz, xm_mid, v_old, v_new, work1, 'remesh_for_RSP2', ierr)
            do k=1,nz
               s% xh(i,k) = v_new(k)
            end do
         end subroutine interpolate1_cell_val
         
         subroutine interpolate1_xa(j)
            integer, intent(in) :: j
            do k=1,nz_old
               v_old(k) = s% xa(j,k)
            end do
            call interpolate_vector_pm( &
               nz_old, xm_mid_old, nz, xm_mid, v_old, v_new, work1, 'remesh_for_RSP2', ierr)
            do k=1,nz
               s% xa(j,k) = v_new(k)
            end do
         end subroutine interpolate1_xa
         
         subroutine rescale_xa
            integer :: k, j
            real(dp) :: sum_xa
            do k=1,nz
               sum_xa = sum(s% xa(1:s% species,k))
               do j=1,s% species
                  s% xa(j,k) = s% xa(j,k)/sum_xa
               end do
            end do
         end subroutine rescale_xa
         
         subroutine revise_lnT_for_QHSE(P_surf, ierr)
            use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results
            use chem_lib
            use chem_def
            use eos_support, only: solve_eos_given_DP
            use eos_def
            use kap_def, only: num_kap_fracs
            use kap_support, only: get_kap
            real(dp), intent(in) :: P_surf
            integer, intent(out) :: ierr
            real(dp) :: logRho, logP, logT_guess, &
               logT_tol, logP_tol, logT, P_m1, P_00, dm_face, &
               kap_fracs(num_kap_fracs), kap, dlnkap_dlnRho, dlnkap_dlnT, &
               old_kap, new_P_surf, new_T_surf
            real(dp), dimension(num_eos_basic_results) :: &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar
            real(dp) :: dres_dxa(num_eos_d_dxa_results,s% species)
            include 'formats'
            ierr = 0
            P_m1 = P_surf
            do k=1,nz
               s% lnT(k) = s% xh(s% i_lnT,k)
               s% lnR(k) = s% xh(s% i_lnR,k)
               s% r(k) = exp(s% lnR(k))
            end do
            !write(*,1) 'before revise_lnT_for_QHSE: logT cntr', s% lnT(nz)/ln10
            do k=1,nz
               if (k < nz) then
                  dm_face = s% dm_bar(k)
               else
                  dm_face = 0.5d0*(s% dm(k-1) + s% dm(k))
               end if
               P_00 = P_m1 + s% cgrav(k)*s% m(k)*dm_face/(4d0*pi*pow4(s% r(k)))
               logP = log10(P_00) ! value for QHSE
               s% lnPeos(k) = logP/ln10
               s% Peos(k) = P_00
               logRho = s% lnd(k)/ln10
               logT_guess = s% lnT(k)/ln10
               logT_tol = 1d-11
               logP_tol = 1d-11
               call solve_eos_given_DP( &
                  s, k, s% xa(:,k), &
                  logRho, logP, logT_guess, logT_tol, logP_tol, &
                  logT, res, d_dlnd, d_dlnT, dres_dxa, ierr)
               if (ierr /= 0) then
                  write(*,2) 'solve_eos_given_DP failed', k
                  write(*,*)
                  write(*,1) 'sum(xa)', sum(s% xa(:,k))
                  do j=1,s% species
                     write(*,4) 'xa(j,k) ' // trim(chem_isos% name(s% chem_id(j))), j, j+s% nvar_hydro, k, s% xa(j,k)
                  end do
                  write(*,1) 'logRho', logRho
                  write(*,1) 'logP', logP
                  write(*,1) 'logT_guess', logT_guess
                  write(*,1) 'logT_tol', logT_tol
                  write(*,1) 'logP_tol', logP_tol
                  write(*,*)
                  stop 'revise_lnT_for_QHSE'
               end if
               s% lnT(k) = logT*ln10
               s% xh(s% i_lnT,k) = s% lnT(k)
               !write(*,2) 'logP dlogT logT logT_guess logRho', k, &
               !   logP, logT - logT_guess, logT, logT_guess, logRho
               P_m1 = P_00
               
               if (k == 1) then ! get opacity and recheck surf BCs
                 call get_kap( & ! assume zbar is set
                     s, k, s% zbar(k), s% xa(:,k), logRho, logT, &
                     res(i_lnfree_e), d_dlnd(i_lnfree_e), d_dlnT(i_lnfree_e), &
                     res(i_eta), d_dlnd(i_eta), d_dlnT(i_eta), &
                     kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, &
                     ierr)
                  if (ierr /= 0) then
                     write(*,2) 'get_kap failed', k
                     stop 'revise_lnT_for_QHSE'
                  end if
                  old_kap = s% opacity(1)
                  s% opacity(1) = kap ! for use by atm surf PT
                  call get_PT_surf(new_P_surf, new_T_surf, ierr)
                  if (ierr /= 0) then
                     write(*,2) 'get_PT_surf failed', k
                     stop 'revise_lnT_for_QHSE'
                  end if
                  write(*,1) 'new old T_surf', new_T_surf, T_surf
                  write(*,1) 'new old P_surf', new_P_surf, P_surf
                  write(*,1) 'new old kap(1)', kap, old_kap
                  !stop 'revise_lnT_for_QHSE'
               end if
               
            end do
            !write(*,1) 'after revise_lnT_for_QHSE: logT cntr', s% lnT(nz)/ln10
            !stop
         end subroutine revise_lnT_for_QHSE

         subroutine set_Hp_face(k)
            use hydro_rsp2, only: get_RSP2_alfa_beta_face_weights
            integer, intent(in) :: k
            real(dp) :: r_00, d_00, Peos_00, Peos_div_rho, Hp_face, &
               d_m1, Peos_m1, alfa, beta
            r_00 = s% r(k)
            d_00 = s% rho(k)
            Peos_00 = s% Peos(k)
            if (k == 1) then
               Peos_div_rho = Peos_00/d_00
               Hp_face = pow2(r_00)*Peos_div_rho/(s% cgrav(k)*s% m(k))
            else
               d_m1 = s% rho(k-1)
               Peos_m1 = s% Peos(k-1)
               call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
               Peos_div_rho = alfa*Peos_00/d_00 + beta*Peos_m1/d_m1
               Hp_face = pow2(r_00)*Peos_div_rho/(s% cgrav(k)*s% m(k))
            end if
            s% Hp_face(k) = Hp_face
            s% xh(s% i_Hp, k) = Hp_face
         end subroutine set_Hp_face

      end subroutine remesh_for_RSP2


      end module hydro_rsp2_support

