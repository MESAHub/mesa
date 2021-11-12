! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

      module hydro_vars

      use star_private_def
      use const_def
      use chem_def, only: chem_isos
      use utils_lib, only: mesa_error, is_bad

      implicit none

      private
      public :: set_vars_if_needed, set_vars, set_final_vars, update_vars, set_cgrav, &
         set_hydro_vars, unpack_xh, set_Teff_info_for_eqns, set_Teff, get_surf_PT, set_grads

      logical, parameter :: dbg = .false.
      logical, parameter :: trace_setvars = .false.


      contains


      subroutine set_vars_if_needed(s, dt, str, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         if (.not. s% need_to_setvars) then
            if (trace_setvars) write(*,*) '** skip set_vars ' // trim(str)
            s% num_skipped_setvars = s% num_skipped_setvars + 1
            return
         end if
         if (trace_setvars) write(*,*) 'set_vars ' // trim(str)
         call set_vars(s, dt, ierr)
      end subroutine set_vars_if_needed


      subroutine set_vars(s, dt, ierr)
         use star_utils, only: reset_starting_vectors
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         logical, parameter :: &
            skip_m_grav_and_grav = .false., &
            skip_net = .false., &
            skip_neu = .false., &
            skip_kap = .false., &
            skip_grads = .false., &
            skip_rotation = .false., &
            skip_brunt = .false., &
            skip_irradiation_heat = .false.
         logical :: skip_set_cz_bdy_mass, skip_mixing_info
         ierr = 0
         s% num_setvars = s% num_setvars + 1
         if (dbg) write(*,*) 'set_vars'
         call reset_starting_vectors(s)
         skip_set_cz_bdy_mass = .not. s% have_new_generation
         skip_mixing_info = .not. s% okay_to_set_mixing_info
         call set_some_vars( &
            s, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat, &
            dt, ierr)
      end subroutine set_vars

      
      subroutine set_final_vars(s, dt, ierr)
         use rates_def, only: num_rvs
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         logical :: &
            skip_basic_vars, &
            skip_micro_vars, &
            skip_m_grav_and_grav, &
            skip_eos, &
            skip_net, &
            skip_neu, &
            skip_kap, &
            skip_grads, &
            skip_rotation, &
            skip_brunt, &
            skip_other_cgrav, &
            skip_mixing_info, &
            skip_set_cz_bdy_mass, &
            skip_mlt
         integer :: nz, ierr1, k, i
         
         include 'formats'
         
         ierr = 0
         nz = s% nz
         
         skip_grads = .false.
         skip_rotation = .false.
         skip_brunt = .false.
         skip_other_cgrav = .false.
         skip_set_cz_bdy_mass = .false.
         skip_m_grav_and_grav = .false.
         skip_mixing_info = .not. s% recalc_mix_info_after_evolve         

         ! only need to do things that were skipped in set_vars_for_solver
         ! i.e., skip what it already did for the last solver iteration
         skip_basic_vars = .not. s% need_to_setvars            
         skip_micro_vars = .not. s% need_to_setvars
         skip_kap = .not. s% need_to_setvars
         skip_neu = .not. s% need_to_setvars
         skip_net = .not. s% need_to_setvars
         skip_eos = .not. s% need_to_setvars
         skip_mlt = .not. s% need_to_setvars
         
         if (s% need_to_setvars) then
            s% num_setvars = s% num_setvars + 1
            if (trace_setvars) write(*,*) 'set_vars in set_final_vars'
         else
            s% num_skipped_setvars = s% num_skipped_setvars + 1
            if (trace_setvars) write(*,*) '** skip set_vars in set_final_vars'
         end if
         if (trace_setvars) write(*,*)
      
         call set_hydro_vars( &
            s, 1, nz, skip_basic_vars, &
            skip_micro_vars, skip_m_grav_and_grav, skip_eos, skip_net, skip_neu, &
            skip_kap, skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_mlt, ierr)
         if (ierr /= 0) then
            if (s% report_ierr .or. dbg) &
               write(*,*) 'update_vars: set_hydro_vars returned ierr', ierr
            return
         end if

         if (s% op_split_burn) then
            do k = 1, nz
               if (s% T_start(k) >= s% op_split_burn_min_T) &
                  s% eps_nuc(k) = s% burn_avg_epsnuc(k)
            end do
         end if

      end subroutine set_final_vars


      subroutine set_some_vars( &
            s, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat, &
            dt, ierr)
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s
         logical, intent(in) :: &
            skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, &
            skip_rotation, skip_brunt, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         integer(8) :: time0, clock_rate
         real(dp) :: total
         logical, parameter :: skip_other_cgrav = .false.
         logical, parameter :: skip_basic_vars = .false.
         logical, parameter :: skip_micro_vars = .false.
         logical, parameter :: skip_mlt = .false.
         logical, parameter :: skip_eos = .false.

         include 'formats'
            
         call update_vars(s, &
            skip_basic_vars, skip_micro_vars, &
            skip_m_grav_and_grav, skip_net, skip_neu, skip_kap, &
            skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat, &
            skip_mlt, skip_eos, dt, ierr)
         if (ierr /= 0) then
            if (s% report_ierr .or. dbg) &
               write(*,*) 'set_some_vars: update_vars returned ierr', ierr
            return
         end if
         
      end subroutine set_some_vars


      subroutine update_vars(s, &
            skip_basic_vars, skip_micro_vars, &
            skip_m_grav_and_grav, skip_net, skip_neu, skip_kap, &
            skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat, &
            skip_mlt, skip_eos, dt, ierr)
         use star_utils, only: eval_irradiation_heat
         type (star_info), pointer :: s
         logical, intent(in) :: &
            skip_basic_vars, skip_micro_vars, &
            skip_m_grav_and_grav, skip_net, skip_neu, skip_kap, skip_grads, &
            skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat, &
            skip_mlt, skip_eos
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         integer :: k, nz
         include 'formats'
         ierr = 0
         nz = s% nz
         if (s% doing_finish_load_model .or. .not. s% RSP_flag) then
            call unpack_xh(s,ierr)
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*,*) 'after unpack ierr', ierr
               return
            end if
            if (.not. skip_mixing_info) then
               s% mixing_type(1:nz) = no_mixing
               s% adjust_mlt_gradT_fraction(1:nz) = -1
            end if            
         end if
         
         call set_hydro_vars( &
            s, 1, nz, skip_basic_vars, &
            skip_micro_vars, skip_m_grav_and_grav, skip_eos, skip_net, skip_neu, &
            skip_kap, skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_mlt, ierr)
         if (ierr /= 0) then
            if (s% report_ierr .or. dbg) &
               write(*,*) 'update_vars: set_hydro_vars returned ierr', ierr
            return
         end if

         if (s% use_other_momentum) then
            call s% other_momentum(s% id, ierr)
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*,*) 'update_vars: other_momentum returned ierr', ierr
               return
            end if
         end if
         
         if (.not. skip_irradiation_heat) then
            if (s% irradiation_flux /= 0) then
               do k=1,nz
                  s% irradiation_heat(k) = eval_irradiation_heat(s,k)
               end do
            else
               s% irradiation_heat(1:nz) = 0
            end if
         end if

      end subroutine update_vars
      
      
      subroutine unpack_xh(s,ierr)
         use star_utils, only: set_qs, set_dm_bar, set_m_and_dm
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: i_lnd, i_lnT, i_lnR, i_w, i_Hp, &
            i_lum, i_v, i_u, i_alpha_RTI, i_Et_RSP, &
            j, k, species, nvar_chem, nz, k_below_just_added
         include 'formats'
         ierr = 0
         nz = s% nz
         k_below_just_added = 1
         species = s% species
         nvar_chem = s% nvar_chem
         i_lnd = s% i_lnd
         i_lnT = s% i_lnT
         i_lnR = s% i_lnR
         i_lum = s% i_lum
         i_w = s% i_w
         i_Hp = s% i_Hp
         i_v = s% i_v
         i_u = s% i_u
         i_alpha_RTI = s% i_alpha_RTI
         i_Et_RSP = s% i_Et_RSP
      
         do j=1,s% nvar_hydro
            if (j == i_lnd) then
               do k=1,nz
                  s% lnd(k) = s% xh(i_lnd,k)
                  s% rho(k) = exp(s% lnd(k))
               end do
               s% dxh_lnd(1:nz) = 0d0
            else if (j == i_lnT) then
               do k=1,nz
                  s% lnT(k) = s% xh(i_lnT,k)
                  s% T(k) = exp(s% lnT(k))
               end do
               s% dxh_lnT(1:nz) = 0d0
            else if (j == i_lnR) then
               do k=1,nz
                  s% lnR(k) = s% xh(i_lnR,k)
                  s% r(k) = exp(s% lnR(k))
               end do
               s% dxh_lnR(1:nz) = 0d0
            else if (j == i_w) then
               do k=1,nz
                  s% w(k) = s% xh(i_w, k)
                  if (s% w(k) < 0d0) then
                     !write(*,4) 'unpack: fix w < 0', k, &
                     !   s% solver_iter, s% model_number, s% w(k)
                     s% w(k) = s% RSP2_w_fix_if_neg
                  end if
               end do
            else if (j == i_Hp) then
               do k=1,nz
                  s% Hp_face(k) = s% xh(i_Hp, k)
               end do
            else if (j == i_lum) then
               do k=1,nz
                  s% L(k) = s% xh(i_lum, k)
               end do
            else if (j == i_v) then
               do k=1,nz
                  s% v(k) = s% xh(i_v,k)
                  s% dxh_v(k) = 0d0
               end do
            else if (j == i_u) then
               do k=1,nz
                  s% u(k) = s% xh(i_u,k)
                  s% dxh_u(k) = 0d0
               end do
            else if (j == i_alpha_RTI) then
               do k=1,nz
                  s% alpha_RTI(k) = max(0d0, s% xh(i_alpha_RTI,k))
                  s% dxh_alpha_RTI(k) = 0d0
               end do
            else if (j == i_Et_RSP) then
               do k=1,nz
                  s% RSP_Et(k) = max(0d0,s% xh(i_Et_RSP,k))
               end do
            end if
         end do

         if (i_lum == 0 .and. .not. s% RSP_flag) s% L(1:nz) = 0d0
         if (i_v == 0) s% v(1:nz) = 0d0
         if (i_u == 0) s% u(1:nz) = 0d0
         if (i_w == 0) s% w(1:nz) = 0d0
         if (i_Hp == 0) s% Hp_face(1:nz) = 0d0

         call set_qs(s, nz, s% q, s% dq, ierr)
         if (ierr /= 0) then
            write(*,*) 'update_vars failed in set_qs'
            return
         end if
         call set_m_and_dm(s)
         call set_dm_bar(s, s% nz, s% dm, s% dm_bar)

      end subroutine unpack_xh


      subroutine set_Teff(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: r_phot, L_surf
         logical, parameter :: skip_partials = .false., &
            need_atm_Psurf = .false., need_atm_Tsurf = .false.
         real(dp) :: Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         ierr = 0
         call set_Teff_info_for_eqns(s, skip_partials, &
            need_atm_Psurf, need_atm_Tsurf, r_phot, L_surf, Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
      end subroutine set_Teff


      subroutine set_Teff_info_for_eqns(s, skip_partials, &
            need_atm_Psurf, need_atm_Tsurf, r_surf, L_surf, Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
         use star_utils, only: set_phot_info
         use atm_lib, only: atm_Teff
         type (star_info), pointer :: s
         logical, intent(in) :: skip_partials, &
            need_atm_Psurf, need_atm_Tsurf
         real(dp), intent(out) :: r_surf, L_surf, Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr

         include 'formats'

         ierr = 0
         
         r_surf = s% r(1)
         L_surf = s% L(1)

         s% P_surf = s% Peos(1)
         s% T_surf = s% T(1)

         call set_phot_info(s) ! sets Teff using L_phot and R_phot
         Teff = s% Teff

         if (s% RSP_flag) then
            lnT_surf = s% lnT(1)
            dlnT_dL = 0d0
            dlnT_dlnR = 0d0
            dlnT_dlnM = 0d0
            dlnT_dlnkap = 0d0
            lnP_surf = s% lnPeos(1)
            dlnP_dL = 0d0
            dlnP_dlnR = 0d0
            dlnP_dlnM = 0d0
            dlnP_dlnkap = 0d0
            return
         end if
         
         if (s% use_other_surface_PT) then
            call s% other_surface_PT( &
               s% id, skip_partials, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
         else
            call get_surf_PT( &
               s, skip_partials, need_atm_Psurf, need_atm_Tsurf, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
         end if
         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*,*) 'error in get_surf_PT'
            end if
            return
         end if
         s% T_surf = exp(lnT_surf)
         s% P_surf = exp(lnP_surf)

         call set_phot_info(s) ! s% T_surf might have changed so call again

      end subroutine set_Teff_info_for_eqns


      subroutine set_hydro_vars( &
            s, nzlo, nzhi, skip_basic_vars, skip_micro_vars, skip_m_grav_and_grav, &
            skip_eos, skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_mlt, ierr)
         use atm_lib, only: atm_eval_T_tau_dq_dtau
         use atm_support, only: get_T_tau_id
         use micro, only: set_micro_vars
         use turb_info, only: set_mlt_vars, check_for_redo_MLT
         use star_utils, only: start_time, update_time, &
            set_m_grav_and_grav, set_scale_height, get_tau, &
            set_abs_du_div_cs, set_conv_time_scales
         use hydro_rotation, only: set_rotation_info, compute_j_fluxes_and_extra_jdot
         use brunt, only: do_brunt_B, do_brunt_N2
         use mix_info, only: set_mixing_info
         use hydro_rsp2, only: set_RSP2_vars

         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: &
            skip_basic_vars, skip_micro_vars, skip_m_grav_and_grav, &
            skip_eos, skip_net, skip_neu, skip_kap, skip_brunt, &
            skip_grads, skip_rotation, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_mlt
         integer, intent(out) :: ierr

         integer :: nz, num_nse, k, T_tau_id
         integer(8) :: time0
         logical, parameter :: dbg = .false.
         real(dp) :: total

         include 'formats'

         if (dbg) write(*,*) 'set_hydro_vars', nzlo, nzhi
         if (s% doing_timing) call start_time(s, time0, total)

         ierr = 0
         nz = s% nz

         if (.not. skip_basic_vars) then
            if (dbg) write(*,*) 'call set_basic_vars'
            call set_basic_vars(s, nzlo, nzhi, ierr)
            if (failed('set_basic_vars')) return
         end if

         if (.not. skip_micro_vars) then
            if (dbg) write(*,*) 'call set_micro_vars'
            call set_micro_vars(s, nzlo, nzhi, &
               skip_eos, skip_net, skip_neu, skip_kap, ierr)
            if (failed('set_micro_vars')) return
         end if

         if (.not. skip_other_cgrav) call set_cgrav(s, ierr)
         
         call get_tau(s, ierr)
         if (failed('get_tau')) return

         if (.not. skip_m_grav_and_grav) then
            ! don't change m_grav or grav during solver iteratons
            if (dbg) write(*,*) 'call set_m_grav_and_grav'
            call set_m_grav_and_grav(s)
         end if

         if (dbg) write(*,*) 'call set_scale_height'
         call set_scale_height(s)

         if (s% rotation_flag .and. (.not. skip_rotation .or. s% w_div_wc_flag)) then
            if (dbg) write(*,*) 'call set_rotation_info'
            call set_rotation_info(s, .true., ierr)
            if (failed('set_rotation_info')) return
         end if

         if (.not. skip_grads) then
            if (dbg) write(*,*) 'call do_brunt_B'
            call do_brunt_B(s, nzlo, nzhi, ierr) ! for unsmoothed_brunt_B
            if (failed('do_brunt_B')) return
            if (dbg) write(*,*) 'call set_grads'
            call set_grads(s, ierr)
            if (failed('set_grads')) return
            call set_conv_time_scales(s) ! uses brunt_B
         end if

         if (.not. skip_mixing_info) then         
            if (.not. s% RSP2_flag) then
               if (dbg) write(*,*) 'call other_adjust_mlt_gradT_fraction'
               call s% other_adjust_mlt_gradT_fraction(s% id,ierr)
               if (failed('other_adjust_mlt_gradT_fraction')) return
            end if         
            if (dbg) write(*,*) 'call set_abs_du_div_cs'
            call set_abs_du_div_cs(s)
         end if
         
         if (.not. skip_mlt .and. .not. s% RSP_flag) then
         
            if (.not. skip_mixing_info) then
               if (s% make_gradr_sticky_in_solver_iters) then
                  s% fixed_gradr_for_rest_of_solver_iters(nzlo:nzhi) = .false.   
               end if         
               s% alpha_mlt(nzlo:nzhi) = s% mixing_length_alpha
               if (s% use_other_alpha_mlt) then
                  call s% other_alpha_mlt(s% id, ierr)
                  if (ierr /= 0) then
                     if (s% report_ierr .or. dbg) &
                        write(*,*) 'other_alpha_mlt returned ierr', ierr
                     return
                  end if
               end if
            end if
            
            if (s% use_other_gradr_factor) then
               if (dbg) write(*,*) 'call other_gradr_factor'
               call s% other_gradr_factor(s% id, ierr)
               if (failed('other_gradr_factor')) return
            else if (s% use_T_tau_gradr_factor) then
               if (dbg) write(*,*) 'call use_T_tau_gradr_factor'
               call get_T_tau_id(s% atm_T_tau_relation, T_tau_id, ierr)
               do k = nzlo, nzhi
                  ! slightly inconsistent to use s% tau, defined at
                  ! faces, with s% gradr, which is a cell average
                  ! but this performs better than using
                  !    taumid = (tau(k)+tau(k+1))/2
                  call atm_eval_T_tau_dq_dtau(T_tau_id, s% tau(k), s% gradr_factor(k), ierr)
                  s% gradr_factor(k) = 1.0_dp + s% gradr_factor(k)
               end do
               if (failed('use_T_tau_gradr_factor')) return
            else
               s% gradr_factor(nzlo:nzhi) = 1d0
            end if
            
            call set_mlt_vars(s, nzlo, nzhi, ierr)
            if (failed('set_mlt_vars')) return
            if (dbg) write(*,*) 'call check_for_redo_MLT'
            
            call check_for_redo_MLT(s, nzlo, nzhi, ierr)
            if (failed('check_for_redo_MLT')) return
            
         end if

         if (.not. skip_brunt) then ! skip_brunt during solver iterations
            if (dbg) write(*,*) 'call do_brunt_N2'
            call do_brunt_N2(s, nzlo, nzhi, ierr)
            if (failed('do_brunt_N2')) return
         end if

         if (.not. skip_mixing_info) then
            if (dbg) write(*,*) 'call set_mixing_info'
            call set_mixing_info(s, skip_set_cz_bdy_mass, ierr)
            if (ierr /= 0) return
         end if

         if (s% j_rot_flag) then
            call compute_j_fluxes_and_extra_jdot(s% id, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in compute_j_fluxes'
            end if
         end if
         
         if (s% RSP2_flag) then
            call set_RSP2_vars(s,ierr)
            if (ierr /= 0) then
               if (len_trim(s% retry_message) == 0) s% retry_message = 'set_RSP2_vars failed'
               if (s% report_ierr) write(*,*) 'ierr from set_RSP2_vars'
               return
            end if
         end if

         if (s% doing_timing) &
            call update_time(s, time0, total, s% time_set_hydro_vars)
         
         s% need_to_setvars = .false.

         
         contains

         logical function failed(str)
            character (len=*), intent(in) :: str
            if (ierr == 0) then
               failed = .false.
               return
            end if
            if (s% report_ierr .or. dbg) &
               write(*,*) 'set_hydro_vars failed in call to ' // trim(str)
            failed = .true.
         end function failed
         
      end subroutine set_hydro_vars


      subroutine check_rs(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k
         logical :: okay
         include 'formats'
         ierr = 0
         okay = .true.
         do k=2, s% nz
            if (s% r(k) > s% r(k-1)) then
               if (s% report_ierr) then
                  write(*,2) 's% r(k) > s% r(k-1)', k, &
                     s% r(k)/Rsun, s% r(k-1)/Rsun, s% r(k)/Rsun-s% r(k-1)/Rsun
               end if
               okay = .false.
            end if
         end do
         if (okay) return
         s% retry_message = 'invalid values for r'
         ierr = -1
         if (s% report_ierr) write(*,*)
      end subroutine check_rs


      subroutine set_basic_vars(s, nzlo, nzhi, ierr)
         use chem_def, only: ini56
         use star_utils, only: set_rv_info, set_rmid
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: j, k, species, nz
         real(dp) :: twoGmrc2, r2, alfa, beta, sum_xa, v, u00, um1, du

         include 'formats'

         if (dbg) write(*,4) 'enter set_basic_vars: nzlo, nzhi, nz', &
            nzlo, nzhi, s% nz
         ierr = 0
         nz = s% nz
         species = s% species
         s% L_phot = s% L(1)/Lsun
         if (.not. s% using_velocity_time_centering) then
            s% d_vc_dv = 1d0
         else
            s% d_vc_dv = 0.5d0
         end if

!$OMP PARALLEL DO PRIVATE(j,k,twoGmrc2,sum_xa) SCHEDULE(dynamic,2)
         do k=nzlo, nzhi
            if (s% lnd_start(k) < -1d90) then
               s% lnd_start(k) = s% lnd(k)
               s% rho_start(k) = s% rho(k)
            end if
            if (s% T_start(k) < 0d0) then
               s% T_start(k) = s% T(k)
               s% lnT_start(k) = s% lnT(k)
            end if
            if (s% v_flag) then
               if (s% v_start(k) < -1d90) s% v_start(k) = s% v(k)
            end if
            if (s% u_flag) then
               if (s% u_start(k) < -1d90) s% u_start(k) = s% u(k)
            end if
            if (s% RTI_flag) then
               if (s% alpha_RTI_start(k) < -1d90) &
                  s% alpha_RTI_start(k) = s% alpha_RTI(k)
            end if
            if (s% RSP_flag) then
               s% RSP_w(k) = sqrt(s% RSP_Et(k))
               if (s% RSP_w_start(k) < -1d90) then
                  s% RSP_w_start(k) = s% RSP_w(k)
               end if
            end if
            if (s% r_start(k) < 0) s% r_start(k) = s% r(k)
            call set_rv_info(s,k)        
            do j=1,species
               s% xa(j,k) = max(0d0, min(1d0, s% xa(j,k)))
            end do
            sum_xa = sum(s% xa(1:species,k))
            if (abs(sum_xa - 1d0) > 1d-12) then
               do j=1,species
                  s% xa(j,k) = s% xa(j,k)/sum_xa
               end do
            end if
         end do
!$OMP END PARALLEL DO

         call set_rmid(s, nzlo, nzhi, ierr)

      end subroutine set_basic_vars
      
      
      subroutine set_cgrav(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         if (s% use_other_cgrav) then
            call s% other_cgrav(s% id, ierr)
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*,*) 'other_cgrav returned ierr', ierr
               return
            end if
         else
            s% cgrav(1:s% nz) = standard_cgrav
         end if
      end subroutine set_cgrav


      subroutine get_surf_PT( &
            s, skip_partials, &
            need_atm_Psurf, need_atm_Tsurf, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)

         use atm_support, only: get_atm_PT
         use atm_def, only: star_debugging_atm_flag, &
            atm_test_partials_val, atm_test_partials_dval_dx
         use chem_def
         use eos_lib, only: Radiation_Pressure

         type (star_info), pointer :: s
         logical, intent(in) :: skip_partials, &
            need_atm_Psurf, need_atm_Tsurf
         real(dp), intent(out) :: &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr

         real(dp) :: L_surf
         real(dp) :: R_surf
         real(dp) :: Teff
         real(dp) :: tau_surf
         real(dp) :: Teff4
         real(dp) :: T_surf4
         real(dp) :: T_surf
         real(dp) :: P_surf_atm
         real(dp) :: P_surf

         include 'formats'

         ! Set up stellar surface parameters
         
         L_surf = s% L(1)
         R_surf = s% r(1)
         Teff = s% Teff
         
         ! Initialize partials
         dlnT_dL = 0._dp; dlnT_dlnR = 0._dp; dlnT_dlnM = 0._dp; dlnT_dlnkap = 0._dp
         dlnP_dL = 0._dp; dlnP_dlnR = 0._dp; dlnP_dlnM = 0._dp; dlnP_dlnkap = 0._dp

         ! Evaluate the surface optical depth

         tau_surf = s% tau_factor*s% tau_base ! tau at outer edge of cell 1
         if (is_bad(tau_surf)) then
            write(*,1) 's% tau_factor', s% tau_factor
            write(*,1) 's% tau_base', s% tau_base
            write(*,1) 'tau_surf', tau_surf
            stop 'bad tau_surf in get_surf_PT'
         end if

         ! Evaluate surface temperature and pressure
             
         if (.not. (need_atm_Psurf .or. need_atm_Tsurf)) then

            ! Special-case boundary condition

            lnP_surf = s% lnPeos_start(1)
            if (is_bad(lnP_surf)) lnP_surf = 0._dp
            T_surf = s% T_start(1)
            lnT_surf = log(T_surf)
            T_surf4 = T_surf*T_surf*T_surf*T_surf
            Teff = pow(four_thirds*T_surf4/(tau_surf + two_thirds), 0.25_dp)

            if (.not. skip_partials) then
               dlnT_dL = 0._dp; dlnT_dlnR = 0._dp; dlnT_dlnM = 0._dp; dlnT_dlnkap = 0._dp
               dlnP_dL = 0._dp; dlnP_dlnR = 0._dp; dlnP_dlnM = 0._dp; dlnP_dlnkap = 0._dp
            endif

         else

            ! Evaluate temperature and pressure based on atm_option
               
            if (.false. .and. 1 == s% solver_test_partials_k .and. &
                  s% solver_iter == s% solver_test_partials_iter_number) then
               star_debugging_atm_flag = .true.
            end if

            call get_atm_PT( &
                 s, tau_surf, L_surf, R_surf, s% m(1), s% cgrav(1), skip_partials, &
                 Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
                 lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
                 ierr)
            if (ierr /= 0) then
               if (s% report_ierr) then
                  write(*,1) 'tau_surf', tau_surf
                  write(*,1) 'L_surf', L_surf
                  write(*,1) 'R_surf', R_surf
                  write(*,1) 'Teff', Teff
                  write(*,1) 's% m(1)', s% m(1)
                  write(*,1) 's% cgrav(1)', s% cgrav(1)
                  write(*,*) 'failed in get_atm_PT'
               end if
               return
            end if

            if (.false. .and. 1 == s% solver_test_partials_k .and. &
                  s% solver_iter == s% solver_test_partials_iter_number) then
               s% solver_test_partials_val = atm_test_partials_val
               s% solver_test_partials_dval_dx = atm_test_partials_dval_dx
            end if

         end if

         ! Check outputs
      
         if (is_bad(lnT_surf) .or. is_bad(lnP_surf)) then
            if (len_trim(s% retry_message) == 0) s% retry_message = 'bad logT surf or logP surf'
            ierr = -1
            write(*,1) 'bad outputs in get_surf_PT'
            write(*,1) 'lnT_surf', lnT_surf
            write(*,1) 'lnP_surf', lnP_surf
            write(*,*) 'atm_option = ', trim(s% atm_option)
            if (s% stop_for_bad_nums) stop 'get PT surf'
         end if

         ! Finish

         return

      end subroutine get_surf_PT


      subroutine set_grads(s, ierr)
         use chem_def, only: chem_isos
         use star_utils, only: smooth, safe_div_val
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, nz, j, cid, max_cid
         real(dp) :: val, max_val, A, Z
         real(dp), pointer, dimension(:) :: dlnP, dlnd, dlnT

         include 'formats'

         ierr = 0
         nz = s% nz
         call do_alloc(ierr)
         if (ierr /= 0) return

         do k = 2, nz
            dlnP(k) = s% lnPeos(k-1) - s% lnPeos(k)
            dlnd(k) = s% lnd(k-1) - s% lnd(k)
            dlnT(k) = s% lnT(k-1) - s% lnT(k)
         end do
         dlnP(1) = dlnP(2)
         dlnd(1) = dlnd(2)
         dlnT(1) = dlnT(2)

         call smooth(dlnP,nz)
         call smooth(dlnd,nz)
         call smooth(dlnT,nz)

         s% grad_density(1) = 0
         s% grad_temperature(1) = 0
         do k = 2, nz
            if (dlnP(k) >= 0) then
               s% grad_density(k) = 0
               s% grad_temperature(k) = 0
            else
               s% grad_density(k) = safe_div_val(s, dlnd(k), dlnP(k))
               s% grad_temperature(k) = safe_div_val(s, dlnT(k), dlnP(k))
            end if
         end do

         call smooth(s% grad_density,nz)
         call smooth(s% grad_temperature,nz)

         if (s% use_Ledoux_criterion .and. s% calculate_Brunt_B) then
            do k=1,nz
               s% gradL_composition_term(k) = s% unsmoothed_brunt_B(k)
            end do
            call smooth_gradL_composition_term
         else
            do k=1,nz
               s% gradL_composition_term(k) = 0d0
            end do
         end if

         call dealloc

         do k=3,nz-2
            max_cid = 0
            max_val = -1d99
            do j=1,s% species
               cid = s% chem_id(j)
               A = dble(chem_isos% Z_plus_N(cid))
               Z = dble(chem_isos% Z(cid))
               val = (s% xa(j,k-2) + s% xa(j,k-1) - s% xa(j,k) - s% xa(j,k+1))*(1d0 + Z)/A
               if (val > max_val) then
                  max_val = val
                  max_cid = cid
               end if
            end do
            s% dominant_iso_for_thermohaline(k) = max_cid
         end do
         s% dominant_iso_for_thermohaline(1:2) = &
            s% dominant_iso_for_thermohaline(3)
         s% dominant_iso_for_thermohaline(nz-1:nz) = &
            s% dominant_iso_for_thermohaline(nz-2)


         contains

         subroutine smooth_gradL_composition_term
            use star_utils, only: weighed_smoothing, threshold_smoothing
            logical, parameter :: preserve_sign = .false.
            real(dp), pointer, dimension(:) :: work
            integer :: k
            include 'formats'
            ierr = 0
            work => dlnd
            if (s% num_cells_for_smooth_gradL_composition_term <= 0) return
            call threshold_smoothing( &
               s% gradL_composition_term, s% threshold_for_smooth_gradL_composition_term, s% nz, &
               s% num_cells_for_smooth_gradL_composition_term, preserve_sign, work)
         end subroutine smooth_gradL_composition_term

         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            call do_work_arrays(.true.,ierr)
         end subroutine do_alloc

         subroutine dealloc
            call do_work_arrays(.false.,ierr)
         end subroutine dealloc

         subroutine do_work_arrays(alloc_flag, ierr)
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               dlnP, nz, nz_alloc_extra, 'mlt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dlnd, nz, nz_alloc_extra, 'mlt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dlnT, nz, nz_alloc_extra, 'mlt', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays

      end subroutine set_grads

   
      end module hydro_vars

