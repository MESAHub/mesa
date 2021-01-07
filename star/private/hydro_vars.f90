! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
      public :: set_vars_if_needed, set_vars, set_final_vars, set_cgrav, &
         set_hydro_vars, set_Teff_info_for_eqns, set_Teff, get_surf_PT

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
            skip_time_derivatives = .false., &
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
            s, skip_time_derivatives, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat, &
            dt, ierr)
      end subroutine set_vars

      
      subroutine set_final_vars(s, dt, ierr)
         use alloc, only: non_crit_get_work_array, non_crit_return_work_array
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
         real(dp), pointer, dimension(:) :: &
            d_mx, cv, nu_st, d_st, d_dsi, d_sh, d_ssi, d_es, d_gsf, rti_info
         
         include 'formats'
         
         ierr = 0
         nz = s% nz

         nullify(d_mx, cv, nu_st, d_st, d_dsi, d_sh, d_ssi, d_es, d_gsf, rti_info)

         ! save and restore mixing coeffs needed for time smoothing
         if (s% rotation_flag) then
            call get_cpy(s% nu_ST, nu_st, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_ST, d_st, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_DSI, d_dsi, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_SH, d_sh, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_SSI, d_ssi, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_ES, d_es, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_GSF, d_gsf, ierr1)
            if (ierr1 /= 0) ierr = -1
         end if

         if (s% RTI_flag) then
            call get_cpy(s% dPdr_dRhodr_info, RTI_info, ierr1)
            if (ierr1 /= 0) ierr = -1
         end if

         call get_cpy(s% D_mix, d_mx, ierr1)
         if (ierr1 /= 0) ierr = -1
         
         if (ierr /= 0) return
         
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

         call set_Teff(s, ierr)
         if (ierr /= 0) then
            if (s% report_ierr .or. dbg) &
               write(*,*) 'update_vars: set_Teff returned ierr', ierr
            return
         end if
         
         if (s% op_split_burn) then
            do k = 1, nz
               if (s% T_start(k) >= s% op_split_burn_min_T) &
                  s% eps_nuc(k) = s% burn_avg_epsnuc(k)
            end do
         end if

         contains


         subroutine get_cpy(src,cpy,ierr)
            real(dp), pointer, dimension(:) :: src, cpy
            integer, intent(out) :: ierr
            integer :: k
            ierr = 0
            call non_crit_get_work_array( &
               s, cpy, nz, nz_alloc_extra, 'set_final_vars', ierr)
            if (ierr /= 0) return
            do k=1,nz
               cpy(k) = src(k)
            end do
         end subroutine get_cpy


         subroutine restore(src,cpy)
            real(dp), pointer, dimension(:) :: src, cpy
            integer :: k
            if (.not. associated(cpy)) return
            do k=1,nz
               src(k) = cpy(k)
            end do
            call non_crit_return_work_array(s, cpy, 'set_final_vars')
         end subroutine restore


      end subroutine set_final_vars


      subroutine set_some_vars( &
            s, skip_time_derivatives, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat, &
            dt, ierr)
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s
         logical, intent(in) :: &
            skip_time_derivatives, skip_m_grav_and_grav, &
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
            skip_time_derivatives, skip_basic_vars, skip_micro_vars, &
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
            skip_time_derivatives, skip_basic_vars, skip_micro_vars, &
            skip_m_grav_and_grav, skip_net, skip_neu, skip_kap, &
            skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat, &
            skip_mlt, skip_eos, dt, ierr)
         use star_utils, only: eval_irradiation_heat, set_qs, set_dm_bar, set_m_and_dm
         type (star_info), pointer :: s
         logical, intent(in) :: &
            skip_time_derivatives, skip_basic_vars, skip_micro_vars, &
            skip_m_grav_and_grav, skip_net, skip_neu, skip_kap, skip_grads, &
            skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_irradiation_heat, &
            skip_mlt, skip_eos
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         integer :: i_lnd, i_lnT, i_lnR, i_eturb, &
            i_lum, i_v, i_u, i_alpha_RTI, i_ln_cvpv0, i_eturb_RSP, &
            j, k, species, nvar_chem, nz, k_below_just_added
         real(dp) :: dt_inv

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
         i_eturb = s% i_eturb
         i_v = s% i_v
         i_u = s% i_u
         i_alpha_RTI = s% i_alpha_RTI
         i_eturb_RSP = s% i_eturb_RSP
         i_ln_cvpv0 = s% i_ln_cvpv0

         if (s% doing_finish_load_model .or. .not. s% RSP_flag) then
            call unpack
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*,*) 'after unpack ierr', ierr
               return
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

         if (s% Teff <= 0 .and. .not. s% doing_finish_load_model) then
            call set_Teff(s, ierr)
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*,*) 'update_vars: set_Teff returned ierr', ierr
               return
            end if
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
         
         contains
         
         subroutine unpack
            integer :: j, k
            include 'formats'
         
            do j=1,s% nvar_hydro
               if (j == i_lnd) then
                  do k=1,nz
                     s% lnd(k) = s% xh(i_lnd,k)
                     s% dxh_lnd(k) = 0d0
                  end do
               else if (j == i_lnT) then
                  do k=1,nz
                     s% lnT(k) = s% xh(i_lnT,k)
                     s% dxh_lnT(k) = 0d0
                  end do
               else if (j == i_lnR) then
                  do k=1,nz
                     s% lnR(k) = s% xh(i_lnR,k)
                     s% dxh_lnR(k) = 0d0
                  end do
               else if (j == i_eturb) then
                  do k=1,nz
                     s% Eturb(k) = max(s% xh(i_eturb, k), min_eturb)
                     s% dxh_eturb(k) = 0d0
                  end do
               else if (j == i_lum) then
                  do k=1,nz
                     s% L(k) = s% xh(i_lum, k)
                  end do
               else if (j == i_v) then
                  do k=1,nz
                     s% v(k) = s% xh(i_v,k)
                  end do
               else if (j == i_u) then
                  do k=1,nz
                     s% u(k) = s% xh(i_u,k)
                  end do
               else if (j == i_alpha_RTI) then
                  do k=1,nz
                     s% alpha_RTI(k) = max(0d0, s% xh(i_alpha_RTI,k))
                  end do
               else if (j == i_eturb_RSP) then
                  do k=1,nz
                     s% Et(k) = max(0d0,s% xh(i_eturb_RSP,k))
                  end do
               else if (j == i_ln_cvpv0) then
                  do k=1,nz
                     s% conv_vel(k) = max(0d0, exp(s% xh(i_ln_cvpv0,k))-s% conv_vel_v0)
                     s% dxh_ln_cvpv0(k) = 0d0
                  end do
               end if
            end do

            if (s% trace_k > 0 .and. s% trace_k <= nz) then
               k = s% trace_k
               if (i_lnd /= 0) &
                  write(*,5) 'update_vars: lnd', k, s% solver_iter, s% solver_adjust_iter, &
                           s% model_number, s% lnd(k)
               if (i_lnT /= 0) &
                  write(*,5) 'update_vars: lnT', k, s% solver_iter, s% solver_adjust_iter, &
                              s% model_number, s% lnT(k)
               if (i_lum /= 0) &
                  write(*,5) 'update_vars: L', k, s% solver_iter, s% solver_adjust_iter, &
                           s% model_number, s% L(k)
               if (i_lnR /= 0) &
                  write(*,5) 'update_vars: lnR', k, s% solver_iter, s% solver_adjust_iter, &
                        s% model_number, s% lnR(k)
               if (i_v /= 0) &
                  write(*,5) 'update_vars: v', k, s% solver_iter, s% solver_adjust_iter, &
                           s% model_number, s% v(k)
               if (i_u /= 0) &
                  write(*,5) 'update_vars: u', k, s% solver_iter, s% solver_adjust_iter, &
                           s% model_number, s% u(k)
               if (i_alpha_RTI /= 0) &
                  write(*,5) 'update_vars: alpha_RTI', k, s% solver_iter, s% solver_adjust_iter, &
                           s% model_number, s% alpha_RTI(k)
               if (i_eturb_RSP /= 0) &
                  write(*,5) 'update_vars: Et', k, s% solver_iter, s% solver_adjust_iter, &
                           s% model_number, s% Et(k)
               if (i_ln_cvpv0 /= 0) &
                  write(*,5) 'update_vars: conv_vel', k, s% solver_iter, s% solver_adjust_iter, &
                           s% model_number, s% conv_vel(k)
            end if

            if (i_lum == 0 .and. .not. (s% RSP_flag .or. s% Eturb_flag)) s% L(1:nz) = 0d0

            if (i_v == 0) s% v(1:nz) = 0d0

            if (i_u == 0) s% u(1:nz) = 0d0

            call set_qs(s, nz, s% q, s% dq, ierr)
            if (ierr /= 0) then
               write(*,*) 'update_vars failed in set_qs'
               return
            end if
            call set_m_and_dm(s)
            call set_dm_bar(s, s% nz, s% dm, s% dm_bar)

            if (.not. skip_mixing_info) then
               s% mixing_type(1:nz) = no_mixing
               s% adjust_mlt_gradT_fraction(1:nz) = -1
            end if

            if (.not. skip_time_derivatives) then

               ! time derivatives at constant q for use by eps_grav and convective velocities
               if (s% generations < 2 .or. dt <= 0 .or. s% nz /= s% nz_old) then

                  if (i_lnd /= 0) s% dlnd_dt_const_q(1:nz) = 0
                  if (i_lnT /= 0) s% dlnT_dt_const_q(1:nz) = 0
                  if (i_ln_cvpv0 /= 0) s% dln_cvpv0_dt_const_q(1:nz) = 0

               else

                  dt_inv = 1/dt
                  if (i_lnd /= 0) then
                     do k=1,nz
                        s% dlnd_dt_const_q(k) = &
                           (s% xh(i_lnd,k) - s% lnd_for_d_dt_const_q(k))*dt_inv
                        if (is_bad(s% dlnd_dt_const_q(k))) then
                           s% retry_message = 'update_vars: bad dlnd_dt_const_q'
                           ierr = -1
                           if (s% report_ierr) &
                              write(*,2) 'update_vars: bad dlnd_dt_const_q', k, &
                                 s% dlnd_dt_const_q(k), &
                                 s% lnd_for_d_dt_const_q(k), s% xh(i_lnd,k)
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad dlnd_dt_const_q', k, s% dlnd_dt_const_q(k)
                              write(*,2) 's% lnd_for_d_dt_const_q(k)', k, s% lnd_for_d_dt_const_q(k)
                              write(*,2) 's% xh(i_lnd,k)', k, s% xh(i_lnd,k)
                              stop 'update_vars'
                           end if
                           return
                           stop 'update_vars'
                        end if
                     end do
                  end if

                  if (i_lnT /= 0) then
                     do k=1,nz
                        s% dlnT_dt_const_q(k) = &
                           (s% xh(i_lnT,k) - s% lnT_for_d_dt_const_q(k))*dt_inv
                        if (is_bad(s% dlnT_dt_const_q(k))) then
                           s% retry_message = 'update_vars: bad dlnT_dt_const_q'
                           ierr = -1
                           if (s% report_ierr) &
                              write(*,2) 'update_vars: bad dlnT_dt_const_q', k, &
                                 s% dlnT_dt_const_q(k), &
                                 s% lnT_for_d_dt_const_q(k), s% xh(i_lnT,k)
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad dlnT_dt_const_q', k, &
                                 s% dlnT_dt_const_q(k), &
                                 s% lnT_for_d_dt_const_q(k), s% xh(i_lnT,k)
                              stop 'update_vars'
                           end if
                           return
                        end if
                     end do
                  end if

                  if (i_ln_cvpv0 /= 0) then
                     do k=1,nz
                        s% dln_cvpv0_dt_const_q(k) = &
                           (s% xh(i_ln_cvpv0,k) - s% ln_cvpv0_for_d_dt_const_q(k))*dt_inv
                        if (is_bad(s% dln_cvpv0_dt_const_q(k))) then
                           ierr = -1
                           s% retry_message = 'update_vars: bad dln_cvpv0_dt_const_q'
                           if (s% report_ierr) &
                              write(*,2) 'update_vars: bad dln_cvpv0_dt_const_q', k, &
                                 s% dln_cvpv0_dt_const_q(k), &
                                 s% ln_cvpv0_for_d_dt_const_q(k), s% xh(i_ln_cvpv0,k)
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad dln_cvpv0_dt_const_q', k, &
                                 s% dln_cvpv0_dt_const_q(k), &
                                 s% ln_cvpv0_for_d_dt_const_q(k), s% xh(i_ln_cvpv0,k)
                              stop 'update_vars'
                           end if
                           return
                        end if
                     end do
                  end if

               end if

               ! time derivatives at constant mass
               if (dt > 0d0) then
                  dt_inv = 1/dt
                  s% dVARdot_dVAR = dt_inv
               else
                  s% dVARdot_dVAR = 0
                  dt_inv = 0
               end if
            
               if (s% generations < 2 .or. dt <= 0 .or. s% nz /= s% nz_old) then

                  if (i_lnd /= 0) s% dlnd_dt(1:nz) = 0
                  if (i_lnT /= 0) s% dlnT_dt(1:nz) = 0
                  if (i_eturb /= 0) s% deturb_dt(1:nz) = 0
                  if (i_lnR /= 0) s% dlnR_dt(1:nz) = 0
                  if (s% v_flag) s% dv_dt(1:nz) = 0
                  if (s% u_flag) s% du_dt(1:nz) = 0
                  if (s% RTI_flag) s% dalpha_RTI_dt(1:nz) = 0
                  if (s% RSP_flag) s% dEt_dt(1:nz) = 0
                  if (s% conv_vel_flag) s% dln_cvpv0_dt(1:nz) = 0

               else

                  do k=1,k_below_just_added-1
                     if (i_lnR /= 0) s% dlnR_dt(k) = 0
                     if (i_v /= 0) s% dv_dt(k) = 0
                     if (i_u /= 0) s% du_dt(k) = 0
                     if (i_lnd /= 0) s% dlnd_dt(k) = 0
                     if (i_lnT /= 0) s% dlnT_dt(k) = 0
                     if (i_eturb /= 0) s% deturb_dt(k) = 0
                     if (i_eturb_RSP /= 0) s% dEt_dt(k) = 0
                     if (i_alpha_RTI /= 0) s% dalpha_RTI_dt(k) = 0
                     if (i_ln_cvpv0 /= 0) s% dln_cvpv0_dt(k) = 0
                  end do

                  if (i_lnd /= 0) then
                     do k=k_below_just_added,nz
                        s% dlnd_dt(k) = &
                           (s% xh(i_lnd,k) - s% lnd_for_d_dt_const_m(k))*dt_inv
                        if (is_bad(s% dlnd_dt(k))) then
                           s% retry_message = 'update_vars: bad dlnd_dt'
                           ierr = -1
                           if (s% report_ierr) &
                              write(*,2) 'update_vars: bad dlnd_dt', k, &
                                 s% dlnd_dt(k), s% lnd_for_d_dt_const_m(k), s% xh(i_lnd,k)
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad dlnd_dt', k, &
                                 s% dlnd_dt(k)
                              stop 'update_vars'
                           end if
                           return
                           stop 'update_vars'
                        end if
                     end do
                  end if

                  if (i_lnT /= 0) then
                     do k=k_below_just_added,nz
                        s% dlnT_dt(k) = (s% xh(i_lnT,k) - s% lnT_for_d_dt_const_m(k))*dt_inv
                        if (is_bad(s% dlnT_dt(k))) then
                           s% retry_message = 'update_vars: bad dlnT_dt'
                           ierr = -1
                           if (s% report_ierr) &
                              write(*,2) 'update_vars: bad dlnT_dt', k, s% dlnT_dt(k), &
                                 s% xh(i_lnT,k), s% lnT_for_d_dt_const_m(k), dt_inv
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad dlnT_dt', k, &
                                 s% dlnT_dt(k)
                              stop 'update_vars'
                           end if
                        end if
                     end do
                  end if

                  if (s% v_flag) then
                     do k=k_below_just_added,nz
                        s% dv_dt(k) = (s% xh(i_v,k) - s% v_for_d_dt_const_m(k))*dt_inv
                        if (is_bad(s% dv_dt(k))) then
                           s% retry_message = 'update_vars: bad dv_dt'
                           ierr = -1
                           if (s% report_ierr) &
                              write(*,2) 'update_vars: bad dv_dt', k, s% dv_dt(k)
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad dv_dt', k, &
                                 s% dv_dt(k), s% xh(i_v,k), s% v_for_d_dt_const_m(k)
                              stop 'update_vars'
                           end if
                        end if
                     end do
                  end if

                  if (s% u_flag) then
                     do k=k_below_just_added,nz
                        s% du_dt(k) = (s% xh(i_u,k) - s% u_for_d_dt_const_m(k))*dt_inv
                        if (is_bad(s% du_dt(k))) then
                           ierr = -1
                           s% retry_message = 'update_vars: bad du_dt'
                           if (s% report_ierr) &
                              write(*,2) 'update_vars: bad du_dt', k, s% du_dt(k)
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad du_dt', k, &
                                 s% du_dt(k)
                              stop 'update_vars'
                           end if
                        end if
                     end do
                  end if

                  if (s% Eturb_flag) then
                     do k=k_below_just_added,nz
                        s% deturb_dt(k) = (s% xh(i_eturb,k) - s% eturb_for_d_dt_const_m(k))*dt_inv
                        if (is_bad(s% deturb_dt(k))) then
                           ierr = -1
                           s% retry_message = 'update_vars: bad deturb_dt'
                           if (s% report_ierr) &
                              write(*,2) 'update_vars: bad deturb_dt', k, s% deturb_dt(k)
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad deturb_dt', k, &
                                 s% du_dt(k)
                              stop 'update_vars'
                           end if
                        end if
                     end do
                  end if

                  if (s% RTI_flag) then
                     do k=k_below_just_added,nz
                        s% dalpha_RTI_dt(k) = &
                           (s% xh(i_alpha_RTI,k) - s% alpha_RTI_for_d_dt_const_m(k))*dt_inv
                        if (is_bad(s% dalpha_RTI_dt(k))) then
                           ierr = -1
                           s% retry_message = 'update_vars: bad dalpha_RTI_dt'
                           if (s% report_ierr) &
                              write(*,3) 'update_vars: bad dalpha_RTI_dt', k, k_below_just_added, &
                                 s% alpha_RTI_for_d_dt_const_m(k), s% dalpha_RTI_dt(k), &
                                 s% xh(i_alpha_RTI,k)
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad dalpha_RTI_dt', k, &
                                 s% dalpha_RTI_dt(k)
                              stop 'update_vars'
                           end if
                        end if
                     end do
                  end if

                  if (s% conv_vel_flag) then
                     do k=k_below_just_added,nz
                        s% dln_cvpv0_dt(k) = &
                           (s% xh(i_ln_cvpv0,k) - s% ln_cvpv0_for_d_dt_const_m(k))*dt_inv
                        if (is_bad(s% dln_cvpv0_dt(k))) then
                           ierr = -1
                           s% retry_message = 'update_vars: bad dln_cvpv0_dt'
                           if (s% report_ierr) &
                              write(*,3) 'update_vars: bad dln_cvpv0_dt', k, k_below_just_added, &
                                 s% ln_cvpv0_for_d_dt_const_m(k), s% dln_cvpv0_dt(k), &
                                 s% xh(i_ln_cvpv0,k), s% ln_cvpv0_for_d_dt_const_m(k)
                           if (s% stop_for_bad_nums) then
                              write(*,2) 'update_vars: bad dln_cvpv0_dt', k, &
                                 s% dln_cvpv0_dt(k)
                              stop 'update_vars'
                           end if
                        end if
                     end do
                  end if

                  do k=k_below_just_added,nz
                     s% dlnR_dt(k) = (s% xh(i_lnR,k) - s% lnR_for_d_dt_const_m(k))*dt_inv
                     if (is_bad(s% dlnR_dt(k))) then
                        ierr = -1
                        s% retry_message = 'update_vars: bad dlnR_dt'
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dlnR_dt', k, s% dlnR_dt(k)
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'update_vars: bad dlnR_dt', k, &
                              s% dlnR_dt(k)
                           stop 'update_vars'
                        end if
                     end if
                  end do

               end if

            end if

         end subroutine unpack

      end subroutine update_vars


      subroutine set_Teff(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: r_phot, L_surf
         logical, parameter :: skip_partials = .true.
         real(dp) :: Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         ierr = 0
         call set_Teff_info_for_eqns(s, skip_partials, r_phot, L_surf, Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
      end subroutine set_Teff


      subroutine set_Teff_info_for_eqns(s, skip_partials, r_surf, L_surf, Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
         use star_utils, only: get_phot_info
         type (star_info), pointer :: s
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: r_surf, L_surf, Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr

         integer :: k_phot, k
         real(dp) :: Tm1, T00, T4_m1, T4_00, P_rad_m1, P_rad_00, y_phot, &
            alfa, beta, A, opacity_face, r_phot, m_phot, v_phot, &
            L_phot, T_phot, cs_phot, kap_phot, logg_phot, dL_dlnR, dL_dlnT
         real(qp) :: q1

         include 'formats'

         ierr = 0
         
         ! Set surface values

         L_surf = s% L(1)
         r_surf = s% r(1)

         if (s% RSP_flag .and. .not. s% RSP_use_atm_grey_with_kap_for_Psurf) then

            call get_phot_info( &
                 s, r_phot, m_phot, v_phot, L_phot, T_phot, &
                 cs_phot, kap_phot, logg_phot, y_phot, k_phot)

            Teff = T_phot
            lnT_surf = s% lnT(1)
            dlnT_dL = 0d0
            dlnT_dlnR = 0d0
            dlnT_dlnM = 0d0
            dlnT_dlnkap = 0d0
            lnP_surf = s% lnP(1)
            dlnP_dL = 0d0
            dlnP_dlnR = 0d0
            dlnP_dlnM = 0d0
            dlnP_dlnkap = 0d0
            s% Teff = Teff
            s% L_phot = L_phot
            s% photosphere_L = s% L_phot
            s% photosphere_r = r_phot/Rsun
            return
         end if

         if (s% use_other_surface_PT) then
            call s% other_surface_PT( &
               s% id, skip_partials, &
               Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
         else
            call get_surf_PT( &
               s, skip_partials, &
               Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
         end if

         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*,*) 'error in get_surf_PT'
            end if
            return
         end if

         s% Teff = Teff
         
         ! Calculate and store photosphere (tau=2/3) values; these
         ! aren't actually used to set up surface values

         call get_phot_info( &
              s, r_phot, m_phot, v_phot, L_phot, T_phot, &
              cs_phot, kap_phot, logg_phot, y_phot, k_phot)

         s% L_phot = L_phot/Lsun
         s% photosphere_L = s% L_phot
         s% photosphere_r = r_phot/Rsun

      end subroutine set_Teff_info_for_eqns


      subroutine set_hydro_vars( &
            s, nzlo, nzhi, skip_basic_vars, skip_micro_vars, skip_m_grav_and_grav, &
            skip_eos, skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_mlt, ierr)
         use atm_lib, only: atm_eval_T_tau_dq_dtau
         use atm_support, only: get_T_tau_id
         use micro, only: set_micro_vars
         use mlt_info, only: set_mlt_vars, check_for_redo_MLT, set_grads
         use star_utils, only: set_k_CpTMdot_lt_L, start_time, update_time, &
            set_m_grav_and_grav, set_scale_height, get_tau, &
            set_abs_du_div_cs
         use hydro_rotation, only: set_rotation_info, compute_j_fluxes_and_extra_jdot
         use hydro_eturb, only: reset_eturb_using_L
         use brunt, only: do_brunt_B, do_brunt_N2
         use mix_info, only: set_mixing_info

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

         if (dbg) write(*,*) 'call set_k_CpTMdot_lt_L'
         call set_k_CpTMdot_lt_L(s)

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
         end if

         if (.not. skip_mixing_info) then
         
            if (.not. s% Eturb_flag) then
               if (dbg) write(*,*) 'call other_adjust_mlt_gradT_fraction'
               call s% other_adjust_mlt_gradT_fraction(s% id,ierr)
               if (failed('other_adjust_mlt_gradT_fraction')) return
            end if
         
            if (s% u_flag) then
               if (dbg) write(*,*) 'call set_abs_du_div_cs'
               call set_abs_du_div_cs(s)
            end if

         end if
         
         if (.not. skip_mlt .and. .not. s% RSP_flag .and. .not. s% Eturb_flag) then
         
            if (.not. skip_mixing_info) then
               if (s% make_gradr_sticky_in_solver_iters) &
                  s% fixed_gradr_for_rest_of_solver_iters(nzlo:nzhi) = .false.            
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
         
         if (s% need_to_reset_eturb) then
            stop 'need_to_reset_eturb is not ready for use'
            if (dbg) write(*,*) 'call set_mlt_vars'
            call set_mlt_vars(s, 1, s% nz, ierr)
            if (failed('set_mlt_vars')) return
            call reset_eturb_using_L(s,ierr)
            if (failed('reset_eturb_using_L')) return
            s% need_to_reset_eturb = .false.
         end if

         if (.not. skip_brunt) then
            if (dbg) write(*,*) 'call do_brunt_N2'
            call do_brunt_N2(s, nzlo, nzhi, ierr)
            if (failed('do_brunt_N2')) return
         end if

         if (.not. skip_mixing_info) then
            if (dbg) write(*,*) 'call set_mixing_info'
            call set_mixing_info(s, skip_set_cz_bdy_mass, ierr)
            if (ierr /= 0) return
            call set_photosphere_start_info
         end if

         if (s% j_rot_flag) then
            call compute_j_fluxes_and_extra_jdot(s% id, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in compute_j_fluxes'
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
         
         subroutine set_photosphere_start_info
            use star_utils, only: get_phot_kap
            real(dp) :: kap
            include 'formats'
            kap = get_phot_kap(s)
            s% photosphere_opacity_start = kap
         end subroutine set_photosphere_start_info

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
         if (.not. s% using_Fraley_time_centering) then
            s% d_vc_dv = 1d0
         else
            s% d_vc_dv = 0.5d0
         end if

!$OMP PARALLEL DO PRIVATE(j,k,twoGmrc2,sum_xa) SCHEDULE(dynamic,2)
         do k=nzlo, nzhi
            s% rho(k) = exp(s% lnd(k))
            if (s% lnd_start(k) < -1d90) then
               s% lnd_start(k) = s% lnd(k)
               s% rho_start(k) = s% rho(k)
            end if
            s% T(k) = exp(s% lnT(k))
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
            if (s% conv_vel_flag) then
               s% conv_vel_start(k) = s% conv_vel(k)
            end if
            if (s% RSP_flag) then
               s% w(k) = sqrt(s% Et(k))
               if (s% w_start(k) < -1d90) then
                  s% w_start(k) = s% w(k)
               end if
            end if
            s% r(k) = exp(s% lnR(k))
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
            Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)

         use atm_support, only: get_atm_PT
         use atm_def, only: star_debugging_atm_flag, &
            atm_test_partials_val, atm_test_partials_dval_dx
         use chem_def
         use eos_lib, only: Radiation_Pressure

         type (star_info), pointer :: s
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr

         real(dp) :: L_surf
         real(dp) :: R_surf
         real(dp) :: tau_surf
         real(dp) :: Teff4
         real(dp) :: T_surf4
         real(dp) :: T_surf
         real(dp) :: P_surf_atm
         real(dp) :: P_surf
         logical :: do_not_need_atm_Psurf, do_not_need_atm_Tsurf

         include 'formats'
         
         do_not_need_atm_Psurf = &
            s% use_compression_outer_BC .or. &
            s% use_zero_Pgas_outer_BC .or. &
            s% use_fixed_vsurf_outer_BC .or. &
            (s% use_momentum_outer_BC .and. trim(s% atm_option) == 'fixed_Psurf')
            
         do_not_need_atm_Tsurf = &
            s% i_lum == 0 .or. &
            ((s% use_fixed_L_for_BB_outer_BC .or. s% tau_for_L_BB > 0d0) .and. &
             (s% use_T_black_body_outer_BC))

         ! Set up stellar surface parameters

         L_surf = s% L(1)
         R_surf = s% r(1)
         
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
             
         if (do_not_need_atm_Psurf .and. do_not_need_atm_Tsurf) then

            ! Special-case boundary conditions

            lnP_surf = s% lnP_start(1)
            if (is_bad(lnP_surf)) lnP_surf = 0._dp
            T_surf = s% T_start(1)
            lnT_surf = log(T_surf)
            T_surf4 = T_surf*T_surf*T_surf*T_surf
            Teff = pow(4._dp/3._dp*T_surf4/(tau_surf + 2._dp/3._dp), 0.25_dp)

            if (.not. skip_partials) then
               dlnT_dL = 0._dp; dlnT_dlnR = 0._dp; dlnT_dlnM = 0._dp; dlnT_dlnkap = 0._dp
               dlnP_dL = 0._dp; dlnP_dlnR = 0._dp; dlnP_dlnM = 0._dp; dlnP_dlnkap = 0._dp
            endif

         else

            ! Evaluate temperature and pressure based on atm_option
            ! (yes, we might want to translate atm_option into an integer,
            ! but these string comparisons are cheap)

            ! The first few are special, 'trivial-atmosphere' options

            select case (s% atm_option)

            case ('fixed_Teff')

               ! set Tsurf from Eddington T-tau relation
               !     for current surface tau and Teff = `atm_fixed_Teff`.
               ! set Psurf = Radiation_Pressure(Tsurf)

               Teff = s% atm_fixed_Teff
               Teff4 = Teff*Teff*Teff*Teff
               T_surf = pow(3._dp/4._dp*Teff4*(tau_surf + 2._dp/3._dp), 0.25_dp)
               lnT_surf = log(T_surf)
               lnP_surf = Radiation_Pressure(T_surf)

               if (.not. skip_partials) then
                  dlnT_dL = 0._dp; dlnT_dlnR = 0._dp; dlnT_dlnM = 0._dp; dlnT_dlnkap = 0._dp
                  dlnP_dL = 0._dp; dlnP_dlnR = 0._dp; dlnP_dlnM = 0._dp; dlnP_dlnkap = 0._dp
               endif
               
            case ('fixed_Tsurf')

               ! set Teff from Eddington T-tau relation for given
               ! Tsurf and tau=2/3 set Psurf =
               ! Radiation_Pressure(Tsurf)

               T_surf = s% atm_fixed_Tsurf
               lnT_surf = log(T_surf)
               T_surf4 = T_surf*T_surf*T_surf*T_surf
               Teff = pow(4._dp/3._dp*T_surf4/(tau_surf + 2._dp/3._dp), 0.25_dp)
               lnP_surf = Radiation_Pressure(T_surf)

               if (.not. skip_partials) then
                  dlnT_dL = 0._dp; dlnT_dlnR = 0._dp; dlnT_dlnM = 0._dp; dlnT_dlnkap = 0._dp
                  dlnP_dL = 0._dp; dlnP_dlnR = 0._dp; dlnP_dlnM = 0._dp; dlnP_dlnkap = 0._dp
               endif

            case ('fixed_Psurf')

               ! set Teff from L and R using L = 4*pi*R^2*boltz_sigma*T^4.
               ! set Tsurf using Eddington T-tau relation

               if (L_surf < 0._dp) then
                  if (s% report_ierr) then
                     write(*,2) 'get_surf_PT: L_surf <= 0', s% model_number, L_surf
                     call mesa_error(__FILE__,__LINE__)
                  end if
                  s% retry_message = 'L_surf < 0'
                  ierr = -1
                  return
               end if

               lnP_surf = safe_log(s% atm_fixed_Psurf)
               if (R_surf <= 0._dp) then
                  T_surf4 = 1._dp
                  T_surf = 1._dp
                  lnT_surf = 0._dp
                  if (.not. skip_partials) then
                     dlnT_dlnR = 0._dp
                     dlnT_dL = 0._dp
                  endif
                  Teff = s% T(1)
               else
                  Teff = pow(L_surf/(4._dp*pi*R_surf*R_surf*boltz_sigma), 0.25_dp)
                  T_surf4 = 3d0/4d0*pow(Teff, 4.d0)*(tau_surf + 2._dp/3._dp)
                  T_surf = pow(T_surf4, 0.25_dp)
                  lnT_surf = log(T_surf)
                  if (.not. skip_partials) then
                     dlnT_dlnR = -0.5_dp
                     dlnT_dL = 1._dp/(4._dp*L_surf)
                  endif
               end if

               if (.not. skip_partials) then
                  dlnT_dlnM = 0._dp; dlnT_dlnkap = 0._dp
                  dlnP_dL = 0._dp; dlnP_dlnR = 0._dp; dlnP_dlnM = 0._dp; dlnP_dlnkap = 0._dp
               endif

            case ('fixed_Psurf_and_Tsurf')

               lnP_surf = safe_log(s% atm_fixed_Psurf)
               T_surf = s% atm_fixed_Tsurf
               lnT_surf = log(T_surf)
               T_surf4 = T_surf*T_surf*T_surf*T_surf
               Teff = pow(4d0/3d0*T_surf4/(tau_surf + 2d0/3d0), 0.25d0)

               if (.not. skip_partials) then
                  dlnT_dL = 0; dlnT_dlnR = 0; dlnT_dlnM = 0; dlnT_dlnkap = 0
                  dlnP_dL = 0; dlnP_dlnR = 0; dlnP_dlnM = 0; dlnP_dlnkap = 0
               endif

            case default

               ! Everything else -- the 'non-trivial atmospheres' ---
               ! gets passed to atm_support
               
               if (.false. .and. 1 == s% solver_test_partials_k .and. &
                     s% solver_iter == s% solver_test_partials_iter_number) then
                  star_debugging_atm_flag = .true.
               end if

               call get_atm_PT( &
                    s, tau_surf, skip_partials, &
                    Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
                    lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
                    ierr)
               if (ierr /= 0) then
                  if (s% report_ierr) write(*,*) 'failed in get_atm_PT'
                  return
               end if

               if (.false. .and. 1 == s% solver_test_partials_k .and. &
                     s% solver_iter == s% solver_test_partials_iter_number) then
                  s% solver_test_partials_val = atm_test_partials_val
                  s% solver_test_partials_dval_dx = atm_test_partials_dval_dx
               end if

            end select

         end if

         ! Add in extra surface pressure

         if (s% surface_extra_Pgas /= 0._dp) then
            P_surf_atm = exp(lnP_surf)
            P_surf = P_surf_atm + s% surface_extra_Pgas
            if (P_surf < 1E-50_dp) then
               lnP_surf = -50*ln10
               if (.not. skip_partials) then
                  dlnP_dL = 0._dp
                  dlnP_dlnR = 0._dp
                  dlnP_dlnM = 0._dp
                  dlnP_dlnkap = 0._dp
               endif
            else
               lnP_surf = log(P_surf)
               if (.not. skip_partials) then
                  dlnP_dL = dlnP_dL*P_surf_atm/P_surf
                  dlnP_dlnR = dlnP_dlnR*P_surf_atm/P_surf
                  dlnP_dlnM = dlnP_dlnM*P_surf_atm/P_surf
                  dlnP_dlnkap = dlnP_dlnkap*P_surf_atm/P_surf
               endif
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

   
      end module hydro_vars

