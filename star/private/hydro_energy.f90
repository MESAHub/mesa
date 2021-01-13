! ***********************************************************************
!
!   Copyright (C) 2018-2019  Bill Paxton & The MESA Team
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


      module hydro_energy

      use star_private_def
      use const_def
      use utils_lib, only: mesa_error, is_bad
      use auto_diff
      use auto_diff_support
      use star_utils, only: em1, e00, ep1

      implicit none

      private
      public :: do1_energy_eqn

      contains
      

      subroutine do1_energy_eqn( & ! energy conservation
            s, k, xscale, equ, skip_partials, do_chem, nvar, ierr)
         use star_utils, only: store_partials
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials, do_chem
         integer, intent(out) :: ierr         
         real(dp), dimension(nvar) :: d_dm1, d_d00, d_dp1      
         include 'formats'         
         call get1_energy_eqn( &
            s, k, xscale, equ, skip_partials, do_chem, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'ierr /= 0 for get1_energy_eqn', k
            return
         end if         
         if (skip_partials) return
         call store_partials(s, k, xscale, s% i_dlnE_dt, nvar, d_dm1, d_d00, d_dp1)
      end subroutine do1_energy_eqn


      subroutine get1_energy_eqn( &
            s, k, xscale, equ, skip_partials, do_chem, nvar, &
            d_dm1, d_d00, d_dp1, ierr)

         use eos_def, only: i_grad_ad, i_lnPgas, i_lnE
         use eps_grav, only: eval_eps_grav_and_partials
         use accurate_sum_auto_diff_18var_order1
         use auto_diff_support
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials, do_chem
         real(dp), intent(out), dimension(nvar) :: d_dm1, d_d00, d_dp1
         integer, intent(out) :: ierr
         
         type(auto_diff_real_18var_order1) :: resid_18, &
            dL_dm_18, sources_18, others_18, dEturb_dt_18, dwork_dm_18, &
            eps_grav_18, dke_dt_18, dpe_dt_18, de_dt_18, P_dV_dt_18
         type(accurate_auto_diff_real_18var_order1) :: esum_18
         real(dp) :: cell_energy_fraction_start, residual, dm, dt, scal
         real(dp), dimension(s% species) :: &
            d_dwork_dxam1, d_dwork_dxa00, d_dwork_dxap1
         integer :: nz, i_dlnE_dt, i_lnd, i_lnT, i_lnR, i_lum, i_v, i_eturb
         logical :: test_partials, doing_op_split_burn, eps_grav_form
                    
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         ierr = 0
         call init
            
         call setup_eps_grav(ierr); if (ierr /= 0) return ! do this first - it sets eps_grav_form         
         call setup_de_dt_and_friends(ierr); if (ierr /= 0) return   
         call setup_dwork_dm(ierr); if (ierr /= 0) return         
         call setup_dL_dm(ierr); if (ierr /= 0) return         
         call setup_sources_and_others(ierr); if (ierr /= 0) return         
         call setup_dEturb_dt(ierr); if (ierr /= 0) return         
         call setup_scal(ierr); if (ierr /= 0) return
         
         s% eps_grav_form_for_energy_eqn(k) = eps_grav_form
         s% dL_dm(k) = dL_dm_18%val
         s% dwork_dm(k) = dwork_dm_18%val
         s% energy_sources(k) = sources_18%val 
            ! nuclear heating, non_nuc_neu_cooling, irradiation heating, extra_heat, eps_mdot
         s% energy_others(k) = others_18%val
            ! eps_WD_sedimentation, eps_diffusion, eps_pre_mix
         s% PdVdt(k) = 0d0
         
         ! sum terms in esum_18 using accurate_auto_diff_real_18var_order1
         if (eps_grav_form) then ! for this case, dwork_dm doesn't include work by P since that is in eps_grav
            esum_18 = - dL_dm_18 + sources_18 + others_18 - dEturb_dt_18 - dwork_dm_18 + eps_grav_18
         else if (s% use_dedt_form_with_total_energy_conservation) then
            esum_18 = - dL_dm_18 + sources_18 + others_18 - dEturb_dt_18 - dwork_dm_18 - dke_dt_18 - dpe_dt_18 - de_dt_18
         else
            esum_18 = - dL_dm_18 + sources_18 + others_18 - dEturb_dt_18 - dwork_dm_18 - P_dV_dt_18 - de_dt_18
            s% PdVdt(k) = P_dV_dt_18%val
         end if
         resid_18 = esum_18 ! convert back to auto_diff_real_18var_order1
         s% ergs_error(k) = -dm*dt*resid_18%val ! save ergs_error before scaling
         resid_18 = scal*resid_18
         residual = resid_18%val
         equ(i_dlnE_dt, k) = residual
         s% E_residual(k) = residual

         if (is_bad(residual)) then
!$omp critical (hydro_equ_l_crit1)
            write(*,2) 'energy eqn residual', k, residual
            stop 'get1_energy_eqn'
!$omp end critical (hydro_equ_l_crit1)
         end if
         if (test_partials) then
            s% solver_test_partials_val = residual
         end if
         if (skip_partials) return
         call unpack_res18(resid_18)

         if (test_partials) then  
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = d_d00(s% solver_test_partials_var)  
            write(*,*) 'get1_energy_eqn', s% solver_test_partials_var
            if (eps_grav_form) write(*,*) 'eps_grav_form', eps_grav_form
            !if (.false. .and. s% solver_iter == s% solver_test_partials_iter_number) then
            if (.true.) then
               write(*,2) 'scal', k, scal
               write(*,2) 'residual', k, residual
               write(*,2) 'sources*scal', k, sources_18%val*scal
               write(*,2) '-dL_dm*scal', k, -dL_dm_18%val*scal
               write(*,2) '-dEturb_dt*scal', k, -dEturb_dt_18%val*scal
               write(*,2) '-dwork_dm*scal', k, -dwork_dm_18%val*scal
               write(*,2) '-dke_dt*scal', k, -dke_dt_18%val*scal
               write(*,2) '-dpe_dt*scal', k, -dpe_dt_18%val*scal
               write(*,2) 'gradT', k, s% gradT(k)
               write(*,2) 'opacity', k, s% opacity(k)
               write(*,2) 'logT', k, s% lnT(k)/ln10
               write(*,2) 'logRho', k, s% lnd(k)/ln10
               write(*,2) 'X', k, s% X(k)
               write(*,2) 'Z', k, s% Z(k)
            end if
            write(*,*)
         end if
         
         contains
         
         subroutine init
            i_dlnE_dt = s% i_dlnE_dt
            i_lnd = s% i_lnd
            i_lnT = s% i_lnT
            i_lnR = s% i_lnR
            i_lum = s% i_lum
            i_v = s% i_v
            i_eturb = s% i_eturb
            nz = s% nz
            dt = s% dt
            dm = s% dm(k)
            cell_energy_fraction_start = &
               s% energy_start(k)*s% dm(k)/s% total_internal_energy_old                    
            doing_op_split_burn = s% op_split_burn .and. &
               s% T_start(k) >= s% op_split_burn_min_T
            d_dm1 = 0d0; d_d00 = 0d0; d_dp1 = 0d0
         end subroutine init
      
         subroutine setup_dwork_dm(ierr)
            integer, intent(out) :: ierr
            real(dp) :: dwork
            logical :: skip_P
            include 'formats'
            ierr = 0
            skip_P = eps_grav_form .or. .not. s% use_dedt_form_with_total_energy_conservation
            ! NOTE: if skip_P then dwork is only that done by turbulence and artificial viscosity
            call eval_dwork(s, k, skip_P, dwork_dm_18, dwork, &
               d_dwork_dxam1, d_dwork_dxa00, d_dwork_dxap1, ierr) 
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'failed in eval_dwork', k
               return
            end if
            dwork_dm_18 = dwork_dm_18/dm
         end subroutine setup_dwork_dm
         
         subroutine setup_dL_dm(ierr)
            integer, intent(out) :: ierr
            type(auto_diff_real_18var_order1) :: L00_18, Lp1_18
            include 'formats'
            ierr = 0         
            L00_18 = wrap_L_00(s, k)
            Lp1_18 = wrap_L_p1(s, k)
            if (s% using_Fraley_time_centering .and. &
                     s% include_L_in_Fraley_time_centering) then
               L00_18 = 0.5d0*(L00_18 + s% L_start(k))
               if (k < s% nz) Lp1_18 = 0.5d0*(Lp1_18 + s% L_start(k+1))
            end if
            dL_dm_18 = (L00_18 - Lp1_18)/dm
         end subroutine setup_dL_dm

         subroutine setup_sources_and_others(ierr) ! sources_18, others_18
            use hydro_Eturb, only: calc_Eq_18
            integer, intent(out) :: ierr
            type(auto_diff_real_18var_order1) :: &
               eps_nuc_18, non_nuc_neu_18, extra_heat_18, Eq_18, RTI_diffusion_18
            include 'formats'
            ierr = 0
         
            if (s% eps_nuc_factor == 0d0 .or. s% nonlocal_NiCo_decay_heat) then
               eps_nuc_18 = 0 ! get eps_nuc from extra_heat instead
            else if (s% op_split_burn .and. s% T_start(k) >= s% op_split_burn_min_T) then
               eps_nuc_18 = 0d0
               eps_nuc_18%val = s% burn_avg_epsnuc(k)
            else
               eps_nuc_18 = 0d0
               eps_nuc_18%val = s% eps_nuc(k)
               eps_nuc_18%d1Array(i_lnd_00) = s% d_epsnuc_dlnd(k)
               eps_nuc_18%d1Array(i_lnT_00) = s% d_epsnuc_dlnT(k)
            end if
            
            non_nuc_neu_18 = 0d0
            ! for reasons lost in the past, we always time center non_nuc_neu
            non_nuc_neu_18%val = 0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))
            non_nuc_neu_18%d1Array(i_lnd_00) = 0.5d0*s% d_nonnucneu_dlnd(k)
            non_nuc_neu_18%d1Array(i_lnT_00) = 0.5d0*s% d_nonnucneu_dlnT(k)
            
            call wrap(extra_heat_18, s% extra_heat(k), &
               s% d_extra_heat_dlndm1(k), s% d_extra_heat_dlnd00(k), s% d_extra_heat_dlndp1(k), &
               s% d_extra_heat_dlnTm1(k), s% d_extra_heat_dlnT00(k), s% d_extra_heat_dlnTp1(k), &
               0d0, 0d0, 0d0, &
               0d0, s% d_extra_heat_dlnR00(k), s% d_extra_heat_dlnRp1(k), &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0)
            
            ! other = eps_WD_sedimentation + eps_diffusion + eps_pre_mix
            ! no partials for any of these
            others_18 = 0d0 
            if (s% do_element_diffusion) then
               if (s% do_WD_sedimentation_heating) then
                  others_18%val = others_18%val + s% eps_WD_sedimentation(k)
               else if (s% do_diffusion_heating) then
                  others_18%val = others_18%val + s% eps_diffusion(k)
               end if
            end if
            if (s% do_conv_premix .and. s% do_premix_heating) &
               others_18%val = others_18%val + s% eps_pre_mix(k)
            
            Eq_18 = 0d0
            if (s% Eturb_flag) then             
               call calc_Eq_18(s, k, Eq_18, ierr)
               if (ierr /= 0) return
            end if   
            
            call setup_RTI_diffusion(RTI_diffusion_18)

            sources_18 = eps_nuc_18 - non_nuc_neu_18 + extra_heat_18 + Eq_18 + RTI_diffusion_18

            sources_18%val = sources_18%val + s% irradiation_heat(k)
            
            if (s% mstar_dot /= 0d0) sources_18%val = sources_18%val + s% eps_mdot(k)

         end subroutine setup_sources_and_others
         
         subroutine setup_RTI_diffusion(diffusion_eps_18)
            type(auto_diff_real_18var_order1), intent(out) :: diffusion_eps_18
            real(dp) :: diffusion_factor, emin_start, sigp1, sig00
            logical :: do_diffusion
            type(auto_diff_real_18var_order1) :: &
               e_m1, e_00, e_p1, diffusion_eps_in, diffusion_eps_out
            include 'formats'
            diffusion_factor = s% dedt_RTI_diffusion_factor
            do_diffusion = s% RTI_flag .and. diffusion_factor > 0d0
            if (.not. do_diffusion) then
               diffusion_eps_18 = 0d0
            else
               if (k < s% nz) then
                  if (s% alpha_RTI(k) > 1d-10 .and. k > 1) then
                     emin_start = min( &
                        s% energy_start(k+1), s% energy_start(k), s% energy_start(k-1))
                     if (emin_start < 5d0*s% RTI_energy_floor) then
                        diffusion_factor = diffusion_factor* &
                           (1d0 + (5d0*s% RTI_energy_floor - emin_start)/emin_start)
                     end if
                  end if
                  sigp1 = diffusion_factor*s% sig_RTI(k+1)
                  e_p1 = wrap_e_p1(s,k)
               else
                  sigp1 = 0
                  e_p1 = 0d0
               end if
               if (k > 1) then
                  sig00 = diffusion_factor*s% sig_RTI(k)
                  e_m1 = wrap_e_m1(s,k)
               else
                  sig00 = 0
                  e_m1 = 0
               end if
               e_00 = wrap_e_00(s,k)
               diffusion_eps_in = sigp1*(e_p1 - e_00)/dm
               diffusion_eps_out = sig00*(e_00 - e_m1)/dm
               diffusion_eps_18 = diffusion_eps_in - diffusion_eps_out
            end if
            s% dedt_RTI(k) = diffusion_eps_18%val
         end subroutine setup_RTI_diffusion
         
         subroutine setup_dEturb_dt(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            dEturb_dt_18 = 0d0
            if (s% Eturb_flag) then
               dEturb_dt_18%val = s% dxh_eturb(k)/dt ! Eturb = Eturb_start + dxh_eturb
               dEturb_dt_18%d1Array(i_eturb_00) = 1d0/dt
            end if
         end subroutine setup_dEturb_dt
         
         subroutine setup_eps_grav(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            
            if (s% u_flag) then ! for now, assume u_flag means no eps_grav 
               eps_grav_form = .false.
               return
            end if

            eps_grav_form = .not. s% use_dedt_form_of_energy_eqn
         
            if (eps_grav_form) then ! check if want it false         
               if (s% always_use_dedt_form_of_energy_eqn) eps_grav_form = .false.            
               if (s% steps_before_always_use_dedt_form_of_energy_eqn >= 0 .and. &
                   s% model_number > s% steps_before_always_use_dedt_form_of_energy_eqn) eps_grav_form = .false.                
            end if
         
            if (.not. eps_grav_form) then ! check if want it true
               if (s% always_use_eps_grav_form_of_energy_eqn) eps_grav_form = .true.             
               if (s% eta_start(k) > s% max_eta_for_dedt_form_of_energy_eqn) eps_grav_form = .true.         
               if (s% doing_relax .and. s% no_dedt_form_during_relax) eps_grav_form = .true.         
               if (s% max_gamma_for_dedt_form_of_energy_eqn > 0d0 .and. .not. eps_grav_form) then
                  if (s% gam_start(k) > s% max_gamma_for_dedt_form_of_energy_eqn) then
                     !write(*,3) 'use eps_grav_form because gamma > max', k, s% model_number, gamma, &
                     !  s% max_gamma_for_dedt_form_of_energy_eqn
                     eps_grav_form = .true.
                  end if
               end if
            end if

            if (eps_grav_form) then
               if (s% Eturb_flag) then
                  stop 'cannot use eps_grav with Eturb yet.  fix energy eqn.'
               end if
               call eval_eps_grav_and_partials(s, k, ierr) ! get eps_grav info
               if (ierr /= 0) then
                  if (s% report_ierr) write(*,2) 'failed in eval_eps_grav_and_partials', k
                  return
               end if
               call wrap(eps_grav_18, s% eps_grav(k), &
                  s% d_eps_grav_dlndm1(k), s% d_eps_grav_dlnd00(k), s% d_eps_grav_dlndp1(k), &
                  s% d_eps_grav_dlnTm1(k), s% d_eps_grav_dlnT00(k), s% d_eps_grav_dlnTp1(k), &
                  0d0, 0d0, 0d0, &
                  0d0, s% d_eps_grav_dlnR00(k), s% d_eps_grav_dlnRp1(k), &
                  0d0, s% d_eps_grav_dv00(k), s% d_eps_grav_dvp1(k), &
                  0d0, s% d_eps_grav_dL00(k), s% d_eps_grav_dLp1(k))
            end if
            
         end subroutine setup_eps_grav

         subroutine setup_de_dt_and_friends(ierr)
            use star_utils, only: get_dke_dt_dpe_dt
            integer, intent(out) :: ierr
            real(dp) :: P_dV
            real(dp) :: dke_dt, d_dkedt_dv00, d_dkedt_dvp1, &
               dpe_dt, d_dpedt_dlnR00, d_dpedt_dlnRp1, &
               de_dt, d_de_dt_dlnd, d_de_dt_dlnT, d_PdV_dlnd, d_PdV_dlnT
            include 'formats'
            ierr = 0

            dke_dt = 0d0; d_dkedt_dv00 = 0d0; d_dkedt_dvp1 = 0d0
            dpe_dt = 0d0; d_dpedt_dlnR00 = 0d0; d_dpedt_dlnRp1 = 0d0
            de_dt = 0d0; d_de_dt_dlnd = 0d0; d_de_dt_dlnT = 0d0
            P_dV = 0d0; d_PdV_dlnd = 0d0; d_PdV_dlnT = 0d0

            if (.not. eps_grav_form) then
               de_dt = (s% energy(k) - s% energy_start(k))/dt
               d_de_dt_dlnd = s% dE_dRho_for_partials(k)*s% rho(k)/dt
               d_de_dt_dlnT = s% Cv_for_partials(k)*s% T(k)/dt
               if (s% use_dedt_form_with_total_energy_conservation) then
                  call get_dke_dt_dpe_dt(s, k, dt, &
                     dke_dt, d_dkedt_dv00, d_dkedt_dvp1, &
                     dpe_dt, d_dpedt_dlnR00, d_dpedt_dlnRp1, ierr)      
                  if (ierr /= 0) then
                     if (s% report_ierr) write(*,2) 'failed in get_dke_dt_dpe_dt', k
                     return
                  end if
                  dke_dt_18 = 0d0
                  dke_dt_18%val = dke_dt
                  dke_dt_18%d1Array(i_v_00) = d_dkedt_dv00
                  dke_dt_18%d1Array(i_v_p1) = d_dkedt_dvp1
                  dpe_dt_18 = 0d0
                  dpe_dt_18%val = dpe_dt
                  dpe_dt_18%d1Array(i_lnR_00) = d_dpedt_dlnR00
                  dpe_dt_18%d1Array(i_lnR_p1) = d_dpedt_dlnRp1
                  de_dt_18 = 0d0
                  de_dt_18%val = de_dt
                  de_dt_18%d1Array(i_lnd_00) = d_de_dt_dlnd
                  de_dt_18%d1Array(i_lnT_00) = d_de_dt_dlnT
               else
                  call get_P_dV(s, k, P_dV, d_PdV_dlnd, d_PdV_dlnT, ierr)
                  if (ierr /= 0) then
                     if (s% report_ierr) write(*,2) 'failed in get_P_dV', k
                     return
                  P_dV_dt_18 = 0d0
                  P_dV_dt_18%val = P_dV/dt
                  P_dV_dt_18%d1Array(i_lnd_00) = d_PdV_dlnd/dt
                  P_dV_dt_18%d1Array(i_lnT_00) = d_PdV_dlnT/dt
                  de_dt_18 = 0d0
                  de_dt_18%val = de_dt
                  de_dt_18%d1Array(i_lnd_00) = d_de_dt_dlnd
                  de_dt_18%d1Array(i_lnT_00) = d_de_dt_dlnT
                  end if
               end if
            end if
            
            s% dkedt(k) = dke_dt
            s% dpedt(k) = dpe_dt
            s% dkedt(k) = dke_dt
            s% dedt(k) = de_dt
         
         end subroutine setup_de_dt_and_friends
         
         subroutine setup_scal(ierr)
            integer, intent(out) :: ierr
            real(dp) :: scal_qp, dt_qp, e0_dq
            include 'formats'
            ierr = 0
            if (k > 1) then
               scal = 1d0
            else
               scal = 1d-6
            end if
            if (s% dedt_eqn_r_scale > 0d0) &
               scal = min(scal, cell_energy_fraction_start*s% dedt_eqn_r_scale) 
            scal_qp = scal              
            dt_qp = dt
            e0_dq = s% energy_start(k)
            scal_qp = scal_qp*dt_qp/e0_dq ! tests show that need qp for this (14700)
            scal = scal_qp
         end subroutine setup_scal
         
         subroutine unpack_res18(res18)
            use star_utils, only: unpack_res18_partials
            type(auto_diff_real_18var_order1) :: res18
            integer :: j

            include 'formats'
            
            call unpack_res18_partials(s, k, nvar, xscale, i_dlnE_dt, &
               res18, d_dm1, d_d00, d_dp1)
            
            ! do partials wrt composition

            if (.not. (s% nonlocal_NiCo_decay_heat .or. doing_op_split_burn)) then
               if (do_chem .and. s% dxdt_nuc_factor > 0d0) then
                  do j=1,s% species
                     call e00(s, xscale, i_dlnE_dt, j+s% nvar_hydro, k, nvar, &
                              scal*s% d_epsnuc_dx(j,k))
                  end do
               end if
            end if

            if (.not. eps_grav_form) then          
               do j=1,s% species
                  call e00(s, xscale, i_dlnE_dt, j+s% nvar_hydro, k, nvar, &
                           -scal*(s%energy(k)/dt)*s% dlnE_dxa_for_partials(j,k))
               end do                           
            else if (do_chem .and. (.not. doing_op_split_burn) .and. &
                     (s% dxdt_nuc_factor > 0d0 .or. s% mix_factor > 0d0)) then                     
               do j=1,s% species
                  call e00(s, xscale, i_dlnE_dt, j+s% nvar_hydro, k, nvar, &
                     scal*s% d_eps_grav_dx(j,k))
               end do               
            end if
            
            do j=1,s% species
               call e00(s, xscale, i_dlnE_dt, j+s% nvar_hydro, k, nvar, -scal*d_dwork_dxa00(j)/dm)
            end do
            if (k > 1) then 
               do j=1,s% species
                  call em1(s, xscale, i_dlnE_dt, j+s% nvar_hydro, k, nvar, -scal*d_dwork_dxam1(j)/dm)
               end do
            end if
            if (k < nz) then
               do j=1,s% species
                  call ep1(s, xscale, i_dlnE_dt, j+s% nvar_hydro, k, nvar, -scal*d_dwork_dxap1(j)/dm)
               end do
            end if         
            
         end subroutine unpack_res18
         
         subroutine unpack1(j, dvar_m1, dvar_00, dvar_p1)
            integer, intent(in) :: j
            real(dp), intent(in) :: dvar_m1, dvar_00, dvar_p1
            d_dm1(j) = dvar_m1
            d_d00(j) = dvar_00
            d_dp1(j) = dvar_p1
         end subroutine unpack1

      end subroutine get1_energy_eqn


      subroutine eval_dwork(s, k, skip_P, dwork_18, dwork, &
            d_dwork_dxam1, d_dwork_dxa00, d_dwork_dxap1, ierr) 
         use accurate_sum_auto_diff_18var_order1
         use auto_diff_support
         use star_utils, only: calc_XP_18_tw
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         logical, intent(in) :: skip_P
         type(auto_diff_real_18var_order1), intent(out) :: dwork_18
         real(dp), intent(out) :: dwork
         real(dp), intent(out), dimension(s% species) :: &
            d_dwork_dxam1, d_dwork_dxa00, d_dwork_dxap1
         integer, intent(out) :: ierr
            
         real(dp) :: work_00, work_p1, dm, dV
         real(dp), dimension(s% species) :: &
            d_work_00_dxa00, d_work_00_dxam1, &
            d_work_p1_dxap1, d_work_p1_dxa00, d_XP_dxa
         type(auto_diff_real_18var_order1) :: work_00_18, work_p1_18, &
            XP_18, dV_18, rho_18
         logical :: test_partials
         integer :: j
         include 'formats'
         ierr = 0

         call eval1_work(s, k, skip_P, &
            work_00_18, work_00, d_work_00_dxa00, d_work_00_dxam1, ierr)
         if (ierr /= 0) return
         call eval1_work(s, k+1, skip_P, &
            work_p1_18, work_p1, d_work_p1_dxap1, d_work_p1_dxa00, ierr)
         if (ierr /= 0) return
         work_p1_18 = shift_p1(work_p1_18) ! shift the partials         
         dwork_18 = work_00_18 - work_p1_18
         dwork = dwork_18%val
         do j=1,s% species
            d_dwork_dxam1(j) = d_work_00_dxam1(j)
            d_dwork_dxa00(j) = d_work_00_dxa00(j) - d_work_p1_dxa00(j)
            d_dwork_dxap1(j) = -d_work_p1_dxap1(j)
         end do         

         !test_partials = (k == s% solver_test_partials_k) 
         test_partials = .false.
            
         if (test_partials) then
            s% solver_test_partials_val = work_00
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = 0d0
            write(*,*) 'eval_dwork', s% solver_test_partials_var
         end if

      end subroutine eval_dwork
      

      ! ergs/s at face(k)
      subroutine eval1_work(s, k, skip_P, &
            work_18, work, d_work_dxa00, d_work_dxam1, ierr)
         use star_utils, only: get_avQ_18, calc_Pt_18_tw
         use accurate_sum_auto_diff_18var_order1
         use auto_diff_support
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         logical, intent(in) :: skip_P
         type(auto_diff_real_18var_order1), intent(out) :: work_18
         real(dp), intent(out) :: work
         real(dp), dimension(s% species), intent(out) :: &
            d_work_dxa00, d_work_dxam1
         integer, intent(out) :: ierr
         real(dp) :: alfa, beta, A, P_face, u_face, theta, P_theta
         real(dp), dimension(s% species) :: d_Pface_dxa00, d_Pface_dxam1
         type(auto_diff_real_18var_order1) :: &
            A_18, P_face_18, u_face_18, mlt_Pturb_18, &
            PtR_18, PtL_18, avQL_18, avQR_18, PL_18, PR_18
         logical :: test_partials
         integer :: j
         include 'formats'
         ierr = 0

         if (k > s% nz .or. (s% dt <= 0d0 .and. .not. (s% v_flag .or. s% u_flag))) then
            work_18 = 0d0
            if (k == s% nz+1) then
               work = 4*pi*s% r_center*s% r_center*s% P_start(s% nz)*s% v_center
               s% work_inward_at_center = work
            end if
            work_18%val = work
            d_work_dxa00 = 0d0
            d_work_dxam1 = 0d0
            return    
         end if
         
         A_18 = 0d0
         A_18%val = 4d0*pi*s% R2(k)
         A_18%d1Array(i_lnR_00) = 4d0*pi*s% d_R2_dlnR(k)
         A = A_18%val

         if (s% u_flag) then ! keep it simple for now.
         
            call wrap_P_face(s, P_face_18, k) 
            call wrap_u_face(s, u_face_18, k) 
            d_Pface_dxa00 = 0d0
            d_Pface_dxam1 = 0d0
         
         else

            theta = 1d0
            P_theta = 1d0
            if (s% using_Fraley_time_centering) then
               theta = 0.5d0
               if (s% include_P_in_Fraley_time_centering) P_theta = 0.5d0
            end if

            PR_18 = 0d0
            avQR_18 = 0d0
            PtR_18 = 0d0         
            if (k > 1) then 
               if (s% use_avQ_art_visc) then
                  call get_avQ_18(s, k-1, avQR_18, ierr)
                  if (ierr /= 0) return
                  avQR_18 = shift_m1(avQR_18)
                  avQR_18 = 0.5d0*(avQR_18 + s% avQ_start(k-1))
               end if            
               if (s% Eturb_flag) then
                  call calc_Pt_18_tw(s, k-1, PtR_18, ierr)
                  if (ierr /= 0) return
                  PtR_18 = shift_m1(PtR_18)
               end if            
               if (.not. skip_P) then
                  PR_18 = wrap_p_m1(s, k)
                  if (s% using_Fraley_time_centering .and. &
                           s% include_P_in_Fraley_time_centering) &
                     PR_18 = 0.5d0*(PR_18 + s% P_start(k))
               end if            
               if (s% use_other_pressure) PR_18%val = PR_18%val + s% extra_pressure(k-1)            
            end if
      
            avQL_18 = 0d0
            if (s% use_avQ_art_visc) then
               call get_avQ_18(s, k, avQL_18, ierr)
               if (ierr /= 0) return
               avQL_18 = 0.5d0*(avQL_18 + s% avQ_start(k))
            end if
         
            PtL_18 = 0d0
            if (s% Eturb_flag) then
               call calc_Pt_18_tw(s, k, PtL_18, ierr)
               if (ierr /= 0) return
            end if
         
            if (skip_P) then
               PL_18 = 0d0
            else
               PL_18 = wrap_p_00(s, k)
               if (s% using_Fraley_time_centering .and. &
                        s% include_P_in_Fraley_time_centering) &
                  PL_18 = 0.5d0*(PL_18 + s% P_start(k))
            end if
            if (s% use_other_pressure) PL_18%val = PL_18%val + s% extra_pressure(k)
         
            d_Pface_dxa00 = 0d0
            d_Pface_dxam1 = 0d0
         
            if (k > 1) then         
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
               beta = 1d0 - alfa
               P_face_18 = alfa*(PL_18 + avQL_18 + PtL_18) + beta*(PR_18 + avQR_18 + PtR_18)            
               if (.not. skip_P) then
                  do j=1,s% species
                     d_Pface_dxa00(j) = d_Pface_dxa00(j) + &
                        alfa*s% dlnP_dxa_for_partials(j,k)*P_theta*s% P(k)
                  end do
                  do j=1,s% species
                     d_Pface_dxam1(j) = d_Pface_dxam1(j) + &
                        beta*s% dlnP_dxa_for_partials(j,k-1)*P_theta*s% P(k-1)
                  end do
               end if         
               if (s% mlt_Pturb_factor > 0d0 .and. s% mlt_vc_start(k) > 0d0 .and. k > 1) then
                  mlt_Pturb_18 = 0d0
                  mlt_Pturb_18%val = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*(s% rho(k-1) + s% rho(k))/6d0
                  mlt_Pturb_18%d1Array(i_lnd_m1) = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*s% rho(k-1)/6d0
                  mlt_Pturb_18%d1Array(i_lnd_00) = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*s% rho(k)/6d0
                  P_face_18 = P_face_18 + mlt_Pturb_18
               end if            
            else ! k == 1
               P_face_18 = PL_18 + avQL_18 + PtL_18            
               if (.not. skip_P) then
                  do j=1,s% species
                     d_Pface_dxa00(j) = d_Pface_dxa00(j) + &
                        s% dlnP_dxa_for_partials(j,k)*P_theta*s% P(k)
                  end do
               end if            
            end if

            u_face_18 = 0d0
            if (s% v_flag) then
               u_face_18%val = s% vc(k)
               u_face_18%d1Array(i_v_00) = s% d_vc_dv
            else
               u_face_18%val = theta*(s% r(k) - s% r_start(k))/s% dt
               u_face_18%d1Array(i_lnR_00) = theta*s% r(k)/s% dt
            end if
            
         end if
         
         work_18 = A_18*P_face_18*u_face_18
         P_face = P_face_18%val
         u_face = u_face_18%val
         work = work_18%val
         if (k == 1) s% work_outward_at_surface = work
         
         do j=1,s% species
            d_work_dxa00(j) = A*d_Pface_dxa00(j)*u_face
            d_work_dxam1(j) = A*d_Pface_dxam1(j)*u_face
         end do

         if (is_bad(work)) then
!$omp critical (hydro_equ_l_crit2)
            write(*,2) 'work', k, work
            stop 'eval1_work'
!$omp end critical (hydro_equ_l_crit2)
         end if

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
            
         if (test_partials) then
            s% solver_test_partials_val = u_face
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = 0d0
            write(*,*) 'eval1_work', s% solver_test_partials_var
         end if
         
      end subroutine eval1_work


      subroutine get_P_dV(s, k, P_dV, d_PdV_dlnd, d_PdV_dlnT, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: P_dV
         real(dp), intent(out) :: d_PdV_dlnd, d_PdV_dlnT
         integer, intent(out) :: ierr
         real(dp) :: rho, dlnd, theta, P, d_P_dlnd00, d_P_dlnT00, dV, d_dV_dlnd
         logical :: test_partials
         include 'formats'         
         ierr = 0         
         
         if (s% using_Fraley_time_centering .and. &
               s% include_P_in_Fraley_time_centering) then
            theta = 0.5d0
            P = 0.5d0*(s% P(k) + s% P_start(k))
         else
            theta = 1d0
            P = s% P(k)
         end if
         d_P_dlnd00 = theta*s% P(k)*s% chiRho_for_partials(k)
         d_P_dlnT00 = theta*s% P(k)*s% chiT_for_partials(k)

         rho = s% rho(k)
         dlnd = s% dxh_lnd(k) ! solver value used for lnd = lnd_start + dxh_lnd

         ! dV = 1/rho - 1/rho_start = 
         ! -(rho/rho_start - 1)/rho = 
         ! -(exp(lnd)/exp(lnd_start) - 1)/rho = 
         ! -(exp(lnd - lnd_start) - 1)/rho =
         ! -(exp(dlnd) - 1)/rho = 
         ! -expm1(dlnd)/rho
         
         dV = -expm1(dlnd)/rho
         d_dV_dlnd = -1d0/rho ! dV = 1/rho - 1/rho_start         
         
         P_dV = P*dV
         d_PdV_dlnd = d_P_dlnd00*dV + P*d_dV_dlnd
         d_PdV_dlnT = d_P_dlnT00*dV

         !test_partials = (k == s% solver_test_partials_k) 
         test_partials = .false.            
         if (test_partials) then
            s% solver_test_partials_val = P_dV
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = d_PdV_dlnd
            write(*,*) 'get_P_dV', s% solver_test_partials_var
         end if

      end subroutine get_P_dV

      
      end module hydro_energy

