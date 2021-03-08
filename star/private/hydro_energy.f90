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
      use star_utils, only: em1, e00, ep1, set_energy_eqn_scal

      implicit none

      private
      public :: do1_energy_eqn

      contains
      

      subroutine do1_energy_eqn( & ! energy conservation
            s, k, skip_partials, do_chem, nvar, ierr)
         use star_utils, only: store_partials
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials, do_chem
         integer, intent(out) :: ierr         
         real(dp), dimension(nvar) :: d_dm1, d_d00, d_dp1      
         include 'formats'
         call get1_energy_eqn( &
            s, k, skip_partials, do_chem, nvar, &
            d_dm1, d_d00, d_dp1, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'ierr /= 0 for get1_energy_eqn', k
            return
         end if         
         if (skip_partials) return
         call store_partials( &
            s, k, s% i_dlnE_dt, nvar, d_dm1, d_d00, d_dp1, 'do1_energy_eqn', ierr)
      end subroutine do1_energy_eqn


      subroutine get1_energy_eqn( &
            s, k, skip_partials, do_chem, nvar, &
            d_dm1, d_d00, d_dp1, ierr)

         use eos_def, only: i_grad_ad, i_lnPgas, i_lnE
         use eps_grav, only: eval_eps_grav_and_partials
         use accurate_sum_auto_diff_star_order1
         use auto_diff_support
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials, do_chem
         real(dp), intent(out), dimension(nvar) :: d_dm1, d_d00, d_dp1
         integer, intent(out) :: ierr
         
         type(auto_diff_real_star_order1) :: resid_ad, &
            dL_dm_ad, sources_ad, others_ad, d_turbulent_energy_dt_ad, &
            dwork_dm_ad, eps_grav_ad, dke_dt_ad, dpe_dt_ad, de_dt_ad
         type(accurate_auto_diff_real_star_order1) :: esum_ad
         real(dp) :: residual, dm, dt, scal
         real(dp), dimension(s% species) :: &
            d_dwork_dxam1, d_dwork_dxa00, d_dwork_dxap1
         integer :: nz, i_dlnE_dt, i_lum, i_v
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
         call setup_d_turbulent_energy_dt(ierr); if (ierr /= 0) return         
         call set_energy_eqn_scal(s, k, scal, ierr); if (ierr /= 0) return
         
         s% eps_grav_form_for_energy_eqn(k) = eps_grav_form
         s% dL_dm(k) = dL_dm_ad%val
         s% dwork_dm(k) = dwork_dm_ad%val
         s% energy_sources(k) = sources_ad%val 
            ! nuclear heating, non_nuc_neu_cooling, irradiation heating, extra_heat, eps_mdot
         s% energy_others(k) = others_ad%val
            ! eps_WD_sedimentation, eps_diffusion, eps_pre_mix
         ! sum terms in esum_ad using accurate_auto_diff_real_star_order1
         if (eps_grav_form) then ! for this case, dwork_dm doesn't include work by P since that is in eps_grav
            esum_ad = - dL_dm_ad + sources_ad + &
               others_ad - d_turbulent_energy_dt_ad - dwork_dm_ad + eps_grav_ad
         else if (s% using_velocity_time_centering .and. &
                s% use_P_d_1_div_rho_form_of_work_when_time_centering_velocity) then
            esum_ad = - dL_dm_ad + sources_ad + &
               others_ad - d_turbulent_energy_dt_ad - dwork_dm_ad - de_dt_ad
         else
            esum_ad = - dL_dm_ad + sources_ad + &
               others_ad - d_turbulent_energy_dt_ad - dwork_dm_ad - dke_dt_ad - dpe_dt_ad - de_dt_ad
         end if
         resid_ad = esum_ad ! convert back to auto_diff_real_star_order1
         s% ergs_error(k) = -dm*dt*resid_ad%val ! save ergs_error before scaling
         resid_ad = scal*resid_ad
         residual = resid_ad%val
         s% equ(i_dlnE_dt, k) = residual
         
         if (test_partials) then
            s% solver_test_partials_val = residual
         end if
         if (skip_partials) return
         call unpack_res18(s% species, resid_ad)

         if (test_partials) then  
            s% solver_test_partials_var = s% i_u
            s% solver_test_partials_dval_dx = d_d00(s% solver_test_partials_var)  
            write(*,*) 'get1_energy_eqn', s% solver_test_partials_var
            if (eps_grav_form) write(*,*) 'eps_grav_form', eps_grav_form
            !if (.false. .and. s% solver_iter == s% solver_test_partials_iter_number) then
            if (.true.) then
               write(*,2) 'scal', k, scal
               write(*,2) 'residual', k, residual
               write(*,2) 'sources*scal', k, sources_ad%val*scal
               write(*,2) '-dL_dm*scal', k, -dL_dm_ad%val*scal
               write(*,2) '-d_turbulent_energy_dt*scal', k, -d_turbulent_energy_dt_ad%val*scal
               write(*,2) '-dwork_dm*scal', k, -dwork_dm_ad%val*scal
               write(*,2) '-dke_dt*scal', k, -dke_dt_ad%val*scal
               write(*,2) '-dpe_dt*scal', k, -dpe_dt_ad%val*scal
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
            i_lum = s% i_lum
            i_v = s% i_v
            nz = s% nz
            dt = s% dt
            dm = s% dm(k)
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
            skip_P = eps_grav_form
            if (s% using_velocity_time_centering .and. &
                s% use_P_d_1_div_rho_form_of_work_when_time_centering_velocity) then
               call eval_Fraley_PdV_work(s, k, skip_P, dwork_dm_ad, dwork, &
                  d_dwork_dxa00, ierr) 
               d_dwork_dxam1 = 0
               d_dwork_dxap1 = 0
            else
               call eval_dwork(s, k, skip_P, dwork_dm_ad, dwork, &
                  d_dwork_dxam1, d_dwork_dxa00, d_dwork_dxap1, ierr) 
            end if
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'failed in setup_dwork_dm', k
               return
            end if
            dwork_dm_ad = dwork_dm_ad/dm
         end subroutine setup_dwork_dm
         
         subroutine setup_dL_dm(ierr)  
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: &
               L00_ad, Lp1_ad, unused
            include 'formats'
            ierr = 0         
            L00_ad = wrap_L_00(s, k)
            Lp1_ad = wrap_L_p1(s, k)
            if (s% using_velocity_time_centering .and. &
                     s% include_L_in_velocity_time_centering) then
               L00_ad = 0.5d0*(L00_ad + s% L_start(k))
               if (k < s% nz) Lp1_ad = 0.5d0*(Lp1_ad + s% L_start(k+1))
            end if
            dL_dm_ad = (L00_ad - Lp1_ad)/dm
         end subroutine setup_dL_dm

         subroutine setup_sources_and_others(ierr) ! sources_ad, others_ad
            use hydro_tdc, only: compute_Eq_cell
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: &
               eps_nuc_ad, non_nuc_neu_ad, extra_heat_ad, Eq_ad, RTI_diffusion_ad
            real(dp) :: d_extra_heat_dlnR00, d_extra_heat_dlnRp1, &
               d_extra_heat_dlnTm1, d_extra_heat_dlnT00, d_extra_heat_dlnTp1, &
               d_extra_heat_dlndm1, d_extra_heat_dlnd00, d_extra_heat_dlndp1
            include 'formats'
            ierr = 0
         
            if (s% eps_nuc_factor == 0d0 .or. s% nonlocal_NiCo_decay_heat) then
               eps_nuc_ad = 0 ! get eps_nuc from extra_heat instead
            else if (s% op_split_burn .and. s% T_start(k) >= s% op_split_burn_min_T) then
               eps_nuc_ad = 0d0
               eps_nuc_ad%val = s% burn_avg_epsnuc(k)
            else
               eps_nuc_ad = 0d0
               eps_nuc_ad%val = s% eps_nuc(k)
               eps_nuc_ad%d1Array(i_lnd_00) = s% d_epsnuc_dlnd(k)
               eps_nuc_ad%d1Array(i_lnT_00) = s% d_epsnuc_dlnT(k)
            end if
            
            non_nuc_neu_ad = 0d0
            ! for reasons lost in the past, we always time center non_nuc_neu
            non_nuc_neu_ad%val = 0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))
            non_nuc_neu_ad%d1Array(i_lnd_00) = 0.5d0*s% d_nonnucneu_dlnd(k)
            non_nuc_neu_ad%d1Array(i_lnT_00) = 0.5d0*s% d_nonnucneu_dlnT(k)
            
            d_extra_heat_dlnR00 = s% d_extra_heat_dlnR00(k)
            d_extra_heat_dlnRp1 = s% d_extra_heat_dlnRp1(k)
            d_extra_heat_dlnTm1 = s% d_extra_heat_dlnTm1(k)
            d_extra_heat_dlnT00 = s% d_extra_heat_dlnT00(k)
            d_extra_heat_dlnTp1 = s% d_extra_heat_dlnTp1(k)
            d_extra_heat_dlndm1 = s% d_extra_heat_dlndm1(k)
            d_extra_heat_dlnd00 = s% d_extra_heat_dlnd00(k)
            d_extra_heat_dlndp1 = s% d_extra_heat_dlndp1(k)
            call wrap(extra_heat_ad, s% extra_heat(k), &
               d_extra_heat_dlndm1, d_extra_heat_dlnd00, d_extra_heat_dlndp1, &
               d_extra_heat_dlnTm1, d_extra_heat_dlnT00, d_extra_heat_dlnTp1, &
               0d0, 0d0, 0d0, &
               0d0, d_extra_heat_dlnR00, d_extra_heat_dlnRp1, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0)
            
            ! other = eps_WD_sedimentation + eps_diffusion + eps_pre_mix
            ! no partials for any of these
            others_ad = 0d0 
            if (s% do_element_diffusion) then
               if (s% do_WD_sedimentation_heating) then
                  others_ad%val = others_ad%val + s% eps_WD_sedimentation(k)
               else if (s% do_diffusion_heating) then
                  others_ad%val = others_ad%val + s% eps_diffusion(k)
               end if
            end if
            if (s% do_conv_premix .and. s% do_premix_heating) &
               others_ad%val = others_ad%val + s% eps_pre_mix(k)
            
            Eq_ad = 0d0
            if (s% TDC_flag) then             
               Eq_ad = compute_Eq_cell(s, k, ierr)
               if (ierr /= 0) return
            end if   
            
            call setup_RTI_diffusion(RTI_diffusion_ad)

            sources_ad = eps_nuc_ad - non_nuc_neu_ad + extra_heat_ad + Eq_ad + RTI_diffusion_ad

            sources_ad%val = sources_ad%val + s% irradiation_heat(k)
            
            if (s% mstar_dot /= 0d0) sources_ad%val = sources_ad%val + s% eps_mdot(k)

         end subroutine setup_sources_and_others
         
         subroutine setup_RTI_diffusion(diffusion_eps_ad)
            type(auto_diff_real_star_order1), intent(out) :: diffusion_eps_ad
            real(dp) :: diffusion_factor, emin_start, sigp1, sig00
            logical :: do_diffusion
            type(auto_diff_real_star_order1) :: &
               e_m1, e_00, e_p1, diffusion_eps_in, diffusion_eps_out
            include 'formats'
            diffusion_factor = s% dedt_RTI_diffusion_factor
            do_diffusion = s% RTI_flag .and. diffusion_factor > 0d0
            if (.not. do_diffusion) then
               diffusion_eps_ad = 0d0
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
               diffusion_eps_ad = diffusion_eps_in - diffusion_eps_out
            end if
            s% dedt_RTI(k) = diffusion_eps_ad%val
         end subroutine setup_RTI_diffusion
         
         subroutine setup_d_turbulent_energy_dt(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            if (s% TDC_flag) then
               d_turbulent_energy_dt_ad = wrap_dxh_etrb(s,k)/dt
            else
               d_turbulent_energy_dt_ad = 0d0
            end if
            s% detrbdt(k) = d_turbulent_energy_dt_ad%val
         end subroutine setup_d_turbulent_energy_dt
         
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
               if (s% TDC_flag) then
                  stop 'cannot use eps_grav with et yet.  fix energy eqn.'
               end if
               call eval_eps_grav_and_partials(s, k, ierr) ! get eps_grav info
               if (ierr /= 0) then
                  if (s% report_ierr) write(*,2) 'failed in eval_eps_grav_and_partials', k
                  return
               end if
               eps_grav_ad = s% eps_grav_ad(k)
            end if
            
         end subroutine setup_eps_grav

         subroutine setup_de_dt_and_friends(ierr)
            use star_utils, only: get_dke_dt_dpe_dt
            integer, intent(out) :: ierr
            real(dp) :: dke_dt, d_dkedt_dv00, d_dkedt_dvp1, &
               dpe_dt, d_dpedt_dlnR00, d_dpedt_dlnRp1, &
               de_dt, d_de_dt_dlnd, d_de_dt_dlnT
            include 'formats'
            ierr = 0

            dke_dt = 0d0; d_dkedt_dv00 = 0d0; d_dkedt_dvp1 = 0d0
            dpe_dt = 0d0; d_dpedt_dlnR00 = 0d0; d_dpedt_dlnRp1 = 0d0
            de_dt = 0d0; d_de_dt_dlnd = 0d0; d_de_dt_dlnT = 0d0

            if (.not. eps_grav_form) then
            
               de_dt = (s% energy(k) - s% energy_start(k))/dt
               d_de_dt_dlnd = s% dE_dRho_for_partials(k)*s% rho(k)/dt
               d_de_dt_dlnT = s% Cv_for_partials(k)*s% T(k)/dt
               de_dt_ad = 0d0
               de_dt_ad%val = de_dt
               de_dt_ad%d1Array(i_lnd_00) = d_de_dt_dlnd
               de_dt_ad%d1Array(i_lnT_00) = d_de_dt_dlnT
               
               call get_dke_dt_dpe_dt(s, k, dt, &
                  dke_dt, d_dkedt_dv00, d_dkedt_dvp1, &
                  dpe_dt, d_dpedt_dlnR00, d_dpedt_dlnRp1, ierr)      
               if (ierr /= 0) then
                  if (s% report_ierr) write(*,2) 'failed in get_dke_dt_dpe_dt', k
                  return
               end if
               dke_dt_ad = 0d0
               dke_dt_ad%val = dke_dt
               dke_dt_ad%d1Array(i_v_00) = d_dkedt_dv00
               dke_dt_ad%d1Array(i_v_p1) = d_dkedt_dvp1
               
               dpe_dt_ad = 0d0
               dpe_dt_ad%val = dpe_dt
               dpe_dt_ad%d1Array(i_lnR_00) = d_dpedt_dlnR00
               dpe_dt_ad%d1Array(i_lnR_p1) = d_dpedt_dlnRp1
               
            end if
            
            s% dkedt(k) = dke_dt
            s% dpedt(k) = dpe_dt
            s% dkedt(k) = dke_dt
            s% dedt(k) = de_dt
         
         end subroutine setup_de_dt_and_friends
         
         subroutine unpack_res18(species,res18)
            use star_utils, only: save_eqn_dxa_partials, unpack_residual_partials
            type(auto_diff_real_star_order1) :: res18
            integer, intent(in) :: species
            real(dp) :: dequ
            integer :: j
            real(dp), dimension(species) :: dxam1, dxa00, dxap1
            logical, parameter :: checking = .true.
            include 'formats'
            
            ! do partials wrt composition
            dxam1 = 0d0; dxa00 = 0d0; dxap1 = 0d0
            if (.not. (s% nonlocal_NiCo_decay_heat .or. doing_op_split_burn)) then
               if (do_chem .and. s% dxdt_nuc_factor > 0d0) then
                  do j=1,s% species
                     dequ = scal*s% d_epsnuc_dx(j,k)
                     if (checking) call check_dequ(dequ,'d_epsnuc_dx')
                     dxa00(j) = dxa00(j) + dequ
                  end do
               end if
            end if

            if (.not. eps_grav_form) then          
               do j=1,s% species
                  dequ = -scal*(s%energy(k)/dt)*s% dlnE_dxa_for_partials(j,k)
                  if (checking) call check_dequ(dequ,'dlnE_dxa_for_partials')
                  dxa00(j) = dxa00(j) + dequ
               end do                           
            else if (do_chem .and. (.not. doing_op_split_burn) .and. &
                     (s% dxdt_nuc_factor > 0d0 .or. s% mix_factor > 0d0)) then                     
               do j=1,s% species
                  dequ = scal*s% d_eps_grav_dx(j,k)
                  if (checking) call check_dequ(dequ,'d_eps_grav_dx')
                  dxa00(j) = dxa00(j) + dequ
               end do               
            end if
            
            do j=1,s% species
               dequ = -scal*d_dwork_dxa00(j)/dm
               if (checking) call check_dequ(dequ,'d_dwork_dxa00')
               dxa00(j) = dxa00(j) + dequ
            end do
            if (k > 1) then 
               do j=1,s% species
                  dequ = -scal*d_dwork_dxam1(j)/dm
                  if (checking) call check_dequ(dequ,'d_dwork_dxam1')
                  dxam1(j) = dxam1(j) + dequ
               end do
            end if
            if (k < nz) then
               do j=1,s% species
                  dequ = -scal*d_dwork_dxap1(j)/dm
                  if (checking) call check_dequ(dequ,'d_dwork_dxap1')
                  dxap1(j) = dxap1(j) + dequ
               end do
            end if         

            call save_eqn_dxa_partials(&
               s, k, nvar, i_dlnE_dt, species, dxam1, dxa00, dxap1, 'get1_energy_eqn', ierr)
            
            call unpack_residual_partials(s, k, nvar, i_dlnE_dt, &
               res18, d_dm1, d_d00, d_dp1)
            
         end subroutine unpack_res18

         subroutine check_dequ(dequ, str)
            real(dp), intent(in) :: dequ
            character (len=*), intent(in) :: str
            include 'formats'
            if (is_bad(dequ)) then
!$omp critical (hydro_energy_crit2)
               ierr = -1
               if (s% report_ierr) then
                  write(*,2) 'get1_energy_eqn: bad ' // trim(str), k, dequ
               end if
               if (s% stop_for_bad_nums) stop 'get1_energy_eqn'
!$omp end critical (hydro_energy_crit2)
               return
            end if
         end subroutine check_dequ
         
         subroutine unpack1(j, dvar_m1, dvar_00, dvar_p1)
            integer, intent(in) :: j
            real(dp), intent(in) :: dvar_m1, dvar_00, dvar_p1
            d_dm1(j) = dvar_m1
            d_d00(j) = dvar_00
            d_dp1(j) = dvar_p1
         end subroutine unpack1

      end subroutine get1_energy_eqn


      subroutine eval_dwork(s, k, skip_P, dwork_ad, dwork, &
            d_dwork_dxam1, d_dwork_dxa00, d_dwork_dxap1, ierr) 
         use accurate_sum_auto_diff_star_order1
         use auto_diff_support
         use star_utils, only: calc_Ptot_ad_tw
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         logical, intent(in) :: skip_P
         type(auto_diff_real_star_order1), intent(out) :: dwork_ad
         real(dp), intent(out) :: dwork
         real(dp), intent(out), dimension(s% species) :: &
            d_dwork_dxam1, d_dwork_dxa00, d_dwork_dxap1
         integer, intent(out) :: ierr
            
         real(dp) :: work_00, work_p1, dm, dV
         real(dp), dimension(s% species) :: &
            d_work_00_dxa00, d_work_00_dxam1, &
            d_work_p1_dxap1, d_work_p1_dxa00, d_Ptot_dxa
         type(auto_diff_real_star_order1) :: work_00_ad, work_p1_ad, &
            Ptot_ad, dV_ad, rho_ad
         logical :: test_partials
         integer :: j
         include 'formats'
         ierr = 0

         call eval1_work(s, k, skip_P, &
            work_00_ad, work_00, d_work_00_dxa00, d_work_00_dxam1, ierr)
         if (ierr /= 0) return
         call eval1_work(s, k+1, skip_P, &
            work_p1_ad, work_p1, d_work_p1_dxap1, d_work_p1_dxa00, ierr)
         if (ierr /= 0) return
         work_p1_ad = shift_p1(work_p1_ad) ! shift the partials         
         dwork_ad = work_00_ad - work_p1_ad
         dwork = dwork_ad%val
         do j=1,s% species
            d_dwork_dxam1(j) = d_work_00_dxam1(j)
            d_dwork_dxa00(j) = d_work_00_dxa00(j) - d_work_p1_dxa00(j)
            d_dwork_dxap1(j) = -d_work_p1_dxap1(j)
         end do         

         !test_partials = (k == s% solver_test_partials_k) 
         test_partials = .false.
            
         if (test_partials) then
            s% solver_test_partials_val = 0
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'eval_dwork', s% solver_test_partials_var
         end if

      end subroutine eval_dwork
      

      ! ergs/s at face(k)
      subroutine eval1_work(s, k, skip_Peos, &
            work_ad, work, d_work_dxa00, d_work_dxam1, ierr)
         use star_utils, only: get_Pvsc_ad, calc_Ptrb_ad_tw
         use accurate_sum_auto_diff_star_order1
         use auto_diff_support
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         logical, intent(in) :: skip_Peos
         type(auto_diff_real_star_order1), intent(out) :: work_ad
         real(dp), intent(out) :: work
         real(dp), dimension(s% species), intent(out) :: &
            d_work_dxa00, d_work_dxam1
         integer, intent(out) :: ierr
         real(dp) :: alfa, beta, P_theta, extra_P, Peos_face, Av_face
         real(dp), dimension(s% species) :: d_Pface_dxa00, d_Pface_dxam1
         type(auto_diff_real_star_order1) :: &
            P_face_ad, A_times_v_face_ad, mlt_Pturb_ad, &
            PtrbR_ad, PtrbL_ad, PvscL_ad, PvscR_ad, PL_ad, PR_ad, &
            Peos_ad, Ptrb_ad, Pvsc_ad, inv_R2
         logical :: test_partials
         integer :: j
         include 'formats'
         ierr = 0

         d_work_dxa00 = 0d0
         d_work_dxam1 = 0d0
         if (k > s% nz .or. (s% dt <= 0d0 .and. .not. (s% v_flag .or. s% u_flag))) then
            work_ad = 0d0
            if (k == s% nz+1) then
               work = pi4*s% r_center*s% r_center*s% Peos_start(s% nz)*s% v_center
               s% work_inward_at_center = work
            end if
            work_ad%val = work
            return    
         end if
         
         call eval1_A_times_v_face_ad(s, k, A_times_v_face_ad, ierr)
         if (ierr /= 0) return

         if (k > 1) then         
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         else
            alfa = 1d0
         end if
         beta = 1d0 - alfa

         if (s% u_flag) then
            P_face_ad = s% P_face_ad(k)
            if (s% using_velocity_time_centering .and. &
                     s% include_P_in_velocity_time_centering) &
               P_face_ad = 0.5d0*(P_face_ad + s% P_face_start(k))
            d_Pface_dxa00 = 0d0
            d_Pface_dxam1 = 0d0
         else ! set P_ad
            d_Pface_dxa00 = 0d0
            d_Pface_dxam1 = 0d0
            if (skip_Peos) then
               Peos_ad = 0d0
            else
               if (k > 1) then 
                  PR_ad = wrap_Peos_m1(s,k)
                  if (s% using_velocity_time_centering .and. &
                           s% include_P_in_velocity_time_centering) &
                     PR_ad = 0.5d0*(PR_ad + s% Peos_start(k-1))
               else
                  PR_ad = 0d0
               end if
               PL_ad = wrap_Peos_00(s,k)
               if (s% using_velocity_time_centering .and. &
                        s% include_P_in_velocity_time_centering) &
                  PL_ad = 0.5d0*(PL_ad + s% Peos_start(k))
               Peos_ad = alfa*PL_ad + beta*PR_ad
               if (s% using_velocity_time_centering .and. &
                        s% include_P_in_velocity_time_centering) then
                  P_theta = 0.5d0
               else
                  P_theta = 1d0
               end if
               if (k > 1) then         
                  do j=1,s% species
                     d_Pface_dxa00(j) = &
                        alfa*s% dlnPeos_dxa_for_partials(j,k)*P_theta*s% Peos(k)
                  end do
                  do j=1,s% species
                     d_Pface_dxam1(j) = &
                        beta*s% dlnPeos_dxa_for_partials(j,k-1)*P_theta*s% Peos(k-1)
                  end do
               else ! k == 1
                  do j=1,s% species
                     d_Pface_dxa00(j) = &
                        s% dlnPeos_dxa_for_partials(j,k)*P_theta*s% Peos(k)
                  end do
               end if
            end if
         
            ! set Pvsc_ad
            if (.not. s% use_Pvsc_art_visc) then
               Pvsc_ad = 0d0
            else
               if (k > 1) then 
                  call get_Pvsc_ad(s, k-1, PvscR_ad, ierr)
                  if (ierr /= 0) return
                  PvscR_ad = shift_m1(PvscR_ad)
                  ! always time center
                  PvscR_ad = 0.5d0*(PvscR_ad + s% Pvsc_start(k-1))
               else
                  PvscR_ad = 0d0
               end if
               call get_Pvsc_ad(s, k, PvscL_ad, ierr)
               if (ierr /= 0) return
               ! always time center
               PvscL_ad = 0.5d0*(PvscL_ad + s% Pvsc_start(k))
               Pvsc_ad = alfa*PvscL_ad + beta*PvscR_ad
            end if
         
            ! set Ptrb_ad
            if (.not. s% TDC_flag) then
               Ptrb_ad = 0d0
            else
               if (k > 1) then 
                  call calc_Ptrb_ad_tw(s, k-1, PtrbR_ad, ierr)
                  if (ierr /= 0) return
                  PtrbR_ad = shift_m1(PtrbR_ad)
               else
                  PtrbR_ad = 0d0
               end if
               call calc_Ptrb_ad_tw(s, k, PtrbL_ad, ierr)
               if (ierr /= 0) return
               Ptrb_ad = alfa*PtrbL_ad + beta*PtrbR_ad
            end if
         
            ! set extra_P
            if (.not. s% use_other_pressure) then
               extra_P = 0d0
            else if (k > 1) then 
               extra_P = alfa*s% extra_pressure(k) + beta*s% extra_pressure(k-1) 
            else
               extra_P = s% extra_pressure(k)
            end if
         
            ! set mlt_Pturb_ad
            mlt_Pturb_ad = 0d0
            if (s% mlt_Pturb_factor > 0d0 .and. s% mlt_vc_start(k) > 0d0 .and. k > 1) then
               mlt_Pturb_ad%val = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*(s% rho(k-1) + s% rho(k))/6d0
               mlt_Pturb_ad%d1Array(i_lnd_m1) = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*s% rho(k-1)/6d0
               mlt_Pturb_ad%d1Array(i_lnd_00) = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*s% rho(k)/6d0
            end if            
         
            P_face_ad = Peos_ad + Pvsc_ad + Ptrb_ad + mlt_Pturb_ad + extra_P
         
         end if
         
         work_ad = A_times_v_face_ad*P_face_ad
         work = work_ad%val
         
         if (k == 1) s% work_outward_at_surface = work
         
         Av_face = A_times_v_face_ad%val
         do j=1,s% species
            d_work_dxa00(j) = Av_face*d_Pface_dxa00(j)
            d_work_dxam1(j) = Av_face*d_Pface_dxam1(j)
         end do

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
            
         if (test_partials) then
            s% solver_test_partials_val = 0
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'eval1_work', s% solver_test_partials_var
         end if
         
      end subroutine eval1_work
      
      
      subroutine eval1_A_times_v_face_ad(s, k, A_times_v_face_ad, ierr)
         use star_utils, only: get_area_info
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: A_times_v_face_ad
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: A_ad, inv_R2, u_face_ad
         include 'formats'

         ierr = 0
         call get_area_info(s, k, A_ad, inv_R2, ierr)
         if (ierr /= 0) return
         
         u_face_ad = 0d0
         if (s% v_flag) then
            u_face_ad%val = s% vc(k)
            u_face_ad%d1Array(i_v_00) = s% d_vc_dv
         else if (s% u_flag) then
            u_face_ad = s% u_face_ad(k)
            if (s% using_velocity_time_centering) &
               u_face_ad = 0.5d0*(u_face_ad + s% u_face_start(k))
         else if (s% using_velocity_time_centering) then
            u_face_ad%val = 0.5d0*(s% r(k) - s% r_start(k))/s% dt
            u_face_ad%d1Array(i_lnR_00) = 0.5d0*s% r(k)/s% dt
         else
            u_face_ad%val = (s% r(k) - s% r_start(k))/s% dt
            u_face_ad%d1Array(i_lnR_00) = s% r(k)/s% dt
         end if
         
         A_times_v_face_ad = A_ad*u_face_ad
      
      end subroutine eval1_A_times_v_face_ad


      subroutine eval_Fraley_PdV_work( &
            s, k, skip_P, dwork_ad, dwork, d_dwork_dxa00, ierr) 
         use accurate_sum_auto_diff_star_order1
         use auto_diff_support
         use star_utils, only: calc_Ptot_ad_tw
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         logical, intent(in) :: skip_P
         type(auto_diff_real_star_order1), intent(out) :: dwork_ad
         real(dp), intent(out) :: dwork
         real(dp), intent(out), dimension(s% species) :: d_dwork_dxa00
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: &
            Av_face00_ad, Av_facep1_ad, Ptot_ad, dV
         real(dp), dimension(s% species) :: d_Ptot_dxa
         real(dp) :: Av_face00, Av_facep1
         logical :: include_mlt_Pturb
         integer :: j

         include 'formats'
         ierr = 0
         
         ! dV = 1/rho - 1/rho_start 
         call eval1_A_times_v_face_ad(s, k, Av_face00_ad, ierr)
         if (ierr /= 0) return
         if (k < s% nz) then
            call eval1_A_times_v_face_ad(s, k+1, Av_facep1_ad, ierr)
            if (ierr /= 0) return
            Av_facep1_ad = shift_p1(Av_facep1_ad)
         else
            Av_facep1_ad = 0d0
            Av_facep1_ad%val = 4*pi*pow2(s% r_center)*s% v_center
         end if
         Av_face00 = Av_face00_ad%val
         Av_facep1 = Av_facep1_ad%val
         dV = Av_face00_ad - Av_facep1_ad

         include_mlt_Pturb = s% mlt_Pturb_factor > 0d0 &
            .and. s% mlt_vc_start(k) > 0d0 .and. k > 1
         
         call calc_Ptot_ad_tw( &
            s, k, skip_P, .not. include_mlt_Pturb, Ptot_ad, d_Ptot_dxa, ierr)
         if (ierr /= 0) return
         
         do j=1,s% species
            d_dwork_dxa00(j) = d_Ptot_dxa(j)*(Av_face00 - Av_facep1)
         end do
         if (k == 1) s% work_outward_at_surface = Ptot_ad%val*Av_face00
         
         dwork_ad = Ptot_ad*dV
         dwork = dwork_ad%val

      end subroutine eval_Fraley_PdV_work

      
      end module hydro_energy

