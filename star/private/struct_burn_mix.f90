! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module struct_burn_mix

      use star_private_def
      use const_def
      use utils_lib, only: is_bad

      implicit none

      private
      public :: do_struct_burn_mix

      contains


      integer function do_struct_burn_mix(s, skip_global_corr_coeff_limit)
         use mix_info, only: set_mixing_info, get_convection_sigmas
         use rates_def, only: num_rvs
         use hydro_vars, only: set_vars_if_needed
         use star_utils, only: start_time, update_time

         type (star_info), pointer :: s
         logical, intent(in) :: skip_global_corr_coeff_limit

         integer :: nz, nvar, species, ierr, j, k, k_bad
         integer(8) :: time0
         logical :: do_chem
         real(dp) :: dt, tol_correction_norm, tol_max_correction, total

         include 'formats'

         ierr = 0
         s% non_epsnuc_energy_change_from_split_burn = 0d0
         dt = s% dt

         if (s% rsp_flag) then
            do_struct_burn_mix = do_rsp_step(s,dt)    
            s% total_num_solver_iterations = &
               s% total_num_solver_iterations + s% num_solver_iterations
            s% total_num_solver_calls_made = s% total_num_solver_calls_made + 1
            if (do_struct_burn_mix == keep_going) &
               s% total_num_solver_calls_converged = &
                  s% total_num_solver_calls_converged + 1
            return
         end if        

         if (s% use_other_before_struct_burn_mix) then
            call s% other_before_struct_burn_mix(s% id, dt, do_struct_burn_mix)
            if (do_struct_burn_mix /= keep_going) return
         end if

         if (s% doing_timing) call start_time(s, time0, total)

         s% doing_struct_burn_mix = .true.
         nz = s% nz

         species = s% species
         s% num_solver_iterations = 0

         do_struct_burn_mix = retry
         
         s% dVARdot_dVAR = 1d0/dt     

         s% do_burn = (s% dxdt_nuc_factor > 0d0)
         s% do_mix = (s% mix_factor > 0d0)
         
         if (s% op_split_burn) then
            do k=1,nz
               if (s% T(k) >= s% op_split_burn_min_T) then
                  s% burn_num_iters(k) = 0
                  s% eps_nuc(k) = 0d0
                  s% d_epsnuc_dlnd(k) = 0d0
                  s% d_epsnuc_dlnT(k) = 0d0
                  s% d_epsnuc_dx(:,k) = 0d0
                  s% eps_nuc_categories(:,k) = 0d0
                  s% dxdt_nuc(:,k) =  0d0
                  s% d_dxdt_nuc_dRho(:,k) =  0d0
                  s% d_dxdt_nuc_dT(:,k) =  0d0
                  s% d_dxdt_nuc_dx(:,:,k) =  0d0
                  s% eps_nuc_neu_total(k) = 0d0
               end if
            end do
         end if
         
         if (s% do_burn .and. s% op_split_burn) then
            total = 0
            do k=1,s% nz
               total = total - s% energy(k)*s% dm(k)
            end do
            if (s% trace_evolve) write(*,*) 'call do_burn'
            do_struct_burn_mix = do_burn(s, dt)
            if (do_struct_burn_mix /= keep_going) then
               write(*,2) 'failed in do_burn', s% model_number
               stop 'do_struct_burn_mix'
               return
            end if            
            call set_vars_if_needed(s, s% dt, 'after do_burn', ierr)
            if (ierr /= 0) return            
            do k=1,s% nz
               total = total + s% energy(k)*s% dm(k)
            end do
            s% non_epsnuc_energy_change_from_split_burn = total
            if (s% trace_evolve) write(*,*) 'done do_burn'
         end if
         
         if (s% doing_first_model_of_run) then
            if (s% i_lum /= 0) then
               s% L_phot_old = s% xh(s% i_lum,1)/Lsun
            else
               s% L_phot_old = 0
            end if
         end if

         if (s% do_mix) then
            call get_convection_sigmas(s, dt, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'get_convection_sigmas failed'
               return
            end if
         else
            s% sig(1:s% nz) = 0
         end if

         do_chem = (s% do_burn .or. s% do_mix)
         if (do_chem) then ! include abundances
            nvar = s% nvar_total
         else ! no chem => just do structure
            nvar = s% nvar_hydro
         end if

         call set_tol_correction(s, maxval(s% T(1:s% nz)), &
            tol_correction_norm, tol_max_correction)
         call set_surf_info(s, nvar)

         if (s% w_div_wc_flag) then
            s% xh(s% i_w_div_wc,:s% nz) = s% w_div_w_crit_roche(:s% nz)
         end if
         
         if (s% j_rot_flag) then
            s% xh(s% i_j_rot,:s% nz) = s% j_rot(:s% nz)
            s% j_rot_start(:s% nz) = s% j_rot(:s% nz)
            s% i_rot_start(:s% nz) = s% i_rot(:s% nz)
            s% total_abs_angular_momentum = dot_product(abs(s% j_rot(:s% nz)),s% dm_bar(:s% nz))
         end if

         call save_start_values(s, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'save_start_values failed'
            return
         end if
                     
         if (s% trace_evolve) write(*,*) 'call solver'
         do_struct_burn_mix = do_solver_converge( &
            s, nvar, skip_global_corr_coeff_limit, &
            tol_correction_norm, tol_max_correction)
         if (s% trace_evolve) write(*,*) 'done solver'

         s% total_num_solver_iterations = &
            s% total_num_solver_iterations + s% num_solver_iterations
         s% total_num_solver_calls_made = s% total_num_solver_calls_made + 1
         if (do_struct_burn_mix == keep_going) &
            s% total_num_solver_calls_converged = &
               s% total_num_solver_calls_converged + 1
         if (s% doing_relax) then
            s% total_num_solver_relax_iterations = &
               s% total_num_solver_relax_iterations + s% num_solver_iterations
            s% total_num_solver_relax_calls_made = s% total_num_solver_relax_calls_made + 1
            if (do_struct_burn_mix == keep_going) &
               s% total_num_solver_relax_calls_converged = &
                  s% total_num_solver_relax_calls_converged + 1
         end if

         if (do_struct_burn_mix /= keep_going) return

         if (.not. s% j_rot_flag) &
            do_struct_burn_mix = do_mix_omega(s,dt)

         if (s% use_other_after_struct_burn_mix) &
            call s% other_after_struct_burn_mix(s% id, dt, do_struct_burn_mix)

         s% solver_iter = 0 ! to indicate that no longer doing solver iterations
         s% doing_struct_burn_mix = .false.
         if (s% doing_timing) call update_time(s, time0, total, s% time_struct_burn_mix)
         
         contains
         
         subroutine test(str)
            use chem_def, only: category_name
            character (len=*), intent(in) :: str
            include 'formats'
            integer :: k, i
            k = s% nz
            i = maxloc(s% eps_nuc_categories(1:num_categories,k),dim=1)
            write(*,3) trim(str) // ' eps_nuc_cat ' // trim(category_name(i)), &
               k, s% model_number, s% eps_nuc_categories(i,k)
         end subroutine test

      end function do_struct_burn_mix


      integer function do_mix_omega(s, dt)
         use solve_omega_mix, only: do_solve_omega_mix
         use hydro_rotation, only: set_i_rot, set_omega

         type (star_info), pointer :: s
         real(dp), intent(in) :: dt

         integer :: ierr

         do_mix_omega = keep_going

         if (s% rotation_flag) then
            ! After solver is done with the structure, recompute moments of inertia and
            ! omega before angular momentum mix
            call set_i_rot(s, .false.)
            call set_omega(s, 'struct_burn_mix')
            if (s% premix_omega) then
               do_mix_omega = do_solve_omega_mix(s, 0.5d0*dt)
            else
               do_mix_omega = do_solve_omega_mix(s, dt)
            end if
            if (do_mix_omega /= keep_going) return
         end if

      end function do_mix_omega


      integer function do_rsp_step(s,dt)
         ! return keep_going, retry, or terminate
         use rsp, only: rsp_one_step
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer :: ierr
         do_rsp_step = keep_going
         ierr = 0
         call rsp_one_step(s,ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'ierr from rsp_one_step'
            do_rsp_step = retry
         end if
      end function do_rsp_step


      subroutine save_start_values(s, ierr)
         use chem_def, only: num_categories
         use hydro_rsp2, only: set_etrb_start_vars
         use star_utils, only: eval_total_energy_integrals, set_luminosity_by_category
         use chem_def, only: ih1
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k, j, i_h1
         include 'formats'
         ierr = 0

         call set_luminosity_by_category(s)

         do k=1,s% nz
            do j=1,s% species
               s% dxdt_nuc_start(j,k) = s% dxdt_nuc(j,k)
            end do
            do j=1,num_categories
               s% luminosity_by_category_start(j,k) = &
                  s% luminosity_by_category(j,k)
            end do
         end do

         do k=1,s% nz
            !s% lnT_start(k) = s% lnT(k)
            !s% T_start(k) set elsewhere
            !s% lnd_start(k) = s% lnd(k) set elsewhere
            !s% rho_start(k) = s% rho(k)
            !s% r_start(k) set elsewhere
            !s% rmid_start(k) set elsewhere
            !s% v_start(k) set elsewhere
            !s% csound_start(k) set elsewhere
            s% lnPeos_start(k) = s% lnPeos(k)
            s% Peos_start(k) = s% Peos(k)
            s% lnPgas_start(k) = s% lnPgas(k)
            s% lnE_start(k) = s% lnE(k)
            s% energy_start(k) = s% energy(k)
            s% lnR_start(k) = s% lnR(k)
            s% u_start(k) = s% u(k)
            s% u_face_start(k) = 0d0 ! s% u_face_ad(k)%val
            s% P_face_start(k) = -1d0 ! mark as unset s% P_face_ad(k)%val
            s% L_start(k) = s% L(k)
            s% omega_start(k) = s% omega(k)
            s% ye_start(k) = s% ye(k)
            s% X_start(k) = s% X(k)
            s% Y_start(k) = s% Y(k)
            s% Z_start(k) = s% Z(k)
            s% j_rot_start(k) = s% j_rot(k)
            s% eps_nuc_start(k) = s% eps_nuc(k)
            s% non_nuc_neu_start(k) = s% non_nuc_neu(k)
            s% mass_correction_start(k) = s% mass_correction(k)
            s% P_div_rho_start(k) = s% Peos(k)/s% rho(k)
            s% Pvsc_start(k) = -1d99
            s% scale_height_start(k) = s% scale_height(k)
            s% gradT_start(k) = s% gradT(k)
            s% gradL_start(k) = s% gradL(k)
            s% grada_start(k) = s% grada(k)
            s% gradr_start(k) = s% gradr(k)
            s% grada_face_start(k) = s% grada_face(k)
            s% chiT_start(k) = s% chiT(k)
            s% chiRho_start(k) = s% chiRho(k)
            s% cp_start(k) = s% cp(k)
            s% Cv_start(k) = s% Cv(k)
            s% dE_dRho_start(k) = s% dE_dRho(k)
            s% gam_start(k) = s% gam(k)
            s% lnS_start(k) = s% lnS(k)
            s% eta_start(k) = s% eta(k)
            s% abar_start(k) = s% abar(k)
            s% zbar_start(k) = s% zbar(k)
            s% z53bar_start(k) = s% z53bar(k)
            s% mu_start(k) = s% mu(k)
            s% phase_start(k) = s% phase(k)
            s% latent_ddlnT_start(k) = s% latent_ddlnT(k)
            s% latent_ddlnRho_start(k) = s% latent_ddlnRho(k)
            s% eps_nuc_start(k) = s% eps_nuc(k)
            s% opacity_start(k) = s% opacity(k)
            s% mlt_mixing_length_start(k) = s% mlt_mixing_length(k)
            s% mlt_mixing_type_start(k) = s% mlt_mixing_type(k)
            s% mlt_D_start(k) = s% mlt_D(k)
            s% mlt_Gamma_start(k) = s% mlt_Gamma(k)
         end do
         
         if (s% using_RSP2) then
            call set_etrb_start_vars(s,ierr)
         end if

         do k=1,s% nz
            do j=1,s% nvar_hydro
               s% xh_start(j,k) = s% xh(j,k)
            end do
         end do
         
         do k=1,s% nz
            do j=1,s% species
               s% xa_start(j,k) = s% xa(j,k)
            end do
         end do
         
         s% start_H_envelope_base_k = s% nz+1
         i_h1 = s% net_iso(ih1)
         if (i_h1 > 0) then
            ! start_H_envelope_base_k = outermost cell where H1 is not most abundant species
            do k=1,s% nz
               j = maxloc(s% xa(1:s% species,k), dim=1)
               if (j /= i_h1) then
                  s% start_H_envelope_base_k = k
                  exit
               end if
            end do
         end if
         if (s% start_H_envelope_base_k > s% nz) s% start_H_envelope_base_k = 0

         call eval_total_energy_integrals(s, &
            s% total_internal_energy_start, &
            s% total_gravitational_energy_start, &
            s% total_radial_kinetic_energy_start, &
            s% total_rotational_kinetic_energy_start, &
            s% total_turbulent_energy_start, &
            s% total_energy_start)
         
      end subroutine save_start_values


      integer function do_solver_converge( &
            s, nvar, skip_global_corr_coeff_limit, &
            tol_correction_norm, tol_max_correction)
         ! return keep_going, retry, or terminate
         use mtx_lib
         use mtx_def
         use num_def
         use star_utils, only: start_time, update_time

         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         logical, intent(in) :: skip_global_corr_coeff_limit
         real(dp), intent(in) :: tol_correction_norm, tol_max_correction

         integer :: ierr, nz, k, n, solver_lwork, solver_liwork
         logical :: report

         include 'formats'

         if (s% dt <= 0d0) then
            do_solver_converge = keep_going
            return
         end if

         do_solver_converge = terminate

         ierr = 0

         nz = s% nz
         n = nz*nvar
         call work_sizes_for_solver(ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'do_solver_converge: work_sizes_for_solver failed'
            do_solver_converge = retry
            s% result_reason = nonzero_ierr
            s% termination_code = t_solve_hydro
            return
         end if

         call alloc_for_solver(ierr)
         if (ierr /= 0) then
            s% termination_code = t_solve_hydro
            s% result_reason = nonzero_ierr
            return
         end if

         report = (s% report_solver_progress .or. s% report_ierr)
         
         s% solver_call_number = s% solver_call_number + 1

         do_solver_converge = do_solver( &
            s, skip_global_corr_coeff_limit, &
            tol_correction_norm, tol_max_correction, &
            report, nz, nvar, s% solver_work, solver_lwork, &
            s% solver_iwork, solver_liwork)


         contains


         subroutine alloc_for_solver(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            if (.not. associated(s% solver_iwork)) then
               allocate(s% solver_iwork(solver_liwork))
            else if (size(s% solver_iwork, dim=1) < solver_liwork) then
               deallocate(s% solver_iwork)
               allocate(s% solver_iwork(int(1.3d0*solver_liwork)+100))
            end if
            if (.not. associated(s% solver_work)) then
               allocate(s% solver_work(solver_lwork))
            else if (size(s% solver_work, dim=1) < solver_lwork) then
               deallocate(s% solver_work)
               allocate(s% solver_work(int(1.3d0*solver_lwork)+100))
            end if
         end subroutine alloc_for_solver

         subroutine work_sizes_for_solver(ierr)
            use star_solver, only: get_solver_work_sizes
            integer, intent(out) :: ierr
            call get_solver_work_sizes(s, nvar, nz, solver_lwork, solver_liwork, ierr)
         end subroutine work_sizes_for_solver

      end function do_solver_converge


      subroutine set_surf_info(s, nvar) ! set to values at start of step
         use star_utils, only: get_lnd_from_xh, get_lnT_from_xh, get_lnR_from_xh
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         if (s% i_lnd > 0 .and. s% i_lnd <= nvar) &
            s% surf_lnd = get_lnd_from_xh(s, 1)
         if (s% i_lnT > 0 .and. s% i_lnT <= nvar) &
            s% surf_lnT = get_lnT_from_xh(s, 1)
         if (s% i_lnR > 0 .and. s% i_lnR <= nvar) &
            s% surf_lnR = get_lnR_from_xh(s, 1)
         if (s% i_v > 0 .and. s% i_v <= nvar) &
            s% surf_v = s% xh(s% i_v,1)
         if (s% i_u > 0 .and. s% i_u <= nvar) &
            s% surf_v = s% xh(s% i_u,1)
         s% surf_lnS = s% lnS(1)
         s% num_surf_revisions = 0
      end subroutine set_surf_info


      subroutine set_xh(s,nvar) ! set xh using current structure info
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         integer :: j1, k, nz
         include 'formats'
         nz = s% nz
         do j1 = 1, min(nvar,s% nvar_hydro)
            if (j1 == s% i_lnd .and. s% i_lnd <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% lnd(k)
               end do
            else if (j1 == s% i_lnT .and. s% i_lnT <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% lnT(k)
               end do
            else if (j1 == s% i_lnR .and. s% i_lnR <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% lnR(k)
               end do
            else if (j1 == s% i_lum .and. s% i_lum <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% L(k)
               end do
            else if (j1 == s% i_w .and. s% i_w <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% w(k)
               end do
            else if (j1 == s% i_Hp .and. s% i_Hp <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% Hp_face(k)
               end do
            else if (j1 == s% i_v .and. s% i_v <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% v(k)
               end do
            else if (j1 == s% i_u .and. s% i_u <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% u(k)
               end do
            else if (j1 == s% i_alpha_RTI .and. s% i_alpha_RTI <= nvar) then
               do k = 1, nz
                  s% xh(j1,k) = s% alpha_RTI(k)
               end do
            else if (j1 == s% i_ln_cvpv0 .and. s% i_ln_cvpv0 <= nvar) then
                  stop 'pablo needs to revise this'
!               do k = 1, nz
!                  ! create a rough first guess using mlt_vc_start and conv_vel when
!                  ! mlt_vc is larger than the starting conv_vel
!                  if (s% mlt_vc_start(k) > 0d0 .and. s% mlt_vc_start(k) > s% conv_vel(k)) then
!                     s% conv_vel(k) = s% conv_vel(k) + &
!                        (s% mlt_vc_start(k) -s% conv_vel(k)) * &
!                        min(1d0, s% dt*s% mlt_vc_start(k)/(s% scale_height_start(k)*s% mixing_length_alpha))
!                  end if
!                  s% xh(j1,k) = log(s% conv_vel(k)+s% conv_vel_v0)
!               end do
            end if
         end do
      end subroutine set_xh


      subroutine set_tol_correction( &
            s, T_max, tol_correction_norm, tol_max_correction)
         type (star_info), pointer :: s
         real(dp), intent(in) :: T_max
         real(dp), intent(out) :: tol_correction_norm, tol_max_correction
         include 'formats'
         if (T_max >= s% tol_correction_extreme_T_limit) then
            tol_correction_norm = s% tol_correction_norm_extreme_T
            tol_max_correction = s% tol_max_correction_extreme_T
         else if (T_max >= s% tol_correction_high_T_limit) then
            tol_correction_norm = s% tol_correction_norm_high_T
            tol_max_correction = s% tol_max_correction_high_T
         else
            tol_correction_norm = s% tol_correction_norm
            tol_max_correction = s% tol_max_correction
         end if
      end subroutine set_tol_correction


      integer function do_solver( &
            s, skip_global_corr_coeff_limit, &
            tol_correction_norm, tol_max_correction, &
            report, nz, nvar, solver_work, solver_lwork, &
            solver_iwork, solver_liwork)
         ! return keep_going, retry, or terminate

         ! when using solver for hydro step,
         ! do not require that functions have been evaluated for starting configuration.
         ! when finish, will have functions evaluated for the final set of primary variables.
         ! for example, the reaction rates will have been computed, so they can be used
         ! as initial values in the following burn and mix.

         use num_def
         use alloc

         type (star_info), pointer :: s
         integer, intent(in) :: nz, nvar
         logical, intent(in) :: skip_global_corr_coeff_limit, report
         real(dp), intent(in) :: tol_correction_norm, tol_max_correction

         integer, intent(in) :: solver_lwork, solver_liwork
         real(dp), pointer :: solver_work(:) ! (solver_lwork)
         integer, pointer :: solver_iwork(:) ! (solver_liwork)
         logical :: converged
         integer :: i, k, species, ierr, alph, j1, j2, gold_tolerances_level
         real(dp) :: varscale, r003, rp13, dV, frac, maxT

         include 'formats'

         species = s% species
         do_solver = keep_going
         s% using_gold_tolerances = .false.
         gold_tolerances_level = 0
         
         if ((s% use_gold2_tolerances .and. s% steps_before_use_gold2_tolerances < 0) .or. &
             (s% steps_before_use_gold2_tolerances >= 0 .and. &
                s% model_number > s% steps_before_use_gold2_tolerances + max(0,s% init_model_number))) then
            s% using_gold_tolerances = .true.
            gold_tolerances_level = 2
         else if ((s% use_gold_tolerances .and. s% steps_before_use_gold_tolerances < 0) .or. &
             (s% steps_before_use_gold_tolerances >= 0 .and. &
                s% model_number > s% steps_before_use_gold_tolerances + max(0,s% init_model_number))) then
            if (s% maxT_for_gold_tolerances > 0) then
               maxT = maxval(s% T(1:nz))
            else
               maxT = -1d0
            end if
            if (maxT > s% maxT_for_gold_tolerances) then 
               !write(*,2) 'exceed maxT_for_gold_tolerances', &
               !   s% model_number, maxT, s% maxT_for_gold_tolerances
            else ! okay for maxT, so check if also ok for eosPC_frac
               s% using_gold_tolerances = .true.
               gold_tolerances_level = 1
            end if
         end if

         call set_xh(s, nvar) ! set xh using current structure info

         do k = 1, nz
            do j1 = 1, min(nvar, s% nvar_hydro)
               s% solver_dx(j1,k) = s% xh(j1,k) - s% xh_start(j1,k)
            end do
         end do

         if (nvar >= s% nvar_hydro+1) then
            do k = 1, nz
               j2 = 1
               do j1 = s% nvar_hydro+1, nvar
                  s% xa_sub_xa_start(j2,k) = s% xa(j2,k) - s% xa_start(j2,k)
                  s% solver_dx(j1,k) = s% xa_sub_xa_start(j2,k)
                  j2 = j2+1
               end do
            end do
         end if

         converged = .false.
         call hydro_solver_step( &
            s, nz, s% nvar_hydro, nvar, skip_global_corr_coeff_limit, &
            gold_tolerances_level, tol_max_correction, tol_correction_norm, &
            solver_work, solver_lwork, &
            solver_iwork, solver_liwork, &
            converged, ierr)

         if (ierr /= 0) then
            if (report) then
               write(*, *) 'hydro_solver_step returned ierr', ierr
               write(*, *) 's% model_number', s% model_number
               write(*, *) 'nz', nz
               write(*, *) 's% num_retries', s% num_retries
               write(*, *)
            end if
            do_solver = retry
            s% result_reason = nonzero_ierr
            s% dt_why_retry_count(Tlim_solver) = &
               s% dt_why_retry_count(Tlim_solver) + 1
            return
         end if

         if (converged) then ! sanity checks before accept it
            converged = check_after_converge(s, report, ierr)
            if (converged .and. s% RTI_flag) & ! special checks
               converged = RTI_check_after_converge(s, report, ierr)
         end if

         if (.not. converged) then
            do_solver = retry
            s% result_reason = hydro_failed_to_converge
            s% dt_why_retry_count(Tlim_solver) = &
               s% dt_why_retry_count(Tlim_solver) + 1
            if (report) then
               write(*,2) 'solver rejected trial model'
               write(*,2) 's% model_number', s% model_number
               write(*,2) 's% solver_call_number', s% solver_call_number
               write(*,2) 'nz', nz
               write(*,2) 's% num_retries', s% num_retries
               write(*,1) 'dt', s% dt
               write(*,1) 'log dt/secyer', log10(s% dt/secyer)
               write(*, *)
            end if
            return
         end if

      end function do_solver


      logical function RTI_check_after_converge(s, report, ierr) result(converged)
         use mesh_adjust, only: set_lnT_for_energy
         use micro, only: do_eos_for_cell
         use chem_def, only: ih1, ihe3, ihe4
         use star_utils, only: store_lnT_in_xh, get_T_and_lnT_from_xh
         type (star_info), pointer :: s
         logical, intent(in) :: report
         integer, intent(out) :: ierr
         integer :: k, nz
         real(dp) :: old_energy, old_IE, new_IE, old_KE, new_KE, new_u, new_v, &
            revised_energy, new_lnT
         include 'formats'
         ierr = 0
         nz = s% nz
         converged = .true.
         !return
         do k=1,nz
            if (k < nz .and. s% alpha_RTI(k) < 1d-10) cycle
            old_energy = s% energy(k)
            old_IE = old_energy*s% dm(k)
            if (s% energy(k) < s% RTI_energy_floor) then
               ! try to take from KE to give to IE
               ! else just bump energy and take hit to energy conservation
               s% energy(k) = s% RTI_energy_floor
               s% lnE(k) = log(s% energy(k))
               call set_lnT_for_energy(s, k, &
                  s% net_iso(ih1), s% net_iso(ihe3), s% net_iso(ihe4), &
                  s% species, s% xa(:,k), &
                  s% rho(k), s% lnd(k)/ln10, s% energy(k), s% lnT(k), &
                  new_lnT, revised_energy, ierr)
               if (ierr /= 0) return ! stop 'do_merge failed in set_lnT_for_energy'
               call store_lnT_in_xh(s, k, new_lnT)
               call get_T_and_lnT_from_xh(s, k, s% T(k), s% lnT(k))
            end if
            new_IE = s% energy(k)*s% dm(k)
            if (s% u_flag) then
               old_KE = 0.5d0*s% dm(k)*s% u(k)*s% u(k)
               new_KE = max(0d0, old_KE + old_IE - new_IE)
               new_u = sqrt(new_KE/(0.5d0*s% dm(k)))
               if (s% u(k) > 0d0) then
                  s% u(k) = new_u
               else
                  s% u(k) = -new_u
               end if
               s% xh(s% i_u, k) = s% u(k)
            else if (s% v_flag) then ! only rough approximation possible here
               old_KE = 0.5d0*s% dm_bar(k)*s% v(k)*s% v(k)
               new_KE = max(0d0, old_KE + old_IE - new_IE)
               new_v = sqrt(max(0d0,new_KE)/(0.5d0*s% dm_bar(k)))
               if (s% v(k) > 0d0) then
                  s% v(k) = new_v
               else
                  s% v(k) = -new_v
               end if
               s% xh(s% i_v, k) = s% v(k)
            end if
         end do
      end function RTI_check_after_converge


      logical function check_after_converge(s, report, ierr) result(converged)
         type (star_info), pointer :: s
         logical, intent(in) :: report
         integer, intent(out) :: ierr
         integer :: k, nz
         include 'formats'
         ierr = 0
         nz = s% nz
         converged = .true.
         if (s% R_center > 0) then
            if (s% R_center > exp(s% lnR(nz))) then
               if (report) &
                  write(*,2) 'volume < 0 in cell nz', nz, &
                     s% R_center - exp(s% lnR(nz)), s% R_center, exp(s% lnR(nz)), &
                     s% dm(nz), s% rho(nz), s% dq(nz)
               converged = .false.
               return
            end if
         end if
         do k=1,nz-1
            if (s% lnR(k) <= s% lnR(k+1)) then
               if (report) write(*,2) 'after hydro, negative cell volume in cell k', &
                     k, s% lnR(k) - s% lnR(k+1), s% lnR(k), s% lnR(k+1), &
                     s% lnR_start(k) - s% lnR_start(k+1), s% lnR_start(k), s% lnR_start(k+1)
               converged = .false.; exit
               stop 'check_after_converge'
            else
               if (s% lnT(k) > ln10*12) then
                  if (report) write(*,2) 'after hydro, logT > 12 in cell k', k, s% lnT(k)
                  converged = .false.!; exit
               else if (s% lnT(k) < ln10) then
                  if (report) write(*,*) 'after hydro, logT < 1 in cell k', k
                  converged = .false.!; exit
               else if (s% lnd(k) > ln10*12) then
                  if (report) write(*,*) 'after hydro, logRho > 12 in cell k', k
                  converged = .false.!; exit
               else if (s% lnd(k) < -ln10*30) then
                  if (report) write(*,*) 'after hydro, logRho < -30 in cell k', k
                  converged = .false.!; exit
               end if
            end if
         end do
      end function check_after_converge


      subroutine hydro_solver_step( &
            s, nz, nvar_hydro, nvar, skip_global_corr_coeff_limit, &
            gold_tolerances_level, tol_max_correction, tol_correction_norm, &
            solver_work, solver_lwork, &
            solver_iwork, solver_liwork, &
            converged, ierr)
         use num_def
         use chem_def
         use mtx_lib
         use mtx_def
         use alloc

         type (star_info), pointer :: s
         integer, intent(in) :: nz, nvar_hydro, nvar
         logical, intent(in) :: skip_global_corr_coeff_limit
         real(dp), intent(in) :: tol_max_correction, tol_correction_norm
         integer, intent(in) :: gold_tolerances_level
         integer, intent(in) :: solver_lwork, solver_liwork
         real(dp), intent(inout), pointer :: solver_work(:) ! (solver_lwork)
         integer, intent(inout), pointer :: solver_iwork(:) ! (solver_liwork)
         logical, intent(out) :: converged
         integer, intent(out) :: ierr

         integer :: i, k, j, matrix_type, neq
         logical :: failure
         real(dp) :: varscale
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         neq = nvar*nz

         if (dbg) write(*, *) 'enter hydro_solver_step'

         s% used_extra_iter_in_solver_for_accretion = .false.

         call check_sizes(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'check_sizes failed'
            return
         end if

         if (dbg) write(*, *) 'call solver'
         call newt(ierr)
         if (ierr /= 0 .and. s% report_ierr) then
            write(*,*) 'solver failed for hydro'
         end if

         converged = (ierr == 0) .and. (.not. failure)
         if (converged) then
            do k=1,nz
               do j=1,min(nvar,nvar_hydro)
                  s% xh(j,k) = s% xh_start(j,k) + s% solver_dx(j,k)
               end do
            end do
            ! s% xa has already been updated by final call to set_solver_vars from solver
         end if


         contains


         subroutine newt(ierr)
            use star_solver, only: solver
            use rates_def, only: warn_rates_for_high_temp
            integer, intent(out) :: ierr
            integer :: k, j
            logical :: save_warn_rates_flag
            include 'formats'
            s% doing_solver_iterations = .true.
            save_warn_rates_flag = warn_rates_for_high_temp
            warn_rates_for_high_temp = .false.        
            call solver( &
               s, nvar, skip_global_corr_coeff_limit, &
               gold_tolerances_level, tol_max_correction, tol_correction_norm, &
               solver_work, solver_lwork, &
               solver_iwork, solver_liwork, &
               failure, ierr)
            s% doing_solver_iterations = .false.
            warn_rates_for_high_temp = save_warn_rates_flag
         end subroutine newt


      end subroutine hydro_solver_step


      integer function do_burn(s, dt)
         use net_lib, only: net_work_size, net_1_zone_burn_const_density_work_size, &
            net_1_zone_burn_work_size
         use star_utils, only: start_time, update_time
         use net, only: get_screening_mode
         use chem_def
         use micro, only: do_eos_for_cell

         type (star_info), pointer :: s
         real(dp), intent(in) :: dt

         integer :: &
            k_bad, net_lwork, ierr, max_num_iters_k, nz, op_err, &
            i, j, k, num_iters, species, max_num_iters_used, &
            screening_mode, burn_lwork, burn_lwork_const_density, kmin
         integer(8) :: time0, clock_rate
         real(dp) :: total, avg_epsnuc, min_T_for_const_density_solver
         logical :: trace, dbg, okay, skip_burn
         logical, parameter :: burn_dbg = .false.

         include 'formats'

         trace = .false.
         
         min_T_for_const_density_solver = s% op_split_burn_min_T_for_variable_T_solver

         do_burn = keep_going
         ierr = 0
         nz = s% nz
         species = s% species

         if (s% eps_nuc_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) then
            do k = 1, nz
               s% eps_nuc(k) = 0d0
               s% burn_num_iters(k) = 0
               s% burn_avg_epsnuc(k) = 0d0
               s% max_burn_correction(k) = 0d0
            end do
            return
         end if

         if (dt <= 0d0) return
         
         net_lwork = net_work_size(s% net_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_burn failed in net_work_size'
            do_burn = terminate
            s% termination_code = t_solve_burn
            s% result_reason = nonzero_ierr
            return
         end if
         
         burn_lwork_const_density = net_1_zone_burn_const_density_work_size(s% net_handle,ierr)
         if (ierr /= 0) then
            write(*,*) 'do_burn failed in net_1_zone_burn_const_density_work_size'
            do_burn = terminate
            s% termination_code = t_solve_burn
            s% result_reason = nonzero_ierr
            return
         end if
         burn_lwork = net_1_zone_burn_work_size(s% net_handle,ierr)
         if (ierr /= 0) then
            write(*,*) 'do_burn failed in net_1_zone_burn_work_size'
            do_burn = terminate
            s% termination_code = t_solve_burn
            s% result_reason = nonzero_ierr
            return
         end if
         
         burn_lwork = max(burn_lwork, burn_lwork_const_density)

         max_num_iters_used = 0
         max_num_iters_k = 0
         k_bad = 0
         
         screening_mode = get_screening_mode(s,ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,*) 'unknown string for screening_mode: ' // trim(s% screening_mode)
            return
            stop 'do1_net'
         end if

         dbg = .false. ! (s% model_number == 1137)

         kmin = nz+1
         do k=1,nz
            if (s% T_start(k) < s% op_split_burn_min_T) then
                ! We get here if we have an off center ignition,
                ! the arrays wont have been initialised earlier as they stop at the 
                ! first temperature that exceeds op_split_burn_min_T
               s% burn_num_iters(k) = 0
               s% burn_avg_epsnuc(k) = 0d0
               cycle
            end if
            kmin = k
            exit
         end do
         
         if (kmin > nz) return

         !skip_burn = s% fe_core_infall > s% op_split_burn_eps_nuc_infall_limit
         skip_burn = (minval(s% v_start(1:s% nz)) < -s% op_split_burn_eps_nuc_infall_limit)

         if (s% doing_timing) call start_time(s, time0, total)

!$OMP PARALLEL DO PRIVATE(k,op_err,num_iters,avg_epsnuc) SCHEDULE(dynamic,2)
         do k = kmin, nz
            if (k_bad /= 0) cycle
            if (s% T_start(k) < s% op_split_burn_min_T) then
               ! We get here if we have an off center ignition,
               ! the arrays wont have been initialised earlier as they stop at the 
               ! first temperature that exceeds op_split_burn_min_T
               s% burn_num_iters(k) = 0
               s% burn_avg_epsnuc(k) = 0d0
               cycle
            end if
            s% max_burn_correction(k) = 0d0
            op_err = 0
            call burn1_zone( &
               s, k, species, min_T_for_const_density_solver, skip_burn, &
               net_lwork, burn_lwork, screening_mode, &
               dt, num_iters, avg_epsnuc, burn_dbg, op_err)
            if (op_err /= 0) then
               ierr = -1
               k_bad = k
               cycle
            end if
            call do_eos_for_cell(s,k,op_err)
            if (op_err /= 0) then
               write(*,2) 'do_burn failed in do_eos_for_cell', k
               ierr = -1
               k_bad = k
               cycle
            end if        
            !write(*,3) 'num_iters', k, num_iters
            s% burn_num_iters(k) = num_iters
            s% burn_avg_epsnuc(k) = avg_epsnuc
            if (num_iters > max_num_iters_used) then
               max_num_iters_used = num_iters
               max_num_iters_k = k
            end if
         end do
!$OMP END PARALLEL DO
         
         s% need_to_setvars = .true.
         
         if (s% doing_timing) &
            call update_time(s, time0, total, s% time_solve_burn)
            
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'do_burn failed', k_bad
            return
            stop 'do_burn'
         
         
            do_burn = retry
            if (trace .or. s% report_ierr) then
               write(*,*) 'do_burn ierr'
               !stop 'do_burn'
            end if
            call restore
            return
         end if

         if (dbg) write(*,2) 'done do_burn'


         contains

         subroutine restore
            integer :: j, k
            do k = 1, nz
               do j=1,species
                  s% xa(j,k) = s% xa_start(j,k)
               end do
            end do
         end subroutine restore

      end function do_burn


      subroutine burn1_zone( &
            s, k, species, min_T_for_const_density_solver, skip_burn, &
            net_lwork, burn_lwork, screening_mode, &
            dt, num_iters_out, avg_epsnuc, dbg_in, ierr)
         use net_lib, only: net_1_zone_burn_const_density, net_1_zone_burn, &
            show_net_reactions_and_info
         use rates_def, only: std_reaction_Qs, std_reaction_neuQs
         use chem_def, only: chem_isos, num_categories, category_name
         use net, only: do1_net
         use star_utils, only: store_lnT_in_xh, get_T_and_lnT_from_xh
         type (star_info), pointer :: s
         integer, intent(in) :: k, species, &
            net_lwork, burn_lwork, screening_mode
         real(dp), intent(in) :: dt, min_T_for_const_density_solver
         logical, intent(in) :: skip_burn, dbg_in
         real(dp), intent(out) :: avg_epsnuc
         integer, intent(out) :: num_iters_out, ierr
         
         real(dp), target :: net_work_ary(net_lwork), burn_work_ary(burn_lwork), xa_start_ary(species)
         real(dp), pointer :: net_work(:), burn_work(:), xa_start(:)
         
         real(dp) :: stptry, eps, odescal, &
            starting_log10T, ending_log10T, ending_eps_neu_total, &
            Cv0, eta0, substep_start_time
         integer :: i, max_steps, nfcn, njac, ntry, naccpt, nrejct
         integer, parameter :: num_times = 1
         real(dp), target, dimension(4*num_times) :: log10Ts_ary, log10Rhos_ary, etas_ary
         real(dp), pointer, dimension(:) :: log10Ts_f1, log10Rhos_f1, etas_f1, &
            dxdt_source_term, times
         logical :: use_pivoting, trace, burn_dbg
         
         include 'formats'

         ierr = 0
         num_iters_out = 0
         
         if (skip_burn) then
            avg_epsnuc = 0d0
            s% eps_nuc(k) = 0d0
            s% d_epsnuc_dlnd(k) = 0d0
            s% d_epsnuc_dlnT(k) = 0d0
            s% d_epsnuc_dx(:,k) = 0d0
            s% dxdt_nuc(:,k) = 0d0
            s% eps_nuc_categories(:,k) = 0d0
            s% d_dxdt_nuc_dRho(:,k) =  0d0
            s% d_dxdt_nuc_dT(:,k) =  0d0
            s% d_dxdt_nuc_dx(:,:,k) =  0d0
            s% eps_nuc_neu_total(k) = 0d0
            return
         end if
         
         log10Ts_f1 => log10Ts_ary
         log10Rhos_f1 => log10Rhos_ary
         etas_f1 => etas_ary
         
         nullify(dxdt_source_term, times)
         
         net_work => net_work_ary
         burn_work => burn_work_ary
         xa_start => xa_start_ary

         stptry = 0d0
         eps = s% op_split_burn_eps
         odescal = s% op_split_burn_odescal         
         max_steps = s% burn_steps_hard_limit         
         use_pivoting = .false. ! .true.
         trace = .false.
         burn_dbg = .false.
         starting_log10T = s% lnT(k)/ln10
         
         do i=1,species
            xa_start(i) = s% xa(i,k)
         end do
         
         substep_start_time = 0d0
         
         if (s% T(k) >= min_T_for_const_density_solver) then
            Cv0 = s% Cv(k)
            eta0 = s% eta(k)
            call net_1_zone_burn_const_density( &
               s% net_handle, s% eos_handle, species, s% num_reactions, 0d0, dt, &
               xa_start, starting_log10T, s% lnd(k)/ln10, &
               get_eos_info_for_burn_at_const_density, &
               s% rate_factors, s% weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, &
               screening_mode, &
               stptry, max_steps, eps, odescal, &
               use_pivoting, trace, burn_dbg, burn_finish_substep, &
               burn_lwork, burn_work, net_lwork, net_work, s% xa(1:species,k), &
               s% eps_nuc_categories(:,k), &
               ending_log10T, avg_epsnuc, ending_eps_neu_total, &
               nfcn, njac, ntry, naccpt, nrejct, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,2) 'net_1_zone_burn_const_density failed', k
               return
               stop 'burn1_zone'
            end if
            ! restore temperature
            call store_lnT_in_xh(s, k, starting_log10T*ln10)
            call get_T_and_lnT_from_xh(s, k, s% T(k), s% lnT(k))
         else
            log10Ts_f1 => log10Ts_ary
            log10Rhos_f1 => log10Rhos_ary
            etas_f1 => etas_ary
            nullify(dxdt_source_term, times)
            log10Ts_f1(1) = s% lnT(k)/ln10
            log10Rhos_f1(1) = s% lnd(k)/ln10
            etas_f1(1) = s% eta(k)
            call net_1_zone_burn( &
               s% net_handle, s% eos_handle, species, s% num_reactions, 0d0, dt, xa_start, &
               num_times, times, log10Ts_f1, log10Rhos_f1, etas_f1, dxdt_source_term, &
               s% rate_factors, s% weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, &
               screening_mode,  &
               stptry, max_steps, eps, odescal, &
               use_pivoting, trace, burn_dbg, burn_finish_substep, &
               burn_lwork, burn_work, net_lwork, net_work, s% xa(1:species,k), &
               s% eps_nuc_categories(:,k), &
               avg_epsnuc, ending_eps_neu_total, &
               nfcn, njac, ntry, naccpt, nrejct, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,2) 'net_1_zone_burn failed', k
               return
               stop 'burn1_zone'
            end if
         end if
         
         num_iters_out = naccpt
         
         ! make extra call to get eps_nuc_categories
         call do1_net(s, k, s% species, s% num_reactions, &
            net_lwork, .false., ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,2) 'net_1_zone_burn final call to do1_net failed', k
            return
            stop 'burn1_zone'
         end if
               
         s% eps_nuc(k) = 0d0
         s% d_epsnuc_dlnd(k) = 0d0
         s% d_epsnuc_dlnT(k) = 0d0
         s% d_epsnuc_dx(:,k) = 0d0
         s% dxdt_nuc(:,k) = 0d0
         !s% eps_nuc_categories(:,k) = 0d0
         s% d_dxdt_nuc_dRho(:,k) =  0d0
         s% d_dxdt_nuc_dT(:,k) =  0d0
         s% d_dxdt_nuc_dx(:,:,k) =  0d0
         s% eps_nuc_neu_total(k) = 0d0
         
         do i=1,species ! for use by dX_nuc_drop timestep limiter
            s% dxdt_nuc(i,k) = (s% xa(i,k)-xa_start(i))/dt
         end do
         
         contains
         
         subroutine get_eos_info_for_burn_at_const_density( &
               eos_handle, species, chem_id, net_iso, xa, &
               Rho, logRho, T, logT, &
               Cv, d_Cv_dlnT, eta, d_eta_dlnT, ierr)
            use eos_lib, only: eosDT_get
            use eos_def
            integer, intent(in) :: eos_handle, species
            integer, pointer :: chem_id(:) ! maps species to chem id
            integer, pointer :: net_iso(:) ! maps chem id to species number
            real(dp), intent(in) :: &
               xa(:), rho, logRho, T, logT
            real(dp), intent(out) :: &
               Cv, d_Cv_dlnT, eta, d_eta_dlnT
            integer, intent(out) :: ierr

            real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
            real(dp) :: d_dxa(num_eos_d_dxa_results,species)

            include 'formats'
            ierr = 0
            
            call eosDT_get( &
               eos_handle, species, chem_id, net_iso, xa, &
               Rho, logRho, T, logT, &
               res, d_dlnd, d_dlnT, d_dxa, ierr)

            if (ierr /= 0) then
               write(*,*) 'failed in eosDT_get'
               return
            end if

            Cv = res(i_cv)
            d_Cv_dlnT = d_dlnT(i_cv)

            eta = res(i_eta)
            d_eta_dlnT = d_dlnT(i_eta)
         
         end subroutine get_eos_info_for_burn_at_const_density


         subroutine burn_finish_substep(nstp, time, y, ierr)
            use chem_def, only: category_name
            integer,intent(in) :: nstp
            real(dp), intent(in) :: time, y(:)
            integer, intent(out) :: ierr
            real(dp) :: frac, step_time
            integer :: j, i
            include 'formats'
            ierr = 0
            !step_time = time - substep_start_time
            !if (step_time <= 0d0) return
            !frac = step_time/dt
            !do j = 1, num_categories
            !   s% eps_nuc_categories(j,k) = &
            !      s% eps_nuc_categories(j,k) + frac*eps_nuc_cat(j)
            !end do
            !if (.false. .and. k == s% nz) then
            !   i = maxloc(eps_nuc_cat(1:num_categories),dim=1)
            !   write(*,3) 'frac time/dt eps_nuc_cat ' // trim(category_name(i)), &
            !      i, k, frac, time/dt, eps_nuc_cat(i), s% eps_nuc_categories(i,k)
            !end if
            !substep_start_time = time
         end subroutine burn_finish_substep

      end subroutine burn1_zone


      end module struct_burn_mix


