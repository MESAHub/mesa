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
         use solve_hydro, only: &
            set_surf_info, set_tol_correction, do_solver_converge
         use solve_burn, only: do_burn
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
         use solve_hydro, only: set_luminosity_by_category
         use chem_def, only: num_categories
         use hydro_tdc, only: set_etrb_start_vars
         use star_utils, only: eval_total_energy_integrals
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
            s% eps_nuc_start(k) = s% eps_nuc(k)
            s% opacity_start(k) = s% opacity(k)
            s% mlt_mixing_length_start(k) = s% mlt_mixing_length(k)
            s% mlt_mixing_type_start(k) = s% mlt_mixing_type(k)
            s% mlt_D_start(k) = s% mlt_D(k)
            s% mlt_vc_start(k) = s% mlt_vc(k)
            s% mlt_Gamma_start(k) = s% mlt_Gamma(k)
         end do
         
         if (s% TDC_flag) then
            call set_etrb_start_vars(s,ierr)
            !do k=1,s% nz ! DEBUGGING INFO
            !   s% xtra1_array(k) = s% rho(k)
            !   s% xtra2_array(k) = s% T(k)
            !   s% xtra3_array(k) = abs(s% w(k)) + 1d0
            !   s% xtra4_array(k) = abs(s% v(k)) + 1d0
            !   s% xtra5_array(k) = s% r(k)
            !end do
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


      end module struct_burn_mix


