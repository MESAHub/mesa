! ***********************************************************************
!
!   Copyright (C) 2013-2019  Bill Paxton & The MESA Team
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

      module solve_burn
      use utils_lib, only: is_bad
      use star_private_def
      use const_def

      implicit none

      private
      public :: do_burn


      contains


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
         logical :: okay_to_reuse_rate_screened, use_pivoting, trace, burn_dbg
         
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
         okay_to_reuse_rate_screened = .false.
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
               okay_to_reuse_rate_screened, &
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
         call do1_net(s, k, s% species, .false., s% num_reactions, &
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


      end module solve_burn

