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

         if ((s% eps_nuc_factor == 0d0 .and. s% dxdt_nuc_factor == 0d0) .or. &
               s% gamma_law_hydro > 0d0) then
            do k = 1, nz
               s% eps_nuc(k) = 0d0
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
            if (k == s% trace_k) &
               write(*,1) 'call do1_burn from do_burn loop'
            if (.false.) then
               call burn1_BE( &
                  s, k, species, net_lwork, dt, &
                  num_iters, avg_epsnuc, burn_dbg, op_err)
            else
               call burn1_zone( &
                  s, k, species, min_T_for_const_density_solver, skip_burn, &
                  net_lwork, burn_lwork, screening_mode, &
                  dt, num_iters, avg_epsnuc, burn_dbg, op_err)
            end if
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
            Cv0 = s% dE_dT(k)
            eta0 = s% eta(k)
            call net_1_zone_burn_const_density( &
               s% net_handle, species, s% num_reactions, 0d0, dt, &
               xa_start, starting_log10T, s% lnd(k)/ln10, &
               get_eos_info_for_burn_at_const_density, &
               s% rate_factors, s% weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, &
               screening_mode, s% theta_e(k),  &
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
            s% xh(s% i_lnT,k) = starting_log10T*ln10
            s% lnT(k) = s% xh(s% i_lnT,k)
            s% T(k) = exp(s% lnT(k))
         else
            log10Ts_f1 => log10Ts_ary
            log10Rhos_f1 => log10Rhos_ary
            etas_f1 => etas_ary
            nullify(dxdt_source_term, times)
            log10Ts_f1(1) = s% lnT(k)/ln10
            log10Rhos_f1(1) = s% lnd(k)/ln10
            etas_f1(1) = s% eta(k)
            call net_1_zone_burn( &
               s% net_handle, species, s% num_reactions, 0d0, dt, xa_start, &
               num_times, times, log10Ts_f1, log10Rhos_f1, etas_f1, dxdt_source_term, &
               s% rate_factors, s% weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, &
               screening_mode, s% theta_e(k),  &
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
               z, xh, abar, zbar, xa, rho, logRho, T, logT, &
               Cv, d_Cv_dlnT, eta, d_eta_dlnT, ierr)
            use eos_lib, only: eos_get_helm_results
            use eos_def
            real(dp), intent(in) :: &
               z, xh, abar, zbar, xa(:), rho, logRho, T, logT
            real(dp), intent(out) :: &
               Cv, d_Cv_dlnT, eta, d_eta_dlnT
            integer, intent(out) :: ierr            
            real(dp) :: res(num_helm_results)
            real(dp), parameter :: coulomb_temp_cut = 1d6, coulomb_den_cut = 1d3
            include 'formats'            
            ierr = 0
            
            Cv = Cv0
            d_Cv_dlnT = 0d0
            
            eta = eta0
            d_eta_dlnT = 0d0
            
            !return
            
            call eos_get_helm_results( &
               xh, abar, zbar, rho, logRho, T, logT, &
               coulomb_temp_cut, coulomb_den_cut, &
               .true., .false., .false., &
               s% eos_rq% logT_ion_HELM, s% eos_rq% logT_neutral_HELM, &
               res, ierr)
            if (ierr /= 0) then
               !write(*,2) 'burn1_zone failed in eos_get_helm_results', k
               return
            end if            
            Cv = res(h_cv)
            d_Cv_dlnT = res(h_dcvdt)*T
            eta = res(h_etaele)
            d_eta_dlnT = res(h_detat)*T
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



      subroutine do_solve_radiative_diffusion(s, dt_total, kmin, kmax, ierr)
         use mesh_adjust, only: set_lnT_for_energy
         use chem_def, only: ih1, ihe3, ihe4
         use hydro_vars, only: set_vars
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt_total
         integer, intent(in) :: kmin, kmax
         integer, intent(out) :: ierr

         integer :: n, i, j, k, max_iters_per_substep, &
            max_iters_total, total_num_iters, num_iters, op_err, kbad
         integer :: steps_used, max_steps, min_steps
         real(qp) :: remaining_time, total_time, time, dt, dt1, dt_limit, dE, &
            max_del, avg_del, tol_correction_max, tol_correction_norm, T, crad_qp
         real(dp) :: total, dt_dble, alfa, kap_face, area, new_energy, &
            new_T, new_lnT, revised_energy, dxdt, dm
         real(dp), pointer, dimension(:) :: sig
         real(qp), pointer, dimension(:) :: &
            du, d, dl, x, b, bp, vp, xp, dX, X_0, X_1, rhs, del, erad_start
         logical :: okay
         real(qp), parameter :: xl0 = 0, xl1 = 1
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         if (dbg) write(*,3) 'do_solve_radiative_diffusion', kmin, kmax

         total_time = dt_total
         time = 0
         max_steps = 20
         min_steps = 4
         max_iters_per_substep = 4
         max_iters_total = 40
         total_num_iters = 0
         tol_correction_max = 1d-4
         tol_correction_norm = 1d-7

         if (dt_total <= 0d0) return
         
         if (kmin <= 1) then
            if (s% report_ierr) write(*,*) 'kmin must be > 1 for do_solve_radiative_diffusion'
            ierr = -1
            return
         end if
         
         if (kmax > s% nz) then
            if (s% report_ierr) write(*,*) 'kmax must be <= nz for do_solve_radiative_diffusion'
            ierr = -1
            return
         end if
         
         n = kmax - kmin + 1

         call do_alloc(ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'allocate failed in do_solve_radiative_diffusion'
            return
         end if
         
         crad_qp = crad
         do i = 1, n
            k = i + kmin - 1
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            kap_face = alfa*s% opacity(k) + (1d0-alfa)*s% opacity(k-1)
            area = 4d0*pi*s% r(k)*s% r(k)
            sig(i) = area*area*clight/(3d0*kap_face*s% dm_bar(k))
            T = s% T(k)
            erad_start(i) = crad_qp*T*T*T*T ! erad(k)
            X_1(i) = erad_start(i)
         end do

         okay = .true.
         steps_used = 0
         dt_limit = 1d99
         do i = 1, n
            k = i + kmin - 1
            dt1 = 1d20*s% dm(k)/sig(i)
            if (dt1 < dt_limit) dt_limit = dt1
         end do
         dt = 0.5d0*dt_limit
         if (dbg) write(*,1) 'dt/dt_limit', dt/dt_limit, dt, dt_limit

      step_loop: do while &
               (total_time - time > 1d-10*total_time .and. &
                  steps_used < max_steps)

            steps_used = steps_used + 1
            remaining_time = total_time - time
            dt = max(dt, 1d-6*remaining_time)
            dt = min(dt, real(dt_total/min_steps,kind=qp))
            if (dt >= remaining_time) then
               dt = remaining_time
            else
               dt = min(dt, 0.5d0*remaining_time)
            end if
            if (steps_used >= max_steps) dt = remaining_time ! just go for it
            if (dbg) write(*,2) 'step dt', steps_used, dt, dt/remaining_time

            ! setup for new step
            do i = 1, n
               X_0(i) = X_1(i)
               dX(i) = 0d0
            end do

         solve_loop: do num_iters = 1, max_iters_per_substep

               if (total_num_iters >= max_iters_total) then
                  ierr = -1
                  if (s% report_ierr) write(*,*) 'failed to converge in allowed number of steps'
                  exit step_loop
               end if

               total_num_iters = total_num_iters+1
               
               do i = 1, n
                  k = i + kmin - 1
                  dm = s% dm(k)
                  if (i == 1) then
                     dxdt = sig(i+1)*(X_1(i+1) - X_1(i))/dm
                  else if (i == n) then
                     dxdt = -sig(i)*(X_1(i) - X_1(i-1))/dm
                  else
                     dxdt = (sig(i+1)*(X_1(i+1) - X_1(i)) - sig(i)*(X_1(i) - X_1(i-1)))/dm
                  end if
                  rhs(i) = dt*dxdt - dX(i)
                  if (i > 1) dl(i-1) = -dt*sig(i)/dm
                  if (i < n) then
                     du(i) = -dt*sig(i+1)/dm
                     d(i) = xl1 + dt*(sig(i) + sig(i+1))/dm
                  else
                     du(i) = xl0
                     d(i) = xl1 + dt*sig(i)/dm
                  end if
                  if (is_bad(rhs(i) + d(i))) then
                     if (s% report_ierr) write(*,2) 'rhs d', i, rhs(i), d(i)
                     return
                     stop 'solve_radiative_diffusion'
                  end if
               end do
               dl(n) = xl0

               call solve_tridiag(dl, d, du, rhs(1:n), del(1:n), n, ierr)
               if (ierr /= 0) then
                  if (s% report_ierr) write(*,*) 'matrix solve failed in solve mix'
                  exit step_loop
               end if

               ! apply the correction dX = dX + del
               ! X_1 = X_0 + dX
               ! X_0 is erad at start of substep
               ! X_1 is candidate for erad at end of substep
               do i = 1, n
                  dX(i) = dX(i) + del(i)
                  X_1(i) = X_0(i) + dX(i)
                  if (is_bad(X_1(i))) then
                     if (s% report_ierr) write(*,2) 'new X_1(i)', i, X_1(i), del(i)
                     return
                     stop 'solve_radiative_diffusion'
                  end if
               end do

               ! if correction small enough, exit solve_loop
               max_del = maxval(abs(del(1:n)))
               avg_del = sum(abs(del(1:n)))/n
               if (max_del <= tol_correction_max .and. avg_del <= tol_correction_norm) then
                  if (dbg) &
                     write(*,2) 'substep converged: iters max_del avg_del dt/total', &
                        num_iters, max_del, avg_del, dt/total_time
                  exit solve_loop ! this substep is done
               end if

               if (dbg) &
                  write(*,2) 'substep not converged: iters max_del avg_del', &
                     num_iters, max_del, avg_del

               if (num_iters == max_iters_per_substep) then
                  if (s% report_ierr) write(*,*) 'num_iters hit max_iters_per_substep'
                  exit step_loop
               end if

            end do solve_loop

            time = time + dt

         end do step_loop

         if (total_time - time > 1d-10*total_time) then
            ierr = -1
            if (s% report_ierr) write(*,*) 'failed in solve_radiative_diffusion'
         end if

         if (dbg) write(*,*)
         
         if (ierr == 0) then ! update lnT for new energy
            kbad = 0
!x$OMP PARALLEL DO PRIVATE(i,k,dE,new_energy,new_T,new_lnT,revised_energy,op_err) SCHEDULE(dynamic,2)
            do i = 1, n
               if (ierr /= 0) cycle
               k = i + kmin - 1
               dE = X_1(i) - erad_start(i)
               if (abs(dE) < s% energy(k)*1d-8) then
                  new_T = s% T(k) + dE/s% dE_dT(k)
                  new_lnT = log(new_T)
                  if (dbg) write(*,2) 'dlogT new old dE/E', k, &
                     (new_lnT - s% lnT(k))/ln10, new_lnT/ln10, s% lnT(k)/ln10, dE/s% energy(k)
                  s% xh(s% i_lnT,k) = new_lnT
                  s% lnT(k) = new_lnT
                  s% T(k) = new_T
               else
                  new_energy = s% energy(k) + dE
                  op_err = 0
                  call set_lnT_for_energy( &
                     s, k, &
                     s% net_iso(ih1), s% net_iso(ihe3), s% net_iso(ihe4), &
                     s% species, s% xa(:,k), &
                     s% rho(k), s% lnd(k)/ln10, new_energy, s% lnT(k), &
                     new_lnT, revised_energy, op_err)
                  if (op_err /= 0) then
                     kbad = k; ierr = -1; cycle
                  end if
                  if (dbg) write(*,2) 'dlogT new old dE/E', k, &
                     (new_lnT - s% lnT(k))/ln10, new_lnT/ln10, s% lnT(k)/ln10, dE/s% energy(k)
                  s% xh(s% i_lnT,k) = new_lnT
                  s% lnT(k) = new_lnT
                  s% T(k) = exp(new_lnT)
               end if
            end do
!x$OMP END PARALLEL DO
         end if

         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,2) 'set_lnT_for_energy failed in call solve_radiative_diffusion', kbad
         end if
         
         if (ierr == 0) then
            call set_vars(s, dt_total, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) &
                  write(*,2) 'set_vars failed in call solve_radiative_diffusion', kbad
            end if
         end if
         
         if (dbg) stop 'solve_radiative_diffusion'

         call dealloc


         contains


         subroutine do_alloc(ierr)
            use alloc
            integer, intent(out) :: ierr

            call non_crit_get_work_array(s, sig, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return

            call non_crit_get_quad_array(s, erad_start, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, du, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, d, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, dl, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, x, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, b, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, bp, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, vp, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, xp, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return

            call non_crit_get_quad_array(s, dX, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, X_0, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, X_1, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, rhs, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return
            call non_crit_get_quad_array(s, del, n, 0, 'solve_radiative_diffusion', ierr)
            if (ierr /= 0) return


         end subroutine do_alloc


         subroutine dealloc
            use alloc
            call non_crit_return_work_array(s, sig, 'solve_radiative_diffusion')
            
            call non_crit_return_quad_array(s, erad_start, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, du, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, d, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, dl, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, x, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, b, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, bp, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, vp, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, xp, 'solve_radiative_diffusion')

            call non_crit_return_quad_array(s, dX, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, X_0, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, X_1, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, rhs, 'solve_radiative_diffusion')
            call non_crit_return_quad_array(s, del, 'solve_radiative_diffusion')
         end subroutine dealloc


         subroutine solve_tridiag(sub, diag, sup, rhs, x, n, ierr)
            !      sub - sub-diagonal
            !      diag - the main diagonal
            !      sup - sup-diagonal
            !      rhs - right hand side
            !      x - the answer
            !      n - number of equations
            integer, intent(in) :: n
            real(qp), dimension(:), intent(in) :: sup, diag, sub
            real(qp), dimension(:), intent(in) :: rhs
            real(qp), dimension(:), intent(out) :: x
            integer, intent(out) :: ierr

            real(qp) :: m
            integer i

            ierr = 0

            bp(1) = diag(1)
            vp(1) = rhs(1)

            do i = 2,n
               m = sub(i-1)/bp(i-1)
               bp(i) = diag(i) - m*sup(i-1)
               vp(i) = rhs(i) - m*vp(i-1)
            end do

            xp(n) = vp(n)/bp(n)
            x(n) = xp(n)
            do i = n-1, 1, -1
               xp(i) = (vp(i) - sup(i)*xp(i+1))/bp(i)
               x(i) = xp(i)
            end do

         end subroutine solve_tridiag


      end subroutine do_solve_radiative_diffusion
   

      ! "BE" = Backward Euler.    no substeps in current version.
      subroutine burn1_BE( &
            s, k, species, net_lwork, dt_total, &
            num_iters_out, avg_epsnuc, dbg_in, ierr)

         use chem_def, only: chem_isos, ih1, ihe3, ihe4
         use chem_lib, only: basic_composition_info
         use net, only: do1_net
         use net_lib, only: get_reaction_id_table_ptr
         use net_def, only: Net_Info
         use eos_lib, only: eos_get_helm_results
         use eos_def
         use utils_lib, only: is_bad
         use rates_def, only: num_rvs, reaction_Name
         use mesh_adjust, only: set_lnT_for_energy
         use micro, only: do_eos_for_cell

         type (star_info), pointer :: s
         integer, intent(in) :: k, species, net_lwork
         real(dp), intent(in) :: dt_total
         logical, intent(in) :: dbg_in
         integer, intent(out) :: num_iters_out
         real(dp), intent(out) :: avg_epsnuc
         integer, intent(out) :: ierr

         logical :: did_refactor, reuse_given_rates, dbg, unchanged

         integer, target :: ipivot_array(species+1)
         integer, pointer :: ipivot(:)
         real(dp), target, dimension((species+1)*(species+1)) :: mtx_array
         real(dp), pointer, dimension(:,:) :: mtx
         real(dp), target, dimension(species+1) :: &
            del_array, x0_array, x1_array, dx_array
         real(dp), pointer, dimension(:) :: del, x0, x1, dx

         real(dp), parameter :: one = 1, zero = 0
         real(dp) :: &
            tol_max_corr, tol_avg_corr, max_resid, total_burn, time, dt, remaining_time, &
            lambda, lambda0, min_lambda, Cv0, lnT0, dCv_dlnT, Cv, T, prev_max_correction, &
            new_energy, revised_energy, new_lnT

         real(dp) :: sumx

         integer, pointer :: reaction_id(:) ! maps net reaction number to reaction id
         integer, pointer :: chem_id(:)
         integer :: i, j, nz, max_iters, nvar, num_substeps, num_iters, i_lnT, &
            max_steps, steps_used, num_tries, max_tries

         real(dp), dimension(species) :: d_dxdt_dRho, d_dxdt_dT
         real(dp) :: d_dxdt_nuc_dx(species,species)

         include 'formats'

         dbg = (k == 1008) ! dbg_in

         ierr = 0
         avg_epsnuc = 0
         nvar = species+1
         chem_id => s% chem_id
         i_lnT = s% i_lnT
         
         do j=1,species
            s% xa_start(j,k) = s% xa(j,k)
         end do

         call basic_composition_info( &
            species, chem_id, s% xa(1:species,k), s% X(k), s% Y(k), s% Z(k), &
            s% abar(k), s% zbar(k), s% z2bar(k), s% z53bar(k), s% ye(k), &
            s% mass_correction(k), sumx)

         reuse_given_rates = .false.
         call do1_net(s, k, species, &
            reuse_given_rates, s% num_reactions, net_lwork, .false., ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'burn1_BE: initial do1_net failed', k
            return
         end if

         if (dt <= 0d0) return

         reuse_given_rates = .false. ! s% burn_reuse_given_rates

         ipivot => ipivot_array
         del => del_array
         x0 => x0_array
         x1 => x1_array
         dx => dx_array
         mtx(1:nvar,1:nvar) => mtx_array(1:nvar*nvar)

         min_lambda = 0.1d0 ! s% op_split_burn_min_lambda

         tol_avg_corr = 1d-4 ! s% op_split_burn_tol_avg_correction
         tol_max_corr = 1d-5 ! s% op_split_burn_tol_max_correction
         
         Cv0 = s% dE_dT(k)
         lnT0 = s% lnT(k)
         dCv_dlnT = s% d_eos_dlnT(i_Cv,k)
         Cv = Cv0
         T = s% T(k)

         do j=1,species
            x0(j) = s% xa(j,k)
            x1(j) = x0(j)
            dx(j) = 0d0
         end do

         x0(nvar) = lnT0
         x1(nvar) = lnT0
         dx(nvar) = 0d0

         total_burn = 0d0
         time = 0
         max_steps = 100
         steps_used = 0
         num_iters = 0
         max_tries = 10
         num_iters = 0
         max_iters = 50 ! s% op_split_burn_max_iterations

         dt = dt_total*1d-5
         lambda = 1d0
         unchanged = .true.

         if (dbg) write(*,*)
         
         step_loop: do while &
               (dt_total - time > 1d-10*dt_total .and. &
                  steps_used < max_steps)
         
            steps_used = steps_used + 1
            
            retry_loop: do num_tries = 1, max_tries
            
               remaining_time = dt_total - time
               if (num_tries == 1) then
                  if (lambda == 1d0 .and. num_iters < max_iters/2) dt = 1.5d0*dt
                  if (dt >= remaining_time) then
                     dt = remaining_time
                  else
                     dt = min(dt, 0.5d0*remaining_time)
                  end if
                  if (steps_used >= max_steps) dt = remaining_time ! just go for it
               end if
               if (dbg) write(*,4) 'burn_BE k, step, try, dt, dt/remaining', &
                     k, steps_used, num_tries, dt, dt/remaining_time
            
               prev_max_correction = 1d99
            
               solve_loop: do num_iters = 1, max_iters

                  do j=1,species
                     del(j) = dt*s% dxdt_nuc(j,k) - dx(j) ! residual
                     do i=1,species
                        mtx(j,i) = -dt*s% d_dxdt_nuc_dx(j,i,k)
                     end do
                     mtx(j,j) = one + mtx(j,j)
                     mtx(nvar,j) = -dt*s% d_epsnuc_dx(j,k)/(Cv*T)
                     mtx(j,nvar) = -dt*s% d_dxdt_nuc_dT(j,k)*T
                  end do
                  del(nvar) = dt*s% eps_nuc(k)/(Cv*T) - dx(nvar)
                  mtx(nvar,nvar) = one + dt/(Cv*Cv*T)* &
                     (s% eps_nuc(k)*dCv_dlnT + Cv*(s% eps_nuc(k) - s% d_epsnuc_dlnT(k)))

                  ! solve mtx*del = rhs (rhs is in del at start)
                  call solve_dp(s, k, nvar, mtx, ipivot, del, nvar, ierr)
                  if (ierr /= 0) then
                     write(*,3) 'solve_dp failed', k, s% model_number
                     call dealloc
                     return
                  end if

                  ! set lambda for positivity
                  call positivity(s, k, species, min_lambda, x1, del, lambda, 'BE', num_iters, dbg)

                  lambda0 = lambda
                  if (lambda < min_lambda) then
                     lambda = min_lambda
                  end if

                  do j=1,species ! apply the correction for abundances
                     dx(j) = dx(j) + del(j)*lambda
                     x1(j) = x0(j) + dx(j)
                     s% xa(j,k) = x1(j)
                     if (abs(s% xa(j,k) - s% xa_start(j,k)) > 1d-14) unchanged = .false.
                  end do
            
                  dx(nvar) = dx(nvar) + del(nvar)*lambda
                  x1(nvar) = x0(nvar) + dx(nvar)
                  s% lnT(k) = x1(nvar)
                  s% xh(i_lnT,k) = s% lnT(k)
                  s% T(k) = exp(s% lnT(k))
                  T = s% T(k)

                  if (.not. unchanged) then
                     
                     call do_eos_for_cell(s, k, ierr)
                     if (ierr /= 0) then
                        if (s% report_ierr) write(*,2) 'burn1_BE failed in do_eos_for_cell', k
                        exit step_loop
                     end if
                     Cv = s% dE_dT(k)
                     dCv_dlnT = s% d_eos_dlnT(i_Cv,k)

                     call do1_net(s, k, species, &
                        reuse_given_rates, s% num_reactions, net_lwork, .false., ierr)
                     if (ierr /= 0) then
                        if (s% report_ierr) write(*,2) 'burn1_BE failed in do1_net', k
                        exit step_loop
                     end if

                  end if

                  s% max_burn_correction(k) = maxval(abs(del(1:nvar)))
                  s% avg_burn_correction(k) = sum(abs(del(1:nvar)))/nvar
               
                  if (dbg) write(*,2) 'iter, max corr, avg corr, lambda', num_iters, &
                     s% max_burn_correction(k), s% avg_burn_correction(k), lambda
                  
                  if (s% max_burn_correction(k) > prev_max_correction) then
                     if (dbg) write(*,1) 'increase in max correction: retry', &
                        s% max_burn_correction(k), prev_max_correction
                     exit solve_loop
                  end if
                  prev_max_correction = s% max_burn_correction(k)

                  if (lambda0 == one) then

                     ! if magnitude of max correction is small enough, consider converged.
                     if (is_bad(s% max_burn_correction(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'max_burn_correction', k, s% max_burn_correction(k)
                        if (.not. (dbg .or. dbg_in .or. s% report_ierr)) return

                        do j=1,species
                           write(*,2) 'del', j, del(j)
                        end do
                        write(*,3) 'max_burn_correction', k, num_iters, s% max_burn_correction(k)
                        stop 'burn1_BE'
                     end if

                     if (s% max_burn_correction(k) <= tol_max_corr .and. &
                         s% avg_burn_correction(k) <= tol_avg_corr) exit retry_loop

                  end if

                  if (num_iters == max_iters) then
                     if (dbg) write(*,2) 'no more iters allowed', max_iters
                     exit solve_loop
                  end if

               end do solve_loop ! do another iter
               
               if (num_tries == max_tries) then
                  if (dbg) write(*,2) 'no more tries allowed', max_tries
                  exit step_loop
               end if
               
               if (dbg) write(*,*) 'retry with smaller dt'
               do j=1,species
                  x1(j) = x0(j)
                  s% xa(j,k) = x1(j)
                  dx(j) = 0d0
               end do
               s% lnT(k) = x0(nvar)
               s% xh(i_lnT,k) = s% lnT(k)
               s% T(k) = exp(s% lnT(k))
               T = s% T(k)

               call do_eos_for_cell(s, k, ierr)
               if (ierr /= 0) then
                  if (s% report_ierr) write(*,2) 'retry burn1_BE failed in do_eos_for_cell', k
                  exit step_loop
               end if
               Cv = s% dE_dT(k)
               dCv_dlnT = s% d_eos_dlnT(i_Cv,k)

               call do1_net(s, k, species, &
                  reuse_given_rates, s% num_reactions, net_lwork, .false., ierr)
               if (ierr /= 0) then
                  if (dbg .or. s% report_ierr) write(*,2) 'do1_net failed in retry for BE_burn', k
                  exit step_loop
               end if

               dt = 0.5d0*dt

            end do retry_loop
            
            if (steps_used >= max_steps) exit step_loop

            do j=1,nvar
               x0(j) = x1(j)
               dx(j) = 0d0
            end do
            
            time = time + dt
            total_burn = total_burn + dt*s% eps_nuc(k)

            if (dbg) write(*,*)
         
         end do step_loop
         
         if (dt_total - time > 1d-10*dt_total) then
            if (dbg .or. s% report_ierr) write(*,2) 'burn1_BE ran out of steps', k
            ierr = -1
         end if
         
         avg_epsnuc = total_burn/dt_total
         num_iters_out = steps_used

         if (ierr == 0) then
            new_energy = s% energy(k) + total_burn
            call set_lnT_for_energy( &
               s, k, &
               s% net_iso(ih1), s% net_iso(ihe3), s% net_iso(ihe4), &
               s% species, s% xa(:,k), &
               s% rho(k), s% lnd(k)/ln10, new_energy, s% lnT(k), &
               new_lnT, revised_energy, ierr)
            if (ierr /= 0) then
               if (dbg .or. s% report_ierr) write(*,2) 'burn1_BE failed in call set_lnT_for_energy', k
            else
               if (dbg) write(*,2) 'dlogT new old dE/E', k, &
                  (new_lnT - s% lnT(k))/ln10, new_lnT/ln10, s% lnT(k)/ln10, total_burn/s% energy(k)
               s% xh(s% i_lnT,k) = new_lnT
               s% lnT(k) = new_lnT
               s% T(k) = exp(new_lnT)
            end if
         end if
         
         if (dbg) write(*,3) 'burn1_BE k steps_used avg_epsnuc', k, steps_used, avg_epsnuc

         call dealloc
         
         !if (dbg) stop 'burn1_BE'


         contains


         subroutine dealloc
         end subroutine dealloc


         subroutine report_failure(k, num_iters)
            integer, intent(in) :: k, num_iters
            include 'formats'
            if (.not. s% report_ierr) return
            write(*,4) 'BE failed to converge: k, num_iters, model, T', &
                  k, num_iters, s% model_number, s% T(k)
         end subroutine report_failure

         subroutine show_stuff(str)
            character (len=*), intent(in) :: str
            include 'formats'
            write(*,4) trim(str) // ' T, logT', &
               num_iters, k, s% model_number, s% T(k), s% lnT(k)/ln10
            write(*,4) trim(str) // ' rho, logRho', &
               num_iters, k, s% model_number, s% rho(k), s% lnd(k)/ln10
            write(*,4) trim(str) // ' Ye, eta', &
               num_iters, k, s% model_number, s% Ye(k), s% eta(k)
            write(*,4) trim(str) // ' dt', &
               num_iters, k, s% model_number, dt
            write(*,4) trim(str) // ' sum(x0)', &
               num_iters, k, s% model_number, sum(x0(1:species))
            write(*,4) trim(str) // ' sum(del)', &
               num_iters, k, s% model_number, sum(del(1:species))
            write(*,4) trim(str) // ' sum(xa)', &
               num_iters, k, s% model_number, sum(s% xa(1:species,k))

            write(*,*)
            write(*,2) 'del'
            do j=1,species
               !if (s% xa(j,k) > 1d-16) &
                  write(*,3) 'del ' // trim(chem_isos% name(s% chem_id(j))), &
                     j, k, del(j)
            end do

            write(*,*)
            write(*,2) 'xa'
            do j=1,species
               !if (s% xa(j,k) > 1d-16) &
                  write(*,3) 'xa ' // trim(chem_isos% name(s% chem_id(j))), &
                     j, k, s% xa(j,k)
            end do


            write(*,*)
            write(*,*)
            write(*,2) 'dt', k, dt
            write(*,2) 'sum dxdt_nuc', k, sum(s% dxdt_nuc(:,k))
            write(*,2) 'dt*(sum dxdt_nuc)', k, dt*sum(s% dxdt_nuc(:,k))
            write(*,2) 'sum del', k, sum(del(:))
            write(*,2) 'sum dx', k, sum(dx(:))
            write(*,2) 'sum xa', k, sum(s% xa(:,k))
            write(*,*)


            write(*,1) 'large magnitude entries'
            write(*,*)
            write(*,2) 'abs(dt*dxdt_nuc) > 1d-4', k
            do j=1,species
               if (dt*abs(s% dxdt_nuc(j,k)) > 1d-4) &
                  write(*,3) 'dt*dxdt_nuc(j,k) ' // trim(chem_isos% name(s% chem_id(j))), &
                     j, k, dt*s% dxdt_nuc(j,k)
            end do
            write(*,*)
            write(*,2) 'abs(dt*d_dxdt_nuc_dx) > 1d12', k
            do j=1,species
               do i=1,species
                  if (dt*abs(s% d_dxdt_nuc_dx(i,j,k)) > 1d12) &
                     write(*,4) 'dt*d_dxdt_nuc_dx(i,j,k) ' // trim(chem_isos% name(s% chem_id(j))) &
                        // ' wrt ' // trim(chem_isos% name(s% chem_id(i))), &
                        i, j, k, dt*s% d_dxdt_nuc_dx(i,j,k)
               end do
            end do
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,1) 'nonzero entries'
            write(*,*)
            do j=1,species
               if (s% dxdt_nuc(j,k) /= 0d0) &
                  write(*,3) 'dxdt_nuc(j,k) ' // trim(chem_isos% name(s% chem_id(j))), &
                     j, k, s% dxdt_nuc(j,k)
            end do
            write(*,*)
            write(*,2) 'd_dxdt_nuc_dx', k
            do j=1,species
               do i=1,species
                  if (s% d_dxdt_nuc_dx(i,j,k) /= 0d0) &
                     write(*,4) 'd_dxdt_nuc_dx(i,j,k) ' // trim(chem_isos% name(s% chem_id(j))) &
                        // ' wrt ' // trim(chem_isos% name(s% chem_id(i))), &
                        i, j, k, s% d_dxdt_nuc_dx(i,j,k)
               end do
            end do
            write(*,2) 'done d_dxdt_nuc_dx', k
            write(*,*)
            do j=1,species
               write(*,3) 'sum d_dxdt_nuc_dx(:,j,k) ' // trim(chem_isos% name(s% chem_id(j))), &
                  j, k, sum(s% d_dxdt_nuc_dx(:,j,k))
            end do
            write(*,*)
            do i=1,species
               write(*,3) 'sum d_dxdt_nuc_dx(i,:,k) ' // trim(chem_isos% name(s% chem_id(i))), &
                  j, k, sum(s% d_dxdt_nuc_dx(i,:,k))
            end do

            !stop 'burn1_BE'
         end subroutine show_stuff


      end subroutine burn1_BE


      subroutine positivity( &
            s, k, species, min_lambda, x, del, lambda, str, iter, dbg)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: k, species, iter
         logical, intent(in) :: dbg
         real(dp), intent(in) :: min_lambda
         real(dp), intent(in), dimension(:) :: x, del
         character(len=*), intent(in) :: str
         real(dp), intent(out) :: lambda

         integer :: i, j, bad_j
         real(dp) :: alpha, new_xa, old_xa, dxa, eps

         include 'formats'

         lambda = 1d0
         eps = 1d-6 ! s% op_split_burn_tol_max_correction ! allow this amount below 0
         bad_j = 0
         do j=1,species
            old_xa = x(j)
            if (old_xa < 1d-99) cycle
            dxa = del(j)
            new_xa = old_xa + dxa
            if (new_xa >= 0) cycle
            alpha = -(old_xa + eps)/dxa
            ! so dxa*alpha = -old_xa - eps,
            ! and therefore old_xa + alpha*dxa = -eps
            if (alpha < lambda) then
               lambda = alpha
               bad_j = j
            end if
         end do
         if (lambda < min_lambda) lambda = min_lambda
         if (dbg .and. lambda < 1) then
            j = bad_j
            write(*,3) 'lambda ' // trim(chem_isos% name(s% chem_id(j))), &
               iter, j, lambda, x(j) + del(j), x(j) + lambda*del(j), x(j), del(j)
         end if

      end subroutine positivity


      subroutine solve_dp(s, k, nvar, mtx, ipivot, del, species, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar, species
         integer :: ipivot(:)
         real(dp) :: mtx(:,:), del(:)
         integer, intent(out) :: ierr
         include 'formats'
         call my_getf2(nvar, mtx, nvar, ipivot, ierr)
         if (ierr /= 0) return
         call my_getrs1(nvar, mtx, nvar, ipivot, del, species, ierr)
         if (ierr /= 0) return
         ! solution is now in del
      end subroutine solve_dp


      include "mtx_solve_routines.inc"


      end module solve_burn

