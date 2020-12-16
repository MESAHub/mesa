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
!   MERCHANT\ILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
!
! ***********************************************************************

      module diffusion

      use const_def
      use chem_def
      use star_private_def
      use diffusion_support
      use diffusion_procs

      implicit none

      integer, parameter :: diffusion_min_nc = 4 ! minimum number of classes


      logical, parameter :: use_dcoeff_dX = .true.


      contains


      subroutine do_solve_diffusion( &
            s, nz, species, nc, m, class, class_chem_id, net_iso, chem_id, &
            abar, ye, free_e, dm_bar_in, dm_in, &
            T, lnT, rho, lnd, r_mid, dlnP_dm, dlnT_dm, dlnRho_dm, &
            L_face, r_face_in, dlnP_dm_face, dlnT_dm_face, dlnRho_dm_face, &
            pure_Coulomb, total_time, dt_div_timescale, &
            max_steps, max_iters_total, max_iters_per_substep, &
            calculate_ionization, typical_charge, nsmooth_typical_charge, &
            min_T_for_radaccel, max_T_for_radaccel, &
            min_Z_for_radaccel, max_Z_for_radaccel, &
            screening_for_radaccel, op_mono_data_path, op_mono_data_cache_filename, &
            v_advection_max, R_center, &
            gamma, gamma_full_on, gamma_full_off, &
            T_full_on, T_full_off, diffusion_factor, &
            xa, steps_used, total_num_iters, total_num_retries, &
            nzlo, nzhi, X_init, X_final, &
            D_self_face, v_advection_face, v_total_face, &
            vlnP_face, vlnT_face, v_rad_face, g_rad_face, &
            E_field_face, g_field_face, &
            ierr)

         use diffusion_procs, only: fixup
         use kap_lib, only: load_op_mono_data

         type (star_info), pointer :: s
         integer, intent(in) :: nz, species, nc, m, &
            min_Z_for_radaccel, max_Z_for_radaccel
            ! nc = number of classes of species
            ! m = nc+1
         integer, intent(in) :: class(:) ! (species)
         integer, intent(in) :: class_chem_id(:) ! (nc)
         integer, intent(in) :: net_iso(:), chem_id(:)
            ! class(i) = class number for species i. class numbers from 1 to nc
            ! class_chem_id(j) = isotope id number from chem_def for "typical" member of class j
         real(dp), intent(in), dimension(:) :: &
            abar, ye, free_e, gamma, dm_bar_in, dm_in, &
            T, lnT, rho, lnd, r_mid, dlnP_dm, dlnT_dm, dlnRho_dm, &
            L_face, r_face_in, dlnP_dm_face, dlnT_dm_face, dlnRho_dm_face
         logical, intent(in) :: pure_Coulomb, screening_for_radaccel
         character (len=*), intent(in) :: &
            op_mono_data_path, op_mono_data_cache_filename
         real(dp), intent(in) :: &
            total_time, dt_div_timescale, &
            min_T_for_radaccel, max_T_for_radaccel, &
            v_advection_max, R_center, gamma_full_on, gamma_full_off, &
            T_full_on, T_full_off, diffusion_factor(nc)
         integer, intent(in) :: max_steps, max_iters_total, max_iters_per_substep
         logical, intent(in) :: calculate_ionization
         real(dp), intent(inout) :: typical_charge(:,:) ! (nc,nz)
         integer, intent(in) :: nsmooth_typical_charge
         real(dp), intent(inout) :: xa(:,:) ! (species,nz) ! mass fractions
         real(dp), intent(out), dimension(:,:) :: X_init, X_final, &
            D_self_face, v_advection_face, v_total_face, &
            vlnP_face, vlnT_face, v_rad_face, g_rad_face ! (nc,nz)
         real(dp), intent(inout) :: E_field_face(:) ! (nz)
         real(dp), intent(inout) :: g_field_face(:) ! (nz)
         integer, intent(out) :: steps_used, total_num_iters, total_num_retries, ierr
         integer, intent(inout) ::  nzlo, nzhi !upper and lower bounds on region

         integer :: i, j, k, num_iters, idiag, n, neqs, ku, kl, &
            ldab, ldafb, ldb, ldx, kmax, h1, he4, nbound, bad_j, bad_k, &
            retry_count, max_retries, ierr_dealloc
         real(dp) :: &
            dt, dt_next, min_dt, time, mstar, mtotal, sum_mass_nzlo_nzhi, xtotal_init(nc), &
            xtotal(nc), mass_init(nc), frac, err, bad_X, bad_sum, bad_Xsum, &
            dx_max, dx_avg, tol_correction_max, tol_correction_norm, remaining_time, &
            sumx, tmp, timescale, min_lambda, upwind_limit, lambda, max_del, avg_del

         integer, parameter :: min_nz_lo = 5
         real(dp), parameter :: &
            dt_retry_factor = 0.5d0, dt_max_factor = 2d0, dt_min_factor = 0.75d0, &
            tiny_X = 1d-50, tiny_C = 1d-50, max_sum_abs = 10, X_cleanup_tol = 1d-2, &
            max_flow_frac = 1d0, max_flow_X_limit = 1d-5, dx_avg_target = 0.975d0

         real(dp), dimension(:), pointer :: cell_dm, dm_bar, A, sum_new_mass
         real(dp), dimension(:,:), pointer :: &
            Z, X, ending_dX_dm, C, dC_dr, C_div_X, &
            rhs, new_mass, del, X_0, X_1, dX, dX_dt
         real(dp), dimension(:,:,:), pointer :: em1, e00, ep1

         real(dp), dimension(:), pointer :: &
            r_face, alfa_face, rho_face, T_face, four_pi_r2_rho_face, &
            dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, sum_starting_mass, &
            limit_coeffs_face, &
            e_ap, e_at, e_ar, e_ax1, &
            g_ap, g_at, g_ar, g_ax1

         real(dp), dimension(:), pointer :: &
            Z1, X1, ending_dX_dm1, C1, dC_dr1, C_div_X1, new_mass1, &
            X_face1, C_div_X_face1, GT_face1, del1, X_01, X_11, dX1, dX_dt1, &
            ending_mass1, starting_mass1, starting_dX_dm1, SIG_face1, sigma_lnC1, &
            AD_face, rad_accel_face1, log10_g_rad1, &
            X_theta1, X_minus_theta1, xm_face, sum_ending_mass
         real(dp), dimension(:,:), pointer :: &
            X_face, C_div_X_face, GT_face, &
            rad_accel_face, log10_g_rad, &
            ending_mass, starting_mass, starting_dX_dm, &
            X_theta, X_minus_theta, e_ax, g_ax
         real(dp), dimension(:,:,:), pointer :: &
            SIG_face, sigma_lnC

         real(dp), pointer, dimension(:) :: rhs1, lblk1, dblk1, ublk1, b1
         integer, pointer :: ipiv1(:)
         integer :: lrd, lid, j_bad, k_bad, kmax_rad_accel, min_num_substeps, &
            iter_dbg, j_dbg, k_dbg, k_max
         integer(8) :: time0, time1, clock_rate
         integer, pointer :: ipar_decsol(:)
         real(dp), pointer :: rpar_decsol(:)
         real(dp), dimension(species) :: xa_total_before, xa_total_after
         real(dp), dimension(m) :: C_face, Z_face
         real(dp) :: X_total_atol, X_total_rtol, b_bad, flow_out, flow_in, &
            vc_target, vc, vc_old, dt_old, &
            max_timestep_factor, min_timestep_factor
         logical :: have_changed_matrix_coeffs, trace, last_step, &
            solved, do_timing, use_isolve


         logical, parameter :: dbg = .false.
         !logical, parameter :: dbg = .true.


         include 'formats'

         steps_used = 0
         total_num_iters = 0
         total_num_retries = 0

         use_isolve = s% diffusion_use_isolve

         min_num_substeps = s% diffusion_min_num_substeps
         trace = s% show_diffusion_substep_info
         do_timing = s% show_diffusion_timing
         nullify(limit_coeffs_face)

         if (dbg) then
            write(*,2) 'nzlo', nzlo
            write(*,2) 'nzhi', nzhi
            write(*,2) 'nz', nz
            write(*,*)
         end if


         min_lambda = 1d-2

         if (do_timing) then
            call system_clock(time0,clock_rate)
         else
            time0 = 0
         endif

         ierr = 0
         ierr_dealloc = 0
         if (m /= nc+1) then
            ierr = -1
            write(*,*) 'm /= nc+1'
            return
         end if

         if (T(1) <= max_T_for_radaccel) then
            if (dbg) write(*,*) 'call load_op_mono_data'
            call load_op_mono_data( &
               op_mono_data_path, op_mono_data_cache_filename, ierr)
            if (dbg) write(*,*) 'done load_op_mono_data'
            if (ierr /= 0) then
               write(*,*) 'error while loading OP data, ierr = ',ierr
               return
            end if
         end if

         X_total_atol = s% diffusion_X_total_atol
         X_total_rtol = s% diffusion_X_total_rtol

         h1 = net_iso(ih1)
         if (h1 == 0) then
            ierr = -1; write(*,*) 'isos must include h1 for diffusion'; return
         end if

         he4 = net_iso(ihe4)
         if (he4 == 0) then
            ierr = -1; write(*,*) 'isos must include he4 for diffusion'; return
         end if

         !reset nzlo if necessary
         nbound=nzlo
         do k=nzlo,nzhi
            if (T(k) > T_full_off) then
               nbound=k; exit
            endif
         enddo
         nzlo=nbound
         
         call do_alloc1(ierr)
         if (ierr /= 0) then
            if (dbg .or. s% report_ierr) write(*,*) 'diffusion failed in do_alloc1'
            return
         end if
         
         call get_limit_coeffs( &
            s, nz, nzlo, nzhi, &
            gamma_full_on, gamma_full_off, &
            T_full_on, T_full_off, &
            gamma, T, limit_coeffs_face, k_max)
         if (dbg) write(*,4) 'k_max for limit_coeffs_face /= 0', &
               k_max, nzhi, nzhi - k_max
         nzhi = k_max ! max k s.t. limit_coeffs_face(k) > 0

         n = nzhi-nzlo+1

         if (dbg) write(*,2) '  nz', nz
         if (dbg) write(*,2) 'nzlo, q', nzlo, s% q(nzlo)
         if (dbg) write(*,2) 'nzhi, q', nzhi, s% q(nzhi)
         if (dbg) write(*,2) '   n', n
         if (n <= 1) then
            call do_dealloc1
            return
         end if

         neqs = n*nc
         idiag = 2*nc
         ku = 2*nc - 1
         kl = ku
         ldab = ku + kl + 1
         ldafb = kl + ldab
         ldb = neqs
         ldx = neqs

         upwind_limit = s% diffusion_upwind_abs_v_limit

         mstar = sum(dm_in(1:nz))
         do i=1,species
            xa_total_before(i) = &
               dot_product(dm_in(1:nz),xa(i,1:nz))/mstar
         end do

         call do_alloc(ierr)
         if (ierr /= 0) then
            if (dbg .or. s% report_ierr) write(*,*) 'diffusion failed in do_alloc'
            call do_dealloc1
            return
         end if

         do k=2,nz
            alfa_face(k) = dm_in(k-1)/(dm_in(k) + dm_in(k-1))
         end do
         alfa_face(1) = 1d0
         ! e.g., T_face(k) = alfa_face(k)*T(k) + (1d0-alfa_face(k))*T(k-1)

         do k=nzlo+1,nzhi
            r_face(k) = r_face_in(k)
         end do
         r_face(nzlo) = r_face_in(1)
         r_face(nzhi+1) = R_center

         ! create classes
         A(1:nc) = chem_isos% Z_plus_N(class_chem_id(1:nc))
         A(m) = me/amu
         do k=1,nz
            cell_dm(k) = dm_in(k)
            dm_bar(k) = dm_bar_in(k)
            do j=1,nc
               X(j,k) = 0
            end do
            do j=1,species
               i = class(j)
               X(i,k) = X(i,k) + max(tiny_X, xa(j,k))
            end do
            do i=1,nc
               if (X(i,k) <= 0) X(i,k) = tiny_X
            end do
            tmp = sum(X(1:nc,k))
            do j=1,nc
               X(j,k) = X(j,k) / tmp
               X_init(j,k) = X(j,k)
            end do
            X(m,k) = 0d0 ! this will be set later
         end do

         if (nzlo > 1) then ! combine cells 1 to nzlo
            mtotal = sum(cell_dm(1:nzlo))
            do j=1,species ! change to species average of 1:nzlo
               xa(j,nzlo) = dot_product(cell_dm(1:nzlo),xa(j,1:nzlo))/mtotal
            end do
            do j=1,nc ! change to class average
               X(j,nzlo) = dot_product(cell_dm(1:nzlo),X(j,1:nzlo))/mtotal
            end do
            cell_dm(nzlo) = mtotal
            sumx = sum(X(1:nc,nzlo))
            do j=1,nc
               X(j,nzlo) = X(j,nzlo)/sumx
               X_init(j,nzlo) = X(j,nzlo)
            end do
         end if
         dm_bar(nzlo+1) = 0.5d0*(cell_dm(nzlo) + cell_dm(nzlo+1))

         do i=1,nc ! save for species conservation checks
            mass_init(i) = dot_product(cell_dm(nzlo:nzhi),X(i,nzlo:nzhi))
         end do
         sum_mass_nzlo_nzhi = sum(cell_dm(nzlo:nzhi))

         ! get info that is held constant during solve
         if (dbg) write(*,*) 'call setup_struct_info'
         call setup_struct_info( &
            s, nz, nzlo, nzhi, species, nc, m, X, A, tiny_X,  &
            dlnP_dm_face, dlnT_dm_face, dlnRho_dm_face, cell_dm, dm_in, &
            abar, free_e, T, lnT, rho, lnd, L_face, r_face, alfa_face, &
            class, class_chem_id, calculate_ionization, nsmooth_typical_charge, &
            min_T_for_radaccel, max_T_for_radaccel, &
            min_Z_for_radaccel, max_Z_for_radaccel, &
            screening_for_radaccel, rho_face, T_face, four_pi_r2_rho_face, &
            dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, &
            Z, typical_charge, xm_face, &
            rad_accel_face, log10_g_rad, g_rad_face, &
            kmax_rad_accel, ierr)
         if (dbg) write(*,*) 'done setup_struct_info'
         if (failed('setup_struct_info')) then
            call dealloc
            return
         end if

         mtotal = sum(cell_dm(nzlo:nzhi))
         do j=1,nc
            xtotal_init(j) = &
               dot_product(cell_dm(nzlo:nzhi),X(j,nzlo:nzhi))/mtotal
         end do

         tol_correction_max = s% diffusion_tol_correction_max
         tol_correction_norm = s% diffusion_tol_correction_norm
         dt = 0
         min_dt = total_time/(1000*max_steps)
         time = 0
         steps_used = 0
         total_num_retries = 0
         total_num_iters = 0

         if (dbg) write(*,*) '1st time call fixup'
         call fixup(s,  &
            nz, nc, m, nzlo, nzhi, &
            total_num_iters, s% diffusion_min_X_hard_limit, X_total_atol, X_total_rtol, &
            cell_dm, mtotal, xtotal_init, X, &
            lnT, sum_ending_mass, ending_mass, ending_dX_dm, &
            bad_j, bad_k, bad_X, bad_sum, bad_Xsum, ierr)
         if (failed('fixup')) then
            call dealloc
            return
         end if
         have_changed_matrix_coeffs = .false.

         vc_old = 0
         dt_old = 0
         avg_del = 0

         if (use_isolve) then
            call solve_with_isolve(ierr)
         else
            call do_step_loop(ierr)
         end if
         if (ierr /= 0) then
            call dealloc
            return
         end if

         if (total_time - time > 1d-10*total_time .or. steps_used == 0) then
            ierr = -1 ! failed to finish
         end if

         if (ierr == 0) then

            call set_new_xa(nz, nzlo, nzhi, species, nc, m, class, X_init, X, cell_dm, xa)

            do k=nzlo,nzhi
               do j=1,nc
                  X_final(j,k) = X(j,k)
               end do
            end do

            do k=1,nzlo-1 ! fully mixed
               do i=1,species
                  xa(i,k) = xa(i,nzlo)
               end do
               do j=1,nc
                  X_final(j,k) = X(j,nzlo)
               end do
            end do

         end if

         call dealloc

         if (dbg .and. ierr /= 0) then
            stop 'failed in diffusion'
         end if

         if (do_timing) then
            call system_clock(time1,clock_rate)
            write(*,2) 'diffusion model timing', s% model_number, &
               dble(time1 - time0)/clock_rate
         end if


         contains


         subroutine solve_with_isolve(ierr)
            use num_def
            use num_lib
            use mtx_def
            use mtx_lib

            integer, intent(out) :: ierr

            integer :: &
               i, k, nsteps, lout, iout, idid, ijac, max_steps, &
               imas, mlmas, mumas, itol, &
               nzmax, lrd, lid, isparse, &
               liwork, lwork, caller_id, which_solver, which_decsol
            real(dp) :: h, max_step_size, x_min, x_max

            integer, pointer :: iwork(:) !(liwork)
            real(dp), pointer :: work(:) !(lwork)
            integer, pointer :: ipar_decsol(:) !(lid)
            real(dp), pointer :: rpar_decsol(:) !(lrd)

            real(dp) :: atol(1), rtol(1), t, tend
            real(dp), dimension(:), pointer :: lblk, dblk, ublk ! =(nc,nc,n)
            real(dp), dimension(:), pointer :: uf_lblk, uf_dblk, uf_ublk ! =(nc,nc,n)

            integer, parameter :: lrpar = 1, lipar = 1
            real(dp), target :: rpar_ary(lrpar)
            integer, target :: ipar_ary(lipar)
            real(dp), pointer :: rpar(:)
            integer, pointer :: ipar(:)

            real(dp), pointer :: vars1(:), vars(:,:)

            include 'formats'
            ierr = 0

            rpar => rpar_ary
            ipar => ipar_ary

            which_decsol = bcyclic_dble

            which_solver = solver_option(s% diffusion_isolve_solver, ierr)
            if (ierr /= 0) then
               write(*,*) 'unknown solver name for diffusion_isolve_solver' // &
                  trim(s% diffusion_isolve_solver)
               return
            end if

            allocate(vars1(neqs))
            if (which_decsol == bcyclic_dble) then
               allocate(lblk(nc*nc*n), dblk(nc*nc*n), ublk(nc*nc*n))
               allocate(uf_lblk(nc*nc*n), uf_dblk(nc*nc*n), uf_ublk(nc*nc*n))
            else
               nullify(lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk)
            end if

            itol = 0
            rtol(1) = s% diffusion_rtol_for_isolve
            atol(1) = s% diffusion_atol_for_isolve

            iout = 1
            max_steps = s% diffusion_maxsteps_for_isolve
            isparse = 0
            lout = 6
            caller_id = 0

            x_min = s% diffusion_min_X_hard_limit
            x_max = 1d0 - s% diffusion_min_X_hard_limit

            ijac = 1 ! analytic jacobian

            imas = 0
            mlmas = 0
            mumas = 0

            ipar = 0
            rpar = 0

            lid = 0; lrd = 0
            if (which_decsol == lapack) then
               stop 'not ready for lapack in diff'
               nzmax = 0
               call lapack_work_sizes(neqs,lrd,lid)
            else if (which_decsol == bcyclic_dble) then
               nzmax = 0
               call bcyclic_dble_work_sizes(nc,n,lrd,lid)
            else
               write(*,*) 'unknown value for which_decsol', which_decsol
               call mesa_error(__FILE__,__LINE__)
            end if

            call isolve_work_sizes(neqs,nzmax,imas,kl,ku,mlmas,mumas,liwork,lwork)

            allocate(iwork(liwork),work(lwork),ipar_decsol(lid),rpar_decsol(lrd),stat=ierr)
            if (ierr /= 0) then
               write(*,*) 'allocate ierr', ierr
               call mesa_error(__FILE__,__LINE__) ! test_int_support
            end if

            iwork = 0
            work = 0

            t = 0
            tend = total_time
            h = total_time ! initial step size
            max_step_size = total_time

            ! set vars
            vars(1:nc,1:n) => vars1(1:neqs)
            do i=1,n
               k = i+nzlo-1
               do j=1,nc
                  vars(j,i) = X(j,k)
               end do
            end do

            if (which_decsol == lapack) then
               stop 'not supported'
               call isolve( &
                  which_solver, neqs, fcn, t, vars1, tend, &
                  h, max_step_size, max_steps, &
                  rtol, atol, itol, x_min, x_max, &
                  jac, ijac, null_sjac, nzmax, isparse, kl, ku, &
                  null_mas, imas, mlmas, mumas, &
                  solout, iout, &
                  lapack_decsol, null_decsols, null_decsolblk, &
                  lrd, rpar_decsol, lid, ipar_decsol, &
                  caller_id, 0, 0, &
                  lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, &
                  null_fcn_blk_dble, null_jac_blk_dble, &
                  work, lwork, iwork, liwork, &
                  lrpar, rpar, lipar, ipar, &
                  lout, idid)
            else if (which_decsol == bcyclic_dble) then
               call isolve( &
                  which_solver, neqs, null_fcn, t, vars1, tend, &
                  h, max_step_size, max_steps, &
                  rtol, atol, itol, x_min, x_max, &
                  null_jac, ijac, null_sjac, nzmax, isparse, kl, ku, &
                  null_mas, imas, mlmas, mumas, &
                  solout, iout, &
                  null_decsol, null_decsols, bcyclic_dble_decsolblk, &
                  lrd, rpar_decsol, lid, ipar_decsol, &
                  caller_id, nc, n, &
                  lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, &
                  fcn_blk_dble, jac_blk_dble, &
                  work, lwork, iwork, liwork, &
                  lrpar, rpar, lipar, ipar, &
                  lout, idid)
            else
               stop 'bad which_decsol in mod_diffusion'
            end if

            if (idid /= 1) ierr = -1
            !if (ierr /= 0) then
            !   write(*,*) 'solver returned ierr /= 0', idid
            !   call mesa_error(__FILE__,__LINE__)
            !end if
            steps_used = iwork(17) ! number of accepted steps
            total_num_retries = iwork(16) - iwork(17) ! num computed - num accepted
            time = tend

            do i=1,n
               k = i+nzlo-1
               do j=1,nc
                  X(j,k) = vars(j,i)
               end do
            end do

            if (which_decsol == bcyclic_dble) &
               deallocate(lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk)
            deallocate(vars1,iwork,work,ipar_decsol,rpar_decsol)


         end subroutine solve_with_isolve


         subroutine fcn_blk_dble(&
               n_blk_dble,caller_id,nvar_blk_dble,nz_blk_dble,&
               time,h,vars,f,lrpar,rpar,lipar,ipar,ierr)
            use const_def, only: dp
            integer, intent(in) :: n_blk_dble, caller_id, &
               nvar_blk_dble, nz_blk_dble, lrpar, lipar
            real(dp), intent(in) :: time,h
            real(dp), intent(inout), pointer :: vars(:) ! (n)
            real(dp), intent(inout), pointer :: f(:) ! (n) dvars/dt
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.

            integer :: i, j, k
            real(dp), pointer :: f2(:,:)

            include 'formats'

            ierr = 0

            if (n_blk_dble /= neqs .or. nvar_blk_dble /= nc .or. nz_blk_dble /= n) &
               stop 'bad args for fcn_blk_dble'

            call update_for_new_vars(vars, time, ierr)
            if (ierr /= 0) then
               !stop 'fcn_blk_dble'
               return
            end if

            call get_dX_dt( &
               nz, nzlo, nzhi, m, nc, X, X_face, C_div_X, GT_face, AD_face, SIG_face, &
               alfa_face, cell_dm, v_advection_face, upwind_limit, dX_dt)

            f2(1:nc,1:n) => f(1:neqs)
            do i=1,n
               k = i+nzlo-1
               do j=1,nc
                  f2(j,i) = dX_dt(j,k)
               end do
            end do

         end subroutine fcn_blk_dble


         subroutine jac_blk_dble( &
               n_blk_dble,caller_id,nvar_blk_dble,nz_blk_dble,&
               time,h,vars,f,lblk1,dblk1,ublk1,lrpar,rpar,lipar,ipar,ierr)
            use const_def,only: dp
            integer,intent(in) :: n_blk_dble,caller_id,nvar_blk_dble,nz_blk_dble,lrpar,lipar
            real(dp),intent(in) :: time,h
            real(dp),intent(inout), pointer :: vars(:) ! (n)
            real(dp),intent(inout), pointer :: f(:) ! (n) dvars/dt
            real(dp),dimension(:),pointer,intent(inout) :: lblk1,dblk1,ublk1
            integer,intent(inout),pointer :: ipar(:) ! (lipar)
            real(dp),intent(inout),pointer :: rpar(:) ! (lrpar)
            integer,intent(out) :: ierr ! nonzero means terminate integration

            integer :: i, j, jj, k
            real(dp), pointer :: f2(:,:)
            real(dp),dimension(:,:,:),pointer :: lblk,dblk,ublk

            include 'formats'

            ierr = 0

            if (n_blk_dble /= neqs .or. nvar_blk_dble /= nc .or. nz_blk_dble /= n) &
               stop 'bad args for jac_blk_dble'

            call update_for_new_vars(vars, time, ierr)
            if (ierr /= 0) then
               !stop 'jac_blk_dble'
               return
            end if

            call get_eqn_matrix_entries( &
               nz, nzlo, nzhi, m, nc, .true., X, X_face, C_div_X, &
               GT_face, AD_face, SIG_face, &
               alfa_face, cell_dm, v_advection_face, upwind_limit, &
               tiny_X, 0d0, dX, dX_dt, del, em1, e00, ep1)

            f2(1:nc,1:n) => f(1:neqs)
            lblk(1:nc,1:nc,1:n) => lblk1(1:nc*nc*n)
            dblk(1:nc,1:nc,1:n) => dblk1(1:nc*nc*n)
            ublk(1:nc,1:nc,1:n) => ublk1(1:nc*nc*n)
            do i=1,n
               k = i+nzlo-1
               do j=1,nc
                  f2(j,i) = dX_dt(j,k)
                  do jj=1,nc
                     lblk(jj,j,i) = em1(jj,j,k)
                     dblk(jj,j,i) = e00(jj,j,k)
                     ublk(jj,j,i) = ep1(jj,j,k)
                  end do
               end do
            end do

         end subroutine jac_blk_dble


         subroutine update_for_new_vars(vars, time, ierr)
            real(dp),intent(inout), pointer :: vars(:)
            real(dp), intent(in) :: time
            integer,intent(out) :: ierr

            integer :: i, j, k
            real(dp), pointer :: vars2(:,:)


            include 'formats'

            ierr = 0
            vars2(1:nc,1:n) => vars(1:neqs)

            ! copy vars to X
            do i=1,n
               k = i+nzlo-1
               do j=1,nc
                  X(j,k) = vars2(j,i)
               end do
            end do

            if (.true.) then

               call fixup(s,  &
                  nz, nc, m, nzlo, nzhi, total_num_iters, &
                  s% diffusion_min_X_hard_limit, X_total_atol, X_total_rtol, &
                  cell_dm, mtotal, xtotal_init, X, &
                  lnT, sum_ending_mass, ending_mass, ending_dX_dm, &
                  bad_j, bad_k, bad_X, bad_sum, bad_Xsum, ierr)
               if (ierr /= 0) then
                  s% retry_message = 'element diffusion failed in fixup'
                  if (s% report_ierr) write(*, *) s% retry_message
                  return
               end if
               if (failed('fixup')) return

               ! copy X to vars
               do i=1,n
                  k = i+nzlo-1
                  do j=1,nc
                     vars2(j,i) = X(j,k)
                  end do
               end do

            end if

            call get_matrix_coeffs( &
               s, nz, nc, m, nzlo, nzhi, class(h1), class(he4), &
               pure_Coulomb, dt, v_advection_max, tiny_C, diffusion_factor, &
               A, X, Z, rho_face, T_face, four_pi_r2_rho_face, &
               xm_face, cell_dm, dm_bar, dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, &
               r_face, r_mid, limit_coeffs_face, alfa_face, &
               rad_accel_face, log10_g_rad, g_rad_face, &
               min_T_for_radaccel, max_T_for_radaccel, &
               X_init, X_face, C, C_div_X, C_div_X_face, &
               e_ap, e_at, e_ar, e_ax, E_field_face, &
               g_ap, g_at, g_ar, g_ax, g_field_face, &
               v_advection_face, v_total_face, vlnP_face, vlnT_face, v_rad_face, &
               GT_face, D_self_face, AD_face, SIG_face, sigma_lnC, ierr)
            if (ierr /= 0) then
               write(*,1) 'failed in get_matrix_coeffs', time
               return
            end if
            if (failed('get_matrix_coeffs')) return

         end subroutine update_for_new_vars


         subroutine fcn(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
            use const_def, only: dp
            integer, intent(in) :: n, lrpar, lipar
            real(dp), intent(in) :: x, h
            real(dp), intent(inout) :: y(:)
            real(dp), intent(inout) :: f(:) ! dy/dx
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
            ierr = 0
            f = 0
         end subroutine fcn


         subroutine jac(n, x, h, y, f, dfdy, ldfy, lrpar, rpar, lipar, ipar, ierr)
            use const_def, only: dp
            integer, intent(in) :: n, ldfy, lrpar, lipar
            real(dp), intent(in) :: x, h
            real(dp), intent(inout) :: y(:)
            real(dp), intent(inout) :: f(:) ! dy/dx
            real(dp), intent(inout) :: dfdy(:,:)
               ! dense: dfdy(i, j) = partial f(i) / partial y(j)
               ! banded: dfdy(i-j+mujac+1, j) = partial f(i) / partial y(j)
                  ! uses rows 1 to mljac+mujac+1 of dfdy.
                  ! The j-th column of the square matrix is stored in the j-th column of the
                  ! array dfdy as follows:
                  ! dfdy(mujac+1+i-j, j) = partial f(i) / partial y(j)
                  ! for max(1, j-mujac)<=i<=min(N, j+mljac)
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr ! nonzero means terminate integration
            ierr = 0
            f = 0
            dfdy = 0
         end subroutine jac


         subroutine solout(nr, xold, time, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
            ! nr is the step number.
            ! x is the current x value; xold is the previous x value.
            ! y is the current y value.
            ! irtrn negative means terminate integration.
            ! rwork_y and iwork_y hold info for interp_y
            ! note that these are not the same as the rwork and iwork arrays for the solver.
            use const_def, only: dp
            integer, intent(in) :: nr, n, lrpar, lipar
            real(dp), intent(in) :: xold, time
            real(dp), intent(inout) :: y(:)
            ! y can be modified if necessary to keep it in valid range of possible solutions.
            real(dp), intent(inout), target :: rwork_y(*)
            integer, intent(inout), target :: iwork_y(*)
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            interface
               include 'num_interp_y.dek'
            end interface
            integer, intent(out) :: irtrn ! < 0 causes solver to return to calling program.
            include 'formats'
            if (s% show_diffusion_substep_info .and. mod(nr,10) == 0) write(*,2) 'nr', nr, time
            irtrn = 0
         end subroutine solout



         subroutine do_step_loop(ierr)
            integer, intent(out) :: ierr

            include 'formats'
            ierr = 0

      step_loop: do while &
               (total_time - time > 1d-10*total_time .and. &
                  steps_used < max_steps)

            steps_used = steps_used + 1

            if (dbg) write(*,*) 'call update_coeffs'
            call update_coeffs((kmax_rad_accel > 0 .and. steps_used > 1), ierr)
            if (failed('update_coeffs')) return

            remaining_time = total_time - time

            if (dbg) write(*,*) 'call get_timescale'
            call get_timescale( &
               s, nz, nzlo, nzhi, m, nc, -s% diffusion_min_X_hard_limit*0.5d0, v_advection_face, &
               upwind_limit, X, X_face, C_div_X, SIG_face, GT_face, AD_face, cell_dm, r_face, &
               steps_used, total_num_iters, dbg, iter_dbg, j_dbg, k_dbg, &
               j, k, class_chem_id, total_time, timescale)

            dt = dt_div_timescale*timescale
            dt = max(dt, 1d-4*remaining_time)
            if (min_num_substeps > 0) dt = min(dt, total_time/min_num_substeps)

            if (dbg) write(*,3) 'timescale remaining_time dt', &
               steps_used, total_num_iters, timescale, &
                  remaining_time, dt, dt_div_timescale

            if (dt >= remaining_time) then
               if (dbg) write(*,3) 'dt >= remaining_time', &
                  steps_used, total_num_iters, dt/remaining_time
               dt = remaining_time
            else if (dt > 0.5d0*remaining_time) then
               if (dbg) write(*,3) 'dt >= 0.5d0*remaining_time', &
                  steps_used, total_num_iters, dt/remaining_time
               dt = 0.5d0*remaining_time
            else
               if (dbg) write(*,3) 'dt < 0.5d0*remaining_time', &
                  steps_used, total_num_iters, dt/remaining_time
            end if

            last_step = time + dt >= total_time*(1d0 - 1d-10)

            if (dbg) write(*,3) &
               'dt, ../total, ../remaining, remaining/total', &
               steps_used, total_num_iters, dt, dt/total_time, &
               dt/remaining_time, remaining_time/total_time
            if (dbg) write(*,*)

            max_retries = s% diffusion_max_retries_per_substep
         retry_loop: do retry_count = 1, max_retries

               if (dt < min_dt) then
                  if (dbg) write(*,1) 'quit because dt < min_dt: min_dt', min_dt
                  exit step_loop
               end if

               if (retry_count == 1) then ! save starting X in X_0
                  if (trace) then
                     write(*,*)
                     write(*,2) '1st try for substep', steps_used
                  end if
                  do k=nzlo,nzhi
                     do j=1,nc
                        X_0(j,k) = X(j,k)
                        X_1(j,k) = X_0(j,k)
                        dX(j,k) = 0d0
                     end do
                  end do
               else ! prepare for retry of solve_loop
                  if (trace) then
                     write(*,*)
                     write(*,2) 'retry substep', steps_used
                  end if
                  dt = dt*dt_retry_factor
                  total_num_retries = total_num_retries+1
                  do k=nzlo,nzhi ! restore X from X_0
                     do j=1,nc
                        X(j,k) = X_0(j,k)
                        X_1(j,k) = X_0(j,k)
                        dX(j,k) = 0d0
                     end do
                  end do
                  if (have_changed_matrix_coeffs) then ! must recalculate them
                     if (dbg) write(*,*) 'call get_matrix_coeffs'
                     call get_matrix_coeffs( &
                        s, nz, nc, m, nzlo, nzhi, class(h1), class(he4), &
                        pure_Coulomb, dt, v_advection_max, tiny_C, diffusion_factor, &
                        A, X, Z, rho_face, T_face, four_pi_r2_rho_face, &
                        xm_face, cell_dm, dm_bar, dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, &
                        r_face, r_mid, limit_coeffs_face, alfa_face, &
                        rad_accel_face, log10_g_rad, g_rad_face, &
                        min_T_for_radaccel, max_T_for_radaccel, &
                        X_init, X_face, C, C_div_X, C_div_X_face, &
                        e_ap, e_at, e_ar, e_ax, E_field_face, &
                        g_ap, g_at, g_ar, g_ax, g_field_face, &
                        v_advection_face, v_total_face, vlnP_face, vlnT_face, v_rad_face, &
                        GT_face, D_self_face, AD_face, SIG_face, sigma_lnC, ierr)
                     if (failed('get_matrix_coeffs')) return
                  end if
               end if

               have_changed_matrix_coeffs = .false.

               if (dbg) write(*,*) 'call BE_solve'
               solved = BE_solve(1d0, last_step, ierr)
               if (ierr /= 0) exit step_loop

               if (solved) then ! this step is done
                  time = time + dt
                  cycle step_loop
               end if

            end do retry_loop

            exit step_loop

         end do step_loop

         end subroutine do_step_loop


         logical function BE_solve(theta, do_smooth, ierr) ! Backwards Euler
            real(dp), intent(in) :: theta ! multiplier for dt
            logical, intent(in) :: do_smooth
            integer, intent(out) :: ierr

            include 'formats'

            BE_solve = .false.

         solve_loop: do num_iters = 1, max_iters_per_substep

               if (total_num_iters >= max_iters_total) then
                  ierr = -1
                  return
               end if

               total_num_iters = total_num_iters+1

               call get_eqn_matrix_entries( &
                  nz, nzlo, nzhi, m, nc, .false., X, X_face, C_div_X, &
                  GT_face, AD_face, SIG_face, &
                  alfa_face, cell_dm, v_advection_face, upwind_limit, &
                  tiny_X, theta*dt, dX, dX_dt, del, em1, e00, ep1)

               call factor_and_solve_matrix( &
                  nz, nzlo, nzhi, nc, n, neqs, &
                  lblk1, dblk1, ublk1, ipiv1, del1, &
                  lrd, rpar_decsol, lid, ipar_decsol, ierr)
               if (ierr /= 0) then
                  if (dbg .or. trace) write(*,*) 'retry because failed to solve matrix eqn'
                  ierr = 0; return
               end if

               ! set lambda for positivity
               call positivity( &
                  nc, nzlo, nzhi, X_1, -s% diffusion_min_X_hard_limit, &
                  del, lambda, num_iters)
               if (dbg) write(*,2) 'lambda', num_iters, lambda

               if (lambda < min_lambda) then
                  lambda = min_lambda
                  if (dbg) write(*,2) 'min lambda', num_iters, lambda
               end if

               do k=nzlo,nzhi ! apply the correction
                  do j=1,nc
                     dX(j,k) = dX(j,k) + lambda*del(j,k)
                     X_1(j,k) = X_0(j,k) + dX(j,k)
                     X(j,k) = X_1(j,k)
                  end do
               end do

               call fixup(s,  &
                  nz, nc, m, nzlo, nzhi, total_num_iters, &
                  s% diffusion_min_X_hard_limit, X_total_atol, X_total_rtol, &
                  cell_dm, mtotal, xtotal_init, X, &
                  lnT, sum_ending_mass, ending_mass, ending_dX_dm, &
                  bad_j, bad_k, bad_X, bad_sum, bad_Xsum, ierr)
               if (ierr /= 0) then
                  if (dbg .or. trace) then
                     write(*,'(a,6i5,3f16.8,1p,99e16.8)') &
                        '(fixup failed) retry: step try iter total_iters j k X sum Xsum xm', &
                           steps_used, retry_count, num_iters, total_num_iters, &
                           bad_j, bad_k, bad_X, bad_sum, bad_Xsum, xm_face(bad_k)/Msun
                  end if

                  if (dbg .and. bad_sum /= 0) then
                     stop 'bad sum from fixup'

                  end if
                  ierr = 0; return
               end if

               if (lambda == 1d0) then
                  ! if magnitude of correction is small enough, consider converged.
                  max_del = maxval(abs(del(1:nc,nzlo:nzhi)))
                  avg_del = sum(abs(del(1:nc,nzlo:nzhi)))/(nc*n)
                  if (max_del <= tol_correction_max .and. avg_del <= tol_correction_norm) then
                     if (dbg .or. trace) &
                        write(*,3) 'substep converged: iters max_del avg_del dt/total remain/total', &
                           steps_used, num_iters, max_del, avg_del, dt/total_time, &
                           remaining_time/total_time
                     BE_solve = .true.
                     return
                  end if
                  if (dbg .or. trace) then
                     if (max_del > tol_correction_max) &
                        write(*,2) 'max_del > tol_correction_max', num_iters, max_del, tol_correction_max
                     if (avg_del > tol_correction_norm) &
                        write(*,2) 'avg_del > tol_correction_norm', num_iters, avg_del, tol_correction_norm
                  end if
               end if

               if (num_iters == max_iters_per_substep) then
                  if (dbg .or. trace) then
                        write(*,3) 'retry because num_iters == max_iters_per_substep', &
                           num_iters, max_iters_per_substep
                  end if
                  return
               end if

               ! prepare for next iteration of solve_loop
               call get_matrix_coeffs( &
                  s, nz, nc, m, nzlo, nzhi, class(h1), class(he4), &
                  pure_Coulomb, dt, v_advection_max, tiny_C, diffusion_factor, &
                  A, X, Z, rho_face, T_face, four_pi_r2_rho_face, &
                  xm_face, cell_dm, dm_bar, dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, &
                  r_face, r_mid, limit_coeffs_face, alfa_face, &
                  rad_accel_face, log10_g_rad, g_rad_face, &
                  min_T_for_radaccel, max_T_for_radaccel, &
                  X_init, X_face, C, C_div_X, C_div_X_face, &
                  e_ap, e_at, e_ar, e_ax, E_field_face, &
                  g_ap, g_at, g_ar, g_ax, g_field_face, &
                  v_advection_face, v_total_face, vlnP_face, vlnT_face, v_rad_face, &
                  GT_face, D_self_face, AD_face, SIG_face, sigma_lnC, ierr)
               if (failed('get_matrix_coeffs')) return
               have_changed_matrix_coeffs = .true.

            end do solve_loop

            stop 'diffusion bug -- should not fall through solve_loop'

         end function BE_solve


         subroutine update_coeffs(update_g_rad, ierr)
            logical, intent(in) :: update_g_rad
            integer, intent(out) :: ierr

            include 'formats'

            ierr = 0

            if (update_g_rad) then
               call update_rad_accel_face( &
                  nzlo, nzhi, nc, m, A, X_init, X, &
                  log10_g_rad, g_rad_face, &
                  rad_accel_face, kmax_rad_accel)
            end if

            call get_matrix_coeffs( &
               s, nz, nc, m, nzlo, nzhi, class(h1), class(he4), &
               pure_Coulomb, dt, v_advection_max, tiny_C, diffusion_factor, &
               A, X, Z, rho_face, T_face, four_pi_r2_rho_face, &
               xm_face, cell_dm, dm_bar, dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, &
               r_face, r_mid, limit_coeffs_face, alfa_face, &
               rad_accel_face, log10_g_rad, g_rad_face, &
               min_T_for_radaccel, max_T_for_radaccel, &
               X_init, X_face, C, C_div_X, C_div_X_face, &
               e_ap, e_at, e_ar, e_ax, E_field_face, &
               g_ap, g_at, g_ar, g_ax, g_field_face, &
               v_advection_face, v_total_face, vlnP_face, vlnT_face, v_rad_face, &
               GT_face, D_self_face, AD_face, SIG_face, sigma_lnC, ierr)

         end subroutine update_coeffs


         subroutine do_alloc1(ierr)
            use alloc, only: work_array
            integer, intent(out) :: ierr
            ierr = 0
            call work_array(s, .true., .false., &
               limit_coeffs_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
         end subroutine do_alloc1

         subroutine do_dealloc1
            use alloc, only: work_array
            call work_array(s, .false., .false., &
               limit_coeffs_face, nz, nz_alloc_extra, 'mod_diffusion', ierr_dealloc)
         end subroutine do_dealloc1

         subroutine do_alloc(ierr)
            use mtx_lib, only: bcyclic_dble_work_sizes
            use alloc, only: get_integer_work_array
            integer, intent(out) :: ierr
            ierr = 0
            call bcyclic_dble_work_sizes(nc, nz, lrd, lid)
            
            call get_integer_work_array(s, ipiv1, nc*nz, nc*nz_alloc_extra, ierr)
            if (ierr /= 0) return
            call get_integer_work_array(s, ipar_decsol, lid, 0, ierr)
            if (ierr /= 0) return
            call do_work_arrays(.true.,ierr)
            if (ierr /= 0) return            
            
            g_ax(1:m,1:nz) => g_ax1(1:m*nz)
            Z(1:m,1:nz) => Z1(1:m*nz)
            X(1:m,1:nz) => X1(1:m*nz)
            C(1:m,1:nz) => C1(1:m*nz)
            dC_dr(1:m,1:nz) => dC_dr1(1:m*nz)
            C_div_X(1:m,1:nz) => C_div_X1(1:m*nz)
            X_face(1:m,1:nz) => X_face1(1:m*nz)
            C_div_X_face(1:m,1:nz) => C_div_X_face1(1:m*nz)
            rad_accel_face(1:m,1:nz) => rad_accel_face1(1:m*nz)
            log10_g_rad(1:m,1:nz) => log10_g_rad1(1:m*nz)
            new_mass(1:nc,1:nz) => new_mass1(1:nc*nz)
            starting_dX_dm(1:nc,1:nz) => starting_dX_dm1(1:nc*nz)
            starting_mass(1:nc,1:nz) => starting_mass1(1:nc*nz)
            ending_mass(1:nc,1:nz) => ending_mass1(1:nc*nz)
            ending_dX_dm(1:nc,1:nz) => ending_dX_dm1(1:nc*nz)
            e_ax(1:m,1:nz) => e_ax1(1:m*nz)
            GT_face(1:nc,1:nz) => GT_face1(1:nc*nz)
            X_theta(1:nc,1:nz) => X_theta1(1:nc*nz)
            X_minus_theta(1:nc,1:nz) => X_minus_theta1(1:nc*nz)
            del(1:nc,1:nz) => del1(1:nc*nz)
            X_0(1:nc,1:nz) => X_01(1:nc*nz)
            X_1(1:nc,1:nz) => X_11(1:nc*nz)
            dX(1:nc,1:nz) => dX1(1:nc*nz)
            dX_dt(1:nc,1:nz) => dX_dt1(1:nc*nz)
            SIG_face(1:nc,1:nc,1:nz) => SIG_face1(1:nc*nc*nz)
            sigma_lnC(1:nc,1:nc,1:nz) => sigma_lnC1(1:nc*nc*nz)
            em1(1:nc,1:nc,1:nz) => lblk1(1:nc*nc*nz)
            e00(1:nc,1:nc,1:nz) => dblk1(1:nc*nc*nz)
            ep1(1:nc,1:nc,1:nz) => ublk1(1:nc*nc*nz)
            rhs(1:nc,1:nz) => rhs1(1:nc*nz)
            
         end subroutine do_alloc


         subroutine dealloc
            use alloc, only: return_integer_work_array
            call do_dealloc1
            call return_integer_work_array(s, ipiv1)
            call return_integer_work_array(s, ipar_decsol)
            call do_work_arrays(.false.,ierr_dealloc)
         end subroutine dealloc


         subroutine do_work_arrays(alloc_flag, ierr)
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0

            call work_array(s, alloc_flag, crit, &
               Z1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               X1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               C1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dC_dr1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               C_div_X1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               X_face1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               C_div_X_face1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               rad_accel_face1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               log10_g_rad1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               sum_ending_mass, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               cell_dm, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dm_bar, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               sum_new_mass, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               sum_starting_mass, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               alfa_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               rho_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               T_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               four_pi_r2_rho_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dlnP_dr_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dlnT_dr_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dlnRho_dr_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               r_face, nz+1, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               xm_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               del1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               X_01, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               X_11, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dX1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dX_dt1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               rhs1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               b1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               lblk1, nc*nc*nz, nc*nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dblk1, nc*nc*nz, nc*nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               ublk1, nc*nc*nz, nc*nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               rpar_decsol, lrd, 0, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               A, m, 0, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            ! for e field
            call work_array(s, alloc_flag, crit, &
               e_ap, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               e_at, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               e_ar, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               e_ax1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            ! for g field
            call work_array(s, alloc_flag, crit, &
               g_ap, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               g_at, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               g_ar, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               g_ax1, m*nz, m*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               new_mass1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               starting_dX_dm1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               starting_mass1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               ending_mass1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               ending_dX_dm1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               GT_face1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               AD_face, nz, nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               X_theta1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               X_minus_theta1, nc*nz, nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

            call work_array(s, alloc_flag, crit, &
               SIG_face1, nc*nc*nz, nc*nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               sigma_lnC1, nc*nc*nz, nc*nc*nz_alloc_extra, 'mod_diffusion', ierr)
            if (ierr /= 0) return

         end subroutine do_work_arrays


         logical function failed(str)
            character (len=*) :: str
            failed = .false.
            if (ierr == 0) return
            if (dbg .or. s% report_ierr) write(*,*) 'failed in ' // trim(str), ierr
            call dealloc
            failed = .true.
         end function failed

      end subroutine do_solve_diffusion


      subroutine get_timescale( &
            s, nz, nzlo, nzhi, m, nc, eps, v_advection_face, upwind_limit, &
            X, X_face, C_div_X, SIG_face, GT_face, AD_face, cell_dm, r_face, &
            steps_used, iter, dbg, iter_dbg, i_dbg, k_dbg, i_t, k_t, &
            class_chem_id, total_time, dt)
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nzlo, nzhi, m, nc
         real(dp), intent(in) :: eps, upwind_limit, total_time
         real(dp), intent(in), dimension(:,:) :: X, X_face, C_div_X, v_advection_face
         real(dp), intent(in) :: SIG_face(:,:,:) ! (nc,nc,nz)
         real(dp), intent(in) :: GT_face(:,:) ! (nc,nz)
         real(dp), intent(in) :: AD_face(:) ! (nz)
         real(dp), intent(in) :: cell_dm(:) ! (nz)
         real(dp), intent(in) :: r_face(:) ! (nz)
         integer, intent(in) :: steps_used, iter, class_chem_id(:)
         logical, intent(in) :: dbg
         integer, intent(in) :: iter_dbg, i_dbg, k_dbg
         integer, intent(out) :: i_t, k_t
         real(dp), intent(out) :: dt
         integer :: i, j, k
         real(dp) :: xdm, f, flow_in, flow_out, dt1, dxdt, coeff, &
            flow_in_GT, flow_out_GT, flow_in_SIG, flow_out_SIG, &
            flow_in_max, flow_out_max, &
            flow_in_GT_max, flow_out_GT_max, &
            flow_in_SIG_max, flow_out_SIG_max

         include 'formats'

         dt = 1d99
         flow_in_max = 0
         flow_out_max = 0
         flow_in_GT_max = 0
         flow_out_GT_max = 0
         flow_in_SIG_max = 0
         flow_out_SIG_max = 0

         do k=nzlo+1,nzhi-1

            do i=1,nc

               if (X(i,k) < eps) cycle

               if (abs(v_advection_face(i,k)) >= upwind_limit) then
                  if (GT_face(i,k) > 0) then
                     flow_out = GT_face(i,k)*X(i,k)
                  else
                     flow_out = GT_face(i,k)*X(i,k-1)
                  end if
               else
                  flow_out = GT_face(i,k)*X_face(i,k)
               end if

               if (abs(v_advection_face(i,k+1)) >= upwind_limit) then
                  if (GT_face(i,k+1) > 0) then
                     flow_in = GT_face(i,k+1)*X(i,k+1)
                  else
                     flow_in = GT_face(i,k+1)*X(i,k)
                  end if
               else
                  flow_in = GT_face(i,k+1)*X_face(i,k+1)
               end if

               flow_in_GT = flow_in
               flow_out_GT = flow_out
               flow_in_SIG = 0
               flow_out_SIG = 0

               do j=1,nc
                  coeff = SIG_face(i,j,k+1)*X_face(i,k+1)/max(smallX,X_face(j,k+1))
                  if (j == i) coeff = coeff + AD_face(k+1)
                  f = coeff*(C_div_X(j,k+1)*X(j,k+1) - C_div_X(j,k)*X(j,k))
                  flow_in = flow_in + f
                  flow_in_SIG = flow_in_SIG + f
                  coeff = SIG_face(i,j,k)*X_face(i,k)/max(smallX,X_face(j,k))
                  if (j == i) coeff = coeff + AD_face(k)
                  f = coeff*(C_div_X(j,k)*X(j,k) - C_div_X(j,k-1)*X(j,k-1))
                  flow_out = flow_out + f
                  flow_out_SIG = flow_out_SIG + f
               end do

               if (dbg .and. i == i_dbg .and. k == k_dbg .and. iter == iter_dbg) then
                  write(*,*)
                  write(*,*) 'dgb get_timescale'
                  call show_stuff
                  write(*,*)
               end if

               dxdt = (flow_in - flow_out)/cell_dm(k)
               if (dxdt >= 0d0) cycle ! only consider decreasing abundances
               dt1 = (X(i,k) + eps)/max(1d-50,-dxdt)

               if (dt1 < dt) then
                  dt = dt1
                  k_t = k
                  i_t = i
                  flow_in_max = flow_in
                  flow_out_max = flow_out
                  flow_in_GT_max = flow_in_GT
                  flow_out_GT_max = flow_out_GT
                  flow_in_SIG_max = flow_in_SIG
                  flow_out_SIG_max = flow_out_SIG
               end if

            end do

         end do

         if (dbg) then

            i = i_t
            k = k_t
            flow_in = flow_in_max
            flow_out = flow_out_max
            flow_in_GT = flow_in_GT_max
            flow_out_GT = flow_out_GT_max
            flow_in_SIG = flow_in_SIG_max
            flow_out_SIG = flow_out_SIG_max

            call show_stuff

            stop 'get_timescale'

         end if


         contains


         subroutine show_stuff
            use chem_def, only: chem_isos

            real(dp) :: dr

            include 'formats'
            dr = 0.5d0*(r_face(k-1) - r_face(k+1))

            write(*,*)
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'timescale total_time', &
               steps_used, iter, i, k, dt, total_time
            write(*,'(a30,4i6)') trim(chem_isos% name(class_chem_id(i))) // &
               ' nzlo, nzhi, nz', nzlo, nzhi, nz
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'dX', &
               steps_used, iter, i, k, &
               dt*(flow_in - flow_out)/cell_dm(k)
            write(*,'(a30,4i6,99e14.5)') 'dX GT', &
               steps_used, iter, i, k, &
               dt*(flow_in_GT - flow_out_GT)/cell_dm(k)
            write(*,'(a30,4i6,99e14.5)') 'dX SIG', &
               steps_used, iter, i, k, &
               dt*(flow_in_SIG - flow_out_SIG)/cell_dm(k)
            write(*,*)

            write(*,'(a30,4i6,99e14.5)') 'dX in', &
               steps_used, iter, i, k, &
               dt*flow_in/cell_dm(k)
            write(*,'(a30,4i6,99e14.5)') 'dX out', &
               steps_used, iter, i, k, &
               dt*flow_out/cell_dm(k)
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'rho', &
               steps_used, iter, i, k+1, s% rho(k+1)
            write(*,'(a30,4i6,99e14.5)') 'rho', &
               steps_used, iter, i, k, s% rho(k)
            write(*,'(a30,4i6,99e14.5)') 'rho', &
               steps_used, iter, i, k-1, s% rho(k-1)
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'C', &
               steps_used, iter, i, k+1, X(i,k+1)*C_div_X(i,k+1)
            write(*,'(a30,4i6,99e14.5)') 'C', &
               steps_used, iter, i, k, X(i,k)*C_div_X(i,k)
            write(*,'(a30,4i6,99e14.5)') 'C', &
               steps_used, iter, i, k-1, X(i,k-1)*C_div_X(i,k-1)
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'X', &
               steps_used, iter, i, k+1, X(i,k+1)
            write(*,'(a30,4i6,99e14.5)') 'X', &
               steps_used, iter, i, k, X(i,k)
            write(*,'(a30,4i6,99e14.5)') 'X', &
               steps_used, iter, i, k-1, X(i,k-1)
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'X_face', &
               steps_used, iter, i, k+1, X_face(i,k+1)
            write(*,'(a30,4i6,99e14.5)') 'X_face', &
               steps_used, iter, i, k, X_face(i,k)
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'dt*v_advection_face/dr', &
               steps_used, iter, i, k+1, dt*v_advection_face(i,k+1)/dr
            write(*,'(a30,4i6,99e14.5)') 'dt*v_advection_face/dr', &
               steps_used, iter, i, k, dt*v_advection_face(i,k)/dr
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'dt*v_advection_face', &
               steps_used, iter, i, k+1, dt*v_advection_face(i,k+1)
            write(*,'(a30,4i6,99e14.5)') 'dt*v_advection_face', &
               steps_used, iter, i, k, dt*v_advection_face(i,k)
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'v_advection_face', &
               steps_used, iter, i, k+1, v_advection_face(i,k+1)
            write(*,'(a30,4i6,99e14.5)') 'v_advection_face', &
               steps_used, iter, i, k, v_advection_face(i,k)
            write(*,*)
            write(*,'(a30,4i6,99e14.5)') 'GT_face', &
               steps_used, iter, i, k+1, GT_face(i,k+1)
            write(*,'(a30,4i6,99e14.5)') 'GT_face', &
               steps_used, iter, i, k, GT_face(i,k)
            write(*,*)

            write(*,2) 'i_dbg', i_dbg
            write(*,2) 'k_dbg', k_dbg
            write(*,2) 'iter_dbg', iter_dbg
            write(*,*) 'i == i_dbg', i == i_dbg
            write(*,*) 'k == k_dbg', k == k_dbg
            write(*,*) 'iter == iter_dbg', iter == iter_dbg

         end subroutine show_stuff


      end subroutine get_timescale


      subroutine get_dX_dt( &
            nz, nzlo, nzhi, m, nc, X, X_face, C_div_X, GT_face, AD_face, SIG_face, &
            alfa_face, cell_dm, v_advection_face, upwind_limit, dX_dt)
         integer, intent(in) :: nz, nzlo, nzhi, m, nc
         real(dp), dimension(:,:), intent(in) :: X, X_face, C_div_X ! (nc,nz)
         real(dp), intent(in) :: upwind_limit
         real(dp), intent(in), dimension(:) :: alfa_face, cell_dm, AD_face
         real(dp), intent(in), dimension(:,:) :: GT_face, v_advection_face
         real(dp), intent(in), dimension(:,:,:) :: SIG_face
         real(dp), intent(inout), dimension(:,:) :: dX_dt ! (nc,nz)
         integer :: k
!$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,2)
         do k=nzlo,nzhi
            call get1_dX_dt( &
               k, nz, nzlo, nzhi, m, nc, X, X_face, C_div_X, &
               GT_face, AD_face, SIG_face, alfa_face, cell_dm, v_advection_face, &
               upwind_limit, dX_dt(:,k))
         end do
!$OMP END PARALLEL DO
      end subroutine get_dX_dt


      subroutine get1_dX_dt( &
            k, nz, nzlo, nzhi, m, nc, X, X_face, C_div_X, &
            GT_face, AD_face, SIG_face, alfa_face, cell_dm, v_advection_face, &
            upwind_limit, dX_dt)

         integer, intent(in) :: k, nz, nzlo, nzhi, m, nc
         real(dp), dimension(:,:), intent(in) :: X, X_face, C_div_X ! (nc,nz)
         real(dp), intent(in) :: upwind_limit
         real(dp), intent(in), dimension(:) :: alfa_face, cell_dm, AD_face
         real(dp), intent(in), dimension(:,:) :: GT_face, v_advection_face
         real(dp), intent(in), dimension(:,:,:) :: SIG_face
         real(dp), intent(inout), dimension(:) :: dX_dt ! (nc)
         integer :: i, j
         real(dp) :: alfa, beta, c, coeff, dC
         include 'formats'

         dX_dt(:) = 0

         if (k > nzlo) then ! do face k

            alfa = alfa_face(k)
            beta = 1d0 - alfa
            do i=1,nc

               c = GT_face(i,k)
               if (abs(v_advection_face(i,k)) >= upwind_limit) then
                  if (c > 0) then
                     dX_dt(i) = dX_dt(i) - X(i,k)*c
                  else
                     dX_dt(i) = dX_dt(i) - X(i,k-1)*c
                  end if
               else
                  dX_dt(i) = dX_dt(i) - X_face(i,k)*c
               end if

               do j=1,nc

                  c = SIG_face(i,j,k)/max(smallX,X_face(j,k))
                  coeff = X_face(i,k)*c
                  if (j == i) coeff = coeff + AD_face(k)
                  dC = C_div_X(j,k)*X(j,k) - C_div_X(j,k-1)*X(j,k-1)
                  dX_dt(i) = dX_dt(i) - coeff*dC

               end do

            end do

         end if

         if (k < nzhi) then ! do face k+1

            alfa = alfa_face(k+1)
            beta = 1d0 - alfa

            do i=1,nc

               c = GT_face(i,k+1)
               if (abs(v_advection_face(i,k+1)) >= upwind_limit) then
                  if (c > 0) then
                     dX_dt(i) = dX_dt(i) + X(i,k+1)*c
                  else
                     dX_dt(i) = dX_dt(i) + X(i,k)*c
                  end if
               else
                  dX_dt(i) = dX_dt(i) + X_face(i,k+1)*c
               end if

               do j=1,nc

                  c = SIG_face(i,j,k+1)/max(smallX,X_face(j,k+1))
                  coeff = X_face(i,k+1)*c
                  if (j == i) coeff = coeff + AD_face(k+1)
                  dC = C_div_X(j,k+1)*X(j,k+1) - C_div_X(j,k)*X(j,k)
                  dX_dt(i) = dX_dt(i) + coeff*dC

               end do

            end do

         end if

         do i=1,nc
            dX_dt(i) = dX_dt(i)/cell_dm(k)
         end do

      end subroutine get1_dX_dt


      subroutine get_eqn_matrix_entries( &
            nz, nzlo, nzhi, m, nc, jac_only, X, X_face, C_div_X, GT_face, AD_face, SIG_face, &
            alfa_face, cell_dm, v_advection_face, upwind_limit, &
            tiny_X, dt, dX, dX_dt, rhs, em1, e00, ep1)
         integer, intent(in) :: nz, nzlo, nzhi, m, nc
         logical, intent(in) :: jac_only
         real(dp), dimension(:,:), intent(in) :: X, X_face, C_div_X ! (nc,nz)
         real(dp), intent(in) :: upwind_limit, tiny_X, dt
         real(dp), intent(in), dimension(:) :: alfa_face, cell_dm, AD_face
         real(dp), intent(in), dimension(:,:) :: GT_face, v_advection_face
         real(dp), intent(in), dimension(:,:,:) :: SIG_face
         real(dp), intent(in), dimension(:,:) :: dX ! (nc,nz)
         real(dp), intent(inout), dimension(:,:) :: dX_dt ! (nc,nz)
         real(dp), intent(inout), dimension(:,:) :: rhs ! (nc,nz)
         real(dp), intent(inout), dimension(:,:,:) :: em1, e00, ep1 ! (nc,nc,nz)
         integer :: k, j
         ! lhs(i,k) := X(i,k) - (flow(i,k) - flow(i,k-1))*dt/cell_dm(k)
         ! rhs(i,k) := X_prev(i,k)
         ! em1(i,j,k) = d(lhs(i,k))/d(X(j,k-1))
         ! e00(i,j,k) = d(lhs(i,k))/d(X(j,k))
         ! ep1(i,j,k) = d(lhs(i,k))/d(X(j,k+1))
         include 'formats'
!$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,2)
         do k=nzlo,nzhi
            call get1_eqn_matrix_entries( &
               k, nz, nzlo, nzhi, m, nc, jac_only, X, X_face, C_div_X, &
               GT_face, AD_face, SIG_face, alfa_face, cell_dm, v_advection_face, &
               upwind_limit, tiny_X, dt, dX, dX_dt(:,k), rhs, em1, e00, ep1)
         end do
!$OMP END PARALLEL DO
      end subroutine get_eqn_matrix_entries


      subroutine get1_eqn_matrix_entries( &
            k, nz, nzlo, nzhi, m, nc, jac_only, X, X_face, C_div_X, &
            GT_face, AD_face, SIG_face, alfa_face, cell_dm, v_advection_face, &
            upwind_limit, tiny_X, dt, dX, dX_dt, rhs, em1, e00, ep1)

         integer, intent(in) :: k, nz, nzlo, nzhi, m, nc
         logical, intent(in) :: jac_only
         real(dp), dimension(:,:), intent(in) :: X, X_face, C_div_X ! (nc,nz)
         real(dp), intent(in) :: upwind_limit, tiny_X, dt
         real(dp), intent(in), dimension(:) :: alfa_face, cell_dm, AD_face
         real(dp), intent(in), dimension(:,:) :: GT_face, v_advection_face
         real(dp), intent(in), dimension(:,:,:) :: SIG_face
         real(dp), intent(in), dimension(:,:) :: dX ! (nc,nz)
         real(dp), intent(inout), dimension(:) :: dX_dt ! (nc)
         real(dp), intent(inout), dimension(:,:) :: rhs ! (nc,nz)
         real(dp), intent(inout), dimension(:,:,:) :: em1, e00, ep1 ! (nc,nc,nz)
         integer :: i, j, jj
         real(dp) :: alfa, beta, c, coeff, dC_dXj00, dC_dXjm1, dC_dXjp1, &
            dC, dcoeff_dXjm1, dcoeff_dXj00, dcoeff_dXjp1, max_coeff, &
            dt_div_dm

         include 'formats'

         dX_dt(:) = 0; em1(:,:,k) = 0; e00(:,:,k) = 0; ep1(:,:,k) = 0

         if (k > nzlo) then ! do face k

            alfa = alfa_face(k)
            beta = 1d0 - alfa
            do i=1,nc

               c = GT_face(i,k)
               if (abs(v_advection_face(i,k)) >= upwind_limit) then
                  if (c > 0) then
                     dX_dt(i) = dX_dt(i) - X(i,k)*c
                     e00(i,i,k) = e00(i,i,k) - c
                  else
                     dX_dt(i) = dX_dt(i) - X(i,k-1)*c
                     em1(i,i,k) = em1(i,i,k) - c
                  end if
               else
                  dX_dt(i) = dX_dt(i) - X_face(i,k)*c
                  em1(i,i,k) = em1(i,i,k) - beta*c
                  e00(i,i,k) = e00(i,i,k) - alfa*c
               end if

               do j=1,nc

                  c = SIG_face(i,j,k)/max(smallX,X_face(j,k))
                  coeff = X_face(i,k)*c
                  if (j == i) coeff = coeff + AD_face(k)
                  dC = C_div_X(j,k)*X(j,k) - C_div_X(j,k-1)*X(j,k-1)
                  dX_dt(i) = dX_dt(i) - coeff*dC

                  if (j == i) then
                     em1(i,j,k) = em1(i,j,k) - beta*c*dC
                     e00(i,j,k) = e00(i,j,k) - alfa*c*dC
                  end if

                  dC_dXjm1 = -C_div_X(j,k-1)
                  dC_dXj00 = C_div_X(j,k)
                  em1(i,j,k) = em1(i,j,k) - coeff*dC_dXjm1
                  e00(i,j,k) = e00(i,j,k) - coeff*dC_dXj00

                  if (use_dcoeff_dX .and. X_face(j,k) > smallX) then
                     dcoeff_dXjm1 = -beta*coeff/X_face(j,k)
                     dcoeff_dXj00 = -alfa*coeff/X_face(j,k)
                     em1(i,j,k) = em1(i,j,k) - dcoeff_dXjm1*dC
                     e00(i,j,k) = e00(i,j,k) - dcoeff_dXj00*dC
                  end if

               end do

            end do

         end if

         if (k < nzhi) then ! do face k+1

            alfa = alfa_face(k+1)
            beta = 1d0 - alfa

            do i=1,nc

               c = GT_face(i,k+1)
               if (abs(v_advection_face(i,k+1)) >= upwind_limit) then
                  if (c > 0) then
                     dX_dt(i) = dX_dt(i) + X(i,k+1)*c
                     ep1(i,i,k) = ep1(i,i,k) + c
                  else
                     dX_dt(i) = dX_dt(i) + X(i,k)*c
                     e00(i,i,k) = e00(i,i,k) + c
                  end if
               else
                  dX_dt(i) = dX_dt(i) + X_face(i,k+1)*c
                  e00(i,i,k) = e00(i,i,k) + beta*c
                  ep1(i,i,k) = ep1(i,i,k) + alfa*c
               end if

               do j=1,nc

                  c = SIG_face(i,j,k+1)/max(smallX,X_face(j,k+1))
                  coeff = X_face(i,k+1)*c
                  if (j == i) coeff = coeff + AD_face(k+1)
                  dC = C_div_X(j,k+1)*X(j,k+1) - C_div_X(j,k)*X(j,k)
                  dX_dt(i) = dX_dt(i) + coeff*dC

                  if (j == i) then
                     e00(i,j,k) = e00(i,j,k) + beta*c*dC
                     ep1(i,j,k) = ep1(i,j,k) + alfa*c*dC
                  end if

                  dC_dXj00 = -C_div_X(j,k)
                  dC_dXjp1 = C_div_X(j,k+1)
                  e00(i,j,k) = e00(i,j,k) + coeff*dC_dXj00
                  ep1(i,j,k) = ep1(i,j,k) + coeff*dC_dXjp1

                  if (use_dcoeff_dX .and. X_face(j,k+1) > smallX) then
                     dcoeff_dXj00 = -beta*coeff/X_face(j,k+1)
                     dcoeff_dXjp1 = -alfa*coeff/X_face(j,k+1)
                     e00(i,j,k) = e00(i,j,k) + dcoeff_dXj00*dC
                     ep1(i,j,k) = ep1(i,j,k) + dcoeff_dXjp1*dC
                  end if

               end do

            end do

         end if

         if (jac_only) then
            do j=1,nc
               dX_dt(j) = dX_dt(j)/cell_dm(k)
               do i=1,nc
                  em1(i,j,k) = em1(i,j,k)/cell_dm(k)
                  e00(i,j,k) = e00(i,j,k)/cell_dm(k)
                  ep1(i,j,k) = ep1(i,j,k)/cell_dm(k)
               end do
            end do
            return
         end if

         dt_div_dm = dt/cell_dm(k)
         do j=1,nc
            dX_dt(j) = dX_dt(j)/cell_dm(k)
            rhs(j,k) = dX_dt(j)*dt - dX(j,k)
            do i=1,nc
               em1(i,j,k) = -em1(i,j,k)*dt_div_dm
               e00(i,j,k) = -e00(i,j,k)*dt_div_dm
               ep1(i,j,k) = -ep1(i,j,k)*dt_div_dm
            end do
            e00(j,j,k) = e00(j,j,k) + 1d0
         end do

      end subroutine get1_eqn_matrix_entries


      subroutine factor_and_solve_matrix( &
            nz, nzlo, nzhi, nc, n, neqs, &
            lblk1_nz, dblk1_nz, ublk1_nz, ipiv1_nz, del1_nz, &
            lrd, rpar_decsol, lid, ipar_decsol, ierr)

         use mtx_lib, only: bcyclic_dble_decsolblk, block_dble_mv
         integer, intent(in) :: nz, nzlo, nzhi, nc, n, neqs
         real(dp), pointer, dimension(:) :: &
            lblk1_nz, dblk1_nz, ublk1_nz, del1_nz
         integer, pointer :: ipiv1_nz(:)

         integer :: lrd, lid
         integer, pointer :: ipar_decsol(:)
         real(dp), pointer :: rpar_decsol(:)
         integer, intent(out) :: ierr

         integer :: i, j, k, caller_id, ierr2
         integer, pointer :: ipiv1_n(:)
         real(dp), pointer, dimension(:) :: rhs1_n, lblk1_n, dblk1_n, ublk1_n
         real(dp), pointer, dimension(:,:,:) :: lblk, dblk, ublk

         ierr2 = 0
         caller_id = 0

         rhs1_n(1:nc*n) => del1_nz(1+nc*(nzlo-1):nc*nzhi)
         ipiv1_n(1:nc*n) => ipiv1_nz(1+nc*(nzlo-1):nc*nzhi)

         lblk1_n(1:nc*nc*n) => lblk1_nz(1+nc*nc*(nzlo-1):nc*nc*nzhi)
         dblk1_n(1:nc*nc*n) => dblk1_nz(1+nc*nc*(nzlo-1):nc*nc*nzhi)
         ublk1_n(1:nc*nc*n) => ublk1_nz(1+nc*nc*(nzlo-1):nc*nc*nzhi)

         lblk(1:nc,1:nc,1:n) => lblk1_n(1:nc*nc*n)
         dblk(1:nc,1:nc,1:n) => dblk1_n(1:nc*nc*n)
         ublk(1:nc,1:nc,1:n) => ublk1_n(1:nc*nc*n)


         ! factor
         call bcyclic_dble_decsolblk( &
            0, caller_id, nc, n, lblk1_n, dblk1_n, ublk1_n, rhs1_n, ipiv1_n, &
            lrd, rpar_decsol, lid, ipar_decsol, ierr)
         if (ierr /= 0) then
            call bcyclic_dble_decsolblk( &
               2, caller_id, nc, n, lblk1_n, dblk1_n, ublk1_n, rhs1_n, ipiv1_n, &
               lrd, rpar_decsol, lid, ipar_decsol, ierr2)
            return
         end if

         ! solve
         call bcyclic_dble_decsolblk( &
            1, caller_id, nc, n, lblk1_n, dblk1_n, ublk1_n, rhs1_n, ipiv1_n, &
            lrd, rpar_decsol, lid, ipar_decsol, ierr)
         ! deallocate
         call bcyclic_dble_decsolblk( &
            2, caller_id, nc, n, lblk1_n, dblk1_n, ublk1_n, rhs1_n, ipiv1_n, &
            lrd, rpar_decsol, lid, ipar_decsol, ierr2)

      end subroutine factor_and_solve_matrix


      subroutine positivity(nc, nzlo, nzhi, X, eps, del, lambda, num_iters)
         integer, intent(in) :: nc, nzlo, nzhi, num_iters
         real(dp), intent(in) :: eps
         real(dp), intent(in), dimension(:,:) :: X, del
         real(dp), intent(out) :: lambda

         integer :: j, k, bad_j, bad_k
         real(dp) :: alpha, new_x, old_x, dx

         include 'formats'

         lambda = 1d0
         bad_j = 0
         do k=nzlo,nzhi
            do j=1,nc
               old_x = X(j,k)
               if (old_x < 1d-99) cycle
               dx = del(j,k)
               new_x = old_x + dx
               if (new_x >= 0) cycle
               alpha = -(old_x + eps)/dx
               ! so dx*alpha = -old_x - eps,
               ! and therefore old_x + alpha*dx = -eps
               if (alpha < lambda) then
                  lambda = alpha
                  bad_k = k
                  bad_j = j
               end if
            end do
         end do
         if (lambda < 1) lambda = 0.8d0*lambda ! under correct

      end subroutine positivity


      subroutine write_plot_data( &
            nzlo, nzhi, nc, class_chem_id, &
            X1, X2, v)
         use chem_def, only: chem_isos
         use utils_lib, only : mkdir
         integer, intent(in) :: nzlo, nzhi, nc, class_chem_id(:)
         real(dp), intent(in), dimension(:,:) :: X1, X2, v

         integer :: k, j, ierr
         character (len=strlen) :: fname
         integer, parameter :: io = 33

         include 'formats'

         write(*,2) 'nzlo', nzlo
         write(*,2) 'nzhi', nzhi

         call mkdir('plot_data/')
         fname = 'plot_data/diffusion_plot.data'
         ierr = 0
         open(io, file=trim(fname), iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname)
            stop 'write_plot_data'
         end if

         write(io,'(a6)',advance='no') 'k'
         do j=1,nc
            write(io,'(a20)',advance='no') &
               'X1_' // trim(chem_isos% name(class_chem_id(j)))
            write(io,'(a20)',advance='no') &
               'X2_' // trim(chem_isos% name(class_chem_id(j)))
         end do
         do j=1,nc
            write(io,'(a20)',advance='no') &
               'logX1_' // trim(chem_isos% name(class_chem_id(j)))
            write(io,'(a20)',advance='no') &
               'logX2_' // trim(chem_isos% name(class_chem_id(j)))
         end do
         do j=1,nc
            write(io,'(a20)',advance='no') &
               'v_' // trim(chem_isos% name(class_chem_id(j)))
         end do
         write(io,*)

         do k=nzlo,nzhi
            write(io,'(i6)',advance='no') k
            do j=1,nc
               write(io,'(1pe20.12)',advance='no') X1(j,k)
               write(io,'(1pe20.12)',advance='no') X2(j,k)
            end do
            do j=1,nc
               write(io,'(1pe20.12)',advance='no') safe_log10(X1(j,k))
               write(io,'(1pe20.12)',advance='no') safe_log10(X2(j,k))
            end do
            do j=1,nc
               write(io,'(1pe20.12)',advance='no') v(j,k)
            end do
            write(io,*)
         end do

         close(io)

      end subroutine write_plot_data



      ! this routine sets up tables for 4 classes: h, he, o, and fe.
      subroutine get_min_number_of_classes( &
            species, chem_id, class, class_chem_id, class_name)
         use chem_def
         integer, parameter :: nc = diffusion_min_nc
         integer, intent(in) :: species
         integer, intent(in) :: chem_id(:) ! (species)
         integer, intent(out) :: class(:) ! (species)
         integer, intent(out) :: class_chem_id(:) ! (nc)
         character (len=8), intent(out) :: class_name(:) ! (nc)
         real(dp) :: A
         integer :: i, j
         integer, parameter :: c_h = 1, c_he = 2, c_o = 3, c_fe = 4
         class_name(c_h) = 'c_h'
         class_name(c_he) = 'c_he'
         class_name(c_o) = 'c_o'
         class_name(c_fe) = 'c_fe'
         class_chem_id(c_h) = ih1
         class_chem_id(c_he) = ihe4
         class_chem_id(c_o) = io16
         class_chem_id(c_fe) = ife56
         do i=1,species
            A = chem_isos% Z_plus_N(chem_id(i))
            if (A < 3) then
               class(i) = c_h
            else if (A < 12) then
               class(i) = c_he
            else if (A < 20) then
               class(i) = c_o
            else
               class(i) = c_fe
            end if
         end do
      end subroutine get_min_number_of_classes


      ! this routine sets up tables with a separate class for each isotope.
      subroutine get_max_number_of_classes( &
            species, chem_id, class, class_chem_id, class_name)
         use chem_def
         integer, intent(in) :: species, chem_id(:) ! (species)
         integer, intent(out) :: class(:) ! (species)
         integer, intent(out) :: class_chem_id(:) ! (species)
         character (len=8), intent(out) :: class_name(:) ! (species)
         integer :: i, ci
         do i=1,species
            ci = chem_id(i)
            class_name(i) = 'c_' // trim(chem_isos% name(ci))
            class(i) = i
            class_chem_id(i) = ci
         end do
      end subroutine get_max_number_of_classes


      subroutine get_diffusion_classes( &
            nc, species, chem_id, class_chem_id, class_A_cutoff, &
            class, class_name)
         use chem_def
         integer, intent(in) :: nc, species, chem_id(:) ! (species)
         integer, intent(in) :: class_chem_id(:) ! (nc)
         real(dp), intent(in) :: class_A_cutoff(:) ! (nc)
         integer, intent(out) :: class(:) ! (species)
         character (len=8), intent(out) :: class_name(:) ! (nc)
         real(dp) :: A
         integer :: i, j
         do i=1,species
            A = chem_isos% Z_plus_N(chem_id(i))
            class(i) = nc
            do j=1,nc-1
               if (A < class_A_cutoff(j)) then
                  class(i) = j
                  exit
               end if
            end do
         end do
         do j=1,nc
            class_name(j) = 'c_' // trim(chem_isos% name(class_chem_id(j)))
         end do
      end subroutine get_diffusion_classes


      subroutine set_diffusion_classes( &
            nc, species, chem_id, class_chem_id, class_A_max, use_full_net, &
            class, class_name)
         use chem_def
         integer, intent(in) :: nc, species, chem_id(:) ! (species)
         integer, intent(in) :: class_chem_id(:) ! (nc)
         real(dp), intent(in) :: class_A_max(:) ! (nc)
         logical, intent(in) :: use_full_net
         integer, intent(out) :: class(:) ! (species)
         character (len=8), intent(out) :: class_name(:) ! (nc)
         real(dp) :: A
         integer :: i, j
         if(use_full_net) then
            do i=1,species
               class(i) = i
            end do
         else
            do i=1,species
               A = chem_isos% Z(chem_id(i)) + chem_isos% N(chem_id(i))
               class(i) = nc
               do j=1,nc-1
                  if (A <= class_A_max(j)) then
                     class(i) = j
                     exit
                  end if
               end do
            end do
         end if
         do j=1,nc
            class_name(j) = 'c_' // trim(chem_isos% name(class_chem_id(j)))
         end do
      end subroutine set_diffusion_classes


      real(dp) function adjust_timestep( &
            vc_target, vc, vc_old, dt, dt_old, &
            max_timestep_factor, min_timestep_factor) result(dt_next)
         ! H211b "low pass" controller.
         ! Soderlind Wang, J of Computational and Applied Math 185 (2006) 225 – 243.
         use num_def
         real(dp), intent(in) :: vc_target, vc, vc_old, dt, dt_old, &
            max_timestep_factor, min_timestep_factor

         real(dp), parameter :: order = 1d0
         real(dp), parameter :: beta1 = 0.25d0/order
         real(dp), parameter :: beta2 = 0.25d0/order
         real(dp), parameter :: alpha2 = 0.25d0
         real(dp) :: vcratio, vcratio_prev, limtr

         include 'formats'

         dt_next = 1d99

         if (vc_old > 0 .and. dt_old > 0) then ! use 2 values to do "low pass" controller
            vcratio = limiter(vc_target/vc)
            vcratio_prev = limiter(vc_target/vc_old)
            limtr = limiter(pow(vcratio,beta1) * pow(vcratio_prev,beta2) * pow(dt/dt_old,-alpha2))
            dt_next = min(dt_next, dt*limtr)
         else ! no history available, so fall back to the basic controller
            dt_next = min(dt_next, dt*vc_target/vc)
         end if

         if (max_timestep_factor > 0 .and. dt_next > max_timestep_factor*dt) then
            dt_next = max_timestep_factor*dt
         end if

         if (min_timestep_factor > 0 .and. dt_next < min_timestep_factor*dt) then
            dt_next = min_timestep_factor*dt
         end if

         contains

         real(dp) function limiter(x)
            real(dp), intent(in) :: x
            real(dp), parameter :: kappa = 2
            ! for x >= 0 and kappa = 2, limiter value is between 0.07 and 4.14
            ! for x = 1, limiter = 1
            limiter = 1 + kappa*atan((x-1)/kappa)
         end function limiter

      end function adjust_timestep



      end module diffusion

