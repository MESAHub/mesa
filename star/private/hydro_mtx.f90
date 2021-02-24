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
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module hydro_mtx

      use star_private_def
      use const_def

      use num_def

      implicit none

      logical, parameter :: dbg = .false.

      integer, parameter :: ipar_id = 1
      integer, parameter :: ipar_first_call = 2
      integer, parameter :: hydro_lipar = ipar_first_call

      integer, parameter :: rpar_dt = 1
      integer, parameter :: hydro_lrpar = 1


      ! for inspectB debugging
      real(dp), pointer :: debug_previous_data(:,:)

      ! for residual debugging
      real(dp), pointer :: debug_previous_equ_data(:,:)


      contains


      subroutine set_solver_vars(s, nvar, dt, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         s% num_solver_setvars = s% num_solver_setvars + 1
         call set_vars_for_solver(s, nvar, 1, s% nz, dt, ierr)
      end subroutine set_solver_vars


      subroutine set_vars_for_solver(s, nvar, nzlo, nzhi, dt, ierr)
         use const_def, only: secyer, Msun, Lsun, Rsun
         use star_utils, only: set_rmid, set_dm_bar, set_m_and_dm, set_rv_info
         use star_utils, only: current_min_xa_hard_limit, current_sum_xa_hard_limit, &
            lookup_nameofvar
         use hydro_vars, only: set_hydro_vars
         use hydro_rotation, only: set_i_rot, set_omega
         use mix_info, only: get_convection_sigmas
         use chem_def
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nzlo, nzhi
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         logical, parameter :: &
            skip_basic_vars = .true., &
            skip_micro_vars = .false., &
            skip_kap = .false., &
            skip_neu = .false., &
            skip_net = .false., &
            skip_eos = .false., &
            skip_mlt = .false., &
            skip_grads = .true., &
            skip_rotation = .true., &
            skip_m_grav_and_grav = .true., &
            skip_brunt = .true., &
            skip_mixing_info = .true., &
            skip_set_cz_bdy_mass = .true., &
            skip_other_cgrav = .true.
         logical :: do_chem, do_struct, try_again, do_edit_lnR, report_dx
         integer :: i, j, k, kk, klo, khi, i_var, &
            i_lnd, i_lnT, i_lnR, i_lum, i_w, i_v, &
            i_u, i_alpha_RTI, i_ln_cvpv0, i_w_div_wc, i_j_rot, &
            fe56, nvar_chem, species, i_chem1, nz, nvar_hydro
         real(dp), dimension(:, :), pointer :: xh_start, xa_start
         integer :: op_err, kbad, &
            cnt, max_fixes, loc(2), k_lo, k_hi, k_const_mass
         real(dp) :: r2, xavg, du, u00, um1, dx_for_i_var, x_for_i_var, &
            dq_sum, xa_err_norm, d_dxdt_dx, min_xa_hard_limit, sum_xa_hard_limit
         logical :: do_lnd, do_lnT, do_lnR, do_lum, do_w, &
            do_u, do_v, do_alpha_RTI, do_conv_vel, do_w_div_wc, do_j_rot

         include 'formats'

         ierr = 0

         nz = s% nz
         nvar_chem = s% nvar_chem
         species = s% species
         nvar_hydro = s% nvar_hydro
         d_dxdt_dx = s% dVARdot_dVAR
         i_chem1 = s% i_chem1

         xh_start => s% xh_start
         xa_start => s% xa_start
         
         report_dx = &
            s% solver_test_partials_dx_0 > 0d0 .and. &
            s% solver_test_partials_k > 0 .and. &
            s% solver_call_number == s% solver_test_partials_call_number .and. &
            s% solver_test_partials_iter_number == s% solver_iter .and. &
            len_trim(s% solver_test_partials_show_dx_var_name) > 0
            
         if (report_dx) then
            k = s% solver_test_partials_k
            i_var = lookup_nameofvar(s, s% solver_test_partials_show_dx_var_name)            
            if (i_var > 0) then
               if (i_var > nvar_hydro) then
                  dx_for_i_var = s% solver_dx(i_var,k)
                  x_for_i_var = xa_start(i_var-nvar_hydro,k) + s% solver_dx(i_var,k)
               else
                  dx_for_i_var = s% solver_dx(i_var,k)
                  x_for_i_var = xh_start(i_var,k) + s% solver_dx(i_var,k)
               end if
               write(*,3) 'dx, x for var name ' // &
                  trim(s% solver_test_partials_show_dx_var_name), &
                  k, s% solver_iter, dx_for_i_var, x_for_i_var
            end if
         end if


         do_chem = (s% do_burn .or. s% do_mix)
         do_struct = (s% do_struct_hydro .or. s% do_struct_thermo)

         if (dbg .and. .not. skip_mixing_info) write(*,2) 'redo mix info iter', s% solver_iter

         min_xa_hard_limit = current_min_xa_hard_limit(s)
         sum_xa_hard_limit = current_sum_xa_hard_limit(s)

         i_lnd = s% i_lnd
         i_lnT = s% i_lnT
         i_lnR = s% i_lnR
         i_lum = s% i_lum
         i_w = s% i_w
         i_v = s% i_v
         i_u = s% i_u
         i_alpha_RTI = s% i_alpha_RTI
         i_ln_cvpv0 = s% i_ln_cvpv0
         i_w_div_wc = s% i_w_div_wc
         i_j_rot = s% i_j_rot

         do_lnd = i_lnd > 0 .and. i_lnd <= nvar
         do_lnT = i_lnT > 0 .and. i_lnT <= nvar
         do_lnR = i_lnR > 0 .and. i_lnR <= nvar
         do_lum = i_lum > 0 .and. i_lum <= nvar
         do_w = i_w > 0 .and. i_w <= nvar
         do_v = i_v > 0 .and. i_v <= nvar
         do_u = i_u > 0 .and. i_u <= nvar
         do_alpha_RTI = i_alpha_RTI > 0 .and. i_alpha_RTI <= nvar
         do_conv_vel = i_ln_cvpv0 > 0 .and. i_ln_cvpv0 <= nvar
         do_w_div_wc = i_w_div_wc > 0 .and. i_w_div_wc <= nvar
         do_j_rot = i_j_rot > 0 .and. i_j_rot <= nvar

         if (s% trace_k > 0 .and. s% trace_k <= nz) then
            k = s% trace_k
            if (i_lnd /= 0) write(*,3) 'set_vars_for_solver: lnd dx', &
               k, s% solver_iter, xh_start(i_lnd,k), s% solver_dx(i_lnd,k)
         end if

         fe56 = s% net_iso(ife56)
         if (fe56 /= 0) fe56 = i_chem1+fe56-1

         if (nvar > nvar_hydro) then
            do k=1,nz
               do j=1,species
                  s% xa_sub_xa_start(j,k) = s% solver_dx(j+nvar_hydro,k)
                  s% xa(j,k) = xa_start(j,k) + s% solver_dx(j+nvar_hydro,k)
               end do
            end do
            max_fixes = 0 !5
            do cnt=1,max_fixes
               loc = minloc(s% xa(1:species,1:nz))
               j = loc(1)
               k = loc(2)
               if (s% xa(j,k) >= 1d-3*min_xa_hard_limit) exit ! too good to fix
               if (s% xa(j,k) < min_xa_hard_limit) then
                  if (s% report_ierr) then
                     khi = nz
                     do kk=k+1,nz
                        if (s% xa(j,kk) < min_xa_hard_limit) cycle
                        khi = kk-1; exit
                     end do
                     klo = 1
                     do kk=k-1,1,-1
                        if (s% xa(j,kk) < min_xa_hard_limit) cycle
                        klo = kk+1; exit
                     end do
                     do k=klo,khi
                        write(*,2) &
                           'negative ' // trim(chem_isos% name(s% chem_id(j))), &
                           k, s% xa(j,k), xa_start(j,k), s% solver_dx(nvar_hydro+j,k), s% m(k)/Msun
                     end do
                  end if
                  s% retry_message = 'some abundance < min_xa_hard_limit'
                  ierr = -1
                  return
                  exit ! too bad to fix
               end if
               if (k == 1) then
                  k_lo = 1; k_hi = 2
               else if (k == nz) then
                  k_lo = nz-1; k_hi = nz
               else if (s% sig(k) > s% sig(k+1)) then
                  k_lo = k-1; k_hi = k
               else if (s% sig(k+1) > 0) then
                  k_lo = k; k_hi = k+1
               else
                  exit
               end if
               try_again = .true.
               do while (try_again .and. sum(s% xa(j,k_lo:k_hi)*s% dq(k_lo:k_hi)) < 0)
                  try_again = .false.
                  if (k_lo > 1) then
                     if (s% sig(k_lo) > 0) then
                        k_lo = k_lo-1
                        try_again = .true.
                     end if
                  end if
                  if (k_hi < nz) then
                     if (s% sig(k_hi+1) > 0) then
                        k_hi = k_hi+1
                        try_again = .true.
                     end if
                  end if
               end do
               !write(*,3) 'no extend', k_lo, k_hi
               if (.not. try_again) exit
               dq_sum = sum(s% dq(k_lo:k_hi))
               if (s% report_ierr) then
                  write(*,5) 'fix xa(j,k_lo:k_hi)', j, k_lo, k, k_hi, s% xa(j,k), &
                     sum(s% xa(j,k_lo:k_hi)*s% dq(k_lo:k_hi))/dq_sum, dq_sum
                  !stop
                end if
               do j=1,species
                  xavg = sum(s% xa(j,k_lo:k_hi)*s% dq(k_lo:k_hi))/dq_sum
                  do kk=k_lo,k_hi
                     s% xa(j,kk) = xavg
                  end do
               end do
            end do
         end if
         
         if (s% solver_test_partials_k > 0 .and. &
             s% solver_test_partials_k <= nz) then
            k = s% solver_test_partials_k
            !write(*,2) 'set_vars_for_solver iter before', s% solver_iter
            !write(*,2) 'rho', k, exp(s% lnd(k))
            !write(*,2) 'T', k, exp(s% lnT(k))
            !write(*,2) 'r', k, exp(s% lnR(k))
            !write(*,2) 'Et', k, s% Et(k)
            !write(*,2) 'v', k, s% v(k)
            !write(*,*)
         end if

         kbad = 0
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k=1,nz
            op_err = 0
            call set1(k,.false.,op_err)
            if (op_err /= 0) then
               kbad = k; ierr = op_err
            end if
         end do
!$OMP END PARALLEL DO

         if (ierr /= 0) then
            if (s% report_ierr) then
               do k=1,nz ! report the errors sequentially
                  call set1(k,.true.,op_err)
               end do
               write(*,3) 'set_vars_for_solver failed: model, nz', &
                  s% model_number, nz ! kbad depends on num threads
            end if
            return
         end if
         
         if (s% solver_test_partials_k > 0 .and. &
             s% solver_test_partials_k <= nz) then
            k = s% solver_test_partials_k
            !write(*,2) 'set_vars_for_solver iter after', s% solver_iter
            !write(*,2) 'rho', k, exp(s% lnd(k))
            !write(*,2) 'T', k, exp(s% lnT(k))
            !write(*,2) 'T_start', k, s% T_start(k)
            !write(*,2) 'T - T_start', k, exp(s% lnT(k)) - s% T_start(k)
            !write(*,2) 'r', k, exp(s% lnR(k))
            !write(*,2) 'Et', k, s% Et(k)
            !write(*,2) 'v', k, s% v(k)
            !write(*,*)
            !stop 'set_vars_for_solver'
         end if

         if (do_lum) then
            do k=1,nz
               if (is_bad_num(s% L(k))) then
                  if (s% report_ierr) write(*,2) 'set_vars_for_solver L', k, s% L(k), &
                     xh_start(i_lum,k) + s% solver_dx(i_lum,k), &
                     xh_start(i_lum,k), s% solver_dx(i_lum,k)
                  s% retry_message = 'bad num for some L'
                  ierr = -1
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver L', k, s% L(k)
                     stop 'set_vars_for_solver'
                  end if
                  return
                  stop
               end if
            end do
         end if

         if (ierr /= 0) then
            if (s% report_ierr) then

               do k=1,nz
                  if (abs(1d0 - sum(s% xa(:,k))) > 1d-3) then
                     write(*,2) 'set_vars_for_solver: bad xa sum', k, &
                        sum(s% xa(:,k)), sum(xa_start(:,k)), sum(s% solver_dx(i_chem1:nvar,k))
                     write(*,'(51x,a)') 'xa, xa_start+dx, xa_start, dx'
                     do j=1,species
                        write(*,2) trim(chem_isos% name(s% chem_id(j))), k, &
                           s% xa(j,k), xa_start(j,k) + s% solver_dx(i_chem1-1+j,k), &
                           xa_start(j,k), s% solver_dx(i_chem1-1+j,k)
                     end do

                     exit

                  end if
               end do
               write(*,*)

            end if
            return
         end if

         if (do_struct) then
            do_edit_lnR = do_lnR .and. .not. (s% doing_check_partials)
            if (do_edit_lnR) call edit_lnR(s, xh_start, s% solver_dx)
!$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,2)
            do k=1,nz
               if (do_edit_lnR) s% r(k) = exp(s% lnR(k))
               call set_rv_info(s,k)
               ! note: m_grav is held constant during solver iterations
               s% grav(k) = s% cgrav(k)*s% m_grav(k)/(s% r(k)*s% r(k))
            end do
!$OMP END PARALLEL DO
         end if

         if (do_lnR) then
            call set_rmid(s, 1, nz, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*, *) 'set_rmid returned ierr', ierr
               return
            end if
         end if

         if (dbg) write(*, *) 'call set_hydro_vars'
         call set_hydro_vars( &
            s, 1, nz, skip_basic_vars, skip_micro_vars, &
            skip_m_grav_and_grav, skip_eos, skip_net, skip_neu, skip_kap, &
            skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_set_cz_bdy_mass, skip_mlt, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'set_hydro_vars returned ierr', ierr
            return
         end if

         if (s% rotation_flag) then
            ! Moments of inertia are kept fixed as those at the start of the hydro solver
            ! neither omege nor irot enter the equations directly
         !   call set_i_rot(s)
            call set_omega(s, 'hydro_mtx')
         end if


         contains


         subroutine set1(k,report,ierr)
            use chem_def, only: chem_isos
            integer, intent(in) :: k
            logical, intent(in) :: report
            integer, intent(out) :: ierr

            real(dp) :: x(nvar)
            ! setting x = xh_start + dx is necessary because of numerical issues.
            ! we want to ensure that we've calculated the variables using exactly the
            ! same values for x as will be returned as the final result.
            real(dp) :: r2, sum_xa
            integer :: j, i, k_below_just_added
            real(dp) :: del_t, starting_value, alfa, beta, v, theta

            include 'formats'
            ierr = 0
            v = 0
            
            k_below_just_added = 1

            if (do_struct) then

               do j=1,min(nvar, nvar_hydro)
                  x(j) = xh_start(j,k) + s% solver_dx(j,k)
                  !write(*,2) 'new ' // s% nameofvar(j), k, x(j)
               end do

               if (do_lnT) then

                  s% lnT(k) = x(i_lnT)
                  s% dxh_lnT(k) = s% solver_dx(i_lnT,k)
                  if (abs(s% lnT(k) - s% lnT_start(k)) > &
                          ln10*s% hydro_mtx_max_allowed_abs_dlogT .and. &
                       s% min_logT_for_hydro_mtx_max_allowed < &
                        ln10*min(s% lnT(k),s% lnT_start(k))) then
                     if (report) &
                        write(*,4) 'hydro_mtx: change too large, dlogT, logT, logT_start', &
                           s% model_number, k, s% solver_iter, &
                           (s% lnT(k) - s% lnT_start(k))/ln10, &
                           s% lnT(k)/ln10, s% lnT_start(k)/ln10
                     write(s% retry_message, *) 'abs(dlogT) > hydro_mtx_max_allowed_abs_dlogT', k
                     ierr = -1
                     return
                  end if
                  if (s% lnT(k) > ln10*s% hydro_mtx_max_allowed_logT .and. &
                       s% min_logT_for_hydro_mtx_max_allowed < &
                        ln10*min(s% lnT(k),s% lnT_start(k))) then
                     if (report) &
                        write(*,4) 'hydro_mtx: logT too large', &
                           s% model_number, k, s% solver_iter, &
                           s% lnT(k)/ln10, s% lnT_start(k)/ln10
                     write(s% retry_message, *) 'logT > hydro_mtx_max_allowed_logT', k
                     ierr = -1
                     return
                  end if
                  if (s% lnT(k) < ln10*s% hydro_mtx_min_allowed_logT) then
                     if (report) &
                        write(*,4) 'hydro_mtx: logT too small', &
                           s% model_number, k, s% solver_iter, &
                           s% lnT(k)/ln10, s% lnT_start(k)/ln10
                     write(s% retry_message, *) 'logT < hydro_mtx_min_allowed_logT', k
                     ierr = -1
                     return
                  end if
                  s% T(k) = exp(s% lnT(k))
                  if (is_bad_num(s% T(k))) then
                     s% retry_message = 'bad num for T'
                     if (report) write(*,2) 'bad num T', k, s% T(k)
                     if (s% stop_for_bad_nums) then
                        write(*,2) 'set_vars_for_solver T', k, s% T(k)
                        stop 'set_vars_for_solver'
                     end if
                     ierr = -1
                  end if

               end if

               if (do_lum) then
                  s% L(k) = x(i_lum)
                  if (is_bad_num(s% L(k))) then
                     s% retry_message = 'bad num for L'
                     if (report) write(*,2) 'bad num L', k, s% L(k)
                     ierr = -1
                     if (s% stop_for_bad_nums) then
                        write(*,2) 'set_vars_for_solver L', k, s% L(k)
                        stop 'set_vars_for_solver'
                     end if
                  end if

               end if

               if (do_w) then
                  s% w(k) = max(x(i_w),0d0)
                  s% dxh_w(k) = s% solver_dx(i_w,k)
                  if (s% w(k) < 0d0 .or. is_bad_num(s% w(k))) then
                     s% retry_message = 'bad num for et'
                     if (report) write(*,2) 'bad num et', k, s% w(k)
                     ierr = -1
                     if (s% stop_for_bad_nums) then
!$omp critical (set_vars_for_solver_crit1)
                        write(*,2) 'set_vars_for_solver et', k, s% w(k)
                        write(*,2) 'set_vars_for_solver et_start', k, s% w_start(k)
                        write(*,2) 'set_vars_for_solver xh_start', k, xh_start(i_w,k)
                        write(*,2) 'set_vars_for_solver dx', k, s% solver_dx(i_w,k)
                        stop 'set_vars_for_solver'
!$omp end critical (set_vars_for_solver_crit1)
                     end if
                  end if

               end if

               if (s% do_struct_hydro) then

                  if (do_v) then
                     s% v(k) = x(i_v)
                     if (is_bad_num(s% v(k))) then
                        s% retry_message = 'bad num for v'
                        if (report) write(*,2) 'bad num v', k, s% v(k)
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'set_vars_for_solver v', k, s% v(k)
                           stop 'set_vars_for_solver'
                        end if
                        ierr = -1
                     end if
                  end if

                  if (do_u) then
                     s% u(k) = x(i_u)
                     if (is_bad_num(s% u(k))) then
                        s% retry_message = 'bad num for u'
                        if (report) write(*,2) 'bad num u', k, s% u(k)
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'set_vars_for_solver u', k, s% u(k)
                           stop 'set_vars_for_solver'
                        end if
                        ierr = -1
                     end if
                  end if

                  if (do_alpha_RTI) then
                     s% alpha_RTI(k) = max(0d0, x(i_alpha_RTI))
                     if (is_bad_num(s% alpha_RTI(k))) then
                        s% retry_message = 'bad num for alpha_RTI'
                        if (report) write(*,2) 'bad num alpha_RTI', k, s% alpha_RTI(k)
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'set_vars_for_solver alpha_RTI', k, s% alpha_RTI(k)
                           stop 'set_vars_for_solver'
                        end if
                        ierr = -1
                     end if
                     if (s% alpha_RTI(k) < 0d0) then
                        s% alpha_RTI(k) = 0d0
                        x(i_alpha_RTI) = 0d0
                     end if
                  end if

                  if (do_conv_vel) then
                     s% dxh_ln_cvpv0(k) = s% solver_dx(i_ln_cvpv0,k)
                     s% conv_vel(k) = max(0d0, exp(x(i_ln_cvpv0))-s% conv_vel_v0)
                     if (s% conv_vel(k) > 1d90 .or. is_bad_num(s% conv_vel(k))) then
                        s% retry_message = 'bad num for conv_vel'
                        if (report) write(*,2) 'bad num conv_vel', k, s% conv_vel(k)
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'set_vars_for_solver conv_vel', k, s% conv_vel(k)
                           stop 'set_vars_for_solver'
                        end if
                        ierr = -1
                     end if
                     if (s% conv_vel(k) < 0d0) then
                        s% conv_vel(k) = 0d0
                        x(i_ln_cvpv0) = log(s% conv_vel_v0)
                        s% dxh_ln_cvpv0(k) = x(i_ln_cvpv0) - xh_start(i_ln_cvpv0,k)
                     end if
                  end if

                  if (do_w_div_wc) then
                     s% w_div_w_crit_roche(k) = x(i_w_div_wc)
                     if (s% w_div_w_crit_roche(k) > 0.99d0) then
                        s% w_div_w_crit_roche(k) = 0.99d0
                     end if
                     if (s% w_div_w_crit_roche(k) < -0.99d0) then
                        s% w_div_w_crit_roche(k) = -0.99d0
                     end if
                     if (is_bad_num(s% w_div_w_crit_roche(k))) then
                        s% retry_message = 'bad num for w_div_w_crit_roche'
                        if (report) write(*,2) 'bad num w_div_w_crit_roche', k, s% w_div_w_crit_roche(k)
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'set_vars_for_solver w_div_w_crit_roche', k, s% w_div_w_crit_roche(k)
                           stop 'set_vars_for_solver'
                        end if
                        ierr = -1
                     end if
                  end if

                  if (do_j_rot) then
                     s% j_rot(k) = x(i_j_rot)
                     if (is_bad_num(s% j_rot(k))) then
                        s% retry_message = 'bad num for j_rot'
                        if (report) write(*,2) 'bad num j_rot', k, s% j_rot(k)
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'set_vars_for_solver j_rot', k, s% j_rot(k)
                           stop 'set_vars_for_solver'
                        end if
                        ierr = -1
                     end if
                  end if

                  if (do_lnR) then
                     s% lnR(k) = x(i_lnR)
                     s% dxh_lnR(k) = s% solver_dx(i_lnR,k)
                     if (is_bad_num(s% lnR(k))) then
                        s% retry_message = 'bad num for lnR'
                        if (report) write(*,2) 'bad num lnR', k, s% lnR(k)
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'set_vars_for_solver lnR', k, s% lnR(k)
                           stop 'set_vars_for_solver'
                        end if
                        ierr = -1
                     end if
                     s% r(k) = exp(s% lnR(k))
                  end if

                  if (do_lnd) then
                     s% lnd(k) = x(i_lnd)
                     s% dxh_lnd(k) = s% solver_dx(i_lnd,k)
                     if (s% lnd(k) < ln10*s% hydro_mtx_min_allowed_logRho) then
                        write(s% retry_message, *) 'logRho < hydro_mtx_min_allowed_logRho', k
                        if (report) &
                           write(*,4) 'hydro_mtx: logRho too small', &
                              s% model_number, k, s% solver_iter, &
                              s% lnd(k)/ln10, s% lnd_start(k)/ln10
                        ierr = -1
                        return
                     end if
                     if (s% lnd(k) > ln10*s% hydro_mtx_max_allowed_logRho .and. &
                          s% min_logT_for_hydro_mtx_max_allowed < &
                           ln10*min(s% lnT(k),s% lnT_start(k))) then
                        write(s% retry_message, *) 'logRho > hydro_mtx_max_allowed_logRho', k
                        if (s% report_ierr .or. report) &
                           write(*,4) 'hydro_mtx: logRho too large', &
                              s% model_number, k, s% solver_iter, &
                              s% lnd(k)/ln10, s% lnd_start(k)/ln10
                        ierr = -1
                        return
                     end if
                     if (abs(s% lnd(k) - s% lnd_start(k)) > &
                           ln10*s% hydro_mtx_max_allowed_abs_dlogRho .and. &
                          s% min_logT_for_hydro_mtx_max_allowed < &
                           ln10*min(s% lnT(k),s% lnT_start(k))) then
                        write(s% retry_message, *) 'abs(dlogRho) > hydro_mtx_max_allowed_abs_dlogRho', k
                        if (s% report_ierr .or. report) &
                           write(*,4) &
                              'hydro_mtx: dlogRho, logRho, logRho_start', &
                              s% model_number, k, s% solver_iter, &
                              (s% lnd(k) - s% lnd_start(k))/ln10, &
                              s% lnd(k)/ln10, s% lnd_start(k)/ln10
                        ierr = -1
                        return
                     end if
                     s% rho(k) = exp(s% lnd(k))
                     if (is_bad_num(s% rho(k))) then
                        write(s% retry_message, *) 'bad num for rho', k
                        if (report) write(*,2) 'bad num rho', k, s% rho(k)
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'set_vars_for_solver rho', k, s% rho(k)
                           stop 'set_vars_for_solver'
                        end if
                        ierr = -1
                     end if
                  end if

               end if

               if (k == s% trace_k) then
                  if (i_lnd /= 0) &
                     write(*,4) 'hydro_mtx: lgd', &
                        k, s% solver_iter, s% model_number, &
                        s% lnd(k)/ln10, xh_start(i_lnd,k), s% solver_dx(i_lnd,k)
                  if (i_lnT /= 0) &
                     write(*,4) 'hydro_mtx: lgT', k, s% solver_iter, &
                        s% model_number, s% lnT(k)/ln10, xh_start(i_lnT,k), s% solver_dx(i_lnT,k)
                  if (i_lum /= 0) &
                     write(*,4) 'hydro_mtx: L', k, s% solver_iter, &
                        s% model_number, s% L(k), xh_start(i_lum,k), s% solver_dx(i_lum,k)
                  write(*,4) 'hydro_mtx: lgR', k, s% solver_iter, &
                        s% model_number, s% lnR(k)/ln10, xh_start(i_lnR,k), s% solver_dx(i_lnR,k)
                  if (i_v /= 0) &
                     write(*,4) 'hydro_mtx: v', k, s% solver_iter, &
                        s% model_number, s% v(k), xh_start(i_v,k), s% solver_dx(i_v,k)
                  if (i_u /= 0) &
                     write(*,4) 'hydro_mtx: u', k, s% solver_iter, &
                        s% model_number, s% u(k), xh_start(i_u,k), s% solver_dx(i_u,k)
               end if

               ! set time derivatives at constant q -- only need the ones for eps_grav
               if (dt == 0) then
                  s% dlnT_dt_const_q(k) = 0
                  if (s% do_struct_hydro) then
                     s% dlnd_dt_const_q(k) = 0
                     if (i_ln_cvpv0 /= 0) s% dln_cvpv0_dt_const_q(k) = 0
                  end if
               else if (k < s% k_below_const_q) then
                  ! use dx to get better accuracy
                  if (i_lnT /= 0) s% dlnT_dt_const_q(k) = s% solver_dx(i_lnT,k)*d_dxdt_dx
                  if (s% do_struct_hydro) then
                     if (i_lnd /= 0) s% dlnd_dt_const_q(k) = s% solver_dx(i_lnd,k)*d_dxdt_dx
                     if (i_ln_cvpv0 /= 0) &
                        s% dln_cvpv0_dt_const_q(k) = s% solver_dx(i_ln_cvpv0,k)*d_dxdt_dx
                  end if
               else
                  if (i_lnT /= 0) s% dlnT_dt_const_q(k) = &
                     (x(i_lnT) - s% lnT_for_d_dt_const_q(k))*d_dxdt_dx
                  if (s% do_struct_hydro) then
                     if (i_lnd /= 0) &
                        s% dlnd_dt_const_q(k) = &
                           (x(i_lnd) - s% lnd_for_d_dt_const_q(k))*d_dxdt_dx
                     if (i_ln_cvpv0 /= 0) &
                        s% dln_cvpv0_dt_const_q(k) = &
                           (x(i_ln_cvpv0) - s% ln_cvpv0_for_d_dt_const_q(k))*d_dxdt_dx
                  end if
               end if

               ! set time derivatives at constant mass
               if (dt == 0) then

                  s% dlnT_dt(k) = 0
                  s% dw_dt(k) = 0
                  if (s% do_struct_hydro) then
                     s% dlnd_dt(k) = 0
                     s% dlnR_dt(k) = 0
                     if (i_v /= 0) s% dv_dt(k) = 0
                     if (i_u /= 0) s% du_dt(k) = 0
                     if (i_alpha_RTI /= 0) s% dalpha_RTI_dt(k) = 0
                     if (i_ln_cvpv0 /= 0) s% dln_cvpv0_dt(k) = 0
                     if (i_j_rot /= 0) s% dj_rot_dt(k) = 0
                  end if

               else if (k >= s% k_const_mass) then
                  ! use dx to get better accuracy

                  if (do_lnT) s% dlnT_dt(k) = s% solver_dx(i_lnT,k)*d_dxdt_dx
                  if (do_w) s% dw_dt(k) = s% solver_dx(i_w,k)*d_dxdt_dx

                  if (s% do_struct_hydro) then
                     if (do_lnd) s% dlnd_dt(k) = s% solver_dx(i_lnd,k)*d_dxdt_dx
                     if (do_v) s% dv_dt(k) = s% solver_dx(i_v,k)*d_dxdt_dx
                     if (do_u) s% du_dt(k) = s% solver_dx(i_u,k)*d_dxdt_dx
                     if (do_alpha_RTI) s% dalpha_RTI_dt(k) = s% solver_dx(i_alpha_RTI,k)*d_dxdt_dx
                     if (do_lnR) s% dlnR_dt(k) = s% solver_dx(i_lnR,k)*d_dxdt_dx
                  end if

                  if (k == s% trace_k) then
                     write(*,4) 'd_dxdt_dx', k, s% solver_iter, &
                        s% model_number, d_dxdt_dx
                     if (do_lnT) write(*,4) 's% dlnT_dt(k)', k, s% solver_iter, &
                        s% model_number, s% dlnT_dt(k)
                     write(*,2) 'd_dxdt_dx', k, d_dxdt_dx
                     write(*,*)
                  end if

               else if (k >= k_below_just_added) then

                  if (do_lnT) &
                     s% dlnT_dt(k) = (x(i_lnT) - s% lnT_for_d_dt_const_m(k))*d_dxdt_dx
                  if (do_w) &
                     s% dw_dt(k) = (x(i_w) - s% w_for_d_dt_const_m(k))*d_dxdt_dx
                  if (s% do_struct_hydro) then
                     s% dlnR_dt(k) = (x(i_lnR) - s% lnR_for_d_dt_const_m(k))*d_dxdt_dx
                     if (do_lnd) &
                        s% dlnd_dt(k) = (x(i_lnd) - s% lnd_for_d_dt_const_m(k))*d_dxdt_dx
                     if (do_v) s% dv_dt(k) = &
                        (x(i_v) - s% v_for_d_dt_const_m(k))*d_dxdt_dx
                     if (do_u) s% du_dt(k) = &
                        (x(i_u) - s% u_for_d_dt_const_m(k))*d_dxdt_dx
                     if (do_alpha_RTI) s% dalpha_RTI_dt(k) = &
                        (x(i_alpha_RTI) - s% alpha_RTI_for_d_dt_const_m(k))*d_dxdt_dx
                  end if

                  if (k == s% trace_k) then
                     write(*,4) 's% dlnT_dt(k)', k, s% solver_iter, &
                        s% model_number, s% dlnT_dt(k)
                     write(*,4) 'x(i_lnT)', k, s% solver_iter, &
                        s% model_number, x(i_lnT)
                     write(*,4) 's% lnT_for_d_dt_const_m(k)', k, &
                        s% solver_iter, s% model_number, s% lnT_for_d_dt_const_m(k)
                     write(*,4) 'xh_start(i_lnT,1)', 1, s% solver_iter, &
                        s% model_number, xh_start(i_lnT,1)
                     write(*,2) 'd_dxdt_dx', k, d_dxdt_dx
                     write(*,*)
                  end if

               else ! k < s% k_below_just_added, so new surface cell

                  s% dlnT_dt(k) = 0
                  s% dw_dt(k) = 0
                  if (s% do_struct_hydro) then
                     s% dlnR_dt(k) = 0
                     s% dlnd_dt(k) = 0
                     s% dv_dt(k) = 0
                  end if

               end if

               if (s% do_struct_hydro .and. do_conv_vel) then
                  if (k >= s% k_const_mass) then
                     s% dln_cvpv0_dt(k) = s% solver_dx(i_ln_cvpv0,k)*d_dxdt_dx
                  else
                     s% dln_cvpv0_dt(k) = &
                        (x(i_ln_cvpv0) - s% ln_cvpv0_for_d_dt_const_m(k))*d_dxdt_dx
                  end if
               end if

            end if

            if (s% do_struct_hydro .and. do_j_rot) then
               s% dj_rot_dt(k) = s% solver_dx(i_j_rot,k)*d_dxdt_dx
            end if

            if (do_chem) &
               call check1_chem( &
                  s, k, min_xa_hard_limit, sum_xa_hard_limit, report, ierr)

         end subroutine set1


      end subroutine set_vars_for_solver


      subroutine check1_chem( &
            s, k, min_xa_hard_limit, sum_xa_hard_limit, report, ierr)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: min_xa_hard_limit, sum_xa_hard_limit
         logical, intent(in) :: report
         integer, intent(out) :: ierr

         integer :: j, species, jmax
         real(dp) :: sum_xa, xsum
         logical :: okay

         include 'formats'

         ierr = 0
         species = s% species
         okay = .true.
         if (min_xa_hard_limit > -1d50) then
            do j=1,species
               if (s% xa(j,k) < min_xa_hard_limit) then
                  s% retry_message = 'some xa < min_xa_hard_limit'
                  if (report .or. s% report_bad_negative_xa) then
                     write(*,3) &
                        'bad negative xa, min_xa_hard_limit, sig, logT ' // &
                        trim(chem_isos% name(s% chem_id(j))), j, k, &
                        s% xa(j,k), min_xa_hard_limit, s% sig(k), s% lnT(k)/ln10
                     okay = .false.
                  end if
                  ierr = -1
!$omp critical (hydro_mtx_crit1)
                  s% why_Tlim = Tlim_neg_X
!$omp end critical (hydro_mtx_crit1)
                  if (.not. report) return
               end if
            end do
         end if

         do j=1,species
            s% xa(j,k) = max(0d0, min(1d0, s% xa(j,k)))
         end do

         sum_xa = sum(s% xa(1:species,k))
         if (is_bad_num(sum_xa)) then
            if (report) then
               write(*,'(a60,i8,99f20.10)') 'bad num sum X', k, s% m(k)/Msun, sum_xa
            end if
            ierr = -1
            s% retry_message = 'bad num for mass fractions'
            if (s% stop_for_bad_nums) then
               write(*,2) 'set_vars_for_solver sum_xa', k, sum_xa
               stop 'set_vars_for_solver'
            end if
!$omp critical (hydro_mtx_crit2)
            if (s% why_Tlim /= Tlim_neg_X) s% why_Tlim = Tlim_bad_Xsum
!$omp end critical (hydro_mtx_crit2)
            if (.not. report) return
         end if
         if (abs(sum_xa - 1d0) > sum_xa_hard_limit) then
            s% retry_message = 'mass fractions > sum_xa_hard_limit'
            if (report) then
               write(*,2) &
                  'bad sumX', k, &
                  sum_xa - 1d0, sum_xa_hard_limit, &
                  sum(s% xa_start(1:species,k)), &
                  sum_xa - sum(s% xa_start(1:species,k))
            end if
            ierr = -1
!$omp critical (hydro_mtx_crit3)
            if (s% why_Tlim /= Tlim_neg_X) s% why_Tlim = Tlim_bad_Xsum
!$omp end critical (hydro_mtx_crit3)
            okay = .false.
            if (.not. report) return
         end if

         if (abs(sum_xa - 1d0) > 1d-12) then
            !jmax = maxloc(s% xa(1:species,k),dim=1)
            !xsum = sum(s% xa(1:species,k)) - s% xa(jmax,k)
            !if (1d0 > xsum) then
            !   s% xa(jmax,k) = 1d0 - xsum
            !else
               do j=1,species
                  s% xa(j,k) = s% xa(j,k)/sum_xa
               end do
            !end if
         end if

         if (s% xa_clip_limit > 0) then
            do j=1,species
               if (s% xa(j,k) < s% xa_clip_limit) s% xa(j,k) = 0d0
            end do
         end if

      end subroutine check1_chem


      subroutine dump_struct(s)
         type (star_info), pointer :: s
         integer :: k, j, i

         include 'formats'

         do k=1,s% nz
            write(*,2) 'dq', k, s% dq(k)
            write(*,2) 'm', k, s% m(k)
            write(*,2) 'T', k, s% T(k)
            write(*,2) 'rho', k, s% rho(k)
            write(*,2) 'Pgas', k, s% Pgas(k)
            write(*,2) 'L', k, s% L(k)
            write(*,2) 'r', k, s% r(k)
            write(*,2) 'grada', k, s% grada(k)
            write(*,2) 'opacity', k, s% opacity(k)
            write(*,2) 'd_opacity_dlnd', k, s% d_opacity_dlnd(k)
            write(*,2) 'd_opacity_dlnT', k, s% d_opacity_dlnT(k)
            write(*,2) 'eps_nuc', k, s% eps_nuc(k)
            write(*,2) 'd_epsnuc_dlnd', k, s% d_epsnuc_dlnd(k)
            write(*,2) 'd_epsnuc_dlnT', k, s% d_epsnuc_dlnT(k)
            write(*,2) 'non_nuc_neu', k, s% non_nuc_neu(k)
            write(*,2) 'd_nonnucneu_dlnd', k, s% d_nonnucneu_dlnd(k)
            write(*,2) 'd_nonnucneu_dlnT', k, s% d_nonnucneu_dlnT(k)
            write(*,2) 'eps_grav', k, s% eps_grav(k)
            write(*,2) 'gradT', k, s% gradT(k)
            !write(*,2) '', k, s% (k)
            do j=1,s% species
               write(*,3) 'xa(j,k)', j, k, s% xa(j,k)
            end do
         end do


      end subroutine dump_struct


      subroutine edit_lnR(s, xh_start, dx)
         ! uses mass and density to set radius
         type (star_info), pointer :: s
         real(dp), dimension(:, :) :: xh_start
         real(dp), dimension(:, :) :: dx
         real(dp) :: vol00, volp1, cell_vol
         integer :: k, nz
         include 'formats'
         vol00 = four_thirds_pi*s% R_center*s% R_center*s% R_center
         nz = s% nz
         do k=nz, 1, -1
            volp1 = vol00
            cell_vol = s% dm(k)/s% rho(k)
            vol00 = volp1 + cell_vol
            s% lnR(k) = log(vol00/four_thirds_pi)/3
            dx(s% i_lnR,k) = s% lnR(k) - xh_start(s% i_lnR,k)
            if (k >= s% k_below_just_added) &
               s% dlnR_dt(k) = &
                  (s% lnR(k) - s% lnR_for_d_dt_const_m(k))*s% dVARDOT_dVAR
         end do
         call edit_dlnR_dt_above_k_below_just_added(s, xh_start)
      end subroutine edit_lnR


      subroutine edit_dlnR_dt_above_k_below_just_added(s, xh_start)
         type (star_info), pointer :: s
         real(dp), dimension(:, :) :: xh_start
         integer :: k, k_below_just_added
         real(dp) :: lnR_start
         k_below_just_added = s% k_below_just_added
         if (k_below_just_added == 1) return
         lnR_start = xh_start(s% i_lnR,1)
         do k = 1, k_below_just_added - 1
            s% dlnR_dt(k) = 0d0
         end do
      end subroutine edit_dlnR_dt_above_k_below_just_added


      subroutine enter_setmatrix(s, &
            nvar, xder, need_solver_to_eval_jacobian, &
            ldA, A1, ierr)
         use mtx_def, only: lapack
         use rsp_def, only: ABB, LD_ABB, NV, MAX_NZN
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         real(dp), pointer, dimension(:,:) :: xder ! (nvar, nz)
         logical, intent(out) :: need_solver_to_eval_jacobian
         integer, intent(in) :: ldA ! leading dimension of A
         real(dp), pointer, dimension(:) :: A1
         integer, intent(out) :: ierr

         real(dp), pointer, dimension(:,:) :: A ! (ldA, neqns)
         integer :: i, j, k, cnt, i_lnR, neqns, nz, nzlo, nzhi
         real(dp) :: dt, lnR00, lnRm1, lnRp1, dlnR_prev, ddx, ddx_limit, &
            epsder_struct, epsder_chem
         integer :: id, nvar_hydro
         logical :: dbg_enter_setmatrix, do_chem
         real(dp), pointer :: blk3(:, :, :, :)

         include 'formats'

         dbg_enter_setmatrix = dbg

         ierr = 0
         nz = s% nz
         neqns = nvar*nz
         
         if (dbg_enter_setmatrix) write(*, '(/,/,/,/,/,/,a)') 'enter_setmatrix'

         if (s% model_number == 1) then
            s% num_solver_iterations = s% num_solver_iterations + 1
            if (s% num_solver_iterations > 60 .and. &
                  mod(s% num_solver_iterations,10) == 0) &
               write(*,*) 'first model is slow to converge: num tries', &
                  s% num_solver_iterations
         end if

         dt = s% dt

         do_chem = (s% do_burn .or. s% do_mix)

         s% jacobian(1:ldA,1:neqns) => A1(1:ldA*neqns)
         A(1:ldA,1:neqns) => A1(1:ldA*neqns)

         do i=1,neqns
            do j=1,lda
               A(j,i) = 0d0
            end do
         end do
         i = nvar*nvar*nz
         if (size(A1,dim=1) < 3*i) then
            write(*,*) 'enter_setmatrix: size(A1,dim=1) < 3*i', size(A1,dim=1), 3*i
            ierr = -1
            return
         end if
         s% ublk(1:nvar,1:nvar,1:nz) => A1(1:i)
         s% dblk(1:nvar,1:nvar,1:nz) => A1(i+1:2*i)
         s% lblk(1:nvar,1:nvar,1:nz) => A1(2*i+1:3*i)

         if (dbg_enter_setmatrix) &
            write(*, *) 'call eval_partials with doing_check_partials = .false.'
         call eval_partials(s, nvar, ierr)
         if (ierr /= 0) return

         call s% other_after_enter_setmatrix(s% id,ierr)

         if (dbg_enter_setmatrix) write(*, *) 'finished enter_setmatrix'
         need_solver_to_eval_jacobian = .false.

      end subroutine enter_setmatrix


      subroutine eval_partials(s, nvar, ierr)
         use hydro_eqns, only: eval_equ
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         logical, parameter :: skip_partials = .false.
         ierr = 0
         call eval_equ(s, nvar, skip_partials, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'eval_partials: eval_equ returned ierr', ierr
            return
         end if
      end subroutine eval_partials


      end module hydro_mtx

