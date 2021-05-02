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
         logical :: do_chem, try_again, do_edit_lnR, report_dx
         integer :: i, j, k, kk, klo, khi, i_var, &
            i_lnd, i_lnT, i_lnR, i_lum, i_w, i_Hp, i_v, &
            i_u, i_alpha_RTI, i_ln_cvpv0, i_w_div_wc, i_j_rot, &
            fe56, nvar_chem, species, nz, nvar_hydro
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
         d_dxdt_dx = 1d0/s% dt

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

         if (dbg .and. .not. skip_mixing_info) write(*,2) 'redo mix info iter', s% solver_iter

         min_xa_hard_limit = current_min_xa_hard_limit(s)
         sum_xa_hard_limit = current_sum_xa_hard_limit(s)

         i_lnd = s% i_lnd
         i_lnT = s% i_lnT
         i_lnR = s% i_lnR
         i_lum = s% i_lum
         i_w = s% i_w
         i_Hp = s% i_Hp
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

         fe56 = s% net_iso(ife56)
         if (fe56 /= 0) fe56 = nvar_hydro+fe56

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
                        sum(s% xa(:,k)), sum(xa_start(:,k)), sum(s% solver_dx(nvar_hydro+1:nvar,k))
                     write(*,'(51x,a)') 'xa, xa_start+dx, xa_start, dx'
                     do j=1,species
                        write(*,2) trim(chem_isos% name(s% chem_id(j))), k, &
                           s% xa(j,k), xa_start(j,k) + s% solver_dx(nvar_hydro+j,k), &
                           xa_start(j,k), s% solver_dx(nvar_hydro+j,k)
                     end do

                     exit

                  end if
               end do
               write(*,*)

            end if
            return
         end if

         do_edit_lnR = do_lnR .and. (.not. s% doing_check_partials)
         if (do_edit_lnR) call edit_lnR(s, xh_start, s% solver_dx)
         do k=1,nz
            if (do_edit_lnR) s% r(k) = exp(s% lnR(k))
            call set_rv_info(s,k)
            ! note: m_grav is held constant during solver iterations
            s% grav(k) = s% cgrav(k)*s% m_grav(k)/(s% r(k)*s% r(k))
         end do

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

            do j=1,min(nvar, nvar_hydro)
               x(j) = xh_start(j,k) + s% solver_dx(j,k)
               !write(*,2) 'new ' // s% nameofvar(j), k, x(j)
            end do

            if (do_lnT) then
               
               s% lnT(k) = x(i_lnT)
               s% T(k) = exp(s% lnT(k))
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
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver T', k, s% T(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num T', k, s% T(k)
                  ierr = -1
               end if

            end if

            if (do_lum) then
               s% L(k) = x(i_lum)
               if (is_bad_num(s% L(k))) then
                  s% retry_message = 'bad num for L'
                  ierr = -1
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver L', k, s% L(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num L', k, s% L(k)
               end if
            end if

            if (do_w) then
               s% w(k) = x(i_w)
               if (s% w(k) < 0d0) s% w(k) = s% RSP2_w_fix_if_neg
               if (is_bad_num(s% w(k))) then
                  s% retry_message = 'bad num for w'
                  ierr = -1
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver w', k, s% w(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num w', k, s% w(k)
               end if
               s% Hp_face(k) = x(i_Hp)
               if (is_bad_num(s% Hp_face(k))) then
                  s% retry_message = 'bad num for Hp_face'
                  ierr = -1
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver Hp_face', k, s% Hp_face(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num Hp_face', k, s% Hp_face(k)
               end if
            end if
            
            if (do_v) then
               s% v(k) = x(i_v)
               s% dxh_v(k) = s% solver_dx(i_v,k)
               if (is_bad_num(s% v(k))) then
                  s% retry_message = 'bad num for v'
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver v', k, s% v(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num v', k, s% v(k)
                  ierr = -1
               end if
            end if

            if (do_u) then
               s% u(k) = x(i_u)
               s% dxh_u(k) = s% solver_dx(i_u,k)
               if (is_bad_num(s% u(k))) then
                  s% retry_message = 'bad num for u'
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver u', k, s% u(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num u', k, s% u(k)
                  ierr = -1
               end if
            end if

            if (do_alpha_RTI) then
               s% alpha_RTI(k) = max(0d0, x(i_alpha_RTI))
               s% dxh_alpha_RTI(k) = s% solver_dx(i_alpha_RTI,k)
               if (is_bad_num(s% alpha_RTI(k))) then
                  s% retry_message = 'bad num for alpha_RTI'
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver alpha_RTI', k, s% alpha_RTI(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num alpha_RTI', k, s% alpha_RTI(k)
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
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver conv_vel', k, s% conv_vel(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num conv_vel', k, s% conv_vel(k)
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
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver w_div_w_crit_roche', k, s% w_div_w_crit_roche(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num w_div_w_crit_roche', k, s% w_div_w_crit_roche(k)
                  ierr = -1
               end if
            end if

            if (do_j_rot) then
               s% j_rot(k) = x(i_j_rot)
               if (is_bad_num(s% j_rot(k))) then
                  s% retry_message = 'bad num for j_rot'
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver j_rot', k, s% j_rot(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num j_rot', k, s% j_rot(k)
                  ierr = -1
               end if
            end if

            if (do_lnR) then
               s% lnR(k) = x(i_lnR)
               s% r(k) = exp(s% lnR(k))
               s% dxh_lnR(k) = s% solver_dx(i_lnR,k)
               if (is_bad_num(s% lnR(k))) then
                  s% retry_message = 'bad num for lnR'
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver r lnR solver_dx', &
                        k, s% r(k), s% lnR(k), s% solver_dx(i_lnR,k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num r lnR solver_dx', &
                     k, s% r(k), s% lnR(k), s% solver_dx(i_lnR,k)
                  ierr = -1
               end if
            end if

            if (do_lnd) then
               s% lnd(k) = x(i_lnd)
               s% rho(k) = exp(s% lnd(k))
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
               if (is_bad_num(s% rho(k))) then
                  write(s% retry_message, *) 'bad num for rho', k
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'set_vars_for_solver rho', k, s% rho(k)
                     stop 'set_vars_for_solver'
                  end if
                  if (report) write(*,2) 'bad num rho', k, s% rho(k)
                  ierr = -1
               end if
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
            ierr = -1
            s% retry_message = 'bad num for mass fractions'
            if (s% stop_for_bad_nums) then
               write(*,2) 'set_vars_for_solver sum_xa', k, sum_xa
               stop 'set_vars_for_solver'
            end if
            if (report) then
               write(*,'(a60,i8,99f20.10)') 'bad num sum X', k, s% m(k)/Msun, sum_xa
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
            write(*,2) 'eps_grav', k, s% eps_grav_ad(k)% val
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
      end subroutine edit_dlnR_dt_above_k_below_just_added


      subroutine prepare_solver_matrix(s, nvar, xder, ldA, A1, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         real(dp), pointer, dimension(:,:) :: xder ! (nvar, nz)
         integer, intent(in) :: ldA ! leading dimension of A
         real(dp), pointer, dimension(:) :: A1
         integer, intent(out) :: ierr

         real(dp), pointer, dimension(:,:) :: A ! (ldA, neqns)
         integer :: i, j, nz, neqns

         include 'formats'

         ierr = 0
         nz = s% nz
         neqns = nvar*nz

         s% jacobian(1:ldA,1:neqns) => A1(1:ldA*neqns)
         A(1:ldA,1:neqns) => A1(1:ldA*neqns)

         do i=1,neqns
            do j=1,lda
               A(j,i) = 0d0
            end do
         end do
         i = nvar*nvar*nz
         if (size(A1,dim=1) < 3*i) then
            write(*,*) 'prepare_solver_matrix: size(A1,dim=1) < 3*i', size(A1,dim=1), 3*i
            ierr = -1
            return
         end if
         s% ublk(1:nvar,1:nvar,1:nz) => A1(1:i)
         s% dblk(1:nvar,1:nvar,1:nz) => A1(i+1:2*i)
         s% lblk(1:nvar,1:nvar,1:nz) => A1(2*i+1:3*i)
         
         ! delete this
         call s% other_after_enter_setmatrix(s% id,ierr)

      end subroutine prepare_solver_matrix


      end module hydro_mtx

