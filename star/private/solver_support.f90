! ***********************************************************************
!
!   Copyright (C) 2012-2019  The MESA Team
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

      module solver_support

      use star_private_def
      use utils_lib, only: is_bad
      use const_def

      use num_def

      implicit none

      logical, parameter :: dbg = .false.

      contains


      subroutine set_xscale_info(s, nvar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr

         integer :: i, j, k, nz, nvar_hydro
         real(dp), parameter :: xscale_min = 1
         real(dp) :: var_scale, lum_scale, vel_scale, omega_scale

         include 'formats'

         ierr = 0

         if (dbg) write(*, *) 'set_xscale'
         nvar_hydro = s% nvar_hydro
         nz = s% nz
         do k=1,nz
            do i=1,nvar
               if (i <= nvar_hydro) then ! structure variable
                  if (i == s% i_j_rot) then
                     s% x_scale(i,k) = 10d0*sqrt(s% cgrav(k)*s% m(k)*s% r_start(k))
                  else
                     s% x_scale(i,k) = max(xscale_min, abs(s% xh_start(i,k)))
                  end if
               else ! abundance variable
                  s% x_scale(i,k) = max(s% xa_scale, s% xa_start(i-nvar_hydro,k))
               end if
            end do
         end do

         contains

         subroutine dump_xscale
            integer :: k, j, k0, k1
            include 'formats'
            !write(*,1) 's% xa_scale', s% xa_scale
            do k=1,s% nz
               do j=1,nvar
                  write(*,2) 'xscale ' // trim(s% nameofvar(j)), k, s% x_scale(j,k)
               end do
               write(*,*)
            end do
            stop 'set_xscale'
         end subroutine dump_xscale

      end subroutine set_xscale_info


      subroutine eval_equations(s,  nvar, ierr)
         use hydro_eqns, only: eval_equ
         use mix_info, only: set_dxdt_mix
         use star_utils, only: update_time, total_times
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr

         integer :: cnt, i, j, k, nz
         integer :: id
         real(dp) :: dt, theta_dt
         include 'formats'

         ierr = 0
         nz = s% nz

         if (dbg) write(*, *) 'eval_equations'

         dt = s% dt

         if (s% solver_iter == 0) then
            if (dbg) write(*, *) 'skip set_solver_vars on call before 1st iter'
         else
            if (dbg) write(*, *) 'call set_solver_vars'
            call set_solver_vars(s, nvar, dt, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) &
                  write(*,2) 'eval_equations: set_solver_vars returned ierr', ierr
               return
            end if
         end if

         call set_dxdt_mix(s)

         if (ierr == 0) then
            do k=1,nz
               do j=1,nvar
                  s% equ(j,k) = 0d0
                  s% residual_weight(j,k) = 1d0
                  s% correction_weight(j,k) = 1d0
               end do
            end do
            if (dbg) write(*, *) 'call eval_equ'
            call eval_equ(s, nvar, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) &
                  write(*, *) 'eval_equations: eval_equ returned ierr', ierr
               return
            end if
         end if

         if (ierr /= 0) return
         
         if (.not. s% stop_for_bad_nums) return

         cnt = 0
         do i=1,nz
            do j=1, nvar
               if (is_bad_num(s% equ(j, i))) then
                  cnt = cnt + 1
                  s% retry_message = 'eval_equations: equ has a bad num'
                  if (s% report_ierr) then
                     write(*,4) 'eval_equations: equ has a bad num ' // trim(s% nameofequ(j)), &
                        j, i, nvar, s% equ(j, i)
                     write(*,2) 'cell', i
                     write(*,2) 'nz', s% nz
                  end if
                  if (s% stop_for_bad_nums) stop 'eval_equations'
               end if
            end do
         end do
         
         if (cnt > 0) then
            ierr = -1
            return
         end if


         contains


         subroutine dump_eval_equ
            integer :: k, j, k0, k1
            include 'formats'
            do k=1,s% nz
               do j=1,nvar
                  write(*,2) 'dx ' // trim(s% nameofvar(j)), k, s% solver_dx(j, k)
               end do
               write(*,*)
            end do
            stop 'dump_eval_equ'
         end subroutine dump_eval_equ


      end subroutine eval_equations


      subroutine sizequ(s, nvar, equ_norm, equ_max, k_max, j_max, ierr) ! equ = residuals
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         real(dp), intent(out) :: equ_norm, equ_max
         integer, intent(out) :: k_max, j_max, ierr

         integer :: j, k, num_terms, n, nz, nvar_hydro, nvar_chem, &
            max_loc, skip_eqn1, skip_eqn2, skip_eqn3
         real(dp) :: sumequ, absq, max_energy_resid, avg_energy_resid
         
         logical :: dbg

         include 'formats'

         ierr = 0
         
         equ_norm = 0d0
         equ_max = 0d0
         k_max = 0
         j_max = 0
         
         dbg = s% solver_check_everything

         nvar_hydro = min(nvar, s% nvar_hydro)
         nvar_chem = s% nvar_chem
         
         nz = s% nz
         n = nz
         num_terms = 0
         sumequ = 0
         skip_eqn1 = 0
         skip_eqn2 = 0
         skip_eqn3 = 0
         if (s% convergence_ignore_equL_residuals) skip_eqn1 = s% i_equL
         if (s% convergence_ignore_alpha_RTI_residuals) skip_eqn2 = s% i_dalpha_RTI_dt
         if (s% do_burn .or. s% do_mix) then
            num_terms = num_terms + nvar*nz
            if (skip_eqn1 > 0) num_terms = num_terms - nz
            if (skip_eqn2 > 0) num_terms = num_terms - nz
            if (skip_eqn3 > 0) num_terms = num_terms - nz
            do k = 1, nz
               do j = 1, nvar
                  if (j == skip_eqn1 .or. j == skip_eqn2 .or. j == skip_eqn3) cycle
                  absq = abs(s% equ(j,k)*s% residual_weight(j,k))
                  sumequ = sumequ + absq
                  if (absq > equ_max) then
                     equ_max = absq
                     j_max = j
                     k_max = k
                  end if
               end do
            end do
         else
            if (skip_eqn1 == 0 .and. skip_eqn2 == 0) then
               num_terms = num_terms + nvar_hydro*nz
            else if (skip_eqn1 > 0 .and. skip_eqn2 > 0) then
               num_terms = num_terms + (nvar_hydro-2)*nz
            else
               num_terms = num_terms + (nvar_hydro-1)*nz
            end if
            do k = 1, nz
               do j = 1, nvar_hydro
                  if (j == skip_eqn1 .or. j == skip_eqn2) cycle
                  absq = abs(s% equ(j,k)*s% residual_weight(j,k))
                  sumequ = sumequ + absq
                  if (is_bad(sumequ)) then
                     if (dbg) then
                        write(*,3) trim(s% nameofequ(j)) // ' sumequ', j, k, sumequ
                        stop 'sizeq 1'
                     end if
                     ierr = -1
                     if (s% report_ierr) &
                        write(*,3) 'bad equ(j,k)*s% residual_weight(j,k) ' // trim(s% nameofequ(j)), &
                           j, k, s% equ(j,k)*s% residual_weight(j,k)
                     if (s% stop_for_bad_nums) stop 'sizeq 2'
                     return
                  end if
                  if (absq > equ_max) then
                     equ_max = absq
                     j_max = j
                     k_max = k
                  end if
               end do
            end do
         end if
         if (s% do_burn .or. s% do_mix) then
            num_terms = num_terms + nvar_chem*nz
            do k = 1, nz
               do j = nvar_hydro+1, nvar
                  absq = abs(s% equ(j,k)*s% residual_weight(j,k))
                  sumequ = sumequ + absq
                  if (absq > equ_max) then
                     equ_max = absq
                     j_max = j
                     k_max = k
                  end if
               end do
            end do
         end if
         if (s% conv_vel_flag) then
            do k = 1, nz
               j = s% i_dln_cvpv0_dt
               absq = abs(s% equ(j,k)*s% residual_weight(j,k))
            end do
         end if

         equ_norm = sumequ/num_terms
         if (dbg) write(*,4) trim(s% nameofequ(j_max)) // ' sizequ equ_max norm', &
            k_max, s% solver_iter, s% model_number, equ_max, equ_norm
         
         if (dbg) call dump_equ
         
         return
         call dump_equ
         stop 'sizequ'
         
         contains

         subroutine dump_equ
            integer :: k, j, k0, k1
            include 'formats'
            do k=1,s% nz
               do j=1,nvar
                  write(*,3) 'equ ' // trim(s% nameofequ(j)), &
                     k, s% solver_iter, s% equ(j, k)
               end do
               write(*,*)
               !if (k == 6) exit
            end do
         end subroutine dump_equ

      end subroutine sizequ


      subroutine sizeB(s, nvar, B, max_correction, correction_norm, max_zone, max_var, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         real(dp), pointer, dimension(:,:) :: B ! (nvar, nz)
         real(dp), intent(out) :: correction_norm ! a measure of the average correction
         real(dp), intent(out) :: max_correction ! magnitude of the max correction
         integer, intent(out) :: max_zone, max_var, ierr

         integer :: k, i, nz, num_terms, j, n, nvar_hydro, jmax, num_xa_terms, &
            skip1, skip2, skip3, skip4, skip5
         real(dp) :: abs_corr, sum_corr, sum_xa_corr, x_limit, &
            max_abs_correction, max_abs_correction_cv, max_abs_corr_for_k, max_abs_xa_corr_for_k
         logical :: found_NaN, found_bad_num, report
         real(dp), parameter :: frac = 0.1d0
         logical, parameter :: dbg = .false.
         logical, parameter :: check_for_bad_nums = .true.
         logical, parameter :: save_max_abs_corr_for_k = .true.

         include 'formats'

         if (dbg) write(*, *) 'enter sizeB'

         ierr = 0
         nz = s% nz
         n = nz
         nvar_hydro = min(nvar, s% nvar_hydro)

         if (s% include_L_in_correction_limits) then
            skip1 = 0
            do k=1,nz
               s% correction_weight(s% i_lum,k) = 1d0/(frac*s% L_start(1) + abs(s% L(k)))
            end do
         else
            skip1 = s% i_lum
         end if

         if (s% u_flag .and. s% include_u_in_correction_limits) then
            skip2 = 0
            do k=1,nz
               s% correction_weight(s% i_u,k) = 1d0/(frac*s% csound_start(k) + abs(s% u(k)))
            end do
         else if (s% v_flag .and. s% include_v_in_correction_limits) then
            skip2 = 0
            do k=1,nz
               s% correction_weight(s% i_v,k) = 1d0/(frac*s% csound_start(k) + abs(s% v(k)))
            end do
         else if (s% u_flag) then
            skip2 = s% i_u
         else
            skip2 = s% i_v
         end if

         if (s% using_RSP2 .and. s% include_w_in_correction_limits) then
            skip3 = 0
            do k=1,nz
               if (abs(s% w(k)) < 1d0) then
                  s% correction_weight(s% i_w,k) = 0d0
               else
                  s% correction_weight(s% i_w,k) = 1d0/(1d3 + abs(s% w(k))) ! don't sweat the small stuff
               end if
            end do
         else
            skip3 = s% i_w
         end if
         skip4 = s% i_Hp
                 
         skip5 = 0

         max_zone = 0
         max_var = 0
         num_terms = 0
         num_xa_terms = 0
         sum_corr = 0
         sum_xa_corr = 0
         max_correction = 0
         max_abs_correction = 0
         x_limit = s% correction_xa_limit
         found_NaN = .false.
         found_bad_num = .false.
         report = s% report_ierr
         cell_loop: do k = 1, nz
            max_abs_corr_for_k = 0
            max_abs_xa_corr_for_k = 0

            if (s% do_burn .or. s% do_mix) then
               jmax = nvar
            else
               jmax = nvar_hydro
            end if
            var_loop: do j = 1, jmax
               if (j == skip1 .or. &
                   j == skip2 .or. &
                   j == skip3 .or. &
                   j == skip4 .or. &
                   j == skip5 .or. &
                   j == s% i_alpha_RTI) cycle
               if (check_for_bad_nums) then
                  if (is_bad_num(B(j,k)*s% correction_weight(j,k))) then
                     found_bad_num = .true.
                     if (report) write(*,2) 'sizeB: bad num for correction ' // &
                        s% nameofvar(j), k, B(j,k)*s% correction_weight(j,k)
                     if (s% stop_for_bad_nums) then
                        found_NaN = .true.
                        write(*,3) s% nameofvar(j) // ' B(j,k)*s% correction_weight(j,k)', &
                           j, k, B(j,k)*s% correction_weight(j,k)
                        stop 'sizeB'
                     end if
                     
                     max_zone = k
                     max_var = j
                     exit cell_loop
                     
                     cycle
                  end if
               end if
               if (j > nvar_hydro) then
                  if (s% xa_start(j-nvar_hydro,k) < x_limit) cycle
               end if

               abs_corr = abs(B(j,k)*s% correction_weight(j,k))
               if (is_bad_num(abs_corr)) then
                  found_bad_num = .true.
                  if (report) write(*,3) 'B(j,k)*s% correction_weight(j,k)', &
                     j, k, B(j,k)*s% correction_weight(j,k)
                  if (s% stop_for_bad_nums) found_NaN = .true.
               end if
               if (abs_corr > max_abs_corr_for_k &
                  .and. .not. (j > nvar_hydro .and. s% ignore_species_in_max_correction)) &
                     max_abs_corr_for_k = abs_corr
               if (abs_corr > max_abs_correction &
                  .and. .not. (j > nvar_hydro .and. s% ignore_species_in_max_correction)) then
                  max_correction = B(j,k)*s% correction_weight(j,k)
                  max_abs_correction = abs_corr
                  max_zone = k
                  max_var = j
               end if
               if (j > nvar_hydro) then
                  num_xa_terms = num_xa_terms + 1
                  sum_xa_corr = sum_xa_corr + abs_corr
                  if (abs_corr > max_abs_xa_corr_for_k) &
                     max_abs_xa_corr_for_k = abs_corr
               else
                  num_terms = num_terms + 1
                  sum_corr = sum_corr + abs_corr
               end if
            end do var_loop
            if (num_xa_terms > 0) then
               num_terms = num_terms + 1
               sum_corr = sum_corr + sum_xa_corr/num_xa_terms
            end if
            if (s% do_burn .or. s% do_mix) then
               species_loop: do j = nvar_hydro+1, nvar
                  i = j - s% nvar_hydro
                  if (check_for_bad_nums) then
                     if (is_bad_num(B(j,k)*s% correction_weight(j,k))) then
                        found_bad_num = .true.
                        if (report) write(*,3) 'chem B(j,k)*s% correction_weight(j,k)', j, k, &
                           B(j,k)*s% correction_weight(j,k)
                        if (s% stop_for_bad_nums) then
                           found_NaN = .true.
                           write(*,3) 'chem B(j,k)*s% correction_weight(j,k)', &
                              j, k, B(j,k)*s% correction_weight(j,k)
                           stop 'sizeB'
                        max_zone = k
                        max_var = j
                        exit cell_loop
                        end if
                     end if
                  end if
                  ! recall that correction dx = B*xscale, so B is a relative correction
                  if (s% xa_start(i,k) >= x_limit) then
                     abs_corr = abs(B(j,k)*s% correction_weight(j,k))
                     if (abs_corr > max_abs_corr_for_k) max_abs_corr_for_k = abs_corr
                     if (abs_corr > max_abs_correction) then
                        max_abs_correction = abs_corr
                        max_correction = B(j,k)*s% correction_weight(j,k)
                        max_zone = k
                        max_var = j
                     end if
                     sum_corr = sum_corr + abs_corr
                     num_terms = num_terms + 1
                  end if
               end do species_loop
            end if
            s% max_abs_xa_corr(k) = max_abs_xa_corr_for_k
         end do cell_loop

         if (found_bad_num) then
            ierr = -1
            if (found_NaN .and. s% stop_for_bad_nums) then
               write(*,*) 'found bad num'
               stop 'sizeB'
            end if
            if (.not. dbg) return
         end if

         if (is_bad_num(sum_corr)) then
            ierr = -1
            if (s% stop_for_bad_nums) then
               if (report) write(*,*) 'sum_corr', sum_corr
               stop 'sizeB'
            end if
            if (.not. dbg) return
            write(*,*) 'sum_corr', sum_corr
            stop 'sizeB'
         end if

         correction_norm = sum_corr/num_terms  !sqrt(sum_corr/num_terms)
         if (dbg) then
            write(*,2) 'sizeB: iter, correction_norm, max_correction', &
               s% solver_iter, correction_norm, max_correction
            if (max_correction > 1d50 .or. is_bad_num(correction_norm)) then
               call show_stuff
               stop 'sizeB'
            end if
         end if

         if (s% solver_show_correction_info) call show_stuff

         abs_corr = max_abs_correction

         s% abs_max_corr2 = s% abs_max_corr1; s% abs_max_corr1 = abs_corr
         s% max_var2 = s% max_var1; s% max_var1 = max_var
         s% max_zone2 = s% max_zone1; s% max_zone1 = max_zone

         if (ierr /= 0) stop 'ierr in sizeB'

         if (is_bad_num(max_correction)) then
            ierr = -1
            if (s% stop_for_bad_nums) then
               if (report) write(*,*) 'max_correction', max_correction
               stop 'sizeB'
            end if
            if (.not. dbg) return
            write(*,*) 'max_correction', max_correction
            stop 'sizeB'
         end if

         if (s% solver_iter < 3) return
         ! check for flailing
         if ( &
             abs_corr > s% tol_max_correction .and. &
             abs_corr > s% abs_max_corr1 .and. s% abs_max_corr1 > s% abs_max_corr2 .and. &
             max_zone == s% max_zone1 .and. s% max_zone1 == s% max_zone2 .and. &
             max_var == s% max_var1 .and. s% max_var1 == s% max_var2) then
            if (s% solver_show_correction_info) then
               write(*,*) 'give up because diverging'
            end if
            max_correction = 1d99
         end if


         contains


         subroutine show_stuff
            integer :: j, k
            real(dp) :: dx, prev, new
            include 'formats'
            if (s% solver_iter == 1) then
               write(*,*)
               write(*,'(4a7,12a16,99a13)') &
                  'model', 'iter', 'var', 'zone', &
                  'corr norm', 'max corr', 'xscale', &
                  'dx', 'new-prev', 'new', 'prev', &
                  'mass loc', 'log dt/yr', 'lgE', 'lgT', 'lgRho'
            end if
            k = max_zone
            j = max_var
            if (j > nvar_hydro) then
               prev = s% xa_start(j - nvar_hydro,k)
            else
               prev = s% xh_start(j,k)
            end if
            dx = B(j,k)*s% correction_weight(j,k)*s% x_scale(j,k)
            new = prev + dx
            write(*,'(2i7,a7,i7,12e16.8,99f13.8)') &
               s% model_number, s% solver_iter, trim(s% nameofvar(max_var)), k, &
               correction_norm, B(j,k)*s% correction_weight(j,k), s% x_scale(j,k), &
               dx, new - prev, new, prev, &
               s% m(k)/Msun, log10(s% dt/secyer), &
               s% lnE(k)/ln10, s% lnT(k)/ln10, &
               s% lnd(k)/ln10
         end subroutine show_stuff


         subroutine dump_B
            integer :: k, j, k0, k1
            include 'formats'
            do k=1,s% nz
               do j=1,nvar
                  write(*,2) 'B ' // trim(s% nameofequ(j)), k, B(j, k)
               end do
               write(*,*)
            end do
            stop 'dump_equ'
         end subroutine dump_B


      end subroutine sizeB


      ! the proposed change to dx is B*xscale*correction_factor
      ! edit correction_factor and/or B as necessary so that the new dx will be valid.
      ! set ierr nonzero if things are beyond repair.
      subroutine Bdomain(s, nvar, B, correction_factor, ierr)
         use const_def, only: dp
         use chem_def, only: chem_isos
         use star_utils, only: current_min_xa_hard_limit, rand
         use rsp_def, only: EFL0
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         real(dp), pointer, dimension(:,:) :: B ! (nvar, nz)
         real(dp), intent(inout) :: correction_factor
         integer, intent(out) :: ierr
         integer :: id, i, j, k, nz, species, bad_j, bad_k, &
            i_alpha_RTI, i_ln_cvpv0, i_w_div_wc
         real(dp) :: alpha, min_alpha, new_xa, old_xa, dxa, eps, min_xa_hard_limit, &
            old_E, dE, new_E, old_lnd, dlnd, new_lnd, dw, new_w, &
            dw_div_wc, old_w_div_wc, new_w_div_wc, dconv_vel, old_conv_vel, new_conv_vel, &
            dalpha_RTI, new_alpha_RTI, old_alpha_RTI, log_conv_vel_v0, &
            dlum_surf, old_lum_surf, new_lum_surf
         include 'formats'
         ierr = 0
         min_alpha = 1d0
         nz = s% nz

         if (s% using_RSP2) & ! clip change in w to maintain non-negativity.
            call clip_so_non_negative(s% i_w, 0d0)

         if (s% RTI_flag) & ! clip change in alpha_RTI to maintain non-negativity.
            call clip_so_non_negative(s% i_alpha_RTI, 0d0)

         if (s% conv_vel_flag) then ! clip change in conv_vel to maintain non-negativity.
            log_conv_vel_v0 = log(s% conv_vel_v0)
            i_ln_cvpv0 = s% i_ln_cvpv0
            !note that dconv_vel and others refers to changes in ln(conv_vel+v0)
            do k = 1, s% nz
               dconv_vel = B(i_ln_cvpv0,k)*s% x_scale(i_ln_cvpv0,k)*correction_factor
               old_conv_vel = s% xh_start(i_ln_cvpv0,k) + s% solver_dx(i_ln_cvpv0,k)
               new_conv_vel = old_conv_vel + dconv_vel
               if (new_conv_vel >= log_conv_vel_v0) cycle
               dconv_vel = -old_conv_vel + log_conv_vel_v0
               B(i_ln_cvpv0,k) = dconv_vel/(s% x_scale(i_ln_cvpv0,k)*correction_factor)
            end do
         end if

         if (s% i_lum>=0 .and. s% scale_max_correction_for_negative_surf_lum) then
            !ensure surface luminosity does not become negative
            dlum_surf = B(s% i_lum,1)*s% x_scale(s% i_lum,1)
            old_lum_surf = s% xh_start(s% i_lum,1) + s% solver_dx(s% i_lum,1)
            new_lum_surf = old_lum_surf + dlum_surf
            if (new_lum_surf < 0d0 .and. old_lum_surf > 0d0) then
               correction_factor = min(correction_factor, &
                  -s% max_frac_for_negative_surf_lum*old_lum_surf/dlum_surf)
            end if
         end if

         !if (s% w_div_wc_flag) then
         !   do k=1,nz
         !      old_xa = s% xh_start(s% i_w_div_wc,k) + dx(s% i_w_div_wc,k)
         !      dxa = B(s% i_w_div_wc,k)*xscale(s% i_w_div_wc,k)*correction_factor
         !      if (dxa > 0.9d0*(1-old_xa)) then
         !         if(k==91)write(*,*) "problems", old_xa, dxa, old_xa+dxa
         !         correction_factor = min(correction_factor,(0.9d0*(1-old_xa))/B(s% i_w_div_wc,k)*xscale(s%i_w_div_wc,k))
         !         dxa = B(s% i_w_div_wc,k)*xscale(s% i_w_div_wc,k)*correction_factor
         !         if(k==91)write(*,*) "need to reduce correction", k, old_xa+dxa
         !      end if
         !   end do
         !end if

         if ((.not. s% do_solver_damping_for_neg_xa) .or. &
             (.not. (s% do_mix .or. s% do_burn))) then
            if (min_alpha < 1d0) then
               min_alpha = max(min_alpha, s% corr_coeff_limit)
               correction_factor = min_alpha*correction_factor
            end if
            return
         end if

         bad_j = 0
         bad_k = 0
         if (nvar > s% nvar_hydro) then
            species = s% species
            min_xa_hard_limit = current_min_xa_hard_limit(s)
            eps = -0.5d0*min_xa_hard_limit ! allow xa to be this much below 0d0
            do k=1,nz
               do j=1,species
                  i = j + s% nvar_hydro
                  old_xa = s% xa_start(j,k) + s% solver_dx(i,k)
                  if (old_xa <= 1d-90) cycle
                  dxa = B(i,k)*s% x_scale(i,k)*correction_factor
                  new_xa = old_xa + dxa
                  if (new_xa >= 0d0) cycle
                  alpha = -(old_xa + eps)/dxa
                  ! so dxa*alpha = -old_xa - eps,
                  ! and therefore old_xa + alpha*dxa = -eps = 0.5*min_xa_hard_limit
                  if (alpha < min_alpha) then
                     min_alpha = alpha
                     bad_j = j
                     bad_k = k
                  end if
               end do
            end do
         end if

         min_alpha = max(min_alpha, s% corr_coeff_limit)
         correction_factor = min_alpha*correction_factor
         if (s% trace_solver_damping .and. min_alpha < 1d0 .and. bad_j > 0) then
            write(*,4) 'solver damping to avoid negative mass fractions: ' // &
               trim(chem_isos% name(s% chem_id(bad_j))), bad_k, &
               s% model_number, s% solver_iter, min_alpha
         end if
         
         contains
         
         subroutine clip_so_non_negative(i,minval)
            integer, intent(in) :: i
            real(dp), intent(in) :: minval
            real(dp) :: dval, old_val, new_val
            do k = 1, s% nz
               dval = B(i,k)*s% x_scale(i,k)*correction_factor
               old_val = s% xh_start(i,k) + s% solver_dx(i,k)
               new_val = old_val + dval
               if (dval >= 0) cycle
               if (new_val >= 0d0) cycle
               dval = minval - old_val
               B(i,k) = dval/(s% x_scale(i,k)*correction_factor)
            end do
         end subroutine clip_so_non_negative

      end subroutine Bdomain


      subroutine inspectB(s, nvar, B, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         real(dp), pointer, dimension(:,:) :: B ! (nvar, nz)
         integer, intent(out) :: ierr

         integer :: id
         integer, parameter :: inspectB_iter_stop = -1
         include 'formats'

         if (dbg) write(*, *) 'inspectB', s% solver_iter
         ierr = 0
         if (s% solver_iter == inspectB_iter_stop) then
            call dumpB
            stop 'debug: inspectB'
         end if

         contains

         subroutine dumpB
            integer :: k, j, k0, k1
            include 'formats'
            do k=1,s% nz
               do j=1,nvar
                  write(*,2) 'B ' // trim(s% nameofvar(j)), k, B(j, k)
                  write(*,2) 'xscale ' // trim(s% nameofvar(j)), k, s% x_scale(j, k)
                  write(*,2) 'dx ' // trim(s% nameofvar(j)), k, s% solver_dx(j, k)
               end do
               write(*,*)
            end do
            stop 'dumpB'
         end subroutine dumpB

      end subroutine inspectB


      ! about to declare victory... but may want to do another iteration
      ! 1 means force another iteration
      ! 0 means don't need to force another
      ! -1 means failure. solver returns with non-convergence.
      integer function force_another_iteration(s, iter, itermin)
         type (star_info), pointer :: s
         integer, intent(in) :: iter ! have finished this many iterations and have converged
         integer, intent(in) :: itermin ! this is the requested minimum.  iter may be < itermin.

         integer :: k, res

         include 'formats'

         if (iter < itermin) then
            force_another_iteration = 1
            return
         end if
         force_another_iteration = 0

         if (s% k_below_just_added > 1 .and. &
               s% num_surf_revisions < s% max_num_surf_revisions .and. &
               abs(s% lnS(1) - s% surf_lnS) > &
                  s% max_abs_rel_change_surf_lnS*max(s% lnS(1),s% surf_lnS)) then
            s% surf_lnT = s% lnT(1)
            s% surf_lnR = s% lnR(1)
            s% surf_lnd = s% lnd(1)
            if (s% i_v /= 0) s% surf_v = s% v(1)
            if (s% i_u /= 0) s% surf_v = s% u_face_ad(1)%val
            s% surf_lnS = s% lnS(1)
            s% num_surf_revisions = s% num_surf_revisions + 1
            force_another_iteration = 1
            s% used_extra_iter_in_solver_for_accretion = .true.
            return
         end if

      end function force_another_iteration


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
               if (s% w(k) < 0d0) then
                  !write(*,4) 'set_vars_for_solver: fix w < 0', k, &
                  !   s% solver_iter, s% model_number, s% w(k)
                  s% w(k) = s% RSP2_w_fix_if_neg
               end if
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


      end module solver_support

