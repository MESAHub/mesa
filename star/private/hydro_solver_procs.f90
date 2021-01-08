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

      module hydro_solver_procs
      use hydro_mtx

      use star_private_def
      use utils_lib, only: is_bad
      use const_def

      use num_def

      implicit none


      contains


      subroutine set_xscale_info(s, nvar, nz, xscale, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nz
         real(dp), pointer :: xscale(:,:) ! (nvar, nz)
         integer, intent(out) :: ierr

         integer :: i, j, k, nvar_hydro
         real(dp), parameter :: xscale_min = 1
         real(dp) :: var_scale, lum_scale, vel_scale, omega_scale

         include 'formats'

         ierr = 0

         if (dbg) write(*, *) 'set_xscale'
         nvar_hydro = s% nvar_hydro

         do k=1,nz
            do i=1,nvar
               if (i <= nvar_hydro) then ! structure variable
                  if (i /= s% i_j_rot) then
                     xscale(i,k) = max(xscale_min, abs(s% xh_start(i,k)))
                  else
                     xscale(i,k) = 10d0*sqrt(s% cgrav(k)*s% m(k)*s% r_start(k))
                  end if
               else ! abundance variable
                  xscale(i,k) = max(s% xa_scale, s% xa_start(i-nvar_hydro,k))
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
                  write(*,2) 'xscale ' // trim(s% nameofvar(j)), k, xscale(j,k)
               end do
               write(*,*)
            end do
            stop 'set_xscale'
         end subroutine dump_xscale

      end subroutine set_xscale_info


      subroutine eval_equations( &
            iter, nvar, nz, dx, xscale, equ_in, lrpar, rpar, lipar, ipar, ierr)
         use hydro_eqns, only: eval_equ
         use mix_info, only: set_dxdt_mix
         use star_utils, only: update_time, total_times
         integer, intent(in) :: iter, nvar, nz
         real(dp), pointer, dimension(:,:) :: dx, xscale, equ_in ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr

         integer :: cnt, i, j, k
         type (star_info), pointer :: s
         integer :: id
         real(dp) :: dt, theta_dt
         real(dp), pointer :: equ(:,:)

         logical, parameter :: skip_partials = .true.

         include 'formats'

         ierr = 0

         id = ipar(ipar_id)
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         equ(1:nvar,1:nz) => s% equ1(1:nvar*nz)

         if (dbg) write(*, *) 'eval_equations'

         dt = rpar(rpar_dt)

         if (ipar(ipar_first_call) /= 0) then
            ipar(ipar_first_call) = 0
            if (dbg) write(*, *) 'skip set_solver_vars on call before 1st iter'
         else
            if (dbg) write(*, *) 'call set_solver_vars'
            call set_solver_vars(s, iter, nvar, dx, xscale, dt, ierr)
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
                  equ(j,k) = 0d0
                  s% residual_weight(j,k) = 1d0
                  s% correction_weight(j,k) = 1d0
               end do
            end do
            if (dbg) write(*, *) 'call eval_equ'
            call eval_equ(s, nvar, skip_partials, xscale, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) &
                  write(*, *) 'eval_equations: eval_equ returned ierr', ierr
               return
            end if
         end if

         if (ierr /= 0) return

         cnt = 0
         do i=1,nz
            do j=1, nvar
               if (is_bad_num(equ(j, i))) then
                  cnt = cnt + 1
                  s% retry_message = 'eval_equations: equ has a bad num'
                  if (s% report_ierr) then
                     write(*,4) 'eval_equations: equ has a bad num ' // trim(s% nameofequ(j)), &
                        j, i, nvar, equ(j, i)
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
                  write(*,2) 'dx ' // trim(s% nameofvar(j)), k, dx(j, k)
               end do
               write(*,*)
            end do
            stop 'dump_eval_equ'
         end subroutine dump_eval_equ


      end subroutine eval_equations


      subroutine sizequ(s, &
            iter, nvar, nz, equ, &
            equ_norm, equ_max, k_max, j_max, &
            lrpar, rpar, lipar, ipar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: iter, nvar, nz
         real(dp), pointer :: equ(:,:) ! (nvar, nz)
         real(dp), intent(out) :: equ_norm, equ_max
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: k_max, j_max, ierr

         integer :: j, k, num_terms, n, i_chem1, nvar_hydro, nvar_chem, &
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
         n = nz
         num_terms = 0
         sumequ = 0
         skip_eqn1 = 0
         skip_eqn2 = 0
         skip_eqn3 = 0
         if (s% convergence_ignore_equL_residuals) skip_eqn1 = s% i_equL
         if (s% convergence_ignore_alpha_RTI_residuals) skip_eqn2 = s% i_dalpha_RTI_dt
         if (s% do_struct_hydro .or. s% do_struct_thermo) then
            if (s% do_burn .or. s% do_mix) then
               num_terms = num_terms + nvar*nz
               if (skip_eqn1 > 0) num_terms = num_terms - nz
               if (skip_eqn2 > 0) num_terms = num_terms - nz
               if (skip_eqn3 > 0) num_terms = num_terms - nz
               do k = 1, nz
                  do j = 1, nvar
                     if (j == skip_eqn1 .or. j == skip_eqn2 .or. j == skip_eqn3) cycle
                     absq = abs(equ(j,k)*s% residual_weight(j,k))
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
                     absq = abs(equ(j,k)*s% residual_weight(j,k))
                     !write(*,3) 'equ(j,k)*s% residual_weight(j,k)', j, k, equ(j,k)*s% residual_weight(j,k)
                     sumequ = sumequ + absq
                     if (is_bad(sumequ)) then
                        if (dbg) then
                           write(*,3) trim(s% nameofequ(j)) // ' sumequ', j, k, sumequ
                           stop 'sizeq 1'
                        end if
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,3) 'bad equ(j,k)*s% residual_weight(j,k) ' // trim(s% nameofequ(j)), &
                              j, k, equ(j,k)*s% residual_weight(j,k)
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
         end if
         if (s% do_burn .or. s% do_mix) then
            i_chem1 = s% i_chem1
            num_terms = num_terms + nvar_chem*nz
            do k = 1, nz
               do j = i_chem1, nvar
                  absq = abs(equ(j,k)*s% residual_weight(j,k))
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
               absq = abs(equ(j,k)*s% residual_weight(j,k))
            end do
         end if

         equ_norm = sumequ/num_terms
         if (dbg) write(*,4) trim(s% nameofequ(j_max)) // ' sizequ equ_max norm', &
            k_max, iter, s% model_number, equ_max, equ_norm
         
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
                     k, iter, equ(j, k)
               end do
               write(*,*)
               !if (k == 6) exit
            end do
         end subroutine dump_equ

      end subroutine sizequ


      subroutine sizeB(s, &
            iter, nvar, nz, B, xscale, &
            max_correction, correction_norm, max_zone, max_var, &
            lrpar, rpar, lipar, ipar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: iter, nvar, nz
         real(dp), pointer, dimension(:,:) :: B, xscale ! (nvar, nz)
         real(dp), intent(out) :: correction_norm ! a measure of the average correction
         real(dp), intent(out) :: max_correction ! magnitude of the max correction
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: max_zone, max_var, ierr

         integer :: k, i, num_terms, j, n, nvar_hydro, &
            skip1, skip2, skip3, jmax, num_xa_terms, i_alpha_RTI, i_ln_cvpv0
         real(dp) :: abs_corr, sum_corr, sum_xa_corr, x_limit, &
            max_abs_correction, max_abs_correction_cv, max_abs_corr_for_k, max_abs_xa_corr_for_k
         logical :: found_NaN, found_bad_num, report
         logical, parameter :: dbg = .false.
         logical, parameter :: check_for_bad_nums = .true.
         logical, parameter :: save_max_abs_corr_for_k = .true.

         include 'formats'

         if (dbg) write(*, *) 'enter sizeB'

         ierr = 0
         n = nz
         nvar_hydro = min(nvar, s% nvar_hydro)

         if (s% include_L_in_error_est) then
            skip1 = 0
         else
            skip1 = s% i_lum
         end if

         if (s% include_v_in_error_est) then
            skip2 = 0
         else
            skip2 = s% i_v
         end if

         if (s% include_u_in_error_est) then
            skip3 = 0
         else
            skip3 = s% i_u
         end if
         
         i_alpha_RTI = s% i_alpha_RTI
         i_ln_cvpv0 = s% i_ln_cvpv0

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

            if (s% do_struct_hydro .or. s% do_struct_thermo) then
               if (s% do_burn .or. s% do_mix) then
                  jmax = nvar
               else
                  jmax = nvar_hydro
               end if
               var_loop: do j = 1, jmax
                  if (j == skip1 .or. &
                      j == skip2 .or. &
                      j == skip3 .or. &
                      j == i_alpha_RTI) cycle
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
            else if (s% do_burn .or. s% do_mix) then
               species_loop: do j = s% i_chem1, nvar
                  i = j - s% nvar_hydro
                  if (check_for_bad_nums) then
                     if (is_bad_num(B(j,k)*s% correction_weight(j,k))) then
                        found_bad_num = .true.
                        if (report) write(*,3) 'chem B(j,k)*s% correction_weight(j,k)', j, k, B(j,k)*s% correction_weight(j,k)
                        if (s% stop_for_bad_nums) then
                           found_NaN = .true.
                           write(*,3) 'chem B(j,k)*s% correction_weight(j,k)', j, k, B(j,k)*s% correction_weight(j,k)
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
               iter, correction_norm, max_correction
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

         if (iter < 3) return
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
            if (iter == 1) then
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
            dx = B(j,k)*s% correction_weight(j,k)*xscale(j,k)
            new = prev + dx
            write(*,'(2i7,a7,i7,12e16.8,99f13.8)') &
               s% model_number, iter, trim(s% nameofvar(max_var)), k, &
               correction_norm, B(j,k)*s% correction_weight(j,k), xscale(j,k), &
               dx, new - prev, new, prev, &
               s% m(k)/Msun, log10(rpar(rpar_dt)/secyer), &
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
      subroutine Bdomain( &
            iter, nvar, nz, B, dx, xscale, correction_factor, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only: dp
         use chem_def, only: chem_isos
         use star_utils, only: current_min_xa_hard_limit, rand
         use rsp_def, only: EFL0
         integer, intent(in) :: iter, nvar, nz
         real(dp), pointer, dimension(:,:) :: dx, xscale, B ! (nvar, nz)
         real(dp), intent(inout) :: correction_factor
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr
         integer :: id, i, j, k, species, bad_j, bad_k, &
            i_alpha_RTI, i_ln_cvpv0, i_w_div_wc, i_eturb
         real(dp) :: alpha, min_alpha, new_xa, old_xa, dxa, eps, min_xa_hard_limit, &
            old_E, dE, new_E, old_lnd, dlnd, new_lnd, deturb, old_eturb, new_eturb, &
            dw_div_wc, old_w_div_wc, new_w_div_wc, &
            dconv_vel, old_conv_vel, new_conv_vel, dEt, old_Et, new_Et, &
            dalpha_RTI, new_alpha_RTI, old_alpha_RTI, log_conv_vel_v0, &
            dlum_surf, old_lum_surf, new_lum_surf
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         id = ipar(ipar_id)
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         min_alpha = 1d0

         if (s% Eturb_flag) then ! clip change in eturb to maintain non-negativity.
            i_eturb = s% i_eturb
            do k = 1, s% nz
               deturb = B(i_eturb,k)*xscale(i_eturb,k)*correction_factor
               old_eturb = s% xh_start(i_eturb,k) + dx(i_eturb,k)
               new_eturb = old_eturb + deturb
               if (deturb >= 0) cycle
               if (new_eturb >= 0d0) cycle
               deturb = min_eturb*1d-6 - old_eturb
               B(i_eturb,k) = deturb/(xscale(i_eturb,k)*correction_factor)
            end do
         end if

         if (s% RTI_flag) then ! clip change in alpha_RTI to maintain non-negativity.
            i_alpha_RTI = s% i_alpha_RTI
            do k = 1, s% nz
               dalpha_RTI = B(i_alpha_RTI,k)*xscale(i_alpha_RTI,k)*correction_factor
               if (dalpha_RTI >= 0) cycle
               old_alpha_RTI = s% xh_start(i_alpha_RTI,k) + dx(i_alpha_RTI,k)
               new_alpha_RTI = old_alpha_RTI + dalpha_RTI
               if (new_alpha_RTI >= 0d0) cycle
               dalpha_RTI = -old_alpha_RTI
               B(i_alpha_RTI,k) = dalpha_RTI/(xscale(i_alpha_RTI,k)*correction_factor)
            end do
         end if

         if (s% conv_vel_flag) then ! clip change in conv_vel to maintain non-negativity.
            log_conv_vel_v0 = log(s% conv_vel_v0)
            i_ln_cvpv0 = s% i_ln_cvpv0
            !note that dconv_vel and others refers to changes in ln(conv_vel+v0)
            do k = 1, s% nz
               dconv_vel = B(i_ln_cvpv0,k)*xscale(i_ln_cvpv0,k)*correction_factor
               old_conv_vel = s% xh_start(i_ln_cvpv0,k) + dx(i_ln_cvpv0,k)
               new_conv_vel = old_conv_vel + dconv_vel
               if (new_conv_vel >= log_conv_vel_v0) cycle
               dconv_vel = -old_conv_vel + log_conv_vel_v0
               B(i_ln_cvpv0,k) = dconv_vel/(xscale(i_ln_cvpv0,k)*correction_factor)
            end do
         end if

         !if (s% w_div_wc_flag) then ! clip change in w_div_wc to keep it between 0 and 1.
         !   i_w_div_wc = s% i_w_div_wc
         !   do k = 1, s% nz
         !      dw_div_wc = B(i_w_div_wc,k)*xscale(i_w_div_wc,k)*correction_factor
         !      old_w_div_wc = s% xh_start(i_w_div_wc,k) + dx(i_w_div_wc,k)
         !      new_w_div_wc = old_w_div_wc + dw_div_wc
         !      !if (k== 19) write(*,*) "check new_w_div_wc", new_w_div_wc, dw_div_wc, s% j_rot(k), &
         !      !   correction_factor
         !      if (new_w_div_wc > 0.99) then
         !         dw_div_wc = 0.9d0*(0.99d0-old_w_div_wc)
         !         B(i_w_div_wc,k) = dw_div_wc/(xscale(i_w_div_wc,k)*correction_factor)
         !         !if (k== 19) write(*,*) "adjust new_w_div_wc", old_w_div_wc+dw_div_wc, dw_div_wc 
         !      else if (new_w_div_wc < -0.99d0) then
         !         dw_div_wc = 0.9d0*(-0.99d0-old_w_div_wc)
         !         B(i_w_div_wc,k) = dw_div_wc/(xscale(i_w_div_wc,k)*correction_factor)
         !         !if (k== 19) write(*,*) "adjust new_w_div_wc", old_w_div_wc+dw_div_wc, dw_div_wc 
         !      end if
         !   end do
         !end if
         !if (s% w_div_wc_flag) then ! clip change in w_div_wc to keep it between 0 and 1.
         !   i_w_div_wc = s% i_w_div_wc
         !   do k = 1, s% nz
         !      dw_div_wc = B(i_w_div_wc,k)*xscale(i_w_div_wc,k)*correction_factor
         !      old_w_div_wc = s% xh_start(i_w_div_wc,k) + dx(i_w_div_wc,k)
         !      new_w_div_wc = old_w_div_wc + dw_div_wc
         !      !if (k== 19) write(*,*) "check new_w_div_wc", new_w_div_wc, dw_div_wc, s% j_rot(k), &
         !      !   correction_factor
         !      if (dw_div_wc > 0.05d0) then
         !         B(i_w_div_wc,k) = 0.05d0/(xscale(i_w_div_wc,k)*correction_factor)
         !         !if (k== 19) write(*,*) "adjust new_w_div_wc", old_w_div_wc+dw_div_wc, dw_div_wc 
         !      else if (dw_div_wc < -0.05d0) then
         !         B(i_w_div_wc,k) = -0.05d0/(xscale(i_w_div_wc,k)*correction_factor)
         !         !if (k== 19) write(*,*) "adjust new_w_div_wc", old_w_div_wc+dw_div_wc, dw_div_wc 
         !      end if
         !   end do
         !end if


         if (s% i_lum>=0 .and. s% scale_max_correction_for_negative_surf_lum) then
            !ensure surface luminosity does not become negative
            dlum_surf = B(s% i_lum,1)*xscale(s% i_lum,1)
            old_lum_surf = s% xh_start(s% i_lum,1) + dx(s% i_lum,1)
            new_lum_surf = old_lum_surf + dlum_surf
            if (new_lum_surf < 0d0 .and. old_lum_surf > 0d0) then
               correction_factor = min(correction_factor, &
                  -s% max_frac_for_negative_surf_lum*old_lum_surf/dlum_surf)
            end if
         end if
         
         if (s% i_lnT > 0 .and. s% i_lnT <= nvar) &
            call clip1(s% i_lnT, s% solver_clip_dlogT*ln10)
         
         if (s% i_lnd > 0 .and. s% i_lnd <= nvar) &
            call clip1(s% i_lnd, s% solver_clip_dlogRho*ln10)
         
         if (s% i_lnR > 0 .and. s% i_lnR <= nvar) &
            call clip1(s% i_lnR, s% solver_clip_dlogR*ln10)

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
                  old_xa = s% xa_start(j,k) + dx(i,k)
                  if (old_xa <= 1d-90) cycle
                  dxa = B(i,k)*xscale(i,k)*correction_factor
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
               s% model_number, iter, min_alpha
         end if
         
         contains
            
         subroutine clip1(i, clip)
            integer, intent(in) :: i
            real(dp), intent(in) :: clip
            integer :: k
            real(dp) :: old_x, delta, abs_delta, abs_B
            include 'formats'
            if (clip <= 0d0) return
            do k = 1, s% nz
               old_x = s% xh_start(i,k) + dx(i,k) ! value before this iteration
               delta = B(i,k)*xscale(i,k)*correction_factor ! change for this iter
               ! skip if change small enough or if too big to change
               if (abs(delta) <= clip*abs(old_x) .or. is_bad(delta) .or. &
                   abs(old_x) < 1d0 .or. abs(delta) > 10d0*clip*abs(old_x)) cycle
               abs_delta = clip*abs(old_x)
               abs_B = abs_delta/(xscale(i,k)*correction_factor)
               B(i,k) = sign(abs_B,B(i,k))
               write(*,2) 'clip change ' // trim(s% nameofvar(i)), k, delta, old_x
               !stop 'Bdomain'
            end do
         end subroutine clip1

      end subroutine Bdomain


      subroutine inspectB(iter, nvar, nz, dx, B, xscale, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: iter, nvar, nz
         real(dp), pointer, dimension(:,:) :: dx, B, xscale ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: id
         integer, parameter :: inspectB_iter_stop = -1
         include 'formats'

         id = ipar(ipar_id)

         if (dbg) write(*, *) 'inspectB', iter
         ierr = 0

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !if (.not. s% solver_inspectB_flag) 
         return

         if (iter == inspectB_iter_stop) then
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
                  write(*,2) 'xscale ' // trim(s% nameofvar(j)), k, xscale(j, k)
                  write(*,2) 'dx ' // trim(s% nameofvar(j)), k, dx(j, k)
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
      integer function force_another_iteration(iter, itermin, lrpar, rpar, lipar, ipar)
         use hydro_mtx, only: ipar_id
         integer, intent(in) :: iter ! have finished this many iterations and have converged
         integer, intent(in) :: itermin ! this is the requested minimum.  iter may be < itermin.
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)

         type (star_info), pointer :: s
         integer :: id, ierr, k, res

         include 'formats'

         if (iter < itermin) then
            force_another_iteration = 1
            return
         end if
         force_another_iteration = 0

         id = ipar(ipar_id)
         ierr = 0

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then ! DISASTER
            force_another_iteration = -1
            return
         end if

         if (s% k_below_just_added > 1 .and. &
               s% num_surf_revisions < s% max_num_surf_revisions .and. &
               abs(s% lnS(1) - s% surf_lnS) > &
                  s% max_abs_rel_change_surf_lnS*max(s% lnS(1),s% surf_lnS)) then
            s% surf_lnT = s% lnT(1)
            s% surf_lnR = s% lnR(1)
            if (s% i_lnd /= 0) s% surf_lnd = s% lnd(1)
            if (s% i_v /= 0) s% surf_v = s% v(1)
            s% surf_lnS = s% lnS(1)
            s% num_surf_revisions = s% num_surf_revisions + 1
            force_another_iteration = 1
            s% used_extra_iter_in_solver_for_accretion = .true.
            return
         end if

      end function force_another_iteration


      end module hydro_solver_procs

