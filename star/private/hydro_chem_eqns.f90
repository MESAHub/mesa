! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module hydro_chem_eqns

      use star_private_def
      use const_def, only: dp
      use utils_lib

      implicit none

      private
      public :: do_chem_eqns, do1_chem_eqns

      contains

      subroutine do_chem_eqns(s, nvar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         integer :: k, op_err
         include 'formats'
         ierr = 0
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k = 1, s% nz
            if (ierr /= 0) cycle
            call do1_chem_eqns(s, k, nvar, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
      end subroutine do_chem_eqns


      subroutine do1_chem_eqns(s, k, nvar, ierr)

         use chem_def
         use net_lib, only: show_net_reactions, show_net_params
         use rates_def, only: i_rate
         use star_utils, only: em1, e00, ep1
         use auto_diff_support
         use implicit_Dmix, only: implicit_Dmix_debug_selected_cell

         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr

         !integer, pointer :: reaction_id(:) ! maps net reaction number to reaction id
         integer :: nz, j, i, jj, ii, species, debug_max_i_var, debug_max_loc, &
            debug_dxa_max_jj, debug_dxa_max_loc
         real(dp) :: &
            dxdt_expected_dxa, dxdt_expected, dxdt_actual, &
            dq, dm, dequ, dxdt_nuc, dxdt_mix, max_abs_residual, &
            sum_dxdt_nuc, &
            d_dxdt_mix_dx00, d_dxdt_mix_dxm1, d_dxdt_mix_dxp1, &
            sum_dx_burning, sum_dx_mixing, residual, &
            dxdt_factor, eqn_scale, &
            dequ_dlnd, dequ_dlnT, dx00, dxp1, debug_max_dequ_abs, &
            debug_max_dequ_val, debug_dxa_max_abs, debug_dxa_max_val
         character (len=16) :: debug_max_var
         type(auto_diff_real_star_order1) :: &
            sig00_ad, sigp1_ad
         logical :: test_partials, doing_op_split_burn
         logical :: do_implicit_dsig_structure, do_implicit_dsig_dxa
         logical, parameter :: checking = .false.

         include 'formats'

         ierr = 0

         species = s% species
         nz = s% nz
         do_implicit_dsig_structure = s% job% implicit_diffusion_flag .and. &
            s% job% implicit_diffusion_include_dsig_structure
         do_implicit_dsig_dxa = s% job% implicit_diffusion_flag .and. &
            s% use_Ledoux_criterion .and. &
            s% job% implicit_diffusion_include_dsig_dxa

         dq = s% dq(k)
         dm = s% dm(k)

         max_abs_residual = 0
         sum_dxdt_nuc = 0

         if (s% do_mix) then

            d_dxdt_mix_dxm1 = s% d_dxdt_mix_dxm1(k)
            d_dxdt_mix_dx00 = s% d_dxdt_mix_dx00(k)
            d_dxdt_mix_dxp1 = s% d_dxdt_mix_dxp1(k)

         else

            d_dxdt_mix_dxm1 = 0
            d_dxdt_mix_dx00 = 0
            d_dxdt_mix_dxp1 = 0

         end if

         sum_dx_burning = 0
         sum_dx_mixing = 0

         do j=1,species  ! composition equation for species j in cell k

            !test_partials = (k == s% solver_test_partials_k .and. s% net_iso(ihe4) == j)
            test_partials = .false.

            i = s%nvar_hydro+j

            dxdt_actual = s% xa_sub_xa_start(j,k)/s% dt

            doing_op_split_burn = s% op_split_burn .and. s% T_start(k) >= s% op_split_burn_min_T
            if (s% do_burn .and. .not. doing_op_split_burn) then
               dxdt_nuc = s% dxdt_nuc(j,k)
            else
               dxdt_nuc = 0
            end if

            if (s% do_mix) then
               dxdt_mix = s% dxdt_mix(j,k)
            else
               dxdt_mix = 0
            end if

            dxdt_expected = dxdt_mix + dxdt_nuc

            dxdt_factor = 1d0

            eqn_scale = max(s% min_chem_eqn_scale, s% x_scale(i,k)/s% dt)
            residual = (dxdt_expected - dxdt_actual)/eqn_scale
            s% equ(i,k) = residual

            if (abs(residual) > max_abs_residual) &
               max_abs_residual = abs(s% equ(i,k))

            if (is_bad(s% equ(i,k))) then
               s% retry_message = 'bad residual for do1_chem_eqns'
!$OMP critical (star_chem_eqns_bad_num)
               if (s% report_ierr) then
                  write(*,3) 'do1_chem_eqns: equ ' // trim(s% nameofequ(i)), &
                        i, k, s% equ(i,k)
                  write(*,2) 'dxdt_expected', k, dxdt_expected
                  write(*,2) 'dxdt_actual', k, dxdt_actual
                  write(*,2) 'eqn_scale', k, eqn_scale
                  write(*,2) 'dxdt_mix', k, dxdt_mix
                  write(*,2) 'dxdt_nuc', k, dxdt_nuc
               end if
               if (s% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'do1_chem_eqns')
!$OMP end critical (star_chem_eqns_bad_num)
            end if

            call e00(s, i, i, k, nvar, -1d0/eqn_scale/s% dt)

            ! all the rest are jacobian terms for dxdt_expected/eqn_scale

            if (s% do_burn) then

               do jj=1,species
                  ii = s% nvar_hydro+jj
                  dxdt_expected_dxa = s% d_dxdt_nuc_dx(j,jj,k)
                  dequ = dxdt_expected_dxa/eqn_scale
                  if (checking) call check_dequ(dequ,'d_dxdt_nuc_dx')
                  call e00(s, i, ii, k, nvar, dxdt_factor*dequ)
               end do

               dequ_dlnd = s% d_dxdt_nuc_drho(j,k)*s% rho(k)/eqn_scale
               call e00(s, i, s% i_lnd, k, nvar, dxdt_factor*dequ_dlnd)

               dequ_dlnT = s% d_dxdt_nuc_dT(j,k)*s% T(k)/eqn_scale
               call e00(s, i, s% i_lnT, k, nvar, dxdt_factor*dequ_dlnT)

            end if

            if (s% do_mix) then

               dxdt_expected_dxa = d_dxdt_mix_dx00
               dequ = dxdt_expected_dxa/eqn_scale
               if (checking) call check_dequ(dequ,'d_dxdt_mix_dx00')
               call e00(s, i, i, k, nvar, dxdt_factor*dequ)
               if (k > 1) then
                  dxdt_expected_dxa = d_dxdt_mix_dxm1
                  dequ = dxdt_expected_dxa/eqn_scale
                  if (checking) call check_dequ(dequ,'d_dxdt_mix_dxm1')
                  call em1(s, i, i, k, nvar, dxdt_factor*dequ)
               end if
               if (k < nz) then
                  dxdt_expected_dxa = d_dxdt_mix_dxp1
                  dequ = dxdt_expected_dxa/eqn_scale
                  if (checking) call check_dequ(dequ,'d_dxdt_mix_dxp1')
                  call ep1(s, i, i, k, nvar, dxdt_factor*dequ)
               end if

               if (do_implicit_dsig_structure .or. do_implicit_dsig_dxa) then
                  if (k > 1) then
                     dx00 = s% xa(j,k-1) - s% xa(j,k)
                     if (do_implicit_dsig_structure) sig00_ad = s% sig_implicit_ad(k)
                  else
                     dx00 = 0d0
                     if (do_implicit_dsig_structure) sig00_ad = 0d0
                  end if
                  if (k < nz) then
                     dxp1 = s% xa(j,k) - s% xa(j,k+1)
                     if (do_implicit_dsig_structure) &
                        sigp1_ad = shift_p1(s% sig_implicit_ad(k+1))
                  else
                     dxp1 = 0d0
                     if (do_implicit_dsig_structure) sigp1_ad = 0d0
                  end if
                  if (do_implicit_dsig_structure) then
                     call debug_sig_ad_partials( &
                        dx00/(dm*eqn_scale), -dxp1/(dm*eqn_scale))
                     call save_sig_ad_partials( &
                        dx00/(dm*eqn_scale), -dxp1/(dm*eqn_scale))
                  end if
                  if (do_implicit_dsig_dxa) then
                     call debug_sig_dxa_partials
                     do jj=1,species
                        ii = s% nvar_hydro+jj
                        if (k > 1) then
                           dequ = s% d_sig_dxa_m1(jj,k)*dx00/(dm*eqn_scale)
                           if (checking) call check_dequ(dequ,'d_sig_dxa_m1')
                           call em1(s, i, ii, k, nvar, dxdt_factor*dequ)
                        end if
                        dequ = 0d0
                        if (k > 1) dequ = dequ + &
                           s% d_sig_dxa_00(jj,k)*dx00/dm
                        if (k < nz) dequ = dequ - &
                           s% d_sig_dxa_m1(jj,k+1)*dxp1/dm
                        dequ = dequ/eqn_scale
                        if (checking) call check_dequ(dequ,'d_sig_dxa_00')
                        call e00(s, i, ii, k, nvar, dxdt_factor*dequ)
                        if (k < nz) then
                           dequ = -s% d_sig_dxa_00(jj,k+1)*dxp1/(dm*eqn_scale)
                           if (checking) call check_dequ(dequ,'d_sig_dxa_p1')
                           call ep1(s, i, ii, k, nvar, dxdt_factor*dequ)
                        end if
                     end do
                  end if
               end if

            end if

            if (test_partials) then
               s% solver_test_partials_dx_sink = s% net_iso(img24)
               s% solver_test_partials_val = s% dxdt_nuc(j,k)
               s% solver_test_partials_var = s% nvar_hydro + j
               s% solver_test_partials_dval_dx = s% d_dxdt_nuc_dx(j,j,k)
               write(*,*) 'do1_chem_eqns', s% solver_test_partials_var
            end if

         end do

         contains

         subroutine check_dequ(dequ, str)
            real(dp), intent(in) :: dequ
            character (len=*), intent(in) :: str
            include 'formats'
            if (is_bad(dequ)) then
               ierr = -1
               if (s% report_ierr) then
                  write(*,2) 'do1_chem_eqns: bad ' // trim(str), k
               end if
               if (s% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'do1_chem_eqns')
               return
            end if
         end subroutine check_dequ


         subroutine debug_sig_ad_partials(sig00_factor, sigp1_factor)
            real(dp), intent(in) :: sig00_factor, sigp1_factor

            real(dp) :: sig00_d1_max, sigp1_d1_max, xa_m1, xa_p1

            if (.not. implicit_Dmix_debug_selected_cell(s, k)) return

            sig00_d1_max = maxval(abs(sig00_ad%d1Array))
            sigp1_d1_max = maxval(abs(sigp1_ad%d1Array))
            debug_max_dequ_abs = 0d0
            debug_max_dequ_val = 0d0
            debug_max_i_var = 0
            debug_max_loc = 0
            debug_max_var = 'none'

            call debug_sig_ad_var('lnd', s% i_lnd, i_lnd_m1, i_lnd_00, i_lnd_p1, &
               sig00_factor, sigp1_factor)
            call debug_sig_ad_var('lnT', s% i_lnT, i_lnT_m1, i_lnT_00, i_lnT_p1, &
               sig00_factor, sigp1_factor)
            call debug_sig_ad_var('lnR', s% i_lnR, i_lnR_m1, i_lnR_00, i_lnR_p1, &
               sig00_factor, sigp1_factor)
            if (s% i_v /= 0) call debug_sig_ad_var( &
               'v', s% i_v, i_v_m1, i_v_00, i_v_p1, sig00_factor, sigp1_factor)
            if (s% i_u /= 0) call debug_sig_ad_var( &
               'u', s% i_u, i_v_m1, i_v_00, i_v_p1, sig00_factor, sigp1_factor)
            if (s% i_lum /= 0) call debug_sig_ad_var( &
               'lum', s% i_lum, i_L_m1, i_L_00, i_L_p1, sig00_factor, sigp1_factor)
            if (s% i_w /= 0) call debug_sig_ad_var( &
               'w', s% i_w, i_w_m1, i_w_00, i_w_p1, sig00_factor, sigp1_factor)
            if (s% i_Hp /= 0) call debug_sig_ad_var( &
               'Hp', s% i_Hp, i_Hp_m1, i_Hp_00, i_Hp_p1, sig00_factor, sigp1_factor)
            if (s% i_w_div_wc /= 0) call debug_sig_ad_var( &
               'w_div_wc', s% i_w_div_wc, i_w_div_wc_m1, i_w_div_wc_00, &
               i_w_div_wc_p1, sig00_factor, sigp1_factor)
            if (s% i_j_rot /= 0) call debug_sig_ad_var( &
               'j_rot', s% i_j_rot, i_jrot_m1, i_jrot_00, i_jrot_p1, &
               sig00_factor, sigp1_factor)

            if (debug_max_dequ_abs == 0d0) return

            xa_m1 = -99d0
            xa_p1 = -99d0
            if (k > 1) xa_m1 = s% xa(j,k-1)
            if (k < nz) xa_p1 = s% xa(j,k+1)

!$OMP critical (implicit_Dmix_debug_chem)
            write(*,*) 'implicit_Dmix debug chem dsig_struct model', &
               s% model_number, ' iter', s% solver_iter, &
               ' k', k, ' j', j, &
               ' isotope', trim(chem_isos% name(s% chem_id(j))), &
               ' mix_type', s% mixing_type(k), &
               ' mlt_type', s% mlt_mixing_type(k)
            write(*,*) '   residual', residual, ' eqn_scale', eqn_scale, &
               ' dxdt_actual', dxdt_actual, ' dxdt_mix', dxdt_mix, &
               ' dxdt_nuc', dxdt_nuc, &
               ' xa_m1', xa_m1, ' xa_00', s% xa(j,k), &
               ' xa_p1', xa_p1, ' xa_start_00', s% xa_start(j,k)
            write(*,*) '   dx00', dx00, ' dxp1', dxp1, ' dm', dm, &
               ' sig00_factor', sig00_factor, &
               ' sigp1_factor', sigp1_factor, &
               ' sig00_val', sig00_ad%val, ' sigp1_val', sigp1_ad%val
            write(*,*) '   sig00_d1_max', sig00_d1_max, &
               ' sigp1_d1_max', sigp1_d1_max, &
               ' max_dsig_struct_dequ', debug_max_dequ_val, &
               ' max_abs_dsig_struct_dequ', debug_max_dequ_abs, &
               ' max_var', trim(debug_max_var), &
               ' max_i_var', debug_max_i_var, &
               ' max_loc', debug_max_loc
!$OMP end critical (implicit_Dmix_debug_chem)
         end subroutine debug_sig_ad_partials


         subroutine debug_sig_dxa_partials
            integer :: jj
            real(dp) :: dequ_m1, dequ_00, dequ_p1

            if (.not. implicit_Dmix_debug_selected_cell(s, k)) return

            debug_dxa_max_abs = 0d0
            debug_dxa_max_val = 0d0
            debug_dxa_max_jj = 0
            debug_dxa_max_loc = 0

            do jj=1,species
               if (k > 1) then
                  dequ_m1 = s% d_sig_dxa_m1(jj,k)*dx00/(dm*eqn_scale)
                  call debug_sig_dxa_candidate(jj, -1, dequ_m1)
               end if
               dequ_00 = 0d0
               if (k > 1) dequ_00 = dequ_00 + &
                  s% d_sig_dxa_00(jj,k)*dx00/dm
               if (k < nz) dequ_00 = dequ_00 - &
                  s% d_sig_dxa_m1(jj,k+1)*dxp1/dm
               dequ_00 = dequ_00/eqn_scale
               call debug_sig_dxa_candidate(jj, 0, dequ_00)
               if (k < nz) then
                  dequ_p1 = -s% d_sig_dxa_00(jj,k+1)*dxp1/(dm*eqn_scale)
                  call debug_sig_dxa_candidate(jj, 1, dequ_p1)
               end if
            end do

            if (debug_dxa_max_abs == 0d0) return

!$OMP critical (implicit_Dmix_debug_chem)
            write(*,*) 'implicit_Dmix debug chem dsig_dxa model', &
               s% model_number, ' iter', s% solver_iter, &
               ' k', k, ' j', j, &
               ' isotope', trim(chem_isos% name(s% chem_id(j))), &
               ' mix_type', s% mixing_type(k), &
               ' mlt_type', s% mlt_mixing_type(k)
            write(*,*) '   residual', residual, ' eqn_scale', eqn_scale, &
               ' dxdt_actual', dxdt_actual, ' dxdt_mix', dxdt_mix, &
               ' dxdt_nuc', dxdt_nuc, &
               ' dx00', dx00, ' dxp1', dxp1, ' dm', dm
            write(*,*) '   max_dsig_dxa_dequ', debug_dxa_max_val, &
               ' max_abs_dsig_dxa_dequ', debug_dxa_max_abs, &
               ' max_jj', debug_dxa_max_jj, &
               ' max_partial_isotope', &
               trim(chem_isos% name(s% chem_id(debug_dxa_max_jj))), &
               ' max_loc', debug_dxa_max_loc
!$OMP end critical (implicit_Dmix_debug_chem)
         end subroutine debug_sig_dxa_partials


         subroutine debug_sig_dxa_candidate(jj_in, loc, dequ)
            integer, intent(in) :: jj_in, loc
            real(dp), intent(in) :: dequ

            if (abs(dequ) <= debug_dxa_max_abs) return
            debug_dxa_max_abs = abs(dequ)
            debug_dxa_max_val = dequ
            debug_dxa_max_jj = jj_in
            debug_dxa_max_loc = loc
         end subroutine debug_sig_dxa_candidate


         subroutine debug_sig_ad_var(var_name, i_var, i_m1, i_00, i_p1, &
               sig00_factor, sigp1_factor)
            character (len=*), intent(in) :: var_name
            integer, intent(in) :: i_var, i_m1, i_00, i_p1
            real(dp), intent(in) :: sig00_factor, sigp1_factor

            real(dp) :: d_m1, d_00, d_p1

            if (i_var == 0) return

            d_m1 = sig00_factor*sig00_ad%d1Array(i_m1) + &
               sigp1_factor*sigp1_ad%d1Array(i_m1)
            d_00 = sig00_factor*sig00_ad%d1Array(i_00) + &
               sigp1_factor*sigp1_ad%d1Array(i_00)
            d_p1 = sig00_factor*sig00_ad%d1Array(i_p1) + &
               sigp1_factor*sigp1_ad%d1Array(i_p1)

            call debug_sig_ad_candidate(var_name, i_var, -1, d_m1)
            call debug_sig_ad_candidate(var_name, i_var, 0, d_00)
            call debug_sig_ad_candidate(var_name, i_var, 1, d_p1)
         end subroutine debug_sig_ad_var


         subroutine debug_sig_ad_candidate(var_name, i_var, loc, dequ)
            character (len=*), intent(in) :: var_name
            integer, intent(in) :: i_var, loc
            real(dp), intent(in) :: dequ

            if (abs(dequ) <= debug_max_dequ_abs) return
            debug_max_dequ_abs = abs(dequ)
            debug_max_dequ_val = dequ
            debug_max_i_var = i_var
            debug_max_loc = loc
            debug_max_var = var_name
         end subroutine debug_sig_ad_candidate


         subroutine save_sig_ad_partials(sig00_factor, sigp1_factor)
            real(dp), intent(in) :: sig00_factor, sigp1_factor

            call add_sig_ad_var(s% i_lnd, i_lnd_m1, i_lnd_00, i_lnd_p1, &
               sig00_factor, sigp1_factor)
            call add_sig_ad_var(s% i_lnT, i_lnT_m1, i_lnT_00, i_lnT_p1, &
               sig00_factor, sigp1_factor)
            call add_sig_ad_var(s% i_lnR, i_lnR_m1, i_lnR_00, i_lnR_p1, &
               sig00_factor, sigp1_factor)
            if (s% i_v /= 0) call add_sig_ad_var( &
               s% i_v, i_v_m1, i_v_00, i_v_p1, sig00_factor, sigp1_factor)
            if (s% i_u /= 0) call add_sig_ad_var( &
               s% i_u, i_v_m1, i_v_00, i_v_p1, sig00_factor, sigp1_factor)
            if (s% i_lum /= 0) call add_sig_ad_var( &
               s% i_lum, i_L_m1, i_L_00, i_L_p1, sig00_factor, sigp1_factor)
            if (s% i_w /= 0) call add_sig_ad_var( &
               s% i_w, i_w_m1, i_w_00, i_w_p1, sig00_factor, sigp1_factor)
            if (s% i_Hp /= 0) call add_sig_ad_var( &
               s% i_Hp, i_Hp_m1, i_Hp_00, i_Hp_p1, sig00_factor, sigp1_factor)
            if (s% i_w_div_wc /= 0) call add_sig_ad_var( &
               s% i_w_div_wc, i_w_div_wc_m1, i_w_div_wc_00, i_w_div_wc_p1, &
               sig00_factor, sigp1_factor)
            if (s% i_j_rot /= 0) call add_sig_ad_var( &
               s% i_j_rot, i_jrot_m1, i_jrot_00, i_jrot_p1, &
               sig00_factor, sigp1_factor)

         end subroutine save_sig_ad_partials


         subroutine add_sig_ad_var(i_var, i_m1, i_00, i_p1, sig00_factor, sigp1_factor)
            integer, intent(in) :: i_var, i_m1, i_00, i_p1
            real(dp), intent(in) :: sig00_factor, sigp1_factor

            real(dp) :: d_m1, d_00, d_p1

            if (i_var == 0) return

            d_m1 = sig00_factor*sig00_ad%d1Array(i_m1) + &
               sigp1_factor*sigp1_ad%d1Array(i_m1)
            d_00 = sig00_factor*sig00_ad%d1Array(i_00) + &
               sigp1_factor*sigp1_ad%d1Array(i_00)
            d_p1 = sig00_factor*sig00_ad%d1Array(i_p1) + &
               sigp1_factor*sigp1_ad%d1Array(i_p1)

            if (d_m1 /= 0d0) call em1(s, i, i_var, k, nvar, d_m1)
            if (d_00 /= 0d0) call e00(s, i, i_var, k, nvar, d_00)
            if (d_p1 /= 0d0) call ep1(s, i, i_var, k, nvar, d_p1)

         end subroutine add_sig_ad_var

      end subroutine do1_chem_eqns

      end module hydro_chem_eqns
