! ***********************************************************************
!
!   Copyright (C) 2010-2020  The MESA Team
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

      module hydro_rsp2

      use star_private_def
      use const_def, only: dp, boltz_sigma, pi, clight, crad, ln10
      use utils_lib, only: is_bad
      use auto_diff
      use auto_diff_support
      use accurate_sum_auto_diff_star_order1
      use star_utils

      implicit none

      private
      public :: do1_rsp2_L_eqn
      public :: do1_turbulent_energy_eqn
      public :: do1_rsp2_Hp_eqn
      public :: compute_Eq_cell
      public :: compute_Uq_face
      public :: set_RSP2_vars
      public :: Hp_face_for_rsp2_val
      public :: Hp_face_for_rsp2_eqn, set_etrb_start_vars
      public :: RSP2_adjust_vars_before_call_solver
      public :: get_RSP2_alfa_beta_face_weights

      real(dp), parameter :: &
         x_ALFAP = 2.d0/3.d0, &  ! Ptrb
         x_ALFAS = (1.d0/2.d0)*sqrt_2_div_3, &  ! PII_face and Lc
         x_ALFAC = (1.d0/2.d0)*sqrt_2_div_3, &  ! Lc
         x_CEDE  = (8.d0/3.d0)*sqrt_2_div_3, &  ! DAMP
         x_GAMMAR = 2.d0*sqrt(3.d0)  ! DAMPR

      contains

      subroutine set_RSP2_vars(s,ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: x
         integer :: k, op_err
         include 'formats'
         ierr = 0
         op_err = 0
         !$OMP PARALLEL DO PRIVATE(k,op_err,x) SCHEDULE(dynamic,2)
         do k=1,s%nz
            op_err = 0
            ! Hp_face(k) <= 0 means it needs to be set.  e.g., after read file
            if (s% Hp_face(k) <= 0) then
               s% Hp_face(k) = get_scale_height_face_val(s,k)
               s% xh(s% i_Hp,k) = s% Hp_face(k)
            end if
            x = compute_RSP2_gradT(s, k, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if
            s% gradT_ad(k) = x
            s% gradT(k) = x%val
            op_err = 0
            x = compute_Y_face(s, k, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if
            op_err = 0
            x = compute_PII_face(s, k, op_err)
            if (op_err /= 0) ierr = op_err
            !Pvsc           skip?
         end do
         !$OMP END PARALLEL DO
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'failed in set_RSP2_vars loop 1', s% model_number
            return
         end if
         !$OMP PARALLEL DO PRIVATE(k,op_err,x) SCHEDULE(dynamic,2)
         do k=1,s% nz
            op_err = 0
            x = compute_Chi_cell(s, k, op_err)
            if (op_err /= 0) ierr = op_err
            op_err = 0
            x = compute_Eq_cell(s, k, op_err)
            if (op_err /= 0) ierr = op_err
            op_err = 0
            x = compute_C(s, k, op_err)  ! COUPL
            if (op_err /= 0) ierr = op_err
            op_err = 0
            x = compute_L_face(s, k, op_err)  ! Lr, Lt, Lc
            if (op_err /= 0) ierr = op_err
            op_err = 0
            x = compute_Lrad_coeff(s, k, op_err)
            if (op_err /= 0) then
               ierr = op_err
            else if (k > 1 .and. x%val > 0d0) then
               s% gradr_ad(k) = (wrap_L_00(s,k) - s% Lt_ad(k))/x
               s% gradr(k) = s% gradr_ad(k)%val
            end if
         end do
         !$OMP END PARALLEL DO
         if (ierr /= 0) then
            if (s% report_ierr) write(*,2) 'failed in set_RSP2_vars loop 2', s% model_number
            return
         end if
         do k = 1, s% RSP2_num_outermost_cells_forced_nonturbulent
            s% Eq(k) = 0d0; s% Eq_ad(k) = 0d0
            s% Chi(k) = 0d0; s% Chi_ad(k) = 0d0
            s% COUPL(k) = 0d0; s% COUPL_ad(k) = 0d0
            !s% Ptrb(k) = 0d0;
            s% Lc(k) = 0d0; s% Lc_ad(k) = 0d0
            s% Lt(k) = 0d0; s% Lt_ad(k) = 0d0
         end do
         do k = s% nz + 1 - int(s% nz/s% RSP2_nz_div_IBOTOM) , s% nz
            s% Eq(k) = 0d0; s% Eq_ad(k) = 0d0
            s% Chi(k) = 0d0; s% Chi_ad(k) = 0d0
            s% COUPL(k) = 0d0; s% COUPL_ad(k) = 0d0
            !s% Ptrb(k) = 0d0;
            s% Lc(k) = 0d0; s% Lc_ad(k) = 0d0
            s% Lt(k) = 0d0; s% Lt_ad(k) = 0d0
         end do
      end subroutine set_RSP2_vars


      subroutine do1_rsp2_L_eqn(s, k, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) ::  &
            L_expected, L_actual,resid
         type(accurate_auto_diff_real_star_order1) :: L_sum
         real(dp) :: scale, residual
         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         if (.not. s% RSP2_flag) then
            ierr = -1
            return
         end if

         ierr = 0
         if (k /= 1) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'do1_rsp2_L_eqn requires k == 1', k
            return
         end if
         L_sum = s% Lr_ad(k)
         L_sum = L_sum + s% Lc_ad(k)
         L_sum = L_sum + s% Lt_ad(k)
         L_expected = L_sum
         L_actual = wrap_L_00(s, k)
         scale = 1d0/max(abs(L_expected%val),abs(L_actual%val),1d0)
         if (is_bad(scale)) then
            write(*,2) 'do1_rsp2_L_eqn scale', k, scale
            call mesa_error(__FILE__,__LINE__,'do1_rsp2_L_eqn')
         end if
         resid = (L_expected - L_actual)*scale

         residual = resid%val
         s% equ(s% i_equL, k) = residual
         if (test_partials) then
            s% solver_test_partials_val = residual
         end if

         call save_eqn_residual_info(s, k, nvar, s% i_equL, resid, 'do1_rsp2_L_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = resid%d1Array(i_lnR_00)
            write(*,4) 'do1_rsp2_L_eqn', s% solver_test_partials_var
         end if
      end subroutine do1_rsp2_L_eqn


      subroutine do1_rsp2_Hp_eqn(s, k, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) ::  &
            Hp_expected, Hp_actual,resid
         real(dp) :: residual, Hp_start
         logical :: test_partials
         include 'formats'
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         if (.not. s% RSP2_flag) then
            ierr = -1
            return
         end if

         ierr = 0
         Hp_expected = Hp_face_for_rsp2_eqn(s, k, ierr)
         if (ierr /= 0) return
         Hp_actual = wrap_Hp_00(s, k)
         Hp_start = s% Hp_face_start(k)
         resid = (Hp_expected - Hp_actual)/max(Hp_expected,Hp_actual)

         residual = resid%val
         s% equ(s% i_equ_Hp, k) = residual
         if (test_partials) then
            s% solver_test_partials_val = residual
         end if

         if (residual > 1d3) then
         !$omp critical (hydro_rsp2_1)
            write(*,2) 'residual', k, residual
            write(*,2) 'Hp_expected', k, Hp_expected%val
            write(*,2) 'Hp_actual', k, Hp_actual%val
            call mesa_error(__FILE__,__LINE__,'do1_rsp2_Hp_eqn')
         !$omp end critical (hydro_rsp2_1)
         end if

         call save_eqn_residual_info(s, k, nvar, s% i_equ_Hp, resid, 'do1_rsp2_Hp_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = resid%d1Array(i_lnR_00)
            write(*,4) 'do1_rsp2_Hp_eqn', s% solver_test_partials_var
         end if

      end subroutine do1_rsp2_Hp_eqn


      real(dp) function Hp_face_for_rsp2_val(s, k, ierr) result(Hp_face)  ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_face_ad
         ierr = 0
         Hp_face_ad = Hp_face_for_rsp2_eqn(s, k, ierr)
         if (ierr /= 0) return
         Hp_face = Hp_face_ad%val
      end function Hp_face_for_rsp2_val


      function Hp_face_for_rsp2_eqn(s, k, ierr) result(Hp_face)  ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_face
         type(auto_diff_real_star_order1) :: &
            rho_face, area, dlnPeos, &
            r_00, Peos_00, d_00, Peos_m1, d_m1, Peos_div_rho, &
            d_face, Peos_face, alt_Hp_face, A
         real(dp) :: alfa, beta
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            Hp_face = 1d0  ! not used
            return
         end if
         if (k > 1 .and. .not. s% RSP2_assume_HSE) then
            call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
            rho_face = alfa*wrap_d_00(s,k) + beta*wrap_d_m1(s,k)
            area = 4d0*pi*pow2(wrap_r_00(s,k))
            dlnPeos = wrap_lnPeos_m1(s,k) - wrap_lnPeos_00(s,k)
            Hp_face = -s% dm_bar(k)/(area*rho_face*dlnPeos)
         else
            r_00 = wrap_r_00(s, k)  ! not time-centered in RSP
            d_00 = wrap_d_00(s, k)
            Peos_00 = wrap_Peos_00(s, k)
            if (k == 1) then
               Peos_div_rho = Peos_00/d_00
               Hp_face = pow2(r_00)*Peos_div_rho/(s% cgrav(k)*s% m(k))
            else
               d_m1 = wrap_d_m1(s, k)
               Peos_m1 = wrap_Peos_m1(s, k)
               call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
               Peos_div_rho = alfa*Peos_00/d_00 + beta*Peos_m1/d_m1
               Hp_face = pow2(r_00)*Peos_div_rho/(s% cgrav(k)*s% m(k))
               if (k==-104) then
                  write(*,3) 'RSP2 Hp P_div_rho Pdrho_00 Pdrho_m1', k, s% solver_iter, &
                     Hp_face%val, Peos_div_rho%val, Peos_00%val/d_00%val, Peos_m1%val/d_m1%val
                  !write(*,3) 'RSP2 Hp r2_div_Gm r_start r', k, s% solver_iter, &
                  !   Hp_face%val, pow2(r_00%val)/(s% cgrav(k)*s% m(k)), &
                  !   s% r_start(k), r_00%val
               end if
               if (s% alt_scale_height_flag) then
                  call mesa_error(__FILE__,__LINE__,'Hp_face_for_rsp2_eqn: cannot use alt_scale_height_flag')
                  ! consider sound speed*hydro time scale as an alternative scale height
                  d_face = alfa*d_00 + beta*d_m1
                  Peos_face = alfa*Peos_00 + beta*Peos_m1
                  alt_Hp_face = sqrt(Peos_face/s% cgrav(k))/d_face
                  if (alt_Hp_face%val < Hp_face%val) then  ! blend
                     A = pow2(alt_Hp_face/Hp_face)  ! 0 <= A%val < 1
                     Hp_face = A*Hp_face + (1d0 - A)*alt_Hp_face
                  end if
               end if
            end if
         end if
      end function Hp_face_for_rsp2_eqn


      subroutine do1_turbulent_energy_eqn(s, k, nvar, ierr)
         use star_utils, only: calc_Ptrb_ad_tw, set_energy_eqn_scal, save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            d_turbulent_energy_ad, Ptrb_dV_ad, dt_C_ad, dt_Eq_ad, &
            source_div_w_ad, D_div_w_ad, Dr_div_w_ad, Eq_div_w_ad, C_div_w_ad, &
            Ptrb_dV_div_w_ad
         type(auto_diff_real_star_order1) :: w_00, Ptrb_div_etrb, dV_ad
         type(auto_diff_real_star_order1) :: tst, resid_ad, dt_dLt_dm_ad
         type(accurate_auto_diff_real_star_order1) :: esum_ad
         logical :: non_turbulent_cell, positive_branch, test_partials
         real(dp) :: residual, scal, P_theta
         include 'formats'
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         ierr = 0
         w_00 = wrap_w_00(s,k)

         non_turbulent_cell = &
            s% mixing_length_alpha == 0d0 .or. &
            k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
            k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)
         if (.not. s% RSP2_flag) then
            resid_ad = w_00 - s% w_start(k)  ! just hold w constant when not using RSP2
         else if (non_turbulent_cell) then
            resid_ad = w_00/s% csound(k)  ! make w = 0
         else
            call set_energy_eqn_scal(s, k, scal, ierr); if (ierr /= 0) return  ! 1/(erg g^-1 s^-1)
            positive_branch = &
               s% RSP2_source_seed == 0d0 .and. s% RSP2_alfat == 0d0 .and. &
               s% w_start(k) == 0d0
            if (positive_branch) then
               source_div_w_ad = compute_Source_div_w(s, k, ierr)
               if (ierr /= 0) return
               Eq_div_w_ad = compute_Eq_div_w_cell(s, k, ierr)
               if (ierr /= 0) return
               positive_branch = source_div_w_ad%val + Eq_div_w_ad%val > 0d0
            end if

            if (positive_branch) then
               ! The unstable zero branch factors from the local residual as w*G(w).
               call calc_Ptrb_ad_tw(s, k, Ptrb_dV_ad, Ptrb_div_etrb, ierr)
               if (ierr /= 0) return
               dV_ad = 1d0/wrap_d_00(s,k) - 1d0/s% rho_start(k)
               P_theta = 1d0
               if (s% using_velocity_time_centering .and. &
                     s% include_P_in_velocity_time_centering) &
                  P_theta = s% P_theta_for_velocity_time_centering
               Ptrb_dV_div_w_ad = P_theta*Ptrb_div_etrb*w_00*dV_ad
               D_div_w_ad = compute_D_div_w(s, k, ierr)
               if (ierr /= 0) return
               Dr_div_w_ad = compute_Dr_div_w(s, k, ierr)
               if (ierr /= 0) return
               C_div_w_ad = source_div_w_ad - D_div_w_ad - Dr_div_w_ad

               esum_ad = w_00
               esum_ad = esum_ad + Ptrb_dV_div_w_ad
               esum_ad = esum_ad - s% dt*C_div_w_ad
               esum_ad = esum_ad - s% dt*Eq_div_w_ad
               ! Scale G as an energy row without restoring the vanishing factor w.
               resid_ad = esum_ad*s% csound(k)*scal/s% dt
            else
               call setup_d_turbulent_energy(ierr); if (ierr /= 0) return  ! erg g^-1 = cm^2 s^-2
               call setup_Ptrb_dV_ad(ierr); if (ierr /= 0) return  ! erg g^-1
               call setup_dt_dLt_dm_ad(ierr); if (ierr /= 0) return  ! erg g^-1
               call setup_dt_C_ad(ierr); if (ierr /= 0) return  ! erg g^-1
               call setup_dt_Eq_ad(ierr); if (ierr /= 0) return  ! erg g^-1
               ! sum terms in esum_ad using accurate_auto_diff_real_star_order1
               esum_ad = d_turbulent_energy_ad
               esum_ad = esum_ad + Ptrb_dV_ad
               esum_ad = esum_ad + dt_dLt_dm_ad
               esum_ad = esum_ad - dt_C_ad
               esum_ad = esum_ad - dt_Eq_ad  ! erg g^-1
               resid_ad = esum_ad*scal/s%dt
            end if

            if (.not. positive_branch .and. k==-35 .and. s% solver_iter == 1) then
                  write(*,3) 'RSP2 w dEt PdV dtC dtEq', k, s% solver_iter, &
                     w_00%val, d_turbulent_energy_ad%val, Ptrb_dV_ad%val, dt_C_ad%val, dt_Eq_ad%val
            end if

         end if

         residual = resid_ad%val
         s% equ(s% i_detrb_dt, k) = residual

         if (test_partials) then
            tst = residual
            s% solver_test_partials_val = tst%val
            if (s% solver_iter == 12) &
               write(*,*) 'do1_turbulent_energy_eqn', s% solver_test_partials_var, s% lnd(k), tst%val
         end if

         call save_eqn_residual_info(s, k, nvar, s% i_detrb_dt, resid_ad, 'do1_turbulent_energy_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = tst%d1Array(i_lnd_00)     ! xi0 good , xi1 partial 0, xi2 good.  Af horrible.'
            write(*,*) 'do1_turbulent_energy_eqn', s% solver_test_partials_var, s% lnd(k)/ln10, tst%val
         end if

         contains

         subroutine setup_d_turbulent_energy(ierr)  ! erg g^-1
            integer, intent(out) :: ierr
            ierr = 0
            d_turbulent_energy_ad = wrap_etrb_00(s,k) - get_etrb_start(s,k)
         end subroutine setup_d_turbulent_energy

         ! Ptrb_dV_ad = Ptrb_ad*dV_ad
         subroutine setup_Ptrb_dV_ad(ierr)  ! erg g^-1
            use star_utils, only: calc_Ptrb_ad_tw
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Ptrb_ad, PT0, dV_ad, d_00
            call calc_Ptrb_ad_tw(s, k, Ptrb_ad, PT0, ierr)
            if (ierr /= 0) return
            d_00 = wrap_d_00(s,k)
            dV_ad = 1d0/d_00 - 1d0/s% rho_start(k)
            Ptrb_dV_ad = Ptrb_ad*dV_ad  ! erg cm^-3 cm^-3 g^-1 = erg g^-1
         end subroutine setup_Ptrb_dV_ad

         subroutine setup_dt_dLt_dm_ad(ierr)  ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Lt_00, Lt_p1
            real(dp) :: L_theta
            include 'formats'
            ierr = 0
            if (s% using_velocity_time_centering .and. &
                     s% include_L_in_velocity_time_centering) then
               L_theta = s% L_theta_for_velocity_time_centering
            else
               L_theta = 1d0
            end if
            Lt_00 = L_theta*s% Lt_ad(k) + (1d0 - L_theta)*s% Lt_start(k)
            if (k == s% nz) then
               Lt_p1 = 0d0
            else
               Lt_p1 = L_theta*shift_p1(s% Lt_ad(k+1)) + (1d0 - L_theta)*s% Lt_start(k+1)
               if (ierr /= 0) return
            end if
            dt_dLt_dm_ad = (Lt_00 - Lt_p1)*s%dt/s%dm(k)
         end subroutine setup_dt_dLt_dm_ad

         subroutine setup_dt_C_ad(ierr)  ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: C
            C = s% COUPL_ad(k)  ! compute_C(s, k, ierr) ! erg g^-1 s^-1
            if (ierr /= 0) return
            dt_C_ad = s%dt*C
         end subroutine setup_dt_C_ad

         subroutine setup_dt_Eq_ad(ierr)  ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Eq_cell
            Eq_cell = s% Eq_ad(k)  ! compute_Eq_cell(s, k, ierr) ! erg g^-1 s^-1
            if (ierr /= 0) return
            dt_Eq_ad = s%dt*Eq_cell
         end subroutine setup_dt_Eq_ad

      end subroutine do1_turbulent_energy_eqn


      subroutine get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: alfa, beta
         ! face_value = alfa*cell_value(k) + beta*cell_value(k-1)
         if (k == 1) call mesa_error(__FILE__,__LINE__,'bad k==1 for get_RSP2_alfa_beta_face_weights')
         if (s% RSP2_use_mass_interp_face_values) then
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1d0 - alfa
         else
            alfa = 0.5d0
            beta = 0.5d0
         end if
      end subroutine get_RSP2_alfa_beta_face_weights


      function compute_Y_face(s, k, ierr) result(Y_face)  ! superadiabatic gradient [unitless]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Y_face
         type(auto_diff_real_star_order1) :: gradT, grad_ad_face, gradL
         include 'formats'
         ierr = 0

         if (k > s% nz) then
            Y_face = 0d0
            return
         end if

         if (k == 1 .or. s% mixing_length_alpha == 0d0) then
            Y_face = 0d0
            s% Y_face(k) = 0d0
            s% Y_face_ad(k) = 0d0
            s% gradT_sub_grada(k) = 0d0
            return
         end if

         gradT = s% gradT_ad(k)
         grad_ad_face = s% grada_face_ad(k)
         gradL = s% gradL_ad(k)
         Y_face = gradT - gradL

         s% gradT_sub_grada(k) = gradT%val - grad_ad_face%val
         s% Y_face_ad(k) = Y_face
         s% Y_face(k) = Y_face%val

      end function compute_Y_face


      function compute_PII_face(s, k, ierr) result(PII_face)  ! ergs g^-1 K^-1 (like Cp)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: PII_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Y_face
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            PII_face = 0d0
            return
         end if
         if (k == 1 .or. s% mixing_length_alpha == 0d0 .or. &
               k == s% nz) then  ! just skip k == nz to be like RSP
            PII_face = 0d0
            s% PII(k) = 0d0
            s% PII_ad(k) = 0d0
            return
         end if
         Y_face = s% Y_face_ad(k)  ! compute_Y_face(s, k, ierr)
         if (ierr /= 0) return
         PII_face = compute_PII_from_Y(s, k, Y_face, ierr)
         if (ierr /= 0) return

         s% PII(k) = PII_face%val
         s% PII_ad(k) = PII_face
         if (k == -2 .and. s% PII(k) < 0d0) then
            write(*,2) 's% PII(k)', k, s% PII(k)
            write(*,2) 'Y_face', k, Y_face%val
            call mesa_error(__FILE__,__LINE__,'compute_PII_face')
         end if
      end function compute_PII_face


      function compute_PII_from_Y(s, k, Y_face, ierr) result(PII_face)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: Y_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: PII_face
         type(auto_diff_real_star_order1) :: &
            Cp_00, Cp_m1, Cp_face, T_00, T_m1, X, FL, scale, &
            T_face, e_face, Peos_face, rho_face, h_face
         real(dp) :: ALFAS_ALFA, alfa, beta

         ierr = 0
         Cp_00 = wrap_Cp_00(s, k)
         Cp_m1 = wrap_Cp_m1(s, k)
         T_00= wrap_T_00(s, k)
         T_m1 = wrap_T_m1(s, k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         Cp_face = alfa*Cp_00 + beta*Cp_m1  ! ergs g^-1 K^-1
         T_face = alfa*T_00 + beta*T_m1
         rho_face = alfa*wrap_d_00(s,k) + beta*wrap_d_m1(s,k)
         Peos_face = alfa*wrap_Peos_00(s,k) + beta*wrap_Peos_m1(s,k)
         e_face = alfa*wrap_e_00(s,k) + beta*wrap_e_m1(s,k)
         h_face = e_face + Peos_face/rho_face
         ALFAS_ALFA = x_ALFAS*s% mixing_length_alpha
         PII_face = ALFAS_ALFA*Cp_face*Y_face

         scale = 1d0
         if (Y_face > 0d0 .and. s% use_TDC_enthalpy_flux_limiter) then
            ! X = G/F
            X = (Cp_face*T_face/h_face)*ALFAS_ALFA* Y_face / sqrt_2_div_3
            FL = flux_limiter_function(X)
            ! Avoid 0/0 or tiny/tiny; for X near 0, FL near X so scale ~ 1 anyway.
            if (abs(X%val) >= 0.95d0) then
               scale = FL / X
            else
               scale = 1d0
            end if
         end if

         PII_face = PII_face*scale
      end function compute_PII_from_Y

      type(auto_diff_real_star_order1) function flux_limiter_function(X) result(FL) ! should be c2 continuous
        type(auto_diff_real_star_order1), intent(in) :: X
        real(dp), parameter :: X0 = 0.95d0, X1 = 1d0

        type(auto_diff_real_star_order1) :: s, p

        ! Region 1: purely linear, FL = X
        if (X%val < X0) then ! should not be encountered
           FL = X

        ! Region 3: saturated, FL = 1
        else if (X%val >= X1) then
           FL = 1.0_dp

        ! Region 2: smooth C2 transition between the two
        else
           ! Normalized coordinate in [0,1]
           s = (X - X0) / (X1 - X0)

           ! Quintic "smootherstep" polynomial:
           ! p(s) = 10 s^3 - 15 s^4 + 6 s^5
           ! p(0)=0, p(1)=1, p'(0)=p'(1)=0, p''(0)=p''(1)=0
           p = pow3(s) * (10.0_dp + s * (-15.0_dp + 6.0_dp * s))

           ! Blend between line FL=X and flat FL=1
           ! At s=0:  FL = X
           ! At s=1:  FL = 1
           ! Because p', p'' vanish at 0 and 1, FL, FL', FL'' all match.
           FL = X + (1.0_dp - X) * p
        end if
      end function flux_limiter_function


      type(auto_diff_real_star_order1) function flux_limiter_derivative(X) result(dFL_dX)
         type(auto_diff_real_star_order1), intent(in) :: X
         real(dp), parameter :: X0 = 0.95d0, X1 = 1d0
         type(auto_diff_real_star_order1) :: s, p, dp_dX

         if (X%val < X0) then
            dFL_dX = 1d0
         else if (X%val >= X1) then
            dFL_dX = 0d0
         else
            s = (X - X0)/(X1 - X0)
            p = pow3(s)*(10d0 + s*(-15d0 + 6d0*s))
            dp_dX = 30d0*pow2(s)*pow2(1d0 - s)/(X1 - X0)
            dFL_dX = 1d0 - p + (1d0 - X)*dp_dX
         end if
      end function flux_limiter_derivative

      function compute_d_v_div_r(s, k, ierr) result(d_v_div_r)  ! s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: d_v_div_r
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
         include 'formats'
         ierr = 0
         v_00 = wrap_v_00(s,k)
         v_p1 = wrap_v_p1(s,k)
         r_00 = wrap_r_00(s,k)
         r_p1 = wrap_r_p1(s,k)
         if (r_p1%val == 0d0) r_p1 = 1d0
         d_v_div_r = v_00/r_00 - v_p1/r_p1  ! units s^-1
      end function compute_d_v_div_r


      function compute_d_v_div_r_opt_time_center(s, k, ierr) result(d_v_div_r)  ! s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: d_v_div_r
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
         include 'formats'
         ierr = 0
         v_00 = wrap_opt_time_center_v_00(s,k)
         v_p1 = wrap_opt_time_center_v_p1(s,k)
         r_00 = wrap_opt_time_center_r_00(s,k)
         r_p1 = wrap_opt_time_center_r_p1(s,k)
         if (r_p1%val == 0d0) r_p1 = 1d0
         d_v_div_r = v_00/r_00 - v_p1/r_p1  ! units s^-1
      end function compute_d_v_div_r_opt_time_center


      function wrap_Hp_cell(s, k) result(Hp_cell)  ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Hp_cell, grav_cell

         grav_cell = 0.5d0*s% cgrav(k)*s% m_grav(k)/pow2(wrap_r_00(s,k))
         if (k < s% nz) &
            grav_cell = grav_cell + &
               0.5d0*s% cgrav(k+1)*s% m_grav(k+1)/pow2(wrap_r_p1(s,k))
         Hp_cell = wrap_Peos_00(s,k)/(wrap_d_00(s,k)*grav_cell)
      end function wrap_Hp_cell


      function Hp_cell_for_Chi(s, k, ierr) result(Hp_cell)  ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_cell
         include 'formats'
         ierr = 0

         Hp_cell = wrap_Hp_cell(s, k)
      end function Hp_cell_for_Chi


      function compute_Chi_cell(s, k, ierr) result(Chi_cell)
         ! eddy viscosity energy (Kuhfuss 1986) [erg]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Chi_cell
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_div_w, w_00
         ierr = 0
         Chi_div_w = compute_Chi_div_w_cell(s, k, ierr)
         if (ierr /= 0) return
         w_00 = wrap_w_00(s,k)
         Chi_cell = w_00*Chi_div_w
         s% Chi(k) = Chi_cell%val
         s% Chi_ad(k) = Chi_cell

      end function compute_Chi_cell


      function compute_Chi_div_w_cell(s, k, ierr) result(Chi_div_w)
         ! eddy viscosity energy divided by turbulent velocity [erg s cm^-1]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Chi_div_w
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            rho2, r6_cell, d_v_div_r, Hp_cell, d_00, r_00, r_p1
         real(dp) :: f, ALFAM_ALFA
         ierr = 0
         ALFAM_ALFA = s% RSP2_alfam*s% mixing_length_alpha
         if (ALFAM_ALFA == 0d0 .or. &
               k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
               k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            Chi_div_w = 0d0
            return
         end if
         Hp_cell = Hp_cell_for_Chi(s, k, ierr)
         if (ierr /= 0) return
         d_v_div_r = compute_d_v_div_r(s, k, ierr)
         if (ierr /= 0) return
         d_00 = wrap_d_00(s,k)
         f = (16d0/3d0)*pi*ALFAM_ALFA/s% dm(k)
         rho2 = pow2(d_00)
         r_00 = wrap_r_00(s,k)
         r_p1 = wrap_r_p1(s,k)
         r6_cell = 0.5d0*(pow6(r_00) + pow6(r_p1))
         Chi_div_w = f*rho2*r6_cell*d_v_div_r*Hp_cell
      end function compute_Chi_div_w_cell


      function compute_Eq_cell(s, k, ierr) result(Eq_cell)  ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Eq_cell
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Eq_div_w, w_00
         ierr = 0
         Eq_div_w = compute_Eq_div_w_cell(s, k, ierr)
         if (ierr /= 0) return
         w_00 = wrap_w_00(s,k)
         Eq_cell = w_00*Eq_div_w
         s% Eq(k) = Eq_cell%val
         s% Eq_ad(k) = Eq_cell
      end function compute_Eq_cell


      function compute_Eq_div_w_cell(s, k, ierr) result(Eq_div_w)
         ! eddy viscosity work divided by turbulent velocity [cm s^-2]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Eq_div_w
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: d_v_div_r, Chi_div_w
         ierr = 0
         Chi_div_w = compute_Chi_div_w_cell(s, k, ierr)
         if (ierr /= 0) return
         d_v_div_r = compute_d_v_div_r_opt_time_center(s, k, ierr)
         if (ierr /= 0) return
         Eq_div_w = 4d0*pi*Chi_div_w*d_v_div_r/s% dm(k)
      end function compute_Eq_div_w_cell


      function compute_Uq_face(s, k, ierr) result(Uq_face)  ! cm s^-2, acceleration
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Uq_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Chi_00, Chi_m1, r_00
         include 'formats'
         ierr = 0
         if (s% mixing_length_alpha == 0d0 .or. &
             k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            Uq_face = 0d0
         else
            r_00 = wrap_opt_time_center_r_00(s,k)
            Chi_00 = s% Chi_ad(k)  ! compute_Chi_cell(s,k,ierr)
            if (k > 1) then
               !Chi_m1 = shift_m1(compute_Chi_cell(s,k-1,ierr))
               Chi_m1 = shift_m1(s% Chi_ad(k-1))
               if (ierr /= 0) return
            else
               Chi_m1 = 0d0
            end if
            Uq_face = 4d0*pi*(Chi_m1 - Chi_00)/(r_00*s% dm_bar(k))

            if (k==-56) then
               write(*,3) 'RSP2 Uq chi_m1 chi_00 r', k, s% solver_iter, &
                  Uq_face%val, Chi_m1%val, Chi_00%val, r_00%val
            end if

         end if
         ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2, acceleration
         s% Uq(k) = Uq_face%val
      end function compute_Uq_face


      function compute_Source(s, k, ierr) result(Source)  ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Source
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: w_00, Source_div_w
         include 'formats'
         ierr = 0
         w_00 = wrap_w_00(s, k)
         Source_div_w = compute_Source_div_w(s, k, ierr)
         if (ierr /= 0) return
         Source = (w_00 + s% RSP2_source_seed)*Source_div_w

         if (k==-109) then
            write(*,3) 'RSP2 Source w source_div_w', k, s% solver_iter, &
               Source%val, w_00%val, Source_div_w%val
         end if
         s% SOURCE(k) = Source%val

      end function compute_Source


      function compute_Source_div_w(s, k, ierr) result(Source_div_w)  ! cm s^-2
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Source_div_w
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            T_00, d_00, Peos_00, Cp_00, chiT_00, chiRho_00, QQ_00, &
            Hp_face_00, Hp_face_p1, PII_face_00, PII_face_p1, PII_div_Hp_cell, &
            P_QQ_div_Cp
         ierr = 0
         T_00 = wrap_T_00(s, k)
         d_00 = wrap_d_00(s, k)
         Peos_00 = wrap_Peos_00(s, k)
         Cp_00 = wrap_Cp_00(s, k)
         chiT_00 = wrap_chiT_00(s, k)
         chiRho_00 = wrap_chiRho_00(s, k)
         QQ_00 = chiT_00/(d_00*T_00*chiRho_00)

         Hp_face_00 = wrap_Hp_00(s,k)
         PII_face_00 = s% PII_ad(k)  ! compute_PII_face(s, k, ierr)
         if (ierr /= 0) return

         if (k == s% nz) then
            PII_div_Hp_cell = PII_face_00/Hp_face_00
         else
            Hp_face_p1 = wrap_Hp_p1(s,k)
            if (ierr /= 0) return
            !PII_face_p1 = shift_p1(compute_PII_face(s, k+1, ierr))
            PII_face_p1 = shift_p1(s% PII_ad(k+1))
            if (ierr /= 0) return
            PII_div_Hp_cell = 0.5d0*(PII_face_00/Hp_face_00 + PII_face_p1/Hp_face_p1)
         end if

         ! Peos_00*QQ_00/Cp_00 = grad_ad if all perfect.
         !grad_ad_00 = wrap_grad_ad_00(s, k)
         P_QQ_div_Cp = Peos_00*QQ_00/Cp_00  ! use this to be same as RSP
         Source_div_w = PII_div_Hp_cell*T_00*P_QQ_div_Cp
      end function compute_Source_div_w


      function compute_D(s, k, ierr) result(D)  ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: D
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: D_div_w, w_00
         ierr = 0
         w_00 = wrap_w_00(s,k)
         D_div_w = compute_D_div_w(s, k, ierr)
         if (ierr /= 0) return
         D = w_00*D_div_w
         s% DAMP(k) = D%val
      end function compute_D


      function compute_D_div_w(s, k, ierr) result(D_div_w)  ! cm s^-2
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: D_div_w
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_cell, w_00
         ierr = 0
         if (s% mixing_length_alpha == 0d0) then
            D_div_w = 0d0
         else
            Hp_cell = wrap_Hp_cell(s,k)
            w_00 = wrap_w_00(s,k)
            D_div_w = (s% RSP2_alfad*x_CEDE/s% mixing_length_alpha)* &
               pow2(w_00)/Hp_cell
         end if
      end function compute_D_div_w


      function compute_Dr(s, k, ierr) result(Dr)  ! erg g^-1 s^-1 = cm^2 s^-3
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Dr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Dr_div_w, w_00
         ierr = 0
         w_00 = wrap_w_00(s,k)
         Dr_div_w = compute_Dr_div_w(s, k, ierr)
         if (ierr /= 0) return
         Dr = w_00*Dr_div_w
         s% DAMPR(k) = Dr%val
      end function compute_Dr


      function compute_Dr_div_w(s, k, ierr) result(Dr_div_w)  ! cm s^-2
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Dr_div_w
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_00, T_00, d_00, Cp_00, kap_00, Hp_cell, POM2
         real(dp) :: gammar, alpha, POM
         ierr = 0
         alpha = s% mixing_length_alpha
         gammar = s% RSP2_alfar*x_GAMMAR
         if (gammar == 0d0) then
            Dr_div_w = 0d0
            return
         end if
         w_00 = wrap_w_00(s,k)
         T_00 = wrap_T_00(s,k)
         d_00 = wrap_d_00(s,k)
         Cp_00 = wrap_Cp_00(s,k)
         kap_00 = wrap_kap_00(s,k)
         Hp_cell = wrap_Hp_cell(s,k)
         POM = 4d0*boltz_sigma*pow2(gammar/alpha)  ! erg cm^-2 K^-4 s^-1
         POM2 = pow3(T_00)/(pow2(d_00)*Cp_00*kap_00)
         Dr_div_w = w_00*POM*POM2/pow2(Hp_cell)
      end function compute_Dr_div_w


      function compute_C(s, k, ierr) result(C)  ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: C
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Source, D, Dr
         if (s% mixing_length_alpha == 0d0 .or. &
             k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            if (k >= 1 .and. k <= s% nz) then
               s% SOURCE(k) = 0d0
               s% DAMP(k) = 0d0
               s% DAMPR(k) = 0d0
               s% COUPL(k) = 0d0
               s% COUPL_ad(k) = 0d0
            end if
            C = 0d0
            return
         end if
         Source = compute_Source(s, k, ierr)
         if (ierr /= 0) return
         D = compute_D(s, k, ierr)
         if (ierr /= 0) return
         Dr = compute_Dr(s, k, ierr)
         if (ierr /= 0) return
         C = Source - D - Dr
         s% COUPL(k) = C%val
         s% COUPL_ad(k) = C
      end function compute_C


      function compute_L_face(s, k, ierr) result(L_face)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: L_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lr, Lc, Lt
         call compute_L_terms(s, k, L_face, Lr, Lc, Lt, ierr)
      end function compute_L_face


      subroutine compute_L_terms(s, k, L, Lr, Lc, Lt, ierr)
         type (star_info), pointer, intent(in) :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: L, Lr, Lc, Lt
         type(accurate_auto_diff_real_star_order1) :: L_sum
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            L = 0d0
            L%val = s% L_center
            Lr = 0d0
            Lc = 0d0
            Lt = 0d0
            return
         end if
         Lr = compute_Lr(s, k, ierr)
         if (ierr /= 0) return
         if (k == 1) then
            Lc = 0d0
            Lt = 0d0
         else
            Lc = compute_Lc(s, k, ierr)
            if (ierr /= 0) return
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
         end if
         L_sum = Lr
         L_sum = L_sum + Lc
         L_sum = L_sum + Lt
         L = L_sum
         s% Lr_ad(k) = Lr
         s% Lc_ad(k) = Lc
         s% Lt_ad(k) = Lt
      end subroutine compute_L_terms


      function compute_Lr(s, k, ierr) result(Lr)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            r_00, area, T_00, T400, Erad, Lrad_coeff, gradT
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            Lr = s% L_center
         else
            r_00 = wrap_r_00(s,k)  ! not time centered
            area = 4d0*pi*pow2(r_00)
            T_00 = wrap_T_00(s,k)
            T400 = pow4(T_00)
            if (k == 1) then
               if (s% RSP2_use_L_eqn_at_surface) then
                  Erad = crad*T400
                  Lr = s% RSP2_Lsurf_factor*area*clight*Erad
               else
                  Lr = wrap_L_00(s,k)
               end if
            else
               Lrad_coeff = compute_Lrad_coeff(s, k, ierr)
               if (ierr /= 0) return
               gradT = s% gradT_ad(k)
               Lr = Lrad_coeff*gradT
            end if
         end if
         s% Lr(k) = Lr%val
      end function compute_Lr


      function compute_Lrad_coeff(s, k, ierr) result(Lrad_coeff)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lrad_coeff
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            r_00, area, T_face, rho_face, kap_face, Hp_face, krad
         real(dp) :: alfa, beta
         include 'formats'
         ierr = 0

         r_00 = wrap_r_00(s,k)
         area = 4d0*pi*pow2(r_00)
         Hp_face = wrap_Hp_00(s,k)
         if (k == 1) then
            T_face = wrap_T_00(s,k)
            rho_face = wrap_d_00(s,k)
            kap_face = wrap_kap_00(s,k)
         else
            call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
            T_face = alfa*wrap_T_00(s,k) + beta*wrap_T_m1(s,k)
            rho_face = alfa*wrap_d_00(s,k) + beta*wrap_d_m1(s,k)
            kap_face = alfa*wrap_kap_00(s,k) + beta*wrap_kap_m1(s,k)
         end if
         krad = 4d0*crad*clight*pow3(T_face)/(3d0*kap_face*rho_face)
         Lrad_coeff = area*krad*T_face/Hp_face
      end function compute_Lrad_coeff


      function compute_Lc(s, k, ierr) result(Lc)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lc_div_w_face
         Lc = compute_Lc_terms(s, k, Lc_div_w_face, ierr)
         s% Lc(k) = Lc%val
      end function compute_Lc


      function compute_Lc_terms(s, k, Lc_div_w_face, ierr) result(Lc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc, Lc_div_w_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: r_00, area, &
            T_m1, T_00, d_m1, d_00, w_m1, w_00, T_rho_face, PII_face, w_face
         real(dp) :: ALFAC, ALFAS, alfa, beta
         include 'formats'
         ierr = 0
         if (s% mixing_length_alpha == 0d0 .or. &
             k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            Lc = 0d0
            Lc_div_w_face = 1
            return
         end if
         r_00 = wrap_r_00(s, k)
         area = 4d0*pi*pow2(r_00)
         T_m1 = wrap_T_m1(s, k)
         T_00 = wrap_T_00(s, k)
         d_m1 = wrap_d_m1(s, k)
         d_00 = wrap_d_00(s, k)
         w_m1 = wrap_w_m1(s, k)
         w_00 = wrap_w_00(s, k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         T_rho_face = alfa*T_00*d_00 + beta*T_m1*d_m1
         PII_face = s% PII_ad(k)  ! compute_PII_face(s, k, ierr)
         w_face = alfa*w_00 + beta*w_m1
         ALFAC = x_ALFAC
         ALFAS = x_ALFAS
         Lc_div_w_face = area*(ALFAC/ALFAS)*T_rho_face*PII_face
         ! units = cm^2 K g cm^-3 ergs g^-1 K^-1 = ergs cm^-1
         Lc = w_face*Lc_div_w_face
         ! units = cm s^-1 ergs cm^-1 = ergs s^-1
         if (k == -458) then
            write(*,2) 'Lc%val', k, Lc%val
            write(*,2) 'w_face%val', k, w_face%val
            write(*,2) 'Lc_div_w_face', k, Lc_div_w_face%val
            write(*,2) 'PII_face%val', k, PII_face%val
            write(*,2) 'T_rho_face%val', k, T_rho_face%val
            !write(*,2) '', k,
            !write(*,2) '', k,
            call mesa_error(__FILE__,__LINE__,'compute_Lc_terms')
         end if
      end function compute_Lc_terms


      function compute_Lc_gradT_coeff( &
            s, k, limiter_argument_coeff, ierr, w_m1_in, w_00_in) &
            result(Lc_gradT_coeff)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc_gradT_coeff
         type(auto_diff_real_star_order1), intent(out) :: limiter_argument_coeff
         integer, intent(out) :: ierr
         real(dp), intent(in), optional :: w_m1_in, w_00_in
         type(auto_diff_real_star_order1) :: &
            r_00, area, T_m1, T_00, d_m1, d_00, w_m1, w_00, &
            Cp_m1, Cp_00, T_rho_face, Cp_face, w_face, T_face, &
            rho_face, Peos_face, e_face, h_face
         real(dp) :: alfa, beta
         include 'formats'
         ierr = 0
         limiter_argument_coeff = 0d0
         if (k == 1 .or. k == s% nz .or. s% mixing_length_alpha == 0d0 .or. &
             k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            Lc_gradT_coeff = 0d0
            return
         end if
         r_00 = wrap_r_00(s, k)
         area = 4d0*pi*pow2(r_00)
         T_m1 = wrap_T_m1(s, k)
         T_00 = wrap_T_00(s, k)
         d_m1 = wrap_d_m1(s, k)
         d_00 = wrap_d_00(s, k)
         if (present(w_m1_in) .and. present(w_00_in)) then
            w_m1 = w_m1_in
            w_00 = w_00_in
         else
            w_m1 = wrap_w_m1(s, k)
            w_00 = wrap_w_00(s, k)
         end if
         Cp_m1 = wrap_Cp_m1(s, k)
         Cp_00 = wrap_Cp_00(s, k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         T_rho_face = alfa*T_00*d_00 + beta*T_m1*d_m1
         Cp_face = alfa*Cp_00 + beta*Cp_m1
         w_face = alfa*w_00 + beta*w_m1
         T_face = alfa*T_00 + beta*T_m1
         rho_face = alfa*d_00 + beta*d_m1
         Peos_face = alfa*wrap_Peos_00(s,k) + beta*wrap_Peos_m1(s,k)
         e_face = alfa*wrap_e_00(s,k) + beta*wrap_e_m1(s,k)
         h_face = e_face + Peos_face/rho_face
         Lc_gradT_coeff = &
            area*x_ALFAC*s% mixing_length_alpha*T_rho_face*Cp_face*w_face
         limiter_argument_coeff = &
            x_ALFAS*s% mixing_length_alpha*Cp_face*T_face/(sqrt_2_div_3*h_face)
         ! Lc = Lc_gradT_coeff*(gradT - gradL) without the flux limiter.
      end function compute_Lc_gradT_coeff


      function compute_RSP2_gradT(s, k, ierr, w_m1_in, w_00_in) result(gradT)
         use turb_support, only: get_TDC_dynamical_gradL
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: gradT
         integer, intent(out) :: ierr
         real(dp), intent(in), optional :: w_m1_in, w_00_in
         type(auto_diff_real_star_order1) :: &
            Lrad_coeff, Lc_gradT_coeff, thermal_conductance, Lt, &
            grad_ad, gradL, dynamical_gradL, limiter_argument_coeff, &
            limiter_argument, flux_constant, flux_resid, d_flux_resid_d_argument
         real(dp) :: alfa, beta, argument, argument_lo, argument_hi, &
            argument_new, flux_scale
         integer :: iter
         include 'formats'
         ierr = 0
         if (k == 1) then
            gradT = s% gradT_ad(k)
            return
         end if
         Lrad_coeff = compute_Lrad_coeff(s, k, ierr)
         if (ierr /= 0) return
         if (present(w_m1_in) .and. present(w_00_in)) then
            Lc_gradT_coeff = compute_Lc_gradT_coeff( &
               s, k, limiter_argument_coeff, ierr, w_m1_in, w_00_in)
         else
            Lc_gradT_coeff = compute_Lc_gradT_coeff( &
               s, k, limiter_argument_coeff, ierr)
         end if
         if (ierr /= 0) return
         thermal_conductance = Lrad_coeff + Lc_gradT_coeff
         if (thermal_conductance%val <= 0d0 .or. is_bad(thermal_conductance%val)) then
            ierr = -1
            if (s% report_ierr) &
               write(*,2) 'bad RSP2 thermal conductance', k, thermal_conductance%val
            return
         end if
         if (present(w_m1_in) .and. present(w_00_in)) then
            Lt = compute_Lt(s, k, ierr, w_m1_in, w_00_in)
         else
            Lt = compute_Lt(s, k, ierr)
         end if
         if (ierr /= 0) return
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         grad_ad = alfa*wrap_grad_ad_00(s,k) + beta*wrap_grad_ad_m1(s,k)
         gradL = grad_ad
         if (s% use_Ledoux_criterion) &
            gradL = gradL + s% gradL_composition_term(k)
         if (s% TDC_use_dynamical_gradL) then
            call get_TDC_dynamical_gradL(s, k, gradL, dynamical_gradL, ierr)
            if (ierr /= 0) return
            gradL = dynamical_gradL
         end if
         s% grada_face_ad(k) = grad_ad
         s% grada_face(k) = grad_ad%val
         s% gradL_ad(k) = gradL
         s% gradL(k) = gradL%val

         ! L - Lt = Kr*gradT + Kc*(gradT - gradL).
         gradT = &
            (wrap_L_00(s,k) - Lt + Lc_gradT_coeff*gradL)/thermal_conductance
         if (.not. s% use_TDC_enthalpy_flux_limiter .or. &
               Lc_gradT_coeff%val <= 0d0 .or. gradT%val <= gradL%val) return

         if (limiter_argument_coeff%val <= 0d0 .or. &
               is_bad(limiter_argument_coeff%val)) then
            ierr = -1
            if (s% report_ierr) write(*,2) &
               'bad RSP2 flux limiter coefficient', k, limiter_argument_coeff%val
            return
         end if

         limiter_argument = limiter_argument_coeff*(gradT - gradL)
         if (limiter_argument%val <= 0.95d0) return

         ! Kr*X + Kc*FL(X) + b*(Kr*gradL + Lt - L) = 0.
         flux_constant = limiter_argument_coeff*( &
            Lrad_coeff*gradL + Lt - wrap_L_00(s,k))
         limiter_argument = -(Lc_gradT_coeff + flux_constant)/Lrad_coeff
         if (is_bad(limiter_argument%val)) then
            ierr = -1
            if (s% report_ierr) write(*,2) &
               'bad saturated RSP2 flux limiter root', k, limiter_argument%val
            return
         end if
         if (limiter_argument%val >= 1d0) then
            gradT = gradL + limiter_argument/limiter_argument_coeff
            return
         end if

         ! The unlimited and saturated roots bracket the transition root.
         argument_lo = 0.95d0
         argument_hi = 1d0
         argument = 0.5d0*(argument_lo + argument_hi)
         flux_scale = max(abs(Lrad_coeff%val), abs(Lc_gradT_coeff%val), &
            abs(flux_constant%val), 1d0)
         do iter = 1, 50
            limiter_argument = argument
            flux_resid = Lrad_coeff*limiter_argument + &
               Lc_gradT_coeff*flux_limiter_function(limiter_argument) + flux_constant
            if (abs(flux_resid%val) <= 1d-12*flux_scale) exit
            if (flux_resid%val > 0d0) then
               argument_hi = argument
            else
               argument_lo = argument
            end if
            d_flux_resid_d_argument = Lrad_coeff + &
               Lc_gradT_coeff*flux_limiter_derivative(limiter_argument)
            if (d_flux_resid_d_argument%val <= 0d0 .or. &
                  is_bad(d_flux_resid_d_argument%val)) then
               ierr = -1
               if (s% report_ierr) write(*,2) &
                  'bad RSP2 flux limiter derivative', k, d_flux_resid_d_argument%val
               return
            end if
            argument_new = argument - &
               flux_resid%val/d_flux_resid_d_argument%val
            if (argument_new <= argument_lo .or. argument_new >= argument_hi) &
               argument_new = 0.5d0*(argument_lo + argument_hi)
            argument = argument_new
         end do

         limiter_argument = argument
         flux_resid = Lrad_coeff*limiter_argument + &
            Lc_gradT_coeff*flux_limiter_function(limiter_argument) + flux_constant
         if (abs(flux_resid%val) > 1d-10*flux_scale) then
            ierr = -1
            if (s% report_ierr) write(*,2) &
               'RSP2 flux limiter root failed', k, flux_resid%val/flux_scale
            return
         end if
         d_flux_resid_d_argument = Lrad_coeff + &
            Lc_gradT_coeff*flux_limiter_derivative(limiter_argument)
         if (d_flux_resid_d_argument%val <= 0d0 .or. &
               is_bad(d_flux_resid_d_argument%val)) then
            ierr = -1
            if (s% report_ierr) write(*,2) &
               'bad RSP2 flux limiter derivative', k, d_flux_resid_d_argument%val
            return
         end if
         limiter_argument = limiter_argument - flux_resid/d_flux_resid_d_argument
         limiter_argument%val = argument
         gradT = gradL + limiter_argument/limiter_argument_coeff
      end function compute_RSP2_gradT


      function compute_Lt(s, k, ierr, w_m1_in, w_00_in) result(Lt)  ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lt
         integer, intent(out) :: ierr
         real(dp), intent(in), optional :: w_m1_in, w_00_in
         type(auto_diff_real_star_order1) :: r_00, area2, d_m1, d_00, &
            rho2_face, Hp_face, w_m1, w_00, w_face, etrb_m1, etrb_00
         real(dp) :: alpha_alpha_t, alfa, beta
         logical :: store_result
         include 'formats'
         ierr = 0
         store_result = .not. (present(w_m1_in) .and. present(w_00_in))
         if (k > s% nz) then
            Lt = 0d0
            return
         end if
         alpha_alpha_t = s% mixing_length_alpha*s% RSP2_alfat
         if (alpha_alpha_t == 0d0 .or. &
             k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - int(s% nz/s% RSP2_nz_div_IBOTOM)) then
            Lt = 0d0
            if (store_result) s% Lt(k) = 0d0
            return
         end if
         r_00 = wrap_r_00(s,k)
         area2 = pow2(4d0*pi*pow2(r_00))
         d_m1 = wrap_d_m1(s,k)
         d_00 = wrap_d_00(s,k)
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         rho2_face = alfa*pow2(d_00) + beta*pow2(d_m1)
         if (present(w_m1_in) .and. present(w_00_in)) then
            w_m1 = w_m1_in
            w_00 = w_00_in
         else
            w_m1 = wrap_w_m1(s,k)
            w_00 = wrap_w_00(s,k)
         end if
         w_face = alfa*w_00 + beta*w_m1
         etrb_m1 = pow2(w_m1)
         etrb_00 = pow2(w_00)
         Hp_face = wrap_Hp_00(s,k)
         ! Ft = - alpha_t * rho_face * alpha * Hp_face * w_face * detrb/dr (thesis eqn 2.44)
         ! replace dr by dm_bar/(area*rho_face)
         ! Ft = - alpha_alpha_t * rho_face * Hp_face * w_face * (area*rho_face) * detrb/dm_bar
         ! Lt = area * Ft
         ! Lt = -alpha_alpha_t * (area*rho_face)**2 * Hp_face * w_face * (etrb(k-1) - etrb(k))/dm_bar
         Lt = - alpha_alpha_t * area2 * rho2_face * Hp_face * w_face * (etrb_m1 - etrb_00) / s% dm_bar(k)
         ! units = (cm^4) (g^2 cm^-6) (cm) (cm s^-1) (ergs g^-1) g^-1 = erg s^-1
         if (store_result) s% Lt(k) = Lt%val
      end function compute_Lt


      subroutine set_etrb_start_vars(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k
         type(auto_diff_real_star_order1) :: Y_face, Lt
         include 'formats'
         ierr = 0
         do k=1,s%nz
            Y_face = compute_Y_face(s, k, ierr)
            if (ierr /= 0) return
            s% Y_face_start(k) = Y_face%val
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
            s% Lt_start(k) = Lt%val
            s% w_start(k) = s% w(k)
            s% Hp_face_start(k) = s% Hp_face(k)
         end do
      end subroutine set_etrb_start_vars


      subroutine RSP2_adjust_vars_before_call_solver(s, ierr)  ! replaces check_omega in RSP
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: PII_div_Hp, QQ, source_coeff, grav_cell, Hp_cell, &
            damping_coeff, rad_damping_coeff, linear_coeff, discr, soln, w_initial
         integer :: k, k_lo, k_hi
         include 'formats'
         ierr = 0
         if (s% mixing_length_alpha == 0d0 .or. s% dt <= 0d0 .or. &
               s% RSP2_source_seed /= 0d0) return

         k_lo = s% RSP2_num_outermost_cells_forced_nonturbulent + 1
         k_hi = s% nz - max(1,int(s% nz/s% RSP2_nz_div_IBOTOM))
         if (k_lo > k_hi) return

         !$OMP PARALLEL DO PRIVATE(k,PII_div_Hp,QQ,source_coeff,grav_cell,Hp_cell, &
         !$OMP    damping_coeff,rad_damping_coeff,linear_coeff,discr,soln,w_initial) &
         !$OMP SCHEDULE(dynamic,2)
         do k = k_lo, k_hi
            ! Seed only an exactly dormant accepted cell.
            if (s% w_start(k) /= 0d0) cycle

            if (s% Hp_face(k) <= 0d0 .or. s% Hp_face(k+1) <= 0d0) cycle
            PII_div_Hp = 0.5d0*( &
               s% PII(k)/s% Hp_face(k) + s% PII(k+1)/s% Hp_face(k+1))
            QQ = s% chiT(k)/(s% rho(k)*s% T(k)*s% chiRho(k))
            source_coeff = PII_div_Hp*s% T(k)*s% Peos(k)*QQ/s% Cp(k)
            if (source_coeff <= 0d0 .or. is_bad(source_coeff)) cycle

            grav_cell = 0.5d0*( &
               s% cgrav(k)*s% m_grav(k)/pow2(s% r(k)) + &
               s% cgrav(k+1)*s% m_grav(k+1)/pow2(s% r(k+1)))
            if (grav_cell <= 0d0 .or. is_bad(grav_cell)) cycle
            Hp_cell = s% Peos(k)/(s% rho(k)*grav_cell)
            if (Hp_cell <= 0d0 .or. is_bad(Hp_cell)) cycle

            damping_coeff = &
               (s% RSP2_alfad*x_CEDE/s% mixing_length_alpha)/Hp_cell
            rad_damping_coeff = &
               4d0*boltz_sigma*pow2(s% RSP2_alfar*x_GAMMAR/ &
                  s% mixing_length_alpha)*pow3(s% T(k))/ &
               (pow2(s% rho(k))*s% Cp(k)*s% opacity(k)*pow2(Hp_cell))
            if (damping_coeff < 0d0 .or. rad_damping_coeff < 0d0 .or. &
                  is_bad(damping_coeff) .or. is_bad(rad_damping_coeff)) cycle

            linear_coeff = 1d0 + s% dt*rad_damping_coeff
            discr = pow2(linear_coeff) + &
               4d0*pow2(s% dt)*damping_coeff*source_coeff
            if (discr < 0d0 .or. is_bad(discr)) cycle
            soln = 2d0*s% dt*source_coeff/(linear_coeff + sqrt(discr))
            if (soln <= s% w(k) .or. is_bad(soln)) cycle

            w_initial = s% w(k)
            s% w(k) = soln
            if (s% RSP2_report_adjust_w) &
               write(*,3) 'RSP2_adjust_vars_before_call_solver w', &
                  k, s% model_number, w_initial, soln
         end do
         !$OMP END PARALLEL DO

      end subroutine RSP2_adjust_vars_before_call_solver

      end module hydro_rsp2
