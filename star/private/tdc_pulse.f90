! ***********************************************************************
!
!   Copyright (C) 2010-2025  Ebraheem Farag & The MESA Team
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

module tdc_pulse

   use star_private_def
   use const_def, only: dp, boltz_sigma, pi, clight, crad, ln10
   use utils_lib, only: is_bad
   use auto_diff
   use auto_diff_support
   use accurate_sum_auto_diff_star_order1
   use star_utils

   implicit none

   private
   public :: &
      compute_tdc_Eq_cell, compute_tdc_Uq_face, compute_tdc_Eq_face, &
      get_RSP2_alfa_beta_face_weights, set_viscosity_vars_TDC

contains

   ! This routine is called to initialize eq and uq for TDC.
   subroutine set_viscosity_vars_TDC(s, ierr)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: x
      integer :: k, op_err
      include 'formats'
      ierr = 0
      op_err = 0

      !$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
      do k = 1, s%nz
         ! Hp_face(k) <= 0 means it needs to be set.  e.g., after read file
         if (s%Hp_face(k) <= 0) then
            ! this scale height for face is already calculated in TDC
            s%Hp_face(k) = get_scale_height_face_val(s, k) ! because this is called before s% scale_height(k) is updated in mlt_vars.
         end if
      end do
      !$OMP END PARALLEL DO
      if (ierr /= 0) then
         if (s%report_ierr) write (*, 2) 'failed in set_viscosity_vars_TDC loop 1', s%model_number
         return
      end if
      !$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
      do k = 1, s%nz
         x = compute_Chi_cell(s, k, op_err)
         if (op_err /= 0) ierr = op_err
         x = compute_tdc_Eq_cell(s, k, op_err)
         if (op_err /= 0) ierr = op_err
         x = compute_tdc_Uq_face(s, k, op_err)
         if (op_err /= 0) ierr = op_err
      end do
      !$OMP END PARALLEL DO
      if (ierr /= 0) then
         if (s%report_ierr) write (*, 2) 'failed in set_viscosity_vars_TDC loop 2', s%model_number
         return
      end if
      if (.not. (s%v_flag .or. s%u_flag)) then ! set values 0 if not using v_flag or u_flag.
         do k = 1, s%nz
            s%Eq(k) = 0d0; s%Eq_ad(k) = 0d0
            s%Chi(k) = 0d0; s%Chi_ad(k) = 0d0
            s%Uq(k) = 0d0
         end do
      end if
   end subroutine set_viscosity_vars_TDC

   subroutine do1_rsp2_L_eqn(s, k, nvar, ierr)
      use star_utils, only: save_eqn_residual_info
      type(star_info), pointer :: s
      integer, intent(in) :: k, nvar
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: &
         L_expected, L_actual, resid
      real(dp) :: scale, residual, L_start_max
      logical :: test_partials
      include 'formats'

      !test_partials = (k == s% solver_test_partials_k)
      test_partials = .false.
      if (.not. s%RSP2_flag) then
         ierr = -1
         return
      end if

      ierr = 0
      !L_expected = compute_L_face(s, k, ierr)
      !if (ierr /= 0) return
      L_expected = s%Lr_ad(k) + s%Lc_ad(k) + s%Lt_ad(k)
      L_actual = wrap_L_00(s, k)
      L_start_max = maxval(s%L_start(1:s%nz))
      scale = 1d0/L_start_max
      if (is_bad(scale)) then
         write (*, 2) 'do1_rsp2_L_eqn scale', k, scale
         call mesa_error(__FILE__, __LINE__, 'do1_rsp2_L_eqn')
      end if
      resid = (L_expected - L_actual)*scale

      residual = resid%val
      s%equ(s%i_equL, k) = residual
      if (test_partials) then
         s%solver_test_partials_val = residual
      end if

      call save_eqn_residual_info(s, k, nvar, s%i_equL, resid, 'do1_rsp2_L_eqn', ierr)
      if (ierr /= 0) return

      if (test_partials) then
         s%solver_test_partials_var = s%i_lnR
         s%solver_test_partials_dval_dx = resid%d1Array(i_lnR_00)
         write (*, 4) 'do1_rsp2_L_eqn', s%solver_test_partials_var
      end if
   end subroutine do1_rsp2_L_eqn

   subroutine do1_rsp2_Hp_eqn(s, k, nvar, ierr)
      use star_utils, only: save_eqn_residual_info
      type(star_info), pointer :: s
      integer, intent(in) :: k, nvar
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: &
         Hp_expected, Hp_actual, resid
      real(dp) :: residual, Hp_start
      logical :: test_partials
      include 'formats'
      !test_partials = (k == s% solver_test_partials_k)
      test_partials = .false.

      if (.not. s%RSP2_flag) then
         ierr = -1
         return
      end if

      ierr = 0
      Hp_expected = Hp_face_for_rsp2_eqn(s, k, ierr)
      if (ierr /= 0) return
      Hp_actual = wrap_Hp_00(s, k)
      Hp_start = s%Hp_face_start(k)
      resid = (Hp_expected - Hp_actual)/max(Hp_expected, Hp_actual)

      residual = resid%val
      s%equ(s%i_equ_Hp, k) = residual
      if (test_partials) then
         s%solver_test_partials_val = residual
      end if

      if (residual > 1d3) then
         !$omp critical (hydro_rsp2_1)
         write (*, 2) 'residual', k, residual
         write (*, 2) 'Hp_expected', k, Hp_expected%val
         write (*, 2) 'Hp_actual', k, Hp_actual%val
         call mesa_error(__FILE__, __LINE__, 'do1_rsp2_Hp_eqn')
         !$omp end critical (hydro_rsp2_1)
      end if

      call save_eqn_residual_info(s, k, nvar, s%i_equ_Hp, resid, 'do1_rsp2_Hp_eqn', ierr)
      if (ierr /= 0) return

      if (test_partials) then
         s%solver_test_partials_var = s%i_lnR
         s%solver_test_partials_dval_dx = resid%d1Array(i_lnR_00)
         write (*, 4) 'do1_rsp2_Hp_eqn', s%solver_test_partials_var
      end if

   end subroutine do1_rsp2_Hp_eqn

   real(dp) function Hp_face_for_rsp2_val(s, k, ierr) result(Hp_face)  ! cm
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Hp_face_ad
      ierr = 0
      Hp_face_ad = Hp_face_for_rsp2_eqn(s, k, ierr)
      if (ierr /= 0) return
      Hp_face = Hp_face_ad%val
   end function Hp_face_for_rsp2_val

   function Hp_face_for_rsp2_eqn(s, k, ierr) result(Hp_face)  ! cm
      type(star_info), pointer :: s
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
      if (k > s%nz) then
         Hp_face = 1d0  ! not used
         return
      end if
      if (k > 1 .and. .not. s%RSP2_assume_HSE) then
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         rho_face = alfa*wrap_d_00(s, k) + beta*wrap_d_m1(s, k)
         area = 4d0*pi*pow2(wrap_r_00(s, k))
         dlnPeos = wrap_lnPeos_m1(s, k) - wrap_lnPeos_00(s, k)
         Hp_face = -s%dm_bar(k)/(area*rho_face*dlnPeos)
      else
         r_00 = wrap_r_00(s, k)  ! not time-centered in RSP
         d_00 = wrap_d_00(s, k)
         Peos_00 = wrap_Peos_00(s, k)
         if (k == 1) then
            Peos_div_rho = Peos_00/d_00
            Hp_face = pow2(r_00)*Peos_div_rho/(s%cgrav(k)*s%m(k))
         else
            d_m1 = wrap_d_m1(s, k)
            Peos_m1 = wrap_Peos_m1(s, k)
            call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
            Peos_div_rho = alfa*Peos_00/d_00 + beta*Peos_m1/d_m1
            Hp_face = pow2(r_00)*Peos_div_rho/(s%cgrav(k)*s%m(k))
            if (k == -104) then
               write (*, 3) 'RSP2 Hp P_div_rho Pdrho_00 Pdrho_m1', k, s%solver_iter, &
                  Hp_face%val, Peos_div_rho%val, Peos_00%val/d_00%val, Peos_m1%val/d_m1%val
               !write(*,3) 'RSP2 Hp r2_div_Gm r_start r', k, s% solver_iter, &
               !   Hp_face%val, pow2(r_00%val)/(s% cgrav(k)*s% m(k)), &
               !   s% r_start(k), r_00%val
            end if
            if (s%alt_scale_height_flag) then
               call mesa_error(__FILE__, __LINE__, 'Hp_face_for_rsp2_eqn: cannot use alt_scale_height_flag')
               ! consider sound speed*hydro time scale as an alternative scale height
               d_face = alfa*d_00 + beta*d_m1
               Peos_face = alfa*Peos_00 + beta*Peos_m1
               alt_Hp_face = sqrt(Peos_face/s%cgrav(k))/d_face
               if (alt_Hp_face%val < Hp_face%val) then  ! blend
                  A = pow2(alt_Hp_face/Hp_face)  ! 0 <= A%val < 1
                  Hp_face = A*Hp_face + (1d0 - A)*alt_Hp_face
               end if
            end if
         end if
      end if
   end function Hp_face_for_rsp2_eqn

   subroutine do1_turbulent_energy_eqn(s, k, nvar, ierr)
      use star_utils, only: set_energy_eqn_scal, save_eqn_residual_info
      type(star_info), pointer :: s
      integer, intent(in) :: k, nvar
      integer, intent(out) :: ierr
      ! for OLD WAY
      type(auto_diff_real_star_order1) :: &
         d_turbulent_energy_ad, Ptrb_dV_ad, dt_C_ad, dt_Eq_ad
      type(auto_diff_real_star_order1) :: w_00
      type(auto_diff_real_star_order1) :: tst, resid_ad, dt_dLt_dm_ad
      type(accurate_auto_diff_real_star_order1) :: esum_ad
      logical :: non_turbulent_cell, test_partials
      real(dp) :: residual, scal
      include 'formats'
      !test_partials = (k == s% solver_test_partials_k)
      test_partials = .false.

      ierr = 0
      w_00 = wrap_w_00(s, k)

      non_turbulent_cell = &
         s%mixing_length_alpha == 0d0 .or. &
         k <= s%RSP2_num_outermost_cells_forced_nonturbulent .or. &
         k > s%nz - int(s%nz/s%RSP2_nz_div_IBOTOM)
      if (.not. s%RSP2_flag) then
         resid_ad = w_00 - s%w_start(k)  ! just hold w constant when not using RSP2
      else if (non_turbulent_cell) then
         resid_ad = w_00/s%csound(k)  ! make w = 0
      else
         call setup_d_turbulent_energy(ierr); if (ierr /= 0) return  ! erg g^-1 = cm^2 s^-2
         call setup_Ptrb_dV_ad(ierr); if (ierr /= 0) return  ! erg g^-1
         call setup_dt_dLt_dm_ad(ierr); if (ierr /= 0) return  ! erg g^-1
         call setup_dt_C_ad(ierr); if (ierr /= 0) return  ! erg g^-1
         call setup_dt_Eq_ad(ierr); if (ierr /= 0) return  ! erg g^-1
         call set_energy_eqn_scal(s, k, scal, ierr); if (ierr /= 0) return  ! 1/(erg g^-1 s^-1)
         ! sum terms in esum_ad using accurate_auto_diff_real_star_order1
         esum_ad = d_turbulent_energy_ad + Ptrb_dV_ad + dt_dLt_dm_ad - dt_C_ad - dt_Eq_ad  ! erg g^-1
         resid_ad = esum_ad

         if (k == -35 .and. s%solver_iter == 1) then
            write (*, 3) 'RSP2 w dEt PdV dtC dtEq', k, s%solver_iter, &
               w_00%val, d_turbulent_energy_ad%val, Ptrb_dV_ad%val, dt_C_ad%val, dt_Eq_ad%val
         end if

         resid_ad = resid_ad*scal/s%dt  ! to make residual unitless, must cancel out the dt in scal

      end if

      residual = resid_ad%val
      s%equ(s%i_detrb_dt, k) = residual

      if (test_partials) then
         tst = residual
         s%solver_test_partials_val = tst%val
         if (s%solver_iter == 12) &
            write (*, *) 'do1_turbulent_energy_eqn', s%solver_test_partials_var, s%lnd(k), tst%val
      end if

      call save_eqn_residual_info(s, k, nvar, s%i_detrb_dt, resid_ad, 'do1_turbulent_energy_eqn', ierr)
      if (ierr /= 0) return

      if (test_partials) then
         s%solver_test_partials_var = s%i_lnd
         s%solver_test_partials_dval_dx = tst%d1Array(i_lnd_00)     ! xi0 good , xi1 partial 0, xi2 good.  Af horrible.'
         write (*, *) 'do1_turbulent_energy_eqn', s%solver_test_partials_var, s%lnd(k)/ln10, tst%val
      end if

   contains

      subroutine setup_d_turbulent_energy(ierr)  ! erg g^-1
         integer, intent(out) :: ierr
         ierr = 0
         d_turbulent_energy_ad = wrap_etrb_00(s, k) - get_etrb_start(s, k)
      end subroutine setup_d_turbulent_energy

      ! Ptrb_dV_ad = Ptrb_ad*dV_ad
      subroutine setup_Ptrb_dV_ad(ierr)  ! erg g^-1
         use star_utils, only: calc_Ptrb_ad_tw
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Ptrb_ad, PT0, dV_ad, d_00
         call calc_Ptrb_ad_tw(s, k, Ptrb_ad, PT0, ierr)
         if (ierr /= 0) return
         d_00 = wrap_d_00(s, k)
         dV_ad = 1d0/d_00 - 1d0/s%rho_start(k)
         Ptrb_dV_ad = Ptrb_ad*dV_ad  ! erg cm^-3 cm^-3 g^-1 = erg g^-1
      end subroutine setup_Ptrb_dV_ad

      subroutine setup_dt_dLt_dm_ad(ierr)  ! erg g^-1
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lt_00, Lt_p1
         real(dp) :: L_theta
         include 'formats'
         ierr = 0
         if (s%using_velocity_time_centering .and. &
             s%include_L_in_velocity_time_centering) then
            L_theta = s%L_theta_for_velocity_time_centering
         else
            L_theta = 1d0
         end if
         Lt_00 = L_theta*s%Lt_ad(k) + (1d0 - L_theta)*s%Lt_start(k)
         if (k == s%nz) then
            Lt_p1 = 0d0
         else
            Lt_p1 = L_theta*shift_p1(s%Lt_ad(k + 1)) + (1d0 - L_theta)*s%Lt_start(k + 1)
            if (ierr /= 0) return
         end if
         dt_dLt_dm_ad = (Lt_00 - Lt_p1)*s%dt/s%dm(k)
      end subroutine setup_dt_dLt_dm_ad

      subroutine setup_dt_C_ad(ierr)  ! erg g^-1
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: C
         C = s%COUPL_ad(k)  ! compute_C(s, k, ierr) ! erg g^-1 s^-1
         if (ierr /= 0) return
         dt_C_ad = s%dt*C
      end subroutine setup_dt_C_ad

      subroutine setup_dt_Eq_ad(ierr)  ! erg g^-1
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Eq_cell
         Eq_cell = s%Eq_ad(k)  ! compute_tdc_Eq_cell(s, k, ierr) ! erg g^-1 s^-1
         if (ierr /= 0) return
         dt_Eq_ad = s%dt*Eq_cell
      end subroutine setup_dt_Eq_ad

   end subroutine do1_turbulent_energy_eqn

   subroutine get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(out) :: alfa, beta
      ! face_value = alfa*cell_value(k) + beta*cell_value(k-1)
      if (k == 1) call mesa_error(__FILE__, __LINE__, 'bad k==1 for get_RSP2_alfa_beta_face_weights')
      if (s%RSP2_use_mass_interp_face_values) then
         alfa = s%dq(k - 1)/(s%dq(k - 1) + s%dq(k))
         beta = 1d0 - alfa
      else
         alfa = 0.5d0
         beta = 0.5d0
      end if
   end subroutine get_RSP2_alfa_beta_face_weights

   function compute_d_v_div_r(s, k, ierr) result(d_v_div_r)  ! s^-1
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: d_v_div_r
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1, term1, term2
      logical :: dbg
      include 'formats'
      ierr = 0
      dbg = .false.
      v_00 = wrap_v_00(s, k)
      v_p1 = wrap_v_p1(s, k)
      r_00 = wrap_r_00(s, k)
      r_p1 = wrap_r_p1(s, k)
      if (r_p1%val == 0d0) r_p1 = 1d0
      d_v_div_r = v_00/r_00 - v_p1/r_p1 ! units s^-1

      ! Debugging output to trace values
      if (dbg .and. k == -63) then
         write (*, *) 'test d_v_div_r, k:', k
         write (*, *) 'v_00:', v_00%val, 'v_p1:', v_p1%val
         write (*, *) 'r_00:', r_00%val, 'r_p1:', r_p1%val
         write (*, *) 'd_v_div_r:', d_v_div_r%val
      end if
   end function compute_d_v_div_r

   function compute_rho_form_of_d_v_div_r(s, k, ierr) result(d_v_div_r)
      type(star_info), pointer :: s
      integer, intent(in)  :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: d_v_div_r
      type(auto_diff_real_star_order1) :: r_cell, rho_cell, v_cell, dlnrho_dt
      real(dp) :: dm_cell
      ierr = 0

      ! shortcuts
      r_cell = 0.5d0*(wrap_r_00(s, k) + wrap_r_p1(s, k))
      rho_cell = wrap_d_00(s, k)
      v_cell = wrap_v_00(s, k)              ! cell-centred velocity (u_flag)
      dlnrho_dt = wrap_dxh_lnd(s, k)/s%dt    ! (∂/∂t)lnρ
      dm_cell = s%dm(k)                     ! cell mass

      ! Eq. (5)
      d_v_div_r = -dm_cell/(4d0*pi*rho_cell)*(dlnrho_dt/pow3(r_cell) + 3d0*v_cell/pow4(r_cell))

      ! units check:  (g) / (g cm) * (s⁻¹ cm⁻3) = s⁻¹        ✓
   end function compute_rho_form_of_d_v_div_r

   function compute_rho_form_of_d_v_div_r_opt_time_center(s, k, ierr) result(d_v_div_r) ! s^-1
      type(star_info), pointer :: s
      integer, intent(in)  :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: d_v_div_r
      type(auto_diff_real_star_order1) :: r_cell, rho_cell, v_cell, dlnrho_dt
      real(dp) :: dm_cell
      ierr = 0

      ! shortcuts -----------------------------------------------------------
      r_cell = 0.5d0*(wrap_opt_time_center_r_00(s, k) + wrap_opt_time_center_r_p1(s, k))
      rho_cell = wrap_d_00(s, k)
      v_cell = wrap_opt_time_center_v_00(s, k)              ! cell-centred velocity (u_flag)
      dlnrho_dt = wrap_dxh_lnd(s, k)/s%dt    ! (∂/∂t)lnρ
      dm_cell = s%dm(k)                     ! cell mass

      ! Eq. (5)
      d_v_div_r = -dm_cell/(4d0*pi*rho_cell)* &
                  (dlnrho_dt/pow3(r_cell) &
                   + 3d0*v_cell/pow4(r_cell))

      ! units check:  (g) / (g cm) * (s⁻¹ cm⁻3) = s⁻¹        ✓
   end function compute_rho_form_of_d_v_div_r_opt_time_center

   function compute_d_v_div_r_opt_time_center(s, k, ierr) result(d_v_div_r)  ! s^-1
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: d_v_div_r
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
      include 'formats'
      ierr = 0
      v_00 = wrap_opt_time_center_v_00(s, k)
      v_p1 = wrap_opt_time_center_v_p1(s, k)
      r_00 = wrap_opt_time_center_r_00(s, k)
      r_p1 = wrap_opt_time_center_r_p1(s, k)
      if (r_p1%val == 0d0) r_p1 = 1d0
      d_v_div_r = v_00/r_00 - v_p1/r_p1  ! units s^-1
   end function compute_d_v_div_r_opt_time_center

   function wrap_Hp_cell(s, k) result(Hp_cell)  ! cm
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: Hp_cell
      Hp_cell = 0.5d0*(wrap_Hp_00(s, k) + wrap_Hp_p1(s, k))
   end function wrap_Hp_cell

   function Hp_cell_for_Chi(s, k, ierr) result(Hp_cell)  ! cm
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Hp_cell
      type(auto_diff_real_star_order1) :: d_00, Peos_00, rmid
      real(dp) :: mmid, cgrav_mid
      include 'formats'
      ierr = 0

      Hp_cell = wrap_Hp_cell(s, k)
      return

      d_00 = wrap_d_00(s, k)
      Peos_00 = wrap_Peos_00(s, k)
      if (k < s%nz) then
         rmid = 0.5d0*(wrap_r_00(s, k) + wrap_r_p1(s, k))
         mmid = 0.5d0*(s%m(k) + s%m(k + 1))
         cgrav_mid = 0.5d0*(s%cgrav(k) + s%cgrav(k + 1))
      else
         rmid = 0.5d0*(wrap_r_00(s, k) + s%r_center)
         mmid = 0.5d0*(s%m(k) + s%m_center)
         cgrav_mid = s%cgrav(k)
      end if
      Hp_cell = pow2(rmid)*Peos_00/(d_00*cgrav_mid*mmid)
      if (s%alt_scale_height_flag) then
         call mesa_error(__FILE__, __LINE__, 'Hp_cell_for_Chi: cannot use alt_scale_height_flag')
      end if
   end function Hp_cell_for_Chi

   function compute_Chi_cell(s, k, ierr) result(Chi_cell)
      ! eddy viscosity energy (Kuhfuss 1986) [erg]
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: Chi_cell
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: &
         rho2, r6_cell, d_v_div_r, Hp_cell, w_00, d_00, r_00, r_p1
      real(dp) :: f, ALFAM_ALFA
      logical :: dbg
      include 'formats'
      ierr = 0
      dbg = .false.

      ! check where we are getting alfam from.
      if (s%MLT_option == 'TDC' .and. .not. s%RSP2_flag) then
         ALFAM_ALFA = s%alpha_TDC_DampM*s%mixing_length_alpha
      else if (s%RSP2_flag) then
         ALFAM_ALFA = s%RSP2_alfam*s%mixing_length_alpha
      else ! this is for safety, but probably is never called.
         ALFAM_ALFA = 0d0
      end if

      if (ALFAM_ALFA == 0d0 .or. &
          k <= s%RSP2_num_outermost_cells_forced_nonturbulent .or. &
          k > s%nz - int(s%nz/s%RSP2_nz_div_IBOTOM)) then
         Chi_cell = 0d0
         if (k >= 1 .and. k <= s%nz) then
            s%Chi(k) = 0d0
            s%Chi_ad(k) = 0d0
         end if
      else
         Hp_cell = Hp_cell_for_Chi(s, k, ierr)
         if (ierr /= 0) return
         if (s%u_flag .or. s%TDC_use_density_form_for_eddy_viscosity) then
            ! new density derivative term
            d_v_div_r = compute_rho_form_of_d_v_div_r(s, k, ierr)
         else
            d_v_div_r = compute_d_v_div_r(s, k, ierr)
         end if
         if (ierr /= 0) return

         ! don't need to check if mlt_vc > 0 here.
         if (s%MLT_option == 'TDC' .and. .not. s%RSP2_flag) then
            if (s%have_mlt_vc .and. s%okay_to_set_mlt_vc) then
               w_00 = s%mlt_vc_old(k)/sqrt_2_div_3! same as info%A0 from TDC
            else
               w_00 = s%mlt_vc(k)/sqrt_2_div_3! same as info%A0 from TDC
            end if
         else ! normal RSP2
            w_00 = wrap_w_00(s, k)
         end if
         d_00 = wrap_d_00(s, k)
         f = (16d0/3d0)*pi*ALFAM_ALFA/s%dm(k)
         rho2 = pow2(d_00)
         r_00 = wrap_r_00(s, k)
         r_p1 = wrap_r_p1(s, k)
         r6_cell = 0.5d0*(pow6(r_00) + pow6(r_p1))
         Chi_cell = f*rho2*r6_cell*d_v_div_r*Hp_cell*w_00
         ! units = g^-1 cm s^-1 g^2 cm^-6 cm^6 s^-1 cm
         !       = g cm^2 s^-2
         !       = erg

      end if
      s%Chi(k) = Chi_cell%val
      s%Chi_ad(k) = Chi_cell

      if (dbg .and. k == -100) then
         write (*, *) ' s% ALFAM_ALFA', ALFAM_ALFA
         write (*, *) 'Hp_cell', Hp_cell%val
         write (*, *) 'd_v_div_r', d_v_div_r%val
         write (*, *) ' f', f
         write (*, *) 'w_00', w_00%val
         write (*, *) 'd_00 ', d_00%val
         write (*, *) 'rho2 ', rho2%val
         write (*, *) 'r_00', r_00%val
         write (*, *) 'r_p1 ', r_p1%val
         write (*, *) 'r6_cell', r6_cell%val
      end if
   end function compute_Chi_cell

   function compute_tdc_Eq_cell(s, k, ierr) result(Eq_cell)  ! erg g^-1 s^-1
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: Eq_cell
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: d_v_div_r, Chi_cell
      include 'formats'
      ierr = 0
      if (s%mixing_length_alpha == 0d0 .or. &
          k <= s%RSP2_num_outermost_cells_forced_nonturbulent .or. &
          k > s%nz - int(s%nz/s%RSP2_nz_div_IBOTOM)) then
         Eq_cell = 0d0
         if (k >= 1 .and. k <= s%nz) s%Eq_ad(k) = 0d0
      else
         Chi_cell = s%Chi_ad(k)  ! compute_Chi_cell(s,k,ierr)
         if (ierr /= 0) return

         if (s%u_flag .or. s%TDC_use_density_form_for_eddy_viscosity) then
            ! new density derivative term
            d_v_div_r = compute_rho_form_of_d_v_div_r_opt_time_center(s, k, ierr)
         else
            d_v_div_r = compute_d_v_div_r_opt_time_center(s, k, ierr)
         end if

         if (ierr /= 0) return
         Eq_cell = 4d0*pi*Chi_cell*d_v_div_r/s%dm(k)  ! erg s^-1 g^-1
      end if
      s%Eq(k) = Eq_cell%val
      s%Eq_ad(k) = Eq_cell
   end function compute_tdc_Eq_cell

   function compute_tdc_Eq_face(s, k, ierr) result(Eq_face)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Eq_face
      real(dp) :: alfa, beta
      include 'formats'
      ierr = 0
      if (k == 1) then
         Eq_face = 0d0
      else
         call get_RSP2_alfa_beta_face_weights(s, k, alfa, beta)
         Eq_face = alfa*compute_tdc_Eq_cell(s, k, ierr) + beta*compute_tdc_Eq_cell(s, k - 1, ierr) ! should it be k and k+1?
      end if
      if (ierr /= 0) return
   end function compute_tdc_Eq_face

   function compute_tdc_Uq_face(s, k, ierr) result(Uq_face)  ! cm s^-2, acceleration
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: Uq_face
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Chi_00, Chi_m1, r_00
      include 'formats'
      ierr = 0
      if (s%mixing_length_alpha == 0d0 .or. &
          k <= s%RSP2_num_outermost_cells_forced_nonturbulent .or. &
          k > s%nz - int(s%nz/s%RSP2_nz_div_IBOTOM)) then
         Uq_face = 0d0
      else
         r_00 = wrap_opt_time_center_r_00(s, k)

         ! which do we adopt?
         Chi_00 = compute_Chi_cell(s, k, ierr)  ! s% Chi_ad(k) XXX
         !Chi_00 = s% Chi_ad(k)  ! compute_Chi_cell(s,k,ierr)

         if (k > 1) then
            Chi_m1 = shift_m1(compute_Chi_cell(s, k - 1, ierr))
            !Chi_m1 = shift_m1(s% Chi_ad(k-1)) XXX
            if (ierr /= 0) return
         else
            Chi_m1 = 0d0
         end if
         Uq_face = 4d0*pi*(Chi_m1 - Chi_00)/(r_00*s%dm_bar(k))

         if (k == -56) then
            write (*, 3) 'RSP2 Uq chi_m1 chi_00 r', k, s%solver_iter, &
               Uq_face%val, Chi_m1%val, Chi_00%val, r_00%val
         end if

      end if
      ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2, acceleration
      s%Uq(k) = Uq_face%val
   end function compute_tdc_Uq_face

end module tdc_pulse
