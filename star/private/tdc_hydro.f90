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

module tdc_hydro

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
      compute_tdc_Uq_face, compute_tdc_Eq_div_w_face, &
      get_TDC_alfa_beta_face_weights, set_viscosity_vars_TDC, compute_tdc_Uq_dm_cell

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

      if (.not. (s%v_flag .or. s%u_flag)) then ! set values 0 if not using v_flag or u_flag.
         do k = 1, s%nz
            s%Eq(k) = 0d0; s%Eq_ad(k) = 0d0
            s%Chi(k) = 0d0; s%Chi_ad(k) = 0d0
            s%Uq(k) = 0d0
         end do
         return
      end if

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
         x = compute_Chi_div_w_face(s, k, op_err) ! Sets Chi_face
         if (op_err /= 0) ierr = op_err
         x = compute_tdc_Eq_div_w_face(s, k, op_err) ! Sets Eq_face
         if (op_err /= 0) ierr = op_err
         if (s% v_flag) then
            x = compute_tdc_Uq_face(s, k, op_err)
         else if (s% u_flag) then
            x = compute_tdc_Uq_dm_cell(s, k, op_err)
         end if
         if (op_err /= 0) ierr = op_err
      end do
      !$OMP END PARALLEL DO
      if (ierr /= 0) then
         if (s%report_ierr) write (*, 2) 'failed in set_viscosity_vars_TDC loop 2', s%model_number
         return
      end if
   end subroutine set_viscosity_vars_TDC

   subroutine get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
      type(star_info), pointer :: s
      integer, intent(in) :: k
      real(dp), intent(out) :: alfa, beta
      ! face_value = alfa*cell_value(k) + beta*cell_value(k-1)
      if (k == 1) call mesa_error(__FILE__, __LINE__, 'bad k==1 for get_TDC_alfa_beta_face_weights')
      if (s%TDC_hydro_use_mass_interp_face_values) then
         alfa = s%dq(k - 1)/(s%dq(k - 1) + s%dq(k))
         beta = 1d0 - alfa
      else
         alfa = 0.5d0
         beta = 0.5d0
      end if
   end subroutine get_TDC_alfa_beta_face_weights


   function wrap_Hp_cell(s, k) result(Hp_cell)  ! cm , different than rsp2
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: Hp1, Hp0, Hp_cell
      Hp0 = get_scale_height_face(s,k)
      Hp1 = 0d0
      if (k+1 < s%nz) then
         Hp1 = shift_p1(get_scale_height_face(s,k+1))
      end if
      Hp_cell = 0.5d0*(Hp0 + Hp1)
      !0.5d0*(wrap_Hp_00(s, k) + wrap_Hp_p1(s, k))
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
      return ! below is skipped, for now.

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

   ! this function is only called internally in TDC_Uq_face, and for v_flag only.
   function compute_Chi_cell(s, k, ierr) result(Chi_cell) ! does not update s% Chi or Chi_ad
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
         ALFAM_ALFA = s%TDC_alpha_M*s%mixing_length_alpha
      else ! this is for safety, but probably is never called.
         ALFAM_ALFA = 0d0
      end if

      if (ALFAM_ALFA == 0d0 .or. &
          k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
          k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
         Chi_cell = 0d0
      else
         Hp_cell = Hp_cell_for_Chi(s, k, ierr)
         if (ierr /= 0) return
         if (s%TDC_use_density_form_for_eddy_viscosity) then
            ! new density derivative term
            d_v_div_r = compute_rho_form_of_d_v_div_r(s, k, ierr)
         else
            d_v_div_r = compute_d_v_div_r(s, k, ierr)
         end if
         if (ierr /= 0) return

         ! don't need to check if mlt_vc > 0 here.
         if (k < s% nz) then
            if (s% okay_to_set_mlt_vc .and. &
               s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then !add option for explicit mlt_vc, operator split in momentum eq.
               w_00 = 0.5d0*(s% mlt_vc_old(k) + s% mlt_vc_old(k+1))/sqrt_2_div_3! same as info%A0 from TDC
            else
               w_00 = 0.5d0*(s% mlt_vc_ad(k) + shift_p1(s% mlt_vc_ad(k+1)))/sqrt_2_div_3! same as info%A0 from TDC
            end if
         else
            if (s% okay_to_set_mlt_vc .and. &
                s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then !add option for explicit mlt_vc, operator split in momentum eq.
               w_00 = 0.5d0*s% mlt_vc_old(k)/sqrt_2_div_3! same as info%A0 from TDC
            else
               w_00 = 0.5d0*s% mlt_vc_ad(k)/sqrt_2_div_3! same as info%A0 from TDC
            end if
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
      ! this is set in Chi_div_w_face
      !s%Chi(k) = Chi_cell%val
      !s%Chi_ad(k) = Chi_cell

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

  ! face centered variables for tdc update below
   function compute_Chi_div_w_face(s, k, ierr) result(Chi_face)
   ! eddy viscosity energy (Kuhfuss 1986) [erg]
   type(star_info), pointer :: s
   integer, intent(in) :: k
   type(auto_diff_real_star_order1) :: Chi_face
   integer, intent(out) :: ierr
   type(auto_diff_real_star_order1) :: &
   rho2, r6_face, d_v_div_r, Hp_face, w_00, d_00, r_00, r_p1
   real(dp) :: f, ALFAM_ALFA, dmbar
   logical :: dbg
   include 'formats'
   ierr = 0
   dbg = .false.

   ! check where we are getting alfam from.
   if (s%MLT_option == 'TDC' .and. .not. s%RSP2_flag) then
      ALFAM_ALFA = s%TDC_alpha_M*s%mixing_length_alpha
   else ! this is for safety, but probably is never called.
      ALFAM_ALFA = 0d0
   end if

   if (ALFAM_ALFA == 0d0 .or. &
      k > s%nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
      Chi_face = 0d0
   else
      Hp_face = get_scale_height_face(s,k) !Hp_cell_for_Chi(s, k, ierr)
      if (ierr /= 0) return
      if (s%TDC_use_density_form_for_eddy_viscosity) then
         ! new density derivative form
         d_v_div_r = compute_rho_form_of_d_v_div_r_face(s, k, ierr)
      else
         d_v_div_r = compute_d_v_div_r_face(s, k, ierr)
      end if
      if (ierr /= 0) return

      if (k >= 2) then
         dmbar = 0.5d0*(s% dm(k) + s% dm(k-1))
      else
         dmbar = 0.5d0*s% dm(k)
      end if
      d_00 = get_rho_face(s, k)
      f = (16d0/3d0)*pi*ALFAM_ALFA/dmbar
      rho2 = pow2(d_00)
      r_00 = wrap_r_00(s, k)
      !r_p1 = wrap_r_p1(s, k)
      r6_face = pow6(r_00) !0.5d0*(pow6(r_00) + pow6(r_p1))
      Chi_face = f*rho2*r6_face*d_v_div_r*Hp_face!*w_00
      ! units = g^-1 cm s^-1 g^2 cm^-6 cm^6 s^-1 cm * [s/cm] ! [1/w_00] = [s/cm]
      !       = g cm^2 s^-2 * [s/cm]
      !       = erg ! * [s / cm] - > [erg] * [s/cm]

   end if

   ! Chi_cell does not set Chi, we store Chi_face in s% Chi and s% Chi_ad
      if (s% okay_to_set_mlt_vc .and. &
         s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then !add option for explicit mlt_vc, operator split in momentum eq.
         w_00 = s% mlt_vc_old(k)/sqrt_2_div_3! same as info%A0 from TDC
      else
         w_00 = s% mlt_vc_ad(k)/sqrt_2_div_3! same as info%A0 from TDC
      end if
      s%Chi(k) = Chi_face%val*w_00%val
      s%Chi_ad(k) = Chi_face*w_00

      if (dbg .and. k == -100) then
      write (*, *) ' s% ALFAM_ALFA', ALFAM_ALFA
      write (*, *) 'Hp_face', Hp_face%val
      write (*, *) 'd_v_div_r', d_v_div_r%val
      write (*, *) ' f', f
      write (*, *) 'w_00', w_00%val
      write (*, *) 'd_00 ', d_00%val
      write (*, *) 'rho2 ', rho2%val
      write (*, *) 'r_00', r_00%val
      write (*, *) 'r_p1 ', r_p1%val
      write (*, *) 'r6_cell', r6_face%val
      end if
   end function compute_Chi_div_w_face

   function compute_tdc_Eq_div_w_face(s, k, ierr) result(Eq_face)  ! erg g^-1 s^-1 * (cm^-1 s^1)
   type(star_info), pointer :: s
   integer, intent(in) :: k
   type(auto_diff_real_star_order1) :: Eq_face
   integer, intent(out) :: ierr
   type(auto_diff_real_star_order1) :: d_v_div_r, Chi_face, w_00
   real(dp) :: dmbar
   include 'formats'
   ierr = 0
   if (s%mixing_length_alpha == 0d0 .or. &
   k > s%nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
      Eq_face = 0d0
      if (k >= 1 .and. k <= s%nz) s%Eq_ad(k) = 0d0
   else
      Chi_face = compute_Chi_div_w_face(s,k,ierr)
      if (ierr /= 0) return

      if (s%TDC_use_density_form_for_eddy_viscosity) then
         ! new density derivative term
         d_v_div_r = compute_rho_form_of_d_v_div_r_face_opt_time_center(s, k, ierr)
      else
         d_v_div_r = compute_d_v_div_r_opt_time_center_face(s, k, ierr)
      end if

      if (k >= 2) then
         dmbar = 0.5d0*(s% dm(k) + s% dm(k-1))
      else
         dmbar = 0.5d0*s% dm(k)
      end if

      if (ierr /= 0) return
      Eq_face = 4d0*pi*Chi_face*d_v_div_r/dmbar  ! erg s^-1 g^-1 * (cm^-1 s^1)
   end if

   ! only for output, really only used for returning Eq to star pointers.
   if (s% okay_to_set_mlt_vc .and. &
      s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then !add option for explicit mlt_vc, operator split in momentum eq.
      w_00 = s% mlt_vc_old(k)/sqrt_2_div_3! same as info%A0 from TDC
   else
      w_00 = s% mlt_vc_ad(k)/sqrt_2_div_3! same as info%A0 from TDC
   end if

   s%Eq(k) = Eq_face%val * w_00%val
   s%Eq_ad(k) = Eq_face * w_00
   end function compute_tdc_Eq_div_w_face

   ! for v_flag only. face centered Uq for hydro_momentum
   function compute_tdc_Uq_face(s, k, ierr) result(Uq_face) !(v_flag only)  ! cm s^-2, acceleration
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: Uq_face
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Chi_00, Chi_m1, r_00
      include 'formats'
      ierr = 0
      if (s%mixing_length_alpha == 0d0 .or. &
          k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
          k > s%nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
         Uq_face = 0d0
      else
         r_00 = wrap_opt_time_center_r_00(s, k)

         ! which do we adopt?
         Chi_00 = compute_Chi_cell(s, k, ierr)  ! s% Chi_ad(k) XXX

         if (k > 1) then
            Chi_m1 = shift_m1(compute_Chi_cell(s, k-1, ierr))
            if (ierr /= 0) return
         else
            Chi_m1 = 0d0
         end if
         Uq_face = 4d0*pi*(Chi_m1 - Chi_00)/(r_00*s%dm_bar(k))

         if (k == -56) then
            write (*, 3) 'TDC Uq chi_m1 chi_00 r', k, s%solver_iter, &
               Uq_face%val, Chi_m1%val, Chi_00%val, r_00%val
         end if

      end if
      ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2, acceleration
      s%Uq(k) = Uq_face%val
   end function compute_tdc_Uq_face

   ! for u_flag only. cell centered Uq as source for Reimann flux.
   function compute_tdc_Uq_dm_cell(s, k, ierr) result(Uq_cell)  ! cm s^-2, acceleration
      type(star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: Chi_00, Chi_p1, r_00, r_p1, w_00, w_p1, r_cell, Uq_cell
      include 'formats'
      ierr = 0
      if (s%mixing_length_alpha == 0d0 .or. &
          k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
          k > s%nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
         Uq_cell = 0d0
      else
         r_00 = wrap_opt_time_center_r_00(s, k)
         r_p1 = wrap_opt_time_center_r_p1(s, k)
         r_cell = 0.5d0*(r_00+r_p1) ! not staggered unlike terms inside chi_div_w_face

         if (s% okay_to_set_mlt_vc .and. &
            s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then
            w_00 = s% mlt_vc_old(k)/sqrt_2_div_3
         else
            w_00 = s% mlt_vc_ad(k)/sqrt_2_div_3
         end if

         Chi_00 = compute_Chi_div_w_face(s, k, ierr) * w_00

         if (k < s% nz) then
            if (s% okay_to_set_mlt_vc .and. &
               s% TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation) then
               w_p1 = s% mlt_vc_old(k+1)/sqrt_2_div_3
            else
               w_p1 = shift_p1(s% mlt_vc_ad(k+1))/sqrt_2_div_3
            end if

            Chi_p1 = shift_p1(compute_Chi_div_w_face(s, k+1, ierr))*w_p1
            if (ierr /= 0) return
         else
            Chi_p1 = 0d0
            w_p1 = 0d0
         end if

         Uq_cell = 4d0*pi*(Chi_00 - Chi_p1)/(r_cell) ! we have neglected the /dm here, because it is restored in the reimann flux calculation
         ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2 [g], acceleration*mass = Force

         if (k == -56) then
            write (*, 3) 'TDC Uq chi_m1 chi_00 r', k, s%solver_iter, &
               Uq_cell%val, Chi_p1%val, Chi_00%val, r_00%val
         end if

      end if
      s%Uq(k) = Uq_cell%val/ s% dm(k)
   end function compute_tdc_Uq_dm_cell


! all the forms of d(v/r)/dr, below
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

   function compute_rho_form_of_d_v_div_r(s, k, ierr) result(d_v_div_r) ! used in Chi_cell
      type(star_info), pointer :: s
      integer, intent(in)  :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: d_v_div_r, v_00, v_p1
      type(auto_diff_real_star_order1) :: r_cell, rho_cell, v_cell, dlnrho_dt
      real(dp) :: dm_cell
      ierr = 0

      r_cell = 0.5d0*(wrap_r_00(s, k) + wrap_r_p1(s, k))
      rho_cell = wrap_d_00(s, k)
      if (s% u_flag) then
         v_cell = wrap_u_00(s,k)
      else ! v flag
         v_cell = 0.5d0*(wrap_v_00(s, k) + wrap_v_p1(s, k))
      end if
      v_00 = wrap_opt_time_center_v_00(s, k)
      v_p1 = wrap_opt_time_center_v_p1(s, k)
      dlnrho_dt = wrap_dxh_lnd(s, k)/s%dt    ! (∂/∂t)lnρ
      dm_cell = s%dm(k)                     ! cell mass

      ! density form
      d_v_div_r = -dm_cell/(4d0*pi*rho_cell)*(dlnrho_dt/pow3(r_cell) + 3d0*v_cell/pow4(r_cell))

      ! dm_cell*(1/r * du/dm - U/4/pi/rho/r^4), more sensitive to geometry
      !d_v_div_r = ((v_00 - v_p1) - dm_cell*v_cell/(4d0*pi*rho_cell*pow3(r_cell)))/r_cell

   end function compute_rho_form_of_d_v_div_r

   function compute_rho_form_of_d_v_div_r_face(s, k, ierr) result(d_v_div_r)
      type(star_info), pointer :: s
      integer, intent(in)  :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: d_v_div_r
      type(auto_diff_real_star_order1) :: r_face, rho_face, v_face, dlnrho_dt
      real(dp) :: dm_bar, alfa, beta
      ierr = 0

      r_face = wrap_r_00(s, k)
      rho_face = get_rho_face(s, k)
      v_face = wrap_v_00(s, k)   ! face-centered velocity
      if (k >= 2) then
         dm_bar = 0.5d0*(s% dm(k) + s% dm(k-1))
         call get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
         dlnrho_dt = (alfa*wrap_dxh_lnd(s, k) + beta*shift_m1(wrap_dxh_lnd(s, k-1)))/s%dt    ! (∂/∂t)lnρ
      else
         dm_bar = 0.5d0*s% dm(k)
         dlnrho_dt = 0.5d0*wrap_dxh_lnd(s, k)/s%dt    ! (∂/∂t)lnρ
      end if

      ! density form
      d_v_div_r = -dm_bar/(4d0*pi*rho_face)*(dlnrho_dt/pow3(r_face) + 3d0*v_face/pow4(r_face))

      ! dm_bar*(1/r * du/dm - U/4/pi/rho/r^4), more sensitive to geometry
      !d_v_div_r = ((wrap_u_m1(s,k) - wrap_u_00(s,k)) - dm_bar*v_face/(4d0*pi*rho_face*pow3(r_face)))/r_face

   end function compute_rho_form_of_d_v_div_r_face

   function compute_rho_form_of_d_v_div_r_face_opt_time_center(s, k, ierr) result(d_v_div_r) ! s^-1
      type(star_info), pointer :: s
      integer, intent(in)  :: k
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: d_v_div_r
      type(auto_diff_real_star_order1) :: r_face, rho_face, v_face, dlnrho_dt
      real(dp) :: dm_bar, alfa, beta
      ierr = 0

      r_face = wrap_opt_time_center_r_00(s, k)
      rho_face = get_rho_face(s, k)
      v_face = wrap_opt_time_center_v_00(s, k)   ! face-centered velocity
      if (k >= 2) then
         dm_bar = 0.5d0*(s% dm(k) + s% dm(k-1))
         call get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
         dlnrho_dt = (alfa*wrap_dxh_lnd(s, k) + beta*shift_m1(wrap_dxh_lnd(s, k-1)))/s%dt    ! (∂/∂t)lnρ
      else
         dm_bar = 0.5d0*s% dm(k)
         dlnrho_dt = 0.5d0*wrap_dxh_lnd(s, k)/s%dt    ! (∂/∂t)lnρ
      end if

      ! density form
      d_v_div_r = -dm_bar/(4d0*pi*rho_face)*(dlnrho_dt/pow3(r_face) + 3d0*v_face/pow4(r_face))

      ! dm_bar*(1/r * du/dm - U/4/pi/rho/r^4), more sensitive to geometry
      !d_v_div_r = ((wrap_opt_time_center_u_m1(s,k) - wrap_opt_time_center_u_00(s,k)) - dm_bar*v_face/(4d0*pi*rho_face*pow3(r_face)))/r_face

   end function compute_rho_form_of_d_v_div_r_face_opt_time_center

   function compute_d_v_div_r_face(s, k, ierr) result(d_v_div_r)  ! s^-1
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: d_v_div_r
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: v_00, v_m1, r_00, r_m1, term1, term2
      logical :: dbg
      include 'formats'
      ierr = 0
      dbg = .false.

     if (s% v_flag) then
         v_00 = 0.5d0*(wrap_v_00(s, k) + wrap_v_p1(s, k))
         v_m1 = 0.5d0*(wrap_v_00(s, k) + wrap_v_m1(s, k))
     else if(s% u_flag) then
         v_00 = wrap_u_00(s,k)
         v_m1 = wrap_u_m1(s,k)
      end if

      if (s% v_flag) then
         r_00 = 0.5d0*(wrap_r_00(s, k) + wrap_r_p1(s, k))
         r_m1 = 0.5d0*(wrap_r_00(s, k) + wrap_r_m1(s, k))
      else if(s% u_flag) then ! stagger r for u_flag to retain tridiagonality.
         r_00 = wrap_r_00(s, k)
         r_m1 = wrap_r_m1(s, k)
      end if

      if (r_00%val == 0d0) r_00 = 1d0
      if (r_m1%val == 0d0) r_m1 = 1d0
      d_v_div_r = v_m1/r_m1 - v_00/r_00 ! units s^-1

      ! Debugging output to trace values
      if (dbg .and. k == -63) then
         write (*, *) 'test d_v_div_r, k:', k
         write (*, *) 'v_00:', v_00%val, 'v_p1:', v_m1%val
         write (*, *) 'r_00:', r_00%val, 'r_p1:', r_m1%val
         write (*, *) 'd_v_div_r:', d_v_div_r%val
      end if
   end function compute_d_v_div_r_face

   function compute_d_v_div_r_opt_time_center_face(s, k, ierr) result(d_v_div_r)  ! s^-1
      type(star_info), pointer :: s
      integer, intent(in) :: k
      type(auto_diff_real_star_order1) :: d_v_div_r
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: v_00, v_m1, r_00, r_m1, term1, term2
      logical :: dbg
      include 'formats'
      ierr = 0
      dbg = .false.

     if (s% v_flag) then
         v_00 = 0.5d0 *(wrap_opt_time_center_v_00(s, k) + wrap_opt_time_center_v_p1(s, k))
         v_m1 = 0.5d0*(wrap_opt_time_center_v_00(s, k) + wrap_opt_time_center_v_m1(s, k))
     else if(s% u_flag) then
         v_00 = wrap_opt_time_center_u_00(s,k)
         v_m1 = wrap_opt_time_center_u_m1(s,k)
      end if

      if (s% v_flag) then
         r_00 = 0.5d0*(wrap_opt_time_center_r_00(s, k) + wrap_opt_time_center_r_p1(s, k))
         r_m1 = 0.5d0*(wrap_opt_time_center_r_00(s, k) + wrap_opt_time_center_r_m1(s, k))
      else if(s% u_flag) then ! stagger r for u_flag to retain tridiagonality.
         r_00 = wrap_opt_time_center_r_00(s, k)
         r_m1 = wrap_opt_time_center_r_m1(s, k)
      end if

      if (r_00%val == 0d0) r_00 = 1d0
      if (r_m1%val == 0d0) r_m1 = 1d0
      d_v_div_r = v_m1/r_m1 - v_00/r_00 ! units s^-1

      ! Debugging output to trace values
      if (dbg .and. k == -63) then
         write (*, *) 'test d_v_div_r, k:', k
         write (*, *) 'v_00:', v_00%val, 'v_p1:', v_m1%val
         write (*, *) 'r_00:', r_00%val, 'r_p1:', r_m1%val
         write (*, *) 'd_v_div_r:', d_v_div_r%val
      end if
   end function compute_d_v_div_r_opt_time_center_face

end module tdc_hydro
