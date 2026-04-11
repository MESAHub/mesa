! ***********************************************************************
!
!   Copyright (C) 2010-2025  The MESA Team
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

module mlt_tdc_face_support

  use star_private_def
  use const_def, only: dp, ln10, pi4, clight, crad
  use auto_diff
  use kap_support, only: get_kap

  implicit none

  private
  public :: get_mlt_face_state_ad
  public :: get_face_eos_kap_ad
  public :: get_face_scale_height_ad

contains

  ! Returns the MLT/TDC face thermodynamic state as
  ! auto_diff_real_star_order1 quantities, either from recomputed face
  ! EOS/opacity data or from the stored face quantities.
  subroutine get_mlt_face_state_ad( &
       s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, &
       opacity_face, scale_height_face, gradr_face, ierr)
    use star_utils, only: get_T_face, get_Peos_face, get_kap_face, get_rho_face, &
       get_ChiRho_face, get_ChiT_face, get_Cp_face, get_grada_face, get_scale_height_face, get_gradr_face

    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(out) :: &
       T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, &
       opacity_face, scale_height_face, gradr_face
    integer, intent(out) :: ierr

    ierr = 0
    if (s%use_face_values_eos_and_kap_mlt_tdc) then
       call get_face_eos_kap_ad( &
          s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, opacity_face, ierr)
       if (ierr /= 0) return
       call set_scale_height_from_face_state(s, k, P_face, rho_face, scale_height_face)
       call set_gradr_from_face_state(s, k, P_face, opacity_face, T_face, gradr_face)
    else
       T_face = get_T_face(s, k)
       P_face = get_Peos_face(s, k)
       opacity_face = get_kap_face(s, k)
       rho_face = get_rho_face(s, k)
       ChiRho_face = get_ChiRho_face(s, k)
       ChiT_face = get_ChiT_face(s, k)
       Cp_face = get_Cp_face(s, k)
       grada_face = get_grada_face(s, k)
       scale_height_face = get_scale_height_face(s, k)
       gradr_face = get_gradr_face(s, k)
    end if
  end subroutine get_mlt_face_state_ad

  subroutine get_face_composition(s, k, use_starting_comp, zbar_face, xa_face, ierr)
    use star_utils, only: get_face_weights

    type(star_info), pointer :: s
    integer, intent(in) :: k
    logical, intent(in) :: use_starting_comp
    real(dp), intent(out) :: zbar_face
    real(dp), intent(out) :: xa_face(:)
    integer, intent(out) :: ierr

    real(dp) :: alfa, beta, sum_xa

    ierr = 0
    if (k == 1) then
       alfa = 1d0
       beta = 0d0
    else
       call get_face_weights(s, k, alfa, beta)
    end if

    if (use_starting_comp) then
       zbar_face = alfa*s%zbar_start(k) + beta*s%zbar_start(max(1, k-1))
       xa_face(1:s%species) = alfa*s%xa_start(1:s%species, k)
       if (k > 1) xa_face(1:s%species) = xa_face(1:s%species) + beta*s%xa_start(1:s%species, k-1)
    else
       zbar_face = alfa*s%zbar(k) + beta*s%zbar(max(1, k-1))
       xa_face(1:s%species) = alfa*s%xa(1:s%species, k)
       if (k > 1) xa_face(1:s%species) = xa_face(1:s%species) + beta*s%xa(1:s%species, k-1)
    end if

    sum_xa = sum(xa_face)
    if (sum_xa <= 0d0) then
      ierr = -1
      return
    end if
    xa_face = xa_face/sum_xa
  end subroutine get_face_composition


  subroutine get_face_eos_inputs( &
       s, k, T_face, rho_face, eos_res, d_dlnd, d_dlnT, ierr)
    use eos_support, only: get_eos
    use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results
    use star_utils, only: get_rho_face, get_T_face

    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(out) :: T_face, rho_face
    real(dp), intent(out) :: eos_res(num_eos_basic_results), d_dlnd(num_eos_basic_results), d_dlnT(num_eos_basic_results)
    integer, intent(out) :: ierr

    real(dp) :: log10_T, log10_rho, zbar_face
    real(dp) :: eos_xa_face(s%species)
    real(dp) :: d_dxa(num_eos_d_dxa_results, s%species)

    ierr = 0
    T_face = get_T_face(s, k)
    rho_face = get_rho_face(s, k)
    if (T_face%val <= 0d0 .or. rho_face%val <= 0d0) then
       ierr = -1
       return
    end if

    call get_face_composition(s, k, .false., zbar_face, eos_xa_face, ierr)
    if (ierr /= 0) return

    log10_T = log10(T_face%val)
    log10_rho = log10(rho_face%val)

    call get_eos( &
       s, k, eos_xa_face, rho_face%val, log10_rho, T_face%val, log10_T, &
       eos_res, d_dlnd, d_dlnT, d_dxa, ierr)
    if (ierr /= 0) return
  end subroutine get_face_eos_inputs


  subroutine get_face_opacity_factor(s, k, log10_T, opacity_factor_face)
    use star_utils, only: get_face_weights

    type(star_info), pointer :: s
    integer, intent(in) :: k
    real(dp), intent(in) :: log10_T
    real(dp), intent(out) :: opacity_factor_face

    real(dp) :: alfa, beta

    if (k == 1) then
       alfa = 1d0
       beta = 0d0
    else
       call get_face_weights(s, k, alfa, beta)
    end if

    opacity_factor_face = alfa*s%extra_opacity_factor(k)
    if (k > 1) opacity_factor_face = opacity_factor_face + beta*s%extra_opacity_factor(k-1)
    if (s%min_logT_for_opacity_factor_off > 0) then
       if (log10_T >= s%max_logT_for_opacity_factor_off .or. &
             log10_T <= s%min_logT_for_opacity_factor_off) then
          opacity_factor_face = 1d0
       else if (log10_T > s%max_logT_for_opacity_factor_on) then
          opacity_factor_face = 1d0 + (opacity_factor_face - 1d0)* &
             (log10_T - s%max_logT_for_opacity_factor_off)/ &
             (s%max_logT_for_opacity_factor_on - s%max_logT_for_opacity_factor_off)
       else if (log10_T < s%min_logT_for_opacity_factor_on) then
          opacity_factor_face = 1d0 + (opacity_factor_face - 1d0)* &
             (log10_T - s%min_logT_for_opacity_factor_off)/ &
             (s%min_logT_for_opacity_factor_on - s%min_logT_for_opacity_factor_off)
       end if
    end if
  end subroutine get_face_opacity_factor


  subroutine set_face_ad_from_value(value, dvalue_dlnd, dvalue_dlnT, T_face, rho_face, quantity_ad)
    type(auto_diff_real_star_order1), intent(in) :: T_face, rho_face
    real(dp), intent(in) :: value, dvalue_dlnd, dvalue_dlnT
    type(auto_diff_real_star_order1), intent(out) :: quantity_ad
    type(auto_diff_real_star_order1) :: lnT_face, lnd_face
    integer :: j

    lnd_face = log(rho_face)
    lnT_face = log(T_face)
    quantity_ad = 0d0
    quantity_ad%val = value
    do j = 1, size(quantity_ad%d1Array)
       quantity_ad%d1Array(j) = dvalue_dlnd*lnd_face%d1Array(j) + dvalue_dlnT*lnT_face%d1Array(j)
    end do
  end subroutine set_face_ad_from_value


  subroutine set_face_ad_from_log(log_value, dlog_dlnd, dlog_dlnT, T_face, rho_face, quantity_ad)
    type(auto_diff_real_star_order1), intent(in) :: T_face, rho_face
    real(dp), intent(in) :: log_value, dlog_dlnd, dlog_dlnT
    type(auto_diff_real_star_order1), intent(out) :: quantity_ad
    real(dp) :: value

    value = exp(log_value)
    call set_face_ad_from_value(value, value*dlog_dlnd, value*dlog_dlnT, T_face, rho_face, quantity_ad)
  end subroutine set_face_ad_from_log


  ! Recomputes the face EOS and opacity state as
  ! auto_diff_real_star_order1 quantities for the MLT/TDC solve,
  ! instead of using the stored face quantities.
  subroutine get_face_eos_kap_ad( &
       s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, opacity_face, ierr)
    use eos_def, only: num_eos_basic_results, i_lnPgas, i_grad_ad, i_gamma1, i_Cp, i_chiRho, i_chiT, i_eta, i_lnfree_e
    use kap_def, only: num_kap_fracs

    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(out) :: T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, opacity_face
    integer, intent(out) :: ierr

    real(dp) :: log10_T, log10_rho, kap_zbar_face, opacity_factor_face
    real(dp) :: eos_res(num_eos_basic_results), d_dlnd(num_eos_basic_results), d_dlnT(num_eos_basic_results)
    real(dp) :: kap, dlnkap_dlnd, dlnkap_dlnT
    real(dp) :: kap_fracs(num_kap_fracs), kap_xa_face(s%species)
    type(auto_diff_real_star_order1) :: Pgas_face, gamma1_face, mlt_Pturb_ad, alpha

    ierr = 0
    call get_face_eos_inputs( &
       s, k, T_face, rho_face, eos_res, d_dlnd, d_dlnT, ierr)
    if (ierr /= 0) return
    log10_T = log10(T_face%val)
    log10_rho = log10(rho_face%val)

    call set_face_ad_from_log(eos_res(i_lnPgas), d_dlnd(i_lnPgas), d_dlnT(i_lnPgas), T_face, rho_face, Pgas_face)
    P_face = Pgas_face + crad*pow4(T_face)/3d0
    call set_face_ad_from_value(eos_res(i_Cp), d_dlnd(i_Cp), d_dlnT(i_Cp), T_face, rho_face, Cp_face)
    call set_face_ad_from_value(eos_res(i_chiRho), d_dlnd(i_chiRho), d_dlnT(i_chiRho), T_face, rho_face, ChiRho_face)
    call set_face_ad_from_value(eos_res(i_chiT), d_dlnd(i_chiT), d_dlnT(i_chiT), T_face, rho_face, ChiT_face)
    call set_face_ad_from_value(eos_res(i_grad_ad), d_dlnd(i_grad_ad), d_dlnT(i_grad_ad), T_face, rho_face, grada_face)
    call set_face_ad_from_value(eos_res(i_gamma1), d_dlnd(i_gamma1), d_dlnT(i_gamma1), T_face, rho_face, gamma1_face)

    call get_face_composition(s, k, s% use_starting_composition_for_kap, kap_zbar_face, kap_xa_face, ierr)
    if (ierr /= 0) return
    call get_face_opacity_factor(s, k, log10_T, opacity_factor_face)

    call get_kap( &
       s, k, kap_zbar_face, kap_xa_face, log10_rho, log10_T, &
       eos_res(i_lnfree_e), d_dlnd(i_lnfree_e), d_dlnT(i_lnfree_e), &
       eos_res(i_eta), d_dlnd(i_eta), d_dlnT(i_eta), &
       kap_fracs, kap, dlnkap_dlnd, dlnkap_dlnT, ierr)
    if (ierr /= 0 .or. is_bad_num(kap) .or. kap <= 0d0) return

    kap = kap*opacity_factor_face
    if (s%opacity_max > 0d0 .and. kap > s%opacity_max) then
       kap = s%opacity_max
       dlnkap_dlnd = 0d0
       dlnkap_dlnT = 0d0
    end if
    if (s%opacity_min > 0d0 .and. kap < s%opacity_min) then
       kap = s%opacity_min
       dlnkap_dlnd = 0d0
       dlnkap_dlnT = 0d0
    end if
    call set_face_ad_from_value(kap, kap*dlnkap_dlnd, kap*dlnkap_dlnT, T_face, rho_face, opacity_face)

    if (s% have_mlt_vc .and. s% okay_to_set_mlt_vc .and. s% include_mlt_Pturb_in_thermodynamic_gradients &
       .and. s% mlt_Pturb_factor > 0d0 .and. k > 1) then
       mlt_Pturb_ad = s% mlt_Pturb_factor*pow2(s% mlt_vc_old(k))*rho_face/3d0
       alpha = mlt_Pturb_ad/(P_face*gamma1_face)
       grada_face = grada_face*(P_face + mlt_Pturb_ad)/(P_face*(1d0 + alpha))
    end if
  end subroutine get_face_eos_kap_ad


  ! Recomputes the face pressure scale height as an
  ! auto_diff_real_star_order1 quantity from the recomputed face state.
  subroutine get_face_scale_height_ad(s, k, scale_height_face, ierr)
    use eos_def, only: num_eos_basic_results, i_lnPgas

    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(out) :: scale_height_face
    integer, intent(out) :: ierr

    real(dp) :: eos_res(num_eos_basic_results), d_dlnd(num_eos_basic_results), d_dlnT(num_eos_basic_results)
    type(auto_diff_real_star_order1) :: T_face, rho_face, Pgas_face, P_face

    ierr = 0
    call get_face_eos_inputs( &
       s, k, T_face, rho_face, eos_res, d_dlnd, d_dlnT, ierr)
    if (ierr /= 0) return

    call set_face_ad_from_log(eos_res(i_lnPgas), d_dlnd(i_lnPgas), d_dlnT(i_lnPgas), T_face, rho_face, Pgas_face)
    P_face = Pgas_face + crad*pow4(T_face)/3d0

    call set_scale_height_from_face_state(s, k, P_face, rho_face, scale_height_face)
  end subroutine get_face_scale_height_ad


  subroutine set_scale_height_from_face_state(s, k, P_face, rho_face, scale_height_face)
    use auto_diff_support, only: wrap_r_00

    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(in) :: P_face, rho_face
    type(auto_diff_real_star_order1), intent(out) :: scale_height_face

    real(dp) :: G
    type(auto_diff_real_star_order1) :: grav, scale_height2

    G = s%cgrav(k)
    grav = G*s%m_grav(k)/pow2(wrap_r_00(s,k))
    scale_height_face = P_face/(grav*rho_face)
    if (s%alt_scale_height_flag) then
       scale_height2 = sqrt(P_face/G)/rho_face
       if (scale_height2 < scale_height_face) scale_height_face = scale_height2
    end if
  end subroutine set_scale_height_from_face_state


  subroutine set_gradr_from_face_state(s, k, P_face, opacity_face, T_face, gradr_face)
    use auto_diff_support, only: wrap_L_00

    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(in) :: P_face, opacity_face, T_face
    type(auto_diff_real_star_order1), intent(out) :: gradr_face

    real(dp) :: L_theta
    type(auto_diff_real_star_order1) :: L_face, Pr_face

    if (s%include_mlt_in_velocity_time_centering) then
       if (s%using_velocity_time_centering .and. &
             s%include_L_in_velocity_time_centering .and. &
             s%lnT(k) <= s%max_logT_for_include_P_and_L_in_velocity_time_centering*ln10) then
          L_theta = s%L_theta_for_velocity_time_centering
       else
          L_theta = 1d0
       end if
       L_face = L_theta*wrap_L_00(s, k) + (1d0 - L_theta)*s%L_start(k)
    else
       L_face = wrap_L_00(s, k)
    end if

    Pr_face = crad*pow4(T_face)/3d0
    gradr_face = P_face*opacity_face*L_face/(4d0*pi4*clight*s%m_grav(k)*s%cgrav(k)*Pr_face)
  end subroutine set_gradr_from_face_state

end module mlt_tdc_face_support
