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
       call ensure_mlt_tdc_face_state_ad(s, k, ierr)
       if (ierr /= 0) return
       T_face = s%mlt_tdc_T_face_ad(k)
       rho_face = s%mlt_tdc_rho_face_ad(k)
       P_face = s%mlt_tdc_P_face_ad(k)
       Cp_face = s%mlt_tdc_Cp_face_ad(k)
       ChiRho_face = s%mlt_tdc_ChiRho_face_ad(k)
       ChiT_face = s%mlt_tdc_ChiT_face_ad(k)
       grada_face = s%mlt_tdc_grada_face_ad(k)
       opacity_face = s%mlt_tdc_opacity_face_ad(k)
       scale_height_face = s%mlt_tdc_scale_height_face_ad(k)
       gradr_face = s%mlt_tdc_gradr_face_ad(k)
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


  ! Ensures that the recomputed MLT and TDC face thermodynamic quantities have
  ! been assembled and cached for face k.
  subroutine ensure_mlt_tdc_face_state_ad(s, k, ierr)
    type(star_info), pointer :: s
    integer, intent(in) :: k
    integer, intent(out) :: ierr

    type(auto_diff_real_star_order1) :: T_face, rho_face, P_face, Cp_face
    type(auto_diff_real_star_order1) :: ChiRho_face, ChiT_face, grada_face, opacity_face
    type(auto_diff_real_star_order1) :: scale_height_face, gradr_face

    ierr = 0
    if (s%have_mlt_tdc_face_state(k)) return

    call build_mlt_tdc_face_state_ad( &
       s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, &
       opacity_face, scale_height_face, gradr_face, ierr)
    if (ierr /= 0) return

    s%mlt_tdc_T_face_ad(k) = T_face
    s%mlt_tdc_rho_face_ad(k) = rho_face
    s%mlt_tdc_P_face_ad(k) = P_face
    s%mlt_tdc_Cp_face_ad(k) = Cp_face
    s%mlt_tdc_ChiRho_face_ad(k) = ChiRho_face
    s%mlt_tdc_ChiT_face_ad(k) = ChiT_face
    s%mlt_tdc_grada_face_ad(k) = grada_face
    s%mlt_tdc_opacity_face_ad(k) = opacity_face
    s%mlt_tdc_scale_height_face_ad(k) = scale_height_face
    s%mlt_tdc_gradr_face_ad(k) = gradr_face
    s%have_mlt_tdc_face_state(k) = .true.
  end subroutine ensure_mlt_tdc_face_state_ad

  ! Reconstructs the face composition from either the current or the
  ! start-of-step composition and renormalizes xa_face.
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
      if (s%report_ierr) then
         !$OMP critical (mlt_tdc_face_report_ierr)
         write(*,*) 'get_face_composition: sum_xa <= 0 for k', k
         !$OMP end critical (mlt_tdc_face_report_ierr)
      end if
      return
    end if
    xa_face = xa_face/sum_xa
  end subroutine get_face_composition


  ! Builds the recomputed face EOS input state by wrapping T and rho to the
  ! face, reconstructing the face composition, and evaluating the EOS there.
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
       if (s%report_ierr) then
          !$OMP critical (mlt_tdc_face_report_ierr)
          write(*,*) 'get_face_eos_inputs: bad face T or rho for k', k, T_face%val, rho_face%val
          !$OMP end critical (mlt_tdc_face_report_ierr)
       end if
       return
    end if

    call get_face_composition(s, k, .false., zbar_face, eos_xa_face, ierr)
    if (ierr /= 0) return

    log10_T = log10(T_face%val)
    log10_rho = log10(rho_face%val)

    call get_eos( &
       s, k, eos_xa_face, rho_face%val, log10_rho, T_face%val, log10_T, &
       eos_res, d_dlnd, d_dlnT, d_dxa, ierr)
    if (ierr /= 0) then
       if (s%report_ierr) call write_face_eos_call_info(s, k, T_face, rho_face, zbar_face, eos_xa_face)
       return
    end if
  end subroutine get_face_eos_inputs


  ! Interpolates extra_opacity_factor to the face and applies the existing
  ! logT taper used to turn that factor on and off.
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


  ! Returns the cached or newly built face EOS and opacity state as
  ! auto_diff_real_star_order1 quantities for the MLT/TDC solve,
  ! instead of using the stored face quantities.
  subroutine get_face_eos_kap_ad( &
       s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, opacity_face, ierr)
    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(out) :: T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, opacity_face
    integer, intent(out) :: ierr

    ierr = 0
    call ensure_mlt_tdc_face_state_ad(s, k, ierr)
    if (ierr /= 0) return

    T_face = s%mlt_tdc_T_face_ad(k)
    rho_face = s%mlt_tdc_rho_face_ad(k)
    P_face = s%mlt_tdc_P_face_ad(k)
    Cp_face = s%mlt_tdc_Cp_face_ad(k)
    ChiRho_face = s%mlt_tdc_ChiRho_face_ad(k)
    ChiT_face = s%mlt_tdc_ChiT_face_ad(k)
    grada_face = s%mlt_tdc_grada_face_ad(k)
    opacity_face = s%mlt_tdc_opacity_face_ad(k)
  end subroutine get_face_eos_kap_ad


  ! Builds the full set of recomputed face thermodynamic quantities for the
  ! MLT and TDC solve from one EOS call and one opacity call at face k.
  subroutine build_mlt_tdc_face_state_ad( &
       s, k, T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, &
       opacity_face, scale_height_face, gradr_face, ierr)
    use eos_def, only: num_eos_basic_results, i_lnPgas, i_grad_ad, i_gamma1, i_Cp, i_chiRho, i_chiT, i_eta, i_lnfree_e
    use kap_def, only: num_kap_fracs

    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(out) :: &
       T_face, rho_face, P_face, Cp_face, ChiRho_face, ChiT_face, grada_face, &
       opacity_face, scale_height_face, gradr_face
    integer, intent(out) :: ierr

    real(dp) :: log10_T, log10_rho, kap_zbar_face, opacity_factor_face
    real(dp) :: eos_res(num_eos_basic_results), d_dlnd(num_eos_basic_results), d_dlnT(num_eos_basic_results)
    real(dp) :: dlnT_face(auto_diff_star_num_vars), dlnd_face(auto_diff_star_num_vars)
    real(dp) :: kap, dlnkap_dlnd, dlnkap_dlnT
    real(dp) :: kap_fracs(num_kap_fracs), kap_xa_face(s%species)
    type(auto_diff_real_star_order1) :: Pgas_face, gamma1_face, mlt_Pturb_ad, alpha

    ierr = 0
    call get_face_eos_inputs( &
       s, k, T_face, rho_face, eos_res, d_dlnd, d_dlnT, ierr)
    if (ierr /= 0) return
    log10_T = log10(T_face%val)
    log10_rho = log10(rho_face%val)
    call set_face_log_partials(T_face, rho_face, dlnT_face, dlnd_face)

    call set_face_ad_from_log(eos_res(i_lnPgas), d_dlnd(i_lnPgas), d_dlnT(i_lnPgas), dlnd_face, dlnT_face, Pgas_face)
    P_face = Pgas_face + crad*pow4(T_face)/3d0
    call set_face_ad_from_value(eos_res(i_Cp), d_dlnd(i_Cp), d_dlnT(i_Cp), dlnd_face, dlnT_face, Cp_face)
    call set_face_ad_from_value(eos_res(i_chiRho), d_dlnd(i_chiRho), d_dlnT(i_chiRho), dlnd_face, dlnT_face, ChiRho_face)
    call set_face_ad_from_value(eos_res(i_chiT), d_dlnd(i_chiT), d_dlnT(i_chiT), dlnd_face, dlnT_face, ChiT_face)
    call set_face_ad_from_value(eos_res(i_grad_ad), d_dlnd(i_grad_ad), d_dlnT(i_grad_ad), dlnd_face, dlnT_face, grada_face)
    call set_face_ad_from_value(eos_res(i_gamma1), d_dlnd(i_gamma1), d_dlnT(i_gamma1), dlnd_face, dlnT_face, gamma1_face)

    call get_face_composition(s, k, s% use_starting_composition_for_kap, kap_zbar_face, kap_xa_face, ierr)
    if (ierr /= 0) return
    call get_face_opacity_factor(s, k, log10_T, opacity_factor_face)

    call get_kap( &
       s, k, kap_zbar_face, kap_xa_face, log10_rho, log10_T, &
       eos_res(i_lnfree_e), d_dlnd(i_lnfree_e), d_dlnT(i_lnfree_e), &
       eos_res(i_eta), d_dlnd(i_eta), d_dlnT(i_eta), &
       kap_fracs, kap, dlnkap_dlnd, dlnkap_dlnT, ierr)
    if (ierr /= 0) then
       if (s%report_ierr) call write_face_kap_call_info( &
          s, k, T_face, rho_face, kap_zbar_face, kap_xa_face, opacity_factor_face)
       return
    end if
    if (is_bad_num(kap) .or. kap <= 0d0) then
       ierr = -1
       if (s%report_ierr) then
          !$OMP critical (mlt_tdc_face_report_ierr)
          write(*,*) 'get_face_eos_kap_ad: bad face opacity for k', k, kap
          !$OMP end critical (mlt_tdc_face_report_ierr)
          call write_face_kap_call_info(s, k, T_face, rho_face, kap_zbar_face, kap_xa_face, opacity_factor_face)
       end if
       return
    end if

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
    call set_face_ad_from_value(kap, kap*dlnkap_dlnd, kap*dlnkap_dlnT, dlnd_face, dlnT_face, opacity_face)

    if (s% have_mlt_vc .and. s% okay_to_set_mlt_vc .and. s% include_mlt_Pturb_in_thermodynamic_gradients &
       .and. s% mlt_Pturb_factor > 0d0 .and. k > 1) then
       mlt_Pturb_ad = s% mlt_Pturb_factor*pow2(s% mlt_vc_old(k))*rho_face/3d0
       alpha = mlt_Pturb_ad/(P_face*gamma1_face)
       grada_face = grada_face*(P_face + mlt_Pturb_ad)/(P_face*(1d0 + alpha))
    end if

    call set_scale_height_from_face_state(s, k, P_face, rho_face, scale_height_face)
    call set_gradr_from_face_state(s, k, P_face, opacity_face, T_face, gradr_face)
  end subroutine build_mlt_tdc_face_state_ad


  ! Returns the cached or newly built face pressure scale height as an
  ! auto_diff_real_star_order1 quantity from the face EOS state.
  subroutine get_face_scale_height_ad(s, k, scale_height_face, ierr)
    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(out) :: scale_height_face
    integer, intent(out) :: ierr

    ierr = 0
    call ensure_mlt_tdc_face_state_ad(s, k, ierr)
    if (ierr /= 0) return
    scale_height_face = s%mlt_tdc_scale_height_face_ad(k)
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

  ! Precomputes dlnT_face and dlnd_face for converting scalar d/dlnT and
  ! d/dlnd microphysics partials into star-order1 autodiff derivatives.
  subroutine set_face_log_partials(T_face, rho_face, dlnT_face, dlnd_face)
    type(auto_diff_real_star_order1), intent(in) :: T_face, rho_face
    real(dp), intent(out) :: dlnT_face(auto_diff_star_num_vars), dlnd_face(auto_diff_star_num_vars)

    dlnT_face = T_face%d1Array/T_face%val
    dlnd_face = rho_face%d1Array/rho_face%val
  end subroutine set_face_log_partials


  ! Converts a scalar value with d/dlnd and d/dlnT partials into an
  ! auto_diff_real_star_order1 quantity using the face chain rule.
  subroutine set_face_ad_from_value(value, dvalue_dlnd, dvalue_dlnT, dlnd_face, dlnT_face, quantity_ad)
    real(dp), intent(in) :: value, dvalue_dlnd, dvalue_dlnT
    real(dp), intent(in) :: dlnd_face(auto_diff_star_num_vars), dlnT_face(auto_diff_star_num_vars)
    type(auto_diff_real_star_order1), intent(out) :: quantity_ad

    quantity_ad = 0d0
    quantity_ad%val = value
    quantity_ad%d1Array = dvalue_dlnd*dlnd_face + dvalue_dlnT*dlnT_face
  end subroutine set_face_ad_from_value


  ! Same as set_face_ad_from_value, but for a quantity returned in
  ! logarithmic form by the microphysics routine.
  subroutine set_face_ad_from_log(log_value, dlog_dlnd, dlog_dlnT, dlnd_face, dlnT_face, quantity_ad)
    real(dp), intent(in) :: log_value, dlog_dlnd, dlog_dlnT
    real(dp), intent(in) :: dlnd_face(auto_diff_star_num_vars), dlnT_face(auto_diff_star_num_vars)
    type(auto_diff_real_star_order1), intent(out) :: quantity_ad
    real(dp) :: value

    value = exp(log_value)
    call set_face_ad_from_value(value, value*dlog_dlnd, value*dlog_dlnT, dlnd_face, dlnT_face, quantity_ad)
  end subroutine set_face_ad_from_log


  ! Writes the recomputed face EOS inputs that were passed to get_eos.
  subroutine write_face_eos_call_info(s, k, T_face, rho_face, zbar_face, xa_face)
    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(in) :: T_face, rho_face
    real(dp), intent(in) :: zbar_face
    real(dp), intent(in) :: xa_face(:)

    integer :: j
    include 'formats'

    !$OMP critical (mlt_tdc_face_eos_call_info)
    write(*,'(A)')
    write(*,*) 'face EOS input info for k', k
    write(*,1) 'T_face', T_face%val
    write(*,1) 'rho_face', rho_face%val
    write(*,1) 'log10_T_face', log10(T_face%val)
    write(*,1) 'log10_rho_face', log10(rho_face%val)
    write(*,1) 'zbar_face', zbar_face
    write(*,1) 'sum(xa_face)', sum(xa_face)
    do j = 1, s%species
       write(*,2) 'xa_face ' // trim(s%nameofequ(j+s%nvar_hydro)), j, xa_face(j)
    end do
    !$OMP end critical (mlt_tdc_face_eos_call_info)
  end subroutine write_face_eos_call_info


  ! Writes the recomputed face opacity inputs that were passed to get_kap.
  subroutine write_face_kap_call_info(s, k, T_face, rho_face, zbar_face, xa_face, opacity_factor_face)
    type(star_info), pointer :: s
    integer, intent(in) :: k
    type(auto_diff_real_star_order1), intent(in) :: T_face, rho_face
    real(dp), intent(in) :: zbar_face
    real(dp), intent(in) :: xa_face(:)
    real(dp), intent(in) :: opacity_factor_face

    integer :: j
    include 'formats'

    !$OMP critical (mlt_tdc_face_kap_call_info)
    write(*,'(A)')
    write(*,*) 'face opacity input info for k', k
    write(*,1) 'T_face', T_face%val
    write(*,1) 'rho_face', rho_face%val
    write(*,1) 'log10_T_face', log10(T_face%val)
    write(*,1) 'log10_rho_face', log10(rho_face%val)
    write(*,1) 'zbar_face', zbar_face
    write(*,1) 'opacity_factor_face', opacity_factor_face
    write(*,1) 'sum(xa_face)', sum(xa_face)
    do j = 1, s%species
       write(*,2) 'xa_face ' // trim(s%nameofequ(j+s%nvar_hydro)), j, xa_face(j)
    end do
    !$OMP end critical (mlt_tdc_face_kap_call_info)
  end subroutine write_face_kap_call_info

end module mlt_tdc_face_support
