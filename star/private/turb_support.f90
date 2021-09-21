! ***********************************************************************
!
!   Copyright (C) 2010-2021  The MESA Team
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


module turb_support

use star_private_def
use const_def
use num_lib
use utils_lib
use auto_diff_support
use star_utils
use turb

implicit none

private
public :: get_gradT, do1_mlt_eval, Get_results

contains

   !> Determines if it is safe (physically) to use TDC instead of MLT.
   !!
   !! Currently we only know we have to fall back to MLT in cells that get touched
   !! by adjust_mass, because there the convection speeds at the start of the
   !! step can be badly out of whack.
   !!
   !! @param s star pointer
   !! @param k face index
   !! @param fallback False if we can use TDC, True if we can fall back to MLT.
   logical function check_if_must_fall_back_to_MLT(s, k) result(fallback)
      type (star_info), pointer :: s
      integer, intent(in) :: k

      fallback = .false.
      if (abs(s%mstar_dot) > 1d-99 .and. k < s% k_const_mass) then
         fallback = .true.
      end if
   end function check_if_must_fall_back_to_MLT

   subroutine get_gradT(s, MLT_option, & ! used to create models
         r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
         iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
         mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
      type (star_info), pointer :: s
      character (len=*), intent(in) :: MLT_option
      real(dp), intent(in) :: &
         r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
         XH1, cgrav, m, gradL_composition_term, mixing_length_alpha
      integer, intent(in) :: iso
      real(dp), intent(out) :: gradT, Y_face, conv_vel, D, Gamma
      integer, intent(out) :: mixing_type, ierr 
      type(auto_diff_real_star_order1) :: &
         gradr_ad, grada_ad, scale_height_ad, gradT_ad, Y_face_ad, mlt_vc_ad, D_ad, &
         Gamma_ad, r_ad, L_ad, T_ad, P_ad, opacity_ad, rho_ad, dV_ad, chiRho_ad, chiT_ad, Cp_ad
      ierr = 0
      r_ad = r
      L_ad = L
      T_ad = T
      P_ad = P
      opacity_ad = opacity
      rho_ad = rho
      dV_ad = 0d0
      chiRho_ad = chiRho
      chiT_ad = chiT
      Cp_ad = Cp
      gradr_ad = gradr
      grada_ad = grada
      scale_height_ad = scale_height
      call Get_results(s, 0, MLT_option, &
         r_ad, L_ad, T_ad, P_ad, opacity_ad, rho_ad, dV_ad, chiRho_ad, chiT_ad, Cp_ad, &
         gradr_ad, grada_ad, scale_height_ad, &
         iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
         s% alpha_semiconvection, s% thermohaline_coeff, &
         mixing_type, gradT_ad, Y_face_ad, mlt_vc_ad, D_ad, Gamma_ad, ierr)
      gradT = gradT_ad%val
      Y_face = Y_face_ad%val
      conv_vel = mlt_vc_ad%val
      D = D_ad%val
      Gamma = Gamma_ad%val
   end subroutine get_gradT
   
      
   subroutine do1_mlt_eval( &
         s, k, MLT_option, gradL_composition_term, &
         gradr, grada, scale_height, mixing_length_alpha, &
         mixing_type, gradT, Y_face, mlt_vc, D, Gamma, ierr)
      use chem_def, only: ih1
      type (star_info), pointer :: s
      integer, intent(in) :: k
      character (len=*), intent(in) :: MLT_option
      type(auto_diff_real_star_order1), intent(in) :: gradr, grada, scale_height
      real(dp), intent(in) :: gradL_composition_term, mixing_length_alpha
      integer, intent(out) :: mixing_type
      type(auto_diff_real_star_order1), intent(out) :: &
         gradT, Y_face, mlt_vc, D, Gamma
      integer, intent(out) :: ierr 
              
      real(dp) :: cgrav, m, XH1, gradL_old, grada_face_old
      integer :: iso, old_mix_type
      type(auto_diff_real_star_order1) :: r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp
      include 'formats'
      ierr = 0
      
      cgrav = s% cgrav(k)
      m = s% m_grav(k)
      L = wrap_L_00(s,k)
      T = get_T_face(s,k)
      P = get_Peos_face(s,k)
      r = wrap_r_00(s,k)
      opacity = get_kap_face(s,k)
      rho = get_Rho_face(s,k)
      dV = 1d0/rho - 1d0/s% rho_start(k)
      chiRho = get_ChiRho_face(s,k)
      chiT = get_ChiT_face(s,k)
      Cp = get_Cp_face(s,k)
      iso = s% dominant_iso_for_thermohaline(k)
      XH1 = s% xa(s% net_iso(ih1),k)
      
      if (s% use_other_mlt_results) then
         call s% other_mlt_results(s% id, k, MLT_option, &
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            s% alpha_semiconvection, s% thermohaline_coeff, &
            mixing_type, gradT, Y_face, mlt_vc, D, Gamma, ierr)
      else         
         call Get_results(s, k, MLT_option, &
            r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            s% alpha_semiconvection, s% thermohaline_coeff, &
            mixing_type, gradT, Y_face, mlt_vc, D, Gamma, ierr)
      end if

   end subroutine do1_mlt_eval


   subroutine Get_results(s, k, MLT_option, &  ! NOTE: k=0 is a valid arg
         r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height, &
         iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
         alpha_semiconvection, thermohaline_coeff, &
         mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
      use star_utils
      type (star_info), pointer :: s
      integer, intent(in) :: k
      character (len=*), intent(in) :: MLT_option
      type(auto_diff_real_star_order1), intent(in) :: &
         r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height
      integer, intent(in) :: iso
      real(dp), intent(in) :: &
         XH1, cgrav, m, gradL_composition_term, &
         mixing_length_alpha, alpha_semiconvection, thermohaline_coeff
      integer, intent(out) :: mixing_type
      type(auto_diff_real_star_order1), intent(out) :: gradT, Y_face, conv_vel, D, Gamma
      integer, intent(out) :: ierr
      
      type(auto_diff_real_star_order1) :: Pr, Pg, grav, Lambda, gradL, beta
      real(dp) :: conv_vel_start, scale
      character (len=256) :: message        
      logical ::  test_partials, using_TDC
      logical, parameter :: report = .false.
      include 'formats'

      ! Pre-calculate some things. 
      Pr = crad*pow4(T)/3d0
      Pg = P - Pr
      beta = Pg / P
      Lambda = mixing_length_alpha*scale_height
      grav = cgrav*m/pow2(r)   
      if (s% use_Ledoux_criterion) then
         gradL = grada + gradL_composition_term ! Ledoux temperature gradient
      else
         gradL = grada
      end if

      ! Initialize with no mixing
      mixing_type = no_mixing
      gradT = gradr
      Y_face = gradT - gradL
      conv_vel = 0d0
      D = 0d0
      Gamma = 0d0  

      ! Bail if we asked for no mixing, or if parameters are bad.
      if (MLT_option == 'none' .or. beta < 1d-10 .or. mixing_length_alpha <= 0d0 .or. &
            opacity%val < 1d-10 .or. P%val < 1d-20 .or. T%val < 1d-10 .or. Rho%val < 1d-20 &
            .or. m < 1d-10 .or. r%val < 1d-10 .or. cgrav < 1d-10) return

      !test_partials = (k == s% solver_test_partials_k)
      test_partials = .false.
      ierr = 0          
      if (k > 0) then
         s% tdc_num_iters(k) = 0
      end if

      if (report) then
         write(*,*)
         write(*,4) 'enter Get_results k slvr_itr model gradr grada scale_height ' // trim(MLT_option), &
            k, s% solver_iter, s% model_number, gradr%val, grada%val, scale_height%val
      end if



      ! check if this particular k can be done with TDC
      using_TDC = .false.
      if (s% MLT_option == 'TDC') using_TDC = .true.
      if (.not. s% have_mlt_vc) using_TDC = .false.
      if (k <= 0 .or. s%dt <= 0d0) using_TDC = .false.
      if (using_TDC) using_TDC = .not. check_if_must_fall_back_to_MLT(s, k)

      if (using_TDC) then
         if (report) write(*,3) 'call set_TDC', k, s% solver_iter
         if (s% okay_to_set_mlt_vc) then
            conv_vel_start = s% mlt_vc_old(k)
         else
            conv_vel_start = s% mlt_vc(k)
         end if

         ! Set scale for judging the TDC luminosity equation Q(Y)=0.
         ! Q has units of a luminosity, so the scale should be a luminosity.
         if (s% solver_iter == 0) then
            scale = max(abs(s% L(k)), 1d-3*maxval(s% L(1:s% nz)))
         else
            scale = max(abs(s% L_start(k)), 1d-3*maxval(s% L_start(1:s% nz)))
         end if

         call set_TDC(&
            conv_vel_start, mixing_length_alpha, s% alpha_TDC_DAMP, s%alpha_TDC_DAMPR, s%alpha_TDC_PtdVdt, s%dt, cgrav, m, report, &
            mixing_type, scale, L, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, s%tdc_num_iters(k), ierr)
      else if (gradr > gradL) then
         if (report) write(*,3) 'call set_MLT', k, s% solver_iter
         call set_MLT(MLT_option, mixing_length_alpha, report, s% Henyey_MLT_nu_param, s% Henyey_MLT_y_param, &
                        chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                        gradr, grada, gradL, k, &
                        Gamma, gradT, Y_face, conv_vel, D, mixing_type, ierr)
      end if

      ! If we're not convecting, try thermohaline and semiconvection.
      if (mixing_type == no_mixing) then
         if (gradL_composition_term < 0) then
            if (report) write(*,3) 'call set_thermohaline', k, s% solver_iter
            call set_thermohaline(s%thermohaline_option, Lambda, grada, gradr, T, opacity, rho, Cp, gradL_composition_term, &
                              iso, XH1, thermohaline_coeff, &
                              D, gradT, Y_face, conv_vel, mixing_type, ierr)
         else if (gradr > grada) then
            if (report) write(*,3) 'call set_semiconvection', k, s% solver_iter
            call set_semiconvection(L, Lambda, m, T, P, Pr, beta, opacity, rho, alpha_semiconvection, &
                                    s% semiconvection_option, cgrav, Cp, gradr, grada, gradL, &
                                    gradL_composition_term, &
                                    gradT, Y_face, conv_vel, D, mixing_type, ierr)
         end if         
      end if 

      ! If there's too-little mixing to bother, or we hit a bad value, fall back on no mixing.
      if (D%val < s% remove_small_D_limit .or. is_bad(D%val)) then
         if (report) write(*,2) 'D < s% remove_small_D_limit', k, D%val, s% remove_small_D_limit
         mixing_type = no_mixing
         gradT = gradr
         Y_face = gradT - gradL
         conv_vel = 0d0
         D = 0d0
         Gamma = 0d0            
      end if
   end subroutine Get_results


end module turb_support
