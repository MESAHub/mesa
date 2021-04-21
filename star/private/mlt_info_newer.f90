! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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


      module mlt_info_newer

      use star_private_def
      use const_def
      use num_lib
      use utils_lib
      use mlt_get_results_newer, only: do1_mlt_eval_newer
      use auto_diff_support, only: wrap

      implicit none

      private
      public :: &
         set_mlt_vars_newer, &
         do1_mlt_2_newer, &
         set_grads_newer, &
         switch_to_radiative_newer, &
         check_for_redo_MLT_newer, &
         set_gradT_excess_alpha_newer
         

      logical, parameter :: dbg = .false.
      integer, parameter :: kdbg = -1
      
      integer, parameter :: nvbs = num_mlt_partials

      contains


      subroutine set_mlt_vars_newer(s, nzlo, nzhi, ierr)
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: k, op_err
         integer(8) :: time0
         real(dp) :: total, opacity, gamma1, Cv, chiRho, chiT, Cp, &
            grada, P, xh, gradL_composition_term
         logical :: make_gradr_sticky_in_solver_iters
         include 'formats'
         ierr = 0
         if (s% x_integer_ctrl(19) > 0) write(*,3) 'doing set_mlt_vars_newer solver model iter', s% model_number, s% solver_iter
         gradL_composition_term = -1d0
         opacity = -1d0
         chiRho = -1d0
         chiT = -1d0
         Cp = -1d0
         grada = -1d0
         P = -1d0
         xh = -1d0 
         if (s% doing_timing) call start_time(s, time0, total)
         if (s% x_integer_ctrl(19) > 0) write(*,*) 'start set_mlt_vars_newer'
!$OMP PARALLEL DO PRIVATE(k,op_err,make_gradr_sticky_in_solver_iters) SCHEDULE(dynamic,2)
         do k = nzlo, nzhi
            op_err = 0
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,3) 'set_mlt_vars_newer call do1_mlt_2_newer'
            end if
            call do1_mlt_2_newer(s, k, s% alpha_mlt(k), gradL_composition_term, &
               opacity, chiRho, chiT, Cp, grada, P, xh, &
               .false., make_gradr_sticky_in_solver_iters, op_err)
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,3) 'set_mlt_vars_newer done do1_mlt_2_newer', k
            end if
            if (make_gradr_sticky_in_solver_iters .and. s% solver_iter > 3) then
               if (.not. s% fixed_gradr_for_rest_of_solver_iters(k)) then
                  s% fixed_gradr_for_rest_of_solver_iters(k) = &
                     (s% mlt_mixing_type(k) == no_mixing)
               end if
            end if            
         end do
!$OMP END PARALLEL DO
         if (s% x_integer_ctrl(19) > 0) write(*,*) 'done set_mlt_vars_newer'
         if (s% x_integer_ctrl(19) > 0) write(*,*)
         if (s% doing_timing) call update_time(s, time0, total, s% time_mlt)
      end subroutine set_mlt_vars_newer


      subroutine do1_mlt_2_newer(s, k, mixing_length_alpha, gradL_composition_term_in, &
            opacity_face_in, chiRho_face_in, &
            chiT_face_in, Cp_face_in, grada_face_in, P_face_in, xh_face_in, &
            from_do1_mlt, make_gradr_sticky_in_solver_iters, ierr)
         ! get convection info for point k
         use star_utils
         use eos_def
         use chem_def, only: ih1
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha, &
            gradL_composition_term_in, opacity_face_in, &
            chiRho_face_in, chiT_face_in, &
            Cp_face_in, grada_face_in, P_face_in, xh_face_in
         logical, intent(in) :: from_do1_mlt
         logical, intent(out) :: make_gradr_sticky_in_solver_iters
         integer, intent(out) :: ierr

         real(dp) :: max_conv_vel, alfa, beta, v, gradr_factor, d_gradr_factor_dw, &
            f, xh_face, gradL_composition_term, abs_du_div_cs, cs
         real(dp), pointer :: vel(:)
         integer :: i, mixing_type, h1, nz, k_T_max
         real(dp), parameter :: conv_vel_mach_limit = 0.9d0
         character (len=32) :: MLT_option
         logical, parameter :: just_gradr = .false.

         type(auto_diff_real_star_order1) :: &
            opacity_face_ad, chiRho_face_ad, chiT_face_ad, Cp_face_ad, &
            grada_face_ad, P_face_ad, T_face_ad, rho_face_ad, gradr_ad, &
            gradT_ad, Y_face_ad, mlt_vc_ad, Gamma_ad, D_ad, scale_height_ad

         include 'formats'

         ierr = 0
         nz = s% nz

            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,3) 'start do1_mlt_2_newer k iter', k, s% solver_iter
            end if
         
         if (k < 1 .or. k > nz) then
            write(*,3) 'bad k for do1_mlt', k, nz
            ierr = -1
            return
            call mesa_error(__FILE__,__LINE__)
         end if
         
         MLT_option = s% MLT_option

         ! alfa is the fraction coming from k; (1-alfa) from k-1.
         if (k == 1) then
            alfa = 1d0
         else
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         end if
         beta = 1d0 - alfa
         h1 = s% net_iso(ih1)

         gradL_composition_term = gradL_composition_term_in      

         opacity_face_ad = opacity_face_in
         if (opacity_face_in < 0) opacity_face_ad = get_kap_face(s,k)

         chiRho_face_ad = chiRho_face_in
         if (chiRho_face_in < 0) chiRho_face_ad = get_ChiRho_face(s,k)

         chiT_face_ad = chiT_face_in
         if (chiT_face_in < 0) chiT_face_ad = get_ChiT_face(s,k)

         Cp_face_ad = Cp_face_in
         if (Cp_face_in < 0) Cp_face_ad = get_Cp_face(s,k)

         grada_face_ad = grada_face_in
         if (grada_face_in < 0) grada_face_ad = get_grada_face(s,k)
         s% grada_face(k) = grada_face_ad%val         

         P_face_ad = P_face_in
         if (P_face_in < 0) P_face_ad = get_Peos_face(s,k)

         T_face_ad = get_T_face(s,k)
   
         rho_face_ad = get_rho_face(s,k)

         scale_height_ad = get_scale_height_face(s,k)

         xh_face = xh_face_in   
         if (h1 /= 0 .and. xh_face < 0) xh_face = s% xa(h1, k)

         if (s% rotation_flag .and. s% mlt_use_rotation_correction) then
            gradr_factor = s% ft_rot(k)/s% fp_rot(k)*s% gradr_factor(k)
            if (s% w_div_wc_flag) then
               d_gradr_factor_dw = gradr_factor*(s% dft_rot_dw_div_wc(k)/s%ft_rot(k) &
                  -s% dfp_rot_dw_div_wc(k)/s%fp_rot(k) )
            end if
         else
            gradr_factor = s% gradr_factor(k)
            d_gradr_factor_dw = 0d0
         end if
         if (is_bad_num(gradr_factor)) then
            ierr = -1
            return
         end if
         gradr_ad = get_gradr_face(s,k)*gradr_factor
         
         if (k == 1 .and. s% mlt_make_surface_no_mixing) then
            call set_no_mixing('surface_no_mixing')
            return
         end if
         
         if (s% lnT_start(k)/ln10 > s% max_logT_for_mlt) then
            call set_no_mixing('max_logT')
            return
         end if
         
         if (s% no_MLT_below_shock .and. (s%u_flag .or. s%v_flag)) then ! check for outward shock above k
            if (s% u_flag) then
               vel => s% u
            else
               vel => s% v
            end if
            do i=k-1,1,-1
               cs = s% csound(i)
               if (vel(i+1) >= cs .and. vel(i) < cs) then
                  call set_no_mixing('below_shock')
                  return
               end if
            end do
         end if
         
         if (s% no_MLT_below_T_max) then
            k_T_max = maxloc(s% T_start(1:nz),dim=1)
            if (k > k_T_max) then
               call set_no_mixing('below_T_max')
               return
            end if
         end if
         
         make_gradr_sticky_in_solver_iters = s% make_gradr_sticky_in_solver_iters
         if (.not. make_gradr_sticky_in_solver_iters .and. &
               s% min_logT_for_make_gradr_sticky_in_solver_iters < 1d20) then
            k_T_max = maxloc(s% lnT_start(1:nz),dim=1)
            make_gradr_sticky_in_solver_iters = &
               (s% lnT_start(k_T_max)/ln10 >= s% min_logT_for_make_gradr_sticky_in_solver_iters)
         end if
         
         if (make_gradr_sticky_in_solver_iters) then
            if (s% fixed_gradr_for_rest_of_solver_iters(k)) then
               call set_no_mixing('gradr_sticky')
               return
            end if
         end if

         max_conv_vel = s% csound_face(k)*s% max_conv_vel_div_csound
         if (s% dt < s% min_dt_for_increases_in_convection_velocity) then
            max_conv_vel = 1d-2*s% dt*s% cgrav(k)
         end if

         if (s% csound_start(k) > 0d0) then
            if (s% u_flag) then        
               abs_du_div_cs = 0d0
               if (s% u_start(k)/1d5 > s% max_v_for_convection) then
                  if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
                     write(*,2) 'u_start(k)/1d5 > s% max_v_for_convection', k, s% u_start(k)/1d5
                  end if
                  max_conv_vel = 0d0              
               else if (s% q(k) > s% max_q_for_convection_with_hydro_on) then
                  if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
                     write(*,2) 's% q(k) > s% max_q_for_convection_with_hydro_on', k, s% q(k)
                  end if
                  max_conv_vel = 0d0
               else if ((abs(s% u_start(k))) >= &
                     s% csound_start(k)*s% max_v_div_cs_for_convection) then
                  if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
                     write(*,2) 'abs(s% u_start(k)))', k, abs(s% u_start(k))/1d5
                  end if
                  max_conv_vel = 0d0              
               else
                  if (k == 1) then
                     abs_du_div_cs = 1d99
                  else if (k < nz) then
                     abs_du_div_cs = max(abs(s% u_start(k) - s% u_start(k+1)), &
                         abs(s% u_start(k) - s% u_start(k-1))) / s% csound_start(k)
                  end if
                  if (abs_du_div_cs > s% max_abs_du_div_cs_for_convection) then
                     if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
                        write(*,2) 'max_v_div_cs_for_convection', k, s% max_v_div_cs_for_convection
                     end if
                     max_conv_vel = 0d0
                  end if
               end if
            else if (s% v_flag) then
               if (s% v_start(k)/1d5 > s% max_v_for_convection) then
                  if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
                     write(*,2) 's% v_start(k)/1d5 > s% max_v_for_convection', k, s% v_start(k)/1d5
                  end if
                  max_conv_vel = 0d0              
               else if (s% q(k) > s% max_q_for_convection_with_hydro_on) then
                  if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
                     write(*,2) 's% q(k) > s% max_q_for_convection_with_hydro_on', k, s% q(k)
                  end if
                  max_conv_vel = 0d0
               else if ((abs(s% v_start(k))) >= &
                     s% csound_start(k)*s% max_v_div_cs_for_convection) then
                  if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
                     write(*,2) 'max_v_div_cs_for_convection', k, s% max_v_div_cs_for_convection
                  end if
                  max_conv_vel = 0d0
               end if
            end if
            if (max_conv_vel == 0d0) then
               if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
                  write(*,2) 'max_conv_vel == 0d0', k, max_conv_vel
               end if
               MLT_option = 'none'
            end if
         end if

         if (s% use_Ledoux_criterion .and. gradL_composition_term < 0) then
            gradL_composition_term = s% gradL_composition_term(k)
         else
            gradL_composition_term = 0d0
         end if

            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'before call do1_mlt_eval_newer gradr_factor comp_term ' // trim(mlt_option), k, gradr_factor, gradL_composition_term
            end if
            
         call do1_mlt_eval_newer(s, k, MLT_option, just_gradr, gradL_composition_term, &
            gradr_ad, grada_face_ad, scale_height_ad, mixing_length_alpha, &
            mixing_type, gradT_ad, Y_face_ad, mlt_vc_ad, D_ad, Gamma_ad, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*,*) 'ierr in do1_mlt_eval_newer for k', k
            end if
            return
         end if
         

            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'after call do1_mlt_eval_newer gradT d_dlnd_00', k, gradT_ad%val, gradT_ad%d1Array(i_lnd_00)
            end if

         s% mlt_mixing_type(k) = mixing_type
         s% gradT_ad(k) = gradT_ad
         s% gradT(k) = s% gradT_ad(k)%val
         s% Y_face_ad(k) = Y_face_ad
         s% Y_face(k) = s% Y_face_ad(k)%val
         s% mlt_vc_ad(k) = mlt_vc_ad
         if (s% okay_to_set_mlt_vc) s% mlt_vc(k) = s% mlt_vc_ad(k)%val  
         s% mlt_D_ad(k) = D_ad
         s% mlt_D(k) = D_ad%val
         s% mlt_Gamma_ad(k) = Gamma_ad
         s% mlt_Gamma(k) = Gamma_ad%val
         s% gradr_ad(k) = gradr_ad
         s% gradr(k) = s% gradr_ad(k)%val
         s% grada_face_ad(k) = grada_face_ad
         s% grada_face(k) = grada_face_ad%val
         s% scale_height_ad(k) = scale_height_ad
         s% scale_height(k) = scale_height_ad%val             
         s% Lambda_ad(k) = mixing_length_alpha*scale_height_ad
         s% mlt_mixing_length(k) = mixing_length_alpha*s% scale_height(k)    
         s% gradL_ad(k) = grada_face_ad + gradL_composition_term
         s% gradL(k) = s% gradL_ad(k)%val
         
         s% mlt_cdc(k) = s% mlt_D(k)*pow2(pi4*pow2(s%r(k))*rho_face_ad%val)
         
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'do1_mlt_eval_newer before adjust_gradT_fraction gradT d_dlnd_00', k, gradT_ad%val, gradT_ad%d1Array(i_lnd_00)
            end if

         if (s% mlt_gradT_fraction >= 0d0 .and. s% mlt_gradT_fraction <= 1d0) then
            f = s% mlt_gradT_fraction
         else
            f = s% adjust_mlt_gradT_fraction(k)
         end if
         call adjust_gradT_fraction(s, k, f)
         
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'do1_mlt_eval_newer after adjust_gradT_fraction gradT d_dlnd_00', k, gradT_ad%val, gradT_ad%d1Array(i_lnd_00)
            end if

         ! Note that L_conv can be negative when conv_vel /= 0 in a convectively
         ! stable region. This can happen when using conv_vel as a variable
         ! Note make_brown_dwarf can have L=0 which means gradr = 0
         if (mlt_option /= 'none' .and. abs(s% gradr(k))>0d0) then
            s% L_conv(k) = s% L(k) * (1d0 - s% gradT(k)/s% gradr(k)) ! C&G 14.109
         else
            s% L_conv(k) = 0d0
         end if

            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'done do1_mlt_2_newer gradT d_dlnd_00', k, gradT_ad%val, gradT_ad%d1Array(i_lnd_00)
            end if

         contains

         subroutine set_no_mixing(str)
            character (len=*) :: str
            include 'formats'
            s% mlt_mixing_type(k) = no_mixing
            !gradr_ad = gradr_factor*gradr_ad
            gradT_ad = gradr_ad
            s% gradT_ad(k) = gradT_ad
            s% gradT(k) = s% gradT_ad(k)%val
            s% mlt_vc_ad(k) = 0d0
            if (s% okay_to_set_mlt_vc) s% mlt_vc(k) = 0d0
            s% mlt_D_ad(k) = 0d0
            s% mlt_D(k) = 0d0
            s% mlt_Gamma_ad(k) = 0d0
            s% mlt_Gamma(k) = 0d0
            s% gradr_ad(k) = gradr_ad
            s% gradr(k) = s% gradr_ad(k)%val
            s% grada_face_ad(k) = grada_face_ad
            s% grada_face(k) = grada_face_ad%val
            s% scale_height_ad(k) = scale_height_ad
            s% scale_height(k) = scale_height_ad%val             
            s% Lambda_ad(k) = mixing_length_alpha*scale_height_ad
            s% mlt_mixing_length(k) = mixing_length_alpha*s% scale_height(k)    
            s% gradL_ad(k) = 0d0
            s% gradL(k) = 0d0
            Y_face_ad = gradT_ad - grada_face_ad
            s% Y_face_ad(k) = Y_face_ad
            s% Y_face(k) = s% Y_face_ad(k)%val
            s% mlt_cdc(k) = 0d0
            s% actual_gradT(k) = 0d0
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,2) 'do1_mlt_2_newer set_no_mixing gradT ' // trim(str), k, s% gradT(k)
            end if
         end subroutine set_no_mixing

      end subroutine do1_mlt_2_newer


      logical function must_limit_conv_vel(s,k0)
         type (star_info), pointer :: s
         integer, intent(in) :: k0
         real(dp) :: alfa, beta, one_m_f
         integer :: k, wdth
         include 'formats'
         must_limit_conv_vel = .false.
         if (s% q(k0) <= s% max_q_for_limit_conv_vel .or. &
             s% m(k0) <= s% max_mass_in_gm_for_limit_conv_vel .or. &
             s% r(k0) <= s% max_r_in_cm_for_limit_conv_vel) then
            must_limit_conv_vel = .true.
            return
         end if
         if (.not. s% v_flag) return
         wdth = s% width_for_limit_conv_vel
         if (wdth < 0) return
         do k = max(2,k0-wdth),min(s% nz,k0+wdth)
            if (s% csound(k) < s% v(k) .and. s% v(k) <= s% csound(k-1)) then
               must_limit_conv_vel = .true.
               return
            end if
         end do
      end function must_limit_conv_vel


      subroutine adjust_gradT_fraction(s,k,f)
         ! replace gradT by combo of grada_face and gradr
         ! then check excess
         use eos_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: f
         integer, intent(in) :: k

         real(dp) :: alfa, beta, one_m_f

         include 'formats'
         
         if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
            write(*,2) 'adjust_gradT_fraction f', k, f
         end if

         if (f >= 0.0 .and. f <= 1.0) then

            ! alfa is the fraction coming from k; (1-alfa) from k-1.
            if (k == 1) then
               alfa = 1.0d0
            else
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            end if
            beta = 1.0d0 - alfa

            if (f == 0d0) then
               s% gradT(k) = s% gradr(k)
               s% d_gradT_dlnR(k) = s% d_gradr_dlnR(k)
               s% d_gradT_dL(k) = s% d_gradr_dL(k)
               s% d_gradT_dw_div_wc(k) = s% d_gradr_dw_div_wc(k)
               s% d_gradT_dlnd00(k) = s% d_gradr_dlnd00(k)
               s% d_gradT_dlnT00(k) = s% d_gradr_dlnT00(k)
               if (k > 1) then
                  s% d_gradT_dlndm1(k) = s% d_gradr_dlndm1(k)
                  s% d_gradT_dlnTm1(k) = s% d_gradr_dlnTm1(k)
               end if
            else ! mix
               one_m_f = 1.0d0 - f
               s% gradT(k) = f*s% grada_face(k) + one_m_f*s% gradr(k)
               s% d_gradT_dlnR(k) = one_m_f*s% d_gradr_dlnR(k)
               s% d_gradT_dL(k) = one_m_f*s% d_gradr_dL(k)
               s% d_gradT_dw_div_wc(k) = one_m_f*s% d_gradr_dw_div_wc(k)
               s% d_gradT_dlnd00(k) = &
                  f*alfa*s% d_eos_dlnd(i_grad_ad, k) + one_m_f*s% d_gradr_dlnd00(k)
               s% d_gradT_dlnT00(k) = &
                  f*alfa*s% d_eos_dlnT(i_grad_ad, k) + one_m_f*s% d_gradr_dlnT00(k)
               if (k > 1) then
                  s% d_gradT_dlndm1(k) = &
                     f*beta*s% d_eos_dlnd(i_grad_ad, k-1) + one_m_f*s% d_gradr_dlndm1(k)
                  s% d_gradT_dlnTm1(k) = &
                     f*beta*s% d_eos_dlnT(i_grad_ad, k-1) + one_m_f*s% d_gradr_dlnTm1(k)
               end if
            end if
            if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = 0.0d0

         end if

! fxt 02jun2019
! gradT_sub_grada was saved before calling adjust_gradT_excess,
! but adjust_gradT_excess which adjusts gradT(k), making gradT(k) inconsistent with gradT_sub_grada
!         s% gradT_sub_grada(k) = s% gradT(k) - s% grada_face(k) ! gradT_excess
!         call adjust_gradT_excess(s, k)

! just recalculate gradT_sub_grada after calling adjust_gradT_excess
         s% gradT_sub_grada(k) = s% gradT(k) - s% grada_face(k) ! gradT_excess
         call adjust_gradT_excess(s, k)
         s% gradT_sub_grada(k) = s% gradT(k) - s% grada_face(k) ! new gradT_excess from adjusted gradT
         
         if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
            write(*,2) 'after adjust_gradT_excess gradT', k, s% gradT(k)
         end if

      end subroutine adjust_gradT_fraction


      subroutine adjust_gradT_excess(s, k)
         use eos_def
         type (star_info), pointer :: s
         integer, intent(in) :: k

         real(dp) :: alfa, beta, log_tau, gradT_excess_alpha, &
            d_grada_face_dlnd00, d_grada_face_dlnT00, &
            d_grada_face_dlndm1, d_grada_face_dlnTm1

         include 'formats'

         !s% gradT_excess_alpha is calculated at start of step and held constant during iterations
         ! gradT_excess_alpha = 0 means no efficiency boost; = 1 means full efficiency boost

         gradT_excess_alpha = s% gradT_excess_alpha
         s% gradT_excess_effect(k) = 0.0d0

         if (gradT_excess_alpha <= 0.0  .or. &
             s% gradT_sub_grada(k) <= s% gradT_excess_f1) then
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,*) 'gradT_excess_alpha <= 0.0', gradT_excess_alpha <= 0.0
               write(*,*) 's% gradT_sub_grada(k) <= s% gradT_excess_f1', &
                  s% gradT_sub_grada(k) <= s% gradT_excess_f1, &
                  s% gradT_sub_grada(k), s% gradT(k), s% grada_face(k), s% gradT_excess_f1
            end if
            return
         end if

         if (s% lnT(k)/ln10 > s% gradT_excess_max_logT) then
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,*) 's% lnT(k)/ln10 > s% gradT_excess_max_logT', k, s% lnT(k)/ln10, s% gradT_excess_max_logT
            end if
            return
         end if

         log_tau = log10(s% tau(k))
         if (log_tau < s% gradT_excess_max_log_tau_full_off) then
            if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
               write(*,*) 'log_tau < s% gradT_excess_max_log_tau_full_off', k, log_tau, s% gradT_excess_max_log_tau_full_off
            end if
            return
         end if

         ! boost efficiency of energy transport

         ! alfa is the fraction coming from k; (1-alfa) from k-1.
         if (k == 1) then
            alfa = 1.0d0
         else
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         end if
         beta = 1.0d0 - alfa

         ! grada_face = alfa*s% grada(k) + beta*s% grada(k-1)
         d_grada_face_dlnd00 = alfa*s% d_eos_dlnd(i_grad_ad, k)
         d_grada_face_dlnT00 = alfa*s% d_eos_dlnT(i_grad_ad, k)
         if (k > 1) then
            d_grada_face_dlndm1 = beta*s% d_eos_dlnd(i_grad_ad, k-1)
            d_grada_face_dlnTm1 = beta*s% d_eos_dlnT(i_grad_ad, k-1)
         else
            d_grada_face_dlndm1 = 0.0d0
            d_grada_face_dlnTm1 = 0.0d0
         end if

         if (log_tau < s% gradT_excess_min_log_tau_full_on) &
            gradT_excess_alpha = gradT_excess_alpha* &
               (log_tau - s% gradT_excess_max_log_tau_full_off)/ &
               (s% gradT_excess_min_log_tau_full_on - s% gradT_excess_max_log_tau_full_off)

         alfa = s% gradT_excess_f2 ! for full boost, use this fraction of gradT
         if (gradT_excess_alpha < 1) then ! only partial boost, so increase alfa
            ! alfa goes to 1 as gradT_excess_alpha goes to 0
            ! alfa unchanged as gradT_excess_alpha goes to 1
            alfa = alfa + (1d0 - alfa)*(1d0 - gradT_excess_alpha)
         end if
         beta = 1.d0 - alfa

         if (k==s% x_integer_ctrl(19) .and. s% x_integer_ctrl(19) > 0 .and. &
                     s% solver_iter == s% x_integer_ctrl(20) .and. (s% model_number == s% x_integer_ctrl(21) .or. s% x_integer_ctrl(21) == 0)) then
            write(*,2) 'adjust_gradT_excess alfa gradT grada', k, alfa, s% gradT(k), s% grada_face(k)
         end if

         s% gradT(k) = alfa*s% gradT(k) + beta*s% grada_face(k)
         s% gradT_excess_effect(k) = beta

         s% d_gradT_dlnR(k) = alfa*s% d_gradT_dlnR(k)
         s% d_gradT_dL(k) = alfa*s% d_gradT_dL(k)
         s% d_gradT_dw_div_wc(k) = alfa*s% d_gradT_dw_div_wc(k)
         if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = alfa*s% d_gradT_dln_cvpv0(k)

         s% d_gradT_dlnd00(k) = &
            alfa*s% d_gradT_dlnd00(k) + beta*d_grada_face_dlnd00
         s% d_gradT_dlnT00(k) = &
            alfa*s% d_gradT_dlnT00(k) + beta*d_grada_face_dlnT00

         if (k > 1) then
            s% d_gradT_dlndm1(k) = &
               alfa*s% d_gradT_dlndm1(k) + beta*d_grada_face_dlndm1
            s% d_gradT_dlnTm1(k) = &
               alfa*s% d_gradT_dlnTm1(k) + beta*d_grada_face_dlnTm1
         end if

      end subroutine adjust_gradT_excess


      subroutine set_grads_newer(s, ierr)
         use chem_def, only: chem_isos
         use star_utils, only: smooth
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, nz, j, cid, max_cid
         real(dp) :: val, max_val, A, Z
         real(dp), pointer, dimension(:) :: dlnP, dlnd, dlnT

         include 'formats'

         ierr = 0
         nz = s% nz
         call do_alloc(ierr)
         if (ierr /= 0) return

         do k = 2, nz
            dlnP(k) = s% lnPeos(k-1) - s% lnPeos(k)
            dlnd(k) = s% lnd(k-1) - s% lnd(k)
            dlnT(k) = s% lnT(k-1) - s% lnT(k)
         end do
         dlnP(1) = dlnP(2)
         dlnd(1) = dlnd(2)
         dlnT(1) = dlnT(2)

         call smooth(dlnP,nz)
         call smooth(dlnd,nz)
         call smooth(dlnT,nz)

         s% grad_density(1) = 0
         s% grad_temperature(1) = 0
         do k = 2, nz
            if (dlnP(k) >= 0) then
               s% grad_density(k) = 0
               s% grad_temperature(k) = 0
            else
               s% grad_density(k) = dlnd(k)/dlnP(k)
               s% grad_temperature(k) = dlnT(k)/dlnP(k)
            end if
         end do

         call smooth(s% grad_density,nz)
         call smooth(s% grad_temperature,nz)

         if (s% use_Ledoux_criterion .and. s% calculate_Brunt_B) then
            do k=1,nz
               s% gradL_composition_term(k) = s% unsmoothed_brunt_B(k)
            end do
            call smooth_gradL_composition_term
         else
            do k=1,nz
               s% gradL_composition_term(k) = 0d0
            end do
         end if

         call dealloc

         do k=3,nz-2
            max_cid = 0
            max_val = -1d99
            do j=1,s% species
               cid = s% chem_id(j)
               A = dble(chem_isos% Z_plus_N(cid))
               Z = dble(chem_isos% Z(cid))
               val = (s% xa(j,k-2) + s% xa(j,k-1) - s% xa(j,k) - s% xa(j,k+1))*(1d0 + Z)/A
               if (val > max_val) then
                  max_val = val
                  max_cid = cid
               end if
            end do
            s% dominant_iso_for_thermohaline(k) = max_cid
         end do
         s% dominant_iso_for_thermohaline(1:2) = &
            s% dominant_iso_for_thermohaline(3)
         s% dominant_iso_for_thermohaline(nz-1:nz) = &
            s% dominant_iso_for_thermohaline(nz-2)


         contains

         subroutine smooth_gradL_composition_term
            use star_utils, only: weighed_smoothing, threshold_smoothing
            logical, parameter :: preserve_sign = .false.
            real(dp), pointer, dimension(:) :: work
            integer :: k
            include 'formats'
            ierr = 0
            work => dlnd
            if (s% num_cells_for_smooth_gradL_composition_term <= 0) return
            call threshold_smoothing( &
               s% gradL_composition_term, s% threshold_for_smooth_gradL_composition_term, s% nz, &
               s% num_cells_for_smooth_gradL_composition_term, preserve_sign, work)
         end subroutine smooth_gradL_composition_term

         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            call do_work_arrays(.true.,ierr)
         end subroutine do_alloc

         subroutine dealloc
            call do_work_arrays(.false.,ierr)
         end subroutine dealloc

         subroutine do_work_arrays(alloc_flag, ierr)
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               dlnP, nz, nz_alloc_extra, 'mlt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dlnd, nz, nz_alloc_extra, 'mlt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dlnT, nz, nz_alloc_extra, 'mlt', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays

      end subroutine set_grads_newer


      subroutine switch_to_no_mixing(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% mlt_mixing_type(k) = no_mixing
         s% mlt_mixing_length(k) = 0
         s% mlt_D(k) = 0
         s% mlt_cdc(k) = 0d0
         s% mlt_vc(k) = 0
      end subroutine switch_to_no_mixing


      subroutine switch_to_radiative_newer(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         call switch_to_no_mixing(s,k)
         s% gradT(k) = s% gradr(k)
         s% d_gradT_dlnd00(k) = s% d_gradr_dlnd00(k)
         s% d_gradT_dlnT00(k) = s% d_gradr_dlnT00(k)
         s% d_gradT_dlndm1(k) = s% d_gradr_dlndm1(k)
         s% d_gradT_dlnTm1(k) = s% d_gradr_dlnTm1(k)
         s% d_gradT_dlnR(k) = s% d_gradr_dlnR(k)
         s% d_gradT_dL(k) = s% d_gradr_dL(k)
         if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = 0d0
         s% d_gradT_dw_div_wc(k) = s% d_gradr_dw_div_wc(k)
      end subroutine switch_to_radiative_newer


      subroutine switch_to_adiabatic(s,k)
         use eos_def, only: i_grad_ad
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% gradT(k) = s% grada(k)
         s% d_gradT_dlnd00(k) = s% d_eos_dlnd(i_grad_ad, k)
         s% d_gradT_dlnT00(k) = s% d_eos_dlnT(i_grad_ad, k)
         s% d_gradT_dlndm1(k) = 0
         s% d_gradT_dlnTm1(k) = 0
         s% d_gradT_dlnR(k) = 0
         s% d_gradT_dL(k) = 0
         if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = 0d0
         s% d_gradT_dw_div_wc(k) = 0
      end subroutine switch_to_adiabatic


      subroutine set_gradT_excess_alpha_newer(s, ierr)
         use alloc
         use star_utils, only: get_Lrad_div_Ledd, after_C_burn
         use chem_def, only: ih1, ihe4

         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         real(dp) :: beta, lambda, phi, tmp, alpha, alpha2
         real(dp) :: &
            beta_limit, &
            lambda1, &
            beta1, &
            lambda2, &
            beta2, &
            dlambda, &
            dbeta

         integer :: k, k_beta, k_lambda, nz, h1, he4

         logical, parameter :: dbg = .false.

         include 'formats'

         if (dbg) write(*,*) 'enter set_gradT_excess_alpha_newer'

         ierr = 0
         if (.not. s% okay_to_reduce_gradT_excess) then
            s% gradT_excess_alpha = 0
            if (dbg) write(*,1) 'okay_to_reduce_gradT_excess'
            return
         end if

         nz = s% nz

         h1 = s% net_iso(ih1)
         if (h1 /= 0) then
            if (s% xa(h1,nz) > s% gradT_excess_max_center_h1) then
               s% gradT_excess_alpha = 0
               if (dbg) write(*,1) 'gradT_excess_max_center_h1'
               return
            end if
         end if

         he4 = s% net_iso(ihe4)
         if (he4 /= 0) then
            if (s% xa(he4,nz) < s% gradT_excess_min_center_he4) then
               s% gradT_excess_alpha = 0
               if (dbg) write(*,1) 'gradT_excess_min_center_he4'
               return
            end if
         end if

         s% gradT_excess_alpha = 1d0

      end subroutine set_gradT_excess_alpha_newer


      subroutine check_for_redo_MLT_newer(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr

         logical :: in_convective_region
         integer :: k, k_bot
         real(dp) :: bot_Hp, bot_r, top_Hp, top_r, dr
         logical :: dbg

         include 'formats'

         ! check_for_redo_MLT_newer assumes that nzlo = 1, nzhi = nz
         ! that is presently true; make sure that assumption doesn't change
         if (.not. ((nzlo.eq.1).and.(nzhi.eq.s%nz))) then
            write(*,*) 'nzlo != 1 or nzhi != nz'
            call mesa_error(__FILE__,__LINE__)
         endif

         ierr = 0
         dbg = .false.

         bot_Hp = 0; bot_r = 0; top_Hp = 0; top_r = 0; dr = 0

         in_convective_region = (s% mlt_mixing_type(nzhi) == convective_mixing)
         k_bot = nzhi
         bot_r = s% r(k_bot)
         bot_Hp = s% scale_height(k_bot)

         do k=nzhi-1, nzlo+1, -1
            if (in_convective_region) then
               if (s% mlt_mixing_type(k) /= convective_mixing) then
                  call end_of_convective_region
               end if
            else ! in non-convective region
               if (s% mlt_mixing_type(k) == convective_mixing) then ! start of a convective region
                  k_bot = k+1
                  in_convective_region = .true.
                  bot_r = s% r(k_bot)
                  bot_Hp = s% scale_height(k_bot)
               end if
            end if
         end do

         if (in_convective_region) then
            k = 1 ! end at top
            call end_of_convective_region
         end if


         contains


         subroutine end_of_convective_region()
            integer :: kk, op_err, mix_type
            real(dp) :: Hp
            logical :: end_dbg

            9 format(a40, 3i7, 99(1pd26.16))
            include 'formats'

            in_convective_region = .false.

            end_dbg = .false.

            top_r = s% r(k)
            top_Hp = s% scale_height(k)
            dr = top_r - bot_r
            Hp = (bot_Hp + top_Hp)/2

            if (dr < s% alpha_mlt(k)*min(top_Hp, bot_Hp) .and. &
                  s% redo_conv_for_dr_lt_mixing_length) then
!$OMP PARALLEL DO PRIVATE(kk,op_err) SCHEDULE(dynamic,2)
               do kk = k, k_bot
                  op_err = 0
                  call redo1_mlt(s,kk,dr,op_err)
                  if (op_err /= 0) ierr = op_err
               end do
!$OMP END PARALLEL DO
            end if

         end subroutine end_of_convective_region


         subroutine redo1_mlt(s,k,dr,ierr)
            type (star_info), pointer :: s
            integer, intent(in) :: k
            real(dp), intent(in) :: dr
            integer, intent(out) :: ierr
            real(dp) :: Hp, reduced_alpha, &
               gradL_composition_term, opacity_face, chiRho_face, &
               chiT_face, Cp_face, grada_face, P_face, xh_face
            logical :: from_do1_mlt, make_gradr_sticky_in_solver_iters
            include 'formats'
            ierr = 0
            if (dr >= s% mlt_mixing_length(k)) return
            ! if convection zone is smaller than mixing length
            ! redo MLT with reduced alpha so mixing_length = dr
            Hp = s% scale_height(k)
            reduced_alpha = dr/Hp
            gradL_composition_term = -1
            opacity_face = -1
            chiRho_face = -1
            chiT_face = -1
            Cp_face = -1
            grada_face = -1
            P_face = -1
            xh_face = -1
            from_do1_mlt = .false.
            call do1_mlt_2_newer(s, k, reduced_alpha, &
               gradL_composition_term, opacity_face, chiRho_face, &
               chiT_face, Cp_face, grada_face, P_face, xh_face, &
               from_do1_mlt, make_gradr_sticky_in_solver_iters, ierr)
         end subroutine redo1_mlt


      end subroutine check_for_redo_MLT_newer


      end module mlt_info_newer
