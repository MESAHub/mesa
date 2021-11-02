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


      module turb_info

      use star_private_def
      use const_def
      use num_lib
      use utils_lib
      use auto_diff_support

      implicit none

      private
      public :: &
         set_mlt_vars, & ! for hydro_vars and conv_premix
         do1_mlt_2, & ! for predictive_mix
         switch_to_radiative, & ! mix_info
         check_for_redo_MLT, & ! for hydro_vars
         set_gradT_excess_alpha ! for evolve


      contains


      subroutine set_mlt_vars(s, nzlo, nzhi, ierr)
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: k, op_err
         integer(8) :: time0
         real(dp) :: total
         logical :: make_gradr_sticky_in_solver_iters
         include 'formats'
         ierr = 0
         if (s% doing_timing) call start_time(s, time0, total)
!$OMP PARALLEL DO PRIVATE(k,op_err,make_gradr_sticky_in_solver_iters) SCHEDULE(dynamic,2)
         do k = nzlo, nzhi
            op_err = 0
            call do1_mlt_2(s, k, make_gradr_sticky_in_solver_iters, op_err)
            if (make_gradr_sticky_in_solver_iters .and. s% solver_iter > 3) then
               if (.not. s% fixed_gradr_for_rest_of_solver_iters(k)) then
                  s% fixed_gradr_for_rest_of_solver_iters(k) = &
                     (s% mlt_mixing_type(k) == no_mixing)
               end if
            end if            
         end do
!$OMP END PARALLEL DO
         if (s% doing_timing) call update_time(s, time0, total, s% time_mlt)

      end subroutine set_mlt_vars


      subroutine do1_mlt_2(s, k, &
            make_gradr_sticky_in_solver_iters, ierr, &
            mixing_length_alpha_in, gradL_composition_term_in)
         ! get convection info for point k
         use star_utils
         use turb_support, only: do1_mlt_eval
         use eos_def
         use chem_def, only: ih1
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(out) :: make_gradr_sticky_in_solver_iters
         integer, intent(out) :: ierr
         real(dp), intent(in), optional :: &
            mixing_length_alpha_in, gradL_composition_term_in

         type(auto_diff_real_star_order1) :: gradr_factor
         real(dp) :: v, f, xh_face, &
            gradL_composition_term, abs_du_div_cs, cs, mixing_length_alpha
         real(dp), pointer :: vel(:)
         integer :: i, mixing_type, h1, nz, k_T_max
         real(dp), parameter :: conv_vel_mach_limit = 0.9d0
         logical :: no_mix
         type(auto_diff_real_star_order1) :: &
            grada_face_ad, scale_height_ad, gradr_ad, rho_face_ad, &
            gradT_ad, Y_face_ad, mlt_vc_ad, D_ad, Gamma_ad
         include 'formats'

         ierr = 0
         nz = s% nz
         
         if (k < 1 .or. k > nz) then
            write(*,3) 'bad k for do1_mlt', k, nz
            ierr = -1
            return
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (present(mixing_length_alpha_in)) then
            mixing_length_alpha = mixing_length_alpha_in
         else
            mixing_length_alpha = s% mixing_length_alpha
         end if

         if (present(gradL_composition_term_in)) then
            gradL_composition_term = gradL_composition_term_in      
         else if (s% use_Ledoux_criterion) then
            gradL_composition_term = s% gradL_composition_term(k)
         else
            gradL_composition_term = 0d0
         end if

         grada_face_ad = get_grada_face(s,k)
         scale_height_ad = get_scale_height_face(s,k)
         gradr_ad = get_gradr_face(s,k)

         if (s% rotation_flag .and. s% mlt_use_rotation_correction) then
            gradr_factor = s% ft_rot(k)/s% fp_rot(k)*s% gradr_factor(k)
         else
            gradr_factor = s% gradr_factor(k)
         end if
         if (is_bad_num(gradr_factor% val)) then
            ierr = -1
            return
         end if
         gradr_ad = gradr_ad*gradr_factor
         
         ! now can call set_no_mixing if necessary
         
         if (k == 1 .and. s% mlt_make_surface_no_mixing) then
            call set_no_mixing('surface_no_mixing')
            return
         end if

         if (s% lnT_start(k)/ln10 > s% max_logT_for_mlt) then
            call set_no_mixing('max_logT')
            return
         end if

         if (s% phase(k) > 0.5d0) then
            call set_no_mixing('solid_no_mixing')
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

         if (s% csound_start(k) > 0d0 .and. (s% u_flag .or. s% v_flag)) then
            no_mix = .false.
            if (s% u_flag) then
               vel => s% u_start
            else
               vel => s% v_start
            end if
            abs_du_div_cs = 0d0
            if (vel(k)/1d5 > s% max_v_for_convection) then
               no_mix = .true.
            else if (s% q(k) > s% max_q_for_convection_with_hydro_on) then
               no_mix = .true.
            else if ((abs(vel(k))) >= &
                  s% csound_start(k)*s% max_v_div_cs_for_convection) then
               no_mix = .true.              
            else if (s% u_flag) then
               if (k == 1) then
                  abs_du_div_cs = 1d99
               else if (k < nz) then
                  abs_du_div_cs = max(abs(vel(k) - vel(k+1)), &
                      abs(vel(k) - vel(k-1))) / s% csound_start(k)
               end if
               if (abs_du_div_cs > s% max_abs_du_div_cs_for_convection) then
                  no_mix = .true.
               end if
            end if
            if (no_mix) then
               call set_no_mixing('no_mix')
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
         if (make_gradr_sticky_in_solver_iters .and. s% fixed_gradr_for_rest_of_solver_iters(k)) then
            call set_no_mixing('gradr_sticky')
            return
         end if
            
         call do1_mlt_eval(s, k, s% MLT_option, gradL_composition_term, &
            gradr_ad, grada_face_ad, scale_height_ad, mixing_length_alpha, &
            mixing_type, gradT_ad, Y_face_ad, mlt_vc_ad, D_ad, Gamma_ad, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*,*) 'ierr in do1_mlt_eval for k', k
            end if
            return
         end if
         
         call store_results

         if (s% mlt_gradT_fraction >= 0d0 .and. s% mlt_gradT_fraction <= 1d0) then
            f = s% mlt_gradT_fraction
         else
            f = s% adjust_mlt_gradT_fraction(k)
         end if
         call adjust_gradT_fraction(s, k, f)
         
         if (s% mlt_mixing_type(k) == no_mixing .or. abs(s% gradr(k)) < 1d-20) then
            s% L_conv(k) = 0d0
         else
            s% L_conv(k) = s% L(k) * (1d0 - s% gradT(k)/s% gradr(k)) ! C&G 14.109            
         end if

         contains         
         
         subroutine store_results
            s% mlt_mixing_type(k) = mixing_type
         
            s% grada_face_ad(k) = grada_face_ad
            s% grada_face(k) = grada_face_ad%val
         
            s% gradT_ad(k) = gradT_ad
            s% gradT(k) = s% gradT_ad(k)%val
            s% mlt_gradT(k) = s% gradT(k) ! prior to adjustments
         
            s% Y_face_ad(k) = Y_face_ad
            s% Y_face(k) = s% Y_face_ad(k)%val
         
            s% mlt_vc_ad(k) = mlt_vc_ad
            if (s% okay_to_set_mlt_vc) s% mlt_vc(k) = s% mlt_vc_ad(k)%val  
         
            s% mlt_D_ad(k) = D_ad
            s% mlt_D(k) = D_ad%val

            rho_face_ad = get_rho_face(s,k)
            s% mlt_cdc(k) = s% mlt_D(k)*pow2(pi4*pow2(s%r(k))*rho_face_ad%val)
         
            s% mlt_Gamma_ad(k) = Gamma_ad
            s% mlt_Gamma(k) = Gamma_ad%val
         
            s% gradr_ad(k) = gradr_ad
            s% gradr(k) = s% gradr_ad(k)%val
         
            s% gradL_ad(k) = s% grada_face_ad(k) + gradL_composition_term
            s% gradL(k) = s% gradL_ad(k)%val
         
            s% scale_height_ad(k) = scale_height_ad
            s% scale_height(k) = scale_height_ad%val     
                 
            s% Lambda_ad(k) = mixing_length_alpha*scale_height_ad
            s% mlt_mixing_length(k) = s% Lambda_ad(k)%val
            
         end subroutine store_results

         subroutine set_no_mixing(str)
            character (len=*) :: str
            include 'formats'
            
            s% mlt_mixing_type(k) = no_mixing
            
            s% grada_face_ad(k) = grada_face_ad
            s% grada_face(k) = grada_face_ad%val
            
            gradT_ad = gradr_ad
            s% gradT_ad(k) = gradT_ad
            s% gradT(k) = s% gradT_ad(k)%val
            
            Y_face_ad = gradT_ad - grada_face_ad
            s% Y_face_ad(k) = Y_face_ad
            s% Y_face(k) = s% Y_face_ad(k)%val
            
            s% mlt_vc_ad(k) = 0d0
            if (s% okay_to_set_mlt_vc) s% mlt_vc(k) = 0d0
            
            s% mlt_D_ad(k) = 0d0
            s% mlt_D(k) = 0d0
            s% mlt_cdc(k) = 0d0
            
            s% mlt_Gamma_ad(k) = 0d0
            s% mlt_Gamma(k) = 0d0
            
            s% gradr_ad(k) = gradr_ad
            s% gradr(k) = s% gradr_ad(k)%val
               
            s% gradL_ad(k) = 0d0
            s% gradL(k) = 0d0
            
            s% scale_height_ad(k) = scale_height_ad
            s% scale_height(k) = scale_height_ad%val    
                     
            s% Lambda_ad(k) = mixing_length_alpha*scale_height_ad
            s% mlt_mixing_length(k) = s% Lambda_ad(k)%val

            s% L_conv(k) = 0d0
            
         end subroutine set_no_mixing

      end subroutine do1_mlt_2


      subroutine adjust_gradT_fraction(s,k,f)
         ! replace gradT by combo of grada_face and gradr
         ! then check excess
         use eos_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: f
         integer, intent(in) :: k
         include 'formats'         
         if (f >= 0.0 .and. f <= 1.0) then
            if (f == 0d0) then
               s% gradT_ad(k) = s% gradr_ad(k)
            else ! mix
               s% gradT_ad(k) = f*s% grada_face_ad(k) + (1.0d0 - f)*s% gradr_ad(k)
            end if
            s% gradT(k) = s% gradT_ad(k)%val
         end if
         call adjust_gradT_excess(s, k)         
         s% gradT_sub_grada(k) = s% gradT(k) - s% grada_face(k)
      end subroutine adjust_gradT_fraction


      subroutine adjust_gradT_excess(s, k)
         use eos_def
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: alfa, log_tau, gradT_excess_alpha, gradT_sub_grada
         include 'formats'
         !s% gradT_excess_alpha is calculated at start of step and held constant during iterations
         ! gradT_excess_alpha = 0 means no efficiency boost; = 1 means full efficiency boost
         gradT_excess_alpha = s% gradT_excess_alpha
         s% gradT_excess_effect(k) = 0.0d0         
         gradT_sub_grada = s% gradT(k) - s% grada_face(k)
         if (gradT_excess_alpha <= 0.0  .or. &
             gradT_sub_grada <= s% gradT_excess_f1) return
         if (s% lnT(k)/ln10 > s% gradT_excess_max_logT) return
         log_tau = log10(s% tau(k))
         if (log_tau < s% gradT_excess_max_log_tau_full_off) return
         if (log_tau < s% gradT_excess_min_log_tau_full_on) &
            gradT_excess_alpha = gradT_excess_alpha* &
               (log_tau - s% gradT_excess_max_log_tau_full_off)/ &
               (s% gradT_excess_min_log_tau_full_on - s% gradT_excess_max_log_tau_full_off)
         alfa = s% gradT_excess_f2 ! for full boost, use this fraction of gradT
         if (gradT_excess_alpha < 1) & ! only partial boost, so increase alfa
            ! alfa goes to 1 as gradT_excess_alpha goes to 0
            ! alfa unchanged as gradT_excess_alpha goes to 1
            alfa = alfa + (1d0 - alfa)*(1d0 - gradT_excess_alpha)
         s% gradT_ad(k) = alfa*s% gradT_ad(k) + (1d0 - alfa)*s% grada_face_ad(k)
         s% gradT(k) = s% gradT_ad(k)%val
         s% gradT_excess_effect(k) = 1d0 - alfa
      end subroutine adjust_gradT_excess


      subroutine switch_to_radiative(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% mlt_mixing_type(k) = no_mixing
         s% mlt_mixing_length(k) = 0
         s% mlt_D(k) = 0
         s% mlt_cdc(k) = 0d0
         s% mlt_vc(k) = 0
         s% gradT_ad(k) = s% gradr_ad(k)
         s% gradT(k) = s% gradT_ad(k)%val
      end subroutine switch_to_radiative


      subroutine switch_to_adiabatic(s,k)
         use eos_def, only: i_grad_ad
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% gradT_ad(k) = s% grada_face_ad(k)
         s% gradT(k) = s% gradT_ad(k)%val
      end subroutine switch_to_adiabatic


      subroutine set_gradT_excess_alpha(s, ierr)
         use alloc
         use star_utils, only: get_Lrad_div_Ledd, after_C_burn
         use chem_def, only: ih1, ihe4
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: beta, lambda, phi, tmp, alpha, alpha2, &
            beta_limit, lambda1, beta1, lambda2, beta2, dlambda, dbeta
         integer :: k, k_beta, k_lambda, nz, h1, he4
         include 'formats'
         ierr = 0
         if (.not. s% okay_to_reduce_gradT_excess) then
            s% gradT_excess_alpha = 0
            return
         end if
         nz = s% nz
         h1 = s% net_iso(ih1)
         if (h1 /= 0) then
            if (s% xa(h1,nz) > s% gradT_excess_max_center_h1) then
               s% gradT_excess_alpha = 0
               return
            end if
         end if
         he4 = s% net_iso(ihe4)
         if (he4 /= 0) then
            if (s% xa(he4,nz) < s% gradT_excess_min_center_he4) then
               s% gradT_excess_alpha = 0
               return
            end if
         end if
         beta = 1d0 ! beta = min over k of Pgas(k)/Peos(k)
         k_beta = 0
         do k=1,nz
            tmp = s% Pgas(k)/s% Peos(k)
            if (tmp < beta) then
               k_beta = k
               beta = tmp
            end if
         end do
         beta = beta*(1d0 + s% xa(1,nz))
         s% gradT_excess_min_beta = beta
         lambda = 0d0 ! lambda = max over k of Lrad(k)/Ledd(k)
         do k=2,k_beta
            tmp = get_Lrad_div_Ledd(s,k)
            if (tmp > lambda) then
               k_lambda = k
               lambda = tmp
            end if
         end do
         lambda = min(1d0,lambda)
         s% gradT_excess_max_lambda = lambda
         lambda1 = s% gradT_excess_lambda1
         beta1 = s% gradT_excess_beta1
         lambda2 = s% gradT_excess_lambda2
         beta2 = s% gradT_excess_beta2
         dlambda = s% gradT_excess_dlambda
         dbeta = s% gradT_excess_dbeta
         ! alpha is fraction of full boost to apply
         ! depends on location in (beta,lambda) plane
         if (lambda1 < 0) then
            alpha = 1
         else if (lambda >= lambda1) then
            if (beta <= beta1) then
               alpha = 1
            else if (beta < beta1 + dbeta) then
               alpha = (beta1 + dbeta - beta)/dbeta
            else ! beta >= beta1 + dbeta
               alpha = 0
            end if
         else if (lambda >= lambda2) then
            beta_limit = beta2 + &
               (lambda - lambda2)*(beta1 - beta2)/(lambda1 - lambda2)
            if (beta <= beta_limit) then
               alpha = 1
            else if (beta < beta_limit + dbeta) then
               alpha = (beta_limit + dbeta - beta)/dbeta
            else
               alpha = 0
            end if
         else if (lambda > lambda2 - dlambda) then
            if (beta <= beta2) then
               alpha = 1
            else if (beta < beta2 + dbeta) then
               alpha = (lambda - (lambda2 - dlambda))/dlambda
            else ! beta >= beta2 + dbeta
               alpha = 0
            end if
         else ! lambda <= lambda2 - dlambda
            alpha = 0
         end if
         if (s% generations > 1 .and. lambda1 >= 0) then ! time smoothing
            s% gradT_excess_alpha = &
               (1d0 - s% gradT_excess_age_fraction)*alpha + &
               s% gradT_excess_age_fraction*s% gradT_excess_alpha_old
            if (s% gradT_excess_max_change > 0d0) then
               if (s% gradT_excess_alpha > s% gradT_excess_alpha_old) then
                  s% gradT_excess_alpha = min(s% gradT_excess_alpha, s% gradT_excess_alpha_old + &
                     s% gradT_excess_max_change)
               else
                  s% gradT_excess_alpha = max(s% gradT_excess_alpha, s% gradT_excess_alpha_old - &
                     s% gradT_excess_max_change)
               end if
            end if
         else
            s% gradT_excess_alpha = alpha
         end if
         if (s% gradT_excess_alpha < 1d-4) s% gradT_excess_alpha = 0d0
         if (s% gradT_excess_alpha > 0.9999d0) s% gradT_excess_alpha = 1d0
      end subroutine set_gradT_excess_alpha


      subroutine check_for_redo_MLT(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         logical :: in_convective_region
         integer :: k, k_bot
         real(dp) :: bot_Hp, bot_r, top_Hp, top_r, dr
         logical :: dbg
         include 'formats'
         ! check_for_redo_MLT assumes that nzlo = 1, nzhi = nz
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
               if (s% mlt_mixing_type(k) == convective_mixing) then 
                  ! start of a convective region
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

         subroutine redo1_mlt(s, k, dr, ierr)
            type (star_info), pointer :: s
            integer, intent(in) :: k
            real(dp), intent(in) :: dr
            integer, intent(out) :: ierr
            logical :: make_gradr_sticky_in_solver_iters
            include 'formats'
            ierr = 0
            if (dr >= s% mlt_mixing_length(k)) return
            ! if convection zone is smaller than mixing length
            ! redo MLT with reduced alpha so mixing_length = dr
            call do1_mlt_2(s, k, make_gradr_sticky_in_solver_iters, ierr, &
               mixing_length_alpha_in = dr/s% scale_height(k))
         end subroutine redo1_mlt

      end subroutine check_for_redo_MLT


      end module turb_info
