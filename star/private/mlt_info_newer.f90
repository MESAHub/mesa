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






!   logical :: have_mlt_vc
!      set false in init
!      set true by read_model if mlt_vc is in .mod file
!      if true, then mlt_vc included in .mod file by write_model
!      set true in finish_load model (mlt_vc = 0 if not already set)
!      must be true before call mlt_info
!   logical :: okay_to_set_mlt_vc
!      set false in init
!      set false in evolve step at start
!      checked in mlt_info. only set mlt_vc if this is true
!      set true in prepare_for_new_try in evolve after set mlt_vc_old
!   real(dp) :: mlt_vc(:)
!      input value for step containing initial values needed by analytic eqn
!   real(dp) :: mlt_vc_old(:)
!      after mlt_vc has been optionally remeshed, copy it to mlt_vc_old
!      the analytic eqn uses mlt_vc_old as starting value
!      also used for retry/redo to reset mlt_vc
!   type(auto_diff_real_star_order1 :: mlt_vc_ad(:)
!      result from mlt along with gradT_ad and gradr_ad
!      mlt_vc(k) = mlt_vc_ad(k)%val

   


      module mlt_info_newer

      use star_private_def
      use const_def
      use num_lib
      use utils_lib
      use auto_diff_support
      use star_utils

      implicit none

      private
      public :: &
         set_mlt_vars_newer, &            ! conv_premix.f90, hydro_vars.f90
         set_grads_newer, &               ! conv_premix.f90
         do1_mlt_newer, &                 ! predictive_mix.f90
         switch_to_radiative_newer, &     ! mix_info.f90
         check_for_redo_MLT_newer, &      ! hydro_vars.f90
         set_gradT_excess_alpha_newer  ! evolve.f90

      contains


      subroutine set_mlt_vars_newer(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: k, op_err
         integer(8) :: time0
         real(dp) :: total
         logical :: make_gradr_sticky_in_solver_iters
         include 'formats'
         ierr = 0
         !write(*,*) 'set_mlt_vars_newer'
         if (s% doing_timing) call start_time(s, time0, total)
         !call set_gradT_excess_alpha_newer(s, ierr)  evolve set_start_of_step_info -> set_gradT_excess_alpha -> this 
         !if (ierr /= 0) then
         !   if (s% report_ierr) write(*,2) 'set_gradT_excess_alpha failed', k
         !   return
         !end if
         if (s% using_TDC .and. .not. s% have_mlt_vc) then
            write(*,*) 's% using_TDC .and. .not. s% have_mlt_vc'
            stop 'set_mlt_vars'
         end if
!$OMP PARALLEL DO PRIVATE(k,op_err,make_gradr_sticky_in_solver_iters) SCHEDULE(dynamic,2)
         do k = nzlo, nzhi
            op_err = 0
            call do1_mlt_newer(s, k, s% alpha_mlt(k), make_gradr_sticky_in_solver_iters, op_err)
            if (op_err /= 0) then
               ierr = op_err
               if (s% report_ierr) write(*,2) 'set_mlt_vars failed', k
            end if
            if (make_gradr_sticky_in_solver_iters .and. s% solver_iter > 3) then
               if (.not. s% newer_fixed_gradr_for_rest_of_solver_iters(k)) &
                  s% newer_fixed_gradr_for_rest_of_solver_iters(k) = &
                     (s% mlt_mixing_type(k) == no_mixing)
            end if            
         end do
!$OMP END PARALLEL DO
         if (s% doing_timing) call update_time(s, time0, total, s% time_mlt)
      end subroutine set_mlt_vars_newer


      subroutine do1_mlt_newer( &
            s, k, mixing_length_alpha, make_gradr_sticky_in_solver_iters, ierr)
         use mlt_get_results_newer, only: do1_mlt_eval_newer
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha
         logical, intent(out) :: make_gradr_sticky_in_solver_iters
         integer, intent(out) :: ierr
         real(dp), pointer :: vel(:)
         real(dp) :: cs, gradL_composition_term, abs_du_div_cs, gradr_factor
         integer :: i, mixing_type, k_T_max
         logical :: just_gradr
         type(auto_diff_real_star_order1) :: grada_face_ad, rho_face, &
            gradT_ad, Y_face_ad, gradr_ad, mlt_vc_ad, Gamma, D, scale_height
         include 'formats'
         ierr = 0         
         
         if (s% use_Ledoux_criterion) then
            gradL_composition_term = s% gradL_composition_term(k)
         else
            gradL_composition_term = 0d0
         end if
         
         gradr_ad = get_gradr_face(s,k)
         grada_face_ad = get_grada_face_val(s,k)
         scale_height = get_scale_height_face(s,k)
         
         mixing_type = no_mixing
         call set_results
         
         if (k == 1 .and. s% mlt_make_surface_no_mixing) return
         if (s% lnT_start(k)/ln10 > s% max_logT_for_mlt) return         
         
         make_gradr_sticky_in_solver_iters = s% make_gradr_sticky_in_solver_iters
         if (.not. make_gradr_sticky_in_solver_iters .and. &
               s% min_logT_for_make_gradr_sticky_in_solver_iters < 1d20) then
            k_T_max = maxloc(s% lnT_start(1:s%nz),dim=1)
            make_gradr_sticky_in_solver_iters = &
               (s% lnT_start(k_T_max)/ln10 >= s% min_logT_for_make_gradr_sticky_in_solver_iters)
         end if
         
         if (make_gradr_sticky_in_solver_iters) then
            if (s% newer_fixed_gradr_for_rest_of_solver_iters(k)) then
               return
            end if
         end if

         if (s% no_MLT_below_shock .and. (s%u_flag .or. s%v_flag)) then 
            ! check for outward moving shock above k
            if (s% u_flag) then
               vel => s% u
            else
               vel => s% v
            end if
            do i=k-1,1,-1
               cs = s% csound(i)
               if (vel(i+1) >= cs .and. vel(i) < cs) return
            end do
         end if

         gradr_factor = s% gradr_factor(k)
         if (s% rotation_flag .and. s% mlt_use_rotation_correction) &
            gradr_factor = s% ft_rot(k)/s% fp_rot(k)*gradr_factor
         gradr_ad = gradr_factor*gradr_ad
         
         just_gradr = .false.
         call do1_mlt_eval_newer(s, k, s% MLT_option, just_gradr, gradL_composition_term, &
            gradr_ad, grada_face_ad, scale_height, mixing_length_alpha, &
            mixing_type, gradT_ad, Y_face_ad, mlt_vc_ad, D, Gamma, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'ierr in do1_mlt_newer_eval', k
            return
         end if
         
         call set_results

         contains
         
         subroutine set_results
            include 'formats'
            if (mixing_type == no_mixing) then
               gradT_ad = gradr_ad
               Y_face_ad = gradT_ad - grada_face_ad
               mlt_vc_ad = 0d0
               D = 0d0
               Gamma = 0d0
            end if         
            s% mlt_mixing_type(k) = mixing_type
            s% gradT_ad(k) = gradT_ad
            s% gradT(k) = s% gradT_ad(k)%val
            s% Y_face_ad(k) = Y_face_ad
            s% Y_face(k) = s% Y_face_ad(k)%val
            s% mlt_vc_ad(k) = mlt_vc_ad
            if (s% okay_to_set_mlt_vc) s% mlt_vc(k) = s% mlt_vc_ad(k)%val  
            s% mlt_D_ad(k) = D
            s% mlt_D(k) = D%val
            s% mlt_Gamma_ad(k) = Gamma
            s% mlt_Gamma(k) = Gamma%val
            s% gradr_ad(k) = gradr_ad
            s% gradr(k) = s% gradr_ad(k)%val
            s% grada_face_ad(k) = grada_face_ad
            s% grada_face(k) = grada_face_ad%val
            s% scale_height_ad(k) = scale_height
            s% scale_height(k) = scale_height%val             
            s% Lambda_ad(k) = mixing_length_alpha*scale_height
            s% mlt_mixing_length(k) = mixing_length_alpha*s% scale_height(k)    
            s% gradL_ad(k) = grada_face_ad + gradL_composition_term
            s% gradL(k) = s% gradL_ad(k)%val
            
            rho_face = get_rho_face(s,k)
            s% mlt_cdc(k) = s% mlt_D(k)*pow2(pi4*pow2(s%r(k))*rho_face%val)

            call adjust_gradT_fraction_newer(s, k, grada_face_ad, gradr_ad)         
            if (s% mlt_mixing_type(k) /= no_mixing .and. &
                s% mlt_option /= 'none' .and. abs(s% gradr(k)) > 0d0) then
               s% L_conv(k) = s% L(k) * (1d0 - s% gradT(k)/s% gradr(k)) ! C&G 14.109
            else
               s% L_conv(k) = 0d0
            end if         
            if (k == 1 .or. s% MLT_option == 'none') then
               s% grad_superad(k) = 0d0
            else
               s% grad_superad(k) = s% gradT(k) - s% grada_face(k)
               if (abs(s% lnPeos(k-1)-s% lnPeos(k)) < 1d-10) then
                  s% grad_superad_actual(k) = 0
               else
                  s% grad_superad_actual(k) = &
                     (s% lnT(k-1)-s% lnT(k))/(s% lnPeos(k-1)-s% lnPeos(k)) - s% grada_face(k)
               end if
            end if
            
         end subroutine set_results

      end subroutine do1_mlt_newer


      subroutine adjust_gradT_fraction_newer(s, k, grada_face_ad, gradr_ad)
         ! replace gradT by combo of grada_face and gradr
         ! then check excess
         use eos_def
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: grada_face_ad, gradr_ad
         include 'formats'         
         real(dp) :: f
         f = s% adjust_mlt_gradT_fraction(k)
         if (f >= 0d0 .and. f <= 1d0) then
            s% gradT_ad(k) = f*grada_face_ad + (1d0 - f)*gradr_ad
            s% gradT(k) = s% gradT_ad(k)%val
         end if
         s% gradT_sub_grada(k) = s% gradT(k) - s% grada_face(k) ! gradT_excess
         call adjust_gradT_excess_newer(s, k, grada_face_ad)
         s% gradT_sub_grada(k) = s% gradT(k) - s% grada_face(k) ! new gradT_excess from adjusted gradT
      end subroutine adjust_gradT_fraction_newer


      subroutine adjust_gradT_excess_newer(s, k, grada_face_ad)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: grada_face_ad
         real(dp) :: alfa, beta, gradT_excess_alpha, excess
         include 'formats'
         !s% gradT_excess_alpha is calculated at start of step and held constant during iterations
         ! gradT_excess_alpha = 0 means no efficiency boost; = 1 means full efficiency boost
         gradT_excess_alpha = s% gradT_excess_alpha
         s% gradT_excess_effect(k) = 0.0d0
         excess = s% gradT(k) - s% grada_face(k)
         if (gradT_excess_alpha <= 0.0 .or. excess <= s% gradT_excess_f1) return
         ! boost efficiency of energy transport
         alfa = s% gradT_excess_f2 ! for full boost, use this fraction of gradT
         if (gradT_excess_alpha < 1) then ! only partial boost, so increase alfa
            ! alfa goes to 1 as gradT_excess_alpha goes to 0
            ! alfa unchanged as gradT_excess_alpha goes to 1
            alfa = alfa + (1d0 - alfa)*(1d0 - gradT_excess_alpha)
         end if
         beta = 1.d0 - alfa
         s% gradT_ad(k) = alfa*s% gradT_ad(k) + beta*grada_face_ad
         s% gradT(k) = s% gradT_ad(k)%val
         s% gradT_excess_effect(k) = beta
      end subroutine adjust_gradT_excess_newer


      subroutine set_grads_newer(s, ierr)
         use chem_def, only: chem_isos
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


      subroutine switch_to_radiative_newer(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% mlt_mixing_type(k) = no_mixing
         s% mlt_mixing_length(k) = 0
         s% mlt_D_ad(k) = 0d0
         s% mlt_D(k) = 0d0
         s% mlt_cdc(k) = 0d0
         s% mlt_vc_ad(k) = 0d0
         s% mlt_vc(k) = 0d0
         s% mlt_Gamma_ad(k) = 0d0
         s% gradT_ad(k) = s% gradr_ad(k)
         s% gradT(k) = s% gradT_ad(k)%val
         s% Y_face_ad(k) = s% gradT_ad(k) - s% grada_face_ad(k)
      end subroutine switch_to_radiative_newer


      subroutine set_gradT_excess_alpha_newer(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         if (.not. s% okay_to_reduce_gradT_excess) then
            s% gradT_excess_alpha = 0
            return
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
                  call redo1_mlt_newer(s,kk,dr,op_err)
                  if (op_err /= 0) ierr = op_err
               end do
!$OMP END PARALLEL DO
            end if

         end subroutine end_of_convective_region


         subroutine redo1_mlt_newer(s,k,dr,ierr)
            type (star_info), pointer :: s
            integer, intent(in) :: k
            real(dp), intent(in) :: dr
            integer, intent(out) :: ierr
            real(dp) :: Hp, reduced_alpha
            logical :: make_gradr_sticky_in_solver_iters
            include 'formats'
            ierr = 0
            if (dr >= s% mlt_mixing_length(k)) return
            ! if convection zone is smaller than mixing length
            ! redo MLT with reduced alpha so mixing_length = dr
            Hp = s% scale_height(k)
            reduced_alpha = dr/Hp
            call do1_mlt_newer(s, k, reduced_alpha, make_gradr_sticky_in_solver_iters, ierr)
         end subroutine redo1_mlt_newer


      end subroutine check_for_redo_MLT_newer


      end module mlt_info_newer
