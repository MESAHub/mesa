! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
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

      module implicit_Dmix

      use const_def, only: dp, ln10, no_mixing, convective_mixing, &
         rotation_mixing, semiconvective_mixing, thermohaline_mixing
      use star_private_def
      use auto_diff_support

      implicit none

      private
      public :: set_Dmix_components, &
         implicit_Dmix_debug_mlt, &
         implicit_Dmix_debug_thermohaline, &
         implicit_Dmix_debug_sigma, &
         implicit_Dmix_debug_sigma_limiter, &
         implicit_Dmix_debug_profile, &
         implicit_Dmix_debug_summary, &
         implicit_Dmix_debug_selected_cell

      contains

      ! If update_explicit_Dmix is true, this is a full mixing-info pass:
      ! refresh Dmix_explicit and Dmix_implicit from the current MESA state.
      ! If false, this is a solver iteration: leave Dmix_explicit fixed
      ! and update only Dmix_implicit from current local turbulence.
      ! Full passes use the post-cleanup mixing_type.  Solver iterations
      ! use current mlt_mixing_type for local promoted transport and update
      ! matching promoted display flags; nonlocal cleanup remains full-pass.
      subroutine set_Dmix_components(s, update_explicit_Dmix)
         type (star_info), pointer :: s
         logical, intent(in) :: update_explicit_Dmix
         integer :: k, mixing_type_for_implicit
         real(dp) :: D_total, D_implicit
         logical :: use_implicit_Dmix

         if (s% job% implicit_diffusion_flag) then
            do k = 1, s% nz
               D_total = max(0d0, s% D_mix(k))
               if (update_explicit_Dmix) then
                  mixing_type_for_implicit = s% mixing_type(k)
               else
                  mixing_type_for_implicit = s% mlt_mixing_type(k)
               end if
               use_implicit_Dmix = is_implicit_mixing_type(mixing_type_for_implicit)
               if (use_implicit_Dmix) then
                  s% Dmix_implicit(k) = s% mlt_D_ad(k)
                  if (s% Dmix_implicit(k)% val <= 0d0) then
                     s% Dmix_implicit(k) = 0d0
                     use_implicit_Dmix = .false.
                  end if
               else
                  s% Dmix_implicit(k) = 0d0
               end if

               if (update_explicit_Dmix) then
                  if (use_implicit_Dmix) then
                     D_implicit = s% Dmix_implicit(k)% val
                     if (D_implicit > D_total) then
                        ! Later full-pass processing can reduce D_mix; keep the split non-negative.
                        if (D_total == 0d0) then
                           s% Dmix_implicit(k) = 0d0
                        else
                           s% Dmix_implicit(k) = (D_total/D_implicit)*s% Dmix_implicit(k)
                        end if
                        D_implicit = D_total
                     end if
                     ! Full mixing-info pass: preserve MESA's total D_mix.
                     ! Dmix_explicit is the non-implicit part of that total.
                     s% Dmix_explicit(k) = D_total - D_implicit
                  else
                     ! Full mixing-info pass: non-promoted transport keeps the
                     ! ordinary MESA coefficient in Dmix_explicit.
                     s% Dmix_explicit(k) = D_total
                  end if
               end if

               s% D_mix(k) = s% Dmix_implicit(k)% val + s% Dmix_explicit(k)
               if (.not. update_explicit_Dmix) call set_solver_iter_mixing_type
               if (s% rotation_flag) then
                  s% D_mix_non_rotation(k) = max(0d0, s% D_mix(k) - s% D_mix_rotation(k))
               else
                  s% D_mix_non_rotation(k) = s% D_mix(k)
               end if
               if (implicit_Dmix_debug_cell( &
                     s, k, mixing_type_for_implicit, use_implicit_Dmix, D_total)) &
                  call write_Dmix_split_debug(s, update_explicit_Dmix, k, mixing_type_for_implicit, &
                     use_implicit_Dmix, D_total)
            end do
         else
            do k = 1, s% nz
               s% Dmix_implicit(k) = 0d0
               s% Dmix_explicit(k) = s% D_mix(k)
            end do
         end if

         contains

         subroutine set_solver_iter_mixing_type
            if (use_implicit_Dmix) then
               s% mixing_type(k) = mixing_type_for_implicit
            else if (is_implicit_mixing_type(s% mixing_type(k))) then
               if (is_implicit_mixing_type(mixing_type_for_implicit)) then
                  s% mixing_type(k) = no_mixing
               else
                  s% mixing_type(k) = mixing_type_for_implicit
               end if
               if (rotation_mixing_dominates()) then
                  s% mixing_type(k) = rotation_mixing
               end if
            end if
         end subroutine set_solver_iter_mixing_type

         logical function rotation_mixing_dominates()
            rotation_mixing_dominates = s% rotation_flag .and. s% D_mix(k) /= 0d0 .and. &
               max(0d0, s% D_mix(k) - s% D_mix_rotation(k)) < s% D_mix_rotation(k)
         end function rotation_mixing_dominates

      end subroutine set_Dmix_components


      logical function is_implicit_mixing_type(mixing_type)
         integer, intent(in) :: mixing_type

         is_implicit_mixing_type = &
            mixing_type == convective_mixing .or. &
            mixing_type == semiconvective_mixing .or. &
            mixing_type == thermohaline_mixing
      end function is_implicit_mixing_type


      subroutine implicit_Dmix_debug_mlt(s, stage, k, mixing_type, D, &
            gradL_composition_term, gradr, grada, gradT, Y_face, conv_vel)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: stage
         integer, intent(in) :: k, mixing_type
         real(dp), intent(in) :: D, gradL_composition_term, gradr, grada, &
            gradT, Y_face, conv_vel

         if (.not. implicit_Dmix_debug_cell(s, k, mixing_type, &
               is_implicit_mixing_type(mixing_type), D)) return

         write(*,*) 'implicit_Dmix debug mlt ', trim(stage), &
            ' model', s% model_number, ' iter', s% solver_iter, &
            ' k', k, ' mlt_type', mixing_type, &
            ' current_mix_type', s% mixing_type(k), &
            ' tdc_num_iters', s% tdc_num_iters(k)
         write(*,*) '   D_mlt', D, ' Dmix_implicit', s% Dmix_implicit(k)%val, &
            ' Dmix_explicit', s% Dmix_explicit(k), &
            ' D_mix', s% D_mix(k), ' D_nonrot', s% D_mix_non_rotation(k)
         write(*,*) '   gradL_comp', gradL_composition_term, &
            ' brunt_B', s% brunt_B(k), ' gradr', gradr, ' grada', grada, &
            ' gradT', gradT, ' Y_face', Y_face, ' conv_vel', conv_vel
         write(*,*) '   gradr_minus_gradL', gradr - (grada + gradL_composition_term), &
            ' gradr_minus_grada', gradr - grada
      end subroutine implicit_Dmix_debug_mlt


      subroutine implicit_Dmix_debug_thermohaline(s, stage, k, mixing_type, &
            D, gradL_composition_term, gradr, grada, gradL, gradT, Y_face, conv_vel)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: stage
         integer, intent(in) :: k, mixing_type
         real(dp), intent(in) :: D, gradL_composition_term, gradr, grada, gradL, &
            gradT, Y_face, conv_vel

         if (.not. implicit_Dmix_debug_cell(s, k, mixing_type, &
               is_implicit_mixing_type(mixing_type), D)) return

         write(*,*) 'implicit_Dmix debug thermohaline ', trim(stage), &
            ' model', s% model_number, ' iter', s% solver_iter, &
            ' k', k, ' proposed_type', mixing_type, &
            ' mlt_type', s% mlt_mixing_type(k), &
            ' current_mix_type', s% mixing_type(k), &
            ' tdc_num_iters', s% tdc_num_iters(k)
         write(*,*) '   D_eval', D, ' mlt_D', s% mlt_D(k), &
            ' Dmix_implicit', s% Dmix_implicit(k)%val, &
            ' Dmix_explicit', s% Dmix_explicit(k), ' D_mix', s% D_mix(k)
         write(*,*) '   gradL_comp', gradL_composition_term, &
            ' brunt_B', s% brunt_B(k), ' gradr', gradr, ' grada', grada, &
            ' gradL', gradL, ' gradT', gradT, ' Y_face', Y_face, &
            ' conv_vel', conv_vel
         write(*,*) '   gradr_minus_gradL', gradr - gradL, &
            ' gradr_minus_grada', gradr - grada
      end subroutine implicit_Dmix_debug_thermohaline


      subroutine implicit_Dmix_debug_sigma(s, k, sig, dD_dB, dsig_factor, &
            do_implicit_dsig_structure, do_implicit_dsig_dxa)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: sig, dD_dB, dsig_factor
         logical, intent(in) :: do_implicit_dsig_structure, do_implicit_dsig_dxa
         real(dp) :: sig_implicit_val, Dimp_d1_max, dsig_struct_max, dB_dxa_m1_max, &
            dB_dxa_00_max, dsig_dxa_m1_max, dsig_dxa_00_max

         if (.not. implicit_Dmix_debug_cell(s, k, s% mixing_type(k), &
               is_implicit_mixing_type(s% mixing_type(k)), s% D_mix(k))) return

         sig_implicit_val = -99d0
         Dimp_d1_max = maxval(abs(s% Dmix_implicit(k)%d1Array))
         dsig_struct_max = 0d0
         if (do_implicit_dsig_structure) then
            sig_implicit_val = s% sig_implicit_ad(k)%val
            dsig_struct_max = maxval(abs(s% sig_implicit_ad(k)%d1Array))
         end if
         dB_dxa_m1_max = 0d0
         dB_dxa_00_max = 0d0
         dsig_dxa_m1_max = 0d0
         dsig_dxa_00_max = 0d0
         if (do_implicit_dsig_dxa) then
            dB_dxa_m1_max = maxval(abs(s% d_brunt_B_dxa_m1(:,k)))
            dB_dxa_00_max = maxval(abs(s% d_brunt_B_dxa_00(:,k)))
            dsig_dxa_m1_max = maxval(abs(s% d_sig_dxa_m1(:,k)))
            dsig_dxa_00_max = maxval(abs(s% d_sig_dxa_00(:,k)))
         end if
         write(*,*) 'implicit_Dmix debug sigma model', s% model_number, &
            ' iter', s% solver_iter, ' k', k, &
            ' mix_type', s% mixing_type(k), ' mlt_type', s% mlt_mixing_type(k), &
            ' do_structure', do_implicit_dsig_structure, &
            ' do_dxa', do_implicit_dsig_dxa
         write(*,*) '   sig', sig, ' sig_implicit_val', sig_implicit_val, &
            ' Dmix_implicit', s% Dmix_implicit(k)%val, &
            ' Dmix_explicit', s% Dmix_explicit(k), ' D_mix', s% D_mix(k)
         write(*,*) '   Dimp_d1_max', Dimp_d1_max, &
            ' dsig_struct_max', dsig_struct_max, &
            ' dD_dB', dD_dB, ' dsig_factor', dsig_factor, &
            ' dB_dxa_m1_max', dB_dxa_m1_max, &
            ' dB_dxa_00_max', dB_dxa_00_max, &
            ' dsig_dxa_m1_max', dsig_dxa_m1_max, &
            ' dsig_dxa_00_max', dsig_dxa_00_max
         end subroutine implicit_Dmix_debug_sigma


         subroutine implicit_Dmix_debug_sigma_limiter(s, stage, k, &
               sig_before, sig_after, siglim, lim, limiter_applied, &
               derivative_action, do_implicit_dsig_structure, do_implicit_dsig_dxa, &
               jlim, Xlim, dXlim, dmlim, delta_m_to_bdy)
            type (star_info), pointer :: s
            character (len=*), intent(in) :: stage, derivative_action
            integer, intent(in) :: k, jlim
            real(dp), intent(in) :: sig_before, sig_after, siglim, lim, &
               Xlim, dXlim, dmlim, delta_m_to_bdy
            logical, intent(in) :: limiter_applied, do_implicit_dsig_structure, &
               do_implicit_dsig_dxa
            integer :: j, j_min_m1, j_min_00
            real(dp) :: sig_implicit_val, sig_ratio, Dimp_d1_max, dsig_struct_max, &
               dsig_dxa_m1_max, dsig_dxa_00_max, dm_face, dm_m1, dq_m1, &
               min_xa_m1, min_xa_00, xa_jlim_m1, xa_jlim_00

            if (.not. implicit_Dmix_debug_selected_cell(s, k)) return

            sig_implicit_val = -99d0
            Dimp_d1_max = maxval(abs(s% Dmix_implicit(k)%d1Array))
            dsig_struct_max = 0d0
            if (do_implicit_dsig_structure) then
               sig_implicit_val = s% sig_implicit_ad(k)%val
               dsig_struct_max = maxval(abs(s% sig_implicit_ad(k)%d1Array))
            end if
            dsig_dxa_m1_max = 0d0
            dsig_dxa_00_max = 0d0
            if (do_implicit_dsig_dxa) then
               dsig_dxa_m1_max = maxval(abs(s% d_sig_dxa_m1(:,k)))
               dsig_dxa_00_max = maxval(abs(s% d_sig_dxa_00(:,k)))
            end if
            sig_ratio = 0d0
            if (sig_before /= 0d0) sig_ratio = sig_after/sig_before

            dm_face = -99d0
            dm_m1 = -99d0
            dq_m1 = -99d0
            min_xa_m1 = 1d99
            min_xa_00 = 1d99
            j_min_m1 = 0
            j_min_00 = 0
            xa_jlim_m1 = -99d0
            xa_jlim_00 = -99d0
            if (k > 1) then
               dm_m1 = s% dm(k-1)
               dq_m1 = s% dq(k-1)
               dm_face = 0.5d0*(s% dm(k-1) + s% dm(k))
               do j = 1, s% species
                  if (s% xa(j,k-1) < min_xa_m1) then
                     min_xa_m1 = s% xa(j,k-1)
                     j_min_m1 = j
                  end if
               end do
            end if
            do j = 1, s% species
               if (s% xa(j,k) < min_xa_00) then
                  min_xa_00 = s% xa(j,k)
                  j_min_00 = j
               end if
            end do
            if (jlim > 0 .and. jlim <= s% species) then
               if (k > 1) xa_jlim_m1 = s% xa(jlim,k-1)
               xa_jlim_00 = s% xa(jlim,k)
            end if

            write(*,*) 'implicit_Dmix debug sigma limiter ', trim(stage), &
               ' model', s% model_number, ' iter', s% solver_iter, &
               ' k', k, ' mix_type', s% mixing_type(k), &
               ' mlt_type', s% mlt_mixing_type(k), &
               ' applied', limiter_applied, &
               ' derivative_action', trim(derivative_action), &
               ' do_structure', do_implicit_dsig_structure, &
               ' do_dxa', do_implicit_dsig_dxa
            write(*,*) '   sig_raw', s% sig_raw(k), &
               ' sig_before', sig_before, ' sig_after', sig_after, &
               ' siglim', siglim, ' lim', lim, ' sig_ratio', sig_ratio, &
               ' sig_implicit_val', sig_implicit_val
            write(*,*) '   Dmix_implicit', s% Dmix_implicit(k)%val, &
               ' Dmix_explicit', s% Dmix_explicit(k), ' D_mix', s% D_mix(k), &
               ' Dimp_d1_max', Dimp_d1_max, &
               ' dsig_struct_max', dsig_struct_max, &
               ' dsig_dxa_m1_max', dsig_dxa_m1_max, &
               ' dsig_dxa_00_max', dsig_dxa_00_max
            write(*,*) '   dt', s% dt, ' logT', s% lnT(k)/ln10, &
               ' dm_face', dm_face, ' dm_m1', dm_m1, &
               ' dm_00', s% dm(k), ' dq_m1', dq_m1, &
               ' dq_00', s% dq(k)
            write(*,*) '   jlim', jlim, ' Xlim', Xlim, ' dXlim', dXlim, &
               ' dmlim', dmlim, ' delta_m_to_bdy', delta_m_to_bdy, &
               ' xa_jlim_m1', xa_jlim_m1, ' xa_jlim_00', xa_jlim_00, &
               ' j_min_m1', j_min_m1, ' min_xa_m1', min_xa_m1, &
               ' j_min_00', j_min_00, ' min_xa_00', min_xa_00
         end subroutine implicit_Dmix_debug_sigma_limiter


      subroutine implicit_Dmix_debug_profile(s, column_name, k, val)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: column_name
         integer, intent(in) :: k
         real(dp), intent(in) :: val

         if (.not. implicit_Dmix_debug_cell(s, k, s% mixing_type(k), &
               is_implicit_mixing_type(s% mixing_type(k)), s% D_mix(k))) return

         write(*,*) 'implicit_Dmix debug profile ', trim(column_name), &
            ' model', s% model_number, ' iter', s% solver_iter, &
            ' k', k, ' val', val, &
            ' mix_type', s% mixing_type(k), ' mlt_type', s% mlt_mixing_type(k)
         write(*,*) '   Dmix_implicit', s% Dmix_implicit(k)%val, &
            ' Dmix_explicit', s% Dmix_explicit(k), &
            ' D_mix', s% D_mix(k), ' D_nonrot', s% D_mix_non_rotation(k), &
            ' mlt_D', s% mlt_D(k)
      end subroutine implicit_Dmix_debug_profile


      subroutine implicit_Dmix_debug_summary(s, stage)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: stage
         integer :: k, n_mlt_th, n_mix_th, n_implicit, k_max_implicit
         real(dp) :: max_implicit

         if (.not. s% job% implicit_diffusion_debug_thermohaline) return

         n_mlt_th = 0
         n_mix_th = 0
         n_implicit = 0
         k_max_implicit = 0
         max_implicit = 0d0
         do k = 1, s% nz
            if (s% mlt_mixing_type(k) == thermohaline_mixing) n_mlt_th = n_mlt_th + 1
            if (s% mixing_type(k) == thermohaline_mixing) n_mix_th = n_mix_th + 1
            if (s% Dmix_implicit(k)%val > 0d0) then
               n_implicit = n_implicit + 1
               if (s% Dmix_implicit(k)%val > max_implicit) then
                  max_implicit = s% Dmix_implicit(k)%val
                  k_max_implicit = k
               end if
            end if
         end do

         write(*,*) 'implicit_Dmix debug summary ', trim(stage), &
            ' model', s% model_number, ' iter', s% solver_iter, &
            ' n_mlt_th', n_mlt_th, ' n_mix_th', n_mix_th, &
            ' n_implicit', n_implicit, ' k_max_implicit', k_max_implicit, &
            ' max_Dmix_implicit', max_implicit

         do k = 1, s% nz
            if (.not. implicit_Dmix_debug_cell(s, k, s% mixing_type(k), &
                  is_implicit_mixing_type(s% mixing_type(k)), s% D_mix(k))) cycle
            write(*,*) '   cell k', k, ' q', s% q(k), &
               ' mlt_type', s% mlt_mixing_type(k), &
               ' mix_type', s% mixing_type(k), &
               ' gradL_comp', s% gradL_composition_term(k), &
               ' brunt_B', s% brunt_B(k), &
               ' mlt_D', s% mlt_D(k), &
               ' Dimp', s% Dmix_implicit(k)%val, &
               ' Dexp', s% Dmix_explicit(k), &
               ' Dmix', s% D_mix(k), &
               ' Dnonrot', s% D_mix_non_rotation(k)
         end do
      end subroutine implicit_Dmix_debug_summary


      logical function implicit_Dmix_debug_cell(s, k, mixing_type_for_implicit, &
            use_implicit_Dmix, D_total)
         type (star_info), pointer :: s
         integer, intent(in) :: k, mixing_type_for_implicit
         logical, intent(in) :: use_implicit_Dmix
         real(dp), intent(in) :: D_total
         integer :: k_select, k_min, k_max

         implicit_Dmix_debug_cell = .false.
         if (.not. s% job% implicit_diffusion_debug_thermohaline) return

         k_select = s% job% implicit_diffusion_debug_thermohaline_k
         if (k_select > 0) then
            implicit_Dmix_debug_cell = (k == k_select)
         else if (k_select < 0) then
            implicit_Dmix_debug_cell = .true.
         else if (s% job% implicit_diffusion_debug_thermohaline_k_min > 0 .or. &
               s% job% implicit_diffusion_debug_thermohaline_k_max > 0) then
            k_min = s% job% implicit_diffusion_debug_thermohaline_k_min
            k_max = s% job% implicit_diffusion_debug_thermohaline_k_max
            if (k_min <= 0) k_min = 1
            if (k_max <= 0) k_max = s% nz
            if (k_min > k_max) then
               k_select = k_min
               k_min = k_max
               k_max = k_select
            end if
            implicit_Dmix_debug_cell = (k >= k_min .and. k <= k_max)
         else
            implicit_Dmix_debug_cell = &
               mixing_type_for_implicit == thermohaline_mixing .or. &
               s% mlt_mixing_type(k) == thermohaline_mixing .or. &
               s% mixing_type(k) == thermohaline_mixing .or. &
               s% gradL_composition_term(k) < 0d0 .or. &
               use_implicit_Dmix .or. &
               s% Dmix_implicit(k)%val > 0d0 .or. &
               s% mlt_D(k) > 0d0 .or. &
               D_total > 0d0
         end if
      end function implicit_Dmix_debug_cell


      logical function implicit_Dmix_debug_selected_cell(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k

         implicit_Dmix_debug_selected_cell = implicit_Dmix_debug_cell( &
            s, k, s% mixing_type(k), &
            is_implicit_mixing_type(s% mixing_type(k)), s% D_mix(k))
      end function implicit_Dmix_debug_selected_cell


      subroutine write_Dmix_split_debug(s, update_explicit_Dmix, k, &
            mixing_type_for_implicit, use_implicit_Dmix, D_total_before_split)
         type (star_info), pointer :: s
         logical, intent(in) :: update_explicit_Dmix, use_implicit_Dmix
         integer, intent(in) :: k, mixing_type_for_implicit
         real(dp), intent(in) :: D_total_before_split

         write(*,*) 'implicit_Dmix debug split model', s% model_number, &
            ' iter', s% solver_iter, ' k', k, &
            ' update_explicit', update_explicit_Dmix, &
            ' use_implicit', use_implicit_Dmix, &
            ' implicit_type', mixing_type_for_implicit
         write(*,*) '   mlt_type', s% mlt_mixing_type(k), &
            ' mix_type', s% mixing_type(k), &
            ' D_total_before', D_total_before_split, &
            ' mlt_D_ad_val', s% mlt_D_ad(k)%val, &
            ' Dmix_implicit', s% Dmix_implicit(k)%val, &
            ' Dmix_explicit', s% Dmix_explicit(k), &
            ' D_mix', s% D_mix(k), &
            ' D_nonrot', s% D_mix_non_rotation(k)
         write(*,*) '   dDimp_dB', s% Dmix_implicit(k)%d1Array(i_xtra2_00), &
            ' gradL_comp', s% gradL_composition_term(k), &
            ' brunt_B', s% brunt_B(k), &
            ' gradr', s% gradr(k), ' grada_face', s% grada_face(k), &
            ' gradT', s% gradT(k)
      end subroutine write_Dmix_split_debug

      end module implicit_Dmix
