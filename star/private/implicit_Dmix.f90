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

      use const_def, only: dp, no_mixing, convective_mixing, rotation_mixing, &
         semiconvective_mixing, thermohaline_mixing
      use star_private_def
      use auto_diff_support

      implicit none

      private
      public :: set_Dmix_components

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

         logical function is_implicit_mixing_type(mixing_type)
            integer, intent(in) :: mixing_type

            is_implicit_mixing_type = &
               mixing_type == convective_mixing .or. &
               mixing_type == semiconvective_mixing .or. &
               mixing_type == thermohaline_mixing
         end function is_implicit_mixing_type

         logical function rotation_mixing_dominates()
            rotation_mixing_dominates = s% rotation_flag .and. s% D_mix(k) /= 0d0 .and. &
               max(0d0, s% D_mix(k) - s% D_mix_rotation(k)) < s% D_mix_rotation(k)
         end function rotation_mixing_dominates

      end subroutine set_Dmix_components

      end module implicit_Dmix
