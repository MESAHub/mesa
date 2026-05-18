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

      use const_def, only: dp, convective_mixing, semiconvective_mixing
      use star_private_def
      use auto_diff_support

      implicit none

      private
      public :: set_Dmix_components

      contains

      ! If update_explicit_Dmix is true, this is a full mixing-info pass:
      ! refresh Dmix_explicit and Dmix_implicit from the current MESA state.
      ! If false, this is a solver iteration: leave Dmix_explicit fixed
      ! and update only Dmix_implicit from current MLT/TDC.
      ! In both cases, the set of implicit zones follows the
      ! post-cleanup mixing_type from the last full mixing-info pass;
      ! solver iterations update coefficients, not that zone selection.
      subroutine set_Dmix_components(s, update_explicit_Dmix)
         type (star_info), pointer :: s
         logical, intent(in) :: update_explicit_Dmix
         integer :: k
         real(dp) :: D_total, D_implicit
         logical :: use_implicit_Dmix

         if (s% job% implicit_diffusion_flag) then
            do k = 1, s% nz
               D_total = max(0d0, s% D_mix(k))
               use_implicit_Dmix = &
                  s% mixing_type(k) == convective_mixing .or. &
                  s% mixing_type(k) == semiconvective_mixing
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
            end do
         else
            do k = 1, s% nz
               s% Dmix_implicit(k) = 0d0
               s% Dmix_explicit(k) = s% D_mix(k)
            end do
         end if

      end subroutine set_Dmix_components

      end module implicit_Dmix
