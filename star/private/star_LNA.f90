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

      module star_lna

      use star_private_def
      use auto_diff
      use hydro_vars, only: set_vars_if_needed
      use hydro_riemann, only: do_uface_and_Pface
      use star_lna_turbulence_closures, only: tdc_lna_active, tdc_conv_vel_for_star_LNA
      use star_lna_support, only: &
         star_LNA_problem, check_star_LNA_model, setup_star_LNA_problem, &
         report_star_LNA_setup, &
         assemble_density_rows, assemble_radius_rows, assemble_momentum_rows, &
         assemble_energy_rows, assemble_luminosity_rows, &
         assemble_rsp2_turbulent_rows, assemble_tdc_turbulent_rows, &
         write_star_LNA_matrix_summary, solve_dense_star_LNA, &
         free_star_LNA_problem, &
         support_star_LNA_L_conv_ad => star_LNA_L_conv_ad

      implicit none

      private
      public :: do_star_LNA
      public :: star_LNA_L_conv_ad

      contains

      subroutine do_star_LNA(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         type(star_LNA_problem) :: problem
         type(auto_diff_real_star_order1) :: mlt_vc_ad(s% nz)
         integer :: k, op_err
         logical :: Riemann_Pturb

         ierr = 0

         call set_vars_if_needed(s, 0d0, 'star_LNA', ierr)
         if (ierr /= 0) return

         call check_star_LNA_model(s, ierr)
         if (ierr /= 0) return

         Riemann_Pturb = s% u_flag .and. s% mlt_Pturb_factor > 0d0
         if (s% u_flag) then
            ! Establish the nonlinear face and start state before the LNA derivatives.
            call do_uface_and_Pface(s, ierr)
            if (ierr /= 0) return
            if (Riemann_Pturb) then
               do k = 1, s% nz
                  mlt_vc_ad(k) = s% mlt_vc_old(k)
                  if (tdc_lna_active(s)) call tdc_conv_vel_for_star_LNA(s, k, mlt_vc_ad(k))
               end do
               call do_uface_and_Pface(s, ierr, mlt_vc_ad)
            end if
            if (ierr /= 0) goto 100
         end if

         call setup_star_LNA_problem(s, problem, ierr)
         if (ierr /= 0) goto 100

         call report_star_LNA_setup(s, problem% map)

         call assemble_star_LNA_equations(s, problem, ierr)
         if (ierr == 0 .and. s% star_LNA_write_matrix_summary) &
            call write_star_LNA_matrix_summary(s, problem, ierr)
         if (ierr == 0) call solve_dense_star_LNA(s, problem% map, problem% mtx, ierr)

         call free_star_LNA_problem(problem)

100      continue
         if (Riemann_Pturb) then
            ! Restore the nonlinear face derivatives.
            call do_uface_and_Pface(s, op_err)
            if (ierr == 0) ierr = op_err
         end if
      end subroutine do_star_LNA


      subroutine assemble_star_LNA_equations(s, problem, ierr)
         type(star_info), pointer :: s
         type(star_LNA_problem), intent(inout) :: problem
         integer, intent(out) :: ierr

         ierr = 0

         ! rho_k*DeltaV_k = dm_k.
         call assemble_density_rows(s, problem% map, problem% mtx, ierr)
         if (ierr /= 0) return

         ! d lnR_k/dt = v_face,k/r_k.
         call assemble_radius_rows(s, problem% map, problem% mtx, ierr)
         if (ierr /= 0) return

         ! d velocity_k/dt = pressure, gravity, and turbulent-stress acceleration.
         call assemble_momentum_rows(s, problem% map, problem% mtx, ierr)
         if (ierr /= 0) return

         ! dE_k/dt = sources - dL/dm - pressure work.
         call assemble_energy_rows(s, problem% map, problem% mtx, ierr)
         if (ierr /= 0) return

         ! Radiative/convective luminosity closure and surface luminosity BC.
         call assemble_luminosity_rows(s, problem% map, problem% mtx, ierr)
         if (ierr /= 0) return

         ! RSP2 turbulent-energy and independent luminosity closures.
         call assemble_rsp2_turbulent_rows(s, problem% map, problem% mtx, ierr)
         if (ierr /= 0) return

         ! TDC convective-velocity closure.
         call assemble_tdc_turbulent_rows(s, problem% map, problem% mtx, ierr)
         if (ierr /= 0) return
      end subroutine assemble_star_LNA_equations


      function star_LNA_L_conv_ad(s, k) result(L_conv_ad)
         type(star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: L_conv_ad

         L_conv_ad = support_star_LNA_L_conv_ad(s, k)
      end function star_LNA_L_conv_ad

      end module star_lna
