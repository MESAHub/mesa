! ***********************************************************************
!
!   Copyright (C) 2015  Bill Paxton
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


      module hydro_alpha_rti_eqns

      use star_private_def
      use const_def

      implicit none

      private
      public :: do1_dalpha_RTI_dt_eqn


      contains


      subroutine do1_dalpha_RTI_dt_eqn(s, k, skip_partials, nvar, ierr)
         use star_utils, only: em1, e00, ep1
         use chem_def, only: ih1, ihe4

         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         integer, pointer :: reaction_id(:) ! maps net reaction number to reaction id
         integer :: nz, i_alpha_RTI, i_dalpha_RTI_dt, j, kk
         real(dp), pointer, dimension(:) :: sig
         logical :: okay
         real(dp) :: &
            dq, dm, dr, d_dadt_mix_dap1, d_dadt_mix_da00, d_dadt_mix_dam1, dap1, da00, &
            sig00, sigp1, dflux00_dam1, dflux00_da00, dflux00_dap1, &
            dfluxp1_da00, dfluxp1_dap1, r_max, r_min, &
            a00, am1, flux00, ap1, fluxp1, dadt_mix, &
            eqn_scale, dPdr_drhodr, instability2, instability, RTI_B, RTI_D, &
            A_plus_B_div_rho, source_plus, source_minus, ds, dadt_actual, dadt_expected, &
            r00, rp1, ravg_start, dP, drho, rho, rmid, cs, ds_da, &
            dadt_source, dequ_dE_const_Rho, dequ_dlnd, &
            dequ_dlnd_const_E, dequ_dlnPgas_const_T, dequ_dlnT, &
            dequ_dlnT_const_Pgas, dlnT_dlnd_const_E, fac, &
            d_dalpha_00, d_dalpha_m1, d_dalpha_p1
         logical :: test_partials

         include 'formats'

         ierr = 0
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         i_alpha_RTI = s% i_alpha_RTI
         i_dalpha_RTI_dt = s% i_dalpha_RTI_dt
         nz = s% nz

         sig => s% sig_RTI
         dq = s% dq(k)
         dm = s% dm(k)
         a00 = s% alpha_RTI(k)
         rho = s% rho_start(k)
         r00 = s% r_start(k)
         
         d_dadt_mix_dap1 = 0
         d_dadt_mix_dam1 = 0
         d_dadt_mix_da00 = 0
         dap1 = 0
         da00 = 0
         
         fac = s% alpha_RTI_diffusion_factor
          
         sig00 = fac*sig(k)
         if (k < nz) then
            sigp1 = fac*sig(k+1)
         else
            sigp1 = 0
         end if

         if (k > 1) then
            dflux00_dam1 = -sig00
            dflux00_da00 = sig00
         else
            dflux00_dam1 = 0
            dflux00_da00 = 0
         end if

         if (k < nz) then
            dfluxp1_da00 = -sigp1
            dfluxp1_dap1 = sigp1
         else
            dfluxp1_da00 = 0
            dfluxp1_dap1 = 0
         end if

         d_dadt_mix_da00 = 0
         d_dadt_mix_dam1 = 0
         flux00 = 0
         fluxp1 = 0
         dadt_mix = 0
         if (k > 1) then
            am1 = s% alpha_RTI(k-1)
            da00 = am1 - a00
            flux00 = -sig00*da00
         end if
         if (k < nz) then
            ap1 = s% alpha_RTI(k+1)
            dap1 = a00 - ap1
            fluxp1 = -sigp1*dap1
         end if
         dadt_mix = (fluxp1 - flux00)/dm
         d_dadt_mix_dap1 = dfluxp1_dap1/dm
         d_dadt_mix_da00 = (dfluxp1_da00 - dflux00_da00)/dm
         d_dadt_mix_dam1 = -dflux00_dam1/dm
         
         ds_da = 0
         dPdr_drhodr = s% dPdr_dRhodr_info(k)

         if (a00 <= 0d0 .or. s% RTI_D <= 0d0) then
            source_minus = 0d0
         else
            cs = s% csound_start(k)
            if (k < nz) then
               rmid = 0.5d0*(s% r_start(k) + s% r_start(k+1))
            else
               rmid = 0.5d0*(s% r_start(k) + s% R_center)
            end if
            RTI_D = s% RTI_D*max(1d0,a00/s% RTI_max_alpha)
            source_minus = RTI_D*a00*cs/rmid
            ds_da = -RTI_D*cs/rmid
         end if
         
         instability2 = -dPdr_drhodr ! > 0 means Rayleigh-Taylor unstable         
         if (instability2 <= 0d0 .or. &
               s% q(k) > s% alpha_RTI_src_max_q .or. &
               s% q(k) < s% alpha_RTI_src_min_q .or. &
               s% rho(k) < 1d99) then
               !s% rho(k) < s% alpha_RTI_src_min_rho) then
            source_plus = 0d0
            instability2 = 0d0
            instability = 0d0
            A_plus_B_div_rho = 0d0
         else
            RTI_B = s% RTI_B
            instability = sqrt(instability2)
            if (s% alpha_RTI_start(k) < s% RTI_max_alpha) then
               A_plus_B_div_rho = (s% RTI_A + RTI_B*a00)/rho
               ds_da = ds_da + RTI_B*instability/rho
            else ! turn off source when reach max
               A_plus_B_div_rho = 0d0
            end if
            source_plus = A_plus_B_div_rho*instability
         end if

         s% source_minus_alpha_RTI(k) = source_minus
         s% source_plus_alpha_RTI(k) = source_plus

         dadt_source = source_plus - source_minus

         dadt_expected = dadt_mix + dadt_source
         dadt_actual = s% dxh_alpha_RTI(k)/s% dt

         if (dadt_expected == 0d0) then
            s% equ(i_dalpha_RTI_dt,k) = dadt_expected
            eqn_scale = 1d0
         else
            eqn_scale = s% x_scale(i_dalpha_RTI_dt,k)/s% dt
            s% equ(i_dalpha_RTI_dt,k) = (dadt_expected - dadt_actual)/eqn_scale
         end if

         if (test_partials) then
            write(*,2) 'a00', k, a00
            write(*,2) 'dadt_mix', k, dadt_mix
            write(*,2) 'dadt_source', k, dadt_source
            write(*,2) 'source_plus', k, source_plus
            write(*,2) 'source_minus', k, source_minus
            write(*,2) 'A_plus_B_div_rho', k, A_plus_B_div_rho
            write(*,2) 'instability', k, instability
            write(*,2) 's% q(k)', k, s% q(k)
            write(*,2) 's% Peos(k)', k, s% Peos(k)
            write(*,2) 's% rho(k)', k, s% rho(k)
            write(*,2) 's% alpha_RTI_src_max_q', k, s% alpha_RTI_src_max_q
            write(*,2) 'dadt_expected', k, dadt_expected
            write(*,2) 'dadt_actual', k, dadt_actual
            write(*,2) 'eqn_scale', k, eqn_scale
            write(*,2) 's% dt', k, s% dt
            write(*,2) 'equ(i,k)', k, s% equ(i_dalpha_RTI_dt,k)
            s% solver_test_partials_val = s% equ(i_dalpha_RTI_dt,k)
         end if

         if (skip_partials) return
         
         d_dalpha_00 = 0
         d_dalpha_m1 = 0
         d_dalpha_p1 = 0
         
         if (dadt_expected == 0d0) then
            d_dalpha_00 = 1d0
         else
            ! partial of -dadt_actual/eqn_scale
            d_dalpha_00 = -1d0/eqn_scale/s% dt

            if (dadt_mix /= 0d0) then ! partials of dadt_mix/eqn_scale
               d_dalpha_00 = d_dalpha_00 + d_dadt_mix_da00/eqn_scale
               if (k > 1) then
                  d_dalpha_m1 = d_dadt_mix_dam1/eqn_scale
                  call em1(s, i_dalpha_RTI_dt, i_alpha_RTI, k, nvar, d_dalpha_m1)
               end if
               if (k < nz) then
                  d_dalpha_p1 = d_dadt_mix_dap1/eqn_scale
                  call ep1(s, i_dalpha_RTI_dt, i_alpha_RTI, k, nvar, d_dalpha_p1)
               end if
            end if

            ! partials of dadt_source/eqn_scale
            if (dadt_source /= 0d0) d_dalpha_00 = d_dalpha_00 + ds_da/eqn_scale
         
         end if

         
         call e00(s, i_dalpha_RTI_dt, i_alpha_RTI, k, nvar, d_dalpha_00)

         if (test_partials) then   
            s% solver_test_partials_var = i_alpha_RTI
            s% solver_test_partials_dval_dx = d_dalpha_00
            write(*,*) 'do1_dalpha_RTI_dt_eqn', s% solver_test_partials_var
         end if

      end subroutine do1_dalpha_RTI_dt_eqn


      end module hydro_alpha_rti_eqns

