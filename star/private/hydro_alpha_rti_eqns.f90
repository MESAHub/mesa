! ***********************************************************************
!
!   Copyright (C) 2015  The MESA Team
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
      use auto_diff_support

      implicit none

      private
      public :: do1_dalpha_RTI_dt_eqn


      contains


      subroutine do1_dalpha_RTI_dt_eqn(s, k, nvar, ierr)
         use star_utils, only: em1, e00, ep1
         use chem_def, only: ih1, ihe4

         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         integer, intent(out) :: ierr

         type(auto_diff_real_4var_order1) :: a_m1, a_00, a_p1
         type(auto_diff_real_4var_order1) :: da00, flux00, dap1, fluxp1, RTI_D, A_plus_B_div_rho
         type(auto_diff_real_4var_order1) :: source_minus, source_plus, dadt_source, dadt_actual
         type(auto_diff_real_4var_order1) :: dadt_mix, dadt_expected, resid

         integer :: nz, i_alpha_RTI, i_dalpha_RTI_dt, j, kk
         real(dp), pointer, dimension(:) :: sig
         logical :: okay
         real(dp) :: &
            dq, dm, dr, d_dadt_mix_dap1, sig00, sigp1, &
            eqn_scale, dPdr_drhodr, instability2, instability, RTI_B, &
            r00, rp1, drho, rho, rmid, cs, fac
         logical :: test_partials

         include 'formats'

         ! Debugging setup
         ierr = 0
         !test_partials = (k == s% ctrl% solver_test_partials_k)
         test_partials = .false.

         ! Equation setup
         i_alpha_RTI = s% i_alpha_RTI
         i_dalpha_RTI_dt = s% i_dalpha_RTI_dt
         nz = s% nz

         ! Starting/fixed quantities
         sig => s% sig_RTI
         dq = s% dq(k)
         dm = s% dm(k)
         rho = s% rho_start(k)
         r00 = s% r_start(k)        
         fac = s% ctrl% alpha_RTI_diffusion_factor
          
         sig00 = fac*sig(k)
         if (k < nz) then
            sigp1 = fac*sig(k+1)
         else
            sigp1 = 0
         end if

         a_00 = s% alpha_RTI(k)
         a_00%d1val2 = 1d0

         ! alpha's and fluxes
         if (k > 1) then
            a_m1 = s% alpha_RTI(k-1)
            a_m1%d1val1 = 1d0

            da00 = a_m1 - a_00
            flux00 = -sig00*da00
         else
            a_m1 = 0d0
            flux00 = 0d0
         end if

         if (k < nz) then
            a_p1 = s% alpha_RTI(k+1)
            a_p1%d1val3 = 1d0

            dap1 = a_00 - a_p1
            fluxp1 = -sigp1*dap1
         else
            a_p1 = 0d0
            fluxp1 = 0d0
         end if

         ! Flux divergence
         dadt_mix = (fluxp1 - flux00)/dm

         ! Sources and sink s        
         dPdr_drhodr = s% dPdr_dRhodr_info(k)

         if (a_00 <= 0d0 .or. s% ctrl% RTI_D <= 0d0) then
            source_minus = 0d0
         else
            cs = s% csound_start(k)
            if (k < nz) then
               rmid = 0.5d0*(s% r_start(k) + s% r_start(k+1))
            else
               rmid = 0.5d0*(s% r_start(k) + s% R_center)
            end if
            RTI_D = s% ctrl% RTI_D*max(1d0,a_00/s% ctrl% RTI_max_alpha)
            source_minus = RTI_D*a_00*cs/rmid
         end if
         
         instability2 = -dPdr_drhodr ! > 0 means Rayleigh-Taylor unstable         
         if (instability2 <= 0d0 .or. &
               s% q(k) > s% ctrl% alpha_RTI_src_max_q .or. &
               s% q(k) < s% ctrl% alpha_RTI_src_min_q .or. &
               s% rho(k) < 1d99) then
            source_plus = 0d0
            instability2 = 0d0
            instability = 0d0
            A_plus_B_div_rho = 0d0
         else
            RTI_B = s% ctrl% RTI_B
            instability = sqrt(instability2)
            if (s% alpha_RTI_start(k) < s% ctrl% RTI_max_alpha) then
               A_plus_B_div_rho = (s% ctrl% RTI_A + RTI_B*a_00)/rho
            else ! turn off source when reach max
               A_plus_B_div_rho = 0d0
            end if
            source_plus = A_plus_B_div_rho*instability
         end if

         s% source_minus_alpha_RTI(k) = source_minus%val
         s% source_plus_alpha_RTI(k) = source_plus%val

         dadt_source = source_plus - source_minus

         dadt_expected = dadt_mix + dadt_source

         ! We use dxh to avoid subtraction errors.
         ! At least that's what I assume. I just preserved this
         ! choice when re-writing it... - Adam S. Jermyn 6/15/2021
         dadt_actual = s% dxh_alpha_RTI(k)/s% dt
         dadt_actual%d1val2 = 1d0 / s% dt

         eqn_scale = max(1d0, s% x_scale(i_dalpha_RTI_dt,k)/s% dt)
         resid = (dadt_expected - dadt_actual)/eqn_scale

         ! Unpack
         s% equ(i_dalpha_RTI_dt,k) = resid%val
         call em1(s, i_dalpha_RTI_dt, i_alpha_RTI, k, nvar, resid%d1val1)
         call e00(s, i_dalpha_RTI_dt, i_alpha_RTI, k, nvar, resid%d1val2)
         call ep1(s, i_dalpha_RTI_dt, i_alpha_RTI, k, nvar, resid%d1val3)

         ! Debugging
         if (test_partials) then
            write(*,2) 'a_00', k, a_00%val
            write(*,2) 'dadt_mix', k, dadt_mix
            write(*,2) 'dadt_source', k, dadt_source
            write(*,2) 'source_plus', k, source_plus
            write(*,2) 'source_minus', k, source_minus
            write(*,2) 'A_plus_B_div_rho', k, A_plus_B_div_rho
            write(*,2) 'instability', k, instability
            write(*,2) 's% q(k)', k, s% q(k)
            write(*,2) 's% Peos(k)', k, s% Peos(k)
            write(*,2) 's% rho(k)', k, s% rho(k)
            write(*,2) 's% ctrl% alpha_RTI_src_max_q', k, s% ctrl% alpha_RTI_src_max_q
            write(*,2) 'dadt_expected', k, dadt_expected
            write(*,2) 'dadt_actual', k, dadt_actual
            write(*,2) 'eqn_scale', k, eqn_scale
            write(*,2) 's% dt', k, s% dt
            write(*,2) 'equ(i,k)', k, s% equ(i_dalpha_RTI_dt,k)
            s% solver_test_partials_val = s% equ(i_dalpha_RTI_dt,k)
         end if


         if (test_partials) then   
            s% solver_test_partials_var = i_alpha_RTI
            s% solver_test_partials_dval_dx = resid%d1val2
            write(*,*) 'do1_dalpha_RTI_dt_eqn', s% solver_test_partials_var
         end if

      end subroutine do1_dalpha_RTI_dt_eqn


      end module hydro_alpha_rti_eqns

