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
      public :: do1_alpha


      contains


      subroutine do1_alpha( &
            s, xscale, k, nvar, skip_partials, equ, ierr)
         use star_utils, only: em1, e00, ep1
         use chem_def, only: ih1, ihe4

         type (star_info), pointer :: s
         real(dp), pointer :: xscale(:,:) ! (nvar, nz)
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         real(dp), intent(inout) :: equ(:,:)
         integer, intent(out) :: ierr

         integer, pointer :: reaction_id(:) ! maps net reaction number to reaction id
         integer :: nz, i, j, kk, i_lnd, i_lnT, i_lnR
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
            dequ_dlnT_const_Pgas, dlnT_dlnd_const_E, dVARdot_dVAR, fac

         include 'formats'

         ierr = 0

         dVARdot_dVAR = s% dVARdot_dVAR
         i = s% i_alpha_RTI
         nz = s% nz

         i_lnd = s% i_lnd
         i_lnT = s% i_lnT
         i_lnR = s% i_lnR
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
               s% q(k) < s% alpha_RTI_src_min_q) then
            source_plus = 0d0
            instability2 = 0d0
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
         dadt_actual = s% dalpha_RTI_dt(k)

         if (associated(xscale)) then
            eqn_scale = xscale(i,k)*dVARdot_dVAR
            equ(i,k) = (dadt_expected - dadt_actual)/eqn_scale
         else
            equ(i,k) = dadt_expected
            eqn_scale = 1d0
         end if

         if (skip_partials) return

         ! partial of -dadt_actual/eqn_scale
         call e00(s, xscale, i, i, k, nvar, -dVARdot_dVAR/eqn_scale)

         if (dadt_mix /= 0d0) then ! partials of dadt_mix/eqn_scale
            call e00(s, xscale, i, i, k, nvar, d_dadt_mix_da00/eqn_scale)
            if (k > 1) then
               call em1(s, xscale, i, i, k, nvar, d_dadt_mix_dam1/eqn_scale)
            end if
            if (k < nz) then
               call ep1(s, xscale, i, i, k, nvar, d_dadt_mix_dap1/eqn_scale)
            end if
         end if

         if (dadt_source == 0d0) return
         ! partials of dadt_source/eqn_scale

         call e00(s, xscale, i, i, k, nvar, ds_da/eqn_scale)

      end subroutine do1_alpha


      end module hydro_alpha_rti_eqns

