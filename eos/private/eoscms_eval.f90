! ***********************************************************************
!
!   Copyright (C) 2020  The MESA Team
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
!
! ***********************************************************************

module eoscms_eval
   use eos_def
   use num_lib
   use const_def, only : avo, crad, ln10, arg_not_provided, mp, kerg, dp, qp, mesa_data_dir
   use utils_lib, only : is_bad, mesa_error
   use math_lib
   use interp_2d_lib_db
   
   implicit none
   
   logical, parameter :: CMS_cubic_in_X = .false.
   
   integer, parameter :: CMS_num_Xs = 11
   integer, parameter :: min_for_cubic = 2
   integer, parameter :: max_for_cubic = CMS_num_Xs - 2
   
   character(len = 3) :: CMS_Xstr(CMS_num_Xs) = ['000', '010', '020', '030', '040', '050', '060', '070', '080', '090', '100']
   
   real(dp) :: CMS_Xvals(CMS_num_Xs) = [ 0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp, 0.7_dp, 0.8_dp, 0.9_dp, 1.0_dp]
   
   !modeled on EosDT_XZ_Info from eos_def
   type eosCMS_X_Info
      real(dp) :: logRho_min
      real(dp) :: logRho_max
      real(dp) :: delta_logRho
      integer :: num_logRhos
      real(dp) :: logT_min
      real(dp) :: logT_max
      real(dp) :: delta_logT
      integer :: num_logTs
      real(dp), pointer :: logRhos(:), logTs(:)
      real(dp), pointer :: tbl1(:)
      integer :: version
   end type eosCMS_X_Info
   
   type(eosCMS_X_Info), target :: eosCMS_X_data(CMS_num_Xs)
   logical :: eosCMS_X_loaded(CMS_num_Xs) = .false.

contains
   
   subroutine eosCMS_init(ierr)
      integer, intent(out) :: ierr
      eosCMS_X_loaded = .false.
      ierr = 0
   end subroutine eosCMS_init
   
   
   subroutine Get_CMS_alfa(&
      rq, logRho, logT, Z, abar, zbar, &
      alfa, d_alfa_dlogT, d_alfa_dlogRho, &
      ierr)
      use const_def
      use auto_diff
      type (EoS_General_Info), pointer :: rq
      real(dp), intent(in) :: logRho, logT, Z, abar, zbar
      real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
      integer, intent(out) :: ierr
      
      type(auto_diff_real_2var_order1) :: logT_auto, logRho_auto, logQ_auto
      type(auto_diff_real_2var_order1) :: blend, blend_logT, blend_logRho, blend_logQ
      
      include 'formats'
      
      ierr = 0
      
      ! logRho is val1
      logRho_auto% val = logRho
      logRho_auto% d1val1 = 1d0
      logRho_auto% d1val2 = 0d0
      
      ! logT is val2
      logT_auto% val = logT
      logT_auto% d1val1 = 0d0
      logT_auto% d1val2 = 1d0
      
      logQ_auto = logRho_auto - 2d0 * logT_auto + 12d0
      
      ! for blend variables 1 is CMS, 0 is other
      ! (this is the opposite of the final alfa)
      
      ! logT blend
      if (logT_auto < rq% logT_min_for_any_CMS) then
         blend_logT = 0d0
      else if (logT_auto < rq% logT_min_for_all_CMS) then
         blend_logT = (logT_auto - rQ% logT_min_for_any_CMS) / (rq% logT_min_for_all_CMS - rq% logT_min_for_any_CMS)
      else if (logT_auto < rq% logT_max_for_all_CMS) then
         blend_logT = 1d0
      else if (logT_auto < rq% logT_max_for_any_CMS) then
         blend_logT = (logT_auto - rQ% logT_max_for_any_CMS) / (rq% logT_max_for_all_CMS - rq% logT_max_for_any_CMS)
      else
         blend_logT = 0
      end if
      
      
      ! logRho blend
      if (logRho_auto < rq% logRho_min_for_any_CMS) then
         blend_logRho = 0d0
      else if (logRho_auto < rq% logRho_min_for_all_CMS) then
         blend_logRho = (logRho_auto - rQ% logRho_min_for_any_CMS) / (rq% logRho_min_for_all_CMS - rq% logRho_min_for_any_CMS)
      else if (logRho_auto < rq% logRho_max_for_all_CMS) then
         blend_logRho = 1d0
      else if (logRho_auto < rq% logRho_max_for_any_CMS) then
         blend_logRho = (logRho_auto - rQ% logRho_max_for_any_CMS) / (rq% logRho_max_for_all_CMS - rq% logRho_max_for_any_CMS)
      else
         blend_logRho = 0
      end if
      
      
      ! logQ blend
      if (logQ_auto < rq% logQ_min_for_any_CMS) then
         blend_logQ = 0d0
      else if (logQ_auto < rq% logQ_min_for_all_CMS) then
         blend_logQ = (logQ_auto - rQ% logQ_min_for_any_CMS) / (rq% logQ_min_for_all_CMS - rq% logQ_min_for_any_CMS)
      else if (logQ_auto < rq% logQ_max_for_all_CMS) then
         blend_logQ = 1d0
      else if (logQ_auto < rq% logQ_max_for_any_CMS) then
         blend_logQ = (logQ_auto - rQ% logQ_max_for_any_CMS) / (rq% logQ_max_for_all_CMS - rq% logQ_max_for_any_CMS)
      else
         blend_logQ = 0
      end if
      
      ! combine blends
      blend = blend_logRho * blend_logT * blend_logQ
      
      alfa = 1d0 - blend% val
      d_alfa_dlogRho = -blend% d1val1
      d_alfa_dlogT = -blend% d1val2
   
   end subroutine Get_CMS_alfa
   
   
   subroutine get_CMS_for_eosdt(&
      handle, dbg, Z, X, abar, zbar, &
      species, chem_id, net_iso, xa, &
      rho, logRho, T, logT, remaining_fraction, &
      res, d_dlnd, d_dlnT, d_dxa, &
      skip, ierr)
      use chem_def, only : chem_isos
      use interp_1d_lib, only : interp_4pt_pm
      integer, intent(in) :: handle
      logical, intent(in) :: dbg
      real(dp), intent(in) :: Z, X, abar, zbar, remaining_fraction
      integer, intent(in) :: species
      integer, pointer :: chem_id(:), net_iso(:)
      real(dp), intent(in) :: xa(:)
      real(dp), intent(in) :: rho, logRho, T, logT
      real(dp), intent(inout), dimension(nv) :: &
         res, d_dlnd, d_dlnT
      real(dp), intent(inout), dimension(nv, species) :: d_dxa
      logical, intent(out) :: skip
      integer, intent(out) :: ierr
      integer :: iX, i
      type (EoS_General_Info), pointer :: rq
      real(dp) :: res1(nv), res2(nv), res3(nv), res4(nv)
      real(dp) :: dres1_dlnT(nv), dres2_dlnT(nv), dres3_dlnT(nv), dres4_dlnT(nv)
      real(dp) :: dres1_dlnRho(nv), dres2_dlnRho(nv), dres3_dlnRho(nv), dres4_dlnRho(nv)
      real(dp) :: d_dX(nv)
      real(dp) :: alfa, beta, dbeta_dX, dalfa_dX, xx(4), y(4), a(3), dw, dX
      rq => eos_handles(handle)
      
      if(rq% CMS_use_fixed_composition)then !do fixed composition (one table only)
         
         if(rq% CMS_fixed_composition_index < 0 .or. rq% CMS_fixed_composition_index > 10)then
            write(*, *) 'invalid value for CMS_fixed_composition_index.  See eos.defaults.'
            ierr = -1
            return
         endif
         
         iX = rq% CMS_fixed_composition_index + 1
         
         call eval_eosCMS_fixed_X(iX, logRho, logT, res, d_dlnT, d_dlnd, ierr)
         
         if(ierr/=0) then
            write(*, *) 'failed in get_CMS_for_eosdt'
            return
         endif
         
         ! composition derivatives; here composition is constant so no change
         d_dxa = 0
      
      else !do full composition
         !locate X values in the tables such that Xvals(iX) <= X < Xvals(iX+1)
         if (X <= CMS_Xvals(1)) then
            iX = 1
            if (X < 0) write(*, *) 'warning: X < 0 in eosCMS'
         else if (X >= CMS_Xvals(CMS_num_Xs - 1)) then
            iX = CMS_num_Xs - 1
            if (X > 1) write(*, *) 'warning: X > 1 in eosCMS'
         else
            do i = 2, CMS_num_Xs - 1
               if (X < CMS_Xvals(i))then
                  iX = i - 1; exit
               endif
            enddo
         endif
         
         !interpolation always bicubic in logRho and logT
         if(CMS_cubic_in_X .and. iX > 2 .and. iX < CMS_num_Xs - 2)then
            !do cubic interpolation in X
            call eval_eosCMS_fixed_X(iX - 1, logRho, logT, res1, dres1_dlnT, dres1_dlnRho, ierr)
            call eval_eosCMS_fixed_X(iX, logRho, logT, res2, dres2_dlnT, dres2_dlnRho, ierr)
            call eval_eosCMS_fixed_X(iX + 1, logRho, logT, res3, dres3_dlnT, dres3_dlnRho, ierr)
            call eval_eosCMS_fixed_X(iX + 2, logRho, logT, res4, dres4_dlnT, dres4_dlnRho, ierr)
            if(ierr/=0) then
               write(*, *) 'failed in get_CMS_for_eosdt'
               return
            endif
            XX(1:4) = CMS_Xvals(iX - 1:iX + 2)
            dX = X - CMS_Xvals(iX) ! assumes fixed dX spacing
            do i = 1, nv
               !result
               y(1:4) = [res1(i), res2(i), res3(i), res4(i)]
               call interp_4pt_pm(XX, y, a)
               res(i) = y(2) + dX * (a(1) + dX * (a(2) + dX * a(3)))
               !dres/dX
               d_dX(i) = a(1) + dX * (2 * a(2) + 3 * dX * a(3))
               !dres/dlnT
               y(1:4) = [dres1_dlnT(i), dres2_dlnT(i), dres3_dlnT(i), dres4_dlnT(i)]
               call interp_4pt_pm(XX, y, a)
               d_dlnT(i) = y(2) + dX * (a(1) + dX * (a(2) + dX * a(3)))
               !dres/dlnRho
               y(1:4) = [dres1_dlnRho(i), dres2_dlnRho(i), dres3_dlnRho(i), dres4_dlnRho(i)]
               call interp_4pt_pm(XX, y, a)
               d_dlnd(i) = y(2) + dX * (a(1) + dX * (a(2) + dX * a(3)))
            enddo
         else !linear interpolation in X
            call eval_eosCMS_fixed_X(iX, logRho, logT, res1, dres1_dlnT, dres1_dlnRho, ierr)
            call eval_eosCMS_fixed_X(iX + 1, logRho, logT, res2, dres2_dlnT, dres2_dlnRho, ierr)
            if(ierr/=0) then
               write(*, *) 'failed in get_CMS_for_eosdt'
               return
            endif
            beta = (X - CMS_Xvals(iX)) / (CMS_Xvals(iX + 1) - CMS_Xvals(iX))
            alfa = 1._dp - beta
            dbeta_dX = 1d0 / (CMS_Xvals(iX + 1) - CMS_Xvals(iX))
            dalfa_dX = -dbeta_dX
            res = alfa * res1 + beta * res2
            d_dX = dalfa_dX * res1 + dbeta_dX * res2
            d_dlnT = alfa * dres1_dlnT + beta * dres2_dlnT
            d_dlnd = alfa * dres1_dlnRho + beta * dres2_dlnRho
         endif
         
         ! composition derivatives
         do i = 1, species
            select case(chem_isos% Z(chem_id(i))) ! charge
            case (1) ! X
               d_dxa(:, i) = d_dX
            case (2) ! Y
               d_dxa(:, i) = 0
            case default ! Z
               d_dxa(:, i) = 0
            end select
         end do
      
      endif
      
      skip = .false.
      
      ! CMS tables do not include radiation.  Add it.
      if (rq% include_radiation) call include_radiation(Z, X, abar, zbar, &
         species, chem_id, net_iso, xa, &
         rho, logRho, T, logT, &
         res, d_dlnd, d_dlnT, d_dxa, &
         ierr)
      
      ! zero phase information
      res(i_phase:i_latent_ddlnRho) = 0d0
      d_dlnT(i_phase:i_latent_ddlnRho) = 0d0
      d_dlnd(i_phase:i_latent_ddlnRho) = 0d0
      
      ! zero all components
      res(i_frac:i_frac + num_eos_frac_results - 1) = 0.0
      d_dlnd(i_frac:i_frac + num_eos_frac_results - 1) = 0.0
      d_dlnT(i_frac:i_frac + num_eos_frac_results - 1) = 0.0
      
      ! mark this one
      res(i_frac_CMS) = 1.0
   
   end subroutine get_CMS_for_eosdt
   
   
   subroutine include_radiation(Z, X, abar, zbar, &
      species, chem_id, net_iso, xa, &
      rho_in, logRho, T_in, logT, &
      res, d_dlnd, d_dlnT, d_dxa, &
      ierr)
      use auto_diff
      real(dp), intent(in) :: Z, X, abar, zbar
      integer, intent(in) :: species
      integer, pointer :: chem_id(:), net_iso(:)
      real(dp), intent(in) :: xa(:)
      real(dp), intent(in) :: rho_in, logRho, T_in, logT
      real(dp), intent(inout), dimension(nv) :: &
         res, d_dlnd, d_dlnT
      real(dp), intent(inout), dimension(nv, species) :: d_dxa
      integer, intent(out) :: ierr
      
      type(auto_diff_real_2var_order1) :: T, Rho, P, Pgas, Prad, lnPgas, lnE, lnS, &
         grad_ad, chiRho, chiT, Cp, Cv, dE_dRho, dS_dT, dS_dRho, mu, lnfree_e, gamma1, gamma3, eta
      
      ierr = 0
      
      T% val = T_in
      T% d1val1 = 1d0
      T% d1val2 = 0d0
      
      Rho% val = Rho_in
      Rho% d1val1 = 0d0
      Rho% d1val2 = 1d0
      
      
      ! unpack results into auto_diff types
      lnPgas % val = res(i_lnPgas)
      lnE % val = res(i_lnE)
      lnS % val = res(i_lnS)
      grad_ad % val = res(i_grad_ad)
      chiRho % val = res(i_chiRho)
      chiT % val = res(i_chiT)
      Cp % val = res(i_Cp)
      Cv % val = res(i_Cv)
      dE_dRho % val = res(i_dE_dRho)
      dS_dT % val = res(i_dS_dT)
      dS_dRho % val = res(i_dS_dRho)
      mu % val = res(i_mu)
      lnfree_e % val = res(i_lnfree_e)
      gamma1 % val = res(i_gamma1)
      gamma3 % val = res(i_gamma3)
      eta % val = res(i_eta)
      
      lnPgas % d1val1 = d_dlnT(i_lnPgas) / T% val
      lnE % d1val1 = d_dlnT(i_lnE) / T% val
      lnS % d1val1 = d_dlnT(i_lnS) / T% val
      grad_ad % d1val1 = d_dlnT(i_grad_ad) / T% val
      chiRho % d1val1 = d_dlnT(i_chiRho) / T% val
      chiT % d1val1 = d_dlnT(i_chiT) / T% val
      Cp % d1val1 = d_dlnT(i_Cp) / T% val
      Cv % d1val1 = d_dlnT(i_Cv) / T% val
      dE_dRho % d1val1 = d_dlnT(i_dE_dRho) / T% val
      dS_dT % d1val1 = d_dlnT(i_dS_dT) / T% val
      dS_dRho % d1val1 = d_dlnT(i_dS_dRho) / T% val
      mu % d1val1 = d_dlnT(i_mu) / T% val
      lnfree_e % d1val1 = d_dlnT(i_lnfree_e) / T% val
      gamma1 % d1val1 = d_dlnT(i_gamma1) / T% val
      gamma3 % d1val1 = d_dlnT(i_gamma3) / T% val
      eta % d1val1 = d_dlnT(i_eta) / T% val
      
      lnPgas % d1val2 = d_dlnd(i_lnPgas) / Rho% val
      lnE % d1val2 = d_dlnd(i_lnE) / Rho% val
      lnS % d1val2 = d_dlnd(i_lnS) / Rho% val
      grad_ad % d1val2 = d_dlnd(i_grad_ad) / Rho% val
      chiRho % d1val2 = d_dlnd(i_chiRho) / Rho% val
      chiT % d1val2 = d_dlnd(i_chiT) / Rho% val
      Cp % d1val2 = d_dlnd(i_Cp) / Rho% val
      Cv % d1val2 = d_dlnd(i_Cv) / Rho% val
      dE_dRho % d1val2 = d_dlnd(i_dE_dRho) / Rho% val
      dS_dT % d1val2 = d_dlnd(i_dS_dT) / Rho% val
      dS_dRho % d1val2 = d_dlnd(i_dS_dRho) / Rho% val
      mu % d1val2 = d_dlnd(i_mu) / Rho% val
      lnfree_e % d1val2 = d_dlnd(i_lnfree_e) / Rho% val
      gamma1 % d1val2 = d_dlnd(i_gamma1) / Rho% val
      gamma3 % d1val2 = d_dlnd(i_gamma3) / Rho% val
      eta % d1val2 = d_dlnd(i_eta) / Rho% val
      
      
      ! add radiation
      Pgas = exp(lnPgas)
      Prad = one_third * crad * pow4(T)
      P = Pgas + Prad
      lnE = log(exp(lnE) + 3d0 * Prad / Rho)
      lnS = log(exp(lnS) + 4d0 * Prad / (Rho * T))
      Cv = Cv + 12d0 * Prad / (Rho * T)
      chiT = chiT * Pgas / P + 4d0 * Prad / P
      chiRho = chiRho * Pgas / P
      gamma3 = 1d0 + P / Rho * chiT / (T * Cv)
      gamma1 = chiT * (gamma3 - 1d0) + chiRho
      grad_ad = (gamma3 - 1d0) / gamma1
      Cp = Cv * gamma1 / chiRho
      dE_dRho = (1d0 - chiT) * P / (Rho * Rho)
      dS_dT = Cv / T
      dS_dRho = -P * chiT / (Rho * Rho * T)
      
      
      ! repack results
      res(i_lnPgas) = lnPgas % val
      res(i_lnE) = lnE % val
      res(i_lnS) = lnS % val
      res(i_grad_ad) = grad_ad % val
      res(i_chiRho) = chiRho % val
      res(i_chiT) = chiT % val
      res(i_Cp) = Cp % val
      res(i_Cv) = Cv % val
      res(i_dE_dRho) = dE_dRho % val
      res(i_dS_dT) = dS_dT % val
      res(i_dS_dRho) = dS_dRho % val
      res(i_mu) = mu % val
      res(i_lnfree_e) = lnfree_e % val
      res(i_gamma1) = gamma1 % val
      res(i_gamma3) = gamma3 % val
      res(i_eta) = eta % val
      
      d_dlnT(i_lnPgas) = lnPgas % d1val1 * T % val
      d_dlnT(i_lnE) = lnE % d1val1 * T % val
      d_dlnT(i_lnS) = lnS % d1val1 * T % val
      d_dlnT(i_grad_ad) = grad_ad % d1val1 * T % val
      d_dlnT(i_chiRho) = chiRho % d1val1 * T % val
      d_dlnT(i_chiT) = chiT % d1val1 * T % val
      d_dlnT(i_Cp) = Cp % d1val1 * T % val
      d_dlnT(i_Cv) = Cv % d1val1 * T % val
      d_dlnT(i_dE_dRho) = dE_dRho % d1val1 * T % val
      d_dlnT(i_dS_dT) = dS_dT % d1val1 * T % val
      d_dlnT(i_dS_dRho) = dS_dRho % d1val1 * T % val
      d_dlnT(i_mu) = mu % d1val1 * T % val
      d_dlnT(i_lnfree_e) = lnfree_e % d1val1 * T % val
      d_dlnT(i_gamma1) = gamma1 % d1val1 * T % val
      d_dlnT(i_gamma3) = gamma3 % d1val1 * T % val
      d_dlnT(i_eta) = eta % d1val1 * T % val
      
      d_dlnd(i_lnPgas) = lnPgas % d1val2 * RHO % val
      d_dlnd(i_lnE) = lnE % d1val2 * RHO % val
      d_dlnd(i_lnS) = lnS % d1val2 * RHO % val
      d_dlnd(i_grad_ad) = grad_ad % d1val2 * RHO % val
      d_dlnd(i_chiRho) = chiRho % d1val2 * RHO % val
      d_dlnd(i_chiT) = chiT % d1val2 * RHO % val
      d_dlnd(i_Cp) = Cp % d1val2 * RHO % val
      d_dlnd(i_Cv) = Cv % d1val2 * RHO % val
      d_dlnd(i_dE_dRho) = dE_dRho % d1val2 * RHO % val
      d_dlnd(i_dS_dT) = dS_dT % d1val2 * RHO % val
      d_dlnd(i_dS_dRho) = dS_dRho % d1val2 * RHO % val
      d_dlnd(i_mu) = mu % d1val2 * RHO % val
      d_dlnd(i_lnfree_e) = lnfree_e % d1val2 * RHO % val
      d_dlnd(i_gamma1) = gamma1 % d1val2 * RHO % val
      d_dlnd(i_gamma3) = gamma3 % d1val2 * RHO % val
      d_dlnd(i_eta) = eta % d1val2 * RHO % val
   
   end subroutine include_radiation
   
   
   subroutine eval_eosCMS_fixed_X(iX, logRho, logT, res, dres_dlnT, dres_dlnRho, ierr)
      use eosdt_support, only : Do_EoS_Interpolations
      integer, intent(in) :: iX
      real(dp), intent(in) :: logRho, logT
      real(dp), intent(out) :: res(nv), dres_dlnT(nv), dres_dlnRho(nv)
      integer, intent(out) :: ierr
      type(eosCMS_X_info), pointer :: c
      integer :: iRho, iT
      real(dp) :: fval(nv), df_dx(nv), df_dy(nv)
      real(dp) :: logT0, logRho0, logT1, logRho1, my_logT, my_logRho
      ierr = 0
      !$OMP CRITICAL(OMP_CRITICAL_IX)
      if(.not.eosCMS_X_loaded(iX)) call load_eosCMS_table(iX, ierr)
      !$OMP END CRITICAL(OMP_CRITICAL_IX)
      
      my_logT = logT
      my_logRho = logRho
      
      c => eosCMS_X_data(iX)
      
      call locate_logRho(c, my_logRho, iRho, logRho0, logRho1)
      call locate_logT  (c, my_logT, iT, logT0, logT1)
      
      call Do_EoS_Interpolations(1, nv, nv, &
         c% num_logTs, c% logTs, c% num_logRhos, c% logRhos, &
         c% tbl1, iT, iRho, logT0, my_logT, logT1, &
         logRho0, my_logRho, logRho1, fval, df_dx, df_dy, ierr)
      
      if(ierr/=0) then
         write(*, *) 'failed in eval_eosCMS_fixed_X'
         return
      endif
      
      res = fval
      dres_dlnT = df_dx * iln10
      dres_dlnRho = df_dy * iln10
      
      ! specific heats were log
      res(i_Cp) = exp(fval(i_Cp))
      dres_dlnT(i_Cp) = dres_dlnT(i_Cp) * res(i_Cp)
      dres_dlnRho(i_Cp) = dres_dlnRho(i_Cp) * res(i_Cp)
      
      res(i_Cv) = exp(fval(i_Cv))
      dres_dlnT(i_Cv) = dres_dlnT(i_Cv) * res(i_Cv)
      dres_dlnRho(i_Cv) = dres_dlnRho(i_Cv) * res(i_Cv)
   
   end subroutine eval_eosCMS_fixed_X
   
   subroutine locate_logT(c, logT, iT, logT0, logT1)
      type(eosCMS_X_info), pointer :: c
      real(dp), intent(inout) :: logT
      integer, intent(out) :: iT
      real(dp), intent(out) :: logT0, logT1
      iT = int((logT - c% logT_min) / c% delta_logT + 0.0001_dp) + 1
      if(iT < 1 .or. iT >= c% num_logTs)then
         if(iT < 1) then
            iT = 1
            logT0 = c% logT_min
            logT1 = logT0 + c% delta_logT
            logT = logT0
         else
            iT = c% num_logTs - 1
            logT0 = c% logT_min + real(iT - 1, kind = dp) * c% delta_logT
            logT1 = logT0 + c% delta_logT
            logT = logT1
         endif
      else
         logT0 = c% logT_min + (iT - 1) * c% delta_logT
         logT1 = logT0 + c% delta_logT
      endif
   end subroutine locate_logT
   
   
   subroutine locate_logRho(c, logRho, iRho, logRho0, logRho1)
      type(eosCMS_X_info), pointer :: c
      real(dp), intent(inout) :: logRho
      integer, intent(out) :: iRho
      real(dp), intent(out) :: logRho0, logRho1
      iRho = int((logRho - c% logRho_min) / c% delta_logRho + 0.0001_dp) + 1
      if(iRho < 1 .or. iRho >= c% num_logRhos)then
         if(iRho < 1) then
            iRho = 1
            logRho0 = c% logRho_min
            logRho1 = logRho0 + c% delta_logRho
            logRho = logRho0
         else
            iRho = c% num_logRhos - 1
            logRho0 = c% logRho_min + real(iRho - 1, kind = dp) * c% delta_logRho
            logRho1 = logRho0 + c% delta_logRho
            logRho = logRho1
         endif
      else
         logRho0 = c% logRho_min + (iRho - 1) * c% delta_logRho
         logRho1 = logRho0 + c% delta_logRho
      endif
   end subroutine locate_logRho
   
   subroutine load_eosCMS_table(iX, ierr)
      integer, intent(in) :: iX
      integer, intent(out) :: ierr
      character(len = 256) :: filename, data_sub_dir
      character(len = 1024) :: message
      integer :: io, ios, i, j, n, v, ili_logTs, ili_logRhos, ibcxmin, ibcxmax, ibcymin, ibcymax
      real(dp), allocatable, target :: f1_ary(:)
      real(dp), pointer :: f1(:), f(:, :, :), vec(:), tbl(:, :, :, :)
      real(dp), target :: vec_ary(50)
      type(eosCMS_X_Info), pointer :: c
      real(dp), allocatable :: bcxmin(:), bcxmax(:), bcymin(:), bcymax(:)
      real(dp) :: X_in, Z_in
      
      ierr = 0
      c => eosCMS_X_data(iX)
      vec => vec_ary
      
      data_sub_dir = '/eosCMS_data/'
      filename = trim(mesa_data_dir) // trim(data_sub_dir) // 'mesa-CMS_' // CMS_Xstr(iX) // 'x.data'
      
      open(newunit = io, file = trim(filename), action = 'read', status = 'old', iostat = ios)
      if(ios/=0) then
         write(*, '(a)') 'failed while reading ' // trim(filename)
         call mesa_error(__FILE__, __LINE__)
         close(io)
         ierr = -1
         return
      endif
      
      read(io, *) !header
      read(io, '(a)') message
      call str_to_vector(message, vec, n, ierr)
      if(ierr/=0) return
      
      c% version = int(vec(1))
      X_in = vec(2)
      Z_in = vec(3)
      c% num_logTs = int(vec(4))
      c% logT_min = vec(5)
      c% logT_Max = vec(6)
      c% delta_logT = vec(7)
      c% num_logRhos = int(vec(8))
      c% logRho_min = vec(9)
      c% logRho_max = vec(10)
      c% delta_logRho = vec(11)
      
      allocate(c% logTs(c% num_logTs))
      allocate(c% logRhos(c% num_logRhos))
      
      c% logTs(1) = c% logT_min
      do i = 2, c% num_logTs
         c% logTs(i) = c% logTs(i - 1) + c% delta_logT
      enddo
      
      c% logRhos(1) = c% logRho_min
      do i = 2, c% num_logRhos
         c% logRhos(i) = c% logRhos(i - 1) + c% delta_logRho
      enddo
      
      !check that the input X value is compatible with the expected value
      if(abs(X_in - CMS_Xvals(iX)) > 0.01_dp) then
         write(*, *) ' eosCMS X value does not match table '
         write(*, *) ' expected X = ', CMS_Xvals(iX)
         write(*, *) ' received X = ', X_in
         call mesa_error(__FILE__, __LINE__)
         ierr = -1
         return
      endif
      
      read(io, *) !header
      read(io, *) !header
      
      allocate(c% tbl1(sz_per_eos_point * nv * c% num_logRhos * c% num_logTs))
      tbl(1:sz_per_eos_point, 1:nv, 1:c% num_logTs, 1:c% num_logRhos) => &
         c% tbl1(1:sz_per_eos_point * nv * c% num_logTs * c% num_logRhos)
      
      allocate(f1_ary(sz_per_eos_point * c% num_logRhos * c% num_logTs))
      f1 => f1_ary
      f(1:sz_per_eos_point, 1:c% num_logTs, 1:c% num_logRhos) => &
         f1(1:sz_per_eos_point * c% num_logTs * c% num_logRhos)
      
      do i = 1, c% num_logTs
         do j = 1, c% num_logRhos
            read(io, '(a)') message
            call str_to_vector(message, vec, n, ierr)
            if(ierr/=0)then
               close(io)
               return
            endif
            tbl(1, i_lnPgas, i, j) = vec(3) * ln10
            tbl(1, i_lnE, i, j) = vec(4) * ln10
            tbl(1, i_lnS, i, j) = vec(5) * ln10
            tbl(1, i_chiRho, i, j) = vec(6)
            tbl(1, i_chiT, i, j) = vec(7)
            tbl(1, i_Cp, i, j) = safe_log(vec(8))
            tbl(1, i_Cv, i, j) = safe_log(vec(9))
            tbl(1, i_dE_dRho, i, j) = vec(10)
            tbl(1, i_dS_dT, i, j) = vec(11)
            tbl(1, i_dS_dRho, i, j) = vec(12)
            tbl(1, i_mu, i, j) = vec(13)
            tbl(1, i_lnfree_e, i, j) = vec(14)
            tbl(1, i_gamma1, i, j) = vec(15)
            tbl(1, i_gamma3, i, j) = vec(16)
            tbl(1, i_grad_ad, i, j) = vec(17)
            tbl(1, i_eta, i, j) = vec(18)
         enddo
      enddo
      
      close(io)
      
      ! logT is "x"
      ! logRho is "y"
      
      ! "not a knot" bc's at edges of tables
      allocate(bcxmin(c% num_logRhos), bcxmax(c% num_logRhos))
      allocate(bcymin(c% num_logTs), bcymax(c% num_logTs))
      ibcxmin = 0; bcxmin(:) = 0
      ibcxmax = 0; bcxmax(:) = 0
      ibcymin = 0; bcymin(:) = 0
      ibcymax = 0; bcymax(:) = 0
      
      !create table for bicubic spline
      do v = 1, nv
         do j = 1, c% num_logRhos
            do i = 1, c% num_logTs
               f(1, i, j) = tbl(1, v, i, j)
            enddo
         enddo
         
         call interp_mkbicub_db(&
            c% logTs, c% num_logTs, c% logRhos, c% num_logRhos, f1, c% num_logTs, &
            ibcxmin, bcxmin, ibcxmax, bcxmax, ibcymin, bcymin, ibcymax, bcymax, &
            ili_logTs, ili_logRhos, ierr)
         
         do j = 1, c% num_logRhos
            do i = 1, c% num_logTs
               tbl(2, v, i, j) = f(2, i, j)
               tbl(3, v, i, j) = f(3, i, j)
               tbl(4, v, i, j) = f(4, i, j)
            enddo
         enddo
      enddo
      
      if(ierr==0) eosCMS_X_loaded(iX) = .true.
   
   end subroutine load_eosCMS_table

end module eoscms_eval

