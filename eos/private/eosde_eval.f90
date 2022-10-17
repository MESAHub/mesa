! ***********************************************************************
!
!   Copyright (C) 2014-2019  The MESA Team
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

module eosDE_eval
   use eos_def
   use const_def
   use math_lib
   
   implicit none


contains
   
   
   subroutine Get_eos_gamma_DE_Results(&
      rq, abar, energy, log10E, rho, log10Rho, gamma, &
      T, log10T, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
      dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
      dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, ierr)
      use utils_lib, only : is_bad
      type (EoS_General_Info), pointer :: rq
      real(dp), intent(in) :: abar, energy, log10E, Rho, log10Rho, gamma
      real(dp), intent(out) :: T, log10T
      real(dp), intent(inout), dimension(:) :: &
         res, d_dlnRho_const_T, d_dlnT_const_Rho
      real(dp), intent(out) :: &
         dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
         dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
      integer, intent(out) :: ierr
      
      real(dp) :: avo_k_div_abar, P, entropy
      include 'formats'
      
      ierr = 0
      
      res(1:nv) = 0d0
      d_dlnRho_const_T(1:nv) = 0d0
      d_dlnT_const_Rho(1:nv) = 0d0
      
      avo_k_div_abar = avo * kerg / abar
      
      P = (gamma - 1d0) * energy * rho
      T = (gamma - 1d0) * energy / avo_k_div_abar
      log10T = log10(T)
      
      res(i_Cv) = avo_k_div_abar / (gamma - 1) ! energy/T
      
      entropy = res(i_Cv) * log(P / pow(rho, gamma)) + 1d9 ! offset to keep it > 0
      
      if (is_bad(entropy) .or. entropy <= 0d0) then
         if (.false.) then
            !$OMP critical (eosde_eval_crit1)
            write(*, 1) 'entropy', entropy
            write(*, 1) 'energy', energy
            write(*, 1) 'T', T
            write(*, 1) 'P', P
            write(*, 1) 'rho', rho
            write(*, 1) 'pow(rho,gamma)', pow(rho, gamma)
            write(*, 1) 'res(i_Cv)', res(i_Cv)
            write(*, 1) 'P/pow(rho,gamma)', P / pow(rho, gamma)
            write(*, 1) 'log(P/pow(rho,gamma)', log(P / pow(rho, gamma))
            write(*, 1) 'gamma', gamma
            call mesa_error(__FILE__, __LINE__, 'bad gamma law entropy')
            !$OMP end critical (eosde_eval_crit1)
         end if
         entropy = 1d-99
      end if
      
      res(i_lnPgas) = log(P) ! treat P as Pgas
      res(i_lnE) = log10E * ln10
      res(i_lnS) = log(entropy)
      res(i_Cp) = gamma * res(i_Cv)
      res(i_grad_ad) = (gamma - 1d0) / gamma ! dlnT_dlnP|S
      res(i_chiRho) = 1d0 ! dlnP_dlnRho|T
      res(i_chiT) = 1d0 ! dlnP_dlnT|Rho
      res(i_dE_dRho) = 0d0
      res(i_dS_dT) = 0d0
      res(i_dS_dRho) = 0d0
      res(i_mu) = 1d0 / abar
      res(i_gamma1) = gamma
      res(i_gamma3) = gamma
      res(i_lnfree_e) = -1d99
      res(i_eta) = 0d0
      
      d_dlnRho_const_T(i_lnPgas) = 1
      
      d_dlnT_const_Rho(i_lnPgas) = 1
      d_dlnT_const_Rho(i_lnE) = 1
      
      dlnT_dlnd_c_E = 0
      dlnPgas_dlnd_c_E = 1
      
      dlnT_dlnE_c_Rho = 1
      dlnPgas_dlnE_c_Rho = 1
   
   end subroutine Get_eos_gamma_DE_Results


end module eosDE_eval
      
