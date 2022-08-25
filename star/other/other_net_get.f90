! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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
 
      module other_net_get

      ! consult star/other/README for general usage instructions
      ! control name: use_other_net_get = .true.
      ! procedure pointer: s% other_net_get => my_routine


      use star_def

      implicit none
      
            
      contains
      
      
      subroutine null_other_net_get(  &
            id, k, net_handle, just_dxdt, n, num_isos, num_reactions,  &
            x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode, &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, ierr)      
         use net_lib, only: net_get
         use net_def, only: Net_Info
         integer, intent(in) :: id ! id for star         
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: net_handle
         logical, intent(in) :: just_dxdt
         type (Net_Info) :: n
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), intent(in)  :: x(:) ! (num_isos)
         real(dp), intent(in)  :: temp, log10temp ! log10 of temp
         real(dp), intent(in)  :: rho, log10rho ! log10 of rho
         real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
         real(dp), intent(in)  :: zbar  ! mean charge per nucleus
         real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
         real(dp), intent(in)  :: ye    
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! electron degeneracy from eos.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         real(dp), intent(out) :: eps_nuc ! ergs/g/s from burning after including losses to reaction neutrinos
         real(dp), intent(out) :: d_eps_nuc_dT
         real(dp), intent(out) :: d_eps_nuc_dRho
         real(dp), intent(inout) :: d_eps_nuc_dx(:) ! (num_isos)       
         real(dp), intent(inout) :: dxdt(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dRho(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dT(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dx(:,:) ! (num_isos, num_isos)            
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
         real(dp), intent(out) :: eps_neu_total ! ergs/g/s neutrinos from weak reactions
         integer, intent(in) :: screening_mode
         integer, intent(in) :: lwork ! size of work >= result from calling net_work_size
         real(dp), pointer :: work(:) ! (lwork)
         integer, intent(out) :: ierr ! 0 means okay
         call net_get( &
            net_handle, just_dxdt, n, num_isos, num_reactions,  &
            x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode, &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, ierr)
      end subroutine null_other_net_get


      end module other_net_get
      
      
      
      
