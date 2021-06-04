! ***********************************************************************
!
!   Copyright (C) 2021 The MESA Team
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

      module other_eos

      ! consult star/other/README for general usage instructions
      ! eos namelist option: use_other = .true.
      ! procedure pointers:
      !   s% eos_rq % other_eos_frac => my_other_eos_frac
      !   s% eos_rq % other_eos_component => my_other_eos_component
      !   s% eos_rq % other_eos_results => my_other_eos_results

      ! the subroutine other_eos_component allows you to provide an
      ! additional component eos.  you must provide a complete set of
      ! eos results.  this eos component has the highest priority in
      ! the blend.

      ! the subroutine other_eos_frac defines the region over which
      ! you want to use other_eos_component.  you can use this to
      ! replace all or part of the MESA eos and control the location
      ! the eos blends associated with this component.

      ! the subroutine other_eos_results lets you modify the final
      ! results from the eos call right before they are returned.
      ! this allows you to make small modifications to the existing
      ! eos results without having to provide a full replacement eos.

      use const_def

      implicit none

      contains

      subroutine null_other_eos_frac( &
              handle, &
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, &
              frac, dfrac_dlogRho, dfrac_dlogT, ierr)

         ! INPUT
         use chem_def, only: num_chem_isos

         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

         integer, intent(in) :: species ! number of species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions

         real(dp), intent(in) :: Rho, log10Rho ! the density
         real(dp), intent(in) :: T, log10T ! the temperature

         ! OUTPUT
         ! this routine must provide a fraction (in [0,1]) of the 'other' eos to use
         ! the remaining fraction (1-frac) will be provided by the standard MESA eos
         real(dp), intent(out) :: frac ! fraction of other_eos to use
         real(dp), intent(out) :: dfrac_dlogRho ! its partial derivative at constant T
         real(dp), intent(out) :: dfrac_dlogT   ! its partial derivative at constant Rho

         integer, intent(out) :: ierr ! 0 means AOK.

         ! default implementation uses other_eos_component everywhere
         frac = 1d0
         dfrac_dlogRho = 0d0
         dfrac_dlogT = 0d0

      end subroutine null_other_eos_frac


      subroutine null_other_eos_component( &
              handle, &
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, &
              res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_const_TRho, ierr)

         ! INPUT
         use chem_def, only: num_chem_isos

         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

         integer, intent(in) :: species ! number of species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions

         real(dp), intent(in) :: Rho, log10Rho ! the density
         real(dp), intent(in) :: T, log10T ! the temperature

         ! OUTPUT

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         ! partial derivatives of the basic results wrt lnd and lnT
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         ! d_dlnRho(i) = d(res(i))/dlnd|T
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         ! d_dlnT(i) = d(res(i))/dlnT|Rho
         real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_basic_results, species)
         ! d_dxa(i,j) = d(res(i))/dxa(j)|T,Rho

         integer, intent(out) :: ierr ! 0 means AOK.

         res = 0
         d_dlnRho_const_T = 0
         d_dlnT_const_Rho = 0
         d_dxa_const_TRho = 0

         write(*,*) 'no implementation for other_eos'
         ierr = -1

      end subroutine null_other_eos_component


      subroutine null_other_eos_results( &
              handle, &
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, &
              res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_const_TRho, ierr)

         ! INPUT
         use chem_def, only: num_chem_isos

         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

         integer, intent(in) :: species ! number of species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions

         real(dp), intent(in) :: Rho, log10Rho ! the density
         real(dp), intent(in) :: T, log10T ! the temperature

         ! OUTPUT

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         ! partial derivatives of the basic results wrt lnd and lnT
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         ! d_dlnRho(i) = d(res(i))/dlnd|T
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         ! d_dlnT(i) = d(res(i))/dlnT|Rho
         real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_basic_results, species)
         ! d_dxa(i,j) = d(res(i))/dxa(j)|T,Rho

         integer, intent(out) :: ierr ! 0 means AOK.

         ! default implementation does not modify results

      end subroutine null_other_eos_results


      end module other_eos
