! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
 
      module other_kap

      ! consult star/other/README for general usage instructions
      ! control name: use_other_kap = .true.
      ! procedure pointers: s% other_kap_get => my_routine
      ! (if using OP MONO)  s% other_kap_get_op_mono => my_routine


      use star_def

      implicit none
      
            
      contains


      subroutine null_other_kap_get( &
            id, k, handle, species, chem_id, net_iso, xa, &
            log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)

         use kap_def, only: num_kap_fracs

         ! INPUT
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! from alloc_kap_handle
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         real(dp), intent(in) :: log10_rho ! density
         real(dp), intent(in) :: log10_T ! temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
            ! eta := electron degeneracy parameter
         
         ! OUTPUT
         real(dp), intent(out) :: kap_fracs(num_kap_fracs)
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
         real(dp), intent(out) :: dln_kap_dxa(:) ! partial derivative w.r.t. to species
         integer, intent(out) :: ierr ! 0 means AOK.
                  
         kap_fracs = 0; kap = 0; dln_kap_dlnRho = 0; dln_kap_dlnT = 0; dln_kap_dxa = 0
         
         write(*,*) 'no implementation for other_kap_get'
         ierr = -1

      end subroutine null_other_kap_get
      
      
      subroutine null_other_kap_get_op_mono( &
            handle, zbar, log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            ! args for op_mono
            use_op_mono_alt_get_kap, &
            nel, izzp, fap, fac, screening, umesh, semesh, ff, rs, &
            ! output
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
         real(dp), intent(in) :: log10_rho ! the density
         real(dp), intent(in) :: log10_T ! the temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         ! args for op_mono_get_kap
         logical, intent(in) :: use_op_mono_alt_get_kap
         integer, intent(in) :: nel
         integer, intent(in) :: izzp(:) ! (nel)
         real(dp), intent(in) :: fap(:) ! (nel) number fractions of elements
         real(dp), intent(in) :: fac(:) ! (nel) scale factors for element opacity
         logical, intent(in) :: screening
         ! work arrays
         real, pointer :: umesh(:), semesh(:), ff(:,:,:,:), rs(:,:,:)
            ! umesh(nptot)
            ! umesh(nptot)
            ! ff(nptot, ipe, 4, 4)
            ! rs(nptot, 4, 4)
            ! ss(nptot, nrad, 4, 4)
         ! output
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
         kap = 0; dlnkap_dlnRho = 0; dlnkap_dlnT = 0
         ierr = -1
      end subroutine null_other_kap_get_op_mono
      


      end module other_kap
      
      
      
      
