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
 
      module other_diffusion_coefficients

      ! consult star/other/README for general usage instructions
      ! control name: use_other_diffusion_coefficients = .true.
      ! procedure pointer: s% other_diffusion_coefficients => my_routine


      use star_def

      implicit none
      
            
      contains



! Compute atomic diffusion coefficients.
! The input parameters are 
! rho          density [gcm^-3]
! T            temperature [K]
! m            number of elements  
! A            mass [amu], note element NN is electrons
! charge       charge [e]
! na           number density [cm^-3]  
! The output are the resistance coefficients in Burgers equations K_ij, z_ij, z'_ij and z''_ij
! It is also possible to output diffusion coefficients D_ij and thermal diffusion coeffcient
!     A_th used in Cowling&Chapman formalism, note Ath(m,i) is Ath_ei
      
      
      
! NOTE: the number of classes of isos = m-1; m is for electrons
! for j from 1 to m-1, you can get the chem id for class j by
!        cid = chem_get_iso_id(s% ctrl% diffusion_class_representative(j))
! e.g., if the representative for class j is he4, then cid will = ihe4 (defined in chem_def)

      
      subroutine null_other_diffusion_coefficients( &
            id, k, nc, m, rho, T, A, X, Z, C, charge, na, &
            Ddiff, Kdiff, Zdiff, Zdiff1, Zdiff2, Ath)
         use const_def, only: dp
         integer, intent(in) :: id, k, nc, m  
         real(dp), intent(in) :: rho, T, charge(m), na(m)
         real(dp), intent(in), dimension(:) :: A, X, Z, C ! (m)
         real(dp), intent(inout), dimension(m,m) :: &
            Ddiff, Kdiff, Zdiff, Zdiff1, Zdiff2, Ath
      end subroutine null_other_diffusion_coefficients


      end module other_diffusion_coefficients
      
      
      
      
