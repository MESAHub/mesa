! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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

module other_diffusion_coefficients

   ! consult star/other/README for general usage instructions
   ! control name: use_other_diffusion_coefficients = .true.
   ! procedure pointer: s% other_diffusion_coefficients => my_routine

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
!        cid = chem_get_iso_id(s% diffusion_class_representative(j))
! e.g., if the representative for class j is he4, then cid will = ihe4 (defined in chem_def)

   subroutine null_other_diffusion_coefficients( &
      id, k, nc, m, rho, T, A, X, Z, C, charge, na, &
      Ddiff, Kdiff, Zdiff, Zdiff1, Zdiff2, Ath)
      use star_def
      use const_def, only: dp
      integer, intent(in) :: id, k, nc, m
      real(dp), intent(in) :: rho, T, charge(m), na(m)
      real(dp), intent(in), dimension(:) :: A, X, Z, C  ! (m)
      real(dp), intent(inout), dimension(m, m) :: &
         Ddiff, Kdiff, Zdiff, Zdiff1, Zdiff2, Ath
   end subroutine null_other_diffusion_coefficients

end module other_diffusion_coefficients

