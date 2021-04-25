! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
 
      module other_mlt_results

      ! consult star/other/README for general usage instructions
      ! control name: use_other_mlt_results = .true.
      ! procedure pointer: s% other_mlt_results => my_routine


      use star_def

      implicit none
      
            
      contains
      
      
      subroutine null_other_mlt_results(id, k, MLT_option, &  ! NOTE: k=0 is a valid arg
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            alpha_semiconvection, thermohaline_coeff, &
            mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
         use const_def, only: dp
         use auto_diff
         integer, intent(in) :: id
         integer, intent(in) :: k
         character (len=*), intent(in) :: MLT_option
         type(auto_diff_real_star_order1), intent(in) :: &
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height
         integer, intent(in) :: iso
         real(dp), intent(in) :: &
            XH1, cgrav, m, gradL_composition_term, &
            mixing_length_alpha, alpha_semiconvection, thermohaline_coeff
         integer, intent(out) :: mixing_type
         type(auto_diff_real_star_order1), intent(out) :: &
            gradT, Y_face, conv_vel, D, Gamma
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_other_mlt_results


      end module other_mlt_results
      
      
      
      
