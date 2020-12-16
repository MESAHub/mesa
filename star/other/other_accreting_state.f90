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
 
      module other_accreting_state

      ! consult star/other/README for general usage instructions
      ! control name: use_other_accreting_state = .true.
      ! procedure pointer: s% other_accreting_state => my_routine

      implicit none
      
            
      contains
      
      ! Note that your routine will be called before many star variables have been set.
      ! If you rely on these, you should call the star_set_vars_in_part1 routine from star_lib
      ! to ensure that they are set.      
      subroutine null_other_accreting_state(id, total_specific_energy, accretion_pressure, accretion_density, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(out) :: total_specific_energy, accretion_pressure, accretion_density
         integer, intent(out) :: ierr
         total_specific_energy = 0d0 ! erg/g
         accretion_pressure = 0d0 ! erg/cm^3
         accretion_density = 0d0 ! g/cm^3
         ierr = 0
      end subroutine null_other_accreting_state


      end module other_accreting_state
      
      
      
      
