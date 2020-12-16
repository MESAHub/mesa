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
 
      module other_eval_i_rot

      ! consult star/other/README for general usage instructions
      ! control name: use_other_eval_i_rot = .true.
      ! procedure pointer: s% other_eval_i_rot => my_routine


      use star_def

      implicit none
      
            
      contains

      subroutine null_other_eval_i_rot(id,ri,r00,ra,w_div_w_crit_roche, i_rot, di_rot_dlnr, di_rot_dw_div_wc)
         integer, intent(in) :: id
         real(dp), intent(in) :: ri,r00,ra,w_div_w_crit_roche
         real(dp), intent(out) :: i_rot, di_rot_dlnr, di_rot_dw_div_wc

         write(*,*) 'no implementation for other_eval_i_rot'
         stop

      end subroutine null_other_eval_i_rot
      

      end module other_eval_i_rot
      
      
      
      
