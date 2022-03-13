! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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
 
      module other_eval_fp_ft

      ! consult star/other/README for general usage instructions
      ! control name: use_other_eval_fp_ft = .true.
      ! procedure pointer: s% other_eval_fp_ft => my_routine


      use star_def

      implicit none
      
            
      contains

      subroutine null_other_eval_fp_ft( &
            id, nz, xm, r, rho, aw, ft, fp, r_polar, r_equatorial, report_ierr, ierr)
         use num_lib
         use star_utils
         use auto_diff_support
         integer, intent(in) :: id
         integer, intent(in) :: nz
         real(dp), intent(in) :: aw(:), r(:), rho(:), xm(:) ! (nz)
         type(auto_diff_real_star_order1), intent(out) :: ft(:), fp(:) ! (nz)
         real(dp), intent(inout) :: r_polar(:), r_equatorial(:) ! (nz)
         logical, intent(in) :: report_ierr
         integer, intent(out) :: ierr

         write(*,*) 'no implementation for other_eval_fp_ft'
         ! must set fp, ft, r_polar and r_equatorial
         ierr = -1

      end subroutine null_other_eval_fp_ft
      

      end module other_eval_fp_ft
      
      
      
      
