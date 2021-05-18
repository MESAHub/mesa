! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
 
      module other_mesh_functions

      ! consult star/other/README for general usage instructions
      ! control name: use_other_mesh_functions = .true.
      ! procedure pointer: s% other_mesh_fcn_data => my_routine


      ! remember to set use_other_mesh_functions = .true. to enable this.
      ! edit the extras_controls routine to set the procedure pointers
      ! e.g.,
         ! s% how_many_other_mesh_fcns => how_many_my_other_mesh_fcns
         ! s% other_mesh_fcn_data => my_other_mesh_fcn_data

      use star_def
      use math_lib

      implicit none
      
            
      contains
      
      
      subroutine null_how_many_other_mesh_fcns(id, n)
         integer, intent(in) :: id
         integer, intent(out) :: n
         n = 0
      end subroutine null_how_many_other_mesh_fcns
      
      
      subroutine null_other_mesh_fcn_data( &
            id, nfcns, names, gval_is_xa_function, vals1, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: nfcns
         character (len=*) :: names(:)
         logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
         real(dp), pointer :: vals1(:) ! =(nz, nfcns)
         integer, intent(out) :: ierr
         gval_is_xa_function(1:nfcns) = .false.
         ierr = 0
      end subroutine null_other_mesh_fcn_data

      
      ! here is an example that adds a mesh function for log(opacity)
      subroutine how_many_other_mesh_fcns(id, n)
         integer, intent(in) :: id
         integer, intent(out) :: n
         n = 1
      end subroutine how_many_other_mesh_fcns
      
      
      subroutine other_mesh_fcn_data( &
            id, nfcns, names, gval_is_xa_function, vals1, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: nfcns
         character (len=*) :: names(:)
         logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
         real(dp), pointer :: vals1(:) ! =(nz, nfcns)
         integer, intent(out) :: ierr
         integer :: nz, k
         real(dp), pointer :: vals(:,:)
         real(dp), parameter :: weight = 20
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = 'kap_function'
         gval_is_xa_function(1) = .false.
         nz = s% nz
         vals(1:nz,1:nfcns) => vals1(1:nz*nfcns)
         do k=1,nz
            vals(k,1) = weight*log10(s% opacity(k))
         end do
      end subroutine other_mesh_fcn_data


      end module other_mesh_functions
      
      
      
      
