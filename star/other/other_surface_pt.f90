! ***********************************************************************
!
!   Copyright (C) 2014  Bill Paxton
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
 
      module other_surface_PT

      ! consult star/other/README for general usage instructions
      ! control name: use_other_surface_PT = .true.
      ! procedure pointer: s% other_surface_PT => my_routine


      use star_def

      implicit none
      
      contains
      
      
      subroutine null_other_surface_PT(id, &
            skip_partials, &
            Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
         use const_def, only: dp
         !use star_lib, only: star_get_surf_PT
         integer, intent(in) :: id
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr
         Teff = 0
         lnT_surf = 0
         dlnT_dL = 0
         dlnT_dlnR = 0
         dlnT_dlnM = 0
         dlnT_dlnkap = 0
         lnP_surf = 0
         dlnP_dL = 0
         dlnP_dlnR = 0
         dlnP_dlnM = 0
         dlnP_dlnkap = 0
         ierr = -1
         
         !call star_get_surf_PT(id, &
         !   skip_partials, &
         !   Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
         !   lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
         
      end subroutine null_other_surface_PT


      end module other_surface_PT
      



