! ***********************************************************************
!
!   Copyright (C) 2014  The MESA Team
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

module other_surface_PT

   ! consult star/other/README for general usage instructions
   ! control name: use_other_surface_PT = .true.
   ! procedure pointer: s% other_surface_PT => my_routine

   implicit none

contains

   ! star_utils:set_phot_info sets s% Teff before this is called
   ! see hydro_vars:set_Teff_info_for_eqns

   subroutine null_other_surface_PT(id, &
                                    skip_partials, &
                                    lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
                                    lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
      use const_def, only: dp
      use star_def
      !use star_lib, only: star_get_surf_PT
      integer, intent(in) :: id
      logical, intent(in) :: skip_partials
      real(dp), intent(out) :: &
         lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
         lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
      integer, intent(out) :: ierr
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

