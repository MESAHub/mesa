! ***********************************************************************
!
!   Copyright (C) 2020-2021 The MESA Team
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

      module other_elect_cond_opacity

      ! consult star/other/README for general usage instructions
      ! kap namelist option: use_other_elect_cond_opacity = .true.
      ! procedure pointers: s% kap_rq % other_elect_cond_opacity => my_routine

      use kap_def

      implicit none


      contains

         subroutine null_other_elect_cond_opacity( &
            handle, &
            zbar, logRho, logT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            use const_def, only: dp
            integer, intent(in) :: handle  ! kap handle; from star, pass s% kap_handle
            real(dp), intent(in) :: zbar  ! average ionic charge (for electron conduction)
            real(dp), intent(in) :: logRho  ! the density
            real(dp), intent(in) :: logT  ! the temperature
            real(dp), intent(out) :: kap  ! electron conduction opacity
            real(dp), intent(out) :: dlnkap_dlnRho  ! partial derivative at constant T
            real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
            integer, intent(out) :: ierr  ! 0 means AOK.

            write(*,*) 'no implementation for other_elect_cond_opacity'
            ierr = -1

            ! can first call kap_lib routine to get standard results, if desired

            ! call kap_get_elect_cond_opacity( &
            !    handle, &
            !    zbar, logRho, logT, &
            !    kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

         end subroutine null_other_elect_cond_opacity

      end module other_elect_cond_opacity
