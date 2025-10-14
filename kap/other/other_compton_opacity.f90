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

      module other_compton_opacity

      ! consult star/other/README for general usage instructions
      ! kap namelist option: use_other_compton_opacity = .true.
      ! procedure pointers: s% kap_rq % other_compton_opacity => my_routine

      use kap_def

      implicit none


      contains

         subroutine null_other_compton_opacity( &
            handle, &
            Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            use const_def, only: dp
            integer, intent(in) :: handle  ! kap handle; from star, pass s% kap_handle
            real(dp), intent(in) :: Rho, T
            real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
            real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
            ! eta := electron degeneracy parameter from eos
            real(dp), intent(out) :: kap  ! electron conduction opacity
            real(dp), intent(out) :: dlnkap_dlnRho, dlnkap_dlnT
            integer, intent(out) :: ierr  ! 0 means AOK.

            write(*,*) 'no implementation for other_compton_opacity'
            ierr = -1

            ! can first call kap_lib routine to get standard results, if desired

            ! subroutine kap_get_compton_opacity( &
            !    handle, &
            !    Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            !    eta, d_eta_dlnRho, d_eta_dlnT, &
            !    kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

         end subroutine null_other_compton_opacity

      end module other_compton_opacity
