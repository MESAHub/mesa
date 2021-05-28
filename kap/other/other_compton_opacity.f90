! ***********************************************************************
!
!   Copyright (C) 2020-2021 The MESA Team
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
            integer, intent(in) :: handle ! kap handle
            real(dp), intent(in) :: Rho, T
            real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
            real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
            ! eta := electron degeneracy parameter from eos
            real(dp), intent(out) :: kap ! electron conduction opacity
            real(dp), intent(out) :: dlnkap_dlnRho, dlnkap_dlnT
            integer, intent(out) :: ierr ! 0 means AOK.

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
