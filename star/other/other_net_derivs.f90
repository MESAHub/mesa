! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
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

module other_net_derivs

   implicit none

   ! set use_other_net_derivs = .true. in your controls namelist

   ! edit the extras_controls routine to set the procedure pointer
   !  other_net_derivs => my_net_derivs

   ! This hook is aimed at being able to modify dydt term.
   ! This means you can add arbitary new reactions that MESA does not know about
   ! by changing the dydt (change in compostion with time)

   ! This hook only works with soft nets (so no approx networks)

contains

   subroutine default_other_net_derivs( &
      n, dydt, eps_nuc_MeV, eta, ye, logtemp, temp, den, abar, zbar, &
      num_reactions, rate_factors, &
      symbolic, just_dydt, ierr)
      use net_def


      type(Net_Info) :: n
      real(qp), pointer, intent(inout) :: dydt(:, :)
      real(qp), intent(out) :: eps_nuc_MeV(:)
      integer, intent(in) :: num_reactions
      real(dp), intent(in) ::eta, ye, logtemp, temp, den, abar, zbar, &
                              rate_factors(:)
      logical, intent(in) :: symbolic, just_dydt
      integer, intent(out) :: ierr

   end subroutine default_other_net_derivs

end module other_net_derivs

