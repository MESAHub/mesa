! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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

module other_neu

   ! consult star/other/README for general usage instructions
   ! control name: use_other_neu = .true.
   ! procedure pointer: s% other_neu => my_routine

   implicit none

contains

   subroutine null_other_neu( &
      id, k, T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
      loss, sources, ierr)
      use star_def
      use neu_lib, only: neu_get
      use neu_def
      integer, intent(in) :: id  ! id for star
      integer, intent(in) :: k  ! cell number or 0 if not for a particular cell
      real(dp), intent(in) :: T  ! temperature
      real(dp), intent(in) :: log10_T  ! log10 of temperature
      real(dp), intent(in) :: Rho  ! density
      real(dp), intent(in) :: log10_Rho  ! log10 of density
      real(dp), intent(in) :: abar  ! mean atomic weight
      real(dp), intent(in) :: zbar  ! mean charge
      real(dp), intent(in) :: log10_Tlim
      logical, intent(inout) :: flags(num_neu_types)  ! true if should include the type of loss
      real(dp), intent(inout) :: loss(num_neu_rvs)  ! total from all sources
      real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
      integer, intent(out) :: ierr
      call neu_get( &
         T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
         loss, sources, ierr)
   end subroutine null_other_neu

end module other_neu

