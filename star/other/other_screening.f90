! ***********************************************************************
!
!   Copyright (C) 2021  The MESA Team
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

module other_screening

   implicit none

   ! set use_other_screening = .true. in your controls namelist

   ! edit the extras_controls routine to set the procedure pointer
   !  other_screening => my_screening

contains

   subroutine default_other_screening(sc, z1, z2, a1, a2, screen, dscreendt, dscreendd, ierr)
      use rates_def


      type(Screen_Info) :: sc  ! See rates_def
      ! This contains lots of useful things like temperature, density etc as well as some precomputed
      ! terms that are useful for screening calculations. The derived type is set in do_screen_set_context (screen.f90)
      real(dp), intent(in) ::    z1, z2      !< charge numbers of reactants
      real(dp), intent(in) ::    a1, a2     !< mass numbers of reactants
      real(dp), intent(out) ::   screen     !< on return, screening factor for this reaction
      real(dp), intent(out) ::   dscreendt     !< on return, temperature derivative of the screening factor
      real(dp), intent(out) ::   dscreendd    !< on return, density derivative of the screening factor
      integer, intent(out) ::   ierr

      screen = 0d0
      dscreendt = 0d0
      dscreendd = 0d0
      ierr = 0
   end subroutine default_other_screening

end module other_screening

