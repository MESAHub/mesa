! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton, Pablo Marchant & The MESA Team
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

      module mod_other_sync_spin_to_orbit

      ! NOTE: remember to set true:
      ! use_other_sync_spin_to_orbit = .true.

      implicit none

      contains

      subroutine null_other_sync_spin_to_orbit(id, nz, osep, qratio, rl, dt_next, Ftid, sync_type, sync_mode, ierr)
         use const_def, only: dp, strlen
         use binary_def, only : binary_info, binary_ptr
         use star_def, only : star_info, star_ptr
         integer, intent(in) :: id
         integer, intent(in) :: nz
         real(dp), intent(in) :: osep  ! orbital separation (cm)
         real(dp), intent(in) :: qratio  ! mass_other_star/mass_this_star
         real(dp), intent(in) :: rl  ! roche lobe radius (cm)
         real(dp), intent(in) :: dt_next  ! next timestep
         real(dp), intent(in) :: Ftid  ! efficiency of tidal synchronization. (time scale / Ftid ).

         character (len=strlen), intent(in) :: sync_type  ! synchronization timescale
         character (len=strlen), intent(in) :: sync_mode  ! where to put/take angular momentum
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: k

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         do k=1,nz
            s% extra_jdot(k) = s% extra_jdot(k) - 0d0
         end do

         write(*,*) "Warning: no implementation for other_sync_spin_to_orbit"
         write(*,*) "Not tidal torques are being included"
      end subroutine null_other_sync_spin_to_orbit

      end module mod_other_sync_spin_to_orbit

