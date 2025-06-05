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

      module mod_other_tsync

      ! NOTE: remember to set true:
      ! use_other_tsync = .true.

      implicit none

      contains

      subroutine null_other_tsync(id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
         use const_def, only: dp, strlen
         use binary_def, only : binary_info, binary_ptr
         use star_def, only : star_info, star_ptr
         integer, intent(in) :: id
         character (len=strlen), intent(in) :: sync_type  ! synchronization timescale
         real(dp), intent(in) :: Ftid  ! efficiency of tidal synchronization. (time scale / Ftid ).
         real(dp), intent(in) :: qratio  ! mass_other_star/mass_this_star
         real(dp), intent(in) :: m
         real(dp), intent(in) :: r_phot
         real(dp), intent(in) :: osep  ! orbital separation (cm)
         real(dp), intent(out) :: t_sync
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         t_sync = 1d99
      end subroutine null_other_tsync

      end module mod_other_tsync

