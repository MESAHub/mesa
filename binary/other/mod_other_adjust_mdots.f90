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

      module mod_other_adjust_mdots

      ! NOTE: remember to set true:
      ! use_other_adjust_mdots = .true.

      implicit none

      contains

      subroutine null_other_adjust_mdots(binary_id, ierr)
         use binary_def, only : binary_info, binary_ptr
         use const_def, only: dp
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% fixed_xfer_fraction = 0d0
         ! note that you should also eval mdot_edd here and the wind mass transfer eff
         ! and angular momentum accretion. You can call the default functions using
         ! the ones provided through binary_lib
         b% mdot_edd = 0d0
         b% mdot_edd_eta = 0d0
         b% mdot_wind_transfer(:) = 0d0
         b% wind_xfer_fraction(:) = 0d0
         b% s_donor% mstar_dot = 0d0
         if (b% point_mass_i == 0) then
            b% component_mdot(b% a_i) = 0d0
         else
            b% s_accretor% mstar_dot = 0d0
         end if
         b% accretion_luminosity = 0d0
         b% mdot_system_transfer(:) = 0d0
         b% mdot_system_cct = 0d0

      end subroutine null_other_adjust_mdots

      end module mod_other_adjust_mdots

