! ***********************************************************************
!
!   Copyright (C) 2010-2019  Meridith Joyce & The MESA Team
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

! Module for handling starspot calculations
! contains parameterized YREC routines
! note: MESA models a pressure contrast rather than temperature contrast

module starspots

   use star_private_def
   use const_def, only: dp, pi4, boltz_sigma, cgas
   use utils_lib

   implicit none

   private
   public :: starspot_tweak_gradr
   public :: starspot_tweak_PT
   public :: starspot_restore_PT

   real(dp) :: L_init

contains

   subroutine starspot_tweak_gradr(s, P, gradr, gradr_spot)
      ! Purpose: adjusts the gradient of the radius to account for starspots.
      !
      ! This subroutine is called at the beginning of Get_results()
      ! in turb_support.f90

      use auto_diff_support
      type(star_info), pointer :: s
      type(auto_diff_real_star_order1), intent(in) :: P
      type(auto_diff_real_star_order1), intent(in) :: gradr
      type(auto_diff_real_star_order1), intent(out) :: gradr_spot

      real(dp) :: mu_ideal_gas, R2, Teff_local, PB_i
      type(auto_diff_real_star_order1) :: xspot_of_r  ! xspot4

      if (.not. s%doing_relax .and. .not. s%doing_first_model_of_run) then
         mu_ideal_gas = s%mu(1)  !1.00794d0 ! for hydrogen, 1 gram per mole
         R2 = pow2(s%R(1))
         Teff_local = pow(s%L(1)/(pi4*boltz_sigma*R2), 0.25d0)
         PB_i = (cgas*s%rho(1)/mu_ideal_gas)*(1.0d0 - s%xspot)*Teff_local
         xspot_of_r = (P - PB_i)/P
         gradr_spot = gradr/(s%fspot*pow4(xspot_of_r) + 1.0d0 - s%fspot)
      else
         gradr_spot = gradr
      end if
   end subroutine starspot_tweak_gradr

   subroutine starspot_tweak_PT(s)
      ! Purpose: saves the surface luminosity and adjusts it and the effective
      ! temperature to account for starspots.
      !
      ! This subroutine is called at the beginning of get_surf_PT()
      ! in hydro_vars.f90

      type(star_info), pointer :: s

      real(dp) :: alp

      alp = 1d0 - s%fspot + s%fspot*pow4(s%xspot)

      ! This is the surface-average value for luminosity
      L_init = s%L(1)

      ! Set the surface L to the unspotted, ambient L
      s%L(1) = s%L(1)/alp

      ! Now, set the Teff. Used in atm table lookup to set boundary conditions
      s%Teff = pow(s%L(1)/(pi4*pow2(s%r(1))*boltz_sigma), 0.25_dp)

   end subroutine starspot_tweak_PT

   subroutine starspot_restore_PT(s)
      ! Purpose: restores the surface luminosity and effective temperature.
      !
      ! Called at the end of get_surf_PT()

      type(star_info), pointer :: s

      s%Teff = pow(L_init/(pi4*pow2(s%r(1))*boltz_sigma), 0.25_dp)
      s%L(1) = L_init

   end subroutine starspot_restore_PT

end module starspots
