! ***********************************************************************
!
!   Copyright (C) 2010-2024  The MESA Team
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

module turb_def

   use const_def, only: dp, no_mixing

   implicit none

   type th_info_t ! Detailed thermohaline info
      real(dp) :: K_therm = 0._dp         ! Thermal conductivity
      real(dp) :: K_T = 0._dp             ! Thermal diffusivity
      real(dp) :: K_C = 0._dp             ! Chemical diffusivity
      real(dp) :: nu = 0._dp              ! Viscosity
      real(dp) :: Pr = 0._dp              ! Prandtl number
      real(dp) :: tau = 0._dp             ! Chemical diffusivity ratio
      real(dp) :: R_0 = 0._dp             ! Density ratio
      real(dp) :: r = 0._dp               ! Reduced density ratio
      real(dp) :: H_B = 0._dp             ! Lorentz force coefficient
      real(dp) :: Pm = 0._dp              ! Magnetic Prandtl number
      real(dp) :: D_B = 0._dp             ! Magnetic diffusivity
      real(dp) :: lam_hat = 0._dp         ! Growth rate of fastest-growing fingering mode
      real(dp) :: l2_hat = 0._dp          ! Horizontal wavenumber squared of fastest-growing fingering mode
      real(dp) :: sigma_max = 0._dp       ! Growth rate of fastest-growing parasitic mode
      real(dp) :: k_z_max = 0._dp         ! Vertical wavenumber of fastest-growing parasitic mode
      real(dp) :: w = 0._dp               ! Saturation flow speed
      real(dp) :: w_HG19 = 0._dp          ! Saturation flow speed in HG19 treatment
      real(dp) :: w_FRG24 = 0._dp         ! Saturation flow speed in FRG24 treatment
      real(dp) :: Nu_C = 0._dp            ! Compositional Nusselt number
      real(dp) :: D_thrm = 0._dp          ! Effective thermohaline mixing diffusivity
      integer  :: mixing_type = no_mixing ! Mixing type
   end type th_info_t

   private
   public :: th_info_t

end module turb_def
