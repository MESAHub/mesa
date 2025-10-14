! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

module atm_def

  use const_def, only: dp

  implicit none

  ! T-tau relations

  integer, parameter :: ATM_T_TAU_EDDINGTON = 1
  integer, parameter :: ATM_T_TAU_SOLAR_HOPF = 2
  integer, parameter :: ATM_T_TAU_KRISHNA_SWAMY = 3
  integer, parameter :: ATM_T_TAU_TRAMPEDACH_SOLAR = 4

  ! Tables

  integer, parameter :: ATM_TABLE_TAU_100 = 100
  integer, parameter :: ATM_TABLE_TAU_10 = 101
  integer, parameter :: ATM_TABLE_TAU_1 = 102
  integer, parameter :: ATM_TABLE_TAU_1M1 = 103
  integer, parameter :: ATM_TABLE_PHOTOSPHERE = 104
  integer, parameter :: ATM_TABLE_WD_TAU_25 = 105
  integer, parameter :: ATM_TABLE_DB_WD_TAU_25 = 106

  integer, parameter :: table_atm_version = 5

  ! Atmosphere structure info

  integer, parameter :: atm_xm = 1  ! mass of atm exterior to this point (g)
  integer, parameter :: atm_delta_r = atm_xm+1  ! radial distance above base of envelope (cm)
  integer, parameter :: atm_lnP = atm_delta_r+1
  integer, parameter :: atm_lnd = atm_lnP+1
  integer, parameter :: atm_lnT = atm_lnd+1
  integer, parameter :: atm_gradT = atm_lnT+1
  integer, parameter :: atm_kap = atm_gradT+1
  integer, parameter :: atm_gamma1 = atm_kap+1
  integer, parameter :: atm_grada = atm_gamma1+1
  integer, parameter :: atm_chiT = atm_grada+1
  integer, parameter :: atm_chiRho = atm_chiT+1
  integer, parameter :: atm_cv = atm_chiRho+1
  integer, parameter :: atm_cp = atm_cv+1
  integer, parameter :: atm_lnfree_e = atm_cp+1
  integer, parameter :: atm_dlnkap_dlnT = atm_lnfree_e+1
  integer, parameter :: atm_dlnkap_dlnd = atm_dlnkap_dlnT+1
  integer, parameter :: atm_lnPgas = atm_dlnkap_dlnd+1
  integer, parameter :: atm_tau = atm_lnPgas+1
  integer, parameter :: atm_gradr = atm_tau+1

  integer, parameter :: num_results_for_build_atm = atm_gradr

  ! Derived-type definitions

  type Atm_Info
     integer :: id
     integer :: nZ
     integer :: ng
     integer :: nT
     integer :: ilinT
     integer :: iling
     real(dp), pointer          :: Teff_array(:)
     real(dp), pointer          :: logg_array(:)
     real(dp), pointer          :: Teff_bound(:)
     real(dp), pointer          :: logZ(:)
     real(dp), pointer          :: alphaFe(:)
     real(dp), pointer          :: Pgas_interp1(:)
     real(dp), pointer          :: T_interp1(:)
     real(dp), pointer          :: Pgas_interp(:,:,:,:)
     real(dp), pointer          :: T_interp(:,:,:,:)
     character(len=8), pointer  :: atm_mix(:)
     character(len=40), pointer :: table_atm_files(:)
     logical, pointer           :: have_atm_table(:)
  end type Atm_Info

  ! Atmosphere tables

  type (Atm_Info), target, save :: &
       ai_two_thirds, ai_100, ai_10, ai_1, &
       ai_1m1, ai_wd_25, ai_db_wd_25

  logical :: table_atm_is_initialized = .false.

  logical :: star_debugging_atm_flag = .false.
  real(dp) :: atm_test_partials_val, atm_test_partials_dval_dx
  real(dp) :: atm_test_partials_L_lo, atm_test_partials_L_hi
  real(dp) :: atm_test_partials_R_lo, atm_test_partials_R_hi
  real(dp) :: atm_test_partials_M_lo, atm_test_partials_M_hi

  abstract interface

     ! Callback routine for EOS evaluation

     subroutine atm_eos_iface( &
          lnP, lnT, &
          lnRho, res, dres_dlnRho, dres_dlnT, &
          ierr)
       use const_def, only: dp
       implicit none
       real(dp), intent(in)  :: lnP
       real(dp), intent(in)  :: lnT
       real(dp), intent(out) :: lnRho
       real(dp), intent(out) :: res(:)
       real(dp), intent(out) :: dres_dlnRho(:)
       real(dp), intent(out) :: dres_dlnT(:)
       integer, intent(out)  :: ierr
     end subroutine atm_eos_iface

     ! Callback routine for opacity evaluation

     subroutine atm_kap_iface( &
          lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
          kap, dlnkap_dlnRho, dlnkap_dlnT, &
          ierr)
       use const_def, only: dp
       implicit none
       real(dp), intent(in)  :: lnRho
       real(dp), intent(in)  :: lnT
       real(dp), intent(in)  :: res(:)
       real(dp), intent(in)  :: dres_dlnRho(:)
       real(dp), intent(in)  :: dres_dlnT(:)
       real(dp), intent(out) :: kap
       real(dp), intent(out) :: dlnkap_dlnRho
       real(dp), intent(out) :: dlnkap_dlnT
       integer, intent(out)  :: ierr
     end subroutine atm_kap_iface

  end interface

end module atm_def
