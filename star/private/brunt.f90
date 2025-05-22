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

      module brunt

      use star_private_def
      use const_def, only: dp, pi4, crad
      use utils_lib

      implicit none

      private
      public :: do_brunt_B
      public :: do_brunt_N2

      logical, parameter :: dbg = .false.

      contains

      ! call this before mlt
      subroutine do_brunt_B(s,nzlo,nzhi,ierr)
         use star_utils, only: get_face_values
         use utils_lib, only: set_nan
         use interp_1d_def

         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: nz, k
         real(dp), allocatable, dimension(:) :: smoothing_array

         include 'formats'

         ierr = 0
         nz = s% nz

         if (.not. s% calculate_Brunt_B) then
            call set_nan(s% brunt_B(1:nz))
            call set_nan(s% unsmoothed_brunt_B(1:nz))
            return
         end if

         if (s% use_other_brunt) then
            call s% other_brunt(s% id, ierr)
            if (ierr /= 0) then
                s% retry_message = 'failed in other_brunt'
                if (s% report_ierr) write(*, *) s% retry_message
            end if

         else if (s% use_eos_partials_for_Brunt) then
            ! ------------------------------------------------------------
            ! chain-rule approach for calculating B from eos partials, see MESA II equation 6.
            ! ------------------------------------------------------------
            call do_brunt_B_eos_partials_form(s, nzlo, nzhi, ierr)
            if (ierr /= 0) then
               s% retry_message = 'failed in do_brunt_B_eos_partials_form'
               if (s% report_ierr) write(*, *) s% retry_message
            end if

         else
            call do_brunt_B_MHM_form(s,nzlo,nzhi,ierr)
            if (ierr /= 0) then
                s% retry_message = 'failed in do_brunt_B_MHM_form'
                if (s% report_ierr) write(*, *) s% retry_message
            end if
         end if
         if (ierr /= 0) return

         ! save unsmoothed brunt_B
         do k=nzlo,nzhi
            s% unsmoothed_brunt_B(k) = s% brunt_B(k)
            if (is_bad(s% unsmoothed_brunt_B(k))) then
               write(*,2) 'unsmoothed_brunt_B(k)', k, s% unsmoothed_brunt_B(k)
               call mesa_error(__FILE__,__LINE__,'brunt')
            end if
         end do

         allocate(smoothing_array(nz))
         call smooth_brunt_B(smoothing_array)
         if (s% use_other_brunt_smoothing) then
            call s% other_brunt_smoothing(s% id, ierr)
            if (ierr /= 0) then
               s% retry_message = 'failed in other_brunt_smoothing'
               if (s% report_ierr) write(*, *) s% retry_message
               return
            end if
         end if

         contains

         subroutine smooth_brunt_B(work)
            use star_utils, only: threshold_smoothing
            real(dp) :: work(:)
            logical, parameter :: preserve_sign = .false.
            if (s% num_cells_for_smooth_brunt_B <= 0) return
            call threshold_smoothing( &
               s% brunt_B, s% threshold_for_smooth_brunt_B, s% nz, &
               s% num_cells_for_smooth_brunt_B, preserve_sign, work)
         end subroutine smooth_brunt_B

      end subroutine do_brunt_B


      ! call this after mlt
      subroutine do_brunt_N2(s,nzlo,nzhi,ierr)
         use star_utils, only: get_face_values

         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr

         real(dp) :: f
         real(dp), allocatable, dimension(:) :: &
            rho_P_chiT_chiRho, rho_P_chiT_chiRho_face

         integer :: nz, k

         include 'formats'

         ierr = 0
         nz = s% nz

         if (.not. (s% calculate_Brunt_B .and. s% calculate_Brunt_N2)) then
            call set_nan(s% brunt_N2(1:nz))
            call set_nan(s% brunt_N2_composition_term(1:nz))
            return
         end if

         allocate(rho_P_chiT_chiRho(nz), rho_P_chiT_chiRho_face(nz))

         do k=1,nz
            rho_P_chiT_chiRho(k) = (s% rho(k)/s% Peos(k))*(s% chiT(k)/s% chiRho(k))
            ! correct for difference between gravitational mass density and baryonic mass density (rho)
            if (s% use_mass_corrections) then
               rho_P_chiT_chiRho(k) = s% mass_correction(k)*rho_P_chiT_chiRho(k)
            end if
         end do

         call get_face_values( &
            s, rho_P_chiT_chiRho, rho_P_chiT_chiRho_face, ierr)
         if (ierr /= 0) then
            s% retry_message = 'failed in get_face_values'
            if (s% report_ierr) write(*, *) s% retry_message
            return
         end if

         do k=1,nz  ! clip B and calculate N^2 from B
            if (abs(s% brunt_B(k)) < s% min_magnitude_brunt_B .or. &
                  s% gradT(k) == 0 .or. is_bad(s% gradT_sub_grada(k))) then
               s% brunt_B(k) = 0
               s% brunt_N2(k) = 0
               s% brunt_N2_composition_term(k) = 0
               cycle
            end if
            f = pow2(s% grav(k))*rho_P_chiT_chiRho_face(k)
            if (is_bad(f) .or. is_bad(s% brunt_B(k)) .or. is_bad(s% gradT_sub_grada(k))) then
               write(*,2) 'f', k, f
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               write(*,2) 's% gradT_sub_grada(k)', k, s% gradT_sub_grada(k)
               write(*,2) 's% gradT(k)', k, s% gradT(k)
               write(*,2) 's% grada_face(k)', k, s% grada_face(k)
               call mesa_error(__FILE__,__LINE__,'brunt')
            end if
            s% brunt_N2(k) = f*(s% brunt_B(k) - s% gradT_sub_grada(k))
            s% brunt_N2_composition_term(k) = f*s% brunt_B(k)
         end do

         if (s% brunt_N2_coefficient /= 1d0) then
            do k=1,nz
               s% brunt_N2(k) = s% brunt_N2_coefficient*s% brunt_N2(k)
               s% brunt_N2_composition_term(k) = &
                  s% brunt_N2_coefficient*s% brunt_N2_composition_term(k)
            end do
         end if

      end subroutine do_brunt_N2


      subroutine do_brunt_B_MHM_form(s, nzlo, nzhi, ierr)
         ! Brassard from Mike Montgomery (MHM)
         use star_utils, only: get_face_values
         use interp_1d_def

         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr


         real(dp), allocatable, dimension(:) :: T_face, rho_face, chiT_face, chiRho_face
         integer :: nz, species, k, op_err
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         nz = s% nz
         species = s% species

         allocate(T_face(nz), rho_face(nz), chiT_face(nz), chiRho_face(nz))

         call get_face_values(s, s% chiT, chiT_face, ierr)
         if (ierr /= 0) return

         call get_face_values(s, s% chiRho, chiRho_face, ierr)
         if (ierr /= 0) return

         call get_face_values(s, s% T, T_face, ierr)
         if (ierr /= 0) return

         call get_face_values(s, s% rho, rho_face, ierr)
         if (ierr /= 0) return

!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k=nzlo,nzhi
            op_err = 0
            call get_brunt_B(&
               s, species, nz, k, T_face(k), rho_face(k), chiT_face(k), chiRho_face(k), op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO

      end subroutine do_brunt_B_MHM_form


      subroutine get_brunt_B(s, species, nz, k, T_face, rho_face, chiT_face, chiRho_face, ierr)
         use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_lnPgas
         use eos_support, only: get_eos

         type (star_info), pointer :: s
         integer, intent(in) :: species, nz, k
         real(dp), intent(in) :: T_face, rho_face, chiT_face, chiRho_face
         integer, intent(out) :: ierr

         real(dp) :: lnP1, lnP2, logRho_face, logT_face, Prad_face, &
            alfa, Ppoint, dlnP_dm, delta_lnP, delta_lnMbar
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_eos_dlnd, d_eos_dlnT
         real(dp) :: d_eos_dxa(num_eos_d_dxa_results,species)

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         s% brunt_B(k) = 0
         if (k <= 1) return

         logT_face = log10(T_face)
         logRho_face = log10(rho_face)
         Prad_face = crad * T_face*T_face*T_face*T_face / 3

         if (is_bad_num(logT_face) .or. is_bad_num(logRho_face)) then
            ierr = -1
            return
            write(*,2) 'logT_face', k, logT_face
            write(*,2) 'logRho_face', k, logRho_face
            write(*,'(A)')
            write(*,2) 's% dq(k-1)', k-1, s% dq(k-1)
            write(*,2) 's% dq(k)', k, s% dq(k)
            write(*,'(A)')
            write(*,2) 's% T(k-1)', k-1, s% T(k-1)
            write(*,2) 'T_face', k, T_face
            write(*,2) 's% T(k)', k, s% T(k)
            write(*,'(A)')
            write(*,2) 's% rho(k-1)', k-1, s% rho(k-1)
            write(*,2) 'rho_face', k, rho_face
            write(*,2) 's% rho(k)', k, s% rho(k)
            write(*,'(A)')
            write(*,2) 's% chiT(k-1)', k-1, s% chiT(k-1)
            write(*,2) 'chiT_face', k, chiT_face
            write(*,2) 's% chiT(k)', k, s% chiT(k)
            write(*,'(A)')
            call mesa_error(__FILE__,__LINE__,'get_brunt_B')
         end if

         call get_eos( &
            s, 0, s% xa(:,k), &
            rho_face, logRho_face, T_face, logT_face, &
            res, d_eos_dlnd, d_eos_dlnT, &
            d_eos_dxa, ierr)
         if (ierr /= 0) return
         lnP1 = log(Prad_face + exp(res(i_lnPgas)))

         call get_eos( &
            s, 0, s% xa(:,k-1), &
            rho_face, logRho_face, T_face, logT_face, &
            res, d_eos_dlnd, d_eos_dlnT, &
            d_eos_dxa, ierr)
         if (ierr /= 0) return
         lnP2 = log(Prad_face + exp(res(i_lnPgas)))

         delta_lnP = s% lnPeos(k-1) - s% lnPeos(k)
         if (delta_lnP > -1d-50) then
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            Ppoint = alfa*s% Peos(k) + (1-alfa)*s% Peos(k-1)
            dlnP_dm = -s% cgrav(k)*s% m(k)/(pi4*pow4(s% r(k))*Ppoint)
            delta_lnP = dlnP_dm*s% dm_bar(k)
         end if

         s% brunt_B(k) = (lnP1 - lnP2)/delta_lnP/chiT_face

         ! add term accounting for the composition-related gradient in gravitational mass
         if (s% use_mass_corrections) then
            delta_lnMbar = log(s% mass_correction(k-1)) - log(s% mass_correction(k))
            s% brunt_B(k) = s% brunt_B(k) - chiRho_face*delta_lnMbar/delta_lnP/chiT_face
         end if

         if (is_bad_num(s% brunt_B(k))) then
            ierr = -1
            s% retry_message = 'bad num for brunt_B'
            if (s% report_ierr) then
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               write(*,2) 'chiT_face', k, chiT_face
               write(*,2) 'delta_lnP', k, delta_lnP
               write(*,2) 's% lnPeos(k)', k, s% lnPeos(k)
               write(*,2) 's% lnPeos(k-1)', k-1, s% lnPeos(k-1)
               write(*,2) 'lnP1', k, lnP1
               write(*,2) 'lnP2', k, lnP2
               write(*,'(A)')
            end if
            if (s% stop_for_bad_nums) then
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               call mesa_error(__FILE__,__LINE__,'do_brunt_B_MHM_form')
            end if
         end if

      end subroutine get_brunt_B


      subroutine do_brunt_B_eos_partials_form(s, nzlo, nzhi, ierr)
         use star_utils, only: get_face_values
         use interp_1d_def

         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr


         real(dp), allocatable, dimension(:) :: T_face, rho_face, chiT_face, chiRho_face
      !   real(dp) :: mass_corr_factor, delta_lnP, delta_lnMbar, B_cell_centered
         real(dp), allocatable, dimension(:,:) :: xa_face! (species, k)
         integer :: nz, species, k, op_err
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         nz = s% nz
         species = s% species

         allocate(T_face(nz), rho_face(nz), chiT_face(nz), chiRho_face(nz), xa_face(species, nz)) ! ,B_cell_centered(nz)

         ! can we paralleize this?
         do k = 1, species
            call get_face_values(s, s% xa(k, :), xa_face(k, :), ierr)
            if (ierr /= 0) return
         end do
         !!!

         call get_face_values(s, s% chiT, chiT_face, ierr)
         if (ierr /= 0) return

         call get_face_values(s, s% chiRho, chiRho_face, ierr)
         if (ierr /= 0) return

         call get_face_values(s, s% T, T_face, ierr)
         if (ierr /= 0) return

         call get_face_values(s, s% rho, rho_face, ierr)
         if (ierr /= 0) return

      !$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k=nzlo,nzhi
            op_err = 0
            call get_brunt_B_from_eos_partials(&
               s, species, nz, k, T_face(k), rho_face(k), chiT_face(k), chiRho_face(k), xa_face(:,:), op_err)
            if (op_err /= 0) ierr = op_err
         end do
      !$OMP END PARALLEL DO


      end subroutine do_brunt_B_eos_partials_form

      subroutine get_brunt_B_from_eos_partials(s, species, nz, k, T_face, rho_face, chiT_face, chiRho_face, xa_face, ierr)
         use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_lnPgas
         use eos_support, only: get_eos

         type (star_info), pointer :: s
         integer, intent(in) :: species, nz, k
         real(dp), intent(in) :: T_face, rho_face, chiT_face, chiRho_face
         integer, intent(out) :: ierr

         real(dp) :: logRho_face, logT_face, Prad_face
         real(dp),  dimension(species,nz), intent(in):: xa_face! (species, nz)
         real(dp), dimension(num_eos_basic_results) :: res, d_eos_dlnd, d_eos_dlnT
         real(dp), dimension(num_eos_d_dxa_results, species) :: d_eos_dxa
         real(dp) :: delta_lnMbar, Ppoint, dlnP_dm, alfa
         real(dp) :: B_term, spatial_derivative_dX_dlnP, chiT
         real(dp) :: delta_lnP_cell, comp, y, t
         integer :: i

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         s% brunt_B(k) = 0
         if (k <= 1) return ! should we add a check for nz?

         logT_face = log10(T_face)
         logRho_face = log10(rho_face)
         Prad_face = crad * pow4(T_face) / 3d0

         if (is_bad_num(logT_face) .or. is_bad_num(logRho_face)) then
            ierr = -1
            return
         end if

         ! Call the EOS to get the required partial derivatives
         call get_eos( &
            s, 0, xa_face(:,k), &
            rho_face, logRho_face, T_face, logT_face, &
            res, d_eos_dlnd, d_eos_dlnT, &
            d_eos_dxa, ierr)
         if (ierr /= 0) return

         ! Extract chiT
         chiT = d_eos_dlnT(i_lnPgas)

         ! Initialize B_term to accumulate the contribution from each species
         B_term = 0d0
         comp   = 0d0          ! compensator

        ! Compute pressure difference across adjacent cells
        delta_lnP_cell  = s%lnPeos(k-1) - s%lnPeos(k) ! center difference
        if (abs(delta_lnP_cell) < tiny(1d0)) then ! is this an okay check?
            ! if the face pressure is flat, then the derivative is numerically zero
            s% brunt_B(k) = 0d0
            return
        end if

         ! Compute the Brunt B composition term
!        do i = 1, species
!            !X_face_p1 = xa_face(i, k)      ! face between zones k and k+1
!            !X_face_m1 = xa_face(i, k-1)    ! face between zones k-1 and k
!            spatial_derivative_dX_dlnP = (s%xa(i,k-1) - s%xa(i,k)) / delta_lnP_cell
!            if (abs(spatial_derivative_dX_dlnP) < 1d-12) spatial_derivative_dX_dlnP = 0d0
!            B_term = B_term - d_eos_dxa(i_lnPgas, i) * spatial_derivative_dX_dlnP
!        end do


         do i = 1, species
             spatial_derivative_dX_dlnP = (s%xa(i,k-1) - s%xa(i,k)) / delta_lnP_cell
             if (abs(spatial_derivative_dX_dlnP) < 1d-12) spatial_derivative_dX_dlnP = 0d0
             y = -d_eos_dxa(i_lnPgas,i) * spatial_derivative_dX_dlnP - comp
             t = B_term + y
             comp = (t - B_term) - y
             B_term = t
         end do

         ! Final calculation of B using chiT
         s% brunt_B(k) = B_term / chiT

        ! this block is always executed and assume HSE, is this an error?
         if (delta_lnP_cell > -1d-50) then
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            Ppoint = alfa*s% Peos(k) + (1-alfa)*s% Peos(k-1)
            dlnP_dm = -s% cgrav(k)*s% m(k)/(pi4*pow4(s% r(k))*Ppoint)
            delta_lnP_cell = dlnP_dm*s% dm_bar(k)
         end if

         ! add term accounting for the composition-related gradient in gravitational mass
         if (s% use_mass_corrections) then
            delta_lnMbar = log(s% mass_correction(k-1)) - log(s% mass_correction(k))
            s% brunt_B(k) = s% brunt_B(k) - chiRho_face*delta_lnMbar/delta_lnP_cell/chiT_face
         end if

         ! Check for bad numbers
         if (is_bad_num(s% brunt_B(k))) then
            ierr = -1
            s% retry_message = 'bad num for brunt_B'
            if (s% report_ierr) then
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               write(*,2) 'chiT', k, chiT
               write(*,2) 'B_term', k, B_term
               call mesa_error(__FILE__, __LINE__, 'get_brunt_B_from_eos_partials')
            end if
         end if
      end subroutine get_brunt_B_from_eos_partials



      end module brunt
