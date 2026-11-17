! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
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

      module implicit_brunt

      use const_def, only: dp, crad, pi4
      use star_private_def
      use utils_lib, only: is_bad_num, mesa_error
      use auto_diff_support

      implicit none

      private
      public :: do_implicit_brunt_B
      public :: set_implicit_gradL_composition_term_ad

      contains

      subroutine do_implicit_brunt_B(s, nzlo, nzhi, ierr)
         use eos_def, only: num_eos_basic_results
         use star_utils, only: get_face_values
         use utils_lib, only: set_nan

         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr

         integer :: k, nz, op_err
         real(dp), allocatable, dimension(:) :: T_face, rho_face, xa_face, xa_path, dcomp_m1, dcomp_00
         real(dp), allocatable, dimension(:,:) :: d_eos_dxa

         include 'formats'

         ierr = 0
         nz = s% nz

         if (.not. s% calculate_Brunt_B) then
            call set_nan(s% brunt_B(1:nz))
            call set_nan(s% unsmoothed_brunt_B(1:nz))
            s% smoothed_brunt_B(1:nz) = 0d0
            s% gradL_composition_term(1:nz) = 0d0
            s% d_brunt_B_dxa_m1(1:s% species,1:nz) = 0d0
            s% d_brunt_B_dxa_00(1:s% species,1:nz) = 0d0
            do k = 1, nz
               s% brunt_B_ad(k) = 0d0
               s% gradL_composition_term_ad(k) = 0d0
            end do
            return
         end if

         allocate(T_face(nz), rho_face(nz))

         call get_face_values(s, s% T, T_face, ierr)
         if (ierr /= 0) return

         call get_face_values(s, s% rho, rho_face, ierr)
         if (ierr /= 0) return

!$OMP PARALLEL PRIVATE(k,op_err,xa_face,xa_path,dcomp_m1,dcomp_00,d_eos_dxa)
         allocate(xa_face(s% species), xa_path(s% species), dcomp_m1(s% species), dcomp_00(s% species), &
            d_eos_dxa(num_eos_basic_results,s% species))
!$OMP DO SCHEDULE(dynamic,2)
         do k = nzlo, nzhi
            op_err = 0
            call get_implicit_brunt_B( &
               s, k, T_face(k), rho_face(k), xa_face, xa_path, dcomp_m1, dcomp_00, d_eos_dxa, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END DO
         deallocate(xa_face, xa_path, dcomp_m1, dcomp_00, d_eos_dxa)
!$OMP END PARALLEL
         if (ierr /= 0) return

         do k = nzlo, nzhi
            s% unsmoothed_brunt_B(k) = s% brunt_B(k)
            if (is_bad_num(s% unsmoothed_brunt_B(k))) then
               write(*,2) 'unsmoothed_brunt_B(k)', k, s% unsmoothed_brunt_B(k)
               call mesa_error(__FILE__,__LINE__,'implicit_brunt')
            end if
         end do

         s% smoothed_brunt_B(1:nz) = s% unsmoothed_brunt_B(1:nz)

         if (s% use_Ledoux_criterion) then
            s% gradL_composition_term(1:nz) = s% brunt_B(1:nz)
         else
            s% gradL_composition_term(1:nz) = 0d0
         end if

         call set_implicit_gradL_composition_term_ad(s)

      end subroutine do_implicit_brunt_B


      subroutine set_implicit_gradL_composition_term_ad(s, use_brunt_ad)
         type (star_info), pointer :: s
         logical, intent(in), optional :: use_brunt_ad
         integer :: k
         logical :: use_ad

         use_ad = .true.
         if (present(use_brunt_ad)) use_ad = use_brunt_ad

         do k = 1, s% nz
            if (s% use_Ledoux_criterion) then
               if (use_ad) then
                  s% gradL_composition_term_ad(k) = s% brunt_B_ad(k)
               else
                  s% gradL_composition_term_ad(k) = &
                     s% gradL_composition_term(k)
               end if
               s% gradL_composition_term_ad(k)%val = &
                  s% gradL_composition_term(k)
               if (use_ad) &
                  s% gradL_composition_term_ad(k)%d1Array(i_xtra2_00) = 1d0
            else
               s% gradL_composition_term_ad(k) = 0d0
            end if
         end do

      end subroutine set_implicit_gradL_composition_term_ad


      subroutine get_implicit_brunt_B( &
            s, k, T_face, rho_face, xa_face, xa_path, dcomp_m1, dcomp_00, d_eos_dxa, ierr)
         use eos_def, only: num_eos_basic_results, i_lnPgas, i_chiRho, i_chiT
         use eos_support, only: get_eos_brunt_dxa_with_moments
         use chem_lib, only: basic_composition_info
         use chem_def, only: chem_isos
         use star_utils, only: get_ChiRho_face, get_ChiT_face, get_rho_face, get_T_face, get_Peos_face

         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: T_face, rho_face
         real(dp), intent(inout) :: xa_face(:), xa_path(:), dcomp_m1(:), dcomp_00(:)
         real(dp), intent(inout) :: d_eos_dxa(:,:)
         integer, intent(out) :: ierr

         integer :: j, species
         real(dp) :: alfa, beta, logRho_face, logT_face, Prad_face, Pgas_face, &
            Peos_face, delta_lnP, dlnP_dm, Ppoint, dlnP_composition, delta_lnMbar
         real(dp) :: lambda, Pgas_path, Peos_path, chiX_path
         real(dp) :: delta_lnP_val, chiT_val, brunt_B_val, &
            ddelta_lnP_m1, ddelta_lnP_00, dchiT_m1, dchiT_00, &
            dmass_corr_m1, dmass_corr_00
         real(dp) :: aion, wion, dlnmc_m1, dlnmc_00
         real(dp) :: X_tmp, Y_tmp, Z_tmp, abar_tmp, zbar_tmp, z2bar_tmp
         real(dp) :: z53bar_tmp, ye_tmp, mass_correction_tmp, sumx_tmp
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_eos_dlnd, d_eos_dlnT
         type(auto_diff_real_star_order1) :: &
            brunt_B_ad, delta_lnP_ad, chiRho_ad, chiT_ad, dlnP_composition_ad, &
            lnd_face_ad, lnT_face_ad, lnP_m1, lnP_00, Ppoint_ad, r_ad
         logical :: use_finite_difference, used_hydro_delta_lnP
         real(dp), parameter :: gauss_offset = 0.28867513459481288225d0
         real(dp), parameter :: gauss_weight = 0.5d0

         include 'formats'

         ierr = 0
         species = s% species
         s% brunt_B(k) = 0d0
         s% brunt_B_ad(k) = 0d0
         s% d_brunt_B_dxa_m1(1:species,k) = 0d0
         s% d_brunt_B_dxa_00(1:species,k) = 0d0
         if (k <= 1) return

         alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         beta = 1d0 - alfa
         xa_face(1:species) = alfa*s% xa(1:species,k) + beta*s% xa(1:species,k-1)

         logT_face = log10(T_face)
         logRho_face = log10(rho_face)
         Prad_face = crad*T_face*T_face*T_face*T_face/3d0

         if (is_bad_num(logT_face) .or. is_bad_num(logRho_face)) then
            ierr = -1
            return
         end if

         dlnP_composition = 0d0
         dcomp_00(1:species) = 0d0
         use_finite_difference = s% job% implicit_diffusion_use_brunt_finite_difference_value
         if (use_finite_difference) then
            lnd_face_ad = log(get_rho_face(s,k))
            lnT_face_ad = log(get_T_face(s,k))
            call eval_lnPeos(s% xa(1:species,k), lnP_00, dcomp_00)
            if (ierr /= 0) return
            call eval_lnPeos(s% xa(1:species,k-1), lnP_m1, dcomp_m1)
            if (ierr /= 0) return
            dlnP_composition_ad = lnP_00 - lnP_m1
            dcomp_m1(1:species) = -dcomp_m1(1:species)
         else if (s% job% implicit_diffusion_use_brunt_gauss_path) then
            lambda = 0.5d0 - gauss_offset
            xa_path(1:species) = (1d0 - lambda)*s% xa(1:species,k-1) + &
               lambda*s% xa(1:species,k)
            call basic_composition_info( &
               species, s% chem_id, xa_path, X_tmp, Y_tmp, Z_tmp, &
               abar_tmp, zbar_tmp, z2bar_tmp, z53bar_tmp, ye_tmp, &
               mass_correction_tmp, sumx_tmp)
            ! cheaper wrapper: path samples need only the lnPgas dxa row.
            call get_eos_brunt_dxa_with_moments( &
               s, 0, xa_path, rho_face, logRho_face, T_face, logT_face, .false., &
               X_tmp, Z_tmp, abar_tmp, zbar_tmp, z2bar_tmp, z53bar_tmp, ye_tmp, &
               mass_correction_tmp, sumx_tmp, &
               res, d_eos_dlnd, d_eos_dlnT, d_eos_dxa, ierr)
            if (ierr /= 0) return
            Pgas_path = exp(res(i_lnPgas))
            Peos_path = Prad_face + Pgas_path
            do j = 1, species
               chiX_path = (Pgas_path/Peos_path)*d_eos_dxa(i_lnPgas,j)
               dcomp_00(j) = dcomp_00(j) + gauss_weight*chiX_path
               dlnP_composition = dlnP_composition + gauss_weight*chiX_path* &
                  (s% xa(j,k) - s% xa(j,k-1))
            end do

            lambda = 0.5d0 + gauss_offset
            xa_path(1:species) = (1d0 - lambda)*s% xa(1:species,k-1) + &
               lambda*s% xa(1:species,k)
            call basic_composition_info( &
               species, s% chem_id, xa_path, X_tmp, Y_tmp, Z_tmp, &
               abar_tmp, zbar_tmp, z2bar_tmp, z53bar_tmp, ye_tmp, &
               mass_correction_tmp, sumx_tmp)
            ! cheaper wrapper: path samples need only the lnPgas dxa row.
            call get_eos_brunt_dxa_with_moments( &
               s, 0, xa_path, rho_face, logRho_face, T_face, logT_face, .false., &
               X_tmp, Z_tmp, abar_tmp, zbar_tmp, z2bar_tmp, z53bar_tmp, ye_tmp, &
               mass_correction_tmp, sumx_tmp, &
               res, d_eos_dlnd, d_eos_dlnT, d_eos_dxa, ierr)
            if (ierr /= 0) return
            Pgas_path = exp(res(i_lnPgas))
            Peos_path = Prad_face + Pgas_path
            do j = 1, species
               chiX_path = (Pgas_path/Peos_path)*d_eos_dxa(i_lnPgas,j)
               dcomp_00(j) = dcomp_00(j) + gauss_weight*chiX_path
               dlnP_composition = dlnP_composition + gauss_weight*chiX_path* &
                  (s% xa(j,k) - s% xa(j,k-1))
            end do
         end if

         call basic_composition_info( &
            species, s% chem_id, xa_face, X_tmp, Y_tmp, Z_tmp, &
            abar_tmp, zbar_tmp, z2bar_tmp, z53bar_tmp, ye_tmp, &
            mass_correction_tmp, sumx_tmp)
         ! cheaper wrapper: the face value needs only brunt's dxa rows.
         call get_eos_brunt_dxa_with_moments( &
            s, 0, xa_face, rho_face, logRho_face, T_face, logT_face, &
            s% job% implicit_diffusion_include_dsig_dxa, &
            X_tmp, Z_tmp, abar_tmp, zbar_tmp, z2bar_tmp, z53bar_tmp, ye_tmp, &
            mass_correction_tmp, sumx_tmp, &
            res, d_eos_dlnd, d_eos_dlnT, d_eos_dxa, ierr)
         if (ierr /= 0) return

         if (.not. use_finite_difference) then
            if (.not. s% job% implicit_diffusion_use_brunt_gauss_path) then
               Pgas_face = exp(res(i_lnPgas))
               Peos_face = Prad_face + Pgas_face
               do j = 1, species
                  dcomp_00(j) = (Pgas_face/Peos_face)*d_eos_dxa(i_lnPgas,j)
                  dlnP_composition = dlnP_composition + dcomp_00(j)*(s% xa(j,k) - s% xa(j,k-1))
               end do
            end if
            ! Hold the face/path pressure partials fixed for the contraction.
            dcomp_m1(1:species) = -dcomp_00(1:species)
            dlnP_composition_ad = dlnP_composition
         end if

         used_hydro_delta_lnP = .false.
         delta_lnP_ad = wrap_lnPeos_m1(s,k) - wrap_lnPeos_00(s,k)
         delta_lnP = delta_lnP_ad% val
         if (delta_lnP > -1d-50) then
            used_hydro_delta_lnP = .true.
            Ppoint = alfa*s% Peos(k) + beta*s% Peos(k-1)
            dlnP_dm = -s% cgrav(k)*s% m(k)/(pi4*s% r(k)*s% r(k)*s% r(k)*s% r(k)*Ppoint)
            delta_lnP = dlnP_dm*s% dm_bar(k)
            delta_lnP_ad = delta_lnP
            if (use_finite_difference) then
               Ppoint_ad = get_Peos_face(s,k)
               r_ad = wrap_r_00(s,k)
               delta_lnP_ad% d1Array = -delta_lnP*(4d0*r_ad% d1Array/r_ad% val + Ppoint_ad% d1Array/Ppoint)
            end if
         end if

         if (use_finite_difference) then
            chiT_ad = res(i_chiT)
            chiT_ad% d1Array = d_eos_dlnd(i_chiT)*lnd_face_ad% d1Array + &
               d_eos_dlnT(i_chiT)*lnT_face_ad% d1Array
         else
            chiT_ad = get_ChiT_face(s,k)
            chiT_ad% val = res(i_chiT)
         end if
         brunt_B_ad = dlnP_composition_ad/(delta_lnP_ad*chiT_ad)

         if (s% use_mass_corrections) then
            if (use_finite_difference) then
               chiRho_ad = res(i_chiRho)
               chiRho_ad% d1Array = d_eos_dlnd(i_chiRho)*lnd_face_ad% d1Array + &
                  d_eos_dlnT(i_chiRho)*lnT_face_ad% d1Array
            else
               chiRho_ad = get_ChiRho_face(s,k)
               chiRho_ad% val = res(i_chiRho)
            end if
            delta_lnMbar = log(s% mass_correction(k-1)) - log(s% mass_correction(k))
            brunt_B_ad = brunt_B_ad - &
               chiRho_ad*delta_lnMbar/(delta_lnP_ad*chiT_ad)
         end if

         s% brunt_B_ad(k) = brunt_B_ad
         s% brunt_B(k) = brunt_B_ad% val
         delta_lnP_val = delta_lnP_ad% val
         chiT_val = chiT_ad% val
         brunt_B_val = brunt_B_ad% val

         if (s% job% implicit_diffusion_include_dsig_dxa) then
            do j = 1, species
               dchiT_m1 = beta*d_eos_dxa(i_chiT,j)
               dchiT_00 = alfa*d_eos_dxa(i_chiT,j)
               if (used_hydro_delta_lnP) then
                  if (use_finite_difference) then
                     ddelta_lnP_m1 = -delta_lnP_val*beta*s% Peos(k-1)/Ppoint*s% dlnPeos_dxa_for_partials(j,k-1)
                     ddelta_lnP_00 = -delta_lnP_val*alfa*s% Peos(k)/Ppoint*s% dlnPeos_dxa_for_partials(j,k)
                  else
                     ddelta_lnP_m1 = 0d0
                     ddelta_lnP_00 = 0d0
                  end if
               else
                  ddelta_lnP_m1 = s% dlnPeos_dxa_for_partials(j,k-1)
                  ddelta_lnP_00 = -s% dlnPeos_dxa_for_partials(j,k)
               end if
               dmass_corr_m1 = 0d0
               dmass_corr_00 = 0d0
               if (s% use_mass_corrections) then
                  aion = dble(chem_isos% Z_plus_N(s% chem_id(j)))
                  wion = chem_isos% W(s% chem_id(j))
                  dlnmc_m1 = (wion/aion - s% mass_correction(k-1))/ &
                     s% mass_correction(k-1)
                  dlnmc_00 = (wion/aion - s% mass_correction(k))/ &
                     s% mass_correction(k)
                  dmass_corr_m1 = delta_lnMbar*beta*d_eos_dxa(i_chiRho,j) + &
                     chiRho_ad% val*dlnmc_m1
                  dmass_corr_00 = delta_lnMbar*alfa*d_eos_dxa(i_chiRho,j) - &
                     chiRho_ad% val*dlnmc_00
               end if
               s% d_brunt_B_dxa_m1(j,k) = &
                  (dcomp_m1(j) - dmass_corr_m1)/(delta_lnP_val*chiT_val) - &
                  brunt_B_val*(ddelta_lnP_m1/delta_lnP_val + dchiT_m1/chiT_val)
               s% d_brunt_B_dxa_00(j,k) = &
                  (dcomp_00(j) - dmass_corr_00)/(delta_lnP_val*chiT_val) - &
                  brunt_B_val*(ddelta_lnP_00/delta_lnP_val + dchiT_00/chiT_val)
            end do
         end if

         if (is_bad_num(s% brunt_B(k))) then
            ierr = -1
            s% retry_message = 'bad num for implicit brunt_B'
            if (s% report_ierr) then
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               write(*,2) 'delta_lnP', k, delta_lnP
               write(*,2) 'res(i_chiT)', k, res(i_chiT)
            end if
         end if

         contains

         subroutine eval_lnPeos(xa, lnPeos, dlnPeos_dxa)
            real(dp), intent(in) :: xa(:)
            type(auto_diff_real_star_order1), intent(out) :: lnPeos
            real(dp), intent(out) :: dlnPeos_dxa(:)
            real(dp) :: Pgas, Peos

            call basic_composition_info( &
               species, s% chem_id, xa, X_tmp, Y_tmp, Z_tmp, &
               abar_tmp, zbar_tmp, z2bar_tmp, z53bar_tmp, ye_tmp, &
               mass_correction_tmp, sumx_tmp)
            ! Evaluate the pressure and its partials at the same endpoint.
            call get_eos_brunt_dxa_with_moments( &
               s, 0, xa, rho_face, logRho_face, T_face, logT_face, .false., &
               X_tmp, Z_tmp, abar_tmp, zbar_tmp, z2bar_tmp, z53bar_tmp, ye_tmp, &
               mass_correction_tmp, sumx_tmp, &
               res, d_eos_dlnd, d_eos_dlnT, d_eos_dxa, ierr)
            if (ierr /= 0) return
            Pgas = exp(res(i_lnPgas))
            Peos = Prad_face + Pgas
            lnPeos = log(Peos)
            ! Convert the gas pressure rows to derivatives of total pressure.
            lnPeos% d1Array = (Pgas*d_eos_dlnd(i_lnPgas)/Peos)*lnd_face_ad% d1Array + &
               ((Pgas*d_eos_dlnT(i_lnPgas) + 4d0*Prad_face)/Peos)*lnT_face_ad% d1Array
            dlnPeos_dxa(1:species) = (Pgas/Peos)*d_eos_dxa(i_lnPgas,1:species)
         end subroutine eval_lnPeos

      end subroutine get_implicit_brunt_B

      end module implicit_brunt
