! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

      module brunt

      use star_private_def
      use const_def
      use utils_lib

      implicit none

      private
      public :: do_brunt_B, do_brunt_N2

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
         integer :: nz, k, i
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

         integer :: nz, k, i

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

         do k=1,nz ! clip B and calculate N^2 from B
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
         real(dp) :: brunt_B
         integer :: nz, species, k, i, op_err
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
               !call mesa_error(__FILE__,__LINE__,'do_brunt_B_MHM_form')
            end if
            if (s% stop_for_bad_nums) then
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               call mesa_error(__FILE__,__LINE__,'do_brunt_B_MHM_form')
            end if
         end if

      end subroutine get_brunt_B


      end module brunt

