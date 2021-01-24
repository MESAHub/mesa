! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
         else if (s% use_brunt_gradmuX_form) then
            call do_brunt_B_gradmuX_form(s,ierr)
            if (ierr /= 0) then
                s% retry_message = 'failed in do_brunt_B_gradmuX_form'
                if (s% report_ierr) write(*, *) s% retry_message
            end if
         else
            call do_brunt_B_MHM_form(s,ierr)
            if (ierr /= 0) then
                s% retry_message = 'failed in do_brunt_B_MHM_form'
                if (s% report_ierr) write(*, *) s% retry_message
            end if
         end if
         if (ierr /= 0) return

         ! save unsmoothed brunt_B
         do k=1,nz
            s% unsmoothed_brunt_B(k) = s% brunt_B(k)
            if (is_bad(s% unsmoothed_brunt_B(k))) then
               write(*,2) 'unsmoothed_brunt_B(k)', k, s% unsmoothed_brunt_B(k)
               stop 'brunt'
            end if
         end do

      end subroutine do_brunt_B


      ! call this after mlt
      subroutine do_brunt_N2(s,nzlo,nzhi,ierr)
         use star_utils, only: get_face_values
         use interp_1d_def

         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr

         real(dp) :: f
         real(dp), pointer, dimension(:) :: &
            work1, f1, rho_P_chiT_chiRho, rho_P_chiT_chiRho_face
         integer, parameter :: nwork = pm_work_size

         integer :: nz, k, i

         include 'formats'

         ierr = 0
         nz = s% nz

         if (.not. (s% calculate_Brunt_B .and. s% calculate_Brunt_N2)) then
            call set_nan(s% brunt_N2(1:nz))
            call set_nan(s% brunt_N2_composition_term(1:nz))
            return
         end if

         call do_alloc(ierr)
         if (ierr /= 0) return

         call smooth_brunt_B(work1)
         if (s% use_other_brunt_smoothing) then
            call s% other_brunt_smoothing(s% id, ierr)
            if (ierr /= 0) then
               s% retry_message = 'failed in other_brunt_smoothing'
               if (s% report_ierr) write(*, *) s% retry_message
               call dealloc
               return
            end if
         end if

         do k=1,nz
            rho_P_chiT_chiRho(k) = (s% rho(k)/s% P(k))*(s% chiT(k)/s% chiRho(k))
         end do
         call get_face_values( &
            s, rho_P_chiT_chiRho, rho_P_chiT_chiRho_face, ierr)
         if (ierr /= 0) then
            s% retry_message = 'failed in get_face_values'
            if (s% report_ierr) write(*, *) s% retry_message
            call dealloc
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
            f = s% grav(k)*s% grav(k)*rho_P_chiT_chiRho_face(k)
            if (is_bad(f) .or. is_bad(s% brunt_B(k)) .or. is_bad(s% gradT_sub_grada(k))) then
               write(*,2) 'f', k, f
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               write(*,2) 's% gradT_sub_grada(k)', k, s% gradT_sub_grada(k)
               write(*,2) 's% gradT(k)', k, s% gradT(k)
               write(*,2) 's% grada_face(k)', k, s% grada_face(k)
               stop 'brunt'
            end if
            s% brunt_N2(k) = f*(s% brunt_B(k) - s% gradT_sub_grada(k))
            s% brunt_N2_composition_term(k) = f*s% brunt_B(k)
         end do

         call dealloc

         if (s% brunt_N2_coefficient /= 1d0) then
            do k=1,nz
               s% brunt_N2(k) = s% brunt_N2_coefficient*s% brunt_N2(k)
               s% brunt_N2_composition_term(k) = &
                  s% brunt_N2_coefficient*s% brunt_N2_composition_term(k)
            end do
         end if


         contains

         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            call do_work_arrays(.true., ierr)
         end subroutine do_alloc

         subroutine dealloc
            call do_work_arrays(.false., ierr)
         end subroutine dealloc
            
         subroutine do_work_arrays(alloc_flag, ierr)
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               rho_P_chiT_chiRho, nz, nz_alloc_extra, 'brunt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               rho_P_chiT_chiRho_face, nz, nz_alloc_extra, 'brunt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               work1, nwork*nz, nz_alloc_extra, 'brunt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               f1, 4*nz, nz_alloc_extra, 'brunt', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays

         subroutine smooth_brunt_B(work)
            use star_utils, only: threshold_smoothing
            real(dp), pointer, dimension(:) :: work
            logical, parameter :: preserve_sign = .false.
            if (s% num_cells_for_smooth_brunt_B <= 0) return
            call threshold_smoothing( &
               s% brunt_B, s% threshold_for_smooth_brunt_B, s% nz, s% num_cells_for_smooth_brunt_B, preserve_sign, work)
         end subroutine smooth_brunt_B


      end subroutine do_brunt_N2


      subroutine do_brunt_B_MHM_form(s, ierr)
         ! Brassard from Mike Montgomery (MHM)
         use star_utils, only: get_face_values
         use interp_1d_def

         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         real(dp), pointer, dimension(:) :: T_face, rho_face, chiT_face
         real(dp) :: brunt_B
         integer :: nz, species, k, i, op_err
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         nz = s% nz
         species = s% species
         
         nullify(T_face, rho_face, chiT_face)

         call do_alloc(ierr)
         if (ierr /= 0) then
            write(*,*) 'allocate failed in do_brunt_MHM_form'
            ierr = -1
            return
         end if

         call get_face_values(s, s% chiT, chiT_face, ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if

         call get_face_values(s, s% T, T_face, ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if

         call get_face_values(s, s% rho, rho_face, ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if

!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k=1,nz
            op_err = 0
            call get_brunt_B(&
               s, species, nz, k, T_face(k), rho_face(k), chiT_face(k), op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
         if (ierr /= 0) then
            call dealloc
            return
         end if

         call dealloc

         contains
            
         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            call do_work_arrays(.true.,ierr)
         end subroutine do_alloc

         subroutine dealloc
            integer :: ierr
            call do_work_arrays(.false.,ierr)
         end subroutine dealloc

         subroutine do_work_arrays(alloc_flag, ierr)
            use interp_1d_def
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               T_face, nz, nz_alloc_extra, 'brunt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               rho_face, nz, nz_alloc_extra, 'brunt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               chiT_face, nz, nz_alloc_extra, 'brunt', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays

      end subroutine do_brunt_B_MHM_form
      

      subroutine get_brunt_B(s, species, nz, k, T_face, rho_face, chiT_face, ierr)
         use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_lnPgas
         use eos_support, only: get_eos

         type (star_info), pointer :: s
         integer, intent(in) :: species, nz, k
         real(dp), intent(in) :: T_face, rho_face, chiT_face
         integer, intent(out) :: ierr

         real(dp) :: lnP1, lnP2, logRho_face, logT_face, Prad_face, &
            alfa, Ppoint, dlnP_dm, delta_lnP
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_eos_dlnd, d_eos_dlnT, d_eos_dabar, d_eos_dzbar
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
            write(*,*)
            write(*,2) 's% dq(k-1)', k-1, s% dq(k-1)
            write(*,2) 's% dq(k)', k, s% dq(k)
            write(*,*)
            write(*,2) 's% T(k-1)', k-1, s% T(k-1)
            write(*,2) 'T_face', k, T_face
            write(*,2) 's% T(k)', k, s% T(k)
            write(*,*)
            write(*,2) 's% rho(k-1)', k-1, s% rho(k-1)
            write(*,2) 'rho_face', k, rho_face
            write(*,2) 's% rho(k)', k, s% rho(k)
            write(*,*)
            write(*,2) 's% chiT(k-1)', k-1, s% chiT(k-1)
            write(*,2) 'chiT_face', k, chiT_face
            write(*,2) 's% chiT(k)', k, s% chiT(k)
            write(*,*)
            stop 'get_brunt_B'
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

         delta_lnP = s% lnP(k-1) - s% lnP(k)
         if (delta_lnP > -1d-50) then
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            Ppoint = alfa*s% P(k) + (1-alfa)*s% P(k-1)
            dlnP_dm = -s% cgrav(k)*s% m(k)/(4*pi*pow4(s% r(k))*Ppoint)
            delta_lnP = dlnP_dm*s% dm_bar(k)
         end if

         s% brunt_B(k) = (lnP1 - lnP2)/delta_lnP/chiT_face
         if (is_bad_num(s% brunt_B(k))) then
            ierr = -1
            s% retry_message = 'bad num for brunt_B'
            if (s% report_ierr) then
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               write(*,2) 'chiT_face', k, chiT_face
               write(*,2) 'delta_lnP', k, delta_lnP
               write(*,2) 's% lnP(k)', k, s% lnP(k)
               write(*,2) 's% lnP(k-1)', k-1, s% lnP(k-1)
               write(*,2) 'lnP1', k, lnP1
               write(*,2) 'lnP2', k, lnP2
               write(*,*)
               !stop 'do_brunt_B_MHM_form'
            end if
            if (s% stop_for_bad_nums) then
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               stop 'do_brunt_B_MHM_form'
            end if
         end if

      end subroutine get_brunt_B


      ! dlnRho_form

      !   A = (1/gamma1)*(dlnP/dlnR) - (dlnRho/dlnR)
      !   N2 = A*g/r
      !
      !   N2 = g^2 * (rho/P)*(chiT/chiRho) * (B - (gradT_sub_grada))
      !
      !   B = gradT_sub_grada + N2/(g^2 * (rho/P)*(chiT/chiRho))

      ! treat gamma1, lnP, and lnd as values at rmid
      ! interpolate to r at edge to get gamma1_edge, lnP_edge, lnd_edge
      ! using lnR as x, interpolate to get slopes dlnP/dlnR and dlnd/dlnR at edges
      ! combine those results to get A.  then N2 from A.
      ! interpolate in m to get (rho/P)*(chiT/chiRho) at edge.
      ! use that with N2 to get B for use with Ledoux

      subroutine do_brunt_B_dlnRho_form(s,rho_P_chiT_chiRho_face,ierr)

         type (star_info), pointer :: s
         real(dp), intent(in) :: rho_P_chiT_chiRho_face(:)
         integer, intent(out) :: ierr

         real(dp) :: alfa, beta, gamma1_face, P_face, dr, dlnd_dr, &
            dlnP_dr, brunt_N2
         integer :: nz, k, species
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         nz = s% nz
         species = s% species

         !   A = (1/gamma1)*(dlnP/dlnR) - (dlnRho/dlnR)
         !   N2 = A*g/r = g*((dlnP/dr)/gamma1 - (dlnRho/dr))
         !   N2 = g^2 * (rho/P)*(chiT/chiRho) * (B - (gradT_sub_grada))
         !   B = gradT_sub_grada + N2/(g^2 * (rho/P)*(chiT/chiRho))

         do k=2,nz
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1 - alfa
            gamma1_face = alfa*s% gamma1(k) + beta*s% gamma1(k-1)
            P_face = alfa*s% P(k) + beta*s% P(k-1)
            dr = s% rmid(k-1) - s% rmid(k)
            dlnd_dr = (s% lnd(k-1) - s% lnd(k))/dr
            dlnP_dr = (s% P(k-1) - s% P(k))/(P_face*dr)
            brunt_N2 = s% grav(k)*(dlnP_dr/gamma1_face - dlnd_dr)
            !s% profile_extra(k,1) = brunt_N2
            s% brunt_B(k) = s% gradT_sub_grada(k) + &
               brunt_N2/(s% grav(k)*s% grav(k)*rho_P_chiT_chiRho_face(k))
         end do
         s% brunt_B(1) = s% brunt_B(2)
         !s% profile_extra(1,1) = 0
         !s% profile_extra_name(1) = 'brunt_N2_dlnRho_form'

      end subroutine do_brunt_B_dlnRho_form


      subroutine do_brunt_B_gradmuX_form(s,ierr)
         use chem_def, only: ih1
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: nz, h1, k
         real(dp) :: A, B, beta_face, beta_factor, dlnmu_X, dlnP, x_face

         include 'formats'

         ierr = 0
         nz = s% nz
         h1 = s% net_iso(ih1)
         if (h1 <= 0) then
            s% brunt_B(1:nz) = 0
            return
         end if

         do k = 2, nz
            A = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            B = 1 - A
            beta_face = A*s% Pgas(k)/s% P(k) + B*s% Pgas(k-1)/s% P(k-1)
            beta_factor = beta_face/(4d0 - 3d0*beta_face)
            x_face = A*s% xa(h1,k) + B*s% xa(h1,k-1)
            dlnmu_X = -(s% xa(h1,k-1) - s% xa(h1,k))/(x_face + 0.6d0)
            dlnP = s% lnP(k-1) - s% lnP(k)
            if (dlnP > -1d-50) then
               !dlnP = -1d-50
               s% brunt_B(k) = 0
               cycle
            end if
            s% brunt_B(k) = beta_factor*dlnmu_X/dlnP
            if (s% brunt_B(k) < -1d10 .and. k == 1496) then
               write(*,2) 's% brunt_B(k)', k, s% brunt_B(k)
               write(*,2) 'beta_factor', k, beta_factor
               write(*,2) 'dlnmu_X', k, dlnmu_X
               write(*,2) 'dlnP', k, dlnP
               write(*,2) 's% xa(h1,k-1)', k-1, s% xa(h1,k-1)
               write(*,2) 's% xa(h1,k)', k, s% xa(h1,k)
               write(*,2) 's% lnP(k-1)', k-1, s% lnP(k-1)
               write(*,2) 's% lnP(k)', k, s% lnP(k)
               write(*,2) 's% lnP(k-1) - s% lnP(k)', k, s% lnP(k-1) - s% lnP(k)
               stop
            end if
         end do
         s% brunt_B(1) = s% brunt_B(2)

      end subroutine do_brunt_B_gradmuX_form


      end module brunt

