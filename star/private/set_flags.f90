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

      module set_flags

      use star_private_def
      use const_def, only: dp
      use utils_lib, only: is_bad
      use alloc

      implicit none

      contains

      subroutine set_v_flag(id, v_flag, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: v_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: nvar_hydro_old, k, nz, i_v, i_u
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% v_flag .eqv. v_flag) return

         nz = s% nz
         s% v_flag = v_flag
         nvar_hydro_old = s% nvar_hydro

         if (.not. v_flag) then  ! remove i_v's
            call del(s% xh)
            call del(s% xh_start)
            if (associated(s% xh_old) .and. s% generations > 1) call del(s% xh_old)
         end if

         call set_var_info(s, ierr)
         if (ierr /= 0) return

         call update_nvar_allocs(s, nvar_hydro_old, s% nvar_chem, ierr)
         if (ierr /= 0) return

         call check_sizes(s, ierr)
         if (ierr /= 0) return

         if (v_flag) then  ! insert i_v's
            i_v = s% i_v
            s% v_center = 0d0
            call insert(s% xh)
            call insert(s% xh_start)
            if (s% u_flag) then
               i_u = s% i_u
               do k=2,nz
                  s% xh(i_v,k) = 0.5d0*(s% xh(i_u,k-1) + s% xh(i_u,k))
               end do
               s% xh(i_v,1) = s% xh(i_u,1)
            else if (s% RSP_flag) then
               s% xh(i_v,1:nz) = 0d0
               s% v(1:nz) = 0d0
            else
               do k=1,nz
                  s% xh(i_v,k) = 0d0
                  if (is_bad(s% xh(i_v,k))) s% xh(i_v,k) = 0d0
                  s% v(k) = s% xh(i_v,k)
               end do
            end if
            if (associated(s% xh_old) .and. s% generations > 1) call insert(s% xh_old)
         end if

         call set_chem_names(s)

         if (v_flag .and. s% u_flag) then  ! turn off u_flag when turn on v_flag
            call set_u_flag(id, .false., ierr)
         end if

         contains

         subroutine del(xs)
            real(dp) :: xs(:,:)
            integer :: j, i_v
            if (size(xs,dim=2) < nz) return
            i_v = s% i_v
            do j = i_v+1, nvar_hydro_old
               xs(j-1,1:nz) = xs(j,1:nz)
            end do
         end subroutine del

         subroutine insert(xs)
            real(dp) :: xs(:,:)
            integer :: j, i_v
            if (size(xs,dim=2) < nz) return
            i_v = s% i_v
            do j = s% nvar_hydro, i_v+1, -1
               xs(j,1:nz) = xs(j-1,1:nz)
            end do
            xs(i_v,1:nz) = 0
         end subroutine insert

      end subroutine set_v_flag


      subroutine set_u_flag(id, u_flag, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: u_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: nvar_hydro_old, k, nz, i_u, i_v
         logical, parameter :: dbg = .false.

         integer :: num_u_vars

         include 'formats'

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% u_flag .eqv. u_flag) return

         nz = s% nz
         s% u_flag = u_flag
         nvar_hydro_old = s% nvar_hydro

         num_u_vars = 1

         if (.not. u_flag) then  ! remove
            call del(s% xh)
            call del(s% xh_start)
            if (associated(s% xh_old) .and. s% generations > 1) call del(s% xh_old)
         end if

         call set_var_info(s, ierr)
         if (ierr /= 0) return

         call update_nvar_allocs(s, nvar_hydro_old, s% nvar_chem, ierr)
         if (ierr /= 0) return

         call check_sizes(s, ierr)
         if (ierr /= 0) return

         if (u_flag) then  ! insert
            i_u = s% i_u
            call insert(s% xh)
            call insert(s% xh_start)
            if (s% v_flag) then  ! use v to initialize u
               i_v = s% i_v
               do k=1,nz-1
                  s% xh(i_u,k) = 0.5d0*(s% xh(i_v,k) + s% xh(i_v,k+1))
               end do
               k = nz
               s% xh(i_u,k) = 0.5d0*(s% xh(i_v,k) + s% v_center)
            else
               do k=1,nz
                  s% xh(i_u,k) = 0d0
               end do
            end if
            if (associated(s% xh_old) .and. s% generations > 1) call insert(s% xh_old)
            call fill_ad_with_zeros(s% u_face_ad,1,-1)
            call fill_ad_with_zeros(s% P_face_ad,1,-1)
         end if

         call set_chem_names(s)

         if (u_flag .and. s% v_flag) then  ! turn off v_flag when turn on u_flag
            call set_v_flag(id, .false., ierr)
         end if


         contains

         subroutine del(xs)
            real(dp) :: xs(:,:)
            integer :: k, j, i_u
            if (size(xs,dim=2) < nz) return
            i_u = s% i_u
            do k = 1, nz
               do j = i_u + num_u_vars, nvar_hydro_old
                  xs(j-num_u_vars,k) = xs(j,k)
               end do
            end do
         end subroutine del

         subroutine insert(xs)
            real(dp) :: xs(:,:)
            integer :: k, j, i_u
            if (size(xs,dim=2) < nz) return
            i_u = s% i_u
            do k = 1, nz
               do j = s% nvar_hydro, i_u + num_u_vars, -1
                  xs(j,k) = xs(j-num_u_vars,k)
               end do
               do j = i_u, i_u + num_u_vars - 1
                  xs(j,k) = 0
               end do
            end do

         end subroutine insert

      end subroutine set_u_flag


      subroutine set_RTI_flag(id, RTI_flag, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: RTI_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: nvar_hydro_old, nz
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% RTI_flag .eqv. RTI_flag) return

         nz = s% nz
         s% RTI_flag = RTI_flag
         nvar_hydro_old = s% nvar_hydro

         if (.not. RTI_flag) then  ! remove i_alpha_RTI's
            call del(s% xh)
            call del(s% xh_start)
            if (associated(s% xh_old) .and. s% generations > 1) call del(s% xh_old)
         end if

         call set_var_info(s, ierr)
         if (ierr /= 0) return

         call update_nvar_allocs(s, nvar_hydro_old, s% nvar_chem, ierr)
         if (ierr /= 0) return

         call check_sizes(s, ierr)
         if (ierr /= 0) return

         if (RTI_flag) then  ! insert i_alpha_RTI's
            call insert(s% xh)
            call insert(s% xh_start)
            s% xh(s% i_alpha_RTI,1:nz) = 0d0
            if (associated(s% xh_old) .and. s% generations > 1) call insert(s% xh_old)
         end if

         call set_chem_names(s)

         contains

         subroutine del(xs)
            real(dp) :: xs(:,:)
            integer :: j, i_alpha_RTI
            if (size(xs,dim=2) < nz) return
            i_alpha_RTI = s% i_alpha_RTI
            do j = i_alpha_RTI+1, nvar_hydro_old
               xs(j-1,1:nz) = xs(j,1:nz)
            end do
         end subroutine del

         subroutine insert(xs)
            real(dp) :: xs(:,:)
            integer :: j, i_alpha_RTI
            if (size(xs,dim=2) < nz) return
            i_alpha_RTI = s% i_alpha_RTI
            do j = s% nvar_hydro, i_alpha_RTI+1, -1
               xs(j,1:nz) = xs(j-1,1:nz)
            end do
            xs(i_alpha_RTI,1:nz) = 0
         end subroutine insert

      end subroutine set_RTI_flag


      subroutine set_pulsation_envelope_mesh(id, ierr)
         use tdc_hydro_support, only: remesh_for_TDC_pulsations
         use hydro_vars, only: set_vars
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% RSP2_flag) then
            write(*,*) 'doing automatic remesh for RSP2'
         else
            write(*,*) 'doing automatic remesh for TDC pulsations'
         end if
         call remesh_for_TDC_pulsations(s,ierr)
         if (ierr /= 0) return
         call set_vars(s, s% dt, ierr)
         if (ierr /= 0) return

      end subroutine set_pulsation_envelope_mesh

      subroutine set_RSP2_flag(id, RSP2_flag, ierr)
         use const_def, only: sqrt_2_div_3
         use hydro_vars, only: set_vars
         use auto_diff_support, only: rsp2_zero_w
         integer, intent(in) :: id
         logical, intent(in) :: RSP2_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp), allocatable :: gradT_old(:)
         real(dp) :: energy, dm_face
         integer :: nvar_hydro_old, i, k, j, first, last, nz
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !write(*,*) 'set_RSP2_flag previous s% RSP2_flag', s% RSP2_flag
         !write(*,*) 'set_RSP2_flag new RSP2_flag', RSP2_flag
         if (s% RSP2_flag .eqv. RSP2_flag) then
            call set_RSP2_3equation_flag(id, RSP2_flag .and. s% RSP2_use_3equation_model, ierr)
            return
         end if
         if (.not. RSP2_flag .and. s% RSP2_3equation_flag) then
            call set_RSP2_3equation_flag(id, .false., ierr)
            if (ierr /= 0) return
         end if

         nz = s% nz
         if (RSP2_flag) gradT_old = s% gradT(1:nz)

         s% RSP2_flag = RSP2_flag
         nvar_hydro_old = s% nvar_hydro

         if (.not. RSP2_flag) then
            call remove1(s% i_Y)
            call remove1(s% i_w)
         end if

         call set_var_info(s, ierr)
         if (ierr /= 0) return

         write(*,*) 'set_RSP2 variables and equations'
         if (.false.) then
            do i=1,s% nvar_hydro
               write(*,'(i3,2a20)') i, trim(s% nameofequ(i)), trim(s% nameofvar(i))
            end do
         end if

         call update_nvar_allocs(s, nvar_hydro_old, s% nvar_chem, ierr)
         if (ierr /= 0) return

         call check_sizes(s, ierr)
         if (ierr /= 0) return

         if (RSP2_flag) then
            call insert1(s% i_w)
            if (s% RSP_flag) then
               ! Deposit each RSP cell energy into its adjoining face volumes.
               first = max(2,s% RSP2_num_outermost_cells_forced_nonturbulent+2)
               last = nz-int(nz/s% RSP2_nz_div_IBOTOM)
               if (s% mixing_length_alpha == 0d0) last = 0
               s% xh(s% i_w,1:nz) = 0d0
               if (first <= last) then
                  do k=1,nz
                     energy = 0.5d0*s% dm(k)*max(0d0,s% xh(s% i_Et_RSP,k))
                     do i=k,k+1
                        j = min(last,max(first,i))
                        s% xh(s% i_w,j) = s% xh(s% i_w,j) + energy
                     end do
                  end do
                  do k=first,last
                     dm_face = 0.5d0*(s% dm(k-1) + s% dm(k))
                     s% xh(s% i_w,k) = sqrt(s% xh(s% i_w,k)/dm_face)
                  end do
               else if (any(s% xh(s% i_Et_RSP,1:nz) > 0d0)) then
                  write(*,*) 'no active face for RSP turbulent energy'
                  ierr = -1
                  return
               end if
            else if (s% have_mlt_vc) then
               do k=1,nz
                  s% xh(s% i_w,k) = s% mlt_vc(k)/sqrt_2_div_3
                  if (rsp2_zero_w(s,k)) s% xh(s% i_w,k) = 0d0
               end do
            else
               write(*,*) 'set_rsp2_flag true requires mlt_vc'
               ierr = -1
               return
            end if
            call insert1(s% i_Y)
         end if

         call set_chem_names(s)

         if (.not. RSP2_flag) return

         if (s% RSP_flag) then  ! turn off RSP_flag when turn on RSP2_flag
            call set_RSP_flag(id, .false., ierr)
            if (ierr /= 0) return
         end if

         if (.not. s% u_flag) then
            call set_v_flag(s% id, .true., ierr)
            if (ierr /= 0) return
         end if

         call set_vars(s, s% dt, ierr)
         if (ierr /= 0) return

         ! Use the target neutral gradient, which is not set by the RSP model.
         s% xh(s% i_Y,1:nz) = gradT_old - s% gradL(1:nz)
         s% xh(s% i_Y,1) = 0d0
         s% xh_start(s% i_Y,1:nz) = s% xh(s% i_Y,1:nz)
         call set_vars(s, s% dt, ierr)
         if (ierr /= 0) return

         if (s% RSP2_remesh_when_load .and. s% nz /= s% TDC_hydro_nz) then
            call set_pulsation_envelope_mesh(id, ierr)
            if (ierr /= 0) return
         end if

         call set_RSP2_3equation_flag(id, s% RSP2_use_3equation_model, ierr)

         contains

         subroutine insert1(i_var)
            integer, intent(in) :: i_var
            include 'formats'
            call insert(s% xh,i_var)
            call insert(s% xh_start,i_var)
            do k=1,nz
               s% xh(i_var,k) = 0d0
            end do
            if (associated(s% xh_old) .and. s% generations > 1) then
               call insert(s% xh_old,i_var)
            end if
         end subroutine insert1

         subroutine remove1(i_remove)
            integer, intent(in) :: i_remove
            call del(s% xh,i_remove)
            call del(s% xh_start,i_remove)
            if (associated(s% xh_old) .and. s% generations > 1) then
               call del(s% xh_old,i_remove)
            end if
         end subroutine remove1

         subroutine del(xs,i_var)
            real(dp) :: xs(:,:)
            integer, intent(in) :: i_var
            integer :: j, k
            if (size(xs,dim=2) < nz) return
            do j = i_var+1, nvar_hydro_old
               do k=1,nz
                  xs(j-1,k) = xs(j,k)
               end do
            end do
         end subroutine del

         subroutine insert(xs,i_var)
            real(dp) :: xs(:,:)
            integer, intent(in) :: i_var
            integer :: j, k
            if (size(xs,dim=2) < nz) return
            do j = s% nvar_hydro, i_var+1, -1
               do k=1,nz
                  xs(j,k) = xs(j-1,k)
               end do
            end do
            xs(i_var,1:nz) = 0d0
         end subroutine insert

      end subroutine set_RSP2_flag


      subroutine set_RSP2_3equation_flag(id, enabled, ierr)
         use hydro_vars, only: set_vars, unpack_xh
         use hydro_rsp2, only: init_rsp2_moments
         integer, intent(in) :: id
         logical, intent(in) :: enabled
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         real(dp), allocatable :: xh_save(:,:), xh_old_save(:,:), Lc_old(:), gradT_old(:)
         integer :: nvar_old, i_Y, i_Pi, i_Phi, j, jj
         logical :: have_old

         call get_star_ptr(id,s,ierr)
         if (ierr /= 0) return
         if (enabled) then
            if (.not. s% RSP2_flag .or. is_bad(s% RSP2_alfa_pi) .or. is_bad(s% RSP2_alfa_phi) .or. &
                  s% RSP2_alfa_pi <= 0d0 .or. &
                  s% RSP2_alfa_phi <= 0d0 .or. s% RSP2_source_seed /= 0d0) then
               write(*,*) 'RSP2 three equation model requires RSP2, positive alfa coefficients and zero source seed'
               ierr = -1
               return
            end if
         end if
         if (s% RSP2_3equation_flag .eqv. enabled) return
         call set_vars(s,s% dt,ierr)
         if (ierr /= 0) return
         nvar_old = s% nvar_hydro
         i_Y = s% i_Y
         i_Pi = s% i_Pi
         i_Phi = s% i_Phi
         xh_save = s% xh(1:nvar_old,1:s% nz)
         Lc_old = s% Lc(1:s% nz)
         gradT_old = s% gradT(1:s% nz)
         have_old = associated(s% xh_old) .and. s% generations > 1
         if (have_old) xh_old_save = s% xh_old(1:nvar_old,1:s% nz_old)

         s% RSP2_3equation_flag = enabled
         call set_var_info(s,ierr)
         if (ierr /= 0) return
         call update_nvar_allocs(s,nvar_old,s% nvar_chem,ierr)
         if (ierr /= 0) return
         s% xh(1:s% nvar_hydro,1:s% nz) = 0d0
         if (have_old) s% xh_old(1:s% nvar_hydro,1:s% nz_old) = 0d0
         do j=1,nvar_old
            jj = j
            if (enabled) then
               if (j > i_Y) jj = j+2
            else
               if (j == i_Pi .or. j == i_Phi) cycle
               if (j > i_Phi) jj = j-2
            end if
            s% xh(jj,1:s% nz) = xh_save(j,1:s% nz)
            if (have_old) s% xh_old(jj,1:s% nz_old) = xh_old_save(j,1:s% nz_old)
         end do
         call set_chem_names(s)
         call set_vars(s,s% dt,ierr)
         if (ierr /= 0) return
         ! Initialize moments using the preserved gradient and its new reference.
         s% xh(s% i_Y,1:s% nz) = gradT_old - s% gradL(1:s% nz)
         s% xh(s% i_Y,1) = 0d0
         call unpack_xh(s,ierr)
         if (ierr /= 0) return
         if (enabled) then
            call init_rsp2_moments(s,Lc_old,ierr)
            if (ierr /= 0) return
         end if
         s% xh_start(1:s% nvar_hydro,1:s% nz) = s% xh(1:s% nvar_hydro,1:s% nz)
         ! Missing moment history cannot be used for temporal extrapolation.
         s% generations = 1
         s% need_to_setvars = .true.
         call set_vars(s,s% dt,ierr)
      end subroutine set_RSP2_3equation_flag


      subroutine set_RSP_flag(id, RSP_flag, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: RSP_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: nvar_hydro_old, k, nz
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% RSP_flag .eqv. RSP_flag) return

         nz = s% nz
         s% RSP_flag = RSP_flag
         nvar_hydro_old = s% nvar_hydro

         if (.not. RSP_flag) then
            call remove1(s% i_Fr_RSP)
            call remove1(s% i_erad_RSP)
            call remove1(s% i_Et_RSP)
         else if (s% i_lum /= 0) then
            call remove1(s% i_lum)
         end if

         call set_var_info(s, ierr)
         if (ierr /= 0) return

         call update_nvar_allocs(s, nvar_hydro_old, s% nvar_chem, ierr)
         if (ierr /= 0) return

         call check_sizes(s, ierr)
         if (ierr /= 0) return

         if (RSP_flag) then
            call insert1(s% i_Et_RSP)
            call insert1(s% i_erad_RSP)
            call insert1(s% i_Fr_RSP)
         else
            call insert1(s% i_lum)
            do k=1,nz
               s% xh(s% i_lum,k) = s% L(k)
            end do
         end if

         call set_chem_names(s)

         if (RSP_flag) call set_v_flag(s% id, .true., ierr)

         contains

         subroutine insert1(i_var)
            integer, intent(in) :: i_var
            call insert(s% xh,i_var)
            call insert(s% xh_start,i_var)
            do k=1,nz
               s% xh(i_var,k) = 0d0
            end do
            if (associated(s% xh_old) .and. s% generations > 1) then
               call insert(s% xh_old,i_var)
            end if
         end subroutine insert1

         subroutine remove1(i_remove)
            integer, intent(in) :: i_remove
            call del(s% xh,i_remove)
            call del(s% xh_start,i_remove)
            if (associated(s% xh_old) .and. s% generations > 1) then
               call del(s% xh_old,i_remove)
            end if
         end subroutine remove1

         subroutine del(xs,i_var)
            real(dp) :: xs(:,:)
            integer, intent(in) :: i_var
            integer :: j, k
            if (size(xs,dim=2) < nz) return
            do j = i_var+1, nvar_hydro_old
               do k=1,nz
                  xs(j-1,k) = xs(j,k)
               end do
            end do
         end subroutine del

         subroutine insert(xs,i_var)
            real(dp) :: xs(:,:)
            integer, intent(in) :: i_var
            integer :: j, k
            if (size(xs,dim=2) < nz) return
            do j = s% nvar_hydro, i_var+1, -1
               do k=1,nz
                  xs(j,k) = xs(j-1,k)
               end do
            end do
            xs(i_var,1:nz) = 0
         end subroutine insert

      end subroutine set_RSP_flag


      subroutine set_w_div_wc_flag(id, w_div_wc_flag, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: w_div_wc_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: nvar_hydro_old, nz
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% w_div_wc_flag .eqv. w_div_wc_flag) return

         nz = s% nz
         s% w_div_wc_flag = w_div_wc_flag
         nvar_hydro_old = s% nvar_hydro

         if (.not. w_div_wc_flag) then  ! remove i_w_div_wc's
            call del(s% xh)
            call del(s% xh_start)
            if (associated(s% xh_old) .and. s% generations > 1) call del(s% xh_old)
         end if

         call set_var_info(s, ierr)
         if (ierr /= 0) return

         call update_nvar_allocs(s, nvar_hydro_old, s% nvar_chem, ierr)
         if (ierr /= 0) return

         call check_sizes(s, ierr)
         if (ierr /= 0) return

         if (w_div_wc_flag) then  ! insert i_w_div_w's
            call insert(s% xh)
            call insert(s% xh_start)
            s% xh(s% i_w_div_wc,1:nz) = 0d0
            if (associated(s% xh_old) .and. s% generations > 1) call insert(s% xh_old)
         end if

         call set_chem_names(s)

         contains

         subroutine del(xs)
            real(dp) :: xs(:,:)
            integer :: j, i_w_div_wc
            if (size(xs,dim=2) < nz) return
            i_w_div_wc = s% i_w_div_wc
            do j = i_w_div_wc+1, nvar_hydro_old
               xs(j-1,1:nz) = xs(j,1:nz)
            end do
         end subroutine del

         subroutine insert(xs)
            real(dp) :: xs(:,:)
            integer :: j, i_w_div_wc
            if (size(xs,dim=2) < nz) return
            i_w_div_wc = s% i_w_div_wc
            do j = s% nvar_hydro, i_w_div_wc+1, -1
               xs(j,1:nz) = xs(j-1,1:nz)
            end do
            xs(i_w_div_wc,1:nz) = 0
         end subroutine insert

      end subroutine set_w_div_wc_flag


      subroutine set_j_rot_flag(id, j_rot_flag, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: j_rot_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: nvar_hydro_old, nz
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% j_rot_flag .eqv. j_rot_flag) return

         nz = s% nz
         s% j_rot_flag = j_rot_flag
         nvar_hydro_old = s% nvar_hydro

         if (.not. j_rot_flag) then  ! remove i_j_rot's
            call del(s% xh)
            call del(s% xh_start)
            if (associated(s% xh_old) .and. s% generations > 1) call del(s% xh_old)
         end if

         call set_var_info(s, ierr)
         if (ierr /= 0) return

         call update_nvar_allocs(s, nvar_hydro_old, s% nvar_chem, ierr)
         if (ierr /= 0) return

         call check_sizes(s, ierr)
         if (ierr /= 0) return

         if (j_rot_flag) then  ! insert i_j_rot's
            call insert(s% xh)
            call insert(s% xh_start)
            s% xh(s% i_j_rot,1:nz) = 0d0
            if (associated(s% xh_old) .and. s% generations > 1) call insert(s% xh_old)
         end if

         call set_chem_names(s)

         contains

         subroutine del(xs)
            real(dp) :: xs(:,:)
            integer :: j, i_j_rot
            if (size(xs,dim=2) < nz) return
            i_j_rot = s% i_j_rot
            do j = i_j_rot+1, nvar_hydro_old
               xs(j-1,1:nz) = xs(j,1:nz)
            end do
         end subroutine del

         subroutine insert(xs)
            real(dp) :: xs(:,:)
            integer :: j, i_j_rot
            if (size(xs,dim=2) < nz) return
            i_j_rot = s% i_j_rot
            do j = s% nvar_hydro, i_j_rot+1, -1
               xs(j,1:nz) = xs(j-1,1:nz)
            end do
            xs(i_j_rot,1:nz) = 0
         end subroutine insert

      end subroutine set_j_rot_flag


      subroutine set_D_omega_flag(id, D_omega_flag, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: D_omega_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% D_omega_flag .eqv. D_omega_flag) return
         s% D_omega_flag = D_omega_flag
         s% D_omega(1:s% nz) = 0
      end subroutine set_D_omega_flag


      subroutine set_am_nu_rot_flag(id, am_nu_rot_flag, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: am_nu_rot_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% am_nu_rot_flag .eqv. am_nu_rot_flag) return
         s% am_nu_rot_flag = am_nu_rot_flag
         s% am_nu_rot(1:s% nz) = 0
      end subroutine set_am_nu_rot_flag


      subroutine set_rotation_flag(id, rotation_flag, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: rotation_flag
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         include 'formats'

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% rotation_flag .eqv. rotation_flag) return

         s% rotation_flag = rotation_flag
         s% omega(1:s% nz) = 0
         s% j_rot(1:s% nz) = 0
         s% D_omega(1:s% nz) = 0
         s% am_nu_rot(1:s% nz) = 0

         if (.not. rotation_flag) then
            call set_w_div_wc_flag(id, .false., ierr)
            if (ierr /= 0) return
            call set_j_rot_flag(id, .false., ierr)
            if (ierr /= 0) return
            return
         end if

         if (s% job% use_w_div_wc_flag_with_rotation) then
            call set_w_div_wc_flag(id, .true., ierr)
            if (ierr /= 0) return
            if (s% job% use_j_rot_flag_with_rotation) then
               call set_j_rot_flag(id, .true., ierr)
               if (ierr /= 0) return
            end if
         end if

         call zero_array(s% nu_ST)
         call zero_array(s% D_ST)
         call zero_array(s% D_DSI)
         call zero_array(s% D_SH)
         call zero_array(s% D_SSI)
         call zero_array(s% D_ES)
         call zero_array(s% D_GSF)

         call zero_array(s% prev_mesh_omega)
         call zero_array(s% prev_mesh_j_rot)


         contains

         subroutine zero_array(d)
            real(dp), pointer :: d(:)
            if (.not. associated(d)) return
            d(:) = 0
         end subroutine zero_array

      end subroutine set_rotation_flag

      end module set_flags
