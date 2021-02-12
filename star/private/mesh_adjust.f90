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


      module mesh_adjust

      use const_def
      use star_private_def
      use chem_def
      use interp_1d_def, only: pm_work_size
      use utils_lib

      implicit none

      private
      public :: do_mesh_adjust, do_prune_mesh_surface, &
         set_lnT_for_energy_with_tol, set_lnT_for_energy

      integer, parameter :: nwork = pm_work_size

      real(dp), parameter :: eta_limit = -1d-6


      logical, parameter :: dbg = .false.


      contains


      subroutine do_mesh_adjust( &
            s, nz, nz_old, xh_old, xa_old, &
            energy_old, eta_old, lnd_old, lnPgas_old, &
            j_rot_old, i_rot_old, omega_old, D_omega_old, &
            conv_vel_old, lnT_old, et_old, specific_PE_old, specific_KE_old, &
            old_m, old_r, old_rho, dPdr_dRhodr_info_old, D_mix_old, &
            cell_type, comes_from, dq_old, xq_old, xh, xa, dq, xq, ierr)
         use interp_1d_def
         use interp_1d_lib
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old
         integer, dimension(:) :: cell_type, comes_from
         real(dp), dimension(:), pointer :: &
            dq_old, xq_old, dq, xq, energy_old, eta_old, &
            lnd_old, lnPgas_old, conv_vel_old, lnT_old, et_old, &
            specific_PE_old, specific_KE_old, &
            old_m, old_r, old_rho, dPdr_dRhodr_info_old, &
            j_rot_old, i_rot_old, omega_old, D_omega_old, D_mix_old
         real(dp), dimension(:,:), pointer :: xh_old, xa_old
         real(dp), dimension(:,:), pointer :: xh, xa
         integer, intent(out) :: ierr

         real(dp) :: dxa, xmstar, mstar, sumx, remove1, remove2, &
            total_internal_energy1, total_internal_energy2, err
         character (len=strlen) :: message
         integer :: k, from_k, j, op_err, nzlo, nzhi, nzlo_old, nzhi_old, species
         logical :: found_bad_one
         real(dp), pointer :: work(:)
         real(dp), dimension(:), allocatable :: &
            dqbar, dqbar_old, new_r, Vol_new, xq_old_plus1, &
            xout_old, xout_new, xq_new, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, &
            energy_new, density_new
         real(dp), dimension(:,:), allocatable :: xa_c0, xa_c1, xa_c2

         include 'formats'

         ierr = 0
         species = s% species

         xmstar = s% xmstar
         mstar = xmstar + s% M_center

         ! check xq's
         do k=1,nz
            if (xq(k) < 0 .or. xq(k) > 1) then
               ierr = -1
               return

               write(*,*) 'k', k
               write(*,*) 'xq(k)', xq(k)
               stop 'debug: do_mesh_adjust'
            end if
         end do

         if (dbg) write(*,*) 'enter do_mesh_adjust'

         nzlo = 0
         do k = 1, nz
            if (cell_type(k) /= unchanged_type) then
               if (dbg) write(*,2) 'nzlo changed', k
               nzlo = k; exit
            end if
         end do
         if (nzlo == 0) then
            if (dbg) write(*,2) 'no cells changed'
            nzlo = nz
         end if

         nzhi = nzlo
         do k = nz, nzlo, -1
            if (cell_type(k) /= unchanged_type) then
               if (dbg) write(*,2) 'nzhi changed', k
               nzhi = k; exit
            end if
         end do

         ! extend range for purposes of interpolation
         if (nzhi < nz) nzhi = nzhi+1
         if (nzlo > 1) nzlo = nzlo-1

         nzlo_old = comes_from(nzlo)
         if (nzhi == nz) then
            nzhi_old = nz_old
         else
            nzhi_old = comes_from(nzhi+1)
         end if

         call do_alloc(ierr)
         if (ierr /= 0) return

         do k=2,nz-1
            dqbar(k) = 0.5d0*(dq(k-1) + dq(k))
         end do
         dqbar(1) = 0.5d0*dq(1)
         dqbar(nz) = 0.5d0*dq(nz-1) + dq(nz)

         do k=2,nz_old-1
            dqbar_old(k) = 0.5d0*(dq_old(k-1) + dq_old(k))
         end do
         dqbar_old(1) = 0.5d0*dq_old(1)
         dqbar_old(nz_old) = 0.5d0*dq_old(nz_old-1) + dq_old(nz_old)

         do k=1,nz_old
            xq_old_plus1(k) = xq_old(k)
         end do
         ! add point at true center so can interpolate xq_new > xq_old(nz_old)
         xq_old_plus1(nz_old+1) = 1

         do k = 1, nzhi - nzlo + 1
            xq_new(k) = xq(nzlo+k-1)
         end do

         xout_old(1) = xq_old(1)
         do k=2,nz_old
            xout_old(k) = xout_old(k-1) + dqbar_old(k-1)
         end do

         xout_new(1) = xq(1)
         do k=2,nz
            xout_new(k) = xout_new(k-1) + dqbar(k-1)
         end do

         if (dbg) write(*,*) 'call do_L'
         call do_L( &
            s, nz, nz_old, nzlo, nzhi, comes_from, &
            xh, xh_old, xq, xq_old_plus1, xq_new, &
            work, tmp1, tmp2, ierr)
         if (failed('do_L')) return

         if (s% RTI_flag) then
            if (dbg) write(*,*) 'call do_alpha_RTI'
            call do_alpha_RTI( &
               s, nz, nz_old, nzlo, nzhi, comes_from, &
               xh, xh_old, xq, xq_old_plus1, xq_new, &
               work, tmp1, tmp2, ierr)
            if (failed('do_alpha_RTI')) return
            call do_interp_pt_val( &
               s, nz, nz_old, nzlo, nzhi, s% dPdr_dRhodr_info, dPdr_dRhodr_info_old, &
               0d0, xq, xq_old_plus1, xq_new, .true., work, tmp1, tmp2, ierr)
            if (failed('dPdr_dRhodr_info')) return
         end if

         if (s% et_flag) then
            if (dbg) write(*,*) 'call do_et'
            call do_et( &
               s, nz, nz_old, nzlo, nzhi, comes_from, &
               xh, xh_old, xq, xq_old_plus1, xq_new, &
               work, tmp1, tmp2, ierr)
            if (failed('do_et')) return
         end if

         if (dbg) write(*,*) 'call do_lnR_and_lnd'
         call do_lnR_and_lnd( &
            s, nz, nz_old, nzlo, nzhi, cell_type, comes_from, &
            xh, xh_old, xmstar, lnd_old, lnPgas_old, &
            dqbar, dqbar_old, old_r, old_m, old_rho, &
            dq, dq_old, xq, xq_old_plus1, density_new, work, &
            tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, ierr)
         if (failed('do_lnR_and_lnd')) return

         if (s% v_flag) then ! calculate new v to conserve kinetic energy
            if (dbg) write(*,*) 'call do_v'
            call do_v( &
               s, nz, nz_old, cell_type, comes_from, &
               xq_old, xq, dq_old, dq, xh, xh_old, &
               xout_old, xout_new, dqbar_old, dqbar, tmp1, ierr)
            if (failed('do_v')) return
         end if

         if (s% u_flag) then ! calculate new u to conserve kinetic energy
            if (dbg) write(*,*) 'call do_u'
            call do_u( &
               s, nz, nz_old, cell_type, comes_from, &
               xq_old, xq, dq_old, dq, xh, xh_old, &
               xout_old, xout_new, tmp1, ierr)
            if (failed('do_u')) return
         end if

         if (s% conv_vel_flag) then
            if (dbg) write(*,*) 'call do_conv_vel'
            call do_conv_vel( &
               s, nz, nz_old, nzlo, nzhi, comes_from, &
               xh, xh_old, xq, xq_old_plus1, xq_new, &
               work, tmp1, tmp2, ierr)
            if (failed('do_conv_vel')) return
         end if

         if (s% rotation_flag) then
            call adjust_omega(s, nz, nz_old, comes_from, &
               xq_old, xq, dq_old, dq, xh, j_rot_old, &
               xout_old, xout_new, dqbar_old, dqbar, ierr)
            if (failed('adjust_omega')) return
            if (s% D_omega_flag) then
               call do_interp_pt_val( &
                  s, nz, nz_old, nzlo, nzhi, s% D_omega, D_omega_old, &
                  0d0, xq, xq_old_plus1, xq_new, .true., work, tmp1, tmp2, ierr)
               if (failed('D_omega')) return
            end if
         end if
         
         call do_interp_pt_val( &
            s, nz, nz_old, nzlo, nzhi, s% D_mix, D_mix_old, &
            0d0, xq, xq_old_plus1, xq_new, .true., work, tmp1, tmp2, ierr)
         if (failed('D_mix')) return

         do k=nzlo_old,nzhi_old  ! 1,nz_old  !

            ! since we must adjust things to make the sum of xa's = 1,
            ! only do linear reconstruction.
            do j=1,species
               call get1_lpp(k, species, nz_old, j, dq_old, xa_old, &
                              .false., xa_c0(:,j), xa_c1(:,j), xa_c2(:,j))
            end do

            sumx = sum(xa_old(1:species,k))
            do j=1,species
               xa_c0(k,j) = xa_old(j,k)/sumx ! make sure that adds to 1
               xa_c2(k,j) = 0 ! no curvature terms
            end do

            ! only reduce magnitude of slopes
            ! so don't risk producing values out of [0..1] range
            if (sum(xa_c1(k,:)) > 0) then
               j = maxloc(xa_c1(k,:), dim=1)
            else
               j = minloc(xa_c1(k,:), dim=1)
            end if
            xa_c1(k,j) = 0
            xa_c1(k,j) = -sum(xa_c1(k,:))
            ! check for valid fractions at boundaries; set slopes to 0 if find a bad one.
            do j=1,species
               dxa = abs(xa_c1(k,j))*dq_old(k)/2
               if (xa_c0(k,j) + dxa > 1 .or. xa_c0(k,j) - dxa < 0) then
                  xa_c1(k,:) = 0
                  exit
               end if
            end do

         end do

         if (failed('adjust_mesh nz_old parallel loop')) return

         if (dbg) write(*,*) 'do xa and lnT'

         total_internal_energy1 = &
            dot_product(dq_old(1:nz_old), energy_old(1:nz_old))

         do k = 1, nz

            op_err = 0

            ! calculate new abundances to conserve species
            call do_xa( &
               s, nz, nz_old, k, species, cell_type, comes_from, xa, xa_old, &
               xa_c0, xa_c1, xa_c2, xq, dq, xq_old, dq_old, &
               s% mesh_adjust_use_quadratic, op_err)
            if (op_err /= 0) then
               write(*,2) 'failed for do_xa', k
               stop
               write(message,*) 'do_xa for k', k
               ierr = op_err
            end if

            ! calculate new temperatures to conserve energy
            call do1_lnT( &
               s, nz_old, k, species, cell_type, comes_from, &
               xa, xh, xh_old, &
               xq, dq, xq_old, dq_old, eta_old, energy_old, lnT_old, &
               specific_PE_old, specific_KE_old, et_old, &
               density_new, energy_new, op_err)
            if (op_err /= 0) then
               write(*,2) 'failed for do1_lnT', k
               stop
               write(message,*) 'do1_lnT for k', k
               ierr = op_err
            end if
            if (is_bad(energy_new(k)) .or. is_bad(dq(k))) then
               write(*,2) 'energy_new', k, energy_new(k)
               write(*,2) 'dq', k, dq(k)
               stop ''
            end if

         end do

         if (failed(message)) return

         total_internal_energy2 = dot_product(dq(1:nz), energy_new(1:nz))
         err = abs(total_internal_energy1 - total_internal_energy2)/ &
               max(abs(total_internal_energy1),abs(total_internal_energy2),1d0)
         s% mesh_adjust_IE_conservation = err

         if (s% trace_mesh_adjust_error_in_conservation) then

            write(*,2) 'mesh adjust error in conservation of IE', &
               s% model_number, err, total_internal_energy2, total_internal_energy1
            if (err > 1d-8) then
               call show_errors
               write(*,*) 'err too large'
               stop 'mesh adjust'
            end if

         end if

         if (dbg) write(*,*) 'call check_species_conservation'
         call check_species_conservation(species,ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in check_species_conservation'
            stop
         end if

         call dealloc

         contains

         subroutine show_errors
            integer :: k0, k0_from, k, k_from, k_outer
            real(dp) :: new_sum, old_sum

            include 'formats'

            k0 = 1
            k0_from = 1
            k_outer = 2
            write(*,*)
            do
               if (cell_type(k_outer) == unchanged_type .and. k_outer /= nz) then
                  k_outer = k_outer + 1
                  cycle
               end if
               k0 = k_outer
               k0_from = comes_from(k_outer)
               do k = k0+1, nz
                  if (cell_type(k) /= unchanged_type .and. k /= nz) cycle
                  new_sum = dot_product(dq(k0:k),energy_new(k0:k))
                  k_from = comes_from(k)
                  old_sum = &
                     dot_product(dq_old(k0_from:k_from),energy_old(k0_from:k_from))
                  write(*,5) 'section err', k0, k, k0_from, k_from, &
                     abs(old_sum - new_sum)/max(abs(old_sum),abs(new_sum),1d0), &
                     old_sum, new_sum
                  k_outer = k
                  exit
               end do
               if (k_outer == nz) exit
               k_outer = k_outer + 1
            end do
            write(*,2) 'nz_old', nz_old
            write(*,2) 'nz', nz
            write(*,*)

         end subroutine show_errors

         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            integer :: sz
            sz = max(nz, nz_old) + 1
            call do_work_arrays(.true.,ierr)
            allocate( &
               dqbar(sz), dqbar_old(sz), new_r(sz), Vol_new(sz), xq_old_plus1(sz), &
               xout_old(sz), xout_new(sz), xq_new(sz), energy_new(sz), density_new(sz), &
               tmp1(sz), tmp2(sz), tmp3(sz), tmp4(sz), tmp5(sz), tmp6(sz), tmp7(sz), &
               xa_c0(sz,species), xa_c1(sz,species), xa_c2(sz,species))            
         end subroutine do_alloc

         subroutine dealloc
            call do_work_arrays(.false.,ierr)
         end subroutine dealloc

         subroutine do_work_arrays(alloc_flag, ierr)
            use interp_1d_def
            use alloc
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
                work, (nz_old+1)*pm_work_size, nz_alloc_extra, 'mesh_adjust', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays


         logical function failed(msg)
            character (len=*) :: msg
            if (ierr == 0) then
               failed = .false.
               return
            end if
            failed = .true.
            if (dbg) write(*, *) 'mesh_revisions failed in ' // trim(msg)
            call dealloc
            return
         end function failed


         subroutine check_species_conservation(species,ierr)
            integer, intent(in) :: species
            integer, intent(out) :: ierr
            integer :: j, k, jbad
            real(dp) :: old_total, new_total
            logical :: okay
            include 'formats'
            ierr = 0
            okay = .true.
            jbad = -1
            do j=1,species
               old_total = dot_product(xa_old(j,1:nz_old),dq_old(1:nz_old))
               if (old_total < 1d-9) cycle
               new_total = dot_product(xa(j,1:nz),dq(1:nz))
               if (abs(new_total - old_total) > 1d-4) then ! check for major problems
                  ierr = -1
                  jbad = j
                  okay = .false.
                  if (dbg) then
                     write(*,*) 'problem with conservation of species ' //  &
                        chem_isos% name(s% chem_id(j))
                     write(*,1) 'new mass fraction', new_total
                     write(*,1) 'old mass fraction', old_total
                     write(*,1) 'new - old', new_total - old_total
                     write(*,1) '(new - old)/old', (new_total - old_total) / old_total
                     write(*,*)
                  end if
               end if
            end do
            if (okay) return
            ierr = -1
            write(*,*)
            do j=1,species
               old_total = dot_product(xa_old(j,1:nz_old),dq_old(1:nz_old))
               if (old_total < 1d-9) cycle
               new_total = dot_product(xa(j,1:nz),dq(1:nz))
               write(*,2) 'new - old mass fraction ' // chem_isos% name(s% chem_id(j)), &
                     j, new_total-old_total
            end do
            write(*,*)
            j = jbad
            do k=2, nz
               if (comes_from(k) == comes_from(k-1)) cycle
               old_total = dot_product(xa_old(j,1:comes_from(k)-1),dq_old(1:comes_from(k)-1))
               if (old_total < 1d-9) cycle
               new_total = dot_product(xa(j,1:k-1),dq(1:k-1))
               write(*,2) 'partial new - old ' // chem_isos% name(s% chem_id(j)), k, &
                  new_total-old_total, new_total, old_total
            end do
            write(*,*)
            do k=415, nz
               write(*,'(a30,99i6)') 'cell_type(k)', k, cell_type(k), comes_from(k)
            end do
            write(*,*)
            write(*,2) 'xq', 439, xq(439)
            write(*,2) 'xq_old', 429, xq_old(429)
            write(*,2) 'dq_old', 429, dq_old(429)
            write(*,2) 'dq', 439, dq(439)
            write(*,*)
            write(*,2) 'xq', 424, xq(424)
            write(*,2) 'xq_old', 428, xq_old(428)
            write(*,2) 'dq_old', 428, dq_old(428)
            write(*,2) 'sum dq', 424, sum(dq(424:438))
            write(*,*)
            write(*,2) 'xq_old + dq_old', 428, xq_old(428) + dq_old(428)
            write(*,2) 'xq_old', 429, xq_old(429)
            write(*,*)
            write(*,2) 'xq + sum dq', 424, xq(424) + sum(dq(424:438))
            write(*,2) 'xq', 439, xq(439)
            write(*,*)
            write(*,1) 'sum dq_old', sum(dq_old(1:nz_old))

            write(*,2) 'dq_old', 427, dq_old(427)
            write(*,2) 'sum new', 416, sum(dq(416:423))
            write(*,2) 'dq_old - sum new', 427, dq_old(427) - sum(dq(416:423))
            write(*,2) 'dq_old', 428, dq_old(428)
            write(*,2) 'sum new', 424, sum(dq(424:438))
            write(*,2) 'dq_old - sum new', 428, dq_old(428) - sum(dq(424:438))
         end subroutine check_species_conservation


      end subroutine do_mesh_adjust


      subroutine do_prune_mesh_surface( &
            s, nz, nz_old, xh_old, xa_old, &
            j_rot_old, i_rot_old, omega_old, D_omega_old, am_nu_rot_old, &
            conv_vel_old, lnT_old, &
            dPdr_dRhodr_info_old, nu_ST_old, D_ST_old, D_DSI_old, D_SH_old, &
            D_SSI_old, D_ES_old, D_GSF_old, D_mix_old, &
            xh, xa, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old
         real(dp), dimension(:), pointer :: &
            j_rot_old, i_rot_old, omega_old, &
            D_omega_old, am_nu_rot_old, conv_vel_old, lnT_old, &
            dPdr_dRhodr_info_old, nu_ST_old, D_ST_old, D_DSI_old, D_SH_old, &
            D_SSI_old, D_ES_old, D_GSF_old, D_mix_old
         real(dp), dimension(:,:), pointer :: xh_old, xa_old
         real(dp), dimension(:,:), pointer :: xh, xa
         integer, intent(out) :: ierr

         integer :: skip, k, k_old, j

         include 'formats'

         ierr = 0
         skip = nz_old - nz
         if (skip < 0) then
            write(*,3) 'bad nz prune_surface do_prune_mesh_surface', nz_old, nz
            ierr = -1
            return
         end if

         do k = 1, nz
            k_old = k + skip
            do j = 1, s% nvar_hydro
               xh(j,k) = xh_old(j,k_old)
            end do
            do j = 1, s% species
               xa(j,k) = xa_old(j,k_old)
            end do
         end do

         call prune1(s% lnT, lnT_old, skip)
         call prune1(s% D_mix, D_mix_old, skip)

         if (s% rotation_flag) then
            call prune1(s% j_rot, j_rot_old, skip)
            call prune1(s% i_rot, i_rot_old, skip)
            call prune1(s% omega, omega_old, skip)
            call prune1(s% nu_ST, nu_ST_old, skip)
            call prune1(s% D_ST, D_ST_old, skip)
            call prune1(s% D_DSI, D_DSI_old, skip)
            call prune1(s% D_SH, D_SH_old, skip)
            call prune1(s% D_SSI, D_SSI_old, skip)
            call prune1(s% D_ES, D_ES_old, skip)
            call prune1(s% D_GSF, D_GSF_old, skip)
         end if

         if (s% D_omega_flag) then
            call prune1(s% D_omega, D_omega_old, skip)
         endif

         if (s% RTI_flag) then
            call prune1(s% dPdr_dRhodr_info, dPdr_dRhodr_info_old, skip)
         endif

         contains

         subroutine prune1(p,p_old,skip)
            real(dp), dimension(:), pointer :: p, p_old
            integer, intent(in) :: skip
            integer :: k
            do k=1,nz
               p(k) = p_old(k+skip)
            end do
         end subroutine prune1

      end subroutine do_prune_mesh_surface


      subroutine do_xh_pt_var( &
            s, i_var, nz, nz_old, nzlo, nzhi, comes_from, xh, xh_old, center_val, &
            xq, xq_old_plus1, xq_new, work, var_old_plus1, var_new, ierr)
         use interp_1d_def
         use interp_1d_lib
         type (star_info), pointer :: s
         integer, intent(in) :: i_var, nz, nz_old, nzlo, nzhi, comes_from(:)
         real(dp),intent(in) :: center_val
         real(dp), dimension(:,:), pointer :: xh, xh_old
         real(dp), dimension(:), pointer :: work
         real(dp), dimension(:) :: &
            xq, xq_old_plus1, var_old_plus1, var_new, xq_new
         integer, intent(out) :: ierr

         integer :: n, k

         include 'formats'

         ierr = 0
         n = nzhi - nzlo + 1

         do k=1,nz_old
            var_old_plus1(k) = xh_old(i_var,k)
         end do
         var_old_plus1(nz_old+1) = center_val

         call interpolate_vector( &
            nz_old+1, xq_old_plus1, n, xq_new, &
            var_old_plus1, var_new, interp_pm, nwork, work, &
            'mesh_adjust do_xh_pt_var', ierr)
         if (ierr /= 0) then
            return
         end if

         do k=nzlo,nzhi
            xh(i_var,k) = max(0d0,var_new(k+1-nzlo))
         end do

         n = nzlo - 1
         if (n > 0) then
            do k=1,n
               xh(i_var,k) = xh_old(i_var,k)
            end do
         end if

         if (nzhi < nz) then
            n = nz - nzhi - 1 ! nz-n = nzhi+1
            do k=0,n
               xh(i_var,nz-k) = xh_old(i_var,nz_old-k)
            end do
         end if

      end subroutine do_xh_pt_var


      subroutine do_L( &
            s, nz, nz_old, nzlo, nzhi, comes_from, xh, xh_old, &
            xq, xq_old_plus1, xq_new, work, L_old_plus1, L_new, ierr)
         use interp_1d_def
         use interp_1d_lib
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old, nzlo, nzhi, comes_from(:)
         real(dp), dimension(:,:), pointer :: xh, xh_old
         real(dp), dimension(:), pointer :: work
         real(dp), dimension(:) :: &
            xq, xq_old_plus1, L_old_plus1, L_new, xq_new
         integer, intent(out) :: ierr

         integer :: n, i_lum, k

         include 'formats'

         ierr = 0
         i_lum = s% i_lum
         if (i_lum == 0) return
         n = nzhi - nzlo + 1

         do k=1,nz_old
            L_old_plus1(k) = xh_old(i_lum,k)
         end do
         L_old_plus1(nz_old+1) = s% L_center

         call interpolate_vector( &
               nz_old+1, xq_old_plus1, n, xq_new, &
               L_old_plus1, L_new, interp_pm, nwork, work, &
               'mesh_adjust do_L', ierr)
         if (ierr /= 0) then
            return

            write(*,*) 'interpolate_vector failed in do_L for remesh'
            stop 'debug: mesh adjust: do_L'
         end if

         do k=nzlo,nzhi
            xh(i_lum,k) = L_new(k+1-nzlo)
         end do

         n = nzlo - 1
         if (n > 0) then
            do k=1,n
               xh(i_lum,k) = xh_old(i_lum,k)
            end do
         end if

         if (nzhi < nz) then
            n = nz - nzhi - 1 ! nz-n = nzhi+1
            do k=0,n
               xh(i_lum,nz-k) = xh_old(i_lum,nz_old-k)
            end do
         end if

      end subroutine do_L


      subroutine do_conv_vel( &
            s, nz, nz_old, nzlo, nzhi, comes_from, xh, xh_old, &
            xq, xq_old_plus1, xq_new, work, cv_old_plus1, cv_new, ierr)
         use interp_1d_def
         use interp_1d_lib
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old, nzlo, nzhi, comes_from(:)
         real(dp), dimension(:,:), pointer :: xh, xh_old
         real(dp), dimension(:), pointer :: work
         real(dp), dimension(:) :: &
            xq, xq_old_plus1, cv_old_plus1, cv_new, xq_new
         integer, intent(out) :: ierr

         integer :: n, i_ln_cvpv0, k

         include 'formats'

         ierr = 0
         i_ln_cvpv0 = s% i_ln_cvpv0
         if (i_ln_cvpv0 == 0) return
         n = nzhi - nzlo + 1

         do k=1,nz_old
            cv_old_plus1(k) = max(0d0,exp(xh_old(i_ln_cvpv0,k))-s% conv_vel_v0)
         end do
         cv_old_plus1(nz_old+1) = cv_old_plus1(nz_old)

         call interpolate_vector( &
               nz_old+1, xq_old_plus1, n, xq_new, &
               cv_old_plus1, cv_new, interp_pm, nwork, work, &
               'mesh_adjust do_conv_vel', ierr)
         if (ierr /= 0) then
            return

            write(*,*) 'interpolate_vector failed in do_conv_vel for remesh'
            stop 'debug: mesh adjust: do_conv_vel'
         end if

         do k=nzlo,nzhi
            xh(i_ln_cvpv0,k) = log(cv_new(k+1-nzlo)+s% conv_vel_v0)
         end do

         n = nzlo - 1
         if (n > 0) then
            do k=1,n
               xh(i_ln_cvpv0,k) = xh_old(i_ln_cvpv0,k)
            end do
         end if

         if (nzhi < nz) then
            n = nz - nzhi - 1 ! nz-n = nzhi+1
            do k=0,n
               xh(i_ln_cvpv0,nz-k) = xh_old(i_ln_cvpv0,nz_old-k)
            end do
         end if

      end subroutine do_conv_vel


      subroutine do_alpha_RTI( &
            s, nz, nz_old, nzlo, nzhi, comes_from, xh, xh_old, &
            xq, xq_old_plus1, xq_new, work, alpha_RTI_old_plus1, alpha_RTI_new, ierr)
         use interp_1d_def
         use interp_1d_lib
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old, nzlo, nzhi, comes_from(:)
         real(dp), dimension(:,:), pointer :: xh, xh_old
         real(dp), dimension(:), pointer :: work
         real(dp), dimension(:) :: &
            xq, xq_old_plus1, alpha_RTI_old_plus1, alpha_RTI_new, xq_new
         integer, intent(out) :: ierr

         integer :: n, i_alpha_RTI, k

         include 'formats'

         ierr = 0
         i_alpha_RTI = s% i_alpha_RTI
         n = nzhi - nzlo + 1

         do k=1,nz_old
            alpha_RTI_old_plus1(k) = xh_old(i_alpha_RTI,k)
         end do
         alpha_RTI_old_plus1(nz_old+1) = 0

         call interpolate_vector( &
               nz_old+1, xq_old_plus1, n, xq_new, &
               alpha_RTI_old_plus1, alpha_RTI_new, interp_pm, nwork, work, &
               'mesh_adjust do_alpha_RTI', ierr)
         if (ierr /= 0) then
            return
         end if

         do k=nzlo,nzhi
            xh(i_alpha_RTI,k) = max(0d0,alpha_RTI_new(k+1-nzlo))
         end do

         n = nzlo - 1
         if (n > 0) then
            do k=1,n
               xh(i_alpha_RTI,k) = xh_old(i_alpha_RTI,k)
            end do
         end if

         if (nzhi < nz) then
            n = nz - nzhi - 1 ! nz-n = nzhi+1
            do k=0,n
               xh(i_alpha_RTI,nz-k) = xh_old(i_alpha_RTI,nz_old-k)
            end do
         end if

      end subroutine do_alpha_RTI


      subroutine do_et( &  
            ! this is not being careful to conserve et.  may need to improve.
            ! similarly not taking et into account in adjusting lnT to conserve energy.
            s, nz, nz_old, nzlo, nzhi, comes_from, xh, xh_old, &
            xq, xq_old_plus1, xq_new, work, et_old_plus1, et_new, ierr)
         use interp_1d_def
         use interp_1d_lib
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old, nzlo, nzhi, comes_from(:)
         real(dp), dimension(:,:), pointer :: xh, xh_old
         real(dp), dimension(:), pointer :: work
         real(dp), dimension(:) :: &
            xq, xq_old_plus1, et_old_plus1, et_new, xq_new
         integer, intent(out) :: ierr

         integer :: n, i_et, k

         include 'formats'

         ierr = 0
         i_et = s% i_et
         n = nzhi - nzlo + 1

         do k=1,nz_old
            et_old_plus1(k) = xh_old(i_et,k)
         end do
         et_old_plus1(nz_old+1) = 0

         call interpolate_vector( &
               nz_old+1, xq_old_plus1, n, xq_new, &
               et_old_plus1, et_new, interp_pm, nwork, work, &
               'mesh_adjust do_et', ierr)
         if (ierr /= 0) then
            return
         end if

         do k=nzlo,nzhi
            xh(i_et,k) = max(0d0,et_new(k+1-nzlo))
         end do

         n = nzlo - 1
         if (n > 0) then
            do k=1,n
               xh(i_et,k) = xh_old(i_et,k)
            end do
         end if

         if (nzhi < nz) then
            n = nz - nzhi - 1 ! nz-n = nzhi+1
            do k=0,n
               xh(i_et,nz-k) = xh_old(i_et,nz_old-k)
            end do
         end if

      end subroutine do_et


      subroutine do_interp_pt_val( &
            s, nz, nz_old, nzlo, nzhi, val, val_old, center_val, &
            xq, xq_old_plus1, xq_new, force_non_negative, &
            work, val_old_plus1, val_new, ierr)
         use interp_1d_def
         use interp_1d_lib
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old, nzlo, nzhi
         real(dp), dimension(:), pointer :: val, val_old
         real(dp), intent(in) :: center_val
         real(dp), dimension(:), pointer :: work
         real(dp), dimension(:) :: &
            xq, xq_old_plus1, xq_new, val_old_plus1, val_new
         logical, intent(in) :: force_non_negative
         integer, intent(out) :: ierr
         integer :: n, k

         include 'formats'

         ierr = 0
         n = nzhi - nzlo + 1

         do k=1,nz_old
            val_old_plus1(k) = val_old(k)
         end do
         val_old_plus1(nz_old+1) = center_val

         call interpolate_vector( &
               nz_old+1, xq_old_plus1, n, xq_new, &
               val_old_plus1, val_new, interp_pm, nwork, work, &
               'mesh_adjust do_interp_pt_val', ierr)
         if (ierr /= 0) then
            return
         end if

         do k=nzlo,nzhi
            val(k) = val_new(k+1-nzlo)
         end do

         n = nzlo - 1
         if (n > 0) then
            do k=1,n
               val(k) = val_old(k)
            end do
         end if

         if (nzhi < nz) then
            n = nz - nzhi - 1 ! nz-n = nzhi+1
            do k=0,n
               val(nz-k) = val_old(nz_old-k)
            end do
         end if

         if (force_non_negative) then
            do k=nzlo,nzhi
               if (val(k) < 0) val(k) = 0
            end do
         end if

      end subroutine do_interp_pt_val


      subroutine do_interp_cell_val( &
            s, nz, nz_old, nzlo, nzhi, val_new_out, val_old, &
            xq, xq_old_plus1, dq, dq_old, work, val_old_plus1, val_new, ierr)
         use interp_1d_def
         use interp_1d_lib
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old, nzlo, nzhi
         real(dp), dimension(:), pointer :: val_new_out, val_old
         real(dp), dimension(:), pointer :: work
         real(dp), dimension(:) :: &
            xq, xq_old_plus1, dq, dq_old, val_old_plus1, val_new
         integer, intent(out) :: ierr

         real(dp), pointer, dimension(:) :: &
            mid_xq_new, mid_xq_old_plus1
         integer :: n, i, j, k

         ierr = 0
         n = nzhi - nzlo + 1
         call do_alloc(ierr)
         if (ierr /= 0) return

         do k=1,nz_old
            val_old_plus1(k) = val_old(k)
            mid_xq_old_plus1(k) = xq_old_plus1(k) + 0.5d0*dq_old(k)
         end do
         val_old_plus1(nz_old+1) = val_old_plus1(nz_old)
         mid_xq_old_plus1(nz_old+1) = 1
         do i=1,n
            mid_xq_new(i) = xq(nzlo+i-1) + 0.5d0*dq(nzlo+i-1)
         end do

         call interpolate_vector( &
               nz_old+1, mid_xq_old_plus1, n, mid_xq_new, &
               val_old_plus1, val_new, interp_pm, nwork, work, &
               'mesh_adjust do_interp_cell_val', ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if

         do i=1,n
            val_new_out(nzlo+i-1) = val_new(i)
         end do

         n = nzlo - 1
         if (n > 0) then
            do i=1,n
               val_new_out(i) = val_old(i)
            end do
         end if

         if (nzhi < nz) then
            n = nz - nzhi - 1 ! nz-n = nzhi+1
            do i=0,n
               val_new_out(nz-i) = val_old(nz_old-i)
            end do
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
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
                mid_xq_old_plus1, nz_old+1, nz_alloc_extra, 'mesh_adjust', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
                mid_xq_new, n, nz_alloc_extra, 'mesh_adjust', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays

      end subroutine do_interp_cell_val


      subroutine do_lnR_and_lnd( &
            s, nz, nz_old, nzlo, nzhi, cell_type, comes_from, &
            xh, xh_old, xmstar, lnd_old, lnPgas_old, &
            dqbar, dqbar_old, old_r, old_m, old_rho, &
            dq, dq_old, xq, xq_old_plus1, density_new, work, &
            Vol_old_plus1, Vol_new, new_r, Vol_init, &
            interp_Vol_new, interp_xq, density_init, ierr)
         use interp_1d_def
         use interp_1d_lib
         use num_lib
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old, nzlo, nzhi, comes_from(:)
         integer :: cell_type(:)
         real(dp), dimension(:,:), pointer :: xh, xh_old
         real(dp), intent(in) :: xmstar
         real(dp), dimension(:), pointer :: work
         real(dp), dimension(:) :: &
            lnd_old, lnPgas_old, dqbar, dqbar_old, old_r, old_m, old_rho, &
            xq, dq, dq_old, xq_old_plus1, density_new, &
            Vol_old_plus1, Vol_new, new_r, Vol_init, &
            interp_Vol_new, interp_xq, density_init
         integer, intent(out) :: ierr

         integer :: k, from_k, kk, n, interp_lo, interp_hi, interp_n, &
            num_revise, i_lnR, i_lnd
         real(dp) :: Vol_min, Vol_max, cell_Vol, Vol_center, Vm1, V00, Vp1

         logical, parameter :: dbg = .false., trace_PE_residual = .false.

         include 'formats'

         ! NOTE: for interpolating volume, need to add point at center

         ierr = 0
         i_lnR = s% i_lnR
         i_lnd = s% i_lnd

         interp_lo = max(1, nzlo-1)
         interp_hi = min(nz, nzhi+1)
         interp_n = interp_hi - interp_lo + 1

         do k=1,nz_old
            Vol_old_plus1(k) = (pi4/3)*old_r(k)*old_r(k)*old_r(k)
         end do
         Vol_center = (pi4/3)*s% R_center*s% R_center*s% R_center
         Vol_old_plus1(nz_old+1) = Vol_center

         ! testing -- check for Vol_old_plus1 strictly decreasing
         do k = 2, nz_old+1
            if (Vol_old_plus1(k) >= Vol_old_plus1(k-1)) then
               ierr = -1
               if (.not. dbg) return
               write(*,3) 'bad old vol', k, nz_old
               write(*,1) 'Vol_old_plus1(k)', Vol_old_plus1(k)
               write(*,1) 'Vol_old_plus1(k-1)', Vol_old_plus1(k-1)
               write(*,*)
               stop 'debug: mesh adjust: do_lnR_and_lnd'
            end if
         end do

         ! testing -- check for q strictly decreasing
         do k = 2, nz
            if (xq(k) <= xq(k-1)) then
               ierr = -1
               if (.not. dbg) return

               write(*,3) 'bad xq', k, nz, xq(k), xq(k-1)
               stop 'debug: mesh adjust: do_lnR_and_lnd'
            end if
         end do

         do k=1,interp_n
            interp_xq(k) = xq(interp_lo+k-1)
         end do

         call interpolate_vector( &
               nz_old+1, xq_old_plus1, interp_n, interp_xq, Vol_old_plus1, &
               interp_Vol_new, interp_pm, nwork, work, 'mesh_adjust do_lnR_and_lnd', ierr)
         if (ierr /= 0) then
            if (.not. dbg) return
            write(*,*) 'failed in interpolate_vector'
            stop 'debug: mesh_adjust'
         end if

         do k=1,interp_n
            Vol_new(interp_lo+k-1) = interp_Vol_new(k)
         end do

         if (Vol_new(interp_lo+1) >= Vol_new(interp_lo)) then
            Vol_new(interp_lo+1) = (Vol_new(interp_lo) + Vol_new(interp_lo+2))/2
            if (dbg) write(*,2) 'fix Vol_new at lo+1', interp_lo+1, Vol_new(interp_lo+1)
            if (Vol_new(interp_lo+1) >= Vol_new(interp_lo)) then
               ierr = -1
               if (.not. dbg) return
               write(*,*) '(Vol_new(interp_lo+1) >= Vol_new(interp_lo))'
               stop 'debug: mesh_adjust'
            end if
         end if

         do k = interp_lo+1, interp_hi-1
            if (Vol_new(k+1) >= Vol_new(k) .or. Vol_new(k) >= Vol_new(k-1)) then
               if (dbg) write(*,2) 'fix interpolated Vol_new', &
                  k, Vol_new(k+1), Vol_new(k), Vol_new(k-1)
               Vol_min = minval(Vol_new(k-1:k+1))
               Vol_max = maxval(Vol_new(k-1:k+1))
               if (Vol_min == Vol_max .or. is_bad(Vol_min) .or. is_bad(Vol_max)) then
                  ierr = -1
                  if (s% stop_for_bad_nums) then
                     write(*,2) 'Vol_min', k, Vol_min
                     write(*,2) 'Vol_max', k, Vol_max
                     stop 'mesh_adjust'
                  end if
                  if (.not. dbg) return
                  write(*,1) 'Vol_min', Vol_min
                  write(*,1) 'Vol_max', Vol_max
                  stop 'debug: mesh_adjust'
               end if
               Vm1 = Vol_new(k-1)
               V00 = Vol_new(k)
               Vp1 = Vol_new(k+1)
               Vol_new(k-1) = Vol_max
               Vol_new(k) = (Vol_max + Vol_min)/2
               Vol_new(k+1) = Vol_min
               if (dbg) write(*,2) 'new Vol_new',  &
                  k, Vol_new(k+1), Vol_new(k), Vol_new(k-1)
               if (Vol_new(k+1) >= Vol_new(k) .or. Vol_new(k) >= Vol_new(k-1)) then
                  ierr = -1
                  if (.not. dbg) return
                  write(*,1) 'Vol_new(k-1)', Vol_new(k-1)
                  write(*,1) 'Vol_new(k)', Vol_new(k)
                  write(*,1) 'Vol_new(k+1)', Vol_new(k+1)
                  stop 'debug: do_lnR_and_lnd in mesh adjust: interpolation gave non-pos volume'
               end if
            end if
         end do

         call set1_new_r(nzlo)
         do k = nzlo, min(nzhi,nz-1)
            if (ierr /= 0) cycle

            call set1_new_r(k+1)

            if (cell_type(k) == unchanged_type) then
               xh(i_lnd,k) = lnd_old(comes_from(k))
               density_new(k) = old_rho(comes_from(k))
               cycle
            end if

            if (new_r(k) <= new_r(k+1)) then
               if (dbg) then
                  write(*,*) 'do_lnR_and_lnd: (new_r(k) <= new_r(k+1))'
                  stop
               end if
               ierr = -1; cycle
            end if

            cell_Vol = Vol_new(k)-Vol_new(k+1)
            if (cell_Vol <= 0) then
               if (dbg) then
                  write(*,2) 'do_lnR_and_lnd: cell_Vol <= 0', k
               end if
               ierr = -1; cycle
            end if
            if (dq(k) <= 0) then
               if (dbg) then
                  write(*,2) 'do_lnR_and_lnd: dq(k) <= 0', k
               end if
               ierr = -1; cycle
            end if
            density_new(k) = xmstar*dq(k)/cell_Vol
            xh(i_lnd,k) = log(density_new(k))

         end do

         if (ierr /= 0) then
            if (.not. dbg) return
            stop 'debug: failed in mesh adjust do_lnR_and_lnd'
         end if

         n = nzlo - 1
         if (n > 0) then
            do k=1,n
               new_r(k) = old_r(k)
               density_new(k) = old_rho(k)
               Vol_new(k) = 4d0/3d0*pi*new_r(k)*new_r(k)*new_r(k)
            end do
         end if

         if (nzhi < nz) then
            n = nz - nzhi - 1 ! nz-n = nzhi+1
            do k=0,n
               new_r(nz-k) = old_r(nz_old-k)
               density_new(nz-k) = old_rho(nz_old-k)
               Vol_new(nz-k) = 4d0/3d0*pi*new_r(nz-k)*new_r(nz-k)*new_r(nz-k)
            end do
         else ! nzhi == nz
            density_new(nz) = xmstar*dq(nz)/(Vol_new(nz) - Vol_center)
            new_r(nz) = pow(Vol_new(nz)/(pi4/3), 1d0/3d0)

            if (dbg) then
               write(*,2) 'old_rho(nz_old)', nz_old, old_rho(nz_old)
               write(*,2) 'density_new(nz)', nz, density_new(nz)
            end if

         end if

         do k=1,nz
            Vol_init(k) = Vol_new(k)
            density_init(k) = density_new(k)
         end do

         do k=1,nz
            from_k = comes_from(k)
            if (new_r(k) == old_r(from_k)) then
               xh(i_lnR,k) = xh_old(i_lnR,from_k)
            else
               xh(i_lnR,k) = log(new_r(k))
            end if
            if (density_new(k) == old_rho(from_k)) then
               xh(i_lnd,k) = lnd_old(from_k)
            else 
               xh(i_lnd,k) = log(density_new(k))
            end if
         end do


         contains

         subroutine set1_new_r(k)
            integer, intent(in) :: k
            include 'formats'
            if (cell_type(k) == unchanged_type) then
               new_r(k) = old_r(comes_from(k))
            else
               new_r(k) = pow(Vol_new(k)/(pi4/3), 1d0/3d0)
            end if
         end subroutine set1_new_r


      end subroutine do_lnR_and_lnd


      subroutine do_xa( &
            s, nz, nz_old, k, species, cell_type, comes_from, &
            xa, xa_old, xa_c0, xa_c1, xa_c2, &
            xq, dq, xq_old,  dq_old, mesh_adjust_use_quadratic, ierr)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old, species, k, cell_type(:), comes_from(:)
         real(dp), dimension(:,:), pointer :: xa, xa_old
         real(dp), dimension(:,:) :: xa_c0, xa_c1, xa_c2
         real(dp), dimension(:), pointer :: xq, dq, xq_old, dq_old
         logical, intent(in) :: mesh_adjust_use_quadratic
         integer, intent(out) :: ierr

         integer :: j, jj, k_old, k_old_last, kdbg, order
         real(dp) :: xq_outer, cell_dq, xa_sum, total(species)
         logical :: dbg_get_integral

         include 'formats'

         ierr = 0

         kdbg = -1074

         if (mesh_adjust_use_quadratic) then
            order = 2
         else
            order = 1
         end if

         if (cell_type(k) == unchanged_type .or. &
               cell_type(k) == revised_type) then
            do j=1,species
               xa(j,k) = xa_old(j,comes_from(k))
            end do
            return
         end if

         xq_outer = xq(k)
         if (k == nz) then
            cell_dq = 1 - xq_outer
         else
            cell_dq = dq(k)
         end if

         k_old = max(comes_from(k)-1,1)

         ! sum the old abundances between xq_outer and xq_inner
         dbg_get_integral = .false.
         total(:) = 0
         do j=1,species
            dbg_get_integral = (k == kdbg) .and. (j == 1) ! h1
            if (dbg_get_integral) write(*,2) trim(chem_isos% name(s% chem_id(j)))
            call get_xq_integral( &
               k_old, nz_old, xq_old, xq_outer, cell_dq, &
               order, xa_c0(:,j), xa_c1(:,j), xa_c2(:,j), &
               total(j), dbg_get_integral, k_old_last, ierr)
         end do

         xa(:,k) = total(:)/cell_dq

         if (k == kdbg) then
            do j=1,species
               write(*,2) 'new ' // trim(chem_isos% name(s% chem_id(j))), k, xa(j,k)
            end do
         end if

         do j=1,species
            if (xa(j,k) > 1 + 1d-8 .or. xa(j,k) < -1d-8) then
               ierr = -1
               return

               do jj=1,species
                  write(*,1) 'xa ' // trim(chem_isos% name(s% chem_id(jj))), xa(jj,k)
               end do
               write(*,*)
               write(*,2) 'sum xa', k, sum(xa(:,k))
               write(*,*)
               write(*,2) 'xa ' // trim(chem_isos% name(s% chem_id(j))), k, xa(j,k)
               write(*,*)
               write(*,2) 'xq_outer', k, xq_outer
               write(*,2) 'xq_inner', k, xq_outer + cell_dq
               write(*,2) 'cell_dq', k, cell_dq
               write(*,*)
               write(*,2) 'xq_old(k_old)', k_old, xq_old(k_old)
               write(*,2) 'xq_inner(k_old)', k_old, xq_old(k_old)+dq_old(k_old)
               write(*,2) 'dq_old(k_old)', k_old, dq_old(k_old)
               write(*,*)
               write(*,2) 'xa_c0(k_old,j)', k_old, xa_c0(k_old,j)
               write(*,2) 'xa_c1(k_old,j)', k_old, xa_c1(k_old,j)
               write(*,2) 'xa_c2(k_old,j)', k_old, xa_c2(k_old,j)
               write(*,*)
               write(*,2) 'old outer', k_old, xa_c0(k_old,j) + xa_c1(k_old,j)*dq_old(k_old)/2
               write(*,2) 'old inner', k_old, xa_c0(k_old,j) - xa_c1(k_old,j)*dq_old(k_old)/2
               write(*,*)
               stop 'debug: mesh adjust: do_xa'
            end if
         end do

         xa_sum = sum(xa(:,k))
         !write(*,1) 'xa_sum', xa_sum

         if (is_bad(xa_sum)) then
            ierr = -1
            if (s% stop_for_bad_nums) then
               write(*,2) 'xa_sum', k, xa_sum
               stop 'mesh adjust: do_xa'
            end if
            return

            write(*,*) 'xa_sum', xa_sum
            write(*,*) 'bug in revise mesh, do_xa bad num: k', k
            stop 'debug: mesh adjust: do_xa'
         end if

         if (abs(1-xa_sum) > 1d-3) then
            ierr = -1
            return
         end if

         xa(:,k) = xa(:,k) / xa_sum

      end subroutine do_xa


      subroutine do1_lnT( &
            s, nz_old, k, &
            species, cell_type, comes_from, &
            xa, xh, xh_old, &
            xq, dq, xq_old, dq_old, eta_old, energy_old, lnT_old, &
            specific_PE_old, specific_KE_old, et_old, &
            density_new, energy_new, ierr)
         use eos_def
         use star_utils, only: set_rmid, cell_specific_PE, cell_specific_KE
         type (star_info), pointer :: s
         integer, intent(in) :: nz_old, k, species, cell_type(:), comes_from(:)
         real(dp), dimension(:,:), pointer :: xa, xh, xh_old
         real(dp), dimension(:) :: &
            xq, dq, xq_old, dq_old, eta_old, energy_old, lnT_old, &
            specific_PE_old, specific_KE_old, et_old, density_new, energy_new
         integer, intent(out) :: ierr

         integer :: k_old, k_old_last, i_lnT, lnT_order, energy_order
         real(dp) :: &
            Rho, logRho, xq_outer, cell_dq, avg_energy, avg_PE, avg_KE, &
            new_PE, new_KE, max_delta_energy, delta_energy, revised_energy, &
            sum_lnT, avg_lnT, new_lnT, sum_energy, new_xa(species), &
            d_dlnR00, d_dlnRp1, d_dv00, d_dvp1
         logical :: dbg_get_integral

         include 'formats'

         ierr = 0
         i_lnT = s% i_lnT
         new_xa(:) = xa(:,k)
         k_old = comes_from(k)

         if (cell_type(k) == unchanged_type) then
            xh(i_lnT, k) = xh_old(i_lnT, k_old)
            energy_new(k) = energy_old(k_old)
            if (is_bad(energy_old(k_old))) then
               write(*,2) 'energy_old(k_old)', k_old, energy_old(k_old)
               stop 'debug: mesh adjust: do1_lnT'
            end if
            return
         end if

         xq_outer = xq(k)
         cell_dq = dq(k)

         if (cell_type(k) == revised_type) then
            avg_lnT = xh_old(i_lnT, k_old)
         else ! find average lnT between xq_outer and xq_inner
            call get_old_value_integral( &
               k, k_old, nz_old, xq_old, dq_old, xq_outer, cell_dq, &
               lnT_old, sum_lnT, dbg, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'get_old_value_integral lnT failed for do1_lnT'
               if (.not. dbg) return
               stop 'debug: mesh adjust: do1_lnT'
            end if
            avg_lnT = sum_lnT/cell_dq
         end if

         if (is_bad(avg_lnT) .or. avg_lnT < 0 .or. avg_lnT > 100) then
            ierr = -1
            if (s% stop_for_bad_nums) then
               write(*,2) 'avg_lnT', k, avg_lnT
               stop 'mesh adjust: do1_lnT'
            end if
            return
         end if

         if (.not. s% mesh_adjust_get_T_from_E) then
            xh(i_lnT, k) = avg_lnT
            energy_new(k) = energy_old(k_old)
            if (is_bad(energy_old(k_old))) then
               write(*,2) 'energy_old(k_old)', k_old, energy_old(k_old)
               stop 'debug: mesh adjust: do1_lnT'
            end if
            return
         end if

         if (eta_old(k_old) >= eta_limit) then
            xh(i_lnT, k) = avg_lnT
            energy_new(k) = energy_old(k_old)
            if (is_bad(energy_old(k_old))) then
               write(*,2) 'eta_old(k_old)', k_old, eta_old(k_old)
               write(*,2) 'energy_old(k_old)', k_old, energy_old(k_old)
               stop 'debug: mesh adjust: do1_lnT'
            end if
            return
         end if

         if (dbg) write(*,2) 'eta_old(k_old)', k_old, eta_old(k_old)

         if (cell_type(k) == revised_type) then
            avg_energy = energy_old(k_old)
         else ! find average internal energy between q_outer and q_inner
            call get_old_value_integral( &
               k, k_old, nz_old, xq_old, dq_old, xq_outer, cell_dq, &
               energy_old, sum_energy, dbg, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'get_old_value_integral failed for do1_lnT'
               if (.not. dbg) return
               stop 'debug: mesh adjust: energy_old do1_lnT'
            end if
            avg_energy = sum_energy/cell_dq
         end if
         
         if (s% max_rel_delta_IE_for_mesh_total_energy_balance == 0d0) then
         
            energy_new(k) = avg_energy
         
         else

            if (cell_type(k) == revised_type) then
               avg_PE = specific_PE_old(k_old)
            else ! find average potential energy between q_outer and q_inner
               call get_old_value_integral( &
                  k, k_old, nz_old, xq_old, dq_old, xq_outer, cell_dq, &
                  specific_PE_old, sum_energy, dbg, ierr)
               if (ierr /= 0) then
                  if (dbg) write(*,*) 'get_old_value_integral failed for do1_lnT'
                  if (.not. dbg) return
                  stop 'debug: mesh adjust: specific_PE_old do1_lnT'
               end if
               avg_PE = sum_energy/cell_dq
            end if

            if (cell_type(k) == revised_type) then
               avg_KE = specific_KE_old(k_old)
            else ! find average kinetic energy between q_outer and q_inner
               call get_old_value_integral( &
                  k, k_old, nz_old, xq_old, dq_old, xq_outer, cell_dq, &
                  specific_KE_old, sum_energy, dbg, ierr)
               if (ierr /= 0) then
                  if (dbg) write(*,*) 'get_old_value_integral failed for do1_lnT'
                  if (.not. dbg) return
                  stop 'debug: mesh adjust: specific_KE_old do1_lnT'
               end if
               avg_KE = sum_energy/cell_dq
            end if
         
            s% r(k) = exp(s% xh(s% i_lnr,k))
            if (k < s% nz) s% r(k+1) = exp(s% xh(s% i_lnr,k+1))
            if (ierr /= 0) return
            new_PE = cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)
            if (s% u_flag) then
               s% u(k) = s% xh(s% i_u,k)
            else if (s% v_flag) then
               s% v(k) = s% xh(s% i_v,k)
               if (k < s% nz) s% v(k+1) = s% xh(s% i_v,k+1)
            end if
            new_KE = cell_specific_KE(s,k,d_dv00,d_dvp1)

            max_delta_energy = avg_energy*s% max_rel_delta_IE_for_mesh_total_energy_balance
            delta_energy = avg_PE + avg_KE - (new_PE + new_KE)
            if (abs(delta_energy) > max_delta_energy) then
               if (s% show_mesh_changes) &
                  write(*,2) 'remesh: delta_energy too large to fix completely', k, &
                     delta_energy, max_delta_energy, &
                     specific_PE_old(k_old), specific_KE_old(k_old), et_old(k_old)
               delta_energy = sign(max_delta_energy,delta_energy)
            end if
            energy_new(k) = avg_energy + delta_energy
            
            if (energy_new(k) <= 0d0) then
               write(*,2) 'energy_new(k) <= 0d0', k, energy_new(k), avg_energy
               energy_new(k) = avg_energy
            end if
            
         end if

         ! call eos to calculate lnT from new internal energy

         Rho = density_new(k)
         logRho = log10(Rho)
         call set_lnT_for_energy( &
            s, k, &
            s% net_iso(ih1), s% net_iso(ihe3), s% net_iso(ihe4), species, new_xa, &
            Rho, logRho, energy_new(k), avg_lnT, new_lnT, revised_energy, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'set_lnT_for_energy failed', k
            new_lnT = avg_lnT
            energy_new(k) = energy_old(k_old)
            ierr = 0
         end if

         xh(i_lnT,k) = new_lnT

         if (ierr /= 0) then
            write(*,2) 'mesh_adjust do1_lnT ierr', ierr
            stop 'do1_lnT'
         end if

      end subroutine do1_lnT


      subroutine get_old_value_integral( &
            k_new, k_old_in, nz_old, xq_old, dq_old, xq_outer, dq_range, &
            value_old, integral, dbg, ierr)
         integer, intent(in) :: k_new, k_old_in, nz_old
         real(dp), intent(in) :: xq_old(:), dq_old(:), xq_outer, dq_range
         real(dp), intent(in), dimension(:) :: value_old
         real(dp), intent(out) :: integral
         logical, intent(in) :: dbg
         integer, intent(out) :: ierr
         real(dp), pointer :: p(:,:)
         nullify(p)
         call get_old_integral( &
            k_new, k_old_in, nz_old, xq_old, dq_old, xq_outer, dq_range, &
            value_old, p, integral, dbg, ierr)
      end subroutine get_old_value_integral


      subroutine get_old_integral( &
            k_new, k_old_in, nz_old, xq_old, dq_old, xq_outer, dq_range, &
            value_old, xh_old, integral, dbg, ierr)
         integer, intent(in) :: k_new, k_old_in, nz_old
         real(dp), intent(in) :: xq_old(:), dq_old(:), xq_outer, dq_range
         real(dp), intent(in) :: value_old(:)
         real(dp), intent(in), pointer :: xh_old(:,:)
         real(dp), intent(out) :: integral
         logical, intent(in) :: dbg
         integer, intent(out) :: ierr

         integer :: k, k_old
         real(dp) :: xq_inner, sum_dqs, old_xq_outer, old_xq_inner, &
            dq_overlap, val

         include 'formats'

         if (dbg) write(*,*)

         ierr = 0
         
         k_old = k_old_in
         ! move starting k_old if necessary
         do
            if (k_old <= 1) exit
            if (xq_old(k_old) <= xq_outer) exit
            k_old = k_old - 1
         end do

         xq_inner = xq_outer + dq_range
         old_xq_inner = xq_old(k_old)
         sum_dqs = 0
         integral = 0d0

         if (dbg) write(*,*)
         if (dbg) write(*,3) 'k_new k_old xq_outer xq_inner dq_range', &
            k_new, k_old, xq_outer, xq_inner, dq_range

         do k = k_old, nz_old

            if (dq_range <= sum_dqs) exit
            old_xq_outer = old_xq_inner
            if (k == nz_old) then
               old_xq_inner = 1
            else
               old_xq_inner = xq_old(k+1)
            end if

            if (dbg) write(*,3) 'k_new k_old old_xq_outer old_xq_inner', &
               k_new, k, old_xq_outer, old_xq_inner

            val = value_old(k)

            if (old_xq_inner <= xq_inner .and. old_xq_outer >= xq_outer) then

               if (dbg) write(*,1) 'entire old cell is in new range'

               sum_dqs = sum_dqs + dq_old(k)
               integral = integral + val*dq_old(k)

            else if (old_xq_inner >= xq_inner .and. old_xq_outer <= xq_outer) then

               if (dbg) write(*,1) 'entire new range is in this old cell'

               sum_dqs = dq_range
               integral = val*dq_range

            else ! only use the part of old cell that is in new range

               if (xq_inner <= old_xq_inner) then

                  if (dbg) write(*,1) 'last part of the new range'                  
                  
                  integral = integral + val*(dq_range - sum_dqs)
                  sum_dqs = dq_range

               else ! partial overlap -- general case

                  dq_overlap = max(0d0, old_xq_inner - xq_outer)
                  sum_dqs = sum_dqs + dq_overlap
                  integral = integral + val*dq_overlap

                  if (dbg) write(*,1) 'partial overlap'

               end if

            end if

         end do

         if (dbg) write(*,2) 'integral/dq_range', &
            k_new, integral/dq_range, integral, dq_range
         if (dbg) write(*,*)

      end subroutine get_old_integral


      subroutine set_lnT_for_energy( &
            s, k, h1, he3, he4, species, xa, &
            Rho, logRho, energy, lnT_guess, lnT, result_energy, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, h1, he3, he4, species
         real(dp), intent(in) :: xa(species), Rho, logRho, energy, lnT_guess
         real(dp), intent(out) :: lnT, result_energy
         integer, intent(out) :: ierr
         call set_lnT_for_energy_with_tol( &
            s, k, h1, he3, he4, species, xa, 1d-11, &
            Rho, logRho, energy, lnT_guess, lnT, result_energy, ierr)
      end subroutine set_lnT_for_energy
      
      
      subroutine set_lnT_for_energy_with_tol( &
            s, k, h1, he3, he4, species, xa, tol, &
            Rho, logRho, energy, lnT_guess, lnT, result_energy, ierr)
         use eos_support, only: solve_eos_given_DE
         use eos_def, only: num_eos_basic_results, i_lnE
         use chem_lib, only: basic_composition_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, h1, he3, he4, species
         real(dp), intent(in) :: xa(species), tol, Rho, logRho, energy, lnT_guess
         real(dp), intent(out) :: lnT, result_energy
         integer, intent(out) :: ierr

         real(dp) :: &
            X, Y, Z, T, logT, res(num_eos_basic_results), &
            d_dlnd(num_eos_basic_results), d_dlnT(num_eos_basic_results), &
            d_dabar(num_eos_basic_results), d_dzbar(num_eos_basic_results), &
            abar, zbar, z53bar, z2bar, ye, mass_correction, sumx, logT_tol, logE_tol
         integer :: j

         include 'formats'

         ierr = 0

         call basic_composition_info(species, s% chem_id, xa(:), X, Y, Z, &
               abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

         logT_tol = tol ! 1d-11
         logE_tol = tol ! 1d-11
         call solve_eos_given_DE( &
            s, k, Z, X, abar, zbar, xa(:), &
            logRho, log10(energy), lnT_guess/ln10, &
            logT_tol, logE_tol, &
            logT, res, d_dlnd, d_dlnT, &
            d_dabar, d_dzbar, &
            ierr)
         lnT = logT*ln10
         
         result_energy = exp(res(i_lnE))

         if (ierr /= 0 .or. is_bad_num(lnT)) then
            ierr = -1
            if (s% stop_for_bad_nums) then
               write(*,2) 'lnT', k, lnT
               stop 'mesh adjust: do1_lnT'
            end if
            return
         end if

         contains

         subroutine show
            include 'formats'
            write(*,*)
            write(*,*) 'set_lnT_for_energy ierr', ierr
            write(*,*) 'k', k
            write(*,*)
            write(*,1) 'lnT =', lnT
            write(*,*)
            write(*,1) 'logRho =', logRho
            write(*,1) 'logT_guess =', lnT_guess/ln10
            write(*,1) 'energy =', energy
            write(*,1) 'Z =', Z
            write(*,1) 'X =', X
            write(*,1) 'abar =', abar
            write(*,1) 'zbar =', zbar
            write(*,*)
            write(*,*)
         end subroutine show

      end subroutine set_lnT_for_energy_with_tol


      ! Stiriba, Youssef, Appl, Numer. Math. 45, 499-511. 2003.

         ! LPP-HARMOD -- local piecewise parabolic reconstruction

         ! interpolant is derived to conserve integral of v in cell k
         ! interpolant slope at cell midpoint is harmonic mean of slopes between adjacent cells
         ! where these slopes between cells are defined as the difference between cell averages
         ! divided by the distance between cell midpoints.
         ! interpolant curvature based on difference between the midpoint slope
         ! and the smaller in magnitude of the slopes between adjacent cells.

         ! interpolant f(dq) = a + b*dq + (c/2)*dq^2, with dq = q - q_midpoint
         ! c0 holds a's, c1 holds b's, and c2 holds c's.


      subroutine get1_lpp(k, ldv, nz, j, dq, v, quad, c0, c1, c2)
         integer, intent(in) :: k, ldv, nz, j
         real(dp), intent(in) :: dq(:) ! (nz)
         real(dp), intent(in) :: v(:,:) ! (ldv,nz)
         logical, intent(in) :: quad
         real(dp), dimension(:) :: c0, c1, c2

         real(dp) :: vbdy1, vbdy2, dqhalf, sm1, s00, sprod
         real(dp), parameter :: rel_curvature_limit = 0.1d0

         logical :: dbg

         include 'formats'

         if (k == 1 .or. k == nz) then
            call set_const
            return
         end if

         dbg = .false.
         !dbg = (k == 30 .and. j == 3) ! .false.

         sm1 = (v(j,k-1) - v(j,k)) / ((dq(k-1) + dq(k))/2)
         s00 = (v(j,k) - v(j,k+1)) / ((dq(k) + dq(k+1))/2)

         sprod = sm1*s00
         if (sprod <= 0) then
            ! at local min or max, so set slope and curvature to 0.
            call set_const
            return
         end if

         if (.not. quad) then
            c0(k) = v(j,k)
            c1(k) = (sm1 + s00)/2 ! use average to smooth abundance transitions
            c2(k) = 0 ! Yan Wang fixed this -- it was left out initially.
         else
            c1(k) = sprod*2/(s00 + sm1) ! harmonic mean slope
            if (abs(sm1) <= abs(s00)) then
               c2(k) = (sm1 - c1(k))/(2*dq(k))
            else
               c2(k) = (c1(k) - s00)/(2*dq(k))
            end if
            c0(k) = v(j,k) - c2(k)*dq(k)*dq(k)/24
         end if

         ! check values at edges for monotonicity
         dqhalf = dq(k)/2
         vbdy1 = c0(k) + c1(k)*dqhalf + c2(k)/2*dqhalf*dqhalf ! value at face(k)
         vbdy2 = c0(k) - c1(k)*dqhalf + c2(k)/2*dqhalf*dqhalf ! value at face(k+1)
         if ((v(j,k-1) - vbdy1)*(vbdy1 - v(j,k)) < 0 .or. &
             (v(j,k) - vbdy2)*(vbdy2 - v(j,k+1)) < 0) then
            if (dbg) then
               write(*,*) 'non-monotonic'
               write(*,2) 'v(j,k-1)', k-1, v(j,k-1)
               write(*,2) 'vbdy1', k, vbdy1
               write(*,2) 'v(j,k)', k, v(j,k)
               write(*,2) 'vbdy2', k, vbdy2
               write(*,2) 'v(j,k+1)', k+1, v(j,k+1)
               write(*,*)
               write(*,2) 'v(j,k-1) - vbdy1', k, v(j,k-1) - vbdy1
               write(*,2) 'vbdy1 - v(j,k+1)', k, vbdy1 - v(j,k+1)
               write(*,*)
               write(*,2) 'v(j,k-1) - vbdy2', k, v(j,k-1) - vbdy2
               write(*,2) 'vbdy2 - v(j,k+1)', k, vbdy2 - v(j,k+1)
               write(*,*)
               write(*,2) 'sm1', k, sm1
               write(*,2) 's00', k, s00
               write(*,*)
               stop 'debug: get1_lpp'
            end if
            c2(k) = 0
            if (abs(sm1) <= abs(s00)) then
               c1(k) = sm1
            else
               c1(k) = s00
            end if
         end if

         contains

         subroutine set_const
            c0(k) = v(j,k)
            c1(k) = 0
            if (quad) c2(k) = 0
         end subroutine set_const

      end subroutine get1_lpp


      subroutine get_xq_integral( &
            k_old_in, nz_old, xq_old, xq_outer, dq, &
            order, c0, c1, c2, integral, dbg, k_old_last, ierr)
         ! integrate val(j,:) from xq_inner to xq_outer, with xq_inner = xq_outer + dq
         integer, intent(in) :: k_old_in, nz_old
         real(dp), intent(in) :: xq_old(:), xq_outer, dq
         integer, intent(in) :: order  ! 0, 1, 2
         real(dp), intent(in), dimension(:) :: c0, c1, c2 ! coefficients
         real(dp), intent(out) :: integral
         logical, intent(in) :: dbg
         integer, intent(out) :: k_old_last, ierr

         integer :: k, k_old
         real(dp) :: a, b, c, old_xq_inner, old_xq_outer, xq_inner, &
            xq_overlap_outer, xq_overlap_inner, dq1, sum_dqs, old_xq_mid, &
            v_overlap_outer, v_overlap_inner, dq_outer, dq_inner, avg_value

         include 'formats'

         if (dbg) write(*,*)

         ierr = 0
         k_old = k_old_in
         ! move starting k_old if necessary
         do
            if (k_old <= 1) exit
            if (xq_old(k_old) <= xq_outer) exit
            k_old = k_old - 1
         end do
         xq_inner = xq_outer + dq
         old_xq_inner = xq_old(k_old)
         sum_dqs = 0
         integral = 0d0

         do k = k_old, nz_old
            if (dq <= sum_dqs) exit
            old_xq_outer = old_xq_inner
            if (k == nz_old) then
               old_xq_inner = 1
            else
               old_xq_inner = xq_old(k+1)
            end if
            xq_overlap_outer = max(xq_outer, old_xq_outer)
            xq_overlap_inner = min(xq_inner, old_xq_inner)

            if (dbg) then
               write(*,2) 'xq_overlap_outer', k, xq_overlap_outer
               write(*,2) 'xq_outer', k, xq_outer
               write(*,2) 'old_xq_outer', k, old_xq_outer
               write(*,2) 'xq_overlap_inner', k, xq_overlap_inner
               write(*,2) 'xq_inner', k, xq_inner
               write(*,2) 'old_xq_inner', k, old_xq_inner
            end if

            if (sum_dqs == 0 .and. xq_overlap_outer == xq_outer .and.  &
               xq_overlap_inner == xq_inner) then
               ! fully contained
               xq_inner = xq_outer + dq
               xq_overlap_inner = xq_inner
               dq1 = dq
            else if (old_xq_inner >= xq_inner) then ! this is the last one
               dq1 = dq - sum_dqs
            else
               dq1 = max(0d0, xq_overlap_inner-xq_overlap_outer)
            end if
            sum_dqs = sum_dqs + dq1
            ! interpolant f(dq) = a + b*dq + (c/2)*dq^2, with dq = q - q_midpoint
            a = c0(k)
            b = c1(k)
            if (order == 2) then
               c = c2(k)
            else
               c = 0
            end if

            if (dq1 == 0 .or. (b==0 .and. c==0)) then
               avg_value = a
            else
               old_xq_mid = (old_xq_outer + old_xq_inner)/2
               dq_outer = old_xq_mid - xq_overlap_outer
               dq_inner = old_xq_mid - xq_overlap_inner

               if (order == 0) then
                  avg_value = a
               else if (order == 1 .or. c == 0) then
                  ! use slope to estimate average value in the region being used
                  if (dbg) write(*,*) 'use slope to estimate average value'
                  v_overlap_outer = a + dq_outer*b
                  v_overlap_inner = a + dq_inner*b
                  avg_value = (v_overlap_outer + v_overlap_inner)/2
               else ! use quadratic reconstruction
                  if (dbg) write(*,*) 'use quadratic reconstruction'
                  avg_value = &
                     a + b*(dq_inner + dq_outer)/2 + &
                        c*(dq_inner*dq_inner + dq_inner*dq_outer + dq_outer*dq_outer)/6
               end if
            end if
            integral = integral + dq1*avg_value
            if (dbg) then
               write(*,2) 'a', k, a
               write(*,2) 'b', k, b
               write(*,2) 'c', k, c
               write(*,2) 'dq1', k, dq1
               write(*,2) 'avg_value', k, avg_value
               write(*,2) 'integral', k, integral
               write(*,*)
            end if
            k_old_last = k
            if (old_xq_inner >= xq_inner) exit ! this is the last one

         end do

         if (dbg) write(*,1) 'integral/dq', integral/dq, integral, dq

      end subroutine get_xq_integral


      subroutine adjust_omega(s, nz, nz_old, comes_from, &
            old_xq, new_xq, old_dq, new_dq, xh, old_j_rot, &
            xout_old, xout_new, old_dqbar, new_dqbar, ierr)
         use alloc
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old
         integer, dimension(:) :: comes_from
         real(dp), dimension(:) :: &
            old_xq, new_xq, old_dq, new_dq, old_j_rot, &
            xout_old, xout_new, old_dqbar, new_dqbar
         real(dp), intent(in) :: xh(:,:)
         integer, intent(out) :: ierr
         integer :: k, op_err, old_k, new_k
         real(dp) :: old_j_tot, new_j_tot
         include 'formats'
         ierr = 0

!$OMP PARALLEL DO PRIVATE(k, op_err) SCHEDULE(dynamic,2)
         do k = 1, nz
            op_err = 0
            call adjust1_omega(s, k, nz, nz_old, comes_from, &
               xout_old, xout_new, old_dqbar, new_dqbar, old_j_rot, xh, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO

      end subroutine adjust_omega


      subroutine adjust1_omega(s, k, nz, nz_old, comes_from, &
            xout_old, xout_new, old_dqbar, new_dqbar, old_j_rot, xh, ierr)
         use hydro_rotation, only: w_div_w_roche_jrot, update1_i_rot_from_xh
         ! set new value for s% omega(k)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nz, nz_old
         integer, dimension(:) :: comes_from
         real(dp), dimension(:), intent(in) :: &
            xout_old, xout_new, old_dqbar, new_dqbar, old_j_rot
         real(dp), intent(in) :: xh(:,:)
         integer, intent(out) :: ierr

         real(dp) :: xq_outer, xq_inner, j_tot, xq0, xq1, new_point_dqbar, dq_sum, dq, r00
         integer :: kk, k_outer, j

         integer, parameter :: k_dbg = -1

         include 'formats'

         ierr = 0
         xq_outer = xout_new(k)
         new_point_dqbar = new_dqbar(k)
         if (k < nz) then
            xq_inner = xq_outer + new_point_dqbar
         else
            xq_inner = 1d0
         end if

         if (k == k_dbg) then
            write(*,2) 'xq_outer', k, xq_outer
            write(*,2) 'xq_inner', k, xq_inner
            write(*,2) 'new_point_dqbar', k, new_point_dqbar
         end if

         dq_sum = 0d0
         j_tot = 0
         if (xq_outer >= xout_old(nz_old)) then
            ! new contained entirely in old center zone
            k_outer = nz_old
            if (k == k_dbg) &
               write(*,2) 'new contained in old center', &
                  k_outer, xout_old(k_outer)
         else if (k == 1) then
            k_outer = 1
         else
            k_outer = comes_from(k-1)
         end if

         do kk = k_outer, nz_old ! loop until reach m_inner

            if (kk == nz_old) then
               xq1 = 1d0
            else
               xq1 = xout_old(kk+1)
            end if
            if (xq1 <= xq_outer) cycle

            xq0 = xout_old(kk)
            if (xq0 >= xq_inner) then
               if (dq_sum < new_point_dqbar .and. kk > 1) then
                  ! need to add a bit more from the previous source
                  dq = new_point_dqbar - dq_sum
                  dq_sum = new_point_dqbar
                  j_tot = j_tot + old_j_rot(kk-1)*dq
                  end if
               exit
            end if

            if (xq1 < xq_outer) then
               ierr = -1
               return
            end if

            if (xq0 >= xq_outer .and. xq1 <= xq_inner) then ! entire old kk is in new k

               dq = old_dqbar(kk)
               dq_sum = dq_sum + dq

               if (dq_sum > new_point_dqbar) then
                  ! dq too large -- numerical roundoff problems
                  dq = dq - (new_point_dqbar - dq_sum)
                  dq_sum = new_point_dqbar
               end if

               j_tot = j_tot + old_j_rot(kk)*dq

            else if (xq0 <= xq_outer .and. xq1 >= xq_inner) then ! entire new k is in old kk

               dq = new_dqbar(k)
               dq_sum = dq_sum + dq
               j_tot = j_tot + old_j_rot(kk)*dq

            else ! only use the part of old kk that is in new k

               if (k == k_dbg) then
                  write(*,*) 'only use the part of old kk that is in new k', xq_inner <= xq1
                  write(*,1) 'xq_outer', xq_outer
                  write(*,1) 'xq_inner', xq_inner
                  write(*,1) 'xq0', xq0
                  write(*,1) 'xq1', xq1
                  write(*,1) 'dq_sum', dq_sum
                  write(*,1) 'new_point_dqbar', new_point_dqbar
                  write(*,1) 'new_point_dqbar - dq_sum', new_point_dqbar - dq_sum
               end if

               if (xq_inner <= xq1) then ! this is the last part of new k

                  if (k == k_dbg) write(*,3) 'this is the last part of new k', k, kk

                  dq = new_point_dqbar - dq_sum
                  dq_sum = new_point_dqbar

               else ! we avoid this case if possible because of numerical roundoff

                  if (k == k_dbg) write(*,3) 'we avoid this case if possible', k, kk

                  dq = max(0d0, xq1 - xq_outer)
                  if (dq_sum + dq > new_point_dqbar) dq = new_point_dqbar - dq_sum
                  dq_sum = dq_sum + dq

               end if

               j_tot = j_tot + old_j_rot(kk)*dq

               if (dq <= 0) then
                  ierr = -1
                  return
               end if

            end if

            if (dq_sum >= new_point_dqbar) then
               if (k == k_dbg) then
                  write(*,2) 'exit for k', k
                  write(*,2) 'dq_sum', kk, dq_sum
                  write(*,2) 'new_point_dqbar', kk, new_point_dqbar
               end if
               exit
            end if

         end do

         s% j_rot(k) = j_tot/dq_sum
         r00 = exp(xh(s% i_lnR, k))
         if (s% fitted_fp_ft_i_rot) then
            s% w_div_w_crit_roche(k) = &
               w_div_w_roche_jrot(r00,s% m(k),s% j_rot(k),s% cgrav(k), &
               s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
         end if
         call update1_i_rot_from_xh(s, k)
         s% omega(k) = s% j_rot(k)/s% i_rot(k)

         if (k_dbg == k) then
            write(*,2) 's% omega(k)', k, s% omega(k)
            write(*,2) 's% j_rot(k)', k, s% j_rot(k)
            write(*,2) 's% i_rot(k)', k, s% i_rot(k)
            if (s% model_number > 1925) stop 'debugging: adjust1_omega'
         end if

      end subroutine adjust1_omega


      ! like adjust_omega.  conserve kinetic energy
      subroutine do_v( &
            s, nz, nz_old, cell_type, comes_from, &
            old_xq, new_xq, old_dq, new_dq, xh, xh_old, &
            xout_old, xout_new, old_dqbar, new_dqbar, old_ke, ierr)
         use alloc
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old
         integer, dimension(:) :: cell_type, comes_from
         real(dp), dimension(:) :: &
            xout_old, xout_new, old_dqbar, new_dqbar, &
            old_xq, new_xq, old_dq, new_dq, old_ke
         real(dp), dimension(:,:) :: xh, xh_old
         integer, intent(out) :: ierr

         integer :: k, j, op_err, old_k, new_k, i_v
         real(dp) :: old_ke_tot, new_ke_tot, xmstar, err

         include 'formats'
         ierr = 0
         i_v = s% i_v
         xmstar = s% xmstar

         old_ke_tot = 0d0
         do k=1,nz_old ! skip common factor 1/2 xmstar in ke
            old_ke(k) = old_dqbar(k)*xh_old(i_v,k)*xh_old(i_v,k)
            old_ke_tot = old_ke_tot + old_ke(k)
         end do

!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k = 1, nz
            op_err = 0
            call adjust1_v( &
               s, k, nz, nz_old, cell_type, comes_from, xout_old, xout_new, &
               old_dqbar, new_dqbar, old_ke, i_v, xh, xh_old, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
         if (ierr /= 0) then
            return
         end if

         new_ke_tot = 0
         do k=1,nz
            new_ke_tot = new_ke_tot + new_dqbar(k)*xh(i_v,k)*xh(i_v,k)
         end do

         if (abs(old_ke_tot - new_ke_tot) > 1d-7*new_ke_tot) then
            ierr = -1
            s% retry_message = 'failed in mesh_adjust do_v'
            if (s% report_ierr) write(*, *) s% retry_message
            !stop 'do_v'
            return
         end if

         err = abs(old_ke_tot - new_ke_tot)/max(new_ke_tot,old_ke_tot,1d0)
         s% mesh_adjust_KE_conservation = err

         if (s% trace_mesh_adjust_error_in_conservation) then
            write(*,2) 'mesh adjust error in conservation of KE', s% model_number, &
               err, new_ke_tot, old_ke_tot
            if (err > 1d-10) then
               write(*,*) 'err too large'
               stop 'do_v'
            end if
         end if

      end subroutine do_v


      subroutine adjust1_v( &
            s, k, nz, nz_old, cell_type, comes_from, xout_old, xout_new, &
            old_dqbar, new_dqbar, old_ke, i_v, xh, xh_old, ierr)
         ! set new value for s% v(k) to conserve kinetic energy
         type (star_info), pointer :: s
         integer, intent(in) :: k, nz, nz_old, i_v
         integer, dimension(:) :: cell_type, comes_from
         real(dp), dimension(:), intent(in) :: &
            xout_old, xout_new, old_dqbar, new_dqbar, old_ke
         real(dp), dimension(:,:) :: xh, xh_old
         integer, intent(out) :: ierr

         real(dp) :: xq_outer, xq_inner, ke_sum, &
            xq0, xq1, new_point_dqbar, dq_sum, dq
         integer :: kk, k_outer, j

         integer, parameter :: k_dbg = -1

         include 'formats'

         ierr = 0

         if (cell_type(k) == unchanged_type .or. &
               cell_type(k) == revised_type) then
            if (k == 1) then
               xh(i_v,k) = xh_old(i_v,comes_from(k))
               return
            end if
            if (cell_type(k-1) == unchanged_type) then
               xh(i_v,k) = xh_old(i_v,comes_from(k))
               return
            end if
         end if

         xq_outer = xout_new(k)
         new_point_dqbar = new_dqbar(k)
         if (k < nz) then
            xq_inner = xq_outer + new_point_dqbar
         else
            xq_inner = 1d0
         end if

         if (k == k_dbg) then
            write(*,2) 'xq_outer', k, xq_outer
            write(*,2) 'xq_inner', k, xq_inner
            write(*,2) 'new_point_dqbar', k, new_point_dqbar
         end if

         dq_sum = 0d0
         ke_sum = 0
         if (xq_outer >= xout_old(nz_old)) then
            ! new contained entirely in old center zone
            k_outer = nz_old
            if (k == k_dbg) &
               write(*,2) 'new contained in old center', &
                  k_outer, xout_old(k_outer)
         else if (k == 1) then
            k_outer = 1
         else
            k_outer = comes_from(k-1)
         end if

         do kk = k_outer, nz_old ! loop until reach xq_inner

            if (kk == nz_old) then
               xq1 = 1d0
            else
               xq1 = xout_old(kk+1)
            end if
            if (xq1 <= xq_outer) cycle

            xq0 = xout_old(kk)
            if (xq0 >= xq_inner) then
               if (dq_sum < new_point_dqbar .and. kk > 1) then
                  ! need to add a bit more from the previous source
                  dq = new_point_dqbar - dq_sum
                  dq_sum = new_point_dqbar
                  ke_sum = ke_sum + old_ke(kk-1)*dq/old_dqbar(kk-1)

                  if (k == k_dbg) &
                     write(*,3) 'new k contains some of old kk-1', &
                        k, kk, old_ke(kk-1)*dq, dq_sum

                  end if
               exit
            end if

            if (xq1 < xq_outer) then
               ierr = -1
               return
            end if

            if (xq0 >= xq_outer .and. xq1 <= xq_inner) then ! entire old kk is in new k

               dq = old_dqbar(kk)
               dq_sum = dq_sum + dq

               if (dq_sum > new_point_dqbar) then
                  ! dq too large -- numerical roundoff problems
                  dq = dq - (new_point_dqbar - dq_sum)
                  dq_sum = new_point_dqbar
               end if

               ke_sum = ke_sum + old_ke(kk)*dq/old_dqbar(kk)

               if (k == k_dbg) &
                  write(*,3) 'new k contains all of old kk', &
                     k, kk, old_ke(kk)*dq, ke_sum

            else if (xq0 <= xq_outer .and. xq1 >= xq_inner) then ! entire new k is in old kk

               dq = new_dqbar(k)
               dq_sum = dq_sum + dq
               ke_sum = ke_sum + old_ke(kk)*dq/old_dqbar(kk)

               if (k == k_dbg) &
                  write(*,3) 'all new k is in old kk', &
                     k, kk, old_ke(kk)*dq, ke_sum

            else ! only use the part of old kk that is in new k

               if (k == k_dbg) then
                  write(*,*) 'only use the part of old kk that is in new k', xq_inner <= xq1
                  write(*,1) 'xq_outer', xq_outer
                  write(*,1) 'xq_inner', xq_inner
                  write(*,1) 'xq0', xq0
                  write(*,1) 'xq1', xq1
                  write(*,1) 'dq_sum', dq_sum
                  write(*,1) 'new_point_dqbar', new_point_dqbar
                  write(*,1) 'new_point_dqbar - dq_sum', new_point_dqbar - dq_sum
               end if

               if (xq_inner <= xq1) then ! this is the last part of new k

                  dq = new_point_dqbar - dq_sum
                  dq_sum = new_point_dqbar

               else ! we avoid this case if possible because of numerical roundoff

                  if (k == k_dbg) write(*,3) 'we avoid this case if possible', k, kk

                  dq = max(0d0, xq1 - xq_outer)
                  if (dq_sum + dq > new_point_dqbar) dq = new_point_dqbar - dq_sum
                  dq_sum = dq_sum + dq

               end if

               if (k == k_dbg) then
                  write(*,3) 'new k use only part of old kk', k, kk
                  write(*,2) 'dq_sum', k, dq_sum
                  write(*,2) 'dq', k, dq
                  write(*,2) 'old_ke(kk)', kk, old_ke(kk)
                  write(*,2) 'old ke_sum', k, ke_sum
                  write(*,2) 'new ke_sum', k, ke_sum + old_ke(kk)*dq
               end if

               ke_sum = ke_sum + old_ke(kk)*dq/old_dqbar(kk)

               if (dq <= 0) then
                  ierr = -1
                  !return
                  write(*,*) 'dq <= 0', dq
                  stop 'debugging: adjust1_v'
               end if

            end if

            if (dq_sum >= new_point_dqbar) then
               exit
            end if

         end do

         xh(i_v,k) = sqrt(ke_sum/new_point_dqbar) ! we have skipped the 1/2 xmstar factor
         if (xh_old(i_v,comes_from(k)) < 0d0) xh(i_v,k) = -xh(i_v,k)

         if (k == k_dbg) then
!$OMP critical (adjust1_v_dbg)
            write(*,2) 'xh(i_v,k) new_dqbar', k, xh(i_v,k), new_dqbar(k)
            write(*,2) 'xh_old(i_v,comes_from(k)) old_dqbar', &
               comes_from(k), xh_old(i_v,comes_from(k)), old_dqbar(comes_from(k))
            if (k == k_dbg) stop 'adjust1_v'
!$OMP end critical (adjust1_v_dbg)
            !stop
         end if

      end subroutine adjust1_v


      subroutine do_u( &
            s, nz, nz_old, cell_type, comes_from, &
            old_xq, new_xq, old_dq, new_dq, xh, xh_old, &
            xout_old, xout_new, old_ke, ierr)
         use alloc
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old
         integer, dimension(:) :: cell_type, comes_from
         real(dp), dimension(:) :: &
            xout_old, xout_new, old_xq, new_xq, old_dq, new_dq, old_ke
         real(dp), dimension(:,:) :: xh, xh_old
         integer, intent(out) :: ierr

         integer :: k, j, op_err, old_k, new_k, i_u
         real(dp) :: old_ke_tot, new_ke_tot, xmstar, err

         include 'formats'
         ierr = 0
         i_u = s% i_u
         xmstar = s% xmstar

         old_ke_tot = 0d0
         do k=1,nz_old ! skip common factor 1/2 xmstar in ke
            old_ke(k) = old_dq(k)*xh_old(i_u,k)*xh_old(i_u,k)
            old_ke_tot = old_ke_tot + old_ke(k)
         end do

!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k = 1, nz
            op_err = 0
            call adjust1_u( &
               s, k, nz, nz_old, cell_type, comes_from, xout_old, xout_new, &
               old_dq, new_dq, old_ke, i_u, xh, xh_old, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
         if (ierr /= 0) then
            return
         end if

         new_ke_tot = 0
         do k=1,nz
            new_ke_tot = new_ke_tot + new_dq(k)*xh(i_u,k)*xh(i_u,k)
         end do

         err = abs(old_ke_tot - new_ke_tot)/max(new_ke_tot,old_ke_tot,1d0)
         s% mesh_adjust_KE_conservation = err

         if (s% trace_mesh_adjust_error_in_conservation) then
            write(*,2) 'mesh adjust error in conservation of KE', s% model_number, &
               err, new_ke_tot, old_ke_tot
            if (err > 1d-10) then
               write(*,*) 'err too large'
               stop 'do_u'
            end if
         end if

      end subroutine do_u


      subroutine adjust1_u( &
            s, k, nz, nz_old, cell_type, comes_from, xout_old, xout_new, &
            old_dq, new_dq, old_ke, i_u, xh, xh_old, ierr)
         ! set new value for s% u(k) to conserve kinetic energy
         type (star_info), pointer :: s
         integer, intent(in) :: k, nz, nz_old, i_u
         integer, dimension(:) :: cell_type, comes_from
         real(dp), dimension(:), intent(in) :: &
            xout_old, xout_new, old_dq, new_dq, old_ke
         real(dp), dimension(:,:) :: xh, xh_old
         integer, intent(out) :: ierr

         real(dp) :: xq_outer, xq_inner, ke_sum, &
            xq0, xq1, new_cell_dq, dq_sum, dq
         integer :: kk, k_outer, j

         integer, parameter :: k_dbg = -1

         include 'formats'

         ierr = 0

         if (cell_type(k) == unchanged_type .or. &
               cell_type(k) == revised_type) then
            if (k == 1) then
               xh(i_u,k) = xh_old(i_u,comes_from(k))
               return
            end if
            if (cell_type(k-1) == unchanged_type) then
               xh(i_u,k) = xh_old(i_u,comes_from(k))
               return
            end if
         end if

         xq_outer = xout_new(k)
         new_cell_dq = new_dq(k)
         if (k < nz) then
            xq_inner = xq_outer + new_cell_dq
         else
            xq_inner = 1d0
         end if

         if (k == k_dbg) then
            write(*,2) 'xq_outer', k, xq_outer
            write(*,2) 'xq_inner', k, xq_inner
            write(*,2) 'new_cell_dq', k, new_cell_dq
         end if

         dq_sum = 0d0
         ke_sum = 0
         if (xq_outer >= xout_old(nz_old)) then
            ! new contained entirely in old center zone
            k_outer = nz_old
            if (k == k_dbg) &
               write(*,2) 'new contained in old center', &
                  k_outer, xout_old(k_outer)
         else if (k == 1) then
            k_outer = 1
         else
            k_outer = comes_from(k-1)
         end if

         do kk = k_outer, nz_old ! loop until reach xq_inner

            if (kk == nz_old) then
               xq1 = 1d0
            else
               xq1 = xout_old(kk+1)
            end if
            if (xq1 <= xq_outer) cycle

            if (xq1 < xq_outer) then
               ierr = -1
               return
            end if

            xq0 = xout_old(kk)
            if (xq0 >= xq_outer .and. xq1 <= xq_inner) then ! entire old kk is in new k

               dq = old_dq(kk)
               dq_sum = dq_sum + dq

               if (dq_sum > new_cell_dq) then
                  ! dq too large -- numerical roundoff problems
                  dq = dq - (new_cell_dq - dq_sum)
                  dq_sum = new_cell_dq
               end if

               ke_sum = ke_sum + old_ke(kk)*dq/old_dq(kk)

               if (k == k_dbg) &
                  write(*,3) 'new k contains all of old kk', &
                     k, kk, old_ke(kk)*dq, ke_sum

            else if (xq0 <= xq_outer .and. xq1 >= xq_inner) then ! entire new k is in old kk

               dq = new_dq(k)
               dq_sum = dq_sum + dq
               ke_sum = ke_sum + old_ke(kk)*dq/old_dq(kk)

               if (k == k_dbg) &
                  write(*,3) 'all new k is in old kk', &
                     k, kk, old_ke(kk)*dq, ke_sum

            else ! only use the part of old kk that is in new k

               if (k == k_dbg) then
                  write(*,*) 'only use the part of old kk that is in new k', xq_inner <= xq1
                  write(*,1) 'xq_outer', xq_outer
                  write(*,1) 'xq_inner', xq_inner
                  write(*,1) 'xq0', xq0
                  write(*,1) 'xq1', xq1
                  write(*,1) 'dq_sum', dq_sum
                  write(*,1) 'new_cell_dq', new_cell_dq
                  write(*,1) 'new_cell_dq - dq_sum', new_cell_dq - dq_sum
               end if

               if (xq_inner <= xq1) then ! this is the last part of new k

                  dq = new_cell_dq - dq_sum
                  dq_sum = new_cell_dq

               else ! we avoid this case if possible because of numerical roundoff

                  if (k == k_dbg) write(*,3) 'we avoid this case if possible', k, kk

                  dq = max(0d0, xq1 - xq_outer)
                  if (dq_sum + dq > new_cell_dq) dq = new_cell_dq - dq_sum
                  dq_sum = dq_sum + dq

               end if

               if (k == k_dbg) then
                  write(*,3) 'new k use only part of old kk', k, kk
                  write(*,2) 'dq_sum', k, dq_sum
                  write(*,2) 'dq', k, dq
                  write(*,2) 'old_ke(kk)', kk, old_ke(kk)
                  write(*,2) 'old ke_sum', k, ke_sum
                  write(*,2) 'new ke_sum', k, ke_sum + old_ke(kk)*dq
               end if

               ke_sum = ke_sum + old_ke(kk)*dq/old_dq(kk)

               if (dq <= 0) then
                  ierr = -1
                  !return
                  write(*,*) 'dq <= 0', dq
                  stop 'debugging: adjust1_u'
               end if

            end if

            if (dq_sum >= new_cell_dq) then
               exit
            end if

         end do

         xh(i_u,k) = sqrt(ke_sum/new_cell_dq) ! we have skipped the 1/2 xmstar factor
         if (xh_old(i_u,comes_from(k)) < 0d0) xh(i_u,k) = -xh(i_u,k)

         if (k == k_dbg) then
!$OMP critical (adjust1_u_dbg)
            write(*,2) 'xh(i_u,k) new_dq', k, xh(i_u,k), new_dq(k)
            write(*,2) 'xh_old(i_u,comes_from(k)) old_dq', &
               comes_from(k), xh_old(i_u,comes_from(k)), old_dq(comes_from(k))
            if (k == k_dbg) stop 'adjust1_u'
!$OMP end critical (adjust1_u_dbg)
            !stop
         end if

      end subroutine adjust1_u



      end module mesh_adjust

