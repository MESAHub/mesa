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

      module adjust_mesh_plot

      use star_private_def
      use const_def
      use adjust_mesh_support

      implicit none

      contains


      subroutine write_plot_data_for_new_mesh( &
            s, nz, nz_old, xh_old, xa_old, &
            D_mix, q_old, &
            xh, xa, dq, q, xq, species, net_iso, &
            num_gvals, gval_names, gvals, &
            which_gval, comes_from, cell_type, &
            delta_gval_max, &
            xmstar, ierr)
         use chem_def
         use interp_1d_def
         use interp_1d_lib
         use utils_lib

         type (star_info), pointer :: s
         integer, intent(in) :: nz, nz_old, species
         integer, dimension(:), pointer :: net_iso
         real(dp), intent(in) :: xmstar
         real(dp), dimension(:), pointer :: &
            q_old, dq, q, xq, D_mix
         real(dp), dimension(:), pointer :: delta_gval_max
         real(dp), dimension(:,:), pointer :: xh_old, xa_old
         real(dp), dimension(:,:), pointer :: xh, xa

         integer, intent(in) :: num_gvals
         character (len=32) :: gval_names(:) ! (num_gvals)  for debugging.
         real(dp), pointer :: gvals(:,:) ! (nz_old, num_gvals)
         integer, pointer :: which_gval(:), comes_from(:), cell_type(:)

         integer, intent(out) :: ierr

         real(dp) :: v, u, ratio, lum, min_ratio, mstar, lgRho, lgT
         real(dp), pointer :: v_old(:), v_new(:,:), work(:)
         integer :: k, nwork, num_vals, iounit, k_min_ratio, d_comes_from_dk, &
            jD_mix, jgval_factor, jgval_max1, &
            jlogR, jFL, jv, ju, j, &
            i_lum, i_lnd, i_lnT, i_v, i_u, i_lnR
         character (len=100) :: filename, name

         ierr = 0

         i_lnd = s% i_lnd
         i_v = s% i_v
         i_u = s% i_u
         i_lum = s% i_lum
         i_lnT = s% i_lnT
         i_lnR = s% i_lnR

         mstar = xmstar + s% M_center

         call mkdir('mesh_plot_data')
         filename = 'mesh_plot_data/new_mesh.data'
         open(newunit=iounit, file=trim(filename), action='write', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            return
         end if

         write(*,*) 'write ' // trim(filename)

         nwork = mp_work_size
         j = 0
         j = j+1; jgval_max1 = j
         j = j+1; jD_mix = j
         j = j+1; jlogR = j
         j = j+1; jFL = j
         j = j+1; jv = j
         j = j+1; ju = j
         j = j+1; jgval_factor = j  ! this one must be last

         num_vals = j + num_gvals

         allocate(v_old(nz_old), v_new(nz, num_vals), work(nz_old*nwork), stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'allocate failed for write new mesh'
            return
         end if

         ! interpolate from old values to new points for comparison

         do k=1,nz_old
            v_old(k) = safe_log10(D_mix(k))
         end do
         call interpolate_vector( &
               nz_old, q_old, nz, q, v_old, v_new(:,jD_mix), interp_m3q, nwork, work, &
               'write_plot_data_for_new_mesh', ierr)
         v_new(:,jD_mix) = max(0d0, v_new(:,jD_mix))

         do k=1,nz_old
            v_old(k) = delta_gval_max(k)
         end do
         call interpolate_vector( &
               nz_old, q_old, nz, q, v_old, v_new(:,jgval_max1), interp_m3q, nwork, work, &
               'write_plot_data_for_new_mesh', ierr)

         do j=1,num_gvals
            do k=1,nz_old
            v_old(k) = gvals(k,j)
            end do
            call interpolate_vector( &
                  nz_old, q_old, nz, q, v_old, v_new(:,jgval_factor+j), interp_m3q, nwork, work, &
                  'write_plot_data_for_new_mesh', ierr)
         end do

         do k=1,nz_old
            v_old(k) = xh_old(i_lnR,k)/ln10
         end do
         call interpolate_vector( &
               nz_old, q_old, nz, q, v_old, v_new(:,jlogR), interp_m3q, nwork, work, &
               'write_plot_data_for_new_mesh', ierr)

         if (i_lum == 0) then
            v_new(:,jFL) = 0
         else
            do k=1,nz_old
               v_old(k) = xh_old(i_lum,k)
            end do
            call interpolate_vector( &
                  nz_old, q_old, nz, q, v_old, v_new(:,jFL), interp_m3q, nwork, work, &
               'write_plot_data_for_new_mesh', ierr)
         end if

         if (i_v == 0) then
            v_new(:,jv) = 0
         else
            do k=1,nz_old
            v_old(k) = xh_old(i_v,k)
            end do
            call interpolate_vector( &
                  nz_old, q_old, nz, q, v_old, v_new(:,jv), interp_m3q, nwork, work, &
               'write_plot_data_for_new_mesh', ierr)
         end if

         if (i_u == 0) then
            v_new(:,ju) = 0
         else
            do k=1,nz_old
            v_old(k) = xh_old(i_u,k)
            end do
            call interpolate_vector( &
                  nz_old, q_old, nz, q, v_old, v_new(:,ju), interp_m3q, nwork, work, &
               'write_plot_data_for_new_mesh', ierr)
         end if

         ierr = 0

         write(iounit, fmt='(a6, 99(a27, 1x))', advance='no') &
            'k', 'q', 'log_dq', 'ratio', 'gval_max', &
            'which_gval', 'comes_from', 'd_comes_from_dk', 'cell_type'
         do j=1,num_gvals
            if (j < 10) then
               write(name,'(a,i1,a)') 'gval_', j, '_' // trim(gval_names(j))
            else
               write(name,'(a,i2,a)') 'gval_', j, '_' // trim(gval_names(j))
            end if
            write(iounit, fmt='(a27, 1x)', advance='no') trim(name)
         end do
         write(iounit, fmt='(99(a27, 1x))', advance='no') &
            'log_D', &
            'mass', 'xq', 'radius', 'logR', 'logT', 'logRho', 'logPgas', &
            'lum', 'v'
         do j=1,num_gvals
            if (j < 10) then
               write(name,'(a,i1,a)') 'gval_', j, '_' // trim(gval_names(j))
            else
               write(name,'(a,i2,a)') 'gval_', j, '_' // trim(gval_names(j))
            end if
            write(iounit, fmt='(a27, 1x)', advance='no') 'd_dk_' // trim(name)
         end do
         write(iounit,*)

         min_ratio = 1d99; k_min_ratio = -1
         do k=1,nz
            if (i_v == 0) then
               v = 0
            else
               v = xh(i_v,k)
            end if
            if (i_lum == 0) then
               lum = 0
            else
               lum = xh(i_lum,k)
            end if
            if (i_lnd /= 0) then
               lgRho = xh(i_lnd,k)/ln10
            else
               lgRho = 0
            end if
            if (i_lnT /= 0) then
               lgT = xh(i_lnT,k)/ln10
            else
               lgT = 0
            end if
            if (k == 1) then
               ratio = 1
            else if (dq(k-1) > 0) then
               ratio = dq(k) / dq(k-1)
            else
               write(*,*) 'k-1', k-1
               write(*,*) 'dq(k-1)', dq(k-1)
               ierr = -1
               return
            end if
            if (ratio < min_ratio) then
               min_ratio = ratio; k_min_ratio = k
            end if
            if (k == 1) then
               d_comes_from_dk = 1
            else
               d_comes_from_dk = comes_from(k) - comes_from(k-1)
            end if
            write(iounit, fmt='(i6, 99(1pes27.16e3,1x))', advance='no') &
               k, q(k), safe_log10(dq(k)), ratio, v_new(k,jgval_max1), dble(which_gval(k)), &
               dble(comes_from(k)), dble(d_comes_from_dk), dble(cell_type(k))
            do j=1,num_gvals
               write(iounit, fmt='(99(1pes27.16e3,1x))', advance='no') v_new(k,jgval_factor+j)
            end do
            write(iounit, fmt='(99(1pes27.16e3,1x))', advance='no') &
               v_new(k,jD_mix), &
               s% m(k)/Msun, xq(k), exp(xh(i_lnR,k))/Rsun, &
               xh(i_lnR,k)/ln10, lgT, lgRho, lum/Lsun, v
            do j=1,num_gvals
               if (k == 1) then
                  write(iounit, fmt='(99(1pes27.16e3,1x))', advance='no') 0d0
               else
                  write(iounit, fmt='(99(1pes27.16e3,1x))', advance='no') &
                     v_new(k,jgval_factor+j) - v_new(k-1,jgval_factor+j)
               end if
            end do
            write(iounit,*)
         end do

         close(iounit)

         contains

         subroutine get_xa(cid)
            integer, intent(in) :: cid
            integer :: k
            if (cid /= 0) then
               do k=1,nz_old
            v_old(k) = xa_old(cid,k)
               end do
            else
               v_old(1:nz_old) = 0
            end if
         end subroutine get_xa

      end subroutine write_plot_data_for_new_mesh


      subroutine write_plot_data_for_mesh_plan( &
            s, nz_old, nz_new, xh_old, xa_old, &
            lnd_old, lnT_old, lnPgas_old, lnE_old, wturb_old, &
            D_mix, mixing_type, &
            dq_old, q_old, xq_old, q_new, &
            species, i_lnR, i_lum, i_v, i_u, comes_from, &
            num_gvals, gval_names, gvals, delta_gval_max, &
            xmstar, ierr)
         use chem_def
         use eos_lib, only: Radiation_Pressure
         use interp_1d_def
         use interp_1d_lib
         use utils_lib

         type (star_info), pointer :: s
         integer, intent(in) :: nz_old, nz_new, species, i_lnR, i_lum, i_v, i_u
         integer, dimension(:), pointer :: mixing_type, comes_from
         real(dp), intent(in) :: xmstar
         real(dp), dimension(:), pointer :: &
            dq_old, q_old, xq_old, q_new, &
            lnd_old, lnT_old, lnPgas_old, lnE_old, wturb_old, &
            D_mix
         real(dp), dimension(:, :), pointer :: xh_old, xa_old

         integer, intent(in) :: num_gvals
         character (len=32) :: gval_names(:) ! (num_gvals)  for debugging.
         real(dp), pointer :: gvals(:,:) ! (nz_old, num_gvals)
         real(dp), pointer :: delta_gval_max(:) ! (nz_old)

         integer, intent(out) :: ierr

         real(dp) :: v, u, lum, mstar, wturb
         integer :: k, iounit, j
         character (len=100) :: filename, name
         
         include 'formats'

         ierr = 0
         mstar = xmstar + s% M_center

         call mkdir('mesh_plot_data')
         filename = 'mesh_plot_data/mesh_plan.data'
         open(newunit=iounit, file=trim(filename), action='write', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            return
         end if

         write(*,*) 'write ' // trim(filename)
         ierr = 0

         write(iounit, fmt='(99(a27, 1x))', advance='no') &
            'k', 'q_old', 'log_dq_old', 'new_in_old', 'log_delta_gval_max'
         do j=1,num_gvals
            if (j < 10) then
               write(name,'(a,i1,a)') 'gval_', j, '_' // trim(gval_names(j))
            else
               write(name,'(a,i2,a)') 'gval_', j, '_' // trim(gval_names(j))
            end if
            write(iounit, fmt='(a27, 1x)', advance='no') trim(name)
         end do
         do j=1,num_gvals
            if (j < 10) then
               write(name,'(a,i1,a)') 'gval_', j, '_' // trim(gval_names(j))
            else
               write(name,'(a,i2,a)') 'gval_', j, '_' // trim(gval_names(j))
            end if
            write(iounit, fmt='(a27, 1x)', advance='no') 'd_dk_' // trim(name)
         end do
         write(iounit, fmt='(99(a27, 1x))', advance='no') &
            'log_D', 'mixing_type', 'mass', &
            'xq', 'radius', 'logR', 'logRho', 'logT', 'logP', 'logPgas', &
            'logE', 'wturb', 'lum', 'v', 'u'
         write(iounit,*)
         wturb = 0
         do k=1,nz_old
            if (i_v == 0) then
               v = 0
            else
               v = xh_old(i_v,k)
            end if
            if (i_u == 0) then
               u = 0
            else
               u = xh_old(i_u,k)
            end if
            if (i_lum == 0) then
               lum = 0
            else
               lum = xh_old(i_lum,k)
            end if
            if (s% w_flag) wturb = wturb_old(k)
            write(iounit, fmt='(i27, 1x, 99(1pes27.16e3,1x))', advance='no') &
               k, q_old(k), safe_log10(dq_old(k)), new_in_old(k), &
               safe_log10(delta_gval_max(k))
            do j=1,num_gvals
               write(iounit, fmt='(99(1pes27.16e3,1x))', advance='no') gvals(k,j)
            end do
            do j=1,num_gvals
               if (k == nz_old) then
                  write(iounit, fmt='(99(1pes27.16e3,1x))', advance='no') 0d0
               else
                  write(iounit, fmt='(99(1pes27.16e3,1x))', advance='no') gvals(k+1,j) - gvals(k,j)
               end if
            end do
            write(iounit, fmt='(99(1pes27.16e3,1x))', advance='no') &
               log10(max(1d0,D_mix(k))), dble(mixing_type(k)), &
               (s% M_center + xmstar*q_old(k))/Msun, xq_old(k), exp(xh_old(i_lnR,k))/Rsun, &
               xh_old(i_lnR,k)/ln10, lnd_old(k)/ln10, lnT_old(k)/ln10, &
               log10(exp(lnPgas_old(k)) + Radiation_Pressure(exp(lnT_old(k)))), &
               lnPgas_old(k)/ln10, lnE_old(k)/ln10, wturb, lum/Lsun, v, u
            write(iounit,*)
         end do

         close(iounit)

         contains

         real(dp) function new_in_old(k_old)
            integer, intent(in) :: k_old
            integer :: k, cnt
            cnt = 0
            do k = 1, nz_new
               if (comes_from(k) == k_old) cnt = cnt+1
            end do
            new_in_old = dble(cnt)
         end function new_in_old


      end subroutine write_plot_data_for_mesh_plan


      end module adjust_mesh_plot


