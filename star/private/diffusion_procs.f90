! ***********************************************************************
!
!   Copyright (C) 2012-2019  The MESA Team
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
!
! ***********************************************************************

      module diffusion_procs

      use const_def
      use chem_def
      use diffusion_support

      implicit none


      contains


      subroutine fixup( &
            s, nz, nc, m, nzlo, nzhi, total_num_iters, &
            min_X_hard_limit, X_total_atol, X_total_rtol, &
            cell_dm, mtotal, xtotal_init, X, &
            lnT, sum_mass, mass, dX_dm, &
            bad_j, bad_k, bad_X, bad_sum, bad_Xsum, ierr)

         type (star_info), pointer :: s
         integer, intent(in) :: &
            nz, nc, m, nzlo, nzhi, total_num_iters
         real(dp), intent(in) :: mtotal, min_X_hard_limit, &
            X_total_atol, X_total_rtol
         real(dp), intent(in), dimension(:) :: cell_dm, xtotal_init, lnT
         real(dp), intent(inout), dimension(:,:) :: X
         real(dp), intent(inout), dimension(:) :: sum_mass
         real(dp), intent(inout), dimension(:,:) :: mass, dX_dm
         real(dp), intent(out) :: bad_X, bad_sum, bad_Xsum
         integer, intent(out) :: bad_j, bad_k, ierr

         integer :: i, j, jj, k, k1, op_err, rep, num_extreme, &
            j_max, k_max, k_lo, k_hi, kk, cnt, maxcnt, nsmooth_x_in_fixup
         real(dp) :: max_lnT_for_smooth, sum_m, m_00, dm, source_mass, frac, err, &
            max_abs_dx, xm1, x00, xp1, m0, sum0, m1, sum1, x_new, xtotal_new, &
            sum_m1, sum_00, sum_p1, m_m1, m_p1, xavg, dm_m1, dm_p1, x1

         logical :: dbg

         include 'formats'

         ierr = 0
         nsmooth_x_in_fixup = 0
         max_lnT_for_smooth = ln10*99 ! max_logT_for_smooth  -- always on for now.
         bad_j = 0
         bad_k = 0
         bad_X = 0d0
         bad_sum = 0d0
         bad_Xsum = 0d0

         dbg = .false. ! total_num_iters == 21

         ! set mass(j,k) using cell_dm and X
         ! note that can have bad X's, so can have bad masses too.
         do k = nzlo,nzhi
            do j=1,nc
               if (X(j,k) < min_X_hard_limit) then
                  bad_j = j
                  bad_k = k
                  bad_X = X(j,k)
                  if (dbg) write(*,3) 'min_X_hard_limit', j, k, X(j,k), min_X_hard_limit
                  if (dbg) call mesa_error(__FILE__,__LINE__,'fixup')
                  ierr = -1
                  return
               end if
               mass(j,k) = cell_dm(k)*X(j,k)
            end do
         end do

         call fix_negative_masses(nz, nc, nzlo, nzhi, cell_dm, mass, ierr)
         if (ierr /= 0) return

         call fix_species_conservation( &
            nz, nc, nzlo, nzhi, mtotal, xtotal_init, X_total_atol, X_total_rtol, &
            mass, sum_mass, bad_Xsum, bad_j, ierr)
         if (ierr /= 0) return

         call fix_single_point_extremes( &
            nz, nc, nzlo, nzhi, max_lnT_for_smooth, lnT, sum_mass, mass, ierr)
         if (ierr /= 0) return

         if (nsmooth_x_in_fixup > 0) then
            call smooth_x( &
               nz, nc, nzlo, nzhi, nsmooth_x_in_fixup, &
               max_lnT_for_smooth, lnT, sum_mass, mass, ierr)
            if (ierr /= 0) return
         end if

!$OMP PARALLEL DO PRIVATE(k, op_err) SCHEDULE(dynamic,2)

         do k = nzlo, nzhi
            call get1_dX_dm( &
               k, nz, nc, nzlo, nzhi, &
               mass, sum_mass, &
               dX_dm, .false., op_err)
            if (op_err /= 0) then
               bad_k = k
               ierr = op_err
            end if
         end do

!$OMP END PARALLEL DO

         if (ierr /= 0) then
            if (dbg) write(*,2) 'failed in get1_dX_dm', bad_k
            if (dbg) call mesa_error(__FILE__,__LINE__,'fixup')
            return
         end if

         if (dbg) write(*,*) 'call redistribute_mass'
         call redistribute_mass( &
            s, nc, nzlo, nzhi, nz, total_num_iters, &
            sum_mass, mass, dX_dm, &
            xtotal_init, cell_dm, X, .false., ierr)
         if (ierr /= 0) then
            if (dbg) write(*,2) 'failed in redistribute_mass'
            if (dbg) call mesa_error(__FILE__,__LINE__,'fixup')
            return
         end if

         if (dbg) call mesa_error(__FILE__,__LINE__,'fixup')

      end subroutine fixup


      subroutine fix_single_point_extremes( &
            nz, nc, nzlo, nzhi, max_lnT_for_smooth, lnT, sum_mass, mass, ierr)

         integer, intent(in) :: nz, nc, nzlo, nzhi
         real(dp), intent(in) :: max_lnT_for_smooth, lnT(:)
         real(dp), intent(inout) :: sum_mass(:), mass(:,:)
         integer, intent(out) :: ierr

         integer :: k, j, k1
         real(dp) :: xm1, x00, xp1, x1, m0, sum0, m1, sum1, dm

         include 'formats'

         ierr = 0

         ! move mass to remove single point extremes in X
         do k=nzlo+1,nzhi-1
            if (lnT(k) > max_lnT_for_smooth) exit
            do j=1,nc
               xm1 = mass(j,k-1)/sum_mass(k-1)
               x00 = mass(j,k)/sum_mass(k)
               xp1 = mass(j,k+1)/sum_mass(k+1)
               if ((xm1-x00)*(x00-xp1) < -1d-24) then
                  if (abs(xm1-x00) < abs(x00-xp1)) then
                     k1 = k-1; x1 = xm1
                  else
                     k1 = k+1; x1 = xp1
                  end if
                  ! make X(j,k)==X(j,k1) while conserving total j mass
                  m0 = mass(j,k)
                  sum0 = sum_mass(k)
                  m1 = mass(j,k1)
                  sum1 = sum_mass(k1)
                  ! find dm s.t. xnew = (m0+dm)/(sum0+dm) = (m1-dm)/(sum1-dm)
                  dm = (m0*sum1 - m1*sum0)/(m0 + m1 - sum0 - sum1)
                  if (dm > 0) then ! moving mass from k1
                     dm = min(dm,0.99999d0*sum1)
                  else ! moving mass from k
                     dm = -min(-dm,0.99999d0*sum0)
                  end if
                  mass(j,k) = m0 + dm
                  sum_mass(k) = sum(mass(1:nc,k)) ! sum0 + dm
                  if (sum_mass(k) < 0d0) then
                     write(*,2) 'sum_mass(k)', k, sum_mass(k)
                     write(*,1) 'sum0', sum0
                     write(*,1) 'sum1', sum1
                     write(*,1) 'm0', m0
                     write(*,1) 'm1', m1
                     write(*,1) 'dm', dm
                     call mesa_error(__FILE__,__LINE__,'fixup extremes')
                  end if
                  mass(j,k1) = m1 - dm
                  sum_mass(k1) = sum(mass(1:nc,k1)) ! sum1 - dm
                  if (sum_mass(k1) < 0d0) then
                     write(*,2) 'sum_mass(k1)', k1, sum_mass(k1)
                     write(*,1) 'sum0', sum0
                     write(*,1) 'sum1', sum1
                     write(*,1) 'm0', m0
                     write(*,1) 'm1', m1
                     write(*,1) 'dm', dm
                     call mesa_error(__FILE__,__LINE__,'fixup extremes')
                  end if
               end if
            end do
         end do

         do k=nzlo,nzhi
            if (sum_mass(k) < 0d0) then
               write(*,2) 'sum_mass(k)', k, sum_mass(k)
               call mesa_error(__FILE__,__LINE__,'fixup 3a')
            end if
            do j=1,nc
               if (mass(j,k) > sum_mass(k)) then
                  write(*,3) 'sum_mass(k)', j, k, mass(j,k), sum_mass(k)
                  call mesa_error(__FILE__,__LINE__,'fixup 3b')
               end if
            end do
         end do

      end subroutine fix_single_point_extremes


      subroutine fix_negative_masses( &
            nz, nc, nzlo, nzhi, cell_dm, mass, ierr)

         integer, intent(in) :: nz, nc, nzlo, nzhi
         real(dp), intent(in), dimension(:) :: cell_dm
         real(dp), intent(inout), dimension(:,:) :: mass
         integer, intent(out) :: ierr

         integer :: k, j, cnt, maxcnt, k_hi, k_lo, kk, jj
         real(dp) :: dm, source_mass, frac, sum_m

         include 'formats'

         ierr = 0

         do k = nzlo,nzhi

            fix1: do j=1,nc

               if (mass(j,k) >= 0d0) cycle
               if (mass(j,k) >= -1d-13*cell_dm(k)) then
                  mass(j,k) = 0d0
                  cycle fix1
               end if

               k_hi = min(k+1,nzhi)
               k_lo = max(k-1,nzlo)
               maxcnt = 2
               do cnt = 1, maxcnt
                  sum_m = sum(mass(j,k_lo:k_hi))
                  if (sum_m >= tiny_mass) exit
                  if (cnt == maxcnt .or. mass(j,k_lo) < 0d0 .or. mass(j,k_hi) < 0d0) then
                     mass(j,k) = 0d0
                     cycle fix1
                  end if
                  k_hi = min(k_hi+1,nzhi)
                  k_lo = max(k_lo-1,nzlo)
               end do

               dm = -mass(j,k) ! dm > 0
               mass(j,k) = 0d0
               ! remove dm from neighbors
               source_mass = sum_m + dm
               frac = sum_m/source_mass
               do kk = k_lo, k_hi
                  mass(j,kk) = mass(j,kk)*frac
               end do
               if (abs(sum_m - sum(mass(j,k_lo:k_hi))) > 1d-12*sum_m) then

                  write(*,5) 'bad (sum_m - sum mass(j,:))/sum_m', j, k, k_lo, k_hi, &
                     (sum_m - sum(mass(j,k_lo:k_hi)))/sum_m
                  write(*,1) 'sum_m - sum(mass(j,k_lo:k_hi))', sum_m - sum(mass(j,k_lo:k_hi))
                  write(*,1) 'sum(mass(j,k_lo:k_hi))', sum(mass(j,k_lo:k_hi))
                  write(*,1) 'sum_m', sum_m
                  write(*,1) 'dm', dm

                  write(*,1) 'sum_m + dm', sum_m + dm
                  write(*,1) 'frac', frac

                  write(*,1) 'dm/cell_dm(k)', dm/cell_dm(k)
                  write(*,2) 'nzlo', nzlo
                  write(*,2) 'k_lo', k_lo
                  write(*,2) 'k', k
                  write(*,2) 'k_hi', k_hi
                  do jj=k_lo,k_hi
                     write(*,3) 'mass', j, jj, mass(j,jj), mass(j,jj)/cell_dm(jj)
                  end do
                  call mesa_error(__FILE__,__LINE__,'fixup')
               end if

            end do fix1

         end do

         do k=nzlo,nzhi
            do j=1,nc
               if (mass(j,k) < 0d0) then
                  write(*,3) 'mass(j,k)', j, k, mass(j,k)
                  call mesa_error(__FILE__,__LINE__,'fix_negative_masses')
               end if
            end do
         end do

      end subroutine fix_negative_masses


      subroutine fix_species_conservation( &
            nz, nc, nzlo, nzhi, mtotal, xtotal_init, X_total_atol, X_total_rtol, &
            mass, sum_mass, bad_Xsum, bad_j, ierr)

         integer, intent(in) :: nz, nc, nzlo, nzhi
         real(dp), intent(in) :: &
            mtotal, xtotal_init(:), X_total_atol, X_total_rtol
         real(dp), intent(inout) :: mass(:,:), sum_mass(:)
         real(dp), intent(out) :: bad_Xsum
         integer, intent(out) :: bad_j, ierr

         integer :: k, j
         real(dp) :: xtotal_new, frac, err

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         bad_j = 0
         bad_Xsum = 0d0

         do j=1,nc
            if (xtotal_init(j) < 1d-50) cycle
            xtotal_new = sum(mass(j,nzlo:nzhi))/mtotal
            frac = xtotal_new/xtotal_init(j)
            err = abs(xtotal_new - xtotal_init(j)) / (X_total_atol + &
               X_total_rtol*max(xtotal_new, xtotal_init(j)))
            if (err > 1d0) then
               if (dbg) write(*,2) 'fixup err', j, err
               if (dbg) write(*,2) 'xtotal_new', j, xtotal_new
               if (dbg) write(*,2) 'xtotal_init(j)', j, xtotal_init(j)
               if (dbg) write(*,2) 'X_total_atol', j, X_total_atol
               if (dbg) write(*,2) 'X_total_rtol', j, X_total_rtol
               if (dbg) write(*,2) 'frac', j, frac
               if (dbg) call mesa_error(__FILE__,__LINE__,'fixup 1')
               bad_j = j
               bad_Xsum = err
               ierr = -1
               return
            end if
            do k=nzlo,nzhi
               mass(j,k) = mass(j,k)/frac
            end do
         end do

         do k=nzlo,nzhi
            sum_mass(k) = sum(mass(1:nc,k))
            if (sum_mass(k) < 0d0) then
               write(*,2) 'sum_mass(k)', k, sum_mass(k)
               call mesa_error(__FILE__,__LINE__,'fixup 2')
            end if
         end do

      end subroutine fix_species_conservation


      subroutine smooth_x( &
            nz, nc, nzlo, nzhi, nsmooth_x_in_fixup, &
            max_lnT_for_smooth, lnT, sum_mass, mass, ierr)

         integer, intent(in) :: nz, nc, nzlo, nzhi, nsmooth_x_in_fixup
         real(dp), intent(in) :: max_lnT_for_smooth, lnT(:)
         real(dp), intent(inout) :: sum_mass(:), mass(:,:)
         integer, intent(out) :: ierr

         integer :: rep, i, k, j
         real(dp) :: sum_m1, sum_00, sum_p1, m_m1, m_00, m_p1, &
            xavg, dm, dm_m1, dm_p1, max_lnT

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         max_lnT = max_lnT_for_smooth

         do rep = 1, nsmooth_x_in_fixup

            do i = 0, 2

               do k = nzlo+1+i, nzhi-1, 3

                  if (lnT(k) > max_lnT) exit

                  sum_m1 = sum_mass(k-1)
                  sum_00 = sum_mass(k)
                  sum_p1 = sum_mass(k+1)

                  do j=1,nc

                     m_m1 = mass(j,k-1)
                     m_00 = mass(j,k)
                     m_p1 = mass(j,k+1)
                     if (m_m1 + m_p1 < tiny_mass) cycle
                     xavg = (m_m1 + m_00 + m_p1)/(sum_m1 + sum_00 + sum_p1)
                     if (abs(xavg - 1d0) < 1d-10) cycle

                     ! find dm s.t. xavg = (m_00 + dm)/(sum_00 + dm)
                     dm = (sum_00*xavg - m_00)/(1d0 - xavg)
                     if (dm >= 0d0) then
                        dm = min(dm, 0.999d0*(m_m1 + m_p1))
                     else ! dm < 0
                        dm = -min(-dm, 0.999d0*m_00)
                     end if

                     mass(j,k) = m_00 + dm
                     if (mass(j,k) < 0d0) then
                        write(*,3) 'mass00', j, k, mass(j,k), m_00, dm
                        call mesa_error(__FILE__,__LINE__,'fixup 5')
                     end if
                     sum_00 = sum(mass(1:nc,k)) ! sum_00 + dm
                     if (sum_00 < 0d0 .or. is_bad_num(sum_00)) then
                        write(*,3) 'sum_00', j, k, sum_00
                        write(*,3) 'dm', j, k, dm
                        write(*,3) 'sum_00 - dm', j, k, sum_00 - dm
                        write(*,3) 'sum_mass(k)', j, k, sum_mass(k)
                        write(*,'(A)')
                        write(*,3) 'sum_m1', j, k, sum_m1
                        write(*,3) 'sum_p1', j, k, sum_p1
                        write(*,'(A)')
                        write(*,3) 'm_m1', j, k, m_m1
                        write(*,3) 'm_00', j, k, m_00
                        write(*,3) 'm_p1', j, k, m_p1
                        write(*,'(A)')
                        write(*,3) 'xavg', j, k, xavg
                        write(*,'(A)')

                        call mesa_error(__FILE__,__LINE__,'fixup')
                     end if
                     sum_mass(k) = sum_00

                     dm_m1 = dm*m_m1/(m_m1 + m_p1)
                     mass(j,k-1) = m_m1 - dm_m1
                     if (mass(j,k-1) < 0d0) then
                        write(*,3) 'massm1', j, k, mass(j,k-1)
                        write(*,3) 'm_m1', j, k-1, m_m1
                        write(*,3) 'm_00', j, k, m_00
                        write(*,3) 'm_p1', j, k+1, m_p1
                        write(*,3) 'xavg', j, k, xavg
                        write(*,3) 'dm', j, k, dm
                        write(*,3) 'm_00 + dm', j, k, m_00 + dm
                        write(*,3) 'dm_m1', j, k, dm_m1
                        write(*,3) 'dm/(m_m1 + m_p1)', j, k, dm/(m_m1 + m_p1)
                        write(*,3) 'm_m1/(m_m1 + m_p1)', j, k, m_m1/(m_m1 + m_p1)
                        write(*,3) 'm_m1 - dm_m1', j, k, m_m1 - dm_m1
                        call mesa_error(__FILE__,__LINE__,'fixup 6')
                     end if
                     sum_m1 = sum(mass(1:nc,k-1)) ! sum_m1 - dm_m1
                     sum_mass(k-1) = sum_m1
                     if (sum_m1 < 0d0 .or. is_bad_num(sum_m1)) then
                        write(*,2) 'sum_m1', j, k, sum_m1, dm_m1
                        call mesa_error(__FILE__,__LINE__,'fixup')
                     end if

                     dm_p1 = dm*m_p1/(m_m1 + m_p1)
                     mass(j,k+1) = m_p1 - dm_p1
                     if (mass(j,k+1) < 0d0) then
                        write(*,3) 'massp1', j, k, mass(j,k+1)
                        write(*,3) 'massm1', j, k, mass(j,k-1)
                        write(*,3) 'm_m1', j, k-1, m_m1
                        write(*,3) 'm_00', j, k, m_00
                        write(*,3) 'm_p1', j, k+1, m_p1
                        write(*,3) 'xavg', j, k, xavg
                        write(*,3) 'dm', j, k, dm
                        write(*,3) 'm_00 + dm', j, k, m_00 + dm
                        write(*,3) 'dm_p1', j, k, dm_m1
                        write(*,3) 'dm/(m_m1 + m_p1)', j, k, dm/(m_m1 + m_p1)
                        write(*,3) 'm_p1/(m_m1 + m_p1)', j, k, m_p1/(m_m1 + m_p1)
                        write(*,3) 'm_p1 - dm_p1', j, k, m_p1 - dm_p1
                        call mesa_error(__FILE__,__LINE__,'fixup 7')
                     end if
                     sum_p1 = sum(mass(1:nc,k+1)) ! sum_p1 - dm_p1
                     sum_mass(k+1) = sum_p1
                     if (sum_p1 < 0d0 .or. is_bad_num(sum_p1)) then
                        write(*,2) 'sum_p1', j, k, sum_p1, dm_p1
                        call mesa_error(__FILE__,__LINE__,'fixup')
                     end if

                     if (abs(xavg - sum(mass(j,k-1:k+1))/sum(sum_mass(k-1:k+1))) > 1d-12*xavg) then
                        write(*,3) 'bad new xavg', j, k, xavg, &
                           sum(mass(j,k-1:k+1))/sum(sum_mass(k-1:k+1))
                        call mesa_error(__FILE__,__LINE__,'fixup 8')
                     end if

                  end do

               end do

            end do

            max_lnT = max_lnT - 0.1d0

         end do

         do k=nzlo,nzhi
            if (sum_mass(k) < 0d0) then
               write(*,2) 'sum_mass(k)', k, sum_mass(k)
               call mesa_error(__FILE__,__LINE__,'fixup 3')
            end if
         end do

      end subroutine smooth_x


      subroutine get1_dX_dm( &
            k, nz, nc, nzlo, nzhi, mass, sum_mass, &
            dX_dm, dbg, ierr)
         integer, intent(in) :: k, nz, nc, nzlo, nzhi
         real(dp), intent(in) :: sum_mass(:), mass(:,:)
         real(dp), intent(inout) :: dX_dm(:,:)
         logical :: dbg
         integer, intent(out) :: ierr

         real(dp) :: slope, sm1, s00, xface_00, xface_p1, &
            dm_half, dm_00, dm_p1, dm_m1, dmbar_p1, dmbar_00, &
            x00, xm1, xp1
         integer :: j
         real(dp), parameter :: tiny_slope = 1d-10

         include 'formats'

         ierr = 0

         dm_00 = sum_mass(k)
         if (dm_00 < tiny_mass) then
            dX_dm(1:nc,k) = 0d0
            return
         end if
         dm_half = 0.5d0*dm_00

         if (nzlo < k .and. k < nzhi) then

            dm_m1 = sum_mass(k-1)
            dm_p1 = sum_mass(k+1)
            if (dm_m1 < tiny_mass .or. dm_p1 < tiny_mass) then
               dX_dm(1:nc,k) = 0d0
               return
            end if

            dmbar_00 = 0.5d0*(dm_00 + dm_m1)
            dmbar_p1 = 0.5d0*(dm_00 + dm_p1)

            do j=1,nc
               xm1 = mass(j,k-1)/dm_m1
               x00 = mass(j,k)/dm_00
               xp1 = mass(j,k+1)/dm_p1
               sm1 = (xm1 - x00)/dmbar_00
               s00 = (x00 - xp1)/dmbar_p1
               slope = 0.5d0*(sm1 + s00)
               if (sm1*s00 <= 0 .or. abs(slope) < tiny_slope) then
                  dX_dm(j,k) = 0d0
               else
                  dX_dm(j,k) = slope
                  xface_00 = x00 + slope*dm_half ! value at face(k)
                  xface_p1 = x00 - slope*dm_half ! value at face(k+1)
                  if (xface_00 > 1d0 .or. xface_00 < 0d0 .or. &
                        (xm1 - xface_00)*(xface_00 - x00) < 0 .or. &
                      xface_p1 > 1d0 .or. xface_p1 < 0d0 .or. &
                        (x00 - xface_p1)*(xface_p1 - xp1) < 0) then
                     if (abs(sm1) <= abs(s00)) then
                        dX_dm(j,k) = sm1
                     else
                        dX_dm(j,k) = s00
                     end if
                  end if
               end if
            end do

         else if (k == nzlo) then

            dm_p1 = sum_mass(k+1)
            if (dm_p1 < tiny_mass) then
               dX_dm(1:nc,k) = 0d0
               return
            end if

            dmbar_p1 = 0.5d0*(dm_00 + dm_p1)

            do j=1,nc
               x00 = mass(j,k)/dm_00
               xp1 = mass(j,k+1)/dm_p1
               slope = (x00 - xp1)/dmbar_p1
               if (abs(slope) < tiny_slope) then
                  dX_dm(j,k) = 0d0
               else
                  dX_dm(j,k) = slope
                  xface_00 = x00 + slope*dm_half ! value at face(k)
                  xface_p1 = x00 - slope*dm_half ! value at face(k+1)
                  if (xface_00 > 1d0 .or. xface_00 < 0d0 .or. &
                        (x00 - xface_p1)*(xface_p1 - xp1) < 0) then
                     dX_dm(j,k) = 0d0
                  end if
               end if
            end do

         else if (k == nzhi) then

            dm_m1 = sum_mass(k-1)
            if (dm_m1 < tiny_mass) then
               dX_dm(1:nc,k) = 0d0
               return
            end if

            dmbar_00 = 0.5d0*(dm_00 + dm_m1)

            do j=1,nc
               xm1 = mass(j,k-1)/dm_m1
               x00 = mass(j,k)/dm_00
               slope = (xm1 - x00)/dmbar_00
               if (abs(slope) < tiny_slope) then
                  dX_dm(j,k) = 0d0
               else
                  dX_dm(j,k) = slope
                  xface_00 = x00 + slope*dm_half ! value at face(k)
                  xface_p1 = x00 - slope*dm_half ! value at face(k+1)
                  if (xface_p1 > 1d0 .or. xface_p1 < 0d0 .or. &
                        (xm1 - xface_00)*(xface_00 - x00) < 0) then
                     dX_dm(j,k) = 0d0
                  end if
               end if
            end do

         else

            write(*,2) 'k bad', k
            call mesa_error(__FILE__,__LINE__,'get1_dX_dm')

         end if

         ! adjust so that sum(dX_dm) = 0
         if (sum(dX_dm(1:nc,k)) > 0d0) then
            j = maxloc(dX_dm(1:nc,k),dim=1)
         else
            j = minloc(dX_dm(1:nc,k),dim=1)
         end if
         dX_dm(j,k) = 0d0 ! remove from sum
         dX_dm(j,k) = -sum(dX_dm(1:nc,k))

         ! recheck for valid values at faces
         do j=1,nc
            x00 = mass(j,k)/dm_00
            slope = dX_dm(j,k)
            xface_00 = x00 + slope*dm_half ! value at face(k)
            xface_p1 = x00 - slope*dm_half ! value at face(k+1)
            if (xface_00 > 1d0 .or. xface_00 < 0d0 .or. &
                xface_p1 > 1d0 .or. xface_p1 < 0d0) then
               if (dbg) then ! .and. abs(slope) > 1d-10) then
                  write(*,3) 'give up on dX_dm', j, k
                  write(*,1) 'slope', slope
                  write(*,1) 'dm_half', dm_half
                  write(*,1) 'xface_00', xface_00
                  write(*,1) 'xface_p1', xface_p1
                  if (k > nzlo) then
                     dm_m1 = sum_mass(k-1)
                  write(*,1) 'xm1', mass(j,k-1)/dm_m1
                  end if
                  write(*,1) 'x00', x00
                  if (k < nzhi) then
                     dm_p1 = sum_mass(k+1)
                     write(*,1) 'xp1', mass(j,k+1)/dm_p1
                  end if
                  write(*,'(A)')
               end if
               dX_dm(1:nc,k) = 0d0
               exit
            end if
         end do

      end subroutine get1_dX_dm


      subroutine redistribute_mass( &
            s, nc, nzlo, nzhi, nz, total_num_iters, sum_mass, mass, dX_dm, &
            xtotal_init, cell_dm, X_new, dbg_in, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nc, nzlo, nzhi, nz, total_num_iters
         real(dp), intent(inout), dimension(:) :: sum_mass
         real(dp), intent(in), dimension(:) :: cell_dm, xtotal_init
         real(dp), intent(in), dimension(:,:) :: mass, dX_dm
         real(dp), intent(inout) :: X_new(:,:)
         logical, intent(in) :: dbg_in
         integer, intent(out) :: ierr

         integer :: k_source, max_iters, k, i, j
         real(dp) :: source_cell_mass, remaining_source_mass, cell_dm_k, &
            remaining_needed_mass, frac, sumX, total_source, total_moved, &
            dm0, dm1, old_sum, new_sum, mtotal, err, dm, diff_dm
         logical :: okay, dbg

         include 'formats'

         ierr = 0
         total_moved = 0d0
         dbg = dbg_in

         ! redistribute mass to make sum(X_new) = 1 for all cells
         ! this is done serially from nzlo to nzhi
         k_source = nzlo
         source_cell_mass = sum_mass(k_source)
         remaining_source_mass = source_cell_mass
         max_iters = 100

         if (dbg) write(*,2) 'nzlo', nzlo
         if (dbg) write(*,2) 'nzhi', nzhi

      cell_loop: do k=nzlo,nzhi

            !dbg = total_num_iters == 26 .and. k == 535

            cell_dm_k = cell_dm(k)
            remaining_needed_mass = cell_dm_k
            X_new(1:nc,k) = 0d0
            if (dbg) write(*,*)

         fill_loop: do i=1,max_iters

               if (dbg) write(*,4) 'remaining_needed_mass/cell_dm_k', &
                  i, k, k_source, remaining_needed_mass/cell_dm_k, &
                  remaining_needed_mass, cell_dm_k, remaining_source_mass
               if (is_bad_num(remaining_needed_mass)) call mesa_error(__FILE__,__LINE__,'redistribute_mass')
               if (is_bad_num(remaining_source_mass)) then
                  write(*,4) 'remaining_source_mass', &
                     i, k, k_source, remaining_source_mass, &
                     remaining_needed_mass, cell_dm_k
                  call mesa_error(__FILE__,__LINE__,'redistribute_mass')
               end if

               if (remaining_needed_mass <= remaining_source_mass) then
                  if (dbg) write(*,1) 'remaining_needed_mass <= remaining_source_mass', &
                     remaining_needed_mass, remaining_source_mass
                  diff_dm = remaining_needed_mass
                  dm0 = remaining_source_mass - diff_dm
                  dm1 = remaining_source_mass
                  if (dm0 < 0 .or. dm1 < 0) then
                     write(*,2) 'dm0', k, dm0
                     write(*,2) 'dm1', k, dm1
                     write(*,2) 'diff_dm', k, diff_dm
                     write(*,2) 'remaining_source_mass', k, remaining_source_mass
                     write(*,2) 'remaining_needed_mass', k, remaining_needed_mass
                     call mesa_error(__FILE__,__LINE__,'redistribute_mass')
                  end if
                  total_moved = total_moved + remaining_needed_mass
                  if (dbg) then
                     write(*,2) 'cell_dm_k', k, cell_dm_k
                     write(*,2) 'remaining_needed_mass', k, remaining_needed_mass
                     write(*,2) &
                        'cell_dm_k - remaining_needed_mass', k, cell_dm_k - remaining_needed_mass
                     write(*,2) 'sum(X_new)*cell_dm_k', k, &
                        sum(X_new(1:nc,k))*cell_dm_k
                  end if
                  do j=1,nc
                     if (dbg) write(*,3) 'init X_new(j,k)', j, k, X_new(j,k)
                     dm = integrate_mass(j,dm0,dm1,diff_dm)
                     if (dm < 0d0) then
                        write(*,3) 'dm', j, k, dm
                        write(*,'(A)')
                        call mesa_error(__FILE__,__LINE__,'redistribute_mass')
                     end if
                     X_new(j,k) = X_new(j,k) + dm/cell_dm_k
                     if (dbg .and. (X_new(j,k) > 1d0 .or. X_new(j,k) < 0d0)) then

                        write(*,3) 'bad X_new(j,k)', j, k, X_new(j,k)
                        write(*,3) 'dm/cell_dm_k', j, k, dm/cell_dm_k
                        write(*,3) 'dm', j, k, dm
                        write(*,3) 'cell_dm_k', j, k, cell_dm_k
                        write(*,3) 'remaining_needed_mass', j, k, remaining_needed_mass
                        write(*,3) 'diff_dm', j, k, diff_dm
                        write(*,3) 'dm0', j, k, dm0
                        write(*,3) 'dm1', j, k, dm1
                        write(*,3) 'remaining_source_mass', j, k, remaining_source_mass
                        write(*,'(A)')
                        call mesa_error(__FILE__,__LINE__,'redistribute_mass')
                     end if
                     remaining_needed_mass = remaining_needed_mass - dm
                  end do
                  if (dbg) write(*,1) 'final remaining_needed_mass/cell_dm_k', &
                     remaining_needed_mass/cell_dm_k, remaining_needed_mass, cell_dm_k
                  remaining_source_mass = dm0
                  remaining_needed_mass = 0d0
                  exit fill_loop
               end if

               if (dbg) write(*,1) 'use all remaining source cell mass'
               ! use all remaining source cell mass
               diff_dm = remaining_source_mass
               dm0 = 0d0
               dm1 = diff_dm
               do j=1,nc
                  dm = integrate_mass(j,dm0,dm1,diff_dm)
                  X_new(j,k) = X_new(j,k) + dm/cell_dm_k
                  if (dbg .and. X_new(j,k) > 1d0 .or. X_new(j,k) < 0d0) then
                     write(*,3) 'bad X_new(j,k)', j, k, X_new(j,k)
                     write(*,1) 'dm', dm
                     write(*,1) 'dm0', dm0
                     write(*,1) 'dm1', dm1
                     write(*,1) 'remaining_source_mass', remaining_source_mass
                     write(*,1) 'remaining_needed_mass', remaining_needed_mass
                     write(*,1) 'cell_dm_k', cell_dm_k
                     call mesa_error(__FILE__,__LINE__,'redistribute_mass')
                  end if
               end do
               if (is_bad_num(remaining_needed_mass)) then
                  write(*,4) 'remaining_needed_mass', i, k, k_source, remaining_needed_mass
                  call mesa_error(__FILE__,__LINE__,'redistribute_mass')
               end if
               total_moved = total_moved + remaining_source_mass
               remaining_needed_mass = remaining_needed_mass - remaining_source_mass
               if (remaining_needed_mass > cell_dm_k) then
                  write(*,2) 'remaining_source_mass', k, remaining_source_mass
                  write(*,2) 'remaining_needed_mass', k, remaining_needed_mass
                  write(*,2) 'cell_dm_k', k, cell_dm_k
                  call mesa_error(__FILE__,__LINE__,'redistribute_mass')
               end if

               ! go to next source cell
               k_source = k_source + 1 ! okay to allow k_source > nzhi; see integrate_mass
               source_cell_mass = sum_mass(min(nzhi,k_source))
               remaining_source_mass = source_cell_mass

               if (source_cell_mass < 0d0 .or. is_bad_num(source_cell_mass)) then
                  ierr = -1
                  s% retry_message = 'bad source cell mass in element diffusion'
                  if (s% report_ierr) then
!$OMP critical (diffusion_source_cell_mass)
                     write(*,4) 'source_cell_mass', &
                        i, k, k_source, remaining_source_mass, source_cell_mass
!$OMP end critical (diffusion_source_cell_mass)
                  end if
                  return

                  write(*,4) 'source_cell_mass', &
                     i, k, k_source, remaining_source_mass, source_cell_mass
                  call mesa_error(__FILE__,__LINE__,'redistribute_mass')
               end if

            end do fill_loop
            if (dbg) write(*,1) 'finished fill_loop'

            if (remaining_needed_mass > 0d0) then
               ierr = -1
               s% retry_message = 'bad remaining mass in element diffusion'
               if (s% report_ierr) then
!$OMP critical (diffusion_dist_cell_mass)
                  write(*,1) 'redistribute_mass: remaining_needed_mass', &
                     remaining_needed_mass
!$OMP end critical (diffusion_dist_cell_mass)
               end if
               return

               write(*,5) 'remaining_needed_mass > 0d0', k, k_source, nzhi, nz, &
                  remaining_needed_mass
               write(*,1) 'source_cell_mass', source_cell_mass
               write(*,1) 'remaining_source_mass', remaining_source_mass
               call mesa_error(__FILE__,__LINE__,'redistribute_mass')
            end if

            sumX = sum(X_new(1:nc,k))
            if (abs(sumX - 1d0) > 1d-10) then
               write(*,1) 'sum(X(k)) - 1d0', sumX - 1d0
               write(*,1) 'sum mass/source_cell_mass', sum(mass(1:nc,k_source))/source_cell_mass
               write(*,1) 'sum dX_dm k_source', sum(dX_dm(1:nc,k_source))
               write(*,2) 'k', k
               write(*,2) 'k_source', k_source
               write(*,2) 'nzlo', nzlo
               write(*,2) 'nzhi', nzhi
               write(*,2) 'total_num_iters', total_num_iters
               call mesa_error(__FILE__,__LINE__,'redistribute_mass')
            end if

            do j=1,nc
               X_new(j,k) = X_new(j,k)/sumX
            end do
            sum_mass(k) = sum_mass(k)/sumX

         end do cell_loop

         if (dbg) write(*,1) 'finished cell_loop'

         total_source = sum(sum_mass(nzlo:nzhi))
         if (abs(total_moved/total_source - 1d0) > 1d-12) then
            write(*,1) 'total_moved/total_source - 1', total_moved/total_source - 1d0
            write(*,1) 'total_source', total_source
            write(*,1) 'total_moved', total_moved
            write(*,1) 'sum cell_dm', sum(cell_dm(nzlo:nzhi))
            write(*,2) 'k_source', k_source
            write(*,2) 'nzlo', nzlo
            write(*,2) 'nzhi', nzhi
            write(*,2) 'nz', nz
            write(*,'(A)')
            call mesa_error(__FILE__,__LINE__,'redistribute_mass')
         end if

         ! check cell sums
         okay = .true.
         do k=nzlo,nzhi
            new_sum = sum(X_new(1:nc,k))
            if (abs(new_sum - 1d0) > 1d-14 .or. is_bad_num(new_sum)) then
               write(*,2) 'redistribute_mass: new_sum-1', k, new_sum-1d0
               okay = .false.
            end if
         end do
         if (.not. okay) call mesa_error(__FILE__,__LINE__,'redistribute_mass')

         ! recheck conservation
         mtotal = sum(cell_dm(nzlo:nzhi))
         okay = .true.
         do j=1,nc
            if (xtotal_init(j) < 1d-20) cycle
            old_sum = mtotal*xtotal_init(j)
            new_sum = dot_product(cell_dm(nzlo:nzhi),X_new(j,nzlo:nzhi))
            err = (new_sum - old_sum)/max(old_sum,new_sum)
            if (abs(err) > 1d-12 .or. is_bad_num(err)) then
               write(*,2) 'redistribute_mass err', j, err, old_sum, new_sum
               okay = .false.
            end if
         end do


         contains


         real(dp) function integrate_mass(j,dm0,dm1,diff)
            integer, intent(in) :: j
            real(dp), intent(in) :: dm0, dm1, diff

            real(dp) :: dm, x, slope, x0, x1, half_dm, xavg
            integer :: k

            include 'formats'

            if (k_source > nzhi) then ! reuse last source cell
               k = nzhi
               slope = 0d0
            else
               k = k_source
               slope = dX_dm(j,k)
            end if
            dm = sum_mass(k)
            if (dm < tiny_mass) then
               integrate_mass = 0d0
               return
            end if
            x = mass(j,k)/dm
            half_dm = 0.5d0*dm
            x0 = x + slope*(dm0 - half_dm)
            x1 = x + slope*(dm1 - half_dm)
            if (dm0 > dm1) then
               write(*,1) 'x0', x0
               write(*,1) 'x1', x1
               write(*,1) 'dm0', dm0
               write(*,1) 'dm1', dm1
               call mesa_error(__FILE__,__LINE__,'integrate_mass')
            end if
            xavg = min(1d0, max(0d0, 0.5d0*(x0 + x1)))
            integrate_mass = xavg*diff

         end function integrate_mass


      end subroutine redistribute_mass


      subroutine get_limit_coeffs( &
            s, nz, nzlo, nzhi, &
            gamma_full_on, gamma_full_off, &
            T_full_on, T_full_off, &
            gamma, T, limit_coeffs_face, k_max)

         type (star_info), pointer :: s
         integer, intent(in) :: nz, nzlo, nzhi
         real(dp), intent(in) :: gamma_full_on, gamma_full_off, &
            T_full_on, T_full_off
         real(dp), dimension(:), intent(in) :: gamma, T
         real(dp), intent(inout) :: limit_coeffs_face(:)
         integer, intent(out) :: k_max

         real(dp) :: lim, gamma_term, T_term, gamma_max, T_min

         integer :: k

         include 'formats'

         k_max = nzhi

         do k=nzhi,nzlo+1,-1

            gamma_max = max(gamma(k),gamma(k-1))
            if (gamma_max >= gamma_full_off) then
               gamma_term = 0
            else if (gamma_max <= gamma_full_on) then
               gamma_term = 1
            else
               gamma_term = (gamma_full_off - gamma_max) / (gamma_full_off - gamma_full_on)
            end if

            T_min = min(T(k),T(k-1))
            if (T_min >= T_full_on) then
               T_term = 1
            else if (T_min <= T_full_off) then
               T_term = 0
            else
               T_term = (T_full_off - T_min) / (T_full_off - T_full_on)
            end if

            lim = gamma_term*T_term*(1d0 - s% phase(k))
            if (lim <= 0d0) then
               limit_coeffs_face(k) = 0d0
               k_max = k-1
            else if (lim >= 1d0) then
               limit_coeffs_face(k) = 1d0
            else
               limit_coeffs_face(k) = 0.5d0*(1d0 - cospi(lim))
            end if

         end do

      end subroutine get_limit_coeffs


      subroutine setup_struct_info( &
            s, nz, nzlo, nzhi, species, nc, m, X, A, tiny_X, &
            dlnP_dm_face, dlnT_dm_face, dlnRho_dm_face, cell_dm, dm_in, &
            abar, free_e, T, lnT, rho, lnd, L_face, r_face, alfa_face, &
            class, class_chem_id, calculate_ionization, nsmooth_typical_charge, &
            min_T_for_radaccel, max_T_for_radaccel, &
            min_Z_for_radaccel, max_Z_for_radaccel, &
            screening, rho_face, T_face, four_pi_r2_rho_face, &
            dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face, &
            Z, typical_charge, xm_face, &
            rad_accel_face, log10_g_rad, g_rad, &
            kmax_rad_accel, ierr)

         type (star_info), pointer :: s
         integer, intent(in) :: nz, nzlo, nzhi, species, nc, m, &
            min_Z_for_radaccel, max_Z_for_radaccel
         real(dp), intent(in) :: X(:,:), A(:)
         real(dp), intent(in) :: tiny_X, &
            min_T_for_radaccel, max_T_for_radaccel
         real(dp), dimension(:), intent(in) :: &
            dlnP_dm_face, dlnT_dm_face, dlnRho_dm_face, cell_dm, dm_in, &
            abar, free_e, T, lnT, rho, lnd, L_face, r_face, alfa_face
         integer, dimension(:), intent(in) :: class, class_chem_id
         logical, intent(in) :: calculate_ionization, screening
         integer, intent(in) :: nsmooth_typical_charge
         real(dp), dimension(:), intent(out) :: &
            rho_face, T_face, four_pi_r2_rho_face, &
            dlnP_dr_face, dlnT_dr_face, dlnRho_dr_face
         real(dp), dimension(:,:), intent(out) :: Z
         real(dp), dimension(:,:), intent(out) :: &
            typical_charge, rad_accel_face, log10_g_rad, g_rad
         real(dp), dimension(:), intent(out) :: xm_face
         integer, intent(out) :: kmax_rad_accel, ierr

         integer :: i, k, j, op_err, kmax

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         if (dbg) write(*,*) 'call do1 for each zone'

         if (calculate_ionization) then
!$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,2)
            do k=nzlo,nzhi
               call do1(k)
            end do
!$OMP END PARALLEL DO
         else
            do k=nzlo,nzhi
               call do1(k)
            end do
         end if

         if (nsmooth_typical_charge > 0) then
            do j=1,nsmooth_typical_charge
               do i=1,nc
                  typical_charge(i,nzlo) = &
                     (2*typical_charge(i,nzlo) + typical_charge(i,nzlo+1))/3
                  do k = nzlo+1, nzhi-2
                     typical_charge(i,k) = &
                        (typical_charge(i,k-1) + typical_charge(i,k) + typical_charge(i,k+1))/3
                  end do
                  typical_charge(i,nzhi-1) = &
                     (2*typical_charge(i,nzhi-1) + typical_charge(i,nzhi-2))/3
                  do k = nzhi-2, nzlo+1, -1
                     typical_charge(i,k) = &
                        (typical_charge(i,k-1) + typical_charge(i,k) + typical_charge(i,k+1))/3
                  end do
                  typical_charge(i,nzhi) = typical_charge(i,nzhi-1)
               end do
            end do
         end if

         xm_face(nzlo) = 0d0
         do k=nzlo,nzhi
            do i=1, nc
               Z(i,k) = typical_charge(i,k)
            end do
            Z(m,k) = -1
            if (k < nzhi) xm_face(k+1) = xm_face(k) + cell_dm(k)
         end do

         if (T_face(nzlo+1) > max_T_for_radaccel) then
            kmax_rad_accel = 0
            return
         end if

         kmax = nzhi
         do k=nzlo+1,nzhi
            if (T_face(k) > max_T_for_radaccel) then
               kmax = k-1
               exit
            end if
         end do
         kmax_rad_accel = kmax

         if (dbg) write(*,*) 'call calc_g_rad'
         call calc_g_rad( &
            nz, nzlo, nzhi, nc, m, kmax_rad_accel, X, A, &
            class_chem_id, s% net_iso, s% op_mono_factors, &
            L_face, rho_face, r_face, T_face, alfa_face, &
            min_T_for_radaccel, max_T_for_radaccel, &
            min_Z_for_radaccel, max_Z_for_radaccel, &
            screening, log10_g_rad, g_rad, &
            rad_accel_face, ierr)
         if (dbg) write(*,*) 'done calc_g_rad'

         return
         do k=nzlo,nzlo+9
            write(*,2) 'alfa_face(k)', k, alfa_face(k)
            write(*,2) 'rho_face(k)', k, rho_face(k)
            write(*,2) 'r_face(k)', k, r_face(k)
            write(*,2) 'four_pi_r2_rho_face(k)', k, four_pi_r2_rho_face(k)
            write(*,2) 'dlnRho_dr_face(k)', k, dlnRho_dr_face(k)
            write(*,2) 'dlnT_dr_face(k)', k, dlnT_dr_face(k)
         end do


         contains


         subroutine do1(k)
            use mod_typical_charge, only: eval_typical_charge
            integer, intent(in) :: k
            integer :: i

            if (k > nzlo) then
               T_face(k) = alfa_face(k)*T(k) + (1d0-alfa_face(k))*T(k-1)
               rho_face(k) = alfa_face(k)*rho(k) + (1d0-alfa_face(k))*rho(k-1)
               four_pi_r2_rho_face(k) = pi4*r_face(k)*r_face(k)*rho_face(k)
               dlnP_dr_face(k) = four_pi_r2_rho_face(k)*dlnP_dm_face(k)
               dlnT_dr_face(k) = four_pi_r2_rho_face(k)*dlnT_dm_face(k)
               dlnRho_dr_face(k) = four_pi_r2_rho_face(k)*dlnRho_dm_face(k)
            end if

            if (calculate_ionization) then
               do i=1, nc
                  typical_charge(i,k) = eval_typical_charge( &
                        class_chem_id(i), abar(k), free_e(k), &
                        T(k), lnT(k)/ln10, rho(k), lnd(k)/ln10)
               end do
            end if

            do i=1,nc
               rad_accel_face(i,k) = 0
               g_rad(i,k) = 0
            end do
            rad_accel_face(m,k) = 0

         end subroutine do1

      end subroutine setup_struct_info


      subroutine calc_g_rad( &
            nz, nzlo, nzhi, nc, m, kmax_rad_accel, X, A, &
            class_chem_id, net_iso, op_mono_factors, &
            L_face, rho_face, r_face, T_face, alfa_face, &
            min_T_for_radaccel, max_T_for_radaccel, &
            min_Z_for_radaccel, max_Z_for_radaccel, &
            screening, log10_g_rad, g_rad, &
            rad_accel_face, ierr)

         use kap_lib, only: get_op_mono_params
         use utils_lib, only: utils_OMP_GET_MAX_THREADS, utils_OMP_GET_THREAD_NUM

         integer, intent(in) :: &
            nz, nzlo, nzhi, nc, m, kmax_rad_accel, class_chem_id(:), net_iso(:), &
            min_Z_for_radaccel, max_Z_for_radaccel
         real(dp), intent(in) :: &
            min_T_for_radaccel, max_T_for_radaccel
         real(dp), dimension(:), intent(in) :: &
            A, L_face, rho_face, r_face, T_face, alfa_face, op_mono_factors
         real(dp), dimension(:,:), intent(in) :: X
         logical, intent(in) :: screening
         real(dp), dimension(:,:), intent(out) :: &
            log10_g_rad, g_rad, rad_accel_face
         integer, intent(out) :: ierr

         integer :: iZ(nc), iZ_rad(nc), n, i, j, k, kk, kmax, op_err, sz, offset, &
            nptot, ipe, nrad, thread_num
         real(dp) :: alfa, beta, X_face(nc)
         real, pointer, dimension(:) :: umesh1, semesh1, ff1, ta1, rs1
         real, pointer :: &
            umesh(:), semesh(:), ff(:,:,:,:), ta(:,:,:,:), rs(:,:,:)

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         kmax = kmax_rad_accel

         g_rad(1:nc,kmax+1:nzhi) = 0
         log10_g_rad(1:nc,kmax+1:nzhi) = 1

         g_rad(1:nc,nzlo) = 0
         log10_g_rad(1:nc,nzlo) = 1

         iZ(1:nc) = chem_isos% Z(class_chem_id(1:nc))

         kk = 0
         do i = 1, nc
            if (iZ(i) >= min_Z_for_radaccel .and. &
                iZ(i) <= max_Z_for_radaccel) then
               kk = kk+1
               iZ_rad(kk) = iZ(i)
            end if
         end do

         if (dbg) write(*,*) 'call get_op_mono_params'
         call get_op_mono_params(nptot, ipe, nrad)
         if (dbg) write(*,*) 'done get_op_mono_params'
         n = utils_OMP_GET_MAX_THREADS()
         if (dbg) write(*,2) 'nptot', nptot
         if (dbg) write(*,2) 'ipe', ipe
         if (dbg) write(*,2) 'nrad', nrad
         if (dbg) write(*,2) 'kmax', kmax
         if (dbg) write(*,2) 'nzlo', nzlo
         if (dbg) write(*,2) 'n', n
         if (dbg) write(*,2) 'nptot*ipe*4*4*n', nptot*ipe*4*4*n
         if (nptot*ipe*4*4*n < 0) call mesa_error(__FILE__,__LINE__,'integer overflow for array size in calc_g_rad')
         allocate( &
            umesh1(nptot*n), semesh1(nptot*n), ff1(nptot*ipe*4*4*n), &
            ta1(nptot*nrad*4*4*n), rs1(nptot*4*4*n), &
            stat=ierr)
         if (ierr /= 0) return

!$OMP PARALLEL DO PRIVATE(i,k,thread_num,op_err,j,alfa,beta,X_face,umesh,ff,ta,rs,sz,offset) SCHEDULE(dynamic,2)
         do k=nzlo+1, kmax
            if (T_face(k) < min_T_for_radaccel) then
               do j = 1, nc
                  rad_accel_face(j,k) = 0d0
               end do
               rad_accel_face(m,k) = 0d0
               cycle
            end if
            i = k - nzlo
            alfa = alfa_face(k)
            beta = 1d0 - alfa
            do j = 1, nc
               X_face(j) = alfa*X(j,k) + beta*X(j,k-1)
            end do
            op_err = 0
            thread_num = utils_OMP_GET_THREAD_NUM()
            sz = nptot; offset = thread_num*sz
            umesh(1:nptot) => umesh1(offset+1:offset+sz)
            semesh(1:nptot) => semesh1(offset+1:offset+sz)
            sz = nptot*ipe*4*4; offset = thread_num*sz
            ff(1:nptot,1:ipe,1:4,1:4) => ff1(offset+1:offset+sz)
            sz = nptot*nrad*4*4; offset = thread_num*sz
            ta(1:nptot,1:nrad,1:4,1:4) => ta1(offset+1:offset+sz)
            sz = nptot*4*4; offset = thread_num*sz
            rs(1:nptot,1:4,1:4) => rs1(offset+1:offset+sz)
            sz = nptot*nrad*4*4; offset = thread_num*sz

            if (dbg) write(*,2) 'call set1_g_rad', k
            if (dbg) stop
            call set1_g_rad( &
               k, nc, iZ, kk, iZ_rad, T_face(k), rho_face(k), &
               L_face(k), r_face(k), A, X_face, screening, &
               min_Z_for_radaccel, max_Z_for_radaccel, &
               class_chem_id, net_iso, op_mono_factors, &
               umesh, semesh, ff, ta, rs, &
               log10_g_rad(:,k), g_rad(:,k), &
               op_err)
            if (op_err /= 0) ierr = op_err
            do j = 1, nc
               rad_accel_face(j,k) = g_rad(j,k)
            end do
            rad_accel_face(m,k) = 0d0
         end do
!$OMP END PARALLEL DO

         deallocate(umesh1, semesh1, ff1, ta1, rs1)

         k = nzlo
         do j = 1, m
            rad_accel_face(j,k) = rad_accel_face(j,k+1) ! for plotting
         end do

      end subroutine calc_g_rad


      subroutine set1_g_rad( &
         k, nc, iZ, kk, iZ_rad, T, rho, L, r, A, X, screening, &
         min_Z_for_radaccel, max_Z_for_radaccel, &
         class_chem_id, net_iso, op_mono_factors, &
         umesh, semesh, ff, ta, rs, &
         log10_g_rad, g_rad, ierr)

         use kap_lib, only: op_mono_get_radacc

         integer, intent(in) :: k, nc, kk, class_chem_id(:), net_iso(:), &
            min_Z_for_radaccel, max_Z_for_radaccel
         integer, dimension(:), intent(in) :: iZ(:) ! (nc)
         integer, dimension(:), intent(in) :: iZ_rad(:) ! (kk)
         real(dp), intent(in) :: T, rho, L, r
         real(dp), dimension(:), intent(in) :: A, X, op_mono_factors
         logical, intent(in) :: screening
         real, pointer :: &
            umesh(:), semesh(:), ff(:,:,:,:), ta(:,:,:,:), rs(:,:,:)
         real(dp), dimension(:), intent(out) :: &
            log10_g_rad, g_rad
         integer, intent(out) :: ierr

         integer :: i, j, ii
         real(dp) :: tot, logT, logRho, flux, logKappa
         real(dp), dimension(nc) :: &
            fa, fac, lgrad

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         tot = 0
         fa = 0
         do i=1,nc
            fa(i) = X(i)/A(i)
            tot = tot + fa(i)
            j = net_iso(class_chem_id(i))
            if (j /= 0) then
               fac(i) = op_mono_factors(j)
            else
               write(*,*) 'bad class_chem_id? in set1_g_rad'
               call mesa_error(__FILE__,__LINE__,'set1_g_rad')
               fac(i) = 1
            end if
         end do
         do i=1,nc
            fa(i) = fa(i)/tot ! number fractions
         end do

         logT = log10(T)
         logRho = log10(rho)
         flux = L/(pi4*r*r)

         if (dbg) write(*,*) 'call op_mono_get_radacc'
         call op_mono_get_radacc( &
            ! input
            kk, iZ_rad, nc, iZ, fa, fac, &
            flux, logT, logRho, screening, &
            ! output
            logKappa, &
            lgrad, &
            ! work arrays
            umesh, semesh, ff, ta, rs, &
            ierr)
         if (dbg) write(*,*) 'done op_mono_get_radacc'

         do i=1,nc
            g_rad(i) = 0d0
            log10_g_rad(i) = 1d0
         end do

         if (ierr == 101) then ! logRho out of range
            ierr = 0
            write(*,2) 'op_mono_get_radacc bad logT logRho', k, logT, logRho
            return
         end if

         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'set1_g_rad')  !return

         do ii = 1, kk
            do i = 1, nc
               if (iZ(i) == iZ_rad(ii)) then
                  log10_g_rad(i) = lgrad(ii)
                  g_rad(i) = exp10(lgrad(ii))
                  exit
               end if
            end do
         end do

      end subroutine set1_g_rad


      subroutine update_rad_accel_face( &
            nzlo, nzhi, nc, m, A, X_init, X, &
            log10_g_rad, g_rad, &
            rad_accel_face, kmax_rad_accel)
         integer, intent(in) :: nzlo, nzhi, nc, m, kmax_rad_accel
         real(dp), intent(in) :: A(:) ! (nc)
         real(dp), dimension(:,:), intent(in) :: &
            X_init, X, log10_g_rad, g_rad
         real(dp), dimension(:,:), intent(inout) :: rad_accel_face

         integer :: j, k, i

         include 'formats'

         do k=nzlo+1,kmax_rad_accel           
            do i=1,nc
               rad_accel_face(i,k) = pow(10d0, log10_g_rad(i,k))
            end do
            rad_accel_face(m,k) = 0d0
         end do

         k = nzlo
         do j = 1, m
            rad_accel_face(j,k) = rad_accel_face(j,k+1) ! for plotting
         end do

      end subroutine update_rad_accel_face


      subroutine set_new_xa( &
            nz, nzlo, nzhi, species, nc, m, class, X_init, X, cell_dm, xa)
         integer, intent(in) :: nz, nzlo, nzhi, species, nc, m, class(:)
         real(dp), intent(in) :: X_init(:,:) ! (nc,nz)
         real(dp), intent(in) :: X(:,:) ! (m,nz)
         real(dp), intent(in) :: cell_dm(:) ! (nz)
         real(dp), intent(inout) :: xa(:,:) ! (species,nz)
         integer :: j, k, i
         real(dp) :: tmp
         include 'formats'
         do k=nzlo,nzhi
            do j=1,species
               i = class(j)
               if (X_init(i,k) <= 0 .or. is_bad_num(X_init(i,k))) then
                  write(*,3) 'X_init(i,k)', i, k, X_init(i,k)
                  call mesa_error(__FILE__,__LINE__,'set_new_xa')
               end if
               xa(j,k) = xa(j,k)*X(i,k)/X_init(i,k)
            end do
            tmp = sum(xa(1:species,k))
            if (tmp <= 0 .or. is_bad_num(tmp)) then
               write(*,2) 'tmp', k, tmp
               call mesa_error(__FILE__,__LINE__,'set_new_xa')
            end if
            do j=1,species
               xa(j,k) = xa(j,k)/tmp
            end do
         end do
      end subroutine set_new_xa


      end module diffusion_procs

