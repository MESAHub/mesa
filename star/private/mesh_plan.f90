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


      module mesh_plan

      use const_def
      use num_lib
      use utils_lib
      use star_private_def

      implicit none

      private
      public :: do_mesh_plan

      logical, parameter :: plan_dbg = .false.
      integer, parameter :: kdbg = -1


      contains


      subroutine do_mesh_plan( &
            s, nz_old, max_allowed_nz, okay_to_merge, D_mix, &
            max_k_old_for_split_in, min_k_old_for_split_in, xq_old, dq_old, &
            min_dq_in, max_dq, min_dq_for_split, mesh_max_allowed_ratio, &
            do_not_split, num_gvals, gval_names, &
            gval_is_xa_function, gval_is_logT_function, gvals, &
            delta_gval_max, max_center_cell_dq, max_surface_cell_dq, &
            max_num_subcells, max_num_merge_cells, &
            nz_new, xq_new, dq_new, which_gval, comes_from, ierr)
         ! return keep_going, or terminate
         use mesh_functions, only: max_allowed_gvals
         ! inputs
         type (star_info), pointer :: s
         integer, intent(in) :: nz_old, max_allowed_nz, max_num_subcells, &
            max_num_merge_cells, max_k_old_for_split_in, min_k_old_for_split_in
         logical, intent(in) :: okay_to_merge
         real(dp), pointer :: D_mix(:) ! (nz_old)
         real(dp), pointer :: xq_old(:) ! (nz_old)
         real(dp), pointer :: dq_old(:) ! (nz_old)
         real(dp), intent(in) :: min_dq_in, max_dq, min_dq_for_split, mesh_max_allowed_ratio
         logical, pointer :: do_not_split(:)
         integer, intent(in) :: num_gvals
         character (len=32) :: gval_names(max_allowed_gvals)
         logical, dimension(max_allowed_gvals) :: gval_is_xa_function, gval_is_logT_function
         real(dp), pointer :: gvals(:,:) ! (nz_old, num_gvals)
         real(dp), pointer :: delta_gval_max(:) ! (nz_old)
         real(dp), intent(in) :: max_center_cell_dq, max_surface_cell_dq
         ! outputs
         integer, intent(out) :: nz_new
         real(dp), pointer :: xq_new(:), dq_new(:) ! (nz_new)
            ! must be allocated on entry; suggested size >= nz_old.
            ! reallocated as necessary if need to enlarge.
            ! size on return is >= nz_new.
         integer, pointer :: which_gval(:) ! (nz_new)  for debugging.
            ! which_gval(k) = gval number that set the size for new cell k.
            ! size may have been reduced below gradient setting by other restrictions.
         integer, pointer :: comes_from(:) ! (nz_new)
            ! xq_old(comes_from(k)+1) > xq_new(k) >= xq_old(comes_from(k)), if comes_from(k) < nz_old.
         integer, intent(out) :: ierr

         integer :: j, k, k_old, k_new, nz, new_capacity, iounit, species, &
            max_num_merge_surface_cells, max_k_old_for_split, min_k_old_for_split
         real(dp) :: D_mix_cutoff, next_xq, next_dq, max_dq_cntr, &
            dq_sum, tmp, min_dq, min_dq_for_xa, min_dq_for_logT

         logical, parameter :: write_plan_debug = .false.

         include 'formats'

         ierr = 0
         
         min_dq = min_dq_in
         
         if (max_k_old_for_split_in < 0) then
            max_k_old_for_split = nz_old + max_k_old_for_split_in
         else
            max_k_old_for_split = max_k_old_for_split_in
         end if
         
         if (min_k_old_for_split_in < 0) then
            min_k_old_for_split = nz_old + min_k_old_for_split_in
         else
            min_k_old_for_split = min_k_old_for_split_in
         end if

         max_num_merge_surface_cells = max_num_merge_cells ! for now

         if (max_dq < min_dq) then
            write(*,1) 'ERROR in controls: max_dq < min_dq', max_dq, min_dq
            ierr = -1
            return
         end if

         if (max_center_cell_dq < min_dq) then
            write(*,1) 'ERROR in controls: max_center_cell_dq < min_dq', max_center_cell_dq, min_dq
            ierr = -1
            return
         end if

         if (max_surface_cell_dq < min_dq) then
            write(*,1) 'ERROR in controls: max_surface_cell_dq < min_dq', max_surface_cell_dq, min_dq
            ierr = -1
            return
         end if

         do k = 2, nz_old
            if (xq_old(k) <= xq_old(k-1)) then
               write(*,3) 'bad xq_old', k, nz_old, xq_old(k), xq_old(k-1)
               ierr = -1
               return
            end if
         end do

         if (s% set_min_D_mix) then
            D_mix_cutoff = s% min_D_mix
         else
            D_mix_cutoff = 0
         end if

         if (write_plan_debug) call open_debug_file

         if (plan_dbg) then
            write(*,*) 'num_gvals', num_gvals
            do j=1,num_gvals
               write(*,*) j, trim(gval_names(j))
            end do
            write(*,'(A)')
         end if

         nz = nz_old
         species = s% species

         new_capacity = min(size(xq_new,dim=1), size(dq_new,dim=1), &
                  size(which_gval,dim=1), size(comes_from,dim=1))
         max_dq_cntr = max(1d-20, max_center_cell_dq, 0.75d0*dq_old(nz_old))

         comes_from(1:nz) = 0
         comes_from(1) = 1

         call pick_new_points(s, ierr)
         if (ierr /= 0) return
         
         if (min_k_old_for_split <= 1) then
            do while (dq_new(1) > max(max_surface_cell_dq,2*min_dq,min_dq_for_split))
               call split1(1, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed trying to split surface cell'
                  stop
                  return
               end if
            end do
         end if

         call smooth_new_points(ierr) ! split as necessary
         if (ierr /= 0) return

         if (write_plan_debug) then
            close(iounit)
         end if

         if (ierr /= 0) return


         contains


         subroutine check_before_smooth(ierr)
            integer, intent(out) :: ierr
            integer :: k, k_old
            real(dp) :: xq0, xqold_start, xqold_end

            include 'formats'
            ierr = 0
            if (comes_from(1) /= 1 .or. xq_old(1) /= xq_new(1)) then
               write(*,*) 'mesh plan check_before_smooth: bad value from k=1'
               ierr = -1
               return
            end if

            do k = 2, nz_new
               xq0 = xq_new(k)
               k_old = comes_from(k)
               xqold_start = xq_old(k_old)
               if (k_old == nz_old) then
                  xqold_end = 1d0
               else
                  xqold_end = xq_old(k_old+1)
               end if
               if (k_old > comes_from(k-1)) then
                  if (xq0 /= xqold_start) then
                     write(*,'(A)')
                     write(*,2) 'nz_new', nz_new
                     write(*,2) 'k', k
                     write(*,2) 'comes_from(k-1)', comes_from(k-1)
                     write(*,2) 'k_old', k_old
                     write(*,2) 'comes_from(k)', comes_from(k)
                     write(*,2) 'nz_old', nz_old
                     write(*,'(A)')
                     write(*,2) 'xq0', k, xq0
                     write(*,2) 'xqold_start', k_old, xqold_start
                     write(*,2) 'xqold_end', k_old, xqold_end
                     write(*,2) 'xqold_end - xqold_start', k_old, xqold_end - xqold_start
                     write(*,2) 'xq0 - xqold_start', k_old, xq0 - xqold_start
                     write(*,'(A)')
                     write(*,2) 'xq_old(comes_from(k))', k, xq_old(comes_from(k))
                     write(*,2) 'xq_old(comes_from(k-1))', k, xq_old(comes_from(k-1))
                     write(*,'(A)')
                     write(*,*) 'k_old > comes_from(k-1)', k_old > comes_from(k-1)
                     write(*,*) 'xq0 /= xqold_start', xq0 /= xqold_start
                     write(*,*) 'mesh plan check_before_smooth'
                     ierr = -1
                     return
                  end if
               else if (k_old == comes_from(k-1)) then
                  if (.not. (xq0 > xqold_start .and. xq0 < xqold_end)) then
                     write(*,'(A)')
                     write(*,*) '(.not. (xq0 > xqold_start .and. xq0 < xqold_end))', k_old
                     write(*,2) 'xq_new(k)', k, xq_new(k)
                     write(*,2) 'xq_old(k_old)', k, xq_old(k_old)
                     write(*,2) 'xq_old(k_old+1)', k, xq_old(k_old+1)
                     write(*,2) 'xqold_end', k, xqold_end
                     write(*,'(A)')
                     write(*,*) 'k_old == comes_from(k-1)', k_old == comes_from(k-1)
                     write(*,*) '.not. (xq0 > xqold_start .and. xq0 < xqold_end)', &
                        .not. (xq0 > xqold_start .and. xq0 < xqold_end)
                     write(*,*) 'mesh plan check_before_smooth'
                     ierr = -1
                     return
                  end if
               else
                  write(*,*) 'comes_from(k) > comes_from(k-1)', k, comes_from(k), comes_from(k-1)
                  write(*,*) 'mesh plan check_before_smooth'
                  ierr = -1
                  return
               end if

               cycle  ! use the following for debugging
               if (dq_new(k-1) < 1d-6*dq_new(k) .or. dq_new(k) < 1d-6*dq_new(k-1)) then
                  write(*,3) 'bad dq_new ratio', k, nz_new, dq_new(k)/dq_new(k-1), dq_new(k), dq_new(k-1)
                  write(*,*) 'check_before_smooth'
                  ierr = -1
                  return
               end if

            end do

         end subroutine check_before_smooth


         subroutine test_new(ierr)
            integer, intent(out) :: ierr
            integer :: k, k_old
            include 'formats'
            ierr = 0
            do k = 1, nz_new
               if (dq_new(k) <= 0) then
                  write(*,3) 'bad dq_new'
                  write(*,2) 'dq_new', k, dq_new(k)
                  write(*,2) 'dq_new', k-1, dq_new(k-1)
                  write(*,2) 'xq_new', k, xq_new(k)
                  write(*,2) 'xq_new', k-1, xq_new(k-1)
                  write(*,3) 'comes_from', k, comes_from(k)
                  write(*,3) 'comes_from', k-1, comes_from(k-1)
                  write(*,2) 'nz_new', nz_new
                  write(*,2) 'nz_old', nz_old
                  write(*,2) 'dq_old(comes_from(k))', comes_from(k), dq_old(comes_from(k))
                  write(*,2) 'dq_old(comes_from(k)-1)', comes_from(k)-1, dq_old(comes_from(k)-1)
                  write(*,2) 'dq_old(comes_from(k)-2)', comes_from(k)-2, dq_old(comes_from(k)-2)
                  write(*,2) 'xq_old(comes_from(k))', comes_from(k), xq_old(comes_from(k))
                  write(*,2) 'xq_old(comes_from(k)-1)', comes_from(k)-1, xq_old(comes_from(k)-1)
                  write(*,2) 'xq_old(comes_from(k)-2)', comes_from(k)-2, xq_old(comes_from(k)-2)
                  write(*,*) 'test_new: mesh plan'
                  ierr = -1
                  return
               end if
               k_old = comes_from(k)
               if (k_old < 1 .or. k_old > nz_old) then
                  write(*,*) 'bad value for comes_from', k, k_old
                  write(*,*) 'test_new: mesh plan'
                  ierr = -1
                  return
               end if
               if (xq_new(k) < xq_old(k_old)) then
                  write(*,3) '(xq_new(k) < xq_old(k_old))', k, k_old, &
                     xq_new(k) - xq_old(k_old), xq_new(k), xq_old(k_old)
                  write(*,*) 'test_new: mesh plan'
                  ierr = -1
                  return
               end if
               if (k_old < nz_old) then
                  if (xq_new(k) > xq_old(k_old+1)) then
                     write(*,3) '(xq_new(k) > xq_old(k_old+1))', k, k_old+1, &
                        xq_new(k) - xq_old(k_old+1), xq_new(k), xq_old(k_old+1)
                     write(*,*) 'test_new: mesh plan'
                     ierr = -1
                     return
                  end if
               end if
            end do
         end subroutine test_new


         subroutine split1(k, ierr)
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            integer :: kk, k_old, k_old_last, split_at_k_old
            real(dp) :: xq_mid, xq_end
            logical :: from_merger, dbg

            include 'formats'
            ierr = 0

            k_old = comes_from(k)
            k_old_last = -1
            xq_end = -1
            from_merger = (xq_new(k) == xq_old(k_old) .and. dq_new(k) > dq_old(k_old) + min_dq*1d-3)

            dbg = (k_old == -1)

            if (dbg) then
               write(*,2) 'start split1 k', k
               write(*,2) 'nz_old', nz_old
               write(*,2) 'nz_new', nz_new
               do kk=1400,nz_new
                  if (dq_new(kk) < 1d-12) then
                     write(*,2) 'dq_new(kk)', kk, dq_new(kk)
                     call mesa_error(__FILE__,__LINE__,'debug: split1')
                  end if
               end do
               do kk = 2, nz_new
                  if (dq_new(kk-1) < 1d-6*dq_new(kk) .or. dq_new(kk) < 1d-6*dq_new(kk-1)) then
                     write(*,3) 'bad dq_new ratio', kk, nz_new, dq_new(kk), dq_new(kk-1)
                     call mesa_error(__FILE__,__LINE__,'debug: split1')
                  end if
               end do
            end if

            if (from_merger) then ! find range of old cells that were merged to form k
               if (k == nz_new) then
                  xq_end = 1
                  k_old_last = nz_old
               else
                  xq_end = xq_new(k+1)
                  k_old_last = 0
                  do kk = k_old+1, nz_old ! find last old cell included in k
                     if (xq_old(kk) == xq_end) then
                        k_old_last = kk-1; exit
                     end if
                     if (xq_old(kk) > xq_end) then
                        write(*,*) 'oops'
                        write(*,2) 'xq_old(kk)', kk, xq_old(kk)
                        write(*,2) 'xq_old(kk-1)', kk-1, xq_old(kk-1)
                        write(*,2) 'xq_end', k, xq_end
                        write(*,*) 'split1'
                        ierr = -1
                        return
                     end if
                  end do
                  if (k_old_last == k_old) then
                     from_merger = .false.
                  else if (k_old_last < k_old) then
                     write(*,*) 'confusion in split1 for k_old_last'
                     write(*,2) 'k_old', k_old
                     write(*,2) 'k_old_last', k_old_last
                     write(*,2) 'nz_old', nz_old
                     write(*,2) 'k', k
                     write(*,2) 'nz_new', nz_new
                     write(*,'(A)')
                     write(*,2) 'dq_new(k)', k, dq_new(k)
                     write(*,2) 'dq_old(k_old)', k_old, dq_old(k_old)
                     write(*,2) 'dq_new(k)-dq_old(k_old)', k_old, dq_new(k)-dq_old(k_old)
                     write(*,'(A)')
                     write(*,2) 'xq_new(k)', k, xq_new(k)
                     write(*,2) 'xq_new(k)+dq_new(k)', k, xq_new(k)+dq_new(k)
                     write(*,2) 'xq_end', k, xq_end
                     write(*,*) 'split1'
                     ierr = -1
                     return
                  end if
               end if
            end if

            if (nz_new == new_capacity) then ! increase allocated size
               new_capacity = (new_capacity*5)/4 + 10
               call realloc(s, nz_new, new_capacity, xq_new, dq_new, which_gval, comes_from, ierr)
               if (ierr /= 0) return
            end if
            nz_new = nz_new + 1
            if (nz_new > max_allowed_nz) then
               write(*,*) 'tried to increase number of mesh points beyond max allowed nz', max_allowed_nz
               ierr = -1
               return
            end if
            do kk = nz_new, k+1, -1
               xq_new(kk) = xq_new(kk-1)
               dq_new(kk) = dq_new(kk-1)
               which_gval(kk) = which_gval(kk-1)
               comes_from(kk) = comes_from(kk-1)
            end do

            if (from_merger) then ! split by breaking up the merger
               xq_mid = xq_new(k) + dq_new(k)*0.5d0
               split_at_k_old = k_old_last
               do kk = k_old+1, k_old_last ! check interior boundaries for closest to xq_mid
                  if (xq_old(kk) >= xq_mid) then
                     if (xq_old(kk) - xq_mid < xq_mid - xq_old(kk-1)) then ! xq_mid closer to xq_old(kk)
                        split_at_k_old = kk
                     else ! xq_mid closer to xq_old(kk-1)
                        split_at_k_old = kk-1
                     end if
                     exit
                  end if
               end do
               comes_from(k+1) = split_at_k_old
               xq_new(k+1) = xq_old(split_at_k_old)
               dq_new(k) = sum(dq_old(k_old:split_at_k_old-1))
               dq_new(k+1) = sum(dq_old(split_at_k_old:k_old_last))
            else
               dq_new(k:k+1) = dq_new(k)*0.5d0
               xq_new(k+1) = xq_new(k) + dq_new(k)
               ! fix up comes_from(k+1)
               comes_from(k+1) = nz_old ! just in case
               do k_old = comes_from(k), nz_old-1
                  if (xq_new(k+1) <= xq_old(k_old+1)) then
                     comes_from(k+1) = k_old
                     exit
                  end if
               end do
            end if

            if (dbg) then
               write(*,2) 'split1 dq_new k', k, dq_new(k)
               write(*,2) 'split1 dq_new k+1', k+1, dq_new(k+1)
               write(*,2) 'comes_from', k_old
               write(*,2) 'nz_old', nz_old
               write(*,2) 'nz_new', nz_new
               do kk=1400,nz_new
                  if (dq_new(kk) < 1d-12) then
                     write(*,2) 'dq_new(kk)', kk, dq_new(kk)
                     call mesa_error(__FILE__,__LINE__,'debug: split1')
                  end if
               end do
            end if

         end subroutine split1


         logical function okay_to_split1(k_old, dq_new, remaining_dq_old)
            integer, intent(in) :: k_old
            real(dp), intent(in) :: dq_new, remaining_dq_old
            real(dp) :: dlnR_old, dr_old, min_dr, rR, rL, dlnR_new
            logical :: dbg

            include 'formats'

            dbg = .false.
            okay_to_split1 = .false.
            if (do_not_split(k_old)) return

            if (max(dq_new,remaining_dq_old) < max(min_dq_for_split,2*min_dq)) return

            if (dbg) then
               write(*,2) 'remaining_dq_old', k_old, remaining_dq_old
               write(*,1) 'dq_new', dq_new
               write(*,1) 'min_dq_for_split', min_dq_for_split
               write(*,1) 'min_dq', min_dq
               write(*,'(A)')
            end if

            if (0d0 < remaining_dq_old .and. dq_new > 0.99d0*remaining_dq_old) then
               if (dbg) then
                  write(*,2) 'dq_old', k_old, remaining_dq_old
                  write(*,1) 'dq_new', dq_new
                  write(*,1) 'dq_new/remaining_dq_old', dq_new/remaining_dq_old
                  write(*,'(A)')
               end if
               return
            end if

            if (k_old < nz_old) then

               rR = s% r(k_old)
               rL = s% r(k_old+1)
               dr_old = rR - rL
               min_dr = s% csound(k_old)*s% mesh_min_dr_div_cs
               if (dr_old*dq_new/dq_old(k_old) < 2*min_dr) then
                  return ! sound crossing time would be too small
               end if

               min_dr = s% mesh_min_dr_div_dRstar*(s% r(1) - s% R_center)
               if (dr_old*dq_new/dq_old(k_old) < 2*min_dr) then
                  return ! new dr would be too small
               end if

               if (s% mesh_min_dlnR > 0d0) then
                  dlnR_new = log((rL + dr_old*dq_new/dq_old(k_old))/rL)
                  if (dlnR_new < 2*s% mesh_min_dlnR) then
                     ! the factor of 2 is a safety margin.
                     return
                  end if
               end if

            end if
            okay_to_split1 = .true.

         end function okay_to_split1


         subroutine smooth_new_points(ierr)
            integer, intent(out) :: ierr

            logical :: dbg, done
            integer :: k, k_old
            real(dp) :: alfa

            include 'formats'

            ierr = 0
            dbg = .false.
            alfa = mesh_max_allowed_ratio

            do ! repeat until nothing left to do
               if (min_k_old_for_split <= 1 .and. &
                     dq_new(1) > alfa*dq_new(2) .and. &
                     okay_to_split1(1,dq_new(1),0d0)) then
                  call split1(1,ierr)
                  if (ierr /= 0) return
                  cycle
               end if
               done = .true.
               k = 2
               do ! check for cell that is too large
                  if (k == nz_new) exit
                  k_old = comes_from(k)
                  if (k_old <= max_k_old_for_split .and. &
                      k_old >= min_k_old_for_split .and. &
                      okay_to_split1(k_old,dq_new(k),0d0) .and. &
                        (dq_new(k) > max_dq .or. &
                           dq_new(k) > alfa*dq_new(k+1) .or. &
                           dq_new(k) > alfa*dq_new(k-1))) then
                     call split1(k,ierr)
                     if (ierr /= 0) return
                     done = .false.
                     ! don't increment k; want to recheck the cell in case need to split again
                  else
                     k = k + 1
                  end if
               end do

               if (dbg) then
                  call test_new(ierr)
                  if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'debug: mesh_plan, smooth_new_points')
               end if

               if (done) exit

               ! now go in opposite direction
               done = .true.
               k = nz_new-1
               do
                  if (k == 1) exit
                  if (okay_to_split1(comes_from(k),dq_new(k),0d0) .and. &
                        (dq_new(k) > alfa*dq_new(k+1) .or. dq_new(k) > alfa*dq_new(k-1))) then
                     call split1(k,ierr)
                     if (ierr /= 0) return
                     done = .false.
                     k = k + 1 ! recheck the same cell
                  else
                     k = k - 1
                  end if
               end do

               if (dbg) then
                  call test_new(ierr)
                  if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'debug: mesh_plan, smooth_new_points')
               end if

               if (done) exit

            end do

         end subroutine smooth_new_points


         subroutine pick_new_points(s, ierr)
            type (star_info), pointer :: s
            integer, intent(out) :: ierr

            logical :: dbg, force_merge_with_one_more
            real(dp) :: dqsum, prev_dq, dq_limit, maxval_delta_xa, next_dq_max, beta_limit, &
               remaining_dq_old, min_dr
            integer :: kk, k_old_init, k_old_next, k_old_next_max, j00, jm1, i, max_merge

            include 'formats'

            beta_limit = 0.1d0

            ierr = 0
            k_old = 1
            k_new = 1
            xq_new(1) = 0
            min_dr = s% mesh_min_dr_div_dRstar*(s% r(1) - s% R_center)

            do ! pick next point location

               dbg = plan_dbg .or. (k_new == kdbg) !.or. (s% mesh_call_number == 2005)

               ! when reach this point,
               ! have set xq_new(k) for k = 1 to k_new, and dq_new(k) for k = 1 to k_new-1.
               ! and have finished using old points from 1 to k_old
               ! i.e., xq_old(k_old+1) > xq_new(k_new) >= xq_old(k_old)

               k_old_init = k_old

               if (k_new == 1) then
                  next_dq_max = max_surface_cell_dq
               else
                  next_dq_max = dq_new(k_new-1)*mesh_max_allowed_ratio
               end if

               ! make initial choice based on gradients. may reduce later.
               if (dbg) then
                  write(*,'(A)')
                  write(*,3) 'call pick_next_dq', k_old, k_new, next_dq_max
               end if
               
               if (s% gradr(k_old) > s% grada(k_old) .and. &
                     s% min_dq_for_xa_convective > 0d0) then
                  min_dq_for_xa = s% min_dq_for_xa_convective
               else
                  min_dq_for_xa = s% min_dq_for_xa
               end if
               
               min_dq_for_logT = s% min_dq_for_logT

               next_dq = pick_next_dq(s, &
                  dbg, next_dq_max, k_old, k_new, nz_old, num_gvals, &
                  xq_new, dq_new, xq_old, dq_old, min_dq, min_dq_for_split, &
                  min_dq_for_xa, min_dq_for_logT, &
                  max_surface_cell_dq, max_dq_cntr, max_num_subcells, &
                  gval_is_xa_function, gval_is_logT_function, gvals, &
                  delta_gval_max, gval_names, which_gval, ierr)
               if (ierr /= 0) return

               if (dbg) &
                  write(*,3) 'pick_next_dq next_dq prev_dq', k_old, k_new, next_dq, dq_old(k_old)

               if (k_new == 1) then
                  max_merge = max_num_merge_surface_cells
               else
                  max_merge = max_num_merge_cells
               end if

               if (k_old < nz_old) then
                  remaining_dq_old = xq_old(k_old+1)-xq_new(k_new)
               else
                  remaining_dq_old = 1d0-xq_new(k_new)
               end if

               if (next_dq < remaining_dq_old) then
                  if (.not. okay_to_split1(k_old, next_dq, remaining_dq_old)) then
                     next_dq = remaining_dq_old ! dq_old(k_old)
                  end if
               end if

               if (next_dq > max_dq) then
                  if (xq_new(k_new) == xq_old(k_old)) then
                     do i = nz_old, k_old+1, -1
                        if (xq_old(i) - xq_new(k_new) <= max_dq) then
                           next_dq = xq_old(i) - xq_new(k_new)
                           exit
                        end if
                     end do
                  end if
               end if

               if (next_dq > max_dq) then
                  next_dq = xq_old(k_old) + dq_old(k_old) - xq_new(k_new)
               end if

               if (k_new == 1 .and. next_dq > max_surface_cell_dq) then
                  if (k_old_init == -1 .and. .true.) write(*,*) 'next_dq > max_surface_cell_dq'
                  next_dq = max_surface_cell_dq
               end if

               next_xq = xq_new(k_new) + next_dq
               if (next_xq > 1 - min_dq) then
                  next_xq = (1 + xq_new(k_new))/2
                  if (k_old < nz_old) then ! make sure don't split current k_old for this case
                     if (xq_old(k_old+1) > next_xq) next_xq = xq_old(k_old+1)
                  end if
                  next_dq = next_xq - xq_new(k_new)
               end if

               if (k_old < nz_old) then

                  if (xq_new(k_new) == xq_old(k_old) .and. &
                        next_dq > dq_old(k_old) - min_dq/2) then

                     if (.not. okay_to_merge) then
                     
                        k_old_next = k_old + 1
                        
                     else if (k_old < min_k_old_for_split .or. &
                              k_old > max_k_old_for_split) then
                     
                        k_old_next = k_old + 1
                              
                     else ! consider doing merge

                        if (next_dq > 1.5d0*dq_old(k_old)) then
                           next_dq = 0.9d0*next_dq ! to avoid split-merge flip-flops
                           next_xq = xq_new(k_new) + next_dq
                        end if

                        k_old_next_max = min(nz_old, k_old + max_merge)
                        k_old_next = k_old_next_max ! will cut this back as necessary
                        do kk=k_old+1,k_old_next_max
                           maxval_delta_xa = maxval(abs(s% xa(:,kk)-s% xa(:,kk-1)))
                           j00 = maxloc(s% xa(:,kk),dim=1)
                           jm1 = maxloc(s% xa(:,kk-1),dim=1)
                           if (maxval_delta_xa > s% max_delta_x_for_merge .or. &
                               j00 /= jm1 .or. is_convective_boundary(kk) .or. &
                               is_crystal_boundary(kk)) then
                              ! don't merge across convective or crystal boundary
                              k_old_next = kk-1
                              exit
                           else if (next_xq <= xq_old(kk) + min_dq/2) then
                              k_old_next = max(k_old+1,kk-1)
                              exit
                           end if
                        end do

                     end if

                     k_old_next = max(k_old_next, k_old+1)
                     next_xq = xq_old(k_old_next)
                     next_dq = next_xq - xq_new(k_new)

                  else if (next_xq >= xq_old(k_old+1) - min_dq/2) then
                     ! this is final subcell of a split, so adjust to finish the parent cell

                     k_old_next = k_old+1
                     next_xq = xq_old(k_old_next)
                     next_dq = next_xq - xq_new(k_new)

                  else ! non-final subcell of split

                     k_old_next = k_old

                  end if

                  if (next_xq == xq_old(k_old_next)) then
                     ! finishing 1 or more old cells and not already at max merge.
                     ! consider forcing a merge with the next cell to make this one larger.
                     force_merge_with_one_more = .false.

                     if (dq_old(k_old) < min_dq) then
                        force_merge_with_one_more = .true.
                     else if (s% merge_if_dlnR_too_small) then
                        if (xq_new(k_new) <= xq_old(k_old) .and. &
                            s% lnR(k_old) - s% lnR(k_old_next) < s% mesh_min_dlnR) then
                           force_merge_with_one_more = .true.
                        else if (k_old_next == nz_old .and. s% R_center > 0) then
                           force_merge_with_one_more = s% lnR(k_old_next) - log(s% R_center) < s% mesh_min_dlnR
                        end if
                     end if

                     if ((.not. force_merge_with_one_more) .and. s% merge_if_dr_div_cs_too_small) then
                        if (xq_new(k_new) <= xq_old(k_old) .and. &
                            s% r(k_old) - s% r(k_old_next) < s% mesh_min_dr_div_cs*s% csound(k_old)) then
                           force_merge_with_one_more = .true.
                        else if (k_old_next == nz_old) then
                           force_merge_with_one_more = (s% r(k_old_next) - s% R_center) < &
                                       s% mesh_min_dr_div_cs*s% csound(k_old_next)
                           if (force_merge_with_one_more .and. dbg) &
                              write(*,3) 'do merge for k_old_next == nz_old', k_old, k_old_next, &
                                 s% r(k_old) - s% R_center, &
                                 s% mesh_min_dr_div_cs*s% csound(k_old_next), &
                                 s% mesh_min_dr_div_cs, s% csound(k_old_next)
                        end if
                     end if

                     if ((.not. force_merge_with_one_more) .and. &
                           s% merge_if_dr_div_dRstar_too_small) then
                        if (xq_new(k_new) <= xq_old(k_old) .and. &
                            s% r(k_old) - s% r(k_old_next) < min_dr) then
                           force_merge_with_one_more = .true.
                        else if (k_old_next == nz_old) then
                           force_merge_with_one_more = (s% r(k_old_next) - s% R_center) < min_dr
                           if (force_merge_with_one_more .and. dbg) &
                              write(*,3) 'do merge for k_old_next == nz_old', k_old, k_old_next, &
                                 s% r(k_old) - s% R_center, min_dr
                        end if
                     end if

                     if (force_merge_with_one_more) then
                        k_old_next = k_old_next + 1
                        if (k_old_next < nz_old) then
                           next_xq = xq_old(k_old_next)
                        else
                           next_xq = 1d0
                        end if
                        next_dq = next_xq - xq_new(k_new)
                        if (dbg) then
                           write(*,3) 'force merge', k_old, k_new, next_xq, next_dq
                           write(*,'(A)')
                        end if
                     end if

                  end if

                  k_old = k_old_next

               end if

               comes_from(k_new) = k_old_init

               ! check if we're done
               if (1 - xq_new(k_new) < max_dq_cntr .or. 1 - next_xq < min_dq) then
                  dq_new(k_new) = 1 - xq_new(k_new)
                  exit
               end if

               dq_new(k_new) = next_dq

               if (k_new == new_capacity) then ! increase allocated size
                  new_capacity = (new_capacity*5)/4 + 10
                  call realloc(s, k_new, new_capacity, xq_new, dq_new, which_gval, comes_from, ierr)
                  if (ierr /= 0) return
               end if

               if (next_xq < xq_new(k_new)) then
                  write(*,*) 'nz_old', nz_old
                  write(*,*) 'k_new', k_new
                  write(*,1) 'next_xq', next_xq
                  write(*,1) 'xq_new(k_new)', xq_new(k_new)
                  write(*,*) 'pick_new_points: next_xq < xq_new(k_new)'
                  ierr = -1
                  return
               end if

               dq_sum = sum(dq_new(1:k_new))

               k_new = k_new + 1
               if (k_new > max_allowed_nz) then
                  write(*,*) 'tried to increase number of mesh points beyond max allowed nz', max_allowed_nz
                  ierr = -1
                  return
               end if

               xq_new(k_new) = next_xq
               if (abs(xq_new(k_new) - dq_sum) > 1d-6) then
                  write(*,'(A)')
                  write(*,*) 'k_new', k_new
                  write(*,1) 'xq_new(k_new) - dq_sum', xq_new(k_new) - dq_sum
                  write(*,1) 'xq_new(k_new)', xq_new(k_new)
                  write(*,1) 'dq_sum', dq_sum
                  write(*,*) 'pick_new_points: abs(xq_new(k_new) - dq_sum) > 1d-6'
                  ierr = -1
                  return
               end if
               !write(*,2) 'xq_new(k_new) - dq_sum', k_new, xq_new(k_new) - dq_sum

               ! increment k_old if necessary
               do while (k_old < nz_old)
                  if (xq_old(k_old+1) > next_xq) exit
                  k_old = k_old + 1
               end do

            end do

            nz_new = k_new

            if (plan_dbg) write(*,2) 'after pick_new_points: nz_new', nz_new

         end subroutine pick_new_points


         logical function is_convective_boundary(kk)
            integer, intent(in) :: kk
            is_convective_boundary = .false.
            if (kk == nz) return
            is_convective_boundary = &
               (s% mixing_type(kk) == convective_mixing .and. &
                s% mixing_type(kk+1) /= convective_mixing) .or. &
               (s% mixing_type(kk+1) == convective_mixing .and. &
                s% mixing_type(kk) /= convective_mixing)
         end function is_convective_boundary

         logical function is_crystal_boundary(kk)
            integer, intent(in) :: kk
            ! replace this with mixing type designation when implemented for crystal
            if(s% m(kk) <= s% crystal_core_boundary_mass .and. &
                 s% m(kk-1) >= s% crystal_core_boundary_mass) then
               is_crystal_boundary = .true.
            else
               is_crystal_boundary = .false.
            end if
         end function is_crystal_boundary

         subroutine open_debug_file
            include 'formats'
            open(newunit=iounit, file=trim('plan_debug.data'), action='write', iostat=ierr)
            if (ierr /= 0) then
               write(*, *) 'open plan_debug.data failed'
               call mesa_error(__FILE__,__LINE__,'debug do_mesh_plan')
            end if
            write(*,*) 'write plan_debug.data'
         end subroutine open_debug_file


      end subroutine do_mesh_plan


      subroutine realloc(s, old_size, new_capacity, xq_new, dq_new, which_gval, comes_from, ierr)
         use alloc
         type (star_info), pointer :: s
         integer, intent(in) :: old_size, new_capacity
         real(dp), pointer :: xq_new(:), dq_new(:)
         integer, pointer :: which_gval(:), comes_from(:)
         integer, intent(out) :: ierr
         integer, parameter :: extra = 100
         call realloc_integer_work_array(s, which_gval, old_size, new_capacity, extra, ierr)
         if (ierr /= 0) return
         call realloc_integer_work_array(s, comes_from, old_size, new_capacity, extra, ierr)
         if (ierr /= 0) return
         call realloc_work_array(s, .false., xq_new, old_size, new_capacity, extra, 'mesh_plan', ierr)
         if (ierr /= 0) return
         call realloc_work_array(s, .false., dq_new, old_size, new_capacity, extra, 'mesh_plan', ierr)
         if (ierr /= 0) return
      end subroutine realloc


      real(dp) function pick_next_dq(s, &
            dbg, next_dq_max, k_old, k_new, nz_old, num_gvals, xq_new, dq_new, &
            xq_old, dq_old, min_dq, min_dq_for_split, min_dq_for_xa, min_dq_for_logT, &
            max_surface_cell_dq, max_dq_cntr, max_num_subcells, &
            gval_is_xa_function, gval_is_logT_function, gvals, delta_gval_max, &
            gval_names, which_gval, ierr)
         use mesh_functions, only: max_allowed_gvals
         type (star_info), pointer :: s
         logical, intent(in) :: dbg
         integer, intent(in) :: k_old, k_new, nz_old, num_gvals, max_num_subcells
         real(dp), pointer :: xq_new(:), dq_new(:) ! (nz)
         real(dp), pointer :: xq_old(:) ! (nz_old)
         real(dp), pointer :: dq_old(:) ! (nz_old)
         logical, dimension(max_allowed_gvals) :: gval_is_xa_function, gval_is_logT_function
         real(dp), pointer :: gvals(:,:) ! (nz_old, num_gvals)
         real(dp), intent(in) :: &
            next_dq_max, min_dq, min_dq_for_split, &
            min_dq_for_xa, min_dq_for_logT, max_surface_cell_dq, max_dq_cntr
         real(dp), pointer :: delta_gval_max(:) ! (nz_old, num_gvals)
         character (len=32) :: gval_names(:) ! (num_gvals)  for debugging.
         integer, pointer :: which_gval(:) ! (nz_new)  for debugging.
         integer, intent(out) :: ierr

         real(dp) :: nxt_dqs(num_gvals), default
         integer :: i, j, jmin, op_err
         logical :: pkdbg

         include 'formats'

         pkdbg = dbg !.or. k_old == 4000
         ierr = 0

         if (pkdbg) write(*,*)
         if (pkdbg) write(*,2) 'dq_old(k_old)', k_old, dq_old(k_old)

         if (k_new == 1) then
            pick_next_dq = min(1d0, sqrt(min_dq))
         else if (k_new <= 20) then
            pick_next_dq = min(1-xq_new(k_new), 10*dq_new(k_new-1))
         else
            pick_next_dq = 1-xq_new(k_new)
         end if

         nxt_dqs(:) = pick_next_dq
         if (k_old == nz_old) then
            which_gval(k_new) = 0
            pick_next_dq = max(min_dq, (1-xq_new(k_new)) - max_dq_cntr)
            do i=1,10
               if (1-xq_new(k_new) <= max_dq_cntr*i) then
                  pick_next_dq = (1-xq_new(k_new))/i
                  exit
               end if
            end do
            return
         end if
         
         default = pick_next_dq ! default size. can be reduced according to gradients of gvals
         do j=1,num_gvals
            nxt_dqs(j) = pick1_dq(s, &
               j, next_dq_max, default, .false., k_old, k_new, nz_old, &
               xq_new, xq_old, dq_old, min_dq, min_dq_for_xa, min_dq_for_logT, max_num_subcells, &
               gval_is_xa_function(j), gval_is_logT_function(j), &
               gvals, delta_gval_max, gval_names, op_err)
            if (op_err /= 0) ierr = op_err
         end do

         jmin = minloc(nxt_dqs(:),dim=1)
         if (pkdbg) write(*,3) 'jmin, k_new, init pick_next_dq', jmin, k_new, pick_next_dq
         which_gval(k_new) = jmin
         pick_next_dq = max(min_dq, min(pick_next_dq, nxt_dqs(jmin)))

      end function pick_next_dq


      real(dp) function pick1_dq(s, &
            j, next_dq_max, default, dbg, k_old, k_new, nz_old, &
            xq_new, xq_old, dq_old, min_dq, min_dq_for_xa, min_dq_for_logT, max_num_subcells, &
            is_xa_function, is_logT_function, gvals, delta_gval_max, gval_names, ierr)
         use num_lib, only: binary_search
         type (star_info), pointer :: s
         integer, intent(in) :: j
         real(dp), intent(in) :: next_dq_max, default, min_dq, min_dq_for_xa, min_dq_for_logT
         logical, intent(in) :: dbg
         integer, intent(in) :: k_old, k_new, nz_old, max_num_subcells
         real(dp), pointer :: xq_new(:) ! (nz)
         real(dp), pointer :: xq_old(:) ! (nz_old)
         real(dp), pointer :: dq_old(:) ! (nz_old)
         logical :: is_xa_function, is_logT_function
         real(dp), pointer :: gvals(:,:) ! (nz_old, num_gvals)
         real(dp), pointer :: delta_gval_max(:) ! (nz_old, num_gvals)
         character (len=32) :: gval_names(:) ! (num_gvals)  for debugging.
         integer, intent(out) :: ierr

         real(dp) :: gnew, dgnew, gmax, gmin, xq, dq_next, dq_sum, sz, dval
         integer :: k
         logical :: dbg1

         include 'formats'

         ierr = 0
         dbg1 = dbg !.or. (k_old == -1)

         pick1_dq = default
         dq_next = -1

         if (dbg1) write(*,*)
         if (dbg1) write(*,*) 'find max allowed dq for gvals(:,j)', j
         ! find max allowed dq for gvals(:,j)

         ! linear interpolate to estimate gvals(:,j) at q_new(k_new)
         dval = gvals(k_old,j) - gvals(k_old+1,j)

         gnew = gvals(k_old+1,j) + dval*(xq_old(k_old+1) - xq_new(k_new))/dq_old(k_old)
         if (dbg1) write(*,*) trim(gval_names(j))

         dgnew = min(delta_gval_max(k_old), delta_gval_max(k_old+1))
         gmax = gnew + dgnew
         gmin = gnew - dgnew

         dq_sum = 0
         do k=k_old+1, nz_old

            if (dbg1) write(*,2) 'gvals(k,j)', k, gvals(k,j), dq_sum, dq_old(k-1)

            if (gvals(k,j) <= gmax .and. gvals(k,j) >= gmin) then
               if (xq_old(k-1) >= xq_new(k_new)) then
                  dq_sum = dq_sum + dq_old(k-1)
               else
                  if (dq_sum /= 0) then
                     write(*,*) '(dq_sum /= 0)'
                     write(*,*) 'pick1_dq'
                     ierr = -1
                     return
                  end if
                  dq_sum = xq_old(k) - xq_new(k_new)
               end if
               if (is_bad(dq_sum)) then
                  write(*,2) 'dq_sum', k, dq_sum
                  write(*,*) 'pick1_dq'
                  ierr = -1
                  if (s% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'mesh plan')
                  return
               end if
               if (dq_sum >= next_dq_max) exit
               if (dq_sum >= default) then
                  dq_sum = default; exit
               end if
               if (k < nz_old) cycle
               ! pick location inside center zone
               if (dbg1) write(*,*) 'pick location inside center zone'
               if (gvals(k-1,j) < gvals(k,j)) then ! see where reach gmax in center cell
                  dq_next = find0(0d0, gvals(k-1,j)-gmax, dq_old(k-1), gvals(k,j)-gmax)
                  if (dbg1) write(*,1) 'gvals(k-1,j) < gvals(k,j)', dq_next
               else if (gvals(k-1,j) > gvals(k,j)) then ! see where reach gmin
                  dq_next = find0(0d0, gvals(k-1,j)-gmin, dq_old(k-1), gvals(k,j)-gmin)
                  if (dbg1) then
                     write(*,1) 'gvals(k-1,j)-gmin', gvals(k-1,j)-gmin
                     write(*,1) 'gvals(k,j)-gmin', gvals(k,j)-gmin
                     write(*,1) 'dq_old(k-1)', dq_old(k-1)
                     write(*,1) 'gvals(k-1,j) > gvals(k,j)', dq_next
                     call mesa_error(__FILE__,__LINE__,'debug pick1_dq')
                  end if
               else ! we're done -- don't need another point for this gval
                  dq_sum = default; exit ! just return the default
               end if
               if (dq_next > 1 - (xq_new(k_new) + min_dq)) then
                  dq_sum = default
               else
                  dq_sum = dq_sum + dq_next
               end if
               exit
            end if

            if (gvals(k,j) > gmax) then ! estimate where = gmax
               dq_next = find0(0d0, gvals(k-1,j)-gmax, dq_old(k-1), gvals(k,j)-gmax)
            else if (gvals(k,j) < gmin) then ! estimate where = gmin
               dq_next = find0(0d0, gvals(k-1,j)-gmin, dq_old(k-1), gvals(k,j)-gmin)
            end if
            if (is_bad(dq_next)) then
               write(*,2) 'dq_next', k, dq_next
               write(*,*) 'gvals(k,j) > gmax', gvals(k,j) > gmax
               write(*,*) 'gvals(k,j) < gmin', gvals(k,j) < gmin
               write(*,3) 'gvals(k-1,j)', k-1, j, gvals(k-1,j)
               write(*,2) 'gmax', k, gmax
               write(*,3) 'gvals(k,j)', k, j, gvals(k,j)
               write(*,2) 'gmin', k, gmin
               write(*,2) 'dq_old(k-1)', k, dq_old(k-1)
               write(*,*) 'pick1_dq'
               ierr = -1
               if (s% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'mesh plan')
               return
            end if
            if (xq_old(k-1) >= xq_new(k_new)) then
               dq_sum = dq_sum + dq_next
            else
               if (dq_sum /= 0) then
                  write(*,*) '(dq_sum /= 0)'
                  write(*,*) 'pick1_dq'
                  ierr = -1
                  return
               end if
               dq_sum = dq_next - (xq_new(k_new) - xq_old(k-1))
            end if
            exit

         end do

         if (dbg1) then
            write(*,1) 'after loop: dq_sum', dq_sum, min_dq, default
         end if

         dq_sum = max(min_dq, dq_sum)
         xq = xq_new(k_new) + dq_sum
         if (dbg1) write(*,1) 'before round_off_xq', xq
         xq = round_off_xq(k_old, nz_old, max_num_subcells, xq, xq_old, dq_old, sz, ierr)
         if (ierr /= 0) return
         if (dbg1) then
            write(*,1) 'after round_off_xq: xq', xq
            write(*,1) 'xq_new(k_new)', xq_new(k_new)
            write(*,1) 'xq - xq_new(k_new)', xq - xq_new(k_new)
            write(*,1) 'sz', sz
            write(*,1) 'min_dq', min_dq
            write(*,1) 'default', default
         end if
         pick1_dq = max(min_dq, sz, xq - xq_new(k_new))
         
         if (is_xa_function .and. pick1_dq < min_dq_for_xa) &
            pick1_dq = min_dq_for_xa
         
         if (is_logT_function .and. pick1_dq < min_dq_for_logT) &
            pick1_dq = min_dq_for_logT

         if (dbg1) then
            write(*,2) 'dq_sum', k_new, dq_sum
            write(*,2) 'xq', k_new, xq
            write(*,2) 'pick1_dq', k_new, pick1_dq
            write(*,2) 'log pick1_dq', k_new, log10(pick1_dq)
         end if

      end function pick1_dq


      real(dp) function round_off_xq(k_old, nz_old, n, xq, xq_old, dq_old, sz, ierr)
         ! adjust to match one of the candidate subcell locations
         ! this prevents generating too many candidate new points
         integer, intent(in) :: k_old, nz_old, n ! n is number of subcells
         real(dp), intent(in) :: xq, xq_old(:), dq_old(:)
         real(dp), intent(out) :: sz ! subcell size at new location
         integer, intent(out) :: ierr

         real(dp) :: dq, tmp
         integer :: i, k, j, knxt
         include 'formats'
         ierr = 0
         knxt = 0
         round_off_xq = -1
         if (xq >= xq_old(nz_old)) then
            knxt = nz_old+1
         else
            j = -1
            do k=k_old,1,-1
               if (xq_old(k) <= xq) then
                  j = k+1; exit
               end if
            end do
            if (j <= 1) then
               write(*,*) 'logic error in mesh plan: round_off_xq'
               ierr = -1
               return
            end if
            do k=j,nz_old
               if (xq_old(k) > xq) then
                  knxt = k; exit
               end if
            end do
         end if
         ! xq is in old cell knxt-1
         ! move location to next subcell boundary
         dq = xq - xq_old(knxt-1)
         sz = dq_old(knxt-1)/dble(n) ! size of subcells
         tmp = dq/sz
         if(tmp>huge(n)) tmp=huge(n)
         i = max(1,floor(tmp))
         if (knxt > nz_old) i = min(i,n/2) ! limit extrapolation into center
         dq = i*sz
         round_off_xq = xq_old(knxt-1) + dq
         if (dq <= 0) then
            write(*,2) 'dq', knxt-1, dq, xq, xq_old(knxt-1), dq_old(knxt-1)
            write(*,3) 'i n sz', i, n, sz
            write(*,*) 'round_off_xq'
            ierr = -1
            return
         end if
      end function round_off_xq


      end module mesh_plan
