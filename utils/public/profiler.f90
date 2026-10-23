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

! Region based wall clock profiler. Regions are opened and closed with
! prof_push/prof_pop and form a call tree: a region is identified by its
! name together with its parent.
!
! The profiler is always compiled in but only records anything when the
! environment variable MESA_PROFILE is set to a value other than '', '0',
! 'no' or 'false'. Otherwise prof_push, prof_pop and prof_report return
! immediately.
!
! Regions are meant to be placed outside of OpenMP parallel regions.
! Calls made from inside a parallel region are ignored, so a routine
! that is called both serially and from parallel loops can be
! instrumented safely; only its serial calls are timed.
!
! TODO(astrobarker): r vs s labeling

module profiler

   use const_def, only: dp, i8
   use utils_lib, only: mesa_error

   implicit none

   private
   public :: prof_init, prof_enabled, prof_push, prof_pop, prof_report

   ! max_regions: max number of distinct (name, parent) regions
   ! max_depth: max nesting depth of open regions
   integer, parameter :: max_regions = 512, max_depth = 64, name_len = 64

   ! per region data, indexed by region id
   character(len=name_len) :: rname(max_regions)
   integer :: rparent(max_regions) = 0  ! 0 for top level regions
   integer(i8) :: rcalls(max_regions) = 0
   integer(i8) :: rincl(max_regions) = 0  ! inclusive clock ticks
   integer(i8) :: rexcl(max_regions) = 0  ! exclusive clock ticks (excluding children)
   integer :: nr = 0

   ! stack of open regions
   integer :: sid(max_depth)  ! region id
   integer(i8) :: sstart(max_depth)  ! clock at push
   integer(i8) :: schild(max_depth)  ! ticks spent in children so far
   integer :: depth = 0

   integer(i8) :: t_first = 0
   logical :: started = .false.

   ! set from MESA_PROFILE by prof_init, on first use
   logical :: initialized = .false., enabled = .false.

contains

   ! read MESA_PROFILE; called automatically on first use, but drivers may
   ! call it earlier
   subroutine prof_init()
      character(len=32) :: val
      integer :: length, status, i

      if (initialized) return
      initialized = .true.

      call get_environment_variable('MESA_PROFILE', val, length, status)
      ! status -1 means truncated, which is still a non-empty value
      if (status > 0 .or. length == 0) return

      ! lower case for comparison
      do i = 1, len_trim(val)
         if (val(i:i) >= 'A' .and. val(i:i) <= 'Z') &
            val(i:i) = achar(iachar(val(i:i)) + 32)
      end do
      enabled = .not. (val == '0' .or. val == 'no' .or. val == 'false')
   end subroutine prof_init

   logical function prof_enabled()
      if (.not. initialized) call prof_init()
      prof_enabled = enabled
   end function prof_enabled

   logical function in_parallel()
!$    use omp_lib, only: omp_in_parallel
      in_parallel = .false.
!$    in_parallel = omp_in_parallel()
   end function in_parallel

   subroutine prof_push(name)
      character(len=*), intent(in) :: name
      integer :: parent, id, i

      if (in_parallel()) return
      if (.not. prof_enabled()) return

      parent = 0
      if (depth > 0) parent = sid(depth)

      ! find the region for (name, parent), adding it if new
      id = 0
      do i = 1, nr
         if (rparent(i) == parent .and. rname(i) == name) then
            id = i
            exit
         end if
      end do

      if (id == 0) then
         if (nr == max_regions) &
            call mesa_error(__FILE__, __LINE__, 'profiler: too many regions (raise max_regions)')
         nr = nr + 1
         id = nr
         rname(id) = name
         rparent(id) = parent
      end if

      if (depth == max_depth) &
         call mesa_error(__FILE__, __LINE__, 'profiler: nesting too deep (missing prof_pop? or raise max_depth)')

      depth = depth + 1
      rcalls(id) = rcalls(id) + 1
      sid(depth) = id
      schild(depth) = 0

      ! read the clock last so the bookkeeping above is not timed
      call system_clock(sstart(depth))
      if (.not. started) then
         t_first = sstart(depth)
         started = .true.
      end if
   end subroutine prof_push

   subroutine prof_pop()
      integer(i8) :: t, elapsed
      integer :: id

      if (in_parallel()) return
      if (.not. prof_enabled()) return

      if (depth == 0) &
         call mesa_error(__FILE__, __LINE__, 'profiler: prof_pop with no open region')

      ! read the clock first so the bookkeeping below is not timed
      call system_clock(t)

      id = sid(depth)
      elapsed = t - sstart(depth)

      ! inclusive time is the full elapsed time; exclusive time removes the
      ! time spent in child regions, which is then credited to the parent
      rincl(id) = rincl(id) + elapsed
      rexcl(id) = rexcl(id) + elapsed - schild(depth)
      depth = depth - 1
      if (depth > 0) schild(depth) = schild(depth) + elapsed
   end subroutine prof_pop

   subroutine prof_report()
      integer(i8) :: rate, t
      real(dp) :: total
      integer :: i

      if (.not. prof_enabled()) return
      if (.not. started) then
         write(*,'(a)') 'Profiler: no regions recorded'
         return
      end if
      call system_clock(t, rate)
      total = max(real(t - t_first, dp)/real(rate, dp), tiny(1d0))

      ! a non-empty stack means a push without a matching pop; the time of
      ! regions that are still open has not been recorded
      if (depth > 0) then
         write(*,'(a,i0,a)') 'WARNING: profiler stack not empty, ', depth, ' region(s) still open:'
         do i = 1, depth
            write(*,'(a,a)') repeat('  ', i), trim(rname(sid(i)))
         end do
      end if

      write(*,'(a)') 'time  %total  %remainder  calls  name'
      write(*,'(a)') '==================='
      call print_tree(0, 0, total, real(rate, dp))
      write(*,'(a)') '==================='
   end subroutine prof_report

   ! %total is inclusive time relative to the time since the first push.
   ! %remainder is exclusive/inclusive time, i.e. the fraction of a region's
   ! time not covered by its instrumented children.
   recursive subroutine print_tree(parent, level, total, rate)
      integer, intent(in) :: parent, level
      real(dp), intent(in) :: total, rate
      integer :: kids(max_regions), nk, i, j, k, tmp, nchild
      character(len=6) :: pct, rem
      real(dp) :: ti

      ! collect children of parent, sorted by inclusive time (descending)
      nk = 0
      do i = 1, nr
         if (rparent(i) /= parent) cycle
         nk = nk + 1
         kids(nk) = i
         j = nk
         do while (j > 1)
            if (rincl(kids(j-1)) >= rincl(kids(j))) exit
            tmp = kids(j-1)
            kids(j-1) = kids(j)
            kids(j) = tmp
            j = j - 1
         end do
      end do

      do j = 1, nk
         k = kids(j)
         ti = real(rincl(k), dp)/rate

         write(pct,'(f5.1)') 1d2*ti/total

         nchild = count(rparent(1:nr) == k)
         if (nchild == 0 .or. rincl(k) == 0) then
            rem = '------'
         else
            write(rem,'(f5.1)') 1d2*real(rexcl(k), dp)/real(rincl(k), dp)
            rem = trim(adjustl(rem))//'%'
         end if

         write(*,'(a,es8.2,a,a,a,a,1x,i0,1x,a)') &
            repeat('|   ', level)//'|-> ', ti, ' sec ', &
            trim(adjustl(pct)), '% ', trim(rem), rcalls(k), trim(rname(k))

         call print_tree(k, level+1, total, rate)
      end do
   end subroutine print_tree

end module profiler
