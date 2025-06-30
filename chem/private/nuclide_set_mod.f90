! ***********************************************************************
!
!   Copyright (C) 2010  Ed Brown & The MESA Team
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

      module nuclide_set_mod

      use chem_def
      use const_def, only: dp

      implicit none

      contains

      function rank_in_set(iso, set)
         character(len=iso_name_length), intent(in) :: iso
         type(nuclide_set), dimension(:), intent(in) :: set
         integer :: rank_in_set
         integer :: low, high, mid, i
         integer, parameter :: max_cycles = 20
         low = 1
         high = size(set)
         if (adjustl(iso) < adjustl(set(low)% nuclide) .or. adjustl(iso) > adjustl(set(high)% nuclide)) then
            rank_in_set = nuclide_not_found
            return
         end if
         do i = 1, max_cycles
            if (high-low <=1) then
               if (adjustl(iso) == adjustl(set(high)% nuclide)) then
                  rank_in_set = set(high)% rank
               else if (adjustl(iso) == adjustl(set(low)% nuclide)) then
                  rank_in_set = set(low)% rank
               else
                  rank_in_set = nuclide_not_found
               end if
               return
            end if
            mid = (high+low)/2
            if (adjustl(iso) <= adjustl(set(mid)% nuclide)) then
               high = mid
            else if (adjustl(iso) > adjustl(set(mid)% nuclide)) then
               low = mid
            end if
         end do
         rank_in_set = nuclide_not_found
      end function rank_in_set


      subroutine sort_nuclide_set(set)
         type(nuclide_set), dimension(:), intent(inout) :: set
         integer :: n, i, ir, j, l
         type(nuclide_set) :: ts

         n = size(set)
         if (size(set) < 2) return
         l = n/2+1
         ir = n
         do
            if (l > 1) then
               l = l-1
               ts = set(l)
            else
               ts = set(ir)
               set(ir) = set(1)
               ir = ir-1
               if (ir == 1) then
                  set(1) = ts
                  return
               end if
            end if
            i = l
            j = l+l
            do
               if (j > ir) exit
               if (j < ir) then
                  if (compare_lt(set(j), set(j+1))) j = j+1
               end if
               if (compare_lt(ts, set(j))) then
                  set(i) = set(j)
                  i = j
                  j = j+j
               else
                  j = ir + 1
               end if
            end do
            set(i) = ts
         end do

         contains

         logical function compare_lt(a, b)
            type(nuclide_set), intent(in) :: a, b
            compare_lt = (adjustl(a% nuclide) < adjustl(b% nuclide))
         end function compare_lt

      end subroutine sort_nuclide_set

      end module nuclide_set_mod
