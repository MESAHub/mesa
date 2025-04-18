! ***********************************************************************
!
!   Copyright (C) 2017 Josiah Schwab & The MESA Team
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

module weak_support

  implicit none

contains

  subroutine parse_weak_rate_name(name, lhs, rhs, ierr)
    use chem_def
    use chem_lib
    character (len=*), intent(in) :: name
    character (len=iso_name_length), intent(out) :: lhs, rhs
    integer, intent(out) :: ierr

    integer :: len, i, j, cid

    ierr = -1
    len = len_trim(name)
    if (name(1:2) /= 'r_') return
    i = 3

    ! get lhs isotope
    call nxt
    cid = chem_get_iso_id(name(i:j))
    if (cid == nuclide_not_found) then
       return
    else
       lhs = name(i:j)
    end if

    ! check middle is wk or wk-minus
    i=j+2
    call nxt
    if (.not. ((name(i:j) == 'wk') .or. (name(i:j) == 'wk-minus'))) then
       return
    end if

    ! get rhs isotope
    i=j+2
    call nxt
    cid = chem_get_iso_id(name(i:j))
    if (cid == nuclide_not_found) then
       return
    else
       rhs = name(i:j)
    end if

    ierr = 0

  contains

    ! calling nxt sets j to last char of token
    subroutine nxt
      j = i
      do
         if (j >= len) return
         j = j+1
         if (name(j:j) == '_') then
            j = j-1; return
         end if
      end do
    end subroutine nxt

  end subroutine parse_weak_rate_name

end module weak_support

