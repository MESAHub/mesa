! ***********************************************************************
!
!   Copyright (C) 2017 Josiah Schwab, Bill Paxton
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

