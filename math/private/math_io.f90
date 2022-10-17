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

module math_io
   
   ! Uses
   
   use const_lib, only : dp
   
   ! No implicit typing
   
   implicit none
   
   ! Access specifiers
   
   private
   
   public :: str_to_vector
   public :: str_to_double
   public :: double_to_str
   
   ! Procedures

contains
   
   subroutine str_to_vector (str, vec, n, ierr)
      
      character (len = *), intent(in) :: str
      real(dp), pointer, intent(inout) :: vec(:)
      integer, intent(out) :: n
      integer, intent(out) :: ierr
      
      integer :: maxlen
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      
      maxlen = size(vec, dim = 1)
      
      n = 0
      ierr = 0
      
      l = len_trim(str)
      
      j = 1
      
      do i = 1, maxlen
         
         do while (j < l .and. str(j:j) == ' ')
            j = j + 1
         end do
         
         k = j
         
         do while (k < l .and. str(k:k) /= ' ')
            k = k + 1
         end do
         
         call str_to_double(str(j:k), vec(i), ierr)
         if (ierr /= 0) then
            return
         end if
         
         n = i
         
         if (k == l) exit
         
         j = k + 1
      
      end do
   
   end subroutine str_to_vector
   
   !****
   
   subroutine str_to_double (str, x, ierr)
      
      character(*), intent(in) :: str
      real(dp), intent(out) :: x
      integer, intent(out) :: ierr
      
      read(str, *, ROUND = 'COMPATIBLE', IOSTAT = ierr) x
   
   end subroutine str_to_double
   
   !****
   
   subroutine double_to_str (x, str)
      
      integer, parameter :: l = 26
      
      real(dp), intent(in) :: x
      character(len = l), intent(out) :: str
      
      write(str, 100, ROUND = 'COMPATIBLE') x
      100 format(1PD26.16)
   
   end subroutine double_to_str

end module math_io
