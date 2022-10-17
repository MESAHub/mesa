! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

module utils_def
   implicit none
   
   integer, parameter :: min_io_unit = 29
   integer, parameter :: max_io_unit = 99
   
   integer, parameter :: eof_token = 1
   integer, parameter :: string_token = 2
   integer, parameter :: name_token = 3
   integer, parameter :: left_paren_token = 4
   integer, parameter :: right_paren_token = 5
   integer, parameter :: comma_token = 6
   
   integer, parameter :: maxlen_key_string = 50
   
   ! see http://en.wikipedia.org/wiki/AVL_tree
   type integer_dict
      character (len = maxlen_key_string) :: key
      integer :: value
      integer :: height
      type (integer_dict), pointer :: left
      type (integer_dict), pointer :: right
      type (hash_entry), pointer :: hash(:)
   end type integer_dict
   
   type hash_entry
      type (integer_dict), pointer :: ptr
   end type hash_entry
   
   type integer_idict
      integer :: key1, key2, value
      integer :: height
      type (integer_idict), pointer :: left
      type (integer_idict), pointer :: right
      type (ihash_entry), pointer :: hash(:)
   end type integer_idict
   
   type ihash_entry
      type (integer_idict), pointer :: ptr
   end type ihash_entry

end module utils_def

