! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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
         character (len=maxlen_key_string) :: key
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
