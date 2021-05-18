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

      module utils_dict
      use utils_def
      
      implicit none            

      contains
      
      
      recursive subroutine do_integer_dict_map(dict, fcn, ierr)
         type (integer_dict), pointer :: dict
         interface
            subroutine fcn(key, value, ierr)
               character (len=*), intent(in) :: key
               integer, intent(in) :: value
               integer, intent(out) :: ierr ! /= 0 means terminate map calls
            end subroutine fcn
         end interface
         type (integer_dict), pointer :: node
         integer :: ierr
         ierr = 0
         if (.not. associated(dict)) return
         node => dict
         do
            if (associated(node% left)) then
               call do_integer_dict_map(node% left, fcn, ierr)
               if (ierr /= 0) return
            end if
            call fcn(node% key, node% value, ierr)
            if (ierr /= 0) return            
            if (.not. associated(node% right)) return
            node => node% right
         end do
      end subroutine do_integer_dict_map
      
      
      subroutine do_get_dict_entries(dict, keys, values)
         type (integer_dict), pointer :: dict
         character (len=maxlen_key_string), pointer :: keys(:)
         integer, pointer :: values(:)
         
         integer :: cnt, ierr, sz
         sz = size_integer_dict(dict)
         sz = min(sz, size(keys,dim=1), size(values,dim=1))
         cnt = 0
         call do_integer_dict_map(dict, fcn, ierr)
         
         contains
         
         subroutine fcn(key, value, ierr)
            character (len=*), intent(in) :: key
            integer, intent(in) :: value
            integer, intent(out) :: ierr ! /= 0 means terminate map calls
            if (cnt >= sz) then
               ierr = -1
               return
            end if
            cnt = cnt+1
            keys(cnt) = key
            values(cnt) = value
         end subroutine fcn
         
      end subroutine do_get_dict_entries
      
      
      recursive subroutine show_key_entries(dict)
         type (integer_dict), pointer :: dict
         type (integer_dict), pointer :: node
         if (.not. associated(dict)) return
         node => dict
         do
            if (associated(node% left)) then
               call show_key_entries(node% left)
            end if
            write(*,fmt=*) trim(node% key), node% value
            if (.not. associated(node% right)) return
            node => node% right
         end do
      end subroutine show_key_entries
      
      
      subroutine find_key_entry(dict, key, node)
         type (integer_dict), pointer :: dict
         character (len=*), intent(in) :: key
         type (integer_dict), pointer :: node ! set null if cannot find key in dict
         type (hash_entry), pointer :: hash(:)
         integer :: i, hash_size, hashkey
         if (.not. associated(dict)) then
            nullify(node); return
         end if
         if (associated(dict% hash)) then
            hash => dict% hash
            hash_size = size(hash)
            hashkey = dict_hashkey(key, hash_size)
            do i=1, hash_size ! find an empty slot
               if (.not. associated(hash(hashkey)% ptr)) exit ! failed to find key
               if (hash(hashkey)% ptr % key == key) then
                  node => hash(hashkey)% ptr
                  return
               end if
               hashkey = hashkey+1
               if (hashkey > hash_size) hashkey = 1
            end do
            nullify(node)
            return ! failed to find key
         end if
         node => dict
         do
            if (node% key == key) return
            if (LLT(node% key,key)) then
               if (.not. associated(node% left)) then
                  nullify(node); return
               end if
               node => node% left
            else
               if (.not. associated(node% right)) then
                  nullify(node); return
               end if
               node => node% right
            end if
         end do      
      end subroutine find_key_entry
      
      
      recursive subroutine insert_node(node, root, duplicate)
         type (integer_dict), pointer :: node ! will be deallocated if a duplicate
         type (integer_dict), pointer :: root
         logical :: duplicate ! true if key was already defined
         
         integer :: height_left, height_right
         logical, parameter :: dbg = .false.

         if (dbg) write(*,*) 'insert ' // trim(node% key) // ' in ' // trim(root% key)
         
         if (node% key == root% key) then
            root% value = node% value
            deallocate(node)
            nullify(node)
            duplicate = .true.
            return
         end if
         
         if (LGT(node% key, root% key)) then ! insert on left
            if (.not. associated(root% left)) then
               root% left => node
            else
               call insert_node(node, root% left, duplicate)
            end if
            height_left = root% left% height
            height_right = height_of_right_branch(root)
            if (height_left - height_right == 2) then ! rebalance
               if (LGT(node% key, root% left% key)) then
                  call single_rotate_with_left(root)
               else
                  call double_rotate_with_left(root)
               end if
            end if
         else ! insert on right
            if (.not. associated(root% right)) then
               root% right => node
            else
               call insert_node(node, root% right, duplicate)
            end if
            height_right = root% right% height
            height_left = height_of_left_branch(root)
            if (height_right - height_left == 2) then ! rebalance
               if (LGT(root% right% key, node% key)) then
                  call single_rotate_with_right(root)
               else
                  call double_rotate_with_right(root)
               end if
            end if
         end if
         
         height_right = height_of_right_branch(root)
         height_left = height_of_left_branch(root)
         root% height = max(height_right, height_left) + 1
         
         if (dbg) write(*,*) 'new root is ' // trim(root% key)
         
         
         contains
         
         
         integer function height_of_left_branch(n)
            type (integer_dict), pointer :: n
            if (.not. associated(n% left)) then
               height_of_left_branch = 0
            else
               height_of_left_branch = n% left% height
            end if
         end function height_of_left_branch
         
         
         integer function height_of_right_branch(n)
            type (integer_dict), pointer :: n
            if (.not. associated(n% right)) then
               height_of_right_branch = 0
            else
               height_of_right_branch = n% right% height
            end if
         end function height_of_right_branch
         
         
         subroutine single_rotate_with_left(k2)
            type (integer_dict), pointer :: k2
            type (integer_dict), pointer :: k1
            k1 => k2% left
            if (.not. associated(k1% right)) then
               nullify(k2% left)
            else
               k2% left => k1% right
            end if
            k1% right => k2
            k2% height = max(height_of_left_branch(k2), height_of_right_branch(k2)) + 1
            k1% height = max(height_of_left_branch(k1), k2% height) + 1
            k2 => k1
         end subroutine single_rotate_with_left
         
         
         subroutine single_rotate_with_right(k1)
            type (integer_dict), pointer :: k1
            type (integer_dict), pointer :: k2
            k2 => k1% right
            if (.not. associated(k2% left)) then
               nullify(k1% right)
            else
               k1% right => k2% left
            end if
            k2% left => k1
            k1% height = max(height_of_right_branch(k1), height_of_left_branch(k1)) + 1
            k2% height = max(height_of_right_branch(k2), k1% height) + 1
            k1 => k2
         end subroutine single_rotate_with_right
         
         
         subroutine double_rotate_with_left(k)
            type (integer_dict), pointer :: k
            call single_rotate_with_right(k% left)
            call single_rotate_with_left(k)
         end subroutine double_rotate_with_left
         
         
         subroutine double_rotate_with_right(k)
            type (integer_dict), pointer :: k
            call single_rotate_with_left(k% right)
            call single_rotate_with_right(k)
         end subroutine double_rotate_with_right
         

      end subroutine insert_node
      
      
      subroutine do_integer_dict_define(dict, key, value, duplicate, ierr)
         type (integer_dict), pointer :: dict ! pass null for empty dict
         character (len=*), intent(in) :: key
         integer, intent(in) :: value
         logical, intent(out) :: duplicate ! true if key was already defined
         integer, intent(out) :: ierr
         type (integer_dict), pointer :: node
         logical, parameter :: dbg = .false.
         ierr = 0
         allocate(node, stat=ierr)
         if (ierr /= 0) return
!$omp critical (dict_define)
         duplicate = .false.
         node% key = key
         node% value = value
         node% height = 1
         nullify(node% left)
         nullify(node% right)
         nullify(node% hash)
         if (dbg) write(*,*) 'insert node ' // trim(key)
         if (.not. associated(dict)) then ! this is the 1st entry
            dict => node
         else
            if (associated(dict% hash)) then
               deallocate(dict% hash)
               nullify(dict% hash)
            end if
            call insert_node(node, dict, duplicate)
         end if
!$omp end critical (dict_define)
         if (dbg) then ! check tree
            write(*,*) 'done insert node ' // trim(key) // ' new root ' // trim(dict% key)
            write(*,*)
            call check_dict(dict, ierr)
            call show_key_entries(dict)
            write(*,*) 'done insert ' // trim(key)
         end if
      end subroutine do_integer_dict_define
      
      
      subroutine do_integer_dict_create_hash(dict, ierr)
         type (integer_dict), pointer :: dict
         integer, intent(out) :: ierr
         
         integer :: cnt, hash_size, i, collisions
         type (hash_entry), pointer :: hash(:)
         
         ierr = 0
         if (.not. associated(dict)) then
            ierr = -1; return
         end if
         if (associated(dict% hash)) return
         
!$omp critical (create_hash)
         if (.not. associated(dict% hash)) then
            cnt = size_integer_dict(dict) ! number of entries         
            if (cnt > 0) then
               hash_size = 4*cnt
               allocate(dict% hash(hash_size), stat=ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in allocate for create hash', hash_size
               else
                  hash => dict% hash
                  do i=1,hash_size
                     nullify(hash(i)% ptr)
                  end do
                  collisions = 0
                  call do_enter_hash(dict, hash, hash_size, collisions)
               end if
            end if
         end if
!$omp end critical (create_hash)
         
      end subroutine do_integer_dict_create_hash

      
      recursive subroutine check_dict(dict, ierr)
         type (integer_dict), pointer :: dict
         integer, intent(out) :: ierr
         integer :: height_left, height_right, height
         if (associated(dict% left)) then
            if (LGT(dict% key, dict% left% key)) then
               write(*,*) 'wrong order dict% key, dict% left% key ' // &
                           trim(dict% key) // ' ' // trim(dict% left% key)
               ierr = -1
               return
            end if
            call check_dict(dict% left, ierr)
            if (ierr /= 0) return
            height_left = dict% left% height
         else
            height_left = 0
         end if
         if (associated(dict% right)) then
            if (LGT(dict% right% key, dict% key)) then
               write(*,*) 'wrong order dict% right% key, dict% key ' // &
                           trim(dict% right% key) // ' ' // trim(dict% key)
               ierr = -1
               return
            end if
            call check_dict(dict% right, ierr)
            if (ierr /= 0) return
            height_right = dict% right% height
         else
            height_right = 0
         end if
         height = max(height_left, height_right) + 1
         if (dict% height /= height) then
            write(*,*) 'bad height for ' // trim(dict% key)
            ierr = -1
         end if
      end subroutine check_dict
      
      
      subroutine do_integer_dict_lookup(dict, key, value, ierr)
         type (integer_dict), pointer :: dict
         character (len=*), intent(in) :: key
         integer, intent(out) :: value
         integer, intent(out) :: ierr ! 0 if found key in dict, -1 if didn't
         type (integer_dict), pointer :: node
         logical, parameter :: dbg = .false.
         if (dbg) then
            call show_key_entries(dict)
            write(*,*)
            write(*,*) 'lookup key ' // trim(key)
            write(*,*)
         end if
         ierr = 0
         value = 0
         call do_integer_dict_create_hash(dict, ierr)
         if (ierr /= 0) return
         call find_key_entry(dict, key, node)
         if (associated(node)) then
            value = node% value
            return
         end if
         ierr = -1 
      end subroutine do_integer_dict_lookup
      
      
      recursive subroutine do_integer_dict_free(dict)
         type (integer_dict), pointer :: dict
         type (integer_dict), pointer :: node, next
         if (.not. associated(dict)) return
         node => dict
         dict => null()
         if (associated(node% hash)) deallocate(node% hash)
         do
            if (associated(node% left)) call do_integer_dict_free(node% left)
            if (.not. associated(node% right)) then
               deallocate(node)
               return
            end if
            next => node% right
            deallocate(node)
            node => next
         end do
      end subroutine do_integer_dict_free
      
      
      recursive function size_integer_dict(dict) result(cnt)
         type (integer_dict), pointer :: dict
         type (integer_dict), pointer :: node, next
         integer :: cnt
         cnt = 0
         if (.not. associated(dict)) return
         node => dict
         do
            cnt = cnt + 1
            if (associated(node% left)) cnt = cnt + size_integer_dict(node% left)
            if (.not. associated(node% right)) return
            next => node% right
            node => next
         end do
      end function size_integer_dict
      
      
      recursive subroutine do_enter_hash(dict, hash, hash_size, collisions)
         type (integer_dict), pointer :: dict
         type (hash_entry), pointer :: hash(:)
         integer, intent(in) :: hash_size
         integer, intent(inout) :: collisions
         type (integer_dict), pointer :: node, next
         integer :: hashkey, size, i
         logical :: okay
         if (.not. associated(dict)) return
         node => dict
         do
            ! enter node in hash
            hashkey = dict_hashkey(node% key, hash_size)
            okay = .false.
            do i=1, hash_size ! find an empty slot
               if (.not. associated(hash(hashkey)% ptr)) then
                  hash(hashkey)% ptr => node
                  okay = .true.
                  exit
               end if
               hashkey = hashkey+1
               collisions = collisions+1
               if (hashkey > hash_size) hashkey = 1
            end do
            if (.not. okay) then
               write(*,*) 'failed in do_enter_hash'
               error stop 1           
         end if
            if (associated(node% left)) &
               call do_enter_hash(node% left, hash, hash_size, collisions)
            if (.not. associated(node% right)) return
            next => node% right
            node => next
         end do
      end subroutine do_enter_hash
      
      
      integer function dict_hashkey(key, hash_size) ! value between 1 and hash_size
         character (len=*) :: key
         integer, intent(in) :: hash_size
         integer:: i, len, new, hash, c
         len = len_trim(key)
         if (len == 0) then
            dict_hashkey = 1
            return
         end if
         ! source: http://www.partow.net/programming/hashfunctions/#APHashFunction
         hash = -1431655766 ! Z'AAAAAAAA'
         do i = 1, len
            c = ichar(key(i:i))
            if (iand(c,1)==1) then
               !new = (hash <<  7) ^ (*str) * (hash >> 3)
               new = ieor(ishft(hash,7), c*ishft(hash,3))
            else
               !new = ~((hash << 11) + (*str) ^ (hash >> 5))
               new = not(ishft(hash,11) + ieor(c,ishft(hash,-5)))
            end if
            hash = ieor(hash, new)
         end do
         dict_hashkey = hash
         if (dict_hashkey < 0) then
            dict_hashkey = -dict_hashkey
         else if (dict_hashkey == 0) then
            dict_hashkey = 1
         end if
         dict_hashkey = 1 + mod(dict_hashkey-1, hash_size)
         if (dict_hashkey <= 0) then
            write(*,*) 'bad dict_hashkey for ' // trim(key), dict_hashkey
            stop 'dict_hashkey'
         end if
      end function dict_hashkey


      end module utils_dict

