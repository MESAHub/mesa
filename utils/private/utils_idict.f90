! ***********************************************************************
!
!   Copyright (C) 2014  The MESA Team
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

      module utils_idict
      use utils_def
      
      implicit none            

      contains
      
      
      recursive subroutine do_integer_idict_map(idict, fcn, ierr)
         type (integer_idict), pointer :: idict
         interface
            subroutine fcn(key1, key2, value, ierr)
               integer, intent(in) :: key1, key2, value
               integer, intent(out) :: ierr ! /= 0 means terminate map calls
            end subroutine fcn
         end interface
         type (integer_idict), pointer :: node
         integer :: ierr
         ierr = 0
         if (.not. associated(idict)) return
         node => idict
         do
            if (associated(node% left)) then
               call do_integer_idict_map(node% left, fcn, ierr)
               if (ierr /= 0) return
            end if
            call fcn(node% key1, node% key2, node% value, ierr)
            if (ierr /= 0) return            
            if (.not. associated(node% right)) return
            node => node% right
         end do
      end subroutine do_integer_idict_map
      
      
      subroutine do_get_idict_entries(idict, key1s, key2s, values)
         type (integer_idict), pointer :: idict
         integer, pointer, dimension(:) :: key1s, key2s, values
         
         integer :: cnt, ierr, sz
         sz = size_integer_idict(idict)
         sz = min(sz, size(key1s,dim=1), size(key2s,dim=1), size(values,dim=1))
         cnt = 0
         call do_integer_idict_map(idict, fcn, ierr)
         
         contains
         
         subroutine fcn(key1, key2, value, ierr)
            integer, intent(in) :: key1, key2, value
            integer, intent(out) :: ierr ! /= 0 means terminate map calls
            if (cnt >= sz) then
               ierr = -1
               return
            end if
            cnt = cnt+1
            key1s(cnt) = key1
            key2s(cnt) = key2
            values(cnt) = value
         end subroutine fcn
         
      end subroutine do_get_idict_entries
      
      
      recursive subroutine show_key1_key2_entries(idict)
         type (integer_idict), pointer :: idict
         type (integer_idict), pointer :: node
         if (.not. associated(idict)) return
         node => idict
         do
            if (associated(node% left)) then
               call show_key1_key2_entries(node% left)
            end if
            write(*,fmt='(3i10)') node% key1, node% key2, node% value
            if (.not. associated(node% right)) return
            node => node% right
         end do
      end subroutine show_key1_key2_entries
      
      
      subroutine find_key1_key2_entry(idict, key1, key2, node)
         type (integer_idict), pointer :: idict
         integer, intent(in) :: key1, key2
         type (integer_idict), pointer :: node ! set null if cannot find key1, key2 in idict
         type (ihash_entry), pointer :: hash(:)
         integer :: i, hash_size, hashkey
         if (.not. associated(idict)) then
            nullify(node); return
         end if
         if (associated(idict% hash)) then
            hash => idict% hash
            hash_size = size(hash)
            hashkey = idict_hashkey(key1, key2, hash_size)
            do i=1, hash_size ! find an empty slot
               if (.not. associated(hash(hashkey)% ptr)) exit
               if (hash(hashkey)% ptr% key1 == key1 .and. &
                   hash(hashkey)% ptr% key2 == key2) then
                  node => hash(hashkey)% ptr
                  return
               end if
               hashkey = hashkey+1
               if (hashkey > hash_size) hashkey = 1
            end do
            nullify(node)
            return ! failed to find key1, key2
         end if
         node => idict
         do
            if (node% key1 == key1 .and. node% key2 == key2) return
            if (node% key1 < key1 .or. &
                  (node% key1 == key1 .and. node% key2 < key2)) then
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
      end subroutine find_key1_key2_entry
      
      
      recursive subroutine insert_node(node, root, duplicate)
         type (integer_idict), pointer :: node ! will be deallocated if a duplicate
         type (integer_idict), pointer :: root
         logical :: duplicate ! true if key was already defined
         
         integer :: height_left, height_right
         logical, parameter :: dbg = .false.
         
         if (node% key1 == root% key1 .and. node% key2 == root% key2) then
            root% value = node% value
            deallocate(node)
            nullify(node)
            duplicate = .true.
            return
         end if
         
         if (node% key1 > root% key1 .or. &
            (node% key1 == root% key1 .and. &
             node% key2 > root% key2)) then ! insert on left
            if (.not. associated(root% left)) then
               root% left => node
            else
               call insert_node(node, root% left, duplicate)
            end if
            height_left = root% left% height
            height_right = height_of_right_branch(root)
            if (height_left - height_right == 2) then ! rebalance
               if (node% key1 > root% left% key1 .or. &
                  (node% key1 == root% left% key1 .and. &
                   node% key2 > root% left% key2)) then ! insert on left
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
               if (root% right% key1 > node% key1 .or. &
                  (root% right% key1 == node% key1 .and. &
                   root% right% key2 > node% key2)) then
                  call single_rotate_with_right(root)
               else
                  call double_rotate_with_right(root)
               end if
            end if
         end if
         
         height_right = height_of_right_branch(root)
         height_left = height_of_left_branch(root)
         root% height = max(height_right, height_left) + 1
                  
         
         contains
         
         
         integer function height_of_left_branch(n)
            type (integer_idict), pointer :: n
            if (.not. associated(n% left)) then
               height_of_left_branch = 0
            else
               height_of_left_branch = n% left% height
            end if
         end function height_of_left_branch
         
         
         integer function height_of_right_branch(n)
            type (integer_idict), pointer :: n
            if (.not. associated(n% right)) then
               height_of_right_branch = 0
            else
               height_of_right_branch = n% right% height
            end if
         end function height_of_right_branch
         
         
         subroutine single_rotate_with_left(k2)
            type (integer_idict), pointer :: k2
            type (integer_idict), pointer :: k1
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
            type (integer_idict), pointer :: k1
            type (integer_idict), pointer :: k2
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
            type (integer_idict), pointer :: k
            call single_rotate_with_right(k% left)
            call single_rotate_with_left(k)
         end subroutine double_rotate_with_left
         
         
         subroutine double_rotate_with_right(k)
            type (integer_idict), pointer :: k
            call single_rotate_with_left(k% right)
            call single_rotate_with_right(k)
         end subroutine double_rotate_with_right
         

      end subroutine insert_node
      
      
      subroutine do_integer_idict_define(idict, key1, key2, value, duplicate, ierr)
         type (integer_idict), pointer :: idict ! pass null for empty idict
         integer, intent(in) :: key1, key2, value
         logical, intent(out) :: duplicate ! true if key was already defined
         integer, intent(out) :: ierr
         type (integer_idict), pointer :: node
         logical, parameter :: dbg = .false.
         ierr = 0
         allocate(node, stat=ierr)
         if (ierr /= 0) return
!$omp critical (idict_define)
         duplicate = .false.
         node% key1 = key1
         node% key2 = key2
         node% value = value
         node% height = 1
         nullify(node% left)
         nullify(node% right)
         nullify(node% hash)
         if (.not. associated(idict)) then ! this is the 1st entry
            idict => node
         else
            if (associated(idict% hash)) then
               deallocate(idict% hash)
               nullify(idict% hash)
            end if
            call insert_node(node, idict, duplicate)
         end if
!$omp end critical (idict_define)
         if (dbg) then ! check tree
            write(*,*)
            call check_idict(idict, ierr)
            call show_key1_key2_entries(idict)
            write(*,*) 'done insert', key1, key2
         end if
      end subroutine do_integer_idict_define
      
      
      subroutine do_integer_idict_create_hash(idict, ierr)
         type (integer_idict), pointer :: idict
         integer, intent(out) :: ierr
         
         integer :: cnt, hash_size, i, collisions
         type (ihash_entry), pointer :: hash(:)
         
         ierr = 0
         if (.not. associated(idict)) then
            ierr = -1; return
         end if
         if (associated(idict% hash)) return
         
!$omp critical (create_hash)
         if (.not. associated(idict% hash)) then
            cnt = size_integer_idict(idict) ! number of entries         
            if (cnt > 0) then
               hash_size = 4*cnt
               allocate(idict% hash(hash_size), stat=ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in allocate for create hash', hash_size
               else
                  hash => idict% hash
                  do i=1,hash_size
                     nullify(hash(i)% ptr)
                  end do
                  collisions = 0
                  call do_enter_hash(idict, hash, hash_size, collisions)
               end if
            end if
         end if
!$omp end critical (create_hash)
         
      end subroutine do_integer_idict_create_hash

      
      recursive subroutine check_idict(idict, ierr)
         type (integer_idict), pointer :: idict
         integer, intent(out) :: ierr
         integer :: height_left, height_right, height
         if (associated(idict% left)) then
            if (idict% key1 > idict% left% key1 .or. &
               (idict% key1 == idict% left% key1 .and. &
                idict% key2 > idict% left% key2)) then
               write(*,*) 'wrong order idict% key1, key2, idict% left% key1, key2', &
                   idict% key1, idict% key2, idict% left% key1, idict% left% key2
               ierr = -1
               return
            end if
            call check_idict(idict% left, ierr)
            if (ierr /= 0) return
            height_left = idict% left% height
         else
            height_left = 0
         end if
         if (associated(idict% right)) then
            if (idict% right% key1 > idict% key1 .or. &
               (idict% right% key1 == idict% key1 .and. &
                idict% right% key2 > idict% key2)) then
               write(*,*) 'wrong order idict% right% key1, key2, idict% key1, key2', &
                   idict% right% key1, idict% right% key2, idict% key1, idict% key2
               ierr = -1
               return
            end if
            call check_idict(idict% right, ierr)
            if (ierr /= 0) return
            height_right = idict% right% height
         else
            height_right = 0
         end if
         height = max(height_left, height_right) + 1
         if (idict% height /= height) then
            write(*,*) 'bad height for', idict% key1, idict% key2
            ierr = -1
         end if
      end subroutine check_idict
      
      
      subroutine do_integer_idict_lookup(idict, key1, key2, value, ierr)
         type (integer_idict), pointer :: idict
         integer, intent(in) :: key1, key2
         integer, intent(out) :: value
         integer, intent(out) :: ierr ! 0 if found key1, key2 in idict, -1 if didn't
         type (integer_idict), pointer :: node
         logical, parameter :: dbg = .false.
         if (dbg) then
            call show_key1_key2_entries(idict)
            write(*,*)
            write(*,*) 'lookup key1, key2', key1, key2
            write(*,*)
         end if
         ierr = 0
         value = 0
         call do_integer_idict_create_hash(idict, ierr)
         if (ierr /= 0) return
         call find_key1_key2_entry(idict, key1, key2, node)
         if (associated(node)) then
            value = node% value
            return
         end if
         ierr = -1 
      end subroutine do_integer_idict_lookup
      
      
      recursive subroutine do_integer_idict_free(idict)
         type (integer_idict), pointer :: idict
         type (integer_idict), pointer :: node, next
         if (.not. associated(idict)) return
         node => idict
         if (associated(node% hash)) deallocate(node% hash)
         do
            if (associated(node% left)) call do_integer_idict_free(node% left)
            if (.not. associated(node% right)) then
               deallocate(node)
               return
            end if
            next => node% right
            deallocate(node)
            node => next
         end do
      end subroutine do_integer_idict_free
      
      
      recursive function size_integer_idict(idict) result(cnt)
         type (integer_idict), pointer :: idict
         type (integer_idict), pointer :: node, next
         integer :: cnt
         cnt = 0
         if (.not. associated(idict)) return
         node => idict
         do
            cnt = cnt + 1
            if (associated(node% left)) cnt = cnt + size_integer_idict(node% left)
            if (.not. associated(node% right)) return
            next => node% right
            node => next
         end do
      end function size_integer_idict
      
      
      recursive subroutine do_enter_hash(idict, hash, hash_size, collisions)
         type (integer_idict), pointer :: idict
         type (ihash_entry), pointer :: hash(:)
         integer, intent(in) :: hash_size
         integer, intent(inout) :: collisions
         type (integer_idict), pointer :: node, next
         integer :: hashkey, size, i
         logical :: okay
         if (.not. associated(idict)) return
         node => idict
         do
            ! enter node in hash
            hashkey = idict_hashkey(node% key1, node% key2, hash_size)
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
      
      
      integer function idict_hashkey(key1, key2, hash_size) ! value between 1 and hash_size
         integer, intent(in) :: key1, key2, hash_size
         integer:: new, hash, c
         ! source: http://www.partow.net/programming/hashfunctions/#APHashFunction
         hash = -1431655766 ! Z'AAAAAAAA'
         c = key1
         if (iand(c,1)==1) then
            !new = (hash <<  7) ^ (*str) * (hash >> 3)
            new = ieor(ishft(hash,7), c*ishft(hash,3))
         else
            !new = ~((hash << 11) + (*str) ^ (hash >> 5))
            new = not(ishft(hash,11) + ieor(c,ishft(hash,-5)))
         end if
         hash = ieor(hash, new)
         c = key2
         if (iand(c,1)==1) then
            !new = (hash <<  7) ^ (*str) * (hash >> 3)
            new = ieor(ishft(hash,7), c*ishft(hash,3))
         else
            !new = ~((hash << 11) + (*str) ^ (hash >> 5))
            new = not(ishft(hash,11) + ieor(c,ishft(hash,-5)))
         end if
         idict_hashkey = ieor(hash, new)
         if (idict_hashkey < 0) then
            idict_hashkey = -idict_hashkey
         else if (idict_hashkey == 0) then
            idict_hashkey = 1
         end if
         idict_hashkey = 1 + mod(idict_hashkey-1, hash_size)
         if (idict_hashkey <= 0) then
            write(*,*) 'bad idict_hashkey for', key1, key2, idict_hashkey
            stop 'idict_hashkey'
         end if
      end function idict_hashkey


      end module utils_idict

