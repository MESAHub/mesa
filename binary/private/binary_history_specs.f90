! ***********************************************************************
!
!   Copyright (C) 2010-2019  Pablo Marchant & The MESA Team
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


      module binary_history_specs

      use const_def
      use star_lib
      use star_def
      use math_lib
      use binary_def
      use binary_private_def

      implicit none

      logical, parameter :: open_close_log = .true.

      contains

      recursive subroutine add_binary_history_columns( &
            b, level, capacity, spec, history_columns_file, ierr)
         use utils_lib
         use utils_def
         use const_def, only: mesa_dir
         type (binary_info), pointer :: b
         integer, intent(in) :: level
         integer, intent(inout) :: capacity
         integer, pointer :: spec(:)
         character (len=*), intent(in) :: history_columns_file
         integer, intent(out) :: ierr

         integer :: iounit, n, i, t, id, j, cnt, ii, nxt_spec
         character (len=256) :: buffer, string, filename
         integer, parameter :: max_level = 20
         logical :: bad_item
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         if (level > max_level) then
            write(*,*) 'too many levels of nesting for binary log column files', level
            ierr = -1
            return
         end if

         ierr = 0

         ! first try local directory
         filename = history_columns_file
         if (len_trim(filename) == 0) filename = 'binary_history_columns.list'
         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then ! if don't find that file, look in binary/defaults
            filename = trim(mesa_dir) // '/binary/defaults/' // trim(filename)
            ierr = 0
            open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then ! fail
               write(*,*) 'failed to open ' // trim(history_columns_file)
               return                  
            end if
         end if
         
         if (dbg) then
            write(*,'(A)')
            write(*,*) 'binary_history_columns_file <' // trim(filename) // '>'
            write(*,'(A)')
         end if

         call count_specs
         
         n = 0
         i = 0
         bad_item = .false.
         
         do
         
            t = token(iounit, n, i, buffer, string)
            if (t == eof_token) exit
            if (t /= name_token) then
               call error; return
            end if
            
            if (string == 'include') then
               t = token(iounit, n, i, buffer, string)
               if (t /= string_token) then
                  call error; return
               end if
               call add_binary_history_columns(b, level+1, capacity, spec, string, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed for included log columns list ' // trim(string)
                  bad_item = .true.
               end if
               if (.not. bad_item) call count_specs
               cycle
            end if
            
            nxt_spec = do1_binary_history_spec(iounit, t, n, i, string, buffer, ierr)
            if (ierr /= 0) bad_item = .true.
            if (.not. bad_item) then
                call insert_spec(nxt_spec, string, ierr)
            end if
                    
         end do
         
         if (dbg) write(*,*) 'finished ' // trim(filename)
         
         close(iounit)
         
         if (bad_item) then
            ierr = -1
            return
         end if
         
         if (dbg) then
            write(*,'(A)')
            write(*,*) 'done add_binary_history_columns ' // trim(filename)
            write(*,'(A)')
         end if
         
         
         contains


         subroutine count_specs
            integer :: i
            j = 1
            do i=1, capacity
               if (spec(i) == 0) then
                  j = i; exit
               end if
            end do
         end subroutine count_specs
         
         
         subroutine make_room(ierr)
            integer, intent(out) :: ierr
            if (j < capacity) return
            capacity = 50 + (3*capacity)/2
            call realloc_integer(spec,capacity,ierr)
            spec(j+1:capacity) = 0
         end subroutine make_room
  
            
         subroutine insert_spec(c, name, ierr)
            integer, intent(in) :: c
            character (len=*) :: name
            integer, intent(out) :: ierr
            integer :: i
            include 'formats'
            do i=1,j-1
               if (spec(i) == c) return
            end do
            call make_room(ierr)
            if (ierr /= 0) return
            spec(j) = c
            if (dbg) write(*,2) trim(name), spec(j)
            j = j+1
         end subroutine insert_spec
         
         
         subroutine error
            ierr = -1
            close(iounit)
         end subroutine error
                  
         
      end subroutine add_binary_history_columns


      integer function do1_binary_history_spec( &
             iounit, t, n, i, string, buffer, ierr) result(spec)
          use utils_lib
          use utils_def
          use chem_lib

          integer :: iounit, t, n, i, j, id
          character (len=*) :: string, buffer
          integer, intent(out) :: ierr

          ierr = 0
          spec = -1

          do j=1, bh_col_id_max
             if (binary_history_column_name(j) == string) then
                 spec = j
                 return
             end if
          end do

          write(*,*) 'bad history list name: ' // trim(string)
          ierr = -1

      end function do1_binary_history_spec

      subroutine set_binary_history_columns(b, history_columns_file, ierr)
         use utils_lib, only: realloc_integer
         type(binary_info), pointer :: b
         character (len=*), intent(in) :: history_columns_file
         integer, intent(out) :: ierr
         integer :: capacity, cnt, i
         logical, parameter :: dbg = .false.
         integer, pointer :: old_history_column_spec(:) => null()
         character (len=strlen) :: fname
         logical :: history_file_exists
         if (dbg) write(*,*) 'set_binary_history_columns'
         ierr = 0
         old_history_column_spec => null()
         if (associated(b% history_column_spec)) old_history_column_spec => b% history_column_spec
         nullify(b% history_column_spec)
         capacity = 100 ! will increase if needed
         allocate(b% history_column_spec(capacity), stat=ierr)
         if (ierr /= 0) return
         b% history_column_spec(:) = 0
         call add_binary_history_columns( &
               b, 1, capacity, b% history_column_spec, history_columns_file, ierr)
         if (ierr /= 0) then
            if (associated(old_history_column_spec)) deallocate(old_history_column_spec)
            return
         end if
         ! delete trailing 0's
         cnt = capacity+1
         do i=1, capacity
            if (b% history_column_spec(i) == 0) then
               cnt = i; exit
            end if
         end do
         capacity = cnt-1
         call realloc_integer(b% history_column_spec, capacity, ierr)
         if (ierr /= 0) return
         if (associated(old_history_column_spec)) then
            ! check that haven't changed the cols specs for an existing log file
            if (ierr /= 0) return
            fname = trim(b% log_directory) // '/' // trim(b% history_name)
            inquire(file=trim(fname), exist=history_file_exists)
            if (history_file_exists) then
               if (capacity /= size(old_history_column_spec)) then
                  ierr = -1
                  write(*,*) 'new size of log col specs', capacity
                  write(*,*) 'old size of log col specs', size(old_history_column_spec)
               else
                  do i=1,capacity
                     if (old_history_column_spec(i) /= b% history_column_spec(i)) then
                        write(*,*) 'change in log col spec', i, old_history_column_spec(i), b% history_column_spec(i)
                        ierr = -1
                     end if
                  end do
               end if
               if (ierr /= 0) then
                  write(*,*) 'ERROR: cannot change binary log columns when have an existing log file'
                  write(*,*) 'please delete the log file or go back to previous log columns list'
               end if
            end if
            deallocate(old_history_column_spec)
            if (ierr /= 0) return
         end if
         if (associated(b% history_column_is_integer)) deallocate(b% history_column_is_integer)
         if (associated(b% history_column_values)) deallocate(b% history_column_values)
         allocate(b% history_column_is_integer(capacity), b% history_column_values(capacity), stat=ierr)
         if (ierr /= 0) return
         b% history_column_is_integer(:) = .false.
         b% history_column_values(:) = 0
         !b% need_to_write_history = .false.
         if (dbg) write(*,*) 'binary num log columns', capacity
         if (dbg) call mesa_error(__FILE__,__LINE__,'debug: set_binary_history_columns')
      end subroutine set_binary_history_columns


      end module binary_history_specs
