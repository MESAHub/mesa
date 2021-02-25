! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module load_weak
      
      use rates_def
      use const_def, only: dp
      use utils_lib, only: mesa_error
      use weaklib_tables, only: weaklib_rate_table
      use suzuki_tables, only: private_load_suzuki_tables

      implicit none
      
      private :: private_load_weak_tables

      contains
      
      
      subroutine load_weak_data(ierr)
         integer, intent(out) :: ierr         
         ierr = 0         
         call private_load_weak_tables(ierr)
         if (ierr /= 0) return

         call load_user_weak_tables(ierr)
         if (ierr /= 0) return
         if (use_suzuki_tables) then
            call private_load_suzuki_tables(ierr)
            if (ierr /= 0) return
         end if

         call load_weak_info_list(ierr)
      end subroutine load_weak_data
      
      
      subroutine load_weak_info_list(ierr)
         use utils_lib
         use math_lib, only: str_to_vector
         integer, intent(out) :: ierr
         
         integer :: iounit, i, nvec
         character (len=256) :: filename, string
         character(len=iso_name_length) :: lhs, rhs
         character(len=2*iso_name_length+1) :: key
         real(dp), target :: vec_ary(2)
         real(dp), pointer :: vec(:)
         integer, parameter :: max_num_weak_info = 1000

         logical, parameter :: dbg = .false.

         include 'formats'
         
         ierr = 0
         vec => vec_ary
         
         filename = trim(weak_data_dir) // '/weak_info.list'
         ierr = 0
         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            return
         end if
         
         if (dbg) then
            write(*,*)
            write(*,*) 'weak info filename <' // trim(filename) // '>'
            write(*,*)
         end if

         do ! skip to line starting with 'from '
            read(iounit,'(a)',iostat=ierr) string
            if (failed('read weak info comments')) return
            if (len_trim(string) > 4) then
               if (string(1:5) == 'from ') exit
            end if
         end do

         nullify(weak_info_list_dict)
         allocate(weak_info_list_halflife(max_num_weak_info))
         allocate(weak_info_list_Qneu(max_num_weak_info))
         num_weak_info_list_reactions = 0
         do i = 1, max_num_weak_info ! keep reading until end of file         
            read(iounit,fmt='(a5,a5,a)',iostat=ierr) lhs, rhs, string
            if (ierr == 0) then
               call str_to_vector(string, vec, nvec, ierr)
               if (nvec < 2) ierr = -1
            end if
            if (ierr /= 0) then
               ierr = 0; exit
            end if
            weak_info_list_halflife(i) = vec(1)
            weak_info_list_Qneu(i) = vec(2)    
            call create_weak_dict_key(lhs, rhs, key)
            !write(*,'(a)') 'weak info list key ' // trim(key)
            call integer_dict_define(weak_info_list_dict, key, i, ierr)
            if (failed('integer_dict_define')) return
            num_weak_info_list_reactions = i
         end do
         
         close(iounit)
         
         if (num_weak_info_list_reactions == 0) then
            ierr = -1
            write(*,*) 'failed trying to read weak_info.list -- no reactions?'
            return
         end if
         
         if (num_weak_info_list_reactions == max_num_weak_info) then
            ierr = -1
            write(*,*) 'failed trying to read weak_info.list -- too many reactions?'
            return
         end if
         
         call integer_dict_create_hash(weak_info_list_dict, ierr)
         if (ierr /= 0) return
         
         call realloc_double(weak_info_list_halflife, num_weak_info_list_reactions, ierr)
         if (ierr /= 0) return
         
         call realloc_double(weak_info_list_Qneu, num_weak_info_list_reactions, ierr)
         if (ierr /= 0) return
         
         
         contains
         
         logical function failed(str)
            character (len=*) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*,*) 'failed: ' // trim(str)
            end if
         end function failed
         
         
      end subroutine load_weak_info_list
      
      
      subroutine private_load_weak_tables(ierr)
         use utils_lib
         use chem_lib, only: chem_get_iso_id
         use chem_def, only: iso_name_length
         integer, intent(out) :: ierr
         
         integer :: iounit, i, ios, id
         character (len=256) :: filename, cache_filename, string
         character(len=iso_name_length) :: lhs1, rhs1, lhs2, rhs2, weak_lhs, weak_rhs
         character(len=2*iso_name_length+1) :: key

         integer, parameter :: weak_num_T9 = 12, weak_num_lYeRho = 11
         real(dp) :: weak_reaction_T9s(weak_num_T9) = &
            (/ 0.01d0, 0.1d0, 0.2d0, 0.4d0, 0.7d0, 1.0d0, 1.5d0, 2.0d0, 3.0d0, 5.0d0, 10.0d0, 30.0d0 /)
         real(dp) :: weak_reaction_lYeRhos(weak_num_lYeRho) = &
            (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0, 10.0d0, 11.0d0 /)

         integer, parameter :: i_ldecay = 1, i_lcapture = 2, i_lneutrino = 3

         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         
         ios = -1
         if (rates_use_cache) then
            cache_filename = trim(rates_cache_dir) // '/weakreactions.bin'
            ios = 0
            open(newunit=iounit,file=trim(cache_filename),action='read', &
               status='old',iostat=ios,form='unformatted')
            if (ios == 0) then ! opened it okay
               call read_weak_cache(iounit,ios)
               close(iounit)
            end if
         end if
         
         if (ios /= 0) then ! need to read data file
         
            filename = trim(weak_data_dir) // '/weakreactions.tables'
            ierr = 0
            open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(filename)
               return
            end if
         
            if (dbg) then
               write(*,*)
               write(*,*) 'weaklib filename <' // trim(filename) // '>'
               write(*,*)
            end if

            do ! skip to after line starting with '='
               read(iounit,'(a)',iostat=ierr) string
               if (failed('read header')) return
               if (len_trim(string) > 0) then
                  if (string(1:1) == '=') exit
               end if
            end do

            if (.not. skip_line()) return
         
            read(iounit,*,iostat=ierr) num_weak_reactions
            if (failed('read num_weak_reactions')) return
         
            if (dbg) write(*,2) 'num_weak_reactions', num_weak_reactions
         
            call alloc
            if (failed('allocate')) return
         
            do i = 1, num_weak_reactions
               if (.not. skip_line()) return
               if (mod(i,2)==1) then ! first of pair
                  if (.not. skip_line()) return
                  if (.not. skip_line()) return
                  read(iounit,fmt='(2a5)',iostat=ierr) lhs1, rhs1
                  if (failed('read lhs1, rhs1')) return
                  if (lhs1 == 'al-6') lhs1 = 'al26-1'
                  if (rhs1 == 'al-6') rhs1 = 'al26-1'
                  if (lhs1 == 'al*6') lhs1 = 'al26-2'
                  if (rhs1 == 'al*6') rhs1 = 'al26-2'
                  read(iounit,fmt='(2a5)',iostat=ierr) lhs2, rhs2
                  if (failed('read lhs2, rhs2')) return
                  if (lhs2 == 'al-6') lhs2 = 'al26-1'
                  if (rhs2 == 'al-6') rhs2 = 'al26-1'
                  if (lhs2 == 'al*6') lhs2 = 'al26-2'
                  if (rhs2 == 'al*6') rhs2 = 'al26-2'
                  if (.not. skip_line()) return
                  if (.not. skip_line()) return
                  weak_lhs = lhs1
                  weak_rhs = rhs1
               else
                  weak_lhs = lhs2
                  weak_rhs = rhs2
               end if
               call adjust_name(weak_lhs)
               call adjust_name(weak_rhs)
               id = chem_get_iso_id(weak_lhs)
               if (id <= 0) then
                  write(*,*) 'weaklib FATAL ERROR: unknown nuclide ' // weak_lhs
                  call mesa_error(__FILE__,__LINE__)
               end if
               weak_lhs_nuclide_id(i) = id
               id = chem_get_iso_id(weak_rhs)
               if (id <= 0) then
                  write(*,*) 'weaklib FATAL ERROR: unknown nuclide ' // weak_rhs
                  call mesa_error(__FILE__,__LINE__)
               end if
               weak_reaclib_id(i) = 0
               weak_rhs_nuclide_id(i) = id
               weak_lhs_nuclide_name(i) = weak_lhs
               weak_rhs_nuclide_name(i) = weak_rhs
               if (.not. skip_line()) return
               call read_table(i,i_ldecay)
               if (failed('read ldecay')) return
               if (.not. skip_line()) return
               call read_table(i,i_lcapture)
               if (failed('read lcapture')) return
               if (.not. skip_line()) return
               call read_table(i,i_lneutrino)
               if (failed('read lneutrino')) return
            end do
         
            close(iounit)
         
            if (rates_use_cache) then
               open(newunit=iounit, file=trim(cache_filename), iostat=ios, &
                     action='write', form='unformatted')
               if (ios == 0) then
                  call write_weak_cache(iounit)
                  close(iounit)
               end if
            end if
         
         end if
         
         nullify(weak_reactions_dict)
         do i = 1, num_weak_reactions
            call create_weak_dict_key(weak_lhs_nuclide_name(i), weak_rhs_nuclide_name(i), key)
            call integer_dict_define(weak_reactions_dict, key, i, ierr)
            if (failed('integer_dict_define')) return
         end do
         
         call integer_dict_create_hash(weak_reactions_dict, ierr)
         if (failed('integer_dict_create_hash')) return

         do i = 1, num_weak_reactions
            associate(t => weak_reactions_tables(i) % t)
              if (ierr == 0) call t% setup(ierr)
            end associate
            if (failed('setup')) return
         end do

         if (dbg) write(*,*) 'finished load_weak_tables'
         
         
         contains


         subroutine read_weak_cache(iounit,ios)
            integer, intent(in) :: iounit
            integer, intent(out) :: ios
            integer :: n, i
            
            include 'formats'
         
            read(iounit,iostat=ios) num_weak_reactions
            if (ios /= 0) return
         
            if (dbg) write(*,2) 'num_weak_reactions', num_weak_reactions
         
            call alloc
            if (failed('allocate')) return
            
            n = num_weak_reactions
            
            read(iounit,iostat=ios) &
               weak_lhs_nuclide_id(1:n), &
               weak_rhs_nuclide_id(1:n), &
               weak_reaclib_id(1:n), &
               weak_lhs_nuclide_name(1:n), &
               weak_rhs_nuclide_name(1:n)

            do i = 1, n
               read(iounit, iostat=ios) &
                    weak_reactions_tables(i) % t % data(1,1:weak_num_T9,1:weak_num_lYeRho,1:3)
            end do
               
         end subroutine read_weak_cache


         subroutine write_weak_cache(iounit)
            integer, intent(in) :: iounit
            integer :: n, i
            
            include 'formats'
         
            write(iounit) num_weak_reactions
            
            n = num_weak_reactions
            
            write(iounit) &
               weak_lhs_nuclide_id(1:n), &
               weak_rhs_nuclide_id(1:n), &
               weak_reaclib_id(1:n), &
               weak_lhs_nuclide_name(1:n), &
               weak_rhs_nuclide_name(1:n)

            do i = 1, n
               write(iounit, iostat=ios) &
                    weak_reactions_tables(i) % t % data(1,1:weak_num_T9,1:weak_num_lYeRho,1:3)
            end do
               
         end subroutine write_weak_cache
                        
         
         subroutine alloc

            integer :: i
            type(weaklib_rate_table) :: table
         
            allocate( &
               weak_reaclib_id(num_weak_reactions), &
               weak_lhs_nuclide_name(num_weak_reactions), &
               weak_rhs_nuclide_name(num_weak_reactions), &
               weak_lhs_nuclide_id(num_weak_reactions), &
               weak_rhs_nuclide_id(num_weak_reactions), &
               weak_reactions_tables(num_weak_reactions), &
               stat=ierr)

            do i = 1, num_weak_reactions
               table = weaklib_rate_table(weak_reaction_T9s, weak_reaction_lYeRhos, .false.)
               allocate(weak_reactions_tables(i)% t, source=table)
            end do
               
         end subroutine alloc
         
         
         subroutine adjust_name(nm)
            character(len=iso_name_length) :: nm
            nm = adjustl(nm)
            if (nm == 'p') then
               nm = 'h1'
            else if (nm == 'n') then
               nm = 'neut'
            end if
         end subroutine adjust_name
         
         
         subroutine read_table(i,ii)
            use math_lib, only: str_to_vector
            integer, intent(in) :: i, ii
            integer :: k, j, skip, nvec
            !real :: buffer(weak_num_T9)
            character (len=256) :: buf
            real(dp), target :: vec_ary(50)
            real(dp), pointer :: vec(:)
            logical, parameter :: dbg = .false.
            vec => vec_ary
            skip = -1
            do j = 1, weak_num_lYeRho
               !read(iounit,fmt=*,iostat=ierr) skip, buffer
               read(iounit,fmt='(a)',iostat=ierr) buf
               if (ierr == 0) then
                  call str_to_vector(buf, vec, nvec, ierr)
                  skip = int(vec(1))
                  if (nvec < weak_num_T9+1) ierr = -1
               end if
               if (ierr /= 0 .or. j /= skip) then
                  if (dbg) then
                     write(*,*) 'error in reading table', j, skip
                     write(*,*) 'these are the NEXT lines after the error'
                     do k=1,20
                        read(iounit,fmt='(a)') string
                        write(*,'(a)') trim(string)
                     end do
                     write(*,*)
                     stop 'read_table'
                  end if
                  return
               end if
               do k=1,weak_num_T9
                  weak_reactions_tables(i) % t % data(1,k,j,ii) = vec(k+1)
               end do
               !if (dbg) write(*,'(a,2i6,99f9.3)') 'read_table', j, skip, buffer
            end do
         end subroutine read_table
         
         
         logical function failed(str)
            character (len=*) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*,*) 'failed: ' // trim(str)
            end if
         end function failed
         
         
         logical function skip_line()
            logical, parameter :: dbg = .false.
            if (dbg) then
               read(iounit,fmt='(a)') string
               write(*,'(a)') 'skip line ' // trim(string)
            else
               read(iounit,'(a)',iostat=ierr)
            end if
            skip_line = .not. (failed('skip line'))
         end function skip_line


      end subroutine private_load_weak_tables

      subroutine load_user_weak_tables(ierr)
        use utils_def
        use utils_lib
        use chem_lib, only: chem_get_iso_id
        use chem_def, only: iso_name_length
        use weak_support, only: parse_weak_rate_name

        integer, intent(out) :: ierr

        character (len=256) :: filename
        character(len=iso_name_length) :: lhs, rhs
        character(len=2*iso_name_length+1) :: key

        integer :: i, iounit, n, t, ir, id
        character (len=256) :: dir, rate_name, rate_fname, buffer

        logical, parameter :: dbg = .false.

        include 'formats'

        ierr = 0

        if (dbg) write(*,*) 'load_user_weak_tables'

        ! first try local rate_tables_dir
        dir = rates_table_dir
        filename = trim(dir) // '/weak_rate_list.txt'
        open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
        if (ierr /= 0) then ! if don't find that file, look in rates_dir
           dir = trim(rates_dir) // '/rate_tables'
           filename = trim(dir) // '/weak_rate_list.txt'
           ierr = 0
           open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
           if (ierr /= 0) then
              write(*,*) 'failed to open weak rates list file ' // trim(filename)
              return
           end if
        end if

        n = 0
        i = 0

        if (dbg) write(*,*) 'read rate list file ' // trim(filename)

        rate_loop: do
           t = token(iounit, n, i, buffer, rate_name)
           if (t == eof_token) exit
           if (t /= name_token) then
              call error; return
           end if
           if (dbg) write(*,*) 'use rate table from file for ', trim(rate_name)

           ! first, parse the weak rate id
           call parse_weak_rate_name(rate_name, lhs, rhs, ierr)
           if (dbg) write(*,*) 'parse_weak_rate_name gives ', trim(lhs), ' ', trim(rhs), ierr

           ! check if we already have a rate with this name
           call create_weak_dict_key(lhs, rhs, key)
           if (dbg) write(*,*) 'key is ', trim(key), ierr

           !write(*,'(a)') 'weak info list key ' // trim(key)
           call integer_dict_lookup(weak_reactions_dict, key, ir, ierr)
           if (dbg) write(*,*) ir, ierr


           if (ierr /= 0 .or. ir <= 0) then
              call extend
              ir = num_weak_reactions

              id = chem_get_iso_id(lhs)
              if (id <= 0) then
                 write(*,*) 'weaklib FATAL ERROR: unknown nuclide ' // lhs
                 call mesa_error(__FILE__,__LINE__)
              end if
              weak_lhs_nuclide_id(ir) = id

              id = chem_get_iso_id(rhs)
              if (id <= 0) then
                 write(*,*) 'weaklib FATAL ERROR: unknown nuclide ' // rhs
                 call mesa_error(__FILE__,__LINE__)
              end if
              weak_reaclib_id(ir) = 0
              weak_rhs_nuclide_id(ir) = id
              weak_lhs_nuclide_name(ir) = lhs
              weak_rhs_nuclide_name(ir) = rhs

           end if

           t = token(iounit, n, i, buffer, rate_fname)
           if (t /= string_token) then
              call error; return
           end if
           if (dbg) write(*,*) 'rate_fname ', trim(rate_fname)

           call read_hd5_file
           if (ierr < 0) then
              write(*,*) 'failed to read hdf5 file in load_user_weak_tables'
              call mesa_error(__FILE__,__LINE__)
           end if

           nullify(weak_reactions_dict)
           do i = 1, num_weak_reactions
              call create_weak_dict_key(weak_lhs_nuclide_name(i), weak_rhs_nuclide_name(i), key)
              call integer_dict_define(weak_reactions_dict, key, i, ierr)
              if (failed('integer_dict_define')) return
           end do

           call integer_dict_create_hash(weak_reactions_dict, ierr)
           if (failed('integer_dict_create_hash')) return

        end do rate_loop

        close(iounit)

        if (dbg) write(*,*) 'finished load_weak_tables'

      contains

        subroutine read_hd5_file

          use hdf5io_lib

          character (len=256)                 :: filename
          type(hdf5io_t)                      :: hi
          real(dp), allocatable, dimension(:) :: T9s, lYeRhos
          integer                             :: num_T9, num_lYeRho
          logical                             :: has_cc
          type(weaklib_rate_table)            :: table

          logical, parameter :: dbg = .false.

          filename = trim(dir) // '/' // trim(rate_fname)
          write(*,*) 'reading user weak rate file ', trim(filename)

          ! open file (read-only)

          hi = hdf5io_t(filename, OPEN_FILE)

          ! read axis data

          call hi% alloc_read_dset('T9s', T9s)
          num_T9 = SIZE(T9s)

          call hi% alloc_read_dset('lYeRhos', lYeRhos)
          num_lYeRho = SIZE(lYeRhos)

          ! check if table has coulomb corrections

          has_cc = hi% dset_exists('delta_Q') .AND. hi% dset_exists('Vs')
          
          ! create the table

          table = weaklib_rate_table(T9s, lYeRhos, has_cc)

          ! read data into it

          call hi% read_dset('ldecay', table% data(1, 1:num_T9, 1:num_lYeRho, table% i_ldecay))
          call hi% read_dset('lcapture', table% data(1, 1:num_T9, 1:num_lYeRho, table% i_lcapture))
          call hi% read_dset('lneutrino', table% data(1, 1:num_T9, 1:num_lYeRho, table% i_lneutrino))

          if (has_cc) then
             call hi% read_dset('delta_Q', table% data(1, 1:num_T9, 1:num_lYeRho, table% i_delta_Q))
             call hi% read_dset('Vs', table% data(1, 1:num_T9, 1:num_lYeRho, table% i_Vs))
          end if

          ! store the table

          allocate(weak_reactions_tables(ir)% t, source=table)
          associate(t => weak_reactions_tables(ir) % t)
            if (ierr == 0) call t% setup(ierr)
            if (failed('setup')) return
          end associate

          ! close file

          call hi% final()

        end subroutine read_hd5_file



        subroutine extend
            integer :: i, n
            type(weaklib_rate_table) :: table

            type(table_c), dimension(:), allocatable :: tmp_weak_reactions_tables

            integer, allocatable, dimension(:) :: &
                 tmp_weak_lhs_nuclide_id, tmp_weak_rhs_nuclide_id, tmp_weak_reaclib_id
            character(len=iso_name_length), dimension(:), allocatable :: &
                 tmp_weak_lhs_nuclide_name, tmp_weak_rhs_nuclide_name

            n = num_weak_reactions + 1

            allocate( &
                 tmp_weak_reaclib_id(n), &
                 tmp_weak_lhs_nuclide_name(n), &
                 tmp_weak_rhs_nuclide_name(n), &
                 tmp_weak_lhs_nuclide_id(n), &
                 tmp_weak_rhs_nuclide_id(n), &
                 stat=ierr)

            tmp_weak_reaclib_id(1:num_weak_reactions) = weak_reaclib_id
            tmp_weak_lhs_nuclide_name(1:num_weak_reactions) = weak_lhs_nuclide_name
            tmp_weak_rhs_nuclide_name(1:num_weak_reactions) = weak_rhs_nuclide_name
            tmp_weak_lhs_nuclide_id(1:num_weak_reactions) = weak_lhs_nuclide_id
            tmp_weak_rhs_nuclide_id(1:num_weak_reactions) = weak_rhs_nuclide_id

            allocate( &
                 weak_reaclib_id(n), &
                 weak_lhs_nuclide_name(n), &
                 weak_rhs_nuclide_name(n), &
                 weak_lhs_nuclide_id(n), &
                 weak_rhs_nuclide_id(n), &
                 stat=ierr)

            weak_reaclib_id(1:num_weak_reactions) = tmp_weak_reaclib_id(1:num_weak_reactions)
            weak_lhs_nuclide_name(1:num_weak_reactions) = tmp_weak_lhs_nuclide_name(1:num_weak_reactions)
            weak_rhs_nuclide_name(1:num_weak_reactions) = tmp_weak_rhs_nuclide_name(1:num_weak_reactions)
            weak_lhs_nuclide_id(1:num_weak_reactions) = tmp_weak_lhs_nuclide_id(1:num_weak_reactions)
            weak_rhs_nuclide_id(1:num_weak_reactions) = tmp_weak_rhs_nuclide_id(1:num_weak_reactions)

            deallocate( &
                 tmp_weak_reaclib_id, &
                 tmp_weak_lhs_nuclide_name, &
                 tmp_weak_rhs_nuclide_name, &
                 tmp_weak_lhs_nuclide_id, &
                 tmp_weak_rhs_nuclide_id, &
                 stat=ierr)

            allocate(tmp_weak_reactions_tables(n))
            tmp_weak_reactions_tables(1:num_weak_reactions) = weak_reactions_tables
            call move_alloc(tmp_weak_reactions_tables, weak_reactions_tables)

            num_weak_reactions = n

          end subroutine extend

          subroutine error
            ierr = -1
            close(iounit)
          end subroutine error

          logical function failed(str)
            character (len=*) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*,*) 'failed: ' // trim(str)
            end if
          end function failed


      end subroutine load_user_weak_tables

      end module load_weak

