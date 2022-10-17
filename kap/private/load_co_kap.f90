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

module load_CO_kap
   use kap_def
   use math_lib
   use const_def, only : dp, use_mesa_temp_cache
   use utils_lib, only : mesa_error, mv, switch_str
   
   implicit none
   
   integer, parameter :: min_version = 37
   
   logical, parameter :: CO_dbg = .false.

contains
   
   
   subroutine Setup_Kap_CO_Tables(rq, co_z_tables, use_cache, load_on_demand, ierr)
      type (Kap_General_Info), pointer :: rq
      type (Kap_CO_Z_Table), dimension(:), pointer :: co_z_tables
      logical, intent(in) :: use_cache, load_on_demand
      integer, intent(out) :: ierr
      integer :: iz, ix, CO_option
      include 'formats'
      ierr = 0
      CO_option = rq% kap_CO_option
      if (.not. associated(co_z_tables)) then
         allocate(co_z_tables(num_kap_CO_Zs(CO_option)), STAT = ierr)
         if (ierr /= 0) return
         do iz = 1, num_kap_CO_Zs(CO_option)
            co_z_tables(iz)% Zbase = kap_CO_Zs(iz, CO_option)
            co_z_tables(iz)% log10_Zbase = safe_log10(kap_CO_Zs(iz, CO_option))
            co_z_tables(iz)% Zfrac_C = -1
            co_z_tables(iz)% Zfrac_N = -1
            co_z_tables(iz)% Zfrac_O = -1
            co_z_tables(iz)% Zfrac_Ne = -1
            allocate(co_z_tables(iz)% x_tables(num_kap_CO_Xs(CO_option)), STAT = ierr)
            if (ierr /= 0) return
         end do
         if (show_allocations) write(*, 2) 'setup co_z_tables', &
            num_kap_CO_Zs(CO_option) + num_kap_CO_Zs(CO_option) * num_kap_CO_Xs(CO_option)
      end if
      do iz = 1, num_kap_CO_Zs(CO_option)
         do ix = 1, num_kap_CO_Xs(CO_option)
            call load_one_CO(rq, co_z_tables, iz, ix, load_on_demand, ierr)
            if (ierr /= 0) return
         end do
      end do
   end subroutine Setup_Kap_CO_Tables
   
   
   subroutine load_one_CO(rq, co_z_tables, iz, ix, read_later, ierr)
      use utils_lib
      type (Kap_General_Info), pointer :: rq
      type (Kap_CO_Z_Table), dimension(:), pointer :: co_z_tables
      integer, intent(in) :: iz, ix
      logical, intent(in) :: read_later
      integer, intent(out) :: ierr
      
      character (len = 256) :: fname, filename, cache_filename, temp_cache_filename
      
      ierr = 0
      
      call Get_CO_Filenames(rq, &
         kap_CO_Zs(iz, rq% kap_CO_option), kap_CO_Xs(ix, rq% kap_CO_option), kap_dir, fname, filename, &
         cache_filename, temp_cache_filename, ierr)
      if (ierr /= 0) return
      
      if (rq% show_info) write(*, *) 'read filename <' // trim(filename) // '>'
      
      call Prepare_Kap_CO_X_Table(rq, &
         iz, co_z_tables, co_z_tables(iz)% x_tables, &
         ix, kap_CO_Xs(ix, rq% kap_CO_option), kap_CO_Zs(iz, rq% kap_CO_option), &
         read_later, fname, filename, &
         cache_filename, temp_cache_filename, ierr)
      if (ierr /= 0) return
      
      if (.not. read_later) co_z_tables(iz)% x_tables(ix)% not_loaded_yet = .false.
   
   end subroutine load_one_CO
   
   
   subroutine Prepare_Kap_CO_X_Table(rq, &
      iz, co_z_tables, x_tables, ix, X_in, Z_in, read_later, &
      fname, filename, cache_filename, temp_cache_filename, ierr)
      ! allocates the arrays and stores data
      type (Kap_General_Info), pointer :: rq
      type (Kap_CO_Z_Table), dimension(:), pointer :: co_z_tables
      type (Kap_CO_X_Table), dimension(:), pointer :: x_tables
      integer, intent(in) :: iz, ix
      real(dp), intent(in) :: X_in, Z_in ! expected contents of the Kap file
      logical, intent(in) :: read_later
      character (len = *), intent(in) :: fname, filename, cache_filename, temp_cache_filename
      integer, intent(out) :: ierr
      
      integer :: io_unit, cache_io_unit
      
      real(dp) :: X, Z
      integer :: ios, form, version, n, num_tables, num_logRs, num_logTs, status
      real(dp) :: xin, zz, Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne, &
         logR_min, logR_max, logT_min, logT_max, xErr, zErr
      type (Kap_CO_Table), dimension(:), pointer :: co_tables
      integer :: num_dXC_gt_dXO ! the number of tables with dXC > dXO
      integer :: CO_table_numbers(num_kap_CO_dXs, num_kap_CO_dXs)
      integer :: next_dXO_table(max_num_CO_tables)
      integer :: next_dXC_table(max_num_CO_tables)
      character (len = 256) :: message
      real(dp), target :: vec_ary(30)
      real(dp), pointer :: vec(:)
      integer :: nvec
      real(dp), parameter :: tiny = 1d-6
      include 'formats'
      
      ierr = 0
      Zfrac_C = 0.d0
      Zfrac_N = 0.d0
      Zfrac_O = 0.d0
      Zfrac_Ne = 0.d0
      
      vec => vec_ary
      X = X_in
      Z = Z_in
      form = 0
      nvec = -1
      
      num_dXC_gt_dXO = -1
      CO_table_numbers = -1
      next_dXO_table = -1
      next_dXC_table = -1
      
      open(newunit = io_unit, file = trim(filename), action = 'read', status = 'old', iostat = ios)
      if (ios /= 0) then
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, *) 'NOTICE: missing kap data ' // trim(filename)
         write(*, *)
         write(*, *) 'Please check the validity of the kap_prefix string for this file.'
         write(*, *)
         write(*, *) 'If it is okay, you may need to install new kap data.'
         write(*, *) 'To do that, remove the directory mesa/data/kap_data,'
         write(*, *) 'and rerun the mesa ./install script.'
         write(*, '(A)')
         call mesa_error(__FILE__, __LINE__)
      end if
      
      version = -1
      read(io_unit, *, iostat = ierr) ! skip the 1st line
      if (ierr == 0) then
         read(io_unit, *, iostat = ierr) ! skip the 2nd line
         if (ierr == 0) then
            read(io_unit, '(a)', iostat = ierr) message
            if (ierr == 0) call str_to_vector(message, vec, nvec, ierr)
            if (nvec < 15) ierr = -1
            if (ierr == 0) then
               form = int(vec(1))
               version = int(vec(2))
               num_tables = int(vec(3))
               xin = vec(4)
               zz = vec(5)
               Zfrac_C = vec(6)
               Zfrac_N = vec(7)
               Zfrac_O = vec(8)
               Zfrac_Ne = vec(9)
               num_logRs = int(vec(10))
               logR_min = vec(11)
               logR_max = vec(12)
               num_logTs = int(vec(13))
               logT_min = vec(14)
               logT_max = vec(15)
            end if
         end if
      end if
      
      if (ierr /= 0 .or. version < min_version) then
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, *) 'FATAL ERROR: out-of-date version of kap data ' // trim(filename)
         write(*, *) 'Please update by removing the directory mesa/data/kap_data,'
         write(*, *) 'and rerunning the mesa ./install script.'
         write(*, '(A)')
         call mesa_error(__FILE__, __LINE__)
      end if
      
      if (form /= kap_table_co_enhanced_form) then
         call mesa_error(__FILE__, __LINE__, 'form /= kap_table_co_enhanced_form')
      end if
      
      if (max_num_CO_tables < num_tables) then
         ierr = -1
         return
      end if
      
      if (Zfrac_C < 0 .or. Zfrac_N < 0 .or. Zfrac_O < 0 .or. Zfrac_Ne < 0) then
         ierr = -1
         return
      end if
      
      co_z_tables(iz)% Zfrac_C = Zfrac_C
      co_z_tables(iz)% Zfrac_N = Zfrac_N
      co_z_tables(iz)% Zfrac_O = Zfrac_O
      co_z_tables(iz)% Zfrac_Ne = Zfrac_Ne
      
      call Setup_Kap_CO_X_Table(ierr)
      if (ierr /= 0) return
      
      if (CO_dbg) write(*, *) 'after Setup_Kap_CO_X_Table: associated co_tables', &
         associated(x_tables(ix)% co_tables), ix
      
      if (read_later) then
         close(io_unit)
         return
      end if
      
      if (kap_use_cache) then
         open(newunit = cache_io_unit, file = trim(cache_filename), action = 'read', &
            status = 'old', iostat = ios, form = 'unformatted')
         if (ios == 0) then ! try reading the cached data
            call Read_Kap_CO_X_Table(cache_io_unit, .true., ierr)
            close(cache_io_unit)
            if (ierr == 0) then
               close(io_unit)
               return
            end if
            ierr = 0
         else
            if (CO_dbg) write(*, '(a)') 'failed to open ' // trim(cache_filename)
         end if
      end if
      
      if (CO_dbg) write(*, *) 'before call: associated co_tables', &
         associated(x_tables(ix)% co_tables), ix
      if (show_allocations) write(*, 2) 'loading ' // trim(filename)
      call Read_Kap_CO_X_Table(io_unit, .false., ierr)
      close(io_unit)
      if (ierr /= 0) then
         write(*, *) 'failed in Read_Kap_CO_X_Table ' // trim(filename)
         return
      end if
      
      if (.not. kap_use_cache) return
      
      open(newunit = cache_io_unit, file = trim(switch_str(temp_cache_filename, cache_filename, use_mesa_temp_cache)), &
         iostat = ios, action = 'write', form = 'unformatted')
      if (ios == 0) then
         write(*, '(a)') 'write ' // trim(cache_filename)
         call Write_Kap_CO_X_Table_Cache(&
            x_tables, ix, cache_io_unit, x_tables(ix)% num_logRs, &
            x_tables(ix)% num_logTs, version)
         close(cache_io_unit)
         if(use_mesa_temp_cache) call mv(temp_cache_filename, cache_filename, .true.)
         if (kap_read_after_write_cache) then
            open(newunit = cache_io_unit, file = trim(cache_filename), action = 'read', &
               status = 'old', iostat = ios, form = 'unformatted')
            if (ios == 0) then
               call Read_Kap_CO_X_Table(cache_io_unit, .true., ierr)
               close(cache_io_unit)
            else
               ierr = -1
            end if
         end if
      end if
   
   
   contains
      
      
      subroutine Setup_Kap_CO_X_Table(ierr)
         integer, intent(out) :: ierr
         
         integer :: i
         include 'formats'
         
         xErr = abs(xin - X); zErr = abs(zz - Z)
         if (num_tables <= 0 .or. xErr > tiny .or. zErr > tiny) then
            ierr = -1
            write(*, *) 'bug in file ' // trim(filename), num_tables, xErr, zErr
            write(*, *) 'num_tables <= 0', num_tables <= 0
            write(*, *) 'xErr > tiny', xErr > tiny
            write(*, *) 'zErr > tiny', zErr > tiny
            return
         end if
         
         x_tables(ix)% not_loaded_yet = .true.
         
         x_tables(ix)% X = X
         x_tables(ix)% Z = Z
         
         x_tables(ix)% logR_min = logR_min
         x_tables(ix)% logR_max = logR_max
         x_tables(ix)% num_logRs = num_logRs
         nullify(x_tables(ix)% logRs)
         
         x_tables(ix)% logT_min = logT_min
         x_tables(ix)% logT_max = logT_max
         x_tables(ix)% num_logTs = num_logTs
         nullify(x_tables(ix)% logTs)
         
         x_tables(ix)% num_CO_tables = num_tables
         
         if (show_allocations) write(*, 2) 'co_tables', num_tables
         allocate(co_tables(num_tables), STAT = status)
         if (status .ne. 0) then
            ierr = -1
            if (CO_dbg) write(*, *) 'InsufficientMemory for Prepare_Kap_CO_X_Table', iz, ix
            return
         end if
         x_tables(ix)% co_tables => co_tables
         if (CO_dbg) write(*, *) 'Setup_Kap_CO_X_Table: allocate co_tables', iz, ix, num_tables, status
         if (CO_dbg) write(*, *) 'Setup_Kap_CO_X_Table: associated(x_tables(ix)% co_tables)', &
            associated(x_tables(ix)% co_tables), associated(co_tables), ix
         
         do i = 1, num_tables
            nullify(co_tables(i)% kap1) ! allocate when read the data
         end do
      
      end subroutine Setup_Kap_CO_X_Table
      
      
      subroutine Read_Kap_CO_X_Table(io_unit, reading_cache, ierr)
         integer, intent(in) :: io_unit ! use this for file access
         logical, intent(in) :: reading_cache
         integer, intent(out) :: ierr
         
         character (len = 256) :: message
         character (len = 1) :: char
         integer :: c_num_tables, c_num_logRs, c_num_logTs, c_version
         real(dp) :: c_xin, c_zz, c_logR_min, c_logR_max, c_logT_min, c_logT_max
         
         include 'formats'
         ierr = 0
         
         if (reading_cache) then
            
            ios = 0
            read(io_unit, iostat = ios) c_version
            if (ios /= 0 .or. c_version /= version) then
               ierr = 1
               if (ios /= 0) write(*, *) 'cache failed in read'
               if (c_version /= version) write(*, *) 'cache failed for c_version /= version'
               return
            end if
            
            read(io_unit, iostat = ios) &
               c_xin, c_zz, &
               c_logR_min, c_logR_max, c_num_logRs, x_tables(ix)% ili_logRs, &
               c_logT_min, c_logT_max, c_num_logTs, x_tables(ix)% ili_logTs
            if (ios /= 0 .or. c_num_logRs /= num_logRs .or. c_num_logTs /= num_logTs) then
               ierr = 1
               if (ios /= 0) write(*, *) 'cache failed in read'
               if (c_num_logRs /= num_logRs) write(*, *) 'cache failed for c_num_logRs /= num_logRs'
               if (c_num_logTs /= num_logTs) write(*, *) 'cache failed for c_num_logTs /= num_logTs'
               return
            end if
            
            read(io_unit, iostat = ios) &
               c_num_tables, num_dXC_gt_dXO, &
               CO_table_numbers, next_dXO_table, next_dXC_table
            
            if (ios /= 0) then
               ierr = 1
               if (ios /= 0) write(*, *) 'cache failed in read'
               return
            end if
         
         end if
         
         xErr = abs(xin - X); zErr = abs(zz - Z)
         if (num_tables <= 0 .or. xErr > tiny .or. zErr > tiny) then
            if (reading_cache) then
               if (num_tables <= 0) write(*, *) 'cache failed for num_tables <= 0'
               if (xErr > tiny) write(*, *) 'cache failed for xErr > tiny'
               if (zErr > tiny) write(*, *) 'cache failed for zErr > tiny'
               ierr = 1; return
            end if
            ierr = -1
            write(*, *) 'num_tables <= 0', num_tables <= 0
            write(*, *) 'xErr > tiny', xErr > tiny
            write(*, *) 'zErr > tiny', zErr > tiny
            return
         end if
         
         if (show_allocations) write(*, 2) 'Read_Kap_CO_X_Table logRs logTs', &
            num_logRs + num_logTs
         allocate(x_tables(ix)% logRs(num_logRs), x_tables(ix)% logTs(num_logTs), STAT = status)
         if (status .ne. 0) then
            ierr = -1
            write(*, *) 'InsufficientMemory for Read_Kap_CO_X_Table'
            return
         end if
         
         if (.not. reading_cache) then
            
            read(io_unit, *, iostat = ierr) ! skip line
            if (ierr /= 0) return
            read(io_unit, '(i3)', iostat = ierr) num_dXC_gt_dXO
            if (ierr /= 0) return
            x_tables(ix)% num_dXC_gt_dXO = num_dXC_gt_dXO
            do ! skip to the start of the first table
               read(io_unit, '(a1)', iostat = ierr) char
               if (ierr /= 0) return
               if (char == '-') exit
            end do
            read(io_unit, *, iostat = ierr) ! skip line
            if (ierr /= 0) return
            x_tables(ix)% CO_table_numbers = -1
            x_tables(ix)% next_dXC_table = -1
            x_tables(ix)% next_dXO_table = -1
         
         else
            
            read(io_unit, iostat = ierr) &
               x_tables(ix)% logRs(1:num_logRs), &
               x_tables(ix)% logTs(1:num_logTs)
            if (ierr /= 0) return
            
            x_tables(ix)% num_dXC_gt_dXO = num_dXC_gt_dXO
            x_tables(ix)% CO_table_numbers = CO_table_numbers
            x_tables(ix)% next_dXO_table = next_dXO_table
            x_tables(ix)% next_dXC_table = next_dXC_table
         
         end if
         
         do n = 1, num_tables
            if (CO_dbg) write(*, *) 'call Read_Kap_CO_Table', ix, n, X, Z
            call Read_Kap_CO_Table(x_tables, ix, n, X, Z, &
               num_logRs, num_logTs, io_unit, reading_cache, ierr)
            if (ierr /= 0) return
         end do
         
         if (.not. reading_cache) then
            do n = 1, num_tables
               next_dXC_table(n) = find_next_dXC_table(n)
               next_dXO_table(n) = find_next_dXO_table(n)
            end do
            x_tables(ix)% next_dXO_table = next_dXO_table
            x_tables(ix)% next_dXC_table = next_dXC_table
         end if
      
      end subroutine Read_Kap_CO_X_Table
      
      
      integer function find_next_dXC_table(i)
         integer, intent(in) :: i
         integer :: j
         real(dp) :: dXC, dXO
         dXC = co_tables(i)% dXC
         dXO = co_tables(i)% dXO
         if (dXC > dXO) then
            find_next_dXC_table = i + 1; return
         end if
         do j = 1, num_tables
            if (co_tables(j)% dXC > dXC .and. co_tables(j)% dXO == dXO) then
               find_next_dXC_table = j; return
            end if
         end do
         find_next_dXC_table = -1
      end function find_next_dXC_table
      
      
      integer function find_next_dXO_table(i)
         integer, intent(in) :: i
         integer :: j
         real(dp) :: dXC, dXO
         dXC = co_tables(i)% dXC
         dXO = co_tables(i)% dXO
         if (dXC < dXO) then
            find_next_dXO_table = i + 1; return
         end if
         do j = 1, num_tables
            if (co_tables(j)% dXO > dXO .and. co_tables(j)% dXC == dXC) then
               find_next_dXO_table = j
               return
            end if
         end do
         find_next_dXO_table = -1
      end function find_next_dXO_table
   
   
   end subroutine Prepare_Kap_CO_X_Table
   
   
   subroutine Read_Kap_CO_Table(&
      x_tables, ix, n, X_in, Z_in, num_logRs, num_logTs, io_unit, reading_cache, ierr)
      type (Kap_CO_X_Table), dimension(:), pointer :: x_tables
      integer, intent(in) :: ix ! index in x_tables
      integer, intent(in) :: n ! index in co_tables
      real(dp), intent(in) :: X_in, Z_in
      integer, intent(in) :: io_unit ! use this for file access
      logical, intent(in) :: reading_cache
      integer, intent(in) :: num_logRs, num_logTs
      integer, intent(out) :: ierr ! return nonzero if had trouble
      
      type (Kap_CO_Table), pointer :: co_tables(:)
      integer :: table_num, i, j, ios, status, iXC, iXO
      real(dp) :: X, Z, xin, zz, Y, dXC, dXO, err, logT
      real(dp), allocatable, target :: kap_table(:) ! data & spline coefficients
      real(dp), pointer :: kap(:, :, :)
      real(dp) :: logKs(num_logRs), logRs(num_logRs)
      character (len = 1000) :: message
      real(dp), target :: vec_ary(50)
      real(dp), pointer :: vec(:)
      integer :: nvec
      real(dp), parameter :: tiny = 1d-6
      
      logical :: store_logRs, store_logTs
      include 'formats'
      
      ierr = 0
      vec => vec_ary
      X = X_in; Z = Z_in
      nvec = -1
      
      if (show_allocations) write(*, 2) 'Read_Kap_CO_Table', &
         sz_per_Kap_point * num_logRs * num_logTs
      allocate(kap_table(sz_per_Kap_point * num_logRs * num_logTs))
      
      co_tables => x_tables(ix)% co_tables
      if (CO_dbg) write(*, *) 'Read_Kap_CO_Table associated(co_tables)', &
         associated(co_tables), ix, reading_cache
      
      if (.not. reading_cache) then
         read(io_unit, *, iostat = ierr)
         if (ierr /= 0) return
         read(io_unit, '(a)', iostat = ierr) message
         if (ierr == 0) call str_to_vector(message, vec, nvec, ierr)
         if (ierr /= 0 .or. nvec < 6) then
            write(*, *) 'Read_Kap_CO_Table str_to_vector'
            ierr = 1
            return
         end if
         table_num = int(vec(1))
         xin = vec(2)
         Y = vec(3)
         zz = vec(4)
         dXC = vec(5)
         dXO = vec(6)
         if (table_num /= n) then
            write(*, *) 'wrong num in opacity file for X', X, 'Z', Z, 'table', n
            ierr = -1
            return
         end if
      else ! reading_cache
         read(io_unit, iostat = ios) table_num, xin, Y, zz, dXC, dXO
         if (ios /= 0) then
            ierr = 1
            return
         end if
      end if
      
      if (abs(xin - X) > tiny .or. abs(zz - Z) > tiny) then
         if (reading_cache) then
            ierr = 1
            return
         end if
         write(*, *) 'error in opacity file for X', X, 'Z', Z, 'table', n
         ierr = -1
         return
      end if
      X = xin; Z = zz ! use the values from the file
      
      err = abs(1d0 - (X + Y + Z + dXC + dXO))
      if (err > tiny) then
         if (reading_cache) then
            ierr = 1; return
         end if
         write(*, *) 'abudance error in opacity table for X=', &
            X, 'Z=', Z, 'dXC=', dXC, 'dXO=', dXO
         ierr = -1
         return
      end if
      
      co_tables(n)% table_num = table_num
      co_tables(n)% X = X
      co_tables(n)% Z = Z
      co_tables(n)% dXC = dXC
      co_tables(n)% dXO = dXO
      co_tables(n)% dXC_lookup = get_dX_lookup(dXC, Z)
      co_tables(n)% dXO_lookup = get_dX_lookup(dXO, Z)
      
      if (show_allocations) write(*, 2) 'co_tables', &
         sz_per_Kap_point * num_logRs * num_logTs
      allocate(co_tables(n)% kap1(sz_per_Kap_point * num_logRs * num_logTs), STAT = status)
      if (status .ne. 0) then
         ierr = -1
         return
      end if
      
      if (reading_cache) then
         
         read(io_unit, iostat = ios) kap_table
         do i = 1, sz_per_Kap_point * num_logRs * num_logTs
            co_tables(n)% kap1(i) = kap_table(i)
         end do
         if (ios /= 0) then
            ierr = 4; return
         end if
      
      else
         
         iXC = find_in_dXs(dXC)
         if (iXC > 0) then
            iXO = find_in_dXs(dXO)
            if (iXO > 0) then
               x_tables(ix)% CO_table_numbers(iXC, iXO) = n
            end if
         end if
         
         read(io_unit, *, iostat = ierr)
         if (ierr /= 0) return
         read(io_unit, *, iostat = ierr)
         if (ierr /= 0) return
         
         read(io_unit, '(a)', iostat = ierr) message
         if (ierr == 0) call str_to_vector(message, vec, nvec, ierr)
         if (ierr /= 0) write(*, *) 'str_to_vector ierr', ierr
         if (nvec < num_logRs) ierr = -1
         if (ierr /= 0) then
            write(*, *) 'Read_Kap_CO_Table str_to_vector logRs: nvec, num_logRs', nvec, num_logRs
            write(*, '(a)') trim(message)
            stop
            return
         end if
         do i = 1, num_logRs
            logRs(i) = vec(i)
         end do
         
         if (n == 1) then
            x_tables(ix)% logRs(:) = logRs(:)
         else ! check
            do i = 1, num_logRs
               if (abs(logRs(i) - x_tables(ix)% logRs(i)) > 1d-7) then
                  write(*, *) 'problem with inconsistent logRs'
                  call mesa_error(__FILE__, __LINE__)
               end if
            end do
         end if
         
         read(io_unit, *, iostat = ierr)
         if (ierr /= 0) return
         
         kap(1:sz_per_Kap_point, 1:num_logRs, 1:num_logTs) => &
            co_tables(n)% kap1(1:sz_per_Kap_point * num_logRs * num_logTs)
         
         do i = 1, num_logTs
            
            read(io_unit, '(a)', iostat = ierr) message
            if (ierr == 0) call str_to_vector(message, vec, nvec, ierr)
            if (nvec < 1 + num_logRs) ierr = -1
            if (ierr /= 0) then
               write(*, *) 'Read_Kap_CO_Table str_to_vector logKs'
               return
            end if
            logT = vec(1)
            do j = 1, num_logRs
               logKs(j) = vec(j + 1)
               kap(1, j, i) = logKs(j)
            end do
            
            if (n == 1) then
               x_tables(ix)% logTs(i) = logT
            else ! check
               if (abs(logT - x_tables(ix)% logTs(i)) > 1d-7) then
                  write(*, *) 'problem with inconsistent logTs'
                  call mesa_error(__FILE__, __LINE__)
               end if
            end if
         
         end do
         read(io_unit, *, iostat = ierr)
         if (ierr /= 0) return
         
         call Make_CO_Interpolation_Data(&
            co_tables, n, x_tables(ix)% logRs, num_logRs, &
            x_tables(ix)% logTs, num_logTs, &
            x_tables(ix)% ili_logRs, x_tables(ix)% ili_logTs, ierr)
         
         if (ierr /= 0) then
            write(*, *) 'call Make_CO_Interpolation_Data ierr', ierr
            write(*, *) 'n', n
            write(*, *) 'num_logRs', num_logRs
            write(*, *) 'num_logTs', num_logTs
         end if
      
      end if
   
   end subroutine Read_Kap_CO_Table
   
   
   subroutine Make_CO_Interpolation_Data(&
      co_tables, n, logRs, num_logRs, logTs, num_logTs, ili_logRs, ili_logTs, ierr)
      use interp_2d_lib_db
      type (Kap_CO_Table), dimension(:), pointer :: co_tables
      integer, intent(in) :: n, num_logRs, num_logTs
      real(dp), intent(in), pointer :: logRs(:) ! =(num_logRs)
      real(dp), intent(in), pointer :: logTs(:) ! =(num_logTs)
      integer, intent(out) :: ili_logRs, ili_logTs, ierr
      
      real(dp), target :: table_ary(sz_per_kap_point * num_logRs * num_logTs)
      real(dp), pointer :: table(:, :, :), table1(:), kap(:, :, :)
      character (len = 256) :: message
      integer :: ibcxmin                   ! bc flag for x=xmin
      real(dp) :: bcxmin(num_logTs)               ! bc data vs. y at x=xmin
      integer :: ibcxmax                   ! bc flag for x=xmax
      real(dp) :: bcxmax(num_logTs)               ! bc data vs. y at x=xmax
      integer :: ibcymin                   ! bc flag for y=ymin
      real(dp) :: bcymin(num_logRs)               ! bc data vs. x at y=ymin
      integer :: ibcymax                   ! bc flag for y=ymax
      real(dp) :: bcymax(num_logRs)               ! bc data vs. x at y=ymax
      integer :: ier                       ! =0 on exit if there is no error.
      integer :: i, j
      real(dp) :: tol
      real(dp), pointer :: x_out(:), y_out(:), f_out(:, :)
      
      ! just use "not a knot" bc's at edges of tables
      ibcxmin = 0; bcxmin(1:num_logTs) = 0d0
      ibcxmax = 0; bcxmax(1:num_logTs) = 0d0
      ibcymin = 0; bcymin(1:num_logRs) = 0d0
      ibcymax = 0; bcymax(1:num_logRs) = 0d0
      
      table1 => table_ary
      table(1:sz_per_kap_point, 1:num_logRs, 1:num_logTs) => &
         table_ary(1:sz_per_kap_point * num_logRs * num_logTs)
      kap(1:sz_per_kap_point, 1:num_logRs, 1:num_logTs) => &
         co_tables(n)% kap1(1:sz_per_kap_point * num_logRs * num_logTs)
      
      do j = 1, num_logTs
         do i = 1, num_logRs
            table(1, i, j) = kap(1, i, j)
         end do
      end do
      
      call interp_mkbicub_db(&
         logRs, num_logRs, logTs, num_logTs, table1, num_logRs, &
         ibcxmin, bcxmin, ibcxmax, bcxmax, &
         ibcymin, bcymin, ibcymax, bcymax, &
         ili_logRs, ili_logTs, ier)
      
      if (ier /= 0) then
         write(*, *) 'interp_mkbicub_db error happened for Make_CO_Interpolation_Data for table', n
         ierr = -1
         return
      end if
      
      call Check_Interpolation_Data
      
      do i = 1, sz_per_kap_point * num_logRs * num_logTs
         co_tables(n)% kap1(i) = table1(i)
      end do
      
      ierr = 0
   
   
   contains
      
      subroutine Check_Interpolation_Data
         use utils_lib, only : is_bad
         integer :: i, iR, jtemp
         real(dp) :: val
         
         do i = 1, sz_per_kap_point
            do iR = 1, num_logRs
               do jtemp = 1, num_logTs
                  val = table(i, iR, jtemp)
                  if (is_bad(val)) then
                     if (.true.) then
                        write(*, *) 'bad value in xz', val, i, iR, jtemp
                        write(*, '(99(a15,3x,f15.8,3x))')  &
                           'logR', logRs(iR), 'logT', logTs(jtemp)
                     end if
                     table(i, iR, jtemp) = 0
                  end if
               end do
            end do
         end do
      
      end subroutine Check_Interpolation_Data
   
   
   end subroutine Make_CO_Interpolation_Data
   
   
   subroutine Write_Kap_CO_X_Table_Cache(x_tables, ix, io_unit, num_logRs, num_logTs, version)
      type (Kap_CO_X_Table), dimension(:), pointer :: x_tables
      integer, intent(in) :: ix, io_unit, num_logRs, num_logTs, version
      
      type (Kap_CO_Table), dimension(:), pointer :: co_tables
      
      integer :: num_tables, n, i
      real(dp) :: X, Z, Y, dXC, dXO
      real(dp), allocatable, target :: kap_table(:) ! data & spline coefficients
      include 'formats'
      
      if (show_allocations) write(*, 2) 'Write_Kap_CO_X_Table_Cache', &
         sz_per_Kap_point * num_logRs * num_logTs
      allocate(kap_table(sz_per_Kap_point * num_logRs * num_logTs))
      
      num_tables = x_tables(ix)% num_CO_tables
      X = x_tables(ix)% X
      Z = x_tables(ix)% Z
      co_tables => x_tables(ix)% co_tables
      
      write(io_unit) version
      
      write(io_unit) &
         X, Z, &
         x_tables(ix)% logR_min, &
         x_tables(ix)% logR_max, &
         num_logRs, x_tables(ix)% ili_logRs, &
         x_tables(ix)% logT_min, &
         x_tables(ix)% logT_max, &
         num_logTs, x_tables(ix)% ili_logTs
      
      write(io_unit) &
         num_tables, x_tables(ix)% num_dXC_gt_dXO, &
         x_tables(ix)% CO_table_numbers(:, :), &
         x_tables(ix)% next_dXO_table(:), &
         x_tables(ix)% next_dXC_table(:)
      
      write(io_unit) &
         x_tables(ix)% logRs(1:num_logRs), &
         x_tables(ix)% logTs(1:num_logTs)
      
      do n = 1, num_tables
         dXC = co_tables(n)% dXC
         dXO = co_tables(n)% dXO
         Y = 1 - (X + Z + dXC + dXO)
         do i = 1, sz_per_Kap_point * num_logRs * num_logTs
            kap_table(i) = co_tables(n)% kap1(i)
         end do
         write(io_unit) co_tables(n)% table_num, X, Y, Z, dXC, dXO
         write(io_unit) kap_table
      end do
   
   end subroutine Write_Kap_CO_X_Table_Cache
   
   
   real(dp) function get_dX_lookup(dX, Z)
      real(dp), intent(in) :: dX, Z
      get_dX_lookup = log10(Z + 1d-3 + dX)
   end function get_dX_lookup
   
   
   integer function find_in_dXs(dX)
      real(dp), intent(in) :: dX
      integer :: i
      real(dp), parameter :: tiny = 1d-6
      do i = 1, num_kap_CO_dXs
         if (abs(dX - kap_CO_dXs(i)) < tiny) then
            find_in_dXs = i; return
         end if
      end do
      find_in_dXs = -1
   end function find_in_dXs
   
   
   subroutine Get_CO_Filenames(rq, &
      Z, X, data_dir, fname, filename, cache_filename, temp_cache_filename, ierr)
      type (Kap_General_Info), pointer :: rq
      real(dp), intent(in) :: Z, X
      character (*), intent(in) :: data_dir
      character (*), intent(out) :: fname, filename, cache_filename, temp_cache_filename
      integer, intent(out) :: ierr
      ierr = 0
      call Create_fname(rq, Z, X, fname, ierr)
      filename = trim(data_dir) // '/' // trim(fname) // '.data'
      cache_filename = trim(kap_cache_dir) // '/' // trim(fname) // '.bin'
      temp_cache_filename = trim(kap_temp_cache_dir) // '/' // trim(fname) // '.bin'
   end subroutine Get_CO_Filenames
   
   
   subroutine Create_fname(rq, Z, X, fname, ierr)
      type (Kap_General_Info), pointer :: rq
      real(dp), intent(in) :: Z, X
      character (len = *), intent(out) :: fname
      integer, intent(out) :: ierr
      character (len = 16) :: zstr, xstr
      ierr = 0
      call get_output_string(Z, zstr, ierr)
      call get_output_string(X, xstr, ierr)
      
      write(fname, '(a)') trim(kap_CO_option_str(rq% kap_CO_option)) // '_z' // &
         trim(zstr) // '_x' // trim(xstr)
   
   end subroutine Create_fname

end module load_CO_kap

