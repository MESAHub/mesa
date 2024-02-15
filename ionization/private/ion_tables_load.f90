! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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
!
! ***********************************************************************

      module ion_tables_load
      use ionization_def
      use const_def, only: dp, use_mesa_temp_cache
      use utils_lib, only: mesa_error, mv, switch_str

      implicit none


      contains
            
      
      subroutine Init_ion_tables(file_prefix, Z1_suffix, use_cache, ierr)
         character(*), intent(in) :: file_prefix, Z1_suffix
         logical, intent(in) :: use_cache
         integer, intent(out) :: ierr ! 0 means AOK.
         ierr = 0
         if (ion_root_is_initialized) return
         ion_file_prefix = file_prefix
         ion_Z1_suffix = Z1_suffix
         use_cache_for_ion = use_cache
         call ion_read_sizes(ierr)
         if (ierr /= 0) return
         ion_root_is_initialized = .true.
      end subroutine Init_ion_tables
      
      
      subroutine ion_read_sizes(ierr)
         use utils_lib
         integer, intent(out) :: ierr

         character (len=256) :: message, fname, cache_filename, temp_cache_filename
         integer :: iounit
         real(dp) :: xin, zz
         
         ierr = 0
         
         call Get_ion_Table_Filenames(ion_Zs(1), ion_Xs(1), fname, cache_filename, temp_cache_filename)
         open(newunit=iounit, file=trim(fname), action='read', status='old', iostat=ierr)
         call check_for_error_in_ion_data(ierr, fname)

         read(iounit,*, iostat=ierr)
         call check_for_error_in_ion_data(ierr, fname)

         read(iounit,*, iostat=ierr) &
               ion_version, xin, zz, ion_num_logTs, ion_logT_min, ion_logT_max, ion_del_logT, &
               ion_num_logQs, ion_logQ_min, ion_logQ_max, ion_del_logQ
         call check_for_error_in_ion_data(ierr, fname)

         close(iounit)
         
         if (ion_version < min_version) call request_user_to_reinstall

      end subroutine ion_read_sizes
      
      
      subroutine request_user_to_reinstall
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*) 'NOTICE: you need to install a new version of the ion data.'
         write(*,*) 'Please update by removing the directory mesa/data/ionization_data,'
         write(*,*) 'and rerunning the mesa ./install script.'
         write(*,*)
         write(*,*)
         call mesa_error(__FILE__,__LINE__)
      end subroutine request_user_to_reinstall
      
      
      subroutine check_for_error_in_ion_data(ierr, fname)
         integer, intent(in) :: ierr
         character (len=*) :: fname
         if (ierr == 0) return
         write(*,*) 'load ion tables ' // trim(fname)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,'(a)') 'FATAL ERROR: missing or bad ion data.'
         write(*,'(a)') 'Please update by removing the directories mesa/data/ionization_data,'
         write(*,'(a)') 'and rerunning the mesa ./install script.'
         write(*,*)
         call mesa_error(__FILE__,__LINE__)
      end subroutine check_for_error_in_ion_data
      
      
      subroutine Load_ion_Table(ierr)
         use utils_lib
         integer, intent(out) :: ierr
         
         integer :: iz, ix, i
         
         ierr = 0
!$OMP CRITICAL (load_ionization_table)
         call do_read
!$OMP END CRITICAL (load_ionization_table)
         
         contains
         
         subroutine do_read
            integer :: sz_ion_tbl
            
            if (ion_is_initialized) return
            sz_ion_tbl = sz_per_ion_point*num_ion_vals* &
                                 ion_num_logQs*ion_num_logTs*num_ion_Xs*num_ion_Zs
            allocate(ion_tbl1(sz_ion_tbl),&
                     ion_logQs(ion_num_logQs), ion_logTs(ion_num_logTs),  &
                     STAT=ierr)
            if (ierr /= 0) return
            ion_tbl(1:sz_per_ion_point, 1:num_ion_vals,&
                    1:ion_num_logQs, 1:ion_num_logTs, 1:num_ion_Xs, 1:num_ion_Zs) => &
               ion_tbl1(1:sz_ion_tbl)

            ion_logQs(1) = ion_logQ_min
            do i = 2, ion_num_logQs-1
               ion_logQs(i) = ion_logQs(i-1) + ion_del_logQ
            end do
            ion_logQs(ion_num_logQs) = ion_logQ_max
            
            ion_logTs(1) = ion_logT_min
            do i = 2, ion_num_logTs-1
               ion_logTs(i) = ion_logTs(i-1) + ion_del_logT
            end do
            ion_logTs(ion_num_logTs) = ion_logT_max

            do iz = 1, num_ion_Zs
               do ix = 1, num_ion_Xs
                  if (ion_Zs(iz) + ion_Xs(ix) > 1.0000001d0) cycle
                  call read_one(ix,iz,ierr)
                  if (ierr /= 0) return
               end do
            end do
         
            ion_is_initialized = .true.
         end subroutine do_read
         
         subroutine read_one(ix,iz,ierr)
            integer, intent(in) :: ix, iz
            integer, intent(out) :: ierr

            character (len=256) :: fname, cache_filename, temp_cache_filename
            
            include 'formats'
            
            call Get_ion_Table_Filenames(&
                     ion_Zs(iz), ion_Xs(ix), fname, cache_filename, temp_cache_filename)
            call Load1_ion_Table(&
                  ion_Xs(ix), ion_Zs(iz), ion_tbl(:,:,:,:,ix,iz),&
                  fname, cache_filename, temp_cache_filename, use_cache_for_ion, ierr)
            if (ierr /= 0) then
               write(*,*) 'Load1_ion_Table ierr', ierr, ix, iz, ion_Xs(ix), ion_Zs(iz)
            end if

         end subroutine read_one
         
      end subroutine Load_ion_Table
      
      
      subroutine Get_ion_Table_Filenames(Z, X, fname, cache_filename, temp_cache_filename)
         use const_def, only: mesa_data_dir
         real(dp), intent(in) :: Z, X
         character (len=*), intent(out) :: fname, cache_filename, temp_cache_filename
         character (len=256) :: Zstr, Xstr, suffix
         
         call setstr(Z,Zstr)
         call setstr(X,Xstr)
         if (Zstr == '100') then
            suffix = ion_Z1_suffix
         else
            suffix = ''
         end if
         
         fname = trim(mesa_data_dir) // &
               '/ionization_data/' // trim(ion_file_prefix) // '_' //&
               trim(Zstr) // 'z' // trim(Xstr) // 'x' // trim(suffix) // '.data'
         cache_filename = trim(ionization_cache_dir) // &
               '/' // trim(ion_file_prefix) // '_' //&
               trim(Zstr) // 'z' // trim(Xstr) // 'x' // trim(suffix) // '.bin'
         
         temp_cache_filename = trim(ionization_temp_cache_dir) // &
               '/' // trim(ion_file_prefix) // '_' //&
               trim(Zstr) // 'z' // trim(Xstr) // 'x' // trim(suffix) // '.bin'
         
         contains
         
         subroutine setstr(v,str)
            real(dp), intent(in) :: v
            character (len=*) :: str
            if (v > 0.99999d0) then
               str = '100'
            else if (v > 0.09999d0) then
               write(str, '(i2)') floor(100d0 * v + 0.5d0)
            else
               write(str, '(a,i1)') '0', floor(100d0 * v + 0.5d0)
            end if
         end subroutine setstr
              
      end subroutine Get_ion_Table_Filenames
      
      
      subroutine Load1_ion_Table(&
            X, Z, tbl, filename, cache_filename, temp_cache_filename, use_cache, info)
         real(dp), intent(in) :: X, Z
         real(dp) :: tbl(sz_per_ion_point, num_ion_vals, ion_num_logQs, ion_num_logTs)
         character (*), intent(in) :: filename, cache_filename, temp_cache_filename
         logical, intent(in) :: use_cache
         integer, intent(out) :: info

         integer :: io_unit, cache_io_unit
         
         integer :: num_logQs_in, num_logTs_in, version_in
         real(dp) :: logT_min_in, logT_max_in, del_logT_in, vals(num_ion_vals),&
               logQ_min_in, logQ_max_in, del_logQ_in, X_in, Z_in, logQ, logT
         integer :: j,i,k,iQ,ios,status,line_number
         character (len=500) :: message, input_line
         real(dp), parameter :: tiny = 1e-6
         
         include 'formats'

         info = 0            

         write(message,*) 'open ', trim(filename)
         open(NEWUNIT=io_unit, FILE=trim(filename), ACTION='READ', STATUS='OLD', IOSTAT=ios)
         call check_for_error_in_ion_data(ios, filename)

         line_number = 0
         read(io_unit,*,iostat=info)
         if (info /= 0) return
         line_number = line_number + 1
         read(io_unit,*,iostat=info) &
               version_in, X_in, Z_in, num_logTs_in, logT_min_in, logT_max_in, del_logT_in, &
               num_logQs_in, logQ_min_in, logQ_max_in, del_logQ_in
         if (info /= 0) return
         line_number = line_number + 1
         read(io_unit,*,iostat=info)
         if (info /= 0) return
         line_number = line_number + 1

         if (ion_version > version_in &
            .or. ion_num_logQs /= num_logQs_in &
            .or. ion_num_logTs /= num_logTs_in&
            .or. abs(X-X_in) > tiny&
            .or. abs(Z-Z_in) > tiny&
            .or. abs(ion_logT_min-logT_min_in) > tiny    &
            .or. abs(ion_logT_max-logT_max_in) > tiny    &
            .or. abs(ion_del_logT-del_logT_in) > tiny    &
            .or. abs(ion_logQ_min-logQ_min_in) > tiny    &
            .or. abs(ion_logQ_max-logQ_max_in) > tiny    &
            .or. abs(ion_del_logQ-del_logQ_in) > tiny    &
           ) then
            write(*,*) 'bad header info in ' // trim(filename)
            info = -1
            close(io_unit)
            write(*,'(a50,l1,2i10)') 'ion_version > version_in', ion_version > version_in, ion_version, version_in
            write(*,'(a50,l1)') 'ion_num_logQs /= num_logQs_in', ion_num_logQs /= num_logQs_in
            write(*,'(a50,l1)') 'ion_num_logTs /= num_logTs_in', ion_num_logTs /= num_logTs_in
            write(*,'(a50,l1)') 'abs(X-X_in) > tiny', abs(X-X_in) > tiny
            write(*,'(a50,l1)') 'abs(Z-Z_in) > tiny', abs(Z-Z_in) > tiny
            write(*,'(a50,l1)') 'abs(ion_logT_min-logT_min_in) > tiny', abs(ion_logT_min-logT_min_in) > tiny
            write(*,'(a50,l1)') 'abs(ion_logT_max-logT_max_in) > tiny', abs(ion_logT_max-logT_max_in) > tiny
            write(*,'(a50,l1)') 'abs(ion_del_logT-del_logT_in) > tiny', abs(ion_del_logT-del_logT_in) > tiny
            write(*,'(a50,l1)') 'abs(ion_logQ_min-logQ_min_in) > tiny', abs(ion_logQ_min-logQ_min_in) > tiny
            write(*,'(a50,l1)') 'abs(ion_logQ_max-logQ_max_in) > tiny', abs(ion_logQ_max-logQ_max_in) > tiny
            write(*,'(a50,l1)') 'abs(ion_del_logQ-del_logQ_in) > tiny', abs(ion_del_logQ-del_logQ_in) > tiny
            write(*,*)
            write(*,1) 'ion_logT_max', ion_logT_max
            write(*,1) 'logT_max_in', logT_max_in
            stop 
            return
         end if

         if (use_cache) then
            call Read_ion_Cache(X, Z, tbl, cache_filename, ios)
            if (ios == 0) then
               close(io_unit)
               return
            end if
         end if
         
         do iQ=1,ion_num_logQs
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            line_number = line_number + 1
            read (io_unit,*,iostat=info) logQ
            if (failed('skip line')) return
            line_number = line_number + 1
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            line_number = line_number + 1
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            line_number = line_number + 1
            do i=1,ion_num_logTs
               read(io_unit,'(a)',iostat=info) input_line
               if (failed('read line')) return
               line_number = line_number + 1
               read (input_line,*,iostat=info) logT, vals(1:num_ion_vals)
               if (failed('read tbl')) then
                  write(*,'(a)') trim(input_line)
                  write(*,*) trim(filename)
                  write(*,*) 'iQ, i', iQ, i
                  write(*,*) 'num_ion_vals', num_ion_vals
                  write(*,*) 'logQ', logQ
                  write(*,*) 'bad input line?'
                  call mesa_error(__FILE__,__LINE__)
               end if
               tbl(1,1:num_ion_vals,iQ,i) = vals(1:num_ion_vals)
            enddo
            if(iQ < ion_num_logQs) read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            line_number = line_number + 1
            if(iQ < ion_num_logQs) read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            line_number = line_number + 1
         end do
            
         close(io_unit)
         
         call Make_ion_Interpolation_Data(tbl, info)
         if (failed('Make_ion_Interpolation_Data')) return
         
         call Check_ion_Interpolation_Data(tbl)
         
         if (.not. use_cache) return

         open(NEWUNIT=cache_io_unit, file=trim(switch_str(temp_cache_filename, cache_filename, use_mesa_temp_cache)),&
          iostat=ios, action='write', form='unformatted')

         if (ios == 0) then
            write(*,'(a)') 'write ' // trim(cache_filename)
            write(cache_io_unit) &
               X_in, Z_in, ion_num_logTs, ion_logT_min, ion_logT_max, ion_del_logT, &
               ion_num_logQs, ion_logQ_min, ion_logQ_max, ion_del_logQ, ion_version
            write(cache_io_unit) &
                  tbl(1:sz_per_ion_point, 1:num_ion_vals, 1:ion_num_logQs, 1:ion_num_logTs)
            close(cache_io_unit)
            if(use_mesa_temp_cache) call mv(temp_cache_filename, cache_filename,.true.)
         end if
         
         contains
         
         subroutine Check_ion_Interpolation_Data(tbl)
            use utils_lib,only:is_bad
            real(dp) :: tbl(sz_per_ion_point, num_ion_vals, ion_num_logQs, ion_num_logTs)
         
            ! for logT > 6.8 and logRho < -10, splines can get bogus higher order terms
            ! replace NaN's and Infinities with 0
         
            integer :: i, j, iQ, jtemp
            
            do i = 1, sz_per_ion_point
               do j = 1, num_ion_vals
                  do iQ = 1, ion_num_logQs
                     do jtemp = 1, ion_num_logTs
                        if (is_bad(tbl(i,j,iQ,jtemp))) then
                           tbl(i,j,iQ,jtemp) = 0
                        end if
                     end do
                  end do
               end do
            end do
         
         end subroutine Check_ion_Interpolation_Data
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (info /= 0)
            if (failed) then
               write(*,*)
               write(*,'(a)') trim(filename)
               write(*,'(a,i9)') &
                  ' Load1_ion_Table failed: ' // trim(str) // ' line', line_number
            end if
         end function failed
         

      end subroutine Load1_ion_Table
      
      
      subroutine Make_ion_Interpolation_Data(tbl, info)
         use interp_2d_lib_db
         use const_def, only: crad, ln10

         real(dp) :: tbl(sz_per_ion_point, num_ion_vals, ion_num_logQs, ion_num_logTs)
         integer, intent(out) :: info

         real(dp) :: logQs(ion_num_logQs)              ! x vector, strict ascending
         real(dp) :: logTs(ion_num_logTs)                    ! y vector, strict ascending
         real(dp) :: Ts(ion_num_logTs)
         real(dp), allocatable, target :: f1_ary(:) ! data & spline coefficients
         real(dp), pointer :: f1(:), f(:,:,:)
         integer :: ibcxmin                   ! bc flag for x=xmin
         real(dp) :: bcxmin(ion_num_logTs)    ! bc data vs. y at x=xmin
         integer :: ibcxmax                   ! bc flag for x=xmax
         real(dp) :: bcxmax(ion_num_logTs)     ! bc data vs. y at x=xmax
         integer :: ibcymin                   ! bc flag for y=ymin
         real(dp) :: bcymin(ion_num_logQs)   ! bc data vs. x at y=ymin
         integer :: ibcymax                   ! bc flag for y=ymax
         real(dp) :: bcymax(ion_num_logQs)   ! bc data vs. x at y=ymax
         integer :: ili_logQs    ! =1: logRho grid is "nearly" equally spaced
         integer :: ili_logTs      ! =1: logT grid is "nearly" equally spaced
         integer :: ier            ! =0 on exit if there is no error.
         real(dp) :: logQ, Rho, logRho, T, P, Cv, chiRho, chiT, logT, logT0, logT1, logQ0, logQ1
         real(dp) :: gamma3, gamma1, grad_ad, Prad, E, S
         integer :: iQ, jtemp, ilogT, ilogQ
         real(dp) :: fval(num_ion_vals), df_dx(num_ion_vals), df_dy(num_ion_vals)
         
         integer :: v, vlist(3), var, i, j
         character (len=256) :: message
         
         allocate(f1_ary(sz_per_ion_point*ion_num_logQs*ion_num_logTs))
         
         f1 => f1_ary
         f(1:sz_per_ion_point,1:ion_num_logQs,1:ion_num_logTs) => &
            f1_ary(1:sz_per_ion_point*ion_num_logQs*ion_num_logTs) 

         info = 0

         do iQ = 1, ion_num_logQs
            logQs(iQ) = ion_logQ_min + (iQ-1) * ion_del_logQ
         end do

         do jtemp = 1, ion_num_logTs
            logTs(jtemp) = ion_logT_min + (jtemp-1) * ion_del_logT
         end do

         ! just use "not a knot" bc's at edges of tables
         ibcxmin = 0; bcxmin(:) = 0
         ibcxmax = 0; bcxmax(:) = 0
         ibcymin = 0; bcymin(:) = 0
         ibcymax = 0; bcymax(:) = 0

         ! create tables for bicubic spline interpolation         
         do v = 1, num_ion_vals
         
            f(1,:,:) = tbl(1,v,:,:)
            call interp_mkbicub_db(&
                  logQs,ion_num_logQs,logTs,ion_num_logTs,f1,ion_num_logQs,&
                  ibcxmin,bcxmin,ibcxmax,bcxmax,&
                  ibcymin,bcymin,ibcymax,bcymax,&
                  ili_logQs,ili_logTs,ier)
            if (ier /= 0) then
               write(*,*) 'Make_ion_Interpolation_Data error happened for ion_value', v
               info = 3
               return
            end if
            tbl(2:4,v,:,:) = f(2:4,:,:)
            
         end do
         
      end subroutine Make_ion_Interpolation_Data
      
      
      subroutine Read_ion_Cache(X, Z, tbl, cache_filename, ios)
         real(dp), intent(in) :: X, Z
         real(dp) :: tbl(sz_per_ion_point, num_ion_vals, ion_num_logQs, ion_num_logTs)
         character (*), intent(in) :: cache_filename
         integer, intent(out) :: ios

         integer :: io_unit ! use this for file access

         real(dp) :: X_in, Z_in, logT_min_in, logT_max_in, del_logT_in, &
               logQ_min_in, logQ_max_in, del_logQ_in
         integer :: num_logQs_in, num_logTs_in, version_in
         real(dp), parameter :: tiny = 1d-6
         
         ios = 0
         open(newunit=io_unit,file=trim(cache_filename),action='read',&
               status='old',iostat=ios,form='unformatted')
         if (ios /= 0) return
         
         read(io_unit, iostat=ios) &
               X_in, Z_in, num_logTs_in, logT_min_in, logT_max_in, del_logT_in, &
               num_logQs_in, logQ_min_in, logQ_max_in, del_logQ_in, version_in
         if (ios /= 0) return
         
         if (ion_version /= version_in) then
            ios = 1
            write(*,*) 'read cache failed for version_in'
         end if
         if (ion_num_logQs /= num_logQs_in) then
            ios = 1
            write(*,*) 'read cache failed for ion_num_logQs'
         end if 
         if (ion_num_logTs /= num_logTs_in) then
            ios = 1
            write(*,*) 'read cache failed for ion_num_logTs'
         end if
         if (abs(X-X_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for X_in'
         end if
         if (abs(Z-Z_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for Z_in'
         end if
         if (abs(ion_logT_min-logT_min_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ion_logT_min'
         end if    
         if (abs(ion_logT_max-logT_max_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ion_logT_max'
         end if    
         if (abs(ion_del_logT-del_logT_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ion_del_logT'
         end if    
         if (abs(ion_logQ_min-logQ_min_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ion_logQ_min'
         end if    
         if (abs(ion_logQ_max-logQ_max_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ion_logQ_max'
         end if
         if (abs(ion_del_logQ-del_logQ_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ion_del_logQ'
         end if
         
         if (ios /= 0) then
            close(io_unit); return
         end if

         read(io_unit, iostat=ios) tbl(1:sz_per_ion_point, &
                        1:num_ion_vals, 1:ion_num_logQs, 1:ion_num_logTs)
         if (ios /= 0) then
            close(io_unit); return
         end if
         
         close(io_unit)

      end subroutine Read_ion_Cache
      
      
      end module ion_tables_load
