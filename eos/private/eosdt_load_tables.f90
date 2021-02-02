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
!
! ***********************************************************************

      module eosDT_load_tables
      use eos_def
      use utils_lib, only: is_bad, mesa_error, mv, switch_str
      use const_def, only: dp, use_mesa_temp_cache
      use math_lib

      implicit none
      
      ! the file EOS data
      integer, parameter :: jlogPgas = 1
      integer, parameter :: jlogE = 2
      integer, parameter :: jlogS = 3
      integer, parameter :: jchiRho = 4
      integer, parameter :: jchiT = 5
      integer, parameter :: jCp = 6
      integer, parameter :: jCv = 7
      integer, parameter :: jdE_dRho = 8
      integer, parameter :: jdS_dT = 9
      integer, parameter :: jdS_dRho = 10
      integer, parameter :: jmu = 11
      integer, parameter :: jlogfree_e = 12
      integer, parameter :: jgamma1 = 13
      integer, parameter :: jgamma3 = 14
      integer, parameter :: jgrad_ad = 15
      integer, parameter :: jeta = 16
      integer, parameter :: num_eos_file_vals = 16

      integer, parameter :: file_max_num_logQs = 1000

      

      contains
      
      
      subroutine request_user_to_reinstall
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*) 'NOTICE: you need to install a new verion of the eos data.'
         write(*,*) 'Please update by removing the directory mesa/data/eosDT_data,'
         write(*,*) 'and rerunning the mesa ./install script.'
         write(*,*)
         write(*,*)
         call mesa_error(__FILE__,__LINE__)
      end subroutine request_user_to_reinstall
      
      
      subroutine check_for_error_in_eosDT_data(ierr, fname)
         integer, intent(in) :: ierr
         character (len=*) :: fname
         if (ierr == 0) return
         write(*,*) 'load eos tables ' // trim(fname)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,'(a)') 'FATAL ERROR: missing or bad eos data.'
         write(*,'(a)') 'Please update by removing the directories ' &
            // 'mesa/data/eos*_data' &
            // ' and rerunning the mesa ./install script.'
         write(*,*)
         call mesa_error(__FILE__,__LINE__)
      end subroutine check_for_error_in_eosDT_data

    
    subroutine load_single_eosDT_table_by_id( &
             rq, which_eosdt, ep, ix, iz, ierr)
         use utils_lib
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: which_eosdt
         type (EosDT_XZ_Info), pointer :: ep
         integer,intent(in) :: iz, ix
         integer, intent(out) :: ierr
      
      if (which_eosdt == eosdt_max_FreeEOS) then
         ep => FreeEOS_XZ_data(ix,iz)
         if (FreeEOS_XZ_loaded(ix,iz)) return
      else if (which_eosdt == eosdt_OPAL_SCVH) then
         ep => eosDT_XZ_data(ix,iz)
         if (eosDT_XZ_loaded(ix,iz)) return
      else
         ierr = -1
         return
      end if
      
!$OMP CRITICAL(eosDT_load)
      if (which_eosdt == eosdt_max_FreeEOS) then
         if (.not. FreeEOS_XZ_loaded(ix,iz)) call do_read
      else
         if (.not. eosDT_XZ_loaded(ix,iz)) call do_read
      end if
!$OMP END CRITICAL(eosDT_load)
         
         contains
         
         subroutine do_read
            call read_one(ix,iz,ierr)
            if (ierr /= 0) return
            if (which_eosdt == eosdt_max_FreeEOS) then
               FreeEOS_XZ_loaded(ix,iz) = .true.
            else
               eosDT_XZ_loaded(ix,iz) = .true.
            end if
         end subroutine do_read
         
         subroutine read_one(ix,iz,ierr)
            use const_def, only: mesa_data_dir
            integer, intent(in) :: ix, iz
            integer, intent(out) :: ierr
            character (len=256) :: fname, cache_filename, temp_cache_filename
            integer :: iounit1, iounit2        
            real(dp) :: X, Z    
            type (DT_xz_Info), pointer :: xz
            include 'formats'
            iounit1 = alloc_iounit(ierr); if (ierr /= 0) return
            iounit2 = alloc_iounit(ierr); if (ierr /= 0) return            
            if (which_eosdt == eosdt_max_FreeEOS) then
               xz => FreeEOS_xz_struct
            else
               xz => eosDT_xz_struct
            end if
            call Get_eosDT_Table_Filenames(rq, which_eosdt, xz, &
               ix, iz, mesa_data_dir, fname, cache_filename, temp_cache_filename)
            call Load1_eosDT_Table(rq, which_eosdt, ep, xz, &
               ix, iz, fname, cache_filename, temp_cache_filename, iounit1, iounit2, use_cache_for_eos, ierr)
            if (ierr /= 0) then
               write(*,*) 'Load1_eosDT_Table ierr', ierr, ix, iz, X, Z
            end if
            call free_iounit(iounit2)
            call free_iounit(iounit1)
         end subroutine read_one
         
      end subroutine load_single_eosDT_table_by_id
      
      
      subroutine Get_eosDT_Table_Filenames(rq, which_eosdt, xz, &
            ix, iz, data_dir, fname, cache_filename, temp_cache_filename)
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: which_eosdt
         type (DT_xz_Info), pointer :: xz
         integer, intent(in) :: ix, iz
         character (*), intent(in) :: data_dir
         character (len=*), intent(out) :: fname, cache_filename, temp_cache_filename
         character (len=256) :: Zstr, Xstr, suffix, data_dir_name, data_prefix
         real(dp) :: X, Z
         
         Z = xz% Zs(iz)
         X = xz% Xs_for_Z(ix,iz)
         
         call setstr(Z,Zstr)
         call setstr(X,Xstr)
         suffix = ''         
         if (which_eosdt == eosdt_max_FreeEOS) then
            data_dir_name = '/eosFreeEOS_data/'
            data_prefix = '-FreeEOS_'
            suffix = rq% suffix_for_FreeEOS_Z(iz)
         else
            data_dir_name = '/eosDT_data/'
            data_prefix = '-eosDT_'
         end if
         
         fname = trim(data_dir) //  &
               trim(data_dir_name) // trim(rq% eosDT_file_prefix) // trim(data_prefix) // &
               trim(Zstr) // 'z' // trim(Xstr) // 'x' // trim(suffix) // '.data'
         cache_filename = trim(eosDT_cache_dir) //  &
               '/' // trim(rq% eosDT_file_prefix) // trim(data_prefix) // &
               trim(Zstr) // 'z' // trim(Xstr) // 'x' // trim(suffix) // '.bin'         
         temp_cache_filename = trim(eosDT_temp_cache_dir) //  &
               '/' // trim(rq% eosDT_file_prefix) // trim(data_prefix) // &
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
              
      end subroutine Get_eosDT_Table_Filenames
      
      
      subroutine Load1_eosDT_Table(rq, which_eosdt, ep, xz, &
            ix, iz, filename, cache_filename, temp_cache_filename, &
            io_unit, cache_io_unit, use_cache, info)
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: which_eosdt
         type (EosDT_XZ_Info), pointer :: ep
         type (DT_xz_Info), pointer :: xz
         integer, intent(in) :: ix, iz
         character (*), intent(in) :: filename, cache_filename, temp_cache_filename
         integer, intent(in) :: io_unit, cache_io_unit
         logical, intent(in) :: use_cache
         integer, intent(out) :: info
         
         real(dp) :: X, Z, logQ, logT, X_in, Z_in
         integer :: j, i, k, iQ, ios, status
         character (len=1000) :: message
         real(dp), parameter :: tiny = 1d-10
         real(dp), pointer :: tbl(:,:,:,:) ! => ep% tbl1
         real(dp), pointer :: tbl2_1(:), tbl2(:,:,:)
         real(dp), target :: vec_ary(50)
         real(dp), pointer :: vec(:)
         integer :: n
         
         include 'formats'

         info = 0    
         vec => vec_ary
         Z = xz% Zs(iz)
         X = xz% Xs_for_Z(ix,iz)

         write(message,*) 'open ', trim(filename)
         
         open(UNIT=io_unit, FILE=trim(filename), ACTION='READ', STATUS='OLD', IOSTAT=ios)
         call check_for_error_in_eosDT_data(ios, filename)

         read(io_unit,*,iostat=info)
         if (info /= 0) return
         
         read(io_unit,'(a)',iostat=info) message
         if (info == 0) call str_to_vector(message, vec, n, info)
         if (info /= 0 .or. n < 11) then
            write(*,'(a)') 'failed while reading ' // trim(filename)
            close(io_unit)
            info = -1
            return
         end if
         ep% version = int(vec(1))
         X_in = vec(2)
         Z_in = vec(3)
         ep% num_logTs = int(vec(4))
         ep% logT_min = vec(5)
         ep% logT_max = vec(6)
         ep% del_logT = vec(7)
         ep% num_logQs = int(vec(8))
         ep% logQ_min = vec(9)
         ep% logQ_max = vec(10)
         ep% del_logQ = vec(11)

         read(io_unit,*,iostat=info)
         if (info /= 0) return

         if (abs(X-X_in) > tiny .or. abs(Z-Z_in) > tiny) then
            write(*,*) 'bad header info in ' // trim(filename)
            info = -1
            close(io_unit)
            if (abs(X-X_in) > tiny) then
               write(*,'(a50,l1)') 'abs(X-X_in) > tiny', abs(X-X_in) > tiny
            end if
            if (abs(Z-Z_in) > tiny) then
               write(*,'(a50,l1)') 'abs(Z-Z_in) > tiny', abs(Z-Z_in) > tiny
            end if
            write(*,*)
            call request_user_to_reinstall
            return
         end if
         
         if (show_allocations) write(*,2) 'Load1_eosDT_Table ep% tbl1', &
             sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs + ep% num_logQs + ep% num_logTs
         allocate(ep% tbl1(sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs), &
            ep% logQs(ep% num_logQs), ep% logTs(ep% num_logTs),   &
            STAT=info)
         if (info /= 0) then
            return
         end if
         
         tbl(1:sz_per_eos_point,1:nv,1:ep% num_logQs,1:ep% num_logTs) =>  &
               ep% tbl1(1:sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs)
         
         ep% logQs(1) = ep% logQ_min
         do i = 2, ep% num_logQs-1
            ep% logQs(i) = ep% logQs(i-1) + ep% del_logQ
         end do
         ep% logQs(ep% num_logQs) = ep% logQ_max
         
         ep% logTs(1) = ep% logT_min
         do i = 2, ep% num_logTs-1
            ep% logTs(i) = ep% logTs(i-1) + ep% del_logT
         end do
         ep% logTs(ep% num_logTs) = ep% logT_max

         if (use_cache) then
            call Read_EoS_Cache(X, Z, ep, cache_filename, cache_io_unit, ios)
            if (ios == 0) then
               close(io_unit)
               return
            end if
         end if
                  
         status = 0
         allocate(tbl2_1(num_eos_file_vals*ep% num_logQs*ep% num_logTs), STAT=status)
         if (status .ne. 0) then
            info = -1
            return
         end if
         
         tbl2(1:num_eos_file_vals,1:ep% num_logQs,1:ep% num_logTs) =>  &
               tbl2_1(1:num_eos_file_vals*ep% num_logQs*ep% num_logTs)

         do iQ=1,ep% num_logQs
         
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
                        
            read(io_unit,'(a)',iostat=info) message
            if (info == 0) call str_to_double(message, vec(1), info)
            if (failed('read logQ')) return
            logQ = vec(1)

            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            
            do i=1,ep% num_logTs
            
               read(io_unit,'(a)',iostat=info) message
               if (failed('read line')) then
                  write(*,'(a)') trim(message)
                  write(*,*) trim(filename)
                  write(*,*) 'iQ, i', iQ, i
                  write(*,*) 'logQ', logQ
                  write(*,*) 'bad input line?'
                  call mesa_error(__FILE__,__LINE__)
               end if
               
               call str_to_vector(message, vec, n, info)
               if (info /= 0 .or. n < 1+num_eos_file_vals) then
                  write(*,'(a)') trim(message)
                  write(*,*) trim(filename)
                  write(*,*) 'iQ, i', iQ, i
                  write(*,*) 'logQ', logQ
                  write(*,*) 'bad input line?'
                  call mesa_error(__FILE__,__LINE__)
               end if
               logT = vec(1)
               do j=1,num_eos_file_vals
                  tbl2(j,iQ,i) = vec(1+j)
               end do
               
            enddo
            
            if(iQ == ep% num_logQs) exit
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            
         end do
            
         close(io_unit)
         
         call Make_XEoS_Interpolation_Data(ep, tbl2_1, info)
         deallocate(tbl2_1)
         if (failed('Make_XEoS_Interpolation_Data')) return
         
         call Check_XEoS_Interpolation_Data(ep)
         
         if (.not. use_cache) return

         open(unit=cache_io_unit, &
            file=trim(switch_str(temp_cache_filename, cache_filename, use_mesa_temp_cache)), &
            iostat=ios,action='write', form='unformatted')

         if (ios == 0) then
            write(*,'(a)') 'write ' // trim(cache_filename)
            write(cache_io_unit)  &
               X_in, Z_in, ep% num_logTs, ep% logT_min, ep% logT_max, ep% del_logT,  &
               ep% num_logQs, ep% logQ_min, ep% logQ_max, ep% del_logQ, ep% version
            write(cache_io_unit)  &
               ep% tbl1(&
                  1:sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs)
            close(cache_io_unit)
            if(use_mesa_temp_cache) call mv(temp_cache_filename, cache_filename,.true.)
         end if
         
         
         contains
         
         subroutine Check_XEoS_Interpolation_Data(ep)
            use utils_lib,only:is_bad
            type (EosDT_XZ_Info), pointer :: ep
            ! for logT > 6.8 and logRho < -10, splines can get bogus higher order terms
            ! replace NaN's and Infinities with 0
            integer :: i, j, iQ, jtemp
            do i = 1, sz_per_eos_point
               do j = 1, nv
                  do iQ = 1, ep% num_logQs
                     do jtemp = 1, ep% num_logTs
                        if (is_bad(tbl(i,j,iQ,jtemp))) then
                           tbl(i,j,iQ,jtemp) = 0
                        end if
                     end do
                  end do
               end do
            end do
         end subroutine Check_XEoS_Interpolation_Data
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (info /= 0)
            if (failed) write(*,*) 'Load1_eosDT_Table failed: ' // trim(str)
         end function failed
         

      end subroutine Load1_eosDT_Table
      
      
      subroutine Make_XEoS_Interpolation_Data(ep, tbl2_1, info)
         use interp_2d_lib_db
         use const_def, only: crad, ln10

         type (EosDT_XZ_Info), pointer :: ep
         real(dp), pointer :: tbl2_1(:) ! =(num_eos_file_vals, ep% num_logQs, ep% num_logTs)
         integer, intent(out) :: info

         real(dp) :: logQs(ep% num_logQs)              ! x vector, strict ascending
         real(dp) :: logTs(ep% num_logTs)                    ! y vector, strict ascending
         real(dp) :: Ts(ep% num_logTs)
         real(dp), allocatable, target :: f1_ary(:) ! data & spline coefficients
         real(dp), pointer :: f1(:), f(:,:,:), ep_tbl(:,:,:,:), tbl2(:,:,:)
         integer :: ibcxmin                   ! bc flag for x=xmin
         real(dp) :: bcxmin(ep% num_logTs)    ! bc data vs. y at x=xmin
         integer :: ibcxmax                   ! bc flag for x=xmax
         real(dp) :: bcxmax(ep% num_logTs)     ! bc data vs. y at x=xmax
         integer :: ibcymin                   ! bc flag for y=ymin
         real(dp) :: bcymin(ep% num_logQs)   ! bc data vs. x at y=ymin
         integer :: ibcymax                   ! bc flag for y=ymax
         real(dp) :: bcymax(ep% num_logQs)   ! bc data vs. x at y=ymax
         integer :: ili_logQs    ! =1: logRho grid is "nearly" equally spaced
         integer :: ili_logTs      ! =1: logT grid is "nearly" equally spaced
         integer :: ier            ! =0 on exit if there is no error.
         real(dp) :: logQ, Rho, logRho, T, P, Cv, chiRho, chiT, logT, logT0, logT1, logQ0, logQ1
         real(dp) :: gamma3, gamma1, grad_ad, Prad, E, S
         integer :: iQ, jtemp, ilogT, ilogQ
         real(dp) :: fval(num_eos_file_vals), df_dx(num_eos_file_vals), df_dy(num_eos_file_vals)
         
         real(dp) :: x, y, dlnT, energy, lnE, entropy, lnS, Pgas, lnPgas, dlogT, &
            dlnPgas_dlnd, dlnE_dlnd, dlnS_dlnd, dlnPgas_dlnT, dlnE_dlnT, dlnS_dlnT
         
         integer :: v, vlist(3), var, i, j, num_logQs, num_logTs, ii, jj
         character (len=256) :: message
         
         include 'formats'

         info = 0

         ! just use "not a knot" bc's at edges of tables
         ibcxmin = 0; bcxmin(:) = 0
         ibcxmax = 0; bcxmax(:) = 0
         ibcymin = 0; bcymin(:) = 0
         ibcymax = 0; bcymax(:) = 0
         
         num_logQs = ep% num_logQs
         num_logTs = ep% num_logTs
         
         ep_tbl(1:sz_per_eos_point,1:nv,1:num_logQs,1:num_logTs) =>  &
               ep% tbl1(1:sz_per_eos_point*nv*num_logQs*num_logTs)

         tbl2(1:num_eos_file_vals,1:num_logQs,1:num_logTs) =>  &
               tbl2_1(1:num_eos_file_vals*num_logQs*num_logTs)
         
         allocate(f1_ary(sz_per_eos_point * ep% num_logQs * ep% num_logTs))
         
         f1 => f1_ary
         f(1:sz_per_eos_point,1:num_logQs,1:num_logTs) => &
               f1_ary(1:sz_per_eos_point*num_logQs*num_logTs)
         
         do iQ = 1, ep% num_logQs
            logQs(iQ) = ep% logQ_min + (iQ-1) * ep% del_logQ
         end do

         do jtemp = 1, ep% num_logTs
            logTs(jtemp) = ep% logT_min + (jtemp-1) * ep% del_logT
         end do
         
         ! copy file eos variables to internal eos interpolation tables
         do j=1,num_logTs
            do i=1,num_logQs
               ep_tbl(1,i_lnPgas,i,j) = tbl2(jlogPgas,i,j)*ln10
               ep_tbl(1,i_lnE,i,j) = tbl2(jlogE,i,j)*ln10
               ep_tbl(1,i_lnS,i,j) = tbl2(jlogS,i,j)*ln10
               ep_tbl(1,i_grad_ad,i,j) = tbl2(jgrad_ad,i,j)
               ep_tbl(1,i_chiRho,i,j) = tbl2(jchiRho,i,j)
               ep_tbl(1,i_chiT,i,j) = tbl2(jchiT,i,j)
               ep_tbl(1,i_Cp,i,j) = tbl2(jCp,i,j)
               ep_tbl(1,i_Cv,i,j) = tbl2(jCv,i,j)
               ep_tbl(1,i_dE_dRho,i,j) = tbl2(jdE_dRho,i,j)
               ep_tbl(1,i_dS_dT,i,j) = tbl2(jdS_dT,i,j)
               ep_tbl(1,i_dS_dRho,i,j) = tbl2(jdS_dRho,i,j)
               ep_tbl(1,i_mu,i,j) = tbl2(jmu,i,j)
               ep_tbl(1,i_lnfree_e,i,j) = max(-4d0,tbl2(jlogfree_e,i,j))*ln10
                  ! to protect against non-monotonic interpolation caused by extreme values
               ep_tbl(1,i_gamma1,i,j) = tbl2(jgamma1,i,j)
               ep_tbl(1,i_gamma3,i,j) = tbl2(jgamma3,i,j)
               ep_tbl(1,i_eta,i,j) = tbl2(jeta,i,j)  
            end do             
         end do

         ! create tables for bicubic spline interpolation         
         do v = 1, nv
            do i=1,ep% num_logQs
               do j=1,ep% num_logTs
                  f(1,i,j) = ep_tbl(1,v,i,j)
               end do
            end do
            call interp_mkbicub_db( &
                  logQs,ep% num_logQs,logTs,ep% num_logTs,f1,ep% num_logQs, &
                  ibcxmin,bcxmin,ibcxmax,bcxmax, &
                  ibcymin,bcymin,ibcymax,bcymax, &
                  ili_logQs,ili_logTs,ier)
            if (ier /= 0) then
               write(*,*) 'interp_mkbicub_db error happened for eos_value', v
               info = 3
               return
            end if
            do i=1,ep% num_logQs
               do j=1,ep% num_logTs
                  ep_tbl(2,v,i,j) = f(2,i,j)
                  ep_tbl(3,v,i,j) = f(3,i,j)
                  ep_tbl(4,v,i,j) = f(4,i,j)
               end do
            end do
         end do
         
         
      end subroutine Make_XEoS_Interpolation_Data
      
      
      subroutine Read_EoS_Cache(X, Z, ep, cache_filename, io_unit, ios)
         real(dp), intent(in) :: X, Z
         type (EosDT_XZ_Info), pointer :: ep
         character (*), intent(in) :: cache_filename
         integer, intent(in) :: io_unit ! use this for file access
         integer, intent(out) :: ios

         real(dp) :: X_in, Z_in, logT_min_in, logT_max_in, del_logT_in,  &
               logQ_min_in, logQ_max_in, del_logQ_in
         integer :: num_logQs_in, num_logTs_in, version_in, i, j
         real(dp), parameter :: tiny = 1d-10
         
         include 'formats'
         
         ios = 0
         open(unit=io_unit,file=trim(cache_filename),action='read', &
               status='old',iostat=ios,form='unformatted')
         if (ios /= 0) return
         
         read(io_unit, iostat=ios)  &
               X_in, Z_in, num_logTs_in, logT_min_in, logT_max_in, del_logT_in,  &
               num_logQs_in, logQ_min_in, logQ_max_in, del_logQ_in, version_in
         if (ios /= 0) return
         
         if (ep% version /= version_in) then
            ios = 1
            write(*,*) 'read cache failed for version_in'
         end if
         if (ep% num_logQs /= num_logQs_in) then
            ios = 1
            write(*,*) 'read cache failed for ep% num_logQs'
         end if 
         if (ep% num_logTs /= num_logTs_in) then
            ios = 1
            write(*,*) 'read cache failed for ep% num_logTs'
         end if
         if (abs(X-X_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for X_in'
         end if
         if (abs(Z-Z_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for Z_in'
         end if
         if (abs(ep% logT_min-logT_min_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for eos_logT_min'
         end if    
         if (abs(ep% logT_max-logT_max_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for eos_logT_max'
         end if    
         if (abs(ep% del_logT-del_logT_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for eos_del_logT'
         end if    
         if (abs(ep% logQ_min-logQ_min_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for eos_logQ_min'
         end if    
         if (abs(ep% logQ_max-logQ_max_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for eos_logQ_max'
         end if
         if (abs(ep% del_logQ-del_logQ_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for eos_del_logQ'
         end if
         
         if (ios /= 0) then
            close(io_unit); return
         end if

         read(io_unit, iostat=ios)  &
            ep% tbl1( &
               1:sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs)
         if (ios /= 0) then
            close(io_unit); return
         end if
         
         close(io_unit)

      end subroutine Read_EoS_Cache
      
      
      end module eosDT_load_tables
