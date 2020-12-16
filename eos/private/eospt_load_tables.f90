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

      module eosPT_load_tables
      use eos_def
      use const_def
      use math_lib
      use utils_lib, only: mesa_error, switch_str, mv
     
      implicit none

      
      ! the file eosPT data
      integer, parameter :: jlogRho = 1
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
      integer, parameter :: num_eosPT_file_vals = 16

      integer, parameter :: file_max_num_logWs = 1000
      

      contains
      
      
      subroutine check_for_error_in_eosPT_data(ierr, fname)
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
      end subroutine check_for_error_in_eosPT_data
      
      
      subroutine request_reinstall_eosPT_data
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*) 'NOTICE: you need to install a new verion of the eos data.'
            write(*,*) 'Please update by removing the directory mesa/data/eosPT_data,'
            write(*,*) 'and rerunning the mesa ./install script.'
            write(*,*)
            write(*,*)
            call mesa_error(__FILE__,__LINE__)
      end subroutine request_reinstall_eosPT_data

   
   
      subroutine Load_single_eosPT_Table(rq,ep,x,z,load_ix,load_iz,ierr)
         type (EoS_General_Info), pointer :: rq
         type (EosPT_XZ_Info), pointer :: ep
         real(dp),intent(in) :: x,z
         integer, intent(out) :: load_ix,load_iz,ierr
         
         integer :: iz, ix, i
         include 'formats'
         ierr = 0

         load_iz=0
         do iz = 1, num_eosPT_Zs-1
            if(z.ge.eosPT_Zs(iz) .and. z.le.eosPT_Zs(iz+1))then
               load_iz=iz
            end if
         end do  
         if(z.le.eosPT_Zs(1)) then
            load_iz=1
         end if
         if(z.ge.eosPT_Zs(num_eosPT_Zs)) then
            load_iz=num_eosPT_Zs
         end if
            
         load_ix=0         
         do ix = 1, num_eosPT_Xs-1
            if(x.ge.eosPT_Xs(ix) .and. x.le.eosPT_Xs(ix+1))then
               load_ix=ix
            end if
         end do  
         
         if(x.le.eosPT_Xs(1)) then
            load_ix=1
         end if
         if(x.ge.eosPT_Xs(num_eosPT_Xs))then
            load_ix=num_eosPT_Xs
         end if
         call load_single_eosPT_table_by_id(rq,ep,load_ix,load_iz,ierr)
         
    end subroutine Load_single_eosPT_Table


    subroutine load_single_eosPT_table_by_id(rq,ep,ix,iz,ierr)
         use utils_lib
         type (EoS_General_Info), pointer :: rq
         type (EosPT_XZ_Info), pointer :: ep
         integer,intent(in) :: iz, ix
         integer, intent(out) :: ierr         
      
      if (rq% use_max_SCVH_for_PT) then
         ep => eosSCVH_PT_XZ_data(ix,iz)
         if (eosSCVH_PT_XZ_loaded(ix,iz)) return
      else
         ep => eosPT_XZ_data(ix,iz)
         if (eosPT_XZ_loaded(ix,iz)) return
      end if
       
!$OMP CRITICAL(eosPT_load)
      if (rq% use_max_SCVH_for_PT) then
         if (.not. eosSCVH_PT_XZ_loaded(ix,iz)) call do_read
      else
         if (.not. eosPT_XZ_loaded(ix,iz)) call do_read
      end if
!$OMP END CRITICAL(eosPT_load)
         
         contains
         
         subroutine do_read
            call read_one(ix,iz,ierr)
            if (ierr /= 0) return
            if (rq% use_max_SCVH_for_PT) then
               eosSCVH_PT_XZ_loaded(ix,iz) = .true.
            else
               eosPT_XZ_loaded(ix,iz) = .true.
            end if
         end subroutine do_read
         
         subroutine read_one(ix,iz,ierr)
            use const_def, only: mesa_data_dir
            integer, intent(in) :: ix, iz
            integer, intent(out) :: ierr
            character (len=256) :: fname, cache_filename, temp_cache_filename
            integer :: iounit1, iounit2            
            iounit1 = alloc_iounit(ierr); if (ierr /= 0) return
            iounit2 = alloc_iounit(ierr); if (ierr /= 0) return
            call Get_eosPT_Table_Filenames(rq, &
               eosPT_Zs(iz), eosPT_Xs(ix), mesa_data_dir, &
               fname, cache_filename, temp_cache_filename)
            call Load1_eosPT_Table(rq, ep, &
               ix, iz, fname, cache_filename, temp_cache_filename, iounit1, iounit2, &
               use_cache_for_eos, ierr)
            if (ierr /= 0) then
               write(*,*) 'Load_eosPT_Table ierr', ierr, ix, iz
            end if
            call free_iounit(iounit2)
            call free_iounit(iounit1)        
         end subroutine read_one
         
    end subroutine load_single_eosPT_table_by_id
         
      
      subroutine Get_eosPT_Table_Filenames(rq, &
            Z, X, data_dir, fname, cache_filename, temp_cache_filename)
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X
         character (*), intent(in) :: data_dir
         character (len=*), intent(out) :: fname, cache_filename, temp_cache_filename
         character (len=256) :: Zstr, Xstr, suffix, data_dir_name, data_prefix
         
         call setstr(Z,Zstr)
         call setstr(X,Xstr)
         !if (Zstr == '100') then
         !   suffix = eosPT_Z1_suffix
         !else
            suffix = ''
         !end if
         
         if (rq% use_max_SCVH_for_PT) then
            data_dir_name = '/eosSCVH_PT_data/'
            data_prefix = '-SCVH_PT_'
         else
            data_dir_name = '/eosPT_data/'
            data_prefix = '-eosPT_'
         end if
         
         fname = trim(data_dir) // trim(data_dir_name) // &
            trim(rq% eosPT_file_prefix) // trim(data_prefix) // &
            trim(Zstr) // 'z' // trim(Xstr) // 'x' // trim(suffix) // '.data'
         cache_filename = trim(eosPT_cache_dir) //  &
            '/' // trim(rq% eosPT_file_prefix) // trim(data_prefix) // &
            trim(Zstr) // 'z' // trim(Xstr) // 'x' // trim(suffix) // '.bin'
         temp_cache_filename = trim(eosPT_temp_cache_dir) //  &
               '/' // trim(rq% eosPT_file_prefix) // trim(data_prefix) // &
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

      end subroutine Get_eosPT_Table_Filenames
      
      
      subroutine Load1_eosPT_Table( &
            rq, ep, ix, iz, filename, cache_filename, temp_cache_filename, io_unit, cache_io_unit, &
            use_cache, info)
         type (EoS_General_Info), pointer :: rq
         type (EosPT_XZ_Info), pointer :: ep
         integer, intent(in) :: ix, iz
         character (*), intent(in) :: filename, cache_filename, temp_cache_filename
         integer, intent(in) :: io_unit, cache_io_unit
         logical, intent(in) :: use_cache
         integer, intent(out) :: info
         
         real(dp) :: X, Z
         real(dp) :: X_in, Z_in, logW, logT
         integer :: j, i, k, iW, ios, status
         character (len=1000) :: message
         real(dp), parameter :: tiny = 1d-10
         real(dp), pointer :: tbl(:,:,:,:) ! => ep% tbl1
         real(dp), pointer :: tbl2_1(:), tbl2(:,:,:)
         real(dp), target :: vec_ary(50)
         real(dp), pointer :: vec(:)
         integer :: n

         info = 0            
         vec => vec_ary
         
         X = eosPT_Xs(ix)
         Z = eosPT_Zs(iz)   
         
         write(message,*) 'open ', trim(filename)
         
         open(UNIT=io_unit, FILE=trim(filename), ACTION='READ', STATUS='OLD', IOSTAT=ios)
         call check_for_error_in_eosPT_data(ios, filename)

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
         ep% num_logWs = int(vec(8))
         ep% logW_min = vec(9)
         ep% logW_max = vec(10)
         ep% del_logW = vec(11)

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
            call request_reinstall_eosPT_data
            return
         end if

         allocate(ep% tbl1(sz_per_eos_point*nv*ep% num_logWs*ep% num_logTs),  &
                  ep% logWs(ep% num_logWs), ep% logTs(ep% num_logTs),   &
                  STAT=info)
         if (info /= 0) return
         
         tbl(1:sz_per_eos_point,1:nv,1:ep% num_logWs,1:ep% num_logTs) =>  &
               ep% tbl1(1:sz_per_eos_point*nv*ep% num_logWs*ep% num_logTs)

         ep% logWs(1) = ep% logW_min
         do i = 2, ep% num_logWs-1
            ep% logWs(i) = ep% logWs(i-1) + ep% del_logW
         end do
         ep% logWs(ep% num_logWs) = ep% logW_max
         
         ep% logTs(1) = ep% logT_min
         do i = 2, ep% num_logTs-1
            ep% logTs(i) = ep% logTs(i-1) + ep% del_logT
         end do
         ep% logTs(ep% num_logTs) = ep% logT_max

         if (use_cache) then
            call Read_eosPT_Cache(X, Z, ep, cache_filename, cache_io_unit, ios)
            if (ios == 0) then
               close(io_unit)
               return
            end if
         end if
         
         status = 0
         allocate(tbl2_1(num_eosPT_file_vals*ep% num_logWs*ep% num_logTs), STAT=status)
         if (status .ne. 0) then
            info = -1
            return
         end if
         
         tbl2(1:num_eosPT_file_vals,1:ep% num_logWs,1:ep% num_logTs) =>  &
               tbl2_1(1:num_eosPT_file_vals*ep% num_logWs*ep% num_logTs)
         
         do iW=1,ep% num_logWs
         
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            
            read(io_unit,'(a)',iostat=info) message
            if (info == 0) call str_to_double(message, vec(1), info)
            if (failed('read logW')) return
            logW = vec(1)
            
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return 
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return
            
            do i=1,ep% num_logTs

               read(io_unit,'(a)',iostat=info) message
               if (failed('read line')) then
                  write(*,'(a)') trim(message)
                  write(*,*) trim(filename)
                  write(*,*) 'iW, i', iW, i
                  write(*,*) 'logW', logW
                  write(*,*) 'bad input line?'
                  call mesa_error(__FILE__,__LINE__)
               end if
               
               call str_to_vector(message, vec, n, info)
               if (info /= 0 .or. n < 1+num_eosPT_file_vals) then
                  write(*,'(a)') trim(message)
                  write(*,*) trim(filename)
                  write(*,*) 'iW, i', iW, i
                  write(*,*) 'logW', logW
                  write(*,*) 'bad input line?'
                  call mesa_error(__FILE__,__LINE__)
               end if
               logT = vec(1)
               do j=1,num_eosPT_file_vals
                  tbl2(j,iW,i) = vec(1+j)
               end do
               
            enddo
            
            if (iW == ep% num_logWs) exit
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return 
            read(io_unit,*,iostat=info)
            if (failed('skip line')) return 
            
         end do
            
         close(io_unit)
         
         call Make_XeosPT_Interpolation_Data(ep, tbl2_1, info)
         deallocate(tbl2_1)
         if (info /= 0) return
         
         call Check_XeosPT_Interpolation_Data(ep)
         
         if (.not. use_cache) return

         open(unit=cache_io_unit, &
            file=trim(switch_str(temp_cache_filename, cache_filename, use_mesa_temp_cache)), &
            iostat=ios,action='write', form='unformatted')

         if (ios == 0) then
            write(*,'(a)') 'write ' // trim(cache_filename)
            write(cache_io_unit)  &
               X_in, Z_in, ep% num_logTs, ep% logT_min, ep% logT_max, ep% del_logT,  &
               ep% num_logWs, ep% logW_min, ep% logW_max, ep% del_logW, ep% version
            write(cache_io_unit)  &
               ep% tbl1(1:sz_per_eos_point*nv*ep% num_logWs*ep% num_logTs)
            close(cache_io_unit)
            if(use_mesa_temp_cache) call mv(temp_cache_filename, cache_filename,.true.)
         else
            write(*,*) 'open failed for ' // trim(cache_filename)
         end if
         
         contains
         
         subroutine Check_XeosPT_Interpolation_Data(ep)
            use utils_lib,only:is_bad
            type (EosPT_XZ_Info), pointer :: ep
            ! for logT > 6.8 and logRho < -10, splines can get bogus higher order terms
            ! replace NaN's and Infinities with 0
            integer :: i, j, iW, jtemp
            do i = 1, sz_per_eos_point
               do j = 1, nv
                  do iW = 1, ep% num_logWs
                     do jtemp = 1, ep% num_logTs
                        if (is_bad(tbl(i,j,iW,jtemp))) then
                           tbl(i,j,iW,jtemp) = 0
                        end if
                     end do
                  end do
               end do
            end do
         end subroutine Check_XeosPT_Interpolation_Data
         
         
          logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (info /= 0)
            if (failed) write(*,*) 'Load1_eosPT_Table failed: ' // trim(str)
         end function failed


      end subroutine Load1_eosPT_Table
      
      
      subroutine Make_XeosPT_Interpolation_Data(ep, tbl2_1, info)
         use interp_2d_lib_db

         type (EosPT_XZ_Info), pointer :: ep
         real(dp), pointer :: tbl2_1(:) ! =(num_eos_file_vals, ep% num_logWs, ep% num_logTs)
         integer, intent(out) :: info

         real(dp) :: logWs(ep% num_logWs)              ! x vector, strict ascending
         real(dp) :: logTs(ep% num_logTs)                    ! y vector, strict ascending
         real(dp) :: Ts(ep% num_logTs)
         real(dp), allocatable, target :: f1_ary(:) ! data & spline coefficients
         real(dp), pointer :: f1(:), f(:,:,:), ep_tbl(:,:,:,:), tbl2(:,:,:)
         integer :: ibcxmin                   ! bc flag for x=xmin
         real(dp) :: bcxmin(ep% num_logTs)    ! bc data vs. y at x=xmin
         integer :: ibcxmax                   ! bc flag for x=xmax
         real(dp) :: bcxmax(ep% num_logTs)     ! bc data vs. y at x=xmax
         integer :: ibcymin                   ! bc flag for y=ymin
         real(dp) :: bcymin(ep% num_logWs)   ! bc data vs. x at y=ymin
         integer :: ibcymax                   ! bc flag for y=ymax
         real(dp) :: bcymax(ep% num_logWs)   ! bc data vs. x at y=ymax
         integer :: ili_logWs    ! =1: logRho grid is "nearly" equally spaced
         integer :: ili_logTs      ! =1: logT grid is "nearly" equally spaced
         integer :: ier            ! =0 on exit if there is no error.
         real(dp) :: logW, logPgas, Pgas, Prad, Rho, logRho, T, P, Cv, chiRho, chiT
         real(dp) :: gamma3, gamma1, grad_ad, E, S, logT, logT0, logT1, logW0, logW1
         integer :: iW, jtemp, ilogT, ilogW, num_logWs, num_logTs
         real(dp) :: fval(num_eosPT_file_vals), df_dx(num_eosPT_file_vals), df_dy(num_eosPT_file_vals)
         
         integer :: v, vlist(3), var, i, j
         character (len=256) :: message

         info = 0
         
         num_logWs = ep% num_logWs
         num_logTs = ep% num_logTs
         
         ep_tbl(1:sz_per_eos_point,1:nv,1:num_logWs,1:num_logTs) =>  &
               ep% tbl1(1:sz_per_eos_point*nv*num_logWs*num_logTs)

         tbl2(1:num_eosPT_file_vals,1:num_logWs,1:num_logTs) =>  &
               tbl2_1(1:num_eosPT_file_vals*num_logWs*num_logTs)
         
         allocate(f1_ary(sz_per_eos_point * ep% num_logWs * ep% num_logTs))
         
         f1 => f1_ary
         f(1:sz_per_eos_point,1:num_logWs,1:num_logTs) => &
               f1_ary(1:sz_per_eos_point*num_logWs*num_logTs)
                  
         do iW = 1, num_logWs
            logWs(iW) = ep% logW_min + (iW-1) * ep% del_logW
         end do

         do jtemp = 1, num_logTs
            logTs(jtemp) = ep% logT_min + (jtemp-1) * ep% del_logT
         end do
         
         ! copy file eos variables to internal eos interpolation tables
         do j=1,num_logTs
            do i=1,num_logWs
               ep_tbl(1,eosPT_ilnRho,i,j) = tbl2(jlogRho,i,j)*ln10
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

         ! just use "not a knot" bc's at edges of tables
         ibcxmin = 0; bcxmin(:) = 0
         ibcxmax = 0; bcxmax(:) = 0
         ibcymin = 0; bcymin(:) = 0
         ibcymax = 0; bcymax(:) = 0

         ! create tables for bicubic spline interpolation         
         do v = 1, nv
            do i=1,ep% num_logWs
               do j=1,ep% num_logTs
                  f(1,i,j) = ep_tbl(1,v,i,j)
               end do
            end do
            call interp_mkbicub_db( &
               logWs,ep% num_logWs,logTs,ep% num_logTs,f1,ep% num_logWs, &
               ibcxmin,bcxmin,ibcxmax,bcxmax, &
               ibcymin,bcymin,ibcymax,bcymax, &
               ili_logWs,ili_logTs,ier)
            if (ier /= 0) then
               write(*,*) 'interp_mkbicub_db error happened for ep% value', v
               info = 3
               return
            end if
            do i=1,ep% num_logWs
               do j=1,ep% num_logTs
                  ep_tbl(2,v,i,j) = f(2,i,j)
                  ep_tbl(3,v,i,j) = f(3,i,j)
                  ep_tbl(4,v,i,j) = f(4,i,j)
               end do
            end do
         end do
         
      end subroutine Make_XeosPT_Interpolation_Data
      
      
      subroutine Read_eosPT_Cache(X, Z, ep, cache_filename, io_unit, ios)
         real(dp), intent(in) :: X, Z
         type (EosPT_XZ_Info), pointer :: ep
         character (*), intent(in) :: cache_filename
         integer, intent(in) :: io_unit ! use this for file access
         integer, intent(out) :: ios

         real(dp) :: X_in, Z_in, logT_min_in, logT_max_in, del_logT_in,  &
               logW_min_in, logW_max_in, del_logW_in
         integer :: num_logWs_in, num_logTs_in, version_in, i, j
         real(dp), parameter :: tiny = 1d-10
         
         ios = 0
         open(unit=io_unit,file=trim(cache_filename),action='read', &
               status='old',iostat=ios,form='unformatted')
         if (ios /= 0) return
         
         read(io_unit, iostat=ios)  &
               X_in, Z_in, num_logTs_in, logT_min_in, logT_max_in, del_logT_in,  &
               num_logWs_in, logW_min_in, logW_max_in, del_logW_in, version_in
         if (ios /= 0) return
         
         if (ep% version /= version_in) then
            ios = 1
            write(*,*) 'read cache failed for version_in'
         end if
         if (ep% num_logWs /= num_logWs_in) then
            ios = 1
            write(*,*) 'read cache failed for ep% num_logWs'
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
            write(*,*) 'read cache failed for ep% logT_min'
         end if    
         if (abs(ep% logT_max-logT_max_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ep% logT_max'
         end if    
         if (abs(ep% del_logT-del_logT_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ep% del_logT'
         end if    
         if (abs(ep% logW_min-logW_min_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ep% logW_min'
         end if    
         if (abs(ep% logW_max-logW_max_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ep% logW_max'
         end if
         if (abs(ep% del_logW-del_logW_in) > tiny) then
            ios = 1
            write(*,*) 'read cache failed for ep% del_logW'
         end if
         
         if (ios /= 0) then
            close(io_unit); return
         end if
         
         read(io_unit, iostat=ios)  &
            ep% tbl1(1:sz_per_eos_point*nv*ep% num_logWs*ep% num_logTs)
         if (ios /= 0) then
            write(*,*) 'read failed for ' // trim(cache_filename)
            close(io_unit); return
         end if
         
         close(io_unit)

      end subroutine Read_eosPT_Cache
      
      
      end module eosPT_load_tables
