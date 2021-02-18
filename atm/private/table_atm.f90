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

      module table_atm
      use atm_def
      use const_def, only: dp
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none


      logical, parameter :: dbg = .false.


      contains



      !reads in table_summary file from atm_data, initializes logZ, Teff_array, 
      ! logg_array, and Teff_bound arrays; sets some flags
      subroutine table_atm_init(use_cache, ierr)
         implicit none
         logical, intent(in) :: use_cache
         integer, intent(out) :: ierr
      
         integer :: nZ, ng, nT, i, j, iounit
         integer, pointer :: ibound(:,:), tmp_version(:)
         character(len=256) :: filename
      
         if (table_atm_is_initialized) call table_atm_shutdown()
         
         ierr = 0

         call load_table_summary(ATM_TABLE_PHOTOSPHERE, 'table_summary.txt', ai_two_thirds, ierr)
         if (ierr /= 0) return         
         call load_table_summary(ATM_TABLE_TAU_100, 'table100_summary.txt', ai_100, ierr)
         if (ierr /= 0) return
         call load_table_summary(ATM_TABLE_TAU_10, 'table10_summary.txt', ai_10, ierr)
         if (ierr /= 0) return
         call load_table_summary(ATM_TABLE_TAU_1, 'table1_summary.txt', ai_1, ierr)
         if (ierr /= 0) return
         call load_table_summary(ATM_TABLE_TAU_1M1, 'table1m1_summary.txt', ai_1m1, ierr)
         if (ierr /= 0) return
         call load_table_summary(ATM_TABLE_WD_TAU_25, 'table_wd_25_summary.txt', ai_wd_25, ierr)
         if (ierr /= 0) return
         call load_table_summary(ATM_TABLE_DB_WD_TAU_25, 'table_db_wd_25_summary.txt', ai_db_wd_25, ierr)
         if (ierr /= 0) return
         
         table_atm_is_initialized = .true.
         
         
         contains
      
      
         subroutine load_table_summary(id, fname, ai, ierr)
            use const_def, only: mesa_data_dir
            integer, intent(in) :: id
            character(len=*), intent(in) :: fname
            type(atm_info), intent(inout) :: ai
            integer, intent(out) :: ierr
            
            integer :: nvec
            character (len=500) :: buf
            real(dp), target :: vec_ary(20)
            real(dp), pointer :: vec(:)

            vec => vec_ary
      
            filename = trim(mesa_data_dir)//'/atm_data/' // trim(fname)
            
            if (dbg) write(*,*) 'read ' // trim(filename)
            
            open(newunit=iounit,file=trim(filename),action='read',status='old',iostat=ierr)
            if (ierr/= 0) then
               write(*,*) 'table_atm_init: missing atm data'
               write(*,*) trim(filename)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*) 'FATAL ERROR: missing or bad atm data.'
               write(*,*) 'Please update by removing the directory mesa/data/atm_data,'
               write(*,*) 'and rerunning the mesa ./install script.'
               write(*,*)
               call mesa_error(__FILE__,__LINE__)
            endif

            !read first line and (nZ, nT, ng)
            read(iounit,*)            !first line is text, skip it
            read(iounit,*) nZ, nT, ng
      
            ai% nZ = nZ
            ai% nT = nT
            ai% ng = ng
            ai% id = id

            allocate( &
               ai% Teff_array(nT), ai% logg_array(ng), ai% Teff_bound(ng), &
               ai% logZ(nZ), ai% alphaFe(nZ), &
               ai% Pgas_interp1(4*ng*nT*nZ), ai% T_interp1(4*ng*nT*nZ), &
               ai% have_atm_table(nZ), ai% atm_mix(nZ), ai% table_atm_files(nZ))
            
            ai% Pgas_interp(1:4,1:ng,1:nT,1:nZ) => ai% Pgas_interp1(1:4*ng*nT*nZ)
            ai% T_interp(1:4,1:ng,1:nT,1:nZ) => ai% T_interp1(1:4*ng*nT*nZ)

            allocate(ibound(ng,nZ), tmp_version(nZ))

            ai% have_atm_table(:) = .false.

            !read filenames and headers
            read(iounit,*)            !text
            do i=1,nZ
               read(iounit,'(a)') ai% table_atm_files(i)
               read(iounit,'(14x,i4)') tmp_version(i)
               read(iounit,1) ai% logZ(i), ai% alphaFe(i), ai% atm_mix(i), ibound(1:ng,i)
            enddo

            !read Teff_array
            read(iounit,*)            !text
            read(iounit,2) ai% Teff_array(:)

            !read logg_array
            read(iounit,*)            !text
            read(iounit,3) ai% logg_array(:)
            
            close(iounit)

 1          format(13x,f5.2,8x,f4.1,1x,a8,1x,15x,99i4)
 2          format(13f7.0)
 3          format(13f7.2)

            !determine table boundaries
            do i=1,ng                 ! -- for each logg, smallest Teff at which Pgas->0
               ai% Teff_bound(i) = ai% Teff_array(ibound(i,1))
               do j=2,nZ
                  ai% Teff_bound(i) = min( ai% Teff_bound(i) , ai% Teff_array(ibound(i,j)) )
               enddo
            enddo
            
            
            if (dbg) write(*,*) 'ai% logg_array(:)', ai% logg_array(:)


            deallocate(ibound, tmp_version)
         
         end subroutine load_table_summary


      end subroutine table_atm_init


      subroutine table_atm_shutdown()

         if (.NOT. table_atm_is_initialized) return
      
         call free_table_summary(ai_two_thirds)
         call free_table_summary(ai_100)
         call free_table_summary(ai_10)
         call free_table_summary(ai_1)
         call free_table_summary(ai_1m1)
         call free_table_summary(ai_wd_25)
         call free_table_summary(ai_db_wd_25)

         table_atm_is_initialized = .false.

       contains

         subroutine free_table_summary(ai)

            type(atm_info), intent(inout) :: ai

            deallocate( &
               ai% Teff_array, ai% logg_array, ai% Teff_bound, &
               ai% logZ, ai% alphaFe, &
               ai% Pgas_interp1, ai% T_interp1, &
               ai% have_atm_table, ai% atm_mix, ai% table_atm_files)
            
          end subroutine free_table_summary

      end subroutine table_atm_shutdown

      
      !interpolate in Z, logg, & Teff: 4pt in Z; bicubic spline in Teff,logg
      subroutine get_table_values( &
            id, newZ, newlogg_in, newTeff_in, &
            newPgas, dPgas_dTeff, dPgas_dlogg, &
            newT, dT_dTeff, dT_dlogg, &
            ierr)
         use chem_def,only:GN93_Zsol
         use interp_1d_lib, only: interp_pm, interp_values
         use interp_1d_def
         use interp_2d_lib_db, only: interp_evbicub_db
         use utils_lib, only: is_bad
         
         implicit none
      
         integer, intent(in) :: id
         real(dp), intent(in) :: newZ, newlogg_in, newTeff_in
         real(dp), intent(out) :: newPgas, dPgas_dTeff, dPgas_dlogg
         real(dp), intent(out) :: newT, dT_dTeff, dT_dlogg
         integer, intent(out) :: ierr
         integer :: i, j, numZs, nZ, ng, nT, ict(6), Zindex, Zlo, Zhi
         integer, parameter :: num_res = 3 !number of results (Pgas, dPgas_dTeff, dPgas_dlogg)
         real(dp) :: newlogg, newTeff, newlogZ(1), result_Z(1)
         real(dp), target :: fZ1_ary(4*4)
         real(dp), pointer :: fZ1(:), fZ(:,:)
         real(dp), target :: work_ary(4*mp_work_size)
         real(dp), pointer :: result_2D(:,:), work(:)
         character(len=256) :: filename
         logical :: clip_Teff, clip_logg, gtv_dbg

         type (Atm_Info), pointer :: ai
         
         include 'formats'
         
         gtv_dbg = dbg
         
         fZ1 => fZ1_ary
         fZ(1:4,1:4) => fZ1(1:4*4)
         
         ierr = 0
         work => work_ary

         ! ensure data have been initialized
         if (.not.table_atm_is_initialized) then
            write(*,*) 'get_table_values: table_atm_init required to proceed'
            ierr=-1
            return
         endif

         newlogZ(1) = safe_log10(newZ/GN93_Zsol)

         select case (id)
         case (ATM_TABLE_PHOTOSPHERE)
            ai => ai_two_thirds
         case (ATM_TABLE_TAU_100)
            ai => ai_100
         case (ATM_TABLE_TAU_10)
            ai => ai_10
         case (ATM_TABLE_TAU_1)
            ai => ai_1
         case (ATM_TABLE_TAU_1M1)
            ai => ai_1m1
         case (ATM_TABLE_WD_TAU_25)
            ai => ai_wd_25
         case (ATM_TABLE_DB_WD_TAU_25)
            ai => ai_db_wd_25
         case default
            write(*,*) 'Invalid id in get_table_values:', id
            ierr = -1
            return
         end select
         
         nZ = ai% nZ
         nT = ai% nT
         ng = ai% ng
         
         clip_Teff = .false.
         if (newTeff_in < ai% Teff_array(1)) then
            newTeff = ai% Teff_array(1) ! clip to table in T
            clip_Teff = .true.
         else if (newTeff_in > ai% Teff_array(nT)) then
            newTeff = ai% Teff_array(nT) ! clip to table in T
            clip_Teff = .true.
         else
            newTeff = newTeff_in
         end if
         
         clip_logg = .false.
         if (newlogg_in < ai% logg_array(1)) then
            newlogg = ai% logg_array(1) ! clip to table in logg
            clip_logg = .true.
         else if (newlogg_in > ai% logg_array(ng)) then
            newlogg = ai% logg_array(ng) ! clip to table in logg
            clip_logg = .true.
         else
            newlogg = newlogg_in
         end if

         allocate(result_2D(6,nZ))

         if (gtv_dbg) write(*,*) 'loaded tables: ', ai% have_atm_table(:)

         if (.false.) then !check for Z within range of tables
            if (nz > 1 .and. (newlogZ(1) < ai% logZ(1) .or. newlogZ(1) > ai% logZ(nZ))) then
               write(*,*) 'get_table_values: Z outside range of atm tables'
               ierr=-1
               return
            endif
         else if (newlogZ(1) < ai% logZ(1)) then
            newlogZ(1) = ai% logZ(1)
         else if (newlogZ(1) > ai% logZ(nZ)) then
            newlogZ(1) = ai% logZ(nZ)
         end if

         !identify surrounding Z values and initialize interpolation
         if (newlogZ(1) <= ai% logZ(1)) then
            Zindex = 1
         else if (newlogZ(1) >= ai% logZ(nz)) then
            Zindex = nz
         else
            Zindex = nz
            do j=1,nz-1
               if (ai% logZ(j) <= newlogZ(1)) then
                  Zindex = j
               end if
            end do
         end if

         !identify upper (Zhi) and lower (Zlo) bounds on Z
         Zlo = max( 1,Zindex-1)
         Zhi = min(nZ,Zindex+2)
         if (Zhi-Zlo /= 3) Zlo=max( 1,Zlo-1)
         if (Zhi-Zlo /= 3) Zhi=min(nZ,Zhi+1)
         if (gtv_dbg) write(*,*) 'Zlo, Zindex, Zhi = ', Zlo, Zindex, Zhi
         numZs = Zhi - Zlo + 1

         !do 2D Teff,logg interpolation on each of 4 Z tables (Zlo-Zhi)
                 !  P(g,T) dP_dg, dP_dT, d2P_dg2, d2P_dT2, d2P_dg_dT
                 !    1      2      3       4        5         6
         ict(:) = (/  1,     1,     1,      0,       0,        0  /)

         if (gtv_dbg) write(*,*) 'do_interp for Pgas', id
         call do_interp(ai% Pgas_interp1, newPgas, dPgas_dlogg, dPgas_dTeff, ierr)
         if (ierr /= 0) then
            if (gtv_dbg) write(*,*) 'do_interp failed'
            return
         end if
         if (clip_logg) dPgas_dlogg = 0
         if (clip_Teff) dPgas_dTeff = 0
         
         if (id == ATM_TABLE_PHOTOSPHERE) then
            newT = newTeff
            if (clip_Teff) then
               dT_dTeff = 0
            else
               dT_dTeff = 1
            end if
            dT_dlogg = 0
            return 
         end if

         if (gtv_dbg) write(*,*) 'do_interp for Temp', id
         call do_interp(ai% T_interp1, newT, dT_dlogg, dT_dTeff, ierr)
         if (ierr /= 0) then
            if (gtv_dbg) write(*,*) 'do_interp failed'
            return
         end if
         if (clip_logg) dT_dlogg = 0
         if (clip_Teff) dT_dTeff = 0
         
         if (dbg .or. is_bad(newPgas) .or. is_bad(newT)) then
            write(*,1) 'newPgas', newPgas
            write(*,1) 'dPgas_dTeff', dPgas_dTeff
            write(*,1) 'dPgas_dlogg', dPgas_dlogg
            write(*,*)
            write(*,1) 'newT', newT
            write(*,1) 'dT_dTeff', dT_dTeff
            write(*,1) 'dT_dlogg', dT_dlogg
            write(*,*)
            if (is_bad(newPgas) .or. is_bad(newT)) stop 'get_table_values'
         end if
         
         !if (dbg) write(*,*) 'loaded tables: ', ai% have_atm_table(:)
      
         deallocate(result_2D)

         contains

         subroutine do_interp(f1, newval, dval_dlogg, dval_dTeff, ierr)
            real(dp), dimension(:), pointer :: f1
            real(dp), intent(out) :: newval, dval_dlogg, dval_dTeff
            integer, intent(out) :: ierr
            
            real(dp) :: res(6)
            integer :: j
            real(dp), pointer :: f(:)
            
            include 'formats'
            
            ierr = 0
            
            do i = Zlo, Zhi
               if (.not. ai% have_atm_table(i)) then
                  call load_atm_table(i,ierr) !<-load on demand
               end if
               if (ierr /= 0) then
                  if (gtv_dbg) write(*,*) 'load_atm_table failed'
                  return
               end if
               
               f(1:4*ng*nT) => f1(1+4*ng*nT*(i-1):4*ng*nT*i)
               call interp_evbicub_db(newlogg, newTeff, ai% logg_array, ng, ai% Teff_array, nT, &
                  ai% iling, ai% ilinT, f, ng, ict, res, ierr)
               do j=1,6 
                  result_2D(j,i) = res(j)
               end do
               if (ierr /= 0) then
                  if (gtv_dbg) write(*,*) 'interp_evbicub_db failed'
                  if (newlogg < ai% logg_array(1)) then
                     write(*,1) 'logg too small for atm table', newlogg, ai% logg_array(1)
                  end if
                  if (newlogg > ai% logg_array(ng)) then
                     write(*,1) 'logg too large for atm table', newlogg, ai% logg_array(ng)
                  end if
                  if (newTeff < ai% Teff_array(1)) then
                     write(*,1) 'Teff too small for atm table', newTeff, ai% Teff_array(1)
                  end if
                  if (newTeff > ai% Teff_array(ng)) then
                     write(*,1) 'Teff too large for atm table', newTeff, ai% Teff_array(ng)
                  end if
                  return
               end if
            enddo

            ! now we have val, dval_dTeff, and dval_dlogg in result_2D for each Z
         
            if (numZs == 1) then
               
               newval = result_2D(1,Zlo)
               dval_dlogg = result_2D(2,Zlo)
               dval_dTeff = result_2D(3,Zlo)
            
            else ! Z interpolation

               fZ(1,1:numZs) = result_2D(1,Zlo:Zhi)
               call interp_pm(ai% logZ(Zlo:Zhi),numZs,fZ1,mp_work_size,work,'atm',ierr)
               if (ierr /= 0) then
                  if (gtv_dbg) write(*,*) 'interp_pm failed for Z interpolation'
                  return
               end if
               call interp_values(ai% logZ(Zlo:Zhi),numZs,fZ1,1,newlogZ,result_Z(1:1),ierr)
               if (ierr /= 0) then
                  if (gtv_dbg) write(*,*) 'interp_values failed for Z interpolation'
                  return
               end if
               newval = result_Z(1)

               fZ(1,1:numZs) = result_2D(2,Zlo:Zhi)
               call interp_pm(ai% logZ(Zlo:Zhi),numZs,fZ1,mp_work_size,work,'atm',ierr)
               if (ierr /= 0) then
                  if (gtv_dbg) write(*,*) 'interp_pm failed for Z interpolation of d_dlogg'
                  return
               end if
               call interp_values(ai% logZ(Zlo:Zhi),numZs,fZ1,1,newlogZ,result_Z(1:1),ierr)
               if (ierr /= 0) then
                  if (gtv_dbg) write(*,*) 'interp_values failed for Z interpolation of d_dlogg'
                  return
               end if
               dval_dlogg = result_Z(1)

               fZ(1,1:numZs) = result_2D(3,Zlo:Zhi)
               call interp_pm(ai% logZ(Zlo:Zhi),numZs,fZ1,mp_work_size,work,'atm',ierr)
               if (ierr /= 0) then
                  if (gtv_dbg) write(*,*) 'interp_pm failed for Z interpolation of d_dTeff'
                  return
               end if
               call interp_values(ai% logZ(Zlo:Zhi),numZs,fZ1,1,newlogZ,result_Z(1:1),ierr)
               if (ierr /= 0) then
                  if (gtv_dbg) write(*,*) 'interp_values failed for Z interpolation of d_dTeff'
                  return
               end if
               dval_dTeff = result_Z(1)
            
            end if
            
         end subroutine do_interp



         !reads in the iZth atm pressure table and pre-processes it for 2D spline interpolation
         !optional: reads or writes binary version of the spline table to atm_data/cache
         subroutine load_atm_table(iZ,ierr)
            use utils_lib
            use interp_2D_lib_db, only: interp_mkbicub_db
            use const_def, only: mesa_data_dir
            implicit none
            integer, intent(in)  :: iZ !index of Z table to be loaded
            integer, intent(out) :: ierr
            integer :: iounit, i, j, ibound_tmp(ng), ibcTmin, ibcTmax, ibcgmin, ibcgmax
            integer :: cache_file_version, text_file_version, nvec
            real(dp) :: Teff_tmp(nT), logg_tmp(ng)
            real(dp) :: bcTmin(nT), bcTmax(ng), bcgmin(ng), bcgmax(ng), data_tmp(ng,nT)
            real(dp), target :: vec_ary(100)
            real(dp), pointer :: f1(:), vec(:)
            character(len=2000) :: buf

            ierr = 0
            vec => vec_ary

            !read in and process the text files
            filename = trim(mesa_data_dir) // '/atm_data/' // trim(ai% table_atm_files(iZ))
            open(newunit=iounit,file=trim(filename),action='read',status='old',iostat=ierr)
            if (ierr/= 0) then
               write(*,*) 'load_atm_table: missing atm data:'
               write(*,*) trim(filename)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*) 'FATAL ERROR: missing or bad atm data.'
               write(*,*) 'Please update by removing the directory mesa/data/atm_data,'
               write(*,*) 'and rerunning the mesa ./install script.'
               write(*,*)
               call mesa_error(__FILE__,__LINE__)
            endif
            
            read(iounit,'(14x,i4)',iostat=ierr) text_file_version
            if (failed(1)) return
            if (text_file_version /= table_atm_version) then
               write(*,*) 'load_atm_table: mismatch in table versions'
               write(*,*) 'expected version ', table_atm_version
               write(*,*) 'received version ', text_file_version
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*) 'FATAL ERROR: out-of-date verion of atm data.'
               write(*,*) 'Please update by removing the directory mesa/data/atm_data,'
               write(*,*) 'and rerunning the mesa ./install script.'
               write(*,*)
               write(*,*)
               call mesa_error(__FILE__,__LINE__)
            endif
            
            ibound_tmp = -1
            read(iounit,1,iostat=ierr) ai% logZ(iZ), ai% alphaFe(iZ), ai% atm_mix(iZ), ibound_tmp(1:ng)
            if (ierr /= 0) then
               write(*,*) 'iZ', iZ
               write(*,*) 'ai% logZ(iZ)', ai% logZ(iZ)
               write(*,*) 'ai% alphaFe(iZ)', ai% alphaFe(iZ)
               write(*,*) 'ai% atm_mix(iZ)', ai% atm_mix(iZ)
               do j=1,ng
                  write(*,*) j, ibound_tmp(j)
               end do
               write(*,*)
            end if
            if (failed(2)) return
            read(iounit,2,iostat=ierr) logg_tmp(:)
            if (failed(3)) return
            data_tmp(:,:) = -1
            do j=1,nT
               read(iounit,'(a)',iostat=ierr) buf
               !write(*,*) 'read ierr', ierr
               if (ierr == 0) then
                  call str_to_vector(buf, vec, nvec, ierr)
                  !write(*,*) 'str_to_vector ierr', ierr
                  if (nvec < ng + 1) ierr = -1
                  !write(*,*) 'nvec ierr', ierr
               end if
               if (ierr /= 0) then
                  write(*,*) 'failed reading Pgas Teff', ai% Teff_array(j)
                  write(*,'(a)') 'buf "' // trim(buf) // '"'
                  write(*,*) 'len_trim(buf)', len_trim(buf)
                  write(*,*) 'nvec', nvec
                  write(*,*) 'ng', ng
                  stop
               end if
               if (failed(4)) return
               Teff_tmp(j) = vec(1)
               do i=1,ng
                  data_tmp(i,j) = vec(i+1)
               end do               
            enddo
            ai% Pgas_interp(1,:,:,iZ) = data_tmp(:,:)
            
            if (ai% id /= ATM_TABLE_PHOTOSPHERE) then ! read T
               read(iounit,2,iostat=ierr) ! skip line
               if (failed(5)) return
               do j=1,nT
                  read(iounit,'(a)',iostat=ierr) buf
                  if (ierr == 0) then
                     call str_to_vector(buf, vec, nvec, ierr)
                     if (nvec < ng + 1) ierr = -1
                  end if
                  if (ierr /= 0) then
                     write(*,*) 'failed reading T Teff', ai% Teff_array(j)
                  end if
                  if (failed(6)) return
                  Teff_tmp(j) = vec(1)
                  do i=1,ng
                     data_tmp(i,j) = vec(i+1)
                  end do 
               enddo
               ai% T_interp(1,:,:,iZ) = data_tmp(:,:)
            end if            

            close(iounit)

 1          format(13x,f5.2,8x,f4.1,1x,a8,1x,15x,99i4)
 2          format(15x,99(9x,f5.2,1x))
 3          format(99e15.7)

            !verify that Teff and logg arrays are the same as globals
            do i=1,nT
               if (abs(Teff_tmp(i) - ai% Teff_array(i)) > 1d-10*Teff_tmp(i)) then
                  ierr = -1
                  if (gtv_dbg) then
                     write(*,'(a30,i6,e24.12)') 'Teff_tmp(i)', i, Teff_tmp(i)
                     write(*,'(a30,i6,e24.12)') 'ai% Teff_array(i)', i, ai% Teff_array(i)
                  end if
                  return
               endif
            enddo

            do i=1,ng
               !write(*,*) 'logg', i, logg_tmp(i), ai% logg_array(i)
               if (abs(logg_tmp(i) - ai% logg_array(i)) > 1d-10*abs(logg_tmp(i)) ) then
                  ierr=-1
                  if (gtv_dbg) then
                     write(*,'(a30,i6,e24.12)') 'logg_tmp(i)', i, logg_tmp(i)
                     write(*,'(a30,i6,e24.12)') 'ai% logg_array(i)', i, ai% logg_array(i)
                  end if
                  return
               endif
            enddo
            !stop

            !make sure boundaries set earlier are still valid
            do i=1,ng
               ai% Teff_bound(i) = min( ai% Teff_bound(i), Teff_tmp(ibound_tmp(i)) )
            enddo
         
            ! use "not a knot" bc's
            ibcTmin = 0; bcTmin(:) = 0d0
            ibcTmax = 0; bcTmax(:) = 0d0
            ibcgmin = 0; bcgmin(:) = 0d0
            ibcgmax = 0; bcgmax(:) = 0d0

            !generate Teff,logg 2D-interpolants
            f1(1:4*ng*nT) => ai% Pgas_interp1(1+4*ng*nT*(iZ-1):4*ng*nT*iZ)
            call interp_mkbicub_db(ai% logg_array, ng, ai% Teff_array, nT, &
               f1, ng, ibcgmin, bcgmin, ibcgmax, bcgmax, &
               ibcTmin, bcTmin, ibcTmax, bcTmax, ai% iling, ai% ilinT, ierr )
            if (ierr /= 0) then
               if (gtv_dbg) write(*,*) 'interp_mkbicub_db failed for Pgas_interp'
               return
            end if
            
            if (ai% id /= ATM_TABLE_PHOTOSPHERE) then
               f1(1:4*ng*nT) => ai% T_interp1(1+4*ng*nT*(iZ-1):4*ng*nT*iZ)
               call interp_mkbicub_db(ai% logg_array, ng, ai% Teff_array, nT, &
                  f1, ng, ibcgmin, bcgmin, ibcgmax, bcgmax, &
                  ibcTmin, bcTmin, ibcTmax, bcTmax, ai% iling, ai% ilinT, ierr )
               if (ierr /= 0) then
                  if (gtv_dbg) write(*,*) 'interp_mkbicub_db failed for T_interp'
                  return
               end if
            end if

            !this file has been loaded and processed
            ai% have_atm_table(iZ) = .true.
         
         end subroutine load_atm_table
         
         
         logical function failed(i)
            integer, intent(in) :: i
            failed = (ierr /= 0)
            if (failed) then
               write(*,'(a,i6)') &
                  'failed in reading atm table ' // trim(filename), i
               stop 'get_table_values'
            end if
         end function failed
         

      end subroutine get_table_values


      end module table_atm
