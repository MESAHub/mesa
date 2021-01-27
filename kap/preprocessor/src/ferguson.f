! ***********************************************************************
!
!   Copyright (C) 2009-2018  Bill Paxton & The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module ferg
      use interp_2D_lib_sg
      use math_lib
      use utils_lib, only: mesa_error

      implicit none
      
        real, pointer :: ferg_logTs(:), ferg_logRs(:)
        real, pointer :: ferg_logKs(:,:,:), ferg_logKs1(:)
      
      real :: ferg_logT_min, ferg_logT_max, ferg_logR_min, ferg_logR_max, ferg_dlogR
      integer :: ferg_num_logTs, ferg_num_logRs
      real :: logT_for_smooth_lowT
      integer :: num_smooth_lowT
      
      !for AF94 tables
      logical :: AF94_lowT
                
      contains

      
      
      subroutine init_ferg(Z, X, data_dir, prefix, Freedman_flag)
         double precision, intent(in) :: Z, X
         character (len=*), intent(in) :: data_dir, prefix
         logical, intent(in) :: Freedman_flag
         character (len=256) :: fname
         
         if (AF94_lowT) then
            ferg_logT_min = 3.0
            ferg_logT_max = 4.1
            ferg_num_logTs = 23
            ferg_logR_min = -7.00
            ferg_logR_max = 1.00
            ferg_dlogR = 0.5
            ferg_num_logRs = 17
         else
            ferg_logT_min = 2.70
            ferg_logT_max = 4.50
            ferg_num_logTs = 85
            ferg_logR_min = -8.00
            ferg_logR_max = 1.00
            ferg_dlogR = 0.5
            ferg_num_logRs = 19
         end if
         allocate( &
            ferg_logTs(ferg_num_logTs), &
            ferg_logRs(ferg_num_logRs), &
            ferg_logKs1(4*ferg_num_logRs*ferg_num_logTs))
         call set_ferg_logTs
         call set_ferg_logRs
         call get_ferg_filename(Z, X, data_dir, prefix, Freedman_flag, fname)
         call read_ferg(fname)
      end subroutine init_ferg
      
      
      subroutine set_ferg_logTs
         integer :: ilgT, cnt
         cnt = 0
         
         if (AF94_lowT) then
            do ilgT = 300, 410, 5
               cnt = cnt+1
               ferg_logTs(cnt) = float(ilgT)/100.0
            end do
            if (cnt /= ferg_num_logTs) then
               write(*,*) 'cnt /= ferg_num_logTs_AF94) in set_ferg_logTs'
               call mesa_error(__FILE__,__LINE__)
            end if
            return
         end if
         
         do ilgT = 270, 285, 5
            cnt = cnt+1
            ferg_logTs(cnt) = float(ilgT)/100.0
         end do
         do ilgT = 290, 349
            cnt = cnt+1
            ferg_logTs(cnt) = float(ilgT)/100.0
         end do
         do ilgT = 350, 450, 5
            cnt = cnt+1
            ferg_logTs(cnt) = float(ilgT)/100.0
         end do
         if (cnt /= ferg_num_logTs) then
            write(*,*) '(cnt /= ferg_num_logTs) in set_ferg_logTs'
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine set_ferg_logTs !'
      
      
      subroutine set_ferg_logRs
         integer :: i
         if (AF94_lowT) then
            do i=1, ferg_num_logRs
               ferg_logRs(i) = ferg_logR_min + (i-1)*ferg_dlogR
            end do
            return
         end if
         do i = 1, ferg_num_logRs
            ferg_logRs(i) = ferg_logR_min + (i-1)*ferg_dlogR
         end do
      end subroutine set_ferg_logRs
      
      
      subroutine read_ferg(fname)
         character (len=*), intent(in) :: fname
      
         integer :: io_logK, i, j, hdr, num_logRs, num_logTs, n1
         real :: logT, lKs(ferg_num_logRs), lRs(ferg_num_logRs), temp_logKs(4, ferg_num_logRs, ferg_num_logTs)
         integer :: ierr                       ! =0 on exit if there is no error.
         character (len=6) :: str
      
         ierr = 0
         
         num_logRs = ferg_num_logRs
         num_logTs = ferg_num_logTs

         ferg_logKs(1:4,1:ferg_num_logRs,1:ferg_num_logTs) => &
             ferg_logKs1(1:4*ferg_num_logRs*ferg_num_logTs)

         io_logK = 40
         open(unit=io_logK, file=trim(fname), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'read_ferg failed to open ', trim(fname)
            call mesa_error(__FILE__,__LINE__)
         end if

         if (AF94_lowT) then
            hdr = 4
         else
            hdr = 3
         end if
   
         do i = 1, hdr ! skip header lines
            read(io_logK, *, iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'read_ferg failed while reading ', trim(fname)
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         
         read(io_logK, '(a6,99(f7.3))', iostat=ierr) str, lRs(1:num_logRs)
         if (ierr /= 0) then
            write(*,*) 'read_ferg failed while reading logRs ', trim(fname)
            call mesa_error(__FILE__,__LINE__)
         end if
         do i = 1, num_logRs
            if (abs(lRs(i) - ferg_logRs(i)) > 0.00001) then
               write(*,*) 'read_ferg: unexpected value for logR', i, lRs(i), ferg_logRs(i)
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         
         if (AF94_lowT) then
            read(io_logK, *, iostat=ierr)
            read(io_logK, *, iostat=ierr)
            if (ierr/=0) write(*,*) 'read_ferg failed while reading ', trim(fname)
         end if
                  
         do i = num_logTs, 1, -1
            if (AF94_lowT) then
               read(io_logK, '(f4.2,99f7.3)', iostat=ierr) logT, lKs(1:num_logRs)            
            else
               read(io_logK, '(f5.3,1x,99f7.3)', iostat=ierr) logT, lKs(1:num_logRs)
            endif
            if (ierr /= 0) then
               write(*,*) 'read_ferg failed while reading logKs ', trim(fname)
               write(*,*) 'i', i
               write(*,*) 'num_logRs', num_logRs
               write(*,*) 'num_logTs', num_logTs
               write(*,*) 'logT', logT             
               call mesa_error(__FILE__,__LINE__)
            end if
            ferg_logKs(1, 1:num_logRs, i) = lKs(1:num_logRs)
            if (abs(logT - ferg_logTs(i)) > 0.00001) then
               write(*,*) 'read_ferg: unexpected value for logT', i, logT, ferg_logTs(i)
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         
         close(io_logK)
         
         ! smooth logKs for logT < logT_for_smooth_lowT
         do j = 1, num_smooth_lowT
            temp_logKs = ferg_logKs
            n1 = max(1, num_logRs-j/2)
            do i = 1, num_logTs-1
               if (ferg_logTs(i) > logT_for_smooth_lowT) exit
               !write(*,*) 'smooth', i, ferg_logTs(i)
               if (i == 1) then
                  ferg_logKs(1, n1:num_logRs, i) = ( &
                     temp_logKs(1, n1:num_logRs, i) + &
                     temp_logKs(1, n1:num_logRs, i+1))/2
               else
                  ferg_logKs(1, n1:num_logRs, i) = ( &
                     temp_logKs(1, n1:num_logRs, i-1) + &
                     temp_logKs(1, n1:num_logRs, i) + &
                     temp_logKs(1, n1:num_logRs, i+1))/3
               end if
            end do
         end do
         
         call interp_mkbipm_sg( &
            ferg_logRs, num_logRs, ferg_logTs, num_logTs, &
            ferg_logKs1, num_logRs, ierr)
         if (ierr /= 0) then
            write(*,*) 'interp_mkbipm_sg error happened for logKs in read_ferg'
            call mesa_error(__FILE__,__LINE__)
         end if
         
      end subroutine read_ferg

      
      subroutine get_ferg_filename(Z, X, data_dir, prefix, Freedman_flag, fname)
         double precision, intent(in) :: Z, X
         character (len=*),intent(in) :: data_dir, prefix
         logical, intent(in) :: Freedman_flag
         character (len=*),intent(out) :: fname       
         character (len=16) :: zstr, xstr       
         integer :: iz, ix
         if (Freedman_flag) then
            iz = 10000
            ix = 70000
         else
            iz = floor(Z*1d5 + 0.1d0)  
            ix = floor(X*1d5 + 0.1d0)
         end if
         select case (iz)
            case (0)
               zstr = '0'
            case (1)
               zstr = '00001'
            case (3)
               zstr = '00003'
            case (10)
               zstr = '0001'
            case (30)
               zstr = '0003'
            case (100)
               zstr = '001'
            case (200)
               zstr = '002'
            case (400)
               zstr = '004'
            case (1000)
               zstr = '01'
            case (2000)
               zstr = '02'
            case (3000)
               zstr = '03'
            case (4000)
               zstr = '04'
            case (5000)
               zstr = '05'
            case (6000)
               zstr = '06'
            case (8000)
               zstr = '08'
            case (10000)
               if (AF94_lowT) then
                  zstr = '10'
               else
                  zstr = '1'
               endif
            case default
               write(*,*) 'get_Zstr: unknown Z value for ferg data', Z
               call mesa_error(__FILE__,__LINE__)
         end select

         select case (ix)
            case (0)
               xstr = '0'
            case (10000)
               xstr = '1'
            case (20000)
               xstr = '2'
            case (35000)
               xstr = '35'
            case (50000)
               xstr = '5'
            case (70000)
               xstr = '7'
            case (80000)
               xstr = '8'
            case (90000)
               xstr = '9'
            case (92000)
               xstr = '92'
            case (94000)
               xstr = '94'
            case (95000)
               xstr = '95'
            case (96000)
               xstr = '96'
            case (97000)
               xstr = '97'
            case (98000)
               xstr = '98'
            case (99000)
               xstr = '99'
            case (99600)
               xstr = '996'
            case (99800)
               xstr = '998'
            case (99900)
               xstr = '999'
            case (99970)
               xstr = '9997'
            case (99990)
               xstr = '9999'
            case (99997)
               xstr = '99997'
            case (99999)
               xstr = '99999'
            case (100000)
               xstr = '99999'
               zstr = '00001'
            case default
               write(*,*) 'get_Xstr: unknown X value for ferg data', X, ix
               call mesa_error(__FILE__,__LINE__)
         end select
      
         !ah, the joys of dealing with opacity tables...
         !AF94 does not have Z=0.06, 0.08 for X>=0.9, it has Z=0.07 instead
         if (AF94_lowT .and. ix >= 90000 .and. (iz == 6000 .or. iz== 8000)) zstr = '07'
         !same rant; AF94 calls these gz*
         if (AF94_lowT .and. ix > 90000) xstr = 'z'

         if ( AF94_lowT ) then
            fname = trim(data_dir) // '/' // trim(prefix) // trim(xstr) &
               // '.' // trim(zstr) // '.tron'
         else
            fname = trim(data_dir) // '/' // trim(prefix) // '.' &
               // trim(xstr) // '.' // trim(zstr) // '.tron'
         end if
         
         !write(*,*) 'read ' // trim(fname)
                        
      end subroutine get_ferg_filename

      
      subroutine eval_ferg(z_in, xh_in, t6_in, r_in, logKap, dbg)
         double precision, intent(in) :: z_in, xh_in, t6_in, r_in
         double precision, intent(out) :: logKap
         logical, intent(in) :: dbg
         
         real :: z, xh, t6, r, logT, logR, f
         integer :: ier, i, j
         real :: min_ferg_logT_for_logR_gt_1
         
         integer :: num_logTs, num_logRs
         real :: logR_min, logR_max, logT_min, logT_max, tiny

         include 'formats'

         num_logTs = ferg_num_logTs
         num_logRs = ferg_num_logRs
         logR_min = ferg_logR_min
         logR_max = ferg_logR_max
         logT_min = ferg_logT_min
         logT_max = ferg_logT_max
         if (AF94_lowT) then
            min_ferg_logT_for_logR_gt_1 = 3.0
         else
            min_ferg_logT_for_logR_gt_1 = 3.5
         endif
      
         z = real(z_in); xh = real(xh_in); t6 = real(t6_in); r = real(r_in)

         logT = min(ferg_logTs(num_logTs), max(ferg_logTs(1), log10(t6*1e6)))
         logR = log10(R)
         logR = max(logR_min, logR)
         if (logR > logR_max) then
            logR = logR_max
            logT = max(logT, min_ferg_logT_for_logR_gt_1)
         end if
         tiny = 1e-6
         logR=min(max(logR, ferg_logRs(1)+tiny), ferg_logRs(num_logRs)-tiny)
         logT=min(max(logT, ferg_logTs(1)+tiny), ferg_logTs(num_logTs)-tiny)
         
         ier = 0
         call interp_evbipm_sg( &
            logR, logT, ferg_logRs, num_logRs, ferg_logTs, num_logTs, &
            ferg_logKs1, num_logRs, f, ier)
            
         if (ier /= 0) then
            write(*,*) 'interp_evbicub_sg error happened in eval_ferg'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         logKap = dble(f)
         
         if (dbg) then
            write(*,*)
            write(*,*) 'eval_ferg'
            write(*,1) 'z_in', z_in
            write(*,1) 'xh_in', xh_in
            write(*,1) 't6_in', t6_in
            write(*,1) 'r_in', r_in
            write(*,*)
            write(*,1) 'logT', logT
            write(*,1) 'logR', logR
            write(*,1) 'logRho', logR + 3*logT - 18
            write(*,*)
            write(*,1) 'logKap', logKap
            write(*,*)
            i = 15  ! logR index
            j = 73  ! logT index
            write(*,2) 'ferg_logRs(i)', i, ferg_logRs(i)
            write(*,2) 'ferg_logTs(j)', i, ferg_logTs(j)
            write(*,1) 'ferg_logKs(1,i,j)', ferg_logKs(1,i,j)
            write(*,*)
            call mesa_error(__FILE__,__LINE__)
         end if

      end subroutine eval_ferg


      end module ferg
