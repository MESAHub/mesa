! ***********************************************************************
!
!   Copyright (C) 2009-2019  Bill Paxton & The MESA Team

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

      module create_fixed_metal_tables
      use const_def
      use kap_support
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none

      contains
      
      
      subroutine Do_Test(data_dir, type1_table)
         character (len=*), intent(in) :: data_dir, type1_table
         real(dp) :: Zbase, X, XC, XO, Rho, logRho, T, logT, logK, logR
         integer :: info
         logical, parameter :: co_enhanced = .false.
         
         1 format(a40,1p e26.16)
         
         Zbase = 0.01d0
         X = 0.70d0
         XC = 0d0
         XO = 0d0
   
         call setup( &
            co_enhanced, data_dir, type1_table, Zbase, X, XC, XO)
         
         ! NOTE: you must have a logR that matches the OPAL tables
         logR = -6d0
         logT = 6d0
         logRho = logR + 3*logT - 18
         T = 10**logT
         rho = 10**logRho
         
         write(*,1) 'logR', logR
         write(*,1) 'logT', logT
         write(*,1) 'logRho', logRho
         write(*,1) 'z', z
         write(*,1) 'xh', x
         write(*,*) 'call Get_Results'
         write(*,*)

         call Get_Results( &
            Zbase, X, XC, XO, Rho, logRho, T, logT,  &
            logK, co_enhanced, data_dir, type1_table, .false., info)
            
         write(*,*) 'info', info
         write(*,1) 'logK', logK
         write(*,1) 'kap', 10**logK
         write(*,*)
         stop 'Do_Test'
         
      end subroutine Do_Test


      subroutine Write_Files( &
            Z_in, which_x, output_dir, data_dir, type1_table, header_info, &
            file_prefix, table_version)
         real(dp), intent(in) :: Z_in, which_x
         character (len=*), intent(in) :: &
            output_dir, data_dir, type1_table, header_info, file_prefix
         integer, intent(in) :: table_version
         integer, parameter :: max_num_Xs = 10
         real(dp) :: Z, Xs(max_num_Xs)
         integer :: ix, iz, num_Xs

         Z = Z_in
         if (which_x < 0) then
            Xs(1:max_num_Xs-1) = &
               (/ 0.00d0, 0.10d0, 0.20d0, 0.35d0, 0.50d0, 0.70d0, 0.80d0, 0.90d0, 0.95d0 /)
               num_Xs = max_num_Xs
         else
            Xs(1) = which_x
            num_Xs = 1
         end if

         if (num_Xs > 1) then
            iz = floor(Z*1d4 + 0.1d0)
            select case (iz)
            case (0)              ! Z = 0.0
               Xs(num_Xs) = 1
            case (1)              ! Z = 0.0001
               Xs(num_Xs) = 0.9999d0
            case (3)              ! Z = 0.0003
               Xs(num_Xs) = 0.9997d0
            case (10)             ! Z = 0.001
               Xs(num_Xs) = 0.999d0
            case (20)             ! Z = 0.002
               Xs(num_Xs) = 0.998d0
            case (40)             ! Z = 0.004
               Xs(num_Xs) = 0.996d0
            case (100)            ! Z = 0.01
               Xs(num_Xs) =  0.99d0
            case (200)            ! Z = 0.02
               Xs(num_Xs) =  0.98d0
            case (300)            ! Z = 0.03
               Xs(num_Xs) =  0.97d0
            case (400)            ! Z = 0.04
               Xs(num_Xs) =  0.96d0
            case (600)            ! Z = 0.06 
               num_Xs = max_num_Xs-1
               Xs(num_Xs) = 0.94d0
            case (800)            ! Z = 0.08
               num_Xs = max_num_Xs-1
               Xs(num_Xs) = 0.92d0
            case (1000)           ! Z = 0.1
               num_Xs = max_num_Xs-2
               Xs(num_Xs) = 0.90d0
            case default
               if (Freedman_flag) then ! use Ferguson Z=0.1 to fill gap for higher T
                  num_Xs = max_num_Xs-2
                  Xs(num_Xs) = 0.90d0
               else
                  write(*,*) 'unknown Z value for fixed metal', Z
                  stop 'Write_Files'
               end if
            end select
         end if

         do ix = 1,num_Xs
            call Do_Table( &
               Z, Xs(ix), output_dir, data_dir, type1_table, header_info, file_prefix, &
               table_version)
         end do
          
      end subroutine Write_Files
      
      
      subroutine Do_Table( &
         Z_in, X, output_dir, data_dir, type1_table, header_info, file_prefix, &
         table_version)
      real(dp), intent(in) :: Z_in, X
      character (len=*), intent(in) :: output_dir, data_dir, type1_table
      character (len=*), intent(in) :: header_info, file_prefix
      integer, intent(in) :: table_version
      integer, parameter :: file_type = 1
      real(dp), parameter :: XC = 0, XO = 0
      logical, parameter :: co_enhanced = .false.
      
      real(dp) :: logR, logT, T, logRho, Rho, logK
      integer :: i, j, io_unit, ios, info
      character (len=256) :: fname
      
      include 'formats'

      Z = Z_in

      io_unit = 34
         
      call create_fname(Z, X, output_dir, file_prefix, fname)
      
      open(unit=io_unit, file=trim(fname), action='write', status='replace', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'fixed_metal Do_Table failed to open ', trim(fname)
         call mesa_error(__FILE__,__LINE__)
      end if
   
      call setup( &
         co_enhanced, data_dir, type1_table, Z, X, XC, XO)

      ! header
      write(io_unit,'(a)') trim(header_info)
      write(io_unit,'(A8,99(2x,A10))') 'form','version', 'X   ','Z   ', 'logRs', 'logR_min', 'logR_max',&
         'logTs', 'logT_min', 'logT_max'
      write(io_unit,advance='no',fmt='(i8,2x,i10,2(2x,F10.6))') file_type, table_version, X, Z
      write(io_unit,advance='no',fmt='(2x,I10,2(2x,F10.6))') logR_points, logR_min, logR_max
      write(io_unit,fmt='(2x,I10,2(2x,F10.6))') logT_points, logT_min, logT_max
      
      ! data
      write(io_unit,'(/,a)') '   logT                       logR = logRho - 3*logT + 18'

      write(io_unit, advance='NO', fmt='(8x)')
      do j=1,logR_points
         logR = logR_min + dlogR*(j-1)
         write(io_unit, advance='NO', fmt='(F8.3)') logR
      enddo
      write(io_unit, *)
      write(io_unit, *)
      
      write(*,*) trim(fname)
         
      do i=1, logT_points
         logT = output_logTs(i)
         T = 10 ** logT
         write(io_unit, advance='NO', fmt='(F8.3)') logT
         do j=1,logR_points
            logR = logR_min + dlogR*(j-1)
            logRho = logR + 3*logT - 18
            Rho = 10**logRho
            info = 0
            call Get_Results(Z, X, XC, XO, Rho, logRho, T, logT,  &
               logK, co_enhanced, data_dir, type1_table, .false., info)
            if (info /= 0 .or. logK > 99d0 .or. logK < -99d0 .or. logK-logK /= 0d0) then
               logK = -99.999d0
            end if
            write(io_unit, advance='NO', fmt='(F8.3)') logK
         enddo
         write(io_unit,*)
      enddo
      write(io_unit,*)
      
      close(io_unit)
      
      end subroutine Do_Table
      
      
      subroutine create_fname(Z, X, output_dir, file_prefix, fname)
         real(dp), intent(in) :: Z, X
         character (len=*), intent(in) :: output_dir, file_prefix
         character (len=*), intent(out) :: fname
         character (len=8) :: zstr, xstr
         integer :: ierr
         
         if (Freedman_flag) then
            call get_output_string(Z,zstr,ierr)
            if (ierr/=0) call mesa_error(__FILE__,__LINE__)
            fname = trim(output_dir) // '/' // trim(file_prefix)// '_z' // &
               trim(zstr) // '.data'
         else     
            call get_output_string(Z,zstr,ierr)
            call get_output_string(X,xstr,ierr)
            if (ierr/=0) call mesa_error(__FILE__,__LINE__)
            fname = trim(output_dir) // '/' // trim(file_prefix)// '_z' // &
               trim(zstr) // '_x' // trim(xstr) // '.data'
         end if
            
      end subroutine create_fname
      
      
      end module create_fixed_metal_tables  

