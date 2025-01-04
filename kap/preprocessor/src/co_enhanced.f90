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

      module create_co_enhanced_tables

      use const_def
      use math_lib
      use kap_support
      use utils_lib, only: mesa_error

      implicit none

      contains
      
      
      subroutine Do_CO_Test(data_dir, type1_table)
      character (len=*), intent(in) :: data_dir, type1_table
      real(dp) :: Zbase, X, dXC, dXO, Rho, logRho, T, logT, logK, opac_rad
      integer :: info
      logical, parameter :: co_enhanced = .true.
      
 1    format(a40,1pe26.16)
      
      Zbase = 0.02d0
      X = 0.70d0
      dXC = 0
      dXO = 0
      
      call setup(co_enhanced, data_dir, type1_table, Zbase, X, dXC, dXO)
      
      logRho = -6d0
      rho = 10**logRho
      logT = 6d0
      T = 10**logT

      call Get_Results( &
      Zbase, X, dXC, dXO, Rho, logRho, T, logT,  &
      logK, co_enhanced, data_dir, type1_table, .false., info)
      
      write(*,*) 'co_enhanced', co_enhanced
      write(*,*) 'info', info
      write(*,1) 'logT', logT
      write(*,1) 'logRho', logRho
      write(*,1) 'z', z
      write(*,1) 'xh', x
      write(*,*)
      write(*,1) 'kap', 10**logK
      write(*,1) 'logK', logK
      write(*,*)
      
      end subroutine Do_CO_Test


      subroutine Write_CO_Files( &
            Zbase, which_x, output_dir, data_dir, type1_table, & 
            header_info, file_prefix, table_version)
         real(dp), intent(in) :: Zbase, which_x
         character (len=*), intent(in) :: &
            output_dir, data_dir, type1_table, header_info, file_prefix
         integer, intent(in) :: table_version
         integer, parameter :: max_num_Xs = 5
         real(dp) :: Xs(max_num_Xs)
         real(dp) :: X
         integer :: ix, num_Xs
         
         if (which_x < 0) then
            num_Xs = 5
            Xs(1:num_Xs) = (/ 0.00d0, 0.03d0, 0.10d0, 0.35d0, 0.70d0 /)
         else
            num_Xs = 1
            Xs(1) = which_x
         end if
         
         do ix = 1,num_Xs
            X = Xs(ix)
            call Do_CO_Tables(Zbase, X, output_dir, data_dir, type1_table, & 
                   header_info, file_prefix, table_version)
         end do
         
      end subroutine Write_CO_Files
      
      
      subroutine Create_Filename(Z, X, output_dir, file_prefix, fname)
         real(dp), intent(in) :: Z, X
         character (len=*),intent(in) :: output_dir, file_prefix
         character (len=*),intent(out) :: fname
         
         character (len=8) :: zstr, xstr
         integer :: i, ierr

         call get_output_string(X,xstr,ierr)
         call get_output_string(Z,zstr,ierr)
         if(ierr/=0) call mesa_error(__FILE__,__LINE__)
         
         write(fname,'(a)') &
            trim(output_dir) // &
            '/' // trim(file_prefix) // '_z' // &
            trim(zstr) // '_x' // trim(xstr) // '.data'
            
      end subroutine Create_Filename
      
      
      subroutine Do_CO_Tables( &
            Zbase, X, output_dir, data_dir, type1_table, header_info, file_prefix, &
            table_version)
         real(dp), intent(in) :: Zbase, X
         character (len=*),intent(in) :: output_dir, data_dir, type1_table, header_info, file_prefix
         integer, intent(in) :: table_version
         integer, parameter :: num_dXCOs = 8
         real(dp) :: dXCOs(num_dXCOs), Y, dXC, dXO, mid!, XC_base, XN_base, XO_base, XNe_base
         integer :: i, j, num_Tables, io_unit, pass, num1, num2, num3, ios
         character (len=256) :: fname
         real(dp), parameter :: tiny = 1d-10
         
         integer, parameter :: file_type = 2

         io_unit = 33
         call Create_Filename(Zbase, X, output_dir, file_prefix, fname)
         open(unit=io_unit, file=trim(fname), action='write', status='replace', iostat=ios)
         if (ios /= 0) then
            write(*,*) 'failed to open ', trim(fname)
            call mesa_error(__FILE__,__LINE__)
         end if      
         write(*,*) 'creating ', trim(fname)
         
         dXCOs(1:num_dXCOs) = (/ 0.00d0, 0.01d0, 0.03d0, 0.10d0, 0.20d0, 0.40d0, 0.60d0, 1.0d0 /)
         Y = 1 - (X + Zbase)    ! this sets upper limit on sum of dXC + dXO         
         mid = Y * 0.5d0
         
         do pass = 1, 3 
           ! count the tables on the 1st pass, write a list of tables on 2nd,
           ! and write the actual tables on the 3rd. 
            
            if (pass == 2) then ! write the header
               write(io_unit,'(a)') trim(header_info)
               write(io_unit,'(A8,2(2x,A8),6(2x,A16),99(2x,A8))') &
                  'form', 'version', 'Tables','X  ','Zbase', &
                  'Zfrac_C', 'Zfrac_N', 'Zfrac_O', 'Zfrac_Ne', &
                  'logRs','min','max','logTs','min','max'
               write(io_unit,advance='no',fmt='(I8,2(2x,I8),6(2x,F16.6))') &
                  file_type, table_version, num_Tables, X, Zbase, &
                  Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne
               write(io_unit,advance='no',fmt='(2x,I8,2(2x,F8.2))') logR_points, logR_min, logR_max
               write(io_unit,fmt='(2x,I8,2(2x,F8.2))') logT_points, logT_min, logT_max
            end if
            
            num_Tables = 0          
            
            if (pass == 2) write(io_unit,'(/,i3,a,/,4(a12))') num1, ' tables with dXC > dXO', 'Num',  &
            'Y   ', 'dXC ', 'dXO '
            ! 1) dXC > dXO
            do j = 1, num_dXCOs
               dXO = dXCOs(j)
               if (dXO >= mid-tiny) exit
               do i = j+1, num_dXCOs
                  dXC = dXCOs(i)
                  if (dXC+dXO > Y) dXC = Y-dXO
                  num_Tables = num_Tables + 1
                  if (pass == 2) then
                     write(io_unit,'(i12,3f12.4)') num_Tables, Y-(dXC+dXO), dXC, dXO
                  else if (pass == 3) then
                     call Do_CO_Table(Zbase,X,dXC,dXO,num_Tables,io_unit,data_dir,type1_table)
                  end if
                  if (dXC+dXO > Y-tiny) exit
               end do
            end do
            if (pass == 1) num1 = num_Tables
            
            if (pass == 2) write(io_unit,'(/,i3,a,/,4(a12))') num2, ' tables with dXC == dXO',  &
                  'Num', 'Y   ', 'dXC ', 'dXO '
            ! 2) dXC == dXO
            do i = 1, num_dXCOs
               if (dXCOs(i) > mid+tiny) exit
               num_Tables = num_Tables + 1
               dXC = dXCOs(i); dXO = dXC
               if (pass == 2) then
                  write(io_unit,'(i12,3f12.4)') num_Tables, Y-(dXC+dXO), dXC, dXO
               else if (pass == 3) then
                  call Do_CO_Table(Zbase, X, dXC, dXO, num_Tables, io_unit, data_dir, type1_table)
               end if
               if (dXCOs(i) >= mid-tiny) exit
            end do
            if (pass == 1) num2 = num_Tables - num1
            
            if (pass == 2) write(io_unit,'(/,i3,a,/,4(a12))') num3,  &
                     ' tables with dXC < dXO', 'Num', 'Y   ', 'dXC ', 'dXO '
            ! 3) dXC < dXO
            do j = 1, num_dXCOs
               dXC = dXCOs(j)
               if (dXC >= mid-tiny) exit
               do i = j+1, num_dXCOs
                  dXO = dXCOs(i)
                  if (dXC+dXO > Y) dXO = Y-dXC
                  num_Tables = num_Tables + 1
                  if (pass == 2) then
                     write(io_unit,'(i12,3f12.4)') num_Tables, Y-(dXC+dXO), dXC, dXO
                  else if (pass == 3) then
                     call Do_CO_Table(Zbase,X,dXC,dXO,num_Tables,io_unit,data_dir,type1_table)
                  end if
                  if (dXC+dXO > Y-tiny) exit
               end do
            end do
            if (pass == 1) num3 = num_Tables - (num1+num2)
            
            if (pass == 2) write(io_unit,'(/,a,/)')  &
            '---------------------------------------------------------------------------'
            
         end do
         
         close(io_unit)
         
      end subroutine Do_CO_Tables
      
      
      subroutine Do_CO_Table(Zbase, X, XC, XO, table_num, io_unit, data_dir, type1_table)
         real(dp), intent(in) :: Zbase, X, XC, XO
         integer, intent(in) :: table_num, io_unit
         character (len=*), intent(in) :: data_dir, type1_table
         
         integer i,j,info
         real(dp) :: opac_rad, logR, logT, T, logRho, Rho, logK
         logical, parameter :: co_enhanced = .true.
         
         call setup(co_enhanced, data_dir, type1_table, Zbase, X, XC, XO)
      
         write(io_unit,'(A10,10x,10(A11))') 'Table', 'X  ',  'Y  ', 'Zbase', 'dXC ', 'dXO '
         write(io_unit,'(6x,I2,17x,5(F6.3,5x))') table_num, X, Y, Zbase, XC, XO
         write(io_unit,'(/,a)') '   logT                       logR = logRho - 3*logT + 18'

         write(io_unit, advance='NO', fmt='(8x)')
         do j=1,logR_points
            logR = logR_min + dlogR*(j-1)
            write(io_unit, advance='NO', fmt='(F8.3)') logR
         enddo
         write(io_unit, *)
         write(io_unit, *)
         
         do i=1, logT_points
            logT = output_logTs(i)
            T = 10 ** logT
            write(io_unit, advance='NO', fmt='(F8.3)') logT
            do j=1,logR_points
               logR = logR_min + dlogR*(j-1)
               logRho = logR + 3*logT - 18
               Rho = 10**logRho
               info = 0
               call Get_Results(Zbase, X, XC, XO, Rho, logRho, T, logT,  &
                  logK, co_enhanced, data_dir, type1_table, .false., info)
               if (info /= 0 .or. logK > 99d0 .or. logK < -99d0 .or. logK-logK /= 0d0) logK = -99.999d0
               
               write(io_unit, advance='NO', fmt='(F8.3)') logK
            enddo
            write(io_unit,*)
         enddo
         write(io_unit,*)
            
      end subroutine Do_CO_Table


      end module create_co_enhanced_tables
