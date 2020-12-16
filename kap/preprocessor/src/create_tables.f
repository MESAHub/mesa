! ***********************************************************************
!
!   Copyright (C) 2009-2019  Bill Paxton & The MESA Team
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

      program create_tables
      
      use create_fixed_metal_tables
      use create_co_enhanced_tables
      use kap_support
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none

      logical :: do_kappa_test = .false.
      !logical :: do_kappa_test = .true.
      
      double precision :: which_z, which_x
      integer :: io_unit, ios, ilgT, cnt, len, status
      character (len=64) :: inlist_fname, whichz_str, whichx_str
      
      status = 0
      call get_command_argument(1, inlist_fname, len, status)
      if (status /= 0) then
         write(*,*) 'failed to get inlist_fname from command line'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call get_command_argument(2, whichz_str, len, status)
      if (status /= 0) then
         write(*,*) 'failed to get whichz from command line'
         call mesa_error(__FILE__,__LINE__)
      end if
      read(whichz_str,*,iostat=status) which_z
      if (status /= 0) then
         write(*,*) 'failed to convert whichz from command line to number'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call get_command_argument(3, whichx_str, len, status)
      if (status /= 0) then
         write(*,*) 'failed to get whichx from command line'
         call mesa_error(__FILE__,__LINE__)
      end if
      read(whichx_str,*,iostat=status) which_x
      if (status /= 0) then
         write(*,*) 'failed to convert whichx from command line to number'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      !write(*,*) 'inlist_fname', trim(inlist_fname)
      !write(*,*) 'which_z', which_z
      !write(*,*) 'which_x', which_x

      call read_namelist(inlist_fname)
      
      call init_preprocessor
      call read_output_logTs

      if (do_kappa_test) then ! testing mode

         call Do_Test(data_dir, type1_table)
         call Do_CO_Test(data_dir, type1_table)
      
      else ! production mode
      
         which_z = min(1d0, max(0d0, which_z))
         
         if (CO_flag) then ! C/O enhanced
            call Write_CO_Files( &
               which_z, which_x, output_dir, data_dir, type1_table, &
               header_info, table_prefix, table_version)
         else
            call Write_Files( &
               which_z, which_x, output_dir, data_dir, type1_table, &
               header_info, table_prefix, table_version)
         end if
         
      end if
      
      
      
      end program create_tables

