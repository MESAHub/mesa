!     ***********************************************************************
!     
!     Copyright (C) 2008-2019  Bill Paxton, Frank Timmes & The MESA Team
!     
!     This file is part of MESA.
!     
!     MESA is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Library Public License as published
!     by the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!     
!     MESA is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU Library General Public License for more details.
!     
!     You should have received a copy of the GNU Library General Public License
!     along with this software; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!     
!     ***********************************************************************

      module kap_support
      
      use const_lib
      use const_def
      use math_lib
      use ferg
      use utils_lib, only: mesa_error
      
      implicit none

      logical :: setup_done = .false.
      

      character (len=256) :: &
      data_dir, type1_table, table_prefix, header_info, &
      logT_file, output_dir, freedman_data_dir
      integer :: table_version
      logical :: OP_file, CO_flag, lowT_flag, Freedman_flag
      real(dp) :: ferg_logT_for_smooth_lowT ! doesn't apply to Freedman data.  ferg only.
      integer :: ferg_num_smooth_lowT
      
      real(dp), pointer :: output_logTs(:) ! (logT_points)

      real(dp) :: logR_min
      real(dp) :: logR_max
      integer :: logR_points

      real(dp) :: Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne
      
      namelist / kappa / &
      data_dir, & 
      type1_table, logT_file, output_dir, freedman_data_dir, &
      table_prefix, header_info, table_version, &
      OP_file, CO_flag, lowT_flag, Freedman_flag, AF94_lowT, &
      logR_min, logR_max, logR_points, &
      Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne, ferg_logT_for_smooth_lowT
      
      integer :: eos_handle
      real(dp) :: logT_min
      real(dp) :: logT_max
      integer :: logT_points
      real(dp) :: dlogR
      
      real(dp) :: T_min, T_max, R_min, R_max
      real(dp) :: Y, Z
      

      contains
      
      
      subroutine init_preprocessor
      use chem_lib
      use chem_def
      use const_lib
      use eos_lib
      
      logical :: use_cache
      integer :: ierr, i
      
      nullify(output_logTs)
      
      call const_init('../..', ierr)
      if (ierr /= 0) then
         write(*,*) 'const_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if

      call math_init()

      call chem_init('isotopes.data', ierr)
      if (ierr /= 0) then
         write(*,*) 'chem_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      end subroutine init_preprocessor
      
      
      subroutine read_namelist(inlist_fname)
         use ferg, only: logT_for_smooth_lowT, num_smooth_lowT
         character (len=*) :: inlist_fname
         integer :: io_unit, ios
         include 'formats'
         
         !set default values of namelist variables
         output_dir = ''
         data_dir = ''     
         type1_table = ''
         freedman_data_dir = ''
         header_info = ''
         lowT_flag = .false.
         logT_file = ''
         Freedman_flag = .false.
         OP_file = .false.
         CO_flag = .false.
         AF94_lowT = .false.
         table_prefix = ''
         table_version = -1
         Zfrac_C = -1
         Zfrac_N = -1
         Zfrac_O = -1
         Zfrac_Ne = -1
         ! build tables at twice the OPAL resolution to reduce spline problems.
         logR_min = 0
         logR_max = 0
         logR_points = 0
         ferg_logT_for_smooth_lowT = 3.3d0
         ferg_num_smooth_lowT = 10

         !read namelist file
         io_unit = 40
         open(io_unit,file=trim(inlist_fname),action='read',delim='quote', &
            status='old',iostat=ios)
         if(ios/=0)then
            write(*,*)
            write(*,*)
            write(*,*) 'Failed to open namelist file: ' // trim(inlist_fname)
            write(*,*)
            write(*,*)
            call mesa_error(__FILE__,__LINE__)
         endif
         read(io_unit, nml=kappa, iostat=ios)
         if(ios/=0)then
            write(*,*)
            write(*,*)
            write(*,*) 'Failed while trying to read namelist file: ' // trim(inlist_fname)
            write(*,*)
            write(*,*)
            call mesa_error(__FILE__,__LINE__)
         endif
         
         R_min = 10**logR_min
         R_max = 10**logR_max
         dlogR = (logR_max - logR_min)/(logR_points - 1)
         
         logT_for_smooth_lowT = real(ferg_logT_for_smooth_lowT)
         num_smooth_lowT = ferg_num_smooth_lowT

         return
         write(*,2) 'logR_points', logR_points
         write(*,1) 'logR_max', logR_max
         write(*,1) 'logR_min', logR_min
         write(*,1) 'dlogR', dlogR
         
         
      end subroutine read_namelist
      
      
      subroutine setup( &
            co_enhanced, data_dir, type1_table, Zbase, X, dXC, dXO)
         use chem_def
         use freedman, only: init_freedman
         logical, intent(in) :: co_enhanced
         character (len=*), intent(in) :: data_dir, type1_table
         real(dp), intent(in) :: Zbase, X, dXC, dXO
         real(dp) :: zalf, ztemp, w(6)
         integer i,iz
         
         include 'formats'
         
         Z = Zbase + dXC + dXO
         if (Z + X > 1d0+1d-9) then
            write(*,*) 'bad composition in setup for kap builder'
            write(*,1) 'Z', Z
            write(*,1) 'X', X
            write(*,1) 'Zbase', Zbase
            write(*,1) 'dXC', dXC
            write(*,1) 'dXO', dXO
            write(*,1) 'Zbase+X+dXC+dXO', Zbase+X+dXC+dXO
            write(*,*)
            call mesa_error(__FILE__,__LINE__)
         end if
         
         Y = max(0d0, min(1d0, 1 - (X + Z)))
         
         if (co_enhanced) then
            call init_opal_type2(Zbase, X)
         else if (.not. lowT_flag) then
            call init_opal_type1(Zbase, X)
         else
            call init_ferg(Zbase, X, data_dir, type1_table, Freedman_flag)
            if (Freedman_flag) call init_freedman(freedman_data_dir, Z)
         end if
      
      end subroutine setup
      
      
      subroutine read_output_logTs
         integer :: iounit, ierr, cnt
         
         include 'formats'
         
         iounit = 33
         open(unit=iounit, file=trim(logT_file), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ', trim(logT_file)
            call mesa_error(__FILE__,__LINE__)
         end if
         read(iounit,*,iostat=ierr) logT_points
         if (ierr /= 0) then
            write(*,*) 'failed to read ', trim(logT_file)
            call mesa_error(__FILE__,__LINE__)
         end if
         if (associated(output_logTs)) deallocate(output_logTs)
         allocate(output_logTs(logT_points))
         do cnt=1,logT_points
            read(iounit,*,iostat=ierr) output_logTs(cnt)
            if (ierr /= 0) then
               write(*,*) 'failed to read ', trim(logT_file), cnt
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         close(iounit)
         logT_min = output_logTs(1)
         logT_max = output_logTs(logT_points)
         T_min = 10**logT_min
         T_max = 10**logT_max
         
         return
         write(*,2) 'logT_points', logT_points
         write(*,1) 'logT_min', logT_min
         write(*,1) 'logT_max', logT_max
         write(*,1) 'T_min', T_min
         write(*,1) 'T_max', T_max
         write(*,*)
         
      end subroutine read_output_logTs

      
      subroutine Get_Results( &
            Zbase, X, dXC_in, dXO_in, Rho_in, logRho_in, T_in, logT_in,  &
            logK, co_enhanced, data_dir, type1_table, dbg, ierr)
         use utils_lib
         real(dp), intent(in) :: Zbase ! initial metals mass fraction
         real(dp), intent(in) :: X ! hydrogen mass fraction
         real(dp), intent(in) :: dXC_in ! enhancement of carbon mass fraction
         real(dp), intent(in) :: dXO_in ! enhancement of oxygen mass mass fraction
         
         real(dp), intent(in) :: Rho_in, logRho_in ! the density
         real(dp), intent(in) :: T_in, logT_in ! the temperature
         logical, intent(in) :: co_enhanced
         character (len=*), intent(in) :: data_dir, type1_table
         logical, intent(in) :: dbg
         
         ! OUTPUT
         
         real(dp), intent(out) :: logK
         integer, intent(out) :: ierr
         
         ! locals
         real(dp) :: dXC, dXO, T6, R, frac, opac, Y, logT, T, logRho, Rho, logR
         character (len=256) :: type1_filename
         
         include 'formats'
         
         ierr = 0
         
         ! clip to table
         
         logT = min(logT_max, max(logT_min, logT_in))
         logRho = logRho_in
         Rho = Rho_in
         
         logR = logRho - 3*logT + 18
         if (logR > logR_max) then
            logR = logR_max
            logRho = logR + 3*logT - 18
            Rho = 10**logRho
         else if (logR < logR_min) then
            logR = logR_min
            logRho = logR + 3*logT - 18
            Rho = 10**logRho
         end if

         T = min(T_max, max(T_min, T_in))

         type1_filename = trim(data_dir) // '/' // trim(type1_table)
                  
         dXC = dXC_in; dXO = dXO_in
         Y = 1 - (X+Zbase+dXC+dXO)
         if (abs(Y) < 1d-6) Y = 0
         if (Y < 0) then ! adjust dXC and dXO
            if (dXC + dXO <= 0) then
               write(*,*) 'Get_Results: dXC + dXO <= 0'
               write(*,1) 'Y', Y
               write(*,1) 'dXC', dXC
               write(*,1) 'dXO', dXO
               write(*,1) 'X', X
               write(*,1) 'Zbase', Zbase
               write(*,1) '(X+Zbase+dXC+dXO)', (X+Zbase+dXC+dXO)
               call mesa_error(__FILE__,__LINE__)
            end if
            frac = (1 - X - Zbase) / (dXC + dXO)
            if (frac <= 0) then
               write(*,*) 'Get_Results: (frac <= 0)'
               write(*,1) 'frac', frac
               write(*,1) 'dXC', dXC
               write(*,1) 'dXO', dXO
               write(*,1) 'X', X
               write(*,1) 'Zbase', Zbase
               write(*,1) '1 - X - Zbase', 1 - X - Zbase
               call mesa_error(__FILE__,__LINE__)
            end if
            dXC = frac*dXC; dXO = frac*dXO
            Y = 1 - (X+Zbase+dXC+dXO)
            if (abs(Y) < 1d-6) Y = 0
            if (Y < 0) then
               write(*,*) 'frac', frac
               write(*,*) 'dXC', dXC
               write(*,*) 'dXO', dXO
               write(*,*) 'Zbase', Zbase
               write(*,*) 'X', X
               write(*,*) '(X+Zbase+dXC+dXO)', (X+Zbase+dXC+dXO)
               write(*,*) 'Y', Y
               write(*,*) 'Get_Results'
               call mesa_error(__FILE__,__LINE__)
            end if
         end if
         
         T6 = T * 1d-6
         R = Rho / (T6*T6*T6)

         if (co_enhanced) then
            call eval_opal_type2 (Zbase, X, dXC, dXO, T6, R, logK, data_dir)
         else if (.not. lowT_flag) then
            if (dbg) write(*,1) 'call eval_opal_type1', Zbase, X, T6, R
            call eval_opal_type1 (Zbase, X, T6, R, logK, type1_filename, OP_file, dbg)
            if (dbg) then ! .or. is_bad(logK)) then
               write(*,1) 'eval_opal_type1 Zbase', Zbase
               write(*,1) 'eval_opal_type1 X', X
               write(*,1) 'eval_opal_type1 T6', T6
               write(*,1) 'eval_opal_type1 R', R
               write(*,1) 'eval_opal_type1 logT', logT
               write(*,1) 'eval_opal_type1 logR', logR
               write(*,1) 'eval_opal_type1 logK', logK
               write(*,*) 'is_bad(logK)', is_bad(logK)
               stop ! dbg
            end if
         else if (Freedman_flag) then
            call get_lowT_blend
         else
            call eval_ferg (Zbase, X, T6, R, logK, dbg)
         end if
         
         if (is_bad(logK)) then
            if (Zbase < 0.99d0) then
               write(*,1) 'bad logK', logK, X, Zbase, logT, logRho, R, T6
               call mesa_error(__FILE__,__LINE__)
            end if
            logK = -2 ! patch for interpolation problems 
         !else
         !   write(*,1) 'good logK', logK, X, Zbase, logT, logRho, R, T6
         end if   
         
         opac = max(1d-99, 10**logK)
         logK = log10(max(1d-99,opac))
         
         if (dbg) write(*,1) 'eval_opal_type1 opac', opac
         if (dbg) write(*,1) 'eval_opal_type1 logK', logK
         
         contains
         
         
         subroutine get_lowT_blend
            use freedman, only: eval_freedman
            real(dp), parameter :: logT_low1 = 3.58d0 ! T ~ 3800
            ! logT_low1 -- below this, all Freedman
            real(dp), parameter :: logT_low2 = 3.60d0 ! T ~ 4000
            ! logT_low2 -- above this, all Ferguson
            real(dp) :: fac, logK_for_blend_with_freedman, logK_freedman
            
            if (logT <= logT_low1 .or. R > 0) then ! all Freedman
               call eval_freedman(T, Rho, logK, dbg)
               return
            end if
            
            if (logT >= logT_low2) then ! all Ferguson
               call eval_ferg (Zbase, X, T6, R, logK, dbg)
               return
            end if
            
            call eval_ferg (Zbase, X, T6, R, logK_for_blend_with_freedman, dbg)
            call eval_freedman(T, Rho, logK_freedman, dbg)
            fac = (logT - logT_low1) / (logT_low2 - logT_low1)
            fac = 0.5d0 * (1.d00 - cos(fac * pi))
            logK = fac * logK_for_blend_with_freedman + (1 - fac) * logK_freedman

         end subroutine get_lowT_blend
         
      end subroutine Get_Results


!borrowed from kap_def
      subroutine get_output_string(x,xstr,ierr) !works with X and Z
         real(dp), intent(in) :: x
         character(len=*), intent(out) :: xstr
         integer, intent(out) :: ierr
         character(len=1), parameter :: c(9)=['1','2','3','4','5','6','7','8','9']
         character(len=8) :: str
         integer :: i, j, k
         
         if(X < 0d0.or.X>1d0) then
            xstr='bad'
            ierr=-1
            return     
         endif
         
         ierr=0
   
         write(str,'(f8.6)') X
         k=0
         do i=1,9
            j=index(str,c(i),back=.true.)
            k=max(k,j)
         enddo
         xstr=str(1:max(k,3))
      end subroutine get_output_string
      
      end module kap_support
