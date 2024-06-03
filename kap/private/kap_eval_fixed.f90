! ***********************************************************************
!
!   Copyright (C) 2012-2019  The MESA Team
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

      module kap_eval_fixed

      use kap_eval_support
      use const_def, only: dp, ln10
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none
      
      contains
      
      
      subroutine Get1_kap_fixed_metal_Results( &
               z_tables, num_Zs, rq, Z, X, Rho, logRho, T, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use const_def
         
         ! INPUT
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         integer, intent(in) :: num_Zs
         type (Kap_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X
         real(dp), intent(inout) :: Rho, logRho ! can be modified to clip to table boundaries
         real(dp), intent(inout) :: T, logT

         ! OUTPUT
         real(dp), intent(out) :: logKap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.
      
         integer :: iz, i
         real(dp) :: Z0, Z1, alfa, beta, lnZ, lnZ0, lnZ1
         real(dp) :: K0, logK0, dlogK0_dlogRho, dlogK0_dlogT
         real(dp) :: logK1, dlogK1_dlogRho, dlogK1_dlogT
         real(dp) :: res
         character (len=256) :: message

         logical :: dbg
         
         include 'formats'
         
         dbg = .false.
         if (dbg) write(*,1) 'Get1_kap_fixed_metal_Results logT', logT

         ierr = 0
         logKap = 0d0
         dlnkap_dlnRho = 0d0
         dlnkap_dlnT = 0d0
         
         if (num_Zs > 1) then
            if (z_tables(1)% Z >= z_tables(2)% Z) then
               ierr = -3
               return
            end if
         end if         
         
         if (num_Zs == 1 .or. Z >= z_tables(num_Zs)% Z) then ! use the largest Z
            call Get_Kap_for_X( &
               z_tables, rq, num_Zs, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            if (dbg) write(*,1) 'logKap logT', logKap, logT
            return
         end if
         
         if (Z <= z_tables(1)% Z) then ! use the smallest Z
            if (dbg) then
               write(*,*) 'use the smallest Z'
               write(*,*) 'Z', Z
               write(*,*) 'z_tables(1)% Z', z_tables(1)% Z
            end if
            call Get_Kap_for_X( &
               z_tables, rq, 1, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            return
         end if
         
         if (dbg) then
            do iz = 1, num_Zs
               write(*,*) 'Z(iz)', iz, z_tables(iz)% Z
            end do
         end if

         do iz = 1, num_Zs-1
            if (Z < z_tables(iz+1)% Z) exit
         end do
      
         Z0 = z_tables(iz)% Z
         Z1 = z_tables(iz+1)% Z
         
         if (dbg) then
            write(*,*) 'Z0', Z0
            write(*,*) 'Z ', Z
            write(*,*) 'Z1', Z1
         end if
         
         if (Z1 <= Z0) then
            ierr = 1
            return
         end if
         
         if (Z <= Z0) then ! use the Z0 table
            if (dbg) then
               write(*,*) 'use the Z0 table', iz, Z0, Z, Z-Z0
               do i = 1, z_tables(iz)% num_Xs
                  write(*,*) 'X', i, z_tables(iz)% x_tables(i)% X
               end do
            end if
            call Get_Kap_for_X( &
               z_tables, rq, iz, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            return
         end if
         
         if (Z >= Z1) then ! use the Z1 table
            if (dbg) write(*,*) 'use the Z1 table', Z1
            call Get_Kap_for_X( &
               z_tables, rq, iz+1, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            return
         end if
         
         if (num_Zs >= 4 .and. rq% cubic_interpolation_in_Z) then
            if (dbg) write(*,*) 'call Get_Kap_for_Z_cubic'
            call Get_Kap_for_Z_cubic(z_tables, num_Zs, rq, iz, Z, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'failed in Get_Kap_for_Z_cubic'
               return
            end if
         else ! linear
            if (dbg) write(*,*) 'call Get_Kap_for_Z_linear'
            call Get_Kap_for_Z_linear(z_tables, rq, iz, Z, Z0, Z1, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'failed in Get_Kap_for_Z_linear'
               return
            end if
         end if

         if (dbg) then
            write(*,1) 'final logK at X and Z', logKap, logT, logRho, X, Z
            write(*,'(A)')
         end if

      end subroutine Get1_kap_fixed_metal_Results
      
      
      ! use tables iz-1 to iz+2 to do piecewise monotonic cubic interpolation in Z
      subroutine Get_Kap_for_Z_cubic( &
            z_tables, num_Zs, rq, iz, Z, X, logRho, logT, &
            logK, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interpolate_vector_autodiff, interp_pm_autodiff
         use auto_diff

         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: num_Zs, iz
         real(dp), intent(in) :: Z, X
         real(dp), intent(inout) :: logK, logRho, logT
         real(dp), intent(out) :: dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr
         
         integer, parameter :: n_old = 4, n_new = 1
         real(dp), dimension(n_old) :: logKs, dlogKs_dlogRho, dlogKs_dlogT
         type(auto_diff_real_2var_order1), dimension(n_old) :: logKs_ad
         type(auto_diff_real_2var_order1) :: logK_ad
         type(auto_diff_real_2var_order1) :: z_old(n_old), z_new(n_new)
         type(auto_diff_real_2var_order1), target :: work_ary(n_old*pm_work_size)
         type(auto_diff_real_2var_order1), pointer :: work(:)
         integer :: i, i1, izz

         logical, parameter :: dbg = .false.
         
         11 format(a40,e20.10)
         
         ierr = 0
         work => work_ary
         
         if (iz+2 > num_Zs) then
            i1 = num_Zs-2
         else if (iz == 1) then
            i1 = 2
         else
            i1 = iz
         end if
         
         if (dbg) write(*,*) 'iz', iz
         
         do i=1,n_old
            izz = i1-2+i
            if (dbg) write(*,*) 'izz', izz
            z_old(i) %val = z_tables(izz)% Z
            z_old(i) % d1val1 = 0d0
            z_old(i) % d1val2 = 0d0
            call Get_Kap_for_X( &
               z_tables, rq, izz, X, logRho, logT, &
               logKs(i), dlogKs_dlogRho(i), dlogKs_dlogT(i), ierr)
            if (ierr /= 0) then
               if (dbg) write(*,11) 'logRho', logRho
               if (dbg) write(*,11) 'logT', logT
               return
            end if
            ! now pack into auto_diff form
            logKs_ad(i) % val = logKs(i)
            logKs_ad(i) % d1val1 = dlogKs_dlogT(i)
            logKs_ad(i) % d1val2 = dlogKs_dlogRho(i)
         end do
         z_new(1) % val = Z
         z_new(1) % d1val1 = 0d0
         z_new(1) % d1val2 = 0d0

         call interp1(logKs_ad, logK_ad, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__,'failed in interp1 for logK')
            return
         end if

         ! unpack auto_diff pack into output reals
         logK = logK_ad % val
         dlnkap_dlnT = logK_ad % d1val1
         dlnkap_dlnRho = logK_ad % d1val2

         if (dbg) then
         
            do i=1,n_old
               write(*,*) 'z_old(i)', z_old(i)
            end do
            write(*,*) 'z_new(1)', z_new(1)
            write(*,'(A)')
         
            do i=1,n_old
               write(*,*) 'logKs(i)', logKs(i)
            end do
            write(*,*) 'logK', logK
            write(*,'(A)')

            do i=1,n_old
               write(*,*) 'dlogKs_dlogRho(i)', dlogKs_dlogRho(i)
            end do
            write(*,*) 'dlnkap_dlnRho', dlnkap_dlnRho
            write(*,'(A)')
         
            do i=1,n_old
               write(*,*) 'dlogKs_dlogT(i)', dlogKs_dlogT(i)
            end do
            write(*,*) 'dlnkap_dlnT', dlnkap_dlnT
            write(*,'(A)')
            
         end if
         
         contains
         
         subroutine interp1(old, new, ierr)
            type(auto_diff_real_2var_order1), intent(in) :: old(n_old)
            type(auto_diff_real_2var_order1), intent(out) :: new
            integer, intent(out) :: ierr
            type(auto_diff_real_2var_order1) :: v_old(n_old), v_new(n_new)
            integer :: i
            do i = 1, n_old
               v_old(i) = old(i)
            end do
            call interpolate_vector_autodiff( &
                  n_old, z_old, n_new, z_new, v_old, v_new, interp_pm_autodiff, pm_work_size, work, &
                  'Get_Kap_for_Z_cubic', ierr)
            new = v_new(1)
         end subroutine interp1
      
      end subroutine Get_Kap_for_Z_cubic
      
      
      subroutine Get_Kap_for_Z_linear(z_tables, rq, iz, Z, Z0, Z1, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz
         real(dp), intent(in) :: Z, Z0, Z1, X
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logKap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr

         real(dp) :: logK0, dlogK0_dlogRho, dlogK0_dlogT
         real(dp) :: logK1, dlogK1_dlogRho, dlogK1_dlogT
         real(dp) :: alfa, beta
      
         logical, parameter :: dbg = .false.

         ierr = 0
         call Get_Kap_for_X( &
            z_tables, rq, iz, X, logRho, logT, &
            logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
         if (ierr /= 0) return
         if (dbg) write(*,*) 'logK0', logK0
      
         call Get_Kap_for_X( &
            z_tables, rq, iz+1, X, logRho, logT, &
            logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) return
         if (dbg) write(*,*) 'logK1', logK1

         ! Z0 result in logK0, Z1 result in logK1
         beta = (Z - Z1) / (Z0 - Z1) ! beta -> 1 as Z -> Z0
         alfa = 1d0 - beta
         
         logKap = beta * logK0 + alfa * logK1
         dlnkap_dlnRho = beta * dlogK0_dlogRho + alfa * dlogK1_dlogRho
         dlnkap_dlnT = beta * dlogK0_dlogT + alfa * dlogK1_dlogT

      end subroutine Get_Kap_for_Z_linear
      
      
      subroutine Get_Kap_for_X(z_tables, rq, iz, X, logRho, logT, &
               logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz
         real(dp), intent(in) :: X
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         type (Kap_X_Table), dimension(:), pointer :: x_tables
         real(dp) :: logK0, dlogK0_dlogRho, dlogK0_dlogT
         real(dp) :: logK1, dlogK1_dlogRho, dlogK1_dlogT
         real(dp) :: X0, X1
         integer :: ix, i, num_Xs
      
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         x_tables => z_tables(iz)% x_tables
         num_Xs = z_tables(iz)% num_Xs
                  
         if (X < 0 .or. X > 1) then
            ierr = -3
            return
         end if
         
         if (num_Xs > 1) then
            if (x_tables(1)% X >= x_tables(2)% X) then
               ierr = -3
               write(*,*) 'x_tables must have increasing X values for Get_Kap_for_X'
               write(*,*) 'lowT', z_tables(iz)% lowT_flag
               write(*,2) 'Z', iz, z_tables(iz)% Z
               write(*,2) 'num_Xs', num_Xs
               do i=1,num_Xs
                  write(*,2) 'X', i, x_tables(i)% X
               end do
               stop
               
               return
            end if
         end if
         
         if (num_Xs == 1 .or. X <= x_tables(1)% X) then ! use the first X
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, 1, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         if (X >= x_tables(num_Xs)% X) then ! use the last X
            if (dbg) write(*,*) 'use the last X: call Get_Kap_for_logRho_logT'
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, num_Xs, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return    
         end if

         ! search for the X
         !write(*,*) 'search for the X'
         ix = num_Xs
         do i = 1, num_Xs-1
            if (X < x_tables(i+1)% X) then
               ix = i; exit
            end if
         end do
         
         if (ix == num_Xs) then
            if (dbg) write(*,*) 'ix == num_Xs: call Get_Kap_for_logRho_logT'
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, num_Xs, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         X0 = x_tables(ix)% X
         X1 = x_tables(ix+1)% X
         
         if (X1 <= X0) then
            ierr = 1
            return
         end if
         
         if (X0 >= X) then ! use the X0 table
            if (dbg) write(*,*) 'use the X0 table'
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, ix, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         if (X1 <= X) then ! use the X1 table
            if (dbg) write(*,*) 'use the X1 table'
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, ix+1, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if

         if (num_Xs >= 4 .and. rq% cubic_interpolation_in_X) then
            !write(*,*) 'call Get_Kap_for_X_cubic'
            call Get_Kap_for_X_cubic( &
               z_tables, rq, iz, ix, num_Xs, x_tables, X, logRho, logT, &
               logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            if (ierr /= 0) then
               !write(*,*) 'failed in Get_Kap_for_X_cubic'
               return
            end if
         else ! linear
            !write(*,*) 'call Get_Kap_for_X_linear'
            call Get_Kap_for_X_linear( &
               z_tables, rq, iz, ix, x_tables, X, X0, X1, logRho, logT, &
               logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            if (ierr /= 0) then
               !write(*,*) 'failed in Get_Kap_for_X_linear'
               return
            end if
         end if
         
         if (.false.) then
            write(*,1) 'logK at X for Z', logK, logT, logRho, X, z_tables(iz)% Z
            write(*,'(A)')
         end if
         
      end subroutine Get_Kap_for_X
      
      
      ! use tables ix-1 to ix+2 to do piecewise monotonic cubic interpolation in X
      subroutine Get_Kap_for_X_cubic( &
            z_tables, rq, iz, ix, num_Xs, x_tables, X, logRho, logT, &
            logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interpolate_vector_autodiff, interp_pm_autodiff
         use auto_diff

         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz, ix, num_Xs
         type (Kap_X_Table), dimension(:), pointer :: x_tables
         real(dp), intent(in) :: X
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         integer, parameter :: n_old = 4, n_new = 1
         real(dp), dimension(n_old) :: logKs, dlogKs_dlogRho, dlogKs_dlogT
         type(auto_diff_real_2var_order1), dimension(n_old) :: logKs_ad
         type(auto_diff_real_2var_order1) :: logK_ad
         type(auto_diff_real_2var_order1) :: x_old(n_old), x_new(n_new)
         type(auto_diff_real_2var_order1), target :: work_ary(n_old*pm_work_size)
         type(auto_diff_real_2var_order1), pointer :: work(:)
         integer :: i, i1, ixx

         logical, parameter :: dbg = .false.
         
         11 format(a40,e20.10)
         
         ierr = 0
         work => work_ary
         
         if (ix+2 > num_Xs) then
            i1 = num_Xs-2
         else if (ix == 1) then
            i1 = 2
         else
            i1 = ix
         end if
         
         if (dbg) write(*,*) 'ix', ix
         
         do i=1,n_old
            ixx = i1-2+i
            if (dbg) write(*,*) 'ixx', ixx
            x_old(i) % val = x_tables(ixx)% X
            x_old(i) % d1val1 = 0d0
            x_old(i) % d1val2 = 0d0
            call Get_Kap_for_logRho_logT( &
                     z_tables, rq, iz, x_tables, ixx, &
                     logRho, logT, logKs(i), dlogKs_dlogRho(i), dlogKs_dlogT(i), ierr)
            if (ierr /= 0) then
               if (dbg) write(*,11) 'logRho', logRho
               if (dbg) write(*,11) 'logT', logT
               return
            end if
            ! now pack into auto_diff form
            logKs_ad(i) % val = logKs(i)
            logKs_ad(i) % d1val1 = dlogKs_dlogT(i)
            logKs_ad(i) % d1val2 = dlogKs_dlogRho(i)
         end do
         x_new(1) % val = X
         x_new(1) % d1val1 = 0d0
         x_new(1) % d1val2 = 0d0
                  
         call interp1(logKs_ad, logK_ad, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__,'failed in interp1 for logK')
            return
         end if

         ! unpack auto_diff pack into output reals
         logK = logK_ad % val
         dlogK_dlogT = logK_ad % d1val1
         dlogK_dlogRho = logK_ad % d1val2
         
         if (dbg) then
         
            do i=1,n_old
               write(*,*) 'x_old(i)', x_old(i)
            end do
            write(*,*) 'x_new(1)', x_new(1)
            write(*,'(A)')
         
            do i=1,n_old
               write(*,*) 'logKs(i)', logKs(i)
            end do
            write(*,*) 'logK', logK
            write(*,'(A)')

            do i=1,n_old
               write(*,*) 'dlogKs_dlogRho(i)', dlogKs_dlogRho(i)
            end do
            write(*,*) 'dlogK_dlogRho', dlogK_dlogRho
            write(*,'(A)')
         
            do i=1,n_old
               write(*,*) 'dlogKs_dlogT(i)', dlogKs_dlogT(i)
            end do
            write(*,*) 'dlogK_dlogT', dlogK_dlogT
            write(*,'(A)')
            
         end if
         
         contains
         
         subroutine interp1(old, new, ierr)
            type(auto_diff_real_2var_order1), intent(in) :: old(n_old)
            type(auto_diff_real_2var_order1), intent(out) :: new
            integer, intent(out) :: ierr
            type(auto_diff_real_2var_order1) :: v_old(n_old), v_new(n_new)
            integer :: i
            do i = 1, n_old
               v_old(i) = old(i)
            end do
            call interpolate_vector_autodiff( &
                  n_old, x_old, n_new, x_new, v_old, v_new, interp_pm_autodiff, pm_work_size, work, &
                  'Get_Kap_for_X_cubic', ierr)
            new = v_new(1)
         end subroutine interp1
          
      end subroutine Get_Kap_for_X_cubic
      
      
      ! use tables ix and ix+1 to do linear interpolation in X
      subroutine Get_Kap_for_X_linear( &
            z_tables, rq, iz, ix, x_tables, X, X0, X1, logRho, logT, &
            logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz, ix
         type (Kap_X_Table), dimension(:), pointer :: x_tables
         real(dp), intent(in) :: X, X0, X1
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         real(dp) :: logK0, dlogK0_dlogRho, dlogK0_dlogT
         real(dp) :: logK1, dlogK1_dlogRho, dlogK1_dlogT
         real(dp) :: alfa, beta
         
         ierr = 0
         
         call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, ix, &
                  logRho, logT, logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
         if (ierr /= 0) return
      
         call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, ix+1, &
                  logRho, logT, logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) return
         
         ! X0 result in logK0, X1 result in logK1
         beta = (X - X1) / (X0 - X1) ! beta -> 1 as X -> X0
         alfa = 1d0 - beta
         
         logK = beta * logK0 + alfa * logK1
         dlogK_dlogRho = beta * dlogK0_dlogRho + alfa * dlogK1_dlogRho
         dlogK_dlogT = beta * dlogK0_dlogT + alfa * dlogK1_dlogT
      
      end subroutine Get_Kap_for_X_linear

      
      subroutine Get_Kap_for_logRho_logT( &
               z_tables, rq, iz, x_tables, ix, &
               logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use load_kap, only: load_one
         use kap_def
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         type (Kap_X_Table), dimension(:), pointer :: x_tables
         integer :: iz, ix
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr

         real(dp) :: logR0, logR1, logT0, logT1, logR, logR_in
         real(dp) :: df_dx, df_dy
         integer :: iR, jtemp, i, num_logRs, num_logTs
         logical :: clipped_logR
         logical, parameter :: read_later = .false.
      
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         if (x_tables(ix)% not_loaded_yet) then ! avoid doing critical section if possible
!$omp critical (load_kap_x_table)
            if (x_tables(ix)% not_loaded_yet) then
               call load_one(rq, &
                  z_tables, iz, ix, x_tables(ix)% X, x_tables(ix)% Z, &
                  read_later, ierr)
            end if
!$omp end critical (load_kap_x_table)
         end if
         if (ierr /= 0) then
            !call mesa_error(__FILE__,__LINE__,'load_one failed in Get_Kap_for_logRho_logT')
            return
         end if

         ! logR from inputs
         logR_in = logRho - 3d0*logT + 18d0
         
         ! blends at higher levels MUST prevent
         ! these tables from being called off their
         ! high/low T and low R edges

         if (logT > x_tables(ix)% logT_max) then
            ierr = -1
            return
         end if

         if (logT < x_tables(ix)% logT_min) then
            ierr = -1
            return
         end if

         if (logR_in < x_tables(ix)% logR_min) then
            ierr = -1
            return
         end if

         
         ! off the high R edge, we use the input temperature
         ! but clip logR to the table edge value

         if (logR_in > x_tables(ix)% logR_max) then
            logR = x_tables(ix)% logR_max
            clipped_logR = .true.
         else
            logR = logR_in
            clipped_logR = .false.
         end if


         num_logRs = x_tables(ix)% num_logRs
         num_logTs = x_tables(ix)% num_logTs

         if (num_logRs <= 0) then
            write(*,*) 'num_logRs', num_logRs
            write(*,*) 'ix', ix
            call mesa_error(__FILE__,__LINE__,'Get_Kap_for_logRho_logT')
         end if

         if (num_logTs <= 0) then
            write(*,*) 'num_logTs', num_logRs
            write(*,*) 'ix', ix
            call mesa_error(__FILE__,__LINE__,'Get_Kap_for_logRho_logT')
         end if

         call Locate_logR( &
            rq, num_logRs, x_tables(ix)% logR_min, x_tables(ix)% logR_max, &
            x_tables(ix)% ili_logRs, x_tables(ix)% logRs, logR, iR, logR0, logR1, ierr)
         if (ierr /= 0) then
            write(*,1) 'x_tables(ix)% logR_min', x_tables(ix)% logR_min
            write(*,1) 'x_tables(ix)% logR_max', x_tables(ix)% logR_max
            write(*,2) 'num_logRs', num_logRs
            write(*,2) 'iR', iR
            write(*,1) 'logR', logR
            write(*,1) 'logR0', logR0
            write(*,1) 'logR1', logR1
            do i=1,num_logRs
               write(*,2) 'logR', i, x_tables(ix)% logRs(i)
            end do
            write(*,*) 'clip_to_kap_table_boundaries', clip_to_kap_table_boundaries
            call mesa_error(__FILE__,__LINE__,'failed in Locate_logR')
            return
         end if

         call Locate_logT( &
            rq, num_logTs, x_tables(ix)% logT_min, x_tables(ix)% logT_max, &
            x_tables(ix)% ili_logTs, x_tables(ix)% logTs, logT, jtemp, logT0, logT1, ierr)
         if (ierr /= 0) return

         call Do_Kap_Interpolations( &
            x_tables(ix)% kap1, num_logRs, num_logTs, &
            iR, jtemp, logR0, logR, logR1, logT0, logT, logT1, &
            logK, df_dx, df_dy)
         if (clipped_logR) df_dx = 0
         
         if (dbg) write(*,1) 'Do_Kap_Interpolations: logK', logK

         ! convert df_dx and df_dy to dlogK_dlogRho, dlogK_dlogT
         dlogK_dlogRho = df_dx
         dlogK_dlogT = df_dy - 3d0*df_dx

         if (dbg) then
            write(*,1) 'logK', logK, logT, logRho, x_tables(ix)% X, z_tables(iz)% Z
         end if
         
      end subroutine Get_Kap_for_logRho_logT
      

      end module kap_eval_fixed
