! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

      module kap_eval_co
      use utils_lib,only: is_bad, mesa_error
      use kap_eval_support
      use const_def, only: dp
      use math_lib
      
      implicit none
            
      contains


      subroutine Get1_kap_CO_Results( &
               rq, Zbase, X, dXC, dXO, Rho, logRho, T, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use const_def
         
         ! INPUT
         type (Kap_General_Info), pointer :: rq
         real(dp), intent(in) :: Zbase, X, dXC, dXO
         real(dp), intent(inout) :: Rho, logRho ! can be modified to clip to table boundaries
         real(dp), intent(inout) :: T, logT

         ! OUTPUT
         real(dp), intent(out) :: logKap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.
         
         integer :: iz, use_iz, num_Zs, CO_option
         real(dp) :: Z0, Z1,log10_Zbase, log10_Z0, log10_Z1
         real(dp) ::  alfa, beta
         real(dp) :: logK0, dlogK0_dlogRho, dlogK0_dlogT, logK1, dlogK1_dlogRho, dlogK1_dlogT
         character (len=256) :: message
      
         logical, parameter :: use_closest_Z = .false.

         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         CO_option = rq% kap_CO_option
         num_Zs = num_kap_CO_Zs(CO_option)

         if (num_Zs > 1) then
            if (kap_co_z_tables(CO_option)% ar(1)% Zbase >= kap_co_z_tables(CO_option)% ar(2)% Zbase) then
               ierr = -3
               return
            end if
         end if         
         
         if (num_Zs == 1 .or. &
               Zbase >= kap_co_z_tables(CO_option)% ar(num_Zs)% Zbase) then ! use the largest Zbase
            if (dbg) write(*,*) 'use the largest Zbase', &
               num_kap_CO_Zs(CO_option), kap_co_z_tables(CO_option)% ar(num_Zs)% Zbase
            call Get_Kap_for_CO_X( &
               rq, dXC, dXO, num_Zs, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            return
         end if
         
         if (Zbase <= kap_co_z_tables(CO_option)% ar(1)% Zbase) then ! use the smallest Zbase
            if (dbg) then
               write(*,*) 'use the smallest Zbase'
               write(*,*) 'Zbase', Zbase
               write(*,*) 'kap_co_z_tables(1)% Zbase', kap_co_z_tables(CO_option)% ar(1)% Zbase
            end if
            call Get_Kap_for_CO_X( &
               rq, dXC, dXO, 1, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            return
         end if

         do iz = 1, num_Zs-1
            if (Zbase < kap_co_z_tables(CO_option)% ar(iz+1)% Zbase) exit
         end do
      
         Z0 = kap_co_z_tables(CO_option)% ar(iz)% Zbase
         Z1 = kap_co_z_tables(CO_option)% ar(iz+1)% Zbase
         
         if (Zbase <= Z0) then ! use the Z0 table
            if (dbg) write(*,*) 'use the Z0 table', Z0
            call Get_Kap_for_CO_X( &
               rq, dXC, dXO, iz, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            return
         end if
         
         if (Zbase >= Z1) then ! use the Z1 table
            if (dbg) write(*,*) 'use the Z1 table', Z1
            call Get_Kap_for_CO_X( &
               rq, dXC, dXO, iz+1, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            return
         end if
         
         if (use_closest_Z) then
            log10_Z0 = kap_co_z_tables(CO_option)% ar(iz)% log10_Zbase
            log10_Z1 = kap_co_z_tables(CO_option)% ar(iz+1)% log10_Zbase
            log10_Zbase = log10(dble(Zbase))
            if (log10_Z1 - log10_Zbase > log10_Zbase - log10_Z0) then ! use the Z0 table
               use_iz = iz
            else
               use_iz = iz+1
            end if
            if (dbg) write(*,*) 'use the Z0 table', Z0
            call Get_Kap_for_CO_X( &
               rq, dXC, dXO, use_iz, X, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            return
         end if
         
         if (dbg) then
            write(*,*) 'iz', iz
            write(*,*) '   Z0', Z0
            write(*,*) 'Zbase', Zbase
            write(*,*) '   Z1', Z1
            write(*,'(A)')
         end if
         
         if (num_Zs >= 4 .and. rq% cubic_interpolation_in_Z) then
            if (dbg) write(*,*) 'call Get_Kap_for_CO_Z_cubic'
            call Get_Kap_for_CO_Z_cubic(rq, iz, Zbase, X, dXC, dXO, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'failed in Get_Kap_for_CO_Z_cubic'
               return
            end if
         else ! linear
            if (dbg) write(*,*) 'call Get_Kap_for_CO_Z_linear'
            call Get_Kap_for_CO_Z_linear(rq, iz, Zbase, Z0, Z1, X, dXC, dXO, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'failed in Get_Kap_for_CO_Z_linear'
               return
            end if
         end if

      end subroutine Get1_kap_CO_Results
      
      
      ! use tables iz-1 to iz+2 to do piecewise monotonic cubic interpolation in Z
      subroutine Get_Kap_for_CO_Z_cubic(rq, iz, Z, X, dXC, dXO, logRho, logT, &
               logK, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interpolate_vector, interp_pm

         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz
         real(dp), intent(in) :: Z, X, dXC, dXO
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr
         
         integer, parameter :: n_old = 4, n_new = 1
         real(dp), dimension(n_old) :: logKs, dlogKs_dlogRho, dlogKs_dlogT
         real(dp) :: z_old(n_old), z_new(n_new)
         real(dp), target :: work_ary(n_old*pm_work_size)
         real(dp), pointer :: work(:)
         integer :: i, i1, izz, num_Zs, CO_option
         
         logical, parameter :: dbg = .false.
         
         11 format(a40,e20.10)
         
         ierr = 0
         work => work_ary
         CO_option = rq% kap_CO_option
         num_Zs = num_kap_CO_Zs(CO_option)
         
         if (iz+2 > num_Zs) then
            i1 = num_Zs-2
         else if (iz == 1) then
            i1 = 2
         else
            i1 = iz
         end if
         
         if (dbg) then
            write(*,*) 'n_old', n_old
            write(*,*) 'i1', i1
            write(*,*) 'iz', iz
            write(*,*) 'Z', Z
            write(*,'(A)')
         end if
         
         do i=1,n_old
            izz = i1-2+i
            z_old(i) = kap_co_z_tables(CO_option)% ar(izz)% Zbase
            if (dbg) then
               write(*,*) 'izz', izz
               write(*,*) 'z_old', i, z_old(i)
            end if
            call Get_Kap_for_CO_X( &
               rq, dXC, dXO, izz, X, logRho, logT, &
               logKs(i), dlogKs_dlogRho(i), dlogKs_dlogT(i), ierr)
            if (dbg) write(*,11) 'logK', logKs(i)
            if (ierr /= 0) then
               return
            end if
         end do
         z_new(1) = Z
         
         call interp1(logKs, logK, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__,'failed in interp1 for logK')
            return
         end if
         
         call interp1(dlogKs_dlogRho, dlnkap_dlnRho, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__,'failed in interp1 for dlogK_dlogRho')
            return
         end if
                  
         call interp1(dlogKs_dlogT, dlnkap_dlnT, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__,'failed in interp1 for dlogK_dlogT')
            return
         end if
         
         if (dbg) then
         
            do i=1,n_old
               write(*,*) 'z_old(i)', z_old(i)
            end do
            write(*,'(A)')
            write(*,*) 'z_new(1)', z_new(1)
            write(*,'(A)')
         
            do i=1,n_old
               write(*,*) 'logK', i, logKs(i)
            end do
            write(*,'(A)')
            write(*,*) 'logK', logK
            write(*,'(A)')

            do i=1,n_old
               write(*,*) 'dlogKs_dlogRho', i, dlogKs_dlogRho(i)
            end do
            write(*,'(A)')
            write(*,*) 'dlnkap_dlnRho', dlnkap_dlnRho
            write(*,'(A)')
         
            do i=1,n_old
               write(*,*) 'dlogKs_dlogT', i, dlogKs_dlogT(i)
            end do
            write(*,'(A)')
            write(*,*) 'dlnkap_dlnT', dlnkap_dlnT
            write(*,'(A)')
            
         end if
         
         contains
         
         subroutine interp1(old, new, ierr)
            real(dp), intent(in) :: old(n_old)
            real(dp), intent(out) :: new
            integer, intent(out) :: ierr
            real(dp) :: v_old(n_old), v_new(n_new)
            v_old(:) = dble(old(:))
            call interpolate_vector( &
               n_old, z_old, n_new, z_new, v_old, v_new, interp_pm, pm_work_size, work, &
               'Get_Kap_for_CO_Z_cubic', ierr)
            new = v_new(1)
         end subroutine interp1
      
      end subroutine Get_Kap_for_CO_Z_cubic
      
      
      subroutine Get_Kap_for_CO_Z_linear( &
               rq, iz, Z, Z0, Z1, X, dXC, dXO, logRho, logT, &
               logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz
         real(dp), intent(in) :: Z, Z0, Z1, X, dXC, dXO
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logKap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr

         real(dp) :: logK0, dlogK0_dlogRho, dlogK0_dlogT, logK1, dlogK1_dlogRho, dlogK1_dlogT
         real(dp) :: alfa, beta
      
         logical, parameter :: dbg = .false.

         ierr = 0
         if (dbg) write(*,*) 'call Get_Kap_for_CO_X'
         call Get_Kap_for_CO_X( &
            rq, dXC, dXO, iz, X, logRho, logT, &
            logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
         if (ierr /= 0) return
         if (dbg) write(*,*) 'logK0', logK0
      
         call Get_Kap_for_CO_X( &
            rq, dXC, dXO, iz+1, X, logRho, logT, &
            logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) return
         if (dbg) write(*,*) 'logK1', logK1

         ! Z0 result in logK0, Z1 result in logK1
         beta = (Z - Z1) / (Z0 - Z1) ! beta -> 1 as Z -> Z0
         alfa = 1d0 - beta
         logKap = beta*logK0 + alfa*logK1
         dlnkap_dlnRho = beta*dlogK0_dlogRho + alfa*dlogK1_dlogRho
         dlnkap_dlnT = beta*dlogK0_dlogT + alfa*dlogK1_dlogT

      end subroutine Get_Kap_for_CO_Z_linear
      
      
      subroutine Get_Kap_for_CO_X(rq, dXC, dXO, iz, X, logRho, logT, &
               logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         ! return opacity from Z table number iz for the given X, dXC, dXO
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz
         real(dp), intent(in) :: X, dXC, dXO
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         type (Kap_CO_X_Table), dimension(:), pointer :: x_tables
         real(dp) :: logK0, dlogK0_dlogRho, dlogK0_dlogT, logK1, dlogK1_dlogRho, dlogK1_dlogT
         real(dp) :: X0, X1
         real(dp) :: alfa, beta
         integer :: ix, i, num_Xs, CO_option
         
         logical, parameter :: dbg = .false.

         CO_option = rq% kap_CO_option
         num_Xs = num_kap_CO_Xs(CO_option)
         x_tables => kap_co_z_tables(CO_option)% ar(iz)% x_tables
                  
         if (X < 0 .or. X > 1) then
            ierr = -3
            return
         end if
         
         if (num_Xs > 1) then
            if (x_tables(1)% X >= x_tables(2)% X) then
               ierr = -3
               return
            end if
         end if

         if (X >= x_tables(num_Xs)% X) then ! use the last X
            if (dbg) write(*,*) 'use the last X'
            call Get_Kap_for_dXCO( &
               rq, iz, x_tables, dXC, dXO, num_Xs, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return    
         end if
         
         if (X <= x_tables(1)% X) then ! use the first X
            if (dbg) write(*,*) 'use the first X'
            call Get_Kap_for_dXCO( &
                  rq, iz, x_tables, dXC, dXO, 1, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if

         ! search for the X
         if (dbg) write(*,*) 'search for the X'
         ix = num_Xs
         do i = 1, num_Xs-1
            if (X < x_tables(i+1)% X) then
               ix = i; exit
            end if
         end do
         
         if (ix == num_Xs) then
            call Get_Kap_for_dXCO( &
                  rq, iz, x_tables, dXC, dXO, num_Xs, &
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
            call Get_Kap_for_dXCO( &
                  rq, iz, x_tables, dXC, dXO, ix, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         if (X1 <= X) then ! use the X1 table
            call Get_Kap_for_dXCO( &
                  rq, iz, x_tables, dXC, dXO, ix+1, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         if (num_Xs >= 4 .and. rq% cubic_interpolation_in_X) then
            call Get_Kap_for_CO_X_cubic( &
                  rq, iz, ix, dXC, dXO, x_tables, X, logRho, logT, &
                  logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         else ! linear
            call Get_Kap_for_CO_X_linear( &
                  rq, iz, ix, dXC, dXO, x_tables, X, X0, X1, logRho, logT, &
                  logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         end if
      
      end subroutine Get_Kap_for_CO_X
      
      
      ! use tables ix-1 to ix+2 to do piecewise monotonic cubic interpolation in X
      subroutine Get_Kap_for_CO_X_cubic( &
            rq, iz, ix, dXC, dXO, x_tables, X, logRho, logT, &
            logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interpolate_vector, interp_pm

         ! return opacity from Z table number iz for the given X, XC, XO
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz, ix
         type (Kap_CO_X_Table), dimension(:), pointer :: x_tables
         real(dp), intent(in) :: X, dXC, dXO
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         integer, parameter :: n_old = 4, n_new = 1
         real(dp), dimension(n_old) :: logKs, dlogKs_dlogRho, dlogKs_dlogT
         real(dp) :: x_old(n_old), x_new(n_new)
         real(dp), target :: work_ary(n_old*pm_work_size)
         real(dp), pointer :: work(:)
         integer :: i, i1, ixx, num_Xs
         
         logical, parameter :: dbg = .false.
         
         11 format(a40,e20.10)
         
         ierr = 0
         work => work_ary

         num_Xs = num_kap_CO_Xs(rq% kap_CO_option)
         
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
            x_old(i) = x_tables(ixx)% X

            call Get_Kap_for_dXCO( &
                  rq, iz, x_tables, dXC, dXO, ixx, &
                  logRho, logT, logKs(i), dlogKs_dlogRho(i), dlogKs_dlogT(i), ierr)
            if (ierr /= 0) then
               if (dbg) write(*,11) 'logRho', logRho
               if (dbg) write(*,11) 'logT', logT
               return
            end if
         end do
         x_new(1) = X
         
         call interp1(logKs, logK, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__,'failed in interp1 for logK')
            return
         end if
         
         call interp1(dlogKs_dlogRho, dlogK_dlogRho, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__,'failed in interp1 for dlogK_dlogRho')
            return
         end if
                  
         call interp1(dlogKs_dlogT, dlogK_dlogT, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__,'failed in interp1 for dlogK_dlogT')
            return
         end if
         
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
            real(dp), intent(in) :: old(n_old)
            real(dp), intent(out) :: new
            integer, intent(out) :: ierr
            real(dp) :: v_old(n_old), v_new(n_new)
            v_old(:) = dble(old(:))
            call interpolate_vector( &
                  n_old, x_old, n_new, x_new, v_old, v_new, interp_pm, pm_work_size, work, &
                  'Get_Kap_for_CO_X_cubic', ierr)
            new = real(v_new(1),kind=dp)
         end subroutine interp1
      
      end subroutine Get_Kap_for_CO_X_cubic
      
      
      subroutine Get_Kap_for_CO_X_linear( &
            rq, iz, ix, dXC, dXO, x_tables, X, X0, X1, logRho, logT, &
            logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         ! return opacity from Z table number iz for the given X, dXC, dXO
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz, ix
         type (Kap_CO_X_Table), dimension(:), pointer :: x_tables
         real(dp), intent(in) :: dXC, dXO, X, X0, X1
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         real(dp) :: logK0, dlogK0_dlogRho, dlogK0_dlogT, logK1, dlogK1_dlogRho, dlogK1_dlogT
         real(dp) :: alfa, beta
         integer :: i
         
         logical, parameter :: dbg = .false.

         ierr = 0
         call Get_Kap_for_dXCO( &
                  rq, iz, x_tables, dXC, dXO, ix, &
                  logRho, logT, logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
         if (ierr /= 0) return
      
         call Get_Kap_for_dXCO( &
                  rq, iz, x_tables, dXC, dXO, ix+1, &
                  logRho, logT, logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) return
         
         ! X0 result in logK0, X1 result in logK1
         beta = (X - X1) / (X0 - X1) ! beta -> 1 as X -> X0
         alfa = 1d0 - beta
         
         logK = beta*logK0 + alfa*logK1
         dlogK_dlogRho = beta*dlogK0_dlogRho + alfa*dlogK1_dlogRho
         dlogK_dlogT = beta*dlogK0_dlogT + alfa*dlogK1_dlogT      
      
      end subroutine Get_Kap_for_CO_X_linear


      subroutine Get_Kap_for_dXCO( &
               rq, iz, x_tables, dXC_in, dXO_in, ix, &
               logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         ! return value from xtable number ix for the given dXC and dXO
         use kap_def
         use load_CO_kap
         type (Kap_General_Info), pointer :: rq
         type (Kap_CO_X_Table), dimension(:), pointer :: x_tables
         integer :: iz, ix
         real(dp), intent(in) :: dXC_in, dXO_in
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         type (Kap_CO_Table), dimension(:), pointer :: co_tables ! stored by table number
         real(dp) :: dXC, dXO, fac, dXCO_max, Z, dXC_lookup, dXO_lookup
         integer :: num_CO_tables, num_dXC_gt_dXO, i1, i2, i3, i4
         real(dp) :: alfa, beta
         real(dp) :: logK1, dlogK1_dlogRho, dlogK1_dlogT, logK2, dlogK2_dlogRho, dlogK2_dlogT
         real(dp) :: logK3, dlogK3_dlogRho, dlogK3_dlogT, logK4, dlogK4_dlogRho, dlogK4_dlogT
         real(dp) :: dXC1_lookup, dXO1_lookup, dXC2_lookup, dXO2_lookup, dXC_2_4_lookup, dXO_2_4_lookup
         real(dp) :: dXC3_lookup, dXO3_lookup, dXC4_lookup, dXO4_lookup, dXC_1_3_lookup, dXO_1_3_lookup
         real(dp) :: logK_2_4, dlogK_2_4_dlogRho, dlogK_2_4_dlogT
         real(dp) :: logK_1_3, dlogK_1_3_dlogRho, dlogK_1_3_dlogT
         logical, parameter :: read_later = .false., dbg = .false.
         
         include 'formats'

         ierr = 0
         
         if (dbg) write(*,1) 'enter Get_Kap_for_dXCO dXC_in dXO_in', dXC_in, dXO_in
         
         dXC = max(0.0_dp, dXC_in)
         dXO = max(0.0_dp, dXO_in)
         
         if (x_tables(ix)% not_loaded_yet) then ! avoid doing critical section if possible
!$omp critical (load_co_table)
            if (x_tables(ix)% not_loaded_yet) then
               call load_one_CO(rq,kap_co_z_tables(rq% kap_CO_option)% ar,iz,ix,read_later,ierr)
            end if
!$omp end critical (load_co_table)
         end if
         if (ierr /= 0) return
         
         co_tables => x_tables(ix)% co_tables
         num_CO_tables = x_tables(ix)% num_CO_tables
         if (dbg) write(*,2) 'num_CO_tables', num_CO_tables
         
         if (num_CO_tables < 1) then
            ierr = -1
            write(*,2) 'num_CO_tables', num_CO_tables
            return
         end if

         num_dXC_gt_dXO = x_tables(ix)% num_dXC_gt_dXO
         Z = x_tables(ix)% Z

         dXCO_max = 1 - ((x_tables(ix)% X) + (x_tables(ix)% Z))
         if (dXC + dXO > dXCO_max) then
            fac = dXCO_max / (dXC + dXO)
            dXC = fac*dXC
            dXO = fac*dXO
         end if
         
         dXC_lookup = get_dX_lookup(dXC, Z)
         dXO_lookup = get_dX_lookup(dXO, Z)
         
         if (dbg) write(*,2) 'call Find_CO_Tables', ix
         call Find_CO_Tables(rq, x_tables, ix, x_tables(ix)% CO_table_numbers,  &
                     x_tables(ix)% next_dXO_table, x_tables(ix)% next_dXC_table,  &
                     co_tables, num_CO_tables, num_dXC_gt_dXO, &
                     dXCO_max, dXC, dXO, dXC_lookup, dXO_lookup, i1, i2, i3, i4,ierr)
         if (ierr /= 0) then
            write(*,*) 'kap failed in Find_CO_Tables'
            return
         endif
         
         if (i1 > 0 .and. i2 <= 0 .and. i3 <= 0 .and. i4 <= 0) then
            call Get_CO_Kap_for_logRho_logT(rq, x_tables, ix, co_tables, i1, logRho, logT,  &
                     logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
                  
         if (i1 <= 0 .or. i2 <= 0 .or. i3 <= 0) call mesa_error(__FILE__,__LINE__,'error in result from Find_CO_Tables')
         
         if (matches_table(i2)) then
            call Get_CO_Kap_for_logRho_logT(rq, x_tables, ix, co_tables, i2, logRho, logT,  &
                     logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         if (matches_table(i3)) then
            call Get_CO_Kap_for_logRho_logT(rq, x_tables, ix, co_tables, i3, logRho, logT,  &
                     logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         if (i4 > 0) then
            if (matches_table(i4)) then
               call Get_CO_Kap_for_logRho_logT(rq, x_tables, ix, co_tables, i4, logRho, logT,  &
                     logK, dlogK_dlogRho, dlogK_dlogT, ierr)
               return
            end if
         end if
         
         call Get_CO_Kap_for_logRho_logT(rq, x_tables, ix, co_tables, i1, logRho, logT,  &
                     logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) return
         dXC1_lookup = co_tables(i1)% dXC_lookup
         dXO1_lookup = co_tables(i1)% dXO_lookup
         
         call Get_CO_Kap_for_logRho_logT(rq, x_tables, ix, co_tables, i2, logRho, logT,  &
                     logK2, dlogK2_dlogRho, dlogK2_dlogT, ierr)
         if (ierr /= 0) return
         dXC2_lookup = co_tables(i2)% dXC_lookup
         dXO2_lookup = co_tables(i2)% dXO_lookup
         
         call Get_CO_Kap_for_logRho_logT(rq, x_tables, ix, co_tables, i3, logRho, logT,  &
                     logK3, dlogK3_dlogRho, dlogK3_dlogT, ierr)
         if (ierr /= 0) return
         dXC3_lookup = co_tables(i3)% dXC_lookup
         dXO3_lookup = co_tables(i3)% dXO_lookup
         
         if (i4 > 0) then
            call Get_CO_Kap_for_logRho_logT(rq, x_tables, ix, co_tables, i4, logRho, logT,  &
                     logK4, dlogK4_dlogRho, dlogK4_dlogT, ierr)
            if (ierr /= 0) return
            dXC4_lookup = co_tables(i4)% dXC_lookup
            dXO4_lookup = co_tables(i4)% dXO_lookup
         else ! copy i3 results
            logK4 = logK3
            dlogK4_dlogRho = dlogK3_dlogRho
            dlogK4_dlogT = dlogK3_dlogT
            dXC4_lookup = dXC3_lookup
            dXO4_lookup = dXO3_lookup
         end if
                     
         if (dXC >= dXO) then ! use values on lines i1-i3 and i2-i4 at dXO
         
            call Get_Kap_at_dXO(dXO_lookup, &
                           dXC2_lookup, dXO2_lookup, logK2, dlogK2_dlogRho, dlogK2_dlogT,  &
                           dXC4_lookup, dXO4_lookup, logK4, dlogK4_dlogRho, dlogK4_dlogT,  &
                           logK_2_4, dlogK_2_4_dlogRho, dlogK_2_4_dlogT, dXC_2_4_lookup)
            
            call Get_Kap_at_dXO(dXO_lookup, &
                           dXC1_lookup, dXO1_lookup, logK1, dlogK1_dlogRho, dlogK1_dlogT,  &
                           dXC3_lookup, dXO3_lookup, logK3, dlogK3_dlogRho, dlogK3_dlogT,  &
                           logK_1_3, dlogK_1_3_dlogRho, dlogK_1_3_dlogT, dXC_1_3_lookup)
            if (dXC_1_3_lookup == dXC_2_4_lookup) then
               alfa = 0d0
            else
               alfa = (dXC_lookup - dXC_2_4_lookup) / (dXC_1_3_lookup - dXC_2_4_lookup)
            end if
            
         else ! use values on lines i1-i3 and i2-i4 at dXC
         
            call Get_Kap_at_dXC(dXC_lookup, &
                           dXC2_lookup, dXO2_lookup, logK2, dlogK2_dlogRho, dlogK2_dlogT,  &
                           dXC4_lookup, dXO4_lookup, logK4, dlogK4_dlogRho, dlogK4_dlogT,  &
                           logK_2_4, dlogK_2_4_dlogRho, dlogK_2_4_dlogT, dXO_2_4_lookup)

            call Get_Kap_at_dXC(dXC_lookup, &
                           dXC1_lookup, dXO1_lookup, logK1, dlogK1_dlogRho, dlogK1_dlogT,  &
                           dXC3_lookup, dXO3_lookup, logK3, dlogK3_dlogRho, dlogK3_dlogT,  &
                           logK_1_3, dlogK_1_3_dlogRho, dlogK_1_3_dlogT, dXO_1_3_lookup)
            if (dXO_1_3_lookup == dXO_2_4_lookup) then
               alfa = 0d0
            else
               alfa = (dXO_lookup - dXO_2_4_lookup) / (dXO_1_3_lookup - dXO_2_4_lookup)
            end if
         
         end if

         beta = 1d0 - alfa

         logK = alfa*logK_1_3 + beta*logK_2_4
         dlogK_dlogRho = alfa*dlogK_1_3_dlogRho + beta*dlogK_2_4_dlogRho
         dlogK_dlogT = alfa*dlogK_1_3_dlogT + beta*dlogK_2_4_dlogT
         
         if (is_bad(logK)) then
            ierr = -1
            return
            write(*,1) 'logK', logK
            write(*,1) 'logK_1_3', logK_1_3
            write(*,1) 'logK_2_4', logK_2_4
            write(*,1) 'alfa', alfa
            write(*,1) 'beta', beta
            write(*,'(A)')
            write(*,2) 'dXC1_lookup', i1, dXC1_lookup
            write(*,2) 'dXO1_lookup', i1, dXO1_lookup
            write(*,2) 'dXC2_lookup', i2, dXC2_lookup
            write(*,2) 'dXO2_lookup', i2, dXO2_lookup
            write(*,2) 'dXC3_lookup', i3, dXC3_lookup
            write(*,2) 'dXO3_lookup', i3, dXO3_lookup
            write(*,2) 'dXC4_lookup', i4, dXC4_lookup
            write(*,2) 'dXO4_lookup', i4, dXO4_lookup
            write(*,1) 'dXC', dXC
            write(*,1) 'dXO', dXO
            call mesa_error(__FILE__,__LINE__,'Get_Kap_for_dXCO')
         end if
         
         contains       
         
         
         logical function matches_table(i)
            integer :: i
            if (i < 1 .or. i > num_CO_tables) then
               write(*,*) 'logRho', logRho
               write(*,*) 'logT', logT
               call mesa_error(__FILE__,__LINE__,'bug in kap_eval_co matches_table')
               matches_table = .false.
            else if (abs(dXC_lookup - co_tables(i)% dXC_lookup) == 0 .and.  &
                     abs(dXO_lookup - co_tables(i)% dXO_lookup) == 0) then
               matches_table = .true.
            else
               matches_table = .false.
            end if
         end function matches_table


         subroutine Get_Kap_at_dXO(dXO_lookup, &
                           dXC_a_lookup, dXO_a_lookup, logK_a, dlogK_a_dlogRho, dlogK_a_dlogT,  &
                           dXC_b_lookup, dXO_b_lookup, logK_b, dlogK_b_dlogRho, dlogK_b_dlogT,  &
                           logK_a_b, dlogK_a_b_dlogRho, dlogK_a_b_dlogT, dXC_a_b_lookup)
            real(dp), intent(in) :: dXO_lookup
            real(dp), intent(in) :: dXC_a_lookup, dXO_a_lookup
            real(dp), intent(in) :: logK_a, dlogK_a_dlogRho, dlogK_a_dlogT
            real(dp), intent(in) :: dXC_b_lookup, dXO_b_lookup
            real(dp), intent(in) :: logK_b, dlogK_b_dlogRho, dlogK_b_dlogT
            real(dp), intent(out) :: logK_a_b, dlogK_a_b_dlogRho, dlogK_a_b_dlogT
            real(dp), intent(out) :: dXC_a_b_lookup
            
            real(dp) :: alfa, beta
            
            if (dXO_a_lookup == dXO_b_lookup) then
               alfa = 0d0
            else
               alfa = (dXO_lookup - dXO_b_lookup) / (dXO_a_lookup - dXO_b_lookup)
            end if
               
            dXC_a_b_lookup = dXC_b_lookup + (dXC_a_lookup - dXC_b_lookup)*alfa
            beta = 1d0 - alfa
            logK_a_b = alfa*logK_a + beta*logK_b
            dlogK_a_b_dlogRho = alfa*dlogK_a_dlogRho + beta*dlogK_b_dlogRho
            dlogK_a_b_dlogT = alfa*dlogK_a_dlogT + beta*dlogK_b_dlogT
            
         end subroutine Get_Kap_at_dXO
         
         
         subroutine Get_Kap_at_dXC(dXC_lookup, &
                           dXC_a_lookup, dXO_a_lookup, logK_a, dlogK_a_dlogRho, dlogK_a_dlogT,  &
                           dXC_b_lookup, dXO_b_lookup, logK_b, dlogK_b_dlogRho, dlogK_b_dlogT,  &
                           logK_a_b, dlogK_a_b_dlogRho, dlogK_a_b_dlogT, dXO_a_b_lookup)
            real(dp), intent(in) :: dXC_lookup
            real(dp), intent(in) :: dXC_a_lookup, dXO_a_lookup
            real(dp), intent(in) :: logK_a, dlogK_a_dlogRho, dlogK_a_dlogT
            real(dp), intent(in) :: dXC_b_lookup, dXO_b_lookup
            real(dp), intent(in) :: logK_b, dlogK_b_dlogRho, dlogK_b_dlogT
            real(dp), intent(out) :: logK_a_b, dlogK_a_b_dlogRho, dlogK_a_b_dlogT
            real(dp), intent(out) :: dXO_a_b_lookup
            
            real(dp) :: alfa, beta
            
            if (dXC_a_lookup == dXC_b_lookup) then
               alfa = 0d0
            else
               alfa = (dXC_lookup - dXC_b_lookup) / (dXC_a_lookup - dXC_b_lookup)
            end if
            
            dXO_a_b_lookup = dXO_b_lookup + (dXO_a_lookup - dXO_b_lookup)*alfa
            beta = 1d0 - alfa
            
            logK_a_b = alfa*logK_a + beta*logK_b
            dlogK_a_b_dlogRho = alfa*dlogK_a_dlogRho + beta*dlogK_b_dlogRho
            dlogK_a_b_dlogT = alfa*dlogK_a_dlogT + beta*dlogK_b_dlogT
            
         end subroutine Get_Kap_at_dXC


      end subroutine Get_Kap_for_dXCO

         
      subroutine Find_CO_Tables( &
                  rq, x_tables, ix, CO_table_numbers, next_dXO_table, next_dXC_table,  &
                  co_tables, num_CO_tables, num_dXC_gt_dXO, &
                  dXCO_max, dXC, dXO, dXC_lookup, dXO_lookup, i1, i2, i3, i4, ierr)
      
         ! for linear interpolation to be smooth, 
         ! must use the smallest convex hull around the given point
         use kap_def
         use load_CO_kap
         
         type (Kap_General_Info), pointer :: rq
         type (Kap_CO_X_Table), dimension(:), pointer :: x_tables
         integer :: ix
         integer, intent(in) :: CO_table_numbers(num_kap_CO_dXs,num_kap_CO_dXs)
         integer, intent(in) :: next_dXO_table(max_num_CO_tables) 
         integer, intent(in) :: next_dXC_table(max_num_CO_tables) 
         type (Kap_CO_Table), dimension(:), pointer :: co_tables
         integer, intent(in) :: num_CO_tables, num_dXC_gt_dXO
         real(dp), intent(in) :: dXCO_max, dXC, dXO, dXC_lookup, dXO_lookup
         integer, intent(out) :: i1, i2, i3, i4
         integer, intent(out) :: ierr

         real(dp) :: dXC2_lookup, dXO2_lookup, dXC4_lookup, dXO4_lookup
         integer :: idXC, idXO
         real(dp), parameter :: tiny = 1d-7
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         if (dbg) write(*,*) 'enter Find_CO_Tables'
         if (dbg) write(*,*) 'associated(co_tables)', associated(co_tables)
         if (dbg) write(*,*) 'size(co_tables,dim=1)', size(co_tables,dim=1)
         if (dbg) write(*,2) 'num_kap_CO_dXs', num_kap_CO_dXs
         if (dbg) write(*,2) 'num_CO_tables', num_CO_tables

         ierr = 0
         
         ! find idXC s.t. kap_CO_dXs(idXC-1) < dXC <= kap_CO_dXs(idXC)
         do idXC = 2, num_kap_CO_dXs
            if (kap_CO_dXs(idXC) >= dXC) exit
         end do
         if (dbg) write(*,2) 'idXC', idXC
         
         ! find idXO s.t. kap_CO_dXs(idXO-1) < dXO <= kap_CO_dXs(idXO)
         do idXO = 2, num_kap_CO_dXs
            if (kap_CO_dXs(idXO) >= dXO) exit
         end do
         if (dbg) write(*,2) 'idXO', idXO
         
         i1 = CO_table_numbers(idXC-1,idXO-1)
         if (dbg) write(*,2) 'i1', i1
         if (matches_table(i1)) then
            i2 = -1; i3 = -1; i4 = -1
            if (dbg) write(*,2) 'matches_table(i1)', i1
            return
         end if
         
         if (dbg) write(*,*) '(dXC >= dXO)', (dXC >= dXO)
         if (dXC >= dXO) then
            i2 = CO_table_numbers(idXC,idXO-1)
            if (dbg) write(*,2) 'i2', i2
            i3 = CO_table_numbers(idXC-1,idXO)
            if (dbg) write(*,2) 'i3', i3
         else
            i2 = CO_table_numbers(idXC-1,idXO)
            if (dbg) write(*,2) 'i2', i2
            i3 = CO_table_numbers(idXC,idXO-1)
            if (dbg) write(*,2) 'i3', i3
         end if
         i4 = CO_table_numbers(idXC,idXO)
         if (dbg) write(*,2) 'i4', i4

         if (i4 > 0) then
            if (i1 <= 0 .or. i2 <= 0 .or. i3 <= 0) then
            
               write(*,2) 'i1', i1
               write(*,2) 'i2', i2
               write(*,2) 'i3', i3
               write(*,2) 'i4', i4
               write(*,2) 'idXC', idXC
               write(*,2) 'idXO', idXO
               
               write(*,1) 'dXCO_max', dble(dXCO_max)
               write(*,1) 'dXC', dble(dXC)
               write(*,1) 'dXO', dble(dXO)
               write(*,1) 'dXC_lookup', dble(dXC_lookup)
               write(*,1) 'dXO_lookup', dble(dXO_lookup)
               
               call mesa_error(__FILE__,__LINE__,'logical failure1 in looking for CO tables')
            end if
            if (matches_table(i2)) then
               i1 = i2; i2 = -1; i3 = -1; i4 = -1; return
            end if
            if (matches_table(i3)) then
               i1 = i3; i2 = -1; i3 = -1; i4 = -1; return
            end if
            if (matches_table(i4)) then
               i1 = i4; i2 = -1; i3 = -1; i4 = -1; return
            end if
            return
         end if
        
         if (on_midline(i1)) then ! middle triangle
            if (dXC >= dXO) then
               i2 = next_dXC_table(i1)
               i3 = next_dXO_table(i1)
            else
               i2 = next_dXO_table(i1)
               i3 = next_dXC_table(i1)
            end if
            return
         end if
         
         ! trapezoid or triangle near the diagonal boundary
            
         if (dXC >= dXO) then
         
            if (i3 <= 0) then
               if (on_diagonal(i1)) then
                  i3 = num_CO_tables
                  if (i2 > 0) then ! bail -- just use i1
                     i2 = -1; i3 = -1; i4 = -1; return
                  end if
               else
                  i3 = num_dXC_gt_dXO
               end if
            end if
            
            if (.not. on_diagonal(i3)) then
               i4 = next_dXC_table(i3)
               if (.not. on_diagonal(i4)) then ! bail -- just use i1
                  i2 = -1; i3 = -1; i4 = -1; return
               end if
            end if
            if (i2 <= 0) i2 = next_dXC_table(i1)
            if (on_diagonal(i2)) return
            
            
            dXC4_lookup = co_tables(i4)% dXC_lookup
            if (dXC_lookup <= dXC4_lookup) return
            ! check if on smaller dXC_lookup side of the line from i2 to i4
            dXC2_lookup = co_tables(i2)% dXC_lookup
            dXO2_lookup = co_tables(i2)% dXO_lookup
            dXO4_lookup = co_tables(i4)% dXO_lookup
            if (dXC_lookup <= dXC4_lookup+(dXC2_lookup-dXC4_lookup)* &
                     (dXO_lookup-dXO4_lookup)/(dXO2_lookup-dXO4_lookup)) return
            ! else we're in the dXC > dXO triangle
            i1 = i2; i2 = next_dXC_table(i1)
            if (.not. on_diagonal(i2)) then ! bail -- just use i1
               i2 = -1; i3 = -1; i4 = -1; return
            end if
            i3 = i4; i4 = -1
            
         else ! dXC < dXO
         
            ! reverse roles of dXC and dXO
         
            if (i3 <= 0) then ! must be in one of the triangles
               if (on_diagonal(i1)) then
                  i3 = num_dXC_gt_dXO
               else
                  i3 = num_CO_tables
                  if (i2 > 0) then ! bail -- just use i1
                     i2 = -1; i3 = -1; i4 = -1; return
                  end if
               end if
            end if
            if (.not. on_diagonal(i3)) then
               i4 = next_dXO_table(i3)
               if (.not. on_diagonal(i4)) then ! bail -- just use i1
                  i2 = -1; i3 = -1; i4 = -1; return
               end if
            end if
            if (i2 <= 0) i2 = next_dXO_table(i1)
            if (on_diagonal(i2)) return
            dXO4_lookup = co_tables(i4)% dXO_lookup
            if (dXO_lookup <= dXO4_lookup) return
            ! check if on smaller dXO_lookup side of the line from i2 to i4
            dXC2_lookup = co_tables(i2)% dXC_lookup
            dXO2_lookup = co_tables(i2)% dXO_lookup
            dXC4_lookup = co_tables(i4)% dXC_lookup
            if (dXO_lookup <= dXO4_lookup+(dXO2_lookup-dXO4_lookup)* &
                     (dXC_lookup-dXC4_lookup)/(dXC2_lookup-dXC4_lookup)) return
            ! else we're in the dXC < dXO triangle
            i1 = i2; i2 = next_dXO_table(i1)
            if (.not. on_diagonal(i2)) then ! bail -- just use i1
               i2 = -1; i3 = -1; i4 = -1; return
            end if
            i3 = num_CO_tables; i4 = -1
         
         end if
         
         
         contains
         
         
         logical function matches_table(i)
            integer :: i
            include 'formats'
            if (dbg) write(*,2) 'matches_table i', i
            if (abs(dXC_lookup - co_tables(i)% dXC_lookup) < tiny .and.  &
                     abs(dXO_lookup - co_tables(i)% dXO_lookup) < tiny) then
               matches_table = .true.
            else
               matches_table = .false.
            end if
            if (dbg) write(*,*) 'matches_table', matches_table
         end function matches_table
         
         
         logical function on_midline(i)
            integer :: i
            if (abs(co_tables(i)% dXC - co_tables(i)% dXO) < tiny) then
               on_midline = .true.
            else
               on_midline = .false.
            end if
         end function on_midline
         
         
         logical function on_diagonal(i)
            integer :: i
            if (abs((co_tables(i)% dXC + co_tables(i)% dXO) - dXCO_max) < tiny) then
               on_diagonal = .true.
            else
               on_diagonal = .false.
            end if
         end function on_diagonal
         

      end subroutine Find_CO_Tables

      
      subroutine Get_CO_Kap_for_logRho_logT( &
               rq, x_tables, ix, co_tables, ico, &
               logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         type (Kap_General_Info), pointer :: rq
         type (Kap_CO_X_Table), dimension(:), pointer :: x_tables
         integer :: ix
         type (Kap_CO_Table), dimension(:), pointer :: co_tables
         integer, intent(in) :: ico
         real(dp), intent(inout) :: logRho, logT
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr

         real(dp) :: logR0, logR1, logT0, logT1, logR, logR_in
         real(dp) :: df_dx, df_dy
         integer :: iR, jtemp, i, num_logRs, num_logTs
         logical :: clipped_logR

         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         if (dbg) write(*,1) 'enter Get_CO_Kap_for_logRho_logT', logRho, logT
         
         ierr = 0

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
            co_tables(ico)% kap1, num_logRs, num_logTs, &
            iR, jtemp, logR0, logR, logR1, logT0, logT, logT1, &
            logK, df_dx, df_dy)
         if (clipped_logR) df_dx = 0
         
         if (dbg) write(*,1) 'Do_Kap_Interpolations: logK', logK

         ! convert df_dx and df_dy to dlogK_dlogRho, dlogK_dlogT
         dlogK_dlogRho = df_dx
         dlogK_dlogT = df_dy - 3d0*df_dx

         if (dbg) write(*,*) 'done Get_CO_Kap_for_logRho_logT'

      end subroutine Get_CO_Kap_for_logRho_logT


      end module kap_eval_co
      
