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

      module freedman
      
      use interp_1d_lib_sg
      use interp_1d_def
      use math_lib
      use utils_lib, only: mesa_error
      use const_def, only: dp

      implicit none
      
      
      integer, parameter :: npoints = 736, num_Ts = 42, max_num_Rhos = 18
      real :: logTs(num_Ts)
      real, dimension(max_num_Rhos,num_Ts) :: logRhos, logKappas
      integer :: num_logRhos_for_logT(num_Ts)
      real, pointer :: f1(:), f(:,:,:)
      real, target :: f_ary(4*max_num_Rhos*num_Ts)
      contains

      
      
      
      subroutine get_Freedman_fname(data_dir, Z, fname)
         real(dp), intent(in) :: Z
         character (len=*),intent(in) :: data_dir
         character (len=*),intent(out) :: fname
         integer :: iz
         iz = floor(Z*1d5 + 0.1d0)
         select case (iz)
         case (1000)
            fname = trim(data_dir) // '/m0.3.txt'
         case (2000)
            fname = trim(data_dir) // '/p0.0.txt'
         case (4000)
            fname = trim(data_dir) // '/p0.3.txt'
         case (10000)
            fname = trim(data_dir) // '/p0.7.txt'
         case (20000)
            fname = trim(data_dir) // '/p1.0.txt'
         case (63000)
            fname = trim(data_dir) // '/p1.5.txt'
         case (100000)
            fname = trim(data_dir) // '/p1.7.txt'
         case default
            write(*,*) 'get_Freedman_fname: unexpected Z value for Freedman data', Z
            call mesa_error(__FILE__,__LINE__)
         end select
      end subroutine get_Freedman_fname


      
      subroutine init_freedman(freedman_data_dir, Z)
         character (len=*),intent(in) :: freedman_data_dir
         real(dp) :: Z
      
         integer :: i, j, k, ierr, io_logK, ii     
         real :: T, P, Rho, kap, logT, logT_prev, logRho, logKappa
         real, pointer :: work1(:)
         real, target :: work_ary(max_num_Rhos*pm_work_size)
         character (len=6) :: str
         character (len=256) :: fname

         include 'formats'

         work1 => work_ary

         ierr = 0
         f(1:4,1:max_num_Rhos,1:num_Ts) => &
             f_ary(1:4*max_num_Rhos*num_Ts)

         call get_Freedman_fname(freedman_data_dir, Z, fname)
         !write(*,*) 'init_freedman: read ' // trim(fname)

         io_logK = 40
         open(unit=io_logK, file=trim(fname), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'init_freedman failed to open ', trim(fname)
            call mesa_error(__FILE__,__LINE__)
         end if
         
         ! skip 2 lines
         do i = 1, 2
            read(io_logK, *, iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'init_freedman failed while reading header ', trim(fname)
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         
         logT_prev = -1
         j = 1
         k = 1
            
         do i = 1, npoints
            read(io_logK, *, iostat=ierr) ii, T, P, Rho, kap
            if (ierr /= 0) then
               write(*,*) 'init_freedman failed while reading ', trim(fname)
               write(*,*) 'i', i
               call mesa_error(__FILE__,__LINE__)
            end if
            logT = log10(T)
            logRho = log10(Rho)
            logKappa = log10(kap)
            !if (logT_prev > 0) write(*,3) 'logT prev', i, j, logT, logT_prev
            if (i == 1) then
               logT_prev = logT
            else if (abs(logT - logT_prev) > 1e-6) then
               num_logRhos_for_logT(j) = k
               j = j + 1
               if (j > num_Ts) then
                  write(*,3) 'init_freedman: too many logTs', i, j, logT
                  call mesa_error(__FILE__,__LINE__)
               end if
               k = 1
               logT_prev = logT
            else
               k = k+1
               if (k > max_num_Rhos) then
                  write(*,2) 'init_freedman: too many rho values: T, Rho', i, T, Rho
                  call mesa_error(__FILE__,__LINE__)
               end if
            end if
            logTs(j) = logT
            logRhos(k,j) = logRho
            logKappas(k,j) = logKappa
            if (i == npoints) then
               num_logRhos_for_logT(j) = k
               if (j /= num_Ts) then
                  write(*,2) 'init_freedman: wrong number of logTs', j
                  call mesa_error(__FILE__,__LINE__)
               end if
            end if
            !write(*,2) 'T Rho kap', i, T, Rho, kap
         end do
         
         close(io_logK)
         
         if (sum(num_logRhos_for_logT) /= npoints) then
            write(*,3) 'init_freedman: bad sum for num logRhos', sum(num_logRhos_for_logT), npoints
            call mesa_error(__FILE__,__LINE__)
         end if

         do j = 1, num_Ts
            f(1,:,j) = logKappas(:,j)
            f1(1:4*max_num_Rhos) => f_ary(1+(j-1)*max_num_Rhos*4:j*max_num_Rhos*4)
            call interp_pm_sg( &
               logRhos(:,j), num_logRhos_for_logT(j), f1, pm_work_size, work1, 'Freedman', ierr)
            if (ierr /= 0) then
               write(*,2) 'init_freedman: failed in interp_pm', j, logTs(j)
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         
       end subroutine init_freedman

      
      subroutine eval_freedman (T,rho,logKap,dbg_in)
         real(dp), intent(in) :: T,rho
         real(dp), intent(out) :: logKap
         logical, intent(in) :: dbg_in
         
         real :: logT, logRho
         integer :: ierr, j, i, j_logT
         integer, parameter :: n_old=4, n_new=1
         real :: x_old(n_old), v_old(n_old), x_new(n_new), v_new(n_new)
         real, pointer :: work1(:)
         real, target :: work_ary(n_old*pm_work_size)
         logical, parameter :: dbg = .false.

         include 'formats'
         
         work1 => work_ary

         logT = max(logTs(1),min(logTs(num_Ts),real(log10(T),kind=kind(logT))))
         logRho = log10(rho)
         ierr = 0
                  
         do j = 1, num_Ts-1
         
            if (logT >= logTs(j) .and. logT <= logTs(j+1)) then
            
               j_logT = min(max(1,j-1),num_Ts-3)
               do i = 1, n_old
                  ! interpolate at logRho for logTs(j_logT)
                  x_new(1) = logRho
                  f1(1:4*max_num_Rhos) => &
                     f_ary(1+(j_logT-1)*max_num_Rhos*4:j_logT*max_num_Rhos*4)
                  call interp_values_sg( &
                     logRhos(:,j_logT), num_logRhos_for_logT(j_logT), &
                     f1, n_new, x_new, v_new, ierr)
                  if (ierr /= 0) then
                     write(*,3) 'eval_freedman: interp failed for logRho', &
                        j, j_logT, logTs(j), logRho
                     call mesa_error(__FILE__,__LINE__)
                  end if
                  x_old(i) = logTs(j_logT)
                  v_old(i) = v_new(1)
                  
                  if (dbg) write(*,2) 'interp: logT logRho logKap', j_logT, logTs(j_logT), logRho, v_new(1)

                  j_logT = j_logT + 1

               end do       
                       
               ! interpolate logKap at logT    
               x_new(1) = logT
               call interpolate_vector_sg( &
                  n_old, x_old, n_new, x_new, v_old, v_new, &
                  interp_pm_sg, pm_work_size, work1, 'Freedman', ierr)
               if (ierr /= 0) then
                  write(*,2) 'eval_freedman: interpolate_vector_sg', &
                     j, logT, logRho
                  do i=1,n_old
                     write(*,2) 'logT', i, x_old(i)
                  end do
                  call mesa_error(__FILE__,__LINE__)
               end if
               logKap = dble(v_new(1))
               
               if (dbg) write(*,2) 'final: logT logRho logKap', 0, logT, logRho, logKap
               if (dbg) stop
               return
               
            end if
            
            if (dbg) write(*,2) 'skip over logT', j, logTs(j), logT
            
         end do

         write(*,1) 'eval_freedman: logT confusion', logT
         call mesa_error(__FILE__,__LINE__)


      end subroutine eval_freedman
      

      end module freedman
