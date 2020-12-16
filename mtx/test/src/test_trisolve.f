! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not,write to the Free Software
!   Foundation,Inc.,59 Temple Place,Suite 330,Boston,MA 02111-1307 USA
!
! ***********************************************************************

#ifdef DBLE
   module test_trisolve_dble
#else
   module test_trisolve_quad
#endif
   
   use mtx_lib
   use mtx_def
   
#ifdef DBLE
      use const_def, only: dp
#else
      use const_def, only: qp
#endif

      use utils_lib, only: is_bad, mesa_error

      implicit none

#ifdef DBLE
      integer, parameter :: fltp = dp
#else
      integer, parameter :: fltp = qp
#endif

   
   contains
   

#ifdef DBLE
   subroutine do_test_trisolve_dble(do_timing, for_release)
#else
   subroutine do_test_trisolve_quad(do_timing, for_release)
#endif
      logical, intent(in) :: do_timing, for_release
      integer :: nvar, nz, ierr, i, j, k, lwork, liwork, n
      real(fltp), dimension(:), pointer :: ublk_tst1, dblk_tst1, lblk_tst1, sol_tst1
      real(fltp), dimension(:,:,:), pointer :: ublk_tst, dblk_tst, lblk_tst
      real(fltp), dimension(:,:), pointer :: sol_tst, x, xcorrect
      real(dp), dimension(:), pointer :: work
      integer, dimension(:), pointer :: iwork
      character (len=255) :: fname
      
      include 'formats.dek'
      
      ierr = 0
      write(*,*) 'do_test_trisolve'
      ! setup test problem
      if (for_release) then
         fname = 'block_tri.data'
      else
         fname = 'block_tri_12.data'
      end if
      call read_testfile(fname, lblk_tst1, dblk_tst1, ublk_tst1)
      lblk_tst(1:nvar,1:nvar,1:nz) => lblk_tst1(1:nvar*nvar*nz)
      dblk_tst(1:nvar,1:nvar,1:nz) => dblk_tst1(1:nvar*nvar*nz)
      ublk_tst(1:nvar,1:nvar,1:nz) => ublk_tst1(1:nvar*nvar*nz)
      
#ifdef DBLE
      call mtx_trisolve_dble_work_sizes(nvar,nz,lwork,liwork)
#else
      call mtx_trisolve_quad_work_sizes(nvar,nz,lwork,liwork)
#endif
      allocate( &
         x(nvar,nz), xcorrect(nvar,nz), sol_tst1(nvar*nz), &
         work(lwork), iwork(liwork))
      
      sol_tst(1:nvar,1:nz) => sol_tst1(1:nvar*nz)

      call set_xcorrect
      call set_brhs(lblk_tst, dblk_tst, ublk_tst, sol_tst)
      
      
      if (.false.) then ! output data in simple format
         open(12,file='test_case.data',status='unknown')
         write(12,*) nvar, nz
         do k=1,nz
            do i=1,nvar
               do j=1,nvar
                  write(12,*) ublk_tst(i,j,k), dblk_tst(i,j,k), lblk_tst(i,j,k)
               end do
            end do
            do i=1,nvar
               write(12,*) sol_tst(i,k), xcorrect(i,k)
            end do
         end do
         close(12)
         stop 'done output to test_case.data'
      end if
      
      
      if (.false.) then ! check data
         write(*,*) trim(fname)
         write(*,3) 'nvar nz', nvar, nz
         do k=1,nz
            do i=1,nvar
               write(*,3) 'sol-xcorrect', i, k, sol_tst(i,k), xcorrect(i,k)
               do j=1,nvar
                  write(*,4) 'u-d-l', i, j, k, &
                     ublk_tst(i,j,k), dblk_tst(i,j,k), lblk_tst(i,j,k)
               end do
            end do
         end do
         stop 'trisolve'
      end if
      
!      sol_tst(1:nvar,1:nz) = 1d0
!      ublk_tst(1:nvar,1:nvar,1:nz) = 1d0
!      dblk_tst(1:nvar,1:nvar,1:nz) = 1d0
!      lblk_tst(1:nvar,1:nvar,1:nz) = 1d0
!      do i=1,nvar
!         do k=1,nz
!            dblk_tst(i,i,k) = 4d0
!         end do
!      end do
      
      
#ifdef DBLE
      call mtx_trisolve_dble( &
#else
      call mtx_trisolve_quad( &
#endif
         nvar, nz, lblk_tst1, dblk_tst1, ublk_tst1, sol_tst1, &
         lwork, work, liwork, iwork, ierr)
      if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
      
      x = sol_tst
      call check_x

!      open(12,file='solution.dat',status='unknown')
!      do n=1,nz
!         do i=1,nvar
!            write(12,*) sol_tst(i,n)
!            write(*,*) sol_tst(i,n)
!         end do
!      end do
      write(*,*)

      deallocate( &
         lblk_tst1, dblk_tst1, ublk_tst1, sol_tst1, x, xcorrect, work, iwork)

      !stop 'trisolve'
      
      
      contains
         
         
      subroutine read_testfile(fname,lblk1,dblk1,ublk1)
         character (len=*), intent(in) :: fname
         real(fltp), dimension(:), pointer :: lblk1,dblk1,ublk1
         integer :: iounit, ierr, i, j, k, line
         iounit = 33; ierr = 0
         open(unit=iounit, file=trim(fname), status='old', action='read', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname)
            call mesa_error(__FILE__,__LINE__)
         end if
#ifdef DBLE
         call mtx_read_block_tridiagonal(iounit,nvar,nz,lblk1,dblk1,ublk1,ierr)
#else
         call mtx_read_quad_block_tridiagonal(iounit,nvar,nz,lblk1,dblk1,ublk1,ierr)
#endif
         if (ierr /= 0) then
            write(*,*) 'failed to read ' // trim(fname)
            call mesa_error(__FILE__,__LINE__)
         end if
         close(iounit)
         
         return
         nz = 9
         write(*,*) 'testing with small nz', nz
         
         
      end subroutine read_testfile
      
      
      subroutine set_xcorrect
         real(fltp) :: cnt
         integer :: k, j
         cnt = 1d0
         do k=1,nz
            do j=1,nvar
               cnt = cnt + 1d-3
               xcorrect(j,k) = cnt
            end do
         end do
      end subroutine set_xcorrect
         
         
      subroutine set_brhs(lblk, dblk, ublk, brhs)
         real(fltp), pointer, dimension(:,:,:) :: lblk, dblk, ublk
         real(fltp), pointer, dimension(:,:) :: brhs
         integer :: k, j
         include 'formats.dek'
         ! set brhs = A*xcorrect
#ifdef DBLE
         call block_dble_mv(lblk, dblk, ublk, xcorrect, brhs)
#else
         call block_quad_mv(lblk, dblk, ublk, xcorrect, brhs)
#endif
      end subroutine set_brhs
      
      
      subroutine check_x
         real(fltp) :: max_err, err, atol, rtol, avg_err
         integer :: i_max, j_max,i, j, rep        
         include 'formats.dek'
         atol = 1d-4
         rtol = 1d-4   
         call check1_x(avg_err, max_err, atol, rtol, i_max, j_max)
         i = i_max; j = j_max
         if (max_err > 1) then
            write(*,3) 'BAD: err, x, xcorrect', i, j, max_err, x(i,j), xcorrect(i,j)
         else
            write(*,3) 'solution matches xcorrect: nvar, nz', nvar, nz
         end if
      end subroutine check_x
         
         
      subroutine check1_x(avg_err, max_err, atol, rtol, i_max, j_max)
         real(fltp), intent(out) :: avg_err, max_err
         real(fltp), intent(in) ::  atol, rtol
         integer, intent(out) :: i_max, j_max
         integer :: i, j
         real(fltp) :: err_sum
         real(fltp) :: err
         include 'formats.dek'      
         max_err = 0; i_max = 0; j_max = 0; err_sum = 0
         do j = 1, nz
            do i = 1, nvar
               if (is_bad(x(i,j))) then
                  write(*,3) 'x xcorrect', i, j, x(i,j), xcorrect(i,j)
                  stop 'check1_x'
               end if
               err = abs(x(i,j) - xcorrect(i,j))/(atol + rtol*max(abs(x(i,j)),abs(xcorrect(i,j))))
               err_sum = err_sum + err
               if (err > max_err) then
                  max_err = err; i_max = i; j_max = j
               end if
               !write(*,3) 'err', i, j, err, x(i,j), xcorrect(i,j)
            end do
         end do
         avg_err = err_sum/(nz*nvar)
      end subroutine check1_x
      
      
      

#ifdef DBLE
   end subroutine do_test_trisolve_dble
#else
   end subroutine do_test_trisolve_quad
#endif


#ifdef DBLE
   end module test_trisolve_dble
#else
   end module test_trisolve_quad
#endif
