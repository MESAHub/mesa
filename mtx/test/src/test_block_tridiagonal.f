! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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
      module test_block_tri_dble
#else
      module test_block_tri_quad
#endif
      
      use mtx_lib
      use mtx_def

#ifdef DBLE
      use const_def, only: dp
#else
      use const_def, only: qp
#endif

      use utils_lib, only: is_bad,mesa_error
      
      implicit none

#ifdef DBLE
      integer, parameter :: fltp = dp
#else
      integer, parameter :: fltp = qp
#endif
      
      integer, parameter :: caller_id = 0
      
      
      contains
      
      
#ifdef DBLE

      subroutine do_test_block_tri_dble
         call test_block(bcyclic_dble, .true.)
         call test_block(block_thomas_dble, .true.)
      end subroutine do_test_block_tri_dble

#else

      subroutine do_test_block_tri_quad
         call test_block(block_thomas_quad, .true.)
      end subroutine do_test_block_tri_quad

#endif
      
      
      subroutine test_block(which_decsol_option, for_release)

         integer, intent(in) :: which_decsol_option
         logical, intent(in) :: for_release
         
         integer :: nrep
         real(fltp), pointer :: lblk1(:), dblk1(:), ublk1(:) ! (nvar,nvar,nz)
         real(fltp), pointer :: lblk(:,:,:), dblk(:,:,:), ublk(:,:,:) ! (nvar,nvar,nz)
         real(fltp), pointer :: l_init(:,:,:), d_init(:,:,:), u_init(:,:,:) ! (nvar,nvar,nz)
         real(fltp), pointer :: x(:,:), xcorrect(:,:), brhs(:,:), work(:,:) ! (nvar,nz)
         real(fltp), pointer :: x1(:) ! =(nvar,nz)
         integer, pointer :: ipiv1(:) ! =(nvar,nz)
         integer, pointer :: ipiv(:,:) ! (nvar,nz)
          
         real(dp), pointer :: rpar_decsol(:) ! (lrd)
         integer, pointer :: ipar_decsol(:) ! (lid)
         real(fltp) :: time_factor, time_solve, time_refine, time_dealloc, sum_times
         real(fltp) :: sum_factor, setup, factor_As, &
            factor_As_thm_fac_avg, factor_As_thm_slv_avg, &
            factor_As_thm_fac_max, factor_As_thm_slv_max, &
            setup_coupling, factor_coupling

         integer :: i, j, k, ierr, lid, lrd, nvar, nz, rep
         logical :: use_given_weights
         character (len=255) :: fname, which_decsol_str
         
         include 'formats.dek'
         
         ierr = 0
         
         call decsol_option_str(which_decsol_option, which_decsol_str, ierr)
         if (ierr /= 0) return
         
         if (for_release) then
            fname = 'block_tri.data'
         else
            fname = 'block_tri_12.data'
         end if
         nrep = 1

         time_factor=0; time_solve=0; time_refine=0; time_dealloc=0
         
         call read_testfile(fname)
         !call xread_testfile(fname)

#ifdef DBLE
         if (which_decsol_option == bcyclic_dble) then
            write(*,*) 'bcyclic_dble'
            call bcyclic_dble_work_sizes(nvar,nz,lrd,lid)
         else
            write(*,*) 'bad value for which_decsol_option in test_block'
            call mesa_error(__FILE__,__LINE__)
         end if
#endif
         
         allocate( &
            rpar_decsol(lrd), ipar_decsol(lid), x1(nvar*nz), xcorrect(nvar,nz), &
            brhs(nvar,nz), ipiv1(nvar*nz), work(nvar,nz), &
            l_init(nvar,nvar,nz), d_init(nvar,nvar,nz), u_init(nvar,nvar,nz), stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in alloc'
            call mesa_error(__FILE__,__LINE__)
         end if
         ipiv(1:nvar,1:nz) => ipiv1(1:nvar*nz)
         x(1:nvar,1:nz) => x1(1:nvar*nz)
         
         
         do k=1,nz
            do i=1,nvar
               do j=1,nvar
                  l_init(j,i,k) = lblk(j,i,k)
                  d_init(j,i,k) = dblk(j,i,k)
                  u_init(j,i,k) = ublk(j,i,k)
               end do
            end do
         end do
         
         call set_xcorrect
         call set_brhs(lblk, dblk, ublk)
      
         if (.false.) then ! check data
            write(*,*) trim(fname)
            write(*,3) 'nvar nz', nvar, nz
            do k=1,nz
               do i=1,nvar
                  write(*,3) 'sol-xcorrect', i, k, brhs(i,k), xcorrect(i,k)
                  do j=1,nvar
                     write(*,4) 'u-d-l', i, j, k, &
                        ublk(i,j,k), dblk(i,j,k), lblk(i,j,k)
                  end do
               end do
            end do
            stop 'test_block_tridiagonal'
         end if

         do rep = 1, nrep
         
            do k=1,nz
               do i=1,nvar
                  do j=1,nvar
                     lblk(j,i,k) = l_init(j,i,k)
                     dblk(j,i,k) = d_init(j,i,k)
                     ublk(j,i,k) = u_init(j,i,k)
                  end do
               end do
            end do
            
            use_given_weights = (rep > 1)

#ifdef DBLE
            if (which_decsol_option == bcyclic_dble) then
               call solve_blocks( &
                  use_given_weights, bcyclic_dble_decsolblk, null_decsolblk_quad)
            else
               write(*,*) 'missing case for which_decsol_option', which_decsol_option
               call mesa_error(__FILE__,__LINE__)
            end if

#endif            
            
            
         end do
         
         call check_x

#ifdef DBLE
         if (which_decsol_option == bcyclic_dble) then
               write(*,*) 'done bcyclic_dble'
         else if (which_decsol_option == block_thomas_dble) then
               write(*,*) 'done block_thomas_dble'
         end if
#else
         if (which_decsol_option == block_thomas_quad) then
               write(*,*) 'done block_thomas_quad'
         end if
#endif
         
         write(*,*)

         deallocate(rpar_decsol, ipar_decsol, x1, xcorrect, work, &
            brhs, ipiv1, lblk1, dblk1, ublk1, l_init, d_init, u_init)
         
         contains
         
         
         subroutine solve_blocks(use_given_weights, decsolblk, decsolblk_quad)
            logical, intent(in) :: use_given_weights
            interface
               include 'mtx_decsolblk_dble.dek'
               include 'mtx_decsolblk_quad.dek'
            end interface
            
            integer :: iop, rep
            real(fltp) :: avg_err, max_err, err, atol, rtol
            integer :: i_max, j_max, i, j, k
         
            include 'formats.dek'         
                     
            iop = 0 ! factor A
#ifdef DBLE
            call decsolblk( &
               iop,caller_id,nvar,nz,lblk1,dblk1,ublk1,x1,ipiv1,lrd,rpar_decsol,lid,ipar_decsol,ierr)
#else
            call decsolblk_quad( &
               iop,caller_id,nvar,nz,lblk1,dblk1,ublk1,x1,ipiv1,lrd,rpar_decsol,lid,ipar_decsol,ierr)
#endif
            if (ierr /= 0) then
               write(*,*) 'decsolblk failed for factor'
               call mesa_error(__FILE__,__LINE__)
            end if
         
            do rep=1,1
            
               iop = 1 ! solve A*x = b
               
               do k=1,nz
                  do j=1,nvar
                     x(j,k) = brhs(j,k)
                  end do
               end do
#ifdef DBLE
               call decsolblk( &
                  iop,caller_id,nvar,nz,lblk1,dblk1,ublk1,x1,ipiv1,lrd,rpar_decsol,lid,ipar_decsol,ierr)
#else
               call decsolblk_quad( &
                  iop,caller_id,nvar,nz,lblk1,dblk1,ublk1,x1,ipiv1,lrd,rpar_decsol,lid,ipar_decsol,ierr)
#endif
               if (ierr /= 0) then
                  write(*,*) 'decsolblk failed for solve'
                  call mesa_error(__FILE__,__LINE__)
               end if
                        
            end do
         
            iop = 2 ! deallocate
#ifdef DBLE
            call decsolblk( &
               iop,caller_id,nvar,nz,lblk1,dblk1,ublk1,x1,ipiv1,lrd,rpar_decsol,lid,ipar_decsol,ierr)
#else
            call decsolblk_quad( &
               iop,caller_id,nvar,nz,lblk1,dblk1,ublk1,x1,ipiv1,lrd,rpar_decsol,lid,ipar_decsol,ierr)
#endif
            if (ierr /= 0) then
               write(*,*) 'decsolblk failed for deallocate'
               call mesa_error(__FILE__,__LINE__)
            end if
                     
         end subroutine solve_blocks
         
         
         subroutine read_testfile(fname)
            character (len=*), intent(in) :: fname
            integer :: iounit, ierr, i, j, k, line
            !write(*,*) 'reading ' // trim(fname)
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
         
            lblk(1:nvar,1:nvar,1:nz) => lblk1(1:nvar*nvar*nz)
            dblk(1:nvar,1:nvar,1:nz) => dblk1(1:nvar*nvar*nz)
            ublk(1:nvar,1:nvar,1:nz) => ublk1(1:nvar*nvar*nz)
            
            !return
            nz = 20
            write(*,*) 'testing with nvar,nz', nvar,nz
            
            
         end subroutine read_testfile
         
         
         subroutine set_brhs(lblk, dblk, ublk)
            real(fltp), pointer, dimension(:,:,:) :: lblk, dblk, ublk
            integer :: k, j
            include 'formats.dek'
            ! set brhs = A*xcorrect
#ifdef DBLE
            call block_dble_mv(nvar, nz, lblk, dblk, ublk, xcorrect, brhs)
#else
            call block_quad_mv(lblk, dblk, ublk, xcorrect, brhs)
#endif
            return
            do k = 1, 2 !nz
               do j = 1, nvar
                  if (brhs(j,k) /= 0) write(*,3) 'brhs xcorrect', j, k, brhs(j,k), xcorrect(j,k)
               end do
            end do
            write(*,*) 'end set_brhs'
            stop
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
               !write(*,3) 'BAD: avg err, max err, x, xcorrect', i, j, avg_err, max_err, x(i,j), xcorrect(i,j)
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
               end do
            end do
            avg_err = err_sum/(nz*nvar)
            !write(*,1) 'avg_err', avg_err
         end subroutine check1_x
         
         
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
           
      end subroutine test_block




#ifdef DBLE
      end module test_block_tri_dble
#else
      end module test_block_tri_quad
#endif
