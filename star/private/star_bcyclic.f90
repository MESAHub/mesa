! ***********************************************************************
! Copyright (C) 2012  The MESA Team
! This file is part of MESA.
! MESA is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Library Public License as published
! by the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
! MESA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General Public License
! along with this software; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
! ***********************************************************************

! derived from BCYCLIC written hirshman et. al.
! S.P.Hirshman, K.S.Perumalla, V.E.Lynch, & R.Sanchez,
! BCYCLIC: A parallel block tridiagonal matrix cyclic solver,
! J. Computational Physics, 229 (2010) 6392-6404.


      module star_bcyclic

      use star_private_def
      use const_def, only: dp, ln10
      use utils_lib, only: set_nan

      implicit none

      private
      public :: bcyclic_factor, bcyclic_solve, clear_storage

      logical, parameter :: dbg = .false.
      logical, parameter :: do_set_nan = .false.


      contains

      subroutine bcyclic_factor ( &
            s, nvar, nz, lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
            B1, row_scale_factors1, col_scale_factors1, &
            equed1, iter, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar ! linear size of each block
         integer, intent(in) :: nz ! number of block rows
         real(dp), pointer, dimension(:) :: &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, &
            B1, row_scale_factors1, col_scale_factors1
         integer, pointer :: ipivot1(:)
         character (len=nz) :: equed1
         integer, intent(in) :: iter ! solver iteration number for debugging output
         integer, intent(out) :: ierr

         integer, pointer :: iptr(:,:), nslevel(:), ipivot(:)
         integer :: neq, ncycle, nstemp, maxlevels, nlevel, i, j, k
         logical :: have_odd_storage
         real(dp), pointer, dimension(:,:) :: dmat, dmatF
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         real(dp) :: min_rcond_from_DGESVX, rpgfac
         integer :: k_min_rcond_from_DGESVX
         
         integer, allocatable :: factored(:)

         include 'formats'
         
         if (s% use_DGESVX_in_bcyclic .and. s% report_min_rcond_from_DGESXV) &
            min_rcond_from_DGESVX = 1d99
         
         allocate(factored(nz))
         do k=1,nz
            factored(k) = 0
         end do

         ierr = 0
         neq = nvar*nz
         !$omp simd
         do i = 1,nvar*neq
            lblkF1(i) = lblk1(i)
            dblkF1(i) = dblk1(i)
            ublkF1(i) = ublk1(i)
         end do

         if (dbg) write(*,*) 'start bcyclic_factor'

         ! compute number of cyclic reduction levels
         ncycle = 1
         maxlevels = 0
         do while (ncycle < nz)
            ncycle = 2*ncycle
            maxlevels = maxlevels+1
         end do
         maxlevels = max(1, maxlevels)

         have_odd_storage = associated(s% bcyclic_odd_storage)
         if (have_odd_storage) then
            if (size(s% bcyclic_odd_storage) < maxlevels) then
               call clear_storage(s)
               have_odd_storage = .false.
            end if
         end if

         if (.not. have_odd_storage) then
            allocate (s% bcyclic_odd_storage(maxlevels+3), stat=ierr)
            if (ierr /= 0) then
               write(*,*) 'alloc failed for odd_storage in bcyclic'
               return
            end if
            do nlevel = 1, size(s% bcyclic_odd_storage)
               s% bcyclic_odd_storage(nlevel)% ul_size = 0
            end do
         end if

         allocate (nslevel(maxlevels), stat=ierr)
         if (ierr /= 0) return

         ncycle = 1
         nstemp = nz
         nlevel = 1

         if (dbg) write(*,*) 'start factor_cycle'

         factor_cycle: do ! perform cyclic-reduction factorization

            nslevel(nlevel) = nstemp

            if (dbg) write(*,2) 'call cycle_onestep', nstemp

            call cycle_onestep( &
               s, nvar, nz, nstemp, ncycle, nlevel, iter, &
               lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
               B1, row_scale_factors1, col_scale_factors1, equed1, factored, &
               min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
               ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if

            if (nstemp == 1) exit

            nstemp = (nstemp+1)/2
            nlevel = nlevel+1
            ncycle = 2*ncycle

            if (nlevel > maxlevels) exit

         end do factor_cycle

         if (dbg) write(*,*) 'done factor_cycle'

         ! factor row 1
         dmat(1:nvar,1:nvar) => dblk1(1:nvar*nvar)
         dmatF(1:nvar,1:nvar) => dblkF1(1:nvar*nvar)
         ipivot(1:nvar) => ipivot1(1:nvar)
         row_scale_factors(1:nvar) => row_scale_factors1(1:nvar)
         col_scale_factors(1:nvar) => col_scale_factors1(1:nvar)
         factored(1) = factored(1) + 1
         call dense_factor(s, 1, nvar, dmat, dmatF, ipivot, &
            row_scale_factors, col_scale_factors, equed, &
            min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
            ierr)
         equed1(1:1) = equed(1:1)
         if (ierr /= 0) then
            write(*,*) 'dense_factor failed'
            call dealloc
            return
         end if
         
         do k=1,nz ! check that every cell factored exactly once
            if (factored(k) /= 1) then
               write(*,3) 'factored /= 1', k, factored(k)
               stop 'bcyclic_factor'
            end if
         end do

         call dealloc
            
         if (s% use_DGESVX_in_bcyclic .and. s% report_min_rcond_from_DGESXV) then
            write(*,4) 'DGESVX: k_min, iter, model, min rcond, rpgfac', &
               k_min_rcond_from_DGESVX, iter, s% model_number, min_rcond_from_DGESVX, rpgfac
         end if

         if (dbg) write(*,*) 'done bcyclic_factor'

         contains

         subroutine dealloc
            deallocate (nslevel)
         end subroutine dealloc


      end subroutine bcyclic_factor


      subroutine cycle_onestep( &
            s, nvar, nz, nblk, ncycle, nlevel, iter, &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
            B1, row_scale_factors1, col_scale_factors1, equed1, factored, &
            min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
            ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nz, nblk, ncycle, nlevel, iter
         real(dp), pointer, dimension(:), intent(inout) :: &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, &
            B1, row_scale_factors1, col_scale_factors1
         character (len=nz) :: equed1
         integer, pointer, intent(inout) :: ipivot1(:)
         integer :: factored(:)
         real(dp) :: min_rcond_from_DGESVX, rpgfac
         integer :: k_min_rcond_from_DGESVX
         integer, intent(out) :: ierr

         integer, pointer :: ipivot(:)
         real(dp), pointer, dimension(:,:) :: dmat, umat, lmat, umat0, lmat0, dmatF
         real(dp), pointer, dimension(:,:) :: lnext, unext, lprev, uprev
         real(dp), pointer, dimension(:) :: mat1
         integer :: i, j, shift, min_sz, new_sz, shift1, shift2, nvar2, &
            ns, op_err, nmin, kcount, k, ii, jj, kk
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed

         include 'formats'

         ierr = 0
         nvar2 = nvar*nvar
         nmin = 1
         kcount = 1+(nblk-nmin)/2
         min_sz = nvar2*kcount
         if (s% bcyclic_odd_storage(nlevel)% ul_size < min_sz) then
            if (s% bcyclic_odd_storage(nlevel)% ul_size > 0) &
               deallocate( &
                  s% bcyclic_odd_storage(nlevel)% umat1, &
                  s% bcyclic_odd_storage(nlevel)% lmat1)
            new_sz = min_sz*1.1 + 100
            s% bcyclic_odd_storage(nlevel)% ul_size = new_sz
            allocate (s% bcyclic_odd_storage(nlevel)% umat1(new_sz), &
                      s% bcyclic_odd_storage(nlevel)% lmat1(new_sz), stat=ierr)
            if (ierr /= 0) then
               write(*,*) 'allocation error in cycle_onestep'
               return
            end if
         end if

!$OMP PARALLEL DO private(ns,kcount,shift,shift2,i) SCHEDULE(static,3)
         do ns = nmin, nblk, 2  ! copy umat and lmat
            kcount = (ns-nmin)/2 + 1
            shift = nvar2*(kcount-1)
            shift2 = nvar2*ncycle*(ns-1)
            do i=1,nvar2
               s% bcyclic_odd_storage(nlevel)% umat1(shift+i) = ublkF1(shift2+i)
               s% bcyclic_odd_storage(nlevel)% lmat1(shift+i) = lblkF1(shift2+i)
            end do
         end do
!$OMP END PARALLEL DO

         if (nvar2*kcount > s% bcyclic_odd_storage(nlevel)% ul_size) then
            write(*,*) 'nvar2*kcount > ul_size in cycle_onestep'
            ierr = -1
            return
         end if

         if (dbg) write(*,*) 'start lu factorization'
         ! compute lu factorization of even diagonal blocks
         nmin = 2
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ipivot,dmat,dmatF,ns,op_err,shift1,shift2,i,j,k,row_scale_factors,col_scale_factors,equed)
         do ns = nmin, nblk, 2

            k = ncycle*(ns-1) + 1
            shift1 = nvar*(k-1)
            shift2 = nvar*shift1
            dmat(1:nvar,1:nvar) => dblk1(shift2+1:shift2+nvar2)
            dmatF(1:nvar,1:nvar) => dblkF1(shift2+1:shift2+nvar2)
            op_err = 0
            ipivot(1:nvar) => ipivot1(shift1+1:shift1+nvar)
            row_scale_factors(1:nvar) => row_scale_factors1(shift1+1:shift1+nvar)
            col_scale_factors(1:nvar) => col_scale_factors1(shift1+1:shift1+nvar)
            factored(k) = factored(k) + 1
            call dense_factor(s, k, nvar, dmat, dmatF, ipivot, &
               row_scale_factors, col_scale_factors, equed, &
               min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
               op_err)
            equed1(k:k) = equed(1:1)
            if (op_err /= 0) then
               ierr = op_err
            end if

         end do
!$OMP END PARALLEL DO
         if (ierr /= 0) then
            !write(*,*) 'factorization failed in bcyclic'
            return
         end if

         if (dbg) write(*,*) 'done lu factorization; start solve'

!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ns,k,shift1,shift2,ipivot,dmat,dmatF,umat,lmat,mat1,i,j,row_scale_factors,col_scale_factors,equed,op_err)
         do ns = nmin, nblk, 2
            ! compute new l=-d[-1]l, u=-d[-1]u for even blocks
            k = ncycle*(ns-1) + 1
            shift1 = nvar*(k-1)
            shift2 = nvar*shift1

            lmat(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)

            dmat(1:nvar,1:nvar) => dblk1(shift2+1:shift2+nvar2)
            dmatF(1:nvar,1:nvar) => dblkF1(shift2+1:shift2+nvar2)
            ipivot(1:nvar) => ipivot1(shift1+1:shift1+nvar)
            row_scale_factors(1:nvar) => row_scale_factors1(shift1+1:shift1+nvar)
            col_scale_factors(1:nvar) => col_scale_factors1(shift1+1:shift1+nvar)
            equed(1:1) = equed1(k:k)
            call dense_solve(s, k, nvar, dmat, dmatF, ipivot, lmat, &
               row_scale_factors, col_scale_factors, equed, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if

            do j=1,nvar
               do i=1,nvar
                  lmat(i,j) = -lmat(i,j)
               end do
            end do

            umat(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)

            call dense_solve(s, k, nvar, dmat, dmatF, ipivot, umat, &
               row_scale_factors, col_scale_factors, equed, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if

            do j=1,nvar
               do i=1,nvar
                  umat(i,j) = -umat(i,j)
               end do
            end do

         end do
!$OMP END PARALLEL DO
         if (dbg) write(*,*) 'done solve'

         if (ierr /= 0) return

         ! compute new odd blocks in terms of even block factors
         ! compute odd hatted matrix elements except at boundaries
         nmin = 1
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(i,ns,shift2,dmat,umat,lmat,lnext,unext,lprev,uprev,kcount,shift,umat0,lmat0,k)
         do i= 1, 3*(1+(nblk-nmin)/2)

            ns = 2*((i-1)/3) + nmin
            k = ncycle*(ns-1) + 1
            if (factored(k) > 0) then
               write(*,2) 'compute new dmat after already factored', k
               stop 'cycle_onestep'
            end if
            shift2 = nvar2*(k-1)
            dmat(1:nvar,1:nvar) => dblkF1(shift2+1:shift2+nvar2)
            umat(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)
            lmat(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)

            if (ns < nblk) then
               shift2 = nvar2*ncycle*ns
               lnext(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)
               unext(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)
            end if

            if (ns > 1) then
               shift2 = nvar2*ncycle*(ns-2)
               lprev(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)
               uprev(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)
            end if

            kcount = 1+(ns-nmin)/2
            shift = nvar2*(kcount-1)
            lmat0(1:nvar,1:nvar) => &
               s% bcyclic_odd_storage(nlevel)% lmat1(shift+1:shift+nvar2)
            umat0(1:nvar,1:nvar) => &
               s% bcyclic_odd_storage(nlevel)% umat1(shift+1:shift+nvar2)

            select case(mod(i-1,3))
            case (0)
               if (ns > 1) then
                  ! lmat = matmul(lmat0, lprev)
                  call my_gemm0_p1(nvar,nvar,nvar,lmat0,nvar,lprev,nvar,lmat,nvar)
               end if
            case (1)
               if (ns < nblk) then
                  ! umat = matmul(umat0, unext)
                  call my_gemm0_p1(nvar,nvar,nvar,umat0,nvar,unext,nvar,umat,nvar)
               end if
            case (2)
               if (ns < nblk) then
                  if (ns > 1) then
                     ! dmat = dmat + matmul(umat0, lnext) + matmul(lmat0,uprev)
                     call my_gemm_plus_mm(nvar,nvar,nvar,umat0,lnext,lmat0,uprev,dmat)
                  else
                     ! dmat = dmat + matmul(umat0, lnext)
                     call my_gemm_p1(nvar,nvar,nvar,umat0,nvar,lnext,nvar,dmat,nvar)
                  end if
               else if (ns > 1) then
                  ! dmat = dmat + matmul(lmat0,uprev)
                  call my_gemm_p1(nvar,nvar,nvar,lmat0,nvar,uprev,nvar,dmat,nvar)
               end if
            end select

         end do
!$OMP END PARALLEL DO
         if (dbg) write(*,*) 'done cycle_onestep'

      end subroutine cycle_onestep


      subroutine cycle_rhs( &
            s, nz, nblk, nvar, ncycle, nlevel, &
            dblk1, dblkF1, soln1, ipivot1, &
            row_scale_factors1, col_scale_factors1, equed1, ierr)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nblk, nvar, ncycle, nlevel
         real(dp), pointer, intent(in), dimension(:) :: &
            dblk1, dblkF1, row_scale_factors1, col_scale_factors1
         real(dp), pointer, intent(inout) :: soln1(:)
         integer, pointer, intent(in) :: ipivot1(:)
         character (len=nz) :: equed1
         integer, intent(out) :: ierr

         integer :: i, k, ns, op_err, nmin, kcount, shift, shift1, shift2, nvar2
         integer, pointer :: ipivot(:)
         real(dp), pointer, dimension(:,:) :: dmatF, dmat, umat, lmat
         real(dp), pointer, dimension(:) :: X, Xprev, Xnext
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         logical :: okay

         include 'formats'

         ierr = 0
         nvar2 = nvar*nvar
         ! compute dblk[-1]*brhs for even indices and store in brhs(even)
         nmin = 2
         op_err = 0
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ns,shift1,ipivot,shift2,k,dmat,dmatF,X,row_scale_factors,col_scale_factors,equed,i,okay,op_err)
         do ns = nmin, nblk, 2
            k = ncycle*(ns-1) + 1
            shift1 = nvar*(k-1)
            shift2 = nvar*shift1
            dmat(1:nvar,1:nvar) => dblk1(shift2+1:shift2+nvar2)
            dmatF(1:nvar,1:nvar) => dblkF1(shift2+1:shift2+nvar2)
            ipivot(1:nvar) => ipivot1(shift1+1:shift1+nvar)
            row_scale_factors(1:nvar) => row_scale_factors1(shift1+1:shift1+nvar)
            col_scale_factors(1:nvar) => col_scale_factors1(shift1+1:shift1+nvar)
            equed(1:1) = equed1(k:k)
            X(1:nvar) => soln1(shift1+1:shift1+nvar)
            call dense_solve1(s, k, nvar, X, dmat, dmatF, ipivot, .true., &
               row_scale_factors, col_scale_factors, equed, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if

         end do
!$OMP END PARALLEL DO

        if (ierr /= 0) return

        ! compute odd (hatted) sources (b-hats) for interior rows
         nmin = 1
         kcount = 0
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ns,shift1,X,kcount,shift,umat,lmat,Xnext,Xprev)
         do ns = nmin, nblk, 2
            shift1 = nvar*ncycle*(ns-1)
            X(1:nvar) => soln1(shift1+1:shift1+nvar)
            kcount = 1+(ns-nmin)/2
            shift = nvar2*(kcount-1)
            umat(1:nvar,1:nvar) => &
               s% bcyclic_odd_storage(nlevel)% umat1(shift+1:shift+nvar2)
            lmat(1:nvar,1:nvar) => &
               s% bcyclic_odd_storage(nlevel)% lmat1(shift+1:shift+nvar2)
            if (ns > 1) then
               shift1 = nvar*ncycle*(ns-2)
               Xprev => soln1(shift1+1:shift1+nvar)
            end if
            if (ns < nblk) then
               shift1 = nvar*ncycle*ns
               Xnext => soln1(shift1+1:shift1+nvar)
               if (ns > 1) then
                  ! bptr = bptr - matmul(umat,bnext) - matmul(lmat,bprev)
                  call my_gemv_mv(nvar,nvar,umat,Xnext,lmat,Xprev,X)
               else
                  ! bptr = bptr - matmul(umat,bnext)
                  call my_gemv(nvar,nvar,umat,nvar,Xnext,X)
               end if
            else if (ns > 1) then
               ! bptr = bptr - matmul(lmat,bprev)
               call my_gemv(nvar,nvar,lmat,nvar,Xprev,X)
            end if
         end do
!$OMP END PARALLEL DO

         if (nvar2*kcount > s% bcyclic_odd_storage(nlevel)% ul_size) then
            write(*,*) 'nvar2*kcount > ul_size in cycle_rhs'
            ierr = -1
            return
         end if

      end subroutine cycle_rhs


      ! computes even index solution from the computed (at previous,higher level)
      ! odd index solutions at this level.
      ! note at this point, the odd brhs values have been replaced (at the highest cycle)
      ! with the solution values (x), at subsequent (lower) cycles, the
      ! odd values are replaced by the even solutions at the next highest cycle. the even
      ! brhs values were multiplied by d[-1] and stored in cycle_rhs
      ! solve for even index values in terms of (computed at this point) odd index values
      subroutine cycle_solve( &
            s, nvar, nz, ncycle, nblk, nlevel, lblk1, ublk1, lblkF1, ublkF1, soln1)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nz, ncycle, nblk, nlevel
         real(dp), pointer, intent(in), dimension(:) :: lblk1, ublk1, lblkF1, ublkF1
         real(dp), pointer, intent(inout) :: soln1(:)

         real(dp), pointer :: umat(:,:), lmat(:,:), bprev(:), bnext(:), bptr(:)
         real(dp), pointer, dimension(:) :: bprevr, bnextr
         integer :: shift1, shift2, nvar2, ns, ierr, nmin, i, j

         include 'formats'

         nvar2 = nvar*nvar
         nmin = 2
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ns,shift1,bptr,shift2,lmat,bprev,umat,bnext)
         do ns = nmin, nblk, 2
            shift1 = ncycle*nvar*(ns-1)
            bptr(1:nvar) => soln1(shift1+1:shift1+nvar)
            shift2 = nvar*shift1
            lmat(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)
            if (ns > 1) then
               shift1 = ncycle*nvar*(ns-2)
               bprev(1:nvar) => soln1(shift1+1:shift1+nvar)
            end if
            if (ns < nblk) then
               umat(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)
               shift1 = ncycle*nvar*ns
               bnext(1:nvar) => soln1(shift1+1:shift1+nvar)
               if (ns > 1) then
                  ! bptr = bptr + matmul(umat,bnext) + matmul(lmat,bprev)
                  call my_gemv_p_mv(nvar,nvar,umat,bnext,lmat,bprev,bptr)
               else
                  ! bptr = bptr + matmul(umat,bnext)
                  call my_gemv_p1(nvar,nvar,umat,nvar,bnext,bptr)
               end if
            else if (ns > 1) then
               ! bptr = bptr + matmul(lmat,bprev)
               call my_gemv_p1(nvar,nvar,lmat,nvar,bprev,bptr)
            end if
         end do
!$OMP END PARALLEL DO

      end subroutine cycle_solve


      subroutine dense_factor(s, k, nvar, mtx, mtxF, ipivot, &
            row_scale_factors, col_scale_factors, equed, &
            min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
            ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer :: mtx(:,:), mtxF(:,:)
         integer, pointer :: ipivot(:)
         real(dp), pointer :: row_scale_factors(:), col_scale_factors(:)
         character (len=1) :: equed
         real(dp) :: min_rcond_from_DGESVX, rpgfac
         integer :: k_min_rcond_from_DGESVX
         integer, intent(out) :: ierr
         logical :: singular
         integer :: i, j
         real(dp), pointer :: work(:)
         integer, pointer :: iwork(:)
         real(dp) :: anorm, rcond
         include 'formats'
         ierr = 0
         
         if (s% use_DGESVX_in_bcyclic) then
            call factor_with_DGESVX
            return
         end if
         
         if (nvar == 4) then
            call my_getf2_n4(mtxF, ipivot, ierr)
         else if (nvar == 5) then
            call my_getf2_n5(mtxF, ipivot, ierr)
         else
            call my_getf2(nvar, mtxF, nvar, ipivot, ierr)
         end if
         
         contains
         
         subroutine factor_with_DGESVX
            character (len=1) :: fact, trans
            integer, parameter :: nrhs = 0
            real(dp) :: rcond
            real(dp) :: a(nvar,nvar), af(nvar,nvar), b(nvar,nrhs), x(nvar,nrhs), &
               r(nvar), c(nvar), ferr(nrhs), berr(nrhs), work(4*nvar)
            integer :: ipiv(nvar), iwork(nvar)
            integer :: i, j
            include 'formats'

            do i=1,nvar
               do j=1,nvar
                  a(i,j) = mtxF(i,j)
               end do
            end do
            
            if (s% use_equilibration_in_DGESVX) then
               fact = 'E' ! matrix A will be equilibrated, then copied to AF and factored
            else
               fact = 'N' ! matrix A will be copied to AF and factored
            end if
            trans = 'N' ! no transpose

!      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
!     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
!     $                   WORK, IWORK, INFO )

            call DGESVX(fact, trans, nvar, nrhs, a, nvar, af, nvar, ipiv, &
                        equed, r, c, b, nvar, x, nvar, rcond, ferr, berr, &
                        work, iwork, ierr)
               
            if (ierr > 0 .and. ierr <= nvar) then ! singular
               write(*,3) 'singular matrix for DGESVX', k, ierr
               stop 'factor_with_DGESVX'
            end if
            if (ierr == nvar+1) then ! means bad rcond, but may not be fatal
               write(*,2) 'DGESVX reports bad matrix conditioning: k, rcond', k, rcond
               ierr = 0
            end if
            
            do i=1,nvar
               do j=1,nvar
                  mtx(i,j) = a(i,j)
                  mtxF(i,j) = af(i,j)
               end do
               row_scale_factors(i) = r(i)
               col_scale_factors(i) = c(i)
               ipivot(i) = ipiv(i)
            end do
            
            if (s% report_min_rcond_from_DGESXV .and. rcond < min_rcond_from_DGESVX) then
               !$OMP critical (bcyclic_dense_factor_crit)
               min_rcond_from_DGESVX = rcond
               k_min_rcond_from_DGESVX = k
               rpgfac = work(1)
               !$OMP end critical (bcyclic_dense_factor_crit)
            end if

         end subroutine factor_with_DGESVX
         
      end subroutine dense_factor
   

      subroutine bcyclic_solve ( &
            s, nvar, nz, lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
            B1, soln1, row_scale_factors1, col_scale_factors1, equed1, &
            iter, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nz, iter
         real(dp), pointer, dimension(:) :: &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, &
            B1, soln1, row_scale_factors1, col_scale_factors1
         integer, pointer :: ipivot1(:)
         character (len=nz) :: equed1
         integer, intent(out) :: ierr

         integer, pointer :: iptr(:,:), nslevel(:), ipivot(:)
         integer :: ncycle, nstemp, maxlevels, nlevel, nvar2, i
         real(dp), pointer, dimension(:,:) :: dmat, dmatF
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         logical :: okay

         include 'formats'


         if (dbg) write(*,*) 'start bcyclic_solve'
         
         ! copy B to soln
         do i=1,nvar*nz
            soln1(i) = B1(i)
         end do

         ierr = 0

         nvar2 = nvar*nvar
         ncycle = 1
         maxlevels = 0
         do while (ncycle < nz)
            ncycle = 2*ncycle
            maxlevels = maxlevels+1
         end do
         maxlevels = max(1, maxlevels)

         allocate (nslevel(maxlevels), stat=ierr)
         if (ierr /= 0) return

         ncycle = 1
         nstemp = nz
         nlevel = 1

         if (dbg) write(*,*) 'start forward_cycle'

         forward_cycle: do

            nslevel(nlevel) = nstemp
            if (dbg) write(*,2) 'call cycle_rhs', nstemp
            call cycle_rhs( &
               s, nz, nstemp, nvar, ncycle, nlevel, &
               dblk1, dblkF1, soln1, ipivot1, &
               row_scale_factors1, col_scale_factors1, equed1, ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if

            if (nstemp == 1) exit

            nstemp = (nstemp+1)/2
            nlevel = nlevel+1
            ncycle = 2*ncycle

            if (nlevel > maxlevels) exit

         end do forward_cycle

         if (dbg) write(*,*) 'done forward_cycle'

         dmat(1:nvar,1:nvar) => dblk1(1:nvar2)
         dmatF(1:nvar,1:nvar) => dblkF1(1:nvar2)
         ipivot(1:nvar) => ipivot1(1:nvar)
         row_scale_factors(1:nvar) => row_scale_factors1(1:nvar)
         col_scale_factors(1:nvar) => col_scale_factors1(1:nvar)
         equed(1:1) = equed1(1:1)
         call dense_solve1(s, 1, nvar, soln1, dmat, dmatF, ipivot, .false., &
            row_scale_factors, col_scale_factors, equed, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in my_getrs1'
            call dealloc
            return
         end if

         ! back solve for even x's
         back_cycle: do while (ncycle > 1)
            ncycle = ncycle/2
            nlevel = nlevel-1
            if (nlevel < 1) then
               ierr = -1
               exit
            end if
            nstemp = nslevel(nlevel)
            call cycle_solve( &
               s, nvar, nz, ncycle, nstemp, nlevel, &
               lblk1, ublk1, lblkF1, ublkF1, soln1)
         end do back_cycle

         call dealloc

         if (dbg) write(*,*) 'done bcyclic_solve'


         contains


         subroutine dealloc
            deallocate (nslevel)
         end subroutine dealloc


      end subroutine bcyclic_solve


      subroutine clear_storage(s)
         type (star_info), pointer :: s
         integer :: nlevel
         nlevel = size(s% bcyclic_odd_storage)
         do while (nlevel > 0)
            if (s% bcyclic_odd_storage(nlevel)% ul_size > 0) then
               deallocate(s% bcyclic_odd_storage(nlevel)% umat1)
               deallocate(s% bcyclic_odd_storage(nlevel)% lmat1)
            end if
            nlevel = nlevel-1
         end do
         deallocate(s% bcyclic_odd_storage)
         nullify(s% bcyclic_odd_storage)
      end subroutine clear_storage


      subroutine dense_solve(s, k, nvar, mtx, mtxF, ipivot, X_mtx, &
            row_scale_factors, col_scale_factors, equed, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer, dimension(:,:) :: mtx, mtxF, X_mtx
         integer, pointer :: ipivot(:)
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         integer, intent(out) :: ierr
         integer :: i
         real(dp), pointer :: X(:)
         ierr = 0
         
         if (s% use_DGESVX_in_bcyclic) then
            call solve_with_DGESVX
            return
         end if

         do i=1,nvar
            X(1:nvar) => X_mtx(1:nvar,i)
            call dense_solve1(s, k, nvar, X, mtx, mtxF, ipivot, .false., &
               row_scale_factors, col_scale_factors, equed, ierr)
            if (ierr /= 0) return
         end do
         
         contains
         
         subroutine solve_with_DGESVX
            character (len=1) :: fact, trans
            real(dp) :: rcond
            real(dp) :: a(nvar,nvar), af(nvar,nvar), b(nvar,nvar), x(nvar,nvar), &
               r(nvar), c(nvar), ferr(nvar), berr(nvar), work(4*nvar)
            integer :: ipiv(nvar), iwork(nvar)
            integer :: i, j, nrhs
            include 'formats'

            nrhs = nvar

            do i=1,nvar
               do j=1,nvar
                  a(i,j) = mtx(i,j)
                  af(i,j) = mtxF(i,j)
                  b(i,j) = X_mtx(i,j)
                  x(i,j) = 0d0
               end do
               r(i) = row_scale_factors(i)
               c(i) = col_scale_factors(i)
               ipiv(i) = ipivot(i)
            end do
            
            fact = 'F' ! factored
            trans = 'N' ! no transpose

!      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
!     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
!     $                   WORK, IWORK, INFO )

            call DGESVX(fact, trans, nvar, nrhs, a, nvar, af, nvar, ipiv, &
                        equed, r, c, b, nvar, x, nvar, rcond, ferr, berr, &
                        work, iwork, ierr)
            if (ierr /= 0) then
               write(*,2) 'solve_with_DGESVX failed', k
            end if
            
            do i=1,nvar
               do j=1,nvar
                  X_mtx(i,j) = x(i,j)
               end do
            end do

         end subroutine solve_with_DGESVX

      end subroutine dense_solve


      subroutine dense_solve1(s, k, nvar, X_vec, mtx, mtxF, ipivot, dbg, &
            row_scale_factors, col_scale_factors, equed, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer :: X_vec(:), mtx(:,:), mtxF(:,:)
         integer, pointer :: ipivot(:)
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         logical, intent(in) :: dbg
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         
         if (s% use_DGESVX_in_bcyclic) then
            call solve1_with_DGESVX
            return
         end if
         
         if (nvar == 4) then
            call my_getrs1_n4(mtxF, ipivot, X_vec, ierr)
         else if (nvar == 5) then
            call my_getrs1_n5(mtxF, ipivot, X_vec, ierr)
         else
            call my_getrs1(nvar, mtxF, nvar, ipivot, X_vec, nvar, ierr)
         end if
         
         contains
         
         subroutine solve1_with_DGESVX
            character (len=1) :: fact, trans
            real(dp) :: rcond
            integer, parameter :: nrhs = 1
            real(dp) :: a(nvar,nvar), af(nvar,nvar), b(nvar,nrhs), x(nvar,nrhs), &
               r(nvar), c(nvar), ferr(nrhs), berr(nrhs), work(4*nvar)
            integer :: ipiv(nvar), iwork(nvar)
            integer :: i, j

            include 'formats'

            do i=1,nvar
               do j=1,nvar
                  a(i,j) = mtx(i,j)
                  af(i,j) = mtxF(i,j)
               end do
               b(i,1) = X_vec(i)
               x(i,1) = 0d0
               r(i) = row_scale_factors(i)
               c(i) = col_scale_factors(i)
               ipiv(i) = ipivot(i)
            end do
            
            fact = 'F' ! factored
            trans = 'N' ! no transpose

!      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
!     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
!     $                   WORK, IWORK, INFO )

            call DGESVX(fact, trans, nvar, nrhs, a, nvar, af, nvar, ipiv, &
                        equed, r, c, b, nvar, x, nvar, rcond, ferr, berr, &
                        work, iwork, ierr)
            
            do i=1,nvar
               X_vec(i) = x(i,1)
            end do

         end subroutine solve1_with_DGESVX

      end subroutine dense_solve1



      subroutine bcyclic_deallocate (s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine bcyclic_deallocate


      include 'mtx_solve_routines.inc'


      end module star_bcyclic
