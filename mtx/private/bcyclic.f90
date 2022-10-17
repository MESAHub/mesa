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

! derived from BCYCLIC as described in hirshman et. al.
! S.P.Hirshman, K.S.Perumalla, V.E.Lynch, & R.Sanchez,
! BCYCLIC: A parallel block tridiagonal matrix cyclic solver,
! J. Computational Physics, 229 (2010) 6392-6404.


module bcyclic
   
   use const_def, only : dp
   use my_lapack95_dble
   use utils_lib, only : set_nan, mesa_error
   
   implicit none
   
   type ulstore
      integer :: ul_size    ! size of umat1 & lmat1 (0 if not allocated)
      real(dp), pointer :: umat1(:), lmat1(:)
   end type ulstore
   
   type(ulstore), pointer :: odd_storage(:) => null()
   
   logical, parameter :: dbg = .false.
   
   logical, parameter :: do_fill_with_NaNs = .false.


contains
   
   
   subroutine bcyclic_factor (&
      lblk1, dblk1, ublk1, ipivot1, brhs1, nvar, nz, sparse, &
      lrd, rpar_decsol, lid, ipar_decsol, ierr)
      real(dp), pointer :: lblk1(:) ! row section of lower block
      real(dp), pointer :: dblk1(:) ! row section of diagonal block
      real(dp), pointer :: ublk1(:) ! row section of upper block
      integer, pointer :: ipivot1(:) ! row section of pivot array for block factorization
      real(dp), pointer :: brhs1(:) ! row section of rhs
      integer, intent(in) :: nvar ! linear size of each block
      integer, intent(in) :: nz ! number of block rows
      logical, intent(in) :: sparse
      integer, intent(in) :: lrd, lid
      real(dp), pointer, intent(inout) :: rpar_decsol(:) ! (lrd)
      integer, pointer, intent(inout) :: ipar_decsol(:) ! (lid)
      integer, intent(out) :: ierr
      
      integer, pointer :: iptr(:, :), nslevel(:), ipivot(:)
      integer :: ncycle, nstemp, maxlevels, nlevel, i, j, k
      logical :: have_odd_storage
      real(dp), pointer, dimension(:, :) :: dmat
      real(dp) :: dlamch, sfmin
      
      include 'formats'
      
      ierr = 0
      
      if (dbg) write(*, *) 'start bcyclic_factor'
      
      ! compute number of cyclic reduction levels
      ncycle = 1
      maxlevels = 0
      do while (ncycle < nz)
         ncycle = 2 * ncycle
         maxlevels = maxlevels + 1
      end do
      maxlevels = max(1, maxlevels)
      
      have_odd_storage = associated(odd_storage)
      if (have_odd_storage) then
         if (size(odd_storage) < maxlevels) then
            call clear_storage
            have_odd_storage = .false.
         end if
      end if
      
      if (.not. have_odd_storage) then
         allocate (odd_storage(maxlevels + 3), stat = ierr)
         if (ierr /= 0) then
            write(*, *) 'alloc failed for odd_storage in bcyclic'
            return
         end if
         do nlevel = 1, size(odd_storage)
            odd_storage(nlevel)% ul_size = 0
         end do
      end if
      
      allocate (nslevel(maxlevels), stat = ierr)
      if (ierr /= 0) return
      
      if (sparse) then
         write(*, *) 'no support for sparse matrix in bcyclic'
         ierr = -1
         return
      end if
      
      ncycle = 1
      nstemp = nz
      nlevel = 1
      
      if (dbg) write(*, *) 'start factor_cycle'
      
      factor_cycle : do ! perform cyclic-reduction factorization
         
         nslevel(nlevel) = nstemp
         
         if (dbg) write(*, 2) 'call cycle_onestep', nstemp
         
         call cycle_onestep(&
            nvar, nz, nstemp, ncycle, nlevel, sparse, &
            lblk1, dblk1, ublk1, ipivot1, ierr)
         if (ierr /= 0) then
            !write(*,*) 'cycle_onestep failed'
            call dealloc
            return
         end if
         
         if (nstemp == 1) exit
         
         nstemp = (nstemp + 1) / 2
         nlevel = nlevel + 1
         ncycle = 2 * ncycle
         
         if (nlevel > maxlevels) exit
      
      end do factor_cycle
      
      if (dbg) write(*, *) 'done factor_cycle'
      
      ! factor row 1
      dmat(1:nvar, 1:nvar) => dblk1(1:nvar * nvar)
      sfmin = dlamch('S')
      ipivot(1:nvar) => ipivot1(1:nvar)
      call my_getf2(nvar, dmat, nvar, ipivot, sfmin, ierr)
      if (ierr /= 0) then
         write(*, *) 'row 1 factor failed in bcyclic_factor'
         call dealloc
         return
      end if
      
      call dealloc
      
      if (dbg) write(*, *) 'done bcyclic_factor'
   
   contains
      
      subroutine dealloc
         deallocate (nslevel)
      end subroutine dealloc
   
   
   end subroutine bcyclic_factor
   
   
   subroutine bcyclic_solve (&
      lblk1, dblk1, ublk1, ipivot1, brhs1, nvar, nz, sparse, &
      lrd, rpar_decsol, lid, ipar_decsol, ierr)
      real(dp), pointer :: lblk1(:) ! row section of lower block
      real(dp), pointer :: dblk1(:) ! row section of diagonal block
      real(dp), pointer :: ublk1(:) ! row section of upper block
      integer, pointer :: ipivot1(:) ! row section of pivot array for block factorization
      real(dp), pointer :: brhs1(:)   ! row section of rhs
      integer, intent(in) :: nvar ! linear size of each block
      integer, intent(in) :: nz     ! number of block rows
      logical, intent(in) :: sparse
      integer, intent(in) :: lrd, lid
      real(dp), pointer, intent(inout) :: rpar_decsol(:) ! (lrd)
      integer, pointer, intent(inout) :: ipar_decsol(:) ! (lid)
      integer, intent(out) :: ierr
      
      integer, pointer :: iptr(:, :), nslevel(:), ipivot(:)
      integer :: ncycle, nstemp, maxlevels, nlevel, nvar2
      real(dp), pointer, dimension(:, :) :: dmat, bptr2
      
      include 'formats'
      
      if (dbg) write(*, *) 'start bcyclic_solve'
      
      ierr = 0
      nvar2 = nvar * nvar
      ncycle = 1
      maxlevels = 0
      do while (ncycle < nz)
         ncycle = 2 * ncycle
         maxlevels = maxlevels + 1
      end do
      maxlevels = max(1, maxlevels)
      
      allocate (nslevel(maxlevels), stat = ierr)
      if (ierr /= 0) return
      
      ncycle = 1
      nstemp = nz
      nlevel = 1
      
      if (dbg) write(*, *) 'start forward_cycle'
      
      forward_cycle : do
         
         nslevel(nlevel) = nstemp
         if (dbg) write(*, 2) 'call cycle_rhs', nstemp
         call cycle_rhs(&
            nstemp, nvar, ncycle, nlevel, sparse, dblk1, brhs1, ipivot1, ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if
         
         if (nstemp == 1) exit
         
         nstemp = (nstemp + 1) / 2
         nlevel = nlevel + 1
         ncycle = 2 * ncycle
         
         if (nlevel > maxlevels) exit
      
      end do forward_cycle
      
      if (dbg) write(*, *) 'done forward_cycle'
      
      ipivot(1:nvar) => ipivot1(1:nvar)
      dmat(1:nvar, 1:nvar) => dblk1(1:nvar2)
      bptr2(1:nvar, 1:1) => brhs1(1:nvar)
      call my_getrs(nvar, 1, dmat, nvar, ipivot, bptr2, nvar, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in bcyclic_solve'
         call dealloc
         return
      end if
      
      ! back solve for even x's
      back_cycle : do while (ncycle > 1)
         ncycle = ncycle / 2
         nlevel = nlevel - 1
         if (nlevel < 1) then
            ierr = -1
            exit
         end if
         nstemp = nslevel(nlevel)
         call cycle_solve(&
            nvar, nz, ncycle, nstemp, nlevel, sparse, lblk1, ublk1, brhs1)
      end do back_cycle
      
      call dealloc
      
      if (dbg) write(*, *) 'done bcyclic_solve'
   
   
   contains
      
      subroutine dealloc
         deallocate (nslevel)
      end subroutine dealloc
   
   
   end subroutine bcyclic_solve
   
   
   subroutine clear_storage
      integer :: nlevel
      nlevel = size(odd_storage)
      do while (nlevel > 0)
         if (odd_storage(nlevel)% ul_size > 0) then
            deallocate(odd_storage(nlevel)% umat1)
            deallocate(odd_storage(nlevel)% lmat1)
         end if
         nlevel = nlevel - 1
      end do
      deallocate(odd_storage)
      nullify(odd_storage)
   end subroutine clear_storage
   
   
   subroutine cycle_onestep(&
      nvar, nz, nblk, ncycle, nlevel, sparse, &
      lblk1, dblk1, ublk1, ipivot1, ierr)
      integer, intent(in) :: nvar, nz, nblk, ncycle, nlevel
      logical, intent(in) :: sparse
      real(dp), pointer, intent(inout) :: lblk1(:), dblk1(:), ublk1(:)
      integer, pointer, intent(inout) :: ipivot1(:)
      integer, intent(out) :: ierr
      
      integer, pointer :: ipivot(:)
      real(dp), pointer, dimension(:, :) :: dmat, umat, lmat, umat0, lmat0
      real(dp), pointer, dimension(:, :) :: lnext, unext, lprev, uprev
      real(dp), pointer :: mat1(:)
      integer :: i, j, shift, min_sz, new_sz, shift1, shift2, nvar2, &
         ns, ierr_loc, nmin, kcount, k, ii, jj, kk
      real(dp) :: dlamch, sfmin
      
      include 'formats'
      
      ierr = 0
      sfmin = dlamch('S')
      nvar2 = nvar * nvar
      nmin = 1
      kcount = 1 + (nblk - nmin) / 2
      min_sz = nvar2 * kcount
      if (odd_storage(nlevel)% ul_size < min_sz) then
         if (odd_storage(nlevel)% ul_size > 0) &
            deallocate(odd_storage(nlevel)% umat1, odd_storage(nlevel)% lmat1)
         new_sz = min_sz * 1.1 + 100
         odd_storage(nlevel)% ul_size = new_sz
         allocate (odd_storage(nlevel)% umat1(new_sz), &
            odd_storage(nlevel)% lmat1(new_sz), stat = ierr)
         if (ierr /= 0) then
            write(*, *) 'allocation error in cycle_onestep'
            return
         end if
      end if
      
      !$omp parallel do private(ns,kcount,shift,shift2,i)
      do ns = nmin, nblk, 2  ! copy umat and lmat
         kcount = (ns - nmin) / 2 + 1
         shift = nvar2 * (kcount - 1)
         shift2 = nvar2 * ncycle * (ns - 1)
         do i = 1, nvar2
            odd_storage(nlevel)% umat1(shift + i) = ublk1(shift2 + i)
            odd_storage(nlevel)% lmat1(shift + i) = lblk1(shift2 + i)
         end do
      end do
      !$omp end parallel do
      
      if (nvar2 * kcount > odd_storage(nlevel)% ul_size) then
         write(*, *) 'nvar2*kcount > ul_size in cycle_onestep'
         ierr = -1
         return
      end if
      
      if (dbg) write(*, *) 'start lu factorization'
      ! compute lu factorization of even diagonal blocks
      nmin = 2
      !$omp parallel do schedule(static,3) &
      !$omp private(ipivot,dmat,ns,ierr_loc,shift1,shift2,k)
      do ns = nmin, nblk, 2
         k = ncycle * (ns - 1) + 1
         shift1 = nvar * (k - 1)
         shift2 = nvar * shift1
         dmat(1:nvar, 1:nvar) => dblk1(shift2 + 1:shift2 + nvar2)
         ierr_loc = 0
         ipivot(1:nvar) => ipivot1(shift1 + 1:shift1 + nvar)
         call my_getf2(nvar, dmat, nvar, ipivot, sfmin, ierr_loc)
         if (ierr_loc /= 0) then
            ierr = ierr_loc
         end if
      end do
      !$omp end parallel do
      if (ierr /= 0) then
         !write(*,*) 'factorization failed in bcyclic'
         return
      end if
      
      if (dbg) write(*, *) 'done lu factorization; start solve'
      
      !$omp parallel do schedule(static,3) &
      !$omp private(ns,k,shift1,shift2,ipivot,dmat,umat,lmat,mat1,ierr_loc)
      do ns = nmin, nblk, 2
         ! compute new l=-d[-1]l, u=-d[-1]u for even blocks
         k = ncycle * (ns - 1) + 1
         shift1 = nvar * (k - 1)
         shift2 = nvar * shift1
         lmat(1:nvar, 1:nvar) => lblk1(shift2 + 1:shift2 + nvar2)
         ipivot(1:nvar) => ipivot1(shift1 + 1:shift1 + nvar)
         dmat(1:nvar, 1:nvar) => dblk1(shift2 + 1:shift2 + nvar2)
         call my_getrs(nvar, nvar, dmat, nvar, ipivot, lmat, nvar, ierr_loc)
         if (ierr_loc /= 0) ierr = ierr_loc
         lmat = -lmat
         umat(1:nvar, 1:nvar) => ublk1(shift2 + 1:shift2 + nvar2)
         call my_getrs(nvar, nvar, dmat, nvar, ipivot, umat, nvar, ierr_loc)
         if (ierr_loc /= 0) ierr = ierr_loc
         umat = -umat
      end do
      !$omp end parallel do
      if (dbg) write(*, *) 'done solve'
      
      if (ierr /= 0) return
      
      ! compute new odd blocks in terms of even block factors
      ! compute odd hatted matrix elements except at boundaries
      nmin = 1
      !$omp parallel do schedule(static,3) &
      !$omp private(i,ns,shift2,dmat,umat,lmat,lnext,unext,lprev,uprev,kcount,shift,umat0,lmat0,k)
      do i = 1, 3 * (1 + (nblk - nmin) / 2)
         
         ns = 2 * ((i - 1) / 3) + nmin
         k = ncycle * (ns - 1) + 1
         shift2 = nvar2 * (k - 1)
         dmat(1:nvar, 1:nvar) => dblk1(shift2 + 1:shift2 + nvar2)
         umat(1:nvar, 1:nvar) => ublk1(shift2 + 1:shift2 + nvar2)
         lmat(1:nvar, 1:nvar) => lblk1(shift2 + 1:shift2 + nvar2)
         
         if (ns < nblk) then
            shift2 = nvar2 * ncycle * ns
            lnext(1:nvar, 1:nvar) => lblk1(shift2 + 1:shift2 + nvar2)
            unext(1:nvar, 1:nvar) => ublk1(shift2 + 1:shift2 + nvar2)
         end if
         
         if (ns > 1) then
            shift2 = nvar2 * ncycle * (ns - 2)
            lprev(1:nvar, 1:nvar) => lblk1(shift2 + 1:shift2 + nvar2)
            uprev(1:nvar, 1:nvar) => ublk1(shift2 + 1:shift2 + nvar2)
         end if
         
         kcount = 1 + (ns - nmin) / 2
         shift = nvar2 * (kcount - 1)
         lmat0(1:nvar, 1:nvar) => odd_storage(nlevel)% lmat1(shift + 1:shift + nvar2)
         umat0(1:nvar, 1:nvar) => odd_storage(nlevel)% umat1(shift + 1:shift + nvar2)
         
         select case(mod(i - 1, 3))
         case (0)
            if (ns > 1) then
               ! lmat = matmul(lmat0, lprev)
               call my_gemm0_p1(nvar, nvar, nvar, lmat0, nvar, lprev, nvar, lmat, nvar)
            end if
         case (1)
            if (ns < nblk) then
               ! umat = matmul(umat0, unext)
               call my_gemm0_p1(nvar, nvar, nvar, umat0, nvar, unext, nvar, umat, nvar)
            end if
         case (2)
            if (ns < nblk) then
               if (ns > 1) then
                  ! dmat = dmat + matmul(umat0, lnext) + matmul(lmat0,uprev)
                  call my_gemm_plus_mm(nvar, nvar, nvar, umat0, lnext, lmat0, uprev, dmat)
               else
                  ! dmat = dmat + matmul(umat0, lnext)
                  call my_gemm_p1(nvar, nvar, nvar, umat0, nvar, lnext, nvar, dmat, nvar)
               end if
            else if (ns > 1) then
               ! dmat = dmat + matmul(lmat0,uprev)
               call my_gemm_p1(nvar, nvar, nvar, lmat0, nvar, uprev, nvar, dmat, nvar)
            end if
         end select
      
      end do
      !$omp end parallel do
      if (dbg) write(*, *) 'done cycle_onestep'
   
   end subroutine cycle_onestep
   
   
   subroutine cycle_rhs(&
      nblk, nvar, ncycle, nlevel, sparse, dblk1, brhs1, ipivot1, ierr)
      integer, intent(in) :: nblk, nvar, ncycle, nlevel
      logical, intent(in) :: sparse
      real(dp), pointer, intent(in) :: dblk1(:)
      real(dp), pointer, intent(inout) :: brhs1(:)
      integer, pointer, intent(in) :: ipivot1(:)
      integer, intent(out) :: ierr
      
      integer :: k, ns, ierr_loc, nmin, kcount, shift, i, shift1, shift2, nvar2
      integer, pointer :: ipivot(:)
      real(dp), pointer, dimension(:, :) :: dmat, umat, lmat, bptr2
      real(dp), pointer, dimension(:) :: bprev, bnext, bptr
      
      include 'formats'
      
      ierr = 0
      nvar2 = nvar * nvar
      ! compute dblk[-1]*brhs for even indices and store in brhs(even)
      nmin = 2
      ierr_loc = 0
      !$omp parallel do schedule(static,3) &
      !$omp private(ns,shift1,ipivot,shift2,k,dmat,bptr2,bptr,ierr_loc)
      do ns = nmin, nblk, 2
         k = ncycle * (ns - 1) + 1
         shift1 = nvar * (k - 1)
         shift2 = nvar * shift1
         ipivot(1:nvar) => ipivot1(shift1 + 1:shift1 + nvar)
         dmat(1:nvar, 1:nvar) => dblk1(shift2 + 1:shift2 + nvar2)
         bptr2(1:nvar, 1:1) => brhs1(shift1 + 1:shift1 + nvar)
         call my_getrs(nvar, 1, dmat, nvar, ipivot, bptr2, nvar, ierr_loc)
         if (ierr_loc /= 0) ierr = ierr_loc
      end do
      !$omp end parallel do
      if (ierr /= 0) return
      
      ! compute odd (hatted) sources (b-hats) for interior rows
      nmin = 1
      kcount = 0
      !$omp parallel do schedule(static,3) &
      !$omp private(ns,shift1,bptr,kcount,shift,umat,lmat,bnext,bprev)
      do ns = nmin, nblk, 2
         shift1 = nvar * ncycle * (ns - 1)
         bptr(1:nvar) => brhs1(shift1 + 1:shift1 + nvar)
         kcount = 1 + (ns - nmin) / 2
         shift = nvar2 * (kcount - 1)
         umat(1:nvar, 1:nvar) => odd_storage(nlevel)% umat1(shift + 1:shift + nvar2)
         lmat(1:nvar, 1:nvar) => odd_storage(nlevel)% lmat1(shift + 1:shift + nvar2)
         if (ns > 1) then
            shift1 = nvar * ncycle * (ns - 2)
            bprev => brhs1(shift1 + 1:shift1 + nvar)
         end if
         if (ns < nblk) then
            shift1 = nvar * ncycle * ns
            bnext => brhs1(shift1 + 1:shift1 + nvar)
            if (ns > 1) then
               ! bptr = bptr - matmul(umat,bnext) - matmul(lmat,bprev)
               call my_gemv_mv(nvar, nvar, umat, bnext, lmat, bprev, bptr)
            else
               ! bptr = bptr - matmul(umat,bnext)
               call my_gemv(nvar, nvar, umat, nvar, bnext, bptr)
            end if
         else if (ns > 1) then
            ! bptr = bptr - matmul(lmat,bprev)
            call my_gemv(nvar, nvar, lmat, nvar, bprev, bptr)
         end if
      end do
      !$omp end parallel do
      
      if (nvar2 * kcount > odd_storage(nlevel)% ul_size) then
         write(*, *) 'nvar2*kcount > ul_size in cycle_rhs'
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
   subroutine cycle_solve(&
      nvar, nz, ncycle, nblk, nlevel, sparse, lblk1, ublk1, brhs1)
      integer, intent(in) :: nvar, nz, ncycle, nblk, nlevel
      logical, intent(in) :: sparse
      real(dp), pointer, intent(in) :: lblk1(:), ublk1(:)
      real(dp), pointer, intent(inout) :: brhs1(:)
      
      real(dp), pointer :: umat(:, :), lmat(:, :), bprev(:), bnext(:), bptr(:)
      real(dp), pointer, dimension(:) :: bprevr, bnextr
      integer :: shift1, shift2, nvar2, ns, ierr, nmin
      
      nvar2 = nvar * nvar
      nmin = 2
      !$omp parallel do schedule(static,3) &
      !$omp private(ns,shift1,bptr,shift2,lmat,bprev,umat,bnext)
      do ns = nmin, nblk, 2
         shift1 = ncycle * nvar * (ns - 1)
         bptr(1:nvar) => brhs1(shift1 + 1:shift1 + nvar)
         shift2 = nvar * shift1
         lmat(1:nvar, 1:nvar) => lblk1(shift2 + 1:shift2 + nvar2)
         if (ns > 1) then
            shift1 = ncycle * nvar * (ns - 2)
            bprev(1:nvar) => brhs1(shift1 + 1:shift1 + nvar)
         end if
         if (ns < nblk) then
            umat(1:nvar, 1:nvar) => ublk1(shift2 + 1:shift2 + nvar2)
            shift1 = ncycle * nvar * ns
            bnext(1:nvar) => brhs1(shift1 + 1:shift1 + nvar)
            if (ns > 1) then
               ! bptr = bptr + matmul(umat,bnext) + matmul(lmat,bprev)
               call my_gemv_p_mv(nvar, nvar, umat, bnext, lmat, bprev, bptr)
            else
               ! bptr = bptr + matmul(umat,bnext)
               call my_gemv_p1(nvar, nvar, umat, nvar, bnext, bptr)
            end if
         else if (ns > 1) then
            ! bptr = bptr + matmul(lmat,bprev)
            call my_gemv_p1(nvar, nvar, lmat, nvar, bprev, bptr)
         end if
      end do
      !$omp end parallel do
   
   end subroutine cycle_solve
   
   
   subroutine bcyclic_deallocate (&
      lblk1, dblk1, ublk1, ipivot1, brhs1, nvar, nz, sparse, &
      lrd, rpar_decsol, lid, ipar_decsol, ierr)
      real(dp), pointer :: lblk1(:) ! row section of lower block
      real(dp), pointer :: dblk1(:) ! row section of diagonal block
      real(dp), pointer :: ublk1(:) ! row section of upper block
      integer, pointer :: ipivot1(:) ! row section of pivot array for block factorization
      real(dp), pointer :: brhs1(:) ! row section of rhs
      integer, intent(in) :: nvar ! linear size of each block
      integer, intent(in) :: nz ! number of block rows
      logical, intent(in) :: sparse
      integer, intent(in) :: lrd, lid
      real(dp), pointer, intent(inout) :: rpar_decsol(:) ! (lrd)
      integer, pointer, intent(inout) :: ipar_decsol(:) ! (lid)
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine bcyclic_deallocate


end module bcyclic
