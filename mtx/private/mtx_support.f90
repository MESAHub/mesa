! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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


      module mtx_support
      
      use const_def, only: dp, qp
      use utils_lib, only: mesa_error
      
      integer, parameter :: num_chunks = 4

      contains


      
      subroutine do_dense_to_band(n,ndim,a,ml,mu,ab,ldab,ierr)
         integer, intent(in) :: n,ndim,ml,mu,ldab
         real(dp), intent(in) :: a(:,:) ! (ndim,n)
         real(dp), intent(inout) :: ab(:,:) ! (ldab,n)
         integer, intent(out) :: ierr 
         integer :: i, j
         if (ml+mu+1 > n) then
            ierr = -1
            return
         end if
         ierr = 0
         ab = 0
         do j=1,n
            do i=max(1,j-mu),min(n,j+ml)
               ab(ldab-ml+i-j,j) = a(i,j)
            end do
         end do
      end subroutine do_dense_to_band

      
      subroutine do_band_to_dense(n,ml,mu,ab,ldab,ndim,a,ierr)
         integer, intent(in) :: n,ndim,ml,mu,ldab
         real(dp), intent(in) :: ab(:,:) ! (ldab,n)
         real(dp), intent(inout) :: a(:,:) ! (ndim,n)
         integer, intent(out) :: ierr 
         integer :: i, j
         if (ml+mu+1 > n) then
            ierr = -1
            return
         end if
         ierr = 0
         a = 0
         do j=1,n
            do i=max(1,j-mu),min(n,j+ml)
               a(i,j) = ab(ldab-ml+i-j,j)
            end do
         end do
      end subroutine do_band_to_dense

      
      subroutine do_band_to_column_sparse(n,ml,mu,ab,ldab,nzmax,nz,colptr,rowind,values,diags,ierr)
         integer, intent(in) :: n,ml,mu,nzmax,ldab
         real(dp), intent(in) :: ab(ldab,n)
         integer, intent(inout) :: colptr(n+1),rowind(nzmax)
         real(dp), intent(inout) :: values(nzmax)
!         real(dp), intent(in) :: ab(:,:) ! (ldab,n)
!         integer, intent(inout) :: colptr(:) ! (n+1)
!         integer, intent(inout) :: rowind(:) ! (nzmax)
!         real(dp), intent(inout) :: values(:) ! (nzmax)
         logical, intent(in) :: diags
         integer, intent(out) :: nz,ierr 
         integer :: i, j
         if (ml+mu+1 > n) then
            ierr = -1
            return
         end if
         ierr = 0
         nz = 0
         do j=1,n
            colptr(j) = nz+1 ! index in values of first entry in this column
            do i=max(1,j-mu),min(n,j+ml)
               if (ab(ldab-ml+i-j,j) == 0) then
                  if (i /= j) cycle ! not a diagonal
                  if (.not. diags) cycle
                  ! else include diagonals even if are == 0
               end if
               nz = nz+1
               if (nz > nzmax) then
                  ierr = j
                  return
               end if
               values(nz) = ab(ldab-ml+i-j,j)
               rowind(nz) = i
            end do
         end do
         colptr(n+1) = nz+1      
      end subroutine do_band_to_column_sparse
      
      
      subroutine do_column_sparse_to_band(n,ml,mu,ab,ldab,nz,colptr,rowind,values,ierr)
         integer, intent(in) :: n,ml,mu,nz,ldab
         
         real(dp), intent(inout) :: ab(ldab,n)
         integer, intent(in) :: colptr(n+1),rowind(nz)
         real(dp), intent(in) :: values(nz)
!         real(dp), intent(inout) :: ab(:,:) ! (ldab,n)
!         integer, intent(inout) :: colptr(:) ! (n+1)
!         integer, intent(inout) :: rowind(:) ! (nzmax)
!         real(dp), intent(in) :: values(:) ! (nz)
         integer, intent(out) :: ierr 
         integer :: i,j,k
         ierr = 0
         ab = 0
         do j=1,n
            do k=colptr(j),colptr(j+1)-1
               i = rowind(k) 
               if (i > j+ml .or. i < j-mu) then
                  ierr = j
                  return
               endif
               ab(ldab-ml+i-j,j) = values(k)
            end do
         end do
      end subroutine do_column_sparse_to_band

      
      subroutine do_band_to_row_sparse(n,ml,mu,ab,ldab,nzmax,nz,rowptr,colind,values,diags,ierr)
         integer, intent(in) :: n,ml,mu,nzmax,ldab
         
         real(dp), intent(in) :: ab(ldab,n)
         integer, intent(inout) :: rowptr(n+1),colind(nzmax)
         real(dp), intent(inout) :: values(nzmax)
!         real(dp), intent(in) :: ab(:,:) ! (ldab,n)
!         integer, intent(out) :: rowptr(:) ! (n+1)
!         integer, intent(out) :: colind(:) ! (nzmax)
!         real(dp), intent(inout) :: values(:) ! (nzmax)
         logical, intent(in) :: diags
         integer, intent(out) :: ierr, nz
         integer :: idiag, op_err, j1, j2, k, i, nz1
         integer, dimension(num_chunks) :: nz_per_chunk, nz_start, nz_max, i_lo, i_hi
         
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'
         
         if (dbg) write(*,*) 'enter do_band_to_row_sparse'
         
         if (ml+mu+1 > n) then
            ierr = -1
            if (dbg) then
               write(*,*) 'do_band_to_row_sparse'
               write(*,*) 'ml+mu+1', ml+mu+1
               write(*,*) 'n', n
            end if
            return
         end if
         
         ierr = 0
         nz = 0
         idiag = ldab - ml
         
         nz_start(1) = 1
         i_lo(1) = 1
         do k = 2, num_chunks
            nz_start(k) = ((k-1)*nzmax)/num_chunks
            nz_max(k-1) = nz_start(k) - 1
            i_lo(k) = ((k-1)*n)/num_chunks
            i_hi(k-1) = i_lo(k) - 1
         end do
         nz_max(num_chunks) = nzmax
         i_hi(num_chunks) = n
         
         if (dbg) write(*,*) 'do_band_to_row_sparse - do chunks'
         op_err = 0
!$OMP PARALLEL DO PRIVATE(k,op_err)
         do k = 1, num_chunks
            call do_chunk_band_to_row_sparse( &
               k, ldab, n, nzmax, nz_per_chunk, nz_start, nz_max, ml, mu, idiag, diags, &
               i_lo, i_hi, ab, rowptr, colind, values, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
         if (dbg) write(*,*) 'do_band_to_row_sparse - done chunks'
         
         if (ierr /= 0) then
         
            write(*,*) 'do_band_to_row_sparse: failed to fit in chunks'
            write(*,*) 'please increase the max fill factor for your sparse matrix'
            call mesa_error(__FILE__,__LINE__)
         
         else 
            
            if (dbg) then
               do k=1,num_chunks
                  write(*,2) 'k', k
                  write(*,2) 'nz_per_chunk(k)', nz_per_chunk(k)
                  write(*,2) 'nz_start(k)', nz_start(k)
                  write(*,2) 'nz_max(k)', nz_max(k)
                  write(*,2) 'i_lo(k)', i_lo(k)
                  write(*,2) 'i_hi(k)', i_hi(k)
                  write(*,*)
               end do
            end if
            
            ! reposition the chunk results
            if (dbg) write(*,*) 'reposition the chunk results'
            i = nz_per_chunk(1)
            do k = 2, num_chunks
               nz1 = nz_per_chunk(k)
               if (dbg) write(*,2) 'i', i
               if (dbg) write(*,2) 'nz1', nz1
               if (dbg) write(*,2) 'k', k
               if (dbg) write(*,2) 'nz_start(k)', nz_start(k)
               j2 = nz_start(k)
               do j1 = i+1, i+nz1 ! avoid ifort segfault
                  values(j1) = values(j2)
                  colind(j1) = colind(j2)
                  j2 = j2+1
               end do
               if (dbg) write(*,2) 'i_lo(k)', i_lo(k)
               if (dbg) write(*,2) 'i_hi(k)', i_hi(k)
               if (dbg) write(*,2) 'nz_start(k)-i-1', nz_start(k)-i-1
               j2 = (nz_start(k)-i-1)
               do j1 = i_lo(k), i_hi(k)
                  rowptr(j1) = rowptr(j1) - j2
               end do
               i = i+nz1
            end do
            nz = i
         end if
         

         rowptr(n+1) = nz+1
         !write(*,*) 'done do_band_to_row_sparse - fill fraction', dble(nz)/dble(nzmax)
                  
      end subroutine do_band_to_row_sparse


      subroutine do_chunk_band_to_row_sparse( &
            num, ldab, n, nzmax, nz_per_chunk, nz_start, nz_max, ml, mu, idiag, diags, &
            i_lo, i_hi, ab, rowptr, colind, values, ierr)
         integer, intent(in) :: num, ldab, n, nzmax, ml, mu, idiag
         logical, intent(in) :: diags

         real(dp), intent(in) :: ab(ldab,n)
         integer, intent(inout) :: rowptr(n+1), colind(nzmax)
         real(dp), intent(inout) :: values(nzmax)
!         real(dp), intent(in) :: ab(:,:) ! (ldab,n)
!         integer, intent(inout) :: rowptr(:) ! (n+1)
!         integer, intent(inout) :: colind(:) ! (nzmax)
!         real(dp), intent(inout) :: values(:) ! (nzmax)
         integer, dimension(num_chunks) :: nz_per_chunk, nz_start, nz_max, i_lo, i_hi
         integer, intent(out) :: ierr
         
         integer :: i, j, nz
         real(dp) :: val
         
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         
         ierr = 0
         nz = nz_start(num) - 1
         do i = i_lo(num), i_hi(num)
            rowptr(i) = nz+1 ! index in values of first entry in this row
            !write(*,*) 'set rowptr(i)', i, rowptr(i)
            do j = max(1,i-ml), min(n,i+mu)
               val = ab(idiag+i-j,j)
               if (val == 0) then
                  if (i /= j) cycle ! not a diagonal, so skip if 0
                  if (.not. diags) cycle
                  ! if (diags) then include diagonals even if are == 0
               end if
               nz = nz+1
               if (nz > nz_max(num)) then
                  ierr = i
                  if (dbg) then
                     write(*,*) 'nz > nz_max(num)', nz, nz_max(num), num
                  end if
                  return
               end if
               values(nz) = val
               colind(nz) = j
            end do
         end do
         nz_per_chunk(num) = nz - nz_start(num) + 1 ! number of non-zero values for this chunk
      end subroutine do_chunk_band_to_row_sparse


      
      subroutine do_row_sparse_to_band(n,ml,mu,ab,ldab,nz,rowptr,colind,values,ierr)
         integer, intent(in) :: n,ml,mu,nz,ldab
         real(dp), intent(inout) :: ab(ldab,n)
         integer, intent(in) :: rowptr(n+1),colind(nz)
         real(dp), intent(in) :: values(nz)
!         real(dp), intent(inout) :: ab(:,:) ! (ldab,n)
!         integer, intent(inout) :: rowptr(:) ! (n+1)
!         integer, intent(inout) :: colind(:) ! (nz)
!         real(dp), intent(in) :: values(:) ! (nz)
         integer, intent(out) :: ierr 
         integer :: i,j,k
         ierr = 0
         ab = 0
         do i=1,n
            do k=rowptr(i),rowptr(i+1)-1
               j = colind(k) 
               if (i > j+ml .or. i < j-mu) then
                  ierr = j
                  return
               endif
               ab(ldab-ml+i-j,j) = values(k)
            end do
         end do
      end subroutine do_row_sparse_to_band


      ! sparse conversion based on similar routines from sparskit_src/formats.f
      
      subroutine do_dense_to_row_sparse(n,ndim,a,nzmax,nz,rowptr,colind,values,diags,ierr)
         integer, intent(in) :: n,ndim,nzmax
         !real(dp), intent(in) :: a(ndim,n)
         !integer, intent(inout) :: rowptr(n+1),colind(nzmax)
         !real(dp), intent(inout) :: values(nzmax)
         real(dp), intent(in) :: a(:,:) ! (ndim,n)
         integer, intent(inout) :: rowptr(:) ! (n+1)
         integer, intent(inout) :: colind(:) ! (nzmax)
         real(dp), intent(inout) :: values(:) ! (nzmax)
         logical, intent(in) :: diags
         integer, intent(out) :: nz,ierr         
         integer :: i,j         
         ierr = 0
         nz = 0
         do i=1,n
            rowptr(i) = nz+1 ! index in values of first entry in this row
            do j=1,n
               if (a(i,j) == 0) then
                  if (i /= j) cycle ! not a diagonal
                  if (.not. diags) cycle
                  ! else include diagonals even if are == 0
               end if
               nz = nz+1
               if (nz > nzmax) then
                  ierr = i
                  return
               end if
               values(nz) = a(i,j)
               colind(nz) = j
            end do
         end do
         rowptr(n+1) = nz+1      
      end subroutine do_dense_to_row_sparse
      
      
      subroutine do_dense_to_row_sparse_0_based( &
            n,ndim,a,nzmax,nz,rowptr,colind,values,diags,ierr)
         integer, intent(in) :: n,ndim,nzmax
         !real(dp), intent(in) :: a(ndim,n)
         !integer, intent(inout) :: rowptr(n+1),colind(nzmax)
         !real(dp), intent(inout) :: values(nzmax)
         real(dp), intent(in) :: a(:,:) ! (ndim,n)
         integer, intent(inout) :: rowptr(:) ! (n+1)
         integer, intent(inout) :: colind(:) ! (nzmax)
         real(dp), intent(inout) :: values(:) ! (nzmax)
         logical, intent(in) :: diags
         integer, intent(out) :: nz,ierr         
         integer :: i,j         
         ierr = 0
         nz = 0
         do i=1,n
            rowptr(i) = nz ! index in values of first entry in this row
            do j=1,n
               if (a(i,j) == 0) then
                  if (i /= j) cycle ! not a diagonal
                  if (.not. diags) cycle
                  ! else include diagonals even if are == 0
               end if
               nz = nz+1
               if (nz > nzmax) then
                  ierr = i
                  return
               end if
               values(nz) = a(i,j)
               colind(nz) = j-1
            end do
         end do
         rowptr(n+1) = nz      
      end subroutine do_dense_to_row_sparse_0_based


      subroutine do_row_sparse_to_dense(n,ndim,a,nz,rowptr,colind,values,ierr) 
         integer, intent(in) :: n,ndim,nz
         real(dp), intent(inout) :: a(ndim,n)
         integer, intent(in) :: rowptr(n+1),colind(nz)
         real(dp), intent(in) :: values(nz)
!         real(dp), intent(inout) :: a(:,:) ! (ndim,n)
!         integer, intent(inout) :: rowptr(:) ! (n+1)
!         integer, intent(inout) :: colind(:) ! (nz)
!         real(dp), intent(inout) :: values(:) ! (nz)
         integer, intent(out) :: ierr     
         integer :: i,j,k
         ierr = 0
         a = 0
         do i=1,n
            do k=rowptr(i),rowptr(i+1)-1
               j = colind(k) 
               if (j > n) then
                  ierr = i
                  return
               endif
               a(i,j) = values(k)
            end do
         end do
      end subroutine do_row_sparse_to_dense


      subroutine do_row_sparse_0_based_to_dense(n,ndim,a,nz,rowptr,colind,values,ierr) 
         integer, intent(in) :: n,ndim,nz
         real(dp), intent(inout) :: a(ndim,n)
         integer, intent(in) :: rowptr(0:n),colind(0:nz-1)
         real(dp), intent(in) :: values(nz)
         integer, intent(out) :: ierr     
         integer :: i,j,k
         ierr = 0
         a = 0
         do i=1,n
            do k=rowptr(i),rowptr(i+1)-1
               j = colind(k) 
               if (j > n) then
                  ierr = i
                  return
               endif
               a(i,j) = values(k)
            end do
         end do
      end subroutine do_row_sparse_0_based_to_dense

      
      subroutine do_dense_to_column_sparse(n,ndim,a,nzmax,nz,colptr,rowind,values,diags,ierr)
         integer, intent(in) :: n,ndim,nzmax
         !real(dp), intent(in) :: a(ndim,n)
         !integer, intent(inout) :: colptr(n+1),rowind(nzmax)
         !real(dp), intent(inout) :: values(nzmax)
         real(dp), intent(in) :: a(:,:) ! (ndim,n)
         integer, intent(inout) :: colptr(:) ! (n+1)
         integer, intent(inout) :: rowind(:) ! (nzmax)
         real(dp), intent(inout) :: values(:) ! (nzmax)
         logical, intent(in) :: diags
         integer, intent(out) :: nz,ierr         
         integer :: i,j         
         ierr = 0
         nz = 0
         do j=1,n
            colptr(j) = nz+1 ! index in values of first entry in this column
            do i=1,n
               if (a(i,j) == 0) then
                  if (i /= j) cycle ! not a diagonal
                  if (.not. diags) cycle
                  ! else include diagonals even if are == 0
               end if
               nz = nz+1
               if (nz > nzmax) then
                  ierr = j
                  return
               end if
               values(nz) = a(i,j)
               rowind(nz) = i
            end do
         end do
         colptr(n+1) = nz+1      
      end subroutine do_dense_to_column_sparse

      
      subroutine do_dense_to_col_sparse_0_based( &
            n,ndim,a,nzmax,nz,colptr,rowind,values,diags,ierr)
         integer, intent(in) :: n,ndim,nzmax
         !real(dp), intent(in) :: a(ndim,n)
         !integer, intent(inout) :: colptr(n+1),rowind(nzmax)
         !real(dp), intent(inout) :: values(nzmax)
         real(dp), intent(in) :: a(:,:) ! (ndim,n)
         integer, intent(inout) :: colptr(:) ! (n+1)
         integer, intent(inout) :: rowind(:) ! (nzmax)
         real(dp), intent(inout) :: values(:) ! (nzmax)
         logical, intent(in) :: diags
         integer, intent(out) :: nz,ierr         
         integer :: i,j         
         ierr = 0
         nz = 0
         do j=1,n
            colptr(j) = nz ! index in values of first entry in this column
            do i=1,n
               if (a(i,j) == 0) then
                  if (i /= j) cycle ! not a diagonal
                  if (.not. diags) cycle
                  ! else include diagonals even if are == 0
               end if
               nz = nz+1
               if (nz > nzmax) then
                  ierr = j
                  return
               end if
               values(nz) = a(i,j)
               rowind(nz) = i-1
            end do
         end do
         colptr(n+1) = nz  
      end subroutine do_dense_to_col_sparse_0_based

      
      subroutine do_dense_to_col_sparse_0_based_qp( &
            n,ndim,a,nzmax,nz,colptr,rowind,values,diags,ierr)
         integer, intent(in) :: n,ndim,nzmax
         !real(dp), intent(in) :: a(ndim,n)
         !integer, intent(inout) :: colptr(n+1),rowind(nzmax)
         !real(dp), intent(inout) :: values(nzmax)
         real(dp), intent(in) :: a(:,:) ! (ndim,n)
         integer, intent(inout) :: colptr(:) ! (n+1)
         integer, intent(inout) :: rowind(:) ! (nzmax)
         real(qp), intent(out) :: values(:) ! (nzmax)
         logical, intent(in) :: diags
         integer, intent(out) :: nz,ierr         
         integer :: i,j         
         ierr = 0
         nz = 0
         do j=1,n
            colptr(j) = nz ! index in values of first entry in this column
            do i=1,n
               if (a(i,j) == 0) then
                  if (i /= j) cycle ! not a diagonal
                  if (.not. diags) cycle
                  ! else include diagonals even if are == 0
               end if
               nz = nz+1
               if (nz > nzmax) then
                  ierr = j
                  return
               end if
               values(nz) = a(i,j)
               rowind(nz) = i-1
            end do
         end do
         colptr(n+1) = nz  
      end subroutine do_dense_to_col_sparse_0_based_qp

      
      subroutine do_quad_dense_to_col_sparse_0_based( &
            n,ndim,a,nzmax,nz,colptr,rowind,values,diags,ierr)
         integer, intent(in) :: n,ndim,nzmax
         real(qp), intent(in) :: a(:,:) ! (ndim,n)
         integer, intent(inout) :: colptr(:) ! (n+1)
         integer, intent(inout) :: rowind(:) ! (nzmax)
         real(qp), intent(out) :: values(:) ! (nzmax)
         logical, intent(in) :: diags
         integer, intent(out) :: nz,ierr         
         integer :: i,j         
         ierr = 0
         nz = 0
         do j=1,n
            colptr(j) = nz ! index in values of first entry in this column
            do i=1,n
               if (a(i,j) == 0) then
                  if (i /= j) cycle ! not a diagonal
                  if (.not. diags) cycle
                  ! else include diagonals even if are == 0
               end if
               nz = nz+1
               if (nz > nzmax) then
                  ierr = j
                  return
               end if
               values(nz) = a(i,j)
               rowind(nz) = i-1
            end do
         end do
         colptr(n+1) = nz  
      end subroutine do_quad_dense_to_col_sparse_0_based


      subroutine do_column_sparse_to_dense(n,ndim,a,nz,colptr,rowind,values,ierr) 
         integer, intent(in) :: n,ndim,nz
         real(dp), intent(inout) :: a(ndim,n)
         integer, intent(in) :: colptr(n+1),rowind(nz)
         real(dp), intent(in) :: values(nz)
!         real(dp), intent(inout) :: a(:,:) ! (ndim,n)
!         integer, intent(in) :: colptr(:) ! (n+1)
!         integer, intent(in) :: rowind(:) ! (nz)
!         real(dp), intent(in) :: values(:) ! (nz)
         integer, intent(out) :: ierr     
         integer :: i,j,k
         ierr = 0
         a = 0
         do j=1,n
            do k=colptr(j),colptr(j+1)-1
               i = rowind(k) 
               if (i > n) then
                  ierr = j
                  return
               endif
               a(i,j) = values(k)
            end do
         end do
      end subroutine do_column_sparse_to_dense


      subroutine do_quad_column_sparse_to_dense(n,ndim,a,nz,colptr,rowind,values,ierr) 
         integer, intent(in) :: n,ndim,nz
         real(qp), intent(out) :: a(ndim,n)
         integer, intent(in) :: colptr(n+1),rowind(nz)
         real(qp), intent(in) :: values(nz)
         integer, intent(out) :: ierr     
         integer :: i,j,k
         ierr = 0
         a = 0
         do j=1,n
            do k=colptr(j),colptr(j+1)-1
               i = rowind(k) 
               if (i > n) then
                  ierr = j
                  return
               endif
               a(i,j) = values(k)
            end do
         end do
      end subroutine do_quad_column_sparse_to_dense


      subroutine do_col_sparse_0_based_to_dense(n,ndim,a,nz,colptr,rowind,values,ierr) 
         integer, intent(in) :: n,ndim,nz
         real(dp), intent(inout) :: a(ndim,n)
         integer, intent(in) :: colptr(n+1),rowind(nz)
         real(dp), intent(in) :: values(nz)
!         real(dp), intent(inout) :: a(:,:) ! (ndim,n)
!         integer, intent(in) :: colptr(:) ! (n+1)
!         integer, intent(in) :: rowind(:) ! (nz)
!         real(dp), intent(in) :: values(:) ! (nz)
         integer, intent(out) :: ierr     
         integer :: i,j,k
         ierr = 0
         a = 0
         do j=1,n
            do k=colptr(j)+1,colptr(j+1)
               i = rowind(k)+1
               if (i > n) then
                  ierr = j
                  return
               endif
               a(i,j) = values(k)
            end do
         end do
      end subroutine do_col_sparse_0_based_to_dense


      subroutine do_block_dble_mv(nvar, nz, lblk, dblk, ublk, b, prod)
         use my_lapack95_dble, only: my_gemv_p1
         integer, intent(in) :: nvar, nz    
         real(dp), pointer, dimension(:,:,:), intent(in) :: lblk, dblk, ublk ! (nvar,nvar,nz)
         real(dp), pointer, dimension(:,:), intent(in) :: b ! (nvar,nz)
         real(dp), pointer, dimension(:,:), intent(inout) :: prod ! (nvar,nz)         
         integer :: k        
         do k = 1, nz
            prod(1:nvar,k) = 0
            call my_gemv_p1(nvar,nvar,dblk(1:nvar,1:nvar,k),nvar,b(1:nvar,k),prod(1:nvar,k))
            if (k > 1) then
               call my_gemv_p1(nvar,nvar,lblk(1:nvar,1:nvar,k),nvar,b(1:nvar,k-1),prod(1:nvar,k))
            end if
            if (k < nz) then
               call my_gemv_p1(nvar,nvar,ublk(1:nvar,1:nvar,k),nvar,b(1:nvar,k+1),prod(1:nvar,k))
            end if
         end do      
      end subroutine do_block_dble_mv                  
      
      
      subroutine do_block_mv_quad(lblk, dblk, ublk, b, prod)
         real(qp), pointer, dimension(:,:,:), intent(in) :: lblk, dblk, ublk ! (nvar,nvar,nz)
         real(qp), pointer, dimension(:,:), intent(in) :: b ! (nvar,nz)
         real(qp), pointer, dimension(:,:), intent(inout) :: prod ! (nvar,nz)         
         integer :: nvar, nz, k        
         include 'formats.dek'
         nvar = size(b,dim=1)
         nz = size(b,dim=2)         
!$OMP PARALLEL DO PRIVATE(k)
         do k = 1, nz
            prod(:,k) = 0
            call qdgemv(nvar,nvar,dblk(:,:,k),nvar,b(:,k),prod(:,k))
            if (k > 1) then
               call qdgemv(nvar,nvar,lblk(:,:,k),nvar,b(:,k-1),prod(:,k))
            end if
            if (k < nz) then
               call qdgemv(nvar,nvar,ublk(:,:,k),nvar,b(:,k+1),prod(:,k))
            end if
         end do      
!$OMP END PARALLEL DO

         contains
         
         subroutine qdgemv(m,n,a,lda,x,y)
            ! y := alpha*a*x + beta*y
         use const_def, only: dp
      
         integer lda,m,n
         real(qp) a(lda,*),x(*),y(*)
         real(qp) :: tmp
         ! trans = 'n'
         ! alpha = 1
         ! beta = 1
         ! incx = 1
         ! incy = 1
         integer :: j, i
           do j = 1,n
               tmp = x(j)
               if (tmp.ne.0d0) then
                   do i = 1,m
                       y(i) = y(i) + tmp*a(i,j)
                   end do
               end if
           end do
         end subroutine qdgemv
         
         
      end subroutine do_block_mv_quad
      
      
      subroutine do_LU_factored_block_dble_mv(lblk, dblk, ublk, b, ipiv, prod)
         real(dp), pointer, dimension(:,:,:), intent(in) :: lblk, dblk, ublk ! (nvar,nvar,nz)
         real(dp), pointer, dimension(:,:), intent(in) :: b ! (nvar,nz)
         integer, intent(in) :: ipiv(:,:) ! (nvar,nz)
         real(dp), pointer, dimension(:,:), intent(inout) :: prod ! (nvar,nz)                  
         integer :: nvar, nz, k, incx, incy         
         nvar = size(b,dim=1)
         nz = size(b,dim=2)         
         incx = 1
         incy = 1         
!$OMP PARALLEL DO PRIVATE(k)
         do k = 1, nz
            call do_LU_factored_square_mv(nvar,dblk(:,:,k),b(:,k),ipiv(:,k),prod(:,k))
            if (k > 1) then
               call dgemv('N',nvar,nvar,1d0,lblk(:,:,k),nvar,b(:,k-1),incx,1d0,prod(:,k),incy)
            end if
            if (k < nz) then
               call dgemv('N',nvar,nvar,1d0,ublk(:,:,k),nvar,b(:,k+1),incx,1d0,prod(:,k),incy)
            end if
         end do      
!$OMP END PARALLEL DO
      end subroutine do_LU_factored_block_dble_mv
      
      
      subroutine do_LU_factored_square_mv(m,a,b,ipiv,x) ! set x = A*b
         ! A factored in LU manner = P*L*U.
         integer, intent(in) :: m
         real(dp), intent(in) :: a(:,:) ! (lda,m), lda >= m
         real(dp), intent(in) :: b(:) ! (m)
         integer, intent(in) :: ipiv(:) ! (m)
         real(dp), intent(inout) :: x(:) ! (m)
         integer :: i, j
         real(dp), dimension(m) :: y
         include 'formats.dek'
         ! y = U*b
         do i=1,m
            y(i) = 0
            do j=i,m
               y(i) = y(i) + a(i,j)*b(j)
            end do
         end do
         ! x = L*y
         do i=1,m
            x(i) = y(i)
            do j=1,i-1
               x(i) = x(i) + a(i,j)*y(j)
            end do
         end do
         ! x = P*x
         call dlaswp(1, x, m, 1, m, ipiv, -1)
      end subroutine do_LU_factored_square_mv
      
      
      subroutine do_LU_factored_square_mm(m,A,B,ipiv,C) ! set C = A*B
         ! A factored in LU manner = P*L*U.
         integer, intent(in) :: m
         real(dp), intent(in) :: A(:,:) ! (m,m)
         real(dp), intent(in) :: B(:,:) ! (m,m)
         integer, intent(in) :: ipiv(:) ! (m)
         real(dp), intent(inout) :: C(:,:) ! (m,m)
         integer :: i, j, k
         real(dp), dimension(m,m) :: Y
         include 'formats.dek'
         ! Y = U*B
         do i=1,m
            do j=1,m
               Y(i,j) = 0
               do k=i,m
                  Y(i,j) = Y(i,j) + A(i,k)*B(k,j)
               end do
            end do
         end do
         ! C = L*Y
         do i=1,m
            do j=1,m
               C(i,j) = Y(i,j)
               do k=1,i-1
                  C(i,j) = C(i,j) + A(i,k)*Y(k,j)
               end do
            end do
         end do
         ! C = P*C
         call dlaswp(m, C, m, 1, m, ipiv, -1)
      end subroutine do_LU_factored_square_mm


      subroutine do_multiply_xa(n, A1, x, b)
         !  calculates b = x*A
         integer, intent(in) :: n
         real(dp), pointer, intent(in) :: A1(:) ! =(n, n)
         real(dp), pointer, intent(in) :: x(:) ! (n)
         real(dp), pointer, intent(inout) :: b(:) ! (n)
         integer :: i, j
         real(dp), pointer :: A(:,:) ! (n, n)
         A(1:n,1:n) => A1(1:n*n)
         do j = 1, n
            b(j) = dot_product(x(1:n),A(1:n,j))
         end do
      end subroutine do_multiply_xa


      subroutine do_quad_multiply_xa(n, A1, x, b)
         !  calculates b = x*A
         integer, intent(in) :: n
         real(qp), pointer, intent(in) :: A1(:) ! =(n, n)
         real(qp), pointer, intent(in) :: x(:) ! (n)
         real(qp), pointer, intent(inout) :: b(:) ! (n)
         integer :: i, j
         real(qp), pointer :: A(:,:) ! (n, n)
         A(1:n,1:n) => A1(1:n*n)
         do j = 1, n
            b(j) = dot_product(x(1:n),A(1:n,j))
         end do
      end subroutine do_quad_multiply_xa


      subroutine do_multiply_xa_plus_c(n, A1, x, c, b)
         !  calculates b = x*A + c
         integer, intent(in) :: n
         real(dp), pointer, intent(in) :: A1(:) ! =(n,n)
         real(dp), pointer, intent(in) :: x(:) ! (n)
         real(dp), pointer, intent(in) :: c(:) ! (n)
         real(dp), pointer, intent(inout) :: b(:) ! (n)
         integer :: i, j
         real(dp), pointer :: A(:,:) ! (n,n)
         A(1:n,1:n) => A1(1:n*n)
         do j = 1, n
            b(j) = dot_product(x(1:n),A(1:n,j)) + c(j)
         end do
      end subroutine do_multiply_xa_plus_c


      subroutine do_quad_multiply_xa_plus_c(n, A1, x, c, b)
         !  calculates b = x*A + c
         integer, intent(in) :: n
         real(qp), pointer, intent(in) :: A1(:) ! =(n, n)
         real(qp), pointer, intent(in) :: x(:) ! (n)
         real(qp), pointer, intent(in) :: c(:) ! (n)
         real(qp), pointer, intent(inout) :: b(:) ! (n)
         integer :: i, j
         real(qp), pointer :: A(:,:) ! (n, n)
         A(1:n,1:n) => A1(1:n*n)
         do j = 1, n
            b(j) = dot_product(x(1:n),A(1:n,j)) + c(j)
         end do
      end subroutine do_quad_multiply_xa_plus_c


      subroutine do_block_multiply_xa(nvar, nz, lblk1, dblk1, ublk1, x1, b1)
         !  calculates b = x*A
         integer, intent(in) :: nvar, nz
         real(dp), dimension(:), intent(in), pointer :: lblk1, dblk1, ublk1 ! =(nvar,nvar,nz)
         real(dp), intent(in), pointer :: x1(:) ! =(nvar,nz)
         real(dp), intent(inout), pointer :: b1(:) ! =(nvar,nz)
         integer :: k, shift, shift2, nvar2
         real(dp), pointer, dimension(:) :: p1, p2, p3, p4
         nvar2 = nvar*nvar
!$OMP PARALLEL DO PRIVATE(k,shift,shift2,p1,p2,p3,p4)
         do k = 1, nz
            shift = nvar*(k-1)
            shift2 = nvar2*(k-1)
            p1(1:nvar2) => dblk1(shift2+1:shift2+nvar2)
            p2(1:nvar) => x1(shift+1:shift+nvar)
            p3(1:nvar) => b1(shift+1:shift+nvar)
            call do_multiply_xa(nvar,p1,p2,p3)
            if (k > 1) then
               p1(1:nvar2) => ublk1(shift2+1:shift2+nvar2)
               p2(1:nvar) => x1(shift+1:shift+nvar)
               p3(1:nvar) => b1(shift+1:shift+nvar)
               p4(1:nvar) => b1(shift+1:shift+nvar)
               call do_multiply_xa_plus_c(nvar,p1,p2,p3,p4)
            end if
            if (k < nz) then
               p1(1:nvar2) => lblk1(shift2+1:shift2+nvar2)
               p2(1:nvar) => x1(shift+1+nvar:shift+2*nvar)
               p3(1:nvar) => b1(shift+1:shift+nvar)
               p4(1:nvar) => b1(shift+1:shift+nvar)
               call do_multiply_xa_plus_c(nvar,p1,p2,p3,p4)
            end if
         end do
!$OMP END PARALLEL DO
      end subroutine do_block_multiply_xa


      subroutine do_quad_block_multiply_xa(nvar, nz, lblk1, dblk1, ublk1, x1, b1)
         !  calculates b = x*A
         integer, intent(in) :: nvar, nz
         real(qp), dimension(:), intent(in), pointer :: lblk1, dblk1, ublk1 ! =(nvar,nvar,nz)
         real(qp), intent(in), pointer :: x1(:) ! =(nvar,nz)
         real(qp), intent(inout), pointer :: b1(:) ! =(nvar,nz)
         integer :: k, shift, shift2, nvar2
         real(qp), pointer, dimension(:) :: p1, p2, p3, p4
         nvar2 = nvar*nvar
!$OMP PARALLEL DO PRIVATE(k,shift,shift2,p1,p2,p3,p4)
         do k = 1, nz
            shift = nvar*(k-1)
            shift2 = nvar2*(k-1)
            p1(1:nvar2) => dblk1(shift2+1:shift2+nvar2)
            p2(1:nvar) => x1(shift+1:shift+nvar)
            p3(1:nvar) => b1(shift+1:shift+nvar)
            call do_quad_multiply_xa(nvar,p1,p2,p3)
            if (k > 1) then
               p1(1:nvar2) => ublk1(shift2+1:shift2+nvar2)
               p2(1:nvar) => x1(shift+1:shift+nvar)
               p3(1:nvar) => b1(shift+1:shift+nvar)
               p4(1:nvar) => b1(shift+1:shift+nvar)
               call do_quad_multiply_xa_plus_c(nvar,p1,p2,p3,p4)
            end if
            if (k < nz) then
               p1(1:nvar2) => lblk1(shift2+1:shift2+nvar2)
               p2(1:nvar) => x1(shift+1+nvar:shift+2*nvar)
               p3(1:nvar) => b1(shift+1:shift+nvar)
               p4(1:nvar) => b1(shift+1:shift+nvar)
               call do_quad_multiply_xa_plus_c(nvar,p1,p2,p3,p4)
            end if
         end do
!$OMP END PARALLEL DO
      end subroutine do_quad_block_multiply_xa


      subroutine do_band_multiply_xa(n, kl, ku, ab1, ldab, x, b)
         !  calculates b = x*a = transpose(a)*x
            integer, intent(in) :: n
         !          the number of linear equations, i.e., the order of the
         !          matrix a.  n >= 0.
            integer, intent(in) :: kl
         !          the number of subdiagonals within the band of a.  kl >= 0.
            integer, intent(in) :: ku
         !          the number of superdiagonals within the band of a.  ku >= 0.
            integer, intent(in) :: ldab
         !          the leading dimension of the array ab.  ldab >= kl+ku+1.
            real(dp), intent(in), pointer :: ab1(:) ! =(ldab, n)
         !          the matrix a in band storage, in rows 1 to kl+ku+1;
         !          the j-th column of a is stored in the j-th column of the
         !          array ab as follows:
         !          ab(ku+1+i-j, j) = a(i, j) for max(1, j-ku)<=i<=min(n, j+kl)
            real(dp), intent(in), pointer :: x(:) ! (n)
         !          the input vector to be multiplied by the matrix.
            real(dp), intent(inout), pointer :: b(:) ! (n)
         !          on exit, set to matrix product of x*a = b
         integer ::  i, j, k
         real(dp), pointer :: ab(:,:)
         ab(1:ldab,1:n) => ab1(1:ldab*n)
         do j = 1, n
            k = ku+1-j
            b(j) = 0
            do i = max(1,j-ku), min(n,j+kl)
               b(j) = b(j) + x(i)*ab(k+i,j)
            end do
         end do
      end subroutine do_band_multiply_xa      


      subroutine do_quad_band_multiply_xa(n, kl, ku, ab, ldab, x, b)
         !  calculates b = x*a = transpose(a)*x
            integer, intent(in) :: n
         !          the number of linear equations, i.e., the order of the
         !          matrix a.  n >= 0.
            integer, intent(in) :: kl
         !          the number of subdiagonals within the band of a.  kl >= 0.
            integer, intent(in) :: ku
         !          the number of superdiagonals within the band of a.  ku >= 0.
            integer, intent(in) :: ldab
         !          the leading dimension of the array ab.  ldab >= kl+ku+1.
            real(qp), intent(in) :: ab(:,:) ! (ldab, n)
         !          the matrix a in band storage, in rows 1 to kl+ku+1;
         !          the j-th column of a is stored in the j-th column of the
         !          array ab as follows:
         !          ab(ku+1+i-j, j) = a(i, j) for max(1, j-ku)<=i<=min(n, j+kl)
            real(qp), intent(in) :: x(:) ! (n)
         !          the input vector to be multiplied by the matrix.
            real(qp), intent(inout) :: b(:) ! (n)
         !          on exit, set to matrix product of x*a = b
         integer ::  i, j, k
         do j = 1, n
            k = ku+1-j
            b(j) = 0
            do i = max(1, j-ku), min(n, j+kl)
               b(j) = b(j) + x(i)*ab(k+i, j)
            end do
         end do
      end subroutine do_quad_band_multiply_xa      
      
      
      
      subroutine do_clip_blocks( &
            mblk, clip_limit, lmat, dmat, umat, dmat_nnz, total_nnz)
         integer, intent(in) :: mblk
         real(dp), intent(in) :: clip_limit
         real(dp), intent(inout) :: lmat(:,:), dmat(:,:), umat(:,:)
         integer, intent(inout) :: dmat_nnz, total_nnz
         integer :: i, j
         dmat_nnz = 0; total_nnz = 0
         do j=1,mblk
            do i=1,mblk
               if (i /= j .and. abs(lmat(i,j)) < clip_limit) lmat(i,j) = 0d0
               if (lmat(i,j) /= 0) total_nnz = total_nnz + 1
               if (i /= j .and. abs(dmat(i,j)) < clip_limit) dmat(i,j) = 0d0
               if (dmat(i,j) /= 0) then
                  total_nnz = total_nnz + 1
                  dmat_nnz = dmat_nnz + 1
               end if
               if (i /= j .and. abs(umat(i,j)) < clip_limit) umat(i,j) = 0d0
               if (umat(i,j) /= 0) total_nnz = total_nnz + 1
            end do
         end do
      end subroutine do_clip_blocks
      
      
      subroutine do_clip_block(mblk, clip_limit, dmat, dmat_nnz)
         integer, intent(in) :: mblk
         real(dp), intent(in) :: clip_limit
         real(dp), intent(inout) :: dmat(:,:)
         integer, intent(inout) :: dmat_nnz
         integer :: i, j
         dmat_nnz = 0
         do j=1,mblk
            do i=1,mblk
               if (i /= j .and. abs(dmat(i,j)) < clip_limit) dmat(i,j) = 0d0
               if (dmat(i,j) /= 0) dmat_nnz = dmat_nnz + 1
            end do
         end do
      end subroutine do_clip_block


      subroutine read_hbcode1(iounit, nrow, ncol, nnzero, values, rowind, colptr, ierr)

      CHARACTER      TITLE*72 , KEY*8    , MXTYPE*3 , &
                     PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20

      INTEGER        TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, &
                     iounit, NROW  , NCOL  , NNZERO, NELTVL

      INTEGER, pointer ::        COLPTR (:), ROWIND (:)

      REAL(dp), pointer ::        VALUES (:)
      integer, intent(out) :: ierr
      
      integer i
      ierr = 0
      READ (iounit, 1000, iostat=ierr ) TITLE , KEY   , &
                           TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, &
                           MXTYPE, NROW  , NCOL  , NNZERO, NELTVL, &
                           PTRFMT, INDFMT, VALFMT, RHSFMT
      if (ierr /= 0) return
 1000 FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )
 
      allocate(VALUES(NNZERO), ROWIND(NNZERO), COLPTR(NCOL+1), stat=ierr)
      if (ierr /= 0) return

      READ (iounit, PTRFMT, iostat=ierr ) ( COLPTR (I), I = 1, NCOL+1 )
      if (ierr /= 0) return

      READ (iounit, INDFMT, iostat=ierr ) ( ROWIND (I), I = 1, NNZERO )
      if (ierr /= 0) return

      IF  ( VALCRD .GT. 0 )  THEN

!         ----------------------
!         ... READ MATRIX VALUES
!         ----------------------

          READ (iounit, VALFMT, iostat=ierr ) ( VALUES (I), I = 1, NNZERO )
          if (ierr /= 0) return

      ENDIF


      end subroutine read_hbcode1


      subroutine read_hbcode1_quad(iounit, nrow, ncol, nnzero, values, rowind, colptr, ierr)

      CHARACTER      TITLE*72 , KEY*8    , MXTYPE*3 , &
                     PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20

      INTEGER        TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, &
                     iounit, NROW  , NCOL  , NNZERO, NELTVL

      INTEGER, pointer ::        COLPTR (:), ROWIND (:)

      REAL(qp), pointer ::        VALUES (:)
      integer, intent(out) :: ierr
      
      integer i
      ierr = 0
      READ (iounit, 1000, iostat=ierr ) TITLE , KEY   , &
                           TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, &
                           MXTYPE, NROW  , NCOL  , NNZERO, NELTVL, &
                           PTRFMT, INDFMT, VALFMT, RHSFMT
      if (ierr /= 0) return
 1000 FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )
 
      allocate(VALUES(NNZERO), ROWIND(NNZERO), COLPTR(NCOL+1), stat=ierr)
      if (ierr /= 0) return

      READ (iounit, PTRFMT, iostat=ierr ) ( COLPTR (I), I = 1, NCOL+1 )
      if (ierr /= 0) return

      READ (iounit, INDFMT, iostat=ierr ) ( ROWIND (I), I = 1, NNZERO )
      if (ierr /= 0) return

      IF  ( VALCRD .GT. 0 )  THEN

!         ----------------------
!         ... READ MATRIX VALUES
!         ----------------------

          READ (iounit, VALFMT, iostat=ierr ) ( VALUES (I), I = 1, NNZERO )
          if (ierr /= 0) return

      ENDIF


      end subroutine read_hbcode1_quad


      subroutine write_hbcode1(iounit, nrow, ncol, nnzero, values, rowind, colptr, ierr)

      CHARACTER      TITLE*72 , KEY*8    , MXTYPE*3 , &
                     PTRFMT*16, INDFMT*16, use_VALFMT*20, VALFMT*20, RHSFMT*20

      INTEGER        TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, &
                     iounit, NROW  , NCOL  , NNZERO, NELTVL

      INTEGER        COLPTR (*), ROWIND (*), ierr

      REAL(dp)         VALUES (*)
      
      integer i
      
      ierr = 0

!    ------------------------
!     ... WRITE HEADER BLOCK
!     ------------------------

!  Line 1 (A72, A8)
!       Col. 1 - 72   Title (TITLE)
!       Col. 73 - 80  Matrix name / identifier (MTRXID)
! 
!  Line 2 (I14, 3(1X, I13))
!       Col. 1 - 14   Total number of lines excluding header (TOTCRD)
!       Col. 16 - 28  Number of lines for pointers (PTRCRD)
!       Col. 30 - 42  Number of lines for row (or variable) indices (INDCRD)
!       Col. 44 - 56  Number of lines for numerical values (VALCRD)
! 
!  Line 3 (A3, 11X, 4(1X, I13))
!       Col. 1 - 3    Matrix type (see below) (MXTYPE)
!       Col. 15 - 28  Compressed Column: Number of rows (NROW)
!                     Elemental: Largest integer used to index variable (MVAR)
!       Col. 30 - 42  Compressed Column: Number of columns (NCOL)
!                     Elemental: Number of element matrices (NELT)
!       Col. 44 - 56  Compressed Column: Number of entries (NNZERO)
!                     Elemental: Number of variable indeces (NVARIX)
!       Col. 58 - 70  Compressed Column: Unused, explicitly zero
!                     Elemental: Number of elemental matrix entries (NELTVL)
! 
!  Line 4 (2A16, A20)
!       Col. 1 - 16   Fortran format for pointers (PTRFMT)
!       Col. 17 - 32  Fortran format for row (or variable) indices (INDFMT)
!       Col. 33 - 52  Fortran format for numerical values of coefficient matrix
!                     (VALFMT)
!                     (blank in the case of matrix patterns)
! 
!  The three character type field on line 3 describes the matrix type.
!  The following table lists the permitted values for each of the three
!  characters. As an example of the type field, RSA denotes that the matrix
!  is real, symmetric, and assembled.
! 
!  First Character:
!       R Real matrix
!       C Complex matrix
!       I integer matrix
!       P Pattern only (no numerical values supplied)
!       Q Pattern only (numerical values supplied in associated auxiliary value
!         file)
! 
!  Second Character:
!       S Symmetric
!       U Unsymmetric
!       H Hermitian
!       Z Skew symmetric
!       R Rectangular
! 
!  Third Character:
!       A Compressed column form
!       E Elemental form
! 



      TITLE = ''
      KEY = ''
      
      PTRFMT = '(10I8)'
      INDFMT = '(12I6)'
      use_VALFMT = '(5(1pE27.16))'
      VALFMT = '(5E27.16)'
      RHSFMT = ''

      PTRCRD = (NCOL+1)/10 + 1 ! number of lines for COLPTR
      INDCRD = NNZERO/12 + 1 ! number of lines for ROWIND
      VALCRD = NNZERO/5 + 1 ! number of lines for VALUES
      RHSCRD = 0
      TOTCRD = 3 + PTRCRD + INDCRD + VALCRD + RHSCRD
      
      MXTYPE = 'RUA'
      NELTVL = 0
      
      WRITE (iounit, 1000 ) TITLE , KEY   , &
                           TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, &
                           MXTYPE, NROW  , NCOL  , NNZERO, NELTVL, &
                           PTRFMT, INDFMT, VALFMT, RHSFMT
 1000 FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )

!     -------------------------
!     ... WRITE MATRIX STRUCTURE
!     -------------------------

      WRITE (iounit, PTRFMT ) ( COLPTR (I), I = 1, NCOL+1 )

      WRITE (iounit, INDFMT ) ( ROWIND (I), I = 1, NNZERO )

      IF  ( VALCRD .GT. 0 )  THEN

!         ----------------------
!         ... WRITE MATRIX VALUES
!         ----------------------

          WRITE (iounit, use_VALFMT ) ( VALUES (I), I = 1, NNZERO )

      ENDIF

      return
      end subroutine write_hbcode1
         
         
      subroutine read_block_tridiagonal(iounit,nvar,nblk,lblk1,dblk1,ublk1,ierr)
         integer, intent(in) :: iounit
         integer, intent(out) :: nvar, nblk
         real(dp), pointer, dimension(:) :: lblk1,dblk1,ublk1 ! will be allocated
         integer, intent(out) :: ierr
         integer :: k
         real(dp), pointer, dimension(:,:,:) :: lblk,dblk,ublk
         ierr = 0
         read(iounit,*,iostat=ierr) nvar, nblk
         if (ierr /= 0) return
         allocate(lblk1(nvar*nvar*nblk), dblk1(nvar*nvar*nblk), ublk1(nvar*nvar*nblk), stat=ierr)
         if (ierr /= 0) return
         lblk(1:nvar,1:nvar,1:nblk) => lblk1(1:nvar*nvar*nblk)
         dblk(1:nvar,1:nvar,1:nblk) => dblk1(1:nvar*nvar*nblk)
         ublk(1:nvar,1:nvar,1:nblk) => ublk1(1:nvar*nvar*nblk)
         do k=1,nblk
            if (k > 1) then
               call read1_sparse_block(iounit, nvar, lblk(:,:,k), ierr)
               if (ierr /= 0) return
            end if
            call read1_sparse_block(iounit, nvar, dblk(:,:,k), ierr)
            if (ierr /= 0) return
            if (k < nblk) then
               call read1_sparse_block(iounit, nvar, ublk(:,:,k), ierr)
               if (ierr /= 0) return
            end if
         end do
         
      end subroutine read_block_tridiagonal
         
         
      subroutine read_quad_block_tridiagonal(iounit,nvar,nblk,lblk1,dblk1,ublk1,ierr)
         integer, intent(in) :: iounit
         integer, intent(out) :: nvar, nblk
         real(qp), pointer, dimension(:) :: lblk1,dblk1,ublk1 ! will be allocated
         integer, intent(out) :: ierr
         integer :: k
         real(qp), pointer, dimension(:,:,:) :: lblk,dblk,ublk
         ierr = 0
         read(iounit,*,iostat=ierr) nvar, nblk
         if (ierr /= 0) return
         allocate(lblk1(nvar*nvar*nblk), dblk1(nvar*nvar*nblk), ublk1(nvar*nvar*nblk), stat=ierr)
         if (ierr /= 0) return
         lblk(1:nvar,1:nvar,1:nblk) => lblk1(1:nvar*nvar*nblk)
         dblk(1:nvar,1:nvar,1:nblk) => dblk1(1:nvar*nvar*nblk)
         ublk(1:nvar,1:nvar,1:nblk) => ublk1(1:nvar*nvar*nblk)
         do k=1,nblk
            if (k > 1) then
               call read1_sparse_block_quad(iounit, nvar, lblk(:,:,k), ierr)
               if (ierr /= 0) return
            end if
            call read1_sparse_block_quad(iounit, nvar, dblk(:,:,k), ierr)
            if (ierr /= 0) return
            if (k < nblk) then
               call read1_sparse_block_quad(iounit, nvar, ublk(:,:,k), ierr)
               if (ierr /= 0) return
            end if
         end do
         
      end subroutine read_quad_block_tridiagonal
      
      
      subroutine read1_sparse_block(iounit, nvar, blk, ierr)
         integer, intent(in) :: iounit, nvar
         real(dp) :: blk(:,:) ! (nvar,nvar)
         integer, intent(out) :: ierr
         integer :: nnz, nrow, ncol
         integer, pointer :: rowind(:), colptr(:)
         real(dp), pointer :: values(:)
         ierr = 0
         call read_hbcode1(iounit, nrow, ncol, nnz, values, rowind, colptr,ierr)
         if (ierr /= 0 .or. nrow /= nvar .or. nrow /= ncol) return
         call do_column_sparse_to_dense(nrow,ncol,blk,nnz,colptr,rowind,values,ierr) 
         deallocate(colptr,rowind,values)
      end subroutine read1_sparse_block
      
      
      subroutine read1_sparse_block_quad(iounit, nvar, blk, ierr)
         integer, intent(in) :: iounit, nvar
         real(qp) :: blk(:,:) ! (nvar,nvar)
         integer, intent(out) :: ierr
         integer :: nnz, nrow, ncol
         integer, pointer :: rowind(:), colptr(:)
         real(qp), pointer :: values(:)
         ierr = 0
         call read_hbcode1_quad(iounit, nrow, ncol, nnz, values, rowind, colptr,ierr)
         if (ierr /= 0 .or. nrow /= nvar .or. nrow /= ncol) return
         call do_quad_column_sparse_to_dense(nrow,ncol,blk,nnz,colptr,rowind,values,ierr) 
         deallocate(colptr,rowind,values)
      end subroutine read1_sparse_block_quad
         
         
      subroutine write_block_tridiagonal(iounit,nvar,nblk,lblk,dblk,ublk,ierr)
         integer, intent(in) :: iounit, nvar, nblk
         real(dp), intent(in), dimension(:,:,:) :: lblk,dblk,ublk
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
         write(iounit,*) nvar, nblk
         do k=1,nblk
            if (k > 1) then
               call write1_sparse_block(iounit, nvar, lblk(:,:,k), ierr)
               if (ierr /= 0) return
            end if
            call write1_sparse_block(iounit, nvar, dblk(:,:,k), ierr)
            if (ierr /= 0) return
            if (k < nblk) then
               call write1_sparse_block(iounit, nvar, ublk(:,:,k), ierr)
               if (ierr /= 0) return
            end if
         end do
      end subroutine write_block_tridiagonal
      
      
      subroutine write1_sparse_block(iounit, nvar, blk, ierr)
         integer, intent(in) :: iounit, nvar
         real(dp), intent(in) :: blk(:,:) ! (nvar,nvar)
         integer, intent(out) :: ierr
         integer :: nnz, rowind(nvar*nvar), colptr(nvar+1)
         real(dp) :: values(nvar*nvar)
         call do_dense_to_column_sparse( &
            nvar, nvar, blk, nvar*nvar, nnz, colptr, rowind, values, .true., ierr)
         if (ierr /= 0) return
         call write_hbcode1(iounit, nvar, nvar, nnz, values, rowind, colptr, ierr)
      end subroutine write1_sparse_block


      end module mtx_support
