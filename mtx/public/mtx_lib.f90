! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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


module mtx_lib
   
   use const_def, only : dp, qp
   
   implicit none


contains
   
   ! mesa includes sources for a subset of BLAS and dble.
   ! you can use those, or, better yet, you can use a package optimized
   ! for your machine such as GotoBLAS or Intel's MKL.
   ! see utils/makefile_header for details.
   
   ! see mtx/blas_src for the subset of BLAS routines included in mtx_lib
   ! see mtx/dble_src for the subset of dble routines included in mtx_lib
   
   ! subroutines for dense and banded matrix decompositions and solves

include "mtx_dble_decsol.dek" ! dble versions
   
   !> Wraps the lapack DGBSVX routine for banded matrices to be used with block tridiagonal matrices.
   !!
   !! @params ublk The upper blocks. Shape is (nvars, nvars, nblocks-1).
   !! @params lblk The lower blocks. Shape is (nvars, nvars, nblocks-1).
   !! @params dblk The diagonal blocks. Shape is (nvars, nvars, nblocks).
   !! @params nblocks The number of blocks on the diagonal.
   !! @params nvars The width and height of each block (e.g. each block is nvars by nvars).
   !! @params x Array of shape (nblocks,nvar). Stores the result of solving Matrix.x = b.
   !! @params b Array of shape (nblocks,nvar). The right-hand side from the definition of x.
   !! @params ierr Integer-valued error code.
   subroutine DGBSVX_block_tridiagonal_padded(ublk, lblk, dblk, x, b, nblocks, nvar, rcond, ierr)
      use DGBSVX_wrapper, only:DGBSVX_tridiagonal_wrapper
      
      ! Inputs
      real(dp), dimension(:, :, :), intent(in) :: ublk, lblk, dblk
   real(dp), dimension(:, :), intent(in) :: b
   integer, intent(in) :: nblocks, nvar
   
   ! Intermediates
   real(dp) :: pre_conditioner(nvar, nblocks)
   
   ! Outputs
   real(dp), dimension(:, :), intent(out) :: x
   real(dp), intent(out) :: rcond
   integer, intent(out) :: ierr
   
   call DGBSVX_tridiagonal_wrapper(ublk(:, :, 1:nblocks-1), lblk(:, :, 2:nblocks), dblk, pre_conditioner, x, b, nblocks, nvar, rcond, ierr)
      
      end subroutine DGBSVX_block_tridiagonal_padded
      
      
      !> Wraps the lapack DGBSVX routine for banded matrices to be used with block tridiagonal matrices.
      !! This version returns the pre-conditioning array which was used to rescale each row.
      !!
      !! @params ublk The upper blocks. Shape is (nvars, nvars, nblocks-1).
      !! @params lblk The lower blocks. Shape is (nvars, nvars, nblocks-1).
      !! @params dblk The diagonal blocks. Shape is (nvars, nvars, nblocks).
      !! @params nblocks The number of blocks on the diagonal.
      !! @params nvars The width and height of each block (e.g. each block is nvars by nvars).
      !! @params x Array of shape (nblocks,nvar). Stores the result of solving Matrix.x = b.
      !! @params b Array of shape (nblocks,nvar). The right-hand side from the definition of x.
      !! @params pre_conditioner Pre-conditioning vector. Computed in this routine.
      !! @params ierr Integer-valued error code.
      subroutine DGBSVX_block_tridiagonal_pc(ublk, lblk, dblk, pre_conditioner, x, b, nblocks, nvar, ierr)
      use DGBSVX_wrapper, only:DGBSVX_tridiagonal_wrapper
      
      ! Inputs
      real(dp), dimension(:, :, :), intent(in) :: ublk, lblk, dblk
   real(dp), dimension(:, :), intent(in) :: b
   integer, intent(in) :: nblocks, nvar
   
   ! Intermediates
   real(dp) :: rcond
   
   ! Outputs
   real(dp), dimension(:, :), intent(out) :: pre_conditioner
   real(dp), dimension(:, :), intent(out) :: x
   integer, intent(out) :: ierr
   
   call DGBSVX_tridiagonal_wrapper(ublk, lblk, dblk, pre_conditioner, x, b, nblocks, nvar, rcond, ierr)
      
      end subroutine DGBSVX_block_tridiagonal_pc
      
      !> Wraps the lapack DGBSVX routine for banded matrices to be used with block tridiagonal matrices.
      !!
      !! @params ublk The upper blocks. Shape is (nvars, nvars, nblocks-1).
      !! @params lblk The lower blocks. Shape is (nvars, nvars, nblocks-1).
      !! @params dblk The diagonal blocks. Shape is (nvars, nvars, nblocks).
      !! @params nblocks The number of blocks on the diagonal.
      !! @params nvars The width and height of each block (e.g. each block is nvars by nvars).
      !! @params x Array of shape (nblocks,nvar). Stores the result of solving Matrix.x = b.
      !! @params b Array of shape (nblocks,nvar). The right-hand side from the definition of x.
      !! @params ierr Integer-valued error code.
      subroutine DGBSVX_block_tridiagonal(ublk, lblk, dblk, x, b, nblocks, nvar, ierr)
   use DGBSVX_wrapper, only:DGBSVX_tridiagonal_wrapper
   
   ! Inputs
   real(dp), dimension(:, :, :), intent(in) :: ublk, lblk, dblk
   real(dp), dimension(:, :), intent(in) :: b
   integer, intent(in) :: nblocks, nvar
   
   ! Intermediates
   real(dp) :: pre_conditioner(nvar, nblocks)
      real(dp) :: rcond
   
   
   ! Outputs
   real(dp), dimension(:, :), intent(out) :: x
   integer, intent(out) :: ierr
   
   call DGBSVX_tridiagonal_wrapper(ublk, lblk, dblk, pre_conditioner, x, b, nblocks, nvar, rcond, ierr)
   
   end subroutine DGBSVX_block_tridiagonal
   
   !> Wraps the lapack DGBSVX routine for banded matrices with a preconditioner
   !! This version returns the pre-conditioning array which was used to rescale each row.
   !!
   !! @param matrix_size The number of rows or columns in the square matrix.
   !! @params n_upper_bands The number of bands above the diagonal.
   !! @params n_lower_bands The number of bands below the diagonal.
   !! @params bands Array of shape (n_upper_bands + n_lower_bands + 1, matrix_size). Stores the bands going from upper-most to lower-most, aligned so the first entry in each band has index 1.
   !! @params x Array of shape (matrix_size). Stores the result of solving Matrix.x = b.
   !! @params b Array of shape (matrix_size). The right-hand side from the definition of x.
   !! @params pre_conditioner Pre-conditioning vector. Computed in this routine.
   !! @params ierr Integer-valued error code.
   subroutine DGBSVX_banded_pc(matrix_size, n_upper_bands, n_lower_bands, bands, x, b, pre_conditioner, ierr)
      use DGBSVX_wrapper, only:DGBSVX_banded_wrapper
      
      ! Inputs
      real(dp), dimension(:, :), intent(in) :: bands
   real(dp), dimension(:), intent(in) :: b
   integer, intent(in) :: matrix_size, n_upper_bands, n_lower_bands
   
   ! Outputs
   real(dp), dimension(:), intent(out) :: x
   real(dp), dimension(:), intent(out) :: pre_conditioner
   integer, intent(out) :: ierr
   
   call DGBSVX_banded_wrapper(matrix_size, n_upper_bands, n_lower_bands, bands, x, b, pre_conditioner, ierr)
      
      end subroutine DGBSVX_banded_pc
      
      
      !> Wraps the lapack DGBSVX routine for banded matrices with a preconditioner
      !!
      !! @param matrix_size The number of rows or columns in the square matrix.
      !! @params n_upper_bands The number of bands above the diagonal.
      !! @params n_lower_bands The number of bands below the diagonal.
      !! @params bands Array of shape (n_upper_bands + n_lower_bands + 1, matrix_size). Stores the bands going from upper-most to lower-most, aligned so the first entry in each band has index 1.
      !! @params x Array of shape (matrix_size). Stores the result of solving Matrix.x = b.
      !! @params b Array of shape (matrix_size). The right-hand side from the definition of x.
      !! @params ierr Integer-valued error code.
      subroutine DGBSVX_banded(matrix_size, n_upper_bands, n_lower_bands, bands, x, b, ierr)
      use DGBSVX_wrapper, only:DGBSVX_banded_wrapper
      
      ! Inputs
      real(dp), dimension(:, :), intent(in) :: bands
   real(dp), dimension(:), intent(in) :: b
   integer, intent(in) :: matrix_size, n_upper_bands, n_lower_bands
   
   ! Intermediates
   real(dp) :: pre_conditioner(matrix_size)
      
      ! Outputs
      real(dp), dimension(:), intent(out) :: x
      integer, intent(out) :: ierr
   
   call DGBSVX_banded_wrapper(matrix_size, n_upper_bands, n_lower_bands, bands, x, b, pre_conditioner, ierr)
      
      end subroutine DGBSVX_banded
      
      
      ! sometimes you just need a null version of a routine
      include "mtx_null_decsol.dek"
      
      ! sometimes you need to debug a jacobian by saving it to plotting data files
      include "mtx_debug_decsol.dek"
         
         ! sparse matrices come in many formats.
         ! for example, compressed row sparse format is used by SPARSKIT,
         ! while compressed column sparse format is used by Super_LU.
         ! here are conversion routines for these two options.
         include "mtx_formats.dek"
      
      subroutine mtx_write_hbcode1(iounit, n, nnzero, values, rowind, colptr, ierr)
      use mtx_support, only:write_hbcode1
      integer, intent(in) :: iounit, n, nnzero
   integer :: rowind(:) ! (nnzero)
   integer :: colptr(:) ! (n+1)
   real(dp) :: values(:) ! (nnzero)
   integer, intent(out) :: ierr
   call write_hbcode1(iounit, n, n, nnzero, values, rowind, colptr, ierr)
      end subroutine mtx_write_hbcode1
      
      
      subroutine mtx_write_block_tridiagonal(iounit, nvar, nblk, lblk, dblk,ublk, ierr)
   use mtx_support, only:write_block_tridiagonal
   integer, intent(in) :: iounit, nvar, nblk
   real(dp), intent(in), dimension(:, :, :) :: lblk, dblk, ublk
   integer, intent(out) :: ierr
   call write_block_tridiagonal(iounit, nvar, nblk, lblk, dblk,ublk, ierr)
      end subroutine mtx_write_block_tridiagonal
      
      subroutine mtx_read_block_tridiagonal(iounit, nvar, nblk, lblk1, dblk1,ublk1, ierr)
      use mtx_support, only:read_block_tridiagonal
      integer, intent(in) :: iounit
      integer, intent(out) :: nvar, nblk
   real(dp), pointer, dimension(:) :: lblk1, dblk1, ublk1 ! =(nvar,nvar,nblk) will be allocated
   integer, intent(out) :: ierr
   call read_block_tridiagonal(iounit, nvar, nblk, lblk1, dblk1,ublk1, ierr)
      end subroutine mtx_read_block_tridiagonal
      
      subroutine mtx_read_quad_block_tridiagonal(iounit, nvar, nblk, lblk1, dblk1,ublk1, ierr)
   use mtx_support, only:read_quad_block_tridiagonal
   integer, intent(in) :: iounit
   integer, intent(out) :: nvar, nblk
   real(qp), pointer, dimension(:) :: lblk1, dblk1, ublk1 ! =(nvar,nvar,nblk) will be allocated
   integer, intent(out) :: ierr
   call read_quad_block_tridiagonal(iounit, nvar, nblk, lblk1, dblk1,ublk1, ierr)
      end subroutine mtx_read_quad_block_tridiagonal
      
      ! BCYCLIC multi-thread block tridiagonal
      include "mtx_bcyclic_dble_decsol.dek"
      ! S.P.Hirshman, K.S.Perumalla, V.E.Lynch, & R.Sanchez,
      ! BCYCLIC: A parallel block tridiagonal matrix cyclic solver,
      ! J. Computational Physics, 229 (2010) 6392-6404.
      
      
      subroutine block_dble_mv(nvar, nz, lblk, dblk, ublk, b, prod)
   ! set prod = A*b with A = block tridiagonal given by lblk, dblk, ublk
   use mtx_support, only:do_block_dble_mv
   integer, intent(in) :: nvar, nz
   real(dp), pointer, dimension(:, :, :), intent(in) :: lblk, dblk, ublk ! (nvar,nvar,nz)
   real(dp), pointer, dimension(:, :), intent(in) :: b ! (nvar,nz)
   real(dp), pointer, dimension(:, :), intent(inout) :: prod ! (nvar,nz)
   call do_block_dble_mv(nvar, nz, lblk, dblk, ublk, b, prod)
      end subroutine block_dble_mv
      
      
      subroutine block_quad_mv(lblk, dblk, ublk, b, prod)
      ! set prod = A*b with A = block tridiagonal given by lblk, dblk, ublk
      use mtx_support, only:do_block_mv_quad
      real(qp), pointer, dimension(:, :, :), intent(in) :: lblk, dblk, ublk ! (nvar,nvar,nz)
   real(qp), pointer, dimension(:, :), intent(in) :: b ! (nvar,nz)
   real(qp), pointer, dimension(:, :), intent(inout) :: prod ! (nvar,nz)
   call do_block_mv_quad(lblk, dblk, ublk, b, prod)
      end subroutine block_quad_mv
      
      
      
      subroutine multiply_xa(n, A1, x, b)
      !  calculates b = x*A
      use mtx_support, only:do_multiply_xa
      integer, intent(in) :: n
      real(dp), pointer, intent(in) :: A1(:) ! =(n, n)
   real(dp), pointer, intent(in) :: x(:) ! (n)
   real(dp), pointer, intent(inout) :: b(:) ! (n)
   call do_multiply_xa(n, A1, x, b)
      end subroutine multiply_xa
      
      
      subroutine block_multiply_xa(nvar, nz, lblk1, dblk1, ublk1, x1, b1)
   !  calculates b = x*A
   use mtx_support, only:do_block_multiply_xa
   integer, intent(in) :: nvar, nz
   real(dp), dimension(:), intent(in), pointer :: lblk1, dblk1, ublk1 ! =(nvar,nvar,nz)
   real(dp), intent(in), pointer :: x1(:) ! =(nvar,nz)
   real(dp), intent(inout), pointer :: b1(:) ! =(nvar,nz)
   call do_block_multiply_xa(nvar, nz, lblk1, dblk1, ublk1, x1, b1)
      end subroutine block_multiply_xa
      
      
      subroutine quad_block_multiply_xa(nvar, nz, lblk1, dblk1, ublk1, x1, b1)
      !  calculates b = x*A
      use mtx_support, only:do_quad_block_multiply_xa
      integer, intent(in) :: nvar, nz
      real(qp), dimension(:), intent(in), pointer :: lblk1, dblk1, ublk1 ! =(nvar,nvar,nz)
   real(qp), intent(in), pointer :: x1(:) ! =(nvar,nz)
   real(qp), intent(inout), pointer :: b1(:) ! =(nvar,nz)
   call do_quad_block_multiply_xa(nvar, nz, lblk1, dblk1, ublk1, x1, b1)
      end subroutine quad_block_multiply_xa
      
      
      subroutine band_multiply_xa(n, kl, ku, ab1, ldab, x, b)
      !  calculates b = x*a = transpose(a)*x
      use mtx_support, only:do_band_multiply_xa
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
   call do_band_multiply_xa(n, kl, ku, ab1, ldab, x, b)
      end subroutine band_multiply_xa
      
      
      include "mtx_lapack95.dek"
      
      
      ! utilities for working with jacobians
      include "mtx_jac.dek"
      
      ! the following call dble routines to estimate matrix condition numbers.
      include "mtx_rcond.dek"
      
      integer function decsol_option(which_decsol_option, ierr)
   use mtx_def
   character (len = *), intent(in) :: which_decsol_option
   integer, intent(out) :: ierr
   character (len = 64) :: option
   ierr = 0
   option = which_decsol_option
      
      if (option == 'lapack') then
   decsol_option = lapack
   
   else if (option == 'bcyclic_dble') then
   decsol_option = bcyclic_dble
   
   else
   ierr = -1
   decsol_option = -1
   end if
      end function decsol_option
      
      
      subroutine decsol_option_str(which_decsol_option, decsol_option, ierr)
      use mtx_def
      integer, intent(in) :: which_decsol_option
      character (len = *), intent(out) :: decsol_option
   integer, intent(out) :: ierr
   ierr = 0
   
   if (which_decsol_option == lapack) then
   decsol_option = 'lapack'
   else if (which_decsol_option == bcyclic_dble) then
   decsol_option = 'bcyclic_dble'
   
   else
   ierr = -1
   decsol_option = ''
   end if
   
   end subroutine decsol_option_str
   
   
   logical function is_block_tridiagonal_decsol(which_decsol_option)
   use mtx_def
   integer, intent(in) :: which_decsol_option
   is_block_tridiagonal_decsol = (which_decsol_option == bcyclic_dble)
end function is_block_tridiagonal_decsol


end module mtx_lib
