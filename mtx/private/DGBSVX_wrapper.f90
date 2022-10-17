module DGBSVX_wrapper
   
   use const_def
   use pre_conditioners
   
   implicit none
   
   private
   public :: DGBSVX_tridiagonal_wrapper, DGBSVX_banded_wrapper

contains
   
   !> Wraps the lapack DGBSVX routine for banded matrices with a preconditioner
   !! and a simpler input format.
   !!
   !! @param matrix_size The number of rows or columns in the square matrix.
   !! @params n_upper_bands The number of bands above the diagonal.
   !! @params n_lower_bands The number of bands below the diagonal.
   !! @params bands Array of shape (n_upper_bands + n_lower_bands + 1, matrix_size). Stores the bands going from upper-most to lower-most, aligned so the first entry in each band has index 1.
   !! @params x Array of shape (matrix_size). Stores the result of solving Matrix.x = b.
   !! @params b Array of shape (matrix_size). The right-hand side from the definition of x.
   !! @params pre_conditioner Pre-conditioning vector. Computed in this routine.
   !! @params ierr Integer-valued error code.
   subroutine DGBSVX_banded_wrapper(matrix_size, n_upper_bands, n_lower_bands, bands, x, b, pre_conditioner, ierr)
      ! Inputs
      real(dp), dimension(:, :), intent(in) :: bands
      real(dp), dimension(:), intent(in) :: b
      integer, intent(in) :: matrix_size, n_upper_bands, n_lower_bands
      
      ! Intermediates
      integer counter, i, j, k, k1, k2
      character(len = 1) :: FACT, TRANS, EQUED
      integer :: NRHS, LDAB, LDAFB, LDB, LDX
      integer :: IWORK(matrix_size), IPIV(matrix_size)
      real(dp) :: WORK(3 * matrix_size), BERR(1), FERR(1)
      real(dp) :: C(matrix_size), R(matrix_size)
      real(dp) :: b_conditioned(matrix_size, 1), x_tmp(matrix_size, 1)
      real(dp) :: AB(n_lower_bands + n_upper_bands + 1, matrix_size)
      real(dp) :: AFB(2 * n_lower_bands + n_upper_bands + 1, matrix_size)
      real(dp) :: RCOND
      
      ! Outputs
      real(dp), dimension(:), intent(out) :: x
      real(dp), dimension(:), intent(out) :: pre_conditioner
      integer, intent(out) :: ierr
      
      ! Fixed arguments
      
      FACT = 'N' ! We've already equilibrated the matrix
      TRANS = 'N' ! This argument is irrelevant if FACT == 'N'
      EQUED = 'N' ! No equilibration
      NRHS = 1                ! Number of columns in equ for A.x=b.
      LDX = matrix_size       ! Number of rows in flattened x
      LDB = matrix_size       ! Number of rows in flattened b
      LDAB = n_lower_bands + n_upper_bands + 1 ! Number of rows in AB format of banded matrix
      LDAFB = 2 * n_lower_bands + n_upper_bands + 1 ! Number of rows in AB format of LU-factored matrix
      
      ! Compute preconditioner
      call compute_band_preconditioner(matrix_size, n_upper_bands, n_lower_bands, bands, pre_conditioner)
      
      ! Pre-condition b
      do i = 1, matrix_size
         b_conditioned(i, 1) = b(i) / pre_conditioner(i)
      end do
      
      ! Put bands into AB form for LAPACK
      ! We store column j in the original matrix in column j of AB.
      ! In that column, row k in the original matrix is stored in row
      ! n_upper_bands+1+k-j of AB. So the mapping from the matrix to AB is
      ! (k,j) -> (n_upper_bands+1+k-j,j)
      AB = 0d0
      
      ! In the upper bands, upper band i at index j (bands(i,j)) corresponds to
      ! position (j, n_upper_bands - i + j + 1) in the matrix. This is then
      ! position (i, n_upper_bands - i + j + 1) in AB.
      do i = 1, n_upper_bands
         do j = 1, matrix_size + i - n_upper_bands - 1
            AB(i, n_upper_bands - i + j + 1) = bands(i, j) / pre_conditioner(j)
         end do
      end do
      
      ! In the lower bands, lower band i+1+n_upper_bands at index j corresponds to
      ! position (i+j, j) in the matrix. This is then
      ! position (n_upper_bands+1+i, j) in AB.
      do i = 1, n_lower_bands
         do j = 1, matrix_size + n_lower_bands - i - 1
            AB(n_upper_bands + 1 + i, j) = bands(i + n_upper_bands + 1, j) / pre_conditioner(i + j)
         end do
      end do
      
      ! In the diagonal band, position k corresponds to index (k,k), which is position
      ! (n_upper_bands+1,k) in AB.
      do i = 1, matrix_size
         AB(n_upper_bands + 1, i) = bands(n_upper_bands + 1, i) / pre_conditioner(i)
      end do
      
      ! Solve
      call DGBSVX(FACT, TRANS, matrix_size, n_lower_bands, n_upper_bands, NRHS, AB, LDAB, AFB, LDAFB, IPIV, &
         EQUED, R, C, b_conditioned, LDB, x_tmp, LDX, RCOND, FERR, BERR, WORK, IWORK, &
         ierr)
      
      ! Check for errors
      if (RCOND < 1d-12 .or. ierr == matrix_size + 1) then
         write(*, *) 'Matrix is singular to machine precision.'
         write(*, *) 'RCOND = ', RCOND
         write(*, *) 'ierr = ', ierr
         write(*, *) 'N = ', matrix_size
         
         open(unit = 10, file = "bands.data")
         do j = 1, LDAB
            do i = 1, matrix_size
               write(10, *) bands(j, i)
            end do
         end do
      else if (ierr /= 0) then
         write(*, *) 'Matrix is exactly singular.'
         write(*, *) 'ierr = ', ierr
         write(*, *) 'N = ', matrix_size
      end if
      
      ! Store output
      do i = 1, matrix_size
         x(i) = x_tmp(i, 1)
      end do
   
   end subroutine DGBSVX_banded_wrapper
   
   !> Wraps the lapack DGBSVX routine for banded matrices to be used with block tridiagonal matrices.
   !! Also includes a pre-conditioning step and a simpler input format.
   !!
   !! @params ublk The upper blocks. Shape is (nvars, nvars, nblocks-1).
   !! @params lblk The lower blocks. Shape is (nvars, nvars, nblocks-1).
   !! @params dblk The diagonal blocks. Shape is (nvars, nvars, nblocks).
   !! @params nblocks The number of blocks on the diagonal.
   !! @params nvars The width and height of each block (e.g. each block is nvars by nvars).
   !! @params x Array of shape (nvar, nblocks). Stores the result of solving Matrix.x = b.
   !! @params b Array of shape (nvar, nblocks). The right-hand side from the definition of x.
   !! @params pre_conditioner Pre-conditioning vector. Computed in this routine.
   !! @params ierr Integer-valued error code.
   subroutine DGBSVX_tridiagonal_wrapper(ublk, lblk, dblk, pre_conditioner, x, b, nblocks, nvar, rcond, ierr)
      ! Inputs
      real(dp), dimension(:, :, :), intent(in) :: ublk, lblk, dblk
      real(dp), dimension(:, :), intent(in) :: b
      integer, intent(in) :: nblocks, nvar
      
      ! Intermediates
      integer counter, i, j, k, k1, k2
      character(len = 1) :: FACT, TRANS, EQUED
      integer :: N, KL, KU, NRHS, LDAB, LDAFB, LDB, LDX
      integer :: IWORK(nvar * nblocks), IPIV(nvar * nblocks)
      real(dp) :: WORK(3 * nvar * nblocks), BERR(1), FERR(1)
      real(dp) :: C(nvar * nblocks), R(nvar * nblocks)
      real(dp) :: b_flat(nvar * nblocks, 1), x_flat(nvar * nblocks, 1)
      real(dp) :: AB(4 * nvar + 1, nvar * nblocks)
      real(dp) :: AFB(6 * nvar + 1, nblocks * nvar)
      real(dp) :: left_pre_conditioner(nvar, nblocks)
      
      ! Outputs
      real(dp), dimension(:, :), intent(out) :: pre_conditioner
      real(dp), dimension(:, :), intent(out) :: x
      real(dp), intent(out) :: rcond
      integer, intent(out) :: ierr
      
      ! Fixed arguments
      
      FACT = 'N' ! We've already equilibrated the matrix
      TRANS = 'N' ! This argument is irrelevant if FACT == 'N'
      EQUED = 'N' ! No equilibration
      N = nblocks * nvar ! Number of equations
      KL = 2 * nvar ! Number of lower bands. We round up to the max number
      ! which could be held by the block triadiagonal matrix.
      ! If performance is ever an issue this solve can be sped
      ! up by tightening KL.
      KU = 2 * nvar ! Number of upper bands. We round up to the max number
      ! which could be held by the block triadiagonal matrix.
      ! If performance is ever an issue this solve can be sped
      ! up by tightening KU.
      NRHS = 1      ! Number of columns in equ for A.x=b.
      LDX = N       ! Number of rows in flattened x
      LDB = N       ! Number of rows in flattened b
      LDAB = KL + KU + 1 ! Number of rows in AB format of banded matrix
      LDAFB = 2 * KL + KU + 1 ! Number of rows in AB format of LU-factored matrix
      
      ! Compute preconditioners
      call compute_block_preconditioner(ublk, lblk, dblk, nblocks, nvar, pre_conditioner)
      call compute_block_left_preconditioner(ublk, lblk, dblk, nblocks, nvar, &
         pre_conditioner, left_pre_conditioner)
      
      ! We solve M_{ij} x_j = b_i
      ! We pre-condition by dividing b_i by p_i, so
      ! M_{ij} x_j / p_i = b_i / p_i
      ! We also pre-condition by dividing x_j by L_j, so
      ! (M_{ij} / (p_i L_j) ) (x_j L_j) = b_i / p_i
      ! We solve for the composite x_j' = x_j L_j, so in the end x_j = x_j' / L_j.
      
      ! Flatten b and precondition
      counter = 1
      do j = 1, nblocks
         do k = 1, nvar
            b_flat(counter, 1) = b(k, j) / pre_conditioner(k, j)
            counter = counter + 1
         end do
      end do
      
      ! Set up banded matrix
      AB = 0d0
      do i = 1, nblocks
         do j = 1, nvar
            do k = 1, nvar
               ! (j,k,i) in dblk corresponds to coordinate ((i-1)*nvar+j,(i-1)*nvar+k)
               ! in the original matrix.
               ! We store column m in the original matrix in column m of AB.
               ! In that column, row l in the original matrix is stored in row
               ! KU+1+l-m of AB. Hence
               ! (j,k,i) goes in spot (KU+1+(i-1)*nvar+j-(i-1)*nvar-k,(i-1)*nvar+k)
               ! = (KU+1+j-k,(i-1)*nvar+k) of AB.
               AB(KU + 1 + j - k, (i - 1) * nvar + k) = dblk(j, k, i) / pre_conditioner(j, i) / left_pre_conditioner(k, i)
            end do
         end do
      end do
      
      do i = 1, nblocks - 1
         do j = 1, nvar
            do k = 1, nvar
               ! (j,k,i) in lblk corresponds to coordinate (i*nvar+j,(i-1)*nvar+k)
               ! in the original matrix.
               ! We store column m in the original matrix in column m of AB.
               ! In that column, row l in the original matrix is stored in row
               ! KU+1+l-m of AB. Hence
               ! (j,k,i) goes in spot (KU+1+i*nvar+j-(i-1)*nvar-k,(i-1)*nvar+k)
               ! = (KU+1+nvar+j-k,(i-1)*nvar+k) of AB.
               AB(KU + 1 + nvar + j - k, (i - 1) * nvar + k) = lblk(j, k, i) / pre_conditioner(j, i + 1) / left_pre_conditioner(k, i)
               
               ! (j,k,i) in ublk corresponds to coordinate ((i-1)*nvar+j,i*nvar+k)
               ! in the original matrix.
               ! We store column m in the original matrix in column m of AB.
               ! In that column, row l in the original matrix is stored in row
               ! KU+1+l-m of AB. Hence
               ! (j,k,i) goes in spot (KU+1+(i-1)*nvar+j-i*nvar-k,i*nvar+k)
               ! = (KU+1-nvar+j-k,(i-1)*nvar+k) of AB.
               AB(KU + 1 - nvar + j - k, i * nvar + k) = ublk(j, k, i) / pre_conditioner(j, i) / left_pre_conditioner(k, i + 1)
            end do
         end do
      end do
      
      ! Solve
      call DGBSVX(FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, &
         EQUED, R, C, b_flat, LDB, x_flat, LDX, RCOND, FERR, BERR, WORK, IWORK, &
         ierr)
      
      ! Check for errors
      if (ierr /= 0 .or. rcond < 1d-10) then
         write(*, *) 'WARNING: Matrix is singular to machine precision.'
         write(*, *) 'RCOND = ', RCOND
         write(*, *) 'ierr = ', ierr
         write(*, *) 'N = ', N
         write(*, *) 'nvars, nblocks', nvar, nblocks
         
         open(unit = 10, file = "lower.data")
         open(unit = 11, file = "upper.data")
         open(unit = 12, file = "diagonal.data")
         do k = 1, nblocks
            do k1 = 1, nvar
               do k2 = 1, nvar
                  if (k < nblocks) then
                     write(10, *) lblk(k1, k2, k)
                     write(11, *) ublk(k1, k2, k)
                  end if
                  write(12, *) dblk(k1, k2, k)
               end do
            end do
         end do
         write(10, *) '---'
         write(11, *) '---'
         write(12, *) '---'
         ierr = N + 1
      else if (ierr /= 0) then
         write(*, *) 'WARNING: Matrix is exactly singular.'
         write(*, *) 'ierr = ', ierr
         write(*, *) 'N = ', N
      end if
      
      ! Unflatten output
      counter = 1
      do j = 1, nblocks
         do k = 1, nvar
            x(k, j) = x_flat(counter, 1) / left_pre_conditioner(k, j)
            counter = counter + 1
         end do
      end do
   
   end subroutine DGBSVX_tridiagonal_wrapper

end module DGBSVX_wrapper
