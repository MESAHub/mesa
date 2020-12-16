module pre_conditioners

      use const_def

      implicit none

      private
      public :: compute_block_preconditioner, compute_band_preconditioner, compute_block_left_preconditioner

      contains

      !> Computes a pre-conditioning array for rows of a block triadiagonal matrix.
      !! Each row in the matrix can then be divided by the corresponding element in this array
      !! to form a usually better-conditioned matrix.
      !! This element is chosen to be the absolute value of the element in that row with the greatest magnitude.
      !!
      !! @params ublk The upper blocks. Shape is (nvars, nvars, nblocks-1).
      !! @params lblk The lower blocks. Shape is (nvars, nvars, nblocks-1).
      !! @params dblk The diagonal blocks. Shape is (nvars, nvars, nblocks).
      !! @params nblocks The number of blocks on the diagonal.
      !! @params nvars The width and height of each block (e.g. each block is nvars by nvars).
      !! @params pre_conditioner Pre-conditioning vector. Computed in this routine.
      subroutine compute_block_preconditioner(ublk, lblk, dblk, nblocks, nvar, pre_conditioner)
         real(dp), dimension(:,:,:), intent(in) :: ublk, lblk, dblk
         integer, intent(in) :: nblocks, nvar

         integer :: j, k

         real(dp), dimension(:,:), intent(out) :: pre_conditioner

         pre_conditioner = 0d0
         do j=1,nblocks
            do k=1,nvar
               pre_conditioner(k,j) = pre_conditioner(k,j) + sum(abs(dblk(k,:,j)))

               if (j < nblocks) then
                  pre_conditioner(k,j) = pre_conditioner(k,j) + sum(abs(ublk(k,:,j))) 
               end if

               if (j > 1) then
                  pre_conditioner(k,j) = pre_conditioner(k,j) + sum(abs(lblk(k,:,j-1))) 
               end if
            end do
         end do

      end subroutine compute_block_preconditioner

      !> Computes a pre-conditioning array for columns of a block triadiagonal matrix.
      !! Each column in the matrix can then be divided by the corresponding element in this array
      !! to form a usually better-conditioned matrix.
      !! This element is chosen to be the absolute value of the element in that column with the greatest magnitude.
      !! Note that the final answer to the linear solve must be divided element-wise by this pre-conditioner.
      !!
      !! @params ublk The upper blocks. Shape is (nvars, nvars, nblocks-1).
      !! @params lblk The lower blocks. Shape is (nvars, nvars, nblocks-1).
      !! @params dblk The diagonal blocks. Shape is (nvars, nvars, nblocks).
      !! @params nblocks The number of blocks on the diagonal.
      !! @params nvars The width and height of each block (e.g. each block is nvars by nvars).
      !! @params pre_conditioner Pre-conditioning vector. Computed in this routine.
      subroutine compute_block_left_preconditioner(ublk, lblk, dblk, nblocks, nvar, &
                                                   right_pre_conditioner, pre_conditioner)
         real(dp), dimension(:,:,:), intent(in) :: ublk, lblk, dblk
         real(dp), dimension(:,:), intent(in) :: right_pre_conditioner
         integer, intent(in) :: nblocks, nvar

         integer :: j, k, k2

         real(dp), dimension(:,:), intent(out) :: pre_conditioner

         pre_conditioner = 0d0
         do j=1,nblocks
            do k=1,nvar
               do k2=1,nvar
                  pre_conditioner(k,j) = pre_conditioner(k,j) + abs(dblk(k2,k,j) / right_pre_conditioner(k2,j))

                  if (j < nblocks) then
                     pre_conditioner(k,j) = pre_conditioner(k,j) + abs(lblk(k2,k,j) / right_pre_conditioner(k2,j+1))
                  end if

                  if (j > 1) then
                     pre_conditioner(k,j) = pre_conditioner(k,j) + abs(ublk(k2,k,j-1) / right_pre_conditioner(k2,j-1))
                  end if
               end do
            end do
         end do

      end subroutine compute_block_left_preconditioner

      !> Computes a pre-conditioning array for rows of a banded matrix.
      !! Each row in the matrix can then be divided by the corresponding element in this array
      !! to form a usually better-conditioned matrix.
      !! This element is chosen to be the absolute value of the element in that row with the greatest magnitude.
      !!
      !! @params n_upper_bands The number of bands above the diagonal.
      !! @params n_lower_bands The number of bands below the diagonal.
      !! @params bands Array of shape (n_upper_bands + n_lower_bands + 1, matrix_size). Stores the bands going from upper-most to lower-most, aligned so the first entry in each band has index 1.
      !! @params pre_conditioner Pre-conditioning vector. Computed in this routine.
      subroutine compute_band_preconditioner(matrix_size, n_upper_bands, n_lower_bands, bands, pre_conditioner)
         ! Inputs
         real(dp), dimension(:,:), intent(in) :: bands
         integer, intent(in) :: matrix_size, n_upper_bands, n_lower_bands

         ! Intermediates
         integer :: i, j, n_bands

         ! Outputs
         real(dp), dimension(:), intent(out) :: pre_conditioner

         n_bands = n_upper_bands + n_lower_bands + 1

         pre_conditioner = 0d0

         ! Upper bands
         do j=1,n_upper_bands
             ! The first upper band has index 1 in bands and runs 1...matrix_size-n_upper_bands.
            do i=1,matrix_size + j - n_upper_bands - 1
               pre_conditioner(i) = max(pre_conditioner(i), abs(bands(j,i)))
            end do
         end do

         ! Diagonal band
         do i=1,matrix_size
            pre_conditioner(i) = max(pre_conditioner(i), abs(bands(n_upper_bands + 1, i)))
         end do

         ! Lower bands
         do j=1,n_lower_bands
            ! The first lower band has index n_upper_bands+2 in bands and runs 1...matrix_size-1
            do i=j+1,matrix_size
               pre_conditioner(i) = max(pre_conditioner(i), abs(bands(j + n_upper_bands + 1, i-j)))
            end do
         end do

      end subroutine compute_band_preconditioner      

end module pre_conditioners