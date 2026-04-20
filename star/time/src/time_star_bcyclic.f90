! ***********************************************************************
!
!   Copyright (C) 2024  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

! This module times the bcyclic algorithm using bcyclic_factor and bcyclic_solve
! from star/private/star_bcyclic.f90

module time_star_bcyclic

   use star_private_def
   use star_bcyclic, only: bcyclic_factor, bcyclic_solve, clear_storage
   use const_def, only: dp
   use omp_lib, only: omp_get_num_threads

   implicit none

contains

   subroutine time_solvers()
      integer, parameter :: n_nz_tests = 4
      integer, parameter :: n_nvar_tests = 4
      integer :: nz_tests(n_nz_tests)
      integer :: nvar_tests(n_nvar_tests)
      integer :: i, j
      integer :: n_threads

      nz_tests = [256, 512, 1024, 2048]
      nvar_tests = [64, 128, 256, 512]

      !$omp parallel
      !$omp master
      n_threads = omp_get_num_threads()
      !$omp end master
      !$omp end parallel
      write(*, '(A, I0)') 'Number of OMP threads: ', n_threads

      write(*, '(A)') '=============================================='
      write(*, '(A)') 'Star BCYCLIC Solver Timing Tests'
      write(*, '(A)') '=============================================='
      write(*, *)

      write(*, '(A)') 'Testing bcyclic_factor and bcyclic_solve from star_bcyclic'
      write(*, '(A)') 'for various block sizes (nvar) and number of blocks (nz)'
      write(*, *)

      write(*, '(A8, A8, A18, A18)') 'nvar', 'nz', 'factor (s)', 'solve (s)'
      write(*, '(A)') '----------------------------------------------'

      do i = 1, n_nvar_tests
         do j = 1, n_nz_tests
            call run_timing_test(nvar_tests(i), nz_tests(j))
         end do
      end do

      write(*, *)
      write(*, '(A)') 'Timing tests complete.'

   end subroutine time_solvers

   subroutine run_timing_test(nvar, nz)
      integer, intent(in) :: nvar, nz

      real(dp) :: time_factor, time_solve

      call time_bcyclic(nvar, nz, time_factor, time_solve)

      write(*, '(I8, I8, ES18.6, ES18.6)') nvar, nz, time_factor, time_solve

   end subroutine run_timing_test

   subroutine time_bcyclic(nvar, nz, time_factor, time_solve)
      integer, intent(in) :: nvar, nz
      real(dp), intent(out) :: time_factor, time_solve

      type(star_info), pointer :: s
      type(star_info), target :: s_data

      real(dp), pointer :: lblk1(:), dblk1(:), ublk1(:)
      real(dp), pointer :: lblkF1(:), dblkF1(:), ublkF1(:)
      real(dp), pointer :: B1(:), soln1(:)
      real(dp), pointer :: row_scale_factors1(:), col_scale_factors1(:)
      integer, pointer :: ipivot1(:)
      character(len=:), allocatable :: equed1

      real(dp), pointer :: lblk(:,:,:), dblk(:,:,:), ublk(:,:,:)
      real(dp), pointer :: brhs(:,:)

      integer :: ierr, nvar2
      integer :: i, j, k
      integer :: count_start, count_end, count_rate

      ierr = 0
      nvar2 = nvar * nvar

      ! Set up minimal star_info structure
      s => s_data
      nullify(s% bcyclic_odd_storage)
      s% use_DGESVX_in_bcyclic = .false.
      s% report_min_rcond_from_DGESXV = .false.
      s% use_equilibration_in_DGESVX = .false.

      ! Allocate arrays
      allocate( &
         lblk1(nvar2*nz), dblk1(nvar2*nz), ublk1(nvar2*nz), &
         lblkF1(nvar2*nz), dblkF1(nvar2*nz), ublkF1(nvar2*nz), &
         B1(nvar*nz), soln1(nvar*nz), &
         ipivot1(nvar*nz), &
         row_scale_factors1(nvar*nz), col_scale_factors1(nvar*nz), &
         stat=ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in alloc for bcyclic'
         time_factor = -1.0_dp
         time_solve = -1.0_dp
         return
      end if

      allocate(character(len=nz) :: equed1)

      ! Set up pointers for convenience
      lblk(1:nvar, 1:nvar, 1:nz) => lblk1(1:nvar2*nz)
      dblk(1:nvar, 1:nvar, 1:nz) => dblk1(1:nvar2*nz)
      ublk(1:nvar, 1:nvar, 1:nz) => ublk1(1:nvar2*nz)
      brhs(1:nvar, 1:nz) => B1(1:nvar*nz)

      call setup_test_matrix(nvar, nz, lblk, dblk, ublk, brhs)

      ! Initialize scale factors
      row_scale_factors1 = 1.0_dp
      col_scale_factors1 = 1.0_dp
      equed1 = repeat('N', nz)

      ! Time factor operation
      call system_clock(count_start)

      call bcyclic_factor( &
         s, nvar, nz, lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
         B1, row_scale_factors1, col_scale_factors1, &
         equed1, 1, ierr)

      call system_clock(count_end, count_rate)
      time_factor = real(count_end - count_start, dp) / real(count_rate, dp)

      if (ierr /= 0) then
         write(*, *) 'bcyclic_factor failed, ierr=', ierr
         time_factor = -1.0_dp
         time_solve = -1.0_dp
         call cleanup()
         return
      end if

      ! Time solve operation
      call system_clock(count_start)

      call bcyclic_solve( &
         s, nvar, nz, lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
         B1, soln1, row_scale_factors1, col_scale_factors1, equed1, &
         1, ierr)

      call system_clock(count_end, count_rate)
      time_solve = real(count_end - count_start, dp) / real(count_rate, dp)

      if (ierr /= 0) then
         write(*, *) 'bcyclic_solve failed, ierr=', ierr
      end if

      call cleanup()

   contains

      subroutine cleanup()
         if (associated(s% bcyclic_odd_storage)) call clear_storage(s)
         deallocate(lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1)
         deallocate(B1, soln1, ipivot1)
         deallocate(row_scale_factors1, col_scale_factors1)
         deallocate(equed1)
      end subroutine cleanup

   end subroutine time_bcyclic

   subroutine setup_test_matrix(nvar, nz, lblk, dblk, ublk, brhs)
      integer, intent(in) :: nvar, nz
      real(dp), pointer, intent(inout) :: lblk(:,:,:), dblk(:,:,:), ublk(:,:,:)
      real(dp), pointer, intent(inout) :: brhs(:,:)
      integer :: i, j, k

      do k = 1, nz
         do j = 1, nvar
            do i = 1, nvar
               if (i == j) then
                  dblk(i, j, k) = 4.0_dp + 0.1_dp * (i + k)
               else
                  dblk(i, j, k) = 0.01_dp * abs(i - j)
               end if
               lblk(i, j, k) = -1.0_dp + 0.001_dp * i
               ublk(i, j, k) = -1.0_dp + 0.001_dp * j
            end do
            brhs(j, k) = 1.0_dp + 0.1_dp * j + 0.01_dp * k
         end do
      end do

      lblk(:, :, 1) = 0.0_dp
      ublk(:, :, nz) = 0.0_dp

   end subroutine setup_test_matrix

end module time_star_bcyclic
