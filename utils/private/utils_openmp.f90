! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module utils_openmp

      implicit none

      private
      public :: eval_OMP_GET_THREAD_NUM
      public :: eval_OMP_GET_MAX_THREADS
      public :: eval_OMP_SET_NUM_THREADS

      integer :: omp_max_threads = -1

      contains

      integer function eval_OMP_GET_THREAD_NUM()
         use omp_lib, only: OMP_GET_THREAD_NUM
         eval_OMP_GET_THREAD_NUM = OMP_GET_THREAD_NUM()
      end function eval_OMP_GET_THREAD_NUM

      integer function eval_OMP_GET_MAX_THREADS()
         use omp_lib, only: OMP_GET_MAX_THREADS
         if (omp_max_threads < 0) &
            omp_max_threads = OMP_GET_MAX_THREADS()
         eval_OMP_GET_MAX_THREADS = omp_max_threads
      end function eval_OMP_GET_MAX_THREADS

      subroutine eval_OMP_SET_NUM_THREADS(threads)
         use omp_lib, only: OMP_SET_NUM_THREADS
         integer, intent(in) :: threads
         call OMP_SET_NUM_THREADS(threads)
      end subroutine eval_OMP_SET_NUM_THREADS

      end module utils_openmp
