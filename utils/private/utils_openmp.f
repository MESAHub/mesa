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

      module utils_openmp
      
      implicit none
      

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

