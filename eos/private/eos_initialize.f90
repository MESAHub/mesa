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
!
! ***********************************************************************

      module eos_initialize

      implicit none

      contains
            
      
      subroutine Init_eos( &
           eosDT_cache_dir_in, &
           use_cache, ierr)
         use eos_def
         use helm_alloc
         use utils_lib, only : mkdir
         use const_def, only: mesa_data_dir, mesa_caches_dir, mesa_temp_caches_dir
         character(*), intent(IN) :: eosDT_cache_dir_in
         logical, intent(in) :: use_cache
         integer, intent(OUT) :: ierr ! 0 means AOK.
         !integer, parameter :: imax = 261, jmax = 101  
            ! dimensions of small version of helm table
         !integer, parameter :: imax = 1081, jmax = 401
            ! dimensions of medium version of helm table; 40 points per decade
         integer, parameter :: imax = 2701, jmax = 1001
            ! dimensions of large version of helm table; 100 points per decade
         ! helm table lives in eosDT_data
         character (len=256) :: eosDT_data_dir
         ierr = 0
         if (eos_root_is_initialized) return
         use_cache_for_eos = use_cache
         eosDT_data_dir = trim(mesa_data_dir) // '/eosDT_data'
         if (use_cache_for_eos) then
            if (len_trim(eosDT_cache_dir_in) > 0) then
               eosDT_cache_dir = eosDT_cache_dir_in
            else if (len_trim(mesa_caches_dir) > 0) then
               eosDT_cache_dir = trim(mesa_caches_dir) // '/eosDT_cache'
            else
               eosDT_cache_dir = trim(eosDT_data_dir) // '/cache'
            end if
            call mkdir(eosDT_cache_dir)
            eosDT_temp_cache_dir = trim(mesa_temp_caches_dir) // '/eosDT_cache'
            if(use_mesa_temp_cache) call mkdir(eosDT_temp_cache_dir)
         end if
         
         call alloc_helm_table(eos_ht, imax, jmax, ierr)
         if (ierr /= 0) return
         
         call read_helm_table(eos_ht, &
            eosDT_data_dir, eosDT_cache_dir, eosDT_temp_cache_dir, use_cache_for_eos, ierr)
         if (ierr /= 0) return

         call eos_def_init
         ! replace defaults from eos_def_init by argument
         
         eos_root_is_initialized = .true.
      
      end subroutine Init_eos
      
      
      end module eos_initialize
