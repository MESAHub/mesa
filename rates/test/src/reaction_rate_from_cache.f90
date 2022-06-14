! ***********************************************************************
!
!   Copyright (C) 2011-2020  Bill Paxton & The MESA Team
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


      program show_rates
      use rates_lib
      use const_lib
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none
      
      character (len=256) :: cache_filename
      integer :: ierr, n
      character (len=32) :: my_mesa_dir
      
      n = COMMAND_ARGUMENT_COUNT()
      if (n /= 1) then
         write(*,*) 'please give full path name of cache file on command line'
         call mesa_error(__FILE__,__LINE__)
      end if
      call GET_COMMAND_ARGUMENT(1, cache_filename)
      write(*,'(a)') '# rates from ' // trim(cache_filename)
      my_mesa_dir = ''         
      call const_init(my_mesa_dir,ierr)     
      if (ierr /= 0) then
         write(*,*) 'const_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if        
      call math_init()
      
      ierr = 0
      call show_reaction_rates_from_cache(cache_filename, ierr) 
      if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
      
      end program show_rates
