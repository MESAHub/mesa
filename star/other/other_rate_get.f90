! ***********************************************************************
!
!   Copyright (C) 2022  The MESA team
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
 
   module other_rate

   implicit none

   ! set use_other_rate_get = .true. in your controls namelist   
        
   ! edit the extras_controls routine to set the procedure pointer
   !  other_rate_get => my_rate_get

   ! This hook is primarily aimed at either modifying an existing mesa reaction or 
   ! providing your own custom rate that depends only on temperature (density factors are handled elsewhere)
   ! This hook can only handle existing reactions. This is not the place to try and add a new reaction.

   ! This is also called early in the setup of MESA so there is no stellar structure yet
   ! Thus we call this multiple times (each with a different temperature) for each reaction.

   ! Note this will effect the cached rate written to data/rates_data/cache 
   ! Becuase of this it will only get called if the rate does NOT already exist in your rates_cache
   ! I recommend you use the rates_cache_dir option to redirect your rates_cache when using this hook
   ! So you dont break your whole mesa install.

   contains


   subroutine default_other_rate_get(ir, temp, tf, raw_rate, ierr)
      use rates_def
      use rates_lib
      implicit none

      integer :: ir ! Rate id
      real(dp),intent(in) ::    temp      !< Temperature
      type (T_Factors) :: tf !< Various temperature factors
      real(dp),intent(inout) ::   raw_rate     !< Unscreened reaction_rate, note this will have the default mesa rate on entry
      integer, intent(out) ::   ierr
   
      ierr = 0

      if (trim(reaction_name(ir)) == 'r_c12_ag_o16') then
         raw_rate = 1.0 * raw_rate
      end if

   end subroutine default_other_rate_get

   end module other_rate
        
        
        
        
  