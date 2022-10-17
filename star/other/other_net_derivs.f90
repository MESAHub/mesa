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

module other_net_derivs
   
   implicit none
   
   ! set use_other_net_derivs = .true. in your controls namelist   
   
   ! edit the extras_controls routine to set the procedure pointer
   !  other_net_derivs => my_net_derivs
   
   ! This hook is aimed at being able to modify dydt term.
   ! This means you can add arbitary new reactions that MESA does not know about
   ! by changing the dydt (change in compostion with time)
   
   ! This hook only works with soft nets (so no approx networks)

contains
   
   
   subroutine default_other_net_derivs(&
      n, dydt, eps_nuc_MeV, eta, ye, logtemp, temp, den, abar, zbar, &
      num_reactions, rate_factors, &
      symbolic, just_dydt, ierr)
      use net_def
      
      implicit none
      
      type(Net_Info) :: n
      real(qp), pointer, intent(inout) :: dydt(:, :)
      real(qp), intent(out) :: eps_nuc_MeV(:)
      integer, intent(in) :: num_reactions
      real(dp), intent(in) :: eta, ye, logtemp, temp, den, abar, zbar, &
         rate_factors(:)
      logical, intent(in) :: symbolic, just_dydt
      integer, intent(out) :: ierr
   
   end subroutine default_other_net_derivs

end module other_net_derivs
        
        
        
