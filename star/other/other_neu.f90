! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
 
      module other_neu

      ! consult star/other/README for general usage instructions
      ! control name: use_other_neu = .true.
      ! procedure pointer: s% other_neu => my_routine


      use star_def

      implicit none
      
            
      contains
      
      
      subroutine null_other_neu(  &
            id, k, T, log10_T, Rho, log10_Rho, abar, zbar, z2bar, log10_Tlim, flags, &
            loss, sources, ierr)
         use neu_lib, only: neu_get
         use neu_def
         integer, intent(in) :: id ! id for star         
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: log10_T ! log10 of temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: log10_Rho ! log10 of density
         real(dp), intent(in) :: abar ! mean atomic weight
         real(dp), intent(in) :: zbar ! mean charge
         real(dp), intent(in) :: z2bar ! mean charge squared
         real(dp), intent(in) :: log10_Tlim 
         logical, intent(inout) :: flags(num_neu_types) ! true if should include the type of loss
         real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
         real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
         integer, intent(out) :: ierr
         call neu_get(  &
            T, log10_T, Rho, log10_Rho, abar, zbar, z2bar, log10_Tlim, flags, &
            loss, sources, ierr)
      end subroutine null_other_neu


      end module other_neu
      
      
      
      
