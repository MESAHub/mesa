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

      module neu_lib
      ! library for calculating neutrino losses from non-nuclear-burning sources
      ! neutrino losses that occur during nuclear reactions are included in the nuclear library
      ! the data interface for the library is defined in neu_def
      
      use const_def, only: dp
      
      implicit none


      contains ! the procedure interface for the library
      ! client programs should only call these routines.

      
      subroutine neu_get(T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
               loss, sources, info)
         use neu_def
         use mod_neu, only : neutrinos

         ! this routine computes neutrino losses from the analytic fits of
         ! itoh et al. apjs 102, 411, 1996, and also returns their derivatives. 
      
         ! provide T or log10_T or both (the code needs both, so pass 'em if you've got 'em!)
         ! same for Rho and log10_Rho
      
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: log10_T ! log10 of temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: log10_Rho ! log10 of density
         real(dp), intent(in) :: abar ! mean atomic weight
         real(dp), intent(in) :: zbar ! mean charge
         real(dp), intent(in) :: log10_Tlim 
         ! log10 of temperature at which begin to cutoff results
         !    NOTE: the Itoh et al data has a lower temperature limit of 10^7
         !    so for T < 10^7, the neutrino losses are simply set to 0
         !    Rather than have an abrupt cutoff, the values are multiplied by a coefficient
         !    that reaches 0 at T = 10^7 and is equal to 1 for log10T > log10_Tlim
         ! log10_Tlim of 7.5 is a reasonable choice.
         logical, intent(in) :: flags(num_neu_types) ! true if should include the type of loss

         real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
         real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
         integer, intent(out) :: info ! 0 means AOK.
         
         call neutrinos(T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim,  &
                  flags, loss, sources, info)
         
      end subroutine neu_get
      

      end module neu_lib

