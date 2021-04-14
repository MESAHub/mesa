! ***********************************************************************
!
!   Copyright (C) 2010-2021  Bill Paxton & The MESA Team
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

      module mod_typical_charge

      use const_def, only: one_third, two_thirds, dp
      use math_lib

      implicit none
      

      logical, parameter :: dbg = .false.
     

      contains

      ! average ionic charges based on
      ! C. Paquette, C. Pelletier, G. Fontaine, G. Michaud,
      ! "Diffusion in White Dwarfs: New Results and Comparative Study",
      ! ApJ Supp. Series, 61:197-217, 1986.

      ! Originally used eqn 21 of above paper for depression of the continuum,
      ! but that expression was missing a rho^1/3, so now corrected to match
      ! eqn 3 of
      ! J. Dupuis, G. Fontaine, C. Pelletier, F. Wesemael,
      ! "A Study of Metal Abundance Patterns in Cool White Dwarfs I.
      ! Time-dependent Calculations of Gravitational Settling",
      ! Apj Supp. Series, 82:505-521, 1992

      ! uses ionization potentials from
      ! Allen, C.W., 1973, "Astrophysical Quantities", 3rd edition, pg 37-38.
      
      ! tables go from H to Zn.
      
      ! the routine is written so that it doesn't ever return 0 as a typical charge.
      ! these are "typical" charges rather than average.  the values are whole numbers.

      real(dp) function eval_typical_charge( &
            cid, abar, free_e, T, log10_T, rho, log10_rho)
         integer, intent(in) :: cid ! chem id such as ihe4.  defined in chem_def.
         real(dp), intent(in) :: abar ! average atomic weight (from chem_lib)
         real(dp), intent(in) :: free_e 
            ! mean number of free electrons per nucleon (from eos_lib)
            ! abar*free_e = (nucleons/particle)*(charge/nucleon) = charge/particle = z1
         real(dp), intent(in) :: T, log10_T, rho, log10_rho
         eval_typical_charge = get_typical_charge( &
            cid, abar, abar*free_e, T, log10_T, rho, log10_rho)    
      end function eval_typical_charge

      
      subroutine chi_info(a1, z1, T, log_T, rho, log_rho, chi, c0, c1, c2)
         real(dp), intent(in) :: a1, z1, T, log_T, rho, log_rho
         real(dp), intent(out) :: chi, c0, c1, c2
         chi = 1.987d-4*T*(-8.392d0 - log_rho + 1.5d0*log_T - log10(z1/a1)) ! eqn 20
         ! coef's used in eqn 21
         c0 = 1.085d-4*rho*T/a1
         c1 = 1.617d4*sqrt(rho*(z1*z1 + z1)/(T*a1))
         c2 = 29.38d0*z1*pow(rho/a1,one_third)
         ! c2 had a typo in eqn 21, now corrected to match Dupuis et al. (1992) eqn 3
      end subroutine chi_info
      
      real(dp) function chi_effective(chi, c0, c1, c2, z1, z2)
         real(dp), intent(in) :: chi, c0, c1, c2, z1, z2
         chi_effective = chi + c0/(z2*z2*z2) + &
            min(c1*z2, c2*(pow(z2/z1,two_thirds) + 0.6d0))
      end function chi_effective
      
      real(dp) function get_typical_charge(cid, a1, z1, T, log_T, rho, log_rho)
         use ionization_potentials
         use chem_def
         integer, intent(in) :: cid
         real(dp), intent(in) :: a1, z1, T, log_T, rho, log_rho      
         real(dp) :: chi, c0, c1, c2, z2, chi_eff
         integer :: i, izmax
         include 'formats'
         if (.not. ionization_tables_okay) then
            call set_ionization_potentials
            ionization_tables_okay = .true.
         end if
         izmax = int(chem_isos% Z(cid))
         get_typical_charge = dble(izmax)
         if (izmax > 30) return
         call chi_info(a1, z1, T, log_T, rho, log_rho, chi, c0, c1, c2)
         do i=1, izmax-1
            z2 = dble(i)
            chi_eff = chi_effective(chi, c0, c1, c2, z1, z2+1)
            if (chi_eff < ip(izmax,i+1)) then
               get_typical_charge = z2
               return
            end if
         end do
      end function get_typical_charge
      
      



      end module mod_typical_charge

