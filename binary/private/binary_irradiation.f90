! ***********************************************************************
!
!   Copyright (C) 2010-2019  Pablo Marchant & The MESA Team
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


module binary_irradiation
   
   use const_def
   use star_lib
   use star_def
   use math_lib, only : safe_log10, pow2
   use binary_def
   
   implicit none

contains
   
   subroutine adjust_irradiation(b)
      type (binary_info), pointer :: b
      real(dp) :: Lx
      integer :: ierr
      type (star_info), pointer :: s
      include 'formats'
      ierr = 0
      s => b% s_donor
      if (b% col_depth_for_eps_extra <= 0) return
      if (b% accretion_powered_irradiation) then
         Lx = b% accretion_luminosity
         s% irradiation_flux = min(b% max_F_irr, Lx / (4 * pi * pow2(b% separation)))
         write(*, 2) 'lg F_irr', s% model_number, safe_log10(s% irradiation_flux)
      else if (b% use_accretor_luminosity_for_irrad) then
         if (.not. b% evolve_both_stars) then
            write(*, *) "Can't use accretor luminosity for irradiation without evolving both stars."
            return
         end if
         Lx = (b% s_accretor% L_phot) * Lsun
         s% irradiation_flux = min(b% max_F_irr, Lx / (4 * pi * pow2(b% separation)))
         write(*, 2) 'lg F_irr', s% model_number, safe_log10(s% irradiation_flux)
      else
         if (b% irrad_flux_at_std_distance <= 0) return
         s% irradiation_flux = b% irrad_flux_at_std_distance * &
            pow2(b% std_distance_for_irradiation / b% separation)
      end if
      s% column_depth_for_irradiation = b% col_depth_for_eps_extra
   end subroutine adjust_irradiation

end module binary_irradiation
