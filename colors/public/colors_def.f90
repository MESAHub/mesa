! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module colors_def
      use const_def, only : strlen, dp
      implicit none
      
      !Public constants for use by clients
      !Have we called colors_init yet?
      logical :: color_is_initialized=.false.
      
      integer, parameter :: max_num_color_files=10
      integer, parameter :: max_num_bcs_per_file=20 
      integer :: bc_total_num_colors
      
      ! color indices are differences in magnitudes in different wavelength bands
      ! as a reminder for non-experts like myself, here's how it goes
      !
      ! msun := apparent magnitude of sun is -26.81
      ! Fsun := solar flux at 1AU is 1.36e6 erg/s/cm^2
      !
      ! "apparent magnitude" m of star with flux F is m = msun - 2.5 log10(F/Fsun)
      ! "absolute magnitude" M for star of apparent magnitude m at distance d is M = m - 5 log(d/d0)
      !     where the standard distance d0 is 10pc.
      !     i.e., absolute magnitude is what the apparent magnitude would be if star were at 10 parsecs.
      !
      ! thus absolute magnitude of sun is about 4.75
      ! 
      ! "bolometric magnitude" = absolute magnitude using flux integrated over all wavelengths
      !     can be derived from the current stellar luminosity using the equation
      !     log(Lstar/Lsun) = (Mbol_sun - Mbol_star)/2.5 using Mbol_sun = 4.75 (LCB)
      !
      ! "visual magnitude" = absolute magnitude only using flux in visible wavelengths
      !      more precisely, this is magnitude as measured with filter centered at 5500A, 890A width.
      !
      ! "bolometric correction" = bolometric magnitude minus visual magnitude
      !      for the sun, the bolometric correction is about -0.11
      !      thus visual magnitude of sun is about 4.86 = Mbol_sun - BC_sun = 4.75 - (-0.11)
      !
      ! in order of increasing wavelength, the "color" magnitudes are as follows:
      !
      ! "U" is the ultraviolet magnitude, center at 365nm.
      ! "B" is the        blue magnitude, center at 440nm.
      ! "V" is the      visual magnitude, center at 550nm.
      ! "R" is the         red magnitude, center at 600nm.
      ! "I" is the   infra-red magnitude, center at 800nm.
      
      ! in addition, longer wavelength "colors" have been defined as well
      ! by order of increasing wavelength, these are J, H, K, L, and M.
      
      ! "color index" is the difference between 2 color magnitudes
      ! for example, B-V is colors_B - colors_V
      ! smaller B-V means larger brightness in blue band compared to visual band, means bluer star.


      ! color magnitude data from Lejeune, Cuisinier, Buser (1998) A&AS 130, 65-75. [LCB]
      ! the coverage is approximately Teff from 50,000K to 2000K, log g 5.5 to -1.02, [Fe/H} 1.0 to -5.0
      !
      ! but not all combination of these are actually represented in the tables.  
      ! the current implementation limits the given arguments to the actual range in the tables.
      ! and it does a simple linear interpolation between tabulated values.

      ! BTW: they use [Fe/H] as a parameter;
      ! the evolution code uses log10(Z/Zsun) as an approximation for this.
      
      ! THE FOLLOWING ARE PRIVATE DEFS -- NOT FOR USE BY CLIENTS
      
      type :: lgz_list ! sorted in decreasing order of lgz ([M/H])
         real(dp) :: lgz ! [Fe_H]
         type (lgz_list), pointer :: nxt => null()
         real(dp),dimension(max_num_bcs_per_file) :: colors
      end type

      type :: lgt_list ! sorted in decreasing order of lgt
         real(dp) :: lgt ! logTeff
         integer :: n_colors
         type (lgt_list), pointer :: nxt => null()
         type (lgg_list), pointer :: glist => null()
      end type

      type :: lgg_list ! sorted in decreasing order of lgg
         real(dp) :: lgg ! log g
         type (lgg_list), pointer :: nxt => null()
         type (lgz_list), pointer :: zlist => null()
      end type
      
      type :: col_list
         !Main data store
         type(lgt_list), pointer :: thead => null()
         CHARACTER(len=strlen),dimension(max_num_bcs_per_file) :: color_names
         integer :: n_colors
      end type

      integer :: num_thead
      type (col_list),dimension(:),pointer :: thead_all => null()
         
      
      end module colors_def

