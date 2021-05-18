! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

      module screen
      use const_def
      use rates_def
      use math_lib
      
      implicit none
      
      contains
      
      subroutine do_screen_set_context( &
            sc, temp, den, logT, logRho, zbar, abar, z2bar, &
            screening_mode, num_isos, y, iso_z158)
         type (Screen_Info) :: sc
         integer, intent(in) :: num_isos
         real(dp), intent(in) ::  &
            temp, den, logT, logRho, zbar, abar, z2bar,  &
            y(:), &
            iso_z158(:) ! Z**1.58
         integer, intent(in) :: screening_mode
      
         real(dp), parameter :: x13   = 1.0d0/3.0d0 
         real(dp), parameter :: x14   = 1.0d0/4.0d0
         real(dp), parameter :: x53   = 5.0d0/3.0d0
         real(dp), parameter :: x532  = 5.0d0/32.0d0
         real(dp), parameter :: x512  = 5.0d0/12.0d0
         real(dp), parameter :: fact  = 1.25992104989487d0 ! the cube root of 2
         real(dp), parameter :: co2   = x13 * 4.248719d3
         real(dp) :: qq
         integer :: j
      
         logical, parameter :: debug = .false.
         !logical, parameter :: debug = .true.
      
         include 'formats'

         if (screening_mode == no_screening .or. zbar == 0d0) return

         sc% temp  = temp
         sc% den   = den
         sc% logT  = logT
         sc% logRho = logRho
         sc% zbar  = zbar
         sc% abar  = abar
         sc% z2bar = z2bar

         ! get the info that depends only on temp, den, and overall composition         

         sc% ytot     = 1.0d0/abar
         sc% rr       = den * sc% ytot
         sc% tempi    = 1.0d0/temp
         sc% dtempi   = -sc% tempi * sc% tempi
         sc% deni     = 1.0d0/den
         qq = 0d0
         do j=1,num_isos
            if (iso_z158(j) == 0d0) cycle
            qq = qq + iso_z158(j) * y(j)
         end do
         sc% z1pt58bar = abar * qq
         sc% zbar13 = pow(zbar,1d0/3d0)
           
         sc% pp       = sqrt(sc% rr * sc% tempi * (z2bar + zbar)) 
         qq            = 0.5d0/(sc% pp) *(z2bar + zbar) 
         sc% dppdt    = qq*sc% rr*sc% dtempi
         sc% dppdd    = qq*sc% ytot*sc% tempi

         sc% qlam0z   = 1.88d8 * sc% tempi * sc% pp
         sc% qlam0zdt = 1.88d8 * (sc% dtempi * sc% pp + sc% tempi * sc% dppdt)
         sc% qlam0zdd = 1.88d8 * sc% tempi * sc% dppdd

         sc% taufac   = co2 * pow(sc% tempi,x13)
         sc% taufacdt = -x13*sc% taufac*sc% tempi

         qq           = sc% rr*zbar
         sc% xni     = pow(qq,x13)
         sc% dxnidd  = x13 * sc% xni * sc% deni

         sc% aa     = 2.27493d5 * sc% tempi * sc% xni
         sc% daadt  = 2.27493d5 * sc% dtempi * sc% xni
         sc% daadd  = 2.27493d5 * sc% tempi * sc% dxnidd
         
         ! ion and electron sphere radii (itoh 1979 eq 1-3)
         sc% ntot  = den / (amu*abar)
         sc% a_e = pow((3.d0 /(pi4 * zbar * sc% ntot)),x13)
      
      end subroutine do_screen_set_context


      end module screen

