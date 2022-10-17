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


module test_screen
   
   use rates_lib
   use rates_def
   use utils_lib, only : mesa_error
   
   implicit none

contains
   
   
   subroutine do_test_screen
      use chem_def
      use chem_lib
      use const_lib
      use math_lib
      
      integer, parameter :: num_isos = 8, max_z_to_cache = 12
      integer :: chem_id(num_isos), i1, i2, jscr, ierr
      integer, pointer :: net_iso(:)
      real(dp) :: xin(num_isos), y(num_isos), iso_z(num_isos), xz, abar, zbar, z2bar, z53bar, ye, sumx, &
         dabar_dx(num_isos), dzbar_dx(num_isos), temp, den, logT, logRho, &
         sc1a, sc1adt, sc1add, xh, xhe, dmc_dx(num_isos), iso_z158(num_isos)
      type (Screen_Info) :: sc
      real(dp) :: zg1, zg2, zg3, zg4
      real(dp) :: zs13, zhat, zhat2, lzav, aznut, zs13inv, mass_correction!approx_abar, approx_zbar
      integer :: screening_mode, i
      character(len = 256) :: scr_option_str
      integer :: h1, he3, he4, c12, n14, o16, ne20, mg24
      character (len = 32) :: my_mesa_dir
      
      include 'formats'
      
      ierr = 0
      my_mesa_dir = '../..'
      call const_init(my_mesa_dir, ierr)
      if (ierr /= 0) then
         write(*, *) 'const_init failed'
         call mesa_error(__FILE__, __LINE__)
      end if
      call math_init()
      
      call chem_init('isotopes.data', ierr)
      if (ierr /= 0) then
         write(*, *) 'chem_init failed'
         call mesa_error(__FILE__, __LINE__)
      end if
      
      h1 = 1
      he3 = 2
      he4 = 3
      c12 = 4
      n14 = 5
      o16 = 6
      ne20 = 7
      mg24 = 8
      
      allocate(net_iso(num_chem_isos))
      
      net_iso = 0
      
      net_iso(ih1) = h1; chem_id(h1) = ih1
      net_iso(ihe3) = he3; chem_id(he3) = ihe3
      net_iso(ihe4) = he4; chem_id(he4) = ihe4
      net_iso(ic12) = c12; chem_id(c12) = ic12
      net_iso(in14) = n14; chem_id(n14) = in14
      net_iso(io16) = o16; chem_id(o16) = io16
      net_iso(ine20) = ne20; chem_id(ne20) = ine20
      net_iso(img24) = mg24; chem_id(mg24) = img24
      
      logT = 7.7110722845770692D+00
      logRho = 4.5306372623742392D+00
      
      xin(net_iso(ih1)) = 9.1649493293186402D-22
      xin(net_iso(ihe3)) = 1.8583723770506327D-24
      xin(net_iso(ihe4)) = 9.8958688392029037D-01
      xin(net_iso(ic12)) = 7.3034226840994307D-04
      xin(net_iso(in14)) = 6.2231918940828723D-03
      xin(net_iso(io16)) = 3.6335428720509150D-04
      xin(net_iso(ine20)) = 1.0527812894109640D-03
      xin(net_iso(img24)) = 2.0434463406007555D-03
      
      i1 = ihe4
      i2 = ic12
      
      if (.false.) then  ! TESTING
         
         xin = 0
         xin(net_iso(ih1)) = 0.72d0
         xin(net_iso(ihe4)) = 0.26d0
         xin(net_iso(in14)) = 0.02d0
         
         i1 = ih1
         i2 = in14
         
         write(*, 1) 'sum(xin)', sum(xin(:))
         
         logT = 7d0
         logRho = 1d0
      
      end if
      
      call composition_info(&
         num_isos, chem_id, xin, xh, xhe, xz, abar, zbar, z2bar, z53bar, &
         ye, mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)
      
      iso_z(:) = chem_isos% Z(chem_id(:))
      
      do i = 1, num_isos
         iso_z158(i) = pow(real(chem_isos% Z(chem_id(i)), kind = dp), 1.58d0)
      end do
      y(:) = xin(:) / chem_isos% Z_plus_N(chem_id(:))
      
      call do1(salpeter_screening, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      
      call do1(extended_screening, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      
      call do1(chugunov_screening, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      
      call do1(no_screening, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      
      deallocate(net_iso)
      
      write(*, *) 'done'
   
   contains
      
      subroutine do1(sc_mode, ierr)
         use math_lib
         integer, intent(in) :: sc_mode
         integer, intent(out) :: ierr
         character (len = 64) :: sc_str
         include 'formats'
         call screening_option_str(sc_mode, sc_str, ierr)
         if (ierr /= 0) return
         write(*, *) trim(sc_str)
         
         temp = exp10(logT)
         den = exp10(logRho)
         
         call screen_init_AZ_info(&
            chem_isos% W(i1), dble(chem_isos% Z(i1)), &
            chem_isos% W(i2), dble(chem_isos% Z(i2)), &
            zs13, &
            zhat, zhat2, lzav, aznut, zs13inv, &
            ierr)
         if (ierr /= 0) return
         
         
         ! set the context for the reactions
         call screen_set_context(&
            sc, temp, den, logT, logRho, zbar, abar, z2bar, &
            sc_mode, num_isos, y, iso_z158)
         
         call screen_pair(&
            sc, &
            chem_isos% W(i1), dble(chem_isos% Z(i1)), &
            chem_isos% W(i2), dble(chem_isos% Z(i2)), &
            sc_mode, &
            zs13, zhat, zhat2, lzav, aznut, zs13inv, rattab_tlo, &
            sc1a, sc1adt, sc1add, ierr)
         if (ierr /= 0) return
         
         write(*, 1) 'logT = ', logT
         write(*, 1) 'logRho = ', logRho
         write(*, 1) 'zbar = ', zbar
         write(*, 1) 'abar = ', abar
         write(*, 1) 'z2bar = ', z2bar
         write(*, 1) 'sc1a = ', sc1a
         write(*, 1) 'sc1adt = ', sc1adt
         write(*, 1) 'sc1add = ', sc1add
         write(*, '(A)')
      
      end subroutine do1
   
   end subroutine do_test_screen

end module test_screen

