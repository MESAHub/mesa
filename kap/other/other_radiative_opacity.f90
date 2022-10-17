! ***********************************************************************
!
!   Copyright (C) 2021 The MESA Team
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

module other_radiative_opacity
   
   ! consult star/other/README for general usage instructions
   ! kap namelist option: use_other_radiative_opacity = .true.
   ! procedure pointers: s% kap_rq % other_radiative_opacity => my_routine
   
   use kap_def
   
   implicit none


contains
   
   subroutine null_other_radiative_opacity(&
      handle, &
      X, Z, XC, XN, XO, XNe, logRho, logT, &
      frac_lowT, frac_highT, frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      
      use const_def, only : dp
      
      ! INPUT
      integer, intent(in) :: handle ! kap handle; from star, pass s% kap_handle
      real(dp), intent(in) :: X, Z, XC, XN, XO, XNe ! composition
      real(dp), intent(in) :: logRho ! density
      real(dp), intent(in) :: logT ! temperature
      
      ! OUTPUT
      real(dp), intent(out) :: frac_lowT, frac_highT, frac_Type2
      real(dp), intent(out) :: kap ! opacity
      real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
      real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
      integer, intent(out) :: ierr ! 0 means AOK.
      
      write(*, *) 'no implementation for other_radiative_opacity'
      ierr = -1
      
      ! can first call kap_lib routine to get standard results, if desired
      
      ! call kap_get_radiative_opacity( &
      !    handle, &
      !    X, Z, XC, XN, XO, XNe, logRho, logT, &
      !    frac_lowT, frac_highT, frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
   
   end subroutine null_other_radiative_opacity

end module other_radiative_opacity
