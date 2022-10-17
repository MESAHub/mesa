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

module math_pown
   
   ! Uses
   
   use const_lib, only : dp
   
   ! No implicit typing
   
   implicit none
   
   ! Generic interfaces
   
   interface powm1
      module procedure powm1_
   end interface powm1
   
   interface pow2
      module procedure pow2_
   end interface pow2
   
   interface pow3
      module procedure pow3_
   end interface pow3
   
   interface pow4
      module procedure pow4_
   end interface pow4
   
   interface pow5
      module procedure pow5_
   end interface pow5
   
   interface pow6
      module procedure pow6_
   end interface pow6
   
   interface pow7
      module procedure pow7_
   end interface pow7
   
   interface pow8
      module procedure pow8_
   end interface pow8
   
   ! Access specifiers
   
   private
   
   public :: powm1
   public :: pow2
   public :: pow3
   public :: pow4
   public :: pow5
   public :: pow6
   public :: pow7
   public :: pow8
   
   ! Procedures

contains
   
   elemental function powm1_ (x) result (powm1_x)
      
      real(dp), intent(in) :: x
      real(dp) :: powm1_x
      
      powm1_x = 1_dp / x
   
   end function powm1_
   
   !****
   
   elemental function pow2_ (x) result (pow2_x)
      
      real(dp), intent(in) :: x
      real(dp) :: pow2_x
      
      pow2_x = x * x
   
   end function pow2_
   
   !****
   
   elemental function pow3_ (x) result (pow3_x)
      
      real(dp), intent(in) :: x
      real(dp) :: pow3_x
      
      pow3_x = x * x * x
   
   end function pow3_
   
   !****
   
   elemental function pow4_ (x) result (pow4_x)
      
      real(dp), intent(in) :: x
      real(dp) :: pow4_x
      
      pow4_x = x * x * x * x
   
   end function pow4_
   
   !****
   
   elemental function pow5_ (x) result (pow5_x)
      
      real(dp), intent(in) :: x
      real(dp) :: pow5_x
      
      pow5_x = x * x * x * x * x
   
   end function pow5_
   
   !****
   
   elemental function pow6_ (x) result (pow6_x)
      
      real(dp), intent(in) :: x
      real(dp) :: pow6_x
      
      pow6_x = x * x * x * x * x * x
   
   end function pow6_
   
   !****
   
   elemental function pow7_ (x) result (pow7_x)
      
      real(dp), intent(in) :: x
      real(dp) :: pow7_x
      
      pow7_x = x * x * x * x * x * x * x
   
   end function pow7_
   
   !****
   
   elemental function pow8_ (x) result (pow8_x)
      
      real(dp), intent(in) :: x
      real(dp) :: pow8_x
      
      pow8_x = x * x * x * x * x * x * x * x
   
   end function pow8_

end module math_pown
