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
!
! ***********************************************************************

module atm_table
   
   ! Uses
   
   use const_def
   use math_lib
   use utils_lib, only : mesa_error
   
   ! No implicit typing
   
   implicit none
   
   ! Access specifiers
   
   private
   
   public :: eval_table
   public :: get_table_alfa_beta
   public :: get_table_base
   
   ! Procedures

contains
   
   subroutine eval_table(&
      L, R, M, cgrav, id, Z, skip_partials, &
      Teff, &
      lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
      lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
      ierr)
      
      use atm_def
      use eos_lib, only : radiation_pressure
      use table_atm, only : get_table_values
      use utils_lib, only : is_bad
      
      real(dp), intent(in) :: L
      real(dp), intent(in) :: R
      real(dp), intent(in) :: M
      real(dp), intent(in) :: cgrav
      integer, intent(in) :: id
      real(dp), intent(in) :: Z
      logical, intent(in) :: skip_partials
      real(dp), intent(in) :: Teff
      real(dp), intent(out) :: lnT
      real(dp), intent(out) :: dlnT_dL
      real(dp), intent(out) :: dlnT_dlnR
      real(dp), intent(out) :: dlnT_dlnM
      real(dp), intent(out) :: dlnT_dlnkap
      real(dp), intent(out) :: lnP
      real(dp), intent(out) :: dlnP_dL
      real(dp), intent(out) :: dlnP_dlnR
      real(dp), intent(out) :: dlnP_dlnM
      real(dp), intent(out) :: dlnP_dlnkap
      integer, intent(out) :: ierr
      
      logical, parameter :: DBG = .FALSE.
      
      real(dp) :: g
      real(dp) :: logg
      real(dp) :: Pgas
      real(dp) :: dPgas_dTeff
      real(dp) :: dPgas_dlogg
      real(dp) :: T
      real(dp) :: dT_dTeff
      real(dp) :: dT_dlogg
      real(dp) :: Prad
      real(dp) :: P
      real(dp) :: dlnTeff_dL
      real(dp) :: dlnTeff_dlnR
      real(dp) :: dTeff_dlnR
      real(dp) :: dlogg_dlnR
      real(dp) :: dlogg_dlnM
      real(dp) :: dPgas_dlnR
      real(dp) :: dPgas_dlnM
      real(dp) :: dPgas_dlnTeff
      real(dp) :: dPrad_dlnTeff
      real(dp) :: dPrad_dlnR
      real(dp) :: dPrad_dL
      real(dp) :: dT_dlnR
      real(dp) :: dT_dlnM
      real(dp) :: dT_dlnT
      real(dp) :: dT_dlnTeff
      
      include 'formats'
      
      ierr = 0
      
      ! Sanity checks
      
      if (L <= 0._dp .OR. R <= 0._dp .OR. M <= 0._dp) then
         ierr = -1
         return
      end if
      
      ! Evaluate the gravity
      
      g = cgrav * M / (R * R)
      logg = log10(g)
      
      ! Perform the table lookup
      
      if (DBG) write(*, *) 'call get_table_values', id
      call get_table_values(&
         id, Z, logg, Teff, &
         Pgas, dPgas_dTeff, dPgas_dlogg, &
         T, dT_dTeff, dT_dlogg, &
         ierr)
      if (ierr /= 0) then
         if (DBG) write(*, *) 'get_table_values(_at_fixed_Z) ierr', ierr
         return
      end if
      
      ! Set up lnT and lnP
      
      lnT = log(T)
      
      Prad = radiation_pressure(T)
      P = Pgas + Prad
      lnP = log(P)
      
      ! Set up partials
      
      if (.NOT. skip_partials) then
         
         dlnTeff_dlnR = -0.5_dp
         dlnTeff_dL = 0.25_dp / L
         dTeff_dlnR = Teff * dlnTeff_dlnR
         
         dlogg_dlnR = -2._dp / ln10
         dlogg_dlnM = 1._dp / ln10
         
         dPgas_dlnR = dPgas_dlogg * dlogg_dlnR + dPgas_dTeff * dTeff_dlnR
         dPgas_dlnM = dPgas_dlogg * dlogg_dlnM
         dPgas_dlnTeff = dPgas_dTeff * Teff
         
         dPrad_dlnTeff = 4._dp * Prad * dT_dTeff * Teff / T
         dPrad_dlnR = dPrad_dlnTeff * dlnTeff_dlnR
         dPrad_dL = dPrad_dlnTeff * dlnTeff_dL
         
         dlnP_dL = (dPgas_dlnTeff + dPrad_dlnTeff) * dlnTeff_dL / P
         dlnP_dlnR = (dPgas_dlnR + dPrad_dlnTeff * dlnTeff_dlnR) / P
         dlnP_dlnM = dPgas_dlnM / P
         dlnP_dlnkap = 0._dp
         
         dT_dlnTeff = dT_dTeff * Teff
         dT_dlnR = dT_dlogg * dlogg_dlnR + dT_dlnTeff * dlnTeff_dlnR
         dT_dlnM = dT_dlogg * dlogg_dlnM
         
         dlnT_dL = dT_dlnTeff * dlnTeff_dL / T
         dlnT_dlnR = dT_dlnR / T
         dlnT_dlnM = dT_dlnM / T
         dlnT_dlnkap = 0._dp
      
      else
         
         dlnP_dL = 0._dp
         dlnP_dlnR = 0._dp
         dlnP_dlnM = 0._dp
         dlnP_dlnkap = 0._dp
         
         dlnT_dL = 0._dp
         dlnT_dlnR = 0._dp
         dlnT_dlnM = 0._dp
         dlnT_dlnkap = 0._dp
      
      endif
      
      if (DBG .or. is_bad(lnP) .or. is_bad(lnT)) then
         write(*, *) 'eval_table'
         write(*, 1) 'Teff', Teff
         write(*, 1) 'T', T
         write(*, 1) 'dT_dTeff', dT_dTeff
         write(*, 1) 'dT_dlogg', dT_dlogg
         write(*, '(A)')
         ierr = -1
         return
         !if (is_bad(lnP) .or. is_bad(lnT)) call mesa_error(__FILE__,__LINE__,'eval_table')
      end if
      
      ! Finish
      
      return
   
   end subroutine eval_table
   
   !****
   
   subroutine get_table_alfa_beta(&
      L, Teff, R, M, cgrav, id, alfa, beta, ierr)
      
      use atm_def
      use table_atm, only : &
         ai_two_thirds, ai_100, ai_10, ai_1, ai_1m1, ai_wd_25, ai_db_wd_25
      
      real(dp), intent(in) :: L
      real(dp), intent(in) :: Teff
      real(dp), intent(in) :: R
      real(dp), intent(in) :: M
      real(dp), intent(in) :: cgrav
      integer, intent(in) :: id
      real(dp), intent(out) :: alfa
      real(dp), intent(out) :: beta
      integer, intent(out) :: ierr
      
      integer, parameter :: PURE_GREY = 1
      integer, parameter :: PURE_TABLE = 2
      integer, parameter :: BLEND_IN_X = 3
      integer, parameter :: BLEND_IN_Y = 4
      integer, parameter :: BLEND_CORNER_OUT = 5
      
      logical, parameter :: DBG = .FALSE.
      
      real(dp) :: g
      real(dp) :: logTeff
      real(dp) :: logg
      type(Atm_Info), pointer :: ai
      real(dp) :: logg_max
      real(dp) :: logg_min
      real(dp) :: logTeff_max
      real(dp) :: logTeff_min
      real(dp) :: Teff_max
      real(dp) :: Teff_min
      real(dp) :: logg_min_margin
      real(dp) :: logg_max_margin
      real(dp) :: logTeff_min_margin
      real(dp) :: logTeff_max_margin
      real(dp) :: logg1
      real(dp) :: logg2
      real(dp) :: logg3
      real(dp) :: logg4
      real(dp) :: logTeff1
      real(dp) :: logTeff2
      real(dp) :: logTeff3
      real(dp) :: logTeff4
      real(dp) :: c_dx
      real(dp) :: c_dy
      integer :: iregion
      integer :: logg_index
      integer :: ng
      integer :: j
      
      include 'formats'
      
      ierr = 0
      
      ! Sanity checks
      
      if (L <= 0._dp .OR. R <= 0._dp .OR. M <= 0._dp) then
         ierr = -1
         return
      end if
      
      ! Evaluate the gravity
      
      g = cgrav * M / (R * R)
      logTeff = log10(Teff)
      logg = log10(g)
      
      ! Set up the table pointer
      
      select case (id)
      case (ATM_TABLE_PHOTOSPHERE)
         ai => ai_two_thirds
      case (ATM_TABLE_TAU_100)
         ai => ai_100
      case (ATM_TABLE_TAU_10)
         ai => ai_10
      case (ATM_TABLE_TAU_1)
         ai => ai_1
      case (ATM_TABLE_TAU_1M1)
         ai => ai_1m1
      case (ATM_TABLE_WD_TAU_25)
         ai => ai_wd_25
      case (ATM_TABLE_DB_WD_TAU_25)
         ai => ai_db_wd_25
      case default
         write(*, *) 'Invalid id in get_table_alfa_beta:', id
         call mesa_error(__FILE__, __LINE__)
      end select
      
      ng = ai% ng
      logg_max = ai% logg_array(ng)
      logg_min = ai% logg_array(1)
      
      ! First, locate current point logg array for use with Teff_bound
      
      if (logg <= logg_min) then
         logg_index = 1
      else if (logg >= logg_max) then
         logg_index = ng
      else
         logg_index = ng
         do j = 1, ng - 1
            if (ai% logg_array(j) <= logg) then
               logg_index = j
            end if
         end do
      end if
      
      Teff_max = ai% Teff_bound(logg_index)
      Teff_min = ai% Teff_array(1)
      logTeff_max = log10(Teff_max)
      logTeff_min = log10(Teff_min)
      
      ! Set up margins for blending
      
      logg_max_margin = 0.01d0
      logg_min_margin = 0.5d0
      logTeff_max_margin = 0.01d0
      logTeff_min_margin = 0.01d0
      
      ! if (id == ATM_TABLE_WD_TAU_25) then
      !    logg_max_margin = 0.5d0
      !    logg_min_margin = 0.5d0
      !    logTeff_max_margin = 0.2d0
      !    logTeff_min_margin = 1d0
      !    logg1 = logg_max + logg_max_margin
      !    logg2 = logg_max
      !    logg3 = logg_min
      !    logg4 = logg_min - logg_min_margin
      !    logTeff1 = logTeff_max + logTeff_max_margin
      !    logTeff2 = logTeff_max
      !    logTeff3 = logTeff_min
      !    logTeff4 = logTeff_min - logTeff_min_margin
      ! else
      logg1 = logg_max
      logg2 = logg_max - logg_max_margin
      logg3 = logg_min + logg_min_margin
      logg4 = logg_min
      logTeff1 = logTeff_max
      logTeff2 = logTeff_max - logTeff_max_margin
      logTeff3 = logTeff_min
      logTeff4 = logTeff_min
      ! end if
      
      ! Decide on what sort of region we're in
      
      if (logg < logg4 .or. logg > logg1 .or. logTeff > logTeff1) then
         iregion = pure_grey
      else if (logTeff > logTeff2) then
         c_dy = (logTeff - logTeff2) / (logTeff1 - logTeff2)
         if (logg > logg2) then
            c_dx = (logg - logg2) / (logg1 - logg2)
            iregion = blend_corner_out
         else if (logg > logg3) then
            iregion = blend_in_y
         else ! logg > logg4
            c_dx = (logg - logg3) / (logg4 - logg3)
            iregion = blend_corner_out
         end if
      else if (logTeff >= logTeff3) then
         if (logg > logg2) then
            c_dx = (logg - logg2) / (logg1 - logg2)
            iregion = blend_in_x
         else if (logg > logg3) then
            iregion = pure_table
         else ! logg > logg4
            c_dx = (logg - logg3) / (logg4 - logg3)
            iregion = blend_in_x
         end if
      else if (logTeff > logTeff4) then
         c_dy = (logTeff - logTeff3) / (logTeff4 - logTeff3)
         if (logg > logg2) then
            c_dx = (logg - logg2) / (logg1 - logg2)
            iregion = blend_corner_out
         else if (logg > logg3) then
            iregion = blend_in_y
         else ! logg > logg4
            c_dx = (logg - logg3) / (logg4 - logg3)
            iregion = blend_corner_out
         end if
      else ! logTeff <= logTeff4
         iregion = pure_grey
      end if
      
      ! Set up alfa and beta
      
      select case (iregion)
      case (pure_grey)
         alfa = 1._dp
         beta = 0._dp
      case (pure_table)
         alfa = 0._dp
         beta = 1._dp
      case (blend_in_y)
         alfa = c_dy
         beta = 1._dp - alfa
      case (blend_in_x)
         alfa = c_dx
         beta = 1._dp - alfa
      case (blend_corner_out)
         alfa = min(1d0, sqrt(c_dx * c_dx + c_dy * c_dy))
         beta = 1 - alfa
      case default
         write(*, *) 'Invalid iregion in get_table_alfa_beta:', iregion
         call mesa_error(__FILE__, __LINE__)
      end select
      
      ! Finish
      
      return
   
   end subroutine get_table_alfa_beta
   
   !****
   
   subroutine get_table_base (id, tau_base, ierr)
      
      use atm_def
      
      integer, intent(in) :: id
      real(dp), intent(out) :: tau_base
      integer, intent(out) :: ierr
      
      ierr = 0
      
      ! Get the base optical depth for the atmosphere
      
      select case (id)
      case (ATM_TABLE_PHOTOSPHERE)
         tau_base = 2._dp / 3._dp
      case (ATM_TABLE_TAU_100)
         tau_base = 100._dp
      case (ATM_TABLE_TAU_10)
         tau_base = 10._dp
      case (ATM_TABLE_TAU_1)
         tau_base = 1._dp
      case (ATM_TABLE_TAU_1M1)
         tau_base = 0.1_dp
      case (ATM_TABLE_WD_TAU_25)
         tau_base = 25.1188_dp
      case (ATM_TABLE_DB_WD_TAU_25)
         tau_base = 25._dp
      case default
         write(*, *) 'Invalid id in get_table_base:', id
         call mesa_error(__FILE__, __LINE__)
      end select
      
      ! Finish
      
      return
   
   end subroutine get_table_base

end module atm_table
